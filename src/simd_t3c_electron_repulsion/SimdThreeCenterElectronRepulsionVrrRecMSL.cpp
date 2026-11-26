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


#include "SimdThreeCenterElectronRepulsionVrrRecMSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_21 = 4.0 / q;
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

    const auto *lsl0_0 = buffer.data(lsl0 + 0);
    const auto *lsl0_3 = buffer.data(lsl0 + 3);
    const auto *lsl0_5 = buffer.data(lsl0 + 5);
    const auto *lsl0_6 = buffer.data(lsl0 + 6);
    const auto *lsl0_9 = buffer.data(lsl0 + 9);
    const auto *lsl0_10 = buffer.data(lsl0 + 10);
    const auto *lsl0_14 = buffer.data(lsl0 + 14);
    const auto *lsl0_15 = buffer.data(lsl0 + 15);
    const auto *lsl0_20 = buffer.data(lsl0 + 20);
    const auto *lsl0_21 = buffer.data(lsl0 + 21);
    const auto *lsl0_27 = buffer.data(lsl0 + 27);
    const auto *lsl0_36 = buffer.data(lsl0 + 36);
    const auto *lsl0_44 = buffer.data(lsl0 + 44);

    const auto *lsk_0 = buffer.data(lsk + 0);
    const auto *lsk_1 = buffer.data(lsk + 1);
    const auto *lsk_2 = buffer.data(lsk + 2);
    const auto *lsk_3 = buffer.data(lsk + 3);
    const auto *lsk_5 = buffer.data(lsk + 5);
    const auto *lsk_6 = buffer.data(lsk + 6);
    const auto *lsk_9 = buffer.data(lsk + 9);
    const auto *lsk_10 = buffer.data(lsk + 10);
    const auto *lsk_14 = buffer.data(lsk + 14);
    const auto *lsk_15 = buffer.data(lsk + 15);
    const auto *lsk_20 = buffer.data(lsk + 20);
    const auto *lsk_28 = buffer.data(lsk + 28);
    const auto *lsk_30 = buffer.data(lsk + 30);
    const auto *lsk_31 = buffer.data(lsk + 31);
    const auto *lsk_32 = buffer.data(lsk + 32);
    const auto *lsk_33 = buffer.data(lsk + 33);
    const auto *lsk_35 = buffer.data(lsk + 35);
    const auto *lsk_64 = buffer.data(lsk + 64);
    const auto *lsk_66 = buffer.data(lsk + 66);
    const auto *lsk_67 = buffer.data(lsk + 67);
    const auto *lsk_68 = buffer.data(lsk + 68);
    const auto *lsk_69 = buffer.data(lsk + 69);
    const auto *lsk_70 = buffer.data(lsk + 70);
    const auto *lsk_71 = buffer.data(lsk + 71);
    const auto *lsk_100 = buffer.data(lsk + 100);
    const auto *lsk_101 = buffer.data(lsk + 101);
    const auto *lsk_102 = buffer.data(lsk + 102);
    const auto *lsk_103 = buffer.data(lsk + 103);
    const auto *lsk_104 = buffer.data(lsk + 104);
    const auto *lsk_105 = buffer.data(lsk + 105);
    const auto *lsk_107 = buffer.data(lsk + 107);

    const auto *lsl1_0 = buffer.data(lsl1 + 0);
    const auto *lsl1_3 = buffer.data(lsl1 + 3);
    const auto *lsl1_5 = buffer.data(lsl1 + 5);
    const auto *lsl1_6 = buffer.data(lsl1 + 6);
    const auto *lsl1_9 = buffer.data(lsl1 + 9);
    const auto *lsl1_10 = buffer.data(lsl1 + 10);
    const auto *lsl1_14 = buffer.data(lsl1 + 14);
    const auto *lsl1_15 = buffer.data(lsl1 + 15);
    const auto *lsl1_20 = buffer.data(lsl1 + 20);
    const auto *lsl1_21 = buffer.data(lsl1 + 21);
    const auto *lsl1_27 = buffer.data(lsl1 + 27);
    const auto *lsl1_36 = buffer.data(lsl1 + 36);
    const auto *lsl1_44 = buffer.data(lsl1 + 44);

    const auto *msi0_0 = buffer.data(msi0 + 0);
    const auto *msi0_1 = buffer.data(msi0 + 1);
    const auto *msi0_2 = buffer.data(msi0 + 2);
    const auto *msi0_3 = buffer.data(msi0 + 3);
    const auto *msi0_5 = buffer.data(msi0 + 5);
    const auto *msi0_6 = buffer.data(msi0 + 6);
    const auto *msi0_8 = buffer.data(msi0 + 8);
    const auto *msi0_9 = buffer.data(msi0 + 9);
    const auto *msi0_10 = buffer.data(msi0 + 10);
    const auto *msi0_12 = buffer.data(msi0 + 12);
    const auto *msi0_13 = buffer.data(msi0 + 13);
    const auto *msi0_14 = buffer.data(msi0 + 14);
    const auto *msi0_21 = buffer.data(msi0 + 21);
    const auto *msi0_23 = buffer.data(msi0 + 23);
    const auto *msi0_24 = buffer.data(msi0 + 24);
    const auto *msi0_25 = buffer.data(msi0 + 25);
    const auto *msi0_26 = buffer.data(msi0 + 26);
    const auto *msi0_27 = buffer.data(msi0 + 27);
    const auto *msi0_31 = buffer.data(msi0 + 31);
    const auto *msi0_34 = buffer.data(msi0 + 34);
    const auto *msi0_35 = buffer.data(msi0 + 35);
    const auto *msi0_38 = buffer.data(msi0 + 38);
    const auto *msi0_39 = buffer.data(msi0 + 39);
    const auto *msi0_40 = buffer.data(msi0 + 40);
    const auto *msi0_49 = buffer.data(msi0 + 49);
    const auto *msi0_50 = buffer.data(msi0 + 50);
    const auto *msi0_51 = buffer.data(msi0 + 51);
    const auto *msi0_52 = buffer.data(msi0 + 52);
    const auto *msi0_53 = buffer.data(msi0 + 53);
    const auto *msi0_58 = buffer.data(msi0 + 58);
    const auto *msi0_60 = buffer.data(msi0 + 60);
    const auto *msi0_61 = buffer.data(msi0 + 61);
    const auto *msi0_63 = buffer.data(msi0 + 63);
    const auto *msi0_64 = buffer.data(msi0 + 64);
    const auto *msi0_65 = buffer.data(msi0 + 65);
    const auto *msi0_67 = buffer.data(msi0 + 67);
    const auto *msi0_68 = buffer.data(msi0 + 68);
    const auto *msi0_69 = buffer.data(msi0 + 69);
    const auto *msi0_70 = buffer.data(msi0 + 70);
    const auto *msi0_78 = buffer.data(msi0 + 78);
    const auto *msi0_79 = buffer.data(msi0 + 79);

    const auto *msi1_0 = buffer.data(msi1 + 0);
    const auto *msi1_1 = buffer.data(msi1 + 1);
    const auto *msi1_2 = buffer.data(msi1 + 2);
    const auto *msi1_3 = buffer.data(msi1 + 3);
    const auto *msi1_5 = buffer.data(msi1 + 5);
    const auto *msi1_6 = buffer.data(msi1 + 6);
    const auto *msi1_8 = buffer.data(msi1 + 8);
    const auto *msi1_9 = buffer.data(msi1 + 9);
    const auto *msi1_10 = buffer.data(msi1 + 10);
    const auto *msi1_12 = buffer.data(msi1 + 12);
    const auto *msi1_13 = buffer.data(msi1 + 13);
    const auto *msi1_14 = buffer.data(msi1 + 14);
    const auto *msi1_21 = buffer.data(msi1 + 21);
    const auto *msi1_23 = buffer.data(msi1 + 23);
    const auto *msi1_24 = buffer.data(msi1 + 24);
    const auto *msi1_25 = buffer.data(msi1 + 25);
    const auto *msi1_26 = buffer.data(msi1 + 26);
    const auto *msi1_27 = buffer.data(msi1 + 27);
    const auto *msi1_31 = buffer.data(msi1 + 31);
    const auto *msi1_34 = buffer.data(msi1 + 34);
    const auto *msi1_35 = buffer.data(msi1 + 35);
    const auto *msi1_38 = buffer.data(msi1 + 38);
    const auto *msi1_39 = buffer.data(msi1 + 39);
    const auto *msi1_40 = buffer.data(msi1 + 40);
    const auto *msi1_49 = buffer.data(msi1 + 49);
    const auto *msi1_50 = buffer.data(msi1 + 50);
    const auto *msi1_51 = buffer.data(msi1 + 51);
    const auto *msi1_52 = buffer.data(msi1 + 52);
    const auto *msi1_53 = buffer.data(msi1 + 53);
    const auto *msi1_58 = buffer.data(msi1 + 58);
    const auto *msi1_60 = buffer.data(msi1 + 60);
    const auto *msi1_61 = buffer.data(msi1 + 61);
    const auto *msi1_63 = buffer.data(msi1 + 63);
    const auto *msi1_64 = buffer.data(msi1 + 64);
    const auto *msi1_65 = buffer.data(msi1 + 65);
    const auto *msi1_67 = buffer.data(msi1 + 67);
    const auto *msi1_68 = buffer.data(msi1 + 68);
    const auto *msi1_69 = buffer.data(msi1 + 69);
    const auto *msi1_70 = buffer.data(msi1 + 70);
    const auto *msi1_78 = buffer.data(msi1 + 78);
    const auto *msi1_79 = buffer.data(msi1 + 79);

    const auto *msk_0 = buffer.data(msk + 0);
    const auto *msk_1 = buffer.data(msk + 1);
    const auto *msk_2 = buffer.data(msk + 2);
    const auto *msk_3 = buffer.data(msk + 3);
    const auto *msk_5 = buffer.data(msk + 5);
    const auto *msk_6 = buffer.data(msk + 6);
    const auto *msk_8 = buffer.data(msk + 8);
    const auto *msk_9 = buffer.data(msk + 9);
    const auto *msk_10 = buffer.data(msk + 10);
    const auto *msk_12 = buffer.data(msk + 12);
    const auto *msk_13 = buffer.data(msk + 13);
    const auto *msk_14 = buffer.data(msk + 14);
    const auto *msk_15 = buffer.data(msk + 15);
    const auto *msk_17 = buffer.data(msk + 17);
    const auto *msk_18 = buffer.data(msk + 18);
    const auto *msk_19 = buffer.data(msk + 19);
    const auto *msk_20 = buffer.data(msk + 20);
    const auto *msk_21 = buffer.data(msk + 21);
    const auto *msk_27 = buffer.data(msk + 27);
    const auto *msk_28 = buffer.data(msk + 28);
    const auto *msk_30 = buffer.data(msk + 30);
    const auto *msk_31 = buffer.data(msk + 31);
    const auto *msk_32 = buffer.data(msk + 32);
    const auto *msk_33 = buffer.data(msk + 33);
    const auto *msk_34 = buffer.data(msk + 34);
    const auto *msk_35 = buffer.data(msk + 35);
    const auto *msk_36 = buffer.data(msk + 36);
    const auto *msk_37 = buffer.data(msk + 37);
    const auto *msk_39 = buffer.data(msk + 39);
    const auto *msk_41 = buffer.data(msk + 41);
    const auto *msk_42 = buffer.data(msk + 42);
    const auto *msk_43 = buffer.data(msk + 43);
    const auto *msk_45 = buffer.data(msk + 45);
    const auto *msk_46 = buffer.data(msk + 46);
    const auto *msk_47 = buffer.data(msk + 47);
    const auto *msk_48 = buffer.data(msk + 48);
    const auto *msk_50 = buffer.data(msk + 50);
    const auto *msk_51 = buffer.data(msk + 51);
    const auto *msk_52 = buffer.data(msk + 52);
    const auto *msk_53 = buffer.data(msk + 53);
    const auto *msk_54 = buffer.data(msk + 54);
    const auto *msk_56 = buffer.data(msk + 56);
    const auto *msk_57 = buffer.data(msk + 57);
    const auto *msk_64 = buffer.data(msk + 64);
    const auto *msk_65 = buffer.data(msk + 65);
    const auto *msk_66 = buffer.data(msk + 66);
    const auto *msk_67 = buffer.data(msk + 67);
    const auto *msk_68 = buffer.data(msk + 68);
    const auto *msk_69 = buffer.data(msk + 69);
    const auto *msk_70 = buffer.data(msk + 70);
    const auto *msk_71 = buffer.data(msk + 71);
    const auto *msk_72 = buffer.data(msk + 72);
    const auto *msk_74 = buffer.data(msk + 74);
    const auto *msk_76 = buffer.data(msk + 76);
    const auto *msk_77 = buffer.data(msk + 77);
    const auto *msk_79 = buffer.data(msk + 79);
    const auto *msk_80 = buffer.data(msk + 80);
    const auto *msk_81 = buffer.data(msk + 81);
    const auto *msk_83 = buffer.data(msk + 83);
    const auto *msk_84 = buffer.data(msk + 84);
    const auto *msk_85 = buffer.data(msk + 85);
    const auto *msk_86 = buffer.data(msk + 86);
    const auto *msk_88 = buffer.data(msk + 88);
    const auto *msk_89 = buffer.data(msk + 89);
    const auto *msk_90 = buffer.data(msk + 90);
    const auto *msk_91 = buffer.data(msk + 91);
    const auto *msk_92 = buffer.data(msk + 92);
    const auto *msk_99 = buffer.data(msk + 99);
    const auto *msk_100 = buffer.data(msk + 100);
    const auto *msk_101 = buffer.data(msk + 101);
    const auto *msk_102 = buffer.data(msk + 102);
    const auto *msk_103 = buffer.data(msk + 103);
    const auto *msk_104 = buffer.data(msk + 104);
    const auto *msk_105 = buffer.data(msk + 105);
    const auto *msk_107 = buffer.data(msk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, lsk_0, msi0_0, \
                         msi1_0, msk_0, msk_1, msk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsk_0[k]
                 + f_1 * msi0_0[k]
                 - f_2 * msi1_0[k]
                 + f_3 * pc_x[k] * msk_0[k];

        t_1[k] = f_3 * pc_y[k] * msk_0[k];

        t_2[k] = f_3 * pc_z[k] * msk_0[k];

        t_3[k] = f_4 * msi0_0[k]
                 - f_5 * msi1_0[k]
                 + f_3 * pc_y[k] * msk_1[k];

        t_4[k] = f_3 * pc_y[k] * msk_2[k];

        t_5[k] = f_4 * msi0_0[k]
                 - f_5 * msi1_0[k]
                 + f_3 * pc_z[k] * msk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, msi0_1, msi0_2, msi0_3, msi1_1, \
                         msi1_2, msi1_3, msk_3, msk_5, msk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * msi0_1[k]
                 - f_7 * msi1_1[k]
                 + f_3 * pc_y[k] * msk_3[k];

        t_7[k] = f_3 * pc_z[k] * msk_3[k];

        t_8[k] = f_3 * pc_y[k] * msk_5[k];

        t_9[k] = f_6 * msi0_2[k]
                 - f_7 * msi1_2[k]
                 + f_3 * pc_z[k] * msk_5[k];

        t_10[k] = f_8 * msi0_3[k]
                  - f_9 * msi1_3[k]
                  + f_3 * pc_y[k] * msk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, msi0_5, msi0_6, \
                         msi1_5, msi1_6, msk_6, msk_8, msk_9, msk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * msk_6[k];

        t_12[k] = f_4 * msi0_5[k]
                  - f_5 * msi1_5[k]
                  + f_3 * pc_y[k] * msk_8[k];

        t_13[k] = f_3 * pc_y[k] * msk_9[k];

        t_14[k] = f_8 * msi0_5[k]
                  - f_9 * msi1_5[k]
                  + f_3 * pc_z[k] * msk_9[k];

        t_15[k] = f_10 * msi0_6[k]
                  - f_11 * msi1_6[k]
                  + f_3 * pc_y[k] * msk_10[k];

        t_16[k] = f_3 * pc_z[k] * msk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, msi0_8, msi0_9, msi1_8, msi1_9, \
                         msk_12, msk_13, msk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * msi0_8[k]
                  - f_7 * msi1_8[k]
                  + f_3 * pc_y[k] * msk_12[k];

        t_18[k] = f_4 * msi0_9[k]
                  - f_5 * msi1_9[k]
                  + f_3 * pc_y[k] * msk_13[k];

        t_19[k] = f_3 * pc_y[k] * msk_14[k];

        t_20[k] = f_10 * msi0_9[k]
                  - f_11 * msi1_9[k]
                  + f_3 * pc_z[k] * msk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, msi0_10, msi0_12, msi0_13, \
                         msi1_10, msi1_12, msi1_13, msk_15, msk_17, \
                         msk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * msi0_10[k]
                  - f_13 * msi1_10[k]
                  + f_3 * pc_y[k] * msk_15[k];

        t_22[k] = f_3 * pc_z[k] * msk_15[k];

        t_23[k] = f_8 * msi0_12[k]
                  - f_9 * msi1_12[k]
                  + f_3 * pc_y[k] * msk_17[k];

        t_24[k] = f_6 * msi0_13[k]
                  - f_7 * msi1_13[k]
                  + f_3 * pc_y[k] * msk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, lsk_28, msi0_14, \
                         msi1_14, msk_19, msk_20, msk_21, msk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * msi0_14[k]
                  - f_5 * msi1_14[k]
                  + f_3 * pc_y[k] * msk_19[k];

        t_26[k] = f_3 * pc_y[k] * msk_20[k];

        t_27[k] = f_12 * msi0_14[k]
                  - f_13 * msi1_14[k]
                  + f_3 * pc_z[k] * msk_20[k];

        t_28[k] = f_0 * lsk_28[k]
                  + f_3 * pc_x[k] * msk_28[k];

        t_29[k] = f_3 * pc_z[k] * msk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, lsk_30, lsk_31, lsk_32, \
                         lsk_33, msk_27, msk_30, msk_31, msk_32, \
                         msk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * lsk_30[k]
                  + f_3 * pc_x[k] * msk_30[k];

        t_31[k] = f_0 * lsk_31[k]
                  + f_3 * pc_x[k] * msk_31[k];

        t_32[k] = f_0 * lsk_32[k]
                  + f_3 * pc_x[k] * msk_32[k];

        t_33[k] = f_0 * lsk_33[k]
                  + f_3 * pc_x[k] * msk_33[k];

        t_34[k] = f_3 * pc_y[k] * msk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, lsk_35, msi0_21, msi0_23, \
                         msi1_21, msi1_23, msk_28, msk_30, msk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * lsk_35[k]
                  + f_3 * pc_x[k] * msk_35[k];

        t_36[k] = f_1 * msi0_21[k]
                  - f_2 * msi1_21[k]
                  + f_3 * pc_y[k] * msk_28[k];

        t_37[k] = f_3 * pc_z[k] * msk_28[k];

        t_38[k] = f_12 * msi0_23[k]
                  - f_13 * msi1_23[k]
                  + f_3 * pc_y[k] * msk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, msi0_24, msi0_25, msi0_26, msi1_24, msi1_25, \
                         msi1_26, msk_31, msk_32, msk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * msi0_24[k]
                  - f_11 * msi1_24[k]
                  + f_3 * pc_y[k] * msk_31[k];

        t_40[k] = f_8 * msi0_25[k]
                  - f_9 * msi1_25[k]
                  + f_3 * pc_y[k] * msk_32[k];

        t_41[k] = f_6 * msi0_26[k]
                  - f_7 * msi1_26[k]
                  + f_3 * pc_y[k] * msk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, lsl0_0, lsk_0, \
                         lsl1_0, msi0_27, msi1_27, msk_34, msk_35, \
                         msk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * msi0_27[k]
                  - f_5 * msi1_27[k]
                  + f_3 * pc_y[k] * msk_34[k];

        t_43[k] = f_3 * pc_y[k] * msk_35[k];

        t_44[k] = f_1 * msi0_27[k]
                  - f_2 * msi1_27[k]
                  + f_3 * pc_z[k] * msk_35[k];

        t_45[k] = pa_y[k] * lsl0_0[k]
                  - f_14 * pc_y[k] * lsl1_0[k];

        t_46[k] = f_15 * lsk_0[k]
                  + f_3 * pc_y[k] * msk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, lsl0_3, lsl0_5, lsk_1, \
                         lsl1_3, lsl1_5, msk_36, msk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * msk_36[k];

        t_48[k] = pa_y[k] * lsl0_3[k]
                  + f_16 * lsk_1[k]
                  - f_14 * pc_y[k] * lsl1_3[k];

        t_49[k] = f_3 * pc_z[k] * msk_37[k];

        t_50[k] = pa_y[k] * lsl0_5[k]
                  - f_14 * pc_y[k] * lsl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, lsl0_6, lsl0_9, lsk_3, \
                         lsk_5, lsl1_6, lsl1_9, msk_39, msk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * lsl0_6[k]
                  + f_17 * lsk_3[k]
                  - f_14 * pc_y[k] * lsl1_6[k];

        t_52[k] = f_3 * pc_z[k] * msk_39[k];

        t_53[k] = f_15 * lsk_5[k]
                  + f_3 * pc_y[k] * msk_41[k];

        t_54[k] = pa_y[k] * lsl0_9[k]
                  - f_14 * pc_y[k] * lsl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, lsl0_10, lsk_6, lsk_9, \
                         lsl1_10, msi0_31, msi1_31, msk_42, msk_43, \
                         msk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * lsl0_10[k]
                  + f_18 * lsk_6[k]
                  - f_14 * pc_y[k] * lsl1_10[k];

        t_56[k] = f_3 * pc_z[k] * msk_42[k];

        t_57[k] = f_4 * msi0_31[k]
                  - f_5 * msi1_31[k]
                  + f_3 * pc_z[k] * msk_43[k];

        t_58[k] = f_15 * lsk_9[k]
                  + f_3 * pc_y[k] * msk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, lsl0_14, lsl0_15, lsk_10, \
                         lsl1_14, lsl1_15, msi0_34, msi1_34, msk_46, \
                         msk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * lsl0_14[k]
                  - f_14 * pc_y[k] * lsl1_14[k];

        t_60[k] = pa_y[k] * lsl0_15[k]
                  + f_19 * lsk_10[k]
                  - f_14 * pc_y[k] * lsl1_15[k];

        t_61[k] = f_3 * pc_z[k] * msk_46[k];

        t_62[k] = f_4 * msi0_34[k]
                  - f_5 * msi1_34[k]
                  + f_3 * pc_z[k] * msk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, lsl0_20, lsk_14, lsl1_20, \
                         msi0_35, msi1_35, msk_48, msk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * msi0_35[k]
                  - f_7 * msi1_35[k]
                  + f_3 * pc_z[k] * msk_48[k];

        t_64[k] = f_15 * lsk_14[k]
                  + f_3 * pc_y[k] * msk_50[k];

        t_65[k] = pa_y[k] * lsl0_20[k]
                  - f_14 * pc_y[k] * lsl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, lsl0_21, lsk_15, lsl1_21, \
                         msi0_38, msi1_38, msk_51, msk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * lsl0_21[k]
                  + f_20 * lsk_15[k]
                  - f_14 * pc_y[k] * lsl1_21[k];

        t_67[k] = f_3 * pc_z[k] * msk_51[k];

        t_68[k] = f_4 * msi0_38[k]
                  - f_5 * msi1_38[k]
                  + f_3 * pc_z[k] * msk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, lsk_20, msi0_39, msi0_40, msi1_39, \
                         msi1_40, msk_53, msk_54, msk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * msi0_39[k]
                  - f_7 * msi1_39[k]
                  + f_3 * pc_z[k] * msk_53[k];

        t_70[k] = f_8 * msi0_40[k]
                  - f_9 * msi1_40[k]
                  + f_3 * pc_z[k] * msk_54[k];

        t_71[k] = f_15 * lsk_20[k]
                  + f_3 * pc_y[k] * msk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, lsl0_27, lsk_64, \
                         lsk_66, lsl1_27, msk_57, msk_64, msk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * lsl0_27[k]
                  - f_14 * pc_y[k] * lsl1_27[k];

        t_73[k] = f_21 * lsk_64[k]
                  + f_3 * pc_x[k] * msk_64[k];

        t_74[k] = f_3 * pc_z[k] * msk_57[k];

        t_75[k] = f_21 * lsk_66[k]
                  + f_3 * pc_x[k] * msk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, lsk_67, lsk_68, lsk_69, lsk_70, \
                         lsk_71, msk_67, msk_68, msk_69, msk_70, \
                         msk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_21 * lsk_67[k]
                  + f_3 * pc_x[k] * msk_67[k];

        t_77[k] = f_21 * lsk_68[k]
                  + f_3 * pc_x[k] * msk_68[k];

        t_78[k] = f_21 * lsk_69[k]
                  + f_3 * pc_x[k] * msk_69[k];

        t_79[k] = f_21 * lsk_70[k]
                  + f_3 * pc_x[k] * msk_70[k];

        t_80[k] = f_21 * lsk_71[k]
                  + f_3 * pc_x[k] * msk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, lsk_28, msi0_49, msi0_50, \
                         msi1_49, msi1_50, msk_64, msk_65, msk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * lsk_28[k]
                  + f_1 * msi0_49[k]
                  - f_2 * msi1_49[k]
                  + f_3 * pc_y[k] * msk_64[k];

        t_82[k] = f_3 * pc_z[k] * msk_64[k];

        t_83[k] = f_4 * msi0_49[k]
                  - f_5 * msi1_49[k]
                  + f_3 * pc_z[k] * msk_65[k];

        t_84[k] = f_6 * msi0_50[k]
                  - f_7 * msi1_50[k]
                  + f_3 * pc_z[k] * msk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, msi0_51, msi0_52, msi0_53, msi1_51, msi1_52, \
                         msi1_53, msk_67, msk_68, msk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * msi0_51[k]
                  - f_9 * msi1_51[k]
                  + f_3 * pc_z[k] * msk_67[k];

        t_86[k] = f_10 * msi0_52[k]
                  - f_11 * msi1_52[k]
                  + f_3 * pc_z[k] * msk_68[k];

        t_87[k] = f_12 * msi0_53[k]
                  - f_13 * msi1_53[k]
                  + f_3 * pc_z[k] * msk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, lsl0_0, lsl0_44, \
                         lsk_35, lsl1_0, lsl1_44, msk_71, msk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * lsk_35[k]
                  + f_3 * pc_y[k] * msk_71[k];

        t_89[k] = pa_y[k] * lsl0_44[k]
                  - f_14 * pc_y[k] * lsl1_44[k];

        t_90[k] = pa_z[k] * lsl0_0[k]
                  - f_14 * pc_z[k] * lsl1_0[k];

        t_91[k] = f_3 * pc_y[k] * msk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, lsl0_3, lsl0_5, lsk_0, \
                         lsk_2, lsl1_3, lsl1_5, msk_72, msk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * lsk_0[k]
                  + f_3 * pc_z[k] * msk_72[k];

        t_93[k] = pa_z[k] * lsl0_3[k]
                  - f_14 * pc_z[k] * lsl1_3[k];

        t_94[k] = f_3 * pc_y[k] * msk_74[k];

        t_95[k] = pa_z[k] * lsl0_5[k]
                  + f_16 * lsk_2[k]
                  - f_14 * pc_z[k] * lsl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, lsl0_6, lsl0_9, lsk_5, \
                         lsl1_6, lsl1_9, msi0_58, msi1_58, msk_76, \
                         msk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * lsl0_6[k]
                  - f_14 * pc_z[k] * lsl1_6[k];

        t_97[k] = f_4 * msi0_58[k]
                  - f_5 * msi1_58[k]
                  + f_3 * pc_y[k] * msk_76[k];

        t_98[k] = f_3 * pc_y[k] * msk_77[k];

        t_99[k] = pa_z[k] * lsl0_9[k]
                  + f_17 * lsk_5[k]
                  - f_14 * pc_z[k] * lsl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, lsl0_10, lsl1_10, \
                         msi0_60, msi0_61, msi1_60, msi1_61, msk_79, msk_80, \
                         msk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * lsl0_10[k]
                   - f_14 * pc_z[k] * lsl1_10[k];

        t_101[k] = f_6 * msi0_60[k]
                   - f_7 * msi1_60[k]
                   + f_3 * pc_y[k] * msk_79[k];

        t_102[k] = f_4 * msi0_61[k]
                   - f_5 * msi1_61[k]
                   + f_3 * pc_y[k] * msk_80[k];

        t_103[k] = f_3 * pc_y[k] * msk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, lsl0_14, lsl0_15, lsk_9, \
                         lsl1_14, lsl1_15, msi0_63, msi1_63, msk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * lsl0_14[k]
                   + f_18 * lsk_9[k]
                   - f_14 * pc_z[k] * lsl1_14[k];

        t_105[k] = pa_z[k] * lsl0_15[k]
                   - f_14 * pc_z[k] * lsl1_15[k];

        t_106[k] = f_8 * msi0_63[k]
                   - f_9 * msi1_63[k]
                   + f_3 * pc_y[k] * msk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, msi0_64, msi0_65, msi1_64, msi1_65, \
                         msk_84, msk_85, msk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * msi0_64[k]
                   - f_7 * msi1_64[k]
                   + f_3 * pc_y[k] * msk_84[k];

        t_108[k] = f_4 * msi0_65[k]
                   - f_5 * msi1_65[k]
                   + f_3 * pc_y[k] * msk_85[k];

        t_109[k] = f_3 * pc_y[k] * msk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, lsl0_20, lsl0_21, lsk_14, \
                         lsl1_20, lsl1_21, msi0_67, msi1_67, msk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * lsl0_20[k]
                   + f_19 * lsk_14[k]
                   - f_14 * pc_z[k] * lsl1_20[k];

        t_111[k] = pa_z[k] * lsl0_21[k]
                   - f_14 * pc_z[k] * lsl1_21[k];

        t_112[k] = f_10 * msi0_67[k]
                   - f_11 * msi1_67[k]
                   + f_3 * pc_y[k] * msk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, msi0_68, msi0_69, msi0_70, msi1_68, \
                         msi1_69, msi1_70, msk_89, msk_90, msk_91, \
                         msk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * msi0_68[k]
                   - f_9 * msi1_68[k]
                   + f_3 * pc_y[k] * msk_89[k];

        t_114[k] = f_6 * msi0_69[k]
                   - f_7 * msi1_69[k]
                   + f_3 * pc_y[k] * msk_90[k];

        t_115[k] = f_4 * msi0_70[k]
                   - f_5 * msi1_70[k]
                   + f_3 * pc_y[k] * msk_91[k];

        t_116[k] = f_3 * pc_y[k] * msk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, lsl0_27, lsk_20, \
                         lsk_100, lsk_101, lsk_102, lsl1_27, msk_100, msk_101, \
                         msk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * lsl0_27[k]
                   + f_20 * lsk_20[k]
                   - f_14 * pc_z[k] * lsl1_27[k];

        t_118[k] = f_21 * lsk_100[k]
                   + f_3 * pc_x[k] * msk_100[k];

        t_119[k] = f_21 * lsk_101[k]
                   + f_3 * pc_x[k] * msk_101[k];

        t_120[k] = f_21 * lsk_102[k]
                   + f_3 * pc_x[k] * msk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, lsk_103, lsk_104, \
                         lsk_105, lsk_107, msk_99, msk_103, msk_104, msk_105, \
                         msk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_21 * lsk_103[k]
                   + f_3 * pc_x[k] * msk_103[k];

        t_122[k] = f_21 * lsk_104[k]
                   + f_3 * pc_x[k] * msk_104[k];

        t_123[k] = f_21 * lsk_105[k]
                   + f_3 * pc_x[k] * msk_105[k];

        t_124[k] = f_3 * pc_y[k] * msk_99[k];

        t_125[k] = f_21 * lsk_107[k]
                   + f_3 * pc_x[k] * msk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, lsl0_36, lsl1_36, msi0_78, \
                         msi0_79, msi1_78, msi1_79, msk_101, msk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * lsl0_36[k]
                   - f_14 * pc_z[k] * lsl1_36[k];

        t_127[k] = f_22 * msi0_78[k]
                   - f_23 * msi1_78[k]
                   + f_3 * pc_y[k] * msk_101[k];

        t_128[k] = f_12 * msi0_79[k]
                   - f_13 * msi1_79[k]
                   + f_3 * pc_y[k] * msk_102[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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
    const auto f_24 = 3.5 / q;

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

    const auto *lsl0_48 = buffer.data(lsl0 + 48);
    const auto *lsl0_51 = buffer.data(lsl0 + 51);
    const auto *lsl0_55 = buffer.data(lsl0 + 55);
    const auto *lsl0_60 = buffer.data(lsl0 + 60);
    const auto *lsl0_66 = buffer.data(lsl0 + 66);
    const auto *lsl0_81 = buffer.data(lsl0 + 81);
    const auto *lsl0_90 = buffer.data(lsl0 + 90);
    const auto *lsl0_95 = buffer.data(lsl0 + 95);
    const auto *lsl0_99 = buffer.data(lsl0 + 99);
    const auto *lsl0_102 = buffer.data(lsl0 + 102);
    const auto *lsl0_104 = buffer.data(lsl0 + 104);
    const auto *lsl0_107 = buffer.data(lsl0 + 107);
    const auto *lsl0_108 = buffer.data(lsl0 + 108);
    const auto *lsl0_110 = buffer.data(lsl0 + 110);
    const auto *lsl0_113 = buffer.data(lsl0 + 113);
    const auto *lsl0_114 = buffer.data(lsl0 + 114);
    const auto *lsl0_115 = buffer.data(lsl0 + 115);
    const auto *lsl0_117 = buffer.data(lsl0 + 117);
    const auto *lsl0_134 = buffer.data(lsl0 + 134);

    const auto *lsk_35 = buffer.data(lsk + 35);
    const auto *lsk_36 = buffer.data(lsk + 36);
    const auto *lsk_39 = buffer.data(lsk + 39);
    const auto *lsk_41 = buffer.data(lsk + 41);
    const auto *lsk_42 = buffer.data(lsk + 42);
    const auto *lsk_45 = buffer.data(lsk + 45);
    const auto *lsk_46 = buffer.data(lsk + 46);
    const auto *lsk_50 = buffer.data(lsk + 50);
    const auto *lsk_51 = buffer.data(lsk + 51);
    const auto *lsk_56 = buffer.data(lsk + 56);
    const auto *lsk_64 = buffer.data(lsk + 64);
    const auto *lsk_71 = buffer.data(lsk + 71);
    const auto *lsk_72 = buffer.data(lsk + 72);
    const auto *lsk_74 = buffer.data(lsk + 74);
    const auto *lsk_77 = buffer.data(lsk + 77);
    const auto *lsk_80 = buffer.data(lsk + 80);
    const auto *lsk_81 = buffer.data(lsk + 81);
    const auto *lsk_84 = buffer.data(lsk + 84);
    const auto *lsk_85 = buffer.data(lsk + 85);
    const auto *lsk_86 = buffer.data(lsk + 86);
    const auto *lsk_89 = buffer.data(lsk + 89);
    const auto *lsk_90 = buffer.data(lsk + 90);
    const auto *lsk_91 = buffer.data(lsk + 91);
    const auto *lsk_92 = buffer.data(lsk + 92);
    const auto *lsk_102 = buffer.data(lsk + 102);
    const auto *lsk_103 = buffer.data(lsk + 103);
    const auto *lsk_104 = buffer.data(lsk + 104);
    const auto *lsk_105 = buffer.data(lsk + 105);
    const auto *lsk_106 = buffer.data(lsk + 106);
    const auto *lsk_107 = buffer.data(lsk + 107);
    const auto *lsk_108 = buffer.data(lsk + 108);
    const auto *lsk_111 = buffer.data(lsk + 111);
    const auto *lsk_114 = buffer.data(lsk + 114);
    const auto *lsk_118 = buffer.data(lsk + 118);
    const auto *lsk_123 = buffer.data(lsk + 123);
    const auto *lsk_129 = buffer.data(lsk + 129);
    const auto *lsk_136 = buffer.data(lsk + 136);
    const auto *lsk_138 = buffer.data(lsk + 138);
    const auto *lsk_139 = buffer.data(lsk + 139);
    const auto *lsk_140 = buffer.data(lsk + 140);
    const auto *lsk_141 = buffer.data(lsk + 141);
    const auto *lsk_142 = buffer.data(lsk + 142);
    const auto *lsk_143 = buffer.data(lsk + 143);
    const auto *lsk_172 = buffer.data(lsk + 172);
    const auto *lsk_173 = buffer.data(lsk + 173);
    const auto *lsk_174 = buffer.data(lsk + 174);
    const auto *lsk_175 = buffer.data(lsk + 175);
    const auto *lsk_176 = buffer.data(lsk + 176);
    const auto *lsk_177 = buffer.data(lsk + 177);
    const auto *lsk_178 = buffer.data(lsk + 178);
    const auto *lsk_179 = buffer.data(lsk + 179);
    const auto *lsk_180 = buffer.data(lsk + 180);
    const auto *lsk_185 = buffer.data(lsk + 185);
    const auto *lsk_189 = buffer.data(lsk + 189);
    const auto *lsk_194 = buffer.data(lsk + 194);
    const auto *lsk_200 = buffer.data(lsk + 200);

    const auto *lsl1_48 = buffer.data(lsl1 + 48);
    const auto *lsl1_51 = buffer.data(lsl1 + 51);
    const auto *lsl1_55 = buffer.data(lsl1 + 55);
    const auto *lsl1_60 = buffer.data(lsl1 + 60);
    const auto *lsl1_66 = buffer.data(lsl1 + 66);
    const auto *lsl1_81 = buffer.data(lsl1 + 81);
    const auto *lsl1_90 = buffer.data(lsl1 + 90);
    const auto *lsl1_95 = buffer.data(lsl1 + 95);
    const auto *lsl1_99 = buffer.data(lsl1 + 99);
    const auto *lsl1_102 = buffer.data(lsl1 + 102);
    const auto *lsl1_104 = buffer.data(lsl1 + 104);
    const auto *lsl1_107 = buffer.data(lsl1 + 107);
    const auto *lsl1_108 = buffer.data(lsl1 + 108);
    const auto *lsl1_110 = buffer.data(lsl1 + 110);
    const auto *lsl1_113 = buffer.data(lsl1 + 113);
    const auto *lsl1_114 = buffer.data(lsl1 + 114);
    const auto *lsl1_115 = buffer.data(lsl1 + 115);
    const auto *lsl1_117 = buffer.data(lsl1 + 117);
    const auto *lsl1_134 = buffer.data(lsl1 + 134);

    const auto *msi0_80 = buffer.data(msi0 + 80);
    const auto *msi0_81 = buffer.data(msi0 + 81);
    const auto *msi0_82 = buffer.data(msi0 + 82);
    const auto *msi0_83 = buffer.data(msi0 + 83);
    const auto *msi0_84 = buffer.data(msi0 + 84);
    const auto *msi0_86 = buffer.data(msi0 + 86);
    const auto *msi0_87 = buffer.data(msi0 + 87);
    const auto *msi0_89 = buffer.data(msi0 + 89);
    const auto *msi0_90 = buffer.data(msi0 + 90);
    const auto *msi0_91 = buffer.data(msi0 + 91);
    const auto *msi0_93 = buffer.data(msi0 + 93);
    const auto *msi0_94 = buffer.data(msi0 + 94);
    const auto *msi0_95 = buffer.data(msi0 + 95);
    const auto *msi0_96 = buffer.data(msi0 + 96);
    const auto *msi0_98 = buffer.data(msi0 + 98);
    const auto *msi0_99 = buffer.data(msi0 + 99);
    const auto *msi0_105 = buffer.data(msi0 + 105);
    const auto *msi0_106 = buffer.data(msi0 + 106);
    const auto *msi0_107 = buffer.data(msi0 + 107);
    const auto *msi0_108 = buffer.data(msi0 + 108);
    const auto *msi0_109 = buffer.data(msi0 + 109);
    const auto *msi0_111 = buffer.data(msi0 + 111);
    const auto *msi0_135 = buffer.data(msi0 + 135);
    const auto *msi0_136 = buffer.data(msi0 + 136);
    const auto *msi0_137 = buffer.data(msi0 + 137);
    const auto *msi0_138 = buffer.data(msi0 + 138);
    const auto *msi0_139 = buffer.data(msi0 + 139);
    const auto *msi0_140 = buffer.data(msi0 + 140);
    const auto *msi0_141 = buffer.data(msi0 + 141);
    const auto *msi0_142 = buffer.data(msi0 + 142);
    const auto *msi0_143 = buffer.data(msi0 + 143);
    const auto *msi0_144 = buffer.data(msi0 + 144);
    const auto *msi0_145 = buffer.data(msi0 + 145);
    const auto *msi0_146 = buffer.data(msi0 + 146);
    const auto *msi0_147 = buffer.data(msi0 + 147);
    const auto *msi0_148 = buffer.data(msi0 + 148);
    const auto *msi0_149 = buffer.data(msi0 + 149);
    const auto *msi0_154 = buffer.data(msi0 + 154);
    const auto *msi0_160 = buffer.data(msi0 + 160);

    const auto *msi1_80 = buffer.data(msi1 + 80);
    const auto *msi1_81 = buffer.data(msi1 + 81);
    const auto *msi1_82 = buffer.data(msi1 + 82);
    const auto *msi1_83 = buffer.data(msi1 + 83);
    const auto *msi1_84 = buffer.data(msi1 + 84);
    const auto *msi1_86 = buffer.data(msi1 + 86);
    const auto *msi1_87 = buffer.data(msi1 + 87);
    const auto *msi1_89 = buffer.data(msi1 + 89);
    const auto *msi1_90 = buffer.data(msi1 + 90);
    const auto *msi1_91 = buffer.data(msi1 + 91);
    const auto *msi1_93 = buffer.data(msi1 + 93);
    const auto *msi1_94 = buffer.data(msi1 + 94);
    const auto *msi1_95 = buffer.data(msi1 + 95);
    const auto *msi1_96 = buffer.data(msi1 + 96);
    const auto *msi1_98 = buffer.data(msi1 + 98);
    const auto *msi1_99 = buffer.data(msi1 + 99);
    const auto *msi1_105 = buffer.data(msi1 + 105);
    const auto *msi1_106 = buffer.data(msi1 + 106);
    const auto *msi1_107 = buffer.data(msi1 + 107);
    const auto *msi1_108 = buffer.data(msi1 + 108);
    const auto *msi1_109 = buffer.data(msi1 + 109);
    const auto *msi1_111 = buffer.data(msi1 + 111);
    const auto *msi1_135 = buffer.data(msi1 + 135);
    const auto *msi1_136 = buffer.data(msi1 + 136);
    const auto *msi1_137 = buffer.data(msi1 + 137);
    const auto *msi1_138 = buffer.data(msi1 + 138);
    const auto *msi1_139 = buffer.data(msi1 + 139);
    const auto *msi1_140 = buffer.data(msi1 + 140);
    const auto *msi1_141 = buffer.data(msi1 + 141);
    const auto *msi1_142 = buffer.data(msi1 + 142);
    const auto *msi1_143 = buffer.data(msi1 + 143);
    const auto *msi1_144 = buffer.data(msi1 + 144);
    const auto *msi1_145 = buffer.data(msi1 + 145);
    const auto *msi1_146 = buffer.data(msi1 + 146);
    const auto *msi1_147 = buffer.data(msi1 + 147);
    const auto *msi1_148 = buffer.data(msi1 + 148);
    const auto *msi1_149 = buffer.data(msi1 + 149);
    const auto *msi1_154 = buffer.data(msi1 + 154);
    const auto *msi1_160 = buffer.data(msi1 + 160);

    const auto *msk_103 = buffer.data(msk + 103);
    const auto *msk_104 = buffer.data(msk + 104);
    const auto *msk_105 = buffer.data(msk + 105);
    const auto *msk_106 = buffer.data(msk + 106);
    const auto *msk_107 = buffer.data(msk + 107);
    const auto *msk_108 = buffer.data(msk + 108);
    const auto *msk_109 = buffer.data(msk + 109);
    const auto *msk_110 = buffer.data(msk + 110);
    const auto *msk_111 = buffer.data(msk + 111);
    const auto *msk_113 = buffer.data(msk + 113);
    const auto *msk_114 = buffer.data(msk + 114);
    const auto *msk_115 = buffer.data(msk + 115);
    const auto *msk_117 = buffer.data(msk + 117);
    const auto *msk_118 = buffer.data(msk + 118);
    const auto *msk_119 = buffer.data(msk + 119);
    const auto *msk_120 = buffer.data(msk + 120);
    const auto *msk_122 = buffer.data(msk + 122);
    const auto *msk_123 = buffer.data(msk + 123);
    const auto *msk_124 = buffer.data(msk + 124);
    const auto *msk_125 = buffer.data(msk + 125);
    const auto *msk_126 = buffer.data(msk + 126);
    const auto *msk_128 = buffer.data(msk + 128);
    const auto *msk_129 = buffer.data(msk + 129);
    const auto *msk_136 = buffer.data(msk + 136);
    const auto *msk_137 = buffer.data(msk + 137);
    const auto *msk_138 = buffer.data(msk + 138);
    const auto *msk_139 = buffer.data(msk + 139);
    const auto *msk_140 = buffer.data(msk + 140);
    const auto *msk_141 = buffer.data(msk + 141);
    const auto *msk_142 = buffer.data(msk + 142);
    const auto *msk_143 = buffer.data(msk + 143);
    const auto *msk_144 = buffer.data(msk + 144);
    const auto *msk_146 = buffer.data(msk + 146);
    const auto *msk_147 = buffer.data(msk + 147);
    const auto *msk_149 = buffer.data(msk + 149);
    const auto *msk_150 = buffer.data(msk + 150);
    const auto *msk_153 = buffer.data(msk + 153);
    const auto *msk_154 = buffer.data(msk + 154);
    const auto *msk_158 = buffer.data(msk + 158);
    const auto *msk_159 = buffer.data(msk + 159);
    const auto *msk_164 = buffer.data(msk + 164);
    const auto *msk_172 = buffer.data(msk + 172);
    const auto *msk_173 = buffer.data(msk + 173);
    const auto *msk_174 = buffer.data(msk + 174);
    const auto *msk_175 = buffer.data(msk + 175);
    const auto *msk_176 = buffer.data(msk + 176);
    const auto *msk_177 = buffer.data(msk + 177);
    const auto *msk_178 = buffer.data(msk + 178);
    const auto *msk_179 = buffer.data(msk + 179);
    const auto *msk_180 = buffer.data(msk + 180);
    const auto *msk_181 = buffer.data(msk + 181);
    const auto *msk_182 = buffer.data(msk + 182);
    const auto *msk_183 = buffer.data(msk + 183);
    const auto *msk_184 = buffer.data(msk + 184);
    const auto *msk_185 = buffer.data(msk + 185);
    const auto *msk_186 = buffer.data(msk + 186);
    const auto *msk_187 = buffer.data(msk + 187);
    const auto *msk_188 = buffer.data(msk + 188);
    const auto *msk_189 = buffer.data(msk + 189);
    const auto *msk_190 = buffer.data(msk + 190);
    const auto *msk_191 = buffer.data(msk + 191);
    const auto *msk_192 = buffer.data(msk + 192);
    const auto *msk_193 = buffer.data(msk + 193);
    const auto *msk_194 = buffer.data(msk + 194);
    const auto *msk_200 = buffer.data(msk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, msi0_80, msi0_81, msi0_82, msi1_80, \
                         msi1_81, msi1_82, msk_103, msk_104, msk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * msi0_80[k]
                   - f_11 * msi1_80[k]
                   + f_3 * pc_y[k] * msk_103[k];

        t_130[k] = f_8 * msi0_81[k]
                   - f_9 * msi1_81[k]
                   + f_3 * pc_y[k] * msk_104[k];

        t_131[k] = f_6 * msi0_82[k]
                   - f_7 * msi1_82[k]
                   + f_3 * pc_y[k] * msk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, lsk_35, lsk_108, \
                         msi0_83, msi0_84, msi1_83, msi1_84, msk_106, msk_107, \
                         msk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * msi0_83[k]
                   - f_5 * msi1_83[k]
                   + f_3 * pc_y[k] * msk_106[k];

        t_133[k] = f_3 * pc_y[k] * msk_107[k];

        t_134[k] = f_15 * lsk_35[k]
                   + f_1 * msi0_83[k]
                   - f_2 * msi1_83[k]
                   + f_3 * pc_z[k] * msk_107[k];

        t_135[k] = f_24 * lsk_108[k]
                   + f_1 * msi0_84[k]
                   - f_2 * msi1_84[k]
                   + f_3 * pc_x[k] * msk_108[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, lsk_36, lsk_111, \
                         msi0_87, msi1_87, msk_108, msk_109, msk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * lsk_36[k]
                   + f_3 * pc_y[k] * msk_108[k];

        t_137[k] = f_3 * pc_z[k] * msk_108[k];

        t_138[k] = f_24 * lsk_111[k]
                   + f_12 * msi0_87[k]
                   - f_13 * msi1_87[k]
                   + f_3 * pc_x[k] * msk_111[k];

        t_139[k] = f_3 * pc_z[k] * msk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_z, lsk_114, msi0_84, msi0_90, msi1_84, \
                         msi1_90, msk_110, msk_111, msk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * msi0_84[k]
                   - f_5 * msi1_84[k]
                   + f_3 * pc_z[k] * msk_110[k];

        t_141[k] = f_24 * lsk_114[k]
                   + f_10 * msi0_90[k]
                   - f_11 * msi1_90[k]
                   + f_3 * pc_x[k] * msk_114[k];

        t_142[k] = f_3 * pc_z[k] * msk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, lsk_41, lsk_118, \
                         msi0_86, msi0_94, msi1_86, msi1_94, msk_113, msk_114, \
                         msk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * lsk_41[k]
                   + f_3 * pc_y[k] * msk_113[k];

        t_144[k] = f_6 * msi0_86[k]
                   - f_7 * msi1_86[k]
                   + f_3 * pc_z[k] * msk_113[k];

        t_145[k] = f_24 * lsk_118[k]
                   + f_8 * msi0_94[k]
                   - f_9 * msi1_94[k]
                   + f_3 * pc_x[k] * msk_118[k];

        t_146[k] = f_3 * pc_z[k] * msk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, lsk_45, msi0_87, msi0_89, msi1_87, \
                         msi1_89, msk_115, msk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * msi0_87[k]
                   - f_5 * msi1_87[k]
                   + f_3 * pc_z[k] * msk_115[k];

        t_148[k] = f_16 * lsk_45[k]
                   + f_3 * pc_y[k] * msk_117[k];

        t_149[k] = f_8 * msi0_89[k]
                   - f_9 * msi1_89[k]
                   + f_3 * pc_z[k] * msk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, lsk_123, msi0_90, msi0_99, msi1_90, \
                         msi1_99, msk_118, msk_119, msk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_24 * lsk_123[k]
                   + f_6 * msi0_99[k]
                   - f_7 * msi1_99[k]
                   + f_3 * pc_x[k] * msk_123[k];

        t_151[k] = f_3 * pc_z[k] * msk_118[k];

        t_152[k] = f_4 * msi0_90[k]
                   - f_5 * msi1_90[k]
                   + f_3 * pc_z[k] * msk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, lsk_50, msi0_91, msi0_93, msi1_91, \
                         msi1_93, msk_120, msk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * msi0_91[k]
                   - f_7 * msi1_91[k]
                   + f_3 * pc_z[k] * msk_120[k];

        t_154[k] = f_16 * lsk_50[k]
                   + f_3 * pc_y[k] * msk_122[k];

        t_155[k] = f_10 * msi0_93[k]
                   - f_11 * msi1_93[k]
                   + f_3 * pc_z[k] * msk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, lsk_129, msi0_94, msi0_105, msi1_94, \
                         msi1_105, msk_123, msk_124, msk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_24 * lsk_129[k]
                   + f_4 * msi0_105[k]
                   - f_5 * msi1_105[k]
                   + f_3 * pc_x[k] * msk_129[k];

        t_157[k] = f_3 * pc_z[k] * msk_123[k];

        t_158[k] = f_4 * msi0_94[k]
                   - f_5 * msi1_94[k]
                   + f_3 * pc_z[k] * msk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, lsk_56, msi0_95, msi0_96, \
                         msi0_98, msi1_95, msi1_96, msi1_98, msk_125, msk_126, \
                         msk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * msi0_95[k]
                   - f_7 * msi1_95[k]
                   + f_3 * pc_z[k] * msk_125[k];

        t_160[k] = f_8 * msi0_96[k]
                   - f_9 * msi1_96[k]
                   + f_3 * pc_z[k] * msk_126[k];

        t_161[k] = f_16 * lsk_56[k]
                   + f_3 * pc_y[k] * msk_128[k];

        t_162[k] = f_12 * msi0_98[k]
                   - f_13 * msi1_98[k]
                   + f_3 * pc_z[k] * msk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, lsk_136, lsk_138, \
                         lsk_139, lsk_140, msk_129, msk_136, msk_138, msk_139, \
                         msk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_24 * lsk_136[k]
                   + f_3 * pc_x[k] * msk_136[k];

        t_164[k] = f_3 * pc_z[k] * msk_129[k];

        t_165[k] = f_24 * lsk_138[k]
                   + f_3 * pc_x[k] * msk_138[k];

        t_166[k] = f_24 * lsk_139[k]
                   + f_3 * pc_x[k] * msk_139[k];

        t_167[k] = f_24 * lsk_140[k]
                   + f_3 * pc_x[k] * msk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, lsk_64, lsk_141, lsk_142, \
                         lsk_143, msi0_105, msi1_105, msk_136, msk_141, msk_142, \
                         msk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_24 * lsk_141[k]
                   + f_3 * pc_x[k] * msk_141[k];

        t_169[k] = f_24 * lsk_142[k]
                   + f_3 * pc_x[k] * msk_142[k];

        t_170[k] = f_24 * lsk_143[k]
                   + f_3 * pc_x[k] * msk_143[k];

        t_171[k] = f_16 * lsk_64[k]
                   + f_1 * msi0_105[k]
                   - f_2 * msi1_105[k]
                   + f_3 * pc_y[k] * msk_136[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_z, msi0_105, msi0_106, msi0_107, \
                         msi1_105, msi1_106, msi1_107, msk_136, msk_137, msk_138, \
                         msk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * msk_136[k];

        t_173[k] = f_4 * msi0_105[k]
                   - f_5 * msi1_105[k]
                   + f_3 * pc_z[k] * msk_137[k];

        t_174[k] = f_6 * msi0_106[k]
                   - f_7 * msi1_106[k]
                   + f_3 * pc_z[k] * msk_138[k];

        t_175[k] = f_8 * msi0_107[k]
                   - f_9 * msi1_107[k]
                   + f_3 * pc_z[k] * msk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, lsk_71, msi0_108, msi0_109, \
                         msi0_111, msi1_108, msi1_109, msi1_111, msk_140, msk_141, \
                         msk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * msi0_108[k]
                   - f_11 * msi1_108[k]
                   + f_3 * pc_z[k] * msk_140[k];

        t_177[k] = f_12 * msi0_109[k]
                   - f_13 * msi1_109[k]
                   + f_3 * pc_z[k] * msk_141[k];

        t_178[k] = f_16 * lsk_71[k]
                   + f_3 * pc_y[k] * msk_143[k];

        t_179[k] = f_1 * msi0_111[k]
                   - f_2 * msi1_111[k]
                   + f_3 * pc_z[k] * msk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, lsl0_48, lsl0_90, \
                         lsk_36, lsk_72, lsl1_48, lsl1_90, msk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * lsl0_90[k]
                   - f_14 * pc_y[k] * lsl1_90[k];

        t_181[k] = f_15 * lsk_72[k]
                   + f_3 * pc_y[k] * msk_144[k];

        t_182[k] = f_15 * lsk_36[k]
                   + f_3 * pc_z[k] * msk_144[k];

        t_183[k] = pa_z[k] * lsl0_48[k]
                   - f_14 * pc_z[k] * lsl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, lsl0_51, lsl0_95, \
                         lsk_39, lsk_74, lsl1_51, lsl1_95, msk_146, \
                         msk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * lsk_74[k]
                   + f_3 * pc_y[k] * msk_146[k];

        t_185[k] = pa_y[k] * lsl0_95[k]
                   - f_14 * pc_y[k] * lsl1_95[k];

        t_186[k] = pa_z[k] * lsl0_51[k]
                   - f_14 * pc_z[k] * lsl1_51[k];

        t_187[k] = f_15 * lsk_39[k]
                   + f_3 * pc_z[k] * msk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, lsl0_55, lsl0_99, \
                         lsk_42, lsk_77, lsl1_55, lsl1_99, msk_149, \
                         msk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * lsk_77[k]
                   + f_3 * pc_y[k] * msk_149[k];

        t_189[k] = pa_y[k] * lsl0_99[k]
                   - f_14 * pc_y[k] * lsl1_99[k];

        t_190[k] = pa_z[k] * lsl0_55[k]
                   - f_14 * pc_z[k] * lsl1_55[k];

        t_191[k] = f_15 * lsk_42[k]
                   + f_3 * pc_z[k] * msk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, lsl0_102, lsl0_104, lsk_80, lsk_81, \
                         lsl1_102, lsl1_104, msk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * lsl0_102[k]
                   + f_16 * lsk_80[k]
                   - f_14 * pc_y[k] * lsl1_102[k];

        t_193[k] = f_15 * lsk_81[k]
                   + f_3 * pc_y[k] * msk_153[k];

        t_194[k] = pa_y[k] * lsl0_104[k]
                   - f_14 * pc_y[k] * lsl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pa_z, pc_y, pc_z, lsl0_60, lsl0_107, \
                         lsk_46, lsk_84, lsl1_60, lsl1_107, msk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * lsl0_60[k]
                   - f_14 * pc_z[k] * lsl1_60[k];

        t_196[k] = f_15 * lsk_46[k]
                   + f_3 * pc_z[k] * msk_154[k];

        t_197[k] = pa_y[k] * lsl0_107[k]
                   + f_17 * lsk_84[k]
                   - f_14 * pc_y[k] * lsl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pc_y, lsl0_108, lsl0_110, lsk_85, lsk_86, \
                         lsl1_108, lsl1_110, msk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_y[k] * lsl0_108[k]
                   + f_16 * lsk_85[k]
                   - f_14 * pc_y[k] * lsl1_108[k];

        t_199[k] = f_15 * lsk_86[k]
                   + f_3 * pc_y[k] * msk_158[k];

        t_200[k] = pa_y[k] * lsl0_110[k]
                   - f_14 * pc_y[k] * lsl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pc_y, pc_z, lsl0_66, lsl0_113, \
                         lsk_51, lsk_89, lsl1_66, lsl1_113, msk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * lsl0_66[k]
                   - f_14 * pc_z[k] * lsl1_66[k];

        t_202[k] = f_15 * lsk_51[k]
                   + f_3 * pc_z[k] * msk_159[k];

        t_203[k] = pa_y[k] * lsl0_113[k]
                   + f_18 * lsk_89[k]
                   - f_14 * pc_y[k] * lsl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_y, lsl0_114, lsl0_115, lsl0_117, \
                         lsk_90, lsk_91, lsk_92, lsl1_114, lsl1_115, lsl1_117, \
                         msk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * lsl0_114[k]
                   + f_17 * lsk_90[k]
                   - f_14 * pc_y[k] * lsl1_114[k];

        t_205[k] = pa_y[k] * lsl0_115[k]
                   + f_16 * lsk_91[k]
                   - f_14 * pc_y[k] * lsl1_115[k];

        t_206[k] = f_15 * lsk_92[k]
                   + f_3 * pc_y[k] * msk_164[k];

        t_207[k] = pa_y[k] * lsl0_117[k]
                   - f_14 * pc_y[k] * lsl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, lsk_172, lsk_173, lsk_174, \
                         lsk_175, lsk_176, msk_172, msk_173, msk_174, msk_175, \
                         msk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_24 * lsk_172[k]
                   + f_3 * pc_x[k] * msk_172[k];

        t_209[k] = f_24 * lsk_173[k]
                   + f_3 * pc_x[k] * msk_173[k];

        t_210[k] = f_24 * lsk_174[k]
                   + f_3 * pc_x[k] * msk_174[k];

        t_211[k] = f_24 * lsk_175[k]
                   + f_3 * pc_x[k] * msk_175[k];

        t_212[k] = f_24 * lsk_176[k]
                   + f_3 * pc_x[k] * msk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, lsl0_81, lsk_177, \
                         lsk_178, lsk_179, lsl1_81, msk_177, msk_178, \
                         msk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_24 * lsk_177[k]
                   + f_3 * pc_x[k] * msk_177[k];

        t_214[k] = f_24 * lsk_178[k]
                   + f_3 * pc_x[k] * msk_178[k];

        t_215[k] = f_24 * lsk_179[k]
                   + f_3 * pc_x[k] * msk_179[k];

        t_216[k] = pa_z[k] * lsl0_81[k]
                   - f_14 * pc_z[k] * lsl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, lsk_64, lsk_102, lsk_103, msi0_135, \
                         msi0_136, msi1_135, msi1_136, msk_172, msk_174, \
                         msk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * lsk_64[k]
                   + f_3 * pc_z[k] * msk_172[k];

        t_218[k] = f_15 * lsk_102[k]
                   + f_12 * msi0_135[k]
                   - f_13 * msi1_135[k]
                   + f_3 * pc_y[k] * msk_174[k];

        t_219[k] = f_15 * lsk_103[k]
                   + f_10 * msi0_136[k]
                   - f_11 * msi1_136[k]
                   + f_3 * pc_y[k] * msk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, lsk_104, lsk_105, lsk_106, msi0_137, \
                         msi0_138, msi0_139, msi1_137, msi1_138, msi1_139, msk_176, msk_177, \
                         msk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * lsk_104[k]
                   + f_8 * msi0_137[k]
                   - f_9 * msi1_137[k]
                   + f_3 * pc_y[k] * msk_176[k];

        t_221[k] = f_15 * lsk_105[k]
                   + f_6 * msi0_138[k]
                   - f_7 * msi1_138[k]
                   + f_3 * pc_y[k] * msk_177[k];

        t_222[k] = f_15 * lsk_106[k]
                   + f_4 * msi0_139[k]
                   - f_5 * msi1_139[k]
                   + f_3 * pc_y[k] * msk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, lsl0_134, lsk_107, \
                         lsk_180, lsl1_134, msi0_140, msi1_140, msk_179, \
                         msk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * lsk_107[k]
                   + f_3 * pc_y[k] * msk_179[k];

        t_224[k] = pa_y[k] * lsl0_134[k]
                   - f_14 * pc_y[k] * lsl1_134[k];

        t_225[k] = f_24 * lsk_180[k]
                   + f_1 * msi0_140[k]
                   - f_2 * msi1_140[k]
                   + f_3 * pc_x[k] * msk_180[k];

        t_226[k] = f_3 * pc_y[k] * msk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, lsk_72, msi0_140, msi1_140, msk_180, \
                         msk_181, msk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * lsk_72[k]
                   + f_3 * pc_z[k] * msk_180[k];

        t_228[k] = f_4 * msi0_140[k]
                   - f_5 * msi1_140[k]
                   + f_3 * pc_y[k] * msk_181[k];

        t_229[k] = f_3 * pc_y[k] * msk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, lsk_185, msi0_141, msi0_142, \
                         msi0_145, msi1_141, msi1_142, msi1_145, msk_183, msk_184, \
                         msk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_24 * lsk_185[k]
                   + f_12 * msi0_145[k]
                   - f_13 * msi1_145[k]
                   + f_3 * pc_x[k] * msk_185[k];

        t_231[k] = f_6 * msi0_141[k]
                   - f_7 * msi1_141[k]
                   + f_3 * pc_y[k] * msk_183[k];

        t_232[k] = f_4 * msi0_142[k]
                   - f_5 * msi1_142[k]
                   + f_3 * pc_y[k] * msk_184[k];

        t_233[k] = f_3 * pc_y[k] * msk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, lsk_189, msi0_143, msi0_144, \
                         msi0_149, msi1_143, msi1_144, msi1_149, msk_186, msk_187, \
                         msk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_24 * lsk_189[k]
                   + f_10 * msi0_149[k]
                   - f_11 * msi1_149[k]
                   + f_3 * pc_x[k] * msk_189[k];

        t_235[k] = f_8 * msi0_143[k]
                   - f_9 * msi1_143[k]
                   + f_3 * pc_y[k] * msk_186[k];

        t_236[k] = f_6 * msi0_144[k]
                   - f_7 * msi1_144[k]
                   + f_3 * pc_y[k] * msk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_y, lsk_194, msi0_145, msi0_154, \
                         msi1_145, msi1_154, msk_188, msk_189, \
                         msk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * msi0_145[k]
                   - f_5 * msi1_145[k]
                   + f_3 * pc_y[k] * msk_188[k];

        t_238[k] = f_3 * pc_y[k] * msk_189[k];

        t_239[k] = f_24 * lsk_194[k]
                   + f_8 * msi0_154[k]
                   - f_9 * msi1_154[k]
                   + f_3 * pc_x[k] * msk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, msi0_146, msi0_147, msi0_148, msi1_146, \
                         msi1_147, msi1_148, msk_190, msk_191, \
                         msk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * msi0_146[k]
                   - f_11 * msi1_146[k]
                   + f_3 * pc_y[k] * msk_190[k];

        t_241[k] = f_8 * msi0_147[k]
                   - f_9 * msi1_147[k]
                   + f_3 * pc_y[k] * msk_191[k];

        t_242[k] = f_6 * msi0_148[k]
                   - f_7 * msi1_148[k]
                   + f_3 * pc_y[k] * msk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_y, lsk_200, msi0_149, msi0_160, \
                         msi1_149, msi1_160, msk_193, msk_194, \
                         msk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * msi0_149[k]
                   - f_5 * msi1_149[k]
                   + f_3 * pc_y[k] * msk_193[k];

        t_244[k] = f_3 * pc_y[k] * msk_194[k];

        t_245[k] = f_24 * lsk_200[k]
                   + f_6 * msi0_160[k]
                   - f_7 * msi1_160[k]
                   + f_3 * pc_x[k] * msk_200[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);
    const auto f_24 = 3.5 / q;

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

    const auto *lsl0_135 = buffer.data(lsl0 + 135);
    const auto *lsl0_138 = buffer.data(lsl0 + 138);
    const auto *lsl0_141 = buffer.data(lsl0 + 141);
    const auto *lsl0_145 = buffer.data(lsl0 + 145);
    const auto *lsl0_147 = buffer.data(lsl0 + 147);
    const auto *lsl0_150 = buffer.data(lsl0 + 150);
    const auto *lsl0_152 = buffer.data(lsl0 + 152);
    const auto *lsl0_153 = buffer.data(lsl0 + 153);
    const auto *lsl0_156 = buffer.data(lsl0 + 156);
    const auto *lsl0_158 = buffer.data(lsl0 + 158);
    const auto *lsl0_159 = buffer.data(lsl0 + 159);
    const auto *lsl0_160 = buffer.data(lsl0 + 160);
    const auto *lsl0_171 = buffer.data(lsl0 + 171);
    const auto *lsl0_225 = buffer.data(lsl0 + 225);

    const auto *lsk_107 = buffer.data(lsk + 107);
    const auto *lsk_108 = buffer.data(lsk + 108);
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
    const auto *lsk_136 = buffer.data(lsk + 136);
    const auto *lsk_143 = buffer.data(lsk + 143);
    const auto *lsk_144 = buffer.data(lsk + 144);
    const auto *lsk_146 = buffer.data(lsk + 146);
    const auto *lsk_149 = buffer.data(lsk + 149);
    const auto *lsk_153 = buffer.data(lsk + 153);
    const auto *lsk_158 = buffer.data(lsk + 158);
    const auto *lsk_164 = buffer.data(lsk + 164);
    const auto *lsk_174 = buffer.data(lsk + 174);
    const auto *lsk_175 = buffer.data(lsk + 175);
    const auto *lsk_176 = buffer.data(lsk + 176);
    const auto *lsk_177 = buffer.data(lsk + 177);
    const auto *lsk_178 = buffer.data(lsk + 178);
    const auto *lsk_179 = buffer.data(lsk + 179);
    const auto *lsk_180 = buffer.data(lsk + 180);
    const auto *lsk_207 = buffer.data(lsk + 207);
    const auto *lsk_208 = buffer.data(lsk + 208);
    const auto *lsk_209 = buffer.data(lsk + 209);
    const auto *lsk_210 = buffer.data(lsk + 210);
    const auto *lsk_211 = buffer.data(lsk + 211);
    const auto *lsk_212 = buffer.data(lsk + 212);
    const auto *lsk_213 = buffer.data(lsk + 213);
    const auto *lsk_215 = buffer.data(lsk + 215);
    const auto *lsk_216 = buffer.data(lsk + 216);
    const auto *lsk_219 = buffer.data(lsk + 219);
    const auto *lsk_222 = buffer.data(lsk + 222);
    const auto *lsk_226 = buffer.data(lsk + 226);
    const auto *lsk_231 = buffer.data(lsk + 231);
    const auto *lsk_237 = buffer.data(lsk + 237);
    const auto *lsk_244 = buffer.data(lsk + 244);
    const auto *lsk_246 = buffer.data(lsk + 246);
    const auto *lsk_247 = buffer.data(lsk + 247);
    const auto *lsk_248 = buffer.data(lsk + 248);
    const auto *lsk_249 = buffer.data(lsk + 249);
    const auto *lsk_250 = buffer.data(lsk + 250);
    const auto *lsk_251 = buffer.data(lsk + 251);
    const auto *lsk_257 = buffer.data(lsk + 257);
    const auto *lsk_261 = buffer.data(lsk + 261);
    const auto *lsk_266 = buffer.data(lsk + 266);
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

    const auto *lsl1_135 = buffer.data(lsl1 + 135);
    const auto *lsl1_138 = buffer.data(lsl1 + 138);
    const auto *lsl1_141 = buffer.data(lsl1 + 141);
    const auto *lsl1_145 = buffer.data(lsl1 + 145);
    const auto *lsl1_147 = buffer.data(lsl1 + 147);
    const auto *lsl1_150 = buffer.data(lsl1 + 150);
    const auto *lsl1_152 = buffer.data(lsl1 + 152);
    const auto *lsl1_153 = buffer.data(lsl1 + 153);
    const auto *lsl1_156 = buffer.data(lsl1 + 156);
    const auto *lsl1_158 = buffer.data(lsl1 + 158);
    const auto *lsl1_159 = buffer.data(lsl1 + 159);
    const auto *lsl1_160 = buffer.data(lsl1 + 160);
    const auto *lsl1_171 = buffer.data(lsl1 + 171);
    const auto *lsl1_225 = buffer.data(lsl1 + 225);

    const auto *msi0_150 = buffer.data(msi0 + 150);
    const auto *msi0_151 = buffer.data(msi0 + 151);
    const auto *msi0_152 = buffer.data(msi0 + 152);
    const auto *msi0_153 = buffer.data(msi0 + 153);
    const auto *msi0_154 = buffer.data(msi0 + 154);
    const auto *msi0_161 = buffer.data(msi0 + 161);
    const auto *msi0_162 = buffer.data(msi0 + 162);
    const auto *msi0_163 = buffer.data(msi0 + 163);
    const auto *msi0_164 = buffer.data(msi0 + 164);
    const auto *msi0_165 = buffer.data(msi0 + 165);
    const auto *msi0_166 = buffer.data(msi0 + 166);
    const auto *msi0_167 = buffer.data(msi0 + 167);
    const auto *msi0_168 = buffer.data(msi0 + 168);
    const auto *msi0_170 = buffer.data(msi0 + 170);
    const auto *msi0_171 = buffer.data(msi0 + 171);
    const auto *msi0_173 = buffer.data(msi0 + 173);
    const auto *msi0_174 = buffer.data(msi0 + 174);
    const auto *msi0_175 = buffer.data(msi0 + 175);
    const auto *msi0_177 = buffer.data(msi0 + 177);
    const auto *msi0_178 = buffer.data(msi0 + 178);
    const auto *msi0_179 = buffer.data(msi0 + 179);
    const auto *msi0_180 = buffer.data(msi0 + 180);
    const auto *msi0_182 = buffer.data(msi0 + 182);
    const auto *msi0_183 = buffer.data(msi0 + 183);
    const auto *msi0_189 = buffer.data(msi0 + 189);
    const auto *msi0_190 = buffer.data(msi0 + 190);
    const auto *msi0_191 = buffer.data(msi0 + 191);
    const auto *msi0_192 = buffer.data(msi0 + 192);
    const auto *msi0_193 = buffer.data(msi0 + 193);
    const auto *msi0_195 = buffer.data(msi0 + 195);
    const auto *msi0_201 = buffer.data(msi0 + 201);
    const auto *msi0_205 = buffer.data(msi0 + 205);
    const auto *msi0_210 = buffer.data(msi0 + 210);
    const auto *msi0_216 = buffer.data(msi0 + 216);
    const auto *msi0_219 = buffer.data(msi0 + 219);
    const auto *msi0_220 = buffer.data(msi0 + 220);
    const auto *msi0_221 = buffer.data(msi0 + 221);
    const auto *msi0_222 = buffer.data(msi0 + 222);
    const auto *msi0_223 = buffer.data(msi0 + 223);

    const auto *msi1_150 = buffer.data(msi1 + 150);
    const auto *msi1_151 = buffer.data(msi1 + 151);
    const auto *msi1_152 = buffer.data(msi1 + 152);
    const auto *msi1_153 = buffer.data(msi1 + 153);
    const auto *msi1_154 = buffer.data(msi1 + 154);
    const auto *msi1_161 = buffer.data(msi1 + 161);
    const auto *msi1_162 = buffer.data(msi1 + 162);
    const auto *msi1_163 = buffer.data(msi1 + 163);
    const auto *msi1_164 = buffer.data(msi1 + 164);
    const auto *msi1_165 = buffer.data(msi1 + 165);
    const auto *msi1_166 = buffer.data(msi1 + 166);
    const auto *msi1_167 = buffer.data(msi1 + 167);
    const auto *msi1_168 = buffer.data(msi1 + 168);
    const auto *msi1_170 = buffer.data(msi1 + 170);
    const auto *msi1_171 = buffer.data(msi1 + 171);
    const auto *msi1_173 = buffer.data(msi1 + 173);
    const auto *msi1_174 = buffer.data(msi1 + 174);
    const auto *msi1_175 = buffer.data(msi1 + 175);
    const auto *msi1_177 = buffer.data(msi1 + 177);
    const auto *msi1_178 = buffer.data(msi1 + 178);
    const auto *msi1_179 = buffer.data(msi1 + 179);
    const auto *msi1_180 = buffer.data(msi1 + 180);
    const auto *msi1_182 = buffer.data(msi1 + 182);
    const auto *msi1_183 = buffer.data(msi1 + 183);
    const auto *msi1_189 = buffer.data(msi1 + 189);
    const auto *msi1_190 = buffer.data(msi1 + 190);
    const auto *msi1_191 = buffer.data(msi1 + 191);
    const auto *msi1_192 = buffer.data(msi1 + 192);
    const auto *msi1_193 = buffer.data(msi1 + 193);
    const auto *msi1_195 = buffer.data(msi1 + 195);
    const auto *msi1_201 = buffer.data(msi1 + 201);
    const auto *msi1_205 = buffer.data(msi1 + 205);
    const auto *msi1_210 = buffer.data(msi1 + 210);
    const auto *msi1_216 = buffer.data(msi1 + 216);
    const auto *msi1_219 = buffer.data(msi1 + 219);
    const auto *msi1_220 = buffer.data(msi1 + 220);
    const auto *msi1_221 = buffer.data(msi1 + 221);
    const auto *msi1_222 = buffer.data(msi1 + 222);
    const auto *msi1_223 = buffer.data(msi1 + 223);

    const auto *msk_195 = buffer.data(msk + 195);
    const auto *msk_196 = buffer.data(msk + 196);
    const auto *msk_197 = buffer.data(msk + 197);
    const auto *msk_198 = buffer.data(msk + 198);
    const auto *msk_199 = buffer.data(msk + 199);
    const auto *msk_200 = buffer.data(msk + 200);
    const auto *msk_207 = buffer.data(msk + 207);
    const auto *msk_208 = buffer.data(msk + 208);
    const auto *msk_209 = buffer.data(msk + 209);
    const auto *msk_210 = buffer.data(msk + 210);
    const auto *msk_211 = buffer.data(msk + 211);
    const auto *msk_212 = buffer.data(msk + 212);
    const auto *msk_213 = buffer.data(msk + 213);
    const auto *msk_214 = buffer.data(msk + 214);
    const auto *msk_215 = buffer.data(msk + 215);
    const auto *msk_216 = buffer.data(msk + 216);
    const auto *msk_217 = buffer.data(msk + 217);
    const auto *msk_218 = buffer.data(msk + 218);
    const auto *msk_219 = buffer.data(msk + 219);
    const auto *msk_221 = buffer.data(msk + 221);
    const auto *msk_222 = buffer.data(msk + 222);
    const auto *msk_223 = buffer.data(msk + 223);
    const auto *msk_225 = buffer.data(msk + 225);
    const auto *msk_226 = buffer.data(msk + 226);
    const auto *msk_227 = buffer.data(msk + 227);
    const auto *msk_228 = buffer.data(msk + 228);
    const auto *msk_230 = buffer.data(msk + 230);
    const auto *msk_231 = buffer.data(msk + 231);
    const auto *msk_232 = buffer.data(msk + 232);
    const auto *msk_233 = buffer.data(msk + 233);
    const auto *msk_234 = buffer.data(msk + 234);
    const auto *msk_236 = buffer.data(msk + 236);
    const auto *msk_237 = buffer.data(msk + 237);
    const auto *msk_244 = buffer.data(msk + 244);
    const auto *msk_245 = buffer.data(msk + 245);
    const auto *msk_246 = buffer.data(msk + 246);
    const auto *msk_247 = buffer.data(msk + 247);
    const auto *msk_248 = buffer.data(msk + 248);
    const auto *msk_249 = buffer.data(msk + 249);
    const auto *msk_250 = buffer.data(msk + 250);
    const auto *msk_251 = buffer.data(msk + 251);
    const auto *msk_252 = buffer.data(msk + 252);
    const auto *msk_254 = buffer.data(msk + 254);
    const auto *msk_255 = buffer.data(msk + 255);
    const auto *msk_257 = buffer.data(msk + 257);
    const auto *msk_258 = buffer.data(msk + 258);
    const auto *msk_261 = buffer.data(msk + 261);
    const auto *msk_262 = buffer.data(msk + 262);
    const auto *msk_266 = buffer.data(msk + 266);
    const auto *msk_267 = buffer.data(msk + 267);
    const auto *msk_272 = buffer.data(msk + 272);
    const auto *msk_279 = buffer.data(msk + 279);
    const auto *msk_280 = buffer.data(msk + 280);
    const auto *msk_281 = buffer.data(msk + 281);
    const auto *msk_282 = buffer.data(msk + 282);
    const auto *msk_283 = buffer.data(msk + 283);
    const auto *msk_284 = buffer.data(msk + 284);
    const auto *msk_285 = buffer.data(msk + 285);
    const auto *msk_286 = buffer.data(msk + 286);
    const auto *msk_287 = buffer.data(msk + 287);
    const auto *msk_288 = buffer.data(msk + 288);

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, msi0_150, msi0_151, msi0_152, msi1_150, \
                         msi1_151, msi1_152, msk_195, msk_196, \
                         msk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * msi0_150[k]
                   - f_13 * msi1_150[k]
                   + f_3 * pc_y[k] * msk_195[k];

        t_247[k] = f_10 * msi0_151[k]
                   - f_11 * msi1_151[k]
                   + f_3 * pc_y[k] * msk_196[k];

        t_248[k] = f_8 * msi0_152[k]
                   - f_9 * msi1_152[k]
                   + f_3 * pc_y[k] * msk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, msi0_153, msi0_154, msi1_153, msi1_154, \
                         msk_198, msk_199, msk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * msi0_153[k]
                   - f_7 * msi1_153[k]
                   + f_3 * pc_y[k] * msk_198[k];

        t_250[k] = f_4 * msi0_154[k]
                   - f_5 * msi1_154[k]
                   + f_3 * pc_y[k] * msk_199[k];

        t_251[k] = f_3 * pc_y[k] * msk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, lsk_207, lsk_208, lsk_209, lsk_210, \
                         msi0_167, msi1_167, msk_207, msk_208, msk_209, \
                         msk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_24 * lsk_207[k]
                   + f_4 * msi0_167[k]
                   - f_5 * msi1_167[k]
                   + f_3 * pc_x[k] * msk_207[k];

        t_253[k] = f_24 * lsk_208[k]
                   + f_3 * pc_x[k] * msk_208[k];

        t_254[k] = f_24 * lsk_209[k]
                   + f_3 * pc_x[k] * msk_209[k];

        t_255[k] = f_24 * lsk_210[k]
                   + f_3 * pc_x[k] * msk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, lsk_211, lsk_212, \
                         lsk_213, lsk_215, msk_207, msk_211, msk_212, msk_213, \
                         msk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_24 * lsk_211[k]
                   + f_3 * pc_x[k] * msk_211[k];

        t_257[k] = f_24 * lsk_212[k]
                   + f_3 * pc_x[k] * msk_212[k];

        t_258[k] = f_24 * lsk_213[k]
                   + f_3 * pc_x[k] * msk_213[k];

        t_259[k] = f_3 * pc_y[k] * msk_207[k];

        t_260[k] = f_24 * lsk_215[k]
                   + f_3 * pc_x[k] * msk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_y, msi0_161, msi0_162, msi0_163, msi1_161, \
                         msi1_162, msi1_163, msk_208, msk_209, \
                         msk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * msi0_161[k]
                   - f_2 * msi1_161[k]
                   + f_3 * pc_y[k] * msk_208[k];

        t_262[k] = f_22 * msi0_162[k]
                   - f_23 * msi1_162[k]
                   + f_3 * pc_y[k] * msk_209[k];

        t_263[k] = f_12 * msi0_163[k]
                   - f_13 * msi1_163[k]
                   + f_3 * pc_y[k] * msk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, msi0_164, msi0_165, msi0_166, msi1_164, \
                         msi1_165, msi1_166, msk_211, msk_212, \
                         msk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * msi0_164[k]
                   - f_11 * msi1_164[k]
                   + f_3 * pc_y[k] * msk_211[k];

        t_265[k] = f_8 * msi0_165[k]
                   - f_9 * msi1_165[k]
                   + f_3 * pc_y[k] * msk_212[k];

        t_266[k] = f_6 * msi0_166[k]
                   - f_7 * msi1_166[k]
                   + f_3 * pc_y[k] * msk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, lsk_107, lsk_216, \
                         msi0_167, msi0_168, msi1_167, msi1_168, msk_214, msk_215, \
                         msk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * msi0_167[k]
                   - f_5 * msi1_167[k]
                   + f_3 * pc_y[k] * msk_214[k];

        t_268[k] = f_3 * pc_y[k] * msk_215[k];

        t_269[k] = f_16 * lsk_107[k]
                   + f_1 * msi0_167[k]
                   - f_2 * msi1_167[k]
                   + f_3 * pc_z[k] * msk_215[k];

        t_270[k] = f_20 * lsk_216[k]
                   + f_1 * msi0_168[k]
                   - f_2 * msi1_168[k]
                   + f_3 * pc_x[k] * msk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, lsk_108, lsk_219, \
                         msi0_171, msi1_171, msk_216, msk_217, \
                         msk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * lsk_108[k]
                   + f_3 * pc_y[k] * msk_216[k];

        t_272[k] = f_3 * pc_z[k] * msk_216[k];

        t_273[k] = f_20 * lsk_219[k]
                   + f_12 * msi0_171[k]
                   - f_13 * msi1_171[k]
                   + f_3 * pc_x[k] * msk_219[k];

        t_274[k] = f_3 * pc_z[k] * msk_217[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, lsk_222, msi0_168, msi0_174, \
                         msi1_168, msi1_174, msk_218, msk_219, \
                         msk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_4 * msi0_168[k]
                   - f_5 * msi1_168[k]
                   + f_3 * pc_z[k] * msk_218[k];

        t_276[k] = f_20 * lsk_222[k]
                   + f_10 * msi0_174[k]
                   - f_11 * msi1_174[k]
                   + f_3 * pc_x[k] * msk_222[k];

        t_277[k] = f_3 * pc_z[k] * msk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pc_x, pc_y, pc_z, lsk_113, lsk_226, \
                         msi0_170, msi0_178, msi1_170, msi1_178, msk_221, msk_222, \
                         msk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * lsk_113[k]
                   + f_3 * pc_y[k] * msk_221[k];

        t_279[k] = f_6 * msi0_170[k]
                   - f_7 * msi1_170[k]
                   + f_3 * pc_z[k] * msk_221[k];

        t_280[k] = f_20 * lsk_226[k]
                   + f_8 * msi0_178[k]
                   - f_9 * msi1_178[k]
                   + f_3 * pc_x[k] * msk_226[k];

        t_281[k] = f_3 * pc_z[k] * msk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pc_y, pc_z, lsk_117, msi0_171, msi0_173, \
                         msi1_171, msi1_173, msk_223, msk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * msi0_171[k]
                   - f_5 * msi1_171[k]
                   + f_3 * pc_z[k] * msk_223[k];

        t_283[k] = f_17 * lsk_117[k]
                   + f_3 * pc_y[k] * msk_225[k];

        t_284[k] = f_8 * msi0_173[k]
                   - f_9 * msi1_173[k]
                   + f_3 * pc_z[k] * msk_225[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_x, pc_z, lsk_231, msi0_174, msi0_183, \
                         msi1_174, msi1_183, msk_226, msk_227, \
                         msk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_20 * lsk_231[k]
                   + f_6 * msi0_183[k]
                   - f_7 * msi1_183[k]
                   + f_3 * pc_x[k] * msk_231[k];

        t_286[k] = f_3 * pc_z[k] * msk_226[k];

        t_287[k] = f_4 * msi0_174[k]
                   - f_5 * msi1_174[k]
                   + f_3 * pc_z[k] * msk_227[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, lsk_122, msi0_175, msi0_177, \
                         msi1_175, msi1_177, msk_228, msk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * msi0_175[k]
                   - f_7 * msi1_175[k]
                   + f_3 * pc_z[k] * msk_228[k];

        t_289[k] = f_17 * lsk_122[k]
                   + f_3 * pc_y[k] * msk_230[k];

        t_290[k] = f_10 * msi0_177[k]
                   - f_11 * msi1_177[k]
                   + f_3 * pc_z[k] * msk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_x, pc_z, lsk_237, msi0_178, msi0_189, \
                         msi1_178, msi1_189, msk_231, msk_232, \
                         msk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_20 * lsk_237[k]
                   + f_4 * msi0_189[k]
                   - f_5 * msi1_189[k]
                   + f_3 * pc_x[k] * msk_237[k];

        t_292[k] = f_3 * pc_z[k] * msk_231[k];

        t_293[k] = f_4 * msi0_178[k]
                   - f_5 * msi1_178[k]
                   + f_3 * pc_z[k] * msk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_y, pc_z, lsk_128, msi0_179, msi0_180, \
                         msi0_182, msi1_179, msi1_180, msi1_182, msk_233, msk_234, \
                         msk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * msi0_179[k]
                   - f_7 * msi1_179[k]
                   + f_3 * pc_z[k] * msk_233[k];

        t_295[k] = f_8 * msi0_180[k]
                   - f_9 * msi1_180[k]
                   + f_3 * pc_z[k] * msk_234[k];

        t_296[k] = f_17 * lsk_128[k]
                   + f_3 * pc_y[k] * msk_236[k];

        t_297[k] = f_12 * msi0_182[k]
                   - f_13 * msi1_182[k]
                   + f_3 * pc_z[k] * msk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_z, lsk_244, lsk_246, \
                         lsk_247, lsk_248, msk_237, msk_244, msk_246, msk_247, \
                         msk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_20 * lsk_244[k]
                   + f_3 * pc_x[k] * msk_244[k];

        t_299[k] = f_3 * pc_z[k] * msk_237[k];

        t_300[k] = f_20 * lsk_246[k]
                   + f_3 * pc_x[k] * msk_246[k];

        t_301[k] = f_20 * lsk_247[k]
                   + f_3 * pc_x[k] * msk_247[k];

        t_302[k] = f_20 * lsk_248[k]
                   + f_3 * pc_x[k] * msk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, pc_y, lsk_136, lsk_249, lsk_250, \
                         lsk_251, msi0_189, msi1_189, msk_244, msk_249, msk_250, \
                         msk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_20 * lsk_249[k]
                   + f_3 * pc_x[k] * msk_249[k];

        t_304[k] = f_20 * lsk_250[k]
                   + f_3 * pc_x[k] * msk_250[k];

        t_305[k] = f_20 * lsk_251[k]
                   + f_3 * pc_x[k] * msk_251[k];

        t_306[k] = f_17 * lsk_136[k]
                   + f_1 * msi0_189[k]
                   - f_2 * msi1_189[k]
                   + f_3 * pc_y[k] * msk_244[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_z, msi0_189, msi0_190, msi0_191, \
                         msi1_189, msi1_190, msi1_191, msk_244, msk_245, msk_246, \
                         msk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_z[k] * msk_244[k];

        t_308[k] = f_4 * msi0_189[k]
                   - f_5 * msi1_189[k]
                   + f_3 * pc_z[k] * msk_245[k];

        t_309[k] = f_6 * msi0_190[k]
                   - f_7 * msi1_190[k]
                   + f_3 * pc_z[k] * msk_246[k];

        t_310[k] = f_8 * msi0_191[k]
                   - f_9 * msi1_191[k]
                   + f_3 * pc_z[k] * msk_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, lsk_143, msi0_192, msi0_193, \
                         msi0_195, msi1_192, msi1_193, msi1_195, msk_248, msk_249, \
                         msk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_10 * msi0_192[k]
                   - f_11 * msi1_192[k]
                   + f_3 * pc_z[k] * msk_248[k];

        t_312[k] = f_12 * msi0_193[k]
                   - f_13 * msi1_193[k]
                   + f_3 * pc_z[k] * msk_249[k];

        t_313[k] = f_17 * lsk_143[k]
                   + f_3 * pc_y[k] * msk_251[k];

        t_314[k] = f_1 * msi0_195[k]
                   - f_2 * msi1_195[k]
                   + f_3 * pc_z[k] * msk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_z, pc_y, pc_z, lsl0_135, lsl0_138, \
                         lsk_108, lsk_144, lsl1_135, lsl1_138, \
                         msk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * lsl0_135[k]
                   - f_14 * pc_z[k] * lsl1_135[k];

        t_316[k] = f_16 * lsk_144[k]
                   + f_3 * pc_y[k] * msk_252[k];

        t_317[k] = f_15 * lsk_108[k]
                   + f_3 * pc_z[k] * msk_252[k];

        t_318[k] = pa_z[k] * lsl0_138[k]
                   - f_14 * pc_z[k] * lsl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_x, pc_y, pc_z, lsl0_141, lsk_146, \
                         lsk_257, lsl1_141, msi0_201, msi1_201, msk_254, \
                         msk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * lsk_146[k]
                   + f_3 * pc_y[k] * msk_254[k];

        t_320[k] = f_20 * lsk_257[k]
                   + f_12 * msi0_201[k]
                   - f_13 * msi1_201[k]
                   + f_3 * pc_x[k] * msk_257[k];

        t_321[k] = pa_z[k] * lsl0_141[k]
                   - f_14 * pc_z[k] * lsl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, lsk_111, lsk_149, lsk_261, \
                         msi0_205, msi1_205, msk_255, msk_257, \
                         msk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * lsk_111[k]
                   + f_3 * pc_z[k] * msk_255[k];

        t_323[k] = f_16 * lsk_149[k]
                   + f_3 * pc_y[k] * msk_257[k];

        t_324[k] = f_20 * lsk_261[k]
                   + f_10 * msi0_205[k]
                   - f_11 * msi1_205[k]
                   + f_3 * pc_x[k] * msk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_z, pc_y, pc_z, lsl0_145, lsl0_147, \
                         lsk_114, lsk_115, lsk_153, lsl1_145, lsl1_147, msk_258, \
                         msk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * lsl0_145[k]
                   - f_14 * pc_z[k] * lsl1_145[k];

        t_326[k] = f_15 * lsk_114[k]
                   + f_3 * pc_z[k] * msk_258[k];

        t_327[k] = pa_z[k] * lsl0_147[k]
                   + f_16 * lsk_115[k]
                   - f_14 * pc_z[k] * lsl1_147[k];

        t_328[k] = f_16 * lsk_153[k]
                   + f_3 * pc_y[k] * msk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pc_x, pc_z, lsl0_150, lsk_118, lsk_266, \
                         lsl1_150, msi0_210, msi1_210, msk_262, \
                         msk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_20 * lsk_266[k]
                   + f_8 * msi0_210[k]
                   - f_9 * msi1_210[k]
                   + f_3 * pc_x[k] * msk_266[k];

        t_330[k] = pa_z[k] * lsl0_150[k]
                   - f_14 * pc_z[k] * lsl1_150[k];

        t_331[k] = f_15 * lsk_118[k]
                   + f_3 * pc_z[k] * msk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pa_z, pc_y, pc_z, lsl0_152, lsl0_153, lsk_119, \
                         lsk_120, lsk_158, lsl1_152, lsl1_153, \
                         msk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_z[k] * lsl0_152[k]
                   + f_16 * lsk_119[k]
                   - f_14 * pc_z[k] * lsl1_152[k];

        t_333[k] = pa_z[k] * lsl0_153[k]
                   + f_17 * lsk_120[k]
                   - f_14 * pc_z[k] * lsl1_153[k];

        t_334[k] = f_16 * lsk_158[k]
                   + f_3 * pc_y[k] * msk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_x, pc_z, lsl0_156, lsk_123, lsk_272, \
                         lsl1_156, msi0_216, msi1_216, msk_267, \
                         msk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_20 * lsk_272[k]
                   + f_6 * msi0_216[k]
                   - f_7 * msi1_216[k]
                   + f_3 * pc_x[k] * msk_272[k];

        t_336[k] = pa_z[k] * lsl0_156[k]
                   - f_14 * pc_z[k] * lsl1_156[k];

        t_337[k] = f_15 * lsk_123[k]
                   + f_3 * pc_z[k] * msk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_z, pc_z, lsl0_158, lsl0_159, lsl0_160, \
                         lsk_124, lsk_125, lsk_126, lsl1_158, lsl1_159, \
                         lsl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_z[k] * lsl0_158[k]
                   + f_16 * lsk_124[k]
                   - f_14 * pc_z[k] * lsl1_158[k];

        t_339[k] = pa_z[k] * lsl0_159[k]
                   + f_17 * lsk_125[k]
                   - f_14 * pc_z[k] * lsl1_159[k];

        t_340[k] = pa_z[k] * lsl0_160[k]
                   + f_18 * lsk_126[k]
                   - f_14 * pc_z[k] * lsl1_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, lsk_164, lsk_279, lsk_280, \
                         lsk_281, msi0_223, msi1_223, msk_272, msk_279, msk_280, \
                         msk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * lsk_164[k]
                   + f_3 * pc_y[k] * msk_272[k];

        t_342[k] = f_20 * lsk_279[k]
                   + f_4 * msi0_223[k]
                   - f_5 * msi1_223[k]
                   + f_3 * pc_x[k] * msk_279[k];

        t_343[k] = f_20 * lsk_280[k]
                   + f_3 * pc_x[k] * msk_280[k];

        t_344[k] = f_20 * lsk_281[k]
                   + f_3 * pc_x[k] * msk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, lsk_282, lsk_283, lsk_284, \
                         lsk_285, lsk_286, msk_282, msk_283, msk_284, msk_285, \
                         msk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_20 * lsk_282[k]
                   + f_3 * pc_x[k] * msk_282[k];

        t_346[k] = f_20 * lsk_283[k]
                   + f_3 * pc_x[k] * msk_283[k];

        t_347[k] = f_20 * lsk_284[k]
                   + f_3 * pc_x[k] * msk_284[k];

        t_348[k] = f_20 * lsk_285[k]
                   + f_3 * pc_x[k] * msk_285[k];

        t_349[k] = f_20 * lsk_286[k]
                   + f_3 * pc_x[k] * msk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pa_z, pc_x, pc_z, lsl0_171, lsk_136, lsk_287, \
                         lsl1_171, msk_280, msk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_20 * lsk_287[k]
                   + f_3 * pc_x[k] * msk_287[k];

        t_351[k] = pa_z[k] * lsl0_171[k]
                   - f_14 * pc_z[k] * lsl1_171[k];

        t_352[k] = f_15 * lsk_136[k]
                   + f_3 * pc_z[k] * msk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, lsk_174, lsk_175, lsk_176, msi0_219, \
                         msi0_220, msi0_221, msi1_219, msi1_220, msi1_221, msk_282, msk_283, \
                         msk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * lsk_174[k]
                   + f_12 * msi0_219[k]
                   - f_13 * msi1_219[k]
                   + f_3 * pc_y[k] * msk_282[k];

        t_354[k] = f_16 * lsk_175[k]
                   + f_10 * msi0_220[k]
                   - f_11 * msi1_220[k]
                   + f_3 * pc_y[k] * msk_283[k];

        t_355[k] = f_16 * lsk_176[k]
                   + f_8 * msi0_221[k]
                   - f_9 * msi1_221[k]
                   + f_3 * pc_y[k] * msk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, lsk_177, lsk_178, lsk_179, msi0_222, \
                         msi0_223, msi1_222, msi1_223, msk_285, msk_286, \
                         msk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * lsk_177[k]
                   + f_6 * msi0_222[k]
                   - f_7 * msi1_222[k]
                   + f_3 * pc_y[k] * msk_285[k];

        t_357[k] = f_16 * lsk_178[k]
                   + f_4 * msi0_223[k]
                   - f_5 * msi1_223[k]
                   + f_3 * pc_y[k] * msk_286[k];

        t_358[k] = f_16 * lsk_179[k]
                   + f_3 * pc_y[k] * msk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_y, pc_y, pc_z, lsl0_225, lsk_143, \
                         lsk_144, lsk_180, lsl1_225, msi0_223, msi1_223, msk_287, \
                         msk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * lsk_143[k]
                   + f_1 * msi0_223[k]
                   - f_2 * msi1_223[k]
                   + f_3 * pc_z[k] * msk_287[k];

        t_360[k] = pa_y[k] * lsl0_225[k]
                   - f_14 * pc_y[k] * lsl1_225[k];

        t_361[k] = f_15 * lsk_180[k]
                   + f_3 * pc_y[k] * msk_288[k];

        t_362[k] = f_16 * lsk_144[k]
                   + f_3 * pc_z[k] * msk_288[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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

    const auto *lsl0_228 = buffer.data(lsl0 + 228);
    const auto *lsl0_230 = buffer.data(lsl0 + 230);
    const auto *lsl0_231 = buffer.data(lsl0 + 231);
    const auto *lsl0_234 = buffer.data(lsl0 + 234);
    const auto *lsl0_235 = buffer.data(lsl0 + 235);
    const auto *lsl0_237 = buffer.data(lsl0 + 237);
    const auto *lsl0_239 = buffer.data(lsl0 + 239);
    const auto *lsl0_240 = buffer.data(lsl0 + 240);
    const auto *lsl0_242 = buffer.data(lsl0 + 242);
    const auto *lsl0_243 = buffer.data(lsl0 + 243);
    const auto *lsl0_245 = buffer.data(lsl0 + 245);
    const auto *lsl0_246 = buffer.data(lsl0 + 246);
    const auto *lsl0_248 = buffer.data(lsl0 + 248);
    const auto *lsl0_249 = buffer.data(lsl0 + 249);
    const auto *lsl0_250 = buffer.data(lsl0 + 250);
    const auto *lsl0_252 = buffer.data(lsl0 + 252);
    const auto *lsl0_269 = buffer.data(lsl0 + 269);

    const auto *lsk_147 = buffer.data(lsk + 147);
    const auto *lsk_150 = buffer.data(lsk + 150);
    const auto *lsk_154 = buffer.data(lsk + 154);
    const auto *lsk_159 = buffer.data(lsk + 159);
    const auto *lsk_172 = buffer.data(lsk + 172);
    const auto *lsk_180 = buffer.data(lsk + 180);
    const auto *lsk_181 = buffer.data(lsk + 181);
    const auto *lsk_182 = buffer.data(lsk + 182);
    const auto *lsk_183 = buffer.data(lsk + 183);
    const auto *lsk_185 = buffer.data(lsk + 185);
    const auto *lsk_186 = buffer.data(lsk + 186);
    const auto *lsk_188 = buffer.data(lsk + 188);
    const auto *lsk_189 = buffer.data(lsk + 189);
    const auto *lsk_190 = buffer.data(lsk + 190);
    const auto *lsk_192 = buffer.data(lsk + 192);
    const auto *lsk_193 = buffer.data(lsk + 193);
    const auto *lsk_194 = buffer.data(lsk + 194);
    const auto *lsk_195 = buffer.data(lsk + 195);
    const auto *lsk_197 = buffer.data(lsk + 197);
    const auto *lsk_198 = buffer.data(lsk + 198);
    const auto *lsk_199 = buffer.data(lsk + 199);
    const auto *lsk_200 = buffer.data(lsk + 200);
    const auto *lsk_208 = buffer.data(lsk + 208);
    const auto *lsk_210 = buffer.data(lsk + 210);
    const auto *lsk_211 = buffer.data(lsk + 211);
    const auto *lsk_212 = buffer.data(lsk + 212);
    const auto *lsk_213 = buffer.data(lsk + 213);
    const auto *lsk_214 = buffer.data(lsk + 214);
    const auto *lsk_215 = buffer.data(lsk + 215);
    const auto *lsk_216 = buffer.data(lsk + 216);
    const auto *lsk_221 = buffer.data(lsk + 221);
    const auto *lsk_225 = buffer.data(lsk + 225);
    const auto *lsk_230 = buffer.data(lsk + 230);
    const auto *lsk_236 = buffer.data(lsk + 236);
    const auto *lsk_316 = buffer.data(lsk + 316);
    const auto *lsk_317 = buffer.data(lsk + 317);
    const auto *lsk_318 = buffer.data(lsk + 318);
    const auto *lsk_319 = buffer.data(lsk + 319);
    const auto *lsk_320 = buffer.data(lsk + 320);
    const auto *lsk_321 = buffer.data(lsk + 321);
    const auto *lsk_322 = buffer.data(lsk + 322);
    const auto *lsk_323 = buffer.data(lsk + 323);
    const auto *lsk_324 = buffer.data(lsk + 324);
    const auto *lsk_329 = buffer.data(lsk + 329);
    const auto *lsk_333 = buffer.data(lsk + 333);
    const auto *lsk_338 = buffer.data(lsk + 338);
    const auto *lsk_344 = buffer.data(lsk + 344);
    const auto *lsk_351 = buffer.data(lsk + 351);
    const auto *lsk_352 = buffer.data(lsk + 352);
    const auto *lsk_353 = buffer.data(lsk + 353);
    const auto *lsk_354 = buffer.data(lsk + 354);
    const auto *lsk_355 = buffer.data(lsk + 355);
    const auto *lsk_356 = buffer.data(lsk + 356);
    const auto *lsk_357 = buffer.data(lsk + 357);
    const auto *lsk_359 = buffer.data(lsk + 359);
    const auto *lsk_360 = buffer.data(lsk + 360);
    const auto *lsk_363 = buffer.data(lsk + 363);
    const auto *lsk_366 = buffer.data(lsk + 366);
    const auto *lsk_370 = buffer.data(lsk + 370);
    const auto *lsk_375 = buffer.data(lsk + 375);
    const auto *lsk_381 = buffer.data(lsk + 381);

    const auto *lsl1_228 = buffer.data(lsl1 + 228);
    const auto *lsl1_230 = buffer.data(lsl1 + 230);
    const auto *lsl1_231 = buffer.data(lsl1 + 231);
    const auto *lsl1_234 = buffer.data(lsl1 + 234);
    const auto *lsl1_235 = buffer.data(lsl1 + 235);
    const auto *lsl1_237 = buffer.data(lsl1 + 237);
    const auto *lsl1_239 = buffer.data(lsl1 + 239);
    const auto *lsl1_240 = buffer.data(lsl1 + 240);
    const auto *lsl1_242 = buffer.data(lsl1 + 242);
    const auto *lsl1_243 = buffer.data(lsl1 + 243);
    const auto *lsl1_245 = buffer.data(lsl1 + 245);
    const auto *lsl1_246 = buffer.data(lsl1 + 246);
    const auto *lsl1_248 = buffer.data(lsl1 + 248);
    const auto *lsl1_249 = buffer.data(lsl1 + 249);
    const auto *lsl1_250 = buffer.data(lsl1 + 250);
    const auto *lsl1_252 = buffer.data(lsl1 + 252);
    const auto *lsl1_269 = buffer.data(lsl1 + 269);

    const auto *msi0_245 = buffer.data(msi0 + 245);
    const auto *msi0_247 = buffer.data(msi0 + 247);
    const auto *msi0_248 = buffer.data(msi0 + 248);
    const auto *msi0_249 = buffer.data(msi0 + 249);
    const auto *msi0_250 = buffer.data(msi0 + 250);
    const auto *msi0_251 = buffer.data(msi0 + 251);
    const auto *msi0_252 = buffer.data(msi0 + 252);
    const auto *msi0_253 = buffer.data(msi0 + 253);
    const auto *msi0_254 = buffer.data(msi0 + 254);
    const auto *msi0_255 = buffer.data(msi0 + 255);
    const auto *msi0_256 = buffer.data(msi0 + 256);
    const auto *msi0_257 = buffer.data(msi0 + 257);
    const auto *msi0_258 = buffer.data(msi0 + 258);
    const auto *msi0_259 = buffer.data(msi0 + 259);
    const auto *msi0_260 = buffer.data(msi0 + 260);
    const auto *msi0_261 = buffer.data(msi0 + 261);
    const auto *msi0_262 = buffer.data(msi0 + 262);
    const auto *msi0_263 = buffer.data(msi0 + 263);
    const auto *msi0_264 = buffer.data(msi0 + 264);
    const auto *msi0_265 = buffer.data(msi0 + 265);
    const auto *msi0_266 = buffer.data(msi0 + 266);
    const auto *msi0_272 = buffer.data(msi0 + 272);
    const auto *msi0_273 = buffer.data(msi0 + 273);
    const auto *msi0_274 = buffer.data(msi0 + 274);
    const auto *msi0_275 = buffer.data(msi0 + 275);
    const auto *msi0_276 = buffer.data(msi0 + 276);
    const auto *msi0_277 = buffer.data(msi0 + 277);
    const auto *msi0_278 = buffer.data(msi0 + 278);
    const auto *msi0_279 = buffer.data(msi0 + 279);
    const auto *msi0_280 = buffer.data(msi0 + 280);
    const auto *msi0_282 = buffer.data(msi0 + 282);
    const auto *msi0_283 = buffer.data(msi0 + 283);
    const auto *msi0_285 = buffer.data(msi0 + 285);
    const auto *msi0_286 = buffer.data(msi0 + 286);
    const auto *msi0_287 = buffer.data(msi0 + 287);
    const auto *msi0_289 = buffer.data(msi0 + 289);
    const auto *msi0_290 = buffer.data(msi0 + 290);
    const auto *msi0_291 = buffer.data(msi0 + 291);
    const auto *msi0_292 = buffer.data(msi0 + 292);
    const auto *msi0_294 = buffer.data(msi0 + 294);
    const auto *msi0_295 = buffer.data(msi0 + 295);
    const auto *msi0_301 = buffer.data(msi0 + 301);

    const auto *msi1_245 = buffer.data(msi1 + 245);
    const auto *msi1_247 = buffer.data(msi1 + 247);
    const auto *msi1_248 = buffer.data(msi1 + 248);
    const auto *msi1_249 = buffer.data(msi1 + 249);
    const auto *msi1_250 = buffer.data(msi1 + 250);
    const auto *msi1_251 = buffer.data(msi1 + 251);
    const auto *msi1_252 = buffer.data(msi1 + 252);
    const auto *msi1_253 = buffer.data(msi1 + 253);
    const auto *msi1_254 = buffer.data(msi1 + 254);
    const auto *msi1_255 = buffer.data(msi1 + 255);
    const auto *msi1_256 = buffer.data(msi1 + 256);
    const auto *msi1_257 = buffer.data(msi1 + 257);
    const auto *msi1_258 = buffer.data(msi1 + 258);
    const auto *msi1_259 = buffer.data(msi1 + 259);
    const auto *msi1_260 = buffer.data(msi1 + 260);
    const auto *msi1_261 = buffer.data(msi1 + 261);
    const auto *msi1_262 = buffer.data(msi1 + 262);
    const auto *msi1_263 = buffer.data(msi1 + 263);
    const auto *msi1_264 = buffer.data(msi1 + 264);
    const auto *msi1_265 = buffer.data(msi1 + 265);
    const auto *msi1_266 = buffer.data(msi1 + 266);
    const auto *msi1_272 = buffer.data(msi1 + 272);
    const auto *msi1_273 = buffer.data(msi1 + 273);
    const auto *msi1_274 = buffer.data(msi1 + 274);
    const auto *msi1_275 = buffer.data(msi1 + 275);
    const auto *msi1_276 = buffer.data(msi1 + 276);
    const auto *msi1_277 = buffer.data(msi1 + 277);
    const auto *msi1_278 = buffer.data(msi1 + 278);
    const auto *msi1_279 = buffer.data(msi1 + 279);
    const auto *msi1_280 = buffer.data(msi1 + 280);
    const auto *msi1_282 = buffer.data(msi1 + 282);
    const auto *msi1_283 = buffer.data(msi1 + 283);
    const auto *msi1_285 = buffer.data(msi1 + 285);
    const auto *msi1_286 = buffer.data(msi1 + 286);
    const auto *msi1_287 = buffer.data(msi1 + 287);
    const auto *msi1_289 = buffer.data(msi1 + 289);
    const auto *msi1_290 = buffer.data(msi1 + 290);
    const auto *msi1_291 = buffer.data(msi1 + 291);
    const auto *msi1_292 = buffer.data(msi1 + 292);
    const auto *msi1_294 = buffer.data(msi1 + 294);
    const auto *msi1_295 = buffer.data(msi1 + 295);
    const auto *msi1_301 = buffer.data(msi1 + 301);

    const auto *msk_290 = buffer.data(msk + 290);
    const auto *msk_291 = buffer.data(msk + 291);
    const auto *msk_293 = buffer.data(msk + 293);
    const auto *msk_294 = buffer.data(msk + 294);
    const auto *msk_297 = buffer.data(msk + 297);
    const auto *msk_298 = buffer.data(msk + 298);
    const auto *msk_302 = buffer.data(msk + 302);
    const auto *msk_303 = buffer.data(msk + 303);
    const auto *msk_308 = buffer.data(msk + 308);
    const auto *msk_316 = buffer.data(msk + 316);
    const auto *msk_317 = buffer.data(msk + 317);
    const auto *msk_318 = buffer.data(msk + 318);
    const auto *msk_319 = buffer.data(msk + 319);
    const auto *msk_320 = buffer.data(msk + 320);
    const auto *msk_321 = buffer.data(msk + 321);
    const auto *msk_322 = buffer.data(msk + 322);
    const auto *msk_323 = buffer.data(msk + 323);
    const auto *msk_324 = buffer.data(msk + 324);
    const auto *msk_325 = buffer.data(msk + 325);
    const auto *msk_326 = buffer.data(msk + 326);
    const auto *msk_327 = buffer.data(msk + 327);
    const auto *msk_328 = buffer.data(msk + 328);
    const auto *msk_329 = buffer.data(msk + 329);
    const auto *msk_330 = buffer.data(msk + 330);
    const auto *msk_331 = buffer.data(msk + 331);
    const auto *msk_332 = buffer.data(msk + 332);
    const auto *msk_333 = buffer.data(msk + 333);
    const auto *msk_334 = buffer.data(msk + 334);
    const auto *msk_335 = buffer.data(msk + 335);
    const auto *msk_336 = buffer.data(msk + 336);
    const auto *msk_337 = buffer.data(msk + 337);
    const auto *msk_338 = buffer.data(msk + 338);
    const auto *msk_339 = buffer.data(msk + 339);
    const auto *msk_340 = buffer.data(msk + 340);
    const auto *msk_341 = buffer.data(msk + 341);
    const auto *msk_342 = buffer.data(msk + 342);
    const auto *msk_343 = buffer.data(msk + 343);
    const auto *msk_344 = buffer.data(msk + 344);
    const auto *msk_351 = buffer.data(msk + 351);
    const auto *msk_352 = buffer.data(msk + 352);
    const auto *msk_353 = buffer.data(msk + 353);
    const auto *msk_354 = buffer.data(msk + 354);
    const auto *msk_355 = buffer.data(msk + 355);
    const auto *msk_356 = buffer.data(msk + 356);
    const auto *msk_357 = buffer.data(msk + 357);
    const auto *msk_358 = buffer.data(msk + 358);
    const auto *msk_359 = buffer.data(msk + 359);
    const auto *msk_360 = buffer.data(msk + 360);
    const auto *msk_361 = buffer.data(msk + 361);
    const auto *msk_362 = buffer.data(msk + 362);
    const auto *msk_363 = buffer.data(msk + 363);
    const auto *msk_365 = buffer.data(msk + 365);
    const auto *msk_366 = buffer.data(msk + 366);
    const auto *msk_367 = buffer.data(msk + 367);
    const auto *msk_369 = buffer.data(msk + 369);
    const auto *msk_370 = buffer.data(msk + 370);
    const auto *msk_371 = buffer.data(msk + 371);
    const auto *msk_372 = buffer.data(msk + 372);
    const auto *msk_374 = buffer.data(msk + 374);
    const auto *msk_375 = buffer.data(msk + 375);
    const auto *msk_376 = buffer.data(msk + 376);
    const auto *msk_377 = buffer.data(msk + 377);
    const auto *msk_378 = buffer.data(msk + 378);
    const auto *msk_380 = buffer.data(msk + 380);
    const auto *msk_381 = buffer.data(msk + 381);

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pc_y, lsl0_228, lsl0_230, lsl0_231, \
                         lsk_181, lsk_182, lsk_183, lsl1_228, lsl1_230, lsl1_231, \
                         msk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_y[k] * lsl0_228[k]
                   + f_16 * lsk_181[k]
                   - f_14 * pc_y[k] * lsl1_228[k];

        t_364[k] = f_15 * lsk_182[k]
                   + f_3 * pc_y[k] * msk_290[k];

        t_365[k] = pa_y[k] * lsl0_230[k]
                   - f_14 * pc_y[k] * lsl1_230[k];

        t_366[k] = pa_y[k] * lsl0_231[k]
                   + f_17 * lsk_183[k]
                   - f_14 * pc_y[k] * lsl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pc_y, pc_z, lsl0_234, lsl0_235, \
                         lsk_147, lsk_185, lsk_186, lsl1_234, lsl1_235, msk_291, \
                         msk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * lsk_147[k]
                   + f_3 * pc_z[k] * msk_291[k];

        t_368[k] = f_15 * lsk_185[k]
                   + f_3 * pc_y[k] * msk_293[k];

        t_369[k] = pa_y[k] * lsl0_234[k]
                   - f_14 * pc_y[k] * lsl1_234[k];

        t_370[k] = pa_y[k] * lsl0_235[k]
                   + f_18 * lsk_186[k]
                   - f_14 * pc_y[k] * lsl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pc_y, pc_z, lsl0_237, lsl0_239, \
                         lsk_150, lsk_188, lsk_189, lsl1_237, lsl1_239, msk_294, \
                         msk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * lsk_150[k]
                   + f_3 * pc_z[k] * msk_294[k];

        t_372[k] = pa_y[k] * lsl0_237[k]
                   + f_16 * lsk_188[k]
                   - f_14 * pc_y[k] * lsl1_237[k];

        t_373[k] = f_15 * lsk_189[k]
                   + f_3 * pc_y[k] * msk_297[k];

        t_374[k] = pa_y[k] * lsl0_239[k]
                   - f_14 * pc_y[k] * lsl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_y, pc_y, pc_z, lsl0_240, lsl0_242, lsk_154, \
                         lsk_190, lsk_192, lsl1_240, lsl1_242, \
                         msk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_y[k] * lsl0_240[k]
                   + f_19 * lsk_190[k]
                   - f_14 * pc_y[k] * lsl1_240[k];

        t_376[k] = f_16 * lsk_154[k]
                   + f_3 * pc_z[k] * msk_298[k];

        t_377[k] = pa_y[k] * lsl0_242[k]
                   + f_17 * lsk_192[k]
                   - f_14 * pc_y[k] * lsl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_y, pc_y, lsl0_243, lsl0_245, lsl0_246, \
                         lsk_193, lsk_194, lsk_195, lsl1_243, lsl1_245, lsl1_246, \
                         msk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * lsl0_243[k]
                   + f_16 * lsk_193[k]
                   - f_14 * pc_y[k] * lsl1_243[k];

        t_379[k] = f_15 * lsk_194[k]
                   + f_3 * pc_y[k] * msk_302[k];

        t_380[k] = pa_y[k] * lsl0_245[k]
                   - f_14 * pc_y[k] * lsl1_245[k];

        t_381[k] = pa_y[k] * lsl0_246[k]
                   + f_20 * lsk_195[k]
                   - f_14 * pc_y[k] * lsl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pa_y, pc_y, pc_z, lsl0_248, lsl0_249, lsk_159, \
                         lsk_197, lsk_198, lsl1_248, lsl1_249, \
                         msk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * lsk_159[k]
                   + f_3 * pc_z[k] * msk_303[k];

        t_383[k] = pa_y[k] * lsl0_248[k]
                   + f_18 * lsk_197[k]
                   - f_14 * pc_y[k] * lsl1_248[k];

        t_384[k] = pa_y[k] * lsl0_249[k]
                   + f_17 * lsk_198[k]
                   - f_14 * pc_y[k] * lsl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pc_x, pc_y, lsl0_250, lsl0_252, \
                         lsk_199, lsk_200, lsk_316, lsl1_250, lsl1_252, msk_308, \
                         msk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pa_y[k] * lsl0_250[k]
                   + f_16 * lsk_199[k]
                   - f_14 * pc_y[k] * lsl1_250[k];

        t_386[k] = f_15 * lsk_200[k]
                   + f_3 * pc_y[k] * msk_308[k];

        t_387[k] = pa_y[k] * lsl0_252[k]
                   - f_14 * pc_y[k] * lsl1_252[k];

        t_388[k] = f_20 * lsk_316[k]
                   + f_3 * pc_x[k] * msk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, lsk_317, lsk_318, lsk_319, \
                         lsk_320, lsk_321, msk_317, msk_318, msk_319, msk_320, \
                         msk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_20 * lsk_317[k]
                   + f_3 * pc_x[k] * msk_317[k];

        t_390[k] = f_20 * lsk_318[k]
                   + f_3 * pc_x[k] * msk_318[k];

        t_391[k] = f_20 * lsk_319[k]
                   + f_3 * pc_x[k] * msk_319[k];

        t_392[k] = f_20 * lsk_320[k]
                   + f_3 * pc_x[k] * msk_320[k];

        t_393[k] = f_20 * lsk_321[k]
                   + f_3 * pc_x[k] * msk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, lsk_172, lsk_208, \
                         lsk_322, lsk_323, msi0_245, msi1_245, msk_316, msk_322, \
                         msk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_20 * lsk_322[k]
                   + f_3 * pc_x[k] * msk_322[k];

        t_395[k] = f_20 * lsk_323[k]
                   + f_3 * pc_x[k] * msk_323[k];

        t_396[k] = f_15 * lsk_208[k]
                   + f_1 * msi0_245[k]
                   - f_2 * msi1_245[k]
                   + f_3 * pc_y[k] * msk_316[k];

        t_397[k] = f_16 * lsk_172[k]
                   + f_3 * pc_z[k] * msk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, lsk_210, lsk_211, lsk_212, msi0_247, \
                         msi0_248, msi0_249, msi1_247, msi1_248, msi1_249, msk_318, msk_319, \
                         msk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * lsk_210[k]
                   + f_12 * msi0_247[k]
                   - f_13 * msi1_247[k]
                   + f_3 * pc_y[k] * msk_318[k];

        t_399[k] = f_15 * lsk_211[k]
                   + f_10 * msi0_248[k]
                   - f_11 * msi1_248[k]
                   + f_3 * pc_y[k] * msk_319[k];

        t_400[k] = f_15 * lsk_212[k]
                   + f_8 * msi0_249[k]
                   - f_9 * msi1_249[k]
                   + f_3 * pc_y[k] * msk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, lsk_213, lsk_214, lsk_215, msi0_250, \
                         msi0_251, msi1_250, msi1_251, msk_321, msk_322, \
                         msk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * lsk_213[k]
                   + f_6 * msi0_250[k]
                   - f_7 * msi1_250[k]
                   + f_3 * pc_y[k] * msk_321[k];

        t_402[k] = f_15 * lsk_214[k]
                   + f_4 * msi0_251[k]
                   - f_5 * msi1_251[k]
                   + f_3 * pc_y[k] * msk_322[k];

        t_403[k] = f_15 * lsk_215[k]
                   + f_3 * pc_y[k] * msk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_y, pc_x, pc_y, pc_z, lsl0_269, \
                         lsk_180, lsk_324, lsl1_269, msi0_252, msi1_252, \
                         msk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * lsl0_269[k]
                   - f_14 * pc_y[k] * lsl1_269[k];

        t_405[k] = f_20 * lsk_324[k]
                   + f_1 * msi0_252[k]
                   - f_2 * msi1_252[k]
                   + f_3 * pc_x[k] * msk_324[k];

        t_406[k] = f_3 * pc_y[k] * msk_324[k];

        t_407[k] = f_17 * lsk_180[k]
                   + f_3 * pc_z[k] * msk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, lsk_329, msi0_252, msi0_257, \
                         msi1_252, msi1_257, msk_325, msk_326, \
                         msk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_4 * msi0_252[k]
                   - f_5 * msi1_252[k]
                   + f_3 * pc_y[k] * msk_325[k];

        t_409[k] = f_3 * pc_y[k] * msk_326[k];

        t_410[k] = f_20 * lsk_329[k]
                   + f_12 * msi0_257[k]
                   - f_13 * msi1_257[k]
                   + f_3 * pc_x[k] * msk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, msi0_253, msi0_254, msi1_253, msi1_254, \
                         msk_327, msk_328, msk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * msi0_253[k]
                   - f_7 * msi1_253[k]
                   + f_3 * pc_y[k] * msk_327[k];

        t_412[k] = f_4 * msi0_254[k]
                   - f_5 * msi1_254[k]
                   + f_3 * pc_y[k] * msk_328[k];

        t_413[k] = f_3 * pc_y[k] * msk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, lsk_333, msi0_255, msi0_256, \
                         msi0_261, msi1_255, msi1_256, msi1_261, msk_330, msk_331, \
                         msk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_20 * lsk_333[k]
                   + f_10 * msi0_261[k]
                   - f_11 * msi1_261[k]
                   + f_3 * pc_x[k] * msk_333[k];

        t_415[k] = f_8 * msi0_255[k]
                   - f_9 * msi1_255[k]
                   + f_3 * pc_y[k] * msk_330[k];

        t_416[k] = f_6 * msi0_256[k]
                   - f_7 * msi1_256[k]
                   + f_3 * pc_y[k] * msk_331[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, lsk_338, msi0_257, msi0_266, \
                         msi1_257, msi1_266, msk_332, msk_333, \
                         msk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * msi0_257[k]
                   - f_5 * msi1_257[k]
                   + f_3 * pc_y[k] * msk_332[k];

        t_418[k] = f_3 * pc_y[k] * msk_333[k];

        t_419[k] = f_20 * lsk_338[k]
                   + f_8 * msi0_266[k]
                   - f_9 * msi1_266[k]
                   + f_3 * pc_x[k] * msk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_y, msi0_258, msi0_259, msi0_260, msi1_258, \
                         msi1_259, msi1_260, msk_334, msk_335, \
                         msk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * msi0_258[k]
                   - f_11 * msi1_258[k]
                   + f_3 * pc_y[k] * msk_334[k];

        t_421[k] = f_8 * msi0_259[k]
                   - f_9 * msi1_259[k]
                   + f_3 * pc_y[k] * msk_335[k];

        t_422[k] = f_6 * msi0_260[k]
                   - f_7 * msi1_260[k]
                   + f_3 * pc_y[k] * msk_336[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, lsk_344, msi0_261, msi0_272, \
                         msi1_261, msi1_272, msk_337, msk_338, \
                         msk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * msi0_261[k]
                   - f_5 * msi1_261[k]
                   + f_3 * pc_y[k] * msk_337[k];

        t_424[k] = f_3 * pc_y[k] * msk_338[k];

        t_425[k] = f_20 * lsk_344[k]
                   + f_6 * msi0_272[k]
                   - f_7 * msi1_272[k]
                   + f_3 * pc_x[k] * msk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_y, msi0_262, msi0_263, msi0_264, msi1_262, \
                         msi1_263, msi1_264, msk_339, msk_340, \
                         msk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * msi0_262[k]
                   - f_13 * msi1_262[k]
                   + f_3 * pc_y[k] * msk_339[k];

        t_427[k] = f_10 * msi0_263[k]
                   - f_11 * msi1_263[k]
                   + f_3 * pc_y[k] * msk_340[k];

        t_428[k] = f_8 * msi0_264[k]
                   - f_9 * msi1_264[k]
                   + f_3 * pc_y[k] * msk_341[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_y, msi0_265, msi0_266, msi1_265, msi1_266, \
                         msk_342, msk_343, msk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * msi0_265[k]
                   - f_7 * msi1_265[k]
                   + f_3 * pc_y[k] * msk_342[k];

        t_430[k] = f_4 * msi0_266[k]
                   - f_5 * msi1_266[k]
                   + f_3 * pc_y[k] * msk_343[k];

        t_431[k] = f_3 * pc_y[k] * msk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, lsk_351, lsk_352, lsk_353, lsk_354, \
                         msi0_279, msi1_279, msk_351, msk_352, msk_353, \
                         msk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_20 * lsk_351[k]
                   + f_4 * msi0_279[k]
                   - f_5 * msi1_279[k]
                   + f_3 * pc_x[k] * msk_351[k];

        t_433[k] = f_20 * lsk_352[k]
                   + f_3 * pc_x[k] * msk_352[k];

        t_434[k] = f_20 * lsk_353[k]
                   + f_3 * pc_x[k] * msk_353[k];

        t_435[k] = f_20 * lsk_354[k]
                   + f_3 * pc_x[k] * msk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_y, lsk_355, lsk_356, \
                         lsk_357, lsk_359, msk_351, msk_355, msk_356, msk_357, \
                         msk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_20 * lsk_355[k]
                   + f_3 * pc_x[k] * msk_355[k];

        t_437[k] = f_20 * lsk_356[k]
                   + f_3 * pc_x[k] * msk_356[k];

        t_438[k] = f_20 * lsk_357[k]
                   + f_3 * pc_x[k] * msk_357[k];

        t_439[k] = f_3 * pc_y[k] * msk_351[k];

        t_440[k] = f_20 * lsk_359[k]
                   + f_3 * pc_x[k] * msk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, msi0_273, msi0_274, msi0_275, msi1_273, \
                         msi1_274, msi1_275, msk_352, msk_353, \
                         msk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * msi0_273[k]
                   - f_2 * msi1_273[k]
                   + f_3 * pc_y[k] * msk_352[k];

        t_442[k] = f_22 * msi0_274[k]
                   - f_23 * msi1_274[k]
                   + f_3 * pc_y[k] * msk_353[k];

        t_443[k] = f_12 * msi0_275[k]
                   - f_13 * msi1_275[k]
                   + f_3 * pc_y[k] * msk_354[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, msi0_276, msi0_277, msi0_278, msi1_276, \
                         msi1_277, msi1_278, msk_355, msk_356, \
                         msk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_10 * msi0_276[k]
                   - f_11 * msi1_276[k]
                   + f_3 * pc_y[k] * msk_355[k];

        t_445[k] = f_8 * msi0_277[k]
                   - f_9 * msi1_277[k]
                   + f_3 * pc_y[k] * msk_356[k];

        t_446[k] = f_6 * msi0_278[k]
                   - f_7 * msi1_278[k]
                   + f_3 * pc_y[k] * msk_357[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pc_x, pc_y, pc_z, lsk_215, lsk_360, \
                         msi0_279, msi0_280, msi1_279, msi1_280, msk_358, msk_359, \
                         msk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * msi0_279[k]
                   - f_5 * msi1_279[k]
                   + f_3 * pc_y[k] * msk_358[k];

        t_448[k] = f_3 * pc_y[k] * msk_359[k];

        t_449[k] = f_17 * lsk_215[k]
                   + f_1 * msi0_279[k]
                   - f_2 * msi1_279[k]
                   + f_3 * pc_z[k] * msk_359[k];

        t_450[k] = f_19 * lsk_360[k]
                   + f_1 * msi0_280[k]
                   - f_2 * msi1_280[k]
                   + f_3 * pc_x[k] * msk_360[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, lsk_216, lsk_363, \
                         msi0_283, msi1_283, msk_360, msk_361, \
                         msk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * lsk_216[k]
                   + f_3 * pc_y[k] * msk_360[k];

        t_452[k] = f_3 * pc_z[k] * msk_360[k];

        t_453[k] = f_19 * lsk_363[k]
                   + f_12 * msi0_283[k]
                   - f_13 * msi1_283[k]
                   + f_3 * pc_x[k] * msk_363[k];

        t_454[k] = f_3 * pc_z[k] * msk_361[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_z, lsk_366, msi0_280, msi0_286, \
                         msi1_280, msi1_286, msk_362, msk_363, \
                         msk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_4 * msi0_280[k]
                   - f_5 * msi1_280[k]
                   + f_3 * pc_z[k] * msk_362[k];

        t_456[k] = f_19 * lsk_366[k]
                   + f_10 * msi0_286[k]
                   - f_11 * msi1_286[k]
                   + f_3 * pc_x[k] * msk_366[k];

        t_457[k] = f_3 * pc_z[k] * msk_363[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pc_x, pc_y, pc_z, lsk_221, lsk_370, \
                         msi0_282, msi0_290, msi1_282, msi1_290, msk_365, msk_366, \
                         msk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_18 * lsk_221[k]
                   + f_3 * pc_y[k] * msk_365[k];

        t_459[k] = f_6 * msi0_282[k]
                   - f_7 * msi1_282[k]
                   + f_3 * pc_z[k] * msk_365[k];

        t_460[k] = f_19 * lsk_370[k]
                   + f_8 * msi0_290[k]
                   - f_9 * msi1_290[k]
                   + f_3 * pc_x[k] * msk_370[k];

        t_461[k] = f_3 * pc_z[k] * msk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, lsk_225, msi0_283, msi0_285, \
                         msi1_283, msi1_285, msk_367, msk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_4 * msi0_283[k]
                   - f_5 * msi1_283[k]
                   + f_3 * pc_z[k] * msk_367[k];

        t_463[k] = f_18 * lsk_225[k]
                   + f_3 * pc_y[k] * msk_369[k];

        t_464[k] = f_8 * msi0_285[k]
                   - f_9 * msi1_285[k]
                   + f_3 * pc_z[k] * msk_369[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, lsk_375, msi0_286, msi0_295, \
                         msi1_286, msi1_295, msk_370, msk_371, \
                         msk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_19 * lsk_375[k]
                   + f_6 * msi0_295[k]
                   - f_7 * msi1_295[k]
                   + f_3 * pc_x[k] * msk_375[k];

        t_466[k] = f_3 * pc_z[k] * msk_370[k];

        t_467[k] = f_4 * msi0_286[k]
                   - f_5 * msi1_286[k]
                   + f_3 * pc_z[k] * msk_371[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_y, pc_z, lsk_230, msi0_287, msi0_289, \
                         msi1_287, msi1_289, msk_372, msk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_6 * msi0_287[k]
                   - f_7 * msi1_287[k]
                   + f_3 * pc_z[k] * msk_372[k];

        t_469[k] = f_18 * lsk_230[k]
                   + f_3 * pc_y[k] * msk_374[k];

        t_470[k] = f_10 * msi0_289[k]
                   - f_11 * msi1_289[k]
                   + f_3 * pc_z[k] * msk_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, lsk_381, msi0_290, msi0_301, \
                         msi1_290, msi1_301, msk_375, msk_376, \
                         msk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_19 * lsk_381[k]
                   + f_4 * msi0_301[k]
                   - f_5 * msi1_301[k]
                   + f_3 * pc_x[k] * msk_381[k];

        t_472[k] = f_3 * pc_z[k] * msk_375[k];

        t_473[k] = f_4 * msi0_290[k]
                   - f_5 * msi1_290[k]
                   + f_3 * pc_z[k] * msk_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_y, pc_z, lsk_236, msi0_291, msi0_292, \
                         msi0_294, msi1_291, msi1_292, msi1_294, msk_377, msk_378, \
                         msk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_6 * msi0_291[k]
                   - f_7 * msi1_291[k]
                   + f_3 * pc_z[k] * msk_377[k];

        t_475[k] = f_8 * msi0_292[k]
                   - f_9 * msi1_292[k]
                   + f_3 * pc_z[k] * msk_378[k];

        t_476[k] = f_18 * lsk_236[k]
                   + f_3 * pc_y[k] * msk_380[k];

        t_477[k] = f_12 * msi0_294[k]
                   - f_13 * msi1_294[k]
                   + f_3 * pc_z[k] * msk_380[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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

    const auto *lsl0_270 = buffer.data(lsl0 + 270);
    const auto *lsl0_273 = buffer.data(lsl0 + 273);
    const auto *lsl0_276 = buffer.data(lsl0 + 276);
    const auto *lsl0_280 = buffer.data(lsl0 + 280);
    const auto *lsl0_282 = buffer.data(lsl0 + 282);
    const auto *lsl0_285 = buffer.data(lsl0 + 285);
    const auto *lsl0_287 = buffer.data(lsl0 + 287);
    const auto *lsl0_288 = buffer.data(lsl0 + 288);
    const auto *lsl0_291 = buffer.data(lsl0 + 291);
    const auto *lsl0_293 = buffer.data(lsl0 + 293);
    const auto *lsl0_294 = buffer.data(lsl0 + 294);
    const auto *lsl0_295 = buffer.data(lsl0 + 295);
    const auto *lsl0_306 = buffer.data(lsl0 + 306);
    const auto *lsl0_405 = buffer.data(lsl0 + 405);

    const auto *lsk_216 = buffer.data(lsk + 216);
    const auto *lsk_219 = buffer.data(lsk + 219);
    const auto *lsk_222 = buffer.data(lsk + 222);
    const auto *lsk_223 = buffer.data(lsk + 223);
    const auto *lsk_226 = buffer.data(lsk + 226);
    const auto *lsk_227 = buffer.data(lsk + 227);
    const auto *lsk_228 = buffer.data(lsk + 228);
    const auto *lsk_231 = buffer.data(lsk + 231);
    const auto *lsk_232 = buffer.data(lsk + 232);
    const auto *lsk_233 = buffer.data(lsk + 233);
    const auto *lsk_234 = buffer.data(lsk + 234);
    const auto *lsk_244 = buffer.data(lsk + 244);
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
    const auto *lsk_280 = buffer.data(lsk + 280);
    const auto *lsk_282 = buffer.data(lsk + 282);
    const auto *lsk_283 = buffer.data(lsk + 283);
    const auto *lsk_284 = buffer.data(lsk + 284);
    const auto *lsk_285 = buffer.data(lsk + 285);
    const auto *lsk_286 = buffer.data(lsk + 286);
    const auto *lsk_287 = buffer.data(lsk + 287);
    const auto *lsk_288 = buffer.data(lsk + 288);
    const auto *lsk_290 = buffer.data(lsk + 290);
    const auto *lsk_293 = buffer.data(lsk + 293);
    const auto *lsk_297 = buffer.data(lsk + 297);
    const auto *lsk_302 = buffer.data(lsk + 302);
    const auto *lsk_308 = buffer.data(lsk + 308);
    const auto *lsk_316 = buffer.data(lsk + 316);
    const auto *lsk_318 = buffer.data(lsk + 318);
    const auto *lsk_319 = buffer.data(lsk + 319);
    const auto *lsk_320 = buffer.data(lsk + 320);
    const auto *lsk_321 = buffer.data(lsk + 321);
    const auto *lsk_322 = buffer.data(lsk + 322);
    const auto *lsk_323 = buffer.data(lsk + 323);
    const auto *lsk_324 = buffer.data(lsk + 324);
    const auto *lsk_388 = buffer.data(lsk + 388);
    const auto *lsk_390 = buffer.data(lsk + 390);
    const auto *lsk_391 = buffer.data(lsk + 391);
    const auto *lsk_392 = buffer.data(lsk + 392);
    const auto *lsk_393 = buffer.data(lsk + 393);
    const auto *lsk_394 = buffer.data(lsk + 394);
    const auto *lsk_395 = buffer.data(lsk + 395);
    const auto *lsk_401 = buffer.data(lsk + 401);
    const auto *lsk_405 = buffer.data(lsk + 405);
    const auto *lsk_410 = buffer.data(lsk + 410);
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

    const auto *lsl1_270 = buffer.data(lsl1 + 270);
    const auto *lsl1_273 = buffer.data(lsl1 + 273);
    const auto *lsl1_276 = buffer.data(lsl1 + 276);
    const auto *lsl1_280 = buffer.data(lsl1 + 280);
    const auto *lsl1_282 = buffer.data(lsl1 + 282);
    const auto *lsl1_285 = buffer.data(lsl1 + 285);
    const auto *lsl1_287 = buffer.data(lsl1 + 287);
    const auto *lsl1_288 = buffer.data(lsl1 + 288);
    const auto *lsl1_291 = buffer.data(lsl1 + 291);
    const auto *lsl1_293 = buffer.data(lsl1 + 293);
    const auto *lsl1_294 = buffer.data(lsl1 + 294);
    const auto *lsl1_295 = buffer.data(lsl1 + 295);
    const auto *lsl1_306 = buffer.data(lsl1 + 306);
    const auto *lsl1_405 = buffer.data(lsl1 + 405);

    const auto *msi0_301 = buffer.data(msi0 + 301);
    const auto *msi0_302 = buffer.data(msi0 + 302);
    const auto *msi0_303 = buffer.data(msi0 + 303);
    const auto *msi0_304 = buffer.data(msi0 + 304);
    const auto *msi0_305 = buffer.data(msi0 + 305);
    const auto *msi0_307 = buffer.data(msi0 + 307);
    const auto *msi0_313 = buffer.data(msi0 + 313);
    const auto *msi0_317 = buffer.data(msi0 + 317);
    const auto *msi0_322 = buffer.data(msi0 + 322);
    const auto *msi0_328 = buffer.data(msi0 + 328);
    const auto *msi0_331 = buffer.data(msi0 + 331);
    const auto *msi0_332 = buffer.data(msi0 + 332);
    const auto *msi0_333 = buffer.data(msi0 + 333);
    const auto *msi0_334 = buffer.data(msi0 + 334);
    const auto *msi0_335 = buffer.data(msi0 + 335);
    const auto *msi0_336 = buffer.data(msi0 + 336);
    const auto *msi0_339 = buffer.data(msi0 + 339);
    const auto *msi0_341 = buffer.data(msi0 + 341);
    const auto *msi0_342 = buffer.data(msi0 + 342);
    const auto *msi0_345 = buffer.data(msi0 + 345);
    const auto *msi0_346 = buffer.data(msi0 + 346);
    const auto *msi0_348 = buffer.data(msi0 + 348);
    const auto *msi0_350 = buffer.data(msi0 + 350);
    const auto *msi0_351 = buffer.data(msi0 + 351);
    const auto *msi0_353 = buffer.data(msi0 + 353);
    const auto *msi0_354 = buffer.data(msi0 + 354);
    const auto *msi0_356 = buffer.data(msi0 + 356);
    const auto *msi0_357 = buffer.data(msi0 + 357);
    const auto *msi0_359 = buffer.data(msi0 + 359);
    const auto *msi0_360 = buffer.data(msi0 + 360);
    const auto *msi0_361 = buffer.data(msi0 + 361);
    const auto *msi0_362 = buffer.data(msi0 + 362);
    const auto *msi0_363 = buffer.data(msi0 + 363);

    const auto *msi1_301 = buffer.data(msi1 + 301);
    const auto *msi1_302 = buffer.data(msi1 + 302);
    const auto *msi1_303 = buffer.data(msi1 + 303);
    const auto *msi1_304 = buffer.data(msi1 + 304);
    const auto *msi1_305 = buffer.data(msi1 + 305);
    const auto *msi1_307 = buffer.data(msi1 + 307);
    const auto *msi1_313 = buffer.data(msi1 + 313);
    const auto *msi1_317 = buffer.data(msi1 + 317);
    const auto *msi1_322 = buffer.data(msi1 + 322);
    const auto *msi1_328 = buffer.data(msi1 + 328);
    const auto *msi1_331 = buffer.data(msi1 + 331);
    const auto *msi1_332 = buffer.data(msi1 + 332);
    const auto *msi1_333 = buffer.data(msi1 + 333);
    const auto *msi1_334 = buffer.data(msi1 + 334);
    const auto *msi1_335 = buffer.data(msi1 + 335);
    const auto *msi1_336 = buffer.data(msi1 + 336);
    const auto *msi1_339 = buffer.data(msi1 + 339);
    const auto *msi1_341 = buffer.data(msi1 + 341);
    const auto *msi1_342 = buffer.data(msi1 + 342);
    const auto *msi1_345 = buffer.data(msi1 + 345);
    const auto *msi1_346 = buffer.data(msi1 + 346);
    const auto *msi1_348 = buffer.data(msi1 + 348);
    const auto *msi1_350 = buffer.data(msi1 + 350);
    const auto *msi1_351 = buffer.data(msi1 + 351);
    const auto *msi1_353 = buffer.data(msi1 + 353);
    const auto *msi1_354 = buffer.data(msi1 + 354);
    const auto *msi1_356 = buffer.data(msi1 + 356);
    const auto *msi1_357 = buffer.data(msi1 + 357);
    const auto *msi1_359 = buffer.data(msi1 + 359);
    const auto *msi1_360 = buffer.data(msi1 + 360);
    const auto *msi1_361 = buffer.data(msi1 + 361);
    const auto *msi1_362 = buffer.data(msi1 + 362);
    const auto *msi1_363 = buffer.data(msi1 + 363);

    const auto *msk_381 = buffer.data(msk + 381);
    const auto *msk_388 = buffer.data(msk + 388);
    const auto *msk_389 = buffer.data(msk + 389);
    const auto *msk_390 = buffer.data(msk + 390);
    const auto *msk_391 = buffer.data(msk + 391);
    const auto *msk_392 = buffer.data(msk + 392);
    const auto *msk_393 = buffer.data(msk + 393);
    const auto *msk_394 = buffer.data(msk + 394);
    const auto *msk_395 = buffer.data(msk + 395);
    const auto *msk_396 = buffer.data(msk + 396);
    const auto *msk_398 = buffer.data(msk + 398);
    const auto *msk_399 = buffer.data(msk + 399);
    const auto *msk_401 = buffer.data(msk + 401);
    const auto *msk_402 = buffer.data(msk + 402);
    const auto *msk_405 = buffer.data(msk + 405);
    const auto *msk_406 = buffer.data(msk + 406);
    const auto *msk_410 = buffer.data(msk + 410);
    const auto *msk_411 = buffer.data(msk + 411);
    const auto *msk_416 = buffer.data(msk + 416);
    const auto *msk_423 = buffer.data(msk + 423);
    const auto *msk_424 = buffer.data(msk + 424);
    const auto *msk_425 = buffer.data(msk + 425);
    const auto *msk_426 = buffer.data(msk + 426);
    const auto *msk_427 = buffer.data(msk + 427);
    const auto *msk_428 = buffer.data(msk + 428);
    const auto *msk_429 = buffer.data(msk + 429);
    const auto *msk_430 = buffer.data(msk + 430);
    const auto *msk_431 = buffer.data(msk + 431);
    const auto *msk_432 = buffer.data(msk + 432);
    const auto *msk_434 = buffer.data(msk + 434);
    const auto *msk_435 = buffer.data(msk + 435);
    const auto *msk_437 = buffer.data(msk + 437);
    const auto *msk_438 = buffer.data(msk + 438);
    const auto *msk_441 = buffer.data(msk + 441);
    const auto *msk_442 = buffer.data(msk + 442);
    const auto *msk_444 = buffer.data(msk + 444);
    const auto *msk_446 = buffer.data(msk + 446);
    const auto *msk_447 = buffer.data(msk + 447);
    const auto *msk_449 = buffer.data(msk + 449);
    const auto *msk_450 = buffer.data(msk + 450);
    const auto *msk_452 = buffer.data(msk + 452);
    const auto *msk_453 = buffer.data(msk + 453);
    const auto *msk_455 = buffer.data(msk + 455);
    const auto *msk_456 = buffer.data(msk + 456);
    const auto *msk_457 = buffer.data(msk + 457);
    const auto *msk_459 = buffer.data(msk + 459);
    const auto *msk_460 = buffer.data(msk + 460);
    const auto *msk_461 = buffer.data(msk + 461);
    const auto *msk_462 = buffer.data(msk + 462);
    const auto *msk_463 = buffer.data(msk + 463);
    const auto *msk_464 = buffer.data(msk + 464);
    const auto *msk_465 = buffer.data(msk + 465);
    const auto *msk_466 = buffer.data(msk + 466);
    const auto *msk_467 = buffer.data(msk + 467);
    const auto *msk_468 = buffer.data(msk + 468);

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pc_x, pc_z, lsk_388, lsk_390, \
                         lsk_391, lsk_392, msk_381, msk_388, msk_390, msk_391, \
                         msk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_19 * lsk_388[k]
                   + f_3 * pc_x[k] * msk_388[k];

        t_479[k] = f_3 * pc_z[k] * msk_381[k];

        t_480[k] = f_19 * lsk_390[k]
                   + f_3 * pc_x[k] * msk_390[k];

        t_481[k] = f_19 * lsk_391[k]
                   + f_3 * pc_x[k] * msk_391[k];

        t_482[k] = f_19 * lsk_392[k]
                   + f_3 * pc_x[k] * msk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pc_x, pc_y, lsk_244, lsk_393, lsk_394, \
                         lsk_395, msi0_301, msi1_301, msk_388, msk_393, msk_394, \
                         msk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_19 * lsk_393[k]
                   + f_3 * pc_x[k] * msk_393[k];

        t_484[k] = f_19 * lsk_394[k]
                   + f_3 * pc_x[k] * msk_394[k];

        t_485[k] = f_19 * lsk_395[k]
                   + f_3 * pc_x[k] * msk_395[k];

        t_486[k] = f_18 * lsk_244[k]
                   + f_1 * msi0_301[k]
                   - f_2 * msi1_301[k]
                   + f_3 * pc_y[k] * msk_388[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pc_z, msi0_301, msi0_302, msi0_303, \
                         msi1_301, msi1_302, msi1_303, msk_388, msk_389, msk_390, \
                         msk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_3 * pc_z[k] * msk_388[k];

        t_488[k] = f_4 * msi0_301[k]
                   - f_5 * msi1_301[k]
                   + f_3 * pc_z[k] * msk_389[k];

        t_489[k] = f_6 * msi0_302[k]
                   - f_7 * msi1_302[k]
                   + f_3 * pc_z[k] * msk_390[k];

        t_490[k] = f_8 * msi0_303[k]
                   - f_9 * msi1_303[k]
                   + f_3 * pc_z[k] * msk_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_y, pc_z, lsk_251, msi0_304, msi0_305, \
                         msi0_307, msi1_304, msi1_305, msi1_307, msk_392, msk_393, \
                         msk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * msi0_304[k]
                   - f_11 * msi1_304[k]
                   + f_3 * pc_z[k] * msk_392[k];

        t_492[k] = f_12 * msi0_305[k]
                   - f_13 * msi1_305[k]
                   + f_3 * pc_z[k] * msk_393[k];

        t_493[k] = f_18 * lsk_251[k]
                   + f_3 * pc_y[k] * msk_395[k];

        t_494[k] = f_1 * msi0_307[k]
                   - f_2 * msi1_307[k]
                   + f_3 * pc_z[k] * msk_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_z, pc_y, pc_z, lsl0_270, lsl0_273, \
                         lsk_216, lsk_252, lsl1_270, lsl1_273, \
                         msk_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * lsl0_270[k]
                   - f_14 * pc_z[k] * lsl1_270[k];

        t_496[k] = f_17 * lsk_252[k]
                   + f_3 * pc_y[k] * msk_396[k];

        t_497[k] = f_15 * lsk_216[k]
                   + f_3 * pc_z[k] * msk_396[k];

        t_498[k] = pa_z[k] * lsl0_273[k]
                   - f_14 * pc_z[k] * lsl1_273[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_z, pc_x, pc_y, pc_z, lsl0_276, lsk_254, \
                         lsk_401, lsl1_276, msi0_313, msi1_313, msk_398, \
                         msk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_17 * lsk_254[k]
                   + f_3 * pc_y[k] * msk_398[k];

        t_500[k] = f_19 * lsk_401[k]
                   + f_12 * msi0_313[k]
                   - f_13 * msi1_313[k]
                   + f_3 * pc_x[k] * msk_401[k];

        t_501[k] = pa_z[k] * lsl0_276[k]
                   - f_14 * pc_z[k] * lsl1_276[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, lsk_219, lsk_257, lsk_405, \
                         msi0_317, msi1_317, msk_399, msk_401, \
                         msk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_15 * lsk_219[k]
                   + f_3 * pc_z[k] * msk_399[k];

        t_503[k] = f_17 * lsk_257[k]
                   + f_3 * pc_y[k] * msk_401[k];

        t_504[k] = f_19 * lsk_405[k]
                   + f_10 * msi0_317[k]
                   - f_11 * msi1_317[k]
                   + f_3 * pc_x[k] * msk_405[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pa_z, pc_y, pc_z, lsl0_280, lsl0_282, \
                         lsk_222, lsk_223, lsk_261, lsl1_280, lsl1_282, msk_402, \
                         msk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_z[k] * lsl0_280[k]
                   - f_14 * pc_z[k] * lsl1_280[k];

        t_506[k] = f_15 * lsk_222[k]
                   + f_3 * pc_z[k] * msk_402[k];

        t_507[k] = pa_z[k] * lsl0_282[k]
                   + f_16 * lsk_223[k]
                   - f_14 * pc_z[k] * lsl1_282[k];

        t_508[k] = f_17 * lsk_261[k]
                   + f_3 * pc_y[k] * msk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_z, pc_x, pc_z, lsl0_285, lsk_226, lsk_410, \
                         lsl1_285, msi0_322, msi1_322, msk_406, \
                         msk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_19 * lsk_410[k]
                   + f_8 * msi0_322[k]
                   - f_9 * msi1_322[k]
                   + f_3 * pc_x[k] * msk_410[k];

        t_510[k] = pa_z[k] * lsl0_285[k]
                   - f_14 * pc_z[k] * lsl1_285[k];

        t_511[k] = f_15 * lsk_226[k]
                   + f_3 * pc_z[k] * msk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_z, pc_y, pc_z, lsl0_287, lsl0_288, lsk_227, \
                         lsk_228, lsk_266, lsl1_287, lsl1_288, \
                         msk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pa_z[k] * lsl0_287[k]
                   + f_16 * lsk_227[k]
                   - f_14 * pc_z[k] * lsl1_287[k];

        t_513[k] = pa_z[k] * lsl0_288[k]
                   + f_17 * lsk_228[k]
                   - f_14 * pc_z[k] * lsl1_288[k];

        t_514[k] = f_17 * lsk_266[k]
                   + f_3 * pc_y[k] * msk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pa_z, pc_x, pc_z, lsl0_291, lsk_231, lsk_416, \
                         lsl1_291, msi0_328, msi1_328, msk_411, \
                         msk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_19 * lsk_416[k]
                   + f_6 * msi0_328[k]
                   - f_7 * msi1_328[k]
                   + f_3 * pc_x[k] * msk_416[k];

        t_516[k] = pa_z[k] * lsl0_291[k]
                   - f_14 * pc_z[k] * lsl1_291[k];

        t_517[k] = f_15 * lsk_231[k]
                   + f_3 * pc_z[k] * msk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pa_z, pc_z, lsl0_293, lsl0_294, lsl0_295, \
                         lsk_232, lsk_233, lsk_234, lsl1_293, lsl1_294, \
                         lsl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pa_z[k] * lsl0_293[k]
                   + f_16 * lsk_232[k]
                   - f_14 * pc_z[k] * lsl1_293[k];

        t_519[k] = pa_z[k] * lsl0_294[k]
                   + f_17 * lsk_233[k]
                   - f_14 * pc_z[k] * lsl1_294[k];

        t_520[k] = pa_z[k] * lsl0_295[k]
                   + f_18 * lsk_234[k]
                   - f_14 * pc_z[k] * lsl1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, lsk_272, lsk_423, lsk_424, \
                         lsk_425, msi0_335, msi1_335, msk_416, msk_423, msk_424, \
                         msk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * lsk_272[k]
                   + f_3 * pc_y[k] * msk_416[k];

        t_522[k] = f_19 * lsk_423[k]
                   + f_4 * msi0_335[k]
                   - f_5 * msi1_335[k]
                   + f_3 * pc_x[k] * msk_423[k];

        t_523[k] = f_19 * lsk_424[k]
                   + f_3 * pc_x[k] * msk_424[k];

        t_524[k] = f_19 * lsk_425[k]
                   + f_3 * pc_x[k] * msk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, lsk_426, lsk_427, lsk_428, \
                         lsk_429, lsk_430, msk_426, msk_427, msk_428, msk_429, \
                         msk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_19 * lsk_426[k]
                   + f_3 * pc_x[k] * msk_426[k];

        t_526[k] = f_19 * lsk_427[k]
                   + f_3 * pc_x[k] * msk_427[k];

        t_527[k] = f_19 * lsk_428[k]
                   + f_3 * pc_x[k] * msk_428[k];

        t_528[k] = f_19 * lsk_429[k]
                   + f_3 * pc_x[k] * msk_429[k];

        t_529[k] = f_19 * lsk_430[k]
                   + f_3 * pc_x[k] * msk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_z, pc_x, pc_z, lsl0_306, lsk_244, lsk_431, \
                         lsl1_306, msk_424, msk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_19 * lsk_431[k]
                   + f_3 * pc_x[k] * msk_431[k];

        t_531[k] = pa_z[k] * lsl0_306[k]
                   - f_14 * pc_z[k] * lsl1_306[k];

        t_532[k] = f_15 * lsk_244[k]
                   + f_3 * pc_z[k] * msk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, lsk_282, lsk_283, lsk_284, msi0_331, \
                         msi0_332, msi0_333, msi1_331, msi1_332, msi1_333, msk_426, msk_427, \
                         msk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * lsk_282[k]
                   + f_12 * msi0_331[k]
                   - f_13 * msi1_331[k]
                   + f_3 * pc_y[k] * msk_426[k];

        t_534[k] = f_17 * lsk_283[k]
                   + f_10 * msi0_332[k]
                   - f_11 * msi1_332[k]
                   + f_3 * pc_y[k] * msk_427[k];

        t_535[k] = f_17 * lsk_284[k]
                   + f_8 * msi0_333[k]
                   - f_9 * msi1_333[k]
                   + f_3 * pc_y[k] * msk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, lsk_285, lsk_286, lsk_287, msi0_334, \
                         msi0_335, msi1_334, msi1_335, msk_429, msk_430, \
                         msk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * lsk_285[k]
                   + f_6 * msi0_334[k]
                   - f_7 * msi1_334[k]
                   + f_3 * pc_y[k] * msk_429[k];

        t_537[k] = f_17 * lsk_286[k]
                   + f_4 * msi0_335[k]
                   - f_5 * msi1_335[k]
                   + f_3 * pc_y[k] * msk_430[k];

        t_538[k] = f_17 * lsk_287[k]
                   + f_3 * pc_y[k] * msk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, lsk_251, lsk_288, lsk_432, \
                         msi0_335, msi0_336, msi1_335, msi1_336, msk_431, \
                         msk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * lsk_251[k]
                   + f_1 * msi0_335[k]
                   - f_2 * msi1_335[k]
                   + f_3 * pc_z[k] * msk_431[k];

        t_540[k] = f_19 * lsk_432[k]
                   + f_1 * msi0_336[k]
                   - f_2 * msi1_336[k]
                   + f_3 * pc_x[k] * msk_432[k];

        t_541[k] = f_16 * lsk_288[k]
                   + f_3 * pc_y[k] * msk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, lsk_252, lsk_290, lsk_435, \
                         msi0_339, msi1_339, msk_432, msk_434, \
                         msk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * lsk_252[k]
                   + f_3 * pc_z[k] * msk_432[k];

        t_543[k] = f_19 * lsk_435[k]
                   + f_12 * msi0_339[k]
                   - f_13 * msi1_339[k]
                   + f_3 * pc_x[k] * msk_435[k];

        t_544[k] = f_16 * lsk_290[k]
                   + f_3 * pc_y[k] * msk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, lsk_255, lsk_437, lsk_438, msi0_341, \
                         msi0_342, msi1_341, msi1_342, msk_435, msk_437, \
                         msk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_19 * lsk_437[k]
                   + f_12 * msi0_341[k]
                   - f_13 * msi1_341[k]
                   + f_3 * pc_x[k] * msk_437[k];

        t_546[k] = f_19 * lsk_438[k]
                   + f_10 * msi0_342[k]
                   - f_11 * msi1_342[k]
                   + f_3 * pc_x[k] * msk_438[k];

        t_547[k] = f_16 * lsk_255[k]
                   + f_3 * pc_z[k] * msk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, lsk_293, lsk_441, lsk_442, msi0_345, \
                         msi0_346, msi1_345, msi1_346, msk_437, msk_441, \
                         msk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * lsk_293[k]
                   + f_3 * pc_y[k] * msk_437[k];

        t_549[k] = f_19 * lsk_441[k]
                   + f_10 * msi0_345[k]
                   - f_11 * msi1_345[k]
                   + f_3 * pc_x[k] * msk_441[k];

        t_550[k] = f_19 * lsk_442[k]
                   + f_8 * msi0_346[k]
                   - f_9 * msi1_346[k]
                   + f_3 * pc_x[k] * msk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, lsk_258, lsk_297, lsk_444, \
                         msi0_348, msi1_348, msk_438, msk_441, \
                         msk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * lsk_258[k]
                   + f_3 * pc_z[k] * msk_438[k];

        t_552[k] = f_19 * lsk_444[k]
                   + f_8 * msi0_348[k]
                   - f_9 * msi1_348[k]
                   + f_3 * pc_x[k] * msk_444[k];

        t_553[k] = f_16 * lsk_297[k]
                   + f_3 * pc_y[k] * msk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, lsk_262, lsk_446, lsk_447, msi0_350, \
                         msi0_351, msi1_350, msi1_351, msk_442, msk_446, \
                         msk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_19 * lsk_446[k]
                   + f_8 * msi0_350[k]
                   - f_9 * msi1_350[k]
                   + f_3 * pc_x[k] * msk_446[k];

        t_555[k] = f_19 * lsk_447[k]
                   + f_6 * msi0_351[k]
                   - f_7 * msi1_351[k]
                   + f_3 * pc_x[k] * msk_447[k];

        t_556[k] = f_16 * lsk_262[k]
                   + f_3 * pc_z[k] * msk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, lsk_302, lsk_449, lsk_450, msi0_353, \
                         msi0_354, msi1_353, msi1_354, msk_446, msk_449, \
                         msk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_19 * lsk_449[k]
                   + f_6 * msi0_353[k]
                   - f_7 * msi1_353[k]
                   + f_3 * pc_x[k] * msk_449[k];

        t_558[k] = f_19 * lsk_450[k]
                   + f_6 * msi0_354[k]
                   - f_7 * msi1_354[k]
                   + f_3 * pc_x[k] * msk_450[k];

        t_559[k] = f_16 * lsk_302[k]
                   + f_3 * pc_y[k] * msk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, lsk_267, lsk_452, lsk_453, msi0_356, \
                         msi0_357, msi1_356, msi1_357, msk_447, msk_452, \
                         msk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_19 * lsk_452[k]
                   + f_6 * msi0_356[k]
                   - f_7 * msi1_356[k]
                   + f_3 * pc_x[k] * msk_452[k];

        t_561[k] = f_19 * lsk_453[k]
                   + f_4 * msi0_357[k]
                   - f_5 * msi1_357[k]
                   + f_3 * pc_x[k] * msk_453[k];

        t_562[k] = f_16 * lsk_267[k]
                   + f_3 * pc_z[k] * msk_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, lsk_455, lsk_456, lsk_457, msi0_359, \
                         msi0_360, msi0_361, msi1_359, msi1_360, msi1_361, msk_455, msk_456, \
                         msk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_19 * lsk_455[k]
                   + f_4 * msi0_359[k]
                   - f_5 * msi1_359[k]
                   + f_3 * pc_x[k] * msk_455[k];

        t_564[k] = f_19 * lsk_456[k]
                   + f_4 * msi0_360[k]
                   - f_5 * msi1_360[k]
                   + f_3 * pc_x[k] * msk_456[k];

        t_565[k] = f_19 * lsk_457[k]
                   + f_4 * msi0_361[k]
                   - f_5 * msi1_361[k]
                   + f_3 * pc_x[k] * msk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, lsk_308, lsk_459, lsk_460, \
                         lsk_461, msi0_363, msi1_363, msk_452, msk_459, msk_460, \
                         msk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * lsk_308[k]
                   + f_3 * pc_y[k] * msk_452[k];

        t_567[k] = f_19 * lsk_459[k]
                   + f_4 * msi0_363[k]
                   - f_5 * msi1_363[k]
                   + f_3 * pc_x[k] * msk_459[k];

        t_568[k] = f_19 * lsk_460[k]
                   + f_3 * pc_x[k] * msk_460[k];

        t_569[k] = f_19 * lsk_461[k]
                   + f_3 * pc_x[k] * msk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, lsk_462, lsk_463, lsk_464, \
                         lsk_465, lsk_466, msk_462, msk_463, msk_464, msk_465, \
                         msk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_19 * lsk_462[k]
                   + f_3 * pc_x[k] * msk_462[k];

        t_571[k] = f_19 * lsk_463[k]
                   + f_3 * pc_x[k] * msk_463[k];

        t_572[k] = f_19 * lsk_464[k]
                   + f_3 * pc_x[k] * msk_464[k];

        t_573[k] = f_19 * lsk_465[k]
                   + f_3 * pc_x[k] * msk_465[k];

        t_574[k] = f_19 * lsk_466[k]
                   + f_3 * pc_x[k] * msk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, lsk_280, lsk_316, lsk_467, \
                         msi0_357, msi1_357, msk_460, msk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_19 * lsk_467[k]
                   + f_3 * pc_x[k] * msk_467[k];

        t_576[k] = f_16 * lsk_316[k]
                   + f_1 * msi0_357[k]
                   - f_2 * msi1_357[k]
                   + f_3 * pc_y[k] * msk_460[k];

        t_577[k] = f_16 * lsk_280[k]
                   + f_3 * pc_z[k] * msk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, lsk_318, lsk_319, lsk_320, msi0_359, \
                         msi0_360, msi0_361, msi1_359, msi1_360, msi1_361, msk_462, msk_463, \
                         msk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * lsk_318[k]
                   + f_12 * msi0_359[k]
                   - f_13 * msi1_359[k]
                   + f_3 * pc_y[k] * msk_462[k];

        t_579[k] = f_16 * lsk_319[k]
                   + f_10 * msi0_360[k]
                   - f_11 * msi1_360[k]
                   + f_3 * pc_y[k] * msk_463[k];

        t_580[k] = f_16 * lsk_320[k]
                   + f_8 * msi0_361[k]
                   - f_9 * msi1_361[k]
                   + f_3 * pc_y[k] * msk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, lsk_321, lsk_322, lsk_323, msi0_362, \
                         msi0_363, msi1_362, msi1_363, msk_465, msk_466, \
                         msk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * lsk_321[k]
                   + f_6 * msi0_362[k]
                   - f_7 * msi1_362[k]
                   + f_3 * pc_y[k] * msk_465[k];

        t_582[k] = f_16 * lsk_322[k]
                   + f_4 * msi0_363[k]
                   - f_5 * msi1_363[k]
                   + f_3 * pc_y[k] * msk_466[k];

        t_583[k] = f_16 * lsk_323[k]
                   + f_3 * pc_y[k] * msk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_y, pc_y, pc_z, lsl0_405, lsk_287, \
                         lsk_288, lsk_324, lsl1_405, msi0_363, msi1_363, msk_467, \
                         msk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * lsk_287[k]
                   + f_1 * msi0_363[k]
                   - f_2 * msi1_363[k]
                   + f_3 * pc_z[k] * msk_467[k];

        t_585[k] = pa_y[k] * lsl0_405[k]
                   - f_14 * pc_y[k] * lsl1_405[k];

        t_586[k] = f_15 * lsk_324[k]
                   + f_3 * pc_y[k] * msk_468[k];

        t_587[k] = f_17 * lsk_288[k]
                   + f_3 * pc_z[k] * msk_468[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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

    const auto *lsl0_408 = buffer.data(lsl0 + 408);
    const auto *lsl0_410 = buffer.data(lsl0 + 410);
    const auto *lsl0_411 = buffer.data(lsl0 + 411);
    const auto *lsl0_414 = buffer.data(lsl0 + 414);
    const auto *lsl0_415 = buffer.data(lsl0 + 415);
    const auto *lsl0_417 = buffer.data(lsl0 + 417);
    const auto *lsl0_419 = buffer.data(lsl0 + 419);
    const auto *lsl0_420 = buffer.data(lsl0 + 420);
    const auto *lsl0_422 = buffer.data(lsl0 + 422);
    const auto *lsl0_423 = buffer.data(lsl0 + 423);
    const auto *lsl0_425 = buffer.data(lsl0 + 425);
    const auto *lsl0_426 = buffer.data(lsl0 + 426);
    const auto *lsl0_428 = buffer.data(lsl0 + 428);
    const auto *lsl0_429 = buffer.data(lsl0 + 429);
    const auto *lsl0_430 = buffer.data(lsl0 + 430);
    const auto *lsl0_432 = buffer.data(lsl0 + 432);
    const auto *lsl0_449 = buffer.data(lsl0 + 449);

    const auto *lsk_291 = buffer.data(lsk + 291);
    const auto *lsk_294 = buffer.data(lsk + 294);
    const auto *lsk_298 = buffer.data(lsk + 298);
    const auto *lsk_303 = buffer.data(lsk + 303);
    const auto *lsk_316 = buffer.data(lsk + 316);
    const auto *lsk_324 = buffer.data(lsk + 324);
    const auto *lsk_325 = buffer.data(lsk + 325);
    const auto *lsk_326 = buffer.data(lsk + 326);
    const auto *lsk_327 = buffer.data(lsk + 327);
    const auto *lsk_329 = buffer.data(lsk + 329);
    const auto *lsk_330 = buffer.data(lsk + 330);
    const auto *lsk_332 = buffer.data(lsk + 332);
    const auto *lsk_333 = buffer.data(lsk + 333);
    const auto *lsk_334 = buffer.data(lsk + 334);
    const auto *lsk_336 = buffer.data(lsk + 336);
    const auto *lsk_337 = buffer.data(lsk + 337);
    const auto *lsk_338 = buffer.data(lsk + 338);
    const auto *lsk_339 = buffer.data(lsk + 339);
    const auto *lsk_341 = buffer.data(lsk + 341);
    const auto *lsk_342 = buffer.data(lsk + 342);
    const auto *lsk_343 = buffer.data(lsk + 343);
    const auto *lsk_344 = buffer.data(lsk + 344);
    const auto *lsk_352 = buffer.data(lsk + 352);
    const auto *lsk_354 = buffer.data(lsk + 354);
    const auto *lsk_355 = buffer.data(lsk + 355);
    const auto *lsk_356 = buffer.data(lsk + 356);
    const auto *lsk_357 = buffer.data(lsk + 357);
    const auto *lsk_358 = buffer.data(lsk + 358);
    const auto *lsk_359 = buffer.data(lsk + 359);
    const auto *lsk_360 = buffer.data(lsk + 360);
    const auto *lsk_365 = buffer.data(lsk + 365);
    const auto *lsk_369 = buffer.data(lsk + 369);
    const auto *lsk_374 = buffer.data(lsk + 374);
    const auto *lsk_380 = buffer.data(lsk + 380);
    const auto *lsk_496 = buffer.data(lsk + 496);
    const auto *lsk_497 = buffer.data(lsk + 497);
    const auto *lsk_498 = buffer.data(lsk + 498);
    const auto *lsk_499 = buffer.data(lsk + 499);
    const auto *lsk_500 = buffer.data(lsk + 500);
    const auto *lsk_501 = buffer.data(lsk + 501);
    const auto *lsk_502 = buffer.data(lsk + 502);
    const auto *lsk_503 = buffer.data(lsk + 503);
    const auto *lsk_504 = buffer.data(lsk + 504);
    const auto *lsk_509 = buffer.data(lsk + 509);
    const auto *lsk_513 = buffer.data(lsk + 513);
    const auto *lsk_518 = buffer.data(lsk + 518);
    const auto *lsk_524 = buffer.data(lsk + 524);
    const auto *lsk_531 = buffer.data(lsk + 531);
    const auto *lsk_532 = buffer.data(lsk + 532);
    const auto *lsk_533 = buffer.data(lsk + 533);
    const auto *lsk_534 = buffer.data(lsk + 534);
    const auto *lsk_535 = buffer.data(lsk + 535);
    const auto *lsk_536 = buffer.data(lsk + 536);
    const auto *lsk_537 = buffer.data(lsk + 537);
    const auto *lsk_539 = buffer.data(lsk + 539);
    const auto *lsk_540 = buffer.data(lsk + 540);
    const auto *lsk_543 = buffer.data(lsk + 543);
    const auto *lsk_546 = buffer.data(lsk + 546);
    const auto *lsk_550 = buffer.data(lsk + 550);
    const auto *lsk_555 = buffer.data(lsk + 555);
    const auto *lsk_561 = buffer.data(lsk + 561);

    const auto *lsl1_408 = buffer.data(lsl1 + 408);
    const auto *lsl1_410 = buffer.data(lsl1 + 410);
    const auto *lsl1_411 = buffer.data(lsl1 + 411);
    const auto *lsl1_414 = buffer.data(lsl1 + 414);
    const auto *lsl1_415 = buffer.data(lsl1 + 415);
    const auto *lsl1_417 = buffer.data(lsl1 + 417);
    const auto *lsl1_419 = buffer.data(lsl1 + 419);
    const auto *lsl1_420 = buffer.data(lsl1 + 420);
    const auto *lsl1_422 = buffer.data(lsl1 + 422);
    const auto *lsl1_423 = buffer.data(lsl1 + 423);
    const auto *lsl1_425 = buffer.data(lsl1 + 425);
    const auto *lsl1_426 = buffer.data(lsl1 + 426);
    const auto *lsl1_428 = buffer.data(lsl1 + 428);
    const auto *lsl1_429 = buffer.data(lsl1 + 429);
    const auto *lsl1_430 = buffer.data(lsl1 + 430);
    const auto *lsl1_432 = buffer.data(lsl1 + 432);
    const auto *lsl1_449 = buffer.data(lsl1 + 449);

    const auto *msi0_385 = buffer.data(msi0 + 385);
    const auto *msi0_387 = buffer.data(msi0 + 387);
    const auto *msi0_388 = buffer.data(msi0 + 388);
    const auto *msi0_389 = buffer.data(msi0 + 389);
    const auto *msi0_390 = buffer.data(msi0 + 390);
    const auto *msi0_391 = buffer.data(msi0 + 391);
    const auto *msi0_392 = buffer.data(msi0 + 392);
    const auto *msi0_393 = buffer.data(msi0 + 393);
    const auto *msi0_394 = buffer.data(msi0 + 394);
    const auto *msi0_395 = buffer.data(msi0 + 395);
    const auto *msi0_396 = buffer.data(msi0 + 396);
    const auto *msi0_397 = buffer.data(msi0 + 397);
    const auto *msi0_398 = buffer.data(msi0 + 398);
    const auto *msi0_399 = buffer.data(msi0 + 399);
    const auto *msi0_400 = buffer.data(msi0 + 400);
    const auto *msi0_401 = buffer.data(msi0 + 401);
    const auto *msi0_402 = buffer.data(msi0 + 402);
    const auto *msi0_403 = buffer.data(msi0 + 403);
    const auto *msi0_404 = buffer.data(msi0 + 404);
    const auto *msi0_405 = buffer.data(msi0 + 405);
    const auto *msi0_406 = buffer.data(msi0 + 406);
    const auto *msi0_412 = buffer.data(msi0 + 412);
    const auto *msi0_413 = buffer.data(msi0 + 413);
    const auto *msi0_414 = buffer.data(msi0 + 414);
    const auto *msi0_415 = buffer.data(msi0 + 415);
    const auto *msi0_416 = buffer.data(msi0 + 416);
    const auto *msi0_417 = buffer.data(msi0 + 417);
    const auto *msi0_418 = buffer.data(msi0 + 418);
    const auto *msi0_419 = buffer.data(msi0 + 419);
    const auto *msi0_420 = buffer.data(msi0 + 420);
    const auto *msi0_422 = buffer.data(msi0 + 422);
    const auto *msi0_423 = buffer.data(msi0 + 423);
    const auto *msi0_425 = buffer.data(msi0 + 425);
    const auto *msi0_426 = buffer.data(msi0 + 426);
    const auto *msi0_427 = buffer.data(msi0 + 427);
    const auto *msi0_429 = buffer.data(msi0 + 429);
    const auto *msi0_430 = buffer.data(msi0 + 430);
    const auto *msi0_431 = buffer.data(msi0 + 431);
    const auto *msi0_432 = buffer.data(msi0 + 432);
    const auto *msi0_434 = buffer.data(msi0 + 434);
    const auto *msi0_435 = buffer.data(msi0 + 435);
    const auto *msi0_441 = buffer.data(msi0 + 441);

    const auto *msi1_385 = buffer.data(msi1 + 385);
    const auto *msi1_387 = buffer.data(msi1 + 387);
    const auto *msi1_388 = buffer.data(msi1 + 388);
    const auto *msi1_389 = buffer.data(msi1 + 389);
    const auto *msi1_390 = buffer.data(msi1 + 390);
    const auto *msi1_391 = buffer.data(msi1 + 391);
    const auto *msi1_392 = buffer.data(msi1 + 392);
    const auto *msi1_393 = buffer.data(msi1 + 393);
    const auto *msi1_394 = buffer.data(msi1 + 394);
    const auto *msi1_395 = buffer.data(msi1 + 395);
    const auto *msi1_396 = buffer.data(msi1 + 396);
    const auto *msi1_397 = buffer.data(msi1 + 397);
    const auto *msi1_398 = buffer.data(msi1 + 398);
    const auto *msi1_399 = buffer.data(msi1 + 399);
    const auto *msi1_400 = buffer.data(msi1 + 400);
    const auto *msi1_401 = buffer.data(msi1 + 401);
    const auto *msi1_402 = buffer.data(msi1 + 402);
    const auto *msi1_403 = buffer.data(msi1 + 403);
    const auto *msi1_404 = buffer.data(msi1 + 404);
    const auto *msi1_405 = buffer.data(msi1 + 405);
    const auto *msi1_406 = buffer.data(msi1 + 406);
    const auto *msi1_412 = buffer.data(msi1 + 412);
    const auto *msi1_413 = buffer.data(msi1 + 413);
    const auto *msi1_414 = buffer.data(msi1 + 414);
    const auto *msi1_415 = buffer.data(msi1 + 415);
    const auto *msi1_416 = buffer.data(msi1 + 416);
    const auto *msi1_417 = buffer.data(msi1 + 417);
    const auto *msi1_418 = buffer.data(msi1 + 418);
    const auto *msi1_419 = buffer.data(msi1 + 419);
    const auto *msi1_420 = buffer.data(msi1 + 420);
    const auto *msi1_422 = buffer.data(msi1 + 422);
    const auto *msi1_423 = buffer.data(msi1 + 423);
    const auto *msi1_425 = buffer.data(msi1 + 425);
    const auto *msi1_426 = buffer.data(msi1 + 426);
    const auto *msi1_427 = buffer.data(msi1 + 427);
    const auto *msi1_429 = buffer.data(msi1 + 429);
    const auto *msi1_430 = buffer.data(msi1 + 430);
    const auto *msi1_431 = buffer.data(msi1 + 431);
    const auto *msi1_432 = buffer.data(msi1 + 432);
    const auto *msi1_434 = buffer.data(msi1 + 434);
    const auto *msi1_435 = buffer.data(msi1 + 435);
    const auto *msi1_441 = buffer.data(msi1 + 441);

    const auto *msk_470 = buffer.data(msk + 470);
    const auto *msk_471 = buffer.data(msk + 471);
    const auto *msk_473 = buffer.data(msk + 473);
    const auto *msk_474 = buffer.data(msk + 474);
    const auto *msk_477 = buffer.data(msk + 477);
    const auto *msk_478 = buffer.data(msk + 478);
    const auto *msk_482 = buffer.data(msk + 482);
    const auto *msk_483 = buffer.data(msk + 483);
    const auto *msk_488 = buffer.data(msk + 488);
    const auto *msk_496 = buffer.data(msk + 496);
    const auto *msk_497 = buffer.data(msk + 497);
    const auto *msk_498 = buffer.data(msk + 498);
    const auto *msk_499 = buffer.data(msk + 499);
    const auto *msk_500 = buffer.data(msk + 500);
    const auto *msk_501 = buffer.data(msk + 501);
    const auto *msk_502 = buffer.data(msk + 502);
    const auto *msk_503 = buffer.data(msk + 503);
    const auto *msk_504 = buffer.data(msk + 504);
    const auto *msk_505 = buffer.data(msk + 505);
    const auto *msk_506 = buffer.data(msk + 506);
    const auto *msk_507 = buffer.data(msk + 507);
    const auto *msk_508 = buffer.data(msk + 508);
    const auto *msk_509 = buffer.data(msk + 509);
    const auto *msk_510 = buffer.data(msk + 510);
    const auto *msk_511 = buffer.data(msk + 511);
    const auto *msk_512 = buffer.data(msk + 512);
    const auto *msk_513 = buffer.data(msk + 513);
    const auto *msk_514 = buffer.data(msk + 514);
    const auto *msk_515 = buffer.data(msk + 515);
    const auto *msk_516 = buffer.data(msk + 516);
    const auto *msk_517 = buffer.data(msk + 517);
    const auto *msk_518 = buffer.data(msk + 518);
    const auto *msk_519 = buffer.data(msk + 519);
    const auto *msk_520 = buffer.data(msk + 520);
    const auto *msk_521 = buffer.data(msk + 521);
    const auto *msk_522 = buffer.data(msk + 522);
    const auto *msk_523 = buffer.data(msk + 523);
    const auto *msk_524 = buffer.data(msk + 524);
    const auto *msk_531 = buffer.data(msk + 531);
    const auto *msk_532 = buffer.data(msk + 532);
    const auto *msk_533 = buffer.data(msk + 533);
    const auto *msk_534 = buffer.data(msk + 534);
    const auto *msk_535 = buffer.data(msk + 535);
    const auto *msk_536 = buffer.data(msk + 536);
    const auto *msk_537 = buffer.data(msk + 537);
    const auto *msk_538 = buffer.data(msk + 538);
    const auto *msk_539 = buffer.data(msk + 539);
    const auto *msk_540 = buffer.data(msk + 540);
    const auto *msk_541 = buffer.data(msk + 541);
    const auto *msk_542 = buffer.data(msk + 542);
    const auto *msk_543 = buffer.data(msk + 543);
    const auto *msk_545 = buffer.data(msk + 545);
    const auto *msk_546 = buffer.data(msk + 546);
    const auto *msk_547 = buffer.data(msk + 547);
    const auto *msk_549 = buffer.data(msk + 549);
    const auto *msk_550 = buffer.data(msk + 550);
    const auto *msk_551 = buffer.data(msk + 551);
    const auto *msk_552 = buffer.data(msk + 552);
    const auto *msk_554 = buffer.data(msk + 554);
    const auto *msk_555 = buffer.data(msk + 555);
    const auto *msk_556 = buffer.data(msk + 556);
    const auto *msk_557 = buffer.data(msk + 557);
    const auto *msk_558 = buffer.data(msk + 558);
    const auto *msk_560 = buffer.data(msk + 560);
    const auto *msk_561 = buffer.data(msk + 561);

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_y, pc_y, lsl0_408, lsl0_410, lsl0_411, \
                         lsk_325, lsk_326, lsk_327, lsl1_408, lsl1_410, lsl1_411, \
                         msk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_y[k] * lsl0_408[k]
                   + f_16 * lsk_325[k]
                   - f_14 * pc_y[k] * lsl1_408[k];

        t_589[k] = f_15 * lsk_326[k]
                   + f_3 * pc_y[k] * msk_470[k];

        t_590[k] = pa_y[k] * lsl0_410[k]
                   - f_14 * pc_y[k] * lsl1_410[k];

        t_591[k] = pa_y[k] * lsl0_411[k]
                   + f_17 * lsk_327[k]
                   - f_14 * pc_y[k] * lsl1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_y, pc_y, pc_z, lsl0_414, lsl0_415, \
                         lsk_291, lsk_329, lsk_330, lsl1_414, lsl1_415, msk_471, \
                         msk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * lsk_291[k]
                   + f_3 * pc_z[k] * msk_471[k];

        t_593[k] = f_15 * lsk_329[k]
                   + f_3 * pc_y[k] * msk_473[k];

        t_594[k] = pa_y[k] * lsl0_414[k]
                   - f_14 * pc_y[k] * lsl1_414[k];

        t_595[k] = pa_y[k] * lsl0_415[k]
                   + f_18 * lsk_330[k]
                   - f_14 * pc_y[k] * lsl1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_y, pc_y, pc_z, lsl0_417, lsl0_419, \
                         lsk_294, lsk_332, lsk_333, lsl1_417, lsl1_419, msk_474, \
                         msk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * lsk_294[k]
                   + f_3 * pc_z[k] * msk_474[k];

        t_597[k] = pa_y[k] * lsl0_417[k]
                   + f_16 * lsk_332[k]
                   - f_14 * pc_y[k] * lsl1_417[k];

        t_598[k] = f_15 * lsk_333[k]
                   + f_3 * pc_y[k] * msk_477[k];

        t_599[k] = pa_y[k] * lsl0_419[k]
                   - f_14 * pc_y[k] * lsl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pa_y, pc_y, pc_z, lsl0_420, lsl0_422, lsk_298, \
                         lsk_334, lsk_336, lsl1_420, lsl1_422, \
                         msk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pa_y[k] * lsl0_420[k]
                   + f_19 * lsk_334[k]
                   - f_14 * pc_y[k] * lsl1_420[k];

        t_601[k] = f_17 * lsk_298[k]
                   + f_3 * pc_z[k] * msk_478[k];

        t_602[k] = pa_y[k] * lsl0_422[k]
                   + f_17 * lsk_336[k]
                   - f_14 * pc_y[k] * lsl1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pa_y, pc_y, lsl0_423, lsl0_425, lsl0_426, \
                         lsk_337, lsk_338, lsk_339, lsl1_423, lsl1_425, lsl1_426, \
                         msk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pa_y[k] * lsl0_423[k]
                   + f_16 * lsk_337[k]
                   - f_14 * pc_y[k] * lsl1_423[k];

        t_604[k] = f_15 * lsk_338[k]
                   + f_3 * pc_y[k] * msk_482[k];

        t_605[k] = pa_y[k] * lsl0_425[k]
                   - f_14 * pc_y[k] * lsl1_425[k];

        t_606[k] = pa_y[k] * lsl0_426[k]
                   + f_20 * lsk_339[k]
                   - f_14 * pc_y[k] * lsl1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pa_y, pc_y, pc_z, lsl0_428, lsl0_429, lsk_303, \
                         lsk_341, lsk_342, lsl1_428, lsl1_429, \
                         msk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * lsk_303[k]
                   + f_3 * pc_z[k] * msk_483[k];

        t_608[k] = pa_y[k] * lsl0_428[k]
                   + f_18 * lsk_341[k]
                   - f_14 * pc_y[k] * lsl1_428[k];

        t_609[k] = pa_y[k] * lsl0_429[k]
                   + f_17 * lsk_342[k]
                   - f_14 * pc_y[k] * lsl1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_y, pc_x, pc_y, lsl0_430, lsl0_432, \
                         lsk_343, lsk_344, lsk_496, lsl1_430, lsl1_432, msk_488, \
                         msk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_y[k] * lsl0_430[k]
                   + f_16 * lsk_343[k]
                   - f_14 * pc_y[k] * lsl1_430[k];

        t_611[k] = f_15 * lsk_344[k]
                   + f_3 * pc_y[k] * msk_488[k];

        t_612[k] = pa_y[k] * lsl0_432[k]
                   - f_14 * pc_y[k] * lsl1_432[k];

        t_613[k] = f_19 * lsk_496[k]
                   + f_3 * pc_x[k] * msk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, lsk_497, lsk_498, lsk_499, \
                         lsk_500, lsk_501, msk_497, msk_498, msk_499, msk_500, \
                         msk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_19 * lsk_497[k]
                   + f_3 * pc_x[k] * msk_497[k];

        t_615[k] = f_19 * lsk_498[k]
                   + f_3 * pc_x[k] * msk_498[k];

        t_616[k] = f_19 * lsk_499[k]
                   + f_3 * pc_x[k] * msk_499[k];

        t_617[k] = f_19 * lsk_500[k]
                   + f_3 * pc_x[k] * msk_500[k];

        t_618[k] = f_19 * lsk_501[k]
                   + f_3 * pc_x[k] * msk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, lsk_316, lsk_352, \
                         lsk_502, lsk_503, msi0_385, msi1_385, msk_496, msk_502, \
                         msk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_19 * lsk_502[k]
                   + f_3 * pc_x[k] * msk_502[k];

        t_620[k] = f_19 * lsk_503[k]
                   + f_3 * pc_x[k] * msk_503[k];

        t_621[k] = f_15 * lsk_352[k]
                   + f_1 * msi0_385[k]
                   - f_2 * msi1_385[k]
                   + f_3 * pc_y[k] * msk_496[k];

        t_622[k] = f_17 * lsk_316[k]
                   + f_3 * pc_z[k] * msk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, lsk_354, lsk_355, lsk_356, msi0_387, \
                         msi0_388, msi0_389, msi1_387, msi1_388, msi1_389, msk_498, msk_499, \
                         msk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * lsk_354[k]
                   + f_12 * msi0_387[k]
                   - f_13 * msi1_387[k]
                   + f_3 * pc_y[k] * msk_498[k];

        t_624[k] = f_15 * lsk_355[k]
                   + f_10 * msi0_388[k]
                   - f_11 * msi1_388[k]
                   + f_3 * pc_y[k] * msk_499[k];

        t_625[k] = f_15 * lsk_356[k]
                   + f_8 * msi0_389[k]
                   - f_9 * msi1_389[k]
                   + f_3 * pc_y[k] * msk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, lsk_357, lsk_358, lsk_359, msi0_390, \
                         msi0_391, msi1_390, msi1_391, msk_501, msk_502, \
                         msk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * lsk_357[k]
                   + f_6 * msi0_390[k]
                   - f_7 * msi1_390[k]
                   + f_3 * pc_y[k] * msk_501[k];

        t_627[k] = f_15 * lsk_358[k]
                   + f_4 * msi0_391[k]
                   - f_5 * msi1_391[k]
                   + f_3 * pc_y[k] * msk_502[k];

        t_628[k] = f_15 * lsk_359[k]
                   + f_3 * pc_y[k] * msk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_y, pc_x, pc_y, pc_z, lsl0_449, \
                         lsk_324, lsk_504, lsl1_449, msi0_392, msi1_392, \
                         msk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_y[k] * lsl0_449[k]
                   - f_14 * pc_y[k] * lsl1_449[k];

        t_630[k] = f_19 * lsk_504[k]
                   + f_1 * msi0_392[k]
                   - f_2 * msi1_392[k]
                   + f_3 * pc_x[k] * msk_504[k];

        t_631[k] = f_3 * pc_y[k] * msk_504[k];

        t_632[k] = f_18 * lsk_324[k]
                   + f_3 * pc_z[k] * msk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, lsk_509, msi0_392, msi0_397, \
                         msi1_392, msi1_397, msk_505, msk_506, \
                         msk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_4 * msi0_392[k]
                   - f_5 * msi1_392[k]
                   + f_3 * pc_y[k] * msk_505[k];

        t_634[k] = f_3 * pc_y[k] * msk_506[k];

        t_635[k] = f_19 * lsk_509[k]
                   + f_12 * msi0_397[k]
                   - f_13 * msi1_397[k]
                   + f_3 * pc_x[k] * msk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_y, msi0_393, msi0_394, msi1_393, msi1_394, \
                         msk_507, msk_508, msk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_6 * msi0_393[k]
                   - f_7 * msi1_393[k]
                   + f_3 * pc_y[k] * msk_507[k];

        t_637[k] = f_4 * msi0_394[k]
                   - f_5 * msi1_394[k]
                   + f_3 * pc_y[k] * msk_508[k];

        t_638[k] = f_3 * pc_y[k] * msk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, lsk_513, msi0_395, msi0_396, \
                         msi0_401, msi1_395, msi1_396, msi1_401, msk_510, msk_511, \
                         msk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_19 * lsk_513[k]
                   + f_10 * msi0_401[k]
                   - f_11 * msi1_401[k]
                   + f_3 * pc_x[k] * msk_513[k];

        t_640[k] = f_8 * msi0_395[k]
                   - f_9 * msi1_395[k]
                   + f_3 * pc_y[k] * msk_510[k];

        t_641[k] = f_6 * msi0_396[k]
                   - f_7 * msi1_396[k]
                   + f_3 * pc_y[k] * msk_511[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, lsk_518, msi0_397, msi0_406, \
                         msi1_397, msi1_406, msk_512, msk_513, \
                         msk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_4 * msi0_397[k]
                   - f_5 * msi1_397[k]
                   + f_3 * pc_y[k] * msk_512[k];

        t_643[k] = f_3 * pc_y[k] * msk_513[k];

        t_644[k] = f_19 * lsk_518[k]
                   + f_8 * msi0_406[k]
                   - f_9 * msi1_406[k]
                   + f_3 * pc_x[k] * msk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_y, msi0_398, msi0_399, msi0_400, msi1_398, \
                         msi1_399, msi1_400, msk_514, msk_515, \
                         msk_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_10 * msi0_398[k]
                   - f_11 * msi1_398[k]
                   + f_3 * pc_y[k] * msk_514[k];

        t_646[k] = f_8 * msi0_399[k]
                   - f_9 * msi1_399[k]
                   + f_3 * pc_y[k] * msk_515[k];

        t_647[k] = f_6 * msi0_400[k]
                   - f_7 * msi1_400[k]
                   + f_3 * pc_y[k] * msk_516[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, lsk_524, msi0_401, msi0_412, \
                         msi1_401, msi1_412, msk_517, msk_518, \
                         msk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_4 * msi0_401[k]
                   - f_5 * msi1_401[k]
                   + f_3 * pc_y[k] * msk_517[k];

        t_649[k] = f_3 * pc_y[k] * msk_518[k];

        t_650[k] = f_19 * lsk_524[k]
                   + f_6 * msi0_412[k]
                   - f_7 * msi1_412[k]
                   + f_3 * pc_x[k] * msk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_y, msi0_402, msi0_403, msi0_404, msi1_402, \
                         msi1_403, msi1_404, msk_519, msk_520, \
                         msk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * msi0_402[k]
                   - f_13 * msi1_402[k]
                   + f_3 * pc_y[k] * msk_519[k];

        t_652[k] = f_10 * msi0_403[k]
                   - f_11 * msi1_403[k]
                   + f_3 * pc_y[k] * msk_520[k];

        t_653[k] = f_8 * msi0_404[k]
                   - f_9 * msi1_404[k]
                   + f_3 * pc_y[k] * msk_521[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_y, msi0_405, msi0_406, msi1_405, msi1_406, \
                         msk_522, msk_523, msk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_6 * msi0_405[k]
                   - f_7 * msi1_405[k]
                   + f_3 * pc_y[k] * msk_522[k];

        t_655[k] = f_4 * msi0_406[k]
                   - f_5 * msi1_406[k]
                   + f_3 * pc_y[k] * msk_523[k];

        t_656[k] = f_3 * pc_y[k] * msk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, lsk_531, lsk_532, lsk_533, lsk_534, \
                         msi0_419, msi1_419, msk_531, msk_532, msk_533, \
                         msk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_19 * lsk_531[k]
                   + f_4 * msi0_419[k]
                   - f_5 * msi1_419[k]
                   + f_3 * pc_x[k] * msk_531[k];

        t_658[k] = f_19 * lsk_532[k]
                   + f_3 * pc_x[k] * msk_532[k];

        t_659[k] = f_19 * lsk_533[k]
                   + f_3 * pc_x[k] * msk_533[k];

        t_660[k] = f_19 * lsk_534[k]
                   + f_3 * pc_x[k] * msk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, pc_y, lsk_535, lsk_536, \
                         lsk_537, lsk_539, msk_531, msk_535, msk_536, msk_537, \
                         msk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_19 * lsk_535[k]
                   + f_3 * pc_x[k] * msk_535[k];

        t_662[k] = f_19 * lsk_536[k]
                   + f_3 * pc_x[k] * msk_536[k];

        t_663[k] = f_19 * lsk_537[k]
                   + f_3 * pc_x[k] * msk_537[k];

        t_664[k] = f_3 * pc_y[k] * msk_531[k];

        t_665[k] = f_19 * lsk_539[k]
                   + f_3 * pc_x[k] * msk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, pc_y, msi0_413, msi0_414, msi0_415, msi1_413, \
                         msi1_414, msi1_415, msk_532, msk_533, \
                         msk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * msi0_413[k]
                   - f_2 * msi1_413[k]
                   + f_3 * pc_y[k] * msk_532[k];

        t_667[k] = f_22 * msi0_414[k]
                   - f_23 * msi1_414[k]
                   + f_3 * pc_y[k] * msk_533[k];

        t_668[k] = f_12 * msi0_415[k]
                   - f_13 * msi1_415[k]
                   + f_3 * pc_y[k] * msk_534[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_y, msi0_416, msi0_417, msi0_418, msi1_416, \
                         msi1_417, msi1_418, msk_535, msk_536, \
                         msk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_10 * msi0_416[k]
                   - f_11 * msi1_416[k]
                   + f_3 * pc_y[k] * msk_535[k];

        t_670[k] = f_8 * msi0_417[k]
                   - f_9 * msi1_417[k]
                   + f_3 * pc_y[k] * msk_536[k];

        t_671[k] = f_6 * msi0_418[k]
                   - f_7 * msi1_418[k]
                   + f_3 * pc_y[k] * msk_537[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, pc_y, pc_z, lsk_359, lsk_540, \
                         msi0_419, msi0_420, msi1_419, msi1_420, msk_538, msk_539, \
                         msk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * msi0_419[k]
                   - f_5 * msi1_419[k]
                   + f_3 * pc_y[k] * msk_538[k];

        t_673[k] = f_3 * pc_y[k] * msk_539[k];

        t_674[k] = f_18 * lsk_359[k]
                   + f_1 * msi0_419[k]
                   - f_2 * msi1_419[k]
                   + f_3 * pc_z[k] * msk_539[k];

        t_675[k] = f_18 * lsk_540[k]
                   + f_1 * msi0_420[k]
                   - f_2 * msi1_420[k]
                   + f_3 * pc_x[k] * msk_540[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pc_x, pc_y, pc_z, lsk_360, lsk_543, \
                         msi0_423, msi1_423, msk_540, msk_541, \
                         msk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_19 * lsk_360[k]
                   + f_3 * pc_y[k] * msk_540[k];

        t_677[k] = f_3 * pc_z[k] * msk_540[k];

        t_678[k] = f_18 * lsk_543[k]
                   + f_12 * msi0_423[k]
                   - f_13 * msi1_423[k]
                   + f_3 * pc_x[k] * msk_543[k];

        t_679[k] = f_3 * pc_z[k] * msk_541[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_z, lsk_546, msi0_420, msi0_426, \
                         msi1_420, msi1_426, msk_542, msk_543, \
                         msk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * msi0_420[k]
                   - f_5 * msi1_420[k]
                   + f_3 * pc_z[k] * msk_542[k];

        t_681[k] = f_18 * lsk_546[k]
                   + f_10 * msi0_426[k]
                   - f_11 * msi1_426[k]
                   + f_3 * pc_x[k] * msk_546[k];

        t_682[k] = f_3 * pc_z[k] * msk_543[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, pc_x, pc_y, pc_z, lsk_365, lsk_550, \
                         msi0_422, msi0_430, msi1_422, msi1_430, msk_545, msk_546, \
                         msk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_19 * lsk_365[k]
                   + f_3 * pc_y[k] * msk_545[k];

        t_684[k] = f_6 * msi0_422[k]
                   - f_7 * msi1_422[k]
                   + f_3 * pc_z[k] * msk_545[k];

        t_685[k] = f_18 * lsk_550[k]
                   + f_8 * msi0_430[k]
                   - f_9 * msi1_430[k]
                   + f_3 * pc_x[k] * msk_550[k];

        t_686[k] = f_3 * pc_z[k] * msk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_y, pc_z, lsk_369, msi0_423, msi0_425, \
                         msi1_423, msi1_425, msk_547, msk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_4 * msi0_423[k]
                   - f_5 * msi1_423[k]
                   + f_3 * pc_z[k] * msk_547[k];

        t_688[k] = f_19 * lsk_369[k]
                   + f_3 * pc_y[k] * msk_549[k];

        t_689[k] = f_8 * msi0_425[k]
                   - f_9 * msi1_425[k]
                   + f_3 * pc_z[k] * msk_549[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_x, pc_z, lsk_555, msi0_426, msi0_435, \
                         msi1_426, msi1_435, msk_550, msk_551, \
                         msk_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_18 * lsk_555[k]
                   + f_6 * msi0_435[k]
                   - f_7 * msi1_435[k]
                   + f_3 * pc_x[k] * msk_555[k];

        t_691[k] = f_3 * pc_z[k] * msk_550[k];

        t_692[k] = f_4 * msi0_426[k]
                   - f_5 * msi1_426[k]
                   + f_3 * pc_z[k] * msk_551[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_y, pc_z, lsk_374, msi0_427, msi0_429, \
                         msi1_427, msi1_429, msk_552, msk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_6 * msi0_427[k]
                   - f_7 * msi1_427[k]
                   + f_3 * pc_z[k] * msk_552[k];

        t_694[k] = f_19 * lsk_374[k]
                   + f_3 * pc_y[k] * msk_554[k];

        t_695[k] = f_10 * msi0_429[k]
                   - f_11 * msi1_429[k]
                   + f_3 * pc_z[k] * msk_554[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_x, pc_z, lsk_561, msi0_430, msi0_441, \
                         msi1_430, msi1_441, msk_555, msk_556, \
                         msk_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_18 * lsk_561[k]
                   + f_4 * msi0_441[k]
                   - f_5 * msi1_441[k]
                   + f_3 * pc_x[k] * msk_561[k];

        t_697[k] = f_3 * pc_z[k] * msk_555[k];

        t_698[k] = f_4 * msi0_430[k]
                   - f_5 * msi1_430[k]
                   + f_3 * pc_z[k] * msk_556[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pc_y, pc_z, lsk_380, msi0_431, msi0_432, \
                         msi0_434, msi1_431, msi1_432, msi1_434, msk_557, msk_558, \
                         msk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_6 * msi0_431[k]
                   - f_7 * msi1_431[k]
                   + f_3 * pc_z[k] * msk_557[k];

        t_700[k] = f_8 * msi0_432[k]
                   - f_9 * msi1_432[k]
                   + f_3 * pc_z[k] * msk_558[k];

        t_701[k] = f_19 * lsk_380[k]
                   + f_3 * pc_y[k] * msk_560[k];

        t_702[k] = f_12 * msi0_434[k]
                   - f_13 * msi1_434[k]
                   + f_3 * pc_z[k] * msk_560[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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

    const auto *lsl0_450 = buffer.data(lsl0 + 450);
    const auto *lsl0_453 = buffer.data(lsl0 + 453);
    const auto *lsl0_456 = buffer.data(lsl0 + 456);
    const auto *lsl0_460 = buffer.data(lsl0 + 460);
    const auto *lsl0_462 = buffer.data(lsl0 + 462);
    const auto *lsl0_465 = buffer.data(lsl0 + 465);
    const auto *lsl0_467 = buffer.data(lsl0 + 467);
    const auto *lsl0_468 = buffer.data(lsl0 + 468);
    const auto *lsl0_471 = buffer.data(lsl0 + 471);
    const auto *lsl0_473 = buffer.data(lsl0 + 473);
    const auto *lsl0_474 = buffer.data(lsl0 + 474);
    const auto *lsl0_475 = buffer.data(lsl0 + 475);
    const auto *lsl0_486 = buffer.data(lsl0 + 486);

    const auto *lsk_360 = buffer.data(lsk + 360);
    const auto *lsk_363 = buffer.data(lsk + 363);
    const auto *lsk_366 = buffer.data(lsk + 366);
    const auto *lsk_367 = buffer.data(lsk + 367);
    const auto *lsk_370 = buffer.data(lsk + 370);
    const auto *lsk_371 = buffer.data(lsk + 371);
    const auto *lsk_372 = buffer.data(lsk + 372);
    const auto *lsk_375 = buffer.data(lsk + 375);
    const auto *lsk_376 = buffer.data(lsk + 376);
    const auto *lsk_377 = buffer.data(lsk + 377);
    const auto *lsk_378 = buffer.data(lsk + 378);
    const auto *lsk_388 = buffer.data(lsk + 388);
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
    const auto *lsk_424 = buffer.data(lsk + 424);
    const auto *lsk_426 = buffer.data(lsk + 426);
    const auto *lsk_427 = buffer.data(lsk + 427);
    const auto *lsk_428 = buffer.data(lsk + 428);
    const auto *lsk_429 = buffer.data(lsk + 429);
    const auto *lsk_430 = buffer.data(lsk + 430);
    const auto *lsk_431 = buffer.data(lsk + 431);
    const auto *lsk_432 = buffer.data(lsk + 432);
    const auto *lsk_434 = buffer.data(lsk + 434);
    const auto *lsk_437 = buffer.data(lsk + 437);
    const auto *lsk_441 = buffer.data(lsk + 441);
    const auto *lsk_446 = buffer.data(lsk + 446);
    const auto *lsk_452 = buffer.data(lsk + 452);
    const auto *lsk_460 = buffer.data(lsk + 460);
    const auto *lsk_462 = buffer.data(lsk + 462);
    const auto *lsk_463 = buffer.data(lsk + 463);
    const auto *lsk_464 = buffer.data(lsk + 464);
    const auto *lsk_465 = buffer.data(lsk + 465);
    const auto *lsk_466 = buffer.data(lsk + 466);
    const auto *lsk_467 = buffer.data(lsk + 467);
    const auto *lsk_468 = buffer.data(lsk + 468);
    const auto *lsk_568 = buffer.data(lsk + 568);
    const auto *lsk_570 = buffer.data(lsk + 570);
    const auto *lsk_571 = buffer.data(lsk + 571);
    const auto *lsk_572 = buffer.data(lsk + 572);
    const auto *lsk_573 = buffer.data(lsk + 573);
    const auto *lsk_574 = buffer.data(lsk + 574);
    const auto *lsk_575 = buffer.data(lsk + 575);
    const auto *lsk_581 = buffer.data(lsk + 581);
    const auto *lsk_585 = buffer.data(lsk + 585);
    const auto *lsk_590 = buffer.data(lsk + 590);
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

    const auto *lsl1_450 = buffer.data(lsl1 + 450);
    const auto *lsl1_453 = buffer.data(lsl1 + 453);
    const auto *lsl1_456 = buffer.data(lsl1 + 456);
    const auto *lsl1_460 = buffer.data(lsl1 + 460);
    const auto *lsl1_462 = buffer.data(lsl1 + 462);
    const auto *lsl1_465 = buffer.data(lsl1 + 465);
    const auto *lsl1_467 = buffer.data(lsl1 + 467);
    const auto *lsl1_468 = buffer.data(lsl1 + 468);
    const auto *lsl1_471 = buffer.data(lsl1 + 471);
    const auto *lsl1_473 = buffer.data(lsl1 + 473);
    const auto *lsl1_474 = buffer.data(lsl1 + 474);
    const auto *lsl1_475 = buffer.data(lsl1 + 475);
    const auto *lsl1_486 = buffer.data(lsl1 + 486);

    const auto *msi0_441 = buffer.data(msi0 + 441);
    const auto *msi0_442 = buffer.data(msi0 + 442);
    const auto *msi0_443 = buffer.data(msi0 + 443);
    const auto *msi0_444 = buffer.data(msi0 + 444);
    const auto *msi0_445 = buffer.data(msi0 + 445);
    const auto *msi0_447 = buffer.data(msi0 + 447);
    const auto *msi0_453 = buffer.data(msi0 + 453);
    const auto *msi0_457 = buffer.data(msi0 + 457);
    const auto *msi0_462 = buffer.data(msi0 + 462);
    const auto *msi0_468 = buffer.data(msi0 + 468);
    const auto *msi0_471 = buffer.data(msi0 + 471);
    const auto *msi0_472 = buffer.data(msi0 + 472);
    const auto *msi0_473 = buffer.data(msi0 + 473);
    const auto *msi0_474 = buffer.data(msi0 + 474);
    const auto *msi0_475 = buffer.data(msi0 + 475);
    const auto *msi0_476 = buffer.data(msi0 + 476);
    const auto *msi0_479 = buffer.data(msi0 + 479);
    const auto *msi0_481 = buffer.data(msi0 + 481);
    const auto *msi0_482 = buffer.data(msi0 + 482);
    const auto *msi0_485 = buffer.data(msi0 + 485);
    const auto *msi0_486 = buffer.data(msi0 + 486);
    const auto *msi0_488 = buffer.data(msi0 + 488);
    const auto *msi0_490 = buffer.data(msi0 + 490);
    const auto *msi0_491 = buffer.data(msi0 + 491);
    const auto *msi0_493 = buffer.data(msi0 + 493);
    const auto *msi0_494 = buffer.data(msi0 + 494);
    const auto *msi0_496 = buffer.data(msi0 + 496);
    const auto *msi0_497 = buffer.data(msi0 + 497);
    const auto *msi0_499 = buffer.data(msi0 + 499);
    const auto *msi0_500 = buffer.data(msi0 + 500);
    const auto *msi0_501 = buffer.data(msi0 + 501);
    const auto *msi0_502 = buffer.data(msi0 + 502);
    const auto *msi0_503 = buffer.data(msi0 + 503);
    const auto *msi0_504 = buffer.data(msi0 + 504);

    const auto *msi1_441 = buffer.data(msi1 + 441);
    const auto *msi1_442 = buffer.data(msi1 + 442);
    const auto *msi1_443 = buffer.data(msi1 + 443);
    const auto *msi1_444 = buffer.data(msi1 + 444);
    const auto *msi1_445 = buffer.data(msi1 + 445);
    const auto *msi1_447 = buffer.data(msi1 + 447);
    const auto *msi1_453 = buffer.data(msi1 + 453);
    const auto *msi1_457 = buffer.data(msi1 + 457);
    const auto *msi1_462 = buffer.data(msi1 + 462);
    const auto *msi1_468 = buffer.data(msi1 + 468);
    const auto *msi1_471 = buffer.data(msi1 + 471);
    const auto *msi1_472 = buffer.data(msi1 + 472);
    const auto *msi1_473 = buffer.data(msi1 + 473);
    const auto *msi1_474 = buffer.data(msi1 + 474);
    const auto *msi1_475 = buffer.data(msi1 + 475);
    const auto *msi1_476 = buffer.data(msi1 + 476);
    const auto *msi1_479 = buffer.data(msi1 + 479);
    const auto *msi1_481 = buffer.data(msi1 + 481);
    const auto *msi1_482 = buffer.data(msi1 + 482);
    const auto *msi1_485 = buffer.data(msi1 + 485);
    const auto *msi1_486 = buffer.data(msi1 + 486);
    const auto *msi1_488 = buffer.data(msi1 + 488);
    const auto *msi1_490 = buffer.data(msi1 + 490);
    const auto *msi1_491 = buffer.data(msi1 + 491);
    const auto *msi1_493 = buffer.data(msi1 + 493);
    const auto *msi1_494 = buffer.data(msi1 + 494);
    const auto *msi1_496 = buffer.data(msi1 + 496);
    const auto *msi1_497 = buffer.data(msi1 + 497);
    const auto *msi1_499 = buffer.data(msi1 + 499);
    const auto *msi1_500 = buffer.data(msi1 + 500);
    const auto *msi1_501 = buffer.data(msi1 + 501);
    const auto *msi1_502 = buffer.data(msi1 + 502);
    const auto *msi1_503 = buffer.data(msi1 + 503);
    const auto *msi1_504 = buffer.data(msi1 + 504);

    const auto *msk_561 = buffer.data(msk + 561);
    const auto *msk_568 = buffer.data(msk + 568);
    const auto *msk_569 = buffer.data(msk + 569);
    const auto *msk_570 = buffer.data(msk + 570);
    const auto *msk_571 = buffer.data(msk + 571);
    const auto *msk_572 = buffer.data(msk + 572);
    const auto *msk_573 = buffer.data(msk + 573);
    const auto *msk_574 = buffer.data(msk + 574);
    const auto *msk_575 = buffer.data(msk + 575);
    const auto *msk_576 = buffer.data(msk + 576);
    const auto *msk_578 = buffer.data(msk + 578);
    const auto *msk_579 = buffer.data(msk + 579);
    const auto *msk_581 = buffer.data(msk + 581);
    const auto *msk_582 = buffer.data(msk + 582);
    const auto *msk_585 = buffer.data(msk + 585);
    const auto *msk_586 = buffer.data(msk + 586);
    const auto *msk_590 = buffer.data(msk + 590);
    const auto *msk_591 = buffer.data(msk + 591);
    const auto *msk_596 = buffer.data(msk + 596);
    const auto *msk_603 = buffer.data(msk + 603);
    const auto *msk_604 = buffer.data(msk + 604);
    const auto *msk_605 = buffer.data(msk + 605);
    const auto *msk_606 = buffer.data(msk + 606);
    const auto *msk_607 = buffer.data(msk + 607);
    const auto *msk_608 = buffer.data(msk + 608);
    const auto *msk_609 = buffer.data(msk + 609);
    const auto *msk_610 = buffer.data(msk + 610);
    const auto *msk_611 = buffer.data(msk + 611);
    const auto *msk_612 = buffer.data(msk + 612);
    const auto *msk_614 = buffer.data(msk + 614);
    const auto *msk_615 = buffer.data(msk + 615);
    const auto *msk_617 = buffer.data(msk + 617);
    const auto *msk_618 = buffer.data(msk + 618);
    const auto *msk_621 = buffer.data(msk + 621);
    const auto *msk_622 = buffer.data(msk + 622);
    const auto *msk_624 = buffer.data(msk + 624);
    const auto *msk_626 = buffer.data(msk + 626);
    const auto *msk_627 = buffer.data(msk + 627);
    const auto *msk_629 = buffer.data(msk + 629);
    const auto *msk_630 = buffer.data(msk + 630);
    const auto *msk_632 = buffer.data(msk + 632);
    const auto *msk_633 = buffer.data(msk + 633);
    const auto *msk_635 = buffer.data(msk + 635);
    const auto *msk_636 = buffer.data(msk + 636);
    const auto *msk_637 = buffer.data(msk + 637);
    const auto *msk_639 = buffer.data(msk + 639);
    const auto *msk_640 = buffer.data(msk + 640);
    const auto *msk_641 = buffer.data(msk + 641);
    const auto *msk_642 = buffer.data(msk + 642);
    const auto *msk_643 = buffer.data(msk + 643);
    const auto *msk_644 = buffer.data(msk + 644);
    const auto *msk_645 = buffer.data(msk + 645);
    const auto *msk_646 = buffer.data(msk + 646);
    const auto *msk_647 = buffer.data(msk + 647);
    const auto *msk_648 = buffer.data(msk + 648);

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, pc_x, pc_z, lsk_568, lsk_570, \
                         lsk_571, lsk_572, msk_561, msk_568, msk_570, msk_571, \
                         msk_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_18 * lsk_568[k]
                   + f_3 * pc_x[k] * msk_568[k];

        t_704[k] = f_3 * pc_z[k] * msk_561[k];

        t_705[k] = f_18 * lsk_570[k]
                   + f_3 * pc_x[k] * msk_570[k];

        t_706[k] = f_18 * lsk_571[k]
                   + f_3 * pc_x[k] * msk_571[k];

        t_707[k] = f_18 * lsk_572[k]
                   + f_3 * pc_x[k] * msk_572[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pc_x, pc_y, lsk_388, lsk_573, lsk_574, \
                         lsk_575, msi0_441, msi1_441, msk_568, msk_573, msk_574, \
                         msk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_18 * lsk_573[k]
                   + f_3 * pc_x[k] * msk_573[k];

        t_709[k] = f_18 * lsk_574[k]
                   + f_3 * pc_x[k] * msk_574[k];

        t_710[k] = f_18 * lsk_575[k]
                   + f_3 * pc_x[k] * msk_575[k];

        t_711[k] = f_19 * lsk_388[k]
                   + f_1 * msi0_441[k]
                   - f_2 * msi1_441[k]
                   + f_3 * pc_y[k] * msk_568[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pc_z, msi0_441, msi0_442, msi0_443, \
                         msi1_441, msi1_442, msi1_443, msk_568, msk_569, msk_570, \
                         msk_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_3 * pc_z[k] * msk_568[k];

        t_713[k] = f_4 * msi0_441[k]
                   - f_5 * msi1_441[k]
                   + f_3 * pc_z[k] * msk_569[k];

        t_714[k] = f_6 * msi0_442[k]
                   - f_7 * msi1_442[k]
                   + f_3 * pc_z[k] * msk_570[k];

        t_715[k] = f_8 * msi0_443[k]
                   - f_9 * msi1_443[k]
                   + f_3 * pc_z[k] * msk_571[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pc_y, pc_z, lsk_395, msi0_444, msi0_445, \
                         msi0_447, msi1_444, msi1_445, msi1_447, msk_572, msk_573, \
                         msk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_10 * msi0_444[k]
                   - f_11 * msi1_444[k]
                   + f_3 * pc_z[k] * msk_572[k];

        t_717[k] = f_12 * msi0_445[k]
                   - f_13 * msi1_445[k]
                   + f_3 * pc_z[k] * msk_573[k];

        t_718[k] = f_19 * lsk_395[k]
                   + f_3 * pc_y[k] * msk_575[k];

        t_719[k] = f_1 * msi0_447[k]
                   - f_2 * msi1_447[k]
                   + f_3 * pc_z[k] * msk_575[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_z, pc_y, pc_z, lsl0_450, lsl0_453, \
                         lsk_360, lsk_396, lsl1_450, lsl1_453, \
                         msk_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_z[k] * lsl0_450[k]
                   - f_14 * pc_z[k] * lsl1_450[k];

        t_721[k] = f_18 * lsk_396[k]
                   + f_3 * pc_y[k] * msk_576[k];

        t_722[k] = f_15 * lsk_360[k]
                   + f_3 * pc_z[k] * msk_576[k];

        t_723[k] = pa_z[k] * lsl0_453[k]
                   - f_14 * pc_z[k] * lsl1_453[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pa_z, pc_x, pc_y, pc_z, lsl0_456, lsk_398, \
                         lsk_581, lsl1_456, msi0_453, msi1_453, msk_578, \
                         msk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_18 * lsk_398[k]
                   + f_3 * pc_y[k] * msk_578[k];

        t_725[k] = f_18 * lsk_581[k]
                   + f_12 * msi0_453[k]
                   - f_13 * msi1_453[k]
                   + f_3 * pc_x[k] * msk_581[k];

        t_726[k] = pa_z[k] * lsl0_456[k]
                   - f_14 * pc_z[k] * lsl1_456[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, lsk_363, lsk_401, lsk_585, \
                         msi0_457, msi1_457, msk_579, msk_581, \
                         msk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_15 * lsk_363[k]
                   + f_3 * pc_z[k] * msk_579[k];

        t_728[k] = f_18 * lsk_401[k]
                   + f_3 * pc_y[k] * msk_581[k];

        t_729[k] = f_18 * lsk_585[k]
                   + f_10 * msi0_457[k]
                   - f_11 * msi1_457[k]
                   + f_3 * pc_x[k] * msk_585[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_z, pc_y, pc_z, lsl0_460, lsl0_462, \
                         lsk_366, lsk_367, lsk_405, lsl1_460, lsl1_462, msk_582, \
                         msk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = pa_z[k] * lsl0_460[k]
                   - f_14 * pc_z[k] * lsl1_460[k];

        t_731[k] = f_15 * lsk_366[k]
                   + f_3 * pc_z[k] * msk_582[k];

        t_732[k] = pa_z[k] * lsl0_462[k]
                   + f_16 * lsk_367[k]
                   - f_14 * pc_z[k] * lsl1_462[k];

        t_733[k] = f_18 * lsk_405[k]
                   + f_3 * pc_y[k] * msk_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pa_z, pc_x, pc_z, lsl0_465, lsk_370, lsk_590, \
                         lsl1_465, msi0_462, msi1_462, msk_586, \
                         msk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_18 * lsk_590[k]
                   + f_8 * msi0_462[k]
                   - f_9 * msi1_462[k]
                   + f_3 * pc_x[k] * msk_590[k];

        t_735[k] = pa_z[k] * lsl0_465[k]
                   - f_14 * pc_z[k] * lsl1_465[k];

        t_736[k] = f_15 * lsk_370[k]
                   + f_3 * pc_z[k] * msk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pa_z, pc_y, pc_z, lsl0_467, lsl0_468, lsk_371, \
                         lsk_372, lsk_410, lsl1_467, lsl1_468, \
                         msk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pa_z[k] * lsl0_467[k]
                   + f_16 * lsk_371[k]
                   - f_14 * pc_z[k] * lsl1_467[k];

        t_738[k] = pa_z[k] * lsl0_468[k]
                   + f_17 * lsk_372[k]
                   - f_14 * pc_z[k] * lsl1_468[k];

        t_739[k] = f_18 * lsk_410[k]
                   + f_3 * pc_y[k] * msk_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pa_z, pc_x, pc_z, lsl0_471, lsk_375, lsk_596, \
                         lsl1_471, msi0_468, msi1_468, msk_591, \
                         msk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_18 * lsk_596[k]
                   + f_6 * msi0_468[k]
                   - f_7 * msi1_468[k]
                   + f_3 * pc_x[k] * msk_596[k];

        t_741[k] = pa_z[k] * lsl0_471[k]
                   - f_14 * pc_z[k] * lsl1_471[k];

        t_742[k] = f_15 * lsk_375[k]
                   + f_3 * pc_z[k] * msk_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pa_z, pc_z, lsl0_473, lsl0_474, lsl0_475, \
                         lsk_376, lsk_377, lsk_378, lsl1_473, lsl1_474, \
                         lsl1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pa_z[k] * lsl0_473[k]
                   + f_16 * lsk_376[k]
                   - f_14 * pc_z[k] * lsl1_473[k];

        t_744[k] = pa_z[k] * lsl0_474[k]
                   + f_17 * lsk_377[k]
                   - f_14 * pc_z[k] * lsl1_474[k];

        t_745[k] = pa_z[k] * lsl0_475[k]
                   + f_18 * lsk_378[k]
                   - f_14 * pc_z[k] * lsl1_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, lsk_416, lsk_603, lsk_604, \
                         lsk_605, msi0_475, msi1_475, msk_596, msk_603, msk_604, \
                         msk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * lsk_416[k]
                   + f_3 * pc_y[k] * msk_596[k];

        t_747[k] = f_18 * lsk_603[k]
                   + f_4 * msi0_475[k]
                   - f_5 * msi1_475[k]
                   + f_3 * pc_x[k] * msk_603[k];

        t_748[k] = f_18 * lsk_604[k]
                   + f_3 * pc_x[k] * msk_604[k];

        t_749[k] = f_18 * lsk_605[k]
                   + f_3 * pc_x[k] * msk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, lsk_606, lsk_607, lsk_608, \
                         lsk_609, lsk_610, msk_606, msk_607, msk_608, msk_609, \
                         msk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_18 * lsk_606[k]
                   + f_3 * pc_x[k] * msk_606[k];

        t_751[k] = f_18 * lsk_607[k]
                   + f_3 * pc_x[k] * msk_607[k];

        t_752[k] = f_18 * lsk_608[k]
                   + f_3 * pc_x[k] * msk_608[k];

        t_753[k] = f_18 * lsk_609[k]
                   + f_3 * pc_x[k] * msk_609[k];

        t_754[k] = f_18 * lsk_610[k]
                   + f_3 * pc_x[k] * msk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pa_z, pc_x, pc_z, lsl0_486, lsk_388, lsk_611, \
                         lsl1_486, msk_604, msk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_18 * lsk_611[k]
                   + f_3 * pc_x[k] * msk_611[k];

        t_756[k] = pa_z[k] * lsl0_486[k]
                   - f_14 * pc_z[k] * lsl1_486[k];

        t_757[k] = f_15 * lsk_388[k]
                   + f_3 * pc_z[k] * msk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, lsk_426, lsk_427, lsk_428, msi0_471, \
                         msi0_472, msi0_473, msi1_471, msi1_472, msi1_473, msk_606, msk_607, \
                         msk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_18 * lsk_426[k]
                   + f_12 * msi0_471[k]
                   - f_13 * msi1_471[k]
                   + f_3 * pc_y[k] * msk_606[k];

        t_759[k] = f_18 * lsk_427[k]
                   + f_10 * msi0_472[k]
                   - f_11 * msi1_472[k]
                   + f_3 * pc_y[k] * msk_607[k];

        t_760[k] = f_18 * lsk_428[k]
                   + f_8 * msi0_473[k]
                   - f_9 * msi1_473[k]
                   + f_3 * pc_y[k] * msk_608[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_y, lsk_429, lsk_430, lsk_431, msi0_474, \
                         msi0_475, msi1_474, msi1_475, msk_609, msk_610, \
                         msk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_18 * lsk_429[k]
                   + f_6 * msi0_474[k]
                   - f_7 * msi1_474[k]
                   + f_3 * pc_y[k] * msk_609[k];

        t_762[k] = f_18 * lsk_430[k]
                   + f_4 * msi0_475[k]
                   - f_5 * msi1_475[k]
                   + f_3 * pc_y[k] * msk_610[k];

        t_763[k] = f_18 * lsk_431[k]
                   + f_3 * pc_y[k] * msk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, pc_z, lsk_395, lsk_432, lsk_612, \
                         msi0_475, msi0_476, msi1_475, msi1_476, msk_611, \
                         msk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * lsk_395[k]
                   + f_1 * msi0_475[k]
                   - f_2 * msi1_475[k]
                   + f_3 * pc_z[k] * msk_611[k];

        t_765[k] = f_18 * lsk_612[k]
                   + f_1 * msi0_476[k]
                   - f_2 * msi1_476[k]
                   + f_3 * pc_x[k] * msk_612[k];

        t_766[k] = f_17 * lsk_432[k]
                   + f_3 * pc_y[k] * msk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_y, pc_z, lsk_396, lsk_434, lsk_615, \
                         msi0_479, msi1_479, msk_612, msk_614, \
                         msk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * lsk_396[k]
                   + f_3 * pc_z[k] * msk_612[k];

        t_768[k] = f_18 * lsk_615[k]
                   + f_12 * msi0_479[k]
                   - f_13 * msi1_479[k]
                   + f_3 * pc_x[k] * msk_615[k];

        t_769[k] = f_17 * lsk_434[k]
                   + f_3 * pc_y[k] * msk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, pc_z, lsk_399, lsk_617, lsk_618, msi0_481, \
                         msi0_482, msi1_481, msi1_482, msk_615, msk_617, \
                         msk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_18 * lsk_617[k]
                   + f_12 * msi0_481[k]
                   - f_13 * msi1_481[k]
                   + f_3 * pc_x[k] * msk_617[k];

        t_771[k] = f_18 * lsk_618[k]
                   + f_10 * msi0_482[k]
                   - f_11 * msi1_482[k]
                   + f_3 * pc_x[k] * msk_618[k];

        t_772[k] = f_16 * lsk_399[k]
                   + f_3 * pc_z[k] * msk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, pc_y, lsk_437, lsk_621, lsk_622, msi0_485, \
                         msi0_486, msi1_485, msi1_486, msk_617, msk_621, \
                         msk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * lsk_437[k]
                   + f_3 * pc_y[k] * msk_617[k];

        t_774[k] = f_18 * lsk_621[k]
                   + f_10 * msi0_485[k]
                   - f_11 * msi1_485[k]
                   + f_3 * pc_x[k] * msk_621[k];

        t_775[k] = f_18 * lsk_622[k]
                   + f_8 * msi0_486[k]
                   - f_9 * msi1_486[k]
                   + f_3 * pc_x[k] * msk_622[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pc_x, pc_y, pc_z, lsk_402, lsk_441, lsk_624, \
                         msi0_488, msi1_488, msk_618, msk_621, \
                         msk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * lsk_402[k]
                   + f_3 * pc_z[k] * msk_618[k];

        t_777[k] = f_18 * lsk_624[k]
                   + f_8 * msi0_488[k]
                   - f_9 * msi1_488[k]
                   + f_3 * pc_x[k] * msk_624[k];

        t_778[k] = f_17 * lsk_441[k]
                   + f_3 * pc_y[k] * msk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_x, pc_z, lsk_406, lsk_626, lsk_627, msi0_490, \
                         msi0_491, msi1_490, msi1_491, msk_622, msk_626, \
                         msk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_18 * lsk_626[k]
                   + f_8 * msi0_490[k]
                   - f_9 * msi1_490[k]
                   + f_3 * pc_x[k] * msk_626[k];

        t_780[k] = f_18 * lsk_627[k]
                   + f_6 * msi0_491[k]
                   - f_7 * msi1_491[k]
                   + f_3 * pc_x[k] * msk_627[k];

        t_781[k] = f_16 * lsk_406[k]
                   + f_3 * pc_z[k] * msk_622[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, lsk_446, lsk_629, lsk_630, msi0_493, \
                         msi0_494, msi1_493, msi1_494, msk_626, msk_629, \
                         msk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_18 * lsk_629[k]
                   + f_6 * msi0_493[k]
                   - f_7 * msi1_493[k]
                   + f_3 * pc_x[k] * msk_629[k];

        t_783[k] = f_18 * lsk_630[k]
                   + f_6 * msi0_494[k]
                   - f_7 * msi1_494[k]
                   + f_3 * pc_x[k] * msk_630[k];

        t_784[k] = f_17 * lsk_446[k]
                   + f_3 * pc_y[k] * msk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_x, pc_z, lsk_411, lsk_632, lsk_633, msi0_496, \
                         msi0_497, msi1_496, msi1_497, msk_627, msk_632, \
                         msk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_18 * lsk_632[k]
                   + f_6 * msi0_496[k]
                   - f_7 * msi1_496[k]
                   + f_3 * pc_x[k] * msk_632[k];

        t_786[k] = f_18 * lsk_633[k]
                   + f_4 * msi0_497[k]
                   - f_5 * msi1_497[k]
                   + f_3 * pc_x[k] * msk_633[k];

        t_787[k] = f_16 * lsk_411[k]
                   + f_3 * pc_z[k] * msk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, lsk_635, lsk_636, lsk_637, msi0_499, \
                         msi0_500, msi0_501, msi1_499, msi1_500, msi1_501, msk_635, msk_636, \
                         msk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_18 * lsk_635[k]
                   + f_4 * msi0_499[k]
                   - f_5 * msi1_499[k]
                   + f_3 * pc_x[k] * msk_635[k];

        t_789[k] = f_18 * lsk_636[k]
                   + f_4 * msi0_500[k]
                   - f_5 * msi1_500[k]
                   + f_3 * pc_x[k] * msk_636[k];

        t_790[k] = f_18 * lsk_637[k]
                   + f_4 * msi0_501[k]
                   - f_5 * msi1_501[k]
                   + f_3 * pc_x[k] * msk_637[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pc_x, pc_y, lsk_452, lsk_639, lsk_640, \
                         lsk_641, msi0_503, msi1_503, msk_632, msk_639, msk_640, \
                         msk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * lsk_452[k]
                   + f_3 * pc_y[k] * msk_632[k];

        t_792[k] = f_18 * lsk_639[k]
                   + f_4 * msi0_503[k]
                   - f_5 * msi1_503[k]
                   + f_3 * pc_x[k] * msk_639[k];

        t_793[k] = f_18 * lsk_640[k]
                   + f_3 * pc_x[k] * msk_640[k];

        t_794[k] = f_18 * lsk_641[k]
                   + f_3 * pc_x[k] * msk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, lsk_642, lsk_643, lsk_644, \
                         lsk_645, lsk_646, msk_642, msk_643, msk_644, msk_645, \
                         msk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_18 * lsk_642[k]
                   + f_3 * pc_x[k] * msk_642[k];

        t_796[k] = f_18 * lsk_643[k]
                   + f_3 * pc_x[k] * msk_643[k];

        t_797[k] = f_18 * lsk_644[k]
                   + f_3 * pc_x[k] * msk_644[k];

        t_798[k] = f_18 * lsk_645[k]
                   + f_3 * pc_x[k] * msk_645[k];

        t_799[k] = f_18 * lsk_646[k]
                   + f_3 * pc_x[k] * msk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, lsk_424, lsk_460, lsk_647, \
                         msi0_497, msi1_497, msk_640, msk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_18 * lsk_647[k]
                   + f_3 * pc_x[k] * msk_647[k];

        t_801[k] = f_17 * lsk_460[k]
                   + f_1 * msi0_497[k]
                   - f_2 * msi1_497[k]
                   + f_3 * pc_y[k] * msk_640[k];

        t_802[k] = f_16 * lsk_424[k]
                   + f_3 * pc_z[k] * msk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, lsk_462, lsk_463, lsk_464, msi0_499, \
                         msi0_500, msi0_501, msi1_499, msi1_500, msi1_501, msk_642, msk_643, \
                         msk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * lsk_462[k]
                   + f_12 * msi0_499[k]
                   - f_13 * msi1_499[k]
                   + f_3 * pc_y[k] * msk_642[k];

        t_804[k] = f_17 * lsk_463[k]
                   + f_10 * msi0_500[k]
                   - f_11 * msi1_500[k]
                   + f_3 * pc_y[k] * msk_643[k];

        t_805[k] = f_17 * lsk_464[k]
                   + f_8 * msi0_501[k]
                   - f_9 * msi1_501[k]
                   + f_3 * pc_y[k] * msk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, lsk_465, lsk_466, lsk_467, msi0_502, \
                         msi0_503, msi1_502, msi1_503, msk_645, msk_646, \
                         msk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * lsk_465[k]
                   + f_6 * msi0_502[k]
                   - f_7 * msi1_502[k]
                   + f_3 * pc_y[k] * msk_645[k];

        t_807[k] = f_17 * lsk_466[k]
                   + f_4 * msi0_503[k]
                   - f_5 * msi1_503[k]
                   + f_3 * pc_y[k] * msk_646[k];

        t_808[k] = f_17 * lsk_467[k]
                   + f_3 * pc_y[k] * msk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_y, pc_z, lsk_431, lsk_468, lsk_648, \
                         msi0_503, msi0_504, msi1_503, msi1_504, msk_647, \
                         msk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * lsk_431[k]
                   + f_1 * msi0_503[k]
                   - f_2 * msi1_503[k]
                   + f_3 * pc_z[k] * msk_647[k];

        t_810[k] = f_18 * lsk_648[k]
                   + f_1 * msi0_504[k]
                   - f_2 * msi1_504[k]
                   + f_3 * pc_x[k] * msk_648[k];

        t_811[k] = f_16 * lsk_468[k]
                   + f_3 * pc_y[k] * msk_648[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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

    const auto *lsl0_630 = buffer.data(lsl0 + 630);
    const auto *lsl0_633 = buffer.data(lsl0 + 633);
    const auto *lsl0_635 = buffer.data(lsl0 + 635);
    const auto *lsl0_636 = buffer.data(lsl0 + 636);
    const auto *lsl0_639 = buffer.data(lsl0 + 639);
    const auto *lsl0_640 = buffer.data(lsl0 + 640);
    const auto *lsl0_642 = buffer.data(lsl0 + 642);
    const auto *lsl0_644 = buffer.data(lsl0 + 644);
    const auto *lsl0_645 = buffer.data(lsl0 + 645);
    const auto *lsl0_647 = buffer.data(lsl0 + 647);
    const auto *lsl0_648 = buffer.data(lsl0 + 648);
    const auto *lsl0_650 = buffer.data(lsl0 + 650);
    const auto *lsl0_651 = buffer.data(lsl0 + 651);
    const auto *lsl0_653 = buffer.data(lsl0 + 653);
    const auto *lsl0_654 = buffer.data(lsl0 + 654);
    const auto *lsl0_655 = buffer.data(lsl0 + 655);
    const auto *lsl0_657 = buffer.data(lsl0 + 657);
    const auto *lsl0_674 = buffer.data(lsl0 + 674);

    const auto *lsk_432 = buffer.data(lsk + 432);
    const auto *lsk_435 = buffer.data(lsk + 435);
    const auto *lsk_438 = buffer.data(lsk + 438);
    const auto *lsk_442 = buffer.data(lsk + 442);
    const auto *lsk_447 = buffer.data(lsk + 447);
    const auto *lsk_460 = buffer.data(lsk + 460);
    const auto *lsk_467 = buffer.data(lsk + 467);
    const auto *lsk_468 = buffer.data(lsk + 468);
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
    const auto *lsk_509 = buffer.data(lsk + 509);
    const auto *lsk_510 = buffer.data(lsk + 510);
    const auto *lsk_512 = buffer.data(lsk + 512);
    const auto *lsk_513 = buffer.data(lsk + 513);
    const auto *lsk_514 = buffer.data(lsk + 514);
    const auto *lsk_516 = buffer.data(lsk + 516);
    const auto *lsk_517 = buffer.data(lsk + 517);
    const auto *lsk_518 = buffer.data(lsk + 518);
    const auto *lsk_519 = buffer.data(lsk + 519);
    const auto *lsk_521 = buffer.data(lsk + 521);
    const auto *lsk_522 = buffer.data(lsk + 522);
    const auto *lsk_523 = buffer.data(lsk + 523);
    const auto *lsk_524 = buffer.data(lsk + 524);
    const auto *lsk_532 = buffer.data(lsk + 532);
    const auto *lsk_534 = buffer.data(lsk + 534);
    const auto *lsk_535 = buffer.data(lsk + 535);
    const auto *lsk_536 = buffer.data(lsk + 536);
    const auto *lsk_537 = buffer.data(lsk + 537);
    const auto *lsk_538 = buffer.data(lsk + 538);
    const auto *lsk_539 = buffer.data(lsk + 539);
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
    const auto *lsk_712 = buffer.data(lsk + 712);
    const auto *lsk_713 = buffer.data(lsk + 713);
    const auto *lsk_714 = buffer.data(lsk + 714);
    const auto *lsk_715 = buffer.data(lsk + 715);
    const auto *lsk_716 = buffer.data(lsk + 716);
    const auto *lsk_717 = buffer.data(lsk + 717);
    const auto *lsk_718 = buffer.data(lsk + 718);
    const auto *lsk_719 = buffer.data(lsk + 719);
    const auto *lsk_720 = buffer.data(lsk + 720);
    const auto *lsk_725 = buffer.data(lsk + 725);
    const auto *lsk_729 = buffer.data(lsk + 729);
    const auto *lsk_734 = buffer.data(lsk + 734);
    const auto *lsk_740 = buffer.data(lsk + 740);

    const auto *lsl1_630 = buffer.data(lsl1 + 630);
    const auto *lsl1_633 = buffer.data(lsl1 + 633);
    const auto *lsl1_635 = buffer.data(lsl1 + 635);
    const auto *lsl1_636 = buffer.data(lsl1 + 636);
    const auto *lsl1_639 = buffer.data(lsl1 + 639);
    const auto *lsl1_640 = buffer.data(lsl1 + 640);
    const auto *lsl1_642 = buffer.data(lsl1 + 642);
    const auto *lsl1_644 = buffer.data(lsl1 + 644);
    const auto *lsl1_645 = buffer.data(lsl1 + 645);
    const auto *lsl1_647 = buffer.data(lsl1 + 647);
    const auto *lsl1_648 = buffer.data(lsl1 + 648);
    const auto *lsl1_650 = buffer.data(lsl1 + 650);
    const auto *lsl1_651 = buffer.data(lsl1 + 651);
    const auto *lsl1_653 = buffer.data(lsl1 + 653);
    const auto *lsl1_654 = buffer.data(lsl1 + 654);
    const auto *lsl1_655 = buffer.data(lsl1 + 655);
    const auto *lsl1_657 = buffer.data(lsl1 + 657);
    const auto *lsl1_674 = buffer.data(lsl1 + 674);

    const auto *msi0_507 = buffer.data(msi0 + 507);
    const auto *msi0_509 = buffer.data(msi0 + 509);
    const auto *msi0_510 = buffer.data(msi0 + 510);
    const auto *msi0_513 = buffer.data(msi0 + 513);
    const auto *msi0_514 = buffer.data(msi0 + 514);
    const auto *msi0_516 = buffer.data(msi0 + 516);
    const auto *msi0_518 = buffer.data(msi0 + 518);
    const auto *msi0_519 = buffer.data(msi0 + 519);
    const auto *msi0_521 = buffer.data(msi0 + 521);
    const auto *msi0_522 = buffer.data(msi0 + 522);
    const auto *msi0_524 = buffer.data(msi0 + 524);
    const auto *msi0_525 = buffer.data(msi0 + 525);
    const auto *msi0_527 = buffer.data(msi0 + 527);
    const auto *msi0_528 = buffer.data(msi0 + 528);
    const auto *msi0_529 = buffer.data(msi0 + 529);
    const auto *msi0_530 = buffer.data(msi0 + 530);
    const auto *msi0_531 = buffer.data(msi0 + 531);
    const auto *msi0_553 = buffer.data(msi0 + 553);
    const auto *msi0_555 = buffer.data(msi0 + 555);
    const auto *msi0_556 = buffer.data(msi0 + 556);
    const auto *msi0_557 = buffer.data(msi0 + 557);
    const auto *msi0_558 = buffer.data(msi0 + 558);
    const auto *msi0_559 = buffer.data(msi0 + 559);
    const auto *msi0_560 = buffer.data(msi0 + 560);
    const auto *msi0_561 = buffer.data(msi0 + 561);
    const auto *msi0_562 = buffer.data(msi0 + 562);
    const auto *msi0_563 = buffer.data(msi0 + 563);
    const auto *msi0_564 = buffer.data(msi0 + 564);
    const auto *msi0_565 = buffer.data(msi0 + 565);
    const auto *msi0_566 = buffer.data(msi0 + 566);
    const auto *msi0_567 = buffer.data(msi0 + 567);
    const auto *msi0_568 = buffer.data(msi0 + 568);
    const auto *msi0_569 = buffer.data(msi0 + 569);
    const auto *msi0_574 = buffer.data(msi0 + 574);
    const auto *msi0_580 = buffer.data(msi0 + 580);

    const auto *msi1_507 = buffer.data(msi1 + 507);
    const auto *msi1_509 = buffer.data(msi1 + 509);
    const auto *msi1_510 = buffer.data(msi1 + 510);
    const auto *msi1_513 = buffer.data(msi1 + 513);
    const auto *msi1_514 = buffer.data(msi1 + 514);
    const auto *msi1_516 = buffer.data(msi1 + 516);
    const auto *msi1_518 = buffer.data(msi1 + 518);
    const auto *msi1_519 = buffer.data(msi1 + 519);
    const auto *msi1_521 = buffer.data(msi1 + 521);
    const auto *msi1_522 = buffer.data(msi1 + 522);
    const auto *msi1_524 = buffer.data(msi1 + 524);
    const auto *msi1_525 = buffer.data(msi1 + 525);
    const auto *msi1_527 = buffer.data(msi1 + 527);
    const auto *msi1_528 = buffer.data(msi1 + 528);
    const auto *msi1_529 = buffer.data(msi1 + 529);
    const auto *msi1_530 = buffer.data(msi1 + 530);
    const auto *msi1_531 = buffer.data(msi1 + 531);
    const auto *msi1_553 = buffer.data(msi1 + 553);
    const auto *msi1_555 = buffer.data(msi1 + 555);
    const auto *msi1_556 = buffer.data(msi1 + 556);
    const auto *msi1_557 = buffer.data(msi1 + 557);
    const auto *msi1_558 = buffer.data(msi1 + 558);
    const auto *msi1_559 = buffer.data(msi1 + 559);
    const auto *msi1_560 = buffer.data(msi1 + 560);
    const auto *msi1_561 = buffer.data(msi1 + 561);
    const auto *msi1_562 = buffer.data(msi1 + 562);
    const auto *msi1_563 = buffer.data(msi1 + 563);
    const auto *msi1_564 = buffer.data(msi1 + 564);
    const auto *msi1_565 = buffer.data(msi1 + 565);
    const auto *msi1_566 = buffer.data(msi1 + 566);
    const auto *msi1_567 = buffer.data(msi1 + 567);
    const auto *msi1_568 = buffer.data(msi1 + 568);
    const auto *msi1_569 = buffer.data(msi1 + 569);
    const auto *msi1_574 = buffer.data(msi1 + 574);
    const auto *msi1_580 = buffer.data(msi1 + 580);

    const auto *msk_648 = buffer.data(msk + 648);
    const auto *msk_650 = buffer.data(msk + 650);
    const auto *msk_651 = buffer.data(msk + 651);
    const auto *msk_653 = buffer.data(msk + 653);
    const auto *msk_654 = buffer.data(msk + 654);
    const auto *msk_657 = buffer.data(msk + 657);
    const auto *msk_658 = buffer.data(msk + 658);
    const auto *msk_660 = buffer.data(msk + 660);
    const auto *msk_662 = buffer.data(msk + 662);
    const auto *msk_663 = buffer.data(msk + 663);
    const auto *msk_665 = buffer.data(msk + 665);
    const auto *msk_666 = buffer.data(msk + 666);
    const auto *msk_668 = buffer.data(msk + 668);
    const auto *msk_669 = buffer.data(msk + 669);
    const auto *msk_671 = buffer.data(msk + 671);
    const auto *msk_672 = buffer.data(msk + 672);
    const auto *msk_673 = buffer.data(msk + 673);
    const auto *msk_675 = buffer.data(msk + 675);
    const auto *msk_676 = buffer.data(msk + 676);
    const auto *msk_677 = buffer.data(msk + 677);
    const auto *msk_678 = buffer.data(msk + 678);
    const auto *msk_679 = buffer.data(msk + 679);
    const auto *msk_680 = buffer.data(msk + 680);
    const auto *msk_681 = buffer.data(msk + 681);
    const auto *msk_682 = buffer.data(msk + 682);
    const auto *msk_683 = buffer.data(msk + 683);
    const auto *msk_684 = buffer.data(msk + 684);
    const auto *msk_686 = buffer.data(msk + 686);
    const auto *msk_687 = buffer.data(msk + 687);
    const auto *msk_689 = buffer.data(msk + 689);
    const auto *msk_690 = buffer.data(msk + 690);
    const auto *msk_693 = buffer.data(msk + 693);
    const auto *msk_694 = buffer.data(msk + 694);
    const auto *msk_698 = buffer.data(msk + 698);
    const auto *msk_699 = buffer.data(msk + 699);
    const auto *msk_704 = buffer.data(msk + 704);
    const auto *msk_712 = buffer.data(msk + 712);
    const auto *msk_713 = buffer.data(msk + 713);
    const auto *msk_714 = buffer.data(msk + 714);
    const auto *msk_715 = buffer.data(msk + 715);
    const auto *msk_716 = buffer.data(msk + 716);
    const auto *msk_717 = buffer.data(msk + 717);
    const auto *msk_718 = buffer.data(msk + 718);
    const auto *msk_719 = buffer.data(msk + 719);
    const auto *msk_720 = buffer.data(msk + 720);
    const auto *msk_721 = buffer.data(msk + 721);
    const auto *msk_722 = buffer.data(msk + 722);
    const auto *msk_723 = buffer.data(msk + 723);
    const auto *msk_724 = buffer.data(msk + 724);
    const auto *msk_725 = buffer.data(msk + 725);
    const auto *msk_726 = buffer.data(msk + 726);
    const auto *msk_727 = buffer.data(msk + 727);
    const auto *msk_728 = buffer.data(msk + 728);
    const auto *msk_729 = buffer.data(msk + 729);
    const auto *msk_730 = buffer.data(msk + 730);
    const auto *msk_731 = buffer.data(msk + 731);
    const auto *msk_732 = buffer.data(msk + 732);
    const auto *msk_733 = buffer.data(msk + 733);
    const auto *msk_734 = buffer.data(msk + 734);
    const auto *msk_740 = buffer.data(msk + 740);

#pragma omp simd aligned(t_812, t_813, t_814, pc_x, pc_y, pc_z, lsk_432, lsk_470, lsk_651, \
                         msi0_507, msi1_507, msk_648, msk_650, \
                         msk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * lsk_432[k]
                   + f_3 * pc_z[k] * msk_648[k];

        t_813[k] = f_18 * lsk_651[k]
                   + f_12 * msi0_507[k]
                   - f_13 * msi1_507[k]
                   + f_3 * pc_x[k] * msk_651[k];

        t_814[k] = f_16 * lsk_470[k]
                   + f_3 * pc_y[k] * msk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_x, pc_z, lsk_435, lsk_653, lsk_654, msi0_509, \
                         msi0_510, msi1_509, msi1_510, msk_651, msk_653, \
                         msk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_18 * lsk_653[k]
                   + f_12 * msi0_509[k]
                   - f_13 * msi1_509[k]
                   + f_3 * pc_x[k] * msk_653[k];

        t_816[k] = f_18 * lsk_654[k]
                   + f_10 * msi0_510[k]
                   - f_11 * msi1_510[k]
                   + f_3 * pc_x[k] * msk_654[k];

        t_817[k] = f_17 * lsk_435[k]
                   + f_3 * pc_z[k] * msk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_y, lsk_473, lsk_657, lsk_658, msi0_513, \
                         msi0_514, msi1_513, msi1_514, msk_653, msk_657, \
                         msk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * lsk_473[k]
                   + f_3 * pc_y[k] * msk_653[k];

        t_819[k] = f_18 * lsk_657[k]
                   + f_10 * msi0_513[k]
                   - f_11 * msi1_513[k]
                   + f_3 * pc_x[k] * msk_657[k];

        t_820[k] = f_18 * lsk_658[k]
                   + f_8 * msi0_514[k]
                   - f_9 * msi1_514[k]
                   + f_3 * pc_x[k] * msk_658[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pc_y, pc_z, lsk_438, lsk_477, lsk_660, \
                         msi0_516, msi1_516, msk_654, msk_657, \
                         msk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * lsk_438[k]
                   + f_3 * pc_z[k] * msk_654[k];

        t_822[k] = f_18 * lsk_660[k]
                   + f_8 * msi0_516[k]
                   - f_9 * msi1_516[k]
                   + f_3 * pc_x[k] * msk_660[k];

        t_823[k] = f_16 * lsk_477[k]
                   + f_3 * pc_y[k] * msk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, lsk_442, lsk_662, lsk_663, msi0_518, \
                         msi0_519, msi1_518, msi1_519, msk_658, msk_662, \
                         msk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_18 * lsk_662[k]
                   + f_8 * msi0_518[k]
                   - f_9 * msi1_518[k]
                   + f_3 * pc_x[k] * msk_662[k];

        t_825[k] = f_18 * lsk_663[k]
                   + f_6 * msi0_519[k]
                   - f_7 * msi1_519[k]
                   + f_3 * pc_x[k] * msk_663[k];

        t_826[k] = f_17 * lsk_442[k]
                   + f_3 * pc_z[k] * msk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pc_x, pc_y, lsk_482, lsk_665, lsk_666, msi0_521, \
                         msi0_522, msi1_521, msi1_522, msk_662, msk_665, \
                         msk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_18 * lsk_665[k]
                   + f_6 * msi0_521[k]
                   - f_7 * msi1_521[k]
                   + f_3 * pc_x[k] * msk_665[k];

        t_828[k] = f_18 * lsk_666[k]
                   + f_6 * msi0_522[k]
                   - f_7 * msi1_522[k]
                   + f_3 * pc_x[k] * msk_666[k];

        t_829[k] = f_16 * lsk_482[k]
                   + f_3 * pc_y[k] * msk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pc_x, pc_z, lsk_447, lsk_668, lsk_669, msi0_524, \
                         msi0_525, msi1_524, msi1_525, msk_663, msk_668, \
                         msk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_18 * lsk_668[k]
                   + f_6 * msi0_524[k]
                   - f_7 * msi1_524[k]
                   + f_3 * pc_x[k] * msk_668[k];

        t_831[k] = f_18 * lsk_669[k]
                   + f_4 * msi0_525[k]
                   - f_5 * msi1_525[k]
                   + f_3 * pc_x[k] * msk_669[k];

        t_832[k] = f_17 * lsk_447[k]
                   + f_3 * pc_z[k] * msk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, lsk_671, lsk_672, lsk_673, msi0_527, \
                         msi0_528, msi0_529, msi1_527, msi1_528, msi1_529, msk_671, msk_672, \
                         msk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_18 * lsk_671[k]
                   + f_4 * msi0_527[k]
                   - f_5 * msi1_527[k]
                   + f_3 * pc_x[k] * msk_671[k];

        t_834[k] = f_18 * lsk_672[k]
                   + f_4 * msi0_528[k]
                   - f_5 * msi1_528[k]
                   + f_3 * pc_x[k] * msk_672[k];

        t_835[k] = f_18 * lsk_673[k]
                   + f_4 * msi0_529[k]
                   - f_5 * msi1_529[k]
                   + f_3 * pc_x[k] * msk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_x, pc_y, lsk_488, lsk_675, lsk_676, \
                         lsk_677, msi0_531, msi1_531, msk_668, msk_675, msk_676, \
                         msk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * lsk_488[k]
                   + f_3 * pc_y[k] * msk_668[k];

        t_837[k] = f_18 * lsk_675[k]
                   + f_4 * msi0_531[k]
                   - f_5 * msi1_531[k]
                   + f_3 * pc_x[k] * msk_675[k];

        t_838[k] = f_18 * lsk_676[k]
                   + f_3 * pc_x[k] * msk_676[k];

        t_839[k] = f_18 * lsk_677[k]
                   + f_3 * pc_x[k] * msk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, lsk_678, lsk_679, lsk_680, \
                         lsk_681, lsk_682, msk_678, msk_679, msk_680, msk_681, \
                         msk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_18 * lsk_678[k]
                   + f_3 * pc_x[k] * msk_678[k];

        t_841[k] = f_18 * lsk_679[k]
                   + f_3 * pc_x[k] * msk_679[k];

        t_842[k] = f_18 * lsk_680[k]
                   + f_3 * pc_x[k] * msk_680[k];

        t_843[k] = f_18 * lsk_681[k]
                   + f_3 * pc_x[k] * msk_681[k];

        t_844[k] = f_18 * lsk_682[k]
                   + f_3 * pc_x[k] * msk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, pc_z, lsk_460, lsk_496, lsk_683, \
                         msi0_525, msi1_525, msk_676, msk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_18 * lsk_683[k]
                   + f_3 * pc_x[k] * msk_683[k];

        t_846[k] = f_16 * lsk_496[k]
                   + f_1 * msi0_525[k]
                   - f_2 * msi1_525[k]
                   + f_3 * pc_y[k] * msk_676[k];

        t_847[k] = f_17 * lsk_460[k]
                   + f_3 * pc_z[k] * msk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, lsk_498, lsk_499, lsk_500, msi0_527, \
                         msi0_528, msi0_529, msi1_527, msi1_528, msi1_529, msk_678, msk_679, \
                         msk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * lsk_498[k]
                   + f_12 * msi0_527[k]
                   - f_13 * msi1_527[k]
                   + f_3 * pc_y[k] * msk_678[k];

        t_849[k] = f_16 * lsk_499[k]
                   + f_10 * msi0_528[k]
                   - f_11 * msi1_528[k]
                   + f_3 * pc_y[k] * msk_679[k];

        t_850[k] = f_16 * lsk_500[k]
                   + f_8 * msi0_529[k]
                   - f_9 * msi1_529[k]
                   + f_3 * pc_y[k] * msk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, lsk_501, lsk_502, lsk_503, msi0_530, \
                         msi0_531, msi1_530, msi1_531, msk_681, msk_682, \
                         msk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * lsk_501[k]
                   + f_6 * msi0_530[k]
                   - f_7 * msi1_530[k]
                   + f_3 * pc_y[k] * msk_681[k];

        t_852[k] = f_16 * lsk_502[k]
                   + f_4 * msi0_531[k]
                   - f_5 * msi1_531[k]
                   + f_3 * pc_y[k] * msk_682[k];

        t_853[k] = f_16 * lsk_503[k]
                   + f_3 * pc_y[k] * msk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_y, pc_y, pc_z, lsl0_630, lsk_467, \
                         lsk_468, lsk_504, lsl1_630, msi0_531, msi1_531, msk_683, \
                         msk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * lsk_467[k]
                   + f_1 * msi0_531[k]
                   - f_2 * msi1_531[k]
                   + f_3 * pc_z[k] * msk_683[k];

        t_855[k] = pa_y[k] * lsl0_630[k]
                   - f_14 * pc_y[k] * lsl1_630[k];

        t_856[k] = f_15 * lsk_504[k]
                   + f_3 * pc_y[k] * msk_684[k];

        t_857[k] = f_18 * lsk_468[k]
                   + f_3 * pc_z[k] * msk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pa_y, pc_y, lsl0_633, lsl0_635, lsl0_636, \
                         lsk_505, lsk_506, lsk_507, lsl1_633, lsl1_635, lsl1_636, \
                         msk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pa_y[k] * lsl0_633[k]
                   + f_16 * lsk_505[k]
                   - f_14 * pc_y[k] * lsl1_633[k];

        t_859[k] = f_15 * lsk_506[k]
                   + f_3 * pc_y[k] * msk_686[k];

        t_860[k] = pa_y[k] * lsl0_635[k]
                   - f_14 * pc_y[k] * lsl1_635[k];

        t_861[k] = pa_y[k] * lsl0_636[k]
                   + f_17 * lsk_507[k]
                   - f_14 * pc_y[k] * lsl1_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pa_y, pc_y, pc_z, lsl0_639, lsl0_640, \
                         lsk_471, lsk_509, lsk_510, lsl1_639, lsl1_640, msk_687, \
                         msk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_18 * lsk_471[k]
                   + f_3 * pc_z[k] * msk_687[k];

        t_863[k] = f_15 * lsk_509[k]
                   + f_3 * pc_y[k] * msk_689[k];

        t_864[k] = pa_y[k] * lsl0_639[k]
                   - f_14 * pc_y[k] * lsl1_639[k];

        t_865[k] = pa_y[k] * lsl0_640[k]
                   + f_18 * lsk_510[k]
                   - f_14 * pc_y[k] * lsl1_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pa_y, pc_y, pc_z, lsl0_642, lsl0_644, \
                         lsk_474, lsk_512, lsk_513, lsl1_642, lsl1_644, msk_690, \
                         msk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * lsk_474[k]
                   + f_3 * pc_z[k] * msk_690[k];

        t_867[k] = pa_y[k] * lsl0_642[k]
                   + f_16 * lsk_512[k]
                   - f_14 * pc_y[k] * lsl1_642[k];

        t_868[k] = f_15 * lsk_513[k]
                   + f_3 * pc_y[k] * msk_693[k];

        t_869[k] = pa_y[k] * lsl0_644[k]
                   - f_14 * pc_y[k] * lsl1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pa_y, pc_y, pc_z, lsl0_645, lsl0_647, lsk_478, \
                         lsk_514, lsk_516, lsl1_645, lsl1_647, \
                         msk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_y[k] * lsl0_645[k]
                   + f_19 * lsk_514[k]
                   - f_14 * pc_y[k] * lsl1_645[k];

        t_871[k] = f_18 * lsk_478[k]
                   + f_3 * pc_z[k] * msk_694[k];

        t_872[k] = pa_y[k] * lsl0_647[k]
                   + f_17 * lsk_516[k]
                   - f_14 * pc_y[k] * lsl1_647[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pa_y, pc_y, lsl0_648, lsl0_650, lsl0_651, \
                         lsk_517, lsk_518, lsk_519, lsl1_648, lsl1_650, lsl1_651, \
                         msk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pa_y[k] * lsl0_648[k]
                   + f_16 * lsk_517[k]
                   - f_14 * pc_y[k] * lsl1_648[k];

        t_874[k] = f_15 * lsk_518[k]
                   + f_3 * pc_y[k] * msk_698[k];

        t_875[k] = pa_y[k] * lsl0_650[k]
                   - f_14 * pc_y[k] * lsl1_650[k];

        t_876[k] = pa_y[k] * lsl0_651[k]
                   + f_20 * lsk_519[k]
                   - f_14 * pc_y[k] * lsl1_651[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pa_y, pc_y, pc_z, lsl0_653, lsl0_654, lsk_483, \
                         lsk_521, lsk_522, lsl1_653, lsl1_654, \
                         msk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * lsk_483[k]
                   + f_3 * pc_z[k] * msk_699[k];

        t_878[k] = pa_y[k] * lsl0_653[k]
                   + f_18 * lsk_521[k]
                   - f_14 * pc_y[k] * lsl1_653[k];

        t_879[k] = pa_y[k] * lsl0_654[k]
                   + f_17 * lsk_522[k]
                   - f_14 * pc_y[k] * lsl1_654[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pa_y, pc_x, pc_y, lsl0_655, lsl0_657, \
                         lsk_523, lsk_524, lsk_712, lsl1_655, lsl1_657, msk_704, \
                         msk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = pa_y[k] * lsl0_655[k]
                   + f_16 * lsk_523[k]
                   - f_14 * pc_y[k] * lsl1_655[k];

        t_881[k] = f_15 * lsk_524[k]
                   + f_3 * pc_y[k] * msk_704[k];

        t_882[k] = pa_y[k] * lsl0_657[k]
                   - f_14 * pc_y[k] * lsl1_657[k];

        t_883[k] = f_18 * lsk_712[k]
                   + f_3 * pc_x[k] * msk_712[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, lsk_713, lsk_714, lsk_715, \
                         lsk_716, lsk_717, msk_713, msk_714, msk_715, msk_716, \
                         msk_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_18 * lsk_713[k]
                   + f_3 * pc_x[k] * msk_713[k];

        t_885[k] = f_18 * lsk_714[k]
                   + f_3 * pc_x[k] * msk_714[k];

        t_886[k] = f_18 * lsk_715[k]
                   + f_3 * pc_x[k] * msk_715[k];

        t_887[k] = f_18 * lsk_716[k]
                   + f_3 * pc_x[k] * msk_716[k];

        t_888[k] = f_18 * lsk_717[k]
                   + f_3 * pc_x[k] * msk_717[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, pc_x, pc_y, pc_z, lsk_496, lsk_532, \
                         lsk_718, lsk_719, msi0_553, msi1_553, msk_712, msk_718, \
                         msk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_18 * lsk_718[k]
                   + f_3 * pc_x[k] * msk_718[k];

        t_890[k] = f_18 * lsk_719[k]
                   + f_3 * pc_x[k] * msk_719[k];

        t_891[k] = f_15 * lsk_532[k]
                   + f_1 * msi0_553[k]
                   - f_2 * msi1_553[k]
                   + f_3 * pc_y[k] * msk_712[k];

        t_892[k] = f_18 * lsk_496[k]
                   + f_3 * pc_z[k] * msk_712[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, lsk_534, lsk_535, lsk_536, msi0_555, \
                         msi0_556, msi0_557, msi1_555, msi1_556, msi1_557, msk_714, msk_715, \
                         msk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * lsk_534[k]
                   + f_12 * msi0_555[k]
                   - f_13 * msi1_555[k]
                   + f_3 * pc_y[k] * msk_714[k];

        t_894[k] = f_15 * lsk_535[k]
                   + f_10 * msi0_556[k]
                   - f_11 * msi1_556[k]
                   + f_3 * pc_y[k] * msk_715[k];

        t_895[k] = f_15 * lsk_536[k]
                   + f_8 * msi0_557[k]
                   - f_9 * msi1_557[k]
                   + f_3 * pc_y[k] * msk_716[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, lsk_537, lsk_538, lsk_539, msi0_558, \
                         msi0_559, msi1_558, msi1_559, msk_717, msk_718, \
                         msk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * lsk_537[k]
                   + f_6 * msi0_558[k]
                   - f_7 * msi1_558[k]
                   + f_3 * pc_y[k] * msk_717[k];

        t_897[k] = f_15 * lsk_538[k]
                   + f_4 * msi0_559[k]
                   - f_5 * msi1_559[k]
                   + f_3 * pc_y[k] * msk_718[k];

        t_898[k] = f_15 * lsk_539[k]
                   + f_3 * pc_y[k] * msk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pa_y, pc_x, pc_y, pc_z, lsl0_674, \
                         lsk_504, lsk_720, lsl1_674, msi0_560, msi1_560, \
                         msk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pa_y[k] * lsl0_674[k]
                   - f_14 * pc_y[k] * lsl1_674[k];

        t_900[k] = f_18 * lsk_720[k]
                   + f_1 * msi0_560[k]
                   - f_2 * msi1_560[k]
                   + f_3 * pc_x[k] * msk_720[k];

        t_901[k] = f_3 * pc_y[k] * msk_720[k];

        t_902[k] = f_19 * lsk_504[k]
                   + f_3 * pc_z[k] * msk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, lsk_725, msi0_560, msi0_565, \
                         msi1_560, msi1_565, msk_721, msk_722, \
                         msk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_4 * msi0_560[k]
                   - f_5 * msi1_560[k]
                   + f_3 * pc_y[k] * msk_721[k];

        t_904[k] = f_3 * pc_y[k] * msk_722[k];

        t_905[k] = f_18 * lsk_725[k]
                   + f_12 * msi0_565[k]
                   - f_13 * msi1_565[k]
                   + f_3 * pc_x[k] * msk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_y, msi0_561, msi0_562, msi1_561, msi1_562, \
                         msk_723, msk_724, msk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_6 * msi0_561[k]
                   - f_7 * msi1_561[k]
                   + f_3 * pc_y[k] * msk_723[k];

        t_907[k] = f_4 * msi0_562[k]
                   - f_5 * msi1_562[k]
                   + f_3 * pc_y[k] * msk_724[k];

        t_908[k] = f_3 * pc_y[k] * msk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_y, lsk_729, msi0_563, msi0_564, \
                         msi0_569, msi1_563, msi1_564, msi1_569, msk_726, msk_727, \
                         msk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_18 * lsk_729[k]
                   + f_10 * msi0_569[k]
                   - f_11 * msi1_569[k]
                   + f_3 * pc_x[k] * msk_729[k];

        t_910[k] = f_8 * msi0_563[k]
                   - f_9 * msi1_563[k]
                   + f_3 * pc_y[k] * msk_726[k];

        t_911[k] = f_6 * msi0_564[k]
                   - f_7 * msi1_564[k]
                   + f_3 * pc_y[k] * msk_727[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, lsk_734, msi0_565, msi0_574, \
                         msi1_565, msi1_574, msk_728, msk_729, \
                         msk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_4 * msi0_565[k]
                   - f_5 * msi1_565[k]
                   + f_3 * pc_y[k] * msk_728[k];

        t_913[k] = f_3 * pc_y[k] * msk_729[k];

        t_914[k] = f_18 * lsk_734[k]
                   + f_8 * msi0_574[k]
                   - f_9 * msi1_574[k]
                   + f_3 * pc_x[k] * msk_734[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_y, msi0_566, msi0_567, msi0_568, msi1_566, \
                         msi1_567, msi1_568, msk_730, msk_731, \
                         msk_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_10 * msi0_566[k]
                   - f_11 * msi1_566[k]
                   + f_3 * pc_y[k] * msk_730[k];

        t_916[k] = f_8 * msi0_567[k]
                   - f_9 * msi1_567[k]
                   + f_3 * pc_y[k] * msk_731[k];

        t_917[k] = f_6 * msi0_568[k]
                   - f_7 * msi1_568[k]
                   + f_3 * pc_y[k] * msk_732[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_x, pc_y, lsk_740, msi0_569, msi0_580, \
                         msi1_569, msi1_580, msk_733, msk_734, \
                         msk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_4 * msi0_569[k]
                   - f_5 * msi1_569[k]
                   + f_3 * pc_y[k] * msk_733[k];

        t_919[k] = f_3 * pc_y[k] * msk_734[k];

        t_920[k] = f_18 * lsk_740[k]
                   + f_6 * msi0_580[k]
                   - f_7 * msi1_580[k]
                   + f_3 * pc_x[k] * msk_740[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsl0,
                                                          const size_t lsk, const size_t lsl1,
                                                          const size_t msi0, const size_t msi1,
                                                          const size_t msk, const size_t ncols,
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

    const auto *lsl0_675 = buffer.data(lsl0 + 675);
    const auto *lsl0_678 = buffer.data(lsl0 + 678);
    const auto *lsl0_681 = buffer.data(lsl0 + 681);
    const auto *lsl0_685 = buffer.data(lsl0 + 685);
    const auto *lsl0_687 = buffer.data(lsl0 + 687);
    const auto *lsl0_690 = buffer.data(lsl0 + 690);
    const auto *lsl0_692 = buffer.data(lsl0 + 692);
    const auto *lsl0_693 = buffer.data(lsl0 + 693);
    const auto *lsl0_696 = buffer.data(lsl0 + 696);
    const auto *lsl0_698 = buffer.data(lsl0 + 698);
    const auto *lsl0_699 = buffer.data(lsl0 + 699);
    const auto *lsl0_700 = buffer.data(lsl0 + 700);
    const auto *lsl0_711 = buffer.data(lsl0 + 711);

    const auto *lsk_539 = buffer.data(lsk + 539);
    const auto *lsk_540 = buffer.data(lsk + 540);
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
    const auto *lsk_568 = buffer.data(lsk + 568);
    const auto *lsk_575 = buffer.data(lsk + 575);
    const auto *lsk_576 = buffer.data(lsk + 576);
    const auto *lsk_578 = buffer.data(lsk + 578);
    const auto *lsk_581 = buffer.data(lsk + 581);
    const auto *lsk_585 = buffer.data(lsk + 585);
    const auto *lsk_590 = buffer.data(lsk + 590);
    const auto *lsk_596 = buffer.data(lsk + 596);
    const auto *lsk_606 = buffer.data(lsk + 606);
    const auto *lsk_607 = buffer.data(lsk + 607);
    const auto *lsk_608 = buffer.data(lsk + 608);
    const auto *lsk_609 = buffer.data(lsk + 609);
    const auto *lsk_610 = buffer.data(lsk + 610);
    const auto *lsk_611 = buffer.data(lsk + 611);
    const auto *lsk_612 = buffer.data(lsk + 612);
    const auto *lsk_747 = buffer.data(lsk + 747);
    const auto *lsk_748 = buffer.data(lsk + 748);
    const auto *lsk_749 = buffer.data(lsk + 749);
    const auto *lsk_750 = buffer.data(lsk + 750);
    const auto *lsk_751 = buffer.data(lsk + 751);
    const auto *lsk_752 = buffer.data(lsk + 752);
    const auto *lsk_753 = buffer.data(lsk + 753);
    const auto *lsk_755 = buffer.data(lsk + 755);
    const auto *lsk_756 = buffer.data(lsk + 756);
    const auto *lsk_759 = buffer.data(lsk + 759);
    const auto *lsk_762 = buffer.data(lsk + 762);
    const auto *lsk_766 = buffer.data(lsk + 766);
    const auto *lsk_771 = buffer.data(lsk + 771);
    const auto *lsk_777 = buffer.data(lsk + 777);
    const auto *lsk_784 = buffer.data(lsk + 784);
    const auto *lsk_786 = buffer.data(lsk + 786);
    const auto *lsk_787 = buffer.data(lsk + 787);
    const auto *lsk_788 = buffer.data(lsk + 788);
    const auto *lsk_789 = buffer.data(lsk + 789);
    const auto *lsk_790 = buffer.data(lsk + 790);
    const auto *lsk_791 = buffer.data(lsk + 791);
    const auto *lsk_797 = buffer.data(lsk + 797);
    const auto *lsk_801 = buffer.data(lsk + 801);
    const auto *lsk_806 = buffer.data(lsk + 806);
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

    const auto *lsl1_675 = buffer.data(lsl1 + 675);
    const auto *lsl1_678 = buffer.data(lsl1 + 678);
    const auto *lsl1_681 = buffer.data(lsl1 + 681);
    const auto *lsl1_685 = buffer.data(lsl1 + 685);
    const auto *lsl1_687 = buffer.data(lsl1 + 687);
    const auto *lsl1_690 = buffer.data(lsl1 + 690);
    const auto *lsl1_692 = buffer.data(lsl1 + 692);
    const auto *lsl1_693 = buffer.data(lsl1 + 693);
    const auto *lsl1_696 = buffer.data(lsl1 + 696);
    const auto *lsl1_698 = buffer.data(lsl1 + 698);
    const auto *lsl1_699 = buffer.data(lsl1 + 699);
    const auto *lsl1_700 = buffer.data(lsl1 + 700);
    const auto *lsl1_711 = buffer.data(lsl1 + 711);

    const auto *msi0_570 = buffer.data(msi0 + 570);
    const auto *msi0_571 = buffer.data(msi0 + 571);
    const auto *msi0_572 = buffer.data(msi0 + 572);
    const auto *msi0_573 = buffer.data(msi0 + 573);
    const auto *msi0_574 = buffer.data(msi0 + 574);
    const auto *msi0_581 = buffer.data(msi0 + 581);
    const auto *msi0_582 = buffer.data(msi0 + 582);
    const auto *msi0_583 = buffer.data(msi0 + 583);
    const auto *msi0_584 = buffer.data(msi0 + 584);
    const auto *msi0_585 = buffer.data(msi0 + 585);
    const auto *msi0_586 = buffer.data(msi0 + 586);
    const auto *msi0_587 = buffer.data(msi0 + 587);
    const auto *msi0_588 = buffer.data(msi0 + 588);
    const auto *msi0_590 = buffer.data(msi0 + 590);
    const auto *msi0_591 = buffer.data(msi0 + 591);
    const auto *msi0_593 = buffer.data(msi0 + 593);
    const auto *msi0_594 = buffer.data(msi0 + 594);
    const auto *msi0_595 = buffer.data(msi0 + 595);
    const auto *msi0_597 = buffer.data(msi0 + 597);
    const auto *msi0_598 = buffer.data(msi0 + 598);
    const auto *msi0_599 = buffer.data(msi0 + 599);
    const auto *msi0_600 = buffer.data(msi0 + 600);
    const auto *msi0_602 = buffer.data(msi0 + 602);
    const auto *msi0_603 = buffer.data(msi0 + 603);
    const auto *msi0_609 = buffer.data(msi0 + 609);
    const auto *msi0_610 = buffer.data(msi0 + 610);
    const auto *msi0_611 = buffer.data(msi0 + 611);
    const auto *msi0_612 = buffer.data(msi0 + 612);
    const auto *msi0_613 = buffer.data(msi0 + 613);
    const auto *msi0_615 = buffer.data(msi0 + 615);
    const auto *msi0_621 = buffer.data(msi0 + 621);
    const auto *msi0_625 = buffer.data(msi0 + 625);
    const auto *msi0_630 = buffer.data(msi0 + 630);
    const auto *msi0_636 = buffer.data(msi0 + 636);
    const auto *msi0_639 = buffer.data(msi0 + 639);
    const auto *msi0_640 = buffer.data(msi0 + 640);
    const auto *msi0_641 = buffer.data(msi0 + 641);
    const auto *msi0_642 = buffer.data(msi0 + 642);
    const auto *msi0_643 = buffer.data(msi0 + 643);
    const auto *msi0_644 = buffer.data(msi0 + 644);

    const auto *msi1_570 = buffer.data(msi1 + 570);
    const auto *msi1_571 = buffer.data(msi1 + 571);
    const auto *msi1_572 = buffer.data(msi1 + 572);
    const auto *msi1_573 = buffer.data(msi1 + 573);
    const auto *msi1_574 = buffer.data(msi1 + 574);
    const auto *msi1_581 = buffer.data(msi1 + 581);
    const auto *msi1_582 = buffer.data(msi1 + 582);
    const auto *msi1_583 = buffer.data(msi1 + 583);
    const auto *msi1_584 = buffer.data(msi1 + 584);
    const auto *msi1_585 = buffer.data(msi1 + 585);
    const auto *msi1_586 = buffer.data(msi1 + 586);
    const auto *msi1_587 = buffer.data(msi1 + 587);
    const auto *msi1_588 = buffer.data(msi1 + 588);
    const auto *msi1_590 = buffer.data(msi1 + 590);
    const auto *msi1_591 = buffer.data(msi1 + 591);
    const auto *msi1_593 = buffer.data(msi1 + 593);
    const auto *msi1_594 = buffer.data(msi1 + 594);
    const auto *msi1_595 = buffer.data(msi1 + 595);
    const auto *msi1_597 = buffer.data(msi1 + 597);
    const auto *msi1_598 = buffer.data(msi1 + 598);
    const auto *msi1_599 = buffer.data(msi1 + 599);
    const auto *msi1_600 = buffer.data(msi1 + 600);
    const auto *msi1_602 = buffer.data(msi1 + 602);
    const auto *msi1_603 = buffer.data(msi1 + 603);
    const auto *msi1_609 = buffer.data(msi1 + 609);
    const auto *msi1_610 = buffer.data(msi1 + 610);
    const auto *msi1_611 = buffer.data(msi1 + 611);
    const auto *msi1_612 = buffer.data(msi1 + 612);
    const auto *msi1_613 = buffer.data(msi1 + 613);
    const auto *msi1_615 = buffer.data(msi1 + 615);
    const auto *msi1_621 = buffer.data(msi1 + 621);
    const auto *msi1_625 = buffer.data(msi1 + 625);
    const auto *msi1_630 = buffer.data(msi1 + 630);
    const auto *msi1_636 = buffer.data(msi1 + 636);
    const auto *msi1_639 = buffer.data(msi1 + 639);
    const auto *msi1_640 = buffer.data(msi1 + 640);
    const auto *msi1_641 = buffer.data(msi1 + 641);
    const auto *msi1_642 = buffer.data(msi1 + 642);
    const auto *msi1_643 = buffer.data(msi1 + 643);
    const auto *msi1_644 = buffer.data(msi1 + 644);

    const auto *msk_735 = buffer.data(msk + 735);
    const auto *msk_736 = buffer.data(msk + 736);
    const auto *msk_737 = buffer.data(msk + 737);
    const auto *msk_738 = buffer.data(msk + 738);
    const auto *msk_739 = buffer.data(msk + 739);
    const auto *msk_740 = buffer.data(msk + 740);
    const auto *msk_747 = buffer.data(msk + 747);
    const auto *msk_748 = buffer.data(msk + 748);
    const auto *msk_749 = buffer.data(msk + 749);
    const auto *msk_750 = buffer.data(msk + 750);
    const auto *msk_751 = buffer.data(msk + 751);
    const auto *msk_752 = buffer.data(msk + 752);
    const auto *msk_753 = buffer.data(msk + 753);
    const auto *msk_754 = buffer.data(msk + 754);
    const auto *msk_755 = buffer.data(msk + 755);
    const auto *msk_756 = buffer.data(msk + 756);
    const auto *msk_757 = buffer.data(msk + 757);
    const auto *msk_758 = buffer.data(msk + 758);
    const auto *msk_759 = buffer.data(msk + 759);
    const auto *msk_761 = buffer.data(msk + 761);
    const auto *msk_762 = buffer.data(msk + 762);
    const auto *msk_763 = buffer.data(msk + 763);
    const auto *msk_765 = buffer.data(msk + 765);
    const auto *msk_766 = buffer.data(msk + 766);
    const auto *msk_767 = buffer.data(msk + 767);
    const auto *msk_768 = buffer.data(msk + 768);
    const auto *msk_770 = buffer.data(msk + 770);
    const auto *msk_771 = buffer.data(msk + 771);
    const auto *msk_772 = buffer.data(msk + 772);
    const auto *msk_773 = buffer.data(msk + 773);
    const auto *msk_774 = buffer.data(msk + 774);
    const auto *msk_776 = buffer.data(msk + 776);
    const auto *msk_777 = buffer.data(msk + 777);
    const auto *msk_784 = buffer.data(msk + 784);
    const auto *msk_785 = buffer.data(msk + 785);
    const auto *msk_786 = buffer.data(msk + 786);
    const auto *msk_787 = buffer.data(msk + 787);
    const auto *msk_788 = buffer.data(msk + 788);
    const auto *msk_789 = buffer.data(msk + 789);
    const auto *msk_790 = buffer.data(msk + 790);
    const auto *msk_791 = buffer.data(msk + 791);
    const auto *msk_792 = buffer.data(msk + 792);
    const auto *msk_794 = buffer.data(msk + 794);
    const auto *msk_795 = buffer.data(msk + 795);
    const auto *msk_797 = buffer.data(msk + 797);
    const auto *msk_798 = buffer.data(msk + 798);
    const auto *msk_801 = buffer.data(msk + 801);
    const auto *msk_802 = buffer.data(msk + 802);
    const auto *msk_806 = buffer.data(msk + 806);
    const auto *msk_807 = buffer.data(msk + 807);
    const auto *msk_812 = buffer.data(msk + 812);
    const auto *msk_819 = buffer.data(msk + 819);
    const auto *msk_820 = buffer.data(msk + 820);
    const auto *msk_821 = buffer.data(msk + 821);
    const auto *msk_822 = buffer.data(msk + 822);
    const auto *msk_823 = buffer.data(msk + 823);
    const auto *msk_824 = buffer.data(msk + 824);
    const auto *msk_825 = buffer.data(msk + 825);
    const auto *msk_826 = buffer.data(msk + 826);
    const auto *msk_827 = buffer.data(msk + 827);
    const auto *msk_828 = buffer.data(msk + 828);

#pragma omp simd aligned(t_921, t_922, t_923, pc_y, msi0_570, msi0_571, msi0_572, msi1_570, \
                         msi1_571, msi1_572, msk_735, msk_736, \
                         msk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_12 * msi0_570[k]
                   - f_13 * msi1_570[k]
                   + f_3 * pc_y[k] * msk_735[k];

        t_922[k] = f_10 * msi0_571[k]
                   - f_11 * msi1_571[k]
                   + f_3 * pc_y[k] * msk_736[k];

        t_923[k] = f_8 * msi0_572[k]
                   - f_9 * msi1_572[k]
                   + f_3 * pc_y[k] * msk_737[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_y, msi0_573, msi0_574, msi1_573, msi1_574, \
                         msk_738, msk_739, msk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_6 * msi0_573[k]
                   - f_7 * msi1_573[k]
                   + f_3 * pc_y[k] * msk_738[k];

        t_925[k] = f_4 * msi0_574[k]
                   - f_5 * msi1_574[k]
                   + f_3 * pc_y[k] * msk_739[k];

        t_926[k] = f_3 * pc_y[k] * msk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pc_x, lsk_747, lsk_748, lsk_749, lsk_750, \
                         msi0_587, msi1_587, msk_747, msk_748, msk_749, \
                         msk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_18 * lsk_747[k]
                   + f_4 * msi0_587[k]
                   - f_5 * msi1_587[k]
                   + f_3 * pc_x[k] * msk_747[k];

        t_928[k] = f_18 * lsk_748[k]
                   + f_3 * pc_x[k] * msk_748[k];

        t_929[k] = f_18 * lsk_749[k]
                   + f_3 * pc_x[k] * msk_749[k];

        t_930[k] = f_18 * lsk_750[k]
                   + f_3 * pc_x[k] * msk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, pc_y, lsk_751, lsk_752, \
                         lsk_753, lsk_755, msk_747, msk_751, msk_752, msk_753, \
                         msk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_18 * lsk_751[k]
                   + f_3 * pc_x[k] * msk_751[k];

        t_932[k] = f_18 * lsk_752[k]
                   + f_3 * pc_x[k] * msk_752[k];

        t_933[k] = f_18 * lsk_753[k]
                   + f_3 * pc_x[k] * msk_753[k];

        t_934[k] = f_3 * pc_y[k] * msk_747[k];

        t_935[k] = f_18 * lsk_755[k]
                   + f_3 * pc_x[k] * msk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, pc_y, msi0_581, msi0_582, msi0_583, msi1_581, \
                         msi1_582, msi1_583, msk_748, msk_749, \
                         msk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * msi0_581[k]
                   - f_2 * msi1_581[k]
                   + f_3 * pc_y[k] * msk_748[k];

        t_937[k] = f_22 * msi0_582[k]
                   - f_23 * msi1_582[k]
                   + f_3 * pc_y[k] * msk_749[k];

        t_938[k] = f_12 * msi0_583[k]
                   - f_13 * msi1_583[k]
                   + f_3 * pc_y[k] * msk_750[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, pc_y, msi0_584, msi0_585, msi0_586, msi1_584, \
                         msi1_585, msi1_586, msk_751, msk_752, \
                         msk_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_10 * msi0_584[k]
                   - f_11 * msi1_584[k]
                   + f_3 * pc_y[k] * msk_751[k];

        t_940[k] = f_8 * msi0_585[k]
                   - f_9 * msi1_585[k]
                   + f_3 * pc_y[k] * msk_752[k];

        t_941[k] = f_6 * msi0_586[k]
                   - f_7 * msi1_586[k]
                   + f_3 * pc_y[k] * msk_753[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, pc_z, lsk_539, lsk_756, \
                         msi0_587, msi0_588, msi1_587, msi1_588, msk_754, msk_755, \
                         msk_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_4 * msi0_587[k]
                   - f_5 * msi1_587[k]
                   + f_3 * pc_y[k] * msk_754[k];

        t_943[k] = f_3 * pc_y[k] * msk_755[k];

        t_944[k] = f_19 * lsk_539[k]
                   + f_1 * msi0_587[k]
                   - f_2 * msi1_587[k]
                   + f_3 * pc_z[k] * msk_755[k];

        t_945[k] = f_17 * lsk_756[k]
                   + f_1 * msi0_588[k]
                   - f_2 * msi1_588[k]
                   + f_3 * pc_x[k] * msk_756[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pc_x, pc_y, pc_z, lsk_540, lsk_759, \
                         msi0_591, msi1_591, msk_756, msk_757, \
                         msk_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_20 * lsk_540[k]
                   + f_3 * pc_y[k] * msk_756[k];

        t_947[k] = f_3 * pc_z[k] * msk_756[k];

        t_948[k] = f_17 * lsk_759[k]
                   + f_12 * msi0_591[k]
                   - f_13 * msi1_591[k]
                   + f_3 * pc_x[k] * msk_759[k];

        t_949[k] = f_3 * pc_z[k] * msk_757[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pc_x, pc_z, lsk_762, msi0_588, msi0_594, \
                         msi1_588, msi1_594, msk_758, msk_759, \
                         msk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_4 * msi0_588[k]
                   - f_5 * msi1_588[k]
                   + f_3 * pc_z[k] * msk_758[k];

        t_951[k] = f_17 * lsk_762[k]
                   + f_10 * msi0_594[k]
                   - f_11 * msi1_594[k]
                   + f_3 * pc_x[k] * msk_762[k];

        t_952[k] = f_3 * pc_z[k] * msk_759[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pc_x, pc_y, pc_z, lsk_545, lsk_766, \
                         msi0_590, msi0_598, msi1_590, msi1_598, msk_761, msk_762, \
                         msk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_20 * lsk_545[k]
                   + f_3 * pc_y[k] * msk_761[k];

        t_954[k] = f_6 * msi0_590[k]
                   - f_7 * msi1_590[k]
                   + f_3 * pc_z[k] * msk_761[k];

        t_955[k] = f_17 * lsk_766[k]
                   + f_8 * msi0_598[k]
                   - f_9 * msi1_598[k]
                   + f_3 * pc_x[k] * msk_766[k];

        t_956[k] = f_3 * pc_z[k] * msk_762[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pc_y, pc_z, lsk_549, msi0_591, msi0_593, \
                         msi1_591, msi1_593, msk_763, msk_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_4 * msi0_591[k]
                   - f_5 * msi1_591[k]
                   + f_3 * pc_z[k] * msk_763[k];

        t_958[k] = f_20 * lsk_549[k]
                   + f_3 * pc_y[k] * msk_765[k];

        t_959[k] = f_8 * msi0_593[k]
                   - f_9 * msi1_593[k]
                   + f_3 * pc_z[k] * msk_765[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pc_x, pc_z, lsk_771, msi0_594, msi0_603, \
                         msi1_594, msi1_603, msk_766, msk_767, \
                         msk_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_17 * lsk_771[k]
                   + f_6 * msi0_603[k]
                   - f_7 * msi1_603[k]
                   + f_3 * pc_x[k] * msk_771[k];

        t_961[k] = f_3 * pc_z[k] * msk_766[k];

        t_962[k] = f_4 * msi0_594[k]
                   - f_5 * msi1_594[k]
                   + f_3 * pc_z[k] * msk_767[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pc_y, pc_z, lsk_554, msi0_595, msi0_597, \
                         msi1_595, msi1_597, msk_768, msk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_6 * msi0_595[k]
                   - f_7 * msi1_595[k]
                   + f_3 * pc_z[k] * msk_768[k];

        t_964[k] = f_20 * lsk_554[k]
                   + f_3 * pc_y[k] * msk_770[k];

        t_965[k] = f_10 * msi0_597[k]
                   - f_11 * msi1_597[k]
                   + f_3 * pc_z[k] * msk_770[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_x, pc_z, lsk_777, msi0_598, msi0_609, \
                         msi1_598, msi1_609, msk_771, msk_772, \
                         msk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_17 * lsk_777[k]
                   + f_4 * msi0_609[k]
                   - f_5 * msi1_609[k]
                   + f_3 * pc_x[k] * msk_777[k];

        t_967[k] = f_3 * pc_z[k] * msk_771[k];

        t_968[k] = f_4 * msi0_598[k]
                   - f_5 * msi1_598[k]
                   + f_3 * pc_z[k] * msk_772[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, pc_y, pc_z, lsk_560, msi0_599, msi0_600, \
                         msi0_602, msi1_599, msi1_600, msi1_602, msk_773, msk_774, \
                         msk_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_6 * msi0_599[k]
                   - f_7 * msi1_599[k]
                   + f_3 * pc_z[k] * msk_773[k];

        t_970[k] = f_8 * msi0_600[k]
                   - f_9 * msi1_600[k]
                   + f_3 * pc_z[k] * msk_774[k];

        t_971[k] = f_20 * lsk_560[k]
                   + f_3 * pc_y[k] * msk_776[k];

        t_972[k] = f_12 * msi0_602[k]
                   - f_13 * msi1_602[k]
                   + f_3 * pc_z[k] * msk_776[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, t_977, pc_x, pc_z, lsk_784, lsk_786, \
                         lsk_787, lsk_788, msk_777, msk_784, msk_786, msk_787, \
                         msk_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_17 * lsk_784[k]
                   + f_3 * pc_x[k] * msk_784[k];

        t_974[k] = f_3 * pc_z[k] * msk_777[k];

        t_975[k] = f_17 * lsk_786[k]
                   + f_3 * pc_x[k] * msk_786[k];

        t_976[k] = f_17 * lsk_787[k]
                   + f_3 * pc_x[k] * msk_787[k];

        t_977[k] = f_17 * lsk_788[k]
                   + f_3 * pc_x[k] * msk_788[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pc_x, pc_y, lsk_568, lsk_789, lsk_790, \
                         lsk_791, msi0_609, msi1_609, msk_784, msk_789, msk_790, \
                         msk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_17 * lsk_789[k]
                   + f_3 * pc_x[k] * msk_789[k];

        t_979[k] = f_17 * lsk_790[k]
                   + f_3 * pc_x[k] * msk_790[k];

        t_980[k] = f_17 * lsk_791[k]
                   + f_3 * pc_x[k] * msk_791[k];

        t_981[k] = f_20 * lsk_568[k]
                   + f_1 * msi0_609[k]
                   - f_2 * msi1_609[k]
                   + f_3 * pc_y[k] * msk_784[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, t_985, pc_z, msi0_609, msi0_610, msi0_611, \
                         msi1_609, msi1_610, msi1_611, msk_784, msk_785, msk_786, \
                         msk_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_3 * pc_z[k] * msk_784[k];

        t_983[k] = f_4 * msi0_609[k]
                   - f_5 * msi1_609[k]
                   + f_3 * pc_z[k] * msk_785[k];

        t_984[k] = f_6 * msi0_610[k]
                   - f_7 * msi1_610[k]
                   + f_3 * pc_z[k] * msk_786[k];

        t_985[k] = f_8 * msi0_611[k]
                   - f_9 * msi1_611[k]
                   + f_3 * pc_z[k] * msk_787[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, pc_y, pc_z, lsk_575, msi0_612, msi0_613, \
                         msi0_615, msi1_612, msi1_613, msi1_615, msk_788, msk_789, \
                         msk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_10 * msi0_612[k]
                   - f_11 * msi1_612[k]
                   + f_3 * pc_z[k] * msk_788[k];

        t_987[k] = f_12 * msi0_613[k]
                   - f_13 * msi1_613[k]
                   + f_3 * pc_z[k] * msk_789[k];

        t_988[k] = f_20 * lsk_575[k]
                   + f_3 * pc_y[k] * msk_791[k];

        t_989[k] = f_1 * msi0_615[k]
                   - f_2 * msi1_615[k]
                   + f_3 * pc_z[k] * msk_791[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, pa_z, pc_y, pc_z, lsl0_675, lsl0_678, \
                         lsk_540, lsk_576, lsl1_675, lsl1_678, \
                         msk_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pa_z[k] * lsl0_675[k]
                   - f_14 * pc_z[k] * lsl1_675[k];

        t_991[k] = f_19 * lsk_576[k]
                   + f_3 * pc_y[k] * msk_792[k];

        t_992[k] = f_15 * lsk_540[k]
                   + f_3 * pc_z[k] * msk_792[k];

        t_993[k] = pa_z[k] * lsl0_678[k]
                   - f_14 * pc_z[k] * lsl1_678[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, pa_z, pc_x, pc_y, pc_z, lsl0_681, lsk_578, \
                         lsk_797, lsl1_681, msi0_621, msi1_621, msk_794, \
                         msk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_19 * lsk_578[k]
                   + f_3 * pc_y[k] * msk_794[k];

        t_995[k] = f_17 * lsk_797[k]
                   + f_12 * msi0_621[k]
                   - f_13 * msi1_621[k]
                   + f_3 * pc_x[k] * msk_797[k];

        t_996[k] = pa_z[k] * lsl0_681[k]
                   - f_14 * pc_z[k] * lsl1_681[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, pc_x, pc_y, pc_z, lsk_543, lsk_581, lsk_801, \
                         msi0_625, msi1_625, msk_795, msk_797, \
                         msk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = f_15 * lsk_543[k]
                   + f_3 * pc_z[k] * msk_795[k];

        t_998[k] = f_19 * lsk_581[k]
                   + f_3 * pc_y[k] * msk_797[k];

        t_999[k] = f_17 * lsk_801[k]
                   + f_10 * msi0_625[k]
                   - f_11 * msi1_625[k]
                   + f_3 * pc_x[k] * msk_801[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pa_z, pc_y, pc_z, lsl0_685, lsl0_687, \
                         lsk_546, lsk_547, lsk_585, lsl1_685, lsl1_687, msk_798, \
                         msk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = pa_z[k] * lsl0_685[k]
                    - f_14 * pc_z[k] * lsl1_685[k];

        t_1001[k] = f_15 * lsk_546[k]
                    + f_3 * pc_z[k] * msk_798[k];

        t_1002[k] = pa_z[k] * lsl0_687[k]
                    + f_16 * lsk_547[k]
                    - f_14 * pc_z[k] * lsl1_687[k];

        t_1003[k] = f_19 * lsk_585[k]
                    + f_3 * pc_y[k] * msk_801[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, pa_z, pc_x, pc_z, lsl0_690, lsk_550, lsk_806, \
                         lsl1_690, msi0_630, msi1_630, msk_802, \
                         msk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_17 * lsk_806[k]
                    + f_8 * msi0_630[k]
                    - f_9 * msi1_630[k]
                    + f_3 * pc_x[k] * msk_806[k];

        t_1005[k] = pa_z[k] * lsl0_690[k]
                    - f_14 * pc_z[k] * lsl1_690[k];

        t_1006[k] = f_15 * lsk_550[k]
                    + f_3 * pc_z[k] * msk_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pa_z, pc_y, pc_z, lsl0_692, lsl0_693, \
                         lsk_551, lsk_552, lsk_590, lsl1_692, lsl1_693, \
                         msk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = pa_z[k] * lsl0_692[k]
                    + f_16 * lsk_551[k]
                    - f_14 * pc_z[k] * lsl1_692[k];

        t_1008[k] = pa_z[k] * lsl0_693[k]
                    + f_17 * lsk_552[k]
                    - f_14 * pc_z[k] * lsl1_693[k];

        t_1009[k] = f_19 * lsk_590[k]
                    + f_3 * pc_y[k] * msk_806[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pa_z, pc_x, pc_z, lsl0_696, lsk_555, lsk_812, \
                         lsl1_696, msi0_636, msi1_636, msk_807, \
                         msk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_17 * lsk_812[k]
                    + f_6 * msi0_636[k]
                    - f_7 * msi1_636[k]
                    + f_3 * pc_x[k] * msk_812[k];

        t_1011[k] = pa_z[k] * lsl0_696[k]
                    - f_14 * pc_z[k] * lsl1_696[k];

        t_1012[k] = f_15 * lsk_555[k]
                    + f_3 * pc_z[k] * msk_807[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pa_z, pc_z, lsl0_698, lsl0_699, lsl0_700, \
                         lsk_556, lsk_557, lsk_558, lsl1_698, lsl1_699, \
                         lsl1_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = pa_z[k] * lsl0_698[k]
                    + f_16 * lsk_556[k]
                    - f_14 * pc_z[k] * lsl1_698[k];

        t_1014[k] = pa_z[k] * lsl0_699[k]
                    + f_17 * lsk_557[k]
                    - f_14 * pc_z[k] * lsl1_699[k];

        t_1015[k] = pa_z[k] * lsl0_700[k]
                    + f_18 * lsk_558[k]
                    - f_14 * pc_z[k] * lsl1_700[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, lsk_596, lsk_819, \
                         lsk_820, lsk_821, msi0_643, msi1_643, msk_812, msk_819, msk_820, \
                         msk_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * lsk_596[k]
                    + f_3 * pc_y[k] * msk_812[k];

        t_1017[k] = f_17 * lsk_819[k]
                    + f_4 * msi0_643[k]
                    - f_5 * msi1_643[k]
                    + f_3 * pc_x[k] * msk_819[k];

        t_1018[k] = f_17 * lsk_820[k]
                    + f_3 * pc_x[k] * msk_820[k];

        t_1019[k] = f_17 * lsk_821[k]
                    + f_3 * pc_x[k] * msk_821[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, pc_x, lsk_822, lsk_823, \
                         lsk_824, lsk_825, lsk_826, msk_822, msk_823, msk_824, msk_825, \
                         msk_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_17 * lsk_822[k]
                    + f_3 * pc_x[k] * msk_822[k];

        t_1021[k] = f_17 * lsk_823[k]
                    + f_3 * pc_x[k] * msk_823[k];

        t_1022[k] = f_17 * lsk_824[k]
                    + f_3 * pc_x[k] * msk_824[k];

        t_1023[k] = f_17 * lsk_825[k]
                    + f_3 * pc_x[k] * msk_825[k];

        t_1024[k] = f_17 * lsk_826[k]
                    + f_3 * pc_x[k] * msk_826[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pa_z, pc_x, pc_z, lsl0_711, lsk_568, lsk_827, \
                         lsl1_711, msk_820, msk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_17 * lsk_827[k]
                    + f_3 * pc_x[k] * msk_827[k];

        t_1026[k] = pa_z[k] * lsl0_711[k]
                    - f_14 * pc_z[k] * lsl1_711[k];

        t_1027[k] = f_15 * lsk_568[k]
                    + f_3 * pc_z[k] * msk_820[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, pc_y, lsk_606, lsk_607, lsk_608, msi0_639, \
                         msi0_640, msi0_641, msi1_639, msi1_640, msi1_641, msk_822, msk_823, \
                         msk_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_19 * lsk_606[k]
                    + f_12 * msi0_639[k]
                    - f_13 * msi1_639[k]
                    + f_3 * pc_y[k] * msk_822[k];

        t_1029[k] = f_19 * lsk_607[k]
                    + f_10 * msi0_640[k]
                    - f_11 * msi1_640[k]
                    + f_3 * pc_y[k] * msk_823[k];

        t_1030[k] = f_19 * lsk_608[k]
                    + f_8 * msi0_641[k]
                    - f_9 * msi1_641[k]
                    + f_3 * pc_y[k] * msk_824[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pc_y, lsk_609, lsk_610, lsk_611, msi0_642, \
                         msi0_643, msi1_642, msi1_643, msk_825, msk_826, \
                         msk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_19 * lsk_609[k]
                    + f_6 * msi0_642[k]
                    - f_7 * msi1_642[k]
                    + f_3 * pc_y[k] * msk_825[k];

        t_1032[k] = f_19 * lsk_610[k]
                    + f_4 * msi0_643[k]
                    - f_5 * msi1_643[k]
                    + f_3 * pc_y[k] * msk_826[k];

        t_1033[k] = f_19 * lsk_611[k]
                    + f_3 * pc_y[k] * msk_827[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_y, pc_z, lsk_575, lsk_612, lsk_828, \
                         msi0_643, msi0_644, msi1_643, msi1_644, msk_827, \
                         msk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_15 * lsk_575[k]
                    + f_1 * msi0_643[k]
                    - f_2 * msi1_643[k]
                    + f_3 * pc_z[k] * msk_827[k];

        t_1035[k] = f_17 * lsk_828[k]
                    + f_1 * msi0_644[k]
                    - f_2 * msi1_644[k]
                    + f_3 * pc_x[k] * msk_828[k];

        t_1036[k] = f_18 * lsk_612[k]
                    + f_3 * pc_y[k] * msk_828[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t lsk, const size_t msi0,
                                                          const size_t msi1, const size_t msk,
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

    const auto *lsk_576 = buffer.data(lsk + 576);
    const auto *lsk_579 = buffer.data(lsk + 579);
    const auto *lsk_582 = buffer.data(lsk + 582);
    const auto *lsk_586 = buffer.data(lsk + 586);
    const auto *lsk_591 = buffer.data(lsk + 591);
    const auto *lsk_604 = buffer.data(lsk + 604);
    const auto *lsk_611 = buffer.data(lsk + 611);
    const auto *lsk_612 = buffer.data(lsk + 612);
    const auto *lsk_614 = buffer.data(lsk + 614);
    const auto *lsk_615 = buffer.data(lsk + 615);
    const auto *lsk_617 = buffer.data(lsk + 617);
    const auto *lsk_618 = buffer.data(lsk + 618);
    const auto *lsk_621 = buffer.data(lsk + 621);
    const auto *lsk_622 = buffer.data(lsk + 622);
    const auto *lsk_626 = buffer.data(lsk + 626);
    const auto *lsk_627 = buffer.data(lsk + 627);
    const auto *lsk_632 = buffer.data(lsk + 632);
    const auto *lsk_640 = buffer.data(lsk + 640);
    const auto *lsk_642 = buffer.data(lsk + 642);
    const auto *lsk_643 = buffer.data(lsk + 643);
    const auto *lsk_644 = buffer.data(lsk + 644);
    const auto *lsk_645 = buffer.data(lsk + 645);
    const auto *lsk_646 = buffer.data(lsk + 646);
    const auto *lsk_647 = buffer.data(lsk + 647);
    const auto *lsk_648 = buffer.data(lsk + 648);
    const auto *lsk_650 = buffer.data(lsk + 650);
    const auto *lsk_651 = buffer.data(lsk + 651);
    const auto *lsk_653 = buffer.data(lsk + 653);
    const auto *lsk_654 = buffer.data(lsk + 654);
    const auto *lsk_657 = buffer.data(lsk + 657);
    const auto *lsk_662 = buffer.data(lsk + 662);
    const auto *lsk_668 = buffer.data(lsk + 668);
    const auto *lsk_676 = buffer.data(lsk + 676);
    const auto *lsk_678 = buffer.data(lsk + 678);
    const auto *lsk_679 = buffer.data(lsk + 679);
    const auto *lsk_680 = buffer.data(lsk + 680);
    const auto *lsk_681 = buffer.data(lsk + 681);
    const auto *lsk_682 = buffer.data(lsk + 682);
    const auto *lsk_683 = buffer.data(lsk + 683);
    const auto *lsk_684 = buffer.data(lsk + 684);
    const auto *lsk_686 = buffer.data(lsk + 686);
    const auto *lsk_689 = buffer.data(lsk + 689);
    const auto *lsk_693 = buffer.data(lsk + 693);
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
    const auto *lsk_903 = buffer.data(lsk + 903);
    const auto *lsk_905 = buffer.data(lsk + 905);
    const auto *lsk_906 = buffer.data(lsk + 906);
    const auto *lsk_909 = buffer.data(lsk + 909);
    const auto *lsk_910 = buffer.data(lsk + 910);
    const auto *lsk_912 = buffer.data(lsk + 912);

    const auto *msi0_647 = buffer.data(msi0 + 647);
    const auto *msi0_649 = buffer.data(msi0 + 649);
    const auto *msi0_650 = buffer.data(msi0 + 650);
    const auto *msi0_653 = buffer.data(msi0 + 653);
    const auto *msi0_654 = buffer.data(msi0 + 654);
    const auto *msi0_656 = buffer.data(msi0 + 656);
    const auto *msi0_658 = buffer.data(msi0 + 658);
    const auto *msi0_659 = buffer.data(msi0 + 659);
    const auto *msi0_661 = buffer.data(msi0 + 661);
    const auto *msi0_662 = buffer.data(msi0 + 662);
    const auto *msi0_664 = buffer.data(msi0 + 664);
    const auto *msi0_665 = buffer.data(msi0 + 665);
    const auto *msi0_667 = buffer.data(msi0 + 667);
    const auto *msi0_668 = buffer.data(msi0 + 668);
    const auto *msi0_669 = buffer.data(msi0 + 669);
    const auto *msi0_670 = buffer.data(msi0 + 670);
    const auto *msi0_671 = buffer.data(msi0 + 671);
    const auto *msi0_672 = buffer.data(msi0 + 672);
    const auto *msi0_675 = buffer.data(msi0 + 675);
    const auto *msi0_677 = buffer.data(msi0 + 677);
    const auto *msi0_678 = buffer.data(msi0 + 678);
    const auto *msi0_681 = buffer.data(msi0 + 681);
    const auto *msi0_682 = buffer.data(msi0 + 682);
    const auto *msi0_684 = buffer.data(msi0 + 684);
    const auto *msi0_686 = buffer.data(msi0 + 686);
    const auto *msi0_687 = buffer.data(msi0 + 687);
    const auto *msi0_689 = buffer.data(msi0 + 689);
    const auto *msi0_690 = buffer.data(msi0 + 690);
    const auto *msi0_692 = buffer.data(msi0 + 692);
    const auto *msi0_693 = buffer.data(msi0 + 693);
    const auto *msi0_695 = buffer.data(msi0 + 695);
    const auto *msi0_696 = buffer.data(msi0 + 696);
    const auto *msi0_697 = buffer.data(msi0 + 697);
    const auto *msi0_698 = buffer.data(msi0 + 698);
    const auto *msi0_699 = buffer.data(msi0 + 699);
    const auto *msi0_700 = buffer.data(msi0 + 700);
    const auto *msi0_703 = buffer.data(msi0 + 703);
    const auto *msi0_705 = buffer.data(msi0 + 705);
    const auto *msi0_706 = buffer.data(msi0 + 706);
    const auto *msi0_709 = buffer.data(msi0 + 709);
    const auto *msi0_710 = buffer.data(msi0 + 710);
    const auto *msi0_712 = buffer.data(msi0 + 712);

    const auto *msi1_647 = buffer.data(msi1 + 647);
    const auto *msi1_649 = buffer.data(msi1 + 649);
    const auto *msi1_650 = buffer.data(msi1 + 650);
    const auto *msi1_653 = buffer.data(msi1 + 653);
    const auto *msi1_654 = buffer.data(msi1 + 654);
    const auto *msi1_656 = buffer.data(msi1 + 656);
    const auto *msi1_658 = buffer.data(msi1 + 658);
    const auto *msi1_659 = buffer.data(msi1 + 659);
    const auto *msi1_661 = buffer.data(msi1 + 661);
    const auto *msi1_662 = buffer.data(msi1 + 662);
    const auto *msi1_664 = buffer.data(msi1 + 664);
    const auto *msi1_665 = buffer.data(msi1 + 665);
    const auto *msi1_667 = buffer.data(msi1 + 667);
    const auto *msi1_668 = buffer.data(msi1 + 668);
    const auto *msi1_669 = buffer.data(msi1 + 669);
    const auto *msi1_670 = buffer.data(msi1 + 670);
    const auto *msi1_671 = buffer.data(msi1 + 671);
    const auto *msi1_672 = buffer.data(msi1 + 672);
    const auto *msi1_675 = buffer.data(msi1 + 675);
    const auto *msi1_677 = buffer.data(msi1 + 677);
    const auto *msi1_678 = buffer.data(msi1 + 678);
    const auto *msi1_681 = buffer.data(msi1 + 681);
    const auto *msi1_682 = buffer.data(msi1 + 682);
    const auto *msi1_684 = buffer.data(msi1 + 684);
    const auto *msi1_686 = buffer.data(msi1 + 686);
    const auto *msi1_687 = buffer.data(msi1 + 687);
    const auto *msi1_689 = buffer.data(msi1 + 689);
    const auto *msi1_690 = buffer.data(msi1 + 690);
    const auto *msi1_692 = buffer.data(msi1 + 692);
    const auto *msi1_693 = buffer.data(msi1 + 693);
    const auto *msi1_695 = buffer.data(msi1 + 695);
    const auto *msi1_696 = buffer.data(msi1 + 696);
    const auto *msi1_697 = buffer.data(msi1 + 697);
    const auto *msi1_698 = buffer.data(msi1 + 698);
    const auto *msi1_699 = buffer.data(msi1 + 699);
    const auto *msi1_700 = buffer.data(msi1 + 700);
    const auto *msi1_703 = buffer.data(msi1 + 703);
    const auto *msi1_705 = buffer.data(msi1 + 705);
    const auto *msi1_706 = buffer.data(msi1 + 706);
    const auto *msi1_709 = buffer.data(msi1 + 709);
    const auto *msi1_710 = buffer.data(msi1 + 710);
    const auto *msi1_712 = buffer.data(msi1 + 712);

    const auto *msk_828 = buffer.data(msk + 828);
    const auto *msk_830 = buffer.data(msk + 830);
    const auto *msk_831 = buffer.data(msk + 831);
    const auto *msk_833 = buffer.data(msk + 833);
    const auto *msk_834 = buffer.data(msk + 834);
    const auto *msk_837 = buffer.data(msk + 837);
    const auto *msk_838 = buffer.data(msk + 838);
    const auto *msk_840 = buffer.data(msk + 840);
    const auto *msk_842 = buffer.data(msk + 842);
    const auto *msk_843 = buffer.data(msk + 843);
    const auto *msk_845 = buffer.data(msk + 845);
    const auto *msk_846 = buffer.data(msk + 846);
    const auto *msk_848 = buffer.data(msk + 848);
    const auto *msk_849 = buffer.data(msk + 849);
    const auto *msk_851 = buffer.data(msk + 851);
    const auto *msk_852 = buffer.data(msk + 852);
    const auto *msk_853 = buffer.data(msk + 853);
    const auto *msk_855 = buffer.data(msk + 855);
    const auto *msk_856 = buffer.data(msk + 856);
    const auto *msk_857 = buffer.data(msk + 857);
    const auto *msk_858 = buffer.data(msk + 858);
    const auto *msk_859 = buffer.data(msk + 859);
    const auto *msk_860 = buffer.data(msk + 860);
    const auto *msk_861 = buffer.data(msk + 861);
    const auto *msk_862 = buffer.data(msk + 862);
    const auto *msk_863 = buffer.data(msk + 863);
    const auto *msk_864 = buffer.data(msk + 864);
    const auto *msk_866 = buffer.data(msk + 866);
    const auto *msk_867 = buffer.data(msk + 867);
    const auto *msk_869 = buffer.data(msk + 869);
    const auto *msk_870 = buffer.data(msk + 870);
    const auto *msk_873 = buffer.data(msk + 873);
    const auto *msk_874 = buffer.data(msk + 874);
    const auto *msk_876 = buffer.data(msk + 876);
    const auto *msk_878 = buffer.data(msk + 878);
    const auto *msk_879 = buffer.data(msk + 879);
    const auto *msk_881 = buffer.data(msk + 881);
    const auto *msk_882 = buffer.data(msk + 882);
    const auto *msk_884 = buffer.data(msk + 884);
    const auto *msk_885 = buffer.data(msk + 885);
    const auto *msk_887 = buffer.data(msk + 887);
    const auto *msk_888 = buffer.data(msk + 888);
    const auto *msk_889 = buffer.data(msk + 889);
    const auto *msk_891 = buffer.data(msk + 891);
    const auto *msk_892 = buffer.data(msk + 892);
    const auto *msk_893 = buffer.data(msk + 893);
    const auto *msk_894 = buffer.data(msk + 894);
    const auto *msk_895 = buffer.data(msk + 895);
    const auto *msk_896 = buffer.data(msk + 896);
    const auto *msk_897 = buffer.data(msk + 897);
    const auto *msk_898 = buffer.data(msk + 898);
    const auto *msk_899 = buffer.data(msk + 899);
    const auto *msk_900 = buffer.data(msk + 900);
    const auto *msk_902 = buffer.data(msk + 902);
    const auto *msk_903 = buffer.data(msk + 903);
    const auto *msk_905 = buffer.data(msk + 905);
    const auto *msk_906 = buffer.data(msk + 906);
    const auto *msk_909 = buffer.data(msk + 909);
    const auto *msk_910 = buffer.data(msk + 910);
    const auto *msk_912 = buffer.data(msk + 912);

#pragma omp simd aligned(t_1037, t_1038, t_1039, pc_x, pc_y, pc_z, lsk_576, lsk_614, lsk_831, \
                         msi0_647, msi1_647, msk_828, msk_830, \
                         msk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_16 * lsk_576[k]
                    + f_3 * pc_z[k] * msk_828[k];

        t_1038[k] = f_17 * lsk_831[k]
                    + f_12 * msi0_647[k]
                    - f_13 * msi1_647[k]
                    + f_3 * pc_x[k] * msk_831[k];

        t_1039[k] = f_18 * lsk_614[k]
                    + f_3 * pc_y[k] * msk_830[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pc_x, pc_z, lsk_579, lsk_833, lsk_834, \
                         msi0_649, msi0_650, msi1_649, msi1_650, msk_831, msk_833, \
                         msk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_17 * lsk_833[k]
                    + f_12 * msi0_649[k]
                    - f_13 * msi1_649[k]
                    + f_3 * pc_x[k] * msk_833[k];

        t_1041[k] = f_17 * lsk_834[k]
                    + f_10 * msi0_650[k]
                    - f_11 * msi1_650[k]
                    + f_3 * pc_x[k] * msk_834[k];

        t_1042[k] = f_16 * lsk_579[k]
                    + f_3 * pc_z[k] * msk_831[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, lsk_617, lsk_837, lsk_838, \
                         msi0_653, msi0_654, msi1_653, msi1_654, msk_833, msk_837, \
                         msk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_18 * lsk_617[k]
                    + f_3 * pc_y[k] * msk_833[k];

        t_1044[k] = f_17 * lsk_837[k]
                    + f_10 * msi0_653[k]
                    - f_11 * msi1_653[k]
                    + f_3 * pc_x[k] * msk_837[k];

        t_1045[k] = f_17 * lsk_838[k]
                    + f_8 * msi0_654[k]
                    - f_9 * msi1_654[k]
                    + f_3 * pc_x[k] * msk_838[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pc_x, pc_y, pc_z, lsk_582, lsk_621, lsk_840, \
                         msi0_656, msi1_656, msk_834, msk_837, \
                         msk_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_16 * lsk_582[k]
                    + f_3 * pc_z[k] * msk_834[k];

        t_1047[k] = f_17 * lsk_840[k]
                    + f_8 * msi0_656[k]
                    - f_9 * msi1_656[k]
                    + f_3 * pc_x[k] * msk_840[k];

        t_1048[k] = f_18 * lsk_621[k]
                    + f_3 * pc_y[k] * msk_837[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pc_x, pc_z, lsk_586, lsk_842, lsk_843, \
                         msi0_658, msi0_659, msi1_658, msi1_659, msk_838, msk_842, \
                         msk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_17 * lsk_842[k]
                    + f_8 * msi0_658[k]
                    - f_9 * msi1_658[k]
                    + f_3 * pc_x[k] * msk_842[k];

        t_1050[k] = f_17 * lsk_843[k]
                    + f_6 * msi0_659[k]
                    - f_7 * msi1_659[k]
                    + f_3 * pc_x[k] * msk_843[k];

        t_1051[k] = f_16 * lsk_586[k]
                    + f_3 * pc_z[k] * msk_838[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pc_x, pc_y, lsk_626, lsk_845, lsk_846, \
                         msi0_661, msi0_662, msi1_661, msi1_662, msk_842, msk_845, \
                         msk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_17 * lsk_845[k]
                    + f_6 * msi0_661[k]
                    - f_7 * msi1_661[k]
                    + f_3 * pc_x[k] * msk_845[k];

        t_1053[k] = f_17 * lsk_846[k]
                    + f_6 * msi0_662[k]
                    - f_7 * msi1_662[k]
                    + f_3 * pc_x[k] * msk_846[k];

        t_1054[k] = f_18 * lsk_626[k]
                    + f_3 * pc_y[k] * msk_842[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, lsk_591, lsk_848, lsk_849, \
                         msi0_664, msi0_665, msi1_664, msi1_665, msk_843, msk_848, \
                         msk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_17 * lsk_848[k]
                    + f_6 * msi0_664[k]
                    - f_7 * msi1_664[k]
                    + f_3 * pc_x[k] * msk_848[k];

        t_1056[k] = f_17 * lsk_849[k]
                    + f_4 * msi0_665[k]
                    - f_5 * msi1_665[k]
                    + f_3 * pc_x[k] * msk_849[k];

        t_1057[k] = f_16 * lsk_591[k]
                    + f_3 * pc_z[k] * msk_843[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_x, lsk_851, lsk_852, lsk_853, msi0_667, \
                         msi0_668, msi0_669, msi1_667, msi1_668, msi1_669, msk_851, msk_852, \
                         msk_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_17 * lsk_851[k]
                    + f_4 * msi0_667[k]
                    - f_5 * msi1_667[k]
                    + f_3 * pc_x[k] * msk_851[k];

        t_1059[k] = f_17 * lsk_852[k]
                    + f_4 * msi0_668[k]
                    - f_5 * msi1_668[k]
                    + f_3 * pc_x[k] * msk_852[k];

        t_1060[k] = f_17 * lsk_853[k]
                    + f_4 * msi0_669[k]
                    - f_5 * msi1_669[k]
                    + f_3 * pc_x[k] * msk_853[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pc_x, pc_y, lsk_632, lsk_855, \
                         lsk_856, lsk_857, msi0_671, msi1_671, msk_848, msk_855, msk_856, \
                         msk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_18 * lsk_632[k]
                    + f_3 * pc_y[k] * msk_848[k];

        t_1062[k] = f_17 * lsk_855[k]
                    + f_4 * msi0_671[k]
                    - f_5 * msi1_671[k]
                    + f_3 * pc_x[k] * msk_855[k];

        t_1063[k] = f_17 * lsk_856[k]
                    + f_3 * pc_x[k] * msk_856[k];

        t_1064[k] = f_17 * lsk_857[k]
                    + f_3 * pc_x[k] * msk_857[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, lsk_858, lsk_859, \
                         lsk_860, lsk_861, lsk_862, msk_858, msk_859, msk_860, msk_861, \
                         msk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_17 * lsk_858[k]
                    + f_3 * pc_x[k] * msk_858[k];

        t_1066[k] = f_17 * lsk_859[k]
                    + f_3 * pc_x[k] * msk_859[k];

        t_1067[k] = f_17 * lsk_860[k]
                    + f_3 * pc_x[k] * msk_860[k];

        t_1068[k] = f_17 * lsk_861[k]
                    + f_3 * pc_x[k] * msk_861[k];

        t_1069[k] = f_17 * lsk_862[k]
                    + f_3 * pc_x[k] * msk_862[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, lsk_604, lsk_640, lsk_863, \
                         msi0_665, msi1_665, msk_856, msk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_17 * lsk_863[k]
                    + f_3 * pc_x[k] * msk_863[k];

        t_1071[k] = f_18 * lsk_640[k]
                    + f_1 * msi0_665[k]
                    - f_2 * msi1_665[k]
                    + f_3 * pc_y[k] * msk_856[k];

        t_1072[k] = f_16 * lsk_604[k]
                    + f_3 * pc_z[k] * msk_856[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_y, lsk_642, lsk_643, lsk_644, msi0_667, \
                         msi0_668, msi0_669, msi1_667, msi1_668, msi1_669, msk_858, msk_859, \
                         msk_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_18 * lsk_642[k]
                    + f_12 * msi0_667[k]
                    - f_13 * msi1_667[k]
                    + f_3 * pc_y[k] * msk_858[k];

        t_1074[k] = f_18 * lsk_643[k]
                    + f_10 * msi0_668[k]
                    - f_11 * msi1_668[k]
                    + f_3 * pc_y[k] * msk_859[k];

        t_1075[k] = f_18 * lsk_644[k]
                    + f_8 * msi0_669[k]
                    - f_9 * msi1_669[k]
                    + f_3 * pc_y[k] * msk_860[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_y, lsk_645, lsk_646, lsk_647, msi0_670, \
                         msi0_671, msi1_670, msi1_671, msk_861, msk_862, \
                         msk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_18 * lsk_645[k]
                    + f_6 * msi0_670[k]
                    - f_7 * msi1_670[k]
                    + f_3 * pc_y[k] * msk_861[k];

        t_1077[k] = f_18 * lsk_646[k]
                    + f_4 * msi0_671[k]
                    - f_5 * msi1_671[k]
                    + f_3 * pc_y[k] * msk_862[k];

        t_1078[k] = f_18 * lsk_647[k]
                    + f_3 * pc_y[k] * msk_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, pc_x, pc_y, pc_z, lsk_611, lsk_648, lsk_864, \
                         msi0_671, msi0_672, msi1_671, msi1_672, msk_863, \
                         msk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_16 * lsk_611[k]
                    + f_1 * msi0_671[k]
                    - f_2 * msi1_671[k]
                    + f_3 * pc_z[k] * msk_863[k];

        t_1080[k] = f_17 * lsk_864[k]
                    + f_1 * msi0_672[k]
                    - f_2 * msi1_672[k]
                    + f_3 * pc_x[k] * msk_864[k];

        t_1081[k] = f_17 * lsk_648[k]
                    + f_3 * pc_y[k] * msk_864[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, lsk_612, lsk_650, lsk_867, \
                         msi0_675, msi1_675, msk_864, msk_866, \
                         msk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_17 * lsk_612[k]
                    + f_3 * pc_z[k] * msk_864[k];

        t_1083[k] = f_17 * lsk_867[k]
                    + f_12 * msi0_675[k]
                    - f_13 * msi1_675[k]
                    + f_3 * pc_x[k] * msk_867[k];

        t_1084[k] = f_17 * lsk_650[k]
                    + f_3 * pc_y[k] * msk_866[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, lsk_615, lsk_869, lsk_870, \
                         msi0_677, msi0_678, msi1_677, msi1_678, msk_867, msk_869, \
                         msk_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_17 * lsk_869[k]
                    + f_12 * msi0_677[k]
                    - f_13 * msi1_677[k]
                    + f_3 * pc_x[k] * msk_869[k];

        t_1086[k] = f_17 * lsk_870[k]
                    + f_10 * msi0_678[k]
                    - f_11 * msi1_678[k]
                    + f_3 * pc_x[k] * msk_870[k];

        t_1087[k] = f_17 * lsk_615[k]
                    + f_3 * pc_z[k] * msk_867[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, lsk_653, lsk_873, lsk_874, \
                         msi0_681, msi0_682, msi1_681, msi1_682, msk_869, msk_873, \
                         msk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * lsk_653[k]
                    + f_3 * pc_y[k] * msk_869[k];

        t_1089[k] = f_17 * lsk_873[k]
                    + f_10 * msi0_681[k]
                    - f_11 * msi1_681[k]
                    + f_3 * pc_x[k] * msk_873[k];

        t_1090[k] = f_17 * lsk_874[k]
                    + f_8 * msi0_682[k]
                    - f_9 * msi1_682[k]
                    + f_3 * pc_x[k] * msk_874[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, lsk_618, lsk_657, lsk_876, \
                         msi0_684, msi1_684, msk_870, msk_873, \
                         msk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_17 * lsk_618[k]
                    + f_3 * pc_z[k] * msk_870[k];

        t_1092[k] = f_17 * lsk_876[k]
                    + f_8 * msi0_684[k]
                    - f_9 * msi1_684[k]
                    + f_3 * pc_x[k] * msk_876[k];

        t_1093[k] = f_17 * lsk_657[k]
                    + f_3 * pc_y[k] * msk_873[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, lsk_622, lsk_878, lsk_879, \
                         msi0_686, msi0_687, msi1_686, msi1_687, msk_874, msk_878, \
                         msk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_17 * lsk_878[k]
                    + f_8 * msi0_686[k]
                    - f_9 * msi1_686[k]
                    + f_3 * pc_x[k] * msk_878[k];

        t_1095[k] = f_17 * lsk_879[k]
                    + f_6 * msi0_687[k]
                    - f_7 * msi1_687[k]
                    + f_3 * pc_x[k] * msk_879[k];

        t_1096[k] = f_17 * lsk_622[k]
                    + f_3 * pc_z[k] * msk_874[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, lsk_662, lsk_881, lsk_882, \
                         msi0_689, msi0_690, msi1_689, msi1_690, msk_878, msk_881, \
                         msk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_17 * lsk_881[k]
                    + f_6 * msi0_689[k]
                    - f_7 * msi1_689[k]
                    + f_3 * pc_x[k] * msk_881[k];

        t_1098[k] = f_17 * lsk_882[k]
                    + f_6 * msi0_690[k]
                    - f_7 * msi1_690[k]
                    + f_3 * pc_x[k] * msk_882[k];

        t_1099[k] = f_17 * lsk_662[k]
                    + f_3 * pc_y[k] * msk_878[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_z, lsk_627, lsk_884, lsk_885, \
                         msi0_692, msi0_693, msi1_692, msi1_693, msk_879, msk_884, \
                         msk_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_17 * lsk_884[k]
                    + f_6 * msi0_692[k]
                    - f_7 * msi1_692[k]
                    + f_3 * pc_x[k] * msk_884[k];

        t_1101[k] = f_17 * lsk_885[k]
                    + f_4 * msi0_693[k]
                    - f_5 * msi1_693[k]
                    + f_3 * pc_x[k] * msk_885[k];

        t_1102[k] = f_17 * lsk_627[k]
                    + f_3 * pc_z[k] * msk_879[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, lsk_887, lsk_888, lsk_889, msi0_695, \
                         msi0_696, msi0_697, msi1_695, msi1_696, msi1_697, msk_887, msk_888, \
                         msk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_17 * lsk_887[k]
                    + f_4 * msi0_695[k]
                    - f_5 * msi1_695[k]
                    + f_3 * pc_x[k] * msk_887[k];

        t_1104[k] = f_17 * lsk_888[k]
                    + f_4 * msi0_696[k]
                    - f_5 * msi1_696[k]
                    + f_3 * pc_x[k] * msk_888[k];

        t_1105[k] = f_17 * lsk_889[k]
                    + f_4 * msi0_697[k]
                    - f_5 * msi1_697[k]
                    + f_3 * pc_x[k] * msk_889[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, pc_y, lsk_668, lsk_891, \
                         lsk_892, lsk_893, msi0_699, msi1_699, msk_884, msk_891, msk_892, \
                         msk_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_17 * lsk_668[k]
                    + f_3 * pc_y[k] * msk_884[k];

        t_1107[k] = f_17 * lsk_891[k]
                    + f_4 * msi0_699[k]
                    - f_5 * msi1_699[k]
                    + f_3 * pc_x[k] * msk_891[k];

        t_1108[k] = f_17 * lsk_892[k]
                    + f_3 * pc_x[k] * msk_892[k];

        t_1109[k] = f_17 * lsk_893[k]
                    + f_3 * pc_x[k] * msk_893[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pc_x, lsk_894, lsk_895, \
                         lsk_896, lsk_897, lsk_898, msk_894, msk_895, msk_896, msk_897, \
                         msk_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_17 * lsk_894[k]
                    + f_3 * pc_x[k] * msk_894[k];

        t_1111[k] = f_17 * lsk_895[k]
                    + f_3 * pc_x[k] * msk_895[k];

        t_1112[k] = f_17 * lsk_896[k]
                    + f_3 * pc_x[k] * msk_896[k];

        t_1113[k] = f_17 * lsk_897[k]
                    + f_3 * pc_x[k] * msk_897[k];

        t_1114[k] = f_17 * lsk_898[k]
                    + f_3 * pc_x[k] * msk_898[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, lsk_640, lsk_676, lsk_899, \
                         msi0_693, msi1_693, msk_892, msk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_17 * lsk_899[k]
                    + f_3 * pc_x[k] * msk_899[k];

        t_1116[k] = f_17 * lsk_676[k]
                    + f_1 * msi0_693[k]
                    - f_2 * msi1_693[k]
                    + f_3 * pc_y[k] * msk_892[k];

        t_1117[k] = f_17 * lsk_640[k]
                    + f_3 * pc_z[k] * msk_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_y, lsk_678, lsk_679, lsk_680, msi0_695, \
                         msi0_696, msi0_697, msi1_695, msi1_696, msi1_697, msk_894, msk_895, \
                         msk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_17 * lsk_678[k]
                    + f_12 * msi0_695[k]
                    - f_13 * msi1_695[k]
                    + f_3 * pc_y[k] * msk_894[k];

        t_1119[k] = f_17 * lsk_679[k]
                    + f_10 * msi0_696[k]
                    - f_11 * msi1_696[k]
                    + f_3 * pc_y[k] * msk_895[k];

        t_1120[k] = f_17 * lsk_680[k]
                    + f_8 * msi0_697[k]
                    - f_9 * msi1_697[k]
                    + f_3 * pc_y[k] * msk_896[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_y, lsk_681, lsk_682, lsk_683, msi0_698, \
                         msi0_699, msi1_698, msi1_699, msk_897, msk_898, \
                         msk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_17 * lsk_681[k]
                    + f_6 * msi0_698[k]
                    - f_7 * msi1_698[k]
                    + f_3 * pc_y[k] * msk_897[k];

        t_1122[k] = f_17 * lsk_682[k]
                    + f_4 * msi0_699[k]
                    - f_5 * msi1_699[k]
                    + f_3 * pc_y[k] * msk_898[k];

        t_1123[k] = f_17 * lsk_683[k]
                    + f_3 * pc_y[k] * msk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, pc_z, lsk_647, lsk_684, lsk_900, \
                         msi0_699, msi0_700, msi1_699, msi1_700, msk_899, \
                         msk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_17 * lsk_647[k]
                    + f_1 * msi0_699[k]
                    - f_2 * msi1_699[k]
                    + f_3 * pc_z[k] * msk_899[k];

        t_1125[k] = f_17 * lsk_900[k]
                    + f_1 * msi0_700[k]
                    - f_2 * msi1_700[k]
                    + f_3 * pc_x[k] * msk_900[k];

        t_1126[k] = f_16 * lsk_684[k]
                    + f_3 * pc_y[k] * msk_900[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, lsk_648, lsk_686, lsk_903, \
                         msi0_703, msi1_703, msk_900, msk_902, \
                         msk_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_18 * lsk_648[k]
                    + f_3 * pc_z[k] * msk_900[k];

        t_1128[k] = f_17 * lsk_903[k]
                    + f_12 * msi0_703[k]
                    - f_13 * msi1_703[k]
                    + f_3 * pc_x[k] * msk_903[k];

        t_1129[k] = f_16 * lsk_686[k]
                    + f_3 * pc_y[k] * msk_902[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, lsk_651, lsk_905, lsk_906, \
                         msi0_705, msi0_706, msi1_705, msi1_706, msk_903, msk_905, \
                         msk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_17 * lsk_905[k]
                    + f_12 * msi0_705[k]
                    - f_13 * msi1_705[k]
                    + f_3 * pc_x[k] * msk_905[k];

        t_1131[k] = f_17 * lsk_906[k]
                    + f_10 * msi0_706[k]
                    - f_11 * msi1_706[k]
                    + f_3 * pc_x[k] * msk_906[k];

        t_1132[k] = f_18 * lsk_651[k]
                    + f_3 * pc_z[k] * msk_903[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, lsk_689, lsk_909, lsk_910, \
                         msi0_709, msi0_710, msi1_709, msi1_710, msk_905, msk_909, \
                         msk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_16 * lsk_689[k]
                    + f_3 * pc_y[k] * msk_905[k];

        t_1134[k] = f_17 * lsk_909[k]
                    + f_10 * msi0_709[k]
                    - f_11 * msi1_709[k]
                    + f_3 * pc_x[k] * msk_909[k];

        t_1135[k] = f_17 * lsk_910[k]
                    + f_8 * msi0_710[k]
                    - f_9 * msi1_710[k]
                    + f_3 * pc_x[k] * msk_910[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, lsk_654, lsk_693, lsk_912, \
                         msi0_712, msi1_712, msk_906, msk_909, \
                         msk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_18 * lsk_654[k]
                    + f_3 * pc_z[k] * msk_906[k];

        t_1137[k] = f_17 * lsk_912[k]
                    + f_8 * msi0_712[k]
                    - f_9 * msi1_712[k]
                    + f_3 * pc_x[k] * msk_912[k];

        t_1138[k] = f_16 * lsk_693[k]
                    + f_3 * pc_y[k] * msk_909[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
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

    const auto *lsl0_900 = buffer.data(lsl0 + 900);
    const auto *lsl0_903 = buffer.data(lsl0 + 903);
    const auto *lsl0_905 = buffer.data(lsl0 + 905);
    const auto *lsl0_906 = buffer.data(lsl0 + 906);
    const auto *lsl0_909 = buffer.data(lsl0 + 909);
    const auto *lsl0_910 = buffer.data(lsl0 + 910);
    const auto *lsl0_912 = buffer.data(lsl0 + 912);
    const auto *lsl0_914 = buffer.data(lsl0 + 914);
    const auto *lsl0_915 = buffer.data(lsl0 + 915);
    const auto *lsl0_917 = buffer.data(lsl0 + 917);
    const auto *lsl0_918 = buffer.data(lsl0 + 918);
    const auto *lsl0_920 = buffer.data(lsl0 + 920);
    const auto *lsl0_921 = buffer.data(lsl0 + 921);
    const auto *lsl0_923 = buffer.data(lsl0 + 923);
    const auto *lsl0_924 = buffer.data(lsl0 + 924);
    const auto *lsl0_925 = buffer.data(lsl0 + 925);
    const auto *lsl0_927 = buffer.data(lsl0 + 927);
    const auto *lsl0_944 = buffer.data(lsl0 + 944);

    const auto *lsk_658 = buffer.data(lsk + 658);
    const auto *lsk_663 = buffer.data(lsk + 663);
    const auto *lsk_676 = buffer.data(lsk + 676);
    const auto *lsk_683 = buffer.data(lsk + 683);
    const auto *lsk_684 = buffer.data(lsk + 684);
    const auto *lsk_687 = buffer.data(lsk + 687);
    const auto *lsk_690 = buffer.data(lsk + 690);
    const auto *lsk_694 = buffer.data(lsk + 694);
    const auto *lsk_698 = buffer.data(lsk + 698);
    const auto *lsk_699 = buffer.data(lsk + 699);
    const auto *lsk_704 = buffer.data(lsk + 704);
    const auto *lsk_712 = buffer.data(lsk + 712);
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
    const auto *lsk_725 = buffer.data(lsk + 725);
    const auto *lsk_726 = buffer.data(lsk + 726);
    const auto *lsk_728 = buffer.data(lsk + 728);
    const auto *lsk_729 = buffer.data(lsk + 729);
    const auto *lsk_730 = buffer.data(lsk + 730);
    const auto *lsk_732 = buffer.data(lsk + 732);
    const auto *lsk_733 = buffer.data(lsk + 733);
    const auto *lsk_734 = buffer.data(lsk + 734);
    const auto *lsk_735 = buffer.data(lsk + 735);
    const auto *lsk_737 = buffer.data(lsk + 737);
    const auto *lsk_738 = buffer.data(lsk + 738);
    const auto *lsk_739 = buffer.data(lsk + 739);
    const auto *lsk_740 = buffer.data(lsk + 740);
    const auto *lsk_748 = buffer.data(lsk + 748);
    const auto *lsk_750 = buffer.data(lsk + 750);
    const auto *lsk_751 = buffer.data(lsk + 751);
    const auto *lsk_752 = buffer.data(lsk + 752);
    const auto *lsk_753 = buffer.data(lsk + 753);
    const auto *lsk_754 = buffer.data(lsk + 754);
    const auto *lsk_755 = buffer.data(lsk + 755);
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
    const auto *lsk_964 = buffer.data(lsk + 964);
    const auto *lsk_965 = buffer.data(lsk + 965);
    const auto *lsk_966 = buffer.data(lsk + 966);
    const auto *lsk_967 = buffer.data(lsk + 967);
    const auto *lsk_968 = buffer.data(lsk + 968);
    const auto *lsk_969 = buffer.data(lsk + 969);
    const auto *lsk_970 = buffer.data(lsk + 970);
    const auto *lsk_971 = buffer.data(lsk + 971);
    const auto *lsk_972 = buffer.data(lsk + 972);
    const auto *lsk_977 = buffer.data(lsk + 977);
    const auto *lsk_981 = buffer.data(lsk + 981);
    const auto *lsk_986 = buffer.data(lsk + 986);
    const auto *lsk_992 = buffer.data(lsk + 992);
    const auto *lsk_999 = buffer.data(lsk + 999);
    const auto *lsk_1000 = buffer.data(lsk + 1000);
    const auto *lsk_1001 = buffer.data(lsk + 1001);
    const auto *lsk_1002 = buffer.data(lsk + 1002);
    const auto *lsk_1003 = buffer.data(lsk + 1003);
    const auto *lsk_1004 = buffer.data(lsk + 1004);
    const auto *lsk_1005 = buffer.data(lsk + 1005);
    const auto *lsk_1007 = buffer.data(lsk + 1007);

    const auto *lsl1_900 = buffer.data(lsl1 + 900);
    const auto *lsl1_903 = buffer.data(lsl1 + 903);
    const auto *lsl1_905 = buffer.data(lsl1 + 905);
    const auto *lsl1_906 = buffer.data(lsl1 + 906);
    const auto *lsl1_909 = buffer.data(lsl1 + 909);
    const auto *lsl1_910 = buffer.data(lsl1 + 910);
    const auto *lsl1_912 = buffer.data(lsl1 + 912);
    const auto *lsl1_914 = buffer.data(lsl1 + 914);
    const auto *lsl1_915 = buffer.data(lsl1 + 915);
    const auto *lsl1_917 = buffer.data(lsl1 + 917);
    const auto *lsl1_918 = buffer.data(lsl1 + 918);
    const auto *lsl1_920 = buffer.data(lsl1 + 920);
    const auto *lsl1_921 = buffer.data(lsl1 + 921);
    const auto *lsl1_923 = buffer.data(lsl1 + 923);
    const auto *lsl1_924 = buffer.data(lsl1 + 924);
    const auto *lsl1_925 = buffer.data(lsl1 + 925);
    const auto *lsl1_927 = buffer.data(lsl1 + 927);
    const auto *lsl1_944 = buffer.data(lsl1 + 944);

    const auto *msi0_714 = buffer.data(msi0 + 714);
    const auto *msi0_715 = buffer.data(msi0 + 715);
    const auto *msi0_717 = buffer.data(msi0 + 717);
    const auto *msi0_718 = buffer.data(msi0 + 718);
    const auto *msi0_720 = buffer.data(msi0 + 720);
    const auto *msi0_721 = buffer.data(msi0 + 721);
    const auto *msi0_723 = buffer.data(msi0 + 723);
    const auto *msi0_724 = buffer.data(msi0 + 724);
    const auto *msi0_725 = buffer.data(msi0 + 725);
    const auto *msi0_726 = buffer.data(msi0 + 726);
    const auto *msi0_727 = buffer.data(msi0 + 727);
    const auto *msi0_749 = buffer.data(msi0 + 749);
    const auto *msi0_751 = buffer.data(msi0 + 751);
    const auto *msi0_752 = buffer.data(msi0 + 752);
    const auto *msi0_753 = buffer.data(msi0 + 753);
    const auto *msi0_754 = buffer.data(msi0 + 754);
    const auto *msi0_755 = buffer.data(msi0 + 755);
    const auto *msi0_756 = buffer.data(msi0 + 756);
    const auto *msi0_757 = buffer.data(msi0 + 757);
    const auto *msi0_758 = buffer.data(msi0 + 758);
    const auto *msi0_759 = buffer.data(msi0 + 759);
    const auto *msi0_760 = buffer.data(msi0 + 760);
    const auto *msi0_761 = buffer.data(msi0 + 761);
    const auto *msi0_762 = buffer.data(msi0 + 762);
    const auto *msi0_763 = buffer.data(msi0 + 763);
    const auto *msi0_764 = buffer.data(msi0 + 764);
    const auto *msi0_765 = buffer.data(msi0 + 765);
    const auto *msi0_766 = buffer.data(msi0 + 766);
    const auto *msi0_767 = buffer.data(msi0 + 767);
    const auto *msi0_768 = buffer.data(msi0 + 768);
    const auto *msi0_769 = buffer.data(msi0 + 769);
    const auto *msi0_770 = buffer.data(msi0 + 770);
    const auto *msi0_776 = buffer.data(msi0 + 776);
    const auto *msi0_783 = buffer.data(msi0 + 783);

    const auto *msi1_714 = buffer.data(msi1 + 714);
    const auto *msi1_715 = buffer.data(msi1 + 715);
    const auto *msi1_717 = buffer.data(msi1 + 717);
    const auto *msi1_718 = buffer.data(msi1 + 718);
    const auto *msi1_720 = buffer.data(msi1 + 720);
    const auto *msi1_721 = buffer.data(msi1 + 721);
    const auto *msi1_723 = buffer.data(msi1 + 723);
    const auto *msi1_724 = buffer.data(msi1 + 724);
    const auto *msi1_725 = buffer.data(msi1 + 725);
    const auto *msi1_726 = buffer.data(msi1 + 726);
    const auto *msi1_727 = buffer.data(msi1 + 727);
    const auto *msi1_749 = buffer.data(msi1 + 749);
    const auto *msi1_751 = buffer.data(msi1 + 751);
    const auto *msi1_752 = buffer.data(msi1 + 752);
    const auto *msi1_753 = buffer.data(msi1 + 753);
    const auto *msi1_754 = buffer.data(msi1 + 754);
    const auto *msi1_755 = buffer.data(msi1 + 755);
    const auto *msi1_756 = buffer.data(msi1 + 756);
    const auto *msi1_757 = buffer.data(msi1 + 757);
    const auto *msi1_758 = buffer.data(msi1 + 758);
    const auto *msi1_759 = buffer.data(msi1 + 759);
    const auto *msi1_760 = buffer.data(msi1 + 760);
    const auto *msi1_761 = buffer.data(msi1 + 761);
    const auto *msi1_762 = buffer.data(msi1 + 762);
    const auto *msi1_763 = buffer.data(msi1 + 763);
    const auto *msi1_764 = buffer.data(msi1 + 764);
    const auto *msi1_765 = buffer.data(msi1 + 765);
    const auto *msi1_766 = buffer.data(msi1 + 766);
    const auto *msi1_767 = buffer.data(msi1 + 767);
    const auto *msi1_768 = buffer.data(msi1 + 768);
    const auto *msi1_769 = buffer.data(msi1 + 769);
    const auto *msi1_770 = buffer.data(msi1 + 770);
    const auto *msi1_776 = buffer.data(msi1 + 776);
    const auto *msi1_783 = buffer.data(msi1 + 783);

    const auto *msk_910 = buffer.data(msk + 910);
    const auto *msk_914 = buffer.data(msk + 914);
    const auto *msk_915 = buffer.data(msk + 915);
    const auto *msk_917 = buffer.data(msk + 917);
    const auto *msk_918 = buffer.data(msk + 918);
    const auto *msk_920 = buffer.data(msk + 920);
    const auto *msk_921 = buffer.data(msk + 921);
    const auto *msk_923 = buffer.data(msk + 923);
    const auto *msk_924 = buffer.data(msk + 924);
    const auto *msk_925 = buffer.data(msk + 925);
    const auto *msk_927 = buffer.data(msk + 927);
    const auto *msk_928 = buffer.data(msk + 928);
    const auto *msk_929 = buffer.data(msk + 929);
    const auto *msk_930 = buffer.data(msk + 930);
    const auto *msk_931 = buffer.data(msk + 931);
    const auto *msk_932 = buffer.data(msk + 932);
    const auto *msk_933 = buffer.data(msk + 933);
    const auto *msk_934 = buffer.data(msk + 934);
    const auto *msk_935 = buffer.data(msk + 935);
    const auto *msk_936 = buffer.data(msk + 936);
    const auto *msk_938 = buffer.data(msk + 938);
    const auto *msk_939 = buffer.data(msk + 939);
    const auto *msk_941 = buffer.data(msk + 941);
    const auto *msk_942 = buffer.data(msk + 942);
    const auto *msk_945 = buffer.data(msk + 945);
    const auto *msk_946 = buffer.data(msk + 946);
    const auto *msk_950 = buffer.data(msk + 950);
    const auto *msk_951 = buffer.data(msk + 951);
    const auto *msk_956 = buffer.data(msk + 956);
    const auto *msk_964 = buffer.data(msk + 964);
    const auto *msk_965 = buffer.data(msk + 965);
    const auto *msk_966 = buffer.data(msk + 966);
    const auto *msk_967 = buffer.data(msk + 967);
    const auto *msk_968 = buffer.data(msk + 968);
    const auto *msk_969 = buffer.data(msk + 969);
    const auto *msk_970 = buffer.data(msk + 970);
    const auto *msk_971 = buffer.data(msk + 971);
    const auto *msk_972 = buffer.data(msk + 972);
    const auto *msk_973 = buffer.data(msk + 973);
    const auto *msk_974 = buffer.data(msk + 974);
    const auto *msk_975 = buffer.data(msk + 975);
    const auto *msk_976 = buffer.data(msk + 976);
    const auto *msk_977 = buffer.data(msk + 977);
    const auto *msk_978 = buffer.data(msk + 978);
    const auto *msk_979 = buffer.data(msk + 979);
    const auto *msk_980 = buffer.data(msk + 980);
    const auto *msk_981 = buffer.data(msk + 981);
    const auto *msk_982 = buffer.data(msk + 982);
    const auto *msk_983 = buffer.data(msk + 983);
    const auto *msk_984 = buffer.data(msk + 984);
    const auto *msk_985 = buffer.data(msk + 985);
    const auto *msk_986 = buffer.data(msk + 986);
    const auto *msk_987 = buffer.data(msk + 987);
    const auto *msk_988 = buffer.data(msk + 988);
    const auto *msk_989 = buffer.data(msk + 989);
    const auto *msk_990 = buffer.data(msk + 990);
    const auto *msk_991 = buffer.data(msk + 991);
    const auto *msk_992 = buffer.data(msk + 992);
    const auto *msk_999 = buffer.data(msk + 999);
    const auto *msk_1000 = buffer.data(msk + 1000);
    const auto *msk_1001 = buffer.data(msk + 1001);
    const auto *msk_1002 = buffer.data(msk + 1002);
    const auto *msk_1003 = buffer.data(msk + 1003);
    const auto *msk_1004 = buffer.data(msk + 1004);
    const auto *msk_1005 = buffer.data(msk + 1005);
    const auto *msk_1007 = buffer.data(msk + 1007);

#pragma omp simd aligned(t_1139, t_1140, t_1141, pc_x, pc_z, lsk_658, lsk_914, lsk_915, \
                         msi0_714, msi0_715, msi1_714, msi1_715, msk_910, msk_914, \
                         msk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_17 * lsk_914[k]
                    + f_8 * msi0_714[k]
                    - f_9 * msi1_714[k]
                    + f_3 * pc_x[k] * msk_914[k];

        t_1140[k] = f_17 * lsk_915[k]
                    + f_6 * msi0_715[k]
                    - f_7 * msi1_715[k]
                    + f_3 * pc_x[k] * msk_915[k];

        t_1141[k] = f_18 * lsk_658[k]
                    + f_3 * pc_z[k] * msk_910[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pc_x, pc_y, lsk_698, lsk_917, lsk_918, \
                         msi0_717, msi0_718, msi1_717, msi1_718, msk_914, msk_917, \
                         msk_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_17 * lsk_917[k]
                    + f_6 * msi0_717[k]
                    - f_7 * msi1_717[k]
                    + f_3 * pc_x[k] * msk_917[k];

        t_1143[k] = f_17 * lsk_918[k]
                    + f_6 * msi0_718[k]
                    - f_7 * msi1_718[k]
                    + f_3 * pc_x[k] * msk_918[k];

        t_1144[k] = f_16 * lsk_698[k]
                    + f_3 * pc_y[k] * msk_914[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pc_x, pc_z, lsk_663, lsk_920, lsk_921, \
                         msi0_720, msi0_721, msi1_720, msi1_721, msk_915, msk_920, \
                         msk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_17 * lsk_920[k]
                    + f_6 * msi0_720[k]
                    - f_7 * msi1_720[k]
                    + f_3 * pc_x[k] * msk_920[k];

        t_1146[k] = f_17 * lsk_921[k]
                    + f_4 * msi0_721[k]
                    - f_5 * msi1_721[k]
                    + f_3 * pc_x[k] * msk_921[k];

        t_1147[k] = f_18 * lsk_663[k]
                    + f_3 * pc_z[k] * msk_915[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pc_x, lsk_923, lsk_924, lsk_925, msi0_723, \
                         msi0_724, msi0_725, msi1_723, msi1_724, msi1_725, msk_923, msk_924, \
                         msk_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_17 * lsk_923[k]
                    + f_4 * msi0_723[k]
                    - f_5 * msi1_723[k]
                    + f_3 * pc_x[k] * msk_923[k];

        t_1149[k] = f_17 * lsk_924[k]
                    + f_4 * msi0_724[k]
                    - f_5 * msi1_724[k]
                    + f_3 * pc_x[k] * msk_924[k];

        t_1150[k] = f_17 * lsk_925[k]
                    + f_4 * msi0_725[k]
                    - f_5 * msi1_725[k]
                    + f_3 * pc_x[k] * msk_925[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pc_x, pc_y, lsk_704, lsk_927, \
                         lsk_928, lsk_929, msi0_727, msi1_727, msk_920, msk_927, msk_928, \
                         msk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_16 * lsk_704[k]
                    + f_3 * pc_y[k] * msk_920[k];

        t_1152[k] = f_17 * lsk_927[k]
                    + f_4 * msi0_727[k]
                    - f_5 * msi1_727[k]
                    + f_3 * pc_x[k] * msk_927[k];

        t_1153[k] = f_17 * lsk_928[k]
                    + f_3 * pc_x[k] * msk_928[k];

        t_1154[k] = f_17 * lsk_929[k]
                    + f_3 * pc_x[k] * msk_929[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pc_x, lsk_930, lsk_931, \
                         lsk_932, lsk_933, lsk_934, msk_930, msk_931, msk_932, msk_933, \
                         msk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_17 * lsk_930[k]
                    + f_3 * pc_x[k] * msk_930[k];

        t_1156[k] = f_17 * lsk_931[k]
                    + f_3 * pc_x[k] * msk_931[k];

        t_1157[k] = f_17 * lsk_932[k]
                    + f_3 * pc_x[k] * msk_932[k];

        t_1158[k] = f_17 * lsk_933[k]
                    + f_3 * pc_x[k] * msk_933[k];

        t_1159[k] = f_17 * lsk_934[k]
                    + f_3 * pc_x[k] * msk_934[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, pc_z, lsk_676, lsk_712, lsk_935, \
                         msi0_721, msi1_721, msk_928, msk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_17 * lsk_935[k]
                    + f_3 * pc_x[k] * msk_935[k];

        t_1161[k] = f_16 * lsk_712[k]
                    + f_1 * msi0_721[k]
                    - f_2 * msi1_721[k]
                    + f_3 * pc_y[k] * msk_928[k];

        t_1162[k] = f_18 * lsk_676[k]
                    + f_3 * pc_z[k] * msk_928[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_y, lsk_714, lsk_715, lsk_716, msi0_723, \
                         msi0_724, msi0_725, msi1_723, msi1_724, msi1_725, msk_930, msk_931, \
                         msk_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * lsk_714[k]
                    + f_12 * msi0_723[k]
                    - f_13 * msi1_723[k]
                    + f_3 * pc_y[k] * msk_930[k];

        t_1164[k] = f_16 * lsk_715[k]
                    + f_10 * msi0_724[k]
                    - f_11 * msi1_724[k]
                    + f_3 * pc_y[k] * msk_931[k];

        t_1165[k] = f_16 * lsk_716[k]
                    + f_8 * msi0_725[k]
                    - f_9 * msi1_725[k]
                    + f_3 * pc_y[k] * msk_932[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_y, lsk_717, lsk_718, lsk_719, msi0_726, \
                         msi0_727, msi1_726, msi1_727, msk_933, msk_934, \
                         msk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_16 * lsk_717[k]
                    + f_6 * msi0_726[k]
                    - f_7 * msi1_726[k]
                    + f_3 * pc_y[k] * msk_933[k];

        t_1167[k] = f_16 * lsk_718[k]
                    + f_4 * msi0_727[k]
                    - f_5 * msi1_727[k]
                    + f_3 * pc_y[k] * msk_934[k];

        t_1168[k] = f_16 * lsk_719[k]
                    + f_3 * pc_y[k] * msk_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pa_y, pc_y, pc_z, lsl0_900, lsk_683, \
                         lsk_684, lsk_720, lsl1_900, msi0_727, msi1_727, msk_935, \
                         msk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_18 * lsk_683[k]
                    + f_1 * msi0_727[k]
                    - f_2 * msi1_727[k]
                    + f_3 * pc_z[k] * msk_935[k];

        t_1170[k] = pa_y[k] * lsl0_900[k]
                    - f_14 * pc_y[k] * lsl1_900[k];

        t_1171[k] = f_15 * lsk_720[k]
                    + f_3 * pc_y[k] * msk_936[k];

        t_1172[k] = f_19 * lsk_684[k]
                    + f_3 * pc_z[k] * msk_936[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, pa_y, pc_y, lsl0_903, lsl0_905, \
                         lsl0_906, lsk_721, lsk_722, lsk_723, lsl1_903, lsl1_905, lsl1_906, \
                         msk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = pa_y[k] * lsl0_903[k]
                    + f_16 * lsk_721[k]
                    - f_14 * pc_y[k] * lsl1_903[k];

        t_1174[k] = f_15 * lsk_722[k]
                    + f_3 * pc_y[k] * msk_938[k];

        t_1175[k] = pa_y[k] * lsl0_905[k]
                    - f_14 * pc_y[k] * lsl1_905[k];

        t_1176[k] = pa_y[k] * lsl0_906[k]
                    + f_17 * lsk_723[k]
                    - f_14 * pc_y[k] * lsl1_906[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, pa_y, pc_y, pc_z, lsl0_909, lsl0_910, \
                         lsk_687, lsk_725, lsk_726, lsl1_909, lsl1_910, msk_939, \
                         msk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_19 * lsk_687[k]
                    + f_3 * pc_z[k] * msk_939[k];

        t_1178[k] = f_15 * lsk_725[k]
                    + f_3 * pc_y[k] * msk_941[k];

        t_1179[k] = pa_y[k] * lsl0_909[k]
                    - f_14 * pc_y[k] * lsl1_909[k];

        t_1180[k] = pa_y[k] * lsl0_910[k]
                    + f_18 * lsk_726[k]
                    - f_14 * pc_y[k] * lsl1_910[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pa_y, pc_y, pc_z, lsl0_912, lsl0_914, \
                         lsk_690, lsk_728, lsk_729, lsl1_912, lsl1_914, msk_942, \
                         msk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_19 * lsk_690[k]
                    + f_3 * pc_z[k] * msk_942[k];

        t_1182[k] = pa_y[k] * lsl0_912[k]
                    + f_16 * lsk_728[k]
                    - f_14 * pc_y[k] * lsl1_912[k];

        t_1183[k] = f_15 * lsk_729[k]
                    + f_3 * pc_y[k] * msk_945[k];

        t_1184[k] = pa_y[k] * lsl0_914[k]
                    - f_14 * pc_y[k] * lsl1_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pa_y, pc_y, pc_z, lsl0_915, lsl0_917, \
                         lsk_694, lsk_730, lsk_732, lsl1_915, lsl1_917, \
                         msk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pa_y[k] * lsl0_915[k]
                    + f_19 * lsk_730[k]
                    - f_14 * pc_y[k] * lsl1_915[k];

        t_1186[k] = f_19 * lsk_694[k]
                    + f_3 * pc_z[k] * msk_946[k];

        t_1187[k] = pa_y[k] * lsl0_917[k]
                    + f_17 * lsk_732[k]
                    - f_14 * pc_y[k] * lsl1_917[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pa_y, pc_y, lsl0_918, lsl0_920, \
                         lsl0_921, lsk_733, lsk_734, lsk_735, lsl1_918, lsl1_920, lsl1_921, \
                         msk_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pa_y[k] * lsl0_918[k]
                    + f_16 * lsk_733[k]
                    - f_14 * pc_y[k] * lsl1_918[k];

        t_1189[k] = f_15 * lsk_734[k]
                    + f_3 * pc_y[k] * msk_950[k];

        t_1190[k] = pa_y[k] * lsl0_920[k]
                    - f_14 * pc_y[k] * lsl1_920[k];

        t_1191[k] = pa_y[k] * lsl0_921[k]
                    + f_20 * lsk_735[k]
                    - f_14 * pc_y[k] * lsl1_921[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, pa_y, pc_y, pc_z, lsl0_923, lsl0_924, \
                         lsk_699, lsk_737, lsk_738, lsl1_923, lsl1_924, \
                         msk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_19 * lsk_699[k]
                    + f_3 * pc_z[k] * msk_951[k];

        t_1193[k] = pa_y[k] * lsl0_923[k]
                    + f_18 * lsk_737[k]
                    - f_14 * pc_y[k] * lsl1_923[k];

        t_1194[k] = pa_y[k] * lsl0_924[k]
                    + f_17 * lsk_738[k]
                    - f_14 * pc_y[k] * lsl1_924[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pa_y, pc_x, pc_y, lsl0_925, lsl0_927, \
                         lsk_739, lsk_740, lsk_964, lsl1_925, lsl1_927, msk_956, \
                         msk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = pa_y[k] * lsl0_925[k]
                    + f_16 * lsk_739[k]
                    - f_14 * pc_y[k] * lsl1_925[k];

        t_1196[k] = f_15 * lsk_740[k]
                    + f_3 * pc_y[k] * msk_956[k];

        t_1197[k] = pa_y[k] * lsl0_927[k]
                    - f_14 * pc_y[k] * lsl1_927[k];

        t_1198[k] = f_17 * lsk_964[k]
                    + f_3 * pc_x[k] * msk_964[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, t_1202, t_1203, pc_x, lsk_965, lsk_966, \
                         lsk_967, lsk_968, lsk_969, msk_965, msk_966, msk_967, msk_968, \
                         msk_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * lsk_965[k]
                    + f_3 * pc_x[k] * msk_965[k];

        t_1200[k] = f_17 * lsk_966[k]
                    + f_3 * pc_x[k] * msk_966[k];

        t_1201[k] = f_17 * lsk_967[k]
                    + f_3 * pc_x[k] * msk_967[k];

        t_1202[k] = f_17 * lsk_968[k]
                    + f_3 * pc_x[k] * msk_968[k];

        t_1203[k] = f_17 * lsk_969[k]
                    + f_3 * pc_x[k] * msk_969[k];
    }

#pragma omp simd aligned(t_1204, t_1205, t_1206, t_1207, pc_x, pc_y, pc_z, lsk_712, lsk_748, \
                         lsk_970, lsk_971, msi0_749, msi1_749, msk_964, msk_970, \
                         msk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1204[k] = f_17 * lsk_970[k]
                    + f_3 * pc_x[k] * msk_970[k];

        t_1205[k] = f_17 * lsk_971[k]
                    + f_3 * pc_x[k] * msk_971[k];

        t_1206[k] = f_15 * lsk_748[k]
                    + f_1 * msi0_749[k]
                    - f_2 * msi1_749[k]
                    + f_3 * pc_y[k] * msk_964[k];

        t_1207[k] = f_19 * lsk_712[k]
                    + f_3 * pc_z[k] * msk_964[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, pc_y, lsk_750, lsk_751, lsk_752, msi0_751, \
                         msi0_752, msi0_753, msi1_751, msi1_752, msi1_753, msk_966, msk_967, \
                         msk_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_15 * lsk_750[k]
                    + f_12 * msi0_751[k]
                    - f_13 * msi1_751[k]
                    + f_3 * pc_y[k] * msk_966[k];

        t_1209[k] = f_15 * lsk_751[k]
                    + f_10 * msi0_752[k]
                    - f_11 * msi1_752[k]
                    + f_3 * pc_y[k] * msk_967[k];

        t_1210[k] = f_15 * lsk_752[k]
                    + f_8 * msi0_753[k]
                    - f_9 * msi1_753[k]
                    + f_3 * pc_y[k] * msk_968[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, pc_y, lsk_753, lsk_754, lsk_755, msi0_754, \
                         msi0_755, msi1_754, msi1_755, msk_969, msk_970, \
                         msk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_15 * lsk_753[k]
                    + f_6 * msi0_754[k]
                    - f_7 * msi1_754[k]
                    + f_3 * pc_y[k] * msk_969[k];

        t_1212[k] = f_15 * lsk_754[k]
                    + f_4 * msi0_755[k]
                    - f_5 * msi1_755[k]
                    + f_3 * pc_y[k] * msk_970[k];

        t_1213[k] = f_15 * lsk_755[k]
                    + f_3 * pc_y[k] * msk_971[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pa_y, pc_x, pc_y, pc_z, lsl0_944, \
                         lsk_720, lsk_972, lsl1_944, msi0_756, msi1_756, \
                         msk_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_y[k] * lsl0_944[k]
                    - f_14 * pc_y[k] * lsl1_944[k];

        t_1215[k] = f_17 * lsk_972[k]
                    + f_1 * msi0_756[k]
                    - f_2 * msi1_756[k]
                    + f_3 * pc_x[k] * msk_972[k];

        t_1216[k] = f_3 * pc_y[k] * msk_972[k];

        t_1217[k] = f_20 * lsk_720[k]
                    + f_3 * pc_z[k] * msk_972[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pc_x, pc_y, lsk_977, msi0_756, msi0_761, \
                         msi1_756, msi1_761, msk_973, msk_974, \
                         msk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_4 * msi0_756[k]
                    - f_5 * msi1_756[k]
                    + f_3 * pc_y[k] * msk_973[k];

        t_1219[k] = f_3 * pc_y[k] * msk_974[k];

        t_1220[k] = f_17 * lsk_977[k]
                    + f_12 * msi0_761[k]
                    - f_13 * msi1_761[k]
                    + f_3 * pc_x[k] * msk_977[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pc_y, msi0_757, msi0_758, msi1_757, msi1_758, \
                         msk_975, msk_976, msk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_6 * msi0_757[k]
                    - f_7 * msi1_757[k]
                    + f_3 * pc_y[k] * msk_975[k];

        t_1222[k] = f_4 * msi0_758[k]
                    - f_5 * msi1_758[k]
                    + f_3 * pc_y[k] * msk_976[k];

        t_1223[k] = f_3 * pc_y[k] * msk_977[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pc_x, pc_y, lsk_981, msi0_759, msi0_760, \
                         msi0_765, msi1_759, msi1_760, msi1_765, msk_978, msk_979, \
                         msk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_17 * lsk_981[k]
                    + f_10 * msi0_765[k]
                    - f_11 * msi1_765[k]
                    + f_3 * pc_x[k] * msk_981[k];

        t_1225[k] = f_8 * msi0_759[k]
                    - f_9 * msi1_759[k]
                    + f_3 * pc_y[k] * msk_978[k];

        t_1226[k] = f_6 * msi0_760[k]
                    - f_7 * msi1_760[k]
                    + f_3 * pc_y[k] * msk_979[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pc_x, pc_y, lsk_986, msi0_761, msi0_770, \
                         msi1_761, msi1_770, msk_980, msk_981, \
                         msk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_4 * msi0_761[k]
                    - f_5 * msi1_761[k]
                    + f_3 * pc_y[k] * msk_980[k];

        t_1228[k] = f_3 * pc_y[k] * msk_981[k];

        t_1229[k] = f_17 * lsk_986[k]
                    + f_8 * msi0_770[k]
                    - f_9 * msi1_770[k]
                    + f_3 * pc_x[k] * msk_986[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pc_y, msi0_762, msi0_763, msi0_764, msi1_762, \
                         msi1_763, msi1_764, msk_982, msk_983, \
                         msk_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_10 * msi0_762[k]
                    - f_11 * msi1_762[k]
                    + f_3 * pc_y[k] * msk_982[k];

        t_1231[k] = f_8 * msi0_763[k]
                    - f_9 * msi1_763[k]
                    + f_3 * pc_y[k] * msk_983[k];

        t_1232[k] = f_6 * msi0_764[k]
                    - f_7 * msi1_764[k]
                    + f_3 * pc_y[k] * msk_984[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, pc_x, pc_y, lsk_992, msi0_765, msi0_776, \
                         msi1_765, msi1_776, msk_985, msk_986, \
                         msk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_4 * msi0_765[k]
                    - f_5 * msi1_765[k]
                    + f_3 * pc_y[k] * msk_985[k];

        t_1234[k] = f_3 * pc_y[k] * msk_986[k];

        t_1235[k] = f_17 * lsk_992[k]
                    + f_6 * msi0_776[k]
                    - f_7 * msi1_776[k]
                    + f_3 * pc_x[k] * msk_992[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pc_y, msi0_766, msi0_767, msi0_768, msi1_766, \
                         msi1_767, msi1_768, msk_987, msk_988, \
                         msk_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = f_12 * msi0_766[k]
                    - f_13 * msi1_766[k]
                    + f_3 * pc_y[k] * msk_987[k];

        t_1237[k] = f_10 * msi0_767[k]
                    - f_11 * msi1_767[k]
                    + f_3 * pc_y[k] * msk_988[k];

        t_1238[k] = f_8 * msi0_768[k]
                    - f_9 * msi1_768[k]
                    + f_3 * pc_y[k] * msk_989[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pc_y, msi0_769, msi0_770, msi1_769, msi1_770, \
                         msk_990, msk_991, msk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_6 * msi0_769[k]
                    - f_7 * msi1_769[k]
                    + f_3 * pc_y[k] * msk_990[k];

        t_1240[k] = f_4 * msi0_770[k]
                    - f_5 * msi1_770[k]
                    + f_3 * pc_y[k] * msk_991[k];

        t_1241[k] = f_3 * pc_y[k] * msk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pc_x, lsk_999, lsk_1000, lsk_1001, \
                         lsk_1002, msi0_783, msi1_783, msk_999, msk_1000, msk_1001, \
                         msk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_17 * lsk_999[k]
                    + f_4 * msi0_783[k]
                    - f_5 * msi1_783[k]
                    + f_3 * pc_x[k] * msk_999[k];

        t_1243[k] = f_17 * lsk_1000[k]
                    + f_3 * pc_x[k] * msk_1000[k];

        t_1244[k] = f_17 * lsk_1001[k]
                    + f_3 * pc_x[k] * msk_1001[k];

        t_1245[k] = f_17 * lsk_1002[k]
                    + f_3 * pc_x[k] * msk_1002[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, pc_x, pc_y, lsk_1003, \
                         lsk_1004, lsk_1005, lsk_1007, msk_999, msk_1003, msk_1004, msk_1005, \
                         msk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_17 * lsk_1003[k]
                    + f_3 * pc_x[k] * msk_1003[k];

        t_1247[k] = f_17 * lsk_1004[k]
                    + f_3 * pc_x[k] * msk_1004[k];

        t_1248[k] = f_17 * lsk_1005[k]
                    + f_3 * pc_x[k] * msk_1005[k];

        t_1249[k] = f_3 * pc_y[k] * msk_999[k];

        t_1250[k] = f_17 * lsk_1007[k]
                    + f_3 * pc_x[k] * msk_1007[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
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
    const auto f_24 = 3.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_945 = buffer.data(lsl0 + 945);
    const auto *lsl0_948 = buffer.data(lsl0 + 948);
    const auto *lsl0_951 = buffer.data(lsl0 + 951);
    const auto *lsl0_955 = buffer.data(lsl0 + 955);
    const auto *lsl0_957 = buffer.data(lsl0 + 957);
    const auto *lsl0_960 = buffer.data(lsl0 + 960);
    const auto *lsl0_962 = buffer.data(lsl0 + 962);
    const auto *lsl0_963 = buffer.data(lsl0 + 963);
    const auto *lsl0_966 = buffer.data(lsl0 + 966);
    const auto *lsl0_968 = buffer.data(lsl0 + 968);
    const auto *lsl0_969 = buffer.data(lsl0 + 969);
    const auto *lsl0_970 = buffer.data(lsl0 + 970);
    const auto *lsl0_981 = buffer.data(lsl0 + 981);

    const auto *lsk_755 = buffer.data(lsk + 755);
    const auto *lsk_756 = buffer.data(lsk + 756);
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
    const auto *lsk_784 = buffer.data(lsk + 784);
    const auto *lsk_791 = buffer.data(lsk + 791);
    const auto *lsk_792 = buffer.data(lsk + 792);
    const auto *lsk_794 = buffer.data(lsk + 794);
    const auto *lsk_795 = buffer.data(lsk + 795);
    const auto *lsk_797 = buffer.data(lsk + 797);
    const auto *lsk_798 = buffer.data(lsk + 798);
    const auto *lsk_801 = buffer.data(lsk + 801);
    const auto *lsk_806 = buffer.data(lsk + 806);
    const auto *lsk_812 = buffer.data(lsk + 812);
    const auto *lsk_822 = buffer.data(lsk + 822);
    const auto *lsk_823 = buffer.data(lsk + 823);
    const auto *lsk_824 = buffer.data(lsk + 824);
    const auto *lsk_825 = buffer.data(lsk + 825);
    const auto *lsk_826 = buffer.data(lsk + 826);
    const auto *lsk_827 = buffer.data(lsk + 827);
    const auto *lsk_828 = buffer.data(lsk + 828);
    const auto *lsk_830 = buffer.data(lsk + 830);
    const auto *lsk_833 = buffer.data(lsk + 833);
    const auto *lsk_837 = buffer.data(lsk + 837);
    const auto *lsk_1008 = buffer.data(lsk + 1008);
    const auto *lsk_1011 = buffer.data(lsk + 1011);
    const auto *lsk_1014 = buffer.data(lsk + 1014);
    const auto *lsk_1018 = buffer.data(lsk + 1018);
    const auto *lsk_1023 = buffer.data(lsk + 1023);
    const auto *lsk_1029 = buffer.data(lsk + 1029);
    const auto *lsk_1036 = buffer.data(lsk + 1036);
    const auto *lsk_1038 = buffer.data(lsk + 1038);
    const auto *lsk_1039 = buffer.data(lsk + 1039);
    const auto *lsk_1040 = buffer.data(lsk + 1040);
    const auto *lsk_1041 = buffer.data(lsk + 1041);
    const auto *lsk_1042 = buffer.data(lsk + 1042);
    const auto *lsk_1043 = buffer.data(lsk + 1043);
    const auto *lsk_1049 = buffer.data(lsk + 1049);
    const auto *lsk_1053 = buffer.data(lsk + 1053);
    const auto *lsk_1058 = buffer.data(lsk + 1058);
    const auto *lsk_1064 = buffer.data(lsk + 1064);
    const auto *lsk_1071 = buffer.data(lsk + 1071);
    const auto *lsk_1072 = buffer.data(lsk + 1072);
    const auto *lsk_1073 = buffer.data(lsk + 1073);
    const auto *lsk_1074 = buffer.data(lsk + 1074);
    const auto *lsk_1075 = buffer.data(lsk + 1075);
    const auto *lsk_1076 = buffer.data(lsk + 1076);
    const auto *lsk_1077 = buffer.data(lsk + 1077);
    const auto *lsk_1078 = buffer.data(lsk + 1078);
    const auto *lsk_1079 = buffer.data(lsk + 1079);
    const auto *lsk_1080 = buffer.data(lsk + 1080);
    const auto *lsk_1083 = buffer.data(lsk + 1083);
    const auto *lsk_1085 = buffer.data(lsk + 1085);
    const auto *lsk_1086 = buffer.data(lsk + 1086);
    const auto *lsk_1089 = buffer.data(lsk + 1089);
    const auto *lsk_1090 = buffer.data(lsk + 1090);
    const auto *lsk_1092 = buffer.data(lsk + 1092);

    const auto *lsl1_945 = buffer.data(lsl1 + 945);
    const auto *lsl1_948 = buffer.data(lsl1 + 948);
    const auto *lsl1_951 = buffer.data(lsl1 + 951);
    const auto *lsl1_955 = buffer.data(lsl1 + 955);
    const auto *lsl1_957 = buffer.data(lsl1 + 957);
    const auto *lsl1_960 = buffer.data(lsl1 + 960);
    const auto *lsl1_962 = buffer.data(lsl1 + 962);
    const auto *lsl1_963 = buffer.data(lsl1 + 963);
    const auto *lsl1_966 = buffer.data(lsl1 + 966);
    const auto *lsl1_968 = buffer.data(lsl1 + 968);
    const auto *lsl1_969 = buffer.data(lsl1 + 969);
    const auto *lsl1_970 = buffer.data(lsl1 + 970);
    const auto *lsl1_981 = buffer.data(lsl1 + 981);

    const auto *msi0_777 = buffer.data(msi0 + 777);
    const auto *msi0_778 = buffer.data(msi0 + 778);
    const auto *msi0_779 = buffer.data(msi0 + 779);
    const auto *msi0_780 = buffer.data(msi0 + 780);
    const auto *msi0_781 = buffer.data(msi0 + 781);
    const auto *msi0_782 = buffer.data(msi0 + 782);
    const auto *msi0_783 = buffer.data(msi0 + 783);
    const auto *msi0_784 = buffer.data(msi0 + 784);
    const auto *msi0_786 = buffer.data(msi0 + 786);
    const auto *msi0_787 = buffer.data(msi0 + 787);
    const auto *msi0_789 = buffer.data(msi0 + 789);
    const auto *msi0_790 = buffer.data(msi0 + 790);
    const auto *msi0_791 = buffer.data(msi0 + 791);
    const auto *msi0_793 = buffer.data(msi0 + 793);
    const auto *msi0_794 = buffer.data(msi0 + 794);
    const auto *msi0_795 = buffer.data(msi0 + 795);
    const auto *msi0_796 = buffer.data(msi0 + 796);
    const auto *msi0_798 = buffer.data(msi0 + 798);
    const auto *msi0_799 = buffer.data(msi0 + 799);
    const auto *msi0_805 = buffer.data(msi0 + 805);
    const auto *msi0_806 = buffer.data(msi0 + 806);
    const auto *msi0_807 = buffer.data(msi0 + 807);
    const auto *msi0_808 = buffer.data(msi0 + 808);
    const auto *msi0_809 = buffer.data(msi0 + 809);
    const auto *msi0_811 = buffer.data(msi0 + 811);
    const auto *msi0_817 = buffer.data(msi0 + 817);
    const auto *msi0_821 = buffer.data(msi0 + 821);
    const auto *msi0_826 = buffer.data(msi0 + 826);
    const auto *msi0_832 = buffer.data(msi0 + 832);
    const auto *msi0_835 = buffer.data(msi0 + 835);
    const auto *msi0_836 = buffer.data(msi0 + 836);
    const auto *msi0_837 = buffer.data(msi0 + 837);
    const auto *msi0_838 = buffer.data(msi0 + 838);
    const auto *msi0_839 = buffer.data(msi0 + 839);
    const auto *msi0_840 = buffer.data(msi0 + 840);
    const auto *msi0_843 = buffer.data(msi0 + 843);
    const auto *msi0_845 = buffer.data(msi0 + 845);
    const auto *msi0_846 = buffer.data(msi0 + 846);
    const auto *msi0_849 = buffer.data(msi0 + 849);
    const auto *msi0_850 = buffer.data(msi0 + 850);
    const auto *msi0_852 = buffer.data(msi0 + 852);

    const auto *msi1_777 = buffer.data(msi1 + 777);
    const auto *msi1_778 = buffer.data(msi1 + 778);
    const auto *msi1_779 = buffer.data(msi1 + 779);
    const auto *msi1_780 = buffer.data(msi1 + 780);
    const auto *msi1_781 = buffer.data(msi1 + 781);
    const auto *msi1_782 = buffer.data(msi1 + 782);
    const auto *msi1_783 = buffer.data(msi1 + 783);
    const auto *msi1_784 = buffer.data(msi1 + 784);
    const auto *msi1_786 = buffer.data(msi1 + 786);
    const auto *msi1_787 = buffer.data(msi1 + 787);
    const auto *msi1_789 = buffer.data(msi1 + 789);
    const auto *msi1_790 = buffer.data(msi1 + 790);
    const auto *msi1_791 = buffer.data(msi1 + 791);
    const auto *msi1_793 = buffer.data(msi1 + 793);
    const auto *msi1_794 = buffer.data(msi1 + 794);
    const auto *msi1_795 = buffer.data(msi1 + 795);
    const auto *msi1_796 = buffer.data(msi1 + 796);
    const auto *msi1_798 = buffer.data(msi1 + 798);
    const auto *msi1_799 = buffer.data(msi1 + 799);
    const auto *msi1_805 = buffer.data(msi1 + 805);
    const auto *msi1_806 = buffer.data(msi1 + 806);
    const auto *msi1_807 = buffer.data(msi1 + 807);
    const auto *msi1_808 = buffer.data(msi1 + 808);
    const auto *msi1_809 = buffer.data(msi1 + 809);
    const auto *msi1_811 = buffer.data(msi1 + 811);
    const auto *msi1_817 = buffer.data(msi1 + 817);
    const auto *msi1_821 = buffer.data(msi1 + 821);
    const auto *msi1_826 = buffer.data(msi1 + 826);
    const auto *msi1_832 = buffer.data(msi1 + 832);
    const auto *msi1_835 = buffer.data(msi1 + 835);
    const auto *msi1_836 = buffer.data(msi1 + 836);
    const auto *msi1_837 = buffer.data(msi1 + 837);
    const auto *msi1_838 = buffer.data(msi1 + 838);
    const auto *msi1_839 = buffer.data(msi1 + 839);
    const auto *msi1_840 = buffer.data(msi1 + 840);
    const auto *msi1_843 = buffer.data(msi1 + 843);
    const auto *msi1_845 = buffer.data(msi1 + 845);
    const auto *msi1_846 = buffer.data(msi1 + 846);
    const auto *msi1_849 = buffer.data(msi1 + 849);
    const auto *msi1_850 = buffer.data(msi1 + 850);
    const auto *msi1_852 = buffer.data(msi1 + 852);

    const auto *msk_1000 = buffer.data(msk + 1000);
    const auto *msk_1001 = buffer.data(msk + 1001);
    const auto *msk_1002 = buffer.data(msk + 1002);
    const auto *msk_1003 = buffer.data(msk + 1003);
    const auto *msk_1004 = buffer.data(msk + 1004);
    const auto *msk_1005 = buffer.data(msk + 1005);
    const auto *msk_1006 = buffer.data(msk + 1006);
    const auto *msk_1007 = buffer.data(msk + 1007);
    const auto *msk_1008 = buffer.data(msk + 1008);
    const auto *msk_1009 = buffer.data(msk + 1009);
    const auto *msk_1010 = buffer.data(msk + 1010);
    const auto *msk_1011 = buffer.data(msk + 1011);
    const auto *msk_1013 = buffer.data(msk + 1013);
    const auto *msk_1014 = buffer.data(msk + 1014);
    const auto *msk_1015 = buffer.data(msk + 1015);
    const auto *msk_1017 = buffer.data(msk + 1017);
    const auto *msk_1018 = buffer.data(msk + 1018);
    const auto *msk_1019 = buffer.data(msk + 1019);
    const auto *msk_1020 = buffer.data(msk + 1020);
    const auto *msk_1022 = buffer.data(msk + 1022);
    const auto *msk_1023 = buffer.data(msk + 1023);
    const auto *msk_1024 = buffer.data(msk + 1024);
    const auto *msk_1025 = buffer.data(msk + 1025);
    const auto *msk_1026 = buffer.data(msk + 1026);
    const auto *msk_1028 = buffer.data(msk + 1028);
    const auto *msk_1029 = buffer.data(msk + 1029);
    const auto *msk_1036 = buffer.data(msk + 1036);
    const auto *msk_1037 = buffer.data(msk + 1037);
    const auto *msk_1038 = buffer.data(msk + 1038);
    const auto *msk_1039 = buffer.data(msk + 1039);
    const auto *msk_1040 = buffer.data(msk + 1040);
    const auto *msk_1041 = buffer.data(msk + 1041);
    const auto *msk_1042 = buffer.data(msk + 1042);
    const auto *msk_1043 = buffer.data(msk + 1043);
    const auto *msk_1044 = buffer.data(msk + 1044);
    const auto *msk_1046 = buffer.data(msk + 1046);
    const auto *msk_1047 = buffer.data(msk + 1047);
    const auto *msk_1049 = buffer.data(msk + 1049);
    const auto *msk_1050 = buffer.data(msk + 1050);
    const auto *msk_1053 = buffer.data(msk + 1053);
    const auto *msk_1054 = buffer.data(msk + 1054);
    const auto *msk_1058 = buffer.data(msk + 1058);
    const auto *msk_1059 = buffer.data(msk + 1059);
    const auto *msk_1064 = buffer.data(msk + 1064);
    const auto *msk_1071 = buffer.data(msk + 1071);
    const auto *msk_1072 = buffer.data(msk + 1072);
    const auto *msk_1073 = buffer.data(msk + 1073);
    const auto *msk_1074 = buffer.data(msk + 1074);
    const auto *msk_1075 = buffer.data(msk + 1075);
    const auto *msk_1076 = buffer.data(msk + 1076);
    const auto *msk_1077 = buffer.data(msk + 1077);
    const auto *msk_1078 = buffer.data(msk + 1078);
    const auto *msk_1079 = buffer.data(msk + 1079);
    const auto *msk_1080 = buffer.data(msk + 1080);
    const auto *msk_1082 = buffer.data(msk + 1082);
    const auto *msk_1083 = buffer.data(msk + 1083);
    const auto *msk_1085 = buffer.data(msk + 1085);
    const auto *msk_1086 = buffer.data(msk + 1086);
    const auto *msk_1089 = buffer.data(msk + 1089);
    const auto *msk_1090 = buffer.data(msk + 1090);
    const auto *msk_1092 = buffer.data(msk + 1092);

#pragma omp simd aligned(t_1251, t_1252, t_1253, pc_y, msi0_777, msi0_778, msi0_779, msi1_777, \
                         msi1_778, msi1_779, msk_1000, msk_1001, \
                         msk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_1 * msi0_777[k]
                    - f_2 * msi1_777[k]
                    + f_3 * pc_y[k] * msk_1000[k];

        t_1252[k] = f_22 * msi0_778[k]
                    - f_23 * msi1_778[k]
                    + f_3 * pc_y[k] * msk_1001[k];

        t_1253[k] = f_12 * msi0_779[k]
                    - f_13 * msi1_779[k]
                    + f_3 * pc_y[k] * msk_1002[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, msi0_780, msi0_781, msi0_782, msi1_780, \
                         msi1_781, msi1_782, msk_1003, msk_1004, \
                         msk_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_10 * msi0_780[k]
                    - f_11 * msi1_780[k]
                    + f_3 * pc_y[k] * msk_1003[k];

        t_1255[k] = f_8 * msi0_781[k]
                    - f_9 * msi1_781[k]
                    + f_3 * pc_y[k] * msk_1004[k];

        t_1256[k] = f_6 * msi0_782[k]
                    - f_7 * msi1_782[k]
                    + f_3 * pc_y[k] * msk_1005[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, t_1260, pc_x, pc_y, pc_z, lsk_755, lsk_1008, \
                         msi0_783, msi0_784, msi1_783, msi1_784, msk_1006, msk_1007, \
                         msk_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_4 * msi0_783[k]
                    - f_5 * msi1_783[k]
                    + f_3 * pc_y[k] * msk_1006[k];

        t_1258[k] = f_3 * pc_y[k] * msk_1007[k];

        t_1259[k] = f_20 * lsk_755[k]
                    + f_1 * msi0_783[k]
                    - f_2 * msi1_783[k]
                    + f_3 * pc_z[k] * msk_1007[k];

        t_1260[k] = f_16 * lsk_1008[k]
                    + f_1 * msi0_784[k]
                    - f_2 * msi1_784[k]
                    + f_3 * pc_x[k] * msk_1008[k];
    }

#pragma omp simd aligned(t_1261, t_1262, t_1263, t_1264, pc_x, pc_y, pc_z, lsk_756, lsk_1011, \
                         msi0_787, msi1_787, msk_1008, msk_1009, \
                         msk_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1261[k] = f_24 * lsk_756[k]
                    + f_3 * pc_y[k] * msk_1008[k];

        t_1262[k] = f_3 * pc_z[k] * msk_1008[k];

        t_1263[k] = f_16 * lsk_1011[k]
                    + f_12 * msi0_787[k]
                    - f_13 * msi1_787[k]
                    + f_3 * pc_x[k] * msk_1011[k];

        t_1264[k] = f_3 * pc_z[k] * msk_1009[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, pc_x, pc_z, lsk_1014, msi0_784, msi0_790, \
                         msi1_784, msi1_790, msk_1010, msk_1011, \
                         msk_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_4 * msi0_784[k]
                    - f_5 * msi1_784[k]
                    + f_3 * pc_z[k] * msk_1010[k];

        t_1266[k] = f_16 * lsk_1014[k]
                    + f_10 * msi0_790[k]
                    - f_11 * msi1_790[k]
                    + f_3 * pc_x[k] * msk_1014[k];

        t_1267[k] = f_3 * pc_z[k] * msk_1011[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pc_x, pc_y, pc_z, lsk_761, lsk_1018, \
                         msi0_786, msi0_794, msi1_786, msi1_794, msk_1013, msk_1014, \
                         msk_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_24 * lsk_761[k]
                    + f_3 * pc_y[k] * msk_1013[k];

        t_1269[k] = f_6 * msi0_786[k]
                    - f_7 * msi1_786[k]
                    + f_3 * pc_z[k] * msk_1013[k];

        t_1270[k] = f_16 * lsk_1018[k]
                    + f_8 * msi0_794[k]
                    - f_9 * msi1_794[k]
                    + f_3 * pc_x[k] * msk_1018[k];

        t_1271[k] = f_3 * pc_z[k] * msk_1014[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pc_y, pc_z, lsk_765, msi0_787, msi0_789, \
                         msi1_787, msi1_789, msk_1015, msk_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * msi0_787[k]
                    - f_5 * msi1_787[k]
                    + f_3 * pc_z[k] * msk_1015[k];

        t_1273[k] = f_24 * lsk_765[k]
                    + f_3 * pc_y[k] * msk_1017[k];

        t_1274[k] = f_8 * msi0_789[k]
                    - f_9 * msi1_789[k]
                    + f_3 * pc_z[k] * msk_1017[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pc_x, pc_z, lsk_1023, msi0_790, msi0_799, \
                         msi1_790, msi1_799, msk_1018, msk_1019, \
                         msk_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_16 * lsk_1023[k]
                    + f_6 * msi0_799[k]
                    - f_7 * msi1_799[k]
                    + f_3 * pc_x[k] * msk_1023[k];

        t_1276[k] = f_3 * pc_z[k] * msk_1018[k];

        t_1277[k] = f_4 * msi0_790[k]
                    - f_5 * msi1_790[k]
                    + f_3 * pc_z[k] * msk_1019[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, pc_y, pc_z, lsk_770, msi0_791, msi0_793, \
                         msi1_791, msi1_793, msk_1020, msk_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_6 * msi0_791[k]
                    - f_7 * msi1_791[k]
                    + f_3 * pc_z[k] * msk_1020[k];

        t_1279[k] = f_24 * lsk_770[k]
                    + f_3 * pc_y[k] * msk_1022[k];

        t_1280[k] = f_10 * msi0_793[k]
                    - f_11 * msi1_793[k]
                    + f_3 * pc_z[k] * msk_1022[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pc_x, pc_z, lsk_1029, msi0_794, msi0_805, \
                         msi1_794, msi1_805, msk_1023, msk_1024, \
                         msk_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_16 * lsk_1029[k]
                    + f_4 * msi0_805[k]
                    - f_5 * msi1_805[k]
                    + f_3 * pc_x[k] * msk_1029[k];

        t_1282[k] = f_3 * pc_z[k] * msk_1023[k];

        t_1283[k] = f_4 * msi0_794[k]
                    - f_5 * msi1_794[k]
                    + f_3 * pc_z[k] * msk_1024[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_y, pc_z, lsk_776, msi0_795, \
                         msi0_796, msi0_798, msi1_795, msi1_796, msi1_798, msk_1025, msk_1026, \
                         msk_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_6 * msi0_795[k]
                    - f_7 * msi1_795[k]
                    + f_3 * pc_z[k] * msk_1025[k];

        t_1285[k] = f_8 * msi0_796[k]
                    - f_9 * msi1_796[k]
                    + f_3 * pc_z[k] * msk_1026[k];

        t_1286[k] = f_24 * lsk_776[k]
                    + f_3 * pc_y[k] * msk_1028[k];

        t_1287[k] = f_12 * msi0_798[k]
                    - f_13 * msi1_798[k]
                    + f_3 * pc_z[k] * msk_1028[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, t_1292, pc_x, pc_z, lsk_1036, \
                         lsk_1038, lsk_1039, lsk_1040, msk_1029, msk_1036, msk_1038, msk_1039, \
                         msk_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = f_16 * lsk_1036[k]
                    + f_3 * pc_x[k] * msk_1036[k];

        t_1289[k] = f_3 * pc_z[k] * msk_1029[k];

        t_1290[k] = f_16 * lsk_1038[k]
                    + f_3 * pc_x[k] * msk_1038[k];

        t_1291[k] = f_16 * lsk_1039[k]
                    + f_3 * pc_x[k] * msk_1039[k];

        t_1292[k] = f_16 * lsk_1040[k]
                    + f_3 * pc_x[k] * msk_1040[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pc_x, pc_y, lsk_784, lsk_1041, \
                         lsk_1042, lsk_1043, msi0_805, msi1_805, msk_1036, msk_1041, msk_1042, \
                         msk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_16 * lsk_1041[k]
                    + f_3 * pc_x[k] * msk_1041[k];

        t_1294[k] = f_16 * lsk_1042[k]
                    + f_3 * pc_x[k] * msk_1042[k];

        t_1295[k] = f_16 * lsk_1043[k]
                    + f_3 * pc_x[k] * msk_1043[k];

        t_1296[k] = f_24 * lsk_784[k]
                    + f_1 * msi0_805[k]
                    - f_2 * msi1_805[k]
                    + f_3 * pc_y[k] * msk_1036[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, t_1300, pc_z, msi0_805, msi0_806, msi0_807, \
                         msi1_805, msi1_806, msi1_807, msk_1036, msk_1037, msk_1038, \
                         msk_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = f_3 * pc_z[k] * msk_1036[k];

        t_1298[k] = f_4 * msi0_805[k]
                    - f_5 * msi1_805[k]
                    + f_3 * pc_z[k] * msk_1037[k];

        t_1299[k] = f_6 * msi0_806[k]
                    - f_7 * msi1_806[k]
                    + f_3 * pc_z[k] * msk_1038[k];

        t_1300[k] = f_8 * msi0_807[k]
                    - f_9 * msi1_807[k]
                    + f_3 * pc_z[k] * msk_1039[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pc_y, pc_z, lsk_791, msi0_808, \
                         msi0_809, msi0_811, msi1_808, msi1_809, msi1_811, msk_1040, msk_1041, \
                         msk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_10 * msi0_808[k]
                    - f_11 * msi1_808[k]
                    + f_3 * pc_z[k] * msk_1040[k];

        t_1302[k] = f_12 * msi0_809[k]
                    - f_13 * msi1_809[k]
                    + f_3 * pc_z[k] * msk_1041[k];

        t_1303[k] = f_24 * lsk_791[k]
                    + f_3 * pc_y[k] * msk_1043[k];

        t_1304[k] = f_1 * msi0_811[k]
                    - f_2 * msi1_811[k]
                    + f_3 * pc_z[k] * msk_1043[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, pa_z, pc_y, pc_z, lsl0_945, lsl0_948, \
                         lsk_756, lsk_792, lsl1_945, lsl1_948, \
                         msk_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = pa_z[k] * lsl0_945[k]
                    - f_14 * pc_z[k] * lsl1_945[k];

        t_1306[k] = f_20 * lsk_792[k]
                    + f_3 * pc_y[k] * msk_1044[k];

        t_1307[k] = f_15 * lsk_756[k]
                    + f_3 * pc_z[k] * msk_1044[k];

        t_1308[k] = pa_z[k] * lsl0_948[k]
                    - f_14 * pc_z[k] * lsl1_948[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, pa_z, pc_x, pc_y, pc_z, lsl0_951, lsk_794, \
                         lsk_1049, lsl1_951, msi0_817, msi1_817, msk_1046, \
                         msk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_20 * lsk_794[k]
                    + f_3 * pc_y[k] * msk_1046[k];

        t_1310[k] = f_16 * lsk_1049[k]
                    + f_12 * msi0_817[k]
                    - f_13 * msi1_817[k]
                    + f_3 * pc_x[k] * msk_1049[k];

        t_1311[k] = pa_z[k] * lsl0_951[k]
                    - f_14 * pc_z[k] * lsl1_951[k];
    }

#pragma omp simd aligned(t_1312, t_1313, t_1314, pc_x, pc_y, pc_z, lsk_759, lsk_797, lsk_1053, \
                         msi0_821, msi1_821, msk_1047, msk_1049, \
                         msk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1312[k] = f_15 * lsk_759[k]
                    + f_3 * pc_z[k] * msk_1047[k];

        t_1313[k] = f_20 * lsk_797[k]
                    + f_3 * pc_y[k] * msk_1049[k];

        t_1314[k] = f_16 * lsk_1053[k]
                    + f_10 * msi0_821[k]
                    - f_11 * msi1_821[k]
                    + f_3 * pc_x[k] * msk_1053[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, pa_z, pc_y, pc_z, lsl0_955, lsl0_957, \
                         lsk_762, lsk_763, lsk_801, lsl1_955, lsl1_957, msk_1050, \
                         msk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = pa_z[k] * lsl0_955[k]
                    - f_14 * pc_z[k] * lsl1_955[k];

        t_1316[k] = f_15 * lsk_762[k]
                    + f_3 * pc_z[k] * msk_1050[k];

        t_1317[k] = pa_z[k] * lsl0_957[k]
                    + f_16 * lsk_763[k]
                    - f_14 * pc_z[k] * lsl1_957[k];

        t_1318[k] = f_20 * lsk_801[k]
                    + f_3 * pc_y[k] * msk_1053[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pa_z, pc_x, pc_z, lsl0_960, lsk_766, \
                         lsk_1058, lsl1_960, msi0_826, msi1_826, msk_1054, \
                         msk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_16 * lsk_1058[k]
                    + f_8 * msi0_826[k]
                    - f_9 * msi1_826[k]
                    + f_3 * pc_x[k] * msk_1058[k];

        t_1320[k] = pa_z[k] * lsl0_960[k]
                    - f_14 * pc_z[k] * lsl1_960[k];

        t_1321[k] = f_15 * lsk_766[k]
                    + f_3 * pc_z[k] * msk_1054[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pa_z, pc_y, pc_z, lsl0_962, lsl0_963, \
                         lsk_767, lsk_768, lsk_806, lsl1_962, lsl1_963, \
                         msk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = pa_z[k] * lsl0_962[k]
                    + f_16 * lsk_767[k]
                    - f_14 * pc_z[k] * lsl1_962[k];

        t_1323[k] = pa_z[k] * lsl0_963[k]
                    + f_17 * lsk_768[k]
                    - f_14 * pc_z[k] * lsl1_963[k];

        t_1324[k] = f_20 * lsk_806[k]
                    + f_3 * pc_y[k] * msk_1058[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pa_z, pc_x, pc_z, lsl0_966, lsk_771, \
                         lsk_1064, lsl1_966, msi0_832, msi1_832, msk_1059, \
                         msk_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_16 * lsk_1064[k]
                    + f_6 * msi0_832[k]
                    - f_7 * msi1_832[k]
                    + f_3 * pc_x[k] * msk_1064[k];

        t_1326[k] = pa_z[k] * lsl0_966[k]
                    - f_14 * pc_z[k] * lsl1_966[k];

        t_1327[k] = f_15 * lsk_771[k]
                    + f_3 * pc_z[k] * msk_1059[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pa_z, pc_z, lsl0_968, lsl0_969, lsl0_970, \
                         lsk_772, lsk_773, lsk_774, lsl1_968, lsl1_969, \
                         lsl1_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = pa_z[k] * lsl0_968[k]
                    + f_16 * lsk_772[k]
                    - f_14 * pc_z[k] * lsl1_968[k];

        t_1329[k] = pa_z[k] * lsl0_969[k]
                    + f_17 * lsk_773[k]
                    - f_14 * pc_z[k] * lsl1_969[k];

        t_1330[k] = pa_z[k] * lsl0_970[k]
                    + f_18 * lsk_774[k]
                    - f_14 * pc_z[k] * lsl1_970[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, pc_x, pc_y, lsk_812, lsk_1071, \
                         lsk_1072, lsk_1073, msi0_839, msi1_839, msk_1064, msk_1071, msk_1072, \
                         msk_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_20 * lsk_812[k]
                    + f_3 * pc_y[k] * msk_1064[k];

        t_1332[k] = f_16 * lsk_1071[k]
                    + f_4 * msi0_839[k]
                    - f_5 * msi1_839[k]
                    + f_3 * pc_x[k] * msk_1071[k];

        t_1333[k] = f_16 * lsk_1072[k]
                    + f_3 * pc_x[k] * msk_1072[k];

        t_1334[k] = f_16 * lsk_1073[k]
                    + f_3 * pc_x[k] * msk_1073[k];
    }

#pragma omp simd aligned(t_1335, t_1336, t_1337, t_1338, t_1339, pc_x, lsk_1074, lsk_1075, \
                         lsk_1076, lsk_1077, lsk_1078, msk_1074, msk_1075, msk_1076, msk_1077, \
                         msk_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1335[k] = f_16 * lsk_1074[k]
                    + f_3 * pc_x[k] * msk_1074[k];

        t_1336[k] = f_16 * lsk_1075[k]
                    + f_3 * pc_x[k] * msk_1075[k];

        t_1337[k] = f_16 * lsk_1076[k]
                    + f_3 * pc_x[k] * msk_1076[k];

        t_1338[k] = f_16 * lsk_1077[k]
                    + f_3 * pc_x[k] * msk_1077[k];

        t_1339[k] = f_16 * lsk_1078[k]
                    + f_3 * pc_x[k] * msk_1078[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pa_z, pc_x, pc_z, lsl0_981, lsk_784, \
                         lsk_1079, lsl1_981, msk_1072, msk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_16 * lsk_1079[k]
                    + f_3 * pc_x[k] * msk_1079[k];

        t_1341[k] = pa_z[k] * lsl0_981[k]
                    - f_14 * pc_z[k] * lsl1_981[k];

        t_1342[k] = f_15 * lsk_784[k]
                    + f_3 * pc_z[k] * msk_1072[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pc_y, lsk_822, lsk_823, lsk_824, msi0_835, \
                         msi0_836, msi0_837, msi1_835, msi1_836, msi1_837, msk_1074, msk_1075, \
                         msk_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_20 * lsk_822[k]
                    + f_12 * msi0_835[k]
                    - f_13 * msi1_835[k]
                    + f_3 * pc_y[k] * msk_1074[k];

        t_1344[k] = f_20 * lsk_823[k]
                    + f_10 * msi0_836[k]
                    - f_11 * msi1_836[k]
                    + f_3 * pc_y[k] * msk_1075[k];

        t_1345[k] = f_20 * lsk_824[k]
                    + f_8 * msi0_837[k]
                    - f_9 * msi1_837[k]
                    + f_3 * pc_y[k] * msk_1076[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pc_y, lsk_825, lsk_826, lsk_827, msi0_838, \
                         msi0_839, msi1_838, msi1_839, msk_1077, msk_1078, \
                         msk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_20 * lsk_825[k]
                    + f_6 * msi0_838[k]
                    - f_7 * msi1_838[k]
                    + f_3 * pc_y[k] * msk_1077[k];

        t_1347[k] = f_20 * lsk_826[k]
                    + f_4 * msi0_839[k]
                    - f_5 * msi1_839[k]
                    + f_3 * pc_y[k] * msk_1078[k];

        t_1348[k] = f_20 * lsk_827[k]
                    + f_3 * pc_y[k] * msk_1079[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pc_x, pc_y, pc_z, lsk_791, lsk_828, lsk_1080, \
                         msi0_839, msi0_840, msi1_839, msi1_840, msk_1079, \
                         msk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_15 * lsk_791[k]
                    + f_1 * msi0_839[k]
                    - f_2 * msi1_839[k]
                    + f_3 * pc_z[k] * msk_1079[k];

        t_1350[k] = f_16 * lsk_1080[k]
                    + f_1 * msi0_840[k]
                    - f_2 * msi1_840[k]
                    + f_3 * pc_x[k] * msk_1080[k];

        t_1351[k] = f_19 * lsk_828[k]
                    + f_3 * pc_y[k] * msk_1080[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pc_x, pc_y, pc_z, lsk_792, lsk_830, lsk_1083, \
                         msi0_843, msi1_843, msk_1080, msk_1082, \
                         msk_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_16 * lsk_792[k]
                    + f_3 * pc_z[k] * msk_1080[k];

        t_1353[k] = f_16 * lsk_1083[k]
                    + f_12 * msi0_843[k]
                    - f_13 * msi1_843[k]
                    + f_3 * pc_x[k] * msk_1083[k];

        t_1354[k] = f_19 * lsk_830[k]
                    + f_3 * pc_y[k] * msk_1082[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pc_x, pc_z, lsk_795, lsk_1085, lsk_1086, \
                         msi0_845, msi0_846, msi1_845, msi1_846, msk_1083, msk_1085, \
                         msk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_16 * lsk_1085[k]
                    + f_12 * msi0_845[k]
                    - f_13 * msi1_845[k]
                    + f_3 * pc_x[k] * msk_1085[k];

        t_1356[k] = f_16 * lsk_1086[k]
                    + f_10 * msi0_846[k]
                    - f_11 * msi1_846[k]
                    + f_3 * pc_x[k] * msk_1086[k];

        t_1357[k] = f_16 * lsk_795[k]
                    + f_3 * pc_z[k] * msk_1083[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, pc_x, pc_y, lsk_833, lsk_1089, lsk_1090, \
                         msi0_849, msi0_850, msi1_849, msi1_850, msk_1085, msk_1089, \
                         msk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_19 * lsk_833[k]
                    + f_3 * pc_y[k] * msk_1085[k];

        t_1359[k] = f_16 * lsk_1089[k]
                    + f_10 * msi0_849[k]
                    - f_11 * msi1_849[k]
                    + f_3 * pc_x[k] * msk_1089[k];

        t_1360[k] = f_16 * lsk_1090[k]
                    + f_8 * msi0_850[k]
                    - f_9 * msi1_850[k]
                    + f_3 * pc_x[k] * msk_1090[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, pc_x, pc_y, pc_z, lsk_798, lsk_837, lsk_1092, \
                         msi0_852, msi1_852, msk_1086, msk_1089, \
                         msk_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = f_16 * lsk_798[k]
                    + f_3 * pc_z[k] * msk_1086[k];

        t_1362[k] = f_16 * lsk_1092[k]
                    + f_8 * msi0_852[k]
                    - f_9 * msi1_852[k]
                    + f_3 * pc_x[k] * msk_1092[k];

        t_1363[k] = f_19 * lsk_837[k]
                    + f_3 * pc_y[k] * msk_1089[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t lsk, const size_t msi0,
                                                           const size_t msi1, const size_t msk,
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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk_802 = buffer.data(lsk + 802);
    const auto *lsk_807 = buffer.data(lsk + 807);
    const auto *lsk_820 = buffer.data(lsk + 820);
    const auto *lsk_827 = buffer.data(lsk + 827);
    const auto *lsk_828 = buffer.data(lsk + 828);
    const auto *lsk_831 = buffer.data(lsk + 831);
    const auto *lsk_834 = buffer.data(lsk + 834);
    const auto *lsk_838 = buffer.data(lsk + 838);
    const auto *lsk_842 = buffer.data(lsk + 842);
    const auto *lsk_843 = buffer.data(lsk + 843);
    const auto *lsk_848 = buffer.data(lsk + 848);
    const auto *lsk_856 = buffer.data(lsk + 856);
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
    const auto *lsk_878 = buffer.data(lsk + 878);
    const auto *lsk_879 = buffer.data(lsk + 879);
    const auto *lsk_884 = buffer.data(lsk + 884);
    const auto *lsk_892 = buffer.data(lsk + 892);
    const auto *lsk_894 = buffer.data(lsk + 894);
    const auto *lsk_895 = buffer.data(lsk + 895);
    const auto *lsk_896 = buffer.data(lsk + 896);
    const auto *lsk_897 = buffer.data(lsk + 897);
    const auto *lsk_898 = buffer.data(lsk + 898);
    const auto *lsk_899 = buffer.data(lsk + 899);
    const auto *lsk_900 = buffer.data(lsk + 900);
    const auto *lsk_902 = buffer.data(lsk + 902);
    const auto *lsk_905 = buffer.data(lsk + 905);
    const auto *lsk_909 = buffer.data(lsk + 909);
    const auto *lsk_914 = buffer.data(lsk + 914);
    const auto *lsk_1094 = buffer.data(lsk + 1094);
    const auto *lsk_1095 = buffer.data(lsk + 1095);
    const auto *lsk_1097 = buffer.data(lsk + 1097);
    const auto *lsk_1098 = buffer.data(lsk + 1098);
    const auto *lsk_1100 = buffer.data(lsk + 1100);
    const auto *lsk_1101 = buffer.data(lsk + 1101);
    const auto *lsk_1103 = buffer.data(lsk + 1103);
    const auto *lsk_1104 = buffer.data(lsk + 1104);
    const auto *lsk_1105 = buffer.data(lsk + 1105);
    const auto *lsk_1107 = buffer.data(lsk + 1107);
    const auto *lsk_1108 = buffer.data(lsk + 1108);
    const auto *lsk_1109 = buffer.data(lsk + 1109);
    const auto *lsk_1110 = buffer.data(lsk + 1110);
    const auto *lsk_1111 = buffer.data(lsk + 1111);
    const auto *lsk_1112 = buffer.data(lsk + 1112);
    const auto *lsk_1113 = buffer.data(lsk + 1113);
    const auto *lsk_1114 = buffer.data(lsk + 1114);
    const auto *lsk_1115 = buffer.data(lsk + 1115);
    const auto *lsk_1116 = buffer.data(lsk + 1116);
    const auto *lsk_1119 = buffer.data(lsk + 1119);
    const auto *lsk_1121 = buffer.data(lsk + 1121);
    const auto *lsk_1122 = buffer.data(lsk + 1122);
    const auto *lsk_1125 = buffer.data(lsk + 1125);
    const auto *lsk_1126 = buffer.data(lsk + 1126);
    const auto *lsk_1128 = buffer.data(lsk + 1128);
    const auto *lsk_1130 = buffer.data(lsk + 1130);
    const auto *lsk_1131 = buffer.data(lsk + 1131);
    const auto *lsk_1133 = buffer.data(lsk + 1133);
    const auto *lsk_1134 = buffer.data(lsk + 1134);
    const auto *lsk_1136 = buffer.data(lsk + 1136);
    const auto *lsk_1137 = buffer.data(lsk + 1137);
    const auto *lsk_1139 = buffer.data(lsk + 1139);
    const auto *lsk_1140 = buffer.data(lsk + 1140);
    const auto *lsk_1141 = buffer.data(lsk + 1141);
    const auto *lsk_1143 = buffer.data(lsk + 1143);
    const auto *lsk_1144 = buffer.data(lsk + 1144);
    const auto *lsk_1145 = buffer.data(lsk + 1145);
    const auto *lsk_1146 = buffer.data(lsk + 1146);
    const auto *lsk_1147 = buffer.data(lsk + 1147);
    const auto *lsk_1148 = buffer.data(lsk + 1148);
    const auto *lsk_1149 = buffer.data(lsk + 1149);
    const auto *lsk_1150 = buffer.data(lsk + 1150);
    const auto *lsk_1151 = buffer.data(lsk + 1151);
    const auto *lsk_1152 = buffer.data(lsk + 1152);
    const auto *lsk_1155 = buffer.data(lsk + 1155);
    const auto *lsk_1157 = buffer.data(lsk + 1157);
    const auto *lsk_1158 = buffer.data(lsk + 1158);
    const auto *lsk_1161 = buffer.data(lsk + 1161);
    const auto *lsk_1162 = buffer.data(lsk + 1162);
    const auto *lsk_1164 = buffer.data(lsk + 1164);
    const auto *lsk_1166 = buffer.data(lsk + 1166);
    const auto *lsk_1167 = buffer.data(lsk + 1167);
    const auto *lsk_1169 = buffer.data(lsk + 1169);
    const auto *lsk_1170 = buffer.data(lsk + 1170);
    const auto *lsk_1172 = buffer.data(lsk + 1172);
    const auto *lsk_1173 = buffer.data(lsk + 1173);
    const auto *lsk_1175 = buffer.data(lsk + 1175);
    const auto *lsk_1176 = buffer.data(lsk + 1176);
    const auto *lsk_1177 = buffer.data(lsk + 1177);

    const auto *msi0_854 = buffer.data(msi0 + 854);
    const auto *msi0_855 = buffer.data(msi0 + 855);
    const auto *msi0_857 = buffer.data(msi0 + 857);
    const auto *msi0_858 = buffer.data(msi0 + 858);
    const auto *msi0_860 = buffer.data(msi0 + 860);
    const auto *msi0_861 = buffer.data(msi0 + 861);
    const auto *msi0_863 = buffer.data(msi0 + 863);
    const auto *msi0_864 = buffer.data(msi0 + 864);
    const auto *msi0_865 = buffer.data(msi0 + 865);
    const auto *msi0_866 = buffer.data(msi0 + 866);
    const auto *msi0_867 = buffer.data(msi0 + 867);
    const auto *msi0_868 = buffer.data(msi0 + 868);
    const auto *msi0_871 = buffer.data(msi0 + 871);
    const auto *msi0_873 = buffer.data(msi0 + 873);
    const auto *msi0_874 = buffer.data(msi0 + 874);
    const auto *msi0_877 = buffer.data(msi0 + 877);
    const auto *msi0_878 = buffer.data(msi0 + 878);
    const auto *msi0_880 = buffer.data(msi0 + 880);
    const auto *msi0_882 = buffer.data(msi0 + 882);
    const auto *msi0_883 = buffer.data(msi0 + 883);
    const auto *msi0_885 = buffer.data(msi0 + 885);
    const auto *msi0_886 = buffer.data(msi0 + 886);
    const auto *msi0_888 = buffer.data(msi0 + 888);
    const auto *msi0_889 = buffer.data(msi0 + 889);
    const auto *msi0_891 = buffer.data(msi0 + 891);
    const auto *msi0_892 = buffer.data(msi0 + 892);
    const auto *msi0_893 = buffer.data(msi0 + 893);
    const auto *msi0_894 = buffer.data(msi0 + 894);
    const auto *msi0_895 = buffer.data(msi0 + 895);
    const auto *msi0_896 = buffer.data(msi0 + 896);
    const auto *msi0_899 = buffer.data(msi0 + 899);
    const auto *msi0_901 = buffer.data(msi0 + 901);
    const auto *msi0_902 = buffer.data(msi0 + 902);
    const auto *msi0_905 = buffer.data(msi0 + 905);
    const auto *msi0_906 = buffer.data(msi0 + 906);
    const auto *msi0_908 = buffer.data(msi0 + 908);
    const auto *msi0_910 = buffer.data(msi0 + 910);
    const auto *msi0_911 = buffer.data(msi0 + 911);
    const auto *msi0_913 = buffer.data(msi0 + 913);
    const auto *msi0_914 = buffer.data(msi0 + 914);
    const auto *msi0_916 = buffer.data(msi0 + 916);
    const auto *msi0_917 = buffer.data(msi0 + 917);
    const auto *msi0_919 = buffer.data(msi0 + 919);
    const auto *msi0_920 = buffer.data(msi0 + 920);
    const auto *msi0_921 = buffer.data(msi0 + 921);

    const auto *msi1_854 = buffer.data(msi1 + 854);
    const auto *msi1_855 = buffer.data(msi1 + 855);
    const auto *msi1_857 = buffer.data(msi1 + 857);
    const auto *msi1_858 = buffer.data(msi1 + 858);
    const auto *msi1_860 = buffer.data(msi1 + 860);
    const auto *msi1_861 = buffer.data(msi1 + 861);
    const auto *msi1_863 = buffer.data(msi1 + 863);
    const auto *msi1_864 = buffer.data(msi1 + 864);
    const auto *msi1_865 = buffer.data(msi1 + 865);
    const auto *msi1_866 = buffer.data(msi1 + 866);
    const auto *msi1_867 = buffer.data(msi1 + 867);
    const auto *msi1_868 = buffer.data(msi1 + 868);
    const auto *msi1_871 = buffer.data(msi1 + 871);
    const auto *msi1_873 = buffer.data(msi1 + 873);
    const auto *msi1_874 = buffer.data(msi1 + 874);
    const auto *msi1_877 = buffer.data(msi1 + 877);
    const auto *msi1_878 = buffer.data(msi1 + 878);
    const auto *msi1_880 = buffer.data(msi1 + 880);
    const auto *msi1_882 = buffer.data(msi1 + 882);
    const auto *msi1_883 = buffer.data(msi1 + 883);
    const auto *msi1_885 = buffer.data(msi1 + 885);
    const auto *msi1_886 = buffer.data(msi1 + 886);
    const auto *msi1_888 = buffer.data(msi1 + 888);
    const auto *msi1_889 = buffer.data(msi1 + 889);
    const auto *msi1_891 = buffer.data(msi1 + 891);
    const auto *msi1_892 = buffer.data(msi1 + 892);
    const auto *msi1_893 = buffer.data(msi1 + 893);
    const auto *msi1_894 = buffer.data(msi1 + 894);
    const auto *msi1_895 = buffer.data(msi1 + 895);
    const auto *msi1_896 = buffer.data(msi1 + 896);
    const auto *msi1_899 = buffer.data(msi1 + 899);
    const auto *msi1_901 = buffer.data(msi1 + 901);
    const auto *msi1_902 = buffer.data(msi1 + 902);
    const auto *msi1_905 = buffer.data(msi1 + 905);
    const auto *msi1_906 = buffer.data(msi1 + 906);
    const auto *msi1_908 = buffer.data(msi1 + 908);
    const auto *msi1_910 = buffer.data(msi1 + 910);
    const auto *msi1_911 = buffer.data(msi1 + 911);
    const auto *msi1_913 = buffer.data(msi1 + 913);
    const auto *msi1_914 = buffer.data(msi1 + 914);
    const auto *msi1_916 = buffer.data(msi1 + 916);
    const auto *msi1_917 = buffer.data(msi1 + 917);
    const auto *msi1_919 = buffer.data(msi1 + 919);
    const auto *msi1_920 = buffer.data(msi1 + 920);
    const auto *msi1_921 = buffer.data(msi1 + 921);

    const auto *msk_1090 = buffer.data(msk + 1090);
    const auto *msk_1094 = buffer.data(msk + 1094);
    const auto *msk_1095 = buffer.data(msk + 1095);
    const auto *msk_1097 = buffer.data(msk + 1097);
    const auto *msk_1098 = buffer.data(msk + 1098);
    const auto *msk_1100 = buffer.data(msk + 1100);
    const auto *msk_1101 = buffer.data(msk + 1101);
    const auto *msk_1103 = buffer.data(msk + 1103);
    const auto *msk_1104 = buffer.data(msk + 1104);
    const auto *msk_1105 = buffer.data(msk + 1105);
    const auto *msk_1107 = buffer.data(msk + 1107);
    const auto *msk_1108 = buffer.data(msk + 1108);
    const auto *msk_1109 = buffer.data(msk + 1109);
    const auto *msk_1110 = buffer.data(msk + 1110);
    const auto *msk_1111 = buffer.data(msk + 1111);
    const auto *msk_1112 = buffer.data(msk + 1112);
    const auto *msk_1113 = buffer.data(msk + 1113);
    const auto *msk_1114 = buffer.data(msk + 1114);
    const auto *msk_1115 = buffer.data(msk + 1115);
    const auto *msk_1116 = buffer.data(msk + 1116);
    const auto *msk_1118 = buffer.data(msk + 1118);
    const auto *msk_1119 = buffer.data(msk + 1119);
    const auto *msk_1121 = buffer.data(msk + 1121);
    const auto *msk_1122 = buffer.data(msk + 1122);
    const auto *msk_1125 = buffer.data(msk + 1125);
    const auto *msk_1126 = buffer.data(msk + 1126);
    const auto *msk_1128 = buffer.data(msk + 1128);
    const auto *msk_1130 = buffer.data(msk + 1130);
    const auto *msk_1131 = buffer.data(msk + 1131);
    const auto *msk_1133 = buffer.data(msk + 1133);
    const auto *msk_1134 = buffer.data(msk + 1134);
    const auto *msk_1136 = buffer.data(msk + 1136);
    const auto *msk_1137 = buffer.data(msk + 1137);
    const auto *msk_1139 = buffer.data(msk + 1139);
    const auto *msk_1140 = buffer.data(msk + 1140);
    const auto *msk_1141 = buffer.data(msk + 1141);
    const auto *msk_1143 = buffer.data(msk + 1143);
    const auto *msk_1144 = buffer.data(msk + 1144);
    const auto *msk_1145 = buffer.data(msk + 1145);
    const auto *msk_1146 = buffer.data(msk + 1146);
    const auto *msk_1147 = buffer.data(msk + 1147);
    const auto *msk_1148 = buffer.data(msk + 1148);
    const auto *msk_1149 = buffer.data(msk + 1149);
    const auto *msk_1150 = buffer.data(msk + 1150);
    const auto *msk_1151 = buffer.data(msk + 1151);
    const auto *msk_1152 = buffer.data(msk + 1152);
    const auto *msk_1154 = buffer.data(msk + 1154);
    const auto *msk_1155 = buffer.data(msk + 1155);
    const auto *msk_1157 = buffer.data(msk + 1157);
    const auto *msk_1158 = buffer.data(msk + 1158);
    const auto *msk_1161 = buffer.data(msk + 1161);
    const auto *msk_1162 = buffer.data(msk + 1162);
    const auto *msk_1164 = buffer.data(msk + 1164);
    const auto *msk_1166 = buffer.data(msk + 1166);
    const auto *msk_1167 = buffer.data(msk + 1167);
    const auto *msk_1169 = buffer.data(msk + 1169);
    const auto *msk_1170 = buffer.data(msk + 1170);
    const auto *msk_1172 = buffer.data(msk + 1172);
    const auto *msk_1173 = buffer.data(msk + 1173);
    const auto *msk_1175 = buffer.data(msk + 1175);
    const auto *msk_1176 = buffer.data(msk + 1176);
    const auto *msk_1177 = buffer.data(msk + 1177);

#pragma omp simd aligned(t_1364, t_1365, t_1366, pc_x, pc_z, lsk_802, lsk_1094, lsk_1095, \
                         msi0_854, msi0_855, msi1_854, msi1_855, msk_1090, msk_1094, \
                         msk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = f_16 * lsk_1094[k]
                    + f_8 * msi0_854[k]
                    - f_9 * msi1_854[k]
                    + f_3 * pc_x[k] * msk_1094[k];

        t_1365[k] = f_16 * lsk_1095[k]
                    + f_6 * msi0_855[k]
                    - f_7 * msi1_855[k]
                    + f_3 * pc_x[k] * msk_1095[k];

        t_1366[k] = f_16 * lsk_802[k]
                    + f_3 * pc_z[k] * msk_1090[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, pc_x, pc_y, lsk_842, lsk_1097, lsk_1098, \
                         msi0_857, msi0_858, msi1_857, msi1_858, msk_1094, msk_1097, \
                         msk_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = f_16 * lsk_1097[k]
                    + f_6 * msi0_857[k]
                    - f_7 * msi1_857[k]
                    + f_3 * pc_x[k] * msk_1097[k];

        t_1368[k] = f_16 * lsk_1098[k]
                    + f_6 * msi0_858[k]
                    - f_7 * msi1_858[k]
                    + f_3 * pc_x[k] * msk_1098[k];

        t_1369[k] = f_19 * lsk_842[k]
                    + f_3 * pc_y[k] * msk_1094[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, pc_x, pc_z, lsk_807, lsk_1100, lsk_1101, \
                         msi0_860, msi0_861, msi1_860, msi1_861, msk_1095, msk_1100, \
                         msk_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_16 * lsk_1100[k]
                    + f_6 * msi0_860[k]
                    - f_7 * msi1_860[k]
                    + f_3 * pc_x[k] * msk_1100[k];

        t_1371[k] = f_16 * lsk_1101[k]
                    + f_4 * msi0_861[k]
                    - f_5 * msi1_861[k]
                    + f_3 * pc_x[k] * msk_1101[k];

        t_1372[k] = f_16 * lsk_807[k]
                    + f_3 * pc_z[k] * msk_1095[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pc_x, lsk_1103, lsk_1104, lsk_1105, msi0_863, \
                         msi0_864, msi0_865, msi1_863, msi1_864, msi1_865, msk_1103, msk_1104, \
                         msk_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = f_16 * lsk_1103[k]
                    + f_4 * msi0_863[k]
                    - f_5 * msi1_863[k]
                    + f_3 * pc_x[k] * msk_1103[k];

        t_1374[k] = f_16 * lsk_1104[k]
                    + f_4 * msi0_864[k]
                    - f_5 * msi1_864[k]
                    + f_3 * pc_x[k] * msk_1104[k];

        t_1375[k] = f_16 * lsk_1105[k]
                    + f_4 * msi0_865[k]
                    - f_5 * msi1_865[k]
                    + f_3 * pc_x[k] * msk_1105[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, t_1379, pc_x, pc_y, lsk_848, lsk_1107, \
                         lsk_1108, lsk_1109, msi0_867, msi1_867, msk_1100, msk_1107, msk_1108, \
                         msk_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_19 * lsk_848[k]
                    + f_3 * pc_y[k] * msk_1100[k];

        t_1377[k] = f_16 * lsk_1107[k]
                    + f_4 * msi0_867[k]
                    - f_5 * msi1_867[k]
                    + f_3 * pc_x[k] * msk_1107[k];

        t_1378[k] = f_16 * lsk_1108[k]
                    + f_3 * pc_x[k] * msk_1108[k];

        t_1379[k] = f_16 * lsk_1109[k]
                    + f_3 * pc_x[k] * msk_1109[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, t_1384, pc_x, lsk_1110, lsk_1111, \
                         lsk_1112, lsk_1113, lsk_1114, msk_1110, msk_1111, msk_1112, msk_1113, \
                         msk_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_16 * lsk_1110[k]
                    + f_3 * pc_x[k] * msk_1110[k];

        t_1381[k] = f_16 * lsk_1111[k]
                    + f_3 * pc_x[k] * msk_1111[k];

        t_1382[k] = f_16 * lsk_1112[k]
                    + f_3 * pc_x[k] * msk_1112[k];

        t_1383[k] = f_16 * lsk_1113[k]
                    + f_3 * pc_x[k] * msk_1113[k];

        t_1384[k] = f_16 * lsk_1114[k]
                    + f_3 * pc_x[k] * msk_1114[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, pc_x, pc_y, pc_z, lsk_820, lsk_856, lsk_1115, \
                         msi0_861, msi1_861, msk_1108, msk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_16 * lsk_1115[k]
                    + f_3 * pc_x[k] * msk_1115[k];

        t_1386[k] = f_19 * lsk_856[k]
                    + f_1 * msi0_861[k]
                    - f_2 * msi1_861[k]
                    + f_3 * pc_y[k] * msk_1108[k];

        t_1387[k] = f_16 * lsk_820[k]
                    + f_3 * pc_z[k] * msk_1108[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, pc_y, lsk_858, lsk_859, lsk_860, msi0_863, \
                         msi0_864, msi0_865, msi1_863, msi1_864, msi1_865, msk_1110, msk_1111, \
                         msk_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_19 * lsk_858[k]
                    + f_12 * msi0_863[k]
                    - f_13 * msi1_863[k]
                    + f_3 * pc_y[k] * msk_1110[k];

        t_1389[k] = f_19 * lsk_859[k]
                    + f_10 * msi0_864[k]
                    - f_11 * msi1_864[k]
                    + f_3 * pc_y[k] * msk_1111[k];

        t_1390[k] = f_19 * lsk_860[k]
                    + f_8 * msi0_865[k]
                    - f_9 * msi1_865[k]
                    + f_3 * pc_y[k] * msk_1112[k];
    }

#pragma omp simd aligned(t_1391, t_1392, t_1393, pc_y, lsk_861, lsk_862, lsk_863, msi0_866, \
                         msi0_867, msi1_866, msi1_867, msk_1113, msk_1114, \
                         msk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1391[k] = f_19 * lsk_861[k]
                    + f_6 * msi0_866[k]
                    - f_7 * msi1_866[k]
                    + f_3 * pc_y[k] * msk_1113[k];

        t_1392[k] = f_19 * lsk_862[k]
                    + f_4 * msi0_867[k]
                    - f_5 * msi1_867[k]
                    + f_3 * pc_y[k] * msk_1114[k];

        t_1393[k] = f_19 * lsk_863[k]
                    + f_3 * pc_y[k] * msk_1115[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, pc_x, pc_y, pc_z, lsk_827, lsk_864, lsk_1116, \
                         msi0_867, msi0_868, msi1_867, msi1_868, msk_1115, \
                         msk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_16 * lsk_827[k]
                    + f_1 * msi0_867[k]
                    - f_2 * msi1_867[k]
                    + f_3 * pc_z[k] * msk_1115[k];

        t_1395[k] = f_16 * lsk_1116[k]
                    + f_1 * msi0_868[k]
                    - f_2 * msi1_868[k]
                    + f_3 * pc_x[k] * msk_1116[k];

        t_1396[k] = f_18 * lsk_864[k]
                    + f_3 * pc_y[k] * msk_1116[k];
    }

#pragma omp simd aligned(t_1397, t_1398, t_1399, pc_x, pc_y, pc_z, lsk_828, lsk_866, lsk_1119, \
                         msi0_871, msi1_871, msk_1116, msk_1118, \
                         msk_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1397[k] = f_17 * lsk_828[k]
                    + f_3 * pc_z[k] * msk_1116[k];

        t_1398[k] = f_16 * lsk_1119[k]
                    + f_12 * msi0_871[k]
                    - f_13 * msi1_871[k]
                    + f_3 * pc_x[k] * msk_1119[k];

        t_1399[k] = f_18 * lsk_866[k]
                    + f_3 * pc_y[k] * msk_1118[k];
    }

#pragma omp simd aligned(t_1400, t_1401, t_1402, pc_x, pc_z, lsk_831, lsk_1121, lsk_1122, \
                         msi0_873, msi0_874, msi1_873, msi1_874, msk_1119, msk_1121, \
                         msk_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1400[k] = f_16 * lsk_1121[k]
                    + f_12 * msi0_873[k]
                    - f_13 * msi1_873[k]
                    + f_3 * pc_x[k] * msk_1121[k];

        t_1401[k] = f_16 * lsk_1122[k]
                    + f_10 * msi0_874[k]
                    - f_11 * msi1_874[k]
                    + f_3 * pc_x[k] * msk_1122[k];

        t_1402[k] = f_17 * lsk_831[k]
                    + f_3 * pc_z[k] * msk_1119[k];
    }

#pragma omp simd aligned(t_1403, t_1404, t_1405, pc_x, pc_y, lsk_869, lsk_1125, lsk_1126, \
                         msi0_877, msi0_878, msi1_877, msi1_878, msk_1121, msk_1125, \
                         msk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1403[k] = f_18 * lsk_869[k]
                    + f_3 * pc_y[k] * msk_1121[k];

        t_1404[k] = f_16 * lsk_1125[k]
                    + f_10 * msi0_877[k]
                    - f_11 * msi1_877[k]
                    + f_3 * pc_x[k] * msk_1125[k];

        t_1405[k] = f_16 * lsk_1126[k]
                    + f_8 * msi0_878[k]
                    - f_9 * msi1_878[k]
                    + f_3 * pc_x[k] * msk_1126[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pc_x, pc_y, pc_z, lsk_834, lsk_873, lsk_1128, \
                         msi0_880, msi1_880, msk_1122, msk_1125, \
                         msk_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_17 * lsk_834[k]
                    + f_3 * pc_z[k] * msk_1122[k];

        t_1407[k] = f_16 * lsk_1128[k]
                    + f_8 * msi0_880[k]
                    - f_9 * msi1_880[k]
                    + f_3 * pc_x[k] * msk_1128[k];

        t_1408[k] = f_18 * lsk_873[k]
                    + f_3 * pc_y[k] * msk_1125[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pc_x, pc_z, lsk_838, lsk_1130, lsk_1131, \
                         msi0_882, msi0_883, msi1_882, msi1_883, msk_1126, msk_1130, \
                         msk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_16 * lsk_1130[k]
                    + f_8 * msi0_882[k]
                    - f_9 * msi1_882[k]
                    + f_3 * pc_x[k] * msk_1130[k];

        t_1410[k] = f_16 * lsk_1131[k]
                    + f_6 * msi0_883[k]
                    - f_7 * msi1_883[k]
                    + f_3 * pc_x[k] * msk_1131[k];

        t_1411[k] = f_17 * lsk_838[k]
                    + f_3 * pc_z[k] * msk_1126[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pc_x, pc_y, lsk_878, lsk_1133, lsk_1134, \
                         msi0_885, msi0_886, msi1_885, msi1_886, msk_1130, msk_1133, \
                         msk_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_16 * lsk_1133[k]
                    + f_6 * msi0_885[k]
                    - f_7 * msi1_885[k]
                    + f_3 * pc_x[k] * msk_1133[k];

        t_1413[k] = f_16 * lsk_1134[k]
                    + f_6 * msi0_886[k]
                    - f_7 * msi1_886[k]
                    + f_3 * pc_x[k] * msk_1134[k];

        t_1414[k] = f_18 * lsk_878[k]
                    + f_3 * pc_y[k] * msk_1130[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pc_x, pc_z, lsk_843, lsk_1136, lsk_1137, \
                         msi0_888, msi0_889, msi1_888, msi1_889, msk_1131, msk_1136, \
                         msk_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_16 * lsk_1136[k]
                    + f_6 * msi0_888[k]
                    - f_7 * msi1_888[k]
                    + f_3 * pc_x[k] * msk_1136[k];

        t_1416[k] = f_16 * lsk_1137[k]
                    + f_4 * msi0_889[k]
                    - f_5 * msi1_889[k]
                    + f_3 * pc_x[k] * msk_1137[k];

        t_1417[k] = f_17 * lsk_843[k]
                    + f_3 * pc_z[k] * msk_1131[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pc_x, lsk_1139, lsk_1140, lsk_1141, msi0_891, \
                         msi0_892, msi0_893, msi1_891, msi1_892, msi1_893, msk_1139, msk_1140, \
                         msk_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_16 * lsk_1139[k]
                    + f_4 * msi0_891[k]
                    - f_5 * msi1_891[k]
                    + f_3 * pc_x[k] * msk_1139[k];

        t_1419[k] = f_16 * lsk_1140[k]
                    + f_4 * msi0_892[k]
                    - f_5 * msi1_892[k]
                    + f_3 * pc_x[k] * msk_1140[k];

        t_1420[k] = f_16 * lsk_1141[k]
                    + f_4 * msi0_893[k]
                    - f_5 * msi1_893[k]
                    + f_3 * pc_x[k] * msk_1141[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, t_1424, pc_x, pc_y, lsk_884, lsk_1143, \
                         lsk_1144, lsk_1145, msi0_895, msi1_895, msk_1136, msk_1143, msk_1144, \
                         msk_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_18 * lsk_884[k]
                    + f_3 * pc_y[k] * msk_1136[k];

        t_1422[k] = f_16 * lsk_1143[k]
                    + f_4 * msi0_895[k]
                    - f_5 * msi1_895[k]
                    + f_3 * pc_x[k] * msk_1143[k];

        t_1423[k] = f_16 * lsk_1144[k]
                    + f_3 * pc_x[k] * msk_1144[k];

        t_1424[k] = f_16 * lsk_1145[k]
                    + f_3 * pc_x[k] * msk_1145[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, t_1428, t_1429, pc_x, lsk_1146, lsk_1147, \
                         lsk_1148, lsk_1149, lsk_1150, msk_1146, msk_1147, msk_1148, msk_1149, \
                         msk_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_16 * lsk_1146[k]
                    + f_3 * pc_x[k] * msk_1146[k];

        t_1426[k] = f_16 * lsk_1147[k]
                    + f_3 * pc_x[k] * msk_1147[k];

        t_1427[k] = f_16 * lsk_1148[k]
                    + f_3 * pc_x[k] * msk_1148[k];

        t_1428[k] = f_16 * lsk_1149[k]
                    + f_3 * pc_x[k] * msk_1149[k];

        t_1429[k] = f_16 * lsk_1150[k]
                    + f_3 * pc_x[k] * msk_1150[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pc_x, pc_y, pc_z, lsk_856, lsk_892, lsk_1151, \
                         msi0_889, msi1_889, msk_1144, msk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_16 * lsk_1151[k]
                    + f_3 * pc_x[k] * msk_1151[k];

        t_1431[k] = f_18 * lsk_892[k]
                    + f_1 * msi0_889[k]
                    - f_2 * msi1_889[k]
                    + f_3 * pc_y[k] * msk_1144[k];

        t_1432[k] = f_17 * lsk_856[k]
                    + f_3 * pc_z[k] * msk_1144[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pc_y, lsk_894, lsk_895, lsk_896, msi0_891, \
                         msi0_892, msi0_893, msi1_891, msi1_892, msi1_893, msk_1146, msk_1147, \
                         msk_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_18 * lsk_894[k]
                    + f_12 * msi0_891[k]
                    - f_13 * msi1_891[k]
                    + f_3 * pc_y[k] * msk_1146[k];

        t_1434[k] = f_18 * lsk_895[k]
                    + f_10 * msi0_892[k]
                    - f_11 * msi1_892[k]
                    + f_3 * pc_y[k] * msk_1147[k];

        t_1435[k] = f_18 * lsk_896[k]
                    + f_8 * msi0_893[k]
                    - f_9 * msi1_893[k]
                    + f_3 * pc_y[k] * msk_1148[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pc_y, lsk_897, lsk_898, lsk_899, msi0_894, \
                         msi0_895, msi1_894, msi1_895, msk_1149, msk_1150, \
                         msk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_18 * lsk_897[k]
                    + f_6 * msi0_894[k]
                    - f_7 * msi1_894[k]
                    + f_3 * pc_y[k] * msk_1149[k];

        t_1437[k] = f_18 * lsk_898[k]
                    + f_4 * msi0_895[k]
                    - f_5 * msi1_895[k]
                    + f_3 * pc_y[k] * msk_1150[k];

        t_1438[k] = f_18 * lsk_899[k]
                    + f_3 * pc_y[k] * msk_1151[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pc_x, pc_y, pc_z, lsk_863, lsk_900, lsk_1152, \
                         msi0_895, msi0_896, msi1_895, msi1_896, msk_1151, \
                         msk_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_17 * lsk_863[k]
                    + f_1 * msi0_895[k]
                    - f_2 * msi1_895[k]
                    + f_3 * pc_z[k] * msk_1151[k];

        t_1440[k] = f_16 * lsk_1152[k]
                    + f_1 * msi0_896[k]
                    - f_2 * msi1_896[k]
                    + f_3 * pc_x[k] * msk_1152[k];

        t_1441[k] = f_17 * lsk_900[k]
                    + f_3 * pc_y[k] * msk_1152[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pc_x, pc_y, pc_z, lsk_864, lsk_902, lsk_1155, \
                         msi0_899, msi1_899, msk_1152, msk_1154, \
                         msk_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_18 * lsk_864[k]
                    + f_3 * pc_z[k] * msk_1152[k];

        t_1443[k] = f_16 * lsk_1155[k]
                    + f_12 * msi0_899[k]
                    - f_13 * msi1_899[k]
                    + f_3 * pc_x[k] * msk_1155[k];

        t_1444[k] = f_17 * lsk_902[k]
                    + f_3 * pc_y[k] * msk_1154[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_x, pc_z, lsk_867, lsk_1157, lsk_1158, \
                         msi0_901, msi0_902, msi1_901, msi1_902, msk_1155, msk_1157, \
                         msk_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_16 * lsk_1157[k]
                    + f_12 * msi0_901[k]
                    - f_13 * msi1_901[k]
                    + f_3 * pc_x[k] * msk_1157[k];

        t_1446[k] = f_16 * lsk_1158[k]
                    + f_10 * msi0_902[k]
                    - f_11 * msi1_902[k]
                    + f_3 * pc_x[k] * msk_1158[k];

        t_1447[k] = f_18 * lsk_867[k]
                    + f_3 * pc_z[k] * msk_1155[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, lsk_905, lsk_1161, lsk_1162, \
                         msi0_905, msi0_906, msi1_905, msi1_906, msk_1157, msk_1161, \
                         msk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_17 * lsk_905[k]
                    + f_3 * pc_y[k] * msk_1157[k];

        t_1449[k] = f_16 * lsk_1161[k]
                    + f_10 * msi0_905[k]
                    - f_11 * msi1_905[k]
                    + f_3 * pc_x[k] * msk_1161[k];

        t_1450[k] = f_16 * lsk_1162[k]
                    + f_8 * msi0_906[k]
                    - f_9 * msi1_906[k]
                    + f_3 * pc_x[k] * msk_1162[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, pc_y, pc_z, lsk_870, lsk_909, lsk_1164, \
                         msi0_908, msi1_908, msk_1158, msk_1161, \
                         msk_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_18 * lsk_870[k]
                    + f_3 * pc_z[k] * msk_1158[k];

        t_1452[k] = f_16 * lsk_1164[k]
                    + f_8 * msi0_908[k]
                    - f_9 * msi1_908[k]
                    + f_3 * pc_x[k] * msk_1164[k];

        t_1453[k] = f_17 * lsk_909[k]
                    + f_3 * pc_y[k] * msk_1161[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_z, lsk_874, lsk_1166, lsk_1167, \
                         msi0_910, msi0_911, msi1_910, msi1_911, msk_1162, msk_1166, \
                         msk_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_16 * lsk_1166[k]
                    + f_8 * msi0_910[k]
                    - f_9 * msi1_910[k]
                    + f_3 * pc_x[k] * msk_1166[k];

        t_1455[k] = f_16 * lsk_1167[k]
                    + f_6 * msi0_911[k]
                    - f_7 * msi1_911[k]
                    + f_3 * pc_x[k] * msk_1167[k];

        t_1456[k] = f_18 * lsk_874[k]
                    + f_3 * pc_z[k] * msk_1162[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, pc_y, lsk_914, lsk_1169, lsk_1170, \
                         msi0_913, msi0_914, msi1_913, msi1_914, msk_1166, msk_1169, \
                         msk_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_16 * lsk_1169[k]
                    + f_6 * msi0_913[k]
                    - f_7 * msi1_913[k]
                    + f_3 * pc_x[k] * msk_1169[k];

        t_1458[k] = f_16 * lsk_1170[k]
                    + f_6 * msi0_914[k]
                    - f_7 * msi1_914[k]
                    + f_3 * pc_x[k] * msk_1170[k];

        t_1459[k] = f_17 * lsk_914[k]
                    + f_3 * pc_y[k] * msk_1166[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pc_x, pc_z, lsk_879, lsk_1172, lsk_1173, \
                         msi0_916, msi0_917, msi1_916, msi1_917, msk_1167, msk_1172, \
                         msk_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_16 * lsk_1172[k]
                    + f_6 * msi0_916[k]
                    - f_7 * msi1_916[k]
                    + f_3 * pc_x[k] * msk_1172[k];

        t_1461[k] = f_16 * lsk_1173[k]
                    + f_4 * msi0_917[k]
                    - f_5 * msi1_917[k]
                    + f_3 * pc_x[k] * msk_1173[k];

        t_1462[k] = f_18 * lsk_879[k]
                    + f_3 * pc_z[k] * msk_1167[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pc_x, lsk_1175, lsk_1176, lsk_1177, msi0_919, \
                         msi0_920, msi0_921, msi1_919, msi1_920, msi1_921, msk_1175, msk_1176, \
                         msk_1177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = f_16 * lsk_1175[k]
                    + f_4 * msi0_919[k]
                    - f_5 * msi1_919[k]
                    + f_3 * pc_x[k] * msk_1175[k];

        t_1464[k] = f_16 * lsk_1176[k]
                    + f_4 * msi0_920[k]
                    - f_5 * msi1_920[k]
                    + f_3 * pc_x[k] * msk_1176[k];

        t_1465[k] = f_16 * lsk_1177[k]
                    + f_4 * msi0_921[k]
                    - f_5 * msi1_921[k]
                    + f_3 * pc_x[k] * msk_1177[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1215 = buffer.data(lsl0 + 1215);
    const auto *lsl0_1218 = buffer.data(lsl0 + 1218);
    const auto *lsl0_1220 = buffer.data(lsl0 + 1220);
    const auto *lsl0_1221 = buffer.data(lsl0 + 1221);
    const auto *lsl0_1224 = buffer.data(lsl0 + 1224);
    const auto *lsl0_1225 = buffer.data(lsl0 + 1225);
    const auto *lsl0_1227 = buffer.data(lsl0 + 1227);
    const auto *lsl0_1229 = buffer.data(lsl0 + 1229);
    const auto *lsl0_1230 = buffer.data(lsl0 + 1230);
    const auto *lsl0_1232 = buffer.data(lsl0 + 1232);
    const auto *lsl0_1233 = buffer.data(lsl0 + 1233);
    const auto *lsl0_1235 = buffer.data(lsl0 + 1235);
    const auto *lsl0_1236 = buffer.data(lsl0 + 1236);
    const auto *lsl0_1238 = buffer.data(lsl0 + 1238);
    const auto *lsl0_1239 = buffer.data(lsl0 + 1239);
    const auto *lsl0_1240 = buffer.data(lsl0 + 1240);
    const auto *lsl0_1242 = buffer.data(lsl0 + 1242);

    const auto *lsk_892 = buffer.data(lsk + 892);
    const auto *lsk_899 = buffer.data(lsk + 899);
    const auto *lsk_900 = buffer.data(lsk + 900);
    const auto *lsk_903 = buffer.data(lsk + 903);
    const auto *lsk_906 = buffer.data(lsk + 906);
    const auto *lsk_910 = buffer.data(lsk + 910);
    const auto *lsk_915 = buffer.data(lsk + 915);
    const auto *lsk_920 = buffer.data(lsk + 920);
    const auto *lsk_928 = buffer.data(lsk + 928);
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
    const auto *lsk_977 = buffer.data(lsk + 977);
    const auto *lsk_978 = buffer.data(lsk + 978);
    const auto *lsk_980 = buffer.data(lsk + 980);
    const auto *lsk_981 = buffer.data(lsk + 981);
    const auto *lsk_982 = buffer.data(lsk + 982);
    const auto *lsk_984 = buffer.data(lsk + 984);
    const auto *lsk_985 = buffer.data(lsk + 985);
    const auto *lsk_986 = buffer.data(lsk + 986);
    const auto *lsk_987 = buffer.data(lsk + 987);
    const auto *lsk_989 = buffer.data(lsk + 989);
    const auto *lsk_990 = buffer.data(lsk + 990);
    const auto *lsk_991 = buffer.data(lsk + 991);
    const auto *lsk_992 = buffer.data(lsk + 992);
    const auto *lsk_1000 = buffer.data(lsk + 1000);
    const auto *lsk_1002 = buffer.data(lsk + 1002);
    const auto *lsk_1003 = buffer.data(lsk + 1003);
    const auto *lsk_1004 = buffer.data(lsk + 1004);
    const auto *lsk_1005 = buffer.data(lsk + 1005);
    const auto *lsk_1006 = buffer.data(lsk + 1006);
    const auto *lsk_1007 = buffer.data(lsk + 1007);
    const auto *lsk_1179 = buffer.data(lsk + 1179);
    const auto *lsk_1180 = buffer.data(lsk + 1180);
    const auto *lsk_1181 = buffer.data(lsk + 1181);
    const auto *lsk_1182 = buffer.data(lsk + 1182);
    const auto *lsk_1183 = buffer.data(lsk + 1183);
    const auto *lsk_1184 = buffer.data(lsk + 1184);
    const auto *lsk_1185 = buffer.data(lsk + 1185);
    const auto *lsk_1186 = buffer.data(lsk + 1186);
    const auto *lsk_1187 = buffer.data(lsk + 1187);
    const auto *lsk_1188 = buffer.data(lsk + 1188);
    const auto *lsk_1191 = buffer.data(lsk + 1191);
    const auto *lsk_1193 = buffer.data(lsk + 1193);
    const auto *lsk_1194 = buffer.data(lsk + 1194);
    const auto *lsk_1197 = buffer.data(lsk + 1197);
    const auto *lsk_1198 = buffer.data(lsk + 1198);
    const auto *lsk_1200 = buffer.data(lsk + 1200);
    const auto *lsk_1202 = buffer.data(lsk + 1202);
    const auto *lsk_1203 = buffer.data(lsk + 1203);
    const auto *lsk_1205 = buffer.data(lsk + 1205);
    const auto *lsk_1206 = buffer.data(lsk + 1206);
    const auto *lsk_1208 = buffer.data(lsk + 1208);
    const auto *lsk_1209 = buffer.data(lsk + 1209);
    const auto *lsk_1211 = buffer.data(lsk + 1211);
    const auto *lsk_1212 = buffer.data(lsk + 1212);
    const auto *lsk_1213 = buffer.data(lsk + 1213);
    const auto *lsk_1215 = buffer.data(lsk + 1215);
    const auto *lsk_1216 = buffer.data(lsk + 1216);
    const auto *lsk_1217 = buffer.data(lsk + 1217);
    const auto *lsk_1218 = buffer.data(lsk + 1218);
    const auto *lsk_1219 = buffer.data(lsk + 1219);
    const auto *lsk_1220 = buffer.data(lsk + 1220);
    const auto *lsk_1221 = buffer.data(lsk + 1221);
    const auto *lsk_1222 = buffer.data(lsk + 1222);
    const auto *lsk_1223 = buffer.data(lsk + 1223);
    const auto *lsk_1252 = buffer.data(lsk + 1252);
    const auto *lsk_1253 = buffer.data(lsk + 1253);
    const auto *lsk_1254 = buffer.data(lsk + 1254);
    const auto *lsk_1255 = buffer.data(lsk + 1255);
    const auto *lsk_1256 = buffer.data(lsk + 1256);
    const auto *lsk_1257 = buffer.data(lsk + 1257);
    const auto *lsk_1258 = buffer.data(lsk + 1258);
    const auto *lsk_1259 = buffer.data(lsk + 1259);

    const auto *lsl1_1215 = buffer.data(lsl1 + 1215);
    const auto *lsl1_1218 = buffer.data(lsl1 + 1218);
    const auto *lsl1_1220 = buffer.data(lsl1 + 1220);
    const auto *lsl1_1221 = buffer.data(lsl1 + 1221);
    const auto *lsl1_1224 = buffer.data(lsl1 + 1224);
    const auto *lsl1_1225 = buffer.data(lsl1 + 1225);
    const auto *lsl1_1227 = buffer.data(lsl1 + 1227);
    const auto *lsl1_1229 = buffer.data(lsl1 + 1229);
    const auto *lsl1_1230 = buffer.data(lsl1 + 1230);
    const auto *lsl1_1232 = buffer.data(lsl1 + 1232);
    const auto *lsl1_1233 = buffer.data(lsl1 + 1233);
    const auto *lsl1_1235 = buffer.data(lsl1 + 1235);
    const auto *lsl1_1236 = buffer.data(lsl1 + 1236);
    const auto *lsl1_1238 = buffer.data(lsl1 + 1238);
    const auto *lsl1_1239 = buffer.data(lsl1 + 1239);
    const auto *lsl1_1240 = buffer.data(lsl1 + 1240);
    const auto *lsl1_1242 = buffer.data(lsl1 + 1242);

    const auto *msi0_917 = buffer.data(msi0 + 917);
    const auto *msi0_919 = buffer.data(msi0 + 919);
    const auto *msi0_920 = buffer.data(msi0 + 920);
    const auto *msi0_921 = buffer.data(msi0 + 921);
    const auto *msi0_922 = buffer.data(msi0 + 922);
    const auto *msi0_923 = buffer.data(msi0 + 923);
    const auto *msi0_924 = buffer.data(msi0 + 924);
    const auto *msi0_927 = buffer.data(msi0 + 927);
    const auto *msi0_929 = buffer.data(msi0 + 929);
    const auto *msi0_930 = buffer.data(msi0 + 930);
    const auto *msi0_933 = buffer.data(msi0 + 933);
    const auto *msi0_934 = buffer.data(msi0 + 934);
    const auto *msi0_936 = buffer.data(msi0 + 936);
    const auto *msi0_938 = buffer.data(msi0 + 938);
    const auto *msi0_939 = buffer.data(msi0 + 939);
    const auto *msi0_941 = buffer.data(msi0 + 941);
    const auto *msi0_942 = buffer.data(msi0 + 942);
    const auto *msi0_944 = buffer.data(msi0 + 944);
    const auto *msi0_945 = buffer.data(msi0 + 945);
    const auto *msi0_947 = buffer.data(msi0 + 947);
    const auto *msi0_948 = buffer.data(msi0 + 948);
    const auto *msi0_949 = buffer.data(msi0 + 949);
    const auto *msi0_950 = buffer.data(msi0 + 950);
    const auto *msi0_951 = buffer.data(msi0 + 951);
    const auto *msi0_973 = buffer.data(msi0 + 973);
    const auto *msi0_975 = buffer.data(msi0 + 975);
    const auto *msi0_976 = buffer.data(msi0 + 976);
    const auto *msi0_977 = buffer.data(msi0 + 977);
    const auto *msi0_978 = buffer.data(msi0 + 978);
    const auto *msi0_979 = buffer.data(msi0 + 979);

    const auto *msi1_917 = buffer.data(msi1 + 917);
    const auto *msi1_919 = buffer.data(msi1 + 919);
    const auto *msi1_920 = buffer.data(msi1 + 920);
    const auto *msi1_921 = buffer.data(msi1 + 921);
    const auto *msi1_922 = buffer.data(msi1 + 922);
    const auto *msi1_923 = buffer.data(msi1 + 923);
    const auto *msi1_924 = buffer.data(msi1 + 924);
    const auto *msi1_927 = buffer.data(msi1 + 927);
    const auto *msi1_929 = buffer.data(msi1 + 929);
    const auto *msi1_930 = buffer.data(msi1 + 930);
    const auto *msi1_933 = buffer.data(msi1 + 933);
    const auto *msi1_934 = buffer.data(msi1 + 934);
    const auto *msi1_936 = buffer.data(msi1 + 936);
    const auto *msi1_938 = buffer.data(msi1 + 938);
    const auto *msi1_939 = buffer.data(msi1 + 939);
    const auto *msi1_941 = buffer.data(msi1 + 941);
    const auto *msi1_942 = buffer.data(msi1 + 942);
    const auto *msi1_944 = buffer.data(msi1 + 944);
    const auto *msi1_945 = buffer.data(msi1 + 945);
    const auto *msi1_947 = buffer.data(msi1 + 947);
    const auto *msi1_948 = buffer.data(msi1 + 948);
    const auto *msi1_949 = buffer.data(msi1 + 949);
    const auto *msi1_950 = buffer.data(msi1 + 950);
    const auto *msi1_951 = buffer.data(msi1 + 951);
    const auto *msi1_973 = buffer.data(msi1 + 973);
    const auto *msi1_975 = buffer.data(msi1 + 975);
    const auto *msi1_976 = buffer.data(msi1 + 976);
    const auto *msi1_977 = buffer.data(msi1 + 977);
    const auto *msi1_978 = buffer.data(msi1 + 978);
    const auto *msi1_979 = buffer.data(msi1 + 979);

    const auto *msk_1172 = buffer.data(msk + 1172);
    const auto *msk_1179 = buffer.data(msk + 1179);
    const auto *msk_1180 = buffer.data(msk + 1180);
    const auto *msk_1181 = buffer.data(msk + 1181);
    const auto *msk_1182 = buffer.data(msk + 1182);
    const auto *msk_1183 = buffer.data(msk + 1183);
    const auto *msk_1184 = buffer.data(msk + 1184);
    const auto *msk_1185 = buffer.data(msk + 1185);
    const auto *msk_1186 = buffer.data(msk + 1186);
    const auto *msk_1187 = buffer.data(msk + 1187);
    const auto *msk_1188 = buffer.data(msk + 1188);
    const auto *msk_1190 = buffer.data(msk + 1190);
    const auto *msk_1191 = buffer.data(msk + 1191);
    const auto *msk_1193 = buffer.data(msk + 1193);
    const auto *msk_1194 = buffer.data(msk + 1194);
    const auto *msk_1197 = buffer.data(msk + 1197);
    const auto *msk_1198 = buffer.data(msk + 1198);
    const auto *msk_1200 = buffer.data(msk + 1200);
    const auto *msk_1202 = buffer.data(msk + 1202);
    const auto *msk_1203 = buffer.data(msk + 1203);
    const auto *msk_1205 = buffer.data(msk + 1205);
    const auto *msk_1206 = buffer.data(msk + 1206);
    const auto *msk_1208 = buffer.data(msk + 1208);
    const auto *msk_1209 = buffer.data(msk + 1209);
    const auto *msk_1211 = buffer.data(msk + 1211);
    const auto *msk_1212 = buffer.data(msk + 1212);
    const auto *msk_1213 = buffer.data(msk + 1213);
    const auto *msk_1215 = buffer.data(msk + 1215);
    const auto *msk_1216 = buffer.data(msk + 1216);
    const auto *msk_1217 = buffer.data(msk + 1217);
    const auto *msk_1218 = buffer.data(msk + 1218);
    const auto *msk_1219 = buffer.data(msk + 1219);
    const auto *msk_1220 = buffer.data(msk + 1220);
    const auto *msk_1221 = buffer.data(msk + 1221);
    const auto *msk_1222 = buffer.data(msk + 1222);
    const auto *msk_1223 = buffer.data(msk + 1223);
    const auto *msk_1224 = buffer.data(msk + 1224);
    const auto *msk_1226 = buffer.data(msk + 1226);
    const auto *msk_1227 = buffer.data(msk + 1227);
    const auto *msk_1229 = buffer.data(msk + 1229);
    const auto *msk_1230 = buffer.data(msk + 1230);
    const auto *msk_1233 = buffer.data(msk + 1233);
    const auto *msk_1234 = buffer.data(msk + 1234);
    const auto *msk_1238 = buffer.data(msk + 1238);
    const auto *msk_1239 = buffer.data(msk + 1239);
    const auto *msk_1244 = buffer.data(msk + 1244);
    const auto *msk_1252 = buffer.data(msk + 1252);
    const auto *msk_1253 = buffer.data(msk + 1253);
    const auto *msk_1254 = buffer.data(msk + 1254);
    const auto *msk_1255 = buffer.data(msk + 1255);
    const auto *msk_1256 = buffer.data(msk + 1256);
    const auto *msk_1257 = buffer.data(msk + 1257);
    const auto *msk_1258 = buffer.data(msk + 1258);
    const auto *msk_1259 = buffer.data(msk + 1259);

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pc_x, pc_y, lsk_920, lsk_1179, \
                         lsk_1180, lsk_1181, msi0_923, msi1_923, msk_1172, msk_1179, msk_1180, \
                         msk_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_17 * lsk_920[k]
                    + f_3 * pc_y[k] * msk_1172[k];

        t_1467[k] = f_16 * lsk_1179[k]
                    + f_4 * msi0_923[k]
                    - f_5 * msi1_923[k]
                    + f_3 * pc_x[k] * msk_1179[k];

        t_1468[k] = f_16 * lsk_1180[k]
                    + f_3 * pc_x[k] * msk_1180[k];

        t_1469[k] = f_16 * lsk_1181[k]
                    + f_3 * pc_x[k] * msk_1181[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, t_1474, pc_x, lsk_1182, lsk_1183, \
                         lsk_1184, lsk_1185, lsk_1186, msk_1182, msk_1183, msk_1184, msk_1185, \
                         msk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_16 * lsk_1182[k]
                    + f_3 * pc_x[k] * msk_1182[k];

        t_1471[k] = f_16 * lsk_1183[k]
                    + f_3 * pc_x[k] * msk_1183[k];

        t_1472[k] = f_16 * lsk_1184[k]
                    + f_3 * pc_x[k] * msk_1184[k];

        t_1473[k] = f_16 * lsk_1185[k]
                    + f_3 * pc_x[k] * msk_1185[k];

        t_1474[k] = f_16 * lsk_1186[k]
                    + f_3 * pc_x[k] * msk_1186[k];
    }

#pragma omp simd aligned(t_1475, t_1476, t_1477, pc_x, pc_y, pc_z, lsk_892, lsk_928, lsk_1187, \
                         msi0_917, msi1_917, msk_1180, msk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1475[k] = f_16 * lsk_1187[k]
                    + f_3 * pc_x[k] * msk_1187[k];

        t_1476[k] = f_17 * lsk_928[k]
                    + f_1 * msi0_917[k]
                    - f_2 * msi1_917[k]
                    + f_3 * pc_y[k] * msk_1180[k];

        t_1477[k] = f_18 * lsk_892[k]
                    + f_3 * pc_z[k] * msk_1180[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_y, lsk_930, lsk_931, lsk_932, msi0_919, \
                         msi0_920, msi0_921, msi1_919, msi1_920, msi1_921, msk_1182, msk_1183, \
                         msk_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_17 * lsk_930[k]
                    + f_12 * msi0_919[k]
                    - f_13 * msi1_919[k]
                    + f_3 * pc_y[k] * msk_1182[k];

        t_1479[k] = f_17 * lsk_931[k]
                    + f_10 * msi0_920[k]
                    - f_11 * msi1_920[k]
                    + f_3 * pc_y[k] * msk_1183[k];

        t_1480[k] = f_17 * lsk_932[k]
                    + f_8 * msi0_921[k]
                    - f_9 * msi1_921[k]
                    + f_3 * pc_y[k] * msk_1184[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, pc_y, lsk_933, lsk_934, lsk_935, msi0_922, \
                         msi0_923, msi1_922, msi1_923, msk_1185, msk_1186, \
                         msk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_17 * lsk_933[k]
                    + f_6 * msi0_922[k]
                    - f_7 * msi1_922[k]
                    + f_3 * pc_y[k] * msk_1185[k];

        t_1482[k] = f_17 * lsk_934[k]
                    + f_4 * msi0_923[k]
                    - f_5 * msi1_923[k]
                    + f_3 * pc_y[k] * msk_1186[k];

        t_1483[k] = f_17 * lsk_935[k]
                    + f_3 * pc_y[k] * msk_1187[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, pc_x, pc_y, pc_z, lsk_899, lsk_936, lsk_1188, \
                         msi0_923, msi0_924, msi1_923, msi1_924, msk_1187, \
                         msk_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_18 * lsk_899[k]
                    + f_1 * msi0_923[k]
                    - f_2 * msi1_923[k]
                    + f_3 * pc_z[k] * msk_1187[k];

        t_1485[k] = f_16 * lsk_1188[k]
                    + f_1 * msi0_924[k]
                    - f_2 * msi1_924[k]
                    + f_3 * pc_x[k] * msk_1188[k];

        t_1486[k] = f_16 * lsk_936[k]
                    + f_3 * pc_y[k] * msk_1188[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pc_x, pc_y, pc_z, lsk_900, lsk_938, lsk_1191, \
                         msi0_927, msi1_927, msk_1188, msk_1190, \
                         msk_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_19 * lsk_900[k]
                    + f_3 * pc_z[k] * msk_1188[k];

        t_1488[k] = f_16 * lsk_1191[k]
                    + f_12 * msi0_927[k]
                    - f_13 * msi1_927[k]
                    + f_3 * pc_x[k] * msk_1191[k];

        t_1489[k] = f_16 * lsk_938[k]
                    + f_3 * pc_y[k] * msk_1190[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pc_x, pc_z, lsk_903, lsk_1193, lsk_1194, \
                         msi0_929, msi0_930, msi1_929, msi1_930, msk_1191, msk_1193, \
                         msk_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_16 * lsk_1193[k]
                    + f_12 * msi0_929[k]
                    - f_13 * msi1_929[k]
                    + f_3 * pc_x[k] * msk_1193[k];

        t_1491[k] = f_16 * lsk_1194[k]
                    + f_10 * msi0_930[k]
                    - f_11 * msi1_930[k]
                    + f_3 * pc_x[k] * msk_1194[k];

        t_1492[k] = f_19 * lsk_903[k]
                    + f_3 * pc_z[k] * msk_1191[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pc_x, pc_y, lsk_941, lsk_1197, lsk_1198, \
                         msi0_933, msi0_934, msi1_933, msi1_934, msk_1193, msk_1197, \
                         msk_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_16 * lsk_941[k]
                    + f_3 * pc_y[k] * msk_1193[k];

        t_1494[k] = f_16 * lsk_1197[k]
                    + f_10 * msi0_933[k]
                    - f_11 * msi1_933[k]
                    + f_3 * pc_x[k] * msk_1197[k];

        t_1495[k] = f_16 * lsk_1198[k]
                    + f_8 * msi0_934[k]
                    - f_9 * msi1_934[k]
                    + f_3 * pc_x[k] * msk_1198[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, pc_x, pc_y, pc_z, lsk_906, lsk_945, lsk_1200, \
                         msi0_936, msi1_936, msk_1194, msk_1197, \
                         msk_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_19 * lsk_906[k]
                    + f_3 * pc_z[k] * msk_1194[k];

        t_1497[k] = f_16 * lsk_1200[k]
                    + f_8 * msi0_936[k]
                    - f_9 * msi1_936[k]
                    + f_3 * pc_x[k] * msk_1200[k];

        t_1498[k] = f_16 * lsk_945[k]
                    + f_3 * pc_y[k] * msk_1197[k];
    }

#pragma omp simd aligned(t_1499, t_1500, t_1501, pc_x, pc_z, lsk_910, lsk_1202, lsk_1203, \
                         msi0_938, msi0_939, msi1_938, msi1_939, msk_1198, msk_1202, \
                         msk_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1499[k] = f_16 * lsk_1202[k]
                    + f_8 * msi0_938[k]
                    - f_9 * msi1_938[k]
                    + f_3 * pc_x[k] * msk_1202[k];

        t_1500[k] = f_16 * lsk_1203[k]
                    + f_6 * msi0_939[k]
                    - f_7 * msi1_939[k]
                    + f_3 * pc_x[k] * msk_1203[k];

        t_1501[k] = f_19 * lsk_910[k]
                    + f_3 * pc_z[k] * msk_1198[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, pc_x, pc_y, lsk_950, lsk_1205, lsk_1206, \
                         msi0_941, msi0_942, msi1_941, msi1_942, msk_1202, msk_1205, \
                         msk_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_16 * lsk_1205[k]
                    + f_6 * msi0_941[k]
                    - f_7 * msi1_941[k]
                    + f_3 * pc_x[k] * msk_1205[k];

        t_1503[k] = f_16 * lsk_1206[k]
                    + f_6 * msi0_942[k]
                    - f_7 * msi1_942[k]
                    + f_3 * pc_x[k] * msk_1206[k];

        t_1504[k] = f_16 * lsk_950[k]
                    + f_3 * pc_y[k] * msk_1202[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pc_x, pc_z, lsk_915, lsk_1208, lsk_1209, \
                         msi0_944, msi0_945, msi1_944, msi1_945, msk_1203, msk_1208, \
                         msk_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_16 * lsk_1208[k]
                    + f_6 * msi0_944[k]
                    - f_7 * msi1_944[k]
                    + f_3 * pc_x[k] * msk_1208[k];

        t_1506[k] = f_16 * lsk_1209[k]
                    + f_4 * msi0_945[k]
                    - f_5 * msi1_945[k]
                    + f_3 * pc_x[k] * msk_1209[k];

        t_1507[k] = f_19 * lsk_915[k]
                    + f_3 * pc_z[k] * msk_1203[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pc_x, lsk_1211, lsk_1212, lsk_1213, msi0_947, \
                         msi0_948, msi0_949, msi1_947, msi1_948, msi1_949, msk_1211, msk_1212, \
                         msk_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_16 * lsk_1211[k]
                    + f_4 * msi0_947[k]
                    - f_5 * msi1_947[k]
                    + f_3 * pc_x[k] * msk_1211[k];

        t_1509[k] = f_16 * lsk_1212[k]
                    + f_4 * msi0_948[k]
                    - f_5 * msi1_948[k]
                    + f_3 * pc_x[k] * msk_1212[k];

        t_1510[k] = f_16 * lsk_1213[k]
                    + f_4 * msi0_949[k]
                    - f_5 * msi1_949[k]
                    + f_3 * pc_x[k] * msk_1213[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pc_x, pc_y, lsk_956, lsk_1215, \
                         lsk_1216, lsk_1217, msi0_951, msi1_951, msk_1208, msk_1215, msk_1216, \
                         msk_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = f_16 * lsk_956[k]
                    + f_3 * pc_y[k] * msk_1208[k];

        t_1512[k] = f_16 * lsk_1215[k]
                    + f_4 * msi0_951[k]
                    - f_5 * msi1_951[k]
                    + f_3 * pc_x[k] * msk_1215[k];

        t_1513[k] = f_16 * lsk_1216[k]
                    + f_3 * pc_x[k] * msk_1216[k];

        t_1514[k] = f_16 * lsk_1217[k]
                    + f_3 * pc_x[k] * msk_1217[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, t_1518, t_1519, pc_x, lsk_1218, lsk_1219, \
                         lsk_1220, lsk_1221, lsk_1222, msk_1218, msk_1219, msk_1220, msk_1221, \
                         msk_1222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_16 * lsk_1218[k]
                    + f_3 * pc_x[k] * msk_1218[k];

        t_1516[k] = f_16 * lsk_1219[k]
                    + f_3 * pc_x[k] * msk_1219[k];

        t_1517[k] = f_16 * lsk_1220[k]
                    + f_3 * pc_x[k] * msk_1220[k];

        t_1518[k] = f_16 * lsk_1221[k]
                    + f_3 * pc_x[k] * msk_1221[k];

        t_1519[k] = f_16 * lsk_1222[k]
                    + f_3 * pc_x[k] * msk_1222[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pc_x, pc_y, pc_z, lsk_928, lsk_964, lsk_1223, \
                         msi0_945, msi1_945, msk_1216, msk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_16 * lsk_1223[k]
                    + f_3 * pc_x[k] * msk_1223[k];

        t_1521[k] = f_16 * lsk_964[k]
                    + f_1 * msi0_945[k]
                    - f_2 * msi1_945[k]
                    + f_3 * pc_y[k] * msk_1216[k];

        t_1522[k] = f_19 * lsk_928[k]
                    + f_3 * pc_z[k] * msk_1216[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_y, lsk_966, lsk_967, lsk_968, msi0_947, \
                         msi0_948, msi0_949, msi1_947, msi1_948, msi1_949, msk_1218, msk_1219, \
                         msk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_16 * lsk_966[k]
                    + f_12 * msi0_947[k]
                    - f_13 * msi1_947[k]
                    + f_3 * pc_y[k] * msk_1218[k];

        t_1524[k] = f_16 * lsk_967[k]
                    + f_10 * msi0_948[k]
                    - f_11 * msi1_948[k]
                    + f_3 * pc_y[k] * msk_1219[k];

        t_1525[k] = f_16 * lsk_968[k]
                    + f_8 * msi0_949[k]
                    - f_9 * msi1_949[k]
                    + f_3 * pc_y[k] * msk_1220[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pc_y, lsk_969, lsk_970, lsk_971, msi0_950, \
                         msi0_951, msi1_950, msi1_951, msk_1221, msk_1222, \
                         msk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_16 * lsk_969[k]
                    + f_6 * msi0_950[k]
                    - f_7 * msi1_950[k]
                    + f_3 * pc_y[k] * msk_1221[k];

        t_1527[k] = f_16 * lsk_970[k]
                    + f_4 * msi0_951[k]
                    - f_5 * msi1_951[k]
                    + f_3 * pc_y[k] * msk_1222[k];

        t_1528[k] = f_16 * lsk_971[k]
                    + f_3 * pc_y[k] * msk_1223[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, pa_y, pc_y, pc_z, lsl0_1215, lsk_935, \
                         lsk_936, lsk_972, lsl1_1215, msi0_951, msi1_951, msk_1223, \
                         msk_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_19 * lsk_935[k]
                    + f_1 * msi0_951[k]
                    - f_2 * msi1_951[k]
                    + f_3 * pc_z[k] * msk_1223[k];

        t_1530[k] = pa_y[k] * lsl0_1215[k]
                    - f_14 * pc_y[k] * lsl1_1215[k];

        t_1531[k] = f_15 * lsk_972[k]
                    + f_3 * pc_y[k] * msk_1224[k];

        t_1532[k] = f_20 * lsk_936[k]
                    + f_3 * pc_z[k] * msk_1224[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, t_1536, pa_y, pc_y, lsl0_1218, lsl0_1220, \
                         lsl0_1221, lsk_973, lsk_974, lsk_975, lsl1_1218, lsl1_1220, \
                         lsl1_1221, msk_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = pa_y[k] * lsl0_1218[k]
                    + f_16 * lsk_973[k]
                    - f_14 * pc_y[k] * lsl1_1218[k];

        t_1534[k] = f_15 * lsk_974[k]
                    + f_3 * pc_y[k] * msk_1226[k];

        t_1535[k] = pa_y[k] * lsl0_1220[k]
                    - f_14 * pc_y[k] * lsl1_1220[k];

        t_1536[k] = pa_y[k] * lsl0_1221[k]
                    + f_17 * lsk_975[k]
                    - f_14 * pc_y[k] * lsl1_1221[k];
    }

#pragma omp simd aligned(t_1537, t_1538, t_1539, t_1540, pa_y, pc_y, pc_z, lsl0_1224, \
                         lsl0_1225, lsk_939, lsk_977, lsk_978, lsl1_1224, lsl1_1225, msk_1227, \
                         msk_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1537[k] = f_20 * lsk_939[k]
                    + f_3 * pc_z[k] * msk_1227[k];

        t_1538[k] = f_15 * lsk_977[k]
                    + f_3 * pc_y[k] * msk_1229[k];

        t_1539[k] = pa_y[k] * lsl0_1224[k]
                    - f_14 * pc_y[k] * lsl1_1224[k];

        t_1540[k] = pa_y[k] * lsl0_1225[k]
                    + f_18 * lsk_978[k]
                    - f_14 * pc_y[k] * lsl1_1225[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, t_1544, pa_y, pc_y, pc_z, lsl0_1227, \
                         lsl0_1229, lsk_942, lsk_980, lsk_981, lsl1_1227, lsl1_1229, msk_1230, \
                         msk_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_20 * lsk_942[k]
                    + f_3 * pc_z[k] * msk_1230[k];

        t_1542[k] = pa_y[k] * lsl0_1227[k]
                    + f_16 * lsk_980[k]
                    - f_14 * pc_y[k] * lsl1_1227[k];

        t_1543[k] = f_15 * lsk_981[k]
                    + f_3 * pc_y[k] * msk_1233[k];

        t_1544[k] = pa_y[k] * lsl0_1229[k]
                    - f_14 * pc_y[k] * lsl1_1229[k];
    }

#pragma omp simd aligned(t_1545, t_1546, t_1547, pa_y, pc_y, pc_z, lsl0_1230, lsl0_1232, \
                         lsk_946, lsk_982, lsk_984, lsl1_1230, lsl1_1232, \
                         msk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1545[k] = pa_y[k] * lsl0_1230[k]
                    + f_19 * lsk_982[k]
                    - f_14 * pc_y[k] * lsl1_1230[k];

        t_1546[k] = f_20 * lsk_946[k]
                    + f_3 * pc_z[k] * msk_1234[k];

        t_1547[k] = pa_y[k] * lsl0_1232[k]
                    + f_17 * lsk_984[k]
                    - f_14 * pc_y[k] * lsl1_1232[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, t_1551, pa_y, pc_y, lsl0_1233, lsl0_1235, \
                         lsl0_1236, lsk_985, lsk_986, lsk_987, lsl1_1233, lsl1_1235, \
                         lsl1_1236, msk_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = pa_y[k] * lsl0_1233[k]
                    + f_16 * lsk_985[k]
                    - f_14 * pc_y[k] * lsl1_1233[k];

        t_1549[k] = f_15 * lsk_986[k]
                    + f_3 * pc_y[k] * msk_1238[k];

        t_1550[k] = pa_y[k] * lsl0_1235[k]
                    - f_14 * pc_y[k] * lsl1_1235[k];

        t_1551[k] = pa_y[k] * lsl0_1236[k]
                    + f_20 * lsk_987[k]
                    - f_14 * pc_y[k] * lsl1_1236[k];
    }

#pragma omp simd aligned(t_1552, t_1553, t_1554, pa_y, pc_y, pc_z, lsl0_1238, lsl0_1239, \
                         lsk_951, lsk_989, lsk_990, lsl1_1238, lsl1_1239, \
                         msk_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1552[k] = f_20 * lsk_951[k]
                    + f_3 * pc_z[k] * msk_1239[k];

        t_1553[k] = pa_y[k] * lsl0_1238[k]
                    + f_18 * lsk_989[k]
                    - f_14 * pc_y[k] * lsl1_1238[k];

        t_1554[k] = pa_y[k] * lsl0_1239[k]
                    + f_17 * lsk_990[k]
                    - f_14 * pc_y[k] * lsl1_1239[k];
    }

#pragma omp simd aligned(t_1555, t_1556, t_1557, t_1558, pa_y, pc_x, pc_y, lsl0_1240, \
                         lsl0_1242, lsk_991, lsk_992, lsk_1252, lsl1_1240, lsl1_1242, \
                         msk_1244, msk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1555[k] = pa_y[k] * lsl0_1240[k]
                    + f_16 * lsk_991[k]
                    - f_14 * pc_y[k] * lsl1_1240[k];

        t_1556[k] = f_15 * lsk_992[k]
                    + f_3 * pc_y[k] * msk_1244[k];

        t_1557[k] = pa_y[k] * lsl0_1242[k]
                    - f_14 * pc_y[k] * lsl1_1242[k];

        t_1558[k] = f_16 * lsk_1252[k]
                    + f_3 * pc_x[k] * msk_1252[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, t_1562, t_1563, pc_x, lsk_1253, lsk_1254, \
                         lsk_1255, lsk_1256, lsk_1257, msk_1253, msk_1254, msk_1255, msk_1256, \
                         msk_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_16 * lsk_1253[k]
                    + f_3 * pc_x[k] * msk_1253[k];

        t_1560[k] = f_16 * lsk_1254[k]
                    + f_3 * pc_x[k] * msk_1254[k];

        t_1561[k] = f_16 * lsk_1255[k]
                    + f_3 * pc_x[k] * msk_1255[k];

        t_1562[k] = f_16 * lsk_1256[k]
                    + f_3 * pc_x[k] * msk_1256[k];

        t_1563[k] = f_16 * lsk_1257[k]
                    + f_3 * pc_x[k] * msk_1257[k];
    }

#pragma omp simd aligned(t_1564, t_1565, t_1566, t_1567, pc_x, pc_y, pc_z, lsk_964, lsk_1000, \
                         lsk_1258, lsk_1259, msi0_973, msi1_973, msk_1252, msk_1258, \
                         msk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1564[k] = f_16 * lsk_1258[k]
                    + f_3 * pc_x[k] * msk_1258[k];

        t_1565[k] = f_16 * lsk_1259[k]
                    + f_3 * pc_x[k] * msk_1259[k];

        t_1566[k] = f_15 * lsk_1000[k]
                    + f_1 * msi0_973[k]
                    - f_2 * msi1_973[k]
                    + f_3 * pc_y[k] * msk_1252[k];

        t_1567[k] = f_20 * lsk_964[k]
                    + f_3 * pc_z[k] * msk_1252[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, pc_y, lsk_1002, lsk_1003, lsk_1004, msi0_975, \
                         msi0_976, msi0_977, msi1_975, msi1_976, msi1_977, msk_1254, msk_1255, \
                         msk_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = f_15 * lsk_1002[k]
                    + f_12 * msi0_975[k]
                    - f_13 * msi1_975[k]
                    + f_3 * pc_y[k] * msk_1254[k];

        t_1569[k] = f_15 * lsk_1003[k]
                    + f_10 * msi0_976[k]
                    - f_11 * msi1_976[k]
                    + f_3 * pc_y[k] * msk_1255[k];

        t_1570[k] = f_15 * lsk_1004[k]
                    + f_8 * msi0_977[k]
                    - f_9 * msi1_977[k]
                    + f_3 * pc_y[k] * msk_1256[k];
    }

#pragma omp simd aligned(t_1571, t_1572, t_1573, pc_y, lsk_1005, lsk_1006, lsk_1007, msi0_978, \
                         msi0_979, msi1_978, msi1_979, msk_1257, msk_1258, \
                         msk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1571[k] = f_15 * lsk_1005[k]
                    + f_6 * msi0_978[k]
                    - f_7 * msi1_978[k]
                    + f_3 * pc_y[k] * msk_1257[k];

        t_1572[k] = f_15 * lsk_1006[k]
                    + f_4 * msi0_979[k]
                    - f_5 * msi1_979[k]
                    + f_3 * pc_y[k] * msk_1258[k];

        t_1573[k] = f_15 * lsk_1007[k]
                    + f_3 * pc_y[k] * msk_1259[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);
    const auto f_24 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1259 = buffer.data(lsl0 + 1259);
    const auto *lsl0_1260 = buffer.data(lsl0 + 1260);
    const auto *lsl0_1263 = buffer.data(lsl0 + 1263);
    const auto *lsl0_1266 = buffer.data(lsl0 + 1266);
    const auto *lsl0_1270 = buffer.data(lsl0 + 1270);
    const auto *lsl0_1275 = buffer.data(lsl0 + 1275);
    const auto *lsl0_1281 = buffer.data(lsl0 + 1281);
    const auto *lsl0_1620 = buffer.data(lsl0 + 1620);
    const auto *lsl0_1623 = buffer.data(lsl0 + 1623);
    const auto *lsl0_1626 = buffer.data(lsl0 + 1626);
    const auto *lsl0_1630 = buffer.data(lsl0 + 1630);
    const auto *lsl0_1635 = buffer.data(lsl0 + 1635);
    const auto *lsl0_1641 = buffer.data(lsl0 + 1641);
    const auto *lsl0_1656 = buffer.data(lsl0 + 1656);
    const auto *lsl0_1658 = buffer.data(lsl0 + 1658);
    const auto *lsl0_1659 = buffer.data(lsl0 + 1659);
    const auto *lsl0_1660 = buffer.data(lsl0 + 1660);
    const auto *lsl0_1661 = buffer.data(lsl0 + 1661);
    const auto *lsl0_1662 = buffer.data(lsl0 + 1662);
    const auto *lsl0_1664 = buffer.data(lsl0 + 1664);
    const auto *lsl0_1670 = buffer.data(lsl0 + 1670);
    const auto *lsl0_1674 = buffer.data(lsl0 + 1674);
    const auto *lsl0_1677 = buffer.data(lsl0 + 1677);
    const auto *lsl0_1679 = buffer.data(lsl0 + 1679);
    const auto *lsl0_1682 = buffer.data(lsl0 + 1682);
    const auto *lsl0_1683 = buffer.data(lsl0 + 1683);
    const auto *lsl0_1685 = buffer.data(lsl0 + 1685);
    const auto *lsl0_1688 = buffer.data(lsl0 + 1688);
    const auto *lsl0_1689 = buffer.data(lsl0 + 1689);
    const auto *lsl0_1690 = buffer.data(lsl0 + 1690);
    const auto *lsl0_1692 = buffer.data(lsl0 + 1692);

    const auto *lsk_972 = buffer.data(lsk + 972);
    const auto *lsk_1007 = buffer.data(lsk + 1007);
    const auto *lsk_1008 = buffer.data(lsk + 1008);
    const auto *lsk_1011 = buffer.data(lsk + 1011);
    const auto *lsk_1013 = buffer.data(lsk + 1013);
    const auto *lsk_1014 = buffer.data(lsk + 1014);
    const auto *lsk_1017 = buffer.data(lsk + 1017);
    const auto *lsk_1018 = buffer.data(lsk + 1018);
    const auto *lsk_1022 = buffer.data(lsk + 1022);
    const auto *lsk_1023 = buffer.data(lsk + 1023);
    const auto *lsk_1028 = buffer.data(lsk + 1028);
    const auto *lsk_1043 = buffer.data(lsk + 1043);
    const auto *lsk_1044 = buffer.data(lsk + 1044);
    const auto *lsk_1046 = buffer.data(lsk + 1046);
    const auto *lsk_1049 = buffer.data(lsk + 1049);
    const auto *lsk_1053 = buffer.data(lsk + 1053);
    const auto *lsk_1058 = buffer.data(lsk + 1058);
    const auto *lsk_1064 = buffer.data(lsk + 1064);
    const auto *lsk_1260 = buffer.data(lsk + 1260);
    const auto *lsk_1265 = buffer.data(lsk + 1265);
    const auto *lsk_1269 = buffer.data(lsk + 1269);
    const auto *lsk_1274 = buffer.data(lsk + 1274);
    const auto *lsk_1280 = buffer.data(lsk + 1280);
    const auto *lsk_1287 = buffer.data(lsk + 1287);
    const auto *lsk_1288 = buffer.data(lsk + 1288);
    const auto *lsk_1289 = buffer.data(lsk + 1289);
    const auto *lsk_1290 = buffer.data(lsk + 1290);
    const auto *lsk_1291 = buffer.data(lsk + 1291);
    const auto *lsk_1292 = buffer.data(lsk + 1292);
    const auto *lsk_1293 = buffer.data(lsk + 1293);
    const auto *lsk_1295 = buffer.data(lsk + 1295);
    const auto *lsk_1296 = buffer.data(lsk + 1296);
    const auto *lsk_1299 = buffer.data(lsk + 1299);
    const auto *lsk_1302 = buffer.data(lsk + 1302);
    const auto *lsk_1306 = buffer.data(lsk + 1306);
    const auto *lsk_1311 = buffer.data(lsk + 1311);
    const auto *lsk_1317 = buffer.data(lsk + 1317);
    const auto *lsk_1324 = buffer.data(lsk + 1324);
    const auto *lsk_1326 = buffer.data(lsk + 1326);
    const auto *lsk_1327 = buffer.data(lsk + 1327);
    const auto *lsk_1328 = buffer.data(lsk + 1328);
    const auto *lsk_1329 = buffer.data(lsk + 1329);
    const auto *lsk_1330 = buffer.data(lsk + 1330);
    const auto *lsk_1331 = buffer.data(lsk + 1331);
    const auto *lsk_1337 = buffer.data(lsk + 1337);
    const auto *lsk_1341 = buffer.data(lsk + 1341);
    const auto *lsk_1344 = buffer.data(lsk + 1344);
    const auto *lsk_1346 = buffer.data(lsk + 1346);
    const auto *lsk_1349 = buffer.data(lsk + 1349);
    const auto *lsk_1350 = buffer.data(lsk + 1350);
    const auto *lsk_1352 = buffer.data(lsk + 1352);
    const auto *lsk_1355 = buffer.data(lsk + 1355);
    const auto *lsk_1356 = buffer.data(lsk + 1356);
    const auto *lsk_1357 = buffer.data(lsk + 1357);
    const auto *lsk_1359 = buffer.data(lsk + 1359);

    const auto *lsl1_1259 = buffer.data(lsl1 + 1259);
    const auto *lsl1_1260 = buffer.data(lsl1 + 1260);
    const auto *lsl1_1263 = buffer.data(lsl1 + 1263);
    const auto *lsl1_1266 = buffer.data(lsl1 + 1266);
    const auto *lsl1_1270 = buffer.data(lsl1 + 1270);
    const auto *lsl1_1275 = buffer.data(lsl1 + 1275);
    const auto *lsl1_1281 = buffer.data(lsl1 + 1281);
    const auto *lsl1_1620 = buffer.data(lsl1 + 1620);
    const auto *lsl1_1623 = buffer.data(lsl1 + 1623);
    const auto *lsl1_1626 = buffer.data(lsl1 + 1626);
    const auto *lsl1_1630 = buffer.data(lsl1 + 1630);
    const auto *lsl1_1635 = buffer.data(lsl1 + 1635);
    const auto *lsl1_1641 = buffer.data(lsl1 + 1641);
    const auto *lsl1_1656 = buffer.data(lsl1 + 1656);
    const auto *lsl1_1658 = buffer.data(lsl1 + 1658);
    const auto *lsl1_1659 = buffer.data(lsl1 + 1659);
    const auto *lsl1_1660 = buffer.data(lsl1 + 1660);
    const auto *lsl1_1661 = buffer.data(lsl1 + 1661);
    const auto *lsl1_1662 = buffer.data(lsl1 + 1662);
    const auto *lsl1_1664 = buffer.data(lsl1 + 1664);
    const auto *lsl1_1670 = buffer.data(lsl1 + 1670);
    const auto *lsl1_1674 = buffer.data(lsl1 + 1674);
    const auto *lsl1_1677 = buffer.data(lsl1 + 1677);
    const auto *lsl1_1679 = buffer.data(lsl1 + 1679);
    const auto *lsl1_1682 = buffer.data(lsl1 + 1682);
    const auto *lsl1_1683 = buffer.data(lsl1 + 1683);
    const auto *lsl1_1685 = buffer.data(lsl1 + 1685);
    const auto *lsl1_1688 = buffer.data(lsl1 + 1688);
    const auto *lsl1_1689 = buffer.data(lsl1 + 1689);
    const auto *lsl1_1690 = buffer.data(lsl1 + 1690);
    const auto *lsl1_1692 = buffer.data(lsl1 + 1692);

    const auto *msi0_980 = buffer.data(msi0 + 980);
    const auto *msi0_981 = buffer.data(msi0 + 981);
    const auto *msi0_982 = buffer.data(msi0 + 982);
    const auto *msi0_983 = buffer.data(msi0 + 983);
    const auto *msi0_984 = buffer.data(msi0 + 984);
    const auto *msi0_985 = buffer.data(msi0 + 985);
    const auto *msi0_986 = buffer.data(msi0 + 986);
    const auto *msi0_987 = buffer.data(msi0 + 987);
    const auto *msi0_988 = buffer.data(msi0 + 988);
    const auto *msi0_989 = buffer.data(msi0 + 989);
    const auto *msi0_990 = buffer.data(msi0 + 990);
    const auto *msi0_991 = buffer.data(msi0 + 991);
    const auto *msi0_992 = buffer.data(msi0 + 992);
    const auto *msi0_993 = buffer.data(msi0 + 993);
    const auto *msi0_994 = buffer.data(msi0 + 994);
    const auto *msi0_1000 = buffer.data(msi0 + 1000);
    const auto *msi0_1001 = buffer.data(msi0 + 1001);
    const auto *msi0_1002 = buffer.data(msi0 + 1002);
    const auto *msi0_1003 = buffer.data(msi0 + 1003);
    const auto *msi0_1004 = buffer.data(msi0 + 1004);
    const auto *msi0_1005 = buffer.data(msi0 + 1005);
    const auto *msi0_1006 = buffer.data(msi0 + 1006);
    const auto *msi0_1007 = buffer.data(msi0 + 1007);
    const auto *msi0_1008 = buffer.data(msi0 + 1008);
    const auto *msi0_1010 = buffer.data(msi0 + 1010);
    const auto *msi0_1011 = buffer.data(msi0 + 1011);
    const auto *msi0_1013 = buffer.data(msi0 + 1013);
    const auto *msi0_1014 = buffer.data(msi0 + 1014);
    const auto *msi0_1015 = buffer.data(msi0 + 1015);
    const auto *msi0_1017 = buffer.data(msi0 + 1017);
    const auto *msi0_1018 = buffer.data(msi0 + 1018);
    const auto *msi0_1019 = buffer.data(msi0 + 1019);
    const auto *msi0_1020 = buffer.data(msi0 + 1020);
    const auto *msi0_1022 = buffer.data(msi0 + 1022);

    const auto *msi1_980 = buffer.data(msi1 + 980);
    const auto *msi1_981 = buffer.data(msi1 + 981);
    const auto *msi1_982 = buffer.data(msi1 + 982);
    const auto *msi1_983 = buffer.data(msi1 + 983);
    const auto *msi1_984 = buffer.data(msi1 + 984);
    const auto *msi1_985 = buffer.data(msi1 + 985);
    const auto *msi1_986 = buffer.data(msi1 + 986);
    const auto *msi1_987 = buffer.data(msi1 + 987);
    const auto *msi1_988 = buffer.data(msi1 + 988);
    const auto *msi1_989 = buffer.data(msi1 + 989);
    const auto *msi1_990 = buffer.data(msi1 + 990);
    const auto *msi1_991 = buffer.data(msi1 + 991);
    const auto *msi1_992 = buffer.data(msi1 + 992);
    const auto *msi1_993 = buffer.data(msi1 + 993);
    const auto *msi1_994 = buffer.data(msi1 + 994);
    const auto *msi1_1000 = buffer.data(msi1 + 1000);
    const auto *msi1_1001 = buffer.data(msi1 + 1001);
    const auto *msi1_1002 = buffer.data(msi1 + 1002);
    const auto *msi1_1003 = buffer.data(msi1 + 1003);
    const auto *msi1_1004 = buffer.data(msi1 + 1004);
    const auto *msi1_1005 = buffer.data(msi1 + 1005);
    const auto *msi1_1006 = buffer.data(msi1 + 1006);
    const auto *msi1_1007 = buffer.data(msi1 + 1007);
    const auto *msi1_1008 = buffer.data(msi1 + 1008);
    const auto *msi1_1010 = buffer.data(msi1 + 1010);
    const auto *msi1_1011 = buffer.data(msi1 + 1011);
    const auto *msi1_1013 = buffer.data(msi1 + 1013);
    const auto *msi1_1014 = buffer.data(msi1 + 1014);
    const auto *msi1_1015 = buffer.data(msi1 + 1015);
    const auto *msi1_1017 = buffer.data(msi1 + 1017);
    const auto *msi1_1018 = buffer.data(msi1 + 1018);
    const auto *msi1_1019 = buffer.data(msi1 + 1019);
    const auto *msi1_1020 = buffer.data(msi1 + 1020);
    const auto *msi1_1022 = buffer.data(msi1 + 1022);

    const auto *msk_1260 = buffer.data(msk + 1260);
    const auto *msk_1261 = buffer.data(msk + 1261);
    const auto *msk_1262 = buffer.data(msk + 1262);
    const auto *msk_1263 = buffer.data(msk + 1263);
    const auto *msk_1264 = buffer.data(msk + 1264);
    const auto *msk_1265 = buffer.data(msk + 1265);
    const auto *msk_1266 = buffer.data(msk + 1266);
    const auto *msk_1267 = buffer.data(msk + 1267);
    const auto *msk_1268 = buffer.data(msk + 1268);
    const auto *msk_1269 = buffer.data(msk + 1269);
    const auto *msk_1270 = buffer.data(msk + 1270);
    const auto *msk_1271 = buffer.data(msk + 1271);
    const auto *msk_1272 = buffer.data(msk + 1272);
    const auto *msk_1273 = buffer.data(msk + 1273);
    const auto *msk_1274 = buffer.data(msk + 1274);
    const auto *msk_1275 = buffer.data(msk + 1275);
    const auto *msk_1276 = buffer.data(msk + 1276);
    const auto *msk_1277 = buffer.data(msk + 1277);
    const auto *msk_1278 = buffer.data(msk + 1278);
    const auto *msk_1279 = buffer.data(msk + 1279);
    const auto *msk_1280 = buffer.data(msk + 1280);
    const auto *msk_1287 = buffer.data(msk + 1287);
    const auto *msk_1288 = buffer.data(msk + 1288);
    const auto *msk_1289 = buffer.data(msk + 1289);
    const auto *msk_1290 = buffer.data(msk + 1290);
    const auto *msk_1291 = buffer.data(msk + 1291);
    const auto *msk_1292 = buffer.data(msk + 1292);
    const auto *msk_1293 = buffer.data(msk + 1293);
    const auto *msk_1294 = buffer.data(msk + 1294);
    const auto *msk_1295 = buffer.data(msk + 1295);
    const auto *msk_1296 = buffer.data(msk + 1296);
    const auto *msk_1297 = buffer.data(msk + 1297);
    const auto *msk_1298 = buffer.data(msk + 1298);
    const auto *msk_1299 = buffer.data(msk + 1299);
    const auto *msk_1301 = buffer.data(msk + 1301);
    const auto *msk_1302 = buffer.data(msk + 1302);
    const auto *msk_1303 = buffer.data(msk + 1303);
    const auto *msk_1305 = buffer.data(msk + 1305);
    const auto *msk_1306 = buffer.data(msk + 1306);
    const auto *msk_1307 = buffer.data(msk + 1307);
    const auto *msk_1308 = buffer.data(msk + 1308);
    const auto *msk_1310 = buffer.data(msk + 1310);
    const auto *msk_1311 = buffer.data(msk + 1311);
    const auto *msk_1312 = buffer.data(msk + 1312);
    const auto *msk_1313 = buffer.data(msk + 1313);
    const auto *msk_1314 = buffer.data(msk + 1314);
    const auto *msk_1316 = buffer.data(msk + 1316);
    const auto *msk_1317 = buffer.data(msk + 1317);
    const auto *msk_1324 = buffer.data(msk + 1324);
    const auto *msk_1326 = buffer.data(msk + 1326);
    const auto *msk_1327 = buffer.data(msk + 1327);
    const auto *msk_1328 = buffer.data(msk + 1328);
    const auto *msk_1329 = buffer.data(msk + 1329);
    const auto *msk_1330 = buffer.data(msk + 1330);
    const auto *msk_1331 = buffer.data(msk + 1331);
    const auto *msk_1332 = buffer.data(msk + 1332);
    const auto *msk_1334 = buffer.data(msk + 1334);
    const auto *msk_1335 = buffer.data(msk + 1335);
    const auto *msk_1337 = buffer.data(msk + 1337);
    const auto *msk_1338 = buffer.data(msk + 1338);
    const auto *msk_1341 = buffer.data(msk + 1341);
    const auto *msk_1342 = buffer.data(msk + 1342);
    const auto *msk_1346 = buffer.data(msk + 1346);
    const auto *msk_1347 = buffer.data(msk + 1347);
    const auto *msk_1352 = buffer.data(msk + 1352);

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pa_y, pc_x, pc_y, pc_z, lsl0_1259, \
                         lsk_972, lsk_1260, lsl1_1259, msi0_980, msi1_980, \
                         msk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = pa_y[k] * lsl0_1259[k]
                    - f_14 * pc_y[k] * lsl1_1259[k];

        t_1575[k] = f_16 * lsk_1260[k]
                    + f_1 * msi0_980[k]
                    - f_2 * msi1_980[k]
                    + f_3 * pc_x[k] * msk_1260[k];

        t_1576[k] = f_3 * pc_y[k] * msk_1260[k];

        t_1577[k] = f_24 * lsk_972[k]
                    + f_3 * pc_z[k] * msk_1260[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pc_x, pc_y, lsk_1265, msi0_980, msi0_985, \
                         msi1_980, msi1_985, msk_1261, msk_1262, \
                         msk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = f_4 * msi0_980[k]
                    - f_5 * msi1_980[k]
                    + f_3 * pc_y[k] * msk_1261[k];

        t_1579[k] = f_3 * pc_y[k] * msk_1262[k];

        t_1580[k] = f_16 * lsk_1265[k]
                    + f_12 * msi0_985[k]
                    - f_13 * msi1_985[k]
                    + f_3 * pc_x[k] * msk_1265[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pc_y, msi0_981, msi0_982, msi1_981, msi1_982, \
                         msk_1263, msk_1264, msk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_6 * msi0_981[k]
                    - f_7 * msi1_981[k]
                    + f_3 * pc_y[k] * msk_1263[k];

        t_1582[k] = f_4 * msi0_982[k]
                    - f_5 * msi1_982[k]
                    + f_3 * pc_y[k] * msk_1264[k];

        t_1583[k] = f_3 * pc_y[k] * msk_1265[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, pc_x, pc_y, lsk_1269, msi0_983, msi0_984, \
                         msi0_989, msi1_983, msi1_984, msi1_989, msk_1266, msk_1267, \
                         msk_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_16 * lsk_1269[k]
                    + f_10 * msi0_989[k]
                    - f_11 * msi1_989[k]
                    + f_3 * pc_x[k] * msk_1269[k];

        t_1585[k] = f_8 * msi0_983[k]
                    - f_9 * msi1_983[k]
                    + f_3 * pc_y[k] * msk_1266[k];

        t_1586[k] = f_6 * msi0_984[k]
                    - f_7 * msi1_984[k]
                    + f_3 * pc_y[k] * msk_1267[k];
    }

#pragma omp simd aligned(t_1587, t_1588, t_1589, pc_x, pc_y, lsk_1274, msi0_985, msi0_994, \
                         msi1_985, msi1_994, msk_1268, msk_1269, \
                         msk_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1587[k] = f_4 * msi0_985[k]
                    - f_5 * msi1_985[k]
                    + f_3 * pc_y[k] * msk_1268[k];

        t_1588[k] = f_3 * pc_y[k] * msk_1269[k];

        t_1589[k] = f_16 * lsk_1274[k]
                    + f_8 * msi0_994[k]
                    - f_9 * msi1_994[k]
                    + f_3 * pc_x[k] * msk_1274[k];
    }

#pragma omp simd aligned(t_1590, t_1591, t_1592, pc_y, msi0_986, msi0_987, msi0_988, msi1_986, \
                         msi1_987, msi1_988, msk_1270, msk_1271, \
                         msk_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1590[k] = f_10 * msi0_986[k]
                    - f_11 * msi1_986[k]
                    + f_3 * pc_y[k] * msk_1270[k];

        t_1591[k] = f_8 * msi0_987[k]
                    - f_9 * msi1_987[k]
                    + f_3 * pc_y[k] * msk_1271[k];

        t_1592[k] = f_6 * msi0_988[k]
                    - f_7 * msi1_988[k]
                    + f_3 * pc_y[k] * msk_1272[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pc_x, pc_y, lsk_1280, msi0_989, msi0_1000, \
                         msi1_989, msi1_1000, msk_1273, msk_1274, \
                         msk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_4 * msi0_989[k]
                    - f_5 * msi1_989[k]
                    + f_3 * pc_y[k] * msk_1273[k];

        t_1594[k] = f_3 * pc_y[k] * msk_1274[k];

        t_1595[k] = f_16 * lsk_1280[k]
                    + f_6 * msi0_1000[k]
                    - f_7 * msi1_1000[k]
                    + f_3 * pc_x[k] * msk_1280[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, pc_y, msi0_990, msi0_991, msi0_992, msi1_990, \
                         msi1_991, msi1_992, msk_1275, msk_1276, \
                         msk_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_12 * msi0_990[k]
                    - f_13 * msi1_990[k]
                    + f_3 * pc_y[k] * msk_1275[k];

        t_1597[k] = f_10 * msi0_991[k]
                    - f_11 * msi1_991[k]
                    + f_3 * pc_y[k] * msk_1276[k];

        t_1598[k] = f_8 * msi0_992[k]
                    - f_9 * msi1_992[k]
                    + f_3 * pc_y[k] * msk_1277[k];
    }

#pragma omp simd aligned(t_1599, t_1600, t_1601, pc_y, msi0_993, msi0_994, msi1_993, msi1_994, \
                         msk_1278, msk_1279, msk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = f_6 * msi0_993[k]
                    - f_7 * msi1_993[k]
                    + f_3 * pc_y[k] * msk_1278[k];

        t_1600[k] = f_4 * msi0_994[k]
                    - f_5 * msi1_994[k]
                    + f_3 * pc_y[k] * msk_1279[k];

        t_1601[k] = f_3 * pc_y[k] * msk_1280[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, t_1605, pc_x, lsk_1287, lsk_1288, lsk_1289, \
                         lsk_1290, msi0_1007, msi1_1007, msk_1287, msk_1288, msk_1289, \
                         msk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = f_16 * lsk_1287[k]
                    + f_4 * msi0_1007[k]
                    - f_5 * msi1_1007[k]
                    + f_3 * pc_x[k] * msk_1287[k];

        t_1603[k] = f_16 * lsk_1288[k]
                    + f_3 * pc_x[k] * msk_1288[k];

        t_1604[k] = f_16 * lsk_1289[k]
                    + f_3 * pc_x[k] * msk_1289[k];

        t_1605[k] = f_16 * lsk_1290[k]
                    + f_3 * pc_x[k] * msk_1290[k];
    }

#pragma omp simd aligned(t_1606, t_1607, t_1608, t_1609, t_1610, pc_x, pc_y, lsk_1291, \
                         lsk_1292, lsk_1293, lsk_1295, msk_1287, msk_1291, msk_1292, msk_1293, \
                         msk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1606[k] = f_16 * lsk_1291[k]
                    + f_3 * pc_x[k] * msk_1291[k];

        t_1607[k] = f_16 * lsk_1292[k]
                    + f_3 * pc_x[k] * msk_1292[k];

        t_1608[k] = f_16 * lsk_1293[k]
                    + f_3 * pc_x[k] * msk_1293[k];

        t_1609[k] = f_3 * pc_y[k] * msk_1287[k];

        t_1610[k] = f_16 * lsk_1295[k]
                    + f_3 * pc_x[k] * msk_1295[k];
    }

#pragma omp simd aligned(t_1611, t_1612, t_1613, pc_y, msi0_1001, msi0_1002, msi0_1003, \
                         msi1_1001, msi1_1002, msi1_1003, msk_1288, msk_1289, \
                         msk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1611[k] = f_1 * msi0_1001[k]
                    - f_2 * msi1_1001[k]
                    + f_3 * pc_y[k] * msk_1288[k];

        t_1612[k] = f_22 * msi0_1002[k]
                    - f_23 * msi1_1002[k]
                    + f_3 * pc_y[k] * msk_1289[k];

        t_1613[k] = f_12 * msi0_1003[k]
                    - f_13 * msi1_1003[k]
                    + f_3 * pc_y[k] * msk_1290[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, pc_y, msi0_1004, msi0_1005, msi0_1006, \
                         msi1_1004, msi1_1005, msi1_1006, msk_1291, msk_1292, \
                         msk_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_10 * msi0_1004[k]
                    - f_11 * msi1_1004[k]
                    + f_3 * pc_y[k] * msk_1291[k];

        t_1615[k] = f_8 * msi0_1005[k]
                    - f_9 * msi1_1005[k]
                    + f_3 * pc_y[k] * msk_1292[k];

        t_1616[k] = f_6 * msi0_1006[k]
                    - f_7 * msi1_1006[k]
                    + f_3 * pc_y[k] * msk_1293[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, t_1620, pa_x, pc_x, pc_y, pc_z, lsl0_1620, \
                         lsk_1007, lsk_1296, lsl1_1620, msi0_1007, msi1_1007, msk_1294, \
                         msk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_4 * msi0_1007[k]
                    - f_5 * msi1_1007[k]
                    + f_3 * pc_y[k] * msk_1294[k];

        t_1618[k] = f_3 * pc_y[k] * msk_1295[k];

        t_1619[k] = f_24 * lsk_1007[k]
                    + f_1 * msi0_1007[k]
                    - f_2 * msi1_1007[k]
                    + f_3 * pc_z[k] * msk_1295[k];

        t_1620[k] = pa_x[k] * lsl0_1620[k]
                    + f_21 * lsk_1296[k]
                    - f_14 * pc_x[k] * lsl1_1620[k];
    }

#pragma omp simd aligned(t_1621, t_1622, t_1623, t_1624, pa_x, pc_x, pc_y, pc_z, lsl0_1623, \
                         lsk_1008, lsk_1299, lsl1_1623, msk_1296, \
                         msk_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1621[k] = f_21 * lsk_1008[k]
                    + f_3 * pc_y[k] * msk_1296[k];

        t_1622[k] = f_3 * pc_z[k] * msk_1296[k];

        t_1623[k] = pa_x[k] * lsl0_1623[k]
                    + f_20 * lsk_1299[k]
                    - f_14 * pc_x[k] * lsl1_1623[k];

        t_1624[k] = f_3 * pc_z[k] * msk_1297[k];
    }

#pragma omp simd aligned(t_1625, t_1626, t_1627, pa_x, pc_x, pc_z, lsl0_1626, lsk_1302, \
                         lsl1_1626, msi0_1008, msi1_1008, msk_1298, \
                         msk_1299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = f_4 * msi0_1008[k]
                    - f_5 * msi1_1008[k]
                    + f_3 * pc_z[k] * msk_1298[k];

        t_1626[k] = pa_x[k] * lsl0_1626[k]
                    + f_19 * lsk_1302[k]
                    - f_14 * pc_x[k] * lsl1_1626[k];

        t_1627[k] = f_3 * pc_z[k] * msk_1299[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, t_1631, pa_x, pc_x, pc_y, pc_z, lsl0_1630, \
                         lsk_1013, lsk_1306, lsl1_1630, msi0_1010, msi1_1010, msk_1301, \
                         msk_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_21 * lsk_1013[k]
                    + f_3 * pc_y[k] * msk_1301[k];

        t_1629[k] = f_6 * msi0_1010[k]
                    - f_7 * msi1_1010[k]
                    + f_3 * pc_z[k] * msk_1301[k];

        t_1630[k] = pa_x[k] * lsl0_1630[k]
                    + f_18 * lsk_1306[k]
                    - f_14 * pc_x[k] * lsl1_1630[k];

        t_1631[k] = f_3 * pc_z[k] * msk_1302[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pc_y, pc_z, lsk_1017, msi0_1011, msi0_1013, \
                         msi1_1011, msi1_1013, msk_1303, msk_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_4 * msi0_1011[k]
                    - f_5 * msi1_1011[k]
                    + f_3 * pc_z[k] * msk_1303[k];

        t_1633[k] = f_21 * lsk_1017[k]
                    + f_3 * pc_y[k] * msk_1305[k];

        t_1634[k] = f_8 * msi0_1013[k]
                    - f_9 * msi1_1013[k]
                    + f_3 * pc_z[k] * msk_1305[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, pa_x, pc_x, pc_z, lsl0_1635, lsk_1311, \
                         lsl1_1635, msi0_1014, msi1_1014, msk_1306, \
                         msk_1307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = pa_x[k] * lsl0_1635[k]
                    + f_17 * lsk_1311[k]
                    - f_14 * pc_x[k] * lsl1_1635[k];

        t_1636[k] = f_3 * pc_z[k] * msk_1306[k];

        t_1637[k] = f_4 * msi0_1014[k]
                    - f_5 * msi1_1014[k]
                    + f_3 * pc_z[k] * msk_1307[k];
    }

#pragma omp simd aligned(t_1638, t_1639, t_1640, pc_y, pc_z, lsk_1022, msi0_1015, msi0_1017, \
                         msi1_1015, msi1_1017, msk_1308, msk_1310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1638[k] = f_6 * msi0_1015[k]
                    - f_7 * msi1_1015[k]
                    + f_3 * pc_z[k] * msk_1308[k];

        t_1639[k] = f_21 * lsk_1022[k]
                    + f_3 * pc_y[k] * msk_1310[k];

        t_1640[k] = f_10 * msi0_1017[k]
                    - f_11 * msi1_1017[k]
                    + f_3 * pc_z[k] * msk_1310[k];
    }

#pragma omp simd aligned(t_1641, t_1642, t_1643, pa_x, pc_x, pc_z, lsl0_1641, lsk_1317, \
                         lsl1_1641, msi0_1018, msi1_1018, msk_1311, \
                         msk_1312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1641[k] = pa_x[k] * lsl0_1641[k]
                    + f_16 * lsk_1317[k]
                    - f_14 * pc_x[k] * lsl1_1641[k];

        t_1642[k] = f_3 * pc_z[k] * msk_1311[k];

        t_1643[k] = f_4 * msi0_1018[k]
                    - f_5 * msi1_1018[k]
                    + f_3 * pc_z[k] * msk_1312[k];
    }

#pragma omp simd aligned(t_1644, t_1645, t_1646, t_1647, pc_y, pc_z, lsk_1028, msi0_1019, \
                         msi0_1020, msi0_1022, msi1_1019, msi1_1020, msi1_1022, msk_1313, \
                         msk_1314, msk_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1644[k] = f_6 * msi0_1019[k]
                    - f_7 * msi1_1019[k]
                    + f_3 * pc_z[k] * msk_1313[k];

        t_1645[k] = f_8 * msi0_1020[k]
                    - f_9 * msi1_1020[k]
                    + f_3 * pc_z[k] * msk_1314[k];

        t_1646[k] = f_21 * lsk_1028[k]
                    + f_3 * pc_y[k] * msk_1316[k];

        t_1647[k] = f_12 * msi0_1022[k]
                    - f_13 * msi1_1022[k]
                    + f_3 * pc_z[k] * msk_1316[k];
    }

#pragma omp simd aligned(t_1648, t_1649, t_1650, t_1651, t_1652, pc_x, pc_z, lsk_1324, \
                         lsk_1326, lsk_1327, lsk_1328, msk_1317, msk_1324, msk_1326, msk_1327, \
                         msk_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1648[k] = f_15 * lsk_1324[k]
                    + f_3 * pc_x[k] * msk_1324[k];

        t_1649[k] = f_3 * pc_z[k] * msk_1317[k];

        t_1650[k] = f_15 * lsk_1326[k]
                    + f_3 * pc_x[k] * msk_1326[k];

        t_1651[k] = f_15 * lsk_1327[k]
                    + f_3 * pc_x[k] * msk_1327[k];

        t_1652[k] = f_15 * lsk_1328[k]
                    + f_3 * pc_x[k] * msk_1328[k];
    }

#pragma omp simd aligned(t_1653, t_1654, t_1655, t_1656, pa_x, pc_x, lsl0_1656, lsk_1329, \
                         lsk_1330, lsk_1331, lsl1_1656, msk_1329, msk_1330, \
                         msk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1653[k] = f_15 * lsk_1329[k]
                    + f_3 * pc_x[k] * msk_1329[k];

        t_1654[k] = f_15 * lsk_1330[k]
                    + f_3 * pc_x[k] * msk_1330[k];

        t_1655[k] = f_15 * lsk_1331[k]
                    + f_3 * pc_x[k] * msk_1331[k];

        t_1656[k] = pa_x[k] * lsl0_1656[k]
                    - f_14 * pc_x[k] * lsl1_1656[k];
    }

#pragma omp simd aligned(t_1657, t_1658, t_1659, t_1660, pa_x, pc_x, pc_z, lsl0_1658, \
                         lsl0_1659, lsl0_1660, lsl1_1658, lsl1_1659, lsl1_1660, \
                         msk_1324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1657[k] = f_3 * pc_z[k] * msk_1324[k];

        t_1658[k] = pa_x[k] * lsl0_1658[k]
                    - f_14 * pc_x[k] * lsl1_1658[k];

        t_1659[k] = pa_x[k] * lsl0_1659[k]
                    - f_14 * pc_x[k] * lsl1_1659[k];

        t_1660[k] = pa_x[k] * lsl0_1660[k]
                    - f_14 * pc_x[k] * lsl1_1660[k];
    }

#pragma omp simd aligned(t_1661, t_1662, t_1663, t_1664, pa_x, pc_x, pc_y, lsl0_1661, \
                         lsl0_1662, lsl0_1664, lsk_1043, lsl1_1661, lsl1_1662, lsl1_1664, \
                         msk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1661[k] = pa_x[k] * lsl0_1661[k]
                    - f_14 * pc_x[k] * lsl1_1661[k];

        t_1662[k] = pa_x[k] * lsl0_1662[k]
                    - f_14 * pc_x[k] * lsl1_1662[k];

        t_1663[k] = f_21 * lsk_1043[k]
                    + f_3 * pc_y[k] * msk_1331[k];

        t_1664[k] = pa_x[k] * lsl0_1664[k]
                    - f_14 * pc_x[k] * lsl1_1664[k];
    }

#pragma omp simd aligned(t_1665, t_1666, t_1667, t_1668, pa_z, pc_y, pc_z, lsl0_1260, \
                         lsl0_1263, lsk_1008, lsk_1044, lsl1_1260, lsl1_1263, \
                         msk_1332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1665[k] = pa_z[k] * lsl0_1260[k]
                    - f_14 * pc_z[k] * lsl1_1260[k];

        t_1666[k] = f_24 * lsk_1044[k]
                    + f_3 * pc_y[k] * msk_1332[k];

        t_1667[k] = f_15 * lsk_1008[k]
                    + f_3 * pc_z[k] * msk_1332[k];

        t_1668[k] = pa_z[k] * lsl0_1263[k]
                    - f_14 * pc_z[k] * lsl1_1263[k];
    }

#pragma omp simd aligned(t_1669, t_1670, t_1671, pa_x, pa_z, pc_x, pc_y, pc_z, lsl0_1266, \
                         lsl0_1670, lsk_1046, lsk_1337, lsl1_1266, lsl1_1670, \
                         msk_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1669[k] = f_24 * lsk_1046[k]
                    + f_3 * pc_y[k] * msk_1334[k];

        t_1670[k] = pa_x[k] * lsl0_1670[k]
                    + f_20 * lsk_1337[k]
                    - f_14 * pc_x[k] * lsl1_1670[k];

        t_1671[k] = pa_z[k] * lsl0_1266[k]
                    - f_14 * pc_z[k] * lsl1_1266[k];
    }

#pragma omp simd aligned(t_1672, t_1673, t_1674, pa_x, pc_x, pc_y, pc_z, lsl0_1674, lsk_1011, \
                         lsk_1049, lsk_1341, lsl1_1674, msk_1335, \
                         msk_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1672[k] = f_15 * lsk_1011[k]
                    + f_3 * pc_z[k] * msk_1335[k];

        t_1673[k] = f_24 * lsk_1049[k]
                    + f_3 * pc_y[k] * msk_1337[k];

        t_1674[k] = pa_x[k] * lsl0_1674[k]
                    + f_19 * lsk_1341[k]
                    - f_14 * pc_x[k] * lsl1_1674[k];
    }

#pragma omp simd aligned(t_1675, t_1676, t_1677, pa_x, pa_z, pc_x, pc_z, lsl0_1270, lsl0_1677, \
                         lsk_1014, lsk_1344, lsl1_1270, lsl1_1677, \
                         msk_1338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1675[k] = pa_z[k] * lsl0_1270[k]
                    - f_14 * pc_z[k] * lsl1_1270[k];

        t_1676[k] = f_15 * lsk_1014[k]
                    + f_3 * pc_z[k] * msk_1338[k];

        t_1677[k] = pa_x[k] * lsl0_1677[k]
                    + f_18 * lsk_1344[k]
                    - f_14 * pc_x[k] * lsl1_1677[k];
    }

#pragma omp simd aligned(t_1678, t_1679, t_1680, pa_x, pa_z, pc_x, pc_y, pc_z, lsl0_1275, \
                         lsl0_1679, lsk_1053, lsk_1346, lsl1_1275, lsl1_1679, \
                         msk_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1678[k] = f_24 * lsk_1053[k]
                    + f_3 * pc_y[k] * msk_1341[k];

        t_1679[k] = pa_x[k] * lsl0_1679[k]
                    + f_18 * lsk_1346[k]
                    - f_14 * pc_x[k] * lsl1_1679[k];

        t_1680[k] = pa_z[k] * lsl0_1275[k]
                    - f_14 * pc_z[k] * lsl1_1275[k];
    }

#pragma omp simd aligned(t_1681, t_1682, t_1683, pa_x, pc_x, pc_z, lsl0_1682, lsl0_1683, \
                         lsk_1018, lsk_1349, lsk_1350, lsl1_1682, lsl1_1683, \
                         msk_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1681[k] = f_15 * lsk_1018[k]
                    + f_3 * pc_z[k] * msk_1342[k];

        t_1682[k] = pa_x[k] * lsl0_1682[k]
                    + f_17 * lsk_1349[k]
                    - f_14 * pc_x[k] * lsl1_1682[k];

        t_1683[k] = pa_x[k] * lsl0_1683[k]
                    + f_17 * lsk_1350[k]
                    - f_14 * pc_x[k] * lsl1_1683[k];
    }

#pragma omp simd aligned(t_1684, t_1685, t_1686, pa_x, pa_z, pc_x, pc_y, pc_z, lsl0_1281, \
                         lsl0_1685, lsk_1058, lsk_1352, lsl1_1281, lsl1_1685, \
                         msk_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1684[k] = f_24 * lsk_1058[k]
                    + f_3 * pc_y[k] * msk_1346[k];

        t_1685[k] = pa_x[k] * lsl0_1685[k]
                    + f_17 * lsk_1352[k]
                    - f_14 * pc_x[k] * lsl1_1685[k];

        t_1686[k] = pa_z[k] * lsl0_1281[k]
                    - f_14 * pc_z[k] * lsl1_1281[k];
    }

#pragma omp simd aligned(t_1687, t_1688, t_1689, pa_x, pc_x, pc_z, lsl0_1688, lsl0_1689, \
                         lsk_1023, lsk_1355, lsk_1356, lsl1_1688, lsl1_1689, \
                         msk_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1687[k] = f_15 * lsk_1023[k]
                    + f_3 * pc_z[k] * msk_1347[k];

        t_1688[k] = pa_x[k] * lsl0_1688[k]
                    + f_16 * lsk_1355[k]
                    - f_14 * pc_x[k] * lsl1_1688[k];

        t_1689[k] = pa_x[k] * lsl0_1689[k]
                    + f_16 * lsk_1356[k]
                    - f_14 * pc_x[k] * lsl1_1689[k];
    }

#pragma omp simd aligned(t_1690, t_1691, t_1692, pa_x, pc_x, pc_y, lsl0_1690, lsl0_1692, \
                         lsk_1064, lsk_1357, lsk_1359, lsl1_1690, lsl1_1692, \
                         msk_1352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1690[k] = pa_x[k] * lsl0_1690[k]
                    + f_16 * lsk_1357[k]
                    - f_14 * pc_x[k] * lsl1_1690[k];

        t_1691[k] = f_24 * lsk_1064[k]
                    + f_3 * pc_y[k] * msk_1352[k];

        t_1692[k] = pa_x[k] * lsl0_1692[k]
                    + f_16 * lsk_1359[k]
                    - f_14 * pc_x[k] * lsl1_1692[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msk, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_24 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1701 = buffer.data(lsl0 + 1701);
    const auto *lsl0_1703 = buffer.data(lsl0 + 1703);
    const auto *lsl0_1704 = buffer.data(lsl0 + 1704);
    const auto *lsl0_1705 = buffer.data(lsl0 + 1705);
    const auto *lsl0_1706 = buffer.data(lsl0 + 1706);
    const auto *lsl0_1707 = buffer.data(lsl0 + 1707);
    const auto *lsl0_1709 = buffer.data(lsl0 + 1709);
    const auto *lsl0_1710 = buffer.data(lsl0 + 1710);
    const auto *lsl0_1713 = buffer.data(lsl0 + 1713);
    const auto *lsl0_1715 = buffer.data(lsl0 + 1715);
    const auto *lsl0_1716 = buffer.data(lsl0 + 1716);
    const auto *lsl0_1719 = buffer.data(lsl0 + 1719);
    const auto *lsl0_1720 = buffer.data(lsl0 + 1720);
    const auto *lsl0_1722 = buffer.data(lsl0 + 1722);
    const auto *lsl0_1724 = buffer.data(lsl0 + 1724);
    const auto *lsl0_1725 = buffer.data(lsl0 + 1725);
    const auto *lsl0_1727 = buffer.data(lsl0 + 1727);
    const auto *lsl0_1728 = buffer.data(lsl0 + 1728);
    const auto *lsl0_1730 = buffer.data(lsl0 + 1730);
    const auto *lsl0_1731 = buffer.data(lsl0 + 1731);
    const auto *lsl0_1733 = buffer.data(lsl0 + 1733);
    const auto *lsl0_1734 = buffer.data(lsl0 + 1734);
    const auto *lsl0_1735 = buffer.data(lsl0 + 1735);
    const auto *lsl0_1737 = buffer.data(lsl0 + 1737);
    const auto *lsl0_1746 = buffer.data(lsl0 + 1746);
    const auto *lsl0_1748 = buffer.data(lsl0 + 1748);
    const auto *lsl0_1749 = buffer.data(lsl0 + 1749);
    const auto *lsl0_1750 = buffer.data(lsl0 + 1750);
    const auto *lsl0_1751 = buffer.data(lsl0 + 1751);
    const auto *lsl0_1752 = buffer.data(lsl0 + 1752);
    const auto *lsl0_1754 = buffer.data(lsl0 + 1754);
    const auto *lsl0_1755 = buffer.data(lsl0 + 1755);
    const auto *lsl0_1758 = buffer.data(lsl0 + 1758);
    const auto *lsl0_1760 = buffer.data(lsl0 + 1760);
    const auto *lsl0_1761 = buffer.data(lsl0 + 1761);
    const auto *lsl0_1764 = buffer.data(lsl0 + 1764);
    const auto *lsl0_1765 = buffer.data(lsl0 + 1765);
    const auto *lsl0_1767 = buffer.data(lsl0 + 1767);
    const auto *lsl0_1769 = buffer.data(lsl0 + 1769);
    const auto *lsl0_1770 = buffer.data(lsl0 + 1770);
    const auto *lsl0_1772 = buffer.data(lsl0 + 1772);
    const auto *lsl0_1773 = buffer.data(lsl0 + 1773);
    const auto *lsl0_1775 = buffer.data(lsl0 + 1775);
    const auto *lsl0_1776 = buffer.data(lsl0 + 1776);
    const auto *lsl0_1778 = buffer.data(lsl0 + 1778);
    const auto *lsl0_1779 = buffer.data(lsl0 + 1779);
    const auto *lsl0_1780 = buffer.data(lsl0 + 1780);
    const auto *lsl0_1782 = buffer.data(lsl0 + 1782);
    const auto *lsl0_1791 = buffer.data(lsl0 + 1791);
    const auto *lsl0_1793 = buffer.data(lsl0 + 1793);
    const auto *lsl0_1794 = buffer.data(lsl0 + 1794);
    const auto *lsl0_1795 = buffer.data(lsl0 + 1795);
    const auto *lsl0_1796 = buffer.data(lsl0 + 1796);
    const auto *lsl0_1797 = buffer.data(lsl0 + 1797);
    const auto *lsl0_1799 = buffer.data(lsl0 + 1799);
    const auto *lsl0_1800 = buffer.data(lsl0 + 1800);
    const auto *lsl0_1803 = buffer.data(lsl0 + 1803);
    const auto *lsl0_1805 = buffer.data(lsl0 + 1805);
    const auto *lsl0_1806 = buffer.data(lsl0 + 1806);
    const auto *lsl0_1809 = buffer.data(lsl0 + 1809);
    const auto *lsl0_1810 = buffer.data(lsl0 + 1810);

    const auto *lsk_1036 = buffer.data(lsk + 1036);
    const auto *lsk_1044 = buffer.data(lsk + 1044);
    const auto *lsk_1047 = buffer.data(lsk + 1047);
    const auto *lsk_1050 = buffer.data(lsk + 1050);
    const auto *lsk_1054 = buffer.data(lsk + 1054);
    const auto *lsk_1059 = buffer.data(lsk + 1059);
    const auto *lsk_1072 = buffer.data(lsk + 1072);
    const auto *lsk_1079 = buffer.data(lsk + 1079);
    const auto *lsk_1080 = buffer.data(lsk + 1080);
    const auto *lsk_1082 = buffer.data(lsk + 1082);
    const auto *lsk_1083 = buffer.data(lsk + 1083);
    const auto *lsk_1085 = buffer.data(lsk + 1085);
    const auto *lsk_1086 = buffer.data(lsk + 1086);
    const auto *lsk_1089 = buffer.data(lsk + 1089);
    const auto *lsk_1090 = buffer.data(lsk + 1090);
    const auto *lsk_1094 = buffer.data(lsk + 1094);
    const auto *lsk_1095 = buffer.data(lsk + 1095);
    const auto *lsk_1100 = buffer.data(lsk + 1100);
    const auto *lsk_1108 = buffer.data(lsk + 1108);
    const auto *lsk_1115 = buffer.data(lsk + 1115);
    const auto *lsk_1116 = buffer.data(lsk + 1116);
    const auto *lsk_1118 = buffer.data(lsk + 1118);
    const auto *lsk_1119 = buffer.data(lsk + 1119);
    const auto *lsk_1121 = buffer.data(lsk + 1121);
    const auto *lsk_1122 = buffer.data(lsk + 1122);
    const auto *lsk_1125 = buffer.data(lsk + 1125);
    const auto *lsk_1130 = buffer.data(lsk + 1130);
    const auto *lsk_1136 = buffer.data(lsk + 1136);
    const auto *lsk_1151 = buffer.data(lsk + 1151);
    const auto *lsk_1152 = buffer.data(lsk + 1152);
    const auto *lsk_1154 = buffer.data(lsk + 1154);
    const auto *lsk_1157 = buffer.data(lsk + 1157);
    const auto *lsk_1360 = buffer.data(lsk + 1360);
    const auto *lsk_1361 = buffer.data(lsk + 1361);
    const auto *lsk_1362 = buffer.data(lsk + 1362);
    const auto *lsk_1363 = buffer.data(lsk + 1363);
    const auto *lsk_1364 = buffer.data(lsk + 1364);
    const auto *lsk_1365 = buffer.data(lsk + 1365);
    const auto *lsk_1366 = buffer.data(lsk + 1366);
    const auto *lsk_1367 = buffer.data(lsk + 1367);
    const auto *lsk_1368 = buffer.data(lsk + 1368);
    const auto *lsk_1371 = buffer.data(lsk + 1371);
    const auto *lsk_1373 = buffer.data(lsk + 1373);
    const auto *lsk_1374 = buffer.data(lsk + 1374);
    const auto *lsk_1377 = buffer.data(lsk + 1377);
    const auto *lsk_1378 = buffer.data(lsk + 1378);
    const auto *lsk_1380 = buffer.data(lsk + 1380);
    const auto *lsk_1382 = buffer.data(lsk + 1382);
    const auto *lsk_1383 = buffer.data(lsk + 1383);
    const auto *lsk_1385 = buffer.data(lsk + 1385);
    const auto *lsk_1386 = buffer.data(lsk + 1386);
    const auto *lsk_1388 = buffer.data(lsk + 1388);
    const auto *lsk_1389 = buffer.data(lsk + 1389);
    const auto *lsk_1391 = buffer.data(lsk + 1391);
    const auto *lsk_1392 = buffer.data(lsk + 1392);
    const auto *lsk_1393 = buffer.data(lsk + 1393);
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
    const auto *lsk_1407 = buffer.data(lsk + 1407);
    const auto *lsk_1409 = buffer.data(lsk + 1409);
    const auto *lsk_1410 = buffer.data(lsk + 1410);
    const auto *lsk_1413 = buffer.data(lsk + 1413);
    const auto *lsk_1414 = buffer.data(lsk + 1414);
    const auto *lsk_1416 = buffer.data(lsk + 1416);
    const auto *lsk_1418 = buffer.data(lsk + 1418);
    const auto *lsk_1419 = buffer.data(lsk + 1419);
    const auto *lsk_1421 = buffer.data(lsk + 1421);
    const auto *lsk_1422 = buffer.data(lsk + 1422);
    const auto *lsk_1424 = buffer.data(lsk + 1424);
    const auto *lsk_1425 = buffer.data(lsk + 1425);
    const auto *lsk_1427 = buffer.data(lsk + 1427);
    const auto *lsk_1428 = buffer.data(lsk + 1428);
    const auto *lsk_1429 = buffer.data(lsk + 1429);
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
    const auto *lsk_1443 = buffer.data(lsk + 1443);
    const auto *lsk_1445 = buffer.data(lsk + 1445);
    const auto *lsk_1446 = buffer.data(lsk + 1446);
    const auto *lsk_1449 = buffer.data(lsk + 1449);
    const auto *lsk_1450 = buffer.data(lsk + 1450);

    const auto *lsl1_1701 = buffer.data(lsl1 + 1701);
    const auto *lsl1_1703 = buffer.data(lsl1 + 1703);
    const auto *lsl1_1704 = buffer.data(lsl1 + 1704);
    const auto *lsl1_1705 = buffer.data(lsl1 + 1705);
    const auto *lsl1_1706 = buffer.data(lsl1 + 1706);
    const auto *lsl1_1707 = buffer.data(lsl1 + 1707);
    const auto *lsl1_1709 = buffer.data(lsl1 + 1709);
    const auto *lsl1_1710 = buffer.data(lsl1 + 1710);
    const auto *lsl1_1713 = buffer.data(lsl1 + 1713);
    const auto *lsl1_1715 = buffer.data(lsl1 + 1715);
    const auto *lsl1_1716 = buffer.data(lsl1 + 1716);
    const auto *lsl1_1719 = buffer.data(lsl1 + 1719);
    const auto *lsl1_1720 = buffer.data(lsl1 + 1720);
    const auto *lsl1_1722 = buffer.data(lsl1 + 1722);
    const auto *lsl1_1724 = buffer.data(lsl1 + 1724);
    const auto *lsl1_1725 = buffer.data(lsl1 + 1725);
    const auto *lsl1_1727 = buffer.data(lsl1 + 1727);
    const auto *lsl1_1728 = buffer.data(lsl1 + 1728);
    const auto *lsl1_1730 = buffer.data(lsl1 + 1730);
    const auto *lsl1_1731 = buffer.data(lsl1 + 1731);
    const auto *lsl1_1733 = buffer.data(lsl1 + 1733);
    const auto *lsl1_1734 = buffer.data(lsl1 + 1734);
    const auto *lsl1_1735 = buffer.data(lsl1 + 1735);
    const auto *lsl1_1737 = buffer.data(lsl1 + 1737);
    const auto *lsl1_1746 = buffer.data(lsl1 + 1746);
    const auto *lsl1_1748 = buffer.data(lsl1 + 1748);
    const auto *lsl1_1749 = buffer.data(lsl1 + 1749);
    const auto *lsl1_1750 = buffer.data(lsl1 + 1750);
    const auto *lsl1_1751 = buffer.data(lsl1 + 1751);
    const auto *lsl1_1752 = buffer.data(lsl1 + 1752);
    const auto *lsl1_1754 = buffer.data(lsl1 + 1754);
    const auto *lsl1_1755 = buffer.data(lsl1 + 1755);
    const auto *lsl1_1758 = buffer.data(lsl1 + 1758);
    const auto *lsl1_1760 = buffer.data(lsl1 + 1760);
    const auto *lsl1_1761 = buffer.data(lsl1 + 1761);
    const auto *lsl1_1764 = buffer.data(lsl1 + 1764);
    const auto *lsl1_1765 = buffer.data(lsl1 + 1765);
    const auto *lsl1_1767 = buffer.data(lsl1 + 1767);
    const auto *lsl1_1769 = buffer.data(lsl1 + 1769);
    const auto *lsl1_1770 = buffer.data(lsl1 + 1770);
    const auto *lsl1_1772 = buffer.data(lsl1 + 1772);
    const auto *lsl1_1773 = buffer.data(lsl1 + 1773);
    const auto *lsl1_1775 = buffer.data(lsl1 + 1775);
    const auto *lsl1_1776 = buffer.data(lsl1 + 1776);
    const auto *lsl1_1778 = buffer.data(lsl1 + 1778);
    const auto *lsl1_1779 = buffer.data(lsl1 + 1779);
    const auto *lsl1_1780 = buffer.data(lsl1 + 1780);
    const auto *lsl1_1782 = buffer.data(lsl1 + 1782);
    const auto *lsl1_1791 = buffer.data(lsl1 + 1791);
    const auto *lsl1_1793 = buffer.data(lsl1 + 1793);
    const auto *lsl1_1794 = buffer.data(lsl1 + 1794);
    const auto *lsl1_1795 = buffer.data(lsl1 + 1795);
    const auto *lsl1_1796 = buffer.data(lsl1 + 1796);
    const auto *lsl1_1797 = buffer.data(lsl1 + 1797);
    const auto *lsl1_1799 = buffer.data(lsl1 + 1799);
    const auto *lsl1_1800 = buffer.data(lsl1 + 1800);
    const auto *lsl1_1803 = buffer.data(lsl1 + 1803);
    const auto *lsl1_1805 = buffer.data(lsl1 + 1805);
    const auto *lsl1_1806 = buffer.data(lsl1 + 1806);
    const auto *lsl1_1809 = buffer.data(lsl1 + 1809);
    const auto *lsl1_1810 = buffer.data(lsl1 + 1810);

    const auto *msk_1360 = buffer.data(msk + 1360);
    const auto *msk_1361 = buffer.data(msk + 1361);
    const auto *msk_1362 = buffer.data(msk + 1362);
    const auto *msk_1363 = buffer.data(msk + 1363);
    const auto *msk_1364 = buffer.data(msk + 1364);
    const auto *msk_1365 = buffer.data(msk + 1365);
    const auto *msk_1366 = buffer.data(msk + 1366);
    const auto *msk_1367 = buffer.data(msk + 1367);
    const auto *msk_1368 = buffer.data(msk + 1368);
    const auto *msk_1370 = buffer.data(msk + 1370);
    const auto *msk_1371 = buffer.data(msk + 1371);
    const auto *msk_1373 = buffer.data(msk + 1373);
    const auto *msk_1374 = buffer.data(msk + 1374);
    const auto *msk_1377 = buffer.data(msk + 1377);
    const auto *msk_1378 = buffer.data(msk + 1378);
    const auto *msk_1382 = buffer.data(msk + 1382);
    const auto *msk_1383 = buffer.data(msk + 1383);
    const auto *msk_1388 = buffer.data(msk + 1388);
    const auto *msk_1396 = buffer.data(msk + 1396);
    const auto *msk_1397 = buffer.data(msk + 1397);
    const auto *msk_1398 = buffer.data(msk + 1398);
    const auto *msk_1399 = buffer.data(msk + 1399);
    const auto *msk_1400 = buffer.data(msk + 1400);
    const auto *msk_1401 = buffer.data(msk + 1401);
    const auto *msk_1402 = buffer.data(msk + 1402);
    const auto *msk_1403 = buffer.data(msk + 1403);
    const auto *msk_1404 = buffer.data(msk + 1404);
    const auto *msk_1406 = buffer.data(msk + 1406);
    const auto *msk_1407 = buffer.data(msk + 1407);
    const auto *msk_1409 = buffer.data(msk + 1409);
    const auto *msk_1410 = buffer.data(msk + 1410);
    const auto *msk_1413 = buffer.data(msk + 1413);
    const auto *msk_1414 = buffer.data(msk + 1414);
    const auto *msk_1418 = buffer.data(msk + 1418);
    const auto *msk_1419 = buffer.data(msk + 1419);
    const auto *msk_1424 = buffer.data(msk + 1424);
    const auto *msk_1432 = buffer.data(msk + 1432);
    const auto *msk_1433 = buffer.data(msk + 1433);
    const auto *msk_1434 = buffer.data(msk + 1434);
    const auto *msk_1435 = buffer.data(msk + 1435);
    const auto *msk_1436 = buffer.data(msk + 1436);
    const auto *msk_1437 = buffer.data(msk + 1437);
    const auto *msk_1438 = buffer.data(msk + 1438);
    const auto *msk_1439 = buffer.data(msk + 1439);
    const auto *msk_1440 = buffer.data(msk + 1440);
    const auto *msk_1442 = buffer.data(msk + 1442);
    const auto *msk_1443 = buffer.data(msk + 1443);
    const auto *msk_1445 = buffer.data(msk + 1445);
    const auto *msk_1446 = buffer.data(msk + 1446);

#pragma omp simd aligned(t_1693, t_1694, t_1695, t_1696, t_1697, pc_x, lsk_1360, lsk_1361, \
                         lsk_1362, lsk_1363, lsk_1364, msk_1360, msk_1361, msk_1362, msk_1363, \
                         msk_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1693[k] = f_15 * lsk_1360[k]
                    + f_3 * pc_x[k] * msk_1360[k];

        t_1694[k] = f_15 * lsk_1361[k]
                    + f_3 * pc_x[k] * msk_1361[k];

        t_1695[k] = f_15 * lsk_1362[k]
                    + f_3 * pc_x[k] * msk_1362[k];

        t_1696[k] = f_15 * lsk_1363[k]
                    + f_3 * pc_x[k] * msk_1363[k];

        t_1697[k] = f_15 * lsk_1364[k]
                    + f_3 * pc_x[k] * msk_1364[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, t_1701, pa_x, pc_x, lsl0_1701, lsk_1365, \
                         lsk_1366, lsk_1367, lsl1_1701, msk_1365, msk_1366, \
                         msk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_15 * lsk_1365[k]
                    + f_3 * pc_x[k] * msk_1365[k];

        t_1699[k] = f_15 * lsk_1366[k]
                    + f_3 * pc_x[k] * msk_1366[k];

        t_1700[k] = f_15 * lsk_1367[k]
                    + f_3 * pc_x[k] * msk_1367[k];

        t_1701[k] = pa_x[k] * lsl0_1701[k]
                    - f_14 * pc_x[k] * lsl1_1701[k];
    }

#pragma omp simd aligned(t_1702, t_1703, t_1704, t_1705, pa_x, pc_x, pc_z, lsl0_1703, \
                         lsl0_1704, lsl0_1705, lsk_1036, lsl1_1703, lsl1_1704, lsl1_1705, \
                         msk_1360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1702[k] = f_15 * lsk_1036[k]
                    + f_3 * pc_z[k] * msk_1360[k];

        t_1703[k] = pa_x[k] * lsl0_1703[k]
                    - f_14 * pc_x[k] * lsl1_1703[k];

        t_1704[k] = pa_x[k] * lsl0_1704[k]
                    - f_14 * pc_x[k] * lsl1_1704[k];

        t_1705[k] = pa_x[k] * lsl0_1705[k]
                    - f_14 * pc_x[k] * lsl1_1705[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, t_1709, pa_x, pc_x, pc_y, lsl0_1706, \
                         lsl0_1707, lsl0_1709, lsk_1079, lsl1_1706, lsl1_1707, lsl1_1709, \
                         msk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = pa_x[k] * lsl0_1706[k]
                    - f_14 * pc_x[k] * lsl1_1706[k];

        t_1707[k] = pa_x[k] * lsl0_1707[k]
                    - f_14 * pc_x[k] * lsl1_1707[k];

        t_1708[k] = f_24 * lsk_1079[k]
                    + f_3 * pc_y[k] * msk_1367[k];

        t_1709[k] = pa_x[k] * lsl0_1709[k]
                    - f_14 * pc_x[k] * lsl1_1709[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, pa_x, pc_x, pc_y, pc_z, lsl0_1710, lsk_1044, \
                         lsk_1080, lsk_1368, lsl1_1710, msk_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = pa_x[k] * lsl0_1710[k]
                    + f_21 * lsk_1368[k]
                    - f_14 * pc_x[k] * lsl1_1710[k];

        t_1711[k] = f_20 * lsk_1080[k]
                    + f_3 * pc_y[k] * msk_1368[k];

        t_1712[k] = f_16 * lsk_1044[k]
                    + f_3 * pc_z[k] * msk_1368[k];
    }

#pragma omp simd aligned(t_1713, t_1714, t_1715, pa_x, pc_x, pc_y, lsl0_1713, lsl0_1715, \
                         lsk_1082, lsk_1371, lsk_1373, lsl1_1713, lsl1_1715, \
                         msk_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1713[k] = pa_x[k] * lsl0_1713[k]
                    + f_20 * lsk_1371[k]
                    - f_14 * pc_x[k] * lsl1_1713[k];

        t_1714[k] = f_20 * lsk_1082[k]
                    + f_3 * pc_y[k] * msk_1370[k];

        t_1715[k] = pa_x[k] * lsl0_1715[k]
                    + f_20 * lsk_1373[k]
                    - f_14 * pc_x[k] * lsl1_1715[k];
    }

#pragma omp simd aligned(t_1716, t_1717, t_1718, pa_x, pc_x, pc_y, pc_z, lsl0_1716, lsk_1047, \
                         lsk_1085, lsk_1374, lsl1_1716, msk_1371, \
                         msk_1373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1716[k] = pa_x[k] * lsl0_1716[k]
                    + f_19 * lsk_1374[k]
                    - f_14 * pc_x[k] * lsl1_1716[k];

        t_1717[k] = f_16 * lsk_1047[k]
                    + f_3 * pc_z[k] * msk_1371[k];

        t_1718[k] = f_20 * lsk_1085[k]
                    + f_3 * pc_y[k] * msk_1373[k];
    }

#pragma omp simd aligned(t_1719, t_1720, t_1721, pa_x, pc_x, pc_z, lsl0_1719, lsl0_1720, \
                         lsk_1050, lsk_1377, lsk_1378, lsl1_1719, lsl1_1720, \
                         msk_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1719[k] = pa_x[k] * lsl0_1719[k]
                    + f_19 * lsk_1377[k]
                    - f_14 * pc_x[k] * lsl1_1719[k];

        t_1720[k] = pa_x[k] * lsl0_1720[k]
                    + f_18 * lsk_1378[k]
                    - f_14 * pc_x[k] * lsl1_1720[k];

        t_1721[k] = f_16 * lsk_1050[k]
                    + f_3 * pc_z[k] * msk_1374[k];
    }

#pragma omp simd aligned(t_1722, t_1723, t_1724, pa_x, pc_x, pc_y, lsl0_1722, lsl0_1724, \
                         lsk_1089, lsk_1380, lsk_1382, lsl1_1722, lsl1_1724, \
                         msk_1377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1722[k] = pa_x[k] * lsl0_1722[k]
                    + f_18 * lsk_1380[k]
                    - f_14 * pc_x[k] * lsl1_1722[k];

        t_1723[k] = f_20 * lsk_1089[k]
                    + f_3 * pc_y[k] * msk_1377[k];

        t_1724[k] = pa_x[k] * lsl0_1724[k]
                    + f_18 * lsk_1382[k]
                    - f_14 * pc_x[k] * lsl1_1724[k];
    }

#pragma omp simd aligned(t_1725, t_1726, t_1727, pa_x, pc_x, pc_z, lsl0_1725, lsl0_1727, \
                         lsk_1054, lsk_1383, lsk_1385, lsl1_1725, lsl1_1727, \
                         msk_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1725[k] = pa_x[k] * lsl0_1725[k]
                    + f_17 * lsk_1383[k]
                    - f_14 * pc_x[k] * lsl1_1725[k];

        t_1726[k] = f_16 * lsk_1054[k]
                    + f_3 * pc_z[k] * msk_1378[k];

        t_1727[k] = pa_x[k] * lsl0_1727[k]
                    + f_17 * lsk_1385[k]
                    - f_14 * pc_x[k] * lsl1_1727[k];
    }

#pragma omp simd aligned(t_1728, t_1729, t_1730, pa_x, pc_x, pc_y, lsl0_1728, lsl0_1730, \
                         lsk_1094, lsk_1386, lsk_1388, lsl1_1728, lsl1_1730, \
                         msk_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1728[k] = pa_x[k] * lsl0_1728[k]
                    + f_17 * lsk_1386[k]
                    - f_14 * pc_x[k] * lsl1_1728[k];

        t_1729[k] = f_20 * lsk_1094[k]
                    + f_3 * pc_y[k] * msk_1382[k];

        t_1730[k] = pa_x[k] * lsl0_1730[k]
                    + f_17 * lsk_1388[k]
                    - f_14 * pc_x[k] * lsl1_1730[k];
    }

#pragma omp simd aligned(t_1731, t_1732, t_1733, pa_x, pc_x, pc_z, lsl0_1731, lsl0_1733, \
                         lsk_1059, lsk_1389, lsk_1391, lsl1_1731, lsl1_1733, \
                         msk_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1731[k] = pa_x[k] * lsl0_1731[k]
                    + f_16 * lsk_1389[k]
                    - f_14 * pc_x[k] * lsl1_1731[k];

        t_1732[k] = f_16 * lsk_1059[k]
                    + f_3 * pc_z[k] * msk_1383[k];

        t_1733[k] = pa_x[k] * lsl0_1733[k]
                    + f_16 * lsk_1391[k]
                    - f_14 * pc_x[k] * lsl1_1733[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, pa_x, pc_x, pc_y, lsl0_1734, lsl0_1735, \
                         lsk_1100, lsk_1392, lsk_1393, lsl1_1734, lsl1_1735, \
                         msk_1388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = pa_x[k] * lsl0_1734[k]
                    + f_16 * lsk_1392[k]
                    - f_14 * pc_x[k] * lsl1_1734[k];

        t_1735[k] = pa_x[k] * lsl0_1735[k]
                    + f_16 * lsk_1393[k]
                    - f_14 * pc_x[k] * lsl1_1735[k];

        t_1736[k] = f_20 * lsk_1100[k]
                    + f_3 * pc_y[k] * msk_1388[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, t_1740, pa_x, pc_x, lsl0_1737, lsk_1395, \
                         lsk_1396, lsk_1397, lsk_1398, lsl1_1737, msk_1396, msk_1397, \
                         msk_1398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = pa_x[k] * lsl0_1737[k]
                    + f_16 * lsk_1395[k]
                    - f_14 * pc_x[k] * lsl1_1737[k];

        t_1738[k] = f_15 * lsk_1396[k]
                    + f_3 * pc_x[k] * msk_1396[k];

        t_1739[k] = f_15 * lsk_1397[k]
                    + f_3 * pc_x[k] * msk_1397[k];

        t_1740[k] = f_15 * lsk_1398[k]
                    + f_3 * pc_x[k] * msk_1398[k];
    }

#pragma omp simd aligned(t_1741, t_1742, t_1743, t_1744, t_1745, pc_x, lsk_1399, lsk_1400, \
                         lsk_1401, lsk_1402, lsk_1403, msk_1399, msk_1400, msk_1401, msk_1402, \
                         msk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1741[k] = f_15 * lsk_1399[k]
                    + f_3 * pc_x[k] * msk_1399[k];

        t_1742[k] = f_15 * lsk_1400[k]
                    + f_3 * pc_x[k] * msk_1400[k];

        t_1743[k] = f_15 * lsk_1401[k]
                    + f_3 * pc_x[k] * msk_1401[k];

        t_1744[k] = f_15 * lsk_1402[k]
                    + f_3 * pc_x[k] * msk_1402[k];

        t_1745[k] = f_15 * lsk_1403[k]
                    + f_3 * pc_x[k] * msk_1403[k];
    }

#pragma omp simd aligned(t_1746, t_1747, t_1748, t_1749, pa_x, pc_x, pc_z, lsl0_1746, \
                         lsl0_1748, lsl0_1749, lsk_1072, lsl1_1746, lsl1_1748, lsl1_1749, \
                         msk_1396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1746[k] = pa_x[k] * lsl0_1746[k]
                    - f_14 * pc_x[k] * lsl1_1746[k];

        t_1747[k] = f_16 * lsk_1072[k]
                    + f_3 * pc_z[k] * msk_1396[k];

        t_1748[k] = pa_x[k] * lsl0_1748[k]
                    - f_14 * pc_x[k] * lsl1_1748[k];

        t_1749[k] = pa_x[k] * lsl0_1749[k]
                    - f_14 * pc_x[k] * lsl1_1749[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, pa_x, pc_x, pc_y, lsl0_1750, \
                         lsl0_1751, lsl0_1752, lsk_1115, lsl1_1750, lsl1_1751, lsl1_1752, \
                         msk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = pa_x[k] * lsl0_1750[k]
                    - f_14 * pc_x[k] * lsl1_1750[k];

        t_1751[k] = pa_x[k] * lsl0_1751[k]
                    - f_14 * pc_x[k] * lsl1_1751[k];

        t_1752[k] = pa_x[k] * lsl0_1752[k]
                    - f_14 * pc_x[k] * lsl1_1752[k];

        t_1753[k] = f_20 * lsk_1115[k]
                    + f_3 * pc_y[k] * msk_1403[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, pa_x, pc_x, pc_y, pc_z, lsl0_1754, \
                         lsl0_1755, lsk_1080, lsk_1116, lsk_1404, lsl1_1754, lsl1_1755, \
                         msk_1404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = pa_x[k] * lsl0_1754[k]
                    - f_14 * pc_x[k] * lsl1_1754[k];

        t_1755[k] = pa_x[k] * lsl0_1755[k]
                    + f_21 * lsk_1404[k]
                    - f_14 * pc_x[k] * lsl1_1755[k];

        t_1756[k] = f_19 * lsk_1116[k]
                    + f_3 * pc_y[k] * msk_1404[k];

        t_1757[k] = f_17 * lsk_1080[k]
                    + f_3 * pc_z[k] * msk_1404[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, pa_x, pc_x, pc_y, lsl0_1758, lsl0_1760, \
                         lsk_1118, lsk_1407, lsk_1409, lsl1_1758, lsl1_1760, \
                         msk_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = pa_x[k] * lsl0_1758[k]
                    + f_20 * lsk_1407[k]
                    - f_14 * pc_x[k] * lsl1_1758[k];

        t_1759[k] = f_19 * lsk_1118[k]
                    + f_3 * pc_y[k] * msk_1406[k];

        t_1760[k] = pa_x[k] * lsl0_1760[k]
                    + f_20 * lsk_1409[k]
                    - f_14 * pc_x[k] * lsl1_1760[k];
    }

#pragma omp simd aligned(t_1761, t_1762, t_1763, pa_x, pc_x, pc_y, pc_z, lsl0_1761, lsk_1083, \
                         lsk_1121, lsk_1410, lsl1_1761, msk_1407, \
                         msk_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1761[k] = pa_x[k] * lsl0_1761[k]
                    + f_19 * lsk_1410[k]
                    - f_14 * pc_x[k] * lsl1_1761[k];

        t_1762[k] = f_17 * lsk_1083[k]
                    + f_3 * pc_z[k] * msk_1407[k];

        t_1763[k] = f_19 * lsk_1121[k]
                    + f_3 * pc_y[k] * msk_1409[k];
    }

#pragma omp simd aligned(t_1764, t_1765, t_1766, pa_x, pc_x, pc_z, lsl0_1764, lsl0_1765, \
                         lsk_1086, lsk_1413, lsk_1414, lsl1_1764, lsl1_1765, \
                         msk_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1764[k] = pa_x[k] * lsl0_1764[k]
                    + f_19 * lsk_1413[k]
                    - f_14 * pc_x[k] * lsl1_1764[k];

        t_1765[k] = pa_x[k] * lsl0_1765[k]
                    + f_18 * lsk_1414[k]
                    - f_14 * pc_x[k] * lsl1_1765[k];

        t_1766[k] = f_17 * lsk_1086[k]
                    + f_3 * pc_z[k] * msk_1410[k];
    }

#pragma omp simd aligned(t_1767, t_1768, t_1769, pa_x, pc_x, pc_y, lsl0_1767, lsl0_1769, \
                         lsk_1125, lsk_1416, lsk_1418, lsl1_1767, lsl1_1769, \
                         msk_1413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1767[k] = pa_x[k] * lsl0_1767[k]
                    + f_18 * lsk_1416[k]
                    - f_14 * pc_x[k] * lsl1_1767[k];

        t_1768[k] = f_19 * lsk_1125[k]
                    + f_3 * pc_y[k] * msk_1413[k];

        t_1769[k] = pa_x[k] * lsl0_1769[k]
                    + f_18 * lsk_1418[k]
                    - f_14 * pc_x[k] * lsl1_1769[k];
    }

#pragma omp simd aligned(t_1770, t_1771, t_1772, pa_x, pc_x, pc_z, lsl0_1770, lsl0_1772, \
                         lsk_1090, lsk_1419, lsk_1421, lsl1_1770, lsl1_1772, \
                         msk_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1770[k] = pa_x[k] * lsl0_1770[k]
                    + f_17 * lsk_1419[k]
                    - f_14 * pc_x[k] * lsl1_1770[k];

        t_1771[k] = f_17 * lsk_1090[k]
                    + f_3 * pc_z[k] * msk_1414[k];

        t_1772[k] = pa_x[k] * lsl0_1772[k]
                    + f_17 * lsk_1421[k]
                    - f_14 * pc_x[k] * lsl1_1772[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pa_x, pc_x, pc_y, lsl0_1773, lsl0_1775, \
                         lsk_1130, lsk_1422, lsk_1424, lsl1_1773, lsl1_1775, \
                         msk_1418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = pa_x[k] * lsl0_1773[k]
                    + f_17 * lsk_1422[k]
                    - f_14 * pc_x[k] * lsl1_1773[k];

        t_1774[k] = f_19 * lsk_1130[k]
                    + f_3 * pc_y[k] * msk_1418[k];

        t_1775[k] = pa_x[k] * lsl0_1775[k]
                    + f_17 * lsk_1424[k]
                    - f_14 * pc_x[k] * lsl1_1775[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pa_x, pc_x, pc_z, lsl0_1776, lsl0_1778, \
                         lsk_1095, lsk_1425, lsk_1427, lsl1_1776, lsl1_1778, \
                         msk_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = pa_x[k] * lsl0_1776[k]
                    + f_16 * lsk_1425[k]
                    - f_14 * pc_x[k] * lsl1_1776[k];

        t_1777[k] = f_17 * lsk_1095[k]
                    + f_3 * pc_z[k] * msk_1419[k];

        t_1778[k] = pa_x[k] * lsl0_1778[k]
                    + f_16 * lsk_1427[k]
                    - f_14 * pc_x[k] * lsl1_1778[k];
    }

#pragma omp simd aligned(t_1779, t_1780, t_1781, pa_x, pc_x, pc_y, lsl0_1779, lsl0_1780, \
                         lsk_1136, lsk_1428, lsk_1429, lsl1_1779, lsl1_1780, \
                         msk_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = pa_x[k] * lsl0_1779[k]
                    + f_16 * lsk_1428[k]
                    - f_14 * pc_x[k] * lsl1_1779[k];

        t_1780[k] = pa_x[k] * lsl0_1780[k]
                    + f_16 * lsk_1429[k]
                    - f_14 * pc_x[k] * lsl1_1780[k];

        t_1781[k] = f_19 * lsk_1136[k]
                    + f_3 * pc_y[k] * msk_1424[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, pa_x, pc_x, lsl0_1782, lsk_1431, \
                         lsk_1432, lsk_1433, lsk_1434, lsl1_1782, msk_1432, msk_1433, \
                         msk_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = pa_x[k] * lsl0_1782[k]
                    + f_16 * lsk_1431[k]
                    - f_14 * pc_x[k] * lsl1_1782[k];

        t_1783[k] = f_15 * lsk_1432[k]
                    + f_3 * pc_x[k] * msk_1432[k];

        t_1784[k] = f_15 * lsk_1433[k]
                    + f_3 * pc_x[k] * msk_1433[k];

        t_1785[k] = f_15 * lsk_1434[k]
                    + f_3 * pc_x[k] * msk_1434[k];
    }

#pragma omp simd aligned(t_1786, t_1787, t_1788, t_1789, t_1790, pc_x, lsk_1435, lsk_1436, \
                         lsk_1437, lsk_1438, lsk_1439, msk_1435, msk_1436, msk_1437, msk_1438, \
                         msk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1786[k] = f_15 * lsk_1435[k]
                    + f_3 * pc_x[k] * msk_1435[k];

        t_1787[k] = f_15 * lsk_1436[k]
                    + f_3 * pc_x[k] * msk_1436[k];

        t_1788[k] = f_15 * lsk_1437[k]
                    + f_3 * pc_x[k] * msk_1437[k];

        t_1789[k] = f_15 * lsk_1438[k]
                    + f_3 * pc_x[k] * msk_1438[k];

        t_1790[k] = f_15 * lsk_1439[k]
                    + f_3 * pc_x[k] * msk_1439[k];
    }

#pragma omp simd aligned(t_1791, t_1792, t_1793, t_1794, pa_x, pc_x, pc_z, lsl0_1791, \
                         lsl0_1793, lsl0_1794, lsk_1108, lsl1_1791, lsl1_1793, lsl1_1794, \
                         msk_1432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1791[k] = pa_x[k] * lsl0_1791[k]
                    - f_14 * pc_x[k] * lsl1_1791[k];

        t_1792[k] = f_17 * lsk_1108[k]
                    + f_3 * pc_z[k] * msk_1432[k];

        t_1793[k] = pa_x[k] * lsl0_1793[k]
                    - f_14 * pc_x[k] * lsl1_1793[k];

        t_1794[k] = pa_x[k] * lsl0_1794[k]
                    - f_14 * pc_x[k] * lsl1_1794[k];
    }

#pragma omp simd aligned(t_1795, t_1796, t_1797, t_1798, pa_x, pc_x, pc_y, lsl0_1795, \
                         lsl0_1796, lsl0_1797, lsk_1151, lsl1_1795, lsl1_1796, lsl1_1797, \
                         msk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1795[k] = pa_x[k] * lsl0_1795[k]
                    - f_14 * pc_x[k] * lsl1_1795[k];

        t_1796[k] = pa_x[k] * lsl0_1796[k]
                    - f_14 * pc_x[k] * lsl1_1796[k];

        t_1797[k] = pa_x[k] * lsl0_1797[k]
                    - f_14 * pc_x[k] * lsl1_1797[k];

        t_1798[k] = f_19 * lsk_1151[k]
                    + f_3 * pc_y[k] * msk_1439[k];
    }

#pragma omp simd aligned(t_1799, t_1800, t_1801, t_1802, pa_x, pc_x, pc_y, pc_z, lsl0_1799, \
                         lsl0_1800, lsk_1116, lsk_1152, lsk_1440, lsl1_1799, lsl1_1800, \
                         msk_1440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1799[k] = pa_x[k] * lsl0_1799[k]
                    - f_14 * pc_x[k] * lsl1_1799[k];

        t_1800[k] = pa_x[k] * lsl0_1800[k]
                    + f_21 * lsk_1440[k]
                    - f_14 * pc_x[k] * lsl1_1800[k];

        t_1801[k] = f_18 * lsk_1152[k]
                    + f_3 * pc_y[k] * msk_1440[k];

        t_1802[k] = f_18 * lsk_1116[k]
                    + f_3 * pc_z[k] * msk_1440[k];
    }

#pragma omp simd aligned(t_1803, t_1804, t_1805, pa_x, pc_x, pc_y, lsl0_1803, lsl0_1805, \
                         lsk_1154, lsk_1443, lsk_1445, lsl1_1803, lsl1_1805, \
                         msk_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1803[k] = pa_x[k] * lsl0_1803[k]
                    + f_20 * lsk_1443[k]
                    - f_14 * pc_x[k] * lsl1_1803[k];

        t_1804[k] = f_18 * lsk_1154[k]
                    + f_3 * pc_y[k] * msk_1442[k];

        t_1805[k] = pa_x[k] * lsl0_1805[k]
                    + f_20 * lsk_1445[k]
                    - f_14 * pc_x[k] * lsl1_1805[k];
    }

#pragma omp simd aligned(t_1806, t_1807, t_1808, pa_x, pc_x, pc_y, pc_z, lsl0_1806, lsk_1119, \
                         lsk_1157, lsk_1446, lsl1_1806, msk_1443, \
                         msk_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1806[k] = pa_x[k] * lsl0_1806[k]
                    + f_19 * lsk_1446[k]
                    - f_14 * pc_x[k] * lsl1_1806[k];

        t_1807[k] = f_18 * lsk_1119[k]
                    + f_3 * pc_z[k] * msk_1443[k];

        t_1808[k] = f_18 * lsk_1157[k]
                    + f_3 * pc_y[k] * msk_1445[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, pa_x, pc_x, pc_z, lsl0_1809, lsl0_1810, \
                         lsk_1122, lsk_1449, lsk_1450, lsl1_1809, lsl1_1810, \
                         msk_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = pa_x[k] * lsl0_1809[k]
                    + f_19 * lsk_1449[k]
                    - f_14 * pc_x[k] * lsl1_1809[k];

        t_1810[k] = pa_x[k] * lsl0_1810[k]
                    + f_18 * lsk_1450[k]
                    - f_14 * pc_x[k] * lsl1_1810[k];

        t_1811[k] = f_18 * lsk_1122[k]
                    + f_3 * pc_z[k] * msk_1446[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msk, const size_t ncols,
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
    const auto f_21 = 4.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1812 = buffer.data(lsl0 + 1812);
    const auto *lsl0_1814 = buffer.data(lsl0 + 1814);
    const auto *lsl0_1815 = buffer.data(lsl0 + 1815);
    const auto *lsl0_1817 = buffer.data(lsl0 + 1817);
    const auto *lsl0_1818 = buffer.data(lsl0 + 1818);
    const auto *lsl0_1820 = buffer.data(lsl0 + 1820);
    const auto *lsl0_1821 = buffer.data(lsl0 + 1821);
    const auto *lsl0_1823 = buffer.data(lsl0 + 1823);
    const auto *lsl0_1824 = buffer.data(lsl0 + 1824);
    const auto *lsl0_1825 = buffer.data(lsl0 + 1825);
    const auto *lsl0_1827 = buffer.data(lsl0 + 1827);
    const auto *lsl0_1836 = buffer.data(lsl0 + 1836);
    const auto *lsl0_1838 = buffer.data(lsl0 + 1838);
    const auto *lsl0_1839 = buffer.data(lsl0 + 1839);
    const auto *lsl0_1840 = buffer.data(lsl0 + 1840);
    const auto *lsl0_1841 = buffer.data(lsl0 + 1841);
    const auto *lsl0_1842 = buffer.data(lsl0 + 1842);
    const auto *lsl0_1844 = buffer.data(lsl0 + 1844);
    const auto *lsl0_1845 = buffer.data(lsl0 + 1845);
    const auto *lsl0_1848 = buffer.data(lsl0 + 1848);
    const auto *lsl0_1850 = buffer.data(lsl0 + 1850);
    const auto *lsl0_1851 = buffer.data(lsl0 + 1851);
    const auto *lsl0_1854 = buffer.data(lsl0 + 1854);
    const auto *lsl0_1855 = buffer.data(lsl0 + 1855);
    const auto *lsl0_1857 = buffer.data(lsl0 + 1857);
    const auto *lsl0_1859 = buffer.data(lsl0 + 1859);
    const auto *lsl0_1860 = buffer.data(lsl0 + 1860);
    const auto *lsl0_1862 = buffer.data(lsl0 + 1862);
    const auto *lsl0_1863 = buffer.data(lsl0 + 1863);
    const auto *lsl0_1865 = buffer.data(lsl0 + 1865);
    const auto *lsl0_1866 = buffer.data(lsl0 + 1866);
    const auto *lsl0_1868 = buffer.data(lsl0 + 1868);
    const auto *lsl0_1869 = buffer.data(lsl0 + 1869);
    const auto *lsl0_1870 = buffer.data(lsl0 + 1870);
    const auto *lsl0_1872 = buffer.data(lsl0 + 1872);
    const auto *lsl0_1881 = buffer.data(lsl0 + 1881);
    const auto *lsl0_1883 = buffer.data(lsl0 + 1883);
    const auto *lsl0_1884 = buffer.data(lsl0 + 1884);
    const auto *lsl0_1885 = buffer.data(lsl0 + 1885);
    const auto *lsl0_1886 = buffer.data(lsl0 + 1886);
    const auto *lsl0_1887 = buffer.data(lsl0 + 1887);
    const auto *lsl0_1889 = buffer.data(lsl0 + 1889);
    const auto *lsl0_1890 = buffer.data(lsl0 + 1890);
    const auto *lsl0_1893 = buffer.data(lsl0 + 1893);
    const auto *lsl0_1895 = buffer.data(lsl0 + 1895);
    const auto *lsl0_1896 = buffer.data(lsl0 + 1896);
    const auto *lsl0_1899 = buffer.data(lsl0 + 1899);
    const auto *lsl0_1900 = buffer.data(lsl0 + 1900);
    const auto *lsl0_1902 = buffer.data(lsl0 + 1902);
    const auto *lsl0_1904 = buffer.data(lsl0 + 1904);
    const auto *lsl0_1905 = buffer.data(lsl0 + 1905);
    const auto *lsl0_1907 = buffer.data(lsl0 + 1907);
    const auto *lsl0_1908 = buffer.data(lsl0 + 1908);
    const auto *lsl0_1910 = buffer.data(lsl0 + 1910);
    const auto *lsl0_1911 = buffer.data(lsl0 + 1911);
    const auto *lsl0_1913 = buffer.data(lsl0 + 1913);
    const auto *lsl0_1914 = buffer.data(lsl0 + 1914);
    const auto *lsl0_1915 = buffer.data(lsl0 + 1915);
    const auto *lsl0_1917 = buffer.data(lsl0 + 1917);
    const auto *lsl0_1926 = buffer.data(lsl0 + 1926);
    const auto *lsl0_1928 = buffer.data(lsl0 + 1928);
    const auto *lsl0_1929 = buffer.data(lsl0 + 1929);

    const auto *lsk_1126 = buffer.data(lsk + 1126);
    const auto *lsk_1131 = buffer.data(lsk + 1131);
    const auto *lsk_1144 = buffer.data(lsk + 1144);
    const auto *lsk_1152 = buffer.data(lsk + 1152);
    const auto *lsk_1155 = buffer.data(lsk + 1155);
    const auto *lsk_1158 = buffer.data(lsk + 1158);
    const auto *lsk_1161 = buffer.data(lsk + 1161);
    const auto *lsk_1162 = buffer.data(lsk + 1162);
    const auto *lsk_1166 = buffer.data(lsk + 1166);
    const auto *lsk_1167 = buffer.data(lsk + 1167);
    const auto *lsk_1172 = buffer.data(lsk + 1172);
    const auto *lsk_1180 = buffer.data(lsk + 1180);
    const auto *lsk_1187 = buffer.data(lsk + 1187);
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
    const auto *lsk_1223 = buffer.data(lsk + 1223);
    const auto *lsk_1224 = buffer.data(lsk + 1224);
    const auto *lsk_1226 = buffer.data(lsk + 1226);
    const auto *lsk_1229 = buffer.data(lsk + 1229);
    const auto *lsk_1233 = buffer.data(lsk + 1233);
    const auto *lsk_1238 = buffer.data(lsk + 1238);
    const auto *lsk_1244 = buffer.data(lsk + 1244);
    const auto *lsk_1452 = buffer.data(lsk + 1452);
    const auto *lsk_1454 = buffer.data(lsk + 1454);
    const auto *lsk_1455 = buffer.data(lsk + 1455);
    const auto *lsk_1457 = buffer.data(lsk + 1457);
    const auto *lsk_1458 = buffer.data(lsk + 1458);
    const auto *lsk_1460 = buffer.data(lsk + 1460);
    const auto *lsk_1461 = buffer.data(lsk + 1461);
    const auto *lsk_1463 = buffer.data(lsk + 1463);
    const auto *lsk_1464 = buffer.data(lsk + 1464);
    const auto *lsk_1465 = buffer.data(lsk + 1465);
    const auto *lsk_1467 = buffer.data(lsk + 1467);
    const auto *lsk_1468 = buffer.data(lsk + 1468);
    const auto *lsk_1469 = buffer.data(lsk + 1469);
    const auto *lsk_1470 = buffer.data(lsk + 1470);
    const auto *lsk_1471 = buffer.data(lsk + 1471);
    const auto *lsk_1472 = buffer.data(lsk + 1472);
    const auto *lsk_1473 = buffer.data(lsk + 1473);
    const auto *lsk_1474 = buffer.data(lsk + 1474);
    const auto *lsk_1475 = buffer.data(lsk + 1475);
    const auto *lsk_1476 = buffer.data(lsk + 1476);
    const auto *lsk_1479 = buffer.data(lsk + 1479);
    const auto *lsk_1481 = buffer.data(lsk + 1481);
    const auto *lsk_1482 = buffer.data(lsk + 1482);
    const auto *lsk_1485 = buffer.data(lsk + 1485);
    const auto *lsk_1486 = buffer.data(lsk + 1486);
    const auto *lsk_1488 = buffer.data(lsk + 1488);
    const auto *lsk_1490 = buffer.data(lsk + 1490);
    const auto *lsk_1491 = buffer.data(lsk + 1491);
    const auto *lsk_1493 = buffer.data(lsk + 1493);
    const auto *lsk_1494 = buffer.data(lsk + 1494);
    const auto *lsk_1496 = buffer.data(lsk + 1496);
    const auto *lsk_1497 = buffer.data(lsk + 1497);
    const auto *lsk_1499 = buffer.data(lsk + 1499);
    const auto *lsk_1500 = buffer.data(lsk + 1500);
    const auto *lsk_1501 = buffer.data(lsk + 1501);
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
    const auto *lsk_1515 = buffer.data(lsk + 1515);
    const auto *lsk_1517 = buffer.data(lsk + 1517);
    const auto *lsk_1518 = buffer.data(lsk + 1518);
    const auto *lsk_1521 = buffer.data(lsk + 1521);
    const auto *lsk_1522 = buffer.data(lsk + 1522);
    const auto *lsk_1524 = buffer.data(lsk + 1524);
    const auto *lsk_1526 = buffer.data(lsk + 1526);
    const auto *lsk_1527 = buffer.data(lsk + 1527);
    const auto *lsk_1529 = buffer.data(lsk + 1529);
    const auto *lsk_1530 = buffer.data(lsk + 1530);
    const auto *lsk_1532 = buffer.data(lsk + 1532);
    const auto *lsk_1533 = buffer.data(lsk + 1533);
    const auto *lsk_1535 = buffer.data(lsk + 1535);
    const auto *lsk_1536 = buffer.data(lsk + 1536);
    const auto *lsk_1537 = buffer.data(lsk + 1537);
    const auto *lsk_1539 = buffer.data(lsk + 1539);
    const auto *lsk_1540 = buffer.data(lsk + 1540);
    const auto *lsk_1541 = buffer.data(lsk + 1541);
    const auto *lsk_1542 = buffer.data(lsk + 1542);
    const auto *lsk_1543 = buffer.data(lsk + 1543);
    const auto *lsk_1544 = buffer.data(lsk + 1544);
    const auto *lsk_1545 = buffer.data(lsk + 1545);
    const auto *lsk_1546 = buffer.data(lsk + 1546);
    const auto *lsk_1547 = buffer.data(lsk + 1547);

    const auto *lsl1_1812 = buffer.data(lsl1 + 1812);
    const auto *lsl1_1814 = buffer.data(lsl1 + 1814);
    const auto *lsl1_1815 = buffer.data(lsl1 + 1815);
    const auto *lsl1_1817 = buffer.data(lsl1 + 1817);
    const auto *lsl1_1818 = buffer.data(lsl1 + 1818);
    const auto *lsl1_1820 = buffer.data(lsl1 + 1820);
    const auto *lsl1_1821 = buffer.data(lsl1 + 1821);
    const auto *lsl1_1823 = buffer.data(lsl1 + 1823);
    const auto *lsl1_1824 = buffer.data(lsl1 + 1824);
    const auto *lsl1_1825 = buffer.data(lsl1 + 1825);
    const auto *lsl1_1827 = buffer.data(lsl1 + 1827);
    const auto *lsl1_1836 = buffer.data(lsl1 + 1836);
    const auto *lsl1_1838 = buffer.data(lsl1 + 1838);
    const auto *lsl1_1839 = buffer.data(lsl1 + 1839);
    const auto *lsl1_1840 = buffer.data(lsl1 + 1840);
    const auto *lsl1_1841 = buffer.data(lsl1 + 1841);
    const auto *lsl1_1842 = buffer.data(lsl1 + 1842);
    const auto *lsl1_1844 = buffer.data(lsl1 + 1844);
    const auto *lsl1_1845 = buffer.data(lsl1 + 1845);
    const auto *lsl1_1848 = buffer.data(lsl1 + 1848);
    const auto *lsl1_1850 = buffer.data(lsl1 + 1850);
    const auto *lsl1_1851 = buffer.data(lsl1 + 1851);
    const auto *lsl1_1854 = buffer.data(lsl1 + 1854);
    const auto *lsl1_1855 = buffer.data(lsl1 + 1855);
    const auto *lsl1_1857 = buffer.data(lsl1 + 1857);
    const auto *lsl1_1859 = buffer.data(lsl1 + 1859);
    const auto *lsl1_1860 = buffer.data(lsl1 + 1860);
    const auto *lsl1_1862 = buffer.data(lsl1 + 1862);
    const auto *lsl1_1863 = buffer.data(lsl1 + 1863);
    const auto *lsl1_1865 = buffer.data(lsl1 + 1865);
    const auto *lsl1_1866 = buffer.data(lsl1 + 1866);
    const auto *lsl1_1868 = buffer.data(lsl1 + 1868);
    const auto *lsl1_1869 = buffer.data(lsl1 + 1869);
    const auto *lsl1_1870 = buffer.data(lsl1 + 1870);
    const auto *lsl1_1872 = buffer.data(lsl1 + 1872);
    const auto *lsl1_1881 = buffer.data(lsl1 + 1881);
    const auto *lsl1_1883 = buffer.data(lsl1 + 1883);
    const auto *lsl1_1884 = buffer.data(lsl1 + 1884);
    const auto *lsl1_1885 = buffer.data(lsl1 + 1885);
    const auto *lsl1_1886 = buffer.data(lsl1 + 1886);
    const auto *lsl1_1887 = buffer.data(lsl1 + 1887);
    const auto *lsl1_1889 = buffer.data(lsl1 + 1889);
    const auto *lsl1_1890 = buffer.data(lsl1 + 1890);
    const auto *lsl1_1893 = buffer.data(lsl1 + 1893);
    const auto *lsl1_1895 = buffer.data(lsl1 + 1895);
    const auto *lsl1_1896 = buffer.data(lsl1 + 1896);
    const auto *lsl1_1899 = buffer.data(lsl1 + 1899);
    const auto *lsl1_1900 = buffer.data(lsl1 + 1900);
    const auto *lsl1_1902 = buffer.data(lsl1 + 1902);
    const auto *lsl1_1904 = buffer.data(lsl1 + 1904);
    const auto *lsl1_1905 = buffer.data(lsl1 + 1905);
    const auto *lsl1_1907 = buffer.data(lsl1 + 1907);
    const auto *lsl1_1908 = buffer.data(lsl1 + 1908);
    const auto *lsl1_1910 = buffer.data(lsl1 + 1910);
    const auto *lsl1_1911 = buffer.data(lsl1 + 1911);
    const auto *lsl1_1913 = buffer.data(lsl1 + 1913);
    const auto *lsl1_1914 = buffer.data(lsl1 + 1914);
    const auto *lsl1_1915 = buffer.data(lsl1 + 1915);
    const auto *lsl1_1917 = buffer.data(lsl1 + 1917);
    const auto *lsl1_1926 = buffer.data(lsl1 + 1926);
    const auto *lsl1_1928 = buffer.data(lsl1 + 1928);
    const auto *lsl1_1929 = buffer.data(lsl1 + 1929);

    const auto *msk_1449 = buffer.data(msk + 1449);
    const auto *msk_1450 = buffer.data(msk + 1450);
    const auto *msk_1454 = buffer.data(msk + 1454);
    const auto *msk_1455 = buffer.data(msk + 1455);
    const auto *msk_1460 = buffer.data(msk + 1460);
    const auto *msk_1468 = buffer.data(msk + 1468);
    const auto *msk_1469 = buffer.data(msk + 1469);
    const auto *msk_1470 = buffer.data(msk + 1470);
    const auto *msk_1471 = buffer.data(msk + 1471);
    const auto *msk_1472 = buffer.data(msk + 1472);
    const auto *msk_1473 = buffer.data(msk + 1473);
    const auto *msk_1474 = buffer.data(msk + 1474);
    const auto *msk_1475 = buffer.data(msk + 1475);
    const auto *msk_1476 = buffer.data(msk + 1476);
    const auto *msk_1478 = buffer.data(msk + 1478);
    const auto *msk_1479 = buffer.data(msk + 1479);
    const auto *msk_1481 = buffer.data(msk + 1481);
    const auto *msk_1482 = buffer.data(msk + 1482);
    const auto *msk_1485 = buffer.data(msk + 1485);
    const auto *msk_1486 = buffer.data(msk + 1486);
    const auto *msk_1490 = buffer.data(msk + 1490);
    const auto *msk_1491 = buffer.data(msk + 1491);
    const auto *msk_1496 = buffer.data(msk + 1496);
    const auto *msk_1504 = buffer.data(msk + 1504);
    const auto *msk_1505 = buffer.data(msk + 1505);
    const auto *msk_1506 = buffer.data(msk + 1506);
    const auto *msk_1507 = buffer.data(msk + 1507);
    const auto *msk_1508 = buffer.data(msk + 1508);
    const auto *msk_1509 = buffer.data(msk + 1509);
    const auto *msk_1510 = buffer.data(msk + 1510);
    const auto *msk_1511 = buffer.data(msk + 1511);
    const auto *msk_1512 = buffer.data(msk + 1512);
    const auto *msk_1514 = buffer.data(msk + 1514);
    const auto *msk_1515 = buffer.data(msk + 1515);
    const auto *msk_1517 = buffer.data(msk + 1517);
    const auto *msk_1518 = buffer.data(msk + 1518);
    const auto *msk_1521 = buffer.data(msk + 1521);
    const auto *msk_1522 = buffer.data(msk + 1522);
    const auto *msk_1526 = buffer.data(msk + 1526);
    const auto *msk_1527 = buffer.data(msk + 1527);
    const auto *msk_1532 = buffer.data(msk + 1532);
    const auto *msk_1540 = buffer.data(msk + 1540);
    const auto *msk_1541 = buffer.data(msk + 1541);
    const auto *msk_1542 = buffer.data(msk + 1542);
    const auto *msk_1543 = buffer.data(msk + 1543);
    const auto *msk_1544 = buffer.data(msk + 1544);
    const auto *msk_1545 = buffer.data(msk + 1545);
    const auto *msk_1546 = buffer.data(msk + 1546);
    const auto *msk_1547 = buffer.data(msk + 1547);

#pragma omp simd aligned(t_1812, t_1813, t_1814, pa_x, pc_x, pc_y, lsl0_1812, lsl0_1814, \
                         lsk_1161, lsk_1452, lsk_1454, lsl1_1812, lsl1_1814, \
                         msk_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1812[k] = pa_x[k] * lsl0_1812[k]
                    + f_18 * lsk_1452[k]
                    - f_14 * pc_x[k] * lsl1_1812[k];

        t_1813[k] = f_18 * lsk_1161[k]
                    + f_3 * pc_y[k] * msk_1449[k];

        t_1814[k] = pa_x[k] * lsl0_1814[k]
                    + f_18 * lsk_1454[k]
                    - f_14 * pc_x[k] * lsl1_1814[k];
    }

#pragma omp simd aligned(t_1815, t_1816, t_1817, pa_x, pc_x, pc_z, lsl0_1815, lsl0_1817, \
                         lsk_1126, lsk_1455, lsk_1457, lsl1_1815, lsl1_1817, \
                         msk_1450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1815[k] = pa_x[k] * lsl0_1815[k]
                    + f_17 * lsk_1455[k]
                    - f_14 * pc_x[k] * lsl1_1815[k];

        t_1816[k] = f_18 * lsk_1126[k]
                    + f_3 * pc_z[k] * msk_1450[k];

        t_1817[k] = pa_x[k] * lsl0_1817[k]
                    + f_17 * lsk_1457[k]
                    - f_14 * pc_x[k] * lsl1_1817[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, pa_x, pc_x, pc_y, lsl0_1818, lsl0_1820, \
                         lsk_1166, lsk_1458, lsk_1460, lsl1_1818, lsl1_1820, \
                         msk_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = pa_x[k] * lsl0_1818[k]
                    + f_17 * lsk_1458[k]
                    - f_14 * pc_x[k] * lsl1_1818[k];

        t_1819[k] = f_18 * lsk_1166[k]
                    + f_3 * pc_y[k] * msk_1454[k];

        t_1820[k] = pa_x[k] * lsl0_1820[k]
                    + f_17 * lsk_1460[k]
                    - f_14 * pc_x[k] * lsl1_1820[k];
    }

#pragma omp simd aligned(t_1821, t_1822, t_1823, pa_x, pc_x, pc_z, lsl0_1821, lsl0_1823, \
                         lsk_1131, lsk_1461, lsk_1463, lsl1_1821, lsl1_1823, \
                         msk_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1821[k] = pa_x[k] * lsl0_1821[k]
                    + f_16 * lsk_1461[k]
                    - f_14 * pc_x[k] * lsl1_1821[k];

        t_1822[k] = f_18 * lsk_1131[k]
                    + f_3 * pc_z[k] * msk_1455[k];

        t_1823[k] = pa_x[k] * lsl0_1823[k]
                    + f_16 * lsk_1463[k]
                    - f_14 * pc_x[k] * lsl1_1823[k];
    }

#pragma omp simd aligned(t_1824, t_1825, t_1826, pa_x, pc_x, pc_y, lsl0_1824, lsl0_1825, \
                         lsk_1172, lsk_1464, lsk_1465, lsl1_1824, lsl1_1825, \
                         msk_1460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1824[k] = pa_x[k] * lsl0_1824[k]
                    + f_16 * lsk_1464[k]
                    - f_14 * pc_x[k] * lsl1_1824[k];

        t_1825[k] = pa_x[k] * lsl0_1825[k]
                    + f_16 * lsk_1465[k]
                    - f_14 * pc_x[k] * lsl1_1825[k];

        t_1826[k] = f_18 * lsk_1172[k]
                    + f_3 * pc_y[k] * msk_1460[k];
    }

#pragma omp simd aligned(t_1827, t_1828, t_1829, t_1830, pa_x, pc_x, lsl0_1827, lsk_1467, \
                         lsk_1468, lsk_1469, lsk_1470, lsl1_1827, msk_1468, msk_1469, \
                         msk_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1827[k] = pa_x[k] * lsl0_1827[k]
                    + f_16 * lsk_1467[k]
                    - f_14 * pc_x[k] * lsl1_1827[k];

        t_1828[k] = f_15 * lsk_1468[k]
                    + f_3 * pc_x[k] * msk_1468[k];

        t_1829[k] = f_15 * lsk_1469[k]
                    + f_3 * pc_x[k] * msk_1469[k];

        t_1830[k] = f_15 * lsk_1470[k]
                    + f_3 * pc_x[k] * msk_1470[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, t_1834, t_1835, pc_x, lsk_1471, lsk_1472, \
                         lsk_1473, lsk_1474, lsk_1475, msk_1471, msk_1472, msk_1473, msk_1474, \
                         msk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_15 * lsk_1471[k]
                    + f_3 * pc_x[k] * msk_1471[k];

        t_1832[k] = f_15 * lsk_1472[k]
                    + f_3 * pc_x[k] * msk_1472[k];

        t_1833[k] = f_15 * lsk_1473[k]
                    + f_3 * pc_x[k] * msk_1473[k];

        t_1834[k] = f_15 * lsk_1474[k]
                    + f_3 * pc_x[k] * msk_1474[k];

        t_1835[k] = f_15 * lsk_1475[k]
                    + f_3 * pc_x[k] * msk_1475[k];
    }

#pragma omp simd aligned(t_1836, t_1837, t_1838, t_1839, pa_x, pc_x, pc_z, lsl0_1836, \
                         lsl0_1838, lsl0_1839, lsk_1144, lsl1_1836, lsl1_1838, lsl1_1839, \
                         msk_1468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1836[k] = pa_x[k] * lsl0_1836[k]
                    - f_14 * pc_x[k] * lsl1_1836[k];

        t_1837[k] = f_18 * lsk_1144[k]
                    + f_3 * pc_z[k] * msk_1468[k];

        t_1838[k] = pa_x[k] * lsl0_1838[k]
                    - f_14 * pc_x[k] * lsl1_1838[k];

        t_1839[k] = pa_x[k] * lsl0_1839[k]
                    - f_14 * pc_x[k] * lsl1_1839[k];
    }

#pragma omp simd aligned(t_1840, t_1841, t_1842, t_1843, pa_x, pc_x, pc_y, lsl0_1840, \
                         lsl0_1841, lsl0_1842, lsk_1187, lsl1_1840, lsl1_1841, lsl1_1842, \
                         msk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1840[k] = pa_x[k] * lsl0_1840[k]
                    - f_14 * pc_x[k] * lsl1_1840[k];

        t_1841[k] = pa_x[k] * lsl0_1841[k]
                    - f_14 * pc_x[k] * lsl1_1841[k];

        t_1842[k] = pa_x[k] * lsl0_1842[k]
                    - f_14 * pc_x[k] * lsl1_1842[k];

        t_1843[k] = f_18 * lsk_1187[k]
                    + f_3 * pc_y[k] * msk_1475[k];
    }

#pragma omp simd aligned(t_1844, t_1845, t_1846, t_1847, pa_x, pc_x, pc_y, pc_z, lsl0_1844, \
                         lsl0_1845, lsk_1152, lsk_1188, lsk_1476, lsl1_1844, lsl1_1845, \
                         msk_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1844[k] = pa_x[k] * lsl0_1844[k]
                    - f_14 * pc_x[k] * lsl1_1844[k];

        t_1845[k] = pa_x[k] * lsl0_1845[k]
                    + f_21 * lsk_1476[k]
                    - f_14 * pc_x[k] * lsl1_1845[k];

        t_1846[k] = f_17 * lsk_1188[k]
                    + f_3 * pc_y[k] * msk_1476[k];

        t_1847[k] = f_19 * lsk_1152[k]
                    + f_3 * pc_z[k] * msk_1476[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, pa_x, pc_x, pc_y, lsl0_1848, lsl0_1850, \
                         lsk_1190, lsk_1479, lsk_1481, lsl1_1848, lsl1_1850, \
                         msk_1478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = pa_x[k] * lsl0_1848[k]
                    + f_20 * lsk_1479[k]
                    - f_14 * pc_x[k] * lsl1_1848[k];

        t_1849[k] = f_17 * lsk_1190[k]
                    + f_3 * pc_y[k] * msk_1478[k];

        t_1850[k] = pa_x[k] * lsl0_1850[k]
                    + f_20 * lsk_1481[k]
                    - f_14 * pc_x[k] * lsl1_1850[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pa_x, pc_x, pc_y, pc_z, lsl0_1851, lsk_1155, \
                         lsk_1193, lsk_1482, lsl1_1851, msk_1479, \
                         msk_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = pa_x[k] * lsl0_1851[k]
                    + f_19 * lsk_1482[k]
                    - f_14 * pc_x[k] * lsl1_1851[k];

        t_1852[k] = f_19 * lsk_1155[k]
                    + f_3 * pc_z[k] * msk_1479[k];

        t_1853[k] = f_17 * lsk_1193[k]
                    + f_3 * pc_y[k] * msk_1481[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, pa_x, pc_x, pc_z, lsl0_1854, lsl0_1855, \
                         lsk_1158, lsk_1485, lsk_1486, lsl1_1854, lsl1_1855, \
                         msk_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = pa_x[k] * lsl0_1854[k]
                    + f_19 * lsk_1485[k]
                    - f_14 * pc_x[k] * lsl1_1854[k];

        t_1855[k] = pa_x[k] * lsl0_1855[k]
                    + f_18 * lsk_1486[k]
                    - f_14 * pc_x[k] * lsl1_1855[k];

        t_1856[k] = f_19 * lsk_1158[k]
                    + f_3 * pc_z[k] * msk_1482[k];
    }

#pragma omp simd aligned(t_1857, t_1858, t_1859, pa_x, pc_x, pc_y, lsl0_1857, lsl0_1859, \
                         lsk_1197, lsk_1488, lsk_1490, lsl1_1857, lsl1_1859, \
                         msk_1485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1857[k] = pa_x[k] * lsl0_1857[k]
                    + f_18 * lsk_1488[k]
                    - f_14 * pc_x[k] * lsl1_1857[k];

        t_1858[k] = f_17 * lsk_1197[k]
                    + f_3 * pc_y[k] * msk_1485[k];

        t_1859[k] = pa_x[k] * lsl0_1859[k]
                    + f_18 * lsk_1490[k]
                    - f_14 * pc_x[k] * lsl1_1859[k];
    }

#pragma omp simd aligned(t_1860, t_1861, t_1862, pa_x, pc_x, pc_z, lsl0_1860, lsl0_1862, \
                         lsk_1162, lsk_1491, lsk_1493, lsl1_1860, lsl1_1862, \
                         msk_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1860[k] = pa_x[k] * lsl0_1860[k]
                    + f_17 * lsk_1491[k]
                    - f_14 * pc_x[k] * lsl1_1860[k];

        t_1861[k] = f_19 * lsk_1162[k]
                    + f_3 * pc_z[k] * msk_1486[k];

        t_1862[k] = pa_x[k] * lsl0_1862[k]
                    + f_17 * lsk_1493[k]
                    - f_14 * pc_x[k] * lsl1_1862[k];
    }

#pragma omp simd aligned(t_1863, t_1864, t_1865, pa_x, pc_x, pc_y, lsl0_1863, lsl0_1865, \
                         lsk_1202, lsk_1494, lsk_1496, lsl1_1863, lsl1_1865, \
                         msk_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1863[k] = pa_x[k] * lsl0_1863[k]
                    + f_17 * lsk_1494[k]
                    - f_14 * pc_x[k] * lsl1_1863[k];

        t_1864[k] = f_17 * lsk_1202[k]
                    + f_3 * pc_y[k] * msk_1490[k];

        t_1865[k] = pa_x[k] * lsl0_1865[k]
                    + f_17 * lsk_1496[k]
                    - f_14 * pc_x[k] * lsl1_1865[k];
    }

#pragma omp simd aligned(t_1866, t_1867, t_1868, pa_x, pc_x, pc_z, lsl0_1866, lsl0_1868, \
                         lsk_1167, lsk_1497, lsk_1499, lsl1_1866, lsl1_1868, \
                         msk_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1866[k] = pa_x[k] * lsl0_1866[k]
                    + f_16 * lsk_1497[k]
                    - f_14 * pc_x[k] * lsl1_1866[k];

        t_1867[k] = f_19 * lsk_1167[k]
                    + f_3 * pc_z[k] * msk_1491[k];

        t_1868[k] = pa_x[k] * lsl0_1868[k]
                    + f_16 * lsk_1499[k]
                    - f_14 * pc_x[k] * lsl1_1868[k];
    }

#pragma omp simd aligned(t_1869, t_1870, t_1871, pa_x, pc_x, pc_y, lsl0_1869, lsl0_1870, \
                         lsk_1208, lsk_1500, lsk_1501, lsl1_1869, lsl1_1870, \
                         msk_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1869[k] = pa_x[k] * lsl0_1869[k]
                    + f_16 * lsk_1500[k]
                    - f_14 * pc_x[k] * lsl1_1869[k];

        t_1870[k] = pa_x[k] * lsl0_1870[k]
                    + f_16 * lsk_1501[k]
                    - f_14 * pc_x[k] * lsl1_1870[k];

        t_1871[k] = f_17 * lsk_1208[k]
                    + f_3 * pc_y[k] * msk_1496[k];
    }

#pragma omp simd aligned(t_1872, t_1873, t_1874, t_1875, pa_x, pc_x, lsl0_1872, lsk_1503, \
                         lsk_1504, lsk_1505, lsk_1506, lsl1_1872, msk_1504, msk_1505, \
                         msk_1506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1872[k] = pa_x[k] * lsl0_1872[k]
                    + f_16 * lsk_1503[k]
                    - f_14 * pc_x[k] * lsl1_1872[k];

        t_1873[k] = f_15 * lsk_1504[k]
                    + f_3 * pc_x[k] * msk_1504[k];

        t_1874[k] = f_15 * lsk_1505[k]
                    + f_3 * pc_x[k] * msk_1505[k];

        t_1875[k] = f_15 * lsk_1506[k]
                    + f_3 * pc_x[k] * msk_1506[k];
    }

#pragma omp simd aligned(t_1876, t_1877, t_1878, t_1879, t_1880, pc_x, lsk_1507, lsk_1508, \
                         lsk_1509, lsk_1510, lsk_1511, msk_1507, msk_1508, msk_1509, msk_1510, \
                         msk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1876[k] = f_15 * lsk_1507[k]
                    + f_3 * pc_x[k] * msk_1507[k];

        t_1877[k] = f_15 * lsk_1508[k]
                    + f_3 * pc_x[k] * msk_1508[k];

        t_1878[k] = f_15 * lsk_1509[k]
                    + f_3 * pc_x[k] * msk_1509[k];

        t_1879[k] = f_15 * lsk_1510[k]
                    + f_3 * pc_x[k] * msk_1510[k];

        t_1880[k] = f_15 * lsk_1511[k]
                    + f_3 * pc_x[k] * msk_1511[k];
    }

#pragma omp simd aligned(t_1881, t_1882, t_1883, t_1884, pa_x, pc_x, pc_z, lsl0_1881, \
                         lsl0_1883, lsl0_1884, lsk_1180, lsl1_1881, lsl1_1883, lsl1_1884, \
                         msk_1504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1881[k] = pa_x[k] * lsl0_1881[k]
                    - f_14 * pc_x[k] * lsl1_1881[k];

        t_1882[k] = f_19 * lsk_1180[k]
                    + f_3 * pc_z[k] * msk_1504[k];

        t_1883[k] = pa_x[k] * lsl0_1883[k]
                    - f_14 * pc_x[k] * lsl1_1883[k];

        t_1884[k] = pa_x[k] * lsl0_1884[k]
                    - f_14 * pc_x[k] * lsl1_1884[k];
    }

#pragma omp simd aligned(t_1885, t_1886, t_1887, t_1888, pa_x, pc_x, pc_y, lsl0_1885, \
                         lsl0_1886, lsl0_1887, lsk_1223, lsl1_1885, lsl1_1886, lsl1_1887, \
                         msk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1885[k] = pa_x[k] * lsl0_1885[k]
                    - f_14 * pc_x[k] * lsl1_1885[k];

        t_1886[k] = pa_x[k] * lsl0_1886[k]
                    - f_14 * pc_x[k] * lsl1_1886[k];

        t_1887[k] = pa_x[k] * lsl0_1887[k]
                    - f_14 * pc_x[k] * lsl1_1887[k];

        t_1888[k] = f_17 * lsk_1223[k]
                    + f_3 * pc_y[k] * msk_1511[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, t_1892, pa_x, pc_x, pc_y, pc_z, lsl0_1889, \
                         lsl0_1890, lsk_1188, lsk_1224, lsk_1512, lsl1_1889, lsl1_1890, \
                         msk_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = pa_x[k] * lsl0_1889[k]
                    - f_14 * pc_x[k] * lsl1_1889[k];

        t_1890[k] = pa_x[k] * lsl0_1890[k]
                    + f_21 * lsk_1512[k]
                    - f_14 * pc_x[k] * lsl1_1890[k];

        t_1891[k] = f_16 * lsk_1224[k]
                    + f_3 * pc_y[k] * msk_1512[k];

        t_1892[k] = f_20 * lsk_1188[k]
                    + f_3 * pc_z[k] * msk_1512[k];
    }

#pragma omp simd aligned(t_1893, t_1894, t_1895, pa_x, pc_x, pc_y, lsl0_1893, lsl0_1895, \
                         lsk_1226, lsk_1515, lsk_1517, lsl1_1893, lsl1_1895, \
                         msk_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1893[k] = pa_x[k] * lsl0_1893[k]
                    + f_20 * lsk_1515[k]
                    - f_14 * pc_x[k] * lsl1_1893[k];

        t_1894[k] = f_16 * lsk_1226[k]
                    + f_3 * pc_y[k] * msk_1514[k];

        t_1895[k] = pa_x[k] * lsl0_1895[k]
                    + f_20 * lsk_1517[k]
                    - f_14 * pc_x[k] * lsl1_1895[k];
    }

#pragma omp simd aligned(t_1896, t_1897, t_1898, pa_x, pc_x, pc_y, pc_z, lsl0_1896, lsk_1191, \
                         lsk_1229, lsk_1518, lsl1_1896, msk_1515, \
                         msk_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1896[k] = pa_x[k] * lsl0_1896[k]
                    + f_19 * lsk_1518[k]
                    - f_14 * pc_x[k] * lsl1_1896[k];

        t_1897[k] = f_20 * lsk_1191[k]
                    + f_3 * pc_z[k] * msk_1515[k];

        t_1898[k] = f_16 * lsk_1229[k]
                    + f_3 * pc_y[k] * msk_1517[k];
    }

#pragma omp simd aligned(t_1899, t_1900, t_1901, pa_x, pc_x, pc_z, lsl0_1899, lsl0_1900, \
                         lsk_1194, lsk_1521, lsk_1522, lsl1_1899, lsl1_1900, \
                         msk_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1899[k] = pa_x[k] * lsl0_1899[k]
                    + f_19 * lsk_1521[k]
                    - f_14 * pc_x[k] * lsl1_1899[k];

        t_1900[k] = pa_x[k] * lsl0_1900[k]
                    + f_18 * lsk_1522[k]
                    - f_14 * pc_x[k] * lsl1_1900[k];

        t_1901[k] = f_20 * lsk_1194[k]
                    + f_3 * pc_z[k] * msk_1518[k];
    }

#pragma omp simd aligned(t_1902, t_1903, t_1904, pa_x, pc_x, pc_y, lsl0_1902, lsl0_1904, \
                         lsk_1233, lsk_1524, lsk_1526, lsl1_1902, lsl1_1904, \
                         msk_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1902[k] = pa_x[k] * lsl0_1902[k]
                    + f_18 * lsk_1524[k]
                    - f_14 * pc_x[k] * lsl1_1902[k];

        t_1903[k] = f_16 * lsk_1233[k]
                    + f_3 * pc_y[k] * msk_1521[k];

        t_1904[k] = pa_x[k] * lsl0_1904[k]
                    + f_18 * lsk_1526[k]
                    - f_14 * pc_x[k] * lsl1_1904[k];
    }

#pragma omp simd aligned(t_1905, t_1906, t_1907, pa_x, pc_x, pc_z, lsl0_1905, lsl0_1907, \
                         lsk_1198, lsk_1527, lsk_1529, lsl1_1905, lsl1_1907, \
                         msk_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1905[k] = pa_x[k] * lsl0_1905[k]
                    + f_17 * lsk_1527[k]
                    - f_14 * pc_x[k] * lsl1_1905[k];

        t_1906[k] = f_20 * lsk_1198[k]
                    + f_3 * pc_z[k] * msk_1522[k];

        t_1907[k] = pa_x[k] * lsl0_1907[k]
                    + f_17 * lsk_1529[k]
                    - f_14 * pc_x[k] * lsl1_1907[k];
    }

#pragma omp simd aligned(t_1908, t_1909, t_1910, pa_x, pc_x, pc_y, lsl0_1908, lsl0_1910, \
                         lsk_1238, lsk_1530, lsk_1532, lsl1_1908, lsl1_1910, \
                         msk_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1908[k] = pa_x[k] * lsl0_1908[k]
                    + f_17 * lsk_1530[k]
                    - f_14 * pc_x[k] * lsl1_1908[k];

        t_1909[k] = f_16 * lsk_1238[k]
                    + f_3 * pc_y[k] * msk_1526[k];

        t_1910[k] = pa_x[k] * lsl0_1910[k]
                    + f_17 * lsk_1532[k]
                    - f_14 * pc_x[k] * lsl1_1910[k];
    }

#pragma omp simd aligned(t_1911, t_1912, t_1913, pa_x, pc_x, pc_z, lsl0_1911, lsl0_1913, \
                         lsk_1203, lsk_1533, lsk_1535, lsl1_1911, lsl1_1913, \
                         msk_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1911[k] = pa_x[k] * lsl0_1911[k]
                    + f_16 * lsk_1533[k]
                    - f_14 * pc_x[k] * lsl1_1911[k];

        t_1912[k] = f_20 * lsk_1203[k]
                    + f_3 * pc_z[k] * msk_1527[k];

        t_1913[k] = pa_x[k] * lsl0_1913[k]
                    + f_16 * lsk_1535[k]
                    - f_14 * pc_x[k] * lsl1_1913[k];
    }

#pragma omp simd aligned(t_1914, t_1915, t_1916, pa_x, pc_x, pc_y, lsl0_1914, lsl0_1915, \
                         lsk_1244, lsk_1536, lsk_1537, lsl1_1914, lsl1_1915, \
                         msk_1532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1914[k] = pa_x[k] * lsl0_1914[k]
                    + f_16 * lsk_1536[k]
                    - f_14 * pc_x[k] * lsl1_1914[k];

        t_1915[k] = pa_x[k] * lsl0_1915[k]
                    + f_16 * lsk_1537[k]
                    - f_14 * pc_x[k] * lsl1_1915[k];

        t_1916[k] = f_16 * lsk_1244[k]
                    + f_3 * pc_y[k] * msk_1532[k];
    }

#pragma omp simd aligned(t_1917, t_1918, t_1919, t_1920, pa_x, pc_x, lsl0_1917, lsk_1539, \
                         lsk_1540, lsk_1541, lsk_1542, lsl1_1917, msk_1540, msk_1541, \
                         msk_1542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1917[k] = pa_x[k] * lsl0_1917[k]
                    + f_16 * lsk_1539[k]
                    - f_14 * pc_x[k] * lsl1_1917[k];

        t_1918[k] = f_15 * lsk_1540[k]
                    + f_3 * pc_x[k] * msk_1540[k];

        t_1919[k] = f_15 * lsk_1541[k]
                    + f_3 * pc_x[k] * msk_1541[k];

        t_1920[k] = f_15 * lsk_1542[k]
                    + f_3 * pc_x[k] * msk_1542[k];
    }

#pragma omp simd aligned(t_1921, t_1922, t_1923, t_1924, t_1925, pc_x, lsk_1543, lsk_1544, \
                         lsk_1545, lsk_1546, lsk_1547, msk_1543, msk_1544, msk_1545, msk_1546, \
                         msk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1921[k] = f_15 * lsk_1543[k]
                    + f_3 * pc_x[k] * msk_1543[k];

        t_1922[k] = f_15 * lsk_1544[k]
                    + f_3 * pc_x[k] * msk_1544[k];

        t_1923[k] = f_15 * lsk_1545[k]
                    + f_3 * pc_x[k] * msk_1545[k];

        t_1924[k] = f_15 * lsk_1546[k]
                    + f_3 * pc_x[k] * msk_1546[k];

        t_1925[k] = f_15 * lsk_1547[k]
                    + f_3 * pc_x[k] * msk_1547[k];
    }

#pragma omp simd aligned(t_1926, t_1927, t_1928, t_1929, pa_x, pc_x, pc_z, lsl0_1926, \
                         lsl0_1928, lsl0_1929, lsk_1216, lsl1_1926, lsl1_1928, lsl1_1929, \
                         msk_1540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1926[k] = pa_x[k] * lsl0_1926[k]
                    - f_14 * pc_x[k] * lsl1_1926[k];

        t_1927[k] = f_20 * lsk_1216[k]
                    + f_3 * pc_z[k] * msk_1540[k];

        t_1928[k] = pa_x[k] * lsl0_1928[k]
                    - f_14 * pc_x[k] * lsl1_1928[k];

        t_1929[k] = pa_x[k] * lsl0_1929[k]
                    - f_14 * pc_x[k] * lsl1_1929[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);
    const auto f_24 = 3.5 / q;

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
    auto *t_2049 = buffer.data(target + 2049);
    auto *t_2050 = buffer.data(target + 2050);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1575 = buffer.data(lsl0 + 1575);
    const auto *lsl0_1580 = buffer.data(lsl0 + 1580);
    const auto *lsl0_1584 = buffer.data(lsl0 + 1584);
    const auto *lsl0_1589 = buffer.data(lsl0 + 1589);
    const auto *lsl0_1595 = buffer.data(lsl0 + 1595);
    const auto *lsl0_1602 = buffer.data(lsl0 + 1602);
    const auto *lsl0_1930 = buffer.data(lsl0 + 1930);
    const auto *lsl0_1931 = buffer.data(lsl0 + 1931);
    const auto *lsl0_1932 = buffer.data(lsl0 + 1932);
    const auto *lsl0_1934 = buffer.data(lsl0 + 1934);
    const auto *lsl0_1938 = buffer.data(lsl0 + 1938);
    const auto *lsl0_1941 = buffer.data(lsl0 + 1941);
    const auto *lsl0_1945 = buffer.data(lsl0 + 1945);
    const auto *lsl0_1947 = buffer.data(lsl0 + 1947);
    const auto *lsl0_1950 = buffer.data(lsl0 + 1950);
    const auto *lsl0_1952 = buffer.data(lsl0 + 1952);
    const auto *lsl0_1953 = buffer.data(lsl0 + 1953);
    const auto *lsl0_1956 = buffer.data(lsl0 + 1956);
    const auto *lsl0_1958 = buffer.data(lsl0 + 1958);
    const auto *lsl0_1959 = buffer.data(lsl0 + 1959);
    const auto *lsl0_1960 = buffer.data(lsl0 + 1960);
    const auto *lsl0_1971 = buffer.data(lsl0 + 1971);
    const auto *lsl0_1973 = buffer.data(lsl0 + 1973);
    const auto *lsl0_1974 = buffer.data(lsl0 + 1974);
    const auto *lsl0_1975 = buffer.data(lsl0 + 1975);
    const auto *lsl0_1976 = buffer.data(lsl0 + 1976);
    const auto *lsl0_1977 = buffer.data(lsl0 + 1977);
    const auto *lsl0_1979 = buffer.data(lsl0 + 1979);
    const auto *lsl0_1980 = buffer.data(lsl0 + 1980);
    const auto *lsl0_1985 = buffer.data(lsl0 + 1985);
    const auto *lsl0_1989 = buffer.data(lsl0 + 1989);
    const auto *lsl0_1994 = buffer.data(lsl0 + 1994);
    const auto *lsl0_2000 = buffer.data(lsl0 + 2000);
    const auto *lsl0_2007 = buffer.data(lsl0 + 2007);
    const auto *lsl0_2016 = buffer.data(lsl0 + 2016);
    const auto *lsl0_2017 = buffer.data(lsl0 + 2017);
    const auto *lsl0_2018 = buffer.data(lsl0 + 2018);
    const auto *lsl0_2019 = buffer.data(lsl0 + 2019);
    const auto *lsl0_2020 = buffer.data(lsl0 + 2020);
    const auto *lsl0_2021 = buffer.data(lsl0 + 2021);
    const auto *lsl0_2022 = buffer.data(lsl0 + 2022);
    const auto *lsl0_2024 = buffer.data(lsl0 + 2024);

    const auto *lsk_1224 = buffer.data(lsk + 1224);
    const auto *lsk_1227 = buffer.data(lsk + 1227);
    const auto *lsk_1230 = buffer.data(lsk + 1230);
    const auto *lsk_1234 = buffer.data(lsk + 1234);
    const auto *lsk_1239 = buffer.data(lsk + 1239);
    const auto *lsk_1252 = buffer.data(lsk + 1252);
    const auto *lsk_1259 = buffer.data(lsk + 1259);
    const auto *lsk_1260 = buffer.data(lsk + 1260);
    const auto *lsk_1262 = buffer.data(lsk + 1262);
    const auto *lsk_1265 = buffer.data(lsk + 1265);
    const auto *lsk_1269 = buffer.data(lsk + 1269);
    const auto *lsk_1274 = buffer.data(lsk + 1274);
    const auto *lsk_1280 = buffer.data(lsk + 1280);
    const auto *lsk_1295 = buffer.data(lsk + 1295);
    const auto *lsk_1551 = buffer.data(lsk + 1551);
    const auto *lsk_1554 = buffer.data(lsk + 1554);
    const auto *lsk_1558 = buffer.data(lsk + 1558);
    const auto *lsk_1560 = buffer.data(lsk + 1560);
    const auto *lsk_1563 = buffer.data(lsk + 1563);
    const auto *lsk_1565 = buffer.data(lsk + 1565);
    const auto *lsk_1566 = buffer.data(lsk + 1566);
    const auto *lsk_1569 = buffer.data(lsk + 1569);
    const auto *lsk_1571 = buffer.data(lsk + 1571);
    const auto *lsk_1572 = buffer.data(lsk + 1572);
    const auto *lsk_1573 = buffer.data(lsk + 1573);
    const auto *lsk_1576 = buffer.data(lsk + 1576);
    const auto *lsk_1577 = buffer.data(lsk + 1577);
    const auto *lsk_1578 = buffer.data(lsk + 1578);
    const auto *lsk_1579 = buffer.data(lsk + 1579);
    const auto *lsk_1580 = buffer.data(lsk + 1580);
    const auto *lsk_1581 = buffer.data(lsk + 1581);
    const auto *lsk_1582 = buffer.data(lsk + 1582);
    const auto *lsk_1583 = buffer.data(lsk + 1583);
    const auto *lsk_1584 = buffer.data(lsk + 1584);
    const auto *lsk_1589 = buffer.data(lsk + 1589);
    const auto *lsk_1593 = buffer.data(lsk + 1593);
    const auto *lsk_1598 = buffer.data(lsk + 1598);
    const auto *lsk_1604 = buffer.data(lsk + 1604);
    const auto *lsk_1611 = buffer.data(lsk + 1611);
    const auto *lsk_1612 = buffer.data(lsk + 1612);
    const auto *lsk_1613 = buffer.data(lsk + 1613);
    const auto *lsk_1614 = buffer.data(lsk + 1614);
    const auto *lsk_1615 = buffer.data(lsk + 1615);
    const auto *lsk_1616 = buffer.data(lsk + 1616);
    const auto *lsk_1617 = buffer.data(lsk + 1617);
    const auto *lsk_1619 = buffer.data(lsk + 1619);

    const auto *lsl1_1575 = buffer.data(lsl1 + 1575);
    const auto *lsl1_1580 = buffer.data(lsl1 + 1580);
    const auto *lsl1_1584 = buffer.data(lsl1 + 1584);
    const auto *lsl1_1589 = buffer.data(lsl1 + 1589);
    const auto *lsl1_1595 = buffer.data(lsl1 + 1595);
    const auto *lsl1_1602 = buffer.data(lsl1 + 1602);
    const auto *lsl1_1930 = buffer.data(lsl1 + 1930);
    const auto *lsl1_1931 = buffer.data(lsl1 + 1931);
    const auto *lsl1_1932 = buffer.data(lsl1 + 1932);
    const auto *lsl1_1934 = buffer.data(lsl1 + 1934);
    const auto *lsl1_1938 = buffer.data(lsl1 + 1938);
    const auto *lsl1_1941 = buffer.data(lsl1 + 1941);
    const auto *lsl1_1945 = buffer.data(lsl1 + 1945);
    const auto *lsl1_1947 = buffer.data(lsl1 + 1947);
    const auto *lsl1_1950 = buffer.data(lsl1 + 1950);
    const auto *lsl1_1952 = buffer.data(lsl1 + 1952);
    const auto *lsl1_1953 = buffer.data(lsl1 + 1953);
    const auto *lsl1_1956 = buffer.data(lsl1 + 1956);
    const auto *lsl1_1958 = buffer.data(lsl1 + 1958);
    const auto *lsl1_1959 = buffer.data(lsl1 + 1959);
    const auto *lsl1_1960 = buffer.data(lsl1 + 1960);
    const auto *lsl1_1971 = buffer.data(lsl1 + 1971);
    const auto *lsl1_1973 = buffer.data(lsl1 + 1973);
    const auto *lsl1_1974 = buffer.data(lsl1 + 1974);
    const auto *lsl1_1975 = buffer.data(lsl1 + 1975);
    const auto *lsl1_1976 = buffer.data(lsl1 + 1976);
    const auto *lsl1_1977 = buffer.data(lsl1 + 1977);
    const auto *lsl1_1979 = buffer.data(lsl1 + 1979);
    const auto *lsl1_1980 = buffer.data(lsl1 + 1980);
    const auto *lsl1_1985 = buffer.data(lsl1 + 1985);
    const auto *lsl1_1989 = buffer.data(lsl1 + 1989);
    const auto *lsl1_1994 = buffer.data(lsl1 + 1994);
    const auto *lsl1_2000 = buffer.data(lsl1 + 2000);
    const auto *lsl1_2007 = buffer.data(lsl1 + 2007);
    const auto *lsl1_2016 = buffer.data(lsl1 + 2016);
    const auto *lsl1_2017 = buffer.data(lsl1 + 2017);
    const auto *lsl1_2018 = buffer.data(lsl1 + 2018);
    const auto *lsl1_2019 = buffer.data(lsl1 + 2019);
    const auto *lsl1_2020 = buffer.data(lsl1 + 2020);
    const auto *lsl1_2021 = buffer.data(lsl1 + 2021);
    const auto *lsl1_2022 = buffer.data(lsl1 + 2022);
    const auto *lsl1_2024 = buffer.data(lsl1 + 2024);

    const auto *msi0_1232 = buffer.data(msi0 + 1232);
    const auto *msi0_1233 = buffer.data(msi0 + 1233);
    const auto *msi0_1234 = buffer.data(msi0 + 1234);
    const auto *msi0_1235 = buffer.data(msi0 + 1235);
    const auto *msi0_1236 = buffer.data(msi0 + 1236);
    const auto *msi0_1237 = buffer.data(msi0 + 1237);
    const auto *msi0_1238 = buffer.data(msi0 + 1238);
    const auto *msi0_1239 = buffer.data(msi0 + 1239);
    const auto *msi0_1240 = buffer.data(msi0 + 1240);
    const auto *msi0_1241 = buffer.data(msi0 + 1241);
    const auto *msi0_1242 = buffer.data(msi0 + 1242);
    const auto *msi0_1243 = buffer.data(msi0 + 1243);
    const auto *msi0_1244 = buffer.data(msi0 + 1244);
    const auto *msi0_1245 = buffer.data(msi0 + 1245);
    const auto *msi0_1246 = buffer.data(msi0 + 1246);
    const auto *msi0_1260 = buffer.data(msi0 + 1260);
    const auto *msi0_1261 = buffer.data(msi0 + 1261);
    const auto *msi0_1263 = buffer.data(msi0 + 1263);
    const auto *msi0_1265 = buffer.data(msi0 + 1265);
    const auto *msi0_1266 = buffer.data(msi0 + 1266);
    const auto *msi0_1268 = buffer.data(msi0 + 1268);
    const auto *msi0_1269 = buffer.data(msi0 + 1269);
    const auto *msi0_1270 = buffer.data(msi0 + 1270);
    const auto *msi0_1272 = buffer.data(msi0 + 1272);
    const auto *msi0_1273 = buffer.data(msi0 + 1273);
    const auto *msi0_1274 = buffer.data(msi0 + 1274);
    const auto *msi0_1275 = buffer.data(msi0 + 1275);
    const auto *msi0_1277 = buffer.data(msi0 + 1277);
    const auto *msi0_1278 = buffer.data(msi0 + 1278);
    const auto *msi0_1279 = buffer.data(msi0 + 1279);
    const auto *msi0_1280 = buffer.data(msi0 + 1280);
    const auto *msi0_1281 = buffer.data(msi0 + 1281);
    const auto *msi0_1283 = buffer.data(msi0 + 1283);
    const auto *msi0_1284 = buffer.data(msi0 + 1284);
    const auto *msi0_1285 = buffer.data(msi0 + 1285);

    const auto *msi1_1232 = buffer.data(msi1 + 1232);
    const auto *msi1_1233 = buffer.data(msi1 + 1233);
    const auto *msi1_1234 = buffer.data(msi1 + 1234);
    const auto *msi1_1235 = buffer.data(msi1 + 1235);
    const auto *msi1_1236 = buffer.data(msi1 + 1236);
    const auto *msi1_1237 = buffer.data(msi1 + 1237);
    const auto *msi1_1238 = buffer.data(msi1 + 1238);
    const auto *msi1_1239 = buffer.data(msi1 + 1239);
    const auto *msi1_1240 = buffer.data(msi1 + 1240);
    const auto *msi1_1241 = buffer.data(msi1 + 1241);
    const auto *msi1_1242 = buffer.data(msi1 + 1242);
    const auto *msi1_1243 = buffer.data(msi1 + 1243);
    const auto *msi1_1244 = buffer.data(msi1 + 1244);
    const auto *msi1_1245 = buffer.data(msi1 + 1245);
    const auto *msi1_1246 = buffer.data(msi1 + 1246);
    const auto *msi1_1260 = buffer.data(msi1 + 1260);
    const auto *msi1_1261 = buffer.data(msi1 + 1261);
    const auto *msi1_1263 = buffer.data(msi1 + 1263);
    const auto *msi1_1265 = buffer.data(msi1 + 1265);
    const auto *msi1_1266 = buffer.data(msi1 + 1266);
    const auto *msi1_1268 = buffer.data(msi1 + 1268);
    const auto *msi1_1269 = buffer.data(msi1 + 1269);
    const auto *msi1_1270 = buffer.data(msi1 + 1270);
    const auto *msi1_1272 = buffer.data(msi1 + 1272);
    const auto *msi1_1273 = buffer.data(msi1 + 1273);
    const auto *msi1_1274 = buffer.data(msi1 + 1274);
    const auto *msi1_1275 = buffer.data(msi1 + 1275);
    const auto *msi1_1277 = buffer.data(msi1 + 1277);
    const auto *msi1_1278 = buffer.data(msi1 + 1278);
    const auto *msi1_1279 = buffer.data(msi1 + 1279);
    const auto *msi1_1280 = buffer.data(msi1 + 1280);
    const auto *msi1_1281 = buffer.data(msi1 + 1281);
    const auto *msi1_1283 = buffer.data(msi1 + 1283);
    const auto *msi1_1284 = buffer.data(msi1 + 1284);
    const auto *msi1_1285 = buffer.data(msi1 + 1285);

    const auto *msk_1547 = buffer.data(msk + 1547);
    const auto *msk_1548 = buffer.data(msk + 1548);
    const auto *msk_1550 = buffer.data(msk + 1550);
    const auto *msk_1551 = buffer.data(msk + 1551);
    const auto *msk_1553 = buffer.data(msk + 1553);
    const auto *msk_1554 = buffer.data(msk + 1554);
    const auto *msk_1557 = buffer.data(msk + 1557);
    const auto *msk_1558 = buffer.data(msk + 1558);
    const auto *msk_1562 = buffer.data(msk + 1562);
    const auto *msk_1563 = buffer.data(msk + 1563);
    const auto *msk_1568 = buffer.data(msk + 1568);
    const auto *msk_1576 = buffer.data(msk + 1576);
    const auto *msk_1577 = buffer.data(msk + 1577);
    const auto *msk_1578 = buffer.data(msk + 1578);
    const auto *msk_1579 = buffer.data(msk + 1579);
    const auto *msk_1580 = buffer.data(msk + 1580);
    const auto *msk_1581 = buffer.data(msk + 1581);
    const auto *msk_1582 = buffer.data(msk + 1582);
    const auto *msk_1583 = buffer.data(msk + 1583);
    const auto *msk_1584 = buffer.data(msk + 1584);
    const auto *msk_1585 = buffer.data(msk + 1585);
    const auto *msk_1586 = buffer.data(msk + 1586);
    const auto *msk_1587 = buffer.data(msk + 1587);
    const auto *msk_1588 = buffer.data(msk + 1588);
    const auto *msk_1589 = buffer.data(msk + 1589);
    const auto *msk_1590 = buffer.data(msk + 1590);
    const auto *msk_1591 = buffer.data(msk + 1591);
    const auto *msk_1592 = buffer.data(msk + 1592);
    const auto *msk_1593 = buffer.data(msk + 1593);
    const auto *msk_1594 = buffer.data(msk + 1594);
    const auto *msk_1595 = buffer.data(msk + 1595);
    const auto *msk_1596 = buffer.data(msk + 1596);
    const auto *msk_1597 = buffer.data(msk + 1597);
    const auto *msk_1598 = buffer.data(msk + 1598);
    const auto *msk_1599 = buffer.data(msk + 1599);
    const auto *msk_1600 = buffer.data(msk + 1600);
    const auto *msk_1601 = buffer.data(msk + 1601);
    const auto *msk_1602 = buffer.data(msk + 1602);
    const auto *msk_1603 = buffer.data(msk + 1603);
    const auto *msk_1604 = buffer.data(msk + 1604);
    const auto *msk_1611 = buffer.data(msk + 1611);
    const auto *msk_1612 = buffer.data(msk + 1612);
    const auto *msk_1613 = buffer.data(msk + 1613);
    const auto *msk_1614 = buffer.data(msk + 1614);
    const auto *msk_1615 = buffer.data(msk + 1615);
    const auto *msk_1616 = buffer.data(msk + 1616);
    const auto *msk_1617 = buffer.data(msk + 1617);
    const auto *msk_1619 = buffer.data(msk + 1619);
    const auto *msk_1620 = buffer.data(msk + 1620);
    const auto *msk_1621 = buffer.data(msk + 1621);
    const auto *msk_1623 = buffer.data(msk + 1623);
    const auto *msk_1625 = buffer.data(msk + 1625);
    const auto *msk_1626 = buffer.data(msk + 1626);
    const auto *msk_1628 = buffer.data(msk + 1628);
    const auto *msk_1629 = buffer.data(msk + 1629);
    const auto *msk_1630 = buffer.data(msk + 1630);
    const auto *msk_1632 = buffer.data(msk + 1632);
    const auto *msk_1633 = buffer.data(msk + 1633);
    const auto *msk_1634 = buffer.data(msk + 1634);
    const auto *msk_1635 = buffer.data(msk + 1635);
    const auto *msk_1637 = buffer.data(msk + 1637);
    const auto *msk_1638 = buffer.data(msk + 1638);
    const auto *msk_1639 = buffer.data(msk + 1639);
    const auto *msk_1640 = buffer.data(msk + 1640);
    const auto *msk_1641 = buffer.data(msk + 1641);
    const auto *msk_1643 = buffer.data(msk + 1643);
    const auto *msk_1644 = buffer.data(msk + 1644);
    const auto *msk_1645 = buffer.data(msk + 1645);

#pragma omp simd aligned(t_1930, t_1931, t_1932, t_1933, pa_x, pc_x, pc_y, lsl0_1930, \
                         lsl0_1931, lsl0_1932, lsk_1259, lsl1_1930, lsl1_1931, lsl1_1932, \
                         msk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1930[k] = pa_x[k] * lsl0_1930[k]
                    - f_14 * pc_x[k] * lsl1_1930[k];

        t_1931[k] = pa_x[k] * lsl0_1931[k]
                    - f_14 * pc_x[k] * lsl1_1931[k];

        t_1932[k] = pa_x[k] * lsl0_1932[k]
                    - f_14 * pc_x[k] * lsl1_1932[k];

        t_1933[k] = f_16 * lsk_1259[k]
                    + f_3 * pc_y[k] * msk_1547[k];
    }

#pragma omp simd aligned(t_1934, t_1935, t_1936, t_1937, pa_x, pa_y, pc_x, pc_y, pc_z, \
                         lsl0_1575, lsl0_1934, lsk_1224, lsk_1260, lsl1_1575, lsl1_1934, \
                         msk_1548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1934[k] = pa_x[k] * lsl0_1934[k]
                    - f_14 * pc_x[k] * lsl1_1934[k];

        t_1935[k] = pa_y[k] * lsl0_1575[k]
                    - f_14 * pc_y[k] * lsl1_1575[k];

        t_1936[k] = f_15 * lsk_1260[k]
                    + f_3 * pc_y[k] * msk_1548[k];

        t_1937[k] = f_24 * lsk_1224[k]
                    + f_3 * pc_z[k] * msk_1548[k];
    }

#pragma omp simd aligned(t_1938, t_1939, t_1940, pa_x, pa_y, pc_x, pc_y, lsl0_1580, lsl0_1938, \
                         lsk_1262, lsk_1551, lsl1_1580, lsl1_1938, \
                         msk_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1938[k] = pa_x[k] * lsl0_1938[k]
                    + f_20 * lsk_1551[k]
                    - f_14 * pc_x[k] * lsl1_1938[k];

        t_1939[k] = f_15 * lsk_1262[k]
                    + f_3 * pc_y[k] * msk_1550[k];

        t_1940[k] = pa_y[k] * lsl0_1580[k]
                    - f_14 * pc_y[k] * lsl1_1580[k];
    }

#pragma omp simd aligned(t_1941, t_1942, t_1943, pa_x, pc_x, pc_y, pc_z, lsl0_1941, lsk_1227, \
                         lsk_1265, lsk_1554, lsl1_1941, msk_1551, \
                         msk_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1941[k] = pa_x[k] * lsl0_1941[k]
                    + f_19 * lsk_1554[k]
                    - f_14 * pc_x[k] * lsl1_1941[k];

        t_1942[k] = f_24 * lsk_1227[k]
                    + f_3 * pc_z[k] * msk_1551[k];

        t_1943[k] = f_15 * lsk_1265[k]
                    + f_3 * pc_y[k] * msk_1553[k];
    }

#pragma omp simd aligned(t_1944, t_1945, t_1946, pa_x, pa_y, pc_x, pc_y, pc_z, lsl0_1584, \
                         lsl0_1945, lsk_1230, lsk_1558, lsl1_1584, lsl1_1945, \
                         msk_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1944[k] = pa_y[k] * lsl0_1584[k]
                    - f_14 * pc_y[k] * lsl1_1584[k];

        t_1945[k] = pa_x[k] * lsl0_1945[k]
                    + f_18 * lsk_1558[k]
                    - f_14 * pc_x[k] * lsl1_1945[k];

        t_1946[k] = f_24 * lsk_1230[k]
                    + f_3 * pc_z[k] * msk_1554[k];
    }

#pragma omp simd aligned(t_1947, t_1948, t_1949, pa_x, pa_y, pc_x, pc_y, lsl0_1589, lsl0_1947, \
                         lsk_1269, lsk_1560, lsl1_1589, lsl1_1947, \
                         msk_1557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1947[k] = pa_x[k] * lsl0_1947[k]
                    + f_18 * lsk_1560[k]
                    - f_14 * pc_x[k] * lsl1_1947[k];

        t_1948[k] = f_15 * lsk_1269[k]
                    + f_3 * pc_y[k] * msk_1557[k];

        t_1949[k] = pa_y[k] * lsl0_1589[k]
                    - f_14 * pc_y[k] * lsl1_1589[k];
    }

#pragma omp simd aligned(t_1950, t_1951, t_1952, pa_x, pc_x, pc_z, lsl0_1950, lsl0_1952, \
                         lsk_1234, lsk_1563, lsk_1565, lsl1_1950, lsl1_1952, \
                         msk_1558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1950[k] = pa_x[k] * lsl0_1950[k]
                    + f_17 * lsk_1563[k]
                    - f_14 * pc_x[k] * lsl1_1950[k];

        t_1951[k] = f_24 * lsk_1234[k]
                    + f_3 * pc_z[k] * msk_1558[k];

        t_1952[k] = pa_x[k] * lsl0_1952[k]
                    + f_17 * lsk_1565[k]
                    - f_14 * pc_x[k] * lsl1_1952[k];
    }

#pragma omp simd aligned(t_1953, t_1954, t_1955, pa_x, pa_y, pc_x, pc_y, lsl0_1595, lsl0_1953, \
                         lsk_1274, lsk_1566, lsl1_1595, lsl1_1953, \
                         msk_1562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1953[k] = pa_x[k] * lsl0_1953[k]
                    + f_17 * lsk_1566[k]
                    - f_14 * pc_x[k] * lsl1_1953[k];

        t_1954[k] = f_15 * lsk_1274[k]
                    + f_3 * pc_y[k] * msk_1562[k];

        t_1955[k] = pa_y[k] * lsl0_1595[k]
                    - f_14 * pc_y[k] * lsl1_1595[k];
    }

#pragma omp simd aligned(t_1956, t_1957, t_1958, pa_x, pc_x, pc_z, lsl0_1956, lsl0_1958, \
                         lsk_1239, lsk_1569, lsk_1571, lsl1_1956, lsl1_1958, \
                         msk_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1956[k] = pa_x[k] * lsl0_1956[k]
                    + f_16 * lsk_1569[k]
                    - f_14 * pc_x[k] * lsl1_1956[k];

        t_1957[k] = f_24 * lsk_1239[k]
                    + f_3 * pc_z[k] * msk_1563[k];

        t_1958[k] = pa_x[k] * lsl0_1958[k]
                    + f_16 * lsk_1571[k]
                    - f_14 * pc_x[k] * lsl1_1958[k];
    }

#pragma omp simd aligned(t_1959, t_1960, t_1961, pa_x, pc_x, pc_y, lsl0_1959, lsl0_1960, \
                         lsk_1280, lsk_1572, lsk_1573, lsl1_1959, lsl1_1960, \
                         msk_1568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1959[k] = pa_x[k] * lsl0_1959[k]
                    + f_16 * lsk_1572[k]
                    - f_14 * pc_x[k] * lsl1_1959[k];

        t_1960[k] = pa_x[k] * lsl0_1960[k]
                    + f_16 * lsk_1573[k]
                    - f_14 * pc_x[k] * lsl1_1960[k];

        t_1961[k] = f_15 * lsk_1280[k]
                    + f_3 * pc_y[k] * msk_1568[k];
    }

#pragma omp simd aligned(t_1962, t_1963, t_1964, t_1965, pa_y, pc_x, pc_y, lsl0_1602, \
                         lsk_1576, lsk_1577, lsk_1578, lsl1_1602, msk_1576, msk_1577, \
                         msk_1578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1962[k] = pa_y[k] * lsl0_1602[k]
                    - f_14 * pc_y[k] * lsl1_1602[k];

        t_1963[k] = f_15 * lsk_1576[k]
                    + f_3 * pc_x[k] * msk_1576[k];

        t_1964[k] = f_15 * lsk_1577[k]
                    + f_3 * pc_x[k] * msk_1577[k];

        t_1965[k] = f_15 * lsk_1578[k]
                    + f_3 * pc_x[k] * msk_1578[k];
    }

#pragma omp simd aligned(t_1966, t_1967, t_1968, t_1969, t_1970, pc_x, lsk_1579, lsk_1580, \
                         lsk_1581, lsk_1582, lsk_1583, msk_1579, msk_1580, msk_1581, msk_1582, \
                         msk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1966[k] = f_15 * lsk_1579[k]
                    + f_3 * pc_x[k] * msk_1579[k];

        t_1967[k] = f_15 * lsk_1580[k]
                    + f_3 * pc_x[k] * msk_1580[k];

        t_1968[k] = f_15 * lsk_1581[k]
                    + f_3 * pc_x[k] * msk_1581[k];

        t_1969[k] = f_15 * lsk_1582[k]
                    + f_3 * pc_x[k] * msk_1582[k];

        t_1970[k] = f_15 * lsk_1583[k]
                    + f_3 * pc_x[k] * msk_1583[k];
    }

#pragma omp simd aligned(t_1971, t_1972, t_1973, t_1974, pa_x, pc_x, pc_z, lsl0_1971, \
                         lsl0_1973, lsl0_1974, lsk_1252, lsl1_1971, lsl1_1973, lsl1_1974, \
                         msk_1576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1971[k] = pa_x[k] * lsl0_1971[k]
                    - f_14 * pc_x[k] * lsl1_1971[k];

        t_1972[k] = f_24 * lsk_1252[k]
                    + f_3 * pc_z[k] * msk_1576[k];

        t_1973[k] = pa_x[k] * lsl0_1973[k]
                    - f_14 * pc_x[k] * lsl1_1973[k];

        t_1974[k] = pa_x[k] * lsl0_1974[k]
                    - f_14 * pc_x[k] * lsl1_1974[k];
    }

#pragma omp simd aligned(t_1975, t_1976, t_1977, t_1978, pa_x, pc_x, pc_y, lsl0_1975, \
                         lsl0_1976, lsl0_1977, lsk_1295, lsl1_1975, lsl1_1976, lsl1_1977, \
                         msk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1975[k] = pa_x[k] * lsl0_1975[k]
                    - f_14 * pc_x[k] * lsl1_1975[k];

        t_1976[k] = pa_x[k] * lsl0_1976[k]
                    - f_14 * pc_x[k] * lsl1_1976[k];

        t_1977[k] = pa_x[k] * lsl0_1977[k]
                    - f_14 * pc_x[k] * lsl1_1977[k];

        t_1978[k] = f_15 * lsk_1295[k]
                    + f_3 * pc_y[k] * msk_1583[k];
    }

#pragma omp simd aligned(t_1979, t_1980, t_1981, t_1982, pa_x, pc_x, pc_y, pc_z, lsl0_1979, \
                         lsl0_1980, lsk_1260, lsk_1584, lsl1_1979, lsl1_1980, \
                         msk_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1979[k] = pa_x[k] * lsl0_1979[k]
                    - f_14 * pc_x[k] * lsl1_1979[k];

        t_1980[k] = pa_x[k] * lsl0_1980[k]
                    + f_21 * lsk_1584[k]
                    - f_14 * pc_x[k] * lsl1_1980[k];

        t_1981[k] = f_3 * pc_y[k] * msk_1584[k];

        t_1982[k] = f_21 * lsk_1260[k]
                    + f_3 * pc_z[k] * msk_1584[k];
    }

#pragma omp simd aligned(t_1983, t_1984, t_1985, pa_x, pc_x, pc_y, lsl0_1985, lsk_1589, \
                         lsl1_1985, msi0_1232, msi1_1232, msk_1585, \
                         msk_1586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1983[k] = f_4 * msi0_1232[k]
                    - f_5 * msi1_1232[k]
                    + f_3 * pc_y[k] * msk_1585[k];

        t_1984[k] = f_3 * pc_y[k] * msk_1586[k];

        t_1985[k] = pa_x[k] * lsl0_1985[k]
                    + f_20 * lsk_1589[k]
                    - f_14 * pc_x[k] * lsl1_1985[k];
    }

#pragma omp simd aligned(t_1986, t_1987, t_1988, pc_y, msi0_1233, msi0_1234, msi1_1233, \
                         msi1_1234, msk_1587, msk_1588, msk_1589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1986[k] = f_6 * msi0_1233[k]
                    - f_7 * msi1_1233[k]
                    + f_3 * pc_y[k] * msk_1587[k];

        t_1987[k] = f_4 * msi0_1234[k]
                    - f_5 * msi1_1234[k]
                    + f_3 * pc_y[k] * msk_1588[k];

        t_1988[k] = f_3 * pc_y[k] * msk_1589[k];
    }

#pragma omp simd aligned(t_1989, t_1990, t_1991, pa_x, pc_x, pc_y, lsl0_1989, lsk_1593, \
                         lsl1_1989, msi0_1235, msi0_1236, msi1_1235, msi1_1236, msk_1590, \
                         msk_1591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1989[k] = pa_x[k] * lsl0_1989[k]
                    + f_19 * lsk_1593[k]
                    - f_14 * pc_x[k] * lsl1_1989[k];

        t_1990[k] = f_8 * msi0_1235[k]
                    - f_9 * msi1_1235[k]
                    + f_3 * pc_y[k] * msk_1590[k];

        t_1991[k] = f_6 * msi0_1236[k]
                    - f_7 * msi1_1236[k]
                    + f_3 * pc_y[k] * msk_1591[k];
    }

#pragma omp simd aligned(t_1992, t_1993, t_1994, pa_x, pc_x, pc_y, lsl0_1994, lsk_1598, \
                         lsl1_1994, msi0_1237, msi1_1237, msk_1592, \
                         msk_1593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1992[k] = f_4 * msi0_1237[k]
                    - f_5 * msi1_1237[k]
                    + f_3 * pc_y[k] * msk_1592[k];

        t_1993[k] = f_3 * pc_y[k] * msk_1593[k];

        t_1994[k] = pa_x[k] * lsl0_1994[k]
                    + f_18 * lsk_1598[k]
                    - f_14 * pc_x[k] * lsl1_1994[k];
    }

#pragma omp simd aligned(t_1995, t_1996, t_1997, pc_y, msi0_1238, msi0_1239, msi0_1240, \
                         msi1_1238, msi1_1239, msi1_1240, msk_1594, msk_1595, \
                         msk_1596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1995[k] = f_10 * msi0_1238[k]
                    - f_11 * msi1_1238[k]
                    + f_3 * pc_y[k] * msk_1594[k];

        t_1996[k] = f_8 * msi0_1239[k]
                    - f_9 * msi1_1239[k]
                    + f_3 * pc_y[k] * msk_1595[k];

        t_1997[k] = f_6 * msi0_1240[k]
                    - f_7 * msi1_1240[k]
                    + f_3 * pc_y[k] * msk_1596[k];
    }

#pragma omp simd aligned(t_1998, t_1999, t_2000, pa_x, pc_x, pc_y, lsl0_2000, lsk_1604, \
                         lsl1_2000, msi0_1241, msi1_1241, msk_1597, \
                         msk_1598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1998[k] = f_4 * msi0_1241[k]
                    - f_5 * msi1_1241[k]
                    + f_3 * pc_y[k] * msk_1597[k];

        t_1999[k] = f_3 * pc_y[k] * msk_1598[k];

        t_2000[k] = pa_x[k] * lsl0_2000[k]
                    + f_17 * lsk_1604[k]
                    - f_14 * pc_x[k] * lsl1_2000[k];
    }

#pragma omp simd aligned(t_2001, t_2002, t_2003, pc_y, msi0_1242, msi0_1243, msi0_1244, \
                         msi1_1242, msi1_1243, msi1_1244, msk_1599, msk_1600, \
                         msk_1601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2001[k] = f_12 * msi0_1242[k]
                    - f_13 * msi1_1242[k]
                    + f_3 * pc_y[k] * msk_1599[k];

        t_2002[k] = f_10 * msi0_1243[k]
                    - f_11 * msi1_1243[k]
                    + f_3 * pc_y[k] * msk_1600[k];

        t_2003[k] = f_8 * msi0_1244[k]
                    - f_9 * msi1_1244[k]
                    + f_3 * pc_y[k] * msk_1601[k];
    }

#pragma omp simd aligned(t_2004, t_2005, t_2006, pc_y, msi0_1245, msi0_1246, msi1_1245, \
                         msi1_1246, msk_1602, msk_1603, msk_1604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2004[k] = f_6 * msi0_1245[k]
                    - f_7 * msi1_1245[k]
                    + f_3 * pc_y[k] * msk_1602[k];

        t_2005[k] = f_4 * msi0_1246[k]
                    - f_5 * msi1_1246[k]
                    + f_3 * pc_y[k] * msk_1603[k];

        t_2006[k] = f_3 * pc_y[k] * msk_1604[k];
    }

#pragma omp simd aligned(t_2007, t_2008, t_2009, t_2010, pa_x, pc_x, lsl0_2007, lsk_1611, \
                         lsk_1612, lsk_1613, lsk_1614, lsl1_2007, msk_1612, msk_1613, \
                         msk_1614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2007[k] = pa_x[k] * lsl0_2007[k]
                    + f_16 * lsk_1611[k]
                    - f_14 * pc_x[k] * lsl1_2007[k];

        t_2008[k] = f_15 * lsk_1612[k]
                    + f_3 * pc_x[k] * msk_1612[k];

        t_2009[k] = f_15 * lsk_1613[k]
                    + f_3 * pc_x[k] * msk_1613[k];

        t_2010[k] = f_15 * lsk_1614[k]
                    + f_3 * pc_x[k] * msk_1614[k];
    }

#pragma omp simd aligned(t_2011, t_2012, t_2013, t_2014, t_2015, pc_x, pc_y, lsk_1615, \
                         lsk_1616, lsk_1617, lsk_1619, msk_1611, msk_1615, msk_1616, msk_1617, \
                         msk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2011[k] = f_15 * lsk_1615[k]
                    + f_3 * pc_x[k] * msk_1615[k];

        t_2012[k] = f_15 * lsk_1616[k]
                    + f_3 * pc_x[k] * msk_1616[k];

        t_2013[k] = f_15 * lsk_1617[k]
                    + f_3 * pc_x[k] * msk_1617[k];

        t_2014[k] = f_3 * pc_y[k] * msk_1611[k];

        t_2015[k] = f_15 * lsk_1619[k]
                    + f_3 * pc_x[k] * msk_1619[k];
    }

#pragma omp simd aligned(t_2016, t_2017, t_2018, t_2019, pa_x, pc_x, lsl0_2016, lsl0_2017, \
                         lsl0_2018, lsl0_2019, lsl1_2016, lsl1_2017, lsl1_2018, \
                         lsl1_2019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2016[k] = pa_x[k] * lsl0_2016[k]
                    - f_14 * pc_x[k] * lsl1_2016[k];

        t_2017[k] = pa_x[k] * lsl0_2017[k]
                    - f_14 * pc_x[k] * lsl1_2017[k];

        t_2018[k] = pa_x[k] * lsl0_2018[k]
                    - f_14 * pc_x[k] * lsl1_2018[k];

        t_2019[k] = pa_x[k] * lsl0_2019[k]
                    - f_14 * pc_x[k] * lsl1_2019[k];
    }

#pragma omp simd aligned(t_2020, t_2021, t_2022, t_2023, pa_x, pc_x, pc_y, lsl0_2020, \
                         lsl0_2021, lsl0_2022, lsl1_2020, lsl1_2021, lsl1_2022, \
                         msk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2020[k] = pa_x[k] * lsl0_2020[k]
                    - f_14 * pc_x[k] * lsl1_2020[k];

        t_2021[k] = pa_x[k] * lsl0_2021[k]
                    - f_14 * pc_x[k] * lsl1_2021[k];

        t_2022[k] = pa_x[k] * lsl0_2022[k]
                    - f_14 * pc_x[k] * lsl1_2022[k];

        t_2023[k] = f_3 * pc_y[k] * msk_1619[k];
    }

#pragma omp simd aligned(t_2024, t_2025, t_2026, t_2027, pa_x, pc_x, pc_z, lsl0_2024, \
                         lsl1_2024, msi0_1260, msi0_1261, msi1_1260, msi1_1261, msk_1620, \
                         msk_1621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2024[k] = pa_x[k] * lsl0_2024[k]
                    - f_14 * pc_x[k] * lsl1_2024[k];

        t_2025[k] = f_1 * msi0_1260[k]
                    - f_2 * msi1_1260[k]
                    + f_3 * pc_x[k] * msk_1620[k];

        t_2026[k] = f_22 * msi0_1261[k]
                    - f_23 * msi1_1261[k]
                    + f_3 * pc_x[k] * msk_1621[k];

        t_2027[k] = f_3 * pc_z[k] * msk_1620[k];
    }

#pragma omp simd aligned(t_2028, t_2029, t_2030, t_2031, pc_x, pc_z, msi0_1263, msi0_1265, \
                         msi0_1266, msi1_1263, msi1_1265, msi1_1266, msk_1621, msk_1623, \
                         msk_1625, msk_1626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2028[k] = f_12 * msi0_1263[k]
                    - f_13 * msi1_1263[k]
                    + f_3 * pc_x[k] * msk_1623[k];

        t_2029[k] = f_3 * pc_z[k] * msk_1621[k];

        t_2030[k] = f_12 * msi0_1265[k]
                    - f_13 * msi1_1265[k]
                    + f_3 * pc_x[k] * msk_1625[k];

        t_2031[k] = f_10 * msi0_1266[k]
                    - f_11 * msi1_1266[k]
                    + f_3 * pc_x[k] * msk_1626[k];
    }

#pragma omp simd aligned(t_2032, t_2033, t_2034, t_2035, pc_x, pc_z, msi0_1268, msi0_1269, \
                         msi0_1270, msi1_1268, msi1_1269, msi1_1270, msk_1623, msk_1628, \
                         msk_1629, msk_1630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2032[k] = f_3 * pc_z[k] * msk_1623[k];

        t_2033[k] = f_10 * msi0_1268[k]
                    - f_11 * msi1_1268[k]
                    + f_3 * pc_x[k] * msk_1628[k];

        t_2034[k] = f_10 * msi0_1269[k]
                    - f_11 * msi1_1269[k]
                    + f_3 * pc_x[k] * msk_1629[k];

        t_2035[k] = f_8 * msi0_1270[k]
                    - f_9 * msi1_1270[k]
                    + f_3 * pc_x[k] * msk_1630[k];
    }

#pragma omp simd aligned(t_2036, t_2037, t_2038, t_2039, pc_x, pc_z, msi0_1272, msi0_1273, \
                         msi0_1274, msi1_1272, msi1_1273, msi1_1274, msk_1626, msk_1632, \
                         msk_1633, msk_1634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2036[k] = f_3 * pc_z[k] * msk_1626[k];

        t_2037[k] = f_8 * msi0_1272[k]
                    - f_9 * msi1_1272[k]
                    + f_3 * pc_x[k] * msk_1632[k];

        t_2038[k] = f_8 * msi0_1273[k]
                    - f_9 * msi1_1273[k]
                    + f_3 * pc_x[k] * msk_1633[k];

        t_2039[k] = f_8 * msi0_1274[k]
                    - f_9 * msi1_1274[k]
                    + f_3 * pc_x[k] * msk_1634[k];
    }

#pragma omp simd aligned(t_2040, t_2041, t_2042, t_2043, pc_x, pc_z, msi0_1275, msi0_1277, \
                         msi0_1278, msi1_1275, msi1_1277, msi1_1278, msk_1630, msk_1635, \
                         msk_1637, msk_1638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2040[k] = f_6 * msi0_1275[k]
                    - f_7 * msi1_1275[k]
                    + f_3 * pc_x[k] * msk_1635[k];

        t_2041[k] = f_3 * pc_z[k] * msk_1630[k];

        t_2042[k] = f_6 * msi0_1277[k]
                    - f_7 * msi1_1277[k]
                    + f_3 * pc_x[k] * msk_1637[k];

        t_2043[k] = f_6 * msi0_1278[k]
                    - f_7 * msi1_1278[k]
                    + f_3 * pc_x[k] * msk_1638[k];
    }

#pragma omp simd aligned(t_2044, t_2045, t_2046, t_2047, pc_x, pc_z, msi0_1279, msi0_1280, \
                         msi0_1281, msi1_1279, msi1_1280, msi1_1281, msk_1635, msk_1639, \
                         msk_1640, msk_1641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2044[k] = f_6 * msi0_1279[k]
                    - f_7 * msi1_1279[k]
                    + f_3 * pc_x[k] * msk_1639[k];

        t_2045[k] = f_6 * msi0_1280[k]
                    - f_7 * msi1_1280[k]
                    + f_3 * pc_x[k] * msk_1640[k];

        t_2046[k] = f_4 * msi0_1281[k]
                    - f_5 * msi1_1281[k]
                    + f_3 * pc_x[k] * msk_1641[k];

        t_2047[k] = f_3 * pc_z[k] * msk_1635[k];
    }

#pragma omp simd aligned(t_2048, t_2049, t_2050, pc_x, msi0_1283, msi0_1284, msi0_1285, \
                         msi1_1283, msi1_1284, msi1_1285, msk_1643, msk_1644, \
                         msk_1645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2048[k] = f_4 * msi0_1283[k]
                    - f_5 * msi1_1283[k]
                    + f_3 * pc_x[k] * msk_1643[k];

        t_2049[k] = f_4 * msi0_1284[k]
                    - f_5 * msi1_1284[k]
                    + f_3 * pc_x[k] * msk_1644[k];

        t_2050[k] = f_4 * msi0_1285[k]
                    - f_5 * msi1_1285[k]
                    + f_3 * pc_x[k] * msk_1645[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);
    const auto f_24 = 3.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1620 = buffer.data(lsl0 + 1620);
    const auto *lsl0_1621 = buffer.data(lsl0 + 1621);
    const auto *lsl0_1623 = buffer.data(lsl0 + 1623);
    const auto *lsl0_1626 = buffer.data(lsl0 + 1626);
    const auto *lsl0_1630 = buffer.data(lsl0 + 1630);
    const auto *lsl0_1635 = buffer.data(lsl0 + 1635);
    const auto *lsl0_1641 = buffer.data(lsl0 + 1641);
    const auto *lsl0_1656 = buffer.data(lsl0 + 1656);
    const auto *lsl0_1658 = buffer.data(lsl0 + 1658);
    const auto *lsl0_1659 = buffer.data(lsl0 + 1659);
    const auto *lsl0_1660 = buffer.data(lsl0 + 1660);
    const auto *lsl0_1661 = buffer.data(lsl0 + 1661);
    const auto *lsl0_1662 = buffer.data(lsl0 + 1662);

    const auto *lsk_1324 = buffer.data(lsk + 1324);
    const auto *lsk_1325 = buffer.data(lsk + 1325);
    const auto *lsk_1326 = buffer.data(lsk + 1326);
    const auto *lsk_1327 = buffer.data(lsk + 1327);
    const auto *lsk_1328 = buffer.data(lsk + 1328);
    const auto *lsk_1329 = buffer.data(lsk + 1329);
    const auto *lsk_1331 = buffer.data(lsk + 1331);
    const auto *lsk_1360 = buffer.data(lsk + 1360);
    const auto *lsk_1367 = buffer.data(lsk + 1367);
    const auto *lsk_1396 = buffer.data(lsk + 1396);
    const auto *lsk_1398 = buffer.data(lsk + 1398);
    const auto *lsk_1399 = buffer.data(lsk + 1399);
    const auto *lsk_1400 = buffer.data(lsk + 1400);
    const auto *lsk_1401 = buffer.data(lsk + 1401);
    const auto *lsk_1402 = buffer.data(lsk + 1402);
    const auto *lsk_1403 = buffer.data(lsk + 1403);

    const auto *lsl1_1620 = buffer.data(lsl1 + 1620);
    const auto *lsl1_1621 = buffer.data(lsl1 + 1621);
    const auto *lsl1_1623 = buffer.data(lsl1 + 1623);
    const auto *lsl1_1626 = buffer.data(lsl1 + 1626);
    const auto *lsl1_1630 = buffer.data(lsl1 + 1630);
    const auto *lsl1_1635 = buffer.data(lsl1 + 1635);
    const auto *lsl1_1641 = buffer.data(lsl1 + 1641);
    const auto *lsl1_1656 = buffer.data(lsl1 + 1656);
    const auto *lsl1_1658 = buffer.data(lsl1 + 1658);
    const auto *lsl1_1659 = buffer.data(lsl1 + 1659);
    const auto *lsl1_1660 = buffer.data(lsl1 + 1660);
    const auto *lsl1_1661 = buffer.data(lsl1 + 1661);
    const auto *lsl1_1662 = buffer.data(lsl1 + 1662);

    const auto *msi0_1281 = buffer.data(msi0 + 1281);
    const auto *msi0_1282 = buffer.data(msi0 + 1282);
    const auto *msi0_1283 = buffer.data(msi0 + 1283);
    const auto *msi0_1284 = buffer.data(msi0 + 1284);
    const auto *msi0_1285 = buffer.data(msi0 + 1285);
    const auto *msi0_1286 = buffer.data(msi0 + 1286);
    const auto *msi0_1287 = buffer.data(msi0 + 1287);
    const auto *msi0_1290 = buffer.data(msi0 + 1290);
    const auto *msi0_1292 = buffer.data(msi0 + 1292);
    const auto *msi0_1293 = buffer.data(msi0 + 1293);
    const auto *msi0_1295 = buffer.data(msi0 + 1295);
    const auto *msi0_1296 = buffer.data(msi0 + 1296);
    const auto *msi0_1297 = buffer.data(msi0 + 1297);
    const auto *msi0_1299 = buffer.data(msi0 + 1299);
    const auto *msi0_1300 = buffer.data(msi0 + 1300);
    const auto *msi0_1301 = buffer.data(msi0 + 1301);
    const auto *msi0_1302 = buffer.data(msi0 + 1302);
    const auto *msi0_1304 = buffer.data(msi0 + 1304);
    const auto *msi0_1305 = buffer.data(msi0 + 1305);
    const auto *msi0_1306 = buffer.data(msi0 + 1306);
    const auto *msi0_1307 = buffer.data(msi0 + 1307);
    const auto *msi0_1308 = buffer.data(msi0 + 1308);
    const auto *msi0_1310 = buffer.data(msi0 + 1310);
    const auto *msi0_1311 = buffer.data(msi0 + 1311);
    const auto *msi0_1312 = buffer.data(msi0 + 1312);
    const auto *msi0_1313 = buffer.data(msi0 + 1313);
    const auto *msi0_1314 = buffer.data(msi0 + 1314);
    const auto *msi0_1315 = buffer.data(msi0 + 1315);
    const auto *msi0_1316 = buffer.data(msi0 + 1316);
    const auto *msi0_1317 = buffer.data(msi0 + 1317);
    const auto *msi0_1318 = buffer.data(msi0 + 1318);
    const auto *msi0_1319 = buffer.data(msi0 + 1319);
    const auto *msi0_1320 = buffer.data(msi0 + 1320);
    const auto *msi0_1321 = buffer.data(msi0 + 1321);
    const auto *msi0_1322 = buffer.data(msi0 + 1322);
    const auto *msi0_1323 = buffer.data(msi0 + 1323);
    const auto *msi0_1324 = buffer.data(msi0 + 1324);
    const auto *msi0_1325 = buffer.data(msi0 + 1325);
    const auto *msi0_1326 = buffer.data(msi0 + 1326);
    const auto *msi0_1327 = buffer.data(msi0 + 1327);
    const auto *msi0_1328 = buffer.data(msi0 + 1328);
    const auto *msi0_1329 = buffer.data(msi0 + 1329);
    const auto *msi0_1330 = buffer.data(msi0 + 1330);
    const auto *msi0_1331 = buffer.data(msi0 + 1331);
    const auto *msi0_1332 = buffer.data(msi0 + 1332);
    const auto *msi0_1333 = buffer.data(msi0 + 1333);
    const auto *msi0_1334 = buffer.data(msi0 + 1334);
    const auto *msi0_1335 = buffer.data(msi0 + 1335);
    const auto *msi0_1336 = buffer.data(msi0 + 1336);
    const auto *msi0_1337 = buffer.data(msi0 + 1337);
    const auto *msi0_1338 = buffer.data(msi0 + 1338);
    const auto *msi0_1339 = buffer.data(msi0 + 1339);
    const auto *msi0_1340 = buffer.data(msi0 + 1340);
    const auto *msi0_1341 = buffer.data(msi0 + 1341);
    const auto *msi0_1342 = buffer.data(msi0 + 1342);
    const auto *msi0_1343 = buffer.data(msi0 + 1343);
    const auto *msi0_1344 = buffer.data(msi0 + 1344);
    const auto *msi0_1345 = buffer.data(msi0 + 1345);
    const auto *msi0_1346 = buffer.data(msi0 + 1346);
    const auto *msi0_1347 = buffer.data(msi0 + 1347);
    const auto *msi0_1348 = buffer.data(msi0 + 1348);
    const auto *msi0_1349 = buffer.data(msi0 + 1349);
    const auto *msi0_1350 = buffer.data(msi0 + 1350);
    const auto *msi0_1351 = buffer.data(msi0 + 1351);
    const auto *msi0_1352 = buffer.data(msi0 + 1352);

    const auto *msi1_1281 = buffer.data(msi1 + 1281);
    const auto *msi1_1282 = buffer.data(msi1 + 1282);
    const auto *msi1_1283 = buffer.data(msi1 + 1283);
    const auto *msi1_1284 = buffer.data(msi1 + 1284);
    const auto *msi1_1285 = buffer.data(msi1 + 1285);
    const auto *msi1_1286 = buffer.data(msi1 + 1286);
    const auto *msi1_1287 = buffer.data(msi1 + 1287);
    const auto *msi1_1290 = buffer.data(msi1 + 1290);
    const auto *msi1_1292 = buffer.data(msi1 + 1292);
    const auto *msi1_1293 = buffer.data(msi1 + 1293);
    const auto *msi1_1295 = buffer.data(msi1 + 1295);
    const auto *msi1_1296 = buffer.data(msi1 + 1296);
    const auto *msi1_1297 = buffer.data(msi1 + 1297);
    const auto *msi1_1299 = buffer.data(msi1 + 1299);
    const auto *msi1_1300 = buffer.data(msi1 + 1300);
    const auto *msi1_1301 = buffer.data(msi1 + 1301);
    const auto *msi1_1302 = buffer.data(msi1 + 1302);
    const auto *msi1_1304 = buffer.data(msi1 + 1304);
    const auto *msi1_1305 = buffer.data(msi1 + 1305);
    const auto *msi1_1306 = buffer.data(msi1 + 1306);
    const auto *msi1_1307 = buffer.data(msi1 + 1307);
    const auto *msi1_1308 = buffer.data(msi1 + 1308);
    const auto *msi1_1310 = buffer.data(msi1 + 1310);
    const auto *msi1_1311 = buffer.data(msi1 + 1311);
    const auto *msi1_1312 = buffer.data(msi1 + 1312);
    const auto *msi1_1313 = buffer.data(msi1 + 1313);
    const auto *msi1_1314 = buffer.data(msi1 + 1314);
    const auto *msi1_1315 = buffer.data(msi1 + 1315);
    const auto *msi1_1316 = buffer.data(msi1 + 1316);
    const auto *msi1_1317 = buffer.data(msi1 + 1317);
    const auto *msi1_1318 = buffer.data(msi1 + 1318);
    const auto *msi1_1319 = buffer.data(msi1 + 1319);
    const auto *msi1_1320 = buffer.data(msi1 + 1320);
    const auto *msi1_1321 = buffer.data(msi1 + 1321);
    const auto *msi1_1322 = buffer.data(msi1 + 1322);
    const auto *msi1_1323 = buffer.data(msi1 + 1323);
    const auto *msi1_1324 = buffer.data(msi1 + 1324);
    const auto *msi1_1325 = buffer.data(msi1 + 1325);
    const auto *msi1_1326 = buffer.data(msi1 + 1326);
    const auto *msi1_1327 = buffer.data(msi1 + 1327);
    const auto *msi1_1328 = buffer.data(msi1 + 1328);
    const auto *msi1_1329 = buffer.data(msi1 + 1329);
    const auto *msi1_1330 = buffer.data(msi1 + 1330);
    const auto *msi1_1331 = buffer.data(msi1 + 1331);
    const auto *msi1_1332 = buffer.data(msi1 + 1332);
    const auto *msi1_1333 = buffer.data(msi1 + 1333);
    const auto *msi1_1334 = buffer.data(msi1 + 1334);
    const auto *msi1_1335 = buffer.data(msi1 + 1335);
    const auto *msi1_1336 = buffer.data(msi1 + 1336);
    const auto *msi1_1337 = buffer.data(msi1 + 1337);
    const auto *msi1_1338 = buffer.data(msi1 + 1338);
    const auto *msi1_1339 = buffer.data(msi1 + 1339);
    const auto *msi1_1340 = buffer.data(msi1 + 1340);
    const auto *msi1_1341 = buffer.data(msi1 + 1341);
    const auto *msi1_1342 = buffer.data(msi1 + 1342);
    const auto *msi1_1343 = buffer.data(msi1 + 1343);
    const auto *msi1_1344 = buffer.data(msi1 + 1344);
    const auto *msi1_1345 = buffer.data(msi1 + 1345);
    const auto *msi1_1346 = buffer.data(msi1 + 1346);
    const auto *msi1_1347 = buffer.data(msi1 + 1347);
    const auto *msi1_1348 = buffer.data(msi1 + 1348);
    const auto *msi1_1349 = buffer.data(msi1 + 1349);
    const auto *msi1_1350 = buffer.data(msi1 + 1350);
    const auto *msi1_1351 = buffer.data(msi1 + 1351);
    const auto *msi1_1352 = buffer.data(msi1 + 1352);

    const auto *msk_1646 = buffer.data(msk + 1646);
    const auto *msk_1647 = buffer.data(msk + 1647);
    const auto *msk_1648 = buffer.data(msk + 1648);
    const auto *msk_1649 = buffer.data(msk + 1649);
    const auto *msk_1650 = buffer.data(msk + 1650);
    const auto *msk_1651 = buffer.data(msk + 1651);
    const auto *msk_1652 = buffer.data(msk + 1652);
    const auto *msk_1653 = buffer.data(msk + 1653);
    const auto *msk_1654 = buffer.data(msk + 1654);
    const auto *msk_1655 = buffer.data(msk + 1655);
    const auto *msk_1658 = buffer.data(msk + 1658);
    const auto *msk_1660 = buffer.data(msk + 1660);
    const auto *msk_1661 = buffer.data(msk + 1661);
    const auto *msk_1663 = buffer.data(msk + 1663);
    const auto *msk_1664 = buffer.data(msk + 1664);
    const auto *msk_1665 = buffer.data(msk + 1665);
    const auto *msk_1667 = buffer.data(msk + 1667);
    const auto *msk_1668 = buffer.data(msk + 1668);
    const auto *msk_1669 = buffer.data(msk + 1669);
    const auto *msk_1670 = buffer.data(msk + 1670);
    const auto *msk_1672 = buffer.data(msk + 1672);
    const auto *msk_1673 = buffer.data(msk + 1673);
    const auto *msk_1674 = buffer.data(msk + 1674);
    const auto *msk_1675 = buffer.data(msk + 1675);
    const auto *msk_1676 = buffer.data(msk + 1676);
    const auto *msk_1678 = buffer.data(msk + 1678);
    const auto *msk_1679 = buffer.data(msk + 1679);
    const auto *msk_1680 = buffer.data(msk + 1680);
    const auto *msk_1681 = buffer.data(msk + 1681);
    const auto *msk_1682 = buffer.data(msk + 1682);
    const auto *msk_1683 = buffer.data(msk + 1683);
    const auto *msk_1684 = buffer.data(msk + 1684);
    const auto *msk_1685 = buffer.data(msk + 1685);
    const auto *msk_1686 = buffer.data(msk + 1686);
    const auto *msk_1687 = buffer.data(msk + 1687);
    const auto *msk_1688 = buffer.data(msk + 1688);
    const auto *msk_1689 = buffer.data(msk + 1689);
    const auto *msk_1690 = buffer.data(msk + 1690);
    const auto *msk_1691 = buffer.data(msk + 1691);
    const auto *msk_1692 = buffer.data(msk + 1692);
    const auto *msk_1693 = buffer.data(msk + 1693);
    const auto *msk_1694 = buffer.data(msk + 1694);
    const auto *msk_1695 = buffer.data(msk + 1695);
    const auto *msk_1696 = buffer.data(msk + 1696);
    const auto *msk_1697 = buffer.data(msk + 1697);
    const auto *msk_1698 = buffer.data(msk + 1698);
    const auto *msk_1699 = buffer.data(msk + 1699);
    const auto *msk_1700 = buffer.data(msk + 1700);
    const auto *msk_1701 = buffer.data(msk + 1701);
    const auto *msk_1702 = buffer.data(msk + 1702);
    const auto *msk_1703 = buffer.data(msk + 1703);
    const auto *msk_1704 = buffer.data(msk + 1704);
    const auto *msk_1705 = buffer.data(msk + 1705);
    const auto *msk_1706 = buffer.data(msk + 1706);
    const auto *msk_1707 = buffer.data(msk + 1707);
    const auto *msk_1708 = buffer.data(msk + 1708);
    const auto *msk_1709 = buffer.data(msk + 1709);
    const auto *msk_1710 = buffer.data(msk + 1710);
    const auto *msk_1711 = buffer.data(msk + 1711);
    const auto *msk_1712 = buffer.data(msk + 1712);
    const auto *msk_1713 = buffer.data(msk + 1713);
    const auto *msk_1714 = buffer.data(msk + 1714);
    const auto *msk_1715 = buffer.data(msk + 1715);
    const auto *msk_1716 = buffer.data(msk + 1716);
    const auto *msk_1717 = buffer.data(msk + 1717);
    const auto *msk_1718 = buffer.data(msk + 1718);
    const auto *msk_1719 = buffer.data(msk + 1719);
    const auto *msk_1720 = buffer.data(msk + 1720);
    const auto *msk_1721 = buffer.data(msk + 1721);
    const auto *msk_1722 = buffer.data(msk + 1722);
    const auto *msk_1723 = buffer.data(msk + 1723);
    const auto *msk_1724 = buffer.data(msk + 1724);
    const auto *msk_1725 = buffer.data(msk + 1725);
    const auto *msk_1726 = buffer.data(msk + 1726);
    const auto *msk_1727 = buffer.data(msk + 1727);
    const auto *msk_1728 = buffer.data(msk + 1728);
    const auto *msk_1729 = buffer.data(msk + 1729);
    const auto *msk_1730 = buffer.data(msk + 1730);
    const auto *msk_1731 = buffer.data(msk + 1731);
    const auto *msk_1732 = buffer.data(msk + 1732);
    const auto *msk_1733 = buffer.data(msk + 1733);
    const auto *msk_1734 = buffer.data(msk + 1734);
    const auto *msk_1735 = buffer.data(msk + 1735);
    const auto *msk_1736 = buffer.data(msk + 1736);

#pragma omp simd aligned(t_2051, t_2052, t_2053, t_2054, t_2055, pc_x, msi0_1286, msi0_1287, \
                         msi1_1286, msi1_1287, msk_1646, msk_1647, msk_1648, msk_1649, \
                         msk_1650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2051[k] = f_4 * msi0_1286[k]
                    - f_5 * msi1_1286[k]
                    + f_3 * pc_x[k] * msk_1646[k];

        t_2052[k] = f_4 * msi0_1287[k]
                    - f_5 * msi1_1287[k]
                    + f_3 * pc_x[k] * msk_1647[k];

        t_2053[k] = f_3 * pc_x[k] * msk_1648[k];

        t_2054[k] = f_3 * pc_x[k] * msk_1649[k];

        t_2055[k] = f_3 * pc_x[k] * msk_1650[k];
    }

#pragma omp simd aligned(t_2056, t_2057, t_2058, t_2059, t_2060, pc_x, msk_1651, msk_1652, \
                         msk_1653, msk_1654, msk_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2056[k] = f_3 * pc_x[k] * msk_1651[k];

        t_2057[k] = f_3 * pc_x[k] * msk_1652[k];

        t_2058[k] = f_3 * pc_x[k] * msk_1653[k];

        t_2059[k] = f_3 * pc_x[k] * msk_1654[k];

        t_2060[k] = f_3 * pc_x[k] * msk_1655[k];
    }

#pragma omp simd aligned(t_2061, t_2062, t_2063, t_2064, pc_y, pc_z, lsk_1324, msi0_1281, \
                         msi0_1282, msi1_1281, msi1_1282, msk_1648, msk_1649, \
                         msk_1650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2061[k] = f_0 * lsk_1324[k]
                    + f_1 * msi0_1281[k]
                    - f_2 * msi1_1281[k]
                    + f_3 * pc_y[k] * msk_1648[k];

        t_2062[k] = f_3 * pc_z[k] * msk_1648[k];

        t_2063[k] = f_4 * msi0_1281[k]
                    - f_5 * msi1_1281[k]
                    + f_3 * pc_z[k] * msk_1649[k];

        t_2064[k] = f_6 * msi0_1282[k]
                    - f_7 * msi1_1282[k]
                    + f_3 * pc_z[k] * msk_1650[k];
    }

#pragma omp simd aligned(t_2065, t_2066, t_2067, pc_z, msi0_1283, msi0_1284, msi0_1285, \
                         msi1_1283, msi1_1284, msi1_1285, msk_1651, msk_1652, \
                         msk_1653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2065[k] = f_8 * msi0_1283[k]
                    - f_9 * msi1_1283[k]
                    + f_3 * pc_z[k] * msk_1651[k];

        t_2066[k] = f_10 * msi0_1284[k]
                    - f_11 * msi1_1284[k]
                    + f_3 * pc_z[k] * msk_1652[k];

        t_2067[k] = f_12 * msi0_1285[k]
                    - f_13 * msi1_1285[k]
                    + f_3 * pc_z[k] * msk_1653[k];
    }

#pragma omp simd aligned(t_2068, t_2069, t_2070, t_2071, pa_z, pc_y, pc_z, lsl0_1620, \
                         lsl0_1621, lsk_1331, lsl1_1620, lsl1_1621, msi0_1287, msi1_1287, \
                         msk_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2068[k] = f_0 * lsk_1331[k]
                    + f_3 * pc_y[k] * msk_1655[k];

        t_2069[k] = f_1 * msi0_1287[k]
                    - f_2 * msi1_1287[k]
                    + f_3 * pc_z[k] * msk_1655[k];

        t_2070[k] = pa_z[k] * lsl0_1620[k]
                    - f_14 * pc_z[k] * lsl1_1620[k];

        t_2071[k] = pa_z[k] * lsl0_1621[k]
                    - f_14 * pc_z[k] * lsl1_1621[k];
    }

#pragma omp simd aligned(t_2072, t_2073, t_2074, pa_z, pc_x, pc_z, lsl0_1623, lsl1_1623, \
                         msi0_1290, msi0_1292, msi1_1290, msi1_1292, msk_1658, \
                         msk_1660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2072[k] = f_22 * msi0_1290[k]
                    - f_23 * msi1_1290[k]
                    + f_3 * pc_x[k] * msk_1658[k];

        t_2073[k] = pa_z[k] * lsl0_1623[k]
                    - f_14 * pc_z[k] * lsl1_1623[k];

        t_2074[k] = f_12 * msi0_1292[k]
                    - f_13 * msi1_1292[k]
                    + f_3 * pc_x[k] * msk_1660[k];
    }

#pragma omp simd aligned(t_2075, t_2076, t_2077, pa_z, pc_x, pc_z, lsl0_1626, lsl1_1626, \
                         msi0_1293, msi0_1295, msi1_1293, msi1_1295, msk_1661, \
                         msk_1663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2075[k] = f_12 * msi0_1293[k]
                    - f_13 * msi1_1293[k]
                    + f_3 * pc_x[k] * msk_1661[k];

        t_2076[k] = pa_z[k] * lsl0_1626[k]
                    - f_14 * pc_z[k] * lsl1_1626[k];

        t_2077[k] = f_10 * msi0_1295[k]
                    - f_11 * msi1_1295[k]
                    + f_3 * pc_x[k] * msk_1663[k];
    }

#pragma omp simd aligned(t_2078, t_2079, t_2080, pa_z, pc_x, pc_z, lsl0_1630, lsl1_1630, \
                         msi0_1296, msi0_1297, msi1_1296, msi1_1297, msk_1664, \
                         msk_1665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2078[k] = f_10 * msi0_1296[k]
                    - f_11 * msi1_1296[k]
                    + f_3 * pc_x[k] * msk_1664[k];

        t_2079[k] = f_10 * msi0_1297[k]
                    - f_11 * msi1_1297[k]
                    + f_3 * pc_x[k] * msk_1665[k];

        t_2080[k] = pa_z[k] * lsl0_1630[k]
                    - f_14 * pc_z[k] * lsl1_1630[k];
    }

#pragma omp simd aligned(t_2081, t_2082, t_2083, pc_x, msi0_1299, msi0_1300, msi0_1301, \
                         msi1_1299, msi1_1300, msi1_1301, msk_1667, msk_1668, \
                         msk_1669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2081[k] = f_8 * msi0_1299[k]
                    - f_9 * msi1_1299[k]
                    + f_3 * pc_x[k] * msk_1667[k];

        t_2082[k] = f_8 * msi0_1300[k]
                    - f_9 * msi1_1300[k]
                    + f_3 * pc_x[k] * msk_1668[k];

        t_2083[k] = f_8 * msi0_1301[k]
                    - f_9 * msi1_1301[k]
                    + f_3 * pc_x[k] * msk_1669[k];
    }

#pragma omp simd aligned(t_2084, t_2085, t_2086, pa_z, pc_x, pc_z, lsl0_1635, lsl1_1635, \
                         msi0_1302, msi0_1304, msi1_1302, msi1_1304, msk_1670, \
                         msk_1672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2084[k] = f_8 * msi0_1302[k]
                    - f_9 * msi1_1302[k]
                    + f_3 * pc_x[k] * msk_1670[k];

        t_2085[k] = pa_z[k] * lsl0_1635[k]
                    - f_14 * pc_z[k] * lsl1_1635[k];

        t_2086[k] = f_6 * msi0_1304[k]
                    - f_7 * msi1_1304[k]
                    + f_3 * pc_x[k] * msk_1672[k];
    }

#pragma omp simd aligned(t_2087, t_2088, t_2089, pc_x, msi0_1305, msi0_1306, msi0_1307, \
                         msi1_1305, msi1_1306, msi1_1307, msk_1673, msk_1674, \
                         msk_1675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2087[k] = f_6 * msi0_1305[k]
                    - f_7 * msi1_1305[k]
                    + f_3 * pc_x[k] * msk_1673[k];

        t_2088[k] = f_6 * msi0_1306[k]
                    - f_7 * msi1_1306[k]
                    + f_3 * pc_x[k] * msk_1674[k];

        t_2089[k] = f_6 * msi0_1307[k]
                    - f_7 * msi1_1307[k]
                    + f_3 * pc_x[k] * msk_1675[k];
    }

#pragma omp simd aligned(t_2090, t_2091, t_2092, pa_z, pc_x, pc_z, lsl0_1641, lsl1_1641, \
                         msi0_1308, msi0_1310, msi1_1308, msi1_1310, msk_1676, \
                         msk_1678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2090[k] = f_6 * msi0_1308[k]
                    - f_7 * msi1_1308[k]
                    + f_3 * pc_x[k] * msk_1676[k];

        t_2091[k] = pa_z[k] * lsl0_1641[k]
                    - f_14 * pc_z[k] * lsl1_1641[k];

        t_2092[k] = f_4 * msi0_1310[k]
                    - f_5 * msi1_1310[k]
                    + f_3 * pc_x[k] * msk_1678[k];
    }

#pragma omp simd aligned(t_2093, t_2094, t_2095, pc_x, msi0_1311, msi0_1312, msi0_1313, \
                         msi1_1311, msi1_1312, msi1_1313, msk_1679, msk_1680, \
                         msk_1681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2093[k] = f_4 * msi0_1311[k]
                    - f_5 * msi1_1311[k]
                    + f_3 * pc_x[k] * msk_1679[k];

        t_2094[k] = f_4 * msi0_1312[k]
                    - f_5 * msi1_1312[k]
                    + f_3 * pc_x[k] * msk_1680[k];

        t_2095[k] = f_4 * msi0_1313[k]
                    - f_5 * msi1_1313[k]
                    + f_3 * pc_x[k] * msk_1681[k];
    }

#pragma omp simd aligned(t_2096, t_2097, t_2098, t_2099, t_2100, pc_x, msi0_1314, msi0_1315, \
                         msi1_1314, msi1_1315, msk_1682, msk_1683, msk_1684, msk_1685, \
                         msk_1686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2096[k] = f_4 * msi0_1314[k]
                    - f_5 * msi1_1314[k]
                    + f_3 * pc_x[k] * msk_1682[k];

        t_2097[k] = f_4 * msi0_1315[k]
                    - f_5 * msi1_1315[k]
                    + f_3 * pc_x[k] * msk_1683[k];

        t_2098[k] = f_3 * pc_x[k] * msk_1684[k];

        t_2099[k] = f_3 * pc_x[k] * msk_1685[k];

        t_2100[k] = f_3 * pc_x[k] * msk_1686[k];
    }

#pragma omp simd aligned(t_2101, t_2102, t_2103, t_2104, t_2105, t_2106, pa_z, pc_x, pc_z, \
                         lsl0_1656, lsl1_1656, msk_1687, msk_1688, msk_1689, msk_1690, \
                         msk_1691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2101[k] = f_3 * pc_x[k] * msk_1687[k];

        t_2102[k] = f_3 * pc_x[k] * msk_1688[k];

        t_2103[k] = f_3 * pc_x[k] * msk_1689[k];

        t_2104[k] = f_3 * pc_x[k] * msk_1690[k];

        t_2105[k] = f_3 * pc_x[k] * msk_1691[k];

        t_2106[k] = pa_z[k] * lsl0_1656[k]
                    - f_14 * pc_z[k] * lsl1_1656[k];
    }

#pragma omp simd aligned(t_2107, t_2108, t_2109, pa_z, pc_z, lsl0_1658, lsl0_1659, lsk_1324, \
                         lsk_1325, lsk_1326, lsl1_1658, lsl1_1659, \
                         msk_1684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2107[k] = f_15 * lsk_1324[k]
                    + f_3 * pc_z[k] * msk_1684[k];

        t_2108[k] = pa_z[k] * lsl0_1658[k]
                    + f_16 * lsk_1325[k]
                    - f_14 * pc_z[k] * lsl1_1658[k];

        t_2109[k] = pa_z[k] * lsl0_1659[k]
                    + f_17 * lsk_1326[k]
                    - f_14 * pc_z[k] * lsl1_1659[k];
    }

#pragma omp simd aligned(t_2110, t_2111, t_2112, pa_z, pc_z, lsl0_1660, lsl0_1661, lsl0_1662, \
                         lsk_1327, lsk_1328, lsk_1329, lsl1_1660, lsl1_1661, \
                         lsl1_1662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2110[k] = pa_z[k] * lsl0_1660[k]
                    + f_18 * lsk_1327[k]
                    - f_14 * pc_z[k] * lsl1_1660[k];

        t_2111[k] = pa_z[k] * lsl0_1661[k]
                    + f_19 * lsk_1328[k]
                    - f_14 * pc_z[k] * lsl1_1661[k];

        t_2112[k] = pa_z[k] * lsl0_1662[k]
                    + f_20 * lsk_1329[k]
                    - f_14 * pc_z[k] * lsl1_1662[k];
    }

#pragma omp simd aligned(t_2113, t_2114, t_2115, pc_x, pc_y, pc_z, lsk_1331, lsk_1367, \
                         msi0_1315, msi0_1316, msi1_1315, msi1_1316, msk_1691, \
                         msk_1692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2113[k] = f_21 * lsk_1367[k]
                    + f_3 * pc_y[k] * msk_1691[k];

        t_2114[k] = f_15 * lsk_1331[k]
                    + f_1 * msi0_1315[k]
                    - f_2 * msi1_1315[k]
                    + f_3 * pc_z[k] * msk_1691[k];

        t_2115[k] = f_1 * msi0_1316[k]
                    - f_2 * msi1_1316[k]
                    + f_3 * pc_x[k] * msk_1692[k];
    }

#pragma omp simd aligned(t_2116, t_2117, t_2118, pc_x, msi0_1317, msi0_1318, msi0_1319, \
                         msi1_1317, msi1_1318, msi1_1319, msk_1693, msk_1694, \
                         msk_1695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2116[k] = f_22 * msi0_1317[k]
                    - f_23 * msi1_1317[k]
                    + f_3 * pc_x[k] * msk_1693[k];

        t_2117[k] = f_22 * msi0_1318[k]
                    - f_23 * msi1_1318[k]
                    + f_3 * pc_x[k] * msk_1694[k];

        t_2118[k] = f_12 * msi0_1319[k]
                    - f_13 * msi1_1319[k]
                    + f_3 * pc_x[k] * msk_1695[k];
    }

#pragma omp simd aligned(t_2119, t_2120, t_2121, pc_x, msi0_1320, msi0_1321, msi0_1322, \
                         msi1_1320, msi1_1321, msi1_1322, msk_1696, msk_1697, \
                         msk_1698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2119[k] = f_12 * msi0_1320[k]
                    - f_13 * msi1_1320[k]
                    + f_3 * pc_x[k] * msk_1696[k];

        t_2120[k] = f_12 * msi0_1321[k]
                    - f_13 * msi1_1321[k]
                    + f_3 * pc_x[k] * msk_1697[k];

        t_2121[k] = f_10 * msi0_1322[k]
                    - f_11 * msi1_1322[k]
                    + f_3 * pc_x[k] * msk_1698[k];
    }

#pragma omp simd aligned(t_2122, t_2123, t_2124, pc_x, msi0_1323, msi0_1324, msi0_1325, \
                         msi1_1323, msi1_1324, msi1_1325, msk_1699, msk_1700, \
                         msk_1701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2122[k] = f_10 * msi0_1323[k]
                    - f_11 * msi1_1323[k]
                    + f_3 * pc_x[k] * msk_1699[k];

        t_2123[k] = f_10 * msi0_1324[k]
                    - f_11 * msi1_1324[k]
                    + f_3 * pc_x[k] * msk_1700[k];

        t_2124[k] = f_10 * msi0_1325[k]
                    - f_11 * msi1_1325[k]
                    + f_3 * pc_x[k] * msk_1701[k];
    }

#pragma omp simd aligned(t_2125, t_2126, t_2127, pc_x, msi0_1326, msi0_1327, msi0_1328, \
                         msi1_1326, msi1_1327, msi1_1328, msk_1702, msk_1703, \
                         msk_1704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2125[k] = f_8 * msi0_1326[k]
                    - f_9 * msi1_1326[k]
                    + f_3 * pc_x[k] * msk_1702[k];

        t_2126[k] = f_8 * msi0_1327[k]
                    - f_9 * msi1_1327[k]
                    + f_3 * pc_x[k] * msk_1703[k];

        t_2127[k] = f_8 * msi0_1328[k]
                    - f_9 * msi1_1328[k]
                    + f_3 * pc_x[k] * msk_1704[k];
    }

#pragma omp simd aligned(t_2128, t_2129, t_2130, pc_x, msi0_1329, msi0_1330, msi0_1331, \
                         msi1_1329, msi1_1330, msi1_1331, msk_1705, msk_1706, \
                         msk_1707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2128[k] = f_8 * msi0_1329[k]
                    - f_9 * msi1_1329[k]
                    + f_3 * pc_x[k] * msk_1705[k];

        t_2129[k] = f_8 * msi0_1330[k]
                    - f_9 * msi1_1330[k]
                    + f_3 * pc_x[k] * msk_1706[k];

        t_2130[k] = f_6 * msi0_1331[k]
                    - f_7 * msi1_1331[k]
                    + f_3 * pc_x[k] * msk_1707[k];
    }

#pragma omp simd aligned(t_2131, t_2132, t_2133, pc_x, msi0_1332, msi0_1333, msi0_1334, \
                         msi1_1332, msi1_1333, msi1_1334, msk_1708, msk_1709, \
                         msk_1710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2131[k] = f_6 * msi0_1332[k]
                    - f_7 * msi1_1332[k]
                    + f_3 * pc_x[k] * msk_1708[k];

        t_2132[k] = f_6 * msi0_1333[k]
                    - f_7 * msi1_1333[k]
                    + f_3 * pc_x[k] * msk_1709[k];

        t_2133[k] = f_6 * msi0_1334[k]
                    - f_7 * msi1_1334[k]
                    + f_3 * pc_x[k] * msk_1710[k];
    }

#pragma omp simd aligned(t_2134, t_2135, t_2136, pc_x, msi0_1335, msi0_1336, msi0_1337, \
                         msi1_1335, msi1_1336, msi1_1337, msk_1711, msk_1712, \
                         msk_1713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2134[k] = f_6 * msi0_1335[k]
                    - f_7 * msi1_1335[k]
                    + f_3 * pc_x[k] * msk_1711[k];

        t_2135[k] = f_6 * msi0_1336[k]
                    - f_7 * msi1_1336[k]
                    + f_3 * pc_x[k] * msk_1712[k];

        t_2136[k] = f_4 * msi0_1337[k]
                    - f_5 * msi1_1337[k]
                    + f_3 * pc_x[k] * msk_1713[k];
    }

#pragma omp simd aligned(t_2137, t_2138, t_2139, pc_x, msi0_1338, msi0_1339, msi0_1340, \
                         msi1_1338, msi1_1339, msi1_1340, msk_1714, msk_1715, \
                         msk_1716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2137[k] = f_4 * msi0_1338[k]
                    - f_5 * msi1_1338[k]
                    + f_3 * pc_x[k] * msk_1714[k];

        t_2138[k] = f_4 * msi0_1339[k]
                    - f_5 * msi1_1339[k]
                    + f_3 * pc_x[k] * msk_1715[k];

        t_2139[k] = f_4 * msi0_1340[k]
                    - f_5 * msi1_1340[k]
                    + f_3 * pc_x[k] * msk_1716[k];
    }

#pragma omp simd aligned(t_2140, t_2141, t_2142, t_2143, pc_x, msi0_1341, msi0_1342, \
                         msi0_1343, msi1_1341, msi1_1342, msi1_1343, msk_1717, msk_1718, \
                         msk_1719, msk_1720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2140[k] = f_4 * msi0_1341[k]
                    - f_5 * msi1_1341[k]
                    + f_3 * pc_x[k] * msk_1717[k];

        t_2141[k] = f_4 * msi0_1342[k]
                    - f_5 * msi1_1342[k]
                    + f_3 * pc_x[k] * msk_1718[k];

        t_2142[k] = f_4 * msi0_1343[k]
                    - f_5 * msi1_1343[k]
                    + f_3 * pc_x[k] * msk_1719[k];

        t_2143[k] = f_3 * pc_x[k] * msk_1720[k];
    }

#pragma omp simd aligned(t_2144, t_2145, t_2146, t_2147, t_2148, t_2149, t_2150, pc_x, \
                         msk_1721, msk_1722, msk_1723, msk_1724, msk_1725, msk_1726, \
                         msk_1727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2144[k] = f_3 * pc_x[k] * msk_1721[k];

        t_2145[k] = f_3 * pc_x[k] * msk_1722[k];

        t_2146[k] = f_3 * pc_x[k] * msk_1723[k];

        t_2147[k] = f_3 * pc_x[k] * msk_1724[k];

        t_2148[k] = f_3 * pc_x[k] * msk_1725[k];

        t_2149[k] = f_3 * pc_x[k] * msk_1726[k];

        t_2150[k] = f_3 * pc_x[k] * msk_1727[k];
    }

#pragma omp simd aligned(t_2151, t_2152, t_2153, pc_y, pc_z, lsk_1360, lsk_1396, lsk_1398, \
                         msi0_1337, msi0_1339, msi1_1337, msi1_1339, msk_1720, \
                         msk_1722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2151[k] = f_24 * lsk_1396[k]
                    + f_1 * msi0_1337[k]
                    - f_2 * msi1_1337[k]
                    + f_3 * pc_y[k] * msk_1720[k];

        t_2152[k] = f_16 * lsk_1360[k]
                    + f_3 * pc_z[k] * msk_1720[k];

        t_2153[k] = f_24 * lsk_1398[k]
                    + f_12 * msi0_1339[k]
                    - f_13 * msi1_1339[k]
                    + f_3 * pc_y[k] * msk_1722[k];
    }

#pragma omp simd aligned(t_2154, t_2155, t_2156, pc_y, lsk_1399, lsk_1400, lsk_1401, \
                         msi0_1340, msi0_1341, msi0_1342, msi1_1340, msi1_1341, msi1_1342, \
                         msk_1723, msk_1724, msk_1725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2154[k] = f_24 * lsk_1399[k]
                    + f_10 * msi0_1340[k]
                    - f_11 * msi1_1340[k]
                    + f_3 * pc_y[k] * msk_1723[k];

        t_2155[k] = f_24 * lsk_1400[k]
                    + f_8 * msi0_1341[k]
                    - f_9 * msi1_1341[k]
                    + f_3 * pc_y[k] * msk_1724[k];

        t_2156[k] = f_24 * lsk_1401[k]
                    + f_6 * msi0_1342[k]
                    - f_7 * msi1_1342[k]
                    + f_3 * pc_y[k] * msk_1725[k];
    }

#pragma omp simd aligned(t_2157, t_2158, t_2159, pc_y, pc_z, lsk_1367, lsk_1402, lsk_1403, \
                         msi0_1343, msi1_1343, msk_1726, msk_1727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2157[k] = f_24 * lsk_1402[k]
                    + f_4 * msi0_1343[k]
                    - f_5 * msi1_1343[k]
                    + f_3 * pc_y[k] * msk_1726[k];

        t_2158[k] = f_24 * lsk_1403[k]
                    + f_3 * pc_y[k] * msk_1727[k];

        t_2159[k] = f_16 * lsk_1367[k]
                    + f_1 * msi0_1343[k]
                    - f_2 * msi1_1343[k]
                    + f_3 * pc_z[k] * msk_1727[k];
    }

#pragma omp simd aligned(t_2160, t_2161, t_2162, pc_x, msi0_1344, msi0_1345, msi0_1346, \
                         msi1_1344, msi1_1345, msi1_1346, msk_1728, msk_1729, \
                         msk_1730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2160[k] = f_1 * msi0_1344[k]
                    - f_2 * msi1_1344[k]
                    + f_3 * pc_x[k] * msk_1728[k];

        t_2161[k] = f_22 * msi0_1345[k]
                    - f_23 * msi1_1345[k]
                    + f_3 * pc_x[k] * msk_1729[k];

        t_2162[k] = f_22 * msi0_1346[k]
                    - f_23 * msi1_1346[k]
                    + f_3 * pc_x[k] * msk_1730[k];
    }

#pragma omp simd aligned(t_2163, t_2164, t_2165, pc_x, msi0_1347, msi0_1348, msi0_1349, \
                         msi1_1347, msi1_1348, msi1_1349, msk_1731, msk_1732, \
                         msk_1733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2163[k] = f_12 * msi0_1347[k]
                    - f_13 * msi1_1347[k]
                    + f_3 * pc_x[k] * msk_1731[k];

        t_2164[k] = f_12 * msi0_1348[k]
                    - f_13 * msi1_1348[k]
                    + f_3 * pc_x[k] * msk_1732[k];

        t_2165[k] = f_12 * msi0_1349[k]
                    - f_13 * msi1_1349[k]
                    + f_3 * pc_x[k] * msk_1733[k];
    }

#pragma omp simd aligned(t_2166, t_2167, t_2168, pc_x, msi0_1350, msi0_1351, msi0_1352, \
                         msi1_1350, msi1_1351, msi1_1352, msk_1734, msk_1735, \
                         msk_1736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2166[k] = f_10 * msi0_1350[k]
                    - f_11 * msi1_1350[k]
                    + f_3 * pc_x[k] * msk_1734[k];

        t_2167[k] = f_10 * msi0_1351[k]
                    - f_11 * msi1_1351[k]
                    + f_3 * pc_x[k] * msk_1735[k];

        t_2168[k] = f_10 * msi0_1352[k]
                    - f_11 * msi1_1352[k]
                    + f_3 * pc_x[k] * msk_1736[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece19(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t lsk, const size_t msi0,
                                                           const size_t msi1, const size_t msk,
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
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

    auto *t_2169 = buffer.data(target + 2169);
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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk_1396 = buffer.data(lsk + 1396);
    const auto *lsk_1403 = buffer.data(lsk + 1403);
    const auto *lsk_1432 = buffer.data(lsk + 1432);
    const auto *lsk_1434 = buffer.data(lsk + 1434);
    const auto *lsk_1435 = buffer.data(lsk + 1435);
    const auto *lsk_1436 = buffer.data(lsk + 1436);
    const auto *lsk_1437 = buffer.data(lsk + 1437);
    const auto *lsk_1438 = buffer.data(lsk + 1438);
    const auto *lsk_1439 = buffer.data(lsk + 1439);
    const auto *lsk_1468 = buffer.data(lsk + 1468);
    const auto *lsk_1470 = buffer.data(lsk + 1470);
    const auto *lsk_1471 = buffer.data(lsk + 1471);
    const auto *lsk_1472 = buffer.data(lsk + 1472);
    const auto *lsk_1473 = buffer.data(lsk + 1473);
    const auto *lsk_1474 = buffer.data(lsk + 1474);
    const auto *lsk_1475 = buffer.data(lsk + 1475);

    const auto *msi0_1353 = buffer.data(msi0 + 1353);
    const auto *msi0_1354 = buffer.data(msi0 + 1354);
    const auto *msi0_1355 = buffer.data(msi0 + 1355);
    const auto *msi0_1356 = buffer.data(msi0 + 1356);
    const auto *msi0_1357 = buffer.data(msi0 + 1357);
    const auto *msi0_1358 = buffer.data(msi0 + 1358);
    const auto *msi0_1359 = buffer.data(msi0 + 1359);
    const auto *msi0_1360 = buffer.data(msi0 + 1360);
    const auto *msi0_1361 = buffer.data(msi0 + 1361);
    const auto *msi0_1362 = buffer.data(msi0 + 1362);
    const auto *msi0_1363 = buffer.data(msi0 + 1363);
    const auto *msi0_1364 = buffer.data(msi0 + 1364);
    const auto *msi0_1365 = buffer.data(msi0 + 1365);
    const auto *msi0_1366 = buffer.data(msi0 + 1366);
    const auto *msi0_1367 = buffer.data(msi0 + 1367);
    const auto *msi0_1368 = buffer.data(msi0 + 1368);
    const auto *msi0_1369 = buffer.data(msi0 + 1369);
    const auto *msi0_1370 = buffer.data(msi0 + 1370);
    const auto *msi0_1371 = buffer.data(msi0 + 1371);
    const auto *msi0_1372 = buffer.data(msi0 + 1372);
    const auto *msi0_1373 = buffer.data(msi0 + 1373);
    const auto *msi0_1374 = buffer.data(msi0 + 1374);
    const auto *msi0_1375 = buffer.data(msi0 + 1375);
    const auto *msi0_1376 = buffer.data(msi0 + 1376);
    const auto *msi0_1377 = buffer.data(msi0 + 1377);
    const auto *msi0_1378 = buffer.data(msi0 + 1378);
    const auto *msi0_1379 = buffer.data(msi0 + 1379);
    const auto *msi0_1380 = buffer.data(msi0 + 1380);
    const auto *msi0_1381 = buffer.data(msi0 + 1381);
    const auto *msi0_1382 = buffer.data(msi0 + 1382);
    const auto *msi0_1383 = buffer.data(msi0 + 1383);
    const auto *msi0_1384 = buffer.data(msi0 + 1384);
    const auto *msi0_1385 = buffer.data(msi0 + 1385);
    const auto *msi0_1386 = buffer.data(msi0 + 1386);
    const auto *msi0_1387 = buffer.data(msi0 + 1387);
    const auto *msi0_1388 = buffer.data(msi0 + 1388);
    const auto *msi0_1389 = buffer.data(msi0 + 1389);
    const auto *msi0_1390 = buffer.data(msi0 + 1390);
    const auto *msi0_1391 = buffer.data(msi0 + 1391);
    const auto *msi0_1392 = buffer.data(msi0 + 1392);
    const auto *msi0_1393 = buffer.data(msi0 + 1393);
    const auto *msi0_1394 = buffer.data(msi0 + 1394);
    const auto *msi0_1395 = buffer.data(msi0 + 1395);
    const auto *msi0_1396 = buffer.data(msi0 + 1396);
    const auto *msi0_1397 = buffer.data(msi0 + 1397);
    const auto *msi0_1398 = buffer.data(msi0 + 1398);
    const auto *msi0_1399 = buffer.data(msi0 + 1399);
    const auto *msi0_1400 = buffer.data(msi0 + 1400);
    const auto *msi0_1401 = buffer.data(msi0 + 1401);
    const auto *msi0_1402 = buffer.data(msi0 + 1402);
    const auto *msi0_1403 = buffer.data(msi0 + 1403);
    const auto *msi0_1404 = buffer.data(msi0 + 1404);
    const auto *msi0_1405 = buffer.data(msi0 + 1405);
    const auto *msi0_1406 = buffer.data(msi0 + 1406);
    const auto *msi0_1407 = buffer.data(msi0 + 1407);
    const auto *msi0_1408 = buffer.data(msi0 + 1408);
    const auto *msi0_1409 = buffer.data(msi0 + 1409);
    const auto *msi0_1410 = buffer.data(msi0 + 1410);
    const auto *msi0_1411 = buffer.data(msi0 + 1411);
    const auto *msi0_1412 = buffer.data(msi0 + 1412);
    const auto *msi0_1413 = buffer.data(msi0 + 1413);
    const auto *msi0_1414 = buffer.data(msi0 + 1414);
    const auto *msi0_1415 = buffer.data(msi0 + 1415);
    const auto *msi0_1416 = buffer.data(msi0 + 1416);
    const auto *msi0_1417 = buffer.data(msi0 + 1417);
    const auto *msi0_1418 = buffer.data(msi0 + 1418);
    const auto *msi0_1419 = buffer.data(msi0 + 1419);
    const auto *msi0_1420 = buffer.data(msi0 + 1420);
    const auto *msi0_1421 = buffer.data(msi0 + 1421);
    const auto *msi0_1422 = buffer.data(msi0 + 1422);
    const auto *msi0_1423 = buffer.data(msi0 + 1423);
    const auto *msi0_1424 = buffer.data(msi0 + 1424);
    const auto *msi0_1425 = buffer.data(msi0 + 1425);
    const auto *msi0_1426 = buffer.data(msi0 + 1426);
    const auto *msi0_1427 = buffer.data(msi0 + 1427);

    const auto *msi1_1353 = buffer.data(msi1 + 1353);
    const auto *msi1_1354 = buffer.data(msi1 + 1354);
    const auto *msi1_1355 = buffer.data(msi1 + 1355);
    const auto *msi1_1356 = buffer.data(msi1 + 1356);
    const auto *msi1_1357 = buffer.data(msi1 + 1357);
    const auto *msi1_1358 = buffer.data(msi1 + 1358);
    const auto *msi1_1359 = buffer.data(msi1 + 1359);
    const auto *msi1_1360 = buffer.data(msi1 + 1360);
    const auto *msi1_1361 = buffer.data(msi1 + 1361);
    const auto *msi1_1362 = buffer.data(msi1 + 1362);
    const auto *msi1_1363 = buffer.data(msi1 + 1363);
    const auto *msi1_1364 = buffer.data(msi1 + 1364);
    const auto *msi1_1365 = buffer.data(msi1 + 1365);
    const auto *msi1_1366 = buffer.data(msi1 + 1366);
    const auto *msi1_1367 = buffer.data(msi1 + 1367);
    const auto *msi1_1368 = buffer.data(msi1 + 1368);
    const auto *msi1_1369 = buffer.data(msi1 + 1369);
    const auto *msi1_1370 = buffer.data(msi1 + 1370);
    const auto *msi1_1371 = buffer.data(msi1 + 1371);
    const auto *msi1_1372 = buffer.data(msi1 + 1372);
    const auto *msi1_1373 = buffer.data(msi1 + 1373);
    const auto *msi1_1374 = buffer.data(msi1 + 1374);
    const auto *msi1_1375 = buffer.data(msi1 + 1375);
    const auto *msi1_1376 = buffer.data(msi1 + 1376);
    const auto *msi1_1377 = buffer.data(msi1 + 1377);
    const auto *msi1_1378 = buffer.data(msi1 + 1378);
    const auto *msi1_1379 = buffer.data(msi1 + 1379);
    const auto *msi1_1380 = buffer.data(msi1 + 1380);
    const auto *msi1_1381 = buffer.data(msi1 + 1381);
    const auto *msi1_1382 = buffer.data(msi1 + 1382);
    const auto *msi1_1383 = buffer.data(msi1 + 1383);
    const auto *msi1_1384 = buffer.data(msi1 + 1384);
    const auto *msi1_1385 = buffer.data(msi1 + 1385);
    const auto *msi1_1386 = buffer.data(msi1 + 1386);
    const auto *msi1_1387 = buffer.data(msi1 + 1387);
    const auto *msi1_1388 = buffer.data(msi1 + 1388);
    const auto *msi1_1389 = buffer.data(msi1 + 1389);
    const auto *msi1_1390 = buffer.data(msi1 + 1390);
    const auto *msi1_1391 = buffer.data(msi1 + 1391);
    const auto *msi1_1392 = buffer.data(msi1 + 1392);
    const auto *msi1_1393 = buffer.data(msi1 + 1393);
    const auto *msi1_1394 = buffer.data(msi1 + 1394);
    const auto *msi1_1395 = buffer.data(msi1 + 1395);
    const auto *msi1_1396 = buffer.data(msi1 + 1396);
    const auto *msi1_1397 = buffer.data(msi1 + 1397);
    const auto *msi1_1398 = buffer.data(msi1 + 1398);
    const auto *msi1_1399 = buffer.data(msi1 + 1399);
    const auto *msi1_1400 = buffer.data(msi1 + 1400);
    const auto *msi1_1401 = buffer.data(msi1 + 1401);
    const auto *msi1_1402 = buffer.data(msi1 + 1402);
    const auto *msi1_1403 = buffer.data(msi1 + 1403);
    const auto *msi1_1404 = buffer.data(msi1 + 1404);
    const auto *msi1_1405 = buffer.data(msi1 + 1405);
    const auto *msi1_1406 = buffer.data(msi1 + 1406);
    const auto *msi1_1407 = buffer.data(msi1 + 1407);
    const auto *msi1_1408 = buffer.data(msi1 + 1408);
    const auto *msi1_1409 = buffer.data(msi1 + 1409);
    const auto *msi1_1410 = buffer.data(msi1 + 1410);
    const auto *msi1_1411 = buffer.data(msi1 + 1411);
    const auto *msi1_1412 = buffer.data(msi1 + 1412);
    const auto *msi1_1413 = buffer.data(msi1 + 1413);
    const auto *msi1_1414 = buffer.data(msi1 + 1414);
    const auto *msi1_1415 = buffer.data(msi1 + 1415);
    const auto *msi1_1416 = buffer.data(msi1 + 1416);
    const auto *msi1_1417 = buffer.data(msi1 + 1417);
    const auto *msi1_1418 = buffer.data(msi1 + 1418);
    const auto *msi1_1419 = buffer.data(msi1 + 1419);
    const auto *msi1_1420 = buffer.data(msi1 + 1420);
    const auto *msi1_1421 = buffer.data(msi1 + 1421);
    const auto *msi1_1422 = buffer.data(msi1 + 1422);
    const auto *msi1_1423 = buffer.data(msi1 + 1423);
    const auto *msi1_1424 = buffer.data(msi1 + 1424);
    const auto *msi1_1425 = buffer.data(msi1 + 1425);
    const auto *msi1_1426 = buffer.data(msi1 + 1426);
    const auto *msi1_1427 = buffer.data(msi1 + 1427);

    const auto *msk_1737 = buffer.data(msk + 1737);
    const auto *msk_1738 = buffer.data(msk + 1738);
    const auto *msk_1739 = buffer.data(msk + 1739);
    const auto *msk_1740 = buffer.data(msk + 1740);
    const auto *msk_1741 = buffer.data(msk + 1741);
    const auto *msk_1742 = buffer.data(msk + 1742);
    const auto *msk_1743 = buffer.data(msk + 1743);
    const auto *msk_1744 = buffer.data(msk + 1744);
    const auto *msk_1745 = buffer.data(msk + 1745);
    const auto *msk_1746 = buffer.data(msk + 1746);
    const auto *msk_1747 = buffer.data(msk + 1747);
    const auto *msk_1748 = buffer.data(msk + 1748);
    const auto *msk_1749 = buffer.data(msk + 1749);
    const auto *msk_1750 = buffer.data(msk + 1750);
    const auto *msk_1751 = buffer.data(msk + 1751);
    const auto *msk_1752 = buffer.data(msk + 1752);
    const auto *msk_1753 = buffer.data(msk + 1753);
    const auto *msk_1754 = buffer.data(msk + 1754);
    const auto *msk_1755 = buffer.data(msk + 1755);
    const auto *msk_1756 = buffer.data(msk + 1756);
    const auto *msk_1757 = buffer.data(msk + 1757);
    const auto *msk_1758 = buffer.data(msk + 1758);
    const auto *msk_1759 = buffer.data(msk + 1759);
    const auto *msk_1760 = buffer.data(msk + 1760);
    const auto *msk_1761 = buffer.data(msk + 1761);
    const auto *msk_1762 = buffer.data(msk + 1762);
    const auto *msk_1763 = buffer.data(msk + 1763);
    const auto *msk_1764 = buffer.data(msk + 1764);
    const auto *msk_1765 = buffer.data(msk + 1765);
    const auto *msk_1766 = buffer.data(msk + 1766);
    const auto *msk_1767 = buffer.data(msk + 1767);
    const auto *msk_1768 = buffer.data(msk + 1768);
    const auto *msk_1769 = buffer.data(msk + 1769);
    const auto *msk_1770 = buffer.data(msk + 1770);
    const auto *msk_1771 = buffer.data(msk + 1771);
    const auto *msk_1772 = buffer.data(msk + 1772);
    const auto *msk_1773 = buffer.data(msk + 1773);
    const auto *msk_1774 = buffer.data(msk + 1774);
    const auto *msk_1775 = buffer.data(msk + 1775);
    const auto *msk_1776 = buffer.data(msk + 1776);
    const auto *msk_1777 = buffer.data(msk + 1777);
    const auto *msk_1778 = buffer.data(msk + 1778);
    const auto *msk_1779 = buffer.data(msk + 1779);
    const auto *msk_1780 = buffer.data(msk + 1780);
    const auto *msk_1781 = buffer.data(msk + 1781);
    const auto *msk_1782 = buffer.data(msk + 1782);
    const auto *msk_1783 = buffer.data(msk + 1783);
    const auto *msk_1784 = buffer.data(msk + 1784);
    const auto *msk_1785 = buffer.data(msk + 1785);
    const auto *msk_1786 = buffer.data(msk + 1786);
    const auto *msk_1787 = buffer.data(msk + 1787);
    const auto *msk_1788 = buffer.data(msk + 1788);
    const auto *msk_1789 = buffer.data(msk + 1789);
    const auto *msk_1790 = buffer.data(msk + 1790);
    const auto *msk_1791 = buffer.data(msk + 1791);
    const auto *msk_1792 = buffer.data(msk + 1792);
    const auto *msk_1793 = buffer.data(msk + 1793);
    const auto *msk_1794 = buffer.data(msk + 1794);
    const auto *msk_1795 = buffer.data(msk + 1795);
    const auto *msk_1796 = buffer.data(msk + 1796);
    const auto *msk_1797 = buffer.data(msk + 1797);
    const auto *msk_1798 = buffer.data(msk + 1798);
    const auto *msk_1799 = buffer.data(msk + 1799);
    const auto *msk_1800 = buffer.data(msk + 1800);
    const auto *msk_1801 = buffer.data(msk + 1801);
    const auto *msk_1802 = buffer.data(msk + 1802);
    const auto *msk_1803 = buffer.data(msk + 1803);
    const auto *msk_1804 = buffer.data(msk + 1804);
    const auto *msk_1805 = buffer.data(msk + 1805);
    const auto *msk_1806 = buffer.data(msk + 1806);
    const auto *msk_1807 = buffer.data(msk + 1807);
    const auto *msk_1808 = buffer.data(msk + 1808);
    const auto *msk_1809 = buffer.data(msk + 1809);
    const auto *msk_1810 = buffer.data(msk + 1810);
    const auto *msk_1811 = buffer.data(msk + 1811);
    const auto *msk_1812 = buffer.data(msk + 1812);
    const auto *msk_1813 = buffer.data(msk + 1813);
    const auto *msk_1814 = buffer.data(msk + 1814);
    const auto *msk_1815 = buffer.data(msk + 1815);
    const auto *msk_1816 = buffer.data(msk + 1816);
    const auto *msk_1817 = buffer.data(msk + 1817);
    const auto *msk_1818 = buffer.data(msk + 1818);
    const auto *msk_1819 = buffer.data(msk + 1819);
    const auto *msk_1820 = buffer.data(msk + 1820);
    const auto *msk_1821 = buffer.data(msk + 1821);
    const auto *msk_1822 = buffer.data(msk + 1822);
    const auto *msk_1823 = buffer.data(msk + 1823);
    const auto *msk_1824 = buffer.data(msk + 1824);
    const auto *msk_1825 = buffer.data(msk + 1825);
    const auto *msk_1826 = buffer.data(msk + 1826);
    const auto *msk_1827 = buffer.data(msk + 1827);
    const auto *msk_1828 = buffer.data(msk + 1828);
    const auto *msk_1829 = buffer.data(msk + 1829);
    const auto *msk_1830 = buffer.data(msk + 1830);
    const auto *msk_1831 = buffer.data(msk + 1831);
    const auto *msk_1832 = buffer.data(msk + 1832);

#pragma omp simd aligned(t_2169, t_2170, t_2171, pc_x, msi0_1353, msi0_1354, msi0_1355, \
                         msi1_1353, msi1_1354, msi1_1355, msk_1737, msk_1738, \
                         msk_1739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2169[k] = f_10 * msi0_1353[k]
                    - f_11 * msi1_1353[k]
                    + f_3 * pc_x[k] * msk_1737[k];

        t_2170[k] = f_8 * msi0_1354[k]
                    - f_9 * msi1_1354[k]
                    + f_3 * pc_x[k] * msk_1738[k];

        t_2171[k] = f_8 * msi0_1355[k]
                    - f_9 * msi1_1355[k]
                    + f_3 * pc_x[k] * msk_1739[k];
    }

#pragma omp simd aligned(t_2172, t_2173, t_2174, pc_x, msi0_1356, msi0_1357, msi0_1358, \
                         msi1_1356, msi1_1357, msi1_1358, msk_1740, msk_1741, \
                         msk_1742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2172[k] = f_8 * msi0_1356[k]
                    - f_9 * msi1_1356[k]
                    + f_3 * pc_x[k] * msk_1740[k];

        t_2173[k] = f_8 * msi0_1357[k]
                    - f_9 * msi1_1357[k]
                    + f_3 * pc_x[k] * msk_1741[k];

        t_2174[k] = f_8 * msi0_1358[k]
                    - f_9 * msi1_1358[k]
                    + f_3 * pc_x[k] * msk_1742[k];
    }

#pragma omp simd aligned(t_2175, t_2176, t_2177, pc_x, msi0_1359, msi0_1360, msi0_1361, \
                         msi1_1359, msi1_1360, msi1_1361, msk_1743, msk_1744, \
                         msk_1745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2175[k] = f_6 * msi0_1359[k]
                    - f_7 * msi1_1359[k]
                    + f_3 * pc_x[k] * msk_1743[k];

        t_2176[k] = f_6 * msi0_1360[k]
                    - f_7 * msi1_1360[k]
                    + f_3 * pc_x[k] * msk_1744[k];

        t_2177[k] = f_6 * msi0_1361[k]
                    - f_7 * msi1_1361[k]
                    + f_3 * pc_x[k] * msk_1745[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, pc_x, msi0_1362, msi0_1363, msi0_1364, \
                         msi1_1362, msi1_1363, msi1_1364, msk_1746, msk_1747, \
                         msk_1748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = f_6 * msi0_1362[k]
                    - f_7 * msi1_1362[k]
                    + f_3 * pc_x[k] * msk_1746[k];

        t_2179[k] = f_6 * msi0_1363[k]
                    - f_7 * msi1_1363[k]
                    + f_3 * pc_x[k] * msk_1747[k];

        t_2180[k] = f_6 * msi0_1364[k]
                    - f_7 * msi1_1364[k]
                    + f_3 * pc_x[k] * msk_1748[k];
    }

#pragma omp simd aligned(t_2181, t_2182, t_2183, pc_x, msi0_1365, msi0_1366, msi0_1367, \
                         msi1_1365, msi1_1366, msi1_1367, msk_1749, msk_1750, \
                         msk_1751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2181[k] = f_4 * msi0_1365[k]
                    - f_5 * msi1_1365[k]
                    + f_3 * pc_x[k] * msk_1749[k];

        t_2182[k] = f_4 * msi0_1366[k]
                    - f_5 * msi1_1366[k]
                    + f_3 * pc_x[k] * msk_1750[k];

        t_2183[k] = f_4 * msi0_1367[k]
                    - f_5 * msi1_1367[k]
                    + f_3 * pc_x[k] * msk_1751[k];
    }

#pragma omp simd aligned(t_2184, t_2185, t_2186, pc_x, msi0_1368, msi0_1369, msi0_1370, \
                         msi1_1368, msi1_1369, msi1_1370, msk_1752, msk_1753, \
                         msk_1754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2184[k] = f_4 * msi0_1368[k]
                    - f_5 * msi1_1368[k]
                    + f_3 * pc_x[k] * msk_1752[k];

        t_2185[k] = f_4 * msi0_1369[k]
                    - f_5 * msi1_1369[k]
                    + f_3 * pc_x[k] * msk_1753[k];

        t_2186[k] = f_4 * msi0_1370[k]
                    - f_5 * msi1_1370[k]
                    + f_3 * pc_x[k] * msk_1754[k];
    }

#pragma omp simd aligned(t_2187, t_2188, t_2189, t_2190, t_2191, t_2192, pc_x, msi0_1371, \
                         msi1_1371, msk_1755, msk_1756, msk_1757, msk_1758, msk_1759, \
                         msk_1760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2187[k] = f_4 * msi0_1371[k]
                    - f_5 * msi1_1371[k]
                    + f_3 * pc_x[k] * msk_1755[k];

        t_2188[k] = f_3 * pc_x[k] * msk_1756[k];

        t_2189[k] = f_3 * pc_x[k] * msk_1757[k];

        t_2190[k] = f_3 * pc_x[k] * msk_1758[k];

        t_2191[k] = f_3 * pc_x[k] * msk_1759[k];

        t_2192[k] = f_3 * pc_x[k] * msk_1760[k];
    }

#pragma omp simd aligned(t_2193, t_2194, t_2195, t_2196, t_2197, pc_x, pc_y, pc_z, lsk_1396, \
                         lsk_1432, msi0_1365, msi1_1365, msk_1756, msk_1761, msk_1762, \
                         msk_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2193[k] = f_3 * pc_x[k] * msk_1761[k];

        t_2194[k] = f_3 * pc_x[k] * msk_1762[k];

        t_2195[k] = f_3 * pc_x[k] * msk_1763[k];

        t_2196[k] = f_20 * lsk_1432[k]
                    + f_1 * msi0_1365[k]
                    - f_2 * msi1_1365[k]
                    + f_3 * pc_y[k] * msk_1756[k];

        t_2197[k] = f_17 * lsk_1396[k]
                    + f_3 * pc_z[k] * msk_1756[k];
    }

#pragma omp simd aligned(t_2198, t_2199, t_2200, pc_y, lsk_1434, lsk_1435, lsk_1436, \
                         msi0_1367, msi0_1368, msi0_1369, msi1_1367, msi1_1368, msi1_1369, \
                         msk_1758, msk_1759, msk_1760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2198[k] = f_20 * lsk_1434[k]
                    + f_12 * msi0_1367[k]
                    - f_13 * msi1_1367[k]
                    + f_3 * pc_y[k] * msk_1758[k];

        t_2199[k] = f_20 * lsk_1435[k]
                    + f_10 * msi0_1368[k]
                    - f_11 * msi1_1368[k]
                    + f_3 * pc_y[k] * msk_1759[k];

        t_2200[k] = f_20 * lsk_1436[k]
                    + f_8 * msi0_1369[k]
                    - f_9 * msi1_1369[k]
                    + f_3 * pc_y[k] * msk_1760[k];
    }

#pragma omp simd aligned(t_2201, t_2202, t_2203, pc_y, lsk_1437, lsk_1438, lsk_1439, \
                         msi0_1370, msi0_1371, msi1_1370, msi1_1371, msk_1761, msk_1762, \
                         msk_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2201[k] = f_20 * lsk_1437[k]
                    + f_6 * msi0_1370[k]
                    - f_7 * msi1_1370[k]
                    + f_3 * pc_y[k] * msk_1761[k];

        t_2202[k] = f_20 * lsk_1438[k]
                    + f_4 * msi0_1371[k]
                    - f_5 * msi1_1371[k]
                    + f_3 * pc_y[k] * msk_1762[k];

        t_2203[k] = f_20 * lsk_1439[k]
                    + f_3 * pc_y[k] * msk_1763[k];
    }

#pragma omp simd aligned(t_2204, t_2205, t_2206, pc_x, pc_z, lsk_1403, msi0_1371, msi0_1372, \
                         msi0_1373, msi1_1371, msi1_1372, msi1_1373, msk_1763, msk_1764, \
                         msk_1765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2204[k] = f_17 * lsk_1403[k]
                    + f_1 * msi0_1371[k]
                    - f_2 * msi1_1371[k]
                    + f_3 * pc_z[k] * msk_1763[k];

        t_2205[k] = f_1 * msi0_1372[k]
                    - f_2 * msi1_1372[k]
                    + f_3 * pc_x[k] * msk_1764[k];

        t_2206[k] = f_22 * msi0_1373[k]
                    - f_23 * msi1_1373[k]
                    + f_3 * pc_x[k] * msk_1765[k];
    }

#pragma omp simd aligned(t_2207, t_2208, t_2209, pc_x, msi0_1374, msi0_1375, msi0_1376, \
                         msi1_1374, msi1_1375, msi1_1376, msk_1766, msk_1767, \
                         msk_1768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2207[k] = f_22 * msi0_1374[k]
                    - f_23 * msi1_1374[k]
                    + f_3 * pc_x[k] * msk_1766[k];

        t_2208[k] = f_12 * msi0_1375[k]
                    - f_13 * msi1_1375[k]
                    + f_3 * pc_x[k] * msk_1767[k];

        t_2209[k] = f_12 * msi0_1376[k]
                    - f_13 * msi1_1376[k]
                    + f_3 * pc_x[k] * msk_1768[k];
    }

#pragma omp simd aligned(t_2210, t_2211, t_2212, pc_x, msi0_1377, msi0_1378, msi0_1379, \
                         msi1_1377, msi1_1378, msi1_1379, msk_1769, msk_1770, \
                         msk_1771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2210[k] = f_12 * msi0_1377[k]
                    - f_13 * msi1_1377[k]
                    + f_3 * pc_x[k] * msk_1769[k];

        t_2211[k] = f_10 * msi0_1378[k]
                    - f_11 * msi1_1378[k]
                    + f_3 * pc_x[k] * msk_1770[k];

        t_2212[k] = f_10 * msi0_1379[k]
                    - f_11 * msi1_1379[k]
                    + f_3 * pc_x[k] * msk_1771[k];
    }

#pragma omp simd aligned(t_2213, t_2214, t_2215, pc_x, msi0_1380, msi0_1381, msi0_1382, \
                         msi1_1380, msi1_1381, msi1_1382, msk_1772, msk_1773, \
                         msk_1774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2213[k] = f_10 * msi0_1380[k]
                    - f_11 * msi1_1380[k]
                    + f_3 * pc_x[k] * msk_1772[k];

        t_2214[k] = f_10 * msi0_1381[k]
                    - f_11 * msi1_1381[k]
                    + f_3 * pc_x[k] * msk_1773[k];

        t_2215[k] = f_8 * msi0_1382[k]
                    - f_9 * msi1_1382[k]
                    + f_3 * pc_x[k] * msk_1774[k];
    }

#pragma omp simd aligned(t_2216, t_2217, t_2218, pc_x, msi0_1383, msi0_1384, msi0_1385, \
                         msi1_1383, msi1_1384, msi1_1385, msk_1775, msk_1776, \
                         msk_1777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2216[k] = f_8 * msi0_1383[k]
                    - f_9 * msi1_1383[k]
                    + f_3 * pc_x[k] * msk_1775[k];

        t_2217[k] = f_8 * msi0_1384[k]
                    - f_9 * msi1_1384[k]
                    + f_3 * pc_x[k] * msk_1776[k];

        t_2218[k] = f_8 * msi0_1385[k]
                    - f_9 * msi1_1385[k]
                    + f_3 * pc_x[k] * msk_1777[k];
    }

#pragma omp simd aligned(t_2219, t_2220, t_2221, pc_x, msi0_1386, msi0_1387, msi0_1388, \
                         msi1_1386, msi1_1387, msi1_1388, msk_1778, msk_1779, \
                         msk_1780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2219[k] = f_8 * msi0_1386[k]
                    - f_9 * msi1_1386[k]
                    + f_3 * pc_x[k] * msk_1778[k];

        t_2220[k] = f_6 * msi0_1387[k]
                    - f_7 * msi1_1387[k]
                    + f_3 * pc_x[k] * msk_1779[k];

        t_2221[k] = f_6 * msi0_1388[k]
                    - f_7 * msi1_1388[k]
                    + f_3 * pc_x[k] * msk_1780[k];
    }

#pragma omp simd aligned(t_2222, t_2223, t_2224, pc_x, msi0_1389, msi0_1390, msi0_1391, \
                         msi1_1389, msi1_1390, msi1_1391, msk_1781, msk_1782, \
                         msk_1783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2222[k] = f_6 * msi0_1389[k]
                    - f_7 * msi1_1389[k]
                    + f_3 * pc_x[k] * msk_1781[k];

        t_2223[k] = f_6 * msi0_1390[k]
                    - f_7 * msi1_1390[k]
                    + f_3 * pc_x[k] * msk_1782[k];

        t_2224[k] = f_6 * msi0_1391[k]
                    - f_7 * msi1_1391[k]
                    + f_3 * pc_x[k] * msk_1783[k];
    }

#pragma omp simd aligned(t_2225, t_2226, t_2227, pc_x, msi0_1392, msi0_1393, msi0_1394, \
                         msi1_1392, msi1_1393, msi1_1394, msk_1784, msk_1785, \
                         msk_1786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2225[k] = f_6 * msi0_1392[k]
                    - f_7 * msi1_1392[k]
                    + f_3 * pc_x[k] * msk_1784[k];

        t_2226[k] = f_4 * msi0_1393[k]
                    - f_5 * msi1_1393[k]
                    + f_3 * pc_x[k] * msk_1785[k];

        t_2227[k] = f_4 * msi0_1394[k]
                    - f_5 * msi1_1394[k]
                    + f_3 * pc_x[k] * msk_1786[k];
    }

#pragma omp simd aligned(t_2228, t_2229, t_2230, pc_x, msi0_1395, msi0_1396, msi0_1397, \
                         msi1_1395, msi1_1396, msi1_1397, msk_1787, msk_1788, \
                         msk_1789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2228[k] = f_4 * msi0_1395[k]
                    - f_5 * msi1_1395[k]
                    + f_3 * pc_x[k] * msk_1787[k];

        t_2229[k] = f_4 * msi0_1396[k]
                    - f_5 * msi1_1396[k]
                    + f_3 * pc_x[k] * msk_1788[k];

        t_2230[k] = f_4 * msi0_1397[k]
                    - f_5 * msi1_1397[k]
                    + f_3 * pc_x[k] * msk_1789[k];
    }

#pragma omp simd aligned(t_2231, t_2232, t_2233, t_2234, t_2235, pc_x, msi0_1398, msi0_1399, \
                         msi1_1398, msi1_1399, msk_1790, msk_1791, msk_1792, msk_1793, \
                         msk_1794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2231[k] = f_4 * msi0_1398[k]
                    - f_5 * msi1_1398[k]
                    + f_3 * pc_x[k] * msk_1790[k];

        t_2232[k] = f_4 * msi0_1399[k]
                    - f_5 * msi1_1399[k]
                    + f_3 * pc_x[k] * msk_1791[k];

        t_2233[k] = f_3 * pc_x[k] * msk_1792[k];

        t_2234[k] = f_3 * pc_x[k] * msk_1793[k];

        t_2235[k] = f_3 * pc_x[k] * msk_1794[k];
    }

#pragma omp simd aligned(t_2236, t_2237, t_2238, t_2239, t_2240, pc_x, msk_1795, msk_1796, \
                         msk_1797, msk_1798, msk_1799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2236[k] = f_3 * pc_x[k] * msk_1795[k];

        t_2237[k] = f_3 * pc_x[k] * msk_1796[k];

        t_2238[k] = f_3 * pc_x[k] * msk_1797[k];

        t_2239[k] = f_3 * pc_x[k] * msk_1798[k];

        t_2240[k] = f_3 * pc_x[k] * msk_1799[k];
    }

#pragma omp simd aligned(t_2241, t_2242, t_2243, pc_y, pc_z, lsk_1432, lsk_1468, lsk_1470, \
                         msi0_1393, msi0_1395, msi1_1393, msi1_1395, msk_1792, \
                         msk_1794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2241[k] = f_19 * lsk_1468[k]
                    + f_1 * msi0_1393[k]
                    - f_2 * msi1_1393[k]
                    + f_3 * pc_y[k] * msk_1792[k];

        t_2242[k] = f_18 * lsk_1432[k]
                    + f_3 * pc_z[k] * msk_1792[k];

        t_2243[k] = f_19 * lsk_1470[k]
                    + f_12 * msi0_1395[k]
                    - f_13 * msi1_1395[k]
                    + f_3 * pc_y[k] * msk_1794[k];
    }

#pragma omp simd aligned(t_2244, t_2245, t_2246, pc_y, lsk_1471, lsk_1472, lsk_1473, \
                         msi0_1396, msi0_1397, msi0_1398, msi1_1396, msi1_1397, msi1_1398, \
                         msk_1795, msk_1796, msk_1797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2244[k] = f_19 * lsk_1471[k]
                    + f_10 * msi0_1396[k]
                    - f_11 * msi1_1396[k]
                    + f_3 * pc_y[k] * msk_1795[k];

        t_2245[k] = f_19 * lsk_1472[k]
                    + f_8 * msi0_1397[k]
                    - f_9 * msi1_1397[k]
                    + f_3 * pc_y[k] * msk_1796[k];

        t_2246[k] = f_19 * lsk_1473[k]
                    + f_6 * msi0_1398[k]
                    - f_7 * msi1_1398[k]
                    + f_3 * pc_y[k] * msk_1797[k];
    }

#pragma omp simd aligned(t_2247, t_2248, t_2249, pc_y, pc_z, lsk_1439, lsk_1474, lsk_1475, \
                         msi0_1399, msi1_1399, msk_1798, msk_1799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2247[k] = f_19 * lsk_1474[k]
                    + f_4 * msi0_1399[k]
                    - f_5 * msi1_1399[k]
                    + f_3 * pc_y[k] * msk_1798[k];

        t_2248[k] = f_19 * lsk_1475[k]
                    + f_3 * pc_y[k] * msk_1799[k];

        t_2249[k] = f_18 * lsk_1439[k]
                    + f_1 * msi0_1399[k]
                    - f_2 * msi1_1399[k]
                    + f_3 * pc_z[k] * msk_1799[k];
    }

#pragma omp simd aligned(t_2250, t_2251, t_2252, pc_x, msi0_1400, msi0_1401, msi0_1402, \
                         msi1_1400, msi1_1401, msi1_1402, msk_1800, msk_1801, \
                         msk_1802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2250[k] = f_1 * msi0_1400[k]
                    - f_2 * msi1_1400[k]
                    + f_3 * pc_x[k] * msk_1800[k];

        t_2251[k] = f_22 * msi0_1401[k]
                    - f_23 * msi1_1401[k]
                    + f_3 * pc_x[k] * msk_1801[k];

        t_2252[k] = f_22 * msi0_1402[k]
                    - f_23 * msi1_1402[k]
                    + f_3 * pc_x[k] * msk_1802[k];
    }

#pragma omp simd aligned(t_2253, t_2254, t_2255, pc_x, msi0_1403, msi0_1404, msi0_1405, \
                         msi1_1403, msi1_1404, msi1_1405, msk_1803, msk_1804, \
                         msk_1805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2253[k] = f_12 * msi0_1403[k]
                    - f_13 * msi1_1403[k]
                    + f_3 * pc_x[k] * msk_1803[k];

        t_2254[k] = f_12 * msi0_1404[k]
                    - f_13 * msi1_1404[k]
                    + f_3 * pc_x[k] * msk_1804[k];

        t_2255[k] = f_12 * msi0_1405[k]
                    - f_13 * msi1_1405[k]
                    + f_3 * pc_x[k] * msk_1805[k];
    }

#pragma omp simd aligned(t_2256, t_2257, t_2258, pc_x, msi0_1406, msi0_1407, msi0_1408, \
                         msi1_1406, msi1_1407, msi1_1408, msk_1806, msk_1807, \
                         msk_1808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2256[k] = f_10 * msi0_1406[k]
                    - f_11 * msi1_1406[k]
                    + f_3 * pc_x[k] * msk_1806[k];

        t_2257[k] = f_10 * msi0_1407[k]
                    - f_11 * msi1_1407[k]
                    + f_3 * pc_x[k] * msk_1807[k];

        t_2258[k] = f_10 * msi0_1408[k]
                    - f_11 * msi1_1408[k]
                    + f_3 * pc_x[k] * msk_1808[k];
    }

#pragma omp simd aligned(t_2259, t_2260, t_2261, pc_x, msi0_1409, msi0_1410, msi0_1411, \
                         msi1_1409, msi1_1410, msi1_1411, msk_1809, msk_1810, \
                         msk_1811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2259[k] = f_10 * msi0_1409[k]
                    - f_11 * msi1_1409[k]
                    + f_3 * pc_x[k] * msk_1809[k];

        t_2260[k] = f_8 * msi0_1410[k]
                    - f_9 * msi1_1410[k]
                    + f_3 * pc_x[k] * msk_1810[k];

        t_2261[k] = f_8 * msi0_1411[k]
                    - f_9 * msi1_1411[k]
                    + f_3 * pc_x[k] * msk_1811[k];
    }

#pragma omp simd aligned(t_2262, t_2263, t_2264, pc_x, msi0_1412, msi0_1413, msi0_1414, \
                         msi1_1412, msi1_1413, msi1_1414, msk_1812, msk_1813, \
                         msk_1814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2262[k] = f_8 * msi0_1412[k]
                    - f_9 * msi1_1412[k]
                    + f_3 * pc_x[k] * msk_1812[k];

        t_2263[k] = f_8 * msi0_1413[k]
                    - f_9 * msi1_1413[k]
                    + f_3 * pc_x[k] * msk_1813[k];

        t_2264[k] = f_8 * msi0_1414[k]
                    - f_9 * msi1_1414[k]
                    + f_3 * pc_x[k] * msk_1814[k];
    }

#pragma omp simd aligned(t_2265, t_2266, t_2267, pc_x, msi0_1415, msi0_1416, msi0_1417, \
                         msi1_1415, msi1_1416, msi1_1417, msk_1815, msk_1816, \
                         msk_1817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2265[k] = f_6 * msi0_1415[k]
                    - f_7 * msi1_1415[k]
                    + f_3 * pc_x[k] * msk_1815[k];

        t_2266[k] = f_6 * msi0_1416[k]
                    - f_7 * msi1_1416[k]
                    + f_3 * pc_x[k] * msk_1816[k];

        t_2267[k] = f_6 * msi0_1417[k]
                    - f_7 * msi1_1417[k]
                    + f_3 * pc_x[k] * msk_1817[k];
    }

#pragma omp simd aligned(t_2268, t_2269, t_2270, pc_x, msi0_1418, msi0_1419, msi0_1420, \
                         msi1_1418, msi1_1419, msi1_1420, msk_1818, msk_1819, \
                         msk_1820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2268[k] = f_6 * msi0_1418[k]
                    - f_7 * msi1_1418[k]
                    + f_3 * pc_x[k] * msk_1818[k];

        t_2269[k] = f_6 * msi0_1419[k]
                    - f_7 * msi1_1419[k]
                    + f_3 * pc_x[k] * msk_1819[k];

        t_2270[k] = f_6 * msi0_1420[k]
                    - f_7 * msi1_1420[k]
                    + f_3 * pc_x[k] * msk_1820[k];
    }

#pragma omp simd aligned(t_2271, t_2272, t_2273, pc_x, msi0_1421, msi0_1422, msi0_1423, \
                         msi1_1421, msi1_1422, msi1_1423, msk_1821, msk_1822, \
                         msk_1823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2271[k] = f_4 * msi0_1421[k]
                    - f_5 * msi1_1421[k]
                    + f_3 * pc_x[k] * msk_1821[k];

        t_2272[k] = f_4 * msi0_1422[k]
                    - f_5 * msi1_1422[k]
                    + f_3 * pc_x[k] * msk_1822[k];

        t_2273[k] = f_4 * msi0_1423[k]
                    - f_5 * msi1_1423[k]
                    + f_3 * pc_x[k] * msk_1823[k];
    }

#pragma omp simd aligned(t_2274, t_2275, t_2276, pc_x, msi0_1424, msi0_1425, msi0_1426, \
                         msi1_1424, msi1_1425, msi1_1426, msk_1824, msk_1825, \
                         msk_1826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2274[k] = f_4 * msi0_1424[k]
                    - f_5 * msi1_1424[k]
                    + f_3 * pc_x[k] * msk_1824[k];

        t_2275[k] = f_4 * msi0_1425[k]
                    - f_5 * msi1_1425[k]
                    + f_3 * pc_x[k] * msk_1825[k];

        t_2276[k] = f_4 * msi0_1426[k]
                    - f_5 * msi1_1426[k]
                    + f_3 * pc_x[k] * msk_1826[k];
    }

#pragma omp simd aligned(t_2277, t_2278, t_2279, t_2280, t_2281, t_2282, pc_x, msi0_1427, \
                         msi1_1427, msk_1827, msk_1828, msk_1829, msk_1830, msk_1831, \
                         msk_1832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2277[k] = f_4 * msi0_1427[k]
                    - f_5 * msi1_1427[k]
                    + f_3 * pc_x[k] * msk_1827[k];

        t_2278[k] = f_3 * pc_x[k] * msk_1828[k];

        t_2279[k] = f_3 * pc_x[k] * msk_1829[k];

        t_2280[k] = f_3 * pc_x[k] * msk_1830[k];

        t_2281[k] = f_3 * pc_x[k] * msk_1831[k];

        t_2282[k] = f_3 * pc_x[k] * msk_1832[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece20(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
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
    const auto f_24 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1980 = buffer.data(lsl0 + 1980);
    const auto *lsl0_1982 = buffer.data(lsl0 + 1982);
    const auto *lsl0_1985 = buffer.data(lsl0 + 1985);
    const auto *lsl0_1989 = buffer.data(lsl0 + 1989);

    const auto *lsk_1468 = buffer.data(lsk + 1468);
    const auto *lsk_1475 = buffer.data(lsk + 1475);
    const auto *lsk_1504 = buffer.data(lsk + 1504);
    const auto *lsk_1506 = buffer.data(lsk + 1506);
    const auto *lsk_1507 = buffer.data(lsk + 1507);
    const auto *lsk_1508 = buffer.data(lsk + 1508);
    const auto *lsk_1509 = buffer.data(lsk + 1509);
    const auto *lsk_1510 = buffer.data(lsk + 1510);
    const auto *lsk_1511 = buffer.data(lsk + 1511);
    const auto *lsk_1540 = buffer.data(lsk + 1540);
    const auto *lsk_1542 = buffer.data(lsk + 1542);
    const auto *lsk_1543 = buffer.data(lsk + 1543);
    const auto *lsk_1544 = buffer.data(lsk + 1544);
    const auto *lsk_1545 = buffer.data(lsk + 1545);
    const auto *lsk_1546 = buffer.data(lsk + 1546);
    const auto *lsk_1547 = buffer.data(lsk + 1547);
    const auto *lsk_1576 = buffer.data(lsk + 1576);
    const auto *lsk_1578 = buffer.data(lsk + 1578);
    const auto *lsk_1579 = buffer.data(lsk + 1579);
    const auto *lsk_1580 = buffer.data(lsk + 1580);
    const auto *lsk_1581 = buffer.data(lsk + 1581);
    const auto *lsk_1582 = buffer.data(lsk + 1582);
    const auto *lsk_1583 = buffer.data(lsk + 1583);

    const auto *lsl1_1980 = buffer.data(lsl1 + 1980);
    const auto *lsl1_1982 = buffer.data(lsl1 + 1982);
    const auto *lsl1_1985 = buffer.data(lsl1 + 1985);
    const auto *lsl1_1989 = buffer.data(lsl1 + 1989);

    const auto *msi0_1421 = buffer.data(msi0 + 1421);
    const auto *msi0_1423 = buffer.data(msi0 + 1423);
    const auto *msi0_1424 = buffer.data(msi0 + 1424);
    const auto *msi0_1425 = buffer.data(msi0 + 1425);
    const auto *msi0_1426 = buffer.data(msi0 + 1426);
    const auto *msi0_1427 = buffer.data(msi0 + 1427);
    const auto *msi0_1428 = buffer.data(msi0 + 1428);
    const auto *msi0_1429 = buffer.data(msi0 + 1429);
    const auto *msi0_1430 = buffer.data(msi0 + 1430);
    const auto *msi0_1431 = buffer.data(msi0 + 1431);
    const auto *msi0_1432 = buffer.data(msi0 + 1432);
    const auto *msi0_1433 = buffer.data(msi0 + 1433);
    const auto *msi0_1434 = buffer.data(msi0 + 1434);
    const auto *msi0_1435 = buffer.data(msi0 + 1435);
    const auto *msi0_1436 = buffer.data(msi0 + 1436);
    const auto *msi0_1437 = buffer.data(msi0 + 1437);
    const auto *msi0_1438 = buffer.data(msi0 + 1438);
    const auto *msi0_1439 = buffer.data(msi0 + 1439);
    const auto *msi0_1440 = buffer.data(msi0 + 1440);
    const auto *msi0_1441 = buffer.data(msi0 + 1441);
    const auto *msi0_1442 = buffer.data(msi0 + 1442);
    const auto *msi0_1443 = buffer.data(msi0 + 1443);
    const auto *msi0_1444 = buffer.data(msi0 + 1444);
    const auto *msi0_1445 = buffer.data(msi0 + 1445);
    const auto *msi0_1446 = buffer.data(msi0 + 1446);
    const auto *msi0_1447 = buffer.data(msi0 + 1447);
    const auto *msi0_1448 = buffer.data(msi0 + 1448);
    const auto *msi0_1449 = buffer.data(msi0 + 1449);
    const auto *msi0_1450 = buffer.data(msi0 + 1450);
    const auto *msi0_1451 = buffer.data(msi0 + 1451);
    const auto *msi0_1452 = buffer.data(msi0 + 1452);
    const auto *msi0_1453 = buffer.data(msi0 + 1453);
    const auto *msi0_1454 = buffer.data(msi0 + 1454);
    const auto *msi0_1455 = buffer.data(msi0 + 1455);
    const auto *msi0_1456 = buffer.data(msi0 + 1456);
    const auto *msi0_1457 = buffer.data(msi0 + 1457);
    const auto *msi0_1458 = buffer.data(msi0 + 1458);
    const auto *msi0_1459 = buffer.data(msi0 + 1459);
    const auto *msi0_1460 = buffer.data(msi0 + 1460);
    const auto *msi0_1461 = buffer.data(msi0 + 1461);
    const auto *msi0_1462 = buffer.data(msi0 + 1462);
    const auto *msi0_1463 = buffer.data(msi0 + 1463);
    const auto *msi0_1464 = buffer.data(msi0 + 1464);
    const auto *msi0_1465 = buffer.data(msi0 + 1465);
    const auto *msi0_1466 = buffer.data(msi0 + 1466);
    const auto *msi0_1467 = buffer.data(msi0 + 1467);
    const auto *msi0_1468 = buffer.data(msi0 + 1468);
    const auto *msi0_1469 = buffer.data(msi0 + 1469);
    const auto *msi0_1470 = buffer.data(msi0 + 1470);
    const auto *msi0_1471 = buffer.data(msi0 + 1471);
    const auto *msi0_1472 = buffer.data(msi0 + 1472);
    const auto *msi0_1473 = buffer.data(msi0 + 1473);
    const auto *msi0_1474 = buffer.data(msi0 + 1474);
    const auto *msi0_1475 = buffer.data(msi0 + 1475);
    const auto *msi0_1476 = buffer.data(msi0 + 1476);
    const auto *msi0_1477 = buffer.data(msi0 + 1477);
    const auto *msi0_1478 = buffer.data(msi0 + 1478);
    const auto *msi0_1479 = buffer.data(msi0 + 1479);
    const auto *msi0_1480 = buffer.data(msi0 + 1480);
    const auto *msi0_1481 = buffer.data(msi0 + 1481);
    const auto *msi0_1482 = buffer.data(msi0 + 1482);
    const auto *msi0_1483 = buffer.data(msi0 + 1483);
    const auto *msi0_1485 = buffer.data(msi0 + 1485);
    const auto *msi0_1487 = buffer.data(msi0 + 1487);
    const auto *msi0_1488 = buffer.data(msi0 + 1488);
    const auto *msi0_1490 = buffer.data(msi0 + 1490);
    const auto *msi0_1491 = buffer.data(msi0 + 1491);
    const auto *msi0_1492 = buffer.data(msi0 + 1492);
    const auto *msi0_1494 = buffer.data(msi0 + 1494);

    const auto *msi1_1421 = buffer.data(msi1 + 1421);
    const auto *msi1_1423 = buffer.data(msi1 + 1423);
    const auto *msi1_1424 = buffer.data(msi1 + 1424);
    const auto *msi1_1425 = buffer.data(msi1 + 1425);
    const auto *msi1_1426 = buffer.data(msi1 + 1426);
    const auto *msi1_1427 = buffer.data(msi1 + 1427);
    const auto *msi1_1428 = buffer.data(msi1 + 1428);
    const auto *msi1_1429 = buffer.data(msi1 + 1429);
    const auto *msi1_1430 = buffer.data(msi1 + 1430);
    const auto *msi1_1431 = buffer.data(msi1 + 1431);
    const auto *msi1_1432 = buffer.data(msi1 + 1432);
    const auto *msi1_1433 = buffer.data(msi1 + 1433);
    const auto *msi1_1434 = buffer.data(msi1 + 1434);
    const auto *msi1_1435 = buffer.data(msi1 + 1435);
    const auto *msi1_1436 = buffer.data(msi1 + 1436);
    const auto *msi1_1437 = buffer.data(msi1 + 1437);
    const auto *msi1_1438 = buffer.data(msi1 + 1438);
    const auto *msi1_1439 = buffer.data(msi1 + 1439);
    const auto *msi1_1440 = buffer.data(msi1 + 1440);
    const auto *msi1_1441 = buffer.data(msi1 + 1441);
    const auto *msi1_1442 = buffer.data(msi1 + 1442);
    const auto *msi1_1443 = buffer.data(msi1 + 1443);
    const auto *msi1_1444 = buffer.data(msi1 + 1444);
    const auto *msi1_1445 = buffer.data(msi1 + 1445);
    const auto *msi1_1446 = buffer.data(msi1 + 1446);
    const auto *msi1_1447 = buffer.data(msi1 + 1447);
    const auto *msi1_1448 = buffer.data(msi1 + 1448);
    const auto *msi1_1449 = buffer.data(msi1 + 1449);
    const auto *msi1_1450 = buffer.data(msi1 + 1450);
    const auto *msi1_1451 = buffer.data(msi1 + 1451);
    const auto *msi1_1452 = buffer.data(msi1 + 1452);
    const auto *msi1_1453 = buffer.data(msi1 + 1453);
    const auto *msi1_1454 = buffer.data(msi1 + 1454);
    const auto *msi1_1455 = buffer.data(msi1 + 1455);
    const auto *msi1_1456 = buffer.data(msi1 + 1456);
    const auto *msi1_1457 = buffer.data(msi1 + 1457);
    const auto *msi1_1458 = buffer.data(msi1 + 1458);
    const auto *msi1_1459 = buffer.data(msi1 + 1459);
    const auto *msi1_1460 = buffer.data(msi1 + 1460);
    const auto *msi1_1461 = buffer.data(msi1 + 1461);
    const auto *msi1_1462 = buffer.data(msi1 + 1462);
    const auto *msi1_1463 = buffer.data(msi1 + 1463);
    const auto *msi1_1464 = buffer.data(msi1 + 1464);
    const auto *msi1_1465 = buffer.data(msi1 + 1465);
    const auto *msi1_1466 = buffer.data(msi1 + 1466);
    const auto *msi1_1467 = buffer.data(msi1 + 1467);
    const auto *msi1_1468 = buffer.data(msi1 + 1468);
    const auto *msi1_1469 = buffer.data(msi1 + 1469);
    const auto *msi1_1470 = buffer.data(msi1 + 1470);
    const auto *msi1_1471 = buffer.data(msi1 + 1471);
    const auto *msi1_1472 = buffer.data(msi1 + 1472);
    const auto *msi1_1473 = buffer.data(msi1 + 1473);
    const auto *msi1_1474 = buffer.data(msi1 + 1474);
    const auto *msi1_1475 = buffer.data(msi1 + 1475);
    const auto *msi1_1476 = buffer.data(msi1 + 1476);
    const auto *msi1_1477 = buffer.data(msi1 + 1477);
    const auto *msi1_1478 = buffer.data(msi1 + 1478);
    const auto *msi1_1479 = buffer.data(msi1 + 1479);
    const auto *msi1_1480 = buffer.data(msi1 + 1480);
    const auto *msi1_1481 = buffer.data(msi1 + 1481);
    const auto *msi1_1482 = buffer.data(msi1 + 1482);
    const auto *msi1_1483 = buffer.data(msi1 + 1483);
    const auto *msi1_1485 = buffer.data(msi1 + 1485);
    const auto *msi1_1487 = buffer.data(msi1 + 1487);
    const auto *msi1_1488 = buffer.data(msi1 + 1488);
    const auto *msi1_1490 = buffer.data(msi1 + 1490);
    const auto *msi1_1491 = buffer.data(msi1 + 1491);
    const auto *msi1_1492 = buffer.data(msi1 + 1492);
    const auto *msi1_1494 = buffer.data(msi1 + 1494);

    const auto *msk_1828 = buffer.data(msk + 1828);
    const auto *msk_1830 = buffer.data(msk + 1830);
    const auto *msk_1831 = buffer.data(msk + 1831);
    const auto *msk_1832 = buffer.data(msk + 1832);
    const auto *msk_1833 = buffer.data(msk + 1833);
    const auto *msk_1834 = buffer.data(msk + 1834);
    const auto *msk_1835 = buffer.data(msk + 1835);
    const auto *msk_1836 = buffer.data(msk + 1836);
    const auto *msk_1837 = buffer.data(msk + 1837);
    const auto *msk_1838 = buffer.data(msk + 1838);
    const auto *msk_1839 = buffer.data(msk + 1839);
    const auto *msk_1840 = buffer.data(msk + 1840);
    const auto *msk_1841 = buffer.data(msk + 1841);
    const auto *msk_1842 = buffer.data(msk + 1842);
    const auto *msk_1843 = buffer.data(msk + 1843);
    const auto *msk_1844 = buffer.data(msk + 1844);
    const auto *msk_1845 = buffer.data(msk + 1845);
    const auto *msk_1846 = buffer.data(msk + 1846);
    const auto *msk_1847 = buffer.data(msk + 1847);
    const auto *msk_1848 = buffer.data(msk + 1848);
    const auto *msk_1849 = buffer.data(msk + 1849);
    const auto *msk_1850 = buffer.data(msk + 1850);
    const auto *msk_1851 = buffer.data(msk + 1851);
    const auto *msk_1852 = buffer.data(msk + 1852);
    const auto *msk_1853 = buffer.data(msk + 1853);
    const auto *msk_1854 = buffer.data(msk + 1854);
    const auto *msk_1855 = buffer.data(msk + 1855);
    const auto *msk_1856 = buffer.data(msk + 1856);
    const auto *msk_1857 = buffer.data(msk + 1857);
    const auto *msk_1858 = buffer.data(msk + 1858);
    const auto *msk_1859 = buffer.data(msk + 1859);
    const auto *msk_1860 = buffer.data(msk + 1860);
    const auto *msk_1861 = buffer.data(msk + 1861);
    const auto *msk_1862 = buffer.data(msk + 1862);
    const auto *msk_1863 = buffer.data(msk + 1863);
    const auto *msk_1864 = buffer.data(msk + 1864);
    const auto *msk_1865 = buffer.data(msk + 1865);
    const auto *msk_1866 = buffer.data(msk + 1866);
    const auto *msk_1867 = buffer.data(msk + 1867);
    const auto *msk_1868 = buffer.data(msk + 1868);
    const auto *msk_1869 = buffer.data(msk + 1869);
    const auto *msk_1870 = buffer.data(msk + 1870);
    const auto *msk_1871 = buffer.data(msk + 1871);
    const auto *msk_1872 = buffer.data(msk + 1872);
    const auto *msk_1873 = buffer.data(msk + 1873);
    const auto *msk_1874 = buffer.data(msk + 1874);
    const auto *msk_1875 = buffer.data(msk + 1875);
    const auto *msk_1876 = buffer.data(msk + 1876);
    const auto *msk_1877 = buffer.data(msk + 1877);
    const auto *msk_1878 = buffer.data(msk + 1878);
    const auto *msk_1879 = buffer.data(msk + 1879);
    const auto *msk_1880 = buffer.data(msk + 1880);
    const auto *msk_1881 = buffer.data(msk + 1881);
    const auto *msk_1882 = buffer.data(msk + 1882);
    const auto *msk_1883 = buffer.data(msk + 1883);
    const auto *msk_1884 = buffer.data(msk + 1884);
    const auto *msk_1885 = buffer.data(msk + 1885);
    const auto *msk_1886 = buffer.data(msk + 1886);
    const auto *msk_1887 = buffer.data(msk + 1887);
    const auto *msk_1888 = buffer.data(msk + 1888);
    const auto *msk_1889 = buffer.data(msk + 1889);
    const auto *msk_1890 = buffer.data(msk + 1890);
    const auto *msk_1891 = buffer.data(msk + 1891);
    const auto *msk_1892 = buffer.data(msk + 1892);
    const auto *msk_1893 = buffer.data(msk + 1893);
    const auto *msk_1894 = buffer.data(msk + 1894);
    const auto *msk_1895 = buffer.data(msk + 1895);
    const auto *msk_1896 = buffer.data(msk + 1896);
    const auto *msk_1897 = buffer.data(msk + 1897);
    const auto *msk_1898 = buffer.data(msk + 1898);
    const auto *msk_1899 = buffer.data(msk + 1899);
    const auto *msk_1900 = buffer.data(msk + 1900);
    const auto *msk_1901 = buffer.data(msk + 1901);
    const auto *msk_1902 = buffer.data(msk + 1902);
    const auto *msk_1903 = buffer.data(msk + 1903);
    const auto *msk_1904 = buffer.data(msk + 1904);
    const auto *msk_1905 = buffer.data(msk + 1905);
    const auto *msk_1906 = buffer.data(msk + 1906);
    const auto *msk_1907 = buffer.data(msk + 1907);
    const auto *msk_1909 = buffer.data(msk + 1909);
    const auto *msk_1911 = buffer.data(msk + 1911);
    const auto *msk_1912 = buffer.data(msk + 1912);
    const auto *msk_1914 = buffer.data(msk + 1914);
    const auto *msk_1915 = buffer.data(msk + 1915);
    const auto *msk_1916 = buffer.data(msk + 1916);
    const auto *msk_1918 = buffer.data(msk + 1918);

#pragma omp simd aligned(t_2283, t_2284, t_2285, t_2286, t_2287, pc_x, pc_y, pc_z, lsk_1468, \
                         lsk_1504, msi0_1421, msi1_1421, msk_1828, msk_1833, msk_1834, \
                         msk_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2283[k] = f_3 * pc_x[k] * msk_1833[k];

        t_2284[k] = f_3 * pc_x[k] * msk_1834[k];

        t_2285[k] = f_3 * pc_x[k] * msk_1835[k];

        t_2286[k] = f_18 * lsk_1504[k]
                    + f_1 * msi0_1421[k]
                    - f_2 * msi1_1421[k]
                    + f_3 * pc_y[k] * msk_1828[k];

        t_2287[k] = f_19 * lsk_1468[k]
                    + f_3 * pc_z[k] * msk_1828[k];
    }

#pragma omp simd aligned(t_2288, t_2289, t_2290, pc_y, lsk_1506, lsk_1507, lsk_1508, \
                         msi0_1423, msi0_1424, msi0_1425, msi1_1423, msi1_1424, msi1_1425, \
                         msk_1830, msk_1831, msk_1832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2288[k] = f_18 * lsk_1506[k]
                    + f_12 * msi0_1423[k]
                    - f_13 * msi1_1423[k]
                    + f_3 * pc_y[k] * msk_1830[k];

        t_2289[k] = f_18 * lsk_1507[k]
                    + f_10 * msi0_1424[k]
                    - f_11 * msi1_1424[k]
                    + f_3 * pc_y[k] * msk_1831[k];

        t_2290[k] = f_18 * lsk_1508[k]
                    + f_8 * msi0_1425[k]
                    - f_9 * msi1_1425[k]
                    + f_3 * pc_y[k] * msk_1832[k];
    }

#pragma omp simd aligned(t_2291, t_2292, t_2293, pc_y, lsk_1509, lsk_1510, lsk_1511, \
                         msi0_1426, msi0_1427, msi1_1426, msi1_1427, msk_1833, msk_1834, \
                         msk_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2291[k] = f_18 * lsk_1509[k]
                    + f_6 * msi0_1426[k]
                    - f_7 * msi1_1426[k]
                    + f_3 * pc_y[k] * msk_1833[k];

        t_2292[k] = f_18 * lsk_1510[k]
                    + f_4 * msi0_1427[k]
                    - f_5 * msi1_1427[k]
                    + f_3 * pc_y[k] * msk_1834[k];

        t_2293[k] = f_18 * lsk_1511[k]
                    + f_3 * pc_y[k] * msk_1835[k];
    }

#pragma omp simd aligned(t_2294, t_2295, t_2296, pc_x, pc_z, lsk_1475, msi0_1427, msi0_1428, \
                         msi0_1429, msi1_1427, msi1_1428, msi1_1429, msk_1835, msk_1836, \
                         msk_1837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2294[k] = f_19 * lsk_1475[k]
                    + f_1 * msi0_1427[k]
                    - f_2 * msi1_1427[k]
                    + f_3 * pc_z[k] * msk_1835[k];

        t_2295[k] = f_1 * msi0_1428[k]
                    - f_2 * msi1_1428[k]
                    + f_3 * pc_x[k] * msk_1836[k];

        t_2296[k] = f_22 * msi0_1429[k]
                    - f_23 * msi1_1429[k]
                    + f_3 * pc_x[k] * msk_1837[k];
    }

#pragma omp simd aligned(t_2297, t_2298, t_2299, pc_x, msi0_1430, msi0_1431, msi0_1432, \
                         msi1_1430, msi1_1431, msi1_1432, msk_1838, msk_1839, \
                         msk_1840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2297[k] = f_22 * msi0_1430[k]
                    - f_23 * msi1_1430[k]
                    + f_3 * pc_x[k] * msk_1838[k];

        t_2298[k] = f_12 * msi0_1431[k]
                    - f_13 * msi1_1431[k]
                    + f_3 * pc_x[k] * msk_1839[k];

        t_2299[k] = f_12 * msi0_1432[k]
                    - f_13 * msi1_1432[k]
                    + f_3 * pc_x[k] * msk_1840[k];
    }

#pragma omp simd aligned(t_2300, t_2301, t_2302, pc_x, msi0_1433, msi0_1434, msi0_1435, \
                         msi1_1433, msi1_1434, msi1_1435, msk_1841, msk_1842, \
                         msk_1843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2300[k] = f_12 * msi0_1433[k]
                    - f_13 * msi1_1433[k]
                    + f_3 * pc_x[k] * msk_1841[k];

        t_2301[k] = f_10 * msi0_1434[k]
                    - f_11 * msi1_1434[k]
                    + f_3 * pc_x[k] * msk_1842[k];

        t_2302[k] = f_10 * msi0_1435[k]
                    - f_11 * msi1_1435[k]
                    + f_3 * pc_x[k] * msk_1843[k];
    }

#pragma omp simd aligned(t_2303, t_2304, t_2305, pc_x, msi0_1436, msi0_1437, msi0_1438, \
                         msi1_1436, msi1_1437, msi1_1438, msk_1844, msk_1845, \
                         msk_1846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2303[k] = f_10 * msi0_1436[k]
                    - f_11 * msi1_1436[k]
                    + f_3 * pc_x[k] * msk_1844[k];

        t_2304[k] = f_10 * msi0_1437[k]
                    - f_11 * msi1_1437[k]
                    + f_3 * pc_x[k] * msk_1845[k];

        t_2305[k] = f_8 * msi0_1438[k]
                    - f_9 * msi1_1438[k]
                    + f_3 * pc_x[k] * msk_1846[k];
    }

#pragma omp simd aligned(t_2306, t_2307, t_2308, pc_x, msi0_1439, msi0_1440, msi0_1441, \
                         msi1_1439, msi1_1440, msi1_1441, msk_1847, msk_1848, \
                         msk_1849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2306[k] = f_8 * msi0_1439[k]
                    - f_9 * msi1_1439[k]
                    + f_3 * pc_x[k] * msk_1847[k];

        t_2307[k] = f_8 * msi0_1440[k]
                    - f_9 * msi1_1440[k]
                    + f_3 * pc_x[k] * msk_1848[k];

        t_2308[k] = f_8 * msi0_1441[k]
                    - f_9 * msi1_1441[k]
                    + f_3 * pc_x[k] * msk_1849[k];
    }

#pragma omp simd aligned(t_2309, t_2310, t_2311, pc_x, msi0_1442, msi0_1443, msi0_1444, \
                         msi1_1442, msi1_1443, msi1_1444, msk_1850, msk_1851, \
                         msk_1852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2309[k] = f_8 * msi0_1442[k]
                    - f_9 * msi1_1442[k]
                    + f_3 * pc_x[k] * msk_1850[k];

        t_2310[k] = f_6 * msi0_1443[k]
                    - f_7 * msi1_1443[k]
                    + f_3 * pc_x[k] * msk_1851[k];

        t_2311[k] = f_6 * msi0_1444[k]
                    - f_7 * msi1_1444[k]
                    + f_3 * pc_x[k] * msk_1852[k];
    }

#pragma omp simd aligned(t_2312, t_2313, t_2314, pc_x, msi0_1445, msi0_1446, msi0_1447, \
                         msi1_1445, msi1_1446, msi1_1447, msk_1853, msk_1854, \
                         msk_1855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2312[k] = f_6 * msi0_1445[k]
                    - f_7 * msi1_1445[k]
                    + f_3 * pc_x[k] * msk_1853[k];

        t_2313[k] = f_6 * msi0_1446[k]
                    - f_7 * msi1_1446[k]
                    + f_3 * pc_x[k] * msk_1854[k];

        t_2314[k] = f_6 * msi0_1447[k]
                    - f_7 * msi1_1447[k]
                    + f_3 * pc_x[k] * msk_1855[k];
    }

#pragma omp simd aligned(t_2315, t_2316, t_2317, pc_x, msi0_1448, msi0_1449, msi0_1450, \
                         msi1_1448, msi1_1449, msi1_1450, msk_1856, msk_1857, \
                         msk_1858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2315[k] = f_6 * msi0_1448[k]
                    - f_7 * msi1_1448[k]
                    + f_3 * pc_x[k] * msk_1856[k];

        t_2316[k] = f_4 * msi0_1449[k]
                    - f_5 * msi1_1449[k]
                    + f_3 * pc_x[k] * msk_1857[k];

        t_2317[k] = f_4 * msi0_1450[k]
                    - f_5 * msi1_1450[k]
                    + f_3 * pc_x[k] * msk_1858[k];
    }

#pragma omp simd aligned(t_2318, t_2319, t_2320, pc_x, msi0_1451, msi0_1452, msi0_1453, \
                         msi1_1451, msi1_1452, msi1_1453, msk_1859, msk_1860, \
                         msk_1861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2318[k] = f_4 * msi0_1451[k]
                    - f_5 * msi1_1451[k]
                    + f_3 * pc_x[k] * msk_1859[k];

        t_2319[k] = f_4 * msi0_1452[k]
                    - f_5 * msi1_1452[k]
                    + f_3 * pc_x[k] * msk_1860[k];

        t_2320[k] = f_4 * msi0_1453[k]
                    - f_5 * msi1_1453[k]
                    + f_3 * pc_x[k] * msk_1861[k];
    }

#pragma omp simd aligned(t_2321, t_2322, t_2323, t_2324, t_2325, pc_x, msi0_1454, msi0_1455, \
                         msi1_1454, msi1_1455, msk_1862, msk_1863, msk_1864, msk_1865, \
                         msk_1866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2321[k] = f_4 * msi0_1454[k]
                    - f_5 * msi1_1454[k]
                    + f_3 * pc_x[k] * msk_1862[k];

        t_2322[k] = f_4 * msi0_1455[k]
                    - f_5 * msi1_1455[k]
                    + f_3 * pc_x[k] * msk_1863[k];

        t_2323[k] = f_3 * pc_x[k] * msk_1864[k];

        t_2324[k] = f_3 * pc_x[k] * msk_1865[k];

        t_2325[k] = f_3 * pc_x[k] * msk_1866[k];
    }

#pragma omp simd aligned(t_2326, t_2327, t_2328, t_2329, t_2330, pc_x, msk_1867, msk_1868, \
                         msk_1869, msk_1870, msk_1871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2326[k] = f_3 * pc_x[k] * msk_1867[k];

        t_2327[k] = f_3 * pc_x[k] * msk_1868[k];

        t_2328[k] = f_3 * pc_x[k] * msk_1869[k];

        t_2329[k] = f_3 * pc_x[k] * msk_1870[k];

        t_2330[k] = f_3 * pc_x[k] * msk_1871[k];
    }

#pragma omp simd aligned(t_2331, t_2332, t_2333, pc_y, pc_z, lsk_1504, lsk_1540, lsk_1542, \
                         msi0_1449, msi0_1451, msi1_1449, msi1_1451, msk_1864, \
                         msk_1866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2331[k] = f_17 * lsk_1540[k]
                    + f_1 * msi0_1449[k]
                    - f_2 * msi1_1449[k]
                    + f_3 * pc_y[k] * msk_1864[k];

        t_2332[k] = f_20 * lsk_1504[k]
                    + f_3 * pc_z[k] * msk_1864[k];

        t_2333[k] = f_17 * lsk_1542[k]
                    + f_12 * msi0_1451[k]
                    - f_13 * msi1_1451[k]
                    + f_3 * pc_y[k] * msk_1866[k];
    }

#pragma omp simd aligned(t_2334, t_2335, t_2336, pc_y, lsk_1543, lsk_1544, lsk_1545, \
                         msi0_1452, msi0_1453, msi0_1454, msi1_1452, msi1_1453, msi1_1454, \
                         msk_1867, msk_1868, msk_1869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2334[k] = f_17 * lsk_1543[k]
                    + f_10 * msi0_1452[k]
                    - f_11 * msi1_1452[k]
                    + f_3 * pc_y[k] * msk_1867[k];

        t_2335[k] = f_17 * lsk_1544[k]
                    + f_8 * msi0_1453[k]
                    - f_9 * msi1_1453[k]
                    + f_3 * pc_y[k] * msk_1868[k];

        t_2336[k] = f_17 * lsk_1545[k]
                    + f_6 * msi0_1454[k]
                    - f_7 * msi1_1454[k]
                    + f_3 * pc_y[k] * msk_1869[k];
    }

#pragma omp simd aligned(t_2337, t_2338, t_2339, pc_y, pc_z, lsk_1511, lsk_1546, lsk_1547, \
                         msi0_1455, msi1_1455, msk_1870, msk_1871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2337[k] = f_17 * lsk_1546[k]
                    + f_4 * msi0_1455[k]
                    - f_5 * msi1_1455[k]
                    + f_3 * pc_y[k] * msk_1870[k];

        t_2338[k] = f_17 * lsk_1547[k]
                    + f_3 * pc_y[k] * msk_1871[k];

        t_2339[k] = f_20 * lsk_1511[k]
                    + f_1 * msi0_1455[k]
                    - f_2 * msi1_1455[k]
                    + f_3 * pc_z[k] * msk_1871[k];
    }

#pragma omp simd aligned(t_2340, t_2341, t_2342, pc_x, msi0_1456, msi0_1457, msi0_1458, \
                         msi1_1456, msi1_1457, msi1_1458, msk_1872, msk_1873, \
                         msk_1874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2340[k] = f_1 * msi0_1456[k]
                    - f_2 * msi1_1456[k]
                    + f_3 * pc_x[k] * msk_1872[k];

        t_2341[k] = f_22 * msi0_1457[k]
                    - f_23 * msi1_1457[k]
                    + f_3 * pc_x[k] * msk_1873[k];

        t_2342[k] = f_22 * msi0_1458[k]
                    - f_23 * msi1_1458[k]
                    + f_3 * pc_x[k] * msk_1874[k];
    }

#pragma omp simd aligned(t_2343, t_2344, t_2345, pc_x, msi0_1459, msi0_1460, msi0_1461, \
                         msi1_1459, msi1_1460, msi1_1461, msk_1875, msk_1876, \
                         msk_1877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2343[k] = f_12 * msi0_1459[k]
                    - f_13 * msi1_1459[k]
                    + f_3 * pc_x[k] * msk_1875[k];

        t_2344[k] = f_12 * msi0_1460[k]
                    - f_13 * msi1_1460[k]
                    + f_3 * pc_x[k] * msk_1876[k];

        t_2345[k] = f_12 * msi0_1461[k]
                    - f_13 * msi1_1461[k]
                    + f_3 * pc_x[k] * msk_1877[k];
    }

#pragma omp simd aligned(t_2346, t_2347, t_2348, pc_x, msi0_1462, msi0_1463, msi0_1464, \
                         msi1_1462, msi1_1463, msi1_1464, msk_1878, msk_1879, \
                         msk_1880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2346[k] = f_10 * msi0_1462[k]
                    - f_11 * msi1_1462[k]
                    + f_3 * pc_x[k] * msk_1878[k];

        t_2347[k] = f_10 * msi0_1463[k]
                    - f_11 * msi1_1463[k]
                    + f_3 * pc_x[k] * msk_1879[k];

        t_2348[k] = f_10 * msi0_1464[k]
                    - f_11 * msi1_1464[k]
                    + f_3 * pc_x[k] * msk_1880[k];
    }

#pragma omp simd aligned(t_2349, t_2350, t_2351, pc_x, msi0_1465, msi0_1466, msi0_1467, \
                         msi1_1465, msi1_1466, msi1_1467, msk_1881, msk_1882, \
                         msk_1883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2349[k] = f_10 * msi0_1465[k]
                    - f_11 * msi1_1465[k]
                    + f_3 * pc_x[k] * msk_1881[k];

        t_2350[k] = f_8 * msi0_1466[k]
                    - f_9 * msi1_1466[k]
                    + f_3 * pc_x[k] * msk_1882[k];

        t_2351[k] = f_8 * msi0_1467[k]
                    - f_9 * msi1_1467[k]
                    + f_3 * pc_x[k] * msk_1883[k];
    }

#pragma omp simd aligned(t_2352, t_2353, t_2354, pc_x, msi0_1468, msi0_1469, msi0_1470, \
                         msi1_1468, msi1_1469, msi1_1470, msk_1884, msk_1885, \
                         msk_1886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2352[k] = f_8 * msi0_1468[k]
                    - f_9 * msi1_1468[k]
                    + f_3 * pc_x[k] * msk_1884[k];

        t_2353[k] = f_8 * msi0_1469[k]
                    - f_9 * msi1_1469[k]
                    + f_3 * pc_x[k] * msk_1885[k];

        t_2354[k] = f_8 * msi0_1470[k]
                    - f_9 * msi1_1470[k]
                    + f_3 * pc_x[k] * msk_1886[k];
    }

#pragma omp simd aligned(t_2355, t_2356, t_2357, pc_x, msi0_1471, msi0_1472, msi0_1473, \
                         msi1_1471, msi1_1472, msi1_1473, msk_1887, msk_1888, \
                         msk_1889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2355[k] = f_6 * msi0_1471[k]
                    - f_7 * msi1_1471[k]
                    + f_3 * pc_x[k] * msk_1887[k];

        t_2356[k] = f_6 * msi0_1472[k]
                    - f_7 * msi1_1472[k]
                    + f_3 * pc_x[k] * msk_1888[k];

        t_2357[k] = f_6 * msi0_1473[k]
                    - f_7 * msi1_1473[k]
                    + f_3 * pc_x[k] * msk_1889[k];
    }

#pragma omp simd aligned(t_2358, t_2359, t_2360, pc_x, msi0_1474, msi0_1475, msi0_1476, \
                         msi1_1474, msi1_1475, msi1_1476, msk_1890, msk_1891, \
                         msk_1892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2358[k] = f_6 * msi0_1474[k]
                    - f_7 * msi1_1474[k]
                    + f_3 * pc_x[k] * msk_1890[k];

        t_2359[k] = f_6 * msi0_1475[k]
                    - f_7 * msi1_1475[k]
                    + f_3 * pc_x[k] * msk_1891[k];

        t_2360[k] = f_6 * msi0_1476[k]
                    - f_7 * msi1_1476[k]
                    + f_3 * pc_x[k] * msk_1892[k];
    }

#pragma omp simd aligned(t_2361, t_2362, t_2363, pc_x, msi0_1477, msi0_1478, msi0_1479, \
                         msi1_1477, msi1_1478, msi1_1479, msk_1893, msk_1894, \
                         msk_1895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2361[k] = f_4 * msi0_1477[k]
                    - f_5 * msi1_1477[k]
                    + f_3 * pc_x[k] * msk_1893[k];

        t_2362[k] = f_4 * msi0_1478[k]
                    - f_5 * msi1_1478[k]
                    + f_3 * pc_x[k] * msk_1894[k];

        t_2363[k] = f_4 * msi0_1479[k]
                    - f_5 * msi1_1479[k]
                    + f_3 * pc_x[k] * msk_1895[k];
    }

#pragma omp simd aligned(t_2364, t_2365, t_2366, pc_x, msi0_1480, msi0_1481, msi0_1482, \
                         msi1_1480, msi1_1481, msi1_1482, msk_1896, msk_1897, \
                         msk_1898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2364[k] = f_4 * msi0_1480[k]
                    - f_5 * msi1_1480[k]
                    + f_3 * pc_x[k] * msk_1896[k];

        t_2365[k] = f_4 * msi0_1481[k]
                    - f_5 * msi1_1481[k]
                    + f_3 * pc_x[k] * msk_1897[k];

        t_2366[k] = f_4 * msi0_1482[k]
                    - f_5 * msi1_1482[k]
                    + f_3 * pc_x[k] * msk_1898[k];
    }

#pragma omp simd aligned(t_2367, t_2368, t_2369, t_2370, t_2371, t_2372, pc_x, msi0_1483, \
                         msi1_1483, msk_1899, msk_1900, msk_1901, msk_1902, msk_1903, \
                         msk_1904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2367[k] = f_4 * msi0_1483[k]
                    - f_5 * msi1_1483[k]
                    + f_3 * pc_x[k] * msk_1899[k];

        t_2368[k] = f_3 * pc_x[k] * msk_1900[k];

        t_2369[k] = f_3 * pc_x[k] * msk_1901[k];

        t_2370[k] = f_3 * pc_x[k] * msk_1902[k];

        t_2371[k] = f_3 * pc_x[k] * msk_1903[k];

        t_2372[k] = f_3 * pc_x[k] * msk_1904[k];
    }

#pragma omp simd aligned(t_2373, t_2374, t_2375, t_2376, t_2377, pc_x, pc_y, pc_z, lsk_1540, \
                         lsk_1576, msi0_1477, msi1_1477, msk_1900, msk_1905, msk_1906, \
                         msk_1907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2373[k] = f_3 * pc_x[k] * msk_1905[k];

        t_2374[k] = f_3 * pc_x[k] * msk_1906[k];

        t_2375[k] = f_3 * pc_x[k] * msk_1907[k];

        t_2376[k] = f_16 * lsk_1576[k]
                    + f_1 * msi0_1477[k]
                    - f_2 * msi1_1477[k]
                    + f_3 * pc_y[k] * msk_1900[k];

        t_2377[k] = f_24 * lsk_1540[k]
                    + f_3 * pc_z[k] * msk_1900[k];
    }

#pragma omp simd aligned(t_2378, t_2379, t_2380, pc_y, lsk_1578, lsk_1579, lsk_1580, \
                         msi0_1479, msi0_1480, msi0_1481, msi1_1479, msi1_1480, msi1_1481, \
                         msk_1902, msk_1903, msk_1904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2378[k] = f_16 * lsk_1578[k]
                    + f_12 * msi0_1479[k]
                    - f_13 * msi1_1479[k]
                    + f_3 * pc_y[k] * msk_1902[k];

        t_2379[k] = f_16 * lsk_1579[k]
                    + f_10 * msi0_1480[k]
                    - f_11 * msi1_1480[k]
                    + f_3 * pc_y[k] * msk_1903[k];

        t_2380[k] = f_16 * lsk_1580[k]
                    + f_8 * msi0_1481[k]
                    - f_9 * msi1_1481[k]
                    + f_3 * pc_y[k] * msk_1904[k];
    }

#pragma omp simd aligned(t_2381, t_2382, t_2383, pc_y, lsk_1581, lsk_1582, lsk_1583, \
                         msi0_1482, msi0_1483, msi1_1482, msi1_1483, msk_1905, msk_1906, \
                         msk_1907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2381[k] = f_16 * lsk_1581[k]
                    + f_6 * msi0_1482[k]
                    - f_7 * msi1_1482[k]
                    + f_3 * pc_y[k] * msk_1905[k];

        t_2382[k] = f_16 * lsk_1582[k]
                    + f_4 * msi0_1483[k]
                    - f_5 * msi1_1483[k]
                    + f_3 * pc_y[k] * msk_1906[k];

        t_2383[k] = f_16 * lsk_1583[k]
                    + f_3 * pc_y[k] * msk_1907[k];
    }

#pragma omp simd aligned(t_2384, t_2385, t_2386, pa_y, pc_x, pc_y, pc_z, lsl0_1980, lsk_1547, \
                         lsl1_1980, msi0_1483, msi0_1485, msi1_1483, msi1_1485, msk_1907, \
                         msk_1909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2384[k] = f_24 * lsk_1547[k]
                    + f_1 * msi0_1483[k]
                    - f_2 * msi1_1483[k]
                    + f_3 * pc_z[k] * msk_1907[k];

        t_2385[k] = pa_y[k] * lsl0_1980[k]
                    - f_14 * pc_y[k] * lsl1_1980[k];

        t_2386[k] = f_22 * msi0_1485[k]
                    - f_23 * msi1_1485[k]
                    + f_3 * pc_x[k] * msk_1909[k];
    }

#pragma omp simd aligned(t_2387, t_2388, t_2389, pa_y, pc_x, pc_y, lsl0_1982, lsl1_1982, \
                         msi0_1487, msi0_1488, msi1_1487, msi1_1488, msk_1911, \
                         msk_1912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2387[k] = pa_y[k] * lsl0_1982[k]
                    - f_14 * pc_y[k] * lsl1_1982[k];

        t_2388[k] = f_12 * msi0_1487[k]
                    - f_13 * msi1_1487[k]
                    + f_3 * pc_x[k] * msk_1911[k];

        t_2389[k] = f_12 * msi0_1488[k]
                    - f_13 * msi1_1488[k]
                    + f_3 * pc_x[k] * msk_1912[k];
    }

#pragma omp simd aligned(t_2390, t_2391, t_2392, pa_y, pc_x, pc_y, lsl0_1985, lsl1_1985, \
                         msi0_1490, msi0_1491, msi1_1490, msi1_1491, msk_1914, \
                         msk_1915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2390[k] = pa_y[k] * lsl0_1985[k]
                    - f_14 * pc_y[k] * lsl1_1985[k];

        t_2391[k] = f_10 * msi0_1490[k]
                    - f_11 * msi1_1490[k]
                    + f_3 * pc_x[k] * msk_1914[k];

        t_2392[k] = f_10 * msi0_1491[k]
                    - f_11 * msi1_1491[k]
                    + f_3 * pc_x[k] * msk_1915[k];
    }

#pragma omp simd aligned(t_2393, t_2394, t_2395, pa_y, pc_x, pc_y, lsl0_1989, lsl1_1989, \
                         msi0_1492, msi0_1494, msi1_1492, msi1_1494, msk_1916, \
                         msk_1918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2393[k] = f_10 * msi0_1492[k]
                    - f_11 * msi1_1492[k]
                    + f_3 * pc_x[k] * msk_1916[k];

        t_2394[k] = pa_y[k] * lsl0_1989[k]
                    - f_14 * pc_y[k] * lsl1_1989[k];

        t_2395[k] = f_8 * msi0_1494[k]
                    - f_9 * msi1_1494[k]
                    + f_3 * pc_x[k] * msk_1918[k];
    }
}

static auto
compute_prim_msl_three_center_electron_repulsion_0_piece21(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsl0,
                                                           const size_t lsk, const size_t lsl1,
                                                           const size_t msi0, const size_t msi1,
                                                           const size_t msk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsl0_1994 = buffer.data(lsl0 + 1994);
    const auto *lsl0_2000 = buffer.data(lsl0 + 2000);
    const auto *lsl0_2007 = buffer.data(lsl0 + 2007);
    const auto *lsl0_2016 = buffer.data(lsl0 + 2016);
    const auto *lsl0_2018 = buffer.data(lsl0 + 2018);
    const auto *lsl0_2019 = buffer.data(lsl0 + 2019);
    const auto *lsl0_2020 = buffer.data(lsl0 + 2020);
    const auto *lsl0_2021 = buffer.data(lsl0 + 2021);
    const auto *lsl0_2022 = buffer.data(lsl0 + 2022);
    const auto *lsl0_2024 = buffer.data(lsl0 + 2024);

    const auto *lsk_1576 = buffer.data(lsk + 1576);
    const auto *lsk_1612 = buffer.data(lsk + 1612);
    const auto *lsk_1614 = buffer.data(lsk + 1614);
    const auto *lsk_1615 = buffer.data(lsk + 1615);
    const auto *lsk_1616 = buffer.data(lsk + 1616);
    const auto *lsk_1617 = buffer.data(lsk + 1617);
    const auto *lsk_1618 = buffer.data(lsk + 1618);
    const auto *lsk_1619 = buffer.data(lsk + 1619);

    const auto *lsl1_1994 = buffer.data(lsl1 + 1994);
    const auto *lsl1_2000 = buffer.data(lsl1 + 2000);
    const auto *lsl1_2007 = buffer.data(lsl1 + 2007);
    const auto *lsl1_2016 = buffer.data(lsl1 + 2016);
    const auto *lsl1_2018 = buffer.data(lsl1 + 2018);
    const auto *lsl1_2019 = buffer.data(lsl1 + 2019);
    const auto *lsl1_2020 = buffer.data(lsl1 + 2020);
    const auto *lsl1_2021 = buffer.data(lsl1 + 2021);
    const auto *lsl1_2022 = buffer.data(lsl1 + 2022);
    const auto *lsl1_2024 = buffer.data(lsl1 + 2024);

    const auto *msi0_1495 = buffer.data(msi0 + 1495);
    const auto *msi0_1496 = buffer.data(msi0 + 1496);
    const auto *msi0_1497 = buffer.data(msi0 + 1497);
    const auto *msi0_1499 = buffer.data(msi0 + 1499);
    const auto *msi0_1500 = buffer.data(msi0 + 1500);
    const auto *msi0_1501 = buffer.data(msi0 + 1501);
    const auto *msi0_1502 = buffer.data(msi0 + 1502);
    const auto *msi0_1503 = buffer.data(msi0 + 1503);
    const auto *msi0_1505 = buffer.data(msi0 + 1505);
    const auto *msi0_1506 = buffer.data(msi0 + 1506);
    const auto *msi0_1507 = buffer.data(msi0 + 1507);
    const auto *msi0_1508 = buffer.data(msi0 + 1508);
    const auto *msi0_1509 = buffer.data(msi0 + 1509);
    const auto *msi0_1510 = buffer.data(msi0 + 1510);
    const auto *msi0_1512 = buffer.data(msi0 + 1512);
    const auto *msi0_1514 = buffer.data(msi0 + 1514);
    const auto *msi0_1515 = buffer.data(msi0 + 1515);
    const auto *msi0_1517 = buffer.data(msi0 + 1517);
    const auto *msi0_1518 = buffer.data(msi0 + 1518);
    const auto *msi0_1519 = buffer.data(msi0 + 1519);
    const auto *msi0_1521 = buffer.data(msi0 + 1521);
    const auto *msi0_1522 = buffer.data(msi0 + 1522);
    const auto *msi0_1523 = buffer.data(msi0 + 1523);
    const auto *msi0_1524 = buffer.data(msi0 + 1524);
    const auto *msi0_1526 = buffer.data(msi0 + 1526);
    const auto *msi0_1527 = buffer.data(msi0 + 1527);
    const auto *msi0_1528 = buffer.data(msi0 + 1528);
    const auto *msi0_1529 = buffer.data(msi0 + 1529);
    const auto *msi0_1530 = buffer.data(msi0 + 1530);
    const auto *msi0_1532 = buffer.data(msi0 + 1532);
    const auto *msi0_1533 = buffer.data(msi0 + 1533);
    const auto *msi0_1534 = buffer.data(msi0 + 1534);
    const auto *msi0_1535 = buffer.data(msi0 + 1535);
    const auto *msi0_1536 = buffer.data(msi0 + 1536);
    const auto *msi0_1537 = buffer.data(msi0 + 1537);
    const auto *msi0_1538 = buffer.data(msi0 + 1538);
    const auto *msi0_1539 = buffer.data(msi0 + 1539);

    const auto *msi1_1495 = buffer.data(msi1 + 1495);
    const auto *msi1_1496 = buffer.data(msi1 + 1496);
    const auto *msi1_1497 = buffer.data(msi1 + 1497);
    const auto *msi1_1499 = buffer.data(msi1 + 1499);
    const auto *msi1_1500 = buffer.data(msi1 + 1500);
    const auto *msi1_1501 = buffer.data(msi1 + 1501);
    const auto *msi1_1502 = buffer.data(msi1 + 1502);
    const auto *msi1_1503 = buffer.data(msi1 + 1503);
    const auto *msi1_1505 = buffer.data(msi1 + 1505);
    const auto *msi1_1506 = buffer.data(msi1 + 1506);
    const auto *msi1_1507 = buffer.data(msi1 + 1507);
    const auto *msi1_1508 = buffer.data(msi1 + 1508);
    const auto *msi1_1509 = buffer.data(msi1 + 1509);
    const auto *msi1_1510 = buffer.data(msi1 + 1510);
    const auto *msi1_1512 = buffer.data(msi1 + 1512);
    const auto *msi1_1514 = buffer.data(msi1 + 1514);
    const auto *msi1_1515 = buffer.data(msi1 + 1515);
    const auto *msi1_1517 = buffer.data(msi1 + 1517);
    const auto *msi1_1518 = buffer.data(msi1 + 1518);
    const auto *msi1_1519 = buffer.data(msi1 + 1519);
    const auto *msi1_1521 = buffer.data(msi1 + 1521);
    const auto *msi1_1522 = buffer.data(msi1 + 1522);
    const auto *msi1_1523 = buffer.data(msi1 + 1523);
    const auto *msi1_1524 = buffer.data(msi1 + 1524);
    const auto *msi1_1526 = buffer.data(msi1 + 1526);
    const auto *msi1_1527 = buffer.data(msi1 + 1527);
    const auto *msi1_1528 = buffer.data(msi1 + 1528);
    const auto *msi1_1529 = buffer.data(msi1 + 1529);
    const auto *msi1_1530 = buffer.data(msi1 + 1530);
    const auto *msi1_1532 = buffer.data(msi1 + 1532);
    const auto *msi1_1533 = buffer.data(msi1 + 1533);
    const auto *msi1_1534 = buffer.data(msi1 + 1534);
    const auto *msi1_1535 = buffer.data(msi1 + 1535);
    const auto *msi1_1536 = buffer.data(msi1 + 1536);
    const auto *msi1_1537 = buffer.data(msi1 + 1537);
    const auto *msi1_1538 = buffer.data(msi1 + 1538);
    const auto *msi1_1539 = buffer.data(msi1 + 1539);

    const auto *msk_1919 = buffer.data(msk + 1919);
    const auto *msk_1920 = buffer.data(msk + 1920);
    const auto *msk_1921 = buffer.data(msk + 1921);
    const auto *msk_1923 = buffer.data(msk + 1923);
    const auto *msk_1924 = buffer.data(msk + 1924);
    const auto *msk_1925 = buffer.data(msk + 1925);
    const auto *msk_1926 = buffer.data(msk + 1926);
    const auto *msk_1927 = buffer.data(msk + 1927);
    const auto *msk_1929 = buffer.data(msk + 1929);
    const auto *msk_1930 = buffer.data(msk + 1930);
    const auto *msk_1931 = buffer.data(msk + 1931);
    const auto *msk_1932 = buffer.data(msk + 1932);
    const auto *msk_1933 = buffer.data(msk + 1933);
    const auto *msk_1934 = buffer.data(msk + 1934);
    const auto *msk_1936 = buffer.data(msk + 1936);
    const auto *msk_1937 = buffer.data(msk + 1937);
    const auto *msk_1938 = buffer.data(msk + 1938);
    const auto *msk_1939 = buffer.data(msk + 1939);
    const auto *msk_1940 = buffer.data(msk + 1940);
    const auto *msk_1941 = buffer.data(msk + 1941);
    const auto *msk_1942 = buffer.data(msk + 1942);
    const auto *msk_1943 = buffer.data(msk + 1943);
    const auto *msk_1944 = buffer.data(msk + 1944);
    const auto *msk_1946 = buffer.data(msk + 1946);
    const auto *msk_1947 = buffer.data(msk + 1947);
    const auto *msk_1949 = buffer.data(msk + 1949);
    const auto *msk_1950 = buffer.data(msk + 1950);
    const auto *msk_1951 = buffer.data(msk + 1951);
    const auto *msk_1953 = buffer.data(msk + 1953);
    const auto *msk_1954 = buffer.data(msk + 1954);
    const auto *msk_1955 = buffer.data(msk + 1955);
    const auto *msk_1956 = buffer.data(msk + 1956);
    const auto *msk_1958 = buffer.data(msk + 1958);
    const auto *msk_1959 = buffer.data(msk + 1959);
    const auto *msk_1960 = buffer.data(msk + 1960);
    const auto *msk_1961 = buffer.data(msk + 1961);
    const auto *msk_1962 = buffer.data(msk + 1962);
    const auto *msk_1964 = buffer.data(msk + 1964);
    const auto *msk_1965 = buffer.data(msk + 1965);
    const auto *msk_1966 = buffer.data(msk + 1966);
    const auto *msk_1967 = buffer.data(msk + 1967);
    const auto *msk_1968 = buffer.data(msk + 1968);
    const auto *msk_1969 = buffer.data(msk + 1969);
    const auto *msk_1971 = buffer.data(msk + 1971);
    const auto *msk_1972 = buffer.data(msk + 1972);
    const auto *msk_1973 = buffer.data(msk + 1973);
    const auto *msk_1974 = buffer.data(msk + 1974);
    const auto *msk_1975 = buffer.data(msk + 1975);
    const auto *msk_1976 = buffer.data(msk + 1976);
    const auto *msk_1977 = buffer.data(msk + 1977);
    const auto *msk_1978 = buffer.data(msk + 1978);
    const auto *msk_1979 = buffer.data(msk + 1979);

#pragma omp simd aligned(t_2396, t_2397, t_2398, pc_x, msi0_1495, msi0_1496, msi0_1497, \
                         msi1_1495, msi1_1496, msi1_1497, msk_1919, msk_1920, \
                         msk_1921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2396[k] = f_8 * msi0_1495[k]
                    - f_9 * msi1_1495[k]
                    + f_3 * pc_x[k] * msk_1919[k];

        t_2397[k] = f_8 * msi0_1496[k]
                    - f_9 * msi1_1496[k]
                    + f_3 * pc_x[k] * msk_1920[k];

        t_2398[k] = f_8 * msi0_1497[k]
                    - f_9 * msi1_1497[k]
                    + f_3 * pc_x[k] * msk_1921[k];
    }

#pragma omp simd aligned(t_2399, t_2400, t_2401, pa_y, pc_x, pc_y, lsl0_1994, lsl1_1994, \
                         msi0_1499, msi0_1500, msi1_1499, msi1_1500, msk_1923, \
                         msk_1924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2399[k] = pa_y[k] * lsl0_1994[k]
                    - f_14 * pc_y[k] * lsl1_1994[k];

        t_2400[k] = f_6 * msi0_1499[k]
                    - f_7 * msi1_1499[k]
                    + f_3 * pc_x[k] * msk_1923[k];

        t_2401[k] = f_6 * msi0_1500[k]
                    - f_7 * msi1_1500[k]
                    + f_3 * pc_x[k] * msk_1924[k];
    }

#pragma omp simd aligned(t_2402, t_2403, t_2404, pc_x, msi0_1501, msi0_1502, msi0_1503, \
                         msi1_1501, msi1_1502, msi1_1503, msk_1925, msk_1926, \
                         msk_1927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2402[k] = f_6 * msi0_1501[k]
                    - f_7 * msi1_1501[k]
                    + f_3 * pc_x[k] * msk_1925[k];

        t_2403[k] = f_6 * msi0_1502[k]
                    - f_7 * msi1_1502[k]
                    + f_3 * pc_x[k] * msk_1926[k];

        t_2404[k] = f_6 * msi0_1503[k]
                    - f_7 * msi1_1503[k]
                    + f_3 * pc_x[k] * msk_1927[k];
    }

#pragma omp simd aligned(t_2405, t_2406, t_2407, pa_y, pc_x, pc_y, lsl0_2000, lsl1_2000, \
                         msi0_1505, msi0_1506, msi1_1505, msi1_1506, msk_1929, \
                         msk_1930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2405[k] = pa_y[k] * lsl0_2000[k]
                    - f_14 * pc_y[k] * lsl1_2000[k];

        t_2406[k] = f_4 * msi0_1505[k]
                    - f_5 * msi1_1505[k]
                    + f_3 * pc_x[k] * msk_1929[k];

        t_2407[k] = f_4 * msi0_1506[k]
                    - f_5 * msi1_1506[k]
                    + f_3 * pc_x[k] * msk_1930[k];
    }

#pragma omp simd aligned(t_2408, t_2409, t_2410, pc_x, msi0_1507, msi0_1508, msi0_1509, \
                         msi1_1507, msi1_1508, msi1_1509, msk_1931, msk_1932, \
                         msk_1933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2408[k] = f_4 * msi0_1507[k]
                    - f_5 * msi1_1507[k]
                    + f_3 * pc_x[k] * msk_1931[k];

        t_2409[k] = f_4 * msi0_1508[k]
                    - f_5 * msi1_1508[k]
                    + f_3 * pc_x[k] * msk_1932[k];

        t_2410[k] = f_4 * msi0_1509[k]
                    - f_5 * msi1_1509[k]
                    + f_3 * pc_x[k] * msk_1933[k];
    }

#pragma omp simd aligned(t_2411, t_2412, t_2413, t_2414, t_2415, pa_y, pc_x, pc_y, lsl0_2007, \
                         lsl1_2007, msi0_1510, msi1_1510, msk_1934, msk_1936, msk_1937, \
                         msk_1938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2411[k] = f_4 * msi0_1510[k]
                    - f_5 * msi1_1510[k]
                    + f_3 * pc_x[k] * msk_1934[k];

        t_2412[k] = pa_y[k] * lsl0_2007[k]
                    - f_14 * pc_y[k] * lsl1_2007[k];

        t_2413[k] = f_3 * pc_x[k] * msk_1936[k];

        t_2414[k] = f_3 * pc_x[k] * msk_1937[k];

        t_2415[k] = f_3 * pc_x[k] * msk_1938[k];
    }

#pragma omp simd aligned(t_2416, t_2417, t_2418, t_2419, t_2420, pc_x, msk_1939, msk_1940, \
                         msk_1941, msk_1942, msk_1943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2416[k] = f_3 * pc_x[k] * msk_1939[k];

        t_2417[k] = f_3 * pc_x[k] * msk_1940[k];

        t_2418[k] = f_3 * pc_x[k] * msk_1941[k];

        t_2419[k] = f_3 * pc_x[k] * msk_1942[k];

        t_2420[k] = f_3 * pc_x[k] * msk_1943[k];
    }

#pragma omp simd aligned(t_2421, t_2422, t_2423, pa_y, pc_y, pc_z, lsl0_2016, lsl0_2018, \
                         lsk_1576, lsk_1612, lsk_1614, lsl1_2016, lsl1_2018, \
                         msk_1936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2421[k] = pa_y[k] * lsl0_2016[k]
                    + f_21 * lsk_1612[k]
                    - f_14 * pc_y[k] * lsl1_2016[k];

        t_2422[k] = f_21 * lsk_1576[k]
                    + f_3 * pc_z[k] * msk_1936[k];

        t_2423[k] = pa_y[k] * lsl0_2018[k]
                    + f_20 * lsk_1614[k]
                    - f_14 * pc_y[k] * lsl1_2018[k];
    }

#pragma omp simd aligned(t_2424, t_2425, t_2426, pa_y, pc_y, lsl0_2019, lsl0_2020, lsl0_2021, \
                         lsk_1615, lsk_1616, lsk_1617, lsl1_2019, lsl1_2020, \
                         lsl1_2021 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2424[k] = pa_y[k] * lsl0_2019[k]
                    + f_19 * lsk_1615[k]
                    - f_14 * pc_y[k] * lsl1_2019[k];

        t_2425[k] = pa_y[k] * lsl0_2020[k]
                    + f_18 * lsk_1616[k]
                    - f_14 * pc_y[k] * lsl1_2020[k];

        t_2426[k] = pa_y[k] * lsl0_2021[k]
                    + f_17 * lsk_1617[k]
                    - f_14 * pc_y[k] * lsl1_2021[k];
    }

#pragma omp simd aligned(t_2427, t_2428, t_2429, pa_y, pc_y, lsl0_2022, lsl0_2024, lsk_1618, \
                         lsk_1619, lsl1_2022, lsl1_2024, msk_1943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2427[k] = pa_y[k] * lsl0_2022[k]
                    + f_16 * lsk_1618[k]
                    - f_14 * pc_y[k] * lsl1_2022[k];

        t_2428[k] = f_15 * lsk_1619[k]
                    + f_3 * pc_y[k] * msk_1943[k];

        t_2429[k] = pa_y[k] * lsl0_2024[k]
                    - f_14 * pc_y[k] * lsl1_2024[k];
    }

#pragma omp simd aligned(t_2430, t_2431, t_2432, t_2433, t_2434, pc_x, pc_y, msi0_1512, \
                         msi0_1514, msi0_1515, msi1_1512, msi1_1514, msi1_1515, msk_1944, \
                         msk_1946, msk_1947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2430[k] = f_1 * msi0_1512[k]
                    - f_2 * msi1_1512[k]
                    + f_3 * pc_x[k] * msk_1944[k];

        t_2431[k] = f_3 * pc_y[k] * msk_1944[k];

        t_2432[k] = f_22 * msi0_1514[k]
                    - f_23 * msi1_1514[k]
                    + f_3 * pc_x[k] * msk_1946[k];

        t_2433[k] = f_12 * msi0_1515[k]
                    - f_13 * msi1_1515[k]
                    + f_3 * pc_x[k] * msk_1947[k];

        t_2434[k] = f_3 * pc_y[k] * msk_1946[k];
    }

#pragma omp simd aligned(t_2435, t_2436, t_2437, t_2438, pc_x, pc_y, msi0_1517, msi0_1518, \
                         msi0_1519, msi1_1517, msi1_1518, msi1_1519, msk_1949, msk_1950, \
                         msk_1951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2435[k] = f_12 * msi0_1517[k]
                    - f_13 * msi1_1517[k]
                    + f_3 * pc_x[k] * msk_1949[k];

        t_2436[k] = f_10 * msi0_1518[k]
                    - f_11 * msi1_1518[k]
                    + f_3 * pc_x[k] * msk_1950[k];

        t_2437[k] = f_10 * msi0_1519[k]
                    - f_11 * msi1_1519[k]
                    + f_3 * pc_x[k] * msk_1951[k];

        t_2438[k] = f_3 * pc_y[k] * msk_1949[k];
    }

#pragma omp simd aligned(t_2439, t_2440, t_2441, pc_x, msi0_1521, msi0_1522, msi0_1523, \
                         msi1_1521, msi1_1522, msi1_1523, msk_1953, msk_1954, \
                         msk_1955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2439[k] = f_10 * msi0_1521[k]
                    - f_11 * msi1_1521[k]
                    + f_3 * pc_x[k] * msk_1953[k];

        t_2440[k] = f_8 * msi0_1522[k]
                    - f_9 * msi1_1522[k]
                    + f_3 * pc_x[k] * msk_1954[k];

        t_2441[k] = f_8 * msi0_1523[k]
                    - f_9 * msi1_1523[k]
                    + f_3 * pc_x[k] * msk_1955[k];
    }

#pragma omp simd aligned(t_2442, t_2443, t_2444, t_2445, pc_x, pc_y, msi0_1524, msi0_1526, \
                         msi0_1527, msi1_1524, msi1_1526, msi1_1527, msk_1953, msk_1956, \
                         msk_1958, msk_1959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2442[k] = f_8 * msi0_1524[k]
                    - f_9 * msi1_1524[k]
                    + f_3 * pc_x[k] * msk_1956[k];

        t_2443[k] = f_3 * pc_y[k] * msk_1953[k];

        t_2444[k] = f_8 * msi0_1526[k]
                    - f_9 * msi1_1526[k]
                    + f_3 * pc_x[k] * msk_1958[k];

        t_2445[k] = f_6 * msi0_1527[k]
                    - f_7 * msi1_1527[k]
                    + f_3 * pc_x[k] * msk_1959[k];
    }

#pragma omp simd aligned(t_2446, t_2447, t_2448, t_2449, pc_x, pc_y, msi0_1528, msi0_1529, \
                         msi0_1530, msi1_1528, msi1_1529, msi1_1530, msk_1958, msk_1960, \
                         msk_1961, msk_1962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2446[k] = f_6 * msi0_1528[k]
                    - f_7 * msi1_1528[k]
                    + f_3 * pc_x[k] * msk_1960[k];

        t_2447[k] = f_6 * msi0_1529[k]
                    - f_7 * msi1_1529[k]
                    + f_3 * pc_x[k] * msk_1961[k];

        t_2448[k] = f_6 * msi0_1530[k]
                    - f_7 * msi1_1530[k]
                    + f_3 * pc_x[k] * msk_1962[k];

        t_2449[k] = f_3 * pc_y[k] * msk_1958[k];
    }

#pragma omp simd aligned(t_2450, t_2451, t_2452, pc_x, msi0_1532, msi0_1533, msi0_1534, \
                         msi1_1532, msi1_1533, msi1_1534, msk_1964, msk_1965, \
                         msk_1966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2450[k] = f_6 * msi0_1532[k]
                    - f_7 * msi1_1532[k]
                    + f_3 * pc_x[k] * msk_1964[k];

        t_2451[k] = f_4 * msi0_1533[k]
                    - f_5 * msi1_1533[k]
                    + f_3 * pc_x[k] * msk_1965[k];

        t_2452[k] = f_4 * msi0_1534[k]
                    - f_5 * msi1_1534[k]
                    + f_3 * pc_x[k] * msk_1966[k];
    }

#pragma omp simd aligned(t_2453, t_2454, t_2455, t_2456, pc_x, pc_y, msi0_1535, msi0_1536, \
                         msi0_1537, msi1_1535, msi1_1536, msi1_1537, msk_1964, msk_1967, \
                         msk_1968, msk_1969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2453[k] = f_4 * msi0_1535[k]
                    - f_5 * msi1_1535[k]
                    + f_3 * pc_x[k] * msk_1967[k];

        t_2454[k] = f_4 * msi0_1536[k]
                    - f_5 * msi1_1536[k]
                    + f_3 * pc_x[k] * msk_1968[k];

        t_2455[k] = f_4 * msi0_1537[k]
                    - f_5 * msi1_1537[k]
                    + f_3 * pc_x[k] * msk_1969[k];

        t_2456[k] = f_3 * pc_y[k] * msk_1964[k];
    }

#pragma omp simd aligned(t_2457, t_2458, t_2459, t_2460, t_2461, t_2462, pc_x, msi0_1539, \
                         msi1_1539, msk_1971, msk_1972, msk_1973, msk_1974, msk_1975, \
                         msk_1976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2457[k] = f_4 * msi0_1539[k]
                    - f_5 * msi1_1539[k]
                    + f_3 * pc_x[k] * msk_1971[k];

        t_2458[k] = f_3 * pc_x[k] * msk_1972[k];

        t_2459[k] = f_3 * pc_x[k] * msk_1973[k];

        t_2460[k] = f_3 * pc_x[k] * msk_1974[k];

        t_2461[k] = f_3 * pc_x[k] * msk_1975[k];

        t_2462[k] = f_3 * pc_x[k] * msk_1976[k];
    }

#pragma omp simd aligned(t_2463, t_2464, t_2465, t_2466, t_2467, pc_x, pc_y, msi0_1533, \
                         msi0_1534, msi1_1533, msi1_1534, msk_1972, msk_1973, msk_1977, \
                         msk_1978, msk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2463[k] = f_3 * pc_x[k] * msk_1977[k];

        t_2464[k] = f_3 * pc_x[k] * msk_1978[k];

        t_2465[k] = f_3 * pc_x[k] * msk_1979[k];

        t_2466[k] = f_1 * msi0_1533[k]
                    - f_2 * msi1_1533[k]
                    + f_3 * pc_y[k] * msk_1972[k];

        t_2467[k] = f_22 * msi0_1534[k]
                    - f_23 * msi1_1534[k]
                    + f_3 * pc_y[k] * msk_1973[k];
    }

#pragma omp simd aligned(t_2468, t_2469, t_2470, pc_y, msi0_1535, msi0_1536, msi0_1537, \
                         msi1_1535, msi1_1536, msi1_1537, msk_1974, msk_1975, \
                         msk_1976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2468[k] = f_12 * msi0_1535[k]
                    - f_13 * msi1_1535[k]
                    + f_3 * pc_y[k] * msk_1974[k];

        t_2469[k] = f_10 * msi0_1536[k]
                    - f_11 * msi1_1536[k]
                    + f_3 * pc_y[k] * msk_1975[k];

        t_2470[k] = f_8 * msi0_1537[k]
                    - f_9 * msi1_1537[k]
                    + f_3 * pc_y[k] * msk_1976[k];
    }

#pragma omp simd aligned(t_2471, t_2472, t_2473, t_2474, pc_y, pc_z, lsk_1619, msi0_1538, \
                         msi0_1539, msi1_1538, msi1_1539, msk_1977, msk_1978, \
                         msk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2471[k] = f_6 * msi0_1538[k]
                    - f_7 * msi1_1538[k]
                    + f_3 * pc_y[k] * msk_1977[k];

        t_2472[k] = f_4 * msi0_1539[k]
                    - f_5 * msi1_1539[k]
                    + f_3 * pc_y[k] * msk_1978[k];

        t_2473[k] = f_3 * pc_y[k] * msk_1979[k];

        t_2474[k] = f_0 * lsk_1619[k]
                    + f_1 * msi0_1539[k]
                    - f_2 * msi1_1539[k]
                    + f_3 * pc_z[k] * msk_1979[k];
    }
}

auto
compute_prim_msl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsl0, const size_t lsk,
                                                   const size_t lsl1, const size_t msi0,
                                                   const size_t msi1, const size_t msk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, lsl0, lsk,
                                                              lsl1, msi0, msi1, msk, ncols,
                                                              gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece9(buffer, target, pc, lsk, msi0,
                                                              msi1, msk, ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece12(buffer, target, pc, lsk, msi0,
                                                               msi1, msk, ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece14(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece15(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msk, ncols, gamma, p,
                                                               q);

    compute_prim_msl_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msk, ncols, gamma, p,
                                                               q);

    compute_prim_msl_three_center_electron_repulsion_0_piece17(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece18(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece19(buffer, target, pc, lsk, msi0,
                                                               msi1, msk, ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece20(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);

    compute_prim_msl_three_center_electron_repulsion_0_piece21(buffer, target, pa, pc, lsl0,
                                                               lsk, lsl1, msi0, msi1, msk,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
