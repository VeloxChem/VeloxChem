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


#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_18 = 4.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.5 / q;

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

    const auto *lsk0_0 = buffer.data(lsk0 + 0);
    const auto *lsk0_3 = buffer.data(lsk0 + 3);
    const auto *lsk0_5 = buffer.data(lsk0 + 5);
    const auto *lsk0_6 = buffer.data(lsk0 + 6);
    const auto *lsk0_9 = buffer.data(lsk0 + 9);
    const auto *lsk0_10 = buffer.data(lsk0 + 10);
    const auto *lsk0_14 = buffer.data(lsk0 + 14);
    const auto *lsk0_15 = buffer.data(lsk0 + 15);
    const auto *lsk0_20 = buffer.data(lsk0 + 20);
    const auto *lsk0_28 = buffer.data(lsk0 + 28);
    const auto *lsk0_35 = buffer.data(lsk0 + 35);

    const auto *lsi_0 = buffer.data(lsi + 0);
    const auto *lsi_1 = buffer.data(lsi + 1);
    const auto *lsi_2 = buffer.data(lsi + 2);
    const auto *lsi_3 = buffer.data(lsi + 3);
    const auto *lsi_5 = buffer.data(lsi + 5);
    const auto *lsi_6 = buffer.data(lsi + 6);
    const auto *lsi_9 = buffer.data(lsi + 9);
    const auto *lsi_10 = buffer.data(lsi + 10);
    const auto *lsi_14 = buffer.data(lsi + 14);
    const auto *lsi_21 = buffer.data(lsi + 21);
    const auto *lsi_23 = buffer.data(lsi + 23);
    const auto *lsi_24 = buffer.data(lsi + 24);
    const auto *lsi_25 = buffer.data(lsi + 25);
    const auto *lsi_27 = buffer.data(lsi + 27);
    const auto *lsi_28 = buffer.data(lsi + 28);
    const auto *lsi_33 = buffer.data(lsi + 33);
    const auto *lsi_37 = buffer.data(lsi + 37);
    const auto *lsi_42 = buffer.data(lsi + 42);
    const auto *lsi_49 = buffer.data(lsi + 49);
    const auto *lsi_51 = buffer.data(lsi + 51);
    const auto *lsi_52 = buffer.data(lsi + 52);
    const auto *lsi_53 = buffer.data(lsi + 53);
    const auto *lsi_54 = buffer.data(lsi + 54);
    const auto *lsi_55 = buffer.data(lsi + 55);
    const auto *lsi_77 = buffer.data(lsi + 77);
    const auto *lsi_78 = buffer.data(lsi + 78);
    const auto *lsi_79 = buffer.data(lsi + 79);
    const auto *lsi_80 = buffer.data(lsi + 80);
    const auto *lsi_81 = buffer.data(lsi + 81);
    const auto *lsi_83 = buffer.data(lsi + 83);
    const auto *lsi_84 = buffer.data(lsi + 84);
    const auto *lsi_87 = buffer.data(lsi + 87);
    const auto *lsi_90 = buffer.data(lsi + 90);
    const auto *lsi_94 = buffer.data(lsi + 94);
    const auto *lsi_99 = buffer.data(lsi + 99);

    const auto *lsk1_0 = buffer.data(lsk1 + 0);
    const auto *lsk1_3 = buffer.data(lsk1 + 3);
    const auto *lsk1_5 = buffer.data(lsk1 + 5);
    const auto *lsk1_6 = buffer.data(lsk1 + 6);
    const auto *lsk1_9 = buffer.data(lsk1 + 9);
    const auto *lsk1_10 = buffer.data(lsk1 + 10);
    const auto *lsk1_14 = buffer.data(lsk1 + 14);
    const auto *lsk1_15 = buffer.data(lsk1 + 15);
    const auto *lsk1_20 = buffer.data(lsk1 + 20);
    const auto *lsk1_28 = buffer.data(lsk1 + 28);
    const auto *lsk1_35 = buffer.data(lsk1 + 35);

    const auto *msh0_0 = buffer.data(msh0 + 0);
    const auto *msh0_1 = buffer.data(msh0 + 1);
    const auto *msh0_2 = buffer.data(msh0 + 2);
    const auto *msh0_3 = buffer.data(msh0 + 3);
    const auto *msh0_5 = buffer.data(msh0 + 5);
    const auto *msh0_6 = buffer.data(msh0 + 6);
    const auto *msh0_8 = buffer.data(msh0 + 8);
    const auto *msh0_9 = buffer.data(msh0 + 9);
    const auto *msh0_15 = buffer.data(msh0 + 15);
    const auto *msh0_17 = buffer.data(msh0 + 17);
    const auto *msh0_18 = buffer.data(msh0 + 18);
    const auto *msh0_19 = buffer.data(msh0 + 19);
    const auto *msh0_20 = buffer.data(msh0 + 20);
    const auto *msh0_24 = buffer.data(msh0 + 24);
    const auto *msh0_27 = buffer.data(msh0 + 27);
    const auto *msh0_28 = buffer.data(msh0 + 28);
    const auto *msh0_36 = buffer.data(msh0 + 36);
    const auto *msh0_37 = buffer.data(msh0 + 37);
    const auto *msh0_38 = buffer.data(msh0 + 38);
    const auto *msh0_39 = buffer.data(msh0 + 39);
    const auto *msh0_44 = buffer.data(msh0 + 44);
    const auto *msh0_46 = buffer.data(msh0 + 46);
    const auto *msh0_47 = buffer.data(msh0 + 47);
    const auto *msh0_49 = buffer.data(msh0 + 49);
    const auto *msh0_50 = buffer.data(msh0 + 50);
    const auto *msh0_51 = buffer.data(msh0 + 51);
    const auto *msh0_58 = buffer.data(msh0 + 58);
    const auto *msh0_59 = buffer.data(msh0 + 59);
    const auto *msh0_60 = buffer.data(msh0 + 60);
    const auto *msh0_61 = buffer.data(msh0 + 61);
    const auto *msh0_62 = buffer.data(msh0 + 62);
    const auto *msh0_63 = buffer.data(msh0 + 63);
    const auto *msh0_65 = buffer.data(msh0 + 65);
    const auto *msh0_66 = buffer.data(msh0 + 66);
    const auto *msh0_68 = buffer.data(msh0 + 68);
    const auto *msh0_69 = buffer.data(msh0 + 69);
    const auto *msh0_70 = buffer.data(msh0 + 70);
    const auto *msh0_72 = buffer.data(msh0 + 72);
    const auto *msh0_73 = buffer.data(msh0 + 73);
    const auto *msh0_78 = buffer.data(msh0 + 78);

    const auto *msh1_0 = buffer.data(msh1 + 0);
    const auto *msh1_1 = buffer.data(msh1 + 1);
    const auto *msh1_2 = buffer.data(msh1 + 2);
    const auto *msh1_3 = buffer.data(msh1 + 3);
    const auto *msh1_5 = buffer.data(msh1 + 5);
    const auto *msh1_6 = buffer.data(msh1 + 6);
    const auto *msh1_8 = buffer.data(msh1 + 8);
    const auto *msh1_9 = buffer.data(msh1 + 9);
    const auto *msh1_15 = buffer.data(msh1 + 15);
    const auto *msh1_17 = buffer.data(msh1 + 17);
    const auto *msh1_18 = buffer.data(msh1 + 18);
    const auto *msh1_19 = buffer.data(msh1 + 19);
    const auto *msh1_20 = buffer.data(msh1 + 20);
    const auto *msh1_24 = buffer.data(msh1 + 24);
    const auto *msh1_27 = buffer.data(msh1 + 27);
    const auto *msh1_28 = buffer.data(msh1 + 28);
    const auto *msh1_36 = buffer.data(msh1 + 36);
    const auto *msh1_37 = buffer.data(msh1 + 37);
    const auto *msh1_38 = buffer.data(msh1 + 38);
    const auto *msh1_39 = buffer.data(msh1 + 39);
    const auto *msh1_44 = buffer.data(msh1 + 44);
    const auto *msh1_46 = buffer.data(msh1 + 46);
    const auto *msh1_47 = buffer.data(msh1 + 47);
    const auto *msh1_49 = buffer.data(msh1 + 49);
    const auto *msh1_50 = buffer.data(msh1 + 50);
    const auto *msh1_51 = buffer.data(msh1 + 51);
    const auto *msh1_58 = buffer.data(msh1 + 58);
    const auto *msh1_59 = buffer.data(msh1 + 59);
    const auto *msh1_60 = buffer.data(msh1 + 60);
    const auto *msh1_61 = buffer.data(msh1 + 61);
    const auto *msh1_62 = buffer.data(msh1 + 62);
    const auto *msh1_63 = buffer.data(msh1 + 63);
    const auto *msh1_65 = buffer.data(msh1 + 65);
    const auto *msh1_66 = buffer.data(msh1 + 66);
    const auto *msh1_68 = buffer.data(msh1 + 68);
    const auto *msh1_69 = buffer.data(msh1 + 69);
    const auto *msh1_70 = buffer.data(msh1 + 70);
    const auto *msh1_72 = buffer.data(msh1 + 72);
    const auto *msh1_73 = buffer.data(msh1 + 73);
    const auto *msh1_78 = buffer.data(msh1 + 78);

    const auto *msi_0 = buffer.data(msi + 0);
    const auto *msi_1 = buffer.data(msi + 1);
    const auto *msi_2 = buffer.data(msi + 2);
    const auto *msi_3 = buffer.data(msi + 3);
    const auto *msi_5 = buffer.data(msi + 5);
    const auto *msi_6 = buffer.data(msi + 6);
    const auto *msi_8 = buffer.data(msi + 8);
    const auto *msi_9 = buffer.data(msi + 9);
    const auto *msi_10 = buffer.data(msi + 10);
    const auto *msi_12 = buffer.data(msi + 12);
    const auto *msi_13 = buffer.data(msi + 13);
    const auto *msi_14 = buffer.data(msi + 14);
    const auto *msi_15 = buffer.data(msi + 15);
    const auto *msi_20 = buffer.data(msi + 20);
    const auto *msi_21 = buffer.data(msi + 21);
    const auto *msi_23 = buffer.data(msi + 23);
    const auto *msi_24 = buffer.data(msi + 24);
    const auto *msi_25 = buffer.data(msi + 25);
    const auto *msi_26 = buffer.data(msi + 26);
    const auto *msi_27 = buffer.data(msi + 27);
    const auto *msi_28 = buffer.data(msi + 28);
    const auto *msi_29 = buffer.data(msi + 29);
    const auto *msi_31 = buffer.data(msi + 31);
    const auto *msi_33 = buffer.data(msi + 33);
    const auto *msi_34 = buffer.data(msi + 34);
    const auto *msi_35 = buffer.data(msi + 35);
    const auto *msi_37 = buffer.data(msi + 37);
    const auto *msi_38 = buffer.data(msi + 38);
    const auto *msi_39 = buffer.data(msi + 39);
    const auto *msi_40 = buffer.data(msi + 40);
    const auto *msi_42 = buffer.data(msi + 42);
    const auto *msi_43 = buffer.data(msi + 43);
    const auto *msi_49 = buffer.data(msi + 49);
    const auto *msi_50 = buffer.data(msi + 50);
    const auto *msi_51 = buffer.data(msi + 51);
    const auto *msi_52 = buffer.data(msi + 52);
    const auto *msi_53 = buffer.data(msi + 53);
    const auto *msi_54 = buffer.data(msi + 54);
    const auto *msi_55 = buffer.data(msi + 55);
    const auto *msi_56 = buffer.data(msi + 56);
    const auto *msi_58 = buffer.data(msi + 58);
    const auto *msi_60 = buffer.data(msi + 60);
    const auto *msi_61 = buffer.data(msi + 61);
    const auto *msi_63 = buffer.data(msi + 63);
    const auto *msi_64 = buffer.data(msi + 64);
    const auto *msi_65 = buffer.data(msi + 65);
    const auto *msi_67 = buffer.data(msi + 67);
    const auto *msi_68 = buffer.data(msi + 68);
    const auto *msi_69 = buffer.data(msi + 69);
    const auto *msi_70 = buffer.data(msi + 70);
    const auto *msi_76 = buffer.data(msi + 76);
    const auto *msi_77 = buffer.data(msi + 77);
    const auto *msi_78 = buffer.data(msi + 78);
    const auto *msi_79 = buffer.data(msi + 79);
    const auto *msi_80 = buffer.data(msi + 80);
    const auto *msi_81 = buffer.data(msi + 81);
    const auto *msi_82 = buffer.data(msi + 82);
    const auto *msi_83 = buffer.data(msi + 83);
    const auto *msi_84 = buffer.data(msi + 84);
    const auto *msi_85 = buffer.data(msi + 85);
    const auto *msi_86 = buffer.data(msi + 86);
    const auto *msi_87 = buffer.data(msi + 87);
    const auto *msi_89 = buffer.data(msi + 89);
    const auto *msi_90 = buffer.data(msi + 90);
    const auto *msi_91 = buffer.data(msi + 91);
    const auto *msi_93 = buffer.data(msi + 93);
    const auto *msi_94 = buffer.data(msi + 94);
    const auto *msi_95 = buffer.data(msi + 95);
    const auto *msi_96 = buffer.data(msi + 96);
    const auto *msi_98 = buffer.data(msi + 98);
    const auto *msi_99 = buffer.data(msi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, lsi_0, msh0_0, \
                         msh1_0, msi_0, msi_1, msi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsi_0[k]
                 + f_1 * msh0_0[k]
                 - f_2 * msh1_0[k]
                 + f_3 * pc_x[k] * msi_0[k];

        t_1[k] = f_3 * pc_y[k] * msi_0[k];

        t_2[k] = f_3 * pc_z[k] * msi_0[k];

        t_3[k] = f_4 * msh0_0[k]
                 - f_5 * msh1_0[k]
                 + f_3 * pc_y[k] * msi_1[k];

        t_4[k] = f_3 * pc_y[k] * msi_2[k];

        t_5[k] = f_4 * msh0_0[k]
                 - f_5 * msh1_0[k]
                 + f_3 * pc_z[k] * msi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, msh0_1, msh0_2, msh0_3, msh1_1, \
                         msh1_2, msh1_3, msi_3, msi_5, msi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * msh0_1[k]
                 - f_7 * msh1_1[k]
                 + f_3 * pc_y[k] * msi_3[k];

        t_7[k] = f_3 * pc_z[k] * msi_3[k];

        t_8[k] = f_3 * pc_y[k] * msi_5[k];

        t_9[k] = f_6 * msh0_2[k]
                 - f_7 * msh1_2[k]
                 + f_3 * pc_z[k] * msi_5[k];

        t_10[k] = f_8 * msh0_3[k]
                  - f_9 * msh1_3[k]
                  + f_3 * pc_y[k] * msi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, msh0_5, msh0_6, \
                         msh1_5, msh1_6, msi_6, msi_8, msi_9, msi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * msi_6[k];

        t_12[k] = f_4 * msh0_5[k]
                  - f_5 * msh1_5[k]
                  + f_3 * pc_y[k] * msi_8[k];

        t_13[k] = f_3 * pc_y[k] * msi_9[k];

        t_14[k] = f_8 * msh0_5[k]
                  - f_9 * msh1_5[k]
                  + f_3 * pc_z[k] * msi_9[k];

        t_15[k] = f_10 * msh0_6[k]
                  - f_11 * msh1_6[k]
                  + f_3 * pc_y[k] * msi_10[k];

        t_16[k] = f_3 * pc_z[k] * msi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, msh0_8, msh0_9, msh1_8, msh1_9, \
                         msi_12, msi_13, msi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * msh0_8[k]
                  - f_7 * msh1_8[k]
                  + f_3 * pc_y[k] * msi_12[k];

        t_18[k] = f_4 * msh0_9[k]
                  - f_5 * msh1_9[k]
                  + f_3 * pc_y[k] * msi_13[k];

        t_19[k] = f_3 * pc_y[k] * msi_14[k];

        t_20[k] = f_10 * msh0_9[k]
                  - f_11 * msh1_9[k]
                  + f_3 * pc_z[k] * msi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, lsi_21, lsi_23, lsi_24, \
                         lsi_25, msi_15, msi_21, msi_23, msi_24, \
                         msi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * lsi_21[k]
                  + f_3 * pc_x[k] * msi_21[k];

        t_22[k] = f_3 * pc_z[k] * msi_15[k];

        t_23[k] = f_0 * lsi_23[k]
                  + f_3 * pc_x[k] * msi_23[k];

        t_24[k] = f_0 * lsi_24[k]
                  + f_3 * pc_x[k] * msi_24[k];

        t_25[k] = f_0 * lsi_25[k]
                  + f_3 * pc_x[k] * msi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, lsi_27, msh0_15, msh1_15, \
                         msi_20, msi_21, msi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * msi_20[k];

        t_27[k] = f_0 * lsi_27[k]
                  + f_3 * pc_x[k] * msi_27[k];

        t_28[k] = f_1 * msh0_15[k]
                  - f_2 * msh1_15[k]
                  + f_3 * pc_y[k] * msi_21[k];

        t_29[k] = f_3 * pc_z[k] * msi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, msh0_17, msh0_18, msh0_19, msh1_17, msh1_18, \
                         msh1_19, msi_23, msi_24, msi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * msh0_17[k]
                  - f_11 * msh1_17[k]
                  + f_3 * pc_y[k] * msi_23[k];

        t_31[k] = f_8 * msh0_18[k]
                  - f_9 * msh1_18[k]
                  + f_3 * pc_y[k] * msi_24[k];

        t_32[k] = f_6 * msh0_19[k]
                  - f_7 * msh1_19[k]
                  + f_3 * pc_y[k] * msi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, lsk0_0, lsi_0, \
                         lsk1_0, msh0_20, msh1_20, msi_26, msi_27, \
                         msi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * msh0_20[k]
                  - f_5 * msh1_20[k]
                  + f_3 * pc_y[k] * msi_26[k];

        t_34[k] = f_3 * pc_y[k] * msi_27[k];

        t_35[k] = f_1 * msh0_20[k]
                  - f_2 * msh1_20[k]
                  + f_3 * pc_z[k] * msi_27[k];

        t_36[k] = pa_y[k] * lsk0_0[k]
                  - f_12 * pc_y[k] * lsk1_0[k];

        t_37[k] = f_13 * lsi_0[k]
                  + f_3 * pc_y[k] * msi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, lsk0_3, lsk0_5, lsi_1, \
                         lsk1_3, lsk1_5, msi_28, msi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * msi_28[k];

        t_39[k] = pa_y[k] * lsk0_3[k]
                  + f_14 * lsi_1[k]
                  - f_12 * pc_y[k] * lsk1_3[k];

        t_40[k] = f_3 * pc_z[k] * msi_29[k];

        t_41[k] = pa_y[k] * lsk0_5[k]
                  - f_12 * pc_y[k] * lsk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, lsk0_6, lsk0_9, lsi_3, \
                         lsi_5, lsk1_6, lsk1_9, msi_31, msi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * lsk0_6[k]
                  + f_15 * lsi_3[k]
                  - f_12 * pc_y[k] * lsk1_6[k];

        t_43[k] = f_3 * pc_z[k] * msi_31[k];

        t_44[k] = f_13 * lsi_5[k]
                  + f_3 * pc_y[k] * msi_33[k];

        t_45[k] = pa_y[k] * lsk0_9[k]
                  - f_12 * pc_y[k] * lsk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, lsk0_10, lsi_6, lsi_9, \
                         lsk1_10, msh0_24, msh1_24, msi_34, msi_35, \
                         msi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * lsk0_10[k]
                  + f_16 * lsi_6[k]
                  - f_12 * pc_y[k] * lsk1_10[k];

        t_47[k] = f_3 * pc_z[k] * msi_34[k];

        t_48[k] = f_4 * msh0_24[k]
                  - f_5 * msh1_24[k]
                  + f_3 * pc_z[k] * msi_35[k];

        t_49[k] = f_13 * lsi_9[k]
                  + f_3 * pc_y[k] * msi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, lsk0_14, lsk0_15, lsi_10, \
                         lsk1_14, lsk1_15, msh0_27, msh1_27, msi_38, \
                         msi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * lsk0_14[k]
                  - f_12 * pc_y[k] * lsk1_14[k];

        t_51[k] = pa_y[k] * lsk0_15[k]
                  + f_17 * lsi_10[k]
                  - f_12 * pc_y[k] * lsk1_15[k];

        t_52[k] = f_3 * pc_z[k] * msi_38[k];

        t_53[k] = f_4 * msh0_27[k]
                  - f_5 * msh1_27[k]
                  + f_3 * pc_z[k] * msi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, lsk0_20, lsi_14, lsk1_20, \
                         msh0_28, msh1_28, msi_40, msi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * msh0_28[k]
                  - f_7 * msh1_28[k]
                  + f_3 * pc_z[k] * msi_40[k];

        t_55[k] = f_13 * lsi_14[k]
                  + f_3 * pc_y[k] * msi_42[k];

        t_56[k] = pa_y[k] * lsk0_20[k]
                  - f_12 * pc_y[k] * lsk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, lsi_49, lsi_51, lsi_52, \
                         lsi_53, msi_43, msi_49, msi_51, msi_52, \
                         msi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_18 * lsi_49[k]
                  + f_3 * pc_x[k] * msi_49[k];

        t_58[k] = f_3 * pc_z[k] * msi_43[k];

        t_59[k] = f_18 * lsi_51[k]
                  + f_3 * pc_x[k] * msi_51[k];

        t_60[k] = f_18 * lsi_52[k]
                  + f_3 * pc_x[k] * msi_52[k];

        t_61[k] = f_18 * lsi_53[k]
                  + f_3 * pc_x[k] * msi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, lsi_21, lsi_54, lsi_55, \
                         msh0_36, msh1_36, msi_49, msi_54, msi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_18 * lsi_54[k]
                  + f_3 * pc_x[k] * msi_54[k];

        t_63[k] = f_18 * lsi_55[k]
                  + f_3 * pc_x[k] * msi_55[k];

        t_64[k] = f_13 * lsi_21[k]
                  + f_1 * msh0_36[k]
                  - f_2 * msh1_36[k]
                  + f_3 * pc_y[k] * msi_49[k];

        t_65[k] = f_3 * pc_z[k] * msi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, msh0_36, msh0_37, msh0_38, msh1_36, msh1_37, \
                         msh1_38, msi_50, msi_51, msi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * msh0_36[k]
                  - f_5 * msh1_36[k]
                  + f_3 * pc_z[k] * msi_50[k];

        t_67[k] = f_6 * msh0_37[k]
                  - f_7 * msh1_37[k]
                  + f_3 * pc_z[k] * msi_51[k];

        t_68[k] = f_8 * msh0_38[k]
                  - f_9 * msh1_38[k]
                  + f_3 * pc_z[k] * msi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, lsk0_35, lsi_27, lsk1_35, \
                         msh0_39, msh1_39, msi_53, msi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * msh0_39[k]
                  - f_11 * msh1_39[k]
                  + f_3 * pc_z[k] * msi_53[k];

        t_70[k] = f_13 * lsi_27[k]
                  + f_3 * pc_y[k] * msi_55[k];

        t_71[k] = pa_y[k] * lsk0_35[k]
                  - f_12 * pc_y[k] * lsk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, lsk0_0, lsk0_3, \
                         lsi_0, lsk1_0, lsk1_3, msi_56, msi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * lsk0_0[k]
                  - f_12 * pc_z[k] * lsk1_0[k];

        t_73[k] = f_3 * pc_y[k] * msi_56[k];

        t_74[k] = f_13 * lsi_0[k]
                  + f_3 * pc_z[k] * msi_56[k];

        t_75[k] = pa_z[k] * lsk0_3[k]
                  - f_12 * pc_z[k] * lsk1_3[k];

        t_76[k] = f_3 * pc_y[k] * msi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, lsk0_5, lsk0_6, lsi_2, \
                         lsk1_5, lsk1_6, msh0_44, msh1_44, msi_60, \
                         msi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * lsk0_5[k]
                  + f_14 * lsi_2[k]
                  - f_12 * pc_z[k] * lsk1_5[k];

        t_78[k] = pa_z[k] * lsk0_6[k]
                  - f_12 * pc_z[k] * lsk1_6[k];

        t_79[k] = f_4 * msh0_44[k]
                  - f_5 * msh1_44[k]
                  + f_3 * pc_y[k] * msi_60[k];

        t_80[k] = f_3 * pc_y[k] * msi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, lsk0_9, lsk0_10, lsi_5, lsk1_9, \
                         lsk1_10, msh0_46, msh1_46, msi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * lsk0_9[k]
                  + f_15 * lsi_5[k]
                  - f_12 * pc_z[k] * lsk1_9[k];

        t_82[k] = pa_z[k] * lsk0_10[k]
                  - f_12 * pc_z[k] * lsk1_10[k];

        t_83[k] = f_6 * msh0_46[k]
                  - f_7 * msh1_46[k]
                  + f_3 * pc_y[k] * msi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, lsk0_14, lsk0_15, lsi_9, \
                         lsk1_14, lsk1_15, msh0_47, msh1_47, msi_64, \
                         msi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * msh0_47[k]
                  - f_5 * msh1_47[k]
                  + f_3 * pc_y[k] * msi_64[k];

        t_85[k] = f_3 * pc_y[k] * msi_65[k];

        t_86[k] = pa_z[k] * lsk0_14[k]
                  + f_16 * lsi_9[k]
                  - f_12 * pc_z[k] * lsk1_14[k];

        t_87[k] = pa_z[k] * lsk0_15[k]
                  - f_12 * pc_z[k] * lsk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, msh0_49, msh0_50, msh0_51, msh1_49, \
                         msh1_50, msh1_51, msi_67, msi_68, msi_69, \
                         msi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * msh0_49[k]
                  - f_9 * msh1_49[k]
                  + f_3 * pc_y[k] * msi_67[k];

        t_89[k] = f_6 * msh0_50[k]
                  - f_7 * msh1_50[k]
                  + f_3 * pc_y[k] * msi_68[k];

        t_90[k] = f_4 * msh0_51[k]
                  - f_5 * msh1_51[k]
                  + f_3 * pc_y[k] * msi_69[k];

        t_91[k] = f_3 * pc_y[k] * msi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, lsk0_20, lsi_14, lsi_77, \
                         lsi_78, lsi_79, lsk1_20, msi_77, msi_78, \
                         msi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * lsk0_20[k]
                  + f_17 * lsi_14[k]
                  - f_12 * pc_z[k] * lsk1_20[k];

        t_93[k] = f_18 * lsi_77[k]
                  + f_3 * pc_x[k] * msi_77[k];

        t_94[k] = f_18 * lsi_78[k]
                  + f_3 * pc_x[k] * msi_78[k];

        t_95[k] = f_18 * lsi_79[k]
                  + f_3 * pc_x[k] * msi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, lsi_80, lsi_81, lsi_83, msi_76, \
                         msi_80, msi_81, msi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_18 * lsi_80[k]
                  + f_3 * pc_x[k] * msi_80[k];

        t_97[k] = f_18 * lsi_81[k]
                  + f_3 * pc_x[k] * msi_81[k];

        t_98[k] = f_3 * pc_y[k] * msi_76[k];

        t_99[k] = f_18 * lsi_83[k]
                  + f_3 * pc_x[k] * msi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, lsk0_28, lsk1_28, msh0_58, \
                         msh0_59, msh1_58, msh1_59, msi_78, msi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * lsk0_28[k]
                   - f_12 * pc_z[k] * lsk1_28[k];

        t_101[k] = f_19 * msh0_58[k]
                   - f_20 * msh1_58[k]
                   + f_3 * pc_y[k] * msi_78[k];

        t_102[k] = f_10 * msh0_59[k]
                   - f_11 * msh1_59[k]
                   + f_3 * pc_y[k] * msi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, msh0_60, msh0_61, msh0_62, msh1_60, \
                         msh1_61, msh1_62, msi_80, msi_81, msi_82, \
                         msi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * msh0_60[k]
                   - f_9 * msh1_60[k]
                   + f_3 * pc_y[k] * msi_80[k];

        t_104[k] = f_6 * msh0_61[k]
                   - f_7 * msh1_61[k]
                   + f_3 * pc_y[k] * msi_81[k];

        t_105[k] = f_4 * msh0_62[k]
                   - f_5 * msh1_62[k]
                   + f_3 * pc_y[k] * msi_82[k];

        t_106[k] = f_3 * pc_y[k] * msi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, lsi_27, lsi_28, lsi_84, \
                         msh0_62, msh0_63, msh1_62, msh1_63, msi_83, \
                         msi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * lsi_27[k]
                   + f_1 * msh0_62[k]
                   - f_2 * msh1_62[k]
                   + f_3 * pc_z[k] * msi_83[k];

        t_108[k] = f_21 * lsi_84[k]
                   + f_1 * msh0_63[k]
                   - f_2 * msh1_63[k]
                   + f_3 * pc_x[k] * msi_84[k];

        t_109[k] = f_14 * lsi_28[k]
                   + f_3 * pc_y[k] * msi_84[k];

        t_110[k] = f_3 * pc_z[k] * msi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, lsi_87, msh0_63, msh0_66, msh1_63, \
                         msh1_66, msi_85, msi_86, msi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_21 * lsi_87[k]
                   + f_10 * msh0_66[k]
                   - f_11 * msh1_66[k]
                   + f_3 * pc_x[k] * msi_87[k];

        t_112[k] = f_3 * pc_z[k] * msi_85[k];

        t_113[k] = f_4 * msh0_63[k]
                   - f_5 * msh1_63[k]
                   + f_3 * pc_z[k] * msi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, lsi_33, lsi_90, \
                         msh0_65, msh0_69, msh1_65, msh1_69, msi_87, msi_89, \
                         msi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_21 * lsi_90[k]
                   + f_8 * msh0_69[k]
                   - f_9 * msh1_69[k]
                   + f_3 * pc_x[k] * msi_90[k];

        t_115[k] = f_3 * pc_z[k] * msi_87[k];

        t_116[k] = f_14 * lsi_33[k]
                   + f_3 * pc_y[k] * msi_89[k];

        t_117[k] = f_6 * msh0_65[k]
                   - f_7 * msh1_65[k]
                   + f_3 * pc_z[k] * msi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, lsi_94, msh0_66, msh0_73, msh1_66, \
                         msh1_73, msi_90, msi_91, msi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_21 * lsi_94[k]
                   + f_6 * msh0_73[k]
                   - f_7 * msh1_73[k]
                   + f_3 * pc_x[k] * msi_94[k];

        t_119[k] = f_3 * pc_z[k] * msi_90[k];

        t_120[k] = f_4 * msh0_66[k]
                   - f_5 * msh1_66[k]
                   + f_3 * pc_z[k] * msi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, lsi_37, lsi_99, \
                         msh0_68, msh0_78, msh1_68, msh1_78, msi_93, msi_94, \
                         msi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * lsi_37[k]
                   + f_3 * pc_y[k] * msi_93[k];

        t_122[k] = f_8 * msh0_68[k]
                   - f_9 * msh1_68[k]
                   + f_3 * pc_z[k] * msi_93[k];

        t_123[k] = f_21 * lsi_99[k]
                   + f_4 * msh0_78[k]
                   - f_5 * msh1_78[k]
                   + f_3 * pc_x[k] * msi_99[k];

        t_124[k] = f_3 * pc_z[k] * msi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, lsi_42, msh0_69, msh0_70, \
                         msh0_72, msh1_69, msh1_70, msh1_72, msi_95, msi_96, \
                         msi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * msh0_69[k]
                   - f_5 * msh1_69[k]
                   + f_3 * pc_z[k] * msi_95[k];

        t_126[k] = f_6 * msh0_70[k]
                   - f_7 * msh1_70[k]
                   + f_3 * pc_z[k] * msi_96[k];

        t_127[k] = f_14 * lsi_42[k]
                   + f_3 * pc_y[k] * msi_98[k];

        t_128[k] = f_10 * msh0_72[k]
                   - f_11 * msh1_72[k]
                   + f_3 * pc_z[k] * msi_98[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / q;

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

    const auto *lsk0_39 = buffer.data(lsk0 + 39);
    const auto *lsk0_42 = buffer.data(lsk0 + 42);
    const auto *lsk0_46 = buffer.data(lsk0 + 46);
    const auto *lsk0_51 = buffer.data(lsk0 + 51);
    const auto *lsk0_64 = buffer.data(lsk0 + 64);
    const auto *lsk0_72 = buffer.data(lsk0 + 72);
    const auto *lsk0_77 = buffer.data(lsk0 + 77);
    const auto *lsk0_81 = buffer.data(lsk0 + 81);
    const auto *lsk0_84 = buffer.data(lsk0 + 84);
    const auto *lsk0_86 = buffer.data(lsk0 + 86);
    const auto *lsk0_89 = buffer.data(lsk0 + 89);
    const auto *lsk0_90 = buffer.data(lsk0 + 90);
    const auto *lsk0_92 = buffer.data(lsk0 + 92);
    const auto *lsk0_107 = buffer.data(lsk0 + 107);

    const auto *lsi_28 = buffer.data(lsi + 28);
    const auto *lsi_31 = buffer.data(lsi + 31);
    const auto *lsi_34 = buffer.data(lsi + 34);
    const auto *lsi_38 = buffer.data(lsi + 38);
    const auto *lsi_49 = buffer.data(lsi + 49);
    const auto *lsi_55 = buffer.data(lsi + 55);
    const auto *lsi_56 = buffer.data(lsi + 56);
    const auto *lsi_58 = buffer.data(lsi + 58);
    const auto *lsi_61 = buffer.data(lsi + 61);
    const auto *lsi_64 = buffer.data(lsi + 64);
    const auto *lsi_65 = buffer.data(lsi + 65);
    const auto *lsi_68 = buffer.data(lsi + 68);
    const auto *lsi_69 = buffer.data(lsi + 69);
    const auto *lsi_70 = buffer.data(lsi + 70);
    const auto *lsi_79 = buffer.data(lsi + 79);
    const auto *lsi_80 = buffer.data(lsi + 80);
    const auto *lsi_81 = buffer.data(lsi + 81);
    const auto *lsi_82 = buffer.data(lsi + 82);
    const auto *lsi_83 = buffer.data(lsi + 83);
    const auto *lsi_84 = buffer.data(lsi + 84);
    const auto *lsi_89 = buffer.data(lsi + 89);
    const auto *lsi_93 = buffer.data(lsi + 93);
    const auto *lsi_98 = buffer.data(lsi + 98);
    const auto *lsi_105 = buffer.data(lsi + 105);
    const auto *lsi_107 = buffer.data(lsi + 107);
    const auto *lsi_108 = buffer.data(lsi + 108);
    const auto *lsi_109 = buffer.data(lsi + 109);
    const auto *lsi_110 = buffer.data(lsi + 110);
    const auto *lsi_111 = buffer.data(lsi + 111);
    const auto *lsi_133 = buffer.data(lsi + 133);
    const auto *lsi_134 = buffer.data(lsi + 134);
    const auto *lsi_135 = buffer.data(lsi + 135);
    const auto *lsi_136 = buffer.data(lsi + 136);
    const auto *lsi_137 = buffer.data(lsi + 137);
    const auto *lsi_138 = buffer.data(lsi + 138);
    const auto *lsi_139 = buffer.data(lsi + 139);
    const auto *lsi_140 = buffer.data(lsi + 140);
    const auto *lsi_145 = buffer.data(lsi + 145);
    const auto *lsi_149 = buffer.data(lsi + 149);
    const auto *lsi_154 = buffer.data(lsi + 154);
    const auto *lsi_160 = buffer.data(lsi + 160);
    const auto *lsi_161 = buffer.data(lsi + 161);
    const auto *lsi_162 = buffer.data(lsi + 162);
    const auto *lsi_163 = buffer.data(lsi + 163);
    const auto *lsi_164 = buffer.data(lsi + 164);
    const auto *lsi_165 = buffer.data(lsi + 165);
    const auto *lsi_167 = buffer.data(lsi + 167);
    const auto *lsi_168 = buffer.data(lsi + 168);
    const auto *lsi_171 = buffer.data(lsi + 171);
    const auto *lsi_174 = buffer.data(lsi + 174);
    const auto *lsi_178 = buffer.data(lsi + 178);
    const auto *lsi_183 = buffer.data(lsi + 183);
    const auto *lsi_189 = buffer.data(lsi + 189);
    const auto *lsi_191 = buffer.data(lsi + 191);
    const auto *lsi_192 = buffer.data(lsi + 192);
    const auto *lsi_193 = buffer.data(lsi + 193);
    const auto *lsi_194 = buffer.data(lsi + 194);
    const auto *lsi_195 = buffer.data(lsi + 195);

    const auto *lsk1_39 = buffer.data(lsk1 + 39);
    const auto *lsk1_42 = buffer.data(lsk1 + 42);
    const auto *lsk1_46 = buffer.data(lsk1 + 46);
    const auto *lsk1_51 = buffer.data(lsk1 + 51);
    const auto *lsk1_64 = buffer.data(lsk1 + 64);
    const auto *lsk1_72 = buffer.data(lsk1 + 72);
    const auto *lsk1_77 = buffer.data(lsk1 + 77);
    const auto *lsk1_81 = buffer.data(lsk1 + 81);
    const auto *lsk1_84 = buffer.data(lsk1 + 84);
    const auto *lsk1_86 = buffer.data(lsk1 + 86);
    const auto *lsk1_89 = buffer.data(lsk1 + 89);
    const auto *lsk1_90 = buffer.data(lsk1 + 90);
    const auto *lsk1_92 = buffer.data(lsk1 + 92);
    const auto *lsk1_107 = buffer.data(lsk1 + 107);

    const auto *msh0_78 = buffer.data(msh0 + 78);
    const auto *msh0_79 = buffer.data(msh0 + 79);
    const auto *msh0_80 = buffer.data(msh0 + 80);
    const auto *msh0_81 = buffer.data(msh0 + 81);
    const auto *msh0_83 = buffer.data(msh0 + 83);
    const auto *msh0_101 = buffer.data(msh0 + 101);
    const auto *msh0_102 = buffer.data(msh0 + 102);
    const auto *msh0_103 = buffer.data(msh0 + 103);
    const auto *msh0_104 = buffer.data(msh0 + 104);
    const auto *msh0_105 = buffer.data(msh0 + 105);
    const auto *msh0_106 = buffer.data(msh0 + 106);
    const auto *msh0_107 = buffer.data(msh0 + 107);
    const auto *msh0_108 = buffer.data(msh0 + 108);
    const auto *msh0_109 = buffer.data(msh0 + 109);
    const auto *msh0_110 = buffer.data(msh0 + 110);
    const auto *msh0_111 = buffer.data(msh0 + 111);
    const auto *msh0_112 = buffer.data(msh0 + 112);
    const auto *msh0_113 = buffer.data(msh0 + 113);
    const auto *msh0_114 = buffer.data(msh0 + 114);
    const auto *msh0_119 = buffer.data(msh0 + 119);
    const auto *msh0_120 = buffer.data(msh0 + 120);
    const auto *msh0_121 = buffer.data(msh0 + 121);
    const auto *msh0_122 = buffer.data(msh0 + 122);
    const auto *msh0_123 = buffer.data(msh0 + 123);
    const auto *msh0_124 = buffer.data(msh0 + 124);
    const auto *msh0_125 = buffer.data(msh0 + 125);
    const auto *msh0_126 = buffer.data(msh0 + 126);
    const auto *msh0_128 = buffer.data(msh0 + 128);
    const auto *msh0_129 = buffer.data(msh0 + 129);
    const auto *msh0_131 = buffer.data(msh0 + 131);
    const auto *msh0_132 = buffer.data(msh0 + 132);
    const auto *msh0_133 = buffer.data(msh0 + 133);
    const auto *msh0_135 = buffer.data(msh0 + 135);
    const auto *msh0_136 = buffer.data(msh0 + 136);
    const auto *msh0_141 = buffer.data(msh0 + 141);
    const auto *msh0_142 = buffer.data(msh0 + 142);
    const auto *msh0_143 = buffer.data(msh0 + 143);
    const auto *msh0_144 = buffer.data(msh0 + 144);

    const auto *msh1_78 = buffer.data(msh1 + 78);
    const auto *msh1_79 = buffer.data(msh1 + 79);
    const auto *msh1_80 = buffer.data(msh1 + 80);
    const auto *msh1_81 = buffer.data(msh1 + 81);
    const auto *msh1_83 = buffer.data(msh1 + 83);
    const auto *msh1_101 = buffer.data(msh1 + 101);
    const auto *msh1_102 = buffer.data(msh1 + 102);
    const auto *msh1_103 = buffer.data(msh1 + 103);
    const auto *msh1_104 = buffer.data(msh1 + 104);
    const auto *msh1_105 = buffer.data(msh1 + 105);
    const auto *msh1_106 = buffer.data(msh1 + 106);
    const auto *msh1_107 = buffer.data(msh1 + 107);
    const auto *msh1_108 = buffer.data(msh1 + 108);
    const auto *msh1_109 = buffer.data(msh1 + 109);
    const auto *msh1_110 = buffer.data(msh1 + 110);
    const auto *msh1_111 = buffer.data(msh1 + 111);
    const auto *msh1_112 = buffer.data(msh1 + 112);
    const auto *msh1_113 = buffer.data(msh1 + 113);
    const auto *msh1_114 = buffer.data(msh1 + 114);
    const auto *msh1_119 = buffer.data(msh1 + 119);
    const auto *msh1_120 = buffer.data(msh1 + 120);
    const auto *msh1_121 = buffer.data(msh1 + 121);
    const auto *msh1_122 = buffer.data(msh1 + 122);
    const auto *msh1_123 = buffer.data(msh1 + 123);
    const auto *msh1_124 = buffer.data(msh1 + 124);
    const auto *msh1_125 = buffer.data(msh1 + 125);
    const auto *msh1_126 = buffer.data(msh1 + 126);
    const auto *msh1_128 = buffer.data(msh1 + 128);
    const auto *msh1_129 = buffer.data(msh1 + 129);
    const auto *msh1_131 = buffer.data(msh1 + 131);
    const auto *msh1_132 = buffer.data(msh1 + 132);
    const auto *msh1_133 = buffer.data(msh1 + 133);
    const auto *msh1_135 = buffer.data(msh1 + 135);
    const auto *msh1_136 = buffer.data(msh1 + 136);
    const auto *msh1_141 = buffer.data(msh1 + 141);
    const auto *msh1_142 = buffer.data(msh1 + 142);
    const auto *msh1_143 = buffer.data(msh1 + 143);
    const auto *msh1_144 = buffer.data(msh1 + 144);

    const auto *msi_99 = buffer.data(msi + 99);
    const auto *msi_105 = buffer.data(msi + 105);
    const auto *msi_106 = buffer.data(msi + 106);
    const auto *msi_107 = buffer.data(msi + 107);
    const auto *msi_108 = buffer.data(msi + 108);
    const auto *msi_109 = buffer.data(msi + 109);
    const auto *msi_110 = buffer.data(msi + 110);
    const auto *msi_111 = buffer.data(msi + 111);
    const auto *msi_112 = buffer.data(msi + 112);
    const auto *msi_114 = buffer.data(msi + 114);
    const auto *msi_115 = buffer.data(msi + 115);
    const auto *msi_117 = buffer.data(msi + 117);
    const auto *msi_118 = buffer.data(msi + 118);
    const auto *msi_121 = buffer.data(msi + 121);
    const auto *msi_122 = buffer.data(msi + 122);
    const auto *msi_126 = buffer.data(msi + 126);
    const auto *msi_133 = buffer.data(msi + 133);
    const auto *msi_134 = buffer.data(msi + 134);
    const auto *msi_135 = buffer.data(msi + 135);
    const auto *msi_136 = buffer.data(msi + 136);
    const auto *msi_137 = buffer.data(msi + 137);
    const auto *msi_138 = buffer.data(msi + 138);
    const auto *msi_139 = buffer.data(msi + 139);
    const auto *msi_140 = buffer.data(msi + 140);
    const auto *msi_141 = buffer.data(msi + 141);
    const auto *msi_142 = buffer.data(msi + 142);
    const auto *msi_143 = buffer.data(msi + 143);
    const auto *msi_144 = buffer.data(msi + 144);
    const auto *msi_145 = buffer.data(msi + 145);
    const auto *msi_146 = buffer.data(msi + 146);
    const auto *msi_147 = buffer.data(msi + 147);
    const auto *msi_148 = buffer.data(msi + 148);
    const auto *msi_149 = buffer.data(msi + 149);
    const auto *msi_150 = buffer.data(msi + 150);
    const auto *msi_151 = buffer.data(msi + 151);
    const auto *msi_152 = buffer.data(msi + 152);
    const auto *msi_153 = buffer.data(msi + 153);
    const auto *msi_154 = buffer.data(msi + 154);
    const auto *msi_160 = buffer.data(msi + 160);
    const auto *msi_161 = buffer.data(msi + 161);
    const auto *msi_162 = buffer.data(msi + 162);
    const auto *msi_163 = buffer.data(msi + 163);
    const auto *msi_164 = buffer.data(msi + 164);
    const auto *msi_165 = buffer.data(msi + 165);
    const auto *msi_166 = buffer.data(msi + 166);
    const auto *msi_167 = buffer.data(msi + 167);
    const auto *msi_168 = buffer.data(msi + 168);
    const auto *msi_169 = buffer.data(msi + 169);
    const auto *msi_170 = buffer.data(msi + 170);
    const auto *msi_171 = buffer.data(msi + 171);
    const auto *msi_173 = buffer.data(msi + 173);
    const auto *msi_174 = buffer.data(msi + 174);
    const auto *msi_175 = buffer.data(msi + 175);
    const auto *msi_177 = buffer.data(msi + 177);
    const auto *msi_178 = buffer.data(msi + 178);
    const auto *msi_179 = buffer.data(msi + 179);
    const auto *msi_180 = buffer.data(msi + 180);
    const auto *msi_182 = buffer.data(msi + 182);
    const auto *msi_183 = buffer.data(msi + 183);
    const auto *msi_189 = buffer.data(msi + 189);
    const auto *msi_190 = buffer.data(msi + 190);
    const auto *msi_191 = buffer.data(msi + 191);
    const auto *msi_192 = buffer.data(msi + 192);
    const auto *msi_193 = buffer.data(msi + 193);
    const auto *msi_194 = buffer.data(msi + 194);
    const auto *msi_195 = buffer.data(msi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, lsi_105, lsi_107, \
                         lsi_108, lsi_109, msi_99, msi_105, msi_107, msi_108, \
                         msi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_21 * lsi_105[k]
                   + f_3 * pc_x[k] * msi_105[k];

        t_130[k] = f_3 * pc_z[k] * msi_99[k];

        t_131[k] = f_21 * lsi_107[k]
                   + f_3 * pc_x[k] * msi_107[k];

        t_132[k] = f_21 * lsi_108[k]
                   + f_3 * pc_x[k] * msi_108[k];

        t_133[k] = f_21 * lsi_109[k]
                   + f_3 * pc_x[k] * msi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, lsi_49, lsi_110, \
                         lsi_111, msh0_78, msh1_78, msi_105, msi_110, \
                         msi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_21 * lsi_110[k]
                   + f_3 * pc_x[k] * msi_110[k];

        t_135[k] = f_21 * lsi_111[k]
                   + f_3 * pc_x[k] * msi_111[k];

        t_136[k] = f_14 * lsi_49[k]
                   + f_1 * msh0_78[k]
                   - f_2 * msh1_78[k]
                   + f_3 * pc_y[k] * msi_105[k];

        t_137[k] = f_3 * pc_z[k] * msi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, msh0_78, msh0_79, msh0_80, msh1_78, \
                         msh1_79, msh1_80, msi_106, msi_107, msi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * msh0_78[k]
                   - f_5 * msh1_78[k]
                   + f_3 * pc_z[k] * msi_106[k];

        t_139[k] = f_6 * msh0_79[k]
                   - f_7 * msh1_79[k]
                   + f_3 * pc_z[k] * msi_107[k];

        t_140[k] = f_8 * msh0_80[k]
                   - f_9 * msh1_80[k]
                   + f_3 * pc_z[k] * msi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, lsk0_72, lsi_55, \
                         lsk1_72, msh0_81, msh0_83, msh1_81, msh1_83, msi_109, \
                         msi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * msh0_81[k]
                   - f_11 * msh1_81[k]
                   + f_3 * pc_z[k] * msi_109[k];

        t_142[k] = f_14 * lsi_55[k]
                   + f_3 * pc_y[k] * msi_111[k];

        t_143[k] = f_1 * msh0_83[k]
                   - f_2 * msh1_83[k]
                   + f_3 * pc_z[k] * msi_111[k];

        t_144[k] = pa_y[k] * lsk0_72[k]
                   - f_12 * pc_y[k] * lsk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, lsk0_39, lsi_28, \
                         lsi_56, lsi_58, lsk1_39, msi_112, msi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * lsi_56[k]
                   + f_3 * pc_y[k] * msi_112[k];

        t_146[k] = f_13 * lsi_28[k]
                   + f_3 * pc_z[k] * msi_112[k];

        t_147[k] = pa_z[k] * lsk0_39[k]
                   - f_12 * pc_z[k] * lsk1_39[k];

        t_148[k] = f_13 * lsi_58[k]
                   + f_3 * pc_y[k] * msi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, lsk0_42, lsk0_77, \
                         lsi_31, lsi_61, lsk1_42, lsk1_77, msi_115, \
                         msi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * lsk0_77[k]
                   - f_12 * pc_y[k] * lsk1_77[k];

        t_150[k] = pa_z[k] * lsk0_42[k]
                   - f_12 * pc_z[k] * lsk1_42[k];

        t_151[k] = f_13 * lsi_31[k]
                   + f_3 * pc_z[k] * msi_115[k];

        t_152[k] = f_13 * lsi_61[k]
                   + f_3 * pc_y[k] * msi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, lsk0_46, lsk0_81, \
                         lsi_34, lsk1_46, lsk1_81, msi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * lsk0_81[k]
                   - f_12 * pc_y[k] * lsk1_81[k];

        t_154[k] = pa_z[k] * lsk0_46[k]
                   - f_12 * pc_z[k] * lsk1_46[k];

        t_155[k] = f_13 * lsi_34[k]
                   + f_3 * pc_z[k] * msi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, lsk0_84, lsk0_86, lsi_64, lsi_65, \
                         lsk1_84, lsk1_86, msi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * lsk0_84[k]
                   + f_14 * lsi_64[k]
                   - f_12 * pc_y[k] * lsk1_84[k];

        t_157[k] = f_13 * lsi_65[k]
                   + f_3 * pc_y[k] * msi_121[k];

        t_158[k] = pa_y[k] * lsk0_86[k]
                   - f_12 * pc_y[k] * lsk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, lsk0_51, lsk0_89, \
                         lsi_38, lsi_68, lsk1_51, lsk1_89, msi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * lsk0_51[k]
                   - f_12 * pc_z[k] * lsk1_51[k];

        t_160[k] = f_13 * lsi_38[k]
                   + f_3 * pc_z[k] * msi_122[k];

        t_161[k] = pa_y[k] * lsk0_89[k]
                   + f_15 * lsi_68[k]
                   - f_12 * pc_y[k] * lsk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, lsk0_90, lsk0_92, \
                         lsi_69, lsi_70, lsi_133, lsk1_90, lsk1_92, msi_126, \
                         msi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * lsk0_90[k]
                   + f_14 * lsi_69[k]
                   - f_12 * pc_y[k] * lsk1_90[k];

        t_163[k] = f_13 * lsi_70[k]
                   + f_3 * pc_y[k] * msi_126[k];

        t_164[k] = pa_y[k] * lsk0_92[k]
                   - f_12 * pc_y[k] * lsk1_92[k];

        t_165[k] = f_21 * lsi_133[k]
                   + f_3 * pc_x[k] * msi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, lsi_134, lsi_135, lsi_136, \
                         lsi_137, lsi_138, msi_134, msi_135, msi_136, msi_137, \
                         msi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_21 * lsi_134[k]
                   + f_3 * pc_x[k] * msi_134[k];

        t_167[k] = f_21 * lsi_135[k]
                   + f_3 * pc_x[k] * msi_135[k];

        t_168[k] = f_21 * lsi_136[k]
                   + f_3 * pc_x[k] * msi_136[k];

        t_169[k] = f_21 * lsi_137[k]
                   + f_3 * pc_x[k] * msi_137[k];

        t_170[k] = f_21 * lsi_138[k]
                   + f_3 * pc_x[k] * msi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, lsk0_64, lsi_49, lsi_139, \
                         lsk1_64, msi_133, msi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_21 * lsi_139[k]
                   + f_3 * pc_x[k] * msi_139[k];

        t_172[k] = pa_z[k] * lsk0_64[k]
                   - f_12 * pc_z[k] * lsk1_64[k];

        t_173[k] = f_13 * lsi_49[k]
                   + f_3 * pc_z[k] * msi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, lsi_79, lsi_80, lsi_81, msh0_101, \
                         msh0_102, msh0_103, msh1_101, msh1_102, msh1_103, msi_135, msi_136, \
                         msi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * lsi_79[k]
                   + f_10 * msh0_101[k]
                   - f_11 * msh1_101[k]
                   + f_3 * pc_y[k] * msi_135[k];

        t_175[k] = f_13 * lsi_80[k]
                   + f_8 * msh0_102[k]
                   - f_9 * msh1_102[k]
                   + f_3 * pc_y[k] * msi_136[k];

        t_176[k] = f_13 * lsi_81[k]
                   + f_6 * msh0_103[k]
                   - f_7 * msh1_103[k]
                   + f_3 * pc_y[k] * msi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, lsk0_107, lsi_82, lsi_83, lsk1_107, \
                         msh0_104, msh1_104, msi_138, msi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * lsi_82[k]
                   + f_4 * msh0_104[k]
                   - f_5 * msh1_104[k]
                   + f_3 * pc_y[k] * msi_138[k];

        t_178[k] = f_13 * lsi_83[k]
                   + f_3 * pc_y[k] * msi_139[k];

        t_179[k] = pa_y[k] * lsk0_107[k]
                   - f_12 * pc_y[k] * lsk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, lsi_56, lsi_140, \
                         msh0_105, msh1_105, msi_140, msi_141, \
                         msi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_21 * lsi_140[k]
                   + f_1 * msh0_105[k]
                   - f_2 * msh1_105[k]
                   + f_3 * pc_x[k] * msi_140[k];

        t_181[k] = f_3 * pc_y[k] * msi_140[k];

        t_182[k] = f_14 * lsi_56[k]
                   + f_3 * pc_z[k] * msi_140[k];

        t_183[k] = f_4 * msh0_105[k]
                   - f_5 * msh1_105[k]
                   + f_3 * pc_y[k] * msi_141[k];

        t_184[k] = f_3 * pc_y[k] * msi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, lsi_145, msh0_106, msh0_107, \
                         msh0_110, msh1_106, msh1_107, msh1_110, msi_143, msi_144, \
                         msi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_21 * lsi_145[k]
                   + f_10 * msh0_110[k]
                   - f_11 * msh1_110[k]
                   + f_3 * pc_x[k] * msi_145[k];

        t_186[k] = f_6 * msh0_106[k]
                   - f_7 * msh1_106[k]
                   + f_3 * pc_y[k] * msi_143[k];

        t_187[k] = f_4 * msh0_107[k]
                   - f_5 * msh1_107[k]
                   + f_3 * pc_y[k] * msi_144[k];

        t_188[k] = f_3 * pc_y[k] * msi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, lsi_149, msh0_108, msh0_109, \
                         msh0_114, msh1_108, msh1_109, msh1_114, msi_146, msi_147, \
                         msi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_21 * lsi_149[k]
                   + f_8 * msh0_114[k]
                   - f_9 * msh1_114[k]
                   + f_3 * pc_x[k] * msi_149[k];

        t_190[k] = f_8 * msh0_108[k]
                   - f_9 * msh1_108[k]
                   + f_3 * pc_y[k] * msi_146[k];

        t_191[k] = f_6 * msh0_109[k]
                   - f_7 * msh1_109[k]
                   + f_3 * pc_y[k] * msi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, lsi_154, msh0_110, msh0_119, \
                         msh1_110, msh1_119, msi_148, msi_149, \
                         msi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * msh0_110[k]
                   - f_5 * msh1_110[k]
                   + f_3 * pc_y[k] * msi_148[k];

        t_193[k] = f_3 * pc_y[k] * msi_149[k];

        t_194[k] = f_21 * lsi_154[k]
                   + f_6 * msh0_119[k]
                   - f_7 * msh1_119[k]
                   + f_3 * pc_x[k] * msi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, msh0_111, msh0_112, msh0_113, msh1_111, \
                         msh1_112, msh1_113, msi_150, msi_151, \
                         msi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * msh0_111[k]
                   - f_11 * msh1_111[k]
                   + f_3 * pc_y[k] * msi_150[k];

        t_196[k] = f_8 * msh0_112[k]
                   - f_9 * msh1_112[k]
                   + f_3 * pc_y[k] * msi_151[k];

        t_197[k] = f_6 * msh0_113[k]
                   - f_7 * msh1_113[k]
                   + f_3 * pc_y[k] * msi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, lsi_160, lsi_161, msh0_114, \
                         msh0_125, msh1_114, msh1_125, msi_153, msi_154, msi_160, \
                         msi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * msh0_114[k]
                   - f_5 * msh1_114[k]
                   + f_3 * pc_y[k] * msi_153[k];

        t_199[k] = f_3 * pc_y[k] * msi_154[k];

        t_200[k] = f_21 * lsi_160[k]
                   + f_4 * msh0_125[k]
                   - f_5 * msh1_125[k]
                   + f_3 * pc_x[k] * msi_160[k];

        t_201[k] = f_21 * lsi_161[k]
                   + f_3 * pc_x[k] * msi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, lsi_162, lsi_163, \
                         lsi_164, lsi_165, msi_160, msi_162, msi_163, msi_164, \
                         msi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_21 * lsi_162[k]
                   + f_3 * pc_x[k] * msi_162[k];

        t_203[k] = f_21 * lsi_163[k]
                   + f_3 * pc_x[k] * msi_163[k];

        t_204[k] = f_21 * lsi_164[k]
                   + f_3 * pc_x[k] * msi_164[k];

        t_205[k] = f_21 * lsi_165[k]
                   + f_3 * pc_x[k] * msi_165[k];

        t_206[k] = f_3 * pc_y[k] * msi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, lsi_167, msh0_120, msh0_121, \
                         msh1_120, msh1_121, msi_161, msi_162, \
                         msi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_21 * lsi_167[k]
                   + f_3 * pc_x[k] * msi_167[k];

        t_208[k] = f_1 * msh0_120[k]
                   - f_2 * msh1_120[k]
                   + f_3 * pc_y[k] * msi_161[k];

        t_209[k] = f_19 * msh0_121[k]
                   - f_20 * msh1_121[k]
                   + f_3 * pc_y[k] * msi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, msh0_122, msh0_123, msh0_124, msh1_122, \
                         msh1_123, msh1_124, msi_163, msi_164, \
                         msi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * msh0_122[k]
                   - f_11 * msh1_122[k]
                   + f_3 * pc_y[k] * msi_163[k];

        t_211[k] = f_8 * msh0_123[k]
                   - f_9 * msh1_123[k]
                   + f_3 * pc_y[k] * msi_164[k];

        t_212[k] = f_6 * msh0_124[k]
                   - f_7 * msh1_124[k]
                   + f_3 * pc_y[k] * msi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, lsi_83, lsi_168, \
                         msh0_125, msh0_126, msh1_125, msh1_126, msi_166, msi_167, \
                         msi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * msh0_125[k]
                   - f_5 * msh1_125[k]
                   + f_3 * pc_y[k] * msi_166[k];

        t_214[k] = f_3 * pc_y[k] * msi_167[k];

        t_215[k] = f_14 * lsi_83[k]
                   + f_1 * msh0_125[k]
                   - f_2 * msh1_125[k]
                   + f_3 * pc_z[k] * msi_167[k];

        t_216[k] = f_22 * lsi_168[k]
                   + f_1 * msh0_126[k]
                   - f_2 * msh1_126[k]
                   + f_3 * pc_x[k] * msi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, lsi_84, lsi_171, \
                         msh0_129, msh1_129, msi_168, msi_169, \
                         msi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * lsi_84[k]
                   + f_3 * pc_y[k] * msi_168[k];

        t_218[k] = f_3 * pc_z[k] * msi_168[k];

        t_219[k] = f_22 * lsi_171[k]
                   + f_10 * msh0_129[k]
                   - f_11 * msh1_129[k]
                   + f_3 * pc_x[k] * msi_171[k];

        t_220[k] = f_3 * pc_z[k] * msi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pc_x, pc_z, lsi_174, msh0_126, msh0_132, \
                         msh1_126, msh1_132, msi_170, msi_171, \
                         msi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * msh0_126[k]
                   - f_5 * msh1_126[k]
                   + f_3 * pc_z[k] * msi_170[k];

        t_222[k] = f_22 * lsi_174[k]
                   + f_8 * msh0_132[k]
                   - f_9 * msh1_132[k]
                   + f_3 * pc_x[k] * msi_174[k];

        t_223[k] = f_3 * pc_z[k] * msi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, lsi_89, lsi_178, \
                         msh0_128, msh0_136, msh1_128, msh1_136, msi_173, msi_174, \
                         msi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * lsi_89[k]
                   + f_3 * pc_y[k] * msi_173[k];

        t_225[k] = f_6 * msh0_128[k]
                   - f_7 * msh1_128[k]
                   + f_3 * pc_z[k] * msi_173[k];

        t_226[k] = f_22 * lsi_178[k]
                   + f_6 * msh0_136[k]
                   - f_7 * msh1_136[k]
                   + f_3 * pc_x[k] * msi_178[k];

        t_227[k] = f_3 * pc_z[k] * msi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, lsi_93, msh0_129, msh0_131, \
                         msh1_129, msh1_131, msi_175, msi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * msh0_129[k]
                   - f_5 * msh1_129[k]
                   + f_3 * pc_z[k] * msi_175[k];

        t_229[k] = f_15 * lsi_93[k]
                   + f_3 * pc_y[k] * msi_177[k];

        t_230[k] = f_8 * msh0_131[k]
                   - f_9 * msh1_131[k]
                   + f_3 * pc_z[k] * msi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_x, pc_z, lsi_183, msh0_132, msh0_141, \
                         msh1_132, msh1_141, msi_178, msi_179, \
                         msi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_22 * lsi_183[k]
                   + f_4 * msh0_141[k]
                   - f_5 * msh1_141[k]
                   + f_3 * pc_x[k] * msi_183[k];

        t_232[k] = f_3 * pc_z[k] * msi_178[k];

        t_233[k] = f_4 * msh0_132[k]
                   - f_5 * msh1_132[k]
                   + f_3 * pc_z[k] * msi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, lsi_98, lsi_189, \
                         msh0_133, msh0_135, msh1_133, msh1_135, msi_180, msi_182, \
                         msi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * msh0_133[k]
                   - f_7 * msh1_133[k]
                   + f_3 * pc_z[k] * msi_180[k];

        t_235[k] = f_15 * lsi_98[k]
                   + f_3 * pc_y[k] * msi_182[k];

        t_236[k] = f_10 * msh0_135[k]
                   - f_11 * msh1_135[k]
                   + f_3 * pc_z[k] * msi_182[k];

        t_237[k] = f_22 * lsi_189[k]
                   + f_3 * pc_x[k] * msi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, lsi_191, lsi_192, \
                         lsi_193, lsi_194, msi_183, msi_191, msi_192, msi_193, \
                         msi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * msi_183[k];

        t_239[k] = f_22 * lsi_191[k]
                   + f_3 * pc_x[k] * msi_191[k];

        t_240[k] = f_22 * lsi_192[k]
                   + f_3 * pc_x[k] * msi_192[k];

        t_241[k] = f_22 * lsi_193[k]
                   + f_3 * pc_x[k] * msi_193[k];

        t_242[k] = f_22 * lsi_194[k]
                   + f_3 * pc_x[k] * msi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, lsi_105, lsi_195, \
                         msh0_141, msh1_141, msi_189, msi_190, \
                         msi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_22 * lsi_195[k]
                   + f_3 * pc_x[k] * msi_195[k];

        t_244[k] = f_15 * lsi_105[k]
                   + f_1 * msh0_141[k]
                   - f_2 * msh1_141[k]
                   + f_3 * pc_y[k] * msi_189[k];

        t_245[k] = f_3 * pc_z[k] * msi_189[k];

        t_246[k] = f_4 * msh0_141[k]
                   - f_5 * msh1_141[k]
                   + f_3 * pc_z[k] * msi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, msh0_142, msh0_143, msh0_144, msh1_142, \
                         msh1_143, msh1_144, msi_191, msi_192, \
                         msi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * msh0_142[k]
                   - f_7 * msh1_142[k]
                   + f_3 * pc_z[k] * msi_191[k];

        t_248[k] = f_8 * msh0_143[k]
                   - f_9 * msh1_143[k]
                   + f_3 * pc_z[k] * msi_192[k];

        t_249[k] = f_10 * msh0_144[k]
                   - f_11 * msh1_144[k]
                   + f_3 * pc_z[k] * msi_193[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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
    const auto f_22 = 3.0 / q;

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

    const auto *lsk0_108 = buffer.data(lsk0 + 108);
    const auto *lsk0_111 = buffer.data(lsk0 + 111);
    const auto *lsk0_114 = buffer.data(lsk0 + 114);
    const auto *lsk0_118 = buffer.data(lsk0 + 118);
    const auto *lsk0_120 = buffer.data(lsk0 + 120);
    const auto *lsk0_123 = buffer.data(lsk0 + 123);
    const auto *lsk0_125 = buffer.data(lsk0 + 125);
    const auto *lsk0_126 = buffer.data(lsk0 + 126);
    const auto *lsk0_136 = buffer.data(lsk0 + 136);
    const auto *lsk0_180 = buffer.data(lsk0 + 180);
    const auto *lsk0_183 = buffer.data(lsk0 + 183);
    const auto *lsk0_185 = buffer.data(lsk0 + 185);
    const auto *lsk0_186 = buffer.data(lsk0 + 186);
    const auto *lsk0_189 = buffer.data(lsk0 + 189);
    const auto *lsk0_190 = buffer.data(lsk0 + 190);
    const auto *lsk0_192 = buffer.data(lsk0 + 192);
    const auto *lsk0_194 = buffer.data(lsk0 + 194);
    const auto *lsk0_195 = buffer.data(lsk0 + 195);
    const auto *lsk0_197 = buffer.data(lsk0 + 197);
    const auto *lsk0_198 = buffer.data(lsk0 + 198);
    const auto *lsk0_200 = buffer.data(lsk0 + 200);
    const auto *lsk0_215 = buffer.data(lsk0 + 215);

    const auto *lsi_84 = buffer.data(lsi + 84);
    const auto *lsi_87 = buffer.data(lsi + 87);
    const auto *lsi_90 = buffer.data(lsi + 90);
    const auto *lsi_91 = buffer.data(lsi + 91);
    const auto *lsi_94 = buffer.data(lsi + 94);
    const auto *lsi_95 = buffer.data(lsi + 95);
    const auto *lsi_96 = buffer.data(lsi + 96);
    const auto *lsi_105 = buffer.data(lsi + 105);
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
    const auto *lsi_135 = buffer.data(lsi + 135);
    const auto *lsi_136 = buffer.data(lsi + 136);
    const auto *lsi_137 = buffer.data(lsi + 137);
    const auto *lsi_138 = buffer.data(lsi + 138);
    const auto *lsi_139 = buffer.data(lsi + 139);
    const auto *lsi_140 = buffer.data(lsi + 140);
    const auto *lsi_141 = buffer.data(lsi + 141);
    const auto *lsi_142 = buffer.data(lsi + 142);
    const auto *lsi_143 = buffer.data(lsi + 143);
    const auto *lsi_145 = buffer.data(lsi + 145);
    const auto *lsi_146 = buffer.data(lsi + 146);
    const auto *lsi_148 = buffer.data(lsi + 148);
    const auto *lsi_149 = buffer.data(lsi + 149);
    const auto *lsi_150 = buffer.data(lsi + 150);
    const auto *lsi_152 = buffer.data(lsi + 152);
    const auto *lsi_153 = buffer.data(lsi + 153);
    const auto *lsi_154 = buffer.data(lsi + 154);
    const auto *lsi_161 = buffer.data(lsi + 161);
    const auto *lsi_163 = buffer.data(lsi + 163);
    const auto *lsi_164 = buffer.data(lsi + 164);
    const auto *lsi_165 = buffer.data(lsi + 165);
    const auto *lsi_166 = buffer.data(lsi + 166);
    const auto *lsi_167 = buffer.data(lsi + 167);
    const auto *lsi_168 = buffer.data(lsi + 168);
    const auto *lsi_201 = buffer.data(lsi + 201);
    const auto *lsi_205 = buffer.data(lsi + 205);
    const auto *lsi_210 = buffer.data(lsi + 210);
    const auto *lsi_216 = buffer.data(lsi + 216);
    const auto *lsi_217 = buffer.data(lsi + 217);
    const auto *lsi_218 = buffer.data(lsi + 218);
    const auto *lsi_219 = buffer.data(lsi + 219);
    const auto *lsi_220 = buffer.data(lsi + 220);
    const auto *lsi_221 = buffer.data(lsi + 221);
    const auto *lsi_222 = buffer.data(lsi + 222);
    const auto *lsi_223 = buffer.data(lsi + 223);
    const auto *lsi_245 = buffer.data(lsi + 245);
    const auto *lsi_246 = buffer.data(lsi + 246);
    const auto *lsi_247 = buffer.data(lsi + 247);
    const auto *lsi_248 = buffer.data(lsi + 248);
    const auto *lsi_249 = buffer.data(lsi + 249);
    const auto *lsi_250 = buffer.data(lsi + 250);
    const auto *lsi_251 = buffer.data(lsi + 251);
    const auto *lsi_252 = buffer.data(lsi + 252);
    const auto *lsi_257 = buffer.data(lsi + 257);
    const auto *lsi_261 = buffer.data(lsi + 261);
    const auto *lsi_266 = buffer.data(lsi + 266);
    const auto *lsi_272 = buffer.data(lsi + 272);
    const auto *lsi_273 = buffer.data(lsi + 273);
    const auto *lsi_274 = buffer.data(lsi + 274);
    const auto *lsi_275 = buffer.data(lsi + 275);
    const auto *lsi_276 = buffer.data(lsi + 276);
    const auto *lsi_277 = buffer.data(lsi + 277);
    const auto *lsi_279 = buffer.data(lsi + 279);
    const auto *lsi_280 = buffer.data(lsi + 280);
    const auto *lsi_283 = buffer.data(lsi + 283);
    const auto *lsi_286 = buffer.data(lsi + 286);

    const auto *lsk1_108 = buffer.data(lsk1 + 108);
    const auto *lsk1_111 = buffer.data(lsk1 + 111);
    const auto *lsk1_114 = buffer.data(lsk1 + 114);
    const auto *lsk1_118 = buffer.data(lsk1 + 118);
    const auto *lsk1_120 = buffer.data(lsk1 + 120);
    const auto *lsk1_123 = buffer.data(lsk1 + 123);
    const auto *lsk1_125 = buffer.data(lsk1 + 125);
    const auto *lsk1_126 = buffer.data(lsk1 + 126);
    const auto *lsk1_136 = buffer.data(lsk1 + 136);
    const auto *lsk1_180 = buffer.data(lsk1 + 180);
    const auto *lsk1_183 = buffer.data(lsk1 + 183);
    const auto *lsk1_185 = buffer.data(lsk1 + 185);
    const auto *lsk1_186 = buffer.data(lsk1 + 186);
    const auto *lsk1_189 = buffer.data(lsk1 + 189);
    const auto *lsk1_190 = buffer.data(lsk1 + 190);
    const auto *lsk1_192 = buffer.data(lsk1 + 192);
    const auto *lsk1_194 = buffer.data(lsk1 + 194);
    const auto *lsk1_195 = buffer.data(lsk1 + 195);
    const auto *lsk1_197 = buffer.data(lsk1 + 197);
    const auto *lsk1_198 = buffer.data(lsk1 + 198);
    const auto *lsk1_200 = buffer.data(lsk1 + 200);
    const auto *lsk1_215 = buffer.data(lsk1 + 215);

    const auto *msh0_146 = buffer.data(msh0 + 146);
    const auto *msh0_152 = buffer.data(msh0 + 152);
    const auto *msh0_156 = buffer.data(msh0 + 156);
    const auto *msh0_161 = buffer.data(msh0 + 161);
    const auto *msh0_164 = buffer.data(msh0 + 164);
    const auto *msh0_165 = buffer.data(msh0 + 165);
    const auto *msh0_166 = buffer.data(msh0 + 166);
    const auto *msh0_167 = buffer.data(msh0 + 167);
    const auto *msh0_183 = buffer.data(msh0 + 183);
    const auto *msh0_185 = buffer.data(msh0 + 185);
    const auto *msh0_186 = buffer.data(msh0 + 186);
    const auto *msh0_187 = buffer.data(msh0 + 187);
    const auto *msh0_188 = buffer.data(msh0 + 188);
    const auto *msh0_189 = buffer.data(msh0 + 189);
    const auto *msh0_190 = buffer.data(msh0 + 190);
    const auto *msh0_191 = buffer.data(msh0 + 191);
    const auto *msh0_192 = buffer.data(msh0 + 192);
    const auto *msh0_193 = buffer.data(msh0 + 193);
    const auto *msh0_194 = buffer.data(msh0 + 194);
    const auto *msh0_195 = buffer.data(msh0 + 195);
    const auto *msh0_196 = buffer.data(msh0 + 196);
    const auto *msh0_197 = buffer.data(msh0 + 197);
    const auto *msh0_198 = buffer.data(msh0 + 198);
    const auto *msh0_203 = buffer.data(msh0 + 203);
    const auto *msh0_204 = buffer.data(msh0 + 204);
    const auto *msh0_205 = buffer.data(msh0 + 205);
    const auto *msh0_206 = buffer.data(msh0 + 206);
    const auto *msh0_207 = buffer.data(msh0 + 207);
    const auto *msh0_208 = buffer.data(msh0 + 208);
    const auto *msh0_209 = buffer.data(msh0 + 209);
    const auto *msh0_210 = buffer.data(msh0 + 210);
    const auto *msh0_213 = buffer.data(msh0 + 213);
    const auto *msh0_216 = buffer.data(msh0 + 216);

    const auto *msh1_146 = buffer.data(msh1 + 146);
    const auto *msh1_152 = buffer.data(msh1 + 152);
    const auto *msh1_156 = buffer.data(msh1 + 156);
    const auto *msh1_161 = buffer.data(msh1 + 161);
    const auto *msh1_164 = buffer.data(msh1 + 164);
    const auto *msh1_165 = buffer.data(msh1 + 165);
    const auto *msh1_166 = buffer.data(msh1 + 166);
    const auto *msh1_167 = buffer.data(msh1 + 167);
    const auto *msh1_183 = buffer.data(msh1 + 183);
    const auto *msh1_185 = buffer.data(msh1 + 185);
    const auto *msh1_186 = buffer.data(msh1 + 186);
    const auto *msh1_187 = buffer.data(msh1 + 187);
    const auto *msh1_188 = buffer.data(msh1 + 188);
    const auto *msh1_189 = buffer.data(msh1 + 189);
    const auto *msh1_190 = buffer.data(msh1 + 190);
    const auto *msh1_191 = buffer.data(msh1 + 191);
    const auto *msh1_192 = buffer.data(msh1 + 192);
    const auto *msh1_193 = buffer.data(msh1 + 193);
    const auto *msh1_194 = buffer.data(msh1 + 194);
    const auto *msh1_195 = buffer.data(msh1 + 195);
    const auto *msh1_196 = buffer.data(msh1 + 196);
    const auto *msh1_197 = buffer.data(msh1 + 197);
    const auto *msh1_198 = buffer.data(msh1 + 198);
    const auto *msh1_203 = buffer.data(msh1 + 203);
    const auto *msh1_204 = buffer.data(msh1 + 204);
    const auto *msh1_205 = buffer.data(msh1 + 205);
    const auto *msh1_206 = buffer.data(msh1 + 206);
    const auto *msh1_207 = buffer.data(msh1 + 207);
    const auto *msh1_208 = buffer.data(msh1 + 208);
    const auto *msh1_209 = buffer.data(msh1 + 209);
    const auto *msh1_210 = buffer.data(msh1 + 210);
    const auto *msh1_213 = buffer.data(msh1 + 213);
    const auto *msh1_216 = buffer.data(msh1 + 216);

    const auto *msi_195 = buffer.data(msi + 195);
    const auto *msi_196 = buffer.data(msi + 196);
    const auto *msi_198 = buffer.data(msi + 198);
    const auto *msi_199 = buffer.data(msi + 199);
    const auto *msi_201 = buffer.data(msi + 201);
    const auto *msi_202 = buffer.data(msi + 202);
    const auto *msi_205 = buffer.data(msi + 205);
    const auto *msi_206 = buffer.data(msi + 206);
    const auto *msi_210 = buffer.data(msi + 210);
    const auto *msi_216 = buffer.data(msi + 216);
    const auto *msi_217 = buffer.data(msi + 217);
    const auto *msi_218 = buffer.data(msi + 218);
    const auto *msi_219 = buffer.data(msi + 219);
    const auto *msi_220 = buffer.data(msi + 220);
    const auto *msi_221 = buffer.data(msi + 221);
    const auto *msi_222 = buffer.data(msi + 222);
    const auto *msi_223 = buffer.data(msi + 223);
    const auto *msi_224 = buffer.data(msi + 224);
    const auto *msi_226 = buffer.data(msi + 226);
    const auto *msi_227 = buffer.data(msi + 227);
    const auto *msi_229 = buffer.data(msi + 229);
    const auto *msi_230 = buffer.data(msi + 230);
    const auto *msi_233 = buffer.data(msi + 233);
    const auto *msi_234 = buffer.data(msi + 234);
    const auto *msi_238 = buffer.data(msi + 238);
    const auto *msi_245 = buffer.data(msi + 245);
    const auto *msi_246 = buffer.data(msi + 246);
    const auto *msi_247 = buffer.data(msi + 247);
    const auto *msi_248 = buffer.data(msi + 248);
    const auto *msi_249 = buffer.data(msi + 249);
    const auto *msi_250 = buffer.data(msi + 250);
    const auto *msi_251 = buffer.data(msi + 251);
    const auto *msi_252 = buffer.data(msi + 252);
    const auto *msi_253 = buffer.data(msi + 253);
    const auto *msi_254 = buffer.data(msi + 254);
    const auto *msi_255 = buffer.data(msi + 255);
    const auto *msi_256 = buffer.data(msi + 256);
    const auto *msi_257 = buffer.data(msi + 257);
    const auto *msi_258 = buffer.data(msi + 258);
    const auto *msi_259 = buffer.data(msi + 259);
    const auto *msi_260 = buffer.data(msi + 260);
    const auto *msi_261 = buffer.data(msi + 261);
    const auto *msi_262 = buffer.data(msi + 262);
    const auto *msi_263 = buffer.data(msi + 263);
    const auto *msi_264 = buffer.data(msi + 264);
    const auto *msi_265 = buffer.data(msi + 265);
    const auto *msi_266 = buffer.data(msi + 266);
    const auto *msi_272 = buffer.data(msi + 272);
    const auto *msi_273 = buffer.data(msi + 273);
    const auto *msi_274 = buffer.data(msi + 274);
    const auto *msi_275 = buffer.data(msi + 275);
    const auto *msi_276 = buffer.data(msi + 276);
    const auto *msi_277 = buffer.data(msi + 277);
    const auto *msi_278 = buffer.data(msi + 278);
    const auto *msi_279 = buffer.data(msi + 279);
    const auto *msi_280 = buffer.data(msi + 280);
    const auto *msi_281 = buffer.data(msi + 281);
    const auto *msi_282 = buffer.data(msi + 282);
    const auto *msi_283 = buffer.data(msi + 283);
    const auto *msi_286 = buffer.data(msi + 286);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, lsk0_108, lsi_111, \
                         lsi_112, lsk1_108, msh0_146, msh1_146, msi_195, \
                         msi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * lsi_111[k]
                   + f_3 * pc_y[k] * msi_195[k];

        t_251[k] = f_1 * msh0_146[k]
                   - f_2 * msh1_146[k]
                   + f_3 * pc_z[k] * msi_195[k];

        t_252[k] = pa_z[k] * lsk0_108[k]
                   - f_12 * pc_z[k] * lsk1_108[k];

        t_253[k] = f_14 * lsi_112[k]
                   + f_3 * pc_y[k] * msi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_y, pc_z, lsk0_111, lsi_84, lsi_114, \
                         lsk1_111, msi_196, msi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * lsi_84[k]
                   + f_3 * pc_z[k] * msi_196[k];

        t_255[k] = pa_z[k] * lsk0_111[k]
                   - f_12 * pc_z[k] * lsk1_111[k];

        t_256[k] = f_14 * lsi_114[k]
                   + f_3 * pc_y[k] * msi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, lsk0_114, lsi_87, lsi_201, \
                         lsk1_114, msh0_152, msh1_152, msi_199, \
                         msi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_22 * lsi_201[k]
                   + f_10 * msh0_152[k]
                   - f_11 * msh1_152[k]
                   + f_3 * pc_x[k] * msi_201[k];

        t_258[k] = pa_z[k] * lsk0_114[k]
                   - f_12 * pc_z[k] * lsk1_114[k];

        t_259[k] = f_13 * lsi_87[k]
                   + f_3 * pc_z[k] * msi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_y, pc_z, lsk0_118, lsi_117, \
                         lsi_205, lsk1_118, msh0_156, msh1_156, msi_201, \
                         msi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * lsi_117[k]
                   + f_3 * pc_y[k] * msi_201[k];

        t_261[k] = f_22 * lsi_205[k]
                   + f_8 * msh0_156[k]
                   - f_9 * msh1_156[k]
                   + f_3 * pc_x[k] * msi_205[k];

        t_262[k] = pa_z[k] * lsk0_118[k]
                   - f_12 * pc_z[k] * lsk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pa_z, pc_y, pc_z, lsk0_120, lsi_90, lsi_91, \
                         lsi_121, lsk1_120, msi_202, msi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * lsi_90[k]
                   + f_3 * pc_z[k] * msi_202[k];

        t_264[k] = pa_z[k] * lsk0_120[k]
                   + f_14 * lsi_91[k]
                   - f_12 * pc_z[k] * lsk1_120[k];

        t_265[k] = f_14 * lsi_121[k]
                   + f_3 * pc_y[k] * msi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, lsk0_123, lsi_94, lsi_210, \
                         lsk1_123, msh0_161, msh1_161, msi_206, \
                         msi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_22 * lsi_210[k]
                   + f_6 * msh0_161[k]
                   - f_7 * msh1_161[k]
                   + f_3 * pc_x[k] * msi_210[k];

        t_267[k] = pa_z[k] * lsk0_123[k]
                   - f_12 * pc_z[k] * lsk1_123[k];

        t_268[k] = f_13 * lsi_94[k]
                   + f_3 * pc_z[k] * msi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_z, pc_y, pc_z, lsk0_125, lsk0_126, lsi_95, \
                         lsi_96, lsi_126, lsk1_125, lsk1_126, msi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_z[k] * lsk0_125[k]
                   + f_14 * lsi_95[k]
                   - f_12 * pc_z[k] * lsk1_125[k];

        t_270[k] = pa_z[k] * lsk0_126[k]
                   + f_15 * lsi_96[k]
                   - f_12 * pc_z[k] * lsk1_126[k];

        t_271[k] = f_14 * lsi_126[k]
                   + f_3 * pc_y[k] * msi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, lsi_216, lsi_217, lsi_218, lsi_219, \
                         msh0_167, msh1_167, msi_216, msi_217, msi_218, \
                         msi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_22 * lsi_216[k]
                   + f_4 * msh0_167[k]
                   - f_5 * msh1_167[k]
                   + f_3 * pc_x[k] * msi_216[k];

        t_273[k] = f_22 * lsi_217[k]
                   + f_3 * pc_x[k] * msi_217[k];

        t_274[k] = f_22 * lsi_218[k]
                   + f_3 * pc_x[k] * msi_218[k];

        t_275[k] = f_22 * lsi_219[k]
                   + f_3 * pc_x[k] * msi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, lsi_220, lsi_221, lsi_222, lsi_223, \
                         msi_220, msi_221, msi_222, msi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_22 * lsi_220[k]
                   + f_3 * pc_x[k] * msi_220[k];

        t_277[k] = f_22 * lsi_221[k]
                   + f_3 * pc_x[k] * msi_221[k];

        t_278[k] = f_22 * lsi_222[k]
                   + f_3 * pc_x[k] * msi_222[k];

        t_279[k] = f_22 * lsi_223[k]
                   + f_3 * pc_x[k] * msi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_z, pc_y, pc_z, lsk0_136, lsi_105, lsi_135, \
                         lsk1_136, msh0_164, msh1_164, msi_217, \
                         msi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * lsk0_136[k]
                   - f_12 * pc_z[k] * lsk1_136[k];

        t_281[k] = f_13 * lsi_105[k]
                   + f_3 * pc_z[k] * msi_217[k];

        t_282[k] = f_14 * lsi_135[k]
                   + f_10 * msh0_164[k]
                   - f_11 * msh1_164[k]
                   + f_3 * pc_y[k] * msi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, lsi_136, lsi_137, lsi_138, msh0_165, \
                         msh0_166, msh0_167, msh1_165, msh1_166, msh1_167, msi_220, msi_221, \
                         msi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * lsi_136[k]
                   + f_8 * msh0_165[k]
                   - f_9 * msh1_165[k]
                   + f_3 * pc_y[k] * msi_220[k];

        t_284[k] = f_14 * lsi_137[k]
                   + f_6 * msh0_166[k]
                   - f_7 * msh1_166[k]
                   + f_3 * pc_y[k] * msi_221[k];

        t_285[k] = f_14 * lsi_138[k]
                   + f_4 * msh0_167[k]
                   - f_5 * msh1_167[k]
                   + f_3 * pc_y[k] * msi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, lsk0_180, lsi_111, \
                         lsi_139, lsi_140, lsk1_180, msh0_167, msh1_167, msi_223, \
                         msi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * lsi_139[k]
                   + f_3 * pc_y[k] * msi_223[k];

        t_287[k] = f_13 * lsi_111[k]
                   + f_1 * msh0_167[k]
                   - f_2 * msh1_167[k]
                   + f_3 * pc_z[k] * msi_223[k];

        t_288[k] = pa_y[k] * lsk0_180[k]
                   - f_12 * pc_y[k] * lsk1_180[k];

        t_289[k] = f_13 * lsi_140[k]
                   + f_3 * pc_y[k] * msi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, pc_z, lsk0_183, lsk0_185, \
                         lsi_112, lsi_141, lsi_142, lsk1_183, lsk1_185, msi_224, \
                         msi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * lsi_112[k]
                   + f_3 * pc_z[k] * msi_224[k];

        t_291[k] = pa_y[k] * lsk0_183[k]
                   + f_14 * lsi_141[k]
                   - f_12 * pc_y[k] * lsk1_183[k];

        t_292[k] = f_13 * lsi_142[k]
                   + f_3 * pc_y[k] * msi_226[k];

        t_293[k] = pa_y[k] * lsk0_185[k]
                   - f_12 * pc_y[k] * lsk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pc_y, pc_z, lsk0_186, lsk0_189, \
                         lsi_115, lsi_143, lsi_145, lsk1_186, lsk1_189, msi_227, \
                         msi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * lsk0_186[k]
                   + f_15 * lsi_143[k]
                   - f_12 * pc_y[k] * lsk1_186[k];

        t_295[k] = f_14 * lsi_115[k]
                   + f_3 * pc_z[k] * msi_227[k];

        t_296[k] = f_13 * lsi_145[k]
                   + f_3 * pc_y[k] * msi_229[k];

        t_297[k] = pa_y[k] * lsk0_189[k]
                   - f_12 * pc_y[k] * lsk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_y, pc_y, pc_z, lsk0_190, lsk0_192, lsi_118, \
                         lsi_146, lsi_148, lsk1_190, lsk1_192, \
                         msi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_y[k] * lsk0_190[k]
                   + f_16 * lsi_146[k]
                   - f_12 * pc_y[k] * lsk1_190[k];

        t_299[k] = f_14 * lsi_118[k]
                   + f_3 * pc_z[k] * msi_230[k];

        t_300[k] = pa_y[k] * lsk0_192[k]
                   + f_14 * lsi_148[k]
                   - f_12 * pc_y[k] * lsk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, pc_y, pc_z, lsk0_194, lsk0_195, \
                         lsi_122, lsi_149, lsi_150, lsk1_194, lsk1_195, msi_233, \
                         msi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * lsi_149[k]
                   + f_3 * pc_y[k] * msi_233[k];

        t_302[k] = pa_y[k] * lsk0_194[k]
                   - f_12 * pc_y[k] * lsk1_194[k];

        t_303[k] = pa_y[k] * lsk0_195[k]
                   + f_17 * lsi_150[k]
                   - f_12 * pc_y[k] * lsk1_195[k];

        t_304[k] = f_14 * lsi_122[k]
                   + f_3 * pc_z[k] * msi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pc_y, lsk0_197, lsk0_198, lsk0_200, \
                         lsi_152, lsi_153, lsi_154, lsk1_197, lsk1_198, lsk1_200, \
                         msi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_y[k] * lsk0_197[k]
                   + f_15 * lsi_152[k]
                   - f_12 * pc_y[k] * lsk1_197[k];

        t_306[k] = pa_y[k] * lsk0_198[k]
                   + f_14 * lsi_153[k]
                   - f_12 * pc_y[k] * lsk1_198[k];

        t_307[k] = f_13 * lsi_154[k]
                   + f_3 * pc_y[k] * msi_238[k];

        t_308[k] = pa_y[k] * lsk0_200[k]
                   - f_12 * pc_y[k] * lsk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, lsi_245, lsi_246, lsi_247, \
                         lsi_248, lsi_249, msi_245, msi_246, msi_247, msi_248, \
                         msi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_22 * lsi_245[k]
                   + f_3 * pc_x[k] * msi_245[k];

        t_310[k] = f_22 * lsi_246[k]
                   + f_3 * pc_x[k] * msi_246[k];

        t_311[k] = f_22 * lsi_247[k]
                   + f_3 * pc_x[k] * msi_247[k];

        t_312[k] = f_22 * lsi_248[k]
                   + f_3 * pc_x[k] * msi_248[k];

        t_313[k] = f_22 * lsi_249[k]
                   + f_3 * pc_x[k] * msi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, lsi_133, lsi_161, \
                         lsi_250, lsi_251, msh0_183, msh1_183, msi_245, msi_250, \
                         msi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_22 * lsi_250[k]
                   + f_3 * pc_x[k] * msi_250[k];

        t_315[k] = f_22 * lsi_251[k]
                   + f_3 * pc_x[k] * msi_251[k];

        t_316[k] = f_13 * lsi_161[k]
                   + f_1 * msh0_183[k]
                   - f_2 * msh1_183[k]
                   + f_3 * pc_y[k] * msi_245[k];

        t_317[k] = f_14 * lsi_133[k]
                   + f_3 * pc_z[k] * msi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, lsi_163, lsi_164, lsi_165, msh0_185, \
                         msh0_186, msh0_187, msh1_185, msh1_186, msh1_187, msi_247, msi_248, \
                         msi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * lsi_163[k]
                   + f_10 * msh0_185[k]
                   - f_11 * msh1_185[k]
                   + f_3 * pc_y[k] * msi_247[k];

        t_319[k] = f_13 * lsi_164[k]
                   + f_8 * msh0_186[k]
                   - f_9 * msh1_186[k]
                   + f_3 * pc_y[k] * msi_248[k];

        t_320[k] = f_13 * lsi_165[k]
                   + f_6 * msh0_187[k]
                   - f_7 * msh1_187[k]
                   + f_3 * pc_y[k] * msi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, lsk0_215, lsi_166, lsi_167, \
                         lsk1_215, msh0_188, msh1_188, msi_250, \
                         msi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * lsi_166[k]
                   + f_4 * msh0_188[k]
                   - f_5 * msh1_188[k]
                   + f_3 * pc_y[k] * msi_250[k];

        t_322[k] = f_13 * lsi_167[k]
                   + f_3 * pc_y[k] * msi_251[k];

        t_323[k] = pa_y[k] * lsk0_215[k]
                   - f_12 * pc_y[k] * lsk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, lsi_140, \
                         lsi_252, msh0_189, msh1_189, msi_252, msi_253, \
                         msi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_22 * lsi_252[k]
                   + f_1 * msh0_189[k]
                   - f_2 * msh1_189[k]
                   + f_3 * pc_x[k] * msi_252[k];

        t_325[k] = f_3 * pc_y[k] * msi_252[k];

        t_326[k] = f_15 * lsi_140[k]
                   + f_3 * pc_z[k] * msi_252[k];

        t_327[k] = f_4 * msh0_189[k]
                   - f_5 * msh1_189[k]
                   + f_3 * pc_y[k] * msi_253[k];

        t_328[k] = f_3 * pc_y[k] * msi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, lsi_257, msh0_190, msh0_191, \
                         msh0_194, msh1_190, msh1_191, msh1_194, msi_255, msi_256, \
                         msi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_22 * lsi_257[k]
                   + f_10 * msh0_194[k]
                   - f_11 * msh1_194[k]
                   + f_3 * pc_x[k] * msi_257[k];

        t_330[k] = f_6 * msh0_190[k]
                   - f_7 * msh1_190[k]
                   + f_3 * pc_y[k] * msi_255[k];

        t_331[k] = f_4 * msh0_191[k]
                   - f_5 * msh1_191[k]
                   + f_3 * pc_y[k] * msi_256[k];

        t_332[k] = f_3 * pc_y[k] * msi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, lsi_261, msh0_192, msh0_193, \
                         msh0_198, msh1_192, msh1_193, msh1_198, msi_258, msi_259, \
                         msi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_22 * lsi_261[k]
                   + f_8 * msh0_198[k]
                   - f_9 * msh1_198[k]
                   + f_3 * pc_x[k] * msi_261[k];

        t_334[k] = f_8 * msh0_192[k]
                   - f_9 * msh1_192[k]
                   + f_3 * pc_y[k] * msi_258[k];

        t_335[k] = f_6 * msh0_193[k]
                   - f_7 * msh1_193[k]
                   + f_3 * pc_y[k] * msi_259[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, lsi_266, msh0_194, msh0_203, \
                         msh1_194, msh1_203, msi_260, msi_261, \
                         msi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_4 * msh0_194[k]
                   - f_5 * msh1_194[k]
                   + f_3 * pc_y[k] * msi_260[k];

        t_337[k] = f_3 * pc_y[k] * msi_261[k];

        t_338[k] = f_22 * lsi_266[k]
                   + f_6 * msh0_203[k]
                   - f_7 * msh1_203[k]
                   + f_3 * pc_x[k] * msi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_y, msh0_195, msh0_196, msh0_197, msh1_195, \
                         msh1_196, msh1_197, msi_262, msi_263, \
                         msi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_10 * msh0_195[k]
                   - f_11 * msh1_195[k]
                   + f_3 * pc_y[k] * msi_262[k];

        t_340[k] = f_8 * msh0_196[k]
                   - f_9 * msh1_196[k]
                   + f_3 * pc_y[k] * msi_263[k];

        t_341[k] = f_6 * msh0_197[k]
                   - f_7 * msh1_197[k]
                   + f_3 * pc_y[k] * msi_264[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_y, lsi_272, lsi_273, msh0_198, \
                         msh0_209, msh1_198, msh1_209, msi_265, msi_266, msi_272, \
                         msi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_4 * msh0_198[k]
                   - f_5 * msh1_198[k]
                   + f_3 * pc_y[k] * msi_265[k];

        t_343[k] = f_3 * pc_y[k] * msi_266[k];

        t_344[k] = f_22 * lsi_272[k]
                   + f_4 * msh0_209[k]
                   - f_5 * msh1_209[k]
                   + f_3 * pc_x[k] * msi_272[k];

        t_345[k] = f_22 * lsi_273[k]
                   + f_3 * pc_x[k] * msi_273[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pc_x, pc_y, lsi_274, lsi_275, \
                         lsi_276, lsi_277, msi_272, msi_274, msi_275, msi_276, \
                         msi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_22 * lsi_274[k]
                   + f_3 * pc_x[k] * msi_274[k];

        t_347[k] = f_22 * lsi_275[k]
                   + f_3 * pc_x[k] * msi_275[k];

        t_348[k] = f_22 * lsi_276[k]
                   + f_3 * pc_x[k] * msi_276[k];

        t_349[k] = f_22 * lsi_277[k]
                   + f_3 * pc_x[k] * msi_277[k];

        t_350[k] = f_3 * pc_y[k] * msi_272[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_x, pc_y, lsi_279, msh0_204, msh0_205, \
                         msh1_204, msh1_205, msi_273, msi_274, \
                         msi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_22 * lsi_279[k]
                   + f_3 * pc_x[k] * msi_279[k];

        t_352[k] = f_1 * msh0_204[k]
                   - f_2 * msh1_204[k]
                   + f_3 * pc_y[k] * msi_273[k];

        t_353[k] = f_19 * msh0_205[k]
                   - f_20 * msh1_205[k]
                   + f_3 * pc_y[k] * msi_274[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_y, msh0_206, msh0_207, msh0_208, msh1_206, \
                         msh1_207, msh1_208, msi_275, msi_276, \
                         msi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_10 * msh0_206[k]
                   - f_11 * msh1_206[k]
                   + f_3 * pc_y[k] * msi_275[k];

        t_355[k] = f_8 * msh0_207[k]
                   - f_9 * msh1_207[k]
                   + f_3 * pc_y[k] * msi_276[k];

        t_356[k] = f_6 * msh0_208[k]
                   - f_7 * msh1_208[k]
                   + f_3 * pc_y[k] * msi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pc_x, pc_y, pc_z, lsi_167, lsi_280, \
                         msh0_209, msh0_210, msh1_209, msh1_210, msi_278, msi_279, \
                         msi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_4 * msh0_209[k]
                   - f_5 * msh1_209[k]
                   + f_3 * pc_y[k] * msi_278[k];

        t_358[k] = f_3 * pc_y[k] * msi_279[k];

        t_359[k] = f_15 * lsi_167[k]
                   + f_1 * msh0_209[k]
                   - f_2 * msh1_209[k]
                   + f_3 * pc_z[k] * msi_279[k];

        t_360[k] = f_17 * lsi_280[k]
                   + f_1 * msh0_210[k]
                   - f_2 * msh1_210[k]
                   + f_3 * pc_x[k] * msi_280[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, lsi_168, lsi_283, \
                         msh0_213, msh1_213, msi_280, msi_281, \
                         msi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_16 * lsi_168[k]
                   + f_3 * pc_y[k] * msi_280[k];

        t_362[k] = f_3 * pc_z[k] * msi_280[k];

        t_363[k] = f_17 * lsi_283[k]
                   + f_10 * msh0_213[k]
                   - f_11 * msh1_213[k]
                   + f_3 * pc_x[k] * msi_283[k];

        t_364[k] = f_3 * pc_z[k] * msi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pc_x, pc_z, lsi_286, msh0_210, msh0_216, \
                         msh1_210, msh1_216, msi_282, msi_283, \
                         msi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * msh0_210[k]
                   - f_5 * msh1_210[k]
                   + f_3 * pc_z[k] * msi_282[k];

        t_366[k] = f_17 * lsi_286[k]
                   + f_8 * msh0_216[k]
                   - f_9 * msh1_216[k]
                   + f_3 * pc_x[k] * msi_286[k];

        t_367[k] = f_3 * pc_z[k] * msi_283[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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

    const auto *lsk0_216 = buffer.data(lsk0 + 216);
    const auto *lsk0_219 = buffer.data(lsk0 + 219);
    const auto *lsk0_222 = buffer.data(lsk0 + 222);
    const auto *lsk0_226 = buffer.data(lsk0 + 226);
    const auto *lsk0_228 = buffer.data(lsk0 + 228);
    const auto *lsk0_231 = buffer.data(lsk0 + 231);
    const auto *lsk0_233 = buffer.data(lsk0 + 233);
    const auto *lsk0_234 = buffer.data(lsk0 + 234);
    const auto *lsk0_244 = buffer.data(lsk0 + 244);
    const auto *lsk0_324 = buffer.data(lsk0 + 324);
    const auto *lsk0_327 = buffer.data(lsk0 + 327);
    const auto *lsk0_329 = buffer.data(lsk0 + 329);
    const auto *lsk0_330 = buffer.data(lsk0 + 330);
    const auto *lsk0_333 = buffer.data(lsk0 + 333);
    const auto *lsk0_334 = buffer.data(lsk0 + 334);
    const auto *lsk0_336 = buffer.data(lsk0 + 336);

    const auto *lsi_168 = buffer.data(lsi + 168);
    const auto *lsi_171 = buffer.data(lsi + 171);
    const auto *lsi_173 = buffer.data(lsi + 173);
    const auto *lsi_174 = buffer.data(lsi + 174);
    const auto *lsi_175 = buffer.data(lsi + 175);
    const auto *lsi_177 = buffer.data(lsi + 177);
    const auto *lsi_178 = buffer.data(lsi + 178);
    const auto *lsi_179 = buffer.data(lsi + 179);
    const auto *lsi_180 = buffer.data(lsi + 180);
    const auto *lsi_182 = buffer.data(lsi + 182);
    const auto *lsi_189 = buffer.data(lsi + 189);
    const auto *lsi_195 = buffer.data(lsi + 195);
    const auto *lsi_196 = buffer.data(lsi + 196);
    const auto *lsi_198 = buffer.data(lsi + 198);
    const auto *lsi_199 = buffer.data(lsi + 199);
    const auto *lsi_201 = buffer.data(lsi + 201);
    const auto *lsi_202 = buffer.data(lsi + 202);
    const auto *lsi_205 = buffer.data(lsi + 205);
    const auto *lsi_206 = buffer.data(lsi + 206);
    const auto *lsi_210 = buffer.data(lsi + 210);
    const auto *lsi_217 = buffer.data(lsi + 217);
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
    const auto *lsi_238 = buffer.data(lsi + 238);
    const auto *lsi_245 = buffer.data(lsi + 245);
    const auto *lsi_247 = buffer.data(lsi + 247);
    const auto *lsi_248 = buffer.data(lsi + 248);
    const auto *lsi_249 = buffer.data(lsi + 249);
    const auto *lsi_250 = buffer.data(lsi + 250);
    const auto *lsi_251 = buffer.data(lsi + 251);
    const auto *lsi_252 = buffer.data(lsi + 252);
    const auto *lsi_253 = buffer.data(lsi + 253);
    const auto *lsi_254 = buffer.data(lsi + 254);
    const auto *lsi_255 = buffer.data(lsi + 255);
    const auto *lsi_257 = buffer.data(lsi + 257);
    const auto *lsi_258 = buffer.data(lsi + 258);
    const auto *lsi_260 = buffer.data(lsi + 260);
    const auto *lsi_290 = buffer.data(lsi + 290);
    const auto *lsi_295 = buffer.data(lsi + 295);
    const auto *lsi_301 = buffer.data(lsi + 301);
    const auto *lsi_303 = buffer.data(lsi + 303);
    const auto *lsi_304 = buffer.data(lsi + 304);
    const auto *lsi_305 = buffer.data(lsi + 305);
    const auto *lsi_306 = buffer.data(lsi + 306);
    const auto *lsi_307 = buffer.data(lsi + 307);
    const auto *lsi_313 = buffer.data(lsi + 313);
    const auto *lsi_317 = buffer.data(lsi + 317);
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

    const auto *lsk1_216 = buffer.data(lsk1 + 216);
    const auto *lsk1_219 = buffer.data(lsk1 + 219);
    const auto *lsk1_222 = buffer.data(lsk1 + 222);
    const auto *lsk1_226 = buffer.data(lsk1 + 226);
    const auto *lsk1_228 = buffer.data(lsk1 + 228);
    const auto *lsk1_231 = buffer.data(lsk1 + 231);
    const auto *lsk1_233 = buffer.data(lsk1 + 233);
    const auto *lsk1_234 = buffer.data(lsk1 + 234);
    const auto *lsk1_244 = buffer.data(lsk1 + 244);
    const auto *lsk1_324 = buffer.data(lsk1 + 324);
    const auto *lsk1_327 = buffer.data(lsk1 + 327);
    const auto *lsk1_329 = buffer.data(lsk1 + 329);
    const auto *lsk1_330 = buffer.data(lsk1 + 330);
    const auto *lsk1_333 = buffer.data(lsk1 + 333);
    const auto *lsk1_334 = buffer.data(lsk1 + 334);
    const auto *lsk1_336 = buffer.data(lsk1 + 336);

    const auto *msh0_212 = buffer.data(msh0 + 212);
    const auto *msh0_213 = buffer.data(msh0 + 213);
    const auto *msh0_215 = buffer.data(msh0 + 215);
    const auto *msh0_216 = buffer.data(msh0 + 216);
    const auto *msh0_217 = buffer.data(msh0 + 217);
    const auto *msh0_219 = buffer.data(msh0 + 219);
    const auto *msh0_220 = buffer.data(msh0 + 220);
    const auto *msh0_225 = buffer.data(msh0 + 225);
    const auto *msh0_226 = buffer.data(msh0 + 226);
    const auto *msh0_227 = buffer.data(msh0 + 227);
    const auto *msh0_228 = buffer.data(msh0 + 228);
    const auto *msh0_230 = buffer.data(msh0 + 230);
    const auto *msh0_236 = buffer.data(msh0 + 236);
    const auto *msh0_240 = buffer.data(msh0 + 240);
    const auto *msh0_245 = buffer.data(msh0 + 245);
    const auto *msh0_248 = buffer.data(msh0 + 248);
    const auto *msh0_249 = buffer.data(msh0 + 249);
    const auto *msh0_250 = buffer.data(msh0 + 250);
    const auto *msh0_251 = buffer.data(msh0 + 251);
    const auto *msh0_252 = buffer.data(msh0 + 252);
    const auto *msh0_255 = buffer.data(msh0 + 255);
    const auto *msh0_257 = buffer.data(msh0 + 257);
    const auto *msh0_258 = buffer.data(msh0 + 258);
    const auto *msh0_261 = buffer.data(msh0 + 261);
    const auto *msh0_262 = buffer.data(msh0 + 262);
    const auto *msh0_264 = buffer.data(msh0 + 264);
    const auto *msh0_266 = buffer.data(msh0 + 266);
    const auto *msh0_267 = buffer.data(msh0 + 267);
    const auto *msh0_269 = buffer.data(msh0 + 269);
    const auto *msh0_270 = buffer.data(msh0 + 270);
    const auto *msh0_271 = buffer.data(msh0 + 271);
    const auto *msh0_272 = buffer.data(msh0 + 272);

    const auto *msh1_212 = buffer.data(msh1 + 212);
    const auto *msh1_213 = buffer.data(msh1 + 213);
    const auto *msh1_215 = buffer.data(msh1 + 215);
    const auto *msh1_216 = buffer.data(msh1 + 216);
    const auto *msh1_217 = buffer.data(msh1 + 217);
    const auto *msh1_219 = buffer.data(msh1 + 219);
    const auto *msh1_220 = buffer.data(msh1 + 220);
    const auto *msh1_225 = buffer.data(msh1 + 225);
    const auto *msh1_226 = buffer.data(msh1 + 226);
    const auto *msh1_227 = buffer.data(msh1 + 227);
    const auto *msh1_228 = buffer.data(msh1 + 228);
    const auto *msh1_230 = buffer.data(msh1 + 230);
    const auto *msh1_236 = buffer.data(msh1 + 236);
    const auto *msh1_240 = buffer.data(msh1 + 240);
    const auto *msh1_245 = buffer.data(msh1 + 245);
    const auto *msh1_248 = buffer.data(msh1 + 248);
    const auto *msh1_249 = buffer.data(msh1 + 249);
    const auto *msh1_250 = buffer.data(msh1 + 250);
    const auto *msh1_251 = buffer.data(msh1 + 251);
    const auto *msh1_252 = buffer.data(msh1 + 252);
    const auto *msh1_255 = buffer.data(msh1 + 255);
    const auto *msh1_257 = buffer.data(msh1 + 257);
    const auto *msh1_258 = buffer.data(msh1 + 258);
    const auto *msh1_261 = buffer.data(msh1 + 261);
    const auto *msh1_262 = buffer.data(msh1 + 262);
    const auto *msh1_264 = buffer.data(msh1 + 264);
    const auto *msh1_266 = buffer.data(msh1 + 266);
    const auto *msh1_267 = buffer.data(msh1 + 267);
    const auto *msh1_269 = buffer.data(msh1 + 269);
    const auto *msh1_270 = buffer.data(msh1 + 270);
    const auto *msh1_271 = buffer.data(msh1 + 271);
    const auto *msh1_272 = buffer.data(msh1 + 272);

    const auto *msi_285 = buffer.data(msi + 285);
    const auto *msi_286 = buffer.data(msi + 286);
    const auto *msi_287 = buffer.data(msi + 287);
    const auto *msi_289 = buffer.data(msi + 289);
    const auto *msi_290 = buffer.data(msi + 290);
    const auto *msi_291 = buffer.data(msi + 291);
    const auto *msi_292 = buffer.data(msi + 292);
    const auto *msi_294 = buffer.data(msi + 294);
    const auto *msi_295 = buffer.data(msi + 295);
    const auto *msi_301 = buffer.data(msi + 301);
    const auto *msi_302 = buffer.data(msi + 302);
    const auto *msi_303 = buffer.data(msi + 303);
    const auto *msi_304 = buffer.data(msi + 304);
    const auto *msi_305 = buffer.data(msi + 305);
    const auto *msi_306 = buffer.data(msi + 306);
    const auto *msi_307 = buffer.data(msi + 307);
    const auto *msi_308 = buffer.data(msi + 308);
    const auto *msi_310 = buffer.data(msi + 310);
    const auto *msi_311 = buffer.data(msi + 311);
    const auto *msi_313 = buffer.data(msi + 313);
    const auto *msi_314 = buffer.data(msi + 314);
    const auto *msi_317 = buffer.data(msi + 317);
    const auto *msi_318 = buffer.data(msi + 318);
    const auto *msi_322 = buffer.data(msi + 322);
    const auto *msi_328 = buffer.data(msi + 328);
    const auto *msi_329 = buffer.data(msi + 329);
    const auto *msi_330 = buffer.data(msi + 330);
    const auto *msi_331 = buffer.data(msi + 331);
    const auto *msi_332 = buffer.data(msi + 332);
    const auto *msi_333 = buffer.data(msi + 333);
    const auto *msi_334 = buffer.data(msi + 334);
    const auto *msi_335 = buffer.data(msi + 335);
    const auto *msi_336 = buffer.data(msi + 336);
    const auto *msi_338 = buffer.data(msi + 338);
    const auto *msi_339 = buffer.data(msi + 339);
    const auto *msi_341 = buffer.data(msi + 341);
    const auto *msi_342 = buffer.data(msi + 342);
    const auto *msi_345 = buffer.data(msi + 345);
    const auto *msi_346 = buffer.data(msi + 346);
    const auto *msi_348 = buffer.data(msi + 348);
    const auto *msi_350 = buffer.data(msi + 350);
    const auto *msi_351 = buffer.data(msi + 351);
    const auto *msi_353 = buffer.data(msi + 353);
    const auto *msi_354 = buffer.data(msi + 354);
    const auto *msi_356 = buffer.data(msi + 356);
    const auto *msi_357 = buffer.data(msi + 357);
    const auto *msi_358 = buffer.data(msi + 358);
    const auto *msi_359 = buffer.data(msi + 359);
    const auto *msi_360 = buffer.data(msi + 360);
    const auto *msi_361 = buffer.data(msi + 361);
    const auto *msi_362 = buffer.data(msi + 362);
    const auto *msi_363 = buffer.data(msi + 363);
    const auto *msi_364 = buffer.data(msi + 364);
    const auto *msi_366 = buffer.data(msi + 366);
    const auto *msi_367 = buffer.data(msi + 367);
    const auto *msi_369 = buffer.data(msi + 369);
    const auto *msi_370 = buffer.data(msi + 370);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pc_x, pc_y, pc_z, lsi_173, lsi_290, \
                         msh0_212, msh0_220, msh1_212, msh1_220, msi_285, msi_286, \
                         msi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_16 * lsi_173[k]
                   + f_3 * pc_y[k] * msi_285[k];

        t_369[k] = f_6 * msh0_212[k]
                   - f_7 * msh1_212[k]
                   + f_3 * pc_z[k] * msi_285[k];

        t_370[k] = f_17 * lsi_290[k]
                   + f_6 * msh0_220[k]
                   - f_7 * msh1_220[k]
                   + f_3 * pc_x[k] * msi_290[k];

        t_371[k] = f_3 * pc_z[k] * msi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, lsi_177, msh0_213, msh0_215, \
                         msh1_213, msh1_215, msi_287, msi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * msh0_213[k]
                   - f_5 * msh1_213[k]
                   + f_3 * pc_z[k] * msi_287[k];

        t_373[k] = f_16 * lsi_177[k]
                   + f_3 * pc_y[k] * msi_289[k];

        t_374[k] = f_8 * msh0_215[k]
                   - f_9 * msh1_215[k]
                   + f_3 * pc_z[k] * msi_289[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_z, lsi_295, msh0_216, msh0_225, \
                         msh1_216, msh1_225, msi_290, msi_291, \
                         msi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_17 * lsi_295[k]
                   + f_4 * msh0_225[k]
                   - f_5 * msh1_225[k]
                   + f_3 * pc_x[k] * msi_295[k];

        t_376[k] = f_3 * pc_z[k] * msi_290[k];

        t_377[k] = f_4 * msh0_216[k]
                   - f_5 * msh1_216[k]
                   + f_3 * pc_z[k] * msi_291[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, lsi_182, lsi_301, \
                         msh0_217, msh0_219, msh1_217, msh1_219, msi_292, msi_294, \
                         msi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * msh0_217[k]
                   - f_7 * msh1_217[k]
                   + f_3 * pc_z[k] * msi_292[k];

        t_379[k] = f_16 * lsi_182[k]
                   + f_3 * pc_y[k] * msi_294[k];

        t_380[k] = f_10 * msh0_219[k]
                   - f_11 * msh1_219[k]
                   + f_3 * pc_z[k] * msi_294[k];

        t_381[k] = f_17 * lsi_301[k]
                   + f_3 * pc_x[k] * msi_301[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_z, lsi_303, lsi_304, \
                         lsi_305, lsi_306, msi_295, msi_303, msi_304, msi_305, \
                         msi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_z[k] * msi_295[k];

        t_383[k] = f_17 * lsi_303[k]
                   + f_3 * pc_x[k] * msi_303[k];

        t_384[k] = f_17 * lsi_304[k]
                   + f_3 * pc_x[k] * msi_304[k];

        t_385[k] = f_17 * lsi_305[k]
                   + f_3 * pc_x[k] * msi_305[k];

        t_386[k] = f_17 * lsi_306[k]
                   + f_3 * pc_x[k] * msi_306[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pc_y, pc_z, lsi_189, lsi_307, \
                         msh0_225, msh1_225, msi_301, msi_302, \
                         msi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_17 * lsi_307[k]
                   + f_3 * pc_x[k] * msi_307[k];

        t_388[k] = f_16 * lsi_189[k]
                   + f_1 * msh0_225[k]
                   - f_2 * msh1_225[k]
                   + f_3 * pc_y[k] * msi_301[k];

        t_389[k] = f_3 * pc_z[k] * msi_301[k];

        t_390[k] = f_4 * msh0_225[k]
                   - f_5 * msh1_225[k]
                   + f_3 * pc_z[k] * msi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_z, msh0_226, msh0_227, msh0_228, msh1_226, \
                         msh1_227, msh1_228, msi_303, msi_304, \
                         msi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_6 * msh0_226[k]
                   - f_7 * msh1_226[k]
                   + f_3 * pc_z[k] * msi_303[k];

        t_392[k] = f_8 * msh0_227[k]
                   - f_9 * msh1_227[k]
                   + f_3 * pc_z[k] * msi_304[k];

        t_393[k] = f_10 * msh0_228[k]
                   - f_11 * msh1_228[k]
                   + f_3 * pc_z[k] * msi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_z, pc_y, pc_z, lsk0_216, lsi_195, \
                         lsi_196, lsk1_216, msh0_230, msh1_230, msi_307, \
                         msi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * lsi_195[k]
                   + f_3 * pc_y[k] * msi_307[k];

        t_395[k] = f_1 * msh0_230[k]
                   - f_2 * msh1_230[k]
                   + f_3 * pc_z[k] * msi_307[k];

        t_396[k] = pa_z[k] * lsk0_216[k]
                   - f_12 * pc_z[k] * lsk1_216[k];

        t_397[k] = f_15 * lsi_196[k]
                   + f_3 * pc_y[k] * msi_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_z, pc_y, pc_z, lsk0_219, lsi_168, lsi_198, \
                         lsk1_219, msi_308, msi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * lsi_168[k]
                   + f_3 * pc_z[k] * msi_308[k];

        t_399[k] = pa_z[k] * lsk0_219[k]
                   - f_12 * pc_z[k] * lsk1_219[k];

        t_400[k] = f_15 * lsi_198[k]
                   + f_3 * pc_y[k] * msi_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_z, pc_x, pc_z, lsk0_222, lsi_171, lsi_313, \
                         lsk1_222, msh0_236, msh1_236, msi_311, \
                         msi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_17 * lsi_313[k]
                   + f_10 * msh0_236[k]
                   - f_11 * msh1_236[k]
                   + f_3 * pc_x[k] * msi_313[k];

        t_402[k] = pa_z[k] * lsk0_222[k]
                   - f_12 * pc_z[k] * lsk1_222[k];

        t_403[k] = f_13 * lsi_171[k]
                   + f_3 * pc_z[k] * msi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_z, pc_x, pc_y, pc_z, lsk0_226, lsi_201, \
                         lsi_317, lsk1_226, msh0_240, msh1_240, msi_313, \
                         msi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * lsi_201[k]
                   + f_3 * pc_y[k] * msi_313[k];

        t_405[k] = f_17 * lsi_317[k]
                   + f_8 * msh0_240[k]
                   - f_9 * msh1_240[k]
                   + f_3 * pc_x[k] * msi_317[k];

        t_406[k] = pa_z[k] * lsk0_226[k]
                   - f_12 * pc_z[k] * lsk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_z, pc_y, pc_z, lsk0_228, lsi_174, lsi_175, \
                         lsi_205, lsk1_228, msi_314, msi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * lsi_174[k]
                   + f_3 * pc_z[k] * msi_314[k];

        t_408[k] = pa_z[k] * lsk0_228[k]
                   + f_14 * lsi_175[k]
                   - f_12 * pc_z[k] * lsk1_228[k];

        t_409[k] = f_15 * lsi_205[k]
                   + f_3 * pc_y[k] * msi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pa_z, pc_x, pc_z, lsk0_231, lsi_178, lsi_322, \
                         lsk1_231, msh0_245, msh1_245, msi_318, \
                         msi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_17 * lsi_322[k]
                   + f_6 * msh0_245[k]
                   - f_7 * msh1_245[k]
                   + f_3 * pc_x[k] * msi_322[k];

        t_411[k] = pa_z[k] * lsk0_231[k]
                   - f_12 * pc_z[k] * lsk1_231[k];

        t_412[k] = f_13 * lsi_178[k]
                   + f_3 * pc_z[k] * msi_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pa_z, pc_y, pc_z, lsk0_233, lsk0_234, lsi_179, \
                         lsi_180, lsi_210, lsk1_233, lsk1_234, \
                         msi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_z[k] * lsk0_233[k]
                   + f_14 * lsi_179[k]
                   - f_12 * pc_z[k] * lsk1_233[k];

        t_414[k] = pa_z[k] * lsk0_234[k]
                   + f_15 * lsi_180[k]
                   - f_12 * pc_z[k] * lsk1_234[k];

        t_415[k] = f_15 * lsi_210[k]
                   + f_3 * pc_y[k] * msi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, lsi_328, lsi_329, lsi_330, lsi_331, \
                         msh0_251, msh1_251, msi_328, msi_329, msi_330, \
                         msi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_17 * lsi_328[k]
                   + f_4 * msh0_251[k]
                   - f_5 * msh1_251[k]
                   + f_3 * pc_x[k] * msi_328[k];

        t_417[k] = f_17 * lsi_329[k]
                   + f_3 * pc_x[k] * msi_329[k];

        t_418[k] = f_17 * lsi_330[k]
                   + f_3 * pc_x[k] * msi_330[k];

        t_419[k] = f_17 * lsi_331[k]
                   + f_3 * pc_x[k] * msi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, lsi_332, lsi_333, lsi_334, lsi_335, \
                         msi_332, msi_333, msi_334, msi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_17 * lsi_332[k]
                   + f_3 * pc_x[k] * msi_332[k];

        t_421[k] = f_17 * lsi_333[k]
                   + f_3 * pc_x[k] * msi_333[k];

        t_422[k] = f_17 * lsi_334[k]
                   + f_3 * pc_x[k] * msi_334[k];

        t_423[k] = f_17 * lsi_335[k]
                   + f_3 * pc_x[k] * msi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_z, pc_y, pc_z, lsk0_244, lsi_189, lsi_219, \
                         lsk1_244, msh0_248, msh1_248, msi_329, \
                         msi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * lsk0_244[k]
                   - f_12 * pc_z[k] * lsk1_244[k];

        t_425[k] = f_13 * lsi_189[k]
                   + f_3 * pc_z[k] * msi_329[k];

        t_426[k] = f_15 * lsi_219[k]
                   + f_10 * msh0_248[k]
                   - f_11 * msh1_248[k]
                   + f_3 * pc_y[k] * msi_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, lsi_220, lsi_221, lsi_222, msh0_249, \
                         msh0_250, msh0_251, msh1_249, msh1_250, msh1_251, msi_332, msi_333, \
                         msi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * lsi_220[k]
                   + f_8 * msh0_249[k]
                   - f_9 * msh1_249[k]
                   + f_3 * pc_y[k] * msi_332[k];

        t_428[k] = f_15 * lsi_221[k]
                   + f_6 * msh0_250[k]
                   - f_7 * msh1_250[k]
                   + f_3 * pc_y[k] * msi_333[k];

        t_429[k] = f_15 * lsi_222[k]
                   + f_4 * msh0_251[k]
                   - f_5 * msh1_251[k]
                   + f_3 * pc_y[k] * msi_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, lsi_195, lsi_223, lsi_336, \
                         msh0_251, msh0_252, msh1_251, msh1_252, msi_335, \
                         msi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * lsi_223[k]
                   + f_3 * pc_y[k] * msi_335[k];

        t_431[k] = f_13 * lsi_195[k]
                   + f_1 * msh0_251[k]
                   - f_2 * msh1_251[k]
                   + f_3 * pc_z[k] * msi_335[k];

        t_432[k] = f_17 * lsi_336[k]
                   + f_1 * msh0_252[k]
                   - f_2 * msh1_252[k]
                   + f_3 * pc_x[k] * msi_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, lsi_196, lsi_224, \
                         lsi_226, lsi_339, msh0_255, msh1_255, msi_336, msi_338, \
                         msi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * lsi_224[k]
                   + f_3 * pc_y[k] * msi_336[k];

        t_434[k] = f_14 * lsi_196[k]
                   + f_3 * pc_z[k] * msi_336[k];

        t_435[k] = f_17 * lsi_339[k]
                   + f_10 * msh0_255[k]
                   - f_11 * msh1_255[k]
                   + f_3 * pc_x[k] * msi_339[k];

        t_436[k] = f_14 * lsi_226[k]
                   + f_3 * pc_y[k] * msi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, lsi_199, lsi_341, lsi_342, msh0_257, \
                         msh0_258, msh1_257, msh1_258, msi_339, msi_341, \
                         msi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_17 * lsi_341[k]
                   + f_10 * msh0_257[k]
                   - f_11 * msh1_257[k]
                   + f_3 * pc_x[k] * msi_341[k];

        t_438[k] = f_17 * lsi_342[k]
                   + f_8 * msh0_258[k]
                   - f_9 * msh1_258[k]
                   + f_3 * pc_x[k] * msi_342[k];

        t_439[k] = f_14 * lsi_199[k]
                   + f_3 * pc_z[k] * msi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, lsi_229, lsi_345, lsi_346, msh0_261, \
                         msh0_262, msh1_261, msh1_262, msi_341, msi_345, \
                         msi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * lsi_229[k]
                   + f_3 * pc_y[k] * msi_341[k];

        t_441[k] = f_17 * lsi_345[k]
                   + f_8 * msh0_261[k]
                   - f_9 * msh1_261[k]
                   + f_3 * pc_x[k] * msi_345[k];

        t_442[k] = f_17 * lsi_346[k]
                   + f_6 * msh0_262[k]
                   - f_7 * msh1_262[k]
                   + f_3 * pc_x[k] * msi_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, lsi_202, lsi_233, lsi_348, \
                         msh0_264, msh1_264, msi_342, msi_345, \
                         msi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * lsi_202[k]
                   + f_3 * pc_z[k] * msi_342[k];

        t_444[k] = f_17 * lsi_348[k]
                   + f_6 * msh0_264[k]
                   - f_7 * msh1_264[k]
                   + f_3 * pc_x[k] * msi_348[k];

        t_445[k] = f_14 * lsi_233[k]
                   + f_3 * pc_y[k] * msi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, lsi_206, lsi_350, lsi_351, msh0_266, \
                         msh0_267, msh1_266, msh1_267, msi_346, msi_350, \
                         msi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_17 * lsi_350[k]
                   + f_6 * msh0_266[k]
                   - f_7 * msh1_266[k]
                   + f_3 * pc_x[k] * msi_350[k];

        t_447[k] = f_17 * lsi_351[k]
                   + f_4 * msh0_267[k]
                   - f_5 * msh1_267[k]
                   + f_3 * pc_x[k] * msi_351[k];

        t_448[k] = f_14 * lsi_206[k]
                   + f_3 * pc_z[k] * msi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, lsi_238, lsi_353, lsi_354, msh0_269, \
                         msh0_270, msh1_269, msh1_270, msi_350, msi_353, \
                         msi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * lsi_353[k]
                   + f_4 * msh0_269[k]
                   - f_5 * msh1_269[k]
                   + f_3 * pc_x[k] * msi_353[k];

        t_450[k] = f_17 * lsi_354[k]
                   + f_4 * msh0_270[k]
                   - f_5 * msh1_270[k]
                   + f_3 * pc_x[k] * msi_354[k];

        t_451[k] = f_14 * lsi_238[k]
                   + f_3 * pc_y[k] * msi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, lsi_356, lsi_357, lsi_358, lsi_359, \
                         msh0_272, msh1_272, msi_356, msi_357, msi_358, \
                         msi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_17 * lsi_356[k]
                   + f_4 * msh0_272[k]
                   - f_5 * msh1_272[k]
                   + f_3 * pc_x[k] * msi_356[k];

        t_453[k] = f_17 * lsi_357[k]
                   + f_3 * pc_x[k] * msi_357[k];

        t_454[k] = f_17 * lsi_358[k]
                   + f_3 * pc_x[k] * msi_358[k];

        t_455[k] = f_17 * lsi_359[k]
                   + f_3 * pc_x[k] * msi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, lsi_360, lsi_361, lsi_362, lsi_363, \
                         msi_360, msi_361, msi_362, msi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_17 * lsi_360[k]
                   + f_3 * pc_x[k] * msi_360[k];

        t_457[k] = f_17 * lsi_361[k]
                   + f_3 * pc_x[k] * msi_361[k];

        t_458[k] = f_17 * lsi_362[k]
                   + f_3 * pc_x[k] * msi_362[k];

        t_459[k] = f_17 * lsi_363[k]
                   + f_3 * pc_x[k] * msi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, lsi_217, lsi_245, lsi_247, msh0_267, \
                         msh0_269, msh1_267, msh1_269, msi_357, \
                         msi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * lsi_245[k]
                   + f_1 * msh0_267[k]
                   - f_2 * msh1_267[k]
                   + f_3 * pc_y[k] * msi_357[k];

        t_461[k] = f_14 * lsi_217[k]
                   + f_3 * pc_z[k] * msi_357[k];

        t_462[k] = f_14 * lsi_247[k]
                   + f_10 * msh0_269[k]
                   - f_11 * msh1_269[k]
                   + f_3 * pc_y[k] * msi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, lsi_248, lsi_249, lsi_250, msh0_270, \
                         msh0_271, msh0_272, msh1_270, msh1_271, msh1_272, msi_360, msi_361, \
                         msi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * lsi_248[k]
                   + f_8 * msh0_270[k]
                   - f_9 * msh1_270[k]
                   + f_3 * pc_y[k] * msi_360[k];

        t_464[k] = f_14 * lsi_249[k]
                   + f_6 * msh0_271[k]
                   - f_7 * msh1_271[k]
                   + f_3 * pc_y[k] * msi_361[k];

        t_465[k] = f_14 * lsi_250[k]
                   + f_4 * msh0_272[k]
                   - f_5 * msh1_272[k]
                   + f_3 * pc_y[k] * msi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_y, pc_y, pc_z, lsk0_324, lsi_223, \
                         lsi_251, lsi_252, lsk1_324, msh0_272, msh1_272, msi_363, \
                         msi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * lsi_251[k]
                   + f_3 * pc_y[k] * msi_363[k];

        t_467[k] = f_14 * lsi_223[k]
                   + f_1 * msh0_272[k]
                   - f_2 * msh1_272[k]
                   + f_3 * pc_z[k] * msi_363[k];

        t_468[k] = pa_y[k] * lsk0_324[k]
                   - f_12 * pc_y[k] * lsk1_324[k];

        t_469[k] = f_13 * lsi_252[k]
                   + f_3 * pc_y[k] * msi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_y, pc_y, pc_z, lsk0_327, lsk0_329, \
                         lsi_224, lsi_253, lsi_254, lsk1_327, lsk1_329, msi_364, \
                         msi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * lsi_224[k]
                   + f_3 * pc_z[k] * msi_364[k];

        t_471[k] = pa_y[k] * lsk0_327[k]
                   + f_14 * lsi_253[k]
                   - f_12 * pc_y[k] * lsk1_327[k];

        t_472[k] = f_13 * lsi_254[k]
                   + f_3 * pc_y[k] * msi_366[k];

        t_473[k] = pa_y[k] * lsk0_329[k]
                   - f_12 * pc_y[k] * lsk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_y, pc_y, pc_z, lsk0_330, lsk0_333, \
                         lsi_227, lsi_255, lsi_257, lsk1_330, lsk1_333, msi_367, \
                         msi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_y[k] * lsk0_330[k]
                   + f_15 * lsi_255[k]
                   - f_12 * pc_y[k] * lsk1_330[k];

        t_475[k] = f_15 * lsi_227[k]
                   + f_3 * pc_z[k] * msi_367[k];

        t_476[k] = f_13 * lsi_257[k]
                   + f_3 * pc_y[k] * msi_369[k];

        t_477[k] = pa_y[k] * lsk0_333[k]
                   - f_12 * pc_y[k] * lsk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_y, pc_y, pc_z, lsk0_334, lsk0_336, lsi_230, \
                         lsi_258, lsi_260, lsk1_334, lsk1_336, \
                         msi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pa_y[k] * lsk0_334[k]
                   + f_16 * lsi_258[k]
                   - f_12 * pc_y[k] * lsk1_334[k];

        t_479[k] = f_15 * lsi_230[k]
                   + f_3 * pc_z[k] * msi_370[k];

        t_480[k] = pa_y[k] * lsk0_336[k]
                   + f_14 * lsi_260[k]
                   - f_12 * pc_y[k] * lsk1_336[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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

    const auto *lsk0_338 = buffer.data(lsk0 + 338);
    const auto *lsk0_339 = buffer.data(lsk0 + 339);
    const auto *lsk0_341 = buffer.data(lsk0 + 341);
    const auto *lsk0_342 = buffer.data(lsk0 + 342);
    const auto *lsk0_344 = buffer.data(lsk0 + 344);
    const auto *lsk0_359 = buffer.data(lsk0 + 359);
    const auto *lsk0_360 = buffer.data(lsk0 + 360);
    const auto *lsk0_363 = buffer.data(lsk0 + 363);
    const auto *lsk0_366 = buffer.data(lsk0 + 366);
    const auto *lsk0_370 = buffer.data(lsk0 + 370);
    const auto *lsk0_372 = buffer.data(lsk0 + 372);
    const auto *lsk0_375 = buffer.data(lsk0 + 375);
    const auto *lsk0_377 = buffer.data(lsk0 + 377);
    const auto *lsk0_378 = buffer.data(lsk0 + 378);

    const auto *lsi_234 = buffer.data(lsi + 234);
    const auto *lsi_245 = buffer.data(lsi + 245);
    const auto *lsi_252 = buffer.data(lsi + 252);
    const auto *lsi_261 = buffer.data(lsi + 261);
    const auto *lsi_262 = buffer.data(lsi + 262);
    const auto *lsi_264 = buffer.data(lsi + 264);
    const auto *lsi_265 = buffer.data(lsi + 265);
    const auto *lsi_266 = buffer.data(lsi + 266);
    const auto *lsi_273 = buffer.data(lsi + 273);
    const auto *lsi_275 = buffer.data(lsi + 275);
    const auto *lsi_276 = buffer.data(lsi + 276);
    const auto *lsi_277 = buffer.data(lsi + 277);
    const auto *lsi_278 = buffer.data(lsi + 278);
    const auto *lsi_279 = buffer.data(lsi + 279);
    const auto *lsi_280 = buffer.data(lsi + 280);
    const auto *lsi_283 = buffer.data(lsi + 283);
    const auto *lsi_285 = buffer.data(lsi + 285);
    const auto *lsi_286 = buffer.data(lsi + 286);
    const auto *lsi_287 = buffer.data(lsi + 287);
    const auto *lsi_289 = buffer.data(lsi + 289);
    const auto *lsi_290 = buffer.data(lsi + 290);
    const auto *lsi_291 = buffer.data(lsi + 291);
    const auto *lsi_292 = buffer.data(lsi + 292);
    const auto *lsi_294 = buffer.data(lsi + 294);
    const auto *lsi_301 = buffer.data(lsi + 301);
    const auto *lsi_307 = buffer.data(lsi + 307);
    const auto *lsi_308 = buffer.data(lsi + 308);
    const auto *lsi_310 = buffer.data(lsi + 310);
    const auto *lsi_313 = buffer.data(lsi + 313);
    const auto *lsi_317 = buffer.data(lsi + 317);
    const auto *lsi_322 = buffer.data(lsi + 322);
    const auto *lsi_385 = buffer.data(lsi + 385);
    const auto *lsi_386 = buffer.data(lsi + 386);
    const auto *lsi_387 = buffer.data(lsi + 387);
    const auto *lsi_388 = buffer.data(lsi + 388);
    const auto *lsi_389 = buffer.data(lsi + 389);
    const auto *lsi_390 = buffer.data(lsi + 390);
    const auto *lsi_391 = buffer.data(lsi + 391);
    const auto *lsi_392 = buffer.data(lsi + 392);
    const auto *lsi_397 = buffer.data(lsi + 397);
    const auto *lsi_401 = buffer.data(lsi + 401);
    const auto *lsi_406 = buffer.data(lsi + 406);
    const auto *lsi_412 = buffer.data(lsi + 412);
    const auto *lsi_413 = buffer.data(lsi + 413);
    const auto *lsi_414 = buffer.data(lsi + 414);
    const auto *lsi_415 = buffer.data(lsi + 415);
    const auto *lsi_416 = buffer.data(lsi + 416);
    const auto *lsi_417 = buffer.data(lsi + 417);
    const auto *lsi_419 = buffer.data(lsi + 419);
    const auto *lsi_420 = buffer.data(lsi + 420);
    const auto *lsi_423 = buffer.data(lsi + 423);
    const auto *lsi_426 = buffer.data(lsi + 426);
    const auto *lsi_430 = buffer.data(lsi + 430);
    const auto *lsi_435 = buffer.data(lsi + 435);
    const auto *lsi_441 = buffer.data(lsi + 441);
    const auto *lsi_443 = buffer.data(lsi + 443);
    const auto *lsi_444 = buffer.data(lsi + 444);
    const auto *lsi_445 = buffer.data(lsi + 445);
    const auto *lsi_446 = buffer.data(lsi + 446);
    const auto *lsi_447 = buffer.data(lsi + 447);
    const auto *lsi_453 = buffer.data(lsi + 453);
    const auto *lsi_457 = buffer.data(lsi + 457);
    const auto *lsi_462 = buffer.data(lsi + 462);
    const auto *lsi_468 = buffer.data(lsi + 468);
    const auto *lsi_469 = buffer.data(lsi + 469);
    const auto *lsi_470 = buffer.data(lsi + 470);
    const auto *lsi_471 = buffer.data(lsi + 471);

    const auto *lsk1_338 = buffer.data(lsk1 + 338);
    const auto *lsk1_339 = buffer.data(lsk1 + 339);
    const auto *lsk1_341 = buffer.data(lsk1 + 341);
    const auto *lsk1_342 = buffer.data(lsk1 + 342);
    const auto *lsk1_344 = buffer.data(lsk1 + 344);
    const auto *lsk1_359 = buffer.data(lsk1 + 359);
    const auto *lsk1_360 = buffer.data(lsk1 + 360);
    const auto *lsk1_363 = buffer.data(lsk1 + 363);
    const auto *lsk1_366 = buffer.data(lsk1 + 366);
    const auto *lsk1_370 = buffer.data(lsk1 + 370);
    const auto *lsk1_372 = buffer.data(lsk1 + 372);
    const auto *lsk1_375 = buffer.data(lsk1 + 375);
    const auto *lsk1_377 = buffer.data(lsk1 + 377);
    const auto *lsk1_378 = buffer.data(lsk1 + 378);

    const auto *msh0_288 = buffer.data(msh0 + 288);
    const auto *msh0_290 = buffer.data(msh0 + 290);
    const auto *msh0_291 = buffer.data(msh0 + 291);
    const auto *msh0_292 = buffer.data(msh0 + 292);
    const auto *msh0_293 = buffer.data(msh0 + 293);
    const auto *msh0_294 = buffer.data(msh0 + 294);
    const auto *msh0_295 = buffer.data(msh0 + 295);
    const auto *msh0_296 = buffer.data(msh0 + 296);
    const auto *msh0_297 = buffer.data(msh0 + 297);
    const auto *msh0_298 = buffer.data(msh0 + 298);
    const auto *msh0_299 = buffer.data(msh0 + 299);
    const auto *msh0_300 = buffer.data(msh0 + 300);
    const auto *msh0_301 = buffer.data(msh0 + 301);
    const auto *msh0_302 = buffer.data(msh0 + 302);
    const auto *msh0_303 = buffer.data(msh0 + 303);
    const auto *msh0_308 = buffer.data(msh0 + 308);
    const auto *msh0_309 = buffer.data(msh0 + 309);
    const auto *msh0_310 = buffer.data(msh0 + 310);
    const auto *msh0_311 = buffer.data(msh0 + 311);
    const auto *msh0_312 = buffer.data(msh0 + 312);
    const auto *msh0_313 = buffer.data(msh0 + 313);
    const auto *msh0_314 = buffer.data(msh0 + 314);
    const auto *msh0_315 = buffer.data(msh0 + 315);
    const auto *msh0_317 = buffer.data(msh0 + 317);
    const auto *msh0_318 = buffer.data(msh0 + 318);
    const auto *msh0_320 = buffer.data(msh0 + 320);
    const auto *msh0_321 = buffer.data(msh0 + 321);
    const auto *msh0_322 = buffer.data(msh0 + 322);
    const auto *msh0_324 = buffer.data(msh0 + 324);
    const auto *msh0_325 = buffer.data(msh0 + 325);
    const auto *msh0_330 = buffer.data(msh0 + 330);
    const auto *msh0_331 = buffer.data(msh0 + 331);
    const auto *msh0_332 = buffer.data(msh0 + 332);
    const auto *msh0_333 = buffer.data(msh0 + 333);
    const auto *msh0_335 = buffer.data(msh0 + 335);
    const auto *msh0_341 = buffer.data(msh0 + 341);
    const auto *msh0_345 = buffer.data(msh0 + 345);
    const auto *msh0_350 = buffer.data(msh0 + 350);
    const auto *msh0_356 = buffer.data(msh0 + 356);

    const auto *msh1_288 = buffer.data(msh1 + 288);
    const auto *msh1_290 = buffer.data(msh1 + 290);
    const auto *msh1_291 = buffer.data(msh1 + 291);
    const auto *msh1_292 = buffer.data(msh1 + 292);
    const auto *msh1_293 = buffer.data(msh1 + 293);
    const auto *msh1_294 = buffer.data(msh1 + 294);
    const auto *msh1_295 = buffer.data(msh1 + 295);
    const auto *msh1_296 = buffer.data(msh1 + 296);
    const auto *msh1_297 = buffer.data(msh1 + 297);
    const auto *msh1_298 = buffer.data(msh1 + 298);
    const auto *msh1_299 = buffer.data(msh1 + 299);
    const auto *msh1_300 = buffer.data(msh1 + 300);
    const auto *msh1_301 = buffer.data(msh1 + 301);
    const auto *msh1_302 = buffer.data(msh1 + 302);
    const auto *msh1_303 = buffer.data(msh1 + 303);
    const auto *msh1_308 = buffer.data(msh1 + 308);
    const auto *msh1_309 = buffer.data(msh1 + 309);
    const auto *msh1_310 = buffer.data(msh1 + 310);
    const auto *msh1_311 = buffer.data(msh1 + 311);
    const auto *msh1_312 = buffer.data(msh1 + 312);
    const auto *msh1_313 = buffer.data(msh1 + 313);
    const auto *msh1_314 = buffer.data(msh1 + 314);
    const auto *msh1_315 = buffer.data(msh1 + 315);
    const auto *msh1_317 = buffer.data(msh1 + 317);
    const auto *msh1_318 = buffer.data(msh1 + 318);
    const auto *msh1_320 = buffer.data(msh1 + 320);
    const auto *msh1_321 = buffer.data(msh1 + 321);
    const auto *msh1_322 = buffer.data(msh1 + 322);
    const auto *msh1_324 = buffer.data(msh1 + 324);
    const auto *msh1_325 = buffer.data(msh1 + 325);
    const auto *msh1_330 = buffer.data(msh1 + 330);
    const auto *msh1_331 = buffer.data(msh1 + 331);
    const auto *msh1_332 = buffer.data(msh1 + 332);
    const auto *msh1_333 = buffer.data(msh1 + 333);
    const auto *msh1_335 = buffer.data(msh1 + 335);
    const auto *msh1_341 = buffer.data(msh1 + 341);
    const auto *msh1_345 = buffer.data(msh1 + 345);
    const auto *msh1_350 = buffer.data(msh1 + 350);
    const auto *msh1_356 = buffer.data(msh1 + 356);

    const auto *msi_373 = buffer.data(msi + 373);
    const auto *msi_374 = buffer.data(msi + 374);
    const auto *msi_378 = buffer.data(msi + 378);
    const auto *msi_385 = buffer.data(msi + 385);
    const auto *msi_386 = buffer.data(msi + 386);
    const auto *msi_387 = buffer.data(msi + 387);
    const auto *msi_388 = buffer.data(msi + 388);
    const auto *msi_389 = buffer.data(msi + 389);
    const auto *msi_390 = buffer.data(msi + 390);
    const auto *msi_391 = buffer.data(msi + 391);
    const auto *msi_392 = buffer.data(msi + 392);
    const auto *msi_393 = buffer.data(msi + 393);
    const auto *msi_394 = buffer.data(msi + 394);
    const auto *msi_395 = buffer.data(msi + 395);
    const auto *msi_396 = buffer.data(msi + 396);
    const auto *msi_397 = buffer.data(msi + 397);
    const auto *msi_398 = buffer.data(msi + 398);
    const auto *msi_399 = buffer.data(msi + 399);
    const auto *msi_400 = buffer.data(msi + 400);
    const auto *msi_401 = buffer.data(msi + 401);
    const auto *msi_402 = buffer.data(msi + 402);
    const auto *msi_403 = buffer.data(msi + 403);
    const auto *msi_404 = buffer.data(msi + 404);
    const auto *msi_405 = buffer.data(msi + 405);
    const auto *msi_406 = buffer.data(msi + 406);
    const auto *msi_412 = buffer.data(msi + 412);
    const auto *msi_413 = buffer.data(msi + 413);
    const auto *msi_414 = buffer.data(msi + 414);
    const auto *msi_415 = buffer.data(msi + 415);
    const auto *msi_416 = buffer.data(msi + 416);
    const auto *msi_417 = buffer.data(msi + 417);
    const auto *msi_418 = buffer.data(msi + 418);
    const auto *msi_419 = buffer.data(msi + 419);
    const auto *msi_420 = buffer.data(msi + 420);
    const auto *msi_421 = buffer.data(msi + 421);
    const auto *msi_422 = buffer.data(msi + 422);
    const auto *msi_423 = buffer.data(msi + 423);
    const auto *msi_425 = buffer.data(msi + 425);
    const auto *msi_426 = buffer.data(msi + 426);
    const auto *msi_427 = buffer.data(msi + 427);
    const auto *msi_429 = buffer.data(msi + 429);
    const auto *msi_430 = buffer.data(msi + 430);
    const auto *msi_431 = buffer.data(msi + 431);
    const auto *msi_432 = buffer.data(msi + 432);
    const auto *msi_434 = buffer.data(msi + 434);
    const auto *msi_435 = buffer.data(msi + 435);
    const auto *msi_441 = buffer.data(msi + 441);
    const auto *msi_442 = buffer.data(msi + 442);
    const auto *msi_443 = buffer.data(msi + 443);
    const auto *msi_444 = buffer.data(msi + 444);
    const auto *msi_445 = buffer.data(msi + 445);
    const auto *msi_446 = buffer.data(msi + 446);
    const auto *msi_447 = buffer.data(msi + 447);
    const auto *msi_448 = buffer.data(msi + 448);
    const auto *msi_450 = buffer.data(msi + 450);
    const auto *msi_451 = buffer.data(msi + 451);
    const auto *msi_453 = buffer.data(msi + 453);
    const auto *msi_454 = buffer.data(msi + 454);
    const auto *msi_457 = buffer.data(msi + 457);
    const auto *msi_458 = buffer.data(msi + 458);
    const auto *msi_462 = buffer.data(msi + 462);
    const auto *msi_468 = buffer.data(msi + 468);
    const auto *msi_469 = buffer.data(msi + 469);
    const auto *msi_470 = buffer.data(msi + 470);
    const auto *msi_471 = buffer.data(msi + 471);

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_y, pc_y, pc_z, lsk0_338, lsk0_339, \
                         lsi_234, lsi_261, lsi_262, lsk1_338, lsk1_339, msi_373, \
                         msi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * lsi_261[k]
                   + f_3 * pc_y[k] * msi_373[k];

        t_482[k] = pa_y[k] * lsk0_338[k]
                   - f_12 * pc_y[k] * lsk1_338[k];

        t_483[k] = pa_y[k] * lsk0_339[k]
                   + f_17 * lsi_262[k]
                   - f_12 * pc_y[k] * lsk1_339[k];

        t_484[k] = f_15 * lsi_234[k]
                   + f_3 * pc_z[k] * msi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_y, pc_y, lsk0_341, lsk0_342, lsk0_344, \
                         lsi_264, lsi_265, lsi_266, lsk1_341, lsk1_342, lsk1_344, \
                         msi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_y[k] * lsk0_341[k]
                   + f_15 * lsi_264[k]
                   - f_12 * pc_y[k] * lsk1_341[k];

        t_486[k] = pa_y[k] * lsk0_342[k]
                   + f_14 * lsi_265[k]
                   - f_12 * pc_y[k] * lsk1_342[k];

        t_487[k] = f_13 * lsi_266[k]
                   + f_3 * pc_y[k] * msi_378[k];

        t_488[k] = pa_y[k] * lsk0_344[k]
                   - f_12 * pc_y[k] * lsk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, lsi_385, lsi_386, lsi_387, \
                         lsi_388, lsi_389, msi_385, msi_386, msi_387, msi_388, \
                         msi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_17 * lsi_385[k]
                   + f_3 * pc_x[k] * msi_385[k];

        t_490[k] = f_17 * lsi_386[k]
                   + f_3 * pc_x[k] * msi_386[k];

        t_491[k] = f_17 * lsi_387[k]
                   + f_3 * pc_x[k] * msi_387[k];

        t_492[k] = f_17 * lsi_388[k]
                   + f_3 * pc_x[k] * msi_388[k];

        t_493[k] = f_17 * lsi_389[k]
                   + f_3 * pc_x[k] * msi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, lsi_245, lsi_273, \
                         lsi_390, lsi_391, msh0_288, msh1_288, msi_385, msi_390, \
                         msi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_17 * lsi_390[k]
                   + f_3 * pc_x[k] * msi_390[k];

        t_495[k] = f_17 * lsi_391[k]
                   + f_3 * pc_x[k] * msi_391[k];

        t_496[k] = f_13 * lsi_273[k]
                   + f_1 * msh0_288[k]
                   - f_2 * msh1_288[k]
                   + f_3 * pc_y[k] * msi_385[k];

        t_497[k] = f_15 * lsi_245[k]
                   + f_3 * pc_z[k] * msi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, lsi_275, lsi_276, lsi_277, msh0_290, \
                         msh0_291, msh0_292, msh1_290, msh1_291, msh1_292, msi_387, msi_388, \
                         msi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * lsi_275[k]
                   + f_10 * msh0_290[k]
                   - f_11 * msh1_290[k]
                   + f_3 * pc_y[k] * msi_387[k];

        t_499[k] = f_13 * lsi_276[k]
                   + f_8 * msh0_291[k]
                   - f_9 * msh1_291[k]
                   + f_3 * pc_y[k] * msi_388[k];

        t_500[k] = f_13 * lsi_277[k]
                   + f_6 * msh0_292[k]
                   - f_7 * msh1_292[k]
                   + f_3 * pc_y[k] * msi_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pc_y, lsk0_359, lsi_278, lsi_279, \
                         lsk1_359, msh0_293, msh1_293, msi_390, \
                         msi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * lsi_278[k]
                   + f_4 * msh0_293[k]
                   - f_5 * msh1_293[k]
                   + f_3 * pc_y[k] * msi_390[k];

        t_502[k] = f_13 * lsi_279[k]
                   + f_3 * pc_y[k] * msi_391[k];

        t_503[k] = pa_y[k] * lsk0_359[k]
                   - f_12 * pc_y[k] * lsk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, lsi_252, \
                         lsi_392, msh0_294, msh1_294, msi_392, msi_393, \
                         msi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_17 * lsi_392[k]
                   + f_1 * msh0_294[k]
                   - f_2 * msh1_294[k]
                   + f_3 * pc_x[k] * msi_392[k];

        t_505[k] = f_3 * pc_y[k] * msi_392[k];

        t_506[k] = f_16 * lsi_252[k]
                   + f_3 * pc_z[k] * msi_392[k];

        t_507[k] = f_4 * msh0_294[k]
                   - f_5 * msh1_294[k]
                   + f_3 * pc_y[k] * msi_393[k];

        t_508[k] = f_3 * pc_y[k] * msi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, lsi_397, msh0_295, msh0_296, \
                         msh0_299, msh1_295, msh1_296, msh1_299, msi_395, msi_396, \
                         msi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_17 * lsi_397[k]
                   + f_10 * msh0_299[k]
                   - f_11 * msh1_299[k]
                   + f_3 * pc_x[k] * msi_397[k];

        t_510[k] = f_6 * msh0_295[k]
                   - f_7 * msh1_295[k]
                   + f_3 * pc_y[k] * msi_395[k];

        t_511[k] = f_4 * msh0_296[k]
                   - f_5 * msh1_296[k]
                   + f_3 * pc_y[k] * msi_396[k];

        t_512[k] = f_3 * pc_y[k] * msi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_y, lsi_401, msh0_297, msh0_298, \
                         msh0_303, msh1_297, msh1_298, msh1_303, msi_398, msi_399, \
                         msi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_17 * lsi_401[k]
                   + f_8 * msh0_303[k]
                   - f_9 * msh1_303[k]
                   + f_3 * pc_x[k] * msi_401[k];

        t_514[k] = f_8 * msh0_297[k]
                   - f_9 * msh1_297[k]
                   + f_3 * pc_y[k] * msi_398[k];

        t_515[k] = f_6 * msh0_298[k]
                   - f_7 * msh1_298[k]
                   + f_3 * pc_y[k] * msi_399[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_x, pc_y, lsi_406, msh0_299, msh0_308, \
                         msh1_299, msh1_308, msi_400, msi_401, \
                         msi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_4 * msh0_299[k]
                   - f_5 * msh1_299[k]
                   + f_3 * pc_y[k] * msi_400[k];

        t_517[k] = f_3 * pc_y[k] * msi_401[k];

        t_518[k] = f_17 * lsi_406[k]
                   + f_6 * msh0_308[k]
                   - f_7 * msh1_308[k]
                   + f_3 * pc_x[k] * msi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, msh0_300, msh0_301, msh0_302, msh1_300, \
                         msh1_301, msh1_302, msi_402, msi_403, \
                         msi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_10 * msh0_300[k]
                   - f_11 * msh1_300[k]
                   + f_3 * pc_y[k] * msi_402[k];

        t_520[k] = f_8 * msh0_301[k]
                   - f_9 * msh1_301[k]
                   + f_3 * pc_y[k] * msi_403[k];

        t_521[k] = f_6 * msh0_302[k]
                   - f_7 * msh1_302[k]
                   + f_3 * pc_y[k] * msi_404[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, lsi_412, lsi_413, msh0_303, \
                         msh0_314, msh1_303, msh1_314, msi_405, msi_406, msi_412, \
                         msi_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_4 * msh0_303[k]
                   - f_5 * msh1_303[k]
                   + f_3 * pc_y[k] * msi_405[k];

        t_523[k] = f_3 * pc_y[k] * msi_406[k];

        t_524[k] = f_17 * lsi_412[k]
                   + f_4 * msh0_314[k]
                   - f_5 * msh1_314[k]
                   + f_3 * pc_x[k] * msi_412[k];

        t_525[k] = f_17 * lsi_413[k]
                   + f_3 * pc_x[k] * msi_413[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, pc_x, pc_y, lsi_414, lsi_415, \
                         lsi_416, lsi_417, msi_412, msi_414, msi_415, msi_416, \
                         msi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_17 * lsi_414[k]
                   + f_3 * pc_x[k] * msi_414[k];

        t_527[k] = f_17 * lsi_415[k]
                   + f_3 * pc_x[k] * msi_415[k];

        t_528[k] = f_17 * lsi_416[k]
                   + f_3 * pc_x[k] * msi_416[k];

        t_529[k] = f_17 * lsi_417[k]
                   + f_3 * pc_x[k] * msi_417[k];

        t_530[k] = f_3 * pc_y[k] * msi_412[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pc_x, pc_y, lsi_419, msh0_309, msh0_310, \
                         msh1_309, msh1_310, msi_413, msi_414, \
                         msi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_17 * lsi_419[k]
                   + f_3 * pc_x[k] * msi_419[k];

        t_532[k] = f_1 * msh0_309[k]
                   - f_2 * msh1_309[k]
                   + f_3 * pc_y[k] * msi_413[k];

        t_533[k] = f_19 * msh0_310[k]
                   - f_20 * msh1_310[k]
                   + f_3 * pc_y[k] * msi_414[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, msh0_311, msh0_312, msh0_313, msh1_311, \
                         msh1_312, msh1_313, msi_415, msi_416, \
                         msi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * msh0_311[k]
                   - f_11 * msh1_311[k]
                   + f_3 * pc_y[k] * msi_415[k];

        t_535[k] = f_8 * msh0_312[k]
                   - f_9 * msh1_312[k]
                   + f_3 * pc_y[k] * msi_416[k];

        t_536[k] = f_6 * msh0_313[k]
                   - f_7 * msh1_313[k]
                   + f_3 * pc_y[k] * msi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, lsi_279, lsi_420, \
                         msh0_314, msh0_315, msh1_314, msh1_315, msi_418, msi_419, \
                         msi_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * msh0_314[k]
                   - f_5 * msh1_314[k]
                   + f_3 * pc_y[k] * msi_418[k];

        t_538[k] = f_3 * pc_y[k] * msi_419[k];

        t_539[k] = f_16 * lsi_279[k]
                   + f_1 * msh0_314[k]
                   - f_2 * msh1_314[k]
                   + f_3 * pc_z[k] * msi_419[k];

        t_540[k] = f_16 * lsi_420[k]
                   + f_1 * msh0_315[k]
                   - f_2 * msh1_315[k]
                   + f_3 * pc_x[k] * msi_420[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, lsi_280, lsi_423, \
                         msh0_318, msh1_318, msi_420, msi_421, \
                         msi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_17 * lsi_280[k]
                   + f_3 * pc_y[k] * msi_420[k];

        t_542[k] = f_3 * pc_z[k] * msi_420[k];

        t_543[k] = f_16 * lsi_423[k]
                   + f_10 * msh0_318[k]
                   - f_11 * msh1_318[k]
                   + f_3 * pc_x[k] * msi_423[k];

        t_544[k] = f_3 * pc_z[k] * msi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, lsi_426, msh0_315, msh0_321, \
                         msh1_315, msh1_321, msi_422, msi_423, \
                         msi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * msh0_315[k]
                   - f_5 * msh1_315[k]
                   + f_3 * pc_z[k] * msi_422[k];

        t_546[k] = f_16 * lsi_426[k]
                   + f_8 * msh0_321[k]
                   - f_9 * msh1_321[k]
                   + f_3 * pc_x[k] * msi_426[k];

        t_547[k] = f_3 * pc_z[k] * msi_423[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pc_x, pc_y, pc_z, lsi_285, lsi_430, \
                         msh0_317, msh0_325, msh1_317, msh1_325, msi_425, msi_426, \
                         msi_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * lsi_285[k]
                   + f_3 * pc_y[k] * msi_425[k];

        t_549[k] = f_6 * msh0_317[k]
                   - f_7 * msh1_317[k]
                   + f_3 * pc_z[k] * msi_425[k];

        t_550[k] = f_16 * lsi_430[k]
                   + f_6 * msh0_325[k]
                   - f_7 * msh1_325[k]
                   + f_3 * pc_x[k] * msi_430[k];

        t_551[k] = f_3 * pc_z[k] * msi_426[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_y, pc_z, lsi_289, msh0_318, msh0_320, \
                         msh1_318, msh1_320, msi_427, msi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_4 * msh0_318[k]
                   - f_5 * msh1_318[k]
                   + f_3 * pc_z[k] * msi_427[k];

        t_553[k] = f_17 * lsi_289[k]
                   + f_3 * pc_y[k] * msi_429[k];

        t_554[k] = f_8 * msh0_320[k]
                   - f_9 * msh1_320[k]
                   + f_3 * pc_z[k] * msi_429[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pc_x, pc_z, lsi_435, msh0_321, msh0_330, \
                         msh1_321, msh1_330, msi_430, msi_431, \
                         msi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_16 * lsi_435[k]
                   + f_4 * msh0_330[k]
                   - f_5 * msh1_330[k]
                   + f_3 * pc_x[k] * msi_435[k];

        t_556[k] = f_3 * pc_z[k] * msi_430[k];

        t_557[k] = f_4 * msh0_321[k]
                   - f_5 * msh1_321[k]
                   + f_3 * pc_z[k] * msi_431[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pc_x, pc_y, pc_z, lsi_294, lsi_441, \
                         msh0_322, msh0_324, msh1_322, msh1_324, msi_432, msi_434, \
                         msi_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_6 * msh0_322[k]
                   - f_7 * msh1_322[k]
                   + f_3 * pc_z[k] * msi_432[k];

        t_559[k] = f_17 * lsi_294[k]
                   + f_3 * pc_y[k] * msi_434[k];

        t_560[k] = f_10 * msh0_324[k]
                   - f_11 * msh1_324[k]
                   + f_3 * pc_z[k] * msi_434[k];

        t_561[k] = f_16 * lsi_441[k]
                   + f_3 * pc_x[k] * msi_441[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pc_x, pc_z, lsi_443, lsi_444, \
                         lsi_445, lsi_446, msi_435, msi_443, msi_444, msi_445, \
                         msi_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_3 * pc_z[k] * msi_435[k];

        t_563[k] = f_16 * lsi_443[k]
                   + f_3 * pc_x[k] * msi_443[k];

        t_564[k] = f_16 * lsi_444[k]
                   + f_3 * pc_x[k] * msi_444[k];

        t_565[k] = f_16 * lsi_445[k]
                   + f_3 * pc_x[k] * msi_445[k];

        t_566[k] = f_16 * lsi_446[k]
                   + f_3 * pc_x[k] * msi_446[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, lsi_301, lsi_447, \
                         msh0_330, msh1_330, msi_441, msi_442, \
                         msi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_16 * lsi_447[k]
                   + f_3 * pc_x[k] * msi_447[k];

        t_568[k] = f_17 * lsi_301[k]
                   + f_1 * msh0_330[k]
                   - f_2 * msh1_330[k]
                   + f_3 * pc_y[k] * msi_441[k];

        t_569[k] = f_3 * pc_z[k] * msi_441[k];

        t_570[k] = f_4 * msh0_330[k]
                   - f_5 * msh1_330[k]
                   + f_3 * pc_z[k] * msi_442[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_z, msh0_331, msh0_332, msh0_333, msh1_331, \
                         msh1_332, msh1_333, msi_443, msi_444, \
                         msi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_6 * msh0_331[k]
                   - f_7 * msh1_331[k]
                   + f_3 * pc_z[k] * msi_443[k];

        t_572[k] = f_8 * msh0_332[k]
                   - f_9 * msh1_332[k]
                   + f_3 * pc_z[k] * msi_444[k];

        t_573[k] = f_10 * msh0_333[k]
                   - f_11 * msh1_333[k]
                   + f_3 * pc_z[k] * msi_445[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_z, pc_y, pc_z, lsk0_360, lsi_307, \
                         lsi_308, lsk1_360, msh0_335, msh1_335, msi_447, \
                         msi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * lsi_307[k]
                   + f_3 * pc_y[k] * msi_447[k];

        t_575[k] = f_1 * msh0_335[k]
                   - f_2 * msh1_335[k]
                   + f_3 * pc_z[k] * msi_447[k];

        t_576[k] = pa_z[k] * lsk0_360[k]
                   - f_12 * pc_z[k] * lsk1_360[k];

        t_577[k] = f_16 * lsi_308[k]
                   + f_3 * pc_y[k] * msi_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_z, pc_y, pc_z, lsk0_363, lsi_280, lsi_310, \
                         lsk1_363, msi_448, msi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * lsi_280[k]
                   + f_3 * pc_z[k] * msi_448[k];

        t_579[k] = pa_z[k] * lsk0_363[k]
                   - f_12 * pc_z[k] * lsk1_363[k];

        t_580[k] = f_16 * lsi_310[k]
                   + f_3 * pc_y[k] * msi_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_z, pc_x, pc_z, lsk0_366, lsi_283, lsi_453, \
                         lsk1_366, msh0_341, msh1_341, msi_451, \
                         msi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * lsi_453[k]
                   + f_10 * msh0_341[k]
                   - f_11 * msh1_341[k]
                   + f_3 * pc_x[k] * msi_453[k];

        t_582[k] = pa_z[k] * lsk0_366[k]
                   - f_12 * pc_z[k] * lsk1_366[k];

        t_583[k] = f_13 * lsi_283[k]
                   + f_3 * pc_z[k] * msi_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pa_z, pc_x, pc_y, pc_z, lsk0_370, lsi_313, \
                         lsi_457, lsk1_370, msh0_345, msh1_345, msi_453, \
                         msi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * lsi_313[k]
                   + f_3 * pc_y[k] * msi_453[k];

        t_585[k] = f_16 * lsi_457[k]
                   + f_8 * msh0_345[k]
                   - f_9 * msh1_345[k]
                   + f_3 * pc_x[k] * msi_457[k];

        t_586[k] = pa_z[k] * lsk0_370[k]
                   - f_12 * pc_z[k] * lsk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pa_z, pc_y, pc_z, lsk0_372, lsi_286, lsi_287, \
                         lsi_317, lsk1_372, msi_454, msi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * lsi_286[k]
                   + f_3 * pc_z[k] * msi_454[k];

        t_588[k] = pa_z[k] * lsk0_372[k]
                   + f_14 * lsi_287[k]
                   - f_12 * pc_z[k] * lsk1_372[k];

        t_589[k] = f_16 * lsi_317[k]
                   + f_3 * pc_y[k] * msi_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pa_z, pc_x, pc_z, lsk0_375, lsi_290, lsi_462, \
                         lsk1_375, msh0_350, msh1_350, msi_458, \
                         msi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_16 * lsi_462[k]
                   + f_6 * msh0_350[k]
                   - f_7 * msh1_350[k]
                   + f_3 * pc_x[k] * msi_462[k];

        t_591[k] = pa_z[k] * lsk0_375[k]
                   - f_12 * pc_z[k] * lsk1_375[k];

        t_592[k] = f_13 * lsi_290[k]
                   + f_3 * pc_z[k] * msi_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pa_z, pc_y, pc_z, lsk0_377, lsk0_378, lsi_291, \
                         lsi_292, lsi_322, lsk1_377, lsk1_378, \
                         msi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pa_z[k] * lsk0_377[k]
                   + f_14 * lsi_291[k]
                   - f_12 * pc_z[k] * lsk1_377[k];

        t_594[k] = pa_z[k] * lsk0_378[k]
                   + f_15 * lsi_292[k]
                   - f_12 * pc_z[k] * lsk1_378[k];

        t_595[k] = f_16 * lsi_322[k]
                   + f_3 * pc_y[k] * msi_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, lsi_468, lsi_469, lsi_470, lsi_471, \
                         msh0_356, msh1_356, msi_468, msi_469, msi_470, \
                         msi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_16 * lsi_468[k]
                   + f_4 * msh0_356[k]
                   - f_5 * msh1_356[k]
                   + f_3 * pc_x[k] * msi_468[k];

        t_597[k] = f_16 * lsi_469[k]
                   + f_3 * pc_x[k] * msi_469[k];

        t_598[k] = f_16 * lsi_470[k]
                   + f_3 * pc_x[k] * msi_470[k];

        t_599[k] = f_16 * lsi_471[k]
                   + f_3 * pc_x[k] * msi_471[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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

    const auto *lsk0_388 = buffer.data(lsk0 + 388);
    const auto *lsk0_504 = buffer.data(lsk0 + 504);
    const auto *lsk0_507 = buffer.data(lsk0 + 507);
    const auto *lsk0_509 = buffer.data(lsk0 + 509);
    const auto *lsk0_510 = buffer.data(lsk0 + 510);
    const auto *lsk0_513 = buffer.data(lsk0 + 513);
    const auto *lsk0_514 = buffer.data(lsk0 + 514);
    const auto *lsk0_516 = buffer.data(lsk0 + 516);
    const auto *lsk0_518 = buffer.data(lsk0 + 518);
    const auto *lsk0_519 = buffer.data(lsk0 + 519);
    const auto *lsk0_521 = buffer.data(lsk0 + 521);
    const auto *lsk0_522 = buffer.data(lsk0 + 522);
    const auto *lsk0_524 = buffer.data(lsk0 + 524);

    const auto *lsi_301 = buffer.data(lsi + 301);
    const auto *lsi_307 = buffer.data(lsi + 307);
    const auto *lsi_308 = buffer.data(lsi + 308);
    const auto *lsi_311 = buffer.data(lsi + 311);
    const auto *lsi_314 = buffer.data(lsi + 314);
    const auto *lsi_318 = buffer.data(lsi + 318);
    const auto *lsi_329 = buffer.data(lsi + 329);
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
    const auto *lsi_350 = buffer.data(lsi + 350);
    const auto *lsi_357 = buffer.data(lsi + 357);
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
    const auto *lsi_373 = buffer.data(lsi + 373);
    const auto *lsi_374 = buffer.data(lsi + 374);
    const auto *lsi_378 = buffer.data(lsi + 378);
    const auto *lsi_385 = buffer.data(lsi + 385);
    const auto *lsi_387 = buffer.data(lsi + 387);
    const auto *lsi_388 = buffer.data(lsi + 388);
    const auto *lsi_389 = buffer.data(lsi + 389);
    const auto *lsi_390 = buffer.data(lsi + 390);
    const auto *lsi_391 = buffer.data(lsi + 391);
    const auto *lsi_392 = buffer.data(lsi + 392);
    const auto *lsi_393 = buffer.data(lsi + 393);
    const auto *lsi_394 = buffer.data(lsi + 394);
    const auto *lsi_395 = buffer.data(lsi + 395);
    const auto *lsi_397 = buffer.data(lsi + 397);
    const auto *lsi_398 = buffer.data(lsi + 398);
    const auto *lsi_400 = buffer.data(lsi + 400);
    const auto *lsi_401 = buffer.data(lsi + 401);
    const auto *lsi_402 = buffer.data(lsi + 402);
    const auto *lsi_404 = buffer.data(lsi + 404);
    const auto *lsi_405 = buffer.data(lsi + 405);
    const auto *lsi_406 = buffer.data(lsi + 406);
    const auto *lsi_472 = buffer.data(lsi + 472);
    const auto *lsi_473 = buffer.data(lsi + 473);
    const auto *lsi_474 = buffer.data(lsi + 474);
    const auto *lsi_475 = buffer.data(lsi + 475);
    const auto *lsi_476 = buffer.data(lsi + 476);
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
    const auto *lsi_553 = buffer.data(lsi + 553);
    const auto *lsi_554 = buffer.data(lsi + 554);
    const auto *lsi_555 = buffer.data(lsi + 555);
    const auto *lsi_556 = buffer.data(lsi + 556);
    const auto *lsi_557 = buffer.data(lsi + 557);

    const auto *lsk1_388 = buffer.data(lsk1 + 388);
    const auto *lsk1_504 = buffer.data(lsk1 + 504);
    const auto *lsk1_507 = buffer.data(lsk1 + 507);
    const auto *lsk1_509 = buffer.data(lsk1 + 509);
    const auto *lsk1_510 = buffer.data(lsk1 + 510);
    const auto *lsk1_513 = buffer.data(lsk1 + 513);
    const auto *lsk1_514 = buffer.data(lsk1 + 514);
    const auto *lsk1_516 = buffer.data(lsk1 + 516);
    const auto *lsk1_518 = buffer.data(lsk1 + 518);
    const auto *lsk1_519 = buffer.data(lsk1 + 519);
    const auto *lsk1_521 = buffer.data(lsk1 + 521);
    const auto *lsk1_522 = buffer.data(lsk1 + 522);
    const auto *lsk1_524 = buffer.data(lsk1 + 524);

    const auto *msh0_353 = buffer.data(msh0 + 353);
    const auto *msh0_354 = buffer.data(msh0 + 354);
    const auto *msh0_355 = buffer.data(msh0 + 355);
    const auto *msh0_356 = buffer.data(msh0 + 356);
    const auto *msh0_357 = buffer.data(msh0 + 357);
    const auto *msh0_360 = buffer.data(msh0 + 360);
    const auto *msh0_362 = buffer.data(msh0 + 362);
    const auto *msh0_363 = buffer.data(msh0 + 363);
    const auto *msh0_366 = buffer.data(msh0 + 366);
    const auto *msh0_367 = buffer.data(msh0 + 367);
    const auto *msh0_369 = buffer.data(msh0 + 369);
    const auto *msh0_371 = buffer.data(msh0 + 371);
    const auto *msh0_372 = buffer.data(msh0 + 372);
    const auto *msh0_374 = buffer.data(msh0 + 374);
    const auto *msh0_375 = buffer.data(msh0 + 375);
    const auto *msh0_376 = buffer.data(msh0 + 376);
    const auto *msh0_377 = buffer.data(msh0 + 377);
    const auto *msh0_378 = buffer.data(msh0 + 378);
    const auto *msh0_381 = buffer.data(msh0 + 381);
    const auto *msh0_383 = buffer.data(msh0 + 383);
    const auto *msh0_384 = buffer.data(msh0 + 384);
    const auto *msh0_387 = buffer.data(msh0 + 387);
    const auto *msh0_388 = buffer.data(msh0 + 388);
    const auto *msh0_390 = buffer.data(msh0 + 390);
    const auto *msh0_392 = buffer.data(msh0 + 392);
    const auto *msh0_393 = buffer.data(msh0 + 393);
    const auto *msh0_395 = buffer.data(msh0 + 395);
    const auto *msh0_396 = buffer.data(msh0 + 396);
    const auto *msh0_397 = buffer.data(msh0 + 397);
    const auto *msh0_398 = buffer.data(msh0 + 398);

    const auto *msh1_353 = buffer.data(msh1 + 353);
    const auto *msh1_354 = buffer.data(msh1 + 354);
    const auto *msh1_355 = buffer.data(msh1 + 355);
    const auto *msh1_356 = buffer.data(msh1 + 356);
    const auto *msh1_357 = buffer.data(msh1 + 357);
    const auto *msh1_360 = buffer.data(msh1 + 360);
    const auto *msh1_362 = buffer.data(msh1 + 362);
    const auto *msh1_363 = buffer.data(msh1 + 363);
    const auto *msh1_366 = buffer.data(msh1 + 366);
    const auto *msh1_367 = buffer.data(msh1 + 367);
    const auto *msh1_369 = buffer.data(msh1 + 369);
    const auto *msh1_371 = buffer.data(msh1 + 371);
    const auto *msh1_372 = buffer.data(msh1 + 372);
    const auto *msh1_374 = buffer.data(msh1 + 374);
    const auto *msh1_375 = buffer.data(msh1 + 375);
    const auto *msh1_376 = buffer.data(msh1 + 376);
    const auto *msh1_377 = buffer.data(msh1 + 377);
    const auto *msh1_378 = buffer.data(msh1 + 378);
    const auto *msh1_381 = buffer.data(msh1 + 381);
    const auto *msh1_383 = buffer.data(msh1 + 383);
    const auto *msh1_384 = buffer.data(msh1 + 384);
    const auto *msh1_387 = buffer.data(msh1 + 387);
    const auto *msh1_388 = buffer.data(msh1 + 388);
    const auto *msh1_390 = buffer.data(msh1 + 390);
    const auto *msh1_392 = buffer.data(msh1 + 392);
    const auto *msh1_393 = buffer.data(msh1 + 393);
    const auto *msh1_395 = buffer.data(msh1 + 395);
    const auto *msh1_396 = buffer.data(msh1 + 396);
    const auto *msh1_397 = buffer.data(msh1 + 397);
    const auto *msh1_398 = buffer.data(msh1 + 398);

    const auto *msi_469 = buffer.data(msi + 469);
    const auto *msi_471 = buffer.data(msi + 471);
    const auto *msi_472 = buffer.data(msi + 472);
    const auto *msi_473 = buffer.data(msi + 473);
    const auto *msi_474 = buffer.data(msi + 474);
    const auto *msi_475 = buffer.data(msi + 475);
    const auto *msi_476 = buffer.data(msi + 476);
    const auto *msi_478 = buffer.data(msi + 478);
    const auto *msi_479 = buffer.data(msi + 479);
    const auto *msi_481 = buffer.data(msi + 481);
    const auto *msi_482 = buffer.data(msi + 482);
    const auto *msi_485 = buffer.data(msi + 485);
    const auto *msi_486 = buffer.data(msi + 486);
    const auto *msi_488 = buffer.data(msi + 488);
    const auto *msi_490 = buffer.data(msi + 490);
    const auto *msi_491 = buffer.data(msi + 491);
    const auto *msi_493 = buffer.data(msi + 493);
    const auto *msi_494 = buffer.data(msi + 494);
    const auto *msi_496 = buffer.data(msi + 496);
    const auto *msi_497 = buffer.data(msi + 497);
    const auto *msi_498 = buffer.data(msi + 498);
    const auto *msi_499 = buffer.data(msi + 499);
    const auto *msi_500 = buffer.data(msi + 500);
    const auto *msi_501 = buffer.data(msi + 501);
    const auto *msi_502 = buffer.data(msi + 502);
    const auto *msi_503 = buffer.data(msi + 503);
    const auto *msi_504 = buffer.data(msi + 504);
    const auto *msi_506 = buffer.data(msi + 506);
    const auto *msi_507 = buffer.data(msi + 507);
    const auto *msi_509 = buffer.data(msi + 509);
    const auto *msi_510 = buffer.data(msi + 510);
    const auto *msi_513 = buffer.data(msi + 513);
    const auto *msi_514 = buffer.data(msi + 514);
    const auto *msi_516 = buffer.data(msi + 516);
    const auto *msi_518 = buffer.data(msi + 518);
    const auto *msi_519 = buffer.data(msi + 519);
    const auto *msi_521 = buffer.data(msi + 521);
    const auto *msi_522 = buffer.data(msi + 522);
    const auto *msi_524 = buffer.data(msi + 524);
    const auto *msi_525 = buffer.data(msi + 525);
    const auto *msi_526 = buffer.data(msi + 526);
    const auto *msi_527 = buffer.data(msi + 527);
    const auto *msi_528 = buffer.data(msi + 528);
    const auto *msi_529 = buffer.data(msi + 529);
    const auto *msi_530 = buffer.data(msi + 530);
    const auto *msi_531 = buffer.data(msi + 531);
    const auto *msi_532 = buffer.data(msi + 532);
    const auto *msi_534 = buffer.data(msi + 534);
    const auto *msi_535 = buffer.data(msi + 535);
    const auto *msi_537 = buffer.data(msi + 537);
    const auto *msi_538 = buffer.data(msi + 538);
    const auto *msi_541 = buffer.data(msi + 541);
    const auto *msi_542 = buffer.data(msi + 542);
    const auto *msi_546 = buffer.data(msi + 546);
    const auto *msi_553 = buffer.data(msi + 553);
    const auto *msi_554 = buffer.data(msi + 554);
    const auto *msi_555 = buffer.data(msi + 555);
    const auto *msi_556 = buffer.data(msi + 556);
    const auto *msi_557 = buffer.data(msi + 557);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, lsi_472, lsi_473, lsi_474, lsi_475, \
                         msi_472, msi_473, msi_474, msi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_16 * lsi_472[k]
                   + f_3 * pc_x[k] * msi_472[k];

        t_601[k] = f_16 * lsi_473[k]
                   + f_3 * pc_x[k] * msi_473[k];

        t_602[k] = f_16 * lsi_474[k]
                   + f_3 * pc_x[k] * msi_474[k];

        t_603[k] = f_16 * lsi_475[k]
                   + f_3 * pc_x[k] * msi_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pa_z, pc_y, pc_z, lsk0_388, lsi_301, lsi_331, \
                         lsk1_388, msh0_353, msh1_353, msi_469, \
                         msi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_z[k] * lsk0_388[k]
                   - f_12 * pc_z[k] * lsk1_388[k];

        t_605[k] = f_13 * lsi_301[k]
                   + f_3 * pc_z[k] * msi_469[k];

        t_606[k] = f_16 * lsi_331[k]
                   + f_10 * msh0_353[k]
                   - f_11 * msh1_353[k]
                   + f_3 * pc_y[k] * msi_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, lsi_332, lsi_333, lsi_334, msh0_354, \
                         msh0_355, msh0_356, msh1_354, msh1_355, msh1_356, msi_472, msi_473, \
                         msi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * lsi_332[k]
                   + f_8 * msh0_354[k]
                   - f_9 * msh1_354[k]
                   + f_3 * pc_y[k] * msi_472[k];

        t_608[k] = f_16 * lsi_333[k]
                   + f_6 * msh0_355[k]
                   - f_7 * msh1_355[k]
                   + f_3 * pc_y[k] * msi_473[k];

        t_609[k] = f_16 * lsi_334[k]
                   + f_4 * msh0_356[k]
                   - f_5 * msh1_356[k]
                   + f_3 * pc_y[k] * msi_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, lsi_307, lsi_335, lsi_476, \
                         msh0_356, msh0_357, msh1_356, msh1_357, msi_475, \
                         msi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * lsi_335[k]
                   + f_3 * pc_y[k] * msi_475[k];

        t_611[k] = f_13 * lsi_307[k]
                   + f_1 * msh0_356[k]
                   - f_2 * msh1_356[k]
                   + f_3 * pc_z[k] * msi_475[k];

        t_612[k] = f_16 * lsi_476[k]
                   + f_1 * msh0_357[k]
                   - f_2 * msh1_357[k]
                   + f_3 * pc_x[k] * msi_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, lsi_308, lsi_336, \
                         lsi_338, lsi_479, msh0_360, msh1_360, msi_476, msi_478, \
                         msi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * lsi_336[k]
                   + f_3 * pc_y[k] * msi_476[k];

        t_614[k] = f_14 * lsi_308[k]
                   + f_3 * pc_z[k] * msi_476[k];

        t_615[k] = f_16 * lsi_479[k]
                   + f_10 * msh0_360[k]
                   - f_11 * msh1_360[k]
                   + f_3 * pc_x[k] * msi_479[k];

        t_616[k] = f_15 * lsi_338[k]
                   + f_3 * pc_y[k] * msi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, lsi_311, lsi_481, lsi_482, msh0_362, \
                         msh0_363, msh1_362, msh1_363, msi_479, msi_481, \
                         msi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_16 * lsi_481[k]
                   + f_10 * msh0_362[k]
                   - f_11 * msh1_362[k]
                   + f_3 * pc_x[k] * msi_481[k];

        t_618[k] = f_16 * lsi_482[k]
                   + f_8 * msh0_363[k]
                   - f_9 * msh1_363[k]
                   + f_3 * pc_x[k] * msi_482[k];

        t_619[k] = f_14 * lsi_311[k]
                   + f_3 * pc_z[k] * msi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, lsi_341, lsi_485, lsi_486, msh0_366, \
                         msh0_367, msh1_366, msh1_367, msi_481, msi_485, \
                         msi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * lsi_341[k]
                   + f_3 * pc_y[k] * msi_481[k];

        t_621[k] = f_16 * lsi_485[k]
                   + f_8 * msh0_366[k]
                   - f_9 * msh1_366[k]
                   + f_3 * pc_x[k] * msi_485[k];

        t_622[k] = f_16 * lsi_486[k]
                   + f_6 * msh0_367[k]
                   - f_7 * msh1_367[k]
                   + f_3 * pc_x[k] * msi_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, lsi_314, lsi_345, lsi_488, \
                         msh0_369, msh1_369, msi_482, msi_485, \
                         msi_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * lsi_314[k]
                   + f_3 * pc_z[k] * msi_482[k];

        t_624[k] = f_16 * lsi_488[k]
                   + f_6 * msh0_369[k]
                   - f_7 * msh1_369[k]
                   + f_3 * pc_x[k] * msi_488[k];

        t_625[k] = f_15 * lsi_345[k]
                   + f_3 * pc_y[k] * msi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, lsi_318, lsi_490, lsi_491, msh0_371, \
                         msh0_372, msh1_371, msh1_372, msi_486, msi_490, \
                         msi_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_16 * lsi_490[k]
                   + f_6 * msh0_371[k]
                   - f_7 * msh1_371[k]
                   + f_3 * pc_x[k] * msi_490[k];

        t_627[k] = f_16 * lsi_491[k]
                   + f_4 * msh0_372[k]
                   - f_5 * msh1_372[k]
                   + f_3 * pc_x[k] * msi_491[k];

        t_628[k] = f_14 * lsi_318[k]
                   + f_3 * pc_z[k] * msi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, lsi_350, lsi_493, lsi_494, msh0_374, \
                         msh0_375, msh1_374, msh1_375, msi_490, msi_493, \
                         msi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_16 * lsi_493[k]
                   + f_4 * msh0_374[k]
                   - f_5 * msh1_374[k]
                   + f_3 * pc_x[k] * msi_493[k];

        t_630[k] = f_16 * lsi_494[k]
                   + f_4 * msh0_375[k]
                   - f_5 * msh1_375[k]
                   + f_3 * pc_x[k] * msi_494[k];

        t_631[k] = f_15 * lsi_350[k]
                   + f_3 * pc_y[k] * msi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, lsi_496, lsi_497, lsi_498, lsi_499, \
                         msh0_377, msh1_377, msi_496, msi_497, msi_498, \
                         msi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_16 * lsi_496[k]
                   + f_4 * msh0_377[k]
                   - f_5 * msh1_377[k]
                   + f_3 * pc_x[k] * msi_496[k];

        t_633[k] = f_16 * lsi_497[k]
                   + f_3 * pc_x[k] * msi_497[k];

        t_634[k] = f_16 * lsi_498[k]
                   + f_3 * pc_x[k] * msi_498[k];

        t_635[k] = f_16 * lsi_499[k]
                   + f_3 * pc_x[k] * msi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, lsi_500, lsi_501, lsi_502, lsi_503, \
                         msi_500, msi_501, msi_502, msi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_16 * lsi_500[k]
                   + f_3 * pc_x[k] * msi_500[k];

        t_637[k] = f_16 * lsi_501[k]
                   + f_3 * pc_x[k] * msi_501[k];

        t_638[k] = f_16 * lsi_502[k]
                   + f_3 * pc_x[k] * msi_502[k];

        t_639[k] = f_16 * lsi_503[k]
                   + f_3 * pc_x[k] * msi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, lsi_329, lsi_357, lsi_359, msh0_372, \
                         msh0_374, msh1_372, msh1_374, msi_497, \
                         msi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * lsi_357[k]
                   + f_1 * msh0_372[k]
                   - f_2 * msh1_372[k]
                   + f_3 * pc_y[k] * msi_497[k];

        t_641[k] = f_14 * lsi_329[k]
                   + f_3 * pc_z[k] * msi_497[k];

        t_642[k] = f_15 * lsi_359[k]
                   + f_10 * msh0_374[k]
                   - f_11 * msh1_374[k]
                   + f_3 * pc_y[k] * msi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, lsi_360, lsi_361, lsi_362, msh0_375, \
                         msh0_376, msh0_377, msh1_375, msh1_376, msh1_377, msi_500, msi_501, \
                         msi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * lsi_360[k]
                   + f_8 * msh0_375[k]
                   - f_9 * msh1_375[k]
                   + f_3 * pc_y[k] * msi_500[k];

        t_644[k] = f_15 * lsi_361[k]
                   + f_6 * msh0_376[k]
                   - f_7 * msh1_376[k]
                   + f_3 * pc_y[k] * msi_501[k];

        t_645[k] = f_15 * lsi_362[k]
                   + f_4 * msh0_377[k]
                   - f_5 * msh1_377[k]
                   + f_3 * pc_y[k] * msi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, lsi_335, lsi_363, lsi_504, \
                         msh0_377, msh0_378, msh1_377, msh1_378, msi_503, \
                         msi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * lsi_363[k]
                   + f_3 * pc_y[k] * msi_503[k];

        t_647[k] = f_14 * lsi_335[k]
                   + f_1 * msh0_377[k]
                   - f_2 * msh1_377[k]
                   + f_3 * pc_z[k] * msi_503[k];

        t_648[k] = f_16 * lsi_504[k]
                   + f_1 * msh0_378[k]
                   - f_2 * msh1_378[k]
                   + f_3 * pc_x[k] * msi_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, lsi_336, lsi_364, \
                         lsi_366, lsi_507, msh0_381, msh1_381, msi_504, msi_506, \
                         msi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * lsi_364[k]
                   + f_3 * pc_y[k] * msi_504[k];

        t_650[k] = f_15 * lsi_336[k]
                   + f_3 * pc_z[k] * msi_504[k];

        t_651[k] = f_16 * lsi_507[k]
                   + f_10 * msh0_381[k]
                   - f_11 * msh1_381[k]
                   + f_3 * pc_x[k] * msi_507[k];

        t_652[k] = f_14 * lsi_366[k]
                   + f_3 * pc_y[k] * msi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, lsi_339, lsi_509, lsi_510, msh0_383, \
                         msh0_384, msh1_383, msh1_384, msi_507, msi_509, \
                         msi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_16 * lsi_509[k]
                   + f_10 * msh0_383[k]
                   - f_11 * msh1_383[k]
                   + f_3 * pc_x[k] * msi_509[k];

        t_654[k] = f_16 * lsi_510[k]
                   + f_8 * msh0_384[k]
                   - f_9 * msh1_384[k]
                   + f_3 * pc_x[k] * msi_510[k];

        t_655[k] = f_15 * lsi_339[k]
                   + f_3 * pc_z[k] * msi_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, lsi_369, lsi_513, lsi_514, msh0_387, \
                         msh0_388, msh1_387, msh1_388, msi_509, msi_513, \
                         msi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * lsi_369[k]
                   + f_3 * pc_y[k] * msi_509[k];

        t_657[k] = f_16 * lsi_513[k]
                   + f_8 * msh0_387[k]
                   - f_9 * msh1_387[k]
                   + f_3 * pc_x[k] * msi_513[k];

        t_658[k] = f_16 * lsi_514[k]
                   + f_6 * msh0_388[k]
                   - f_7 * msh1_388[k]
                   + f_3 * pc_x[k] * msi_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, lsi_342, lsi_373, lsi_516, \
                         msh0_390, msh1_390, msi_510, msi_513, \
                         msi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * lsi_342[k]
                   + f_3 * pc_z[k] * msi_510[k];

        t_660[k] = f_16 * lsi_516[k]
                   + f_6 * msh0_390[k]
                   - f_7 * msh1_390[k]
                   + f_3 * pc_x[k] * msi_516[k];

        t_661[k] = f_14 * lsi_373[k]
                   + f_3 * pc_y[k] * msi_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, lsi_346, lsi_518, lsi_519, msh0_392, \
                         msh0_393, msh1_392, msh1_393, msi_514, msi_518, \
                         msi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_16 * lsi_518[k]
                   + f_6 * msh0_392[k]
                   - f_7 * msh1_392[k]
                   + f_3 * pc_x[k] * msi_518[k];

        t_663[k] = f_16 * lsi_519[k]
                   + f_4 * msh0_393[k]
                   - f_5 * msh1_393[k]
                   + f_3 * pc_x[k] * msi_519[k];

        t_664[k] = f_15 * lsi_346[k]
                   + f_3 * pc_z[k] * msi_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, lsi_378, lsi_521, lsi_522, msh0_395, \
                         msh0_396, msh1_395, msh1_396, msi_518, msi_521, \
                         msi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_16 * lsi_521[k]
                   + f_4 * msh0_395[k]
                   - f_5 * msh1_395[k]
                   + f_3 * pc_x[k] * msi_521[k];

        t_666[k] = f_16 * lsi_522[k]
                   + f_4 * msh0_396[k]
                   - f_5 * msh1_396[k]
                   + f_3 * pc_x[k] * msi_522[k];

        t_667[k] = f_14 * lsi_378[k]
                   + f_3 * pc_y[k] * msi_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, lsi_524, lsi_525, lsi_526, lsi_527, \
                         msh0_398, msh1_398, msi_524, msi_525, msi_526, \
                         msi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_16 * lsi_524[k]
                   + f_4 * msh0_398[k]
                   - f_5 * msh1_398[k]
                   + f_3 * pc_x[k] * msi_524[k];

        t_669[k] = f_16 * lsi_525[k]
                   + f_3 * pc_x[k] * msi_525[k];

        t_670[k] = f_16 * lsi_526[k]
                   + f_3 * pc_x[k] * msi_526[k];

        t_671[k] = f_16 * lsi_527[k]
                   + f_3 * pc_x[k] * msi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, lsi_528, lsi_529, lsi_530, lsi_531, \
                         msi_528, msi_529, msi_530, msi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_16 * lsi_528[k]
                   + f_3 * pc_x[k] * msi_528[k];

        t_673[k] = f_16 * lsi_529[k]
                   + f_3 * pc_x[k] * msi_529[k];

        t_674[k] = f_16 * lsi_530[k]
                   + f_3 * pc_x[k] * msi_530[k];

        t_675[k] = f_16 * lsi_531[k]
                   + f_3 * pc_x[k] * msi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, lsi_357, lsi_385, lsi_387, msh0_393, \
                         msh0_395, msh1_393, msh1_395, msi_525, \
                         msi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * lsi_385[k]
                   + f_1 * msh0_393[k]
                   - f_2 * msh1_393[k]
                   + f_3 * pc_y[k] * msi_525[k];

        t_677[k] = f_15 * lsi_357[k]
                   + f_3 * pc_z[k] * msi_525[k];

        t_678[k] = f_14 * lsi_387[k]
                   + f_10 * msh0_395[k]
                   - f_11 * msh1_395[k]
                   + f_3 * pc_y[k] * msi_527[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, lsi_388, lsi_389, lsi_390, msh0_396, \
                         msh0_397, msh0_398, msh1_396, msh1_397, msh1_398, msi_528, msi_529, \
                         msi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * lsi_388[k]
                   + f_8 * msh0_396[k]
                   - f_9 * msh1_396[k]
                   + f_3 * pc_y[k] * msi_528[k];

        t_680[k] = f_14 * lsi_389[k]
                   + f_6 * msh0_397[k]
                   - f_7 * msh1_397[k]
                   + f_3 * pc_y[k] * msi_529[k];

        t_681[k] = f_14 * lsi_390[k]
                   + f_4 * msh0_398[k]
                   - f_5 * msh1_398[k]
                   + f_3 * pc_y[k] * msi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_y, pc_y, pc_z, lsk0_504, lsi_363, \
                         lsi_391, lsi_392, lsk1_504, msh0_398, msh1_398, msi_531, \
                         msi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * lsi_391[k]
                   + f_3 * pc_y[k] * msi_531[k];

        t_683[k] = f_15 * lsi_363[k]
                   + f_1 * msh0_398[k]
                   - f_2 * msh1_398[k]
                   + f_3 * pc_z[k] * msi_531[k];

        t_684[k] = pa_y[k] * lsk0_504[k]
                   - f_12 * pc_y[k] * lsk1_504[k];

        t_685[k] = f_13 * lsi_392[k]
                   + f_3 * pc_y[k] * msi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_y, pc_y, pc_z, lsk0_507, lsk0_509, \
                         lsi_364, lsi_393, lsi_394, lsk1_507, lsk1_509, msi_532, \
                         msi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * lsi_364[k]
                   + f_3 * pc_z[k] * msi_532[k];

        t_687[k] = pa_y[k] * lsk0_507[k]
                   + f_14 * lsi_393[k]
                   - f_12 * pc_y[k] * lsk1_507[k];

        t_688[k] = f_13 * lsi_394[k]
                   + f_3 * pc_y[k] * msi_534[k];

        t_689[k] = pa_y[k] * lsk0_509[k]
                   - f_12 * pc_y[k] * lsk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_y, pc_y, pc_z, lsk0_510, lsk0_513, \
                         lsi_367, lsi_395, lsi_397, lsk1_510, lsk1_513, msi_535, \
                         msi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_y[k] * lsk0_510[k]
                   + f_15 * lsi_395[k]
                   - f_12 * pc_y[k] * lsk1_510[k];

        t_691[k] = f_16 * lsi_367[k]
                   + f_3 * pc_z[k] * msi_535[k];

        t_692[k] = f_13 * lsi_397[k]
                   + f_3 * pc_y[k] * msi_537[k];

        t_693[k] = pa_y[k] * lsk0_513[k]
                   - f_12 * pc_y[k] * lsk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pa_y, pc_y, pc_z, lsk0_514, lsk0_516, lsi_370, \
                         lsi_398, lsi_400, lsk1_514, lsk1_516, \
                         msi_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pa_y[k] * lsk0_514[k]
                   + f_16 * lsi_398[k]
                   - f_12 * pc_y[k] * lsk1_514[k];

        t_695[k] = f_16 * lsi_370[k]
                   + f_3 * pc_z[k] * msi_538[k];

        t_696[k] = pa_y[k] * lsk0_516[k]
                   + f_14 * lsi_400[k]
                   - f_12 * pc_y[k] * lsk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pa_y, pc_y, pc_z, lsk0_518, lsk0_519, \
                         lsi_374, lsi_401, lsi_402, lsk1_518, lsk1_519, msi_541, \
                         msi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * lsi_401[k]
                   + f_3 * pc_y[k] * msi_541[k];

        t_698[k] = pa_y[k] * lsk0_518[k]
                   - f_12 * pc_y[k] * lsk1_518[k];

        t_699[k] = pa_y[k] * lsk0_519[k]
                   + f_17 * lsi_402[k]
                   - f_12 * pc_y[k] * lsk1_519[k];

        t_700[k] = f_16 * lsi_374[k]
                   + f_3 * pc_z[k] * msi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pa_y, pc_y, lsk0_521, lsk0_522, lsk0_524, \
                         lsi_404, lsi_405, lsi_406, lsk1_521, lsk1_522, lsk1_524, \
                         msi_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pa_y[k] * lsk0_521[k]
                   + f_15 * lsi_404[k]
                   - f_12 * pc_y[k] * lsk1_521[k];

        t_702[k] = pa_y[k] * lsk0_522[k]
                   + f_14 * lsi_405[k]
                   - f_12 * pc_y[k] * lsk1_522[k];

        t_703[k] = f_13 * lsi_406[k]
                   + f_3 * pc_y[k] * msi_546[k];

        t_704[k] = pa_y[k] * lsk0_524[k]
                   - f_12 * pc_y[k] * lsk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, lsi_553, lsi_554, lsi_555, \
                         lsi_556, lsi_557, msi_553, msi_554, msi_555, msi_556, \
                         msi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_16 * lsi_553[k]
                   + f_3 * pc_x[k] * msi_553[k];

        t_706[k] = f_16 * lsi_554[k]
                   + f_3 * pc_x[k] * msi_554[k];

        t_707[k] = f_16 * lsi_555[k]
                   + f_3 * pc_x[k] * msi_555[k];

        t_708[k] = f_16 * lsi_556[k]
                   + f_3 * pc_x[k] * msi_556[k];

        t_709[k] = f_16 * lsi_557[k]
                   + f_3 * pc_x[k] * msi_557[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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
    const auto f_22 = 3.0 / q;

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

    const auto *lsk0_539 = buffer.data(lsk0 + 539);
    const auto *lsk0_540 = buffer.data(lsk0 + 540);
    const auto *lsk0_543 = buffer.data(lsk0 + 543);
    const auto *lsk0_546 = buffer.data(lsk0 + 546);
    const auto *lsk0_550 = buffer.data(lsk0 + 550);
    const auto *lsk0_552 = buffer.data(lsk0 + 552);
    const auto *lsk0_555 = buffer.data(lsk0 + 555);
    const auto *lsk0_557 = buffer.data(lsk0 + 557);
    const auto *lsk0_558 = buffer.data(lsk0 + 558);
    const auto *lsk0_568 = buffer.data(lsk0 + 568);

    const auto *lsi_385 = buffer.data(lsi + 385);
    const auto *lsi_392 = buffer.data(lsi + 392);
    const auto *lsi_413 = buffer.data(lsi + 413);
    const auto *lsi_415 = buffer.data(lsi + 415);
    const auto *lsi_416 = buffer.data(lsi + 416);
    const auto *lsi_417 = buffer.data(lsi + 417);
    const auto *lsi_418 = buffer.data(lsi + 418);
    const auto *lsi_419 = buffer.data(lsi + 419);
    const auto *lsi_420 = buffer.data(lsi + 420);
    const auto *lsi_423 = buffer.data(lsi + 423);
    const auto *lsi_425 = buffer.data(lsi + 425);
    const auto *lsi_426 = buffer.data(lsi + 426);
    const auto *lsi_427 = buffer.data(lsi + 427);
    const auto *lsi_429 = buffer.data(lsi + 429);
    const auto *lsi_430 = buffer.data(lsi + 430);
    const auto *lsi_431 = buffer.data(lsi + 431);
    const auto *lsi_432 = buffer.data(lsi + 432);
    const auto *lsi_434 = buffer.data(lsi + 434);
    const auto *lsi_441 = buffer.data(lsi + 441);
    const auto *lsi_447 = buffer.data(lsi + 447);
    const auto *lsi_448 = buffer.data(lsi + 448);
    const auto *lsi_450 = buffer.data(lsi + 450);
    const auto *lsi_453 = buffer.data(lsi + 453);
    const auto *lsi_457 = buffer.data(lsi + 457);
    const auto *lsi_462 = buffer.data(lsi + 462);
    const auto *lsi_471 = buffer.data(lsi + 471);
    const auto *lsi_472 = buffer.data(lsi + 472);
    const auto *lsi_473 = buffer.data(lsi + 473);
    const auto *lsi_474 = buffer.data(lsi + 474);
    const auto *lsi_558 = buffer.data(lsi + 558);
    const auto *lsi_559 = buffer.data(lsi + 559);
    const auto *lsi_560 = buffer.data(lsi + 560);
    const auto *lsi_565 = buffer.data(lsi + 565);
    const auto *lsi_569 = buffer.data(lsi + 569);
    const auto *lsi_574 = buffer.data(lsi + 574);
    const auto *lsi_580 = buffer.data(lsi + 580);
    const auto *lsi_581 = buffer.data(lsi + 581);
    const auto *lsi_582 = buffer.data(lsi + 582);
    const auto *lsi_583 = buffer.data(lsi + 583);
    const auto *lsi_584 = buffer.data(lsi + 584);
    const auto *lsi_585 = buffer.data(lsi + 585);
    const auto *lsi_587 = buffer.data(lsi + 587);
    const auto *lsi_588 = buffer.data(lsi + 588);
    const auto *lsi_591 = buffer.data(lsi + 591);
    const auto *lsi_594 = buffer.data(lsi + 594);
    const auto *lsi_598 = buffer.data(lsi + 598);
    const auto *lsi_603 = buffer.data(lsi + 603);
    const auto *lsi_609 = buffer.data(lsi + 609);
    const auto *lsi_611 = buffer.data(lsi + 611);
    const auto *lsi_612 = buffer.data(lsi + 612);
    const auto *lsi_613 = buffer.data(lsi + 613);
    const auto *lsi_614 = buffer.data(lsi + 614);
    const auto *lsi_615 = buffer.data(lsi + 615);
    const auto *lsi_621 = buffer.data(lsi + 621);
    const auto *lsi_625 = buffer.data(lsi + 625);
    const auto *lsi_630 = buffer.data(lsi + 630);
    const auto *lsi_636 = buffer.data(lsi + 636);
    const auto *lsi_637 = buffer.data(lsi + 637);
    const auto *lsi_638 = buffer.data(lsi + 638);
    const auto *lsi_639 = buffer.data(lsi + 639);
    const auto *lsi_640 = buffer.data(lsi + 640);
    const auto *lsi_641 = buffer.data(lsi + 641);
    const auto *lsi_642 = buffer.data(lsi + 642);
    const auto *lsi_643 = buffer.data(lsi + 643);

    const auto *lsk1_539 = buffer.data(lsk1 + 539);
    const auto *lsk1_540 = buffer.data(lsk1 + 540);
    const auto *lsk1_543 = buffer.data(lsk1 + 543);
    const auto *lsk1_546 = buffer.data(lsk1 + 546);
    const auto *lsk1_550 = buffer.data(lsk1 + 550);
    const auto *lsk1_552 = buffer.data(lsk1 + 552);
    const auto *lsk1_555 = buffer.data(lsk1 + 555);
    const auto *lsk1_557 = buffer.data(lsk1 + 557);
    const auto *lsk1_558 = buffer.data(lsk1 + 558);
    const auto *lsk1_568 = buffer.data(lsk1 + 568);

    const auto *msh0_414 = buffer.data(msh0 + 414);
    const auto *msh0_416 = buffer.data(msh0 + 416);
    const auto *msh0_417 = buffer.data(msh0 + 417);
    const auto *msh0_418 = buffer.data(msh0 + 418);
    const auto *msh0_419 = buffer.data(msh0 + 419);
    const auto *msh0_420 = buffer.data(msh0 + 420);
    const auto *msh0_421 = buffer.data(msh0 + 421);
    const auto *msh0_422 = buffer.data(msh0 + 422);
    const auto *msh0_423 = buffer.data(msh0 + 423);
    const auto *msh0_424 = buffer.data(msh0 + 424);
    const auto *msh0_425 = buffer.data(msh0 + 425);
    const auto *msh0_426 = buffer.data(msh0 + 426);
    const auto *msh0_427 = buffer.data(msh0 + 427);
    const auto *msh0_428 = buffer.data(msh0 + 428);
    const auto *msh0_429 = buffer.data(msh0 + 429);
    const auto *msh0_434 = buffer.data(msh0 + 434);
    const auto *msh0_435 = buffer.data(msh0 + 435);
    const auto *msh0_436 = buffer.data(msh0 + 436);
    const auto *msh0_437 = buffer.data(msh0 + 437);
    const auto *msh0_438 = buffer.data(msh0 + 438);
    const auto *msh0_439 = buffer.data(msh0 + 439);
    const auto *msh0_440 = buffer.data(msh0 + 440);
    const auto *msh0_441 = buffer.data(msh0 + 441);
    const auto *msh0_443 = buffer.data(msh0 + 443);
    const auto *msh0_444 = buffer.data(msh0 + 444);
    const auto *msh0_446 = buffer.data(msh0 + 446);
    const auto *msh0_447 = buffer.data(msh0 + 447);
    const auto *msh0_448 = buffer.data(msh0 + 448);
    const auto *msh0_450 = buffer.data(msh0 + 450);
    const auto *msh0_451 = buffer.data(msh0 + 451);
    const auto *msh0_456 = buffer.data(msh0 + 456);
    const auto *msh0_457 = buffer.data(msh0 + 457);
    const auto *msh0_458 = buffer.data(msh0 + 458);
    const auto *msh0_459 = buffer.data(msh0 + 459);
    const auto *msh0_461 = buffer.data(msh0 + 461);
    const auto *msh0_467 = buffer.data(msh0 + 467);
    const auto *msh0_471 = buffer.data(msh0 + 471);
    const auto *msh0_476 = buffer.data(msh0 + 476);
    const auto *msh0_479 = buffer.data(msh0 + 479);
    const auto *msh0_480 = buffer.data(msh0 + 480);
    const auto *msh0_481 = buffer.data(msh0 + 481);
    const auto *msh0_482 = buffer.data(msh0 + 482);

    const auto *msh1_414 = buffer.data(msh1 + 414);
    const auto *msh1_416 = buffer.data(msh1 + 416);
    const auto *msh1_417 = buffer.data(msh1 + 417);
    const auto *msh1_418 = buffer.data(msh1 + 418);
    const auto *msh1_419 = buffer.data(msh1 + 419);
    const auto *msh1_420 = buffer.data(msh1 + 420);
    const auto *msh1_421 = buffer.data(msh1 + 421);
    const auto *msh1_422 = buffer.data(msh1 + 422);
    const auto *msh1_423 = buffer.data(msh1 + 423);
    const auto *msh1_424 = buffer.data(msh1 + 424);
    const auto *msh1_425 = buffer.data(msh1 + 425);
    const auto *msh1_426 = buffer.data(msh1 + 426);
    const auto *msh1_427 = buffer.data(msh1 + 427);
    const auto *msh1_428 = buffer.data(msh1 + 428);
    const auto *msh1_429 = buffer.data(msh1 + 429);
    const auto *msh1_434 = buffer.data(msh1 + 434);
    const auto *msh1_435 = buffer.data(msh1 + 435);
    const auto *msh1_436 = buffer.data(msh1 + 436);
    const auto *msh1_437 = buffer.data(msh1 + 437);
    const auto *msh1_438 = buffer.data(msh1 + 438);
    const auto *msh1_439 = buffer.data(msh1 + 439);
    const auto *msh1_440 = buffer.data(msh1 + 440);
    const auto *msh1_441 = buffer.data(msh1 + 441);
    const auto *msh1_443 = buffer.data(msh1 + 443);
    const auto *msh1_444 = buffer.data(msh1 + 444);
    const auto *msh1_446 = buffer.data(msh1 + 446);
    const auto *msh1_447 = buffer.data(msh1 + 447);
    const auto *msh1_448 = buffer.data(msh1 + 448);
    const auto *msh1_450 = buffer.data(msh1 + 450);
    const auto *msh1_451 = buffer.data(msh1 + 451);
    const auto *msh1_456 = buffer.data(msh1 + 456);
    const auto *msh1_457 = buffer.data(msh1 + 457);
    const auto *msh1_458 = buffer.data(msh1 + 458);
    const auto *msh1_459 = buffer.data(msh1 + 459);
    const auto *msh1_461 = buffer.data(msh1 + 461);
    const auto *msh1_467 = buffer.data(msh1 + 467);
    const auto *msh1_471 = buffer.data(msh1 + 471);
    const auto *msh1_476 = buffer.data(msh1 + 476);
    const auto *msh1_479 = buffer.data(msh1 + 479);
    const auto *msh1_480 = buffer.data(msh1 + 480);
    const auto *msh1_481 = buffer.data(msh1 + 481);
    const auto *msh1_482 = buffer.data(msh1 + 482);

    const auto *msi_553 = buffer.data(msi + 553);
    const auto *msi_555 = buffer.data(msi + 555);
    const auto *msi_556 = buffer.data(msi + 556);
    const auto *msi_557 = buffer.data(msi + 557);
    const auto *msi_558 = buffer.data(msi + 558);
    const auto *msi_559 = buffer.data(msi + 559);
    const auto *msi_560 = buffer.data(msi + 560);
    const auto *msi_561 = buffer.data(msi + 561);
    const auto *msi_562 = buffer.data(msi + 562);
    const auto *msi_563 = buffer.data(msi + 563);
    const auto *msi_564 = buffer.data(msi + 564);
    const auto *msi_565 = buffer.data(msi + 565);
    const auto *msi_566 = buffer.data(msi + 566);
    const auto *msi_567 = buffer.data(msi + 567);
    const auto *msi_568 = buffer.data(msi + 568);
    const auto *msi_569 = buffer.data(msi + 569);
    const auto *msi_570 = buffer.data(msi + 570);
    const auto *msi_571 = buffer.data(msi + 571);
    const auto *msi_572 = buffer.data(msi + 572);
    const auto *msi_573 = buffer.data(msi + 573);
    const auto *msi_574 = buffer.data(msi + 574);
    const auto *msi_580 = buffer.data(msi + 580);
    const auto *msi_581 = buffer.data(msi + 581);
    const auto *msi_582 = buffer.data(msi + 582);
    const auto *msi_583 = buffer.data(msi + 583);
    const auto *msi_584 = buffer.data(msi + 584);
    const auto *msi_585 = buffer.data(msi + 585);
    const auto *msi_586 = buffer.data(msi + 586);
    const auto *msi_587 = buffer.data(msi + 587);
    const auto *msi_588 = buffer.data(msi + 588);
    const auto *msi_589 = buffer.data(msi + 589);
    const auto *msi_590 = buffer.data(msi + 590);
    const auto *msi_591 = buffer.data(msi + 591);
    const auto *msi_593 = buffer.data(msi + 593);
    const auto *msi_594 = buffer.data(msi + 594);
    const auto *msi_595 = buffer.data(msi + 595);
    const auto *msi_597 = buffer.data(msi + 597);
    const auto *msi_598 = buffer.data(msi + 598);
    const auto *msi_599 = buffer.data(msi + 599);
    const auto *msi_600 = buffer.data(msi + 600);
    const auto *msi_602 = buffer.data(msi + 602);
    const auto *msi_603 = buffer.data(msi + 603);
    const auto *msi_609 = buffer.data(msi + 609);
    const auto *msi_610 = buffer.data(msi + 610);
    const auto *msi_611 = buffer.data(msi + 611);
    const auto *msi_612 = buffer.data(msi + 612);
    const auto *msi_613 = buffer.data(msi + 613);
    const auto *msi_614 = buffer.data(msi + 614);
    const auto *msi_615 = buffer.data(msi + 615);
    const auto *msi_616 = buffer.data(msi + 616);
    const auto *msi_618 = buffer.data(msi + 618);
    const auto *msi_619 = buffer.data(msi + 619);
    const auto *msi_621 = buffer.data(msi + 621);
    const auto *msi_622 = buffer.data(msi + 622);
    const auto *msi_625 = buffer.data(msi + 625);
    const auto *msi_626 = buffer.data(msi + 626);
    const auto *msi_630 = buffer.data(msi + 630);
    const auto *msi_636 = buffer.data(msi + 636);
    const auto *msi_637 = buffer.data(msi + 637);
    const auto *msi_638 = buffer.data(msi + 638);
    const auto *msi_639 = buffer.data(msi + 639);
    const auto *msi_640 = buffer.data(msi + 640);
    const auto *msi_641 = buffer.data(msi + 641);
    const auto *msi_642 = buffer.data(msi + 642);
    const auto *msi_643 = buffer.data(msi + 643);

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, lsi_385, lsi_413, \
                         lsi_558, lsi_559, msh0_414, msh1_414, msi_553, msi_558, \
                         msi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_16 * lsi_558[k]
                   + f_3 * pc_x[k] * msi_558[k];

        t_711[k] = f_16 * lsi_559[k]
                   + f_3 * pc_x[k] * msi_559[k];

        t_712[k] = f_13 * lsi_413[k]
                   + f_1 * msh0_414[k]
                   - f_2 * msh1_414[k]
                   + f_3 * pc_y[k] * msi_553[k];

        t_713[k] = f_16 * lsi_385[k]
                   + f_3 * pc_z[k] * msi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, lsi_415, lsi_416, lsi_417, msh0_416, \
                         msh0_417, msh0_418, msh1_416, msh1_417, msh1_418, msi_555, msi_556, \
                         msi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * lsi_415[k]
                   + f_10 * msh0_416[k]
                   - f_11 * msh1_416[k]
                   + f_3 * pc_y[k] * msi_555[k];

        t_715[k] = f_13 * lsi_416[k]
                   + f_8 * msh0_417[k]
                   - f_9 * msh1_417[k]
                   + f_3 * pc_y[k] * msi_556[k];

        t_716[k] = f_13 * lsi_417[k]
                   + f_6 * msh0_418[k]
                   - f_7 * msh1_418[k]
                   + f_3 * pc_y[k] * msi_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_y, pc_y, lsk0_539, lsi_418, lsi_419, \
                         lsk1_539, msh0_419, msh1_419, msi_558, \
                         msi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * lsi_418[k]
                   + f_4 * msh0_419[k]
                   - f_5 * msh1_419[k]
                   + f_3 * pc_y[k] * msi_558[k];

        t_718[k] = f_13 * lsi_419[k]
                   + f_3 * pc_y[k] * msi_559[k];

        t_719[k] = pa_y[k] * lsk0_539[k]
                   - f_12 * pc_y[k] * lsk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pc_x, pc_y, pc_z, lsi_392, \
                         lsi_560, msh0_420, msh1_420, msi_560, msi_561, \
                         msi_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_16 * lsi_560[k]
                   + f_1 * msh0_420[k]
                   - f_2 * msh1_420[k]
                   + f_3 * pc_x[k] * msi_560[k];

        t_721[k] = f_3 * pc_y[k] * msi_560[k];

        t_722[k] = f_17 * lsi_392[k]
                   + f_3 * pc_z[k] * msi_560[k];

        t_723[k] = f_4 * msh0_420[k]
                   - f_5 * msh1_420[k]
                   + f_3 * pc_y[k] * msi_561[k];

        t_724[k] = f_3 * pc_y[k] * msi_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, pc_y, lsi_565, msh0_421, msh0_422, \
                         msh0_425, msh1_421, msh1_422, msh1_425, msi_563, msi_564, \
                         msi_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_16 * lsi_565[k]
                   + f_10 * msh0_425[k]
                   - f_11 * msh1_425[k]
                   + f_3 * pc_x[k] * msi_565[k];

        t_726[k] = f_6 * msh0_421[k]
                   - f_7 * msh1_421[k]
                   + f_3 * pc_y[k] * msi_563[k];

        t_727[k] = f_4 * msh0_422[k]
                   - f_5 * msh1_422[k]
                   + f_3 * pc_y[k] * msi_564[k];

        t_728[k] = f_3 * pc_y[k] * msi_565[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, lsi_569, msh0_423, msh0_424, \
                         msh0_429, msh1_423, msh1_424, msh1_429, msi_566, msi_567, \
                         msi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_16 * lsi_569[k]
                   + f_8 * msh0_429[k]
                   - f_9 * msh1_429[k]
                   + f_3 * pc_x[k] * msi_569[k];

        t_730[k] = f_8 * msh0_423[k]
                   - f_9 * msh1_423[k]
                   + f_3 * pc_y[k] * msi_566[k];

        t_731[k] = f_6 * msh0_424[k]
                   - f_7 * msh1_424[k]
                   + f_3 * pc_y[k] * msi_567[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_x, pc_y, lsi_574, msh0_425, msh0_434, \
                         msh1_425, msh1_434, msi_568, msi_569, \
                         msi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_4 * msh0_425[k]
                   - f_5 * msh1_425[k]
                   + f_3 * pc_y[k] * msi_568[k];

        t_733[k] = f_3 * pc_y[k] * msi_569[k];

        t_734[k] = f_16 * lsi_574[k]
                   + f_6 * msh0_434[k]
                   - f_7 * msh1_434[k]
                   + f_3 * pc_x[k] * msi_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pc_y, msh0_426, msh0_427, msh0_428, msh1_426, \
                         msh1_427, msh1_428, msi_570, msi_571, \
                         msi_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_10 * msh0_426[k]
                   - f_11 * msh1_426[k]
                   + f_3 * pc_y[k] * msi_570[k];

        t_736[k] = f_8 * msh0_427[k]
                   - f_9 * msh1_427[k]
                   + f_3 * pc_y[k] * msi_571[k];

        t_737[k] = f_6 * msh0_428[k]
                   - f_7 * msh1_428[k]
                   + f_3 * pc_y[k] * msi_572[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, pc_x, pc_y, lsi_580, lsi_581, msh0_429, \
                         msh0_440, msh1_429, msh1_440, msi_573, msi_574, msi_580, \
                         msi_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_4 * msh0_429[k]
                   - f_5 * msh1_429[k]
                   + f_3 * pc_y[k] * msi_573[k];

        t_739[k] = f_3 * pc_y[k] * msi_574[k];

        t_740[k] = f_16 * lsi_580[k]
                   + f_4 * msh0_440[k]
                   - f_5 * msh1_440[k]
                   + f_3 * pc_x[k] * msi_580[k];

        t_741[k] = f_16 * lsi_581[k]
                   + f_3 * pc_x[k] * msi_581[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, pc_x, pc_y, lsi_582, lsi_583, \
                         lsi_584, lsi_585, msi_580, msi_582, msi_583, msi_584, \
                         msi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_16 * lsi_582[k]
                   + f_3 * pc_x[k] * msi_582[k];

        t_743[k] = f_16 * lsi_583[k]
                   + f_3 * pc_x[k] * msi_583[k];

        t_744[k] = f_16 * lsi_584[k]
                   + f_3 * pc_x[k] * msi_584[k];

        t_745[k] = f_16 * lsi_585[k]
                   + f_3 * pc_x[k] * msi_585[k];

        t_746[k] = f_3 * pc_y[k] * msi_580[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_x, pc_y, lsi_587, msh0_435, msh0_436, \
                         msh1_435, msh1_436, msi_581, msi_582, \
                         msi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_16 * lsi_587[k]
                   + f_3 * pc_x[k] * msi_587[k];

        t_748[k] = f_1 * msh0_435[k]
                   - f_2 * msh1_435[k]
                   + f_3 * pc_y[k] * msi_581[k];

        t_749[k] = f_19 * msh0_436[k]
                   - f_20 * msh1_436[k]
                   + f_3 * pc_y[k] * msi_582[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, msh0_437, msh0_438, msh0_439, msh1_437, \
                         msh1_438, msh1_439, msi_583, msi_584, \
                         msi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_10 * msh0_437[k]
                   - f_11 * msh1_437[k]
                   + f_3 * pc_y[k] * msi_583[k];

        t_751[k] = f_8 * msh0_438[k]
                   - f_9 * msh1_438[k]
                   + f_3 * pc_y[k] * msi_584[k];

        t_752[k] = f_6 * msh0_439[k]
                   - f_7 * msh1_439[k]
                   + f_3 * pc_y[k] * msi_585[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pc_x, pc_y, pc_z, lsi_419, lsi_588, \
                         msh0_440, msh0_441, msh1_440, msh1_441, msi_586, msi_587, \
                         msi_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * msh0_440[k]
                   - f_5 * msh1_440[k]
                   + f_3 * pc_y[k] * msi_586[k];

        t_754[k] = f_3 * pc_y[k] * msi_587[k];

        t_755[k] = f_17 * lsi_419[k]
                   + f_1 * msh0_440[k]
                   - f_2 * msh1_440[k]
                   + f_3 * pc_z[k] * msi_587[k];

        t_756[k] = f_15 * lsi_588[k]
                   + f_1 * msh0_441[k]
                   - f_2 * msh1_441[k]
                   + f_3 * pc_x[k] * msi_588[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pc_x, pc_y, pc_z, lsi_420, lsi_591, \
                         msh0_444, msh1_444, msi_588, msi_589, \
                         msi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_22 * lsi_420[k]
                   + f_3 * pc_y[k] * msi_588[k];

        t_758[k] = f_3 * pc_z[k] * msi_588[k];

        t_759[k] = f_15 * lsi_591[k]
                   + f_10 * msh0_444[k]
                   - f_11 * msh1_444[k]
                   + f_3 * pc_x[k] * msi_591[k];

        t_760[k] = f_3 * pc_z[k] * msi_589[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_x, pc_z, lsi_594, msh0_441, msh0_447, \
                         msh1_441, msh1_447, msi_590, msi_591, \
                         msi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_4 * msh0_441[k]
                   - f_5 * msh1_441[k]
                   + f_3 * pc_z[k] * msi_590[k];

        t_762[k] = f_15 * lsi_594[k]
                   + f_8 * msh0_447[k]
                   - f_9 * msh1_447[k]
                   + f_3 * pc_x[k] * msi_594[k];

        t_763[k] = f_3 * pc_z[k] * msi_591[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, lsi_425, lsi_598, \
                         msh0_443, msh0_451, msh1_443, msh1_451, msi_593, msi_594, \
                         msi_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_22 * lsi_425[k]
                   + f_3 * pc_y[k] * msi_593[k];

        t_765[k] = f_6 * msh0_443[k]
                   - f_7 * msh1_443[k]
                   + f_3 * pc_z[k] * msi_593[k];

        t_766[k] = f_15 * lsi_598[k]
                   + f_6 * msh0_451[k]
                   - f_7 * msh1_451[k]
                   + f_3 * pc_x[k] * msi_598[k];

        t_767[k] = f_3 * pc_z[k] * msi_594[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_y, pc_z, lsi_429, msh0_444, msh0_446, \
                         msh1_444, msh1_446, msi_595, msi_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * msh0_444[k]
                   - f_5 * msh1_444[k]
                   + f_3 * pc_z[k] * msi_595[k];

        t_769[k] = f_22 * lsi_429[k]
                   + f_3 * pc_y[k] * msi_597[k];

        t_770[k] = f_8 * msh0_446[k]
                   - f_9 * msh1_446[k]
                   + f_3 * pc_z[k] * msi_597[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, pc_x, pc_z, lsi_603, msh0_447, msh0_456, \
                         msh1_447, msh1_456, msi_598, msi_599, \
                         msi_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_15 * lsi_603[k]
                   + f_4 * msh0_456[k]
                   - f_5 * msh1_456[k]
                   + f_3 * pc_x[k] * msi_603[k];

        t_772[k] = f_3 * pc_z[k] * msi_598[k];

        t_773[k] = f_4 * msh0_447[k]
                   - f_5 * msh1_447[k]
                   + f_3 * pc_z[k] * msi_599[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pc_x, pc_y, pc_z, lsi_434, lsi_609, \
                         msh0_448, msh0_450, msh1_448, msh1_450, msi_600, msi_602, \
                         msi_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_6 * msh0_448[k]
                   - f_7 * msh1_448[k]
                   + f_3 * pc_z[k] * msi_600[k];

        t_775[k] = f_22 * lsi_434[k]
                   + f_3 * pc_y[k] * msi_602[k];

        t_776[k] = f_10 * msh0_450[k]
                   - f_11 * msh1_450[k]
                   + f_3 * pc_z[k] * msi_602[k];

        t_777[k] = f_15 * lsi_609[k]
                   + f_3 * pc_x[k] * msi_609[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, pc_x, pc_z, lsi_611, lsi_612, \
                         lsi_613, lsi_614, msi_603, msi_611, msi_612, msi_613, \
                         msi_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_3 * pc_z[k] * msi_603[k];

        t_779[k] = f_15 * lsi_611[k]
                   + f_3 * pc_x[k] * msi_611[k];

        t_780[k] = f_15 * lsi_612[k]
                   + f_3 * pc_x[k] * msi_612[k];

        t_781[k] = f_15 * lsi_613[k]
                   + f_3 * pc_x[k] * msi_613[k];

        t_782[k] = f_15 * lsi_614[k]
                   + f_3 * pc_x[k] * msi_614[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pc_x, pc_y, pc_z, lsi_441, lsi_615, \
                         msh0_456, msh1_456, msi_609, msi_610, \
                         msi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_15 * lsi_615[k]
                   + f_3 * pc_x[k] * msi_615[k];

        t_784[k] = f_22 * lsi_441[k]
                   + f_1 * msh0_456[k]
                   - f_2 * msh1_456[k]
                   + f_3 * pc_y[k] * msi_609[k];

        t_785[k] = f_3 * pc_z[k] * msi_609[k];

        t_786[k] = f_4 * msh0_456[k]
                   - f_5 * msh1_456[k]
                   + f_3 * pc_z[k] * msi_610[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_z, msh0_457, msh0_458, msh0_459, msh1_457, \
                         msh1_458, msh1_459, msi_611, msi_612, \
                         msi_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_6 * msh0_457[k]
                   - f_7 * msh1_457[k]
                   + f_3 * pc_z[k] * msi_611[k];

        t_788[k] = f_8 * msh0_458[k]
                   - f_9 * msh1_458[k]
                   + f_3 * pc_z[k] * msi_612[k];

        t_789[k] = f_10 * msh0_459[k]
                   - f_11 * msh1_459[k]
                   + f_3 * pc_z[k] * msi_613[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_y, pc_z, lsk0_540, lsi_447, \
                         lsi_448, lsk1_540, msh0_461, msh1_461, msi_615, \
                         msi_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_22 * lsi_447[k]
                   + f_3 * pc_y[k] * msi_615[k];

        t_791[k] = f_1 * msh0_461[k]
                   - f_2 * msh1_461[k]
                   + f_3 * pc_z[k] * msi_615[k];

        t_792[k] = pa_z[k] * lsk0_540[k]
                   - f_12 * pc_z[k] * lsk1_540[k];

        t_793[k] = f_17 * lsi_448[k]
                   + f_3 * pc_y[k] * msi_616[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pa_z, pc_y, pc_z, lsk0_543, lsi_420, lsi_450, \
                         lsk1_543, msi_616, msi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * lsi_420[k]
                   + f_3 * pc_z[k] * msi_616[k];

        t_795[k] = pa_z[k] * lsk0_543[k]
                   - f_12 * pc_z[k] * lsk1_543[k];

        t_796[k] = f_17 * lsi_450[k]
                   + f_3 * pc_y[k] * msi_618[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pa_z, pc_x, pc_z, lsk0_546, lsi_423, lsi_621, \
                         lsk1_546, msh0_467, msh1_467, msi_619, \
                         msi_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_15 * lsi_621[k]
                   + f_10 * msh0_467[k]
                   - f_11 * msh1_467[k]
                   + f_3 * pc_x[k] * msi_621[k];

        t_798[k] = pa_z[k] * lsk0_546[k]
                   - f_12 * pc_z[k] * lsk1_546[k];

        t_799[k] = f_13 * lsi_423[k]
                   + f_3 * pc_z[k] * msi_619[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pa_z, pc_x, pc_y, pc_z, lsk0_550, lsi_453, \
                         lsi_625, lsk1_550, msh0_471, msh1_471, msi_621, \
                         msi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * lsi_453[k]
                   + f_3 * pc_y[k] * msi_621[k];

        t_801[k] = f_15 * lsi_625[k]
                   + f_8 * msh0_471[k]
                   - f_9 * msh1_471[k]
                   + f_3 * pc_x[k] * msi_625[k];

        t_802[k] = pa_z[k] * lsk0_550[k]
                   - f_12 * pc_z[k] * lsk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pa_z, pc_y, pc_z, lsk0_552, lsi_426, lsi_427, \
                         lsi_457, lsk1_552, msi_622, msi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_13 * lsi_426[k]
                   + f_3 * pc_z[k] * msi_622[k];

        t_804[k] = pa_z[k] * lsk0_552[k]
                   + f_14 * lsi_427[k]
                   - f_12 * pc_z[k] * lsk1_552[k];

        t_805[k] = f_17 * lsi_457[k]
                   + f_3 * pc_y[k] * msi_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pa_z, pc_x, pc_z, lsk0_555, lsi_430, lsi_630, \
                         lsk1_555, msh0_476, msh1_476, msi_626, \
                         msi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_15 * lsi_630[k]
                   + f_6 * msh0_476[k]
                   - f_7 * msh1_476[k]
                   + f_3 * pc_x[k] * msi_630[k];

        t_807[k] = pa_z[k] * lsk0_555[k]
                   - f_12 * pc_z[k] * lsk1_555[k];

        t_808[k] = f_13 * lsi_430[k]
                   + f_3 * pc_z[k] * msi_626[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pa_z, pc_y, pc_z, lsk0_557, lsk0_558, lsi_431, \
                         lsi_432, lsi_462, lsk1_557, lsk1_558, \
                         msi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pa_z[k] * lsk0_557[k]
                   + f_14 * lsi_431[k]
                   - f_12 * pc_z[k] * lsk1_557[k];

        t_810[k] = pa_z[k] * lsk0_558[k]
                   + f_15 * lsi_432[k]
                   - f_12 * pc_z[k] * lsk1_558[k];

        t_811[k] = f_17 * lsi_462[k]
                   + f_3 * pc_y[k] * msi_630[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, lsi_636, lsi_637, lsi_638, lsi_639, \
                         msh0_482, msh1_482, msi_636, msi_637, msi_638, \
                         msi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_15 * lsi_636[k]
                   + f_4 * msh0_482[k]
                   - f_5 * msh1_482[k]
                   + f_3 * pc_x[k] * msi_636[k];

        t_813[k] = f_15 * lsi_637[k]
                   + f_3 * pc_x[k] * msi_637[k];

        t_814[k] = f_15 * lsi_638[k]
                   + f_3 * pc_x[k] * msi_638[k];

        t_815[k] = f_15 * lsi_639[k]
                   + f_3 * pc_x[k] * msi_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, pc_x, lsi_640, lsi_641, lsi_642, lsi_643, \
                         msi_640, msi_641, msi_642, msi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_15 * lsi_640[k]
                   + f_3 * pc_x[k] * msi_640[k];

        t_817[k] = f_15 * lsi_641[k]
                   + f_3 * pc_x[k] * msi_641[k];

        t_818[k] = f_15 * lsi_642[k]
                   + f_3 * pc_x[k] * msi_642[k];

        t_819[k] = f_15 * lsi_643[k]
                   + f_3 * pc_x[k] * msi_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_y, pc_z, lsk0_568, lsi_441, lsi_471, \
                         lsk1_568, msh0_479, msh1_479, msi_637, \
                         msi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pa_z[k] * lsk0_568[k]
                   - f_12 * pc_z[k] * lsk1_568[k];

        t_821[k] = f_13 * lsi_441[k]
                   + f_3 * pc_z[k] * msi_637[k];

        t_822[k] = f_17 * lsi_471[k]
                   + f_10 * msh0_479[k]
                   - f_11 * msh1_479[k]
                   + f_3 * pc_y[k] * msi_639[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_y, lsi_472, lsi_473, lsi_474, msh0_480, \
                         msh0_481, msh0_482, msh1_480, msh1_481, msh1_482, msi_640, msi_641, \
                         msi_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_17 * lsi_472[k]
                   + f_8 * msh0_480[k]
                   - f_9 * msh1_480[k]
                   + f_3 * pc_y[k] * msi_640[k];

        t_824[k] = f_17 * lsi_473[k]
                   + f_6 * msh0_481[k]
                   - f_7 * msh1_481[k]
                   + f_3 * pc_y[k] * msi_641[k];

        t_825[k] = f_17 * lsi_474[k]
                   + f_4 * msh0_482[k]
                   - f_5 * msh1_482[k]
                   + f_3 * pc_y[k] * msi_642[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t lsi, const size_t msh0,
                                                          const size_t msh1, const size_t msi,
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

    const auto *lsi_447 = buffer.data(lsi + 447);
    const auto *lsi_448 = buffer.data(lsi + 448);
    const auto *lsi_451 = buffer.data(lsi + 451);
    const auto *lsi_454 = buffer.data(lsi + 454);
    const auto *lsi_458 = buffer.data(lsi + 458);
    const auto *lsi_469 = buffer.data(lsi + 469);
    const auto *lsi_475 = buffer.data(lsi + 475);
    const auto *lsi_476 = buffer.data(lsi + 476);
    const auto *lsi_478 = buffer.data(lsi + 478);
    const auto *lsi_479 = buffer.data(lsi + 479);
    const auto *lsi_481 = buffer.data(lsi + 481);
    const auto *lsi_482 = buffer.data(lsi + 482);
    const auto *lsi_485 = buffer.data(lsi + 485);
    const auto *lsi_486 = buffer.data(lsi + 486);
    const auto *lsi_490 = buffer.data(lsi + 490);
    const auto *lsi_497 = buffer.data(lsi + 497);
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
    const auto *lsi_518 = buffer.data(lsi + 518);
    const auto *lsi_525 = buffer.data(lsi + 525);
    const auto *lsi_527 = buffer.data(lsi + 527);
    const auto *lsi_528 = buffer.data(lsi + 528);
    const auto *lsi_529 = buffer.data(lsi + 529);
    const auto *lsi_530 = buffer.data(lsi + 530);
    const auto *lsi_531 = buffer.data(lsi + 531);
    const auto *lsi_532 = buffer.data(lsi + 532);
    const auto *lsi_534 = buffer.data(lsi + 534);
    const auto *lsi_537 = buffer.data(lsi + 537);
    const auto *lsi_541 = buffer.data(lsi + 541);
    const auto *lsi_546 = buffer.data(lsi + 546);
    const auto *lsi_553 = buffer.data(lsi + 553);
    const auto *lsi_555 = buffer.data(lsi + 555);
    const auto *lsi_644 = buffer.data(lsi + 644);
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

    const auto *msh0_482 = buffer.data(msh0 + 482);
    const auto *msh0_483 = buffer.data(msh0 + 483);
    const auto *msh0_486 = buffer.data(msh0 + 486);
    const auto *msh0_488 = buffer.data(msh0 + 488);
    const auto *msh0_489 = buffer.data(msh0 + 489);
    const auto *msh0_492 = buffer.data(msh0 + 492);
    const auto *msh0_493 = buffer.data(msh0 + 493);
    const auto *msh0_495 = buffer.data(msh0 + 495);
    const auto *msh0_497 = buffer.data(msh0 + 497);
    const auto *msh0_498 = buffer.data(msh0 + 498);
    const auto *msh0_500 = buffer.data(msh0 + 500);
    const auto *msh0_501 = buffer.data(msh0 + 501);
    const auto *msh0_502 = buffer.data(msh0 + 502);
    const auto *msh0_503 = buffer.data(msh0 + 503);
    const auto *msh0_504 = buffer.data(msh0 + 504);
    const auto *msh0_507 = buffer.data(msh0 + 507);
    const auto *msh0_509 = buffer.data(msh0 + 509);
    const auto *msh0_510 = buffer.data(msh0 + 510);
    const auto *msh0_513 = buffer.data(msh0 + 513);
    const auto *msh0_514 = buffer.data(msh0 + 514);
    const auto *msh0_516 = buffer.data(msh0 + 516);
    const auto *msh0_518 = buffer.data(msh0 + 518);
    const auto *msh0_519 = buffer.data(msh0 + 519);
    const auto *msh0_521 = buffer.data(msh0 + 521);
    const auto *msh0_522 = buffer.data(msh0 + 522);
    const auto *msh0_523 = buffer.data(msh0 + 523);
    const auto *msh0_524 = buffer.data(msh0 + 524);
    const auto *msh0_525 = buffer.data(msh0 + 525);
    const auto *msh0_528 = buffer.data(msh0 + 528);
    const auto *msh0_530 = buffer.data(msh0 + 530);
    const auto *msh0_531 = buffer.data(msh0 + 531);
    const auto *msh0_534 = buffer.data(msh0 + 534);
    const auto *msh0_535 = buffer.data(msh0 + 535);
    const auto *msh0_537 = buffer.data(msh0 + 537);
    const auto *msh0_539 = buffer.data(msh0 + 539);
    const auto *msh0_540 = buffer.data(msh0 + 540);
    const auto *msh0_542 = buffer.data(msh0 + 542);
    const auto *msh0_543 = buffer.data(msh0 + 543);
    const auto *msh0_545 = buffer.data(msh0 + 545);

    const auto *msh1_482 = buffer.data(msh1 + 482);
    const auto *msh1_483 = buffer.data(msh1 + 483);
    const auto *msh1_486 = buffer.data(msh1 + 486);
    const auto *msh1_488 = buffer.data(msh1 + 488);
    const auto *msh1_489 = buffer.data(msh1 + 489);
    const auto *msh1_492 = buffer.data(msh1 + 492);
    const auto *msh1_493 = buffer.data(msh1 + 493);
    const auto *msh1_495 = buffer.data(msh1 + 495);
    const auto *msh1_497 = buffer.data(msh1 + 497);
    const auto *msh1_498 = buffer.data(msh1 + 498);
    const auto *msh1_500 = buffer.data(msh1 + 500);
    const auto *msh1_501 = buffer.data(msh1 + 501);
    const auto *msh1_502 = buffer.data(msh1 + 502);
    const auto *msh1_503 = buffer.data(msh1 + 503);
    const auto *msh1_504 = buffer.data(msh1 + 504);
    const auto *msh1_507 = buffer.data(msh1 + 507);
    const auto *msh1_509 = buffer.data(msh1 + 509);
    const auto *msh1_510 = buffer.data(msh1 + 510);
    const auto *msh1_513 = buffer.data(msh1 + 513);
    const auto *msh1_514 = buffer.data(msh1 + 514);
    const auto *msh1_516 = buffer.data(msh1 + 516);
    const auto *msh1_518 = buffer.data(msh1 + 518);
    const auto *msh1_519 = buffer.data(msh1 + 519);
    const auto *msh1_521 = buffer.data(msh1 + 521);
    const auto *msh1_522 = buffer.data(msh1 + 522);
    const auto *msh1_523 = buffer.data(msh1 + 523);
    const auto *msh1_524 = buffer.data(msh1 + 524);
    const auto *msh1_525 = buffer.data(msh1 + 525);
    const auto *msh1_528 = buffer.data(msh1 + 528);
    const auto *msh1_530 = buffer.data(msh1 + 530);
    const auto *msh1_531 = buffer.data(msh1 + 531);
    const auto *msh1_534 = buffer.data(msh1 + 534);
    const auto *msh1_535 = buffer.data(msh1 + 535);
    const auto *msh1_537 = buffer.data(msh1 + 537);
    const auto *msh1_539 = buffer.data(msh1 + 539);
    const auto *msh1_540 = buffer.data(msh1 + 540);
    const auto *msh1_542 = buffer.data(msh1 + 542);
    const auto *msh1_543 = buffer.data(msh1 + 543);
    const auto *msh1_545 = buffer.data(msh1 + 545);

    const auto *msi_643 = buffer.data(msi + 643);
    const auto *msi_644 = buffer.data(msi + 644);
    const auto *msi_646 = buffer.data(msi + 646);
    const auto *msi_647 = buffer.data(msi + 647);
    const auto *msi_649 = buffer.data(msi + 649);
    const auto *msi_650 = buffer.data(msi + 650);
    const auto *msi_653 = buffer.data(msi + 653);
    const auto *msi_654 = buffer.data(msi + 654);
    const auto *msi_656 = buffer.data(msi + 656);
    const auto *msi_658 = buffer.data(msi + 658);
    const auto *msi_659 = buffer.data(msi + 659);
    const auto *msi_661 = buffer.data(msi + 661);
    const auto *msi_662 = buffer.data(msi + 662);
    const auto *msi_664 = buffer.data(msi + 664);
    const auto *msi_665 = buffer.data(msi + 665);
    const auto *msi_666 = buffer.data(msi + 666);
    const auto *msi_667 = buffer.data(msi + 667);
    const auto *msi_668 = buffer.data(msi + 668);
    const auto *msi_669 = buffer.data(msi + 669);
    const auto *msi_670 = buffer.data(msi + 670);
    const auto *msi_671 = buffer.data(msi + 671);
    const auto *msi_672 = buffer.data(msi + 672);
    const auto *msi_674 = buffer.data(msi + 674);
    const auto *msi_675 = buffer.data(msi + 675);
    const auto *msi_677 = buffer.data(msi + 677);
    const auto *msi_678 = buffer.data(msi + 678);
    const auto *msi_681 = buffer.data(msi + 681);
    const auto *msi_682 = buffer.data(msi + 682);
    const auto *msi_684 = buffer.data(msi + 684);
    const auto *msi_686 = buffer.data(msi + 686);
    const auto *msi_687 = buffer.data(msi + 687);
    const auto *msi_689 = buffer.data(msi + 689);
    const auto *msi_690 = buffer.data(msi + 690);
    const auto *msi_692 = buffer.data(msi + 692);
    const auto *msi_693 = buffer.data(msi + 693);
    const auto *msi_694 = buffer.data(msi + 694);
    const auto *msi_695 = buffer.data(msi + 695);
    const auto *msi_696 = buffer.data(msi + 696);
    const auto *msi_697 = buffer.data(msi + 697);
    const auto *msi_698 = buffer.data(msi + 698);
    const auto *msi_699 = buffer.data(msi + 699);
    const auto *msi_700 = buffer.data(msi + 700);
    const auto *msi_702 = buffer.data(msi + 702);
    const auto *msi_703 = buffer.data(msi + 703);
    const auto *msi_705 = buffer.data(msi + 705);
    const auto *msi_706 = buffer.data(msi + 706);
    const auto *msi_709 = buffer.data(msi + 709);
    const auto *msi_710 = buffer.data(msi + 710);
    const auto *msi_712 = buffer.data(msi + 712);
    const auto *msi_714 = buffer.data(msi + 714);
    const auto *msi_715 = buffer.data(msi + 715);
    const auto *msi_717 = buffer.data(msi + 717);
    const auto *msi_718 = buffer.data(msi + 718);
    const auto *msi_720 = buffer.data(msi + 720);
    const auto *msi_721 = buffer.data(msi + 721);
    const auto *msi_722 = buffer.data(msi + 722);
    const auto *msi_723 = buffer.data(msi + 723);
    const auto *msi_724 = buffer.data(msi + 724);
    const auto *msi_725 = buffer.data(msi + 725);
    const auto *msi_726 = buffer.data(msi + 726);
    const auto *msi_727 = buffer.data(msi + 727);

#pragma omp simd aligned(t_826, t_827, t_828, pc_x, pc_y, pc_z, lsi_447, lsi_475, lsi_644, \
                         msh0_482, msh0_483, msh1_482, msh1_483, msi_643, \
                         msi_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * lsi_475[k]
                   + f_3 * pc_y[k] * msi_643[k];

        t_827[k] = f_13 * lsi_447[k]
                   + f_1 * msh0_482[k]
                   - f_2 * msh1_482[k]
                   + f_3 * pc_z[k] * msi_643[k];

        t_828[k] = f_15 * lsi_644[k]
                   + f_1 * msh0_483[k]
                   - f_2 * msh1_483[k]
                   + f_3 * pc_x[k] * msi_644[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, pc_y, pc_z, lsi_448, lsi_476, \
                         lsi_478, lsi_647, msh0_486, msh1_486, msi_644, msi_646, \
                         msi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_16 * lsi_476[k]
                   + f_3 * pc_y[k] * msi_644[k];

        t_830[k] = f_14 * lsi_448[k]
                   + f_3 * pc_z[k] * msi_644[k];

        t_831[k] = f_15 * lsi_647[k]
                   + f_10 * msh0_486[k]
                   - f_11 * msh1_486[k]
                   + f_3 * pc_x[k] * msi_647[k];

        t_832[k] = f_16 * lsi_478[k]
                   + f_3 * pc_y[k] * msi_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, pc_z, lsi_451, lsi_649, lsi_650, msh0_488, \
                         msh0_489, msh1_488, msh1_489, msi_647, msi_649, \
                         msi_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_15 * lsi_649[k]
                   + f_10 * msh0_488[k]
                   - f_11 * msh1_488[k]
                   + f_3 * pc_x[k] * msi_649[k];

        t_834[k] = f_15 * lsi_650[k]
                   + f_8 * msh0_489[k]
                   - f_9 * msh1_489[k]
                   + f_3 * pc_x[k] * msi_650[k];

        t_835[k] = f_14 * lsi_451[k]
                   + f_3 * pc_z[k] * msi_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_x, pc_y, lsi_481, lsi_653, lsi_654, msh0_492, \
                         msh0_493, msh1_492, msh1_493, msi_649, msi_653, \
                         msi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * lsi_481[k]
                   + f_3 * pc_y[k] * msi_649[k];

        t_837[k] = f_15 * lsi_653[k]
                   + f_8 * msh0_492[k]
                   - f_9 * msh1_492[k]
                   + f_3 * pc_x[k] * msi_653[k];

        t_838[k] = f_15 * lsi_654[k]
                   + f_6 * msh0_493[k]
                   - f_7 * msh1_493[k]
                   + f_3 * pc_x[k] * msi_654[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_y, pc_z, lsi_454, lsi_485, lsi_656, \
                         msh0_495, msh1_495, msi_650, msi_653, \
                         msi_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * lsi_454[k]
                   + f_3 * pc_z[k] * msi_650[k];

        t_840[k] = f_15 * lsi_656[k]
                   + f_6 * msh0_495[k]
                   - f_7 * msh1_495[k]
                   + f_3 * pc_x[k] * msi_656[k];

        t_841[k] = f_16 * lsi_485[k]
                   + f_3 * pc_y[k] * msi_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_z, lsi_458, lsi_658, lsi_659, msh0_497, \
                         msh0_498, msh1_497, msh1_498, msi_654, msi_658, \
                         msi_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_15 * lsi_658[k]
                   + f_6 * msh0_497[k]
                   - f_7 * msh1_497[k]
                   + f_3 * pc_x[k] * msi_658[k];

        t_843[k] = f_15 * lsi_659[k]
                   + f_4 * msh0_498[k]
                   - f_5 * msh1_498[k]
                   + f_3 * pc_x[k] * msi_659[k];

        t_844[k] = f_14 * lsi_458[k]
                   + f_3 * pc_z[k] * msi_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, lsi_490, lsi_661, lsi_662, msh0_500, \
                         msh0_501, msh1_500, msh1_501, msi_658, msi_661, \
                         msi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_15 * lsi_661[k]
                   + f_4 * msh0_500[k]
                   - f_5 * msh1_500[k]
                   + f_3 * pc_x[k] * msi_661[k];

        t_846[k] = f_15 * lsi_662[k]
                   + f_4 * msh0_501[k]
                   - f_5 * msh1_501[k]
                   + f_3 * pc_x[k] * msi_662[k];

        t_847[k] = f_16 * lsi_490[k]
                   + f_3 * pc_y[k] * msi_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, lsi_664, lsi_665, lsi_666, lsi_667, \
                         msh0_503, msh1_503, msi_664, msi_665, msi_666, \
                         msi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_15 * lsi_664[k]
                   + f_4 * msh0_503[k]
                   - f_5 * msh1_503[k]
                   + f_3 * pc_x[k] * msi_664[k];

        t_849[k] = f_15 * lsi_665[k]
                   + f_3 * pc_x[k] * msi_665[k];

        t_850[k] = f_15 * lsi_666[k]
                   + f_3 * pc_x[k] * msi_666[k];

        t_851[k] = f_15 * lsi_667[k]
                   + f_3 * pc_x[k] * msi_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, lsi_668, lsi_669, lsi_670, lsi_671, \
                         msi_668, msi_669, msi_670, msi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_15 * lsi_668[k]
                   + f_3 * pc_x[k] * msi_668[k];

        t_853[k] = f_15 * lsi_669[k]
                   + f_3 * pc_x[k] * msi_669[k];

        t_854[k] = f_15 * lsi_670[k]
                   + f_3 * pc_x[k] * msi_670[k];

        t_855[k] = f_15 * lsi_671[k]
                   + f_3 * pc_x[k] * msi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, lsi_469, lsi_497, lsi_499, msh0_498, \
                         msh0_500, msh1_498, msh1_500, msi_665, \
                         msi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * lsi_497[k]
                   + f_1 * msh0_498[k]
                   - f_2 * msh1_498[k]
                   + f_3 * pc_y[k] * msi_665[k];

        t_857[k] = f_14 * lsi_469[k]
                   + f_3 * pc_z[k] * msi_665[k];

        t_858[k] = f_16 * lsi_499[k]
                   + f_10 * msh0_500[k]
                   - f_11 * msh1_500[k]
                   + f_3 * pc_y[k] * msi_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, lsi_500, lsi_501, lsi_502, msh0_501, \
                         msh0_502, msh0_503, msh1_501, msh1_502, msh1_503, msi_668, msi_669, \
                         msi_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * lsi_500[k]
                   + f_8 * msh0_501[k]
                   - f_9 * msh1_501[k]
                   + f_3 * pc_y[k] * msi_668[k];

        t_860[k] = f_16 * lsi_501[k]
                   + f_6 * msh0_502[k]
                   - f_7 * msh1_502[k]
                   + f_3 * pc_y[k] * msi_669[k];

        t_861[k] = f_16 * lsi_502[k]
                   + f_4 * msh0_503[k]
                   - f_5 * msh1_503[k]
                   + f_3 * pc_y[k] * msi_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, lsi_475, lsi_503, lsi_672, \
                         msh0_503, msh0_504, msh1_503, msh1_504, msi_671, \
                         msi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * lsi_503[k]
                   + f_3 * pc_y[k] * msi_671[k];

        t_863[k] = f_14 * lsi_475[k]
                   + f_1 * msh0_503[k]
                   - f_2 * msh1_503[k]
                   + f_3 * pc_z[k] * msi_671[k];

        t_864[k] = f_15 * lsi_672[k]
                   + f_1 * msh0_504[k]
                   - f_2 * msh1_504[k]
                   + f_3 * pc_x[k] * msi_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, lsi_476, lsi_504, \
                         lsi_506, lsi_675, msh0_507, msh1_507, msi_672, msi_674, \
                         msi_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_15 * lsi_504[k]
                   + f_3 * pc_y[k] * msi_672[k];

        t_866[k] = f_15 * lsi_476[k]
                   + f_3 * pc_z[k] * msi_672[k];

        t_867[k] = f_15 * lsi_675[k]
                   + f_10 * msh0_507[k]
                   - f_11 * msh1_507[k]
                   + f_3 * pc_x[k] * msi_675[k];

        t_868[k] = f_15 * lsi_506[k]
                   + f_3 * pc_y[k] * msi_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, lsi_479, lsi_677, lsi_678, msh0_509, \
                         msh0_510, msh1_509, msh1_510, msi_675, msi_677, \
                         msi_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_15 * lsi_677[k]
                   + f_10 * msh0_509[k]
                   - f_11 * msh1_509[k]
                   + f_3 * pc_x[k] * msi_677[k];

        t_870[k] = f_15 * lsi_678[k]
                   + f_8 * msh0_510[k]
                   - f_9 * msh1_510[k]
                   + f_3 * pc_x[k] * msi_678[k];

        t_871[k] = f_15 * lsi_479[k]
                   + f_3 * pc_z[k] * msi_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, pc_y, lsi_509, lsi_681, lsi_682, msh0_513, \
                         msh0_514, msh1_513, msh1_514, msi_677, msi_681, \
                         msi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * lsi_509[k]
                   + f_3 * pc_y[k] * msi_677[k];

        t_873[k] = f_15 * lsi_681[k]
                   + f_8 * msh0_513[k]
                   - f_9 * msh1_513[k]
                   + f_3 * pc_x[k] * msi_681[k];

        t_874[k] = f_15 * lsi_682[k]
                   + f_6 * msh0_514[k]
                   - f_7 * msh1_514[k]
                   + f_3 * pc_x[k] * msi_682[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, pc_y, pc_z, lsi_482, lsi_513, lsi_684, \
                         msh0_516, msh1_516, msi_678, msi_681, \
                         msi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * lsi_482[k]
                   + f_3 * pc_z[k] * msi_678[k];

        t_876[k] = f_15 * lsi_684[k]
                   + f_6 * msh0_516[k]
                   - f_7 * msh1_516[k]
                   + f_3 * pc_x[k] * msi_684[k];

        t_877[k] = f_15 * lsi_513[k]
                   + f_3 * pc_y[k] * msi_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, pc_z, lsi_486, lsi_686, lsi_687, msh0_518, \
                         msh0_519, msh1_518, msh1_519, msi_682, msi_686, \
                         msi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_15 * lsi_686[k]
                   + f_6 * msh0_518[k]
                   - f_7 * msh1_518[k]
                   + f_3 * pc_x[k] * msi_686[k];

        t_879[k] = f_15 * lsi_687[k]
                   + f_4 * msh0_519[k]
                   - f_5 * msh1_519[k]
                   + f_3 * pc_x[k] * msi_687[k];

        t_880[k] = f_15 * lsi_486[k]
                   + f_3 * pc_z[k] * msi_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_y, lsi_518, lsi_689, lsi_690, msh0_521, \
                         msh0_522, msh1_521, msh1_522, msi_686, msi_689, \
                         msi_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_15 * lsi_689[k]
                   + f_4 * msh0_521[k]
                   - f_5 * msh1_521[k]
                   + f_3 * pc_x[k] * msi_689[k];

        t_882[k] = f_15 * lsi_690[k]
                   + f_4 * msh0_522[k]
                   - f_5 * msh1_522[k]
                   + f_3 * pc_x[k] * msi_690[k];

        t_883[k] = f_15 * lsi_518[k]
                   + f_3 * pc_y[k] * msi_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pc_x, lsi_692, lsi_693, lsi_694, lsi_695, \
                         msh0_524, msh1_524, msi_692, msi_693, msi_694, \
                         msi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_15 * lsi_692[k]
                   + f_4 * msh0_524[k]
                   - f_5 * msh1_524[k]
                   + f_3 * pc_x[k] * msi_692[k];

        t_885[k] = f_15 * lsi_693[k]
                   + f_3 * pc_x[k] * msi_693[k];

        t_886[k] = f_15 * lsi_694[k]
                   + f_3 * pc_x[k] * msi_694[k];

        t_887[k] = f_15 * lsi_695[k]
                   + f_3 * pc_x[k] * msi_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, lsi_696, lsi_697, lsi_698, lsi_699, \
                         msi_696, msi_697, msi_698, msi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_15 * lsi_696[k]
                   + f_3 * pc_x[k] * msi_696[k];

        t_889[k] = f_15 * lsi_697[k]
                   + f_3 * pc_x[k] * msi_697[k];

        t_890[k] = f_15 * lsi_698[k]
                   + f_3 * pc_x[k] * msi_698[k];

        t_891[k] = f_15 * lsi_699[k]
                   + f_3 * pc_x[k] * msi_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, lsi_497, lsi_525, lsi_527, msh0_519, \
                         msh0_521, msh1_519, msh1_521, msi_693, \
                         msi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * lsi_525[k]
                   + f_1 * msh0_519[k]
                   - f_2 * msh1_519[k]
                   + f_3 * pc_y[k] * msi_693[k];

        t_893[k] = f_15 * lsi_497[k]
                   + f_3 * pc_z[k] * msi_693[k];

        t_894[k] = f_15 * lsi_527[k]
                   + f_10 * msh0_521[k]
                   - f_11 * msh1_521[k]
                   + f_3 * pc_y[k] * msi_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, lsi_528, lsi_529, lsi_530, msh0_522, \
                         msh0_523, msh0_524, msh1_522, msh1_523, msh1_524, msi_696, msi_697, \
                         msi_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * lsi_528[k]
                   + f_8 * msh0_522[k]
                   - f_9 * msh1_522[k]
                   + f_3 * pc_y[k] * msi_696[k];

        t_896[k] = f_15 * lsi_529[k]
                   + f_6 * msh0_523[k]
                   - f_7 * msh1_523[k]
                   + f_3 * pc_y[k] * msi_697[k];

        t_897[k] = f_15 * lsi_530[k]
                   + f_4 * msh0_524[k]
                   - f_5 * msh1_524[k]
                   + f_3 * pc_y[k] * msi_698[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_x, pc_y, pc_z, lsi_503, lsi_531, lsi_700, \
                         msh0_524, msh0_525, msh1_524, msh1_525, msi_699, \
                         msi_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * lsi_531[k]
                   + f_3 * pc_y[k] * msi_699[k];

        t_899[k] = f_15 * lsi_503[k]
                   + f_1 * msh0_524[k]
                   - f_2 * msh1_524[k]
                   + f_3 * pc_z[k] * msi_699[k];

        t_900[k] = f_15 * lsi_700[k]
                   + f_1 * msh0_525[k]
                   - f_2 * msh1_525[k]
                   + f_3 * pc_x[k] * msi_700[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pc_x, pc_y, pc_z, lsi_504, lsi_532, \
                         lsi_534, lsi_703, msh0_528, msh1_528, msi_700, msi_702, \
                         msi_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_14 * lsi_532[k]
                   + f_3 * pc_y[k] * msi_700[k];

        t_902[k] = f_16 * lsi_504[k]
                   + f_3 * pc_z[k] * msi_700[k];

        t_903[k] = f_15 * lsi_703[k]
                   + f_10 * msh0_528[k]
                   - f_11 * msh1_528[k]
                   + f_3 * pc_x[k] * msi_703[k];

        t_904[k] = f_14 * lsi_534[k]
                   + f_3 * pc_y[k] * msi_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, lsi_507, lsi_705, lsi_706, msh0_530, \
                         msh0_531, msh1_530, msh1_531, msi_703, msi_705, \
                         msi_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_15 * lsi_705[k]
                   + f_10 * msh0_530[k]
                   - f_11 * msh1_530[k]
                   + f_3 * pc_x[k] * msi_705[k];

        t_906[k] = f_15 * lsi_706[k]
                   + f_8 * msh0_531[k]
                   - f_9 * msh1_531[k]
                   + f_3 * pc_x[k] * msi_706[k];

        t_907[k] = f_16 * lsi_507[k]
                   + f_3 * pc_z[k] * msi_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, lsi_537, lsi_709, lsi_710, msh0_534, \
                         msh0_535, msh1_534, msh1_535, msi_705, msi_709, \
                         msi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * lsi_537[k]
                   + f_3 * pc_y[k] * msi_705[k];

        t_909[k] = f_15 * lsi_709[k]
                   + f_8 * msh0_534[k]
                   - f_9 * msh1_534[k]
                   + f_3 * pc_x[k] * msi_709[k];

        t_910[k] = f_15 * lsi_710[k]
                   + f_6 * msh0_535[k]
                   - f_7 * msh1_535[k]
                   + f_3 * pc_x[k] * msi_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, pc_y, pc_z, lsi_510, lsi_541, lsi_712, \
                         msh0_537, msh1_537, msi_706, msi_709, \
                         msi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * lsi_510[k]
                   + f_3 * pc_z[k] * msi_706[k];

        t_912[k] = f_15 * lsi_712[k]
                   + f_6 * msh0_537[k]
                   - f_7 * msh1_537[k]
                   + f_3 * pc_x[k] * msi_712[k];

        t_913[k] = f_14 * lsi_541[k]
                   + f_3 * pc_y[k] * msi_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, lsi_514, lsi_714, lsi_715, msh0_539, \
                         msh0_540, msh1_539, msh1_540, msi_710, msi_714, \
                         msi_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_15 * lsi_714[k]
                   + f_6 * msh0_539[k]
                   - f_7 * msh1_539[k]
                   + f_3 * pc_x[k] * msi_714[k];

        t_915[k] = f_15 * lsi_715[k]
                   + f_4 * msh0_540[k]
                   - f_5 * msh1_540[k]
                   + f_3 * pc_x[k] * msi_715[k];

        t_916[k] = f_16 * lsi_514[k]
                   + f_3 * pc_z[k] * msi_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, pc_y, lsi_546, lsi_717, lsi_718, msh0_542, \
                         msh0_543, msh1_542, msh1_543, msi_714, msi_717, \
                         msi_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_15 * lsi_717[k]
                   + f_4 * msh0_542[k]
                   - f_5 * msh1_542[k]
                   + f_3 * pc_x[k] * msi_717[k];

        t_918[k] = f_15 * lsi_718[k]
                   + f_4 * msh0_543[k]
                   - f_5 * msh1_543[k]
                   + f_3 * pc_x[k] * msi_718[k];

        t_919[k] = f_14 * lsi_546[k]
                   + f_3 * pc_y[k] * msi_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, lsi_720, lsi_721, lsi_722, lsi_723, \
                         msh0_545, msh1_545, msi_720, msi_721, msi_722, \
                         msi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_15 * lsi_720[k]
                   + f_4 * msh0_545[k]
                   - f_5 * msh1_545[k]
                   + f_3 * pc_x[k] * msi_720[k];

        t_921[k] = f_15 * lsi_721[k]
                   + f_3 * pc_x[k] * msi_721[k];

        t_922[k] = f_15 * lsi_722[k]
                   + f_3 * pc_x[k] * msi_722[k];

        t_923[k] = f_15 * lsi_723[k]
                   + f_3 * pc_x[k] * msi_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, lsi_724, lsi_725, lsi_726, lsi_727, \
                         msi_724, msi_725, msi_726, msi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_15 * lsi_724[k]
                   + f_3 * pc_x[k] * msi_724[k];

        t_925[k] = f_15 * lsi_725[k]
                   + f_3 * pc_x[k] * msi_725[k];

        t_926[k] = f_15 * lsi_726[k]
                   + f_3 * pc_x[k] * msi_726[k];

        t_927[k] = f_15 * lsi_727[k]
                   + f_3 * pc_x[k] * msi_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, lsi_525, lsi_553, lsi_555, msh0_540, \
                         msh0_542, msh1_540, msh1_542, msi_721, \
                         msi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * lsi_553[k]
                   + f_1 * msh0_540[k]
                   - f_2 * msh1_540[k]
                   + f_3 * pc_y[k] * msi_721[k];

        t_929[k] = f_16 * lsi_525[k]
                   + f_3 * pc_z[k] * msi_721[k];

        t_930[k] = f_14 * lsi_555[k]
                   + f_10 * msh0_542[k]
                   - f_11 * msh1_542[k]
                   + f_3 * pc_y[k] * msi_723[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_720 = buffer.data(lsk0 + 720);
    const auto *lsk0_723 = buffer.data(lsk0 + 723);
    const auto *lsk0_725 = buffer.data(lsk0 + 725);
    const auto *lsk0_726 = buffer.data(lsk0 + 726);
    const auto *lsk0_729 = buffer.data(lsk0 + 729);
    const auto *lsk0_730 = buffer.data(lsk0 + 730);
    const auto *lsk0_732 = buffer.data(lsk0 + 732);
    const auto *lsk0_734 = buffer.data(lsk0 + 734);
    const auto *lsk0_735 = buffer.data(lsk0 + 735);
    const auto *lsk0_737 = buffer.data(lsk0 + 737);
    const auto *lsk0_738 = buffer.data(lsk0 + 738);
    const auto *lsk0_740 = buffer.data(lsk0 + 740);
    const auto *lsk0_755 = buffer.data(lsk0 + 755);
    const auto *lsk0_756 = buffer.data(lsk0 + 756);
    const auto *lsk0_759 = buffer.data(lsk0 + 759);

    const auto *lsi_531 = buffer.data(lsi + 531);
    const auto *lsi_532 = buffer.data(lsi + 532);
    const auto *lsi_535 = buffer.data(lsi + 535);
    const auto *lsi_538 = buffer.data(lsi + 538);
    const auto *lsi_542 = buffer.data(lsi + 542);
    const auto *lsi_553 = buffer.data(lsi + 553);
    const auto *lsi_556 = buffer.data(lsi + 556);
    const auto *lsi_557 = buffer.data(lsi + 557);
    const auto *lsi_558 = buffer.data(lsi + 558);
    const auto *lsi_559 = buffer.data(lsi + 559);
    const auto *lsi_560 = buffer.data(lsi + 560);
    const auto *lsi_561 = buffer.data(lsi + 561);
    const auto *lsi_562 = buffer.data(lsi + 562);
    const auto *lsi_563 = buffer.data(lsi + 563);
    const auto *lsi_565 = buffer.data(lsi + 565);
    const auto *lsi_566 = buffer.data(lsi + 566);
    const auto *lsi_568 = buffer.data(lsi + 568);
    const auto *lsi_569 = buffer.data(lsi + 569);
    const auto *lsi_570 = buffer.data(lsi + 570);
    const auto *lsi_572 = buffer.data(lsi + 572);
    const auto *lsi_573 = buffer.data(lsi + 573);
    const auto *lsi_574 = buffer.data(lsi + 574);
    const auto *lsi_581 = buffer.data(lsi + 581);
    const auto *lsi_583 = buffer.data(lsi + 583);
    const auto *lsi_584 = buffer.data(lsi + 584);
    const auto *lsi_585 = buffer.data(lsi + 585);
    const auto *lsi_586 = buffer.data(lsi + 586);
    const auto *lsi_587 = buffer.data(lsi + 587);
    const auto *lsi_588 = buffer.data(lsi + 588);
    const auto *lsi_593 = buffer.data(lsi + 593);
    const auto *lsi_597 = buffer.data(lsi + 597);
    const auto *lsi_602 = buffer.data(lsi + 602);
    const auto *lsi_609 = buffer.data(lsi + 609);
    const auto *lsi_615 = buffer.data(lsi + 615);
    const auto *lsi_616 = buffer.data(lsi + 616);
    const auto *lsi_618 = buffer.data(lsi + 618);
    const auto *lsi_749 = buffer.data(lsi + 749);
    const auto *lsi_750 = buffer.data(lsi + 750);
    const auto *lsi_751 = buffer.data(lsi + 751);
    const auto *lsi_752 = buffer.data(lsi + 752);
    const auto *lsi_753 = buffer.data(lsi + 753);
    const auto *lsi_754 = buffer.data(lsi + 754);
    const auto *lsi_755 = buffer.data(lsi + 755);
    const auto *lsi_756 = buffer.data(lsi + 756);
    const auto *lsi_761 = buffer.data(lsi + 761);
    const auto *lsi_765 = buffer.data(lsi + 765);
    const auto *lsi_770 = buffer.data(lsi + 770);
    const auto *lsi_776 = buffer.data(lsi + 776);
    const auto *lsi_777 = buffer.data(lsi + 777);
    const auto *lsi_778 = buffer.data(lsi + 778);
    const auto *lsi_779 = buffer.data(lsi + 779);
    const auto *lsi_780 = buffer.data(lsi + 780);
    const auto *lsi_781 = buffer.data(lsi + 781);
    const auto *lsi_783 = buffer.data(lsi + 783);
    const auto *lsi_784 = buffer.data(lsi + 784);
    const auto *lsi_787 = buffer.data(lsi + 787);
    const auto *lsi_790 = buffer.data(lsi + 790);
    const auto *lsi_794 = buffer.data(lsi + 794);
    const auto *lsi_799 = buffer.data(lsi + 799);
    const auto *lsi_805 = buffer.data(lsi + 805);
    const auto *lsi_807 = buffer.data(lsi + 807);
    const auto *lsi_808 = buffer.data(lsi + 808);
    const auto *lsi_809 = buffer.data(lsi + 809);
    const auto *lsi_810 = buffer.data(lsi + 810);
    const auto *lsi_811 = buffer.data(lsi + 811);

    const auto *lsk1_720 = buffer.data(lsk1 + 720);
    const auto *lsk1_723 = buffer.data(lsk1 + 723);
    const auto *lsk1_725 = buffer.data(lsk1 + 725);
    const auto *lsk1_726 = buffer.data(lsk1 + 726);
    const auto *lsk1_729 = buffer.data(lsk1 + 729);
    const auto *lsk1_730 = buffer.data(lsk1 + 730);
    const auto *lsk1_732 = buffer.data(lsk1 + 732);
    const auto *lsk1_734 = buffer.data(lsk1 + 734);
    const auto *lsk1_735 = buffer.data(lsk1 + 735);
    const auto *lsk1_737 = buffer.data(lsk1 + 737);
    const auto *lsk1_738 = buffer.data(lsk1 + 738);
    const auto *lsk1_740 = buffer.data(lsk1 + 740);
    const auto *lsk1_755 = buffer.data(lsk1 + 755);
    const auto *lsk1_756 = buffer.data(lsk1 + 756);
    const auto *lsk1_759 = buffer.data(lsk1 + 759);

    const auto *msh0_543 = buffer.data(msh0 + 543);
    const auto *msh0_544 = buffer.data(msh0 + 544);
    const auto *msh0_545 = buffer.data(msh0 + 545);
    const auto *msh0_561 = buffer.data(msh0 + 561);
    const auto *msh0_563 = buffer.data(msh0 + 563);
    const auto *msh0_564 = buffer.data(msh0 + 564);
    const auto *msh0_565 = buffer.data(msh0 + 565);
    const auto *msh0_566 = buffer.data(msh0 + 566);
    const auto *msh0_567 = buffer.data(msh0 + 567);
    const auto *msh0_568 = buffer.data(msh0 + 568);
    const auto *msh0_569 = buffer.data(msh0 + 569);
    const auto *msh0_570 = buffer.data(msh0 + 570);
    const auto *msh0_571 = buffer.data(msh0 + 571);
    const auto *msh0_572 = buffer.data(msh0 + 572);
    const auto *msh0_573 = buffer.data(msh0 + 573);
    const auto *msh0_574 = buffer.data(msh0 + 574);
    const auto *msh0_575 = buffer.data(msh0 + 575);
    const auto *msh0_576 = buffer.data(msh0 + 576);
    const auto *msh0_581 = buffer.data(msh0 + 581);
    const auto *msh0_582 = buffer.data(msh0 + 582);
    const auto *msh0_583 = buffer.data(msh0 + 583);
    const auto *msh0_584 = buffer.data(msh0 + 584);
    const auto *msh0_585 = buffer.data(msh0 + 585);
    const auto *msh0_586 = buffer.data(msh0 + 586);
    const auto *msh0_587 = buffer.data(msh0 + 587);
    const auto *msh0_588 = buffer.data(msh0 + 588);
    const auto *msh0_590 = buffer.data(msh0 + 590);
    const auto *msh0_591 = buffer.data(msh0 + 591);
    const auto *msh0_593 = buffer.data(msh0 + 593);
    const auto *msh0_594 = buffer.data(msh0 + 594);
    const auto *msh0_595 = buffer.data(msh0 + 595);
    const auto *msh0_597 = buffer.data(msh0 + 597);
    const auto *msh0_598 = buffer.data(msh0 + 598);
    const auto *msh0_603 = buffer.data(msh0 + 603);
    const auto *msh0_604 = buffer.data(msh0 + 604);
    const auto *msh0_605 = buffer.data(msh0 + 605);
    const auto *msh0_606 = buffer.data(msh0 + 606);
    const auto *msh0_608 = buffer.data(msh0 + 608);

    const auto *msh1_543 = buffer.data(msh1 + 543);
    const auto *msh1_544 = buffer.data(msh1 + 544);
    const auto *msh1_545 = buffer.data(msh1 + 545);
    const auto *msh1_561 = buffer.data(msh1 + 561);
    const auto *msh1_563 = buffer.data(msh1 + 563);
    const auto *msh1_564 = buffer.data(msh1 + 564);
    const auto *msh1_565 = buffer.data(msh1 + 565);
    const auto *msh1_566 = buffer.data(msh1 + 566);
    const auto *msh1_567 = buffer.data(msh1 + 567);
    const auto *msh1_568 = buffer.data(msh1 + 568);
    const auto *msh1_569 = buffer.data(msh1 + 569);
    const auto *msh1_570 = buffer.data(msh1 + 570);
    const auto *msh1_571 = buffer.data(msh1 + 571);
    const auto *msh1_572 = buffer.data(msh1 + 572);
    const auto *msh1_573 = buffer.data(msh1 + 573);
    const auto *msh1_574 = buffer.data(msh1 + 574);
    const auto *msh1_575 = buffer.data(msh1 + 575);
    const auto *msh1_576 = buffer.data(msh1 + 576);
    const auto *msh1_581 = buffer.data(msh1 + 581);
    const auto *msh1_582 = buffer.data(msh1 + 582);
    const auto *msh1_583 = buffer.data(msh1 + 583);
    const auto *msh1_584 = buffer.data(msh1 + 584);
    const auto *msh1_585 = buffer.data(msh1 + 585);
    const auto *msh1_586 = buffer.data(msh1 + 586);
    const auto *msh1_587 = buffer.data(msh1 + 587);
    const auto *msh1_588 = buffer.data(msh1 + 588);
    const auto *msh1_590 = buffer.data(msh1 + 590);
    const auto *msh1_591 = buffer.data(msh1 + 591);
    const auto *msh1_593 = buffer.data(msh1 + 593);
    const auto *msh1_594 = buffer.data(msh1 + 594);
    const auto *msh1_595 = buffer.data(msh1 + 595);
    const auto *msh1_597 = buffer.data(msh1 + 597);
    const auto *msh1_598 = buffer.data(msh1 + 598);
    const auto *msh1_603 = buffer.data(msh1 + 603);
    const auto *msh1_604 = buffer.data(msh1 + 604);
    const auto *msh1_605 = buffer.data(msh1 + 605);
    const auto *msh1_606 = buffer.data(msh1 + 606);
    const auto *msh1_608 = buffer.data(msh1 + 608);

    const auto *msi_724 = buffer.data(msi + 724);
    const auto *msi_725 = buffer.data(msi + 725);
    const auto *msi_726 = buffer.data(msi + 726);
    const auto *msi_727 = buffer.data(msi + 727);
    const auto *msi_728 = buffer.data(msi + 728);
    const auto *msi_730 = buffer.data(msi + 730);
    const auto *msi_731 = buffer.data(msi + 731);
    const auto *msi_733 = buffer.data(msi + 733);
    const auto *msi_734 = buffer.data(msi + 734);
    const auto *msi_737 = buffer.data(msi + 737);
    const auto *msi_738 = buffer.data(msi + 738);
    const auto *msi_742 = buffer.data(msi + 742);
    const auto *msi_749 = buffer.data(msi + 749);
    const auto *msi_750 = buffer.data(msi + 750);
    const auto *msi_751 = buffer.data(msi + 751);
    const auto *msi_752 = buffer.data(msi + 752);
    const auto *msi_753 = buffer.data(msi + 753);
    const auto *msi_754 = buffer.data(msi + 754);
    const auto *msi_755 = buffer.data(msi + 755);
    const auto *msi_756 = buffer.data(msi + 756);
    const auto *msi_757 = buffer.data(msi + 757);
    const auto *msi_758 = buffer.data(msi + 758);
    const auto *msi_759 = buffer.data(msi + 759);
    const auto *msi_760 = buffer.data(msi + 760);
    const auto *msi_761 = buffer.data(msi + 761);
    const auto *msi_762 = buffer.data(msi + 762);
    const auto *msi_763 = buffer.data(msi + 763);
    const auto *msi_764 = buffer.data(msi + 764);
    const auto *msi_765 = buffer.data(msi + 765);
    const auto *msi_766 = buffer.data(msi + 766);
    const auto *msi_767 = buffer.data(msi + 767);
    const auto *msi_768 = buffer.data(msi + 768);
    const auto *msi_769 = buffer.data(msi + 769);
    const auto *msi_770 = buffer.data(msi + 770);
    const auto *msi_776 = buffer.data(msi + 776);
    const auto *msi_777 = buffer.data(msi + 777);
    const auto *msi_778 = buffer.data(msi + 778);
    const auto *msi_779 = buffer.data(msi + 779);
    const auto *msi_780 = buffer.data(msi + 780);
    const auto *msi_781 = buffer.data(msi + 781);
    const auto *msi_782 = buffer.data(msi + 782);
    const auto *msi_783 = buffer.data(msi + 783);
    const auto *msi_784 = buffer.data(msi + 784);
    const auto *msi_785 = buffer.data(msi + 785);
    const auto *msi_786 = buffer.data(msi + 786);
    const auto *msi_787 = buffer.data(msi + 787);
    const auto *msi_789 = buffer.data(msi + 789);
    const auto *msi_790 = buffer.data(msi + 790);
    const auto *msi_791 = buffer.data(msi + 791);
    const auto *msi_793 = buffer.data(msi + 793);
    const auto *msi_794 = buffer.data(msi + 794);
    const auto *msi_795 = buffer.data(msi + 795);
    const auto *msi_796 = buffer.data(msi + 796);
    const auto *msi_798 = buffer.data(msi + 798);
    const auto *msi_799 = buffer.data(msi + 799);
    const auto *msi_805 = buffer.data(msi + 805);
    const auto *msi_806 = buffer.data(msi + 806);
    const auto *msi_807 = buffer.data(msi + 807);
    const auto *msi_808 = buffer.data(msi + 808);
    const auto *msi_809 = buffer.data(msi + 809);
    const auto *msi_810 = buffer.data(msi + 810);
    const auto *msi_811 = buffer.data(msi + 811);
    const auto *msi_812 = buffer.data(msi + 812);
    const auto *msi_814 = buffer.data(msi + 814);

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, lsi_556, lsi_557, lsi_558, msh0_543, \
                         msh0_544, msh0_545, msh1_543, msh1_544, msh1_545, msi_724, msi_725, \
                         msi_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * lsi_556[k]
                   + f_8 * msh0_543[k]
                   - f_9 * msh1_543[k]
                   + f_3 * pc_y[k] * msi_724[k];

        t_932[k] = f_14 * lsi_557[k]
                   + f_6 * msh0_544[k]
                   - f_7 * msh1_544[k]
                   + f_3 * pc_y[k] * msi_725[k];

        t_933[k] = f_14 * lsi_558[k]
                   + f_4 * msh0_545[k]
                   - f_5 * msh1_545[k]
                   + f_3 * pc_y[k] * msi_726[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pa_y, pc_y, pc_z, lsk0_720, lsi_531, \
                         lsi_559, lsi_560, lsk1_720, msh0_545, msh1_545, msi_727, \
                         msi_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * lsi_559[k]
                   + f_3 * pc_y[k] * msi_727[k];

        t_935[k] = f_16 * lsi_531[k]
                   + f_1 * msh0_545[k]
                   - f_2 * msh1_545[k]
                   + f_3 * pc_z[k] * msi_727[k];

        t_936[k] = pa_y[k] * lsk0_720[k]
                   - f_12 * pc_y[k] * lsk1_720[k];

        t_937[k] = f_13 * lsi_560[k]
                   + f_3 * pc_y[k] * msi_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pa_y, pc_y, pc_z, lsk0_723, lsk0_725, \
                         lsi_532, lsi_561, lsi_562, lsk1_723, lsk1_725, msi_728, \
                         msi_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * lsi_532[k]
                   + f_3 * pc_z[k] * msi_728[k];

        t_939[k] = pa_y[k] * lsk0_723[k]
                   + f_14 * lsi_561[k]
                   - f_12 * pc_y[k] * lsk1_723[k];

        t_940[k] = f_13 * lsi_562[k]
                   + f_3 * pc_y[k] * msi_730[k];

        t_941[k] = pa_y[k] * lsk0_725[k]
                   - f_12 * pc_y[k] * lsk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_y, pc_y, pc_z, lsk0_726, lsk0_729, \
                         lsi_535, lsi_563, lsi_565, lsk1_726, lsk1_729, msi_731, \
                         msi_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pa_y[k] * lsk0_726[k]
                   + f_15 * lsi_563[k]
                   - f_12 * pc_y[k] * lsk1_726[k];

        t_943[k] = f_17 * lsi_535[k]
                   + f_3 * pc_z[k] * msi_731[k];

        t_944[k] = f_13 * lsi_565[k]
                   + f_3 * pc_y[k] * msi_733[k];

        t_945[k] = pa_y[k] * lsk0_729[k]
                   - f_12 * pc_y[k] * lsk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pa_y, pc_y, pc_z, lsk0_730, lsk0_732, lsi_538, \
                         lsi_566, lsi_568, lsk1_730, lsk1_732, \
                         msi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pa_y[k] * lsk0_730[k]
                   + f_16 * lsi_566[k]
                   - f_12 * pc_y[k] * lsk1_730[k];

        t_947[k] = f_17 * lsi_538[k]
                   + f_3 * pc_z[k] * msi_734[k];

        t_948[k] = pa_y[k] * lsk0_732[k]
                   + f_14 * lsi_568[k]
                   - f_12 * pc_y[k] * lsk1_732[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, lsk0_734, lsk0_735, \
                         lsi_542, lsi_569, lsi_570, lsk1_734, lsk1_735, msi_737, \
                         msi_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_13 * lsi_569[k]
                   + f_3 * pc_y[k] * msi_737[k];

        t_950[k] = pa_y[k] * lsk0_734[k]
                   - f_12 * pc_y[k] * lsk1_734[k];

        t_951[k] = pa_y[k] * lsk0_735[k]
                   + f_17 * lsi_570[k]
                   - f_12 * pc_y[k] * lsk1_735[k];

        t_952[k] = f_17 * lsi_542[k]
                   + f_3 * pc_z[k] * msi_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, lsk0_737, lsk0_738, lsk0_740, \
                         lsi_572, lsi_573, lsi_574, lsk1_737, lsk1_738, lsk1_740, \
                         msi_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pa_y[k] * lsk0_737[k]
                   + f_15 * lsi_572[k]
                   - f_12 * pc_y[k] * lsk1_737[k];

        t_954[k] = pa_y[k] * lsk0_738[k]
                   + f_14 * lsi_573[k]
                   - f_12 * pc_y[k] * lsk1_738[k];

        t_955[k] = f_13 * lsi_574[k]
                   + f_3 * pc_y[k] * msi_742[k];

        t_956[k] = pa_y[k] * lsk0_740[k]
                   - f_12 * pc_y[k] * lsk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, lsi_749, lsi_750, lsi_751, \
                         lsi_752, lsi_753, msi_749, msi_750, msi_751, msi_752, \
                         msi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_15 * lsi_749[k]
                   + f_3 * pc_x[k] * msi_749[k];

        t_958[k] = f_15 * lsi_750[k]
                   + f_3 * pc_x[k] * msi_750[k];

        t_959[k] = f_15 * lsi_751[k]
                   + f_3 * pc_x[k] * msi_751[k];

        t_960[k] = f_15 * lsi_752[k]
                   + f_3 * pc_x[k] * msi_752[k];

        t_961[k] = f_15 * lsi_753[k]
                   + f_3 * pc_x[k] * msi_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pc_x, pc_y, pc_z, lsi_553, lsi_581, \
                         lsi_754, lsi_755, msh0_561, msh1_561, msi_749, msi_754, \
                         msi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_15 * lsi_754[k]
                   + f_3 * pc_x[k] * msi_754[k];

        t_963[k] = f_15 * lsi_755[k]
                   + f_3 * pc_x[k] * msi_755[k];

        t_964[k] = f_13 * lsi_581[k]
                   + f_1 * msh0_561[k]
                   - f_2 * msh1_561[k]
                   + f_3 * pc_y[k] * msi_749[k];

        t_965[k] = f_17 * lsi_553[k]
                   + f_3 * pc_z[k] * msi_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_y, lsi_583, lsi_584, lsi_585, msh0_563, \
                         msh0_564, msh0_565, msh1_563, msh1_564, msh1_565, msi_751, msi_752, \
                         msi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_13 * lsi_583[k]
                   + f_10 * msh0_563[k]
                   - f_11 * msh1_563[k]
                   + f_3 * pc_y[k] * msi_751[k];

        t_967[k] = f_13 * lsi_584[k]
                   + f_8 * msh0_564[k]
                   - f_9 * msh1_564[k]
                   + f_3 * pc_y[k] * msi_752[k];

        t_968[k] = f_13 * lsi_585[k]
                   + f_6 * msh0_565[k]
                   - f_7 * msh1_565[k]
                   + f_3 * pc_y[k] * msi_753[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_y, pc_y, lsk0_755, lsi_586, lsi_587, \
                         lsk1_755, msh0_566, msh1_566, msi_754, \
                         msi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_13 * lsi_586[k]
                   + f_4 * msh0_566[k]
                   - f_5 * msh1_566[k]
                   + f_3 * pc_y[k] * msi_754[k];

        t_970[k] = f_13 * lsi_587[k]
                   + f_3 * pc_y[k] * msi_755[k];

        t_971[k] = pa_y[k] * lsk0_755[k]
                   - f_12 * pc_y[k] * lsk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, pc_x, pc_y, pc_z, lsi_560, \
                         lsi_756, msh0_567, msh1_567, msi_756, msi_757, \
                         msi_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_15 * lsi_756[k]
                   + f_1 * msh0_567[k]
                   - f_2 * msh1_567[k]
                   + f_3 * pc_x[k] * msi_756[k];

        t_973[k] = f_3 * pc_y[k] * msi_756[k];

        t_974[k] = f_22 * lsi_560[k]
                   + f_3 * pc_z[k] * msi_756[k];

        t_975[k] = f_4 * msh0_567[k]
                   - f_5 * msh1_567[k]
                   + f_3 * pc_y[k] * msi_757[k];

        t_976[k] = f_3 * pc_y[k] * msi_758[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, pc_x, pc_y, lsi_761, msh0_568, msh0_569, \
                         msh0_572, msh1_568, msh1_569, msh1_572, msi_759, msi_760, \
                         msi_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_15 * lsi_761[k]
                   + f_10 * msh0_572[k]
                   - f_11 * msh1_572[k]
                   + f_3 * pc_x[k] * msi_761[k];

        t_978[k] = f_6 * msh0_568[k]
                   - f_7 * msh1_568[k]
                   + f_3 * pc_y[k] * msi_759[k];

        t_979[k] = f_4 * msh0_569[k]
                   - f_5 * msh1_569[k]
                   + f_3 * pc_y[k] * msi_760[k];

        t_980[k] = f_3 * pc_y[k] * msi_761[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_x, pc_y, lsi_765, msh0_570, msh0_571, \
                         msh0_576, msh1_570, msh1_571, msh1_576, msi_762, msi_763, \
                         msi_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_15 * lsi_765[k]
                   + f_8 * msh0_576[k]
                   - f_9 * msh1_576[k]
                   + f_3 * pc_x[k] * msi_765[k];

        t_982[k] = f_8 * msh0_570[k]
                   - f_9 * msh1_570[k]
                   + f_3 * pc_y[k] * msi_762[k];

        t_983[k] = f_6 * msh0_571[k]
                   - f_7 * msh1_571[k]
                   + f_3 * pc_y[k] * msi_763[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_x, pc_y, lsi_770, msh0_572, msh0_581, \
                         msh1_572, msh1_581, msi_764, msi_765, \
                         msi_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_4 * msh0_572[k]
                   - f_5 * msh1_572[k]
                   + f_3 * pc_y[k] * msi_764[k];

        t_985[k] = f_3 * pc_y[k] * msi_765[k];

        t_986[k] = f_15 * lsi_770[k]
                   + f_6 * msh0_581[k]
                   - f_7 * msh1_581[k]
                   + f_3 * pc_x[k] * msi_770[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pc_y, msh0_573, msh0_574, msh0_575, msh1_573, \
                         msh1_574, msh1_575, msi_766, msi_767, \
                         msi_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_10 * msh0_573[k]
                   - f_11 * msh1_573[k]
                   + f_3 * pc_y[k] * msi_766[k];

        t_988[k] = f_8 * msh0_574[k]
                   - f_9 * msh1_574[k]
                   + f_3 * pc_y[k] * msi_767[k];

        t_989[k] = f_6 * msh0_575[k]
                   - f_7 * msh1_575[k]
                   + f_3 * pc_y[k] * msi_768[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, pc_x, pc_y, lsi_776, lsi_777, msh0_576, \
                         msh0_587, msh1_576, msh1_587, msi_769, msi_770, msi_776, \
                         msi_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = f_4 * msh0_576[k]
                   - f_5 * msh1_576[k]
                   + f_3 * pc_y[k] * msi_769[k];

        t_991[k] = f_3 * pc_y[k] * msi_770[k];

        t_992[k] = f_15 * lsi_776[k]
                   + f_4 * msh0_587[k]
                   - f_5 * msh1_587[k]
                   + f_3 * pc_x[k] * msi_776[k];

        t_993[k] = f_15 * lsi_777[k]
                   + f_3 * pc_x[k] * msi_777[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, t_998, pc_x, pc_y, lsi_778, lsi_779, \
                         lsi_780, lsi_781, msi_776, msi_778, msi_779, msi_780, \
                         msi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_15 * lsi_778[k]
                   + f_3 * pc_x[k] * msi_778[k];

        t_995[k] = f_15 * lsi_779[k]
                   + f_3 * pc_x[k] * msi_779[k];

        t_996[k] = f_15 * lsi_780[k]
                   + f_3 * pc_x[k] * msi_780[k];

        t_997[k] = f_15 * lsi_781[k]
                   + f_3 * pc_x[k] * msi_781[k];

        t_998[k] = f_3 * pc_y[k] * msi_776[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, pc_x, pc_y, lsi_783, msh0_582, msh0_583, \
                         msh1_582, msh1_583, msi_777, msi_778, \
                         msi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_15 * lsi_783[k]
                   + f_3 * pc_x[k] * msi_783[k];

        t_1000[k] = f_1 * msh0_582[k]
                    - f_2 * msh1_582[k]
                    + f_3 * pc_y[k] * msi_777[k];

        t_1001[k] = f_19 * msh0_583[k]
                    - f_20 * msh1_583[k]
                    + f_3 * pc_y[k] * msi_778[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, msh0_584, msh0_585, msh0_586, msh1_584, \
                         msh1_585, msh1_586, msi_779, msi_780, \
                         msi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_10 * msh0_584[k]
                    - f_11 * msh1_584[k]
                    + f_3 * pc_y[k] * msi_779[k];

        t_1003[k] = f_8 * msh0_585[k]
                    - f_9 * msh1_585[k]
                    + f_3 * pc_y[k] * msi_780[k];

        t_1004[k] = f_6 * msh0_586[k]
                    - f_7 * msh1_586[k]
                    + f_3 * pc_y[k] * msi_781[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, pc_x, pc_y, pc_z, lsi_587, lsi_784, \
                         msh0_587, msh0_588, msh1_587, msh1_588, msi_782, msi_783, \
                         msi_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_4 * msh0_587[k]
                    - f_5 * msh1_587[k]
                    + f_3 * pc_y[k] * msi_782[k];

        t_1006[k] = f_3 * pc_y[k] * msi_783[k];

        t_1007[k] = f_22 * lsi_587[k]
                    + f_1 * msh0_587[k]
                    - f_2 * msh1_587[k]
                    + f_3 * pc_z[k] * msi_783[k];

        t_1008[k] = f_14 * lsi_784[k]
                    + f_1 * msh0_588[k]
                    - f_2 * msh1_588[k]
                    + f_3 * pc_x[k] * msi_784[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pc_x, pc_y, pc_z, lsi_588, lsi_787, \
                         msh0_591, msh1_591, msi_784, msi_785, \
                         msi_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_21 * lsi_588[k]
                    + f_3 * pc_y[k] * msi_784[k];

        t_1010[k] = f_3 * pc_z[k] * msi_784[k];

        t_1011[k] = f_14 * lsi_787[k]
                    + f_10 * msh0_591[k]
                    - f_11 * msh1_591[k]
                    + f_3 * pc_x[k] * msi_787[k];

        t_1012[k] = f_3 * pc_z[k] * msi_785[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, pc_z, lsi_790, msh0_588, msh0_594, \
                         msh1_588, msh1_594, msi_786, msi_787, \
                         msi_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_4 * msh0_588[k]
                    - f_5 * msh1_588[k]
                    + f_3 * pc_z[k] * msi_786[k];

        t_1014[k] = f_14 * lsi_790[k]
                    + f_8 * msh0_594[k]
                    - f_9 * msh1_594[k]
                    + f_3 * pc_x[k] * msi_790[k];

        t_1015[k] = f_3 * pc_z[k] * msi_787[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, pc_z, lsi_593, lsi_794, \
                         msh0_590, msh0_598, msh1_590, msh1_598, msi_789, msi_790, \
                         msi_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_21 * lsi_593[k]
                    + f_3 * pc_y[k] * msi_789[k];

        t_1017[k] = f_6 * msh0_590[k]
                    - f_7 * msh1_590[k]
                    + f_3 * pc_z[k] * msi_789[k];

        t_1018[k] = f_14 * lsi_794[k]
                    + f_6 * msh0_598[k]
                    - f_7 * msh1_598[k]
                    + f_3 * pc_x[k] * msi_794[k];

        t_1019[k] = f_3 * pc_z[k] * msi_790[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, pc_y, pc_z, lsi_597, msh0_591, msh0_593, \
                         msh1_591, msh1_593, msi_791, msi_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * msh0_591[k]
                    - f_5 * msh1_591[k]
                    + f_3 * pc_z[k] * msi_791[k];

        t_1021[k] = f_21 * lsi_597[k]
                    + f_3 * pc_y[k] * msi_793[k];

        t_1022[k] = f_8 * msh0_593[k]
                    - f_9 * msh1_593[k]
                    + f_3 * pc_z[k] * msi_793[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, pc_x, pc_z, lsi_799, msh0_594, msh0_603, \
                         msh1_594, msh1_603, msi_794, msi_795, \
                         msi_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_14 * lsi_799[k]
                    + f_4 * msh0_603[k]
                    - f_5 * msh1_603[k]
                    + f_3 * pc_x[k] * msi_799[k];

        t_1024[k] = f_3 * pc_z[k] * msi_794[k];

        t_1025[k] = f_4 * msh0_594[k]
                    - f_5 * msh1_594[k]
                    + f_3 * pc_z[k] * msi_795[k];
    }

#pragma omp simd aligned(t_1026, t_1027, t_1028, t_1029, pc_x, pc_y, pc_z, lsi_602, lsi_805, \
                         msh0_595, msh0_597, msh1_595, msh1_597, msi_796, msi_798, \
                         msi_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1026[k] = f_6 * msh0_595[k]
                    - f_7 * msh1_595[k]
                    + f_3 * pc_z[k] * msi_796[k];

        t_1027[k] = f_21 * lsi_602[k]
                    + f_3 * pc_y[k] * msi_798[k];

        t_1028[k] = f_10 * msh0_597[k]
                    - f_11 * msh1_597[k]
                    + f_3 * pc_z[k] * msi_798[k];

        t_1029[k] = f_14 * lsi_805[k]
                    + f_3 * pc_x[k] * msi_805[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, t_1034, pc_x, pc_z, lsi_807, lsi_808, \
                         lsi_809, lsi_810, msi_799, msi_807, msi_808, msi_809, \
                         msi_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_3 * pc_z[k] * msi_799[k];

        t_1031[k] = f_14 * lsi_807[k]
                    + f_3 * pc_x[k] * msi_807[k];

        t_1032[k] = f_14 * lsi_808[k]
                    + f_3 * pc_x[k] * msi_808[k];

        t_1033[k] = f_14 * lsi_809[k]
                    + f_3 * pc_x[k] * msi_809[k];

        t_1034[k] = f_14 * lsi_810[k]
                    + f_3 * pc_x[k] * msi_810[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, pc_x, pc_y, pc_z, lsi_609, lsi_811, \
                         msh0_603, msh1_603, msi_805, msi_806, \
                         msi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_14 * lsi_811[k]
                    + f_3 * pc_x[k] * msi_811[k];

        t_1036[k] = f_21 * lsi_609[k]
                    + f_1 * msh0_603[k]
                    - f_2 * msh1_603[k]
                    + f_3 * pc_y[k] * msi_805[k];

        t_1037[k] = f_3 * pc_z[k] * msi_805[k];

        t_1038[k] = f_4 * msh0_603[k]
                    - f_5 * msh1_603[k]
                    + f_3 * pc_z[k] * msi_806[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_z, msh0_604, msh0_605, msh0_606, msh1_604, \
                         msh1_605, msh1_606, msi_807, msi_808, \
                         msi_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_6 * msh0_604[k]
                    - f_7 * msh1_604[k]
                    + f_3 * pc_z[k] * msi_807[k];

        t_1040[k] = f_8 * msh0_605[k]
                    - f_9 * msh1_605[k]
                    + f_3 * pc_z[k] * msi_808[k];

        t_1041[k] = f_10 * msh0_606[k]
                    - f_11 * msh1_606[k]
                    + f_3 * pc_z[k] * msi_809[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pa_z, pc_y, pc_z, lsk0_756, lsi_615, \
                         lsi_616, lsk1_756, msh0_608, msh1_608, msi_811, \
                         msi_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_21 * lsi_615[k]
                    + f_3 * pc_y[k] * msi_811[k];

        t_1043[k] = f_1 * msh0_608[k]
                    - f_2 * msh1_608[k]
                    + f_3 * pc_z[k] * msi_811[k];

        t_1044[k] = pa_z[k] * lsk0_756[k]
                    - f_12 * pc_z[k] * lsk1_756[k];

        t_1045[k] = f_22 * lsi_616[k]
                    + f_3 * pc_y[k] * msi_812[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pa_z, pc_y, pc_z, lsk0_759, lsi_588, lsi_618, \
                         lsk1_759, msi_812, msi_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_13 * lsi_588[k]
                    + f_3 * pc_z[k] * msi_812[k];

        t_1047[k] = pa_z[k] * lsk0_759[k]
                    - f_12 * pc_z[k] * lsk1_759[k];

        t_1048[k] = f_22 * lsi_618[k]
                    + f_3 * pc_y[k] * msi_814[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsk0,
                                                          const size_t lsi, const size_t lsk1,
                                                          const size_t msh0, const size_t msh1,
                                                          const size_t msi, const size_t ncols,
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
    const auto f_22 = 3.0 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_762 = buffer.data(lsk0 + 762);
    const auto *lsk0_766 = buffer.data(lsk0 + 766);
    const auto *lsk0_768 = buffer.data(lsk0 + 768);
    const auto *lsk0_771 = buffer.data(lsk0 + 771);
    const auto *lsk0_773 = buffer.data(lsk0 + 773);
    const auto *lsk0_774 = buffer.data(lsk0 + 774);
    const auto *lsk0_784 = buffer.data(lsk0 + 784);

    const auto *lsi_591 = buffer.data(lsi + 591);
    const auto *lsi_594 = buffer.data(lsi + 594);
    const auto *lsi_595 = buffer.data(lsi + 595);
    const auto *lsi_598 = buffer.data(lsi + 598);
    const auto *lsi_599 = buffer.data(lsi + 599);
    const auto *lsi_600 = buffer.data(lsi + 600);
    const auto *lsi_609 = buffer.data(lsi + 609);
    const auto *lsi_615 = buffer.data(lsi + 615);
    const auto *lsi_616 = buffer.data(lsi + 616);
    const auto *lsi_619 = buffer.data(lsi + 619);
    const auto *lsi_621 = buffer.data(lsi + 621);
    const auto *lsi_622 = buffer.data(lsi + 622);
    const auto *lsi_625 = buffer.data(lsi + 625);
    const auto *lsi_626 = buffer.data(lsi + 626);
    const auto *lsi_630 = buffer.data(lsi + 630);
    const auto *lsi_637 = buffer.data(lsi + 637);
    const auto *lsi_639 = buffer.data(lsi + 639);
    const auto *lsi_640 = buffer.data(lsi + 640);
    const auto *lsi_641 = buffer.data(lsi + 641);
    const auto *lsi_642 = buffer.data(lsi + 642);
    const auto *lsi_643 = buffer.data(lsi + 643);
    const auto *lsi_644 = buffer.data(lsi + 644);
    const auto *lsi_646 = buffer.data(lsi + 646);
    const auto *lsi_647 = buffer.data(lsi + 647);
    const auto *lsi_649 = buffer.data(lsi + 649);
    const auto *lsi_650 = buffer.data(lsi + 650);
    const auto *lsi_653 = buffer.data(lsi + 653);
    const auto *lsi_654 = buffer.data(lsi + 654);
    const auto *lsi_658 = buffer.data(lsi + 658);
    const auto *lsi_665 = buffer.data(lsi + 665);
    const auto *lsi_667 = buffer.data(lsi + 667);
    const auto *lsi_668 = buffer.data(lsi + 668);
    const auto *lsi_669 = buffer.data(lsi + 669);
    const auto *lsi_670 = buffer.data(lsi + 670);
    const auto *lsi_671 = buffer.data(lsi + 671);
    const auto *lsi_672 = buffer.data(lsi + 672);
    const auto *lsi_674 = buffer.data(lsi + 674);
    const auto *lsi_677 = buffer.data(lsi + 677);
    const auto *lsi_681 = buffer.data(lsi + 681);
    const auto *lsi_686 = buffer.data(lsi + 686);
    const auto *lsi_693 = buffer.data(lsi + 693);
    const auto *lsi_695 = buffer.data(lsi + 695);
    const auto *lsi_696 = buffer.data(lsi + 696);
    const auto *lsi_697 = buffer.data(lsi + 697);
    const auto *lsi_698 = buffer.data(lsi + 698);
    const auto *lsi_699 = buffer.data(lsi + 699);
    const auto *lsi_817 = buffer.data(lsi + 817);
    const auto *lsi_821 = buffer.data(lsi + 821);
    const auto *lsi_826 = buffer.data(lsi + 826);
    const auto *lsi_832 = buffer.data(lsi + 832);
    const auto *lsi_833 = buffer.data(lsi + 833);
    const auto *lsi_834 = buffer.data(lsi + 834);
    const auto *lsi_835 = buffer.data(lsi + 835);
    const auto *lsi_836 = buffer.data(lsi + 836);
    const auto *lsi_837 = buffer.data(lsi + 837);
    const auto *lsi_838 = buffer.data(lsi + 838);
    const auto *lsi_839 = buffer.data(lsi + 839);
    const auto *lsi_840 = buffer.data(lsi + 840);
    const auto *lsi_843 = buffer.data(lsi + 843);
    const auto *lsi_845 = buffer.data(lsi + 845);
    const auto *lsi_846 = buffer.data(lsi + 846);
    const auto *lsi_849 = buffer.data(lsi + 849);
    const auto *lsi_850 = buffer.data(lsi + 850);
    const auto *lsi_852 = buffer.data(lsi + 852);
    const auto *lsi_854 = buffer.data(lsi + 854);
    const auto *lsi_855 = buffer.data(lsi + 855);
    const auto *lsi_857 = buffer.data(lsi + 857);
    const auto *lsi_858 = buffer.data(lsi + 858);
    const auto *lsi_860 = buffer.data(lsi + 860);
    const auto *lsi_861 = buffer.data(lsi + 861);
    const auto *lsi_862 = buffer.data(lsi + 862);
    const auto *lsi_863 = buffer.data(lsi + 863);
    const auto *lsi_864 = buffer.data(lsi + 864);
    const auto *lsi_865 = buffer.data(lsi + 865);
    const auto *lsi_866 = buffer.data(lsi + 866);
    const auto *lsi_867 = buffer.data(lsi + 867);
    const auto *lsi_868 = buffer.data(lsi + 868);
    const auto *lsi_871 = buffer.data(lsi + 871);
    const auto *lsi_873 = buffer.data(lsi + 873);
    const auto *lsi_874 = buffer.data(lsi + 874);
    const auto *lsi_877 = buffer.data(lsi + 877);
    const auto *lsi_878 = buffer.data(lsi + 878);
    const auto *lsi_880 = buffer.data(lsi + 880);
    const auto *lsi_882 = buffer.data(lsi + 882);
    const auto *lsi_883 = buffer.data(lsi + 883);
    const auto *lsi_885 = buffer.data(lsi + 885);
    const auto *lsi_886 = buffer.data(lsi + 886);
    const auto *lsi_888 = buffer.data(lsi + 888);
    const auto *lsi_889 = buffer.data(lsi + 889);
    const auto *lsi_890 = buffer.data(lsi + 890);
    const auto *lsi_891 = buffer.data(lsi + 891);
    const auto *lsi_892 = buffer.data(lsi + 892);
    const auto *lsi_893 = buffer.data(lsi + 893);
    const auto *lsi_894 = buffer.data(lsi + 894);
    const auto *lsi_895 = buffer.data(lsi + 895);
    const auto *lsi_896 = buffer.data(lsi + 896);

    const auto *lsk1_762 = buffer.data(lsk1 + 762);
    const auto *lsk1_766 = buffer.data(lsk1 + 766);
    const auto *lsk1_768 = buffer.data(lsk1 + 768);
    const auto *lsk1_771 = buffer.data(lsk1 + 771);
    const auto *lsk1_773 = buffer.data(lsk1 + 773);
    const auto *lsk1_774 = buffer.data(lsk1 + 774);
    const auto *lsk1_784 = buffer.data(lsk1 + 784);

    const auto *msh0_614 = buffer.data(msh0 + 614);
    const auto *msh0_618 = buffer.data(msh0 + 618);
    const auto *msh0_623 = buffer.data(msh0 + 623);
    const auto *msh0_626 = buffer.data(msh0 + 626);
    const auto *msh0_627 = buffer.data(msh0 + 627);
    const auto *msh0_628 = buffer.data(msh0 + 628);
    const auto *msh0_629 = buffer.data(msh0 + 629);
    const auto *msh0_630 = buffer.data(msh0 + 630);
    const auto *msh0_633 = buffer.data(msh0 + 633);
    const auto *msh0_635 = buffer.data(msh0 + 635);
    const auto *msh0_636 = buffer.data(msh0 + 636);
    const auto *msh0_639 = buffer.data(msh0 + 639);
    const auto *msh0_640 = buffer.data(msh0 + 640);
    const auto *msh0_642 = buffer.data(msh0 + 642);
    const auto *msh0_644 = buffer.data(msh0 + 644);
    const auto *msh0_645 = buffer.data(msh0 + 645);
    const auto *msh0_647 = buffer.data(msh0 + 647);
    const auto *msh0_648 = buffer.data(msh0 + 648);
    const auto *msh0_649 = buffer.data(msh0 + 649);
    const auto *msh0_650 = buffer.data(msh0 + 650);
    const auto *msh0_651 = buffer.data(msh0 + 651);
    const auto *msh0_654 = buffer.data(msh0 + 654);
    const auto *msh0_656 = buffer.data(msh0 + 656);
    const auto *msh0_657 = buffer.data(msh0 + 657);
    const auto *msh0_660 = buffer.data(msh0 + 660);
    const auto *msh0_661 = buffer.data(msh0 + 661);
    const auto *msh0_663 = buffer.data(msh0 + 663);
    const auto *msh0_665 = buffer.data(msh0 + 665);
    const auto *msh0_666 = buffer.data(msh0 + 666);
    const auto *msh0_668 = buffer.data(msh0 + 668);
    const auto *msh0_669 = buffer.data(msh0 + 669);
    const auto *msh0_670 = buffer.data(msh0 + 670);
    const auto *msh0_671 = buffer.data(msh0 + 671);
    const auto *msh0_672 = buffer.data(msh0 + 672);

    const auto *msh1_614 = buffer.data(msh1 + 614);
    const auto *msh1_618 = buffer.data(msh1 + 618);
    const auto *msh1_623 = buffer.data(msh1 + 623);
    const auto *msh1_626 = buffer.data(msh1 + 626);
    const auto *msh1_627 = buffer.data(msh1 + 627);
    const auto *msh1_628 = buffer.data(msh1 + 628);
    const auto *msh1_629 = buffer.data(msh1 + 629);
    const auto *msh1_630 = buffer.data(msh1 + 630);
    const auto *msh1_633 = buffer.data(msh1 + 633);
    const auto *msh1_635 = buffer.data(msh1 + 635);
    const auto *msh1_636 = buffer.data(msh1 + 636);
    const auto *msh1_639 = buffer.data(msh1 + 639);
    const auto *msh1_640 = buffer.data(msh1 + 640);
    const auto *msh1_642 = buffer.data(msh1 + 642);
    const auto *msh1_644 = buffer.data(msh1 + 644);
    const auto *msh1_645 = buffer.data(msh1 + 645);
    const auto *msh1_647 = buffer.data(msh1 + 647);
    const auto *msh1_648 = buffer.data(msh1 + 648);
    const auto *msh1_649 = buffer.data(msh1 + 649);
    const auto *msh1_650 = buffer.data(msh1 + 650);
    const auto *msh1_651 = buffer.data(msh1 + 651);
    const auto *msh1_654 = buffer.data(msh1 + 654);
    const auto *msh1_656 = buffer.data(msh1 + 656);
    const auto *msh1_657 = buffer.data(msh1 + 657);
    const auto *msh1_660 = buffer.data(msh1 + 660);
    const auto *msh1_661 = buffer.data(msh1 + 661);
    const auto *msh1_663 = buffer.data(msh1 + 663);
    const auto *msh1_665 = buffer.data(msh1 + 665);
    const auto *msh1_666 = buffer.data(msh1 + 666);
    const auto *msh1_668 = buffer.data(msh1 + 668);
    const auto *msh1_669 = buffer.data(msh1 + 669);
    const auto *msh1_670 = buffer.data(msh1 + 670);
    const auto *msh1_671 = buffer.data(msh1 + 671);
    const auto *msh1_672 = buffer.data(msh1 + 672);

    const auto *msi_815 = buffer.data(msi + 815);
    const auto *msi_817 = buffer.data(msi + 817);
    const auto *msi_818 = buffer.data(msi + 818);
    const auto *msi_821 = buffer.data(msi + 821);
    const auto *msi_822 = buffer.data(msi + 822);
    const auto *msi_826 = buffer.data(msi + 826);
    const auto *msi_832 = buffer.data(msi + 832);
    const auto *msi_833 = buffer.data(msi + 833);
    const auto *msi_834 = buffer.data(msi + 834);
    const auto *msi_835 = buffer.data(msi + 835);
    const auto *msi_836 = buffer.data(msi + 836);
    const auto *msi_837 = buffer.data(msi + 837);
    const auto *msi_838 = buffer.data(msi + 838);
    const auto *msi_839 = buffer.data(msi + 839);
    const auto *msi_840 = buffer.data(msi + 840);
    const auto *msi_842 = buffer.data(msi + 842);
    const auto *msi_843 = buffer.data(msi + 843);
    const auto *msi_845 = buffer.data(msi + 845);
    const auto *msi_846 = buffer.data(msi + 846);
    const auto *msi_849 = buffer.data(msi + 849);
    const auto *msi_850 = buffer.data(msi + 850);
    const auto *msi_852 = buffer.data(msi + 852);
    const auto *msi_854 = buffer.data(msi + 854);
    const auto *msi_855 = buffer.data(msi + 855);
    const auto *msi_857 = buffer.data(msi + 857);
    const auto *msi_858 = buffer.data(msi + 858);
    const auto *msi_860 = buffer.data(msi + 860);
    const auto *msi_861 = buffer.data(msi + 861);
    const auto *msi_862 = buffer.data(msi + 862);
    const auto *msi_863 = buffer.data(msi + 863);
    const auto *msi_864 = buffer.data(msi + 864);
    const auto *msi_865 = buffer.data(msi + 865);
    const auto *msi_866 = buffer.data(msi + 866);
    const auto *msi_867 = buffer.data(msi + 867);
    const auto *msi_868 = buffer.data(msi + 868);
    const auto *msi_870 = buffer.data(msi + 870);
    const auto *msi_871 = buffer.data(msi + 871);
    const auto *msi_873 = buffer.data(msi + 873);
    const auto *msi_874 = buffer.data(msi + 874);
    const auto *msi_877 = buffer.data(msi + 877);
    const auto *msi_878 = buffer.data(msi + 878);
    const auto *msi_880 = buffer.data(msi + 880);
    const auto *msi_882 = buffer.data(msi + 882);
    const auto *msi_883 = buffer.data(msi + 883);
    const auto *msi_885 = buffer.data(msi + 885);
    const auto *msi_886 = buffer.data(msi + 886);
    const auto *msi_888 = buffer.data(msi + 888);
    const auto *msi_889 = buffer.data(msi + 889);
    const auto *msi_890 = buffer.data(msi + 890);
    const auto *msi_891 = buffer.data(msi + 891);
    const auto *msi_892 = buffer.data(msi + 892);
    const auto *msi_893 = buffer.data(msi + 893);
    const auto *msi_894 = buffer.data(msi + 894);
    const auto *msi_895 = buffer.data(msi + 895);
    const auto *msi_896 = buffer.data(msi + 896);

#pragma omp simd aligned(t_1049, t_1050, t_1051, pa_z, pc_x, pc_z, lsk0_762, lsi_591, lsi_817, \
                         lsk1_762, msh0_614, msh1_614, msi_815, \
                         msi_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_14 * lsi_817[k]
                    + f_10 * msh0_614[k]
                    - f_11 * msh1_614[k]
                    + f_3 * pc_x[k] * msi_817[k];

        t_1050[k] = pa_z[k] * lsk0_762[k]
                    - f_12 * pc_z[k] * lsk1_762[k];

        t_1051[k] = f_13 * lsi_591[k]
                    + f_3 * pc_z[k] * msi_815[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pa_z, pc_x, pc_y, pc_z, lsk0_766, lsi_621, \
                         lsi_821, lsk1_766, msh0_618, msh1_618, msi_817, \
                         msi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_22 * lsi_621[k]
                    + f_3 * pc_y[k] * msi_817[k];

        t_1053[k] = f_14 * lsi_821[k]
                    + f_8 * msh0_618[k]
                    - f_9 * msh1_618[k]
                    + f_3 * pc_x[k] * msi_821[k];

        t_1054[k] = pa_z[k] * lsk0_766[k]
                    - f_12 * pc_z[k] * lsk1_766[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pa_z, pc_y, pc_z, lsk0_768, lsi_594, lsi_595, \
                         lsi_625, lsk1_768, msi_818, msi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_13 * lsi_594[k]
                    + f_3 * pc_z[k] * msi_818[k];

        t_1056[k] = pa_z[k] * lsk0_768[k]
                    + f_14 * lsi_595[k]
                    - f_12 * pc_z[k] * lsk1_768[k];

        t_1057[k] = f_22 * lsi_625[k]
                    + f_3 * pc_y[k] * msi_821[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pa_z, pc_x, pc_z, lsk0_771, lsi_598, lsi_826, \
                         lsk1_771, msh0_623, msh1_623, msi_822, \
                         msi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_14 * lsi_826[k]
                    + f_6 * msh0_623[k]
                    - f_7 * msh1_623[k]
                    + f_3 * pc_x[k] * msi_826[k];

        t_1059[k] = pa_z[k] * lsk0_771[k]
                    - f_12 * pc_z[k] * lsk1_771[k];

        t_1060[k] = f_13 * lsi_598[k]
                    + f_3 * pc_z[k] * msi_822[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pa_z, pc_y, pc_z, lsk0_773, lsk0_774, \
                         lsi_599, lsi_600, lsi_630, lsk1_773, lsk1_774, \
                         msi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pa_z[k] * lsk0_773[k]
                    + f_14 * lsi_599[k]
                    - f_12 * pc_z[k] * lsk1_773[k];

        t_1062[k] = pa_z[k] * lsk0_774[k]
                    + f_15 * lsi_600[k]
                    - f_12 * pc_z[k] * lsk1_774[k];

        t_1063[k] = f_22 * lsi_630[k]
                    + f_3 * pc_y[k] * msi_826[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, pc_x, lsi_832, lsi_833, lsi_834, \
                         lsi_835, msh0_629, msh1_629, msi_832, msi_833, msi_834, \
                         msi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_14 * lsi_832[k]
                    + f_4 * msh0_629[k]
                    - f_5 * msh1_629[k]
                    + f_3 * pc_x[k] * msi_832[k];

        t_1065[k] = f_14 * lsi_833[k]
                    + f_3 * pc_x[k] * msi_833[k];

        t_1066[k] = f_14 * lsi_834[k]
                    + f_3 * pc_x[k] * msi_834[k];

        t_1067[k] = f_14 * lsi_835[k]
                    + f_3 * pc_x[k] * msi_835[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, pc_x, lsi_836, lsi_837, lsi_838, \
                         lsi_839, msi_836, msi_837, msi_838, msi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_14 * lsi_836[k]
                    + f_3 * pc_x[k] * msi_836[k];

        t_1069[k] = f_14 * lsi_837[k]
                    + f_3 * pc_x[k] * msi_837[k];

        t_1070[k] = f_14 * lsi_838[k]
                    + f_3 * pc_x[k] * msi_838[k];

        t_1071[k] = f_14 * lsi_839[k]
                    + f_3 * pc_x[k] * msi_839[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, pa_z, pc_y, pc_z, lsk0_784, lsi_609, lsi_639, \
                         lsk1_784, msh0_626, msh1_626, msi_833, \
                         msi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = pa_z[k] * lsk0_784[k]
                    - f_12 * pc_z[k] * lsk1_784[k];

        t_1073[k] = f_13 * lsi_609[k]
                    + f_3 * pc_z[k] * msi_833[k];

        t_1074[k] = f_22 * lsi_639[k]
                    + f_10 * msh0_626[k]
                    - f_11 * msh1_626[k]
                    + f_3 * pc_y[k] * msi_835[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, pc_y, lsi_640, lsi_641, lsi_642, msh0_627, \
                         msh0_628, msh0_629, msh1_627, msh1_628, msh1_629, msi_836, msi_837, \
                         msi_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_22 * lsi_640[k]
                    + f_8 * msh0_627[k]
                    - f_9 * msh1_627[k]
                    + f_3 * pc_y[k] * msi_836[k];

        t_1076[k] = f_22 * lsi_641[k]
                    + f_6 * msh0_628[k]
                    - f_7 * msh1_628[k]
                    + f_3 * pc_y[k] * msi_837[k];

        t_1077[k] = f_22 * lsi_642[k]
                    + f_4 * msh0_629[k]
                    - f_5 * msh1_629[k]
                    + f_3 * pc_y[k] * msi_838[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, pc_x, pc_y, pc_z, lsi_615, lsi_643, lsi_840, \
                         msh0_629, msh0_630, msh1_629, msh1_630, msi_839, \
                         msi_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_22 * lsi_643[k]
                    + f_3 * pc_y[k] * msi_839[k];

        t_1079[k] = f_13 * lsi_615[k]
                    + f_1 * msh0_629[k]
                    - f_2 * msh1_629[k]
                    + f_3 * pc_z[k] * msi_839[k];

        t_1080[k] = f_14 * lsi_840[k]
                    + f_1 * msh0_630[k]
                    - f_2 * msh1_630[k]
                    + f_3 * pc_x[k] * msi_840[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, lsi_616, lsi_644, \
                         lsi_646, lsi_843, msh0_633, msh1_633, msi_840, msi_842, \
                         msi_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_17 * lsi_644[k]
                    + f_3 * pc_y[k] * msi_840[k];

        t_1082[k] = f_14 * lsi_616[k]
                    + f_3 * pc_z[k] * msi_840[k];

        t_1083[k] = f_14 * lsi_843[k]
                    + f_10 * msh0_633[k]
                    - f_11 * msh1_633[k]
                    + f_3 * pc_x[k] * msi_843[k];

        t_1084[k] = f_17 * lsi_646[k]
                    + f_3 * pc_y[k] * msi_842[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, lsi_619, lsi_845, lsi_846, \
                         msh0_635, msh0_636, msh1_635, msh1_636, msi_843, msi_845, \
                         msi_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_14 * lsi_845[k]
                    + f_10 * msh0_635[k]
                    - f_11 * msh1_635[k]
                    + f_3 * pc_x[k] * msi_845[k];

        t_1086[k] = f_14 * lsi_846[k]
                    + f_8 * msh0_636[k]
                    - f_9 * msh1_636[k]
                    + f_3 * pc_x[k] * msi_846[k];

        t_1087[k] = f_14 * lsi_619[k]
                    + f_3 * pc_z[k] * msi_843[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, lsi_649, lsi_849, lsi_850, \
                         msh0_639, msh0_640, msh1_639, msh1_640, msi_845, msi_849, \
                         msi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * lsi_649[k]
                    + f_3 * pc_y[k] * msi_845[k];

        t_1089[k] = f_14 * lsi_849[k]
                    + f_8 * msh0_639[k]
                    - f_9 * msh1_639[k]
                    + f_3 * pc_x[k] * msi_849[k];

        t_1090[k] = f_14 * lsi_850[k]
                    + f_6 * msh0_640[k]
                    - f_7 * msh1_640[k]
                    + f_3 * pc_x[k] * msi_850[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, lsi_622, lsi_653, lsi_852, \
                         msh0_642, msh1_642, msi_846, msi_849, \
                         msi_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_14 * lsi_622[k]
                    + f_3 * pc_z[k] * msi_846[k];

        t_1092[k] = f_14 * lsi_852[k]
                    + f_6 * msh0_642[k]
                    - f_7 * msh1_642[k]
                    + f_3 * pc_x[k] * msi_852[k];

        t_1093[k] = f_17 * lsi_653[k]
                    + f_3 * pc_y[k] * msi_849[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, lsi_626, lsi_854, lsi_855, \
                         msh0_644, msh0_645, msh1_644, msh1_645, msi_850, msi_854, \
                         msi_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_14 * lsi_854[k]
                    + f_6 * msh0_644[k]
                    - f_7 * msh1_644[k]
                    + f_3 * pc_x[k] * msi_854[k];

        t_1095[k] = f_14 * lsi_855[k]
                    + f_4 * msh0_645[k]
                    - f_5 * msh1_645[k]
                    + f_3 * pc_x[k] * msi_855[k];

        t_1096[k] = f_14 * lsi_626[k]
                    + f_3 * pc_z[k] * msi_850[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, lsi_658, lsi_857, lsi_858, \
                         msh0_647, msh0_648, msh1_647, msh1_648, msi_854, msi_857, \
                         msi_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_14 * lsi_857[k]
                    + f_4 * msh0_647[k]
                    - f_5 * msh1_647[k]
                    + f_3 * pc_x[k] * msi_857[k];

        t_1098[k] = f_14 * lsi_858[k]
                    + f_4 * msh0_648[k]
                    - f_5 * msh1_648[k]
                    + f_3 * pc_x[k] * msi_858[k];

        t_1099[k] = f_17 * lsi_658[k]
                    + f_3 * pc_y[k] * msi_854[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, lsi_860, lsi_861, lsi_862, \
                         lsi_863, msh0_650, msh1_650, msi_860, msi_861, msi_862, \
                         msi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_14 * lsi_860[k]
                    + f_4 * msh0_650[k]
                    - f_5 * msh1_650[k]
                    + f_3 * pc_x[k] * msi_860[k];

        t_1101[k] = f_14 * lsi_861[k]
                    + f_3 * pc_x[k] * msi_861[k];

        t_1102[k] = f_14 * lsi_862[k]
                    + f_3 * pc_x[k] * msi_862[k];

        t_1103[k] = f_14 * lsi_863[k]
                    + f_3 * pc_x[k] * msi_863[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, lsi_864, lsi_865, lsi_866, \
                         lsi_867, msi_864, msi_865, msi_866, msi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_14 * lsi_864[k]
                    + f_3 * pc_x[k] * msi_864[k];

        t_1105[k] = f_14 * lsi_865[k]
                    + f_3 * pc_x[k] * msi_865[k];

        t_1106[k] = f_14 * lsi_866[k]
                    + f_3 * pc_x[k] * msi_866[k];

        t_1107[k] = f_14 * lsi_867[k]
                    + f_3 * pc_x[k] * msi_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, lsi_637, lsi_665, lsi_667, \
                         msh0_645, msh0_647, msh1_645, msh1_647, msi_861, \
                         msi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * lsi_665[k]
                    + f_1 * msh0_645[k]
                    - f_2 * msh1_645[k]
                    + f_3 * pc_y[k] * msi_861[k];

        t_1109[k] = f_14 * lsi_637[k]
                    + f_3 * pc_z[k] * msi_861[k];

        t_1110[k] = f_17 * lsi_667[k]
                    + f_10 * msh0_647[k]
                    - f_11 * msh1_647[k]
                    + f_3 * pc_y[k] * msi_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pc_y, lsi_668, lsi_669, lsi_670, msh0_648, \
                         msh0_649, msh0_650, msh1_648, msh1_649, msh1_650, msi_864, msi_865, \
                         msi_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_17 * lsi_668[k]
                    + f_8 * msh0_648[k]
                    - f_9 * msh1_648[k]
                    + f_3 * pc_y[k] * msi_864[k];

        t_1112[k] = f_17 * lsi_669[k]
                    + f_6 * msh0_649[k]
                    - f_7 * msh1_649[k]
                    + f_3 * pc_y[k] * msi_865[k];

        t_1113[k] = f_17 * lsi_670[k]
                    + f_4 * msh0_650[k]
                    - f_5 * msh1_650[k]
                    + f_3 * pc_y[k] * msi_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_x, pc_y, pc_z, lsi_643, lsi_671, lsi_868, \
                         msh0_650, msh0_651, msh1_650, msh1_651, msi_867, \
                         msi_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * lsi_671[k]
                    + f_3 * pc_y[k] * msi_867[k];

        t_1115[k] = f_14 * lsi_643[k]
                    + f_1 * msh0_650[k]
                    - f_2 * msh1_650[k]
                    + f_3 * pc_z[k] * msi_867[k];

        t_1116[k] = f_14 * lsi_868[k]
                    + f_1 * msh0_651[k]
                    - f_2 * msh1_651[k]
                    + f_3 * pc_x[k] * msi_868[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, pc_x, pc_y, pc_z, lsi_644, lsi_672, \
                         lsi_674, lsi_871, msh0_654, msh1_654, msi_868, msi_870, \
                         msi_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_16 * lsi_672[k]
                    + f_3 * pc_y[k] * msi_868[k];

        t_1118[k] = f_15 * lsi_644[k]
                    + f_3 * pc_z[k] * msi_868[k];

        t_1119[k] = f_14 * lsi_871[k]
                    + f_10 * msh0_654[k]
                    - f_11 * msh1_654[k]
                    + f_3 * pc_x[k] * msi_871[k];

        t_1120[k] = f_16 * lsi_674[k]
                    + f_3 * pc_y[k] * msi_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_x, pc_z, lsi_647, lsi_873, lsi_874, \
                         msh0_656, msh0_657, msh1_656, msh1_657, msi_871, msi_873, \
                         msi_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_14 * lsi_873[k]
                    + f_10 * msh0_656[k]
                    - f_11 * msh1_656[k]
                    + f_3 * pc_x[k] * msi_873[k];

        t_1122[k] = f_14 * lsi_874[k]
                    + f_8 * msh0_657[k]
                    - f_9 * msh1_657[k]
                    + f_3 * pc_x[k] * msi_874[k];

        t_1123[k] = f_15 * lsi_647[k]
                    + f_3 * pc_z[k] * msi_871[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, lsi_677, lsi_877, lsi_878, \
                         msh0_660, msh0_661, msh1_660, msh1_661, msi_873, msi_877, \
                         msi_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_16 * lsi_677[k]
                    + f_3 * pc_y[k] * msi_873[k];

        t_1125[k] = f_14 * lsi_877[k]
                    + f_8 * msh0_660[k]
                    - f_9 * msh1_660[k]
                    + f_3 * pc_x[k] * msi_877[k];

        t_1126[k] = f_14 * lsi_878[k]
                    + f_6 * msh0_661[k]
                    - f_7 * msh1_661[k]
                    + f_3 * pc_x[k] * msi_878[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, lsi_650, lsi_681, lsi_880, \
                         msh0_663, msh1_663, msi_874, msi_877, \
                         msi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_15 * lsi_650[k]
                    + f_3 * pc_z[k] * msi_874[k];

        t_1128[k] = f_14 * lsi_880[k]
                    + f_6 * msh0_663[k]
                    - f_7 * msh1_663[k]
                    + f_3 * pc_x[k] * msi_880[k];

        t_1129[k] = f_16 * lsi_681[k]
                    + f_3 * pc_y[k] * msi_877[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, lsi_654, lsi_882, lsi_883, \
                         msh0_665, msh0_666, msh1_665, msh1_666, msi_878, msi_882, \
                         msi_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_14 * lsi_882[k]
                    + f_6 * msh0_665[k]
                    - f_7 * msh1_665[k]
                    + f_3 * pc_x[k] * msi_882[k];

        t_1131[k] = f_14 * lsi_883[k]
                    + f_4 * msh0_666[k]
                    - f_5 * msh1_666[k]
                    + f_3 * pc_x[k] * msi_883[k];

        t_1132[k] = f_15 * lsi_654[k]
                    + f_3 * pc_z[k] * msi_878[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, lsi_686, lsi_885, lsi_886, \
                         msh0_668, msh0_669, msh1_668, msh1_669, msi_882, msi_885, \
                         msi_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_14 * lsi_885[k]
                    + f_4 * msh0_668[k]
                    - f_5 * msh1_668[k]
                    + f_3 * pc_x[k] * msi_885[k];

        t_1134[k] = f_14 * lsi_886[k]
                    + f_4 * msh0_669[k]
                    - f_5 * msh1_669[k]
                    + f_3 * pc_x[k] * msi_886[k];

        t_1135[k] = f_16 * lsi_686[k]
                    + f_3 * pc_y[k] * msi_882[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, pc_x, lsi_888, lsi_889, lsi_890, \
                         lsi_891, msh0_671, msh1_671, msi_888, msi_889, msi_890, \
                         msi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_14 * lsi_888[k]
                    + f_4 * msh0_671[k]
                    - f_5 * msh1_671[k]
                    + f_3 * pc_x[k] * msi_888[k];

        t_1137[k] = f_14 * lsi_889[k]
                    + f_3 * pc_x[k] * msi_889[k];

        t_1138[k] = f_14 * lsi_890[k]
                    + f_3 * pc_x[k] * msi_890[k];

        t_1139[k] = f_14 * lsi_891[k]
                    + f_3 * pc_x[k] * msi_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pc_x, lsi_892, lsi_893, lsi_894, \
                         lsi_895, msi_892, msi_893, msi_894, msi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_14 * lsi_892[k]
                    + f_3 * pc_x[k] * msi_892[k];

        t_1141[k] = f_14 * lsi_893[k]
                    + f_3 * pc_x[k] * msi_893[k];

        t_1142[k] = f_14 * lsi_894[k]
                    + f_3 * pc_x[k] * msi_894[k];

        t_1143[k] = f_14 * lsi_895[k]
                    + f_3 * pc_x[k] * msi_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pc_y, pc_z, lsi_665, lsi_693, lsi_695, \
                         msh0_666, msh0_668, msh1_666, msh1_668, msi_889, \
                         msi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_16 * lsi_693[k]
                    + f_1 * msh0_666[k]
                    - f_2 * msh1_666[k]
                    + f_3 * pc_y[k] * msi_889[k];

        t_1145[k] = f_15 * lsi_665[k]
                    + f_3 * pc_z[k] * msi_889[k];

        t_1146[k] = f_16 * lsi_695[k]
                    + f_10 * msh0_668[k]
                    - f_11 * msh1_668[k]
                    + f_3 * pc_y[k] * msi_891[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pc_y, lsi_696, lsi_697, lsi_698, msh0_669, \
                         msh0_670, msh0_671, msh1_669, msh1_670, msh1_671, msi_892, msi_893, \
                         msi_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_16 * lsi_696[k]
                    + f_8 * msh0_669[k]
                    - f_9 * msh1_669[k]
                    + f_3 * pc_y[k] * msi_892[k];

        t_1148[k] = f_16 * lsi_697[k]
                    + f_6 * msh0_670[k]
                    - f_7 * msh1_670[k]
                    + f_3 * pc_y[k] * msi_893[k];

        t_1149[k] = f_16 * lsi_698[k]
                    + f_4 * msh0_671[k]
                    - f_5 * msh1_671[k]
                    + f_3 * pc_y[k] * msi_894[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, lsi_671, lsi_699, lsi_896, \
                         msh0_671, msh0_672, msh1_671, msh1_672, msi_895, \
                         msi_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_16 * lsi_699[k]
                    + f_3 * pc_y[k] * msi_895[k];

        t_1151[k] = f_15 * lsi_671[k]
                    + f_1 * msh0_671[k]
                    - f_2 * msh1_671[k]
                    + f_3 * pc_z[k] * msi_895[k];

        t_1152[k] = f_14 * lsi_896[k]
                    + f_1 * msh0_672[k]
                    - f_2 * msh1_672[k]
                    + f_3 * pc_x[k] * msi_896[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsk0,
                                                           const size_t lsi, const size_t lsk1,
                                                           const size_t msh0, const size_t msh1,
                                                           const size_t msi, const size_t ncols,
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
    const auto f_22 = 3.0 / q;

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

    const auto *lsk0_972 = buffer.data(lsk0 + 972);
    const auto *lsk0_975 = buffer.data(lsk0 + 975);
    const auto *lsk0_977 = buffer.data(lsk0 + 977);
    const auto *lsk0_978 = buffer.data(lsk0 + 978);
    const auto *lsk0_981 = buffer.data(lsk0 + 981);
    const auto *lsk0_982 = buffer.data(lsk0 + 982);
    const auto *lsk0_984 = buffer.data(lsk0 + 984);
    const auto *lsk0_986 = buffer.data(lsk0 + 986);
    const auto *lsk0_987 = buffer.data(lsk0 + 987);
    const auto *lsk0_989 = buffer.data(lsk0 + 989);
    const auto *lsk0_990 = buffer.data(lsk0 + 990);
    const auto *lsk0_992 = buffer.data(lsk0 + 992);
    const auto *lsk0_1007 = buffer.data(lsk0 + 1007);

    const auto *lsi_672 = buffer.data(lsi + 672);
    const auto *lsi_675 = buffer.data(lsi + 675);
    const auto *lsi_678 = buffer.data(lsi + 678);
    const auto *lsi_682 = buffer.data(lsi + 682);
    const auto *lsi_693 = buffer.data(lsi + 693);
    const auto *lsi_699 = buffer.data(lsi + 699);
    const auto *lsi_700 = buffer.data(lsi + 700);
    const auto *lsi_702 = buffer.data(lsi + 702);
    const auto *lsi_703 = buffer.data(lsi + 703);
    const auto *lsi_705 = buffer.data(lsi + 705);
    const auto *lsi_706 = buffer.data(lsi + 706);
    const auto *lsi_709 = buffer.data(lsi + 709);
    const auto *lsi_710 = buffer.data(lsi + 710);
    const auto *lsi_714 = buffer.data(lsi + 714);
    const auto *lsi_721 = buffer.data(lsi + 721);
    const auto *lsi_723 = buffer.data(lsi + 723);
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
    const auto *lsi_751 = buffer.data(lsi + 751);
    const auto *lsi_752 = buffer.data(lsi + 752);
    const auto *lsi_753 = buffer.data(lsi + 753);
    const auto *lsi_754 = buffer.data(lsi + 754);
    const auto *lsi_755 = buffer.data(lsi + 755);
    const auto *lsi_756 = buffer.data(lsi + 756);
    const auto *lsi_757 = buffer.data(lsi + 757);
    const auto *lsi_758 = buffer.data(lsi + 758);
    const auto *lsi_759 = buffer.data(lsi + 759);
    const auto *lsi_761 = buffer.data(lsi + 761);
    const auto *lsi_762 = buffer.data(lsi + 762);
    const auto *lsi_764 = buffer.data(lsi + 764);
    const auto *lsi_765 = buffer.data(lsi + 765);
    const auto *lsi_766 = buffer.data(lsi + 766);
    const auto *lsi_768 = buffer.data(lsi + 768);
    const auto *lsi_769 = buffer.data(lsi + 769);
    const auto *lsi_770 = buffer.data(lsi + 770);
    const auto *lsi_777 = buffer.data(lsi + 777);
    const auto *lsi_779 = buffer.data(lsi + 779);
    const auto *lsi_780 = buffer.data(lsi + 780);
    const auto *lsi_781 = buffer.data(lsi + 781);
    const auto *lsi_782 = buffer.data(lsi + 782);
    const auto *lsi_783 = buffer.data(lsi + 783);
    const auto *lsi_899 = buffer.data(lsi + 899);
    const auto *lsi_901 = buffer.data(lsi + 901);
    const auto *lsi_902 = buffer.data(lsi + 902);
    const auto *lsi_905 = buffer.data(lsi + 905);
    const auto *lsi_906 = buffer.data(lsi + 906);
    const auto *lsi_908 = buffer.data(lsi + 908);
    const auto *lsi_910 = buffer.data(lsi + 910);
    const auto *lsi_911 = buffer.data(lsi + 911);
    const auto *lsi_913 = buffer.data(lsi + 913);
    const auto *lsi_914 = buffer.data(lsi + 914);
    const auto *lsi_916 = buffer.data(lsi + 916);
    const auto *lsi_917 = buffer.data(lsi + 917);
    const auto *lsi_918 = buffer.data(lsi + 918);
    const auto *lsi_919 = buffer.data(lsi + 919);
    const auto *lsi_920 = buffer.data(lsi + 920);
    const auto *lsi_921 = buffer.data(lsi + 921);
    const auto *lsi_922 = buffer.data(lsi + 922);
    const auto *lsi_923 = buffer.data(lsi + 923);
    const auto *lsi_924 = buffer.data(lsi + 924);
    const auto *lsi_927 = buffer.data(lsi + 927);
    const auto *lsi_929 = buffer.data(lsi + 929);
    const auto *lsi_930 = buffer.data(lsi + 930);
    const auto *lsi_933 = buffer.data(lsi + 933);
    const auto *lsi_934 = buffer.data(lsi + 934);
    const auto *lsi_936 = buffer.data(lsi + 936);
    const auto *lsi_938 = buffer.data(lsi + 938);
    const auto *lsi_939 = buffer.data(lsi + 939);
    const auto *lsi_941 = buffer.data(lsi + 941);
    const auto *lsi_942 = buffer.data(lsi + 942);
    const auto *lsi_944 = buffer.data(lsi + 944);
    const auto *lsi_945 = buffer.data(lsi + 945);
    const auto *lsi_946 = buffer.data(lsi + 946);
    const auto *lsi_947 = buffer.data(lsi + 947);
    const auto *lsi_948 = buffer.data(lsi + 948);
    const auto *lsi_949 = buffer.data(lsi + 949);
    const auto *lsi_950 = buffer.data(lsi + 950);
    const auto *lsi_951 = buffer.data(lsi + 951);
    const auto *lsi_973 = buffer.data(lsi + 973);
    const auto *lsi_974 = buffer.data(lsi + 974);
    const auto *lsi_975 = buffer.data(lsi + 975);
    const auto *lsi_976 = buffer.data(lsi + 976);
    const auto *lsi_977 = buffer.data(lsi + 977);
    const auto *lsi_978 = buffer.data(lsi + 978);
    const auto *lsi_979 = buffer.data(lsi + 979);

    const auto *lsk1_972 = buffer.data(lsk1 + 972);
    const auto *lsk1_975 = buffer.data(lsk1 + 975);
    const auto *lsk1_977 = buffer.data(lsk1 + 977);
    const auto *lsk1_978 = buffer.data(lsk1 + 978);
    const auto *lsk1_981 = buffer.data(lsk1 + 981);
    const auto *lsk1_982 = buffer.data(lsk1 + 982);
    const auto *lsk1_984 = buffer.data(lsk1 + 984);
    const auto *lsk1_986 = buffer.data(lsk1 + 986);
    const auto *lsk1_987 = buffer.data(lsk1 + 987);
    const auto *lsk1_989 = buffer.data(lsk1 + 989);
    const auto *lsk1_990 = buffer.data(lsk1 + 990);
    const auto *lsk1_992 = buffer.data(lsk1 + 992);
    const auto *lsk1_1007 = buffer.data(lsk1 + 1007);

    const auto *msh0_675 = buffer.data(msh0 + 675);
    const auto *msh0_677 = buffer.data(msh0 + 677);
    const auto *msh0_678 = buffer.data(msh0 + 678);
    const auto *msh0_681 = buffer.data(msh0 + 681);
    const auto *msh0_682 = buffer.data(msh0 + 682);
    const auto *msh0_684 = buffer.data(msh0 + 684);
    const auto *msh0_686 = buffer.data(msh0 + 686);
    const auto *msh0_687 = buffer.data(msh0 + 687);
    const auto *msh0_689 = buffer.data(msh0 + 689);
    const auto *msh0_690 = buffer.data(msh0 + 690);
    const auto *msh0_691 = buffer.data(msh0 + 691);
    const auto *msh0_692 = buffer.data(msh0 + 692);
    const auto *msh0_693 = buffer.data(msh0 + 693);
    const auto *msh0_696 = buffer.data(msh0 + 696);
    const auto *msh0_698 = buffer.data(msh0 + 698);
    const auto *msh0_699 = buffer.data(msh0 + 699);
    const auto *msh0_702 = buffer.data(msh0 + 702);
    const auto *msh0_703 = buffer.data(msh0 + 703);
    const auto *msh0_705 = buffer.data(msh0 + 705);
    const auto *msh0_707 = buffer.data(msh0 + 707);
    const auto *msh0_708 = buffer.data(msh0 + 708);
    const auto *msh0_710 = buffer.data(msh0 + 710);
    const auto *msh0_711 = buffer.data(msh0 + 711);
    const auto *msh0_712 = buffer.data(msh0 + 712);
    const auto *msh0_713 = buffer.data(msh0 + 713);
    const auto *msh0_729 = buffer.data(msh0 + 729);
    const auto *msh0_731 = buffer.data(msh0 + 731);
    const auto *msh0_732 = buffer.data(msh0 + 732);
    const auto *msh0_733 = buffer.data(msh0 + 733);
    const auto *msh0_734 = buffer.data(msh0 + 734);

    const auto *msh1_675 = buffer.data(msh1 + 675);
    const auto *msh1_677 = buffer.data(msh1 + 677);
    const auto *msh1_678 = buffer.data(msh1 + 678);
    const auto *msh1_681 = buffer.data(msh1 + 681);
    const auto *msh1_682 = buffer.data(msh1 + 682);
    const auto *msh1_684 = buffer.data(msh1 + 684);
    const auto *msh1_686 = buffer.data(msh1 + 686);
    const auto *msh1_687 = buffer.data(msh1 + 687);
    const auto *msh1_689 = buffer.data(msh1 + 689);
    const auto *msh1_690 = buffer.data(msh1 + 690);
    const auto *msh1_691 = buffer.data(msh1 + 691);
    const auto *msh1_692 = buffer.data(msh1 + 692);
    const auto *msh1_693 = buffer.data(msh1 + 693);
    const auto *msh1_696 = buffer.data(msh1 + 696);
    const auto *msh1_698 = buffer.data(msh1 + 698);
    const auto *msh1_699 = buffer.data(msh1 + 699);
    const auto *msh1_702 = buffer.data(msh1 + 702);
    const auto *msh1_703 = buffer.data(msh1 + 703);
    const auto *msh1_705 = buffer.data(msh1 + 705);
    const auto *msh1_707 = buffer.data(msh1 + 707);
    const auto *msh1_708 = buffer.data(msh1 + 708);
    const auto *msh1_710 = buffer.data(msh1 + 710);
    const auto *msh1_711 = buffer.data(msh1 + 711);
    const auto *msh1_712 = buffer.data(msh1 + 712);
    const auto *msh1_713 = buffer.data(msh1 + 713);
    const auto *msh1_729 = buffer.data(msh1 + 729);
    const auto *msh1_731 = buffer.data(msh1 + 731);
    const auto *msh1_732 = buffer.data(msh1 + 732);
    const auto *msh1_733 = buffer.data(msh1 + 733);
    const auto *msh1_734 = buffer.data(msh1 + 734);

    const auto *msi_896 = buffer.data(msi + 896);
    const auto *msi_898 = buffer.data(msi + 898);
    const auto *msi_899 = buffer.data(msi + 899);
    const auto *msi_901 = buffer.data(msi + 901);
    const auto *msi_902 = buffer.data(msi + 902);
    const auto *msi_905 = buffer.data(msi + 905);
    const auto *msi_906 = buffer.data(msi + 906);
    const auto *msi_908 = buffer.data(msi + 908);
    const auto *msi_910 = buffer.data(msi + 910);
    const auto *msi_911 = buffer.data(msi + 911);
    const auto *msi_913 = buffer.data(msi + 913);
    const auto *msi_914 = buffer.data(msi + 914);
    const auto *msi_916 = buffer.data(msi + 916);
    const auto *msi_917 = buffer.data(msi + 917);
    const auto *msi_918 = buffer.data(msi + 918);
    const auto *msi_919 = buffer.data(msi + 919);
    const auto *msi_920 = buffer.data(msi + 920);
    const auto *msi_921 = buffer.data(msi + 921);
    const auto *msi_922 = buffer.data(msi + 922);
    const auto *msi_923 = buffer.data(msi + 923);
    const auto *msi_924 = buffer.data(msi + 924);
    const auto *msi_926 = buffer.data(msi + 926);
    const auto *msi_927 = buffer.data(msi + 927);
    const auto *msi_929 = buffer.data(msi + 929);
    const auto *msi_930 = buffer.data(msi + 930);
    const auto *msi_933 = buffer.data(msi + 933);
    const auto *msi_934 = buffer.data(msi + 934);
    const auto *msi_936 = buffer.data(msi + 936);
    const auto *msi_938 = buffer.data(msi + 938);
    const auto *msi_939 = buffer.data(msi + 939);
    const auto *msi_941 = buffer.data(msi + 941);
    const auto *msi_942 = buffer.data(msi + 942);
    const auto *msi_944 = buffer.data(msi + 944);
    const auto *msi_945 = buffer.data(msi + 945);
    const auto *msi_946 = buffer.data(msi + 946);
    const auto *msi_947 = buffer.data(msi + 947);
    const auto *msi_948 = buffer.data(msi + 948);
    const auto *msi_949 = buffer.data(msi + 949);
    const auto *msi_950 = buffer.data(msi + 950);
    const auto *msi_951 = buffer.data(msi + 951);
    const auto *msi_952 = buffer.data(msi + 952);
    const auto *msi_954 = buffer.data(msi + 954);
    const auto *msi_955 = buffer.data(msi + 955);
    const auto *msi_957 = buffer.data(msi + 957);
    const auto *msi_958 = buffer.data(msi + 958);
    const auto *msi_961 = buffer.data(msi + 961);
    const auto *msi_962 = buffer.data(msi + 962);
    const auto *msi_966 = buffer.data(msi + 966);
    const auto *msi_973 = buffer.data(msi + 973);
    const auto *msi_974 = buffer.data(msi + 974);
    const auto *msi_975 = buffer.data(msi + 975);
    const auto *msi_976 = buffer.data(msi + 976);
    const auto *msi_977 = buffer.data(msi + 977);
    const auto *msi_978 = buffer.data(msi + 978);
    const auto *msi_979 = buffer.data(msi + 979);

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pc_x, pc_y, pc_z, lsi_672, lsi_700, \
                         lsi_702, lsi_899, msh0_675, msh1_675, msi_896, msi_898, \
                         msi_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_15 * lsi_700[k]
                    + f_3 * pc_y[k] * msi_896[k];

        t_1154[k] = f_16 * lsi_672[k]
                    + f_3 * pc_z[k] * msi_896[k];

        t_1155[k] = f_14 * lsi_899[k]
                    + f_10 * msh0_675[k]
                    - f_11 * msh1_675[k]
                    + f_3 * pc_x[k] * msi_899[k];

        t_1156[k] = f_15 * lsi_702[k]
                    + f_3 * pc_y[k] * msi_898[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, pc_z, lsi_675, lsi_901, lsi_902, \
                         msh0_677, msh0_678, msh1_677, msh1_678, msi_899, msi_901, \
                         msi_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_14 * lsi_901[k]
                    + f_10 * msh0_677[k]
                    - f_11 * msh1_677[k]
                    + f_3 * pc_x[k] * msi_901[k];

        t_1158[k] = f_14 * lsi_902[k]
                    + f_8 * msh0_678[k]
                    - f_9 * msh1_678[k]
                    + f_3 * pc_x[k] * msi_902[k];

        t_1159[k] = f_16 * lsi_675[k]
                    + f_3 * pc_z[k] * msi_899[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, lsi_705, lsi_905, lsi_906, \
                         msh0_681, msh0_682, msh1_681, msh1_682, msi_901, msi_905, \
                         msi_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * lsi_705[k]
                    + f_3 * pc_y[k] * msi_901[k];

        t_1161[k] = f_14 * lsi_905[k]
                    + f_8 * msh0_681[k]
                    - f_9 * msh1_681[k]
                    + f_3 * pc_x[k] * msi_905[k];

        t_1162[k] = f_14 * lsi_906[k]
                    + f_6 * msh0_682[k]
                    - f_7 * msh1_682[k]
                    + f_3 * pc_x[k] * msi_906[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_x, pc_y, pc_z, lsi_678, lsi_709, lsi_908, \
                         msh0_684, msh1_684, msi_902, msi_905, \
                         msi_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * lsi_678[k]
                    + f_3 * pc_z[k] * msi_902[k];

        t_1164[k] = f_14 * lsi_908[k]
                    + f_6 * msh0_684[k]
                    - f_7 * msh1_684[k]
                    + f_3 * pc_x[k] * msi_908[k];

        t_1165[k] = f_15 * lsi_709[k]
                    + f_3 * pc_y[k] * msi_905[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_x, pc_z, lsi_682, lsi_910, lsi_911, \
                         msh0_686, msh0_687, msh1_686, msh1_687, msi_906, msi_910, \
                         msi_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_14 * lsi_910[k]
                    + f_6 * msh0_686[k]
                    - f_7 * msh1_686[k]
                    + f_3 * pc_x[k] * msi_910[k];

        t_1167[k] = f_14 * lsi_911[k]
                    + f_4 * msh0_687[k]
                    - f_5 * msh1_687[k]
                    + f_3 * pc_x[k] * msi_911[k];

        t_1168[k] = f_16 * lsi_682[k]
                    + f_3 * pc_z[k] * msi_906[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_x, pc_y, lsi_714, lsi_913, lsi_914, \
                         msh0_689, msh0_690, msh1_689, msh1_690, msi_910, msi_913, \
                         msi_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_14 * lsi_913[k]
                    + f_4 * msh0_689[k]
                    - f_5 * msh1_689[k]
                    + f_3 * pc_x[k] * msi_913[k];

        t_1170[k] = f_14 * lsi_914[k]
                    + f_4 * msh0_690[k]
                    - f_5 * msh1_690[k]
                    + f_3 * pc_x[k] * msi_914[k];

        t_1171[k] = f_15 * lsi_714[k]
                    + f_3 * pc_y[k] * msi_910[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, pc_x, lsi_916, lsi_917, lsi_918, \
                         lsi_919, msh0_692, msh1_692, msi_916, msi_917, msi_918, \
                         msi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_14 * lsi_916[k]
                    + f_4 * msh0_692[k]
                    - f_5 * msh1_692[k]
                    + f_3 * pc_x[k] * msi_916[k];

        t_1173[k] = f_14 * lsi_917[k]
                    + f_3 * pc_x[k] * msi_917[k];

        t_1174[k] = f_14 * lsi_918[k]
                    + f_3 * pc_x[k] * msi_918[k];

        t_1175[k] = f_14 * lsi_919[k]
                    + f_3 * pc_x[k] * msi_919[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pc_x, lsi_920, lsi_921, lsi_922, \
                         lsi_923, msi_920, msi_921, msi_922, msi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_14 * lsi_920[k]
                    + f_3 * pc_x[k] * msi_920[k];

        t_1177[k] = f_14 * lsi_921[k]
                    + f_3 * pc_x[k] * msi_921[k];

        t_1178[k] = f_14 * lsi_922[k]
                    + f_3 * pc_x[k] * msi_922[k];

        t_1179[k] = f_14 * lsi_923[k]
                    + f_3 * pc_x[k] * msi_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_y, pc_z, lsi_693, lsi_721, lsi_723, \
                         msh0_687, msh0_689, msh1_687, msh1_689, msi_917, \
                         msi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_15 * lsi_721[k]
                    + f_1 * msh0_687[k]
                    - f_2 * msh1_687[k]
                    + f_3 * pc_y[k] * msi_917[k];

        t_1181[k] = f_16 * lsi_693[k]
                    + f_3 * pc_z[k] * msi_917[k];

        t_1182[k] = f_15 * lsi_723[k]
                    + f_10 * msh0_689[k]
                    - f_11 * msh1_689[k]
                    + f_3 * pc_y[k] * msi_919[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pc_y, lsi_724, lsi_725, lsi_726, msh0_690, \
                         msh0_691, msh0_692, msh1_690, msh1_691, msh1_692, msi_920, msi_921, \
                         msi_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_15 * lsi_724[k]
                    + f_8 * msh0_690[k]
                    - f_9 * msh1_690[k]
                    + f_3 * pc_y[k] * msi_920[k];

        t_1184[k] = f_15 * lsi_725[k]
                    + f_6 * msh0_691[k]
                    - f_7 * msh1_691[k]
                    + f_3 * pc_y[k] * msi_921[k];

        t_1185[k] = f_15 * lsi_726[k]
                    + f_4 * msh0_692[k]
                    - f_5 * msh1_692[k]
                    + f_3 * pc_y[k] * msi_922[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pc_x, pc_y, pc_z, lsi_699, lsi_727, lsi_924, \
                         msh0_692, msh0_693, msh1_692, msh1_693, msi_923, \
                         msi_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_15 * lsi_727[k]
                    + f_3 * pc_y[k] * msi_923[k];

        t_1187[k] = f_16 * lsi_699[k]
                    + f_1 * msh0_692[k]
                    - f_2 * msh1_692[k]
                    + f_3 * pc_z[k] * msi_923[k];

        t_1188[k] = f_14 * lsi_924[k]
                    + f_1 * msh0_693[k]
                    - f_2 * msh1_693[k]
                    + f_3 * pc_x[k] * msi_924[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pc_x, pc_y, pc_z, lsi_700, lsi_728, \
                         lsi_730, lsi_927, msh0_696, msh1_696, msi_924, msi_926, \
                         msi_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_14 * lsi_728[k]
                    + f_3 * pc_y[k] * msi_924[k];

        t_1190[k] = f_17 * lsi_700[k]
                    + f_3 * pc_z[k] * msi_924[k];

        t_1191[k] = f_14 * lsi_927[k]
                    + f_10 * msh0_696[k]
                    - f_11 * msh1_696[k]
                    + f_3 * pc_x[k] * msi_927[k];

        t_1192[k] = f_14 * lsi_730[k]
                    + f_3 * pc_y[k] * msi_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pc_x, pc_z, lsi_703, lsi_929, lsi_930, \
                         msh0_698, msh0_699, msh1_698, msh1_699, msi_927, msi_929, \
                         msi_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_14 * lsi_929[k]
                    + f_10 * msh0_698[k]
                    - f_11 * msh1_698[k]
                    + f_3 * pc_x[k] * msi_929[k];

        t_1194[k] = f_14 * lsi_930[k]
                    + f_8 * msh0_699[k]
                    - f_9 * msh1_699[k]
                    + f_3 * pc_x[k] * msi_930[k];

        t_1195[k] = f_17 * lsi_703[k]
                    + f_3 * pc_z[k] * msi_927[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, lsi_733, lsi_933, lsi_934, \
                         msh0_702, msh0_703, msh1_702, msh1_703, msi_929, msi_933, \
                         msi_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_14 * lsi_733[k]
                    + f_3 * pc_y[k] * msi_929[k];

        t_1197[k] = f_14 * lsi_933[k]
                    + f_8 * msh0_702[k]
                    - f_9 * msh1_702[k]
                    + f_3 * pc_x[k] * msi_933[k];

        t_1198[k] = f_14 * lsi_934[k]
                    + f_6 * msh0_703[k]
                    - f_7 * msh1_703[k]
                    + f_3 * pc_x[k] * msi_934[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_y, pc_z, lsi_706, lsi_737, lsi_936, \
                         msh0_705, msh1_705, msi_930, msi_933, \
                         msi_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * lsi_706[k]
                    + f_3 * pc_z[k] * msi_930[k];

        t_1200[k] = f_14 * lsi_936[k]
                    + f_6 * msh0_705[k]
                    - f_7 * msh1_705[k]
                    + f_3 * pc_x[k] * msi_936[k];

        t_1201[k] = f_14 * lsi_737[k]
                    + f_3 * pc_y[k] * msi_933[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, pc_z, lsi_710, lsi_938, lsi_939, \
                         msh0_707, msh0_708, msh1_707, msh1_708, msi_934, msi_938, \
                         msi_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_14 * lsi_938[k]
                    + f_6 * msh0_707[k]
                    - f_7 * msh1_707[k]
                    + f_3 * pc_x[k] * msi_938[k];

        t_1203[k] = f_14 * lsi_939[k]
                    + f_4 * msh0_708[k]
                    - f_5 * msh1_708[k]
                    + f_3 * pc_x[k] * msi_939[k];

        t_1204[k] = f_17 * lsi_710[k]
                    + f_3 * pc_z[k] * msi_934[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pc_x, pc_y, lsi_742, lsi_941, lsi_942, \
                         msh0_710, msh0_711, msh1_710, msh1_711, msi_938, msi_941, \
                         msi_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_14 * lsi_941[k]
                    + f_4 * msh0_710[k]
                    - f_5 * msh1_710[k]
                    + f_3 * pc_x[k] * msi_941[k];

        t_1206[k] = f_14 * lsi_942[k]
                    + f_4 * msh0_711[k]
                    - f_5 * msh1_711[k]
                    + f_3 * pc_x[k] * msi_942[k];

        t_1207[k] = f_14 * lsi_742[k]
                    + f_3 * pc_y[k] * msi_938[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, pc_x, lsi_944, lsi_945, lsi_946, \
                         lsi_947, msh0_713, msh1_713, msi_944, msi_945, msi_946, \
                         msi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_14 * lsi_944[k]
                    + f_4 * msh0_713[k]
                    - f_5 * msh1_713[k]
                    + f_3 * pc_x[k] * msi_944[k];

        t_1209[k] = f_14 * lsi_945[k]
                    + f_3 * pc_x[k] * msi_945[k];

        t_1210[k] = f_14 * lsi_946[k]
                    + f_3 * pc_x[k] * msi_946[k];

        t_1211[k] = f_14 * lsi_947[k]
                    + f_3 * pc_x[k] * msi_947[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, pc_x, lsi_948, lsi_949, lsi_950, \
                         lsi_951, msi_948, msi_949, msi_950, msi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_14 * lsi_948[k]
                    + f_3 * pc_x[k] * msi_948[k];

        t_1213[k] = f_14 * lsi_949[k]
                    + f_3 * pc_x[k] * msi_949[k];

        t_1214[k] = f_14 * lsi_950[k]
                    + f_3 * pc_x[k] * msi_950[k];

        t_1215[k] = f_14 * lsi_951[k]
                    + f_3 * pc_x[k] * msi_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_y, pc_z, lsi_721, lsi_749, lsi_751, \
                         msh0_708, msh0_710, msh1_708, msh1_710, msi_945, \
                         msi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_14 * lsi_749[k]
                    + f_1 * msh0_708[k]
                    - f_2 * msh1_708[k]
                    + f_3 * pc_y[k] * msi_945[k];

        t_1217[k] = f_17 * lsi_721[k]
                    + f_3 * pc_z[k] * msi_945[k];

        t_1218[k] = f_14 * lsi_751[k]
                    + f_10 * msh0_710[k]
                    - f_11 * msh1_710[k]
                    + f_3 * pc_y[k] * msi_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pc_y, lsi_752, lsi_753, lsi_754, msh0_711, \
                         msh0_712, msh0_713, msh1_711, msh1_712, msh1_713, msi_948, msi_949, \
                         msi_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_14 * lsi_752[k]
                    + f_8 * msh0_711[k]
                    - f_9 * msh1_711[k]
                    + f_3 * pc_y[k] * msi_948[k];

        t_1220[k] = f_14 * lsi_753[k]
                    + f_6 * msh0_712[k]
                    - f_7 * msh1_712[k]
                    + f_3 * pc_y[k] * msi_949[k];

        t_1221[k] = f_14 * lsi_754[k]
                    + f_4 * msh0_713[k]
                    - f_5 * msh1_713[k]
                    + f_3 * pc_y[k] * msi_950[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pa_y, pc_y, pc_z, lsk0_972, lsi_727, \
                         lsi_755, lsi_756, lsk1_972, msh0_713, msh1_713, msi_951, \
                         msi_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * lsi_755[k]
                    + f_3 * pc_y[k] * msi_951[k];

        t_1223[k] = f_17 * lsi_727[k]
                    + f_1 * msh0_713[k]
                    - f_2 * msh1_713[k]
                    + f_3 * pc_z[k] * msi_951[k];

        t_1224[k] = pa_y[k] * lsk0_972[k]
                    - f_12 * pc_y[k] * lsk1_972[k];

        t_1225[k] = f_13 * lsi_756[k]
                    + f_3 * pc_y[k] * msi_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pa_y, pc_y, pc_z, lsk0_975, lsk0_977, \
                         lsi_728, lsi_757, lsi_758, lsk1_975, lsk1_977, msi_952, \
                         msi_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_22 * lsi_728[k]
                    + f_3 * pc_z[k] * msi_952[k];

        t_1227[k] = pa_y[k] * lsk0_975[k]
                    + f_14 * lsi_757[k]
                    - f_12 * pc_y[k] * lsk1_975[k];

        t_1228[k] = f_13 * lsi_758[k]
                    + f_3 * pc_y[k] * msi_954[k];

        t_1229[k] = pa_y[k] * lsk0_977[k]
                    - f_12 * pc_y[k] * lsk1_977[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pa_y, pc_y, pc_z, lsk0_978, lsk0_981, \
                         lsi_731, lsi_759, lsi_761, lsk1_978, lsk1_981, msi_955, \
                         msi_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = pa_y[k] * lsk0_978[k]
                    + f_15 * lsi_759[k]
                    - f_12 * pc_y[k] * lsk1_978[k];

        t_1231[k] = f_22 * lsi_731[k]
                    + f_3 * pc_z[k] * msi_955[k];

        t_1232[k] = f_13 * lsi_761[k]
                    + f_3 * pc_y[k] * msi_957[k];

        t_1233[k] = pa_y[k] * lsk0_981[k]
                    - f_12 * pc_y[k] * lsk1_981[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pa_y, pc_y, pc_z, lsk0_982, lsk0_984, \
                         lsi_734, lsi_762, lsi_764, lsk1_982, lsk1_984, \
                         msi_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = pa_y[k] * lsk0_982[k]
                    + f_16 * lsi_762[k]
                    - f_12 * pc_y[k] * lsk1_982[k];

        t_1235[k] = f_22 * lsi_734[k]
                    + f_3 * pc_z[k] * msi_958[k];

        t_1236[k] = pa_y[k] * lsk0_984[k]
                    + f_14 * lsi_764[k]
                    - f_12 * pc_y[k] * lsk1_984[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pa_y, pc_y, pc_z, lsk0_986, lsk0_987, \
                         lsi_738, lsi_765, lsi_766, lsk1_986, lsk1_987, msi_961, \
                         msi_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_13 * lsi_765[k]
                    + f_3 * pc_y[k] * msi_961[k];

        t_1238[k] = pa_y[k] * lsk0_986[k]
                    - f_12 * pc_y[k] * lsk1_986[k];

        t_1239[k] = pa_y[k] * lsk0_987[k]
                    + f_17 * lsi_766[k]
                    - f_12 * pc_y[k] * lsk1_987[k];

        t_1240[k] = f_22 * lsi_738[k]
                    + f_3 * pc_z[k] * msi_962[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, pa_y, pc_y, lsk0_989, lsk0_990, \
                         lsk0_992, lsi_768, lsi_769, lsi_770, lsk1_989, lsk1_990, lsk1_992, \
                         msi_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = pa_y[k] * lsk0_989[k]
                    + f_15 * lsi_768[k]
                    - f_12 * pc_y[k] * lsk1_989[k];

        t_1242[k] = pa_y[k] * lsk0_990[k]
                    + f_14 * lsi_769[k]
                    - f_12 * pc_y[k] * lsk1_990[k];

        t_1243[k] = f_13 * lsi_770[k]
                    + f_3 * pc_y[k] * msi_966[k];

        t_1244[k] = pa_y[k] * lsk0_992[k]
                    - f_12 * pc_y[k] * lsk1_992[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, lsi_973, lsi_974, \
                         lsi_975, lsi_976, lsi_977, msi_973, msi_974, msi_975, msi_976, \
                         msi_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_14 * lsi_973[k]
                    + f_3 * pc_x[k] * msi_973[k];

        t_1246[k] = f_14 * lsi_974[k]
                    + f_3 * pc_x[k] * msi_974[k];

        t_1247[k] = f_14 * lsi_975[k]
                    + f_3 * pc_x[k] * msi_975[k];

        t_1248[k] = f_14 * lsi_976[k]
                    + f_3 * pc_x[k] * msi_976[k];

        t_1249[k] = f_14 * lsi_977[k]
                    + f_3 * pc_x[k] * msi_977[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, pc_z, lsi_749, lsi_777, \
                         lsi_978, lsi_979, msh0_729, msh1_729, msi_973, msi_978, \
                         msi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_14 * lsi_978[k]
                    + f_3 * pc_x[k] * msi_978[k];

        t_1251[k] = f_14 * lsi_979[k]
                    + f_3 * pc_x[k] * msi_979[k];

        t_1252[k] = f_13 * lsi_777[k]
                    + f_1 * msh0_729[k]
                    - f_2 * msh1_729[k]
                    + f_3 * pc_y[k] * msi_973[k];

        t_1253[k] = f_22 * lsi_749[k]
                    + f_3 * pc_z[k] * msi_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, lsi_779, lsi_780, lsi_781, msh0_731, \
                         msh0_732, msh0_733, msh1_731, msh1_732, msh1_733, msi_975, msi_976, \
                         msi_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_13 * lsi_779[k]
                    + f_10 * msh0_731[k]
                    - f_11 * msh1_731[k]
                    + f_3 * pc_y[k] * msi_975[k];

        t_1255[k] = f_13 * lsi_780[k]
                    + f_8 * msh0_732[k]
                    - f_9 * msh1_732[k]
                    + f_3 * pc_y[k] * msi_976[k];

        t_1256[k] = f_13 * lsi_781[k]
                    + f_6 * msh0_733[k]
                    - f_7 * msh1_733[k]
                    + f_3 * pc_y[k] * msi_977[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pa_y, pc_y, lsk0_1007, lsi_782, lsi_783, \
                         lsk1_1007, msh0_734, msh1_734, msi_978, \
                         msi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_13 * lsi_782[k]
                    + f_4 * msh0_734[k]
                    - f_5 * msh1_734[k]
                    + f_3 * pc_y[k] * msi_978[k];

        t_1258[k] = f_13 * lsi_783[k]
                    + f_3 * pc_y[k] * msi_979[k];

        t_1259[k] = pa_y[k] * lsk0_1007[k]
                    - f_12 * pc_y[k] * lsk1_1007[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsk0,
                                                           const size_t lsi, const size_t lsk1,
                                                           const size_t msh0, const size_t msh1,
                                                           const size_t msi, const size_t ncols,
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
    const auto f_18 = 4.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_1008 = buffer.data(lsk0 + 1008);
    const auto *lsk0_1011 = buffer.data(lsk0 + 1011);
    const auto *lsk0_1014 = buffer.data(lsk0 + 1014);
    const auto *lsk0_1018 = buffer.data(lsk0 + 1018);
    const auto *lsk0_1023 = buffer.data(lsk0 + 1023);
    const auto *lsk0_1296 = buffer.data(lsk0 + 1296);
    const auto *lsk0_1299 = buffer.data(lsk0 + 1299);
    const auto *lsk0_1302 = buffer.data(lsk0 + 1302);
    const auto *lsk0_1306 = buffer.data(lsk0 + 1306);
    const auto *lsk0_1311 = buffer.data(lsk0 + 1311);
    const auto *lsk0_1324 = buffer.data(lsk0 + 1324);
    const auto *lsk0_1326 = buffer.data(lsk0 + 1326);
    const auto *lsk0_1327 = buffer.data(lsk0 + 1327);
    const auto *lsk0_1328 = buffer.data(lsk0 + 1328);
    const auto *lsk0_1329 = buffer.data(lsk0 + 1329);
    const auto *lsk0_1331 = buffer.data(lsk0 + 1331);
    const auto *lsk0_1337 = buffer.data(lsk0 + 1337);
    const auto *lsk0_1341 = buffer.data(lsk0 + 1341);
    const auto *lsk0_1344 = buffer.data(lsk0 + 1344);
    const auto *lsk0_1346 = buffer.data(lsk0 + 1346);
    const auto *lsk0_1349 = buffer.data(lsk0 + 1349);
    const auto *lsk0_1350 = buffer.data(lsk0 + 1350);
    const auto *lsk0_1352 = buffer.data(lsk0 + 1352);
    const auto *lsk0_1360 = buffer.data(lsk0 + 1360);
    const auto *lsk0_1362 = buffer.data(lsk0 + 1362);
    const auto *lsk0_1363 = buffer.data(lsk0 + 1363);
    const auto *lsk0_1364 = buffer.data(lsk0 + 1364);
    const auto *lsk0_1365 = buffer.data(lsk0 + 1365);
    const auto *lsk0_1367 = buffer.data(lsk0 + 1367);
    const auto *lsk0_1368 = buffer.data(lsk0 + 1368);
    const auto *lsk0_1371 = buffer.data(lsk0 + 1371);
    const auto *lsk0_1373 = buffer.data(lsk0 + 1373);
    const auto *lsk0_1374 = buffer.data(lsk0 + 1374);
    const auto *lsk0_1377 = buffer.data(lsk0 + 1377);
    const auto *lsk0_1378 = buffer.data(lsk0 + 1378);
    const auto *lsk0_1380 = buffer.data(lsk0 + 1380);

    const auto *lsi_756 = buffer.data(lsi + 756);
    const auto *lsi_783 = buffer.data(lsi + 783);
    const auto *lsi_784 = buffer.data(lsi + 784);
    const auto *lsi_787 = buffer.data(lsi + 787);
    const auto *lsi_789 = buffer.data(lsi + 789);
    const auto *lsi_790 = buffer.data(lsi + 790);
    const auto *lsi_793 = buffer.data(lsi + 793);
    const auto *lsi_794 = buffer.data(lsi + 794);
    const auto *lsi_798 = buffer.data(lsi + 798);
    const auto *lsi_805 = buffer.data(lsi + 805);
    const auto *lsi_811 = buffer.data(lsi + 811);
    const auto *lsi_812 = buffer.data(lsi + 812);
    const auto *lsi_814 = buffer.data(lsi + 814);
    const auto *lsi_815 = buffer.data(lsi + 815);
    const auto *lsi_817 = buffer.data(lsi + 817);
    const auto *lsi_818 = buffer.data(lsi + 818);
    const auto *lsi_821 = buffer.data(lsi + 821);
    const auto *lsi_826 = buffer.data(lsi + 826);
    const auto *lsi_839 = buffer.data(lsi + 839);
    const auto *lsi_840 = buffer.data(lsi + 840);
    const auto *lsi_842 = buffer.data(lsi + 842);
    const auto *lsi_845 = buffer.data(lsi + 845);
    const auto *lsi_849 = buffer.data(lsi + 849);
    const auto *lsi_980 = buffer.data(lsi + 980);
    const auto *lsi_985 = buffer.data(lsi + 985);
    const auto *lsi_989 = buffer.data(lsi + 989);
    const auto *lsi_994 = buffer.data(lsi + 994);
    const auto *lsi_1000 = buffer.data(lsi + 1000);
    const auto *lsi_1001 = buffer.data(lsi + 1001);
    const auto *lsi_1002 = buffer.data(lsi + 1002);
    const auto *lsi_1003 = buffer.data(lsi + 1003);
    const auto *lsi_1004 = buffer.data(lsi + 1004);
    const auto *lsi_1005 = buffer.data(lsi + 1005);
    const auto *lsi_1007 = buffer.data(lsi + 1007);
    const auto *lsi_1008 = buffer.data(lsi + 1008);
    const auto *lsi_1011 = buffer.data(lsi + 1011);
    const auto *lsi_1014 = buffer.data(lsi + 1014);
    const auto *lsi_1018 = buffer.data(lsi + 1018);
    const auto *lsi_1023 = buffer.data(lsi + 1023);
    const auto *lsi_1029 = buffer.data(lsi + 1029);
    const auto *lsi_1031 = buffer.data(lsi + 1031);
    const auto *lsi_1032 = buffer.data(lsi + 1032);
    const auto *lsi_1033 = buffer.data(lsi + 1033);
    const auto *lsi_1034 = buffer.data(lsi + 1034);
    const auto *lsi_1035 = buffer.data(lsi + 1035);
    const auto *lsi_1041 = buffer.data(lsi + 1041);
    const auto *lsi_1045 = buffer.data(lsi + 1045);
    const auto *lsi_1048 = buffer.data(lsi + 1048);
    const auto *lsi_1050 = buffer.data(lsi + 1050);
    const auto *lsi_1053 = buffer.data(lsi + 1053);
    const auto *lsi_1054 = buffer.data(lsi + 1054);
    const auto *lsi_1056 = buffer.data(lsi + 1056);
    const auto *lsi_1057 = buffer.data(lsi + 1057);
    const auto *lsi_1058 = buffer.data(lsi + 1058);
    const auto *lsi_1059 = buffer.data(lsi + 1059);
    const auto *lsi_1060 = buffer.data(lsi + 1060);
    const auto *lsi_1061 = buffer.data(lsi + 1061);
    const auto *lsi_1062 = buffer.data(lsi + 1062);
    const auto *lsi_1063 = buffer.data(lsi + 1063);
    const auto *lsi_1064 = buffer.data(lsi + 1064);
    const auto *lsi_1067 = buffer.data(lsi + 1067);
    const auto *lsi_1069 = buffer.data(lsi + 1069);
    const auto *lsi_1070 = buffer.data(lsi + 1070);
    const auto *lsi_1073 = buffer.data(lsi + 1073);
    const auto *lsi_1074 = buffer.data(lsi + 1074);
    const auto *lsi_1076 = buffer.data(lsi + 1076);

    const auto *lsk1_1008 = buffer.data(lsk1 + 1008);
    const auto *lsk1_1011 = buffer.data(lsk1 + 1011);
    const auto *lsk1_1014 = buffer.data(lsk1 + 1014);
    const auto *lsk1_1018 = buffer.data(lsk1 + 1018);
    const auto *lsk1_1023 = buffer.data(lsk1 + 1023);
    const auto *lsk1_1296 = buffer.data(lsk1 + 1296);
    const auto *lsk1_1299 = buffer.data(lsk1 + 1299);
    const auto *lsk1_1302 = buffer.data(lsk1 + 1302);
    const auto *lsk1_1306 = buffer.data(lsk1 + 1306);
    const auto *lsk1_1311 = buffer.data(lsk1 + 1311);
    const auto *lsk1_1324 = buffer.data(lsk1 + 1324);
    const auto *lsk1_1326 = buffer.data(lsk1 + 1326);
    const auto *lsk1_1327 = buffer.data(lsk1 + 1327);
    const auto *lsk1_1328 = buffer.data(lsk1 + 1328);
    const auto *lsk1_1329 = buffer.data(lsk1 + 1329);
    const auto *lsk1_1331 = buffer.data(lsk1 + 1331);
    const auto *lsk1_1337 = buffer.data(lsk1 + 1337);
    const auto *lsk1_1341 = buffer.data(lsk1 + 1341);
    const auto *lsk1_1344 = buffer.data(lsk1 + 1344);
    const auto *lsk1_1346 = buffer.data(lsk1 + 1346);
    const auto *lsk1_1349 = buffer.data(lsk1 + 1349);
    const auto *lsk1_1350 = buffer.data(lsk1 + 1350);
    const auto *lsk1_1352 = buffer.data(lsk1 + 1352);
    const auto *lsk1_1360 = buffer.data(lsk1 + 1360);
    const auto *lsk1_1362 = buffer.data(lsk1 + 1362);
    const auto *lsk1_1363 = buffer.data(lsk1 + 1363);
    const auto *lsk1_1364 = buffer.data(lsk1 + 1364);
    const auto *lsk1_1365 = buffer.data(lsk1 + 1365);
    const auto *lsk1_1367 = buffer.data(lsk1 + 1367);
    const auto *lsk1_1368 = buffer.data(lsk1 + 1368);
    const auto *lsk1_1371 = buffer.data(lsk1 + 1371);
    const auto *lsk1_1373 = buffer.data(lsk1 + 1373);
    const auto *lsk1_1374 = buffer.data(lsk1 + 1374);
    const auto *lsk1_1377 = buffer.data(lsk1 + 1377);
    const auto *lsk1_1378 = buffer.data(lsk1 + 1378);
    const auto *lsk1_1380 = buffer.data(lsk1 + 1380);

    const auto *msh0_735 = buffer.data(msh0 + 735);
    const auto *msh0_736 = buffer.data(msh0 + 736);
    const auto *msh0_737 = buffer.data(msh0 + 737);
    const auto *msh0_738 = buffer.data(msh0 + 738);
    const auto *msh0_739 = buffer.data(msh0 + 739);
    const auto *msh0_740 = buffer.data(msh0 + 740);
    const auto *msh0_741 = buffer.data(msh0 + 741);
    const auto *msh0_742 = buffer.data(msh0 + 742);
    const auto *msh0_743 = buffer.data(msh0 + 743);
    const auto *msh0_744 = buffer.data(msh0 + 744);
    const auto *msh0_749 = buffer.data(msh0 + 749);
    const auto *msh0_750 = buffer.data(msh0 + 750);
    const auto *msh0_751 = buffer.data(msh0 + 751);
    const auto *msh0_752 = buffer.data(msh0 + 752);
    const auto *msh0_753 = buffer.data(msh0 + 753);
    const auto *msh0_754 = buffer.data(msh0 + 754);
    const auto *msh0_755 = buffer.data(msh0 + 755);
    const auto *msh0_756 = buffer.data(msh0 + 756);
    const auto *msh0_758 = buffer.data(msh0 + 758);
    const auto *msh0_759 = buffer.data(msh0 + 759);
    const auto *msh0_761 = buffer.data(msh0 + 761);
    const auto *msh0_762 = buffer.data(msh0 + 762);
    const auto *msh0_763 = buffer.data(msh0 + 763);
    const auto *msh0_765 = buffer.data(msh0 + 765);

    const auto *msh1_735 = buffer.data(msh1 + 735);
    const auto *msh1_736 = buffer.data(msh1 + 736);
    const auto *msh1_737 = buffer.data(msh1 + 737);
    const auto *msh1_738 = buffer.data(msh1 + 738);
    const auto *msh1_739 = buffer.data(msh1 + 739);
    const auto *msh1_740 = buffer.data(msh1 + 740);
    const auto *msh1_741 = buffer.data(msh1 + 741);
    const auto *msh1_742 = buffer.data(msh1 + 742);
    const auto *msh1_743 = buffer.data(msh1 + 743);
    const auto *msh1_744 = buffer.data(msh1 + 744);
    const auto *msh1_749 = buffer.data(msh1 + 749);
    const auto *msh1_750 = buffer.data(msh1 + 750);
    const auto *msh1_751 = buffer.data(msh1 + 751);
    const auto *msh1_752 = buffer.data(msh1 + 752);
    const auto *msh1_753 = buffer.data(msh1 + 753);
    const auto *msh1_754 = buffer.data(msh1 + 754);
    const auto *msh1_755 = buffer.data(msh1 + 755);
    const auto *msh1_756 = buffer.data(msh1 + 756);
    const auto *msh1_758 = buffer.data(msh1 + 758);
    const auto *msh1_759 = buffer.data(msh1 + 759);
    const auto *msh1_761 = buffer.data(msh1 + 761);
    const auto *msh1_762 = buffer.data(msh1 + 762);
    const auto *msh1_763 = buffer.data(msh1 + 763);
    const auto *msh1_765 = buffer.data(msh1 + 765);

    const auto *msi_980 = buffer.data(msi + 980);
    const auto *msi_981 = buffer.data(msi + 981);
    const auto *msi_982 = buffer.data(msi + 982);
    const auto *msi_983 = buffer.data(msi + 983);
    const auto *msi_984 = buffer.data(msi + 984);
    const auto *msi_985 = buffer.data(msi + 985);
    const auto *msi_986 = buffer.data(msi + 986);
    const auto *msi_987 = buffer.data(msi + 987);
    const auto *msi_988 = buffer.data(msi + 988);
    const auto *msi_989 = buffer.data(msi + 989);
    const auto *msi_990 = buffer.data(msi + 990);
    const auto *msi_991 = buffer.data(msi + 991);
    const auto *msi_992 = buffer.data(msi + 992);
    const auto *msi_993 = buffer.data(msi + 993);
    const auto *msi_994 = buffer.data(msi + 994);
    const auto *msi_1000 = buffer.data(msi + 1000);
    const auto *msi_1001 = buffer.data(msi + 1001);
    const auto *msi_1002 = buffer.data(msi + 1002);
    const auto *msi_1003 = buffer.data(msi + 1003);
    const auto *msi_1004 = buffer.data(msi + 1004);
    const auto *msi_1005 = buffer.data(msi + 1005);
    const auto *msi_1006 = buffer.data(msi + 1006);
    const auto *msi_1007 = buffer.data(msi + 1007);
    const auto *msi_1008 = buffer.data(msi + 1008);
    const auto *msi_1009 = buffer.data(msi + 1009);
    const auto *msi_1010 = buffer.data(msi + 1010);
    const auto *msi_1011 = buffer.data(msi + 1011);
    const auto *msi_1013 = buffer.data(msi + 1013);
    const auto *msi_1014 = buffer.data(msi + 1014);
    const auto *msi_1015 = buffer.data(msi + 1015);
    const auto *msi_1017 = buffer.data(msi + 1017);
    const auto *msi_1018 = buffer.data(msi + 1018);
    const auto *msi_1019 = buffer.data(msi + 1019);
    const auto *msi_1020 = buffer.data(msi + 1020);
    const auto *msi_1022 = buffer.data(msi + 1022);
    const auto *msi_1023 = buffer.data(msi + 1023);
    const auto *msi_1029 = buffer.data(msi + 1029);
    const auto *msi_1031 = buffer.data(msi + 1031);
    const auto *msi_1032 = buffer.data(msi + 1032);
    const auto *msi_1033 = buffer.data(msi + 1033);
    const auto *msi_1034 = buffer.data(msi + 1034);
    const auto *msi_1035 = buffer.data(msi + 1035);
    const auto *msi_1036 = buffer.data(msi + 1036);
    const auto *msi_1038 = buffer.data(msi + 1038);
    const auto *msi_1039 = buffer.data(msi + 1039);
    const auto *msi_1041 = buffer.data(msi + 1041);
    const auto *msi_1042 = buffer.data(msi + 1042);
    const auto *msi_1045 = buffer.data(msi + 1045);
    const auto *msi_1046 = buffer.data(msi + 1046);
    const auto *msi_1050 = buffer.data(msi + 1050);
    const auto *msi_1057 = buffer.data(msi + 1057);
    const auto *msi_1058 = buffer.data(msi + 1058);
    const auto *msi_1059 = buffer.data(msi + 1059);
    const auto *msi_1060 = buffer.data(msi + 1060);
    const auto *msi_1061 = buffer.data(msi + 1061);
    const auto *msi_1062 = buffer.data(msi + 1062);
    const auto *msi_1063 = buffer.data(msi + 1063);
    const auto *msi_1064 = buffer.data(msi + 1064);
    const auto *msi_1066 = buffer.data(msi + 1066);
    const auto *msi_1067 = buffer.data(msi + 1067);
    const auto *msi_1069 = buffer.data(msi + 1069);
    const auto *msi_1070 = buffer.data(msi + 1070);
    const auto *msi_1073 = buffer.data(msi + 1073);

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, t_1264, pc_x, pc_y, pc_z, lsi_756, \
                         lsi_980, msh0_735, msh1_735, msi_980, msi_981, \
                         msi_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_14 * lsi_980[k]
                    + f_1 * msh0_735[k]
                    - f_2 * msh1_735[k]
                    + f_3 * pc_x[k] * msi_980[k];

        t_1261[k] = f_3 * pc_y[k] * msi_980[k];

        t_1262[k] = f_21 * lsi_756[k]
                    + f_3 * pc_z[k] * msi_980[k];

        t_1263[k] = f_4 * msh0_735[k]
                    - f_5 * msh1_735[k]
                    + f_3 * pc_y[k] * msi_981[k];

        t_1264[k] = f_3 * pc_y[k] * msi_982[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, t_1268, pc_x, pc_y, lsi_985, msh0_736, \
                         msh0_737, msh0_740, msh1_736, msh1_737, msh1_740, msi_983, msi_984, \
                         msi_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_14 * lsi_985[k]
                    + f_10 * msh0_740[k]
                    - f_11 * msh1_740[k]
                    + f_3 * pc_x[k] * msi_985[k];

        t_1266[k] = f_6 * msh0_736[k]
                    - f_7 * msh1_736[k]
                    + f_3 * pc_y[k] * msi_983[k];

        t_1267[k] = f_4 * msh0_737[k]
                    - f_5 * msh1_737[k]
                    + f_3 * pc_y[k] * msi_984[k];

        t_1268[k] = f_3 * pc_y[k] * msi_985[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, pc_x, pc_y, lsi_989, msh0_738, msh0_739, \
                         msh0_744, msh1_738, msh1_739, msh1_744, msi_986, msi_987, \
                         msi_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = f_14 * lsi_989[k]
                    + f_8 * msh0_744[k]
                    - f_9 * msh1_744[k]
                    + f_3 * pc_x[k] * msi_989[k];

        t_1270[k] = f_8 * msh0_738[k]
                    - f_9 * msh1_738[k]
                    + f_3 * pc_y[k] * msi_986[k];

        t_1271[k] = f_6 * msh0_739[k]
                    - f_7 * msh1_739[k]
                    + f_3 * pc_y[k] * msi_987[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pc_x, pc_y, lsi_994, msh0_740, msh0_749, \
                         msh1_740, msh1_749, msi_988, msi_989, \
                         msi_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * msh0_740[k]
                    - f_5 * msh1_740[k]
                    + f_3 * pc_y[k] * msi_988[k];

        t_1273[k] = f_3 * pc_y[k] * msi_989[k];

        t_1274[k] = f_14 * lsi_994[k]
                    + f_6 * msh0_749[k]
                    - f_7 * msh1_749[k]
                    + f_3 * pc_x[k] * msi_994[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pc_y, msh0_741, msh0_742, msh0_743, msh1_741, \
                         msh1_742, msh1_743, msi_990, msi_991, \
                         msi_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_10 * msh0_741[k]
                    - f_11 * msh1_741[k]
                    + f_3 * pc_y[k] * msi_990[k];

        t_1276[k] = f_8 * msh0_742[k]
                    - f_9 * msh1_742[k]
                    + f_3 * pc_y[k] * msi_991[k];

        t_1277[k] = f_6 * msh0_743[k]
                    - f_7 * msh1_743[k]
                    + f_3 * pc_y[k] * msi_992[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, t_1281, pc_x, pc_y, lsi_1000, lsi_1001, \
                         msh0_744, msh0_755, msh1_744, msh1_755, msi_993, msi_994, msi_1000, \
                         msi_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_4 * msh0_744[k]
                    - f_5 * msh1_744[k]
                    + f_3 * pc_y[k] * msi_993[k];

        t_1279[k] = f_3 * pc_y[k] * msi_994[k];

        t_1280[k] = f_14 * lsi_1000[k]
                    + f_4 * msh0_755[k]
                    - f_5 * msh1_755[k]
                    + f_3 * pc_x[k] * msi_1000[k];

        t_1281[k] = f_14 * lsi_1001[k]
                    + f_3 * pc_x[k] * msi_1001[k];
    }

#pragma omp simd aligned(t_1282, t_1283, t_1284, t_1285, t_1286, pc_x, pc_y, lsi_1002, \
                         lsi_1003, lsi_1004, lsi_1005, msi_1000, msi_1002, msi_1003, msi_1004, \
                         msi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1282[k] = f_14 * lsi_1002[k]
                    + f_3 * pc_x[k] * msi_1002[k];

        t_1283[k] = f_14 * lsi_1003[k]
                    + f_3 * pc_x[k] * msi_1003[k];

        t_1284[k] = f_14 * lsi_1004[k]
                    + f_3 * pc_x[k] * msi_1004[k];

        t_1285[k] = f_14 * lsi_1005[k]
                    + f_3 * pc_x[k] * msi_1005[k];

        t_1286[k] = f_3 * pc_y[k] * msi_1000[k];
    }

#pragma omp simd aligned(t_1287, t_1288, t_1289, pc_x, pc_y, lsi_1007, msh0_750, msh0_751, \
                         msh1_750, msh1_751, msi_1001, msi_1002, \
                         msi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1287[k] = f_14 * lsi_1007[k]
                    + f_3 * pc_x[k] * msi_1007[k];

        t_1288[k] = f_1 * msh0_750[k]
                    - f_2 * msh1_750[k]
                    + f_3 * pc_y[k] * msi_1001[k];

        t_1289[k] = f_19 * msh0_751[k]
                    - f_20 * msh1_751[k]
                    + f_3 * pc_y[k] * msi_1002[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, pc_y, msh0_752, msh0_753, msh0_754, msh1_752, \
                         msh1_753, msh1_754, msi_1003, msi_1004, \
                         msi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_10 * msh0_752[k]
                    - f_11 * msh1_752[k]
                    + f_3 * pc_y[k] * msi_1003[k];

        t_1291[k] = f_8 * msh0_753[k]
                    - f_9 * msh1_753[k]
                    + f_3 * pc_y[k] * msi_1004[k];

        t_1292[k] = f_6 * msh0_754[k]
                    - f_7 * msh1_754[k]
                    + f_3 * pc_y[k] * msi_1005[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pa_x, pc_x, pc_y, pc_z, lsk0_1296, \
                         lsi_783, lsi_1008, lsk1_1296, msh0_755, msh1_755, msi_1006, \
                         msi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_4 * msh0_755[k]
                    - f_5 * msh1_755[k]
                    + f_3 * pc_y[k] * msi_1006[k];

        t_1294[k] = f_3 * pc_y[k] * msi_1007[k];

        t_1295[k] = f_21 * lsi_783[k]
                    + f_1 * msh0_755[k]
                    - f_2 * msh1_755[k]
                    + f_3 * pc_z[k] * msi_1007[k];

        t_1296[k] = pa_x[k] * lsk0_1296[k]
                    + f_21 * lsi_1008[k]
                    - f_12 * pc_x[k] * lsk1_1296[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, t_1300, pa_x, pc_x, pc_y, pc_z, lsk0_1299, \
                         lsi_784, lsi_1011, lsk1_1299, msi_1008, \
                         msi_1009 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = f_18 * lsi_784[k]
                    + f_3 * pc_y[k] * msi_1008[k];

        t_1298[k] = f_3 * pc_z[k] * msi_1008[k];

        t_1299[k] = pa_x[k] * lsk0_1299[k]
                    + f_17 * lsi_1011[k]
                    - f_12 * pc_x[k] * lsk1_1299[k];

        t_1300[k] = f_3 * pc_z[k] * msi_1009[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, pa_x, pc_x, pc_z, lsk0_1302, lsi_1014, \
                         lsk1_1302, msh0_756, msh1_756, msi_1010, \
                         msi_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_4 * msh0_756[k]
                    - f_5 * msh1_756[k]
                    + f_3 * pc_z[k] * msi_1010[k];

        t_1302[k] = pa_x[k] * lsk0_1302[k]
                    + f_16 * lsi_1014[k]
                    - f_12 * pc_x[k] * lsk1_1302[k];

        t_1303[k] = f_3 * pc_z[k] * msi_1011[k];
    }

#pragma omp simd aligned(t_1304, t_1305, t_1306, t_1307, pa_x, pc_x, pc_y, pc_z, lsk0_1306, \
                         lsi_789, lsi_1018, lsk1_1306, msh0_758, msh1_758, msi_1013, \
                         msi_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1304[k] = f_18 * lsi_789[k]
                    + f_3 * pc_y[k] * msi_1013[k];

        t_1305[k] = f_6 * msh0_758[k]
                    - f_7 * msh1_758[k]
                    + f_3 * pc_z[k] * msi_1013[k];

        t_1306[k] = pa_x[k] * lsk0_1306[k]
                    + f_15 * lsi_1018[k]
                    - f_12 * pc_x[k] * lsk1_1306[k];

        t_1307[k] = f_3 * pc_z[k] * msi_1014[k];
    }

#pragma omp simd aligned(t_1308, t_1309, t_1310, pc_y, pc_z, lsi_793, msh0_759, msh0_761, \
                         msh1_759, msh1_761, msi_1015, msi_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1308[k] = f_4 * msh0_759[k]
                    - f_5 * msh1_759[k]
                    + f_3 * pc_z[k] * msi_1015[k];

        t_1309[k] = f_18 * lsi_793[k]
                    + f_3 * pc_y[k] * msi_1017[k];

        t_1310[k] = f_8 * msh0_761[k]
                    - f_9 * msh1_761[k]
                    + f_3 * pc_z[k] * msi_1017[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, pa_x, pc_x, pc_z, lsk0_1311, lsi_1023, \
                         lsk1_1311, msh0_762, msh1_762, msi_1018, \
                         msi_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = pa_x[k] * lsk0_1311[k]
                    + f_14 * lsi_1023[k]
                    - f_12 * pc_x[k] * lsk1_1311[k];

        t_1312[k] = f_3 * pc_z[k] * msi_1018[k];

        t_1313[k] = f_4 * msh0_762[k]
                    - f_5 * msh1_762[k]
                    + f_3 * pc_z[k] * msi_1019[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, t_1317, pc_x, pc_y, pc_z, lsi_798, lsi_1029, \
                         msh0_763, msh0_765, msh1_763, msh1_765, msi_1020, msi_1022, \
                         msi_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_6 * msh0_763[k]
                    - f_7 * msh1_763[k]
                    + f_3 * pc_z[k] * msi_1020[k];

        t_1315[k] = f_18 * lsi_798[k]
                    + f_3 * pc_y[k] * msi_1022[k];

        t_1316[k] = f_10 * msh0_765[k]
                    - f_11 * msh1_765[k]
                    + f_3 * pc_z[k] * msi_1022[k];

        t_1317[k] = f_13 * lsi_1029[k]
                    + f_3 * pc_x[k] * msi_1029[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, t_1321, t_1322, pc_x, pc_z, lsi_1031, \
                         lsi_1032, lsi_1033, lsi_1034, msi_1023, msi_1031, msi_1032, msi_1033, \
                         msi_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_3 * pc_z[k] * msi_1023[k];

        t_1319[k] = f_13 * lsi_1031[k]
                    + f_3 * pc_x[k] * msi_1031[k];

        t_1320[k] = f_13 * lsi_1032[k]
                    + f_3 * pc_x[k] * msi_1032[k];

        t_1321[k] = f_13 * lsi_1033[k]
                    + f_3 * pc_x[k] * msi_1033[k];

        t_1322[k] = f_13 * lsi_1034[k]
                    + f_3 * pc_x[k] * msi_1034[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, t_1326, pa_x, pc_x, pc_z, lsk0_1324, \
                         lsk0_1326, lsi_1035, lsk1_1324, lsk1_1326, msi_1029, \
                         msi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = f_13 * lsi_1035[k]
                    + f_3 * pc_x[k] * msi_1035[k];

        t_1324[k] = pa_x[k] * lsk0_1324[k]
                    - f_12 * pc_x[k] * lsk1_1324[k];

        t_1325[k] = f_3 * pc_z[k] * msi_1029[k];

        t_1326[k] = pa_x[k] * lsk0_1326[k]
                    - f_12 * pc_x[k] * lsk1_1326[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, t_1330, pa_x, pc_x, pc_y, lsk0_1327, \
                         lsk0_1328, lsk0_1329, lsi_811, lsk1_1327, lsk1_1328, lsk1_1329, \
                         msi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = pa_x[k] * lsk0_1327[k]
                    - f_12 * pc_x[k] * lsk1_1327[k];

        t_1328[k] = pa_x[k] * lsk0_1328[k]
                    - f_12 * pc_x[k] * lsk1_1328[k];

        t_1329[k] = pa_x[k] * lsk0_1329[k]
                    - f_12 * pc_x[k] * lsk1_1329[k];

        t_1330[k] = f_18 * lsi_811[k]
                    + f_3 * pc_y[k] * msi_1035[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, pa_x, pa_z, pc_x, pc_y, pc_z, \
                         lsk0_1008, lsk0_1331, lsi_784, lsi_812, lsk1_1008, lsk1_1331, \
                         msi_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = pa_x[k] * lsk0_1331[k]
                    - f_12 * pc_x[k] * lsk1_1331[k];

        t_1332[k] = pa_z[k] * lsk0_1008[k]
                    - f_12 * pc_z[k] * lsk1_1008[k];

        t_1333[k] = f_21 * lsi_812[k]
                    + f_3 * pc_y[k] * msi_1036[k];

        t_1334[k] = f_13 * lsi_784[k]
                    + f_3 * pc_z[k] * msi_1036[k];
    }

#pragma omp simd aligned(t_1335, t_1336, t_1337, pa_x, pa_z, pc_x, pc_y, pc_z, lsk0_1011, \
                         lsk0_1337, lsi_814, lsi_1041, lsk1_1011, lsk1_1337, \
                         msi_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1335[k] = pa_z[k] * lsk0_1011[k]
                    - f_12 * pc_z[k] * lsk1_1011[k];

        t_1336[k] = f_21 * lsi_814[k]
                    + f_3 * pc_y[k] * msi_1038[k];

        t_1337[k] = pa_x[k] * lsk0_1337[k]
                    + f_17 * lsi_1041[k]
                    - f_12 * pc_x[k] * lsk1_1337[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, pa_z, pc_y, pc_z, lsk0_1014, lsi_787, \
                         lsi_817, lsk1_1014, msi_1039, msi_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = pa_z[k] * lsk0_1014[k]
                    - f_12 * pc_z[k] * lsk1_1014[k];

        t_1339[k] = f_13 * lsi_787[k]
                    + f_3 * pc_z[k] * msi_1039[k];

        t_1340[k] = f_21 * lsi_817[k]
                    + f_3 * pc_y[k] * msi_1041[k];
    }

#pragma omp simd aligned(t_1341, t_1342, t_1343, pa_x, pa_z, pc_x, pc_z, lsk0_1018, lsk0_1341, \
                         lsi_790, lsi_1045, lsk1_1018, lsk1_1341, \
                         msi_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1341[k] = pa_x[k] * lsk0_1341[k]
                    + f_16 * lsi_1045[k]
                    - f_12 * pc_x[k] * lsk1_1341[k];

        t_1342[k] = pa_z[k] * lsk0_1018[k]
                    - f_12 * pc_z[k] * lsk1_1018[k];

        t_1343[k] = f_13 * lsi_790[k]
                    + f_3 * pc_z[k] * msi_1042[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, pa_x, pc_x, pc_y, lsk0_1344, lsk0_1346, \
                         lsi_821, lsi_1048, lsi_1050, lsk1_1344, lsk1_1346, \
                         msi_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = pa_x[k] * lsk0_1344[k]
                    + f_15 * lsi_1048[k]
                    - f_12 * pc_x[k] * lsk1_1344[k];

        t_1345[k] = f_21 * lsi_821[k]
                    + f_3 * pc_y[k] * msi_1045[k];

        t_1346[k] = pa_x[k] * lsk0_1346[k]
                    + f_15 * lsi_1050[k]
                    - f_12 * pc_x[k] * lsk1_1346[k];
    }

#pragma omp simd aligned(t_1347, t_1348, t_1349, pa_x, pa_z, pc_x, pc_z, lsk0_1023, lsk0_1349, \
                         lsi_794, lsi_1053, lsk1_1023, lsk1_1349, \
                         msi_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1347[k] = pa_z[k] * lsk0_1023[k]
                    - f_12 * pc_z[k] * lsk1_1023[k];

        t_1348[k] = f_13 * lsi_794[k]
                    + f_3 * pc_z[k] * msi_1046[k];

        t_1349[k] = pa_x[k] * lsk0_1349[k]
                    + f_14 * lsi_1053[k]
                    - f_12 * pc_x[k] * lsk1_1349[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, pa_x, pc_x, pc_y, lsk0_1350, lsk0_1352, \
                         lsi_826, lsi_1054, lsi_1056, lsk1_1350, lsk1_1352, \
                         msi_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = pa_x[k] * lsk0_1350[k]
                    + f_14 * lsi_1054[k]
                    - f_12 * pc_x[k] * lsk1_1350[k];

        t_1351[k] = f_21 * lsi_826[k]
                    + f_3 * pc_y[k] * msi_1050[k];

        t_1352[k] = pa_x[k] * lsk0_1352[k]
                    + f_14 * lsi_1056[k]
                    - f_12 * pc_x[k] * lsk1_1352[k];
    }

#pragma omp simd aligned(t_1353, t_1354, t_1355, t_1356, t_1357, pc_x, lsi_1057, lsi_1058, \
                         lsi_1059, lsi_1060, lsi_1061, msi_1057, msi_1058, msi_1059, msi_1060, \
                         msi_1061 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1353[k] = f_13 * lsi_1057[k]
                    + f_3 * pc_x[k] * msi_1057[k];

        t_1354[k] = f_13 * lsi_1058[k]
                    + f_3 * pc_x[k] * msi_1058[k];

        t_1355[k] = f_13 * lsi_1059[k]
                    + f_3 * pc_x[k] * msi_1059[k];

        t_1356[k] = f_13 * lsi_1060[k]
                    + f_3 * pc_x[k] * msi_1060[k];

        t_1357[k] = f_13 * lsi_1061[k]
                    + f_3 * pc_x[k] * msi_1061[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pa_x, pc_x, pc_z, lsk0_1360, lsi_805, \
                         lsi_1062, lsi_1063, lsk1_1360, msi_1057, msi_1062, \
                         msi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_13 * lsi_1062[k]
                    + f_3 * pc_x[k] * msi_1062[k];

        t_1359[k] = f_13 * lsi_1063[k]
                    + f_3 * pc_x[k] * msi_1063[k];

        t_1360[k] = pa_x[k] * lsk0_1360[k]
                    - f_12 * pc_x[k] * lsk1_1360[k];

        t_1361[k] = f_13 * lsi_805[k]
                    + f_3 * pc_z[k] * msi_1057[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, pa_x, pc_x, lsk0_1362, lsk0_1363, \
                         lsk0_1364, lsk0_1365, lsk1_1362, lsk1_1363, lsk1_1364, \
                         lsk1_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = pa_x[k] * lsk0_1362[k]
                    - f_12 * pc_x[k] * lsk1_1362[k];

        t_1363[k] = pa_x[k] * lsk0_1363[k]
                    - f_12 * pc_x[k] * lsk1_1363[k];

        t_1364[k] = pa_x[k] * lsk0_1364[k]
                    - f_12 * pc_x[k] * lsk1_1364[k];

        t_1365[k] = pa_x[k] * lsk0_1365[k]
                    - f_12 * pc_x[k] * lsk1_1365[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, t_1369, pa_x, pc_x, pc_y, lsk0_1367, \
                         lsk0_1368, lsi_839, lsi_840, lsi_1064, lsk1_1367, lsk1_1368, \
                         msi_1063, msi_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_21 * lsi_839[k]
                    + f_3 * pc_y[k] * msi_1063[k];

        t_1367[k] = pa_x[k] * lsk0_1367[k]
                    - f_12 * pc_x[k] * lsk1_1367[k];

        t_1368[k] = pa_x[k] * lsk0_1368[k]
                    + f_21 * lsi_1064[k]
                    - f_12 * pc_x[k] * lsk1_1368[k];

        t_1369[k] = f_22 * lsi_840[k]
                    + f_3 * pc_y[k] * msi_1064[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, pa_x, pc_x, pc_y, pc_z, lsk0_1371, lsi_812, \
                         lsi_842, lsi_1067, lsk1_1371, msi_1064, \
                         msi_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_14 * lsi_812[k]
                    + f_3 * pc_z[k] * msi_1064[k];

        t_1371[k] = pa_x[k] * lsk0_1371[k]
                    + f_17 * lsi_1067[k]
                    - f_12 * pc_x[k] * lsk1_1371[k];

        t_1372[k] = f_22 * lsi_842[k]
                    + f_3 * pc_y[k] * msi_1066[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pa_x, pc_x, pc_z, lsk0_1373, lsk0_1374, \
                         lsi_815, lsi_1069, lsi_1070, lsk1_1373, lsk1_1374, \
                         msi_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = pa_x[k] * lsk0_1373[k]
                    + f_17 * lsi_1069[k]
                    - f_12 * pc_x[k] * lsk1_1373[k];

        t_1374[k] = pa_x[k] * lsk0_1374[k]
                    + f_16 * lsi_1070[k]
                    - f_12 * pc_x[k] * lsk1_1374[k];

        t_1375[k] = f_14 * lsi_815[k]
                    + f_3 * pc_z[k] * msi_1067[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, pa_x, pc_x, pc_y, lsk0_1377, lsk0_1378, \
                         lsi_845, lsi_1073, lsi_1074, lsk1_1377, lsk1_1378, \
                         msi_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_22 * lsi_845[k]
                    + f_3 * pc_y[k] * msi_1069[k];

        t_1377[k] = pa_x[k] * lsk0_1377[k]
                    + f_16 * lsi_1073[k]
                    - f_12 * pc_x[k] * lsk1_1377[k];

        t_1378[k] = pa_x[k] * lsk0_1378[k]
                    + f_15 * lsi_1074[k]
                    - f_12 * pc_x[k] * lsk1_1378[k];
    }

#pragma omp simd aligned(t_1379, t_1380, t_1381, pa_x, pc_x, pc_y, pc_z, lsk0_1380, lsi_818, \
                         lsi_849, lsi_1076, lsk1_1380, msi_1070, \
                         msi_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1379[k] = f_14 * lsi_818[k]
                    + f_3 * pc_z[k] * msi_1070[k];

        t_1380[k] = pa_x[k] * lsk0_1380[k]
                    + f_15 * lsi_1076[k]
                    - f_12 * pc_x[k] * lsk1_1380[k];

        t_1381[k] = f_22 * lsi_849[k]
                    + f_3 * pc_y[k] * msi_1073[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsk0,
                                                           const size_t lsi, const size_t lsk1,
                                                           const size_t msi, const size_t ncols,
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
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_1382 = buffer.data(lsk0 + 1382);
    const auto *lsk0_1383 = buffer.data(lsk0 + 1383);
    const auto *lsk0_1385 = buffer.data(lsk0 + 1385);
    const auto *lsk0_1386 = buffer.data(lsk0 + 1386);
    const auto *lsk0_1388 = buffer.data(lsk0 + 1388);
    const auto *lsk0_1396 = buffer.data(lsk0 + 1396);
    const auto *lsk0_1398 = buffer.data(lsk0 + 1398);
    const auto *lsk0_1399 = buffer.data(lsk0 + 1399);
    const auto *lsk0_1400 = buffer.data(lsk0 + 1400);
    const auto *lsk0_1401 = buffer.data(lsk0 + 1401);
    const auto *lsk0_1403 = buffer.data(lsk0 + 1403);
    const auto *lsk0_1404 = buffer.data(lsk0 + 1404);
    const auto *lsk0_1407 = buffer.data(lsk0 + 1407);
    const auto *lsk0_1409 = buffer.data(lsk0 + 1409);
    const auto *lsk0_1410 = buffer.data(lsk0 + 1410);
    const auto *lsk0_1413 = buffer.data(lsk0 + 1413);
    const auto *lsk0_1414 = buffer.data(lsk0 + 1414);
    const auto *lsk0_1416 = buffer.data(lsk0 + 1416);
    const auto *lsk0_1418 = buffer.data(lsk0 + 1418);
    const auto *lsk0_1419 = buffer.data(lsk0 + 1419);
    const auto *lsk0_1421 = buffer.data(lsk0 + 1421);
    const auto *lsk0_1422 = buffer.data(lsk0 + 1422);
    const auto *lsk0_1424 = buffer.data(lsk0 + 1424);
    const auto *lsk0_1432 = buffer.data(lsk0 + 1432);
    const auto *lsk0_1434 = buffer.data(lsk0 + 1434);
    const auto *lsk0_1435 = buffer.data(lsk0 + 1435);
    const auto *lsk0_1436 = buffer.data(lsk0 + 1436);
    const auto *lsk0_1437 = buffer.data(lsk0 + 1437);
    const auto *lsk0_1439 = buffer.data(lsk0 + 1439);
    const auto *lsk0_1440 = buffer.data(lsk0 + 1440);
    const auto *lsk0_1443 = buffer.data(lsk0 + 1443);
    const auto *lsk0_1445 = buffer.data(lsk0 + 1445);
    const auto *lsk0_1446 = buffer.data(lsk0 + 1446);
    const auto *lsk0_1449 = buffer.data(lsk0 + 1449);
    const auto *lsk0_1450 = buffer.data(lsk0 + 1450);
    const auto *lsk0_1452 = buffer.data(lsk0 + 1452);
    const auto *lsk0_1454 = buffer.data(lsk0 + 1454);
    const auto *lsk0_1455 = buffer.data(lsk0 + 1455);
    const auto *lsk0_1457 = buffer.data(lsk0 + 1457);
    const auto *lsk0_1458 = buffer.data(lsk0 + 1458);
    const auto *lsk0_1460 = buffer.data(lsk0 + 1460);
    const auto *lsk0_1468 = buffer.data(lsk0 + 1468);
    const auto *lsk0_1470 = buffer.data(lsk0 + 1470);
    const auto *lsk0_1471 = buffer.data(lsk0 + 1471);
    const auto *lsk0_1472 = buffer.data(lsk0 + 1472);
    const auto *lsk0_1473 = buffer.data(lsk0 + 1473);
    const auto *lsk0_1475 = buffer.data(lsk0 + 1475);
    const auto *lsk0_1476 = buffer.data(lsk0 + 1476);
    const auto *lsk0_1479 = buffer.data(lsk0 + 1479);
    const auto *lsk0_1481 = buffer.data(lsk0 + 1481);
    const auto *lsk0_1482 = buffer.data(lsk0 + 1482);
    const auto *lsk0_1485 = buffer.data(lsk0 + 1485);
    const auto *lsk0_1486 = buffer.data(lsk0 + 1486);
    const auto *lsk0_1488 = buffer.data(lsk0 + 1488);
    const auto *lsk0_1490 = buffer.data(lsk0 + 1490);
    const auto *lsk0_1491 = buffer.data(lsk0 + 1491);
    const auto *lsk0_1493 = buffer.data(lsk0 + 1493);
    const auto *lsk0_1494 = buffer.data(lsk0 + 1494);
    const auto *lsk0_1496 = buffer.data(lsk0 + 1496);

    const auto *lsi_822 = buffer.data(lsi + 822);
    const auto *lsi_833 = buffer.data(lsi + 833);
    const auto *lsi_840 = buffer.data(lsi + 840);
    const auto *lsi_843 = buffer.data(lsi + 843);
    const auto *lsi_846 = buffer.data(lsi + 846);
    const auto *lsi_850 = buffer.data(lsi + 850);
    const auto *lsi_854 = buffer.data(lsi + 854);
    const auto *lsi_861 = buffer.data(lsi + 861);
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
    const auto *lsi_895 = buffer.data(lsi + 895);
    const auto *lsi_896 = buffer.data(lsi + 896);
    const auto *lsi_898 = buffer.data(lsi + 898);
    const auto *lsi_899 = buffer.data(lsi + 899);
    const auto *lsi_901 = buffer.data(lsi + 901);
    const auto *lsi_902 = buffer.data(lsi + 902);
    const auto *lsi_905 = buffer.data(lsi + 905);
    const auto *lsi_906 = buffer.data(lsi + 906);
    const auto *lsi_910 = buffer.data(lsi + 910);
    const auto *lsi_923 = buffer.data(lsi + 923);
    const auto *lsi_924 = buffer.data(lsi + 924);
    const auto *lsi_926 = buffer.data(lsi + 926);
    const auto *lsi_929 = buffer.data(lsi + 929);
    const auto *lsi_933 = buffer.data(lsi + 933);
    const auto *lsi_938 = buffer.data(lsi + 938);
    const auto *lsi_1078 = buffer.data(lsi + 1078);
    const auto *lsi_1079 = buffer.data(lsi + 1079);
    const auto *lsi_1081 = buffer.data(lsi + 1081);
    const auto *lsi_1082 = buffer.data(lsi + 1082);
    const auto *lsi_1084 = buffer.data(lsi + 1084);
    const auto *lsi_1085 = buffer.data(lsi + 1085);
    const auto *lsi_1086 = buffer.data(lsi + 1086);
    const auto *lsi_1087 = buffer.data(lsi + 1087);
    const auto *lsi_1088 = buffer.data(lsi + 1088);
    const auto *lsi_1089 = buffer.data(lsi + 1089);
    const auto *lsi_1090 = buffer.data(lsi + 1090);
    const auto *lsi_1091 = buffer.data(lsi + 1091);
    const auto *lsi_1092 = buffer.data(lsi + 1092);
    const auto *lsi_1095 = buffer.data(lsi + 1095);
    const auto *lsi_1097 = buffer.data(lsi + 1097);
    const auto *lsi_1098 = buffer.data(lsi + 1098);
    const auto *lsi_1101 = buffer.data(lsi + 1101);
    const auto *lsi_1102 = buffer.data(lsi + 1102);
    const auto *lsi_1104 = buffer.data(lsi + 1104);
    const auto *lsi_1106 = buffer.data(lsi + 1106);
    const auto *lsi_1107 = buffer.data(lsi + 1107);
    const auto *lsi_1109 = buffer.data(lsi + 1109);
    const auto *lsi_1110 = buffer.data(lsi + 1110);
    const auto *lsi_1112 = buffer.data(lsi + 1112);
    const auto *lsi_1113 = buffer.data(lsi + 1113);
    const auto *lsi_1114 = buffer.data(lsi + 1114);
    const auto *lsi_1115 = buffer.data(lsi + 1115);
    const auto *lsi_1116 = buffer.data(lsi + 1116);
    const auto *lsi_1117 = buffer.data(lsi + 1117);
    const auto *lsi_1118 = buffer.data(lsi + 1118);
    const auto *lsi_1119 = buffer.data(lsi + 1119);
    const auto *lsi_1120 = buffer.data(lsi + 1120);
    const auto *lsi_1123 = buffer.data(lsi + 1123);
    const auto *lsi_1125 = buffer.data(lsi + 1125);
    const auto *lsi_1126 = buffer.data(lsi + 1126);
    const auto *lsi_1129 = buffer.data(lsi + 1129);
    const auto *lsi_1130 = buffer.data(lsi + 1130);
    const auto *lsi_1132 = buffer.data(lsi + 1132);
    const auto *lsi_1134 = buffer.data(lsi + 1134);
    const auto *lsi_1135 = buffer.data(lsi + 1135);
    const auto *lsi_1137 = buffer.data(lsi + 1137);
    const auto *lsi_1138 = buffer.data(lsi + 1138);
    const auto *lsi_1140 = buffer.data(lsi + 1140);
    const auto *lsi_1141 = buffer.data(lsi + 1141);
    const auto *lsi_1142 = buffer.data(lsi + 1142);
    const auto *lsi_1143 = buffer.data(lsi + 1143);
    const auto *lsi_1144 = buffer.data(lsi + 1144);
    const auto *lsi_1145 = buffer.data(lsi + 1145);
    const auto *lsi_1146 = buffer.data(lsi + 1146);
    const auto *lsi_1147 = buffer.data(lsi + 1147);
    const auto *lsi_1148 = buffer.data(lsi + 1148);
    const auto *lsi_1151 = buffer.data(lsi + 1151);
    const auto *lsi_1153 = buffer.data(lsi + 1153);
    const auto *lsi_1154 = buffer.data(lsi + 1154);
    const auto *lsi_1157 = buffer.data(lsi + 1157);
    const auto *lsi_1158 = buffer.data(lsi + 1158);
    const auto *lsi_1160 = buffer.data(lsi + 1160);
    const auto *lsi_1162 = buffer.data(lsi + 1162);
    const auto *lsi_1163 = buffer.data(lsi + 1163);
    const auto *lsi_1165 = buffer.data(lsi + 1165);
    const auto *lsi_1166 = buffer.data(lsi + 1166);
    const auto *lsi_1168 = buffer.data(lsi + 1168);

    const auto *lsk1_1382 = buffer.data(lsk1 + 1382);
    const auto *lsk1_1383 = buffer.data(lsk1 + 1383);
    const auto *lsk1_1385 = buffer.data(lsk1 + 1385);
    const auto *lsk1_1386 = buffer.data(lsk1 + 1386);
    const auto *lsk1_1388 = buffer.data(lsk1 + 1388);
    const auto *lsk1_1396 = buffer.data(lsk1 + 1396);
    const auto *lsk1_1398 = buffer.data(lsk1 + 1398);
    const auto *lsk1_1399 = buffer.data(lsk1 + 1399);
    const auto *lsk1_1400 = buffer.data(lsk1 + 1400);
    const auto *lsk1_1401 = buffer.data(lsk1 + 1401);
    const auto *lsk1_1403 = buffer.data(lsk1 + 1403);
    const auto *lsk1_1404 = buffer.data(lsk1 + 1404);
    const auto *lsk1_1407 = buffer.data(lsk1 + 1407);
    const auto *lsk1_1409 = buffer.data(lsk1 + 1409);
    const auto *lsk1_1410 = buffer.data(lsk1 + 1410);
    const auto *lsk1_1413 = buffer.data(lsk1 + 1413);
    const auto *lsk1_1414 = buffer.data(lsk1 + 1414);
    const auto *lsk1_1416 = buffer.data(lsk1 + 1416);
    const auto *lsk1_1418 = buffer.data(lsk1 + 1418);
    const auto *lsk1_1419 = buffer.data(lsk1 + 1419);
    const auto *lsk1_1421 = buffer.data(lsk1 + 1421);
    const auto *lsk1_1422 = buffer.data(lsk1 + 1422);
    const auto *lsk1_1424 = buffer.data(lsk1 + 1424);
    const auto *lsk1_1432 = buffer.data(lsk1 + 1432);
    const auto *lsk1_1434 = buffer.data(lsk1 + 1434);
    const auto *lsk1_1435 = buffer.data(lsk1 + 1435);
    const auto *lsk1_1436 = buffer.data(lsk1 + 1436);
    const auto *lsk1_1437 = buffer.data(lsk1 + 1437);
    const auto *lsk1_1439 = buffer.data(lsk1 + 1439);
    const auto *lsk1_1440 = buffer.data(lsk1 + 1440);
    const auto *lsk1_1443 = buffer.data(lsk1 + 1443);
    const auto *lsk1_1445 = buffer.data(lsk1 + 1445);
    const auto *lsk1_1446 = buffer.data(lsk1 + 1446);
    const auto *lsk1_1449 = buffer.data(lsk1 + 1449);
    const auto *lsk1_1450 = buffer.data(lsk1 + 1450);
    const auto *lsk1_1452 = buffer.data(lsk1 + 1452);
    const auto *lsk1_1454 = buffer.data(lsk1 + 1454);
    const auto *lsk1_1455 = buffer.data(lsk1 + 1455);
    const auto *lsk1_1457 = buffer.data(lsk1 + 1457);
    const auto *lsk1_1458 = buffer.data(lsk1 + 1458);
    const auto *lsk1_1460 = buffer.data(lsk1 + 1460);
    const auto *lsk1_1468 = buffer.data(lsk1 + 1468);
    const auto *lsk1_1470 = buffer.data(lsk1 + 1470);
    const auto *lsk1_1471 = buffer.data(lsk1 + 1471);
    const auto *lsk1_1472 = buffer.data(lsk1 + 1472);
    const auto *lsk1_1473 = buffer.data(lsk1 + 1473);
    const auto *lsk1_1475 = buffer.data(lsk1 + 1475);
    const auto *lsk1_1476 = buffer.data(lsk1 + 1476);
    const auto *lsk1_1479 = buffer.data(lsk1 + 1479);
    const auto *lsk1_1481 = buffer.data(lsk1 + 1481);
    const auto *lsk1_1482 = buffer.data(lsk1 + 1482);
    const auto *lsk1_1485 = buffer.data(lsk1 + 1485);
    const auto *lsk1_1486 = buffer.data(lsk1 + 1486);
    const auto *lsk1_1488 = buffer.data(lsk1 + 1488);
    const auto *lsk1_1490 = buffer.data(lsk1 + 1490);
    const auto *lsk1_1491 = buffer.data(lsk1 + 1491);
    const auto *lsk1_1493 = buffer.data(lsk1 + 1493);
    const auto *lsk1_1494 = buffer.data(lsk1 + 1494);
    const auto *lsk1_1496 = buffer.data(lsk1 + 1496);

    const auto *msi_1074 = buffer.data(msi + 1074);
    const auto *msi_1078 = buffer.data(msi + 1078);
    const auto *msi_1085 = buffer.data(msi + 1085);
    const auto *msi_1086 = buffer.data(msi + 1086);
    const auto *msi_1087 = buffer.data(msi + 1087);
    const auto *msi_1088 = buffer.data(msi + 1088);
    const auto *msi_1089 = buffer.data(msi + 1089);
    const auto *msi_1090 = buffer.data(msi + 1090);
    const auto *msi_1091 = buffer.data(msi + 1091);
    const auto *msi_1092 = buffer.data(msi + 1092);
    const auto *msi_1094 = buffer.data(msi + 1094);
    const auto *msi_1095 = buffer.data(msi + 1095);
    const auto *msi_1097 = buffer.data(msi + 1097);
    const auto *msi_1098 = buffer.data(msi + 1098);
    const auto *msi_1101 = buffer.data(msi + 1101);
    const auto *msi_1102 = buffer.data(msi + 1102);
    const auto *msi_1106 = buffer.data(msi + 1106);
    const auto *msi_1113 = buffer.data(msi + 1113);
    const auto *msi_1114 = buffer.data(msi + 1114);
    const auto *msi_1115 = buffer.data(msi + 1115);
    const auto *msi_1116 = buffer.data(msi + 1116);
    const auto *msi_1117 = buffer.data(msi + 1117);
    const auto *msi_1118 = buffer.data(msi + 1118);
    const auto *msi_1119 = buffer.data(msi + 1119);
    const auto *msi_1120 = buffer.data(msi + 1120);
    const auto *msi_1122 = buffer.data(msi + 1122);
    const auto *msi_1123 = buffer.data(msi + 1123);
    const auto *msi_1125 = buffer.data(msi + 1125);
    const auto *msi_1126 = buffer.data(msi + 1126);
    const auto *msi_1129 = buffer.data(msi + 1129);
    const auto *msi_1130 = buffer.data(msi + 1130);
    const auto *msi_1134 = buffer.data(msi + 1134);
    const auto *msi_1141 = buffer.data(msi + 1141);
    const auto *msi_1142 = buffer.data(msi + 1142);
    const auto *msi_1143 = buffer.data(msi + 1143);
    const auto *msi_1144 = buffer.data(msi + 1144);
    const auto *msi_1145 = buffer.data(msi + 1145);
    const auto *msi_1146 = buffer.data(msi + 1146);
    const auto *msi_1147 = buffer.data(msi + 1147);
    const auto *msi_1148 = buffer.data(msi + 1148);
    const auto *msi_1150 = buffer.data(msi + 1150);
    const auto *msi_1151 = buffer.data(msi + 1151);
    const auto *msi_1153 = buffer.data(msi + 1153);
    const auto *msi_1154 = buffer.data(msi + 1154);
    const auto *msi_1157 = buffer.data(msi + 1157);
    const auto *msi_1158 = buffer.data(msi + 1158);
    const auto *msi_1162 = buffer.data(msi + 1162);

#pragma omp simd aligned(t_1382, t_1383, t_1384, pa_x, pc_x, pc_z, lsk0_1382, lsk0_1383, \
                         lsi_822, lsi_1078, lsi_1079, lsk1_1382, lsk1_1383, \
                         msi_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1382[k] = pa_x[k] * lsk0_1382[k]
                    + f_15 * lsi_1078[k]
                    - f_12 * pc_x[k] * lsk1_1382[k];

        t_1383[k] = pa_x[k] * lsk0_1383[k]
                    + f_14 * lsi_1079[k]
                    - f_12 * pc_x[k] * lsk1_1383[k];

        t_1384[k] = f_14 * lsi_822[k]
                    + f_3 * pc_z[k] * msi_1074[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, pa_x, pc_x, pc_y, lsk0_1385, lsk0_1386, \
                         lsi_854, lsi_1081, lsi_1082, lsk1_1385, lsk1_1386, \
                         msi_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = pa_x[k] * lsk0_1385[k]
                    + f_14 * lsi_1081[k]
                    - f_12 * pc_x[k] * lsk1_1385[k];

        t_1386[k] = pa_x[k] * lsk0_1386[k]
                    + f_14 * lsi_1082[k]
                    - f_12 * pc_x[k] * lsk1_1386[k];

        t_1387[k] = f_22 * lsi_854[k]
                    + f_3 * pc_y[k] * msi_1078[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, t_1391, pa_x, pc_x, lsk0_1388, lsi_1084, \
                         lsi_1085, lsi_1086, lsi_1087, lsk1_1388, msi_1085, msi_1086, \
                         msi_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = pa_x[k] * lsk0_1388[k]
                    + f_14 * lsi_1084[k]
                    - f_12 * pc_x[k] * lsk1_1388[k];

        t_1389[k] = f_13 * lsi_1085[k]
                    + f_3 * pc_x[k] * msi_1085[k];

        t_1390[k] = f_13 * lsi_1086[k]
                    + f_3 * pc_x[k] * msi_1086[k];

        t_1391[k] = f_13 * lsi_1087[k]
                    + f_3 * pc_x[k] * msi_1087[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, t_1395, pc_x, lsi_1088, lsi_1089, lsi_1090, \
                         lsi_1091, msi_1088, msi_1089, msi_1090, \
                         msi_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_13 * lsi_1088[k]
                    + f_3 * pc_x[k] * msi_1088[k];

        t_1393[k] = f_13 * lsi_1089[k]
                    + f_3 * pc_x[k] * msi_1089[k];

        t_1394[k] = f_13 * lsi_1090[k]
                    + f_3 * pc_x[k] * msi_1090[k];

        t_1395[k] = f_13 * lsi_1091[k]
                    + f_3 * pc_x[k] * msi_1091[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, t_1399, pa_x, pc_x, pc_z, lsk0_1396, \
                         lsk0_1398, lsk0_1399, lsi_833, lsk1_1396, lsk1_1398, lsk1_1399, \
                         msi_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = pa_x[k] * lsk0_1396[k]
                    - f_12 * pc_x[k] * lsk1_1396[k];

        t_1397[k] = f_14 * lsi_833[k]
                    + f_3 * pc_z[k] * msi_1085[k];

        t_1398[k] = pa_x[k] * lsk0_1398[k]
                    - f_12 * pc_x[k] * lsk1_1398[k];

        t_1399[k] = pa_x[k] * lsk0_1399[k]
                    - f_12 * pc_x[k] * lsk1_1399[k];
    }

#pragma omp simd aligned(t_1400, t_1401, t_1402, t_1403, pa_x, pc_x, pc_y, lsk0_1400, \
                         lsk0_1401, lsk0_1403, lsi_867, lsk1_1400, lsk1_1401, lsk1_1403, \
                         msi_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1400[k] = pa_x[k] * lsk0_1400[k]
                    - f_12 * pc_x[k] * lsk1_1400[k];

        t_1401[k] = pa_x[k] * lsk0_1401[k]
                    - f_12 * pc_x[k] * lsk1_1401[k];

        t_1402[k] = f_22 * lsi_867[k]
                    + f_3 * pc_y[k] * msi_1091[k];

        t_1403[k] = pa_x[k] * lsk0_1403[k]
                    - f_12 * pc_x[k] * lsk1_1403[k];
    }

#pragma omp simd aligned(t_1404, t_1405, t_1406, pa_x, pc_x, pc_y, pc_z, lsk0_1404, lsi_840, \
                         lsi_868, lsi_1092, lsk1_1404, msi_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1404[k] = pa_x[k] * lsk0_1404[k]
                    + f_21 * lsi_1092[k]
                    - f_12 * pc_x[k] * lsk1_1404[k];

        t_1405[k] = f_17 * lsi_868[k]
                    + f_3 * pc_y[k] * msi_1092[k];

        t_1406[k] = f_15 * lsi_840[k]
                    + f_3 * pc_z[k] * msi_1092[k];
    }

#pragma omp simd aligned(t_1407, t_1408, t_1409, pa_x, pc_x, pc_y, lsk0_1407, lsk0_1409, \
                         lsi_870, lsi_1095, lsi_1097, lsk1_1407, lsk1_1409, \
                         msi_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1407[k] = pa_x[k] * lsk0_1407[k]
                    + f_17 * lsi_1095[k]
                    - f_12 * pc_x[k] * lsk1_1407[k];

        t_1408[k] = f_17 * lsi_870[k]
                    + f_3 * pc_y[k] * msi_1094[k];

        t_1409[k] = pa_x[k] * lsk0_1409[k]
                    + f_17 * lsi_1097[k]
                    - f_12 * pc_x[k] * lsk1_1409[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pa_x, pc_x, pc_y, pc_z, lsk0_1410, lsi_843, \
                         lsi_873, lsi_1098, lsk1_1410, msi_1095, \
                         msi_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = pa_x[k] * lsk0_1410[k]
                    + f_16 * lsi_1098[k]
                    - f_12 * pc_x[k] * lsk1_1410[k];

        t_1411[k] = f_15 * lsi_843[k]
                    + f_3 * pc_z[k] * msi_1095[k];

        t_1412[k] = f_17 * lsi_873[k]
                    + f_3 * pc_y[k] * msi_1097[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, pa_x, pc_x, pc_z, lsk0_1413, lsk0_1414, \
                         lsi_846, lsi_1101, lsi_1102, lsk1_1413, lsk1_1414, \
                         msi_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = pa_x[k] * lsk0_1413[k]
                    + f_16 * lsi_1101[k]
                    - f_12 * pc_x[k] * lsk1_1413[k];

        t_1414[k] = pa_x[k] * lsk0_1414[k]
                    + f_15 * lsi_1102[k]
                    - f_12 * pc_x[k] * lsk1_1414[k];

        t_1415[k] = f_15 * lsi_846[k]
                    + f_3 * pc_z[k] * msi_1098[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, pa_x, pc_x, pc_y, lsk0_1416, lsk0_1418, \
                         lsi_877, lsi_1104, lsi_1106, lsk1_1416, lsk1_1418, \
                         msi_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = pa_x[k] * lsk0_1416[k]
                    + f_15 * lsi_1104[k]
                    - f_12 * pc_x[k] * lsk1_1416[k];

        t_1417[k] = f_17 * lsi_877[k]
                    + f_3 * pc_y[k] * msi_1101[k];

        t_1418[k] = pa_x[k] * lsk0_1418[k]
                    + f_15 * lsi_1106[k]
                    - f_12 * pc_x[k] * lsk1_1418[k];
    }

#pragma omp simd aligned(t_1419, t_1420, t_1421, pa_x, pc_x, pc_z, lsk0_1419, lsk0_1421, \
                         lsi_850, lsi_1107, lsi_1109, lsk1_1419, lsk1_1421, \
                         msi_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1419[k] = pa_x[k] * lsk0_1419[k]
                    + f_14 * lsi_1107[k]
                    - f_12 * pc_x[k] * lsk1_1419[k];

        t_1420[k] = f_15 * lsi_850[k]
                    + f_3 * pc_z[k] * msi_1102[k];

        t_1421[k] = pa_x[k] * lsk0_1421[k]
                    + f_14 * lsi_1109[k]
                    - f_12 * pc_x[k] * lsk1_1421[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, pa_x, pc_x, pc_y, lsk0_1422, lsk0_1424, \
                         lsi_882, lsi_1110, lsi_1112, lsk1_1422, lsk1_1424, \
                         msi_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = pa_x[k] * lsk0_1422[k]
                    + f_14 * lsi_1110[k]
                    - f_12 * pc_x[k] * lsk1_1422[k];

        t_1423[k] = f_17 * lsi_882[k]
                    + f_3 * pc_y[k] * msi_1106[k];

        t_1424[k] = pa_x[k] * lsk0_1424[k]
                    + f_14 * lsi_1112[k]
                    - f_12 * pc_x[k] * lsk1_1424[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, t_1428, t_1429, pc_x, lsi_1113, lsi_1114, \
                         lsi_1115, lsi_1116, lsi_1117, msi_1113, msi_1114, msi_1115, msi_1116, \
                         msi_1117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_13 * lsi_1113[k]
                    + f_3 * pc_x[k] * msi_1113[k];

        t_1426[k] = f_13 * lsi_1114[k]
                    + f_3 * pc_x[k] * msi_1114[k];

        t_1427[k] = f_13 * lsi_1115[k]
                    + f_3 * pc_x[k] * msi_1115[k];

        t_1428[k] = f_13 * lsi_1116[k]
                    + f_3 * pc_x[k] * msi_1116[k];

        t_1429[k] = f_13 * lsi_1117[k]
                    + f_3 * pc_x[k] * msi_1117[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, t_1433, pa_x, pc_x, pc_z, lsk0_1432, lsi_861, \
                         lsi_1118, lsi_1119, lsk1_1432, msi_1113, msi_1118, \
                         msi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_13 * lsi_1118[k]
                    + f_3 * pc_x[k] * msi_1118[k];

        t_1431[k] = f_13 * lsi_1119[k]
                    + f_3 * pc_x[k] * msi_1119[k];

        t_1432[k] = pa_x[k] * lsk0_1432[k]
                    - f_12 * pc_x[k] * lsk1_1432[k];

        t_1433[k] = f_15 * lsi_861[k]
                    + f_3 * pc_z[k] * msi_1113[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, t_1437, pa_x, pc_x, lsk0_1434, lsk0_1435, \
                         lsk0_1436, lsk0_1437, lsk1_1434, lsk1_1435, lsk1_1436, \
                         lsk1_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = pa_x[k] * lsk0_1434[k]
                    - f_12 * pc_x[k] * lsk1_1434[k];

        t_1435[k] = pa_x[k] * lsk0_1435[k]
                    - f_12 * pc_x[k] * lsk1_1435[k];

        t_1436[k] = pa_x[k] * lsk0_1436[k]
                    - f_12 * pc_x[k] * lsk1_1436[k];

        t_1437[k] = pa_x[k] * lsk0_1437[k]
                    - f_12 * pc_x[k] * lsk1_1437[k];
    }

#pragma omp simd aligned(t_1438, t_1439, t_1440, t_1441, pa_x, pc_x, pc_y, lsk0_1439, \
                         lsk0_1440, lsi_895, lsi_896, lsi_1120, lsk1_1439, lsk1_1440, \
                         msi_1119, msi_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1438[k] = f_17 * lsi_895[k]
                    + f_3 * pc_y[k] * msi_1119[k];

        t_1439[k] = pa_x[k] * lsk0_1439[k]
                    - f_12 * pc_x[k] * lsk1_1439[k];

        t_1440[k] = pa_x[k] * lsk0_1440[k]
                    + f_21 * lsi_1120[k]
                    - f_12 * pc_x[k] * lsk1_1440[k];

        t_1441[k] = f_16 * lsi_896[k]
                    + f_3 * pc_y[k] * msi_1120[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pa_x, pc_x, pc_y, pc_z, lsk0_1443, lsi_868, \
                         lsi_898, lsi_1123, lsk1_1443, msi_1120, \
                         msi_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_16 * lsi_868[k]
                    + f_3 * pc_z[k] * msi_1120[k];

        t_1443[k] = pa_x[k] * lsk0_1443[k]
                    + f_17 * lsi_1123[k]
                    - f_12 * pc_x[k] * lsk1_1443[k];

        t_1444[k] = f_16 * lsi_898[k]
                    + f_3 * pc_y[k] * msi_1122[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pa_x, pc_x, pc_z, lsk0_1445, lsk0_1446, \
                         lsi_871, lsi_1125, lsi_1126, lsk1_1445, lsk1_1446, \
                         msi_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = pa_x[k] * lsk0_1445[k]
                    + f_17 * lsi_1125[k]
                    - f_12 * pc_x[k] * lsk1_1445[k];

        t_1446[k] = pa_x[k] * lsk0_1446[k]
                    + f_16 * lsi_1126[k]
                    - f_12 * pc_x[k] * lsk1_1446[k];

        t_1447[k] = f_16 * lsi_871[k]
                    + f_3 * pc_z[k] * msi_1123[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pa_x, pc_x, pc_y, lsk0_1449, lsk0_1450, \
                         lsi_901, lsi_1129, lsi_1130, lsk1_1449, lsk1_1450, \
                         msi_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_16 * lsi_901[k]
                    + f_3 * pc_y[k] * msi_1125[k];

        t_1449[k] = pa_x[k] * lsk0_1449[k]
                    + f_16 * lsi_1129[k]
                    - f_12 * pc_x[k] * lsk1_1449[k];

        t_1450[k] = pa_x[k] * lsk0_1450[k]
                    + f_15 * lsi_1130[k]
                    - f_12 * pc_x[k] * lsk1_1450[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pa_x, pc_x, pc_y, pc_z, lsk0_1452, lsi_874, \
                         lsi_905, lsi_1132, lsk1_1452, msi_1126, \
                         msi_1129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_16 * lsi_874[k]
                    + f_3 * pc_z[k] * msi_1126[k];

        t_1452[k] = pa_x[k] * lsk0_1452[k]
                    + f_15 * lsi_1132[k]
                    - f_12 * pc_x[k] * lsk1_1452[k];

        t_1453[k] = f_16 * lsi_905[k]
                    + f_3 * pc_y[k] * msi_1129[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pa_x, pc_x, pc_z, lsk0_1454, lsk0_1455, \
                         lsi_878, lsi_1134, lsi_1135, lsk1_1454, lsk1_1455, \
                         msi_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = pa_x[k] * lsk0_1454[k]
                    + f_15 * lsi_1134[k]
                    - f_12 * pc_x[k] * lsk1_1454[k];

        t_1455[k] = pa_x[k] * lsk0_1455[k]
                    + f_14 * lsi_1135[k]
                    - f_12 * pc_x[k] * lsk1_1455[k];

        t_1456[k] = f_16 * lsi_878[k]
                    + f_3 * pc_z[k] * msi_1130[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pa_x, pc_x, pc_y, lsk0_1457, lsk0_1458, \
                         lsi_910, lsi_1137, lsi_1138, lsk1_1457, lsk1_1458, \
                         msi_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = pa_x[k] * lsk0_1457[k]
                    + f_14 * lsi_1137[k]
                    - f_12 * pc_x[k] * lsk1_1457[k];

        t_1458[k] = pa_x[k] * lsk0_1458[k]
                    + f_14 * lsi_1138[k]
                    - f_12 * pc_x[k] * lsk1_1458[k];

        t_1459[k] = f_16 * lsi_910[k]
                    + f_3 * pc_y[k] * msi_1134[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, t_1463, pa_x, pc_x, lsk0_1460, lsi_1140, \
                         lsi_1141, lsi_1142, lsi_1143, lsk1_1460, msi_1141, msi_1142, \
                         msi_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = pa_x[k] * lsk0_1460[k]
                    + f_14 * lsi_1140[k]
                    - f_12 * pc_x[k] * lsk1_1460[k];

        t_1461[k] = f_13 * lsi_1141[k]
                    + f_3 * pc_x[k] * msi_1141[k];

        t_1462[k] = f_13 * lsi_1142[k]
                    + f_3 * pc_x[k] * msi_1142[k];

        t_1463[k] = f_13 * lsi_1143[k]
                    + f_3 * pc_x[k] * msi_1143[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, t_1467, pc_x, lsi_1144, lsi_1145, lsi_1146, \
                         lsi_1147, msi_1144, msi_1145, msi_1146, \
                         msi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_13 * lsi_1144[k]
                    + f_3 * pc_x[k] * msi_1144[k];

        t_1465[k] = f_13 * lsi_1145[k]
                    + f_3 * pc_x[k] * msi_1145[k];

        t_1466[k] = f_13 * lsi_1146[k]
                    + f_3 * pc_x[k] * msi_1146[k];

        t_1467[k] = f_13 * lsi_1147[k]
                    + f_3 * pc_x[k] * msi_1147[k];
    }

#pragma omp simd aligned(t_1468, t_1469, t_1470, t_1471, pa_x, pc_x, pc_z, lsk0_1468, \
                         lsk0_1470, lsk0_1471, lsi_889, lsk1_1468, lsk1_1470, lsk1_1471, \
                         msi_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1468[k] = pa_x[k] * lsk0_1468[k]
                    - f_12 * pc_x[k] * lsk1_1468[k];

        t_1469[k] = f_16 * lsi_889[k]
                    + f_3 * pc_z[k] * msi_1141[k];

        t_1470[k] = pa_x[k] * lsk0_1470[k]
                    - f_12 * pc_x[k] * lsk1_1470[k];

        t_1471[k] = pa_x[k] * lsk0_1471[k]
                    - f_12 * pc_x[k] * lsk1_1471[k];
    }

#pragma omp simd aligned(t_1472, t_1473, t_1474, t_1475, pa_x, pc_x, pc_y, lsk0_1472, \
                         lsk0_1473, lsk0_1475, lsi_923, lsk1_1472, lsk1_1473, lsk1_1475, \
                         msi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1472[k] = pa_x[k] * lsk0_1472[k]
                    - f_12 * pc_x[k] * lsk1_1472[k];

        t_1473[k] = pa_x[k] * lsk0_1473[k]
                    - f_12 * pc_x[k] * lsk1_1473[k];

        t_1474[k] = f_16 * lsi_923[k]
                    + f_3 * pc_y[k] * msi_1147[k];

        t_1475[k] = pa_x[k] * lsk0_1475[k]
                    - f_12 * pc_x[k] * lsk1_1475[k];
    }

#pragma omp simd aligned(t_1476, t_1477, t_1478, pa_x, pc_x, pc_y, pc_z, lsk0_1476, lsi_896, \
                         lsi_924, lsi_1148, lsk1_1476, msi_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1476[k] = pa_x[k] * lsk0_1476[k]
                    + f_21 * lsi_1148[k]
                    - f_12 * pc_x[k] * lsk1_1476[k];

        t_1477[k] = f_15 * lsi_924[k]
                    + f_3 * pc_y[k] * msi_1148[k];

        t_1478[k] = f_17 * lsi_896[k]
                    + f_3 * pc_z[k] * msi_1148[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, pa_x, pc_x, pc_y, lsk0_1479, lsk0_1481, \
                         lsi_926, lsi_1151, lsi_1153, lsk1_1479, lsk1_1481, \
                         msi_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = pa_x[k] * lsk0_1479[k]
                    + f_17 * lsi_1151[k]
                    - f_12 * pc_x[k] * lsk1_1479[k];

        t_1480[k] = f_15 * lsi_926[k]
                    + f_3 * pc_y[k] * msi_1150[k];

        t_1481[k] = pa_x[k] * lsk0_1481[k]
                    + f_17 * lsi_1153[k]
                    - f_12 * pc_x[k] * lsk1_1481[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, pa_x, pc_x, pc_y, pc_z, lsk0_1482, lsi_899, \
                         lsi_929, lsi_1154, lsk1_1482, msi_1151, \
                         msi_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = pa_x[k] * lsk0_1482[k]
                    + f_16 * lsi_1154[k]
                    - f_12 * pc_x[k] * lsk1_1482[k];

        t_1483[k] = f_17 * lsi_899[k]
                    + f_3 * pc_z[k] * msi_1151[k];

        t_1484[k] = f_15 * lsi_929[k]
                    + f_3 * pc_y[k] * msi_1153[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, pa_x, pc_x, pc_z, lsk0_1485, lsk0_1486, \
                         lsi_902, lsi_1157, lsi_1158, lsk1_1485, lsk1_1486, \
                         msi_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = pa_x[k] * lsk0_1485[k]
                    + f_16 * lsi_1157[k]
                    - f_12 * pc_x[k] * lsk1_1485[k];

        t_1486[k] = pa_x[k] * lsk0_1486[k]
                    + f_15 * lsi_1158[k]
                    - f_12 * pc_x[k] * lsk1_1486[k];

        t_1487[k] = f_17 * lsi_902[k]
                    + f_3 * pc_z[k] * msi_1154[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pa_x, pc_x, pc_y, lsk0_1488, lsk0_1490, \
                         lsi_933, lsi_1160, lsi_1162, lsk1_1488, lsk1_1490, \
                         msi_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = pa_x[k] * lsk0_1488[k]
                    + f_15 * lsi_1160[k]
                    - f_12 * pc_x[k] * lsk1_1488[k];

        t_1489[k] = f_15 * lsi_933[k]
                    + f_3 * pc_y[k] * msi_1157[k];

        t_1490[k] = pa_x[k] * lsk0_1490[k]
                    + f_15 * lsi_1162[k]
                    - f_12 * pc_x[k] * lsk1_1490[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pa_x, pc_x, pc_z, lsk0_1491, lsk0_1493, \
                         lsi_906, lsi_1163, lsi_1165, lsk1_1491, lsk1_1493, \
                         msi_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = pa_x[k] * lsk0_1491[k]
                    + f_14 * lsi_1163[k]
                    - f_12 * pc_x[k] * lsk1_1491[k];

        t_1492[k] = f_17 * lsi_906[k]
                    + f_3 * pc_z[k] * msi_1158[k];

        t_1493[k] = pa_x[k] * lsk0_1493[k]
                    + f_14 * lsi_1165[k]
                    - f_12 * pc_x[k] * lsk1_1493[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, pa_x, pc_x, pc_y, lsk0_1494, lsk0_1496, \
                         lsi_938, lsi_1166, lsi_1168, lsk1_1494, lsk1_1496, \
                         msi_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = pa_x[k] * lsk0_1494[k]
                    + f_14 * lsi_1166[k]
                    - f_12 * pc_x[k] * lsk1_1494[k];

        t_1495[k] = f_15 * lsi_938[k]
                    + f_3 * pc_y[k] * msi_1162[k];

        t_1496[k] = pa_x[k] * lsk0_1496[k]
                    + f_14 * lsi_1168[k]
                    - f_12 * pc_x[k] * lsk1_1496[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsk0,
                                                           const size_t lsi, const size_t lsk1,
                                                           const size_t msh0, const size_t msh1,
                                                           const size_t msi, const size_t ncols,
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
    const auto f_18 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_1260 = buffer.data(lsk0 + 1260);
    const auto *lsk0_1265 = buffer.data(lsk0 + 1265);
    const auto *lsk0_1269 = buffer.data(lsk0 + 1269);
    const auto *lsk0_1274 = buffer.data(lsk0 + 1274);
    const auto *lsk0_1280 = buffer.data(lsk0 + 1280);
    const auto *lsk0_1504 = buffer.data(lsk0 + 1504);
    const auto *lsk0_1506 = buffer.data(lsk0 + 1506);
    const auto *lsk0_1507 = buffer.data(lsk0 + 1507);
    const auto *lsk0_1508 = buffer.data(lsk0 + 1508);
    const auto *lsk0_1509 = buffer.data(lsk0 + 1509);
    const auto *lsk0_1511 = buffer.data(lsk0 + 1511);
    const auto *lsk0_1512 = buffer.data(lsk0 + 1512);
    const auto *lsk0_1515 = buffer.data(lsk0 + 1515);
    const auto *lsk0_1517 = buffer.data(lsk0 + 1517);
    const auto *lsk0_1518 = buffer.data(lsk0 + 1518);
    const auto *lsk0_1521 = buffer.data(lsk0 + 1521);
    const auto *lsk0_1522 = buffer.data(lsk0 + 1522);
    const auto *lsk0_1524 = buffer.data(lsk0 + 1524);
    const auto *lsk0_1526 = buffer.data(lsk0 + 1526);
    const auto *lsk0_1527 = buffer.data(lsk0 + 1527);
    const auto *lsk0_1529 = buffer.data(lsk0 + 1529);
    const auto *lsk0_1530 = buffer.data(lsk0 + 1530);
    const auto *lsk0_1532 = buffer.data(lsk0 + 1532);
    const auto *lsk0_1540 = buffer.data(lsk0 + 1540);
    const auto *lsk0_1542 = buffer.data(lsk0 + 1542);
    const auto *lsk0_1543 = buffer.data(lsk0 + 1543);
    const auto *lsk0_1544 = buffer.data(lsk0 + 1544);
    const auto *lsk0_1545 = buffer.data(lsk0 + 1545);
    const auto *lsk0_1547 = buffer.data(lsk0 + 1547);
    const auto *lsk0_1551 = buffer.data(lsk0 + 1551);
    const auto *lsk0_1554 = buffer.data(lsk0 + 1554);
    const auto *lsk0_1558 = buffer.data(lsk0 + 1558);
    const auto *lsk0_1560 = buffer.data(lsk0 + 1560);
    const auto *lsk0_1563 = buffer.data(lsk0 + 1563);
    const auto *lsk0_1565 = buffer.data(lsk0 + 1565);
    const auto *lsk0_1566 = buffer.data(lsk0 + 1566);
    const auto *lsk0_1576 = buffer.data(lsk0 + 1576);
    const auto *lsk0_1578 = buffer.data(lsk0 + 1578);
    const auto *lsk0_1579 = buffer.data(lsk0 + 1579);
    const auto *lsk0_1580 = buffer.data(lsk0 + 1580);
    const auto *lsk0_1581 = buffer.data(lsk0 + 1581);
    const auto *lsk0_1583 = buffer.data(lsk0 + 1583);
    const auto *lsk0_1584 = buffer.data(lsk0 + 1584);
    const auto *lsk0_1589 = buffer.data(lsk0 + 1589);
    const auto *lsk0_1593 = buffer.data(lsk0 + 1593);
    const auto *lsk0_1598 = buffer.data(lsk0 + 1598);
    const auto *lsk0_1604 = buffer.data(lsk0 + 1604);
    const auto *lsk0_1612 = buffer.data(lsk0 + 1612);
    const auto *lsk0_1613 = buffer.data(lsk0 + 1613);
    const auto *lsk0_1614 = buffer.data(lsk0 + 1614);
    const auto *lsk0_1615 = buffer.data(lsk0 + 1615);
    const auto *lsk0_1616 = buffer.data(lsk0 + 1616);
    const auto *lsk0_1617 = buffer.data(lsk0 + 1617);
    const auto *lsk0_1619 = buffer.data(lsk0 + 1619);

    const auto *lsi_917 = buffer.data(lsi + 917);
    const auto *lsi_924 = buffer.data(lsi + 924);
    const auto *lsi_927 = buffer.data(lsi + 927);
    const auto *lsi_930 = buffer.data(lsi + 930);
    const auto *lsi_934 = buffer.data(lsi + 934);
    const auto *lsi_945 = buffer.data(lsi + 945);
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
    const auto *lsi_979 = buffer.data(lsi + 979);
    const auto *lsi_980 = buffer.data(lsi + 980);
    const auto *lsi_982 = buffer.data(lsi + 982);
    const auto *lsi_985 = buffer.data(lsi + 985);
    const auto *lsi_989 = buffer.data(lsi + 989);
    const auto *lsi_994 = buffer.data(lsi + 994);
    const auto *lsi_1007 = buffer.data(lsi + 1007);
    const auto *lsi_1169 = buffer.data(lsi + 1169);
    const auto *lsi_1170 = buffer.data(lsi + 1170);
    const auto *lsi_1171 = buffer.data(lsi + 1171);
    const auto *lsi_1172 = buffer.data(lsi + 1172);
    const auto *lsi_1173 = buffer.data(lsi + 1173);
    const auto *lsi_1174 = buffer.data(lsi + 1174);
    const auto *lsi_1175 = buffer.data(lsi + 1175);
    const auto *lsi_1176 = buffer.data(lsi + 1176);
    const auto *lsi_1179 = buffer.data(lsi + 1179);
    const auto *lsi_1181 = buffer.data(lsi + 1181);
    const auto *lsi_1182 = buffer.data(lsi + 1182);
    const auto *lsi_1185 = buffer.data(lsi + 1185);
    const auto *lsi_1186 = buffer.data(lsi + 1186);
    const auto *lsi_1188 = buffer.data(lsi + 1188);
    const auto *lsi_1190 = buffer.data(lsi + 1190);
    const auto *lsi_1191 = buffer.data(lsi + 1191);
    const auto *lsi_1193 = buffer.data(lsi + 1193);
    const auto *lsi_1194 = buffer.data(lsi + 1194);
    const auto *lsi_1196 = buffer.data(lsi + 1196);
    const auto *lsi_1197 = buffer.data(lsi + 1197);
    const auto *lsi_1198 = buffer.data(lsi + 1198);
    const auto *lsi_1199 = buffer.data(lsi + 1199);
    const auto *lsi_1200 = buffer.data(lsi + 1200);
    const auto *lsi_1201 = buffer.data(lsi + 1201);
    const auto *lsi_1202 = buffer.data(lsi + 1202);
    const auto *lsi_1203 = buffer.data(lsi + 1203);
    const auto *lsi_1207 = buffer.data(lsi + 1207);
    const auto *lsi_1210 = buffer.data(lsi + 1210);
    const auto *lsi_1214 = buffer.data(lsi + 1214);
    const auto *lsi_1216 = buffer.data(lsi + 1216);
    const auto *lsi_1219 = buffer.data(lsi + 1219);
    const auto *lsi_1221 = buffer.data(lsi + 1221);
    const auto *lsi_1222 = buffer.data(lsi + 1222);
    const auto *lsi_1225 = buffer.data(lsi + 1225);
    const auto *lsi_1226 = buffer.data(lsi + 1226);
    const auto *lsi_1227 = buffer.data(lsi + 1227);
    const auto *lsi_1228 = buffer.data(lsi + 1228);
    const auto *lsi_1229 = buffer.data(lsi + 1229);
    const auto *lsi_1230 = buffer.data(lsi + 1230);
    const auto *lsi_1231 = buffer.data(lsi + 1231);
    const auto *lsi_1232 = buffer.data(lsi + 1232);
    const auto *lsi_1237 = buffer.data(lsi + 1237);
    const auto *lsi_1241 = buffer.data(lsi + 1241);
    const auto *lsi_1246 = buffer.data(lsi + 1246);
    const auto *lsi_1252 = buffer.data(lsi + 1252);
    const auto *lsi_1253 = buffer.data(lsi + 1253);
    const auto *lsi_1254 = buffer.data(lsi + 1254);
    const auto *lsi_1255 = buffer.data(lsi + 1255);
    const auto *lsi_1256 = buffer.data(lsi + 1256);
    const auto *lsi_1257 = buffer.data(lsi + 1257);
    const auto *lsi_1259 = buffer.data(lsi + 1259);

    const auto *lsk1_1260 = buffer.data(lsk1 + 1260);
    const auto *lsk1_1265 = buffer.data(lsk1 + 1265);
    const auto *lsk1_1269 = buffer.data(lsk1 + 1269);
    const auto *lsk1_1274 = buffer.data(lsk1 + 1274);
    const auto *lsk1_1280 = buffer.data(lsk1 + 1280);
    const auto *lsk1_1504 = buffer.data(lsk1 + 1504);
    const auto *lsk1_1506 = buffer.data(lsk1 + 1506);
    const auto *lsk1_1507 = buffer.data(lsk1 + 1507);
    const auto *lsk1_1508 = buffer.data(lsk1 + 1508);
    const auto *lsk1_1509 = buffer.data(lsk1 + 1509);
    const auto *lsk1_1511 = buffer.data(lsk1 + 1511);
    const auto *lsk1_1512 = buffer.data(lsk1 + 1512);
    const auto *lsk1_1515 = buffer.data(lsk1 + 1515);
    const auto *lsk1_1517 = buffer.data(lsk1 + 1517);
    const auto *lsk1_1518 = buffer.data(lsk1 + 1518);
    const auto *lsk1_1521 = buffer.data(lsk1 + 1521);
    const auto *lsk1_1522 = buffer.data(lsk1 + 1522);
    const auto *lsk1_1524 = buffer.data(lsk1 + 1524);
    const auto *lsk1_1526 = buffer.data(lsk1 + 1526);
    const auto *lsk1_1527 = buffer.data(lsk1 + 1527);
    const auto *lsk1_1529 = buffer.data(lsk1 + 1529);
    const auto *lsk1_1530 = buffer.data(lsk1 + 1530);
    const auto *lsk1_1532 = buffer.data(lsk1 + 1532);
    const auto *lsk1_1540 = buffer.data(lsk1 + 1540);
    const auto *lsk1_1542 = buffer.data(lsk1 + 1542);
    const auto *lsk1_1543 = buffer.data(lsk1 + 1543);
    const auto *lsk1_1544 = buffer.data(lsk1 + 1544);
    const auto *lsk1_1545 = buffer.data(lsk1 + 1545);
    const auto *lsk1_1547 = buffer.data(lsk1 + 1547);
    const auto *lsk1_1551 = buffer.data(lsk1 + 1551);
    const auto *lsk1_1554 = buffer.data(lsk1 + 1554);
    const auto *lsk1_1558 = buffer.data(lsk1 + 1558);
    const auto *lsk1_1560 = buffer.data(lsk1 + 1560);
    const auto *lsk1_1563 = buffer.data(lsk1 + 1563);
    const auto *lsk1_1565 = buffer.data(lsk1 + 1565);
    const auto *lsk1_1566 = buffer.data(lsk1 + 1566);
    const auto *lsk1_1576 = buffer.data(lsk1 + 1576);
    const auto *lsk1_1578 = buffer.data(lsk1 + 1578);
    const auto *lsk1_1579 = buffer.data(lsk1 + 1579);
    const auto *lsk1_1580 = buffer.data(lsk1 + 1580);
    const auto *lsk1_1581 = buffer.data(lsk1 + 1581);
    const auto *lsk1_1583 = buffer.data(lsk1 + 1583);
    const auto *lsk1_1584 = buffer.data(lsk1 + 1584);
    const auto *lsk1_1589 = buffer.data(lsk1 + 1589);
    const auto *lsk1_1593 = buffer.data(lsk1 + 1593);
    const auto *lsk1_1598 = buffer.data(lsk1 + 1598);
    const auto *lsk1_1604 = buffer.data(lsk1 + 1604);
    const auto *lsk1_1612 = buffer.data(lsk1 + 1612);
    const auto *lsk1_1613 = buffer.data(lsk1 + 1613);
    const auto *lsk1_1614 = buffer.data(lsk1 + 1614);
    const auto *lsk1_1615 = buffer.data(lsk1 + 1615);
    const auto *lsk1_1616 = buffer.data(lsk1 + 1616);
    const auto *lsk1_1617 = buffer.data(lsk1 + 1617);
    const auto *lsk1_1619 = buffer.data(lsk1 + 1619);

    const auto *msh0_924 = buffer.data(msh0 + 924);
    const auto *msh0_925 = buffer.data(msh0 + 925);
    const auto *msh0_926 = buffer.data(msh0 + 926);
    const auto *msh0_927 = buffer.data(msh0 + 927);
    const auto *msh0_928 = buffer.data(msh0 + 928);
    const auto *msh0_929 = buffer.data(msh0 + 929);
    const auto *msh0_930 = buffer.data(msh0 + 930);
    const auto *msh0_931 = buffer.data(msh0 + 931);
    const auto *msh0_932 = buffer.data(msh0 + 932);
    const auto *msh0_933 = buffer.data(msh0 + 933);

    const auto *msh1_924 = buffer.data(msh1 + 924);
    const auto *msh1_925 = buffer.data(msh1 + 925);
    const auto *msh1_926 = buffer.data(msh1 + 926);
    const auto *msh1_927 = buffer.data(msh1 + 927);
    const auto *msh1_928 = buffer.data(msh1 + 928);
    const auto *msh1_929 = buffer.data(msh1 + 929);
    const auto *msh1_930 = buffer.data(msh1 + 930);
    const auto *msh1_931 = buffer.data(msh1 + 931);
    const auto *msh1_932 = buffer.data(msh1 + 932);
    const auto *msh1_933 = buffer.data(msh1 + 933);

    const auto *msi_1169 = buffer.data(msi + 1169);
    const auto *msi_1170 = buffer.data(msi + 1170);
    const auto *msi_1171 = buffer.data(msi + 1171);
    const auto *msi_1172 = buffer.data(msi + 1172);
    const auto *msi_1173 = buffer.data(msi + 1173);
    const auto *msi_1174 = buffer.data(msi + 1174);
    const auto *msi_1175 = buffer.data(msi + 1175);
    const auto *msi_1176 = buffer.data(msi + 1176);
    const auto *msi_1178 = buffer.data(msi + 1178);
    const auto *msi_1179 = buffer.data(msi + 1179);
    const auto *msi_1181 = buffer.data(msi + 1181);
    const auto *msi_1182 = buffer.data(msi + 1182);
    const auto *msi_1185 = buffer.data(msi + 1185);
    const auto *msi_1186 = buffer.data(msi + 1186);
    const auto *msi_1190 = buffer.data(msi + 1190);
    const auto *msi_1197 = buffer.data(msi + 1197);
    const auto *msi_1198 = buffer.data(msi + 1198);
    const auto *msi_1199 = buffer.data(msi + 1199);
    const auto *msi_1200 = buffer.data(msi + 1200);
    const auto *msi_1201 = buffer.data(msi + 1201);
    const auto *msi_1202 = buffer.data(msi + 1202);
    const auto *msi_1203 = buffer.data(msi + 1203);
    const auto *msi_1204 = buffer.data(msi + 1204);
    const auto *msi_1206 = buffer.data(msi + 1206);
    const auto *msi_1207 = buffer.data(msi + 1207);
    const auto *msi_1209 = buffer.data(msi + 1209);
    const auto *msi_1210 = buffer.data(msi + 1210);
    const auto *msi_1213 = buffer.data(msi + 1213);
    const auto *msi_1214 = buffer.data(msi + 1214);
    const auto *msi_1218 = buffer.data(msi + 1218);
    const auto *msi_1225 = buffer.data(msi + 1225);
    const auto *msi_1226 = buffer.data(msi + 1226);
    const auto *msi_1227 = buffer.data(msi + 1227);
    const auto *msi_1228 = buffer.data(msi + 1228);
    const auto *msi_1229 = buffer.data(msi + 1229);
    const auto *msi_1230 = buffer.data(msi + 1230);
    const auto *msi_1231 = buffer.data(msi + 1231);
    const auto *msi_1232 = buffer.data(msi + 1232);
    const auto *msi_1233 = buffer.data(msi + 1233);
    const auto *msi_1234 = buffer.data(msi + 1234);
    const auto *msi_1235 = buffer.data(msi + 1235);
    const auto *msi_1236 = buffer.data(msi + 1236);
    const auto *msi_1237 = buffer.data(msi + 1237);
    const auto *msi_1238 = buffer.data(msi + 1238);
    const auto *msi_1239 = buffer.data(msi + 1239);
    const auto *msi_1240 = buffer.data(msi + 1240);
    const auto *msi_1241 = buffer.data(msi + 1241);
    const auto *msi_1242 = buffer.data(msi + 1242);
    const auto *msi_1243 = buffer.data(msi + 1243);
    const auto *msi_1244 = buffer.data(msi + 1244);
    const auto *msi_1245 = buffer.data(msi + 1245);
    const auto *msi_1246 = buffer.data(msi + 1246);
    const auto *msi_1252 = buffer.data(msi + 1252);
    const auto *msi_1253 = buffer.data(msi + 1253);
    const auto *msi_1254 = buffer.data(msi + 1254);
    const auto *msi_1255 = buffer.data(msi + 1255);
    const auto *msi_1256 = buffer.data(msi + 1256);
    const auto *msi_1257 = buffer.data(msi + 1257);
    const auto *msi_1259 = buffer.data(msi + 1259);

#pragma omp simd aligned(t_1497, t_1498, t_1499, t_1500, t_1501, pc_x, lsi_1169, lsi_1170, \
                         lsi_1171, lsi_1172, lsi_1173, msi_1169, msi_1170, msi_1171, msi_1172, \
                         msi_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1497[k] = f_13 * lsi_1169[k]
                    + f_3 * pc_x[k] * msi_1169[k];

        t_1498[k] = f_13 * lsi_1170[k]
                    + f_3 * pc_x[k] * msi_1170[k];

        t_1499[k] = f_13 * lsi_1171[k]
                    + f_3 * pc_x[k] * msi_1171[k];

        t_1500[k] = f_13 * lsi_1172[k]
                    + f_3 * pc_x[k] * msi_1172[k];

        t_1501[k] = f_13 * lsi_1173[k]
                    + f_3 * pc_x[k] * msi_1173[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, t_1505, pa_x, pc_x, pc_z, lsk0_1504, lsi_917, \
                         lsi_1174, lsi_1175, lsk1_1504, msi_1169, msi_1174, \
                         msi_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_13 * lsi_1174[k]
                    + f_3 * pc_x[k] * msi_1174[k];

        t_1503[k] = f_13 * lsi_1175[k]
                    + f_3 * pc_x[k] * msi_1175[k];

        t_1504[k] = pa_x[k] * lsk0_1504[k]
                    - f_12 * pc_x[k] * lsk1_1504[k];

        t_1505[k] = f_17 * lsi_917[k]
                    + f_3 * pc_z[k] * msi_1169[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, t_1509, pa_x, pc_x, lsk0_1506, lsk0_1507, \
                         lsk0_1508, lsk0_1509, lsk1_1506, lsk1_1507, lsk1_1508, \
                         lsk1_1509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = pa_x[k] * lsk0_1506[k]
                    - f_12 * pc_x[k] * lsk1_1506[k];

        t_1507[k] = pa_x[k] * lsk0_1507[k]
                    - f_12 * pc_x[k] * lsk1_1507[k];

        t_1508[k] = pa_x[k] * lsk0_1508[k]
                    - f_12 * pc_x[k] * lsk1_1508[k];

        t_1509[k] = pa_x[k] * lsk0_1509[k]
                    - f_12 * pc_x[k] * lsk1_1509[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, t_1513, pa_x, pc_x, pc_y, lsk0_1511, \
                         lsk0_1512, lsi_951, lsi_952, lsi_1176, lsk1_1511, lsk1_1512, \
                         msi_1175, msi_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_15 * lsi_951[k]
                    + f_3 * pc_y[k] * msi_1175[k];

        t_1511[k] = pa_x[k] * lsk0_1511[k]
                    - f_12 * pc_x[k] * lsk1_1511[k];

        t_1512[k] = pa_x[k] * lsk0_1512[k]
                    + f_21 * lsi_1176[k]
                    - f_12 * pc_x[k] * lsk1_1512[k];

        t_1513[k] = f_14 * lsi_952[k]
                    + f_3 * pc_y[k] * msi_1176[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, pa_x, pc_x, pc_y, pc_z, lsk0_1515, lsi_924, \
                         lsi_954, lsi_1179, lsk1_1515, msi_1176, \
                         msi_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_22 * lsi_924[k]
                    + f_3 * pc_z[k] * msi_1176[k];

        t_1515[k] = pa_x[k] * lsk0_1515[k]
                    + f_17 * lsi_1179[k]
                    - f_12 * pc_x[k] * lsk1_1515[k];

        t_1516[k] = f_14 * lsi_954[k]
                    + f_3 * pc_y[k] * msi_1178[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, pa_x, pc_x, pc_z, lsk0_1517, lsk0_1518, \
                         lsi_927, lsi_1181, lsi_1182, lsk1_1517, lsk1_1518, \
                         msi_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = pa_x[k] * lsk0_1517[k]
                    + f_17 * lsi_1181[k]
                    - f_12 * pc_x[k] * lsk1_1517[k];

        t_1518[k] = pa_x[k] * lsk0_1518[k]
                    + f_16 * lsi_1182[k]
                    - f_12 * pc_x[k] * lsk1_1518[k];

        t_1519[k] = f_22 * lsi_927[k]
                    + f_3 * pc_z[k] * msi_1179[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pa_x, pc_x, pc_y, lsk0_1521, lsk0_1522, \
                         lsi_957, lsi_1185, lsi_1186, lsk1_1521, lsk1_1522, \
                         msi_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_14 * lsi_957[k]
                    + f_3 * pc_y[k] * msi_1181[k];

        t_1521[k] = pa_x[k] * lsk0_1521[k]
                    + f_16 * lsi_1185[k]
                    - f_12 * pc_x[k] * lsk1_1521[k];

        t_1522[k] = pa_x[k] * lsk0_1522[k]
                    + f_15 * lsi_1186[k]
                    - f_12 * pc_x[k] * lsk1_1522[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pa_x, pc_x, pc_y, pc_z, lsk0_1524, lsi_930, \
                         lsi_961, lsi_1188, lsk1_1524, msi_1182, \
                         msi_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_22 * lsi_930[k]
                    + f_3 * pc_z[k] * msi_1182[k];

        t_1524[k] = pa_x[k] * lsk0_1524[k]
                    + f_15 * lsi_1188[k]
                    - f_12 * pc_x[k] * lsk1_1524[k];

        t_1525[k] = f_14 * lsi_961[k]
                    + f_3 * pc_y[k] * msi_1185[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pa_x, pc_x, pc_z, lsk0_1526, lsk0_1527, \
                         lsi_934, lsi_1190, lsi_1191, lsk1_1526, lsk1_1527, \
                         msi_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = pa_x[k] * lsk0_1526[k]
                    + f_15 * lsi_1190[k]
                    - f_12 * pc_x[k] * lsk1_1526[k];

        t_1527[k] = pa_x[k] * lsk0_1527[k]
                    + f_14 * lsi_1191[k]
                    - f_12 * pc_x[k] * lsk1_1527[k];

        t_1528[k] = f_22 * lsi_934[k]
                    + f_3 * pc_z[k] * msi_1186[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, pa_x, pc_x, pc_y, lsk0_1529, lsk0_1530, \
                         lsi_966, lsi_1193, lsi_1194, lsk1_1529, lsk1_1530, \
                         msi_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = pa_x[k] * lsk0_1529[k]
                    + f_14 * lsi_1193[k]
                    - f_12 * pc_x[k] * lsk1_1529[k];

        t_1530[k] = pa_x[k] * lsk0_1530[k]
                    + f_14 * lsi_1194[k]
                    - f_12 * pc_x[k] * lsk1_1530[k];

        t_1531[k] = f_14 * lsi_966[k]
                    + f_3 * pc_y[k] * msi_1190[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, t_1535, pa_x, pc_x, lsk0_1532, lsi_1196, \
                         lsi_1197, lsi_1198, lsi_1199, lsk1_1532, msi_1197, msi_1198, \
                         msi_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = pa_x[k] * lsk0_1532[k]
                    + f_14 * lsi_1196[k]
                    - f_12 * pc_x[k] * lsk1_1532[k];

        t_1533[k] = f_13 * lsi_1197[k]
                    + f_3 * pc_x[k] * msi_1197[k];

        t_1534[k] = f_13 * lsi_1198[k]
                    + f_3 * pc_x[k] * msi_1198[k];

        t_1535[k] = f_13 * lsi_1199[k]
                    + f_3 * pc_x[k] * msi_1199[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pc_x, lsi_1200, lsi_1201, lsi_1202, \
                         lsi_1203, msi_1200, msi_1201, msi_1202, \
                         msi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_13 * lsi_1200[k]
                    + f_3 * pc_x[k] * msi_1200[k];

        t_1537[k] = f_13 * lsi_1201[k]
                    + f_3 * pc_x[k] * msi_1201[k];

        t_1538[k] = f_13 * lsi_1202[k]
                    + f_3 * pc_x[k] * msi_1202[k];

        t_1539[k] = f_13 * lsi_1203[k]
                    + f_3 * pc_x[k] * msi_1203[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, t_1543, pa_x, pc_x, pc_z, lsk0_1540, \
                         lsk0_1542, lsk0_1543, lsi_945, lsk1_1540, lsk1_1542, lsk1_1543, \
                         msi_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = pa_x[k] * lsk0_1540[k]
                    - f_12 * pc_x[k] * lsk1_1540[k];

        t_1541[k] = f_22 * lsi_945[k]
                    + f_3 * pc_z[k] * msi_1197[k];

        t_1542[k] = pa_x[k] * lsk0_1542[k]
                    - f_12 * pc_x[k] * lsk1_1542[k];

        t_1543[k] = pa_x[k] * lsk0_1543[k]
                    - f_12 * pc_x[k] * lsk1_1543[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, t_1547, pa_x, pc_x, pc_y, lsk0_1544, \
                         lsk0_1545, lsk0_1547, lsi_979, lsk1_1544, lsk1_1545, lsk1_1547, \
                         msi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = pa_x[k] * lsk0_1544[k]
                    - f_12 * pc_x[k] * lsk1_1544[k];

        t_1545[k] = pa_x[k] * lsk0_1545[k]
                    - f_12 * pc_x[k] * lsk1_1545[k];

        t_1546[k] = f_14 * lsi_979[k]
                    + f_3 * pc_y[k] * msi_1203[k];

        t_1547[k] = pa_x[k] * lsk0_1547[k]
                    - f_12 * pc_x[k] * lsk1_1547[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, pa_y, pc_y, pc_z, lsk0_1260, lsi_952, \
                         lsi_980, lsk1_1260, msi_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = pa_y[k] * lsk0_1260[k]
                    - f_12 * pc_y[k] * lsk1_1260[k];

        t_1549[k] = f_13 * lsi_980[k]
                    + f_3 * pc_y[k] * msi_1204[k];

        t_1550[k] = f_21 * lsi_952[k]
                    + f_3 * pc_z[k] * msi_1204[k];
    }

#pragma omp simd aligned(t_1551, t_1552, t_1553, pa_x, pa_y, pc_x, pc_y, lsk0_1265, lsk0_1551, \
                         lsi_982, lsi_1207, lsk1_1265, lsk1_1551, \
                         msi_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1551[k] = pa_x[k] * lsk0_1551[k]
                    + f_17 * lsi_1207[k]
                    - f_12 * pc_x[k] * lsk1_1551[k];

        t_1552[k] = f_13 * lsi_982[k]
                    + f_3 * pc_y[k] * msi_1206[k];

        t_1553[k] = pa_y[k] * lsk0_1265[k]
                    - f_12 * pc_y[k] * lsk1_1265[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, pa_x, pc_x, pc_y, pc_z, lsk0_1554, lsi_955, \
                         lsi_985, lsi_1210, lsk1_1554, msi_1207, \
                         msi_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = pa_x[k] * lsk0_1554[k]
                    + f_16 * lsi_1210[k]
                    - f_12 * pc_x[k] * lsk1_1554[k];

        t_1555[k] = f_21 * lsi_955[k]
                    + f_3 * pc_z[k] * msi_1207[k];

        t_1556[k] = f_13 * lsi_985[k]
                    + f_3 * pc_y[k] * msi_1209[k];
    }

#pragma omp simd aligned(t_1557, t_1558, t_1559, pa_x, pa_y, pc_x, pc_y, pc_z, lsk0_1269, \
                         lsk0_1558, lsi_958, lsi_1214, lsk1_1269, lsk1_1558, \
                         msi_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1557[k] = pa_y[k] * lsk0_1269[k]
                    - f_12 * pc_y[k] * lsk1_1269[k];

        t_1558[k] = pa_x[k] * lsk0_1558[k]
                    + f_15 * lsi_1214[k]
                    - f_12 * pc_x[k] * lsk1_1558[k];

        t_1559[k] = f_21 * lsi_958[k]
                    + f_3 * pc_z[k] * msi_1210[k];
    }

#pragma omp simd aligned(t_1560, t_1561, t_1562, pa_x, pa_y, pc_x, pc_y, lsk0_1274, lsk0_1560, \
                         lsi_989, lsi_1216, lsk1_1274, lsk1_1560, \
                         msi_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1560[k] = pa_x[k] * lsk0_1560[k]
                    + f_15 * lsi_1216[k]
                    - f_12 * pc_x[k] * lsk1_1560[k];

        t_1561[k] = f_13 * lsi_989[k]
                    + f_3 * pc_y[k] * msi_1213[k];

        t_1562[k] = pa_y[k] * lsk0_1274[k]
                    - f_12 * pc_y[k] * lsk1_1274[k];
    }

#pragma omp simd aligned(t_1563, t_1564, t_1565, pa_x, pc_x, pc_z, lsk0_1563, lsk0_1565, \
                         lsi_962, lsi_1219, lsi_1221, lsk1_1563, lsk1_1565, \
                         msi_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1563[k] = pa_x[k] * lsk0_1563[k]
                    + f_14 * lsi_1219[k]
                    - f_12 * pc_x[k] * lsk1_1563[k];

        t_1564[k] = f_21 * lsi_962[k]
                    + f_3 * pc_z[k] * msi_1214[k];

        t_1565[k] = pa_x[k] * lsk0_1565[k]
                    + f_14 * lsi_1221[k]
                    - f_12 * pc_x[k] * lsk1_1565[k];
    }

#pragma omp simd aligned(t_1566, t_1567, t_1568, pa_x, pa_y, pc_x, pc_y, lsk0_1280, lsk0_1566, \
                         lsi_994, lsi_1222, lsk1_1280, lsk1_1566, \
                         msi_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1566[k] = pa_x[k] * lsk0_1566[k]
                    + f_14 * lsi_1222[k]
                    - f_12 * pc_x[k] * lsk1_1566[k];

        t_1567[k] = f_13 * lsi_994[k]
                    + f_3 * pc_y[k] * msi_1218[k];

        t_1568[k] = pa_y[k] * lsk0_1280[k]
                    - f_12 * pc_y[k] * lsk1_1280[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, t_1572, t_1573, pc_x, lsi_1225, lsi_1226, \
                         lsi_1227, lsi_1228, lsi_1229, msi_1225, msi_1226, msi_1227, msi_1228, \
                         msi_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = f_13 * lsi_1225[k]
                    + f_3 * pc_x[k] * msi_1225[k];

        t_1570[k] = f_13 * lsi_1226[k]
                    + f_3 * pc_x[k] * msi_1226[k];

        t_1571[k] = f_13 * lsi_1227[k]
                    + f_3 * pc_x[k] * msi_1227[k];

        t_1572[k] = f_13 * lsi_1228[k]
                    + f_3 * pc_x[k] * msi_1228[k];

        t_1573[k] = f_13 * lsi_1229[k]
                    + f_3 * pc_x[k] * msi_1229[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pa_x, pc_x, pc_z, lsk0_1576, lsi_973, \
                         lsi_1230, lsi_1231, lsk1_1576, msi_1225, msi_1230, \
                         msi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = f_13 * lsi_1230[k]
                    + f_3 * pc_x[k] * msi_1230[k];

        t_1575[k] = f_13 * lsi_1231[k]
                    + f_3 * pc_x[k] * msi_1231[k];

        t_1576[k] = pa_x[k] * lsk0_1576[k]
                    - f_12 * pc_x[k] * lsk1_1576[k];

        t_1577[k] = f_21 * lsi_973[k]
                    + f_3 * pc_z[k] * msi_1225[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, t_1581, pa_x, pc_x, lsk0_1578, lsk0_1579, \
                         lsk0_1580, lsk0_1581, lsk1_1578, lsk1_1579, lsk1_1580, \
                         lsk1_1581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pa_x[k] * lsk0_1578[k]
                    - f_12 * pc_x[k] * lsk1_1578[k];

        t_1579[k] = pa_x[k] * lsk0_1579[k]
                    - f_12 * pc_x[k] * lsk1_1579[k];

        t_1580[k] = pa_x[k] * lsk0_1580[k]
                    - f_12 * pc_x[k] * lsk1_1580[k];

        t_1581[k] = pa_x[k] * lsk0_1581[k]
                    - f_12 * pc_x[k] * lsk1_1581[k];
    }

#pragma omp simd aligned(t_1582, t_1583, t_1584, t_1585, pa_x, pc_x, pc_y, lsk0_1583, \
                         lsk0_1584, lsi_1007, lsi_1232, lsk1_1583, lsk1_1584, msi_1231, \
                         msi_1232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1582[k] = f_13 * lsi_1007[k]
                    + f_3 * pc_y[k] * msi_1231[k];

        t_1583[k] = pa_x[k] * lsk0_1583[k]
                    - f_12 * pc_x[k] * lsk1_1583[k];

        t_1584[k] = pa_x[k] * lsk0_1584[k]
                    + f_21 * lsi_1232[k]
                    - f_12 * pc_x[k] * lsk1_1584[k];

        t_1585[k] = f_3 * pc_y[k] * msi_1232[k];
    }

#pragma omp simd aligned(t_1586, t_1587, t_1588, pc_y, pc_z, lsi_980, msh0_924, msh1_924, \
                         msi_1232, msi_1233, msi_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1586[k] = f_18 * lsi_980[k]
                    + f_3 * pc_z[k] * msi_1232[k];

        t_1587[k] = f_4 * msh0_924[k]
                    - f_5 * msh1_924[k]
                    + f_3 * pc_y[k] * msi_1233[k];

        t_1588[k] = f_3 * pc_y[k] * msi_1234[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, pa_x, pc_x, pc_y, lsk0_1589, lsi_1237, \
                         lsk1_1589, msh0_925, msh0_926, msh1_925, msh1_926, msi_1235, \
                         msi_1236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = pa_x[k] * lsk0_1589[k]
                    + f_17 * lsi_1237[k]
                    - f_12 * pc_x[k] * lsk1_1589[k];

        t_1590[k] = f_6 * msh0_925[k]
                    - f_7 * msh1_925[k]
                    + f_3 * pc_y[k] * msi_1235[k];

        t_1591[k] = f_4 * msh0_926[k]
                    - f_5 * msh1_926[k]
                    + f_3 * pc_y[k] * msi_1236[k];
    }

#pragma omp simd aligned(t_1592, t_1593, t_1594, pa_x, pc_x, pc_y, lsk0_1593, lsi_1241, \
                         lsk1_1593, msh0_927, msh1_927, msi_1237, \
                         msi_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1592[k] = f_3 * pc_y[k] * msi_1237[k];

        t_1593[k] = pa_x[k] * lsk0_1593[k]
                    + f_16 * lsi_1241[k]
                    - f_12 * pc_x[k] * lsk1_1593[k];

        t_1594[k] = f_8 * msh0_927[k]
                    - f_9 * msh1_927[k]
                    + f_3 * pc_y[k] * msi_1238[k];
    }

#pragma omp simd aligned(t_1595, t_1596, t_1597, pc_y, msh0_928, msh0_929, msh1_928, msh1_929, \
                         msi_1239, msi_1240, msi_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1595[k] = f_6 * msh0_928[k]
                    - f_7 * msh1_928[k]
                    + f_3 * pc_y[k] * msi_1239[k];

        t_1596[k] = f_4 * msh0_929[k]
                    - f_5 * msh1_929[k]
                    + f_3 * pc_y[k] * msi_1240[k];

        t_1597[k] = f_3 * pc_y[k] * msi_1241[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, pa_x, pc_x, pc_y, lsk0_1598, lsi_1246, \
                         lsk1_1598, msh0_930, msh0_931, msh1_930, msh1_931, msi_1242, \
                         msi_1243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = pa_x[k] * lsk0_1598[k]
                    + f_15 * lsi_1246[k]
                    - f_12 * pc_x[k] * lsk1_1598[k];

        t_1599[k] = f_10 * msh0_930[k]
                    - f_11 * msh1_930[k]
                    + f_3 * pc_y[k] * msi_1242[k];

        t_1600[k] = f_8 * msh0_931[k]
                    - f_9 * msh1_931[k]
                    + f_3 * pc_y[k] * msi_1243[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, pc_y, msh0_932, msh0_933, msh1_932, msh1_933, \
                         msi_1244, msi_1245, msi_1246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = f_6 * msh0_932[k]
                    - f_7 * msh1_932[k]
                    + f_3 * pc_y[k] * msi_1244[k];

        t_1602[k] = f_4 * msh0_933[k]
                    - f_5 * msh1_933[k]
                    + f_3 * pc_y[k] * msi_1245[k];

        t_1603[k] = f_3 * pc_y[k] * msi_1246[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, t_1607, pa_x, pc_x, lsk0_1604, lsi_1252, \
                         lsi_1253, lsi_1254, lsi_1255, lsk1_1604, msi_1253, msi_1254, \
                         msi_1255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = pa_x[k] * lsk0_1604[k]
                    + f_14 * lsi_1252[k]
                    - f_12 * pc_x[k] * lsk1_1604[k];

        t_1605[k] = f_13 * lsi_1253[k]
                    + f_3 * pc_x[k] * msi_1253[k];

        t_1606[k] = f_13 * lsi_1254[k]
                    + f_3 * pc_x[k] * msi_1254[k];

        t_1607[k] = f_13 * lsi_1255[k]
                    + f_3 * pc_x[k] * msi_1255[k];
    }

#pragma omp simd aligned(t_1608, t_1609, t_1610, t_1611, pc_x, pc_y, lsi_1256, lsi_1257, \
                         lsi_1259, msi_1252, msi_1256, msi_1257, \
                         msi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1608[k] = f_13 * lsi_1256[k]
                    + f_3 * pc_x[k] * msi_1256[k];

        t_1609[k] = f_13 * lsi_1257[k]
                    + f_3 * pc_x[k] * msi_1257[k];

        t_1610[k] = f_3 * pc_y[k] * msi_1252[k];

        t_1611[k] = f_13 * lsi_1259[k]
                    + f_3 * pc_x[k] * msi_1259[k];
    }

#pragma omp simd aligned(t_1612, t_1613, t_1614, t_1615, pa_x, pc_x, lsk0_1612, lsk0_1613, \
                         lsk0_1614, lsk0_1615, lsk1_1612, lsk1_1613, lsk1_1614, \
                         lsk1_1615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1612[k] = pa_x[k] * lsk0_1612[k]
                    - f_12 * pc_x[k] * lsk1_1612[k];

        t_1613[k] = pa_x[k] * lsk0_1613[k]
                    - f_12 * pc_x[k] * lsk1_1613[k];

        t_1614[k] = pa_x[k] * lsk0_1614[k]
                    - f_12 * pc_x[k] * lsk1_1614[k];

        t_1615[k] = pa_x[k] * lsk0_1615[k]
                    - f_12 * pc_x[k] * lsk1_1615[k];
    }

#pragma omp simd aligned(t_1616, t_1617, t_1618, t_1619, pa_x, pc_x, pc_y, lsk0_1616, \
                         lsk0_1617, lsk0_1619, lsk1_1616, lsk1_1617, lsk1_1619, \
                         msi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1616[k] = pa_x[k] * lsk0_1616[k]
                    - f_12 * pc_x[k] * lsk1_1616[k];

        t_1617[k] = pa_x[k] * lsk0_1617[k]
                    - f_12 * pc_x[k] * lsk1_1617[k];

        t_1618[k] = f_3 * pc_y[k] * msi_1259[k];

        t_1619[k] = pa_x[k] * lsk0_1619[k]
                    - f_12 * pc_x[k] * lsk1_1619[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsk0,
                                                           const size_t lsi, const size_t lsk1,
                                                           const size_t msh0, const size_t msh1,
                                                           const size_t msi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_18 = 4.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_1296 = buffer.data(lsk0 + 1296);
    const auto *lsk0_1297 = buffer.data(lsk0 + 1297);
    const auto *lsk0_1299 = buffer.data(lsk0 + 1299);
    const auto *lsk0_1302 = buffer.data(lsk0 + 1302);
    const auto *lsk0_1306 = buffer.data(lsk0 + 1306);
    const auto *lsk0_1311 = buffer.data(lsk0 + 1311);
    const auto *lsk0_1324 = buffer.data(lsk0 + 1324);
    const auto *lsk0_1326 = buffer.data(lsk0 + 1326);
    const auto *lsk0_1327 = buffer.data(lsk0 + 1327);
    const auto *lsk0_1328 = buffer.data(lsk0 + 1328);
    const auto *lsk0_1329 = buffer.data(lsk0 + 1329);

    const auto *lsi_1029 = buffer.data(lsi + 1029);
    const auto *lsi_1030 = buffer.data(lsi + 1030);
    const auto *lsi_1031 = buffer.data(lsi + 1031);
    const auto *lsi_1032 = buffer.data(lsi + 1032);
    const auto *lsi_1033 = buffer.data(lsi + 1033);
    const auto *lsi_1035 = buffer.data(lsi + 1035);
    const auto *lsi_1057 = buffer.data(lsi + 1057);
    const auto *lsi_1063 = buffer.data(lsi + 1063);
    const auto *lsi_1085 = buffer.data(lsi + 1085);
    const auto *lsi_1087 = buffer.data(lsi + 1087);
    const auto *lsi_1088 = buffer.data(lsi + 1088);
    const auto *lsi_1089 = buffer.data(lsi + 1089);
    const auto *lsi_1090 = buffer.data(lsi + 1090);
    const auto *lsi_1091 = buffer.data(lsi + 1091);

    const auto *lsk1_1296 = buffer.data(lsk1 + 1296);
    const auto *lsk1_1297 = buffer.data(lsk1 + 1297);
    const auto *lsk1_1299 = buffer.data(lsk1 + 1299);
    const auto *lsk1_1302 = buffer.data(lsk1 + 1302);
    const auto *lsk1_1306 = buffer.data(lsk1 + 1306);
    const auto *lsk1_1311 = buffer.data(lsk1 + 1311);
    const auto *lsk1_1324 = buffer.data(lsk1 + 1324);
    const auto *lsk1_1326 = buffer.data(lsk1 + 1326);
    const auto *lsk1_1327 = buffer.data(lsk1 + 1327);
    const auto *lsk1_1328 = buffer.data(lsk1 + 1328);
    const auto *lsk1_1329 = buffer.data(lsk1 + 1329);

    const auto *msh0_945 = buffer.data(msh0 + 945);
    const auto *msh0_946 = buffer.data(msh0 + 946);
    const auto *msh0_948 = buffer.data(msh0 + 948);
    const auto *msh0_950 = buffer.data(msh0 + 950);
    const auto *msh0_951 = buffer.data(msh0 + 951);
    const auto *msh0_953 = buffer.data(msh0 + 953);
    const auto *msh0_954 = buffer.data(msh0 + 954);
    const auto *msh0_955 = buffer.data(msh0 + 955);
    const auto *msh0_957 = buffer.data(msh0 + 957);
    const auto *msh0_958 = buffer.data(msh0 + 958);
    const auto *msh0_959 = buffer.data(msh0 + 959);
    const auto *msh0_960 = buffer.data(msh0 + 960);
    const auto *msh0_961 = buffer.data(msh0 + 961);
    const auto *msh0_962 = buffer.data(msh0 + 962);
    const auto *msh0_963 = buffer.data(msh0 + 963);
    const auto *msh0_964 = buffer.data(msh0 + 964);
    const auto *msh0_965 = buffer.data(msh0 + 965);
    const auto *msh0_968 = buffer.data(msh0 + 968);
    const auto *msh0_970 = buffer.data(msh0 + 970);
    const auto *msh0_971 = buffer.data(msh0 + 971);
    const auto *msh0_973 = buffer.data(msh0 + 973);
    const auto *msh0_974 = buffer.data(msh0 + 974);
    const auto *msh0_975 = buffer.data(msh0 + 975);
    const auto *msh0_977 = buffer.data(msh0 + 977);
    const auto *msh0_978 = buffer.data(msh0 + 978);
    const auto *msh0_979 = buffer.data(msh0 + 979);
    const auto *msh0_980 = buffer.data(msh0 + 980);
    const auto *msh0_982 = buffer.data(msh0 + 982);
    const auto *msh0_983 = buffer.data(msh0 + 983);
    const auto *msh0_984 = buffer.data(msh0 + 984);
    const auto *msh0_985 = buffer.data(msh0 + 985);
    const auto *msh0_986 = buffer.data(msh0 + 986);
    const auto *msh0_987 = buffer.data(msh0 + 987);
    const auto *msh0_988 = buffer.data(msh0 + 988);
    const auto *msh0_989 = buffer.data(msh0 + 989);
    const auto *msh0_990 = buffer.data(msh0 + 990);
    const auto *msh0_991 = buffer.data(msh0 + 991);
    const auto *msh0_992 = buffer.data(msh0 + 992);
    const auto *msh0_993 = buffer.data(msh0 + 993);
    const auto *msh0_994 = buffer.data(msh0 + 994);
    const auto *msh0_995 = buffer.data(msh0 + 995);
    const auto *msh0_996 = buffer.data(msh0 + 996);
    const auto *msh0_997 = buffer.data(msh0 + 997);
    const auto *msh0_998 = buffer.data(msh0 + 998);
    const auto *msh0_999 = buffer.data(msh0 + 999);
    const auto *msh0_1000 = buffer.data(msh0 + 1000);
    const auto *msh0_1001 = buffer.data(msh0 + 1001);
    const auto *msh0_1002 = buffer.data(msh0 + 1002);
    const auto *msh0_1003 = buffer.data(msh0 + 1003);
    const auto *msh0_1004 = buffer.data(msh0 + 1004);
    const auto *msh0_1005 = buffer.data(msh0 + 1005);
    const auto *msh0_1006 = buffer.data(msh0 + 1006);
    const auto *msh0_1007 = buffer.data(msh0 + 1007);
    const auto *msh0_1008 = buffer.data(msh0 + 1008);
    const auto *msh0_1009 = buffer.data(msh0 + 1009);
    const auto *msh0_1010 = buffer.data(msh0 + 1010);
    const auto *msh0_1011 = buffer.data(msh0 + 1011);
    const auto *msh0_1012 = buffer.data(msh0 + 1012);
    const auto *msh0_1013 = buffer.data(msh0 + 1013);
    const auto *msh0_1014 = buffer.data(msh0 + 1014);
    const auto *msh0_1015 = buffer.data(msh0 + 1015);
    const auto *msh0_1016 = buffer.data(msh0 + 1016);
    const auto *msh0_1017 = buffer.data(msh0 + 1017);

    const auto *msh1_945 = buffer.data(msh1 + 945);
    const auto *msh1_946 = buffer.data(msh1 + 946);
    const auto *msh1_948 = buffer.data(msh1 + 948);
    const auto *msh1_950 = buffer.data(msh1 + 950);
    const auto *msh1_951 = buffer.data(msh1 + 951);
    const auto *msh1_953 = buffer.data(msh1 + 953);
    const auto *msh1_954 = buffer.data(msh1 + 954);
    const auto *msh1_955 = buffer.data(msh1 + 955);
    const auto *msh1_957 = buffer.data(msh1 + 957);
    const auto *msh1_958 = buffer.data(msh1 + 958);
    const auto *msh1_959 = buffer.data(msh1 + 959);
    const auto *msh1_960 = buffer.data(msh1 + 960);
    const auto *msh1_961 = buffer.data(msh1 + 961);
    const auto *msh1_962 = buffer.data(msh1 + 962);
    const auto *msh1_963 = buffer.data(msh1 + 963);
    const auto *msh1_964 = buffer.data(msh1 + 964);
    const auto *msh1_965 = buffer.data(msh1 + 965);
    const auto *msh1_968 = buffer.data(msh1 + 968);
    const auto *msh1_970 = buffer.data(msh1 + 970);
    const auto *msh1_971 = buffer.data(msh1 + 971);
    const auto *msh1_973 = buffer.data(msh1 + 973);
    const auto *msh1_974 = buffer.data(msh1 + 974);
    const auto *msh1_975 = buffer.data(msh1 + 975);
    const auto *msh1_977 = buffer.data(msh1 + 977);
    const auto *msh1_978 = buffer.data(msh1 + 978);
    const auto *msh1_979 = buffer.data(msh1 + 979);
    const auto *msh1_980 = buffer.data(msh1 + 980);
    const auto *msh1_982 = buffer.data(msh1 + 982);
    const auto *msh1_983 = buffer.data(msh1 + 983);
    const auto *msh1_984 = buffer.data(msh1 + 984);
    const auto *msh1_985 = buffer.data(msh1 + 985);
    const auto *msh1_986 = buffer.data(msh1 + 986);
    const auto *msh1_987 = buffer.data(msh1 + 987);
    const auto *msh1_988 = buffer.data(msh1 + 988);
    const auto *msh1_989 = buffer.data(msh1 + 989);
    const auto *msh1_990 = buffer.data(msh1 + 990);
    const auto *msh1_991 = buffer.data(msh1 + 991);
    const auto *msh1_992 = buffer.data(msh1 + 992);
    const auto *msh1_993 = buffer.data(msh1 + 993);
    const auto *msh1_994 = buffer.data(msh1 + 994);
    const auto *msh1_995 = buffer.data(msh1 + 995);
    const auto *msh1_996 = buffer.data(msh1 + 996);
    const auto *msh1_997 = buffer.data(msh1 + 997);
    const auto *msh1_998 = buffer.data(msh1 + 998);
    const auto *msh1_999 = buffer.data(msh1 + 999);
    const auto *msh1_1000 = buffer.data(msh1 + 1000);
    const auto *msh1_1001 = buffer.data(msh1 + 1001);
    const auto *msh1_1002 = buffer.data(msh1 + 1002);
    const auto *msh1_1003 = buffer.data(msh1 + 1003);
    const auto *msh1_1004 = buffer.data(msh1 + 1004);
    const auto *msh1_1005 = buffer.data(msh1 + 1005);
    const auto *msh1_1006 = buffer.data(msh1 + 1006);
    const auto *msh1_1007 = buffer.data(msh1 + 1007);
    const auto *msh1_1008 = buffer.data(msh1 + 1008);
    const auto *msh1_1009 = buffer.data(msh1 + 1009);
    const auto *msh1_1010 = buffer.data(msh1 + 1010);
    const auto *msh1_1011 = buffer.data(msh1 + 1011);
    const auto *msh1_1012 = buffer.data(msh1 + 1012);
    const auto *msh1_1013 = buffer.data(msh1 + 1013);
    const auto *msh1_1014 = buffer.data(msh1 + 1014);
    const auto *msh1_1015 = buffer.data(msh1 + 1015);
    const auto *msh1_1016 = buffer.data(msh1 + 1016);
    const auto *msh1_1017 = buffer.data(msh1 + 1017);

    const auto *msi_1260 = buffer.data(msi + 1260);
    const auto *msi_1261 = buffer.data(msi + 1261);
    const auto *msi_1263 = buffer.data(msi + 1263);
    const auto *msi_1265 = buffer.data(msi + 1265);
    const auto *msi_1266 = buffer.data(msi + 1266);
    const auto *msi_1268 = buffer.data(msi + 1268);
    const auto *msi_1269 = buffer.data(msi + 1269);
    const auto *msi_1270 = buffer.data(msi + 1270);
    const auto *msi_1272 = buffer.data(msi + 1272);
    const auto *msi_1273 = buffer.data(msi + 1273);
    const auto *msi_1274 = buffer.data(msi + 1274);
    const auto *msi_1275 = buffer.data(msi + 1275);
    const auto *msi_1277 = buffer.data(msi + 1277);
    const auto *msi_1278 = buffer.data(msi + 1278);
    const auto *msi_1279 = buffer.data(msi + 1279);
    const auto *msi_1280 = buffer.data(msi + 1280);
    const auto *msi_1281 = buffer.data(msi + 1281);
    const auto *msi_1282 = buffer.data(msi + 1282);
    const auto *msi_1283 = buffer.data(msi + 1283);
    const auto *msi_1284 = buffer.data(msi + 1284);
    const auto *msi_1285 = buffer.data(msi + 1285);
    const auto *msi_1286 = buffer.data(msi + 1286);
    const auto *msi_1287 = buffer.data(msi + 1287);
    const auto *msi_1290 = buffer.data(msi + 1290);
    const auto *msi_1292 = buffer.data(msi + 1292);
    const auto *msi_1293 = buffer.data(msi + 1293);
    const auto *msi_1295 = buffer.data(msi + 1295);
    const auto *msi_1296 = buffer.data(msi + 1296);
    const auto *msi_1297 = buffer.data(msi + 1297);
    const auto *msi_1299 = buffer.data(msi + 1299);
    const auto *msi_1300 = buffer.data(msi + 1300);
    const auto *msi_1301 = buffer.data(msi + 1301);
    const auto *msi_1302 = buffer.data(msi + 1302);
    const auto *msi_1304 = buffer.data(msi + 1304);
    const auto *msi_1305 = buffer.data(msi + 1305);
    const auto *msi_1306 = buffer.data(msi + 1306);
    const auto *msi_1307 = buffer.data(msi + 1307);
    const auto *msi_1308 = buffer.data(msi + 1308);
    const auto *msi_1309 = buffer.data(msi + 1309);
    const auto *msi_1310 = buffer.data(msi + 1310);
    const auto *msi_1311 = buffer.data(msi + 1311);
    const auto *msi_1312 = buffer.data(msi + 1312);
    const auto *msi_1313 = buffer.data(msi + 1313);
    const auto *msi_1314 = buffer.data(msi + 1314);
    const auto *msi_1315 = buffer.data(msi + 1315);
    const auto *msi_1316 = buffer.data(msi + 1316);
    const auto *msi_1317 = buffer.data(msi + 1317);
    const auto *msi_1318 = buffer.data(msi + 1318);
    const auto *msi_1319 = buffer.data(msi + 1319);
    const auto *msi_1320 = buffer.data(msi + 1320);
    const auto *msi_1321 = buffer.data(msi + 1321);
    const auto *msi_1322 = buffer.data(msi + 1322);
    const auto *msi_1323 = buffer.data(msi + 1323);
    const auto *msi_1324 = buffer.data(msi + 1324);
    const auto *msi_1325 = buffer.data(msi + 1325);
    const auto *msi_1326 = buffer.data(msi + 1326);
    const auto *msi_1327 = buffer.data(msi + 1327);
    const auto *msi_1328 = buffer.data(msi + 1328);
    const auto *msi_1329 = buffer.data(msi + 1329);
    const auto *msi_1330 = buffer.data(msi + 1330);
    const auto *msi_1331 = buffer.data(msi + 1331);
    const auto *msi_1332 = buffer.data(msi + 1332);
    const auto *msi_1333 = buffer.data(msi + 1333);
    const auto *msi_1334 = buffer.data(msi + 1334);
    const auto *msi_1335 = buffer.data(msi + 1335);
    const auto *msi_1336 = buffer.data(msi + 1336);
    const auto *msi_1337 = buffer.data(msi + 1337);
    const auto *msi_1338 = buffer.data(msi + 1338);
    const auto *msi_1339 = buffer.data(msi + 1339);
    const auto *msi_1340 = buffer.data(msi + 1340);
    const auto *msi_1341 = buffer.data(msi + 1341);
    const auto *msi_1342 = buffer.data(msi + 1342);
    const auto *msi_1343 = buffer.data(msi + 1343);
    const auto *msi_1344 = buffer.data(msi + 1344);
    const auto *msi_1345 = buffer.data(msi + 1345);
    const auto *msi_1346 = buffer.data(msi + 1346);
    const auto *msi_1347 = buffer.data(msi + 1347);
    const auto *msi_1348 = buffer.data(msi + 1348);
    const auto *msi_1349 = buffer.data(msi + 1349);
    const auto *msi_1350 = buffer.data(msi + 1350);
    const auto *msi_1351 = buffer.data(msi + 1351);
    const auto *msi_1352 = buffer.data(msi + 1352);
    const auto *msi_1353 = buffer.data(msi + 1353);

#pragma omp simd aligned(t_1620, t_1621, t_1622, t_1623, t_1624, pc_x, pc_z, msh0_945, \
                         msh0_946, msh0_948, msh1_945, msh1_946, msh1_948, msi_1260, msi_1261, \
                         msi_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1620[k] = f_1 * msh0_945[k]
                    - f_2 * msh1_945[k]
                    + f_3 * pc_x[k] * msi_1260[k];

        t_1621[k] = f_19 * msh0_946[k]
                    - f_20 * msh1_946[k]
                    + f_3 * pc_x[k] * msi_1261[k];

        t_1622[k] = f_3 * pc_z[k] * msi_1260[k];

        t_1623[k] = f_10 * msh0_948[k]
                    - f_11 * msh1_948[k]
                    + f_3 * pc_x[k] * msi_1263[k];

        t_1624[k] = f_3 * pc_z[k] * msi_1261[k];
    }

#pragma omp simd aligned(t_1625, t_1626, t_1627, t_1628, pc_x, pc_z, msh0_950, msh0_951, \
                         msh0_953, msh1_950, msh1_951, msh1_953, msi_1263, msi_1265, msi_1266, \
                         msi_1268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = f_10 * msh0_950[k]
                    - f_11 * msh1_950[k]
                    + f_3 * pc_x[k] * msi_1265[k];

        t_1626[k] = f_8 * msh0_951[k]
                    - f_9 * msh1_951[k]
                    + f_3 * pc_x[k] * msi_1266[k];

        t_1627[k] = f_3 * pc_z[k] * msi_1263[k];

        t_1628[k] = f_8 * msh0_953[k]
                    - f_9 * msh1_953[k]
                    + f_3 * pc_x[k] * msi_1268[k];
    }

#pragma omp simd aligned(t_1629, t_1630, t_1631, t_1632, pc_x, pc_z, msh0_954, msh0_955, \
                         msh0_957, msh1_954, msh1_955, msh1_957, msi_1266, msi_1269, msi_1270, \
                         msi_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1629[k] = f_8 * msh0_954[k]
                    - f_9 * msh1_954[k]
                    + f_3 * pc_x[k] * msi_1269[k];

        t_1630[k] = f_6 * msh0_955[k]
                    - f_7 * msh1_955[k]
                    + f_3 * pc_x[k] * msi_1270[k];

        t_1631[k] = f_3 * pc_z[k] * msi_1266[k];

        t_1632[k] = f_6 * msh0_957[k]
                    - f_7 * msh1_957[k]
                    + f_3 * pc_x[k] * msi_1272[k];
    }

#pragma omp simd aligned(t_1633, t_1634, t_1635, t_1636, pc_x, pc_z, msh0_958, msh0_959, \
                         msh0_960, msh1_958, msh1_959, msh1_960, msi_1270, msi_1273, msi_1274, \
                         msi_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1633[k] = f_6 * msh0_958[k]
                    - f_7 * msh1_958[k]
                    + f_3 * pc_x[k] * msi_1273[k];

        t_1634[k] = f_6 * msh0_959[k]
                    - f_7 * msh1_959[k]
                    + f_3 * pc_x[k] * msi_1274[k];

        t_1635[k] = f_4 * msh0_960[k]
                    - f_5 * msh1_960[k]
                    + f_3 * pc_x[k] * msi_1275[k];

        t_1636[k] = f_3 * pc_z[k] * msi_1270[k];
    }

#pragma omp simd aligned(t_1637, t_1638, t_1639, pc_x, msh0_962, msh0_963, msh0_964, msh1_962, \
                         msh1_963, msh1_964, msi_1277, msi_1278, \
                         msi_1279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1637[k] = f_4 * msh0_962[k]
                    - f_5 * msh1_962[k]
                    + f_3 * pc_x[k] * msi_1277[k];

        t_1638[k] = f_4 * msh0_963[k]
                    - f_5 * msh1_963[k]
                    + f_3 * pc_x[k] * msi_1278[k];

        t_1639[k] = f_4 * msh0_964[k]
                    - f_5 * msh1_964[k]
                    + f_3 * pc_x[k] * msi_1279[k];
    }

#pragma omp simd aligned(t_1640, t_1641, t_1642, t_1643, t_1644, t_1645, pc_x, msh0_965, \
                         msh1_965, msi_1280, msi_1281, msi_1282, msi_1283, msi_1284, \
                         msi_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1640[k] = f_4 * msh0_965[k]
                    - f_5 * msh1_965[k]
                    + f_3 * pc_x[k] * msi_1280[k];

        t_1641[k] = f_3 * pc_x[k] * msi_1281[k];

        t_1642[k] = f_3 * pc_x[k] * msi_1282[k];

        t_1643[k] = f_3 * pc_x[k] * msi_1283[k];

        t_1644[k] = f_3 * pc_x[k] * msi_1284[k];

        t_1645[k] = f_3 * pc_x[k] * msi_1285[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, t_1649, t_1650, pc_x, pc_y, pc_z, lsi_1029, \
                         msh0_960, msh1_960, msi_1281, msi_1282, msi_1286, \
                         msi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = f_3 * pc_x[k] * msi_1286[k];

        t_1647[k] = f_3 * pc_x[k] * msi_1287[k];

        t_1648[k] = f_0 * lsi_1029[k]
                    + f_1 * msh0_960[k]
                    - f_2 * msh1_960[k]
                    + f_3 * pc_y[k] * msi_1281[k];

        t_1649[k] = f_3 * pc_z[k] * msi_1281[k];

        t_1650[k] = f_4 * msh0_960[k]
                    - f_5 * msh1_960[k]
                    + f_3 * pc_z[k] * msi_1282[k];
    }

#pragma omp simd aligned(t_1651, t_1652, t_1653, pc_z, msh0_961, msh0_962, msh0_963, msh1_961, \
                         msh1_962, msh1_963, msi_1283, msi_1284, \
                         msi_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1651[k] = f_6 * msh0_961[k]
                    - f_7 * msh1_961[k]
                    + f_3 * pc_z[k] * msi_1283[k];

        t_1652[k] = f_8 * msh0_962[k]
                    - f_9 * msh1_962[k]
                    + f_3 * pc_z[k] * msi_1284[k];

        t_1653[k] = f_10 * msh0_963[k]
                    - f_11 * msh1_963[k]
                    + f_3 * pc_z[k] * msi_1285[k];
    }

#pragma omp simd aligned(t_1654, t_1655, t_1656, t_1657, pa_z, pc_y, pc_z, lsk0_1296, \
                         lsk0_1297, lsi_1035, lsk1_1296, lsk1_1297, msh0_965, msh1_965, \
                         msi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1654[k] = f_0 * lsi_1035[k]
                    + f_3 * pc_y[k] * msi_1287[k];

        t_1655[k] = f_1 * msh0_965[k]
                    - f_2 * msh1_965[k]
                    + f_3 * pc_z[k] * msi_1287[k];

        t_1656[k] = pa_z[k] * lsk0_1296[k]
                    - f_12 * pc_z[k] * lsk1_1296[k];

        t_1657[k] = pa_z[k] * lsk0_1297[k]
                    - f_12 * pc_z[k] * lsk1_1297[k];
    }

#pragma omp simd aligned(t_1658, t_1659, t_1660, pa_z, pc_x, pc_z, lsk0_1299, lsk1_1299, \
                         msh0_968, msh0_970, msh1_968, msh1_970, msi_1290, \
                         msi_1292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1658[k] = f_19 * msh0_968[k]
                    - f_20 * msh1_968[k]
                    + f_3 * pc_x[k] * msi_1290[k];

        t_1659[k] = pa_z[k] * lsk0_1299[k]
                    - f_12 * pc_z[k] * lsk1_1299[k];

        t_1660[k] = f_10 * msh0_970[k]
                    - f_11 * msh1_970[k]
                    + f_3 * pc_x[k] * msi_1292[k];
    }

#pragma omp simd aligned(t_1661, t_1662, t_1663, pa_z, pc_x, pc_z, lsk0_1302, lsk1_1302, \
                         msh0_971, msh0_973, msh1_971, msh1_973, msi_1293, \
                         msi_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1661[k] = f_10 * msh0_971[k]
                    - f_11 * msh1_971[k]
                    + f_3 * pc_x[k] * msi_1293[k];

        t_1662[k] = pa_z[k] * lsk0_1302[k]
                    - f_12 * pc_z[k] * lsk1_1302[k];

        t_1663[k] = f_8 * msh0_973[k]
                    - f_9 * msh1_973[k]
                    + f_3 * pc_x[k] * msi_1295[k];
    }

#pragma omp simd aligned(t_1664, t_1665, t_1666, pa_z, pc_x, pc_z, lsk0_1306, lsk1_1306, \
                         msh0_974, msh0_975, msh1_974, msh1_975, msi_1296, \
                         msi_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1664[k] = f_8 * msh0_974[k]
                    - f_9 * msh1_974[k]
                    + f_3 * pc_x[k] * msi_1296[k];

        t_1665[k] = f_8 * msh0_975[k]
                    - f_9 * msh1_975[k]
                    + f_3 * pc_x[k] * msi_1297[k];

        t_1666[k] = pa_z[k] * lsk0_1306[k]
                    - f_12 * pc_z[k] * lsk1_1306[k];
    }

#pragma omp simd aligned(t_1667, t_1668, t_1669, pc_x, msh0_977, msh0_978, msh0_979, msh1_977, \
                         msh1_978, msh1_979, msi_1299, msi_1300, \
                         msi_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1667[k] = f_6 * msh0_977[k]
                    - f_7 * msh1_977[k]
                    + f_3 * pc_x[k] * msi_1299[k];

        t_1668[k] = f_6 * msh0_978[k]
                    - f_7 * msh1_978[k]
                    + f_3 * pc_x[k] * msi_1300[k];

        t_1669[k] = f_6 * msh0_979[k]
                    - f_7 * msh1_979[k]
                    + f_3 * pc_x[k] * msi_1301[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, pa_z, pc_x, pc_z, lsk0_1311, lsk1_1311, \
                         msh0_980, msh0_982, msh1_980, msh1_982, msi_1302, \
                         msi_1304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_6 * msh0_980[k]
                    - f_7 * msh1_980[k]
                    + f_3 * pc_x[k] * msi_1302[k];

        t_1671[k] = pa_z[k] * lsk0_1311[k]
                    - f_12 * pc_z[k] * lsk1_1311[k];

        t_1672[k] = f_4 * msh0_982[k]
                    - f_5 * msh1_982[k]
                    + f_3 * pc_x[k] * msi_1304[k];
    }

#pragma omp simd aligned(t_1673, t_1674, t_1675, pc_x, msh0_983, msh0_984, msh0_985, msh1_983, \
                         msh1_984, msh1_985, msi_1305, msi_1306, \
                         msi_1307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1673[k] = f_4 * msh0_983[k]
                    - f_5 * msh1_983[k]
                    + f_3 * pc_x[k] * msi_1305[k];

        t_1674[k] = f_4 * msh0_984[k]
                    - f_5 * msh1_984[k]
                    + f_3 * pc_x[k] * msi_1306[k];

        t_1675[k] = f_4 * msh0_985[k]
                    - f_5 * msh1_985[k]
                    + f_3 * pc_x[k] * msi_1307[k];
    }

#pragma omp simd aligned(t_1676, t_1677, t_1678, t_1679, t_1680, t_1681, pc_x, msh0_986, \
                         msh1_986, msi_1308, msi_1309, msi_1310, msi_1311, msi_1312, \
                         msi_1313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1676[k] = f_4 * msh0_986[k]
                    - f_5 * msh1_986[k]
                    + f_3 * pc_x[k] * msi_1308[k];

        t_1677[k] = f_3 * pc_x[k] * msi_1309[k];

        t_1678[k] = f_3 * pc_x[k] * msi_1310[k];

        t_1679[k] = f_3 * pc_x[k] * msi_1311[k];

        t_1680[k] = f_3 * pc_x[k] * msi_1312[k];

        t_1681[k] = f_3 * pc_x[k] * msi_1313[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, t_1685, pa_z, pc_x, pc_z, lsk0_1324, \
                         lsi_1029, lsk1_1324, msi_1309, msi_1314, \
                         msi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_3 * pc_x[k] * msi_1314[k];

        t_1683[k] = f_3 * pc_x[k] * msi_1315[k];

        t_1684[k] = pa_z[k] * lsk0_1324[k]
                    - f_12 * pc_z[k] * lsk1_1324[k];

        t_1685[k] = f_13 * lsi_1029[k]
                    + f_3 * pc_z[k] * msi_1309[k];
    }

#pragma omp simd aligned(t_1686, t_1687, t_1688, pa_z, pc_z, lsk0_1326, lsk0_1327, lsk0_1328, \
                         lsi_1030, lsi_1031, lsi_1032, lsk1_1326, lsk1_1327, \
                         lsk1_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1686[k] = pa_z[k] * lsk0_1326[k]
                    + f_14 * lsi_1030[k]
                    - f_12 * pc_z[k] * lsk1_1326[k];

        t_1687[k] = pa_z[k] * lsk0_1327[k]
                    + f_15 * lsi_1031[k]
                    - f_12 * pc_z[k] * lsk1_1327[k];

        t_1688[k] = pa_z[k] * lsk0_1328[k]
                    + f_16 * lsi_1032[k]
                    - f_12 * pc_z[k] * lsk1_1328[k];
    }

#pragma omp simd aligned(t_1689, t_1690, t_1691, pa_z, pc_y, pc_z, lsk0_1329, lsi_1033, \
                         lsi_1035, lsi_1063, lsk1_1329, msh0_986, msh1_986, \
                         msi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1689[k] = pa_z[k] * lsk0_1329[k]
                    + f_17 * lsi_1033[k]
                    - f_12 * pc_z[k] * lsk1_1329[k];

        t_1690[k] = f_18 * lsi_1063[k]
                    + f_3 * pc_y[k] * msi_1315[k];

        t_1691[k] = f_13 * lsi_1035[k]
                    + f_1 * msh0_986[k]
                    - f_2 * msh1_986[k]
                    + f_3 * pc_z[k] * msi_1315[k];
    }

#pragma omp simd aligned(t_1692, t_1693, t_1694, pc_x, msh0_987, msh0_988, msh0_989, msh1_987, \
                         msh1_988, msh1_989, msi_1316, msi_1317, \
                         msi_1318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1692[k] = f_1 * msh0_987[k]
                    - f_2 * msh1_987[k]
                    + f_3 * pc_x[k] * msi_1316[k];

        t_1693[k] = f_19 * msh0_988[k]
                    - f_20 * msh1_988[k]
                    + f_3 * pc_x[k] * msi_1317[k];

        t_1694[k] = f_19 * msh0_989[k]
                    - f_20 * msh1_989[k]
                    + f_3 * pc_x[k] * msi_1318[k];
    }

#pragma omp simd aligned(t_1695, t_1696, t_1697, pc_x, msh0_990, msh0_991, msh0_992, msh1_990, \
                         msh1_991, msh1_992, msi_1319, msi_1320, \
                         msi_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1695[k] = f_10 * msh0_990[k]
                    - f_11 * msh1_990[k]
                    + f_3 * pc_x[k] * msi_1319[k];

        t_1696[k] = f_10 * msh0_991[k]
                    - f_11 * msh1_991[k]
                    + f_3 * pc_x[k] * msi_1320[k];

        t_1697[k] = f_10 * msh0_992[k]
                    - f_11 * msh1_992[k]
                    + f_3 * pc_x[k] * msi_1321[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, pc_x, msh0_993, msh0_994, msh0_995, msh1_993, \
                         msh1_994, msh1_995, msi_1322, msi_1323, \
                         msi_1324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_8 * msh0_993[k]
                    - f_9 * msh1_993[k]
                    + f_3 * pc_x[k] * msi_1322[k];

        t_1699[k] = f_8 * msh0_994[k]
                    - f_9 * msh1_994[k]
                    + f_3 * pc_x[k] * msi_1323[k];

        t_1700[k] = f_8 * msh0_995[k]
                    - f_9 * msh1_995[k]
                    + f_3 * pc_x[k] * msi_1324[k];
    }

#pragma omp simd aligned(t_1701, t_1702, t_1703, pc_x, msh0_996, msh0_997, msh0_998, msh1_996, \
                         msh1_997, msh1_998, msi_1325, msi_1326, \
                         msi_1327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1701[k] = f_8 * msh0_996[k]
                    - f_9 * msh1_996[k]
                    + f_3 * pc_x[k] * msi_1325[k];

        t_1702[k] = f_6 * msh0_997[k]
                    - f_7 * msh1_997[k]
                    + f_3 * pc_x[k] * msi_1326[k];

        t_1703[k] = f_6 * msh0_998[k]
                    - f_7 * msh1_998[k]
                    + f_3 * pc_x[k] * msi_1327[k];
    }

#pragma omp simd aligned(t_1704, t_1705, t_1706, pc_x, msh0_999, msh0_1000, msh0_1001, \
                         msh1_999, msh1_1000, msh1_1001, msi_1328, msi_1329, \
                         msi_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1704[k] = f_6 * msh0_999[k]
                    - f_7 * msh1_999[k]
                    + f_3 * pc_x[k] * msi_1328[k];

        t_1705[k] = f_6 * msh0_1000[k]
                    - f_7 * msh1_1000[k]
                    + f_3 * pc_x[k] * msi_1329[k];

        t_1706[k] = f_6 * msh0_1001[k]
                    - f_7 * msh1_1001[k]
                    + f_3 * pc_x[k] * msi_1330[k];
    }

#pragma omp simd aligned(t_1707, t_1708, t_1709, pc_x, msh0_1002, msh0_1003, msh0_1004, \
                         msh1_1002, msh1_1003, msh1_1004, msi_1331, msi_1332, \
                         msi_1333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1707[k] = f_4 * msh0_1002[k]
                    - f_5 * msh1_1002[k]
                    + f_3 * pc_x[k] * msi_1331[k];

        t_1708[k] = f_4 * msh0_1003[k]
                    - f_5 * msh1_1003[k]
                    + f_3 * pc_x[k] * msi_1332[k];

        t_1709[k] = f_4 * msh0_1004[k]
                    - f_5 * msh1_1004[k]
                    + f_3 * pc_x[k] * msi_1333[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, t_1713, pc_x, msh0_1005, msh0_1006, \
                         msh0_1007, msh1_1005, msh1_1006, msh1_1007, msi_1334, msi_1335, \
                         msi_1336, msi_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = f_4 * msh0_1005[k]
                    - f_5 * msh1_1005[k]
                    + f_3 * pc_x[k] * msi_1334[k];

        t_1711[k] = f_4 * msh0_1006[k]
                    - f_5 * msh1_1006[k]
                    + f_3 * pc_x[k] * msi_1335[k];

        t_1712[k] = f_4 * msh0_1007[k]
                    - f_5 * msh1_1007[k]
                    + f_3 * pc_x[k] * msi_1336[k];

        t_1713[k] = f_3 * pc_x[k] * msi_1337[k];
    }

#pragma omp simd aligned(t_1714, t_1715, t_1716, t_1717, t_1718, t_1719, pc_x, msi_1338, \
                         msi_1339, msi_1340, msi_1341, msi_1342, \
                         msi_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = f_3 * pc_x[k] * msi_1338[k];

        t_1715[k] = f_3 * pc_x[k] * msi_1339[k];

        t_1716[k] = f_3 * pc_x[k] * msi_1340[k];

        t_1717[k] = f_3 * pc_x[k] * msi_1341[k];

        t_1718[k] = f_3 * pc_x[k] * msi_1342[k];

        t_1719[k] = f_3 * pc_x[k] * msi_1343[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, pc_y, pc_z, lsi_1057, lsi_1085, lsi_1087, \
                         msh0_1002, msh0_1004, msh1_1002, msh1_1004, msi_1337, \
                         msi_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = f_21 * lsi_1085[k]
                    + f_1 * msh0_1002[k]
                    - f_2 * msh1_1002[k]
                    + f_3 * pc_y[k] * msi_1337[k];

        t_1721[k] = f_14 * lsi_1057[k]
                    + f_3 * pc_z[k] * msi_1337[k];

        t_1722[k] = f_21 * lsi_1087[k]
                    + f_10 * msh0_1004[k]
                    - f_11 * msh1_1004[k]
                    + f_3 * pc_y[k] * msi_1339[k];
    }

#pragma omp simd aligned(t_1723, t_1724, t_1725, pc_y, lsi_1088, lsi_1089, lsi_1090, \
                         msh0_1005, msh0_1006, msh0_1007, msh1_1005, msh1_1006, msh1_1007, \
                         msi_1340, msi_1341, msi_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1723[k] = f_21 * lsi_1088[k]
                    + f_8 * msh0_1005[k]
                    - f_9 * msh1_1005[k]
                    + f_3 * pc_y[k] * msi_1340[k];

        t_1724[k] = f_21 * lsi_1089[k]
                    + f_6 * msh0_1006[k]
                    - f_7 * msh1_1006[k]
                    + f_3 * pc_y[k] * msi_1341[k];

        t_1725[k] = f_21 * lsi_1090[k]
                    + f_4 * msh0_1007[k]
                    - f_5 * msh1_1007[k]
                    + f_3 * pc_y[k] * msi_1342[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, pc_x, pc_y, pc_z, lsi_1063, lsi_1091, \
                         msh0_1007, msh0_1008, msh1_1007, msh1_1008, msi_1343, \
                         msi_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_21 * lsi_1091[k]
                    + f_3 * pc_y[k] * msi_1343[k];

        t_1727[k] = f_14 * lsi_1063[k]
                    + f_1 * msh0_1007[k]
                    - f_2 * msh1_1007[k]
                    + f_3 * pc_z[k] * msi_1343[k];

        t_1728[k] = f_1 * msh0_1008[k]
                    - f_2 * msh1_1008[k]
                    + f_3 * pc_x[k] * msi_1344[k];
    }

#pragma omp simd aligned(t_1729, t_1730, t_1731, pc_x, msh0_1009, msh0_1010, msh0_1011, \
                         msh1_1009, msh1_1010, msh1_1011, msi_1345, msi_1346, \
                         msi_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1729[k] = f_19 * msh0_1009[k]
                    - f_20 * msh1_1009[k]
                    + f_3 * pc_x[k] * msi_1345[k];

        t_1730[k] = f_19 * msh0_1010[k]
                    - f_20 * msh1_1010[k]
                    + f_3 * pc_x[k] * msi_1346[k];

        t_1731[k] = f_10 * msh0_1011[k]
                    - f_11 * msh1_1011[k]
                    + f_3 * pc_x[k] * msi_1347[k];
    }

#pragma omp simd aligned(t_1732, t_1733, t_1734, pc_x, msh0_1012, msh0_1013, msh0_1014, \
                         msh1_1012, msh1_1013, msh1_1014, msi_1348, msi_1349, \
                         msi_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1732[k] = f_10 * msh0_1012[k]
                    - f_11 * msh1_1012[k]
                    + f_3 * pc_x[k] * msi_1348[k];

        t_1733[k] = f_10 * msh0_1013[k]
                    - f_11 * msh1_1013[k]
                    + f_3 * pc_x[k] * msi_1349[k];

        t_1734[k] = f_8 * msh0_1014[k]
                    - f_9 * msh1_1014[k]
                    + f_3 * pc_x[k] * msi_1350[k];
    }

#pragma omp simd aligned(t_1735, t_1736, t_1737, pc_x, msh0_1015, msh0_1016, msh0_1017, \
                         msh1_1015, msh1_1016, msh1_1017, msi_1351, msi_1352, \
                         msi_1353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1735[k] = f_8 * msh0_1015[k]
                    - f_9 * msh1_1015[k]
                    + f_3 * pc_x[k] * msi_1351[k];

        t_1736[k] = f_8 * msh0_1016[k]
                    - f_9 * msh1_1016[k]
                    + f_3 * pc_x[k] * msi_1352[k];

        t_1737[k] = f_8 * msh0_1017[k]
                    - f_9 * msh1_1017[k]
                    + f_3 * pc_x[k] * msi_1353[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t lsi, const size_t msh0,
                                                           const size_t msh1, const size_t msi,
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
    const auto f_22 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi_1085 = buffer.data(lsi + 1085);
    const auto *lsi_1091 = buffer.data(lsi + 1091);
    const auto *lsi_1113 = buffer.data(lsi + 1113);
    const auto *lsi_1115 = buffer.data(lsi + 1115);
    const auto *lsi_1116 = buffer.data(lsi + 1116);
    const auto *lsi_1117 = buffer.data(lsi + 1117);
    const auto *lsi_1118 = buffer.data(lsi + 1118);
    const auto *lsi_1119 = buffer.data(lsi + 1119);
    const auto *lsi_1141 = buffer.data(lsi + 1141);
    const auto *lsi_1143 = buffer.data(lsi + 1143);
    const auto *lsi_1144 = buffer.data(lsi + 1144);
    const auto *lsi_1145 = buffer.data(lsi + 1145);
    const auto *lsi_1146 = buffer.data(lsi + 1146);
    const auto *lsi_1147 = buffer.data(lsi + 1147);
    const auto *lsi_1169 = buffer.data(lsi + 1169);
    const auto *lsi_1171 = buffer.data(lsi + 1171);
    const auto *lsi_1172 = buffer.data(lsi + 1172);
    const auto *lsi_1173 = buffer.data(lsi + 1173);
    const auto *lsi_1174 = buffer.data(lsi + 1174);
    const auto *lsi_1175 = buffer.data(lsi + 1175);

    const auto *msh0_1018 = buffer.data(msh0 + 1018);
    const auto *msh0_1019 = buffer.data(msh0 + 1019);
    const auto *msh0_1020 = buffer.data(msh0 + 1020);
    const auto *msh0_1021 = buffer.data(msh0 + 1021);
    const auto *msh0_1022 = buffer.data(msh0 + 1022);
    const auto *msh0_1023 = buffer.data(msh0 + 1023);
    const auto *msh0_1024 = buffer.data(msh0 + 1024);
    const auto *msh0_1025 = buffer.data(msh0 + 1025);
    const auto *msh0_1026 = buffer.data(msh0 + 1026);
    const auto *msh0_1027 = buffer.data(msh0 + 1027);
    const auto *msh0_1028 = buffer.data(msh0 + 1028);
    const auto *msh0_1029 = buffer.data(msh0 + 1029);
    const auto *msh0_1030 = buffer.data(msh0 + 1030);
    const auto *msh0_1031 = buffer.data(msh0 + 1031);
    const auto *msh0_1032 = buffer.data(msh0 + 1032);
    const auto *msh0_1033 = buffer.data(msh0 + 1033);
    const auto *msh0_1034 = buffer.data(msh0 + 1034);
    const auto *msh0_1035 = buffer.data(msh0 + 1035);
    const auto *msh0_1036 = buffer.data(msh0 + 1036);
    const auto *msh0_1037 = buffer.data(msh0 + 1037);
    const auto *msh0_1038 = buffer.data(msh0 + 1038);
    const auto *msh0_1039 = buffer.data(msh0 + 1039);
    const auto *msh0_1040 = buffer.data(msh0 + 1040);
    const auto *msh0_1041 = buffer.data(msh0 + 1041);
    const auto *msh0_1042 = buffer.data(msh0 + 1042);
    const auto *msh0_1043 = buffer.data(msh0 + 1043);
    const auto *msh0_1044 = buffer.data(msh0 + 1044);
    const auto *msh0_1045 = buffer.data(msh0 + 1045);
    const auto *msh0_1046 = buffer.data(msh0 + 1046);
    const auto *msh0_1047 = buffer.data(msh0 + 1047);
    const auto *msh0_1048 = buffer.data(msh0 + 1048);
    const auto *msh0_1049 = buffer.data(msh0 + 1049);
    const auto *msh0_1050 = buffer.data(msh0 + 1050);
    const auto *msh0_1051 = buffer.data(msh0 + 1051);
    const auto *msh0_1052 = buffer.data(msh0 + 1052);
    const auto *msh0_1053 = buffer.data(msh0 + 1053);
    const auto *msh0_1054 = buffer.data(msh0 + 1054);
    const auto *msh0_1055 = buffer.data(msh0 + 1055);
    const auto *msh0_1056 = buffer.data(msh0 + 1056);
    const auto *msh0_1057 = buffer.data(msh0 + 1057);
    const auto *msh0_1058 = buffer.data(msh0 + 1058);
    const auto *msh0_1059 = buffer.data(msh0 + 1059);
    const auto *msh0_1060 = buffer.data(msh0 + 1060);
    const auto *msh0_1061 = buffer.data(msh0 + 1061);
    const auto *msh0_1062 = buffer.data(msh0 + 1062);
    const auto *msh0_1063 = buffer.data(msh0 + 1063);
    const auto *msh0_1064 = buffer.data(msh0 + 1064);
    const auto *msh0_1065 = buffer.data(msh0 + 1065);
    const auto *msh0_1066 = buffer.data(msh0 + 1066);
    const auto *msh0_1067 = buffer.data(msh0 + 1067);
    const auto *msh0_1068 = buffer.data(msh0 + 1068);
    const auto *msh0_1069 = buffer.data(msh0 + 1069);
    const auto *msh0_1070 = buffer.data(msh0 + 1070);
    const auto *msh0_1071 = buffer.data(msh0 + 1071);
    const auto *msh0_1072 = buffer.data(msh0 + 1072);
    const auto *msh0_1073 = buffer.data(msh0 + 1073);
    const auto *msh0_1074 = buffer.data(msh0 + 1074);
    const auto *msh0_1075 = buffer.data(msh0 + 1075);
    const auto *msh0_1076 = buffer.data(msh0 + 1076);
    const auto *msh0_1077 = buffer.data(msh0 + 1077);
    const auto *msh0_1078 = buffer.data(msh0 + 1078);
    const auto *msh0_1079 = buffer.data(msh0 + 1079);
    const auto *msh0_1080 = buffer.data(msh0 + 1080);
    const auto *msh0_1081 = buffer.data(msh0 + 1081);
    const auto *msh0_1082 = buffer.data(msh0 + 1082);
    const auto *msh0_1083 = buffer.data(msh0 + 1083);
    const auto *msh0_1084 = buffer.data(msh0 + 1084);
    const auto *msh0_1085 = buffer.data(msh0 + 1085);
    const auto *msh0_1086 = buffer.data(msh0 + 1086);

    const auto *msh1_1018 = buffer.data(msh1 + 1018);
    const auto *msh1_1019 = buffer.data(msh1 + 1019);
    const auto *msh1_1020 = buffer.data(msh1 + 1020);
    const auto *msh1_1021 = buffer.data(msh1 + 1021);
    const auto *msh1_1022 = buffer.data(msh1 + 1022);
    const auto *msh1_1023 = buffer.data(msh1 + 1023);
    const auto *msh1_1024 = buffer.data(msh1 + 1024);
    const auto *msh1_1025 = buffer.data(msh1 + 1025);
    const auto *msh1_1026 = buffer.data(msh1 + 1026);
    const auto *msh1_1027 = buffer.data(msh1 + 1027);
    const auto *msh1_1028 = buffer.data(msh1 + 1028);
    const auto *msh1_1029 = buffer.data(msh1 + 1029);
    const auto *msh1_1030 = buffer.data(msh1 + 1030);
    const auto *msh1_1031 = buffer.data(msh1 + 1031);
    const auto *msh1_1032 = buffer.data(msh1 + 1032);
    const auto *msh1_1033 = buffer.data(msh1 + 1033);
    const auto *msh1_1034 = buffer.data(msh1 + 1034);
    const auto *msh1_1035 = buffer.data(msh1 + 1035);
    const auto *msh1_1036 = buffer.data(msh1 + 1036);
    const auto *msh1_1037 = buffer.data(msh1 + 1037);
    const auto *msh1_1038 = buffer.data(msh1 + 1038);
    const auto *msh1_1039 = buffer.data(msh1 + 1039);
    const auto *msh1_1040 = buffer.data(msh1 + 1040);
    const auto *msh1_1041 = buffer.data(msh1 + 1041);
    const auto *msh1_1042 = buffer.data(msh1 + 1042);
    const auto *msh1_1043 = buffer.data(msh1 + 1043);
    const auto *msh1_1044 = buffer.data(msh1 + 1044);
    const auto *msh1_1045 = buffer.data(msh1 + 1045);
    const auto *msh1_1046 = buffer.data(msh1 + 1046);
    const auto *msh1_1047 = buffer.data(msh1 + 1047);
    const auto *msh1_1048 = buffer.data(msh1 + 1048);
    const auto *msh1_1049 = buffer.data(msh1 + 1049);
    const auto *msh1_1050 = buffer.data(msh1 + 1050);
    const auto *msh1_1051 = buffer.data(msh1 + 1051);
    const auto *msh1_1052 = buffer.data(msh1 + 1052);
    const auto *msh1_1053 = buffer.data(msh1 + 1053);
    const auto *msh1_1054 = buffer.data(msh1 + 1054);
    const auto *msh1_1055 = buffer.data(msh1 + 1055);
    const auto *msh1_1056 = buffer.data(msh1 + 1056);
    const auto *msh1_1057 = buffer.data(msh1 + 1057);
    const auto *msh1_1058 = buffer.data(msh1 + 1058);
    const auto *msh1_1059 = buffer.data(msh1 + 1059);
    const auto *msh1_1060 = buffer.data(msh1 + 1060);
    const auto *msh1_1061 = buffer.data(msh1 + 1061);
    const auto *msh1_1062 = buffer.data(msh1 + 1062);
    const auto *msh1_1063 = buffer.data(msh1 + 1063);
    const auto *msh1_1064 = buffer.data(msh1 + 1064);
    const auto *msh1_1065 = buffer.data(msh1 + 1065);
    const auto *msh1_1066 = buffer.data(msh1 + 1066);
    const auto *msh1_1067 = buffer.data(msh1 + 1067);
    const auto *msh1_1068 = buffer.data(msh1 + 1068);
    const auto *msh1_1069 = buffer.data(msh1 + 1069);
    const auto *msh1_1070 = buffer.data(msh1 + 1070);
    const auto *msh1_1071 = buffer.data(msh1 + 1071);
    const auto *msh1_1072 = buffer.data(msh1 + 1072);
    const auto *msh1_1073 = buffer.data(msh1 + 1073);
    const auto *msh1_1074 = buffer.data(msh1 + 1074);
    const auto *msh1_1075 = buffer.data(msh1 + 1075);
    const auto *msh1_1076 = buffer.data(msh1 + 1076);
    const auto *msh1_1077 = buffer.data(msh1 + 1077);
    const auto *msh1_1078 = buffer.data(msh1 + 1078);
    const auto *msh1_1079 = buffer.data(msh1 + 1079);
    const auto *msh1_1080 = buffer.data(msh1 + 1080);
    const auto *msh1_1081 = buffer.data(msh1 + 1081);
    const auto *msh1_1082 = buffer.data(msh1 + 1082);
    const auto *msh1_1083 = buffer.data(msh1 + 1083);
    const auto *msh1_1084 = buffer.data(msh1 + 1084);
    const auto *msh1_1085 = buffer.data(msh1 + 1085);
    const auto *msh1_1086 = buffer.data(msh1 + 1086);

    const auto *msi_1354 = buffer.data(msi + 1354);
    const auto *msi_1355 = buffer.data(msi + 1355);
    const auto *msi_1356 = buffer.data(msi + 1356);
    const auto *msi_1357 = buffer.data(msi + 1357);
    const auto *msi_1358 = buffer.data(msi + 1358);
    const auto *msi_1359 = buffer.data(msi + 1359);
    const auto *msi_1360 = buffer.data(msi + 1360);
    const auto *msi_1361 = buffer.data(msi + 1361);
    const auto *msi_1362 = buffer.data(msi + 1362);
    const auto *msi_1363 = buffer.data(msi + 1363);
    const auto *msi_1364 = buffer.data(msi + 1364);
    const auto *msi_1365 = buffer.data(msi + 1365);
    const auto *msi_1366 = buffer.data(msi + 1366);
    const auto *msi_1367 = buffer.data(msi + 1367);
    const auto *msi_1368 = buffer.data(msi + 1368);
    const auto *msi_1369 = buffer.data(msi + 1369);
    const auto *msi_1370 = buffer.data(msi + 1370);
    const auto *msi_1371 = buffer.data(msi + 1371);
    const auto *msi_1372 = buffer.data(msi + 1372);
    const auto *msi_1373 = buffer.data(msi + 1373);
    const auto *msi_1374 = buffer.data(msi + 1374);
    const auto *msi_1375 = buffer.data(msi + 1375);
    const auto *msi_1376 = buffer.data(msi + 1376);
    const auto *msi_1377 = buffer.data(msi + 1377);
    const auto *msi_1378 = buffer.data(msi + 1378);
    const auto *msi_1379 = buffer.data(msi + 1379);
    const auto *msi_1380 = buffer.data(msi + 1380);
    const auto *msi_1381 = buffer.data(msi + 1381);
    const auto *msi_1382 = buffer.data(msi + 1382);
    const auto *msi_1383 = buffer.data(msi + 1383);
    const auto *msi_1384 = buffer.data(msi + 1384);
    const auto *msi_1385 = buffer.data(msi + 1385);
    const auto *msi_1386 = buffer.data(msi + 1386);
    const auto *msi_1387 = buffer.data(msi + 1387);
    const auto *msi_1388 = buffer.data(msi + 1388);
    const auto *msi_1389 = buffer.data(msi + 1389);
    const auto *msi_1390 = buffer.data(msi + 1390);
    const auto *msi_1391 = buffer.data(msi + 1391);
    const auto *msi_1392 = buffer.data(msi + 1392);
    const auto *msi_1393 = buffer.data(msi + 1393);
    const auto *msi_1394 = buffer.data(msi + 1394);
    const auto *msi_1395 = buffer.data(msi + 1395);
    const auto *msi_1396 = buffer.data(msi + 1396);
    const auto *msi_1397 = buffer.data(msi + 1397);
    const auto *msi_1398 = buffer.data(msi + 1398);
    const auto *msi_1399 = buffer.data(msi + 1399);
    const auto *msi_1400 = buffer.data(msi + 1400);
    const auto *msi_1401 = buffer.data(msi + 1401);
    const auto *msi_1402 = buffer.data(msi + 1402);
    const auto *msi_1403 = buffer.data(msi + 1403);
    const auto *msi_1404 = buffer.data(msi + 1404);
    const auto *msi_1405 = buffer.data(msi + 1405);
    const auto *msi_1406 = buffer.data(msi + 1406);
    const auto *msi_1407 = buffer.data(msi + 1407);
    const auto *msi_1408 = buffer.data(msi + 1408);
    const auto *msi_1409 = buffer.data(msi + 1409);
    const auto *msi_1410 = buffer.data(msi + 1410);
    const auto *msi_1411 = buffer.data(msi + 1411);
    const auto *msi_1412 = buffer.data(msi + 1412);
    const auto *msi_1413 = buffer.data(msi + 1413);
    const auto *msi_1414 = buffer.data(msi + 1414);
    const auto *msi_1415 = buffer.data(msi + 1415);
    const auto *msi_1416 = buffer.data(msi + 1416);
    const auto *msi_1417 = buffer.data(msi + 1417);
    const auto *msi_1418 = buffer.data(msi + 1418);
    const auto *msi_1419 = buffer.data(msi + 1419);
    const auto *msi_1420 = buffer.data(msi + 1420);
    const auto *msi_1421 = buffer.data(msi + 1421);
    const auto *msi_1422 = buffer.data(msi + 1422);
    const auto *msi_1423 = buffer.data(msi + 1423);
    const auto *msi_1424 = buffer.data(msi + 1424);
    const auto *msi_1425 = buffer.data(msi + 1425);
    const auto *msi_1426 = buffer.data(msi + 1426);
    const auto *msi_1427 = buffer.data(msi + 1427);
    const auto *msi_1428 = buffer.data(msi + 1428);
    const auto *msi_1429 = buffer.data(msi + 1429);
    const auto *msi_1430 = buffer.data(msi + 1430);
    const auto *msi_1431 = buffer.data(msi + 1431);
    const auto *msi_1432 = buffer.data(msi + 1432);
    const auto *msi_1433 = buffer.data(msi + 1433);
    const auto *msi_1434 = buffer.data(msi + 1434);
    const auto *msi_1435 = buffer.data(msi + 1435);
    const auto *msi_1436 = buffer.data(msi + 1436);
    const auto *msi_1437 = buffer.data(msi + 1437);
    const auto *msi_1438 = buffer.data(msi + 1438);
    const auto *msi_1439 = buffer.data(msi + 1439);
    const auto *msi_1440 = buffer.data(msi + 1440);
    const auto *msi_1441 = buffer.data(msi + 1441);
    const auto *msi_1442 = buffer.data(msi + 1442);
    const auto *msi_1443 = buffer.data(msi + 1443);

#pragma omp simd aligned(t_1738, t_1739, t_1740, pc_x, msh0_1018, msh0_1019, msh0_1020, \
                         msh1_1018, msh1_1019, msh1_1020, msi_1354, msi_1355, \
                         msi_1356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1738[k] = f_6 * msh0_1018[k]
                    - f_7 * msh1_1018[k]
                    + f_3 * pc_x[k] * msi_1354[k];

        t_1739[k] = f_6 * msh0_1019[k]
                    - f_7 * msh1_1019[k]
                    + f_3 * pc_x[k] * msi_1355[k];

        t_1740[k] = f_6 * msh0_1020[k]
                    - f_7 * msh1_1020[k]
                    + f_3 * pc_x[k] * msi_1356[k];
    }

#pragma omp simd aligned(t_1741, t_1742, t_1743, pc_x, msh0_1021, msh0_1022, msh0_1023, \
                         msh1_1021, msh1_1022, msh1_1023, msi_1357, msi_1358, \
                         msi_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1741[k] = f_6 * msh0_1021[k]
                    - f_7 * msh1_1021[k]
                    + f_3 * pc_x[k] * msi_1357[k];

        t_1742[k] = f_6 * msh0_1022[k]
                    - f_7 * msh1_1022[k]
                    + f_3 * pc_x[k] * msi_1358[k];

        t_1743[k] = f_4 * msh0_1023[k]
                    - f_5 * msh1_1023[k]
                    + f_3 * pc_x[k] * msi_1359[k];
    }

#pragma omp simd aligned(t_1744, t_1745, t_1746, pc_x, msh0_1024, msh0_1025, msh0_1026, \
                         msh1_1024, msh1_1025, msh1_1026, msi_1360, msi_1361, \
                         msi_1362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1744[k] = f_4 * msh0_1024[k]
                    - f_5 * msh1_1024[k]
                    + f_3 * pc_x[k] * msi_1360[k];

        t_1745[k] = f_4 * msh0_1025[k]
                    - f_5 * msh1_1025[k]
                    + f_3 * pc_x[k] * msi_1361[k];

        t_1746[k] = f_4 * msh0_1026[k]
                    - f_5 * msh1_1026[k]
                    + f_3 * pc_x[k] * msi_1362[k];
    }

#pragma omp simd aligned(t_1747, t_1748, t_1749, t_1750, t_1751, pc_x, msh0_1027, msh0_1028, \
                         msh1_1027, msh1_1028, msi_1363, msi_1364, msi_1365, msi_1366, \
                         msi_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1747[k] = f_4 * msh0_1027[k]
                    - f_5 * msh1_1027[k]
                    + f_3 * pc_x[k] * msi_1363[k];

        t_1748[k] = f_4 * msh0_1028[k]
                    - f_5 * msh1_1028[k]
                    + f_3 * pc_x[k] * msi_1364[k];

        t_1749[k] = f_3 * pc_x[k] * msi_1365[k];

        t_1750[k] = f_3 * pc_x[k] * msi_1366[k];

        t_1751[k] = f_3 * pc_x[k] * msi_1367[k];
    }

#pragma omp simd aligned(t_1752, t_1753, t_1754, t_1755, t_1756, pc_x, pc_y, lsi_1113, \
                         msh0_1023, msh1_1023, msi_1365, msi_1368, msi_1369, msi_1370, \
                         msi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1752[k] = f_3 * pc_x[k] * msi_1368[k];

        t_1753[k] = f_3 * pc_x[k] * msi_1369[k];

        t_1754[k] = f_3 * pc_x[k] * msi_1370[k];

        t_1755[k] = f_3 * pc_x[k] * msi_1371[k];

        t_1756[k] = f_22 * lsi_1113[k]
                    + f_1 * msh0_1023[k]
                    - f_2 * msh1_1023[k]
                    + f_3 * pc_y[k] * msi_1365[k];
    }

#pragma omp simd aligned(t_1757, t_1758, t_1759, pc_y, pc_z, lsi_1085, lsi_1115, lsi_1116, \
                         msh0_1025, msh0_1026, msh1_1025, msh1_1026, msi_1365, msi_1367, \
                         msi_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1757[k] = f_15 * lsi_1085[k]
                    + f_3 * pc_z[k] * msi_1365[k];

        t_1758[k] = f_22 * lsi_1115[k]
                    + f_10 * msh0_1025[k]
                    - f_11 * msh1_1025[k]
                    + f_3 * pc_y[k] * msi_1367[k];

        t_1759[k] = f_22 * lsi_1116[k]
                    + f_8 * msh0_1026[k]
                    - f_9 * msh1_1026[k]
                    + f_3 * pc_y[k] * msi_1368[k];
    }

#pragma omp simd aligned(t_1760, t_1761, t_1762, pc_y, lsi_1117, lsi_1118, lsi_1119, \
                         msh0_1027, msh0_1028, msh1_1027, msh1_1028, msi_1369, msi_1370, \
                         msi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1760[k] = f_22 * lsi_1117[k]
                    + f_6 * msh0_1027[k]
                    - f_7 * msh1_1027[k]
                    + f_3 * pc_y[k] * msi_1369[k];

        t_1761[k] = f_22 * lsi_1118[k]
                    + f_4 * msh0_1028[k]
                    - f_5 * msh1_1028[k]
                    + f_3 * pc_y[k] * msi_1370[k];

        t_1762[k] = f_22 * lsi_1119[k]
                    + f_3 * pc_y[k] * msi_1371[k];
    }

#pragma omp simd aligned(t_1763, t_1764, t_1765, pc_x, pc_z, lsi_1091, msh0_1028, msh0_1029, \
                         msh0_1030, msh1_1028, msh1_1029, msh1_1030, msi_1371, msi_1372, \
                         msi_1373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1763[k] = f_15 * lsi_1091[k]
                    + f_1 * msh0_1028[k]
                    - f_2 * msh1_1028[k]
                    + f_3 * pc_z[k] * msi_1371[k];

        t_1764[k] = f_1 * msh0_1029[k]
                    - f_2 * msh1_1029[k]
                    + f_3 * pc_x[k] * msi_1372[k];

        t_1765[k] = f_19 * msh0_1030[k]
                    - f_20 * msh1_1030[k]
                    + f_3 * pc_x[k] * msi_1373[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pc_x, msh0_1031, msh0_1032, msh0_1033, \
                         msh1_1031, msh1_1032, msh1_1033, msi_1374, msi_1375, \
                         msi_1376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_19 * msh0_1031[k]
                    - f_20 * msh1_1031[k]
                    + f_3 * pc_x[k] * msi_1374[k];

        t_1767[k] = f_10 * msh0_1032[k]
                    - f_11 * msh1_1032[k]
                    + f_3 * pc_x[k] * msi_1375[k];

        t_1768[k] = f_10 * msh0_1033[k]
                    - f_11 * msh1_1033[k]
                    + f_3 * pc_x[k] * msi_1376[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, pc_x, msh0_1034, msh0_1035, msh0_1036, \
                         msh1_1034, msh1_1035, msh1_1036, msi_1377, msi_1378, \
                         msi_1379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = f_10 * msh0_1034[k]
                    - f_11 * msh1_1034[k]
                    + f_3 * pc_x[k] * msi_1377[k];

        t_1770[k] = f_8 * msh0_1035[k]
                    - f_9 * msh1_1035[k]
                    + f_3 * pc_x[k] * msi_1378[k];

        t_1771[k] = f_8 * msh0_1036[k]
                    - f_9 * msh1_1036[k]
                    + f_3 * pc_x[k] * msi_1379[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, pc_x, msh0_1037, msh0_1038, msh0_1039, \
                         msh1_1037, msh1_1038, msh1_1039, msi_1380, msi_1381, \
                         msi_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = f_8 * msh0_1037[k]
                    - f_9 * msh1_1037[k]
                    + f_3 * pc_x[k] * msi_1380[k];

        t_1773[k] = f_8 * msh0_1038[k]
                    - f_9 * msh1_1038[k]
                    + f_3 * pc_x[k] * msi_1381[k];

        t_1774[k] = f_6 * msh0_1039[k]
                    - f_7 * msh1_1039[k]
                    + f_3 * pc_x[k] * msi_1382[k];
    }

#pragma omp simd aligned(t_1775, t_1776, t_1777, pc_x, msh0_1040, msh0_1041, msh0_1042, \
                         msh1_1040, msh1_1041, msh1_1042, msi_1383, msi_1384, \
                         msi_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1775[k] = f_6 * msh0_1040[k]
                    - f_7 * msh1_1040[k]
                    + f_3 * pc_x[k] * msi_1383[k];

        t_1776[k] = f_6 * msh0_1041[k]
                    - f_7 * msh1_1041[k]
                    + f_3 * pc_x[k] * msi_1384[k];

        t_1777[k] = f_6 * msh0_1042[k]
                    - f_7 * msh1_1042[k]
                    + f_3 * pc_x[k] * msi_1385[k];
    }

#pragma omp simd aligned(t_1778, t_1779, t_1780, pc_x, msh0_1043, msh0_1044, msh0_1045, \
                         msh1_1043, msh1_1044, msh1_1045, msi_1386, msi_1387, \
                         msi_1388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1778[k] = f_6 * msh0_1043[k]
                    - f_7 * msh1_1043[k]
                    + f_3 * pc_x[k] * msi_1386[k];

        t_1779[k] = f_4 * msh0_1044[k]
                    - f_5 * msh1_1044[k]
                    + f_3 * pc_x[k] * msi_1387[k];

        t_1780[k] = f_4 * msh0_1045[k]
                    - f_5 * msh1_1045[k]
                    + f_3 * pc_x[k] * msi_1388[k];
    }

#pragma omp simd aligned(t_1781, t_1782, t_1783, pc_x, msh0_1046, msh0_1047, msh0_1048, \
                         msh1_1046, msh1_1047, msh1_1048, msi_1389, msi_1390, \
                         msi_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1781[k] = f_4 * msh0_1046[k]
                    - f_5 * msh1_1046[k]
                    + f_3 * pc_x[k] * msi_1389[k];

        t_1782[k] = f_4 * msh0_1047[k]
                    - f_5 * msh1_1047[k]
                    + f_3 * pc_x[k] * msi_1390[k];

        t_1783[k] = f_4 * msh0_1048[k]
                    - f_5 * msh1_1048[k]
                    + f_3 * pc_x[k] * msi_1391[k];
    }

#pragma omp simd aligned(t_1784, t_1785, t_1786, t_1787, t_1788, t_1789, pc_x, msh0_1049, \
                         msh1_1049, msi_1392, msi_1393, msi_1394, msi_1395, msi_1396, \
                         msi_1397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1784[k] = f_4 * msh0_1049[k]
                    - f_5 * msh1_1049[k]
                    + f_3 * pc_x[k] * msi_1392[k];

        t_1785[k] = f_3 * pc_x[k] * msi_1393[k];

        t_1786[k] = f_3 * pc_x[k] * msi_1394[k];

        t_1787[k] = f_3 * pc_x[k] * msi_1395[k];

        t_1788[k] = f_3 * pc_x[k] * msi_1396[k];

        t_1789[k] = f_3 * pc_x[k] * msi_1397[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, t_1793, pc_x, pc_y, pc_z, lsi_1113, lsi_1141, \
                         msh0_1044, msh1_1044, msi_1393, msi_1398, \
                         msi_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_3 * pc_x[k] * msi_1398[k];

        t_1791[k] = f_3 * pc_x[k] * msi_1399[k];

        t_1792[k] = f_17 * lsi_1141[k]
                    + f_1 * msh0_1044[k]
                    - f_2 * msh1_1044[k]
                    + f_3 * pc_y[k] * msi_1393[k];

        t_1793[k] = f_16 * lsi_1113[k]
                    + f_3 * pc_z[k] * msi_1393[k];
    }

#pragma omp simd aligned(t_1794, t_1795, t_1796, pc_y, lsi_1143, lsi_1144, lsi_1145, \
                         msh0_1046, msh0_1047, msh0_1048, msh1_1046, msh1_1047, msh1_1048, \
                         msi_1395, msi_1396, msi_1397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1794[k] = f_17 * lsi_1143[k]
                    + f_10 * msh0_1046[k]
                    - f_11 * msh1_1046[k]
                    + f_3 * pc_y[k] * msi_1395[k];

        t_1795[k] = f_17 * lsi_1144[k]
                    + f_8 * msh0_1047[k]
                    - f_9 * msh1_1047[k]
                    + f_3 * pc_y[k] * msi_1396[k];

        t_1796[k] = f_17 * lsi_1145[k]
                    + f_6 * msh0_1048[k]
                    - f_7 * msh1_1048[k]
                    + f_3 * pc_y[k] * msi_1397[k];
    }

#pragma omp simd aligned(t_1797, t_1798, t_1799, pc_y, pc_z, lsi_1119, lsi_1146, lsi_1147, \
                         msh0_1049, msh1_1049, msi_1398, msi_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1797[k] = f_17 * lsi_1146[k]
                    + f_4 * msh0_1049[k]
                    - f_5 * msh1_1049[k]
                    + f_3 * pc_y[k] * msi_1398[k];

        t_1798[k] = f_17 * lsi_1147[k]
                    + f_3 * pc_y[k] * msi_1399[k];

        t_1799[k] = f_16 * lsi_1119[k]
                    + f_1 * msh0_1049[k]
                    - f_2 * msh1_1049[k]
                    + f_3 * pc_z[k] * msi_1399[k];
    }

#pragma omp simd aligned(t_1800, t_1801, t_1802, pc_x, msh0_1050, msh0_1051, msh0_1052, \
                         msh1_1050, msh1_1051, msh1_1052, msi_1400, msi_1401, \
                         msi_1402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1800[k] = f_1 * msh0_1050[k]
                    - f_2 * msh1_1050[k]
                    + f_3 * pc_x[k] * msi_1400[k];

        t_1801[k] = f_19 * msh0_1051[k]
                    - f_20 * msh1_1051[k]
                    + f_3 * pc_x[k] * msi_1401[k];

        t_1802[k] = f_19 * msh0_1052[k]
                    - f_20 * msh1_1052[k]
                    + f_3 * pc_x[k] * msi_1402[k];
    }

#pragma omp simd aligned(t_1803, t_1804, t_1805, pc_x, msh0_1053, msh0_1054, msh0_1055, \
                         msh1_1053, msh1_1054, msh1_1055, msi_1403, msi_1404, \
                         msi_1405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1803[k] = f_10 * msh0_1053[k]
                    - f_11 * msh1_1053[k]
                    + f_3 * pc_x[k] * msi_1403[k];

        t_1804[k] = f_10 * msh0_1054[k]
                    - f_11 * msh1_1054[k]
                    + f_3 * pc_x[k] * msi_1404[k];

        t_1805[k] = f_10 * msh0_1055[k]
                    - f_11 * msh1_1055[k]
                    + f_3 * pc_x[k] * msi_1405[k];
    }

#pragma omp simd aligned(t_1806, t_1807, t_1808, pc_x, msh0_1056, msh0_1057, msh0_1058, \
                         msh1_1056, msh1_1057, msh1_1058, msi_1406, msi_1407, \
                         msi_1408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1806[k] = f_8 * msh0_1056[k]
                    - f_9 * msh1_1056[k]
                    + f_3 * pc_x[k] * msi_1406[k];

        t_1807[k] = f_8 * msh0_1057[k]
                    - f_9 * msh1_1057[k]
                    + f_3 * pc_x[k] * msi_1407[k];

        t_1808[k] = f_8 * msh0_1058[k]
                    - f_9 * msh1_1058[k]
                    + f_3 * pc_x[k] * msi_1408[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, pc_x, msh0_1059, msh0_1060, msh0_1061, \
                         msh1_1059, msh1_1060, msh1_1061, msi_1409, msi_1410, \
                         msi_1411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = f_8 * msh0_1059[k]
                    - f_9 * msh1_1059[k]
                    + f_3 * pc_x[k] * msi_1409[k];

        t_1810[k] = f_6 * msh0_1060[k]
                    - f_7 * msh1_1060[k]
                    + f_3 * pc_x[k] * msi_1410[k];

        t_1811[k] = f_6 * msh0_1061[k]
                    - f_7 * msh1_1061[k]
                    + f_3 * pc_x[k] * msi_1411[k];
    }

#pragma omp simd aligned(t_1812, t_1813, t_1814, pc_x, msh0_1062, msh0_1063, msh0_1064, \
                         msh1_1062, msh1_1063, msh1_1064, msi_1412, msi_1413, \
                         msi_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1812[k] = f_6 * msh0_1062[k]
                    - f_7 * msh1_1062[k]
                    + f_3 * pc_x[k] * msi_1412[k];

        t_1813[k] = f_6 * msh0_1063[k]
                    - f_7 * msh1_1063[k]
                    + f_3 * pc_x[k] * msi_1413[k];

        t_1814[k] = f_6 * msh0_1064[k]
                    - f_7 * msh1_1064[k]
                    + f_3 * pc_x[k] * msi_1414[k];
    }

#pragma omp simd aligned(t_1815, t_1816, t_1817, pc_x, msh0_1065, msh0_1066, msh0_1067, \
                         msh1_1065, msh1_1066, msh1_1067, msi_1415, msi_1416, \
                         msi_1417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1815[k] = f_4 * msh0_1065[k]
                    - f_5 * msh1_1065[k]
                    + f_3 * pc_x[k] * msi_1415[k];

        t_1816[k] = f_4 * msh0_1066[k]
                    - f_5 * msh1_1066[k]
                    + f_3 * pc_x[k] * msi_1416[k];

        t_1817[k] = f_4 * msh0_1067[k]
                    - f_5 * msh1_1067[k]
                    + f_3 * pc_x[k] * msi_1417[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, t_1821, pc_x, msh0_1068, msh0_1069, \
                         msh0_1070, msh1_1068, msh1_1069, msh1_1070, msi_1418, msi_1419, \
                         msi_1420, msi_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = f_4 * msh0_1068[k]
                    - f_5 * msh1_1068[k]
                    + f_3 * pc_x[k] * msi_1418[k];

        t_1819[k] = f_4 * msh0_1069[k]
                    - f_5 * msh1_1069[k]
                    + f_3 * pc_x[k] * msi_1419[k];

        t_1820[k] = f_4 * msh0_1070[k]
                    - f_5 * msh1_1070[k]
                    + f_3 * pc_x[k] * msi_1420[k];

        t_1821[k] = f_3 * pc_x[k] * msi_1421[k];
    }

#pragma omp simd aligned(t_1822, t_1823, t_1824, t_1825, t_1826, t_1827, pc_x, msi_1422, \
                         msi_1423, msi_1424, msi_1425, msi_1426, \
                         msi_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1822[k] = f_3 * pc_x[k] * msi_1422[k];

        t_1823[k] = f_3 * pc_x[k] * msi_1423[k];

        t_1824[k] = f_3 * pc_x[k] * msi_1424[k];

        t_1825[k] = f_3 * pc_x[k] * msi_1425[k];

        t_1826[k] = f_3 * pc_x[k] * msi_1426[k];

        t_1827[k] = f_3 * pc_x[k] * msi_1427[k];
    }

#pragma omp simd aligned(t_1828, t_1829, t_1830, pc_y, pc_z, lsi_1141, lsi_1169, lsi_1171, \
                         msh0_1065, msh0_1067, msh1_1065, msh1_1067, msi_1421, \
                         msi_1423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1828[k] = f_16 * lsi_1169[k]
                    + f_1 * msh0_1065[k]
                    - f_2 * msh1_1065[k]
                    + f_3 * pc_y[k] * msi_1421[k];

        t_1829[k] = f_17 * lsi_1141[k]
                    + f_3 * pc_z[k] * msi_1421[k];

        t_1830[k] = f_16 * lsi_1171[k]
                    + f_10 * msh0_1067[k]
                    - f_11 * msh1_1067[k]
                    + f_3 * pc_y[k] * msi_1423[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, pc_y, lsi_1172, lsi_1173, lsi_1174, \
                         msh0_1068, msh0_1069, msh0_1070, msh1_1068, msh1_1069, msh1_1070, \
                         msi_1424, msi_1425, msi_1426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_16 * lsi_1172[k]
                    + f_8 * msh0_1068[k]
                    - f_9 * msh1_1068[k]
                    + f_3 * pc_y[k] * msi_1424[k];

        t_1832[k] = f_16 * lsi_1173[k]
                    + f_6 * msh0_1069[k]
                    - f_7 * msh1_1069[k]
                    + f_3 * pc_y[k] * msi_1425[k];

        t_1833[k] = f_16 * lsi_1174[k]
                    + f_4 * msh0_1070[k]
                    - f_5 * msh1_1070[k]
                    + f_3 * pc_y[k] * msi_1426[k];
    }

#pragma omp simd aligned(t_1834, t_1835, t_1836, pc_x, pc_y, pc_z, lsi_1147, lsi_1175, \
                         msh0_1070, msh0_1071, msh1_1070, msh1_1071, msi_1427, \
                         msi_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = f_16 * lsi_1175[k]
                    + f_3 * pc_y[k] * msi_1427[k];

        t_1835[k] = f_17 * lsi_1147[k]
                    + f_1 * msh0_1070[k]
                    - f_2 * msh1_1070[k]
                    + f_3 * pc_z[k] * msi_1427[k];

        t_1836[k] = f_1 * msh0_1071[k]
                    - f_2 * msh1_1071[k]
                    + f_3 * pc_x[k] * msi_1428[k];
    }

#pragma omp simd aligned(t_1837, t_1838, t_1839, pc_x, msh0_1072, msh0_1073, msh0_1074, \
                         msh1_1072, msh1_1073, msh1_1074, msi_1429, msi_1430, \
                         msi_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1837[k] = f_19 * msh0_1072[k]
                    - f_20 * msh1_1072[k]
                    + f_3 * pc_x[k] * msi_1429[k];

        t_1838[k] = f_19 * msh0_1073[k]
                    - f_20 * msh1_1073[k]
                    + f_3 * pc_x[k] * msi_1430[k];

        t_1839[k] = f_10 * msh0_1074[k]
                    - f_11 * msh1_1074[k]
                    + f_3 * pc_x[k] * msi_1431[k];
    }

#pragma omp simd aligned(t_1840, t_1841, t_1842, pc_x, msh0_1075, msh0_1076, msh0_1077, \
                         msh1_1075, msh1_1076, msh1_1077, msi_1432, msi_1433, \
                         msi_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1840[k] = f_10 * msh0_1075[k]
                    - f_11 * msh1_1075[k]
                    + f_3 * pc_x[k] * msi_1432[k];

        t_1841[k] = f_10 * msh0_1076[k]
                    - f_11 * msh1_1076[k]
                    + f_3 * pc_x[k] * msi_1433[k];

        t_1842[k] = f_8 * msh0_1077[k]
                    - f_9 * msh1_1077[k]
                    + f_3 * pc_x[k] * msi_1434[k];
    }

#pragma omp simd aligned(t_1843, t_1844, t_1845, pc_x, msh0_1078, msh0_1079, msh0_1080, \
                         msh1_1078, msh1_1079, msh1_1080, msi_1435, msi_1436, \
                         msi_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1843[k] = f_8 * msh0_1078[k]
                    - f_9 * msh1_1078[k]
                    + f_3 * pc_x[k] * msi_1435[k];

        t_1844[k] = f_8 * msh0_1079[k]
                    - f_9 * msh1_1079[k]
                    + f_3 * pc_x[k] * msi_1436[k];

        t_1845[k] = f_8 * msh0_1080[k]
                    - f_9 * msh1_1080[k]
                    + f_3 * pc_x[k] * msi_1437[k];
    }

#pragma omp simd aligned(t_1846, t_1847, t_1848, pc_x, msh0_1081, msh0_1082, msh0_1083, \
                         msh1_1081, msh1_1082, msh1_1083, msi_1438, msi_1439, \
                         msi_1440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1846[k] = f_6 * msh0_1081[k]
                    - f_7 * msh1_1081[k]
                    + f_3 * pc_x[k] * msi_1438[k];

        t_1847[k] = f_6 * msh0_1082[k]
                    - f_7 * msh1_1082[k]
                    + f_3 * pc_x[k] * msi_1439[k];

        t_1848[k] = f_6 * msh0_1083[k]
                    - f_7 * msh1_1083[k]
                    + f_3 * pc_x[k] * msi_1440[k];
    }

#pragma omp simd aligned(t_1849, t_1850, t_1851, pc_x, msh0_1084, msh0_1085, msh0_1086, \
                         msh1_1084, msh1_1085, msh1_1086, msi_1441, msi_1442, \
                         msi_1443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1849[k] = f_6 * msh0_1084[k]
                    - f_7 * msh1_1084[k]
                    + f_3 * pc_x[k] * msi_1441[k];

        t_1850[k] = f_6 * msh0_1085[k]
                    - f_7 * msh1_1085[k]
                    + f_3 * pc_x[k] * msi_1442[k];

        t_1851[k] = f_4 * msh0_1086[k]
                    - f_5 * msh1_1086[k]
                    + f_3 * pc_x[k] * msi_1443[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsk0,
                                                           const size_t lsi, const size_t lsk1,
                                                           const size_t msh0, const size_t msh1,
                                                           const size_t msi, const size_t ncols,
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
    const auto f_18 = 4.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsk0_1584 = buffer.data(lsk0 + 1584);
    const auto *lsk0_1586 = buffer.data(lsk0 + 1586);
    const auto *lsk0_1589 = buffer.data(lsk0 + 1589);
    const auto *lsk0_1593 = buffer.data(lsk0 + 1593);
    const auto *lsk0_1598 = buffer.data(lsk0 + 1598);
    const auto *lsk0_1604 = buffer.data(lsk0 + 1604);
    const auto *lsk0_1612 = buffer.data(lsk0 + 1612);
    const auto *lsk0_1614 = buffer.data(lsk0 + 1614);
    const auto *lsk0_1615 = buffer.data(lsk0 + 1615);
    const auto *lsk0_1616 = buffer.data(lsk0 + 1616);
    const auto *lsk0_1617 = buffer.data(lsk0 + 1617);
    const auto *lsk0_1619 = buffer.data(lsk0 + 1619);

    const auto *lsi_1169 = buffer.data(lsi + 1169);
    const auto *lsi_1175 = buffer.data(lsi + 1175);
    const auto *lsi_1197 = buffer.data(lsi + 1197);
    const auto *lsi_1199 = buffer.data(lsi + 1199);
    const auto *lsi_1200 = buffer.data(lsi + 1200);
    const auto *lsi_1201 = buffer.data(lsi + 1201);
    const auto *lsi_1202 = buffer.data(lsi + 1202);
    const auto *lsi_1203 = buffer.data(lsi + 1203);
    const auto *lsi_1225 = buffer.data(lsi + 1225);
    const auto *lsi_1227 = buffer.data(lsi + 1227);
    const auto *lsi_1228 = buffer.data(lsi + 1228);
    const auto *lsi_1229 = buffer.data(lsi + 1229);
    const auto *lsi_1230 = buffer.data(lsi + 1230);
    const auto *lsi_1231 = buffer.data(lsi + 1231);
    const auto *lsi_1253 = buffer.data(lsi + 1253);
    const auto *lsi_1255 = buffer.data(lsi + 1255);
    const auto *lsi_1256 = buffer.data(lsi + 1256);
    const auto *lsi_1257 = buffer.data(lsi + 1257);
    const auto *lsi_1258 = buffer.data(lsi + 1258);
    const auto *lsi_1259 = buffer.data(lsi + 1259);

    const auto *lsk1_1584 = buffer.data(lsk1 + 1584);
    const auto *lsk1_1586 = buffer.data(lsk1 + 1586);
    const auto *lsk1_1589 = buffer.data(lsk1 + 1589);
    const auto *lsk1_1593 = buffer.data(lsk1 + 1593);
    const auto *lsk1_1598 = buffer.data(lsk1 + 1598);
    const auto *lsk1_1604 = buffer.data(lsk1 + 1604);
    const auto *lsk1_1612 = buffer.data(lsk1 + 1612);
    const auto *lsk1_1614 = buffer.data(lsk1 + 1614);
    const auto *lsk1_1615 = buffer.data(lsk1 + 1615);
    const auto *lsk1_1616 = buffer.data(lsk1 + 1616);
    const auto *lsk1_1617 = buffer.data(lsk1 + 1617);
    const auto *lsk1_1619 = buffer.data(lsk1 + 1619);

    const auto *msh0_1086 = buffer.data(msh0 + 1086);
    const auto *msh0_1087 = buffer.data(msh0 + 1087);
    const auto *msh0_1088 = buffer.data(msh0 + 1088);
    const auto *msh0_1089 = buffer.data(msh0 + 1089);
    const auto *msh0_1090 = buffer.data(msh0 + 1090);
    const auto *msh0_1091 = buffer.data(msh0 + 1091);
    const auto *msh0_1092 = buffer.data(msh0 + 1092);
    const auto *msh0_1093 = buffer.data(msh0 + 1093);
    const auto *msh0_1094 = buffer.data(msh0 + 1094);
    const auto *msh0_1095 = buffer.data(msh0 + 1095);
    const auto *msh0_1096 = buffer.data(msh0 + 1096);
    const auto *msh0_1097 = buffer.data(msh0 + 1097);
    const auto *msh0_1098 = buffer.data(msh0 + 1098);
    const auto *msh0_1099 = buffer.data(msh0 + 1099);
    const auto *msh0_1100 = buffer.data(msh0 + 1100);
    const auto *msh0_1101 = buffer.data(msh0 + 1101);
    const auto *msh0_1102 = buffer.data(msh0 + 1102);
    const auto *msh0_1103 = buffer.data(msh0 + 1103);
    const auto *msh0_1104 = buffer.data(msh0 + 1104);
    const auto *msh0_1105 = buffer.data(msh0 + 1105);
    const auto *msh0_1106 = buffer.data(msh0 + 1106);
    const auto *msh0_1107 = buffer.data(msh0 + 1107);
    const auto *msh0_1108 = buffer.data(msh0 + 1108);
    const auto *msh0_1109 = buffer.data(msh0 + 1109);
    const auto *msh0_1110 = buffer.data(msh0 + 1110);
    const auto *msh0_1111 = buffer.data(msh0 + 1111);
    const auto *msh0_1112 = buffer.data(msh0 + 1112);
    const auto *msh0_1114 = buffer.data(msh0 + 1114);
    const auto *msh0_1116 = buffer.data(msh0 + 1116);
    const auto *msh0_1117 = buffer.data(msh0 + 1117);
    const auto *msh0_1119 = buffer.data(msh0 + 1119);
    const auto *msh0_1120 = buffer.data(msh0 + 1120);
    const auto *msh0_1121 = buffer.data(msh0 + 1121);
    const auto *msh0_1123 = buffer.data(msh0 + 1123);
    const auto *msh0_1124 = buffer.data(msh0 + 1124);
    const auto *msh0_1125 = buffer.data(msh0 + 1125);
    const auto *msh0_1126 = buffer.data(msh0 + 1126);
    const auto *msh0_1128 = buffer.data(msh0 + 1128);
    const auto *msh0_1129 = buffer.data(msh0 + 1129);
    const auto *msh0_1130 = buffer.data(msh0 + 1130);
    const auto *msh0_1131 = buffer.data(msh0 + 1131);
    const auto *msh0_1132 = buffer.data(msh0 + 1132);
    const auto *msh0_1134 = buffer.data(msh0 + 1134);
    const auto *msh0_1136 = buffer.data(msh0 + 1136);
    const auto *msh0_1137 = buffer.data(msh0 + 1137);
    const auto *msh0_1139 = buffer.data(msh0 + 1139);
    const auto *msh0_1140 = buffer.data(msh0 + 1140);
    const auto *msh0_1141 = buffer.data(msh0 + 1141);
    const auto *msh0_1143 = buffer.data(msh0 + 1143);
    const auto *msh0_1144 = buffer.data(msh0 + 1144);
    const auto *msh0_1145 = buffer.data(msh0 + 1145);
    const auto *msh0_1146 = buffer.data(msh0 + 1146);
    const auto *msh0_1148 = buffer.data(msh0 + 1148);
    const auto *msh0_1149 = buffer.data(msh0 + 1149);
    const auto *msh0_1150 = buffer.data(msh0 + 1150);
    const auto *msh0_1151 = buffer.data(msh0 + 1151);
    const auto *msh0_1152 = buffer.data(msh0 + 1152);
    const auto *msh0_1154 = buffer.data(msh0 + 1154);

    const auto *msh1_1086 = buffer.data(msh1 + 1086);
    const auto *msh1_1087 = buffer.data(msh1 + 1087);
    const auto *msh1_1088 = buffer.data(msh1 + 1088);
    const auto *msh1_1089 = buffer.data(msh1 + 1089);
    const auto *msh1_1090 = buffer.data(msh1 + 1090);
    const auto *msh1_1091 = buffer.data(msh1 + 1091);
    const auto *msh1_1092 = buffer.data(msh1 + 1092);
    const auto *msh1_1093 = buffer.data(msh1 + 1093);
    const auto *msh1_1094 = buffer.data(msh1 + 1094);
    const auto *msh1_1095 = buffer.data(msh1 + 1095);
    const auto *msh1_1096 = buffer.data(msh1 + 1096);
    const auto *msh1_1097 = buffer.data(msh1 + 1097);
    const auto *msh1_1098 = buffer.data(msh1 + 1098);
    const auto *msh1_1099 = buffer.data(msh1 + 1099);
    const auto *msh1_1100 = buffer.data(msh1 + 1100);
    const auto *msh1_1101 = buffer.data(msh1 + 1101);
    const auto *msh1_1102 = buffer.data(msh1 + 1102);
    const auto *msh1_1103 = buffer.data(msh1 + 1103);
    const auto *msh1_1104 = buffer.data(msh1 + 1104);
    const auto *msh1_1105 = buffer.data(msh1 + 1105);
    const auto *msh1_1106 = buffer.data(msh1 + 1106);
    const auto *msh1_1107 = buffer.data(msh1 + 1107);
    const auto *msh1_1108 = buffer.data(msh1 + 1108);
    const auto *msh1_1109 = buffer.data(msh1 + 1109);
    const auto *msh1_1110 = buffer.data(msh1 + 1110);
    const auto *msh1_1111 = buffer.data(msh1 + 1111);
    const auto *msh1_1112 = buffer.data(msh1 + 1112);
    const auto *msh1_1114 = buffer.data(msh1 + 1114);
    const auto *msh1_1116 = buffer.data(msh1 + 1116);
    const auto *msh1_1117 = buffer.data(msh1 + 1117);
    const auto *msh1_1119 = buffer.data(msh1 + 1119);
    const auto *msh1_1120 = buffer.data(msh1 + 1120);
    const auto *msh1_1121 = buffer.data(msh1 + 1121);
    const auto *msh1_1123 = buffer.data(msh1 + 1123);
    const auto *msh1_1124 = buffer.data(msh1 + 1124);
    const auto *msh1_1125 = buffer.data(msh1 + 1125);
    const auto *msh1_1126 = buffer.data(msh1 + 1126);
    const auto *msh1_1128 = buffer.data(msh1 + 1128);
    const auto *msh1_1129 = buffer.data(msh1 + 1129);
    const auto *msh1_1130 = buffer.data(msh1 + 1130);
    const auto *msh1_1131 = buffer.data(msh1 + 1131);
    const auto *msh1_1132 = buffer.data(msh1 + 1132);
    const auto *msh1_1134 = buffer.data(msh1 + 1134);
    const auto *msh1_1136 = buffer.data(msh1 + 1136);
    const auto *msh1_1137 = buffer.data(msh1 + 1137);
    const auto *msh1_1139 = buffer.data(msh1 + 1139);
    const auto *msh1_1140 = buffer.data(msh1 + 1140);
    const auto *msh1_1141 = buffer.data(msh1 + 1141);
    const auto *msh1_1143 = buffer.data(msh1 + 1143);
    const auto *msh1_1144 = buffer.data(msh1 + 1144);
    const auto *msh1_1145 = buffer.data(msh1 + 1145);
    const auto *msh1_1146 = buffer.data(msh1 + 1146);
    const auto *msh1_1148 = buffer.data(msh1 + 1148);
    const auto *msh1_1149 = buffer.data(msh1 + 1149);
    const auto *msh1_1150 = buffer.data(msh1 + 1150);
    const auto *msh1_1151 = buffer.data(msh1 + 1151);
    const auto *msh1_1152 = buffer.data(msh1 + 1152);
    const auto *msh1_1154 = buffer.data(msh1 + 1154);

    const auto *msi_1444 = buffer.data(msi + 1444);
    const auto *msi_1445 = buffer.data(msi + 1445);
    const auto *msi_1446 = buffer.data(msi + 1446);
    const auto *msi_1447 = buffer.data(msi + 1447);
    const auto *msi_1448 = buffer.data(msi + 1448);
    const auto *msi_1449 = buffer.data(msi + 1449);
    const auto *msi_1450 = buffer.data(msi + 1450);
    const auto *msi_1451 = buffer.data(msi + 1451);
    const auto *msi_1452 = buffer.data(msi + 1452);
    const auto *msi_1453 = buffer.data(msi + 1453);
    const auto *msi_1454 = buffer.data(msi + 1454);
    const auto *msi_1455 = buffer.data(msi + 1455);
    const auto *msi_1456 = buffer.data(msi + 1456);
    const auto *msi_1457 = buffer.data(msi + 1457);
    const auto *msi_1458 = buffer.data(msi + 1458);
    const auto *msi_1459 = buffer.data(msi + 1459);
    const auto *msi_1460 = buffer.data(msi + 1460);
    const auto *msi_1461 = buffer.data(msi + 1461);
    const auto *msi_1462 = buffer.data(msi + 1462);
    const auto *msi_1463 = buffer.data(msi + 1463);
    const auto *msi_1464 = buffer.data(msi + 1464);
    const auto *msi_1465 = buffer.data(msi + 1465);
    const auto *msi_1466 = buffer.data(msi + 1466);
    const auto *msi_1467 = buffer.data(msi + 1467);
    const auto *msi_1468 = buffer.data(msi + 1468);
    const auto *msi_1469 = buffer.data(msi + 1469);
    const auto *msi_1470 = buffer.data(msi + 1470);
    const auto *msi_1471 = buffer.data(msi + 1471);
    const auto *msi_1472 = buffer.data(msi + 1472);
    const auto *msi_1473 = buffer.data(msi + 1473);
    const auto *msi_1474 = buffer.data(msi + 1474);
    const auto *msi_1475 = buffer.data(msi + 1475);
    const auto *msi_1476 = buffer.data(msi + 1476);
    const auto *msi_1477 = buffer.data(msi + 1477);
    const auto *msi_1478 = buffer.data(msi + 1478);
    const auto *msi_1479 = buffer.data(msi + 1479);
    const auto *msi_1480 = buffer.data(msi + 1480);
    const auto *msi_1481 = buffer.data(msi + 1481);
    const auto *msi_1482 = buffer.data(msi + 1482);
    const auto *msi_1483 = buffer.data(msi + 1483);
    const auto *msi_1485 = buffer.data(msi + 1485);
    const auto *msi_1487 = buffer.data(msi + 1487);
    const auto *msi_1488 = buffer.data(msi + 1488);
    const auto *msi_1490 = buffer.data(msi + 1490);
    const auto *msi_1491 = buffer.data(msi + 1491);
    const auto *msi_1492 = buffer.data(msi + 1492);
    const auto *msi_1494 = buffer.data(msi + 1494);
    const auto *msi_1495 = buffer.data(msi + 1495);
    const auto *msi_1496 = buffer.data(msi + 1496);
    const auto *msi_1497 = buffer.data(msi + 1497);
    const auto *msi_1499 = buffer.data(msi + 1499);
    const auto *msi_1500 = buffer.data(msi + 1500);
    const auto *msi_1501 = buffer.data(msi + 1501);
    const auto *msi_1502 = buffer.data(msi + 1502);
    const auto *msi_1503 = buffer.data(msi + 1503);
    const auto *msi_1505 = buffer.data(msi + 1505);
    const auto *msi_1506 = buffer.data(msi + 1506);
    const auto *msi_1507 = buffer.data(msi + 1507);
    const auto *msi_1508 = buffer.data(msi + 1508);
    const auto *msi_1509 = buffer.data(msi + 1509);
    const auto *msi_1510 = buffer.data(msi + 1510);
    const auto *msi_1511 = buffer.data(msi + 1511);
    const auto *msi_1512 = buffer.data(msi + 1512);
    const auto *msi_1514 = buffer.data(msi + 1514);
    const auto *msi_1515 = buffer.data(msi + 1515);
    const auto *msi_1517 = buffer.data(msi + 1517);
    const auto *msi_1518 = buffer.data(msi + 1518);
    const auto *msi_1519 = buffer.data(msi + 1519);
    const auto *msi_1521 = buffer.data(msi + 1521);
    const auto *msi_1522 = buffer.data(msi + 1522);
    const auto *msi_1523 = buffer.data(msi + 1523);
    const auto *msi_1524 = buffer.data(msi + 1524);
    const auto *msi_1526 = buffer.data(msi + 1526);
    const auto *msi_1527 = buffer.data(msi + 1527);
    const auto *msi_1528 = buffer.data(msi + 1528);
    const auto *msi_1529 = buffer.data(msi + 1529);
    const auto *msi_1530 = buffer.data(msi + 1530);
    const auto *msi_1532 = buffer.data(msi + 1532);
    const auto *msi_1533 = buffer.data(msi + 1533);
    const auto *msi_1534 = buffer.data(msi + 1534);
    const auto *msi_1535 = buffer.data(msi + 1535);
    const auto *msi_1536 = buffer.data(msi + 1536);
    const auto *msi_1537 = buffer.data(msi + 1537);
    const auto *msi_1538 = buffer.data(msi + 1538);
    const auto *msi_1539 = buffer.data(msi + 1539);

#pragma omp simd aligned(t_1852, t_1853, t_1854, pc_x, msh0_1087, msh0_1088, msh0_1089, \
                         msh1_1087, msh1_1088, msh1_1089, msi_1444, msi_1445, \
                         msi_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1852[k] = f_4 * msh0_1087[k]
                    - f_5 * msh1_1087[k]
                    + f_3 * pc_x[k] * msi_1444[k];

        t_1853[k] = f_4 * msh0_1088[k]
                    - f_5 * msh1_1088[k]
                    + f_3 * pc_x[k] * msi_1445[k];

        t_1854[k] = f_4 * msh0_1089[k]
                    - f_5 * msh1_1089[k]
                    + f_3 * pc_x[k] * msi_1446[k];
    }

#pragma omp simd aligned(t_1855, t_1856, t_1857, t_1858, t_1859, pc_x, msh0_1090, msh0_1091, \
                         msh1_1090, msh1_1091, msi_1447, msi_1448, msi_1449, msi_1450, \
                         msi_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1855[k] = f_4 * msh0_1090[k]
                    - f_5 * msh1_1090[k]
                    + f_3 * pc_x[k] * msi_1447[k];

        t_1856[k] = f_4 * msh0_1091[k]
                    - f_5 * msh1_1091[k]
                    + f_3 * pc_x[k] * msi_1448[k];

        t_1857[k] = f_3 * pc_x[k] * msi_1449[k];

        t_1858[k] = f_3 * pc_x[k] * msi_1450[k];

        t_1859[k] = f_3 * pc_x[k] * msi_1451[k];
    }

#pragma omp simd aligned(t_1860, t_1861, t_1862, t_1863, t_1864, pc_x, pc_y, lsi_1197, \
                         msh0_1086, msh1_1086, msi_1449, msi_1452, msi_1453, msi_1454, \
                         msi_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1860[k] = f_3 * pc_x[k] * msi_1452[k];

        t_1861[k] = f_3 * pc_x[k] * msi_1453[k];

        t_1862[k] = f_3 * pc_x[k] * msi_1454[k];

        t_1863[k] = f_3 * pc_x[k] * msi_1455[k];

        t_1864[k] = f_15 * lsi_1197[k]
                    + f_1 * msh0_1086[k]
                    - f_2 * msh1_1086[k]
                    + f_3 * pc_y[k] * msi_1449[k];
    }

#pragma omp simd aligned(t_1865, t_1866, t_1867, pc_y, pc_z, lsi_1169, lsi_1199, lsi_1200, \
                         msh0_1088, msh0_1089, msh1_1088, msh1_1089, msi_1449, msi_1451, \
                         msi_1452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1865[k] = f_22 * lsi_1169[k]
                    + f_3 * pc_z[k] * msi_1449[k];

        t_1866[k] = f_15 * lsi_1199[k]
                    + f_10 * msh0_1088[k]
                    - f_11 * msh1_1088[k]
                    + f_3 * pc_y[k] * msi_1451[k];

        t_1867[k] = f_15 * lsi_1200[k]
                    + f_8 * msh0_1089[k]
                    - f_9 * msh1_1089[k]
                    + f_3 * pc_y[k] * msi_1452[k];
    }

#pragma omp simd aligned(t_1868, t_1869, t_1870, pc_y, lsi_1201, lsi_1202, lsi_1203, \
                         msh0_1090, msh0_1091, msh1_1090, msh1_1091, msi_1453, msi_1454, \
                         msi_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1868[k] = f_15 * lsi_1201[k]
                    + f_6 * msh0_1090[k]
                    - f_7 * msh1_1090[k]
                    + f_3 * pc_y[k] * msi_1453[k];

        t_1869[k] = f_15 * lsi_1202[k]
                    + f_4 * msh0_1091[k]
                    - f_5 * msh1_1091[k]
                    + f_3 * pc_y[k] * msi_1454[k];

        t_1870[k] = f_15 * lsi_1203[k]
                    + f_3 * pc_y[k] * msi_1455[k];
    }

#pragma omp simd aligned(t_1871, t_1872, t_1873, pc_x, pc_z, lsi_1175, msh0_1091, msh0_1092, \
                         msh0_1093, msh1_1091, msh1_1092, msh1_1093, msi_1455, msi_1456, \
                         msi_1457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1871[k] = f_22 * lsi_1175[k]
                    + f_1 * msh0_1091[k]
                    - f_2 * msh1_1091[k]
                    + f_3 * pc_z[k] * msi_1455[k];

        t_1872[k] = f_1 * msh0_1092[k]
                    - f_2 * msh1_1092[k]
                    + f_3 * pc_x[k] * msi_1456[k];

        t_1873[k] = f_19 * msh0_1093[k]
                    - f_20 * msh1_1093[k]
                    + f_3 * pc_x[k] * msi_1457[k];
    }

#pragma omp simd aligned(t_1874, t_1875, t_1876, pc_x, msh0_1094, msh0_1095, msh0_1096, \
                         msh1_1094, msh1_1095, msh1_1096, msi_1458, msi_1459, \
                         msi_1460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1874[k] = f_19 * msh0_1094[k]
                    - f_20 * msh1_1094[k]
                    + f_3 * pc_x[k] * msi_1458[k];

        t_1875[k] = f_10 * msh0_1095[k]
                    - f_11 * msh1_1095[k]
                    + f_3 * pc_x[k] * msi_1459[k];

        t_1876[k] = f_10 * msh0_1096[k]
                    - f_11 * msh1_1096[k]
                    + f_3 * pc_x[k] * msi_1460[k];
    }

#pragma omp simd aligned(t_1877, t_1878, t_1879, pc_x, msh0_1097, msh0_1098, msh0_1099, \
                         msh1_1097, msh1_1098, msh1_1099, msi_1461, msi_1462, \
                         msi_1463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1877[k] = f_10 * msh0_1097[k]
                    - f_11 * msh1_1097[k]
                    + f_3 * pc_x[k] * msi_1461[k];

        t_1878[k] = f_8 * msh0_1098[k]
                    - f_9 * msh1_1098[k]
                    + f_3 * pc_x[k] * msi_1462[k];

        t_1879[k] = f_8 * msh0_1099[k]
                    - f_9 * msh1_1099[k]
                    + f_3 * pc_x[k] * msi_1463[k];
    }

#pragma omp simd aligned(t_1880, t_1881, t_1882, pc_x, msh0_1100, msh0_1101, msh0_1102, \
                         msh1_1100, msh1_1101, msh1_1102, msi_1464, msi_1465, \
                         msi_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1880[k] = f_8 * msh0_1100[k]
                    - f_9 * msh1_1100[k]
                    + f_3 * pc_x[k] * msi_1464[k];

        t_1881[k] = f_8 * msh0_1101[k]
                    - f_9 * msh1_1101[k]
                    + f_3 * pc_x[k] * msi_1465[k];

        t_1882[k] = f_6 * msh0_1102[k]
                    - f_7 * msh1_1102[k]
                    + f_3 * pc_x[k] * msi_1466[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pc_x, msh0_1103, msh0_1104, msh0_1105, \
                         msh1_1103, msh1_1104, msh1_1105, msi_1467, msi_1468, \
                         msi_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_6 * msh0_1103[k]
                    - f_7 * msh1_1103[k]
                    + f_3 * pc_x[k] * msi_1467[k];

        t_1884[k] = f_6 * msh0_1104[k]
                    - f_7 * msh1_1104[k]
                    + f_3 * pc_x[k] * msi_1468[k];

        t_1885[k] = f_6 * msh0_1105[k]
                    - f_7 * msh1_1105[k]
                    + f_3 * pc_x[k] * msi_1469[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pc_x, msh0_1106, msh0_1107, msh0_1108, \
                         msh1_1106, msh1_1107, msh1_1108, msi_1470, msi_1471, \
                         msi_1472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = f_6 * msh0_1106[k]
                    - f_7 * msh1_1106[k]
                    + f_3 * pc_x[k] * msi_1470[k];

        t_1887[k] = f_4 * msh0_1107[k]
                    - f_5 * msh1_1107[k]
                    + f_3 * pc_x[k] * msi_1471[k];

        t_1888[k] = f_4 * msh0_1108[k]
                    - f_5 * msh1_1108[k]
                    + f_3 * pc_x[k] * msi_1472[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, pc_x, msh0_1109, msh0_1110, msh0_1111, \
                         msh1_1109, msh1_1110, msh1_1111, msi_1473, msi_1474, \
                         msi_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = f_4 * msh0_1109[k]
                    - f_5 * msh1_1109[k]
                    + f_3 * pc_x[k] * msi_1473[k];

        t_1890[k] = f_4 * msh0_1110[k]
                    - f_5 * msh1_1110[k]
                    + f_3 * pc_x[k] * msi_1474[k];

        t_1891[k] = f_4 * msh0_1111[k]
                    - f_5 * msh1_1111[k]
                    + f_3 * pc_x[k] * msi_1475[k];
    }

#pragma omp simd aligned(t_1892, t_1893, t_1894, t_1895, t_1896, t_1897, pc_x, msh0_1112, \
                         msh1_1112, msi_1476, msi_1477, msi_1478, msi_1479, msi_1480, \
                         msi_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1892[k] = f_4 * msh0_1112[k]
                    - f_5 * msh1_1112[k]
                    + f_3 * pc_x[k] * msi_1476[k];

        t_1893[k] = f_3 * pc_x[k] * msi_1477[k];

        t_1894[k] = f_3 * pc_x[k] * msi_1478[k];

        t_1895[k] = f_3 * pc_x[k] * msi_1479[k];

        t_1896[k] = f_3 * pc_x[k] * msi_1480[k];

        t_1897[k] = f_3 * pc_x[k] * msi_1481[k];
    }

#pragma omp simd aligned(t_1898, t_1899, t_1900, t_1901, pc_x, pc_y, pc_z, lsi_1197, lsi_1225, \
                         msh0_1107, msh1_1107, msi_1477, msi_1482, \
                         msi_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1898[k] = f_3 * pc_x[k] * msi_1482[k];

        t_1899[k] = f_3 * pc_x[k] * msi_1483[k];

        t_1900[k] = f_14 * lsi_1225[k]
                    + f_1 * msh0_1107[k]
                    - f_2 * msh1_1107[k]
                    + f_3 * pc_y[k] * msi_1477[k];

        t_1901[k] = f_21 * lsi_1197[k]
                    + f_3 * pc_z[k] * msi_1477[k];
    }

#pragma omp simd aligned(t_1902, t_1903, t_1904, pc_y, lsi_1227, lsi_1228, lsi_1229, \
                         msh0_1109, msh0_1110, msh0_1111, msh1_1109, msh1_1110, msh1_1111, \
                         msi_1479, msi_1480, msi_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1902[k] = f_14 * lsi_1227[k]
                    + f_10 * msh0_1109[k]
                    - f_11 * msh1_1109[k]
                    + f_3 * pc_y[k] * msi_1479[k];

        t_1903[k] = f_14 * lsi_1228[k]
                    + f_8 * msh0_1110[k]
                    - f_9 * msh1_1110[k]
                    + f_3 * pc_y[k] * msi_1480[k];

        t_1904[k] = f_14 * lsi_1229[k]
                    + f_6 * msh0_1111[k]
                    - f_7 * msh1_1111[k]
                    + f_3 * pc_y[k] * msi_1481[k];
    }

#pragma omp simd aligned(t_1905, t_1906, t_1907, t_1908, pa_y, pc_y, pc_z, lsk0_1584, \
                         lsi_1203, lsi_1230, lsi_1231, lsk1_1584, msh0_1112, msh1_1112, \
                         msi_1482, msi_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1905[k] = f_14 * lsi_1230[k]
                    + f_4 * msh0_1112[k]
                    - f_5 * msh1_1112[k]
                    + f_3 * pc_y[k] * msi_1482[k];

        t_1906[k] = f_14 * lsi_1231[k]
                    + f_3 * pc_y[k] * msi_1483[k];

        t_1907[k] = f_21 * lsi_1203[k]
                    + f_1 * msh0_1112[k]
                    - f_2 * msh1_1112[k]
                    + f_3 * pc_z[k] * msi_1483[k];

        t_1908[k] = pa_y[k] * lsk0_1584[k]
                    - f_12 * pc_y[k] * lsk1_1584[k];
    }

#pragma omp simd aligned(t_1909, t_1910, t_1911, pa_y, pc_x, pc_y, lsk0_1586, lsk1_1586, \
                         msh0_1114, msh0_1116, msh1_1114, msh1_1116, msi_1485, \
                         msi_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1909[k] = f_19 * msh0_1114[k]
                    - f_20 * msh1_1114[k]
                    + f_3 * pc_x[k] * msi_1485[k];

        t_1910[k] = pa_y[k] * lsk0_1586[k]
                    - f_12 * pc_y[k] * lsk1_1586[k];

        t_1911[k] = f_10 * msh0_1116[k]
                    - f_11 * msh1_1116[k]
                    + f_3 * pc_x[k] * msi_1487[k];
    }

#pragma omp simd aligned(t_1912, t_1913, t_1914, pa_y, pc_x, pc_y, lsk0_1589, lsk1_1589, \
                         msh0_1117, msh0_1119, msh1_1117, msh1_1119, msi_1488, \
                         msi_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1912[k] = f_10 * msh0_1117[k]
                    - f_11 * msh1_1117[k]
                    + f_3 * pc_x[k] * msi_1488[k];

        t_1913[k] = pa_y[k] * lsk0_1589[k]
                    - f_12 * pc_y[k] * lsk1_1589[k];

        t_1914[k] = f_8 * msh0_1119[k]
                    - f_9 * msh1_1119[k]
                    + f_3 * pc_x[k] * msi_1490[k];
    }

#pragma omp simd aligned(t_1915, t_1916, t_1917, pa_y, pc_x, pc_y, lsk0_1593, lsk1_1593, \
                         msh0_1120, msh0_1121, msh1_1120, msh1_1121, msi_1491, \
                         msi_1492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1915[k] = f_8 * msh0_1120[k]
                    - f_9 * msh1_1120[k]
                    + f_3 * pc_x[k] * msi_1491[k];

        t_1916[k] = f_8 * msh0_1121[k]
                    - f_9 * msh1_1121[k]
                    + f_3 * pc_x[k] * msi_1492[k];

        t_1917[k] = pa_y[k] * lsk0_1593[k]
                    - f_12 * pc_y[k] * lsk1_1593[k];
    }

#pragma omp simd aligned(t_1918, t_1919, t_1920, pc_x, msh0_1123, msh0_1124, msh0_1125, \
                         msh1_1123, msh1_1124, msh1_1125, msi_1494, msi_1495, \
                         msi_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1918[k] = f_6 * msh0_1123[k]
                    - f_7 * msh1_1123[k]
                    + f_3 * pc_x[k] * msi_1494[k];

        t_1919[k] = f_6 * msh0_1124[k]
                    - f_7 * msh1_1124[k]
                    + f_3 * pc_x[k] * msi_1495[k];

        t_1920[k] = f_6 * msh0_1125[k]
                    - f_7 * msh1_1125[k]
                    + f_3 * pc_x[k] * msi_1496[k];
    }

#pragma omp simd aligned(t_1921, t_1922, t_1923, pa_y, pc_x, pc_y, lsk0_1598, lsk1_1598, \
                         msh0_1126, msh0_1128, msh1_1126, msh1_1128, msi_1497, \
                         msi_1499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1921[k] = f_6 * msh0_1126[k]
                    - f_7 * msh1_1126[k]
                    + f_3 * pc_x[k] * msi_1497[k];

        t_1922[k] = pa_y[k] * lsk0_1598[k]
                    - f_12 * pc_y[k] * lsk1_1598[k];

        t_1923[k] = f_4 * msh0_1128[k]
                    - f_5 * msh1_1128[k]
                    + f_3 * pc_x[k] * msi_1499[k];
    }

#pragma omp simd aligned(t_1924, t_1925, t_1926, pc_x, msh0_1129, msh0_1130, msh0_1131, \
                         msh1_1129, msh1_1130, msh1_1131, msi_1500, msi_1501, \
                         msi_1502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1924[k] = f_4 * msh0_1129[k]
                    - f_5 * msh1_1129[k]
                    + f_3 * pc_x[k] * msi_1500[k];

        t_1925[k] = f_4 * msh0_1130[k]
                    - f_5 * msh1_1130[k]
                    + f_3 * pc_x[k] * msi_1501[k];

        t_1926[k] = f_4 * msh0_1131[k]
                    - f_5 * msh1_1131[k]
                    + f_3 * pc_x[k] * msi_1502[k];
    }

#pragma omp simd aligned(t_1927, t_1928, t_1929, t_1930, t_1931, pa_y, pc_x, pc_y, lsk0_1604, \
                         lsk1_1604, msh0_1132, msh1_1132, msi_1503, msi_1505, msi_1506, \
                         msi_1507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1927[k] = f_4 * msh0_1132[k]
                    - f_5 * msh1_1132[k]
                    + f_3 * pc_x[k] * msi_1503[k];

        t_1928[k] = pa_y[k] * lsk0_1604[k]
                    - f_12 * pc_y[k] * lsk1_1604[k];

        t_1929[k] = f_3 * pc_x[k] * msi_1505[k];

        t_1930[k] = f_3 * pc_x[k] * msi_1506[k];

        t_1931[k] = f_3 * pc_x[k] * msi_1507[k];
    }

#pragma omp simd aligned(t_1932, t_1933, t_1934, t_1935, t_1936, pa_y, pc_x, pc_y, lsk0_1612, \
                         lsi_1253, lsk1_1612, msi_1508, msi_1509, msi_1510, \
                         msi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1932[k] = f_3 * pc_x[k] * msi_1508[k];

        t_1933[k] = f_3 * pc_x[k] * msi_1509[k];

        t_1934[k] = f_3 * pc_x[k] * msi_1510[k];

        t_1935[k] = f_3 * pc_x[k] * msi_1511[k];

        t_1936[k] = pa_y[k] * lsk0_1612[k]
                    + f_21 * lsi_1253[k]
                    - f_12 * pc_y[k] * lsk1_1612[k];
    }

#pragma omp simd aligned(t_1937, t_1938, t_1939, pa_y, pc_y, pc_z, lsk0_1614, lsk0_1615, \
                         lsi_1225, lsi_1255, lsi_1256, lsk1_1614, lsk1_1615, \
                         msi_1505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1937[k] = f_18 * lsi_1225[k]
                    + f_3 * pc_z[k] * msi_1505[k];

        t_1938[k] = pa_y[k] * lsk0_1614[k]
                    + f_17 * lsi_1255[k]
                    - f_12 * pc_y[k] * lsk1_1614[k];

        t_1939[k] = pa_y[k] * lsk0_1615[k]
                    + f_16 * lsi_1256[k]
                    - f_12 * pc_y[k] * lsk1_1615[k];
    }

#pragma omp simd aligned(t_1940, t_1941, t_1942, t_1943, pa_y, pc_y, lsk0_1616, lsk0_1617, \
                         lsk0_1619, lsi_1257, lsi_1258, lsi_1259, lsk1_1616, lsk1_1617, \
                         lsk1_1619, msi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1940[k] = pa_y[k] * lsk0_1616[k]
                    + f_15 * lsi_1257[k]
                    - f_12 * pc_y[k] * lsk1_1616[k];

        t_1941[k] = pa_y[k] * lsk0_1617[k]
                    + f_14 * lsi_1258[k]
                    - f_12 * pc_y[k] * lsk1_1617[k];

        t_1942[k] = f_13 * lsi_1259[k]
                    + f_3 * pc_y[k] * msi_1511[k];

        t_1943[k] = pa_y[k] * lsk0_1619[k]
                    - f_12 * pc_y[k] * lsk1_1619[k];
    }

#pragma omp simd aligned(t_1944, t_1945, t_1946, t_1947, t_1948, pc_x, pc_y, msh0_1134, \
                         msh0_1136, msh0_1137, msh1_1134, msh1_1136, msh1_1137, msi_1512, \
                         msi_1514, msi_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1944[k] = f_1 * msh0_1134[k]
                    - f_2 * msh1_1134[k]
                    + f_3 * pc_x[k] * msi_1512[k];

        t_1945[k] = f_3 * pc_y[k] * msi_1512[k];

        t_1946[k] = f_19 * msh0_1136[k]
                    - f_20 * msh1_1136[k]
                    + f_3 * pc_x[k] * msi_1514[k];

        t_1947[k] = f_10 * msh0_1137[k]
                    - f_11 * msh1_1137[k]
                    + f_3 * pc_x[k] * msi_1515[k];

        t_1948[k] = f_3 * pc_y[k] * msi_1514[k];
    }

#pragma omp simd aligned(t_1949, t_1950, t_1951, t_1952, pc_x, pc_y, msh0_1139, msh0_1140, \
                         msh0_1141, msh1_1139, msh1_1140, msh1_1141, msi_1517, msi_1518, \
                         msi_1519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1949[k] = f_10 * msh0_1139[k]
                    - f_11 * msh1_1139[k]
                    + f_3 * pc_x[k] * msi_1517[k];

        t_1950[k] = f_8 * msh0_1140[k]
                    - f_9 * msh1_1140[k]
                    + f_3 * pc_x[k] * msi_1518[k];

        t_1951[k] = f_8 * msh0_1141[k]
                    - f_9 * msh1_1141[k]
                    + f_3 * pc_x[k] * msi_1519[k];

        t_1952[k] = f_3 * pc_y[k] * msi_1517[k];
    }

#pragma omp simd aligned(t_1953, t_1954, t_1955, pc_x, msh0_1143, msh0_1144, msh0_1145, \
                         msh1_1143, msh1_1144, msh1_1145, msi_1521, msi_1522, \
                         msi_1523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1953[k] = f_8 * msh0_1143[k]
                    - f_9 * msh1_1143[k]
                    + f_3 * pc_x[k] * msi_1521[k];

        t_1954[k] = f_6 * msh0_1144[k]
                    - f_7 * msh1_1144[k]
                    + f_3 * pc_x[k] * msi_1522[k];

        t_1955[k] = f_6 * msh0_1145[k]
                    - f_7 * msh1_1145[k]
                    + f_3 * pc_x[k] * msi_1523[k];
    }

#pragma omp simd aligned(t_1956, t_1957, t_1958, t_1959, pc_x, pc_y, msh0_1146, msh0_1148, \
                         msh0_1149, msh1_1146, msh1_1148, msh1_1149, msi_1521, msi_1524, \
                         msi_1526, msi_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1956[k] = f_6 * msh0_1146[k]
                    - f_7 * msh1_1146[k]
                    + f_3 * pc_x[k] * msi_1524[k];

        t_1957[k] = f_3 * pc_y[k] * msi_1521[k];

        t_1958[k] = f_6 * msh0_1148[k]
                    - f_7 * msh1_1148[k]
                    + f_3 * pc_x[k] * msi_1526[k];

        t_1959[k] = f_4 * msh0_1149[k]
                    - f_5 * msh1_1149[k]
                    + f_3 * pc_x[k] * msi_1527[k];
    }

#pragma omp simd aligned(t_1960, t_1961, t_1962, t_1963, pc_x, pc_y, msh0_1150, msh0_1151, \
                         msh0_1152, msh1_1150, msh1_1151, msh1_1152, msi_1526, msi_1528, \
                         msi_1529, msi_1530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1960[k] = f_4 * msh0_1150[k]
                    - f_5 * msh1_1150[k]
                    + f_3 * pc_x[k] * msi_1528[k];

        t_1961[k] = f_4 * msh0_1151[k]
                    - f_5 * msh1_1151[k]
                    + f_3 * pc_x[k] * msi_1529[k];

        t_1962[k] = f_4 * msh0_1152[k]
                    - f_5 * msh1_1152[k]
                    + f_3 * pc_x[k] * msi_1530[k];

        t_1963[k] = f_3 * pc_y[k] * msi_1526[k];
    }

#pragma omp simd aligned(t_1964, t_1965, t_1966, t_1967, t_1968, t_1969, pc_x, msh0_1154, \
                         msh1_1154, msi_1532, msi_1533, msi_1534, msi_1535, msi_1536, \
                         msi_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1964[k] = f_4 * msh0_1154[k]
                    - f_5 * msh1_1154[k]
                    + f_3 * pc_x[k] * msi_1532[k];

        t_1965[k] = f_3 * pc_x[k] * msi_1533[k];

        t_1966[k] = f_3 * pc_x[k] * msi_1534[k];

        t_1967[k] = f_3 * pc_x[k] * msi_1535[k];

        t_1968[k] = f_3 * pc_x[k] * msi_1536[k];

        t_1969[k] = f_3 * pc_x[k] * msi_1537[k];
    }

#pragma omp simd aligned(t_1970, t_1971, t_1972, t_1973, pc_x, pc_y, msh0_1149, msh0_1150, \
                         msh1_1149, msh1_1150, msi_1533, msi_1534, msi_1538, \
                         msi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1970[k] = f_3 * pc_x[k] * msi_1538[k];

        t_1971[k] = f_3 * pc_x[k] * msi_1539[k];

        t_1972[k] = f_1 * msh0_1149[k]
                    - f_2 * msh1_1149[k]
                    + f_3 * pc_y[k] * msi_1533[k];

        t_1973[k] = f_19 * msh0_1150[k]
                    - f_20 * msh1_1150[k]
                    + f_3 * pc_y[k] * msi_1534[k];
    }
}

static auto
compute_prim_msk_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t lsi, const size_t msh0,
                                                           const size_t msh1, const size_t msi,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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

    auto *t_1974 = buffer.data(target + 1974);
    auto *t_1975 = buffer.data(target + 1975);
    auto *t_1976 = buffer.data(target + 1976);
    auto *t_1977 = buffer.data(target + 1977);
    auto *t_1978 = buffer.data(target + 1978);
    auto *t_1979 = buffer.data(target + 1979);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi_1259 = buffer.data(lsi + 1259);

    const auto *msh0_1151 = buffer.data(msh0 + 1151);
    const auto *msh0_1152 = buffer.data(msh0 + 1152);
    const auto *msh0_1153 = buffer.data(msh0 + 1153);
    const auto *msh0_1154 = buffer.data(msh0 + 1154);

    const auto *msh1_1151 = buffer.data(msh1 + 1151);
    const auto *msh1_1152 = buffer.data(msh1 + 1152);
    const auto *msh1_1153 = buffer.data(msh1 + 1153);
    const auto *msh1_1154 = buffer.data(msh1 + 1154);

    const auto *msi_1535 = buffer.data(msi + 1535);
    const auto *msi_1536 = buffer.data(msi + 1536);
    const auto *msi_1537 = buffer.data(msi + 1537);
    const auto *msi_1538 = buffer.data(msi + 1538);
    const auto *msi_1539 = buffer.data(msi + 1539);

#pragma omp simd aligned(t_1974, t_1975, t_1976, pc_y, msh0_1151, msh0_1152, msh0_1153, \
                         msh1_1151, msh1_1152, msh1_1153, msi_1535, msi_1536, \
                         msi_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1974[k] = f_10 * msh0_1151[k]
                    - f_11 * msh1_1151[k]
                    + f_3 * pc_y[k] * msi_1535[k];

        t_1975[k] = f_8 * msh0_1152[k]
                    - f_9 * msh1_1152[k]
                    + f_3 * pc_y[k] * msi_1536[k];

        t_1976[k] = f_6 * msh0_1153[k]
                    - f_7 * msh1_1153[k]
                    + f_3 * pc_y[k] * msi_1537[k];
    }

#pragma omp simd aligned(t_1977, t_1978, t_1979, pc_y, pc_z, lsi_1259, msh0_1154, msh1_1154, \
                         msi_1538, msi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1977[k] = f_4 * msh0_1154[k]
                    - f_5 * msh1_1154[k]
                    + f_3 * pc_y[k] * msi_1538[k];

        t_1978[k] = f_3 * pc_y[k] * msi_1539[k];

        t_1979[k] = f_0 * lsi_1259[k]
                    + f_1 * msh0_1154[k]
                    - f_2 * msh1_1154[k]
                    + f_3 * pc_z[k] * msi_1539[k];
    }
}

auto
compute_prim_msk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsk0, const size_t lsi,
                                                   const size_t lsk1, const size_t msh0,
                                                   const size_t msh1, const size_t msi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece7(buffer, target, pc, lsi, msh0,
                                                              msh1, msi, ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, lsk0, lsi,
                                                              lsk1, msh0, msh1, msi, ncols,
                                                              gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, lsk0,
                                                               lsi, lsk1, msh0, msh1, msi,
                                                               ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, lsk0,
                                                               lsi, lsk1, msh0, msh1, msi,
                                                               ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece12(buffer, target, pa, pc, lsk0,
                                                               lsi, lsk1, msi, ncols, gamma, p,
                                                               q);

    compute_prim_msk_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, lsk0,
                                                               lsi, lsk1, msh0, msh1, msi,
                                                               ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece14(buffer, target, pa, pc, lsk0,
                                                               lsi, lsk1, msh0, msh1, msi,
                                                               ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece15(buffer, target, pc, lsi, msh0,
                                                               msh1, msi, ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, lsk0,
                                                               lsi, lsk1, msh0, msh1, msi,
                                                               ncols, gamma, p, q);

    compute_prim_msk_three_center_electron_repulsion_0_piece17(buffer, target, pc, lsi, msh0,
                                                               msh1, msi, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
