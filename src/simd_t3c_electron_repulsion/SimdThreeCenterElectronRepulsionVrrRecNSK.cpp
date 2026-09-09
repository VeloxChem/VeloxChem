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


#include "SimdThreeCenterElectronRepulsionVrrRecNSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 4.0 / q;

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

    const auto *msk0_0 = buffer.data(msk0 + 0);
    const auto *msk0_3 = buffer.data(msk0 + 3);
    const auto *msk0_5 = buffer.data(msk0 + 5);
    const auto *msk0_6 = buffer.data(msk0 + 6);
    const auto *msk0_9 = buffer.data(msk0 + 9);
    const auto *msk0_10 = buffer.data(msk0 + 10);
    const auto *msk0_14 = buffer.data(msk0 + 14);
    const auto *msk0_15 = buffer.data(msk0 + 15);
    const auto *msk0_20 = buffer.data(msk0 + 20);
    const auto *msk0_28 = buffer.data(msk0 + 28);
    const auto *msk0_35 = buffer.data(msk0 + 35);

    const auto *msi_0 = buffer.data(msi + 0);
    const auto *msi_1 = buffer.data(msi + 1);
    const auto *msi_2 = buffer.data(msi + 2);
    const auto *msi_3 = buffer.data(msi + 3);
    const auto *msi_5 = buffer.data(msi + 5);
    const auto *msi_6 = buffer.data(msi + 6);
    const auto *msi_9 = buffer.data(msi + 9);
    const auto *msi_10 = buffer.data(msi + 10);
    const auto *msi_14 = buffer.data(msi + 14);
    const auto *msi_21 = buffer.data(msi + 21);
    const auto *msi_23 = buffer.data(msi + 23);
    const auto *msi_24 = buffer.data(msi + 24);
    const auto *msi_25 = buffer.data(msi + 25);
    const auto *msi_27 = buffer.data(msi + 27);
    const auto *msi_28 = buffer.data(msi + 28);
    const auto *msi_33 = buffer.data(msi + 33);
    const auto *msi_37 = buffer.data(msi + 37);
    const auto *msi_42 = buffer.data(msi + 42);
    const auto *msi_49 = buffer.data(msi + 49);
    const auto *msi_51 = buffer.data(msi + 51);
    const auto *msi_52 = buffer.data(msi + 52);
    const auto *msi_53 = buffer.data(msi + 53);
    const auto *msi_54 = buffer.data(msi + 54);
    const auto *msi_55 = buffer.data(msi + 55);
    const auto *msi_77 = buffer.data(msi + 77);
    const auto *msi_78 = buffer.data(msi + 78);
    const auto *msi_79 = buffer.data(msi + 79);
    const auto *msi_80 = buffer.data(msi + 80);
    const auto *msi_81 = buffer.data(msi + 81);
    const auto *msi_83 = buffer.data(msi + 83);
    const auto *msi_84 = buffer.data(msi + 84);
    const auto *msi_87 = buffer.data(msi + 87);
    const auto *msi_90 = buffer.data(msi + 90);
    const auto *msi_94 = buffer.data(msi + 94);
    const auto *msi_99 = buffer.data(msi + 99);

    const auto *msk1_0 = buffer.data(msk1 + 0);
    const auto *msk1_3 = buffer.data(msk1 + 3);
    const auto *msk1_5 = buffer.data(msk1 + 5);
    const auto *msk1_6 = buffer.data(msk1 + 6);
    const auto *msk1_9 = buffer.data(msk1 + 9);
    const auto *msk1_10 = buffer.data(msk1 + 10);
    const auto *msk1_14 = buffer.data(msk1 + 14);
    const auto *msk1_15 = buffer.data(msk1 + 15);
    const auto *msk1_20 = buffer.data(msk1 + 20);
    const auto *msk1_28 = buffer.data(msk1 + 28);
    const auto *msk1_35 = buffer.data(msk1 + 35);

    const auto *nsh0_0 = buffer.data(nsh0 + 0);
    const auto *nsh0_1 = buffer.data(nsh0 + 1);
    const auto *nsh0_2 = buffer.data(nsh0 + 2);
    const auto *nsh0_3 = buffer.data(nsh0 + 3);
    const auto *nsh0_5 = buffer.data(nsh0 + 5);
    const auto *nsh0_6 = buffer.data(nsh0 + 6);
    const auto *nsh0_8 = buffer.data(nsh0 + 8);
    const auto *nsh0_9 = buffer.data(nsh0 + 9);
    const auto *nsh0_15 = buffer.data(nsh0 + 15);
    const auto *nsh0_17 = buffer.data(nsh0 + 17);
    const auto *nsh0_18 = buffer.data(nsh0 + 18);
    const auto *nsh0_19 = buffer.data(nsh0 + 19);
    const auto *nsh0_20 = buffer.data(nsh0 + 20);
    const auto *nsh0_24 = buffer.data(nsh0 + 24);
    const auto *nsh0_27 = buffer.data(nsh0 + 27);
    const auto *nsh0_28 = buffer.data(nsh0 + 28);
    const auto *nsh0_36 = buffer.data(nsh0 + 36);
    const auto *nsh0_37 = buffer.data(nsh0 + 37);
    const auto *nsh0_38 = buffer.data(nsh0 + 38);
    const auto *nsh0_39 = buffer.data(nsh0 + 39);
    const auto *nsh0_44 = buffer.data(nsh0 + 44);
    const auto *nsh0_46 = buffer.data(nsh0 + 46);
    const auto *nsh0_47 = buffer.data(nsh0 + 47);
    const auto *nsh0_49 = buffer.data(nsh0 + 49);
    const auto *nsh0_50 = buffer.data(nsh0 + 50);
    const auto *nsh0_51 = buffer.data(nsh0 + 51);
    const auto *nsh0_58 = buffer.data(nsh0 + 58);
    const auto *nsh0_59 = buffer.data(nsh0 + 59);
    const auto *nsh0_60 = buffer.data(nsh0 + 60);
    const auto *nsh0_61 = buffer.data(nsh0 + 61);
    const auto *nsh0_62 = buffer.data(nsh0 + 62);
    const auto *nsh0_63 = buffer.data(nsh0 + 63);
    const auto *nsh0_65 = buffer.data(nsh0 + 65);
    const auto *nsh0_66 = buffer.data(nsh0 + 66);
    const auto *nsh0_68 = buffer.data(nsh0 + 68);
    const auto *nsh0_69 = buffer.data(nsh0 + 69);
    const auto *nsh0_70 = buffer.data(nsh0 + 70);
    const auto *nsh0_72 = buffer.data(nsh0 + 72);
    const auto *nsh0_73 = buffer.data(nsh0 + 73);
    const auto *nsh0_78 = buffer.data(nsh0 + 78);

    const auto *nsh1_0 = buffer.data(nsh1 + 0);
    const auto *nsh1_1 = buffer.data(nsh1 + 1);
    const auto *nsh1_2 = buffer.data(nsh1 + 2);
    const auto *nsh1_3 = buffer.data(nsh1 + 3);
    const auto *nsh1_5 = buffer.data(nsh1 + 5);
    const auto *nsh1_6 = buffer.data(nsh1 + 6);
    const auto *nsh1_8 = buffer.data(nsh1 + 8);
    const auto *nsh1_9 = buffer.data(nsh1 + 9);
    const auto *nsh1_15 = buffer.data(nsh1 + 15);
    const auto *nsh1_17 = buffer.data(nsh1 + 17);
    const auto *nsh1_18 = buffer.data(nsh1 + 18);
    const auto *nsh1_19 = buffer.data(nsh1 + 19);
    const auto *nsh1_20 = buffer.data(nsh1 + 20);
    const auto *nsh1_24 = buffer.data(nsh1 + 24);
    const auto *nsh1_27 = buffer.data(nsh1 + 27);
    const auto *nsh1_28 = buffer.data(nsh1 + 28);
    const auto *nsh1_36 = buffer.data(nsh1 + 36);
    const auto *nsh1_37 = buffer.data(nsh1 + 37);
    const auto *nsh1_38 = buffer.data(nsh1 + 38);
    const auto *nsh1_39 = buffer.data(nsh1 + 39);
    const auto *nsh1_44 = buffer.data(nsh1 + 44);
    const auto *nsh1_46 = buffer.data(nsh1 + 46);
    const auto *nsh1_47 = buffer.data(nsh1 + 47);
    const auto *nsh1_49 = buffer.data(nsh1 + 49);
    const auto *nsh1_50 = buffer.data(nsh1 + 50);
    const auto *nsh1_51 = buffer.data(nsh1 + 51);
    const auto *nsh1_58 = buffer.data(nsh1 + 58);
    const auto *nsh1_59 = buffer.data(nsh1 + 59);
    const auto *nsh1_60 = buffer.data(nsh1 + 60);
    const auto *nsh1_61 = buffer.data(nsh1 + 61);
    const auto *nsh1_62 = buffer.data(nsh1 + 62);
    const auto *nsh1_63 = buffer.data(nsh1 + 63);
    const auto *nsh1_65 = buffer.data(nsh1 + 65);
    const auto *nsh1_66 = buffer.data(nsh1 + 66);
    const auto *nsh1_68 = buffer.data(nsh1 + 68);
    const auto *nsh1_69 = buffer.data(nsh1 + 69);
    const auto *nsh1_70 = buffer.data(nsh1 + 70);
    const auto *nsh1_72 = buffer.data(nsh1 + 72);
    const auto *nsh1_73 = buffer.data(nsh1 + 73);
    const auto *nsh1_78 = buffer.data(nsh1 + 78);

    const auto *nsi_0 = buffer.data(nsi + 0);
    const auto *nsi_1 = buffer.data(nsi + 1);
    const auto *nsi_2 = buffer.data(nsi + 2);
    const auto *nsi_3 = buffer.data(nsi + 3);
    const auto *nsi_5 = buffer.data(nsi + 5);
    const auto *nsi_6 = buffer.data(nsi + 6);
    const auto *nsi_8 = buffer.data(nsi + 8);
    const auto *nsi_9 = buffer.data(nsi + 9);
    const auto *nsi_10 = buffer.data(nsi + 10);
    const auto *nsi_12 = buffer.data(nsi + 12);
    const auto *nsi_13 = buffer.data(nsi + 13);
    const auto *nsi_14 = buffer.data(nsi + 14);
    const auto *nsi_15 = buffer.data(nsi + 15);
    const auto *nsi_20 = buffer.data(nsi + 20);
    const auto *nsi_21 = buffer.data(nsi + 21);
    const auto *nsi_23 = buffer.data(nsi + 23);
    const auto *nsi_24 = buffer.data(nsi + 24);
    const auto *nsi_25 = buffer.data(nsi + 25);
    const auto *nsi_26 = buffer.data(nsi + 26);
    const auto *nsi_27 = buffer.data(nsi + 27);
    const auto *nsi_28 = buffer.data(nsi + 28);
    const auto *nsi_29 = buffer.data(nsi + 29);
    const auto *nsi_31 = buffer.data(nsi + 31);
    const auto *nsi_33 = buffer.data(nsi + 33);
    const auto *nsi_34 = buffer.data(nsi + 34);
    const auto *nsi_35 = buffer.data(nsi + 35);
    const auto *nsi_37 = buffer.data(nsi + 37);
    const auto *nsi_38 = buffer.data(nsi + 38);
    const auto *nsi_39 = buffer.data(nsi + 39);
    const auto *nsi_40 = buffer.data(nsi + 40);
    const auto *nsi_42 = buffer.data(nsi + 42);
    const auto *nsi_43 = buffer.data(nsi + 43);
    const auto *nsi_49 = buffer.data(nsi + 49);
    const auto *nsi_50 = buffer.data(nsi + 50);
    const auto *nsi_51 = buffer.data(nsi + 51);
    const auto *nsi_52 = buffer.data(nsi + 52);
    const auto *nsi_53 = buffer.data(nsi + 53);
    const auto *nsi_54 = buffer.data(nsi + 54);
    const auto *nsi_55 = buffer.data(nsi + 55);
    const auto *nsi_56 = buffer.data(nsi + 56);
    const auto *nsi_58 = buffer.data(nsi + 58);
    const auto *nsi_60 = buffer.data(nsi + 60);
    const auto *nsi_61 = buffer.data(nsi + 61);
    const auto *nsi_63 = buffer.data(nsi + 63);
    const auto *nsi_64 = buffer.data(nsi + 64);
    const auto *nsi_65 = buffer.data(nsi + 65);
    const auto *nsi_67 = buffer.data(nsi + 67);
    const auto *nsi_68 = buffer.data(nsi + 68);
    const auto *nsi_69 = buffer.data(nsi + 69);
    const auto *nsi_70 = buffer.data(nsi + 70);
    const auto *nsi_76 = buffer.data(nsi + 76);
    const auto *nsi_77 = buffer.data(nsi + 77);
    const auto *nsi_78 = buffer.data(nsi + 78);
    const auto *nsi_79 = buffer.data(nsi + 79);
    const auto *nsi_80 = buffer.data(nsi + 80);
    const auto *nsi_81 = buffer.data(nsi + 81);
    const auto *nsi_82 = buffer.data(nsi + 82);
    const auto *nsi_83 = buffer.data(nsi + 83);
    const auto *nsi_84 = buffer.data(nsi + 84);
    const auto *nsi_85 = buffer.data(nsi + 85);
    const auto *nsi_86 = buffer.data(nsi + 86);
    const auto *nsi_87 = buffer.data(nsi + 87);
    const auto *nsi_89 = buffer.data(nsi + 89);
    const auto *nsi_90 = buffer.data(nsi + 90);
    const auto *nsi_91 = buffer.data(nsi + 91);
    const auto *nsi_93 = buffer.data(nsi + 93);
    const auto *nsi_94 = buffer.data(nsi + 94);
    const auto *nsi_95 = buffer.data(nsi + 95);
    const auto *nsi_96 = buffer.data(nsi + 96);
    const auto *nsi_98 = buffer.data(nsi + 98);
    const auto *nsi_99 = buffer.data(nsi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, msi_0, nsh0_0, \
                         nsh1_0, nsi_0, nsi_1, nsi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * msi_0[k]
                 + f_1 * nsh0_0[k]
                 - f_2 * nsh1_0[k]
                 + f_3 * pc_x[k] * nsi_0[k];

        t_1[k] = f_3 * pc_y[k] * nsi_0[k];

        t_2[k] = f_3 * pc_z[k] * nsi_0[k];

        t_3[k] = f_4 * nsh0_0[k]
                 - f_5 * nsh1_0[k]
                 + f_3 * pc_y[k] * nsi_1[k];

        t_4[k] = f_3 * pc_y[k] * nsi_2[k];

        t_5[k] = f_4 * nsh0_0[k]
                 - f_5 * nsh1_0[k]
                 + f_3 * pc_z[k] * nsi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, nsh0_1, nsh0_2, nsh0_3, nsh1_1, \
                         nsh1_2, nsh1_3, nsi_3, nsi_5, nsi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * nsh0_1[k]
                 - f_7 * nsh1_1[k]
                 + f_3 * pc_y[k] * nsi_3[k];

        t_7[k] = f_3 * pc_z[k] * nsi_3[k];

        t_8[k] = f_3 * pc_y[k] * nsi_5[k];

        t_9[k] = f_6 * nsh0_2[k]
                 - f_7 * nsh1_2[k]
                 + f_3 * pc_z[k] * nsi_5[k];

        t_10[k] = f_8 * nsh0_3[k]
                  - f_9 * nsh1_3[k]
                  + f_3 * pc_y[k] * nsi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, nsh0_5, nsh0_6, \
                         nsh1_5, nsh1_6, nsi_6, nsi_8, nsi_9, nsi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * nsi_6[k];

        t_12[k] = f_4 * nsh0_5[k]
                  - f_5 * nsh1_5[k]
                  + f_3 * pc_y[k] * nsi_8[k];

        t_13[k] = f_3 * pc_y[k] * nsi_9[k];

        t_14[k] = f_8 * nsh0_5[k]
                  - f_9 * nsh1_5[k]
                  + f_3 * pc_z[k] * nsi_9[k];

        t_15[k] = f_10 * nsh0_6[k]
                  - f_11 * nsh1_6[k]
                  + f_3 * pc_y[k] * nsi_10[k];

        t_16[k] = f_3 * pc_z[k] * nsi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, nsh0_8, nsh0_9, nsh1_8, nsh1_9, \
                         nsi_12, nsi_13, nsi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * nsh0_8[k]
                  - f_7 * nsh1_8[k]
                  + f_3 * pc_y[k] * nsi_12[k];

        t_18[k] = f_4 * nsh0_9[k]
                  - f_5 * nsh1_9[k]
                  + f_3 * pc_y[k] * nsi_13[k];

        t_19[k] = f_3 * pc_y[k] * nsi_14[k];

        t_20[k] = f_10 * nsh0_9[k]
                  - f_11 * nsh1_9[k]
                  + f_3 * pc_z[k] * nsi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, msi_21, msi_23, msi_24, \
                         msi_25, nsi_15, nsi_21, nsi_23, nsi_24, \
                         nsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * msi_21[k]
                  + f_3 * pc_x[k] * nsi_21[k];

        t_22[k] = f_3 * pc_z[k] * nsi_15[k];

        t_23[k] = f_0 * msi_23[k]
                  + f_3 * pc_x[k] * nsi_23[k];

        t_24[k] = f_0 * msi_24[k]
                  + f_3 * pc_x[k] * nsi_24[k];

        t_25[k] = f_0 * msi_25[k]
                  + f_3 * pc_x[k] * nsi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, msi_27, nsh0_15, nsh1_15, \
                         nsi_20, nsi_21, nsi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * nsi_20[k];

        t_27[k] = f_0 * msi_27[k]
                  + f_3 * pc_x[k] * nsi_27[k];

        t_28[k] = f_1 * nsh0_15[k]
                  - f_2 * nsh1_15[k]
                  + f_3 * pc_y[k] * nsi_21[k];

        t_29[k] = f_3 * pc_z[k] * nsi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, nsh0_17, nsh0_18, nsh0_19, nsh1_17, nsh1_18, \
                         nsh1_19, nsi_23, nsi_24, nsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * nsh0_17[k]
                  - f_11 * nsh1_17[k]
                  + f_3 * pc_y[k] * nsi_23[k];

        t_31[k] = f_8 * nsh0_18[k]
                  - f_9 * nsh1_18[k]
                  + f_3 * pc_y[k] * nsi_24[k];

        t_32[k] = f_6 * nsh0_19[k]
                  - f_7 * nsh1_19[k]
                  + f_3 * pc_y[k] * nsi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, msk0_0, msi_0, \
                         msk1_0, nsh0_20, nsh1_20, nsi_26, nsi_27, \
                         nsi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * nsh0_20[k]
                  - f_5 * nsh1_20[k]
                  + f_3 * pc_y[k] * nsi_26[k];

        t_34[k] = f_3 * pc_y[k] * nsi_27[k];

        t_35[k] = f_1 * nsh0_20[k]
                  - f_2 * nsh1_20[k]
                  + f_3 * pc_z[k] * nsi_27[k];

        t_36[k] = pa_y[k] * msk0_0[k]
                  - f_12 * pc_y[k] * msk1_0[k];

        t_37[k] = f_13 * msi_0[k]
                  + f_3 * pc_y[k] * nsi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, msk0_3, msk0_5, msi_1, \
                         msk1_3, msk1_5, nsi_28, nsi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * nsi_28[k];

        t_39[k] = pa_y[k] * msk0_3[k]
                  + f_14 * msi_1[k]
                  - f_12 * pc_y[k] * msk1_3[k];

        t_40[k] = f_3 * pc_z[k] * nsi_29[k];

        t_41[k] = pa_y[k] * msk0_5[k]
                  - f_12 * pc_y[k] * msk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, msk0_6, msk0_9, msi_3, \
                         msi_5, msk1_6, msk1_9, nsi_31, nsi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * msk0_6[k]
                  + f_15 * msi_3[k]
                  - f_12 * pc_y[k] * msk1_6[k];

        t_43[k] = f_3 * pc_z[k] * nsi_31[k];

        t_44[k] = f_13 * msi_5[k]
                  + f_3 * pc_y[k] * nsi_33[k];

        t_45[k] = pa_y[k] * msk0_9[k]
                  - f_12 * pc_y[k] * msk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, msk0_10, msi_6, msi_9, \
                         msk1_10, nsh0_24, nsh1_24, nsi_34, nsi_35, \
                         nsi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * msk0_10[k]
                  + f_16 * msi_6[k]
                  - f_12 * pc_y[k] * msk1_10[k];

        t_47[k] = f_3 * pc_z[k] * nsi_34[k];

        t_48[k] = f_4 * nsh0_24[k]
                  - f_5 * nsh1_24[k]
                  + f_3 * pc_z[k] * nsi_35[k];

        t_49[k] = f_13 * msi_9[k]
                  + f_3 * pc_y[k] * nsi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, msk0_14, msk0_15, msi_10, \
                         msk1_14, msk1_15, nsh0_27, nsh1_27, nsi_38, \
                         nsi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * msk0_14[k]
                  - f_12 * pc_y[k] * msk1_14[k];

        t_51[k] = pa_y[k] * msk0_15[k]
                  + f_17 * msi_10[k]
                  - f_12 * pc_y[k] * msk1_15[k];

        t_52[k] = f_3 * pc_z[k] * nsi_38[k];

        t_53[k] = f_4 * nsh0_27[k]
                  - f_5 * nsh1_27[k]
                  + f_3 * pc_z[k] * nsi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, msk0_20, msi_14, msk1_20, \
                         nsh0_28, nsh1_28, nsi_40, nsi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * nsh0_28[k]
                  - f_7 * nsh1_28[k]
                  + f_3 * pc_z[k] * nsi_40[k];

        t_55[k] = f_13 * msi_14[k]
                  + f_3 * pc_y[k] * nsi_42[k];

        t_56[k] = pa_y[k] * msk0_20[k]
                  - f_12 * pc_y[k] * msk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, msi_49, msi_51, msi_52, \
                         msi_53, nsi_43, nsi_49, nsi_51, nsi_52, \
                         nsi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_18 * msi_49[k]
                  + f_3 * pc_x[k] * nsi_49[k];

        t_58[k] = f_3 * pc_z[k] * nsi_43[k];

        t_59[k] = f_18 * msi_51[k]
                  + f_3 * pc_x[k] * nsi_51[k];

        t_60[k] = f_18 * msi_52[k]
                  + f_3 * pc_x[k] * nsi_52[k];

        t_61[k] = f_18 * msi_53[k]
                  + f_3 * pc_x[k] * nsi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, msi_21, msi_54, msi_55, \
                         nsh0_36, nsh1_36, nsi_49, nsi_54, nsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_18 * msi_54[k]
                  + f_3 * pc_x[k] * nsi_54[k];

        t_63[k] = f_18 * msi_55[k]
                  + f_3 * pc_x[k] * nsi_55[k];

        t_64[k] = f_13 * msi_21[k]
                  + f_1 * nsh0_36[k]
                  - f_2 * nsh1_36[k]
                  + f_3 * pc_y[k] * nsi_49[k];

        t_65[k] = f_3 * pc_z[k] * nsi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, nsh0_36, nsh0_37, nsh0_38, nsh1_36, nsh1_37, \
                         nsh1_38, nsi_50, nsi_51, nsi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * nsh0_36[k]
                  - f_5 * nsh1_36[k]
                  + f_3 * pc_z[k] * nsi_50[k];

        t_67[k] = f_6 * nsh0_37[k]
                  - f_7 * nsh1_37[k]
                  + f_3 * pc_z[k] * nsi_51[k];

        t_68[k] = f_8 * nsh0_38[k]
                  - f_9 * nsh1_38[k]
                  + f_3 * pc_z[k] * nsi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, msk0_35, msi_27, msk1_35, \
                         nsh0_39, nsh1_39, nsi_53, nsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * nsh0_39[k]
                  - f_11 * nsh1_39[k]
                  + f_3 * pc_z[k] * nsi_53[k];

        t_70[k] = f_13 * msi_27[k]
                  + f_3 * pc_y[k] * nsi_55[k];

        t_71[k] = pa_y[k] * msk0_35[k]
                  - f_12 * pc_y[k] * msk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, msk0_0, msk0_3, \
                         msi_0, msk1_0, msk1_3, nsi_56, nsi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * msk0_0[k]
                  - f_12 * pc_z[k] * msk1_0[k];

        t_73[k] = f_3 * pc_y[k] * nsi_56[k];

        t_74[k] = f_13 * msi_0[k]
                  + f_3 * pc_z[k] * nsi_56[k];

        t_75[k] = pa_z[k] * msk0_3[k]
                  - f_12 * pc_z[k] * msk1_3[k];

        t_76[k] = f_3 * pc_y[k] * nsi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, msk0_5, msk0_6, msi_2, \
                         msk1_5, msk1_6, nsh0_44, nsh1_44, nsi_60, \
                         nsi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * msk0_5[k]
                  + f_14 * msi_2[k]
                  - f_12 * pc_z[k] * msk1_5[k];

        t_78[k] = pa_z[k] * msk0_6[k]
                  - f_12 * pc_z[k] * msk1_6[k];

        t_79[k] = f_4 * nsh0_44[k]
                  - f_5 * nsh1_44[k]
                  + f_3 * pc_y[k] * nsi_60[k];

        t_80[k] = f_3 * pc_y[k] * nsi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, msk0_9, msk0_10, msi_5, msk1_9, \
                         msk1_10, nsh0_46, nsh1_46, nsi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * msk0_9[k]
                  + f_15 * msi_5[k]
                  - f_12 * pc_z[k] * msk1_9[k];

        t_82[k] = pa_z[k] * msk0_10[k]
                  - f_12 * pc_z[k] * msk1_10[k];

        t_83[k] = f_6 * nsh0_46[k]
                  - f_7 * nsh1_46[k]
                  + f_3 * pc_y[k] * nsi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, msk0_14, msk0_15, msi_9, \
                         msk1_14, msk1_15, nsh0_47, nsh1_47, nsi_64, \
                         nsi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * nsh0_47[k]
                  - f_5 * nsh1_47[k]
                  + f_3 * pc_y[k] * nsi_64[k];

        t_85[k] = f_3 * pc_y[k] * nsi_65[k];

        t_86[k] = pa_z[k] * msk0_14[k]
                  + f_16 * msi_9[k]
                  - f_12 * pc_z[k] * msk1_14[k];

        t_87[k] = pa_z[k] * msk0_15[k]
                  - f_12 * pc_z[k] * msk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, nsh0_49, nsh0_50, nsh0_51, nsh1_49, \
                         nsh1_50, nsh1_51, nsi_67, nsi_68, nsi_69, \
                         nsi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * nsh0_49[k]
                  - f_9 * nsh1_49[k]
                  + f_3 * pc_y[k] * nsi_67[k];

        t_89[k] = f_6 * nsh0_50[k]
                  - f_7 * nsh1_50[k]
                  + f_3 * pc_y[k] * nsi_68[k];

        t_90[k] = f_4 * nsh0_51[k]
                  - f_5 * nsh1_51[k]
                  + f_3 * pc_y[k] * nsi_69[k];

        t_91[k] = f_3 * pc_y[k] * nsi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, msk0_20, msi_14, msi_77, \
                         msi_78, msi_79, msk1_20, nsi_77, nsi_78, \
                         nsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * msk0_20[k]
                  + f_17 * msi_14[k]
                  - f_12 * pc_z[k] * msk1_20[k];

        t_93[k] = f_18 * msi_77[k]
                  + f_3 * pc_x[k] * nsi_77[k];

        t_94[k] = f_18 * msi_78[k]
                  + f_3 * pc_x[k] * nsi_78[k];

        t_95[k] = f_18 * msi_79[k]
                  + f_3 * pc_x[k] * nsi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, msi_80, msi_81, msi_83, nsi_76, \
                         nsi_80, nsi_81, nsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_18 * msi_80[k]
                  + f_3 * pc_x[k] * nsi_80[k];

        t_97[k] = f_18 * msi_81[k]
                  + f_3 * pc_x[k] * nsi_81[k];

        t_98[k] = f_3 * pc_y[k] * nsi_76[k];

        t_99[k] = f_18 * msi_83[k]
                  + f_3 * pc_x[k] * nsi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, msk0_28, msk1_28, nsh0_58, \
                         nsh0_59, nsh1_58, nsh1_59, nsi_78, nsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * msk0_28[k]
                   - f_12 * pc_z[k] * msk1_28[k];

        t_101[k] = f_19 * nsh0_58[k]
                   - f_20 * nsh1_58[k]
                   + f_3 * pc_y[k] * nsi_78[k];

        t_102[k] = f_10 * nsh0_59[k]
                   - f_11 * nsh1_59[k]
                   + f_3 * pc_y[k] * nsi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, nsh0_60, nsh0_61, nsh0_62, nsh1_60, \
                         nsh1_61, nsh1_62, nsi_80, nsi_81, nsi_82, \
                         nsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * nsh0_60[k]
                   - f_9 * nsh1_60[k]
                   + f_3 * pc_y[k] * nsi_80[k];

        t_104[k] = f_6 * nsh0_61[k]
                   - f_7 * nsh1_61[k]
                   + f_3 * pc_y[k] * nsi_81[k];

        t_105[k] = f_4 * nsh0_62[k]
                   - f_5 * nsh1_62[k]
                   + f_3 * pc_y[k] * nsi_82[k];

        t_106[k] = f_3 * pc_y[k] * nsi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, msi_27, msi_28, msi_84, \
                         nsh0_62, nsh0_63, nsh1_62, nsh1_63, nsi_83, \
                         nsi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * msi_27[k]
                   + f_1 * nsh0_62[k]
                   - f_2 * nsh1_62[k]
                   + f_3 * pc_z[k] * nsi_83[k];

        t_108[k] = f_21 * msi_84[k]
                   + f_1 * nsh0_63[k]
                   - f_2 * nsh1_63[k]
                   + f_3 * pc_x[k] * nsi_84[k];

        t_109[k] = f_14 * msi_28[k]
                   + f_3 * pc_y[k] * nsi_84[k];

        t_110[k] = f_3 * pc_z[k] * nsi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, msi_87, nsh0_63, nsh0_66, nsh1_63, \
                         nsh1_66, nsi_85, nsi_86, nsi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_21 * msi_87[k]
                   + f_10 * nsh0_66[k]
                   - f_11 * nsh1_66[k]
                   + f_3 * pc_x[k] * nsi_87[k];

        t_112[k] = f_3 * pc_z[k] * nsi_85[k];

        t_113[k] = f_4 * nsh0_63[k]
                   - f_5 * nsh1_63[k]
                   + f_3 * pc_z[k] * nsi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, msi_33, msi_90, \
                         nsh0_65, nsh0_69, nsh1_65, nsh1_69, nsi_87, nsi_89, \
                         nsi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_21 * msi_90[k]
                   + f_8 * nsh0_69[k]
                   - f_9 * nsh1_69[k]
                   + f_3 * pc_x[k] * nsi_90[k];

        t_115[k] = f_3 * pc_z[k] * nsi_87[k];

        t_116[k] = f_14 * msi_33[k]
                   + f_3 * pc_y[k] * nsi_89[k];

        t_117[k] = f_6 * nsh0_65[k]
                   - f_7 * nsh1_65[k]
                   + f_3 * pc_z[k] * nsi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, msi_94, nsh0_66, nsh0_73, nsh1_66, \
                         nsh1_73, nsi_90, nsi_91, nsi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_21 * msi_94[k]
                   + f_6 * nsh0_73[k]
                   - f_7 * nsh1_73[k]
                   + f_3 * pc_x[k] * nsi_94[k];

        t_119[k] = f_3 * pc_z[k] * nsi_90[k];

        t_120[k] = f_4 * nsh0_66[k]
                   - f_5 * nsh1_66[k]
                   + f_3 * pc_z[k] * nsi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, msi_37, msi_99, \
                         nsh0_68, nsh0_78, nsh1_68, nsh1_78, nsi_93, nsi_94, \
                         nsi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * msi_37[k]
                   + f_3 * pc_y[k] * nsi_93[k];

        t_122[k] = f_8 * nsh0_68[k]
                   - f_9 * nsh1_68[k]
                   + f_3 * pc_z[k] * nsi_93[k];

        t_123[k] = f_21 * msi_99[k]
                   + f_4 * nsh0_78[k]
                   - f_5 * nsh1_78[k]
                   + f_3 * pc_x[k] * nsi_99[k];

        t_124[k] = f_3 * pc_z[k] * nsi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, msi_42, nsh0_69, nsh0_70, \
                         nsh0_72, nsh1_69, nsh1_70, nsh1_72, nsi_95, nsi_96, \
                         nsi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * nsh0_69[k]
                   - f_5 * nsh1_69[k]
                   + f_3 * pc_z[k] * nsi_95[k];

        t_126[k] = f_6 * nsh0_70[k]
                   - f_7 * nsh1_70[k]
                   + f_3 * pc_z[k] * nsi_96[k];

        t_127[k] = f_14 * msi_42[k]
                   + f_3 * pc_y[k] * nsi_98[k];

        t_128[k] = f_10 * nsh0_72[k]
                   - f_11 * nsh1_72[k]
                   + f_3 * pc_z[k] * nsi_98[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *msk0_39 = buffer.data(msk0 + 39);
    const auto *msk0_42 = buffer.data(msk0 + 42);
    const auto *msk0_46 = buffer.data(msk0 + 46);
    const auto *msk0_51 = buffer.data(msk0 + 51);
    const auto *msk0_64 = buffer.data(msk0 + 64);
    const auto *msk0_72 = buffer.data(msk0 + 72);
    const auto *msk0_77 = buffer.data(msk0 + 77);
    const auto *msk0_81 = buffer.data(msk0 + 81);
    const auto *msk0_84 = buffer.data(msk0 + 84);
    const auto *msk0_86 = buffer.data(msk0 + 86);
    const auto *msk0_89 = buffer.data(msk0 + 89);
    const auto *msk0_90 = buffer.data(msk0 + 90);
    const auto *msk0_92 = buffer.data(msk0 + 92);
    const auto *msk0_107 = buffer.data(msk0 + 107);

    const auto *msi_28 = buffer.data(msi + 28);
    const auto *msi_31 = buffer.data(msi + 31);
    const auto *msi_34 = buffer.data(msi + 34);
    const auto *msi_38 = buffer.data(msi + 38);
    const auto *msi_49 = buffer.data(msi + 49);
    const auto *msi_55 = buffer.data(msi + 55);
    const auto *msi_56 = buffer.data(msi + 56);
    const auto *msi_58 = buffer.data(msi + 58);
    const auto *msi_61 = buffer.data(msi + 61);
    const auto *msi_64 = buffer.data(msi + 64);
    const auto *msi_65 = buffer.data(msi + 65);
    const auto *msi_68 = buffer.data(msi + 68);
    const auto *msi_69 = buffer.data(msi + 69);
    const auto *msi_70 = buffer.data(msi + 70);
    const auto *msi_79 = buffer.data(msi + 79);
    const auto *msi_80 = buffer.data(msi + 80);
    const auto *msi_81 = buffer.data(msi + 81);
    const auto *msi_82 = buffer.data(msi + 82);
    const auto *msi_83 = buffer.data(msi + 83);
    const auto *msi_84 = buffer.data(msi + 84);
    const auto *msi_89 = buffer.data(msi + 89);
    const auto *msi_93 = buffer.data(msi + 93);
    const auto *msi_98 = buffer.data(msi + 98);
    const auto *msi_105 = buffer.data(msi + 105);
    const auto *msi_107 = buffer.data(msi + 107);
    const auto *msi_108 = buffer.data(msi + 108);
    const auto *msi_109 = buffer.data(msi + 109);
    const auto *msi_110 = buffer.data(msi + 110);
    const auto *msi_111 = buffer.data(msi + 111);
    const auto *msi_133 = buffer.data(msi + 133);
    const auto *msi_134 = buffer.data(msi + 134);
    const auto *msi_135 = buffer.data(msi + 135);
    const auto *msi_136 = buffer.data(msi + 136);
    const auto *msi_137 = buffer.data(msi + 137);
    const auto *msi_138 = buffer.data(msi + 138);
    const auto *msi_139 = buffer.data(msi + 139);
    const auto *msi_140 = buffer.data(msi + 140);
    const auto *msi_145 = buffer.data(msi + 145);
    const auto *msi_149 = buffer.data(msi + 149);
    const auto *msi_154 = buffer.data(msi + 154);
    const auto *msi_160 = buffer.data(msi + 160);
    const auto *msi_161 = buffer.data(msi + 161);
    const auto *msi_162 = buffer.data(msi + 162);
    const auto *msi_163 = buffer.data(msi + 163);
    const auto *msi_164 = buffer.data(msi + 164);
    const auto *msi_165 = buffer.data(msi + 165);
    const auto *msi_167 = buffer.data(msi + 167);
    const auto *msi_168 = buffer.data(msi + 168);
    const auto *msi_171 = buffer.data(msi + 171);
    const auto *msi_174 = buffer.data(msi + 174);
    const auto *msi_178 = buffer.data(msi + 178);
    const auto *msi_183 = buffer.data(msi + 183);
    const auto *msi_189 = buffer.data(msi + 189);
    const auto *msi_191 = buffer.data(msi + 191);
    const auto *msi_192 = buffer.data(msi + 192);
    const auto *msi_193 = buffer.data(msi + 193);
    const auto *msi_194 = buffer.data(msi + 194);
    const auto *msi_195 = buffer.data(msi + 195);

    const auto *msk1_39 = buffer.data(msk1 + 39);
    const auto *msk1_42 = buffer.data(msk1 + 42);
    const auto *msk1_46 = buffer.data(msk1 + 46);
    const auto *msk1_51 = buffer.data(msk1 + 51);
    const auto *msk1_64 = buffer.data(msk1 + 64);
    const auto *msk1_72 = buffer.data(msk1 + 72);
    const auto *msk1_77 = buffer.data(msk1 + 77);
    const auto *msk1_81 = buffer.data(msk1 + 81);
    const auto *msk1_84 = buffer.data(msk1 + 84);
    const auto *msk1_86 = buffer.data(msk1 + 86);
    const auto *msk1_89 = buffer.data(msk1 + 89);
    const auto *msk1_90 = buffer.data(msk1 + 90);
    const auto *msk1_92 = buffer.data(msk1 + 92);
    const auto *msk1_107 = buffer.data(msk1 + 107);

    const auto *nsh0_78 = buffer.data(nsh0 + 78);
    const auto *nsh0_79 = buffer.data(nsh0 + 79);
    const auto *nsh0_80 = buffer.data(nsh0 + 80);
    const auto *nsh0_81 = buffer.data(nsh0 + 81);
    const auto *nsh0_83 = buffer.data(nsh0 + 83);
    const auto *nsh0_101 = buffer.data(nsh0 + 101);
    const auto *nsh0_102 = buffer.data(nsh0 + 102);
    const auto *nsh0_103 = buffer.data(nsh0 + 103);
    const auto *nsh0_104 = buffer.data(nsh0 + 104);
    const auto *nsh0_105 = buffer.data(nsh0 + 105);
    const auto *nsh0_106 = buffer.data(nsh0 + 106);
    const auto *nsh0_107 = buffer.data(nsh0 + 107);
    const auto *nsh0_108 = buffer.data(nsh0 + 108);
    const auto *nsh0_109 = buffer.data(nsh0 + 109);
    const auto *nsh0_110 = buffer.data(nsh0 + 110);
    const auto *nsh0_111 = buffer.data(nsh0 + 111);
    const auto *nsh0_112 = buffer.data(nsh0 + 112);
    const auto *nsh0_113 = buffer.data(nsh0 + 113);
    const auto *nsh0_114 = buffer.data(nsh0 + 114);
    const auto *nsh0_119 = buffer.data(nsh0 + 119);
    const auto *nsh0_120 = buffer.data(nsh0 + 120);
    const auto *nsh0_121 = buffer.data(nsh0 + 121);
    const auto *nsh0_122 = buffer.data(nsh0 + 122);
    const auto *nsh0_123 = buffer.data(nsh0 + 123);
    const auto *nsh0_124 = buffer.data(nsh0 + 124);
    const auto *nsh0_125 = buffer.data(nsh0 + 125);
    const auto *nsh0_126 = buffer.data(nsh0 + 126);
    const auto *nsh0_128 = buffer.data(nsh0 + 128);
    const auto *nsh0_129 = buffer.data(nsh0 + 129);
    const auto *nsh0_131 = buffer.data(nsh0 + 131);
    const auto *nsh0_132 = buffer.data(nsh0 + 132);
    const auto *nsh0_133 = buffer.data(nsh0 + 133);
    const auto *nsh0_135 = buffer.data(nsh0 + 135);
    const auto *nsh0_136 = buffer.data(nsh0 + 136);
    const auto *nsh0_141 = buffer.data(nsh0 + 141);
    const auto *nsh0_142 = buffer.data(nsh0 + 142);
    const auto *nsh0_143 = buffer.data(nsh0 + 143);
    const auto *nsh0_144 = buffer.data(nsh0 + 144);

    const auto *nsh1_78 = buffer.data(nsh1 + 78);
    const auto *nsh1_79 = buffer.data(nsh1 + 79);
    const auto *nsh1_80 = buffer.data(nsh1 + 80);
    const auto *nsh1_81 = buffer.data(nsh1 + 81);
    const auto *nsh1_83 = buffer.data(nsh1 + 83);
    const auto *nsh1_101 = buffer.data(nsh1 + 101);
    const auto *nsh1_102 = buffer.data(nsh1 + 102);
    const auto *nsh1_103 = buffer.data(nsh1 + 103);
    const auto *nsh1_104 = buffer.data(nsh1 + 104);
    const auto *nsh1_105 = buffer.data(nsh1 + 105);
    const auto *nsh1_106 = buffer.data(nsh1 + 106);
    const auto *nsh1_107 = buffer.data(nsh1 + 107);
    const auto *nsh1_108 = buffer.data(nsh1 + 108);
    const auto *nsh1_109 = buffer.data(nsh1 + 109);
    const auto *nsh1_110 = buffer.data(nsh1 + 110);
    const auto *nsh1_111 = buffer.data(nsh1 + 111);
    const auto *nsh1_112 = buffer.data(nsh1 + 112);
    const auto *nsh1_113 = buffer.data(nsh1 + 113);
    const auto *nsh1_114 = buffer.data(nsh1 + 114);
    const auto *nsh1_119 = buffer.data(nsh1 + 119);
    const auto *nsh1_120 = buffer.data(nsh1 + 120);
    const auto *nsh1_121 = buffer.data(nsh1 + 121);
    const auto *nsh1_122 = buffer.data(nsh1 + 122);
    const auto *nsh1_123 = buffer.data(nsh1 + 123);
    const auto *nsh1_124 = buffer.data(nsh1 + 124);
    const auto *nsh1_125 = buffer.data(nsh1 + 125);
    const auto *nsh1_126 = buffer.data(nsh1 + 126);
    const auto *nsh1_128 = buffer.data(nsh1 + 128);
    const auto *nsh1_129 = buffer.data(nsh1 + 129);
    const auto *nsh1_131 = buffer.data(nsh1 + 131);
    const auto *nsh1_132 = buffer.data(nsh1 + 132);
    const auto *nsh1_133 = buffer.data(nsh1 + 133);
    const auto *nsh1_135 = buffer.data(nsh1 + 135);
    const auto *nsh1_136 = buffer.data(nsh1 + 136);
    const auto *nsh1_141 = buffer.data(nsh1 + 141);
    const auto *nsh1_142 = buffer.data(nsh1 + 142);
    const auto *nsh1_143 = buffer.data(nsh1 + 143);
    const auto *nsh1_144 = buffer.data(nsh1 + 144);

    const auto *nsi_99 = buffer.data(nsi + 99);
    const auto *nsi_105 = buffer.data(nsi + 105);
    const auto *nsi_106 = buffer.data(nsi + 106);
    const auto *nsi_107 = buffer.data(nsi + 107);
    const auto *nsi_108 = buffer.data(nsi + 108);
    const auto *nsi_109 = buffer.data(nsi + 109);
    const auto *nsi_110 = buffer.data(nsi + 110);
    const auto *nsi_111 = buffer.data(nsi + 111);
    const auto *nsi_112 = buffer.data(nsi + 112);
    const auto *nsi_114 = buffer.data(nsi + 114);
    const auto *nsi_115 = buffer.data(nsi + 115);
    const auto *nsi_117 = buffer.data(nsi + 117);
    const auto *nsi_118 = buffer.data(nsi + 118);
    const auto *nsi_121 = buffer.data(nsi + 121);
    const auto *nsi_122 = buffer.data(nsi + 122);
    const auto *nsi_126 = buffer.data(nsi + 126);
    const auto *nsi_133 = buffer.data(nsi + 133);
    const auto *nsi_134 = buffer.data(nsi + 134);
    const auto *nsi_135 = buffer.data(nsi + 135);
    const auto *nsi_136 = buffer.data(nsi + 136);
    const auto *nsi_137 = buffer.data(nsi + 137);
    const auto *nsi_138 = buffer.data(nsi + 138);
    const auto *nsi_139 = buffer.data(nsi + 139);
    const auto *nsi_140 = buffer.data(nsi + 140);
    const auto *nsi_141 = buffer.data(nsi + 141);
    const auto *nsi_142 = buffer.data(nsi + 142);
    const auto *nsi_143 = buffer.data(nsi + 143);
    const auto *nsi_144 = buffer.data(nsi + 144);
    const auto *nsi_145 = buffer.data(nsi + 145);
    const auto *nsi_146 = buffer.data(nsi + 146);
    const auto *nsi_147 = buffer.data(nsi + 147);
    const auto *nsi_148 = buffer.data(nsi + 148);
    const auto *nsi_149 = buffer.data(nsi + 149);
    const auto *nsi_150 = buffer.data(nsi + 150);
    const auto *nsi_151 = buffer.data(nsi + 151);
    const auto *nsi_152 = buffer.data(nsi + 152);
    const auto *nsi_153 = buffer.data(nsi + 153);
    const auto *nsi_154 = buffer.data(nsi + 154);
    const auto *nsi_160 = buffer.data(nsi + 160);
    const auto *nsi_161 = buffer.data(nsi + 161);
    const auto *nsi_162 = buffer.data(nsi + 162);
    const auto *nsi_163 = buffer.data(nsi + 163);
    const auto *nsi_164 = buffer.data(nsi + 164);
    const auto *nsi_165 = buffer.data(nsi + 165);
    const auto *nsi_166 = buffer.data(nsi + 166);
    const auto *nsi_167 = buffer.data(nsi + 167);
    const auto *nsi_168 = buffer.data(nsi + 168);
    const auto *nsi_169 = buffer.data(nsi + 169);
    const auto *nsi_170 = buffer.data(nsi + 170);
    const auto *nsi_171 = buffer.data(nsi + 171);
    const auto *nsi_173 = buffer.data(nsi + 173);
    const auto *nsi_174 = buffer.data(nsi + 174);
    const auto *nsi_175 = buffer.data(nsi + 175);
    const auto *nsi_177 = buffer.data(nsi + 177);
    const auto *nsi_178 = buffer.data(nsi + 178);
    const auto *nsi_179 = buffer.data(nsi + 179);
    const auto *nsi_180 = buffer.data(nsi + 180);
    const auto *nsi_182 = buffer.data(nsi + 182);
    const auto *nsi_183 = buffer.data(nsi + 183);
    const auto *nsi_189 = buffer.data(nsi + 189);
    const auto *nsi_190 = buffer.data(nsi + 190);
    const auto *nsi_191 = buffer.data(nsi + 191);
    const auto *nsi_192 = buffer.data(nsi + 192);
    const auto *nsi_193 = buffer.data(nsi + 193);
    const auto *nsi_194 = buffer.data(nsi + 194);
    const auto *nsi_195 = buffer.data(nsi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, msi_105, msi_107, \
                         msi_108, msi_109, nsi_99, nsi_105, nsi_107, nsi_108, \
                         nsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_21 * msi_105[k]
                   + f_3 * pc_x[k] * nsi_105[k];

        t_130[k] = f_3 * pc_z[k] * nsi_99[k];

        t_131[k] = f_21 * msi_107[k]
                   + f_3 * pc_x[k] * nsi_107[k];

        t_132[k] = f_21 * msi_108[k]
                   + f_3 * pc_x[k] * nsi_108[k];

        t_133[k] = f_21 * msi_109[k]
                   + f_3 * pc_x[k] * nsi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, msi_49, msi_110, \
                         msi_111, nsh0_78, nsh1_78, nsi_105, nsi_110, \
                         nsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_21 * msi_110[k]
                   + f_3 * pc_x[k] * nsi_110[k];

        t_135[k] = f_21 * msi_111[k]
                   + f_3 * pc_x[k] * nsi_111[k];

        t_136[k] = f_14 * msi_49[k]
                   + f_1 * nsh0_78[k]
                   - f_2 * nsh1_78[k]
                   + f_3 * pc_y[k] * nsi_105[k];

        t_137[k] = f_3 * pc_z[k] * nsi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, nsh0_78, nsh0_79, nsh0_80, nsh1_78, \
                         nsh1_79, nsh1_80, nsi_106, nsi_107, nsi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * nsh0_78[k]
                   - f_5 * nsh1_78[k]
                   + f_3 * pc_z[k] * nsi_106[k];

        t_139[k] = f_6 * nsh0_79[k]
                   - f_7 * nsh1_79[k]
                   + f_3 * pc_z[k] * nsi_107[k];

        t_140[k] = f_8 * nsh0_80[k]
                   - f_9 * nsh1_80[k]
                   + f_3 * pc_z[k] * nsi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, msk0_72, msi_55, \
                         msk1_72, nsh0_81, nsh0_83, nsh1_81, nsh1_83, nsi_109, \
                         nsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * nsh0_81[k]
                   - f_11 * nsh1_81[k]
                   + f_3 * pc_z[k] * nsi_109[k];

        t_142[k] = f_14 * msi_55[k]
                   + f_3 * pc_y[k] * nsi_111[k];

        t_143[k] = f_1 * nsh0_83[k]
                   - f_2 * nsh1_83[k]
                   + f_3 * pc_z[k] * nsi_111[k];

        t_144[k] = pa_y[k] * msk0_72[k]
                   - f_12 * pc_y[k] * msk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, msk0_39, msi_28, \
                         msi_56, msi_58, msk1_39, nsi_112, nsi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * msi_56[k]
                   + f_3 * pc_y[k] * nsi_112[k];

        t_146[k] = f_13 * msi_28[k]
                   + f_3 * pc_z[k] * nsi_112[k];

        t_147[k] = pa_z[k] * msk0_39[k]
                   - f_12 * pc_z[k] * msk1_39[k];

        t_148[k] = f_13 * msi_58[k]
                   + f_3 * pc_y[k] * nsi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, msk0_42, msk0_77, \
                         msi_31, msi_61, msk1_42, msk1_77, nsi_115, \
                         nsi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * msk0_77[k]
                   - f_12 * pc_y[k] * msk1_77[k];

        t_150[k] = pa_z[k] * msk0_42[k]
                   - f_12 * pc_z[k] * msk1_42[k];

        t_151[k] = f_13 * msi_31[k]
                   + f_3 * pc_z[k] * nsi_115[k];

        t_152[k] = f_13 * msi_61[k]
                   + f_3 * pc_y[k] * nsi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, msk0_46, msk0_81, \
                         msi_34, msk1_46, msk1_81, nsi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * msk0_81[k]
                   - f_12 * pc_y[k] * msk1_81[k];

        t_154[k] = pa_z[k] * msk0_46[k]
                   - f_12 * pc_z[k] * msk1_46[k];

        t_155[k] = f_13 * msi_34[k]
                   + f_3 * pc_z[k] * nsi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, msk0_84, msk0_86, msi_64, msi_65, \
                         msk1_84, msk1_86, nsi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * msk0_84[k]
                   + f_14 * msi_64[k]
                   - f_12 * pc_y[k] * msk1_84[k];

        t_157[k] = f_13 * msi_65[k]
                   + f_3 * pc_y[k] * nsi_121[k];

        t_158[k] = pa_y[k] * msk0_86[k]
                   - f_12 * pc_y[k] * msk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, msk0_51, msk0_89, \
                         msi_38, msi_68, msk1_51, msk1_89, nsi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * msk0_51[k]
                   - f_12 * pc_z[k] * msk1_51[k];

        t_160[k] = f_13 * msi_38[k]
                   + f_3 * pc_z[k] * nsi_122[k];

        t_161[k] = pa_y[k] * msk0_89[k]
                   + f_15 * msi_68[k]
                   - f_12 * pc_y[k] * msk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, msk0_90, msk0_92, \
                         msi_69, msi_70, msi_133, msk1_90, msk1_92, nsi_126, \
                         nsi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * msk0_90[k]
                   + f_14 * msi_69[k]
                   - f_12 * pc_y[k] * msk1_90[k];

        t_163[k] = f_13 * msi_70[k]
                   + f_3 * pc_y[k] * nsi_126[k];

        t_164[k] = pa_y[k] * msk0_92[k]
                   - f_12 * pc_y[k] * msk1_92[k];

        t_165[k] = f_21 * msi_133[k]
                   + f_3 * pc_x[k] * nsi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, msi_134, msi_135, msi_136, \
                         msi_137, msi_138, nsi_134, nsi_135, nsi_136, nsi_137, \
                         nsi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_21 * msi_134[k]
                   + f_3 * pc_x[k] * nsi_134[k];

        t_167[k] = f_21 * msi_135[k]
                   + f_3 * pc_x[k] * nsi_135[k];

        t_168[k] = f_21 * msi_136[k]
                   + f_3 * pc_x[k] * nsi_136[k];

        t_169[k] = f_21 * msi_137[k]
                   + f_3 * pc_x[k] * nsi_137[k];

        t_170[k] = f_21 * msi_138[k]
                   + f_3 * pc_x[k] * nsi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, msk0_64, msi_49, msi_139, \
                         msk1_64, nsi_133, nsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_21 * msi_139[k]
                   + f_3 * pc_x[k] * nsi_139[k];

        t_172[k] = pa_z[k] * msk0_64[k]
                   - f_12 * pc_z[k] * msk1_64[k];

        t_173[k] = f_13 * msi_49[k]
                   + f_3 * pc_z[k] * nsi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, msi_79, msi_80, msi_81, nsh0_101, \
                         nsh0_102, nsh0_103, nsh1_101, nsh1_102, nsh1_103, nsi_135, nsi_136, \
                         nsi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * msi_79[k]
                   + f_10 * nsh0_101[k]
                   - f_11 * nsh1_101[k]
                   + f_3 * pc_y[k] * nsi_135[k];

        t_175[k] = f_13 * msi_80[k]
                   + f_8 * nsh0_102[k]
                   - f_9 * nsh1_102[k]
                   + f_3 * pc_y[k] * nsi_136[k];

        t_176[k] = f_13 * msi_81[k]
                   + f_6 * nsh0_103[k]
                   - f_7 * nsh1_103[k]
                   + f_3 * pc_y[k] * nsi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, msk0_107, msi_82, msi_83, msk1_107, \
                         nsh0_104, nsh1_104, nsi_138, nsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * msi_82[k]
                   + f_4 * nsh0_104[k]
                   - f_5 * nsh1_104[k]
                   + f_3 * pc_y[k] * nsi_138[k];

        t_178[k] = f_13 * msi_83[k]
                   + f_3 * pc_y[k] * nsi_139[k];

        t_179[k] = pa_y[k] * msk0_107[k]
                   - f_12 * pc_y[k] * msk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, msi_56, msi_140, \
                         nsh0_105, nsh1_105, nsi_140, nsi_141, \
                         nsi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_21 * msi_140[k]
                   + f_1 * nsh0_105[k]
                   - f_2 * nsh1_105[k]
                   + f_3 * pc_x[k] * nsi_140[k];

        t_181[k] = f_3 * pc_y[k] * nsi_140[k];

        t_182[k] = f_14 * msi_56[k]
                   + f_3 * pc_z[k] * nsi_140[k];

        t_183[k] = f_4 * nsh0_105[k]
                   - f_5 * nsh1_105[k]
                   + f_3 * pc_y[k] * nsi_141[k];

        t_184[k] = f_3 * pc_y[k] * nsi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, msi_145, nsh0_106, nsh0_107, \
                         nsh0_110, nsh1_106, nsh1_107, nsh1_110, nsi_143, nsi_144, \
                         nsi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_21 * msi_145[k]
                   + f_10 * nsh0_110[k]
                   - f_11 * nsh1_110[k]
                   + f_3 * pc_x[k] * nsi_145[k];

        t_186[k] = f_6 * nsh0_106[k]
                   - f_7 * nsh1_106[k]
                   + f_3 * pc_y[k] * nsi_143[k];

        t_187[k] = f_4 * nsh0_107[k]
                   - f_5 * nsh1_107[k]
                   + f_3 * pc_y[k] * nsi_144[k];

        t_188[k] = f_3 * pc_y[k] * nsi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, msi_149, nsh0_108, nsh0_109, \
                         nsh0_114, nsh1_108, nsh1_109, nsh1_114, nsi_146, nsi_147, \
                         nsi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_21 * msi_149[k]
                   + f_8 * nsh0_114[k]
                   - f_9 * nsh1_114[k]
                   + f_3 * pc_x[k] * nsi_149[k];

        t_190[k] = f_8 * nsh0_108[k]
                   - f_9 * nsh1_108[k]
                   + f_3 * pc_y[k] * nsi_146[k];

        t_191[k] = f_6 * nsh0_109[k]
                   - f_7 * nsh1_109[k]
                   + f_3 * pc_y[k] * nsi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, msi_154, nsh0_110, nsh0_119, \
                         nsh1_110, nsh1_119, nsi_148, nsi_149, \
                         nsi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * nsh0_110[k]
                   - f_5 * nsh1_110[k]
                   + f_3 * pc_y[k] * nsi_148[k];

        t_193[k] = f_3 * pc_y[k] * nsi_149[k];

        t_194[k] = f_21 * msi_154[k]
                   + f_6 * nsh0_119[k]
                   - f_7 * nsh1_119[k]
                   + f_3 * pc_x[k] * nsi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, nsh0_111, nsh0_112, nsh0_113, nsh1_111, \
                         nsh1_112, nsh1_113, nsi_150, nsi_151, \
                         nsi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * nsh0_111[k]
                   - f_11 * nsh1_111[k]
                   + f_3 * pc_y[k] * nsi_150[k];

        t_196[k] = f_8 * nsh0_112[k]
                   - f_9 * nsh1_112[k]
                   + f_3 * pc_y[k] * nsi_151[k];

        t_197[k] = f_6 * nsh0_113[k]
                   - f_7 * nsh1_113[k]
                   + f_3 * pc_y[k] * nsi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, msi_160, msi_161, nsh0_114, \
                         nsh0_125, nsh1_114, nsh1_125, nsi_153, nsi_154, nsi_160, \
                         nsi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * nsh0_114[k]
                   - f_5 * nsh1_114[k]
                   + f_3 * pc_y[k] * nsi_153[k];

        t_199[k] = f_3 * pc_y[k] * nsi_154[k];

        t_200[k] = f_21 * msi_160[k]
                   + f_4 * nsh0_125[k]
                   - f_5 * nsh1_125[k]
                   + f_3 * pc_x[k] * nsi_160[k];

        t_201[k] = f_21 * msi_161[k]
                   + f_3 * pc_x[k] * nsi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, msi_162, msi_163, \
                         msi_164, msi_165, nsi_160, nsi_162, nsi_163, nsi_164, \
                         nsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_21 * msi_162[k]
                   + f_3 * pc_x[k] * nsi_162[k];

        t_203[k] = f_21 * msi_163[k]
                   + f_3 * pc_x[k] * nsi_163[k];

        t_204[k] = f_21 * msi_164[k]
                   + f_3 * pc_x[k] * nsi_164[k];

        t_205[k] = f_21 * msi_165[k]
                   + f_3 * pc_x[k] * nsi_165[k];

        t_206[k] = f_3 * pc_y[k] * nsi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, msi_167, nsh0_120, nsh0_121, \
                         nsh1_120, nsh1_121, nsi_161, nsi_162, \
                         nsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_21 * msi_167[k]
                   + f_3 * pc_x[k] * nsi_167[k];

        t_208[k] = f_1 * nsh0_120[k]
                   - f_2 * nsh1_120[k]
                   + f_3 * pc_y[k] * nsi_161[k];

        t_209[k] = f_19 * nsh0_121[k]
                   - f_20 * nsh1_121[k]
                   + f_3 * pc_y[k] * nsi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, nsh0_122, nsh0_123, nsh0_124, nsh1_122, \
                         nsh1_123, nsh1_124, nsi_163, nsi_164, \
                         nsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * nsh0_122[k]
                   - f_11 * nsh1_122[k]
                   + f_3 * pc_y[k] * nsi_163[k];

        t_211[k] = f_8 * nsh0_123[k]
                   - f_9 * nsh1_123[k]
                   + f_3 * pc_y[k] * nsi_164[k];

        t_212[k] = f_6 * nsh0_124[k]
                   - f_7 * nsh1_124[k]
                   + f_3 * pc_y[k] * nsi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, msi_83, msi_168, \
                         nsh0_125, nsh0_126, nsh1_125, nsh1_126, nsi_166, nsi_167, \
                         nsi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * nsh0_125[k]
                   - f_5 * nsh1_125[k]
                   + f_3 * pc_y[k] * nsi_166[k];

        t_214[k] = f_3 * pc_y[k] * nsi_167[k];

        t_215[k] = f_14 * msi_83[k]
                   + f_1 * nsh0_125[k]
                   - f_2 * nsh1_125[k]
                   + f_3 * pc_z[k] * nsi_167[k];

        t_216[k] = f_22 * msi_168[k]
                   + f_1 * nsh0_126[k]
                   - f_2 * nsh1_126[k]
                   + f_3 * pc_x[k] * nsi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, msi_84, msi_171, \
                         nsh0_129, nsh1_129, nsi_168, nsi_169, \
                         nsi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * msi_84[k]
                   + f_3 * pc_y[k] * nsi_168[k];

        t_218[k] = f_3 * pc_z[k] * nsi_168[k];

        t_219[k] = f_22 * msi_171[k]
                   + f_10 * nsh0_129[k]
                   - f_11 * nsh1_129[k]
                   + f_3 * pc_x[k] * nsi_171[k];

        t_220[k] = f_3 * pc_z[k] * nsi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pc_x, pc_z, msi_174, nsh0_126, nsh0_132, \
                         nsh1_126, nsh1_132, nsi_170, nsi_171, \
                         nsi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * nsh0_126[k]
                   - f_5 * nsh1_126[k]
                   + f_3 * pc_z[k] * nsi_170[k];

        t_222[k] = f_22 * msi_174[k]
                   + f_8 * nsh0_132[k]
                   - f_9 * nsh1_132[k]
                   + f_3 * pc_x[k] * nsi_174[k];

        t_223[k] = f_3 * pc_z[k] * nsi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, msi_89, msi_178, \
                         nsh0_128, nsh0_136, nsh1_128, nsh1_136, nsi_173, nsi_174, \
                         nsi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * msi_89[k]
                   + f_3 * pc_y[k] * nsi_173[k];

        t_225[k] = f_6 * nsh0_128[k]
                   - f_7 * nsh1_128[k]
                   + f_3 * pc_z[k] * nsi_173[k];

        t_226[k] = f_22 * msi_178[k]
                   + f_6 * nsh0_136[k]
                   - f_7 * nsh1_136[k]
                   + f_3 * pc_x[k] * nsi_178[k];

        t_227[k] = f_3 * pc_z[k] * nsi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, msi_93, nsh0_129, nsh0_131, \
                         nsh1_129, nsh1_131, nsi_175, nsi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * nsh0_129[k]
                   - f_5 * nsh1_129[k]
                   + f_3 * pc_z[k] * nsi_175[k];

        t_229[k] = f_15 * msi_93[k]
                   + f_3 * pc_y[k] * nsi_177[k];

        t_230[k] = f_8 * nsh0_131[k]
                   - f_9 * nsh1_131[k]
                   + f_3 * pc_z[k] * nsi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_x, pc_z, msi_183, nsh0_132, nsh0_141, \
                         nsh1_132, nsh1_141, nsi_178, nsi_179, \
                         nsi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_22 * msi_183[k]
                   + f_4 * nsh0_141[k]
                   - f_5 * nsh1_141[k]
                   + f_3 * pc_x[k] * nsi_183[k];

        t_232[k] = f_3 * pc_z[k] * nsi_178[k];

        t_233[k] = f_4 * nsh0_132[k]
                   - f_5 * nsh1_132[k]
                   + f_3 * pc_z[k] * nsi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, msi_98, msi_189, \
                         nsh0_133, nsh0_135, nsh1_133, nsh1_135, nsi_180, nsi_182, \
                         nsi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * nsh0_133[k]
                   - f_7 * nsh1_133[k]
                   + f_3 * pc_z[k] * nsi_180[k];

        t_235[k] = f_15 * msi_98[k]
                   + f_3 * pc_y[k] * nsi_182[k];

        t_236[k] = f_10 * nsh0_135[k]
                   - f_11 * nsh1_135[k]
                   + f_3 * pc_z[k] * nsi_182[k];

        t_237[k] = f_22 * msi_189[k]
                   + f_3 * pc_x[k] * nsi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, msi_191, msi_192, \
                         msi_193, msi_194, nsi_183, nsi_191, nsi_192, nsi_193, \
                         nsi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * nsi_183[k];

        t_239[k] = f_22 * msi_191[k]
                   + f_3 * pc_x[k] * nsi_191[k];

        t_240[k] = f_22 * msi_192[k]
                   + f_3 * pc_x[k] * nsi_192[k];

        t_241[k] = f_22 * msi_193[k]
                   + f_3 * pc_x[k] * nsi_193[k];

        t_242[k] = f_22 * msi_194[k]
                   + f_3 * pc_x[k] * nsi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, msi_105, msi_195, \
                         nsh0_141, nsh1_141, nsi_189, nsi_190, \
                         nsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_22 * msi_195[k]
                   + f_3 * pc_x[k] * nsi_195[k];

        t_244[k] = f_15 * msi_105[k]
                   + f_1 * nsh0_141[k]
                   - f_2 * nsh1_141[k]
                   + f_3 * pc_y[k] * nsi_189[k];

        t_245[k] = f_3 * pc_z[k] * nsi_189[k];

        t_246[k] = f_4 * nsh0_141[k]
                   - f_5 * nsh1_141[k]
                   + f_3 * pc_z[k] * nsi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, nsh0_142, nsh0_143, nsh0_144, nsh1_142, \
                         nsh1_143, nsh1_144, nsi_191, nsi_192, \
                         nsi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * nsh0_142[k]
                   - f_7 * nsh1_142[k]
                   + f_3 * pc_z[k] * nsi_191[k];

        t_248[k] = f_8 * nsh0_143[k]
                   - f_9 * nsh1_143[k]
                   + f_3 * pc_z[k] * nsi_192[k];

        t_249[k] = f_10 * nsh0_144[k]
                   - f_11 * nsh1_144[k]
                   + f_3 * pc_z[k] * nsi_193[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_108 = buffer.data(msk0 + 108);
    const auto *msk0_111 = buffer.data(msk0 + 111);
    const auto *msk0_114 = buffer.data(msk0 + 114);
    const auto *msk0_118 = buffer.data(msk0 + 118);
    const auto *msk0_120 = buffer.data(msk0 + 120);
    const auto *msk0_123 = buffer.data(msk0 + 123);
    const auto *msk0_125 = buffer.data(msk0 + 125);
    const auto *msk0_126 = buffer.data(msk0 + 126);
    const auto *msk0_136 = buffer.data(msk0 + 136);
    const auto *msk0_180 = buffer.data(msk0 + 180);
    const auto *msk0_183 = buffer.data(msk0 + 183);
    const auto *msk0_185 = buffer.data(msk0 + 185);
    const auto *msk0_186 = buffer.data(msk0 + 186);
    const auto *msk0_189 = buffer.data(msk0 + 189);
    const auto *msk0_190 = buffer.data(msk0 + 190);
    const auto *msk0_192 = buffer.data(msk0 + 192);
    const auto *msk0_194 = buffer.data(msk0 + 194);
    const auto *msk0_195 = buffer.data(msk0 + 195);
    const auto *msk0_197 = buffer.data(msk0 + 197);
    const auto *msk0_198 = buffer.data(msk0 + 198);
    const auto *msk0_200 = buffer.data(msk0 + 200);
    const auto *msk0_215 = buffer.data(msk0 + 215);

    const auto *msi_84 = buffer.data(msi + 84);
    const auto *msi_87 = buffer.data(msi + 87);
    const auto *msi_90 = buffer.data(msi + 90);
    const auto *msi_91 = buffer.data(msi + 91);
    const auto *msi_94 = buffer.data(msi + 94);
    const auto *msi_95 = buffer.data(msi + 95);
    const auto *msi_96 = buffer.data(msi + 96);
    const auto *msi_105 = buffer.data(msi + 105);
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
    const auto *msi_135 = buffer.data(msi + 135);
    const auto *msi_136 = buffer.data(msi + 136);
    const auto *msi_137 = buffer.data(msi + 137);
    const auto *msi_138 = buffer.data(msi + 138);
    const auto *msi_139 = buffer.data(msi + 139);
    const auto *msi_140 = buffer.data(msi + 140);
    const auto *msi_141 = buffer.data(msi + 141);
    const auto *msi_142 = buffer.data(msi + 142);
    const auto *msi_143 = buffer.data(msi + 143);
    const auto *msi_145 = buffer.data(msi + 145);
    const auto *msi_146 = buffer.data(msi + 146);
    const auto *msi_148 = buffer.data(msi + 148);
    const auto *msi_149 = buffer.data(msi + 149);
    const auto *msi_150 = buffer.data(msi + 150);
    const auto *msi_152 = buffer.data(msi + 152);
    const auto *msi_153 = buffer.data(msi + 153);
    const auto *msi_154 = buffer.data(msi + 154);
    const auto *msi_161 = buffer.data(msi + 161);
    const auto *msi_163 = buffer.data(msi + 163);
    const auto *msi_164 = buffer.data(msi + 164);
    const auto *msi_165 = buffer.data(msi + 165);
    const auto *msi_166 = buffer.data(msi + 166);
    const auto *msi_167 = buffer.data(msi + 167);
    const auto *msi_168 = buffer.data(msi + 168);
    const auto *msi_201 = buffer.data(msi + 201);
    const auto *msi_205 = buffer.data(msi + 205);
    const auto *msi_210 = buffer.data(msi + 210);
    const auto *msi_216 = buffer.data(msi + 216);
    const auto *msi_217 = buffer.data(msi + 217);
    const auto *msi_218 = buffer.data(msi + 218);
    const auto *msi_219 = buffer.data(msi + 219);
    const auto *msi_220 = buffer.data(msi + 220);
    const auto *msi_221 = buffer.data(msi + 221);
    const auto *msi_222 = buffer.data(msi + 222);
    const auto *msi_223 = buffer.data(msi + 223);
    const auto *msi_245 = buffer.data(msi + 245);
    const auto *msi_246 = buffer.data(msi + 246);
    const auto *msi_247 = buffer.data(msi + 247);
    const auto *msi_248 = buffer.data(msi + 248);
    const auto *msi_249 = buffer.data(msi + 249);
    const auto *msi_250 = buffer.data(msi + 250);
    const auto *msi_251 = buffer.data(msi + 251);
    const auto *msi_252 = buffer.data(msi + 252);
    const auto *msi_257 = buffer.data(msi + 257);
    const auto *msi_261 = buffer.data(msi + 261);
    const auto *msi_266 = buffer.data(msi + 266);
    const auto *msi_272 = buffer.data(msi + 272);
    const auto *msi_273 = buffer.data(msi + 273);
    const auto *msi_274 = buffer.data(msi + 274);
    const auto *msi_275 = buffer.data(msi + 275);
    const auto *msi_276 = buffer.data(msi + 276);
    const auto *msi_277 = buffer.data(msi + 277);
    const auto *msi_279 = buffer.data(msi + 279);
    const auto *msi_280 = buffer.data(msi + 280);
    const auto *msi_283 = buffer.data(msi + 283);
    const auto *msi_286 = buffer.data(msi + 286);

    const auto *msk1_108 = buffer.data(msk1 + 108);
    const auto *msk1_111 = buffer.data(msk1 + 111);
    const auto *msk1_114 = buffer.data(msk1 + 114);
    const auto *msk1_118 = buffer.data(msk1 + 118);
    const auto *msk1_120 = buffer.data(msk1 + 120);
    const auto *msk1_123 = buffer.data(msk1 + 123);
    const auto *msk1_125 = buffer.data(msk1 + 125);
    const auto *msk1_126 = buffer.data(msk1 + 126);
    const auto *msk1_136 = buffer.data(msk1 + 136);
    const auto *msk1_180 = buffer.data(msk1 + 180);
    const auto *msk1_183 = buffer.data(msk1 + 183);
    const auto *msk1_185 = buffer.data(msk1 + 185);
    const auto *msk1_186 = buffer.data(msk1 + 186);
    const auto *msk1_189 = buffer.data(msk1 + 189);
    const auto *msk1_190 = buffer.data(msk1 + 190);
    const auto *msk1_192 = buffer.data(msk1 + 192);
    const auto *msk1_194 = buffer.data(msk1 + 194);
    const auto *msk1_195 = buffer.data(msk1 + 195);
    const auto *msk1_197 = buffer.data(msk1 + 197);
    const auto *msk1_198 = buffer.data(msk1 + 198);
    const auto *msk1_200 = buffer.data(msk1 + 200);
    const auto *msk1_215 = buffer.data(msk1 + 215);

    const auto *nsh0_146 = buffer.data(nsh0 + 146);
    const auto *nsh0_152 = buffer.data(nsh0 + 152);
    const auto *nsh0_156 = buffer.data(nsh0 + 156);
    const auto *nsh0_161 = buffer.data(nsh0 + 161);
    const auto *nsh0_164 = buffer.data(nsh0 + 164);
    const auto *nsh0_165 = buffer.data(nsh0 + 165);
    const auto *nsh0_166 = buffer.data(nsh0 + 166);
    const auto *nsh0_167 = buffer.data(nsh0 + 167);
    const auto *nsh0_183 = buffer.data(nsh0 + 183);
    const auto *nsh0_185 = buffer.data(nsh0 + 185);
    const auto *nsh0_186 = buffer.data(nsh0 + 186);
    const auto *nsh0_187 = buffer.data(nsh0 + 187);
    const auto *nsh0_188 = buffer.data(nsh0 + 188);
    const auto *nsh0_189 = buffer.data(nsh0 + 189);
    const auto *nsh0_190 = buffer.data(nsh0 + 190);
    const auto *nsh0_191 = buffer.data(nsh0 + 191);
    const auto *nsh0_192 = buffer.data(nsh0 + 192);
    const auto *nsh0_193 = buffer.data(nsh0 + 193);
    const auto *nsh0_194 = buffer.data(nsh0 + 194);
    const auto *nsh0_195 = buffer.data(nsh0 + 195);
    const auto *nsh0_196 = buffer.data(nsh0 + 196);
    const auto *nsh0_197 = buffer.data(nsh0 + 197);
    const auto *nsh0_198 = buffer.data(nsh0 + 198);
    const auto *nsh0_203 = buffer.data(nsh0 + 203);
    const auto *nsh0_204 = buffer.data(nsh0 + 204);
    const auto *nsh0_205 = buffer.data(nsh0 + 205);
    const auto *nsh0_206 = buffer.data(nsh0 + 206);
    const auto *nsh0_207 = buffer.data(nsh0 + 207);
    const auto *nsh0_208 = buffer.data(nsh0 + 208);
    const auto *nsh0_209 = buffer.data(nsh0 + 209);
    const auto *nsh0_210 = buffer.data(nsh0 + 210);
    const auto *nsh0_213 = buffer.data(nsh0 + 213);
    const auto *nsh0_216 = buffer.data(nsh0 + 216);

    const auto *nsh1_146 = buffer.data(nsh1 + 146);
    const auto *nsh1_152 = buffer.data(nsh1 + 152);
    const auto *nsh1_156 = buffer.data(nsh1 + 156);
    const auto *nsh1_161 = buffer.data(nsh1 + 161);
    const auto *nsh1_164 = buffer.data(nsh1 + 164);
    const auto *nsh1_165 = buffer.data(nsh1 + 165);
    const auto *nsh1_166 = buffer.data(nsh1 + 166);
    const auto *nsh1_167 = buffer.data(nsh1 + 167);
    const auto *nsh1_183 = buffer.data(nsh1 + 183);
    const auto *nsh1_185 = buffer.data(nsh1 + 185);
    const auto *nsh1_186 = buffer.data(nsh1 + 186);
    const auto *nsh1_187 = buffer.data(nsh1 + 187);
    const auto *nsh1_188 = buffer.data(nsh1 + 188);
    const auto *nsh1_189 = buffer.data(nsh1 + 189);
    const auto *nsh1_190 = buffer.data(nsh1 + 190);
    const auto *nsh1_191 = buffer.data(nsh1 + 191);
    const auto *nsh1_192 = buffer.data(nsh1 + 192);
    const auto *nsh1_193 = buffer.data(nsh1 + 193);
    const auto *nsh1_194 = buffer.data(nsh1 + 194);
    const auto *nsh1_195 = buffer.data(nsh1 + 195);
    const auto *nsh1_196 = buffer.data(nsh1 + 196);
    const auto *nsh1_197 = buffer.data(nsh1 + 197);
    const auto *nsh1_198 = buffer.data(nsh1 + 198);
    const auto *nsh1_203 = buffer.data(nsh1 + 203);
    const auto *nsh1_204 = buffer.data(nsh1 + 204);
    const auto *nsh1_205 = buffer.data(nsh1 + 205);
    const auto *nsh1_206 = buffer.data(nsh1 + 206);
    const auto *nsh1_207 = buffer.data(nsh1 + 207);
    const auto *nsh1_208 = buffer.data(nsh1 + 208);
    const auto *nsh1_209 = buffer.data(nsh1 + 209);
    const auto *nsh1_210 = buffer.data(nsh1 + 210);
    const auto *nsh1_213 = buffer.data(nsh1 + 213);
    const auto *nsh1_216 = buffer.data(nsh1 + 216);

    const auto *nsi_195 = buffer.data(nsi + 195);
    const auto *nsi_196 = buffer.data(nsi + 196);
    const auto *nsi_198 = buffer.data(nsi + 198);
    const auto *nsi_199 = buffer.data(nsi + 199);
    const auto *nsi_201 = buffer.data(nsi + 201);
    const auto *nsi_202 = buffer.data(nsi + 202);
    const auto *nsi_205 = buffer.data(nsi + 205);
    const auto *nsi_206 = buffer.data(nsi + 206);
    const auto *nsi_210 = buffer.data(nsi + 210);
    const auto *nsi_216 = buffer.data(nsi + 216);
    const auto *nsi_217 = buffer.data(nsi + 217);
    const auto *nsi_218 = buffer.data(nsi + 218);
    const auto *nsi_219 = buffer.data(nsi + 219);
    const auto *nsi_220 = buffer.data(nsi + 220);
    const auto *nsi_221 = buffer.data(nsi + 221);
    const auto *nsi_222 = buffer.data(nsi + 222);
    const auto *nsi_223 = buffer.data(nsi + 223);
    const auto *nsi_224 = buffer.data(nsi + 224);
    const auto *nsi_226 = buffer.data(nsi + 226);
    const auto *nsi_227 = buffer.data(nsi + 227);
    const auto *nsi_229 = buffer.data(nsi + 229);
    const auto *nsi_230 = buffer.data(nsi + 230);
    const auto *nsi_233 = buffer.data(nsi + 233);
    const auto *nsi_234 = buffer.data(nsi + 234);
    const auto *nsi_238 = buffer.data(nsi + 238);
    const auto *nsi_245 = buffer.data(nsi + 245);
    const auto *nsi_246 = buffer.data(nsi + 246);
    const auto *nsi_247 = buffer.data(nsi + 247);
    const auto *nsi_248 = buffer.data(nsi + 248);
    const auto *nsi_249 = buffer.data(nsi + 249);
    const auto *nsi_250 = buffer.data(nsi + 250);
    const auto *nsi_251 = buffer.data(nsi + 251);
    const auto *nsi_252 = buffer.data(nsi + 252);
    const auto *nsi_253 = buffer.data(nsi + 253);
    const auto *nsi_254 = buffer.data(nsi + 254);
    const auto *nsi_255 = buffer.data(nsi + 255);
    const auto *nsi_256 = buffer.data(nsi + 256);
    const auto *nsi_257 = buffer.data(nsi + 257);
    const auto *nsi_258 = buffer.data(nsi + 258);
    const auto *nsi_259 = buffer.data(nsi + 259);
    const auto *nsi_260 = buffer.data(nsi + 260);
    const auto *nsi_261 = buffer.data(nsi + 261);
    const auto *nsi_262 = buffer.data(nsi + 262);
    const auto *nsi_263 = buffer.data(nsi + 263);
    const auto *nsi_264 = buffer.data(nsi + 264);
    const auto *nsi_265 = buffer.data(nsi + 265);
    const auto *nsi_266 = buffer.data(nsi + 266);
    const auto *nsi_272 = buffer.data(nsi + 272);
    const auto *nsi_273 = buffer.data(nsi + 273);
    const auto *nsi_274 = buffer.data(nsi + 274);
    const auto *nsi_275 = buffer.data(nsi + 275);
    const auto *nsi_276 = buffer.data(nsi + 276);
    const auto *nsi_277 = buffer.data(nsi + 277);
    const auto *nsi_278 = buffer.data(nsi + 278);
    const auto *nsi_279 = buffer.data(nsi + 279);
    const auto *nsi_280 = buffer.data(nsi + 280);
    const auto *nsi_281 = buffer.data(nsi + 281);
    const auto *nsi_282 = buffer.data(nsi + 282);
    const auto *nsi_283 = buffer.data(nsi + 283);
    const auto *nsi_286 = buffer.data(nsi + 286);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, msk0_108, msi_111, \
                         msi_112, msk1_108, nsh0_146, nsh1_146, nsi_195, \
                         nsi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * msi_111[k]
                   + f_3 * pc_y[k] * nsi_195[k];

        t_251[k] = f_1 * nsh0_146[k]
                   - f_2 * nsh1_146[k]
                   + f_3 * pc_z[k] * nsi_195[k];

        t_252[k] = pa_z[k] * msk0_108[k]
                   - f_12 * pc_z[k] * msk1_108[k];

        t_253[k] = f_14 * msi_112[k]
                   + f_3 * pc_y[k] * nsi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_y, pc_z, msk0_111, msi_84, msi_114, \
                         msk1_111, nsi_196, nsi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * msi_84[k]
                   + f_3 * pc_z[k] * nsi_196[k];

        t_255[k] = pa_z[k] * msk0_111[k]
                   - f_12 * pc_z[k] * msk1_111[k];

        t_256[k] = f_14 * msi_114[k]
                   + f_3 * pc_y[k] * nsi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, msk0_114, msi_87, msi_201, \
                         msk1_114, nsh0_152, nsh1_152, nsi_199, \
                         nsi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_22 * msi_201[k]
                   + f_10 * nsh0_152[k]
                   - f_11 * nsh1_152[k]
                   + f_3 * pc_x[k] * nsi_201[k];

        t_258[k] = pa_z[k] * msk0_114[k]
                   - f_12 * pc_z[k] * msk1_114[k];

        t_259[k] = f_13 * msi_87[k]
                   + f_3 * pc_z[k] * nsi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_y, pc_z, msk0_118, msi_117, \
                         msi_205, msk1_118, nsh0_156, nsh1_156, nsi_201, \
                         nsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * msi_117[k]
                   + f_3 * pc_y[k] * nsi_201[k];

        t_261[k] = f_22 * msi_205[k]
                   + f_8 * nsh0_156[k]
                   - f_9 * nsh1_156[k]
                   + f_3 * pc_x[k] * nsi_205[k];

        t_262[k] = pa_z[k] * msk0_118[k]
                   - f_12 * pc_z[k] * msk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pa_z, pc_y, pc_z, msk0_120, msi_90, msi_91, \
                         msi_121, msk1_120, nsi_202, nsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * msi_90[k]
                   + f_3 * pc_z[k] * nsi_202[k];

        t_264[k] = pa_z[k] * msk0_120[k]
                   + f_14 * msi_91[k]
                   - f_12 * pc_z[k] * msk1_120[k];

        t_265[k] = f_14 * msi_121[k]
                   + f_3 * pc_y[k] * nsi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, msk0_123, msi_94, msi_210, \
                         msk1_123, nsh0_161, nsh1_161, nsi_206, \
                         nsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_22 * msi_210[k]
                   + f_6 * nsh0_161[k]
                   - f_7 * nsh1_161[k]
                   + f_3 * pc_x[k] * nsi_210[k];

        t_267[k] = pa_z[k] * msk0_123[k]
                   - f_12 * pc_z[k] * msk1_123[k];

        t_268[k] = f_13 * msi_94[k]
                   + f_3 * pc_z[k] * nsi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_z, pc_y, pc_z, msk0_125, msk0_126, msi_95, \
                         msi_96, msi_126, msk1_125, msk1_126, nsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_z[k] * msk0_125[k]
                   + f_14 * msi_95[k]
                   - f_12 * pc_z[k] * msk1_125[k];

        t_270[k] = pa_z[k] * msk0_126[k]
                   + f_15 * msi_96[k]
                   - f_12 * pc_z[k] * msk1_126[k];

        t_271[k] = f_14 * msi_126[k]
                   + f_3 * pc_y[k] * nsi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, msi_216, msi_217, msi_218, msi_219, \
                         nsh0_167, nsh1_167, nsi_216, nsi_217, nsi_218, \
                         nsi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_22 * msi_216[k]
                   + f_4 * nsh0_167[k]
                   - f_5 * nsh1_167[k]
                   + f_3 * pc_x[k] * nsi_216[k];

        t_273[k] = f_22 * msi_217[k]
                   + f_3 * pc_x[k] * nsi_217[k];

        t_274[k] = f_22 * msi_218[k]
                   + f_3 * pc_x[k] * nsi_218[k];

        t_275[k] = f_22 * msi_219[k]
                   + f_3 * pc_x[k] * nsi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, msi_220, msi_221, msi_222, msi_223, \
                         nsi_220, nsi_221, nsi_222, nsi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_22 * msi_220[k]
                   + f_3 * pc_x[k] * nsi_220[k];

        t_277[k] = f_22 * msi_221[k]
                   + f_3 * pc_x[k] * nsi_221[k];

        t_278[k] = f_22 * msi_222[k]
                   + f_3 * pc_x[k] * nsi_222[k];

        t_279[k] = f_22 * msi_223[k]
                   + f_3 * pc_x[k] * nsi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_z, pc_y, pc_z, msk0_136, msi_105, msi_135, \
                         msk1_136, nsh0_164, nsh1_164, nsi_217, \
                         nsi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * msk0_136[k]
                   - f_12 * pc_z[k] * msk1_136[k];

        t_281[k] = f_13 * msi_105[k]
                   + f_3 * pc_z[k] * nsi_217[k];

        t_282[k] = f_14 * msi_135[k]
                   + f_10 * nsh0_164[k]
                   - f_11 * nsh1_164[k]
                   + f_3 * pc_y[k] * nsi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, msi_136, msi_137, msi_138, nsh0_165, \
                         nsh0_166, nsh0_167, nsh1_165, nsh1_166, nsh1_167, nsi_220, nsi_221, \
                         nsi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * msi_136[k]
                   + f_8 * nsh0_165[k]
                   - f_9 * nsh1_165[k]
                   + f_3 * pc_y[k] * nsi_220[k];

        t_284[k] = f_14 * msi_137[k]
                   + f_6 * nsh0_166[k]
                   - f_7 * nsh1_166[k]
                   + f_3 * pc_y[k] * nsi_221[k];

        t_285[k] = f_14 * msi_138[k]
                   + f_4 * nsh0_167[k]
                   - f_5 * nsh1_167[k]
                   + f_3 * pc_y[k] * nsi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, msk0_180, msi_111, \
                         msi_139, msi_140, msk1_180, nsh0_167, nsh1_167, nsi_223, \
                         nsi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * msi_139[k]
                   + f_3 * pc_y[k] * nsi_223[k];

        t_287[k] = f_13 * msi_111[k]
                   + f_1 * nsh0_167[k]
                   - f_2 * nsh1_167[k]
                   + f_3 * pc_z[k] * nsi_223[k];

        t_288[k] = pa_y[k] * msk0_180[k]
                   - f_12 * pc_y[k] * msk1_180[k];

        t_289[k] = f_13 * msi_140[k]
                   + f_3 * pc_y[k] * nsi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, pc_z, msk0_183, msk0_185, \
                         msi_112, msi_141, msi_142, msk1_183, msk1_185, nsi_224, \
                         nsi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * msi_112[k]
                   + f_3 * pc_z[k] * nsi_224[k];

        t_291[k] = pa_y[k] * msk0_183[k]
                   + f_14 * msi_141[k]
                   - f_12 * pc_y[k] * msk1_183[k];

        t_292[k] = f_13 * msi_142[k]
                   + f_3 * pc_y[k] * nsi_226[k];

        t_293[k] = pa_y[k] * msk0_185[k]
                   - f_12 * pc_y[k] * msk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pc_y, pc_z, msk0_186, msk0_189, \
                         msi_115, msi_143, msi_145, msk1_186, msk1_189, nsi_227, \
                         nsi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * msk0_186[k]
                   + f_15 * msi_143[k]
                   - f_12 * pc_y[k] * msk1_186[k];

        t_295[k] = f_14 * msi_115[k]
                   + f_3 * pc_z[k] * nsi_227[k];

        t_296[k] = f_13 * msi_145[k]
                   + f_3 * pc_y[k] * nsi_229[k];

        t_297[k] = pa_y[k] * msk0_189[k]
                   - f_12 * pc_y[k] * msk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_y, pc_y, pc_z, msk0_190, msk0_192, msi_118, \
                         msi_146, msi_148, msk1_190, msk1_192, \
                         nsi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_y[k] * msk0_190[k]
                   + f_16 * msi_146[k]
                   - f_12 * pc_y[k] * msk1_190[k];

        t_299[k] = f_14 * msi_118[k]
                   + f_3 * pc_z[k] * nsi_230[k];

        t_300[k] = pa_y[k] * msk0_192[k]
                   + f_14 * msi_148[k]
                   - f_12 * pc_y[k] * msk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, pc_y, pc_z, msk0_194, msk0_195, \
                         msi_122, msi_149, msi_150, msk1_194, msk1_195, nsi_233, \
                         nsi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * msi_149[k]
                   + f_3 * pc_y[k] * nsi_233[k];

        t_302[k] = pa_y[k] * msk0_194[k]
                   - f_12 * pc_y[k] * msk1_194[k];

        t_303[k] = pa_y[k] * msk0_195[k]
                   + f_17 * msi_150[k]
                   - f_12 * pc_y[k] * msk1_195[k];

        t_304[k] = f_14 * msi_122[k]
                   + f_3 * pc_z[k] * nsi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pc_y, msk0_197, msk0_198, msk0_200, \
                         msi_152, msi_153, msi_154, msk1_197, msk1_198, msk1_200, \
                         nsi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_y[k] * msk0_197[k]
                   + f_15 * msi_152[k]
                   - f_12 * pc_y[k] * msk1_197[k];

        t_306[k] = pa_y[k] * msk0_198[k]
                   + f_14 * msi_153[k]
                   - f_12 * pc_y[k] * msk1_198[k];

        t_307[k] = f_13 * msi_154[k]
                   + f_3 * pc_y[k] * nsi_238[k];

        t_308[k] = pa_y[k] * msk0_200[k]
                   - f_12 * pc_y[k] * msk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, msi_245, msi_246, msi_247, \
                         msi_248, msi_249, nsi_245, nsi_246, nsi_247, nsi_248, \
                         nsi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_22 * msi_245[k]
                   + f_3 * pc_x[k] * nsi_245[k];

        t_310[k] = f_22 * msi_246[k]
                   + f_3 * pc_x[k] * nsi_246[k];

        t_311[k] = f_22 * msi_247[k]
                   + f_3 * pc_x[k] * nsi_247[k];

        t_312[k] = f_22 * msi_248[k]
                   + f_3 * pc_x[k] * nsi_248[k];

        t_313[k] = f_22 * msi_249[k]
                   + f_3 * pc_x[k] * nsi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, msi_133, msi_161, \
                         msi_250, msi_251, nsh0_183, nsh1_183, nsi_245, nsi_250, \
                         nsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_22 * msi_250[k]
                   + f_3 * pc_x[k] * nsi_250[k];

        t_315[k] = f_22 * msi_251[k]
                   + f_3 * pc_x[k] * nsi_251[k];

        t_316[k] = f_13 * msi_161[k]
                   + f_1 * nsh0_183[k]
                   - f_2 * nsh1_183[k]
                   + f_3 * pc_y[k] * nsi_245[k];

        t_317[k] = f_14 * msi_133[k]
                   + f_3 * pc_z[k] * nsi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, msi_163, msi_164, msi_165, nsh0_185, \
                         nsh0_186, nsh0_187, nsh1_185, nsh1_186, nsh1_187, nsi_247, nsi_248, \
                         nsi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * msi_163[k]
                   + f_10 * nsh0_185[k]
                   - f_11 * nsh1_185[k]
                   + f_3 * pc_y[k] * nsi_247[k];

        t_319[k] = f_13 * msi_164[k]
                   + f_8 * nsh0_186[k]
                   - f_9 * nsh1_186[k]
                   + f_3 * pc_y[k] * nsi_248[k];

        t_320[k] = f_13 * msi_165[k]
                   + f_6 * nsh0_187[k]
                   - f_7 * nsh1_187[k]
                   + f_3 * pc_y[k] * nsi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, msk0_215, msi_166, msi_167, \
                         msk1_215, nsh0_188, nsh1_188, nsi_250, \
                         nsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * msi_166[k]
                   + f_4 * nsh0_188[k]
                   - f_5 * nsh1_188[k]
                   + f_3 * pc_y[k] * nsi_250[k];

        t_322[k] = f_13 * msi_167[k]
                   + f_3 * pc_y[k] * nsi_251[k];

        t_323[k] = pa_y[k] * msk0_215[k]
                   - f_12 * pc_y[k] * msk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, msi_140, \
                         msi_252, nsh0_189, nsh1_189, nsi_252, nsi_253, \
                         nsi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_22 * msi_252[k]
                   + f_1 * nsh0_189[k]
                   - f_2 * nsh1_189[k]
                   + f_3 * pc_x[k] * nsi_252[k];

        t_325[k] = f_3 * pc_y[k] * nsi_252[k];

        t_326[k] = f_15 * msi_140[k]
                   + f_3 * pc_z[k] * nsi_252[k];

        t_327[k] = f_4 * nsh0_189[k]
                   - f_5 * nsh1_189[k]
                   + f_3 * pc_y[k] * nsi_253[k];

        t_328[k] = f_3 * pc_y[k] * nsi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, msi_257, nsh0_190, nsh0_191, \
                         nsh0_194, nsh1_190, nsh1_191, nsh1_194, nsi_255, nsi_256, \
                         nsi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_22 * msi_257[k]
                   + f_10 * nsh0_194[k]
                   - f_11 * nsh1_194[k]
                   + f_3 * pc_x[k] * nsi_257[k];

        t_330[k] = f_6 * nsh0_190[k]
                   - f_7 * nsh1_190[k]
                   + f_3 * pc_y[k] * nsi_255[k];

        t_331[k] = f_4 * nsh0_191[k]
                   - f_5 * nsh1_191[k]
                   + f_3 * pc_y[k] * nsi_256[k];

        t_332[k] = f_3 * pc_y[k] * nsi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, msi_261, nsh0_192, nsh0_193, \
                         nsh0_198, nsh1_192, nsh1_193, nsh1_198, nsi_258, nsi_259, \
                         nsi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_22 * msi_261[k]
                   + f_8 * nsh0_198[k]
                   - f_9 * nsh1_198[k]
                   + f_3 * pc_x[k] * nsi_261[k];

        t_334[k] = f_8 * nsh0_192[k]
                   - f_9 * nsh1_192[k]
                   + f_3 * pc_y[k] * nsi_258[k];

        t_335[k] = f_6 * nsh0_193[k]
                   - f_7 * nsh1_193[k]
                   + f_3 * pc_y[k] * nsi_259[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, msi_266, nsh0_194, nsh0_203, \
                         nsh1_194, nsh1_203, nsi_260, nsi_261, \
                         nsi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_4 * nsh0_194[k]
                   - f_5 * nsh1_194[k]
                   + f_3 * pc_y[k] * nsi_260[k];

        t_337[k] = f_3 * pc_y[k] * nsi_261[k];

        t_338[k] = f_22 * msi_266[k]
                   + f_6 * nsh0_203[k]
                   - f_7 * nsh1_203[k]
                   + f_3 * pc_x[k] * nsi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_y, nsh0_195, nsh0_196, nsh0_197, nsh1_195, \
                         nsh1_196, nsh1_197, nsi_262, nsi_263, \
                         nsi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_10 * nsh0_195[k]
                   - f_11 * nsh1_195[k]
                   + f_3 * pc_y[k] * nsi_262[k];

        t_340[k] = f_8 * nsh0_196[k]
                   - f_9 * nsh1_196[k]
                   + f_3 * pc_y[k] * nsi_263[k];

        t_341[k] = f_6 * nsh0_197[k]
                   - f_7 * nsh1_197[k]
                   + f_3 * pc_y[k] * nsi_264[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_y, msi_272, msi_273, nsh0_198, \
                         nsh0_209, nsh1_198, nsh1_209, nsi_265, nsi_266, nsi_272, \
                         nsi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_4 * nsh0_198[k]
                   - f_5 * nsh1_198[k]
                   + f_3 * pc_y[k] * nsi_265[k];

        t_343[k] = f_3 * pc_y[k] * nsi_266[k];

        t_344[k] = f_22 * msi_272[k]
                   + f_4 * nsh0_209[k]
                   - f_5 * nsh1_209[k]
                   + f_3 * pc_x[k] * nsi_272[k];

        t_345[k] = f_22 * msi_273[k]
                   + f_3 * pc_x[k] * nsi_273[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pc_x, pc_y, msi_274, msi_275, \
                         msi_276, msi_277, nsi_272, nsi_274, nsi_275, nsi_276, \
                         nsi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_22 * msi_274[k]
                   + f_3 * pc_x[k] * nsi_274[k];

        t_347[k] = f_22 * msi_275[k]
                   + f_3 * pc_x[k] * nsi_275[k];

        t_348[k] = f_22 * msi_276[k]
                   + f_3 * pc_x[k] * nsi_276[k];

        t_349[k] = f_22 * msi_277[k]
                   + f_3 * pc_x[k] * nsi_277[k];

        t_350[k] = f_3 * pc_y[k] * nsi_272[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_x, pc_y, msi_279, nsh0_204, nsh0_205, \
                         nsh1_204, nsh1_205, nsi_273, nsi_274, \
                         nsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_22 * msi_279[k]
                   + f_3 * pc_x[k] * nsi_279[k];

        t_352[k] = f_1 * nsh0_204[k]
                   - f_2 * nsh1_204[k]
                   + f_3 * pc_y[k] * nsi_273[k];

        t_353[k] = f_19 * nsh0_205[k]
                   - f_20 * nsh1_205[k]
                   + f_3 * pc_y[k] * nsi_274[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_y, nsh0_206, nsh0_207, nsh0_208, nsh1_206, \
                         nsh1_207, nsh1_208, nsi_275, nsi_276, \
                         nsi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_10 * nsh0_206[k]
                   - f_11 * nsh1_206[k]
                   + f_3 * pc_y[k] * nsi_275[k];

        t_355[k] = f_8 * nsh0_207[k]
                   - f_9 * nsh1_207[k]
                   + f_3 * pc_y[k] * nsi_276[k];

        t_356[k] = f_6 * nsh0_208[k]
                   - f_7 * nsh1_208[k]
                   + f_3 * pc_y[k] * nsi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pc_x, pc_y, pc_z, msi_167, msi_280, \
                         nsh0_209, nsh0_210, nsh1_209, nsh1_210, nsi_278, nsi_279, \
                         nsi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_4 * nsh0_209[k]
                   - f_5 * nsh1_209[k]
                   + f_3 * pc_y[k] * nsi_278[k];

        t_358[k] = f_3 * pc_y[k] * nsi_279[k];

        t_359[k] = f_15 * msi_167[k]
                   + f_1 * nsh0_209[k]
                   - f_2 * nsh1_209[k]
                   + f_3 * pc_z[k] * nsi_279[k];

        t_360[k] = f_23 * msi_280[k]
                   + f_1 * nsh0_210[k]
                   - f_2 * nsh1_210[k]
                   + f_3 * pc_x[k] * nsi_280[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, msi_168, msi_283, \
                         nsh0_213, nsh1_213, nsi_280, nsi_281, \
                         nsi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_16 * msi_168[k]
                   + f_3 * pc_y[k] * nsi_280[k];

        t_362[k] = f_3 * pc_z[k] * nsi_280[k];

        t_363[k] = f_23 * msi_283[k]
                   + f_10 * nsh0_213[k]
                   - f_11 * nsh1_213[k]
                   + f_3 * pc_x[k] * nsi_283[k];

        t_364[k] = f_3 * pc_z[k] * nsi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pc_x, pc_z, msi_286, nsh0_210, nsh0_216, \
                         nsh1_210, nsh1_216, nsi_282, nsi_283, \
                         nsi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * nsh0_210[k]
                   - f_5 * nsh1_210[k]
                   + f_3 * pc_z[k] * nsi_282[k];

        t_366[k] = f_23 * msi_286[k]
                   + f_8 * nsh0_216[k]
                   - f_9 * nsh1_216[k]
                   + f_3 * pc_x[k] * nsi_286[k];

        t_367[k] = f_3 * pc_z[k] * nsi_283[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_216 = buffer.data(msk0 + 216);
    const auto *msk0_219 = buffer.data(msk0 + 219);
    const auto *msk0_222 = buffer.data(msk0 + 222);
    const auto *msk0_226 = buffer.data(msk0 + 226);
    const auto *msk0_228 = buffer.data(msk0 + 228);
    const auto *msk0_231 = buffer.data(msk0 + 231);
    const auto *msk0_233 = buffer.data(msk0 + 233);
    const auto *msk0_234 = buffer.data(msk0 + 234);
    const auto *msk0_244 = buffer.data(msk0 + 244);
    const auto *msk0_324 = buffer.data(msk0 + 324);
    const auto *msk0_327 = buffer.data(msk0 + 327);
    const auto *msk0_329 = buffer.data(msk0 + 329);
    const auto *msk0_330 = buffer.data(msk0 + 330);
    const auto *msk0_333 = buffer.data(msk0 + 333);
    const auto *msk0_334 = buffer.data(msk0 + 334);
    const auto *msk0_336 = buffer.data(msk0 + 336);

    const auto *msi_168 = buffer.data(msi + 168);
    const auto *msi_171 = buffer.data(msi + 171);
    const auto *msi_173 = buffer.data(msi + 173);
    const auto *msi_174 = buffer.data(msi + 174);
    const auto *msi_175 = buffer.data(msi + 175);
    const auto *msi_177 = buffer.data(msi + 177);
    const auto *msi_178 = buffer.data(msi + 178);
    const auto *msi_179 = buffer.data(msi + 179);
    const auto *msi_180 = buffer.data(msi + 180);
    const auto *msi_182 = buffer.data(msi + 182);
    const auto *msi_189 = buffer.data(msi + 189);
    const auto *msi_195 = buffer.data(msi + 195);
    const auto *msi_196 = buffer.data(msi + 196);
    const auto *msi_198 = buffer.data(msi + 198);
    const auto *msi_199 = buffer.data(msi + 199);
    const auto *msi_201 = buffer.data(msi + 201);
    const auto *msi_202 = buffer.data(msi + 202);
    const auto *msi_205 = buffer.data(msi + 205);
    const auto *msi_206 = buffer.data(msi + 206);
    const auto *msi_210 = buffer.data(msi + 210);
    const auto *msi_217 = buffer.data(msi + 217);
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
    const auto *msi_238 = buffer.data(msi + 238);
    const auto *msi_245 = buffer.data(msi + 245);
    const auto *msi_247 = buffer.data(msi + 247);
    const auto *msi_248 = buffer.data(msi + 248);
    const auto *msi_249 = buffer.data(msi + 249);
    const auto *msi_250 = buffer.data(msi + 250);
    const auto *msi_251 = buffer.data(msi + 251);
    const auto *msi_252 = buffer.data(msi + 252);
    const auto *msi_253 = buffer.data(msi + 253);
    const auto *msi_254 = buffer.data(msi + 254);
    const auto *msi_255 = buffer.data(msi + 255);
    const auto *msi_257 = buffer.data(msi + 257);
    const auto *msi_258 = buffer.data(msi + 258);
    const auto *msi_260 = buffer.data(msi + 260);
    const auto *msi_290 = buffer.data(msi + 290);
    const auto *msi_295 = buffer.data(msi + 295);
    const auto *msi_301 = buffer.data(msi + 301);
    const auto *msi_303 = buffer.data(msi + 303);
    const auto *msi_304 = buffer.data(msi + 304);
    const auto *msi_305 = buffer.data(msi + 305);
    const auto *msi_306 = buffer.data(msi + 306);
    const auto *msi_307 = buffer.data(msi + 307);
    const auto *msi_313 = buffer.data(msi + 313);
    const auto *msi_317 = buffer.data(msi + 317);
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

    const auto *msk1_216 = buffer.data(msk1 + 216);
    const auto *msk1_219 = buffer.data(msk1 + 219);
    const auto *msk1_222 = buffer.data(msk1 + 222);
    const auto *msk1_226 = buffer.data(msk1 + 226);
    const auto *msk1_228 = buffer.data(msk1 + 228);
    const auto *msk1_231 = buffer.data(msk1 + 231);
    const auto *msk1_233 = buffer.data(msk1 + 233);
    const auto *msk1_234 = buffer.data(msk1 + 234);
    const auto *msk1_244 = buffer.data(msk1 + 244);
    const auto *msk1_324 = buffer.data(msk1 + 324);
    const auto *msk1_327 = buffer.data(msk1 + 327);
    const auto *msk1_329 = buffer.data(msk1 + 329);
    const auto *msk1_330 = buffer.data(msk1 + 330);
    const auto *msk1_333 = buffer.data(msk1 + 333);
    const auto *msk1_334 = buffer.data(msk1 + 334);
    const auto *msk1_336 = buffer.data(msk1 + 336);

    const auto *nsh0_212 = buffer.data(nsh0 + 212);
    const auto *nsh0_213 = buffer.data(nsh0 + 213);
    const auto *nsh0_215 = buffer.data(nsh0 + 215);
    const auto *nsh0_216 = buffer.data(nsh0 + 216);
    const auto *nsh0_217 = buffer.data(nsh0 + 217);
    const auto *nsh0_219 = buffer.data(nsh0 + 219);
    const auto *nsh0_220 = buffer.data(nsh0 + 220);
    const auto *nsh0_225 = buffer.data(nsh0 + 225);
    const auto *nsh0_226 = buffer.data(nsh0 + 226);
    const auto *nsh0_227 = buffer.data(nsh0 + 227);
    const auto *nsh0_228 = buffer.data(nsh0 + 228);
    const auto *nsh0_230 = buffer.data(nsh0 + 230);
    const auto *nsh0_236 = buffer.data(nsh0 + 236);
    const auto *nsh0_240 = buffer.data(nsh0 + 240);
    const auto *nsh0_245 = buffer.data(nsh0 + 245);
    const auto *nsh0_248 = buffer.data(nsh0 + 248);
    const auto *nsh0_249 = buffer.data(nsh0 + 249);
    const auto *nsh0_250 = buffer.data(nsh0 + 250);
    const auto *nsh0_251 = buffer.data(nsh0 + 251);
    const auto *nsh0_252 = buffer.data(nsh0 + 252);
    const auto *nsh0_255 = buffer.data(nsh0 + 255);
    const auto *nsh0_257 = buffer.data(nsh0 + 257);
    const auto *nsh0_258 = buffer.data(nsh0 + 258);
    const auto *nsh0_261 = buffer.data(nsh0 + 261);
    const auto *nsh0_262 = buffer.data(nsh0 + 262);
    const auto *nsh0_264 = buffer.data(nsh0 + 264);
    const auto *nsh0_266 = buffer.data(nsh0 + 266);
    const auto *nsh0_267 = buffer.data(nsh0 + 267);
    const auto *nsh0_269 = buffer.data(nsh0 + 269);
    const auto *nsh0_270 = buffer.data(nsh0 + 270);
    const auto *nsh0_271 = buffer.data(nsh0 + 271);
    const auto *nsh0_272 = buffer.data(nsh0 + 272);

    const auto *nsh1_212 = buffer.data(nsh1 + 212);
    const auto *nsh1_213 = buffer.data(nsh1 + 213);
    const auto *nsh1_215 = buffer.data(nsh1 + 215);
    const auto *nsh1_216 = buffer.data(nsh1 + 216);
    const auto *nsh1_217 = buffer.data(nsh1 + 217);
    const auto *nsh1_219 = buffer.data(nsh1 + 219);
    const auto *nsh1_220 = buffer.data(nsh1 + 220);
    const auto *nsh1_225 = buffer.data(nsh1 + 225);
    const auto *nsh1_226 = buffer.data(nsh1 + 226);
    const auto *nsh1_227 = buffer.data(nsh1 + 227);
    const auto *nsh1_228 = buffer.data(nsh1 + 228);
    const auto *nsh1_230 = buffer.data(nsh1 + 230);
    const auto *nsh1_236 = buffer.data(nsh1 + 236);
    const auto *nsh1_240 = buffer.data(nsh1 + 240);
    const auto *nsh1_245 = buffer.data(nsh1 + 245);
    const auto *nsh1_248 = buffer.data(nsh1 + 248);
    const auto *nsh1_249 = buffer.data(nsh1 + 249);
    const auto *nsh1_250 = buffer.data(nsh1 + 250);
    const auto *nsh1_251 = buffer.data(nsh1 + 251);
    const auto *nsh1_252 = buffer.data(nsh1 + 252);
    const auto *nsh1_255 = buffer.data(nsh1 + 255);
    const auto *nsh1_257 = buffer.data(nsh1 + 257);
    const auto *nsh1_258 = buffer.data(nsh1 + 258);
    const auto *nsh1_261 = buffer.data(nsh1 + 261);
    const auto *nsh1_262 = buffer.data(nsh1 + 262);
    const auto *nsh1_264 = buffer.data(nsh1 + 264);
    const auto *nsh1_266 = buffer.data(nsh1 + 266);
    const auto *nsh1_267 = buffer.data(nsh1 + 267);
    const auto *nsh1_269 = buffer.data(nsh1 + 269);
    const auto *nsh1_270 = buffer.data(nsh1 + 270);
    const auto *nsh1_271 = buffer.data(nsh1 + 271);
    const auto *nsh1_272 = buffer.data(nsh1 + 272);

    const auto *nsi_285 = buffer.data(nsi + 285);
    const auto *nsi_286 = buffer.data(nsi + 286);
    const auto *nsi_287 = buffer.data(nsi + 287);
    const auto *nsi_289 = buffer.data(nsi + 289);
    const auto *nsi_290 = buffer.data(nsi + 290);
    const auto *nsi_291 = buffer.data(nsi + 291);
    const auto *nsi_292 = buffer.data(nsi + 292);
    const auto *nsi_294 = buffer.data(nsi + 294);
    const auto *nsi_295 = buffer.data(nsi + 295);
    const auto *nsi_301 = buffer.data(nsi + 301);
    const auto *nsi_302 = buffer.data(nsi + 302);
    const auto *nsi_303 = buffer.data(nsi + 303);
    const auto *nsi_304 = buffer.data(nsi + 304);
    const auto *nsi_305 = buffer.data(nsi + 305);
    const auto *nsi_306 = buffer.data(nsi + 306);
    const auto *nsi_307 = buffer.data(nsi + 307);
    const auto *nsi_308 = buffer.data(nsi + 308);
    const auto *nsi_310 = buffer.data(nsi + 310);
    const auto *nsi_311 = buffer.data(nsi + 311);
    const auto *nsi_313 = buffer.data(nsi + 313);
    const auto *nsi_314 = buffer.data(nsi + 314);
    const auto *nsi_317 = buffer.data(nsi + 317);
    const auto *nsi_318 = buffer.data(nsi + 318);
    const auto *nsi_322 = buffer.data(nsi + 322);
    const auto *nsi_328 = buffer.data(nsi + 328);
    const auto *nsi_329 = buffer.data(nsi + 329);
    const auto *nsi_330 = buffer.data(nsi + 330);
    const auto *nsi_331 = buffer.data(nsi + 331);
    const auto *nsi_332 = buffer.data(nsi + 332);
    const auto *nsi_333 = buffer.data(nsi + 333);
    const auto *nsi_334 = buffer.data(nsi + 334);
    const auto *nsi_335 = buffer.data(nsi + 335);
    const auto *nsi_336 = buffer.data(nsi + 336);
    const auto *nsi_338 = buffer.data(nsi + 338);
    const auto *nsi_339 = buffer.data(nsi + 339);
    const auto *nsi_341 = buffer.data(nsi + 341);
    const auto *nsi_342 = buffer.data(nsi + 342);
    const auto *nsi_345 = buffer.data(nsi + 345);
    const auto *nsi_346 = buffer.data(nsi + 346);
    const auto *nsi_348 = buffer.data(nsi + 348);
    const auto *nsi_350 = buffer.data(nsi + 350);
    const auto *nsi_351 = buffer.data(nsi + 351);
    const auto *nsi_353 = buffer.data(nsi + 353);
    const auto *nsi_354 = buffer.data(nsi + 354);
    const auto *nsi_356 = buffer.data(nsi + 356);
    const auto *nsi_357 = buffer.data(nsi + 357);
    const auto *nsi_358 = buffer.data(nsi + 358);
    const auto *nsi_359 = buffer.data(nsi + 359);
    const auto *nsi_360 = buffer.data(nsi + 360);
    const auto *nsi_361 = buffer.data(nsi + 361);
    const auto *nsi_362 = buffer.data(nsi + 362);
    const auto *nsi_363 = buffer.data(nsi + 363);
    const auto *nsi_364 = buffer.data(nsi + 364);
    const auto *nsi_366 = buffer.data(nsi + 366);
    const auto *nsi_367 = buffer.data(nsi + 367);
    const auto *nsi_369 = buffer.data(nsi + 369);
    const auto *nsi_370 = buffer.data(nsi + 370);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pc_x, pc_y, pc_z, msi_173, msi_290, \
                         nsh0_212, nsh0_220, nsh1_212, nsh1_220, nsi_285, nsi_286, \
                         nsi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_16 * msi_173[k]
                   + f_3 * pc_y[k] * nsi_285[k];

        t_369[k] = f_6 * nsh0_212[k]
                   - f_7 * nsh1_212[k]
                   + f_3 * pc_z[k] * nsi_285[k];

        t_370[k] = f_23 * msi_290[k]
                   + f_6 * nsh0_220[k]
                   - f_7 * nsh1_220[k]
                   + f_3 * pc_x[k] * nsi_290[k];

        t_371[k] = f_3 * pc_z[k] * nsi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, msi_177, nsh0_213, nsh0_215, \
                         nsh1_213, nsh1_215, nsi_287, nsi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * nsh0_213[k]
                   - f_5 * nsh1_213[k]
                   + f_3 * pc_z[k] * nsi_287[k];

        t_373[k] = f_16 * msi_177[k]
                   + f_3 * pc_y[k] * nsi_289[k];

        t_374[k] = f_8 * nsh0_215[k]
                   - f_9 * nsh1_215[k]
                   + f_3 * pc_z[k] * nsi_289[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_z, msi_295, nsh0_216, nsh0_225, \
                         nsh1_216, nsh1_225, nsi_290, nsi_291, \
                         nsi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_23 * msi_295[k]
                   + f_4 * nsh0_225[k]
                   - f_5 * nsh1_225[k]
                   + f_3 * pc_x[k] * nsi_295[k];

        t_376[k] = f_3 * pc_z[k] * nsi_290[k];

        t_377[k] = f_4 * nsh0_216[k]
                   - f_5 * nsh1_216[k]
                   + f_3 * pc_z[k] * nsi_291[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, msi_182, msi_301, \
                         nsh0_217, nsh0_219, nsh1_217, nsh1_219, nsi_292, nsi_294, \
                         nsi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * nsh0_217[k]
                   - f_7 * nsh1_217[k]
                   + f_3 * pc_z[k] * nsi_292[k];

        t_379[k] = f_16 * msi_182[k]
                   + f_3 * pc_y[k] * nsi_294[k];

        t_380[k] = f_10 * nsh0_219[k]
                   - f_11 * nsh1_219[k]
                   + f_3 * pc_z[k] * nsi_294[k];

        t_381[k] = f_23 * msi_301[k]
                   + f_3 * pc_x[k] * nsi_301[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_z, msi_303, msi_304, \
                         msi_305, msi_306, nsi_295, nsi_303, nsi_304, nsi_305, \
                         nsi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_z[k] * nsi_295[k];

        t_383[k] = f_23 * msi_303[k]
                   + f_3 * pc_x[k] * nsi_303[k];

        t_384[k] = f_23 * msi_304[k]
                   + f_3 * pc_x[k] * nsi_304[k];

        t_385[k] = f_23 * msi_305[k]
                   + f_3 * pc_x[k] * nsi_305[k];

        t_386[k] = f_23 * msi_306[k]
                   + f_3 * pc_x[k] * nsi_306[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pc_y, pc_z, msi_189, msi_307, \
                         nsh0_225, nsh1_225, nsi_301, nsi_302, \
                         nsi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_23 * msi_307[k]
                   + f_3 * pc_x[k] * nsi_307[k];

        t_388[k] = f_16 * msi_189[k]
                   + f_1 * nsh0_225[k]
                   - f_2 * nsh1_225[k]
                   + f_3 * pc_y[k] * nsi_301[k];

        t_389[k] = f_3 * pc_z[k] * nsi_301[k];

        t_390[k] = f_4 * nsh0_225[k]
                   - f_5 * nsh1_225[k]
                   + f_3 * pc_z[k] * nsi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_z, nsh0_226, nsh0_227, nsh0_228, nsh1_226, \
                         nsh1_227, nsh1_228, nsi_303, nsi_304, \
                         nsi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_6 * nsh0_226[k]
                   - f_7 * nsh1_226[k]
                   + f_3 * pc_z[k] * nsi_303[k];

        t_392[k] = f_8 * nsh0_227[k]
                   - f_9 * nsh1_227[k]
                   + f_3 * pc_z[k] * nsi_304[k];

        t_393[k] = f_10 * nsh0_228[k]
                   - f_11 * nsh1_228[k]
                   + f_3 * pc_z[k] * nsi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_z, pc_y, pc_z, msk0_216, msi_195, \
                         msi_196, msk1_216, nsh0_230, nsh1_230, nsi_307, \
                         nsi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * msi_195[k]
                   + f_3 * pc_y[k] * nsi_307[k];

        t_395[k] = f_1 * nsh0_230[k]
                   - f_2 * nsh1_230[k]
                   + f_3 * pc_z[k] * nsi_307[k];

        t_396[k] = pa_z[k] * msk0_216[k]
                   - f_12 * pc_z[k] * msk1_216[k];

        t_397[k] = f_15 * msi_196[k]
                   + f_3 * pc_y[k] * nsi_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_z, pc_y, pc_z, msk0_219, msi_168, msi_198, \
                         msk1_219, nsi_308, nsi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * msi_168[k]
                   + f_3 * pc_z[k] * nsi_308[k];

        t_399[k] = pa_z[k] * msk0_219[k]
                   - f_12 * pc_z[k] * msk1_219[k];

        t_400[k] = f_15 * msi_198[k]
                   + f_3 * pc_y[k] * nsi_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_z, pc_x, pc_z, msk0_222, msi_171, msi_313, \
                         msk1_222, nsh0_236, nsh1_236, nsi_311, \
                         nsi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_23 * msi_313[k]
                   + f_10 * nsh0_236[k]
                   - f_11 * nsh1_236[k]
                   + f_3 * pc_x[k] * nsi_313[k];

        t_402[k] = pa_z[k] * msk0_222[k]
                   - f_12 * pc_z[k] * msk1_222[k];

        t_403[k] = f_13 * msi_171[k]
                   + f_3 * pc_z[k] * nsi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_z, pc_x, pc_y, pc_z, msk0_226, msi_201, \
                         msi_317, msk1_226, nsh0_240, nsh1_240, nsi_313, \
                         nsi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * msi_201[k]
                   + f_3 * pc_y[k] * nsi_313[k];

        t_405[k] = f_23 * msi_317[k]
                   + f_8 * nsh0_240[k]
                   - f_9 * nsh1_240[k]
                   + f_3 * pc_x[k] * nsi_317[k];

        t_406[k] = pa_z[k] * msk0_226[k]
                   - f_12 * pc_z[k] * msk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_z, pc_y, pc_z, msk0_228, msi_174, msi_175, \
                         msi_205, msk1_228, nsi_314, nsi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * msi_174[k]
                   + f_3 * pc_z[k] * nsi_314[k];

        t_408[k] = pa_z[k] * msk0_228[k]
                   + f_14 * msi_175[k]
                   - f_12 * pc_z[k] * msk1_228[k];

        t_409[k] = f_15 * msi_205[k]
                   + f_3 * pc_y[k] * nsi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pa_z, pc_x, pc_z, msk0_231, msi_178, msi_322, \
                         msk1_231, nsh0_245, nsh1_245, nsi_318, \
                         nsi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_23 * msi_322[k]
                   + f_6 * nsh0_245[k]
                   - f_7 * nsh1_245[k]
                   + f_3 * pc_x[k] * nsi_322[k];

        t_411[k] = pa_z[k] * msk0_231[k]
                   - f_12 * pc_z[k] * msk1_231[k];

        t_412[k] = f_13 * msi_178[k]
                   + f_3 * pc_z[k] * nsi_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pa_z, pc_y, pc_z, msk0_233, msk0_234, msi_179, \
                         msi_180, msi_210, msk1_233, msk1_234, \
                         nsi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_z[k] * msk0_233[k]
                   + f_14 * msi_179[k]
                   - f_12 * pc_z[k] * msk1_233[k];

        t_414[k] = pa_z[k] * msk0_234[k]
                   + f_15 * msi_180[k]
                   - f_12 * pc_z[k] * msk1_234[k];

        t_415[k] = f_15 * msi_210[k]
                   + f_3 * pc_y[k] * nsi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, msi_328, msi_329, msi_330, msi_331, \
                         nsh0_251, nsh1_251, nsi_328, nsi_329, nsi_330, \
                         nsi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_23 * msi_328[k]
                   + f_4 * nsh0_251[k]
                   - f_5 * nsh1_251[k]
                   + f_3 * pc_x[k] * nsi_328[k];

        t_417[k] = f_23 * msi_329[k]
                   + f_3 * pc_x[k] * nsi_329[k];

        t_418[k] = f_23 * msi_330[k]
                   + f_3 * pc_x[k] * nsi_330[k];

        t_419[k] = f_23 * msi_331[k]
                   + f_3 * pc_x[k] * nsi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, msi_332, msi_333, msi_334, msi_335, \
                         nsi_332, nsi_333, nsi_334, nsi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_23 * msi_332[k]
                   + f_3 * pc_x[k] * nsi_332[k];

        t_421[k] = f_23 * msi_333[k]
                   + f_3 * pc_x[k] * nsi_333[k];

        t_422[k] = f_23 * msi_334[k]
                   + f_3 * pc_x[k] * nsi_334[k];

        t_423[k] = f_23 * msi_335[k]
                   + f_3 * pc_x[k] * nsi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_z, pc_y, pc_z, msk0_244, msi_189, msi_219, \
                         msk1_244, nsh0_248, nsh1_248, nsi_329, \
                         nsi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * msk0_244[k]
                   - f_12 * pc_z[k] * msk1_244[k];

        t_425[k] = f_13 * msi_189[k]
                   + f_3 * pc_z[k] * nsi_329[k];

        t_426[k] = f_15 * msi_219[k]
                   + f_10 * nsh0_248[k]
                   - f_11 * nsh1_248[k]
                   + f_3 * pc_y[k] * nsi_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, msi_220, msi_221, msi_222, nsh0_249, \
                         nsh0_250, nsh0_251, nsh1_249, nsh1_250, nsh1_251, nsi_332, nsi_333, \
                         nsi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * msi_220[k]
                   + f_8 * nsh0_249[k]
                   - f_9 * nsh1_249[k]
                   + f_3 * pc_y[k] * nsi_332[k];

        t_428[k] = f_15 * msi_221[k]
                   + f_6 * nsh0_250[k]
                   - f_7 * nsh1_250[k]
                   + f_3 * pc_y[k] * nsi_333[k];

        t_429[k] = f_15 * msi_222[k]
                   + f_4 * nsh0_251[k]
                   - f_5 * nsh1_251[k]
                   + f_3 * pc_y[k] * nsi_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, msi_195, msi_223, msi_336, \
                         nsh0_251, nsh0_252, nsh1_251, nsh1_252, nsi_335, \
                         nsi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * msi_223[k]
                   + f_3 * pc_y[k] * nsi_335[k];

        t_431[k] = f_13 * msi_195[k]
                   + f_1 * nsh0_251[k]
                   - f_2 * nsh1_251[k]
                   + f_3 * pc_z[k] * nsi_335[k];

        t_432[k] = f_23 * msi_336[k]
                   + f_1 * nsh0_252[k]
                   - f_2 * nsh1_252[k]
                   + f_3 * pc_x[k] * nsi_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, msi_196, msi_224, \
                         msi_226, msi_339, nsh0_255, nsh1_255, nsi_336, nsi_338, \
                         nsi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * msi_224[k]
                   + f_3 * pc_y[k] * nsi_336[k];

        t_434[k] = f_14 * msi_196[k]
                   + f_3 * pc_z[k] * nsi_336[k];

        t_435[k] = f_23 * msi_339[k]
                   + f_10 * nsh0_255[k]
                   - f_11 * nsh1_255[k]
                   + f_3 * pc_x[k] * nsi_339[k];

        t_436[k] = f_14 * msi_226[k]
                   + f_3 * pc_y[k] * nsi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, msi_199, msi_341, msi_342, nsh0_257, \
                         nsh0_258, nsh1_257, nsh1_258, nsi_339, nsi_341, \
                         nsi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_23 * msi_341[k]
                   + f_10 * nsh0_257[k]
                   - f_11 * nsh1_257[k]
                   + f_3 * pc_x[k] * nsi_341[k];

        t_438[k] = f_23 * msi_342[k]
                   + f_8 * nsh0_258[k]
                   - f_9 * nsh1_258[k]
                   + f_3 * pc_x[k] * nsi_342[k];

        t_439[k] = f_14 * msi_199[k]
                   + f_3 * pc_z[k] * nsi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, msi_229, msi_345, msi_346, nsh0_261, \
                         nsh0_262, nsh1_261, nsh1_262, nsi_341, nsi_345, \
                         nsi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * msi_229[k]
                   + f_3 * pc_y[k] * nsi_341[k];

        t_441[k] = f_23 * msi_345[k]
                   + f_8 * nsh0_261[k]
                   - f_9 * nsh1_261[k]
                   + f_3 * pc_x[k] * nsi_345[k];

        t_442[k] = f_23 * msi_346[k]
                   + f_6 * nsh0_262[k]
                   - f_7 * nsh1_262[k]
                   + f_3 * pc_x[k] * nsi_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, msi_202, msi_233, msi_348, \
                         nsh0_264, nsh1_264, nsi_342, nsi_345, \
                         nsi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * msi_202[k]
                   + f_3 * pc_z[k] * nsi_342[k];

        t_444[k] = f_23 * msi_348[k]
                   + f_6 * nsh0_264[k]
                   - f_7 * nsh1_264[k]
                   + f_3 * pc_x[k] * nsi_348[k];

        t_445[k] = f_14 * msi_233[k]
                   + f_3 * pc_y[k] * nsi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, msi_206, msi_350, msi_351, nsh0_266, \
                         nsh0_267, nsh1_266, nsh1_267, nsi_346, nsi_350, \
                         nsi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_23 * msi_350[k]
                   + f_6 * nsh0_266[k]
                   - f_7 * nsh1_266[k]
                   + f_3 * pc_x[k] * nsi_350[k];

        t_447[k] = f_23 * msi_351[k]
                   + f_4 * nsh0_267[k]
                   - f_5 * nsh1_267[k]
                   + f_3 * pc_x[k] * nsi_351[k];

        t_448[k] = f_14 * msi_206[k]
                   + f_3 * pc_z[k] * nsi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, msi_238, msi_353, msi_354, nsh0_269, \
                         nsh0_270, nsh1_269, nsh1_270, nsi_350, nsi_353, \
                         nsi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_23 * msi_353[k]
                   + f_4 * nsh0_269[k]
                   - f_5 * nsh1_269[k]
                   + f_3 * pc_x[k] * nsi_353[k];

        t_450[k] = f_23 * msi_354[k]
                   + f_4 * nsh0_270[k]
                   - f_5 * nsh1_270[k]
                   + f_3 * pc_x[k] * nsi_354[k];

        t_451[k] = f_14 * msi_238[k]
                   + f_3 * pc_y[k] * nsi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, msi_356, msi_357, msi_358, msi_359, \
                         nsh0_272, nsh1_272, nsi_356, nsi_357, nsi_358, \
                         nsi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_23 * msi_356[k]
                   + f_4 * nsh0_272[k]
                   - f_5 * nsh1_272[k]
                   + f_3 * pc_x[k] * nsi_356[k];

        t_453[k] = f_23 * msi_357[k]
                   + f_3 * pc_x[k] * nsi_357[k];

        t_454[k] = f_23 * msi_358[k]
                   + f_3 * pc_x[k] * nsi_358[k];

        t_455[k] = f_23 * msi_359[k]
                   + f_3 * pc_x[k] * nsi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, msi_360, msi_361, msi_362, msi_363, \
                         nsi_360, nsi_361, nsi_362, nsi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_23 * msi_360[k]
                   + f_3 * pc_x[k] * nsi_360[k];

        t_457[k] = f_23 * msi_361[k]
                   + f_3 * pc_x[k] * nsi_361[k];

        t_458[k] = f_23 * msi_362[k]
                   + f_3 * pc_x[k] * nsi_362[k];

        t_459[k] = f_23 * msi_363[k]
                   + f_3 * pc_x[k] * nsi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, msi_217, msi_245, msi_247, nsh0_267, \
                         nsh0_269, nsh1_267, nsh1_269, nsi_357, \
                         nsi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * msi_245[k]
                   + f_1 * nsh0_267[k]
                   - f_2 * nsh1_267[k]
                   + f_3 * pc_y[k] * nsi_357[k];

        t_461[k] = f_14 * msi_217[k]
                   + f_3 * pc_z[k] * nsi_357[k];

        t_462[k] = f_14 * msi_247[k]
                   + f_10 * nsh0_269[k]
                   - f_11 * nsh1_269[k]
                   + f_3 * pc_y[k] * nsi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, msi_248, msi_249, msi_250, nsh0_270, \
                         nsh0_271, nsh0_272, nsh1_270, nsh1_271, nsh1_272, nsi_360, nsi_361, \
                         nsi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * msi_248[k]
                   + f_8 * nsh0_270[k]
                   - f_9 * nsh1_270[k]
                   + f_3 * pc_y[k] * nsi_360[k];

        t_464[k] = f_14 * msi_249[k]
                   + f_6 * nsh0_271[k]
                   - f_7 * nsh1_271[k]
                   + f_3 * pc_y[k] * nsi_361[k];

        t_465[k] = f_14 * msi_250[k]
                   + f_4 * nsh0_272[k]
                   - f_5 * nsh1_272[k]
                   + f_3 * pc_y[k] * nsi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_y, pc_y, pc_z, msk0_324, msi_223, \
                         msi_251, msi_252, msk1_324, nsh0_272, nsh1_272, nsi_363, \
                         nsi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * msi_251[k]
                   + f_3 * pc_y[k] * nsi_363[k];

        t_467[k] = f_14 * msi_223[k]
                   + f_1 * nsh0_272[k]
                   - f_2 * nsh1_272[k]
                   + f_3 * pc_z[k] * nsi_363[k];

        t_468[k] = pa_y[k] * msk0_324[k]
                   - f_12 * pc_y[k] * msk1_324[k];

        t_469[k] = f_13 * msi_252[k]
                   + f_3 * pc_y[k] * nsi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_y, pc_y, pc_z, msk0_327, msk0_329, \
                         msi_224, msi_253, msi_254, msk1_327, msk1_329, nsi_364, \
                         nsi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * msi_224[k]
                   + f_3 * pc_z[k] * nsi_364[k];

        t_471[k] = pa_y[k] * msk0_327[k]
                   + f_14 * msi_253[k]
                   - f_12 * pc_y[k] * msk1_327[k];

        t_472[k] = f_13 * msi_254[k]
                   + f_3 * pc_y[k] * nsi_366[k];

        t_473[k] = pa_y[k] * msk0_329[k]
                   - f_12 * pc_y[k] * msk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_y, pc_y, pc_z, msk0_330, msk0_333, \
                         msi_227, msi_255, msi_257, msk1_330, msk1_333, nsi_367, \
                         nsi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_y[k] * msk0_330[k]
                   + f_15 * msi_255[k]
                   - f_12 * pc_y[k] * msk1_330[k];

        t_475[k] = f_15 * msi_227[k]
                   + f_3 * pc_z[k] * nsi_367[k];

        t_476[k] = f_13 * msi_257[k]
                   + f_3 * pc_y[k] * nsi_369[k];

        t_477[k] = pa_y[k] * msk0_333[k]
                   - f_12 * pc_y[k] * msk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_y, pc_y, pc_z, msk0_334, msk0_336, msi_230, \
                         msi_258, msi_260, msk1_334, msk1_336, \
                         nsi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pa_y[k] * msk0_334[k]
                   + f_16 * msi_258[k]
                   - f_12 * pc_y[k] * msk1_334[k];

        t_479[k] = f_15 * msi_230[k]
                   + f_3 * pc_z[k] * nsi_370[k];

        t_480[k] = pa_y[k] * msk0_336[k]
                   + f_14 * msi_260[k]
                   - f_12 * pc_y[k] * msk1_336[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_338 = buffer.data(msk0 + 338);
    const auto *msk0_339 = buffer.data(msk0 + 339);
    const auto *msk0_341 = buffer.data(msk0 + 341);
    const auto *msk0_342 = buffer.data(msk0 + 342);
    const auto *msk0_344 = buffer.data(msk0 + 344);
    const auto *msk0_359 = buffer.data(msk0 + 359);
    const auto *msk0_360 = buffer.data(msk0 + 360);
    const auto *msk0_363 = buffer.data(msk0 + 363);
    const auto *msk0_366 = buffer.data(msk0 + 366);
    const auto *msk0_370 = buffer.data(msk0 + 370);
    const auto *msk0_372 = buffer.data(msk0 + 372);
    const auto *msk0_375 = buffer.data(msk0 + 375);
    const auto *msk0_377 = buffer.data(msk0 + 377);
    const auto *msk0_378 = buffer.data(msk0 + 378);

    const auto *msi_234 = buffer.data(msi + 234);
    const auto *msi_245 = buffer.data(msi + 245);
    const auto *msi_252 = buffer.data(msi + 252);
    const auto *msi_261 = buffer.data(msi + 261);
    const auto *msi_262 = buffer.data(msi + 262);
    const auto *msi_264 = buffer.data(msi + 264);
    const auto *msi_265 = buffer.data(msi + 265);
    const auto *msi_266 = buffer.data(msi + 266);
    const auto *msi_273 = buffer.data(msi + 273);
    const auto *msi_275 = buffer.data(msi + 275);
    const auto *msi_276 = buffer.data(msi + 276);
    const auto *msi_277 = buffer.data(msi + 277);
    const auto *msi_278 = buffer.data(msi + 278);
    const auto *msi_279 = buffer.data(msi + 279);
    const auto *msi_280 = buffer.data(msi + 280);
    const auto *msi_283 = buffer.data(msi + 283);
    const auto *msi_285 = buffer.data(msi + 285);
    const auto *msi_286 = buffer.data(msi + 286);
    const auto *msi_287 = buffer.data(msi + 287);
    const auto *msi_289 = buffer.data(msi + 289);
    const auto *msi_290 = buffer.data(msi + 290);
    const auto *msi_291 = buffer.data(msi + 291);
    const auto *msi_292 = buffer.data(msi + 292);
    const auto *msi_294 = buffer.data(msi + 294);
    const auto *msi_301 = buffer.data(msi + 301);
    const auto *msi_307 = buffer.data(msi + 307);
    const auto *msi_308 = buffer.data(msi + 308);
    const auto *msi_310 = buffer.data(msi + 310);
    const auto *msi_313 = buffer.data(msi + 313);
    const auto *msi_317 = buffer.data(msi + 317);
    const auto *msi_322 = buffer.data(msi + 322);
    const auto *msi_385 = buffer.data(msi + 385);
    const auto *msi_386 = buffer.data(msi + 386);
    const auto *msi_387 = buffer.data(msi + 387);
    const auto *msi_388 = buffer.data(msi + 388);
    const auto *msi_389 = buffer.data(msi + 389);
    const auto *msi_390 = buffer.data(msi + 390);
    const auto *msi_391 = buffer.data(msi + 391);
    const auto *msi_392 = buffer.data(msi + 392);
    const auto *msi_397 = buffer.data(msi + 397);
    const auto *msi_401 = buffer.data(msi + 401);
    const auto *msi_406 = buffer.data(msi + 406);
    const auto *msi_412 = buffer.data(msi + 412);
    const auto *msi_413 = buffer.data(msi + 413);
    const auto *msi_414 = buffer.data(msi + 414);
    const auto *msi_415 = buffer.data(msi + 415);
    const auto *msi_416 = buffer.data(msi + 416);
    const auto *msi_417 = buffer.data(msi + 417);
    const auto *msi_419 = buffer.data(msi + 419);
    const auto *msi_420 = buffer.data(msi + 420);
    const auto *msi_423 = buffer.data(msi + 423);
    const auto *msi_426 = buffer.data(msi + 426);
    const auto *msi_430 = buffer.data(msi + 430);
    const auto *msi_435 = buffer.data(msi + 435);
    const auto *msi_441 = buffer.data(msi + 441);
    const auto *msi_443 = buffer.data(msi + 443);
    const auto *msi_444 = buffer.data(msi + 444);
    const auto *msi_445 = buffer.data(msi + 445);
    const auto *msi_446 = buffer.data(msi + 446);
    const auto *msi_447 = buffer.data(msi + 447);
    const auto *msi_453 = buffer.data(msi + 453);
    const auto *msi_457 = buffer.data(msi + 457);
    const auto *msi_462 = buffer.data(msi + 462);
    const auto *msi_468 = buffer.data(msi + 468);
    const auto *msi_469 = buffer.data(msi + 469);
    const auto *msi_470 = buffer.data(msi + 470);
    const auto *msi_471 = buffer.data(msi + 471);

    const auto *msk1_338 = buffer.data(msk1 + 338);
    const auto *msk1_339 = buffer.data(msk1 + 339);
    const auto *msk1_341 = buffer.data(msk1 + 341);
    const auto *msk1_342 = buffer.data(msk1 + 342);
    const auto *msk1_344 = buffer.data(msk1 + 344);
    const auto *msk1_359 = buffer.data(msk1 + 359);
    const auto *msk1_360 = buffer.data(msk1 + 360);
    const auto *msk1_363 = buffer.data(msk1 + 363);
    const auto *msk1_366 = buffer.data(msk1 + 366);
    const auto *msk1_370 = buffer.data(msk1 + 370);
    const auto *msk1_372 = buffer.data(msk1 + 372);
    const auto *msk1_375 = buffer.data(msk1 + 375);
    const auto *msk1_377 = buffer.data(msk1 + 377);
    const auto *msk1_378 = buffer.data(msk1 + 378);

    const auto *nsh0_288 = buffer.data(nsh0 + 288);
    const auto *nsh0_290 = buffer.data(nsh0 + 290);
    const auto *nsh0_291 = buffer.data(nsh0 + 291);
    const auto *nsh0_292 = buffer.data(nsh0 + 292);
    const auto *nsh0_293 = buffer.data(nsh0 + 293);
    const auto *nsh0_294 = buffer.data(nsh0 + 294);
    const auto *nsh0_295 = buffer.data(nsh0 + 295);
    const auto *nsh0_296 = buffer.data(nsh0 + 296);
    const auto *nsh0_297 = buffer.data(nsh0 + 297);
    const auto *nsh0_298 = buffer.data(nsh0 + 298);
    const auto *nsh0_299 = buffer.data(nsh0 + 299);
    const auto *nsh0_300 = buffer.data(nsh0 + 300);
    const auto *nsh0_301 = buffer.data(nsh0 + 301);
    const auto *nsh0_302 = buffer.data(nsh0 + 302);
    const auto *nsh0_303 = buffer.data(nsh0 + 303);
    const auto *nsh0_308 = buffer.data(nsh0 + 308);
    const auto *nsh0_309 = buffer.data(nsh0 + 309);
    const auto *nsh0_310 = buffer.data(nsh0 + 310);
    const auto *nsh0_311 = buffer.data(nsh0 + 311);
    const auto *nsh0_312 = buffer.data(nsh0 + 312);
    const auto *nsh0_313 = buffer.data(nsh0 + 313);
    const auto *nsh0_314 = buffer.data(nsh0 + 314);
    const auto *nsh0_315 = buffer.data(nsh0 + 315);
    const auto *nsh0_317 = buffer.data(nsh0 + 317);
    const auto *nsh0_318 = buffer.data(nsh0 + 318);
    const auto *nsh0_320 = buffer.data(nsh0 + 320);
    const auto *nsh0_321 = buffer.data(nsh0 + 321);
    const auto *nsh0_322 = buffer.data(nsh0 + 322);
    const auto *nsh0_324 = buffer.data(nsh0 + 324);
    const auto *nsh0_325 = buffer.data(nsh0 + 325);
    const auto *nsh0_330 = buffer.data(nsh0 + 330);
    const auto *nsh0_331 = buffer.data(nsh0 + 331);
    const auto *nsh0_332 = buffer.data(nsh0 + 332);
    const auto *nsh0_333 = buffer.data(nsh0 + 333);
    const auto *nsh0_335 = buffer.data(nsh0 + 335);
    const auto *nsh0_341 = buffer.data(nsh0 + 341);
    const auto *nsh0_345 = buffer.data(nsh0 + 345);
    const auto *nsh0_350 = buffer.data(nsh0 + 350);
    const auto *nsh0_356 = buffer.data(nsh0 + 356);

    const auto *nsh1_288 = buffer.data(nsh1 + 288);
    const auto *nsh1_290 = buffer.data(nsh1 + 290);
    const auto *nsh1_291 = buffer.data(nsh1 + 291);
    const auto *nsh1_292 = buffer.data(nsh1 + 292);
    const auto *nsh1_293 = buffer.data(nsh1 + 293);
    const auto *nsh1_294 = buffer.data(nsh1 + 294);
    const auto *nsh1_295 = buffer.data(nsh1 + 295);
    const auto *nsh1_296 = buffer.data(nsh1 + 296);
    const auto *nsh1_297 = buffer.data(nsh1 + 297);
    const auto *nsh1_298 = buffer.data(nsh1 + 298);
    const auto *nsh1_299 = buffer.data(nsh1 + 299);
    const auto *nsh1_300 = buffer.data(nsh1 + 300);
    const auto *nsh1_301 = buffer.data(nsh1 + 301);
    const auto *nsh1_302 = buffer.data(nsh1 + 302);
    const auto *nsh1_303 = buffer.data(nsh1 + 303);
    const auto *nsh1_308 = buffer.data(nsh1 + 308);
    const auto *nsh1_309 = buffer.data(nsh1 + 309);
    const auto *nsh1_310 = buffer.data(nsh1 + 310);
    const auto *nsh1_311 = buffer.data(nsh1 + 311);
    const auto *nsh1_312 = buffer.data(nsh1 + 312);
    const auto *nsh1_313 = buffer.data(nsh1 + 313);
    const auto *nsh1_314 = buffer.data(nsh1 + 314);
    const auto *nsh1_315 = buffer.data(nsh1 + 315);
    const auto *nsh1_317 = buffer.data(nsh1 + 317);
    const auto *nsh1_318 = buffer.data(nsh1 + 318);
    const auto *nsh1_320 = buffer.data(nsh1 + 320);
    const auto *nsh1_321 = buffer.data(nsh1 + 321);
    const auto *nsh1_322 = buffer.data(nsh1 + 322);
    const auto *nsh1_324 = buffer.data(nsh1 + 324);
    const auto *nsh1_325 = buffer.data(nsh1 + 325);
    const auto *nsh1_330 = buffer.data(nsh1 + 330);
    const auto *nsh1_331 = buffer.data(nsh1 + 331);
    const auto *nsh1_332 = buffer.data(nsh1 + 332);
    const auto *nsh1_333 = buffer.data(nsh1 + 333);
    const auto *nsh1_335 = buffer.data(nsh1 + 335);
    const auto *nsh1_341 = buffer.data(nsh1 + 341);
    const auto *nsh1_345 = buffer.data(nsh1 + 345);
    const auto *nsh1_350 = buffer.data(nsh1 + 350);
    const auto *nsh1_356 = buffer.data(nsh1 + 356);

    const auto *nsi_373 = buffer.data(nsi + 373);
    const auto *nsi_374 = buffer.data(nsi + 374);
    const auto *nsi_378 = buffer.data(nsi + 378);
    const auto *nsi_385 = buffer.data(nsi + 385);
    const auto *nsi_386 = buffer.data(nsi + 386);
    const auto *nsi_387 = buffer.data(nsi + 387);
    const auto *nsi_388 = buffer.data(nsi + 388);
    const auto *nsi_389 = buffer.data(nsi + 389);
    const auto *nsi_390 = buffer.data(nsi + 390);
    const auto *nsi_391 = buffer.data(nsi + 391);
    const auto *nsi_392 = buffer.data(nsi + 392);
    const auto *nsi_393 = buffer.data(nsi + 393);
    const auto *nsi_394 = buffer.data(nsi + 394);
    const auto *nsi_395 = buffer.data(nsi + 395);
    const auto *nsi_396 = buffer.data(nsi + 396);
    const auto *nsi_397 = buffer.data(nsi + 397);
    const auto *nsi_398 = buffer.data(nsi + 398);
    const auto *nsi_399 = buffer.data(nsi + 399);
    const auto *nsi_400 = buffer.data(nsi + 400);
    const auto *nsi_401 = buffer.data(nsi + 401);
    const auto *nsi_402 = buffer.data(nsi + 402);
    const auto *nsi_403 = buffer.data(nsi + 403);
    const auto *nsi_404 = buffer.data(nsi + 404);
    const auto *nsi_405 = buffer.data(nsi + 405);
    const auto *nsi_406 = buffer.data(nsi + 406);
    const auto *nsi_412 = buffer.data(nsi + 412);
    const auto *nsi_413 = buffer.data(nsi + 413);
    const auto *nsi_414 = buffer.data(nsi + 414);
    const auto *nsi_415 = buffer.data(nsi + 415);
    const auto *nsi_416 = buffer.data(nsi + 416);
    const auto *nsi_417 = buffer.data(nsi + 417);
    const auto *nsi_418 = buffer.data(nsi + 418);
    const auto *nsi_419 = buffer.data(nsi + 419);
    const auto *nsi_420 = buffer.data(nsi + 420);
    const auto *nsi_421 = buffer.data(nsi + 421);
    const auto *nsi_422 = buffer.data(nsi + 422);
    const auto *nsi_423 = buffer.data(nsi + 423);
    const auto *nsi_425 = buffer.data(nsi + 425);
    const auto *nsi_426 = buffer.data(nsi + 426);
    const auto *nsi_427 = buffer.data(nsi + 427);
    const auto *nsi_429 = buffer.data(nsi + 429);
    const auto *nsi_430 = buffer.data(nsi + 430);
    const auto *nsi_431 = buffer.data(nsi + 431);
    const auto *nsi_432 = buffer.data(nsi + 432);
    const auto *nsi_434 = buffer.data(nsi + 434);
    const auto *nsi_435 = buffer.data(nsi + 435);
    const auto *nsi_441 = buffer.data(nsi + 441);
    const auto *nsi_442 = buffer.data(nsi + 442);
    const auto *nsi_443 = buffer.data(nsi + 443);
    const auto *nsi_444 = buffer.data(nsi + 444);
    const auto *nsi_445 = buffer.data(nsi + 445);
    const auto *nsi_446 = buffer.data(nsi + 446);
    const auto *nsi_447 = buffer.data(nsi + 447);
    const auto *nsi_448 = buffer.data(nsi + 448);
    const auto *nsi_450 = buffer.data(nsi + 450);
    const auto *nsi_451 = buffer.data(nsi + 451);
    const auto *nsi_453 = buffer.data(nsi + 453);
    const auto *nsi_454 = buffer.data(nsi + 454);
    const auto *nsi_457 = buffer.data(nsi + 457);
    const auto *nsi_458 = buffer.data(nsi + 458);
    const auto *nsi_462 = buffer.data(nsi + 462);
    const auto *nsi_468 = buffer.data(nsi + 468);
    const auto *nsi_469 = buffer.data(nsi + 469);
    const auto *nsi_470 = buffer.data(nsi + 470);
    const auto *nsi_471 = buffer.data(nsi + 471);

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_y, pc_y, pc_z, msk0_338, msk0_339, \
                         msi_234, msi_261, msi_262, msk1_338, msk1_339, nsi_373, \
                         nsi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * msi_261[k]
                   + f_3 * pc_y[k] * nsi_373[k];

        t_482[k] = pa_y[k] * msk0_338[k]
                   - f_12 * pc_y[k] * msk1_338[k];

        t_483[k] = pa_y[k] * msk0_339[k]
                   + f_17 * msi_262[k]
                   - f_12 * pc_y[k] * msk1_339[k];

        t_484[k] = f_15 * msi_234[k]
                   + f_3 * pc_z[k] * nsi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_y, pc_y, msk0_341, msk0_342, msk0_344, \
                         msi_264, msi_265, msi_266, msk1_341, msk1_342, msk1_344, \
                         nsi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_y[k] * msk0_341[k]
                   + f_15 * msi_264[k]
                   - f_12 * pc_y[k] * msk1_341[k];

        t_486[k] = pa_y[k] * msk0_342[k]
                   + f_14 * msi_265[k]
                   - f_12 * pc_y[k] * msk1_342[k];

        t_487[k] = f_13 * msi_266[k]
                   + f_3 * pc_y[k] * nsi_378[k];

        t_488[k] = pa_y[k] * msk0_344[k]
                   - f_12 * pc_y[k] * msk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, msi_385, msi_386, msi_387, \
                         msi_388, msi_389, nsi_385, nsi_386, nsi_387, nsi_388, \
                         nsi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_23 * msi_385[k]
                   + f_3 * pc_x[k] * nsi_385[k];

        t_490[k] = f_23 * msi_386[k]
                   + f_3 * pc_x[k] * nsi_386[k];

        t_491[k] = f_23 * msi_387[k]
                   + f_3 * pc_x[k] * nsi_387[k];

        t_492[k] = f_23 * msi_388[k]
                   + f_3 * pc_x[k] * nsi_388[k];

        t_493[k] = f_23 * msi_389[k]
                   + f_3 * pc_x[k] * nsi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, msi_245, msi_273, \
                         msi_390, msi_391, nsh0_288, nsh1_288, nsi_385, nsi_390, \
                         nsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_23 * msi_390[k]
                   + f_3 * pc_x[k] * nsi_390[k];

        t_495[k] = f_23 * msi_391[k]
                   + f_3 * pc_x[k] * nsi_391[k];

        t_496[k] = f_13 * msi_273[k]
                   + f_1 * nsh0_288[k]
                   - f_2 * nsh1_288[k]
                   + f_3 * pc_y[k] * nsi_385[k];

        t_497[k] = f_15 * msi_245[k]
                   + f_3 * pc_z[k] * nsi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, msi_275, msi_276, msi_277, nsh0_290, \
                         nsh0_291, nsh0_292, nsh1_290, nsh1_291, nsh1_292, nsi_387, nsi_388, \
                         nsi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * msi_275[k]
                   + f_10 * nsh0_290[k]
                   - f_11 * nsh1_290[k]
                   + f_3 * pc_y[k] * nsi_387[k];

        t_499[k] = f_13 * msi_276[k]
                   + f_8 * nsh0_291[k]
                   - f_9 * nsh1_291[k]
                   + f_3 * pc_y[k] * nsi_388[k];

        t_500[k] = f_13 * msi_277[k]
                   + f_6 * nsh0_292[k]
                   - f_7 * nsh1_292[k]
                   + f_3 * pc_y[k] * nsi_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pc_y, msk0_359, msi_278, msi_279, \
                         msk1_359, nsh0_293, nsh1_293, nsi_390, \
                         nsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * msi_278[k]
                   + f_4 * nsh0_293[k]
                   - f_5 * nsh1_293[k]
                   + f_3 * pc_y[k] * nsi_390[k];

        t_502[k] = f_13 * msi_279[k]
                   + f_3 * pc_y[k] * nsi_391[k];

        t_503[k] = pa_y[k] * msk0_359[k]
                   - f_12 * pc_y[k] * msk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, msi_252, \
                         msi_392, nsh0_294, nsh1_294, nsi_392, nsi_393, \
                         nsi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_23 * msi_392[k]
                   + f_1 * nsh0_294[k]
                   - f_2 * nsh1_294[k]
                   + f_3 * pc_x[k] * nsi_392[k];

        t_505[k] = f_3 * pc_y[k] * nsi_392[k];

        t_506[k] = f_16 * msi_252[k]
                   + f_3 * pc_z[k] * nsi_392[k];

        t_507[k] = f_4 * nsh0_294[k]
                   - f_5 * nsh1_294[k]
                   + f_3 * pc_y[k] * nsi_393[k];

        t_508[k] = f_3 * pc_y[k] * nsi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, msi_397, nsh0_295, nsh0_296, \
                         nsh0_299, nsh1_295, nsh1_296, nsh1_299, nsi_395, nsi_396, \
                         nsi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_23 * msi_397[k]
                   + f_10 * nsh0_299[k]
                   - f_11 * nsh1_299[k]
                   + f_3 * pc_x[k] * nsi_397[k];

        t_510[k] = f_6 * nsh0_295[k]
                   - f_7 * nsh1_295[k]
                   + f_3 * pc_y[k] * nsi_395[k];

        t_511[k] = f_4 * nsh0_296[k]
                   - f_5 * nsh1_296[k]
                   + f_3 * pc_y[k] * nsi_396[k];

        t_512[k] = f_3 * pc_y[k] * nsi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_y, msi_401, nsh0_297, nsh0_298, \
                         nsh0_303, nsh1_297, nsh1_298, nsh1_303, nsi_398, nsi_399, \
                         nsi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_23 * msi_401[k]
                   + f_8 * nsh0_303[k]
                   - f_9 * nsh1_303[k]
                   + f_3 * pc_x[k] * nsi_401[k];

        t_514[k] = f_8 * nsh0_297[k]
                   - f_9 * nsh1_297[k]
                   + f_3 * pc_y[k] * nsi_398[k];

        t_515[k] = f_6 * nsh0_298[k]
                   - f_7 * nsh1_298[k]
                   + f_3 * pc_y[k] * nsi_399[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_x, pc_y, msi_406, nsh0_299, nsh0_308, \
                         nsh1_299, nsh1_308, nsi_400, nsi_401, \
                         nsi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_4 * nsh0_299[k]
                   - f_5 * nsh1_299[k]
                   + f_3 * pc_y[k] * nsi_400[k];

        t_517[k] = f_3 * pc_y[k] * nsi_401[k];

        t_518[k] = f_23 * msi_406[k]
                   + f_6 * nsh0_308[k]
                   - f_7 * nsh1_308[k]
                   + f_3 * pc_x[k] * nsi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, nsh0_300, nsh0_301, nsh0_302, nsh1_300, \
                         nsh1_301, nsh1_302, nsi_402, nsi_403, \
                         nsi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_10 * nsh0_300[k]
                   - f_11 * nsh1_300[k]
                   + f_3 * pc_y[k] * nsi_402[k];

        t_520[k] = f_8 * nsh0_301[k]
                   - f_9 * nsh1_301[k]
                   + f_3 * pc_y[k] * nsi_403[k];

        t_521[k] = f_6 * nsh0_302[k]
                   - f_7 * nsh1_302[k]
                   + f_3 * pc_y[k] * nsi_404[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, msi_412, msi_413, nsh0_303, \
                         nsh0_314, nsh1_303, nsh1_314, nsi_405, nsi_406, nsi_412, \
                         nsi_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_4 * nsh0_303[k]
                   - f_5 * nsh1_303[k]
                   + f_3 * pc_y[k] * nsi_405[k];

        t_523[k] = f_3 * pc_y[k] * nsi_406[k];

        t_524[k] = f_23 * msi_412[k]
                   + f_4 * nsh0_314[k]
                   - f_5 * nsh1_314[k]
                   + f_3 * pc_x[k] * nsi_412[k];

        t_525[k] = f_23 * msi_413[k]
                   + f_3 * pc_x[k] * nsi_413[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, pc_x, pc_y, msi_414, msi_415, \
                         msi_416, msi_417, nsi_412, nsi_414, nsi_415, nsi_416, \
                         nsi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_23 * msi_414[k]
                   + f_3 * pc_x[k] * nsi_414[k];

        t_527[k] = f_23 * msi_415[k]
                   + f_3 * pc_x[k] * nsi_415[k];

        t_528[k] = f_23 * msi_416[k]
                   + f_3 * pc_x[k] * nsi_416[k];

        t_529[k] = f_23 * msi_417[k]
                   + f_3 * pc_x[k] * nsi_417[k];

        t_530[k] = f_3 * pc_y[k] * nsi_412[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pc_x, pc_y, msi_419, nsh0_309, nsh0_310, \
                         nsh1_309, nsh1_310, nsi_413, nsi_414, \
                         nsi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_23 * msi_419[k]
                   + f_3 * pc_x[k] * nsi_419[k];

        t_532[k] = f_1 * nsh0_309[k]
                   - f_2 * nsh1_309[k]
                   + f_3 * pc_y[k] * nsi_413[k];

        t_533[k] = f_19 * nsh0_310[k]
                   - f_20 * nsh1_310[k]
                   + f_3 * pc_y[k] * nsi_414[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, nsh0_311, nsh0_312, nsh0_313, nsh1_311, \
                         nsh1_312, nsh1_313, nsi_415, nsi_416, \
                         nsi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * nsh0_311[k]
                   - f_11 * nsh1_311[k]
                   + f_3 * pc_y[k] * nsi_415[k];

        t_535[k] = f_8 * nsh0_312[k]
                   - f_9 * nsh1_312[k]
                   + f_3 * pc_y[k] * nsi_416[k];

        t_536[k] = f_6 * nsh0_313[k]
                   - f_7 * nsh1_313[k]
                   + f_3 * pc_y[k] * nsi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, msi_279, msi_420, \
                         nsh0_314, nsh0_315, nsh1_314, nsh1_315, nsi_418, nsi_419, \
                         nsi_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * nsh0_314[k]
                   - f_5 * nsh1_314[k]
                   + f_3 * pc_y[k] * nsi_418[k];

        t_538[k] = f_3 * pc_y[k] * nsi_419[k];

        t_539[k] = f_16 * msi_279[k]
                   + f_1 * nsh0_314[k]
                   - f_2 * nsh1_314[k]
                   + f_3 * pc_z[k] * nsi_419[k];

        t_540[k] = f_17 * msi_420[k]
                   + f_1 * nsh0_315[k]
                   - f_2 * nsh1_315[k]
                   + f_3 * pc_x[k] * nsi_420[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, msi_280, msi_423, \
                         nsh0_318, nsh1_318, nsi_420, nsi_421, \
                         nsi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_17 * msi_280[k]
                   + f_3 * pc_y[k] * nsi_420[k];

        t_542[k] = f_3 * pc_z[k] * nsi_420[k];

        t_543[k] = f_17 * msi_423[k]
                   + f_10 * nsh0_318[k]
                   - f_11 * nsh1_318[k]
                   + f_3 * pc_x[k] * nsi_423[k];

        t_544[k] = f_3 * pc_z[k] * nsi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, msi_426, nsh0_315, nsh0_321, \
                         nsh1_315, nsh1_321, nsi_422, nsi_423, \
                         nsi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * nsh0_315[k]
                   - f_5 * nsh1_315[k]
                   + f_3 * pc_z[k] * nsi_422[k];

        t_546[k] = f_17 * msi_426[k]
                   + f_8 * nsh0_321[k]
                   - f_9 * nsh1_321[k]
                   + f_3 * pc_x[k] * nsi_426[k];

        t_547[k] = f_3 * pc_z[k] * nsi_423[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pc_x, pc_y, pc_z, msi_285, msi_430, \
                         nsh0_317, nsh0_325, nsh1_317, nsh1_325, nsi_425, nsi_426, \
                         nsi_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * msi_285[k]
                   + f_3 * pc_y[k] * nsi_425[k];

        t_549[k] = f_6 * nsh0_317[k]
                   - f_7 * nsh1_317[k]
                   + f_3 * pc_z[k] * nsi_425[k];

        t_550[k] = f_17 * msi_430[k]
                   + f_6 * nsh0_325[k]
                   - f_7 * nsh1_325[k]
                   + f_3 * pc_x[k] * nsi_430[k];

        t_551[k] = f_3 * pc_z[k] * nsi_426[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_y, pc_z, msi_289, nsh0_318, nsh0_320, \
                         nsh1_318, nsh1_320, nsi_427, nsi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_4 * nsh0_318[k]
                   - f_5 * nsh1_318[k]
                   + f_3 * pc_z[k] * nsi_427[k];

        t_553[k] = f_17 * msi_289[k]
                   + f_3 * pc_y[k] * nsi_429[k];

        t_554[k] = f_8 * nsh0_320[k]
                   - f_9 * nsh1_320[k]
                   + f_3 * pc_z[k] * nsi_429[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pc_x, pc_z, msi_435, nsh0_321, nsh0_330, \
                         nsh1_321, nsh1_330, nsi_430, nsi_431, \
                         nsi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_17 * msi_435[k]
                   + f_4 * nsh0_330[k]
                   - f_5 * nsh1_330[k]
                   + f_3 * pc_x[k] * nsi_435[k];

        t_556[k] = f_3 * pc_z[k] * nsi_430[k];

        t_557[k] = f_4 * nsh0_321[k]
                   - f_5 * nsh1_321[k]
                   + f_3 * pc_z[k] * nsi_431[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pc_x, pc_y, pc_z, msi_294, msi_441, \
                         nsh0_322, nsh0_324, nsh1_322, nsh1_324, nsi_432, nsi_434, \
                         nsi_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_6 * nsh0_322[k]
                   - f_7 * nsh1_322[k]
                   + f_3 * pc_z[k] * nsi_432[k];

        t_559[k] = f_17 * msi_294[k]
                   + f_3 * pc_y[k] * nsi_434[k];

        t_560[k] = f_10 * nsh0_324[k]
                   - f_11 * nsh1_324[k]
                   + f_3 * pc_z[k] * nsi_434[k];

        t_561[k] = f_17 * msi_441[k]
                   + f_3 * pc_x[k] * nsi_441[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pc_x, pc_z, msi_443, msi_444, \
                         msi_445, msi_446, nsi_435, nsi_443, nsi_444, nsi_445, \
                         nsi_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_3 * pc_z[k] * nsi_435[k];

        t_563[k] = f_17 * msi_443[k]
                   + f_3 * pc_x[k] * nsi_443[k];

        t_564[k] = f_17 * msi_444[k]
                   + f_3 * pc_x[k] * nsi_444[k];

        t_565[k] = f_17 * msi_445[k]
                   + f_3 * pc_x[k] * nsi_445[k];

        t_566[k] = f_17 * msi_446[k]
                   + f_3 * pc_x[k] * nsi_446[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, msi_301, msi_447, \
                         nsh0_330, nsh1_330, nsi_441, nsi_442, \
                         nsi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_17 * msi_447[k]
                   + f_3 * pc_x[k] * nsi_447[k];

        t_568[k] = f_17 * msi_301[k]
                   + f_1 * nsh0_330[k]
                   - f_2 * nsh1_330[k]
                   + f_3 * pc_y[k] * nsi_441[k];

        t_569[k] = f_3 * pc_z[k] * nsi_441[k];

        t_570[k] = f_4 * nsh0_330[k]
                   - f_5 * nsh1_330[k]
                   + f_3 * pc_z[k] * nsi_442[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_z, nsh0_331, nsh0_332, nsh0_333, nsh1_331, \
                         nsh1_332, nsh1_333, nsi_443, nsi_444, \
                         nsi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_6 * nsh0_331[k]
                   - f_7 * nsh1_331[k]
                   + f_3 * pc_z[k] * nsi_443[k];

        t_572[k] = f_8 * nsh0_332[k]
                   - f_9 * nsh1_332[k]
                   + f_3 * pc_z[k] * nsi_444[k];

        t_573[k] = f_10 * nsh0_333[k]
                   - f_11 * nsh1_333[k]
                   + f_3 * pc_z[k] * nsi_445[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_z, pc_y, pc_z, msk0_360, msi_307, \
                         msi_308, msk1_360, nsh0_335, nsh1_335, nsi_447, \
                         nsi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * msi_307[k]
                   + f_3 * pc_y[k] * nsi_447[k];

        t_575[k] = f_1 * nsh0_335[k]
                   - f_2 * nsh1_335[k]
                   + f_3 * pc_z[k] * nsi_447[k];

        t_576[k] = pa_z[k] * msk0_360[k]
                   - f_12 * pc_z[k] * msk1_360[k];

        t_577[k] = f_16 * msi_308[k]
                   + f_3 * pc_y[k] * nsi_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_z, pc_y, pc_z, msk0_363, msi_280, msi_310, \
                         msk1_363, nsi_448, nsi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * msi_280[k]
                   + f_3 * pc_z[k] * nsi_448[k];

        t_579[k] = pa_z[k] * msk0_363[k]
                   - f_12 * pc_z[k] * msk1_363[k];

        t_580[k] = f_16 * msi_310[k]
                   + f_3 * pc_y[k] * nsi_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_z, pc_x, pc_z, msk0_366, msi_283, msi_453, \
                         msk1_366, nsh0_341, nsh1_341, nsi_451, \
                         nsi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_17 * msi_453[k]
                   + f_10 * nsh0_341[k]
                   - f_11 * nsh1_341[k]
                   + f_3 * pc_x[k] * nsi_453[k];

        t_582[k] = pa_z[k] * msk0_366[k]
                   - f_12 * pc_z[k] * msk1_366[k];

        t_583[k] = f_13 * msi_283[k]
                   + f_3 * pc_z[k] * nsi_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pa_z, pc_x, pc_y, pc_z, msk0_370, msi_313, \
                         msi_457, msk1_370, nsh0_345, nsh1_345, nsi_453, \
                         nsi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * msi_313[k]
                   + f_3 * pc_y[k] * nsi_453[k];

        t_585[k] = f_17 * msi_457[k]
                   + f_8 * nsh0_345[k]
                   - f_9 * nsh1_345[k]
                   + f_3 * pc_x[k] * nsi_457[k];

        t_586[k] = pa_z[k] * msk0_370[k]
                   - f_12 * pc_z[k] * msk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pa_z, pc_y, pc_z, msk0_372, msi_286, msi_287, \
                         msi_317, msk1_372, nsi_454, nsi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * msi_286[k]
                   + f_3 * pc_z[k] * nsi_454[k];

        t_588[k] = pa_z[k] * msk0_372[k]
                   + f_14 * msi_287[k]
                   - f_12 * pc_z[k] * msk1_372[k];

        t_589[k] = f_16 * msi_317[k]
                   + f_3 * pc_y[k] * nsi_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pa_z, pc_x, pc_z, msk0_375, msi_290, msi_462, \
                         msk1_375, nsh0_350, nsh1_350, nsi_458, \
                         nsi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_17 * msi_462[k]
                   + f_6 * nsh0_350[k]
                   - f_7 * nsh1_350[k]
                   + f_3 * pc_x[k] * nsi_462[k];

        t_591[k] = pa_z[k] * msk0_375[k]
                   - f_12 * pc_z[k] * msk1_375[k];

        t_592[k] = f_13 * msi_290[k]
                   + f_3 * pc_z[k] * nsi_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pa_z, pc_y, pc_z, msk0_377, msk0_378, msi_291, \
                         msi_292, msi_322, msk1_377, msk1_378, \
                         nsi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pa_z[k] * msk0_377[k]
                   + f_14 * msi_291[k]
                   - f_12 * pc_z[k] * msk1_377[k];

        t_594[k] = pa_z[k] * msk0_378[k]
                   + f_15 * msi_292[k]
                   - f_12 * pc_z[k] * msk1_378[k];

        t_595[k] = f_16 * msi_322[k]
                   + f_3 * pc_y[k] * nsi_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, msi_468, msi_469, msi_470, msi_471, \
                         nsh0_356, nsh1_356, nsi_468, nsi_469, nsi_470, \
                         nsi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * msi_468[k]
                   + f_4 * nsh0_356[k]
                   - f_5 * nsh1_356[k]
                   + f_3 * pc_x[k] * nsi_468[k];

        t_597[k] = f_17 * msi_469[k]
                   + f_3 * pc_x[k] * nsi_469[k];

        t_598[k] = f_17 * msi_470[k]
                   + f_3 * pc_x[k] * nsi_470[k];

        t_599[k] = f_17 * msi_471[k]
                   + f_3 * pc_x[k] * nsi_471[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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

    const auto *msk0_388 = buffer.data(msk0 + 388);
    const auto *msk0_504 = buffer.data(msk0 + 504);
    const auto *msk0_507 = buffer.data(msk0 + 507);
    const auto *msk0_509 = buffer.data(msk0 + 509);
    const auto *msk0_510 = buffer.data(msk0 + 510);
    const auto *msk0_513 = buffer.data(msk0 + 513);
    const auto *msk0_514 = buffer.data(msk0 + 514);
    const auto *msk0_516 = buffer.data(msk0 + 516);
    const auto *msk0_518 = buffer.data(msk0 + 518);
    const auto *msk0_519 = buffer.data(msk0 + 519);
    const auto *msk0_521 = buffer.data(msk0 + 521);
    const auto *msk0_522 = buffer.data(msk0 + 522);
    const auto *msk0_524 = buffer.data(msk0 + 524);

    const auto *msi_301 = buffer.data(msi + 301);
    const auto *msi_307 = buffer.data(msi + 307);
    const auto *msi_308 = buffer.data(msi + 308);
    const auto *msi_311 = buffer.data(msi + 311);
    const auto *msi_314 = buffer.data(msi + 314);
    const auto *msi_318 = buffer.data(msi + 318);
    const auto *msi_329 = buffer.data(msi + 329);
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
    const auto *msi_350 = buffer.data(msi + 350);
    const auto *msi_357 = buffer.data(msi + 357);
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
    const auto *msi_373 = buffer.data(msi + 373);
    const auto *msi_374 = buffer.data(msi + 374);
    const auto *msi_378 = buffer.data(msi + 378);
    const auto *msi_385 = buffer.data(msi + 385);
    const auto *msi_387 = buffer.data(msi + 387);
    const auto *msi_388 = buffer.data(msi + 388);
    const auto *msi_389 = buffer.data(msi + 389);
    const auto *msi_390 = buffer.data(msi + 390);
    const auto *msi_391 = buffer.data(msi + 391);
    const auto *msi_392 = buffer.data(msi + 392);
    const auto *msi_393 = buffer.data(msi + 393);
    const auto *msi_394 = buffer.data(msi + 394);
    const auto *msi_395 = buffer.data(msi + 395);
    const auto *msi_397 = buffer.data(msi + 397);
    const auto *msi_398 = buffer.data(msi + 398);
    const auto *msi_400 = buffer.data(msi + 400);
    const auto *msi_401 = buffer.data(msi + 401);
    const auto *msi_402 = buffer.data(msi + 402);
    const auto *msi_404 = buffer.data(msi + 404);
    const auto *msi_405 = buffer.data(msi + 405);
    const auto *msi_406 = buffer.data(msi + 406);
    const auto *msi_472 = buffer.data(msi + 472);
    const auto *msi_473 = buffer.data(msi + 473);
    const auto *msi_474 = buffer.data(msi + 474);
    const auto *msi_475 = buffer.data(msi + 475);
    const auto *msi_476 = buffer.data(msi + 476);
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
    const auto *msi_553 = buffer.data(msi + 553);
    const auto *msi_554 = buffer.data(msi + 554);
    const auto *msi_555 = buffer.data(msi + 555);
    const auto *msi_556 = buffer.data(msi + 556);
    const auto *msi_557 = buffer.data(msi + 557);

    const auto *msk1_388 = buffer.data(msk1 + 388);
    const auto *msk1_504 = buffer.data(msk1 + 504);
    const auto *msk1_507 = buffer.data(msk1 + 507);
    const auto *msk1_509 = buffer.data(msk1 + 509);
    const auto *msk1_510 = buffer.data(msk1 + 510);
    const auto *msk1_513 = buffer.data(msk1 + 513);
    const auto *msk1_514 = buffer.data(msk1 + 514);
    const auto *msk1_516 = buffer.data(msk1 + 516);
    const auto *msk1_518 = buffer.data(msk1 + 518);
    const auto *msk1_519 = buffer.data(msk1 + 519);
    const auto *msk1_521 = buffer.data(msk1 + 521);
    const auto *msk1_522 = buffer.data(msk1 + 522);
    const auto *msk1_524 = buffer.data(msk1 + 524);

    const auto *nsh0_353 = buffer.data(nsh0 + 353);
    const auto *nsh0_354 = buffer.data(nsh0 + 354);
    const auto *nsh0_355 = buffer.data(nsh0 + 355);
    const auto *nsh0_356 = buffer.data(nsh0 + 356);
    const auto *nsh0_357 = buffer.data(nsh0 + 357);
    const auto *nsh0_360 = buffer.data(nsh0 + 360);
    const auto *nsh0_362 = buffer.data(nsh0 + 362);
    const auto *nsh0_363 = buffer.data(nsh0 + 363);
    const auto *nsh0_366 = buffer.data(nsh0 + 366);
    const auto *nsh0_367 = buffer.data(nsh0 + 367);
    const auto *nsh0_369 = buffer.data(nsh0 + 369);
    const auto *nsh0_371 = buffer.data(nsh0 + 371);
    const auto *nsh0_372 = buffer.data(nsh0 + 372);
    const auto *nsh0_374 = buffer.data(nsh0 + 374);
    const auto *nsh0_375 = buffer.data(nsh0 + 375);
    const auto *nsh0_376 = buffer.data(nsh0 + 376);
    const auto *nsh0_377 = buffer.data(nsh0 + 377);
    const auto *nsh0_378 = buffer.data(nsh0 + 378);
    const auto *nsh0_381 = buffer.data(nsh0 + 381);
    const auto *nsh0_383 = buffer.data(nsh0 + 383);
    const auto *nsh0_384 = buffer.data(nsh0 + 384);
    const auto *nsh0_387 = buffer.data(nsh0 + 387);
    const auto *nsh0_388 = buffer.data(nsh0 + 388);
    const auto *nsh0_390 = buffer.data(nsh0 + 390);
    const auto *nsh0_392 = buffer.data(nsh0 + 392);
    const auto *nsh0_393 = buffer.data(nsh0 + 393);
    const auto *nsh0_395 = buffer.data(nsh0 + 395);
    const auto *nsh0_396 = buffer.data(nsh0 + 396);
    const auto *nsh0_397 = buffer.data(nsh0 + 397);
    const auto *nsh0_398 = buffer.data(nsh0 + 398);

    const auto *nsh1_353 = buffer.data(nsh1 + 353);
    const auto *nsh1_354 = buffer.data(nsh1 + 354);
    const auto *nsh1_355 = buffer.data(nsh1 + 355);
    const auto *nsh1_356 = buffer.data(nsh1 + 356);
    const auto *nsh1_357 = buffer.data(nsh1 + 357);
    const auto *nsh1_360 = buffer.data(nsh1 + 360);
    const auto *nsh1_362 = buffer.data(nsh1 + 362);
    const auto *nsh1_363 = buffer.data(nsh1 + 363);
    const auto *nsh1_366 = buffer.data(nsh1 + 366);
    const auto *nsh1_367 = buffer.data(nsh1 + 367);
    const auto *nsh1_369 = buffer.data(nsh1 + 369);
    const auto *nsh1_371 = buffer.data(nsh1 + 371);
    const auto *nsh1_372 = buffer.data(nsh1 + 372);
    const auto *nsh1_374 = buffer.data(nsh1 + 374);
    const auto *nsh1_375 = buffer.data(nsh1 + 375);
    const auto *nsh1_376 = buffer.data(nsh1 + 376);
    const auto *nsh1_377 = buffer.data(nsh1 + 377);
    const auto *nsh1_378 = buffer.data(nsh1 + 378);
    const auto *nsh1_381 = buffer.data(nsh1 + 381);
    const auto *nsh1_383 = buffer.data(nsh1 + 383);
    const auto *nsh1_384 = buffer.data(nsh1 + 384);
    const auto *nsh1_387 = buffer.data(nsh1 + 387);
    const auto *nsh1_388 = buffer.data(nsh1 + 388);
    const auto *nsh1_390 = buffer.data(nsh1 + 390);
    const auto *nsh1_392 = buffer.data(nsh1 + 392);
    const auto *nsh1_393 = buffer.data(nsh1 + 393);
    const auto *nsh1_395 = buffer.data(nsh1 + 395);
    const auto *nsh1_396 = buffer.data(nsh1 + 396);
    const auto *nsh1_397 = buffer.data(nsh1 + 397);
    const auto *nsh1_398 = buffer.data(nsh1 + 398);

    const auto *nsi_469 = buffer.data(nsi + 469);
    const auto *nsi_471 = buffer.data(nsi + 471);
    const auto *nsi_472 = buffer.data(nsi + 472);
    const auto *nsi_473 = buffer.data(nsi + 473);
    const auto *nsi_474 = buffer.data(nsi + 474);
    const auto *nsi_475 = buffer.data(nsi + 475);
    const auto *nsi_476 = buffer.data(nsi + 476);
    const auto *nsi_478 = buffer.data(nsi + 478);
    const auto *nsi_479 = buffer.data(nsi + 479);
    const auto *nsi_481 = buffer.data(nsi + 481);
    const auto *nsi_482 = buffer.data(nsi + 482);
    const auto *nsi_485 = buffer.data(nsi + 485);
    const auto *nsi_486 = buffer.data(nsi + 486);
    const auto *nsi_488 = buffer.data(nsi + 488);
    const auto *nsi_490 = buffer.data(nsi + 490);
    const auto *nsi_491 = buffer.data(nsi + 491);
    const auto *nsi_493 = buffer.data(nsi + 493);
    const auto *nsi_494 = buffer.data(nsi + 494);
    const auto *nsi_496 = buffer.data(nsi + 496);
    const auto *nsi_497 = buffer.data(nsi + 497);
    const auto *nsi_498 = buffer.data(nsi + 498);
    const auto *nsi_499 = buffer.data(nsi + 499);
    const auto *nsi_500 = buffer.data(nsi + 500);
    const auto *nsi_501 = buffer.data(nsi + 501);
    const auto *nsi_502 = buffer.data(nsi + 502);
    const auto *nsi_503 = buffer.data(nsi + 503);
    const auto *nsi_504 = buffer.data(nsi + 504);
    const auto *nsi_506 = buffer.data(nsi + 506);
    const auto *nsi_507 = buffer.data(nsi + 507);
    const auto *nsi_509 = buffer.data(nsi + 509);
    const auto *nsi_510 = buffer.data(nsi + 510);
    const auto *nsi_513 = buffer.data(nsi + 513);
    const auto *nsi_514 = buffer.data(nsi + 514);
    const auto *nsi_516 = buffer.data(nsi + 516);
    const auto *nsi_518 = buffer.data(nsi + 518);
    const auto *nsi_519 = buffer.data(nsi + 519);
    const auto *nsi_521 = buffer.data(nsi + 521);
    const auto *nsi_522 = buffer.data(nsi + 522);
    const auto *nsi_524 = buffer.data(nsi + 524);
    const auto *nsi_525 = buffer.data(nsi + 525);
    const auto *nsi_526 = buffer.data(nsi + 526);
    const auto *nsi_527 = buffer.data(nsi + 527);
    const auto *nsi_528 = buffer.data(nsi + 528);
    const auto *nsi_529 = buffer.data(nsi + 529);
    const auto *nsi_530 = buffer.data(nsi + 530);
    const auto *nsi_531 = buffer.data(nsi + 531);
    const auto *nsi_532 = buffer.data(nsi + 532);
    const auto *nsi_534 = buffer.data(nsi + 534);
    const auto *nsi_535 = buffer.data(nsi + 535);
    const auto *nsi_537 = buffer.data(nsi + 537);
    const auto *nsi_538 = buffer.data(nsi + 538);
    const auto *nsi_541 = buffer.data(nsi + 541);
    const auto *nsi_542 = buffer.data(nsi + 542);
    const auto *nsi_546 = buffer.data(nsi + 546);
    const auto *nsi_553 = buffer.data(nsi + 553);
    const auto *nsi_554 = buffer.data(nsi + 554);
    const auto *nsi_555 = buffer.data(nsi + 555);
    const auto *nsi_556 = buffer.data(nsi + 556);
    const auto *nsi_557 = buffer.data(nsi + 557);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, msi_472, msi_473, msi_474, msi_475, \
                         nsi_472, nsi_473, nsi_474, nsi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_17 * msi_472[k]
                   + f_3 * pc_x[k] * nsi_472[k];

        t_601[k] = f_17 * msi_473[k]
                   + f_3 * pc_x[k] * nsi_473[k];

        t_602[k] = f_17 * msi_474[k]
                   + f_3 * pc_x[k] * nsi_474[k];

        t_603[k] = f_17 * msi_475[k]
                   + f_3 * pc_x[k] * nsi_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pa_z, pc_y, pc_z, msk0_388, msi_301, msi_331, \
                         msk1_388, nsh0_353, nsh1_353, nsi_469, \
                         nsi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_z[k] * msk0_388[k]
                   - f_12 * pc_z[k] * msk1_388[k];

        t_605[k] = f_13 * msi_301[k]
                   + f_3 * pc_z[k] * nsi_469[k];

        t_606[k] = f_16 * msi_331[k]
                   + f_10 * nsh0_353[k]
                   - f_11 * nsh1_353[k]
                   + f_3 * pc_y[k] * nsi_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, msi_332, msi_333, msi_334, nsh0_354, \
                         nsh0_355, nsh0_356, nsh1_354, nsh1_355, nsh1_356, nsi_472, nsi_473, \
                         nsi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * msi_332[k]
                   + f_8 * nsh0_354[k]
                   - f_9 * nsh1_354[k]
                   + f_3 * pc_y[k] * nsi_472[k];

        t_608[k] = f_16 * msi_333[k]
                   + f_6 * nsh0_355[k]
                   - f_7 * nsh1_355[k]
                   + f_3 * pc_y[k] * nsi_473[k];

        t_609[k] = f_16 * msi_334[k]
                   + f_4 * nsh0_356[k]
                   - f_5 * nsh1_356[k]
                   + f_3 * pc_y[k] * nsi_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, msi_307, msi_335, msi_476, \
                         nsh0_356, nsh0_357, nsh1_356, nsh1_357, nsi_475, \
                         nsi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * msi_335[k]
                   + f_3 * pc_y[k] * nsi_475[k];

        t_611[k] = f_13 * msi_307[k]
                   + f_1 * nsh0_356[k]
                   - f_2 * nsh1_356[k]
                   + f_3 * pc_z[k] * nsi_475[k];

        t_612[k] = f_17 * msi_476[k]
                   + f_1 * nsh0_357[k]
                   - f_2 * nsh1_357[k]
                   + f_3 * pc_x[k] * nsi_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, msi_308, msi_336, \
                         msi_338, msi_479, nsh0_360, nsh1_360, nsi_476, nsi_478, \
                         nsi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * msi_336[k]
                   + f_3 * pc_y[k] * nsi_476[k];

        t_614[k] = f_14 * msi_308[k]
                   + f_3 * pc_z[k] * nsi_476[k];

        t_615[k] = f_17 * msi_479[k]
                   + f_10 * nsh0_360[k]
                   - f_11 * nsh1_360[k]
                   + f_3 * pc_x[k] * nsi_479[k];

        t_616[k] = f_15 * msi_338[k]
                   + f_3 * pc_y[k] * nsi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, msi_311, msi_481, msi_482, nsh0_362, \
                         nsh0_363, nsh1_362, nsh1_363, nsi_479, nsi_481, \
                         nsi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_17 * msi_481[k]
                   + f_10 * nsh0_362[k]
                   - f_11 * nsh1_362[k]
                   + f_3 * pc_x[k] * nsi_481[k];

        t_618[k] = f_17 * msi_482[k]
                   + f_8 * nsh0_363[k]
                   - f_9 * nsh1_363[k]
                   + f_3 * pc_x[k] * nsi_482[k];

        t_619[k] = f_14 * msi_311[k]
                   + f_3 * pc_z[k] * nsi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, msi_341, msi_485, msi_486, nsh0_366, \
                         nsh0_367, nsh1_366, nsh1_367, nsi_481, nsi_485, \
                         nsi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * msi_341[k]
                   + f_3 * pc_y[k] * nsi_481[k];

        t_621[k] = f_17 * msi_485[k]
                   + f_8 * nsh0_366[k]
                   - f_9 * nsh1_366[k]
                   + f_3 * pc_x[k] * nsi_485[k];

        t_622[k] = f_17 * msi_486[k]
                   + f_6 * nsh0_367[k]
                   - f_7 * nsh1_367[k]
                   + f_3 * pc_x[k] * nsi_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, msi_314, msi_345, msi_488, \
                         nsh0_369, nsh1_369, nsi_482, nsi_485, \
                         nsi_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * msi_314[k]
                   + f_3 * pc_z[k] * nsi_482[k];

        t_624[k] = f_17 * msi_488[k]
                   + f_6 * nsh0_369[k]
                   - f_7 * nsh1_369[k]
                   + f_3 * pc_x[k] * nsi_488[k];

        t_625[k] = f_15 * msi_345[k]
                   + f_3 * pc_y[k] * nsi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, msi_318, msi_490, msi_491, nsh0_371, \
                         nsh0_372, nsh1_371, nsh1_372, nsi_486, nsi_490, \
                         nsi_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_17 * msi_490[k]
                   + f_6 * nsh0_371[k]
                   - f_7 * nsh1_371[k]
                   + f_3 * pc_x[k] * nsi_490[k];

        t_627[k] = f_17 * msi_491[k]
                   + f_4 * nsh0_372[k]
                   - f_5 * nsh1_372[k]
                   + f_3 * pc_x[k] * nsi_491[k];

        t_628[k] = f_14 * msi_318[k]
                   + f_3 * pc_z[k] * nsi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, msi_350, msi_493, msi_494, nsh0_374, \
                         nsh0_375, nsh1_374, nsh1_375, nsi_490, nsi_493, \
                         nsi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_17 * msi_493[k]
                   + f_4 * nsh0_374[k]
                   - f_5 * nsh1_374[k]
                   + f_3 * pc_x[k] * nsi_493[k];

        t_630[k] = f_17 * msi_494[k]
                   + f_4 * nsh0_375[k]
                   - f_5 * nsh1_375[k]
                   + f_3 * pc_x[k] * nsi_494[k];

        t_631[k] = f_15 * msi_350[k]
                   + f_3 * pc_y[k] * nsi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, msi_496, msi_497, msi_498, msi_499, \
                         nsh0_377, nsh1_377, nsi_496, nsi_497, nsi_498, \
                         nsi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_17 * msi_496[k]
                   + f_4 * nsh0_377[k]
                   - f_5 * nsh1_377[k]
                   + f_3 * pc_x[k] * nsi_496[k];

        t_633[k] = f_17 * msi_497[k]
                   + f_3 * pc_x[k] * nsi_497[k];

        t_634[k] = f_17 * msi_498[k]
                   + f_3 * pc_x[k] * nsi_498[k];

        t_635[k] = f_17 * msi_499[k]
                   + f_3 * pc_x[k] * nsi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, msi_500, msi_501, msi_502, msi_503, \
                         nsi_500, nsi_501, nsi_502, nsi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_17 * msi_500[k]
                   + f_3 * pc_x[k] * nsi_500[k];

        t_637[k] = f_17 * msi_501[k]
                   + f_3 * pc_x[k] * nsi_501[k];

        t_638[k] = f_17 * msi_502[k]
                   + f_3 * pc_x[k] * nsi_502[k];

        t_639[k] = f_17 * msi_503[k]
                   + f_3 * pc_x[k] * nsi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, msi_329, msi_357, msi_359, nsh0_372, \
                         nsh0_374, nsh1_372, nsh1_374, nsi_497, \
                         nsi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * msi_357[k]
                   + f_1 * nsh0_372[k]
                   - f_2 * nsh1_372[k]
                   + f_3 * pc_y[k] * nsi_497[k];

        t_641[k] = f_14 * msi_329[k]
                   + f_3 * pc_z[k] * nsi_497[k];

        t_642[k] = f_15 * msi_359[k]
                   + f_10 * nsh0_374[k]
                   - f_11 * nsh1_374[k]
                   + f_3 * pc_y[k] * nsi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, msi_360, msi_361, msi_362, nsh0_375, \
                         nsh0_376, nsh0_377, nsh1_375, nsh1_376, nsh1_377, nsi_500, nsi_501, \
                         nsi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * msi_360[k]
                   + f_8 * nsh0_375[k]
                   - f_9 * nsh1_375[k]
                   + f_3 * pc_y[k] * nsi_500[k];

        t_644[k] = f_15 * msi_361[k]
                   + f_6 * nsh0_376[k]
                   - f_7 * nsh1_376[k]
                   + f_3 * pc_y[k] * nsi_501[k];

        t_645[k] = f_15 * msi_362[k]
                   + f_4 * nsh0_377[k]
                   - f_5 * nsh1_377[k]
                   + f_3 * pc_y[k] * nsi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, msi_335, msi_363, msi_504, \
                         nsh0_377, nsh0_378, nsh1_377, nsh1_378, nsi_503, \
                         nsi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * msi_363[k]
                   + f_3 * pc_y[k] * nsi_503[k];

        t_647[k] = f_14 * msi_335[k]
                   + f_1 * nsh0_377[k]
                   - f_2 * nsh1_377[k]
                   + f_3 * pc_z[k] * nsi_503[k];

        t_648[k] = f_17 * msi_504[k]
                   + f_1 * nsh0_378[k]
                   - f_2 * nsh1_378[k]
                   + f_3 * pc_x[k] * nsi_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, msi_336, msi_364, \
                         msi_366, msi_507, nsh0_381, nsh1_381, nsi_504, nsi_506, \
                         nsi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * msi_364[k]
                   + f_3 * pc_y[k] * nsi_504[k];

        t_650[k] = f_15 * msi_336[k]
                   + f_3 * pc_z[k] * nsi_504[k];

        t_651[k] = f_17 * msi_507[k]
                   + f_10 * nsh0_381[k]
                   - f_11 * nsh1_381[k]
                   + f_3 * pc_x[k] * nsi_507[k];

        t_652[k] = f_14 * msi_366[k]
                   + f_3 * pc_y[k] * nsi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, msi_339, msi_509, msi_510, nsh0_383, \
                         nsh0_384, nsh1_383, nsh1_384, nsi_507, nsi_509, \
                         nsi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_17 * msi_509[k]
                   + f_10 * nsh0_383[k]
                   - f_11 * nsh1_383[k]
                   + f_3 * pc_x[k] * nsi_509[k];

        t_654[k] = f_17 * msi_510[k]
                   + f_8 * nsh0_384[k]
                   - f_9 * nsh1_384[k]
                   + f_3 * pc_x[k] * nsi_510[k];

        t_655[k] = f_15 * msi_339[k]
                   + f_3 * pc_z[k] * nsi_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, msi_369, msi_513, msi_514, nsh0_387, \
                         nsh0_388, nsh1_387, nsh1_388, nsi_509, nsi_513, \
                         nsi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * msi_369[k]
                   + f_3 * pc_y[k] * nsi_509[k];

        t_657[k] = f_17 * msi_513[k]
                   + f_8 * nsh0_387[k]
                   - f_9 * nsh1_387[k]
                   + f_3 * pc_x[k] * nsi_513[k];

        t_658[k] = f_17 * msi_514[k]
                   + f_6 * nsh0_388[k]
                   - f_7 * nsh1_388[k]
                   + f_3 * pc_x[k] * nsi_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, msi_342, msi_373, msi_516, \
                         nsh0_390, nsh1_390, nsi_510, nsi_513, \
                         nsi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * msi_342[k]
                   + f_3 * pc_z[k] * nsi_510[k];

        t_660[k] = f_17 * msi_516[k]
                   + f_6 * nsh0_390[k]
                   - f_7 * nsh1_390[k]
                   + f_3 * pc_x[k] * nsi_516[k];

        t_661[k] = f_14 * msi_373[k]
                   + f_3 * pc_y[k] * nsi_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, msi_346, msi_518, msi_519, nsh0_392, \
                         nsh0_393, nsh1_392, nsh1_393, nsi_514, nsi_518, \
                         nsi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_17 * msi_518[k]
                   + f_6 * nsh0_392[k]
                   - f_7 * nsh1_392[k]
                   + f_3 * pc_x[k] * nsi_518[k];

        t_663[k] = f_17 * msi_519[k]
                   + f_4 * nsh0_393[k]
                   - f_5 * nsh1_393[k]
                   + f_3 * pc_x[k] * nsi_519[k];

        t_664[k] = f_15 * msi_346[k]
                   + f_3 * pc_z[k] * nsi_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, msi_378, msi_521, msi_522, nsh0_395, \
                         nsh0_396, nsh1_395, nsh1_396, nsi_518, nsi_521, \
                         nsi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_17 * msi_521[k]
                   + f_4 * nsh0_395[k]
                   - f_5 * nsh1_395[k]
                   + f_3 * pc_x[k] * nsi_521[k];

        t_666[k] = f_17 * msi_522[k]
                   + f_4 * nsh0_396[k]
                   - f_5 * nsh1_396[k]
                   + f_3 * pc_x[k] * nsi_522[k];

        t_667[k] = f_14 * msi_378[k]
                   + f_3 * pc_y[k] * nsi_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, msi_524, msi_525, msi_526, msi_527, \
                         nsh0_398, nsh1_398, nsi_524, nsi_525, nsi_526, \
                         nsi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_17 * msi_524[k]
                   + f_4 * nsh0_398[k]
                   - f_5 * nsh1_398[k]
                   + f_3 * pc_x[k] * nsi_524[k];

        t_669[k] = f_17 * msi_525[k]
                   + f_3 * pc_x[k] * nsi_525[k];

        t_670[k] = f_17 * msi_526[k]
                   + f_3 * pc_x[k] * nsi_526[k];

        t_671[k] = f_17 * msi_527[k]
                   + f_3 * pc_x[k] * nsi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, msi_528, msi_529, msi_530, msi_531, \
                         nsi_528, nsi_529, nsi_530, nsi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_17 * msi_528[k]
                   + f_3 * pc_x[k] * nsi_528[k];

        t_673[k] = f_17 * msi_529[k]
                   + f_3 * pc_x[k] * nsi_529[k];

        t_674[k] = f_17 * msi_530[k]
                   + f_3 * pc_x[k] * nsi_530[k];

        t_675[k] = f_17 * msi_531[k]
                   + f_3 * pc_x[k] * nsi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, msi_357, msi_385, msi_387, nsh0_393, \
                         nsh0_395, nsh1_393, nsh1_395, nsi_525, \
                         nsi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * msi_385[k]
                   + f_1 * nsh0_393[k]
                   - f_2 * nsh1_393[k]
                   + f_3 * pc_y[k] * nsi_525[k];

        t_677[k] = f_15 * msi_357[k]
                   + f_3 * pc_z[k] * nsi_525[k];

        t_678[k] = f_14 * msi_387[k]
                   + f_10 * nsh0_395[k]
                   - f_11 * nsh1_395[k]
                   + f_3 * pc_y[k] * nsi_527[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, msi_388, msi_389, msi_390, nsh0_396, \
                         nsh0_397, nsh0_398, nsh1_396, nsh1_397, nsh1_398, nsi_528, nsi_529, \
                         nsi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * msi_388[k]
                   + f_8 * nsh0_396[k]
                   - f_9 * nsh1_396[k]
                   + f_3 * pc_y[k] * nsi_528[k];

        t_680[k] = f_14 * msi_389[k]
                   + f_6 * nsh0_397[k]
                   - f_7 * nsh1_397[k]
                   + f_3 * pc_y[k] * nsi_529[k];

        t_681[k] = f_14 * msi_390[k]
                   + f_4 * nsh0_398[k]
                   - f_5 * nsh1_398[k]
                   + f_3 * pc_y[k] * nsi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_y, pc_y, pc_z, msk0_504, msi_363, \
                         msi_391, msi_392, msk1_504, nsh0_398, nsh1_398, nsi_531, \
                         nsi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * msi_391[k]
                   + f_3 * pc_y[k] * nsi_531[k];

        t_683[k] = f_15 * msi_363[k]
                   + f_1 * nsh0_398[k]
                   - f_2 * nsh1_398[k]
                   + f_3 * pc_z[k] * nsi_531[k];

        t_684[k] = pa_y[k] * msk0_504[k]
                   - f_12 * pc_y[k] * msk1_504[k];

        t_685[k] = f_13 * msi_392[k]
                   + f_3 * pc_y[k] * nsi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_y, pc_y, pc_z, msk0_507, msk0_509, \
                         msi_364, msi_393, msi_394, msk1_507, msk1_509, nsi_532, \
                         nsi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * msi_364[k]
                   + f_3 * pc_z[k] * nsi_532[k];

        t_687[k] = pa_y[k] * msk0_507[k]
                   + f_14 * msi_393[k]
                   - f_12 * pc_y[k] * msk1_507[k];

        t_688[k] = f_13 * msi_394[k]
                   + f_3 * pc_y[k] * nsi_534[k];

        t_689[k] = pa_y[k] * msk0_509[k]
                   - f_12 * pc_y[k] * msk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_y, pc_y, pc_z, msk0_510, msk0_513, \
                         msi_367, msi_395, msi_397, msk1_510, msk1_513, nsi_535, \
                         nsi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_y[k] * msk0_510[k]
                   + f_15 * msi_395[k]
                   - f_12 * pc_y[k] * msk1_510[k];

        t_691[k] = f_16 * msi_367[k]
                   + f_3 * pc_z[k] * nsi_535[k];

        t_692[k] = f_13 * msi_397[k]
                   + f_3 * pc_y[k] * nsi_537[k];

        t_693[k] = pa_y[k] * msk0_513[k]
                   - f_12 * pc_y[k] * msk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pa_y, pc_y, pc_z, msk0_514, msk0_516, msi_370, \
                         msi_398, msi_400, msk1_514, msk1_516, \
                         nsi_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pa_y[k] * msk0_514[k]
                   + f_16 * msi_398[k]
                   - f_12 * pc_y[k] * msk1_514[k];

        t_695[k] = f_16 * msi_370[k]
                   + f_3 * pc_z[k] * nsi_538[k];

        t_696[k] = pa_y[k] * msk0_516[k]
                   + f_14 * msi_400[k]
                   - f_12 * pc_y[k] * msk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pa_y, pc_y, pc_z, msk0_518, msk0_519, \
                         msi_374, msi_401, msi_402, msk1_518, msk1_519, nsi_541, \
                         nsi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * msi_401[k]
                   + f_3 * pc_y[k] * nsi_541[k];

        t_698[k] = pa_y[k] * msk0_518[k]
                   - f_12 * pc_y[k] * msk1_518[k];

        t_699[k] = pa_y[k] * msk0_519[k]
                   + f_17 * msi_402[k]
                   - f_12 * pc_y[k] * msk1_519[k];

        t_700[k] = f_16 * msi_374[k]
                   + f_3 * pc_z[k] * nsi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pa_y, pc_y, msk0_521, msk0_522, msk0_524, \
                         msi_404, msi_405, msi_406, msk1_521, msk1_522, msk1_524, \
                         nsi_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pa_y[k] * msk0_521[k]
                   + f_15 * msi_404[k]
                   - f_12 * pc_y[k] * msk1_521[k];

        t_702[k] = pa_y[k] * msk0_522[k]
                   + f_14 * msi_405[k]
                   - f_12 * pc_y[k] * msk1_522[k];

        t_703[k] = f_13 * msi_406[k]
                   + f_3 * pc_y[k] * nsi_546[k];

        t_704[k] = pa_y[k] * msk0_524[k]
                   - f_12 * pc_y[k] * msk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, msi_553, msi_554, msi_555, \
                         msi_556, msi_557, nsi_553, nsi_554, nsi_555, nsi_556, \
                         nsi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_17 * msi_553[k]
                   + f_3 * pc_x[k] * nsi_553[k];

        t_706[k] = f_17 * msi_554[k]
                   + f_3 * pc_x[k] * nsi_554[k];

        t_707[k] = f_17 * msi_555[k]
                   + f_3 * pc_x[k] * nsi_555[k];

        t_708[k] = f_17 * msi_556[k]
                   + f_3 * pc_x[k] * nsi_556[k];

        t_709[k] = f_17 * msi_557[k]
                   + f_3 * pc_x[k] * nsi_557[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_539 = buffer.data(msk0 + 539);
    const auto *msk0_540 = buffer.data(msk0 + 540);
    const auto *msk0_543 = buffer.data(msk0 + 543);
    const auto *msk0_546 = buffer.data(msk0 + 546);
    const auto *msk0_550 = buffer.data(msk0 + 550);
    const auto *msk0_552 = buffer.data(msk0 + 552);
    const auto *msk0_555 = buffer.data(msk0 + 555);
    const auto *msk0_557 = buffer.data(msk0 + 557);
    const auto *msk0_558 = buffer.data(msk0 + 558);
    const auto *msk0_568 = buffer.data(msk0 + 568);

    const auto *msi_385 = buffer.data(msi + 385);
    const auto *msi_392 = buffer.data(msi + 392);
    const auto *msi_413 = buffer.data(msi + 413);
    const auto *msi_415 = buffer.data(msi + 415);
    const auto *msi_416 = buffer.data(msi + 416);
    const auto *msi_417 = buffer.data(msi + 417);
    const auto *msi_418 = buffer.data(msi + 418);
    const auto *msi_419 = buffer.data(msi + 419);
    const auto *msi_420 = buffer.data(msi + 420);
    const auto *msi_423 = buffer.data(msi + 423);
    const auto *msi_425 = buffer.data(msi + 425);
    const auto *msi_426 = buffer.data(msi + 426);
    const auto *msi_427 = buffer.data(msi + 427);
    const auto *msi_429 = buffer.data(msi + 429);
    const auto *msi_430 = buffer.data(msi + 430);
    const auto *msi_431 = buffer.data(msi + 431);
    const auto *msi_432 = buffer.data(msi + 432);
    const auto *msi_434 = buffer.data(msi + 434);
    const auto *msi_441 = buffer.data(msi + 441);
    const auto *msi_447 = buffer.data(msi + 447);
    const auto *msi_448 = buffer.data(msi + 448);
    const auto *msi_450 = buffer.data(msi + 450);
    const auto *msi_453 = buffer.data(msi + 453);
    const auto *msi_457 = buffer.data(msi + 457);
    const auto *msi_462 = buffer.data(msi + 462);
    const auto *msi_471 = buffer.data(msi + 471);
    const auto *msi_472 = buffer.data(msi + 472);
    const auto *msi_473 = buffer.data(msi + 473);
    const auto *msi_474 = buffer.data(msi + 474);
    const auto *msi_558 = buffer.data(msi + 558);
    const auto *msi_559 = buffer.data(msi + 559);
    const auto *msi_560 = buffer.data(msi + 560);
    const auto *msi_565 = buffer.data(msi + 565);
    const auto *msi_569 = buffer.data(msi + 569);
    const auto *msi_574 = buffer.data(msi + 574);
    const auto *msi_580 = buffer.data(msi + 580);
    const auto *msi_581 = buffer.data(msi + 581);
    const auto *msi_582 = buffer.data(msi + 582);
    const auto *msi_583 = buffer.data(msi + 583);
    const auto *msi_584 = buffer.data(msi + 584);
    const auto *msi_585 = buffer.data(msi + 585);
    const auto *msi_587 = buffer.data(msi + 587);
    const auto *msi_588 = buffer.data(msi + 588);
    const auto *msi_591 = buffer.data(msi + 591);
    const auto *msi_594 = buffer.data(msi + 594);
    const auto *msi_598 = buffer.data(msi + 598);
    const auto *msi_603 = buffer.data(msi + 603);
    const auto *msi_609 = buffer.data(msi + 609);
    const auto *msi_611 = buffer.data(msi + 611);
    const auto *msi_612 = buffer.data(msi + 612);
    const auto *msi_613 = buffer.data(msi + 613);
    const auto *msi_614 = buffer.data(msi + 614);
    const auto *msi_615 = buffer.data(msi + 615);
    const auto *msi_621 = buffer.data(msi + 621);
    const auto *msi_625 = buffer.data(msi + 625);
    const auto *msi_630 = buffer.data(msi + 630);
    const auto *msi_636 = buffer.data(msi + 636);
    const auto *msi_637 = buffer.data(msi + 637);
    const auto *msi_638 = buffer.data(msi + 638);
    const auto *msi_639 = buffer.data(msi + 639);
    const auto *msi_640 = buffer.data(msi + 640);
    const auto *msi_641 = buffer.data(msi + 641);
    const auto *msi_642 = buffer.data(msi + 642);
    const auto *msi_643 = buffer.data(msi + 643);

    const auto *msk1_539 = buffer.data(msk1 + 539);
    const auto *msk1_540 = buffer.data(msk1 + 540);
    const auto *msk1_543 = buffer.data(msk1 + 543);
    const auto *msk1_546 = buffer.data(msk1 + 546);
    const auto *msk1_550 = buffer.data(msk1 + 550);
    const auto *msk1_552 = buffer.data(msk1 + 552);
    const auto *msk1_555 = buffer.data(msk1 + 555);
    const auto *msk1_557 = buffer.data(msk1 + 557);
    const auto *msk1_558 = buffer.data(msk1 + 558);
    const auto *msk1_568 = buffer.data(msk1 + 568);

    const auto *nsh0_414 = buffer.data(nsh0 + 414);
    const auto *nsh0_416 = buffer.data(nsh0 + 416);
    const auto *nsh0_417 = buffer.data(nsh0 + 417);
    const auto *nsh0_418 = buffer.data(nsh0 + 418);
    const auto *nsh0_419 = buffer.data(nsh0 + 419);
    const auto *nsh0_420 = buffer.data(nsh0 + 420);
    const auto *nsh0_421 = buffer.data(nsh0 + 421);
    const auto *nsh0_422 = buffer.data(nsh0 + 422);
    const auto *nsh0_423 = buffer.data(nsh0 + 423);
    const auto *nsh0_424 = buffer.data(nsh0 + 424);
    const auto *nsh0_425 = buffer.data(nsh0 + 425);
    const auto *nsh0_426 = buffer.data(nsh0 + 426);
    const auto *nsh0_427 = buffer.data(nsh0 + 427);
    const auto *nsh0_428 = buffer.data(nsh0 + 428);
    const auto *nsh0_429 = buffer.data(nsh0 + 429);
    const auto *nsh0_434 = buffer.data(nsh0 + 434);
    const auto *nsh0_435 = buffer.data(nsh0 + 435);
    const auto *nsh0_436 = buffer.data(nsh0 + 436);
    const auto *nsh0_437 = buffer.data(nsh0 + 437);
    const auto *nsh0_438 = buffer.data(nsh0 + 438);
    const auto *nsh0_439 = buffer.data(nsh0 + 439);
    const auto *nsh0_440 = buffer.data(nsh0 + 440);
    const auto *nsh0_441 = buffer.data(nsh0 + 441);
    const auto *nsh0_443 = buffer.data(nsh0 + 443);
    const auto *nsh0_444 = buffer.data(nsh0 + 444);
    const auto *nsh0_446 = buffer.data(nsh0 + 446);
    const auto *nsh0_447 = buffer.data(nsh0 + 447);
    const auto *nsh0_448 = buffer.data(nsh0 + 448);
    const auto *nsh0_450 = buffer.data(nsh0 + 450);
    const auto *nsh0_451 = buffer.data(nsh0 + 451);
    const auto *nsh0_456 = buffer.data(nsh0 + 456);
    const auto *nsh0_457 = buffer.data(nsh0 + 457);
    const auto *nsh0_458 = buffer.data(nsh0 + 458);
    const auto *nsh0_459 = buffer.data(nsh0 + 459);
    const auto *nsh0_461 = buffer.data(nsh0 + 461);
    const auto *nsh0_467 = buffer.data(nsh0 + 467);
    const auto *nsh0_471 = buffer.data(nsh0 + 471);
    const auto *nsh0_476 = buffer.data(nsh0 + 476);
    const auto *nsh0_479 = buffer.data(nsh0 + 479);
    const auto *nsh0_480 = buffer.data(nsh0 + 480);
    const auto *nsh0_481 = buffer.data(nsh0 + 481);
    const auto *nsh0_482 = buffer.data(nsh0 + 482);

    const auto *nsh1_414 = buffer.data(nsh1 + 414);
    const auto *nsh1_416 = buffer.data(nsh1 + 416);
    const auto *nsh1_417 = buffer.data(nsh1 + 417);
    const auto *nsh1_418 = buffer.data(nsh1 + 418);
    const auto *nsh1_419 = buffer.data(nsh1 + 419);
    const auto *nsh1_420 = buffer.data(nsh1 + 420);
    const auto *nsh1_421 = buffer.data(nsh1 + 421);
    const auto *nsh1_422 = buffer.data(nsh1 + 422);
    const auto *nsh1_423 = buffer.data(nsh1 + 423);
    const auto *nsh1_424 = buffer.data(nsh1 + 424);
    const auto *nsh1_425 = buffer.data(nsh1 + 425);
    const auto *nsh1_426 = buffer.data(nsh1 + 426);
    const auto *nsh1_427 = buffer.data(nsh1 + 427);
    const auto *nsh1_428 = buffer.data(nsh1 + 428);
    const auto *nsh1_429 = buffer.data(nsh1 + 429);
    const auto *nsh1_434 = buffer.data(nsh1 + 434);
    const auto *nsh1_435 = buffer.data(nsh1 + 435);
    const auto *nsh1_436 = buffer.data(nsh1 + 436);
    const auto *nsh1_437 = buffer.data(nsh1 + 437);
    const auto *nsh1_438 = buffer.data(nsh1 + 438);
    const auto *nsh1_439 = buffer.data(nsh1 + 439);
    const auto *nsh1_440 = buffer.data(nsh1 + 440);
    const auto *nsh1_441 = buffer.data(nsh1 + 441);
    const auto *nsh1_443 = buffer.data(nsh1 + 443);
    const auto *nsh1_444 = buffer.data(nsh1 + 444);
    const auto *nsh1_446 = buffer.data(nsh1 + 446);
    const auto *nsh1_447 = buffer.data(nsh1 + 447);
    const auto *nsh1_448 = buffer.data(nsh1 + 448);
    const auto *nsh1_450 = buffer.data(nsh1 + 450);
    const auto *nsh1_451 = buffer.data(nsh1 + 451);
    const auto *nsh1_456 = buffer.data(nsh1 + 456);
    const auto *nsh1_457 = buffer.data(nsh1 + 457);
    const auto *nsh1_458 = buffer.data(nsh1 + 458);
    const auto *nsh1_459 = buffer.data(nsh1 + 459);
    const auto *nsh1_461 = buffer.data(nsh1 + 461);
    const auto *nsh1_467 = buffer.data(nsh1 + 467);
    const auto *nsh1_471 = buffer.data(nsh1 + 471);
    const auto *nsh1_476 = buffer.data(nsh1 + 476);
    const auto *nsh1_479 = buffer.data(nsh1 + 479);
    const auto *nsh1_480 = buffer.data(nsh1 + 480);
    const auto *nsh1_481 = buffer.data(nsh1 + 481);
    const auto *nsh1_482 = buffer.data(nsh1 + 482);

    const auto *nsi_553 = buffer.data(nsi + 553);
    const auto *nsi_555 = buffer.data(nsi + 555);
    const auto *nsi_556 = buffer.data(nsi + 556);
    const auto *nsi_557 = buffer.data(nsi + 557);
    const auto *nsi_558 = buffer.data(nsi + 558);
    const auto *nsi_559 = buffer.data(nsi + 559);
    const auto *nsi_560 = buffer.data(nsi + 560);
    const auto *nsi_561 = buffer.data(nsi + 561);
    const auto *nsi_562 = buffer.data(nsi + 562);
    const auto *nsi_563 = buffer.data(nsi + 563);
    const auto *nsi_564 = buffer.data(nsi + 564);
    const auto *nsi_565 = buffer.data(nsi + 565);
    const auto *nsi_566 = buffer.data(nsi + 566);
    const auto *nsi_567 = buffer.data(nsi + 567);
    const auto *nsi_568 = buffer.data(nsi + 568);
    const auto *nsi_569 = buffer.data(nsi + 569);
    const auto *nsi_570 = buffer.data(nsi + 570);
    const auto *nsi_571 = buffer.data(nsi + 571);
    const auto *nsi_572 = buffer.data(nsi + 572);
    const auto *nsi_573 = buffer.data(nsi + 573);
    const auto *nsi_574 = buffer.data(nsi + 574);
    const auto *nsi_580 = buffer.data(nsi + 580);
    const auto *nsi_581 = buffer.data(nsi + 581);
    const auto *nsi_582 = buffer.data(nsi + 582);
    const auto *nsi_583 = buffer.data(nsi + 583);
    const auto *nsi_584 = buffer.data(nsi + 584);
    const auto *nsi_585 = buffer.data(nsi + 585);
    const auto *nsi_586 = buffer.data(nsi + 586);
    const auto *nsi_587 = buffer.data(nsi + 587);
    const auto *nsi_588 = buffer.data(nsi + 588);
    const auto *nsi_589 = buffer.data(nsi + 589);
    const auto *nsi_590 = buffer.data(nsi + 590);
    const auto *nsi_591 = buffer.data(nsi + 591);
    const auto *nsi_593 = buffer.data(nsi + 593);
    const auto *nsi_594 = buffer.data(nsi + 594);
    const auto *nsi_595 = buffer.data(nsi + 595);
    const auto *nsi_597 = buffer.data(nsi + 597);
    const auto *nsi_598 = buffer.data(nsi + 598);
    const auto *nsi_599 = buffer.data(nsi + 599);
    const auto *nsi_600 = buffer.data(nsi + 600);
    const auto *nsi_602 = buffer.data(nsi + 602);
    const auto *nsi_603 = buffer.data(nsi + 603);
    const auto *nsi_609 = buffer.data(nsi + 609);
    const auto *nsi_610 = buffer.data(nsi + 610);
    const auto *nsi_611 = buffer.data(nsi + 611);
    const auto *nsi_612 = buffer.data(nsi + 612);
    const auto *nsi_613 = buffer.data(nsi + 613);
    const auto *nsi_614 = buffer.data(nsi + 614);
    const auto *nsi_615 = buffer.data(nsi + 615);
    const auto *nsi_616 = buffer.data(nsi + 616);
    const auto *nsi_618 = buffer.data(nsi + 618);
    const auto *nsi_619 = buffer.data(nsi + 619);
    const auto *nsi_621 = buffer.data(nsi + 621);
    const auto *nsi_622 = buffer.data(nsi + 622);
    const auto *nsi_625 = buffer.data(nsi + 625);
    const auto *nsi_626 = buffer.data(nsi + 626);
    const auto *nsi_630 = buffer.data(nsi + 630);
    const auto *nsi_636 = buffer.data(nsi + 636);
    const auto *nsi_637 = buffer.data(nsi + 637);
    const auto *nsi_638 = buffer.data(nsi + 638);
    const auto *nsi_639 = buffer.data(nsi + 639);
    const auto *nsi_640 = buffer.data(nsi + 640);
    const auto *nsi_641 = buffer.data(nsi + 641);
    const auto *nsi_642 = buffer.data(nsi + 642);
    const auto *nsi_643 = buffer.data(nsi + 643);

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, msi_385, msi_413, \
                         msi_558, msi_559, nsh0_414, nsh1_414, nsi_553, nsi_558, \
                         nsi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_17 * msi_558[k]
                   + f_3 * pc_x[k] * nsi_558[k];

        t_711[k] = f_17 * msi_559[k]
                   + f_3 * pc_x[k] * nsi_559[k];

        t_712[k] = f_13 * msi_413[k]
                   + f_1 * nsh0_414[k]
                   - f_2 * nsh1_414[k]
                   + f_3 * pc_y[k] * nsi_553[k];

        t_713[k] = f_16 * msi_385[k]
                   + f_3 * pc_z[k] * nsi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, msi_415, msi_416, msi_417, nsh0_416, \
                         nsh0_417, nsh0_418, nsh1_416, nsh1_417, nsh1_418, nsi_555, nsi_556, \
                         nsi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * msi_415[k]
                   + f_10 * nsh0_416[k]
                   - f_11 * nsh1_416[k]
                   + f_3 * pc_y[k] * nsi_555[k];

        t_715[k] = f_13 * msi_416[k]
                   + f_8 * nsh0_417[k]
                   - f_9 * nsh1_417[k]
                   + f_3 * pc_y[k] * nsi_556[k];

        t_716[k] = f_13 * msi_417[k]
                   + f_6 * nsh0_418[k]
                   - f_7 * nsh1_418[k]
                   + f_3 * pc_y[k] * nsi_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_y, pc_y, msk0_539, msi_418, msi_419, \
                         msk1_539, nsh0_419, nsh1_419, nsi_558, \
                         nsi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * msi_418[k]
                   + f_4 * nsh0_419[k]
                   - f_5 * nsh1_419[k]
                   + f_3 * pc_y[k] * nsi_558[k];

        t_718[k] = f_13 * msi_419[k]
                   + f_3 * pc_y[k] * nsi_559[k];

        t_719[k] = pa_y[k] * msk0_539[k]
                   - f_12 * pc_y[k] * msk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pc_x, pc_y, pc_z, msi_392, \
                         msi_560, nsh0_420, nsh1_420, nsi_560, nsi_561, \
                         nsi_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_17 * msi_560[k]
                   + f_1 * nsh0_420[k]
                   - f_2 * nsh1_420[k]
                   + f_3 * pc_x[k] * nsi_560[k];

        t_721[k] = f_3 * pc_y[k] * nsi_560[k];

        t_722[k] = f_17 * msi_392[k]
                   + f_3 * pc_z[k] * nsi_560[k];

        t_723[k] = f_4 * nsh0_420[k]
                   - f_5 * nsh1_420[k]
                   + f_3 * pc_y[k] * nsi_561[k];

        t_724[k] = f_3 * pc_y[k] * nsi_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, pc_y, msi_565, nsh0_421, nsh0_422, \
                         nsh0_425, nsh1_421, nsh1_422, nsh1_425, nsi_563, nsi_564, \
                         nsi_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_17 * msi_565[k]
                   + f_10 * nsh0_425[k]
                   - f_11 * nsh1_425[k]
                   + f_3 * pc_x[k] * nsi_565[k];

        t_726[k] = f_6 * nsh0_421[k]
                   - f_7 * nsh1_421[k]
                   + f_3 * pc_y[k] * nsi_563[k];

        t_727[k] = f_4 * nsh0_422[k]
                   - f_5 * nsh1_422[k]
                   + f_3 * pc_y[k] * nsi_564[k];

        t_728[k] = f_3 * pc_y[k] * nsi_565[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, msi_569, nsh0_423, nsh0_424, \
                         nsh0_429, nsh1_423, nsh1_424, nsh1_429, nsi_566, nsi_567, \
                         nsi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_17 * msi_569[k]
                   + f_8 * nsh0_429[k]
                   - f_9 * nsh1_429[k]
                   + f_3 * pc_x[k] * nsi_569[k];

        t_730[k] = f_8 * nsh0_423[k]
                   - f_9 * nsh1_423[k]
                   + f_3 * pc_y[k] * nsi_566[k];

        t_731[k] = f_6 * nsh0_424[k]
                   - f_7 * nsh1_424[k]
                   + f_3 * pc_y[k] * nsi_567[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_x, pc_y, msi_574, nsh0_425, nsh0_434, \
                         nsh1_425, nsh1_434, nsi_568, nsi_569, \
                         nsi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_4 * nsh0_425[k]
                   - f_5 * nsh1_425[k]
                   + f_3 * pc_y[k] * nsi_568[k];

        t_733[k] = f_3 * pc_y[k] * nsi_569[k];

        t_734[k] = f_17 * msi_574[k]
                   + f_6 * nsh0_434[k]
                   - f_7 * nsh1_434[k]
                   + f_3 * pc_x[k] * nsi_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pc_y, nsh0_426, nsh0_427, nsh0_428, nsh1_426, \
                         nsh1_427, nsh1_428, nsi_570, nsi_571, \
                         nsi_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_10 * nsh0_426[k]
                   - f_11 * nsh1_426[k]
                   + f_3 * pc_y[k] * nsi_570[k];

        t_736[k] = f_8 * nsh0_427[k]
                   - f_9 * nsh1_427[k]
                   + f_3 * pc_y[k] * nsi_571[k];

        t_737[k] = f_6 * nsh0_428[k]
                   - f_7 * nsh1_428[k]
                   + f_3 * pc_y[k] * nsi_572[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, pc_x, pc_y, msi_580, msi_581, nsh0_429, \
                         nsh0_440, nsh1_429, nsh1_440, nsi_573, nsi_574, nsi_580, \
                         nsi_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_4 * nsh0_429[k]
                   - f_5 * nsh1_429[k]
                   + f_3 * pc_y[k] * nsi_573[k];

        t_739[k] = f_3 * pc_y[k] * nsi_574[k];

        t_740[k] = f_17 * msi_580[k]
                   + f_4 * nsh0_440[k]
                   - f_5 * nsh1_440[k]
                   + f_3 * pc_x[k] * nsi_580[k];

        t_741[k] = f_17 * msi_581[k]
                   + f_3 * pc_x[k] * nsi_581[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, pc_x, pc_y, msi_582, msi_583, \
                         msi_584, msi_585, nsi_580, nsi_582, nsi_583, nsi_584, \
                         nsi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_17 * msi_582[k]
                   + f_3 * pc_x[k] * nsi_582[k];

        t_743[k] = f_17 * msi_583[k]
                   + f_3 * pc_x[k] * nsi_583[k];

        t_744[k] = f_17 * msi_584[k]
                   + f_3 * pc_x[k] * nsi_584[k];

        t_745[k] = f_17 * msi_585[k]
                   + f_3 * pc_x[k] * nsi_585[k];

        t_746[k] = f_3 * pc_y[k] * nsi_580[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_x, pc_y, msi_587, nsh0_435, nsh0_436, \
                         nsh1_435, nsh1_436, nsi_581, nsi_582, \
                         nsi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_17 * msi_587[k]
                   + f_3 * pc_x[k] * nsi_587[k];

        t_748[k] = f_1 * nsh0_435[k]
                   - f_2 * nsh1_435[k]
                   + f_3 * pc_y[k] * nsi_581[k];

        t_749[k] = f_19 * nsh0_436[k]
                   - f_20 * nsh1_436[k]
                   + f_3 * pc_y[k] * nsi_582[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, nsh0_437, nsh0_438, nsh0_439, nsh1_437, \
                         nsh1_438, nsh1_439, nsi_583, nsi_584, \
                         nsi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_10 * nsh0_437[k]
                   - f_11 * nsh1_437[k]
                   + f_3 * pc_y[k] * nsi_583[k];

        t_751[k] = f_8 * nsh0_438[k]
                   - f_9 * nsh1_438[k]
                   + f_3 * pc_y[k] * nsi_584[k];

        t_752[k] = f_6 * nsh0_439[k]
                   - f_7 * nsh1_439[k]
                   + f_3 * pc_y[k] * nsi_585[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pc_x, pc_y, pc_z, msi_419, msi_588, \
                         nsh0_440, nsh0_441, nsh1_440, nsh1_441, nsi_586, nsi_587, \
                         nsi_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * nsh0_440[k]
                   - f_5 * nsh1_440[k]
                   + f_3 * pc_y[k] * nsi_586[k];

        t_754[k] = f_3 * pc_y[k] * nsi_587[k];

        t_755[k] = f_17 * msi_419[k]
                   + f_1 * nsh0_440[k]
                   - f_2 * nsh1_440[k]
                   + f_3 * pc_z[k] * nsi_587[k];

        t_756[k] = f_16 * msi_588[k]
                   + f_1 * nsh0_441[k]
                   - f_2 * nsh1_441[k]
                   + f_3 * pc_x[k] * nsi_588[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pc_x, pc_y, pc_z, msi_420, msi_591, \
                         nsh0_444, nsh1_444, nsi_588, nsi_589, \
                         nsi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_23 * msi_420[k]
                   + f_3 * pc_y[k] * nsi_588[k];

        t_758[k] = f_3 * pc_z[k] * nsi_588[k];

        t_759[k] = f_16 * msi_591[k]
                   + f_10 * nsh0_444[k]
                   - f_11 * nsh1_444[k]
                   + f_3 * pc_x[k] * nsi_591[k];

        t_760[k] = f_3 * pc_z[k] * nsi_589[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_x, pc_z, msi_594, nsh0_441, nsh0_447, \
                         nsh1_441, nsh1_447, nsi_590, nsi_591, \
                         nsi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_4 * nsh0_441[k]
                   - f_5 * nsh1_441[k]
                   + f_3 * pc_z[k] * nsi_590[k];

        t_762[k] = f_16 * msi_594[k]
                   + f_8 * nsh0_447[k]
                   - f_9 * nsh1_447[k]
                   + f_3 * pc_x[k] * nsi_594[k];

        t_763[k] = f_3 * pc_z[k] * nsi_591[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, msi_425, msi_598, \
                         nsh0_443, nsh0_451, nsh1_443, nsh1_451, nsi_593, nsi_594, \
                         nsi_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_23 * msi_425[k]
                   + f_3 * pc_y[k] * nsi_593[k];

        t_765[k] = f_6 * nsh0_443[k]
                   - f_7 * nsh1_443[k]
                   + f_3 * pc_z[k] * nsi_593[k];

        t_766[k] = f_16 * msi_598[k]
                   + f_6 * nsh0_451[k]
                   - f_7 * nsh1_451[k]
                   + f_3 * pc_x[k] * nsi_598[k];

        t_767[k] = f_3 * pc_z[k] * nsi_594[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_y, pc_z, msi_429, nsh0_444, nsh0_446, \
                         nsh1_444, nsh1_446, nsi_595, nsi_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * nsh0_444[k]
                   - f_5 * nsh1_444[k]
                   + f_3 * pc_z[k] * nsi_595[k];

        t_769[k] = f_23 * msi_429[k]
                   + f_3 * pc_y[k] * nsi_597[k];

        t_770[k] = f_8 * nsh0_446[k]
                   - f_9 * nsh1_446[k]
                   + f_3 * pc_z[k] * nsi_597[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, pc_x, pc_z, msi_603, nsh0_447, nsh0_456, \
                         nsh1_447, nsh1_456, nsi_598, nsi_599, \
                         nsi_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_16 * msi_603[k]
                   + f_4 * nsh0_456[k]
                   - f_5 * nsh1_456[k]
                   + f_3 * pc_x[k] * nsi_603[k];

        t_772[k] = f_3 * pc_z[k] * nsi_598[k];

        t_773[k] = f_4 * nsh0_447[k]
                   - f_5 * nsh1_447[k]
                   + f_3 * pc_z[k] * nsi_599[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pc_x, pc_y, pc_z, msi_434, msi_609, \
                         nsh0_448, nsh0_450, nsh1_448, nsh1_450, nsi_600, nsi_602, \
                         nsi_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_6 * nsh0_448[k]
                   - f_7 * nsh1_448[k]
                   + f_3 * pc_z[k] * nsi_600[k];

        t_775[k] = f_23 * msi_434[k]
                   + f_3 * pc_y[k] * nsi_602[k];

        t_776[k] = f_10 * nsh0_450[k]
                   - f_11 * nsh1_450[k]
                   + f_3 * pc_z[k] * nsi_602[k];

        t_777[k] = f_16 * msi_609[k]
                   + f_3 * pc_x[k] * nsi_609[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, pc_x, pc_z, msi_611, msi_612, \
                         msi_613, msi_614, nsi_603, nsi_611, nsi_612, nsi_613, \
                         nsi_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_3 * pc_z[k] * nsi_603[k];

        t_779[k] = f_16 * msi_611[k]
                   + f_3 * pc_x[k] * nsi_611[k];

        t_780[k] = f_16 * msi_612[k]
                   + f_3 * pc_x[k] * nsi_612[k];

        t_781[k] = f_16 * msi_613[k]
                   + f_3 * pc_x[k] * nsi_613[k];

        t_782[k] = f_16 * msi_614[k]
                   + f_3 * pc_x[k] * nsi_614[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pc_x, pc_y, pc_z, msi_441, msi_615, \
                         nsh0_456, nsh1_456, nsi_609, nsi_610, \
                         nsi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_16 * msi_615[k]
                   + f_3 * pc_x[k] * nsi_615[k];

        t_784[k] = f_23 * msi_441[k]
                   + f_1 * nsh0_456[k]
                   - f_2 * nsh1_456[k]
                   + f_3 * pc_y[k] * nsi_609[k];

        t_785[k] = f_3 * pc_z[k] * nsi_609[k];

        t_786[k] = f_4 * nsh0_456[k]
                   - f_5 * nsh1_456[k]
                   + f_3 * pc_z[k] * nsi_610[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_z, nsh0_457, nsh0_458, nsh0_459, nsh1_457, \
                         nsh1_458, nsh1_459, nsi_611, nsi_612, \
                         nsi_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_6 * nsh0_457[k]
                   - f_7 * nsh1_457[k]
                   + f_3 * pc_z[k] * nsi_611[k];

        t_788[k] = f_8 * nsh0_458[k]
                   - f_9 * nsh1_458[k]
                   + f_3 * pc_z[k] * nsi_612[k];

        t_789[k] = f_10 * nsh0_459[k]
                   - f_11 * nsh1_459[k]
                   + f_3 * pc_z[k] * nsi_613[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_y, pc_z, msk0_540, msi_447, \
                         msi_448, msk1_540, nsh0_461, nsh1_461, nsi_615, \
                         nsi_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_23 * msi_447[k]
                   + f_3 * pc_y[k] * nsi_615[k];

        t_791[k] = f_1 * nsh0_461[k]
                   - f_2 * nsh1_461[k]
                   + f_3 * pc_z[k] * nsi_615[k];

        t_792[k] = pa_z[k] * msk0_540[k]
                   - f_12 * pc_z[k] * msk1_540[k];

        t_793[k] = f_17 * msi_448[k]
                   + f_3 * pc_y[k] * nsi_616[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pa_z, pc_y, pc_z, msk0_543, msi_420, msi_450, \
                         msk1_543, nsi_616, nsi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * msi_420[k]
                   + f_3 * pc_z[k] * nsi_616[k];

        t_795[k] = pa_z[k] * msk0_543[k]
                   - f_12 * pc_z[k] * msk1_543[k];

        t_796[k] = f_17 * msi_450[k]
                   + f_3 * pc_y[k] * nsi_618[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pa_z, pc_x, pc_z, msk0_546, msi_423, msi_621, \
                         msk1_546, nsh0_467, nsh1_467, nsi_619, \
                         nsi_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_16 * msi_621[k]
                   + f_10 * nsh0_467[k]
                   - f_11 * nsh1_467[k]
                   + f_3 * pc_x[k] * nsi_621[k];

        t_798[k] = pa_z[k] * msk0_546[k]
                   - f_12 * pc_z[k] * msk1_546[k];

        t_799[k] = f_13 * msi_423[k]
                   + f_3 * pc_z[k] * nsi_619[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pa_z, pc_x, pc_y, pc_z, msk0_550, msi_453, \
                         msi_625, msk1_550, nsh0_471, nsh1_471, nsi_621, \
                         nsi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * msi_453[k]
                   + f_3 * pc_y[k] * nsi_621[k];

        t_801[k] = f_16 * msi_625[k]
                   + f_8 * nsh0_471[k]
                   - f_9 * nsh1_471[k]
                   + f_3 * pc_x[k] * nsi_625[k];

        t_802[k] = pa_z[k] * msk0_550[k]
                   - f_12 * pc_z[k] * msk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pa_z, pc_y, pc_z, msk0_552, msi_426, msi_427, \
                         msi_457, msk1_552, nsi_622, nsi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_13 * msi_426[k]
                   + f_3 * pc_z[k] * nsi_622[k];

        t_804[k] = pa_z[k] * msk0_552[k]
                   + f_14 * msi_427[k]
                   - f_12 * pc_z[k] * msk1_552[k];

        t_805[k] = f_17 * msi_457[k]
                   + f_3 * pc_y[k] * nsi_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pa_z, pc_x, pc_z, msk0_555, msi_430, msi_630, \
                         msk1_555, nsh0_476, nsh1_476, nsi_626, \
                         nsi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_16 * msi_630[k]
                   + f_6 * nsh0_476[k]
                   - f_7 * nsh1_476[k]
                   + f_3 * pc_x[k] * nsi_630[k];

        t_807[k] = pa_z[k] * msk0_555[k]
                   - f_12 * pc_z[k] * msk1_555[k];

        t_808[k] = f_13 * msi_430[k]
                   + f_3 * pc_z[k] * nsi_626[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pa_z, pc_y, pc_z, msk0_557, msk0_558, msi_431, \
                         msi_432, msi_462, msk1_557, msk1_558, \
                         nsi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pa_z[k] * msk0_557[k]
                   + f_14 * msi_431[k]
                   - f_12 * pc_z[k] * msk1_557[k];

        t_810[k] = pa_z[k] * msk0_558[k]
                   + f_15 * msi_432[k]
                   - f_12 * pc_z[k] * msk1_558[k];

        t_811[k] = f_17 * msi_462[k]
                   + f_3 * pc_y[k] * nsi_630[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, msi_636, msi_637, msi_638, msi_639, \
                         nsh0_482, nsh1_482, nsi_636, nsi_637, nsi_638, \
                         nsi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_16 * msi_636[k]
                   + f_4 * nsh0_482[k]
                   - f_5 * nsh1_482[k]
                   + f_3 * pc_x[k] * nsi_636[k];

        t_813[k] = f_16 * msi_637[k]
                   + f_3 * pc_x[k] * nsi_637[k];

        t_814[k] = f_16 * msi_638[k]
                   + f_3 * pc_x[k] * nsi_638[k];

        t_815[k] = f_16 * msi_639[k]
                   + f_3 * pc_x[k] * nsi_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, pc_x, msi_640, msi_641, msi_642, msi_643, \
                         nsi_640, nsi_641, nsi_642, nsi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_16 * msi_640[k]
                   + f_3 * pc_x[k] * nsi_640[k];

        t_817[k] = f_16 * msi_641[k]
                   + f_3 * pc_x[k] * nsi_641[k];

        t_818[k] = f_16 * msi_642[k]
                   + f_3 * pc_x[k] * nsi_642[k];

        t_819[k] = f_16 * msi_643[k]
                   + f_3 * pc_x[k] * nsi_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_y, pc_z, msk0_568, msi_441, msi_471, \
                         msk1_568, nsh0_479, nsh1_479, nsi_637, \
                         nsi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pa_z[k] * msk0_568[k]
                   - f_12 * pc_z[k] * msk1_568[k];

        t_821[k] = f_13 * msi_441[k]
                   + f_3 * pc_z[k] * nsi_637[k];

        t_822[k] = f_17 * msi_471[k]
                   + f_10 * nsh0_479[k]
                   - f_11 * nsh1_479[k]
                   + f_3 * pc_y[k] * nsi_639[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_y, msi_472, msi_473, msi_474, nsh0_480, \
                         nsh0_481, nsh0_482, nsh1_480, nsh1_481, nsh1_482, nsi_640, nsi_641, \
                         nsi_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_17 * msi_472[k]
                   + f_8 * nsh0_480[k]
                   - f_9 * nsh1_480[k]
                   + f_3 * pc_y[k] * nsi_640[k];

        t_824[k] = f_17 * msi_473[k]
                   + f_6 * nsh0_481[k]
                   - f_7 * nsh1_481[k]
                   + f_3 * pc_y[k] * nsi_641[k];

        t_825[k] = f_17 * msi_474[k]
                   + f_4 * nsh0_482[k]
                   - f_5 * nsh1_482[k]
                   + f_3 * pc_y[k] * nsi_642[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t msi, const size_t nsh0,
                                                          const size_t nsh1, const size_t nsi,
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

    const auto *msi_447 = buffer.data(msi + 447);
    const auto *msi_448 = buffer.data(msi + 448);
    const auto *msi_451 = buffer.data(msi + 451);
    const auto *msi_454 = buffer.data(msi + 454);
    const auto *msi_458 = buffer.data(msi + 458);
    const auto *msi_469 = buffer.data(msi + 469);
    const auto *msi_475 = buffer.data(msi + 475);
    const auto *msi_476 = buffer.data(msi + 476);
    const auto *msi_478 = buffer.data(msi + 478);
    const auto *msi_479 = buffer.data(msi + 479);
    const auto *msi_481 = buffer.data(msi + 481);
    const auto *msi_482 = buffer.data(msi + 482);
    const auto *msi_485 = buffer.data(msi + 485);
    const auto *msi_486 = buffer.data(msi + 486);
    const auto *msi_490 = buffer.data(msi + 490);
    const auto *msi_497 = buffer.data(msi + 497);
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
    const auto *msi_518 = buffer.data(msi + 518);
    const auto *msi_525 = buffer.data(msi + 525);
    const auto *msi_527 = buffer.data(msi + 527);
    const auto *msi_528 = buffer.data(msi + 528);
    const auto *msi_529 = buffer.data(msi + 529);
    const auto *msi_530 = buffer.data(msi + 530);
    const auto *msi_531 = buffer.data(msi + 531);
    const auto *msi_532 = buffer.data(msi + 532);
    const auto *msi_534 = buffer.data(msi + 534);
    const auto *msi_537 = buffer.data(msi + 537);
    const auto *msi_541 = buffer.data(msi + 541);
    const auto *msi_546 = buffer.data(msi + 546);
    const auto *msi_553 = buffer.data(msi + 553);
    const auto *msi_555 = buffer.data(msi + 555);
    const auto *msi_644 = buffer.data(msi + 644);
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

    const auto *nsh0_482 = buffer.data(nsh0 + 482);
    const auto *nsh0_483 = buffer.data(nsh0 + 483);
    const auto *nsh0_486 = buffer.data(nsh0 + 486);
    const auto *nsh0_488 = buffer.data(nsh0 + 488);
    const auto *nsh0_489 = buffer.data(nsh0 + 489);
    const auto *nsh0_492 = buffer.data(nsh0 + 492);
    const auto *nsh0_493 = buffer.data(nsh0 + 493);
    const auto *nsh0_495 = buffer.data(nsh0 + 495);
    const auto *nsh0_497 = buffer.data(nsh0 + 497);
    const auto *nsh0_498 = buffer.data(nsh0 + 498);
    const auto *nsh0_500 = buffer.data(nsh0 + 500);
    const auto *nsh0_501 = buffer.data(nsh0 + 501);
    const auto *nsh0_502 = buffer.data(nsh0 + 502);
    const auto *nsh0_503 = buffer.data(nsh0 + 503);
    const auto *nsh0_504 = buffer.data(nsh0 + 504);
    const auto *nsh0_507 = buffer.data(nsh0 + 507);
    const auto *nsh0_509 = buffer.data(nsh0 + 509);
    const auto *nsh0_510 = buffer.data(nsh0 + 510);
    const auto *nsh0_513 = buffer.data(nsh0 + 513);
    const auto *nsh0_514 = buffer.data(nsh0 + 514);
    const auto *nsh0_516 = buffer.data(nsh0 + 516);
    const auto *nsh0_518 = buffer.data(nsh0 + 518);
    const auto *nsh0_519 = buffer.data(nsh0 + 519);
    const auto *nsh0_521 = buffer.data(nsh0 + 521);
    const auto *nsh0_522 = buffer.data(nsh0 + 522);
    const auto *nsh0_523 = buffer.data(nsh0 + 523);
    const auto *nsh0_524 = buffer.data(nsh0 + 524);
    const auto *nsh0_525 = buffer.data(nsh0 + 525);
    const auto *nsh0_528 = buffer.data(nsh0 + 528);
    const auto *nsh0_530 = buffer.data(nsh0 + 530);
    const auto *nsh0_531 = buffer.data(nsh0 + 531);
    const auto *nsh0_534 = buffer.data(nsh0 + 534);
    const auto *nsh0_535 = buffer.data(nsh0 + 535);
    const auto *nsh0_537 = buffer.data(nsh0 + 537);
    const auto *nsh0_539 = buffer.data(nsh0 + 539);
    const auto *nsh0_540 = buffer.data(nsh0 + 540);
    const auto *nsh0_542 = buffer.data(nsh0 + 542);
    const auto *nsh0_543 = buffer.data(nsh0 + 543);
    const auto *nsh0_545 = buffer.data(nsh0 + 545);

    const auto *nsh1_482 = buffer.data(nsh1 + 482);
    const auto *nsh1_483 = buffer.data(nsh1 + 483);
    const auto *nsh1_486 = buffer.data(nsh1 + 486);
    const auto *nsh1_488 = buffer.data(nsh1 + 488);
    const auto *nsh1_489 = buffer.data(nsh1 + 489);
    const auto *nsh1_492 = buffer.data(nsh1 + 492);
    const auto *nsh1_493 = buffer.data(nsh1 + 493);
    const auto *nsh1_495 = buffer.data(nsh1 + 495);
    const auto *nsh1_497 = buffer.data(nsh1 + 497);
    const auto *nsh1_498 = buffer.data(nsh1 + 498);
    const auto *nsh1_500 = buffer.data(nsh1 + 500);
    const auto *nsh1_501 = buffer.data(nsh1 + 501);
    const auto *nsh1_502 = buffer.data(nsh1 + 502);
    const auto *nsh1_503 = buffer.data(nsh1 + 503);
    const auto *nsh1_504 = buffer.data(nsh1 + 504);
    const auto *nsh1_507 = buffer.data(nsh1 + 507);
    const auto *nsh1_509 = buffer.data(nsh1 + 509);
    const auto *nsh1_510 = buffer.data(nsh1 + 510);
    const auto *nsh1_513 = buffer.data(nsh1 + 513);
    const auto *nsh1_514 = buffer.data(nsh1 + 514);
    const auto *nsh1_516 = buffer.data(nsh1 + 516);
    const auto *nsh1_518 = buffer.data(nsh1 + 518);
    const auto *nsh1_519 = buffer.data(nsh1 + 519);
    const auto *nsh1_521 = buffer.data(nsh1 + 521);
    const auto *nsh1_522 = buffer.data(nsh1 + 522);
    const auto *nsh1_523 = buffer.data(nsh1 + 523);
    const auto *nsh1_524 = buffer.data(nsh1 + 524);
    const auto *nsh1_525 = buffer.data(nsh1 + 525);
    const auto *nsh1_528 = buffer.data(nsh1 + 528);
    const auto *nsh1_530 = buffer.data(nsh1 + 530);
    const auto *nsh1_531 = buffer.data(nsh1 + 531);
    const auto *nsh1_534 = buffer.data(nsh1 + 534);
    const auto *nsh1_535 = buffer.data(nsh1 + 535);
    const auto *nsh1_537 = buffer.data(nsh1 + 537);
    const auto *nsh1_539 = buffer.data(nsh1 + 539);
    const auto *nsh1_540 = buffer.data(nsh1 + 540);
    const auto *nsh1_542 = buffer.data(nsh1 + 542);
    const auto *nsh1_543 = buffer.data(nsh1 + 543);
    const auto *nsh1_545 = buffer.data(nsh1 + 545);

    const auto *nsi_643 = buffer.data(nsi + 643);
    const auto *nsi_644 = buffer.data(nsi + 644);
    const auto *nsi_646 = buffer.data(nsi + 646);
    const auto *nsi_647 = buffer.data(nsi + 647);
    const auto *nsi_649 = buffer.data(nsi + 649);
    const auto *nsi_650 = buffer.data(nsi + 650);
    const auto *nsi_653 = buffer.data(nsi + 653);
    const auto *nsi_654 = buffer.data(nsi + 654);
    const auto *nsi_656 = buffer.data(nsi + 656);
    const auto *nsi_658 = buffer.data(nsi + 658);
    const auto *nsi_659 = buffer.data(nsi + 659);
    const auto *nsi_661 = buffer.data(nsi + 661);
    const auto *nsi_662 = buffer.data(nsi + 662);
    const auto *nsi_664 = buffer.data(nsi + 664);
    const auto *nsi_665 = buffer.data(nsi + 665);
    const auto *nsi_666 = buffer.data(nsi + 666);
    const auto *nsi_667 = buffer.data(nsi + 667);
    const auto *nsi_668 = buffer.data(nsi + 668);
    const auto *nsi_669 = buffer.data(nsi + 669);
    const auto *nsi_670 = buffer.data(nsi + 670);
    const auto *nsi_671 = buffer.data(nsi + 671);
    const auto *nsi_672 = buffer.data(nsi + 672);
    const auto *nsi_674 = buffer.data(nsi + 674);
    const auto *nsi_675 = buffer.data(nsi + 675);
    const auto *nsi_677 = buffer.data(nsi + 677);
    const auto *nsi_678 = buffer.data(nsi + 678);
    const auto *nsi_681 = buffer.data(nsi + 681);
    const auto *nsi_682 = buffer.data(nsi + 682);
    const auto *nsi_684 = buffer.data(nsi + 684);
    const auto *nsi_686 = buffer.data(nsi + 686);
    const auto *nsi_687 = buffer.data(nsi + 687);
    const auto *nsi_689 = buffer.data(nsi + 689);
    const auto *nsi_690 = buffer.data(nsi + 690);
    const auto *nsi_692 = buffer.data(nsi + 692);
    const auto *nsi_693 = buffer.data(nsi + 693);
    const auto *nsi_694 = buffer.data(nsi + 694);
    const auto *nsi_695 = buffer.data(nsi + 695);
    const auto *nsi_696 = buffer.data(nsi + 696);
    const auto *nsi_697 = buffer.data(nsi + 697);
    const auto *nsi_698 = buffer.data(nsi + 698);
    const auto *nsi_699 = buffer.data(nsi + 699);
    const auto *nsi_700 = buffer.data(nsi + 700);
    const auto *nsi_702 = buffer.data(nsi + 702);
    const auto *nsi_703 = buffer.data(nsi + 703);
    const auto *nsi_705 = buffer.data(nsi + 705);
    const auto *nsi_706 = buffer.data(nsi + 706);
    const auto *nsi_709 = buffer.data(nsi + 709);
    const auto *nsi_710 = buffer.data(nsi + 710);
    const auto *nsi_712 = buffer.data(nsi + 712);
    const auto *nsi_714 = buffer.data(nsi + 714);
    const auto *nsi_715 = buffer.data(nsi + 715);
    const auto *nsi_717 = buffer.data(nsi + 717);
    const auto *nsi_718 = buffer.data(nsi + 718);
    const auto *nsi_720 = buffer.data(nsi + 720);
    const auto *nsi_721 = buffer.data(nsi + 721);
    const auto *nsi_722 = buffer.data(nsi + 722);
    const auto *nsi_723 = buffer.data(nsi + 723);
    const auto *nsi_724 = buffer.data(nsi + 724);
    const auto *nsi_725 = buffer.data(nsi + 725);
    const auto *nsi_726 = buffer.data(nsi + 726);
    const auto *nsi_727 = buffer.data(nsi + 727);

#pragma omp simd aligned(t_826, t_827, t_828, pc_x, pc_y, pc_z, msi_447, msi_475, msi_644, \
                         nsh0_482, nsh0_483, nsh1_482, nsh1_483, nsi_643, \
                         nsi_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * msi_475[k]
                   + f_3 * pc_y[k] * nsi_643[k];

        t_827[k] = f_13 * msi_447[k]
                   + f_1 * nsh0_482[k]
                   - f_2 * nsh1_482[k]
                   + f_3 * pc_z[k] * nsi_643[k];

        t_828[k] = f_16 * msi_644[k]
                   + f_1 * nsh0_483[k]
                   - f_2 * nsh1_483[k]
                   + f_3 * pc_x[k] * nsi_644[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, pc_y, pc_z, msi_448, msi_476, \
                         msi_478, msi_647, nsh0_486, nsh1_486, nsi_644, nsi_646, \
                         nsi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_16 * msi_476[k]
                   + f_3 * pc_y[k] * nsi_644[k];

        t_830[k] = f_14 * msi_448[k]
                   + f_3 * pc_z[k] * nsi_644[k];

        t_831[k] = f_16 * msi_647[k]
                   + f_10 * nsh0_486[k]
                   - f_11 * nsh1_486[k]
                   + f_3 * pc_x[k] * nsi_647[k];

        t_832[k] = f_16 * msi_478[k]
                   + f_3 * pc_y[k] * nsi_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, pc_z, msi_451, msi_649, msi_650, nsh0_488, \
                         nsh0_489, nsh1_488, nsh1_489, nsi_647, nsi_649, \
                         nsi_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_16 * msi_649[k]
                   + f_10 * nsh0_488[k]
                   - f_11 * nsh1_488[k]
                   + f_3 * pc_x[k] * nsi_649[k];

        t_834[k] = f_16 * msi_650[k]
                   + f_8 * nsh0_489[k]
                   - f_9 * nsh1_489[k]
                   + f_3 * pc_x[k] * nsi_650[k];

        t_835[k] = f_14 * msi_451[k]
                   + f_3 * pc_z[k] * nsi_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_x, pc_y, msi_481, msi_653, msi_654, nsh0_492, \
                         nsh0_493, nsh1_492, nsh1_493, nsi_649, nsi_653, \
                         nsi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * msi_481[k]
                   + f_3 * pc_y[k] * nsi_649[k];

        t_837[k] = f_16 * msi_653[k]
                   + f_8 * nsh0_492[k]
                   - f_9 * nsh1_492[k]
                   + f_3 * pc_x[k] * nsi_653[k];

        t_838[k] = f_16 * msi_654[k]
                   + f_6 * nsh0_493[k]
                   - f_7 * nsh1_493[k]
                   + f_3 * pc_x[k] * nsi_654[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_y, pc_z, msi_454, msi_485, msi_656, \
                         nsh0_495, nsh1_495, nsi_650, nsi_653, \
                         nsi_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * msi_454[k]
                   + f_3 * pc_z[k] * nsi_650[k];

        t_840[k] = f_16 * msi_656[k]
                   + f_6 * nsh0_495[k]
                   - f_7 * nsh1_495[k]
                   + f_3 * pc_x[k] * nsi_656[k];

        t_841[k] = f_16 * msi_485[k]
                   + f_3 * pc_y[k] * nsi_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_z, msi_458, msi_658, msi_659, nsh0_497, \
                         nsh0_498, nsh1_497, nsh1_498, nsi_654, nsi_658, \
                         nsi_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_16 * msi_658[k]
                   + f_6 * nsh0_497[k]
                   - f_7 * nsh1_497[k]
                   + f_3 * pc_x[k] * nsi_658[k];

        t_843[k] = f_16 * msi_659[k]
                   + f_4 * nsh0_498[k]
                   - f_5 * nsh1_498[k]
                   + f_3 * pc_x[k] * nsi_659[k];

        t_844[k] = f_14 * msi_458[k]
                   + f_3 * pc_z[k] * nsi_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, msi_490, msi_661, msi_662, nsh0_500, \
                         nsh0_501, nsh1_500, nsh1_501, nsi_658, nsi_661, \
                         nsi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_16 * msi_661[k]
                   + f_4 * nsh0_500[k]
                   - f_5 * nsh1_500[k]
                   + f_3 * pc_x[k] * nsi_661[k];

        t_846[k] = f_16 * msi_662[k]
                   + f_4 * nsh0_501[k]
                   - f_5 * nsh1_501[k]
                   + f_3 * pc_x[k] * nsi_662[k];

        t_847[k] = f_16 * msi_490[k]
                   + f_3 * pc_y[k] * nsi_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, msi_664, msi_665, msi_666, msi_667, \
                         nsh0_503, nsh1_503, nsi_664, nsi_665, nsi_666, \
                         nsi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * msi_664[k]
                   + f_4 * nsh0_503[k]
                   - f_5 * nsh1_503[k]
                   + f_3 * pc_x[k] * nsi_664[k];

        t_849[k] = f_16 * msi_665[k]
                   + f_3 * pc_x[k] * nsi_665[k];

        t_850[k] = f_16 * msi_666[k]
                   + f_3 * pc_x[k] * nsi_666[k];

        t_851[k] = f_16 * msi_667[k]
                   + f_3 * pc_x[k] * nsi_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, msi_668, msi_669, msi_670, msi_671, \
                         nsi_668, nsi_669, nsi_670, nsi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_16 * msi_668[k]
                   + f_3 * pc_x[k] * nsi_668[k];

        t_853[k] = f_16 * msi_669[k]
                   + f_3 * pc_x[k] * nsi_669[k];

        t_854[k] = f_16 * msi_670[k]
                   + f_3 * pc_x[k] * nsi_670[k];

        t_855[k] = f_16 * msi_671[k]
                   + f_3 * pc_x[k] * nsi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, msi_469, msi_497, msi_499, nsh0_498, \
                         nsh0_500, nsh1_498, nsh1_500, nsi_665, \
                         nsi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * msi_497[k]
                   + f_1 * nsh0_498[k]
                   - f_2 * nsh1_498[k]
                   + f_3 * pc_y[k] * nsi_665[k];

        t_857[k] = f_14 * msi_469[k]
                   + f_3 * pc_z[k] * nsi_665[k];

        t_858[k] = f_16 * msi_499[k]
                   + f_10 * nsh0_500[k]
                   - f_11 * nsh1_500[k]
                   + f_3 * pc_y[k] * nsi_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, msi_500, msi_501, msi_502, nsh0_501, \
                         nsh0_502, nsh0_503, nsh1_501, nsh1_502, nsh1_503, nsi_668, nsi_669, \
                         nsi_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * msi_500[k]
                   + f_8 * nsh0_501[k]
                   - f_9 * nsh1_501[k]
                   + f_3 * pc_y[k] * nsi_668[k];

        t_860[k] = f_16 * msi_501[k]
                   + f_6 * nsh0_502[k]
                   - f_7 * nsh1_502[k]
                   + f_3 * pc_y[k] * nsi_669[k];

        t_861[k] = f_16 * msi_502[k]
                   + f_4 * nsh0_503[k]
                   - f_5 * nsh1_503[k]
                   + f_3 * pc_y[k] * nsi_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, msi_475, msi_503, msi_672, \
                         nsh0_503, nsh0_504, nsh1_503, nsh1_504, nsi_671, \
                         nsi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * msi_503[k]
                   + f_3 * pc_y[k] * nsi_671[k];

        t_863[k] = f_14 * msi_475[k]
                   + f_1 * nsh0_503[k]
                   - f_2 * nsh1_503[k]
                   + f_3 * pc_z[k] * nsi_671[k];

        t_864[k] = f_16 * msi_672[k]
                   + f_1 * nsh0_504[k]
                   - f_2 * nsh1_504[k]
                   + f_3 * pc_x[k] * nsi_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, msi_476, msi_504, \
                         msi_506, msi_675, nsh0_507, nsh1_507, nsi_672, nsi_674, \
                         nsi_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_15 * msi_504[k]
                   + f_3 * pc_y[k] * nsi_672[k];

        t_866[k] = f_15 * msi_476[k]
                   + f_3 * pc_z[k] * nsi_672[k];

        t_867[k] = f_16 * msi_675[k]
                   + f_10 * nsh0_507[k]
                   - f_11 * nsh1_507[k]
                   + f_3 * pc_x[k] * nsi_675[k];

        t_868[k] = f_15 * msi_506[k]
                   + f_3 * pc_y[k] * nsi_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, msi_479, msi_677, msi_678, nsh0_509, \
                         nsh0_510, nsh1_509, nsh1_510, nsi_675, nsi_677, \
                         nsi_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_16 * msi_677[k]
                   + f_10 * nsh0_509[k]
                   - f_11 * nsh1_509[k]
                   + f_3 * pc_x[k] * nsi_677[k];

        t_870[k] = f_16 * msi_678[k]
                   + f_8 * nsh0_510[k]
                   - f_9 * nsh1_510[k]
                   + f_3 * pc_x[k] * nsi_678[k];

        t_871[k] = f_15 * msi_479[k]
                   + f_3 * pc_z[k] * nsi_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, pc_y, msi_509, msi_681, msi_682, nsh0_513, \
                         nsh0_514, nsh1_513, nsh1_514, nsi_677, nsi_681, \
                         nsi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * msi_509[k]
                   + f_3 * pc_y[k] * nsi_677[k];

        t_873[k] = f_16 * msi_681[k]
                   + f_8 * nsh0_513[k]
                   - f_9 * nsh1_513[k]
                   + f_3 * pc_x[k] * nsi_681[k];

        t_874[k] = f_16 * msi_682[k]
                   + f_6 * nsh0_514[k]
                   - f_7 * nsh1_514[k]
                   + f_3 * pc_x[k] * nsi_682[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, pc_y, pc_z, msi_482, msi_513, msi_684, \
                         nsh0_516, nsh1_516, nsi_678, nsi_681, \
                         nsi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * msi_482[k]
                   + f_3 * pc_z[k] * nsi_678[k];

        t_876[k] = f_16 * msi_684[k]
                   + f_6 * nsh0_516[k]
                   - f_7 * nsh1_516[k]
                   + f_3 * pc_x[k] * nsi_684[k];

        t_877[k] = f_15 * msi_513[k]
                   + f_3 * pc_y[k] * nsi_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, pc_z, msi_486, msi_686, msi_687, nsh0_518, \
                         nsh0_519, nsh1_518, nsh1_519, nsi_682, nsi_686, \
                         nsi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_16 * msi_686[k]
                   + f_6 * nsh0_518[k]
                   - f_7 * nsh1_518[k]
                   + f_3 * pc_x[k] * nsi_686[k];

        t_879[k] = f_16 * msi_687[k]
                   + f_4 * nsh0_519[k]
                   - f_5 * nsh1_519[k]
                   + f_3 * pc_x[k] * nsi_687[k];

        t_880[k] = f_15 * msi_486[k]
                   + f_3 * pc_z[k] * nsi_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_y, msi_518, msi_689, msi_690, nsh0_521, \
                         nsh0_522, nsh1_521, nsh1_522, nsi_686, nsi_689, \
                         nsi_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_16 * msi_689[k]
                   + f_4 * nsh0_521[k]
                   - f_5 * nsh1_521[k]
                   + f_3 * pc_x[k] * nsi_689[k];

        t_882[k] = f_16 * msi_690[k]
                   + f_4 * nsh0_522[k]
                   - f_5 * nsh1_522[k]
                   + f_3 * pc_x[k] * nsi_690[k];

        t_883[k] = f_15 * msi_518[k]
                   + f_3 * pc_y[k] * nsi_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pc_x, msi_692, msi_693, msi_694, msi_695, \
                         nsh0_524, nsh1_524, nsi_692, nsi_693, nsi_694, \
                         nsi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_16 * msi_692[k]
                   + f_4 * nsh0_524[k]
                   - f_5 * nsh1_524[k]
                   + f_3 * pc_x[k] * nsi_692[k];

        t_885[k] = f_16 * msi_693[k]
                   + f_3 * pc_x[k] * nsi_693[k];

        t_886[k] = f_16 * msi_694[k]
                   + f_3 * pc_x[k] * nsi_694[k];

        t_887[k] = f_16 * msi_695[k]
                   + f_3 * pc_x[k] * nsi_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, msi_696, msi_697, msi_698, msi_699, \
                         nsi_696, nsi_697, nsi_698, nsi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_16 * msi_696[k]
                   + f_3 * pc_x[k] * nsi_696[k];

        t_889[k] = f_16 * msi_697[k]
                   + f_3 * pc_x[k] * nsi_697[k];

        t_890[k] = f_16 * msi_698[k]
                   + f_3 * pc_x[k] * nsi_698[k];

        t_891[k] = f_16 * msi_699[k]
                   + f_3 * pc_x[k] * nsi_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, msi_497, msi_525, msi_527, nsh0_519, \
                         nsh0_521, nsh1_519, nsh1_521, nsi_693, \
                         nsi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * msi_525[k]
                   + f_1 * nsh0_519[k]
                   - f_2 * nsh1_519[k]
                   + f_3 * pc_y[k] * nsi_693[k];

        t_893[k] = f_15 * msi_497[k]
                   + f_3 * pc_z[k] * nsi_693[k];

        t_894[k] = f_15 * msi_527[k]
                   + f_10 * nsh0_521[k]
                   - f_11 * nsh1_521[k]
                   + f_3 * pc_y[k] * nsi_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, msi_528, msi_529, msi_530, nsh0_522, \
                         nsh0_523, nsh0_524, nsh1_522, nsh1_523, nsh1_524, nsi_696, nsi_697, \
                         nsi_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * msi_528[k]
                   + f_8 * nsh0_522[k]
                   - f_9 * nsh1_522[k]
                   + f_3 * pc_y[k] * nsi_696[k];

        t_896[k] = f_15 * msi_529[k]
                   + f_6 * nsh0_523[k]
                   - f_7 * nsh1_523[k]
                   + f_3 * pc_y[k] * nsi_697[k];

        t_897[k] = f_15 * msi_530[k]
                   + f_4 * nsh0_524[k]
                   - f_5 * nsh1_524[k]
                   + f_3 * pc_y[k] * nsi_698[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_x, pc_y, pc_z, msi_503, msi_531, msi_700, \
                         nsh0_524, nsh0_525, nsh1_524, nsh1_525, nsi_699, \
                         nsi_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * msi_531[k]
                   + f_3 * pc_y[k] * nsi_699[k];

        t_899[k] = f_15 * msi_503[k]
                   + f_1 * nsh0_524[k]
                   - f_2 * nsh1_524[k]
                   + f_3 * pc_z[k] * nsi_699[k];

        t_900[k] = f_16 * msi_700[k]
                   + f_1 * nsh0_525[k]
                   - f_2 * nsh1_525[k]
                   + f_3 * pc_x[k] * nsi_700[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pc_x, pc_y, pc_z, msi_504, msi_532, \
                         msi_534, msi_703, nsh0_528, nsh1_528, nsi_700, nsi_702, \
                         nsi_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_14 * msi_532[k]
                   + f_3 * pc_y[k] * nsi_700[k];

        t_902[k] = f_16 * msi_504[k]
                   + f_3 * pc_z[k] * nsi_700[k];

        t_903[k] = f_16 * msi_703[k]
                   + f_10 * nsh0_528[k]
                   - f_11 * nsh1_528[k]
                   + f_3 * pc_x[k] * nsi_703[k];

        t_904[k] = f_14 * msi_534[k]
                   + f_3 * pc_y[k] * nsi_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, msi_507, msi_705, msi_706, nsh0_530, \
                         nsh0_531, nsh1_530, nsh1_531, nsi_703, nsi_705, \
                         nsi_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_16 * msi_705[k]
                   + f_10 * nsh0_530[k]
                   - f_11 * nsh1_530[k]
                   + f_3 * pc_x[k] * nsi_705[k];

        t_906[k] = f_16 * msi_706[k]
                   + f_8 * nsh0_531[k]
                   - f_9 * nsh1_531[k]
                   + f_3 * pc_x[k] * nsi_706[k];

        t_907[k] = f_16 * msi_507[k]
                   + f_3 * pc_z[k] * nsi_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, msi_537, msi_709, msi_710, nsh0_534, \
                         nsh0_535, nsh1_534, nsh1_535, nsi_705, nsi_709, \
                         nsi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * msi_537[k]
                   + f_3 * pc_y[k] * nsi_705[k];

        t_909[k] = f_16 * msi_709[k]
                   + f_8 * nsh0_534[k]
                   - f_9 * nsh1_534[k]
                   + f_3 * pc_x[k] * nsi_709[k];

        t_910[k] = f_16 * msi_710[k]
                   + f_6 * nsh0_535[k]
                   - f_7 * nsh1_535[k]
                   + f_3 * pc_x[k] * nsi_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, pc_y, pc_z, msi_510, msi_541, msi_712, \
                         nsh0_537, nsh1_537, nsi_706, nsi_709, \
                         nsi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * msi_510[k]
                   + f_3 * pc_z[k] * nsi_706[k];

        t_912[k] = f_16 * msi_712[k]
                   + f_6 * nsh0_537[k]
                   - f_7 * nsh1_537[k]
                   + f_3 * pc_x[k] * nsi_712[k];

        t_913[k] = f_14 * msi_541[k]
                   + f_3 * pc_y[k] * nsi_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, msi_514, msi_714, msi_715, nsh0_539, \
                         nsh0_540, nsh1_539, nsh1_540, nsi_710, nsi_714, \
                         nsi_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_16 * msi_714[k]
                   + f_6 * nsh0_539[k]
                   - f_7 * nsh1_539[k]
                   + f_3 * pc_x[k] * nsi_714[k];

        t_915[k] = f_16 * msi_715[k]
                   + f_4 * nsh0_540[k]
                   - f_5 * nsh1_540[k]
                   + f_3 * pc_x[k] * nsi_715[k];

        t_916[k] = f_16 * msi_514[k]
                   + f_3 * pc_z[k] * nsi_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, pc_y, msi_546, msi_717, msi_718, nsh0_542, \
                         nsh0_543, nsh1_542, nsh1_543, nsi_714, nsi_717, \
                         nsi_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_16 * msi_717[k]
                   + f_4 * nsh0_542[k]
                   - f_5 * nsh1_542[k]
                   + f_3 * pc_x[k] * nsi_717[k];

        t_918[k] = f_16 * msi_718[k]
                   + f_4 * nsh0_543[k]
                   - f_5 * nsh1_543[k]
                   + f_3 * pc_x[k] * nsi_718[k];

        t_919[k] = f_14 * msi_546[k]
                   + f_3 * pc_y[k] * nsi_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, msi_720, msi_721, msi_722, msi_723, \
                         nsh0_545, nsh1_545, nsi_720, nsi_721, nsi_722, \
                         nsi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_16 * msi_720[k]
                   + f_4 * nsh0_545[k]
                   - f_5 * nsh1_545[k]
                   + f_3 * pc_x[k] * nsi_720[k];

        t_921[k] = f_16 * msi_721[k]
                   + f_3 * pc_x[k] * nsi_721[k];

        t_922[k] = f_16 * msi_722[k]
                   + f_3 * pc_x[k] * nsi_722[k];

        t_923[k] = f_16 * msi_723[k]
                   + f_3 * pc_x[k] * nsi_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, msi_724, msi_725, msi_726, msi_727, \
                         nsi_724, nsi_725, nsi_726, nsi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_16 * msi_724[k]
                   + f_3 * pc_x[k] * nsi_724[k];

        t_925[k] = f_16 * msi_725[k]
                   + f_3 * pc_x[k] * nsi_725[k];

        t_926[k] = f_16 * msi_726[k]
                   + f_3 * pc_x[k] * nsi_726[k];

        t_927[k] = f_16 * msi_727[k]
                   + f_3 * pc_x[k] * nsi_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, msi_525, msi_553, msi_555, nsh0_540, \
                         nsh0_542, nsh1_540, nsh1_542, nsi_721, \
                         nsi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * msi_553[k]
                   + f_1 * nsh0_540[k]
                   - f_2 * nsh1_540[k]
                   + f_3 * pc_y[k] * nsi_721[k];

        t_929[k] = f_16 * msi_525[k]
                   + f_3 * pc_z[k] * nsi_721[k];

        t_930[k] = f_14 * msi_555[k]
                   + f_10 * nsh0_542[k]
                   - f_11 * nsh1_542[k]
                   + f_3 * pc_y[k] * nsi_723[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_720 = buffer.data(msk0 + 720);
    const auto *msk0_723 = buffer.data(msk0 + 723);
    const auto *msk0_725 = buffer.data(msk0 + 725);
    const auto *msk0_726 = buffer.data(msk0 + 726);
    const auto *msk0_729 = buffer.data(msk0 + 729);
    const auto *msk0_730 = buffer.data(msk0 + 730);
    const auto *msk0_732 = buffer.data(msk0 + 732);
    const auto *msk0_734 = buffer.data(msk0 + 734);
    const auto *msk0_735 = buffer.data(msk0 + 735);
    const auto *msk0_737 = buffer.data(msk0 + 737);
    const auto *msk0_738 = buffer.data(msk0 + 738);
    const auto *msk0_740 = buffer.data(msk0 + 740);
    const auto *msk0_755 = buffer.data(msk0 + 755);
    const auto *msk0_756 = buffer.data(msk0 + 756);
    const auto *msk0_759 = buffer.data(msk0 + 759);

    const auto *msi_531 = buffer.data(msi + 531);
    const auto *msi_532 = buffer.data(msi + 532);
    const auto *msi_535 = buffer.data(msi + 535);
    const auto *msi_538 = buffer.data(msi + 538);
    const auto *msi_542 = buffer.data(msi + 542);
    const auto *msi_553 = buffer.data(msi + 553);
    const auto *msi_556 = buffer.data(msi + 556);
    const auto *msi_557 = buffer.data(msi + 557);
    const auto *msi_558 = buffer.data(msi + 558);
    const auto *msi_559 = buffer.data(msi + 559);
    const auto *msi_560 = buffer.data(msi + 560);
    const auto *msi_561 = buffer.data(msi + 561);
    const auto *msi_562 = buffer.data(msi + 562);
    const auto *msi_563 = buffer.data(msi + 563);
    const auto *msi_565 = buffer.data(msi + 565);
    const auto *msi_566 = buffer.data(msi + 566);
    const auto *msi_568 = buffer.data(msi + 568);
    const auto *msi_569 = buffer.data(msi + 569);
    const auto *msi_570 = buffer.data(msi + 570);
    const auto *msi_572 = buffer.data(msi + 572);
    const auto *msi_573 = buffer.data(msi + 573);
    const auto *msi_574 = buffer.data(msi + 574);
    const auto *msi_581 = buffer.data(msi + 581);
    const auto *msi_583 = buffer.data(msi + 583);
    const auto *msi_584 = buffer.data(msi + 584);
    const auto *msi_585 = buffer.data(msi + 585);
    const auto *msi_586 = buffer.data(msi + 586);
    const auto *msi_587 = buffer.data(msi + 587);
    const auto *msi_588 = buffer.data(msi + 588);
    const auto *msi_593 = buffer.data(msi + 593);
    const auto *msi_597 = buffer.data(msi + 597);
    const auto *msi_602 = buffer.data(msi + 602);
    const auto *msi_609 = buffer.data(msi + 609);
    const auto *msi_615 = buffer.data(msi + 615);
    const auto *msi_616 = buffer.data(msi + 616);
    const auto *msi_618 = buffer.data(msi + 618);
    const auto *msi_749 = buffer.data(msi + 749);
    const auto *msi_750 = buffer.data(msi + 750);
    const auto *msi_751 = buffer.data(msi + 751);
    const auto *msi_752 = buffer.data(msi + 752);
    const auto *msi_753 = buffer.data(msi + 753);
    const auto *msi_754 = buffer.data(msi + 754);
    const auto *msi_755 = buffer.data(msi + 755);
    const auto *msi_756 = buffer.data(msi + 756);
    const auto *msi_761 = buffer.data(msi + 761);
    const auto *msi_765 = buffer.data(msi + 765);
    const auto *msi_770 = buffer.data(msi + 770);
    const auto *msi_776 = buffer.data(msi + 776);
    const auto *msi_777 = buffer.data(msi + 777);
    const auto *msi_778 = buffer.data(msi + 778);
    const auto *msi_779 = buffer.data(msi + 779);
    const auto *msi_780 = buffer.data(msi + 780);
    const auto *msi_781 = buffer.data(msi + 781);
    const auto *msi_783 = buffer.data(msi + 783);
    const auto *msi_784 = buffer.data(msi + 784);
    const auto *msi_787 = buffer.data(msi + 787);
    const auto *msi_790 = buffer.data(msi + 790);
    const auto *msi_794 = buffer.data(msi + 794);
    const auto *msi_799 = buffer.data(msi + 799);
    const auto *msi_805 = buffer.data(msi + 805);
    const auto *msi_807 = buffer.data(msi + 807);
    const auto *msi_808 = buffer.data(msi + 808);
    const auto *msi_809 = buffer.data(msi + 809);
    const auto *msi_810 = buffer.data(msi + 810);
    const auto *msi_811 = buffer.data(msi + 811);

    const auto *msk1_720 = buffer.data(msk1 + 720);
    const auto *msk1_723 = buffer.data(msk1 + 723);
    const auto *msk1_725 = buffer.data(msk1 + 725);
    const auto *msk1_726 = buffer.data(msk1 + 726);
    const auto *msk1_729 = buffer.data(msk1 + 729);
    const auto *msk1_730 = buffer.data(msk1 + 730);
    const auto *msk1_732 = buffer.data(msk1 + 732);
    const auto *msk1_734 = buffer.data(msk1 + 734);
    const auto *msk1_735 = buffer.data(msk1 + 735);
    const auto *msk1_737 = buffer.data(msk1 + 737);
    const auto *msk1_738 = buffer.data(msk1 + 738);
    const auto *msk1_740 = buffer.data(msk1 + 740);
    const auto *msk1_755 = buffer.data(msk1 + 755);
    const auto *msk1_756 = buffer.data(msk1 + 756);
    const auto *msk1_759 = buffer.data(msk1 + 759);

    const auto *nsh0_543 = buffer.data(nsh0 + 543);
    const auto *nsh0_544 = buffer.data(nsh0 + 544);
    const auto *nsh0_545 = buffer.data(nsh0 + 545);
    const auto *nsh0_561 = buffer.data(nsh0 + 561);
    const auto *nsh0_563 = buffer.data(nsh0 + 563);
    const auto *nsh0_564 = buffer.data(nsh0 + 564);
    const auto *nsh0_565 = buffer.data(nsh0 + 565);
    const auto *nsh0_566 = buffer.data(nsh0 + 566);
    const auto *nsh0_567 = buffer.data(nsh0 + 567);
    const auto *nsh0_568 = buffer.data(nsh0 + 568);
    const auto *nsh0_569 = buffer.data(nsh0 + 569);
    const auto *nsh0_570 = buffer.data(nsh0 + 570);
    const auto *nsh0_571 = buffer.data(nsh0 + 571);
    const auto *nsh0_572 = buffer.data(nsh0 + 572);
    const auto *nsh0_573 = buffer.data(nsh0 + 573);
    const auto *nsh0_574 = buffer.data(nsh0 + 574);
    const auto *nsh0_575 = buffer.data(nsh0 + 575);
    const auto *nsh0_576 = buffer.data(nsh0 + 576);
    const auto *nsh0_581 = buffer.data(nsh0 + 581);
    const auto *nsh0_582 = buffer.data(nsh0 + 582);
    const auto *nsh0_583 = buffer.data(nsh0 + 583);
    const auto *nsh0_584 = buffer.data(nsh0 + 584);
    const auto *nsh0_585 = buffer.data(nsh0 + 585);
    const auto *nsh0_586 = buffer.data(nsh0 + 586);
    const auto *nsh0_587 = buffer.data(nsh0 + 587);
    const auto *nsh0_588 = buffer.data(nsh0 + 588);
    const auto *nsh0_590 = buffer.data(nsh0 + 590);
    const auto *nsh0_591 = buffer.data(nsh0 + 591);
    const auto *nsh0_593 = buffer.data(nsh0 + 593);
    const auto *nsh0_594 = buffer.data(nsh0 + 594);
    const auto *nsh0_595 = buffer.data(nsh0 + 595);
    const auto *nsh0_597 = buffer.data(nsh0 + 597);
    const auto *nsh0_598 = buffer.data(nsh0 + 598);
    const auto *nsh0_603 = buffer.data(nsh0 + 603);
    const auto *nsh0_604 = buffer.data(nsh0 + 604);
    const auto *nsh0_605 = buffer.data(nsh0 + 605);
    const auto *nsh0_606 = buffer.data(nsh0 + 606);
    const auto *nsh0_608 = buffer.data(nsh0 + 608);

    const auto *nsh1_543 = buffer.data(nsh1 + 543);
    const auto *nsh1_544 = buffer.data(nsh1 + 544);
    const auto *nsh1_545 = buffer.data(nsh1 + 545);
    const auto *nsh1_561 = buffer.data(nsh1 + 561);
    const auto *nsh1_563 = buffer.data(nsh1 + 563);
    const auto *nsh1_564 = buffer.data(nsh1 + 564);
    const auto *nsh1_565 = buffer.data(nsh1 + 565);
    const auto *nsh1_566 = buffer.data(nsh1 + 566);
    const auto *nsh1_567 = buffer.data(nsh1 + 567);
    const auto *nsh1_568 = buffer.data(nsh1 + 568);
    const auto *nsh1_569 = buffer.data(nsh1 + 569);
    const auto *nsh1_570 = buffer.data(nsh1 + 570);
    const auto *nsh1_571 = buffer.data(nsh1 + 571);
    const auto *nsh1_572 = buffer.data(nsh1 + 572);
    const auto *nsh1_573 = buffer.data(nsh1 + 573);
    const auto *nsh1_574 = buffer.data(nsh1 + 574);
    const auto *nsh1_575 = buffer.data(nsh1 + 575);
    const auto *nsh1_576 = buffer.data(nsh1 + 576);
    const auto *nsh1_581 = buffer.data(nsh1 + 581);
    const auto *nsh1_582 = buffer.data(nsh1 + 582);
    const auto *nsh1_583 = buffer.data(nsh1 + 583);
    const auto *nsh1_584 = buffer.data(nsh1 + 584);
    const auto *nsh1_585 = buffer.data(nsh1 + 585);
    const auto *nsh1_586 = buffer.data(nsh1 + 586);
    const auto *nsh1_587 = buffer.data(nsh1 + 587);
    const auto *nsh1_588 = buffer.data(nsh1 + 588);
    const auto *nsh1_590 = buffer.data(nsh1 + 590);
    const auto *nsh1_591 = buffer.data(nsh1 + 591);
    const auto *nsh1_593 = buffer.data(nsh1 + 593);
    const auto *nsh1_594 = buffer.data(nsh1 + 594);
    const auto *nsh1_595 = buffer.data(nsh1 + 595);
    const auto *nsh1_597 = buffer.data(nsh1 + 597);
    const auto *nsh1_598 = buffer.data(nsh1 + 598);
    const auto *nsh1_603 = buffer.data(nsh1 + 603);
    const auto *nsh1_604 = buffer.data(nsh1 + 604);
    const auto *nsh1_605 = buffer.data(nsh1 + 605);
    const auto *nsh1_606 = buffer.data(nsh1 + 606);
    const auto *nsh1_608 = buffer.data(nsh1 + 608);

    const auto *nsi_724 = buffer.data(nsi + 724);
    const auto *nsi_725 = buffer.data(nsi + 725);
    const auto *nsi_726 = buffer.data(nsi + 726);
    const auto *nsi_727 = buffer.data(nsi + 727);
    const auto *nsi_728 = buffer.data(nsi + 728);
    const auto *nsi_730 = buffer.data(nsi + 730);
    const auto *nsi_731 = buffer.data(nsi + 731);
    const auto *nsi_733 = buffer.data(nsi + 733);
    const auto *nsi_734 = buffer.data(nsi + 734);
    const auto *nsi_737 = buffer.data(nsi + 737);
    const auto *nsi_738 = buffer.data(nsi + 738);
    const auto *nsi_742 = buffer.data(nsi + 742);
    const auto *nsi_749 = buffer.data(nsi + 749);
    const auto *nsi_750 = buffer.data(nsi + 750);
    const auto *nsi_751 = buffer.data(nsi + 751);
    const auto *nsi_752 = buffer.data(nsi + 752);
    const auto *nsi_753 = buffer.data(nsi + 753);
    const auto *nsi_754 = buffer.data(nsi + 754);
    const auto *nsi_755 = buffer.data(nsi + 755);
    const auto *nsi_756 = buffer.data(nsi + 756);
    const auto *nsi_757 = buffer.data(nsi + 757);
    const auto *nsi_758 = buffer.data(nsi + 758);
    const auto *nsi_759 = buffer.data(nsi + 759);
    const auto *nsi_760 = buffer.data(nsi + 760);
    const auto *nsi_761 = buffer.data(nsi + 761);
    const auto *nsi_762 = buffer.data(nsi + 762);
    const auto *nsi_763 = buffer.data(nsi + 763);
    const auto *nsi_764 = buffer.data(nsi + 764);
    const auto *nsi_765 = buffer.data(nsi + 765);
    const auto *nsi_766 = buffer.data(nsi + 766);
    const auto *nsi_767 = buffer.data(nsi + 767);
    const auto *nsi_768 = buffer.data(nsi + 768);
    const auto *nsi_769 = buffer.data(nsi + 769);
    const auto *nsi_770 = buffer.data(nsi + 770);
    const auto *nsi_776 = buffer.data(nsi + 776);
    const auto *nsi_777 = buffer.data(nsi + 777);
    const auto *nsi_778 = buffer.data(nsi + 778);
    const auto *nsi_779 = buffer.data(nsi + 779);
    const auto *nsi_780 = buffer.data(nsi + 780);
    const auto *nsi_781 = buffer.data(nsi + 781);
    const auto *nsi_782 = buffer.data(nsi + 782);
    const auto *nsi_783 = buffer.data(nsi + 783);
    const auto *nsi_784 = buffer.data(nsi + 784);
    const auto *nsi_785 = buffer.data(nsi + 785);
    const auto *nsi_786 = buffer.data(nsi + 786);
    const auto *nsi_787 = buffer.data(nsi + 787);
    const auto *nsi_789 = buffer.data(nsi + 789);
    const auto *nsi_790 = buffer.data(nsi + 790);
    const auto *nsi_791 = buffer.data(nsi + 791);
    const auto *nsi_793 = buffer.data(nsi + 793);
    const auto *nsi_794 = buffer.data(nsi + 794);
    const auto *nsi_795 = buffer.data(nsi + 795);
    const auto *nsi_796 = buffer.data(nsi + 796);
    const auto *nsi_798 = buffer.data(nsi + 798);
    const auto *nsi_799 = buffer.data(nsi + 799);
    const auto *nsi_805 = buffer.data(nsi + 805);
    const auto *nsi_806 = buffer.data(nsi + 806);
    const auto *nsi_807 = buffer.data(nsi + 807);
    const auto *nsi_808 = buffer.data(nsi + 808);
    const auto *nsi_809 = buffer.data(nsi + 809);
    const auto *nsi_810 = buffer.data(nsi + 810);
    const auto *nsi_811 = buffer.data(nsi + 811);
    const auto *nsi_812 = buffer.data(nsi + 812);
    const auto *nsi_814 = buffer.data(nsi + 814);

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, msi_556, msi_557, msi_558, nsh0_543, \
                         nsh0_544, nsh0_545, nsh1_543, nsh1_544, nsh1_545, nsi_724, nsi_725, \
                         nsi_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * msi_556[k]
                   + f_8 * nsh0_543[k]
                   - f_9 * nsh1_543[k]
                   + f_3 * pc_y[k] * nsi_724[k];

        t_932[k] = f_14 * msi_557[k]
                   + f_6 * nsh0_544[k]
                   - f_7 * nsh1_544[k]
                   + f_3 * pc_y[k] * nsi_725[k];

        t_933[k] = f_14 * msi_558[k]
                   + f_4 * nsh0_545[k]
                   - f_5 * nsh1_545[k]
                   + f_3 * pc_y[k] * nsi_726[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pa_y, pc_y, pc_z, msk0_720, msi_531, \
                         msi_559, msi_560, msk1_720, nsh0_545, nsh1_545, nsi_727, \
                         nsi_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * msi_559[k]
                   + f_3 * pc_y[k] * nsi_727[k];

        t_935[k] = f_16 * msi_531[k]
                   + f_1 * nsh0_545[k]
                   - f_2 * nsh1_545[k]
                   + f_3 * pc_z[k] * nsi_727[k];

        t_936[k] = pa_y[k] * msk0_720[k]
                   - f_12 * pc_y[k] * msk1_720[k];

        t_937[k] = f_13 * msi_560[k]
                   + f_3 * pc_y[k] * nsi_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pa_y, pc_y, pc_z, msk0_723, msk0_725, \
                         msi_532, msi_561, msi_562, msk1_723, msk1_725, nsi_728, \
                         nsi_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * msi_532[k]
                   + f_3 * pc_z[k] * nsi_728[k];

        t_939[k] = pa_y[k] * msk0_723[k]
                   + f_14 * msi_561[k]
                   - f_12 * pc_y[k] * msk1_723[k];

        t_940[k] = f_13 * msi_562[k]
                   + f_3 * pc_y[k] * nsi_730[k];

        t_941[k] = pa_y[k] * msk0_725[k]
                   - f_12 * pc_y[k] * msk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_y, pc_y, pc_z, msk0_726, msk0_729, \
                         msi_535, msi_563, msi_565, msk1_726, msk1_729, nsi_731, \
                         nsi_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pa_y[k] * msk0_726[k]
                   + f_15 * msi_563[k]
                   - f_12 * pc_y[k] * msk1_726[k];

        t_943[k] = f_17 * msi_535[k]
                   + f_3 * pc_z[k] * nsi_731[k];

        t_944[k] = f_13 * msi_565[k]
                   + f_3 * pc_y[k] * nsi_733[k];

        t_945[k] = pa_y[k] * msk0_729[k]
                   - f_12 * pc_y[k] * msk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pa_y, pc_y, pc_z, msk0_730, msk0_732, msi_538, \
                         msi_566, msi_568, msk1_730, msk1_732, \
                         nsi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pa_y[k] * msk0_730[k]
                   + f_16 * msi_566[k]
                   - f_12 * pc_y[k] * msk1_730[k];

        t_947[k] = f_17 * msi_538[k]
                   + f_3 * pc_z[k] * nsi_734[k];

        t_948[k] = pa_y[k] * msk0_732[k]
                   + f_14 * msi_568[k]
                   - f_12 * pc_y[k] * msk1_732[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, msk0_734, msk0_735, \
                         msi_542, msi_569, msi_570, msk1_734, msk1_735, nsi_737, \
                         nsi_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_13 * msi_569[k]
                   + f_3 * pc_y[k] * nsi_737[k];

        t_950[k] = pa_y[k] * msk0_734[k]
                   - f_12 * pc_y[k] * msk1_734[k];

        t_951[k] = pa_y[k] * msk0_735[k]
                   + f_17 * msi_570[k]
                   - f_12 * pc_y[k] * msk1_735[k];

        t_952[k] = f_17 * msi_542[k]
                   + f_3 * pc_z[k] * nsi_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, msk0_737, msk0_738, msk0_740, \
                         msi_572, msi_573, msi_574, msk1_737, msk1_738, msk1_740, \
                         nsi_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pa_y[k] * msk0_737[k]
                   + f_15 * msi_572[k]
                   - f_12 * pc_y[k] * msk1_737[k];

        t_954[k] = pa_y[k] * msk0_738[k]
                   + f_14 * msi_573[k]
                   - f_12 * pc_y[k] * msk1_738[k];

        t_955[k] = f_13 * msi_574[k]
                   + f_3 * pc_y[k] * nsi_742[k];

        t_956[k] = pa_y[k] * msk0_740[k]
                   - f_12 * pc_y[k] * msk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, msi_749, msi_750, msi_751, \
                         msi_752, msi_753, nsi_749, nsi_750, nsi_751, nsi_752, \
                         nsi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_16 * msi_749[k]
                   + f_3 * pc_x[k] * nsi_749[k];

        t_958[k] = f_16 * msi_750[k]
                   + f_3 * pc_x[k] * nsi_750[k];

        t_959[k] = f_16 * msi_751[k]
                   + f_3 * pc_x[k] * nsi_751[k];

        t_960[k] = f_16 * msi_752[k]
                   + f_3 * pc_x[k] * nsi_752[k];

        t_961[k] = f_16 * msi_753[k]
                   + f_3 * pc_x[k] * nsi_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pc_x, pc_y, pc_z, msi_553, msi_581, \
                         msi_754, msi_755, nsh0_561, nsh1_561, nsi_749, nsi_754, \
                         nsi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_16 * msi_754[k]
                   + f_3 * pc_x[k] * nsi_754[k];

        t_963[k] = f_16 * msi_755[k]
                   + f_3 * pc_x[k] * nsi_755[k];

        t_964[k] = f_13 * msi_581[k]
                   + f_1 * nsh0_561[k]
                   - f_2 * nsh1_561[k]
                   + f_3 * pc_y[k] * nsi_749[k];

        t_965[k] = f_17 * msi_553[k]
                   + f_3 * pc_z[k] * nsi_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_y, msi_583, msi_584, msi_585, nsh0_563, \
                         nsh0_564, nsh0_565, nsh1_563, nsh1_564, nsh1_565, nsi_751, nsi_752, \
                         nsi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_13 * msi_583[k]
                   + f_10 * nsh0_563[k]
                   - f_11 * nsh1_563[k]
                   + f_3 * pc_y[k] * nsi_751[k];

        t_967[k] = f_13 * msi_584[k]
                   + f_8 * nsh0_564[k]
                   - f_9 * nsh1_564[k]
                   + f_3 * pc_y[k] * nsi_752[k];

        t_968[k] = f_13 * msi_585[k]
                   + f_6 * nsh0_565[k]
                   - f_7 * nsh1_565[k]
                   + f_3 * pc_y[k] * nsi_753[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_y, pc_y, msk0_755, msi_586, msi_587, \
                         msk1_755, nsh0_566, nsh1_566, nsi_754, \
                         nsi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_13 * msi_586[k]
                   + f_4 * nsh0_566[k]
                   - f_5 * nsh1_566[k]
                   + f_3 * pc_y[k] * nsi_754[k];

        t_970[k] = f_13 * msi_587[k]
                   + f_3 * pc_y[k] * nsi_755[k];

        t_971[k] = pa_y[k] * msk0_755[k]
                   - f_12 * pc_y[k] * msk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, pc_x, pc_y, pc_z, msi_560, \
                         msi_756, nsh0_567, nsh1_567, nsi_756, nsi_757, \
                         nsi_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_16 * msi_756[k]
                   + f_1 * nsh0_567[k]
                   - f_2 * nsh1_567[k]
                   + f_3 * pc_x[k] * nsi_756[k];

        t_973[k] = f_3 * pc_y[k] * nsi_756[k];

        t_974[k] = f_23 * msi_560[k]
                   + f_3 * pc_z[k] * nsi_756[k];

        t_975[k] = f_4 * nsh0_567[k]
                   - f_5 * nsh1_567[k]
                   + f_3 * pc_y[k] * nsi_757[k];

        t_976[k] = f_3 * pc_y[k] * nsi_758[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, pc_x, pc_y, msi_761, nsh0_568, nsh0_569, \
                         nsh0_572, nsh1_568, nsh1_569, nsh1_572, nsi_759, nsi_760, \
                         nsi_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_16 * msi_761[k]
                   + f_10 * nsh0_572[k]
                   - f_11 * nsh1_572[k]
                   + f_3 * pc_x[k] * nsi_761[k];

        t_978[k] = f_6 * nsh0_568[k]
                   - f_7 * nsh1_568[k]
                   + f_3 * pc_y[k] * nsi_759[k];

        t_979[k] = f_4 * nsh0_569[k]
                   - f_5 * nsh1_569[k]
                   + f_3 * pc_y[k] * nsi_760[k];

        t_980[k] = f_3 * pc_y[k] * nsi_761[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_x, pc_y, msi_765, nsh0_570, nsh0_571, \
                         nsh0_576, nsh1_570, nsh1_571, nsh1_576, nsi_762, nsi_763, \
                         nsi_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_16 * msi_765[k]
                   + f_8 * nsh0_576[k]
                   - f_9 * nsh1_576[k]
                   + f_3 * pc_x[k] * nsi_765[k];

        t_982[k] = f_8 * nsh0_570[k]
                   - f_9 * nsh1_570[k]
                   + f_3 * pc_y[k] * nsi_762[k];

        t_983[k] = f_6 * nsh0_571[k]
                   - f_7 * nsh1_571[k]
                   + f_3 * pc_y[k] * nsi_763[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_x, pc_y, msi_770, nsh0_572, nsh0_581, \
                         nsh1_572, nsh1_581, nsi_764, nsi_765, \
                         nsi_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_4 * nsh0_572[k]
                   - f_5 * nsh1_572[k]
                   + f_3 * pc_y[k] * nsi_764[k];

        t_985[k] = f_3 * pc_y[k] * nsi_765[k];

        t_986[k] = f_16 * msi_770[k]
                   + f_6 * nsh0_581[k]
                   - f_7 * nsh1_581[k]
                   + f_3 * pc_x[k] * nsi_770[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pc_y, nsh0_573, nsh0_574, nsh0_575, nsh1_573, \
                         nsh1_574, nsh1_575, nsi_766, nsi_767, \
                         nsi_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_10 * nsh0_573[k]
                   - f_11 * nsh1_573[k]
                   + f_3 * pc_y[k] * nsi_766[k];

        t_988[k] = f_8 * nsh0_574[k]
                   - f_9 * nsh1_574[k]
                   + f_3 * pc_y[k] * nsi_767[k];

        t_989[k] = f_6 * nsh0_575[k]
                   - f_7 * nsh1_575[k]
                   + f_3 * pc_y[k] * nsi_768[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, pc_x, pc_y, msi_776, msi_777, nsh0_576, \
                         nsh0_587, nsh1_576, nsh1_587, nsi_769, nsi_770, nsi_776, \
                         nsi_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = f_4 * nsh0_576[k]
                   - f_5 * nsh1_576[k]
                   + f_3 * pc_y[k] * nsi_769[k];

        t_991[k] = f_3 * pc_y[k] * nsi_770[k];

        t_992[k] = f_16 * msi_776[k]
                   + f_4 * nsh0_587[k]
                   - f_5 * nsh1_587[k]
                   + f_3 * pc_x[k] * nsi_776[k];

        t_993[k] = f_16 * msi_777[k]
                   + f_3 * pc_x[k] * nsi_777[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, t_998, pc_x, pc_y, msi_778, msi_779, \
                         msi_780, msi_781, nsi_776, nsi_778, nsi_779, nsi_780, \
                         nsi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_16 * msi_778[k]
                   + f_3 * pc_x[k] * nsi_778[k];

        t_995[k] = f_16 * msi_779[k]
                   + f_3 * pc_x[k] * nsi_779[k];

        t_996[k] = f_16 * msi_780[k]
                   + f_3 * pc_x[k] * nsi_780[k];

        t_997[k] = f_16 * msi_781[k]
                   + f_3 * pc_x[k] * nsi_781[k];

        t_998[k] = f_3 * pc_y[k] * nsi_776[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, pc_x, pc_y, msi_783, nsh0_582, nsh0_583, \
                         nsh1_582, nsh1_583, nsi_777, nsi_778, \
                         nsi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_16 * msi_783[k]
                   + f_3 * pc_x[k] * nsi_783[k];

        t_1000[k] = f_1 * nsh0_582[k]
                    - f_2 * nsh1_582[k]
                    + f_3 * pc_y[k] * nsi_777[k];

        t_1001[k] = f_19 * nsh0_583[k]
                    - f_20 * nsh1_583[k]
                    + f_3 * pc_y[k] * nsi_778[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, nsh0_584, nsh0_585, nsh0_586, nsh1_584, \
                         nsh1_585, nsh1_586, nsi_779, nsi_780, \
                         nsi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_10 * nsh0_584[k]
                    - f_11 * nsh1_584[k]
                    + f_3 * pc_y[k] * nsi_779[k];

        t_1003[k] = f_8 * nsh0_585[k]
                    - f_9 * nsh1_585[k]
                    + f_3 * pc_y[k] * nsi_780[k];

        t_1004[k] = f_6 * nsh0_586[k]
                    - f_7 * nsh1_586[k]
                    + f_3 * pc_y[k] * nsi_781[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, pc_x, pc_y, pc_z, msi_587, msi_784, \
                         nsh0_587, nsh0_588, nsh1_587, nsh1_588, nsi_782, nsi_783, \
                         nsi_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_4 * nsh0_587[k]
                    - f_5 * nsh1_587[k]
                    + f_3 * pc_y[k] * nsi_782[k];

        t_1006[k] = f_3 * pc_y[k] * nsi_783[k];

        t_1007[k] = f_23 * msi_587[k]
                    + f_1 * nsh0_587[k]
                    - f_2 * nsh1_587[k]
                    + f_3 * pc_z[k] * nsi_783[k];

        t_1008[k] = f_15 * msi_784[k]
                    + f_1 * nsh0_588[k]
                    - f_2 * nsh1_588[k]
                    + f_3 * pc_x[k] * nsi_784[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pc_x, pc_y, pc_z, msi_588, msi_787, \
                         nsh0_591, nsh1_591, nsi_784, nsi_785, \
                         nsi_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_22 * msi_588[k]
                    + f_3 * pc_y[k] * nsi_784[k];

        t_1010[k] = f_3 * pc_z[k] * nsi_784[k];

        t_1011[k] = f_15 * msi_787[k]
                    + f_10 * nsh0_591[k]
                    - f_11 * nsh1_591[k]
                    + f_3 * pc_x[k] * nsi_787[k];

        t_1012[k] = f_3 * pc_z[k] * nsi_785[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, pc_z, msi_790, nsh0_588, nsh0_594, \
                         nsh1_588, nsh1_594, nsi_786, nsi_787, \
                         nsi_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_4 * nsh0_588[k]
                    - f_5 * nsh1_588[k]
                    + f_3 * pc_z[k] * nsi_786[k];

        t_1014[k] = f_15 * msi_790[k]
                    + f_8 * nsh0_594[k]
                    - f_9 * nsh1_594[k]
                    + f_3 * pc_x[k] * nsi_790[k];

        t_1015[k] = f_3 * pc_z[k] * nsi_787[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, pc_z, msi_593, msi_794, \
                         nsh0_590, nsh0_598, nsh1_590, nsh1_598, nsi_789, nsi_790, \
                         nsi_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_22 * msi_593[k]
                    + f_3 * pc_y[k] * nsi_789[k];

        t_1017[k] = f_6 * nsh0_590[k]
                    - f_7 * nsh1_590[k]
                    + f_3 * pc_z[k] * nsi_789[k];

        t_1018[k] = f_15 * msi_794[k]
                    + f_6 * nsh0_598[k]
                    - f_7 * nsh1_598[k]
                    + f_3 * pc_x[k] * nsi_794[k];

        t_1019[k] = f_3 * pc_z[k] * nsi_790[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, pc_y, pc_z, msi_597, nsh0_591, nsh0_593, \
                         nsh1_591, nsh1_593, nsi_791, nsi_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * nsh0_591[k]
                    - f_5 * nsh1_591[k]
                    + f_3 * pc_z[k] * nsi_791[k];

        t_1021[k] = f_22 * msi_597[k]
                    + f_3 * pc_y[k] * nsi_793[k];

        t_1022[k] = f_8 * nsh0_593[k]
                    - f_9 * nsh1_593[k]
                    + f_3 * pc_z[k] * nsi_793[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, pc_x, pc_z, msi_799, nsh0_594, nsh0_603, \
                         nsh1_594, nsh1_603, nsi_794, nsi_795, \
                         nsi_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_15 * msi_799[k]
                    + f_4 * nsh0_603[k]
                    - f_5 * nsh1_603[k]
                    + f_3 * pc_x[k] * nsi_799[k];

        t_1024[k] = f_3 * pc_z[k] * nsi_794[k];

        t_1025[k] = f_4 * nsh0_594[k]
                    - f_5 * nsh1_594[k]
                    + f_3 * pc_z[k] * nsi_795[k];
    }

#pragma omp simd aligned(t_1026, t_1027, t_1028, t_1029, pc_x, pc_y, pc_z, msi_602, msi_805, \
                         nsh0_595, nsh0_597, nsh1_595, nsh1_597, nsi_796, nsi_798, \
                         nsi_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1026[k] = f_6 * nsh0_595[k]
                    - f_7 * nsh1_595[k]
                    + f_3 * pc_z[k] * nsi_796[k];

        t_1027[k] = f_22 * msi_602[k]
                    + f_3 * pc_y[k] * nsi_798[k];

        t_1028[k] = f_10 * nsh0_597[k]
                    - f_11 * nsh1_597[k]
                    + f_3 * pc_z[k] * nsi_798[k];

        t_1029[k] = f_15 * msi_805[k]
                    + f_3 * pc_x[k] * nsi_805[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, t_1034, pc_x, pc_z, msi_807, msi_808, \
                         msi_809, msi_810, nsi_799, nsi_807, nsi_808, nsi_809, \
                         nsi_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_3 * pc_z[k] * nsi_799[k];

        t_1031[k] = f_15 * msi_807[k]
                    + f_3 * pc_x[k] * nsi_807[k];

        t_1032[k] = f_15 * msi_808[k]
                    + f_3 * pc_x[k] * nsi_808[k];

        t_1033[k] = f_15 * msi_809[k]
                    + f_3 * pc_x[k] * nsi_809[k];

        t_1034[k] = f_15 * msi_810[k]
                    + f_3 * pc_x[k] * nsi_810[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, pc_x, pc_y, pc_z, msi_609, msi_811, \
                         nsh0_603, nsh1_603, nsi_805, nsi_806, \
                         nsi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_15 * msi_811[k]
                    + f_3 * pc_x[k] * nsi_811[k];

        t_1036[k] = f_22 * msi_609[k]
                    + f_1 * nsh0_603[k]
                    - f_2 * nsh1_603[k]
                    + f_3 * pc_y[k] * nsi_805[k];

        t_1037[k] = f_3 * pc_z[k] * nsi_805[k];

        t_1038[k] = f_4 * nsh0_603[k]
                    - f_5 * nsh1_603[k]
                    + f_3 * pc_z[k] * nsi_806[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_z, nsh0_604, nsh0_605, nsh0_606, nsh1_604, \
                         nsh1_605, nsh1_606, nsi_807, nsi_808, \
                         nsi_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_6 * nsh0_604[k]
                    - f_7 * nsh1_604[k]
                    + f_3 * pc_z[k] * nsi_807[k];

        t_1040[k] = f_8 * nsh0_605[k]
                    - f_9 * nsh1_605[k]
                    + f_3 * pc_z[k] * nsi_808[k];

        t_1041[k] = f_10 * nsh0_606[k]
                    - f_11 * nsh1_606[k]
                    + f_3 * pc_z[k] * nsi_809[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pa_z, pc_y, pc_z, msk0_756, msi_615, \
                         msi_616, msk1_756, nsh0_608, nsh1_608, nsi_811, \
                         nsi_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_22 * msi_615[k]
                    + f_3 * pc_y[k] * nsi_811[k];

        t_1043[k] = f_1 * nsh0_608[k]
                    - f_2 * nsh1_608[k]
                    + f_3 * pc_z[k] * nsi_811[k];

        t_1044[k] = pa_z[k] * msk0_756[k]
                    - f_12 * pc_z[k] * msk1_756[k];

        t_1045[k] = f_23 * msi_616[k]
                    + f_3 * pc_y[k] * nsi_812[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pa_z, pc_y, pc_z, msk0_759, msi_588, msi_618, \
                         msk1_759, nsi_812, nsi_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_13 * msi_588[k]
                    + f_3 * pc_z[k] * nsi_812[k];

        t_1047[k] = pa_z[k] * msk0_759[k]
                    - f_12 * pc_z[k] * msk1_759[k];

        t_1048[k] = f_23 * msi_618[k]
                    + f_3 * pc_y[k] * nsi_814[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msk0,
                                                          const size_t msi, const size_t msk1,
                                                          const size_t nsh0, const size_t nsh1,
                                                          const size_t nsi, const size_t ncols,
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
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_762 = buffer.data(msk0 + 762);
    const auto *msk0_766 = buffer.data(msk0 + 766);
    const auto *msk0_768 = buffer.data(msk0 + 768);
    const auto *msk0_771 = buffer.data(msk0 + 771);
    const auto *msk0_773 = buffer.data(msk0 + 773);
    const auto *msk0_774 = buffer.data(msk0 + 774);
    const auto *msk0_784 = buffer.data(msk0 + 784);

    const auto *msi_591 = buffer.data(msi + 591);
    const auto *msi_594 = buffer.data(msi + 594);
    const auto *msi_595 = buffer.data(msi + 595);
    const auto *msi_598 = buffer.data(msi + 598);
    const auto *msi_599 = buffer.data(msi + 599);
    const auto *msi_600 = buffer.data(msi + 600);
    const auto *msi_609 = buffer.data(msi + 609);
    const auto *msi_615 = buffer.data(msi + 615);
    const auto *msi_616 = buffer.data(msi + 616);
    const auto *msi_619 = buffer.data(msi + 619);
    const auto *msi_621 = buffer.data(msi + 621);
    const auto *msi_622 = buffer.data(msi + 622);
    const auto *msi_625 = buffer.data(msi + 625);
    const auto *msi_626 = buffer.data(msi + 626);
    const auto *msi_630 = buffer.data(msi + 630);
    const auto *msi_637 = buffer.data(msi + 637);
    const auto *msi_639 = buffer.data(msi + 639);
    const auto *msi_640 = buffer.data(msi + 640);
    const auto *msi_641 = buffer.data(msi + 641);
    const auto *msi_642 = buffer.data(msi + 642);
    const auto *msi_643 = buffer.data(msi + 643);
    const auto *msi_644 = buffer.data(msi + 644);
    const auto *msi_646 = buffer.data(msi + 646);
    const auto *msi_647 = buffer.data(msi + 647);
    const auto *msi_649 = buffer.data(msi + 649);
    const auto *msi_650 = buffer.data(msi + 650);
    const auto *msi_653 = buffer.data(msi + 653);
    const auto *msi_654 = buffer.data(msi + 654);
    const auto *msi_658 = buffer.data(msi + 658);
    const auto *msi_665 = buffer.data(msi + 665);
    const auto *msi_667 = buffer.data(msi + 667);
    const auto *msi_668 = buffer.data(msi + 668);
    const auto *msi_669 = buffer.data(msi + 669);
    const auto *msi_670 = buffer.data(msi + 670);
    const auto *msi_671 = buffer.data(msi + 671);
    const auto *msi_672 = buffer.data(msi + 672);
    const auto *msi_674 = buffer.data(msi + 674);
    const auto *msi_677 = buffer.data(msi + 677);
    const auto *msi_681 = buffer.data(msi + 681);
    const auto *msi_686 = buffer.data(msi + 686);
    const auto *msi_693 = buffer.data(msi + 693);
    const auto *msi_695 = buffer.data(msi + 695);
    const auto *msi_696 = buffer.data(msi + 696);
    const auto *msi_697 = buffer.data(msi + 697);
    const auto *msi_698 = buffer.data(msi + 698);
    const auto *msi_699 = buffer.data(msi + 699);
    const auto *msi_817 = buffer.data(msi + 817);
    const auto *msi_821 = buffer.data(msi + 821);
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

    const auto *msk1_762 = buffer.data(msk1 + 762);
    const auto *msk1_766 = buffer.data(msk1 + 766);
    const auto *msk1_768 = buffer.data(msk1 + 768);
    const auto *msk1_771 = buffer.data(msk1 + 771);
    const auto *msk1_773 = buffer.data(msk1 + 773);
    const auto *msk1_774 = buffer.data(msk1 + 774);
    const auto *msk1_784 = buffer.data(msk1 + 784);

    const auto *nsh0_614 = buffer.data(nsh0 + 614);
    const auto *nsh0_618 = buffer.data(nsh0 + 618);
    const auto *nsh0_623 = buffer.data(nsh0 + 623);
    const auto *nsh0_626 = buffer.data(nsh0 + 626);
    const auto *nsh0_627 = buffer.data(nsh0 + 627);
    const auto *nsh0_628 = buffer.data(nsh0 + 628);
    const auto *nsh0_629 = buffer.data(nsh0 + 629);
    const auto *nsh0_630 = buffer.data(nsh0 + 630);
    const auto *nsh0_633 = buffer.data(nsh0 + 633);
    const auto *nsh0_635 = buffer.data(nsh0 + 635);
    const auto *nsh0_636 = buffer.data(nsh0 + 636);
    const auto *nsh0_639 = buffer.data(nsh0 + 639);
    const auto *nsh0_640 = buffer.data(nsh0 + 640);
    const auto *nsh0_642 = buffer.data(nsh0 + 642);
    const auto *nsh0_644 = buffer.data(nsh0 + 644);
    const auto *nsh0_645 = buffer.data(nsh0 + 645);
    const auto *nsh0_647 = buffer.data(nsh0 + 647);
    const auto *nsh0_648 = buffer.data(nsh0 + 648);
    const auto *nsh0_649 = buffer.data(nsh0 + 649);
    const auto *nsh0_650 = buffer.data(nsh0 + 650);
    const auto *nsh0_651 = buffer.data(nsh0 + 651);
    const auto *nsh0_654 = buffer.data(nsh0 + 654);
    const auto *nsh0_656 = buffer.data(nsh0 + 656);
    const auto *nsh0_657 = buffer.data(nsh0 + 657);
    const auto *nsh0_660 = buffer.data(nsh0 + 660);
    const auto *nsh0_661 = buffer.data(nsh0 + 661);
    const auto *nsh0_663 = buffer.data(nsh0 + 663);
    const auto *nsh0_665 = buffer.data(nsh0 + 665);
    const auto *nsh0_666 = buffer.data(nsh0 + 666);
    const auto *nsh0_668 = buffer.data(nsh0 + 668);
    const auto *nsh0_669 = buffer.data(nsh0 + 669);
    const auto *nsh0_670 = buffer.data(nsh0 + 670);
    const auto *nsh0_671 = buffer.data(nsh0 + 671);
    const auto *nsh0_672 = buffer.data(nsh0 + 672);

    const auto *nsh1_614 = buffer.data(nsh1 + 614);
    const auto *nsh1_618 = buffer.data(nsh1 + 618);
    const auto *nsh1_623 = buffer.data(nsh1 + 623);
    const auto *nsh1_626 = buffer.data(nsh1 + 626);
    const auto *nsh1_627 = buffer.data(nsh1 + 627);
    const auto *nsh1_628 = buffer.data(nsh1 + 628);
    const auto *nsh1_629 = buffer.data(nsh1 + 629);
    const auto *nsh1_630 = buffer.data(nsh1 + 630);
    const auto *nsh1_633 = buffer.data(nsh1 + 633);
    const auto *nsh1_635 = buffer.data(nsh1 + 635);
    const auto *nsh1_636 = buffer.data(nsh1 + 636);
    const auto *nsh1_639 = buffer.data(nsh1 + 639);
    const auto *nsh1_640 = buffer.data(nsh1 + 640);
    const auto *nsh1_642 = buffer.data(nsh1 + 642);
    const auto *nsh1_644 = buffer.data(nsh1 + 644);
    const auto *nsh1_645 = buffer.data(nsh1 + 645);
    const auto *nsh1_647 = buffer.data(nsh1 + 647);
    const auto *nsh1_648 = buffer.data(nsh1 + 648);
    const auto *nsh1_649 = buffer.data(nsh1 + 649);
    const auto *nsh1_650 = buffer.data(nsh1 + 650);
    const auto *nsh1_651 = buffer.data(nsh1 + 651);
    const auto *nsh1_654 = buffer.data(nsh1 + 654);
    const auto *nsh1_656 = buffer.data(nsh1 + 656);
    const auto *nsh1_657 = buffer.data(nsh1 + 657);
    const auto *nsh1_660 = buffer.data(nsh1 + 660);
    const auto *nsh1_661 = buffer.data(nsh1 + 661);
    const auto *nsh1_663 = buffer.data(nsh1 + 663);
    const auto *nsh1_665 = buffer.data(nsh1 + 665);
    const auto *nsh1_666 = buffer.data(nsh1 + 666);
    const auto *nsh1_668 = buffer.data(nsh1 + 668);
    const auto *nsh1_669 = buffer.data(nsh1 + 669);
    const auto *nsh1_670 = buffer.data(nsh1 + 670);
    const auto *nsh1_671 = buffer.data(nsh1 + 671);
    const auto *nsh1_672 = buffer.data(nsh1 + 672);

    const auto *nsi_815 = buffer.data(nsi + 815);
    const auto *nsi_817 = buffer.data(nsi + 817);
    const auto *nsi_818 = buffer.data(nsi + 818);
    const auto *nsi_821 = buffer.data(nsi + 821);
    const auto *nsi_822 = buffer.data(nsi + 822);
    const auto *nsi_826 = buffer.data(nsi + 826);
    const auto *nsi_832 = buffer.data(nsi + 832);
    const auto *nsi_833 = buffer.data(nsi + 833);
    const auto *nsi_834 = buffer.data(nsi + 834);
    const auto *nsi_835 = buffer.data(nsi + 835);
    const auto *nsi_836 = buffer.data(nsi + 836);
    const auto *nsi_837 = buffer.data(nsi + 837);
    const auto *nsi_838 = buffer.data(nsi + 838);
    const auto *nsi_839 = buffer.data(nsi + 839);
    const auto *nsi_840 = buffer.data(nsi + 840);
    const auto *nsi_842 = buffer.data(nsi + 842);
    const auto *nsi_843 = buffer.data(nsi + 843);
    const auto *nsi_845 = buffer.data(nsi + 845);
    const auto *nsi_846 = buffer.data(nsi + 846);
    const auto *nsi_849 = buffer.data(nsi + 849);
    const auto *nsi_850 = buffer.data(nsi + 850);
    const auto *nsi_852 = buffer.data(nsi + 852);
    const auto *nsi_854 = buffer.data(nsi + 854);
    const auto *nsi_855 = buffer.data(nsi + 855);
    const auto *nsi_857 = buffer.data(nsi + 857);
    const auto *nsi_858 = buffer.data(nsi + 858);
    const auto *nsi_860 = buffer.data(nsi + 860);
    const auto *nsi_861 = buffer.data(nsi + 861);
    const auto *nsi_862 = buffer.data(nsi + 862);
    const auto *nsi_863 = buffer.data(nsi + 863);
    const auto *nsi_864 = buffer.data(nsi + 864);
    const auto *nsi_865 = buffer.data(nsi + 865);
    const auto *nsi_866 = buffer.data(nsi + 866);
    const auto *nsi_867 = buffer.data(nsi + 867);
    const auto *nsi_868 = buffer.data(nsi + 868);
    const auto *nsi_870 = buffer.data(nsi + 870);
    const auto *nsi_871 = buffer.data(nsi + 871);
    const auto *nsi_873 = buffer.data(nsi + 873);
    const auto *nsi_874 = buffer.data(nsi + 874);
    const auto *nsi_877 = buffer.data(nsi + 877);
    const auto *nsi_878 = buffer.data(nsi + 878);
    const auto *nsi_880 = buffer.data(nsi + 880);
    const auto *nsi_882 = buffer.data(nsi + 882);
    const auto *nsi_883 = buffer.data(nsi + 883);
    const auto *nsi_885 = buffer.data(nsi + 885);
    const auto *nsi_886 = buffer.data(nsi + 886);
    const auto *nsi_888 = buffer.data(nsi + 888);
    const auto *nsi_889 = buffer.data(nsi + 889);
    const auto *nsi_890 = buffer.data(nsi + 890);
    const auto *nsi_891 = buffer.data(nsi + 891);
    const auto *nsi_892 = buffer.data(nsi + 892);
    const auto *nsi_893 = buffer.data(nsi + 893);
    const auto *nsi_894 = buffer.data(nsi + 894);
    const auto *nsi_895 = buffer.data(nsi + 895);
    const auto *nsi_896 = buffer.data(nsi + 896);

#pragma omp simd aligned(t_1049, t_1050, t_1051, pa_z, pc_x, pc_z, msk0_762, msi_591, msi_817, \
                         msk1_762, nsh0_614, nsh1_614, nsi_815, \
                         nsi_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_15 * msi_817[k]
                    + f_10 * nsh0_614[k]
                    - f_11 * nsh1_614[k]
                    + f_3 * pc_x[k] * nsi_817[k];

        t_1050[k] = pa_z[k] * msk0_762[k]
                    - f_12 * pc_z[k] * msk1_762[k];

        t_1051[k] = f_13 * msi_591[k]
                    + f_3 * pc_z[k] * nsi_815[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pa_z, pc_x, pc_y, pc_z, msk0_766, msi_621, \
                         msi_821, msk1_766, nsh0_618, nsh1_618, nsi_817, \
                         nsi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_23 * msi_621[k]
                    + f_3 * pc_y[k] * nsi_817[k];

        t_1053[k] = f_15 * msi_821[k]
                    + f_8 * nsh0_618[k]
                    - f_9 * nsh1_618[k]
                    + f_3 * pc_x[k] * nsi_821[k];

        t_1054[k] = pa_z[k] * msk0_766[k]
                    - f_12 * pc_z[k] * msk1_766[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pa_z, pc_y, pc_z, msk0_768, msi_594, msi_595, \
                         msi_625, msk1_768, nsi_818, nsi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_13 * msi_594[k]
                    + f_3 * pc_z[k] * nsi_818[k];

        t_1056[k] = pa_z[k] * msk0_768[k]
                    + f_14 * msi_595[k]
                    - f_12 * pc_z[k] * msk1_768[k];

        t_1057[k] = f_23 * msi_625[k]
                    + f_3 * pc_y[k] * nsi_821[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pa_z, pc_x, pc_z, msk0_771, msi_598, msi_826, \
                         msk1_771, nsh0_623, nsh1_623, nsi_822, \
                         nsi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_15 * msi_826[k]
                    + f_6 * nsh0_623[k]
                    - f_7 * nsh1_623[k]
                    + f_3 * pc_x[k] * nsi_826[k];

        t_1059[k] = pa_z[k] * msk0_771[k]
                    - f_12 * pc_z[k] * msk1_771[k];

        t_1060[k] = f_13 * msi_598[k]
                    + f_3 * pc_z[k] * nsi_822[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pa_z, pc_y, pc_z, msk0_773, msk0_774, \
                         msi_599, msi_600, msi_630, msk1_773, msk1_774, \
                         nsi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pa_z[k] * msk0_773[k]
                    + f_14 * msi_599[k]
                    - f_12 * pc_z[k] * msk1_773[k];

        t_1062[k] = pa_z[k] * msk0_774[k]
                    + f_15 * msi_600[k]
                    - f_12 * pc_z[k] * msk1_774[k];

        t_1063[k] = f_23 * msi_630[k]
                    + f_3 * pc_y[k] * nsi_826[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, pc_x, msi_832, msi_833, msi_834, \
                         msi_835, nsh0_629, nsh1_629, nsi_832, nsi_833, nsi_834, \
                         nsi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_15 * msi_832[k]
                    + f_4 * nsh0_629[k]
                    - f_5 * nsh1_629[k]
                    + f_3 * pc_x[k] * nsi_832[k];

        t_1065[k] = f_15 * msi_833[k]
                    + f_3 * pc_x[k] * nsi_833[k];

        t_1066[k] = f_15 * msi_834[k]
                    + f_3 * pc_x[k] * nsi_834[k];

        t_1067[k] = f_15 * msi_835[k]
                    + f_3 * pc_x[k] * nsi_835[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, pc_x, msi_836, msi_837, msi_838, \
                         msi_839, nsi_836, nsi_837, nsi_838, nsi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_15 * msi_836[k]
                    + f_3 * pc_x[k] * nsi_836[k];

        t_1069[k] = f_15 * msi_837[k]
                    + f_3 * pc_x[k] * nsi_837[k];

        t_1070[k] = f_15 * msi_838[k]
                    + f_3 * pc_x[k] * nsi_838[k];

        t_1071[k] = f_15 * msi_839[k]
                    + f_3 * pc_x[k] * nsi_839[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, pa_z, pc_y, pc_z, msk0_784, msi_609, msi_639, \
                         msk1_784, nsh0_626, nsh1_626, nsi_833, \
                         nsi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = pa_z[k] * msk0_784[k]
                    - f_12 * pc_z[k] * msk1_784[k];

        t_1073[k] = f_13 * msi_609[k]
                    + f_3 * pc_z[k] * nsi_833[k];

        t_1074[k] = f_23 * msi_639[k]
                    + f_10 * nsh0_626[k]
                    - f_11 * nsh1_626[k]
                    + f_3 * pc_y[k] * nsi_835[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, pc_y, msi_640, msi_641, msi_642, nsh0_627, \
                         nsh0_628, nsh0_629, nsh1_627, nsh1_628, nsh1_629, nsi_836, nsi_837, \
                         nsi_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_23 * msi_640[k]
                    + f_8 * nsh0_627[k]
                    - f_9 * nsh1_627[k]
                    + f_3 * pc_y[k] * nsi_836[k];

        t_1076[k] = f_23 * msi_641[k]
                    + f_6 * nsh0_628[k]
                    - f_7 * nsh1_628[k]
                    + f_3 * pc_y[k] * nsi_837[k];

        t_1077[k] = f_23 * msi_642[k]
                    + f_4 * nsh0_629[k]
                    - f_5 * nsh1_629[k]
                    + f_3 * pc_y[k] * nsi_838[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, pc_x, pc_y, pc_z, msi_615, msi_643, msi_840, \
                         nsh0_629, nsh0_630, nsh1_629, nsh1_630, nsi_839, \
                         nsi_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_23 * msi_643[k]
                    + f_3 * pc_y[k] * nsi_839[k];

        t_1079[k] = f_13 * msi_615[k]
                    + f_1 * nsh0_629[k]
                    - f_2 * nsh1_629[k]
                    + f_3 * pc_z[k] * nsi_839[k];

        t_1080[k] = f_15 * msi_840[k]
                    + f_1 * nsh0_630[k]
                    - f_2 * nsh1_630[k]
                    + f_3 * pc_x[k] * nsi_840[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, msi_616, msi_644, \
                         msi_646, msi_843, nsh0_633, nsh1_633, nsi_840, nsi_842, \
                         nsi_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_17 * msi_644[k]
                    + f_3 * pc_y[k] * nsi_840[k];

        t_1082[k] = f_14 * msi_616[k]
                    + f_3 * pc_z[k] * nsi_840[k];

        t_1083[k] = f_15 * msi_843[k]
                    + f_10 * nsh0_633[k]
                    - f_11 * nsh1_633[k]
                    + f_3 * pc_x[k] * nsi_843[k];

        t_1084[k] = f_17 * msi_646[k]
                    + f_3 * pc_y[k] * nsi_842[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, msi_619, msi_845, msi_846, \
                         nsh0_635, nsh0_636, nsh1_635, nsh1_636, nsi_843, nsi_845, \
                         nsi_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_15 * msi_845[k]
                    + f_10 * nsh0_635[k]
                    - f_11 * nsh1_635[k]
                    + f_3 * pc_x[k] * nsi_845[k];

        t_1086[k] = f_15 * msi_846[k]
                    + f_8 * nsh0_636[k]
                    - f_9 * nsh1_636[k]
                    + f_3 * pc_x[k] * nsi_846[k];

        t_1087[k] = f_14 * msi_619[k]
                    + f_3 * pc_z[k] * nsi_843[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, msi_649, msi_849, msi_850, \
                         nsh0_639, nsh0_640, nsh1_639, nsh1_640, nsi_845, nsi_849, \
                         nsi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * msi_649[k]
                    + f_3 * pc_y[k] * nsi_845[k];

        t_1089[k] = f_15 * msi_849[k]
                    + f_8 * nsh0_639[k]
                    - f_9 * nsh1_639[k]
                    + f_3 * pc_x[k] * nsi_849[k];

        t_1090[k] = f_15 * msi_850[k]
                    + f_6 * nsh0_640[k]
                    - f_7 * nsh1_640[k]
                    + f_3 * pc_x[k] * nsi_850[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, msi_622, msi_653, msi_852, \
                         nsh0_642, nsh1_642, nsi_846, nsi_849, \
                         nsi_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_14 * msi_622[k]
                    + f_3 * pc_z[k] * nsi_846[k];

        t_1092[k] = f_15 * msi_852[k]
                    + f_6 * nsh0_642[k]
                    - f_7 * nsh1_642[k]
                    + f_3 * pc_x[k] * nsi_852[k];

        t_1093[k] = f_17 * msi_653[k]
                    + f_3 * pc_y[k] * nsi_849[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, msi_626, msi_854, msi_855, \
                         nsh0_644, nsh0_645, nsh1_644, nsh1_645, nsi_850, nsi_854, \
                         nsi_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_15 * msi_854[k]
                    + f_6 * nsh0_644[k]
                    - f_7 * nsh1_644[k]
                    + f_3 * pc_x[k] * nsi_854[k];

        t_1095[k] = f_15 * msi_855[k]
                    + f_4 * nsh0_645[k]
                    - f_5 * nsh1_645[k]
                    + f_3 * pc_x[k] * nsi_855[k];

        t_1096[k] = f_14 * msi_626[k]
                    + f_3 * pc_z[k] * nsi_850[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, msi_658, msi_857, msi_858, \
                         nsh0_647, nsh0_648, nsh1_647, nsh1_648, nsi_854, nsi_857, \
                         nsi_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_15 * msi_857[k]
                    + f_4 * nsh0_647[k]
                    - f_5 * nsh1_647[k]
                    + f_3 * pc_x[k] * nsi_857[k];

        t_1098[k] = f_15 * msi_858[k]
                    + f_4 * nsh0_648[k]
                    - f_5 * nsh1_648[k]
                    + f_3 * pc_x[k] * nsi_858[k];

        t_1099[k] = f_17 * msi_658[k]
                    + f_3 * pc_y[k] * nsi_854[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, msi_860, msi_861, msi_862, \
                         msi_863, nsh0_650, nsh1_650, nsi_860, nsi_861, nsi_862, \
                         nsi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_15 * msi_860[k]
                    + f_4 * nsh0_650[k]
                    - f_5 * nsh1_650[k]
                    + f_3 * pc_x[k] * nsi_860[k];

        t_1101[k] = f_15 * msi_861[k]
                    + f_3 * pc_x[k] * nsi_861[k];

        t_1102[k] = f_15 * msi_862[k]
                    + f_3 * pc_x[k] * nsi_862[k];

        t_1103[k] = f_15 * msi_863[k]
                    + f_3 * pc_x[k] * nsi_863[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, msi_864, msi_865, msi_866, \
                         msi_867, nsi_864, nsi_865, nsi_866, nsi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_15 * msi_864[k]
                    + f_3 * pc_x[k] * nsi_864[k];

        t_1105[k] = f_15 * msi_865[k]
                    + f_3 * pc_x[k] * nsi_865[k];

        t_1106[k] = f_15 * msi_866[k]
                    + f_3 * pc_x[k] * nsi_866[k];

        t_1107[k] = f_15 * msi_867[k]
                    + f_3 * pc_x[k] * nsi_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, msi_637, msi_665, msi_667, \
                         nsh0_645, nsh0_647, nsh1_645, nsh1_647, nsi_861, \
                         nsi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * msi_665[k]
                    + f_1 * nsh0_645[k]
                    - f_2 * nsh1_645[k]
                    + f_3 * pc_y[k] * nsi_861[k];

        t_1109[k] = f_14 * msi_637[k]
                    + f_3 * pc_z[k] * nsi_861[k];

        t_1110[k] = f_17 * msi_667[k]
                    + f_10 * nsh0_647[k]
                    - f_11 * nsh1_647[k]
                    + f_3 * pc_y[k] * nsi_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pc_y, msi_668, msi_669, msi_670, nsh0_648, \
                         nsh0_649, nsh0_650, nsh1_648, nsh1_649, nsh1_650, nsi_864, nsi_865, \
                         nsi_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_17 * msi_668[k]
                    + f_8 * nsh0_648[k]
                    - f_9 * nsh1_648[k]
                    + f_3 * pc_y[k] * nsi_864[k];

        t_1112[k] = f_17 * msi_669[k]
                    + f_6 * nsh0_649[k]
                    - f_7 * nsh1_649[k]
                    + f_3 * pc_y[k] * nsi_865[k];

        t_1113[k] = f_17 * msi_670[k]
                    + f_4 * nsh0_650[k]
                    - f_5 * nsh1_650[k]
                    + f_3 * pc_y[k] * nsi_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_x, pc_y, pc_z, msi_643, msi_671, msi_868, \
                         nsh0_650, nsh0_651, nsh1_650, nsh1_651, nsi_867, \
                         nsi_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * msi_671[k]
                    + f_3 * pc_y[k] * nsi_867[k];

        t_1115[k] = f_14 * msi_643[k]
                    + f_1 * nsh0_650[k]
                    - f_2 * nsh1_650[k]
                    + f_3 * pc_z[k] * nsi_867[k];

        t_1116[k] = f_15 * msi_868[k]
                    + f_1 * nsh0_651[k]
                    - f_2 * nsh1_651[k]
                    + f_3 * pc_x[k] * nsi_868[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, pc_x, pc_y, pc_z, msi_644, msi_672, \
                         msi_674, msi_871, nsh0_654, nsh1_654, nsi_868, nsi_870, \
                         nsi_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_16 * msi_672[k]
                    + f_3 * pc_y[k] * nsi_868[k];

        t_1118[k] = f_15 * msi_644[k]
                    + f_3 * pc_z[k] * nsi_868[k];

        t_1119[k] = f_15 * msi_871[k]
                    + f_10 * nsh0_654[k]
                    - f_11 * nsh1_654[k]
                    + f_3 * pc_x[k] * nsi_871[k];

        t_1120[k] = f_16 * msi_674[k]
                    + f_3 * pc_y[k] * nsi_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_x, pc_z, msi_647, msi_873, msi_874, \
                         nsh0_656, nsh0_657, nsh1_656, nsh1_657, nsi_871, nsi_873, \
                         nsi_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_15 * msi_873[k]
                    + f_10 * nsh0_656[k]
                    - f_11 * nsh1_656[k]
                    + f_3 * pc_x[k] * nsi_873[k];

        t_1122[k] = f_15 * msi_874[k]
                    + f_8 * nsh0_657[k]
                    - f_9 * nsh1_657[k]
                    + f_3 * pc_x[k] * nsi_874[k];

        t_1123[k] = f_15 * msi_647[k]
                    + f_3 * pc_z[k] * nsi_871[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, msi_677, msi_877, msi_878, \
                         nsh0_660, nsh0_661, nsh1_660, nsh1_661, nsi_873, nsi_877, \
                         nsi_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_16 * msi_677[k]
                    + f_3 * pc_y[k] * nsi_873[k];

        t_1125[k] = f_15 * msi_877[k]
                    + f_8 * nsh0_660[k]
                    - f_9 * nsh1_660[k]
                    + f_3 * pc_x[k] * nsi_877[k];

        t_1126[k] = f_15 * msi_878[k]
                    + f_6 * nsh0_661[k]
                    - f_7 * nsh1_661[k]
                    + f_3 * pc_x[k] * nsi_878[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, msi_650, msi_681, msi_880, \
                         nsh0_663, nsh1_663, nsi_874, nsi_877, \
                         nsi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_15 * msi_650[k]
                    + f_3 * pc_z[k] * nsi_874[k];

        t_1128[k] = f_15 * msi_880[k]
                    + f_6 * nsh0_663[k]
                    - f_7 * nsh1_663[k]
                    + f_3 * pc_x[k] * nsi_880[k];

        t_1129[k] = f_16 * msi_681[k]
                    + f_3 * pc_y[k] * nsi_877[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, msi_654, msi_882, msi_883, \
                         nsh0_665, nsh0_666, nsh1_665, nsh1_666, nsi_878, nsi_882, \
                         nsi_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_15 * msi_882[k]
                    + f_6 * nsh0_665[k]
                    - f_7 * nsh1_665[k]
                    + f_3 * pc_x[k] * nsi_882[k];

        t_1131[k] = f_15 * msi_883[k]
                    + f_4 * nsh0_666[k]
                    - f_5 * nsh1_666[k]
                    + f_3 * pc_x[k] * nsi_883[k];

        t_1132[k] = f_15 * msi_654[k]
                    + f_3 * pc_z[k] * nsi_878[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, msi_686, msi_885, msi_886, \
                         nsh0_668, nsh0_669, nsh1_668, nsh1_669, nsi_882, nsi_885, \
                         nsi_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_15 * msi_885[k]
                    + f_4 * nsh0_668[k]
                    - f_5 * nsh1_668[k]
                    + f_3 * pc_x[k] * nsi_885[k];

        t_1134[k] = f_15 * msi_886[k]
                    + f_4 * nsh0_669[k]
                    - f_5 * nsh1_669[k]
                    + f_3 * pc_x[k] * nsi_886[k];

        t_1135[k] = f_16 * msi_686[k]
                    + f_3 * pc_y[k] * nsi_882[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, pc_x, msi_888, msi_889, msi_890, \
                         msi_891, nsh0_671, nsh1_671, nsi_888, nsi_889, nsi_890, \
                         nsi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_15 * msi_888[k]
                    + f_4 * nsh0_671[k]
                    - f_5 * nsh1_671[k]
                    + f_3 * pc_x[k] * nsi_888[k];

        t_1137[k] = f_15 * msi_889[k]
                    + f_3 * pc_x[k] * nsi_889[k];

        t_1138[k] = f_15 * msi_890[k]
                    + f_3 * pc_x[k] * nsi_890[k];

        t_1139[k] = f_15 * msi_891[k]
                    + f_3 * pc_x[k] * nsi_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pc_x, msi_892, msi_893, msi_894, \
                         msi_895, nsi_892, nsi_893, nsi_894, nsi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_15 * msi_892[k]
                    + f_3 * pc_x[k] * nsi_892[k];

        t_1141[k] = f_15 * msi_893[k]
                    + f_3 * pc_x[k] * nsi_893[k];

        t_1142[k] = f_15 * msi_894[k]
                    + f_3 * pc_x[k] * nsi_894[k];

        t_1143[k] = f_15 * msi_895[k]
                    + f_3 * pc_x[k] * nsi_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pc_y, pc_z, msi_665, msi_693, msi_695, \
                         nsh0_666, nsh0_668, nsh1_666, nsh1_668, nsi_889, \
                         nsi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_16 * msi_693[k]
                    + f_1 * nsh0_666[k]
                    - f_2 * nsh1_666[k]
                    + f_3 * pc_y[k] * nsi_889[k];

        t_1145[k] = f_15 * msi_665[k]
                    + f_3 * pc_z[k] * nsi_889[k];

        t_1146[k] = f_16 * msi_695[k]
                    + f_10 * nsh0_668[k]
                    - f_11 * nsh1_668[k]
                    + f_3 * pc_y[k] * nsi_891[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pc_y, msi_696, msi_697, msi_698, nsh0_669, \
                         nsh0_670, nsh0_671, nsh1_669, nsh1_670, nsh1_671, nsi_892, nsi_893, \
                         nsi_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_16 * msi_696[k]
                    + f_8 * nsh0_669[k]
                    - f_9 * nsh1_669[k]
                    + f_3 * pc_y[k] * nsi_892[k];

        t_1148[k] = f_16 * msi_697[k]
                    + f_6 * nsh0_670[k]
                    - f_7 * nsh1_670[k]
                    + f_3 * pc_y[k] * nsi_893[k];

        t_1149[k] = f_16 * msi_698[k]
                    + f_4 * nsh0_671[k]
                    - f_5 * nsh1_671[k]
                    + f_3 * pc_y[k] * nsi_894[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, msi_671, msi_699, msi_896, \
                         nsh0_671, nsh0_672, nsh1_671, nsh1_672, nsi_895, \
                         nsi_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_16 * msi_699[k]
                    + f_3 * pc_y[k] * nsi_895[k];

        t_1151[k] = f_15 * msi_671[k]
                    + f_1 * nsh0_671[k]
                    - f_2 * nsh1_671[k]
                    + f_3 * pc_z[k] * nsi_895[k];

        t_1152[k] = f_15 * msi_896[k]
                    + f_1 * nsh0_672[k]
                    - f_2 * nsh1_672[k]
                    + f_3 * pc_x[k] * nsi_896[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsh0, const size_t nsh1,
                                                           const size_t nsi, const size_t ncols,
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
    const auto f_23 = 3.0 / q;

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

    const auto *msk0_972 = buffer.data(msk0 + 972);
    const auto *msk0_975 = buffer.data(msk0 + 975);
    const auto *msk0_977 = buffer.data(msk0 + 977);
    const auto *msk0_978 = buffer.data(msk0 + 978);
    const auto *msk0_981 = buffer.data(msk0 + 981);
    const auto *msk0_982 = buffer.data(msk0 + 982);
    const auto *msk0_984 = buffer.data(msk0 + 984);
    const auto *msk0_986 = buffer.data(msk0 + 986);
    const auto *msk0_987 = buffer.data(msk0 + 987);
    const auto *msk0_989 = buffer.data(msk0 + 989);
    const auto *msk0_990 = buffer.data(msk0 + 990);
    const auto *msk0_992 = buffer.data(msk0 + 992);
    const auto *msk0_1007 = buffer.data(msk0 + 1007);

    const auto *msi_672 = buffer.data(msi + 672);
    const auto *msi_675 = buffer.data(msi + 675);
    const auto *msi_678 = buffer.data(msi + 678);
    const auto *msi_682 = buffer.data(msi + 682);
    const auto *msi_693 = buffer.data(msi + 693);
    const auto *msi_699 = buffer.data(msi + 699);
    const auto *msi_700 = buffer.data(msi + 700);
    const auto *msi_702 = buffer.data(msi + 702);
    const auto *msi_703 = buffer.data(msi + 703);
    const auto *msi_705 = buffer.data(msi + 705);
    const auto *msi_706 = buffer.data(msi + 706);
    const auto *msi_709 = buffer.data(msi + 709);
    const auto *msi_710 = buffer.data(msi + 710);
    const auto *msi_714 = buffer.data(msi + 714);
    const auto *msi_721 = buffer.data(msi + 721);
    const auto *msi_723 = buffer.data(msi + 723);
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
    const auto *msi_751 = buffer.data(msi + 751);
    const auto *msi_752 = buffer.data(msi + 752);
    const auto *msi_753 = buffer.data(msi + 753);
    const auto *msi_754 = buffer.data(msi + 754);
    const auto *msi_755 = buffer.data(msi + 755);
    const auto *msi_756 = buffer.data(msi + 756);
    const auto *msi_757 = buffer.data(msi + 757);
    const auto *msi_758 = buffer.data(msi + 758);
    const auto *msi_759 = buffer.data(msi + 759);
    const auto *msi_761 = buffer.data(msi + 761);
    const auto *msi_762 = buffer.data(msi + 762);
    const auto *msi_764 = buffer.data(msi + 764);
    const auto *msi_765 = buffer.data(msi + 765);
    const auto *msi_766 = buffer.data(msi + 766);
    const auto *msi_768 = buffer.data(msi + 768);
    const auto *msi_769 = buffer.data(msi + 769);
    const auto *msi_770 = buffer.data(msi + 770);
    const auto *msi_777 = buffer.data(msi + 777);
    const auto *msi_779 = buffer.data(msi + 779);
    const auto *msi_780 = buffer.data(msi + 780);
    const auto *msi_781 = buffer.data(msi + 781);
    const auto *msi_782 = buffer.data(msi + 782);
    const auto *msi_783 = buffer.data(msi + 783);
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
    const auto *msi_973 = buffer.data(msi + 973);
    const auto *msi_974 = buffer.data(msi + 974);
    const auto *msi_975 = buffer.data(msi + 975);
    const auto *msi_976 = buffer.data(msi + 976);
    const auto *msi_977 = buffer.data(msi + 977);
    const auto *msi_978 = buffer.data(msi + 978);
    const auto *msi_979 = buffer.data(msi + 979);

    const auto *msk1_972 = buffer.data(msk1 + 972);
    const auto *msk1_975 = buffer.data(msk1 + 975);
    const auto *msk1_977 = buffer.data(msk1 + 977);
    const auto *msk1_978 = buffer.data(msk1 + 978);
    const auto *msk1_981 = buffer.data(msk1 + 981);
    const auto *msk1_982 = buffer.data(msk1 + 982);
    const auto *msk1_984 = buffer.data(msk1 + 984);
    const auto *msk1_986 = buffer.data(msk1 + 986);
    const auto *msk1_987 = buffer.data(msk1 + 987);
    const auto *msk1_989 = buffer.data(msk1 + 989);
    const auto *msk1_990 = buffer.data(msk1 + 990);
    const auto *msk1_992 = buffer.data(msk1 + 992);
    const auto *msk1_1007 = buffer.data(msk1 + 1007);

    const auto *nsh0_675 = buffer.data(nsh0 + 675);
    const auto *nsh0_677 = buffer.data(nsh0 + 677);
    const auto *nsh0_678 = buffer.data(nsh0 + 678);
    const auto *nsh0_681 = buffer.data(nsh0 + 681);
    const auto *nsh0_682 = buffer.data(nsh0 + 682);
    const auto *nsh0_684 = buffer.data(nsh0 + 684);
    const auto *nsh0_686 = buffer.data(nsh0 + 686);
    const auto *nsh0_687 = buffer.data(nsh0 + 687);
    const auto *nsh0_689 = buffer.data(nsh0 + 689);
    const auto *nsh0_690 = buffer.data(nsh0 + 690);
    const auto *nsh0_691 = buffer.data(nsh0 + 691);
    const auto *nsh0_692 = buffer.data(nsh0 + 692);
    const auto *nsh0_693 = buffer.data(nsh0 + 693);
    const auto *nsh0_696 = buffer.data(nsh0 + 696);
    const auto *nsh0_698 = buffer.data(nsh0 + 698);
    const auto *nsh0_699 = buffer.data(nsh0 + 699);
    const auto *nsh0_702 = buffer.data(nsh0 + 702);
    const auto *nsh0_703 = buffer.data(nsh0 + 703);
    const auto *nsh0_705 = buffer.data(nsh0 + 705);
    const auto *nsh0_707 = buffer.data(nsh0 + 707);
    const auto *nsh0_708 = buffer.data(nsh0 + 708);
    const auto *nsh0_710 = buffer.data(nsh0 + 710);
    const auto *nsh0_711 = buffer.data(nsh0 + 711);
    const auto *nsh0_712 = buffer.data(nsh0 + 712);
    const auto *nsh0_713 = buffer.data(nsh0 + 713);
    const auto *nsh0_729 = buffer.data(nsh0 + 729);
    const auto *nsh0_731 = buffer.data(nsh0 + 731);
    const auto *nsh0_732 = buffer.data(nsh0 + 732);
    const auto *nsh0_733 = buffer.data(nsh0 + 733);
    const auto *nsh0_734 = buffer.data(nsh0 + 734);

    const auto *nsh1_675 = buffer.data(nsh1 + 675);
    const auto *nsh1_677 = buffer.data(nsh1 + 677);
    const auto *nsh1_678 = buffer.data(nsh1 + 678);
    const auto *nsh1_681 = buffer.data(nsh1 + 681);
    const auto *nsh1_682 = buffer.data(nsh1 + 682);
    const auto *nsh1_684 = buffer.data(nsh1 + 684);
    const auto *nsh1_686 = buffer.data(nsh1 + 686);
    const auto *nsh1_687 = buffer.data(nsh1 + 687);
    const auto *nsh1_689 = buffer.data(nsh1 + 689);
    const auto *nsh1_690 = buffer.data(nsh1 + 690);
    const auto *nsh1_691 = buffer.data(nsh1 + 691);
    const auto *nsh1_692 = buffer.data(nsh1 + 692);
    const auto *nsh1_693 = buffer.data(nsh1 + 693);
    const auto *nsh1_696 = buffer.data(nsh1 + 696);
    const auto *nsh1_698 = buffer.data(nsh1 + 698);
    const auto *nsh1_699 = buffer.data(nsh1 + 699);
    const auto *nsh1_702 = buffer.data(nsh1 + 702);
    const auto *nsh1_703 = buffer.data(nsh1 + 703);
    const auto *nsh1_705 = buffer.data(nsh1 + 705);
    const auto *nsh1_707 = buffer.data(nsh1 + 707);
    const auto *nsh1_708 = buffer.data(nsh1 + 708);
    const auto *nsh1_710 = buffer.data(nsh1 + 710);
    const auto *nsh1_711 = buffer.data(nsh1 + 711);
    const auto *nsh1_712 = buffer.data(nsh1 + 712);
    const auto *nsh1_713 = buffer.data(nsh1 + 713);
    const auto *nsh1_729 = buffer.data(nsh1 + 729);
    const auto *nsh1_731 = buffer.data(nsh1 + 731);
    const auto *nsh1_732 = buffer.data(nsh1 + 732);
    const auto *nsh1_733 = buffer.data(nsh1 + 733);
    const auto *nsh1_734 = buffer.data(nsh1 + 734);

    const auto *nsi_896 = buffer.data(nsi + 896);
    const auto *nsi_898 = buffer.data(nsi + 898);
    const auto *nsi_899 = buffer.data(nsi + 899);
    const auto *nsi_901 = buffer.data(nsi + 901);
    const auto *nsi_902 = buffer.data(nsi + 902);
    const auto *nsi_905 = buffer.data(nsi + 905);
    const auto *nsi_906 = buffer.data(nsi + 906);
    const auto *nsi_908 = buffer.data(nsi + 908);
    const auto *nsi_910 = buffer.data(nsi + 910);
    const auto *nsi_911 = buffer.data(nsi + 911);
    const auto *nsi_913 = buffer.data(nsi + 913);
    const auto *nsi_914 = buffer.data(nsi + 914);
    const auto *nsi_916 = buffer.data(nsi + 916);
    const auto *nsi_917 = buffer.data(nsi + 917);
    const auto *nsi_918 = buffer.data(nsi + 918);
    const auto *nsi_919 = buffer.data(nsi + 919);
    const auto *nsi_920 = buffer.data(nsi + 920);
    const auto *nsi_921 = buffer.data(nsi + 921);
    const auto *nsi_922 = buffer.data(nsi + 922);
    const auto *nsi_923 = buffer.data(nsi + 923);
    const auto *nsi_924 = buffer.data(nsi + 924);
    const auto *nsi_926 = buffer.data(nsi + 926);
    const auto *nsi_927 = buffer.data(nsi + 927);
    const auto *nsi_929 = buffer.data(nsi + 929);
    const auto *nsi_930 = buffer.data(nsi + 930);
    const auto *nsi_933 = buffer.data(nsi + 933);
    const auto *nsi_934 = buffer.data(nsi + 934);
    const auto *nsi_936 = buffer.data(nsi + 936);
    const auto *nsi_938 = buffer.data(nsi + 938);
    const auto *nsi_939 = buffer.data(nsi + 939);
    const auto *nsi_941 = buffer.data(nsi + 941);
    const auto *nsi_942 = buffer.data(nsi + 942);
    const auto *nsi_944 = buffer.data(nsi + 944);
    const auto *nsi_945 = buffer.data(nsi + 945);
    const auto *nsi_946 = buffer.data(nsi + 946);
    const auto *nsi_947 = buffer.data(nsi + 947);
    const auto *nsi_948 = buffer.data(nsi + 948);
    const auto *nsi_949 = buffer.data(nsi + 949);
    const auto *nsi_950 = buffer.data(nsi + 950);
    const auto *nsi_951 = buffer.data(nsi + 951);
    const auto *nsi_952 = buffer.data(nsi + 952);
    const auto *nsi_954 = buffer.data(nsi + 954);
    const auto *nsi_955 = buffer.data(nsi + 955);
    const auto *nsi_957 = buffer.data(nsi + 957);
    const auto *nsi_958 = buffer.data(nsi + 958);
    const auto *nsi_961 = buffer.data(nsi + 961);
    const auto *nsi_962 = buffer.data(nsi + 962);
    const auto *nsi_966 = buffer.data(nsi + 966);
    const auto *nsi_973 = buffer.data(nsi + 973);
    const auto *nsi_974 = buffer.data(nsi + 974);
    const auto *nsi_975 = buffer.data(nsi + 975);
    const auto *nsi_976 = buffer.data(nsi + 976);
    const auto *nsi_977 = buffer.data(nsi + 977);
    const auto *nsi_978 = buffer.data(nsi + 978);
    const auto *nsi_979 = buffer.data(nsi + 979);

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pc_x, pc_y, pc_z, msi_672, msi_700, \
                         msi_702, msi_899, nsh0_675, nsh1_675, nsi_896, nsi_898, \
                         nsi_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_15 * msi_700[k]
                    + f_3 * pc_y[k] * nsi_896[k];

        t_1154[k] = f_16 * msi_672[k]
                    + f_3 * pc_z[k] * nsi_896[k];

        t_1155[k] = f_15 * msi_899[k]
                    + f_10 * nsh0_675[k]
                    - f_11 * nsh1_675[k]
                    + f_3 * pc_x[k] * nsi_899[k];

        t_1156[k] = f_15 * msi_702[k]
                    + f_3 * pc_y[k] * nsi_898[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, pc_z, msi_675, msi_901, msi_902, \
                         nsh0_677, nsh0_678, nsh1_677, nsh1_678, nsi_899, nsi_901, \
                         nsi_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_15 * msi_901[k]
                    + f_10 * nsh0_677[k]
                    - f_11 * nsh1_677[k]
                    + f_3 * pc_x[k] * nsi_901[k];

        t_1158[k] = f_15 * msi_902[k]
                    + f_8 * nsh0_678[k]
                    - f_9 * nsh1_678[k]
                    + f_3 * pc_x[k] * nsi_902[k];

        t_1159[k] = f_16 * msi_675[k]
                    + f_3 * pc_z[k] * nsi_899[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, msi_705, msi_905, msi_906, \
                         nsh0_681, nsh0_682, nsh1_681, nsh1_682, nsi_901, nsi_905, \
                         nsi_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * msi_705[k]
                    + f_3 * pc_y[k] * nsi_901[k];

        t_1161[k] = f_15 * msi_905[k]
                    + f_8 * nsh0_681[k]
                    - f_9 * nsh1_681[k]
                    + f_3 * pc_x[k] * nsi_905[k];

        t_1162[k] = f_15 * msi_906[k]
                    + f_6 * nsh0_682[k]
                    - f_7 * nsh1_682[k]
                    + f_3 * pc_x[k] * nsi_906[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_x, pc_y, pc_z, msi_678, msi_709, msi_908, \
                         nsh0_684, nsh1_684, nsi_902, nsi_905, \
                         nsi_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * msi_678[k]
                    + f_3 * pc_z[k] * nsi_902[k];

        t_1164[k] = f_15 * msi_908[k]
                    + f_6 * nsh0_684[k]
                    - f_7 * nsh1_684[k]
                    + f_3 * pc_x[k] * nsi_908[k];

        t_1165[k] = f_15 * msi_709[k]
                    + f_3 * pc_y[k] * nsi_905[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_x, pc_z, msi_682, msi_910, msi_911, \
                         nsh0_686, nsh0_687, nsh1_686, nsh1_687, nsi_906, nsi_910, \
                         nsi_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_15 * msi_910[k]
                    + f_6 * nsh0_686[k]
                    - f_7 * nsh1_686[k]
                    + f_3 * pc_x[k] * nsi_910[k];

        t_1167[k] = f_15 * msi_911[k]
                    + f_4 * nsh0_687[k]
                    - f_5 * nsh1_687[k]
                    + f_3 * pc_x[k] * nsi_911[k];

        t_1168[k] = f_16 * msi_682[k]
                    + f_3 * pc_z[k] * nsi_906[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_x, pc_y, msi_714, msi_913, msi_914, \
                         nsh0_689, nsh0_690, nsh1_689, nsh1_690, nsi_910, nsi_913, \
                         nsi_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_15 * msi_913[k]
                    + f_4 * nsh0_689[k]
                    - f_5 * nsh1_689[k]
                    + f_3 * pc_x[k] * nsi_913[k];

        t_1170[k] = f_15 * msi_914[k]
                    + f_4 * nsh0_690[k]
                    - f_5 * nsh1_690[k]
                    + f_3 * pc_x[k] * nsi_914[k];

        t_1171[k] = f_15 * msi_714[k]
                    + f_3 * pc_y[k] * nsi_910[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, pc_x, msi_916, msi_917, msi_918, \
                         msi_919, nsh0_692, nsh1_692, nsi_916, nsi_917, nsi_918, \
                         nsi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_15 * msi_916[k]
                    + f_4 * nsh0_692[k]
                    - f_5 * nsh1_692[k]
                    + f_3 * pc_x[k] * nsi_916[k];

        t_1173[k] = f_15 * msi_917[k]
                    + f_3 * pc_x[k] * nsi_917[k];

        t_1174[k] = f_15 * msi_918[k]
                    + f_3 * pc_x[k] * nsi_918[k];

        t_1175[k] = f_15 * msi_919[k]
                    + f_3 * pc_x[k] * nsi_919[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pc_x, msi_920, msi_921, msi_922, \
                         msi_923, nsi_920, nsi_921, nsi_922, nsi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_15 * msi_920[k]
                    + f_3 * pc_x[k] * nsi_920[k];

        t_1177[k] = f_15 * msi_921[k]
                    + f_3 * pc_x[k] * nsi_921[k];

        t_1178[k] = f_15 * msi_922[k]
                    + f_3 * pc_x[k] * nsi_922[k];

        t_1179[k] = f_15 * msi_923[k]
                    + f_3 * pc_x[k] * nsi_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_y, pc_z, msi_693, msi_721, msi_723, \
                         nsh0_687, nsh0_689, nsh1_687, nsh1_689, nsi_917, \
                         nsi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_15 * msi_721[k]
                    + f_1 * nsh0_687[k]
                    - f_2 * nsh1_687[k]
                    + f_3 * pc_y[k] * nsi_917[k];

        t_1181[k] = f_16 * msi_693[k]
                    + f_3 * pc_z[k] * nsi_917[k];

        t_1182[k] = f_15 * msi_723[k]
                    + f_10 * nsh0_689[k]
                    - f_11 * nsh1_689[k]
                    + f_3 * pc_y[k] * nsi_919[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pc_y, msi_724, msi_725, msi_726, nsh0_690, \
                         nsh0_691, nsh0_692, nsh1_690, nsh1_691, nsh1_692, nsi_920, nsi_921, \
                         nsi_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_15 * msi_724[k]
                    + f_8 * nsh0_690[k]
                    - f_9 * nsh1_690[k]
                    + f_3 * pc_y[k] * nsi_920[k];

        t_1184[k] = f_15 * msi_725[k]
                    + f_6 * nsh0_691[k]
                    - f_7 * nsh1_691[k]
                    + f_3 * pc_y[k] * nsi_921[k];

        t_1185[k] = f_15 * msi_726[k]
                    + f_4 * nsh0_692[k]
                    - f_5 * nsh1_692[k]
                    + f_3 * pc_y[k] * nsi_922[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pc_x, pc_y, pc_z, msi_699, msi_727, msi_924, \
                         nsh0_692, nsh0_693, nsh1_692, nsh1_693, nsi_923, \
                         nsi_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_15 * msi_727[k]
                    + f_3 * pc_y[k] * nsi_923[k];

        t_1187[k] = f_16 * msi_699[k]
                    + f_1 * nsh0_692[k]
                    - f_2 * nsh1_692[k]
                    + f_3 * pc_z[k] * nsi_923[k];

        t_1188[k] = f_15 * msi_924[k]
                    + f_1 * nsh0_693[k]
                    - f_2 * nsh1_693[k]
                    + f_3 * pc_x[k] * nsi_924[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pc_x, pc_y, pc_z, msi_700, msi_728, \
                         msi_730, msi_927, nsh0_696, nsh1_696, nsi_924, nsi_926, \
                         nsi_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_14 * msi_728[k]
                    + f_3 * pc_y[k] * nsi_924[k];

        t_1190[k] = f_17 * msi_700[k]
                    + f_3 * pc_z[k] * nsi_924[k];

        t_1191[k] = f_15 * msi_927[k]
                    + f_10 * nsh0_696[k]
                    - f_11 * nsh1_696[k]
                    + f_3 * pc_x[k] * nsi_927[k];

        t_1192[k] = f_14 * msi_730[k]
                    + f_3 * pc_y[k] * nsi_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pc_x, pc_z, msi_703, msi_929, msi_930, \
                         nsh0_698, nsh0_699, nsh1_698, nsh1_699, nsi_927, nsi_929, \
                         nsi_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_15 * msi_929[k]
                    + f_10 * nsh0_698[k]
                    - f_11 * nsh1_698[k]
                    + f_3 * pc_x[k] * nsi_929[k];

        t_1194[k] = f_15 * msi_930[k]
                    + f_8 * nsh0_699[k]
                    - f_9 * nsh1_699[k]
                    + f_3 * pc_x[k] * nsi_930[k];

        t_1195[k] = f_17 * msi_703[k]
                    + f_3 * pc_z[k] * nsi_927[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, msi_733, msi_933, msi_934, \
                         nsh0_702, nsh0_703, nsh1_702, nsh1_703, nsi_929, nsi_933, \
                         nsi_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_14 * msi_733[k]
                    + f_3 * pc_y[k] * nsi_929[k];

        t_1197[k] = f_15 * msi_933[k]
                    + f_8 * nsh0_702[k]
                    - f_9 * nsh1_702[k]
                    + f_3 * pc_x[k] * nsi_933[k];

        t_1198[k] = f_15 * msi_934[k]
                    + f_6 * nsh0_703[k]
                    - f_7 * nsh1_703[k]
                    + f_3 * pc_x[k] * nsi_934[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_y, pc_z, msi_706, msi_737, msi_936, \
                         nsh0_705, nsh1_705, nsi_930, nsi_933, \
                         nsi_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * msi_706[k]
                    + f_3 * pc_z[k] * nsi_930[k];

        t_1200[k] = f_15 * msi_936[k]
                    + f_6 * nsh0_705[k]
                    - f_7 * nsh1_705[k]
                    + f_3 * pc_x[k] * nsi_936[k];

        t_1201[k] = f_14 * msi_737[k]
                    + f_3 * pc_y[k] * nsi_933[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, pc_z, msi_710, msi_938, msi_939, \
                         nsh0_707, nsh0_708, nsh1_707, nsh1_708, nsi_934, nsi_938, \
                         nsi_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_15 * msi_938[k]
                    + f_6 * nsh0_707[k]
                    - f_7 * nsh1_707[k]
                    + f_3 * pc_x[k] * nsi_938[k];

        t_1203[k] = f_15 * msi_939[k]
                    + f_4 * nsh0_708[k]
                    - f_5 * nsh1_708[k]
                    + f_3 * pc_x[k] * nsi_939[k];

        t_1204[k] = f_17 * msi_710[k]
                    + f_3 * pc_z[k] * nsi_934[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pc_x, pc_y, msi_742, msi_941, msi_942, \
                         nsh0_710, nsh0_711, nsh1_710, nsh1_711, nsi_938, nsi_941, \
                         nsi_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_15 * msi_941[k]
                    + f_4 * nsh0_710[k]
                    - f_5 * nsh1_710[k]
                    + f_3 * pc_x[k] * nsi_941[k];

        t_1206[k] = f_15 * msi_942[k]
                    + f_4 * nsh0_711[k]
                    - f_5 * nsh1_711[k]
                    + f_3 * pc_x[k] * nsi_942[k];

        t_1207[k] = f_14 * msi_742[k]
                    + f_3 * pc_y[k] * nsi_938[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, pc_x, msi_944, msi_945, msi_946, \
                         msi_947, nsh0_713, nsh1_713, nsi_944, nsi_945, nsi_946, \
                         nsi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_15 * msi_944[k]
                    + f_4 * nsh0_713[k]
                    - f_5 * nsh1_713[k]
                    + f_3 * pc_x[k] * nsi_944[k];

        t_1209[k] = f_15 * msi_945[k]
                    + f_3 * pc_x[k] * nsi_945[k];

        t_1210[k] = f_15 * msi_946[k]
                    + f_3 * pc_x[k] * nsi_946[k];

        t_1211[k] = f_15 * msi_947[k]
                    + f_3 * pc_x[k] * nsi_947[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, pc_x, msi_948, msi_949, msi_950, \
                         msi_951, nsi_948, nsi_949, nsi_950, nsi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_15 * msi_948[k]
                    + f_3 * pc_x[k] * nsi_948[k];

        t_1213[k] = f_15 * msi_949[k]
                    + f_3 * pc_x[k] * nsi_949[k];

        t_1214[k] = f_15 * msi_950[k]
                    + f_3 * pc_x[k] * nsi_950[k];

        t_1215[k] = f_15 * msi_951[k]
                    + f_3 * pc_x[k] * nsi_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_y, pc_z, msi_721, msi_749, msi_751, \
                         nsh0_708, nsh0_710, nsh1_708, nsh1_710, nsi_945, \
                         nsi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_14 * msi_749[k]
                    + f_1 * nsh0_708[k]
                    - f_2 * nsh1_708[k]
                    + f_3 * pc_y[k] * nsi_945[k];

        t_1217[k] = f_17 * msi_721[k]
                    + f_3 * pc_z[k] * nsi_945[k];

        t_1218[k] = f_14 * msi_751[k]
                    + f_10 * nsh0_710[k]
                    - f_11 * nsh1_710[k]
                    + f_3 * pc_y[k] * nsi_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pc_y, msi_752, msi_753, msi_754, nsh0_711, \
                         nsh0_712, nsh0_713, nsh1_711, nsh1_712, nsh1_713, nsi_948, nsi_949, \
                         nsi_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_14 * msi_752[k]
                    + f_8 * nsh0_711[k]
                    - f_9 * nsh1_711[k]
                    + f_3 * pc_y[k] * nsi_948[k];

        t_1220[k] = f_14 * msi_753[k]
                    + f_6 * nsh0_712[k]
                    - f_7 * nsh1_712[k]
                    + f_3 * pc_y[k] * nsi_949[k];

        t_1221[k] = f_14 * msi_754[k]
                    + f_4 * nsh0_713[k]
                    - f_5 * nsh1_713[k]
                    + f_3 * pc_y[k] * nsi_950[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pa_y, pc_y, pc_z, msk0_972, msi_727, \
                         msi_755, msi_756, msk1_972, nsh0_713, nsh1_713, nsi_951, \
                         nsi_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * msi_755[k]
                    + f_3 * pc_y[k] * nsi_951[k];

        t_1223[k] = f_17 * msi_727[k]
                    + f_1 * nsh0_713[k]
                    - f_2 * nsh1_713[k]
                    + f_3 * pc_z[k] * nsi_951[k];

        t_1224[k] = pa_y[k] * msk0_972[k]
                    - f_12 * pc_y[k] * msk1_972[k];

        t_1225[k] = f_13 * msi_756[k]
                    + f_3 * pc_y[k] * nsi_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pa_y, pc_y, pc_z, msk0_975, msk0_977, \
                         msi_728, msi_757, msi_758, msk1_975, msk1_977, nsi_952, \
                         nsi_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_23 * msi_728[k]
                    + f_3 * pc_z[k] * nsi_952[k];

        t_1227[k] = pa_y[k] * msk0_975[k]
                    + f_14 * msi_757[k]
                    - f_12 * pc_y[k] * msk1_975[k];

        t_1228[k] = f_13 * msi_758[k]
                    + f_3 * pc_y[k] * nsi_954[k];

        t_1229[k] = pa_y[k] * msk0_977[k]
                    - f_12 * pc_y[k] * msk1_977[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pa_y, pc_y, pc_z, msk0_978, msk0_981, \
                         msi_731, msi_759, msi_761, msk1_978, msk1_981, nsi_955, \
                         nsi_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = pa_y[k] * msk0_978[k]
                    + f_15 * msi_759[k]
                    - f_12 * pc_y[k] * msk1_978[k];

        t_1231[k] = f_23 * msi_731[k]
                    + f_3 * pc_z[k] * nsi_955[k];

        t_1232[k] = f_13 * msi_761[k]
                    + f_3 * pc_y[k] * nsi_957[k];

        t_1233[k] = pa_y[k] * msk0_981[k]
                    - f_12 * pc_y[k] * msk1_981[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pa_y, pc_y, pc_z, msk0_982, msk0_984, \
                         msi_734, msi_762, msi_764, msk1_982, msk1_984, \
                         nsi_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = pa_y[k] * msk0_982[k]
                    + f_16 * msi_762[k]
                    - f_12 * pc_y[k] * msk1_982[k];

        t_1235[k] = f_23 * msi_734[k]
                    + f_3 * pc_z[k] * nsi_958[k];

        t_1236[k] = pa_y[k] * msk0_984[k]
                    + f_14 * msi_764[k]
                    - f_12 * pc_y[k] * msk1_984[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pa_y, pc_y, pc_z, msk0_986, msk0_987, \
                         msi_738, msi_765, msi_766, msk1_986, msk1_987, nsi_961, \
                         nsi_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_13 * msi_765[k]
                    + f_3 * pc_y[k] * nsi_961[k];

        t_1238[k] = pa_y[k] * msk0_986[k]
                    - f_12 * pc_y[k] * msk1_986[k];

        t_1239[k] = pa_y[k] * msk0_987[k]
                    + f_17 * msi_766[k]
                    - f_12 * pc_y[k] * msk1_987[k];

        t_1240[k] = f_23 * msi_738[k]
                    + f_3 * pc_z[k] * nsi_962[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, pa_y, pc_y, msk0_989, msk0_990, \
                         msk0_992, msi_768, msi_769, msi_770, msk1_989, msk1_990, msk1_992, \
                         nsi_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = pa_y[k] * msk0_989[k]
                    + f_15 * msi_768[k]
                    - f_12 * pc_y[k] * msk1_989[k];

        t_1242[k] = pa_y[k] * msk0_990[k]
                    + f_14 * msi_769[k]
                    - f_12 * pc_y[k] * msk1_990[k];

        t_1243[k] = f_13 * msi_770[k]
                    + f_3 * pc_y[k] * nsi_966[k];

        t_1244[k] = pa_y[k] * msk0_992[k]
                    - f_12 * pc_y[k] * msk1_992[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, msi_973, msi_974, \
                         msi_975, msi_976, msi_977, nsi_973, nsi_974, nsi_975, nsi_976, \
                         nsi_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_15 * msi_973[k]
                    + f_3 * pc_x[k] * nsi_973[k];

        t_1246[k] = f_15 * msi_974[k]
                    + f_3 * pc_x[k] * nsi_974[k];

        t_1247[k] = f_15 * msi_975[k]
                    + f_3 * pc_x[k] * nsi_975[k];

        t_1248[k] = f_15 * msi_976[k]
                    + f_3 * pc_x[k] * nsi_976[k];

        t_1249[k] = f_15 * msi_977[k]
                    + f_3 * pc_x[k] * nsi_977[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, pc_z, msi_749, msi_777, \
                         msi_978, msi_979, nsh0_729, nsh1_729, nsi_973, nsi_978, \
                         nsi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_15 * msi_978[k]
                    + f_3 * pc_x[k] * nsi_978[k];

        t_1251[k] = f_15 * msi_979[k]
                    + f_3 * pc_x[k] * nsi_979[k];

        t_1252[k] = f_13 * msi_777[k]
                    + f_1 * nsh0_729[k]
                    - f_2 * nsh1_729[k]
                    + f_3 * pc_y[k] * nsi_973[k];

        t_1253[k] = f_23 * msi_749[k]
                    + f_3 * pc_z[k] * nsi_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, msi_779, msi_780, msi_781, nsh0_731, \
                         nsh0_732, nsh0_733, nsh1_731, nsh1_732, nsh1_733, nsi_975, nsi_976, \
                         nsi_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_13 * msi_779[k]
                    + f_10 * nsh0_731[k]
                    - f_11 * nsh1_731[k]
                    + f_3 * pc_y[k] * nsi_975[k];

        t_1255[k] = f_13 * msi_780[k]
                    + f_8 * nsh0_732[k]
                    - f_9 * nsh1_732[k]
                    + f_3 * pc_y[k] * nsi_976[k];

        t_1256[k] = f_13 * msi_781[k]
                    + f_6 * nsh0_733[k]
                    - f_7 * nsh1_733[k]
                    + f_3 * pc_y[k] * nsi_977[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pa_y, pc_y, msk0_1007, msi_782, msi_783, \
                         msk1_1007, nsh0_734, nsh1_734, nsi_978, \
                         nsi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_13 * msi_782[k]
                    + f_4 * nsh0_734[k]
                    - f_5 * nsh1_734[k]
                    + f_3 * pc_y[k] * nsi_978[k];

        t_1258[k] = f_13 * msi_783[k]
                    + f_3 * pc_y[k] * nsi_979[k];

        t_1259[k] = pa_y[k] * msk0_1007[k]
                    - f_12 * pc_y[k] * msk1_1007[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsh0, const size_t nsh1,
                                                           const size_t nsi, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1008 = buffer.data(msk0 + 1008);
    const auto *msk0_1011 = buffer.data(msk0 + 1011);
    const auto *msk0_1014 = buffer.data(msk0 + 1014);
    const auto *msk0_1018 = buffer.data(msk0 + 1018);
    const auto *msk0_1020 = buffer.data(msk0 + 1020);
    const auto *msk0_1023 = buffer.data(msk0 + 1023);
    const auto *msk0_1025 = buffer.data(msk0 + 1025);
    const auto *msk0_1026 = buffer.data(msk0 + 1026);
    const auto *msk0_1036 = buffer.data(msk0 + 1036);

    const auto *msi_756 = buffer.data(msi + 756);
    const auto *msi_783 = buffer.data(msi + 783);
    const auto *msi_784 = buffer.data(msi + 784);
    const auto *msi_787 = buffer.data(msi + 787);
    const auto *msi_789 = buffer.data(msi + 789);
    const auto *msi_790 = buffer.data(msi + 790);
    const auto *msi_791 = buffer.data(msi + 791);
    const auto *msi_793 = buffer.data(msi + 793);
    const auto *msi_794 = buffer.data(msi + 794);
    const auto *msi_795 = buffer.data(msi + 795);
    const auto *msi_796 = buffer.data(msi + 796);
    const auto *msi_798 = buffer.data(msi + 798);
    const auto *msi_805 = buffer.data(msi + 805);
    const auto *msi_811 = buffer.data(msi + 811);
    const auto *msi_812 = buffer.data(msi + 812);
    const auto *msi_814 = buffer.data(msi + 814);
    const auto *msi_815 = buffer.data(msi + 815);
    const auto *msi_817 = buffer.data(msi + 817);
    const auto *msi_821 = buffer.data(msi + 821);
    const auto *msi_826 = buffer.data(msi + 826);
    const auto *msi_835 = buffer.data(msi + 835);
    const auto *msi_836 = buffer.data(msi + 836);
    const auto *msi_837 = buffer.data(msi + 837);
    const auto *msi_838 = buffer.data(msi + 838);
    const auto *msi_839 = buffer.data(msi + 839);
    const auto *msi_840 = buffer.data(msi + 840);
    const auto *msi_842 = buffer.data(msi + 842);
    const auto *msi_980 = buffer.data(msi + 980);
    const auto *msi_985 = buffer.data(msi + 985);
    const auto *msi_989 = buffer.data(msi + 989);
    const auto *msi_994 = buffer.data(msi + 994);
    const auto *msi_1000 = buffer.data(msi + 1000);
    const auto *msi_1001 = buffer.data(msi + 1001);
    const auto *msi_1002 = buffer.data(msi + 1002);
    const auto *msi_1003 = buffer.data(msi + 1003);
    const auto *msi_1004 = buffer.data(msi + 1004);
    const auto *msi_1005 = buffer.data(msi + 1005);
    const auto *msi_1007 = buffer.data(msi + 1007);
    const auto *msi_1008 = buffer.data(msi + 1008);
    const auto *msi_1011 = buffer.data(msi + 1011);
    const auto *msi_1014 = buffer.data(msi + 1014);
    const auto *msi_1018 = buffer.data(msi + 1018);
    const auto *msi_1023 = buffer.data(msi + 1023);
    const auto *msi_1029 = buffer.data(msi + 1029);
    const auto *msi_1031 = buffer.data(msi + 1031);
    const auto *msi_1032 = buffer.data(msi + 1032);
    const auto *msi_1033 = buffer.data(msi + 1033);
    const auto *msi_1034 = buffer.data(msi + 1034);
    const auto *msi_1035 = buffer.data(msi + 1035);
    const auto *msi_1041 = buffer.data(msi + 1041);
    const auto *msi_1045 = buffer.data(msi + 1045);
    const auto *msi_1050 = buffer.data(msi + 1050);
    const auto *msi_1056 = buffer.data(msi + 1056);
    const auto *msi_1057 = buffer.data(msi + 1057);
    const auto *msi_1058 = buffer.data(msi + 1058);
    const auto *msi_1059 = buffer.data(msi + 1059);
    const auto *msi_1060 = buffer.data(msi + 1060);
    const auto *msi_1061 = buffer.data(msi + 1061);
    const auto *msi_1062 = buffer.data(msi + 1062);
    const auto *msi_1063 = buffer.data(msi + 1063);
    const auto *msi_1064 = buffer.data(msi + 1064);
    const auto *msi_1067 = buffer.data(msi + 1067);
    const auto *msi_1069 = buffer.data(msi + 1069);
    const auto *msi_1070 = buffer.data(msi + 1070);

    const auto *msk1_1008 = buffer.data(msk1 + 1008);
    const auto *msk1_1011 = buffer.data(msk1 + 1011);
    const auto *msk1_1014 = buffer.data(msk1 + 1014);
    const auto *msk1_1018 = buffer.data(msk1 + 1018);
    const auto *msk1_1020 = buffer.data(msk1 + 1020);
    const auto *msk1_1023 = buffer.data(msk1 + 1023);
    const auto *msk1_1025 = buffer.data(msk1 + 1025);
    const auto *msk1_1026 = buffer.data(msk1 + 1026);
    const auto *msk1_1036 = buffer.data(msk1 + 1036);

    const auto *nsh0_735 = buffer.data(nsh0 + 735);
    const auto *nsh0_736 = buffer.data(nsh0 + 736);
    const auto *nsh0_737 = buffer.data(nsh0 + 737);
    const auto *nsh0_738 = buffer.data(nsh0 + 738);
    const auto *nsh0_739 = buffer.data(nsh0 + 739);
    const auto *nsh0_740 = buffer.data(nsh0 + 740);
    const auto *nsh0_741 = buffer.data(nsh0 + 741);
    const auto *nsh0_742 = buffer.data(nsh0 + 742);
    const auto *nsh0_743 = buffer.data(nsh0 + 743);
    const auto *nsh0_744 = buffer.data(nsh0 + 744);
    const auto *nsh0_749 = buffer.data(nsh0 + 749);
    const auto *nsh0_750 = buffer.data(nsh0 + 750);
    const auto *nsh0_751 = buffer.data(nsh0 + 751);
    const auto *nsh0_752 = buffer.data(nsh0 + 752);
    const auto *nsh0_753 = buffer.data(nsh0 + 753);
    const auto *nsh0_754 = buffer.data(nsh0 + 754);
    const auto *nsh0_755 = buffer.data(nsh0 + 755);
    const auto *nsh0_756 = buffer.data(nsh0 + 756);
    const auto *nsh0_758 = buffer.data(nsh0 + 758);
    const auto *nsh0_759 = buffer.data(nsh0 + 759);
    const auto *nsh0_761 = buffer.data(nsh0 + 761);
    const auto *nsh0_762 = buffer.data(nsh0 + 762);
    const auto *nsh0_763 = buffer.data(nsh0 + 763);
    const auto *nsh0_765 = buffer.data(nsh0 + 765);
    const auto *nsh0_766 = buffer.data(nsh0 + 766);
    const auto *nsh0_771 = buffer.data(nsh0 + 771);
    const auto *nsh0_772 = buffer.data(nsh0 + 772);
    const auto *nsh0_773 = buffer.data(nsh0 + 773);
    const auto *nsh0_774 = buffer.data(nsh0 + 774);
    const auto *nsh0_776 = buffer.data(nsh0 + 776);
    const auto *nsh0_782 = buffer.data(nsh0 + 782);
    const auto *nsh0_786 = buffer.data(nsh0 + 786);
    const auto *nsh0_791 = buffer.data(nsh0 + 791);
    const auto *nsh0_794 = buffer.data(nsh0 + 794);
    const auto *nsh0_795 = buffer.data(nsh0 + 795);
    const auto *nsh0_796 = buffer.data(nsh0 + 796);
    const auto *nsh0_797 = buffer.data(nsh0 + 797);
    const auto *nsh0_798 = buffer.data(nsh0 + 798);
    const auto *nsh0_801 = buffer.data(nsh0 + 801);
    const auto *nsh0_803 = buffer.data(nsh0 + 803);
    const auto *nsh0_804 = buffer.data(nsh0 + 804);

    const auto *nsh1_735 = buffer.data(nsh1 + 735);
    const auto *nsh1_736 = buffer.data(nsh1 + 736);
    const auto *nsh1_737 = buffer.data(nsh1 + 737);
    const auto *nsh1_738 = buffer.data(nsh1 + 738);
    const auto *nsh1_739 = buffer.data(nsh1 + 739);
    const auto *nsh1_740 = buffer.data(nsh1 + 740);
    const auto *nsh1_741 = buffer.data(nsh1 + 741);
    const auto *nsh1_742 = buffer.data(nsh1 + 742);
    const auto *nsh1_743 = buffer.data(nsh1 + 743);
    const auto *nsh1_744 = buffer.data(nsh1 + 744);
    const auto *nsh1_749 = buffer.data(nsh1 + 749);
    const auto *nsh1_750 = buffer.data(nsh1 + 750);
    const auto *nsh1_751 = buffer.data(nsh1 + 751);
    const auto *nsh1_752 = buffer.data(nsh1 + 752);
    const auto *nsh1_753 = buffer.data(nsh1 + 753);
    const auto *nsh1_754 = buffer.data(nsh1 + 754);
    const auto *nsh1_755 = buffer.data(nsh1 + 755);
    const auto *nsh1_756 = buffer.data(nsh1 + 756);
    const auto *nsh1_758 = buffer.data(nsh1 + 758);
    const auto *nsh1_759 = buffer.data(nsh1 + 759);
    const auto *nsh1_761 = buffer.data(nsh1 + 761);
    const auto *nsh1_762 = buffer.data(nsh1 + 762);
    const auto *nsh1_763 = buffer.data(nsh1 + 763);
    const auto *nsh1_765 = buffer.data(nsh1 + 765);
    const auto *nsh1_766 = buffer.data(nsh1 + 766);
    const auto *nsh1_771 = buffer.data(nsh1 + 771);
    const auto *nsh1_772 = buffer.data(nsh1 + 772);
    const auto *nsh1_773 = buffer.data(nsh1 + 773);
    const auto *nsh1_774 = buffer.data(nsh1 + 774);
    const auto *nsh1_776 = buffer.data(nsh1 + 776);
    const auto *nsh1_782 = buffer.data(nsh1 + 782);
    const auto *nsh1_786 = buffer.data(nsh1 + 786);
    const auto *nsh1_791 = buffer.data(nsh1 + 791);
    const auto *nsh1_794 = buffer.data(nsh1 + 794);
    const auto *nsh1_795 = buffer.data(nsh1 + 795);
    const auto *nsh1_796 = buffer.data(nsh1 + 796);
    const auto *nsh1_797 = buffer.data(nsh1 + 797);
    const auto *nsh1_798 = buffer.data(nsh1 + 798);
    const auto *nsh1_801 = buffer.data(nsh1 + 801);
    const auto *nsh1_803 = buffer.data(nsh1 + 803);
    const auto *nsh1_804 = buffer.data(nsh1 + 804);

    const auto *nsi_980 = buffer.data(nsi + 980);
    const auto *nsi_981 = buffer.data(nsi + 981);
    const auto *nsi_982 = buffer.data(nsi + 982);
    const auto *nsi_983 = buffer.data(nsi + 983);
    const auto *nsi_984 = buffer.data(nsi + 984);
    const auto *nsi_985 = buffer.data(nsi + 985);
    const auto *nsi_986 = buffer.data(nsi + 986);
    const auto *nsi_987 = buffer.data(nsi + 987);
    const auto *nsi_988 = buffer.data(nsi + 988);
    const auto *nsi_989 = buffer.data(nsi + 989);
    const auto *nsi_990 = buffer.data(nsi + 990);
    const auto *nsi_991 = buffer.data(nsi + 991);
    const auto *nsi_992 = buffer.data(nsi + 992);
    const auto *nsi_993 = buffer.data(nsi + 993);
    const auto *nsi_994 = buffer.data(nsi + 994);
    const auto *nsi_1000 = buffer.data(nsi + 1000);
    const auto *nsi_1001 = buffer.data(nsi + 1001);
    const auto *nsi_1002 = buffer.data(nsi + 1002);
    const auto *nsi_1003 = buffer.data(nsi + 1003);
    const auto *nsi_1004 = buffer.data(nsi + 1004);
    const auto *nsi_1005 = buffer.data(nsi + 1005);
    const auto *nsi_1006 = buffer.data(nsi + 1006);
    const auto *nsi_1007 = buffer.data(nsi + 1007);
    const auto *nsi_1008 = buffer.data(nsi + 1008);
    const auto *nsi_1009 = buffer.data(nsi + 1009);
    const auto *nsi_1010 = buffer.data(nsi + 1010);
    const auto *nsi_1011 = buffer.data(nsi + 1011);
    const auto *nsi_1013 = buffer.data(nsi + 1013);
    const auto *nsi_1014 = buffer.data(nsi + 1014);
    const auto *nsi_1015 = buffer.data(nsi + 1015);
    const auto *nsi_1017 = buffer.data(nsi + 1017);
    const auto *nsi_1018 = buffer.data(nsi + 1018);
    const auto *nsi_1019 = buffer.data(nsi + 1019);
    const auto *nsi_1020 = buffer.data(nsi + 1020);
    const auto *nsi_1022 = buffer.data(nsi + 1022);
    const auto *nsi_1023 = buffer.data(nsi + 1023);
    const auto *nsi_1029 = buffer.data(nsi + 1029);
    const auto *nsi_1030 = buffer.data(nsi + 1030);
    const auto *nsi_1031 = buffer.data(nsi + 1031);
    const auto *nsi_1032 = buffer.data(nsi + 1032);
    const auto *nsi_1033 = buffer.data(nsi + 1033);
    const auto *nsi_1034 = buffer.data(nsi + 1034);
    const auto *nsi_1035 = buffer.data(nsi + 1035);
    const auto *nsi_1036 = buffer.data(nsi + 1036);
    const auto *nsi_1038 = buffer.data(nsi + 1038);
    const auto *nsi_1039 = buffer.data(nsi + 1039);
    const auto *nsi_1041 = buffer.data(nsi + 1041);
    const auto *nsi_1042 = buffer.data(nsi + 1042);
    const auto *nsi_1045 = buffer.data(nsi + 1045);
    const auto *nsi_1046 = buffer.data(nsi + 1046);
    const auto *nsi_1050 = buffer.data(nsi + 1050);
    const auto *nsi_1056 = buffer.data(nsi + 1056);
    const auto *nsi_1057 = buffer.data(nsi + 1057);
    const auto *nsi_1058 = buffer.data(nsi + 1058);
    const auto *nsi_1059 = buffer.data(nsi + 1059);
    const auto *nsi_1060 = buffer.data(nsi + 1060);
    const auto *nsi_1061 = buffer.data(nsi + 1061);
    const auto *nsi_1062 = buffer.data(nsi + 1062);
    const auto *nsi_1063 = buffer.data(nsi + 1063);
    const auto *nsi_1064 = buffer.data(nsi + 1064);
    const auto *nsi_1066 = buffer.data(nsi + 1066);
    const auto *nsi_1067 = buffer.data(nsi + 1067);
    const auto *nsi_1069 = buffer.data(nsi + 1069);
    const auto *nsi_1070 = buffer.data(nsi + 1070);

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, t_1264, pc_x, pc_y, pc_z, msi_756, \
                         msi_980, nsh0_735, nsh1_735, nsi_980, nsi_981, \
                         nsi_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_15 * msi_980[k]
                    + f_1 * nsh0_735[k]
                    - f_2 * nsh1_735[k]
                    + f_3 * pc_x[k] * nsi_980[k];

        t_1261[k] = f_3 * pc_y[k] * nsi_980[k];

        t_1262[k] = f_22 * msi_756[k]
                    + f_3 * pc_z[k] * nsi_980[k];

        t_1263[k] = f_4 * nsh0_735[k]
                    - f_5 * nsh1_735[k]
                    + f_3 * pc_y[k] * nsi_981[k];

        t_1264[k] = f_3 * pc_y[k] * nsi_982[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, t_1268, pc_x, pc_y, msi_985, nsh0_736, \
                         nsh0_737, nsh0_740, nsh1_736, nsh1_737, nsh1_740, nsi_983, nsi_984, \
                         nsi_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_15 * msi_985[k]
                    + f_10 * nsh0_740[k]
                    - f_11 * nsh1_740[k]
                    + f_3 * pc_x[k] * nsi_985[k];

        t_1266[k] = f_6 * nsh0_736[k]
                    - f_7 * nsh1_736[k]
                    + f_3 * pc_y[k] * nsi_983[k];

        t_1267[k] = f_4 * nsh0_737[k]
                    - f_5 * nsh1_737[k]
                    + f_3 * pc_y[k] * nsi_984[k];

        t_1268[k] = f_3 * pc_y[k] * nsi_985[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, pc_x, pc_y, msi_989, nsh0_738, nsh0_739, \
                         nsh0_744, nsh1_738, nsh1_739, nsh1_744, nsi_986, nsi_987, \
                         nsi_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = f_15 * msi_989[k]
                    + f_8 * nsh0_744[k]
                    - f_9 * nsh1_744[k]
                    + f_3 * pc_x[k] * nsi_989[k];

        t_1270[k] = f_8 * nsh0_738[k]
                    - f_9 * nsh1_738[k]
                    + f_3 * pc_y[k] * nsi_986[k];

        t_1271[k] = f_6 * nsh0_739[k]
                    - f_7 * nsh1_739[k]
                    + f_3 * pc_y[k] * nsi_987[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pc_x, pc_y, msi_994, nsh0_740, nsh0_749, \
                         nsh1_740, nsh1_749, nsi_988, nsi_989, \
                         nsi_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * nsh0_740[k]
                    - f_5 * nsh1_740[k]
                    + f_3 * pc_y[k] * nsi_988[k];

        t_1273[k] = f_3 * pc_y[k] * nsi_989[k];

        t_1274[k] = f_15 * msi_994[k]
                    + f_6 * nsh0_749[k]
                    - f_7 * nsh1_749[k]
                    + f_3 * pc_x[k] * nsi_994[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pc_y, nsh0_741, nsh0_742, nsh0_743, nsh1_741, \
                         nsh1_742, nsh1_743, nsi_990, nsi_991, \
                         nsi_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_10 * nsh0_741[k]
                    - f_11 * nsh1_741[k]
                    + f_3 * pc_y[k] * nsi_990[k];

        t_1276[k] = f_8 * nsh0_742[k]
                    - f_9 * nsh1_742[k]
                    + f_3 * pc_y[k] * nsi_991[k];

        t_1277[k] = f_6 * nsh0_743[k]
                    - f_7 * nsh1_743[k]
                    + f_3 * pc_y[k] * nsi_992[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, t_1281, pc_x, pc_y, msi_1000, msi_1001, \
                         nsh0_744, nsh0_755, nsh1_744, nsh1_755, nsi_993, nsi_994, nsi_1000, \
                         nsi_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_4 * nsh0_744[k]
                    - f_5 * nsh1_744[k]
                    + f_3 * pc_y[k] * nsi_993[k];

        t_1279[k] = f_3 * pc_y[k] * nsi_994[k];

        t_1280[k] = f_15 * msi_1000[k]
                    + f_4 * nsh0_755[k]
                    - f_5 * nsh1_755[k]
                    + f_3 * pc_x[k] * nsi_1000[k];

        t_1281[k] = f_15 * msi_1001[k]
                    + f_3 * pc_x[k] * nsi_1001[k];
    }

#pragma omp simd aligned(t_1282, t_1283, t_1284, t_1285, t_1286, pc_x, pc_y, msi_1002, \
                         msi_1003, msi_1004, msi_1005, nsi_1000, nsi_1002, nsi_1003, nsi_1004, \
                         nsi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1282[k] = f_15 * msi_1002[k]
                    + f_3 * pc_x[k] * nsi_1002[k];

        t_1283[k] = f_15 * msi_1003[k]
                    + f_3 * pc_x[k] * nsi_1003[k];

        t_1284[k] = f_15 * msi_1004[k]
                    + f_3 * pc_x[k] * nsi_1004[k];

        t_1285[k] = f_15 * msi_1005[k]
                    + f_3 * pc_x[k] * nsi_1005[k];

        t_1286[k] = f_3 * pc_y[k] * nsi_1000[k];
    }

#pragma omp simd aligned(t_1287, t_1288, t_1289, pc_x, pc_y, msi_1007, nsh0_750, nsh0_751, \
                         nsh1_750, nsh1_751, nsi_1001, nsi_1002, \
                         nsi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1287[k] = f_15 * msi_1007[k]
                    + f_3 * pc_x[k] * nsi_1007[k];

        t_1288[k] = f_1 * nsh0_750[k]
                    - f_2 * nsh1_750[k]
                    + f_3 * pc_y[k] * nsi_1001[k];

        t_1289[k] = f_19 * nsh0_751[k]
                    - f_20 * nsh1_751[k]
                    + f_3 * pc_y[k] * nsi_1002[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, pc_y, nsh0_752, nsh0_753, nsh0_754, nsh1_752, \
                         nsh1_753, nsh1_754, nsi_1003, nsi_1004, \
                         nsi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_10 * nsh0_752[k]
                    - f_11 * nsh1_752[k]
                    + f_3 * pc_y[k] * nsi_1003[k];

        t_1291[k] = f_8 * nsh0_753[k]
                    - f_9 * nsh1_753[k]
                    + f_3 * pc_y[k] * nsi_1004[k];

        t_1292[k] = f_6 * nsh0_754[k]
                    - f_7 * nsh1_754[k]
                    + f_3 * pc_y[k] * nsi_1005[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pc_x, pc_y, pc_z, msi_783, msi_1008, \
                         nsh0_755, nsh0_756, nsh1_755, nsh1_756, nsi_1006, nsi_1007, \
                         nsi_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_4 * nsh0_755[k]
                    - f_5 * nsh1_755[k]
                    + f_3 * pc_y[k] * nsi_1006[k];

        t_1294[k] = f_3 * pc_y[k] * nsi_1007[k];

        t_1295[k] = f_22 * msi_783[k]
                    + f_1 * nsh0_755[k]
                    - f_2 * nsh1_755[k]
                    + f_3 * pc_z[k] * nsi_1007[k];

        t_1296[k] = f_14 * msi_1008[k]
                    + f_1 * nsh0_756[k]
                    - f_2 * nsh1_756[k]
                    + f_3 * pc_x[k] * nsi_1008[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, t_1300, pc_x, pc_y, pc_z, msi_784, msi_1011, \
                         nsh0_759, nsh1_759, nsi_1008, nsi_1009, \
                         nsi_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = f_21 * msi_784[k]
                    + f_3 * pc_y[k] * nsi_1008[k];

        t_1298[k] = f_3 * pc_z[k] * nsi_1008[k];

        t_1299[k] = f_14 * msi_1011[k]
                    + f_10 * nsh0_759[k]
                    - f_11 * nsh1_759[k]
                    + f_3 * pc_x[k] * nsi_1011[k];

        t_1300[k] = f_3 * pc_z[k] * nsi_1009[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, pc_x, pc_z, msi_1014, nsh0_756, nsh0_762, \
                         nsh1_756, nsh1_762, nsi_1010, nsi_1011, \
                         nsi_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_4 * nsh0_756[k]
                    - f_5 * nsh1_756[k]
                    + f_3 * pc_z[k] * nsi_1010[k];

        t_1302[k] = f_14 * msi_1014[k]
                    + f_8 * nsh0_762[k]
                    - f_9 * nsh1_762[k]
                    + f_3 * pc_x[k] * nsi_1014[k];

        t_1303[k] = f_3 * pc_z[k] * nsi_1011[k];
    }

#pragma omp simd aligned(t_1304, t_1305, t_1306, t_1307, pc_x, pc_y, pc_z, msi_789, msi_1018, \
                         nsh0_758, nsh0_766, nsh1_758, nsh1_766, nsi_1013, nsi_1014, \
                         nsi_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1304[k] = f_21 * msi_789[k]
                    + f_3 * pc_y[k] * nsi_1013[k];

        t_1305[k] = f_6 * nsh0_758[k]
                    - f_7 * nsh1_758[k]
                    + f_3 * pc_z[k] * nsi_1013[k];

        t_1306[k] = f_14 * msi_1018[k]
                    + f_6 * nsh0_766[k]
                    - f_7 * nsh1_766[k]
                    + f_3 * pc_x[k] * nsi_1018[k];

        t_1307[k] = f_3 * pc_z[k] * nsi_1014[k];
    }

#pragma omp simd aligned(t_1308, t_1309, t_1310, pc_y, pc_z, msi_793, nsh0_759, nsh0_761, \
                         nsh1_759, nsh1_761, nsi_1015, nsi_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1308[k] = f_4 * nsh0_759[k]
                    - f_5 * nsh1_759[k]
                    + f_3 * pc_z[k] * nsi_1015[k];

        t_1309[k] = f_21 * msi_793[k]
                    + f_3 * pc_y[k] * nsi_1017[k];

        t_1310[k] = f_8 * nsh0_761[k]
                    - f_9 * nsh1_761[k]
                    + f_3 * pc_z[k] * nsi_1017[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, pc_x, pc_z, msi_1023, nsh0_762, nsh0_771, \
                         nsh1_762, nsh1_771, nsi_1018, nsi_1019, \
                         nsi_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = f_14 * msi_1023[k]
                    + f_4 * nsh0_771[k]
                    - f_5 * nsh1_771[k]
                    + f_3 * pc_x[k] * nsi_1023[k];

        t_1312[k] = f_3 * pc_z[k] * nsi_1018[k];

        t_1313[k] = f_4 * nsh0_762[k]
                    - f_5 * nsh1_762[k]
                    + f_3 * pc_z[k] * nsi_1019[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, t_1317, pc_x, pc_y, pc_z, msi_798, msi_1029, \
                         nsh0_763, nsh0_765, nsh1_763, nsh1_765, nsi_1020, nsi_1022, \
                         nsi_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_6 * nsh0_763[k]
                    - f_7 * nsh1_763[k]
                    + f_3 * pc_z[k] * nsi_1020[k];

        t_1315[k] = f_21 * msi_798[k]
                    + f_3 * pc_y[k] * nsi_1022[k];

        t_1316[k] = f_10 * nsh0_765[k]
                    - f_11 * nsh1_765[k]
                    + f_3 * pc_z[k] * nsi_1022[k];

        t_1317[k] = f_14 * msi_1029[k]
                    + f_3 * pc_x[k] * nsi_1029[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, t_1321, t_1322, pc_x, pc_z, msi_1031, \
                         msi_1032, msi_1033, msi_1034, nsi_1023, nsi_1031, nsi_1032, nsi_1033, \
                         nsi_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_3 * pc_z[k] * nsi_1023[k];

        t_1319[k] = f_14 * msi_1031[k]
                    + f_3 * pc_x[k] * nsi_1031[k];

        t_1320[k] = f_14 * msi_1032[k]
                    + f_3 * pc_x[k] * nsi_1032[k];

        t_1321[k] = f_14 * msi_1033[k]
                    + f_3 * pc_x[k] * nsi_1033[k];

        t_1322[k] = f_14 * msi_1034[k]
                    + f_3 * pc_x[k] * nsi_1034[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, t_1326, pc_x, pc_y, pc_z, msi_805, msi_1035, \
                         nsh0_771, nsh1_771, nsi_1029, nsi_1030, \
                         nsi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = f_14 * msi_1035[k]
                    + f_3 * pc_x[k] * nsi_1035[k];

        t_1324[k] = f_21 * msi_805[k]
                    + f_1 * nsh0_771[k]
                    - f_2 * nsh1_771[k]
                    + f_3 * pc_y[k] * nsi_1029[k];

        t_1325[k] = f_3 * pc_z[k] * nsi_1029[k];

        t_1326[k] = f_4 * nsh0_771[k]
                    - f_5 * nsh1_771[k]
                    + f_3 * pc_z[k] * nsi_1030[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, pc_z, nsh0_772, nsh0_773, nsh0_774, nsh1_772, \
                         nsh1_773, nsh1_774, nsi_1031, nsi_1032, \
                         nsi_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_6 * nsh0_772[k]
                    - f_7 * nsh1_772[k]
                    + f_3 * pc_z[k] * nsi_1031[k];

        t_1328[k] = f_8 * nsh0_773[k]
                    - f_9 * nsh1_773[k]
                    + f_3 * pc_z[k] * nsi_1032[k];

        t_1329[k] = f_10 * nsh0_774[k]
                    - f_11 * nsh1_774[k]
                    + f_3 * pc_z[k] * nsi_1033[k];
    }

#pragma omp simd aligned(t_1330, t_1331, t_1332, t_1333, pa_z, pc_y, pc_z, msk0_1008, msi_811, \
                         msi_812, msk1_1008, nsh0_776, nsh1_776, nsi_1035, \
                         nsi_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1330[k] = f_21 * msi_811[k]
                    + f_3 * pc_y[k] * nsi_1035[k];

        t_1331[k] = f_1 * nsh0_776[k]
                    - f_2 * nsh1_776[k]
                    + f_3 * pc_z[k] * nsi_1035[k];

        t_1332[k] = pa_z[k] * msk0_1008[k]
                    - f_12 * pc_z[k] * msk1_1008[k];

        t_1333[k] = f_22 * msi_812[k]
                    + f_3 * pc_y[k] * nsi_1036[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, pa_z, pc_y, pc_z, msk0_1011, msi_784, \
                         msi_814, msk1_1011, nsi_1036, nsi_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_13 * msi_784[k]
                    + f_3 * pc_z[k] * nsi_1036[k];

        t_1335[k] = pa_z[k] * msk0_1011[k]
                    - f_12 * pc_z[k] * msk1_1011[k];

        t_1336[k] = f_22 * msi_814[k]
                    + f_3 * pc_y[k] * nsi_1038[k];
    }

#pragma omp simd aligned(t_1337, t_1338, t_1339, pa_z, pc_x, pc_z, msk0_1014, msi_787, \
                         msi_1041, msk1_1014, nsh0_782, nsh1_782, nsi_1039, \
                         nsi_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1337[k] = f_14 * msi_1041[k]
                    + f_10 * nsh0_782[k]
                    - f_11 * nsh1_782[k]
                    + f_3 * pc_x[k] * nsi_1041[k];

        t_1338[k] = pa_z[k] * msk0_1014[k]
                    - f_12 * pc_z[k] * msk1_1014[k];

        t_1339[k] = f_13 * msi_787[k]
                    + f_3 * pc_z[k] * nsi_1039[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pa_z, pc_x, pc_y, pc_z, msk0_1018, msi_817, \
                         msi_1045, msk1_1018, nsh0_786, nsh1_786, nsi_1041, \
                         nsi_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_22 * msi_817[k]
                    + f_3 * pc_y[k] * nsi_1041[k];

        t_1341[k] = f_14 * msi_1045[k]
                    + f_8 * nsh0_786[k]
                    - f_9 * nsh1_786[k]
                    + f_3 * pc_x[k] * nsi_1045[k];

        t_1342[k] = pa_z[k] * msk0_1018[k]
                    - f_12 * pc_z[k] * msk1_1018[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pa_z, pc_y, pc_z, msk0_1020, msi_790, \
                         msi_791, msi_821, msk1_1020, nsi_1042, \
                         nsi_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_13 * msi_790[k]
                    + f_3 * pc_z[k] * nsi_1042[k];

        t_1344[k] = pa_z[k] * msk0_1020[k]
                    + f_14 * msi_791[k]
                    - f_12 * pc_z[k] * msk1_1020[k];

        t_1345[k] = f_22 * msi_821[k]
                    + f_3 * pc_y[k] * nsi_1045[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pa_z, pc_x, pc_z, msk0_1023, msi_794, \
                         msi_1050, msk1_1023, nsh0_791, nsh1_791, nsi_1046, \
                         nsi_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_14 * msi_1050[k]
                    + f_6 * nsh0_791[k]
                    - f_7 * nsh1_791[k]
                    + f_3 * pc_x[k] * nsi_1050[k];

        t_1347[k] = pa_z[k] * msk0_1023[k]
                    - f_12 * pc_z[k] * msk1_1023[k];

        t_1348[k] = f_13 * msi_794[k]
                    + f_3 * pc_z[k] * nsi_1046[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pa_z, pc_y, pc_z, msk0_1025, msk0_1026, \
                         msi_795, msi_796, msi_826, msk1_1025, msk1_1026, \
                         nsi_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = pa_z[k] * msk0_1025[k]
                    + f_14 * msi_795[k]
                    - f_12 * pc_z[k] * msk1_1025[k];

        t_1350[k] = pa_z[k] * msk0_1026[k]
                    + f_15 * msi_796[k]
                    - f_12 * pc_z[k] * msk1_1026[k];

        t_1351[k] = f_22 * msi_826[k]
                    + f_3 * pc_y[k] * nsi_1050[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, t_1355, pc_x, msi_1056, msi_1057, msi_1058, \
                         msi_1059, nsh0_797, nsh1_797, nsi_1056, nsi_1057, nsi_1058, \
                         nsi_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_14 * msi_1056[k]
                    + f_4 * nsh0_797[k]
                    - f_5 * nsh1_797[k]
                    + f_3 * pc_x[k] * nsi_1056[k];

        t_1353[k] = f_14 * msi_1057[k]
                    + f_3 * pc_x[k] * nsi_1057[k];

        t_1354[k] = f_14 * msi_1058[k]
                    + f_3 * pc_x[k] * nsi_1058[k];

        t_1355[k] = f_14 * msi_1059[k]
                    + f_3 * pc_x[k] * nsi_1059[k];
    }

#pragma omp simd aligned(t_1356, t_1357, t_1358, t_1359, pc_x, msi_1060, msi_1061, msi_1062, \
                         msi_1063, nsi_1060, nsi_1061, nsi_1062, \
                         nsi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1356[k] = f_14 * msi_1060[k]
                    + f_3 * pc_x[k] * nsi_1060[k];

        t_1357[k] = f_14 * msi_1061[k]
                    + f_3 * pc_x[k] * nsi_1061[k];

        t_1358[k] = f_14 * msi_1062[k]
                    + f_3 * pc_x[k] * nsi_1062[k];

        t_1359[k] = f_14 * msi_1063[k]
                    + f_3 * pc_x[k] * nsi_1063[k];
    }

#pragma omp simd aligned(t_1360, t_1361, t_1362, pa_z, pc_y, pc_z, msk0_1036, msi_805, \
                         msi_835, msk1_1036, nsh0_794, nsh1_794, nsi_1057, \
                         nsi_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1360[k] = pa_z[k] * msk0_1036[k]
                    - f_12 * pc_z[k] * msk1_1036[k];

        t_1361[k] = f_13 * msi_805[k]
                    + f_3 * pc_z[k] * nsi_1057[k];

        t_1362[k] = f_22 * msi_835[k]
                    + f_10 * nsh0_794[k]
                    - f_11 * nsh1_794[k]
                    + f_3 * pc_y[k] * nsi_1059[k];
    }

#pragma omp simd aligned(t_1363, t_1364, t_1365, pc_y, msi_836, msi_837, msi_838, nsh0_795, \
                         nsh0_796, nsh0_797, nsh1_795, nsh1_796, nsh1_797, nsi_1060, nsi_1061, \
                         nsi_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1363[k] = f_22 * msi_836[k]
                    + f_8 * nsh0_795[k]
                    - f_9 * nsh1_795[k]
                    + f_3 * pc_y[k] * nsi_1060[k];

        t_1364[k] = f_22 * msi_837[k]
                    + f_6 * nsh0_796[k]
                    - f_7 * nsh1_796[k]
                    + f_3 * pc_y[k] * nsi_1061[k];

        t_1365[k] = f_22 * msi_838[k]
                    + f_4 * nsh0_797[k]
                    - f_5 * nsh1_797[k]
                    + f_3 * pc_y[k] * nsi_1062[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, pc_x, pc_y, pc_z, msi_811, msi_839, msi_1064, \
                         nsh0_797, nsh0_798, nsh1_797, nsh1_798, nsi_1063, \
                         nsi_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_22 * msi_839[k]
                    + f_3 * pc_y[k] * nsi_1063[k];

        t_1367[k] = f_13 * msi_811[k]
                    + f_1 * nsh0_797[k]
                    - f_2 * nsh1_797[k]
                    + f_3 * pc_z[k] * nsi_1063[k];

        t_1368[k] = f_14 * msi_1064[k]
                    + f_1 * nsh0_798[k]
                    - f_2 * nsh1_798[k]
                    + f_3 * pc_x[k] * nsi_1064[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, t_1372, pc_x, pc_y, pc_z, msi_812, msi_840, \
                         msi_842, msi_1067, nsh0_801, nsh1_801, nsi_1064, nsi_1066, \
                         nsi_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_23 * msi_840[k]
                    + f_3 * pc_y[k] * nsi_1064[k];

        t_1370[k] = f_14 * msi_812[k]
                    + f_3 * pc_z[k] * nsi_1064[k];

        t_1371[k] = f_14 * msi_1067[k]
                    + f_10 * nsh0_801[k]
                    - f_11 * nsh1_801[k]
                    + f_3 * pc_x[k] * nsi_1067[k];

        t_1372[k] = f_23 * msi_842[k]
                    + f_3 * pc_y[k] * nsi_1066[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pc_x, pc_z, msi_815, msi_1069, msi_1070, \
                         nsh0_803, nsh0_804, nsh1_803, nsh1_804, nsi_1067, nsi_1069, \
                         nsi_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = f_14 * msi_1069[k]
                    + f_10 * nsh0_803[k]
                    - f_11 * nsh1_803[k]
                    + f_3 * pc_x[k] * nsi_1069[k];

        t_1374[k] = f_14 * msi_1070[k]
                    + f_8 * nsh0_804[k]
                    - f_9 * nsh1_804[k]
                    + f_3 * pc_x[k] * nsi_1070[k];

        t_1375[k] = f_14 * msi_815[k]
                    + f_3 * pc_z[k] * nsi_1067[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t msi, const size_t nsh0,
                                                           const size_t nsh1, const size_t nsi,
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
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi_818 = buffer.data(msi + 818);
    const auto *msi_822 = buffer.data(msi + 822);
    const auto *msi_833 = buffer.data(msi + 833);
    const auto *msi_839 = buffer.data(msi + 839);
    const auto *msi_840 = buffer.data(msi + 840);
    const auto *msi_843 = buffer.data(msi + 843);
    const auto *msi_845 = buffer.data(msi + 845);
    const auto *msi_846 = buffer.data(msi + 846);
    const auto *msi_849 = buffer.data(msi + 849);
    const auto *msi_850 = buffer.data(msi + 850);
    const auto *msi_854 = buffer.data(msi + 854);
    const auto *msi_861 = buffer.data(msi + 861);
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
    const auto *msi_882 = buffer.data(msi + 882);
    const auto *msi_889 = buffer.data(msi + 889);
    const auto *msi_891 = buffer.data(msi + 891);
    const auto *msi_892 = buffer.data(msi + 892);
    const auto *msi_893 = buffer.data(msi + 893);
    const auto *msi_894 = buffer.data(msi + 894);
    const auto *msi_895 = buffer.data(msi + 895);
    const auto *msi_896 = buffer.data(msi + 896);
    const auto *msi_898 = buffer.data(msi + 898);
    const auto *msi_901 = buffer.data(msi + 901);
    const auto *msi_905 = buffer.data(msi + 905);
    const auto *msi_910 = buffer.data(msi + 910);
    const auto *msi_917 = buffer.data(msi + 917);
    const auto *msi_919 = buffer.data(msi + 919);
    const auto *msi_920 = buffer.data(msi + 920);
    const auto *msi_921 = buffer.data(msi + 921);
    const auto *msi_922 = buffer.data(msi + 922);
    const auto *msi_923 = buffer.data(msi + 923);
    const auto *msi_924 = buffer.data(msi + 924);
    const auto *msi_926 = buffer.data(msi + 926);
    const auto *msi_1073 = buffer.data(msi + 1073);
    const auto *msi_1074 = buffer.data(msi + 1074);
    const auto *msi_1076 = buffer.data(msi + 1076);
    const auto *msi_1078 = buffer.data(msi + 1078);
    const auto *msi_1079 = buffer.data(msi + 1079);
    const auto *msi_1081 = buffer.data(msi + 1081);
    const auto *msi_1082 = buffer.data(msi + 1082);
    const auto *msi_1084 = buffer.data(msi + 1084);
    const auto *msi_1085 = buffer.data(msi + 1085);
    const auto *msi_1086 = buffer.data(msi + 1086);
    const auto *msi_1087 = buffer.data(msi + 1087);
    const auto *msi_1088 = buffer.data(msi + 1088);
    const auto *msi_1089 = buffer.data(msi + 1089);
    const auto *msi_1090 = buffer.data(msi + 1090);
    const auto *msi_1091 = buffer.data(msi + 1091);
    const auto *msi_1092 = buffer.data(msi + 1092);
    const auto *msi_1095 = buffer.data(msi + 1095);
    const auto *msi_1097 = buffer.data(msi + 1097);
    const auto *msi_1098 = buffer.data(msi + 1098);
    const auto *msi_1101 = buffer.data(msi + 1101);
    const auto *msi_1102 = buffer.data(msi + 1102);
    const auto *msi_1104 = buffer.data(msi + 1104);
    const auto *msi_1106 = buffer.data(msi + 1106);
    const auto *msi_1107 = buffer.data(msi + 1107);
    const auto *msi_1109 = buffer.data(msi + 1109);
    const auto *msi_1110 = buffer.data(msi + 1110);
    const auto *msi_1112 = buffer.data(msi + 1112);
    const auto *msi_1113 = buffer.data(msi + 1113);
    const auto *msi_1114 = buffer.data(msi + 1114);
    const auto *msi_1115 = buffer.data(msi + 1115);
    const auto *msi_1116 = buffer.data(msi + 1116);
    const auto *msi_1117 = buffer.data(msi + 1117);
    const auto *msi_1118 = buffer.data(msi + 1118);
    const auto *msi_1119 = buffer.data(msi + 1119);
    const auto *msi_1120 = buffer.data(msi + 1120);
    const auto *msi_1123 = buffer.data(msi + 1123);
    const auto *msi_1125 = buffer.data(msi + 1125);
    const auto *msi_1126 = buffer.data(msi + 1126);
    const auto *msi_1129 = buffer.data(msi + 1129);
    const auto *msi_1130 = buffer.data(msi + 1130);
    const auto *msi_1132 = buffer.data(msi + 1132);
    const auto *msi_1134 = buffer.data(msi + 1134);
    const auto *msi_1135 = buffer.data(msi + 1135);
    const auto *msi_1137 = buffer.data(msi + 1137);
    const auto *msi_1138 = buffer.data(msi + 1138);
    const auto *msi_1140 = buffer.data(msi + 1140);
    const auto *msi_1141 = buffer.data(msi + 1141);
    const auto *msi_1142 = buffer.data(msi + 1142);
    const auto *msi_1143 = buffer.data(msi + 1143);
    const auto *msi_1144 = buffer.data(msi + 1144);
    const auto *msi_1145 = buffer.data(msi + 1145);
    const auto *msi_1146 = buffer.data(msi + 1146);
    const auto *msi_1147 = buffer.data(msi + 1147);
    const auto *msi_1148 = buffer.data(msi + 1148);
    const auto *msi_1151 = buffer.data(msi + 1151);

    const auto *nsh0_807 = buffer.data(nsh0 + 807);
    const auto *nsh0_808 = buffer.data(nsh0 + 808);
    const auto *nsh0_810 = buffer.data(nsh0 + 810);
    const auto *nsh0_812 = buffer.data(nsh0 + 812);
    const auto *nsh0_813 = buffer.data(nsh0 + 813);
    const auto *nsh0_815 = buffer.data(nsh0 + 815);
    const auto *nsh0_816 = buffer.data(nsh0 + 816);
    const auto *nsh0_817 = buffer.data(nsh0 + 817);
    const auto *nsh0_818 = buffer.data(nsh0 + 818);
    const auto *nsh0_819 = buffer.data(nsh0 + 819);
    const auto *nsh0_822 = buffer.data(nsh0 + 822);
    const auto *nsh0_824 = buffer.data(nsh0 + 824);
    const auto *nsh0_825 = buffer.data(nsh0 + 825);
    const auto *nsh0_828 = buffer.data(nsh0 + 828);
    const auto *nsh0_829 = buffer.data(nsh0 + 829);
    const auto *nsh0_831 = buffer.data(nsh0 + 831);
    const auto *nsh0_833 = buffer.data(nsh0 + 833);
    const auto *nsh0_834 = buffer.data(nsh0 + 834);
    const auto *nsh0_836 = buffer.data(nsh0 + 836);
    const auto *nsh0_837 = buffer.data(nsh0 + 837);
    const auto *nsh0_838 = buffer.data(nsh0 + 838);
    const auto *nsh0_839 = buffer.data(nsh0 + 839);
    const auto *nsh0_840 = buffer.data(nsh0 + 840);
    const auto *nsh0_843 = buffer.data(nsh0 + 843);
    const auto *nsh0_845 = buffer.data(nsh0 + 845);
    const auto *nsh0_846 = buffer.data(nsh0 + 846);
    const auto *nsh0_849 = buffer.data(nsh0 + 849);
    const auto *nsh0_850 = buffer.data(nsh0 + 850);
    const auto *nsh0_852 = buffer.data(nsh0 + 852);
    const auto *nsh0_854 = buffer.data(nsh0 + 854);
    const auto *nsh0_855 = buffer.data(nsh0 + 855);
    const auto *nsh0_857 = buffer.data(nsh0 + 857);
    const auto *nsh0_858 = buffer.data(nsh0 + 858);
    const auto *nsh0_859 = buffer.data(nsh0 + 859);
    const auto *nsh0_860 = buffer.data(nsh0 + 860);
    const auto *nsh0_861 = buffer.data(nsh0 + 861);
    const auto *nsh0_864 = buffer.data(nsh0 + 864);

    const auto *nsh1_807 = buffer.data(nsh1 + 807);
    const auto *nsh1_808 = buffer.data(nsh1 + 808);
    const auto *nsh1_810 = buffer.data(nsh1 + 810);
    const auto *nsh1_812 = buffer.data(nsh1 + 812);
    const auto *nsh1_813 = buffer.data(nsh1 + 813);
    const auto *nsh1_815 = buffer.data(nsh1 + 815);
    const auto *nsh1_816 = buffer.data(nsh1 + 816);
    const auto *nsh1_817 = buffer.data(nsh1 + 817);
    const auto *nsh1_818 = buffer.data(nsh1 + 818);
    const auto *nsh1_819 = buffer.data(nsh1 + 819);
    const auto *nsh1_822 = buffer.data(nsh1 + 822);
    const auto *nsh1_824 = buffer.data(nsh1 + 824);
    const auto *nsh1_825 = buffer.data(nsh1 + 825);
    const auto *nsh1_828 = buffer.data(nsh1 + 828);
    const auto *nsh1_829 = buffer.data(nsh1 + 829);
    const auto *nsh1_831 = buffer.data(nsh1 + 831);
    const auto *nsh1_833 = buffer.data(nsh1 + 833);
    const auto *nsh1_834 = buffer.data(nsh1 + 834);
    const auto *nsh1_836 = buffer.data(nsh1 + 836);
    const auto *nsh1_837 = buffer.data(nsh1 + 837);
    const auto *nsh1_838 = buffer.data(nsh1 + 838);
    const auto *nsh1_839 = buffer.data(nsh1 + 839);
    const auto *nsh1_840 = buffer.data(nsh1 + 840);
    const auto *nsh1_843 = buffer.data(nsh1 + 843);
    const auto *nsh1_845 = buffer.data(nsh1 + 845);
    const auto *nsh1_846 = buffer.data(nsh1 + 846);
    const auto *nsh1_849 = buffer.data(nsh1 + 849);
    const auto *nsh1_850 = buffer.data(nsh1 + 850);
    const auto *nsh1_852 = buffer.data(nsh1 + 852);
    const auto *nsh1_854 = buffer.data(nsh1 + 854);
    const auto *nsh1_855 = buffer.data(nsh1 + 855);
    const auto *nsh1_857 = buffer.data(nsh1 + 857);
    const auto *nsh1_858 = buffer.data(nsh1 + 858);
    const auto *nsh1_859 = buffer.data(nsh1 + 859);
    const auto *nsh1_860 = buffer.data(nsh1 + 860);
    const auto *nsh1_861 = buffer.data(nsh1 + 861);
    const auto *nsh1_864 = buffer.data(nsh1 + 864);

    const auto *nsi_1069 = buffer.data(nsi + 1069);
    const auto *nsi_1070 = buffer.data(nsi + 1070);
    const auto *nsi_1073 = buffer.data(nsi + 1073);
    const auto *nsi_1074 = buffer.data(nsi + 1074);
    const auto *nsi_1076 = buffer.data(nsi + 1076);
    const auto *nsi_1078 = buffer.data(nsi + 1078);
    const auto *nsi_1079 = buffer.data(nsi + 1079);
    const auto *nsi_1081 = buffer.data(nsi + 1081);
    const auto *nsi_1082 = buffer.data(nsi + 1082);
    const auto *nsi_1084 = buffer.data(nsi + 1084);
    const auto *nsi_1085 = buffer.data(nsi + 1085);
    const auto *nsi_1086 = buffer.data(nsi + 1086);
    const auto *nsi_1087 = buffer.data(nsi + 1087);
    const auto *nsi_1088 = buffer.data(nsi + 1088);
    const auto *nsi_1089 = buffer.data(nsi + 1089);
    const auto *nsi_1090 = buffer.data(nsi + 1090);
    const auto *nsi_1091 = buffer.data(nsi + 1091);
    const auto *nsi_1092 = buffer.data(nsi + 1092);
    const auto *nsi_1094 = buffer.data(nsi + 1094);
    const auto *nsi_1095 = buffer.data(nsi + 1095);
    const auto *nsi_1097 = buffer.data(nsi + 1097);
    const auto *nsi_1098 = buffer.data(nsi + 1098);
    const auto *nsi_1101 = buffer.data(nsi + 1101);
    const auto *nsi_1102 = buffer.data(nsi + 1102);
    const auto *nsi_1104 = buffer.data(nsi + 1104);
    const auto *nsi_1106 = buffer.data(nsi + 1106);
    const auto *nsi_1107 = buffer.data(nsi + 1107);
    const auto *nsi_1109 = buffer.data(nsi + 1109);
    const auto *nsi_1110 = buffer.data(nsi + 1110);
    const auto *nsi_1112 = buffer.data(nsi + 1112);
    const auto *nsi_1113 = buffer.data(nsi + 1113);
    const auto *nsi_1114 = buffer.data(nsi + 1114);
    const auto *nsi_1115 = buffer.data(nsi + 1115);
    const auto *nsi_1116 = buffer.data(nsi + 1116);
    const auto *nsi_1117 = buffer.data(nsi + 1117);
    const auto *nsi_1118 = buffer.data(nsi + 1118);
    const auto *nsi_1119 = buffer.data(nsi + 1119);
    const auto *nsi_1120 = buffer.data(nsi + 1120);
    const auto *nsi_1122 = buffer.data(nsi + 1122);
    const auto *nsi_1123 = buffer.data(nsi + 1123);
    const auto *nsi_1125 = buffer.data(nsi + 1125);
    const auto *nsi_1126 = buffer.data(nsi + 1126);
    const auto *nsi_1129 = buffer.data(nsi + 1129);
    const auto *nsi_1130 = buffer.data(nsi + 1130);
    const auto *nsi_1132 = buffer.data(nsi + 1132);
    const auto *nsi_1134 = buffer.data(nsi + 1134);
    const auto *nsi_1135 = buffer.data(nsi + 1135);
    const auto *nsi_1137 = buffer.data(nsi + 1137);
    const auto *nsi_1138 = buffer.data(nsi + 1138);
    const auto *nsi_1140 = buffer.data(nsi + 1140);
    const auto *nsi_1141 = buffer.data(nsi + 1141);
    const auto *nsi_1142 = buffer.data(nsi + 1142);
    const auto *nsi_1143 = buffer.data(nsi + 1143);
    const auto *nsi_1144 = buffer.data(nsi + 1144);
    const auto *nsi_1145 = buffer.data(nsi + 1145);
    const auto *nsi_1146 = buffer.data(nsi + 1146);
    const auto *nsi_1147 = buffer.data(nsi + 1147);
    const auto *nsi_1148 = buffer.data(nsi + 1148);
    const auto *nsi_1150 = buffer.data(nsi + 1150);
    const auto *nsi_1151 = buffer.data(nsi + 1151);

#pragma omp simd aligned(t_1376, t_1377, t_1378, pc_x, pc_y, msi_845, msi_1073, msi_1074, \
                         nsh0_807, nsh0_808, nsh1_807, nsh1_808, nsi_1069, nsi_1073, \
                         nsi_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_23 * msi_845[k]
                    + f_3 * pc_y[k] * nsi_1069[k];

        t_1377[k] = f_14 * msi_1073[k]
                    + f_8 * nsh0_807[k]
                    - f_9 * nsh1_807[k]
                    + f_3 * pc_x[k] * nsi_1073[k];

        t_1378[k] = f_14 * msi_1074[k]
                    + f_6 * nsh0_808[k]
                    - f_7 * nsh1_808[k]
                    + f_3 * pc_x[k] * nsi_1074[k];
    }

#pragma omp simd aligned(t_1379, t_1380, t_1381, pc_x, pc_y, pc_z, msi_818, msi_849, msi_1076, \
                         nsh0_810, nsh1_810, nsi_1070, nsi_1073, \
                         nsi_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1379[k] = f_14 * msi_818[k]
                    + f_3 * pc_z[k] * nsi_1070[k];

        t_1380[k] = f_14 * msi_1076[k]
                    + f_6 * nsh0_810[k]
                    - f_7 * nsh1_810[k]
                    + f_3 * pc_x[k] * nsi_1076[k];

        t_1381[k] = f_23 * msi_849[k]
                    + f_3 * pc_y[k] * nsi_1073[k];
    }

#pragma omp simd aligned(t_1382, t_1383, t_1384, pc_x, pc_z, msi_822, msi_1078, msi_1079, \
                         nsh0_812, nsh0_813, nsh1_812, nsh1_813, nsi_1074, nsi_1078, \
                         nsi_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1382[k] = f_14 * msi_1078[k]
                    + f_6 * nsh0_812[k]
                    - f_7 * nsh1_812[k]
                    + f_3 * pc_x[k] * nsi_1078[k];

        t_1383[k] = f_14 * msi_1079[k]
                    + f_4 * nsh0_813[k]
                    - f_5 * nsh1_813[k]
                    + f_3 * pc_x[k] * nsi_1079[k];

        t_1384[k] = f_14 * msi_822[k]
                    + f_3 * pc_z[k] * nsi_1074[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, pc_x, pc_y, msi_854, msi_1081, msi_1082, \
                         nsh0_815, nsh0_816, nsh1_815, nsh1_816, nsi_1078, nsi_1081, \
                         nsi_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_14 * msi_1081[k]
                    + f_4 * nsh0_815[k]
                    - f_5 * nsh1_815[k]
                    + f_3 * pc_x[k] * nsi_1081[k];

        t_1386[k] = f_14 * msi_1082[k]
                    + f_4 * nsh0_816[k]
                    - f_5 * nsh1_816[k]
                    + f_3 * pc_x[k] * nsi_1082[k];

        t_1387[k] = f_23 * msi_854[k]
                    + f_3 * pc_y[k] * nsi_1078[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, t_1391, pc_x, msi_1084, msi_1085, msi_1086, \
                         msi_1087, nsh0_818, nsh1_818, nsi_1084, nsi_1085, nsi_1086, \
                         nsi_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_14 * msi_1084[k]
                    + f_4 * nsh0_818[k]
                    - f_5 * nsh1_818[k]
                    + f_3 * pc_x[k] * nsi_1084[k];

        t_1389[k] = f_14 * msi_1085[k]
                    + f_3 * pc_x[k] * nsi_1085[k];

        t_1390[k] = f_14 * msi_1086[k]
                    + f_3 * pc_x[k] * nsi_1086[k];

        t_1391[k] = f_14 * msi_1087[k]
                    + f_3 * pc_x[k] * nsi_1087[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, t_1395, pc_x, msi_1088, msi_1089, msi_1090, \
                         msi_1091, nsi_1088, nsi_1089, nsi_1090, \
                         nsi_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_14 * msi_1088[k]
                    + f_3 * pc_x[k] * nsi_1088[k];

        t_1393[k] = f_14 * msi_1089[k]
                    + f_3 * pc_x[k] * nsi_1089[k];

        t_1394[k] = f_14 * msi_1090[k]
                    + f_3 * pc_x[k] * nsi_1090[k];

        t_1395[k] = f_14 * msi_1091[k]
                    + f_3 * pc_x[k] * nsi_1091[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, pc_y, pc_z, msi_833, msi_861, msi_863, \
                         nsh0_813, nsh0_815, nsh1_813, nsh1_815, nsi_1085, \
                         nsi_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_23 * msi_861[k]
                    + f_1 * nsh0_813[k]
                    - f_2 * nsh1_813[k]
                    + f_3 * pc_y[k] * nsi_1085[k];

        t_1397[k] = f_14 * msi_833[k]
                    + f_3 * pc_z[k] * nsi_1085[k];

        t_1398[k] = f_23 * msi_863[k]
                    + f_10 * nsh0_815[k]
                    - f_11 * nsh1_815[k]
                    + f_3 * pc_y[k] * nsi_1087[k];
    }

#pragma omp simd aligned(t_1399, t_1400, t_1401, pc_y, msi_864, msi_865, msi_866, nsh0_816, \
                         nsh0_817, nsh0_818, nsh1_816, nsh1_817, nsh1_818, nsi_1088, nsi_1089, \
                         nsi_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1399[k] = f_23 * msi_864[k]
                    + f_8 * nsh0_816[k]
                    - f_9 * nsh1_816[k]
                    + f_3 * pc_y[k] * nsi_1088[k];

        t_1400[k] = f_23 * msi_865[k]
                    + f_6 * nsh0_817[k]
                    - f_7 * nsh1_817[k]
                    + f_3 * pc_y[k] * nsi_1089[k];

        t_1401[k] = f_23 * msi_866[k]
                    + f_4 * nsh0_818[k]
                    - f_5 * nsh1_818[k]
                    + f_3 * pc_y[k] * nsi_1090[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, pc_x, pc_y, pc_z, msi_839, msi_867, msi_1092, \
                         nsh0_818, nsh0_819, nsh1_818, nsh1_819, nsi_1091, \
                         nsi_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_23 * msi_867[k]
                    + f_3 * pc_y[k] * nsi_1091[k];

        t_1403[k] = f_14 * msi_839[k]
                    + f_1 * nsh0_818[k]
                    - f_2 * nsh1_818[k]
                    + f_3 * pc_z[k] * nsi_1091[k];

        t_1404[k] = f_14 * msi_1092[k]
                    + f_1 * nsh0_819[k]
                    - f_2 * nsh1_819[k]
                    + f_3 * pc_x[k] * nsi_1092[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, t_1408, pc_x, pc_y, pc_z, msi_840, msi_868, \
                         msi_870, msi_1095, nsh0_822, nsh1_822, nsi_1092, nsi_1094, \
                         nsi_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_17 * msi_868[k]
                    + f_3 * pc_y[k] * nsi_1092[k];

        t_1406[k] = f_15 * msi_840[k]
                    + f_3 * pc_z[k] * nsi_1092[k];

        t_1407[k] = f_14 * msi_1095[k]
                    + f_10 * nsh0_822[k]
                    - f_11 * nsh1_822[k]
                    + f_3 * pc_x[k] * nsi_1095[k];

        t_1408[k] = f_17 * msi_870[k]
                    + f_3 * pc_y[k] * nsi_1094[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pc_x, pc_z, msi_843, msi_1097, msi_1098, \
                         nsh0_824, nsh0_825, nsh1_824, nsh1_825, nsi_1095, nsi_1097, \
                         nsi_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_14 * msi_1097[k]
                    + f_10 * nsh0_824[k]
                    - f_11 * nsh1_824[k]
                    + f_3 * pc_x[k] * nsi_1097[k];

        t_1410[k] = f_14 * msi_1098[k]
                    + f_8 * nsh0_825[k]
                    - f_9 * nsh1_825[k]
                    + f_3 * pc_x[k] * nsi_1098[k];

        t_1411[k] = f_15 * msi_843[k]
                    + f_3 * pc_z[k] * nsi_1095[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pc_x, pc_y, msi_873, msi_1101, msi_1102, \
                         nsh0_828, nsh0_829, nsh1_828, nsh1_829, nsi_1097, nsi_1101, \
                         nsi_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_17 * msi_873[k]
                    + f_3 * pc_y[k] * nsi_1097[k];

        t_1413[k] = f_14 * msi_1101[k]
                    + f_8 * nsh0_828[k]
                    - f_9 * nsh1_828[k]
                    + f_3 * pc_x[k] * nsi_1101[k];

        t_1414[k] = f_14 * msi_1102[k]
                    + f_6 * nsh0_829[k]
                    - f_7 * nsh1_829[k]
                    + f_3 * pc_x[k] * nsi_1102[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pc_x, pc_y, pc_z, msi_846, msi_877, msi_1104, \
                         nsh0_831, nsh1_831, nsi_1098, nsi_1101, \
                         nsi_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_15 * msi_846[k]
                    + f_3 * pc_z[k] * nsi_1098[k];

        t_1416[k] = f_14 * msi_1104[k]
                    + f_6 * nsh0_831[k]
                    - f_7 * nsh1_831[k]
                    + f_3 * pc_x[k] * nsi_1104[k];

        t_1417[k] = f_17 * msi_877[k]
                    + f_3 * pc_y[k] * nsi_1101[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pc_x, pc_z, msi_850, msi_1106, msi_1107, \
                         nsh0_833, nsh0_834, nsh1_833, nsh1_834, nsi_1102, nsi_1106, \
                         nsi_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_14 * msi_1106[k]
                    + f_6 * nsh0_833[k]
                    - f_7 * nsh1_833[k]
                    + f_3 * pc_x[k] * nsi_1106[k];

        t_1419[k] = f_14 * msi_1107[k]
                    + f_4 * nsh0_834[k]
                    - f_5 * nsh1_834[k]
                    + f_3 * pc_x[k] * nsi_1107[k];

        t_1420[k] = f_15 * msi_850[k]
                    + f_3 * pc_z[k] * nsi_1102[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, pc_x, pc_y, msi_882, msi_1109, msi_1110, \
                         nsh0_836, nsh0_837, nsh1_836, nsh1_837, nsi_1106, nsi_1109, \
                         nsi_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_14 * msi_1109[k]
                    + f_4 * nsh0_836[k]
                    - f_5 * nsh1_836[k]
                    + f_3 * pc_x[k] * nsi_1109[k];

        t_1422[k] = f_14 * msi_1110[k]
                    + f_4 * nsh0_837[k]
                    - f_5 * nsh1_837[k]
                    + f_3 * pc_x[k] * nsi_1110[k];

        t_1423[k] = f_17 * msi_882[k]
                    + f_3 * pc_y[k] * nsi_1106[k];
    }

#pragma omp simd aligned(t_1424, t_1425, t_1426, t_1427, pc_x, msi_1112, msi_1113, msi_1114, \
                         msi_1115, nsh0_839, nsh1_839, nsi_1112, nsi_1113, nsi_1114, \
                         nsi_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1424[k] = f_14 * msi_1112[k]
                    + f_4 * nsh0_839[k]
                    - f_5 * nsh1_839[k]
                    + f_3 * pc_x[k] * nsi_1112[k];

        t_1425[k] = f_14 * msi_1113[k]
                    + f_3 * pc_x[k] * nsi_1113[k];

        t_1426[k] = f_14 * msi_1114[k]
                    + f_3 * pc_x[k] * nsi_1114[k];

        t_1427[k] = f_14 * msi_1115[k]
                    + f_3 * pc_x[k] * nsi_1115[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, t_1431, pc_x, msi_1116, msi_1117, msi_1118, \
                         msi_1119, nsi_1116, nsi_1117, nsi_1118, \
                         nsi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_14 * msi_1116[k]
                    + f_3 * pc_x[k] * nsi_1116[k];

        t_1429[k] = f_14 * msi_1117[k]
                    + f_3 * pc_x[k] * nsi_1117[k];

        t_1430[k] = f_14 * msi_1118[k]
                    + f_3 * pc_x[k] * nsi_1118[k];

        t_1431[k] = f_14 * msi_1119[k]
                    + f_3 * pc_x[k] * nsi_1119[k];
    }

#pragma omp simd aligned(t_1432, t_1433, t_1434, pc_y, pc_z, msi_861, msi_889, msi_891, \
                         nsh0_834, nsh0_836, nsh1_834, nsh1_836, nsi_1113, \
                         nsi_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1432[k] = f_17 * msi_889[k]
                    + f_1 * nsh0_834[k]
                    - f_2 * nsh1_834[k]
                    + f_3 * pc_y[k] * nsi_1113[k];

        t_1433[k] = f_15 * msi_861[k]
                    + f_3 * pc_z[k] * nsi_1113[k];

        t_1434[k] = f_17 * msi_891[k]
                    + f_10 * nsh0_836[k]
                    - f_11 * nsh1_836[k]
                    + f_3 * pc_y[k] * nsi_1115[k];
    }

#pragma omp simd aligned(t_1435, t_1436, t_1437, pc_y, msi_892, msi_893, msi_894, nsh0_837, \
                         nsh0_838, nsh0_839, nsh1_837, nsh1_838, nsh1_839, nsi_1116, nsi_1117, \
                         nsi_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1435[k] = f_17 * msi_892[k]
                    + f_8 * nsh0_837[k]
                    - f_9 * nsh1_837[k]
                    + f_3 * pc_y[k] * nsi_1116[k];

        t_1436[k] = f_17 * msi_893[k]
                    + f_6 * nsh0_838[k]
                    - f_7 * nsh1_838[k]
                    + f_3 * pc_y[k] * nsi_1117[k];

        t_1437[k] = f_17 * msi_894[k]
                    + f_4 * nsh0_839[k]
                    - f_5 * nsh1_839[k]
                    + f_3 * pc_y[k] * nsi_1118[k];
    }

#pragma omp simd aligned(t_1438, t_1439, t_1440, pc_x, pc_y, pc_z, msi_867, msi_895, msi_1120, \
                         nsh0_839, nsh0_840, nsh1_839, nsh1_840, nsi_1119, \
                         nsi_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1438[k] = f_17 * msi_895[k]
                    + f_3 * pc_y[k] * nsi_1119[k];

        t_1439[k] = f_15 * msi_867[k]
                    + f_1 * nsh0_839[k]
                    - f_2 * nsh1_839[k]
                    + f_3 * pc_z[k] * nsi_1119[k];

        t_1440[k] = f_14 * msi_1120[k]
                    + f_1 * nsh0_840[k]
                    - f_2 * nsh1_840[k]
                    + f_3 * pc_x[k] * nsi_1120[k];
    }

#pragma omp simd aligned(t_1441, t_1442, t_1443, t_1444, pc_x, pc_y, pc_z, msi_868, msi_896, \
                         msi_898, msi_1123, nsh0_843, nsh1_843, nsi_1120, nsi_1122, \
                         nsi_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1441[k] = f_16 * msi_896[k]
                    + f_3 * pc_y[k] * nsi_1120[k];

        t_1442[k] = f_16 * msi_868[k]
                    + f_3 * pc_z[k] * nsi_1120[k];

        t_1443[k] = f_14 * msi_1123[k]
                    + f_10 * nsh0_843[k]
                    - f_11 * nsh1_843[k]
                    + f_3 * pc_x[k] * nsi_1123[k];

        t_1444[k] = f_16 * msi_898[k]
                    + f_3 * pc_y[k] * nsi_1122[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_x, pc_z, msi_871, msi_1125, msi_1126, \
                         nsh0_845, nsh0_846, nsh1_845, nsh1_846, nsi_1123, nsi_1125, \
                         nsi_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_14 * msi_1125[k]
                    + f_10 * nsh0_845[k]
                    - f_11 * nsh1_845[k]
                    + f_3 * pc_x[k] * nsi_1125[k];

        t_1446[k] = f_14 * msi_1126[k]
                    + f_8 * nsh0_846[k]
                    - f_9 * nsh1_846[k]
                    + f_3 * pc_x[k] * nsi_1126[k];

        t_1447[k] = f_16 * msi_871[k]
                    + f_3 * pc_z[k] * nsi_1123[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, msi_901, msi_1129, msi_1130, \
                         nsh0_849, nsh0_850, nsh1_849, nsh1_850, nsi_1125, nsi_1129, \
                         nsi_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_16 * msi_901[k]
                    + f_3 * pc_y[k] * nsi_1125[k];

        t_1449[k] = f_14 * msi_1129[k]
                    + f_8 * nsh0_849[k]
                    - f_9 * nsh1_849[k]
                    + f_3 * pc_x[k] * nsi_1129[k];

        t_1450[k] = f_14 * msi_1130[k]
                    + f_6 * nsh0_850[k]
                    - f_7 * nsh1_850[k]
                    + f_3 * pc_x[k] * nsi_1130[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, pc_y, pc_z, msi_874, msi_905, msi_1132, \
                         nsh0_852, nsh1_852, nsi_1126, nsi_1129, \
                         nsi_1132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_16 * msi_874[k]
                    + f_3 * pc_z[k] * nsi_1126[k];

        t_1452[k] = f_14 * msi_1132[k]
                    + f_6 * nsh0_852[k]
                    - f_7 * nsh1_852[k]
                    + f_3 * pc_x[k] * nsi_1132[k];

        t_1453[k] = f_16 * msi_905[k]
                    + f_3 * pc_y[k] * nsi_1129[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_z, msi_878, msi_1134, msi_1135, \
                         nsh0_854, nsh0_855, nsh1_854, nsh1_855, nsi_1130, nsi_1134, \
                         nsi_1135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_14 * msi_1134[k]
                    + f_6 * nsh0_854[k]
                    - f_7 * nsh1_854[k]
                    + f_3 * pc_x[k] * nsi_1134[k];

        t_1455[k] = f_14 * msi_1135[k]
                    + f_4 * nsh0_855[k]
                    - f_5 * nsh1_855[k]
                    + f_3 * pc_x[k] * nsi_1135[k];

        t_1456[k] = f_16 * msi_878[k]
                    + f_3 * pc_z[k] * nsi_1130[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, pc_y, msi_910, msi_1137, msi_1138, \
                         nsh0_857, nsh0_858, nsh1_857, nsh1_858, nsi_1134, nsi_1137, \
                         nsi_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_14 * msi_1137[k]
                    + f_4 * nsh0_857[k]
                    - f_5 * nsh1_857[k]
                    + f_3 * pc_x[k] * nsi_1137[k];

        t_1458[k] = f_14 * msi_1138[k]
                    + f_4 * nsh0_858[k]
                    - f_5 * nsh1_858[k]
                    + f_3 * pc_x[k] * nsi_1138[k];

        t_1459[k] = f_16 * msi_910[k]
                    + f_3 * pc_y[k] * nsi_1134[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, t_1463, pc_x, msi_1140, msi_1141, msi_1142, \
                         msi_1143, nsh0_860, nsh1_860, nsi_1140, nsi_1141, nsi_1142, \
                         nsi_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_14 * msi_1140[k]
                    + f_4 * nsh0_860[k]
                    - f_5 * nsh1_860[k]
                    + f_3 * pc_x[k] * nsi_1140[k];

        t_1461[k] = f_14 * msi_1141[k]
                    + f_3 * pc_x[k] * nsi_1141[k];

        t_1462[k] = f_14 * msi_1142[k]
                    + f_3 * pc_x[k] * nsi_1142[k];

        t_1463[k] = f_14 * msi_1143[k]
                    + f_3 * pc_x[k] * nsi_1143[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, t_1467, pc_x, msi_1144, msi_1145, msi_1146, \
                         msi_1147, nsi_1144, nsi_1145, nsi_1146, \
                         nsi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_14 * msi_1144[k]
                    + f_3 * pc_x[k] * nsi_1144[k];

        t_1465[k] = f_14 * msi_1145[k]
                    + f_3 * pc_x[k] * nsi_1145[k];

        t_1466[k] = f_14 * msi_1146[k]
                    + f_3 * pc_x[k] * nsi_1146[k];

        t_1467[k] = f_14 * msi_1147[k]
                    + f_3 * pc_x[k] * nsi_1147[k];
    }

#pragma omp simd aligned(t_1468, t_1469, t_1470, pc_y, pc_z, msi_889, msi_917, msi_919, \
                         nsh0_855, nsh0_857, nsh1_855, nsh1_857, nsi_1141, \
                         nsi_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1468[k] = f_16 * msi_917[k]
                    + f_1 * nsh0_855[k]
                    - f_2 * nsh1_855[k]
                    + f_3 * pc_y[k] * nsi_1141[k];

        t_1469[k] = f_16 * msi_889[k]
                    + f_3 * pc_z[k] * nsi_1141[k];

        t_1470[k] = f_16 * msi_919[k]
                    + f_10 * nsh0_857[k]
                    - f_11 * nsh1_857[k]
                    + f_3 * pc_y[k] * nsi_1143[k];
    }

#pragma omp simd aligned(t_1471, t_1472, t_1473, pc_y, msi_920, msi_921, msi_922, nsh0_858, \
                         nsh0_859, nsh0_860, nsh1_858, nsh1_859, nsh1_860, nsi_1144, nsi_1145, \
                         nsi_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1471[k] = f_16 * msi_920[k]
                    + f_8 * nsh0_858[k]
                    - f_9 * nsh1_858[k]
                    + f_3 * pc_y[k] * nsi_1144[k];

        t_1472[k] = f_16 * msi_921[k]
                    + f_6 * nsh0_859[k]
                    - f_7 * nsh1_859[k]
                    + f_3 * pc_y[k] * nsi_1145[k];

        t_1473[k] = f_16 * msi_922[k]
                    + f_4 * nsh0_860[k]
                    - f_5 * nsh1_860[k]
                    + f_3 * pc_y[k] * nsi_1146[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, pc_x, pc_y, pc_z, msi_895, msi_923, msi_1148, \
                         nsh0_860, nsh0_861, nsh1_860, nsh1_861, nsi_1147, \
                         nsi_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_16 * msi_923[k]
                    + f_3 * pc_y[k] * nsi_1147[k];

        t_1475[k] = f_16 * msi_895[k]
                    + f_1 * nsh0_860[k]
                    - f_2 * nsh1_860[k]
                    + f_3 * pc_z[k] * nsi_1147[k];

        t_1476[k] = f_14 * msi_1148[k]
                    + f_1 * nsh0_861[k]
                    - f_2 * nsh1_861[k]
                    + f_3 * pc_x[k] * nsi_1148[k];
    }

#pragma omp simd aligned(t_1477, t_1478, t_1479, t_1480, pc_x, pc_y, pc_z, msi_896, msi_924, \
                         msi_926, msi_1151, nsh0_864, nsh1_864, nsi_1148, nsi_1150, \
                         nsi_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1477[k] = f_15 * msi_924[k]
                    + f_3 * pc_y[k] * nsi_1148[k];

        t_1478[k] = f_17 * msi_896[k]
                    + f_3 * pc_z[k] * nsi_1148[k];

        t_1479[k] = f_14 * msi_1151[k]
                    + f_10 * nsh0_864[k]
                    - f_11 * nsh1_864[k]
                    + f_3 * pc_x[k] * nsi_1151[k];

        t_1480[k] = f_15 * msi_926[k]
                    + f_3 * pc_y[k] * nsi_1150[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsh0, const size_t nsh1,
                                                           const size_t nsi, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1260 = buffer.data(msk0 + 1260);
    const auto *msk0_1263 = buffer.data(msk0 + 1263);
    const auto *msk0_1265 = buffer.data(msk0 + 1265);
    const auto *msk0_1266 = buffer.data(msk0 + 1266);
    const auto *msk0_1269 = buffer.data(msk0 + 1269);
    const auto *msk0_1270 = buffer.data(msk0 + 1270);
    const auto *msk0_1272 = buffer.data(msk0 + 1272);
    const auto *msk0_1274 = buffer.data(msk0 + 1274);
    const auto *msk0_1275 = buffer.data(msk0 + 1275);
    const auto *msk0_1277 = buffer.data(msk0 + 1277);
    const auto *msk0_1278 = buffer.data(msk0 + 1278);
    const auto *msk0_1280 = buffer.data(msk0 + 1280);
    const auto *msk0_1295 = buffer.data(msk0 + 1295);

    const auto *msi_899 = buffer.data(msi + 899);
    const auto *msi_902 = buffer.data(msi + 902);
    const auto *msi_906 = buffer.data(msi + 906);
    const auto *msi_917 = buffer.data(msi + 917);
    const auto *msi_923 = buffer.data(msi + 923);
    const auto *msi_924 = buffer.data(msi + 924);
    const auto *msi_927 = buffer.data(msi + 927);
    const auto *msi_929 = buffer.data(msi + 929);
    const auto *msi_930 = buffer.data(msi + 930);
    const auto *msi_933 = buffer.data(msi + 933);
    const auto *msi_934 = buffer.data(msi + 934);
    const auto *msi_938 = buffer.data(msi + 938);
    const auto *msi_945 = buffer.data(msi + 945);
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
    const auto *msi_975 = buffer.data(msi + 975);
    const auto *msi_976 = buffer.data(msi + 976);
    const auto *msi_977 = buffer.data(msi + 977);
    const auto *msi_978 = buffer.data(msi + 978);
    const auto *msi_979 = buffer.data(msi + 979);
    const auto *msi_980 = buffer.data(msi + 980);
    const auto *msi_981 = buffer.data(msi + 981);
    const auto *msi_982 = buffer.data(msi + 982);
    const auto *msi_983 = buffer.data(msi + 983);
    const auto *msi_985 = buffer.data(msi + 985);
    const auto *msi_986 = buffer.data(msi + 986);
    const auto *msi_988 = buffer.data(msi + 988);
    const auto *msi_989 = buffer.data(msi + 989);
    const auto *msi_990 = buffer.data(msi + 990);
    const auto *msi_992 = buffer.data(msi + 992);
    const auto *msi_993 = buffer.data(msi + 993);
    const auto *msi_994 = buffer.data(msi + 994);
    const auto *msi_1001 = buffer.data(msi + 1001);
    const auto *msi_1003 = buffer.data(msi + 1003);
    const auto *msi_1004 = buffer.data(msi + 1004);
    const auto *msi_1005 = buffer.data(msi + 1005);
    const auto *msi_1006 = buffer.data(msi + 1006);
    const auto *msi_1007 = buffer.data(msi + 1007);
    const auto *msi_1153 = buffer.data(msi + 1153);
    const auto *msi_1154 = buffer.data(msi + 1154);
    const auto *msi_1157 = buffer.data(msi + 1157);
    const auto *msi_1158 = buffer.data(msi + 1158);
    const auto *msi_1160 = buffer.data(msi + 1160);
    const auto *msi_1162 = buffer.data(msi + 1162);
    const auto *msi_1163 = buffer.data(msi + 1163);
    const auto *msi_1165 = buffer.data(msi + 1165);
    const auto *msi_1166 = buffer.data(msi + 1166);
    const auto *msi_1168 = buffer.data(msi + 1168);
    const auto *msi_1169 = buffer.data(msi + 1169);
    const auto *msi_1170 = buffer.data(msi + 1170);
    const auto *msi_1171 = buffer.data(msi + 1171);
    const auto *msi_1172 = buffer.data(msi + 1172);
    const auto *msi_1173 = buffer.data(msi + 1173);
    const auto *msi_1174 = buffer.data(msi + 1174);
    const auto *msi_1175 = buffer.data(msi + 1175);
    const auto *msi_1176 = buffer.data(msi + 1176);
    const auto *msi_1179 = buffer.data(msi + 1179);
    const auto *msi_1181 = buffer.data(msi + 1181);
    const auto *msi_1182 = buffer.data(msi + 1182);
    const auto *msi_1185 = buffer.data(msi + 1185);
    const auto *msi_1186 = buffer.data(msi + 1186);
    const auto *msi_1188 = buffer.data(msi + 1188);
    const auto *msi_1190 = buffer.data(msi + 1190);
    const auto *msi_1191 = buffer.data(msi + 1191);
    const auto *msi_1193 = buffer.data(msi + 1193);
    const auto *msi_1194 = buffer.data(msi + 1194);
    const auto *msi_1196 = buffer.data(msi + 1196);
    const auto *msi_1197 = buffer.data(msi + 1197);
    const auto *msi_1198 = buffer.data(msi + 1198);
    const auto *msi_1199 = buffer.data(msi + 1199);
    const auto *msi_1200 = buffer.data(msi + 1200);
    const auto *msi_1201 = buffer.data(msi + 1201);
    const auto *msi_1202 = buffer.data(msi + 1202);
    const auto *msi_1203 = buffer.data(msi + 1203);
    const auto *msi_1225 = buffer.data(msi + 1225);
    const auto *msi_1226 = buffer.data(msi + 1226);
    const auto *msi_1227 = buffer.data(msi + 1227);
    const auto *msi_1228 = buffer.data(msi + 1228);
    const auto *msi_1229 = buffer.data(msi + 1229);
    const auto *msi_1230 = buffer.data(msi + 1230);
    const auto *msi_1231 = buffer.data(msi + 1231);
    const auto *msi_1232 = buffer.data(msi + 1232);

    const auto *msk1_1260 = buffer.data(msk1 + 1260);
    const auto *msk1_1263 = buffer.data(msk1 + 1263);
    const auto *msk1_1265 = buffer.data(msk1 + 1265);
    const auto *msk1_1266 = buffer.data(msk1 + 1266);
    const auto *msk1_1269 = buffer.data(msk1 + 1269);
    const auto *msk1_1270 = buffer.data(msk1 + 1270);
    const auto *msk1_1272 = buffer.data(msk1 + 1272);
    const auto *msk1_1274 = buffer.data(msk1 + 1274);
    const auto *msk1_1275 = buffer.data(msk1 + 1275);
    const auto *msk1_1277 = buffer.data(msk1 + 1277);
    const auto *msk1_1278 = buffer.data(msk1 + 1278);
    const auto *msk1_1280 = buffer.data(msk1 + 1280);
    const auto *msk1_1295 = buffer.data(msk1 + 1295);

    const auto *nsh0_866 = buffer.data(nsh0 + 866);
    const auto *nsh0_867 = buffer.data(nsh0 + 867);
    const auto *nsh0_870 = buffer.data(nsh0 + 870);
    const auto *nsh0_871 = buffer.data(nsh0 + 871);
    const auto *nsh0_873 = buffer.data(nsh0 + 873);
    const auto *nsh0_875 = buffer.data(nsh0 + 875);
    const auto *nsh0_876 = buffer.data(nsh0 + 876);
    const auto *nsh0_878 = buffer.data(nsh0 + 878);
    const auto *nsh0_879 = buffer.data(nsh0 + 879);
    const auto *nsh0_880 = buffer.data(nsh0 + 880);
    const auto *nsh0_881 = buffer.data(nsh0 + 881);
    const auto *nsh0_882 = buffer.data(nsh0 + 882);
    const auto *nsh0_885 = buffer.data(nsh0 + 885);
    const auto *nsh0_887 = buffer.data(nsh0 + 887);
    const auto *nsh0_888 = buffer.data(nsh0 + 888);
    const auto *nsh0_891 = buffer.data(nsh0 + 891);
    const auto *nsh0_892 = buffer.data(nsh0 + 892);
    const auto *nsh0_894 = buffer.data(nsh0 + 894);
    const auto *nsh0_896 = buffer.data(nsh0 + 896);
    const auto *nsh0_897 = buffer.data(nsh0 + 897);
    const auto *nsh0_899 = buffer.data(nsh0 + 899);
    const auto *nsh0_900 = buffer.data(nsh0 + 900);
    const auto *nsh0_901 = buffer.data(nsh0 + 901);
    const auto *nsh0_902 = buffer.data(nsh0 + 902);
    const auto *nsh0_918 = buffer.data(nsh0 + 918);
    const auto *nsh0_920 = buffer.data(nsh0 + 920);
    const auto *nsh0_921 = buffer.data(nsh0 + 921);
    const auto *nsh0_922 = buffer.data(nsh0 + 922);
    const auto *nsh0_923 = buffer.data(nsh0 + 923);
    const auto *nsh0_924 = buffer.data(nsh0 + 924);

    const auto *nsh1_866 = buffer.data(nsh1 + 866);
    const auto *nsh1_867 = buffer.data(nsh1 + 867);
    const auto *nsh1_870 = buffer.data(nsh1 + 870);
    const auto *nsh1_871 = buffer.data(nsh1 + 871);
    const auto *nsh1_873 = buffer.data(nsh1 + 873);
    const auto *nsh1_875 = buffer.data(nsh1 + 875);
    const auto *nsh1_876 = buffer.data(nsh1 + 876);
    const auto *nsh1_878 = buffer.data(nsh1 + 878);
    const auto *nsh1_879 = buffer.data(nsh1 + 879);
    const auto *nsh1_880 = buffer.data(nsh1 + 880);
    const auto *nsh1_881 = buffer.data(nsh1 + 881);
    const auto *nsh1_882 = buffer.data(nsh1 + 882);
    const auto *nsh1_885 = buffer.data(nsh1 + 885);
    const auto *nsh1_887 = buffer.data(nsh1 + 887);
    const auto *nsh1_888 = buffer.data(nsh1 + 888);
    const auto *nsh1_891 = buffer.data(nsh1 + 891);
    const auto *nsh1_892 = buffer.data(nsh1 + 892);
    const auto *nsh1_894 = buffer.data(nsh1 + 894);
    const auto *nsh1_896 = buffer.data(nsh1 + 896);
    const auto *nsh1_897 = buffer.data(nsh1 + 897);
    const auto *nsh1_899 = buffer.data(nsh1 + 899);
    const auto *nsh1_900 = buffer.data(nsh1 + 900);
    const auto *nsh1_901 = buffer.data(nsh1 + 901);
    const auto *nsh1_902 = buffer.data(nsh1 + 902);
    const auto *nsh1_918 = buffer.data(nsh1 + 918);
    const auto *nsh1_920 = buffer.data(nsh1 + 920);
    const auto *nsh1_921 = buffer.data(nsh1 + 921);
    const auto *nsh1_922 = buffer.data(nsh1 + 922);
    const auto *nsh1_923 = buffer.data(nsh1 + 923);
    const auto *nsh1_924 = buffer.data(nsh1 + 924);

    const auto *nsi_1151 = buffer.data(nsi + 1151);
    const auto *nsi_1153 = buffer.data(nsi + 1153);
    const auto *nsi_1154 = buffer.data(nsi + 1154);
    const auto *nsi_1157 = buffer.data(nsi + 1157);
    const auto *nsi_1158 = buffer.data(nsi + 1158);
    const auto *nsi_1160 = buffer.data(nsi + 1160);
    const auto *nsi_1162 = buffer.data(nsi + 1162);
    const auto *nsi_1163 = buffer.data(nsi + 1163);
    const auto *nsi_1165 = buffer.data(nsi + 1165);
    const auto *nsi_1166 = buffer.data(nsi + 1166);
    const auto *nsi_1168 = buffer.data(nsi + 1168);
    const auto *nsi_1169 = buffer.data(nsi + 1169);
    const auto *nsi_1170 = buffer.data(nsi + 1170);
    const auto *nsi_1171 = buffer.data(nsi + 1171);
    const auto *nsi_1172 = buffer.data(nsi + 1172);
    const auto *nsi_1173 = buffer.data(nsi + 1173);
    const auto *nsi_1174 = buffer.data(nsi + 1174);
    const auto *nsi_1175 = buffer.data(nsi + 1175);
    const auto *nsi_1176 = buffer.data(nsi + 1176);
    const auto *nsi_1178 = buffer.data(nsi + 1178);
    const auto *nsi_1179 = buffer.data(nsi + 1179);
    const auto *nsi_1181 = buffer.data(nsi + 1181);
    const auto *nsi_1182 = buffer.data(nsi + 1182);
    const auto *nsi_1185 = buffer.data(nsi + 1185);
    const auto *nsi_1186 = buffer.data(nsi + 1186);
    const auto *nsi_1188 = buffer.data(nsi + 1188);
    const auto *nsi_1190 = buffer.data(nsi + 1190);
    const auto *nsi_1191 = buffer.data(nsi + 1191);
    const auto *nsi_1193 = buffer.data(nsi + 1193);
    const auto *nsi_1194 = buffer.data(nsi + 1194);
    const auto *nsi_1196 = buffer.data(nsi + 1196);
    const auto *nsi_1197 = buffer.data(nsi + 1197);
    const auto *nsi_1198 = buffer.data(nsi + 1198);
    const auto *nsi_1199 = buffer.data(nsi + 1199);
    const auto *nsi_1200 = buffer.data(nsi + 1200);
    const auto *nsi_1201 = buffer.data(nsi + 1201);
    const auto *nsi_1202 = buffer.data(nsi + 1202);
    const auto *nsi_1203 = buffer.data(nsi + 1203);
    const auto *nsi_1204 = buffer.data(nsi + 1204);
    const auto *nsi_1206 = buffer.data(nsi + 1206);
    const auto *nsi_1207 = buffer.data(nsi + 1207);
    const auto *nsi_1209 = buffer.data(nsi + 1209);
    const auto *nsi_1210 = buffer.data(nsi + 1210);
    const auto *nsi_1213 = buffer.data(nsi + 1213);
    const auto *nsi_1214 = buffer.data(nsi + 1214);
    const auto *nsi_1218 = buffer.data(nsi + 1218);
    const auto *nsi_1225 = buffer.data(nsi + 1225);
    const auto *nsi_1226 = buffer.data(nsi + 1226);
    const auto *nsi_1227 = buffer.data(nsi + 1227);
    const auto *nsi_1228 = buffer.data(nsi + 1228);
    const auto *nsi_1229 = buffer.data(nsi + 1229);
    const auto *nsi_1230 = buffer.data(nsi + 1230);
    const auto *nsi_1231 = buffer.data(nsi + 1231);
    const auto *nsi_1232 = buffer.data(nsi + 1232);
    const auto *nsi_1233 = buffer.data(nsi + 1233);
    const auto *nsi_1234 = buffer.data(nsi + 1234);

#pragma omp simd aligned(t_1481, t_1482, t_1483, pc_x, pc_z, msi_899, msi_1153, msi_1154, \
                         nsh0_866, nsh0_867, nsh1_866, nsh1_867, nsi_1151, nsi_1153, \
                         nsi_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_14 * msi_1153[k]
                    + f_10 * nsh0_866[k]
                    - f_11 * nsh1_866[k]
                    + f_3 * pc_x[k] * nsi_1153[k];

        t_1482[k] = f_14 * msi_1154[k]
                    + f_8 * nsh0_867[k]
                    - f_9 * nsh1_867[k]
                    + f_3 * pc_x[k] * nsi_1154[k];

        t_1483[k] = f_17 * msi_899[k]
                    + f_3 * pc_z[k] * nsi_1151[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, pc_x, pc_y, msi_929, msi_1157, msi_1158, \
                         nsh0_870, nsh0_871, nsh1_870, nsh1_871, nsi_1153, nsi_1157, \
                         nsi_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_15 * msi_929[k]
                    + f_3 * pc_y[k] * nsi_1153[k];

        t_1485[k] = f_14 * msi_1157[k]
                    + f_8 * nsh0_870[k]
                    - f_9 * nsh1_870[k]
                    + f_3 * pc_x[k] * nsi_1157[k];

        t_1486[k] = f_14 * msi_1158[k]
                    + f_6 * nsh0_871[k]
                    - f_7 * nsh1_871[k]
                    + f_3 * pc_x[k] * nsi_1158[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pc_x, pc_y, pc_z, msi_902, msi_933, msi_1160, \
                         nsh0_873, nsh1_873, nsi_1154, nsi_1157, \
                         nsi_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_17 * msi_902[k]
                    + f_3 * pc_z[k] * nsi_1154[k];

        t_1488[k] = f_14 * msi_1160[k]
                    + f_6 * nsh0_873[k]
                    - f_7 * nsh1_873[k]
                    + f_3 * pc_x[k] * nsi_1160[k];

        t_1489[k] = f_15 * msi_933[k]
                    + f_3 * pc_y[k] * nsi_1157[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pc_x, pc_z, msi_906, msi_1162, msi_1163, \
                         nsh0_875, nsh0_876, nsh1_875, nsh1_876, nsi_1158, nsi_1162, \
                         nsi_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_14 * msi_1162[k]
                    + f_6 * nsh0_875[k]
                    - f_7 * nsh1_875[k]
                    + f_3 * pc_x[k] * nsi_1162[k];

        t_1491[k] = f_14 * msi_1163[k]
                    + f_4 * nsh0_876[k]
                    - f_5 * nsh1_876[k]
                    + f_3 * pc_x[k] * nsi_1163[k];

        t_1492[k] = f_17 * msi_906[k]
                    + f_3 * pc_z[k] * nsi_1158[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pc_x, pc_y, msi_938, msi_1165, msi_1166, \
                         nsh0_878, nsh0_879, nsh1_878, nsh1_879, nsi_1162, nsi_1165, \
                         nsi_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_14 * msi_1165[k]
                    + f_4 * nsh0_878[k]
                    - f_5 * nsh1_878[k]
                    + f_3 * pc_x[k] * nsi_1165[k];

        t_1494[k] = f_14 * msi_1166[k]
                    + f_4 * nsh0_879[k]
                    - f_5 * nsh1_879[k]
                    + f_3 * pc_x[k] * nsi_1166[k];

        t_1495[k] = f_15 * msi_938[k]
                    + f_3 * pc_y[k] * nsi_1162[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pc_x, msi_1168, msi_1169, msi_1170, \
                         msi_1171, nsh0_881, nsh1_881, nsi_1168, nsi_1169, nsi_1170, \
                         nsi_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_14 * msi_1168[k]
                    + f_4 * nsh0_881[k]
                    - f_5 * nsh1_881[k]
                    + f_3 * pc_x[k] * nsi_1168[k];

        t_1497[k] = f_14 * msi_1169[k]
                    + f_3 * pc_x[k] * nsi_1169[k];

        t_1498[k] = f_14 * msi_1170[k]
                    + f_3 * pc_x[k] * nsi_1170[k];

        t_1499[k] = f_14 * msi_1171[k]
                    + f_3 * pc_x[k] * nsi_1171[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, pc_x, msi_1172, msi_1173, msi_1174, \
                         msi_1175, nsi_1172, nsi_1173, nsi_1174, \
                         nsi_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_14 * msi_1172[k]
                    + f_3 * pc_x[k] * nsi_1172[k];

        t_1501[k] = f_14 * msi_1173[k]
                    + f_3 * pc_x[k] * nsi_1173[k];

        t_1502[k] = f_14 * msi_1174[k]
                    + f_3 * pc_x[k] * nsi_1174[k];

        t_1503[k] = f_14 * msi_1175[k]
                    + f_3 * pc_x[k] * nsi_1175[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, pc_y, pc_z, msi_917, msi_945, msi_947, \
                         nsh0_876, nsh0_878, nsh1_876, nsh1_878, nsi_1169, \
                         nsi_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = f_15 * msi_945[k]
                    + f_1 * nsh0_876[k]
                    - f_2 * nsh1_876[k]
                    + f_3 * pc_y[k] * nsi_1169[k];

        t_1505[k] = f_17 * msi_917[k]
                    + f_3 * pc_z[k] * nsi_1169[k];

        t_1506[k] = f_15 * msi_947[k]
                    + f_10 * nsh0_878[k]
                    - f_11 * nsh1_878[k]
                    + f_3 * pc_y[k] * nsi_1171[k];
    }

#pragma omp simd aligned(t_1507, t_1508, t_1509, pc_y, msi_948, msi_949, msi_950, nsh0_879, \
                         nsh0_880, nsh0_881, nsh1_879, nsh1_880, nsh1_881, nsi_1172, nsi_1173, \
                         nsi_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1507[k] = f_15 * msi_948[k]
                    + f_8 * nsh0_879[k]
                    - f_9 * nsh1_879[k]
                    + f_3 * pc_y[k] * nsi_1172[k];

        t_1508[k] = f_15 * msi_949[k]
                    + f_6 * nsh0_880[k]
                    - f_7 * nsh1_880[k]
                    + f_3 * pc_y[k] * nsi_1173[k];

        t_1509[k] = f_15 * msi_950[k]
                    + f_4 * nsh0_881[k]
                    - f_5 * nsh1_881[k]
                    + f_3 * pc_y[k] * nsi_1174[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, pc_x, pc_y, pc_z, msi_923, msi_951, msi_1176, \
                         nsh0_881, nsh0_882, nsh1_881, nsh1_882, nsi_1175, \
                         nsi_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_15 * msi_951[k]
                    + f_3 * pc_y[k] * nsi_1175[k];

        t_1511[k] = f_17 * msi_923[k]
                    + f_1 * nsh0_881[k]
                    - f_2 * nsh1_881[k]
                    + f_3 * pc_z[k] * nsi_1175[k];

        t_1512[k] = f_14 * msi_1176[k]
                    + f_1 * nsh0_882[k]
                    - f_2 * nsh1_882[k]
                    + f_3 * pc_x[k] * nsi_1176[k];
    }

#pragma omp simd aligned(t_1513, t_1514, t_1515, t_1516, pc_x, pc_y, pc_z, msi_924, msi_952, \
                         msi_954, msi_1179, nsh0_885, nsh1_885, nsi_1176, nsi_1178, \
                         nsi_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1513[k] = f_14 * msi_952[k]
                    + f_3 * pc_y[k] * nsi_1176[k];

        t_1514[k] = f_23 * msi_924[k]
                    + f_3 * pc_z[k] * nsi_1176[k];

        t_1515[k] = f_14 * msi_1179[k]
                    + f_10 * nsh0_885[k]
                    - f_11 * nsh1_885[k]
                    + f_3 * pc_x[k] * nsi_1179[k];

        t_1516[k] = f_14 * msi_954[k]
                    + f_3 * pc_y[k] * nsi_1178[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, pc_x, pc_z, msi_927, msi_1181, msi_1182, \
                         nsh0_887, nsh0_888, nsh1_887, nsh1_888, nsi_1179, nsi_1181, \
                         nsi_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_14 * msi_1181[k]
                    + f_10 * nsh0_887[k]
                    - f_11 * nsh1_887[k]
                    + f_3 * pc_x[k] * nsi_1181[k];

        t_1518[k] = f_14 * msi_1182[k]
                    + f_8 * nsh0_888[k]
                    - f_9 * nsh1_888[k]
                    + f_3 * pc_x[k] * nsi_1182[k];

        t_1519[k] = f_23 * msi_927[k]
                    + f_3 * pc_z[k] * nsi_1179[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pc_x, pc_y, msi_957, msi_1185, msi_1186, \
                         nsh0_891, nsh0_892, nsh1_891, nsh1_892, nsi_1181, nsi_1185, \
                         nsi_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_14 * msi_957[k]
                    + f_3 * pc_y[k] * nsi_1181[k];

        t_1521[k] = f_14 * msi_1185[k]
                    + f_8 * nsh0_891[k]
                    - f_9 * nsh1_891[k]
                    + f_3 * pc_x[k] * nsi_1185[k];

        t_1522[k] = f_14 * msi_1186[k]
                    + f_6 * nsh0_892[k]
                    - f_7 * nsh1_892[k]
                    + f_3 * pc_x[k] * nsi_1186[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_x, pc_y, pc_z, msi_930, msi_961, msi_1188, \
                         nsh0_894, nsh1_894, nsi_1182, nsi_1185, \
                         nsi_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_23 * msi_930[k]
                    + f_3 * pc_z[k] * nsi_1182[k];

        t_1524[k] = f_14 * msi_1188[k]
                    + f_6 * nsh0_894[k]
                    - f_7 * nsh1_894[k]
                    + f_3 * pc_x[k] * nsi_1188[k];

        t_1525[k] = f_14 * msi_961[k]
                    + f_3 * pc_y[k] * nsi_1185[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pc_x, pc_z, msi_934, msi_1190, msi_1191, \
                         nsh0_896, nsh0_897, nsh1_896, nsh1_897, nsi_1186, nsi_1190, \
                         nsi_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_14 * msi_1190[k]
                    + f_6 * nsh0_896[k]
                    - f_7 * nsh1_896[k]
                    + f_3 * pc_x[k] * nsi_1190[k];

        t_1527[k] = f_14 * msi_1191[k]
                    + f_4 * nsh0_897[k]
                    - f_5 * nsh1_897[k]
                    + f_3 * pc_x[k] * nsi_1191[k];

        t_1528[k] = f_23 * msi_934[k]
                    + f_3 * pc_z[k] * nsi_1186[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, pc_x, pc_y, msi_966, msi_1193, msi_1194, \
                         nsh0_899, nsh0_900, nsh1_899, nsh1_900, nsi_1190, nsi_1193, \
                         nsi_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_14 * msi_1193[k]
                    + f_4 * nsh0_899[k]
                    - f_5 * nsh1_899[k]
                    + f_3 * pc_x[k] * nsi_1193[k];

        t_1530[k] = f_14 * msi_1194[k]
                    + f_4 * nsh0_900[k]
                    - f_5 * nsh1_900[k]
                    + f_3 * pc_x[k] * nsi_1194[k];

        t_1531[k] = f_14 * msi_966[k]
                    + f_3 * pc_y[k] * nsi_1190[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, t_1535, pc_x, msi_1196, msi_1197, msi_1198, \
                         msi_1199, nsh0_902, nsh1_902, nsi_1196, nsi_1197, nsi_1198, \
                         nsi_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = f_14 * msi_1196[k]
                    + f_4 * nsh0_902[k]
                    - f_5 * nsh1_902[k]
                    + f_3 * pc_x[k] * nsi_1196[k];

        t_1533[k] = f_14 * msi_1197[k]
                    + f_3 * pc_x[k] * nsi_1197[k];

        t_1534[k] = f_14 * msi_1198[k]
                    + f_3 * pc_x[k] * nsi_1198[k];

        t_1535[k] = f_14 * msi_1199[k]
                    + f_3 * pc_x[k] * nsi_1199[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pc_x, msi_1200, msi_1201, msi_1202, \
                         msi_1203, nsi_1200, nsi_1201, nsi_1202, \
                         nsi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_14 * msi_1200[k]
                    + f_3 * pc_x[k] * nsi_1200[k];

        t_1537[k] = f_14 * msi_1201[k]
                    + f_3 * pc_x[k] * nsi_1201[k];

        t_1538[k] = f_14 * msi_1202[k]
                    + f_3 * pc_x[k] * nsi_1202[k];

        t_1539[k] = f_14 * msi_1203[k]
                    + f_3 * pc_x[k] * nsi_1203[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, pc_y, pc_z, msi_945, msi_973, msi_975, \
                         nsh0_897, nsh0_899, nsh1_897, nsh1_899, nsi_1197, \
                         nsi_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = f_14 * msi_973[k]
                    + f_1 * nsh0_897[k]
                    - f_2 * nsh1_897[k]
                    + f_3 * pc_y[k] * nsi_1197[k];

        t_1541[k] = f_23 * msi_945[k]
                    + f_3 * pc_z[k] * nsi_1197[k];

        t_1542[k] = f_14 * msi_975[k]
                    + f_10 * nsh0_899[k]
                    - f_11 * nsh1_899[k]
                    + f_3 * pc_y[k] * nsi_1199[k];
    }

#pragma omp simd aligned(t_1543, t_1544, t_1545, pc_y, msi_976, msi_977, msi_978, nsh0_900, \
                         nsh0_901, nsh0_902, nsh1_900, nsh1_901, nsh1_902, nsi_1200, nsi_1201, \
                         nsi_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1543[k] = f_14 * msi_976[k]
                    + f_8 * nsh0_900[k]
                    - f_9 * nsh1_900[k]
                    + f_3 * pc_y[k] * nsi_1200[k];

        t_1544[k] = f_14 * msi_977[k]
                    + f_6 * nsh0_901[k]
                    - f_7 * nsh1_901[k]
                    + f_3 * pc_y[k] * nsi_1201[k];

        t_1545[k] = f_14 * msi_978[k]
                    + f_4 * nsh0_902[k]
                    - f_5 * nsh1_902[k]
                    + f_3 * pc_y[k] * nsi_1202[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pa_y, pc_y, pc_z, msk0_1260, msi_951, \
                         msi_979, msi_980, msk1_1260, nsh0_902, nsh1_902, nsi_1203, \
                         nsi_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_14 * msi_979[k]
                    + f_3 * pc_y[k] * nsi_1203[k];

        t_1547[k] = f_23 * msi_951[k]
                    + f_1 * nsh0_902[k]
                    - f_2 * nsh1_902[k]
                    + f_3 * pc_z[k] * nsi_1203[k];

        t_1548[k] = pa_y[k] * msk0_1260[k]
                    - f_12 * pc_y[k] * msk1_1260[k];

        t_1549[k] = f_13 * msi_980[k]
                    + f_3 * pc_y[k] * nsi_1204[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, t_1553, pa_y, pc_y, pc_z, msk0_1263, \
                         msk0_1265, msi_952, msi_981, msi_982, msk1_1263, msk1_1265, nsi_1204, \
                         nsi_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_22 * msi_952[k]
                    + f_3 * pc_z[k] * nsi_1204[k];

        t_1551[k] = pa_y[k] * msk0_1263[k]
                    + f_14 * msi_981[k]
                    - f_12 * pc_y[k] * msk1_1263[k];

        t_1552[k] = f_13 * msi_982[k]
                    + f_3 * pc_y[k] * nsi_1206[k];

        t_1553[k] = pa_y[k] * msk0_1265[k]
                    - f_12 * pc_y[k] * msk1_1265[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, t_1557, pa_y, pc_y, pc_z, msk0_1266, \
                         msk0_1269, msi_955, msi_983, msi_985, msk1_1266, msk1_1269, nsi_1207, \
                         nsi_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = pa_y[k] * msk0_1266[k]
                    + f_15 * msi_983[k]
                    - f_12 * pc_y[k] * msk1_1266[k];

        t_1555[k] = f_22 * msi_955[k]
                    + f_3 * pc_z[k] * nsi_1207[k];

        t_1556[k] = f_13 * msi_985[k]
                    + f_3 * pc_y[k] * nsi_1209[k];

        t_1557[k] = pa_y[k] * msk0_1269[k]
                    - f_12 * pc_y[k] * msk1_1269[k];
    }

#pragma omp simd aligned(t_1558, t_1559, t_1560, pa_y, pc_y, pc_z, msk0_1270, msk0_1272, \
                         msi_958, msi_986, msi_988, msk1_1270, msk1_1272, \
                         nsi_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1558[k] = pa_y[k] * msk0_1270[k]
                    + f_16 * msi_986[k]
                    - f_12 * pc_y[k] * msk1_1270[k];

        t_1559[k] = f_22 * msi_958[k]
                    + f_3 * pc_z[k] * nsi_1210[k];

        t_1560[k] = pa_y[k] * msk0_1272[k]
                    + f_14 * msi_988[k]
                    - f_12 * pc_y[k] * msk1_1272[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, pa_y, pc_y, pc_z, msk0_1274, \
                         msk0_1275, msi_962, msi_989, msi_990, msk1_1274, msk1_1275, nsi_1213, \
                         nsi_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = f_13 * msi_989[k]
                    + f_3 * pc_y[k] * nsi_1213[k];

        t_1562[k] = pa_y[k] * msk0_1274[k]
                    - f_12 * pc_y[k] * msk1_1274[k];

        t_1563[k] = pa_y[k] * msk0_1275[k]
                    + f_17 * msi_990[k]
                    - f_12 * pc_y[k] * msk1_1275[k];

        t_1564[k] = f_22 * msi_962[k]
                    + f_3 * pc_z[k] * nsi_1214[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, t_1568, pa_y, pc_y, msk0_1277, msk0_1278, \
                         msk0_1280, msi_992, msi_993, msi_994, msk1_1277, msk1_1278, \
                         msk1_1280, nsi_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = pa_y[k] * msk0_1277[k]
                    + f_15 * msi_992[k]
                    - f_12 * pc_y[k] * msk1_1277[k];

        t_1566[k] = pa_y[k] * msk0_1278[k]
                    + f_14 * msi_993[k]
                    - f_12 * pc_y[k] * msk1_1278[k];

        t_1567[k] = f_13 * msi_994[k]
                    + f_3 * pc_y[k] * nsi_1218[k];

        t_1568[k] = pa_y[k] * msk0_1280[k]
                    - f_12 * pc_y[k] * msk1_1280[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, t_1572, t_1573, pc_x, msi_1225, msi_1226, \
                         msi_1227, msi_1228, msi_1229, nsi_1225, nsi_1226, nsi_1227, nsi_1228, \
                         nsi_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = f_14 * msi_1225[k]
                    + f_3 * pc_x[k] * nsi_1225[k];

        t_1570[k] = f_14 * msi_1226[k]
                    + f_3 * pc_x[k] * nsi_1226[k];

        t_1571[k] = f_14 * msi_1227[k]
                    + f_3 * pc_x[k] * nsi_1227[k];

        t_1572[k] = f_14 * msi_1228[k]
                    + f_3 * pc_x[k] * nsi_1228[k];

        t_1573[k] = f_14 * msi_1229[k]
                    + f_3 * pc_x[k] * nsi_1229[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pc_x, pc_y, pc_z, msi_973, msi_1001, \
                         msi_1230, msi_1231, nsh0_918, nsh1_918, nsi_1225, nsi_1230, \
                         nsi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = f_14 * msi_1230[k]
                    + f_3 * pc_x[k] * nsi_1230[k];

        t_1575[k] = f_14 * msi_1231[k]
                    + f_3 * pc_x[k] * nsi_1231[k];

        t_1576[k] = f_13 * msi_1001[k]
                    + f_1 * nsh0_918[k]
                    - f_2 * nsh1_918[k]
                    + f_3 * pc_y[k] * nsi_1225[k];

        t_1577[k] = f_22 * msi_973[k]
                    + f_3 * pc_z[k] * nsi_1225[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pc_y, msi_1003, msi_1004, msi_1005, nsh0_920, \
                         nsh0_921, nsh0_922, nsh1_920, nsh1_921, nsh1_922, nsi_1227, nsi_1228, \
                         nsi_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = f_13 * msi_1003[k]
                    + f_10 * nsh0_920[k]
                    - f_11 * nsh1_920[k]
                    + f_3 * pc_y[k] * nsi_1227[k];

        t_1579[k] = f_13 * msi_1004[k]
                    + f_8 * nsh0_921[k]
                    - f_9 * nsh1_921[k]
                    + f_3 * pc_y[k] * nsi_1228[k];

        t_1580[k] = f_13 * msi_1005[k]
                    + f_6 * nsh0_922[k]
                    - f_7 * nsh1_922[k]
                    + f_3 * pc_y[k] * nsi_1229[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pa_y, pc_y, msk0_1295, msi_1006, msi_1007, \
                         msk1_1295, nsh0_923, nsh1_923, nsi_1230, \
                         nsi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_13 * msi_1006[k]
                    + f_4 * nsh0_923[k]
                    - f_5 * nsh1_923[k]
                    + f_3 * pc_y[k] * nsi_1230[k];

        t_1582[k] = f_13 * msi_1007[k]
                    + f_3 * pc_y[k] * nsi_1231[k];

        t_1583[k] = pa_y[k] * msk0_1295[k]
                    - f_12 * pc_y[k] * msk1_1295[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, t_1587, t_1588, pc_x, pc_y, pc_z, msi_980, \
                         msi_1232, nsh0_924, nsh1_924, nsi_1232, nsi_1233, \
                         nsi_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_14 * msi_1232[k]
                    + f_1 * nsh0_924[k]
                    - f_2 * nsh1_924[k]
                    + f_3 * pc_x[k] * nsi_1232[k];

        t_1585[k] = f_3 * pc_y[k] * nsi_1232[k];

        t_1586[k] = f_21 * msi_980[k]
                    + f_3 * pc_z[k] * nsi_1232[k];

        t_1587[k] = f_4 * nsh0_924[k]
                    - f_5 * nsh1_924[k]
                    + f_3 * pc_y[k] * nsi_1233[k];

        t_1588[k] = f_3 * pc_y[k] * nsi_1234[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsh0, const size_t nsh1,
                                                           const size_t nsi, const size_t ncols,
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1296 = buffer.data(msk0 + 1296);
    const auto *msk0_1299 = buffer.data(msk0 + 1299);
    const auto *msk0_1302 = buffer.data(msk0 + 1302);
    const auto *msk0_1306 = buffer.data(msk0 + 1306);
    const auto *msk0_1311 = buffer.data(msk0 + 1311);
    const auto *msk0_1620 = buffer.data(msk0 + 1620);
    const auto *msk0_1623 = buffer.data(msk0 + 1623);
    const auto *msk0_1626 = buffer.data(msk0 + 1626);
    const auto *msk0_1630 = buffer.data(msk0 + 1630);
    const auto *msk0_1635 = buffer.data(msk0 + 1635);
    const auto *msk0_1648 = buffer.data(msk0 + 1648);
    const auto *msk0_1650 = buffer.data(msk0 + 1650);
    const auto *msk0_1651 = buffer.data(msk0 + 1651);
    const auto *msk0_1652 = buffer.data(msk0 + 1652);
    const auto *msk0_1653 = buffer.data(msk0 + 1653);
    const auto *msk0_1655 = buffer.data(msk0 + 1655);
    const auto *msk0_1661 = buffer.data(msk0 + 1661);
    const auto *msk0_1665 = buffer.data(msk0 + 1665);
    const auto *msk0_1668 = buffer.data(msk0 + 1668);
    const auto *msk0_1670 = buffer.data(msk0 + 1670);
    const auto *msk0_1673 = buffer.data(msk0 + 1673);
    const auto *msk0_1674 = buffer.data(msk0 + 1674);
    const auto *msk0_1676 = buffer.data(msk0 + 1676);
    const auto *msk0_1684 = buffer.data(msk0 + 1684);
    const auto *msk0_1686 = buffer.data(msk0 + 1686);
    const auto *msk0_1687 = buffer.data(msk0 + 1687);
    const auto *msk0_1688 = buffer.data(msk0 + 1688);
    const auto *msk0_1689 = buffer.data(msk0 + 1689);
    const auto *msk0_1691 = buffer.data(msk0 + 1691);
    const auto *msk0_1692 = buffer.data(msk0 + 1692);
    const auto *msk0_1695 = buffer.data(msk0 + 1695);
    const auto *msk0_1697 = buffer.data(msk0 + 1697);
    const auto *msk0_1698 = buffer.data(msk0 + 1698);
    const auto *msk0_1701 = buffer.data(msk0 + 1701);
    const auto *msk0_1702 = buffer.data(msk0 + 1702);
    const auto *msk0_1704 = buffer.data(msk0 + 1704);
    const auto *msk0_1706 = buffer.data(msk0 + 1706);
    const auto *msk0_1707 = buffer.data(msk0 + 1707);

    const auto *msi_1007 = buffer.data(msi + 1007);
    const auto *msi_1008 = buffer.data(msi + 1008);
    const auto *msi_1011 = buffer.data(msi + 1011);
    const auto *msi_1013 = buffer.data(msi + 1013);
    const auto *msi_1014 = buffer.data(msi + 1014);
    const auto *msi_1017 = buffer.data(msi + 1017);
    const auto *msi_1018 = buffer.data(msi + 1018);
    const auto *msi_1022 = buffer.data(msi + 1022);
    const auto *msi_1029 = buffer.data(msi + 1029);
    const auto *msi_1035 = buffer.data(msi + 1035);
    const auto *msi_1036 = buffer.data(msi + 1036);
    const auto *msi_1038 = buffer.data(msi + 1038);
    const auto *msi_1039 = buffer.data(msi + 1039);
    const auto *msi_1041 = buffer.data(msi + 1041);
    const auto *msi_1042 = buffer.data(msi + 1042);
    const auto *msi_1045 = buffer.data(msi + 1045);
    const auto *msi_1046 = buffer.data(msi + 1046);
    const auto *msi_1050 = buffer.data(msi + 1050);
    const auto *msi_1063 = buffer.data(msi + 1063);
    const auto *msi_1064 = buffer.data(msi + 1064);
    const auto *msi_1066 = buffer.data(msi + 1066);
    const auto *msi_1069 = buffer.data(msi + 1069);
    const auto *msi_1073 = buffer.data(msi + 1073);
    const auto *msi_1237 = buffer.data(msi + 1237);
    const auto *msi_1241 = buffer.data(msi + 1241);
    const auto *msi_1246 = buffer.data(msi + 1246);
    const auto *msi_1252 = buffer.data(msi + 1252);
    const auto *msi_1253 = buffer.data(msi + 1253);
    const auto *msi_1254 = buffer.data(msi + 1254);
    const auto *msi_1255 = buffer.data(msi + 1255);
    const auto *msi_1256 = buffer.data(msi + 1256);
    const auto *msi_1257 = buffer.data(msi + 1257);
    const auto *msi_1259 = buffer.data(msi + 1259);
    const auto *msi_1260 = buffer.data(msi + 1260);
    const auto *msi_1263 = buffer.data(msi + 1263);
    const auto *msi_1266 = buffer.data(msi + 1266);
    const auto *msi_1270 = buffer.data(msi + 1270);
    const auto *msi_1275 = buffer.data(msi + 1275);
    const auto *msi_1281 = buffer.data(msi + 1281);
    const auto *msi_1283 = buffer.data(msi + 1283);
    const auto *msi_1284 = buffer.data(msi + 1284);
    const auto *msi_1285 = buffer.data(msi + 1285);
    const auto *msi_1286 = buffer.data(msi + 1286);
    const auto *msi_1287 = buffer.data(msi + 1287);
    const auto *msi_1293 = buffer.data(msi + 1293);
    const auto *msi_1297 = buffer.data(msi + 1297);
    const auto *msi_1300 = buffer.data(msi + 1300);
    const auto *msi_1302 = buffer.data(msi + 1302);
    const auto *msi_1305 = buffer.data(msi + 1305);
    const auto *msi_1306 = buffer.data(msi + 1306);
    const auto *msi_1308 = buffer.data(msi + 1308);
    const auto *msi_1309 = buffer.data(msi + 1309);
    const auto *msi_1310 = buffer.data(msi + 1310);
    const auto *msi_1311 = buffer.data(msi + 1311);
    const auto *msi_1312 = buffer.data(msi + 1312);
    const auto *msi_1313 = buffer.data(msi + 1313);
    const auto *msi_1314 = buffer.data(msi + 1314);
    const auto *msi_1315 = buffer.data(msi + 1315);
    const auto *msi_1316 = buffer.data(msi + 1316);
    const auto *msi_1319 = buffer.data(msi + 1319);
    const auto *msi_1321 = buffer.data(msi + 1321);
    const auto *msi_1322 = buffer.data(msi + 1322);
    const auto *msi_1325 = buffer.data(msi + 1325);
    const auto *msi_1326 = buffer.data(msi + 1326);
    const auto *msi_1328 = buffer.data(msi + 1328);
    const auto *msi_1330 = buffer.data(msi + 1330);
    const auto *msi_1331 = buffer.data(msi + 1331);

    const auto *msk1_1296 = buffer.data(msk1 + 1296);
    const auto *msk1_1299 = buffer.data(msk1 + 1299);
    const auto *msk1_1302 = buffer.data(msk1 + 1302);
    const auto *msk1_1306 = buffer.data(msk1 + 1306);
    const auto *msk1_1311 = buffer.data(msk1 + 1311);
    const auto *msk1_1620 = buffer.data(msk1 + 1620);
    const auto *msk1_1623 = buffer.data(msk1 + 1623);
    const auto *msk1_1626 = buffer.data(msk1 + 1626);
    const auto *msk1_1630 = buffer.data(msk1 + 1630);
    const auto *msk1_1635 = buffer.data(msk1 + 1635);
    const auto *msk1_1648 = buffer.data(msk1 + 1648);
    const auto *msk1_1650 = buffer.data(msk1 + 1650);
    const auto *msk1_1651 = buffer.data(msk1 + 1651);
    const auto *msk1_1652 = buffer.data(msk1 + 1652);
    const auto *msk1_1653 = buffer.data(msk1 + 1653);
    const auto *msk1_1655 = buffer.data(msk1 + 1655);
    const auto *msk1_1661 = buffer.data(msk1 + 1661);
    const auto *msk1_1665 = buffer.data(msk1 + 1665);
    const auto *msk1_1668 = buffer.data(msk1 + 1668);
    const auto *msk1_1670 = buffer.data(msk1 + 1670);
    const auto *msk1_1673 = buffer.data(msk1 + 1673);
    const auto *msk1_1674 = buffer.data(msk1 + 1674);
    const auto *msk1_1676 = buffer.data(msk1 + 1676);
    const auto *msk1_1684 = buffer.data(msk1 + 1684);
    const auto *msk1_1686 = buffer.data(msk1 + 1686);
    const auto *msk1_1687 = buffer.data(msk1 + 1687);
    const auto *msk1_1688 = buffer.data(msk1 + 1688);
    const auto *msk1_1689 = buffer.data(msk1 + 1689);
    const auto *msk1_1691 = buffer.data(msk1 + 1691);
    const auto *msk1_1692 = buffer.data(msk1 + 1692);
    const auto *msk1_1695 = buffer.data(msk1 + 1695);
    const auto *msk1_1697 = buffer.data(msk1 + 1697);
    const auto *msk1_1698 = buffer.data(msk1 + 1698);
    const auto *msk1_1701 = buffer.data(msk1 + 1701);
    const auto *msk1_1702 = buffer.data(msk1 + 1702);
    const auto *msk1_1704 = buffer.data(msk1 + 1704);
    const auto *msk1_1706 = buffer.data(msk1 + 1706);
    const auto *msk1_1707 = buffer.data(msk1 + 1707);

    const auto *nsh0_925 = buffer.data(nsh0 + 925);
    const auto *nsh0_926 = buffer.data(nsh0 + 926);
    const auto *nsh0_927 = buffer.data(nsh0 + 927);
    const auto *nsh0_928 = buffer.data(nsh0 + 928);
    const auto *nsh0_929 = buffer.data(nsh0 + 929);
    const auto *nsh0_930 = buffer.data(nsh0 + 930);
    const auto *nsh0_931 = buffer.data(nsh0 + 931);
    const auto *nsh0_932 = buffer.data(nsh0 + 932);
    const auto *nsh0_933 = buffer.data(nsh0 + 933);
    const auto *nsh0_938 = buffer.data(nsh0 + 938);
    const auto *nsh0_939 = buffer.data(nsh0 + 939);
    const auto *nsh0_940 = buffer.data(nsh0 + 940);
    const auto *nsh0_941 = buffer.data(nsh0 + 941);
    const auto *nsh0_942 = buffer.data(nsh0 + 942);
    const auto *nsh0_943 = buffer.data(nsh0 + 943);
    const auto *nsh0_944 = buffer.data(nsh0 + 944);
    const auto *nsh0_945 = buffer.data(nsh0 + 945);
    const auto *nsh0_947 = buffer.data(nsh0 + 947);
    const auto *nsh0_948 = buffer.data(nsh0 + 948);
    const auto *nsh0_950 = buffer.data(nsh0 + 950);
    const auto *nsh0_951 = buffer.data(nsh0 + 951);
    const auto *nsh0_952 = buffer.data(nsh0 + 952);
    const auto *nsh0_954 = buffer.data(nsh0 + 954);

    const auto *nsh1_925 = buffer.data(nsh1 + 925);
    const auto *nsh1_926 = buffer.data(nsh1 + 926);
    const auto *nsh1_927 = buffer.data(nsh1 + 927);
    const auto *nsh1_928 = buffer.data(nsh1 + 928);
    const auto *nsh1_929 = buffer.data(nsh1 + 929);
    const auto *nsh1_930 = buffer.data(nsh1 + 930);
    const auto *nsh1_931 = buffer.data(nsh1 + 931);
    const auto *nsh1_932 = buffer.data(nsh1 + 932);
    const auto *nsh1_933 = buffer.data(nsh1 + 933);
    const auto *nsh1_938 = buffer.data(nsh1 + 938);
    const auto *nsh1_939 = buffer.data(nsh1 + 939);
    const auto *nsh1_940 = buffer.data(nsh1 + 940);
    const auto *nsh1_941 = buffer.data(nsh1 + 941);
    const auto *nsh1_942 = buffer.data(nsh1 + 942);
    const auto *nsh1_943 = buffer.data(nsh1 + 943);
    const auto *nsh1_944 = buffer.data(nsh1 + 944);
    const auto *nsh1_945 = buffer.data(nsh1 + 945);
    const auto *nsh1_947 = buffer.data(nsh1 + 947);
    const auto *nsh1_948 = buffer.data(nsh1 + 948);
    const auto *nsh1_950 = buffer.data(nsh1 + 950);
    const auto *nsh1_951 = buffer.data(nsh1 + 951);
    const auto *nsh1_952 = buffer.data(nsh1 + 952);
    const auto *nsh1_954 = buffer.data(nsh1 + 954);

    const auto *nsi_1235 = buffer.data(nsi + 1235);
    const auto *nsi_1236 = buffer.data(nsi + 1236);
    const auto *nsi_1237 = buffer.data(nsi + 1237);
    const auto *nsi_1238 = buffer.data(nsi + 1238);
    const auto *nsi_1239 = buffer.data(nsi + 1239);
    const auto *nsi_1240 = buffer.data(nsi + 1240);
    const auto *nsi_1241 = buffer.data(nsi + 1241);
    const auto *nsi_1242 = buffer.data(nsi + 1242);
    const auto *nsi_1243 = buffer.data(nsi + 1243);
    const auto *nsi_1244 = buffer.data(nsi + 1244);
    const auto *nsi_1245 = buffer.data(nsi + 1245);
    const auto *nsi_1246 = buffer.data(nsi + 1246);
    const auto *nsi_1252 = buffer.data(nsi + 1252);
    const auto *nsi_1253 = buffer.data(nsi + 1253);
    const auto *nsi_1254 = buffer.data(nsi + 1254);
    const auto *nsi_1255 = buffer.data(nsi + 1255);
    const auto *nsi_1256 = buffer.data(nsi + 1256);
    const auto *nsi_1257 = buffer.data(nsi + 1257);
    const auto *nsi_1258 = buffer.data(nsi + 1258);
    const auto *nsi_1259 = buffer.data(nsi + 1259);
    const auto *nsi_1260 = buffer.data(nsi + 1260);
    const auto *nsi_1261 = buffer.data(nsi + 1261);
    const auto *nsi_1262 = buffer.data(nsi + 1262);
    const auto *nsi_1263 = buffer.data(nsi + 1263);
    const auto *nsi_1265 = buffer.data(nsi + 1265);
    const auto *nsi_1266 = buffer.data(nsi + 1266);
    const auto *nsi_1267 = buffer.data(nsi + 1267);
    const auto *nsi_1269 = buffer.data(nsi + 1269);
    const auto *nsi_1270 = buffer.data(nsi + 1270);
    const auto *nsi_1271 = buffer.data(nsi + 1271);
    const auto *nsi_1272 = buffer.data(nsi + 1272);
    const auto *nsi_1274 = buffer.data(nsi + 1274);
    const auto *nsi_1275 = buffer.data(nsi + 1275);
    const auto *nsi_1281 = buffer.data(nsi + 1281);
    const auto *nsi_1283 = buffer.data(nsi + 1283);
    const auto *nsi_1284 = buffer.data(nsi + 1284);
    const auto *nsi_1285 = buffer.data(nsi + 1285);
    const auto *nsi_1286 = buffer.data(nsi + 1286);
    const auto *nsi_1287 = buffer.data(nsi + 1287);
    const auto *nsi_1288 = buffer.data(nsi + 1288);
    const auto *nsi_1290 = buffer.data(nsi + 1290);
    const auto *nsi_1291 = buffer.data(nsi + 1291);
    const auto *nsi_1293 = buffer.data(nsi + 1293);
    const auto *nsi_1294 = buffer.data(nsi + 1294);
    const auto *nsi_1297 = buffer.data(nsi + 1297);
    const auto *nsi_1298 = buffer.data(nsi + 1298);
    const auto *nsi_1302 = buffer.data(nsi + 1302);
    const auto *nsi_1309 = buffer.data(nsi + 1309);
    const auto *nsi_1310 = buffer.data(nsi + 1310);
    const auto *nsi_1311 = buffer.data(nsi + 1311);
    const auto *nsi_1312 = buffer.data(nsi + 1312);
    const auto *nsi_1313 = buffer.data(nsi + 1313);
    const auto *nsi_1314 = buffer.data(nsi + 1314);
    const auto *nsi_1315 = buffer.data(nsi + 1315);
    const auto *nsi_1316 = buffer.data(nsi + 1316);
    const auto *nsi_1318 = buffer.data(nsi + 1318);
    const auto *nsi_1319 = buffer.data(nsi + 1319);
    const auto *nsi_1321 = buffer.data(nsi + 1321);
    const auto *nsi_1322 = buffer.data(nsi + 1322);
    const auto *nsi_1325 = buffer.data(nsi + 1325);
    const auto *nsi_1326 = buffer.data(nsi + 1326);

#pragma omp simd aligned(t_1589, t_1590, t_1591, t_1592, pc_x, pc_y, msi_1237, nsh0_925, \
                         nsh0_926, nsh0_929, nsh1_925, nsh1_926, nsh1_929, nsi_1235, nsi_1236, \
                         nsi_1237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = f_14 * msi_1237[k]
                    + f_10 * nsh0_929[k]
                    - f_11 * nsh1_929[k]
                    + f_3 * pc_x[k] * nsi_1237[k];

        t_1590[k] = f_6 * nsh0_925[k]
                    - f_7 * nsh1_925[k]
                    + f_3 * pc_y[k] * nsi_1235[k];

        t_1591[k] = f_4 * nsh0_926[k]
                    - f_5 * nsh1_926[k]
                    + f_3 * pc_y[k] * nsi_1236[k];

        t_1592[k] = f_3 * pc_y[k] * nsi_1237[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pc_x, pc_y, msi_1241, nsh0_927, nsh0_928, \
                         nsh0_933, nsh1_927, nsh1_928, nsh1_933, nsi_1238, nsi_1239, \
                         nsi_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_14 * msi_1241[k]
                    + f_8 * nsh0_933[k]
                    - f_9 * nsh1_933[k]
                    + f_3 * pc_x[k] * nsi_1241[k];

        t_1594[k] = f_8 * nsh0_927[k]
                    - f_9 * nsh1_927[k]
                    + f_3 * pc_y[k] * nsi_1238[k];

        t_1595[k] = f_6 * nsh0_928[k]
                    - f_7 * nsh1_928[k]
                    + f_3 * pc_y[k] * nsi_1239[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, pc_x, pc_y, msi_1246, nsh0_929, nsh0_938, \
                         nsh1_929, nsh1_938, nsi_1240, nsi_1241, \
                         nsi_1246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_4 * nsh0_929[k]
                    - f_5 * nsh1_929[k]
                    + f_3 * pc_y[k] * nsi_1240[k];

        t_1597[k] = f_3 * pc_y[k] * nsi_1241[k];

        t_1598[k] = f_14 * msi_1246[k]
                    + f_6 * nsh0_938[k]
                    - f_7 * nsh1_938[k]
                    + f_3 * pc_x[k] * nsi_1246[k];
    }

#pragma omp simd aligned(t_1599, t_1600, t_1601, pc_y, nsh0_930, nsh0_931, nsh0_932, nsh1_930, \
                         nsh1_931, nsh1_932, nsi_1242, nsi_1243, \
                         nsi_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = f_10 * nsh0_930[k]
                    - f_11 * nsh1_930[k]
                    + f_3 * pc_y[k] * nsi_1242[k];

        t_1600[k] = f_8 * nsh0_931[k]
                    - f_9 * nsh1_931[k]
                    + f_3 * pc_y[k] * nsi_1243[k];

        t_1601[k] = f_6 * nsh0_932[k]
                    - f_7 * nsh1_932[k]
                    + f_3 * pc_y[k] * nsi_1244[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, t_1605, pc_x, pc_y, msi_1252, msi_1253, \
                         nsh0_933, nsh0_944, nsh1_933, nsh1_944, nsi_1245, nsi_1246, nsi_1252, \
                         nsi_1253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = f_4 * nsh0_933[k]
                    - f_5 * nsh1_933[k]
                    + f_3 * pc_y[k] * nsi_1245[k];

        t_1603[k] = f_3 * pc_y[k] * nsi_1246[k];

        t_1604[k] = f_14 * msi_1252[k]
                    + f_4 * nsh0_944[k]
                    - f_5 * nsh1_944[k]
                    + f_3 * pc_x[k] * nsi_1252[k];

        t_1605[k] = f_14 * msi_1253[k]
                    + f_3 * pc_x[k] * nsi_1253[k];
    }

#pragma omp simd aligned(t_1606, t_1607, t_1608, t_1609, t_1610, pc_x, pc_y, msi_1254, \
                         msi_1255, msi_1256, msi_1257, nsi_1252, nsi_1254, nsi_1255, nsi_1256, \
                         nsi_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1606[k] = f_14 * msi_1254[k]
                    + f_3 * pc_x[k] * nsi_1254[k];

        t_1607[k] = f_14 * msi_1255[k]
                    + f_3 * pc_x[k] * nsi_1255[k];

        t_1608[k] = f_14 * msi_1256[k]
                    + f_3 * pc_x[k] * nsi_1256[k];

        t_1609[k] = f_14 * msi_1257[k]
                    + f_3 * pc_x[k] * nsi_1257[k];

        t_1610[k] = f_3 * pc_y[k] * nsi_1252[k];
    }

#pragma omp simd aligned(t_1611, t_1612, t_1613, pc_x, pc_y, msi_1259, nsh0_939, nsh0_940, \
                         nsh1_939, nsh1_940, nsi_1253, nsi_1254, \
                         nsi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1611[k] = f_14 * msi_1259[k]
                    + f_3 * pc_x[k] * nsi_1259[k];

        t_1612[k] = f_1 * nsh0_939[k]
                    - f_2 * nsh1_939[k]
                    + f_3 * pc_y[k] * nsi_1253[k];

        t_1613[k] = f_19 * nsh0_940[k]
                    - f_20 * nsh1_940[k]
                    + f_3 * pc_y[k] * nsi_1254[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, pc_y, nsh0_941, nsh0_942, nsh0_943, nsh1_941, \
                         nsh1_942, nsh1_943, nsi_1255, nsi_1256, \
                         nsi_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_10 * nsh0_941[k]
                    - f_11 * nsh1_941[k]
                    + f_3 * pc_y[k] * nsi_1255[k];

        t_1615[k] = f_8 * nsh0_942[k]
                    - f_9 * nsh1_942[k]
                    + f_3 * pc_y[k] * nsi_1256[k];

        t_1616[k] = f_6 * nsh0_943[k]
                    - f_7 * nsh1_943[k]
                    + f_3 * pc_y[k] * nsi_1257[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, t_1620, pa_x, pc_x, pc_y, pc_z, msk0_1620, \
                         msi_1007, msi_1260, msk1_1620, nsh0_944, nsh1_944, nsi_1258, \
                         nsi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_4 * nsh0_944[k]
                    - f_5 * nsh1_944[k]
                    + f_3 * pc_y[k] * nsi_1258[k];

        t_1618[k] = f_3 * pc_y[k] * nsi_1259[k];

        t_1619[k] = f_21 * msi_1007[k]
                    + f_1 * nsh0_944[k]
                    - f_2 * nsh1_944[k]
                    + f_3 * pc_z[k] * nsi_1259[k];

        t_1620[k] = pa_x[k] * msk0_1620[k]
                    + f_22 * msi_1260[k]
                    - f_12 * pc_x[k] * msk1_1620[k];
    }

#pragma omp simd aligned(t_1621, t_1622, t_1623, t_1624, pa_x, pc_x, pc_y, pc_z, msk0_1623, \
                         msi_1008, msi_1263, msk1_1623, nsi_1260, \
                         nsi_1261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1621[k] = f_18 * msi_1008[k]
                    + f_3 * pc_y[k] * nsi_1260[k];

        t_1622[k] = f_3 * pc_z[k] * nsi_1260[k];

        t_1623[k] = pa_x[k] * msk0_1623[k]
                    + f_17 * msi_1263[k]
                    - f_12 * pc_x[k] * msk1_1623[k];

        t_1624[k] = f_3 * pc_z[k] * nsi_1261[k];
    }

#pragma omp simd aligned(t_1625, t_1626, t_1627, pa_x, pc_x, pc_z, msk0_1626, msi_1266, \
                         msk1_1626, nsh0_945, nsh1_945, nsi_1262, \
                         nsi_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = f_4 * nsh0_945[k]
                    - f_5 * nsh1_945[k]
                    + f_3 * pc_z[k] * nsi_1262[k];

        t_1626[k] = pa_x[k] * msk0_1626[k]
                    + f_16 * msi_1266[k]
                    - f_12 * pc_x[k] * msk1_1626[k];

        t_1627[k] = f_3 * pc_z[k] * nsi_1263[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, t_1631, pa_x, pc_x, pc_y, pc_z, msk0_1630, \
                         msi_1013, msi_1270, msk1_1630, nsh0_947, nsh1_947, nsi_1265, \
                         nsi_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_18 * msi_1013[k]
                    + f_3 * pc_y[k] * nsi_1265[k];

        t_1629[k] = f_6 * nsh0_947[k]
                    - f_7 * nsh1_947[k]
                    + f_3 * pc_z[k] * nsi_1265[k];

        t_1630[k] = pa_x[k] * msk0_1630[k]
                    + f_15 * msi_1270[k]
                    - f_12 * pc_x[k] * msk1_1630[k];

        t_1631[k] = f_3 * pc_z[k] * nsi_1266[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pc_y, pc_z, msi_1017, nsh0_948, nsh0_950, \
                         nsh1_948, nsh1_950, nsi_1267, nsi_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_4 * nsh0_948[k]
                    - f_5 * nsh1_948[k]
                    + f_3 * pc_z[k] * nsi_1267[k];

        t_1633[k] = f_18 * msi_1017[k]
                    + f_3 * pc_y[k] * nsi_1269[k];

        t_1634[k] = f_8 * nsh0_950[k]
                    - f_9 * nsh1_950[k]
                    + f_3 * pc_z[k] * nsi_1269[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, pa_x, pc_x, pc_z, msk0_1635, msi_1275, \
                         msk1_1635, nsh0_951, nsh1_951, nsi_1270, \
                         nsi_1271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = pa_x[k] * msk0_1635[k]
                    + f_14 * msi_1275[k]
                    - f_12 * pc_x[k] * msk1_1635[k];

        t_1636[k] = f_3 * pc_z[k] * nsi_1270[k];

        t_1637[k] = f_4 * nsh0_951[k]
                    - f_5 * nsh1_951[k]
                    + f_3 * pc_z[k] * nsi_1271[k];
    }

#pragma omp simd aligned(t_1638, t_1639, t_1640, t_1641, pc_x, pc_y, pc_z, msi_1022, msi_1281, \
                         nsh0_952, nsh0_954, nsh1_952, nsh1_954, nsi_1272, nsi_1274, \
                         nsi_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1638[k] = f_6 * nsh0_952[k]
                    - f_7 * nsh1_952[k]
                    + f_3 * pc_z[k] * nsi_1272[k];

        t_1639[k] = f_18 * msi_1022[k]
                    + f_3 * pc_y[k] * nsi_1274[k];

        t_1640[k] = f_10 * nsh0_954[k]
                    - f_11 * nsh1_954[k]
                    + f_3 * pc_z[k] * nsi_1274[k];

        t_1641[k] = f_13 * msi_1281[k]
                    + f_3 * pc_x[k] * nsi_1281[k];
    }

#pragma omp simd aligned(t_1642, t_1643, t_1644, t_1645, t_1646, pc_x, pc_z, msi_1283, \
                         msi_1284, msi_1285, msi_1286, nsi_1275, nsi_1283, nsi_1284, nsi_1285, \
                         nsi_1286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1642[k] = f_3 * pc_z[k] * nsi_1275[k];

        t_1643[k] = f_13 * msi_1283[k]
                    + f_3 * pc_x[k] * nsi_1283[k];

        t_1644[k] = f_13 * msi_1284[k]
                    + f_3 * pc_x[k] * nsi_1284[k];

        t_1645[k] = f_13 * msi_1285[k]
                    + f_3 * pc_x[k] * nsi_1285[k];

        t_1646[k] = f_13 * msi_1286[k]
                    + f_3 * pc_x[k] * nsi_1286[k];
    }

#pragma omp simd aligned(t_1647, t_1648, t_1649, t_1650, pa_x, pc_x, pc_z, msk0_1648, \
                         msk0_1650, msi_1287, msk1_1648, msk1_1650, nsi_1281, \
                         nsi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1647[k] = f_13 * msi_1287[k]
                    + f_3 * pc_x[k] * nsi_1287[k];

        t_1648[k] = pa_x[k] * msk0_1648[k]
                    - f_12 * pc_x[k] * msk1_1648[k];

        t_1649[k] = f_3 * pc_z[k] * nsi_1281[k];

        t_1650[k] = pa_x[k] * msk0_1650[k]
                    - f_12 * pc_x[k] * msk1_1650[k];
    }

#pragma omp simd aligned(t_1651, t_1652, t_1653, t_1654, pa_x, pc_x, pc_y, msk0_1651, \
                         msk0_1652, msk0_1653, msi_1035, msk1_1651, msk1_1652, msk1_1653, \
                         nsi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1651[k] = pa_x[k] * msk0_1651[k]
                    - f_12 * pc_x[k] * msk1_1651[k];

        t_1652[k] = pa_x[k] * msk0_1652[k]
                    - f_12 * pc_x[k] * msk1_1652[k];

        t_1653[k] = pa_x[k] * msk0_1653[k]
                    - f_12 * pc_x[k] * msk1_1653[k];

        t_1654[k] = f_18 * msi_1035[k]
                    + f_3 * pc_y[k] * nsi_1287[k];
    }

#pragma omp simd aligned(t_1655, t_1656, t_1657, t_1658, pa_x, pa_z, pc_x, pc_y, pc_z, \
                         msk0_1296, msk0_1655, msi_1008, msi_1036, msk1_1296, msk1_1655, \
                         nsi_1288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1655[k] = pa_x[k] * msk0_1655[k]
                    - f_12 * pc_x[k] * msk1_1655[k];

        t_1656[k] = pa_z[k] * msk0_1296[k]
                    - f_12 * pc_z[k] * msk1_1296[k];

        t_1657[k] = f_21 * msi_1036[k]
                    + f_3 * pc_y[k] * nsi_1288[k];

        t_1658[k] = f_13 * msi_1008[k]
                    + f_3 * pc_z[k] * nsi_1288[k];
    }

#pragma omp simd aligned(t_1659, t_1660, t_1661, pa_x, pa_z, pc_x, pc_y, pc_z, msk0_1299, \
                         msk0_1661, msi_1038, msi_1293, msk1_1299, msk1_1661, \
                         nsi_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1659[k] = pa_z[k] * msk0_1299[k]
                    - f_12 * pc_z[k] * msk1_1299[k];

        t_1660[k] = f_21 * msi_1038[k]
                    + f_3 * pc_y[k] * nsi_1290[k];

        t_1661[k] = pa_x[k] * msk0_1661[k]
                    + f_17 * msi_1293[k]
                    - f_12 * pc_x[k] * msk1_1661[k];
    }

#pragma omp simd aligned(t_1662, t_1663, t_1664, pa_z, pc_y, pc_z, msk0_1302, msi_1011, \
                         msi_1041, msk1_1302, nsi_1291, nsi_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1662[k] = pa_z[k] * msk0_1302[k]
                    - f_12 * pc_z[k] * msk1_1302[k];

        t_1663[k] = f_13 * msi_1011[k]
                    + f_3 * pc_z[k] * nsi_1291[k];

        t_1664[k] = f_21 * msi_1041[k]
                    + f_3 * pc_y[k] * nsi_1293[k];
    }

#pragma omp simd aligned(t_1665, t_1666, t_1667, pa_x, pa_z, pc_x, pc_z, msk0_1306, msk0_1665, \
                         msi_1014, msi_1297, msk1_1306, msk1_1665, \
                         nsi_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1665[k] = pa_x[k] * msk0_1665[k]
                    + f_16 * msi_1297[k]
                    - f_12 * pc_x[k] * msk1_1665[k];

        t_1666[k] = pa_z[k] * msk0_1306[k]
                    - f_12 * pc_z[k] * msk1_1306[k];

        t_1667[k] = f_13 * msi_1014[k]
                    + f_3 * pc_z[k] * nsi_1294[k];
    }

#pragma omp simd aligned(t_1668, t_1669, t_1670, pa_x, pc_x, pc_y, msk0_1668, msk0_1670, \
                         msi_1045, msi_1300, msi_1302, msk1_1668, msk1_1670, \
                         nsi_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1668[k] = pa_x[k] * msk0_1668[k]
                    + f_15 * msi_1300[k]
                    - f_12 * pc_x[k] * msk1_1668[k];

        t_1669[k] = f_21 * msi_1045[k]
                    + f_3 * pc_y[k] * nsi_1297[k];

        t_1670[k] = pa_x[k] * msk0_1670[k]
                    + f_15 * msi_1302[k]
                    - f_12 * pc_x[k] * msk1_1670[k];
    }

#pragma omp simd aligned(t_1671, t_1672, t_1673, pa_x, pa_z, pc_x, pc_z, msk0_1311, msk0_1673, \
                         msi_1018, msi_1305, msk1_1311, msk1_1673, \
                         nsi_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1671[k] = pa_z[k] * msk0_1311[k]
                    - f_12 * pc_z[k] * msk1_1311[k];

        t_1672[k] = f_13 * msi_1018[k]
                    + f_3 * pc_z[k] * nsi_1298[k];

        t_1673[k] = pa_x[k] * msk0_1673[k]
                    + f_14 * msi_1305[k]
                    - f_12 * pc_x[k] * msk1_1673[k];
    }

#pragma omp simd aligned(t_1674, t_1675, t_1676, pa_x, pc_x, pc_y, msk0_1674, msk0_1676, \
                         msi_1050, msi_1306, msi_1308, msk1_1674, msk1_1676, \
                         nsi_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1674[k] = pa_x[k] * msk0_1674[k]
                    + f_14 * msi_1306[k]
                    - f_12 * pc_x[k] * msk1_1674[k];

        t_1675[k] = f_21 * msi_1050[k]
                    + f_3 * pc_y[k] * nsi_1302[k];

        t_1676[k] = pa_x[k] * msk0_1676[k]
                    + f_14 * msi_1308[k]
                    - f_12 * pc_x[k] * msk1_1676[k];
    }

#pragma omp simd aligned(t_1677, t_1678, t_1679, t_1680, t_1681, pc_x, msi_1309, msi_1310, \
                         msi_1311, msi_1312, msi_1313, nsi_1309, nsi_1310, nsi_1311, nsi_1312, \
                         nsi_1313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1677[k] = f_13 * msi_1309[k]
                    + f_3 * pc_x[k] * nsi_1309[k];

        t_1678[k] = f_13 * msi_1310[k]
                    + f_3 * pc_x[k] * nsi_1310[k];

        t_1679[k] = f_13 * msi_1311[k]
                    + f_3 * pc_x[k] * nsi_1311[k];

        t_1680[k] = f_13 * msi_1312[k]
                    + f_3 * pc_x[k] * nsi_1312[k];

        t_1681[k] = f_13 * msi_1313[k]
                    + f_3 * pc_x[k] * nsi_1313[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, t_1685, pa_x, pc_x, pc_z, msk0_1684, \
                         msi_1029, msi_1314, msi_1315, msk1_1684, nsi_1309, nsi_1314, \
                         nsi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_13 * msi_1314[k]
                    + f_3 * pc_x[k] * nsi_1314[k];

        t_1683[k] = f_13 * msi_1315[k]
                    + f_3 * pc_x[k] * nsi_1315[k];

        t_1684[k] = pa_x[k] * msk0_1684[k]
                    - f_12 * pc_x[k] * msk1_1684[k];

        t_1685[k] = f_13 * msi_1029[k]
                    + f_3 * pc_z[k] * nsi_1309[k];
    }

#pragma omp simd aligned(t_1686, t_1687, t_1688, t_1689, pa_x, pc_x, msk0_1686, msk0_1687, \
                         msk0_1688, msk0_1689, msk1_1686, msk1_1687, msk1_1688, \
                         msk1_1689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1686[k] = pa_x[k] * msk0_1686[k]
                    - f_12 * pc_x[k] * msk1_1686[k];

        t_1687[k] = pa_x[k] * msk0_1687[k]
                    - f_12 * pc_x[k] * msk1_1687[k];

        t_1688[k] = pa_x[k] * msk0_1688[k]
                    - f_12 * pc_x[k] * msk1_1688[k];

        t_1689[k] = pa_x[k] * msk0_1689[k]
                    - f_12 * pc_x[k] * msk1_1689[k];
    }

#pragma omp simd aligned(t_1690, t_1691, t_1692, t_1693, pa_x, pc_x, pc_y, msk0_1691, \
                         msk0_1692, msi_1063, msi_1064, msi_1316, msk1_1691, msk1_1692, \
                         nsi_1315, nsi_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1690[k] = f_21 * msi_1063[k]
                    + f_3 * pc_y[k] * nsi_1315[k];

        t_1691[k] = pa_x[k] * msk0_1691[k]
                    - f_12 * pc_x[k] * msk1_1691[k];

        t_1692[k] = pa_x[k] * msk0_1692[k]
                    + f_22 * msi_1316[k]
                    - f_12 * pc_x[k] * msk1_1692[k];

        t_1693[k] = f_22 * msi_1064[k]
                    + f_3 * pc_y[k] * nsi_1316[k];
    }

#pragma omp simd aligned(t_1694, t_1695, t_1696, pa_x, pc_x, pc_y, pc_z, msk0_1695, msi_1036, \
                         msi_1066, msi_1319, msk1_1695, nsi_1316, \
                         nsi_1318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1694[k] = f_14 * msi_1036[k]
                    + f_3 * pc_z[k] * nsi_1316[k];

        t_1695[k] = pa_x[k] * msk0_1695[k]
                    + f_17 * msi_1319[k]
                    - f_12 * pc_x[k] * msk1_1695[k];

        t_1696[k] = f_22 * msi_1066[k]
                    + f_3 * pc_y[k] * nsi_1318[k];
    }

#pragma omp simd aligned(t_1697, t_1698, t_1699, pa_x, pc_x, pc_z, msk0_1697, msk0_1698, \
                         msi_1039, msi_1321, msi_1322, msk1_1697, msk1_1698, \
                         nsi_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1697[k] = pa_x[k] * msk0_1697[k]
                    + f_17 * msi_1321[k]
                    - f_12 * pc_x[k] * msk1_1697[k];

        t_1698[k] = pa_x[k] * msk0_1698[k]
                    + f_16 * msi_1322[k]
                    - f_12 * pc_x[k] * msk1_1698[k];

        t_1699[k] = f_14 * msi_1039[k]
                    + f_3 * pc_z[k] * nsi_1319[k];
    }

#pragma omp simd aligned(t_1700, t_1701, t_1702, pa_x, pc_x, pc_y, msk0_1701, msk0_1702, \
                         msi_1069, msi_1325, msi_1326, msk1_1701, msk1_1702, \
                         nsi_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1700[k] = f_22 * msi_1069[k]
                    + f_3 * pc_y[k] * nsi_1321[k];

        t_1701[k] = pa_x[k] * msk0_1701[k]
                    + f_16 * msi_1325[k]
                    - f_12 * pc_x[k] * msk1_1701[k];

        t_1702[k] = pa_x[k] * msk0_1702[k]
                    + f_15 * msi_1326[k]
                    - f_12 * pc_x[k] * msk1_1702[k];
    }

#pragma omp simd aligned(t_1703, t_1704, t_1705, pa_x, pc_x, pc_y, pc_z, msk0_1704, msi_1042, \
                         msi_1073, msi_1328, msk1_1704, nsi_1322, \
                         nsi_1325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1703[k] = f_14 * msi_1042[k]
                    + f_3 * pc_z[k] * nsi_1322[k];

        t_1704[k] = pa_x[k] * msk0_1704[k]
                    + f_15 * msi_1328[k]
                    - f_12 * pc_x[k] * msk1_1704[k];

        t_1705[k] = f_22 * msi_1073[k]
                    + f_3 * pc_y[k] * nsi_1325[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, pa_x, pc_x, pc_z, msk0_1706, msk0_1707, \
                         msi_1046, msi_1330, msi_1331, msk1_1706, msk1_1707, \
                         nsi_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = pa_x[k] * msk0_1706[k]
                    + f_15 * msi_1330[k]
                    - f_12 * pc_x[k] * msk1_1706[k];

        t_1707[k] = pa_x[k] * msk0_1707[k]
                    + f_14 * msi_1331[k]
                    - f_12 * pc_x[k] * msk1_1707[k];

        t_1708[k] = f_14 * msi_1046[k]
                    + f_3 * pc_z[k] * nsi_1326[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsi, const size_t ncols,
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
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1709 = buffer.data(msk0 + 1709);
    const auto *msk0_1710 = buffer.data(msk0 + 1710);
    const auto *msk0_1712 = buffer.data(msk0 + 1712);
    const auto *msk0_1720 = buffer.data(msk0 + 1720);
    const auto *msk0_1722 = buffer.data(msk0 + 1722);
    const auto *msk0_1723 = buffer.data(msk0 + 1723);
    const auto *msk0_1724 = buffer.data(msk0 + 1724);
    const auto *msk0_1725 = buffer.data(msk0 + 1725);
    const auto *msk0_1727 = buffer.data(msk0 + 1727);
    const auto *msk0_1728 = buffer.data(msk0 + 1728);
    const auto *msk0_1731 = buffer.data(msk0 + 1731);
    const auto *msk0_1733 = buffer.data(msk0 + 1733);
    const auto *msk0_1734 = buffer.data(msk0 + 1734);
    const auto *msk0_1737 = buffer.data(msk0 + 1737);
    const auto *msk0_1738 = buffer.data(msk0 + 1738);
    const auto *msk0_1740 = buffer.data(msk0 + 1740);
    const auto *msk0_1742 = buffer.data(msk0 + 1742);
    const auto *msk0_1743 = buffer.data(msk0 + 1743);
    const auto *msk0_1745 = buffer.data(msk0 + 1745);
    const auto *msk0_1746 = buffer.data(msk0 + 1746);
    const auto *msk0_1748 = buffer.data(msk0 + 1748);
    const auto *msk0_1756 = buffer.data(msk0 + 1756);
    const auto *msk0_1758 = buffer.data(msk0 + 1758);
    const auto *msk0_1759 = buffer.data(msk0 + 1759);
    const auto *msk0_1760 = buffer.data(msk0 + 1760);
    const auto *msk0_1761 = buffer.data(msk0 + 1761);
    const auto *msk0_1763 = buffer.data(msk0 + 1763);
    const auto *msk0_1764 = buffer.data(msk0 + 1764);
    const auto *msk0_1767 = buffer.data(msk0 + 1767);
    const auto *msk0_1769 = buffer.data(msk0 + 1769);
    const auto *msk0_1770 = buffer.data(msk0 + 1770);
    const auto *msk0_1773 = buffer.data(msk0 + 1773);
    const auto *msk0_1774 = buffer.data(msk0 + 1774);
    const auto *msk0_1776 = buffer.data(msk0 + 1776);
    const auto *msk0_1778 = buffer.data(msk0 + 1778);
    const auto *msk0_1779 = buffer.data(msk0 + 1779);
    const auto *msk0_1781 = buffer.data(msk0 + 1781);
    const auto *msk0_1782 = buffer.data(msk0 + 1782);
    const auto *msk0_1784 = buffer.data(msk0 + 1784);
    const auto *msk0_1792 = buffer.data(msk0 + 1792);
    const auto *msk0_1794 = buffer.data(msk0 + 1794);
    const auto *msk0_1795 = buffer.data(msk0 + 1795);
    const auto *msk0_1796 = buffer.data(msk0 + 1796);
    const auto *msk0_1797 = buffer.data(msk0 + 1797);
    const auto *msk0_1799 = buffer.data(msk0 + 1799);
    const auto *msk0_1800 = buffer.data(msk0 + 1800);
    const auto *msk0_1803 = buffer.data(msk0 + 1803);
    const auto *msk0_1805 = buffer.data(msk0 + 1805);
    const auto *msk0_1806 = buffer.data(msk0 + 1806);
    const auto *msk0_1809 = buffer.data(msk0 + 1809);
    const auto *msk0_1810 = buffer.data(msk0 + 1810);
    const auto *msk0_1812 = buffer.data(msk0 + 1812);
    const auto *msk0_1814 = buffer.data(msk0 + 1814);
    const auto *msk0_1815 = buffer.data(msk0 + 1815);
    const auto *msk0_1817 = buffer.data(msk0 + 1817);
    const auto *msk0_1818 = buffer.data(msk0 + 1818);
    const auto *msk0_1820 = buffer.data(msk0 + 1820);

    const auto *msi_1057 = buffer.data(msi + 1057);
    const auto *msi_1064 = buffer.data(msi + 1064);
    const auto *msi_1067 = buffer.data(msi + 1067);
    const auto *msi_1070 = buffer.data(msi + 1070);
    const auto *msi_1074 = buffer.data(msi + 1074);
    const auto *msi_1078 = buffer.data(msi + 1078);
    const auto *msi_1085 = buffer.data(msi + 1085);
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
    const auto *msi_1119 = buffer.data(msi + 1119);
    const auto *msi_1120 = buffer.data(msi + 1120);
    const auto *msi_1122 = buffer.data(msi + 1122);
    const auto *msi_1123 = buffer.data(msi + 1123);
    const auto *msi_1125 = buffer.data(msi + 1125);
    const auto *msi_1126 = buffer.data(msi + 1126);
    const auto *msi_1129 = buffer.data(msi + 1129);
    const auto *msi_1130 = buffer.data(msi + 1130);
    const auto *msi_1134 = buffer.data(msi + 1134);
    const auto *msi_1147 = buffer.data(msi + 1147);
    const auto *msi_1148 = buffer.data(msi + 1148);
    const auto *msi_1150 = buffer.data(msi + 1150);
    const auto *msi_1153 = buffer.data(msi + 1153);
    const auto *msi_1157 = buffer.data(msi + 1157);
    const auto *msi_1162 = buffer.data(msi + 1162);
    const auto *msi_1333 = buffer.data(msi + 1333);
    const auto *msi_1334 = buffer.data(msi + 1334);
    const auto *msi_1336 = buffer.data(msi + 1336);
    const auto *msi_1337 = buffer.data(msi + 1337);
    const auto *msi_1338 = buffer.data(msi + 1338);
    const auto *msi_1339 = buffer.data(msi + 1339);
    const auto *msi_1340 = buffer.data(msi + 1340);
    const auto *msi_1341 = buffer.data(msi + 1341);
    const auto *msi_1342 = buffer.data(msi + 1342);
    const auto *msi_1343 = buffer.data(msi + 1343);
    const auto *msi_1344 = buffer.data(msi + 1344);
    const auto *msi_1347 = buffer.data(msi + 1347);
    const auto *msi_1349 = buffer.data(msi + 1349);
    const auto *msi_1350 = buffer.data(msi + 1350);
    const auto *msi_1353 = buffer.data(msi + 1353);
    const auto *msi_1354 = buffer.data(msi + 1354);
    const auto *msi_1356 = buffer.data(msi + 1356);
    const auto *msi_1358 = buffer.data(msi + 1358);
    const auto *msi_1359 = buffer.data(msi + 1359);
    const auto *msi_1361 = buffer.data(msi + 1361);
    const auto *msi_1362 = buffer.data(msi + 1362);
    const auto *msi_1364 = buffer.data(msi + 1364);
    const auto *msi_1365 = buffer.data(msi + 1365);
    const auto *msi_1366 = buffer.data(msi + 1366);
    const auto *msi_1367 = buffer.data(msi + 1367);
    const auto *msi_1368 = buffer.data(msi + 1368);
    const auto *msi_1369 = buffer.data(msi + 1369);
    const auto *msi_1370 = buffer.data(msi + 1370);
    const auto *msi_1371 = buffer.data(msi + 1371);
    const auto *msi_1372 = buffer.data(msi + 1372);
    const auto *msi_1375 = buffer.data(msi + 1375);
    const auto *msi_1377 = buffer.data(msi + 1377);
    const auto *msi_1378 = buffer.data(msi + 1378);
    const auto *msi_1381 = buffer.data(msi + 1381);
    const auto *msi_1382 = buffer.data(msi + 1382);
    const auto *msi_1384 = buffer.data(msi + 1384);
    const auto *msi_1386 = buffer.data(msi + 1386);
    const auto *msi_1387 = buffer.data(msi + 1387);
    const auto *msi_1389 = buffer.data(msi + 1389);
    const auto *msi_1390 = buffer.data(msi + 1390);
    const auto *msi_1392 = buffer.data(msi + 1392);
    const auto *msi_1393 = buffer.data(msi + 1393);
    const auto *msi_1394 = buffer.data(msi + 1394);
    const auto *msi_1395 = buffer.data(msi + 1395);
    const auto *msi_1396 = buffer.data(msi + 1396);
    const auto *msi_1397 = buffer.data(msi + 1397);
    const auto *msi_1398 = buffer.data(msi + 1398);
    const auto *msi_1399 = buffer.data(msi + 1399);
    const auto *msi_1400 = buffer.data(msi + 1400);
    const auto *msi_1403 = buffer.data(msi + 1403);
    const auto *msi_1405 = buffer.data(msi + 1405);
    const auto *msi_1406 = buffer.data(msi + 1406);
    const auto *msi_1409 = buffer.data(msi + 1409);
    const auto *msi_1410 = buffer.data(msi + 1410);
    const auto *msi_1412 = buffer.data(msi + 1412);
    const auto *msi_1414 = buffer.data(msi + 1414);
    const auto *msi_1415 = buffer.data(msi + 1415);
    const auto *msi_1417 = buffer.data(msi + 1417);
    const auto *msi_1418 = buffer.data(msi + 1418);
    const auto *msi_1420 = buffer.data(msi + 1420);
    const auto *msi_1421 = buffer.data(msi + 1421);
    const auto *msi_1422 = buffer.data(msi + 1422);
    const auto *msi_1423 = buffer.data(msi + 1423);
    const auto *msi_1424 = buffer.data(msi + 1424);
    const auto *msi_1425 = buffer.data(msi + 1425);

    const auto *msk1_1709 = buffer.data(msk1 + 1709);
    const auto *msk1_1710 = buffer.data(msk1 + 1710);
    const auto *msk1_1712 = buffer.data(msk1 + 1712);
    const auto *msk1_1720 = buffer.data(msk1 + 1720);
    const auto *msk1_1722 = buffer.data(msk1 + 1722);
    const auto *msk1_1723 = buffer.data(msk1 + 1723);
    const auto *msk1_1724 = buffer.data(msk1 + 1724);
    const auto *msk1_1725 = buffer.data(msk1 + 1725);
    const auto *msk1_1727 = buffer.data(msk1 + 1727);
    const auto *msk1_1728 = buffer.data(msk1 + 1728);
    const auto *msk1_1731 = buffer.data(msk1 + 1731);
    const auto *msk1_1733 = buffer.data(msk1 + 1733);
    const auto *msk1_1734 = buffer.data(msk1 + 1734);
    const auto *msk1_1737 = buffer.data(msk1 + 1737);
    const auto *msk1_1738 = buffer.data(msk1 + 1738);
    const auto *msk1_1740 = buffer.data(msk1 + 1740);
    const auto *msk1_1742 = buffer.data(msk1 + 1742);
    const auto *msk1_1743 = buffer.data(msk1 + 1743);
    const auto *msk1_1745 = buffer.data(msk1 + 1745);
    const auto *msk1_1746 = buffer.data(msk1 + 1746);
    const auto *msk1_1748 = buffer.data(msk1 + 1748);
    const auto *msk1_1756 = buffer.data(msk1 + 1756);
    const auto *msk1_1758 = buffer.data(msk1 + 1758);
    const auto *msk1_1759 = buffer.data(msk1 + 1759);
    const auto *msk1_1760 = buffer.data(msk1 + 1760);
    const auto *msk1_1761 = buffer.data(msk1 + 1761);
    const auto *msk1_1763 = buffer.data(msk1 + 1763);
    const auto *msk1_1764 = buffer.data(msk1 + 1764);
    const auto *msk1_1767 = buffer.data(msk1 + 1767);
    const auto *msk1_1769 = buffer.data(msk1 + 1769);
    const auto *msk1_1770 = buffer.data(msk1 + 1770);
    const auto *msk1_1773 = buffer.data(msk1 + 1773);
    const auto *msk1_1774 = buffer.data(msk1 + 1774);
    const auto *msk1_1776 = buffer.data(msk1 + 1776);
    const auto *msk1_1778 = buffer.data(msk1 + 1778);
    const auto *msk1_1779 = buffer.data(msk1 + 1779);
    const auto *msk1_1781 = buffer.data(msk1 + 1781);
    const auto *msk1_1782 = buffer.data(msk1 + 1782);
    const auto *msk1_1784 = buffer.data(msk1 + 1784);
    const auto *msk1_1792 = buffer.data(msk1 + 1792);
    const auto *msk1_1794 = buffer.data(msk1 + 1794);
    const auto *msk1_1795 = buffer.data(msk1 + 1795);
    const auto *msk1_1796 = buffer.data(msk1 + 1796);
    const auto *msk1_1797 = buffer.data(msk1 + 1797);
    const auto *msk1_1799 = buffer.data(msk1 + 1799);
    const auto *msk1_1800 = buffer.data(msk1 + 1800);
    const auto *msk1_1803 = buffer.data(msk1 + 1803);
    const auto *msk1_1805 = buffer.data(msk1 + 1805);
    const auto *msk1_1806 = buffer.data(msk1 + 1806);
    const auto *msk1_1809 = buffer.data(msk1 + 1809);
    const auto *msk1_1810 = buffer.data(msk1 + 1810);
    const auto *msk1_1812 = buffer.data(msk1 + 1812);
    const auto *msk1_1814 = buffer.data(msk1 + 1814);
    const auto *msk1_1815 = buffer.data(msk1 + 1815);
    const auto *msk1_1817 = buffer.data(msk1 + 1817);
    const auto *msk1_1818 = buffer.data(msk1 + 1818);
    const auto *msk1_1820 = buffer.data(msk1 + 1820);

    const auto *nsi_1330 = buffer.data(nsi + 1330);
    const auto *nsi_1337 = buffer.data(nsi + 1337);
    const auto *nsi_1338 = buffer.data(nsi + 1338);
    const auto *nsi_1339 = buffer.data(nsi + 1339);
    const auto *nsi_1340 = buffer.data(nsi + 1340);
    const auto *nsi_1341 = buffer.data(nsi + 1341);
    const auto *nsi_1342 = buffer.data(nsi + 1342);
    const auto *nsi_1343 = buffer.data(nsi + 1343);
    const auto *nsi_1344 = buffer.data(nsi + 1344);
    const auto *nsi_1346 = buffer.data(nsi + 1346);
    const auto *nsi_1347 = buffer.data(nsi + 1347);
    const auto *nsi_1349 = buffer.data(nsi + 1349);
    const auto *nsi_1350 = buffer.data(nsi + 1350);
    const auto *nsi_1353 = buffer.data(nsi + 1353);
    const auto *nsi_1354 = buffer.data(nsi + 1354);
    const auto *nsi_1358 = buffer.data(nsi + 1358);
    const auto *nsi_1365 = buffer.data(nsi + 1365);
    const auto *nsi_1366 = buffer.data(nsi + 1366);
    const auto *nsi_1367 = buffer.data(nsi + 1367);
    const auto *nsi_1368 = buffer.data(nsi + 1368);
    const auto *nsi_1369 = buffer.data(nsi + 1369);
    const auto *nsi_1370 = buffer.data(nsi + 1370);
    const auto *nsi_1371 = buffer.data(nsi + 1371);
    const auto *nsi_1372 = buffer.data(nsi + 1372);
    const auto *nsi_1374 = buffer.data(nsi + 1374);
    const auto *nsi_1375 = buffer.data(nsi + 1375);
    const auto *nsi_1377 = buffer.data(nsi + 1377);
    const auto *nsi_1378 = buffer.data(nsi + 1378);
    const auto *nsi_1381 = buffer.data(nsi + 1381);
    const auto *nsi_1382 = buffer.data(nsi + 1382);
    const auto *nsi_1386 = buffer.data(nsi + 1386);
    const auto *nsi_1393 = buffer.data(nsi + 1393);
    const auto *nsi_1394 = buffer.data(nsi + 1394);
    const auto *nsi_1395 = buffer.data(nsi + 1395);
    const auto *nsi_1396 = buffer.data(nsi + 1396);
    const auto *nsi_1397 = buffer.data(nsi + 1397);
    const auto *nsi_1398 = buffer.data(nsi + 1398);
    const auto *nsi_1399 = buffer.data(nsi + 1399);
    const auto *nsi_1400 = buffer.data(nsi + 1400);
    const auto *nsi_1402 = buffer.data(nsi + 1402);
    const auto *nsi_1403 = buffer.data(nsi + 1403);
    const auto *nsi_1405 = buffer.data(nsi + 1405);
    const auto *nsi_1406 = buffer.data(nsi + 1406);
    const auto *nsi_1409 = buffer.data(nsi + 1409);
    const auto *nsi_1410 = buffer.data(nsi + 1410);
    const auto *nsi_1414 = buffer.data(nsi + 1414);
    const auto *nsi_1421 = buffer.data(nsi + 1421);
    const auto *nsi_1422 = buffer.data(nsi + 1422);
    const auto *nsi_1423 = buffer.data(nsi + 1423);
    const auto *nsi_1424 = buffer.data(nsi + 1424);
    const auto *nsi_1425 = buffer.data(nsi + 1425);

#pragma omp simd aligned(t_1709, t_1710, t_1711, pa_x, pc_x, pc_y, msk0_1709, msk0_1710, \
                         msi_1078, msi_1333, msi_1334, msk1_1709, msk1_1710, \
                         nsi_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1709[k] = pa_x[k] * msk0_1709[k]
                    + f_14 * msi_1333[k]
                    - f_12 * pc_x[k] * msk1_1709[k];

        t_1710[k] = pa_x[k] * msk0_1710[k]
                    + f_14 * msi_1334[k]
                    - f_12 * pc_x[k] * msk1_1710[k];

        t_1711[k] = f_22 * msi_1078[k]
                    + f_3 * pc_y[k] * nsi_1330[k];
    }

#pragma omp simd aligned(t_1712, t_1713, t_1714, t_1715, pa_x, pc_x, msk0_1712, msi_1336, \
                         msi_1337, msi_1338, msi_1339, msk1_1712, nsi_1337, nsi_1338, \
                         nsi_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1712[k] = pa_x[k] * msk0_1712[k]
                    + f_14 * msi_1336[k]
                    - f_12 * pc_x[k] * msk1_1712[k];

        t_1713[k] = f_13 * msi_1337[k]
                    + f_3 * pc_x[k] * nsi_1337[k];

        t_1714[k] = f_13 * msi_1338[k]
                    + f_3 * pc_x[k] * nsi_1338[k];

        t_1715[k] = f_13 * msi_1339[k]
                    + f_3 * pc_x[k] * nsi_1339[k];
    }

#pragma omp simd aligned(t_1716, t_1717, t_1718, t_1719, pc_x, msi_1340, msi_1341, msi_1342, \
                         msi_1343, nsi_1340, nsi_1341, nsi_1342, \
                         nsi_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1716[k] = f_13 * msi_1340[k]
                    + f_3 * pc_x[k] * nsi_1340[k];

        t_1717[k] = f_13 * msi_1341[k]
                    + f_3 * pc_x[k] * nsi_1341[k];

        t_1718[k] = f_13 * msi_1342[k]
                    + f_3 * pc_x[k] * nsi_1342[k];

        t_1719[k] = f_13 * msi_1343[k]
                    + f_3 * pc_x[k] * nsi_1343[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, t_1723, pa_x, pc_x, pc_z, msk0_1720, \
                         msk0_1722, msk0_1723, msi_1057, msk1_1720, msk1_1722, msk1_1723, \
                         nsi_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = pa_x[k] * msk0_1720[k]
                    - f_12 * pc_x[k] * msk1_1720[k];

        t_1721[k] = f_14 * msi_1057[k]
                    + f_3 * pc_z[k] * nsi_1337[k];

        t_1722[k] = pa_x[k] * msk0_1722[k]
                    - f_12 * pc_x[k] * msk1_1722[k];

        t_1723[k] = pa_x[k] * msk0_1723[k]
                    - f_12 * pc_x[k] * msk1_1723[k];
    }

#pragma omp simd aligned(t_1724, t_1725, t_1726, t_1727, pa_x, pc_x, pc_y, msk0_1724, \
                         msk0_1725, msk0_1727, msi_1091, msk1_1724, msk1_1725, msk1_1727, \
                         nsi_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1724[k] = pa_x[k] * msk0_1724[k]
                    - f_12 * pc_x[k] * msk1_1724[k];

        t_1725[k] = pa_x[k] * msk0_1725[k]
                    - f_12 * pc_x[k] * msk1_1725[k];

        t_1726[k] = f_22 * msi_1091[k]
                    + f_3 * pc_y[k] * nsi_1343[k];

        t_1727[k] = pa_x[k] * msk0_1727[k]
                    - f_12 * pc_x[k] * msk1_1727[k];
    }

#pragma omp simd aligned(t_1728, t_1729, t_1730, pa_x, pc_x, pc_y, pc_z, msk0_1728, msi_1064, \
                         msi_1092, msi_1344, msk1_1728, nsi_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1728[k] = pa_x[k] * msk0_1728[k]
                    + f_22 * msi_1344[k]
                    - f_12 * pc_x[k] * msk1_1728[k];

        t_1729[k] = f_23 * msi_1092[k]
                    + f_3 * pc_y[k] * nsi_1344[k];

        t_1730[k] = f_15 * msi_1064[k]
                    + f_3 * pc_z[k] * nsi_1344[k];
    }

#pragma omp simd aligned(t_1731, t_1732, t_1733, pa_x, pc_x, pc_y, msk0_1731, msk0_1733, \
                         msi_1094, msi_1347, msi_1349, msk1_1731, msk1_1733, \
                         nsi_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1731[k] = pa_x[k] * msk0_1731[k]
                    + f_17 * msi_1347[k]
                    - f_12 * pc_x[k] * msk1_1731[k];

        t_1732[k] = f_23 * msi_1094[k]
                    + f_3 * pc_y[k] * nsi_1346[k];

        t_1733[k] = pa_x[k] * msk0_1733[k]
                    + f_17 * msi_1349[k]
                    - f_12 * pc_x[k] * msk1_1733[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, pa_x, pc_x, pc_y, pc_z, msk0_1734, msi_1067, \
                         msi_1097, msi_1350, msk1_1734, nsi_1347, \
                         nsi_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = pa_x[k] * msk0_1734[k]
                    + f_16 * msi_1350[k]
                    - f_12 * pc_x[k] * msk1_1734[k];

        t_1735[k] = f_15 * msi_1067[k]
                    + f_3 * pc_z[k] * nsi_1347[k];

        t_1736[k] = f_23 * msi_1097[k]
                    + f_3 * pc_y[k] * nsi_1349[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, pa_x, pc_x, pc_z, msk0_1737, msk0_1738, \
                         msi_1070, msi_1353, msi_1354, msk1_1737, msk1_1738, \
                         nsi_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = pa_x[k] * msk0_1737[k]
                    + f_16 * msi_1353[k]
                    - f_12 * pc_x[k] * msk1_1737[k];

        t_1738[k] = pa_x[k] * msk0_1738[k]
                    + f_15 * msi_1354[k]
                    - f_12 * pc_x[k] * msk1_1738[k];

        t_1739[k] = f_15 * msi_1070[k]
                    + f_3 * pc_z[k] * nsi_1350[k];
    }

#pragma omp simd aligned(t_1740, t_1741, t_1742, pa_x, pc_x, pc_y, msk0_1740, msk0_1742, \
                         msi_1101, msi_1356, msi_1358, msk1_1740, msk1_1742, \
                         nsi_1353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1740[k] = pa_x[k] * msk0_1740[k]
                    + f_15 * msi_1356[k]
                    - f_12 * pc_x[k] * msk1_1740[k];

        t_1741[k] = f_23 * msi_1101[k]
                    + f_3 * pc_y[k] * nsi_1353[k];

        t_1742[k] = pa_x[k] * msk0_1742[k]
                    + f_15 * msi_1358[k]
                    - f_12 * pc_x[k] * msk1_1742[k];
    }

#pragma omp simd aligned(t_1743, t_1744, t_1745, pa_x, pc_x, pc_z, msk0_1743, msk0_1745, \
                         msi_1074, msi_1359, msi_1361, msk1_1743, msk1_1745, \
                         nsi_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1743[k] = pa_x[k] * msk0_1743[k]
                    + f_14 * msi_1359[k]
                    - f_12 * pc_x[k] * msk1_1743[k];

        t_1744[k] = f_15 * msi_1074[k]
                    + f_3 * pc_z[k] * nsi_1354[k];

        t_1745[k] = pa_x[k] * msk0_1745[k]
                    + f_14 * msi_1361[k]
                    - f_12 * pc_x[k] * msk1_1745[k];
    }

#pragma omp simd aligned(t_1746, t_1747, t_1748, pa_x, pc_x, pc_y, msk0_1746, msk0_1748, \
                         msi_1106, msi_1362, msi_1364, msk1_1746, msk1_1748, \
                         nsi_1358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1746[k] = pa_x[k] * msk0_1746[k]
                    + f_14 * msi_1362[k]
                    - f_12 * pc_x[k] * msk1_1746[k];

        t_1747[k] = f_23 * msi_1106[k]
                    + f_3 * pc_y[k] * nsi_1358[k];

        t_1748[k] = pa_x[k] * msk0_1748[k]
                    + f_14 * msi_1364[k]
                    - f_12 * pc_x[k] * msk1_1748[k];
    }

#pragma omp simd aligned(t_1749, t_1750, t_1751, t_1752, t_1753, pc_x, msi_1365, msi_1366, \
                         msi_1367, msi_1368, msi_1369, nsi_1365, nsi_1366, nsi_1367, nsi_1368, \
                         nsi_1369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1749[k] = f_13 * msi_1365[k]
                    + f_3 * pc_x[k] * nsi_1365[k];

        t_1750[k] = f_13 * msi_1366[k]
                    + f_3 * pc_x[k] * nsi_1366[k];

        t_1751[k] = f_13 * msi_1367[k]
                    + f_3 * pc_x[k] * nsi_1367[k];

        t_1752[k] = f_13 * msi_1368[k]
                    + f_3 * pc_x[k] * nsi_1368[k];

        t_1753[k] = f_13 * msi_1369[k]
                    + f_3 * pc_x[k] * nsi_1369[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, pa_x, pc_x, pc_z, msk0_1756, \
                         msi_1085, msi_1370, msi_1371, msk1_1756, nsi_1365, nsi_1370, \
                         nsi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = f_13 * msi_1370[k]
                    + f_3 * pc_x[k] * nsi_1370[k];

        t_1755[k] = f_13 * msi_1371[k]
                    + f_3 * pc_x[k] * nsi_1371[k];

        t_1756[k] = pa_x[k] * msk0_1756[k]
                    - f_12 * pc_x[k] * msk1_1756[k];

        t_1757[k] = f_15 * msi_1085[k]
                    + f_3 * pc_z[k] * nsi_1365[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, t_1761, pa_x, pc_x, msk0_1758, msk0_1759, \
                         msk0_1760, msk0_1761, msk1_1758, msk1_1759, msk1_1760, \
                         msk1_1761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = pa_x[k] * msk0_1758[k]
                    - f_12 * pc_x[k] * msk1_1758[k];

        t_1759[k] = pa_x[k] * msk0_1759[k]
                    - f_12 * pc_x[k] * msk1_1759[k];

        t_1760[k] = pa_x[k] * msk0_1760[k]
                    - f_12 * pc_x[k] * msk1_1760[k];

        t_1761[k] = pa_x[k] * msk0_1761[k]
                    - f_12 * pc_x[k] * msk1_1761[k];
    }

#pragma omp simd aligned(t_1762, t_1763, t_1764, t_1765, pa_x, pc_x, pc_y, msk0_1763, \
                         msk0_1764, msi_1119, msi_1120, msi_1372, msk1_1763, msk1_1764, \
                         nsi_1371, nsi_1372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1762[k] = f_23 * msi_1119[k]
                    + f_3 * pc_y[k] * nsi_1371[k];

        t_1763[k] = pa_x[k] * msk0_1763[k]
                    - f_12 * pc_x[k] * msk1_1763[k];

        t_1764[k] = pa_x[k] * msk0_1764[k]
                    + f_22 * msi_1372[k]
                    - f_12 * pc_x[k] * msk1_1764[k];

        t_1765[k] = f_17 * msi_1120[k]
                    + f_3 * pc_y[k] * nsi_1372[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pa_x, pc_x, pc_y, pc_z, msk0_1767, msi_1092, \
                         msi_1122, msi_1375, msk1_1767, nsi_1372, \
                         nsi_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_16 * msi_1092[k]
                    + f_3 * pc_z[k] * nsi_1372[k];

        t_1767[k] = pa_x[k] * msk0_1767[k]
                    + f_17 * msi_1375[k]
                    - f_12 * pc_x[k] * msk1_1767[k];

        t_1768[k] = f_17 * msi_1122[k]
                    + f_3 * pc_y[k] * nsi_1374[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, pa_x, pc_x, pc_z, msk0_1769, msk0_1770, \
                         msi_1095, msi_1377, msi_1378, msk1_1769, msk1_1770, \
                         nsi_1375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = pa_x[k] * msk0_1769[k]
                    + f_17 * msi_1377[k]
                    - f_12 * pc_x[k] * msk1_1769[k];

        t_1770[k] = pa_x[k] * msk0_1770[k]
                    + f_16 * msi_1378[k]
                    - f_12 * pc_x[k] * msk1_1770[k];

        t_1771[k] = f_16 * msi_1095[k]
                    + f_3 * pc_z[k] * nsi_1375[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, pa_x, pc_x, pc_y, msk0_1773, msk0_1774, \
                         msi_1125, msi_1381, msi_1382, msk1_1773, msk1_1774, \
                         nsi_1377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = f_17 * msi_1125[k]
                    + f_3 * pc_y[k] * nsi_1377[k];

        t_1773[k] = pa_x[k] * msk0_1773[k]
                    + f_16 * msi_1381[k]
                    - f_12 * pc_x[k] * msk1_1773[k];

        t_1774[k] = pa_x[k] * msk0_1774[k]
                    + f_15 * msi_1382[k]
                    - f_12 * pc_x[k] * msk1_1774[k];
    }

#pragma omp simd aligned(t_1775, t_1776, t_1777, pa_x, pc_x, pc_y, pc_z, msk0_1776, msi_1098, \
                         msi_1129, msi_1384, msk1_1776, nsi_1378, \
                         nsi_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1775[k] = f_16 * msi_1098[k]
                    + f_3 * pc_z[k] * nsi_1378[k];

        t_1776[k] = pa_x[k] * msk0_1776[k]
                    + f_15 * msi_1384[k]
                    - f_12 * pc_x[k] * msk1_1776[k];

        t_1777[k] = f_17 * msi_1129[k]
                    + f_3 * pc_y[k] * nsi_1381[k];
    }

#pragma omp simd aligned(t_1778, t_1779, t_1780, pa_x, pc_x, pc_z, msk0_1778, msk0_1779, \
                         msi_1102, msi_1386, msi_1387, msk1_1778, msk1_1779, \
                         nsi_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1778[k] = pa_x[k] * msk0_1778[k]
                    + f_15 * msi_1386[k]
                    - f_12 * pc_x[k] * msk1_1778[k];

        t_1779[k] = pa_x[k] * msk0_1779[k]
                    + f_14 * msi_1387[k]
                    - f_12 * pc_x[k] * msk1_1779[k];

        t_1780[k] = f_16 * msi_1102[k]
                    + f_3 * pc_z[k] * nsi_1382[k];
    }

#pragma omp simd aligned(t_1781, t_1782, t_1783, pa_x, pc_x, pc_y, msk0_1781, msk0_1782, \
                         msi_1134, msi_1389, msi_1390, msk1_1781, msk1_1782, \
                         nsi_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1781[k] = pa_x[k] * msk0_1781[k]
                    + f_14 * msi_1389[k]
                    - f_12 * pc_x[k] * msk1_1781[k];

        t_1782[k] = pa_x[k] * msk0_1782[k]
                    + f_14 * msi_1390[k]
                    - f_12 * pc_x[k] * msk1_1782[k];

        t_1783[k] = f_17 * msi_1134[k]
                    + f_3 * pc_y[k] * nsi_1386[k];
    }

#pragma omp simd aligned(t_1784, t_1785, t_1786, t_1787, pa_x, pc_x, msk0_1784, msi_1392, \
                         msi_1393, msi_1394, msi_1395, msk1_1784, nsi_1393, nsi_1394, \
                         nsi_1395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1784[k] = pa_x[k] * msk0_1784[k]
                    + f_14 * msi_1392[k]
                    - f_12 * pc_x[k] * msk1_1784[k];

        t_1785[k] = f_13 * msi_1393[k]
                    + f_3 * pc_x[k] * nsi_1393[k];

        t_1786[k] = f_13 * msi_1394[k]
                    + f_3 * pc_x[k] * nsi_1394[k];

        t_1787[k] = f_13 * msi_1395[k]
                    + f_3 * pc_x[k] * nsi_1395[k];
    }

#pragma omp simd aligned(t_1788, t_1789, t_1790, t_1791, pc_x, msi_1396, msi_1397, msi_1398, \
                         msi_1399, nsi_1396, nsi_1397, nsi_1398, \
                         nsi_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1788[k] = f_13 * msi_1396[k]
                    + f_3 * pc_x[k] * nsi_1396[k];

        t_1789[k] = f_13 * msi_1397[k]
                    + f_3 * pc_x[k] * nsi_1397[k];

        t_1790[k] = f_13 * msi_1398[k]
                    + f_3 * pc_x[k] * nsi_1398[k];

        t_1791[k] = f_13 * msi_1399[k]
                    + f_3 * pc_x[k] * nsi_1399[k];
    }

#pragma omp simd aligned(t_1792, t_1793, t_1794, t_1795, pa_x, pc_x, pc_z, msk0_1792, \
                         msk0_1794, msk0_1795, msi_1113, msk1_1792, msk1_1794, msk1_1795, \
                         nsi_1393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1792[k] = pa_x[k] * msk0_1792[k]
                    - f_12 * pc_x[k] * msk1_1792[k];

        t_1793[k] = f_16 * msi_1113[k]
                    + f_3 * pc_z[k] * nsi_1393[k];

        t_1794[k] = pa_x[k] * msk0_1794[k]
                    - f_12 * pc_x[k] * msk1_1794[k];

        t_1795[k] = pa_x[k] * msk0_1795[k]
                    - f_12 * pc_x[k] * msk1_1795[k];
    }

#pragma omp simd aligned(t_1796, t_1797, t_1798, t_1799, pa_x, pc_x, pc_y, msk0_1796, \
                         msk0_1797, msk0_1799, msi_1147, msk1_1796, msk1_1797, msk1_1799, \
                         nsi_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1796[k] = pa_x[k] * msk0_1796[k]
                    - f_12 * pc_x[k] * msk1_1796[k];

        t_1797[k] = pa_x[k] * msk0_1797[k]
                    - f_12 * pc_x[k] * msk1_1797[k];

        t_1798[k] = f_17 * msi_1147[k]
                    + f_3 * pc_y[k] * nsi_1399[k];

        t_1799[k] = pa_x[k] * msk0_1799[k]
                    - f_12 * pc_x[k] * msk1_1799[k];
    }

#pragma omp simd aligned(t_1800, t_1801, t_1802, pa_x, pc_x, pc_y, pc_z, msk0_1800, msi_1120, \
                         msi_1148, msi_1400, msk1_1800, nsi_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1800[k] = pa_x[k] * msk0_1800[k]
                    + f_22 * msi_1400[k]
                    - f_12 * pc_x[k] * msk1_1800[k];

        t_1801[k] = f_16 * msi_1148[k]
                    + f_3 * pc_y[k] * nsi_1400[k];

        t_1802[k] = f_17 * msi_1120[k]
                    + f_3 * pc_z[k] * nsi_1400[k];
    }

#pragma omp simd aligned(t_1803, t_1804, t_1805, pa_x, pc_x, pc_y, msk0_1803, msk0_1805, \
                         msi_1150, msi_1403, msi_1405, msk1_1803, msk1_1805, \
                         nsi_1402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1803[k] = pa_x[k] * msk0_1803[k]
                    + f_17 * msi_1403[k]
                    - f_12 * pc_x[k] * msk1_1803[k];

        t_1804[k] = f_16 * msi_1150[k]
                    + f_3 * pc_y[k] * nsi_1402[k];

        t_1805[k] = pa_x[k] * msk0_1805[k]
                    + f_17 * msi_1405[k]
                    - f_12 * pc_x[k] * msk1_1805[k];
    }

#pragma omp simd aligned(t_1806, t_1807, t_1808, pa_x, pc_x, pc_y, pc_z, msk0_1806, msi_1123, \
                         msi_1153, msi_1406, msk1_1806, nsi_1403, \
                         nsi_1405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1806[k] = pa_x[k] * msk0_1806[k]
                    + f_16 * msi_1406[k]
                    - f_12 * pc_x[k] * msk1_1806[k];

        t_1807[k] = f_17 * msi_1123[k]
                    + f_3 * pc_z[k] * nsi_1403[k];

        t_1808[k] = f_16 * msi_1153[k]
                    + f_3 * pc_y[k] * nsi_1405[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, pa_x, pc_x, pc_z, msk0_1809, msk0_1810, \
                         msi_1126, msi_1409, msi_1410, msk1_1809, msk1_1810, \
                         nsi_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = pa_x[k] * msk0_1809[k]
                    + f_16 * msi_1409[k]
                    - f_12 * pc_x[k] * msk1_1809[k];

        t_1810[k] = pa_x[k] * msk0_1810[k]
                    + f_15 * msi_1410[k]
                    - f_12 * pc_x[k] * msk1_1810[k];

        t_1811[k] = f_17 * msi_1126[k]
                    + f_3 * pc_z[k] * nsi_1406[k];
    }

#pragma omp simd aligned(t_1812, t_1813, t_1814, pa_x, pc_x, pc_y, msk0_1812, msk0_1814, \
                         msi_1157, msi_1412, msi_1414, msk1_1812, msk1_1814, \
                         nsi_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1812[k] = pa_x[k] * msk0_1812[k]
                    + f_15 * msi_1412[k]
                    - f_12 * pc_x[k] * msk1_1812[k];

        t_1813[k] = f_16 * msi_1157[k]
                    + f_3 * pc_y[k] * nsi_1409[k];

        t_1814[k] = pa_x[k] * msk0_1814[k]
                    + f_15 * msi_1414[k]
                    - f_12 * pc_x[k] * msk1_1814[k];
    }

#pragma omp simd aligned(t_1815, t_1816, t_1817, pa_x, pc_x, pc_z, msk0_1815, msk0_1817, \
                         msi_1130, msi_1415, msi_1417, msk1_1815, msk1_1817, \
                         nsi_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1815[k] = pa_x[k] * msk0_1815[k]
                    + f_14 * msi_1415[k]
                    - f_12 * pc_x[k] * msk1_1815[k];

        t_1816[k] = f_17 * msi_1130[k]
                    + f_3 * pc_z[k] * nsi_1410[k];

        t_1817[k] = pa_x[k] * msk0_1817[k]
                    + f_14 * msi_1417[k]
                    - f_12 * pc_x[k] * msk1_1817[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, pa_x, pc_x, pc_y, msk0_1818, msk0_1820, \
                         msi_1162, msi_1418, msi_1420, msk1_1818, msk1_1820, \
                         nsi_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = pa_x[k] * msk0_1818[k]
                    + f_14 * msi_1418[k]
                    - f_12 * pc_x[k] * msk1_1818[k];

        t_1819[k] = f_16 * msi_1162[k]
                    + f_3 * pc_y[k] * nsi_1414[k];

        t_1820[k] = pa_x[k] * msk0_1820[k]
                    + f_14 * msi_1420[k]
                    - f_12 * pc_x[k] * msk1_1820[k];
    }

#pragma omp simd aligned(t_1821, t_1822, t_1823, t_1824, t_1825, pc_x, msi_1421, msi_1422, \
                         msi_1423, msi_1424, msi_1425, nsi_1421, nsi_1422, nsi_1423, nsi_1424, \
                         nsi_1425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1821[k] = f_13 * msi_1421[k]
                    + f_3 * pc_x[k] * nsi_1421[k];

        t_1822[k] = f_13 * msi_1422[k]
                    + f_3 * pc_x[k] * nsi_1422[k];

        t_1823[k] = f_13 * msi_1423[k]
                    + f_3 * pc_x[k] * nsi_1423[k];

        t_1824[k] = f_13 * msi_1424[k]
                    + f_3 * pc_x[k] * nsi_1424[k];

        t_1825[k] = f_13 * msi_1425[k]
                    + f_3 * pc_x[k] * nsi_1425[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsi, const size_t ncols,
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
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1584 = buffer.data(msk0 + 1584);
    const auto *msk0_1589 = buffer.data(msk0 + 1589);
    const auto *msk0_1593 = buffer.data(msk0 + 1593);
    const auto *msk0_1598 = buffer.data(msk0 + 1598);
    const auto *msk0_1604 = buffer.data(msk0 + 1604);
    const auto *msk0_1828 = buffer.data(msk0 + 1828);
    const auto *msk0_1830 = buffer.data(msk0 + 1830);
    const auto *msk0_1831 = buffer.data(msk0 + 1831);
    const auto *msk0_1832 = buffer.data(msk0 + 1832);
    const auto *msk0_1833 = buffer.data(msk0 + 1833);
    const auto *msk0_1835 = buffer.data(msk0 + 1835);
    const auto *msk0_1836 = buffer.data(msk0 + 1836);
    const auto *msk0_1839 = buffer.data(msk0 + 1839);
    const auto *msk0_1841 = buffer.data(msk0 + 1841);
    const auto *msk0_1842 = buffer.data(msk0 + 1842);
    const auto *msk0_1845 = buffer.data(msk0 + 1845);
    const auto *msk0_1846 = buffer.data(msk0 + 1846);
    const auto *msk0_1848 = buffer.data(msk0 + 1848);
    const auto *msk0_1850 = buffer.data(msk0 + 1850);
    const auto *msk0_1851 = buffer.data(msk0 + 1851);
    const auto *msk0_1853 = buffer.data(msk0 + 1853);
    const auto *msk0_1854 = buffer.data(msk0 + 1854);
    const auto *msk0_1856 = buffer.data(msk0 + 1856);
    const auto *msk0_1864 = buffer.data(msk0 + 1864);
    const auto *msk0_1866 = buffer.data(msk0 + 1866);
    const auto *msk0_1867 = buffer.data(msk0 + 1867);
    const auto *msk0_1868 = buffer.data(msk0 + 1868);
    const auto *msk0_1869 = buffer.data(msk0 + 1869);
    const auto *msk0_1871 = buffer.data(msk0 + 1871);
    const auto *msk0_1872 = buffer.data(msk0 + 1872);
    const auto *msk0_1875 = buffer.data(msk0 + 1875);
    const auto *msk0_1877 = buffer.data(msk0 + 1877);
    const auto *msk0_1878 = buffer.data(msk0 + 1878);
    const auto *msk0_1881 = buffer.data(msk0 + 1881);
    const auto *msk0_1882 = buffer.data(msk0 + 1882);
    const auto *msk0_1884 = buffer.data(msk0 + 1884);
    const auto *msk0_1886 = buffer.data(msk0 + 1886);
    const auto *msk0_1887 = buffer.data(msk0 + 1887);
    const auto *msk0_1889 = buffer.data(msk0 + 1889);
    const auto *msk0_1890 = buffer.data(msk0 + 1890);
    const auto *msk0_1892 = buffer.data(msk0 + 1892);
    const auto *msk0_1900 = buffer.data(msk0 + 1900);
    const auto *msk0_1902 = buffer.data(msk0 + 1902);
    const auto *msk0_1903 = buffer.data(msk0 + 1903);
    const auto *msk0_1904 = buffer.data(msk0 + 1904);
    const auto *msk0_1905 = buffer.data(msk0 + 1905);
    const auto *msk0_1907 = buffer.data(msk0 + 1907);
    const auto *msk0_1911 = buffer.data(msk0 + 1911);
    const auto *msk0_1914 = buffer.data(msk0 + 1914);
    const auto *msk0_1918 = buffer.data(msk0 + 1918);
    const auto *msk0_1920 = buffer.data(msk0 + 1920);
    const auto *msk0_1923 = buffer.data(msk0 + 1923);
    const auto *msk0_1925 = buffer.data(msk0 + 1925);
    const auto *msk0_1926 = buffer.data(msk0 + 1926);
    const auto *msk0_1936 = buffer.data(msk0 + 1936);
    const auto *msk0_1938 = buffer.data(msk0 + 1938);
    const auto *msk0_1939 = buffer.data(msk0 + 1939);
    const auto *msk0_1940 = buffer.data(msk0 + 1940);
    const auto *msk0_1941 = buffer.data(msk0 + 1941);
    const auto *msk0_1943 = buffer.data(msk0 + 1943);

    const auto *msi_1141 = buffer.data(msi + 1141);
    const auto *msi_1148 = buffer.data(msi + 1148);
    const auto *msi_1151 = buffer.data(msi + 1151);
    const auto *msi_1154 = buffer.data(msi + 1154);
    const auto *msi_1158 = buffer.data(msi + 1158);
    const auto *msi_1169 = buffer.data(msi + 1169);
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
    const auto *msi_1231 = buffer.data(msi + 1231);
    const auto *msi_1232 = buffer.data(msi + 1232);
    const auto *msi_1234 = buffer.data(msi + 1234);
    const auto *msi_1237 = buffer.data(msi + 1237);
    const auto *msi_1241 = buffer.data(msi + 1241);
    const auto *msi_1246 = buffer.data(msi + 1246);
    const auto *msi_1259 = buffer.data(msi + 1259);
    const auto *msi_1426 = buffer.data(msi + 1426);
    const auto *msi_1427 = buffer.data(msi + 1427);
    const auto *msi_1428 = buffer.data(msi + 1428);
    const auto *msi_1431 = buffer.data(msi + 1431);
    const auto *msi_1433 = buffer.data(msi + 1433);
    const auto *msi_1434 = buffer.data(msi + 1434);
    const auto *msi_1437 = buffer.data(msi + 1437);
    const auto *msi_1438 = buffer.data(msi + 1438);
    const auto *msi_1440 = buffer.data(msi + 1440);
    const auto *msi_1442 = buffer.data(msi + 1442);
    const auto *msi_1443 = buffer.data(msi + 1443);
    const auto *msi_1445 = buffer.data(msi + 1445);
    const auto *msi_1446 = buffer.data(msi + 1446);
    const auto *msi_1448 = buffer.data(msi + 1448);
    const auto *msi_1449 = buffer.data(msi + 1449);
    const auto *msi_1450 = buffer.data(msi + 1450);
    const auto *msi_1451 = buffer.data(msi + 1451);
    const auto *msi_1452 = buffer.data(msi + 1452);
    const auto *msi_1453 = buffer.data(msi + 1453);
    const auto *msi_1454 = buffer.data(msi + 1454);
    const auto *msi_1455 = buffer.data(msi + 1455);
    const auto *msi_1456 = buffer.data(msi + 1456);
    const auto *msi_1459 = buffer.data(msi + 1459);
    const auto *msi_1461 = buffer.data(msi + 1461);
    const auto *msi_1462 = buffer.data(msi + 1462);
    const auto *msi_1465 = buffer.data(msi + 1465);
    const auto *msi_1466 = buffer.data(msi + 1466);
    const auto *msi_1468 = buffer.data(msi + 1468);
    const auto *msi_1470 = buffer.data(msi + 1470);
    const auto *msi_1471 = buffer.data(msi + 1471);
    const auto *msi_1473 = buffer.data(msi + 1473);
    const auto *msi_1474 = buffer.data(msi + 1474);
    const auto *msi_1476 = buffer.data(msi + 1476);
    const auto *msi_1477 = buffer.data(msi + 1477);
    const auto *msi_1478 = buffer.data(msi + 1478);
    const auto *msi_1479 = buffer.data(msi + 1479);
    const auto *msi_1480 = buffer.data(msi + 1480);
    const auto *msi_1481 = buffer.data(msi + 1481);
    const auto *msi_1482 = buffer.data(msi + 1482);
    const auto *msi_1483 = buffer.data(msi + 1483);
    const auto *msi_1487 = buffer.data(msi + 1487);
    const auto *msi_1490 = buffer.data(msi + 1490);
    const auto *msi_1494 = buffer.data(msi + 1494);
    const auto *msi_1496 = buffer.data(msi + 1496);
    const auto *msi_1499 = buffer.data(msi + 1499);
    const auto *msi_1501 = buffer.data(msi + 1501);
    const auto *msi_1502 = buffer.data(msi + 1502);
    const auto *msi_1505 = buffer.data(msi + 1505);
    const auto *msi_1506 = buffer.data(msi + 1506);
    const auto *msi_1507 = buffer.data(msi + 1507);
    const auto *msi_1508 = buffer.data(msi + 1508);
    const auto *msi_1509 = buffer.data(msi + 1509);
    const auto *msi_1510 = buffer.data(msi + 1510);
    const auto *msi_1511 = buffer.data(msi + 1511);

    const auto *msk1_1584 = buffer.data(msk1 + 1584);
    const auto *msk1_1589 = buffer.data(msk1 + 1589);
    const auto *msk1_1593 = buffer.data(msk1 + 1593);
    const auto *msk1_1598 = buffer.data(msk1 + 1598);
    const auto *msk1_1604 = buffer.data(msk1 + 1604);
    const auto *msk1_1828 = buffer.data(msk1 + 1828);
    const auto *msk1_1830 = buffer.data(msk1 + 1830);
    const auto *msk1_1831 = buffer.data(msk1 + 1831);
    const auto *msk1_1832 = buffer.data(msk1 + 1832);
    const auto *msk1_1833 = buffer.data(msk1 + 1833);
    const auto *msk1_1835 = buffer.data(msk1 + 1835);
    const auto *msk1_1836 = buffer.data(msk1 + 1836);
    const auto *msk1_1839 = buffer.data(msk1 + 1839);
    const auto *msk1_1841 = buffer.data(msk1 + 1841);
    const auto *msk1_1842 = buffer.data(msk1 + 1842);
    const auto *msk1_1845 = buffer.data(msk1 + 1845);
    const auto *msk1_1846 = buffer.data(msk1 + 1846);
    const auto *msk1_1848 = buffer.data(msk1 + 1848);
    const auto *msk1_1850 = buffer.data(msk1 + 1850);
    const auto *msk1_1851 = buffer.data(msk1 + 1851);
    const auto *msk1_1853 = buffer.data(msk1 + 1853);
    const auto *msk1_1854 = buffer.data(msk1 + 1854);
    const auto *msk1_1856 = buffer.data(msk1 + 1856);
    const auto *msk1_1864 = buffer.data(msk1 + 1864);
    const auto *msk1_1866 = buffer.data(msk1 + 1866);
    const auto *msk1_1867 = buffer.data(msk1 + 1867);
    const auto *msk1_1868 = buffer.data(msk1 + 1868);
    const auto *msk1_1869 = buffer.data(msk1 + 1869);
    const auto *msk1_1871 = buffer.data(msk1 + 1871);
    const auto *msk1_1872 = buffer.data(msk1 + 1872);
    const auto *msk1_1875 = buffer.data(msk1 + 1875);
    const auto *msk1_1877 = buffer.data(msk1 + 1877);
    const auto *msk1_1878 = buffer.data(msk1 + 1878);
    const auto *msk1_1881 = buffer.data(msk1 + 1881);
    const auto *msk1_1882 = buffer.data(msk1 + 1882);
    const auto *msk1_1884 = buffer.data(msk1 + 1884);
    const auto *msk1_1886 = buffer.data(msk1 + 1886);
    const auto *msk1_1887 = buffer.data(msk1 + 1887);
    const auto *msk1_1889 = buffer.data(msk1 + 1889);
    const auto *msk1_1890 = buffer.data(msk1 + 1890);
    const auto *msk1_1892 = buffer.data(msk1 + 1892);
    const auto *msk1_1900 = buffer.data(msk1 + 1900);
    const auto *msk1_1902 = buffer.data(msk1 + 1902);
    const auto *msk1_1903 = buffer.data(msk1 + 1903);
    const auto *msk1_1904 = buffer.data(msk1 + 1904);
    const auto *msk1_1905 = buffer.data(msk1 + 1905);
    const auto *msk1_1907 = buffer.data(msk1 + 1907);
    const auto *msk1_1911 = buffer.data(msk1 + 1911);
    const auto *msk1_1914 = buffer.data(msk1 + 1914);
    const auto *msk1_1918 = buffer.data(msk1 + 1918);
    const auto *msk1_1920 = buffer.data(msk1 + 1920);
    const auto *msk1_1923 = buffer.data(msk1 + 1923);
    const auto *msk1_1925 = buffer.data(msk1 + 1925);
    const auto *msk1_1926 = buffer.data(msk1 + 1926);
    const auto *msk1_1936 = buffer.data(msk1 + 1936);
    const auto *msk1_1938 = buffer.data(msk1 + 1938);
    const auto *msk1_1939 = buffer.data(msk1 + 1939);
    const auto *msk1_1940 = buffer.data(msk1 + 1940);
    const auto *msk1_1941 = buffer.data(msk1 + 1941);
    const auto *msk1_1943 = buffer.data(msk1 + 1943);

    const auto *nsi_1421 = buffer.data(nsi + 1421);
    const auto *nsi_1426 = buffer.data(nsi + 1426);
    const auto *nsi_1427 = buffer.data(nsi + 1427);
    const auto *nsi_1428 = buffer.data(nsi + 1428);
    const auto *nsi_1430 = buffer.data(nsi + 1430);
    const auto *nsi_1431 = buffer.data(nsi + 1431);
    const auto *nsi_1433 = buffer.data(nsi + 1433);
    const auto *nsi_1434 = buffer.data(nsi + 1434);
    const auto *nsi_1437 = buffer.data(nsi + 1437);
    const auto *nsi_1438 = buffer.data(nsi + 1438);
    const auto *nsi_1442 = buffer.data(nsi + 1442);
    const auto *nsi_1449 = buffer.data(nsi + 1449);
    const auto *nsi_1450 = buffer.data(nsi + 1450);
    const auto *nsi_1451 = buffer.data(nsi + 1451);
    const auto *nsi_1452 = buffer.data(nsi + 1452);
    const auto *nsi_1453 = buffer.data(nsi + 1453);
    const auto *nsi_1454 = buffer.data(nsi + 1454);
    const auto *nsi_1455 = buffer.data(nsi + 1455);
    const auto *nsi_1456 = buffer.data(nsi + 1456);
    const auto *nsi_1458 = buffer.data(nsi + 1458);
    const auto *nsi_1459 = buffer.data(nsi + 1459);
    const auto *nsi_1461 = buffer.data(nsi + 1461);
    const auto *nsi_1462 = buffer.data(nsi + 1462);
    const auto *nsi_1465 = buffer.data(nsi + 1465);
    const auto *nsi_1466 = buffer.data(nsi + 1466);
    const auto *nsi_1470 = buffer.data(nsi + 1470);
    const auto *nsi_1477 = buffer.data(nsi + 1477);
    const auto *nsi_1478 = buffer.data(nsi + 1478);
    const auto *nsi_1479 = buffer.data(nsi + 1479);
    const auto *nsi_1480 = buffer.data(nsi + 1480);
    const auto *nsi_1481 = buffer.data(nsi + 1481);
    const auto *nsi_1482 = buffer.data(nsi + 1482);
    const auto *nsi_1483 = buffer.data(nsi + 1483);
    const auto *nsi_1484 = buffer.data(nsi + 1484);
    const auto *nsi_1486 = buffer.data(nsi + 1486);
    const auto *nsi_1487 = buffer.data(nsi + 1487);
    const auto *nsi_1489 = buffer.data(nsi + 1489);
    const auto *nsi_1490 = buffer.data(nsi + 1490);
    const auto *nsi_1493 = buffer.data(nsi + 1493);
    const auto *nsi_1494 = buffer.data(nsi + 1494);
    const auto *nsi_1498 = buffer.data(nsi + 1498);
    const auto *nsi_1505 = buffer.data(nsi + 1505);
    const auto *nsi_1506 = buffer.data(nsi + 1506);
    const auto *nsi_1507 = buffer.data(nsi + 1507);
    const auto *nsi_1508 = buffer.data(nsi + 1508);
    const auto *nsi_1509 = buffer.data(nsi + 1509);
    const auto *nsi_1510 = buffer.data(nsi + 1510);
    const auto *nsi_1511 = buffer.data(nsi + 1511);

#pragma omp simd aligned(t_1826, t_1827, t_1828, t_1829, pa_x, pc_x, pc_z, msk0_1828, \
                         msi_1141, msi_1426, msi_1427, msk1_1828, nsi_1421, nsi_1426, \
                         nsi_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1826[k] = f_13 * msi_1426[k]
                    + f_3 * pc_x[k] * nsi_1426[k];

        t_1827[k] = f_13 * msi_1427[k]
                    + f_3 * pc_x[k] * nsi_1427[k];

        t_1828[k] = pa_x[k] * msk0_1828[k]
                    - f_12 * pc_x[k] * msk1_1828[k];

        t_1829[k] = f_17 * msi_1141[k]
                    + f_3 * pc_z[k] * nsi_1421[k];
    }

#pragma omp simd aligned(t_1830, t_1831, t_1832, t_1833, pa_x, pc_x, msk0_1830, msk0_1831, \
                         msk0_1832, msk0_1833, msk1_1830, msk1_1831, msk1_1832, \
                         msk1_1833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1830[k] = pa_x[k] * msk0_1830[k]
                    - f_12 * pc_x[k] * msk1_1830[k];

        t_1831[k] = pa_x[k] * msk0_1831[k]
                    - f_12 * pc_x[k] * msk1_1831[k];

        t_1832[k] = pa_x[k] * msk0_1832[k]
                    - f_12 * pc_x[k] * msk1_1832[k];

        t_1833[k] = pa_x[k] * msk0_1833[k]
                    - f_12 * pc_x[k] * msk1_1833[k];
    }

#pragma omp simd aligned(t_1834, t_1835, t_1836, t_1837, pa_x, pc_x, pc_y, msk0_1835, \
                         msk0_1836, msi_1175, msi_1176, msi_1428, msk1_1835, msk1_1836, \
                         nsi_1427, nsi_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = f_16 * msi_1175[k]
                    + f_3 * pc_y[k] * nsi_1427[k];

        t_1835[k] = pa_x[k] * msk0_1835[k]
                    - f_12 * pc_x[k] * msk1_1835[k];

        t_1836[k] = pa_x[k] * msk0_1836[k]
                    + f_22 * msi_1428[k]
                    - f_12 * pc_x[k] * msk1_1836[k];

        t_1837[k] = f_15 * msi_1176[k]
                    + f_3 * pc_y[k] * nsi_1428[k];
    }

#pragma omp simd aligned(t_1838, t_1839, t_1840, pa_x, pc_x, pc_y, pc_z, msk0_1839, msi_1148, \
                         msi_1178, msi_1431, msk1_1839, nsi_1428, \
                         nsi_1430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1838[k] = f_23 * msi_1148[k]
                    + f_3 * pc_z[k] * nsi_1428[k];

        t_1839[k] = pa_x[k] * msk0_1839[k]
                    + f_17 * msi_1431[k]
                    - f_12 * pc_x[k] * msk1_1839[k];

        t_1840[k] = f_15 * msi_1178[k]
                    + f_3 * pc_y[k] * nsi_1430[k];
    }

#pragma omp simd aligned(t_1841, t_1842, t_1843, pa_x, pc_x, pc_z, msk0_1841, msk0_1842, \
                         msi_1151, msi_1433, msi_1434, msk1_1841, msk1_1842, \
                         nsi_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1841[k] = pa_x[k] * msk0_1841[k]
                    + f_17 * msi_1433[k]
                    - f_12 * pc_x[k] * msk1_1841[k];

        t_1842[k] = pa_x[k] * msk0_1842[k]
                    + f_16 * msi_1434[k]
                    - f_12 * pc_x[k] * msk1_1842[k];

        t_1843[k] = f_23 * msi_1151[k]
                    + f_3 * pc_z[k] * nsi_1431[k];
    }

#pragma omp simd aligned(t_1844, t_1845, t_1846, pa_x, pc_x, pc_y, msk0_1845, msk0_1846, \
                         msi_1181, msi_1437, msi_1438, msk1_1845, msk1_1846, \
                         nsi_1433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1844[k] = f_15 * msi_1181[k]
                    + f_3 * pc_y[k] * nsi_1433[k];

        t_1845[k] = pa_x[k] * msk0_1845[k]
                    + f_16 * msi_1437[k]
                    - f_12 * pc_x[k] * msk1_1845[k];

        t_1846[k] = pa_x[k] * msk0_1846[k]
                    + f_15 * msi_1438[k]
                    - f_12 * pc_x[k] * msk1_1846[k];
    }

#pragma omp simd aligned(t_1847, t_1848, t_1849, pa_x, pc_x, pc_y, pc_z, msk0_1848, msi_1154, \
                         msi_1185, msi_1440, msk1_1848, nsi_1434, \
                         nsi_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1847[k] = f_23 * msi_1154[k]
                    + f_3 * pc_z[k] * nsi_1434[k];

        t_1848[k] = pa_x[k] * msk0_1848[k]
                    + f_15 * msi_1440[k]
                    - f_12 * pc_x[k] * msk1_1848[k];

        t_1849[k] = f_15 * msi_1185[k]
                    + f_3 * pc_y[k] * nsi_1437[k];
    }

#pragma omp simd aligned(t_1850, t_1851, t_1852, pa_x, pc_x, pc_z, msk0_1850, msk0_1851, \
                         msi_1158, msi_1442, msi_1443, msk1_1850, msk1_1851, \
                         nsi_1438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1850[k] = pa_x[k] * msk0_1850[k]
                    + f_15 * msi_1442[k]
                    - f_12 * pc_x[k] * msk1_1850[k];

        t_1851[k] = pa_x[k] * msk0_1851[k]
                    + f_14 * msi_1443[k]
                    - f_12 * pc_x[k] * msk1_1851[k];

        t_1852[k] = f_23 * msi_1158[k]
                    + f_3 * pc_z[k] * nsi_1438[k];
    }

#pragma omp simd aligned(t_1853, t_1854, t_1855, pa_x, pc_x, pc_y, msk0_1853, msk0_1854, \
                         msi_1190, msi_1445, msi_1446, msk1_1853, msk1_1854, \
                         nsi_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1853[k] = pa_x[k] * msk0_1853[k]
                    + f_14 * msi_1445[k]
                    - f_12 * pc_x[k] * msk1_1853[k];

        t_1854[k] = pa_x[k] * msk0_1854[k]
                    + f_14 * msi_1446[k]
                    - f_12 * pc_x[k] * msk1_1854[k];

        t_1855[k] = f_15 * msi_1190[k]
                    + f_3 * pc_y[k] * nsi_1442[k];
    }

#pragma omp simd aligned(t_1856, t_1857, t_1858, t_1859, pa_x, pc_x, msk0_1856, msi_1448, \
                         msi_1449, msi_1450, msi_1451, msk1_1856, nsi_1449, nsi_1450, \
                         nsi_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1856[k] = pa_x[k] * msk0_1856[k]
                    + f_14 * msi_1448[k]
                    - f_12 * pc_x[k] * msk1_1856[k];

        t_1857[k] = f_13 * msi_1449[k]
                    + f_3 * pc_x[k] * nsi_1449[k];

        t_1858[k] = f_13 * msi_1450[k]
                    + f_3 * pc_x[k] * nsi_1450[k];

        t_1859[k] = f_13 * msi_1451[k]
                    + f_3 * pc_x[k] * nsi_1451[k];
    }

#pragma omp simd aligned(t_1860, t_1861, t_1862, t_1863, pc_x, msi_1452, msi_1453, msi_1454, \
                         msi_1455, nsi_1452, nsi_1453, nsi_1454, \
                         nsi_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1860[k] = f_13 * msi_1452[k]
                    + f_3 * pc_x[k] * nsi_1452[k];

        t_1861[k] = f_13 * msi_1453[k]
                    + f_3 * pc_x[k] * nsi_1453[k];

        t_1862[k] = f_13 * msi_1454[k]
                    + f_3 * pc_x[k] * nsi_1454[k];

        t_1863[k] = f_13 * msi_1455[k]
                    + f_3 * pc_x[k] * nsi_1455[k];
    }

#pragma omp simd aligned(t_1864, t_1865, t_1866, t_1867, pa_x, pc_x, pc_z, msk0_1864, \
                         msk0_1866, msk0_1867, msi_1169, msk1_1864, msk1_1866, msk1_1867, \
                         nsi_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1864[k] = pa_x[k] * msk0_1864[k]
                    - f_12 * pc_x[k] * msk1_1864[k];

        t_1865[k] = f_23 * msi_1169[k]
                    + f_3 * pc_z[k] * nsi_1449[k];

        t_1866[k] = pa_x[k] * msk0_1866[k]
                    - f_12 * pc_x[k] * msk1_1866[k];

        t_1867[k] = pa_x[k] * msk0_1867[k]
                    - f_12 * pc_x[k] * msk1_1867[k];
    }

#pragma omp simd aligned(t_1868, t_1869, t_1870, t_1871, pa_x, pc_x, pc_y, msk0_1868, \
                         msk0_1869, msk0_1871, msi_1203, msk1_1868, msk1_1869, msk1_1871, \
                         nsi_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1868[k] = pa_x[k] * msk0_1868[k]
                    - f_12 * pc_x[k] * msk1_1868[k];

        t_1869[k] = pa_x[k] * msk0_1869[k]
                    - f_12 * pc_x[k] * msk1_1869[k];

        t_1870[k] = f_15 * msi_1203[k]
                    + f_3 * pc_y[k] * nsi_1455[k];

        t_1871[k] = pa_x[k] * msk0_1871[k]
                    - f_12 * pc_x[k] * msk1_1871[k];
    }

#pragma omp simd aligned(t_1872, t_1873, t_1874, pa_x, pc_x, pc_y, pc_z, msk0_1872, msi_1176, \
                         msi_1204, msi_1456, msk1_1872, nsi_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1872[k] = pa_x[k] * msk0_1872[k]
                    + f_22 * msi_1456[k]
                    - f_12 * pc_x[k] * msk1_1872[k];

        t_1873[k] = f_14 * msi_1204[k]
                    + f_3 * pc_y[k] * nsi_1456[k];

        t_1874[k] = f_22 * msi_1176[k]
                    + f_3 * pc_z[k] * nsi_1456[k];
    }

#pragma omp simd aligned(t_1875, t_1876, t_1877, pa_x, pc_x, pc_y, msk0_1875, msk0_1877, \
                         msi_1206, msi_1459, msi_1461, msk1_1875, msk1_1877, \
                         nsi_1458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1875[k] = pa_x[k] * msk0_1875[k]
                    + f_17 * msi_1459[k]
                    - f_12 * pc_x[k] * msk1_1875[k];

        t_1876[k] = f_14 * msi_1206[k]
                    + f_3 * pc_y[k] * nsi_1458[k];

        t_1877[k] = pa_x[k] * msk0_1877[k]
                    + f_17 * msi_1461[k]
                    - f_12 * pc_x[k] * msk1_1877[k];
    }

#pragma omp simd aligned(t_1878, t_1879, t_1880, pa_x, pc_x, pc_y, pc_z, msk0_1878, msi_1179, \
                         msi_1209, msi_1462, msk1_1878, nsi_1459, \
                         nsi_1461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1878[k] = pa_x[k] * msk0_1878[k]
                    + f_16 * msi_1462[k]
                    - f_12 * pc_x[k] * msk1_1878[k];

        t_1879[k] = f_22 * msi_1179[k]
                    + f_3 * pc_z[k] * nsi_1459[k];

        t_1880[k] = f_14 * msi_1209[k]
                    + f_3 * pc_y[k] * nsi_1461[k];
    }

#pragma omp simd aligned(t_1881, t_1882, t_1883, pa_x, pc_x, pc_z, msk0_1881, msk0_1882, \
                         msi_1182, msi_1465, msi_1466, msk1_1881, msk1_1882, \
                         nsi_1462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1881[k] = pa_x[k] * msk0_1881[k]
                    + f_16 * msi_1465[k]
                    - f_12 * pc_x[k] * msk1_1881[k];

        t_1882[k] = pa_x[k] * msk0_1882[k]
                    + f_15 * msi_1466[k]
                    - f_12 * pc_x[k] * msk1_1882[k];

        t_1883[k] = f_22 * msi_1182[k]
                    + f_3 * pc_z[k] * nsi_1462[k];
    }

#pragma omp simd aligned(t_1884, t_1885, t_1886, pa_x, pc_x, pc_y, msk0_1884, msk0_1886, \
                         msi_1213, msi_1468, msi_1470, msk1_1884, msk1_1886, \
                         nsi_1465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1884[k] = pa_x[k] * msk0_1884[k]
                    + f_15 * msi_1468[k]
                    - f_12 * pc_x[k] * msk1_1884[k];

        t_1885[k] = f_14 * msi_1213[k]
                    + f_3 * pc_y[k] * nsi_1465[k];

        t_1886[k] = pa_x[k] * msk0_1886[k]
                    + f_15 * msi_1470[k]
                    - f_12 * pc_x[k] * msk1_1886[k];
    }

#pragma omp simd aligned(t_1887, t_1888, t_1889, pa_x, pc_x, pc_z, msk0_1887, msk0_1889, \
                         msi_1186, msi_1471, msi_1473, msk1_1887, msk1_1889, \
                         nsi_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1887[k] = pa_x[k] * msk0_1887[k]
                    + f_14 * msi_1471[k]
                    - f_12 * pc_x[k] * msk1_1887[k];

        t_1888[k] = f_22 * msi_1186[k]
                    + f_3 * pc_z[k] * nsi_1466[k];

        t_1889[k] = pa_x[k] * msk0_1889[k]
                    + f_14 * msi_1473[k]
                    - f_12 * pc_x[k] * msk1_1889[k];
    }

#pragma omp simd aligned(t_1890, t_1891, t_1892, pa_x, pc_x, pc_y, msk0_1890, msk0_1892, \
                         msi_1218, msi_1474, msi_1476, msk1_1890, msk1_1892, \
                         nsi_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1890[k] = pa_x[k] * msk0_1890[k]
                    + f_14 * msi_1474[k]
                    - f_12 * pc_x[k] * msk1_1890[k];

        t_1891[k] = f_14 * msi_1218[k]
                    + f_3 * pc_y[k] * nsi_1470[k];

        t_1892[k] = pa_x[k] * msk0_1892[k]
                    + f_14 * msi_1476[k]
                    - f_12 * pc_x[k] * msk1_1892[k];
    }

#pragma omp simd aligned(t_1893, t_1894, t_1895, t_1896, t_1897, pc_x, msi_1477, msi_1478, \
                         msi_1479, msi_1480, msi_1481, nsi_1477, nsi_1478, nsi_1479, nsi_1480, \
                         nsi_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1893[k] = f_13 * msi_1477[k]
                    + f_3 * pc_x[k] * nsi_1477[k];

        t_1894[k] = f_13 * msi_1478[k]
                    + f_3 * pc_x[k] * nsi_1478[k];

        t_1895[k] = f_13 * msi_1479[k]
                    + f_3 * pc_x[k] * nsi_1479[k];

        t_1896[k] = f_13 * msi_1480[k]
                    + f_3 * pc_x[k] * nsi_1480[k];

        t_1897[k] = f_13 * msi_1481[k]
                    + f_3 * pc_x[k] * nsi_1481[k];
    }

#pragma omp simd aligned(t_1898, t_1899, t_1900, t_1901, pa_x, pc_x, pc_z, msk0_1900, \
                         msi_1197, msi_1482, msi_1483, msk1_1900, nsi_1477, nsi_1482, \
                         nsi_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1898[k] = f_13 * msi_1482[k]
                    + f_3 * pc_x[k] * nsi_1482[k];

        t_1899[k] = f_13 * msi_1483[k]
                    + f_3 * pc_x[k] * nsi_1483[k];

        t_1900[k] = pa_x[k] * msk0_1900[k]
                    - f_12 * pc_x[k] * msk1_1900[k];

        t_1901[k] = f_22 * msi_1197[k]
                    + f_3 * pc_z[k] * nsi_1477[k];
    }

#pragma omp simd aligned(t_1902, t_1903, t_1904, t_1905, pa_x, pc_x, msk0_1902, msk0_1903, \
                         msk0_1904, msk0_1905, msk1_1902, msk1_1903, msk1_1904, \
                         msk1_1905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1902[k] = pa_x[k] * msk0_1902[k]
                    - f_12 * pc_x[k] * msk1_1902[k];

        t_1903[k] = pa_x[k] * msk0_1903[k]
                    - f_12 * pc_x[k] * msk1_1903[k];

        t_1904[k] = pa_x[k] * msk0_1904[k]
                    - f_12 * pc_x[k] * msk1_1904[k];

        t_1905[k] = pa_x[k] * msk0_1905[k]
                    - f_12 * pc_x[k] * msk1_1905[k];
    }

#pragma omp simd aligned(t_1906, t_1907, t_1908, t_1909, pa_x, pa_y, pc_x, pc_y, msk0_1584, \
                         msk0_1907, msi_1231, msi_1232, msk1_1584, msk1_1907, nsi_1483, \
                         nsi_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1906[k] = f_14 * msi_1231[k]
                    + f_3 * pc_y[k] * nsi_1483[k];

        t_1907[k] = pa_x[k] * msk0_1907[k]
                    - f_12 * pc_x[k] * msk1_1907[k];

        t_1908[k] = pa_y[k] * msk0_1584[k]
                    - f_12 * pc_y[k] * msk1_1584[k];

        t_1909[k] = f_13 * msi_1232[k]
                    + f_3 * pc_y[k] * nsi_1484[k];
    }

#pragma omp simd aligned(t_1910, t_1911, t_1912, pa_x, pc_x, pc_y, pc_z, msk0_1911, msi_1204, \
                         msi_1234, msi_1487, msk1_1911, nsi_1484, \
                         nsi_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1910[k] = f_21 * msi_1204[k]
                    + f_3 * pc_z[k] * nsi_1484[k];

        t_1911[k] = pa_x[k] * msk0_1911[k]
                    + f_17 * msi_1487[k]
                    - f_12 * pc_x[k] * msk1_1911[k];

        t_1912[k] = f_13 * msi_1234[k]
                    + f_3 * pc_y[k] * nsi_1486[k];
    }

#pragma omp simd aligned(t_1913, t_1914, t_1915, pa_x, pa_y, pc_x, pc_y, pc_z, msk0_1589, \
                         msk0_1914, msi_1207, msi_1490, msk1_1589, msk1_1914, \
                         nsi_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1913[k] = pa_y[k] * msk0_1589[k]
                    - f_12 * pc_y[k] * msk1_1589[k];

        t_1914[k] = pa_x[k] * msk0_1914[k]
                    + f_16 * msi_1490[k]
                    - f_12 * pc_x[k] * msk1_1914[k];

        t_1915[k] = f_21 * msi_1207[k]
                    + f_3 * pc_z[k] * nsi_1487[k];
    }

#pragma omp simd aligned(t_1916, t_1917, t_1918, pa_x, pa_y, pc_x, pc_y, msk0_1593, msk0_1918, \
                         msi_1237, msi_1494, msk1_1593, msk1_1918, \
                         nsi_1489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1916[k] = f_13 * msi_1237[k]
                    + f_3 * pc_y[k] * nsi_1489[k];

        t_1917[k] = pa_y[k] * msk0_1593[k]
                    - f_12 * pc_y[k] * msk1_1593[k];

        t_1918[k] = pa_x[k] * msk0_1918[k]
                    + f_15 * msi_1494[k]
                    - f_12 * pc_x[k] * msk1_1918[k];
    }

#pragma omp simd aligned(t_1919, t_1920, t_1921, pa_x, pc_x, pc_y, pc_z, msk0_1920, msi_1210, \
                         msi_1241, msi_1496, msk1_1920, nsi_1490, \
                         nsi_1493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1919[k] = f_21 * msi_1210[k]
                    + f_3 * pc_z[k] * nsi_1490[k];

        t_1920[k] = pa_x[k] * msk0_1920[k]
                    + f_15 * msi_1496[k]
                    - f_12 * pc_x[k] * msk1_1920[k];

        t_1921[k] = f_13 * msi_1241[k]
                    + f_3 * pc_y[k] * nsi_1493[k];
    }

#pragma omp simd aligned(t_1922, t_1923, t_1924, pa_x, pa_y, pc_x, pc_y, pc_z, msk0_1598, \
                         msk0_1923, msi_1214, msi_1499, msk1_1598, msk1_1923, \
                         nsi_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1922[k] = pa_y[k] * msk0_1598[k]
                    - f_12 * pc_y[k] * msk1_1598[k];

        t_1923[k] = pa_x[k] * msk0_1923[k]
                    + f_14 * msi_1499[k]
                    - f_12 * pc_x[k] * msk1_1923[k];

        t_1924[k] = f_21 * msi_1214[k]
                    + f_3 * pc_z[k] * nsi_1494[k];
    }

#pragma omp simd aligned(t_1925, t_1926, t_1927, pa_x, pc_x, pc_y, msk0_1925, msk0_1926, \
                         msi_1246, msi_1501, msi_1502, msk1_1925, msk1_1926, \
                         nsi_1498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1925[k] = pa_x[k] * msk0_1925[k]
                    + f_14 * msi_1501[k]
                    - f_12 * pc_x[k] * msk1_1925[k];

        t_1926[k] = pa_x[k] * msk0_1926[k]
                    + f_14 * msi_1502[k]
                    - f_12 * pc_x[k] * msk1_1926[k];

        t_1927[k] = f_13 * msi_1246[k]
                    + f_3 * pc_y[k] * nsi_1498[k];
    }

#pragma omp simd aligned(t_1928, t_1929, t_1930, t_1931, pa_y, pc_x, pc_y, msk0_1604, \
                         msi_1505, msi_1506, msi_1507, msk1_1604, nsi_1505, nsi_1506, \
                         nsi_1507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1928[k] = pa_y[k] * msk0_1604[k]
                    - f_12 * pc_y[k] * msk1_1604[k];

        t_1929[k] = f_13 * msi_1505[k]
                    + f_3 * pc_x[k] * nsi_1505[k];

        t_1930[k] = f_13 * msi_1506[k]
                    + f_3 * pc_x[k] * nsi_1506[k];

        t_1931[k] = f_13 * msi_1507[k]
                    + f_3 * pc_x[k] * nsi_1507[k];
    }

#pragma omp simd aligned(t_1932, t_1933, t_1934, t_1935, pc_x, msi_1508, msi_1509, msi_1510, \
                         msi_1511, nsi_1508, nsi_1509, nsi_1510, \
                         nsi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1932[k] = f_13 * msi_1508[k]
                    + f_3 * pc_x[k] * nsi_1508[k];

        t_1933[k] = f_13 * msi_1509[k]
                    + f_3 * pc_x[k] * nsi_1509[k];

        t_1934[k] = f_13 * msi_1510[k]
                    + f_3 * pc_x[k] * nsi_1510[k];

        t_1935[k] = f_13 * msi_1511[k]
                    + f_3 * pc_x[k] * nsi_1511[k];
    }

#pragma omp simd aligned(t_1936, t_1937, t_1938, t_1939, pa_x, pc_x, pc_z, msk0_1936, \
                         msk0_1938, msk0_1939, msi_1225, msk1_1936, msk1_1938, msk1_1939, \
                         nsi_1505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1936[k] = pa_x[k] * msk0_1936[k]
                    - f_12 * pc_x[k] * msk1_1936[k];

        t_1937[k] = f_21 * msi_1225[k]
                    + f_3 * pc_z[k] * nsi_1505[k];

        t_1938[k] = pa_x[k] * msk0_1938[k]
                    - f_12 * pc_x[k] * msk1_1938[k];

        t_1939[k] = pa_x[k] * msk0_1939[k]
                    - f_12 * pc_x[k] * msk1_1939[k];
    }

#pragma omp simd aligned(t_1940, t_1941, t_1942, t_1943, pa_x, pc_x, pc_y, msk0_1940, \
                         msk0_1941, msk0_1943, msi_1259, msk1_1940, msk1_1941, msk1_1943, \
                         nsi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1940[k] = pa_x[k] * msk0_1940[k]
                    - f_12 * pc_x[k] * msk1_1940[k];

        t_1941[k] = pa_x[k] * msk0_1941[k]
                    - f_12 * pc_x[k] * msk1_1941[k];

        t_1942[k] = f_13 * msi_1259[k]
                    + f_3 * pc_y[k] * nsi_1511[k];

        t_1943[k] = pa_x[k] * msk0_1943[k]
                    - f_12 * pc_x[k] * msk1_1943[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsh0, const size_t nsh1,
                                                           const size_t nsi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_22 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1620 = buffer.data(msk0 + 1620);
    const auto *msk0_1621 = buffer.data(msk0 + 1621);
    const auto *msk0_1623 = buffer.data(msk0 + 1623);
    const auto *msk0_1626 = buffer.data(msk0 + 1626);
    const auto *msk0_1630 = buffer.data(msk0 + 1630);
    const auto *msk0_1635 = buffer.data(msk0 + 1635);
    const auto *msk0_1648 = buffer.data(msk0 + 1648);
    const auto *msk0_1650 = buffer.data(msk0 + 1650);
    const auto *msk0_1651 = buffer.data(msk0 + 1651);
    const auto *msk0_1652 = buffer.data(msk0 + 1652);
    const auto *msk0_1653 = buffer.data(msk0 + 1653);
    const auto *msk0_1944 = buffer.data(msk0 + 1944);
    const auto *msk0_1949 = buffer.data(msk0 + 1949);
    const auto *msk0_1953 = buffer.data(msk0 + 1953);
    const auto *msk0_1958 = buffer.data(msk0 + 1958);
    const auto *msk0_1964 = buffer.data(msk0 + 1964);
    const auto *msk0_1972 = buffer.data(msk0 + 1972);
    const auto *msk0_1973 = buffer.data(msk0 + 1973);
    const auto *msk0_1974 = buffer.data(msk0 + 1974);
    const auto *msk0_1975 = buffer.data(msk0 + 1975);
    const auto *msk0_1976 = buffer.data(msk0 + 1976);
    const auto *msk0_1977 = buffer.data(msk0 + 1977);
    const auto *msk0_1979 = buffer.data(msk0 + 1979);

    const auto *msi_1232 = buffer.data(msi + 1232);
    const auto *msi_1281 = buffer.data(msi + 1281);
    const auto *msi_1282 = buffer.data(msi + 1282);
    const auto *msi_1283 = buffer.data(msi + 1283);
    const auto *msi_1284 = buffer.data(msi + 1284);
    const auto *msi_1285 = buffer.data(msi + 1285);
    const auto *msi_1287 = buffer.data(msi + 1287);
    const auto *msi_1315 = buffer.data(msi + 1315);
    const auto *msi_1512 = buffer.data(msi + 1512);
    const auto *msi_1517 = buffer.data(msi + 1517);
    const auto *msi_1521 = buffer.data(msi + 1521);
    const auto *msi_1526 = buffer.data(msi + 1526);
    const auto *msi_1532 = buffer.data(msi + 1532);
    const auto *msi_1533 = buffer.data(msi + 1533);
    const auto *msi_1534 = buffer.data(msi + 1534);
    const auto *msi_1535 = buffer.data(msi + 1535);
    const auto *msi_1536 = buffer.data(msi + 1536);
    const auto *msi_1537 = buffer.data(msi + 1537);
    const auto *msi_1539 = buffer.data(msi + 1539);

    const auto *msk1_1620 = buffer.data(msk1 + 1620);
    const auto *msk1_1621 = buffer.data(msk1 + 1621);
    const auto *msk1_1623 = buffer.data(msk1 + 1623);
    const auto *msk1_1626 = buffer.data(msk1 + 1626);
    const auto *msk1_1630 = buffer.data(msk1 + 1630);
    const auto *msk1_1635 = buffer.data(msk1 + 1635);
    const auto *msk1_1648 = buffer.data(msk1 + 1648);
    const auto *msk1_1650 = buffer.data(msk1 + 1650);
    const auto *msk1_1651 = buffer.data(msk1 + 1651);
    const auto *msk1_1652 = buffer.data(msk1 + 1652);
    const auto *msk1_1653 = buffer.data(msk1 + 1653);
    const auto *msk1_1944 = buffer.data(msk1 + 1944);
    const auto *msk1_1949 = buffer.data(msk1 + 1949);
    const auto *msk1_1953 = buffer.data(msk1 + 1953);
    const auto *msk1_1958 = buffer.data(msk1 + 1958);
    const auto *msk1_1964 = buffer.data(msk1 + 1964);
    const auto *msk1_1972 = buffer.data(msk1 + 1972);
    const auto *msk1_1973 = buffer.data(msk1 + 1973);
    const auto *msk1_1974 = buffer.data(msk1 + 1974);
    const auto *msk1_1975 = buffer.data(msk1 + 1975);
    const auto *msk1_1976 = buffer.data(msk1 + 1976);
    const auto *msk1_1977 = buffer.data(msk1 + 1977);
    const auto *msk1_1979 = buffer.data(msk1 + 1979);

    const auto *nsh0_1134 = buffer.data(nsh0 + 1134);
    const auto *nsh0_1135 = buffer.data(nsh0 + 1135);
    const auto *nsh0_1136 = buffer.data(nsh0 + 1136);
    const auto *nsh0_1137 = buffer.data(nsh0 + 1137);
    const auto *nsh0_1138 = buffer.data(nsh0 + 1138);
    const auto *nsh0_1139 = buffer.data(nsh0 + 1139);
    const auto *nsh0_1140 = buffer.data(nsh0 + 1140);
    const auto *nsh0_1141 = buffer.data(nsh0 + 1141);
    const auto *nsh0_1142 = buffer.data(nsh0 + 1142);
    const auto *nsh0_1143 = buffer.data(nsh0 + 1143);
    const auto *nsh0_1155 = buffer.data(nsh0 + 1155);
    const auto *nsh0_1156 = buffer.data(nsh0 + 1156);
    const auto *nsh0_1158 = buffer.data(nsh0 + 1158);
    const auto *nsh0_1160 = buffer.data(nsh0 + 1160);
    const auto *nsh0_1161 = buffer.data(nsh0 + 1161);
    const auto *nsh0_1163 = buffer.data(nsh0 + 1163);
    const auto *nsh0_1164 = buffer.data(nsh0 + 1164);
    const auto *nsh0_1165 = buffer.data(nsh0 + 1165);
    const auto *nsh0_1167 = buffer.data(nsh0 + 1167);
    const auto *nsh0_1168 = buffer.data(nsh0 + 1168);
    const auto *nsh0_1169 = buffer.data(nsh0 + 1169);
    const auto *nsh0_1170 = buffer.data(nsh0 + 1170);
    const auto *nsh0_1171 = buffer.data(nsh0 + 1171);
    const auto *nsh0_1172 = buffer.data(nsh0 + 1172);
    const auto *nsh0_1173 = buffer.data(nsh0 + 1173);
    const auto *nsh0_1174 = buffer.data(nsh0 + 1174);
    const auto *nsh0_1175 = buffer.data(nsh0 + 1175);
    const auto *nsh0_1178 = buffer.data(nsh0 + 1178);
    const auto *nsh0_1180 = buffer.data(nsh0 + 1180);
    const auto *nsh0_1181 = buffer.data(nsh0 + 1181);
    const auto *nsh0_1183 = buffer.data(nsh0 + 1183);
    const auto *nsh0_1184 = buffer.data(nsh0 + 1184);
    const auto *nsh0_1185 = buffer.data(nsh0 + 1185);
    const auto *nsh0_1187 = buffer.data(nsh0 + 1187);
    const auto *nsh0_1188 = buffer.data(nsh0 + 1188);
    const auto *nsh0_1189 = buffer.data(nsh0 + 1189);
    const auto *nsh0_1190 = buffer.data(nsh0 + 1190);
    const auto *nsh0_1192 = buffer.data(nsh0 + 1192);
    const auto *nsh0_1193 = buffer.data(nsh0 + 1193);
    const auto *nsh0_1194 = buffer.data(nsh0 + 1194);
    const auto *nsh0_1195 = buffer.data(nsh0 + 1195);
    const auto *nsh0_1196 = buffer.data(nsh0 + 1196);
    const auto *nsh0_1197 = buffer.data(nsh0 + 1197);
    const auto *nsh0_1198 = buffer.data(nsh0 + 1198);
    const auto *nsh0_1199 = buffer.data(nsh0 + 1199);
    const auto *nsh0_1200 = buffer.data(nsh0 + 1200);
    const auto *nsh0_1201 = buffer.data(nsh0 + 1201);
    const auto *nsh0_1202 = buffer.data(nsh0 + 1202);
    const auto *nsh0_1203 = buffer.data(nsh0 + 1203);
    const auto *nsh0_1204 = buffer.data(nsh0 + 1204);
    const auto *nsh0_1205 = buffer.data(nsh0 + 1205);
    const auto *nsh0_1206 = buffer.data(nsh0 + 1206);
    const auto *nsh0_1207 = buffer.data(nsh0 + 1207);
    const auto *nsh0_1208 = buffer.data(nsh0 + 1208);
    const auto *nsh0_1209 = buffer.data(nsh0 + 1209);
    const auto *nsh0_1210 = buffer.data(nsh0 + 1210);
    const auto *nsh0_1211 = buffer.data(nsh0 + 1211);

    const auto *nsh1_1134 = buffer.data(nsh1 + 1134);
    const auto *nsh1_1135 = buffer.data(nsh1 + 1135);
    const auto *nsh1_1136 = buffer.data(nsh1 + 1136);
    const auto *nsh1_1137 = buffer.data(nsh1 + 1137);
    const auto *nsh1_1138 = buffer.data(nsh1 + 1138);
    const auto *nsh1_1139 = buffer.data(nsh1 + 1139);
    const auto *nsh1_1140 = buffer.data(nsh1 + 1140);
    const auto *nsh1_1141 = buffer.data(nsh1 + 1141);
    const auto *nsh1_1142 = buffer.data(nsh1 + 1142);
    const auto *nsh1_1143 = buffer.data(nsh1 + 1143);
    const auto *nsh1_1155 = buffer.data(nsh1 + 1155);
    const auto *nsh1_1156 = buffer.data(nsh1 + 1156);
    const auto *nsh1_1158 = buffer.data(nsh1 + 1158);
    const auto *nsh1_1160 = buffer.data(nsh1 + 1160);
    const auto *nsh1_1161 = buffer.data(nsh1 + 1161);
    const auto *nsh1_1163 = buffer.data(nsh1 + 1163);
    const auto *nsh1_1164 = buffer.data(nsh1 + 1164);
    const auto *nsh1_1165 = buffer.data(nsh1 + 1165);
    const auto *nsh1_1167 = buffer.data(nsh1 + 1167);
    const auto *nsh1_1168 = buffer.data(nsh1 + 1168);
    const auto *nsh1_1169 = buffer.data(nsh1 + 1169);
    const auto *nsh1_1170 = buffer.data(nsh1 + 1170);
    const auto *nsh1_1171 = buffer.data(nsh1 + 1171);
    const auto *nsh1_1172 = buffer.data(nsh1 + 1172);
    const auto *nsh1_1173 = buffer.data(nsh1 + 1173);
    const auto *nsh1_1174 = buffer.data(nsh1 + 1174);
    const auto *nsh1_1175 = buffer.data(nsh1 + 1175);
    const auto *nsh1_1178 = buffer.data(nsh1 + 1178);
    const auto *nsh1_1180 = buffer.data(nsh1 + 1180);
    const auto *nsh1_1181 = buffer.data(nsh1 + 1181);
    const auto *nsh1_1183 = buffer.data(nsh1 + 1183);
    const auto *nsh1_1184 = buffer.data(nsh1 + 1184);
    const auto *nsh1_1185 = buffer.data(nsh1 + 1185);
    const auto *nsh1_1187 = buffer.data(nsh1 + 1187);
    const auto *nsh1_1188 = buffer.data(nsh1 + 1188);
    const auto *nsh1_1189 = buffer.data(nsh1 + 1189);
    const auto *nsh1_1190 = buffer.data(nsh1 + 1190);
    const auto *nsh1_1192 = buffer.data(nsh1 + 1192);
    const auto *nsh1_1193 = buffer.data(nsh1 + 1193);
    const auto *nsh1_1194 = buffer.data(nsh1 + 1194);
    const auto *nsh1_1195 = buffer.data(nsh1 + 1195);
    const auto *nsh1_1196 = buffer.data(nsh1 + 1196);
    const auto *nsh1_1197 = buffer.data(nsh1 + 1197);
    const auto *nsh1_1198 = buffer.data(nsh1 + 1198);
    const auto *nsh1_1199 = buffer.data(nsh1 + 1199);
    const auto *nsh1_1200 = buffer.data(nsh1 + 1200);
    const auto *nsh1_1201 = buffer.data(nsh1 + 1201);
    const auto *nsh1_1202 = buffer.data(nsh1 + 1202);
    const auto *nsh1_1203 = buffer.data(nsh1 + 1203);
    const auto *nsh1_1204 = buffer.data(nsh1 + 1204);
    const auto *nsh1_1205 = buffer.data(nsh1 + 1205);
    const auto *nsh1_1206 = buffer.data(nsh1 + 1206);
    const auto *nsh1_1207 = buffer.data(nsh1 + 1207);
    const auto *nsh1_1208 = buffer.data(nsh1 + 1208);
    const auto *nsh1_1209 = buffer.data(nsh1 + 1209);
    const auto *nsh1_1210 = buffer.data(nsh1 + 1210);
    const auto *nsh1_1211 = buffer.data(nsh1 + 1211);

    const auto *nsi_1512 = buffer.data(nsi + 1512);
    const auto *nsi_1513 = buffer.data(nsi + 1513);
    const auto *nsi_1514 = buffer.data(nsi + 1514);
    const auto *nsi_1515 = buffer.data(nsi + 1515);
    const auto *nsi_1516 = buffer.data(nsi + 1516);
    const auto *nsi_1517 = buffer.data(nsi + 1517);
    const auto *nsi_1518 = buffer.data(nsi + 1518);
    const auto *nsi_1519 = buffer.data(nsi + 1519);
    const auto *nsi_1520 = buffer.data(nsi + 1520);
    const auto *nsi_1521 = buffer.data(nsi + 1521);
    const auto *nsi_1522 = buffer.data(nsi + 1522);
    const auto *nsi_1523 = buffer.data(nsi + 1523);
    const auto *nsi_1524 = buffer.data(nsi + 1524);
    const auto *nsi_1525 = buffer.data(nsi + 1525);
    const auto *nsi_1526 = buffer.data(nsi + 1526);
    const auto *nsi_1532 = buffer.data(nsi + 1532);
    const auto *nsi_1533 = buffer.data(nsi + 1533);
    const auto *nsi_1534 = buffer.data(nsi + 1534);
    const auto *nsi_1535 = buffer.data(nsi + 1535);
    const auto *nsi_1536 = buffer.data(nsi + 1536);
    const auto *nsi_1537 = buffer.data(nsi + 1537);
    const auto *nsi_1539 = buffer.data(nsi + 1539);
    const auto *nsi_1540 = buffer.data(nsi + 1540);
    const auto *nsi_1541 = buffer.data(nsi + 1541);
    const auto *nsi_1543 = buffer.data(nsi + 1543);
    const auto *nsi_1545 = buffer.data(nsi + 1545);
    const auto *nsi_1546 = buffer.data(nsi + 1546);
    const auto *nsi_1548 = buffer.data(nsi + 1548);
    const auto *nsi_1549 = buffer.data(nsi + 1549);
    const auto *nsi_1550 = buffer.data(nsi + 1550);
    const auto *nsi_1552 = buffer.data(nsi + 1552);
    const auto *nsi_1553 = buffer.data(nsi + 1553);
    const auto *nsi_1554 = buffer.data(nsi + 1554);
    const auto *nsi_1555 = buffer.data(nsi + 1555);
    const auto *nsi_1557 = buffer.data(nsi + 1557);
    const auto *nsi_1558 = buffer.data(nsi + 1558);
    const auto *nsi_1559 = buffer.data(nsi + 1559);
    const auto *nsi_1560 = buffer.data(nsi + 1560);
    const auto *nsi_1561 = buffer.data(nsi + 1561);
    const auto *nsi_1562 = buffer.data(nsi + 1562);
    const auto *nsi_1563 = buffer.data(nsi + 1563);
    const auto *nsi_1564 = buffer.data(nsi + 1564);
    const auto *nsi_1565 = buffer.data(nsi + 1565);
    const auto *nsi_1566 = buffer.data(nsi + 1566);
    const auto *nsi_1567 = buffer.data(nsi + 1567);
    const auto *nsi_1570 = buffer.data(nsi + 1570);
    const auto *nsi_1572 = buffer.data(nsi + 1572);
    const auto *nsi_1573 = buffer.data(nsi + 1573);
    const auto *nsi_1575 = buffer.data(nsi + 1575);
    const auto *nsi_1576 = buffer.data(nsi + 1576);
    const auto *nsi_1577 = buffer.data(nsi + 1577);
    const auto *nsi_1579 = buffer.data(nsi + 1579);
    const auto *nsi_1580 = buffer.data(nsi + 1580);
    const auto *nsi_1581 = buffer.data(nsi + 1581);
    const auto *nsi_1582 = buffer.data(nsi + 1582);
    const auto *nsi_1584 = buffer.data(nsi + 1584);
    const auto *nsi_1585 = buffer.data(nsi + 1585);
    const auto *nsi_1586 = buffer.data(nsi + 1586);
    const auto *nsi_1587 = buffer.data(nsi + 1587);
    const auto *nsi_1588 = buffer.data(nsi + 1588);
    const auto *nsi_1589 = buffer.data(nsi + 1589);
    const auto *nsi_1590 = buffer.data(nsi + 1590);
    const auto *nsi_1591 = buffer.data(nsi + 1591);
    const auto *nsi_1592 = buffer.data(nsi + 1592);
    const auto *nsi_1593 = buffer.data(nsi + 1593);
    const auto *nsi_1594 = buffer.data(nsi + 1594);
    const auto *nsi_1595 = buffer.data(nsi + 1595);
    const auto *nsi_1596 = buffer.data(nsi + 1596);
    const auto *nsi_1597 = buffer.data(nsi + 1597);
    const auto *nsi_1598 = buffer.data(nsi + 1598);
    const auto *nsi_1599 = buffer.data(nsi + 1599);
    const auto *nsi_1600 = buffer.data(nsi + 1600);
    const auto *nsi_1601 = buffer.data(nsi + 1601);
    const auto *nsi_1602 = buffer.data(nsi + 1602);
    const auto *nsi_1603 = buffer.data(nsi + 1603);
    const auto *nsi_1604 = buffer.data(nsi + 1604);
    const auto *nsi_1605 = buffer.data(nsi + 1605);
    const auto *nsi_1606 = buffer.data(nsi + 1606);
    const auto *nsi_1607 = buffer.data(nsi + 1607);
    const auto *nsi_1608 = buffer.data(nsi + 1608);
    const auto *nsi_1609 = buffer.data(nsi + 1609);
    const auto *nsi_1610 = buffer.data(nsi + 1610);

#pragma omp simd aligned(t_1944, t_1945, t_1946, t_1947, pa_x, pc_x, pc_y, pc_z, msk0_1944, \
                         msi_1232, msi_1512, msk1_1944, nsh0_1134, nsh1_1134, nsi_1512, \
                         nsi_1513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1944[k] = pa_x[k] * msk0_1944[k]
                    + f_22 * msi_1512[k]
                    - f_12 * pc_x[k] * msk1_1944[k];

        t_1945[k] = f_3 * pc_y[k] * nsi_1512[k];

        t_1946[k] = f_18 * msi_1232[k]
                    + f_3 * pc_z[k] * nsi_1512[k];

        t_1947[k] = f_4 * nsh0_1134[k]
                    - f_5 * nsh1_1134[k]
                    + f_3 * pc_y[k] * nsi_1513[k];
    }

#pragma omp simd aligned(t_1948, t_1949, t_1950, pa_x, pc_x, pc_y, msk0_1949, msi_1517, \
                         msk1_1949, nsh0_1135, nsh1_1135, nsi_1514, \
                         nsi_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1948[k] = f_3 * pc_y[k] * nsi_1514[k];

        t_1949[k] = pa_x[k] * msk0_1949[k]
                    + f_17 * msi_1517[k]
                    - f_12 * pc_x[k] * msk1_1949[k];

        t_1950[k] = f_6 * nsh0_1135[k]
                    - f_7 * nsh1_1135[k]
                    + f_3 * pc_y[k] * nsi_1515[k];
    }

#pragma omp simd aligned(t_1951, t_1952, t_1953, pa_x, pc_x, pc_y, msk0_1953, msi_1521, \
                         msk1_1953, nsh0_1136, nsh1_1136, nsi_1516, \
                         nsi_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1951[k] = f_4 * nsh0_1136[k]
                    - f_5 * nsh1_1136[k]
                    + f_3 * pc_y[k] * nsi_1516[k];

        t_1952[k] = f_3 * pc_y[k] * nsi_1517[k];

        t_1953[k] = pa_x[k] * msk0_1953[k]
                    + f_16 * msi_1521[k]
                    - f_12 * pc_x[k] * msk1_1953[k];
    }

#pragma omp simd aligned(t_1954, t_1955, t_1956, t_1957, pc_y, nsh0_1137, nsh0_1138, \
                         nsh0_1139, nsh1_1137, nsh1_1138, nsh1_1139, nsi_1518, nsi_1519, \
                         nsi_1520, nsi_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1954[k] = f_8 * nsh0_1137[k]
                    - f_9 * nsh1_1137[k]
                    + f_3 * pc_y[k] * nsi_1518[k];

        t_1955[k] = f_6 * nsh0_1138[k]
                    - f_7 * nsh1_1138[k]
                    + f_3 * pc_y[k] * nsi_1519[k];

        t_1956[k] = f_4 * nsh0_1139[k]
                    - f_5 * nsh1_1139[k]
                    + f_3 * pc_y[k] * nsi_1520[k];

        t_1957[k] = f_3 * pc_y[k] * nsi_1521[k];
    }

#pragma omp simd aligned(t_1958, t_1959, t_1960, pa_x, pc_x, pc_y, msk0_1958, msi_1526, \
                         msk1_1958, nsh0_1140, nsh0_1141, nsh1_1140, nsh1_1141, nsi_1522, \
                         nsi_1523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1958[k] = pa_x[k] * msk0_1958[k]
                    + f_15 * msi_1526[k]
                    - f_12 * pc_x[k] * msk1_1958[k];

        t_1959[k] = f_10 * nsh0_1140[k]
                    - f_11 * nsh1_1140[k]
                    + f_3 * pc_y[k] * nsi_1522[k];

        t_1960[k] = f_8 * nsh0_1141[k]
                    - f_9 * nsh1_1141[k]
                    + f_3 * pc_y[k] * nsi_1523[k];
    }

#pragma omp simd aligned(t_1961, t_1962, t_1963, pc_y, nsh0_1142, nsh0_1143, nsh1_1142, \
                         nsh1_1143, nsi_1524, nsi_1525, nsi_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1961[k] = f_6 * nsh0_1142[k]
                    - f_7 * nsh1_1142[k]
                    + f_3 * pc_y[k] * nsi_1524[k];

        t_1962[k] = f_4 * nsh0_1143[k]
                    - f_5 * nsh1_1143[k]
                    + f_3 * pc_y[k] * nsi_1525[k];

        t_1963[k] = f_3 * pc_y[k] * nsi_1526[k];
    }

#pragma omp simd aligned(t_1964, t_1965, t_1966, t_1967, pa_x, pc_x, msk0_1964, msi_1532, \
                         msi_1533, msi_1534, msi_1535, msk1_1964, nsi_1533, nsi_1534, \
                         nsi_1535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1964[k] = pa_x[k] * msk0_1964[k]
                    + f_14 * msi_1532[k]
                    - f_12 * pc_x[k] * msk1_1964[k];

        t_1965[k] = f_13 * msi_1533[k]
                    + f_3 * pc_x[k] * nsi_1533[k];

        t_1966[k] = f_13 * msi_1534[k]
                    + f_3 * pc_x[k] * nsi_1534[k];

        t_1967[k] = f_13 * msi_1535[k]
                    + f_3 * pc_x[k] * nsi_1535[k];
    }

#pragma omp simd aligned(t_1968, t_1969, t_1970, t_1971, pc_x, pc_y, msi_1536, msi_1537, \
                         msi_1539, nsi_1532, nsi_1536, nsi_1537, \
                         nsi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1968[k] = f_13 * msi_1536[k]
                    + f_3 * pc_x[k] * nsi_1536[k];

        t_1969[k] = f_13 * msi_1537[k]
                    + f_3 * pc_x[k] * nsi_1537[k];

        t_1970[k] = f_3 * pc_y[k] * nsi_1532[k];

        t_1971[k] = f_13 * msi_1539[k]
                    + f_3 * pc_x[k] * nsi_1539[k];
    }

#pragma omp simd aligned(t_1972, t_1973, t_1974, t_1975, pa_x, pc_x, msk0_1972, msk0_1973, \
                         msk0_1974, msk0_1975, msk1_1972, msk1_1973, msk1_1974, \
                         msk1_1975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1972[k] = pa_x[k] * msk0_1972[k]
                    - f_12 * pc_x[k] * msk1_1972[k];

        t_1973[k] = pa_x[k] * msk0_1973[k]
                    - f_12 * pc_x[k] * msk1_1973[k];

        t_1974[k] = pa_x[k] * msk0_1974[k]
                    - f_12 * pc_x[k] * msk1_1974[k];

        t_1975[k] = pa_x[k] * msk0_1975[k]
                    - f_12 * pc_x[k] * msk1_1975[k];
    }

#pragma omp simd aligned(t_1976, t_1977, t_1978, t_1979, pa_x, pc_x, pc_y, msk0_1976, \
                         msk0_1977, msk0_1979, msk1_1976, msk1_1977, msk1_1979, \
                         nsi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1976[k] = pa_x[k] * msk0_1976[k]
                    - f_12 * pc_x[k] * msk1_1976[k];

        t_1977[k] = pa_x[k] * msk0_1977[k]
                    - f_12 * pc_x[k] * msk1_1977[k];

        t_1978[k] = f_3 * pc_y[k] * nsi_1539[k];

        t_1979[k] = pa_x[k] * msk0_1979[k]
                    - f_12 * pc_x[k] * msk1_1979[k];
    }

#pragma omp simd aligned(t_1980, t_1981, t_1982, t_1983, t_1984, pc_x, pc_z, nsh0_1155, \
                         nsh0_1156, nsh0_1158, nsh1_1155, nsh1_1156, nsh1_1158, nsi_1540, \
                         nsi_1541, nsi_1543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1980[k] = f_1 * nsh0_1155[k]
                    - f_2 * nsh1_1155[k]
                    + f_3 * pc_x[k] * nsi_1540[k];

        t_1981[k] = f_19 * nsh0_1156[k]
                    - f_20 * nsh1_1156[k]
                    + f_3 * pc_x[k] * nsi_1541[k];

        t_1982[k] = f_3 * pc_z[k] * nsi_1540[k];

        t_1983[k] = f_10 * nsh0_1158[k]
                    - f_11 * nsh1_1158[k]
                    + f_3 * pc_x[k] * nsi_1543[k];

        t_1984[k] = f_3 * pc_z[k] * nsi_1541[k];
    }

#pragma omp simd aligned(t_1985, t_1986, t_1987, t_1988, pc_x, pc_z, nsh0_1160, nsh0_1161, \
                         nsh0_1163, nsh1_1160, nsh1_1161, nsh1_1163, nsi_1543, nsi_1545, \
                         nsi_1546, nsi_1548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1985[k] = f_10 * nsh0_1160[k]
                    - f_11 * nsh1_1160[k]
                    + f_3 * pc_x[k] * nsi_1545[k];

        t_1986[k] = f_8 * nsh0_1161[k]
                    - f_9 * nsh1_1161[k]
                    + f_3 * pc_x[k] * nsi_1546[k];

        t_1987[k] = f_3 * pc_z[k] * nsi_1543[k];

        t_1988[k] = f_8 * nsh0_1163[k]
                    - f_9 * nsh1_1163[k]
                    + f_3 * pc_x[k] * nsi_1548[k];
    }

#pragma omp simd aligned(t_1989, t_1990, t_1991, t_1992, pc_x, pc_z, nsh0_1164, nsh0_1165, \
                         nsh0_1167, nsh1_1164, nsh1_1165, nsh1_1167, nsi_1546, nsi_1549, \
                         nsi_1550, nsi_1552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1989[k] = f_8 * nsh0_1164[k]
                    - f_9 * nsh1_1164[k]
                    + f_3 * pc_x[k] * nsi_1549[k];

        t_1990[k] = f_6 * nsh0_1165[k]
                    - f_7 * nsh1_1165[k]
                    + f_3 * pc_x[k] * nsi_1550[k];

        t_1991[k] = f_3 * pc_z[k] * nsi_1546[k];

        t_1992[k] = f_6 * nsh0_1167[k]
                    - f_7 * nsh1_1167[k]
                    + f_3 * pc_x[k] * nsi_1552[k];
    }

#pragma omp simd aligned(t_1993, t_1994, t_1995, t_1996, pc_x, pc_z, nsh0_1168, nsh0_1169, \
                         nsh0_1170, nsh1_1168, nsh1_1169, nsh1_1170, nsi_1550, nsi_1553, \
                         nsi_1554, nsi_1555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1993[k] = f_6 * nsh0_1168[k]
                    - f_7 * nsh1_1168[k]
                    + f_3 * pc_x[k] * nsi_1553[k];

        t_1994[k] = f_6 * nsh0_1169[k]
                    - f_7 * nsh1_1169[k]
                    + f_3 * pc_x[k] * nsi_1554[k];

        t_1995[k] = f_4 * nsh0_1170[k]
                    - f_5 * nsh1_1170[k]
                    + f_3 * pc_x[k] * nsi_1555[k];

        t_1996[k] = f_3 * pc_z[k] * nsi_1550[k];
    }

#pragma omp simd aligned(t_1997, t_1998, t_1999, pc_x, nsh0_1172, nsh0_1173, nsh0_1174, \
                         nsh1_1172, nsh1_1173, nsh1_1174, nsi_1557, nsi_1558, \
                         nsi_1559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1997[k] = f_4 * nsh0_1172[k]
                    - f_5 * nsh1_1172[k]
                    + f_3 * pc_x[k] * nsi_1557[k];

        t_1998[k] = f_4 * nsh0_1173[k]
                    - f_5 * nsh1_1173[k]
                    + f_3 * pc_x[k] * nsi_1558[k];

        t_1999[k] = f_4 * nsh0_1174[k]
                    - f_5 * nsh1_1174[k]
                    + f_3 * pc_x[k] * nsi_1559[k];
    }

#pragma omp simd aligned(t_2000, t_2001, t_2002, t_2003, t_2004, t_2005, pc_x, nsh0_1175, \
                         nsh1_1175, nsi_1560, nsi_1561, nsi_1562, nsi_1563, nsi_1564, \
                         nsi_1565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2000[k] = f_4 * nsh0_1175[k]
                    - f_5 * nsh1_1175[k]
                    + f_3 * pc_x[k] * nsi_1560[k];

        t_2001[k] = f_3 * pc_x[k] * nsi_1561[k];

        t_2002[k] = f_3 * pc_x[k] * nsi_1562[k];

        t_2003[k] = f_3 * pc_x[k] * nsi_1563[k];

        t_2004[k] = f_3 * pc_x[k] * nsi_1564[k];

        t_2005[k] = f_3 * pc_x[k] * nsi_1565[k];
    }

#pragma omp simd aligned(t_2006, t_2007, t_2008, t_2009, t_2010, pc_x, pc_y, pc_z, msi_1281, \
                         nsh0_1170, nsh1_1170, nsi_1561, nsi_1562, nsi_1566, \
                         nsi_1567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2006[k] = f_3 * pc_x[k] * nsi_1566[k];

        t_2007[k] = f_3 * pc_x[k] * nsi_1567[k];

        t_2008[k] = f_0 * msi_1281[k]
                    + f_1 * nsh0_1170[k]
                    - f_2 * nsh1_1170[k]
                    + f_3 * pc_y[k] * nsi_1561[k];

        t_2009[k] = f_3 * pc_z[k] * nsi_1561[k];

        t_2010[k] = f_4 * nsh0_1170[k]
                    - f_5 * nsh1_1170[k]
                    + f_3 * pc_z[k] * nsi_1562[k];
    }

#pragma omp simd aligned(t_2011, t_2012, t_2013, pc_z, nsh0_1171, nsh0_1172, nsh0_1173, \
                         nsh1_1171, nsh1_1172, nsh1_1173, nsi_1563, nsi_1564, \
                         nsi_1565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2011[k] = f_6 * nsh0_1171[k]
                    - f_7 * nsh1_1171[k]
                    + f_3 * pc_z[k] * nsi_1563[k];

        t_2012[k] = f_8 * nsh0_1172[k]
                    - f_9 * nsh1_1172[k]
                    + f_3 * pc_z[k] * nsi_1564[k];

        t_2013[k] = f_10 * nsh0_1173[k]
                    - f_11 * nsh1_1173[k]
                    + f_3 * pc_z[k] * nsi_1565[k];
    }

#pragma omp simd aligned(t_2014, t_2015, t_2016, t_2017, pa_z, pc_y, pc_z, msk0_1620, \
                         msk0_1621, msi_1287, msk1_1620, msk1_1621, nsh0_1175, nsh1_1175, \
                         nsi_1567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2014[k] = f_0 * msi_1287[k]
                    + f_3 * pc_y[k] * nsi_1567[k];

        t_2015[k] = f_1 * nsh0_1175[k]
                    - f_2 * nsh1_1175[k]
                    + f_3 * pc_z[k] * nsi_1567[k];

        t_2016[k] = pa_z[k] * msk0_1620[k]
                    - f_12 * pc_z[k] * msk1_1620[k];

        t_2017[k] = pa_z[k] * msk0_1621[k]
                    - f_12 * pc_z[k] * msk1_1621[k];
    }

#pragma omp simd aligned(t_2018, t_2019, t_2020, pa_z, pc_x, pc_z, msk0_1623, msk1_1623, \
                         nsh0_1178, nsh0_1180, nsh1_1178, nsh1_1180, nsi_1570, \
                         nsi_1572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2018[k] = f_19 * nsh0_1178[k]
                    - f_20 * nsh1_1178[k]
                    + f_3 * pc_x[k] * nsi_1570[k];

        t_2019[k] = pa_z[k] * msk0_1623[k]
                    - f_12 * pc_z[k] * msk1_1623[k];

        t_2020[k] = f_10 * nsh0_1180[k]
                    - f_11 * nsh1_1180[k]
                    + f_3 * pc_x[k] * nsi_1572[k];
    }

#pragma omp simd aligned(t_2021, t_2022, t_2023, pa_z, pc_x, pc_z, msk0_1626, msk1_1626, \
                         nsh0_1181, nsh0_1183, nsh1_1181, nsh1_1183, nsi_1573, \
                         nsi_1575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2021[k] = f_10 * nsh0_1181[k]
                    - f_11 * nsh1_1181[k]
                    + f_3 * pc_x[k] * nsi_1573[k];

        t_2022[k] = pa_z[k] * msk0_1626[k]
                    - f_12 * pc_z[k] * msk1_1626[k];

        t_2023[k] = f_8 * nsh0_1183[k]
                    - f_9 * nsh1_1183[k]
                    + f_3 * pc_x[k] * nsi_1575[k];
    }

#pragma omp simd aligned(t_2024, t_2025, t_2026, pa_z, pc_x, pc_z, msk0_1630, msk1_1630, \
                         nsh0_1184, nsh0_1185, nsh1_1184, nsh1_1185, nsi_1576, \
                         nsi_1577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2024[k] = f_8 * nsh0_1184[k]
                    - f_9 * nsh1_1184[k]
                    + f_3 * pc_x[k] * nsi_1576[k];

        t_2025[k] = f_8 * nsh0_1185[k]
                    - f_9 * nsh1_1185[k]
                    + f_3 * pc_x[k] * nsi_1577[k];

        t_2026[k] = pa_z[k] * msk0_1630[k]
                    - f_12 * pc_z[k] * msk1_1630[k];
    }

#pragma omp simd aligned(t_2027, t_2028, t_2029, pc_x, nsh0_1187, nsh0_1188, nsh0_1189, \
                         nsh1_1187, nsh1_1188, nsh1_1189, nsi_1579, nsi_1580, \
                         nsi_1581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2027[k] = f_6 * nsh0_1187[k]
                    - f_7 * nsh1_1187[k]
                    + f_3 * pc_x[k] * nsi_1579[k];

        t_2028[k] = f_6 * nsh0_1188[k]
                    - f_7 * nsh1_1188[k]
                    + f_3 * pc_x[k] * nsi_1580[k];

        t_2029[k] = f_6 * nsh0_1189[k]
                    - f_7 * nsh1_1189[k]
                    + f_3 * pc_x[k] * nsi_1581[k];
    }

#pragma omp simd aligned(t_2030, t_2031, t_2032, pa_z, pc_x, pc_z, msk0_1635, msk1_1635, \
                         nsh0_1190, nsh0_1192, nsh1_1190, nsh1_1192, nsi_1582, \
                         nsi_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2030[k] = f_6 * nsh0_1190[k]
                    - f_7 * nsh1_1190[k]
                    + f_3 * pc_x[k] * nsi_1582[k];

        t_2031[k] = pa_z[k] * msk0_1635[k]
                    - f_12 * pc_z[k] * msk1_1635[k];

        t_2032[k] = f_4 * nsh0_1192[k]
                    - f_5 * nsh1_1192[k]
                    + f_3 * pc_x[k] * nsi_1584[k];
    }

#pragma omp simd aligned(t_2033, t_2034, t_2035, pc_x, nsh0_1193, nsh0_1194, nsh0_1195, \
                         nsh1_1193, nsh1_1194, nsh1_1195, nsi_1585, nsi_1586, \
                         nsi_1587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2033[k] = f_4 * nsh0_1193[k]
                    - f_5 * nsh1_1193[k]
                    + f_3 * pc_x[k] * nsi_1585[k];

        t_2034[k] = f_4 * nsh0_1194[k]
                    - f_5 * nsh1_1194[k]
                    + f_3 * pc_x[k] * nsi_1586[k];

        t_2035[k] = f_4 * nsh0_1195[k]
                    - f_5 * nsh1_1195[k]
                    + f_3 * pc_x[k] * nsi_1587[k];
    }

#pragma omp simd aligned(t_2036, t_2037, t_2038, t_2039, t_2040, t_2041, pc_x, nsh0_1196, \
                         nsh1_1196, nsi_1588, nsi_1589, nsi_1590, nsi_1591, nsi_1592, \
                         nsi_1593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2036[k] = f_4 * nsh0_1196[k]
                    - f_5 * nsh1_1196[k]
                    + f_3 * pc_x[k] * nsi_1588[k];

        t_2037[k] = f_3 * pc_x[k] * nsi_1589[k];

        t_2038[k] = f_3 * pc_x[k] * nsi_1590[k];

        t_2039[k] = f_3 * pc_x[k] * nsi_1591[k];

        t_2040[k] = f_3 * pc_x[k] * nsi_1592[k];

        t_2041[k] = f_3 * pc_x[k] * nsi_1593[k];
    }

#pragma omp simd aligned(t_2042, t_2043, t_2044, t_2045, pa_z, pc_x, pc_z, msk0_1648, \
                         msi_1281, msk1_1648, nsi_1589, nsi_1594, \
                         nsi_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2042[k] = f_3 * pc_x[k] * nsi_1594[k];

        t_2043[k] = f_3 * pc_x[k] * nsi_1595[k];

        t_2044[k] = pa_z[k] * msk0_1648[k]
                    - f_12 * pc_z[k] * msk1_1648[k];

        t_2045[k] = f_13 * msi_1281[k]
                    + f_3 * pc_z[k] * nsi_1589[k];
    }

#pragma omp simd aligned(t_2046, t_2047, t_2048, pa_z, pc_z, msk0_1650, msk0_1651, msk0_1652, \
                         msi_1282, msi_1283, msi_1284, msk1_1650, msk1_1651, \
                         msk1_1652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2046[k] = pa_z[k] * msk0_1650[k]
                    + f_14 * msi_1282[k]
                    - f_12 * pc_z[k] * msk1_1650[k];

        t_2047[k] = pa_z[k] * msk0_1651[k]
                    + f_15 * msi_1283[k]
                    - f_12 * pc_z[k] * msk1_1651[k];

        t_2048[k] = pa_z[k] * msk0_1652[k]
                    + f_16 * msi_1284[k]
                    - f_12 * pc_z[k] * msk1_1652[k];
    }

#pragma omp simd aligned(t_2049, t_2050, t_2051, pa_z, pc_y, pc_z, msk0_1653, msi_1285, \
                         msi_1287, msi_1315, msk1_1653, nsh0_1196, nsh1_1196, \
                         nsi_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2049[k] = pa_z[k] * msk0_1653[k]
                    + f_17 * msi_1285[k]
                    - f_12 * pc_z[k] * msk1_1653[k];

        t_2050[k] = f_18 * msi_1315[k]
                    + f_3 * pc_y[k] * nsi_1595[k];

        t_2051[k] = f_13 * msi_1287[k]
                    + f_1 * nsh0_1196[k]
                    - f_2 * nsh1_1196[k]
                    + f_3 * pc_z[k] * nsi_1595[k];
    }

#pragma omp simd aligned(t_2052, t_2053, t_2054, pc_x, nsh0_1197, nsh0_1198, nsh0_1199, \
                         nsh1_1197, nsh1_1198, nsh1_1199, nsi_1596, nsi_1597, \
                         nsi_1598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2052[k] = f_1 * nsh0_1197[k]
                    - f_2 * nsh1_1197[k]
                    + f_3 * pc_x[k] * nsi_1596[k];

        t_2053[k] = f_19 * nsh0_1198[k]
                    - f_20 * nsh1_1198[k]
                    + f_3 * pc_x[k] * nsi_1597[k];

        t_2054[k] = f_19 * nsh0_1199[k]
                    - f_20 * nsh1_1199[k]
                    + f_3 * pc_x[k] * nsi_1598[k];
    }

#pragma omp simd aligned(t_2055, t_2056, t_2057, pc_x, nsh0_1200, nsh0_1201, nsh0_1202, \
                         nsh1_1200, nsh1_1201, nsh1_1202, nsi_1599, nsi_1600, \
                         nsi_1601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2055[k] = f_10 * nsh0_1200[k]
                    - f_11 * nsh1_1200[k]
                    + f_3 * pc_x[k] * nsi_1599[k];

        t_2056[k] = f_10 * nsh0_1201[k]
                    - f_11 * nsh1_1201[k]
                    + f_3 * pc_x[k] * nsi_1600[k];

        t_2057[k] = f_10 * nsh0_1202[k]
                    - f_11 * nsh1_1202[k]
                    + f_3 * pc_x[k] * nsi_1601[k];
    }

#pragma omp simd aligned(t_2058, t_2059, t_2060, pc_x, nsh0_1203, nsh0_1204, nsh0_1205, \
                         nsh1_1203, nsh1_1204, nsh1_1205, nsi_1602, nsi_1603, \
                         nsi_1604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2058[k] = f_8 * nsh0_1203[k]
                    - f_9 * nsh1_1203[k]
                    + f_3 * pc_x[k] * nsi_1602[k];

        t_2059[k] = f_8 * nsh0_1204[k]
                    - f_9 * nsh1_1204[k]
                    + f_3 * pc_x[k] * nsi_1603[k];

        t_2060[k] = f_8 * nsh0_1205[k]
                    - f_9 * nsh1_1205[k]
                    + f_3 * pc_x[k] * nsi_1604[k];
    }

#pragma omp simd aligned(t_2061, t_2062, t_2063, pc_x, nsh0_1206, nsh0_1207, nsh0_1208, \
                         nsh1_1206, nsh1_1207, nsh1_1208, nsi_1605, nsi_1606, \
                         nsi_1607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2061[k] = f_8 * nsh0_1206[k]
                    - f_9 * nsh1_1206[k]
                    + f_3 * pc_x[k] * nsi_1605[k];

        t_2062[k] = f_6 * nsh0_1207[k]
                    - f_7 * nsh1_1207[k]
                    + f_3 * pc_x[k] * nsi_1606[k];

        t_2063[k] = f_6 * nsh0_1208[k]
                    - f_7 * nsh1_1208[k]
                    + f_3 * pc_x[k] * nsi_1607[k];
    }

#pragma omp simd aligned(t_2064, t_2065, t_2066, pc_x, nsh0_1209, nsh0_1210, nsh0_1211, \
                         nsh1_1209, nsh1_1210, nsh1_1211, nsi_1608, nsi_1609, \
                         nsi_1610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2064[k] = f_6 * nsh0_1209[k]
                    - f_7 * nsh1_1209[k]
                    + f_3 * pc_x[k] * nsi_1608[k];

        t_2065[k] = f_6 * nsh0_1210[k]
                    - f_7 * nsh1_1210[k]
                    + f_3 * pc_x[k] * nsi_1609[k];

        t_2066[k] = f_6 * nsh0_1211[k]
                    - f_7 * nsh1_1211[k]
                    + f_3 * pc_x[k] * nsi_1610[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t msi, const size_t nsh0,
                                                           const size_t nsh1, const size_t nsi,
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
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi_1309 = buffer.data(msi + 1309);
    const auto *msi_1315 = buffer.data(msi + 1315);
    const auto *msi_1337 = buffer.data(msi + 1337);
    const auto *msi_1339 = buffer.data(msi + 1339);
    const auto *msi_1340 = buffer.data(msi + 1340);
    const auto *msi_1341 = buffer.data(msi + 1341);
    const auto *msi_1342 = buffer.data(msi + 1342);
    const auto *msi_1343 = buffer.data(msi + 1343);
    const auto *msi_1365 = buffer.data(msi + 1365);
    const auto *msi_1367 = buffer.data(msi + 1367);
    const auto *msi_1368 = buffer.data(msi + 1368);
    const auto *msi_1369 = buffer.data(msi + 1369);
    const auto *msi_1370 = buffer.data(msi + 1370);
    const auto *msi_1371 = buffer.data(msi + 1371);
    const auto *msi_1393 = buffer.data(msi + 1393);
    const auto *msi_1395 = buffer.data(msi + 1395);
    const auto *msi_1396 = buffer.data(msi + 1396);
    const auto *msi_1397 = buffer.data(msi + 1397);
    const auto *msi_1398 = buffer.data(msi + 1398);
    const auto *msi_1399 = buffer.data(msi + 1399);

    const auto *nsh0_1212 = buffer.data(nsh0 + 1212);
    const auto *nsh0_1213 = buffer.data(nsh0 + 1213);
    const auto *nsh0_1214 = buffer.data(nsh0 + 1214);
    const auto *nsh0_1215 = buffer.data(nsh0 + 1215);
    const auto *nsh0_1216 = buffer.data(nsh0 + 1216);
    const auto *nsh0_1217 = buffer.data(nsh0 + 1217);
    const auto *nsh0_1218 = buffer.data(nsh0 + 1218);
    const auto *nsh0_1219 = buffer.data(nsh0 + 1219);
    const auto *nsh0_1220 = buffer.data(nsh0 + 1220);
    const auto *nsh0_1221 = buffer.data(nsh0 + 1221);
    const auto *nsh0_1222 = buffer.data(nsh0 + 1222);
    const auto *nsh0_1223 = buffer.data(nsh0 + 1223);
    const auto *nsh0_1224 = buffer.data(nsh0 + 1224);
    const auto *nsh0_1225 = buffer.data(nsh0 + 1225);
    const auto *nsh0_1226 = buffer.data(nsh0 + 1226);
    const auto *nsh0_1227 = buffer.data(nsh0 + 1227);
    const auto *nsh0_1228 = buffer.data(nsh0 + 1228);
    const auto *nsh0_1229 = buffer.data(nsh0 + 1229);
    const auto *nsh0_1230 = buffer.data(nsh0 + 1230);
    const auto *nsh0_1231 = buffer.data(nsh0 + 1231);
    const auto *nsh0_1232 = buffer.data(nsh0 + 1232);
    const auto *nsh0_1233 = buffer.data(nsh0 + 1233);
    const auto *nsh0_1234 = buffer.data(nsh0 + 1234);
    const auto *nsh0_1235 = buffer.data(nsh0 + 1235);
    const auto *nsh0_1236 = buffer.data(nsh0 + 1236);
    const auto *nsh0_1237 = buffer.data(nsh0 + 1237);
    const auto *nsh0_1238 = buffer.data(nsh0 + 1238);
    const auto *nsh0_1239 = buffer.data(nsh0 + 1239);
    const auto *nsh0_1240 = buffer.data(nsh0 + 1240);
    const auto *nsh0_1241 = buffer.data(nsh0 + 1241);
    const auto *nsh0_1242 = buffer.data(nsh0 + 1242);
    const auto *nsh0_1243 = buffer.data(nsh0 + 1243);
    const auto *nsh0_1244 = buffer.data(nsh0 + 1244);
    const auto *nsh0_1245 = buffer.data(nsh0 + 1245);
    const auto *nsh0_1246 = buffer.data(nsh0 + 1246);
    const auto *nsh0_1247 = buffer.data(nsh0 + 1247);
    const auto *nsh0_1248 = buffer.data(nsh0 + 1248);
    const auto *nsh0_1249 = buffer.data(nsh0 + 1249);
    const auto *nsh0_1250 = buffer.data(nsh0 + 1250);
    const auto *nsh0_1251 = buffer.data(nsh0 + 1251);
    const auto *nsh0_1252 = buffer.data(nsh0 + 1252);
    const auto *nsh0_1253 = buffer.data(nsh0 + 1253);
    const auto *nsh0_1254 = buffer.data(nsh0 + 1254);
    const auto *nsh0_1255 = buffer.data(nsh0 + 1255);
    const auto *nsh0_1256 = buffer.data(nsh0 + 1256);
    const auto *nsh0_1257 = buffer.data(nsh0 + 1257);
    const auto *nsh0_1258 = buffer.data(nsh0 + 1258);
    const auto *nsh0_1259 = buffer.data(nsh0 + 1259);
    const auto *nsh0_1260 = buffer.data(nsh0 + 1260);
    const auto *nsh0_1261 = buffer.data(nsh0 + 1261);
    const auto *nsh0_1262 = buffer.data(nsh0 + 1262);
    const auto *nsh0_1263 = buffer.data(nsh0 + 1263);
    const auto *nsh0_1264 = buffer.data(nsh0 + 1264);
    const auto *nsh0_1265 = buffer.data(nsh0 + 1265);
    const auto *nsh0_1266 = buffer.data(nsh0 + 1266);
    const auto *nsh0_1267 = buffer.data(nsh0 + 1267);
    const auto *nsh0_1268 = buffer.data(nsh0 + 1268);
    const auto *nsh0_1269 = buffer.data(nsh0 + 1269);
    const auto *nsh0_1270 = buffer.data(nsh0 + 1270);
    const auto *nsh0_1271 = buffer.data(nsh0 + 1271);
    const auto *nsh0_1272 = buffer.data(nsh0 + 1272);
    const auto *nsh0_1273 = buffer.data(nsh0 + 1273);
    const auto *nsh0_1274 = buffer.data(nsh0 + 1274);
    const auto *nsh0_1275 = buffer.data(nsh0 + 1275);
    const auto *nsh0_1276 = buffer.data(nsh0 + 1276);
    const auto *nsh0_1277 = buffer.data(nsh0 + 1277);
    const auto *nsh0_1278 = buffer.data(nsh0 + 1278);
    const auto *nsh0_1279 = buffer.data(nsh0 + 1279);
    const auto *nsh0_1280 = buffer.data(nsh0 + 1280);

    const auto *nsh1_1212 = buffer.data(nsh1 + 1212);
    const auto *nsh1_1213 = buffer.data(nsh1 + 1213);
    const auto *nsh1_1214 = buffer.data(nsh1 + 1214);
    const auto *nsh1_1215 = buffer.data(nsh1 + 1215);
    const auto *nsh1_1216 = buffer.data(nsh1 + 1216);
    const auto *nsh1_1217 = buffer.data(nsh1 + 1217);
    const auto *nsh1_1218 = buffer.data(nsh1 + 1218);
    const auto *nsh1_1219 = buffer.data(nsh1 + 1219);
    const auto *nsh1_1220 = buffer.data(nsh1 + 1220);
    const auto *nsh1_1221 = buffer.data(nsh1 + 1221);
    const auto *nsh1_1222 = buffer.data(nsh1 + 1222);
    const auto *nsh1_1223 = buffer.data(nsh1 + 1223);
    const auto *nsh1_1224 = buffer.data(nsh1 + 1224);
    const auto *nsh1_1225 = buffer.data(nsh1 + 1225);
    const auto *nsh1_1226 = buffer.data(nsh1 + 1226);
    const auto *nsh1_1227 = buffer.data(nsh1 + 1227);
    const auto *nsh1_1228 = buffer.data(nsh1 + 1228);
    const auto *nsh1_1229 = buffer.data(nsh1 + 1229);
    const auto *nsh1_1230 = buffer.data(nsh1 + 1230);
    const auto *nsh1_1231 = buffer.data(nsh1 + 1231);
    const auto *nsh1_1232 = buffer.data(nsh1 + 1232);
    const auto *nsh1_1233 = buffer.data(nsh1 + 1233);
    const auto *nsh1_1234 = buffer.data(nsh1 + 1234);
    const auto *nsh1_1235 = buffer.data(nsh1 + 1235);
    const auto *nsh1_1236 = buffer.data(nsh1 + 1236);
    const auto *nsh1_1237 = buffer.data(nsh1 + 1237);
    const auto *nsh1_1238 = buffer.data(nsh1 + 1238);
    const auto *nsh1_1239 = buffer.data(nsh1 + 1239);
    const auto *nsh1_1240 = buffer.data(nsh1 + 1240);
    const auto *nsh1_1241 = buffer.data(nsh1 + 1241);
    const auto *nsh1_1242 = buffer.data(nsh1 + 1242);
    const auto *nsh1_1243 = buffer.data(nsh1 + 1243);
    const auto *nsh1_1244 = buffer.data(nsh1 + 1244);
    const auto *nsh1_1245 = buffer.data(nsh1 + 1245);
    const auto *nsh1_1246 = buffer.data(nsh1 + 1246);
    const auto *nsh1_1247 = buffer.data(nsh1 + 1247);
    const auto *nsh1_1248 = buffer.data(nsh1 + 1248);
    const auto *nsh1_1249 = buffer.data(nsh1 + 1249);
    const auto *nsh1_1250 = buffer.data(nsh1 + 1250);
    const auto *nsh1_1251 = buffer.data(nsh1 + 1251);
    const auto *nsh1_1252 = buffer.data(nsh1 + 1252);
    const auto *nsh1_1253 = buffer.data(nsh1 + 1253);
    const auto *nsh1_1254 = buffer.data(nsh1 + 1254);
    const auto *nsh1_1255 = buffer.data(nsh1 + 1255);
    const auto *nsh1_1256 = buffer.data(nsh1 + 1256);
    const auto *nsh1_1257 = buffer.data(nsh1 + 1257);
    const auto *nsh1_1258 = buffer.data(nsh1 + 1258);
    const auto *nsh1_1259 = buffer.data(nsh1 + 1259);
    const auto *nsh1_1260 = buffer.data(nsh1 + 1260);
    const auto *nsh1_1261 = buffer.data(nsh1 + 1261);
    const auto *nsh1_1262 = buffer.data(nsh1 + 1262);
    const auto *nsh1_1263 = buffer.data(nsh1 + 1263);
    const auto *nsh1_1264 = buffer.data(nsh1 + 1264);
    const auto *nsh1_1265 = buffer.data(nsh1 + 1265);
    const auto *nsh1_1266 = buffer.data(nsh1 + 1266);
    const auto *nsh1_1267 = buffer.data(nsh1 + 1267);
    const auto *nsh1_1268 = buffer.data(nsh1 + 1268);
    const auto *nsh1_1269 = buffer.data(nsh1 + 1269);
    const auto *nsh1_1270 = buffer.data(nsh1 + 1270);
    const auto *nsh1_1271 = buffer.data(nsh1 + 1271);
    const auto *nsh1_1272 = buffer.data(nsh1 + 1272);
    const auto *nsh1_1273 = buffer.data(nsh1 + 1273);
    const auto *nsh1_1274 = buffer.data(nsh1 + 1274);
    const auto *nsh1_1275 = buffer.data(nsh1 + 1275);
    const auto *nsh1_1276 = buffer.data(nsh1 + 1276);
    const auto *nsh1_1277 = buffer.data(nsh1 + 1277);
    const auto *nsh1_1278 = buffer.data(nsh1 + 1278);
    const auto *nsh1_1279 = buffer.data(nsh1 + 1279);
    const auto *nsh1_1280 = buffer.data(nsh1 + 1280);

    const auto *nsi_1611 = buffer.data(nsi + 1611);
    const auto *nsi_1612 = buffer.data(nsi + 1612);
    const auto *nsi_1613 = buffer.data(nsi + 1613);
    const auto *nsi_1614 = buffer.data(nsi + 1614);
    const auto *nsi_1615 = buffer.data(nsi + 1615);
    const auto *nsi_1616 = buffer.data(nsi + 1616);
    const auto *nsi_1617 = buffer.data(nsi + 1617);
    const auto *nsi_1618 = buffer.data(nsi + 1618);
    const auto *nsi_1619 = buffer.data(nsi + 1619);
    const auto *nsi_1620 = buffer.data(nsi + 1620);
    const auto *nsi_1621 = buffer.data(nsi + 1621);
    const auto *nsi_1622 = buffer.data(nsi + 1622);
    const auto *nsi_1623 = buffer.data(nsi + 1623);
    const auto *nsi_1624 = buffer.data(nsi + 1624);
    const auto *nsi_1625 = buffer.data(nsi + 1625);
    const auto *nsi_1626 = buffer.data(nsi + 1626);
    const auto *nsi_1627 = buffer.data(nsi + 1627);
    const auto *nsi_1628 = buffer.data(nsi + 1628);
    const auto *nsi_1629 = buffer.data(nsi + 1629);
    const auto *nsi_1630 = buffer.data(nsi + 1630);
    const auto *nsi_1631 = buffer.data(nsi + 1631);
    const auto *nsi_1632 = buffer.data(nsi + 1632);
    const auto *nsi_1633 = buffer.data(nsi + 1633);
    const auto *nsi_1634 = buffer.data(nsi + 1634);
    const auto *nsi_1635 = buffer.data(nsi + 1635);
    const auto *nsi_1636 = buffer.data(nsi + 1636);
    const auto *nsi_1637 = buffer.data(nsi + 1637);
    const auto *nsi_1638 = buffer.data(nsi + 1638);
    const auto *nsi_1639 = buffer.data(nsi + 1639);
    const auto *nsi_1640 = buffer.data(nsi + 1640);
    const auto *nsi_1641 = buffer.data(nsi + 1641);
    const auto *nsi_1642 = buffer.data(nsi + 1642);
    const auto *nsi_1643 = buffer.data(nsi + 1643);
    const auto *nsi_1644 = buffer.data(nsi + 1644);
    const auto *nsi_1645 = buffer.data(nsi + 1645);
    const auto *nsi_1646 = buffer.data(nsi + 1646);
    const auto *nsi_1647 = buffer.data(nsi + 1647);
    const auto *nsi_1648 = buffer.data(nsi + 1648);
    const auto *nsi_1649 = buffer.data(nsi + 1649);
    const auto *nsi_1650 = buffer.data(nsi + 1650);
    const auto *nsi_1651 = buffer.data(nsi + 1651);
    const auto *nsi_1652 = buffer.data(nsi + 1652);
    const auto *nsi_1653 = buffer.data(nsi + 1653);
    const auto *nsi_1654 = buffer.data(nsi + 1654);
    const auto *nsi_1655 = buffer.data(nsi + 1655);
    const auto *nsi_1656 = buffer.data(nsi + 1656);
    const auto *nsi_1657 = buffer.data(nsi + 1657);
    const auto *nsi_1658 = buffer.data(nsi + 1658);
    const auto *nsi_1659 = buffer.data(nsi + 1659);
    const auto *nsi_1660 = buffer.data(nsi + 1660);
    const auto *nsi_1661 = buffer.data(nsi + 1661);
    const auto *nsi_1662 = buffer.data(nsi + 1662);
    const auto *nsi_1663 = buffer.data(nsi + 1663);
    const auto *nsi_1664 = buffer.data(nsi + 1664);
    const auto *nsi_1665 = buffer.data(nsi + 1665);
    const auto *nsi_1666 = buffer.data(nsi + 1666);
    const auto *nsi_1667 = buffer.data(nsi + 1667);
    const auto *nsi_1668 = buffer.data(nsi + 1668);
    const auto *nsi_1669 = buffer.data(nsi + 1669);
    const auto *nsi_1670 = buffer.data(nsi + 1670);
    const auto *nsi_1671 = buffer.data(nsi + 1671);
    const auto *nsi_1672 = buffer.data(nsi + 1672);
    const auto *nsi_1673 = buffer.data(nsi + 1673);
    const auto *nsi_1674 = buffer.data(nsi + 1674);
    const auto *nsi_1675 = buffer.data(nsi + 1675);
    const auto *nsi_1676 = buffer.data(nsi + 1676);
    const auto *nsi_1677 = buffer.data(nsi + 1677);
    const auto *nsi_1678 = buffer.data(nsi + 1678);
    const auto *nsi_1679 = buffer.data(nsi + 1679);
    const auto *nsi_1680 = buffer.data(nsi + 1680);
    const auto *nsi_1681 = buffer.data(nsi + 1681);
    const auto *nsi_1682 = buffer.data(nsi + 1682);
    const auto *nsi_1683 = buffer.data(nsi + 1683);
    const auto *nsi_1684 = buffer.data(nsi + 1684);
    const auto *nsi_1685 = buffer.data(nsi + 1685);
    const auto *nsi_1686 = buffer.data(nsi + 1686);
    const auto *nsi_1687 = buffer.data(nsi + 1687);
    const auto *nsi_1688 = buffer.data(nsi + 1688);
    const auto *nsi_1689 = buffer.data(nsi + 1689);
    const auto *nsi_1690 = buffer.data(nsi + 1690);
    const auto *nsi_1691 = buffer.data(nsi + 1691);
    const auto *nsi_1692 = buffer.data(nsi + 1692);
    const auto *nsi_1693 = buffer.data(nsi + 1693);
    const auto *nsi_1694 = buffer.data(nsi + 1694);
    const auto *nsi_1695 = buffer.data(nsi + 1695);
    const auto *nsi_1696 = buffer.data(nsi + 1696);
    const auto *nsi_1697 = buffer.data(nsi + 1697);
    const auto *nsi_1698 = buffer.data(nsi + 1698);
    const auto *nsi_1699 = buffer.data(nsi + 1699);
    const auto *nsi_1700 = buffer.data(nsi + 1700);
    const auto *nsi_1701 = buffer.data(nsi + 1701);

#pragma omp simd aligned(t_2067, t_2068, t_2069, pc_x, nsh0_1212, nsh0_1213, nsh0_1214, \
                         nsh1_1212, nsh1_1213, nsh1_1214, nsi_1611, nsi_1612, \
                         nsi_1613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2067[k] = f_4 * nsh0_1212[k]
                    - f_5 * nsh1_1212[k]
                    + f_3 * pc_x[k] * nsi_1611[k];

        t_2068[k] = f_4 * nsh0_1213[k]
                    - f_5 * nsh1_1213[k]
                    + f_3 * pc_x[k] * nsi_1612[k];

        t_2069[k] = f_4 * nsh0_1214[k]
                    - f_5 * nsh1_1214[k]
                    + f_3 * pc_x[k] * nsi_1613[k];
    }

#pragma omp simd aligned(t_2070, t_2071, t_2072, t_2073, pc_x, nsh0_1215, nsh0_1216, \
                         nsh0_1217, nsh1_1215, nsh1_1216, nsh1_1217, nsi_1614, nsi_1615, \
                         nsi_1616, nsi_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2070[k] = f_4 * nsh0_1215[k]
                    - f_5 * nsh1_1215[k]
                    + f_3 * pc_x[k] * nsi_1614[k];

        t_2071[k] = f_4 * nsh0_1216[k]
                    - f_5 * nsh1_1216[k]
                    + f_3 * pc_x[k] * nsi_1615[k];

        t_2072[k] = f_4 * nsh0_1217[k]
                    - f_5 * nsh1_1217[k]
                    + f_3 * pc_x[k] * nsi_1616[k];

        t_2073[k] = f_3 * pc_x[k] * nsi_1617[k];
    }

#pragma omp simd aligned(t_2074, t_2075, t_2076, t_2077, t_2078, t_2079, pc_x, nsi_1618, \
                         nsi_1619, nsi_1620, nsi_1621, nsi_1622, \
                         nsi_1623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2074[k] = f_3 * pc_x[k] * nsi_1618[k];

        t_2075[k] = f_3 * pc_x[k] * nsi_1619[k];

        t_2076[k] = f_3 * pc_x[k] * nsi_1620[k];

        t_2077[k] = f_3 * pc_x[k] * nsi_1621[k];

        t_2078[k] = f_3 * pc_x[k] * nsi_1622[k];

        t_2079[k] = f_3 * pc_x[k] * nsi_1623[k];
    }

#pragma omp simd aligned(t_2080, t_2081, t_2082, pc_y, pc_z, msi_1309, msi_1337, msi_1339, \
                         nsh0_1212, nsh0_1214, nsh1_1212, nsh1_1214, nsi_1617, \
                         nsi_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2080[k] = f_21 * msi_1337[k]
                    + f_1 * nsh0_1212[k]
                    - f_2 * nsh1_1212[k]
                    + f_3 * pc_y[k] * nsi_1617[k];

        t_2081[k] = f_14 * msi_1309[k]
                    + f_3 * pc_z[k] * nsi_1617[k];

        t_2082[k] = f_21 * msi_1339[k]
                    + f_10 * nsh0_1214[k]
                    - f_11 * nsh1_1214[k]
                    + f_3 * pc_y[k] * nsi_1619[k];
    }

#pragma omp simd aligned(t_2083, t_2084, t_2085, pc_y, msi_1340, msi_1341, msi_1342, \
                         nsh0_1215, nsh0_1216, nsh0_1217, nsh1_1215, nsh1_1216, nsh1_1217, \
                         nsi_1620, nsi_1621, nsi_1622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2083[k] = f_21 * msi_1340[k]
                    + f_8 * nsh0_1215[k]
                    - f_9 * nsh1_1215[k]
                    + f_3 * pc_y[k] * nsi_1620[k];

        t_2084[k] = f_21 * msi_1341[k]
                    + f_6 * nsh0_1216[k]
                    - f_7 * nsh1_1216[k]
                    + f_3 * pc_y[k] * nsi_1621[k];

        t_2085[k] = f_21 * msi_1342[k]
                    + f_4 * nsh0_1217[k]
                    - f_5 * nsh1_1217[k]
                    + f_3 * pc_y[k] * nsi_1622[k];
    }

#pragma omp simd aligned(t_2086, t_2087, t_2088, pc_x, pc_y, pc_z, msi_1315, msi_1343, \
                         nsh0_1217, nsh0_1218, nsh1_1217, nsh1_1218, nsi_1623, \
                         nsi_1624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2086[k] = f_21 * msi_1343[k]
                    + f_3 * pc_y[k] * nsi_1623[k];

        t_2087[k] = f_14 * msi_1315[k]
                    + f_1 * nsh0_1217[k]
                    - f_2 * nsh1_1217[k]
                    + f_3 * pc_z[k] * nsi_1623[k];

        t_2088[k] = f_1 * nsh0_1218[k]
                    - f_2 * nsh1_1218[k]
                    + f_3 * pc_x[k] * nsi_1624[k];
    }

#pragma omp simd aligned(t_2089, t_2090, t_2091, pc_x, nsh0_1219, nsh0_1220, nsh0_1221, \
                         nsh1_1219, nsh1_1220, nsh1_1221, nsi_1625, nsi_1626, \
                         nsi_1627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2089[k] = f_19 * nsh0_1219[k]
                    - f_20 * nsh1_1219[k]
                    + f_3 * pc_x[k] * nsi_1625[k];

        t_2090[k] = f_19 * nsh0_1220[k]
                    - f_20 * nsh1_1220[k]
                    + f_3 * pc_x[k] * nsi_1626[k];

        t_2091[k] = f_10 * nsh0_1221[k]
                    - f_11 * nsh1_1221[k]
                    + f_3 * pc_x[k] * nsi_1627[k];
    }

#pragma omp simd aligned(t_2092, t_2093, t_2094, pc_x, nsh0_1222, nsh0_1223, nsh0_1224, \
                         nsh1_1222, nsh1_1223, nsh1_1224, nsi_1628, nsi_1629, \
                         nsi_1630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2092[k] = f_10 * nsh0_1222[k]
                    - f_11 * nsh1_1222[k]
                    + f_3 * pc_x[k] * nsi_1628[k];

        t_2093[k] = f_10 * nsh0_1223[k]
                    - f_11 * nsh1_1223[k]
                    + f_3 * pc_x[k] * nsi_1629[k];

        t_2094[k] = f_8 * nsh0_1224[k]
                    - f_9 * nsh1_1224[k]
                    + f_3 * pc_x[k] * nsi_1630[k];
    }

#pragma omp simd aligned(t_2095, t_2096, t_2097, pc_x, nsh0_1225, nsh0_1226, nsh0_1227, \
                         nsh1_1225, nsh1_1226, nsh1_1227, nsi_1631, nsi_1632, \
                         nsi_1633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2095[k] = f_8 * nsh0_1225[k]
                    - f_9 * nsh1_1225[k]
                    + f_3 * pc_x[k] * nsi_1631[k];

        t_2096[k] = f_8 * nsh0_1226[k]
                    - f_9 * nsh1_1226[k]
                    + f_3 * pc_x[k] * nsi_1632[k];

        t_2097[k] = f_8 * nsh0_1227[k]
                    - f_9 * nsh1_1227[k]
                    + f_3 * pc_x[k] * nsi_1633[k];
    }

#pragma omp simd aligned(t_2098, t_2099, t_2100, pc_x, nsh0_1228, nsh0_1229, nsh0_1230, \
                         nsh1_1228, nsh1_1229, nsh1_1230, nsi_1634, nsi_1635, \
                         nsi_1636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2098[k] = f_6 * nsh0_1228[k]
                    - f_7 * nsh1_1228[k]
                    + f_3 * pc_x[k] * nsi_1634[k];

        t_2099[k] = f_6 * nsh0_1229[k]
                    - f_7 * nsh1_1229[k]
                    + f_3 * pc_x[k] * nsi_1635[k];

        t_2100[k] = f_6 * nsh0_1230[k]
                    - f_7 * nsh1_1230[k]
                    + f_3 * pc_x[k] * nsi_1636[k];
    }

#pragma omp simd aligned(t_2101, t_2102, t_2103, pc_x, nsh0_1231, nsh0_1232, nsh0_1233, \
                         nsh1_1231, nsh1_1232, nsh1_1233, nsi_1637, nsi_1638, \
                         nsi_1639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2101[k] = f_6 * nsh0_1231[k]
                    - f_7 * nsh1_1231[k]
                    + f_3 * pc_x[k] * nsi_1637[k];

        t_2102[k] = f_6 * nsh0_1232[k]
                    - f_7 * nsh1_1232[k]
                    + f_3 * pc_x[k] * nsi_1638[k];

        t_2103[k] = f_4 * nsh0_1233[k]
                    - f_5 * nsh1_1233[k]
                    + f_3 * pc_x[k] * nsi_1639[k];
    }

#pragma omp simd aligned(t_2104, t_2105, t_2106, pc_x, nsh0_1234, nsh0_1235, nsh0_1236, \
                         nsh1_1234, nsh1_1235, nsh1_1236, nsi_1640, nsi_1641, \
                         nsi_1642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2104[k] = f_4 * nsh0_1234[k]
                    - f_5 * nsh1_1234[k]
                    + f_3 * pc_x[k] * nsi_1640[k];

        t_2105[k] = f_4 * nsh0_1235[k]
                    - f_5 * nsh1_1235[k]
                    + f_3 * pc_x[k] * nsi_1641[k];

        t_2106[k] = f_4 * nsh0_1236[k]
                    - f_5 * nsh1_1236[k]
                    + f_3 * pc_x[k] * nsi_1642[k];
    }

#pragma omp simd aligned(t_2107, t_2108, t_2109, t_2110, t_2111, pc_x, nsh0_1237, nsh0_1238, \
                         nsh1_1237, nsh1_1238, nsi_1643, nsi_1644, nsi_1645, nsi_1646, \
                         nsi_1647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2107[k] = f_4 * nsh0_1237[k]
                    - f_5 * nsh1_1237[k]
                    + f_3 * pc_x[k] * nsi_1643[k];

        t_2108[k] = f_4 * nsh0_1238[k]
                    - f_5 * nsh1_1238[k]
                    + f_3 * pc_x[k] * nsi_1644[k];

        t_2109[k] = f_3 * pc_x[k] * nsi_1645[k];

        t_2110[k] = f_3 * pc_x[k] * nsi_1646[k];

        t_2111[k] = f_3 * pc_x[k] * nsi_1647[k];
    }

#pragma omp simd aligned(t_2112, t_2113, t_2114, t_2115, t_2116, pc_x, pc_y, msi_1365, \
                         nsh0_1233, nsh1_1233, nsi_1645, nsi_1648, nsi_1649, nsi_1650, \
                         nsi_1651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2112[k] = f_3 * pc_x[k] * nsi_1648[k];

        t_2113[k] = f_3 * pc_x[k] * nsi_1649[k];

        t_2114[k] = f_3 * pc_x[k] * nsi_1650[k];

        t_2115[k] = f_3 * pc_x[k] * nsi_1651[k];

        t_2116[k] = f_22 * msi_1365[k]
                    + f_1 * nsh0_1233[k]
                    - f_2 * nsh1_1233[k]
                    + f_3 * pc_y[k] * nsi_1645[k];
    }

#pragma omp simd aligned(t_2117, t_2118, t_2119, pc_y, pc_z, msi_1337, msi_1367, msi_1368, \
                         nsh0_1235, nsh0_1236, nsh1_1235, nsh1_1236, nsi_1645, nsi_1647, \
                         nsi_1648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2117[k] = f_15 * msi_1337[k]
                    + f_3 * pc_z[k] * nsi_1645[k];

        t_2118[k] = f_22 * msi_1367[k]
                    + f_10 * nsh0_1235[k]
                    - f_11 * nsh1_1235[k]
                    + f_3 * pc_y[k] * nsi_1647[k];

        t_2119[k] = f_22 * msi_1368[k]
                    + f_8 * nsh0_1236[k]
                    - f_9 * nsh1_1236[k]
                    + f_3 * pc_y[k] * nsi_1648[k];
    }

#pragma omp simd aligned(t_2120, t_2121, t_2122, pc_y, msi_1369, msi_1370, msi_1371, \
                         nsh0_1237, nsh0_1238, nsh1_1237, nsh1_1238, nsi_1649, nsi_1650, \
                         nsi_1651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2120[k] = f_22 * msi_1369[k]
                    + f_6 * nsh0_1237[k]
                    - f_7 * nsh1_1237[k]
                    + f_3 * pc_y[k] * nsi_1649[k];

        t_2121[k] = f_22 * msi_1370[k]
                    + f_4 * nsh0_1238[k]
                    - f_5 * nsh1_1238[k]
                    + f_3 * pc_y[k] * nsi_1650[k];

        t_2122[k] = f_22 * msi_1371[k]
                    + f_3 * pc_y[k] * nsi_1651[k];
    }

#pragma omp simd aligned(t_2123, t_2124, t_2125, pc_x, pc_z, msi_1343, nsh0_1238, nsh0_1239, \
                         nsh0_1240, nsh1_1238, nsh1_1239, nsh1_1240, nsi_1651, nsi_1652, \
                         nsi_1653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2123[k] = f_15 * msi_1343[k]
                    + f_1 * nsh0_1238[k]
                    - f_2 * nsh1_1238[k]
                    + f_3 * pc_z[k] * nsi_1651[k];

        t_2124[k] = f_1 * nsh0_1239[k]
                    - f_2 * nsh1_1239[k]
                    + f_3 * pc_x[k] * nsi_1652[k];

        t_2125[k] = f_19 * nsh0_1240[k]
                    - f_20 * nsh1_1240[k]
                    + f_3 * pc_x[k] * nsi_1653[k];
    }

#pragma omp simd aligned(t_2126, t_2127, t_2128, pc_x, nsh0_1241, nsh0_1242, nsh0_1243, \
                         nsh1_1241, nsh1_1242, nsh1_1243, nsi_1654, nsi_1655, \
                         nsi_1656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2126[k] = f_19 * nsh0_1241[k]
                    - f_20 * nsh1_1241[k]
                    + f_3 * pc_x[k] * nsi_1654[k];

        t_2127[k] = f_10 * nsh0_1242[k]
                    - f_11 * nsh1_1242[k]
                    + f_3 * pc_x[k] * nsi_1655[k];

        t_2128[k] = f_10 * nsh0_1243[k]
                    - f_11 * nsh1_1243[k]
                    + f_3 * pc_x[k] * nsi_1656[k];
    }

#pragma omp simd aligned(t_2129, t_2130, t_2131, pc_x, nsh0_1244, nsh0_1245, nsh0_1246, \
                         nsh1_1244, nsh1_1245, nsh1_1246, nsi_1657, nsi_1658, \
                         nsi_1659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2129[k] = f_10 * nsh0_1244[k]
                    - f_11 * nsh1_1244[k]
                    + f_3 * pc_x[k] * nsi_1657[k];

        t_2130[k] = f_8 * nsh0_1245[k]
                    - f_9 * nsh1_1245[k]
                    + f_3 * pc_x[k] * nsi_1658[k];

        t_2131[k] = f_8 * nsh0_1246[k]
                    - f_9 * nsh1_1246[k]
                    + f_3 * pc_x[k] * nsi_1659[k];
    }

#pragma omp simd aligned(t_2132, t_2133, t_2134, pc_x, nsh0_1247, nsh0_1248, nsh0_1249, \
                         nsh1_1247, nsh1_1248, nsh1_1249, nsi_1660, nsi_1661, \
                         nsi_1662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2132[k] = f_8 * nsh0_1247[k]
                    - f_9 * nsh1_1247[k]
                    + f_3 * pc_x[k] * nsi_1660[k];

        t_2133[k] = f_8 * nsh0_1248[k]
                    - f_9 * nsh1_1248[k]
                    + f_3 * pc_x[k] * nsi_1661[k];

        t_2134[k] = f_6 * nsh0_1249[k]
                    - f_7 * nsh1_1249[k]
                    + f_3 * pc_x[k] * nsi_1662[k];
    }

#pragma omp simd aligned(t_2135, t_2136, t_2137, pc_x, nsh0_1250, nsh0_1251, nsh0_1252, \
                         nsh1_1250, nsh1_1251, nsh1_1252, nsi_1663, nsi_1664, \
                         nsi_1665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2135[k] = f_6 * nsh0_1250[k]
                    - f_7 * nsh1_1250[k]
                    + f_3 * pc_x[k] * nsi_1663[k];

        t_2136[k] = f_6 * nsh0_1251[k]
                    - f_7 * nsh1_1251[k]
                    + f_3 * pc_x[k] * nsi_1664[k];

        t_2137[k] = f_6 * nsh0_1252[k]
                    - f_7 * nsh1_1252[k]
                    + f_3 * pc_x[k] * nsi_1665[k];
    }

#pragma omp simd aligned(t_2138, t_2139, t_2140, pc_x, nsh0_1253, nsh0_1254, nsh0_1255, \
                         nsh1_1253, nsh1_1254, nsh1_1255, nsi_1666, nsi_1667, \
                         nsi_1668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2138[k] = f_6 * nsh0_1253[k]
                    - f_7 * nsh1_1253[k]
                    + f_3 * pc_x[k] * nsi_1666[k];

        t_2139[k] = f_4 * nsh0_1254[k]
                    - f_5 * nsh1_1254[k]
                    + f_3 * pc_x[k] * nsi_1667[k];

        t_2140[k] = f_4 * nsh0_1255[k]
                    - f_5 * nsh1_1255[k]
                    + f_3 * pc_x[k] * nsi_1668[k];
    }

#pragma omp simd aligned(t_2141, t_2142, t_2143, pc_x, nsh0_1256, nsh0_1257, nsh0_1258, \
                         nsh1_1256, nsh1_1257, nsh1_1258, nsi_1669, nsi_1670, \
                         nsi_1671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2141[k] = f_4 * nsh0_1256[k]
                    - f_5 * nsh1_1256[k]
                    + f_3 * pc_x[k] * nsi_1669[k];

        t_2142[k] = f_4 * nsh0_1257[k]
                    - f_5 * nsh1_1257[k]
                    + f_3 * pc_x[k] * nsi_1670[k];

        t_2143[k] = f_4 * nsh0_1258[k]
                    - f_5 * nsh1_1258[k]
                    + f_3 * pc_x[k] * nsi_1671[k];
    }

#pragma omp simd aligned(t_2144, t_2145, t_2146, t_2147, t_2148, t_2149, pc_x, nsh0_1259, \
                         nsh1_1259, nsi_1672, nsi_1673, nsi_1674, nsi_1675, nsi_1676, \
                         nsi_1677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2144[k] = f_4 * nsh0_1259[k]
                    - f_5 * nsh1_1259[k]
                    + f_3 * pc_x[k] * nsi_1672[k];

        t_2145[k] = f_3 * pc_x[k] * nsi_1673[k];

        t_2146[k] = f_3 * pc_x[k] * nsi_1674[k];

        t_2147[k] = f_3 * pc_x[k] * nsi_1675[k];

        t_2148[k] = f_3 * pc_x[k] * nsi_1676[k];

        t_2149[k] = f_3 * pc_x[k] * nsi_1677[k];
    }

#pragma omp simd aligned(t_2150, t_2151, t_2152, t_2153, pc_x, pc_y, pc_z, msi_1365, msi_1393, \
                         nsh0_1254, nsh1_1254, nsi_1673, nsi_1678, \
                         nsi_1679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2150[k] = f_3 * pc_x[k] * nsi_1678[k];

        t_2151[k] = f_3 * pc_x[k] * nsi_1679[k];

        t_2152[k] = f_23 * msi_1393[k]
                    + f_1 * nsh0_1254[k]
                    - f_2 * nsh1_1254[k]
                    + f_3 * pc_y[k] * nsi_1673[k];

        t_2153[k] = f_16 * msi_1365[k]
                    + f_3 * pc_z[k] * nsi_1673[k];
    }

#pragma omp simd aligned(t_2154, t_2155, t_2156, pc_y, msi_1395, msi_1396, msi_1397, \
                         nsh0_1256, nsh0_1257, nsh0_1258, nsh1_1256, nsh1_1257, nsh1_1258, \
                         nsi_1675, nsi_1676, nsi_1677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2154[k] = f_23 * msi_1395[k]
                    + f_10 * nsh0_1256[k]
                    - f_11 * nsh1_1256[k]
                    + f_3 * pc_y[k] * nsi_1675[k];

        t_2155[k] = f_23 * msi_1396[k]
                    + f_8 * nsh0_1257[k]
                    - f_9 * nsh1_1257[k]
                    + f_3 * pc_y[k] * nsi_1676[k];

        t_2156[k] = f_23 * msi_1397[k]
                    + f_6 * nsh0_1258[k]
                    - f_7 * nsh1_1258[k]
                    + f_3 * pc_y[k] * nsi_1677[k];
    }

#pragma omp simd aligned(t_2157, t_2158, t_2159, pc_y, pc_z, msi_1371, msi_1398, msi_1399, \
                         nsh0_1259, nsh1_1259, nsi_1678, nsi_1679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2157[k] = f_23 * msi_1398[k]
                    + f_4 * nsh0_1259[k]
                    - f_5 * nsh1_1259[k]
                    + f_3 * pc_y[k] * nsi_1678[k];

        t_2158[k] = f_23 * msi_1399[k]
                    + f_3 * pc_y[k] * nsi_1679[k];

        t_2159[k] = f_16 * msi_1371[k]
                    + f_1 * nsh0_1259[k]
                    - f_2 * nsh1_1259[k]
                    + f_3 * pc_z[k] * nsi_1679[k];
    }

#pragma omp simd aligned(t_2160, t_2161, t_2162, pc_x, nsh0_1260, nsh0_1261, nsh0_1262, \
                         nsh1_1260, nsh1_1261, nsh1_1262, nsi_1680, nsi_1681, \
                         nsi_1682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2160[k] = f_1 * nsh0_1260[k]
                    - f_2 * nsh1_1260[k]
                    + f_3 * pc_x[k] * nsi_1680[k];

        t_2161[k] = f_19 * nsh0_1261[k]
                    - f_20 * nsh1_1261[k]
                    + f_3 * pc_x[k] * nsi_1681[k];

        t_2162[k] = f_19 * nsh0_1262[k]
                    - f_20 * nsh1_1262[k]
                    + f_3 * pc_x[k] * nsi_1682[k];
    }

#pragma omp simd aligned(t_2163, t_2164, t_2165, pc_x, nsh0_1263, nsh0_1264, nsh0_1265, \
                         nsh1_1263, nsh1_1264, nsh1_1265, nsi_1683, nsi_1684, \
                         nsi_1685 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2163[k] = f_10 * nsh0_1263[k]
                    - f_11 * nsh1_1263[k]
                    + f_3 * pc_x[k] * nsi_1683[k];

        t_2164[k] = f_10 * nsh0_1264[k]
                    - f_11 * nsh1_1264[k]
                    + f_3 * pc_x[k] * nsi_1684[k];

        t_2165[k] = f_10 * nsh0_1265[k]
                    - f_11 * nsh1_1265[k]
                    + f_3 * pc_x[k] * nsi_1685[k];
    }

#pragma omp simd aligned(t_2166, t_2167, t_2168, pc_x, nsh0_1266, nsh0_1267, nsh0_1268, \
                         nsh1_1266, nsh1_1267, nsh1_1268, nsi_1686, nsi_1687, \
                         nsi_1688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2166[k] = f_8 * nsh0_1266[k]
                    - f_9 * nsh1_1266[k]
                    + f_3 * pc_x[k] * nsi_1686[k];

        t_2167[k] = f_8 * nsh0_1267[k]
                    - f_9 * nsh1_1267[k]
                    + f_3 * pc_x[k] * nsi_1687[k];

        t_2168[k] = f_8 * nsh0_1268[k]
                    - f_9 * nsh1_1268[k]
                    + f_3 * pc_x[k] * nsi_1688[k];
    }

#pragma omp simd aligned(t_2169, t_2170, t_2171, pc_x, nsh0_1269, nsh0_1270, nsh0_1271, \
                         nsh1_1269, nsh1_1270, nsh1_1271, nsi_1689, nsi_1690, \
                         nsi_1691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2169[k] = f_8 * nsh0_1269[k]
                    - f_9 * nsh1_1269[k]
                    + f_3 * pc_x[k] * nsi_1689[k];

        t_2170[k] = f_6 * nsh0_1270[k]
                    - f_7 * nsh1_1270[k]
                    + f_3 * pc_x[k] * nsi_1690[k];

        t_2171[k] = f_6 * nsh0_1271[k]
                    - f_7 * nsh1_1271[k]
                    + f_3 * pc_x[k] * nsi_1691[k];
    }

#pragma omp simd aligned(t_2172, t_2173, t_2174, pc_x, nsh0_1272, nsh0_1273, nsh0_1274, \
                         nsh1_1272, nsh1_1273, nsh1_1274, nsi_1692, nsi_1693, \
                         nsi_1694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2172[k] = f_6 * nsh0_1272[k]
                    - f_7 * nsh1_1272[k]
                    + f_3 * pc_x[k] * nsi_1692[k];

        t_2173[k] = f_6 * nsh0_1273[k]
                    - f_7 * nsh1_1273[k]
                    + f_3 * pc_x[k] * nsi_1693[k];

        t_2174[k] = f_6 * nsh0_1274[k]
                    - f_7 * nsh1_1274[k]
                    + f_3 * pc_x[k] * nsi_1694[k];
    }

#pragma omp simd aligned(t_2175, t_2176, t_2177, pc_x, nsh0_1275, nsh0_1276, nsh0_1277, \
                         nsh1_1275, nsh1_1276, nsh1_1277, nsi_1695, nsi_1696, \
                         nsi_1697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2175[k] = f_4 * nsh0_1275[k]
                    - f_5 * nsh1_1275[k]
                    + f_3 * pc_x[k] * nsi_1695[k];

        t_2176[k] = f_4 * nsh0_1276[k]
                    - f_5 * nsh1_1276[k]
                    + f_3 * pc_x[k] * nsi_1696[k];

        t_2177[k] = f_4 * nsh0_1277[k]
                    - f_5 * nsh1_1277[k]
                    + f_3 * pc_x[k] * nsi_1697[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, t_2181, pc_x, nsh0_1278, nsh0_1279, \
                         nsh0_1280, nsh1_1278, nsh1_1279, nsh1_1280, nsi_1698, nsi_1699, \
                         nsi_1700, nsi_1701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = f_4 * nsh0_1278[k]
                    - f_5 * nsh1_1278[k]
                    + f_3 * pc_x[k] * nsi_1698[k];

        t_2179[k] = f_4 * nsh0_1279[k]
                    - f_5 * nsh1_1279[k]
                    + f_3 * pc_x[k] * nsi_1699[k];

        t_2180[k] = f_4 * nsh0_1280[k]
                    - f_5 * nsh1_1280[k]
                    + f_3 * pc_x[k] * nsi_1700[k];

        t_2181[k] = f_3 * pc_x[k] * nsi_1701[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece19(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t msi, const size_t nsh0,
                                                           const size_t nsh1, const size_t nsi,
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
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi_1393 = buffer.data(msi + 1393);
    const auto *msi_1399 = buffer.data(msi + 1399);
    const auto *msi_1421 = buffer.data(msi + 1421);
    const auto *msi_1423 = buffer.data(msi + 1423);
    const auto *msi_1424 = buffer.data(msi + 1424);
    const auto *msi_1425 = buffer.data(msi + 1425);
    const auto *msi_1426 = buffer.data(msi + 1426);
    const auto *msi_1427 = buffer.data(msi + 1427);
    const auto *msi_1449 = buffer.data(msi + 1449);
    const auto *msi_1451 = buffer.data(msi + 1451);
    const auto *msi_1452 = buffer.data(msi + 1452);
    const auto *msi_1453 = buffer.data(msi + 1453);
    const auto *msi_1454 = buffer.data(msi + 1454);
    const auto *msi_1455 = buffer.data(msi + 1455);
    const auto *msi_1477 = buffer.data(msi + 1477);
    const auto *msi_1479 = buffer.data(msi + 1479);
    const auto *msi_1480 = buffer.data(msi + 1480);
    const auto *msi_1481 = buffer.data(msi + 1481);
    const auto *msi_1482 = buffer.data(msi + 1482);
    const auto *msi_1483 = buffer.data(msi + 1483);
    const auto *msi_1505 = buffer.data(msi + 1505);
    const auto *msi_1507 = buffer.data(msi + 1507);

    const auto *nsh0_1275 = buffer.data(nsh0 + 1275);
    const auto *nsh0_1277 = buffer.data(nsh0 + 1277);
    const auto *nsh0_1278 = buffer.data(nsh0 + 1278);
    const auto *nsh0_1279 = buffer.data(nsh0 + 1279);
    const auto *nsh0_1280 = buffer.data(nsh0 + 1280);
    const auto *nsh0_1281 = buffer.data(nsh0 + 1281);
    const auto *nsh0_1282 = buffer.data(nsh0 + 1282);
    const auto *nsh0_1283 = buffer.data(nsh0 + 1283);
    const auto *nsh0_1284 = buffer.data(nsh0 + 1284);
    const auto *nsh0_1285 = buffer.data(nsh0 + 1285);
    const auto *nsh0_1286 = buffer.data(nsh0 + 1286);
    const auto *nsh0_1287 = buffer.data(nsh0 + 1287);
    const auto *nsh0_1288 = buffer.data(nsh0 + 1288);
    const auto *nsh0_1289 = buffer.data(nsh0 + 1289);
    const auto *nsh0_1290 = buffer.data(nsh0 + 1290);
    const auto *nsh0_1291 = buffer.data(nsh0 + 1291);
    const auto *nsh0_1292 = buffer.data(nsh0 + 1292);
    const auto *nsh0_1293 = buffer.data(nsh0 + 1293);
    const auto *nsh0_1294 = buffer.data(nsh0 + 1294);
    const auto *nsh0_1295 = buffer.data(nsh0 + 1295);
    const auto *nsh0_1296 = buffer.data(nsh0 + 1296);
    const auto *nsh0_1297 = buffer.data(nsh0 + 1297);
    const auto *nsh0_1298 = buffer.data(nsh0 + 1298);
    const auto *nsh0_1299 = buffer.data(nsh0 + 1299);
    const auto *nsh0_1300 = buffer.data(nsh0 + 1300);
    const auto *nsh0_1301 = buffer.data(nsh0 + 1301);
    const auto *nsh0_1302 = buffer.data(nsh0 + 1302);
    const auto *nsh0_1303 = buffer.data(nsh0 + 1303);
    const auto *nsh0_1304 = buffer.data(nsh0 + 1304);
    const auto *nsh0_1305 = buffer.data(nsh0 + 1305);
    const auto *nsh0_1306 = buffer.data(nsh0 + 1306);
    const auto *nsh0_1307 = buffer.data(nsh0 + 1307);
    const auto *nsh0_1308 = buffer.data(nsh0 + 1308);
    const auto *nsh0_1309 = buffer.data(nsh0 + 1309);
    const auto *nsh0_1310 = buffer.data(nsh0 + 1310);
    const auto *nsh0_1311 = buffer.data(nsh0 + 1311);
    const auto *nsh0_1312 = buffer.data(nsh0 + 1312);
    const auto *nsh0_1313 = buffer.data(nsh0 + 1313);
    const auto *nsh0_1314 = buffer.data(nsh0 + 1314);
    const auto *nsh0_1315 = buffer.data(nsh0 + 1315);
    const auto *nsh0_1316 = buffer.data(nsh0 + 1316);
    const auto *nsh0_1317 = buffer.data(nsh0 + 1317);
    const auto *nsh0_1318 = buffer.data(nsh0 + 1318);
    const auto *nsh0_1319 = buffer.data(nsh0 + 1319);
    const auto *nsh0_1320 = buffer.data(nsh0 + 1320);
    const auto *nsh0_1321 = buffer.data(nsh0 + 1321);
    const auto *nsh0_1322 = buffer.data(nsh0 + 1322);
    const auto *nsh0_1323 = buffer.data(nsh0 + 1323);
    const auto *nsh0_1324 = buffer.data(nsh0 + 1324);
    const auto *nsh0_1325 = buffer.data(nsh0 + 1325);
    const auto *nsh0_1326 = buffer.data(nsh0 + 1326);
    const auto *nsh0_1327 = buffer.data(nsh0 + 1327);
    const auto *nsh0_1328 = buffer.data(nsh0 + 1328);
    const auto *nsh0_1329 = buffer.data(nsh0 + 1329);
    const auto *nsh0_1330 = buffer.data(nsh0 + 1330);
    const auto *nsh0_1331 = buffer.data(nsh0 + 1331);
    const auto *nsh0_1332 = buffer.data(nsh0 + 1332);
    const auto *nsh0_1333 = buffer.data(nsh0 + 1333);
    const auto *nsh0_1334 = buffer.data(nsh0 + 1334);
    const auto *nsh0_1335 = buffer.data(nsh0 + 1335);
    const auto *nsh0_1336 = buffer.data(nsh0 + 1336);
    const auto *nsh0_1337 = buffer.data(nsh0 + 1337);
    const auto *nsh0_1338 = buffer.data(nsh0 + 1338);
    const auto *nsh0_1339 = buffer.data(nsh0 + 1339);
    const auto *nsh0_1340 = buffer.data(nsh0 + 1340);
    const auto *nsh0_1341 = buffer.data(nsh0 + 1341);
    const auto *nsh0_1342 = buffer.data(nsh0 + 1342);
    const auto *nsh0_1343 = buffer.data(nsh0 + 1343);

    const auto *nsh1_1275 = buffer.data(nsh1 + 1275);
    const auto *nsh1_1277 = buffer.data(nsh1 + 1277);
    const auto *nsh1_1278 = buffer.data(nsh1 + 1278);
    const auto *nsh1_1279 = buffer.data(nsh1 + 1279);
    const auto *nsh1_1280 = buffer.data(nsh1 + 1280);
    const auto *nsh1_1281 = buffer.data(nsh1 + 1281);
    const auto *nsh1_1282 = buffer.data(nsh1 + 1282);
    const auto *nsh1_1283 = buffer.data(nsh1 + 1283);
    const auto *nsh1_1284 = buffer.data(nsh1 + 1284);
    const auto *nsh1_1285 = buffer.data(nsh1 + 1285);
    const auto *nsh1_1286 = buffer.data(nsh1 + 1286);
    const auto *nsh1_1287 = buffer.data(nsh1 + 1287);
    const auto *nsh1_1288 = buffer.data(nsh1 + 1288);
    const auto *nsh1_1289 = buffer.data(nsh1 + 1289);
    const auto *nsh1_1290 = buffer.data(nsh1 + 1290);
    const auto *nsh1_1291 = buffer.data(nsh1 + 1291);
    const auto *nsh1_1292 = buffer.data(nsh1 + 1292);
    const auto *nsh1_1293 = buffer.data(nsh1 + 1293);
    const auto *nsh1_1294 = buffer.data(nsh1 + 1294);
    const auto *nsh1_1295 = buffer.data(nsh1 + 1295);
    const auto *nsh1_1296 = buffer.data(nsh1 + 1296);
    const auto *nsh1_1297 = buffer.data(nsh1 + 1297);
    const auto *nsh1_1298 = buffer.data(nsh1 + 1298);
    const auto *nsh1_1299 = buffer.data(nsh1 + 1299);
    const auto *nsh1_1300 = buffer.data(nsh1 + 1300);
    const auto *nsh1_1301 = buffer.data(nsh1 + 1301);
    const auto *nsh1_1302 = buffer.data(nsh1 + 1302);
    const auto *nsh1_1303 = buffer.data(nsh1 + 1303);
    const auto *nsh1_1304 = buffer.data(nsh1 + 1304);
    const auto *nsh1_1305 = buffer.data(nsh1 + 1305);
    const auto *nsh1_1306 = buffer.data(nsh1 + 1306);
    const auto *nsh1_1307 = buffer.data(nsh1 + 1307);
    const auto *nsh1_1308 = buffer.data(nsh1 + 1308);
    const auto *nsh1_1309 = buffer.data(nsh1 + 1309);
    const auto *nsh1_1310 = buffer.data(nsh1 + 1310);
    const auto *nsh1_1311 = buffer.data(nsh1 + 1311);
    const auto *nsh1_1312 = buffer.data(nsh1 + 1312);
    const auto *nsh1_1313 = buffer.data(nsh1 + 1313);
    const auto *nsh1_1314 = buffer.data(nsh1 + 1314);
    const auto *nsh1_1315 = buffer.data(nsh1 + 1315);
    const auto *nsh1_1316 = buffer.data(nsh1 + 1316);
    const auto *nsh1_1317 = buffer.data(nsh1 + 1317);
    const auto *nsh1_1318 = buffer.data(nsh1 + 1318);
    const auto *nsh1_1319 = buffer.data(nsh1 + 1319);
    const auto *nsh1_1320 = buffer.data(nsh1 + 1320);
    const auto *nsh1_1321 = buffer.data(nsh1 + 1321);
    const auto *nsh1_1322 = buffer.data(nsh1 + 1322);
    const auto *nsh1_1323 = buffer.data(nsh1 + 1323);
    const auto *nsh1_1324 = buffer.data(nsh1 + 1324);
    const auto *nsh1_1325 = buffer.data(nsh1 + 1325);
    const auto *nsh1_1326 = buffer.data(nsh1 + 1326);
    const auto *nsh1_1327 = buffer.data(nsh1 + 1327);
    const auto *nsh1_1328 = buffer.data(nsh1 + 1328);
    const auto *nsh1_1329 = buffer.data(nsh1 + 1329);
    const auto *nsh1_1330 = buffer.data(nsh1 + 1330);
    const auto *nsh1_1331 = buffer.data(nsh1 + 1331);
    const auto *nsh1_1332 = buffer.data(nsh1 + 1332);
    const auto *nsh1_1333 = buffer.data(nsh1 + 1333);
    const auto *nsh1_1334 = buffer.data(nsh1 + 1334);
    const auto *nsh1_1335 = buffer.data(nsh1 + 1335);
    const auto *nsh1_1336 = buffer.data(nsh1 + 1336);
    const auto *nsh1_1337 = buffer.data(nsh1 + 1337);
    const auto *nsh1_1338 = buffer.data(nsh1 + 1338);
    const auto *nsh1_1339 = buffer.data(nsh1 + 1339);
    const auto *nsh1_1340 = buffer.data(nsh1 + 1340);
    const auto *nsh1_1341 = buffer.data(nsh1 + 1341);
    const auto *nsh1_1342 = buffer.data(nsh1 + 1342);
    const auto *nsh1_1343 = buffer.data(nsh1 + 1343);

    const auto *nsi_1701 = buffer.data(nsi + 1701);
    const auto *nsi_1702 = buffer.data(nsi + 1702);
    const auto *nsi_1703 = buffer.data(nsi + 1703);
    const auto *nsi_1704 = buffer.data(nsi + 1704);
    const auto *nsi_1705 = buffer.data(nsi + 1705);
    const auto *nsi_1706 = buffer.data(nsi + 1706);
    const auto *nsi_1707 = buffer.data(nsi + 1707);
    const auto *nsi_1708 = buffer.data(nsi + 1708);
    const auto *nsi_1709 = buffer.data(nsi + 1709);
    const auto *nsi_1710 = buffer.data(nsi + 1710);
    const auto *nsi_1711 = buffer.data(nsi + 1711);
    const auto *nsi_1712 = buffer.data(nsi + 1712);
    const auto *nsi_1713 = buffer.data(nsi + 1713);
    const auto *nsi_1714 = buffer.data(nsi + 1714);
    const auto *nsi_1715 = buffer.data(nsi + 1715);
    const auto *nsi_1716 = buffer.data(nsi + 1716);
    const auto *nsi_1717 = buffer.data(nsi + 1717);
    const auto *nsi_1718 = buffer.data(nsi + 1718);
    const auto *nsi_1719 = buffer.data(nsi + 1719);
    const auto *nsi_1720 = buffer.data(nsi + 1720);
    const auto *nsi_1721 = buffer.data(nsi + 1721);
    const auto *nsi_1722 = buffer.data(nsi + 1722);
    const auto *nsi_1723 = buffer.data(nsi + 1723);
    const auto *nsi_1724 = buffer.data(nsi + 1724);
    const auto *nsi_1725 = buffer.data(nsi + 1725);
    const auto *nsi_1726 = buffer.data(nsi + 1726);
    const auto *nsi_1727 = buffer.data(nsi + 1727);
    const auto *nsi_1728 = buffer.data(nsi + 1728);
    const auto *nsi_1729 = buffer.data(nsi + 1729);
    const auto *nsi_1730 = buffer.data(nsi + 1730);
    const auto *nsi_1731 = buffer.data(nsi + 1731);
    const auto *nsi_1732 = buffer.data(nsi + 1732);
    const auto *nsi_1733 = buffer.data(nsi + 1733);
    const auto *nsi_1734 = buffer.data(nsi + 1734);
    const auto *nsi_1735 = buffer.data(nsi + 1735);
    const auto *nsi_1736 = buffer.data(nsi + 1736);
    const auto *nsi_1737 = buffer.data(nsi + 1737);
    const auto *nsi_1738 = buffer.data(nsi + 1738);
    const auto *nsi_1739 = buffer.data(nsi + 1739);
    const auto *nsi_1740 = buffer.data(nsi + 1740);
    const auto *nsi_1741 = buffer.data(nsi + 1741);
    const auto *nsi_1742 = buffer.data(nsi + 1742);
    const auto *nsi_1743 = buffer.data(nsi + 1743);
    const auto *nsi_1744 = buffer.data(nsi + 1744);
    const auto *nsi_1745 = buffer.data(nsi + 1745);
    const auto *nsi_1746 = buffer.data(nsi + 1746);
    const auto *nsi_1747 = buffer.data(nsi + 1747);
    const auto *nsi_1748 = buffer.data(nsi + 1748);
    const auto *nsi_1749 = buffer.data(nsi + 1749);
    const auto *nsi_1750 = buffer.data(nsi + 1750);
    const auto *nsi_1751 = buffer.data(nsi + 1751);
    const auto *nsi_1752 = buffer.data(nsi + 1752);
    const auto *nsi_1753 = buffer.data(nsi + 1753);
    const auto *nsi_1754 = buffer.data(nsi + 1754);
    const auto *nsi_1755 = buffer.data(nsi + 1755);
    const auto *nsi_1756 = buffer.data(nsi + 1756);
    const auto *nsi_1757 = buffer.data(nsi + 1757);
    const auto *nsi_1758 = buffer.data(nsi + 1758);
    const auto *nsi_1759 = buffer.data(nsi + 1759);
    const auto *nsi_1760 = buffer.data(nsi + 1760);
    const auto *nsi_1761 = buffer.data(nsi + 1761);
    const auto *nsi_1762 = buffer.data(nsi + 1762);
    const auto *nsi_1763 = buffer.data(nsi + 1763);
    const auto *nsi_1764 = buffer.data(nsi + 1764);
    const auto *nsi_1765 = buffer.data(nsi + 1765);
    const auto *nsi_1766 = buffer.data(nsi + 1766);
    const auto *nsi_1767 = buffer.data(nsi + 1767);
    const auto *nsi_1768 = buffer.data(nsi + 1768);
    const auto *nsi_1769 = buffer.data(nsi + 1769);
    const auto *nsi_1770 = buffer.data(nsi + 1770);
    const auto *nsi_1771 = buffer.data(nsi + 1771);
    const auto *nsi_1772 = buffer.data(nsi + 1772);
    const auto *nsi_1773 = buffer.data(nsi + 1773);
    const auto *nsi_1774 = buffer.data(nsi + 1774);
    const auto *nsi_1775 = buffer.data(nsi + 1775);
    const auto *nsi_1776 = buffer.data(nsi + 1776);
    const auto *nsi_1777 = buffer.data(nsi + 1777);
    const auto *nsi_1778 = buffer.data(nsi + 1778);
    const auto *nsi_1779 = buffer.data(nsi + 1779);
    const auto *nsi_1780 = buffer.data(nsi + 1780);
    const auto *nsi_1781 = buffer.data(nsi + 1781);
    const auto *nsi_1782 = buffer.data(nsi + 1782);
    const auto *nsi_1783 = buffer.data(nsi + 1783);
    const auto *nsi_1784 = buffer.data(nsi + 1784);
    const auto *nsi_1785 = buffer.data(nsi + 1785);
    const auto *nsi_1786 = buffer.data(nsi + 1786);
    const auto *nsi_1787 = buffer.data(nsi + 1787);
    const auto *nsi_1788 = buffer.data(nsi + 1788);
    const auto *nsi_1789 = buffer.data(nsi + 1789);
    const auto *nsi_1790 = buffer.data(nsi + 1790);
    const auto *nsi_1791 = buffer.data(nsi + 1791);

#pragma omp simd aligned(t_2182, t_2183, t_2184, t_2185, t_2186, t_2187, pc_x, nsi_1702, \
                         nsi_1703, nsi_1704, nsi_1705, nsi_1706, \
                         nsi_1707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2182[k] = f_3 * pc_x[k] * nsi_1702[k];

        t_2183[k] = f_3 * pc_x[k] * nsi_1703[k];

        t_2184[k] = f_3 * pc_x[k] * nsi_1704[k];

        t_2185[k] = f_3 * pc_x[k] * nsi_1705[k];

        t_2186[k] = f_3 * pc_x[k] * nsi_1706[k];

        t_2187[k] = f_3 * pc_x[k] * nsi_1707[k];
    }

#pragma omp simd aligned(t_2188, t_2189, t_2190, pc_y, pc_z, msi_1393, msi_1421, msi_1423, \
                         nsh0_1275, nsh0_1277, nsh1_1275, nsh1_1277, nsi_1701, \
                         nsi_1703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2188[k] = f_17 * msi_1421[k]
                    + f_1 * nsh0_1275[k]
                    - f_2 * nsh1_1275[k]
                    + f_3 * pc_y[k] * nsi_1701[k];

        t_2189[k] = f_17 * msi_1393[k]
                    + f_3 * pc_z[k] * nsi_1701[k];

        t_2190[k] = f_17 * msi_1423[k]
                    + f_10 * nsh0_1277[k]
                    - f_11 * nsh1_1277[k]
                    + f_3 * pc_y[k] * nsi_1703[k];
    }

#pragma omp simd aligned(t_2191, t_2192, t_2193, pc_y, msi_1424, msi_1425, msi_1426, \
                         nsh0_1278, nsh0_1279, nsh0_1280, nsh1_1278, nsh1_1279, nsh1_1280, \
                         nsi_1704, nsi_1705, nsi_1706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2191[k] = f_17 * msi_1424[k]
                    + f_8 * nsh0_1278[k]
                    - f_9 * nsh1_1278[k]
                    + f_3 * pc_y[k] * nsi_1704[k];

        t_2192[k] = f_17 * msi_1425[k]
                    + f_6 * nsh0_1279[k]
                    - f_7 * nsh1_1279[k]
                    + f_3 * pc_y[k] * nsi_1705[k];

        t_2193[k] = f_17 * msi_1426[k]
                    + f_4 * nsh0_1280[k]
                    - f_5 * nsh1_1280[k]
                    + f_3 * pc_y[k] * nsi_1706[k];
    }

#pragma omp simd aligned(t_2194, t_2195, t_2196, pc_x, pc_y, pc_z, msi_1399, msi_1427, \
                         nsh0_1280, nsh0_1281, nsh1_1280, nsh1_1281, nsi_1707, \
                         nsi_1708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2194[k] = f_17 * msi_1427[k]
                    + f_3 * pc_y[k] * nsi_1707[k];

        t_2195[k] = f_17 * msi_1399[k]
                    + f_1 * nsh0_1280[k]
                    - f_2 * nsh1_1280[k]
                    + f_3 * pc_z[k] * nsi_1707[k];

        t_2196[k] = f_1 * nsh0_1281[k]
                    - f_2 * nsh1_1281[k]
                    + f_3 * pc_x[k] * nsi_1708[k];
    }

#pragma omp simd aligned(t_2197, t_2198, t_2199, pc_x, nsh0_1282, nsh0_1283, nsh0_1284, \
                         nsh1_1282, nsh1_1283, nsh1_1284, nsi_1709, nsi_1710, \
                         nsi_1711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2197[k] = f_19 * nsh0_1282[k]
                    - f_20 * nsh1_1282[k]
                    + f_3 * pc_x[k] * nsi_1709[k];

        t_2198[k] = f_19 * nsh0_1283[k]
                    - f_20 * nsh1_1283[k]
                    + f_3 * pc_x[k] * nsi_1710[k];

        t_2199[k] = f_10 * nsh0_1284[k]
                    - f_11 * nsh1_1284[k]
                    + f_3 * pc_x[k] * nsi_1711[k];
    }

#pragma omp simd aligned(t_2200, t_2201, t_2202, pc_x, nsh0_1285, nsh0_1286, nsh0_1287, \
                         nsh1_1285, nsh1_1286, nsh1_1287, nsi_1712, nsi_1713, \
                         nsi_1714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2200[k] = f_10 * nsh0_1285[k]
                    - f_11 * nsh1_1285[k]
                    + f_3 * pc_x[k] * nsi_1712[k];

        t_2201[k] = f_10 * nsh0_1286[k]
                    - f_11 * nsh1_1286[k]
                    + f_3 * pc_x[k] * nsi_1713[k];

        t_2202[k] = f_8 * nsh0_1287[k]
                    - f_9 * nsh1_1287[k]
                    + f_3 * pc_x[k] * nsi_1714[k];
    }

#pragma omp simd aligned(t_2203, t_2204, t_2205, pc_x, nsh0_1288, nsh0_1289, nsh0_1290, \
                         nsh1_1288, nsh1_1289, nsh1_1290, nsi_1715, nsi_1716, \
                         nsi_1717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2203[k] = f_8 * nsh0_1288[k]
                    - f_9 * nsh1_1288[k]
                    + f_3 * pc_x[k] * nsi_1715[k];

        t_2204[k] = f_8 * nsh0_1289[k]
                    - f_9 * nsh1_1289[k]
                    + f_3 * pc_x[k] * nsi_1716[k];

        t_2205[k] = f_8 * nsh0_1290[k]
                    - f_9 * nsh1_1290[k]
                    + f_3 * pc_x[k] * nsi_1717[k];
    }

#pragma omp simd aligned(t_2206, t_2207, t_2208, pc_x, nsh0_1291, nsh0_1292, nsh0_1293, \
                         nsh1_1291, nsh1_1292, nsh1_1293, nsi_1718, nsi_1719, \
                         nsi_1720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2206[k] = f_6 * nsh0_1291[k]
                    - f_7 * nsh1_1291[k]
                    + f_3 * pc_x[k] * nsi_1718[k];

        t_2207[k] = f_6 * nsh0_1292[k]
                    - f_7 * nsh1_1292[k]
                    + f_3 * pc_x[k] * nsi_1719[k];

        t_2208[k] = f_6 * nsh0_1293[k]
                    - f_7 * nsh1_1293[k]
                    + f_3 * pc_x[k] * nsi_1720[k];
    }

#pragma omp simd aligned(t_2209, t_2210, t_2211, pc_x, nsh0_1294, nsh0_1295, nsh0_1296, \
                         nsh1_1294, nsh1_1295, nsh1_1296, nsi_1721, nsi_1722, \
                         nsi_1723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2209[k] = f_6 * nsh0_1294[k]
                    - f_7 * nsh1_1294[k]
                    + f_3 * pc_x[k] * nsi_1721[k];

        t_2210[k] = f_6 * nsh0_1295[k]
                    - f_7 * nsh1_1295[k]
                    + f_3 * pc_x[k] * nsi_1722[k];

        t_2211[k] = f_4 * nsh0_1296[k]
                    - f_5 * nsh1_1296[k]
                    + f_3 * pc_x[k] * nsi_1723[k];
    }

#pragma omp simd aligned(t_2212, t_2213, t_2214, pc_x, nsh0_1297, nsh0_1298, nsh0_1299, \
                         nsh1_1297, nsh1_1298, nsh1_1299, nsi_1724, nsi_1725, \
                         nsi_1726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2212[k] = f_4 * nsh0_1297[k]
                    - f_5 * nsh1_1297[k]
                    + f_3 * pc_x[k] * nsi_1724[k];

        t_2213[k] = f_4 * nsh0_1298[k]
                    - f_5 * nsh1_1298[k]
                    + f_3 * pc_x[k] * nsi_1725[k];

        t_2214[k] = f_4 * nsh0_1299[k]
                    - f_5 * nsh1_1299[k]
                    + f_3 * pc_x[k] * nsi_1726[k];
    }

#pragma omp simd aligned(t_2215, t_2216, t_2217, t_2218, t_2219, pc_x, nsh0_1300, nsh0_1301, \
                         nsh1_1300, nsh1_1301, nsi_1727, nsi_1728, nsi_1729, nsi_1730, \
                         nsi_1731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2215[k] = f_4 * nsh0_1300[k]
                    - f_5 * nsh1_1300[k]
                    + f_3 * pc_x[k] * nsi_1727[k];

        t_2216[k] = f_4 * nsh0_1301[k]
                    - f_5 * nsh1_1301[k]
                    + f_3 * pc_x[k] * nsi_1728[k];

        t_2217[k] = f_3 * pc_x[k] * nsi_1729[k];

        t_2218[k] = f_3 * pc_x[k] * nsi_1730[k];

        t_2219[k] = f_3 * pc_x[k] * nsi_1731[k];
    }

#pragma omp simd aligned(t_2220, t_2221, t_2222, t_2223, t_2224, pc_x, pc_y, msi_1449, \
                         nsh0_1296, nsh1_1296, nsi_1729, nsi_1732, nsi_1733, nsi_1734, \
                         nsi_1735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2220[k] = f_3 * pc_x[k] * nsi_1732[k];

        t_2221[k] = f_3 * pc_x[k] * nsi_1733[k];

        t_2222[k] = f_3 * pc_x[k] * nsi_1734[k];

        t_2223[k] = f_3 * pc_x[k] * nsi_1735[k];

        t_2224[k] = f_16 * msi_1449[k]
                    + f_1 * nsh0_1296[k]
                    - f_2 * nsh1_1296[k]
                    + f_3 * pc_y[k] * nsi_1729[k];
    }

#pragma omp simd aligned(t_2225, t_2226, t_2227, pc_y, pc_z, msi_1421, msi_1451, msi_1452, \
                         nsh0_1298, nsh0_1299, nsh1_1298, nsh1_1299, nsi_1729, nsi_1731, \
                         nsi_1732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2225[k] = f_23 * msi_1421[k]
                    + f_3 * pc_z[k] * nsi_1729[k];

        t_2226[k] = f_16 * msi_1451[k]
                    + f_10 * nsh0_1298[k]
                    - f_11 * nsh1_1298[k]
                    + f_3 * pc_y[k] * nsi_1731[k];

        t_2227[k] = f_16 * msi_1452[k]
                    + f_8 * nsh0_1299[k]
                    - f_9 * nsh1_1299[k]
                    + f_3 * pc_y[k] * nsi_1732[k];
    }

#pragma omp simd aligned(t_2228, t_2229, t_2230, pc_y, msi_1453, msi_1454, msi_1455, \
                         nsh0_1300, nsh0_1301, nsh1_1300, nsh1_1301, nsi_1733, nsi_1734, \
                         nsi_1735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2228[k] = f_16 * msi_1453[k]
                    + f_6 * nsh0_1300[k]
                    - f_7 * nsh1_1300[k]
                    + f_3 * pc_y[k] * nsi_1733[k];

        t_2229[k] = f_16 * msi_1454[k]
                    + f_4 * nsh0_1301[k]
                    - f_5 * nsh1_1301[k]
                    + f_3 * pc_y[k] * nsi_1734[k];

        t_2230[k] = f_16 * msi_1455[k]
                    + f_3 * pc_y[k] * nsi_1735[k];
    }

#pragma omp simd aligned(t_2231, t_2232, t_2233, pc_x, pc_z, msi_1427, nsh0_1301, nsh0_1302, \
                         nsh0_1303, nsh1_1301, nsh1_1302, nsh1_1303, nsi_1735, nsi_1736, \
                         nsi_1737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2231[k] = f_23 * msi_1427[k]
                    + f_1 * nsh0_1301[k]
                    - f_2 * nsh1_1301[k]
                    + f_3 * pc_z[k] * nsi_1735[k];

        t_2232[k] = f_1 * nsh0_1302[k]
                    - f_2 * nsh1_1302[k]
                    + f_3 * pc_x[k] * nsi_1736[k];

        t_2233[k] = f_19 * nsh0_1303[k]
                    - f_20 * nsh1_1303[k]
                    + f_3 * pc_x[k] * nsi_1737[k];
    }

#pragma omp simd aligned(t_2234, t_2235, t_2236, pc_x, nsh0_1304, nsh0_1305, nsh0_1306, \
                         nsh1_1304, nsh1_1305, nsh1_1306, nsi_1738, nsi_1739, \
                         nsi_1740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2234[k] = f_19 * nsh0_1304[k]
                    - f_20 * nsh1_1304[k]
                    + f_3 * pc_x[k] * nsi_1738[k];

        t_2235[k] = f_10 * nsh0_1305[k]
                    - f_11 * nsh1_1305[k]
                    + f_3 * pc_x[k] * nsi_1739[k];

        t_2236[k] = f_10 * nsh0_1306[k]
                    - f_11 * nsh1_1306[k]
                    + f_3 * pc_x[k] * nsi_1740[k];
    }

#pragma omp simd aligned(t_2237, t_2238, t_2239, pc_x, nsh0_1307, nsh0_1308, nsh0_1309, \
                         nsh1_1307, nsh1_1308, nsh1_1309, nsi_1741, nsi_1742, \
                         nsi_1743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2237[k] = f_10 * nsh0_1307[k]
                    - f_11 * nsh1_1307[k]
                    + f_3 * pc_x[k] * nsi_1741[k];

        t_2238[k] = f_8 * nsh0_1308[k]
                    - f_9 * nsh1_1308[k]
                    + f_3 * pc_x[k] * nsi_1742[k];

        t_2239[k] = f_8 * nsh0_1309[k]
                    - f_9 * nsh1_1309[k]
                    + f_3 * pc_x[k] * nsi_1743[k];
    }

#pragma omp simd aligned(t_2240, t_2241, t_2242, pc_x, nsh0_1310, nsh0_1311, nsh0_1312, \
                         nsh1_1310, nsh1_1311, nsh1_1312, nsi_1744, nsi_1745, \
                         nsi_1746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2240[k] = f_8 * nsh0_1310[k]
                    - f_9 * nsh1_1310[k]
                    + f_3 * pc_x[k] * nsi_1744[k];

        t_2241[k] = f_8 * nsh0_1311[k]
                    - f_9 * nsh1_1311[k]
                    + f_3 * pc_x[k] * nsi_1745[k];

        t_2242[k] = f_6 * nsh0_1312[k]
                    - f_7 * nsh1_1312[k]
                    + f_3 * pc_x[k] * nsi_1746[k];
    }

#pragma omp simd aligned(t_2243, t_2244, t_2245, pc_x, nsh0_1313, nsh0_1314, nsh0_1315, \
                         nsh1_1313, nsh1_1314, nsh1_1315, nsi_1747, nsi_1748, \
                         nsi_1749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2243[k] = f_6 * nsh0_1313[k]
                    - f_7 * nsh1_1313[k]
                    + f_3 * pc_x[k] * nsi_1747[k];

        t_2244[k] = f_6 * nsh0_1314[k]
                    - f_7 * nsh1_1314[k]
                    + f_3 * pc_x[k] * nsi_1748[k];

        t_2245[k] = f_6 * nsh0_1315[k]
                    - f_7 * nsh1_1315[k]
                    + f_3 * pc_x[k] * nsi_1749[k];
    }

#pragma omp simd aligned(t_2246, t_2247, t_2248, pc_x, nsh0_1316, nsh0_1317, nsh0_1318, \
                         nsh1_1316, nsh1_1317, nsh1_1318, nsi_1750, nsi_1751, \
                         nsi_1752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2246[k] = f_6 * nsh0_1316[k]
                    - f_7 * nsh1_1316[k]
                    + f_3 * pc_x[k] * nsi_1750[k];

        t_2247[k] = f_4 * nsh0_1317[k]
                    - f_5 * nsh1_1317[k]
                    + f_3 * pc_x[k] * nsi_1751[k];

        t_2248[k] = f_4 * nsh0_1318[k]
                    - f_5 * nsh1_1318[k]
                    + f_3 * pc_x[k] * nsi_1752[k];
    }

#pragma omp simd aligned(t_2249, t_2250, t_2251, pc_x, nsh0_1319, nsh0_1320, nsh0_1321, \
                         nsh1_1319, nsh1_1320, nsh1_1321, nsi_1753, nsi_1754, \
                         nsi_1755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2249[k] = f_4 * nsh0_1319[k]
                    - f_5 * nsh1_1319[k]
                    + f_3 * pc_x[k] * nsi_1753[k];

        t_2250[k] = f_4 * nsh0_1320[k]
                    - f_5 * nsh1_1320[k]
                    + f_3 * pc_x[k] * nsi_1754[k];

        t_2251[k] = f_4 * nsh0_1321[k]
                    - f_5 * nsh1_1321[k]
                    + f_3 * pc_x[k] * nsi_1755[k];
    }

#pragma omp simd aligned(t_2252, t_2253, t_2254, t_2255, t_2256, t_2257, pc_x, nsh0_1322, \
                         nsh1_1322, nsi_1756, nsi_1757, nsi_1758, nsi_1759, nsi_1760, \
                         nsi_1761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2252[k] = f_4 * nsh0_1322[k]
                    - f_5 * nsh1_1322[k]
                    + f_3 * pc_x[k] * nsi_1756[k];

        t_2253[k] = f_3 * pc_x[k] * nsi_1757[k];

        t_2254[k] = f_3 * pc_x[k] * nsi_1758[k];

        t_2255[k] = f_3 * pc_x[k] * nsi_1759[k];

        t_2256[k] = f_3 * pc_x[k] * nsi_1760[k];

        t_2257[k] = f_3 * pc_x[k] * nsi_1761[k];
    }

#pragma omp simd aligned(t_2258, t_2259, t_2260, t_2261, pc_x, pc_y, pc_z, msi_1449, msi_1477, \
                         nsh0_1317, nsh1_1317, nsi_1757, nsi_1762, \
                         nsi_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2258[k] = f_3 * pc_x[k] * nsi_1762[k];

        t_2259[k] = f_3 * pc_x[k] * nsi_1763[k];

        t_2260[k] = f_15 * msi_1477[k]
                    + f_1 * nsh0_1317[k]
                    - f_2 * nsh1_1317[k]
                    + f_3 * pc_y[k] * nsi_1757[k];

        t_2261[k] = f_22 * msi_1449[k]
                    + f_3 * pc_z[k] * nsi_1757[k];
    }

#pragma omp simd aligned(t_2262, t_2263, t_2264, pc_y, msi_1479, msi_1480, msi_1481, \
                         nsh0_1319, nsh0_1320, nsh0_1321, nsh1_1319, nsh1_1320, nsh1_1321, \
                         nsi_1759, nsi_1760, nsi_1761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2262[k] = f_15 * msi_1479[k]
                    + f_10 * nsh0_1319[k]
                    - f_11 * nsh1_1319[k]
                    + f_3 * pc_y[k] * nsi_1759[k];

        t_2263[k] = f_15 * msi_1480[k]
                    + f_8 * nsh0_1320[k]
                    - f_9 * nsh1_1320[k]
                    + f_3 * pc_y[k] * nsi_1760[k];

        t_2264[k] = f_15 * msi_1481[k]
                    + f_6 * nsh0_1321[k]
                    - f_7 * nsh1_1321[k]
                    + f_3 * pc_y[k] * nsi_1761[k];
    }

#pragma omp simd aligned(t_2265, t_2266, t_2267, pc_y, pc_z, msi_1455, msi_1482, msi_1483, \
                         nsh0_1322, nsh1_1322, nsi_1762, nsi_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2265[k] = f_15 * msi_1482[k]
                    + f_4 * nsh0_1322[k]
                    - f_5 * nsh1_1322[k]
                    + f_3 * pc_y[k] * nsi_1762[k];

        t_2266[k] = f_15 * msi_1483[k]
                    + f_3 * pc_y[k] * nsi_1763[k];

        t_2267[k] = f_22 * msi_1455[k]
                    + f_1 * nsh0_1322[k]
                    - f_2 * nsh1_1322[k]
                    + f_3 * pc_z[k] * nsi_1763[k];
    }

#pragma omp simd aligned(t_2268, t_2269, t_2270, pc_x, nsh0_1323, nsh0_1324, nsh0_1325, \
                         nsh1_1323, nsh1_1324, nsh1_1325, nsi_1764, nsi_1765, \
                         nsi_1766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2268[k] = f_1 * nsh0_1323[k]
                    - f_2 * nsh1_1323[k]
                    + f_3 * pc_x[k] * nsi_1764[k];

        t_2269[k] = f_19 * nsh0_1324[k]
                    - f_20 * nsh1_1324[k]
                    + f_3 * pc_x[k] * nsi_1765[k];

        t_2270[k] = f_19 * nsh0_1325[k]
                    - f_20 * nsh1_1325[k]
                    + f_3 * pc_x[k] * nsi_1766[k];
    }

#pragma omp simd aligned(t_2271, t_2272, t_2273, pc_x, nsh0_1326, nsh0_1327, nsh0_1328, \
                         nsh1_1326, nsh1_1327, nsh1_1328, nsi_1767, nsi_1768, \
                         nsi_1769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2271[k] = f_10 * nsh0_1326[k]
                    - f_11 * nsh1_1326[k]
                    + f_3 * pc_x[k] * nsi_1767[k];

        t_2272[k] = f_10 * nsh0_1327[k]
                    - f_11 * nsh1_1327[k]
                    + f_3 * pc_x[k] * nsi_1768[k];

        t_2273[k] = f_10 * nsh0_1328[k]
                    - f_11 * nsh1_1328[k]
                    + f_3 * pc_x[k] * nsi_1769[k];
    }

#pragma omp simd aligned(t_2274, t_2275, t_2276, pc_x, nsh0_1329, nsh0_1330, nsh0_1331, \
                         nsh1_1329, nsh1_1330, nsh1_1331, nsi_1770, nsi_1771, \
                         nsi_1772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2274[k] = f_8 * nsh0_1329[k]
                    - f_9 * nsh1_1329[k]
                    + f_3 * pc_x[k] * nsi_1770[k];

        t_2275[k] = f_8 * nsh0_1330[k]
                    - f_9 * nsh1_1330[k]
                    + f_3 * pc_x[k] * nsi_1771[k];

        t_2276[k] = f_8 * nsh0_1331[k]
                    - f_9 * nsh1_1331[k]
                    + f_3 * pc_x[k] * nsi_1772[k];
    }

#pragma omp simd aligned(t_2277, t_2278, t_2279, pc_x, nsh0_1332, nsh0_1333, nsh0_1334, \
                         nsh1_1332, nsh1_1333, nsh1_1334, nsi_1773, nsi_1774, \
                         nsi_1775 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2277[k] = f_8 * nsh0_1332[k]
                    - f_9 * nsh1_1332[k]
                    + f_3 * pc_x[k] * nsi_1773[k];

        t_2278[k] = f_6 * nsh0_1333[k]
                    - f_7 * nsh1_1333[k]
                    + f_3 * pc_x[k] * nsi_1774[k];

        t_2279[k] = f_6 * nsh0_1334[k]
                    - f_7 * nsh1_1334[k]
                    + f_3 * pc_x[k] * nsi_1775[k];
    }

#pragma omp simd aligned(t_2280, t_2281, t_2282, pc_x, nsh0_1335, nsh0_1336, nsh0_1337, \
                         nsh1_1335, nsh1_1336, nsh1_1337, nsi_1776, nsi_1777, \
                         nsi_1778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2280[k] = f_6 * nsh0_1335[k]
                    - f_7 * nsh1_1335[k]
                    + f_3 * pc_x[k] * nsi_1776[k];

        t_2281[k] = f_6 * nsh0_1336[k]
                    - f_7 * nsh1_1336[k]
                    + f_3 * pc_x[k] * nsi_1777[k];

        t_2282[k] = f_6 * nsh0_1337[k]
                    - f_7 * nsh1_1337[k]
                    + f_3 * pc_x[k] * nsi_1778[k];
    }

#pragma omp simd aligned(t_2283, t_2284, t_2285, pc_x, nsh0_1338, nsh0_1339, nsh0_1340, \
                         nsh1_1338, nsh1_1339, nsh1_1340, nsi_1779, nsi_1780, \
                         nsi_1781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2283[k] = f_4 * nsh0_1338[k]
                    - f_5 * nsh1_1338[k]
                    + f_3 * pc_x[k] * nsi_1779[k];

        t_2284[k] = f_4 * nsh0_1339[k]
                    - f_5 * nsh1_1339[k]
                    + f_3 * pc_x[k] * nsi_1780[k];

        t_2285[k] = f_4 * nsh0_1340[k]
                    - f_5 * nsh1_1340[k]
                    + f_3 * pc_x[k] * nsi_1781[k];
    }

#pragma omp simd aligned(t_2286, t_2287, t_2288, t_2289, pc_x, nsh0_1341, nsh0_1342, \
                         nsh0_1343, nsh1_1341, nsh1_1342, nsh1_1343, nsi_1782, nsi_1783, \
                         nsi_1784, nsi_1785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2286[k] = f_4 * nsh0_1341[k]
                    - f_5 * nsh1_1341[k]
                    + f_3 * pc_x[k] * nsi_1782[k];

        t_2287[k] = f_4 * nsh0_1342[k]
                    - f_5 * nsh1_1342[k]
                    + f_3 * pc_x[k] * nsi_1783[k];

        t_2288[k] = f_4 * nsh0_1343[k]
                    - f_5 * nsh1_1343[k]
                    + f_3 * pc_x[k] * nsi_1784[k];

        t_2289[k] = f_3 * pc_x[k] * nsi_1785[k];
    }

#pragma omp simd aligned(t_2290, t_2291, t_2292, t_2293, t_2294, t_2295, pc_x, nsi_1786, \
                         nsi_1787, nsi_1788, nsi_1789, nsi_1790, \
                         nsi_1791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2290[k] = f_3 * pc_x[k] * nsi_1786[k];

        t_2291[k] = f_3 * pc_x[k] * nsi_1787[k];

        t_2292[k] = f_3 * pc_x[k] * nsi_1788[k];

        t_2293[k] = f_3 * pc_x[k] * nsi_1789[k];

        t_2294[k] = f_3 * pc_x[k] * nsi_1790[k];

        t_2295[k] = f_3 * pc_x[k] * nsi_1791[k];
    }

#pragma omp simd aligned(t_2296, t_2297, t_2298, pc_y, pc_z, msi_1477, msi_1505, msi_1507, \
                         nsh0_1338, nsh0_1340, nsh1_1338, nsh1_1340, nsi_1785, \
                         nsi_1787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2296[k] = f_14 * msi_1505[k]
                    + f_1 * nsh0_1338[k]
                    - f_2 * nsh1_1338[k]
                    + f_3 * pc_y[k] * nsi_1785[k];

        t_2297[k] = f_21 * msi_1477[k]
                    + f_3 * pc_z[k] * nsi_1785[k];

        t_2298[k] = f_14 * msi_1507[k]
                    + f_10 * nsh0_1340[k]
                    - f_11 * nsh1_1340[k]
                    + f_3 * pc_y[k] * nsi_1787[k];
    }
}

static auto
compute_prim_nsk_three_center_electron_repulsion_0_piece20(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msk0,
                                                           const size_t msi, const size_t msk1,
                                                           const size_t nsh0, const size_t nsh1,
                                                           const size_t nsi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msk0_1944 = buffer.data(msk0 + 1944);
    const auto *msk0_1946 = buffer.data(msk0 + 1946);
    const auto *msk0_1949 = buffer.data(msk0 + 1949);
    const auto *msk0_1953 = buffer.data(msk0 + 1953);
    const auto *msk0_1958 = buffer.data(msk0 + 1958);
    const auto *msk0_1964 = buffer.data(msk0 + 1964);
    const auto *msk0_1972 = buffer.data(msk0 + 1972);
    const auto *msk0_1974 = buffer.data(msk0 + 1974);
    const auto *msk0_1975 = buffer.data(msk0 + 1975);
    const auto *msk0_1976 = buffer.data(msk0 + 1976);
    const auto *msk0_1977 = buffer.data(msk0 + 1977);
    const auto *msk0_1979 = buffer.data(msk0 + 1979);

    const auto *msi_1483 = buffer.data(msi + 1483);
    const auto *msi_1505 = buffer.data(msi + 1505);
    const auto *msi_1508 = buffer.data(msi + 1508);
    const auto *msi_1509 = buffer.data(msi + 1509);
    const auto *msi_1510 = buffer.data(msi + 1510);
    const auto *msi_1511 = buffer.data(msi + 1511);
    const auto *msi_1533 = buffer.data(msi + 1533);
    const auto *msi_1535 = buffer.data(msi + 1535);
    const auto *msi_1536 = buffer.data(msi + 1536);
    const auto *msi_1537 = buffer.data(msi + 1537);
    const auto *msi_1538 = buffer.data(msi + 1538);
    const auto *msi_1539 = buffer.data(msi + 1539);

    const auto *msk1_1944 = buffer.data(msk1 + 1944);
    const auto *msk1_1946 = buffer.data(msk1 + 1946);
    const auto *msk1_1949 = buffer.data(msk1 + 1949);
    const auto *msk1_1953 = buffer.data(msk1 + 1953);
    const auto *msk1_1958 = buffer.data(msk1 + 1958);
    const auto *msk1_1964 = buffer.data(msk1 + 1964);
    const auto *msk1_1972 = buffer.data(msk1 + 1972);
    const auto *msk1_1974 = buffer.data(msk1 + 1974);
    const auto *msk1_1975 = buffer.data(msk1 + 1975);
    const auto *msk1_1976 = buffer.data(msk1 + 1976);
    const auto *msk1_1977 = buffer.data(msk1 + 1977);
    const auto *msk1_1979 = buffer.data(msk1 + 1979);

    const auto *nsh0_1341 = buffer.data(nsh0 + 1341);
    const auto *nsh0_1342 = buffer.data(nsh0 + 1342);
    const auto *nsh0_1343 = buffer.data(nsh0 + 1343);
    const auto *nsh0_1345 = buffer.data(nsh0 + 1345);
    const auto *nsh0_1347 = buffer.data(nsh0 + 1347);
    const auto *nsh0_1348 = buffer.data(nsh0 + 1348);
    const auto *nsh0_1350 = buffer.data(nsh0 + 1350);
    const auto *nsh0_1351 = buffer.data(nsh0 + 1351);
    const auto *nsh0_1352 = buffer.data(nsh0 + 1352);
    const auto *nsh0_1354 = buffer.data(nsh0 + 1354);
    const auto *nsh0_1355 = buffer.data(nsh0 + 1355);
    const auto *nsh0_1356 = buffer.data(nsh0 + 1356);
    const auto *nsh0_1357 = buffer.data(nsh0 + 1357);
    const auto *nsh0_1359 = buffer.data(nsh0 + 1359);
    const auto *nsh0_1360 = buffer.data(nsh0 + 1360);
    const auto *nsh0_1361 = buffer.data(nsh0 + 1361);
    const auto *nsh0_1362 = buffer.data(nsh0 + 1362);
    const auto *nsh0_1363 = buffer.data(nsh0 + 1363);
    const auto *nsh0_1365 = buffer.data(nsh0 + 1365);
    const auto *nsh0_1367 = buffer.data(nsh0 + 1367);
    const auto *nsh0_1368 = buffer.data(nsh0 + 1368);
    const auto *nsh0_1370 = buffer.data(nsh0 + 1370);
    const auto *nsh0_1371 = buffer.data(nsh0 + 1371);
    const auto *nsh0_1372 = buffer.data(nsh0 + 1372);
    const auto *nsh0_1374 = buffer.data(nsh0 + 1374);
    const auto *nsh0_1375 = buffer.data(nsh0 + 1375);
    const auto *nsh0_1376 = buffer.data(nsh0 + 1376);
    const auto *nsh0_1377 = buffer.data(nsh0 + 1377);
    const auto *nsh0_1379 = buffer.data(nsh0 + 1379);
    const auto *nsh0_1380 = buffer.data(nsh0 + 1380);
    const auto *nsh0_1381 = buffer.data(nsh0 + 1381);
    const auto *nsh0_1382 = buffer.data(nsh0 + 1382);
    const auto *nsh0_1383 = buffer.data(nsh0 + 1383);
    const auto *nsh0_1384 = buffer.data(nsh0 + 1384);
    const auto *nsh0_1385 = buffer.data(nsh0 + 1385);

    const auto *nsh1_1341 = buffer.data(nsh1 + 1341);
    const auto *nsh1_1342 = buffer.data(nsh1 + 1342);
    const auto *nsh1_1343 = buffer.data(nsh1 + 1343);
    const auto *nsh1_1345 = buffer.data(nsh1 + 1345);
    const auto *nsh1_1347 = buffer.data(nsh1 + 1347);
    const auto *nsh1_1348 = buffer.data(nsh1 + 1348);
    const auto *nsh1_1350 = buffer.data(nsh1 + 1350);
    const auto *nsh1_1351 = buffer.data(nsh1 + 1351);
    const auto *nsh1_1352 = buffer.data(nsh1 + 1352);
    const auto *nsh1_1354 = buffer.data(nsh1 + 1354);
    const auto *nsh1_1355 = buffer.data(nsh1 + 1355);
    const auto *nsh1_1356 = buffer.data(nsh1 + 1356);
    const auto *nsh1_1357 = buffer.data(nsh1 + 1357);
    const auto *nsh1_1359 = buffer.data(nsh1 + 1359);
    const auto *nsh1_1360 = buffer.data(nsh1 + 1360);
    const auto *nsh1_1361 = buffer.data(nsh1 + 1361);
    const auto *nsh1_1362 = buffer.data(nsh1 + 1362);
    const auto *nsh1_1363 = buffer.data(nsh1 + 1363);
    const auto *nsh1_1365 = buffer.data(nsh1 + 1365);
    const auto *nsh1_1367 = buffer.data(nsh1 + 1367);
    const auto *nsh1_1368 = buffer.data(nsh1 + 1368);
    const auto *nsh1_1370 = buffer.data(nsh1 + 1370);
    const auto *nsh1_1371 = buffer.data(nsh1 + 1371);
    const auto *nsh1_1372 = buffer.data(nsh1 + 1372);
    const auto *nsh1_1374 = buffer.data(nsh1 + 1374);
    const auto *nsh1_1375 = buffer.data(nsh1 + 1375);
    const auto *nsh1_1376 = buffer.data(nsh1 + 1376);
    const auto *nsh1_1377 = buffer.data(nsh1 + 1377);
    const auto *nsh1_1379 = buffer.data(nsh1 + 1379);
    const auto *nsh1_1380 = buffer.data(nsh1 + 1380);
    const auto *nsh1_1381 = buffer.data(nsh1 + 1381);
    const auto *nsh1_1382 = buffer.data(nsh1 + 1382);
    const auto *nsh1_1383 = buffer.data(nsh1 + 1383);
    const auto *nsh1_1384 = buffer.data(nsh1 + 1384);
    const auto *nsh1_1385 = buffer.data(nsh1 + 1385);

    const auto *nsi_1788 = buffer.data(nsi + 1788);
    const auto *nsi_1789 = buffer.data(nsi + 1789);
    const auto *nsi_1790 = buffer.data(nsi + 1790);
    const auto *nsi_1791 = buffer.data(nsi + 1791);
    const auto *nsi_1793 = buffer.data(nsi + 1793);
    const auto *nsi_1795 = buffer.data(nsi + 1795);
    const auto *nsi_1796 = buffer.data(nsi + 1796);
    const auto *nsi_1798 = buffer.data(nsi + 1798);
    const auto *nsi_1799 = buffer.data(nsi + 1799);
    const auto *nsi_1800 = buffer.data(nsi + 1800);
    const auto *nsi_1802 = buffer.data(nsi + 1802);
    const auto *nsi_1803 = buffer.data(nsi + 1803);
    const auto *nsi_1804 = buffer.data(nsi + 1804);
    const auto *nsi_1805 = buffer.data(nsi + 1805);
    const auto *nsi_1807 = buffer.data(nsi + 1807);
    const auto *nsi_1808 = buffer.data(nsi + 1808);
    const auto *nsi_1809 = buffer.data(nsi + 1809);
    const auto *nsi_1810 = buffer.data(nsi + 1810);
    const auto *nsi_1811 = buffer.data(nsi + 1811);
    const auto *nsi_1813 = buffer.data(nsi + 1813);
    const auto *nsi_1814 = buffer.data(nsi + 1814);
    const auto *nsi_1815 = buffer.data(nsi + 1815);
    const auto *nsi_1816 = buffer.data(nsi + 1816);
    const auto *nsi_1817 = buffer.data(nsi + 1817);
    const auto *nsi_1818 = buffer.data(nsi + 1818);
    const auto *nsi_1819 = buffer.data(nsi + 1819);
    const auto *nsi_1820 = buffer.data(nsi + 1820);
    const auto *nsi_1822 = buffer.data(nsi + 1822);
    const auto *nsi_1823 = buffer.data(nsi + 1823);
    const auto *nsi_1825 = buffer.data(nsi + 1825);
    const auto *nsi_1826 = buffer.data(nsi + 1826);
    const auto *nsi_1827 = buffer.data(nsi + 1827);
    const auto *nsi_1829 = buffer.data(nsi + 1829);
    const auto *nsi_1830 = buffer.data(nsi + 1830);
    const auto *nsi_1831 = buffer.data(nsi + 1831);
    const auto *nsi_1832 = buffer.data(nsi + 1832);
    const auto *nsi_1834 = buffer.data(nsi + 1834);
    const auto *nsi_1835 = buffer.data(nsi + 1835);
    const auto *nsi_1836 = buffer.data(nsi + 1836);
    const auto *nsi_1837 = buffer.data(nsi + 1837);
    const auto *nsi_1838 = buffer.data(nsi + 1838);
    const auto *nsi_1840 = buffer.data(nsi + 1840);
    const auto *nsi_1841 = buffer.data(nsi + 1841);
    const auto *nsi_1842 = buffer.data(nsi + 1842);
    const auto *nsi_1843 = buffer.data(nsi + 1843);
    const auto *nsi_1844 = buffer.data(nsi + 1844);
    const auto *nsi_1845 = buffer.data(nsi + 1845);
    const auto *nsi_1846 = buffer.data(nsi + 1846);
    const auto *nsi_1847 = buffer.data(nsi + 1847);

#pragma omp simd aligned(t_2299, t_2300, t_2301, pc_y, msi_1508, msi_1509, msi_1510, \
                         nsh0_1341, nsh0_1342, nsh0_1343, nsh1_1341, nsh1_1342, nsh1_1343, \
                         nsi_1788, nsi_1789, nsi_1790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2299[k] = f_14 * msi_1508[k]
                    + f_8 * nsh0_1341[k]
                    - f_9 * nsh1_1341[k]
                    + f_3 * pc_y[k] * nsi_1788[k];

        t_2300[k] = f_14 * msi_1509[k]
                    + f_6 * nsh0_1342[k]
                    - f_7 * nsh1_1342[k]
                    + f_3 * pc_y[k] * nsi_1789[k];

        t_2301[k] = f_14 * msi_1510[k]
                    + f_4 * nsh0_1343[k]
                    - f_5 * nsh1_1343[k]
                    + f_3 * pc_y[k] * nsi_1790[k];
    }

#pragma omp simd aligned(t_2302, t_2303, t_2304, pa_y, pc_y, pc_z, msk0_1944, msi_1483, \
                         msi_1511, msk1_1944, nsh0_1343, nsh1_1343, \
                         nsi_1791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2302[k] = f_14 * msi_1511[k]
                    + f_3 * pc_y[k] * nsi_1791[k];

        t_2303[k] = f_21 * msi_1483[k]
                    + f_1 * nsh0_1343[k]
                    - f_2 * nsh1_1343[k]
                    + f_3 * pc_z[k] * nsi_1791[k];

        t_2304[k] = pa_y[k] * msk0_1944[k]
                    - f_12 * pc_y[k] * msk1_1944[k];
    }

#pragma omp simd aligned(t_2305, t_2306, t_2307, pa_y, pc_x, pc_y, msk0_1946, msk1_1946, \
                         nsh0_1345, nsh0_1347, nsh1_1345, nsh1_1347, nsi_1793, \
                         nsi_1795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2305[k] = f_19 * nsh0_1345[k]
                    - f_20 * nsh1_1345[k]
                    + f_3 * pc_x[k] * nsi_1793[k];

        t_2306[k] = pa_y[k] * msk0_1946[k]
                    - f_12 * pc_y[k] * msk1_1946[k];

        t_2307[k] = f_10 * nsh0_1347[k]
                    - f_11 * nsh1_1347[k]
                    + f_3 * pc_x[k] * nsi_1795[k];
    }

#pragma omp simd aligned(t_2308, t_2309, t_2310, pa_y, pc_x, pc_y, msk0_1949, msk1_1949, \
                         nsh0_1348, nsh0_1350, nsh1_1348, nsh1_1350, nsi_1796, \
                         nsi_1798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2308[k] = f_10 * nsh0_1348[k]
                    - f_11 * nsh1_1348[k]
                    + f_3 * pc_x[k] * nsi_1796[k];

        t_2309[k] = pa_y[k] * msk0_1949[k]
                    - f_12 * pc_y[k] * msk1_1949[k];

        t_2310[k] = f_8 * nsh0_1350[k]
                    - f_9 * nsh1_1350[k]
                    + f_3 * pc_x[k] * nsi_1798[k];
    }

#pragma omp simd aligned(t_2311, t_2312, t_2313, pa_y, pc_x, pc_y, msk0_1953, msk1_1953, \
                         nsh0_1351, nsh0_1352, nsh1_1351, nsh1_1352, nsi_1799, \
                         nsi_1800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2311[k] = f_8 * nsh0_1351[k]
                    - f_9 * nsh1_1351[k]
                    + f_3 * pc_x[k] * nsi_1799[k];

        t_2312[k] = f_8 * nsh0_1352[k]
                    - f_9 * nsh1_1352[k]
                    + f_3 * pc_x[k] * nsi_1800[k];

        t_2313[k] = pa_y[k] * msk0_1953[k]
                    - f_12 * pc_y[k] * msk1_1953[k];
    }

#pragma omp simd aligned(t_2314, t_2315, t_2316, pc_x, nsh0_1354, nsh0_1355, nsh0_1356, \
                         nsh1_1354, nsh1_1355, nsh1_1356, nsi_1802, nsi_1803, \
                         nsi_1804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2314[k] = f_6 * nsh0_1354[k]
                    - f_7 * nsh1_1354[k]
                    + f_3 * pc_x[k] * nsi_1802[k];

        t_2315[k] = f_6 * nsh0_1355[k]
                    - f_7 * nsh1_1355[k]
                    + f_3 * pc_x[k] * nsi_1803[k];

        t_2316[k] = f_6 * nsh0_1356[k]
                    - f_7 * nsh1_1356[k]
                    + f_3 * pc_x[k] * nsi_1804[k];
    }

#pragma omp simd aligned(t_2317, t_2318, t_2319, pa_y, pc_x, pc_y, msk0_1958, msk1_1958, \
                         nsh0_1357, nsh0_1359, nsh1_1357, nsh1_1359, nsi_1805, \
                         nsi_1807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2317[k] = f_6 * nsh0_1357[k]
                    - f_7 * nsh1_1357[k]
                    + f_3 * pc_x[k] * nsi_1805[k];

        t_2318[k] = pa_y[k] * msk0_1958[k]
                    - f_12 * pc_y[k] * msk1_1958[k];

        t_2319[k] = f_4 * nsh0_1359[k]
                    - f_5 * nsh1_1359[k]
                    + f_3 * pc_x[k] * nsi_1807[k];
    }

#pragma omp simd aligned(t_2320, t_2321, t_2322, pc_x, nsh0_1360, nsh0_1361, nsh0_1362, \
                         nsh1_1360, nsh1_1361, nsh1_1362, nsi_1808, nsi_1809, \
                         nsi_1810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2320[k] = f_4 * nsh0_1360[k]
                    - f_5 * nsh1_1360[k]
                    + f_3 * pc_x[k] * nsi_1808[k];

        t_2321[k] = f_4 * nsh0_1361[k]
                    - f_5 * nsh1_1361[k]
                    + f_3 * pc_x[k] * nsi_1809[k];

        t_2322[k] = f_4 * nsh0_1362[k]
                    - f_5 * nsh1_1362[k]
                    + f_3 * pc_x[k] * nsi_1810[k];
    }

#pragma omp simd aligned(t_2323, t_2324, t_2325, t_2326, t_2327, pa_y, pc_x, pc_y, msk0_1964, \
                         msk1_1964, nsh0_1363, nsh1_1363, nsi_1811, nsi_1813, nsi_1814, \
                         nsi_1815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2323[k] = f_4 * nsh0_1363[k]
                    - f_5 * nsh1_1363[k]
                    + f_3 * pc_x[k] * nsi_1811[k];

        t_2324[k] = pa_y[k] * msk0_1964[k]
                    - f_12 * pc_y[k] * msk1_1964[k];

        t_2325[k] = f_3 * pc_x[k] * nsi_1813[k];

        t_2326[k] = f_3 * pc_x[k] * nsi_1814[k];

        t_2327[k] = f_3 * pc_x[k] * nsi_1815[k];
    }

#pragma omp simd aligned(t_2328, t_2329, t_2330, t_2331, t_2332, pa_y, pc_x, pc_y, msk0_1972, \
                         msi_1533, msk1_1972, nsi_1816, nsi_1817, nsi_1818, \
                         nsi_1819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2328[k] = f_3 * pc_x[k] * nsi_1816[k];

        t_2329[k] = f_3 * pc_x[k] * nsi_1817[k];

        t_2330[k] = f_3 * pc_x[k] * nsi_1818[k];

        t_2331[k] = f_3 * pc_x[k] * nsi_1819[k];

        t_2332[k] = pa_y[k] * msk0_1972[k]
                    + f_22 * msi_1533[k]
                    - f_12 * pc_y[k] * msk1_1972[k];
    }

#pragma omp simd aligned(t_2333, t_2334, t_2335, pa_y, pc_y, pc_z, msk0_1974, msk0_1975, \
                         msi_1505, msi_1535, msi_1536, msk1_1974, msk1_1975, \
                         nsi_1813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2333[k] = f_18 * msi_1505[k]
                    + f_3 * pc_z[k] * nsi_1813[k];

        t_2334[k] = pa_y[k] * msk0_1974[k]
                    + f_17 * msi_1535[k]
                    - f_12 * pc_y[k] * msk1_1974[k];

        t_2335[k] = pa_y[k] * msk0_1975[k]
                    + f_16 * msi_1536[k]
                    - f_12 * pc_y[k] * msk1_1975[k];
    }

#pragma omp simd aligned(t_2336, t_2337, t_2338, t_2339, pa_y, pc_y, msk0_1976, msk0_1977, \
                         msk0_1979, msi_1537, msi_1538, msi_1539, msk1_1976, msk1_1977, \
                         msk1_1979, nsi_1819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2336[k] = pa_y[k] * msk0_1976[k]
                    + f_15 * msi_1537[k]
                    - f_12 * pc_y[k] * msk1_1976[k];

        t_2337[k] = pa_y[k] * msk0_1977[k]
                    + f_14 * msi_1538[k]
                    - f_12 * pc_y[k] * msk1_1977[k];

        t_2338[k] = f_13 * msi_1539[k]
                    + f_3 * pc_y[k] * nsi_1819[k];

        t_2339[k] = pa_y[k] * msk0_1979[k]
                    - f_12 * pc_y[k] * msk1_1979[k];
    }

#pragma omp simd aligned(t_2340, t_2341, t_2342, t_2343, t_2344, pc_x, pc_y, nsh0_1365, \
                         nsh0_1367, nsh0_1368, nsh1_1365, nsh1_1367, nsh1_1368, nsi_1820, \
                         nsi_1822, nsi_1823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2340[k] = f_1 * nsh0_1365[k]
                    - f_2 * nsh1_1365[k]
                    + f_3 * pc_x[k] * nsi_1820[k];

        t_2341[k] = f_3 * pc_y[k] * nsi_1820[k];

        t_2342[k] = f_19 * nsh0_1367[k]
                    - f_20 * nsh1_1367[k]
                    + f_3 * pc_x[k] * nsi_1822[k];

        t_2343[k] = f_10 * nsh0_1368[k]
                    - f_11 * nsh1_1368[k]
                    + f_3 * pc_x[k] * nsi_1823[k];

        t_2344[k] = f_3 * pc_y[k] * nsi_1822[k];
    }

#pragma omp simd aligned(t_2345, t_2346, t_2347, t_2348, pc_x, pc_y, nsh0_1370, nsh0_1371, \
                         nsh0_1372, nsh1_1370, nsh1_1371, nsh1_1372, nsi_1825, nsi_1826, \
                         nsi_1827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2345[k] = f_10 * nsh0_1370[k]
                    - f_11 * nsh1_1370[k]
                    + f_3 * pc_x[k] * nsi_1825[k];

        t_2346[k] = f_8 * nsh0_1371[k]
                    - f_9 * nsh1_1371[k]
                    + f_3 * pc_x[k] * nsi_1826[k];

        t_2347[k] = f_8 * nsh0_1372[k]
                    - f_9 * nsh1_1372[k]
                    + f_3 * pc_x[k] * nsi_1827[k];

        t_2348[k] = f_3 * pc_y[k] * nsi_1825[k];
    }

#pragma omp simd aligned(t_2349, t_2350, t_2351, pc_x, nsh0_1374, nsh0_1375, nsh0_1376, \
                         nsh1_1374, nsh1_1375, nsh1_1376, nsi_1829, nsi_1830, \
                         nsi_1831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2349[k] = f_8 * nsh0_1374[k]
                    - f_9 * nsh1_1374[k]
                    + f_3 * pc_x[k] * nsi_1829[k];

        t_2350[k] = f_6 * nsh0_1375[k]
                    - f_7 * nsh1_1375[k]
                    + f_3 * pc_x[k] * nsi_1830[k];

        t_2351[k] = f_6 * nsh0_1376[k]
                    - f_7 * nsh1_1376[k]
                    + f_3 * pc_x[k] * nsi_1831[k];
    }

#pragma omp simd aligned(t_2352, t_2353, t_2354, t_2355, pc_x, pc_y, nsh0_1377, nsh0_1379, \
                         nsh0_1380, nsh1_1377, nsh1_1379, nsh1_1380, nsi_1829, nsi_1832, \
                         nsi_1834, nsi_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2352[k] = f_6 * nsh0_1377[k]
                    - f_7 * nsh1_1377[k]
                    + f_3 * pc_x[k] * nsi_1832[k];

        t_2353[k] = f_3 * pc_y[k] * nsi_1829[k];

        t_2354[k] = f_6 * nsh0_1379[k]
                    - f_7 * nsh1_1379[k]
                    + f_3 * pc_x[k] * nsi_1834[k];

        t_2355[k] = f_4 * nsh0_1380[k]
                    - f_5 * nsh1_1380[k]
                    + f_3 * pc_x[k] * nsi_1835[k];
    }

#pragma omp simd aligned(t_2356, t_2357, t_2358, t_2359, pc_x, pc_y, nsh0_1381, nsh0_1382, \
                         nsh0_1383, nsh1_1381, nsh1_1382, nsh1_1383, nsi_1834, nsi_1836, \
                         nsi_1837, nsi_1838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2356[k] = f_4 * nsh0_1381[k]
                    - f_5 * nsh1_1381[k]
                    + f_3 * pc_x[k] * nsi_1836[k];

        t_2357[k] = f_4 * nsh0_1382[k]
                    - f_5 * nsh1_1382[k]
                    + f_3 * pc_x[k] * nsi_1837[k];

        t_2358[k] = f_4 * nsh0_1383[k]
                    - f_5 * nsh1_1383[k]
                    + f_3 * pc_x[k] * nsi_1838[k];

        t_2359[k] = f_3 * pc_y[k] * nsi_1834[k];
    }

#pragma omp simd aligned(t_2360, t_2361, t_2362, t_2363, t_2364, t_2365, pc_x, nsh0_1385, \
                         nsh1_1385, nsi_1840, nsi_1841, nsi_1842, nsi_1843, nsi_1844, \
                         nsi_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2360[k] = f_4 * nsh0_1385[k]
                    - f_5 * nsh1_1385[k]
                    + f_3 * pc_x[k] * nsi_1840[k];

        t_2361[k] = f_3 * pc_x[k] * nsi_1841[k];

        t_2362[k] = f_3 * pc_x[k] * nsi_1842[k];

        t_2363[k] = f_3 * pc_x[k] * nsi_1843[k];

        t_2364[k] = f_3 * pc_x[k] * nsi_1844[k];

        t_2365[k] = f_3 * pc_x[k] * nsi_1845[k];
    }

#pragma omp simd aligned(t_2366, t_2367, t_2368, t_2369, pc_x, pc_y, nsh0_1380, nsh0_1381, \
                         nsh1_1380, nsh1_1381, nsi_1841, nsi_1842, nsi_1846, \
                         nsi_1847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2366[k] = f_3 * pc_x[k] * nsi_1846[k];

        t_2367[k] = f_3 * pc_x[k] * nsi_1847[k];

        t_2368[k] = f_1 * nsh0_1380[k]
                    - f_2 * nsh1_1380[k]
                    + f_3 * pc_y[k] * nsi_1841[k];

        t_2369[k] = f_19 * nsh0_1381[k]
                    - f_20 * nsh1_1381[k]
                    + f_3 * pc_y[k] * nsi_1842[k];
    }

#pragma omp simd aligned(t_2370, t_2371, t_2372, pc_y, nsh0_1382, nsh0_1383, nsh0_1384, \
                         nsh1_1382, nsh1_1383, nsh1_1384, nsi_1843, nsi_1844, \
                         nsi_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2370[k] = f_10 * nsh0_1382[k]
                    - f_11 * nsh1_1382[k]
                    + f_3 * pc_y[k] * nsi_1843[k];

        t_2371[k] = f_8 * nsh0_1383[k]
                    - f_9 * nsh1_1383[k]
                    + f_3 * pc_y[k] * nsi_1844[k];

        t_2372[k] = f_6 * nsh0_1384[k]
                    - f_7 * nsh1_1384[k]
                    + f_3 * pc_y[k] * nsi_1845[k];
    }

#pragma omp simd aligned(t_2373, t_2374, t_2375, pc_y, pc_z, msi_1539, nsh0_1385, nsh1_1385, \
                         nsi_1846, nsi_1847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2373[k] = f_4 * nsh0_1385[k]
                    - f_5 * nsh1_1385[k]
                    + f_3 * pc_y[k] * nsi_1846[k];

        t_2374[k] = f_3 * pc_y[k] * nsi_1847[k];

        t_2375[k] = f_0 * msi_1539[k]
                    + f_1 * nsh0_1385[k]
                    - f_2 * nsh1_1385[k]
                    + f_3 * pc_z[k] * nsi_1847[k];
    }
}

auto
compute_prim_nsk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t msk0, const size_t msi,
                                                   const size_t msk1, const size_t nsh0,
                                                   const size_t nsh1, const size_t nsi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_nsk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece7(buffer, target, pc, msi, nsh0,
                                                              nsh1, nsi, ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, msk0, msi,
                                                              msk1, nsh0, nsh1, nsi, ncols,
                                                              gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsh0, nsh1, nsi,
                                                               ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsh0, nsh1, nsi,
                                                               ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece12(buffer, target, pc, msi, nsh0,
                                                               nsh1, nsi, ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsh0, nsh1, nsi,
                                                               ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece14(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsh0, nsh1, nsi,
                                                               ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece15(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsi, ncols, gamma, p,
                                                               q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsi, ncols, gamma, p,
                                                               q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece17(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsh0, nsh1, nsi,
                                                               ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece18(buffer, target, pc, msi, nsh0,
                                                               nsh1, nsi, ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece19(buffer, target, pc, msi, nsh0,
                                                               nsh1, nsi, ncols, gamma, p, q);

    compute_prim_nsk_three_center_electron_repulsion_0_piece20(buffer, target, pa, pc, msk0,
                                                               msi, msk1, nsh0, nsh1, nsi,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
