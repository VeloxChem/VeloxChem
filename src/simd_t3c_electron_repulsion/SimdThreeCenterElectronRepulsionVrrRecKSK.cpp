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


#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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

    const auto *isk0_0 = buffer.data(isk0 + 0);
    const auto *isk0_3 = buffer.data(isk0 + 3);
    const auto *isk0_5 = buffer.data(isk0 + 5);
    const auto *isk0_6 = buffer.data(isk0 + 6);
    const auto *isk0_9 = buffer.data(isk0 + 9);
    const auto *isk0_10 = buffer.data(isk0 + 10);
    const auto *isk0_14 = buffer.data(isk0 + 14);
    const auto *isk0_15 = buffer.data(isk0 + 15);
    const auto *isk0_20 = buffer.data(isk0 + 20);
    const auto *isk0_28 = buffer.data(isk0 + 28);
    const auto *isk0_35 = buffer.data(isk0 + 35);

    const auto *isi_0 = buffer.data(isi + 0);
    const auto *isi_1 = buffer.data(isi + 1);
    const auto *isi_2 = buffer.data(isi + 2);
    const auto *isi_3 = buffer.data(isi + 3);
    const auto *isi_5 = buffer.data(isi + 5);
    const auto *isi_6 = buffer.data(isi + 6);
    const auto *isi_9 = buffer.data(isi + 9);
    const auto *isi_10 = buffer.data(isi + 10);
    const auto *isi_14 = buffer.data(isi + 14);
    const auto *isi_21 = buffer.data(isi + 21);
    const auto *isi_23 = buffer.data(isi + 23);
    const auto *isi_24 = buffer.data(isi + 24);
    const auto *isi_25 = buffer.data(isi + 25);
    const auto *isi_27 = buffer.data(isi + 27);
    const auto *isi_28 = buffer.data(isi + 28);
    const auto *isi_33 = buffer.data(isi + 33);
    const auto *isi_37 = buffer.data(isi + 37);
    const auto *isi_42 = buffer.data(isi + 42);
    const auto *isi_49 = buffer.data(isi + 49);
    const auto *isi_51 = buffer.data(isi + 51);
    const auto *isi_52 = buffer.data(isi + 52);
    const auto *isi_53 = buffer.data(isi + 53);
    const auto *isi_54 = buffer.data(isi + 54);
    const auto *isi_55 = buffer.data(isi + 55);
    const auto *isi_77 = buffer.data(isi + 77);
    const auto *isi_78 = buffer.data(isi + 78);
    const auto *isi_79 = buffer.data(isi + 79);
    const auto *isi_80 = buffer.data(isi + 80);
    const auto *isi_81 = buffer.data(isi + 81);
    const auto *isi_83 = buffer.data(isi + 83);
    const auto *isi_84 = buffer.data(isi + 84);
    const auto *isi_87 = buffer.data(isi + 87);
    const auto *isi_90 = buffer.data(isi + 90);
    const auto *isi_94 = buffer.data(isi + 94);
    const auto *isi_99 = buffer.data(isi + 99);

    const auto *isk1_0 = buffer.data(isk1 + 0);
    const auto *isk1_3 = buffer.data(isk1 + 3);
    const auto *isk1_5 = buffer.data(isk1 + 5);
    const auto *isk1_6 = buffer.data(isk1 + 6);
    const auto *isk1_9 = buffer.data(isk1 + 9);
    const auto *isk1_10 = buffer.data(isk1 + 10);
    const auto *isk1_14 = buffer.data(isk1 + 14);
    const auto *isk1_15 = buffer.data(isk1 + 15);
    const auto *isk1_20 = buffer.data(isk1 + 20);
    const auto *isk1_28 = buffer.data(isk1 + 28);
    const auto *isk1_35 = buffer.data(isk1 + 35);

    const auto *ksh0_0 = buffer.data(ksh0 + 0);
    const auto *ksh0_1 = buffer.data(ksh0 + 1);
    const auto *ksh0_2 = buffer.data(ksh0 + 2);
    const auto *ksh0_3 = buffer.data(ksh0 + 3);
    const auto *ksh0_5 = buffer.data(ksh0 + 5);
    const auto *ksh0_6 = buffer.data(ksh0 + 6);
    const auto *ksh0_8 = buffer.data(ksh0 + 8);
    const auto *ksh0_9 = buffer.data(ksh0 + 9);
    const auto *ksh0_15 = buffer.data(ksh0 + 15);
    const auto *ksh0_17 = buffer.data(ksh0 + 17);
    const auto *ksh0_18 = buffer.data(ksh0 + 18);
    const auto *ksh0_19 = buffer.data(ksh0 + 19);
    const auto *ksh0_20 = buffer.data(ksh0 + 20);
    const auto *ksh0_24 = buffer.data(ksh0 + 24);
    const auto *ksh0_27 = buffer.data(ksh0 + 27);
    const auto *ksh0_28 = buffer.data(ksh0 + 28);
    const auto *ksh0_36 = buffer.data(ksh0 + 36);
    const auto *ksh0_37 = buffer.data(ksh0 + 37);
    const auto *ksh0_38 = buffer.data(ksh0 + 38);
    const auto *ksh0_39 = buffer.data(ksh0 + 39);
    const auto *ksh0_44 = buffer.data(ksh0 + 44);
    const auto *ksh0_46 = buffer.data(ksh0 + 46);
    const auto *ksh0_47 = buffer.data(ksh0 + 47);
    const auto *ksh0_49 = buffer.data(ksh0 + 49);
    const auto *ksh0_50 = buffer.data(ksh0 + 50);
    const auto *ksh0_51 = buffer.data(ksh0 + 51);
    const auto *ksh0_58 = buffer.data(ksh0 + 58);
    const auto *ksh0_59 = buffer.data(ksh0 + 59);
    const auto *ksh0_60 = buffer.data(ksh0 + 60);
    const auto *ksh0_61 = buffer.data(ksh0 + 61);
    const auto *ksh0_62 = buffer.data(ksh0 + 62);
    const auto *ksh0_63 = buffer.data(ksh0 + 63);
    const auto *ksh0_65 = buffer.data(ksh0 + 65);
    const auto *ksh0_66 = buffer.data(ksh0 + 66);
    const auto *ksh0_68 = buffer.data(ksh0 + 68);
    const auto *ksh0_69 = buffer.data(ksh0 + 69);
    const auto *ksh0_70 = buffer.data(ksh0 + 70);
    const auto *ksh0_72 = buffer.data(ksh0 + 72);
    const auto *ksh0_73 = buffer.data(ksh0 + 73);
    const auto *ksh0_78 = buffer.data(ksh0 + 78);

    const auto *ksh1_0 = buffer.data(ksh1 + 0);
    const auto *ksh1_1 = buffer.data(ksh1 + 1);
    const auto *ksh1_2 = buffer.data(ksh1 + 2);
    const auto *ksh1_3 = buffer.data(ksh1 + 3);
    const auto *ksh1_5 = buffer.data(ksh1 + 5);
    const auto *ksh1_6 = buffer.data(ksh1 + 6);
    const auto *ksh1_8 = buffer.data(ksh1 + 8);
    const auto *ksh1_9 = buffer.data(ksh1 + 9);
    const auto *ksh1_15 = buffer.data(ksh1 + 15);
    const auto *ksh1_17 = buffer.data(ksh1 + 17);
    const auto *ksh1_18 = buffer.data(ksh1 + 18);
    const auto *ksh1_19 = buffer.data(ksh1 + 19);
    const auto *ksh1_20 = buffer.data(ksh1 + 20);
    const auto *ksh1_24 = buffer.data(ksh1 + 24);
    const auto *ksh1_27 = buffer.data(ksh1 + 27);
    const auto *ksh1_28 = buffer.data(ksh1 + 28);
    const auto *ksh1_36 = buffer.data(ksh1 + 36);
    const auto *ksh1_37 = buffer.data(ksh1 + 37);
    const auto *ksh1_38 = buffer.data(ksh1 + 38);
    const auto *ksh1_39 = buffer.data(ksh1 + 39);
    const auto *ksh1_44 = buffer.data(ksh1 + 44);
    const auto *ksh1_46 = buffer.data(ksh1 + 46);
    const auto *ksh1_47 = buffer.data(ksh1 + 47);
    const auto *ksh1_49 = buffer.data(ksh1 + 49);
    const auto *ksh1_50 = buffer.data(ksh1 + 50);
    const auto *ksh1_51 = buffer.data(ksh1 + 51);
    const auto *ksh1_58 = buffer.data(ksh1 + 58);
    const auto *ksh1_59 = buffer.data(ksh1 + 59);
    const auto *ksh1_60 = buffer.data(ksh1 + 60);
    const auto *ksh1_61 = buffer.data(ksh1 + 61);
    const auto *ksh1_62 = buffer.data(ksh1 + 62);
    const auto *ksh1_63 = buffer.data(ksh1 + 63);
    const auto *ksh1_65 = buffer.data(ksh1 + 65);
    const auto *ksh1_66 = buffer.data(ksh1 + 66);
    const auto *ksh1_68 = buffer.data(ksh1 + 68);
    const auto *ksh1_69 = buffer.data(ksh1 + 69);
    const auto *ksh1_70 = buffer.data(ksh1 + 70);
    const auto *ksh1_72 = buffer.data(ksh1 + 72);
    const auto *ksh1_73 = buffer.data(ksh1 + 73);
    const auto *ksh1_78 = buffer.data(ksh1 + 78);

    const auto *ksi_0 = buffer.data(ksi + 0);
    const auto *ksi_1 = buffer.data(ksi + 1);
    const auto *ksi_2 = buffer.data(ksi + 2);
    const auto *ksi_3 = buffer.data(ksi + 3);
    const auto *ksi_5 = buffer.data(ksi + 5);
    const auto *ksi_6 = buffer.data(ksi + 6);
    const auto *ksi_8 = buffer.data(ksi + 8);
    const auto *ksi_9 = buffer.data(ksi + 9);
    const auto *ksi_10 = buffer.data(ksi + 10);
    const auto *ksi_12 = buffer.data(ksi + 12);
    const auto *ksi_13 = buffer.data(ksi + 13);
    const auto *ksi_14 = buffer.data(ksi + 14);
    const auto *ksi_15 = buffer.data(ksi + 15);
    const auto *ksi_20 = buffer.data(ksi + 20);
    const auto *ksi_21 = buffer.data(ksi + 21);
    const auto *ksi_23 = buffer.data(ksi + 23);
    const auto *ksi_24 = buffer.data(ksi + 24);
    const auto *ksi_25 = buffer.data(ksi + 25);
    const auto *ksi_26 = buffer.data(ksi + 26);
    const auto *ksi_27 = buffer.data(ksi + 27);
    const auto *ksi_28 = buffer.data(ksi + 28);
    const auto *ksi_29 = buffer.data(ksi + 29);
    const auto *ksi_31 = buffer.data(ksi + 31);
    const auto *ksi_33 = buffer.data(ksi + 33);
    const auto *ksi_34 = buffer.data(ksi + 34);
    const auto *ksi_35 = buffer.data(ksi + 35);
    const auto *ksi_37 = buffer.data(ksi + 37);
    const auto *ksi_38 = buffer.data(ksi + 38);
    const auto *ksi_39 = buffer.data(ksi + 39);
    const auto *ksi_40 = buffer.data(ksi + 40);
    const auto *ksi_42 = buffer.data(ksi + 42);
    const auto *ksi_43 = buffer.data(ksi + 43);
    const auto *ksi_49 = buffer.data(ksi + 49);
    const auto *ksi_50 = buffer.data(ksi + 50);
    const auto *ksi_51 = buffer.data(ksi + 51);
    const auto *ksi_52 = buffer.data(ksi + 52);
    const auto *ksi_53 = buffer.data(ksi + 53);
    const auto *ksi_54 = buffer.data(ksi + 54);
    const auto *ksi_55 = buffer.data(ksi + 55);
    const auto *ksi_56 = buffer.data(ksi + 56);
    const auto *ksi_58 = buffer.data(ksi + 58);
    const auto *ksi_60 = buffer.data(ksi + 60);
    const auto *ksi_61 = buffer.data(ksi + 61);
    const auto *ksi_63 = buffer.data(ksi + 63);
    const auto *ksi_64 = buffer.data(ksi + 64);
    const auto *ksi_65 = buffer.data(ksi + 65);
    const auto *ksi_67 = buffer.data(ksi + 67);
    const auto *ksi_68 = buffer.data(ksi + 68);
    const auto *ksi_69 = buffer.data(ksi + 69);
    const auto *ksi_70 = buffer.data(ksi + 70);
    const auto *ksi_76 = buffer.data(ksi + 76);
    const auto *ksi_77 = buffer.data(ksi + 77);
    const auto *ksi_78 = buffer.data(ksi + 78);
    const auto *ksi_79 = buffer.data(ksi + 79);
    const auto *ksi_80 = buffer.data(ksi + 80);
    const auto *ksi_81 = buffer.data(ksi + 81);
    const auto *ksi_82 = buffer.data(ksi + 82);
    const auto *ksi_83 = buffer.data(ksi + 83);
    const auto *ksi_84 = buffer.data(ksi + 84);
    const auto *ksi_85 = buffer.data(ksi + 85);
    const auto *ksi_86 = buffer.data(ksi + 86);
    const auto *ksi_87 = buffer.data(ksi + 87);
    const auto *ksi_89 = buffer.data(ksi + 89);
    const auto *ksi_90 = buffer.data(ksi + 90);
    const auto *ksi_91 = buffer.data(ksi + 91);
    const auto *ksi_93 = buffer.data(ksi + 93);
    const auto *ksi_94 = buffer.data(ksi + 94);
    const auto *ksi_95 = buffer.data(ksi + 95);
    const auto *ksi_96 = buffer.data(ksi + 96);
    const auto *ksi_98 = buffer.data(ksi + 98);
    const auto *ksi_99 = buffer.data(ksi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, isi_0, ksh0_0, \
                         ksh1_0, ksi_0, ksi_1, ksi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * isi_0[k]
                 + f_1 * ksh0_0[k]
                 - f_2 * ksh1_0[k]
                 + f_3 * pc_x[k] * ksi_0[k];

        t_1[k] = f_3 * pc_y[k] * ksi_0[k];

        t_2[k] = f_3 * pc_z[k] * ksi_0[k];

        t_3[k] = f_4 * ksh0_0[k]
                 - f_5 * ksh1_0[k]
                 + f_3 * pc_y[k] * ksi_1[k];

        t_4[k] = f_3 * pc_y[k] * ksi_2[k];

        t_5[k] = f_4 * ksh0_0[k]
                 - f_5 * ksh1_0[k]
                 + f_3 * pc_z[k] * ksi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, ksh0_1, ksh0_2, ksh0_3, ksh1_1, \
                         ksh1_2, ksh1_3, ksi_3, ksi_5, ksi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ksh0_1[k]
                 - f_7 * ksh1_1[k]
                 + f_3 * pc_y[k] * ksi_3[k];

        t_7[k] = f_3 * pc_z[k] * ksi_3[k];

        t_8[k] = f_3 * pc_y[k] * ksi_5[k];

        t_9[k] = f_6 * ksh0_2[k]
                 - f_7 * ksh1_2[k]
                 + f_3 * pc_z[k] * ksi_5[k];

        t_10[k] = f_8 * ksh0_3[k]
                  - f_9 * ksh1_3[k]
                  + f_3 * pc_y[k] * ksi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, ksh0_5, ksh0_6, \
                         ksh1_5, ksh1_6, ksi_6, ksi_8, ksi_9, ksi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * ksi_6[k];

        t_12[k] = f_4 * ksh0_5[k]
                  - f_5 * ksh1_5[k]
                  + f_3 * pc_y[k] * ksi_8[k];

        t_13[k] = f_3 * pc_y[k] * ksi_9[k];

        t_14[k] = f_8 * ksh0_5[k]
                  - f_9 * ksh1_5[k]
                  + f_3 * pc_z[k] * ksi_9[k];

        t_15[k] = f_10 * ksh0_6[k]
                  - f_11 * ksh1_6[k]
                  + f_3 * pc_y[k] * ksi_10[k];

        t_16[k] = f_3 * pc_z[k] * ksi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, ksh0_8, ksh0_9, ksh1_8, ksh1_9, \
                         ksi_12, ksi_13, ksi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * ksh0_8[k]
                  - f_7 * ksh1_8[k]
                  + f_3 * pc_y[k] * ksi_12[k];

        t_18[k] = f_4 * ksh0_9[k]
                  - f_5 * ksh1_9[k]
                  + f_3 * pc_y[k] * ksi_13[k];

        t_19[k] = f_3 * pc_y[k] * ksi_14[k];

        t_20[k] = f_10 * ksh0_9[k]
                  - f_11 * ksh1_9[k]
                  + f_3 * pc_z[k] * ksi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, isi_21, isi_23, isi_24, \
                         isi_25, ksi_15, ksi_21, ksi_23, ksi_24, \
                         ksi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * isi_21[k]
                  + f_3 * pc_x[k] * ksi_21[k];

        t_22[k] = f_3 * pc_z[k] * ksi_15[k];

        t_23[k] = f_0 * isi_23[k]
                  + f_3 * pc_x[k] * ksi_23[k];

        t_24[k] = f_0 * isi_24[k]
                  + f_3 * pc_x[k] * ksi_24[k];

        t_25[k] = f_0 * isi_25[k]
                  + f_3 * pc_x[k] * ksi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, isi_27, ksh0_15, ksh1_15, \
                         ksi_20, ksi_21, ksi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * ksi_20[k];

        t_27[k] = f_0 * isi_27[k]
                  + f_3 * pc_x[k] * ksi_27[k];

        t_28[k] = f_1 * ksh0_15[k]
                  - f_2 * ksh1_15[k]
                  + f_3 * pc_y[k] * ksi_21[k];

        t_29[k] = f_3 * pc_z[k] * ksi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, ksh0_17, ksh0_18, ksh0_19, ksh1_17, ksh1_18, \
                         ksh1_19, ksi_23, ksi_24, ksi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * ksh0_17[k]
                  - f_11 * ksh1_17[k]
                  + f_3 * pc_y[k] * ksi_23[k];

        t_31[k] = f_8 * ksh0_18[k]
                  - f_9 * ksh1_18[k]
                  + f_3 * pc_y[k] * ksi_24[k];

        t_32[k] = f_6 * ksh0_19[k]
                  - f_7 * ksh1_19[k]
                  + f_3 * pc_y[k] * ksi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, isk0_0, isi_0, \
                         isk1_0, ksh0_20, ksh1_20, ksi_26, ksi_27, \
                         ksi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * ksh0_20[k]
                  - f_5 * ksh1_20[k]
                  + f_3 * pc_y[k] * ksi_26[k];

        t_34[k] = f_3 * pc_y[k] * ksi_27[k];

        t_35[k] = f_1 * ksh0_20[k]
                  - f_2 * ksh1_20[k]
                  + f_3 * pc_z[k] * ksi_27[k];

        t_36[k] = pa_y[k] * isk0_0[k]
                  - f_12 * pc_y[k] * isk1_0[k];

        t_37[k] = f_13 * isi_0[k]
                  + f_3 * pc_y[k] * ksi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, isk0_3, isk0_5, isi_1, \
                         isk1_3, isk1_5, ksi_28, ksi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * ksi_28[k];

        t_39[k] = pa_y[k] * isk0_3[k]
                  + f_14 * isi_1[k]
                  - f_12 * pc_y[k] * isk1_3[k];

        t_40[k] = f_3 * pc_z[k] * ksi_29[k];

        t_41[k] = pa_y[k] * isk0_5[k]
                  - f_12 * pc_y[k] * isk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, isk0_6, isk0_9, isi_3, \
                         isi_5, isk1_6, isk1_9, ksi_31, ksi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * isk0_6[k]
                  + f_15 * isi_3[k]
                  - f_12 * pc_y[k] * isk1_6[k];

        t_43[k] = f_3 * pc_z[k] * ksi_31[k];

        t_44[k] = f_13 * isi_5[k]
                  + f_3 * pc_y[k] * ksi_33[k];

        t_45[k] = pa_y[k] * isk0_9[k]
                  - f_12 * pc_y[k] * isk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, isk0_10, isi_6, isi_9, \
                         isk1_10, ksh0_24, ksh1_24, ksi_34, ksi_35, \
                         ksi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * isk0_10[k]
                  + f_16 * isi_6[k]
                  - f_12 * pc_y[k] * isk1_10[k];

        t_47[k] = f_3 * pc_z[k] * ksi_34[k];

        t_48[k] = f_4 * ksh0_24[k]
                  - f_5 * ksh1_24[k]
                  + f_3 * pc_z[k] * ksi_35[k];

        t_49[k] = f_13 * isi_9[k]
                  + f_3 * pc_y[k] * ksi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, isk0_14, isk0_15, isi_10, \
                         isk1_14, isk1_15, ksh0_27, ksh1_27, ksi_38, \
                         ksi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * isk0_14[k]
                  - f_12 * pc_y[k] * isk1_14[k];

        t_51[k] = pa_y[k] * isk0_15[k]
                  + f_17 * isi_10[k]
                  - f_12 * pc_y[k] * isk1_15[k];

        t_52[k] = f_3 * pc_z[k] * ksi_38[k];

        t_53[k] = f_4 * ksh0_27[k]
                  - f_5 * ksh1_27[k]
                  + f_3 * pc_z[k] * ksi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, isk0_20, isi_14, isk1_20, \
                         ksh0_28, ksh1_28, ksi_40, ksi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * ksh0_28[k]
                  - f_7 * ksh1_28[k]
                  + f_3 * pc_z[k] * ksi_40[k];

        t_55[k] = f_13 * isi_14[k]
                  + f_3 * pc_y[k] * ksi_42[k];

        t_56[k] = pa_y[k] * isk0_20[k]
                  - f_12 * pc_y[k] * isk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, isi_49, isi_51, isi_52, \
                         isi_53, ksi_43, ksi_49, ksi_51, ksi_52, \
                         ksi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_18 * isi_49[k]
                  + f_3 * pc_x[k] * ksi_49[k];

        t_58[k] = f_3 * pc_z[k] * ksi_43[k];

        t_59[k] = f_18 * isi_51[k]
                  + f_3 * pc_x[k] * ksi_51[k];

        t_60[k] = f_18 * isi_52[k]
                  + f_3 * pc_x[k] * ksi_52[k];

        t_61[k] = f_18 * isi_53[k]
                  + f_3 * pc_x[k] * ksi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, isi_21, isi_54, isi_55, \
                         ksh0_36, ksh1_36, ksi_49, ksi_54, ksi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_18 * isi_54[k]
                  + f_3 * pc_x[k] * ksi_54[k];

        t_63[k] = f_18 * isi_55[k]
                  + f_3 * pc_x[k] * ksi_55[k];

        t_64[k] = f_13 * isi_21[k]
                  + f_1 * ksh0_36[k]
                  - f_2 * ksh1_36[k]
                  + f_3 * pc_y[k] * ksi_49[k];

        t_65[k] = f_3 * pc_z[k] * ksi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, ksh0_36, ksh0_37, ksh0_38, ksh1_36, ksh1_37, \
                         ksh1_38, ksi_50, ksi_51, ksi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * ksh0_36[k]
                  - f_5 * ksh1_36[k]
                  + f_3 * pc_z[k] * ksi_50[k];

        t_67[k] = f_6 * ksh0_37[k]
                  - f_7 * ksh1_37[k]
                  + f_3 * pc_z[k] * ksi_51[k];

        t_68[k] = f_8 * ksh0_38[k]
                  - f_9 * ksh1_38[k]
                  + f_3 * pc_z[k] * ksi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, isk0_35, isi_27, isk1_35, \
                         ksh0_39, ksh1_39, ksi_53, ksi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ksh0_39[k]
                  - f_11 * ksh1_39[k]
                  + f_3 * pc_z[k] * ksi_53[k];

        t_70[k] = f_13 * isi_27[k]
                  + f_3 * pc_y[k] * ksi_55[k];

        t_71[k] = pa_y[k] * isk0_35[k]
                  - f_12 * pc_y[k] * isk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, isk0_0, isk0_3, \
                         isi_0, isk1_0, isk1_3, ksi_56, ksi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * isk0_0[k]
                  - f_12 * pc_z[k] * isk1_0[k];

        t_73[k] = f_3 * pc_y[k] * ksi_56[k];

        t_74[k] = f_13 * isi_0[k]
                  + f_3 * pc_z[k] * ksi_56[k];

        t_75[k] = pa_z[k] * isk0_3[k]
                  - f_12 * pc_z[k] * isk1_3[k];

        t_76[k] = f_3 * pc_y[k] * ksi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, isk0_5, isk0_6, isi_2, \
                         isk1_5, isk1_6, ksh0_44, ksh1_44, ksi_60, \
                         ksi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * isk0_5[k]
                  + f_14 * isi_2[k]
                  - f_12 * pc_z[k] * isk1_5[k];

        t_78[k] = pa_z[k] * isk0_6[k]
                  - f_12 * pc_z[k] * isk1_6[k];

        t_79[k] = f_4 * ksh0_44[k]
                  - f_5 * ksh1_44[k]
                  + f_3 * pc_y[k] * ksi_60[k];

        t_80[k] = f_3 * pc_y[k] * ksi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, isk0_9, isk0_10, isi_5, isk1_9, \
                         isk1_10, ksh0_46, ksh1_46, ksi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * isk0_9[k]
                  + f_15 * isi_5[k]
                  - f_12 * pc_z[k] * isk1_9[k];

        t_82[k] = pa_z[k] * isk0_10[k]
                  - f_12 * pc_z[k] * isk1_10[k];

        t_83[k] = f_6 * ksh0_46[k]
                  - f_7 * ksh1_46[k]
                  + f_3 * pc_y[k] * ksi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, isk0_14, isk0_15, isi_9, \
                         isk1_14, isk1_15, ksh0_47, ksh1_47, ksi_64, \
                         ksi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * ksh0_47[k]
                  - f_5 * ksh1_47[k]
                  + f_3 * pc_y[k] * ksi_64[k];

        t_85[k] = f_3 * pc_y[k] * ksi_65[k];

        t_86[k] = pa_z[k] * isk0_14[k]
                  + f_16 * isi_9[k]
                  - f_12 * pc_z[k] * isk1_14[k];

        t_87[k] = pa_z[k] * isk0_15[k]
                  - f_12 * pc_z[k] * isk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, ksh0_49, ksh0_50, ksh0_51, ksh1_49, \
                         ksh1_50, ksh1_51, ksi_67, ksi_68, ksi_69, \
                         ksi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * ksh0_49[k]
                  - f_9 * ksh1_49[k]
                  + f_3 * pc_y[k] * ksi_67[k];

        t_89[k] = f_6 * ksh0_50[k]
                  - f_7 * ksh1_50[k]
                  + f_3 * pc_y[k] * ksi_68[k];

        t_90[k] = f_4 * ksh0_51[k]
                  - f_5 * ksh1_51[k]
                  + f_3 * pc_y[k] * ksi_69[k];

        t_91[k] = f_3 * pc_y[k] * ksi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, isk0_20, isi_14, isi_77, \
                         isi_78, isi_79, isk1_20, ksi_77, ksi_78, \
                         ksi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * isk0_20[k]
                  + f_17 * isi_14[k]
                  - f_12 * pc_z[k] * isk1_20[k];

        t_93[k] = f_18 * isi_77[k]
                  + f_3 * pc_x[k] * ksi_77[k];

        t_94[k] = f_18 * isi_78[k]
                  + f_3 * pc_x[k] * ksi_78[k];

        t_95[k] = f_18 * isi_79[k]
                  + f_3 * pc_x[k] * ksi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, isi_80, isi_81, isi_83, ksi_76, \
                         ksi_80, ksi_81, ksi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_18 * isi_80[k]
                  + f_3 * pc_x[k] * ksi_80[k];

        t_97[k] = f_18 * isi_81[k]
                  + f_3 * pc_x[k] * ksi_81[k];

        t_98[k] = f_3 * pc_y[k] * ksi_76[k];

        t_99[k] = f_18 * isi_83[k]
                  + f_3 * pc_x[k] * ksi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, isk0_28, isk1_28, ksh0_58, \
                         ksh0_59, ksh1_58, ksh1_59, ksi_78, ksi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * isk0_28[k]
                   - f_12 * pc_z[k] * isk1_28[k];

        t_101[k] = f_19 * ksh0_58[k]
                   - f_20 * ksh1_58[k]
                   + f_3 * pc_y[k] * ksi_78[k];

        t_102[k] = f_10 * ksh0_59[k]
                   - f_11 * ksh1_59[k]
                   + f_3 * pc_y[k] * ksi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, ksh0_60, ksh0_61, ksh0_62, ksh1_60, \
                         ksh1_61, ksh1_62, ksi_80, ksi_81, ksi_82, \
                         ksi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * ksh0_60[k]
                   - f_9 * ksh1_60[k]
                   + f_3 * pc_y[k] * ksi_80[k];

        t_104[k] = f_6 * ksh0_61[k]
                   - f_7 * ksh1_61[k]
                   + f_3 * pc_y[k] * ksi_81[k];

        t_105[k] = f_4 * ksh0_62[k]
                   - f_5 * ksh1_62[k]
                   + f_3 * pc_y[k] * ksi_82[k];

        t_106[k] = f_3 * pc_y[k] * ksi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, isi_27, isi_28, isi_84, \
                         ksh0_62, ksh0_63, ksh1_62, ksh1_63, ksi_83, \
                         ksi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * isi_27[k]
                   + f_1 * ksh0_62[k]
                   - f_2 * ksh1_62[k]
                   + f_3 * pc_z[k] * ksi_83[k];

        t_108[k] = f_17 * isi_84[k]
                   + f_1 * ksh0_63[k]
                   - f_2 * ksh1_63[k]
                   + f_3 * pc_x[k] * ksi_84[k];

        t_109[k] = f_14 * isi_28[k]
                   + f_3 * pc_y[k] * ksi_84[k];

        t_110[k] = f_3 * pc_z[k] * ksi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, isi_87, ksh0_63, ksh0_66, ksh1_63, \
                         ksh1_66, ksi_85, ksi_86, ksi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_17 * isi_87[k]
                   + f_10 * ksh0_66[k]
                   - f_11 * ksh1_66[k]
                   + f_3 * pc_x[k] * ksi_87[k];

        t_112[k] = f_3 * pc_z[k] * ksi_85[k];

        t_113[k] = f_4 * ksh0_63[k]
                   - f_5 * ksh1_63[k]
                   + f_3 * pc_z[k] * ksi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, isi_33, isi_90, \
                         ksh0_65, ksh0_69, ksh1_65, ksh1_69, ksi_87, ksi_89, \
                         ksi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_17 * isi_90[k]
                   + f_8 * ksh0_69[k]
                   - f_9 * ksh1_69[k]
                   + f_3 * pc_x[k] * ksi_90[k];

        t_115[k] = f_3 * pc_z[k] * ksi_87[k];

        t_116[k] = f_14 * isi_33[k]
                   + f_3 * pc_y[k] * ksi_89[k];

        t_117[k] = f_6 * ksh0_65[k]
                   - f_7 * ksh1_65[k]
                   + f_3 * pc_z[k] * ksi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, isi_94, ksh0_66, ksh0_73, ksh1_66, \
                         ksh1_73, ksi_90, ksi_91, ksi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_17 * isi_94[k]
                   + f_6 * ksh0_73[k]
                   - f_7 * ksh1_73[k]
                   + f_3 * pc_x[k] * ksi_94[k];

        t_119[k] = f_3 * pc_z[k] * ksi_90[k];

        t_120[k] = f_4 * ksh0_66[k]
                   - f_5 * ksh1_66[k]
                   + f_3 * pc_z[k] * ksi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, isi_37, isi_99, \
                         ksh0_68, ksh0_78, ksh1_68, ksh1_78, ksi_93, ksi_94, \
                         ksi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * isi_37[k]
                   + f_3 * pc_y[k] * ksi_93[k];

        t_122[k] = f_8 * ksh0_68[k]
                   - f_9 * ksh1_68[k]
                   + f_3 * pc_z[k] * ksi_93[k];

        t_123[k] = f_17 * isi_99[k]
                   + f_4 * ksh0_78[k]
                   - f_5 * ksh1_78[k]
                   + f_3 * pc_x[k] * ksi_99[k];

        t_124[k] = f_3 * pc_z[k] * ksi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, isi_42, ksh0_69, ksh0_70, \
                         ksh0_72, ksh1_69, ksh1_70, ksh1_72, ksi_95, ksi_96, \
                         ksi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * ksh0_69[k]
                   - f_5 * ksh1_69[k]
                   + f_3 * pc_z[k] * ksi_95[k];

        t_126[k] = f_6 * ksh0_70[k]
                   - f_7 * ksh1_70[k]
                   + f_3 * pc_z[k] * ksi_96[k];

        t_127[k] = f_14 * isi_42[k]
                   + f_3 * pc_y[k] * ksi_98[k];

        t_128[k] = f_10 * ksh0_72[k]
                   - f_11 * ksh1_72[k]
                   + f_3 * pc_z[k] * ksi_98[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
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

    const auto *isk0_39 = buffer.data(isk0 + 39);
    const auto *isk0_42 = buffer.data(isk0 + 42);
    const auto *isk0_46 = buffer.data(isk0 + 46);
    const auto *isk0_51 = buffer.data(isk0 + 51);
    const auto *isk0_64 = buffer.data(isk0 + 64);
    const auto *isk0_72 = buffer.data(isk0 + 72);
    const auto *isk0_77 = buffer.data(isk0 + 77);
    const auto *isk0_81 = buffer.data(isk0 + 81);
    const auto *isk0_84 = buffer.data(isk0 + 84);
    const auto *isk0_86 = buffer.data(isk0 + 86);
    const auto *isk0_89 = buffer.data(isk0 + 89);
    const auto *isk0_90 = buffer.data(isk0 + 90);
    const auto *isk0_92 = buffer.data(isk0 + 92);
    const auto *isk0_107 = buffer.data(isk0 + 107);

    const auto *isi_28 = buffer.data(isi + 28);
    const auto *isi_31 = buffer.data(isi + 31);
    const auto *isi_34 = buffer.data(isi + 34);
    const auto *isi_38 = buffer.data(isi + 38);
    const auto *isi_49 = buffer.data(isi + 49);
    const auto *isi_55 = buffer.data(isi + 55);
    const auto *isi_56 = buffer.data(isi + 56);
    const auto *isi_58 = buffer.data(isi + 58);
    const auto *isi_61 = buffer.data(isi + 61);
    const auto *isi_64 = buffer.data(isi + 64);
    const auto *isi_65 = buffer.data(isi + 65);
    const auto *isi_68 = buffer.data(isi + 68);
    const auto *isi_69 = buffer.data(isi + 69);
    const auto *isi_70 = buffer.data(isi + 70);
    const auto *isi_79 = buffer.data(isi + 79);
    const auto *isi_80 = buffer.data(isi + 80);
    const auto *isi_81 = buffer.data(isi + 81);
    const auto *isi_82 = buffer.data(isi + 82);
    const auto *isi_83 = buffer.data(isi + 83);
    const auto *isi_84 = buffer.data(isi + 84);
    const auto *isi_89 = buffer.data(isi + 89);
    const auto *isi_93 = buffer.data(isi + 93);
    const auto *isi_98 = buffer.data(isi + 98);
    const auto *isi_105 = buffer.data(isi + 105);
    const auto *isi_107 = buffer.data(isi + 107);
    const auto *isi_108 = buffer.data(isi + 108);
    const auto *isi_109 = buffer.data(isi + 109);
    const auto *isi_110 = buffer.data(isi + 110);
    const auto *isi_111 = buffer.data(isi + 111);
    const auto *isi_133 = buffer.data(isi + 133);
    const auto *isi_134 = buffer.data(isi + 134);
    const auto *isi_135 = buffer.data(isi + 135);
    const auto *isi_136 = buffer.data(isi + 136);
    const auto *isi_137 = buffer.data(isi + 137);
    const auto *isi_138 = buffer.data(isi + 138);
    const auto *isi_139 = buffer.data(isi + 139);
    const auto *isi_140 = buffer.data(isi + 140);
    const auto *isi_145 = buffer.data(isi + 145);
    const auto *isi_149 = buffer.data(isi + 149);
    const auto *isi_154 = buffer.data(isi + 154);
    const auto *isi_160 = buffer.data(isi + 160);
    const auto *isi_161 = buffer.data(isi + 161);
    const auto *isi_162 = buffer.data(isi + 162);
    const auto *isi_163 = buffer.data(isi + 163);
    const auto *isi_164 = buffer.data(isi + 164);
    const auto *isi_165 = buffer.data(isi + 165);
    const auto *isi_167 = buffer.data(isi + 167);
    const auto *isi_168 = buffer.data(isi + 168);
    const auto *isi_171 = buffer.data(isi + 171);
    const auto *isi_174 = buffer.data(isi + 174);
    const auto *isi_178 = buffer.data(isi + 178);
    const auto *isi_183 = buffer.data(isi + 183);
    const auto *isi_189 = buffer.data(isi + 189);
    const auto *isi_191 = buffer.data(isi + 191);
    const auto *isi_192 = buffer.data(isi + 192);
    const auto *isi_193 = buffer.data(isi + 193);
    const auto *isi_194 = buffer.data(isi + 194);
    const auto *isi_195 = buffer.data(isi + 195);

    const auto *isk1_39 = buffer.data(isk1 + 39);
    const auto *isk1_42 = buffer.data(isk1 + 42);
    const auto *isk1_46 = buffer.data(isk1 + 46);
    const auto *isk1_51 = buffer.data(isk1 + 51);
    const auto *isk1_64 = buffer.data(isk1 + 64);
    const auto *isk1_72 = buffer.data(isk1 + 72);
    const auto *isk1_77 = buffer.data(isk1 + 77);
    const auto *isk1_81 = buffer.data(isk1 + 81);
    const auto *isk1_84 = buffer.data(isk1 + 84);
    const auto *isk1_86 = buffer.data(isk1 + 86);
    const auto *isk1_89 = buffer.data(isk1 + 89);
    const auto *isk1_90 = buffer.data(isk1 + 90);
    const auto *isk1_92 = buffer.data(isk1 + 92);
    const auto *isk1_107 = buffer.data(isk1 + 107);

    const auto *ksh0_78 = buffer.data(ksh0 + 78);
    const auto *ksh0_79 = buffer.data(ksh0 + 79);
    const auto *ksh0_80 = buffer.data(ksh0 + 80);
    const auto *ksh0_81 = buffer.data(ksh0 + 81);
    const auto *ksh0_83 = buffer.data(ksh0 + 83);
    const auto *ksh0_101 = buffer.data(ksh0 + 101);
    const auto *ksh0_102 = buffer.data(ksh0 + 102);
    const auto *ksh0_103 = buffer.data(ksh0 + 103);
    const auto *ksh0_104 = buffer.data(ksh0 + 104);
    const auto *ksh0_105 = buffer.data(ksh0 + 105);
    const auto *ksh0_106 = buffer.data(ksh0 + 106);
    const auto *ksh0_107 = buffer.data(ksh0 + 107);
    const auto *ksh0_108 = buffer.data(ksh0 + 108);
    const auto *ksh0_109 = buffer.data(ksh0 + 109);
    const auto *ksh0_110 = buffer.data(ksh0 + 110);
    const auto *ksh0_111 = buffer.data(ksh0 + 111);
    const auto *ksh0_112 = buffer.data(ksh0 + 112);
    const auto *ksh0_113 = buffer.data(ksh0 + 113);
    const auto *ksh0_114 = buffer.data(ksh0 + 114);
    const auto *ksh0_119 = buffer.data(ksh0 + 119);
    const auto *ksh0_120 = buffer.data(ksh0 + 120);
    const auto *ksh0_121 = buffer.data(ksh0 + 121);
    const auto *ksh0_122 = buffer.data(ksh0 + 122);
    const auto *ksh0_123 = buffer.data(ksh0 + 123);
    const auto *ksh0_124 = buffer.data(ksh0 + 124);
    const auto *ksh0_125 = buffer.data(ksh0 + 125);
    const auto *ksh0_126 = buffer.data(ksh0 + 126);
    const auto *ksh0_128 = buffer.data(ksh0 + 128);
    const auto *ksh0_129 = buffer.data(ksh0 + 129);
    const auto *ksh0_131 = buffer.data(ksh0 + 131);
    const auto *ksh0_132 = buffer.data(ksh0 + 132);
    const auto *ksh0_133 = buffer.data(ksh0 + 133);
    const auto *ksh0_135 = buffer.data(ksh0 + 135);
    const auto *ksh0_136 = buffer.data(ksh0 + 136);
    const auto *ksh0_141 = buffer.data(ksh0 + 141);
    const auto *ksh0_142 = buffer.data(ksh0 + 142);
    const auto *ksh0_143 = buffer.data(ksh0 + 143);
    const auto *ksh0_144 = buffer.data(ksh0 + 144);

    const auto *ksh1_78 = buffer.data(ksh1 + 78);
    const auto *ksh1_79 = buffer.data(ksh1 + 79);
    const auto *ksh1_80 = buffer.data(ksh1 + 80);
    const auto *ksh1_81 = buffer.data(ksh1 + 81);
    const auto *ksh1_83 = buffer.data(ksh1 + 83);
    const auto *ksh1_101 = buffer.data(ksh1 + 101);
    const auto *ksh1_102 = buffer.data(ksh1 + 102);
    const auto *ksh1_103 = buffer.data(ksh1 + 103);
    const auto *ksh1_104 = buffer.data(ksh1 + 104);
    const auto *ksh1_105 = buffer.data(ksh1 + 105);
    const auto *ksh1_106 = buffer.data(ksh1 + 106);
    const auto *ksh1_107 = buffer.data(ksh1 + 107);
    const auto *ksh1_108 = buffer.data(ksh1 + 108);
    const auto *ksh1_109 = buffer.data(ksh1 + 109);
    const auto *ksh1_110 = buffer.data(ksh1 + 110);
    const auto *ksh1_111 = buffer.data(ksh1 + 111);
    const auto *ksh1_112 = buffer.data(ksh1 + 112);
    const auto *ksh1_113 = buffer.data(ksh1 + 113);
    const auto *ksh1_114 = buffer.data(ksh1 + 114);
    const auto *ksh1_119 = buffer.data(ksh1 + 119);
    const auto *ksh1_120 = buffer.data(ksh1 + 120);
    const auto *ksh1_121 = buffer.data(ksh1 + 121);
    const auto *ksh1_122 = buffer.data(ksh1 + 122);
    const auto *ksh1_123 = buffer.data(ksh1 + 123);
    const auto *ksh1_124 = buffer.data(ksh1 + 124);
    const auto *ksh1_125 = buffer.data(ksh1 + 125);
    const auto *ksh1_126 = buffer.data(ksh1 + 126);
    const auto *ksh1_128 = buffer.data(ksh1 + 128);
    const auto *ksh1_129 = buffer.data(ksh1 + 129);
    const auto *ksh1_131 = buffer.data(ksh1 + 131);
    const auto *ksh1_132 = buffer.data(ksh1 + 132);
    const auto *ksh1_133 = buffer.data(ksh1 + 133);
    const auto *ksh1_135 = buffer.data(ksh1 + 135);
    const auto *ksh1_136 = buffer.data(ksh1 + 136);
    const auto *ksh1_141 = buffer.data(ksh1 + 141);
    const auto *ksh1_142 = buffer.data(ksh1 + 142);
    const auto *ksh1_143 = buffer.data(ksh1 + 143);
    const auto *ksh1_144 = buffer.data(ksh1 + 144);

    const auto *ksi_99 = buffer.data(ksi + 99);
    const auto *ksi_105 = buffer.data(ksi + 105);
    const auto *ksi_106 = buffer.data(ksi + 106);
    const auto *ksi_107 = buffer.data(ksi + 107);
    const auto *ksi_108 = buffer.data(ksi + 108);
    const auto *ksi_109 = buffer.data(ksi + 109);
    const auto *ksi_110 = buffer.data(ksi + 110);
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
    const auto *ksi_134 = buffer.data(ksi + 134);
    const auto *ksi_135 = buffer.data(ksi + 135);
    const auto *ksi_136 = buffer.data(ksi + 136);
    const auto *ksi_137 = buffer.data(ksi + 137);
    const auto *ksi_138 = buffer.data(ksi + 138);
    const auto *ksi_139 = buffer.data(ksi + 139);
    const auto *ksi_140 = buffer.data(ksi + 140);
    const auto *ksi_141 = buffer.data(ksi + 141);
    const auto *ksi_142 = buffer.data(ksi + 142);
    const auto *ksi_143 = buffer.data(ksi + 143);
    const auto *ksi_144 = buffer.data(ksi + 144);
    const auto *ksi_145 = buffer.data(ksi + 145);
    const auto *ksi_146 = buffer.data(ksi + 146);
    const auto *ksi_147 = buffer.data(ksi + 147);
    const auto *ksi_148 = buffer.data(ksi + 148);
    const auto *ksi_149 = buffer.data(ksi + 149);
    const auto *ksi_150 = buffer.data(ksi + 150);
    const auto *ksi_151 = buffer.data(ksi + 151);
    const auto *ksi_152 = buffer.data(ksi + 152);
    const auto *ksi_153 = buffer.data(ksi + 153);
    const auto *ksi_154 = buffer.data(ksi + 154);
    const auto *ksi_160 = buffer.data(ksi + 160);
    const auto *ksi_161 = buffer.data(ksi + 161);
    const auto *ksi_162 = buffer.data(ksi + 162);
    const auto *ksi_163 = buffer.data(ksi + 163);
    const auto *ksi_164 = buffer.data(ksi + 164);
    const auto *ksi_165 = buffer.data(ksi + 165);
    const auto *ksi_166 = buffer.data(ksi + 166);
    const auto *ksi_167 = buffer.data(ksi + 167);
    const auto *ksi_168 = buffer.data(ksi + 168);
    const auto *ksi_169 = buffer.data(ksi + 169);
    const auto *ksi_170 = buffer.data(ksi + 170);
    const auto *ksi_171 = buffer.data(ksi + 171);
    const auto *ksi_173 = buffer.data(ksi + 173);
    const auto *ksi_174 = buffer.data(ksi + 174);
    const auto *ksi_175 = buffer.data(ksi + 175);
    const auto *ksi_177 = buffer.data(ksi + 177);
    const auto *ksi_178 = buffer.data(ksi + 178);
    const auto *ksi_179 = buffer.data(ksi + 179);
    const auto *ksi_180 = buffer.data(ksi + 180);
    const auto *ksi_182 = buffer.data(ksi + 182);
    const auto *ksi_183 = buffer.data(ksi + 183);
    const auto *ksi_189 = buffer.data(ksi + 189);
    const auto *ksi_190 = buffer.data(ksi + 190);
    const auto *ksi_191 = buffer.data(ksi + 191);
    const auto *ksi_192 = buffer.data(ksi + 192);
    const auto *ksi_193 = buffer.data(ksi + 193);
    const auto *ksi_194 = buffer.data(ksi + 194);
    const auto *ksi_195 = buffer.data(ksi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, isi_105, isi_107, \
                         isi_108, isi_109, ksi_99, ksi_105, ksi_107, ksi_108, \
                         ksi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_17 * isi_105[k]
                   + f_3 * pc_x[k] * ksi_105[k];

        t_130[k] = f_3 * pc_z[k] * ksi_99[k];

        t_131[k] = f_17 * isi_107[k]
                   + f_3 * pc_x[k] * ksi_107[k];

        t_132[k] = f_17 * isi_108[k]
                   + f_3 * pc_x[k] * ksi_108[k];

        t_133[k] = f_17 * isi_109[k]
                   + f_3 * pc_x[k] * ksi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, isi_49, isi_110, \
                         isi_111, ksh0_78, ksh1_78, ksi_105, ksi_110, \
                         ksi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_17 * isi_110[k]
                   + f_3 * pc_x[k] * ksi_110[k];

        t_135[k] = f_17 * isi_111[k]
                   + f_3 * pc_x[k] * ksi_111[k];

        t_136[k] = f_14 * isi_49[k]
                   + f_1 * ksh0_78[k]
                   - f_2 * ksh1_78[k]
                   + f_3 * pc_y[k] * ksi_105[k];

        t_137[k] = f_3 * pc_z[k] * ksi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, ksh0_78, ksh0_79, ksh0_80, ksh1_78, \
                         ksh1_79, ksh1_80, ksi_106, ksi_107, ksi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * ksh0_78[k]
                   - f_5 * ksh1_78[k]
                   + f_3 * pc_z[k] * ksi_106[k];

        t_139[k] = f_6 * ksh0_79[k]
                   - f_7 * ksh1_79[k]
                   + f_3 * pc_z[k] * ksi_107[k];

        t_140[k] = f_8 * ksh0_80[k]
                   - f_9 * ksh1_80[k]
                   + f_3 * pc_z[k] * ksi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, isk0_72, isi_55, \
                         isk1_72, ksh0_81, ksh0_83, ksh1_81, ksh1_83, ksi_109, \
                         ksi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * ksh0_81[k]
                   - f_11 * ksh1_81[k]
                   + f_3 * pc_z[k] * ksi_109[k];

        t_142[k] = f_14 * isi_55[k]
                   + f_3 * pc_y[k] * ksi_111[k];

        t_143[k] = f_1 * ksh0_83[k]
                   - f_2 * ksh1_83[k]
                   + f_3 * pc_z[k] * ksi_111[k];

        t_144[k] = pa_y[k] * isk0_72[k]
                   - f_12 * pc_y[k] * isk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, isk0_39, isi_28, \
                         isi_56, isi_58, isk1_39, ksi_112, ksi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * isi_56[k]
                   + f_3 * pc_y[k] * ksi_112[k];

        t_146[k] = f_13 * isi_28[k]
                   + f_3 * pc_z[k] * ksi_112[k];

        t_147[k] = pa_z[k] * isk0_39[k]
                   - f_12 * pc_z[k] * isk1_39[k];

        t_148[k] = f_13 * isi_58[k]
                   + f_3 * pc_y[k] * ksi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, isk0_42, isk0_77, \
                         isi_31, isi_61, isk1_42, isk1_77, ksi_115, \
                         ksi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * isk0_77[k]
                   - f_12 * pc_y[k] * isk1_77[k];

        t_150[k] = pa_z[k] * isk0_42[k]
                   - f_12 * pc_z[k] * isk1_42[k];

        t_151[k] = f_13 * isi_31[k]
                   + f_3 * pc_z[k] * ksi_115[k];

        t_152[k] = f_13 * isi_61[k]
                   + f_3 * pc_y[k] * ksi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, isk0_46, isk0_81, \
                         isi_34, isk1_46, isk1_81, ksi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * isk0_81[k]
                   - f_12 * pc_y[k] * isk1_81[k];

        t_154[k] = pa_z[k] * isk0_46[k]
                   - f_12 * pc_z[k] * isk1_46[k];

        t_155[k] = f_13 * isi_34[k]
                   + f_3 * pc_z[k] * ksi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, isk0_84, isk0_86, isi_64, isi_65, \
                         isk1_84, isk1_86, ksi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * isk0_84[k]
                   + f_14 * isi_64[k]
                   - f_12 * pc_y[k] * isk1_84[k];

        t_157[k] = f_13 * isi_65[k]
                   + f_3 * pc_y[k] * ksi_121[k];

        t_158[k] = pa_y[k] * isk0_86[k]
                   - f_12 * pc_y[k] * isk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, isk0_51, isk0_89, \
                         isi_38, isi_68, isk1_51, isk1_89, ksi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * isk0_51[k]
                   - f_12 * pc_z[k] * isk1_51[k];

        t_160[k] = f_13 * isi_38[k]
                   + f_3 * pc_z[k] * ksi_122[k];

        t_161[k] = pa_y[k] * isk0_89[k]
                   + f_15 * isi_68[k]
                   - f_12 * pc_y[k] * isk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, isk0_90, isk0_92, \
                         isi_69, isi_70, isi_133, isk1_90, isk1_92, ksi_126, \
                         ksi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * isk0_90[k]
                   + f_14 * isi_69[k]
                   - f_12 * pc_y[k] * isk1_90[k];

        t_163[k] = f_13 * isi_70[k]
                   + f_3 * pc_y[k] * ksi_126[k];

        t_164[k] = pa_y[k] * isk0_92[k]
                   - f_12 * pc_y[k] * isk1_92[k];

        t_165[k] = f_17 * isi_133[k]
                   + f_3 * pc_x[k] * ksi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, isi_134, isi_135, isi_136, \
                         isi_137, isi_138, ksi_134, ksi_135, ksi_136, ksi_137, \
                         ksi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_17 * isi_134[k]
                   + f_3 * pc_x[k] * ksi_134[k];

        t_167[k] = f_17 * isi_135[k]
                   + f_3 * pc_x[k] * ksi_135[k];

        t_168[k] = f_17 * isi_136[k]
                   + f_3 * pc_x[k] * ksi_136[k];

        t_169[k] = f_17 * isi_137[k]
                   + f_3 * pc_x[k] * ksi_137[k];

        t_170[k] = f_17 * isi_138[k]
                   + f_3 * pc_x[k] * ksi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, isk0_64, isi_49, isi_139, \
                         isk1_64, ksi_133, ksi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_17 * isi_139[k]
                   + f_3 * pc_x[k] * ksi_139[k];

        t_172[k] = pa_z[k] * isk0_64[k]
                   - f_12 * pc_z[k] * isk1_64[k];

        t_173[k] = f_13 * isi_49[k]
                   + f_3 * pc_z[k] * ksi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, isi_79, isi_80, isi_81, ksh0_101, \
                         ksh0_102, ksh0_103, ksh1_101, ksh1_102, ksh1_103, ksi_135, ksi_136, \
                         ksi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * isi_79[k]
                   + f_10 * ksh0_101[k]
                   - f_11 * ksh1_101[k]
                   + f_3 * pc_y[k] * ksi_135[k];

        t_175[k] = f_13 * isi_80[k]
                   + f_8 * ksh0_102[k]
                   - f_9 * ksh1_102[k]
                   + f_3 * pc_y[k] * ksi_136[k];

        t_176[k] = f_13 * isi_81[k]
                   + f_6 * ksh0_103[k]
                   - f_7 * ksh1_103[k]
                   + f_3 * pc_y[k] * ksi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, isk0_107, isi_82, isi_83, isk1_107, \
                         ksh0_104, ksh1_104, ksi_138, ksi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * isi_82[k]
                   + f_4 * ksh0_104[k]
                   - f_5 * ksh1_104[k]
                   + f_3 * pc_y[k] * ksi_138[k];

        t_178[k] = f_13 * isi_83[k]
                   + f_3 * pc_y[k] * ksi_139[k];

        t_179[k] = pa_y[k] * isk0_107[k]
                   - f_12 * pc_y[k] * isk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, isi_56, isi_140, \
                         ksh0_105, ksh1_105, ksi_140, ksi_141, \
                         ksi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_17 * isi_140[k]
                   + f_1 * ksh0_105[k]
                   - f_2 * ksh1_105[k]
                   + f_3 * pc_x[k] * ksi_140[k];

        t_181[k] = f_3 * pc_y[k] * ksi_140[k];

        t_182[k] = f_14 * isi_56[k]
                   + f_3 * pc_z[k] * ksi_140[k];

        t_183[k] = f_4 * ksh0_105[k]
                   - f_5 * ksh1_105[k]
                   + f_3 * pc_y[k] * ksi_141[k];

        t_184[k] = f_3 * pc_y[k] * ksi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, isi_145, ksh0_106, ksh0_107, \
                         ksh0_110, ksh1_106, ksh1_107, ksh1_110, ksi_143, ksi_144, \
                         ksi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_17 * isi_145[k]
                   + f_10 * ksh0_110[k]
                   - f_11 * ksh1_110[k]
                   + f_3 * pc_x[k] * ksi_145[k];

        t_186[k] = f_6 * ksh0_106[k]
                   - f_7 * ksh1_106[k]
                   + f_3 * pc_y[k] * ksi_143[k];

        t_187[k] = f_4 * ksh0_107[k]
                   - f_5 * ksh1_107[k]
                   + f_3 * pc_y[k] * ksi_144[k];

        t_188[k] = f_3 * pc_y[k] * ksi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, isi_149, ksh0_108, ksh0_109, \
                         ksh0_114, ksh1_108, ksh1_109, ksh1_114, ksi_146, ksi_147, \
                         ksi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_17 * isi_149[k]
                   + f_8 * ksh0_114[k]
                   - f_9 * ksh1_114[k]
                   + f_3 * pc_x[k] * ksi_149[k];

        t_190[k] = f_8 * ksh0_108[k]
                   - f_9 * ksh1_108[k]
                   + f_3 * pc_y[k] * ksi_146[k];

        t_191[k] = f_6 * ksh0_109[k]
                   - f_7 * ksh1_109[k]
                   + f_3 * pc_y[k] * ksi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, isi_154, ksh0_110, ksh0_119, \
                         ksh1_110, ksh1_119, ksi_148, ksi_149, \
                         ksi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * ksh0_110[k]
                   - f_5 * ksh1_110[k]
                   + f_3 * pc_y[k] * ksi_148[k];

        t_193[k] = f_3 * pc_y[k] * ksi_149[k];

        t_194[k] = f_17 * isi_154[k]
                   + f_6 * ksh0_119[k]
                   - f_7 * ksh1_119[k]
                   + f_3 * pc_x[k] * ksi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, ksh0_111, ksh0_112, ksh0_113, ksh1_111, \
                         ksh1_112, ksh1_113, ksi_150, ksi_151, \
                         ksi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * ksh0_111[k]
                   - f_11 * ksh1_111[k]
                   + f_3 * pc_y[k] * ksi_150[k];

        t_196[k] = f_8 * ksh0_112[k]
                   - f_9 * ksh1_112[k]
                   + f_3 * pc_y[k] * ksi_151[k];

        t_197[k] = f_6 * ksh0_113[k]
                   - f_7 * ksh1_113[k]
                   + f_3 * pc_y[k] * ksi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, isi_160, isi_161, ksh0_114, \
                         ksh0_125, ksh1_114, ksh1_125, ksi_153, ksi_154, ksi_160, \
                         ksi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * ksh0_114[k]
                   - f_5 * ksh1_114[k]
                   + f_3 * pc_y[k] * ksi_153[k];

        t_199[k] = f_3 * pc_y[k] * ksi_154[k];

        t_200[k] = f_17 * isi_160[k]
                   + f_4 * ksh0_125[k]
                   - f_5 * ksh1_125[k]
                   + f_3 * pc_x[k] * ksi_160[k];

        t_201[k] = f_17 * isi_161[k]
                   + f_3 * pc_x[k] * ksi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, isi_162, isi_163, \
                         isi_164, isi_165, ksi_160, ksi_162, ksi_163, ksi_164, \
                         ksi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_17 * isi_162[k]
                   + f_3 * pc_x[k] * ksi_162[k];

        t_203[k] = f_17 * isi_163[k]
                   + f_3 * pc_x[k] * ksi_163[k];

        t_204[k] = f_17 * isi_164[k]
                   + f_3 * pc_x[k] * ksi_164[k];

        t_205[k] = f_17 * isi_165[k]
                   + f_3 * pc_x[k] * ksi_165[k];

        t_206[k] = f_3 * pc_y[k] * ksi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, isi_167, ksh0_120, ksh0_121, \
                         ksh1_120, ksh1_121, ksi_161, ksi_162, \
                         ksi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_17 * isi_167[k]
                   + f_3 * pc_x[k] * ksi_167[k];

        t_208[k] = f_1 * ksh0_120[k]
                   - f_2 * ksh1_120[k]
                   + f_3 * pc_y[k] * ksi_161[k];

        t_209[k] = f_19 * ksh0_121[k]
                   - f_20 * ksh1_121[k]
                   + f_3 * pc_y[k] * ksi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, ksh0_122, ksh0_123, ksh0_124, ksh1_122, \
                         ksh1_123, ksh1_124, ksi_163, ksi_164, \
                         ksi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * ksh0_122[k]
                   - f_11 * ksh1_122[k]
                   + f_3 * pc_y[k] * ksi_163[k];

        t_211[k] = f_8 * ksh0_123[k]
                   - f_9 * ksh1_123[k]
                   + f_3 * pc_y[k] * ksi_164[k];

        t_212[k] = f_6 * ksh0_124[k]
                   - f_7 * ksh1_124[k]
                   + f_3 * pc_y[k] * ksi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, isi_83, isi_168, \
                         ksh0_125, ksh0_126, ksh1_125, ksh1_126, ksi_166, ksi_167, \
                         ksi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * ksh0_125[k]
                   - f_5 * ksh1_125[k]
                   + f_3 * pc_y[k] * ksi_166[k];

        t_214[k] = f_3 * pc_y[k] * ksi_167[k];

        t_215[k] = f_14 * isi_83[k]
                   + f_1 * ksh0_125[k]
                   - f_2 * ksh1_125[k]
                   + f_3 * pc_z[k] * ksi_167[k];

        t_216[k] = f_16 * isi_168[k]
                   + f_1 * ksh0_126[k]
                   - f_2 * ksh1_126[k]
                   + f_3 * pc_x[k] * ksi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, isi_84, isi_171, \
                         ksh0_129, ksh1_129, ksi_168, ksi_169, \
                         ksi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * isi_84[k]
                   + f_3 * pc_y[k] * ksi_168[k];

        t_218[k] = f_3 * pc_z[k] * ksi_168[k];

        t_219[k] = f_16 * isi_171[k]
                   + f_10 * ksh0_129[k]
                   - f_11 * ksh1_129[k]
                   + f_3 * pc_x[k] * ksi_171[k];

        t_220[k] = f_3 * pc_z[k] * ksi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pc_x, pc_z, isi_174, ksh0_126, ksh0_132, \
                         ksh1_126, ksh1_132, ksi_170, ksi_171, \
                         ksi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * ksh0_126[k]
                   - f_5 * ksh1_126[k]
                   + f_3 * pc_z[k] * ksi_170[k];

        t_222[k] = f_16 * isi_174[k]
                   + f_8 * ksh0_132[k]
                   - f_9 * ksh1_132[k]
                   + f_3 * pc_x[k] * ksi_174[k];

        t_223[k] = f_3 * pc_z[k] * ksi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, isi_89, isi_178, \
                         ksh0_128, ksh0_136, ksh1_128, ksh1_136, ksi_173, ksi_174, \
                         ksi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * isi_89[k]
                   + f_3 * pc_y[k] * ksi_173[k];

        t_225[k] = f_6 * ksh0_128[k]
                   - f_7 * ksh1_128[k]
                   + f_3 * pc_z[k] * ksi_173[k];

        t_226[k] = f_16 * isi_178[k]
                   + f_6 * ksh0_136[k]
                   - f_7 * ksh1_136[k]
                   + f_3 * pc_x[k] * ksi_178[k];

        t_227[k] = f_3 * pc_z[k] * ksi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, isi_93, ksh0_129, ksh0_131, \
                         ksh1_129, ksh1_131, ksi_175, ksi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * ksh0_129[k]
                   - f_5 * ksh1_129[k]
                   + f_3 * pc_z[k] * ksi_175[k];

        t_229[k] = f_15 * isi_93[k]
                   + f_3 * pc_y[k] * ksi_177[k];

        t_230[k] = f_8 * ksh0_131[k]
                   - f_9 * ksh1_131[k]
                   + f_3 * pc_z[k] * ksi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_x, pc_z, isi_183, ksh0_132, ksh0_141, \
                         ksh1_132, ksh1_141, ksi_178, ksi_179, \
                         ksi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_16 * isi_183[k]
                   + f_4 * ksh0_141[k]
                   - f_5 * ksh1_141[k]
                   + f_3 * pc_x[k] * ksi_183[k];

        t_232[k] = f_3 * pc_z[k] * ksi_178[k];

        t_233[k] = f_4 * ksh0_132[k]
                   - f_5 * ksh1_132[k]
                   + f_3 * pc_z[k] * ksi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, isi_98, isi_189, \
                         ksh0_133, ksh0_135, ksh1_133, ksh1_135, ksi_180, ksi_182, \
                         ksi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * ksh0_133[k]
                   - f_7 * ksh1_133[k]
                   + f_3 * pc_z[k] * ksi_180[k];

        t_235[k] = f_15 * isi_98[k]
                   + f_3 * pc_y[k] * ksi_182[k];

        t_236[k] = f_10 * ksh0_135[k]
                   - f_11 * ksh1_135[k]
                   + f_3 * pc_z[k] * ksi_182[k];

        t_237[k] = f_16 * isi_189[k]
                   + f_3 * pc_x[k] * ksi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, isi_191, isi_192, \
                         isi_193, isi_194, ksi_183, ksi_191, ksi_192, ksi_193, \
                         ksi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * ksi_183[k];

        t_239[k] = f_16 * isi_191[k]
                   + f_3 * pc_x[k] * ksi_191[k];

        t_240[k] = f_16 * isi_192[k]
                   + f_3 * pc_x[k] * ksi_192[k];

        t_241[k] = f_16 * isi_193[k]
                   + f_3 * pc_x[k] * ksi_193[k];

        t_242[k] = f_16 * isi_194[k]
                   + f_3 * pc_x[k] * ksi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, isi_105, isi_195, \
                         ksh0_141, ksh1_141, ksi_189, ksi_190, \
                         ksi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_16 * isi_195[k]
                   + f_3 * pc_x[k] * ksi_195[k];

        t_244[k] = f_15 * isi_105[k]
                   + f_1 * ksh0_141[k]
                   - f_2 * ksh1_141[k]
                   + f_3 * pc_y[k] * ksi_189[k];

        t_245[k] = f_3 * pc_z[k] * ksi_189[k];

        t_246[k] = f_4 * ksh0_141[k]
                   - f_5 * ksh1_141[k]
                   + f_3 * pc_z[k] * ksi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, ksh0_142, ksh0_143, ksh0_144, ksh1_142, \
                         ksh1_143, ksh1_144, ksi_191, ksi_192, \
                         ksi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * ksh0_142[k]
                   - f_7 * ksh1_142[k]
                   + f_3 * pc_z[k] * ksi_191[k];

        t_248[k] = f_8 * ksh0_143[k]
                   - f_9 * ksh1_143[k]
                   + f_3 * pc_z[k] * ksi_192[k];

        t_249[k] = f_10 * ksh0_144[k]
                   - f_11 * ksh1_144[k]
                   + f_3 * pc_z[k] * ksi_193[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
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

    const auto *isk0_108 = buffer.data(isk0 + 108);
    const auto *isk0_111 = buffer.data(isk0 + 111);
    const auto *isk0_114 = buffer.data(isk0 + 114);
    const auto *isk0_118 = buffer.data(isk0 + 118);
    const auto *isk0_120 = buffer.data(isk0 + 120);
    const auto *isk0_123 = buffer.data(isk0 + 123);
    const auto *isk0_125 = buffer.data(isk0 + 125);
    const auto *isk0_126 = buffer.data(isk0 + 126);
    const auto *isk0_136 = buffer.data(isk0 + 136);
    const auto *isk0_180 = buffer.data(isk0 + 180);
    const auto *isk0_183 = buffer.data(isk0 + 183);
    const auto *isk0_185 = buffer.data(isk0 + 185);
    const auto *isk0_186 = buffer.data(isk0 + 186);
    const auto *isk0_189 = buffer.data(isk0 + 189);
    const auto *isk0_190 = buffer.data(isk0 + 190);
    const auto *isk0_192 = buffer.data(isk0 + 192);
    const auto *isk0_194 = buffer.data(isk0 + 194);
    const auto *isk0_195 = buffer.data(isk0 + 195);
    const auto *isk0_197 = buffer.data(isk0 + 197);
    const auto *isk0_198 = buffer.data(isk0 + 198);
    const auto *isk0_200 = buffer.data(isk0 + 200);
    const auto *isk0_215 = buffer.data(isk0 + 215);

    const auto *isi_84 = buffer.data(isi + 84);
    const auto *isi_87 = buffer.data(isi + 87);
    const auto *isi_90 = buffer.data(isi + 90);
    const auto *isi_91 = buffer.data(isi + 91);
    const auto *isi_94 = buffer.data(isi + 94);
    const auto *isi_95 = buffer.data(isi + 95);
    const auto *isi_96 = buffer.data(isi + 96);
    const auto *isi_105 = buffer.data(isi + 105);
    const auto *isi_111 = buffer.data(isi + 111);
    const auto *isi_112 = buffer.data(isi + 112);
    const auto *isi_114 = buffer.data(isi + 114);
    const auto *isi_115 = buffer.data(isi + 115);
    const auto *isi_117 = buffer.data(isi + 117);
    const auto *isi_118 = buffer.data(isi + 118);
    const auto *isi_121 = buffer.data(isi + 121);
    const auto *isi_122 = buffer.data(isi + 122);
    const auto *isi_126 = buffer.data(isi + 126);
    const auto *isi_133 = buffer.data(isi + 133);
    const auto *isi_135 = buffer.data(isi + 135);
    const auto *isi_136 = buffer.data(isi + 136);
    const auto *isi_137 = buffer.data(isi + 137);
    const auto *isi_138 = buffer.data(isi + 138);
    const auto *isi_139 = buffer.data(isi + 139);
    const auto *isi_140 = buffer.data(isi + 140);
    const auto *isi_141 = buffer.data(isi + 141);
    const auto *isi_142 = buffer.data(isi + 142);
    const auto *isi_143 = buffer.data(isi + 143);
    const auto *isi_145 = buffer.data(isi + 145);
    const auto *isi_146 = buffer.data(isi + 146);
    const auto *isi_148 = buffer.data(isi + 148);
    const auto *isi_149 = buffer.data(isi + 149);
    const auto *isi_150 = buffer.data(isi + 150);
    const auto *isi_152 = buffer.data(isi + 152);
    const auto *isi_153 = buffer.data(isi + 153);
    const auto *isi_154 = buffer.data(isi + 154);
    const auto *isi_161 = buffer.data(isi + 161);
    const auto *isi_163 = buffer.data(isi + 163);
    const auto *isi_164 = buffer.data(isi + 164);
    const auto *isi_165 = buffer.data(isi + 165);
    const auto *isi_166 = buffer.data(isi + 166);
    const auto *isi_167 = buffer.data(isi + 167);
    const auto *isi_168 = buffer.data(isi + 168);
    const auto *isi_201 = buffer.data(isi + 201);
    const auto *isi_205 = buffer.data(isi + 205);
    const auto *isi_210 = buffer.data(isi + 210);
    const auto *isi_216 = buffer.data(isi + 216);
    const auto *isi_217 = buffer.data(isi + 217);
    const auto *isi_218 = buffer.data(isi + 218);
    const auto *isi_219 = buffer.data(isi + 219);
    const auto *isi_220 = buffer.data(isi + 220);
    const auto *isi_221 = buffer.data(isi + 221);
    const auto *isi_222 = buffer.data(isi + 222);
    const auto *isi_223 = buffer.data(isi + 223);
    const auto *isi_245 = buffer.data(isi + 245);
    const auto *isi_246 = buffer.data(isi + 246);
    const auto *isi_247 = buffer.data(isi + 247);
    const auto *isi_248 = buffer.data(isi + 248);
    const auto *isi_249 = buffer.data(isi + 249);
    const auto *isi_250 = buffer.data(isi + 250);
    const auto *isi_251 = buffer.data(isi + 251);
    const auto *isi_252 = buffer.data(isi + 252);
    const auto *isi_257 = buffer.data(isi + 257);
    const auto *isi_261 = buffer.data(isi + 261);
    const auto *isi_266 = buffer.data(isi + 266);
    const auto *isi_272 = buffer.data(isi + 272);
    const auto *isi_273 = buffer.data(isi + 273);
    const auto *isi_274 = buffer.data(isi + 274);
    const auto *isi_275 = buffer.data(isi + 275);
    const auto *isi_276 = buffer.data(isi + 276);
    const auto *isi_277 = buffer.data(isi + 277);
    const auto *isi_279 = buffer.data(isi + 279);
    const auto *isi_280 = buffer.data(isi + 280);
    const auto *isi_283 = buffer.data(isi + 283);
    const auto *isi_286 = buffer.data(isi + 286);

    const auto *isk1_108 = buffer.data(isk1 + 108);
    const auto *isk1_111 = buffer.data(isk1 + 111);
    const auto *isk1_114 = buffer.data(isk1 + 114);
    const auto *isk1_118 = buffer.data(isk1 + 118);
    const auto *isk1_120 = buffer.data(isk1 + 120);
    const auto *isk1_123 = buffer.data(isk1 + 123);
    const auto *isk1_125 = buffer.data(isk1 + 125);
    const auto *isk1_126 = buffer.data(isk1 + 126);
    const auto *isk1_136 = buffer.data(isk1 + 136);
    const auto *isk1_180 = buffer.data(isk1 + 180);
    const auto *isk1_183 = buffer.data(isk1 + 183);
    const auto *isk1_185 = buffer.data(isk1 + 185);
    const auto *isk1_186 = buffer.data(isk1 + 186);
    const auto *isk1_189 = buffer.data(isk1 + 189);
    const auto *isk1_190 = buffer.data(isk1 + 190);
    const auto *isk1_192 = buffer.data(isk1 + 192);
    const auto *isk1_194 = buffer.data(isk1 + 194);
    const auto *isk1_195 = buffer.data(isk1 + 195);
    const auto *isk1_197 = buffer.data(isk1 + 197);
    const auto *isk1_198 = buffer.data(isk1 + 198);
    const auto *isk1_200 = buffer.data(isk1 + 200);
    const auto *isk1_215 = buffer.data(isk1 + 215);

    const auto *ksh0_146 = buffer.data(ksh0 + 146);
    const auto *ksh0_152 = buffer.data(ksh0 + 152);
    const auto *ksh0_156 = buffer.data(ksh0 + 156);
    const auto *ksh0_161 = buffer.data(ksh0 + 161);
    const auto *ksh0_164 = buffer.data(ksh0 + 164);
    const auto *ksh0_165 = buffer.data(ksh0 + 165);
    const auto *ksh0_166 = buffer.data(ksh0 + 166);
    const auto *ksh0_167 = buffer.data(ksh0 + 167);
    const auto *ksh0_183 = buffer.data(ksh0 + 183);
    const auto *ksh0_185 = buffer.data(ksh0 + 185);
    const auto *ksh0_186 = buffer.data(ksh0 + 186);
    const auto *ksh0_187 = buffer.data(ksh0 + 187);
    const auto *ksh0_188 = buffer.data(ksh0 + 188);
    const auto *ksh0_189 = buffer.data(ksh0 + 189);
    const auto *ksh0_190 = buffer.data(ksh0 + 190);
    const auto *ksh0_191 = buffer.data(ksh0 + 191);
    const auto *ksh0_192 = buffer.data(ksh0 + 192);
    const auto *ksh0_193 = buffer.data(ksh0 + 193);
    const auto *ksh0_194 = buffer.data(ksh0 + 194);
    const auto *ksh0_195 = buffer.data(ksh0 + 195);
    const auto *ksh0_196 = buffer.data(ksh0 + 196);
    const auto *ksh0_197 = buffer.data(ksh0 + 197);
    const auto *ksh0_198 = buffer.data(ksh0 + 198);
    const auto *ksh0_203 = buffer.data(ksh0 + 203);
    const auto *ksh0_204 = buffer.data(ksh0 + 204);
    const auto *ksh0_205 = buffer.data(ksh0 + 205);
    const auto *ksh0_206 = buffer.data(ksh0 + 206);
    const auto *ksh0_207 = buffer.data(ksh0 + 207);
    const auto *ksh0_208 = buffer.data(ksh0 + 208);
    const auto *ksh0_209 = buffer.data(ksh0 + 209);
    const auto *ksh0_210 = buffer.data(ksh0 + 210);
    const auto *ksh0_213 = buffer.data(ksh0 + 213);
    const auto *ksh0_216 = buffer.data(ksh0 + 216);

    const auto *ksh1_146 = buffer.data(ksh1 + 146);
    const auto *ksh1_152 = buffer.data(ksh1 + 152);
    const auto *ksh1_156 = buffer.data(ksh1 + 156);
    const auto *ksh1_161 = buffer.data(ksh1 + 161);
    const auto *ksh1_164 = buffer.data(ksh1 + 164);
    const auto *ksh1_165 = buffer.data(ksh1 + 165);
    const auto *ksh1_166 = buffer.data(ksh1 + 166);
    const auto *ksh1_167 = buffer.data(ksh1 + 167);
    const auto *ksh1_183 = buffer.data(ksh1 + 183);
    const auto *ksh1_185 = buffer.data(ksh1 + 185);
    const auto *ksh1_186 = buffer.data(ksh1 + 186);
    const auto *ksh1_187 = buffer.data(ksh1 + 187);
    const auto *ksh1_188 = buffer.data(ksh1 + 188);
    const auto *ksh1_189 = buffer.data(ksh1 + 189);
    const auto *ksh1_190 = buffer.data(ksh1 + 190);
    const auto *ksh1_191 = buffer.data(ksh1 + 191);
    const auto *ksh1_192 = buffer.data(ksh1 + 192);
    const auto *ksh1_193 = buffer.data(ksh1 + 193);
    const auto *ksh1_194 = buffer.data(ksh1 + 194);
    const auto *ksh1_195 = buffer.data(ksh1 + 195);
    const auto *ksh1_196 = buffer.data(ksh1 + 196);
    const auto *ksh1_197 = buffer.data(ksh1 + 197);
    const auto *ksh1_198 = buffer.data(ksh1 + 198);
    const auto *ksh1_203 = buffer.data(ksh1 + 203);
    const auto *ksh1_204 = buffer.data(ksh1 + 204);
    const auto *ksh1_205 = buffer.data(ksh1 + 205);
    const auto *ksh1_206 = buffer.data(ksh1 + 206);
    const auto *ksh1_207 = buffer.data(ksh1 + 207);
    const auto *ksh1_208 = buffer.data(ksh1 + 208);
    const auto *ksh1_209 = buffer.data(ksh1 + 209);
    const auto *ksh1_210 = buffer.data(ksh1 + 210);
    const auto *ksh1_213 = buffer.data(ksh1 + 213);
    const auto *ksh1_216 = buffer.data(ksh1 + 216);

    const auto *ksi_195 = buffer.data(ksi + 195);
    const auto *ksi_196 = buffer.data(ksi + 196);
    const auto *ksi_198 = buffer.data(ksi + 198);
    const auto *ksi_199 = buffer.data(ksi + 199);
    const auto *ksi_201 = buffer.data(ksi + 201);
    const auto *ksi_202 = buffer.data(ksi + 202);
    const auto *ksi_205 = buffer.data(ksi + 205);
    const auto *ksi_206 = buffer.data(ksi + 206);
    const auto *ksi_210 = buffer.data(ksi + 210);
    const auto *ksi_216 = buffer.data(ksi + 216);
    const auto *ksi_217 = buffer.data(ksi + 217);
    const auto *ksi_218 = buffer.data(ksi + 218);
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
    const auto *ksi_234 = buffer.data(ksi + 234);
    const auto *ksi_238 = buffer.data(ksi + 238);
    const auto *ksi_245 = buffer.data(ksi + 245);
    const auto *ksi_246 = buffer.data(ksi + 246);
    const auto *ksi_247 = buffer.data(ksi + 247);
    const auto *ksi_248 = buffer.data(ksi + 248);
    const auto *ksi_249 = buffer.data(ksi + 249);
    const auto *ksi_250 = buffer.data(ksi + 250);
    const auto *ksi_251 = buffer.data(ksi + 251);
    const auto *ksi_252 = buffer.data(ksi + 252);
    const auto *ksi_253 = buffer.data(ksi + 253);
    const auto *ksi_254 = buffer.data(ksi + 254);
    const auto *ksi_255 = buffer.data(ksi + 255);
    const auto *ksi_256 = buffer.data(ksi + 256);
    const auto *ksi_257 = buffer.data(ksi + 257);
    const auto *ksi_258 = buffer.data(ksi + 258);
    const auto *ksi_259 = buffer.data(ksi + 259);
    const auto *ksi_260 = buffer.data(ksi + 260);
    const auto *ksi_261 = buffer.data(ksi + 261);
    const auto *ksi_262 = buffer.data(ksi + 262);
    const auto *ksi_263 = buffer.data(ksi + 263);
    const auto *ksi_264 = buffer.data(ksi + 264);
    const auto *ksi_265 = buffer.data(ksi + 265);
    const auto *ksi_266 = buffer.data(ksi + 266);
    const auto *ksi_272 = buffer.data(ksi + 272);
    const auto *ksi_273 = buffer.data(ksi + 273);
    const auto *ksi_274 = buffer.data(ksi + 274);
    const auto *ksi_275 = buffer.data(ksi + 275);
    const auto *ksi_276 = buffer.data(ksi + 276);
    const auto *ksi_277 = buffer.data(ksi + 277);
    const auto *ksi_278 = buffer.data(ksi + 278);
    const auto *ksi_279 = buffer.data(ksi + 279);
    const auto *ksi_280 = buffer.data(ksi + 280);
    const auto *ksi_281 = buffer.data(ksi + 281);
    const auto *ksi_282 = buffer.data(ksi + 282);
    const auto *ksi_283 = buffer.data(ksi + 283);
    const auto *ksi_286 = buffer.data(ksi + 286);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, isk0_108, isi_111, \
                         isi_112, isk1_108, ksh0_146, ksh1_146, ksi_195, \
                         ksi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * isi_111[k]
                   + f_3 * pc_y[k] * ksi_195[k];

        t_251[k] = f_1 * ksh0_146[k]
                   - f_2 * ksh1_146[k]
                   + f_3 * pc_z[k] * ksi_195[k];

        t_252[k] = pa_z[k] * isk0_108[k]
                   - f_12 * pc_z[k] * isk1_108[k];

        t_253[k] = f_14 * isi_112[k]
                   + f_3 * pc_y[k] * ksi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_y, pc_z, isk0_111, isi_84, isi_114, \
                         isk1_111, ksi_196, ksi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * isi_84[k]
                   + f_3 * pc_z[k] * ksi_196[k];

        t_255[k] = pa_z[k] * isk0_111[k]
                   - f_12 * pc_z[k] * isk1_111[k];

        t_256[k] = f_14 * isi_114[k]
                   + f_3 * pc_y[k] * ksi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, isk0_114, isi_87, isi_201, \
                         isk1_114, ksh0_152, ksh1_152, ksi_199, \
                         ksi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_16 * isi_201[k]
                   + f_10 * ksh0_152[k]
                   - f_11 * ksh1_152[k]
                   + f_3 * pc_x[k] * ksi_201[k];

        t_258[k] = pa_z[k] * isk0_114[k]
                   - f_12 * pc_z[k] * isk1_114[k];

        t_259[k] = f_13 * isi_87[k]
                   + f_3 * pc_z[k] * ksi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_y, pc_z, isk0_118, isi_117, \
                         isi_205, isk1_118, ksh0_156, ksh1_156, ksi_201, \
                         ksi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * isi_117[k]
                   + f_3 * pc_y[k] * ksi_201[k];

        t_261[k] = f_16 * isi_205[k]
                   + f_8 * ksh0_156[k]
                   - f_9 * ksh1_156[k]
                   + f_3 * pc_x[k] * ksi_205[k];

        t_262[k] = pa_z[k] * isk0_118[k]
                   - f_12 * pc_z[k] * isk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pa_z, pc_y, pc_z, isk0_120, isi_90, isi_91, \
                         isi_121, isk1_120, ksi_202, ksi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * isi_90[k]
                   + f_3 * pc_z[k] * ksi_202[k];

        t_264[k] = pa_z[k] * isk0_120[k]
                   + f_14 * isi_91[k]
                   - f_12 * pc_z[k] * isk1_120[k];

        t_265[k] = f_14 * isi_121[k]
                   + f_3 * pc_y[k] * ksi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, isk0_123, isi_94, isi_210, \
                         isk1_123, ksh0_161, ksh1_161, ksi_206, \
                         ksi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_16 * isi_210[k]
                   + f_6 * ksh0_161[k]
                   - f_7 * ksh1_161[k]
                   + f_3 * pc_x[k] * ksi_210[k];

        t_267[k] = pa_z[k] * isk0_123[k]
                   - f_12 * pc_z[k] * isk1_123[k];

        t_268[k] = f_13 * isi_94[k]
                   + f_3 * pc_z[k] * ksi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_z, pc_y, pc_z, isk0_125, isk0_126, isi_95, \
                         isi_96, isi_126, isk1_125, isk1_126, ksi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_z[k] * isk0_125[k]
                   + f_14 * isi_95[k]
                   - f_12 * pc_z[k] * isk1_125[k];

        t_270[k] = pa_z[k] * isk0_126[k]
                   + f_15 * isi_96[k]
                   - f_12 * pc_z[k] * isk1_126[k];

        t_271[k] = f_14 * isi_126[k]
                   + f_3 * pc_y[k] * ksi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, isi_216, isi_217, isi_218, isi_219, \
                         ksh0_167, ksh1_167, ksi_216, ksi_217, ksi_218, \
                         ksi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_16 * isi_216[k]
                   + f_4 * ksh0_167[k]
                   - f_5 * ksh1_167[k]
                   + f_3 * pc_x[k] * ksi_216[k];

        t_273[k] = f_16 * isi_217[k]
                   + f_3 * pc_x[k] * ksi_217[k];

        t_274[k] = f_16 * isi_218[k]
                   + f_3 * pc_x[k] * ksi_218[k];

        t_275[k] = f_16 * isi_219[k]
                   + f_3 * pc_x[k] * ksi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, isi_220, isi_221, isi_222, isi_223, \
                         ksi_220, ksi_221, ksi_222, ksi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_16 * isi_220[k]
                   + f_3 * pc_x[k] * ksi_220[k];

        t_277[k] = f_16 * isi_221[k]
                   + f_3 * pc_x[k] * ksi_221[k];

        t_278[k] = f_16 * isi_222[k]
                   + f_3 * pc_x[k] * ksi_222[k];

        t_279[k] = f_16 * isi_223[k]
                   + f_3 * pc_x[k] * ksi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_z, pc_y, pc_z, isk0_136, isi_105, isi_135, \
                         isk1_136, ksh0_164, ksh1_164, ksi_217, \
                         ksi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * isk0_136[k]
                   - f_12 * pc_z[k] * isk1_136[k];

        t_281[k] = f_13 * isi_105[k]
                   + f_3 * pc_z[k] * ksi_217[k];

        t_282[k] = f_14 * isi_135[k]
                   + f_10 * ksh0_164[k]
                   - f_11 * ksh1_164[k]
                   + f_3 * pc_y[k] * ksi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, isi_136, isi_137, isi_138, ksh0_165, \
                         ksh0_166, ksh0_167, ksh1_165, ksh1_166, ksh1_167, ksi_220, ksi_221, \
                         ksi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * isi_136[k]
                   + f_8 * ksh0_165[k]
                   - f_9 * ksh1_165[k]
                   + f_3 * pc_y[k] * ksi_220[k];

        t_284[k] = f_14 * isi_137[k]
                   + f_6 * ksh0_166[k]
                   - f_7 * ksh1_166[k]
                   + f_3 * pc_y[k] * ksi_221[k];

        t_285[k] = f_14 * isi_138[k]
                   + f_4 * ksh0_167[k]
                   - f_5 * ksh1_167[k]
                   + f_3 * pc_y[k] * ksi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, isk0_180, isi_111, \
                         isi_139, isi_140, isk1_180, ksh0_167, ksh1_167, ksi_223, \
                         ksi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * isi_139[k]
                   + f_3 * pc_y[k] * ksi_223[k];

        t_287[k] = f_13 * isi_111[k]
                   + f_1 * ksh0_167[k]
                   - f_2 * ksh1_167[k]
                   + f_3 * pc_z[k] * ksi_223[k];

        t_288[k] = pa_y[k] * isk0_180[k]
                   - f_12 * pc_y[k] * isk1_180[k];

        t_289[k] = f_13 * isi_140[k]
                   + f_3 * pc_y[k] * ksi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, pc_z, isk0_183, isk0_185, \
                         isi_112, isi_141, isi_142, isk1_183, isk1_185, ksi_224, \
                         ksi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * isi_112[k]
                   + f_3 * pc_z[k] * ksi_224[k];

        t_291[k] = pa_y[k] * isk0_183[k]
                   + f_14 * isi_141[k]
                   - f_12 * pc_y[k] * isk1_183[k];

        t_292[k] = f_13 * isi_142[k]
                   + f_3 * pc_y[k] * ksi_226[k];

        t_293[k] = pa_y[k] * isk0_185[k]
                   - f_12 * pc_y[k] * isk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pc_y, pc_z, isk0_186, isk0_189, \
                         isi_115, isi_143, isi_145, isk1_186, isk1_189, ksi_227, \
                         ksi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * isk0_186[k]
                   + f_15 * isi_143[k]
                   - f_12 * pc_y[k] * isk1_186[k];

        t_295[k] = f_14 * isi_115[k]
                   + f_3 * pc_z[k] * ksi_227[k];

        t_296[k] = f_13 * isi_145[k]
                   + f_3 * pc_y[k] * ksi_229[k];

        t_297[k] = pa_y[k] * isk0_189[k]
                   - f_12 * pc_y[k] * isk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_y, pc_y, pc_z, isk0_190, isk0_192, isi_118, \
                         isi_146, isi_148, isk1_190, isk1_192, \
                         ksi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_y[k] * isk0_190[k]
                   + f_16 * isi_146[k]
                   - f_12 * pc_y[k] * isk1_190[k];

        t_299[k] = f_14 * isi_118[k]
                   + f_3 * pc_z[k] * ksi_230[k];

        t_300[k] = pa_y[k] * isk0_192[k]
                   + f_14 * isi_148[k]
                   - f_12 * pc_y[k] * isk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, pc_y, pc_z, isk0_194, isk0_195, \
                         isi_122, isi_149, isi_150, isk1_194, isk1_195, ksi_233, \
                         ksi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * isi_149[k]
                   + f_3 * pc_y[k] * ksi_233[k];

        t_302[k] = pa_y[k] * isk0_194[k]
                   - f_12 * pc_y[k] * isk1_194[k];

        t_303[k] = pa_y[k] * isk0_195[k]
                   + f_17 * isi_150[k]
                   - f_12 * pc_y[k] * isk1_195[k];

        t_304[k] = f_14 * isi_122[k]
                   + f_3 * pc_z[k] * ksi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pc_y, isk0_197, isk0_198, isk0_200, \
                         isi_152, isi_153, isi_154, isk1_197, isk1_198, isk1_200, \
                         ksi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_y[k] * isk0_197[k]
                   + f_15 * isi_152[k]
                   - f_12 * pc_y[k] * isk1_197[k];

        t_306[k] = pa_y[k] * isk0_198[k]
                   + f_14 * isi_153[k]
                   - f_12 * pc_y[k] * isk1_198[k];

        t_307[k] = f_13 * isi_154[k]
                   + f_3 * pc_y[k] * ksi_238[k];

        t_308[k] = pa_y[k] * isk0_200[k]
                   - f_12 * pc_y[k] * isk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, isi_245, isi_246, isi_247, \
                         isi_248, isi_249, ksi_245, ksi_246, ksi_247, ksi_248, \
                         ksi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_16 * isi_245[k]
                   + f_3 * pc_x[k] * ksi_245[k];

        t_310[k] = f_16 * isi_246[k]
                   + f_3 * pc_x[k] * ksi_246[k];

        t_311[k] = f_16 * isi_247[k]
                   + f_3 * pc_x[k] * ksi_247[k];

        t_312[k] = f_16 * isi_248[k]
                   + f_3 * pc_x[k] * ksi_248[k];

        t_313[k] = f_16 * isi_249[k]
                   + f_3 * pc_x[k] * ksi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, isi_133, isi_161, \
                         isi_250, isi_251, ksh0_183, ksh1_183, ksi_245, ksi_250, \
                         ksi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_16 * isi_250[k]
                   + f_3 * pc_x[k] * ksi_250[k];

        t_315[k] = f_16 * isi_251[k]
                   + f_3 * pc_x[k] * ksi_251[k];

        t_316[k] = f_13 * isi_161[k]
                   + f_1 * ksh0_183[k]
                   - f_2 * ksh1_183[k]
                   + f_3 * pc_y[k] * ksi_245[k];

        t_317[k] = f_14 * isi_133[k]
                   + f_3 * pc_z[k] * ksi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, isi_163, isi_164, isi_165, ksh0_185, \
                         ksh0_186, ksh0_187, ksh1_185, ksh1_186, ksh1_187, ksi_247, ksi_248, \
                         ksi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * isi_163[k]
                   + f_10 * ksh0_185[k]
                   - f_11 * ksh1_185[k]
                   + f_3 * pc_y[k] * ksi_247[k];

        t_319[k] = f_13 * isi_164[k]
                   + f_8 * ksh0_186[k]
                   - f_9 * ksh1_186[k]
                   + f_3 * pc_y[k] * ksi_248[k];

        t_320[k] = f_13 * isi_165[k]
                   + f_6 * ksh0_187[k]
                   - f_7 * ksh1_187[k]
                   + f_3 * pc_y[k] * ksi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, isk0_215, isi_166, isi_167, \
                         isk1_215, ksh0_188, ksh1_188, ksi_250, \
                         ksi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * isi_166[k]
                   + f_4 * ksh0_188[k]
                   - f_5 * ksh1_188[k]
                   + f_3 * pc_y[k] * ksi_250[k];

        t_322[k] = f_13 * isi_167[k]
                   + f_3 * pc_y[k] * ksi_251[k];

        t_323[k] = pa_y[k] * isk0_215[k]
                   - f_12 * pc_y[k] * isk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, isi_140, \
                         isi_252, ksh0_189, ksh1_189, ksi_252, ksi_253, \
                         ksi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_16 * isi_252[k]
                   + f_1 * ksh0_189[k]
                   - f_2 * ksh1_189[k]
                   + f_3 * pc_x[k] * ksi_252[k];

        t_325[k] = f_3 * pc_y[k] * ksi_252[k];

        t_326[k] = f_15 * isi_140[k]
                   + f_3 * pc_z[k] * ksi_252[k];

        t_327[k] = f_4 * ksh0_189[k]
                   - f_5 * ksh1_189[k]
                   + f_3 * pc_y[k] * ksi_253[k];

        t_328[k] = f_3 * pc_y[k] * ksi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, isi_257, ksh0_190, ksh0_191, \
                         ksh0_194, ksh1_190, ksh1_191, ksh1_194, ksi_255, ksi_256, \
                         ksi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_16 * isi_257[k]
                   + f_10 * ksh0_194[k]
                   - f_11 * ksh1_194[k]
                   + f_3 * pc_x[k] * ksi_257[k];

        t_330[k] = f_6 * ksh0_190[k]
                   - f_7 * ksh1_190[k]
                   + f_3 * pc_y[k] * ksi_255[k];

        t_331[k] = f_4 * ksh0_191[k]
                   - f_5 * ksh1_191[k]
                   + f_3 * pc_y[k] * ksi_256[k];

        t_332[k] = f_3 * pc_y[k] * ksi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, isi_261, ksh0_192, ksh0_193, \
                         ksh0_198, ksh1_192, ksh1_193, ksh1_198, ksi_258, ksi_259, \
                         ksi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_16 * isi_261[k]
                   + f_8 * ksh0_198[k]
                   - f_9 * ksh1_198[k]
                   + f_3 * pc_x[k] * ksi_261[k];

        t_334[k] = f_8 * ksh0_192[k]
                   - f_9 * ksh1_192[k]
                   + f_3 * pc_y[k] * ksi_258[k];

        t_335[k] = f_6 * ksh0_193[k]
                   - f_7 * ksh1_193[k]
                   + f_3 * pc_y[k] * ksi_259[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, isi_266, ksh0_194, ksh0_203, \
                         ksh1_194, ksh1_203, ksi_260, ksi_261, \
                         ksi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_4 * ksh0_194[k]
                   - f_5 * ksh1_194[k]
                   + f_3 * pc_y[k] * ksi_260[k];

        t_337[k] = f_3 * pc_y[k] * ksi_261[k];

        t_338[k] = f_16 * isi_266[k]
                   + f_6 * ksh0_203[k]
                   - f_7 * ksh1_203[k]
                   + f_3 * pc_x[k] * ksi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_y, ksh0_195, ksh0_196, ksh0_197, ksh1_195, \
                         ksh1_196, ksh1_197, ksi_262, ksi_263, \
                         ksi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_10 * ksh0_195[k]
                   - f_11 * ksh1_195[k]
                   + f_3 * pc_y[k] * ksi_262[k];

        t_340[k] = f_8 * ksh0_196[k]
                   - f_9 * ksh1_196[k]
                   + f_3 * pc_y[k] * ksi_263[k];

        t_341[k] = f_6 * ksh0_197[k]
                   - f_7 * ksh1_197[k]
                   + f_3 * pc_y[k] * ksi_264[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_y, isi_272, isi_273, ksh0_198, \
                         ksh0_209, ksh1_198, ksh1_209, ksi_265, ksi_266, ksi_272, \
                         ksi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_4 * ksh0_198[k]
                   - f_5 * ksh1_198[k]
                   + f_3 * pc_y[k] * ksi_265[k];

        t_343[k] = f_3 * pc_y[k] * ksi_266[k];

        t_344[k] = f_16 * isi_272[k]
                   + f_4 * ksh0_209[k]
                   - f_5 * ksh1_209[k]
                   + f_3 * pc_x[k] * ksi_272[k];

        t_345[k] = f_16 * isi_273[k]
                   + f_3 * pc_x[k] * ksi_273[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pc_x, pc_y, isi_274, isi_275, \
                         isi_276, isi_277, ksi_272, ksi_274, ksi_275, ksi_276, \
                         ksi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_16 * isi_274[k]
                   + f_3 * pc_x[k] * ksi_274[k];

        t_347[k] = f_16 * isi_275[k]
                   + f_3 * pc_x[k] * ksi_275[k];

        t_348[k] = f_16 * isi_276[k]
                   + f_3 * pc_x[k] * ksi_276[k];

        t_349[k] = f_16 * isi_277[k]
                   + f_3 * pc_x[k] * ksi_277[k];

        t_350[k] = f_3 * pc_y[k] * ksi_272[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_x, pc_y, isi_279, ksh0_204, ksh0_205, \
                         ksh1_204, ksh1_205, ksi_273, ksi_274, \
                         ksi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_16 * isi_279[k]
                   + f_3 * pc_x[k] * ksi_279[k];

        t_352[k] = f_1 * ksh0_204[k]
                   - f_2 * ksh1_204[k]
                   + f_3 * pc_y[k] * ksi_273[k];

        t_353[k] = f_19 * ksh0_205[k]
                   - f_20 * ksh1_205[k]
                   + f_3 * pc_y[k] * ksi_274[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_y, ksh0_206, ksh0_207, ksh0_208, ksh1_206, \
                         ksh1_207, ksh1_208, ksi_275, ksi_276, \
                         ksi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_10 * ksh0_206[k]
                   - f_11 * ksh1_206[k]
                   + f_3 * pc_y[k] * ksi_275[k];

        t_355[k] = f_8 * ksh0_207[k]
                   - f_9 * ksh1_207[k]
                   + f_3 * pc_y[k] * ksi_276[k];

        t_356[k] = f_6 * ksh0_208[k]
                   - f_7 * ksh1_208[k]
                   + f_3 * pc_y[k] * ksi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pc_x, pc_y, pc_z, isi_167, isi_280, \
                         ksh0_209, ksh0_210, ksh1_209, ksh1_210, ksi_278, ksi_279, \
                         ksi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_4 * ksh0_209[k]
                   - f_5 * ksh1_209[k]
                   + f_3 * pc_y[k] * ksi_278[k];

        t_358[k] = f_3 * pc_y[k] * ksi_279[k];

        t_359[k] = f_15 * isi_167[k]
                   + f_1 * ksh0_209[k]
                   - f_2 * ksh1_209[k]
                   + f_3 * pc_z[k] * ksi_279[k];

        t_360[k] = f_15 * isi_280[k]
                   + f_1 * ksh0_210[k]
                   - f_2 * ksh1_210[k]
                   + f_3 * pc_x[k] * ksi_280[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, isi_168, isi_283, \
                         ksh0_213, ksh1_213, ksi_280, ksi_281, \
                         ksi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_16 * isi_168[k]
                   + f_3 * pc_y[k] * ksi_280[k];

        t_362[k] = f_3 * pc_z[k] * ksi_280[k];

        t_363[k] = f_15 * isi_283[k]
                   + f_10 * ksh0_213[k]
                   - f_11 * ksh1_213[k]
                   + f_3 * pc_x[k] * ksi_283[k];

        t_364[k] = f_3 * pc_z[k] * ksi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pc_x, pc_z, isi_286, ksh0_210, ksh0_216, \
                         ksh1_210, ksh1_216, ksi_282, ksi_283, \
                         ksi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * ksh0_210[k]
                   - f_5 * ksh1_210[k]
                   + f_3 * pc_z[k] * ksi_282[k];

        t_366[k] = f_15 * isi_286[k]
                   + f_8 * ksh0_216[k]
                   - f_9 * ksh1_216[k]
                   + f_3 * pc_x[k] * ksi_286[k];

        t_367[k] = f_3 * pc_z[k] * ksi_283[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
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

    const auto *isk0_216 = buffer.data(isk0 + 216);
    const auto *isk0_219 = buffer.data(isk0 + 219);
    const auto *isk0_222 = buffer.data(isk0 + 222);
    const auto *isk0_226 = buffer.data(isk0 + 226);
    const auto *isk0_228 = buffer.data(isk0 + 228);
    const auto *isk0_231 = buffer.data(isk0 + 231);
    const auto *isk0_233 = buffer.data(isk0 + 233);
    const auto *isk0_234 = buffer.data(isk0 + 234);
    const auto *isk0_244 = buffer.data(isk0 + 244);
    const auto *isk0_324 = buffer.data(isk0 + 324);
    const auto *isk0_327 = buffer.data(isk0 + 327);
    const auto *isk0_329 = buffer.data(isk0 + 329);
    const auto *isk0_330 = buffer.data(isk0 + 330);
    const auto *isk0_333 = buffer.data(isk0 + 333);
    const auto *isk0_334 = buffer.data(isk0 + 334);
    const auto *isk0_336 = buffer.data(isk0 + 336);

    const auto *isi_168 = buffer.data(isi + 168);
    const auto *isi_171 = buffer.data(isi + 171);
    const auto *isi_173 = buffer.data(isi + 173);
    const auto *isi_174 = buffer.data(isi + 174);
    const auto *isi_175 = buffer.data(isi + 175);
    const auto *isi_177 = buffer.data(isi + 177);
    const auto *isi_178 = buffer.data(isi + 178);
    const auto *isi_179 = buffer.data(isi + 179);
    const auto *isi_180 = buffer.data(isi + 180);
    const auto *isi_182 = buffer.data(isi + 182);
    const auto *isi_189 = buffer.data(isi + 189);
    const auto *isi_195 = buffer.data(isi + 195);
    const auto *isi_196 = buffer.data(isi + 196);
    const auto *isi_198 = buffer.data(isi + 198);
    const auto *isi_199 = buffer.data(isi + 199);
    const auto *isi_201 = buffer.data(isi + 201);
    const auto *isi_202 = buffer.data(isi + 202);
    const auto *isi_205 = buffer.data(isi + 205);
    const auto *isi_206 = buffer.data(isi + 206);
    const auto *isi_210 = buffer.data(isi + 210);
    const auto *isi_217 = buffer.data(isi + 217);
    const auto *isi_219 = buffer.data(isi + 219);
    const auto *isi_220 = buffer.data(isi + 220);
    const auto *isi_221 = buffer.data(isi + 221);
    const auto *isi_222 = buffer.data(isi + 222);
    const auto *isi_223 = buffer.data(isi + 223);
    const auto *isi_224 = buffer.data(isi + 224);
    const auto *isi_226 = buffer.data(isi + 226);
    const auto *isi_227 = buffer.data(isi + 227);
    const auto *isi_229 = buffer.data(isi + 229);
    const auto *isi_230 = buffer.data(isi + 230);
    const auto *isi_233 = buffer.data(isi + 233);
    const auto *isi_238 = buffer.data(isi + 238);
    const auto *isi_245 = buffer.data(isi + 245);
    const auto *isi_247 = buffer.data(isi + 247);
    const auto *isi_248 = buffer.data(isi + 248);
    const auto *isi_249 = buffer.data(isi + 249);
    const auto *isi_250 = buffer.data(isi + 250);
    const auto *isi_251 = buffer.data(isi + 251);
    const auto *isi_252 = buffer.data(isi + 252);
    const auto *isi_253 = buffer.data(isi + 253);
    const auto *isi_254 = buffer.data(isi + 254);
    const auto *isi_255 = buffer.data(isi + 255);
    const auto *isi_257 = buffer.data(isi + 257);
    const auto *isi_258 = buffer.data(isi + 258);
    const auto *isi_260 = buffer.data(isi + 260);
    const auto *isi_290 = buffer.data(isi + 290);
    const auto *isi_295 = buffer.data(isi + 295);
    const auto *isi_301 = buffer.data(isi + 301);
    const auto *isi_303 = buffer.data(isi + 303);
    const auto *isi_304 = buffer.data(isi + 304);
    const auto *isi_305 = buffer.data(isi + 305);
    const auto *isi_306 = buffer.data(isi + 306);
    const auto *isi_307 = buffer.data(isi + 307);
    const auto *isi_313 = buffer.data(isi + 313);
    const auto *isi_317 = buffer.data(isi + 317);
    const auto *isi_322 = buffer.data(isi + 322);
    const auto *isi_328 = buffer.data(isi + 328);
    const auto *isi_329 = buffer.data(isi + 329);
    const auto *isi_330 = buffer.data(isi + 330);
    const auto *isi_331 = buffer.data(isi + 331);
    const auto *isi_332 = buffer.data(isi + 332);
    const auto *isi_333 = buffer.data(isi + 333);
    const auto *isi_334 = buffer.data(isi + 334);
    const auto *isi_335 = buffer.data(isi + 335);
    const auto *isi_336 = buffer.data(isi + 336);
    const auto *isi_339 = buffer.data(isi + 339);
    const auto *isi_341 = buffer.data(isi + 341);
    const auto *isi_342 = buffer.data(isi + 342);
    const auto *isi_345 = buffer.data(isi + 345);
    const auto *isi_346 = buffer.data(isi + 346);
    const auto *isi_348 = buffer.data(isi + 348);
    const auto *isi_350 = buffer.data(isi + 350);
    const auto *isi_351 = buffer.data(isi + 351);
    const auto *isi_353 = buffer.data(isi + 353);
    const auto *isi_354 = buffer.data(isi + 354);
    const auto *isi_356 = buffer.data(isi + 356);
    const auto *isi_357 = buffer.data(isi + 357);
    const auto *isi_358 = buffer.data(isi + 358);
    const auto *isi_359 = buffer.data(isi + 359);
    const auto *isi_360 = buffer.data(isi + 360);
    const auto *isi_361 = buffer.data(isi + 361);
    const auto *isi_362 = buffer.data(isi + 362);
    const auto *isi_363 = buffer.data(isi + 363);

    const auto *isk1_216 = buffer.data(isk1 + 216);
    const auto *isk1_219 = buffer.data(isk1 + 219);
    const auto *isk1_222 = buffer.data(isk1 + 222);
    const auto *isk1_226 = buffer.data(isk1 + 226);
    const auto *isk1_228 = buffer.data(isk1 + 228);
    const auto *isk1_231 = buffer.data(isk1 + 231);
    const auto *isk1_233 = buffer.data(isk1 + 233);
    const auto *isk1_234 = buffer.data(isk1 + 234);
    const auto *isk1_244 = buffer.data(isk1 + 244);
    const auto *isk1_324 = buffer.data(isk1 + 324);
    const auto *isk1_327 = buffer.data(isk1 + 327);
    const auto *isk1_329 = buffer.data(isk1 + 329);
    const auto *isk1_330 = buffer.data(isk1 + 330);
    const auto *isk1_333 = buffer.data(isk1 + 333);
    const auto *isk1_334 = buffer.data(isk1 + 334);
    const auto *isk1_336 = buffer.data(isk1 + 336);

    const auto *ksh0_212 = buffer.data(ksh0 + 212);
    const auto *ksh0_213 = buffer.data(ksh0 + 213);
    const auto *ksh0_215 = buffer.data(ksh0 + 215);
    const auto *ksh0_216 = buffer.data(ksh0 + 216);
    const auto *ksh0_217 = buffer.data(ksh0 + 217);
    const auto *ksh0_219 = buffer.data(ksh0 + 219);
    const auto *ksh0_220 = buffer.data(ksh0 + 220);
    const auto *ksh0_225 = buffer.data(ksh0 + 225);
    const auto *ksh0_226 = buffer.data(ksh0 + 226);
    const auto *ksh0_227 = buffer.data(ksh0 + 227);
    const auto *ksh0_228 = buffer.data(ksh0 + 228);
    const auto *ksh0_230 = buffer.data(ksh0 + 230);
    const auto *ksh0_236 = buffer.data(ksh0 + 236);
    const auto *ksh0_240 = buffer.data(ksh0 + 240);
    const auto *ksh0_245 = buffer.data(ksh0 + 245);
    const auto *ksh0_248 = buffer.data(ksh0 + 248);
    const auto *ksh0_249 = buffer.data(ksh0 + 249);
    const auto *ksh0_250 = buffer.data(ksh0 + 250);
    const auto *ksh0_251 = buffer.data(ksh0 + 251);
    const auto *ksh0_252 = buffer.data(ksh0 + 252);
    const auto *ksh0_255 = buffer.data(ksh0 + 255);
    const auto *ksh0_257 = buffer.data(ksh0 + 257);
    const auto *ksh0_258 = buffer.data(ksh0 + 258);
    const auto *ksh0_261 = buffer.data(ksh0 + 261);
    const auto *ksh0_262 = buffer.data(ksh0 + 262);
    const auto *ksh0_264 = buffer.data(ksh0 + 264);
    const auto *ksh0_266 = buffer.data(ksh0 + 266);
    const auto *ksh0_267 = buffer.data(ksh0 + 267);
    const auto *ksh0_269 = buffer.data(ksh0 + 269);
    const auto *ksh0_270 = buffer.data(ksh0 + 270);
    const auto *ksh0_271 = buffer.data(ksh0 + 271);
    const auto *ksh0_272 = buffer.data(ksh0 + 272);

    const auto *ksh1_212 = buffer.data(ksh1 + 212);
    const auto *ksh1_213 = buffer.data(ksh1 + 213);
    const auto *ksh1_215 = buffer.data(ksh1 + 215);
    const auto *ksh1_216 = buffer.data(ksh1 + 216);
    const auto *ksh1_217 = buffer.data(ksh1 + 217);
    const auto *ksh1_219 = buffer.data(ksh1 + 219);
    const auto *ksh1_220 = buffer.data(ksh1 + 220);
    const auto *ksh1_225 = buffer.data(ksh1 + 225);
    const auto *ksh1_226 = buffer.data(ksh1 + 226);
    const auto *ksh1_227 = buffer.data(ksh1 + 227);
    const auto *ksh1_228 = buffer.data(ksh1 + 228);
    const auto *ksh1_230 = buffer.data(ksh1 + 230);
    const auto *ksh1_236 = buffer.data(ksh1 + 236);
    const auto *ksh1_240 = buffer.data(ksh1 + 240);
    const auto *ksh1_245 = buffer.data(ksh1 + 245);
    const auto *ksh1_248 = buffer.data(ksh1 + 248);
    const auto *ksh1_249 = buffer.data(ksh1 + 249);
    const auto *ksh1_250 = buffer.data(ksh1 + 250);
    const auto *ksh1_251 = buffer.data(ksh1 + 251);
    const auto *ksh1_252 = buffer.data(ksh1 + 252);
    const auto *ksh1_255 = buffer.data(ksh1 + 255);
    const auto *ksh1_257 = buffer.data(ksh1 + 257);
    const auto *ksh1_258 = buffer.data(ksh1 + 258);
    const auto *ksh1_261 = buffer.data(ksh1 + 261);
    const auto *ksh1_262 = buffer.data(ksh1 + 262);
    const auto *ksh1_264 = buffer.data(ksh1 + 264);
    const auto *ksh1_266 = buffer.data(ksh1 + 266);
    const auto *ksh1_267 = buffer.data(ksh1 + 267);
    const auto *ksh1_269 = buffer.data(ksh1 + 269);
    const auto *ksh1_270 = buffer.data(ksh1 + 270);
    const auto *ksh1_271 = buffer.data(ksh1 + 271);
    const auto *ksh1_272 = buffer.data(ksh1 + 272);

    const auto *ksi_285 = buffer.data(ksi + 285);
    const auto *ksi_286 = buffer.data(ksi + 286);
    const auto *ksi_287 = buffer.data(ksi + 287);
    const auto *ksi_289 = buffer.data(ksi + 289);
    const auto *ksi_290 = buffer.data(ksi + 290);
    const auto *ksi_291 = buffer.data(ksi + 291);
    const auto *ksi_292 = buffer.data(ksi + 292);
    const auto *ksi_294 = buffer.data(ksi + 294);
    const auto *ksi_295 = buffer.data(ksi + 295);
    const auto *ksi_301 = buffer.data(ksi + 301);
    const auto *ksi_302 = buffer.data(ksi + 302);
    const auto *ksi_303 = buffer.data(ksi + 303);
    const auto *ksi_304 = buffer.data(ksi + 304);
    const auto *ksi_305 = buffer.data(ksi + 305);
    const auto *ksi_306 = buffer.data(ksi + 306);
    const auto *ksi_307 = buffer.data(ksi + 307);
    const auto *ksi_308 = buffer.data(ksi + 308);
    const auto *ksi_310 = buffer.data(ksi + 310);
    const auto *ksi_311 = buffer.data(ksi + 311);
    const auto *ksi_313 = buffer.data(ksi + 313);
    const auto *ksi_314 = buffer.data(ksi + 314);
    const auto *ksi_317 = buffer.data(ksi + 317);
    const auto *ksi_318 = buffer.data(ksi + 318);
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
    const auto *ksi_338 = buffer.data(ksi + 338);
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
    const auto *ksi_364 = buffer.data(ksi + 364);
    const auto *ksi_366 = buffer.data(ksi + 366);
    const auto *ksi_367 = buffer.data(ksi + 367);
    const auto *ksi_369 = buffer.data(ksi + 369);
    const auto *ksi_370 = buffer.data(ksi + 370);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pc_x, pc_y, pc_z, isi_173, isi_290, \
                         ksh0_212, ksh0_220, ksh1_212, ksh1_220, ksi_285, ksi_286, \
                         ksi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_16 * isi_173[k]
                   + f_3 * pc_y[k] * ksi_285[k];

        t_369[k] = f_6 * ksh0_212[k]
                   - f_7 * ksh1_212[k]
                   + f_3 * pc_z[k] * ksi_285[k];

        t_370[k] = f_15 * isi_290[k]
                   + f_6 * ksh0_220[k]
                   - f_7 * ksh1_220[k]
                   + f_3 * pc_x[k] * ksi_290[k];

        t_371[k] = f_3 * pc_z[k] * ksi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, isi_177, ksh0_213, ksh0_215, \
                         ksh1_213, ksh1_215, ksi_287, ksi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * ksh0_213[k]
                   - f_5 * ksh1_213[k]
                   + f_3 * pc_z[k] * ksi_287[k];

        t_373[k] = f_16 * isi_177[k]
                   + f_3 * pc_y[k] * ksi_289[k];

        t_374[k] = f_8 * ksh0_215[k]
                   - f_9 * ksh1_215[k]
                   + f_3 * pc_z[k] * ksi_289[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_z, isi_295, ksh0_216, ksh0_225, \
                         ksh1_216, ksh1_225, ksi_290, ksi_291, \
                         ksi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_15 * isi_295[k]
                   + f_4 * ksh0_225[k]
                   - f_5 * ksh1_225[k]
                   + f_3 * pc_x[k] * ksi_295[k];

        t_376[k] = f_3 * pc_z[k] * ksi_290[k];

        t_377[k] = f_4 * ksh0_216[k]
                   - f_5 * ksh1_216[k]
                   + f_3 * pc_z[k] * ksi_291[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, isi_182, isi_301, \
                         ksh0_217, ksh0_219, ksh1_217, ksh1_219, ksi_292, ksi_294, \
                         ksi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * ksh0_217[k]
                   - f_7 * ksh1_217[k]
                   + f_3 * pc_z[k] * ksi_292[k];

        t_379[k] = f_16 * isi_182[k]
                   + f_3 * pc_y[k] * ksi_294[k];

        t_380[k] = f_10 * ksh0_219[k]
                   - f_11 * ksh1_219[k]
                   + f_3 * pc_z[k] * ksi_294[k];

        t_381[k] = f_15 * isi_301[k]
                   + f_3 * pc_x[k] * ksi_301[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_z, isi_303, isi_304, \
                         isi_305, isi_306, ksi_295, ksi_303, ksi_304, ksi_305, \
                         ksi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_z[k] * ksi_295[k];

        t_383[k] = f_15 * isi_303[k]
                   + f_3 * pc_x[k] * ksi_303[k];

        t_384[k] = f_15 * isi_304[k]
                   + f_3 * pc_x[k] * ksi_304[k];

        t_385[k] = f_15 * isi_305[k]
                   + f_3 * pc_x[k] * ksi_305[k];

        t_386[k] = f_15 * isi_306[k]
                   + f_3 * pc_x[k] * ksi_306[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pc_y, pc_z, isi_189, isi_307, \
                         ksh0_225, ksh1_225, ksi_301, ksi_302, \
                         ksi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_15 * isi_307[k]
                   + f_3 * pc_x[k] * ksi_307[k];

        t_388[k] = f_16 * isi_189[k]
                   + f_1 * ksh0_225[k]
                   - f_2 * ksh1_225[k]
                   + f_3 * pc_y[k] * ksi_301[k];

        t_389[k] = f_3 * pc_z[k] * ksi_301[k];

        t_390[k] = f_4 * ksh0_225[k]
                   - f_5 * ksh1_225[k]
                   + f_3 * pc_z[k] * ksi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_z, ksh0_226, ksh0_227, ksh0_228, ksh1_226, \
                         ksh1_227, ksh1_228, ksi_303, ksi_304, \
                         ksi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_6 * ksh0_226[k]
                   - f_7 * ksh1_226[k]
                   + f_3 * pc_z[k] * ksi_303[k];

        t_392[k] = f_8 * ksh0_227[k]
                   - f_9 * ksh1_227[k]
                   + f_3 * pc_z[k] * ksi_304[k];

        t_393[k] = f_10 * ksh0_228[k]
                   - f_11 * ksh1_228[k]
                   + f_3 * pc_z[k] * ksi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_z, pc_y, pc_z, isk0_216, isi_195, \
                         isi_196, isk1_216, ksh0_230, ksh1_230, ksi_307, \
                         ksi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * isi_195[k]
                   + f_3 * pc_y[k] * ksi_307[k];

        t_395[k] = f_1 * ksh0_230[k]
                   - f_2 * ksh1_230[k]
                   + f_3 * pc_z[k] * ksi_307[k];

        t_396[k] = pa_z[k] * isk0_216[k]
                   - f_12 * pc_z[k] * isk1_216[k];

        t_397[k] = f_15 * isi_196[k]
                   + f_3 * pc_y[k] * ksi_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_z, pc_y, pc_z, isk0_219, isi_168, isi_198, \
                         isk1_219, ksi_308, ksi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * isi_168[k]
                   + f_3 * pc_z[k] * ksi_308[k];

        t_399[k] = pa_z[k] * isk0_219[k]
                   - f_12 * pc_z[k] * isk1_219[k];

        t_400[k] = f_15 * isi_198[k]
                   + f_3 * pc_y[k] * ksi_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_z, pc_x, pc_z, isk0_222, isi_171, isi_313, \
                         isk1_222, ksh0_236, ksh1_236, ksi_311, \
                         ksi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * isi_313[k]
                   + f_10 * ksh0_236[k]
                   - f_11 * ksh1_236[k]
                   + f_3 * pc_x[k] * ksi_313[k];

        t_402[k] = pa_z[k] * isk0_222[k]
                   - f_12 * pc_z[k] * isk1_222[k];

        t_403[k] = f_13 * isi_171[k]
                   + f_3 * pc_z[k] * ksi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_z, pc_x, pc_y, pc_z, isk0_226, isi_201, \
                         isi_317, isk1_226, ksh0_240, ksh1_240, ksi_313, \
                         ksi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * isi_201[k]
                   + f_3 * pc_y[k] * ksi_313[k];

        t_405[k] = f_15 * isi_317[k]
                   + f_8 * ksh0_240[k]
                   - f_9 * ksh1_240[k]
                   + f_3 * pc_x[k] * ksi_317[k];

        t_406[k] = pa_z[k] * isk0_226[k]
                   - f_12 * pc_z[k] * isk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_z, pc_y, pc_z, isk0_228, isi_174, isi_175, \
                         isi_205, isk1_228, ksi_314, ksi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * isi_174[k]
                   + f_3 * pc_z[k] * ksi_314[k];

        t_408[k] = pa_z[k] * isk0_228[k]
                   + f_14 * isi_175[k]
                   - f_12 * pc_z[k] * isk1_228[k];

        t_409[k] = f_15 * isi_205[k]
                   + f_3 * pc_y[k] * ksi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pa_z, pc_x, pc_z, isk0_231, isi_178, isi_322, \
                         isk1_231, ksh0_245, ksh1_245, ksi_318, \
                         ksi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_15 * isi_322[k]
                   + f_6 * ksh0_245[k]
                   - f_7 * ksh1_245[k]
                   + f_3 * pc_x[k] * ksi_322[k];

        t_411[k] = pa_z[k] * isk0_231[k]
                   - f_12 * pc_z[k] * isk1_231[k];

        t_412[k] = f_13 * isi_178[k]
                   + f_3 * pc_z[k] * ksi_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pa_z, pc_y, pc_z, isk0_233, isk0_234, isi_179, \
                         isi_180, isi_210, isk1_233, isk1_234, \
                         ksi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_z[k] * isk0_233[k]
                   + f_14 * isi_179[k]
                   - f_12 * pc_z[k] * isk1_233[k];

        t_414[k] = pa_z[k] * isk0_234[k]
                   + f_15 * isi_180[k]
                   - f_12 * pc_z[k] * isk1_234[k];

        t_415[k] = f_15 * isi_210[k]
                   + f_3 * pc_y[k] * ksi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, isi_328, isi_329, isi_330, isi_331, \
                         ksh0_251, ksh1_251, ksi_328, ksi_329, ksi_330, \
                         ksi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_15 * isi_328[k]
                   + f_4 * ksh0_251[k]
                   - f_5 * ksh1_251[k]
                   + f_3 * pc_x[k] * ksi_328[k];

        t_417[k] = f_15 * isi_329[k]
                   + f_3 * pc_x[k] * ksi_329[k];

        t_418[k] = f_15 * isi_330[k]
                   + f_3 * pc_x[k] * ksi_330[k];

        t_419[k] = f_15 * isi_331[k]
                   + f_3 * pc_x[k] * ksi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, isi_332, isi_333, isi_334, isi_335, \
                         ksi_332, ksi_333, ksi_334, ksi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_15 * isi_332[k]
                   + f_3 * pc_x[k] * ksi_332[k];

        t_421[k] = f_15 * isi_333[k]
                   + f_3 * pc_x[k] * ksi_333[k];

        t_422[k] = f_15 * isi_334[k]
                   + f_3 * pc_x[k] * ksi_334[k];

        t_423[k] = f_15 * isi_335[k]
                   + f_3 * pc_x[k] * ksi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_z, pc_y, pc_z, isk0_244, isi_189, isi_219, \
                         isk1_244, ksh0_248, ksh1_248, ksi_329, \
                         ksi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * isk0_244[k]
                   - f_12 * pc_z[k] * isk1_244[k];

        t_425[k] = f_13 * isi_189[k]
                   + f_3 * pc_z[k] * ksi_329[k];

        t_426[k] = f_15 * isi_219[k]
                   + f_10 * ksh0_248[k]
                   - f_11 * ksh1_248[k]
                   + f_3 * pc_y[k] * ksi_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, isi_220, isi_221, isi_222, ksh0_249, \
                         ksh0_250, ksh0_251, ksh1_249, ksh1_250, ksh1_251, ksi_332, ksi_333, \
                         ksi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * isi_220[k]
                   + f_8 * ksh0_249[k]
                   - f_9 * ksh1_249[k]
                   + f_3 * pc_y[k] * ksi_332[k];

        t_428[k] = f_15 * isi_221[k]
                   + f_6 * ksh0_250[k]
                   - f_7 * ksh1_250[k]
                   + f_3 * pc_y[k] * ksi_333[k];

        t_429[k] = f_15 * isi_222[k]
                   + f_4 * ksh0_251[k]
                   - f_5 * ksh1_251[k]
                   + f_3 * pc_y[k] * ksi_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, isi_195, isi_223, isi_336, \
                         ksh0_251, ksh0_252, ksh1_251, ksh1_252, ksi_335, \
                         ksi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * isi_223[k]
                   + f_3 * pc_y[k] * ksi_335[k];

        t_431[k] = f_13 * isi_195[k]
                   + f_1 * ksh0_251[k]
                   - f_2 * ksh1_251[k]
                   + f_3 * pc_z[k] * ksi_335[k];

        t_432[k] = f_15 * isi_336[k]
                   + f_1 * ksh0_252[k]
                   - f_2 * ksh1_252[k]
                   + f_3 * pc_x[k] * ksi_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, isi_196, isi_224, \
                         isi_226, isi_339, ksh0_255, ksh1_255, ksi_336, ksi_338, \
                         ksi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * isi_224[k]
                   + f_3 * pc_y[k] * ksi_336[k];

        t_434[k] = f_14 * isi_196[k]
                   + f_3 * pc_z[k] * ksi_336[k];

        t_435[k] = f_15 * isi_339[k]
                   + f_10 * ksh0_255[k]
                   - f_11 * ksh1_255[k]
                   + f_3 * pc_x[k] * ksi_339[k];

        t_436[k] = f_14 * isi_226[k]
                   + f_3 * pc_y[k] * ksi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, isi_199, isi_341, isi_342, ksh0_257, \
                         ksh0_258, ksh1_257, ksh1_258, ksi_339, ksi_341, \
                         ksi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_15 * isi_341[k]
                   + f_10 * ksh0_257[k]
                   - f_11 * ksh1_257[k]
                   + f_3 * pc_x[k] * ksi_341[k];

        t_438[k] = f_15 * isi_342[k]
                   + f_8 * ksh0_258[k]
                   - f_9 * ksh1_258[k]
                   + f_3 * pc_x[k] * ksi_342[k];

        t_439[k] = f_14 * isi_199[k]
                   + f_3 * pc_z[k] * ksi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, isi_229, isi_345, isi_346, ksh0_261, \
                         ksh0_262, ksh1_261, ksh1_262, ksi_341, ksi_345, \
                         ksi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * isi_229[k]
                   + f_3 * pc_y[k] * ksi_341[k];

        t_441[k] = f_15 * isi_345[k]
                   + f_8 * ksh0_261[k]
                   - f_9 * ksh1_261[k]
                   + f_3 * pc_x[k] * ksi_345[k];

        t_442[k] = f_15 * isi_346[k]
                   + f_6 * ksh0_262[k]
                   - f_7 * ksh1_262[k]
                   + f_3 * pc_x[k] * ksi_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, isi_202, isi_233, isi_348, \
                         ksh0_264, ksh1_264, ksi_342, ksi_345, \
                         ksi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * isi_202[k]
                   + f_3 * pc_z[k] * ksi_342[k];

        t_444[k] = f_15 * isi_348[k]
                   + f_6 * ksh0_264[k]
                   - f_7 * ksh1_264[k]
                   + f_3 * pc_x[k] * ksi_348[k];

        t_445[k] = f_14 * isi_233[k]
                   + f_3 * pc_y[k] * ksi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, isi_206, isi_350, isi_351, ksh0_266, \
                         ksh0_267, ksh1_266, ksh1_267, ksi_346, ksi_350, \
                         ksi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_15 * isi_350[k]
                   + f_6 * ksh0_266[k]
                   - f_7 * ksh1_266[k]
                   + f_3 * pc_x[k] * ksi_350[k];

        t_447[k] = f_15 * isi_351[k]
                   + f_4 * ksh0_267[k]
                   - f_5 * ksh1_267[k]
                   + f_3 * pc_x[k] * ksi_351[k];

        t_448[k] = f_14 * isi_206[k]
                   + f_3 * pc_z[k] * ksi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, isi_238, isi_353, isi_354, ksh0_269, \
                         ksh0_270, ksh1_269, ksh1_270, ksi_350, ksi_353, \
                         ksi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_15 * isi_353[k]
                   + f_4 * ksh0_269[k]
                   - f_5 * ksh1_269[k]
                   + f_3 * pc_x[k] * ksi_353[k];

        t_450[k] = f_15 * isi_354[k]
                   + f_4 * ksh0_270[k]
                   - f_5 * ksh1_270[k]
                   + f_3 * pc_x[k] * ksi_354[k];

        t_451[k] = f_14 * isi_238[k]
                   + f_3 * pc_y[k] * ksi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, isi_356, isi_357, isi_358, isi_359, \
                         ksh0_272, ksh1_272, ksi_356, ksi_357, ksi_358, \
                         ksi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_15 * isi_356[k]
                   + f_4 * ksh0_272[k]
                   - f_5 * ksh1_272[k]
                   + f_3 * pc_x[k] * ksi_356[k];

        t_453[k] = f_15 * isi_357[k]
                   + f_3 * pc_x[k] * ksi_357[k];

        t_454[k] = f_15 * isi_358[k]
                   + f_3 * pc_x[k] * ksi_358[k];

        t_455[k] = f_15 * isi_359[k]
                   + f_3 * pc_x[k] * ksi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, isi_360, isi_361, isi_362, isi_363, \
                         ksi_360, ksi_361, ksi_362, ksi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_15 * isi_360[k]
                   + f_3 * pc_x[k] * ksi_360[k];

        t_457[k] = f_15 * isi_361[k]
                   + f_3 * pc_x[k] * ksi_361[k];

        t_458[k] = f_15 * isi_362[k]
                   + f_3 * pc_x[k] * ksi_362[k];

        t_459[k] = f_15 * isi_363[k]
                   + f_3 * pc_x[k] * ksi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, isi_217, isi_245, isi_247, ksh0_267, \
                         ksh0_269, ksh1_267, ksh1_269, ksi_357, \
                         ksi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * isi_245[k]
                   + f_1 * ksh0_267[k]
                   - f_2 * ksh1_267[k]
                   + f_3 * pc_y[k] * ksi_357[k];

        t_461[k] = f_14 * isi_217[k]
                   + f_3 * pc_z[k] * ksi_357[k];

        t_462[k] = f_14 * isi_247[k]
                   + f_10 * ksh0_269[k]
                   - f_11 * ksh1_269[k]
                   + f_3 * pc_y[k] * ksi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, isi_248, isi_249, isi_250, ksh0_270, \
                         ksh0_271, ksh0_272, ksh1_270, ksh1_271, ksh1_272, ksi_360, ksi_361, \
                         ksi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * isi_248[k]
                   + f_8 * ksh0_270[k]
                   - f_9 * ksh1_270[k]
                   + f_3 * pc_y[k] * ksi_360[k];

        t_464[k] = f_14 * isi_249[k]
                   + f_6 * ksh0_271[k]
                   - f_7 * ksh1_271[k]
                   + f_3 * pc_y[k] * ksi_361[k];

        t_465[k] = f_14 * isi_250[k]
                   + f_4 * ksh0_272[k]
                   - f_5 * ksh1_272[k]
                   + f_3 * pc_y[k] * ksi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_y, pc_y, pc_z, isk0_324, isi_223, \
                         isi_251, isi_252, isk1_324, ksh0_272, ksh1_272, ksi_363, \
                         ksi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * isi_251[k]
                   + f_3 * pc_y[k] * ksi_363[k];

        t_467[k] = f_14 * isi_223[k]
                   + f_1 * ksh0_272[k]
                   - f_2 * ksh1_272[k]
                   + f_3 * pc_z[k] * ksi_363[k];

        t_468[k] = pa_y[k] * isk0_324[k]
                   - f_12 * pc_y[k] * isk1_324[k];

        t_469[k] = f_13 * isi_252[k]
                   + f_3 * pc_y[k] * ksi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_y, pc_y, pc_z, isk0_327, isk0_329, \
                         isi_224, isi_253, isi_254, isk1_327, isk1_329, ksi_364, \
                         ksi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * isi_224[k]
                   + f_3 * pc_z[k] * ksi_364[k];

        t_471[k] = pa_y[k] * isk0_327[k]
                   + f_14 * isi_253[k]
                   - f_12 * pc_y[k] * isk1_327[k];

        t_472[k] = f_13 * isi_254[k]
                   + f_3 * pc_y[k] * ksi_366[k];

        t_473[k] = pa_y[k] * isk0_329[k]
                   - f_12 * pc_y[k] * isk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_y, pc_y, pc_z, isk0_330, isk0_333, \
                         isi_227, isi_255, isi_257, isk1_330, isk1_333, ksi_367, \
                         ksi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_y[k] * isk0_330[k]
                   + f_15 * isi_255[k]
                   - f_12 * pc_y[k] * isk1_330[k];

        t_475[k] = f_15 * isi_227[k]
                   + f_3 * pc_z[k] * ksi_367[k];

        t_476[k] = f_13 * isi_257[k]
                   + f_3 * pc_y[k] * ksi_369[k];

        t_477[k] = pa_y[k] * isk0_333[k]
                   - f_12 * pc_y[k] * isk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_y, pc_y, pc_z, isk0_334, isk0_336, isi_230, \
                         isi_258, isi_260, isk1_334, isk1_336, \
                         ksi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pa_y[k] * isk0_334[k]
                   + f_16 * isi_258[k]
                   - f_12 * pc_y[k] * isk1_334[k];

        t_479[k] = f_15 * isi_230[k]
                   + f_3 * pc_z[k] * ksi_370[k];

        t_480[k] = pa_y[k] * isk0_336[k]
                   + f_14 * isi_260[k]
                   - f_12 * pc_y[k] * isk1_336[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
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

    const auto *isk0_338 = buffer.data(isk0 + 338);
    const auto *isk0_339 = buffer.data(isk0 + 339);
    const auto *isk0_341 = buffer.data(isk0 + 341);
    const auto *isk0_342 = buffer.data(isk0 + 342);
    const auto *isk0_344 = buffer.data(isk0 + 344);
    const auto *isk0_359 = buffer.data(isk0 + 359);
    const auto *isk0_360 = buffer.data(isk0 + 360);
    const auto *isk0_363 = buffer.data(isk0 + 363);
    const auto *isk0_366 = buffer.data(isk0 + 366);
    const auto *isk0_370 = buffer.data(isk0 + 370);
    const auto *isk0_372 = buffer.data(isk0 + 372);
    const auto *isk0_375 = buffer.data(isk0 + 375);
    const auto *isk0_377 = buffer.data(isk0 + 377);
    const auto *isk0_378 = buffer.data(isk0 + 378);

    const auto *isi_234 = buffer.data(isi + 234);
    const auto *isi_245 = buffer.data(isi + 245);
    const auto *isi_252 = buffer.data(isi + 252);
    const auto *isi_261 = buffer.data(isi + 261);
    const auto *isi_262 = buffer.data(isi + 262);
    const auto *isi_264 = buffer.data(isi + 264);
    const auto *isi_265 = buffer.data(isi + 265);
    const auto *isi_266 = buffer.data(isi + 266);
    const auto *isi_273 = buffer.data(isi + 273);
    const auto *isi_275 = buffer.data(isi + 275);
    const auto *isi_276 = buffer.data(isi + 276);
    const auto *isi_277 = buffer.data(isi + 277);
    const auto *isi_278 = buffer.data(isi + 278);
    const auto *isi_279 = buffer.data(isi + 279);
    const auto *isi_280 = buffer.data(isi + 280);
    const auto *isi_283 = buffer.data(isi + 283);
    const auto *isi_285 = buffer.data(isi + 285);
    const auto *isi_286 = buffer.data(isi + 286);
    const auto *isi_287 = buffer.data(isi + 287);
    const auto *isi_289 = buffer.data(isi + 289);
    const auto *isi_290 = buffer.data(isi + 290);
    const auto *isi_291 = buffer.data(isi + 291);
    const auto *isi_292 = buffer.data(isi + 292);
    const auto *isi_294 = buffer.data(isi + 294);
    const auto *isi_301 = buffer.data(isi + 301);
    const auto *isi_307 = buffer.data(isi + 307);
    const auto *isi_308 = buffer.data(isi + 308);
    const auto *isi_310 = buffer.data(isi + 310);
    const auto *isi_313 = buffer.data(isi + 313);
    const auto *isi_317 = buffer.data(isi + 317);
    const auto *isi_322 = buffer.data(isi + 322);
    const auto *isi_385 = buffer.data(isi + 385);
    const auto *isi_386 = buffer.data(isi + 386);
    const auto *isi_387 = buffer.data(isi + 387);
    const auto *isi_388 = buffer.data(isi + 388);
    const auto *isi_389 = buffer.data(isi + 389);
    const auto *isi_390 = buffer.data(isi + 390);
    const auto *isi_391 = buffer.data(isi + 391);
    const auto *isi_392 = buffer.data(isi + 392);
    const auto *isi_397 = buffer.data(isi + 397);
    const auto *isi_401 = buffer.data(isi + 401);
    const auto *isi_406 = buffer.data(isi + 406);
    const auto *isi_412 = buffer.data(isi + 412);
    const auto *isi_413 = buffer.data(isi + 413);
    const auto *isi_414 = buffer.data(isi + 414);
    const auto *isi_415 = buffer.data(isi + 415);
    const auto *isi_416 = buffer.data(isi + 416);
    const auto *isi_417 = buffer.data(isi + 417);
    const auto *isi_419 = buffer.data(isi + 419);
    const auto *isi_420 = buffer.data(isi + 420);
    const auto *isi_423 = buffer.data(isi + 423);
    const auto *isi_426 = buffer.data(isi + 426);
    const auto *isi_430 = buffer.data(isi + 430);
    const auto *isi_435 = buffer.data(isi + 435);
    const auto *isi_441 = buffer.data(isi + 441);
    const auto *isi_443 = buffer.data(isi + 443);
    const auto *isi_444 = buffer.data(isi + 444);
    const auto *isi_445 = buffer.data(isi + 445);
    const auto *isi_446 = buffer.data(isi + 446);
    const auto *isi_447 = buffer.data(isi + 447);
    const auto *isi_453 = buffer.data(isi + 453);
    const auto *isi_457 = buffer.data(isi + 457);
    const auto *isi_462 = buffer.data(isi + 462);
    const auto *isi_468 = buffer.data(isi + 468);
    const auto *isi_469 = buffer.data(isi + 469);
    const auto *isi_470 = buffer.data(isi + 470);
    const auto *isi_471 = buffer.data(isi + 471);

    const auto *isk1_338 = buffer.data(isk1 + 338);
    const auto *isk1_339 = buffer.data(isk1 + 339);
    const auto *isk1_341 = buffer.data(isk1 + 341);
    const auto *isk1_342 = buffer.data(isk1 + 342);
    const auto *isk1_344 = buffer.data(isk1 + 344);
    const auto *isk1_359 = buffer.data(isk1 + 359);
    const auto *isk1_360 = buffer.data(isk1 + 360);
    const auto *isk1_363 = buffer.data(isk1 + 363);
    const auto *isk1_366 = buffer.data(isk1 + 366);
    const auto *isk1_370 = buffer.data(isk1 + 370);
    const auto *isk1_372 = buffer.data(isk1 + 372);
    const auto *isk1_375 = buffer.data(isk1 + 375);
    const auto *isk1_377 = buffer.data(isk1 + 377);
    const auto *isk1_378 = buffer.data(isk1 + 378);

    const auto *ksh0_288 = buffer.data(ksh0 + 288);
    const auto *ksh0_290 = buffer.data(ksh0 + 290);
    const auto *ksh0_291 = buffer.data(ksh0 + 291);
    const auto *ksh0_292 = buffer.data(ksh0 + 292);
    const auto *ksh0_293 = buffer.data(ksh0 + 293);
    const auto *ksh0_294 = buffer.data(ksh0 + 294);
    const auto *ksh0_295 = buffer.data(ksh0 + 295);
    const auto *ksh0_296 = buffer.data(ksh0 + 296);
    const auto *ksh0_297 = buffer.data(ksh0 + 297);
    const auto *ksh0_298 = buffer.data(ksh0 + 298);
    const auto *ksh0_299 = buffer.data(ksh0 + 299);
    const auto *ksh0_300 = buffer.data(ksh0 + 300);
    const auto *ksh0_301 = buffer.data(ksh0 + 301);
    const auto *ksh0_302 = buffer.data(ksh0 + 302);
    const auto *ksh0_303 = buffer.data(ksh0 + 303);
    const auto *ksh0_308 = buffer.data(ksh0 + 308);
    const auto *ksh0_309 = buffer.data(ksh0 + 309);
    const auto *ksh0_310 = buffer.data(ksh0 + 310);
    const auto *ksh0_311 = buffer.data(ksh0 + 311);
    const auto *ksh0_312 = buffer.data(ksh0 + 312);
    const auto *ksh0_313 = buffer.data(ksh0 + 313);
    const auto *ksh0_314 = buffer.data(ksh0 + 314);
    const auto *ksh0_315 = buffer.data(ksh0 + 315);
    const auto *ksh0_317 = buffer.data(ksh0 + 317);
    const auto *ksh0_318 = buffer.data(ksh0 + 318);
    const auto *ksh0_320 = buffer.data(ksh0 + 320);
    const auto *ksh0_321 = buffer.data(ksh0 + 321);
    const auto *ksh0_322 = buffer.data(ksh0 + 322);
    const auto *ksh0_324 = buffer.data(ksh0 + 324);
    const auto *ksh0_325 = buffer.data(ksh0 + 325);
    const auto *ksh0_330 = buffer.data(ksh0 + 330);
    const auto *ksh0_331 = buffer.data(ksh0 + 331);
    const auto *ksh0_332 = buffer.data(ksh0 + 332);
    const auto *ksh0_333 = buffer.data(ksh0 + 333);
    const auto *ksh0_335 = buffer.data(ksh0 + 335);
    const auto *ksh0_341 = buffer.data(ksh0 + 341);
    const auto *ksh0_345 = buffer.data(ksh0 + 345);
    const auto *ksh0_350 = buffer.data(ksh0 + 350);
    const auto *ksh0_356 = buffer.data(ksh0 + 356);

    const auto *ksh1_288 = buffer.data(ksh1 + 288);
    const auto *ksh1_290 = buffer.data(ksh1 + 290);
    const auto *ksh1_291 = buffer.data(ksh1 + 291);
    const auto *ksh1_292 = buffer.data(ksh1 + 292);
    const auto *ksh1_293 = buffer.data(ksh1 + 293);
    const auto *ksh1_294 = buffer.data(ksh1 + 294);
    const auto *ksh1_295 = buffer.data(ksh1 + 295);
    const auto *ksh1_296 = buffer.data(ksh1 + 296);
    const auto *ksh1_297 = buffer.data(ksh1 + 297);
    const auto *ksh1_298 = buffer.data(ksh1 + 298);
    const auto *ksh1_299 = buffer.data(ksh1 + 299);
    const auto *ksh1_300 = buffer.data(ksh1 + 300);
    const auto *ksh1_301 = buffer.data(ksh1 + 301);
    const auto *ksh1_302 = buffer.data(ksh1 + 302);
    const auto *ksh1_303 = buffer.data(ksh1 + 303);
    const auto *ksh1_308 = buffer.data(ksh1 + 308);
    const auto *ksh1_309 = buffer.data(ksh1 + 309);
    const auto *ksh1_310 = buffer.data(ksh1 + 310);
    const auto *ksh1_311 = buffer.data(ksh1 + 311);
    const auto *ksh1_312 = buffer.data(ksh1 + 312);
    const auto *ksh1_313 = buffer.data(ksh1 + 313);
    const auto *ksh1_314 = buffer.data(ksh1 + 314);
    const auto *ksh1_315 = buffer.data(ksh1 + 315);
    const auto *ksh1_317 = buffer.data(ksh1 + 317);
    const auto *ksh1_318 = buffer.data(ksh1 + 318);
    const auto *ksh1_320 = buffer.data(ksh1 + 320);
    const auto *ksh1_321 = buffer.data(ksh1 + 321);
    const auto *ksh1_322 = buffer.data(ksh1 + 322);
    const auto *ksh1_324 = buffer.data(ksh1 + 324);
    const auto *ksh1_325 = buffer.data(ksh1 + 325);
    const auto *ksh1_330 = buffer.data(ksh1 + 330);
    const auto *ksh1_331 = buffer.data(ksh1 + 331);
    const auto *ksh1_332 = buffer.data(ksh1 + 332);
    const auto *ksh1_333 = buffer.data(ksh1 + 333);
    const auto *ksh1_335 = buffer.data(ksh1 + 335);
    const auto *ksh1_341 = buffer.data(ksh1 + 341);
    const auto *ksh1_345 = buffer.data(ksh1 + 345);
    const auto *ksh1_350 = buffer.data(ksh1 + 350);
    const auto *ksh1_356 = buffer.data(ksh1 + 356);

    const auto *ksi_373 = buffer.data(ksi + 373);
    const auto *ksi_374 = buffer.data(ksi + 374);
    const auto *ksi_378 = buffer.data(ksi + 378);
    const auto *ksi_385 = buffer.data(ksi + 385);
    const auto *ksi_386 = buffer.data(ksi + 386);
    const auto *ksi_387 = buffer.data(ksi + 387);
    const auto *ksi_388 = buffer.data(ksi + 388);
    const auto *ksi_389 = buffer.data(ksi + 389);
    const auto *ksi_390 = buffer.data(ksi + 390);
    const auto *ksi_391 = buffer.data(ksi + 391);
    const auto *ksi_392 = buffer.data(ksi + 392);
    const auto *ksi_393 = buffer.data(ksi + 393);
    const auto *ksi_394 = buffer.data(ksi + 394);
    const auto *ksi_395 = buffer.data(ksi + 395);
    const auto *ksi_396 = buffer.data(ksi + 396);
    const auto *ksi_397 = buffer.data(ksi + 397);
    const auto *ksi_398 = buffer.data(ksi + 398);
    const auto *ksi_399 = buffer.data(ksi + 399);
    const auto *ksi_400 = buffer.data(ksi + 400);
    const auto *ksi_401 = buffer.data(ksi + 401);
    const auto *ksi_402 = buffer.data(ksi + 402);
    const auto *ksi_403 = buffer.data(ksi + 403);
    const auto *ksi_404 = buffer.data(ksi + 404);
    const auto *ksi_405 = buffer.data(ksi + 405);
    const auto *ksi_406 = buffer.data(ksi + 406);
    const auto *ksi_412 = buffer.data(ksi + 412);
    const auto *ksi_413 = buffer.data(ksi + 413);
    const auto *ksi_414 = buffer.data(ksi + 414);
    const auto *ksi_415 = buffer.data(ksi + 415);
    const auto *ksi_416 = buffer.data(ksi + 416);
    const auto *ksi_417 = buffer.data(ksi + 417);
    const auto *ksi_418 = buffer.data(ksi + 418);
    const auto *ksi_419 = buffer.data(ksi + 419);
    const auto *ksi_420 = buffer.data(ksi + 420);
    const auto *ksi_421 = buffer.data(ksi + 421);
    const auto *ksi_422 = buffer.data(ksi + 422);
    const auto *ksi_423 = buffer.data(ksi + 423);
    const auto *ksi_425 = buffer.data(ksi + 425);
    const auto *ksi_426 = buffer.data(ksi + 426);
    const auto *ksi_427 = buffer.data(ksi + 427);
    const auto *ksi_429 = buffer.data(ksi + 429);
    const auto *ksi_430 = buffer.data(ksi + 430);
    const auto *ksi_431 = buffer.data(ksi + 431);
    const auto *ksi_432 = buffer.data(ksi + 432);
    const auto *ksi_434 = buffer.data(ksi + 434);
    const auto *ksi_435 = buffer.data(ksi + 435);
    const auto *ksi_441 = buffer.data(ksi + 441);
    const auto *ksi_442 = buffer.data(ksi + 442);
    const auto *ksi_443 = buffer.data(ksi + 443);
    const auto *ksi_444 = buffer.data(ksi + 444);
    const auto *ksi_445 = buffer.data(ksi + 445);
    const auto *ksi_446 = buffer.data(ksi + 446);
    const auto *ksi_447 = buffer.data(ksi + 447);
    const auto *ksi_448 = buffer.data(ksi + 448);
    const auto *ksi_450 = buffer.data(ksi + 450);
    const auto *ksi_451 = buffer.data(ksi + 451);
    const auto *ksi_453 = buffer.data(ksi + 453);
    const auto *ksi_454 = buffer.data(ksi + 454);
    const auto *ksi_457 = buffer.data(ksi + 457);
    const auto *ksi_458 = buffer.data(ksi + 458);
    const auto *ksi_462 = buffer.data(ksi + 462);
    const auto *ksi_468 = buffer.data(ksi + 468);
    const auto *ksi_469 = buffer.data(ksi + 469);
    const auto *ksi_470 = buffer.data(ksi + 470);
    const auto *ksi_471 = buffer.data(ksi + 471);

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_y, pc_y, pc_z, isk0_338, isk0_339, \
                         isi_234, isi_261, isi_262, isk1_338, isk1_339, ksi_373, \
                         ksi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * isi_261[k]
                   + f_3 * pc_y[k] * ksi_373[k];

        t_482[k] = pa_y[k] * isk0_338[k]
                   - f_12 * pc_y[k] * isk1_338[k];

        t_483[k] = pa_y[k] * isk0_339[k]
                   + f_17 * isi_262[k]
                   - f_12 * pc_y[k] * isk1_339[k];

        t_484[k] = f_15 * isi_234[k]
                   + f_3 * pc_z[k] * ksi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_y, pc_y, isk0_341, isk0_342, isk0_344, \
                         isi_264, isi_265, isi_266, isk1_341, isk1_342, isk1_344, \
                         ksi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_y[k] * isk0_341[k]
                   + f_15 * isi_264[k]
                   - f_12 * pc_y[k] * isk1_341[k];

        t_486[k] = pa_y[k] * isk0_342[k]
                   + f_14 * isi_265[k]
                   - f_12 * pc_y[k] * isk1_342[k];

        t_487[k] = f_13 * isi_266[k]
                   + f_3 * pc_y[k] * ksi_378[k];

        t_488[k] = pa_y[k] * isk0_344[k]
                   - f_12 * pc_y[k] * isk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, isi_385, isi_386, isi_387, \
                         isi_388, isi_389, ksi_385, ksi_386, ksi_387, ksi_388, \
                         ksi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_15 * isi_385[k]
                   + f_3 * pc_x[k] * ksi_385[k];

        t_490[k] = f_15 * isi_386[k]
                   + f_3 * pc_x[k] * ksi_386[k];

        t_491[k] = f_15 * isi_387[k]
                   + f_3 * pc_x[k] * ksi_387[k];

        t_492[k] = f_15 * isi_388[k]
                   + f_3 * pc_x[k] * ksi_388[k];

        t_493[k] = f_15 * isi_389[k]
                   + f_3 * pc_x[k] * ksi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, isi_245, isi_273, \
                         isi_390, isi_391, ksh0_288, ksh1_288, ksi_385, ksi_390, \
                         ksi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_15 * isi_390[k]
                   + f_3 * pc_x[k] * ksi_390[k];

        t_495[k] = f_15 * isi_391[k]
                   + f_3 * pc_x[k] * ksi_391[k];

        t_496[k] = f_13 * isi_273[k]
                   + f_1 * ksh0_288[k]
                   - f_2 * ksh1_288[k]
                   + f_3 * pc_y[k] * ksi_385[k];

        t_497[k] = f_15 * isi_245[k]
                   + f_3 * pc_z[k] * ksi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, isi_275, isi_276, isi_277, ksh0_290, \
                         ksh0_291, ksh0_292, ksh1_290, ksh1_291, ksh1_292, ksi_387, ksi_388, \
                         ksi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * isi_275[k]
                   + f_10 * ksh0_290[k]
                   - f_11 * ksh1_290[k]
                   + f_3 * pc_y[k] * ksi_387[k];

        t_499[k] = f_13 * isi_276[k]
                   + f_8 * ksh0_291[k]
                   - f_9 * ksh1_291[k]
                   + f_3 * pc_y[k] * ksi_388[k];

        t_500[k] = f_13 * isi_277[k]
                   + f_6 * ksh0_292[k]
                   - f_7 * ksh1_292[k]
                   + f_3 * pc_y[k] * ksi_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pc_y, isk0_359, isi_278, isi_279, \
                         isk1_359, ksh0_293, ksh1_293, ksi_390, \
                         ksi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * isi_278[k]
                   + f_4 * ksh0_293[k]
                   - f_5 * ksh1_293[k]
                   + f_3 * pc_y[k] * ksi_390[k];

        t_502[k] = f_13 * isi_279[k]
                   + f_3 * pc_y[k] * ksi_391[k];

        t_503[k] = pa_y[k] * isk0_359[k]
                   - f_12 * pc_y[k] * isk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, isi_252, \
                         isi_392, ksh0_294, ksh1_294, ksi_392, ksi_393, \
                         ksi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_15 * isi_392[k]
                   + f_1 * ksh0_294[k]
                   - f_2 * ksh1_294[k]
                   + f_3 * pc_x[k] * ksi_392[k];

        t_505[k] = f_3 * pc_y[k] * ksi_392[k];

        t_506[k] = f_16 * isi_252[k]
                   + f_3 * pc_z[k] * ksi_392[k];

        t_507[k] = f_4 * ksh0_294[k]
                   - f_5 * ksh1_294[k]
                   + f_3 * pc_y[k] * ksi_393[k];

        t_508[k] = f_3 * pc_y[k] * ksi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, isi_397, ksh0_295, ksh0_296, \
                         ksh0_299, ksh1_295, ksh1_296, ksh1_299, ksi_395, ksi_396, \
                         ksi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_15 * isi_397[k]
                   + f_10 * ksh0_299[k]
                   - f_11 * ksh1_299[k]
                   + f_3 * pc_x[k] * ksi_397[k];

        t_510[k] = f_6 * ksh0_295[k]
                   - f_7 * ksh1_295[k]
                   + f_3 * pc_y[k] * ksi_395[k];

        t_511[k] = f_4 * ksh0_296[k]
                   - f_5 * ksh1_296[k]
                   + f_3 * pc_y[k] * ksi_396[k];

        t_512[k] = f_3 * pc_y[k] * ksi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_y, isi_401, ksh0_297, ksh0_298, \
                         ksh0_303, ksh1_297, ksh1_298, ksh1_303, ksi_398, ksi_399, \
                         ksi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_15 * isi_401[k]
                   + f_8 * ksh0_303[k]
                   - f_9 * ksh1_303[k]
                   + f_3 * pc_x[k] * ksi_401[k];

        t_514[k] = f_8 * ksh0_297[k]
                   - f_9 * ksh1_297[k]
                   + f_3 * pc_y[k] * ksi_398[k];

        t_515[k] = f_6 * ksh0_298[k]
                   - f_7 * ksh1_298[k]
                   + f_3 * pc_y[k] * ksi_399[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_x, pc_y, isi_406, ksh0_299, ksh0_308, \
                         ksh1_299, ksh1_308, ksi_400, ksi_401, \
                         ksi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_4 * ksh0_299[k]
                   - f_5 * ksh1_299[k]
                   + f_3 * pc_y[k] * ksi_400[k];

        t_517[k] = f_3 * pc_y[k] * ksi_401[k];

        t_518[k] = f_15 * isi_406[k]
                   + f_6 * ksh0_308[k]
                   - f_7 * ksh1_308[k]
                   + f_3 * pc_x[k] * ksi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, ksh0_300, ksh0_301, ksh0_302, ksh1_300, \
                         ksh1_301, ksh1_302, ksi_402, ksi_403, \
                         ksi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_10 * ksh0_300[k]
                   - f_11 * ksh1_300[k]
                   + f_3 * pc_y[k] * ksi_402[k];

        t_520[k] = f_8 * ksh0_301[k]
                   - f_9 * ksh1_301[k]
                   + f_3 * pc_y[k] * ksi_403[k];

        t_521[k] = f_6 * ksh0_302[k]
                   - f_7 * ksh1_302[k]
                   + f_3 * pc_y[k] * ksi_404[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, isi_412, isi_413, ksh0_303, \
                         ksh0_314, ksh1_303, ksh1_314, ksi_405, ksi_406, ksi_412, \
                         ksi_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_4 * ksh0_303[k]
                   - f_5 * ksh1_303[k]
                   + f_3 * pc_y[k] * ksi_405[k];

        t_523[k] = f_3 * pc_y[k] * ksi_406[k];

        t_524[k] = f_15 * isi_412[k]
                   + f_4 * ksh0_314[k]
                   - f_5 * ksh1_314[k]
                   + f_3 * pc_x[k] * ksi_412[k];

        t_525[k] = f_15 * isi_413[k]
                   + f_3 * pc_x[k] * ksi_413[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, pc_x, pc_y, isi_414, isi_415, \
                         isi_416, isi_417, ksi_412, ksi_414, ksi_415, ksi_416, \
                         ksi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_15 * isi_414[k]
                   + f_3 * pc_x[k] * ksi_414[k];

        t_527[k] = f_15 * isi_415[k]
                   + f_3 * pc_x[k] * ksi_415[k];

        t_528[k] = f_15 * isi_416[k]
                   + f_3 * pc_x[k] * ksi_416[k];

        t_529[k] = f_15 * isi_417[k]
                   + f_3 * pc_x[k] * ksi_417[k];

        t_530[k] = f_3 * pc_y[k] * ksi_412[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pc_x, pc_y, isi_419, ksh0_309, ksh0_310, \
                         ksh1_309, ksh1_310, ksi_413, ksi_414, \
                         ksi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_15 * isi_419[k]
                   + f_3 * pc_x[k] * ksi_419[k];

        t_532[k] = f_1 * ksh0_309[k]
                   - f_2 * ksh1_309[k]
                   + f_3 * pc_y[k] * ksi_413[k];

        t_533[k] = f_19 * ksh0_310[k]
                   - f_20 * ksh1_310[k]
                   + f_3 * pc_y[k] * ksi_414[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, ksh0_311, ksh0_312, ksh0_313, ksh1_311, \
                         ksh1_312, ksh1_313, ksi_415, ksi_416, \
                         ksi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * ksh0_311[k]
                   - f_11 * ksh1_311[k]
                   + f_3 * pc_y[k] * ksi_415[k];

        t_535[k] = f_8 * ksh0_312[k]
                   - f_9 * ksh1_312[k]
                   + f_3 * pc_y[k] * ksi_416[k];

        t_536[k] = f_6 * ksh0_313[k]
                   - f_7 * ksh1_313[k]
                   + f_3 * pc_y[k] * ksi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, isi_279, isi_420, \
                         ksh0_314, ksh0_315, ksh1_314, ksh1_315, ksi_418, ksi_419, \
                         ksi_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * ksh0_314[k]
                   - f_5 * ksh1_314[k]
                   + f_3 * pc_y[k] * ksi_418[k];

        t_538[k] = f_3 * pc_y[k] * ksi_419[k];

        t_539[k] = f_16 * isi_279[k]
                   + f_1 * ksh0_314[k]
                   - f_2 * ksh1_314[k]
                   + f_3 * pc_z[k] * ksi_419[k];

        t_540[k] = f_14 * isi_420[k]
                   + f_1 * ksh0_315[k]
                   - f_2 * ksh1_315[k]
                   + f_3 * pc_x[k] * ksi_420[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, isi_280, isi_423, \
                         ksh0_318, ksh1_318, ksi_420, ksi_421, \
                         ksi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_17 * isi_280[k]
                   + f_3 * pc_y[k] * ksi_420[k];

        t_542[k] = f_3 * pc_z[k] * ksi_420[k];

        t_543[k] = f_14 * isi_423[k]
                   + f_10 * ksh0_318[k]
                   - f_11 * ksh1_318[k]
                   + f_3 * pc_x[k] * ksi_423[k];

        t_544[k] = f_3 * pc_z[k] * ksi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, isi_426, ksh0_315, ksh0_321, \
                         ksh1_315, ksh1_321, ksi_422, ksi_423, \
                         ksi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * ksh0_315[k]
                   - f_5 * ksh1_315[k]
                   + f_3 * pc_z[k] * ksi_422[k];

        t_546[k] = f_14 * isi_426[k]
                   + f_8 * ksh0_321[k]
                   - f_9 * ksh1_321[k]
                   + f_3 * pc_x[k] * ksi_426[k];

        t_547[k] = f_3 * pc_z[k] * ksi_423[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pc_x, pc_y, pc_z, isi_285, isi_430, \
                         ksh0_317, ksh0_325, ksh1_317, ksh1_325, ksi_425, ksi_426, \
                         ksi_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * isi_285[k]
                   + f_3 * pc_y[k] * ksi_425[k];

        t_549[k] = f_6 * ksh0_317[k]
                   - f_7 * ksh1_317[k]
                   + f_3 * pc_z[k] * ksi_425[k];

        t_550[k] = f_14 * isi_430[k]
                   + f_6 * ksh0_325[k]
                   - f_7 * ksh1_325[k]
                   + f_3 * pc_x[k] * ksi_430[k];

        t_551[k] = f_3 * pc_z[k] * ksi_426[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_y, pc_z, isi_289, ksh0_318, ksh0_320, \
                         ksh1_318, ksh1_320, ksi_427, ksi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_4 * ksh0_318[k]
                   - f_5 * ksh1_318[k]
                   + f_3 * pc_z[k] * ksi_427[k];

        t_553[k] = f_17 * isi_289[k]
                   + f_3 * pc_y[k] * ksi_429[k];

        t_554[k] = f_8 * ksh0_320[k]
                   - f_9 * ksh1_320[k]
                   + f_3 * pc_z[k] * ksi_429[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pc_x, pc_z, isi_435, ksh0_321, ksh0_330, \
                         ksh1_321, ksh1_330, ksi_430, ksi_431, \
                         ksi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_14 * isi_435[k]
                   + f_4 * ksh0_330[k]
                   - f_5 * ksh1_330[k]
                   + f_3 * pc_x[k] * ksi_435[k];

        t_556[k] = f_3 * pc_z[k] * ksi_430[k];

        t_557[k] = f_4 * ksh0_321[k]
                   - f_5 * ksh1_321[k]
                   + f_3 * pc_z[k] * ksi_431[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pc_x, pc_y, pc_z, isi_294, isi_441, \
                         ksh0_322, ksh0_324, ksh1_322, ksh1_324, ksi_432, ksi_434, \
                         ksi_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_6 * ksh0_322[k]
                   - f_7 * ksh1_322[k]
                   + f_3 * pc_z[k] * ksi_432[k];

        t_559[k] = f_17 * isi_294[k]
                   + f_3 * pc_y[k] * ksi_434[k];

        t_560[k] = f_10 * ksh0_324[k]
                   - f_11 * ksh1_324[k]
                   + f_3 * pc_z[k] * ksi_434[k];

        t_561[k] = f_14 * isi_441[k]
                   + f_3 * pc_x[k] * ksi_441[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pc_x, pc_z, isi_443, isi_444, \
                         isi_445, isi_446, ksi_435, ksi_443, ksi_444, ksi_445, \
                         ksi_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_3 * pc_z[k] * ksi_435[k];

        t_563[k] = f_14 * isi_443[k]
                   + f_3 * pc_x[k] * ksi_443[k];

        t_564[k] = f_14 * isi_444[k]
                   + f_3 * pc_x[k] * ksi_444[k];

        t_565[k] = f_14 * isi_445[k]
                   + f_3 * pc_x[k] * ksi_445[k];

        t_566[k] = f_14 * isi_446[k]
                   + f_3 * pc_x[k] * ksi_446[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, isi_301, isi_447, \
                         ksh0_330, ksh1_330, ksi_441, ksi_442, \
                         ksi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_14 * isi_447[k]
                   + f_3 * pc_x[k] * ksi_447[k];

        t_568[k] = f_17 * isi_301[k]
                   + f_1 * ksh0_330[k]
                   - f_2 * ksh1_330[k]
                   + f_3 * pc_y[k] * ksi_441[k];

        t_569[k] = f_3 * pc_z[k] * ksi_441[k];

        t_570[k] = f_4 * ksh0_330[k]
                   - f_5 * ksh1_330[k]
                   + f_3 * pc_z[k] * ksi_442[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_z, ksh0_331, ksh0_332, ksh0_333, ksh1_331, \
                         ksh1_332, ksh1_333, ksi_443, ksi_444, \
                         ksi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_6 * ksh0_331[k]
                   - f_7 * ksh1_331[k]
                   + f_3 * pc_z[k] * ksi_443[k];

        t_572[k] = f_8 * ksh0_332[k]
                   - f_9 * ksh1_332[k]
                   + f_3 * pc_z[k] * ksi_444[k];

        t_573[k] = f_10 * ksh0_333[k]
                   - f_11 * ksh1_333[k]
                   + f_3 * pc_z[k] * ksi_445[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_z, pc_y, pc_z, isk0_360, isi_307, \
                         isi_308, isk1_360, ksh0_335, ksh1_335, ksi_447, \
                         ksi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * isi_307[k]
                   + f_3 * pc_y[k] * ksi_447[k];

        t_575[k] = f_1 * ksh0_335[k]
                   - f_2 * ksh1_335[k]
                   + f_3 * pc_z[k] * ksi_447[k];

        t_576[k] = pa_z[k] * isk0_360[k]
                   - f_12 * pc_z[k] * isk1_360[k];

        t_577[k] = f_16 * isi_308[k]
                   + f_3 * pc_y[k] * ksi_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_z, pc_y, pc_z, isk0_363, isi_280, isi_310, \
                         isk1_363, ksi_448, ksi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * isi_280[k]
                   + f_3 * pc_z[k] * ksi_448[k];

        t_579[k] = pa_z[k] * isk0_363[k]
                   - f_12 * pc_z[k] * isk1_363[k];

        t_580[k] = f_16 * isi_310[k]
                   + f_3 * pc_y[k] * ksi_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_z, pc_x, pc_z, isk0_366, isi_283, isi_453, \
                         isk1_366, ksh0_341, ksh1_341, ksi_451, \
                         ksi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_14 * isi_453[k]
                   + f_10 * ksh0_341[k]
                   - f_11 * ksh1_341[k]
                   + f_3 * pc_x[k] * ksi_453[k];

        t_582[k] = pa_z[k] * isk0_366[k]
                   - f_12 * pc_z[k] * isk1_366[k];

        t_583[k] = f_13 * isi_283[k]
                   + f_3 * pc_z[k] * ksi_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pa_z, pc_x, pc_y, pc_z, isk0_370, isi_313, \
                         isi_457, isk1_370, ksh0_345, ksh1_345, ksi_453, \
                         ksi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * isi_313[k]
                   + f_3 * pc_y[k] * ksi_453[k];

        t_585[k] = f_14 * isi_457[k]
                   + f_8 * ksh0_345[k]
                   - f_9 * ksh1_345[k]
                   + f_3 * pc_x[k] * ksi_457[k];

        t_586[k] = pa_z[k] * isk0_370[k]
                   - f_12 * pc_z[k] * isk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pa_z, pc_y, pc_z, isk0_372, isi_286, isi_287, \
                         isi_317, isk1_372, ksi_454, ksi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * isi_286[k]
                   + f_3 * pc_z[k] * ksi_454[k];

        t_588[k] = pa_z[k] * isk0_372[k]
                   + f_14 * isi_287[k]
                   - f_12 * pc_z[k] * isk1_372[k];

        t_589[k] = f_16 * isi_317[k]
                   + f_3 * pc_y[k] * ksi_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pa_z, pc_x, pc_z, isk0_375, isi_290, isi_462, \
                         isk1_375, ksh0_350, ksh1_350, ksi_458, \
                         ksi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_14 * isi_462[k]
                   + f_6 * ksh0_350[k]
                   - f_7 * ksh1_350[k]
                   + f_3 * pc_x[k] * ksi_462[k];

        t_591[k] = pa_z[k] * isk0_375[k]
                   - f_12 * pc_z[k] * isk1_375[k];

        t_592[k] = f_13 * isi_290[k]
                   + f_3 * pc_z[k] * ksi_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pa_z, pc_y, pc_z, isk0_377, isk0_378, isi_291, \
                         isi_292, isi_322, isk1_377, isk1_378, \
                         ksi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pa_z[k] * isk0_377[k]
                   + f_14 * isi_291[k]
                   - f_12 * pc_z[k] * isk1_377[k];

        t_594[k] = pa_z[k] * isk0_378[k]
                   + f_15 * isi_292[k]
                   - f_12 * pc_z[k] * isk1_378[k];

        t_595[k] = f_16 * isi_322[k]
                   + f_3 * pc_y[k] * ksi_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, isi_468, isi_469, isi_470, isi_471, \
                         ksh0_356, ksh1_356, ksi_468, ksi_469, ksi_470, \
                         ksi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_14 * isi_468[k]
                   + f_4 * ksh0_356[k]
                   - f_5 * ksh1_356[k]
                   + f_3 * pc_x[k] * ksi_468[k];

        t_597[k] = f_14 * isi_469[k]
                   + f_3 * pc_x[k] * ksi_469[k];

        t_598[k] = f_14 * isi_470[k]
                   + f_3 * pc_x[k] * ksi_470[k];

        t_599[k] = f_14 * isi_471[k]
                   + f_3 * pc_x[k] * ksi_471[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
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

    const auto *isk0_388 = buffer.data(isk0 + 388);
    const auto *isk0_504 = buffer.data(isk0 + 504);
    const auto *isk0_507 = buffer.data(isk0 + 507);
    const auto *isk0_509 = buffer.data(isk0 + 509);
    const auto *isk0_510 = buffer.data(isk0 + 510);
    const auto *isk0_513 = buffer.data(isk0 + 513);
    const auto *isk0_514 = buffer.data(isk0 + 514);
    const auto *isk0_516 = buffer.data(isk0 + 516);
    const auto *isk0_518 = buffer.data(isk0 + 518);
    const auto *isk0_519 = buffer.data(isk0 + 519);
    const auto *isk0_521 = buffer.data(isk0 + 521);
    const auto *isk0_522 = buffer.data(isk0 + 522);
    const auto *isk0_524 = buffer.data(isk0 + 524);

    const auto *isi_301 = buffer.data(isi + 301);
    const auto *isi_307 = buffer.data(isi + 307);
    const auto *isi_308 = buffer.data(isi + 308);
    const auto *isi_311 = buffer.data(isi + 311);
    const auto *isi_314 = buffer.data(isi + 314);
    const auto *isi_318 = buffer.data(isi + 318);
    const auto *isi_329 = buffer.data(isi + 329);
    const auto *isi_331 = buffer.data(isi + 331);
    const auto *isi_332 = buffer.data(isi + 332);
    const auto *isi_333 = buffer.data(isi + 333);
    const auto *isi_334 = buffer.data(isi + 334);
    const auto *isi_335 = buffer.data(isi + 335);
    const auto *isi_336 = buffer.data(isi + 336);
    const auto *isi_338 = buffer.data(isi + 338);
    const auto *isi_339 = buffer.data(isi + 339);
    const auto *isi_341 = buffer.data(isi + 341);
    const auto *isi_342 = buffer.data(isi + 342);
    const auto *isi_345 = buffer.data(isi + 345);
    const auto *isi_346 = buffer.data(isi + 346);
    const auto *isi_350 = buffer.data(isi + 350);
    const auto *isi_357 = buffer.data(isi + 357);
    const auto *isi_359 = buffer.data(isi + 359);
    const auto *isi_360 = buffer.data(isi + 360);
    const auto *isi_361 = buffer.data(isi + 361);
    const auto *isi_362 = buffer.data(isi + 362);
    const auto *isi_363 = buffer.data(isi + 363);
    const auto *isi_364 = buffer.data(isi + 364);
    const auto *isi_366 = buffer.data(isi + 366);
    const auto *isi_367 = buffer.data(isi + 367);
    const auto *isi_369 = buffer.data(isi + 369);
    const auto *isi_370 = buffer.data(isi + 370);
    const auto *isi_373 = buffer.data(isi + 373);
    const auto *isi_374 = buffer.data(isi + 374);
    const auto *isi_378 = buffer.data(isi + 378);
    const auto *isi_385 = buffer.data(isi + 385);
    const auto *isi_387 = buffer.data(isi + 387);
    const auto *isi_388 = buffer.data(isi + 388);
    const auto *isi_389 = buffer.data(isi + 389);
    const auto *isi_390 = buffer.data(isi + 390);
    const auto *isi_391 = buffer.data(isi + 391);
    const auto *isi_392 = buffer.data(isi + 392);
    const auto *isi_393 = buffer.data(isi + 393);
    const auto *isi_394 = buffer.data(isi + 394);
    const auto *isi_395 = buffer.data(isi + 395);
    const auto *isi_397 = buffer.data(isi + 397);
    const auto *isi_398 = buffer.data(isi + 398);
    const auto *isi_400 = buffer.data(isi + 400);
    const auto *isi_401 = buffer.data(isi + 401);
    const auto *isi_402 = buffer.data(isi + 402);
    const auto *isi_404 = buffer.data(isi + 404);
    const auto *isi_405 = buffer.data(isi + 405);
    const auto *isi_406 = buffer.data(isi + 406);
    const auto *isi_472 = buffer.data(isi + 472);
    const auto *isi_473 = buffer.data(isi + 473);
    const auto *isi_474 = buffer.data(isi + 474);
    const auto *isi_475 = buffer.data(isi + 475);
    const auto *isi_476 = buffer.data(isi + 476);
    const auto *isi_479 = buffer.data(isi + 479);
    const auto *isi_481 = buffer.data(isi + 481);
    const auto *isi_482 = buffer.data(isi + 482);
    const auto *isi_485 = buffer.data(isi + 485);
    const auto *isi_486 = buffer.data(isi + 486);
    const auto *isi_488 = buffer.data(isi + 488);
    const auto *isi_490 = buffer.data(isi + 490);
    const auto *isi_491 = buffer.data(isi + 491);
    const auto *isi_493 = buffer.data(isi + 493);
    const auto *isi_494 = buffer.data(isi + 494);
    const auto *isi_496 = buffer.data(isi + 496);
    const auto *isi_497 = buffer.data(isi + 497);
    const auto *isi_498 = buffer.data(isi + 498);
    const auto *isi_499 = buffer.data(isi + 499);
    const auto *isi_500 = buffer.data(isi + 500);
    const auto *isi_501 = buffer.data(isi + 501);
    const auto *isi_502 = buffer.data(isi + 502);
    const auto *isi_503 = buffer.data(isi + 503);
    const auto *isi_504 = buffer.data(isi + 504);
    const auto *isi_507 = buffer.data(isi + 507);
    const auto *isi_509 = buffer.data(isi + 509);
    const auto *isi_510 = buffer.data(isi + 510);
    const auto *isi_513 = buffer.data(isi + 513);
    const auto *isi_514 = buffer.data(isi + 514);
    const auto *isi_516 = buffer.data(isi + 516);
    const auto *isi_518 = buffer.data(isi + 518);
    const auto *isi_519 = buffer.data(isi + 519);
    const auto *isi_521 = buffer.data(isi + 521);
    const auto *isi_522 = buffer.data(isi + 522);
    const auto *isi_524 = buffer.data(isi + 524);
    const auto *isi_525 = buffer.data(isi + 525);
    const auto *isi_526 = buffer.data(isi + 526);
    const auto *isi_527 = buffer.data(isi + 527);
    const auto *isi_528 = buffer.data(isi + 528);
    const auto *isi_529 = buffer.data(isi + 529);
    const auto *isi_530 = buffer.data(isi + 530);
    const auto *isi_531 = buffer.data(isi + 531);
    const auto *isi_553 = buffer.data(isi + 553);
    const auto *isi_554 = buffer.data(isi + 554);
    const auto *isi_555 = buffer.data(isi + 555);
    const auto *isi_556 = buffer.data(isi + 556);
    const auto *isi_557 = buffer.data(isi + 557);

    const auto *isk1_388 = buffer.data(isk1 + 388);
    const auto *isk1_504 = buffer.data(isk1 + 504);
    const auto *isk1_507 = buffer.data(isk1 + 507);
    const auto *isk1_509 = buffer.data(isk1 + 509);
    const auto *isk1_510 = buffer.data(isk1 + 510);
    const auto *isk1_513 = buffer.data(isk1 + 513);
    const auto *isk1_514 = buffer.data(isk1 + 514);
    const auto *isk1_516 = buffer.data(isk1 + 516);
    const auto *isk1_518 = buffer.data(isk1 + 518);
    const auto *isk1_519 = buffer.data(isk1 + 519);
    const auto *isk1_521 = buffer.data(isk1 + 521);
    const auto *isk1_522 = buffer.data(isk1 + 522);
    const auto *isk1_524 = buffer.data(isk1 + 524);

    const auto *ksh0_353 = buffer.data(ksh0 + 353);
    const auto *ksh0_354 = buffer.data(ksh0 + 354);
    const auto *ksh0_355 = buffer.data(ksh0 + 355);
    const auto *ksh0_356 = buffer.data(ksh0 + 356);
    const auto *ksh0_357 = buffer.data(ksh0 + 357);
    const auto *ksh0_360 = buffer.data(ksh0 + 360);
    const auto *ksh0_362 = buffer.data(ksh0 + 362);
    const auto *ksh0_363 = buffer.data(ksh0 + 363);
    const auto *ksh0_366 = buffer.data(ksh0 + 366);
    const auto *ksh0_367 = buffer.data(ksh0 + 367);
    const auto *ksh0_369 = buffer.data(ksh0 + 369);
    const auto *ksh0_371 = buffer.data(ksh0 + 371);
    const auto *ksh0_372 = buffer.data(ksh0 + 372);
    const auto *ksh0_374 = buffer.data(ksh0 + 374);
    const auto *ksh0_375 = buffer.data(ksh0 + 375);
    const auto *ksh0_376 = buffer.data(ksh0 + 376);
    const auto *ksh0_377 = buffer.data(ksh0 + 377);
    const auto *ksh0_378 = buffer.data(ksh0 + 378);
    const auto *ksh0_381 = buffer.data(ksh0 + 381);
    const auto *ksh0_383 = buffer.data(ksh0 + 383);
    const auto *ksh0_384 = buffer.data(ksh0 + 384);
    const auto *ksh0_387 = buffer.data(ksh0 + 387);
    const auto *ksh0_388 = buffer.data(ksh0 + 388);
    const auto *ksh0_390 = buffer.data(ksh0 + 390);
    const auto *ksh0_392 = buffer.data(ksh0 + 392);
    const auto *ksh0_393 = buffer.data(ksh0 + 393);
    const auto *ksh0_395 = buffer.data(ksh0 + 395);
    const auto *ksh0_396 = buffer.data(ksh0 + 396);
    const auto *ksh0_397 = buffer.data(ksh0 + 397);
    const auto *ksh0_398 = buffer.data(ksh0 + 398);

    const auto *ksh1_353 = buffer.data(ksh1 + 353);
    const auto *ksh1_354 = buffer.data(ksh1 + 354);
    const auto *ksh1_355 = buffer.data(ksh1 + 355);
    const auto *ksh1_356 = buffer.data(ksh1 + 356);
    const auto *ksh1_357 = buffer.data(ksh1 + 357);
    const auto *ksh1_360 = buffer.data(ksh1 + 360);
    const auto *ksh1_362 = buffer.data(ksh1 + 362);
    const auto *ksh1_363 = buffer.data(ksh1 + 363);
    const auto *ksh1_366 = buffer.data(ksh1 + 366);
    const auto *ksh1_367 = buffer.data(ksh1 + 367);
    const auto *ksh1_369 = buffer.data(ksh1 + 369);
    const auto *ksh1_371 = buffer.data(ksh1 + 371);
    const auto *ksh1_372 = buffer.data(ksh1 + 372);
    const auto *ksh1_374 = buffer.data(ksh1 + 374);
    const auto *ksh1_375 = buffer.data(ksh1 + 375);
    const auto *ksh1_376 = buffer.data(ksh1 + 376);
    const auto *ksh1_377 = buffer.data(ksh1 + 377);
    const auto *ksh1_378 = buffer.data(ksh1 + 378);
    const auto *ksh1_381 = buffer.data(ksh1 + 381);
    const auto *ksh1_383 = buffer.data(ksh1 + 383);
    const auto *ksh1_384 = buffer.data(ksh1 + 384);
    const auto *ksh1_387 = buffer.data(ksh1 + 387);
    const auto *ksh1_388 = buffer.data(ksh1 + 388);
    const auto *ksh1_390 = buffer.data(ksh1 + 390);
    const auto *ksh1_392 = buffer.data(ksh1 + 392);
    const auto *ksh1_393 = buffer.data(ksh1 + 393);
    const auto *ksh1_395 = buffer.data(ksh1 + 395);
    const auto *ksh1_396 = buffer.data(ksh1 + 396);
    const auto *ksh1_397 = buffer.data(ksh1 + 397);
    const auto *ksh1_398 = buffer.data(ksh1 + 398);

    const auto *ksi_469 = buffer.data(ksi + 469);
    const auto *ksi_471 = buffer.data(ksi + 471);
    const auto *ksi_472 = buffer.data(ksi + 472);
    const auto *ksi_473 = buffer.data(ksi + 473);
    const auto *ksi_474 = buffer.data(ksi + 474);
    const auto *ksi_475 = buffer.data(ksi + 475);
    const auto *ksi_476 = buffer.data(ksi + 476);
    const auto *ksi_478 = buffer.data(ksi + 478);
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
    const auto *ksi_506 = buffer.data(ksi + 506);
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
    const auto *ksi_532 = buffer.data(ksi + 532);
    const auto *ksi_534 = buffer.data(ksi + 534);
    const auto *ksi_535 = buffer.data(ksi + 535);
    const auto *ksi_537 = buffer.data(ksi + 537);
    const auto *ksi_538 = buffer.data(ksi + 538);
    const auto *ksi_541 = buffer.data(ksi + 541);
    const auto *ksi_542 = buffer.data(ksi + 542);
    const auto *ksi_546 = buffer.data(ksi + 546);
    const auto *ksi_553 = buffer.data(ksi + 553);
    const auto *ksi_554 = buffer.data(ksi + 554);
    const auto *ksi_555 = buffer.data(ksi + 555);
    const auto *ksi_556 = buffer.data(ksi + 556);
    const auto *ksi_557 = buffer.data(ksi + 557);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, isi_472, isi_473, isi_474, isi_475, \
                         ksi_472, ksi_473, ksi_474, ksi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_14 * isi_472[k]
                   + f_3 * pc_x[k] * ksi_472[k];

        t_601[k] = f_14 * isi_473[k]
                   + f_3 * pc_x[k] * ksi_473[k];

        t_602[k] = f_14 * isi_474[k]
                   + f_3 * pc_x[k] * ksi_474[k];

        t_603[k] = f_14 * isi_475[k]
                   + f_3 * pc_x[k] * ksi_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pa_z, pc_y, pc_z, isk0_388, isi_301, isi_331, \
                         isk1_388, ksh0_353, ksh1_353, ksi_469, \
                         ksi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_z[k] * isk0_388[k]
                   - f_12 * pc_z[k] * isk1_388[k];

        t_605[k] = f_13 * isi_301[k]
                   + f_3 * pc_z[k] * ksi_469[k];

        t_606[k] = f_16 * isi_331[k]
                   + f_10 * ksh0_353[k]
                   - f_11 * ksh1_353[k]
                   + f_3 * pc_y[k] * ksi_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, isi_332, isi_333, isi_334, ksh0_354, \
                         ksh0_355, ksh0_356, ksh1_354, ksh1_355, ksh1_356, ksi_472, ksi_473, \
                         ksi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * isi_332[k]
                   + f_8 * ksh0_354[k]
                   - f_9 * ksh1_354[k]
                   + f_3 * pc_y[k] * ksi_472[k];

        t_608[k] = f_16 * isi_333[k]
                   + f_6 * ksh0_355[k]
                   - f_7 * ksh1_355[k]
                   + f_3 * pc_y[k] * ksi_473[k];

        t_609[k] = f_16 * isi_334[k]
                   + f_4 * ksh0_356[k]
                   - f_5 * ksh1_356[k]
                   + f_3 * pc_y[k] * ksi_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, isi_307, isi_335, isi_476, \
                         ksh0_356, ksh0_357, ksh1_356, ksh1_357, ksi_475, \
                         ksi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * isi_335[k]
                   + f_3 * pc_y[k] * ksi_475[k];

        t_611[k] = f_13 * isi_307[k]
                   + f_1 * ksh0_356[k]
                   - f_2 * ksh1_356[k]
                   + f_3 * pc_z[k] * ksi_475[k];

        t_612[k] = f_14 * isi_476[k]
                   + f_1 * ksh0_357[k]
                   - f_2 * ksh1_357[k]
                   + f_3 * pc_x[k] * ksi_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, isi_308, isi_336, \
                         isi_338, isi_479, ksh0_360, ksh1_360, ksi_476, ksi_478, \
                         ksi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * isi_336[k]
                   + f_3 * pc_y[k] * ksi_476[k];

        t_614[k] = f_14 * isi_308[k]
                   + f_3 * pc_z[k] * ksi_476[k];

        t_615[k] = f_14 * isi_479[k]
                   + f_10 * ksh0_360[k]
                   - f_11 * ksh1_360[k]
                   + f_3 * pc_x[k] * ksi_479[k];

        t_616[k] = f_15 * isi_338[k]
                   + f_3 * pc_y[k] * ksi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, isi_311, isi_481, isi_482, ksh0_362, \
                         ksh0_363, ksh1_362, ksh1_363, ksi_479, ksi_481, \
                         ksi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_14 * isi_481[k]
                   + f_10 * ksh0_362[k]
                   - f_11 * ksh1_362[k]
                   + f_3 * pc_x[k] * ksi_481[k];

        t_618[k] = f_14 * isi_482[k]
                   + f_8 * ksh0_363[k]
                   - f_9 * ksh1_363[k]
                   + f_3 * pc_x[k] * ksi_482[k];

        t_619[k] = f_14 * isi_311[k]
                   + f_3 * pc_z[k] * ksi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, isi_341, isi_485, isi_486, ksh0_366, \
                         ksh0_367, ksh1_366, ksh1_367, ksi_481, ksi_485, \
                         ksi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * isi_341[k]
                   + f_3 * pc_y[k] * ksi_481[k];

        t_621[k] = f_14 * isi_485[k]
                   + f_8 * ksh0_366[k]
                   - f_9 * ksh1_366[k]
                   + f_3 * pc_x[k] * ksi_485[k];

        t_622[k] = f_14 * isi_486[k]
                   + f_6 * ksh0_367[k]
                   - f_7 * ksh1_367[k]
                   + f_3 * pc_x[k] * ksi_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, isi_314, isi_345, isi_488, \
                         ksh0_369, ksh1_369, ksi_482, ksi_485, \
                         ksi_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * isi_314[k]
                   + f_3 * pc_z[k] * ksi_482[k];

        t_624[k] = f_14 * isi_488[k]
                   + f_6 * ksh0_369[k]
                   - f_7 * ksh1_369[k]
                   + f_3 * pc_x[k] * ksi_488[k];

        t_625[k] = f_15 * isi_345[k]
                   + f_3 * pc_y[k] * ksi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, isi_318, isi_490, isi_491, ksh0_371, \
                         ksh0_372, ksh1_371, ksh1_372, ksi_486, ksi_490, \
                         ksi_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_14 * isi_490[k]
                   + f_6 * ksh0_371[k]
                   - f_7 * ksh1_371[k]
                   + f_3 * pc_x[k] * ksi_490[k];

        t_627[k] = f_14 * isi_491[k]
                   + f_4 * ksh0_372[k]
                   - f_5 * ksh1_372[k]
                   + f_3 * pc_x[k] * ksi_491[k];

        t_628[k] = f_14 * isi_318[k]
                   + f_3 * pc_z[k] * ksi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, isi_350, isi_493, isi_494, ksh0_374, \
                         ksh0_375, ksh1_374, ksh1_375, ksi_490, ksi_493, \
                         ksi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_14 * isi_493[k]
                   + f_4 * ksh0_374[k]
                   - f_5 * ksh1_374[k]
                   + f_3 * pc_x[k] * ksi_493[k];

        t_630[k] = f_14 * isi_494[k]
                   + f_4 * ksh0_375[k]
                   - f_5 * ksh1_375[k]
                   + f_3 * pc_x[k] * ksi_494[k];

        t_631[k] = f_15 * isi_350[k]
                   + f_3 * pc_y[k] * ksi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, isi_496, isi_497, isi_498, isi_499, \
                         ksh0_377, ksh1_377, ksi_496, ksi_497, ksi_498, \
                         ksi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_14 * isi_496[k]
                   + f_4 * ksh0_377[k]
                   - f_5 * ksh1_377[k]
                   + f_3 * pc_x[k] * ksi_496[k];

        t_633[k] = f_14 * isi_497[k]
                   + f_3 * pc_x[k] * ksi_497[k];

        t_634[k] = f_14 * isi_498[k]
                   + f_3 * pc_x[k] * ksi_498[k];

        t_635[k] = f_14 * isi_499[k]
                   + f_3 * pc_x[k] * ksi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, isi_500, isi_501, isi_502, isi_503, \
                         ksi_500, ksi_501, ksi_502, ksi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_14 * isi_500[k]
                   + f_3 * pc_x[k] * ksi_500[k];

        t_637[k] = f_14 * isi_501[k]
                   + f_3 * pc_x[k] * ksi_501[k];

        t_638[k] = f_14 * isi_502[k]
                   + f_3 * pc_x[k] * ksi_502[k];

        t_639[k] = f_14 * isi_503[k]
                   + f_3 * pc_x[k] * ksi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, isi_329, isi_357, isi_359, ksh0_372, \
                         ksh0_374, ksh1_372, ksh1_374, ksi_497, \
                         ksi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * isi_357[k]
                   + f_1 * ksh0_372[k]
                   - f_2 * ksh1_372[k]
                   + f_3 * pc_y[k] * ksi_497[k];

        t_641[k] = f_14 * isi_329[k]
                   + f_3 * pc_z[k] * ksi_497[k];

        t_642[k] = f_15 * isi_359[k]
                   + f_10 * ksh0_374[k]
                   - f_11 * ksh1_374[k]
                   + f_3 * pc_y[k] * ksi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, isi_360, isi_361, isi_362, ksh0_375, \
                         ksh0_376, ksh0_377, ksh1_375, ksh1_376, ksh1_377, ksi_500, ksi_501, \
                         ksi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * isi_360[k]
                   + f_8 * ksh0_375[k]
                   - f_9 * ksh1_375[k]
                   + f_3 * pc_y[k] * ksi_500[k];

        t_644[k] = f_15 * isi_361[k]
                   + f_6 * ksh0_376[k]
                   - f_7 * ksh1_376[k]
                   + f_3 * pc_y[k] * ksi_501[k];

        t_645[k] = f_15 * isi_362[k]
                   + f_4 * ksh0_377[k]
                   - f_5 * ksh1_377[k]
                   + f_3 * pc_y[k] * ksi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, isi_335, isi_363, isi_504, \
                         ksh0_377, ksh0_378, ksh1_377, ksh1_378, ksi_503, \
                         ksi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * isi_363[k]
                   + f_3 * pc_y[k] * ksi_503[k];

        t_647[k] = f_14 * isi_335[k]
                   + f_1 * ksh0_377[k]
                   - f_2 * ksh1_377[k]
                   + f_3 * pc_z[k] * ksi_503[k];

        t_648[k] = f_14 * isi_504[k]
                   + f_1 * ksh0_378[k]
                   - f_2 * ksh1_378[k]
                   + f_3 * pc_x[k] * ksi_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, isi_336, isi_364, \
                         isi_366, isi_507, ksh0_381, ksh1_381, ksi_504, ksi_506, \
                         ksi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * isi_364[k]
                   + f_3 * pc_y[k] * ksi_504[k];

        t_650[k] = f_15 * isi_336[k]
                   + f_3 * pc_z[k] * ksi_504[k];

        t_651[k] = f_14 * isi_507[k]
                   + f_10 * ksh0_381[k]
                   - f_11 * ksh1_381[k]
                   + f_3 * pc_x[k] * ksi_507[k];

        t_652[k] = f_14 * isi_366[k]
                   + f_3 * pc_y[k] * ksi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, isi_339, isi_509, isi_510, ksh0_383, \
                         ksh0_384, ksh1_383, ksh1_384, ksi_507, ksi_509, \
                         ksi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_14 * isi_509[k]
                   + f_10 * ksh0_383[k]
                   - f_11 * ksh1_383[k]
                   + f_3 * pc_x[k] * ksi_509[k];

        t_654[k] = f_14 * isi_510[k]
                   + f_8 * ksh0_384[k]
                   - f_9 * ksh1_384[k]
                   + f_3 * pc_x[k] * ksi_510[k];

        t_655[k] = f_15 * isi_339[k]
                   + f_3 * pc_z[k] * ksi_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, isi_369, isi_513, isi_514, ksh0_387, \
                         ksh0_388, ksh1_387, ksh1_388, ksi_509, ksi_513, \
                         ksi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * isi_369[k]
                   + f_3 * pc_y[k] * ksi_509[k];

        t_657[k] = f_14 * isi_513[k]
                   + f_8 * ksh0_387[k]
                   - f_9 * ksh1_387[k]
                   + f_3 * pc_x[k] * ksi_513[k];

        t_658[k] = f_14 * isi_514[k]
                   + f_6 * ksh0_388[k]
                   - f_7 * ksh1_388[k]
                   + f_3 * pc_x[k] * ksi_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, isi_342, isi_373, isi_516, \
                         ksh0_390, ksh1_390, ksi_510, ksi_513, \
                         ksi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * isi_342[k]
                   + f_3 * pc_z[k] * ksi_510[k];

        t_660[k] = f_14 * isi_516[k]
                   + f_6 * ksh0_390[k]
                   - f_7 * ksh1_390[k]
                   + f_3 * pc_x[k] * ksi_516[k];

        t_661[k] = f_14 * isi_373[k]
                   + f_3 * pc_y[k] * ksi_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, isi_346, isi_518, isi_519, ksh0_392, \
                         ksh0_393, ksh1_392, ksh1_393, ksi_514, ksi_518, \
                         ksi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_14 * isi_518[k]
                   + f_6 * ksh0_392[k]
                   - f_7 * ksh1_392[k]
                   + f_3 * pc_x[k] * ksi_518[k];

        t_663[k] = f_14 * isi_519[k]
                   + f_4 * ksh0_393[k]
                   - f_5 * ksh1_393[k]
                   + f_3 * pc_x[k] * ksi_519[k];

        t_664[k] = f_15 * isi_346[k]
                   + f_3 * pc_z[k] * ksi_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, isi_378, isi_521, isi_522, ksh0_395, \
                         ksh0_396, ksh1_395, ksh1_396, ksi_518, ksi_521, \
                         ksi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_14 * isi_521[k]
                   + f_4 * ksh0_395[k]
                   - f_5 * ksh1_395[k]
                   + f_3 * pc_x[k] * ksi_521[k];

        t_666[k] = f_14 * isi_522[k]
                   + f_4 * ksh0_396[k]
                   - f_5 * ksh1_396[k]
                   + f_3 * pc_x[k] * ksi_522[k];

        t_667[k] = f_14 * isi_378[k]
                   + f_3 * pc_y[k] * ksi_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, isi_524, isi_525, isi_526, isi_527, \
                         ksh0_398, ksh1_398, ksi_524, ksi_525, ksi_526, \
                         ksi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_14 * isi_524[k]
                   + f_4 * ksh0_398[k]
                   - f_5 * ksh1_398[k]
                   + f_3 * pc_x[k] * ksi_524[k];

        t_669[k] = f_14 * isi_525[k]
                   + f_3 * pc_x[k] * ksi_525[k];

        t_670[k] = f_14 * isi_526[k]
                   + f_3 * pc_x[k] * ksi_526[k];

        t_671[k] = f_14 * isi_527[k]
                   + f_3 * pc_x[k] * ksi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, isi_528, isi_529, isi_530, isi_531, \
                         ksi_528, ksi_529, ksi_530, ksi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_14 * isi_528[k]
                   + f_3 * pc_x[k] * ksi_528[k];

        t_673[k] = f_14 * isi_529[k]
                   + f_3 * pc_x[k] * ksi_529[k];

        t_674[k] = f_14 * isi_530[k]
                   + f_3 * pc_x[k] * ksi_530[k];

        t_675[k] = f_14 * isi_531[k]
                   + f_3 * pc_x[k] * ksi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, isi_357, isi_385, isi_387, ksh0_393, \
                         ksh0_395, ksh1_393, ksh1_395, ksi_525, \
                         ksi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * isi_385[k]
                   + f_1 * ksh0_393[k]
                   - f_2 * ksh1_393[k]
                   + f_3 * pc_y[k] * ksi_525[k];

        t_677[k] = f_15 * isi_357[k]
                   + f_3 * pc_z[k] * ksi_525[k];

        t_678[k] = f_14 * isi_387[k]
                   + f_10 * ksh0_395[k]
                   - f_11 * ksh1_395[k]
                   + f_3 * pc_y[k] * ksi_527[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, isi_388, isi_389, isi_390, ksh0_396, \
                         ksh0_397, ksh0_398, ksh1_396, ksh1_397, ksh1_398, ksi_528, ksi_529, \
                         ksi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * isi_388[k]
                   + f_8 * ksh0_396[k]
                   - f_9 * ksh1_396[k]
                   + f_3 * pc_y[k] * ksi_528[k];

        t_680[k] = f_14 * isi_389[k]
                   + f_6 * ksh0_397[k]
                   - f_7 * ksh1_397[k]
                   + f_3 * pc_y[k] * ksi_529[k];

        t_681[k] = f_14 * isi_390[k]
                   + f_4 * ksh0_398[k]
                   - f_5 * ksh1_398[k]
                   + f_3 * pc_y[k] * ksi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_y, pc_y, pc_z, isk0_504, isi_363, \
                         isi_391, isi_392, isk1_504, ksh0_398, ksh1_398, ksi_531, \
                         ksi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * isi_391[k]
                   + f_3 * pc_y[k] * ksi_531[k];

        t_683[k] = f_15 * isi_363[k]
                   + f_1 * ksh0_398[k]
                   - f_2 * ksh1_398[k]
                   + f_3 * pc_z[k] * ksi_531[k];

        t_684[k] = pa_y[k] * isk0_504[k]
                   - f_12 * pc_y[k] * isk1_504[k];

        t_685[k] = f_13 * isi_392[k]
                   + f_3 * pc_y[k] * ksi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_y, pc_y, pc_z, isk0_507, isk0_509, \
                         isi_364, isi_393, isi_394, isk1_507, isk1_509, ksi_532, \
                         ksi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * isi_364[k]
                   + f_3 * pc_z[k] * ksi_532[k];

        t_687[k] = pa_y[k] * isk0_507[k]
                   + f_14 * isi_393[k]
                   - f_12 * pc_y[k] * isk1_507[k];

        t_688[k] = f_13 * isi_394[k]
                   + f_3 * pc_y[k] * ksi_534[k];

        t_689[k] = pa_y[k] * isk0_509[k]
                   - f_12 * pc_y[k] * isk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_y, pc_y, pc_z, isk0_510, isk0_513, \
                         isi_367, isi_395, isi_397, isk1_510, isk1_513, ksi_535, \
                         ksi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_y[k] * isk0_510[k]
                   + f_15 * isi_395[k]
                   - f_12 * pc_y[k] * isk1_510[k];

        t_691[k] = f_16 * isi_367[k]
                   + f_3 * pc_z[k] * ksi_535[k];

        t_692[k] = f_13 * isi_397[k]
                   + f_3 * pc_y[k] * ksi_537[k];

        t_693[k] = pa_y[k] * isk0_513[k]
                   - f_12 * pc_y[k] * isk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pa_y, pc_y, pc_z, isk0_514, isk0_516, isi_370, \
                         isi_398, isi_400, isk1_514, isk1_516, \
                         ksi_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pa_y[k] * isk0_514[k]
                   + f_16 * isi_398[k]
                   - f_12 * pc_y[k] * isk1_514[k];

        t_695[k] = f_16 * isi_370[k]
                   + f_3 * pc_z[k] * ksi_538[k];

        t_696[k] = pa_y[k] * isk0_516[k]
                   + f_14 * isi_400[k]
                   - f_12 * pc_y[k] * isk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pa_y, pc_y, pc_z, isk0_518, isk0_519, \
                         isi_374, isi_401, isi_402, isk1_518, isk1_519, ksi_541, \
                         ksi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * isi_401[k]
                   + f_3 * pc_y[k] * ksi_541[k];

        t_698[k] = pa_y[k] * isk0_518[k]
                   - f_12 * pc_y[k] * isk1_518[k];

        t_699[k] = pa_y[k] * isk0_519[k]
                   + f_17 * isi_402[k]
                   - f_12 * pc_y[k] * isk1_519[k];

        t_700[k] = f_16 * isi_374[k]
                   + f_3 * pc_z[k] * ksi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pa_y, pc_y, isk0_521, isk0_522, isk0_524, \
                         isi_404, isi_405, isi_406, isk1_521, isk1_522, isk1_524, \
                         ksi_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pa_y[k] * isk0_521[k]
                   + f_15 * isi_404[k]
                   - f_12 * pc_y[k] * isk1_521[k];

        t_702[k] = pa_y[k] * isk0_522[k]
                   + f_14 * isi_405[k]
                   - f_12 * pc_y[k] * isk1_522[k];

        t_703[k] = f_13 * isi_406[k]
                   + f_3 * pc_y[k] * ksi_546[k];

        t_704[k] = pa_y[k] * isk0_524[k]
                   - f_12 * pc_y[k] * isk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, isi_553, isi_554, isi_555, \
                         isi_556, isi_557, ksi_553, ksi_554, ksi_555, ksi_556, \
                         ksi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_14 * isi_553[k]
                   + f_3 * pc_x[k] * ksi_553[k];

        t_706[k] = f_14 * isi_554[k]
                   + f_3 * pc_x[k] * ksi_554[k];

        t_707[k] = f_14 * isi_555[k]
                   + f_3 * pc_x[k] * ksi_555[k];

        t_708[k] = f_14 * isi_556[k]
                   + f_3 * pc_x[k] * ksi_556[k];

        t_709[k] = f_14 * isi_557[k]
                   + f_3 * pc_x[k] * ksi_557[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isk0_539 = buffer.data(isk0 + 539);
    const auto *isk0_540 = buffer.data(isk0 + 540);
    const auto *isk0_543 = buffer.data(isk0 + 543);
    const auto *isk0_546 = buffer.data(isk0 + 546);
    const auto *isk0_550 = buffer.data(isk0 + 550);
    const auto *isk0_555 = buffer.data(isk0 + 555);
    const auto *isk0_756 = buffer.data(isk0 + 756);
    const auto *isk0_759 = buffer.data(isk0 + 759);
    const auto *isk0_762 = buffer.data(isk0 + 762);
    const auto *isk0_766 = buffer.data(isk0 + 766);
    const auto *isk0_771 = buffer.data(isk0 + 771);
    const auto *isk0_784 = buffer.data(isk0 + 784);
    const auto *isk0_786 = buffer.data(isk0 + 786);
    const auto *isk0_787 = buffer.data(isk0 + 787);
    const auto *isk0_788 = buffer.data(isk0 + 788);
    const auto *isk0_789 = buffer.data(isk0 + 789);
    const auto *isk0_791 = buffer.data(isk0 + 791);
    const auto *isk0_797 = buffer.data(isk0 + 797);
    const auto *isk0_801 = buffer.data(isk0 + 801);
    const auto *isk0_804 = buffer.data(isk0 + 804);
    const auto *isk0_806 = buffer.data(isk0 + 806);
    const auto *isk0_809 = buffer.data(isk0 + 809);
    const auto *isk0_810 = buffer.data(isk0 + 810);
    const auto *isk0_812 = buffer.data(isk0 + 812);
    const auto *isk0_820 = buffer.data(isk0 + 820);
    const auto *isk0_822 = buffer.data(isk0 + 822);
    const auto *isk0_823 = buffer.data(isk0 + 823);
    const auto *isk0_824 = buffer.data(isk0 + 824);
    const auto *isk0_825 = buffer.data(isk0 + 825);
    const auto *isk0_827 = buffer.data(isk0 + 827);
    const auto *isk0_828 = buffer.data(isk0 + 828);

    const auto *isi_385 = buffer.data(isi + 385);
    const auto *isi_392 = buffer.data(isi + 392);
    const auto *isi_413 = buffer.data(isi + 413);
    const auto *isi_415 = buffer.data(isi + 415);
    const auto *isi_416 = buffer.data(isi + 416);
    const auto *isi_417 = buffer.data(isi + 417);
    const auto *isi_418 = buffer.data(isi + 418);
    const auto *isi_419 = buffer.data(isi + 419);
    const auto *isi_420 = buffer.data(isi + 420);
    const auto *isi_423 = buffer.data(isi + 423);
    const auto *isi_425 = buffer.data(isi + 425);
    const auto *isi_426 = buffer.data(isi + 426);
    const auto *isi_429 = buffer.data(isi + 429);
    const auto *isi_430 = buffer.data(isi + 430);
    const auto *isi_434 = buffer.data(isi + 434);
    const auto *isi_441 = buffer.data(isi + 441);
    const auto *isi_447 = buffer.data(isi + 447);
    const auto *isi_448 = buffer.data(isi + 448);
    const auto *isi_450 = buffer.data(isi + 450);
    const auto *isi_453 = buffer.data(isi + 453);
    const auto *isi_457 = buffer.data(isi + 457);
    const auto *isi_462 = buffer.data(isi + 462);
    const auto *isi_475 = buffer.data(isi + 475);
    const auto *isi_476 = buffer.data(isi + 476);
    const auto *isi_558 = buffer.data(isi + 558);
    const auto *isi_559 = buffer.data(isi + 559);
    const auto *isi_560 = buffer.data(isi + 560);
    const auto *isi_565 = buffer.data(isi + 565);
    const auto *isi_569 = buffer.data(isi + 569);
    const auto *isi_574 = buffer.data(isi + 574);
    const auto *isi_580 = buffer.data(isi + 580);
    const auto *isi_581 = buffer.data(isi + 581);
    const auto *isi_582 = buffer.data(isi + 582);
    const auto *isi_583 = buffer.data(isi + 583);
    const auto *isi_584 = buffer.data(isi + 584);
    const auto *isi_585 = buffer.data(isi + 585);
    const auto *isi_587 = buffer.data(isi + 587);
    const auto *isi_588 = buffer.data(isi + 588);
    const auto *isi_591 = buffer.data(isi + 591);
    const auto *isi_594 = buffer.data(isi + 594);
    const auto *isi_598 = buffer.data(isi + 598);
    const auto *isi_603 = buffer.data(isi + 603);
    const auto *isi_609 = buffer.data(isi + 609);
    const auto *isi_611 = buffer.data(isi + 611);
    const auto *isi_612 = buffer.data(isi + 612);
    const auto *isi_613 = buffer.data(isi + 613);
    const auto *isi_614 = buffer.data(isi + 614);
    const auto *isi_615 = buffer.data(isi + 615);
    const auto *isi_621 = buffer.data(isi + 621);
    const auto *isi_625 = buffer.data(isi + 625);
    const auto *isi_628 = buffer.data(isi + 628);
    const auto *isi_630 = buffer.data(isi + 630);
    const auto *isi_633 = buffer.data(isi + 633);
    const auto *isi_634 = buffer.data(isi + 634);
    const auto *isi_636 = buffer.data(isi + 636);
    const auto *isi_637 = buffer.data(isi + 637);
    const auto *isi_638 = buffer.data(isi + 638);
    const auto *isi_639 = buffer.data(isi + 639);
    const auto *isi_640 = buffer.data(isi + 640);
    const auto *isi_641 = buffer.data(isi + 641);
    const auto *isi_642 = buffer.data(isi + 642);
    const auto *isi_643 = buffer.data(isi + 643);
    const auto *isi_644 = buffer.data(isi + 644);

    const auto *isk1_539 = buffer.data(isk1 + 539);
    const auto *isk1_540 = buffer.data(isk1 + 540);
    const auto *isk1_543 = buffer.data(isk1 + 543);
    const auto *isk1_546 = buffer.data(isk1 + 546);
    const auto *isk1_550 = buffer.data(isk1 + 550);
    const auto *isk1_555 = buffer.data(isk1 + 555);
    const auto *isk1_756 = buffer.data(isk1 + 756);
    const auto *isk1_759 = buffer.data(isk1 + 759);
    const auto *isk1_762 = buffer.data(isk1 + 762);
    const auto *isk1_766 = buffer.data(isk1 + 766);
    const auto *isk1_771 = buffer.data(isk1 + 771);
    const auto *isk1_784 = buffer.data(isk1 + 784);
    const auto *isk1_786 = buffer.data(isk1 + 786);
    const auto *isk1_787 = buffer.data(isk1 + 787);
    const auto *isk1_788 = buffer.data(isk1 + 788);
    const auto *isk1_789 = buffer.data(isk1 + 789);
    const auto *isk1_791 = buffer.data(isk1 + 791);
    const auto *isk1_797 = buffer.data(isk1 + 797);
    const auto *isk1_801 = buffer.data(isk1 + 801);
    const auto *isk1_804 = buffer.data(isk1 + 804);
    const auto *isk1_806 = buffer.data(isk1 + 806);
    const auto *isk1_809 = buffer.data(isk1 + 809);
    const auto *isk1_810 = buffer.data(isk1 + 810);
    const auto *isk1_812 = buffer.data(isk1 + 812);
    const auto *isk1_820 = buffer.data(isk1 + 820);
    const auto *isk1_822 = buffer.data(isk1 + 822);
    const auto *isk1_823 = buffer.data(isk1 + 823);
    const auto *isk1_824 = buffer.data(isk1 + 824);
    const auto *isk1_825 = buffer.data(isk1 + 825);
    const auto *isk1_827 = buffer.data(isk1 + 827);
    const auto *isk1_828 = buffer.data(isk1 + 828);

    const auto *ksh0_414 = buffer.data(ksh0 + 414);
    const auto *ksh0_416 = buffer.data(ksh0 + 416);
    const auto *ksh0_417 = buffer.data(ksh0 + 417);
    const auto *ksh0_418 = buffer.data(ksh0 + 418);
    const auto *ksh0_419 = buffer.data(ksh0 + 419);
    const auto *ksh0_420 = buffer.data(ksh0 + 420);
    const auto *ksh0_421 = buffer.data(ksh0 + 421);
    const auto *ksh0_422 = buffer.data(ksh0 + 422);
    const auto *ksh0_423 = buffer.data(ksh0 + 423);
    const auto *ksh0_424 = buffer.data(ksh0 + 424);
    const auto *ksh0_425 = buffer.data(ksh0 + 425);
    const auto *ksh0_426 = buffer.data(ksh0 + 426);
    const auto *ksh0_427 = buffer.data(ksh0 + 427);
    const auto *ksh0_428 = buffer.data(ksh0 + 428);
    const auto *ksh0_429 = buffer.data(ksh0 + 429);
    const auto *ksh0_434 = buffer.data(ksh0 + 434);
    const auto *ksh0_435 = buffer.data(ksh0 + 435);
    const auto *ksh0_436 = buffer.data(ksh0 + 436);
    const auto *ksh0_437 = buffer.data(ksh0 + 437);
    const auto *ksh0_438 = buffer.data(ksh0 + 438);
    const auto *ksh0_439 = buffer.data(ksh0 + 439);
    const auto *ksh0_440 = buffer.data(ksh0 + 440);
    const auto *ksh0_441 = buffer.data(ksh0 + 441);
    const auto *ksh0_443 = buffer.data(ksh0 + 443);
    const auto *ksh0_444 = buffer.data(ksh0 + 444);
    const auto *ksh0_446 = buffer.data(ksh0 + 446);
    const auto *ksh0_447 = buffer.data(ksh0 + 447);
    const auto *ksh0_448 = buffer.data(ksh0 + 448);
    const auto *ksh0_450 = buffer.data(ksh0 + 450);

    const auto *ksh1_414 = buffer.data(ksh1 + 414);
    const auto *ksh1_416 = buffer.data(ksh1 + 416);
    const auto *ksh1_417 = buffer.data(ksh1 + 417);
    const auto *ksh1_418 = buffer.data(ksh1 + 418);
    const auto *ksh1_419 = buffer.data(ksh1 + 419);
    const auto *ksh1_420 = buffer.data(ksh1 + 420);
    const auto *ksh1_421 = buffer.data(ksh1 + 421);
    const auto *ksh1_422 = buffer.data(ksh1 + 422);
    const auto *ksh1_423 = buffer.data(ksh1 + 423);
    const auto *ksh1_424 = buffer.data(ksh1 + 424);
    const auto *ksh1_425 = buffer.data(ksh1 + 425);
    const auto *ksh1_426 = buffer.data(ksh1 + 426);
    const auto *ksh1_427 = buffer.data(ksh1 + 427);
    const auto *ksh1_428 = buffer.data(ksh1 + 428);
    const auto *ksh1_429 = buffer.data(ksh1 + 429);
    const auto *ksh1_434 = buffer.data(ksh1 + 434);
    const auto *ksh1_435 = buffer.data(ksh1 + 435);
    const auto *ksh1_436 = buffer.data(ksh1 + 436);
    const auto *ksh1_437 = buffer.data(ksh1 + 437);
    const auto *ksh1_438 = buffer.data(ksh1 + 438);
    const auto *ksh1_439 = buffer.data(ksh1 + 439);
    const auto *ksh1_440 = buffer.data(ksh1 + 440);
    const auto *ksh1_441 = buffer.data(ksh1 + 441);
    const auto *ksh1_443 = buffer.data(ksh1 + 443);
    const auto *ksh1_444 = buffer.data(ksh1 + 444);
    const auto *ksh1_446 = buffer.data(ksh1 + 446);
    const auto *ksh1_447 = buffer.data(ksh1 + 447);
    const auto *ksh1_448 = buffer.data(ksh1 + 448);
    const auto *ksh1_450 = buffer.data(ksh1 + 450);

    const auto *ksi_553 = buffer.data(ksi + 553);
    const auto *ksi_555 = buffer.data(ksi + 555);
    const auto *ksi_556 = buffer.data(ksi + 556);
    const auto *ksi_557 = buffer.data(ksi + 557);
    const auto *ksi_558 = buffer.data(ksi + 558);
    const auto *ksi_559 = buffer.data(ksi + 559);
    const auto *ksi_560 = buffer.data(ksi + 560);
    const auto *ksi_561 = buffer.data(ksi + 561);
    const auto *ksi_562 = buffer.data(ksi + 562);
    const auto *ksi_563 = buffer.data(ksi + 563);
    const auto *ksi_564 = buffer.data(ksi + 564);
    const auto *ksi_565 = buffer.data(ksi + 565);
    const auto *ksi_566 = buffer.data(ksi + 566);
    const auto *ksi_567 = buffer.data(ksi + 567);
    const auto *ksi_568 = buffer.data(ksi + 568);
    const auto *ksi_569 = buffer.data(ksi + 569);
    const auto *ksi_570 = buffer.data(ksi + 570);
    const auto *ksi_571 = buffer.data(ksi + 571);
    const auto *ksi_572 = buffer.data(ksi + 572);
    const auto *ksi_573 = buffer.data(ksi + 573);
    const auto *ksi_574 = buffer.data(ksi + 574);
    const auto *ksi_580 = buffer.data(ksi + 580);
    const auto *ksi_581 = buffer.data(ksi + 581);
    const auto *ksi_582 = buffer.data(ksi + 582);
    const auto *ksi_583 = buffer.data(ksi + 583);
    const auto *ksi_584 = buffer.data(ksi + 584);
    const auto *ksi_585 = buffer.data(ksi + 585);
    const auto *ksi_586 = buffer.data(ksi + 586);
    const auto *ksi_587 = buffer.data(ksi + 587);
    const auto *ksi_588 = buffer.data(ksi + 588);
    const auto *ksi_589 = buffer.data(ksi + 589);
    const auto *ksi_590 = buffer.data(ksi + 590);
    const auto *ksi_591 = buffer.data(ksi + 591);
    const auto *ksi_593 = buffer.data(ksi + 593);
    const auto *ksi_594 = buffer.data(ksi + 594);
    const auto *ksi_595 = buffer.data(ksi + 595);
    const auto *ksi_597 = buffer.data(ksi + 597);
    const auto *ksi_598 = buffer.data(ksi + 598);
    const auto *ksi_599 = buffer.data(ksi + 599);
    const auto *ksi_600 = buffer.data(ksi + 600);
    const auto *ksi_602 = buffer.data(ksi + 602);
    const auto *ksi_603 = buffer.data(ksi + 603);
    const auto *ksi_609 = buffer.data(ksi + 609);
    const auto *ksi_611 = buffer.data(ksi + 611);
    const auto *ksi_612 = buffer.data(ksi + 612);
    const auto *ksi_613 = buffer.data(ksi + 613);
    const auto *ksi_614 = buffer.data(ksi + 614);
    const auto *ksi_615 = buffer.data(ksi + 615);
    const auto *ksi_616 = buffer.data(ksi + 616);
    const auto *ksi_618 = buffer.data(ksi + 618);
    const auto *ksi_619 = buffer.data(ksi + 619);
    const auto *ksi_621 = buffer.data(ksi + 621);
    const auto *ksi_622 = buffer.data(ksi + 622);
    const auto *ksi_625 = buffer.data(ksi + 625);
    const auto *ksi_626 = buffer.data(ksi + 626);
    const auto *ksi_630 = buffer.data(ksi + 630);
    const auto *ksi_637 = buffer.data(ksi + 637);
    const auto *ksi_638 = buffer.data(ksi + 638);
    const auto *ksi_639 = buffer.data(ksi + 639);
    const auto *ksi_640 = buffer.data(ksi + 640);
    const auto *ksi_641 = buffer.data(ksi + 641);
    const auto *ksi_642 = buffer.data(ksi + 642);
    const auto *ksi_643 = buffer.data(ksi + 643);
    const auto *ksi_644 = buffer.data(ksi + 644);

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, isi_385, isi_413, \
                         isi_558, isi_559, ksh0_414, ksh1_414, ksi_553, ksi_558, \
                         ksi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_14 * isi_558[k]
                   + f_3 * pc_x[k] * ksi_558[k];

        t_711[k] = f_14 * isi_559[k]
                   + f_3 * pc_x[k] * ksi_559[k];

        t_712[k] = f_13 * isi_413[k]
                   + f_1 * ksh0_414[k]
                   - f_2 * ksh1_414[k]
                   + f_3 * pc_y[k] * ksi_553[k];

        t_713[k] = f_16 * isi_385[k]
                   + f_3 * pc_z[k] * ksi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, isi_415, isi_416, isi_417, ksh0_416, \
                         ksh0_417, ksh0_418, ksh1_416, ksh1_417, ksh1_418, ksi_555, ksi_556, \
                         ksi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * isi_415[k]
                   + f_10 * ksh0_416[k]
                   - f_11 * ksh1_416[k]
                   + f_3 * pc_y[k] * ksi_555[k];

        t_715[k] = f_13 * isi_416[k]
                   + f_8 * ksh0_417[k]
                   - f_9 * ksh1_417[k]
                   + f_3 * pc_y[k] * ksi_556[k];

        t_716[k] = f_13 * isi_417[k]
                   + f_6 * ksh0_418[k]
                   - f_7 * ksh1_418[k]
                   + f_3 * pc_y[k] * ksi_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_y, pc_y, isk0_539, isi_418, isi_419, \
                         isk1_539, ksh0_419, ksh1_419, ksi_558, \
                         ksi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * isi_418[k]
                   + f_4 * ksh0_419[k]
                   - f_5 * ksh1_419[k]
                   + f_3 * pc_y[k] * ksi_558[k];

        t_718[k] = f_13 * isi_419[k]
                   + f_3 * pc_y[k] * ksi_559[k];

        t_719[k] = pa_y[k] * isk0_539[k]
                   - f_12 * pc_y[k] * isk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pc_x, pc_y, pc_z, isi_392, \
                         isi_560, ksh0_420, ksh1_420, ksi_560, ksi_561, \
                         ksi_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_14 * isi_560[k]
                   + f_1 * ksh0_420[k]
                   - f_2 * ksh1_420[k]
                   + f_3 * pc_x[k] * ksi_560[k];

        t_721[k] = f_3 * pc_y[k] * ksi_560[k];

        t_722[k] = f_17 * isi_392[k]
                   + f_3 * pc_z[k] * ksi_560[k];

        t_723[k] = f_4 * ksh0_420[k]
                   - f_5 * ksh1_420[k]
                   + f_3 * pc_y[k] * ksi_561[k];

        t_724[k] = f_3 * pc_y[k] * ksi_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, pc_y, isi_565, ksh0_421, ksh0_422, \
                         ksh0_425, ksh1_421, ksh1_422, ksh1_425, ksi_563, ksi_564, \
                         ksi_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_14 * isi_565[k]
                   + f_10 * ksh0_425[k]
                   - f_11 * ksh1_425[k]
                   + f_3 * pc_x[k] * ksi_565[k];

        t_726[k] = f_6 * ksh0_421[k]
                   - f_7 * ksh1_421[k]
                   + f_3 * pc_y[k] * ksi_563[k];

        t_727[k] = f_4 * ksh0_422[k]
                   - f_5 * ksh1_422[k]
                   + f_3 * pc_y[k] * ksi_564[k];

        t_728[k] = f_3 * pc_y[k] * ksi_565[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, isi_569, ksh0_423, ksh0_424, \
                         ksh0_429, ksh1_423, ksh1_424, ksh1_429, ksi_566, ksi_567, \
                         ksi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_14 * isi_569[k]
                   + f_8 * ksh0_429[k]
                   - f_9 * ksh1_429[k]
                   + f_3 * pc_x[k] * ksi_569[k];

        t_730[k] = f_8 * ksh0_423[k]
                   - f_9 * ksh1_423[k]
                   + f_3 * pc_y[k] * ksi_566[k];

        t_731[k] = f_6 * ksh0_424[k]
                   - f_7 * ksh1_424[k]
                   + f_3 * pc_y[k] * ksi_567[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_x, pc_y, isi_574, ksh0_425, ksh0_434, \
                         ksh1_425, ksh1_434, ksi_568, ksi_569, \
                         ksi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_4 * ksh0_425[k]
                   - f_5 * ksh1_425[k]
                   + f_3 * pc_y[k] * ksi_568[k];

        t_733[k] = f_3 * pc_y[k] * ksi_569[k];

        t_734[k] = f_14 * isi_574[k]
                   + f_6 * ksh0_434[k]
                   - f_7 * ksh1_434[k]
                   + f_3 * pc_x[k] * ksi_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pc_y, ksh0_426, ksh0_427, ksh0_428, ksh1_426, \
                         ksh1_427, ksh1_428, ksi_570, ksi_571, \
                         ksi_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_10 * ksh0_426[k]
                   - f_11 * ksh1_426[k]
                   + f_3 * pc_y[k] * ksi_570[k];

        t_736[k] = f_8 * ksh0_427[k]
                   - f_9 * ksh1_427[k]
                   + f_3 * pc_y[k] * ksi_571[k];

        t_737[k] = f_6 * ksh0_428[k]
                   - f_7 * ksh1_428[k]
                   + f_3 * pc_y[k] * ksi_572[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, pc_x, pc_y, isi_580, isi_581, ksh0_429, \
                         ksh0_440, ksh1_429, ksh1_440, ksi_573, ksi_574, ksi_580, \
                         ksi_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_4 * ksh0_429[k]
                   - f_5 * ksh1_429[k]
                   + f_3 * pc_y[k] * ksi_573[k];

        t_739[k] = f_3 * pc_y[k] * ksi_574[k];

        t_740[k] = f_14 * isi_580[k]
                   + f_4 * ksh0_440[k]
                   - f_5 * ksh1_440[k]
                   + f_3 * pc_x[k] * ksi_580[k];

        t_741[k] = f_14 * isi_581[k]
                   + f_3 * pc_x[k] * ksi_581[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, pc_x, pc_y, isi_582, isi_583, \
                         isi_584, isi_585, ksi_580, ksi_582, ksi_583, ksi_584, \
                         ksi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_14 * isi_582[k]
                   + f_3 * pc_x[k] * ksi_582[k];

        t_743[k] = f_14 * isi_583[k]
                   + f_3 * pc_x[k] * ksi_583[k];

        t_744[k] = f_14 * isi_584[k]
                   + f_3 * pc_x[k] * ksi_584[k];

        t_745[k] = f_14 * isi_585[k]
                   + f_3 * pc_x[k] * ksi_585[k];

        t_746[k] = f_3 * pc_y[k] * ksi_580[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_x, pc_y, isi_587, ksh0_435, ksh0_436, \
                         ksh1_435, ksh1_436, ksi_581, ksi_582, \
                         ksi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_14 * isi_587[k]
                   + f_3 * pc_x[k] * ksi_587[k];

        t_748[k] = f_1 * ksh0_435[k]
                   - f_2 * ksh1_435[k]
                   + f_3 * pc_y[k] * ksi_581[k];

        t_749[k] = f_19 * ksh0_436[k]
                   - f_20 * ksh1_436[k]
                   + f_3 * pc_y[k] * ksi_582[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, ksh0_437, ksh0_438, ksh0_439, ksh1_437, \
                         ksh1_438, ksh1_439, ksi_583, ksi_584, \
                         ksi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_10 * ksh0_437[k]
                   - f_11 * ksh1_437[k]
                   + f_3 * pc_y[k] * ksi_583[k];

        t_751[k] = f_8 * ksh0_438[k]
                   - f_9 * ksh1_438[k]
                   + f_3 * pc_y[k] * ksi_584[k];

        t_752[k] = f_6 * ksh0_439[k]
                   - f_7 * ksh1_439[k]
                   + f_3 * pc_y[k] * ksi_585[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pa_x, pc_x, pc_y, pc_z, isk0_756, \
                         isi_419, isi_588, isk1_756, ksh0_440, ksh1_440, ksi_586, \
                         ksi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * ksh0_440[k]
                   - f_5 * ksh1_440[k]
                   + f_3 * pc_y[k] * ksi_586[k];

        t_754[k] = f_3 * pc_y[k] * ksi_587[k];

        t_755[k] = f_17 * isi_419[k]
                   + f_1 * ksh0_440[k]
                   - f_2 * ksh1_440[k]
                   + f_3 * pc_z[k] * ksi_587[k];

        t_756[k] = pa_x[k] * isk0_756[k]
                   + f_0 * isi_588[k]
                   - f_12 * pc_x[k] * isk1_756[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pa_x, pc_x, pc_y, pc_z, isk0_759, \
                         isi_420, isi_591, isk1_759, ksi_588, ksi_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_18 * isi_420[k]
                   + f_3 * pc_y[k] * ksi_588[k];

        t_758[k] = f_3 * pc_z[k] * ksi_588[k];

        t_759[k] = pa_x[k] * isk0_759[k]
                   + f_17 * isi_591[k]
                   - f_12 * pc_x[k] * isk1_759[k];

        t_760[k] = f_3 * pc_z[k] * ksi_589[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pa_x, pc_x, pc_z, isk0_762, isi_594, isk1_762, \
                         ksh0_441, ksh1_441, ksi_590, ksi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_4 * ksh0_441[k]
                   - f_5 * ksh1_441[k]
                   + f_3 * pc_z[k] * ksi_590[k];

        t_762[k] = pa_x[k] * isk0_762[k]
                   + f_16 * isi_594[k]
                   - f_12 * pc_x[k] * isk1_762[k];

        t_763[k] = f_3 * pc_z[k] * ksi_591[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pa_x, pc_x, pc_y, pc_z, isk0_766, \
                         isi_425, isi_598, isk1_766, ksh0_443, ksh1_443, ksi_593, \
                         ksi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_18 * isi_425[k]
                   + f_3 * pc_y[k] * ksi_593[k];

        t_765[k] = f_6 * ksh0_443[k]
                   - f_7 * ksh1_443[k]
                   + f_3 * pc_z[k] * ksi_593[k];

        t_766[k] = pa_x[k] * isk0_766[k]
                   + f_15 * isi_598[k]
                   - f_12 * pc_x[k] * isk1_766[k];

        t_767[k] = f_3 * pc_z[k] * ksi_594[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_y, pc_z, isi_429, ksh0_444, ksh0_446, \
                         ksh1_444, ksh1_446, ksi_595, ksi_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * ksh0_444[k]
                   - f_5 * ksh1_444[k]
                   + f_3 * pc_z[k] * ksi_595[k];

        t_769[k] = f_18 * isi_429[k]
                   + f_3 * pc_y[k] * ksi_597[k];

        t_770[k] = f_8 * ksh0_446[k]
                   - f_9 * ksh1_446[k]
                   + f_3 * pc_z[k] * ksi_597[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, pa_x, pc_x, pc_z, isk0_771, isi_603, isk1_771, \
                         ksh0_447, ksh1_447, ksi_598, ksi_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = pa_x[k] * isk0_771[k]
                   + f_14 * isi_603[k]
                   - f_12 * pc_x[k] * isk1_771[k];

        t_772[k] = f_3 * pc_z[k] * ksi_598[k];

        t_773[k] = f_4 * ksh0_447[k]
                   - f_5 * ksh1_447[k]
                   + f_3 * pc_z[k] * ksi_599[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pc_x, pc_y, pc_z, isi_434, isi_609, \
                         ksh0_448, ksh0_450, ksh1_448, ksh1_450, ksi_600, ksi_602, \
                         ksi_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_6 * ksh0_448[k]
                   - f_7 * ksh1_448[k]
                   + f_3 * pc_z[k] * ksi_600[k];

        t_775[k] = f_18 * isi_434[k]
                   + f_3 * pc_y[k] * ksi_602[k];

        t_776[k] = f_10 * ksh0_450[k]
                   - f_11 * ksh1_450[k]
                   + f_3 * pc_z[k] * ksi_602[k];

        t_777[k] = f_13 * isi_609[k]
                   + f_3 * pc_x[k] * ksi_609[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, pc_x, pc_z, isi_611, isi_612, \
                         isi_613, isi_614, ksi_603, ksi_611, ksi_612, ksi_613, \
                         ksi_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_3 * pc_z[k] * ksi_603[k];

        t_779[k] = f_13 * isi_611[k]
                   + f_3 * pc_x[k] * ksi_611[k];

        t_780[k] = f_13 * isi_612[k]
                   + f_3 * pc_x[k] * ksi_612[k];

        t_781[k] = f_13 * isi_613[k]
                   + f_3 * pc_x[k] * ksi_613[k];

        t_782[k] = f_13 * isi_614[k]
                   + f_3 * pc_x[k] * ksi_614[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pa_x, pc_x, pc_z, isk0_784, isk0_786, \
                         isi_615, isk1_784, isk1_786, ksi_609, \
                         ksi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_13 * isi_615[k]
                   + f_3 * pc_x[k] * ksi_615[k];

        t_784[k] = pa_x[k] * isk0_784[k]
                   - f_12 * pc_x[k] * isk1_784[k];

        t_785[k] = f_3 * pc_z[k] * ksi_609[k];

        t_786[k] = pa_x[k] * isk0_786[k]
                   - f_12 * pc_x[k] * isk1_786[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, t_790, pa_x, pc_x, pc_y, isk0_787, isk0_788, \
                         isk0_789, isi_447, isk1_787, isk1_788, isk1_789, \
                         ksi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = pa_x[k] * isk0_787[k]
                   - f_12 * pc_x[k] * isk1_787[k];

        t_788[k] = pa_x[k] * isk0_788[k]
                   - f_12 * pc_x[k] * isk1_788[k];

        t_789[k] = pa_x[k] * isk0_789[k]
                   - f_12 * pc_x[k] * isk1_789[k];

        t_790[k] = f_18 * isi_447[k]
                   + f_3 * pc_y[k] * ksi_615[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pa_x, pa_z, pc_x, pc_y, pc_z, isk0_540, \
                         isk0_791, isi_420, isi_448, isk1_540, isk1_791, \
                         ksi_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = pa_x[k] * isk0_791[k]
                   - f_12 * pc_x[k] * isk1_791[k];

        t_792[k] = pa_z[k] * isk0_540[k]
                   - f_12 * pc_z[k] * isk1_540[k];

        t_793[k] = f_17 * isi_448[k]
                   + f_3 * pc_y[k] * ksi_616[k];

        t_794[k] = f_13 * isi_420[k]
                   + f_3 * pc_z[k] * ksi_616[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, pa_x, pa_z, pc_x, pc_y, pc_z, isk0_543, \
                         isk0_797, isi_450, isi_621, isk1_543, isk1_797, \
                         ksi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = pa_z[k] * isk0_543[k]
                   - f_12 * pc_z[k] * isk1_543[k];

        t_796[k] = f_17 * isi_450[k]
                   + f_3 * pc_y[k] * ksi_618[k];

        t_797[k] = pa_x[k] * isk0_797[k]
                   + f_17 * isi_621[k]
                   - f_12 * pc_x[k] * isk1_797[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, pa_z, pc_y, pc_z, isk0_546, isi_423, isi_453, \
                         isk1_546, ksi_619, ksi_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = pa_z[k] * isk0_546[k]
                   - f_12 * pc_z[k] * isk1_546[k];

        t_799[k] = f_13 * isi_423[k]
                   + f_3 * pc_z[k] * ksi_619[k];

        t_800[k] = f_17 * isi_453[k]
                   + f_3 * pc_y[k] * ksi_621[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pa_x, pa_z, pc_x, pc_z, isk0_550, isk0_801, \
                         isi_426, isi_625, isk1_550, isk1_801, \
                         ksi_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = pa_x[k] * isk0_801[k]
                   + f_16 * isi_625[k]
                   - f_12 * pc_x[k] * isk1_801[k];

        t_802[k] = pa_z[k] * isk0_550[k]
                   - f_12 * pc_z[k] * isk1_550[k];

        t_803[k] = f_13 * isi_426[k]
                   + f_3 * pc_z[k] * ksi_622[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pa_x, pc_x, pc_y, isk0_804, isk0_806, isi_457, \
                         isi_628, isi_630, isk1_804, isk1_806, \
                         ksi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = pa_x[k] * isk0_804[k]
                   + f_15 * isi_628[k]
                   - f_12 * pc_x[k] * isk1_804[k];

        t_805[k] = f_17 * isi_457[k]
                   + f_3 * pc_y[k] * ksi_625[k];

        t_806[k] = pa_x[k] * isk0_806[k]
                   + f_15 * isi_630[k]
                   - f_12 * pc_x[k] * isk1_806[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, pa_x, pa_z, pc_x, pc_z, isk0_555, isk0_809, \
                         isi_430, isi_633, isk1_555, isk1_809, \
                         ksi_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pa_z[k] * isk0_555[k]
                   - f_12 * pc_z[k] * isk1_555[k];

        t_808[k] = f_13 * isi_430[k]
                   + f_3 * pc_z[k] * ksi_626[k];

        t_809[k] = pa_x[k] * isk0_809[k]
                   + f_14 * isi_633[k]
                   - f_12 * pc_x[k] * isk1_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, pa_x, pc_x, pc_y, isk0_810, isk0_812, isi_462, \
                         isi_634, isi_636, isk1_810, isk1_812, \
                         ksi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = pa_x[k] * isk0_810[k]
                   + f_14 * isi_634[k]
                   - f_12 * pc_x[k] * isk1_810[k];

        t_811[k] = f_17 * isi_462[k]
                   + f_3 * pc_y[k] * ksi_630[k];

        t_812[k] = pa_x[k] * isk0_812[k]
                   + f_14 * isi_636[k]
                   - f_12 * pc_x[k] * isk1_812[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, t_817, pc_x, isi_637, isi_638, isi_639, \
                         isi_640, isi_641, ksi_637, ksi_638, ksi_639, ksi_640, \
                         ksi_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_13 * isi_637[k]
                   + f_3 * pc_x[k] * ksi_637[k];

        t_814[k] = f_13 * isi_638[k]
                   + f_3 * pc_x[k] * ksi_638[k];

        t_815[k] = f_13 * isi_639[k]
                   + f_3 * pc_x[k] * ksi_639[k];

        t_816[k] = f_13 * isi_640[k]
                   + f_3 * pc_x[k] * ksi_640[k];

        t_817[k] = f_13 * isi_641[k]
                   + f_3 * pc_x[k] * ksi_641[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, pa_x, pc_x, pc_z, isk0_820, isi_441, \
                         isi_642, isi_643, isk1_820, ksi_637, ksi_642, \
                         ksi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_13 * isi_642[k]
                   + f_3 * pc_x[k] * ksi_642[k];

        t_819[k] = f_13 * isi_643[k]
                   + f_3 * pc_x[k] * ksi_643[k];

        t_820[k] = pa_x[k] * isk0_820[k]
                   - f_12 * pc_x[k] * isk1_820[k];

        t_821[k] = f_13 * isi_441[k]
                   + f_3 * pc_z[k] * ksi_637[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pa_x, pc_x, isk0_822, isk0_823, isk0_824, \
                         isk0_825, isk1_822, isk1_823, isk1_824, \
                         isk1_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pa_x[k] * isk0_822[k]
                   - f_12 * pc_x[k] * isk1_822[k];

        t_823[k] = pa_x[k] * isk0_823[k]
                   - f_12 * pc_x[k] * isk1_823[k];

        t_824[k] = pa_x[k] * isk0_824[k]
                   - f_12 * pc_x[k] * isk1_824[k];

        t_825[k] = pa_x[k] * isk0_825[k]
                   - f_12 * pc_x[k] * isk1_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pa_x, pc_x, pc_y, isk0_827, isk0_828, \
                         isi_475, isi_476, isi_644, isk1_827, isk1_828, ksi_643, \
                         ksi_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * isi_475[k]
                   + f_3 * pc_y[k] * ksi_643[k];

        t_827[k] = pa_x[k] * isk0_827[k]
                   - f_12 * pc_x[k] * isk1_827[k];

        t_828[k] = pa_x[k] * isk0_828[k]
                   + f_0 * isi_644[k]
                   - f_12 * pc_x[k] * isk1_828[k];

        t_829[k] = f_16 * isi_476[k]
                   + f_3 * pc_y[k] * ksi_644[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_3 = p / q;
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isk0_720 = buffer.data(isk0 + 720);
    const auto *isk0_725 = buffer.data(isk0 + 725);
    const auto *isk0_729 = buffer.data(isk0 + 729);
    const auto *isk0_831 = buffer.data(isk0 + 831);
    const auto *isk0_833 = buffer.data(isk0 + 833);
    const auto *isk0_834 = buffer.data(isk0 + 834);
    const auto *isk0_837 = buffer.data(isk0 + 837);
    const auto *isk0_838 = buffer.data(isk0 + 838);
    const auto *isk0_840 = buffer.data(isk0 + 840);
    const auto *isk0_842 = buffer.data(isk0 + 842);
    const auto *isk0_843 = buffer.data(isk0 + 843);
    const auto *isk0_845 = buffer.data(isk0 + 845);
    const auto *isk0_846 = buffer.data(isk0 + 846);
    const auto *isk0_848 = buffer.data(isk0 + 848);
    const auto *isk0_856 = buffer.data(isk0 + 856);
    const auto *isk0_858 = buffer.data(isk0 + 858);
    const auto *isk0_859 = buffer.data(isk0 + 859);
    const auto *isk0_860 = buffer.data(isk0 + 860);
    const auto *isk0_861 = buffer.data(isk0 + 861);
    const auto *isk0_863 = buffer.data(isk0 + 863);
    const auto *isk0_864 = buffer.data(isk0 + 864);
    const auto *isk0_867 = buffer.data(isk0 + 867);
    const auto *isk0_869 = buffer.data(isk0 + 869);
    const auto *isk0_870 = buffer.data(isk0 + 870);
    const auto *isk0_873 = buffer.data(isk0 + 873);
    const auto *isk0_874 = buffer.data(isk0 + 874);
    const auto *isk0_876 = buffer.data(isk0 + 876);
    const auto *isk0_878 = buffer.data(isk0 + 878);
    const auto *isk0_879 = buffer.data(isk0 + 879);
    const auto *isk0_881 = buffer.data(isk0 + 881);
    const auto *isk0_882 = buffer.data(isk0 + 882);
    const auto *isk0_884 = buffer.data(isk0 + 884);
    const auto *isk0_892 = buffer.data(isk0 + 892);
    const auto *isk0_894 = buffer.data(isk0 + 894);
    const auto *isk0_895 = buffer.data(isk0 + 895);
    const auto *isk0_896 = buffer.data(isk0 + 896);
    const auto *isk0_897 = buffer.data(isk0 + 897);
    const auto *isk0_899 = buffer.data(isk0 + 899);
    const auto *isk0_900 = buffer.data(isk0 + 900);
    const auto *isk0_903 = buffer.data(isk0 + 903);
    const auto *isk0_905 = buffer.data(isk0 + 905);
    const auto *isk0_906 = buffer.data(isk0 + 906);
    const auto *isk0_909 = buffer.data(isk0 + 909);
    const auto *isk0_910 = buffer.data(isk0 + 910);
    const auto *isk0_912 = buffer.data(isk0 + 912);
    const auto *isk0_914 = buffer.data(isk0 + 914);
    const auto *isk0_915 = buffer.data(isk0 + 915);
    const auto *isk0_917 = buffer.data(isk0 + 917);
    const auto *isk0_918 = buffer.data(isk0 + 918);
    const auto *isk0_920 = buffer.data(isk0 + 920);
    const auto *isk0_928 = buffer.data(isk0 + 928);
    const auto *isk0_930 = buffer.data(isk0 + 930);
    const auto *isk0_931 = buffer.data(isk0 + 931);
    const auto *isk0_932 = buffer.data(isk0 + 932);
    const auto *isk0_933 = buffer.data(isk0 + 933);
    const auto *isk0_935 = buffer.data(isk0 + 935);
    const auto *isk0_939 = buffer.data(isk0 + 939);
    const auto *isk0_942 = buffer.data(isk0 + 942);
    const auto *isk0_946 = buffer.data(isk0 + 946);

    const auto *isi_448 = buffer.data(isi + 448);
    const auto *isi_451 = buffer.data(isi + 451);
    const auto *isi_454 = buffer.data(isi + 454);
    const auto *isi_458 = buffer.data(isi + 458);
    const auto *isi_469 = buffer.data(isi + 469);
    const auto *isi_476 = buffer.data(isi + 476);
    const auto *isi_478 = buffer.data(isi + 478);
    const auto *isi_479 = buffer.data(isi + 479);
    const auto *isi_481 = buffer.data(isi + 481);
    const auto *isi_482 = buffer.data(isi + 482);
    const auto *isi_485 = buffer.data(isi + 485);
    const auto *isi_486 = buffer.data(isi + 486);
    const auto *isi_490 = buffer.data(isi + 490);
    const auto *isi_497 = buffer.data(isi + 497);
    const auto *isi_503 = buffer.data(isi + 503);
    const auto *isi_504 = buffer.data(isi + 504);
    const auto *isi_506 = buffer.data(isi + 506);
    const auto *isi_507 = buffer.data(isi + 507);
    const auto *isi_509 = buffer.data(isi + 509);
    const auto *isi_510 = buffer.data(isi + 510);
    const auto *isi_513 = buffer.data(isi + 513);
    const auto *isi_514 = buffer.data(isi + 514);
    const auto *isi_518 = buffer.data(isi + 518);
    const auto *isi_525 = buffer.data(isi + 525);
    const auto *isi_531 = buffer.data(isi + 531);
    const auto *isi_532 = buffer.data(isi + 532);
    const auto *isi_534 = buffer.data(isi + 534);
    const auto *isi_535 = buffer.data(isi + 535);
    const auto *isi_537 = buffer.data(isi + 537);
    const auto *isi_538 = buffer.data(isi + 538);
    const auto *isi_541 = buffer.data(isi + 541);
    const auto *isi_546 = buffer.data(isi + 546);
    const auto *isi_559 = buffer.data(isi + 559);
    const auto *isi_560 = buffer.data(isi + 560);
    const auto *isi_562 = buffer.data(isi + 562);
    const auto *isi_565 = buffer.data(isi + 565);
    const auto *isi_647 = buffer.data(isi + 647);
    const auto *isi_649 = buffer.data(isi + 649);
    const auto *isi_650 = buffer.data(isi + 650);
    const auto *isi_653 = buffer.data(isi + 653);
    const auto *isi_654 = buffer.data(isi + 654);
    const auto *isi_656 = buffer.data(isi + 656);
    const auto *isi_658 = buffer.data(isi + 658);
    const auto *isi_659 = buffer.data(isi + 659);
    const auto *isi_661 = buffer.data(isi + 661);
    const auto *isi_662 = buffer.data(isi + 662);
    const auto *isi_664 = buffer.data(isi + 664);
    const auto *isi_665 = buffer.data(isi + 665);
    const auto *isi_666 = buffer.data(isi + 666);
    const auto *isi_667 = buffer.data(isi + 667);
    const auto *isi_668 = buffer.data(isi + 668);
    const auto *isi_669 = buffer.data(isi + 669);
    const auto *isi_670 = buffer.data(isi + 670);
    const auto *isi_671 = buffer.data(isi + 671);
    const auto *isi_672 = buffer.data(isi + 672);
    const auto *isi_675 = buffer.data(isi + 675);
    const auto *isi_677 = buffer.data(isi + 677);
    const auto *isi_678 = buffer.data(isi + 678);
    const auto *isi_681 = buffer.data(isi + 681);
    const auto *isi_682 = buffer.data(isi + 682);
    const auto *isi_684 = buffer.data(isi + 684);
    const auto *isi_686 = buffer.data(isi + 686);
    const auto *isi_687 = buffer.data(isi + 687);
    const auto *isi_689 = buffer.data(isi + 689);
    const auto *isi_690 = buffer.data(isi + 690);
    const auto *isi_692 = buffer.data(isi + 692);
    const auto *isi_693 = buffer.data(isi + 693);
    const auto *isi_694 = buffer.data(isi + 694);
    const auto *isi_695 = buffer.data(isi + 695);
    const auto *isi_696 = buffer.data(isi + 696);
    const auto *isi_697 = buffer.data(isi + 697);
    const auto *isi_698 = buffer.data(isi + 698);
    const auto *isi_699 = buffer.data(isi + 699);
    const auto *isi_700 = buffer.data(isi + 700);
    const auto *isi_703 = buffer.data(isi + 703);
    const auto *isi_705 = buffer.data(isi + 705);
    const auto *isi_706 = buffer.data(isi + 706);
    const auto *isi_709 = buffer.data(isi + 709);
    const auto *isi_710 = buffer.data(isi + 710);
    const auto *isi_712 = buffer.data(isi + 712);
    const auto *isi_714 = buffer.data(isi + 714);
    const auto *isi_715 = buffer.data(isi + 715);
    const auto *isi_717 = buffer.data(isi + 717);
    const auto *isi_718 = buffer.data(isi + 718);
    const auto *isi_720 = buffer.data(isi + 720);
    const auto *isi_721 = buffer.data(isi + 721);
    const auto *isi_722 = buffer.data(isi + 722);
    const auto *isi_723 = buffer.data(isi + 723);
    const auto *isi_724 = buffer.data(isi + 724);
    const auto *isi_725 = buffer.data(isi + 725);
    const auto *isi_726 = buffer.data(isi + 726);
    const auto *isi_727 = buffer.data(isi + 727);
    const auto *isi_731 = buffer.data(isi + 731);
    const auto *isi_734 = buffer.data(isi + 734);
    const auto *isi_738 = buffer.data(isi + 738);

    const auto *isk1_720 = buffer.data(isk1 + 720);
    const auto *isk1_725 = buffer.data(isk1 + 725);
    const auto *isk1_729 = buffer.data(isk1 + 729);
    const auto *isk1_831 = buffer.data(isk1 + 831);
    const auto *isk1_833 = buffer.data(isk1 + 833);
    const auto *isk1_834 = buffer.data(isk1 + 834);
    const auto *isk1_837 = buffer.data(isk1 + 837);
    const auto *isk1_838 = buffer.data(isk1 + 838);
    const auto *isk1_840 = buffer.data(isk1 + 840);
    const auto *isk1_842 = buffer.data(isk1 + 842);
    const auto *isk1_843 = buffer.data(isk1 + 843);
    const auto *isk1_845 = buffer.data(isk1 + 845);
    const auto *isk1_846 = buffer.data(isk1 + 846);
    const auto *isk1_848 = buffer.data(isk1 + 848);
    const auto *isk1_856 = buffer.data(isk1 + 856);
    const auto *isk1_858 = buffer.data(isk1 + 858);
    const auto *isk1_859 = buffer.data(isk1 + 859);
    const auto *isk1_860 = buffer.data(isk1 + 860);
    const auto *isk1_861 = buffer.data(isk1 + 861);
    const auto *isk1_863 = buffer.data(isk1 + 863);
    const auto *isk1_864 = buffer.data(isk1 + 864);
    const auto *isk1_867 = buffer.data(isk1 + 867);
    const auto *isk1_869 = buffer.data(isk1 + 869);
    const auto *isk1_870 = buffer.data(isk1 + 870);
    const auto *isk1_873 = buffer.data(isk1 + 873);
    const auto *isk1_874 = buffer.data(isk1 + 874);
    const auto *isk1_876 = buffer.data(isk1 + 876);
    const auto *isk1_878 = buffer.data(isk1 + 878);
    const auto *isk1_879 = buffer.data(isk1 + 879);
    const auto *isk1_881 = buffer.data(isk1 + 881);
    const auto *isk1_882 = buffer.data(isk1 + 882);
    const auto *isk1_884 = buffer.data(isk1 + 884);
    const auto *isk1_892 = buffer.data(isk1 + 892);
    const auto *isk1_894 = buffer.data(isk1 + 894);
    const auto *isk1_895 = buffer.data(isk1 + 895);
    const auto *isk1_896 = buffer.data(isk1 + 896);
    const auto *isk1_897 = buffer.data(isk1 + 897);
    const auto *isk1_899 = buffer.data(isk1 + 899);
    const auto *isk1_900 = buffer.data(isk1 + 900);
    const auto *isk1_903 = buffer.data(isk1 + 903);
    const auto *isk1_905 = buffer.data(isk1 + 905);
    const auto *isk1_906 = buffer.data(isk1 + 906);
    const auto *isk1_909 = buffer.data(isk1 + 909);
    const auto *isk1_910 = buffer.data(isk1 + 910);
    const auto *isk1_912 = buffer.data(isk1 + 912);
    const auto *isk1_914 = buffer.data(isk1 + 914);
    const auto *isk1_915 = buffer.data(isk1 + 915);
    const auto *isk1_917 = buffer.data(isk1 + 917);
    const auto *isk1_918 = buffer.data(isk1 + 918);
    const auto *isk1_920 = buffer.data(isk1 + 920);
    const auto *isk1_928 = buffer.data(isk1 + 928);
    const auto *isk1_930 = buffer.data(isk1 + 930);
    const auto *isk1_931 = buffer.data(isk1 + 931);
    const auto *isk1_932 = buffer.data(isk1 + 932);
    const auto *isk1_933 = buffer.data(isk1 + 933);
    const auto *isk1_935 = buffer.data(isk1 + 935);
    const auto *isk1_939 = buffer.data(isk1 + 939);
    const auto *isk1_942 = buffer.data(isk1 + 942);
    const auto *isk1_946 = buffer.data(isk1 + 946);

    const auto *ksi_644 = buffer.data(ksi + 644);
    const auto *ksi_646 = buffer.data(ksi + 646);
    const auto *ksi_647 = buffer.data(ksi + 647);
    const auto *ksi_649 = buffer.data(ksi + 649);
    const auto *ksi_650 = buffer.data(ksi + 650);
    const auto *ksi_653 = buffer.data(ksi + 653);
    const auto *ksi_654 = buffer.data(ksi + 654);
    const auto *ksi_658 = buffer.data(ksi + 658);
    const auto *ksi_665 = buffer.data(ksi + 665);
    const auto *ksi_666 = buffer.data(ksi + 666);
    const auto *ksi_667 = buffer.data(ksi + 667);
    const auto *ksi_668 = buffer.data(ksi + 668);
    const auto *ksi_669 = buffer.data(ksi + 669);
    const auto *ksi_670 = buffer.data(ksi + 670);
    const auto *ksi_671 = buffer.data(ksi + 671);
    const auto *ksi_672 = buffer.data(ksi + 672);
    const auto *ksi_674 = buffer.data(ksi + 674);
    const auto *ksi_675 = buffer.data(ksi + 675);
    const auto *ksi_677 = buffer.data(ksi + 677);
    const auto *ksi_678 = buffer.data(ksi + 678);
    const auto *ksi_681 = buffer.data(ksi + 681);
    const auto *ksi_682 = buffer.data(ksi + 682);
    const auto *ksi_686 = buffer.data(ksi + 686);
    const auto *ksi_693 = buffer.data(ksi + 693);
    const auto *ksi_694 = buffer.data(ksi + 694);
    const auto *ksi_695 = buffer.data(ksi + 695);
    const auto *ksi_696 = buffer.data(ksi + 696);
    const auto *ksi_697 = buffer.data(ksi + 697);
    const auto *ksi_698 = buffer.data(ksi + 698);
    const auto *ksi_699 = buffer.data(ksi + 699);
    const auto *ksi_700 = buffer.data(ksi + 700);
    const auto *ksi_702 = buffer.data(ksi + 702);
    const auto *ksi_703 = buffer.data(ksi + 703);
    const auto *ksi_705 = buffer.data(ksi + 705);
    const auto *ksi_706 = buffer.data(ksi + 706);
    const auto *ksi_709 = buffer.data(ksi + 709);
    const auto *ksi_710 = buffer.data(ksi + 710);
    const auto *ksi_714 = buffer.data(ksi + 714);
    const auto *ksi_721 = buffer.data(ksi + 721);
    const auto *ksi_722 = buffer.data(ksi + 722);
    const auto *ksi_723 = buffer.data(ksi + 723);
    const auto *ksi_724 = buffer.data(ksi + 724);
    const auto *ksi_725 = buffer.data(ksi + 725);
    const auto *ksi_726 = buffer.data(ksi + 726);
    const auto *ksi_727 = buffer.data(ksi + 727);
    const auto *ksi_728 = buffer.data(ksi + 728);
    const auto *ksi_730 = buffer.data(ksi + 730);
    const auto *ksi_731 = buffer.data(ksi + 731);
    const auto *ksi_733 = buffer.data(ksi + 733);
    const auto *ksi_734 = buffer.data(ksi + 734);

#pragma omp simd aligned(t_830, t_831, t_832, pa_x, pc_x, pc_y, pc_z, isk0_831, isi_448, \
                         isi_478, isi_647, isk1_831, ksi_644, ksi_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_14 * isi_448[k]
                   + f_3 * pc_z[k] * ksi_644[k];

        t_831[k] = pa_x[k] * isk0_831[k]
                   + f_17 * isi_647[k]
                   - f_12 * pc_x[k] * isk1_831[k];

        t_832[k] = f_16 * isi_478[k]
                   + f_3 * pc_y[k] * ksi_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pa_x, pc_x, pc_z, isk0_833, isk0_834, isi_451, \
                         isi_649, isi_650, isk1_833, isk1_834, \
                         ksi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = pa_x[k] * isk0_833[k]
                   + f_17 * isi_649[k]
                   - f_12 * pc_x[k] * isk1_833[k];

        t_834[k] = pa_x[k] * isk0_834[k]
                   + f_16 * isi_650[k]
                   - f_12 * pc_x[k] * isk1_834[k];

        t_835[k] = f_14 * isi_451[k]
                   + f_3 * pc_z[k] * ksi_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pa_x, pc_x, pc_y, isk0_837, isk0_838, isi_481, \
                         isi_653, isi_654, isk1_837, isk1_838, \
                         ksi_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * isi_481[k]
                   + f_3 * pc_y[k] * ksi_649[k];

        t_837[k] = pa_x[k] * isk0_837[k]
                   + f_16 * isi_653[k]
                   - f_12 * pc_x[k] * isk1_837[k];

        t_838[k] = pa_x[k] * isk0_838[k]
                   + f_15 * isi_654[k]
                   - f_12 * pc_x[k] * isk1_838[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pa_x, pc_x, pc_y, pc_z, isk0_840, isi_454, \
                         isi_485, isi_656, isk1_840, ksi_650, ksi_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * isi_454[k]
                   + f_3 * pc_z[k] * ksi_650[k];

        t_840[k] = pa_x[k] * isk0_840[k]
                   + f_15 * isi_656[k]
                   - f_12 * pc_x[k] * isk1_840[k];

        t_841[k] = f_16 * isi_485[k]
                   + f_3 * pc_y[k] * ksi_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pa_x, pc_x, pc_z, isk0_842, isk0_843, isi_458, \
                         isi_658, isi_659, isk1_842, isk1_843, \
                         ksi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = pa_x[k] * isk0_842[k]
                   + f_15 * isi_658[k]
                   - f_12 * pc_x[k] * isk1_842[k];

        t_843[k] = pa_x[k] * isk0_843[k]
                   + f_14 * isi_659[k]
                   - f_12 * pc_x[k] * isk1_843[k];

        t_844[k] = f_14 * isi_458[k]
                   + f_3 * pc_z[k] * ksi_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pa_x, pc_x, pc_y, isk0_845, isk0_846, isi_490, \
                         isi_661, isi_662, isk1_845, isk1_846, \
                         ksi_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = pa_x[k] * isk0_845[k]
                   + f_14 * isi_661[k]
                   - f_12 * pc_x[k] * isk1_845[k];

        t_846[k] = pa_x[k] * isk0_846[k]
                   + f_14 * isi_662[k]
                   - f_12 * pc_x[k] * isk1_846[k];

        t_847[k] = f_16 * isi_490[k]
                   + f_3 * pc_y[k] * ksi_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pa_x, pc_x, isk0_848, isi_664, isi_665, \
                         isi_666, isi_667, isk1_848, ksi_665, ksi_666, \
                         ksi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = pa_x[k] * isk0_848[k]
                   + f_14 * isi_664[k]
                   - f_12 * pc_x[k] * isk1_848[k];

        t_849[k] = f_13 * isi_665[k]
                   + f_3 * pc_x[k] * ksi_665[k];

        t_850[k] = f_13 * isi_666[k]
                   + f_3 * pc_x[k] * ksi_666[k];

        t_851[k] = f_13 * isi_667[k]
                   + f_3 * pc_x[k] * ksi_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, isi_668, isi_669, isi_670, isi_671, \
                         ksi_668, ksi_669, ksi_670, ksi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_13 * isi_668[k]
                   + f_3 * pc_x[k] * ksi_668[k];

        t_853[k] = f_13 * isi_669[k]
                   + f_3 * pc_x[k] * ksi_669[k];

        t_854[k] = f_13 * isi_670[k]
                   + f_3 * pc_x[k] * ksi_670[k];

        t_855[k] = f_13 * isi_671[k]
                   + f_3 * pc_x[k] * ksi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, t_859, pa_x, pc_x, pc_z, isk0_856, isk0_858, \
                         isk0_859, isi_469, isk1_856, isk1_858, isk1_859, \
                         ksi_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = pa_x[k] * isk0_856[k]
                   - f_12 * pc_x[k] * isk1_856[k];

        t_857[k] = f_14 * isi_469[k]
                   + f_3 * pc_z[k] * ksi_665[k];

        t_858[k] = pa_x[k] * isk0_858[k]
                   - f_12 * pc_x[k] * isk1_858[k];

        t_859[k] = pa_x[k] * isk0_859[k]
                   - f_12 * pc_x[k] * isk1_859[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, pa_x, pc_x, pc_y, isk0_860, isk0_861, \
                         isk0_863, isi_503, isk1_860, isk1_861, isk1_863, \
                         ksi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = pa_x[k] * isk0_860[k]
                   - f_12 * pc_x[k] * isk1_860[k];

        t_861[k] = pa_x[k] * isk0_861[k]
                   - f_12 * pc_x[k] * isk1_861[k];

        t_862[k] = f_16 * isi_503[k]
                   + f_3 * pc_y[k] * ksi_671[k];

        t_863[k] = pa_x[k] * isk0_863[k]
                   - f_12 * pc_x[k] * isk1_863[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pa_x, pc_x, pc_y, pc_z, isk0_864, isi_476, \
                         isi_504, isi_672, isk1_864, ksi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = pa_x[k] * isk0_864[k]
                   + f_0 * isi_672[k]
                   - f_12 * pc_x[k] * isk1_864[k];

        t_865[k] = f_15 * isi_504[k]
                   + f_3 * pc_y[k] * ksi_672[k];

        t_866[k] = f_15 * isi_476[k]
                   + f_3 * pc_z[k] * ksi_672[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pa_x, pc_x, pc_y, isk0_867, isk0_869, isi_506, \
                         isi_675, isi_677, isk1_867, isk1_869, \
                         ksi_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = pa_x[k] * isk0_867[k]
                   + f_17 * isi_675[k]
                   - f_12 * pc_x[k] * isk1_867[k];

        t_868[k] = f_15 * isi_506[k]
                   + f_3 * pc_y[k] * ksi_674[k];

        t_869[k] = pa_x[k] * isk0_869[k]
                   + f_17 * isi_677[k]
                   - f_12 * pc_x[k] * isk1_869[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pa_x, pc_x, pc_y, pc_z, isk0_870, isi_479, \
                         isi_509, isi_678, isk1_870, ksi_675, ksi_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_x[k] * isk0_870[k]
                   + f_16 * isi_678[k]
                   - f_12 * pc_x[k] * isk1_870[k];

        t_871[k] = f_15 * isi_479[k]
                   + f_3 * pc_z[k] * ksi_675[k];

        t_872[k] = f_15 * isi_509[k]
                   + f_3 * pc_y[k] * ksi_677[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pa_x, pc_x, pc_z, isk0_873, isk0_874, isi_482, \
                         isi_681, isi_682, isk1_873, isk1_874, \
                         ksi_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pa_x[k] * isk0_873[k]
                   + f_16 * isi_681[k]
                   - f_12 * pc_x[k] * isk1_873[k];

        t_874[k] = pa_x[k] * isk0_874[k]
                   + f_15 * isi_682[k]
                   - f_12 * pc_x[k] * isk1_874[k];

        t_875[k] = f_15 * isi_482[k]
                   + f_3 * pc_z[k] * ksi_678[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pa_x, pc_x, pc_y, isk0_876, isk0_878, isi_513, \
                         isi_684, isi_686, isk1_876, isk1_878, \
                         ksi_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = pa_x[k] * isk0_876[k]
                   + f_15 * isi_684[k]
                   - f_12 * pc_x[k] * isk1_876[k];

        t_877[k] = f_15 * isi_513[k]
                   + f_3 * pc_y[k] * ksi_681[k];

        t_878[k] = pa_x[k] * isk0_878[k]
                   + f_15 * isi_686[k]
                   - f_12 * pc_x[k] * isk1_878[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pa_x, pc_x, pc_z, isk0_879, isk0_881, isi_486, \
                         isi_687, isi_689, isk1_879, isk1_881, \
                         ksi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = pa_x[k] * isk0_879[k]
                   + f_14 * isi_687[k]
                   - f_12 * pc_x[k] * isk1_879[k];

        t_880[k] = f_15 * isi_486[k]
                   + f_3 * pc_z[k] * ksi_682[k];

        t_881[k] = pa_x[k] * isk0_881[k]
                   + f_14 * isi_689[k]
                   - f_12 * pc_x[k] * isk1_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, pa_x, pc_x, pc_y, isk0_882, isk0_884, isi_518, \
                         isi_690, isi_692, isk1_882, isk1_884, \
                         ksi_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pa_x[k] * isk0_882[k]
                   + f_14 * isi_690[k]
                   - f_12 * pc_x[k] * isk1_882[k];

        t_883[k] = f_15 * isi_518[k]
                   + f_3 * pc_y[k] * ksi_686[k];

        t_884[k] = pa_x[k] * isk0_884[k]
                   + f_14 * isi_692[k]
                   - f_12 * pc_x[k] * isk1_884[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, pc_x, isi_693, isi_694, isi_695, \
                         isi_696, isi_697, ksi_693, ksi_694, ksi_695, ksi_696, \
                         ksi_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = f_13 * isi_693[k]
                   + f_3 * pc_x[k] * ksi_693[k];

        t_886[k] = f_13 * isi_694[k]
                   + f_3 * pc_x[k] * ksi_694[k];

        t_887[k] = f_13 * isi_695[k]
                   + f_3 * pc_x[k] * ksi_695[k];

        t_888[k] = f_13 * isi_696[k]
                   + f_3 * pc_x[k] * ksi_696[k];

        t_889[k] = f_13 * isi_697[k]
                   + f_3 * pc_x[k] * ksi_697[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pa_x, pc_x, pc_z, isk0_892, isi_497, \
                         isi_698, isi_699, isk1_892, ksi_693, ksi_698, \
                         ksi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * isi_698[k]
                   + f_3 * pc_x[k] * ksi_698[k];

        t_891[k] = f_13 * isi_699[k]
                   + f_3 * pc_x[k] * ksi_699[k];

        t_892[k] = pa_x[k] * isk0_892[k]
                   - f_12 * pc_x[k] * isk1_892[k];

        t_893[k] = f_15 * isi_497[k]
                   + f_3 * pc_z[k] * ksi_693[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pa_x, pc_x, isk0_894, isk0_895, isk0_896, \
                         isk0_897, isk1_894, isk1_895, isk1_896, \
                         isk1_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = pa_x[k] * isk0_894[k]
                   - f_12 * pc_x[k] * isk1_894[k];

        t_895[k] = pa_x[k] * isk0_895[k]
                   - f_12 * pc_x[k] * isk1_895[k];

        t_896[k] = pa_x[k] * isk0_896[k]
                   - f_12 * pc_x[k] * isk1_896[k];

        t_897[k] = pa_x[k] * isk0_897[k]
                   - f_12 * pc_x[k] * isk1_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pa_x, pc_x, pc_y, isk0_899, isk0_900, \
                         isi_531, isi_532, isi_700, isk1_899, isk1_900, ksi_699, \
                         ksi_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * isi_531[k]
                   + f_3 * pc_y[k] * ksi_699[k];

        t_899[k] = pa_x[k] * isk0_899[k]
                   - f_12 * pc_x[k] * isk1_899[k];

        t_900[k] = pa_x[k] * isk0_900[k]
                   + f_0 * isi_700[k]
                   - f_12 * pc_x[k] * isk1_900[k];

        t_901[k] = f_14 * isi_532[k]
                   + f_3 * pc_y[k] * ksi_700[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pa_x, pc_x, pc_y, pc_z, isk0_903, isi_504, \
                         isi_534, isi_703, isk1_903, ksi_700, ksi_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_16 * isi_504[k]
                   + f_3 * pc_z[k] * ksi_700[k];

        t_903[k] = pa_x[k] * isk0_903[k]
                   + f_17 * isi_703[k]
                   - f_12 * pc_x[k] * isk1_903[k];

        t_904[k] = f_14 * isi_534[k]
                   + f_3 * pc_y[k] * ksi_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pa_x, pc_x, pc_z, isk0_905, isk0_906, isi_507, \
                         isi_705, isi_706, isk1_905, isk1_906, \
                         ksi_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = pa_x[k] * isk0_905[k]
                   + f_17 * isi_705[k]
                   - f_12 * pc_x[k] * isk1_905[k];

        t_906[k] = pa_x[k] * isk0_906[k]
                   + f_16 * isi_706[k]
                   - f_12 * pc_x[k] * isk1_906[k];

        t_907[k] = f_16 * isi_507[k]
                   + f_3 * pc_z[k] * ksi_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pa_x, pc_x, pc_y, isk0_909, isk0_910, isi_537, \
                         isi_709, isi_710, isk1_909, isk1_910, \
                         ksi_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * isi_537[k]
                   + f_3 * pc_y[k] * ksi_705[k];

        t_909[k] = pa_x[k] * isk0_909[k]
                   + f_16 * isi_709[k]
                   - f_12 * pc_x[k] * isk1_909[k];

        t_910[k] = pa_x[k] * isk0_910[k]
                   + f_15 * isi_710[k]
                   - f_12 * pc_x[k] * isk1_910[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pa_x, pc_x, pc_y, pc_z, isk0_912, isi_510, \
                         isi_541, isi_712, isk1_912, ksi_706, ksi_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * isi_510[k]
                   + f_3 * pc_z[k] * ksi_706[k];

        t_912[k] = pa_x[k] * isk0_912[k]
                   + f_15 * isi_712[k]
                   - f_12 * pc_x[k] * isk1_912[k];

        t_913[k] = f_14 * isi_541[k]
                   + f_3 * pc_y[k] * ksi_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pa_x, pc_x, pc_z, isk0_914, isk0_915, isi_514, \
                         isi_714, isi_715, isk1_914, isk1_915, \
                         ksi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = pa_x[k] * isk0_914[k]
                   + f_15 * isi_714[k]
                   - f_12 * pc_x[k] * isk1_914[k];

        t_915[k] = pa_x[k] * isk0_915[k]
                   + f_14 * isi_715[k]
                   - f_12 * pc_x[k] * isk1_915[k];

        t_916[k] = f_16 * isi_514[k]
                   + f_3 * pc_z[k] * ksi_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pa_x, pc_x, pc_y, isk0_917, isk0_918, isi_546, \
                         isi_717, isi_718, isk1_917, isk1_918, \
                         ksi_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = pa_x[k] * isk0_917[k]
                   + f_14 * isi_717[k]
                   - f_12 * pc_x[k] * isk1_917[k];

        t_918[k] = pa_x[k] * isk0_918[k]
                   + f_14 * isi_718[k]
                   - f_12 * pc_x[k] * isk1_918[k];

        t_919[k] = f_14 * isi_546[k]
                   + f_3 * pc_y[k] * ksi_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pa_x, pc_x, isk0_920, isi_720, isi_721, \
                         isi_722, isi_723, isk1_920, ksi_721, ksi_722, \
                         ksi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pa_x[k] * isk0_920[k]
                   + f_14 * isi_720[k]
                   - f_12 * pc_x[k] * isk1_920[k];

        t_921[k] = f_13 * isi_721[k]
                   + f_3 * pc_x[k] * ksi_721[k];

        t_922[k] = f_13 * isi_722[k]
                   + f_3 * pc_x[k] * ksi_722[k];

        t_923[k] = f_13 * isi_723[k]
                   + f_3 * pc_x[k] * ksi_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, isi_724, isi_725, isi_726, isi_727, \
                         ksi_724, ksi_725, ksi_726, ksi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_13 * isi_724[k]
                   + f_3 * pc_x[k] * ksi_724[k];

        t_925[k] = f_13 * isi_725[k]
                   + f_3 * pc_x[k] * ksi_725[k];

        t_926[k] = f_13 * isi_726[k]
                   + f_3 * pc_x[k] * ksi_726[k];

        t_927[k] = f_13 * isi_727[k]
                   + f_3 * pc_x[k] * ksi_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, pa_x, pc_x, pc_z, isk0_928, isk0_930, \
                         isk0_931, isi_525, isk1_928, isk1_930, isk1_931, \
                         ksi_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = pa_x[k] * isk0_928[k]
                   - f_12 * pc_x[k] * isk1_928[k];

        t_929[k] = f_16 * isi_525[k]
                   + f_3 * pc_z[k] * ksi_721[k];

        t_930[k] = pa_x[k] * isk0_930[k]
                   - f_12 * pc_x[k] * isk1_930[k];

        t_931[k] = pa_x[k] * isk0_931[k]
                   - f_12 * pc_x[k] * isk1_931[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, t_935, pa_x, pc_x, pc_y, isk0_932, isk0_933, \
                         isk0_935, isi_559, isk1_932, isk1_933, isk1_935, \
                         ksi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = pa_x[k] * isk0_932[k]
                   - f_12 * pc_x[k] * isk1_932[k];

        t_933[k] = pa_x[k] * isk0_933[k]
                   - f_12 * pc_x[k] * isk1_933[k];

        t_934[k] = f_14 * isi_559[k]
                   + f_3 * pc_y[k] * ksi_727[k];

        t_935[k] = pa_x[k] * isk0_935[k]
                   - f_12 * pc_x[k] * isk1_935[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, pa_y, pc_y, pc_z, isk0_720, isi_532, isi_560, \
                         isk1_720, ksi_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = pa_y[k] * isk0_720[k]
                   - f_12 * pc_y[k] * isk1_720[k];

        t_937[k] = f_13 * isi_560[k]
                   + f_3 * pc_y[k] * ksi_728[k];

        t_938[k] = f_17 * isi_532[k]
                   + f_3 * pc_z[k] * ksi_728[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, pa_x, pa_y, pc_x, pc_y, isk0_725, isk0_939, \
                         isi_562, isi_731, isk1_725, isk1_939, \
                         ksi_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = pa_x[k] * isk0_939[k]
                   + f_17 * isi_731[k]
                   - f_12 * pc_x[k] * isk1_939[k];

        t_940[k] = f_13 * isi_562[k]
                   + f_3 * pc_y[k] * ksi_730[k];

        t_941[k] = pa_y[k] * isk0_725[k]
                   - f_12 * pc_y[k] * isk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, pa_x, pc_x, pc_y, pc_z, isk0_942, isi_535, \
                         isi_565, isi_734, isk1_942, ksi_731, ksi_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pa_x[k] * isk0_942[k]
                   + f_16 * isi_734[k]
                   - f_12 * pc_x[k] * isk1_942[k];

        t_943[k] = f_17 * isi_535[k]
                   + f_3 * pc_z[k] * ksi_731[k];

        t_944[k] = f_13 * isi_565[k]
                   + f_3 * pc_y[k] * ksi_733[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, pa_x, pa_y, pc_x, pc_y, pc_z, isk0_729, \
                         isk0_946, isi_538, isi_738, isk1_729, isk1_946, \
                         ksi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = pa_y[k] * isk0_729[k]
                   - f_12 * pc_y[k] * isk1_729[k];

        t_946[k] = pa_x[k] * isk0_946[k]
                   + f_15 * isi_738[k]
                   - f_12 * pc_x[k] * isk1_946[k];

        t_947[k] = f_17 * isi_538[k]
                   + f_3 * pc_z[k] * ksi_734[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isk0_734 = buffer.data(isk0 + 734);
    const auto *isk0_740 = buffer.data(isk0 + 740);
    const auto *isk0_756 = buffer.data(isk0 + 756);
    const auto *isk0_757 = buffer.data(isk0 + 757);
    const auto *isk0_759 = buffer.data(isk0 + 759);
    const auto *isk0_762 = buffer.data(isk0 + 762);
    const auto *isk0_766 = buffer.data(isk0 + 766);
    const auto *isk0_771 = buffer.data(isk0 + 771);
    const auto *isk0_784 = buffer.data(isk0 + 784);
    const auto *isk0_948 = buffer.data(isk0 + 948);
    const auto *isk0_951 = buffer.data(isk0 + 951);
    const auto *isk0_953 = buffer.data(isk0 + 953);
    const auto *isk0_954 = buffer.data(isk0 + 954);
    const auto *isk0_964 = buffer.data(isk0 + 964);
    const auto *isk0_966 = buffer.data(isk0 + 966);
    const auto *isk0_967 = buffer.data(isk0 + 967);
    const auto *isk0_968 = buffer.data(isk0 + 968);
    const auto *isk0_969 = buffer.data(isk0 + 969);
    const auto *isk0_971 = buffer.data(isk0 + 971);
    const auto *isk0_972 = buffer.data(isk0 + 972);
    const auto *isk0_977 = buffer.data(isk0 + 977);
    const auto *isk0_981 = buffer.data(isk0 + 981);
    const auto *isk0_986 = buffer.data(isk0 + 986);
    const auto *isk0_992 = buffer.data(isk0 + 992);
    const auto *isk0_1000 = buffer.data(isk0 + 1000);
    const auto *isk0_1001 = buffer.data(isk0 + 1001);
    const auto *isk0_1002 = buffer.data(isk0 + 1002);
    const auto *isk0_1003 = buffer.data(isk0 + 1003);
    const auto *isk0_1004 = buffer.data(isk0 + 1004);
    const auto *isk0_1005 = buffer.data(isk0 + 1005);
    const auto *isk0_1007 = buffer.data(isk0 + 1007);

    const auto *isi_542 = buffer.data(isi + 542);
    const auto *isi_553 = buffer.data(isi + 553);
    const auto *isi_560 = buffer.data(isi + 560);
    const auto *isi_569 = buffer.data(isi + 569);
    const auto *isi_574 = buffer.data(isi + 574);
    const auto *isi_587 = buffer.data(isi + 587);
    const auto *isi_609 = buffer.data(isi + 609);
    const auto *isi_615 = buffer.data(isi + 615);
    const auto *isi_740 = buffer.data(isi + 740);
    const auto *isi_743 = buffer.data(isi + 743);
    const auto *isi_745 = buffer.data(isi + 745);
    const auto *isi_746 = buffer.data(isi + 746);
    const auto *isi_749 = buffer.data(isi + 749);
    const auto *isi_750 = buffer.data(isi + 750);
    const auto *isi_751 = buffer.data(isi + 751);
    const auto *isi_752 = buffer.data(isi + 752);
    const auto *isi_753 = buffer.data(isi + 753);
    const auto *isi_754 = buffer.data(isi + 754);
    const auto *isi_755 = buffer.data(isi + 755);
    const auto *isi_756 = buffer.data(isi + 756);
    const auto *isi_761 = buffer.data(isi + 761);
    const auto *isi_765 = buffer.data(isi + 765);
    const auto *isi_770 = buffer.data(isi + 770);
    const auto *isi_776 = buffer.data(isi + 776);
    const auto *isi_777 = buffer.data(isi + 777);
    const auto *isi_778 = buffer.data(isi + 778);
    const auto *isi_779 = buffer.data(isi + 779);
    const auto *isi_780 = buffer.data(isi + 780);
    const auto *isi_781 = buffer.data(isi + 781);
    const auto *isi_783 = buffer.data(isi + 783);

    const auto *isk1_734 = buffer.data(isk1 + 734);
    const auto *isk1_740 = buffer.data(isk1 + 740);
    const auto *isk1_756 = buffer.data(isk1 + 756);
    const auto *isk1_757 = buffer.data(isk1 + 757);
    const auto *isk1_759 = buffer.data(isk1 + 759);
    const auto *isk1_762 = buffer.data(isk1 + 762);
    const auto *isk1_766 = buffer.data(isk1 + 766);
    const auto *isk1_771 = buffer.data(isk1 + 771);
    const auto *isk1_784 = buffer.data(isk1 + 784);
    const auto *isk1_948 = buffer.data(isk1 + 948);
    const auto *isk1_951 = buffer.data(isk1 + 951);
    const auto *isk1_953 = buffer.data(isk1 + 953);
    const auto *isk1_954 = buffer.data(isk1 + 954);
    const auto *isk1_964 = buffer.data(isk1 + 964);
    const auto *isk1_966 = buffer.data(isk1 + 966);
    const auto *isk1_967 = buffer.data(isk1 + 967);
    const auto *isk1_968 = buffer.data(isk1 + 968);
    const auto *isk1_969 = buffer.data(isk1 + 969);
    const auto *isk1_971 = buffer.data(isk1 + 971);
    const auto *isk1_972 = buffer.data(isk1 + 972);
    const auto *isk1_977 = buffer.data(isk1 + 977);
    const auto *isk1_981 = buffer.data(isk1 + 981);
    const auto *isk1_986 = buffer.data(isk1 + 986);
    const auto *isk1_992 = buffer.data(isk1 + 992);
    const auto *isk1_1000 = buffer.data(isk1 + 1000);
    const auto *isk1_1001 = buffer.data(isk1 + 1001);
    const auto *isk1_1002 = buffer.data(isk1 + 1002);
    const auto *isk1_1003 = buffer.data(isk1 + 1003);
    const auto *isk1_1004 = buffer.data(isk1 + 1004);
    const auto *isk1_1005 = buffer.data(isk1 + 1005);
    const auto *isk1_1007 = buffer.data(isk1 + 1007);

    const auto *ksh0_567 = buffer.data(ksh0 + 567);
    const auto *ksh0_568 = buffer.data(ksh0 + 568);
    const auto *ksh0_569 = buffer.data(ksh0 + 569);
    const auto *ksh0_570 = buffer.data(ksh0 + 570);
    const auto *ksh0_571 = buffer.data(ksh0 + 571);
    const auto *ksh0_572 = buffer.data(ksh0 + 572);
    const auto *ksh0_573 = buffer.data(ksh0 + 573);
    const auto *ksh0_574 = buffer.data(ksh0 + 574);
    const auto *ksh0_575 = buffer.data(ksh0 + 575);
    const auto *ksh0_576 = buffer.data(ksh0 + 576);
    const auto *ksh0_588 = buffer.data(ksh0 + 588);
    const auto *ksh0_589 = buffer.data(ksh0 + 589);
    const auto *ksh0_591 = buffer.data(ksh0 + 591);
    const auto *ksh0_593 = buffer.data(ksh0 + 593);
    const auto *ksh0_594 = buffer.data(ksh0 + 594);
    const auto *ksh0_596 = buffer.data(ksh0 + 596);
    const auto *ksh0_597 = buffer.data(ksh0 + 597);
    const auto *ksh0_598 = buffer.data(ksh0 + 598);
    const auto *ksh0_600 = buffer.data(ksh0 + 600);
    const auto *ksh0_601 = buffer.data(ksh0 + 601);
    const auto *ksh0_602 = buffer.data(ksh0 + 602);
    const auto *ksh0_603 = buffer.data(ksh0 + 603);
    const auto *ksh0_604 = buffer.data(ksh0 + 604);
    const auto *ksh0_605 = buffer.data(ksh0 + 605);
    const auto *ksh0_606 = buffer.data(ksh0 + 606);
    const auto *ksh0_607 = buffer.data(ksh0 + 607);
    const auto *ksh0_608 = buffer.data(ksh0 + 608);
    const auto *ksh0_611 = buffer.data(ksh0 + 611);
    const auto *ksh0_613 = buffer.data(ksh0 + 613);
    const auto *ksh0_614 = buffer.data(ksh0 + 614);
    const auto *ksh0_616 = buffer.data(ksh0 + 616);
    const auto *ksh0_617 = buffer.data(ksh0 + 617);
    const auto *ksh0_618 = buffer.data(ksh0 + 618);
    const auto *ksh0_620 = buffer.data(ksh0 + 620);
    const auto *ksh0_621 = buffer.data(ksh0 + 621);
    const auto *ksh0_622 = buffer.data(ksh0 + 622);
    const auto *ksh0_623 = buffer.data(ksh0 + 623);
    const auto *ksh0_625 = buffer.data(ksh0 + 625);
    const auto *ksh0_626 = buffer.data(ksh0 + 626);
    const auto *ksh0_627 = buffer.data(ksh0 + 627);
    const auto *ksh0_628 = buffer.data(ksh0 + 628);
    const auto *ksh0_629 = buffer.data(ksh0 + 629);

    const auto *ksh1_567 = buffer.data(ksh1 + 567);
    const auto *ksh1_568 = buffer.data(ksh1 + 568);
    const auto *ksh1_569 = buffer.data(ksh1 + 569);
    const auto *ksh1_570 = buffer.data(ksh1 + 570);
    const auto *ksh1_571 = buffer.data(ksh1 + 571);
    const auto *ksh1_572 = buffer.data(ksh1 + 572);
    const auto *ksh1_573 = buffer.data(ksh1 + 573);
    const auto *ksh1_574 = buffer.data(ksh1 + 574);
    const auto *ksh1_575 = buffer.data(ksh1 + 575);
    const auto *ksh1_576 = buffer.data(ksh1 + 576);
    const auto *ksh1_588 = buffer.data(ksh1 + 588);
    const auto *ksh1_589 = buffer.data(ksh1 + 589);
    const auto *ksh1_591 = buffer.data(ksh1 + 591);
    const auto *ksh1_593 = buffer.data(ksh1 + 593);
    const auto *ksh1_594 = buffer.data(ksh1 + 594);
    const auto *ksh1_596 = buffer.data(ksh1 + 596);
    const auto *ksh1_597 = buffer.data(ksh1 + 597);
    const auto *ksh1_598 = buffer.data(ksh1 + 598);
    const auto *ksh1_600 = buffer.data(ksh1 + 600);
    const auto *ksh1_601 = buffer.data(ksh1 + 601);
    const auto *ksh1_602 = buffer.data(ksh1 + 602);
    const auto *ksh1_603 = buffer.data(ksh1 + 603);
    const auto *ksh1_604 = buffer.data(ksh1 + 604);
    const auto *ksh1_605 = buffer.data(ksh1 + 605);
    const auto *ksh1_606 = buffer.data(ksh1 + 606);
    const auto *ksh1_607 = buffer.data(ksh1 + 607);
    const auto *ksh1_608 = buffer.data(ksh1 + 608);
    const auto *ksh1_611 = buffer.data(ksh1 + 611);
    const auto *ksh1_613 = buffer.data(ksh1 + 613);
    const auto *ksh1_614 = buffer.data(ksh1 + 614);
    const auto *ksh1_616 = buffer.data(ksh1 + 616);
    const auto *ksh1_617 = buffer.data(ksh1 + 617);
    const auto *ksh1_618 = buffer.data(ksh1 + 618);
    const auto *ksh1_620 = buffer.data(ksh1 + 620);
    const auto *ksh1_621 = buffer.data(ksh1 + 621);
    const auto *ksh1_622 = buffer.data(ksh1 + 622);
    const auto *ksh1_623 = buffer.data(ksh1 + 623);
    const auto *ksh1_625 = buffer.data(ksh1 + 625);
    const auto *ksh1_626 = buffer.data(ksh1 + 626);
    const auto *ksh1_627 = buffer.data(ksh1 + 627);
    const auto *ksh1_628 = buffer.data(ksh1 + 628);
    const auto *ksh1_629 = buffer.data(ksh1 + 629);

    const auto *ksi_737 = buffer.data(ksi + 737);
    const auto *ksi_738 = buffer.data(ksi + 738);
    const auto *ksi_742 = buffer.data(ksi + 742);
    const auto *ksi_749 = buffer.data(ksi + 749);
    const auto *ksi_750 = buffer.data(ksi + 750);
    const auto *ksi_751 = buffer.data(ksi + 751);
    const auto *ksi_752 = buffer.data(ksi + 752);
    const auto *ksi_753 = buffer.data(ksi + 753);
    const auto *ksi_754 = buffer.data(ksi + 754);
    const auto *ksi_755 = buffer.data(ksi + 755);
    const auto *ksi_756 = buffer.data(ksi + 756);
    const auto *ksi_757 = buffer.data(ksi + 757);
    const auto *ksi_758 = buffer.data(ksi + 758);
    const auto *ksi_759 = buffer.data(ksi + 759);
    const auto *ksi_760 = buffer.data(ksi + 760);
    const auto *ksi_761 = buffer.data(ksi + 761);
    const auto *ksi_762 = buffer.data(ksi + 762);
    const auto *ksi_763 = buffer.data(ksi + 763);
    const auto *ksi_764 = buffer.data(ksi + 764);
    const auto *ksi_765 = buffer.data(ksi + 765);
    const auto *ksi_766 = buffer.data(ksi + 766);
    const auto *ksi_767 = buffer.data(ksi + 767);
    const auto *ksi_768 = buffer.data(ksi + 768);
    const auto *ksi_769 = buffer.data(ksi + 769);
    const auto *ksi_770 = buffer.data(ksi + 770);
    const auto *ksi_776 = buffer.data(ksi + 776);
    const auto *ksi_777 = buffer.data(ksi + 777);
    const auto *ksi_778 = buffer.data(ksi + 778);
    const auto *ksi_779 = buffer.data(ksi + 779);
    const auto *ksi_780 = buffer.data(ksi + 780);
    const auto *ksi_781 = buffer.data(ksi + 781);
    const auto *ksi_783 = buffer.data(ksi + 783);
    const auto *ksi_784 = buffer.data(ksi + 784);
    const auto *ksi_785 = buffer.data(ksi + 785);
    const auto *ksi_787 = buffer.data(ksi + 787);
    const auto *ksi_789 = buffer.data(ksi + 789);
    const auto *ksi_790 = buffer.data(ksi + 790);
    const auto *ksi_792 = buffer.data(ksi + 792);
    const auto *ksi_793 = buffer.data(ksi + 793);
    const auto *ksi_794 = buffer.data(ksi + 794);
    const auto *ksi_796 = buffer.data(ksi + 796);
    const auto *ksi_797 = buffer.data(ksi + 797);
    const auto *ksi_798 = buffer.data(ksi + 798);
    const auto *ksi_799 = buffer.data(ksi + 799);
    const auto *ksi_801 = buffer.data(ksi + 801);
    const auto *ksi_802 = buffer.data(ksi + 802);
    const auto *ksi_803 = buffer.data(ksi + 803);
    const auto *ksi_804 = buffer.data(ksi + 804);
    const auto *ksi_805 = buffer.data(ksi + 805);
    const auto *ksi_806 = buffer.data(ksi + 806);
    const auto *ksi_807 = buffer.data(ksi + 807);
    const auto *ksi_808 = buffer.data(ksi + 808);
    const auto *ksi_809 = buffer.data(ksi + 809);
    const auto *ksi_810 = buffer.data(ksi + 810);
    const auto *ksi_811 = buffer.data(ksi + 811);
    const auto *ksi_814 = buffer.data(ksi + 814);
    const auto *ksi_816 = buffer.data(ksi + 816);
    const auto *ksi_817 = buffer.data(ksi + 817);
    const auto *ksi_819 = buffer.data(ksi + 819);
    const auto *ksi_820 = buffer.data(ksi + 820);
    const auto *ksi_821 = buffer.data(ksi + 821);
    const auto *ksi_823 = buffer.data(ksi + 823);
    const auto *ksi_824 = buffer.data(ksi + 824);
    const auto *ksi_825 = buffer.data(ksi + 825);
    const auto *ksi_826 = buffer.data(ksi + 826);
    const auto *ksi_828 = buffer.data(ksi + 828);
    const auto *ksi_829 = buffer.data(ksi + 829);
    const auto *ksi_830 = buffer.data(ksi + 830);
    const auto *ksi_831 = buffer.data(ksi + 831);
    const auto *ksi_832 = buffer.data(ksi + 832);
    const auto *ksi_833 = buffer.data(ksi + 833);
    const auto *ksi_834 = buffer.data(ksi + 834);
    const auto *ksi_835 = buffer.data(ksi + 835);
    const auto *ksi_836 = buffer.data(ksi + 836);
    const auto *ksi_837 = buffer.data(ksi + 837);
    const auto *ksi_838 = buffer.data(ksi + 838);
    const auto *ksi_839 = buffer.data(ksi + 839);

#pragma omp simd aligned(t_948, t_949, t_950, pa_x, pa_y, pc_x, pc_y, isk0_734, isk0_948, \
                         isi_569, isi_740, isk1_734, isk1_948, \
                         ksi_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = pa_x[k] * isk0_948[k]
                   + f_15 * isi_740[k]
                   - f_12 * pc_x[k] * isk1_948[k];

        t_949[k] = f_13 * isi_569[k]
                   + f_3 * pc_y[k] * ksi_737[k];

        t_950[k] = pa_y[k] * isk0_734[k]
                   - f_12 * pc_y[k] * isk1_734[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, pa_x, pc_x, pc_z, isk0_951, isk0_953, isi_542, \
                         isi_743, isi_745, isk1_951, isk1_953, \
                         ksi_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = pa_x[k] * isk0_951[k]
                   + f_14 * isi_743[k]
                   - f_12 * pc_x[k] * isk1_951[k];

        t_952[k] = f_17 * isi_542[k]
                   + f_3 * pc_z[k] * ksi_738[k];

        t_953[k] = pa_x[k] * isk0_953[k]
                   + f_14 * isi_745[k]
                   - f_12 * pc_x[k] * isk1_953[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pa_x, pa_y, pc_x, pc_y, isk0_740, isk0_954, \
                         isi_574, isi_746, isk1_740, isk1_954, \
                         ksi_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = pa_x[k] * isk0_954[k]
                   + f_14 * isi_746[k]
                   - f_12 * pc_x[k] * isk1_954[k];

        t_955[k] = f_13 * isi_574[k]
                   + f_3 * pc_y[k] * ksi_742[k];

        t_956[k] = pa_y[k] * isk0_740[k]
                   - f_12 * pc_y[k] * isk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, isi_749, isi_750, isi_751, \
                         isi_752, isi_753, ksi_749, ksi_750, ksi_751, ksi_752, \
                         ksi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_13 * isi_749[k]
                   + f_3 * pc_x[k] * ksi_749[k];

        t_958[k] = f_13 * isi_750[k]
                   + f_3 * pc_x[k] * ksi_750[k];

        t_959[k] = f_13 * isi_751[k]
                   + f_3 * pc_x[k] * ksi_751[k];

        t_960[k] = f_13 * isi_752[k]
                   + f_3 * pc_x[k] * ksi_752[k];

        t_961[k] = f_13 * isi_753[k]
                   + f_3 * pc_x[k] * ksi_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pa_x, pc_x, pc_z, isk0_964, isi_553, \
                         isi_754, isi_755, isk1_964, ksi_749, ksi_754, \
                         ksi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_13 * isi_754[k]
                   + f_3 * pc_x[k] * ksi_754[k];

        t_963[k] = f_13 * isi_755[k]
                   + f_3 * pc_x[k] * ksi_755[k];

        t_964[k] = pa_x[k] * isk0_964[k]
                   - f_12 * pc_x[k] * isk1_964[k];

        t_965[k] = f_17 * isi_553[k]
                   + f_3 * pc_z[k] * ksi_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, t_969, pa_x, pc_x, isk0_966, isk0_967, isk0_968, \
                         isk0_969, isk1_966, isk1_967, isk1_968, \
                         isk1_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = pa_x[k] * isk0_966[k]
                   - f_12 * pc_x[k] * isk1_966[k];

        t_967[k] = pa_x[k] * isk0_967[k]
                   - f_12 * pc_x[k] * isk1_967[k];

        t_968[k] = pa_x[k] * isk0_968[k]
                   - f_12 * pc_x[k] * isk1_968[k];

        t_969[k] = pa_x[k] * isk0_969[k]
                   - f_12 * pc_x[k] * isk1_969[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, pa_x, pc_x, pc_y, isk0_971, isk0_972, \
                         isi_587, isi_756, isk1_971, isk1_972, ksi_755, \
                         ksi_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_13 * isi_587[k]
                   + f_3 * pc_y[k] * ksi_755[k];

        t_971[k] = pa_x[k] * isk0_971[k]
                   - f_12 * pc_x[k] * isk1_971[k];

        t_972[k] = pa_x[k] * isk0_972[k]
                   + f_0 * isi_756[k]
                   - f_12 * pc_x[k] * isk1_972[k];

        t_973[k] = f_3 * pc_y[k] * ksi_756[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, pc_y, pc_z, isi_560, ksh0_567, ksh1_567, \
                         ksi_756, ksi_757, ksi_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_18 * isi_560[k]
                   + f_3 * pc_z[k] * ksi_756[k];

        t_975[k] = f_4 * ksh0_567[k]
                   - f_5 * ksh1_567[k]
                   + f_3 * pc_y[k] * ksi_757[k];

        t_976[k] = f_3 * pc_y[k] * ksi_758[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, pa_x, pc_x, pc_y, isk0_977, isi_761, isk1_977, \
                         ksh0_568, ksh0_569, ksh1_568, ksh1_569, ksi_759, \
                         ksi_760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = pa_x[k] * isk0_977[k]
                   + f_17 * isi_761[k]
                   - f_12 * pc_x[k] * isk1_977[k];

        t_978[k] = f_6 * ksh0_568[k]
                   - f_7 * ksh1_568[k]
                   + f_3 * pc_y[k] * ksi_759[k];

        t_979[k] = f_4 * ksh0_569[k]
                   - f_5 * ksh1_569[k]
                   + f_3 * pc_y[k] * ksi_760[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, pa_x, pc_x, pc_y, isk0_981, isi_765, isk1_981, \
                         ksh0_570, ksh1_570, ksi_761, ksi_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_3 * pc_y[k] * ksi_761[k];

        t_981[k] = pa_x[k] * isk0_981[k]
                   + f_16 * isi_765[k]
                   - f_12 * pc_x[k] * isk1_981[k];

        t_982[k] = f_8 * ksh0_570[k]
                   - f_9 * ksh1_570[k]
                   + f_3 * pc_y[k] * ksi_762[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_y, ksh0_571, ksh0_572, ksh1_571, ksh1_572, \
                         ksi_763, ksi_764, ksi_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_6 * ksh0_571[k]
                   - f_7 * ksh1_571[k]
                   + f_3 * pc_y[k] * ksi_763[k];

        t_984[k] = f_4 * ksh0_572[k]
                   - f_5 * ksh1_572[k]
                   + f_3 * pc_y[k] * ksi_764[k];

        t_985[k] = f_3 * pc_y[k] * ksi_765[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pa_x, pc_x, pc_y, isk0_986, isi_770, isk1_986, \
                         ksh0_573, ksh0_574, ksh1_573, ksh1_574, ksi_766, \
                         ksi_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = pa_x[k] * isk0_986[k]
                   + f_15 * isi_770[k]
                   - f_12 * pc_x[k] * isk1_986[k];

        t_987[k] = f_10 * ksh0_573[k]
                   - f_11 * ksh1_573[k]
                   + f_3 * pc_y[k] * ksi_766[k];

        t_988[k] = f_8 * ksh0_574[k]
                   - f_9 * ksh1_574[k]
                   + f_3 * pc_y[k] * ksi_767[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_y, ksh0_575, ksh0_576, ksh1_575, ksh1_576, \
                         ksi_768, ksi_769, ksi_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_6 * ksh0_575[k]
                   - f_7 * ksh1_575[k]
                   + f_3 * pc_y[k] * ksi_768[k];

        t_990[k] = f_4 * ksh0_576[k]
                   - f_5 * ksh1_576[k]
                   + f_3 * pc_y[k] * ksi_769[k];

        t_991[k] = f_3 * pc_y[k] * ksi_770[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pa_x, pc_x, isk0_992, isi_776, isi_777, \
                         isi_778, isi_779, isk1_992, ksi_777, ksi_778, \
                         ksi_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = pa_x[k] * isk0_992[k]
                   + f_14 * isi_776[k]
                   - f_12 * pc_x[k] * isk1_992[k];

        t_993[k] = f_13 * isi_777[k]
                   + f_3 * pc_x[k] * ksi_777[k];

        t_994[k] = f_13 * isi_778[k]
                   + f_3 * pc_x[k] * ksi_778[k];

        t_995[k] = f_13 * isi_779[k]
                   + f_3 * pc_x[k] * ksi_779[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, pc_x, pc_y, isi_780, isi_781, isi_783, \
                         ksi_776, ksi_780, ksi_781, ksi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_13 * isi_780[k]
                   + f_3 * pc_x[k] * ksi_780[k];

        t_997[k] = f_13 * isi_781[k]
                   + f_3 * pc_x[k] * ksi_781[k];

        t_998[k] = f_3 * pc_y[k] * ksi_776[k];

        t_999[k] = f_13 * isi_783[k]
                   + f_3 * pc_x[k] * ksi_783[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pa_x, pc_x, isk0_1000, isk0_1001, \
                         isk0_1002, isk0_1003, isk1_1000, isk1_1001, isk1_1002, \
                         isk1_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = pa_x[k] * isk0_1000[k]
                    - f_12 * pc_x[k] * isk1_1000[k];

        t_1001[k] = pa_x[k] * isk0_1001[k]
                    - f_12 * pc_x[k] * isk1_1001[k];

        t_1002[k] = pa_x[k] * isk0_1002[k]
                    - f_12 * pc_x[k] * isk1_1002[k];

        t_1003[k] = pa_x[k] * isk0_1003[k]
                    - f_12 * pc_x[k] * isk1_1003[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pa_x, pc_x, pc_y, isk0_1004, \
                         isk0_1005, isk0_1007, isk1_1004, isk1_1005, isk1_1007, \
                         ksi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = pa_x[k] * isk0_1004[k]
                    - f_12 * pc_x[k] * isk1_1004[k];

        t_1005[k] = pa_x[k] * isk0_1005[k]
                    - f_12 * pc_x[k] * isk1_1005[k];

        t_1006[k] = f_3 * pc_y[k] * ksi_783[k];

        t_1007[k] = pa_x[k] * isk0_1007[k]
                    - f_12 * pc_x[k] * isk1_1007[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, t_1012, pc_x, pc_z, ksh0_588, \
                         ksh0_589, ksh0_591, ksh1_588, ksh1_589, ksh1_591, ksi_784, ksi_785, \
                         ksi_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_1 * ksh0_588[k]
                    - f_2 * ksh1_588[k]
                    + f_3 * pc_x[k] * ksi_784[k];

        t_1009[k] = f_19 * ksh0_589[k]
                    - f_20 * ksh1_589[k]
                    + f_3 * pc_x[k] * ksi_785[k];

        t_1010[k] = f_3 * pc_z[k] * ksi_784[k];

        t_1011[k] = f_10 * ksh0_591[k]
                    - f_11 * ksh1_591[k]
                    + f_3 * pc_x[k] * ksi_787[k];

        t_1012[k] = f_3 * pc_z[k] * ksi_785[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, t_1016, pc_x, pc_z, ksh0_593, ksh0_594, \
                         ksh0_596, ksh1_593, ksh1_594, ksh1_596, ksi_787, ksi_789, ksi_790, \
                         ksi_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_10 * ksh0_593[k]
                    - f_11 * ksh1_593[k]
                    + f_3 * pc_x[k] * ksi_789[k];

        t_1014[k] = f_8 * ksh0_594[k]
                    - f_9 * ksh1_594[k]
                    + f_3 * pc_x[k] * ksi_790[k];

        t_1015[k] = f_3 * pc_z[k] * ksi_787[k];

        t_1016[k] = f_8 * ksh0_596[k]
                    - f_9 * ksh1_596[k]
                    + f_3 * pc_x[k] * ksi_792[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, pc_x, pc_z, ksh0_597, ksh0_598, \
                         ksh0_600, ksh1_597, ksh1_598, ksh1_600, ksi_790, ksi_793, ksi_794, \
                         ksi_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = f_8 * ksh0_597[k]
                    - f_9 * ksh1_597[k]
                    + f_3 * pc_x[k] * ksi_793[k];

        t_1018[k] = f_6 * ksh0_598[k]
                    - f_7 * ksh1_598[k]
                    + f_3 * pc_x[k] * ksi_794[k];

        t_1019[k] = f_3 * pc_z[k] * ksi_790[k];

        t_1020[k] = f_6 * ksh0_600[k]
                    - f_7 * ksh1_600[k]
                    + f_3 * pc_x[k] * ksi_796[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pc_x, pc_z, ksh0_601, ksh0_602, \
                         ksh0_603, ksh1_601, ksh1_602, ksh1_603, ksi_794, ksi_797, ksi_798, \
                         ksi_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_6 * ksh0_601[k]
                    - f_7 * ksh1_601[k]
                    + f_3 * pc_x[k] * ksi_797[k];

        t_1022[k] = f_6 * ksh0_602[k]
                    - f_7 * ksh1_602[k]
                    + f_3 * pc_x[k] * ksi_798[k];

        t_1023[k] = f_4 * ksh0_603[k]
                    - f_5 * ksh1_603[k]
                    + f_3 * pc_x[k] * ksi_799[k];

        t_1024[k] = f_3 * pc_z[k] * ksi_794[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pc_x, ksh0_605, ksh0_606, ksh0_607, ksh1_605, \
                         ksh1_606, ksh1_607, ksi_801, ksi_802, \
                         ksi_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_4 * ksh0_605[k]
                    - f_5 * ksh1_605[k]
                    + f_3 * pc_x[k] * ksi_801[k];

        t_1026[k] = f_4 * ksh0_606[k]
                    - f_5 * ksh1_606[k]
                    + f_3 * pc_x[k] * ksi_802[k];

        t_1027[k] = f_4 * ksh0_607[k]
                    - f_5 * ksh1_607[k]
                    + f_3 * pc_x[k] * ksi_803[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, t_1031, t_1032, t_1033, pc_x, ksh0_608, \
                         ksh1_608, ksi_804, ksi_805, ksi_806, ksi_807, ksi_808, \
                         ksi_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_4 * ksh0_608[k]
                    - f_5 * ksh1_608[k]
                    + f_3 * pc_x[k] * ksi_804[k];

        t_1029[k] = f_3 * pc_x[k] * ksi_805[k];

        t_1030[k] = f_3 * pc_x[k] * ksi_806[k];

        t_1031[k] = f_3 * pc_x[k] * ksi_807[k];

        t_1032[k] = f_3 * pc_x[k] * ksi_808[k];

        t_1033[k] = f_3 * pc_x[k] * ksi_809[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, t_1037, t_1038, pc_x, pc_y, pc_z, isi_609, \
                         ksh0_603, ksh1_603, ksi_805, ksi_806, ksi_810, \
                         ksi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_3 * pc_x[k] * ksi_810[k];

        t_1035[k] = f_3 * pc_x[k] * ksi_811[k];

        t_1036[k] = f_0 * isi_609[k]
                    + f_1 * ksh0_603[k]
                    - f_2 * ksh1_603[k]
                    + f_3 * pc_y[k] * ksi_805[k];

        t_1037[k] = f_3 * pc_z[k] * ksi_805[k];

        t_1038[k] = f_4 * ksh0_603[k]
                    - f_5 * ksh1_603[k]
                    + f_3 * pc_z[k] * ksi_806[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_z, ksh0_604, ksh0_605, ksh0_606, ksh1_604, \
                         ksh1_605, ksh1_606, ksi_807, ksi_808, \
                         ksi_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_6 * ksh0_604[k]
                    - f_7 * ksh1_604[k]
                    + f_3 * pc_z[k] * ksi_807[k];

        t_1040[k] = f_8 * ksh0_605[k]
                    - f_9 * ksh1_605[k]
                    + f_3 * pc_z[k] * ksi_808[k];

        t_1041[k] = f_10 * ksh0_606[k]
                    - f_11 * ksh1_606[k]
                    + f_3 * pc_z[k] * ksi_809[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pa_z, pc_y, pc_z, isk0_756, isk0_757, \
                         isi_615, isk1_756, isk1_757, ksh0_608, ksh1_608, \
                         ksi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_0 * isi_615[k]
                    + f_3 * pc_y[k] * ksi_811[k];

        t_1043[k] = f_1 * ksh0_608[k]
                    - f_2 * ksh1_608[k]
                    + f_3 * pc_z[k] * ksi_811[k];

        t_1044[k] = pa_z[k] * isk0_756[k]
                    - f_12 * pc_z[k] * isk1_756[k];

        t_1045[k] = pa_z[k] * isk0_757[k]
                    - f_12 * pc_z[k] * isk1_757[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pa_z, pc_x, pc_z, isk0_759, isk1_759, \
                         ksh0_611, ksh0_613, ksh1_611, ksh1_613, ksi_814, \
                         ksi_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_19 * ksh0_611[k]
                    - f_20 * ksh1_611[k]
                    + f_3 * pc_x[k] * ksi_814[k];

        t_1047[k] = pa_z[k] * isk0_759[k]
                    - f_12 * pc_z[k] * isk1_759[k];

        t_1048[k] = f_10 * ksh0_613[k]
                    - f_11 * ksh1_613[k]
                    + f_3 * pc_x[k] * ksi_816[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pa_z, pc_x, pc_z, isk0_762, isk1_762, \
                         ksh0_614, ksh0_616, ksh1_614, ksh1_616, ksi_817, \
                         ksi_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_10 * ksh0_614[k]
                    - f_11 * ksh1_614[k]
                    + f_3 * pc_x[k] * ksi_817[k];

        t_1050[k] = pa_z[k] * isk0_762[k]
                    - f_12 * pc_z[k] * isk1_762[k];

        t_1051[k] = f_8 * ksh0_616[k]
                    - f_9 * ksh1_616[k]
                    + f_3 * pc_x[k] * ksi_819[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pa_z, pc_x, pc_z, isk0_766, isk1_766, \
                         ksh0_617, ksh0_618, ksh1_617, ksh1_618, ksi_820, \
                         ksi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_8 * ksh0_617[k]
                    - f_9 * ksh1_617[k]
                    + f_3 * pc_x[k] * ksi_820[k];

        t_1053[k] = f_8 * ksh0_618[k]
                    - f_9 * ksh1_618[k]
                    + f_3 * pc_x[k] * ksi_821[k];

        t_1054[k] = pa_z[k] * isk0_766[k]
                    - f_12 * pc_z[k] * isk1_766[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, ksh0_620, ksh0_621, ksh0_622, ksh1_620, \
                         ksh1_621, ksh1_622, ksi_823, ksi_824, \
                         ksi_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_6 * ksh0_620[k]
                    - f_7 * ksh1_620[k]
                    + f_3 * pc_x[k] * ksi_823[k];

        t_1056[k] = f_6 * ksh0_621[k]
                    - f_7 * ksh1_621[k]
                    + f_3 * pc_x[k] * ksi_824[k];

        t_1057[k] = f_6 * ksh0_622[k]
                    - f_7 * ksh1_622[k]
                    + f_3 * pc_x[k] * ksi_825[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pa_z, pc_x, pc_z, isk0_771, isk1_771, \
                         ksh0_623, ksh0_625, ksh1_623, ksh1_625, ksi_826, \
                         ksi_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_6 * ksh0_623[k]
                    - f_7 * ksh1_623[k]
                    + f_3 * pc_x[k] * ksi_826[k];

        t_1059[k] = pa_z[k] * isk0_771[k]
                    - f_12 * pc_z[k] * isk1_771[k];

        t_1060[k] = f_4 * ksh0_625[k]
                    - f_5 * ksh1_625[k]
                    + f_3 * pc_x[k] * ksi_828[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_x, ksh0_626, ksh0_627, ksh0_628, ksh1_626, \
                         ksh1_627, ksh1_628, ksi_829, ksi_830, \
                         ksi_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_4 * ksh0_626[k]
                    - f_5 * ksh1_626[k]
                    + f_3 * pc_x[k] * ksi_829[k];

        t_1062[k] = f_4 * ksh0_627[k]
                    - f_5 * ksh1_627[k]
                    + f_3 * pc_x[k] * ksi_830[k];

        t_1063[k] = f_4 * ksh0_628[k]
                    - f_5 * ksh1_628[k]
                    + f_3 * pc_x[k] * ksi_831[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, ksh0_629, \
                         ksh1_629, ksi_832, ksi_833, ksi_834, ksi_835, ksi_836, \
                         ksi_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_4 * ksh0_629[k]
                    - f_5 * ksh1_629[k]
                    + f_3 * pc_x[k] * ksi_832[k];

        t_1065[k] = f_3 * pc_x[k] * ksi_833[k];

        t_1066[k] = f_3 * pc_x[k] * ksi_834[k];

        t_1067[k] = f_3 * pc_x[k] * ksi_835[k];

        t_1068[k] = f_3 * pc_x[k] * ksi_836[k];

        t_1069[k] = f_3 * pc_x[k] * ksi_837[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, pa_z, pc_x, pc_z, isk0_784, isi_609, \
                         isk1_784, ksi_833, ksi_838, ksi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_3 * pc_x[k] * ksi_838[k];

        t_1071[k] = f_3 * pc_x[k] * ksi_839[k];

        t_1072[k] = pa_z[k] * isk0_784[k]
                    - f_12 * pc_z[k] * isk1_784[k];

        t_1073[k] = f_13 * isi_609[k]
                    + f_3 * pc_z[k] * ksi_833[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isk0,
                                                          const size_t isi, const size_t isk1,
                                                          const size_t ksh0, const size_t ksh1,
                                                          const size_t ksi, const size_t ncols,
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isk0_786 = buffer.data(isk0 + 786);
    const auto *isk0_787 = buffer.data(isk0 + 787);
    const auto *isk0_788 = buffer.data(isk0 + 788);
    const auto *isk0_789 = buffer.data(isk0 + 789);

    const auto *isi_610 = buffer.data(isi + 610);
    const auto *isi_611 = buffer.data(isi + 611);
    const auto *isi_612 = buffer.data(isi + 612);
    const auto *isi_613 = buffer.data(isi + 613);
    const auto *isi_615 = buffer.data(isi + 615);
    const auto *isi_637 = buffer.data(isi + 637);
    const auto *isi_643 = buffer.data(isi + 643);
    const auto *isi_665 = buffer.data(isi + 665);
    const auto *isi_667 = buffer.data(isi + 667);
    const auto *isi_668 = buffer.data(isi + 668);
    const auto *isi_669 = buffer.data(isi + 669);
    const auto *isi_670 = buffer.data(isi + 670);
    const auto *isi_671 = buffer.data(isi + 671);
    const auto *isi_693 = buffer.data(isi + 693);
    const auto *isi_695 = buffer.data(isi + 695);
    const auto *isi_696 = buffer.data(isi + 696);
    const auto *isi_697 = buffer.data(isi + 697);
    const auto *isi_698 = buffer.data(isi + 698);
    const auto *isi_699 = buffer.data(isi + 699);
    const auto *isi_721 = buffer.data(isi + 721);
    const auto *isi_723 = buffer.data(isi + 723);
    const auto *isi_724 = buffer.data(isi + 724);
    const auto *isi_725 = buffer.data(isi + 725);
    const auto *isi_726 = buffer.data(isi + 726);
    const auto *isi_727 = buffer.data(isi + 727);

    const auto *isk1_786 = buffer.data(isk1 + 786);
    const auto *isk1_787 = buffer.data(isk1 + 787);
    const auto *isk1_788 = buffer.data(isk1 + 788);
    const auto *isk1_789 = buffer.data(isk1 + 789);

    const auto *ksh0_629 = buffer.data(ksh0 + 629);
    const auto *ksh0_630 = buffer.data(ksh0 + 630);
    const auto *ksh0_631 = buffer.data(ksh0 + 631);
    const auto *ksh0_632 = buffer.data(ksh0 + 632);
    const auto *ksh0_633 = buffer.data(ksh0 + 633);
    const auto *ksh0_634 = buffer.data(ksh0 + 634);
    const auto *ksh0_635 = buffer.data(ksh0 + 635);
    const auto *ksh0_636 = buffer.data(ksh0 + 636);
    const auto *ksh0_637 = buffer.data(ksh0 + 637);
    const auto *ksh0_638 = buffer.data(ksh0 + 638);
    const auto *ksh0_639 = buffer.data(ksh0 + 639);
    const auto *ksh0_640 = buffer.data(ksh0 + 640);
    const auto *ksh0_641 = buffer.data(ksh0 + 641);
    const auto *ksh0_642 = buffer.data(ksh0 + 642);
    const auto *ksh0_643 = buffer.data(ksh0 + 643);
    const auto *ksh0_644 = buffer.data(ksh0 + 644);
    const auto *ksh0_645 = buffer.data(ksh0 + 645);
    const auto *ksh0_646 = buffer.data(ksh0 + 646);
    const auto *ksh0_647 = buffer.data(ksh0 + 647);
    const auto *ksh0_648 = buffer.data(ksh0 + 648);
    const auto *ksh0_649 = buffer.data(ksh0 + 649);
    const auto *ksh0_650 = buffer.data(ksh0 + 650);
    const auto *ksh0_651 = buffer.data(ksh0 + 651);
    const auto *ksh0_652 = buffer.data(ksh0 + 652);
    const auto *ksh0_653 = buffer.data(ksh0 + 653);
    const auto *ksh0_654 = buffer.data(ksh0 + 654);
    const auto *ksh0_655 = buffer.data(ksh0 + 655);
    const auto *ksh0_656 = buffer.data(ksh0 + 656);
    const auto *ksh0_657 = buffer.data(ksh0 + 657);
    const auto *ksh0_658 = buffer.data(ksh0 + 658);
    const auto *ksh0_659 = buffer.data(ksh0 + 659);
    const auto *ksh0_660 = buffer.data(ksh0 + 660);
    const auto *ksh0_661 = buffer.data(ksh0 + 661);
    const auto *ksh0_662 = buffer.data(ksh0 + 662);
    const auto *ksh0_663 = buffer.data(ksh0 + 663);
    const auto *ksh0_664 = buffer.data(ksh0 + 664);
    const auto *ksh0_665 = buffer.data(ksh0 + 665);
    const auto *ksh0_666 = buffer.data(ksh0 + 666);
    const auto *ksh0_667 = buffer.data(ksh0 + 667);
    const auto *ksh0_668 = buffer.data(ksh0 + 668);
    const auto *ksh0_669 = buffer.data(ksh0 + 669);
    const auto *ksh0_670 = buffer.data(ksh0 + 670);
    const auto *ksh0_671 = buffer.data(ksh0 + 671);
    const auto *ksh0_672 = buffer.data(ksh0 + 672);
    const auto *ksh0_673 = buffer.data(ksh0 + 673);
    const auto *ksh0_674 = buffer.data(ksh0 + 674);
    const auto *ksh0_675 = buffer.data(ksh0 + 675);
    const auto *ksh0_676 = buffer.data(ksh0 + 676);
    const auto *ksh0_677 = buffer.data(ksh0 + 677);
    const auto *ksh0_678 = buffer.data(ksh0 + 678);
    const auto *ksh0_679 = buffer.data(ksh0 + 679);
    const auto *ksh0_680 = buffer.data(ksh0 + 680);
    const auto *ksh0_681 = buffer.data(ksh0 + 681);
    const auto *ksh0_682 = buffer.data(ksh0 + 682);
    const auto *ksh0_683 = buffer.data(ksh0 + 683);
    const auto *ksh0_684 = buffer.data(ksh0 + 684);
    const auto *ksh0_685 = buffer.data(ksh0 + 685);
    const auto *ksh0_686 = buffer.data(ksh0 + 686);
    const auto *ksh0_687 = buffer.data(ksh0 + 687);
    const auto *ksh0_688 = buffer.data(ksh0 + 688);
    const auto *ksh0_689 = buffer.data(ksh0 + 689);
    const auto *ksh0_690 = buffer.data(ksh0 + 690);
    const auto *ksh0_691 = buffer.data(ksh0 + 691);
    const auto *ksh0_692 = buffer.data(ksh0 + 692);

    const auto *ksh1_629 = buffer.data(ksh1 + 629);
    const auto *ksh1_630 = buffer.data(ksh1 + 630);
    const auto *ksh1_631 = buffer.data(ksh1 + 631);
    const auto *ksh1_632 = buffer.data(ksh1 + 632);
    const auto *ksh1_633 = buffer.data(ksh1 + 633);
    const auto *ksh1_634 = buffer.data(ksh1 + 634);
    const auto *ksh1_635 = buffer.data(ksh1 + 635);
    const auto *ksh1_636 = buffer.data(ksh1 + 636);
    const auto *ksh1_637 = buffer.data(ksh1 + 637);
    const auto *ksh1_638 = buffer.data(ksh1 + 638);
    const auto *ksh1_639 = buffer.data(ksh1 + 639);
    const auto *ksh1_640 = buffer.data(ksh1 + 640);
    const auto *ksh1_641 = buffer.data(ksh1 + 641);
    const auto *ksh1_642 = buffer.data(ksh1 + 642);
    const auto *ksh1_643 = buffer.data(ksh1 + 643);
    const auto *ksh1_644 = buffer.data(ksh1 + 644);
    const auto *ksh1_645 = buffer.data(ksh1 + 645);
    const auto *ksh1_646 = buffer.data(ksh1 + 646);
    const auto *ksh1_647 = buffer.data(ksh1 + 647);
    const auto *ksh1_648 = buffer.data(ksh1 + 648);
    const auto *ksh1_649 = buffer.data(ksh1 + 649);
    const auto *ksh1_650 = buffer.data(ksh1 + 650);
    const auto *ksh1_651 = buffer.data(ksh1 + 651);
    const auto *ksh1_652 = buffer.data(ksh1 + 652);
    const auto *ksh1_653 = buffer.data(ksh1 + 653);
    const auto *ksh1_654 = buffer.data(ksh1 + 654);
    const auto *ksh1_655 = buffer.data(ksh1 + 655);
    const auto *ksh1_656 = buffer.data(ksh1 + 656);
    const auto *ksh1_657 = buffer.data(ksh1 + 657);
    const auto *ksh1_658 = buffer.data(ksh1 + 658);
    const auto *ksh1_659 = buffer.data(ksh1 + 659);
    const auto *ksh1_660 = buffer.data(ksh1 + 660);
    const auto *ksh1_661 = buffer.data(ksh1 + 661);
    const auto *ksh1_662 = buffer.data(ksh1 + 662);
    const auto *ksh1_663 = buffer.data(ksh1 + 663);
    const auto *ksh1_664 = buffer.data(ksh1 + 664);
    const auto *ksh1_665 = buffer.data(ksh1 + 665);
    const auto *ksh1_666 = buffer.data(ksh1 + 666);
    const auto *ksh1_667 = buffer.data(ksh1 + 667);
    const auto *ksh1_668 = buffer.data(ksh1 + 668);
    const auto *ksh1_669 = buffer.data(ksh1 + 669);
    const auto *ksh1_670 = buffer.data(ksh1 + 670);
    const auto *ksh1_671 = buffer.data(ksh1 + 671);
    const auto *ksh1_672 = buffer.data(ksh1 + 672);
    const auto *ksh1_673 = buffer.data(ksh1 + 673);
    const auto *ksh1_674 = buffer.data(ksh1 + 674);
    const auto *ksh1_675 = buffer.data(ksh1 + 675);
    const auto *ksh1_676 = buffer.data(ksh1 + 676);
    const auto *ksh1_677 = buffer.data(ksh1 + 677);
    const auto *ksh1_678 = buffer.data(ksh1 + 678);
    const auto *ksh1_679 = buffer.data(ksh1 + 679);
    const auto *ksh1_680 = buffer.data(ksh1 + 680);
    const auto *ksh1_681 = buffer.data(ksh1 + 681);
    const auto *ksh1_682 = buffer.data(ksh1 + 682);
    const auto *ksh1_683 = buffer.data(ksh1 + 683);
    const auto *ksh1_684 = buffer.data(ksh1 + 684);
    const auto *ksh1_685 = buffer.data(ksh1 + 685);
    const auto *ksh1_686 = buffer.data(ksh1 + 686);
    const auto *ksh1_687 = buffer.data(ksh1 + 687);
    const auto *ksh1_688 = buffer.data(ksh1 + 688);
    const auto *ksh1_689 = buffer.data(ksh1 + 689);
    const auto *ksh1_690 = buffer.data(ksh1 + 690);
    const auto *ksh1_691 = buffer.data(ksh1 + 691);
    const auto *ksh1_692 = buffer.data(ksh1 + 692);

    const auto *ksi_839 = buffer.data(ksi + 839);
    const auto *ksi_840 = buffer.data(ksi + 840);
    const auto *ksi_841 = buffer.data(ksi + 841);
    const auto *ksi_842 = buffer.data(ksi + 842);
    const auto *ksi_843 = buffer.data(ksi + 843);
    const auto *ksi_844 = buffer.data(ksi + 844);
    const auto *ksi_845 = buffer.data(ksi + 845);
    const auto *ksi_846 = buffer.data(ksi + 846);
    const auto *ksi_847 = buffer.data(ksi + 847);
    const auto *ksi_848 = buffer.data(ksi + 848);
    const auto *ksi_849 = buffer.data(ksi + 849);
    const auto *ksi_850 = buffer.data(ksi + 850);
    const auto *ksi_851 = buffer.data(ksi + 851);
    const auto *ksi_852 = buffer.data(ksi + 852);
    const auto *ksi_853 = buffer.data(ksi + 853);
    const auto *ksi_854 = buffer.data(ksi + 854);
    const auto *ksi_855 = buffer.data(ksi + 855);
    const auto *ksi_856 = buffer.data(ksi + 856);
    const auto *ksi_857 = buffer.data(ksi + 857);
    const auto *ksi_858 = buffer.data(ksi + 858);
    const auto *ksi_859 = buffer.data(ksi + 859);
    const auto *ksi_860 = buffer.data(ksi + 860);
    const auto *ksi_861 = buffer.data(ksi + 861);
    const auto *ksi_862 = buffer.data(ksi + 862);
    const auto *ksi_863 = buffer.data(ksi + 863);
    const auto *ksi_864 = buffer.data(ksi + 864);
    const auto *ksi_865 = buffer.data(ksi + 865);
    const auto *ksi_866 = buffer.data(ksi + 866);
    const auto *ksi_867 = buffer.data(ksi + 867);
    const auto *ksi_868 = buffer.data(ksi + 868);
    const auto *ksi_869 = buffer.data(ksi + 869);
    const auto *ksi_870 = buffer.data(ksi + 870);
    const auto *ksi_871 = buffer.data(ksi + 871);
    const auto *ksi_872 = buffer.data(ksi + 872);
    const auto *ksi_873 = buffer.data(ksi + 873);
    const auto *ksi_874 = buffer.data(ksi + 874);
    const auto *ksi_875 = buffer.data(ksi + 875);
    const auto *ksi_876 = buffer.data(ksi + 876);
    const auto *ksi_877 = buffer.data(ksi + 877);
    const auto *ksi_878 = buffer.data(ksi + 878);
    const auto *ksi_879 = buffer.data(ksi + 879);
    const auto *ksi_880 = buffer.data(ksi + 880);
    const auto *ksi_881 = buffer.data(ksi + 881);
    const auto *ksi_882 = buffer.data(ksi + 882);
    const auto *ksi_883 = buffer.data(ksi + 883);
    const auto *ksi_884 = buffer.data(ksi + 884);
    const auto *ksi_885 = buffer.data(ksi + 885);
    const auto *ksi_886 = buffer.data(ksi + 886);
    const auto *ksi_887 = buffer.data(ksi + 887);
    const auto *ksi_888 = buffer.data(ksi + 888);
    const auto *ksi_889 = buffer.data(ksi + 889);
    const auto *ksi_890 = buffer.data(ksi + 890);
    const auto *ksi_891 = buffer.data(ksi + 891);
    const auto *ksi_892 = buffer.data(ksi + 892);
    const auto *ksi_893 = buffer.data(ksi + 893);
    const auto *ksi_894 = buffer.data(ksi + 894);
    const auto *ksi_895 = buffer.data(ksi + 895);
    const auto *ksi_896 = buffer.data(ksi + 896);
    const auto *ksi_897 = buffer.data(ksi + 897);
    const auto *ksi_898 = buffer.data(ksi + 898);
    const auto *ksi_899 = buffer.data(ksi + 899);
    const auto *ksi_900 = buffer.data(ksi + 900);
    const auto *ksi_901 = buffer.data(ksi + 901);
    const auto *ksi_902 = buffer.data(ksi + 902);
    const auto *ksi_903 = buffer.data(ksi + 903);
    const auto *ksi_904 = buffer.data(ksi + 904);
    const auto *ksi_905 = buffer.data(ksi + 905);
    const auto *ksi_906 = buffer.data(ksi + 906);
    const auto *ksi_907 = buffer.data(ksi + 907);
    const auto *ksi_908 = buffer.data(ksi + 908);
    const auto *ksi_909 = buffer.data(ksi + 909);
    const auto *ksi_910 = buffer.data(ksi + 910);
    const auto *ksi_911 = buffer.data(ksi + 911);
    const auto *ksi_912 = buffer.data(ksi + 912);
    const auto *ksi_913 = buffer.data(ksi + 913);
    const auto *ksi_914 = buffer.data(ksi + 914);
    const auto *ksi_915 = buffer.data(ksi + 915);
    const auto *ksi_916 = buffer.data(ksi + 916);
    const auto *ksi_917 = buffer.data(ksi + 917);
    const auto *ksi_918 = buffer.data(ksi + 918);
    const auto *ksi_919 = buffer.data(ksi + 919);
    const auto *ksi_920 = buffer.data(ksi + 920);
    const auto *ksi_921 = buffer.data(ksi + 921);
    const auto *ksi_922 = buffer.data(ksi + 922);
    const auto *ksi_923 = buffer.data(ksi + 923);

#pragma omp simd aligned(t_1074, t_1075, t_1076, pa_z, pc_z, isk0_786, isk0_787, isk0_788, \
                         isi_610, isi_611, isi_612, isk1_786, isk1_787, \
                         isk1_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = pa_z[k] * isk0_786[k]
                    + f_14 * isi_610[k]
                    - f_12 * pc_z[k] * isk1_786[k];

        t_1075[k] = pa_z[k] * isk0_787[k]
                    + f_15 * isi_611[k]
                    - f_12 * pc_z[k] * isk1_787[k];

        t_1076[k] = pa_z[k] * isk0_788[k]
                    + f_16 * isi_612[k]
                    - f_12 * pc_z[k] * isk1_788[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, pa_z, pc_y, pc_z, isk0_789, isi_613, isi_615, \
                         isi_643, isk1_789, ksh0_629, ksh1_629, \
                         ksi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = pa_z[k] * isk0_789[k]
                    + f_17 * isi_613[k]
                    - f_12 * pc_z[k] * isk1_789[k];

        t_1078[k] = f_18 * isi_643[k]
                    + f_3 * pc_y[k] * ksi_839[k];

        t_1079[k] = f_13 * isi_615[k]
                    + f_1 * ksh0_629[k]
                    - f_2 * ksh1_629[k]
                    + f_3 * pc_z[k] * ksi_839[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, pc_x, ksh0_630, ksh0_631, ksh0_632, ksh1_630, \
                         ksh1_631, ksh1_632, ksi_840, ksi_841, \
                         ksi_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_1 * ksh0_630[k]
                    - f_2 * ksh1_630[k]
                    + f_3 * pc_x[k] * ksi_840[k];

        t_1081[k] = f_19 * ksh0_631[k]
                    - f_20 * ksh1_631[k]
                    + f_3 * pc_x[k] * ksi_841[k];

        t_1082[k] = f_19 * ksh0_632[k]
                    - f_20 * ksh1_632[k]
                    + f_3 * pc_x[k] * ksi_842[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pc_x, ksh0_633, ksh0_634, ksh0_635, ksh1_633, \
                         ksh1_634, ksh1_635, ksi_843, ksi_844, \
                         ksi_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_10 * ksh0_633[k]
                    - f_11 * ksh1_633[k]
                    + f_3 * pc_x[k] * ksi_843[k];

        t_1084[k] = f_10 * ksh0_634[k]
                    - f_11 * ksh1_634[k]
                    + f_3 * pc_x[k] * ksi_844[k];

        t_1085[k] = f_10 * ksh0_635[k]
                    - f_11 * ksh1_635[k]
                    + f_3 * pc_x[k] * ksi_845[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pc_x, ksh0_636, ksh0_637, ksh0_638, ksh1_636, \
                         ksh1_637, ksh1_638, ksi_846, ksi_847, \
                         ksi_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_8 * ksh0_636[k]
                    - f_9 * ksh1_636[k]
                    + f_3 * pc_x[k] * ksi_846[k];

        t_1087[k] = f_8 * ksh0_637[k]
                    - f_9 * ksh1_637[k]
                    + f_3 * pc_x[k] * ksi_847[k];

        t_1088[k] = f_8 * ksh0_638[k]
                    - f_9 * ksh1_638[k]
                    + f_3 * pc_x[k] * ksi_848[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pc_x, ksh0_639, ksh0_640, ksh0_641, ksh1_639, \
                         ksh1_640, ksh1_641, ksi_849, ksi_850, \
                         ksi_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_8 * ksh0_639[k]
                    - f_9 * ksh1_639[k]
                    + f_3 * pc_x[k] * ksi_849[k];

        t_1090[k] = f_6 * ksh0_640[k]
                    - f_7 * ksh1_640[k]
                    + f_3 * pc_x[k] * ksi_850[k];

        t_1091[k] = f_6 * ksh0_641[k]
                    - f_7 * ksh1_641[k]
                    + f_3 * pc_x[k] * ksi_851[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pc_x, ksh0_642, ksh0_643, ksh0_644, ksh1_642, \
                         ksh1_643, ksh1_644, ksi_852, ksi_853, \
                         ksi_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_6 * ksh0_642[k]
                    - f_7 * ksh1_642[k]
                    + f_3 * pc_x[k] * ksi_852[k];

        t_1093[k] = f_6 * ksh0_643[k]
                    - f_7 * ksh1_643[k]
                    + f_3 * pc_x[k] * ksi_853[k];

        t_1094[k] = f_6 * ksh0_644[k]
                    - f_7 * ksh1_644[k]
                    + f_3 * pc_x[k] * ksi_854[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pc_x, ksh0_645, ksh0_646, ksh0_647, ksh1_645, \
                         ksh1_646, ksh1_647, ksi_855, ksi_856, \
                         ksi_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_4 * ksh0_645[k]
                    - f_5 * ksh1_645[k]
                    + f_3 * pc_x[k] * ksi_855[k];

        t_1096[k] = f_4 * ksh0_646[k]
                    - f_5 * ksh1_646[k]
                    + f_3 * pc_x[k] * ksi_856[k];

        t_1097[k] = f_4 * ksh0_647[k]
                    - f_5 * ksh1_647[k]
                    + f_3 * pc_x[k] * ksi_857[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, t_1101, pc_x, ksh0_648, ksh0_649, ksh0_650, \
                         ksh1_648, ksh1_649, ksh1_650, ksi_858, ksi_859, ksi_860, \
                         ksi_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_4 * ksh0_648[k]
                    - f_5 * ksh1_648[k]
                    + f_3 * pc_x[k] * ksi_858[k];

        t_1099[k] = f_4 * ksh0_649[k]
                    - f_5 * ksh1_649[k]
                    + f_3 * pc_x[k] * ksi_859[k];

        t_1100[k] = f_4 * ksh0_650[k]
                    - f_5 * ksh1_650[k]
                    + f_3 * pc_x[k] * ksi_860[k];

        t_1101[k] = f_3 * pc_x[k] * ksi_861[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, t_1105, t_1106, t_1107, pc_x, ksi_862, \
                         ksi_863, ksi_864, ksi_865, ksi_866, ksi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = f_3 * pc_x[k] * ksi_862[k];

        t_1103[k] = f_3 * pc_x[k] * ksi_863[k];

        t_1104[k] = f_3 * pc_x[k] * ksi_864[k];

        t_1105[k] = f_3 * pc_x[k] * ksi_865[k];

        t_1106[k] = f_3 * pc_x[k] * ksi_866[k];

        t_1107[k] = f_3 * pc_x[k] * ksi_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, isi_637, isi_665, isi_667, \
                         ksh0_645, ksh0_647, ksh1_645, ksh1_647, ksi_861, \
                         ksi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * isi_665[k]
                    + f_1 * ksh0_645[k]
                    - f_2 * ksh1_645[k]
                    + f_3 * pc_y[k] * ksi_861[k];

        t_1109[k] = f_14 * isi_637[k]
                    + f_3 * pc_z[k] * ksi_861[k];

        t_1110[k] = f_17 * isi_667[k]
                    + f_10 * ksh0_647[k]
                    - f_11 * ksh1_647[k]
                    + f_3 * pc_y[k] * ksi_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pc_y, isi_668, isi_669, isi_670, ksh0_648, \
                         ksh0_649, ksh0_650, ksh1_648, ksh1_649, ksh1_650, ksi_864, ksi_865, \
                         ksi_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_17 * isi_668[k]
                    + f_8 * ksh0_648[k]
                    - f_9 * ksh1_648[k]
                    + f_3 * pc_y[k] * ksi_864[k];

        t_1112[k] = f_17 * isi_669[k]
                    + f_6 * ksh0_649[k]
                    - f_7 * ksh1_649[k]
                    + f_3 * pc_y[k] * ksi_865[k];

        t_1113[k] = f_17 * isi_670[k]
                    + f_4 * ksh0_650[k]
                    - f_5 * ksh1_650[k]
                    + f_3 * pc_y[k] * ksi_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_x, pc_y, pc_z, isi_643, isi_671, ksh0_650, \
                         ksh0_651, ksh1_650, ksh1_651, ksi_867, \
                         ksi_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * isi_671[k]
                    + f_3 * pc_y[k] * ksi_867[k];

        t_1115[k] = f_14 * isi_643[k]
                    + f_1 * ksh0_650[k]
                    - f_2 * ksh1_650[k]
                    + f_3 * pc_z[k] * ksi_867[k];

        t_1116[k] = f_1 * ksh0_651[k]
                    - f_2 * ksh1_651[k]
                    + f_3 * pc_x[k] * ksi_868[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, pc_x, ksh0_652, ksh0_653, ksh0_654, ksh1_652, \
                         ksh1_653, ksh1_654, ksi_869, ksi_870, \
                         ksi_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_19 * ksh0_652[k]
                    - f_20 * ksh1_652[k]
                    + f_3 * pc_x[k] * ksi_869[k];

        t_1118[k] = f_19 * ksh0_653[k]
                    - f_20 * ksh1_653[k]
                    + f_3 * pc_x[k] * ksi_870[k];

        t_1119[k] = f_10 * ksh0_654[k]
                    - f_11 * ksh1_654[k]
                    + f_3 * pc_x[k] * ksi_871[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, ksh0_655, ksh0_656, ksh0_657, ksh1_655, \
                         ksh1_656, ksh1_657, ksi_872, ksi_873, \
                         ksi_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_10 * ksh0_655[k]
                    - f_11 * ksh1_655[k]
                    + f_3 * pc_x[k] * ksi_872[k];

        t_1121[k] = f_10 * ksh0_656[k]
                    - f_11 * ksh1_656[k]
                    + f_3 * pc_x[k] * ksi_873[k];

        t_1122[k] = f_8 * ksh0_657[k]
                    - f_9 * ksh1_657[k]
                    + f_3 * pc_x[k] * ksi_874[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, ksh0_658, ksh0_659, ksh0_660, ksh1_658, \
                         ksh1_659, ksh1_660, ksi_875, ksi_876, \
                         ksi_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_8 * ksh0_658[k]
                    - f_9 * ksh1_658[k]
                    + f_3 * pc_x[k] * ksi_875[k];

        t_1124[k] = f_8 * ksh0_659[k]
                    - f_9 * ksh1_659[k]
                    + f_3 * pc_x[k] * ksi_876[k];

        t_1125[k] = f_8 * ksh0_660[k]
                    - f_9 * ksh1_660[k]
                    + f_3 * pc_x[k] * ksi_877[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, ksh0_661, ksh0_662, ksh0_663, ksh1_661, \
                         ksh1_662, ksh1_663, ksi_878, ksi_879, \
                         ksi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_6 * ksh0_661[k]
                    - f_7 * ksh1_661[k]
                    + f_3 * pc_x[k] * ksi_878[k];

        t_1127[k] = f_6 * ksh0_662[k]
                    - f_7 * ksh1_662[k]
                    + f_3 * pc_x[k] * ksi_879[k];

        t_1128[k] = f_6 * ksh0_663[k]
                    - f_7 * ksh1_663[k]
                    + f_3 * pc_x[k] * ksi_880[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, ksh0_664, ksh0_665, ksh0_666, ksh1_664, \
                         ksh1_665, ksh1_666, ksi_881, ksi_882, \
                         ksi_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_6 * ksh0_664[k]
                    - f_7 * ksh1_664[k]
                    + f_3 * pc_x[k] * ksi_881[k];

        t_1130[k] = f_6 * ksh0_665[k]
                    - f_7 * ksh1_665[k]
                    + f_3 * pc_x[k] * ksi_882[k];

        t_1131[k] = f_4 * ksh0_666[k]
                    - f_5 * ksh1_666[k]
                    + f_3 * pc_x[k] * ksi_883[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, pc_x, ksh0_667, ksh0_668, ksh0_669, ksh1_667, \
                         ksh1_668, ksh1_669, ksi_884, ksi_885, \
                         ksi_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_4 * ksh0_667[k]
                    - f_5 * ksh1_667[k]
                    + f_3 * pc_x[k] * ksi_884[k];

        t_1133[k] = f_4 * ksh0_668[k]
                    - f_5 * ksh1_668[k]
                    + f_3 * pc_x[k] * ksi_885[k];

        t_1134[k] = f_4 * ksh0_669[k]
                    - f_5 * ksh1_669[k]
                    + f_3 * pc_x[k] * ksi_886[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, ksh0_670, ksh0_671, \
                         ksh1_670, ksh1_671, ksi_887, ksi_888, ksi_889, ksi_890, \
                         ksi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_4 * ksh0_670[k]
                    - f_5 * ksh1_670[k]
                    + f_3 * pc_x[k] * ksi_887[k];

        t_1136[k] = f_4 * ksh0_671[k]
                    - f_5 * ksh1_671[k]
                    + f_3 * pc_x[k] * ksi_888[k];

        t_1137[k] = f_3 * pc_x[k] * ksi_889[k];

        t_1138[k] = f_3 * pc_x[k] * ksi_890[k];

        t_1139[k] = f_3 * pc_x[k] * ksi_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, t_1144, pc_x, pc_y, isi_693, \
                         ksh0_666, ksh1_666, ksi_889, ksi_892, ksi_893, ksi_894, \
                         ksi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_3 * pc_x[k] * ksi_892[k];

        t_1141[k] = f_3 * pc_x[k] * ksi_893[k];

        t_1142[k] = f_3 * pc_x[k] * ksi_894[k];

        t_1143[k] = f_3 * pc_x[k] * ksi_895[k];

        t_1144[k] = f_16 * isi_693[k]
                    + f_1 * ksh0_666[k]
                    - f_2 * ksh1_666[k]
                    + f_3 * pc_y[k] * ksi_889[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pc_y, pc_z, isi_665, isi_695, isi_696, \
                         ksh0_668, ksh0_669, ksh1_668, ksh1_669, ksi_889, ksi_891, \
                         ksi_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_15 * isi_665[k]
                    + f_3 * pc_z[k] * ksi_889[k];

        t_1146[k] = f_16 * isi_695[k]
                    + f_10 * ksh0_668[k]
                    - f_11 * ksh1_668[k]
                    + f_3 * pc_y[k] * ksi_891[k];

        t_1147[k] = f_16 * isi_696[k]
                    + f_8 * ksh0_669[k]
                    - f_9 * ksh1_669[k]
                    + f_3 * pc_y[k] * ksi_892[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pc_y, isi_697, isi_698, isi_699, ksh0_670, \
                         ksh0_671, ksh1_670, ksh1_671, ksi_893, ksi_894, \
                         ksi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_16 * isi_697[k]
                    + f_6 * ksh0_670[k]
                    - f_7 * ksh1_670[k]
                    + f_3 * pc_y[k] * ksi_893[k];

        t_1149[k] = f_16 * isi_698[k]
                    + f_4 * ksh0_671[k]
                    - f_5 * ksh1_671[k]
                    + f_3 * pc_y[k] * ksi_894[k];

        t_1150[k] = f_16 * isi_699[k]
                    + f_3 * pc_y[k] * ksi_895[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, pc_x, pc_z, isi_671, ksh0_671, ksh0_672, \
                         ksh0_673, ksh1_671, ksh1_672, ksh1_673, ksi_895, ksi_896, \
                         ksi_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_15 * isi_671[k]
                    + f_1 * ksh0_671[k]
                    - f_2 * ksh1_671[k]
                    + f_3 * pc_z[k] * ksi_895[k];

        t_1152[k] = f_1 * ksh0_672[k]
                    - f_2 * ksh1_672[k]
                    + f_3 * pc_x[k] * ksi_896[k];

        t_1153[k] = f_19 * ksh0_673[k]
                    - f_20 * ksh1_673[k]
                    + f_3 * pc_x[k] * ksi_897[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, pc_x, ksh0_674, ksh0_675, ksh0_676, ksh1_674, \
                         ksh1_675, ksh1_676, ksi_898, ksi_899, \
                         ksi_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = f_19 * ksh0_674[k]
                    - f_20 * ksh1_674[k]
                    + f_3 * pc_x[k] * ksi_898[k];

        t_1155[k] = f_10 * ksh0_675[k]
                    - f_11 * ksh1_675[k]
                    + f_3 * pc_x[k] * ksi_899[k];

        t_1156[k] = f_10 * ksh0_676[k]
                    - f_11 * ksh1_676[k]
                    + f_3 * pc_x[k] * ksi_900[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, ksh0_677, ksh0_678, ksh0_679, ksh1_677, \
                         ksh1_678, ksh1_679, ksi_901, ksi_902, \
                         ksi_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_10 * ksh0_677[k]
                    - f_11 * ksh1_677[k]
                    + f_3 * pc_x[k] * ksi_901[k];

        t_1158[k] = f_8 * ksh0_678[k]
                    - f_9 * ksh1_678[k]
                    + f_3 * pc_x[k] * ksi_902[k];

        t_1159[k] = f_8 * ksh0_679[k]
                    - f_9 * ksh1_679[k]
                    + f_3 * pc_x[k] * ksi_903[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, ksh0_680, ksh0_681, ksh0_682, ksh1_680, \
                         ksh1_681, ksh1_682, ksi_904, ksi_905, \
                         ksi_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_8 * ksh0_680[k]
                    - f_9 * ksh1_680[k]
                    + f_3 * pc_x[k] * ksi_904[k];

        t_1161[k] = f_8 * ksh0_681[k]
                    - f_9 * ksh1_681[k]
                    + f_3 * pc_x[k] * ksi_905[k];

        t_1162[k] = f_6 * ksh0_682[k]
                    - f_7 * ksh1_682[k]
                    + f_3 * pc_x[k] * ksi_906[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_x, ksh0_683, ksh0_684, ksh0_685, ksh1_683, \
                         ksh1_684, ksh1_685, ksi_907, ksi_908, \
                         ksi_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_6 * ksh0_683[k]
                    - f_7 * ksh1_683[k]
                    + f_3 * pc_x[k] * ksi_907[k];

        t_1164[k] = f_6 * ksh0_684[k]
                    - f_7 * ksh1_684[k]
                    + f_3 * pc_x[k] * ksi_908[k];

        t_1165[k] = f_6 * ksh0_685[k]
                    - f_7 * ksh1_685[k]
                    + f_3 * pc_x[k] * ksi_909[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_x, ksh0_686, ksh0_687, ksh0_688, ksh1_686, \
                         ksh1_687, ksh1_688, ksi_910, ksi_911, \
                         ksi_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_6 * ksh0_686[k]
                    - f_7 * ksh1_686[k]
                    + f_3 * pc_x[k] * ksi_910[k];

        t_1167[k] = f_4 * ksh0_687[k]
                    - f_5 * ksh1_687[k]
                    + f_3 * pc_x[k] * ksi_911[k];

        t_1168[k] = f_4 * ksh0_688[k]
                    - f_5 * ksh1_688[k]
                    + f_3 * pc_x[k] * ksi_912[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_x, ksh0_689, ksh0_690, ksh0_691, ksh1_689, \
                         ksh1_690, ksh1_691, ksi_913, ksi_914, \
                         ksi_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_4 * ksh0_689[k]
                    - f_5 * ksh1_689[k]
                    + f_3 * pc_x[k] * ksi_913[k];

        t_1170[k] = f_4 * ksh0_690[k]
                    - f_5 * ksh1_690[k]
                    + f_3 * pc_x[k] * ksi_914[k];

        t_1171[k] = f_4 * ksh0_691[k]
                    - f_5 * ksh1_691[k]
                    + f_3 * pc_x[k] * ksi_915[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, t_1176, t_1177, pc_x, ksh0_692, \
                         ksh1_692, ksi_916, ksi_917, ksi_918, ksi_919, ksi_920, \
                         ksi_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_4 * ksh0_692[k]
                    - f_5 * ksh1_692[k]
                    + f_3 * pc_x[k] * ksi_916[k];

        t_1173[k] = f_3 * pc_x[k] * ksi_917[k];

        t_1174[k] = f_3 * pc_x[k] * ksi_918[k];

        t_1175[k] = f_3 * pc_x[k] * ksi_919[k];

        t_1176[k] = f_3 * pc_x[k] * ksi_920[k];

        t_1177[k] = f_3 * pc_x[k] * ksi_921[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, t_1181, pc_x, pc_y, pc_z, isi_693, isi_721, \
                         ksh0_687, ksh1_687, ksi_917, ksi_922, \
                         ksi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_3 * pc_x[k] * ksi_922[k];

        t_1179[k] = f_3 * pc_x[k] * ksi_923[k];

        t_1180[k] = f_15 * isi_721[k]
                    + f_1 * ksh0_687[k]
                    - f_2 * ksh1_687[k]
                    + f_3 * pc_y[k] * ksi_917[k];

        t_1181[k] = f_16 * isi_693[k]
                    + f_3 * pc_z[k] * ksi_917[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_y, isi_723, isi_724, isi_725, ksh0_689, \
                         ksh0_690, ksh0_691, ksh1_689, ksh1_690, ksh1_691, ksi_919, ksi_920, \
                         ksi_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_15 * isi_723[k]
                    + f_10 * ksh0_689[k]
                    - f_11 * ksh1_689[k]
                    + f_3 * pc_y[k] * ksi_919[k];

        t_1183[k] = f_15 * isi_724[k]
                    + f_8 * ksh0_690[k]
                    - f_9 * ksh1_690[k]
                    + f_3 * pc_y[k] * ksi_920[k];

        t_1184[k] = f_15 * isi_725[k]
                    + f_6 * ksh0_691[k]
                    - f_7 * ksh1_691[k]
                    + f_3 * pc_y[k] * ksi_921[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_y, pc_z, isi_699, isi_726, isi_727, \
                         ksh0_692, ksh1_692, ksi_922, ksi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_15 * isi_726[k]
                    + f_4 * ksh0_692[k]
                    - f_5 * ksh1_692[k]
                    + f_3 * pc_y[k] * ksi_922[k];

        t_1186[k] = f_15 * isi_727[k]
                    + f_3 * pc_y[k] * ksi_923[k];

        t_1187[k] = f_16 * isi_699[k]
                    + f_1 * ksh0_692[k]
                    - f_2 * ksh1_692[k]
                    + f_3 * pc_z[k] * ksi_923[k];
    }
}

static auto
compute_prim_ksk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t isk0,
                                                           const size_t isi, const size_t isk1,
                                                           const size_t ksh0, const size_t ksh1,
                                                           const size_t ksi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isk0_972 = buffer.data(isk0 + 972);
    const auto *isk0_974 = buffer.data(isk0 + 974);
    const auto *isk0_977 = buffer.data(isk0 + 977);
    const auto *isk0_981 = buffer.data(isk0 + 981);
    const auto *isk0_986 = buffer.data(isk0 + 986);
    const auto *isk0_992 = buffer.data(isk0 + 992);
    const auto *isk0_1000 = buffer.data(isk0 + 1000);
    const auto *isk0_1002 = buffer.data(isk0 + 1002);
    const auto *isk0_1003 = buffer.data(isk0 + 1003);
    const auto *isk0_1004 = buffer.data(isk0 + 1004);
    const auto *isk0_1005 = buffer.data(isk0 + 1005);
    const auto *isk0_1007 = buffer.data(isk0 + 1007);

    const auto *isi_721 = buffer.data(isi + 721);
    const auto *isi_727 = buffer.data(isi + 727);
    const auto *isi_749 = buffer.data(isi + 749);
    const auto *isi_751 = buffer.data(isi + 751);
    const auto *isi_752 = buffer.data(isi + 752);
    const auto *isi_753 = buffer.data(isi + 753);
    const auto *isi_754 = buffer.data(isi + 754);
    const auto *isi_755 = buffer.data(isi + 755);
    const auto *isi_777 = buffer.data(isi + 777);
    const auto *isi_779 = buffer.data(isi + 779);
    const auto *isi_780 = buffer.data(isi + 780);
    const auto *isi_781 = buffer.data(isi + 781);
    const auto *isi_782 = buffer.data(isi + 782);
    const auto *isi_783 = buffer.data(isi + 783);

    const auto *isk1_972 = buffer.data(isk1 + 972);
    const auto *isk1_974 = buffer.data(isk1 + 974);
    const auto *isk1_977 = buffer.data(isk1 + 977);
    const auto *isk1_981 = buffer.data(isk1 + 981);
    const auto *isk1_986 = buffer.data(isk1 + 986);
    const auto *isk1_992 = buffer.data(isk1 + 992);
    const auto *isk1_1000 = buffer.data(isk1 + 1000);
    const auto *isk1_1002 = buffer.data(isk1 + 1002);
    const auto *isk1_1003 = buffer.data(isk1 + 1003);
    const auto *isk1_1004 = buffer.data(isk1 + 1004);
    const auto *isk1_1005 = buffer.data(isk1 + 1005);
    const auto *isk1_1007 = buffer.data(isk1 + 1007);

    const auto *ksh0_693 = buffer.data(ksh0 + 693);
    const auto *ksh0_694 = buffer.data(ksh0 + 694);
    const auto *ksh0_695 = buffer.data(ksh0 + 695);
    const auto *ksh0_696 = buffer.data(ksh0 + 696);
    const auto *ksh0_697 = buffer.data(ksh0 + 697);
    const auto *ksh0_698 = buffer.data(ksh0 + 698);
    const auto *ksh0_699 = buffer.data(ksh0 + 699);
    const auto *ksh0_700 = buffer.data(ksh0 + 700);
    const auto *ksh0_701 = buffer.data(ksh0 + 701);
    const auto *ksh0_702 = buffer.data(ksh0 + 702);
    const auto *ksh0_703 = buffer.data(ksh0 + 703);
    const auto *ksh0_704 = buffer.data(ksh0 + 704);
    const auto *ksh0_705 = buffer.data(ksh0 + 705);
    const auto *ksh0_706 = buffer.data(ksh0 + 706);
    const auto *ksh0_707 = buffer.data(ksh0 + 707);
    const auto *ksh0_708 = buffer.data(ksh0 + 708);
    const auto *ksh0_709 = buffer.data(ksh0 + 709);
    const auto *ksh0_710 = buffer.data(ksh0 + 710);
    const auto *ksh0_711 = buffer.data(ksh0 + 711);
    const auto *ksh0_712 = buffer.data(ksh0 + 712);
    const auto *ksh0_713 = buffer.data(ksh0 + 713);
    const auto *ksh0_715 = buffer.data(ksh0 + 715);
    const auto *ksh0_717 = buffer.data(ksh0 + 717);
    const auto *ksh0_718 = buffer.data(ksh0 + 718);
    const auto *ksh0_720 = buffer.data(ksh0 + 720);
    const auto *ksh0_721 = buffer.data(ksh0 + 721);
    const auto *ksh0_722 = buffer.data(ksh0 + 722);
    const auto *ksh0_724 = buffer.data(ksh0 + 724);
    const auto *ksh0_725 = buffer.data(ksh0 + 725);
    const auto *ksh0_726 = buffer.data(ksh0 + 726);
    const auto *ksh0_727 = buffer.data(ksh0 + 727);
    const auto *ksh0_729 = buffer.data(ksh0 + 729);
    const auto *ksh0_730 = buffer.data(ksh0 + 730);
    const auto *ksh0_731 = buffer.data(ksh0 + 731);
    const auto *ksh0_732 = buffer.data(ksh0 + 732);
    const auto *ksh0_733 = buffer.data(ksh0 + 733);
    const auto *ksh0_735 = buffer.data(ksh0 + 735);
    const auto *ksh0_737 = buffer.data(ksh0 + 737);
    const auto *ksh0_738 = buffer.data(ksh0 + 738);
    const auto *ksh0_740 = buffer.data(ksh0 + 740);
    const auto *ksh0_741 = buffer.data(ksh0 + 741);
    const auto *ksh0_742 = buffer.data(ksh0 + 742);
    const auto *ksh0_744 = buffer.data(ksh0 + 744);
    const auto *ksh0_745 = buffer.data(ksh0 + 745);
    const auto *ksh0_746 = buffer.data(ksh0 + 746);
    const auto *ksh0_747 = buffer.data(ksh0 + 747);
    const auto *ksh0_749 = buffer.data(ksh0 + 749);
    const auto *ksh0_750 = buffer.data(ksh0 + 750);
    const auto *ksh0_751 = buffer.data(ksh0 + 751);
    const auto *ksh0_752 = buffer.data(ksh0 + 752);
    const auto *ksh0_753 = buffer.data(ksh0 + 753);
    const auto *ksh0_754 = buffer.data(ksh0 + 754);
    const auto *ksh0_755 = buffer.data(ksh0 + 755);

    const auto *ksh1_693 = buffer.data(ksh1 + 693);
    const auto *ksh1_694 = buffer.data(ksh1 + 694);
    const auto *ksh1_695 = buffer.data(ksh1 + 695);
    const auto *ksh1_696 = buffer.data(ksh1 + 696);
    const auto *ksh1_697 = buffer.data(ksh1 + 697);
    const auto *ksh1_698 = buffer.data(ksh1 + 698);
    const auto *ksh1_699 = buffer.data(ksh1 + 699);
    const auto *ksh1_700 = buffer.data(ksh1 + 700);
    const auto *ksh1_701 = buffer.data(ksh1 + 701);
    const auto *ksh1_702 = buffer.data(ksh1 + 702);
    const auto *ksh1_703 = buffer.data(ksh1 + 703);
    const auto *ksh1_704 = buffer.data(ksh1 + 704);
    const auto *ksh1_705 = buffer.data(ksh1 + 705);
    const auto *ksh1_706 = buffer.data(ksh1 + 706);
    const auto *ksh1_707 = buffer.data(ksh1 + 707);
    const auto *ksh1_708 = buffer.data(ksh1 + 708);
    const auto *ksh1_709 = buffer.data(ksh1 + 709);
    const auto *ksh1_710 = buffer.data(ksh1 + 710);
    const auto *ksh1_711 = buffer.data(ksh1 + 711);
    const auto *ksh1_712 = buffer.data(ksh1 + 712);
    const auto *ksh1_713 = buffer.data(ksh1 + 713);
    const auto *ksh1_715 = buffer.data(ksh1 + 715);
    const auto *ksh1_717 = buffer.data(ksh1 + 717);
    const auto *ksh1_718 = buffer.data(ksh1 + 718);
    const auto *ksh1_720 = buffer.data(ksh1 + 720);
    const auto *ksh1_721 = buffer.data(ksh1 + 721);
    const auto *ksh1_722 = buffer.data(ksh1 + 722);
    const auto *ksh1_724 = buffer.data(ksh1 + 724);
    const auto *ksh1_725 = buffer.data(ksh1 + 725);
    const auto *ksh1_726 = buffer.data(ksh1 + 726);
    const auto *ksh1_727 = buffer.data(ksh1 + 727);
    const auto *ksh1_729 = buffer.data(ksh1 + 729);
    const auto *ksh1_730 = buffer.data(ksh1 + 730);
    const auto *ksh1_731 = buffer.data(ksh1 + 731);
    const auto *ksh1_732 = buffer.data(ksh1 + 732);
    const auto *ksh1_733 = buffer.data(ksh1 + 733);
    const auto *ksh1_735 = buffer.data(ksh1 + 735);
    const auto *ksh1_737 = buffer.data(ksh1 + 737);
    const auto *ksh1_738 = buffer.data(ksh1 + 738);
    const auto *ksh1_740 = buffer.data(ksh1 + 740);
    const auto *ksh1_741 = buffer.data(ksh1 + 741);
    const auto *ksh1_742 = buffer.data(ksh1 + 742);
    const auto *ksh1_744 = buffer.data(ksh1 + 744);
    const auto *ksh1_745 = buffer.data(ksh1 + 745);
    const auto *ksh1_746 = buffer.data(ksh1 + 746);
    const auto *ksh1_747 = buffer.data(ksh1 + 747);
    const auto *ksh1_749 = buffer.data(ksh1 + 749);
    const auto *ksh1_750 = buffer.data(ksh1 + 750);
    const auto *ksh1_751 = buffer.data(ksh1 + 751);
    const auto *ksh1_752 = buffer.data(ksh1 + 752);
    const auto *ksh1_753 = buffer.data(ksh1 + 753);
    const auto *ksh1_754 = buffer.data(ksh1 + 754);
    const auto *ksh1_755 = buffer.data(ksh1 + 755);

    const auto *ksi_924 = buffer.data(ksi + 924);
    const auto *ksi_925 = buffer.data(ksi + 925);
    const auto *ksi_926 = buffer.data(ksi + 926);
    const auto *ksi_927 = buffer.data(ksi + 927);
    const auto *ksi_928 = buffer.data(ksi + 928);
    const auto *ksi_929 = buffer.data(ksi + 929);
    const auto *ksi_930 = buffer.data(ksi + 930);
    const auto *ksi_931 = buffer.data(ksi + 931);
    const auto *ksi_932 = buffer.data(ksi + 932);
    const auto *ksi_933 = buffer.data(ksi + 933);
    const auto *ksi_934 = buffer.data(ksi + 934);
    const auto *ksi_935 = buffer.data(ksi + 935);
    const auto *ksi_936 = buffer.data(ksi + 936);
    const auto *ksi_937 = buffer.data(ksi + 937);
    const auto *ksi_938 = buffer.data(ksi + 938);
    const auto *ksi_939 = buffer.data(ksi + 939);
    const auto *ksi_940 = buffer.data(ksi + 940);
    const auto *ksi_941 = buffer.data(ksi + 941);
    const auto *ksi_942 = buffer.data(ksi + 942);
    const auto *ksi_943 = buffer.data(ksi + 943);
    const auto *ksi_944 = buffer.data(ksi + 944);
    const auto *ksi_945 = buffer.data(ksi + 945);
    const auto *ksi_946 = buffer.data(ksi + 946);
    const auto *ksi_947 = buffer.data(ksi + 947);
    const auto *ksi_948 = buffer.data(ksi + 948);
    const auto *ksi_949 = buffer.data(ksi + 949);
    const auto *ksi_950 = buffer.data(ksi + 950);
    const auto *ksi_951 = buffer.data(ksi + 951);
    const auto *ksi_953 = buffer.data(ksi + 953);
    const auto *ksi_955 = buffer.data(ksi + 955);
    const auto *ksi_956 = buffer.data(ksi + 956);
    const auto *ksi_958 = buffer.data(ksi + 958);
    const auto *ksi_959 = buffer.data(ksi + 959);
    const auto *ksi_960 = buffer.data(ksi + 960);
    const auto *ksi_962 = buffer.data(ksi + 962);
    const auto *ksi_963 = buffer.data(ksi + 963);
    const auto *ksi_964 = buffer.data(ksi + 964);
    const auto *ksi_965 = buffer.data(ksi + 965);
    const auto *ksi_967 = buffer.data(ksi + 967);
    const auto *ksi_968 = buffer.data(ksi + 968);
    const auto *ksi_969 = buffer.data(ksi + 969);
    const auto *ksi_970 = buffer.data(ksi + 970);
    const auto *ksi_971 = buffer.data(ksi + 971);
    const auto *ksi_973 = buffer.data(ksi + 973);
    const auto *ksi_974 = buffer.data(ksi + 974);
    const auto *ksi_975 = buffer.data(ksi + 975);
    const auto *ksi_976 = buffer.data(ksi + 976);
    const auto *ksi_977 = buffer.data(ksi + 977);
    const auto *ksi_978 = buffer.data(ksi + 978);
    const auto *ksi_979 = buffer.data(ksi + 979);
    const auto *ksi_980 = buffer.data(ksi + 980);
    const auto *ksi_982 = buffer.data(ksi + 982);
    const auto *ksi_983 = buffer.data(ksi + 983);
    const auto *ksi_985 = buffer.data(ksi + 985);
    const auto *ksi_986 = buffer.data(ksi + 986);
    const auto *ksi_987 = buffer.data(ksi + 987);
    const auto *ksi_989 = buffer.data(ksi + 989);
    const auto *ksi_990 = buffer.data(ksi + 990);
    const auto *ksi_991 = buffer.data(ksi + 991);
    const auto *ksi_992 = buffer.data(ksi + 992);
    const auto *ksi_994 = buffer.data(ksi + 994);
    const auto *ksi_995 = buffer.data(ksi + 995);
    const auto *ksi_996 = buffer.data(ksi + 996);
    const auto *ksi_997 = buffer.data(ksi + 997);
    const auto *ksi_998 = buffer.data(ksi + 998);
    const auto *ksi_1000 = buffer.data(ksi + 1000);
    const auto *ksi_1001 = buffer.data(ksi + 1001);
    const auto *ksi_1002 = buffer.data(ksi + 1002);
    const auto *ksi_1003 = buffer.data(ksi + 1003);
    const auto *ksi_1004 = buffer.data(ksi + 1004);
    const auto *ksi_1005 = buffer.data(ksi + 1005);
    const auto *ksi_1006 = buffer.data(ksi + 1006);
    const auto *ksi_1007 = buffer.data(ksi + 1007);

#pragma omp simd aligned(t_1188, t_1189, t_1190, pc_x, ksh0_693, ksh0_694, ksh0_695, ksh1_693, \
                         ksh1_694, ksh1_695, ksi_924, ksi_925, \
                         ksi_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_1 * ksh0_693[k]
                    - f_2 * ksh1_693[k]
                    + f_3 * pc_x[k] * ksi_924[k];

        t_1189[k] = f_19 * ksh0_694[k]
                    - f_20 * ksh1_694[k]
                    + f_3 * pc_x[k] * ksi_925[k];

        t_1190[k] = f_19 * ksh0_695[k]
                    - f_20 * ksh1_695[k]
                    + f_3 * pc_x[k] * ksi_926[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, pc_x, ksh0_696, ksh0_697, ksh0_698, ksh1_696, \
                         ksh1_697, ksh1_698, ksi_927, ksi_928, \
                         ksi_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_10 * ksh0_696[k]
                    - f_11 * ksh1_696[k]
                    + f_3 * pc_x[k] * ksi_927[k];

        t_1192[k] = f_10 * ksh0_697[k]
                    - f_11 * ksh1_697[k]
                    + f_3 * pc_x[k] * ksi_928[k];

        t_1193[k] = f_10 * ksh0_698[k]
                    - f_11 * ksh1_698[k]
                    + f_3 * pc_x[k] * ksi_929[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, pc_x, ksh0_699, ksh0_700, ksh0_701, ksh1_699, \
                         ksh1_700, ksh1_701, ksi_930, ksi_931, \
                         ksi_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = f_8 * ksh0_699[k]
                    - f_9 * ksh1_699[k]
                    + f_3 * pc_x[k] * ksi_930[k];

        t_1195[k] = f_8 * ksh0_700[k]
                    - f_9 * ksh1_700[k]
                    + f_3 * pc_x[k] * ksi_931[k];

        t_1196[k] = f_8 * ksh0_701[k]
                    - f_9 * ksh1_701[k]
                    + f_3 * pc_x[k] * ksi_932[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pc_x, ksh0_702, ksh0_703, ksh0_704, ksh1_702, \
                         ksh1_703, ksh1_704, ksi_933, ksi_934, \
                         ksi_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = f_8 * ksh0_702[k]
                    - f_9 * ksh1_702[k]
                    + f_3 * pc_x[k] * ksi_933[k];

        t_1198[k] = f_6 * ksh0_703[k]
                    - f_7 * ksh1_703[k]
                    + f_3 * pc_x[k] * ksi_934[k];

        t_1199[k] = f_6 * ksh0_704[k]
                    - f_7 * ksh1_704[k]
                    + f_3 * pc_x[k] * ksi_935[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pc_x, ksh0_705, ksh0_706, ksh0_707, ksh1_705, \
                         ksh1_706, ksh1_707, ksi_936, ksi_937, \
                         ksi_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = f_6 * ksh0_705[k]
                    - f_7 * ksh1_705[k]
                    + f_3 * pc_x[k] * ksi_936[k];

        t_1201[k] = f_6 * ksh0_706[k]
                    - f_7 * ksh1_706[k]
                    + f_3 * pc_x[k] * ksi_937[k];

        t_1202[k] = f_6 * ksh0_707[k]
                    - f_7 * ksh1_707[k]
                    + f_3 * pc_x[k] * ksi_938[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, pc_x, ksh0_708, ksh0_709, ksh0_710, ksh1_708, \
                         ksh1_709, ksh1_710, ksi_939, ksi_940, \
                         ksi_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_4 * ksh0_708[k]
                    - f_5 * ksh1_708[k]
                    + f_3 * pc_x[k] * ksi_939[k];

        t_1204[k] = f_4 * ksh0_709[k]
                    - f_5 * ksh1_709[k]
                    + f_3 * pc_x[k] * ksi_940[k];

        t_1205[k] = f_4 * ksh0_710[k]
                    - f_5 * ksh1_710[k]
                    + f_3 * pc_x[k] * ksi_941[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pc_x, ksh0_711, ksh0_712, ksh0_713, \
                         ksh1_711, ksh1_712, ksh1_713, ksi_942, ksi_943, ksi_944, \
                         ksi_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_4 * ksh0_711[k]
                    - f_5 * ksh1_711[k]
                    + f_3 * pc_x[k] * ksi_942[k];

        t_1207[k] = f_4 * ksh0_712[k]
                    - f_5 * ksh1_712[k]
                    + f_3 * pc_x[k] * ksi_943[k];

        t_1208[k] = f_4 * ksh0_713[k]
                    - f_5 * ksh1_713[k]
                    + f_3 * pc_x[k] * ksi_944[k];

        t_1209[k] = f_3 * pc_x[k] * ksi_945[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, t_1214, t_1215, pc_x, ksi_946, \
                         ksi_947, ksi_948, ksi_949, ksi_950, ksi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_3 * pc_x[k] * ksi_946[k];

        t_1211[k] = f_3 * pc_x[k] * ksi_947[k];

        t_1212[k] = f_3 * pc_x[k] * ksi_948[k];

        t_1213[k] = f_3 * pc_x[k] * ksi_949[k];

        t_1214[k] = f_3 * pc_x[k] * ksi_950[k];

        t_1215[k] = f_3 * pc_x[k] * ksi_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_y, pc_z, isi_721, isi_749, isi_751, \
                         ksh0_708, ksh0_710, ksh1_708, ksh1_710, ksi_945, \
                         ksi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_14 * isi_749[k]
                    + f_1 * ksh0_708[k]
                    - f_2 * ksh1_708[k]
                    + f_3 * pc_y[k] * ksi_945[k];

        t_1217[k] = f_17 * isi_721[k]
                    + f_3 * pc_z[k] * ksi_945[k];

        t_1218[k] = f_14 * isi_751[k]
                    + f_10 * ksh0_710[k]
                    - f_11 * ksh1_710[k]
                    + f_3 * pc_y[k] * ksi_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pc_y, isi_752, isi_753, isi_754, ksh0_711, \
                         ksh0_712, ksh0_713, ksh1_711, ksh1_712, ksh1_713, ksi_948, ksi_949, \
                         ksi_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_14 * isi_752[k]
                    + f_8 * ksh0_711[k]
                    - f_9 * ksh1_711[k]
                    + f_3 * pc_y[k] * ksi_948[k];

        t_1220[k] = f_14 * isi_753[k]
                    + f_6 * ksh0_712[k]
                    - f_7 * ksh1_712[k]
                    + f_3 * pc_y[k] * ksi_949[k];

        t_1221[k] = f_14 * isi_754[k]
                    + f_4 * ksh0_713[k]
                    - f_5 * ksh1_713[k]
                    + f_3 * pc_y[k] * ksi_950[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, pa_y, pc_y, pc_z, isk0_972, isi_727, isi_755, \
                         isk1_972, ksh0_713, ksh1_713, ksi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * isi_755[k]
                    + f_3 * pc_y[k] * ksi_951[k];

        t_1223[k] = f_17 * isi_727[k]
                    + f_1 * ksh0_713[k]
                    - f_2 * ksh1_713[k]
                    + f_3 * pc_z[k] * ksi_951[k];

        t_1224[k] = pa_y[k] * isk0_972[k]
                    - f_12 * pc_y[k] * isk1_972[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pa_y, pc_x, pc_y, isk0_974, isk1_974, \
                         ksh0_715, ksh0_717, ksh1_715, ksh1_717, ksi_953, \
                         ksi_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_19 * ksh0_715[k]
                    - f_20 * ksh1_715[k]
                    + f_3 * pc_x[k] * ksi_953[k];

        t_1226[k] = pa_y[k] * isk0_974[k]
                    - f_12 * pc_y[k] * isk1_974[k];

        t_1227[k] = f_10 * ksh0_717[k]
                    - f_11 * ksh1_717[k]
                    + f_3 * pc_x[k] * ksi_955[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, pa_y, pc_x, pc_y, isk0_977, isk1_977, \
                         ksh0_718, ksh0_720, ksh1_718, ksh1_720, ksi_956, \
                         ksi_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = f_10 * ksh0_718[k]
                    - f_11 * ksh1_718[k]
                    + f_3 * pc_x[k] * ksi_956[k];

        t_1229[k] = pa_y[k] * isk0_977[k]
                    - f_12 * pc_y[k] * isk1_977[k];

        t_1230[k] = f_8 * ksh0_720[k]
                    - f_9 * ksh1_720[k]
                    + f_3 * pc_x[k] * ksi_958[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, pa_y, pc_x, pc_y, isk0_981, isk1_981, \
                         ksh0_721, ksh0_722, ksh1_721, ksh1_722, ksi_959, \
                         ksi_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_8 * ksh0_721[k]
                    - f_9 * ksh1_721[k]
                    + f_3 * pc_x[k] * ksi_959[k];

        t_1232[k] = f_8 * ksh0_722[k]
                    - f_9 * ksh1_722[k]
                    + f_3 * pc_x[k] * ksi_960[k];

        t_1233[k] = pa_y[k] * isk0_981[k]
                    - f_12 * pc_y[k] * isk1_981[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pc_x, ksh0_724, ksh0_725, ksh0_726, ksh1_724, \
                         ksh1_725, ksh1_726, ksi_962, ksi_963, \
                         ksi_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = f_6 * ksh0_724[k]
                    - f_7 * ksh1_724[k]
                    + f_3 * pc_x[k] * ksi_962[k];

        t_1235[k] = f_6 * ksh0_725[k]
                    - f_7 * ksh1_725[k]
                    + f_3 * pc_x[k] * ksi_963[k];

        t_1236[k] = f_6 * ksh0_726[k]
                    - f_7 * ksh1_726[k]
                    + f_3 * pc_x[k] * ksi_964[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, pa_y, pc_x, pc_y, isk0_986, isk1_986, \
                         ksh0_727, ksh0_729, ksh1_727, ksh1_729, ksi_965, \
                         ksi_967 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_6 * ksh0_727[k]
                    - f_7 * ksh1_727[k]
                    + f_3 * pc_x[k] * ksi_965[k];

        t_1238[k] = pa_y[k] * isk0_986[k]
                    - f_12 * pc_y[k] * isk1_986[k];

        t_1239[k] = f_4 * ksh0_729[k]
                    - f_5 * ksh1_729[k]
                    + f_3 * pc_x[k] * ksi_967[k];
    }

#pragma omp simd aligned(t_1240, t_1241, t_1242, pc_x, ksh0_730, ksh0_731, ksh0_732, ksh1_730, \
                         ksh1_731, ksh1_732, ksi_968, ksi_969, \
                         ksi_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1240[k] = f_4 * ksh0_730[k]
                    - f_5 * ksh1_730[k]
                    + f_3 * pc_x[k] * ksi_968[k];

        t_1241[k] = f_4 * ksh0_731[k]
                    - f_5 * ksh1_731[k]
                    + f_3 * pc_x[k] * ksi_969[k];

        t_1242[k] = f_4 * ksh0_732[k]
                    - f_5 * ksh1_732[k]
                    + f_3 * pc_x[k] * ksi_970[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, t_1247, pa_y, pc_x, pc_y, isk0_992, \
                         isk1_992, ksh0_733, ksh1_733, ksi_971, ksi_973, ksi_974, \
                         ksi_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_4 * ksh0_733[k]
                    - f_5 * ksh1_733[k]
                    + f_3 * pc_x[k] * ksi_971[k];

        t_1244[k] = pa_y[k] * isk0_992[k]
                    - f_12 * pc_y[k] * isk1_992[k];

        t_1245[k] = f_3 * pc_x[k] * ksi_973[k];

        t_1246[k] = f_3 * pc_x[k] * ksi_974[k];

        t_1247[k] = f_3 * pc_x[k] * ksi_975[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, pa_y, pc_x, pc_y, isk0_1000, \
                         isi_777, isk1_1000, ksi_976, ksi_977, ksi_978, \
                         ksi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_3 * pc_x[k] * ksi_976[k];

        t_1249[k] = f_3 * pc_x[k] * ksi_977[k];

        t_1250[k] = f_3 * pc_x[k] * ksi_978[k];

        t_1251[k] = f_3 * pc_x[k] * ksi_979[k];

        t_1252[k] = pa_y[k] * isk0_1000[k]
                    + f_0 * isi_777[k]
                    - f_12 * pc_y[k] * isk1_1000[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, pa_y, pc_y, pc_z, isk0_1002, isk0_1003, \
                         isi_749, isi_779, isi_780, isk1_1002, isk1_1003, \
                         ksi_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = f_18 * isi_749[k]
                    + f_3 * pc_z[k] * ksi_973[k];

        t_1254[k] = pa_y[k] * isk0_1002[k]
                    + f_17 * isi_779[k]
                    - f_12 * pc_y[k] * isk1_1002[k];

        t_1255[k] = pa_y[k] * isk0_1003[k]
                    + f_16 * isi_780[k]
                    - f_12 * pc_y[k] * isk1_1003[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pa_y, pc_y, isk0_1004, isk0_1005, \
                         isk0_1007, isi_781, isi_782, isi_783, isk1_1004, isk1_1005, \
                         isk1_1007, ksi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = pa_y[k] * isk0_1004[k]
                    + f_15 * isi_781[k]
                    - f_12 * pc_y[k] * isk1_1004[k];

        t_1257[k] = pa_y[k] * isk0_1005[k]
                    + f_14 * isi_782[k]
                    - f_12 * pc_y[k] * isk1_1005[k];

        t_1258[k] = f_13 * isi_783[k]
                    + f_3 * pc_y[k] * ksi_979[k];

        t_1259[k] = pa_y[k] * isk0_1007[k]
                    - f_12 * pc_y[k] * isk1_1007[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, t_1264, pc_x, pc_y, ksh0_735, \
                         ksh0_737, ksh0_738, ksh1_735, ksh1_737, ksh1_738, ksi_980, ksi_982, \
                         ksi_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_1 * ksh0_735[k]
                    - f_2 * ksh1_735[k]
                    + f_3 * pc_x[k] * ksi_980[k];

        t_1261[k] = f_3 * pc_y[k] * ksi_980[k];

        t_1262[k] = f_19 * ksh0_737[k]
                    - f_20 * ksh1_737[k]
                    + f_3 * pc_x[k] * ksi_982[k];

        t_1263[k] = f_10 * ksh0_738[k]
                    - f_11 * ksh1_738[k]
                    + f_3 * pc_x[k] * ksi_983[k];

        t_1264[k] = f_3 * pc_y[k] * ksi_982[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, t_1268, pc_x, pc_y, ksh0_740, ksh0_741, \
                         ksh0_742, ksh1_740, ksh1_741, ksh1_742, ksi_985, ksi_986, \
                         ksi_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_10 * ksh0_740[k]
                    - f_11 * ksh1_740[k]
                    + f_3 * pc_x[k] * ksi_985[k];

        t_1266[k] = f_8 * ksh0_741[k]
                    - f_9 * ksh1_741[k]
                    + f_3 * pc_x[k] * ksi_986[k];

        t_1267[k] = f_8 * ksh0_742[k]
                    - f_9 * ksh1_742[k]
                    + f_3 * pc_x[k] * ksi_987[k];

        t_1268[k] = f_3 * pc_y[k] * ksi_985[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, pc_x, ksh0_744, ksh0_745, ksh0_746, ksh1_744, \
                         ksh1_745, ksh1_746, ksi_989, ksi_990, \
                         ksi_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = f_8 * ksh0_744[k]
                    - f_9 * ksh1_744[k]
                    + f_3 * pc_x[k] * ksi_989[k];

        t_1270[k] = f_6 * ksh0_745[k]
                    - f_7 * ksh1_745[k]
                    + f_3 * pc_x[k] * ksi_990[k];

        t_1271[k] = f_6 * ksh0_746[k]
                    - f_7 * ksh1_746[k]
                    + f_3 * pc_x[k] * ksi_991[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, ksh0_747, ksh0_749, \
                         ksh0_750, ksh1_747, ksh1_749, ksh1_750, ksi_989, ksi_992, ksi_994, \
                         ksi_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_6 * ksh0_747[k]
                    - f_7 * ksh1_747[k]
                    + f_3 * pc_x[k] * ksi_992[k];

        t_1273[k] = f_3 * pc_y[k] * ksi_989[k];

        t_1274[k] = f_6 * ksh0_749[k]
                    - f_7 * ksh1_749[k]
                    + f_3 * pc_x[k] * ksi_994[k];

        t_1275[k] = f_4 * ksh0_750[k]
                    - f_5 * ksh1_750[k]
                    + f_3 * pc_x[k] * ksi_995[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, t_1279, pc_x, pc_y, ksh0_751, ksh0_752, \
                         ksh0_753, ksh1_751, ksh1_752, ksh1_753, ksi_994, ksi_996, ksi_997, \
                         ksi_998 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_4 * ksh0_751[k]
                    - f_5 * ksh1_751[k]
                    + f_3 * pc_x[k] * ksi_996[k];

        t_1277[k] = f_4 * ksh0_752[k]
                    - f_5 * ksh1_752[k]
                    + f_3 * pc_x[k] * ksi_997[k];

        t_1278[k] = f_4 * ksh0_753[k]
                    - f_5 * ksh1_753[k]
                    + f_3 * pc_x[k] * ksi_998[k];

        t_1279[k] = f_3 * pc_y[k] * ksi_994[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, t_1283, t_1284, t_1285, pc_x, ksh0_755, \
                         ksh1_755, ksi_1000, ksi_1001, ksi_1002, ksi_1003, ksi_1004, \
                         ksi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = f_4 * ksh0_755[k]
                    - f_5 * ksh1_755[k]
                    + f_3 * pc_x[k] * ksi_1000[k];

        t_1281[k] = f_3 * pc_x[k] * ksi_1001[k];

        t_1282[k] = f_3 * pc_x[k] * ksi_1002[k];

        t_1283[k] = f_3 * pc_x[k] * ksi_1003[k];

        t_1284[k] = f_3 * pc_x[k] * ksi_1004[k];

        t_1285[k] = f_3 * pc_x[k] * ksi_1005[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, t_1289, pc_x, pc_y, ksh0_750, ksh0_751, \
                         ksh1_750, ksh1_751, ksi_1001, ksi_1002, ksi_1006, \
                         ksi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_3 * pc_x[k] * ksi_1006[k];

        t_1287[k] = f_3 * pc_x[k] * ksi_1007[k];

        t_1288[k] = f_1 * ksh0_750[k]
                    - f_2 * ksh1_750[k]
                    + f_3 * pc_y[k] * ksi_1001[k];

        t_1289[k] = f_19 * ksh0_751[k]
                    - f_20 * ksh1_751[k]
                    + f_3 * pc_y[k] * ksi_1002[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, pc_y, ksh0_752, ksh0_753, ksh0_754, ksh1_752, \
                         ksh1_753, ksh1_754, ksi_1003, ksi_1004, \
                         ksi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_10 * ksh0_752[k]
                    - f_11 * ksh1_752[k]
                    + f_3 * pc_y[k] * ksi_1003[k];

        t_1291[k] = f_8 * ksh0_753[k]
                    - f_9 * ksh1_753[k]
                    + f_3 * pc_y[k] * ksi_1004[k];

        t_1292[k] = f_6 * ksh0_754[k]
                    - f_7 * ksh1_754[k]
                    + f_3 * pc_y[k] * ksi_1005[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, pc_y, pc_z, isi_783, ksh0_755, ksh1_755, \
                         ksi_1006, ksi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_4 * ksh0_755[k]
                    - f_5 * ksh1_755[k]
                    + f_3 * pc_y[k] * ksi_1006[k];

        t_1294[k] = f_3 * pc_y[k] * ksi_1007[k];

        t_1295[k] = f_0 * isi_783[k]
                    + f_1 * ksh0_755[k]
                    - f_2 * ksh1_755[k]
                    + f_3 * pc_z[k] * ksi_1007[k];
    }
}

auto
compute_prim_ksk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t isk0, const size_t isi,
                                                   const size_t isk1, const size_t ksh0,
                                                   const size_t ksh1, const size_t ksi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksi, ncols, gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, isk0, isi,
                                                              isk1, ksh0, ksh1, ksi, ncols,
                                                              gamma, p, q);

    compute_prim_ksk_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, isk0,
                                                               isi, isk1, ksh0, ksh1, ksi,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
