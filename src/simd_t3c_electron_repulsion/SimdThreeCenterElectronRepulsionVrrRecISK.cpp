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


#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);

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

    const auto *hsk0_0 = buffer.data(hsk0 + 0);
    const auto *hsk0_3 = buffer.data(hsk0 + 3);
    const auto *hsk0_5 = buffer.data(hsk0 + 5);
    const auto *hsk0_6 = buffer.data(hsk0 + 6);
    const auto *hsk0_9 = buffer.data(hsk0 + 9);
    const auto *hsk0_10 = buffer.data(hsk0 + 10);
    const auto *hsk0_14 = buffer.data(hsk0 + 14);
    const auto *hsk0_15 = buffer.data(hsk0 + 15);
    const auto *hsk0_20 = buffer.data(hsk0 + 20);
    const auto *hsk0_28 = buffer.data(hsk0 + 28);
    const auto *hsk0_35 = buffer.data(hsk0 + 35);

    const auto *hsi_0 = buffer.data(hsi + 0);
    const auto *hsi_1 = buffer.data(hsi + 1);
    const auto *hsi_2 = buffer.data(hsi + 2);
    const auto *hsi_3 = buffer.data(hsi + 3);
    const auto *hsi_5 = buffer.data(hsi + 5);
    const auto *hsi_6 = buffer.data(hsi + 6);
    const auto *hsi_9 = buffer.data(hsi + 9);
    const auto *hsi_10 = buffer.data(hsi + 10);
    const auto *hsi_14 = buffer.data(hsi + 14);
    const auto *hsi_21 = buffer.data(hsi + 21);
    const auto *hsi_23 = buffer.data(hsi + 23);
    const auto *hsi_24 = buffer.data(hsi + 24);
    const auto *hsi_25 = buffer.data(hsi + 25);
    const auto *hsi_27 = buffer.data(hsi + 27);
    const auto *hsi_28 = buffer.data(hsi + 28);
    const auto *hsi_33 = buffer.data(hsi + 33);
    const auto *hsi_37 = buffer.data(hsi + 37);
    const auto *hsi_42 = buffer.data(hsi + 42);
    const auto *hsi_49 = buffer.data(hsi + 49);
    const auto *hsi_51 = buffer.data(hsi + 51);
    const auto *hsi_52 = buffer.data(hsi + 52);
    const auto *hsi_53 = buffer.data(hsi + 53);
    const auto *hsi_54 = buffer.data(hsi + 54);
    const auto *hsi_55 = buffer.data(hsi + 55);
    const auto *hsi_77 = buffer.data(hsi + 77);
    const auto *hsi_78 = buffer.data(hsi + 78);
    const auto *hsi_79 = buffer.data(hsi + 79);
    const auto *hsi_80 = buffer.data(hsi + 80);
    const auto *hsi_81 = buffer.data(hsi + 81);
    const auto *hsi_83 = buffer.data(hsi + 83);
    const auto *hsi_84 = buffer.data(hsi + 84);
    const auto *hsi_87 = buffer.data(hsi + 87);
    const auto *hsi_90 = buffer.data(hsi + 90);
    const auto *hsi_94 = buffer.data(hsi + 94);
    const auto *hsi_99 = buffer.data(hsi + 99);

    const auto *hsk1_0 = buffer.data(hsk1 + 0);
    const auto *hsk1_3 = buffer.data(hsk1 + 3);
    const auto *hsk1_5 = buffer.data(hsk1 + 5);
    const auto *hsk1_6 = buffer.data(hsk1 + 6);
    const auto *hsk1_9 = buffer.data(hsk1 + 9);
    const auto *hsk1_10 = buffer.data(hsk1 + 10);
    const auto *hsk1_14 = buffer.data(hsk1 + 14);
    const auto *hsk1_15 = buffer.data(hsk1 + 15);
    const auto *hsk1_20 = buffer.data(hsk1 + 20);
    const auto *hsk1_28 = buffer.data(hsk1 + 28);
    const auto *hsk1_35 = buffer.data(hsk1 + 35);

    const auto *ish0_0 = buffer.data(ish0 + 0);
    const auto *ish0_1 = buffer.data(ish0 + 1);
    const auto *ish0_2 = buffer.data(ish0 + 2);
    const auto *ish0_3 = buffer.data(ish0 + 3);
    const auto *ish0_5 = buffer.data(ish0 + 5);
    const auto *ish0_6 = buffer.data(ish0 + 6);
    const auto *ish0_8 = buffer.data(ish0 + 8);
    const auto *ish0_9 = buffer.data(ish0 + 9);
    const auto *ish0_15 = buffer.data(ish0 + 15);
    const auto *ish0_17 = buffer.data(ish0 + 17);
    const auto *ish0_18 = buffer.data(ish0 + 18);
    const auto *ish0_19 = buffer.data(ish0 + 19);
    const auto *ish0_20 = buffer.data(ish0 + 20);
    const auto *ish0_24 = buffer.data(ish0 + 24);
    const auto *ish0_27 = buffer.data(ish0 + 27);
    const auto *ish0_28 = buffer.data(ish0 + 28);
    const auto *ish0_36 = buffer.data(ish0 + 36);
    const auto *ish0_37 = buffer.data(ish0 + 37);
    const auto *ish0_38 = buffer.data(ish0 + 38);
    const auto *ish0_39 = buffer.data(ish0 + 39);
    const auto *ish0_44 = buffer.data(ish0 + 44);
    const auto *ish0_46 = buffer.data(ish0 + 46);
    const auto *ish0_47 = buffer.data(ish0 + 47);
    const auto *ish0_49 = buffer.data(ish0 + 49);
    const auto *ish0_50 = buffer.data(ish0 + 50);
    const auto *ish0_51 = buffer.data(ish0 + 51);
    const auto *ish0_58 = buffer.data(ish0 + 58);
    const auto *ish0_59 = buffer.data(ish0 + 59);
    const auto *ish0_60 = buffer.data(ish0 + 60);
    const auto *ish0_61 = buffer.data(ish0 + 61);
    const auto *ish0_62 = buffer.data(ish0 + 62);
    const auto *ish0_63 = buffer.data(ish0 + 63);
    const auto *ish0_65 = buffer.data(ish0 + 65);
    const auto *ish0_66 = buffer.data(ish0 + 66);
    const auto *ish0_68 = buffer.data(ish0 + 68);
    const auto *ish0_69 = buffer.data(ish0 + 69);
    const auto *ish0_70 = buffer.data(ish0 + 70);
    const auto *ish0_72 = buffer.data(ish0 + 72);
    const auto *ish0_73 = buffer.data(ish0 + 73);
    const auto *ish0_78 = buffer.data(ish0 + 78);

    const auto *ish1_0 = buffer.data(ish1 + 0);
    const auto *ish1_1 = buffer.data(ish1 + 1);
    const auto *ish1_2 = buffer.data(ish1 + 2);
    const auto *ish1_3 = buffer.data(ish1 + 3);
    const auto *ish1_5 = buffer.data(ish1 + 5);
    const auto *ish1_6 = buffer.data(ish1 + 6);
    const auto *ish1_8 = buffer.data(ish1 + 8);
    const auto *ish1_9 = buffer.data(ish1 + 9);
    const auto *ish1_15 = buffer.data(ish1 + 15);
    const auto *ish1_17 = buffer.data(ish1 + 17);
    const auto *ish1_18 = buffer.data(ish1 + 18);
    const auto *ish1_19 = buffer.data(ish1 + 19);
    const auto *ish1_20 = buffer.data(ish1 + 20);
    const auto *ish1_24 = buffer.data(ish1 + 24);
    const auto *ish1_27 = buffer.data(ish1 + 27);
    const auto *ish1_28 = buffer.data(ish1 + 28);
    const auto *ish1_36 = buffer.data(ish1 + 36);
    const auto *ish1_37 = buffer.data(ish1 + 37);
    const auto *ish1_38 = buffer.data(ish1 + 38);
    const auto *ish1_39 = buffer.data(ish1 + 39);
    const auto *ish1_44 = buffer.data(ish1 + 44);
    const auto *ish1_46 = buffer.data(ish1 + 46);
    const auto *ish1_47 = buffer.data(ish1 + 47);
    const auto *ish1_49 = buffer.data(ish1 + 49);
    const auto *ish1_50 = buffer.data(ish1 + 50);
    const auto *ish1_51 = buffer.data(ish1 + 51);
    const auto *ish1_58 = buffer.data(ish1 + 58);
    const auto *ish1_59 = buffer.data(ish1 + 59);
    const auto *ish1_60 = buffer.data(ish1 + 60);
    const auto *ish1_61 = buffer.data(ish1 + 61);
    const auto *ish1_62 = buffer.data(ish1 + 62);
    const auto *ish1_63 = buffer.data(ish1 + 63);
    const auto *ish1_65 = buffer.data(ish1 + 65);
    const auto *ish1_66 = buffer.data(ish1 + 66);
    const auto *ish1_68 = buffer.data(ish1 + 68);
    const auto *ish1_69 = buffer.data(ish1 + 69);
    const auto *ish1_70 = buffer.data(ish1 + 70);
    const auto *ish1_72 = buffer.data(ish1 + 72);
    const auto *ish1_73 = buffer.data(ish1 + 73);
    const auto *ish1_78 = buffer.data(ish1 + 78);

    const auto *isi_0 = buffer.data(isi + 0);
    const auto *isi_1 = buffer.data(isi + 1);
    const auto *isi_2 = buffer.data(isi + 2);
    const auto *isi_3 = buffer.data(isi + 3);
    const auto *isi_5 = buffer.data(isi + 5);
    const auto *isi_6 = buffer.data(isi + 6);
    const auto *isi_8 = buffer.data(isi + 8);
    const auto *isi_9 = buffer.data(isi + 9);
    const auto *isi_10 = buffer.data(isi + 10);
    const auto *isi_12 = buffer.data(isi + 12);
    const auto *isi_13 = buffer.data(isi + 13);
    const auto *isi_14 = buffer.data(isi + 14);
    const auto *isi_15 = buffer.data(isi + 15);
    const auto *isi_20 = buffer.data(isi + 20);
    const auto *isi_21 = buffer.data(isi + 21);
    const auto *isi_23 = buffer.data(isi + 23);
    const auto *isi_24 = buffer.data(isi + 24);
    const auto *isi_25 = buffer.data(isi + 25);
    const auto *isi_26 = buffer.data(isi + 26);
    const auto *isi_27 = buffer.data(isi + 27);
    const auto *isi_28 = buffer.data(isi + 28);
    const auto *isi_29 = buffer.data(isi + 29);
    const auto *isi_31 = buffer.data(isi + 31);
    const auto *isi_33 = buffer.data(isi + 33);
    const auto *isi_34 = buffer.data(isi + 34);
    const auto *isi_35 = buffer.data(isi + 35);
    const auto *isi_37 = buffer.data(isi + 37);
    const auto *isi_38 = buffer.data(isi + 38);
    const auto *isi_39 = buffer.data(isi + 39);
    const auto *isi_40 = buffer.data(isi + 40);
    const auto *isi_42 = buffer.data(isi + 42);
    const auto *isi_43 = buffer.data(isi + 43);
    const auto *isi_49 = buffer.data(isi + 49);
    const auto *isi_50 = buffer.data(isi + 50);
    const auto *isi_51 = buffer.data(isi + 51);
    const auto *isi_52 = buffer.data(isi + 52);
    const auto *isi_53 = buffer.data(isi + 53);
    const auto *isi_54 = buffer.data(isi + 54);
    const auto *isi_55 = buffer.data(isi + 55);
    const auto *isi_56 = buffer.data(isi + 56);
    const auto *isi_58 = buffer.data(isi + 58);
    const auto *isi_60 = buffer.data(isi + 60);
    const auto *isi_61 = buffer.data(isi + 61);
    const auto *isi_63 = buffer.data(isi + 63);
    const auto *isi_64 = buffer.data(isi + 64);
    const auto *isi_65 = buffer.data(isi + 65);
    const auto *isi_67 = buffer.data(isi + 67);
    const auto *isi_68 = buffer.data(isi + 68);
    const auto *isi_69 = buffer.data(isi + 69);
    const auto *isi_70 = buffer.data(isi + 70);
    const auto *isi_76 = buffer.data(isi + 76);
    const auto *isi_77 = buffer.data(isi + 77);
    const auto *isi_78 = buffer.data(isi + 78);
    const auto *isi_79 = buffer.data(isi + 79);
    const auto *isi_80 = buffer.data(isi + 80);
    const auto *isi_81 = buffer.data(isi + 81);
    const auto *isi_82 = buffer.data(isi + 82);
    const auto *isi_83 = buffer.data(isi + 83);
    const auto *isi_84 = buffer.data(isi + 84);
    const auto *isi_85 = buffer.data(isi + 85);
    const auto *isi_86 = buffer.data(isi + 86);
    const auto *isi_87 = buffer.data(isi + 87);
    const auto *isi_89 = buffer.data(isi + 89);
    const auto *isi_90 = buffer.data(isi + 90);
    const auto *isi_91 = buffer.data(isi + 91);
    const auto *isi_93 = buffer.data(isi + 93);
    const auto *isi_94 = buffer.data(isi + 94);
    const auto *isi_95 = buffer.data(isi + 95);
    const auto *isi_96 = buffer.data(isi + 96);
    const auto *isi_98 = buffer.data(isi + 98);
    const auto *isi_99 = buffer.data(isi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, hsi_0, ish0_0, \
                         ish1_0, isi_0, isi_1, isi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsi_0[k]
                 + f_1 * ish0_0[k]
                 - f_2 * ish1_0[k]
                 + f_3 * pc_x[k] * isi_0[k];

        t_1[k] = f_3 * pc_y[k] * isi_0[k];

        t_2[k] = f_3 * pc_z[k] * isi_0[k];

        t_3[k] = f_4 * ish0_0[k]
                 - f_5 * ish1_0[k]
                 + f_3 * pc_y[k] * isi_1[k];

        t_4[k] = f_3 * pc_y[k] * isi_2[k];

        t_5[k] = f_4 * ish0_0[k]
                 - f_5 * ish1_0[k]
                 + f_3 * pc_z[k] * isi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, ish0_1, ish0_2, ish0_3, ish1_1, \
                         ish1_2, ish1_3, isi_3, isi_5, isi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ish0_1[k]
                 - f_7 * ish1_1[k]
                 + f_3 * pc_y[k] * isi_3[k];

        t_7[k] = f_3 * pc_z[k] * isi_3[k];

        t_8[k] = f_3 * pc_y[k] * isi_5[k];

        t_9[k] = f_6 * ish0_2[k]
                 - f_7 * ish1_2[k]
                 + f_3 * pc_z[k] * isi_5[k];

        t_10[k] = f_8 * ish0_3[k]
                  - f_9 * ish1_3[k]
                  + f_3 * pc_y[k] * isi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, ish0_5, ish0_6, \
                         ish1_5, ish1_6, isi_6, isi_8, isi_9, isi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * isi_6[k];

        t_12[k] = f_4 * ish0_5[k]
                  - f_5 * ish1_5[k]
                  + f_3 * pc_y[k] * isi_8[k];

        t_13[k] = f_3 * pc_y[k] * isi_9[k];

        t_14[k] = f_8 * ish0_5[k]
                  - f_9 * ish1_5[k]
                  + f_3 * pc_z[k] * isi_9[k];

        t_15[k] = f_10 * ish0_6[k]
                  - f_11 * ish1_6[k]
                  + f_3 * pc_y[k] * isi_10[k];

        t_16[k] = f_3 * pc_z[k] * isi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, ish0_8, ish0_9, ish1_8, ish1_9, \
                         isi_12, isi_13, isi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * ish0_8[k]
                  - f_7 * ish1_8[k]
                  + f_3 * pc_y[k] * isi_12[k];

        t_18[k] = f_4 * ish0_9[k]
                  - f_5 * ish1_9[k]
                  + f_3 * pc_y[k] * isi_13[k];

        t_19[k] = f_3 * pc_y[k] * isi_14[k];

        t_20[k] = f_10 * ish0_9[k]
                  - f_11 * ish1_9[k]
                  + f_3 * pc_z[k] * isi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, hsi_21, hsi_23, hsi_24, \
                         hsi_25, isi_15, isi_21, isi_23, isi_24, \
                         isi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * hsi_21[k]
                  + f_3 * pc_x[k] * isi_21[k];

        t_22[k] = f_3 * pc_z[k] * isi_15[k];

        t_23[k] = f_0 * hsi_23[k]
                  + f_3 * pc_x[k] * isi_23[k];

        t_24[k] = f_0 * hsi_24[k]
                  + f_3 * pc_x[k] * isi_24[k];

        t_25[k] = f_0 * hsi_25[k]
                  + f_3 * pc_x[k] * isi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, hsi_27, ish0_15, ish1_15, \
                         isi_20, isi_21, isi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * isi_20[k];

        t_27[k] = f_0 * hsi_27[k]
                  + f_3 * pc_x[k] * isi_27[k];

        t_28[k] = f_1 * ish0_15[k]
                  - f_2 * ish1_15[k]
                  + f_3 * pc_y[k] * isi_21[k];

        t_29[k] = f_3 * pc_z[k] * isi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, ish0_17, ish0_18, ish0_19, ish1_17, ish1_18, \
                         ish1_19, isi_23, isi_24, isi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * ish0_17[k]
                  - f_11 * ish1_17[k]
                  + f_3 * pc_y[k] * isi_23[k];

        t_31[k] = f_8 * ish0_18[k]
                  - f_9 * ish1_18[k]
                  + f_3 * pc_y[k] * isi_24[k];

        t_32[k] = f_6 * ish0_19[k]
                  - f_7 * ish1_19[k]
                  + f_3 * pc_y[k] * isi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, hsk0_0, hsi_0, \
                         hsk1_0, ish0_20, ish1_20, isi_26, isi_27, \
                         isi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * ish0_20[k]
                  - f_5 * ish1_20[k]
                  + f_3 * pc_y[k] * isi_26[k];

        t_34[k] = f_3 * pc_y[k] * isi_27[k];

        t_35[k] = f_1 * ish0_20[k]
                  - f_2 * ish1_20[k]
                  + f_3 * pc_z[k] * isi_27[k];

        t_36[k] = pa_y[k] * hsk0_0[k]
                  - f_12 * pc_y[k] * hsk1_0[k];

        t_37[k] = f_13 * hsi_0[k]
                  + f_3 * pc_y[k] * isi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, hsk0_3, hsk0_5, hsi_1, \
                         hsk1_3, hsk1_5, isi_28, isi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * isi_28[k];

        t_39[k] = pa_y[k] * hsk0_3[k]
                  + f_14 * hsi_1[k]
                  - f_12 * pc_y[k] * hsk1_3[k];

        t_40[k] = f_3 * pc_z[k] * isi_29[k];

        t_41[k] = pa_y[k] * hsk0_5[k]
                  - f_12 * pc_y[k] * hsk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, hsk0_6, hsk0_9, hsi_3, \
                         hsi_5, hsk1_6, hsk1_9, isi_31, isi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * hsk0_6[k]
                  + f_15 * hsi_3[k]
                  - f_12 * pc_y[k] * hsk1_6[k];

        t_43[k] = f_3 * pc_z[k] * isi_31[k];

        t_44[k] = f_13 * hsi_5[k]
                  + f_3 * pc_y[k] * isi_33[k];

        t_45[k] = pa_y[k] * hsk0_9[k]
                  - f_12 * pc_y[k] * hsk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, hsk0_10, hsi_6, hsi_9, \
                         hsk1_10, ish0_24, ish1_24, isi_34, isi_35, \
                         isi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * hsk0_10[k]
                  + f_16 * hsi_6[k]
                  - f_12 * pc_y[k] * hsk1_10[k];

        t_47[k] = f_3 * pc_z[k] * isi_34[k];

        t_48[k] = f_4 * ish0_24[k]
                  - f_5 * ish1_24[k]
                  + f_3 * pc_z[k] * isi_35[k];

        t_49[k] = f_13 * hsi_9[k]
                  + f_3 * pc_y[k] * isi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, hsk0_14, hsk0_15, hsi_10, \
                         hsk1_14, hsk1_15, ish0_27, ish1_27, isi_38, \
                         isi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * hsk0_14[k]
                  - f_12 * pc_y[k] * hsk1_14[k];

        t_51[k] = pa_y[k] * hsk0_15[k]
                  + f_17 * hsi_10[k]
                  - f_12 * pc_y[k] * hsk1_15[k];

        t_52[k] = f_3 * pc_z[k] * isi_38[k];

        t_53[k] = f_4 * ish0_27[k]
                  - f_5 * ish1_27[k]
                  + f_3 * pc_z[k] * isi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, hsk0_20, hsi_14, hsk1_20, \
                         ish0_28, ish1_28, isi_40, isi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * ish0_28[k]
                  - f_7 * ish1_28[k]
                  + f_3 * pc_z[k] * isi_40[k];

        t_55[k] = f_13 * hsi_14[k]
                  + f_3 * pc_y[k] * isi_42[k];

        t_56[k] = pa_y[k] * hsk0_20[k]
                  - f_12 * pc_y[k] * hsk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, hsi_49, hsi_51, hsi_52, \
                         hsi_53, isi_43, isi_49, isi_51, isi_52, \
                         isi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_17 * hsi_49[k]
                  + f_3 * pc_x[k] * isi_49[k];

        t_58[k] = f_3 * pc_z[k] * isi_43[k];

        t_59[k] = f_17 * hsi_51[k]
                  + f_3 * pc_x[k] * isi_51[k];

        t_60[k] = f_17 * hsi_52[k]
                  + f_3 * pc_x[k] * isi_52[k];

        t_61[k] = f_17 * hsi_53[k]
                  + f_3 * pc_x[k] * isi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, hsi_21, hsi_54, hsi_55, \
                         ish0_36, ish1_36, isi_49, isi_54, isi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_17 * hsi_54[k]
                  + f_3 * pc_x[k] * isi_54[k];

        t_63[k] = f_17 * hsi_55[k]
                  + f_3 * pc_x[k] * isi_55[k];

        t_64[k] = f_13 * hsi_21[k]
                  + f_1 * ish0_36[k]
                  - f_2 * ish1_36[k]
                  + f_3 * pc_y[k] * isi_49[k];

        t_65[k] = f_3 * pc_z[k] * isi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, ish0_36, ish0_37, ish0_38, ish1_36, ish1_37, \
                         ish1_38, isi_50, isi_51, isi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * ish0_36[k]
                  - f_5 * ish1_36[k]
                  + f_3 * pc_z[k] * isi_50[k];

        t_67[k] = f_6 * ish0_37[k]
                  - f_7 * ish1_37[k]
                  + f_3 * pc_z[k] * isi_51[k];

        t_68[k] = f_8 * ish0_38[k]
                  - f_9 * ish1_38[k]
                  + f_3 * pc_z[k] * isi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, hsk0_35, hsi_27, hsk1_35, \
                         ish0_39, ish1_39, isi_53, isi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ish0_39[k]
                  - f_11 * ish1_39[k]
                  + f_3 * pc_z[k] * isi_53[k];

        t_70[k] = f_13 * hsi_27[k]
                  + f_3 * pc_y[k] * isi_55[k];

        t_71[k] = pa_y[k] * hsk0_35[k]
                  - f_12 * pc_y[k] * hsk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, hsk0_0, hsk0_3, \
                         hsi_0, hsk1_0, hsk1_3, isi_56, isi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * hsk0_0[k]
                  - f_12 * pc_z[k] * hsk1_0[k];

        t_73[k] = f_3 * pc_y[k] * isi_56[k];

        t_74[k] = f_13 * hsi_0[k]
                  + f_3 * pc_z[k] * isi_56[k];

        t_75[k] = pa_z[k] * hsk0_3[k]
                  - f_12 * pc_z[k] * hsk1_3[k];

        t_76[k] = f_3 * pc_y[k] * isi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, hsk0_5, hsk0_6, hsi_2, \
                         hsk1_5, hsk1_6, ish0_44, ish1_44, isi_60, \
                         isi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * hsk0_5[k]
                  + f_14 * hsi_2[k]
                  - f_12 * pc_z[k] * hsk1_5[k];

        t_78[k] = pa_z[k] * hsk0_6[k]
                  - f_12 * pc_z[k] * hsk1_6[k];

        t_79[k] = f_4 * ish0_44[k]
                  - f_5 * ish1_44[k]
                  + f_3 * pc_y[k] * isi_60[k];

        t_80[k] = f_3 * pc_y[k] * isi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, hsk0_9, hsk0_10, hsi_5, hsk1_9, \
                         hsk1_10, ish0_46, ish1_46, isi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * hsk0_9[k]
                  + f_15 * hsi_5[k]
                  - f_12 * pc_z[k] * hsk1_9[k];

        t_82[k] = pa_z[k] * hsk0_10[k]
                  - f_12 * pc_z[k] * hsk1_10[k];

        t_83[k] = f_6 * ish0_46[k]
                  - f_7 * ish1_46[k]
                  + f_3 * pc_y[k] * isi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, hsk0_14, hsk0_15, hsi_9, \
                         hsk1_14, hsk1_15, ish0_47, ish1_47, isi_64, \
                         isi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * ish0_47[k]
                  - f_5 * ish1_47[k]
                  + f_3 * pc_y[k] * isi_64[k];

        t_85[k] = f_3 * pc_y[k] * isi_65[k];

        t_86[k] = pa_z[k] * hsk0_14[k]
                  + f_16 * hsi_9[k]
                  - f_12 * pc_z[k] * hsk1_14[k];

        t_87[k] = pa_z[k] * hsk0_15[k]
                  - f_12 * pc_z[k] * hsk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, ish0_49, ish0_50, ish0_51, ish1_49, \
                         ish1_50, ish1_51, isi_67, isi_68, isi_69, \
                         isi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * ish0_49[k]
                  - f_9 * ish1_49[k]
                  + f_3 * pc_y[k] * isi_67[k];

        t_89[k] = f_6 * ish0_50[k]
                  - f_7 * ish1_50[k]
                  + f_3 * pc_y[k] * isi_68[k];

        t_90[k] = f_4 * ish0_51[k]
                  - f_5 * ish1_51[k]
                  + f_3 * pc_y[k] * isi_69[k];

        t_91[k] = f_3 * pc_y[k] * isi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, hsk0_20, hsi_14, hsi_77, \
                         hsi_78, hsi_79, hsk1_20, isi_77, isi_78, \
                         isi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * hsk0_20[k]
                  + f_17 * hsi_14[k]
                  - f_12 * pc_z[k] * hsk1_20[k];

        t_93[k] = f_17 * hsi_77[k]
                  + f_3 * pc_x[k] * isi_77[k];

        t_94[k] = f_17 * hsi_78[k]
                  + f_3 * pc_x[k] * isi_78[k];

        t_95[k] = f_17 * hsi_79[k]
                  + f_3 * pc_x[k] * isi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, hsi_80, hsi_81, hsi_83, isi_76, \
                         isi_80, isi_81, isi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_17 * hsi_80[k]
                  + f_3 * pc_x[k] * isi_80[k];

        t_97[k] = f_17 * hsi_81[k]
                  + f_3 * pc_x[k] * isi_81[k];

        t_98[k] = f_3 * pc_y[k] * isi_76[k];

        t_99[k] = f_17 * hsi_83[k]
                  + f_3 * pc_x[k] * isi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, hsk0_28, hsk1_28, ish0_58, \
                         ish0_59, ish1_58, ish1_59, isi_78, isi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * hsk0_28[k]
                   - f_12 * pc_z[k] * hsk1_28[k];

        t_101[k] = f_18 * ish0_58[k]
                   - f_19 * ish1_58[k]
                   + f_3 * pc_y[k] * isi_78[k];

        t_102[k] = f_10 * ish0_59[k]
                   - f_11 * ish1_59[k]
                   + f_3 * pc_y[k] * isi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, ish0_60, ish0_61, ish0_62, ish1_60, \
                         ish1_61, ish1_62, isi_80, isi_81, isi_82, \
                         isi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * ish0_60[k]
                   - f_9 * ish1_60[k]
                   + f_3 * pc_y[k] * isi_80[k];

        t_104[k] = f_6 * ish0_61[k]
                   - f_7 * ish1_61[k]
                   + f_3 * pc_y[k] * isi_81[k];

        t_105[k] = f_4 * ish0_62[k]
                   - f_5 * ish1_62[k]
                   + f_3 * pc_y[k] * isi_82[k];

        t_106[k] = f_3 * pc_y[k] * isi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, hsi_27, hsi_28, hsi_84, \
                         ish0_62, ish0_63, ish1_62, ish1_63, isi_83, \
                         isi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * hsi_27[k]
                   + f_1 * ish0_62[k]
                   - f_2 * ish1_62[k]
                   + f_3 * pc_z[k] * isi_83[k];

        t_108[k] = f_16 * hsi_84[k]
                   + f_1 * ish0_63[k]
                   - f_2 * ish1_63[k]
                   + f_3 * pc_x[k] * isi_84[k];

        t_109[k] = f_14 * hsi_28[k]
                   + f_3 * pc_y[k] * isi_84[k];

        t_110[k] = f_3 * pc_z[k] * isi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, hsi_87, ish0_63, ish0_66, ish1_63, \
                         ish1_66, isi_85, isi_86, isi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_16 * hsi_87[k]
                   + f_10 * ish0_66[k]
                   - f_11 * ish1_66[k]
                   + f_3 * pc_x[k] * isi_87[k];

        t_112[k] = f_3 * pc_z[k] * isi_85[k];

        t_113[k] = f_4 * ish0_63[k]
                   - f_5 * ish1_63[k]
                   + f_3 * pc_z[k] * isi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, hsi_33, hsi_90, \
                         ish0_65, ish0_69, ish1_65, ish1_69, isi_87, isi_89, \
                         isi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_16 * hsi_90[k]
                   + f_8 * ish0_69[k]
                   - f_9 * ish1_69[k]
                   + f_3 * pc_x[k] * isi_90[k];

        t_115[k] = f_3 * pc_z[k] * isi_87[k];

        t_116[k] = f_14 * hsi_33[k]
                   + f_3 * pc_y[k] * isi_89[k];

        t_117[k] = f_6 * ish0_65[k]
                   - f_7 * ish1_65[k]
                   + f_3 * pc_z[k] * isi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, hsi_94, ish0_66, ish0_73, ish1_66, \
                         ish1_73, isi_90, isi_91, isi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_16 * hsi_94[k]
                   + f_6 * ish0_73[k]
                   - f_7 * ish1_73[k]
                   + f_3 * pc_x[k] * isi_94[k];

        t_119[k] = f_3 * pc_z[k] * isi_90[k];

        t_120[k] = f_4 * ish0_66[k]
                   - f_5 * ish1_66[k]
                   + f_3 * pc_z[k] * isi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, hsi_37, hsi_99, \
                         ish0_68, ish0_78, ish1_68, ish1_78, isi_93, isi_94, \
                         isi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * hsi_37[k]
                   + f_3 * pc_y[k] * isi_93[k];

        t_122[k] = f_8 * ish0_68[k]
                   - f_9 * ish1_68[k]
                   + f_3 * pc_z[k] * isi_93[k];

        t_123[k] = f_16 * hsi_99[k]
                   + f_4 * ish0_78[k]
                   - f_5 * ish1_78[k]
                   + f_3 * pc_x[k] * isi_99[k];

        t_124[k] = f_3 * pc_z[k] * isi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, hsi_42, ish0_69, ish0_70, \
                         ish0_72, ish1_69, ish1_70, ish1_72, isi_95, isi_96, \
                         isi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * ish0_69[k]
                   - f_5 * ish1_69[k]
                   + f_3 * pc_z[k] * isi_95[k];

        t_126[k] = f_6 * ish0_70[k]
                   - f_7 * ish1_70[k]
                   + f_3 * pc_z[k] * isi_96[k];

        t_127[k] = f_14 * hsi_42[k]
                   + f_3 * pc_y[k] * isi_98[k];

        t_128[k] = f_10 * ish0_72[k]
                   - f_11 * ish1_72[k]
                   + f_3 * pc_z[k] * isi_98[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
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
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);

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

    const auto *hsk0_39 = buffer.data(hsk0 + 39);
    const auto *hsk0_42 = buffer.data(hsk0 + 42);
    const auto *hsk0_46 = buffer.data(hsk0 + 46);
    const auto *hsk0_51 = buffer.data(hsk0 + 51);
    const auto *hsk0_64 = buffer.data(hsk0 + 64);
    const auto *hsk0_72 = buffer.data(hsk0 + 72);
    const auto *hsk0_77 = buffer.data(hsk0 + 77);
    const auto *hsk0_81 = buffer.data(hsk0 + 81);
    const auto *hsk0_84 = buffer.data(hsk0 + 84);
    const auto *hsk0_86 = buffer.data(hsk0 + 86);
    const auto *hsk0_89 = buffer.data(hsk0 + 89);
    const auto *hsk0_90 = buffer.data(hsk0 + 90);
    const auto *hsk0_92 = buffer.data(hsk0 + 92);
    const auto *hsk0_107 = buffer.data(hsk0 + 107);

    const auto *hsi_28 = buffer.data(hsi + 28);
    const auto *hsi_31 = buffer.data(hsi + 31);
    const auto *hsi_34 = buffer.data(hsi + 34);
    const auto *hsi_38 = buffer.data(hsi + 38);
    const auto *hsi_49 = buffer.data(hsi + 49);
    const auto *hsi_55 = buffer.data(hsi + 55);
    const auto *hsi_56 = buffer.data(hsi + 56);
    const auto *hsi_58 = buffer.data(hsi + 58);
    const auto *hsi_61 = buffer.data(hsi + 61);
    const auto *hsi_64 = buffer.data(hsi + 64);
    const auto *hsi_65 = buffer.data(hsi + 65);
    const auto *hsi_68 = buffer.data(hsi + 68);
    const auto *hsi_69 = buffer.data(hsi + 69);
    const auto *hsi_70 = buffer.data(hsi + 70);
    const auto *hsi_79 = buffer.data(hsi + 79);
    const auto *hsi_80 = buffer.data(hsi + 80);
    const auto *hsi_81 = buffer.data(hsi + 81);
    const auto *hsi_82 = buffer.data(hsi + 82);
    const auto *hsi_83 = buffer.data(hsi + 83);
    const auto *hsi_84 = buffer.data(hsi + 84);
    const auto *hsi_89 = buffer.data(hsi + 89);
    const auto *hsi_93 = buffer.data(hsi + 93);
    const auto *hsi_98 = buffer.data(hsi + 98);
    const auto *hsi_105 = buffer.data(hsi + 105);
    const auto *hsi_107 = buffer.data(hsi + 107);
    const auto *hsi_108 = buffer.data(hsi + 108);
    const auto *hsi_109 = buffer.data(hsi + 109);
    const auto *hsi_110 = buffer.data(hsi + 110);
    const auto *hsi_111 = buffer.data(hsi + 111);
    const auto *hsi_133 = buffer.data(hsi + 133);
    const auto *hsi_134 = buffer.data(hsi + 134);
    const auto *hsi_135 = buffer.data(hsi + 135);
    const auto *hsi_136 = buffer.data(hsi + 136);
    const auto *hsi_137 = buffer.data(hsi + 137);
    const auto *hsi_138 = buffer.data(hsi + 138);
    const auto *hsi_139 = buffer.data(hsi + 139);
    const auto *hsi_140 = buffer.data(hsi + 140);
    const auto *hsi_145 = buffer.data(hsi + 145);
    const auto *hsi_149 = buffer.data(hsi + 149);
    const auto *hsi_154 = buffer.data(hsi + 154);
    const auto *hsi_160 = buffer.data(hsi + 160);
    const auto *hsi_161 = buffer.data(hsi + 161);
    const auto *hsi_162 = buffer.data(hsi + 162);
    const auto *hsi_163 = buffer.data(hsi + 163);
    const auto *hsi_164 = buffer.data(hsi + 164);
    const auto *hsi_165 = buffer.data(hsi + 165);
    const auto *hsi_167 = buffer.data(hsi + 167);
    const auto *hsi_168 = buffer.data(hsi + 168);
    const auto *hsi_171 = buffer.data(hsi + 171);
    const auto *hsi_174 = buffer.data(hsi + 174);
    const auto *hsi_178 = buffer.data(hsi + 178);
    const auto *hsi_183 = buffer.data(hsi + 183);
    const auto *hsi_189 = buffer.data(hsi + 189);
    const auto *hsi_191 = buffer.data(hsi + 191);
    const auto *hsi_192 = buffer.data(hsi + 192);
    const auto *hsi_193 = buffer.data(hsi + 193);
    const auto *hsi_194 = buffer.data(hsi + 194);
    const auto *hsi_195 = buffer.data(hsi + 195);

    const auto *hsk1_39 = buffer.data(hsk1 + 39);
    const auto *hsk1_42 = buffer.data(hsk1 + 42);
    const auto *hsk1_46 = buffer.data(hsk1 + 46);
    const auto *hsk1_51 = buffer.data(hsk1 + 51);
    const auto *hsk1_64 = buffer.data(hsk1 + 64);
    const auto *hsk1_72 = buffer.data(hsk1 + 72);
    const auto *hsk1_77 = buffer.data(hsk1 + 77);
    const auto *hsk1_81 = buffer.data(hsk1 + 81);
    const auto *hsk1_84 = buffer.data(hsk1 + 84);
    const auto *hsk1_86 = buffer.data(hsk1 + 86);
    const auto *hsk1_89 = buffer.data(hsk1 + 89);
    const auto *hsk1_90 = buffer.data(hsk1 + 90);
    const auto *hsk1_92 = buffer.data(hsk1 + 92);
    const auto *hsk1_107 = buffer.data(hsk1 + 107);

    const auto *ish0_78 = buffer.data(ish0 + 78);
    const auto *ish0_79 = buffer.data(ish0 + 79);
    const auto *ish0_80 = buffer.data(ish0 + 80);
    const auto *ish0_81 = buffer.data(ish0 + 81);
    const auto *ish0_83 = buffer.data(ish0 + 83);
    const auto *ish0_101 = buffer.data(ish0 + 101);
    const auto *ish0_102 = buffer.data(ish0 + 102);
    const auto *ish0_103 = buffer.data(ish0 + 103);
    const auto *ish0_104 = buffer.data(ish0 + 104);
    const auto *ish0_105 = buffer.data(ish0 + 105);
    const auto *ish0_106 = buffer.data(ish0 + 106);
    const auto *ish0_107 = buffer.data(ish0 + 107);
    const auto *ish0_108 = buffer.data(ish0 + 108);
    const auto *ish0_109 = buffer.data(ish0 + 109);
    const auto *ish0_110 = buffer.data(ish0 + 110);
    const auto *ish0_111 = buffer.data(ish0 + 111);
    const auto *ish0_112 = buffer.data(ish0 + 112);
    const auto *ish0_113 = buffer.data(ish0 + 113);
    const auto *ish0_114 = buffer.data(ish0 + 114);
    const auto *ish0_119 = buffer.data(ish0 + 119);
    const auto *ish0_120 = buffer.data(ish0 + 120);
    const auto *ish0_121 = buffer.data(ish0 + 121);
    const auto *ish0_122 = buffer.data(ish0 + 122);
    const auto *ish0_123 = buffer.data(ish0 + 123);
    const auto *ish0_124 = buffer.data(ish0 + 124);
    const auto *ish0_125 = buffer.data(ish0 + 125);
    const auto *ish0_126 = buffer.data(ish0 + 126);
    const auto *ish0_128 = buffer.data(ish0 + 128);
    const auto *ish0_129 = buffer.data(ish0 + 129);
    const auto *ish0_131 = buffer.data(ish0 + 131);
    const auto *ish0_132 = buffer.data(ish0 + 132);
    const auto *ish0_133 = buffer.data(ish0 + 133);
    const auto *ish0_135 = buffer.data(ish0 + 135);
    const auto *ish0_136 = buffer.data(ish0 + 136);
    const auto *ish0_141 = buffer.data(ish0 + 141);
    const auto *ish0_142 = buffer.data(ish0 + 142);
    const auto *ish0_143 = buffer.data(ish0 + 143);
    const auto *ish0_144 = buffer.data(ish0 + 144);

    const auto *ish1_78 = buffer.data(ish1 + 78);
    const auto *ish1_79 = buffer.data(ish1 + 79);
    const auto *ish1_80 = buffer.data(ish1 + 80);
    const auto *ish1_81 = buffer.data(ish1 + 81);
    const auto *ish1_83 = buffer.data(ish1 + 83);
    const auto *ish1_101 = buffer.data(ish1 + 101);
    const auto *ish1_102 = buffer.data(ish1 + 102);
    const auto *ish1_103 = buffer.data(ish1 + 103);
    const auto *ish1_104 = buffer.data(ish1 + 104);
    const auto *ish1_105 = buffer.data(ish1 + 105);
    const auto *ish1_106 = buffer.data(ish1 + 106);
    const auto *ish1_107 = buffer.data(ish1 + 107);
    const auto *ish1_108 = buffer.data(ish1 + 108);
    const auto *ish1_109 = buffer.data(ish1 + 109);
    const auto *ish1_110 = buffer.data(ish1 + 110);
    const auto *ish1_111 = buffer.data(ish1 + 111);
    const auto *ish1_112 = buffer.data(ish1 + 112);
    const auto *ish1_113 = buffer.data(ish1 + 113);
    const auto *ish1_114 = buffer.data(ish1 + 114);
    const auto *ish1_119 = buffer.data(ish1 + 119);
    const auto *ish1_120 = buffer.data(ish1 + 120);
    const auto *ish1_121 = buffer.data(ish1 + 121);
    const auto *ish1_122 = buffer.data(ish1 + 122);
    const auto *ish1_123 = buffer.data(ish1 + 123);
    const auto *ish1_124 = buffer.data(ish1 + 124);
    const auto *ish1_125 = buffer.data(ish1 + 125);
    const auto *ish1_126 = buffer.data(ish1 + 126);
    const auto *ish1_128 = buffer.data(ish1 + 128);
    const auto *ish1_129 = buffer.data(ish1 + 129);
    const auto *ish1_131 = buffer.data(ish1 + 131);
    const auto *ish1_132 = buffer.data(ish1 + 132);
    const auto *ish1_133 = buffer.data(ish1 + 133);
    const auto *ish1_135 = buffer.data(ish1 + 135);
    const auto *ish1_136 = buffer.data(ish1 + 136);
    const auto *ish1_141 = buffer.data(ish1 + 141);
    const auto *ish1_142 = buffer.data(ish1 + 142);
    const auto *ish1_143 = buffer.data(ish1 + 143);
    const auto *ish1_144 = buffer.data(ish1 + 144);

    const auto *isi_99 = buffer.data(isi + 99);
    const auto *isi_105 = buffer.data(isi + 105);
    const auto *isi_106 = buffer.data(isi + 106);
    const auto *isi_107 = buffer.data(isi + 107);
    const auto *isi_108 = buffer.data(isi + 108);
    const auto *isi_109 = buffer.data(isi + 109);
    const auto *isi_110 = buffer.data(isi + 110);
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
    const auto *isi_134 = buffer.data(isi + 134);
    const auto *isi_135 = buffer.data(isi + 135);
    const auto *isi_136 = buffer.data(isi + 136);
    const auto *isi_137 = buffer.data(isi + 137);
    const auto *isi_138 = buffer.data(isi + 138);
    const auto *isi_139 = buffer.data(isi + 139);
    const auto *isi_140 = buffer.data(isi + 140);
    const auto *isi_141 = buffer.data(isi + 141);
    const auto *isi_142 = buffer.data(isi + 142);
    const auto *isi_143 = buffer.data(isi + 143);
    const auto *isi_144 = buffer.data(isi + 144);
    const auto *isi_145 = buffer.data(isi + 145);
    const auto *isi_146 = buffer.data(isi + 146);
    const auto *isi_147 = buffer.data(isi + 147);
    const auto *isi_148 = buffer.data(isi + 148);
    const auto *isi_149 = buffer.data(isi + 149);
    const auto *isi_150 = buffer.data(isi + 150);
    const auto *isi_151 = buffer.data(isi + 151);
    const auto *isi_152 = buffer.data(isi + 152);
    const auto *isi_153 = buffer.data(isi + 153);
    const auto *isi_154 = buffer.data(isi + 154);
    const auto *isi_160 = buffer.data(isi + 160);
    const auto *isi_161 = buffer.data(isi + 161);
    const auto *isi_162 = buffer.data(isi + 162);
    const auto *isi_163 = buffer.data(isi + 163);
    const auto *isi_164 = buffer.data(isi + 164);
    const auto *isi_165 = buffer.data(isi + 165);
    const auto *isi_166 = buffer.data(isi + 166);
    const auto *isi_167 = buffer.data(isi + 167);
    const auto *isi_168 = buffer.data(isi + 168);
    const auto *isi_169 = buffer.data(isi + 169);
    const auto *isi_170 = buffer.data(isi + 170);
    const auto *isi_171 = buffer.data(isi + 171);
    const auto *isi_173 = buffer.data(isi + 173);
    const auto *isi_174 = buffer.data(isi + 174);
    const auto *isi_175 = buffer.data(isi + 175);
    const auto *isi_177 = buffer.data(isi + 177);
    const auto *isi_178 = buffer.data(isi + 178);
    const auto *isi_179 = buffer.data(isi + 179);
    const auto *isi_180 = buffer.data(isi + 180);
    const auto *isi_182 = buffer.data(isi + 182);
    const auto *isi_183 = buffer.data(isi + 183);
    const auto *isi_189 = buffer.data(isi + 189);
    const auto *isi_190 = buffer.data(isi + 190);
    const auto *isi_191 = buffer.data(isi + 191);
    const auto *isi_192 = buffer.data(isi + 192);
    const auto *isi_193 = buffer.data(isi + 193);
    const auto *isi_194 = buffer.data(isi + 194);
    const auto *isi_195 = buffer.data(isi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, hsi_105, hsi_107, \
                         hsi_108, hsi_109, isi_99, isi_105, isi_107, isi_108, \
                         isi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_16 * hsi_105[k]
                   + f_3 * pc_x[k] * isi_105[k];

        t_130[k] = f_3 * pc_z[k] * isi_99[k];

        t_131[k] = f_16 * hsi_107[k]
                   + f_3 * pc_x[k] * isi_107[k];

        t_132[k] = f_16 * hsi_108[k]
                   + f_3 * pc_x[k] * isi_108[k];

        t_133[k] = f_16 * hsi_109[k]
                   + f_3 * pc_x[k] * isi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, hsi_49, hsi_110, \
                         hsi_111, ish0_78, ish1_78, isi_105, isi_110, \
                         isi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_16 * hsi_110[k]
                   + f_3 * pc_x[k] * isi_110[k];

        t_135[k] = f_16 * hsi_111[k]
                   + f_3 * pc_x[k] * isi_111[k];

        t_136[k] = f_14 * hsi_49[k]
                   + f_1 * ish0_78[k]
                   - f_2 * ish1_78[k]
                   + f_3 * pc_y[k] * isi_105[k];

        t_137[k] = f_3 * pc_z[k] * isi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, ish0_78, ish0_79, ish0_80, ish1_78, \
                         ish1_79, ish1_80, isi_106, isi_107, isi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * ish0_78[k]
                   - f_5 * ish1_78[k]
                   + f_3 * pc_z[k] * isi_106[k];

        t_139[k] = f_6 * ish0_79[k]
                   - f_7 * ish1_79[k]
                   + f_3 * pc_z[k] * isi_107[k];

        t_140[k] = f_8 * ish0_80[k]
                   - f_9 * ish1_80[k]
                   + f_3 * pc_z[k] * isi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, hsk0_72, hsi_55, \
                         hsk1_72, ish0_81, ish0_83, ish1_81, ish1_83, isi_109, \
                         isi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * ish0_81[k]
                   - f_11 * ish1_81[k]
                   + f_3 * pc_z[k] * isi_109[k];

        t_142[k] = f_14 * hsi_55[k]
                   + f_3 * pc_y[k] * isi_111[k];

        t_143[k] = f_1 * ish0_83[k]
                   - f_2 * ish1_83[k]
                   + f_3 * pc_z[k] * isi_111[k];

        t_144[k] = pa_y[k] * hsk0_72[k]
                   - f_12 * pc_y[k] * hsk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, hsk0_39, hsi_28, \
                         hsi_56, hsi_58, hsk1_39, isi_112, isi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * hsi_56[k]
                   + f_3 * pc_y[k] * isi_112[k];

        t_146[k] = f_13 * hsi_28[k]
                   + f_3 * pc_z[k] * isi_112[k];

        t_147[k] = pa_z[k] * hsk0_39[k]
                   - f_12 * pc_z[k] * hsk1_39[k];

        t_148[k] = f_13 * hsi_58[k]
                   + f_3 * pc_y[k] * isi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, hsk0_42, hsk0_77, \
                         hsi_31, hsi_61, hsk1_42, hsk1_77, isi_115, \
                         isi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * hsk0_77[k]
                   - f_12 * pc_y[k] * hsk1_77[k];

        t_150[k] = pa_z[k] * hsk0_42[k]
                   - f_12 * pc_z[k] * hsk1_42[k];

        t_151[k] = f_13 * hsi_31[k]
                   + f_3 * pc_z[k] * isi_115[k];

        t_152[k] = f_13 * hsi_61[k]
                   + f_3 * pc_y[k] * isi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, hsk0_46, hsk0_81, \
                         hsi_34, hsk1_46, hsk1_81, isi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * hsk0_81[k]
                   - f_12 * pc_y[k] * hsk1_81[k];

        t_154[k] = pa_z[k] * hsk0_46[k]
                   - f_12 * pc_z[k] * hsk1_46[k];

        t_155[k] = f_13 * hsi_34[k]
                   + f_3 * pc_z[k] * isi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, hsk0_84, hsk0_86, hsi_64, hsi_65, \
                         hsk1_84, hsk1_86, isi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * hsk0_84[k]
                   + f_14 * hsi_64[k]
                   - f_12 * pc_y[k] * hsk1_84[k];

        t_157[k] = f_13 * hsi_65[k]
                   + f_3 * pc_y[k] * isi_121[k];

        t_158[k] = pa_y[k] * hsk0_86[k]
                   - f_12 * pc_y[k] * hsk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, hsk0_51, hsk0_89, \
                         hsi_38, hsi_68, hsk1_51, hsk1_89, isi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * hsk0_51[k]
                   - f_12 * pc_z[k] * hsk1_51[k];

        t_160[k] = f_13 * hsi_38[k]
                   + f_3 * pc_z[k] * isi_122[k];

        t_161[k] = pa_y[k] * hsk0_89[k]
                   + f_15 * hsi_68[k]
                   - f_12 * pc_y[k] * hsk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, hsk0_90, hsk0_92, \
                         hsi_69, hsi_70, hsi_133, hsk1_90, hsk1_92, isi_126, \
                         isi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * hsk0_90[k]
                   + f_14 * hsi_69[k]
                   - f_12 * pc_y[k] * hsk1_90[k];

        t_163[k] = f_13 * hsi_70[k]
                   + f_3 * pc_y[k] * isi_126[k];

        t_164[k] = pa_y[k] * hsk0_92[k]
                   - f_12 * pc_y[k] * hsk1_92[k];

        t_165[k] = f_16 * hsi_133[k]
                   + f_3 * pc_x[k] * isi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, hsi_134, hsi_135, hsi_136, \
                         hsi_137, hsi_138, isi_134, isi_135, isi_136, isi_137, \
                         isi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_16 * hsi_134[k]
                   + f_3 * pc_x[k] * isi_134[k];

        t_167[k] = f_16 * hsi_135[k]
                   + f_3 * pc_x[k] * isi_135[k];

        t_168[k] = f_16 * hsi_136[k]
                   + f_3 * pc_x[k] * isi_136[k];

        t_169[k] = f_16 * hsi_137[k]
                   + f_3 * pc_x[k] * isi_137[k];

        t_170[k] = f_16 * hsi_138[k]
                   + f_3 * pc_x[k] * isi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, hsk0_64, hsi_49, hsi_139, \
                         hsk1_64, isi_133, isi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * hsi_139[k]
                   + f_3 * pc_x[k] * isi_139[k];

        t_172[k] = pa_z[k] * hsk0_64[k]
                   - f_12 * pc_z[k] * hsk1_64[k];

        t_173[k] = f_13 * hsi_49[k]
                   + f_3 * pc_z[k] * isi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, hsi_79, hsi_80, hsi_81, ish0_101, \
                         ish0_102, ish0_103, ish1_101, ish1_102, ish1_103, isi_135, isi_136, \
                         isi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * hsi_79[k]
                   + f_10 * ish0_101[k]
                   - f_11 * ish1_101[k]
                   + f_3 * pc_y[k] * isi_135[k];

        t_175[k] = f_13 * hsi_80[k]
                   + f_8 * ish0_102[k]
                   - f_9 * ish1_102[k]
                   + f_3 * pc_y[k] * isi_136[k];

        t_176[k] = f_13 * hsi_81[k]
                   + f_6 * ish0_103[k]
                   - f_7 * ish1_103[k]
                   + f_3 * pc_y[k] * isi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, hsk0_107, hsi_82, hsi_83, hsk1_107, \
                         ish0_104, ish1_104, isi_138, isi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * hsi_82[k]
                   + f_4 * ish0_104[k]
                   - f_5 * ish1_104[k]
                   + f_3 * pc_y[k] * isi_138[k];

        t_178[k] = f_13 * hsi_83[k]
                   + f_3 * pc_y[k] * isi_139[k];

        t_179[k] = pa_y[k] * hsk0_107[k]
                   - f_12 * pc_y[k] * hsk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, hsi_56, hsi_140, \
                         ish0_105, ish1_105, isi_140, isi_141, \
                         isi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_16 * hsi_140[k]
                   + f_1 * ish0_105[k]
                   - f_2 * ish1_105[k]
                   + f_3 * pc_x[k] * isi_140[k];

        t_181[k] = f_3 * pc_y[k] * isi_140[k];

        t_182[k] = f_14 * hsi_56[k]
                   + f_3 * pc_z[k] * isi_140[k];

        t_183[k] = f_4 * ish0_105[k]
                   - f_5 * ish1_105[k]
                   + f_3 * pc_y[k] * isi_141[k];

        t_184[k] = f_3 * pc_y[k] * isi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, hsi_145, ish0_106, ish0_107, \
                         ish0_110, ish1_106, ish1_107, ish1_110, isi_143, isi_144, \
                         isi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_16 * hsi_145[k]
                   + f_10 * ish0_110[k]
                   - f_11 * ish1_110[k]
                   + f_3 * pc_x[k] * isi_145[k];

        t_186[k] = f_6 * ish0_106[k]
                   - f_7 * ish1_106[k]
                   + f_3 * pc_y[k] * isi_143[k];

        t_187[k] = f_4 * ish0_107[k]
                   - f_5 * ish1_107[k]
                   + f_3 * pc_y[k] * isi_144[k];

        t_188[k] = f_3 * pc_y[k] * isi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, hsi_149, ish0_108, ish0_109, \
                         ish0_114, ish1_108, ish1_109, ish1_114, isi_146, isi_147, \
                         isi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * hsi_149[k]
                   + f_8 * ish0_114[k]
                   - f_9 * ish1_114[k]
                   + f_3 * pc_x[k] * isi_149[k];

        t_190[k] = f_8 * ish0_108[k]
                   - f_9 * ish1_108[k]
                   + f_3 * pc_y[k] * isi_146[k];

        t_191[k] = f_6 * ish0_109[k]
                   - f_7 * ish1_109[k]
                   + f_3 * pc_y[k] * isi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, hsi_154, ish0_110, ish0_119, \
                         ish1_110, ish1_119, isi_148, isi_149, \
                         isi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * ish0_110[k]
                   - f_5 * ish1_110[k]
                   + f_3 * pc_y[k] * isi_148[k];

        t_193[k] = f_3 * pc_y[k] * isi_149[k];

        t_194[k] = f_16 * hsi_154[k]
                   + f_6 * ish0_119[k]
                   - f_7 * ish1_119[k]
                   + f_3 * pc_x[k] * isi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, ish0_111, ish0_112, ish0_113, ish1_111, \
                         ish1_112, ish1_113, isi_150, isi_151, \
                         isi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * ish0_111[k]
                   - f_11 * ish1_111[k]
                   + f_3 * pc_y[k] * isi_150[k];

        t_196[k] = f_8 * ish0_112[k]
                   - f_9 * ish1_112[k]
                   + f_3 * pc_y[k] * isi_151[k];

        t_197[k] = f_6 * ish0_113[k]
                   - f_7 * ish1_113[k]
                   + f_3 * pc_y[k] * isi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, hsi_160, hsi_161, ish0_114, \
                         ish0_125, ish1_114, ish1_125, isi_153, isi_154, isi_160, \
                         isi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * ish0_114[k]
                   - f_5 * ish1_114[k]
                   + f_3 * pc_y[k] * isi_153[k];

        t_199[k] = f_3 * pc_y[k] * isi_154[k];

        t_200[k] = f_16 * hsi_160[k]
                   + f_4 * ish0_125[k]
                   - f_5 * ish1_125[k]
                   + f_3 * pc_x[k] * isi_160[k];

        t_201[k] = f_16 * hsi_161[k]
                   + f_3 * pc_x[k] * isi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, hsi_162, hsi_163, \
                         hsi_164, hsi_165, isi_160, isi_162, isi_163, isi_164, \
                         isi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_16 * hsi_162[k]
                   + f_3 * pc_x[k] * isi_162[k];

        t_203[k] = f_16 * hsi_163[k]
                   + f_3 * pc_x[k] * isi_163[k];

        t_204[k] = f_16 * hsi_164[k]
                   + f_3 * pc_x[k] * isi_164[k];

        t_205[k] = f_16 * hsi_165[k]
                   + f_3 * pc_x[k] * isi_165[k];

        t_206[k] = f_3 * pc_y[k] * isi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, hsi_167, ish0_120, ish0_121, \
                         ish1_120, ish1_121, isi_161, isi_162, \
                         isi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_16 * hsi_167[k]
                   + f_3 * pc_x[k] * isi_167[k];

        t_208[k] = f_1 * ish0_120[k]
                   - f_2 * ish1_120[k]
                   + f_3 * pc_y[k] * isi_161[k];

        t_209[k] = f_18 * ish0_121[k]
                   - f_19 * ish1_121[k]
                   + f_3 * pc_y[k] * isi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, ish0_122, ish0_123, ish0_124, ish1_122, \
                         ish1_123, ish1_124, isi_163, isi_164, \
                         isi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * ish0_122[k]
                   - f_11 * ish1_122[k]
                   + f_3 * pc_y[k] * isi_163[k];

        t_211[k] = f_8 * ish0_123[k]
                   - f_9 * ish1_123[k]
                   + f_3 * pc_y[k] * isi_164[k];

        t_212[k] = f_6 * ish0_124[k]
                   - f_7 * ish1_124[k]
                   + f_3 * pc_y[k] * isi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, hsi_83, hsi_168, \
                         ish0_125, ish0_126, ish1_125, ish1_126, isi_166, isi_167, \
                         isi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * ish0_125[k]
                   - f_5 * ish1_125[k]
                   + f_3 * pc_y[k] * isi_166[k];

        t_214[k] = f_3 * pc_y[k] * isi_167[k];

        t_215[k] = f_14 * hsi_83[k]
                   + f_1 * ish0_125[k]
                   - f_2 * ish1_125[k]
                   + f_3 * pc_z[k] * isi_167[k];

        t_216[k] = f_15 * hsi_168[k]
                   + f_1 * ish0_126[k]
                   - f_2 * ish1_126[k]
                   + f_3 * pc_x[k] * isi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, hsi_84, hsi_171, \
                         ish0_129, ish1_129, isi_168, isi_169, \
                         isi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * hsi_84[k]
                   + f_3 * pc_y[k] * isi_168[k];

        t_218[k] = f_3 * pc_z[k] * isi_168[k];

        t_219[k] = f_15 * hsi_171[k]
                   + f_10 * ish0_129[k]
                   - f_11 * ish1_129[k]
                   + f_3 * pc_x[k] * isi_171[k];

        t_220[k] = f_3 * pc_z[k] * isi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pc_x, pc_z, hsi_174, ish0_126, ish0_132, \
                         ish1_126, ish1_132, isi_170, isi_171, \
                         isi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * ish0_126[k]
                   - f_5 * ish1_126[k]
                   + f_3 * pc_z[k] * isi_170[k];

        t_222[k] = f_15 * hsi_174[k]
                   + f_8 * ish0_132[k]
                   - f_9 * ish1_132[k]
                   + f_3 * pc_x[k] * isi_174[k];

        t_223[k] = f_3 * pc_z[k] * isi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, hsi_89, hsi_178, \
                         ish0_128, ish0_136, ish1_128, ish1_136, isi_173, isi_174, \
                         isi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * hsi_89[k]
                   + f_3 * pc_y[k] * isi_173[k];

        t_225[k] = f_6 * ish0_128[k]
                   - f_7 * ish1_128[k]
                   + f_3 * pc_z[k] * isi_173[k];

        t_226[k] = f_15 * hsi_178[k]
                   + f_6 * ish0_136[k]
                   - f_7 * ish1_136[k]
                   + f_3 * pc_x[k] * isi_178[k];

        t_227[k] = f_3 * pc_z[k] * isi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, hsi_93, ish0_129, ish0_131, \
                         ish1_129, ish1_131, isi_175, isi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * ish0_129[k]
                   - f_5 * ish1_129[k]
                   + f_3 * pc_z[k] * isi_175[k];

        t_229[k] = f_15 * hsi_93[k]
                   + f_3 * pc_y[k] * isi_177[k];

        t_230[k] = f_8 * ish0_131[k]
                   - f_9 * ish1_131[k]
                   + f_3 * pc_z[k] * isi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_x, pc_z, hsi_183, ish0_132, ish0_141, \
                         ish1_132, ish1_141, isi_178, isi_179, \
                         isi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_15 * hsi_183[k]
                   + f_4 * ish0_141[k]
                   - f_5 * ish1_141[k]
                   + f_3 * pc_x[k] * isi_183[k];

        t_232[k] = f_3 * pc_z[k] * isi_178[k];

        t_233[k] = f_4 * ish0_132[k]
                   - f_5 * ish1_132[k]
                   + f_3 * pc_z[k] * isi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, hsi_98, hsi_189, \
                         ish0_133, ish0_135, ish1_133, ish1_135, isi_180, isi_182, \
                         isi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * ish0_133[k]
                   - f_7 * ish1_133[k]
                   + f_3 * pc_z[k] * isi_180[k];

        t_235[k] = f_15 * hsi_98[k]
                   + f_3 * pc_y[k] * isi_182[k];

        t_236[k] = f_10 * ish0_135[k]
                   - f_11 * ish1_135[k]
                   + f_3 * pc_z[k] * isi_182[k];

        t_237[k] = f_15 * hsi_189[k]
                   + f_3 * pc_x[k] * isi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, hsi_191, hsi_192, \
                         hsi_193, hsi_194, isi_183, isi_191, isi_192, isi_193, \
                         isi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * isi_183[k];

        t_239[k] = f_15 * hsi_191[k]
                   + f_3 * pc_x[k] * isi_191[k];

        t_240[k] = f_15 * hsi_192[k]
                   + f_3 * pc_x[k] * isi_192[k];

        t_241[k] = f_15 * hsi_193[k]
                   + f_3 * pc_x[k] * isi_193[k];

        t_242[k] = f_15 * hsi_194[k]
                   + f_3 * pc_x[k] * isi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, hsi_105, hsi_195, \
                         ish0_141, ish1_141, isi_189, isi_190, \
                         isi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_15 * hsi_195[k]
                   + f_3 * pc_x[k] * isi_195[k];

        t_244[k] = f_15 * hsi_105[k]
                   + f_1 * ish0_141[k]
                   - f_2 * ish1_141[k]
                   + f_3 * pc_y[k] * isi_189[k];

        t_245[k] = f_3 * pc_z[k] * isi_189[k];

        t_246[k] = f_4 * ish0_141[k]
                   - f_5 * ish1_141[k]
                   + f_3 * pc_z[k] * isi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, ish0_142, ish0_143, ish0_144, ish1_142, \
                         ish1_143, ish1_144, isi_191, isi_192, \
                         isi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * ish0_142[k]
                   - f_7 * ish1_142[k]
                   + f_3 * pc_z[k] * isi_191[k];

        t_248[k] = f_8 * ish0_143[k]
                   - f_9 * ish1_143[k]
                   + f_3 * pc_z[k] * isi_192[k];

        t_249[k] = f_10 * ish0_144[k]
                   - f_11 * ish1_144[k]
                   + f_3 * pc_z[k] * isi_193[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
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
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);

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

    const auto *hsk0_108 = buffer.data(hsk0 + 108);
    const auto *hsk0_111 = buffer.data(hsk0 + 111);
    const auto *hsk0_114 = buffer.data(hsk0 + 114);
    const auto *hsk0_118 = buffer.data(hsk0 + 118);
    const auto *hsk0_120 = buffer.data(hsk0 + 120);
    const auto *hsk0_123 = buffer.data(hsk0 + 123);
    const auto *hsk0_125 = buffer.data(hsk0 + 125);
    const auto *hsk0_126 = buffer.data(hsk0 + 126);
    const auto *hsk0_136 = buffer.data(hsk0 + 136);
    const auto *hsk0_180 = buffer.data(hsk0 + 180);
    const auto *hsk0_183 = buffer.data(hsk0 + 183);
    const auto *hsk0_185 = buffer.data(hsk0 + 185);
    const auto *hsk0_186 = buffer.data(hsk0 + 186);
    const auto *hsk0_189 = buffer.data(hsk0 + 189);
    const auto *hsk0_190 = buffer.data(hsk0 + 190);
    const auto *hsk0_192 = buffer.data(hsk0 + 192);
    const auto *hsk0_194 = buffer.data(hsk0 + 194);
    const auto *hsk0_195 = buffer.data(hsk0 + 195);
    const auto *hsk0_197 = buffer.data(hsk0 + 197);
    const auto *hsk0_198 = buffer.data(hsk0 + 198);
    const auto *hsk0_200 = buffer.data(hsk0 + 200);
    const auto *hsk0_215 = buffer.data(hsk0 + 215);

    const auto *hsi_84 = buffer.data(hsi + 84);
    const auto *hsi_87 = buffer.data(hsi + 87);
    const auto *hsi_90 = buffer.data(hsi + 90);
    const auto *hsi_91 = buffer.data(hsi + 91);
    const auto *hsi_94 = buffer.data(hsi + 94);
    const auto *hsi_95 = buffer.data(hsi + 95);
    const auto *hsi_96 = buffer.data(hsi + 96);
    const auto *hsi_105 = buffer.data(hsi + 105);
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
    const auto *hsi_135 = buffer.data(hsi + 135);
    const auto *hsi_136 = buffer.data(hsi + 136);
    const auto *hsi_137 = buffer.data(hsi + 137);
    const auto *hsi_138 = buffer.data(hsi + 138);
    const auto *hsi_139 = buffer.data(hsi + 139);
    const auto *hsi_140 = buffer.data(hsi + 140);
    const auto *hsi_141 = buffer.data(hsi + 141);
    const auto *hsi_142 = buffer.data(hsi + 142);
    const auto *hsi_143 = buffer.data(hsi + 143);
    const auto *hsi_145 = buffer.data(hsi + 145);
    const auto *hsi_146 = buffer.data(hsi + 146);
    const auto *hsi_148 = buffer.data(hsi + 148);
    const auto *hsi_149 = buffer.data(hsi + 149);
    const auto *hsi_150 = buffer.data(hsi + 150);
    const auto *hsi_152 = buffer.data(hsi + 152);
    const auto *hsi_153 = buffer.data(hsi + 153);
    const auto *hsi_154 = buffer.data(hsi + 154);
    const auto *hsi_161 = buffer.data(hsi + 161);
    const auto *hsi_163 = buffer.data(hsi + 163);
    const auto *hsi_164 = buffer.data(hsi + 164);
    const auto *hsi_165 = buffer.data(hsi + 165);
    const auto *hsi_166 = buffer.data(hsi + 166);
    const auto *hsi_167 = buffer.data(hsi + 167);
    const auto *hsi_168 = buffer.data(hsi + 168);
    const auto *hsi_201 = buffer.data(hsi + 201);
    const auto *hsi_205 = buffer.data(hsi + 205);
    const auto *hsi_210 = buffer.data(hsi + 210);
    const auto *hsi_216 = buffer.data(hsi + 216);
    const auto *hsi_217 = buffer.data(hsi + 217);
    const auto *hsi_218 = buffer.data(hsi + 218);
    const auto *hsi_219 = buffer.data(hsi + 219);
    const auto *hsi_220 = buffer.data(hsi + 220);
    const auto *hsi_221 = buffer.data(hsi + 221);
    const auto *hsi_222 = buffer.data(hsi + 222);
    const auto *hsi_223 = buffer.data(hsi + 223);
    const auto *hsi_245 = buffer.data(hsi + 245);
    const auto *hsi_246 = buffer.data(hsi + 246);
    const auto *hsi_247 = buffer.data(hsi + 247);
    const auto *hsi_248 = buffer.data(hsi + 248);
    const auto *hsi_249 = buffer.data(hsi + 249);
    const auto *hsi_250 = buffer.data(hsi + 250);
    const auto *hsi_251 = buffer.data(hsi + 251);
    const auto *hsi_252 = buffer.data(hsi + 252);
    const auto *hsi_257 = buffer.data(hsi + 257);
    const auto *hsi_261 = buffer.data(hsi + 261);
    const auto *hsi_266 = buffer.data(hsi + 266);
    const auto *hsi_272 = buffer.data(hsi + 272);
    const auto *hsi_273 = buffer.data(hsi + 273);
    const auto *hsi_274 = buffer.data(hsi + 274);
    const auto *hsi_275 = buffer.data(hsi + 275);
    const auto *hsi_276 = buffer.data(hsi + 276);
    const auto *hsi_277 = buffer.data(hsi + 277);
    const auto *hsi_279 = buffer.data(hsi + 279);
    const auto *hsi_280 = buffer.data(hsi + 280);
    const auto *hsi_283 = buffer.data(hsi + 283);
    const auto *hsi_286 = buffer.data(hsi + 286);

    const auto *hsk1_108 = buffer.data(hsk1 + 108);
    const auto *hsk1_111 = buffer.data(hsk1 + 111);
    const auto *hsk1_114 = buffer.data(hsk1 + 114);
    const auto *hsk1_118 = buffer.data(hsk1 + 118);
    const auto *hsk1_120 = buffer.data(hsk1 + 120);
    const auto *hsk1_123 = buffer.data(hsk1 + 123);
    const auto *hsk1_125 = buffer.data(hsk1 + 125);
    const auto *hsk1_126 = buffer.data(hsk1 + 126);
    const auto *hsk1_136 = buffer.data(hsk1 + 136);
    const auto *hsk1_180 = buffer.data(hsk1 + 180);
    const auto *hsk1_183 = buffer.data(hsk1 + 183);
    const auto *hsk1_185 = buffer.data(hsk1 + 185);
    const auto *hsk1_186 = buffer.data(hsk1 + 186);
    const auto *hsk1_189 = buffer.data(hsk1 + 189);
    const auto *hsk1_190 = buffer.data(hsk1 + 190);
    const auto *hsk1_192 = buffer.data(hsk1 + 192);
    const auto *hsk1_194 = buffer.data(hsk1 + 194);
    const auto *hsk1_195 = buffer.data(hsk1 + 195);
    const auto *hsk1_197 = buffer.data(hsk1 + 197);
    const auto *hsk1_198 = buffer.data(hsk1 + 198);
    const auto *hsk1_200 = buffer.data(hsk1 + 200);
    const auto *hsk1_215 = buffer.data(hsk1 + 215);

    const auto *ish0_146 = buffer.data(ish0 + 146);
    const auto *ish0_152 = buffer.data(ish0 + 152);
    const auto *ish0_156 = buffer.data(ish0 + 156);
    const auto *ish0_161 = buffer.data(ish0 + 161);
    const auto *ish0_164 = buffer.data(ish0 + 164);
    const auto *ish0_165 = buffer.data(ish0 + 165);
    const auto *ish0_166 = buffer.data(ish0 + 166);
    const auto *ish0_167 = buffer.data(ish0 + 167);
    const auto *ish0_183 = buffer.data(ish0 + 183);
    const auto *ish0_185 = buffer.data(ish0 + 185);
    const auto *ish0_186 = buffer.data(ish0 + 186);
    const auto *ish0_187 = buffer.data(ish0 + 187);
    const auto *ish0_188 = buffer.data(ish0 + 188);
    const auto *ish0_189 = buffer.data(ish0 + 189);
    const auto *ish0_190 = buffer.data(ish0 + 190);
    const auto *ish0_191 = buffer.data(ish0 + 191);
    const auto *ish0_192 = buffer.data(ish0 + 192);
    const auto *ish0_193 = buffer.data(ish0 + 193);
    const auto *ish0_194 = buffer.data(ish0 + 194);
    const auto *ish0_195 = buffer.data(ish0 + 195);
    const auto *ish0_196 = buffer.data(ish0 + 196);
    const auto *ish0_197 = buffer.data(ish0 + 197);
    const auto *ish0_198 = buffer.data(ish0 + 198);
    const auto *ish0_203 = buffer.data(ish0 + 203);
    const auto *ish0_204 = buffer.data(ish0 + 204);
    const auto *ish0_205 = buffer.data(ish0 + 205);
    const auto *ish0_206 = buffer.data(ish0 + 206);
    const auto *ish0_207 = buffer.data(ish0 + 207);
    const auto *ish0_208 = buffer.data(ish0 + 208);
    const auto *ish0_209 = buffer.data(ish0 + 209);
    const auto *ish0_210 = buffer.data(ish0 + 210);
    const auto *ish0_213 = buffer.data(ish0 + 213);
    const auto *ish0_216 = buffer.data(ish0 + 216);

    const auto *ish1_146 = buffer.data(ish1 + 146);
    const auto *ish1_152 = buffer.data(ish1 + 152);
    const auto *ish1_156 = buffer.data(ish1 + 156);
    const auto *ish1_161 = buffer.data(ish1 + 161);
    const auto *ish1_164 = buffer.data(ish1 + 164);
    const auto *ish1_165 = buffer.data(ish1 + 165);
    const auto *ish1_166 = buffer.data(ish1 + 166);
    const auto *ish1_167 = buffer.data(ish1 + 167);
    const auto *ish1_183 = buffer.data(ish1 + 183);
    const auto *ish1_185 = buffer.data(ish1 + 185);
    const auto *ish1_186 = buffer.data(ish1 + 186);
    const auto *ish1_187 = buffer.data(ish1 + 187);
    const auto *ish1_188 = buffer.data(ish1 + 188);
    const auto *ish1_189 = buffer.data(ish1 + 189);
    const auto *ish1_190 = buffer.data(ish1 + 190);
    const auto *ish1_191 = buffer.data(ish1 + 191);
    const auto *ish1_192 = buffer.data(ish1 + 192);
    const auto *ish1_193 = buffer.data(ish1 + 193);
    const auto *ish1_194 = buffer.data(ish1 + 194);
    const auto *ish1_195 = buffer.data(ish1 + 195);
    const auto *ish1_196 = buffer.data(ish1 + 196);
    const auto *ish1_197 = buffer.data(ish1 + 197);
    const auto *ish1_198 = buffer.data(ish1 + 198);
    const auto *ish1_203 = buffer.data(ish1 + 203);
    const auto *ish1_204 = buffer.data(ish1 + 204);
    const auto *ish1_205 = buffer.data(ish1 + 205);
    const auto *ish1_206 = buffer.data(ish1 + 206);
    const auto *ish1_207 = buffer.data(ish1 + 207);
    const auto *ish1_208 = buffer.data(ish1 + 208);
    const auto *ish1_209 = buffer.data(ish1 + 209);
    const auto *ish1_210 = buffer.data(ish1 + 210);
    const auto *ish1_213 = buffer.data(ish1 + 213);
    const auto *ish1_216 = buffer.data(ish1 + 216);

    const auto *isi_195 = buffer.data(isi + 195);
    const auto *isi_196 = buffer.data(isi + 196);
    const auto *isi_198 = buffer.data(isi + 198);
    const auto *isi_199 = buffer.data(isi + 199);
    const auto *isi_201 = buffer.data(isi + 201);
    const auto *isi_202 = buffer.data(isi + 202);
    const auto *isi_205 = buffer.data(isi + 205);
    const auto *isi_206 = buffer.data(isi + 206);
    const auto *isi_210 = buffer.data(isi + 210);
    const auto *isi_216 = buffer.data(isi + 216);
    const auto *isi_217 = buffer.data(isi + 217);
    const auto *isi_218 = buffer.data(isi + 218);
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
    const auto *isi_234 = buffer.data(isi + 234);
    const auto *isi_238 = buffer.data(isi + 238);
    const auto *isi_245 = buffer.data(isi + 245);
    const auto *isi_246 = buffer.data(isi + 246);
    const auto *isi_247 = buffer.data(isi + 247);
    const auto *isi_248 = buffer.data(isi + 248);
    const auto *isi_249 = buffer.data(isi + 249);
    const auto *isi_250 = buffer.data(isi + 250);
    const auto *isi_251 = buffer.data(isi + 251);
    const auto *isi_252 = buffer.data(isi + 252);
    const auto *isi_253 = buffer.data(isi + 253);
    const auto *isi_254 = buffer.data(isi + 254);
    const auto *isi_255 = buffer.data(isi + 255);
    const auto *isi_256 = buffer.data(isi + 256);
    const auto *isi_257 = buffer.data(isi + 257);
    const auto *isi_258 = buffer.data(isi + 258);
    const auto *isi_259 = buffer.data(isi + 259);
    const auto *isi_260 = buffer.data(isi + 260);
    const auto *isi_261 = buffer.data(isi + 261);
    const auto *isi_262 = buffer.data(isi + 262);
    const auto *isi_263 = buffer.data(isi + 263);
    const auto *isi_264 = buffer.data(isi + 264);
    const auto *isi_265 = buffer.data(isi + 265);
    const auto *isi_266 = buffer.data(isi + 266);
    const auto *isi_272 = buffer.data(isi + 272);
    const auto *isi_273 = buffer.data(isi + 273);
    const auto *isi_274 = buffer.data(isi + 274);
    const auto *isi_275 = buffer.data(isi + 275);
    const auto *isi_276 = buffer.data(isi + 276);
    const auto *isi_277 = buffer.data(isi + 277);
    const auto *isi_278 = buffer.data(isi + 278);
    const auto *isi_279 = buffer.data(isi + 279);
    const auto *isi_280 = buffer.data(isi + 280);
    const auto *isi_281 = buffer.data(isi + 281);
    const auto *isi_282 = buffer.data(isi + 282);
    const auto *isi_283 = buffer.data(isi + 283);
    const auto *isi_286 = buffer.data(isi + 286);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, hsk0_108, hsi_111, \
                         hsi_112, hsk1_108, ish0_146, ish1_146, isi_195, \
                         isi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * hsi_111[k]
                   + f_3 * pc_y[k] * isi_195[k];

        t_251[k] = f_1 * ish0_146[k]
                   - f_2 * ish1_146[k]
                   + f_3 * pc_z[k] * isi_195[k];

        t_252[k] = pa_z[k] * hsk0_108[k]
                   - f_12 * pc_z[k] * hsk1_108[k];

        t_253[k] = f_14 * hsi_112[k]
                   + f_3 * pc_y[k] * isi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_y, pc_z, hsk0_111, hsi_84, hsi_114, \
                         hsk1_111, isi_196, isi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * hsi_84[k]
                   + f_3 * pc_z[k] * isi_196[k];

        t_255[k] = pa_z[k] * hsk0_111[k]
                   - f_12 * pc_z[k] * hsk1_111[k];

        t_256[k] = f_14 * hsi_114[k]
                   + f_3 * pc_y[k] * isi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, hsk0_114, hsi_87, hsi_201, \
                         hsk1_114, ish0_152, ish1_152, isi_199, \
                         isi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * hsi_201[k]
                   + f_10 * ish0_152[k]
                   - f_11 * ish1_152[k]
                   + f_3 * pc_x[k] * isi_201[k];

        t_258[k] = pa_z[k] * hsk0_114[k]
                   - f_12 * pc_z[k] * hsk1_114[k];

        t_259[k] = f_13 * hsi_87[k]
                   + f_3 * pc_z[k] * isi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_y, pc_z, hsk0_118, hsi_117, \
                         hsi_205, hsk1_118, ish0_156, ish1_156, isi_201, \
                         isi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * hsi_117[k]
                   + f_3 * pc_y[k] * isi_201[k];

        t_261[k] = f_15 * hsi_205[k]
                   + f_8 * ish0_156[k]
                   - f_9 * ish1_156[k]
                   + f_3 * pc_x[k] * isi_205[k];

        t_262[k] = pa_z[k] * hsk0_118[k]
                   - f_12 * pc_z[k] * hsk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pa_z, pc_y, pc_z, hsk0_120, hsi_90, hsi_91, \
                         hsi_121, hsk1_120, isi_202, isi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * hsi_90[k]
                   + f_3 * pc_z[k] * isi_202[k];

        t_264[k] = pa_z[k] * hsk0_120[k]
                   + f_14 * hsi_91[k]
                   - f_12 * pc_z[k] * hsk1_120[k];

        t_265[k] = f_14 * hsi_121[k]
                   + f_3 * pc_y[k] * isi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, hsk0_123, hsi_94, hsi_210, \
                         hsk1_123, ish0_161, ish1_161, isi_206, \
                         isi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_15 * hsi_210[k]
                   + f_6 * ish0_161[k]
                   - f_7 * ish1_161[k]
                   + f_3 * pc_x[k] * isi_210[k];

        t_267[k] = pa_z[k] * hsk0_123[k]
                   - f_12 * pc_z[k] * hsk1_123[k];

        t_268[k] = f_13 * hsi_94[k]
                   + f_3 * pc_z[k] * isi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_z, pc_y, pc_z, hsk0_125, hsk0_126, hsi_95, \
                         hsi_96, hsi_126, hsk1_125, hsk1_126, isi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_z[k] * hsk0_125[k]
                   + f_14 * hsi_95[k]
                   - f_12 * pc_z[k] * hsk1_125[k];

        t_270[k] = pa_z[k] * hsk0_126[k]
                   + f_15 * hsi_96[k]
                   - f_12 * pc_z[k] * hsk1_126[k];

        t_271[k] = f_14 * hsi_126[k]
                   + f_3 * pc_y[k] * isi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, hsi_216, hsi_217, hsi_218, hsi_219, \
                         ish0_167, ish1_167, isi_216, isi_217, isi_218, \
                         isi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_15 * hsi_216[k]
                   + f_4 * ish0_167[k]
                   - f_5 * ish1_167[k]
                   + f_3 * pc_x[k] * isi_216[k];

        t_273[k] = f_15 * hsi_217[k]
                   + f_3 * pc_x[k] * isi_217[k];

        t_274[k] = f_15 * hsi_218[k]
                   + f_3 * pc_x[k] * isi_218[k];

        t_275[k] = f_15 * hsi_219[k]
                   + f_3 * pc_x[k] * isi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, hsi_220, hsi_221, hsi_222, hsi_223, \
                         isi_220, isi_221, isi_222, isi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_15 * hsi_220[k]
                   + f_3 * pc_x[k] * isi_220[k];

        t_277[k] = f_15 * hsi_221[k]
                   + f_3 * pc_x[k] * isi_221[k];

        t_278[k] = f_15 * hsi_222[k]
                   + f_3 * pc_x[k] * isi_222[k];

        t_279[k] = f_15 * hsi_223[k]
                   + f_3 * pc_x[k] * isi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_z, pc_y, pc_z, hsk0_136, hsi_105, hsi_135, \
                         hsk1_136, ish0_164, ish1_164, isi_217, \
                         isi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * hsk0_136[k]
                   - f_12 * pc_z[k] * hsk1_136[k];

        t_281[k] = f_13 * hsi_105[k]
                   + f_3 * pc_z[k] * isi_217[k];

        t_282[k] = f_14 * hsi_135[k]
                   + f_10 * ish0_164[k]
                   - f_11 * ish1_164[k]
                   + f_3 * pc_y[k] * isi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, hsi_136, hsi_137, hsi_138, ish0_165, \
                         ish0_166, ish0_167, ish1_165, ish1_166, ish1_167, isi_220, isi_221, \
                         isi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * hsi_136[k]
                   + f_8 * ish0_165[k]
                   - f_9 * ish1_165[k]
                   + f_3 * pc_y[k] * isi_220[k];

        t_284[k] = f_14 * hsi_137[k]
                   + f_6 * ish0_166[k]
                   - f_7 * ish1_166[k]
                   + f_3 * pc_y[k] * isi_221[k];

        t_285[k] = f_14 * hsi_138[k]
                   + f_4 * ish0_167[k]
                   - f_5 * ish1_167[k]
                   + f_3 * pc_y[k] * isi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, hsk0_180, hsi_111, \
                         hsi_139, hsi_140, hsk1_180, ish0_167, ish1_167, isi_223, \
                         isi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * hsi_139[k]
                   + f_3 * pc_y[k] * isi_223[k];

        t_287[k] = f_13 * hsi_111[k]
                   + f_1 * ish0_167[k]
                   - f_2 * ish1_167[k]
                   + f_3 * pc_z[k] * isi_223[k];

        t_288[k] = pa_y[k] * hsk0_180[k]
                   - f_12 * pc_y[k] * hsk1_180[k];

        t_289[k] = f_13 * hsi_140[k]
                   + f_3 * pc_y[k] * isi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, pc_z, hsk0_183, hsk0_185, \
                         hsi_112, hsi_141, hsi_142, hsk1_183, hsk1_185, isi_224, \
                         isi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * hsi_112[k]
                   + f_3 * pc_z[k] * isi_224[k];

        t_291[k] = pa_y[k] * hsk0_183[k]
                   + f_14 * hsi_141[k]
                   - f_12 * pc_y[k] * hsk1_183[k];

        t_292[k] = f_13 * hsi_142[k]
                   + f_3 * pc_y[k] * isi_226[k];

        t_293[k] = pa_y[k] * hsk0_185[k]
                   - f_12 * pc_y[k] * hsk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pc_y, pc_z, hsk0_186, hsk0_189, \
                         hsi_115, hsi_143, hsi_145, hsk1_186, hsk1_189, isi_227, \
                         isi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * hsk0_186[k]
                   + f_15 * hsi_143[k]
                   - f_12 * pc_y[k] * hsk1_186[k];

        t_295[k] = f_14 * hsi_115[k]
                   + f_3 * pc_z[k] * isi_227[k];

        t_296[k] = f_13 * hsi_145[k]
                   + f_3 * pc_y[k] * isi_229[k];

        t_297[k] = pa_y[k] * hsk0_189[k]
                   - f_12 * pc_y[k] * hsk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_y, pc_y, pc_z, hsk0_190, hsk0_192, hsi_118, \
                         hsi_146, hsi_148, hsk1_190, hsk1_192, \
                         isi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_y[k] * hsk0_190[k]
                   + f_16 * hsi_146[k]
                   - f_12 * pc_y[k] * hsk1_190[k];

        t_299[k] = f_14 * hsi_118[k]
                   + f_3 * pc_z[k] * isi_230[k];

        t_300[k] = pa_y[k] * hsk0_192[k]
                   + f_14 * hsi_148[k]
                   - f_12 * pc_y[k] * hsk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, pc_y, pc_z, hsk0_194, hsk0_195, \
                         hsi_122, hsi_149, hsi_150, hsk1_194, hsk1_195, isi_233, \
                         isi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * hsi_149[k]
                   + f_3 * pc_y[k] * isi_233[k];

        t_302[k] = pa_y[k] * hsk0_194[k]
                   - f_12 * pc_y[k] * hsk1_194[k];

        t_303[k] = pa_y[k] * hsk0_195[k]
                   + f_17 * hsi_150[k]
                   - f_12 * pc_y[k] * hsk1_195[k];

        t_304[k] = f_14 * hsi_122[k]
                   + f_3 * pc_z[k] * isi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pc_y, hsk0_197, hsk0_198, hsk0_200, \
                         hsi_152, hsi_153, hsi_154, hsk1_197, hsk1_198, hsk1_200, \
                         isi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_y[k] * hsk0_197[k]
                   + f_15 * hsi_152[k]
                   - f_12 * pc_y[k] * hsk1_197[k];

        t_306[k] = pa_y[k] * hsk0_198[k]
                   + f_14 * hsi_153[k]
                   - f_12 * pc_y[k] * hsk1_198[k];

        t_307[k] = f_13 * hsi_154[k]
                   + f_3 * pc_y[k] * isi_238[k];

        t_308[k] = pa_y[k] * hsk0_200[k]
                   - f_12 * pc_y[k] * hsk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, hsi_245, hsi_246, hsi_247, \
                         hsi_248, hsi_249, isi_245, isi_246, isi_247, isi_248, \
                         isi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_15 * hsi_245[k]
                   + f_3 * pc_x[k] * isi_245[k];

        t_310[k] = f_15 * hsi_246[k]
                   + f_3 * pc_x[k] * isi_246[k];

        t_311[k] = f_15 * hsi_247[k]
                   + f_3 * pc_x[k] * isi_247[k];

        t_312[k] = f_15 * hsi_248[k]
                   + f_3 * pc_x[k] * isi_248[k];

        t_313[k] = f_15 * hsi_249[k]
                   + f_3 * pc_x[k] * isi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, hsi_133, hsi_161, \
                         hsi_250, hsi_251, ish0_183, ish1_183, isi_245, isi_250, \
                         isi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_15 * hsi_250[k]
                   + f_3 * pc_x[k] * isi_250[k];

        t_315[k] = f_15 * hsi_251[k]
                   + f_3 * pc_x[k] * isi_251[k];

        t_316[k] = f_13 * hsi_161[k]
                   + f_1 * ish0_183[k]
                   - f_2 * ish1_183[k]
                   + f_3 * pc_y[k] * isi_245[k];

        t_317[k] = f_14 * hsi_133[k]
                   + f_3 * pc_z[k] * isi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, hsi_163, hsi_164, hsi_165, ish0_185, \
                         ish0_186, ish0_187, ish1_185, ish1_186, ish1_187, isi_247, isi_248, \
                         isi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * hsi_163[k]
                   + f_10 * ish0_185[k]
                   - f_11 * ish1_185[k]
                   + f_3 * pc_y[k] * isi_247[k];

        t_319[k] = f_13 * hsi_164[k]
                   + f_8 * ish0_186[k]
                   - f_9 * ish1_186[k]
                   + f_3 * pc_y[k] * isi_248[k];

        t_320[k] = f_13 * hsi_165[k]
                   + f_6 * ish0_187[k]
                   - f_7 * ish1_187[k]
                   + f_3 * pc_y[k] * isi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, hsk0_215, hsi_166, hsi_167, \
                         hsk1_215, ish0_188, ish1_188, isi_250, \
                         isi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * hsi_166[k]
                   + f_4 * ish0_188[k]
                   - f_5 * ish1_188[k]
                   + f_3 * pc_y[k] * isi_250[k];

        t_322[k] = f_13 * hsi_167[k]
                   + f_3 * pc_y[k] * isi_251[k];

        t_323[k] = pa_y[k] * hsk0_215[k]
                   - f_12 * pc_y[k] * hsk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, hsi_140, \
                         hsi_252, ish0_189, ish1_189, isi_252, isi_253, \
                         isi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_15 * hsi_252[k]
                   + f_1 * ish0_189[k]
                   - f_2 * ish1_189[k]
                   + f_3 * pc_x[k] * isi_252[k];

        t_325[k] = f_3 * pc_y[k] * isi_252[k];

        t_326[k] = f_15 * hsi_140[k]
                   + f_3 * pc_z[k] * isi_252[k];

        t_327[k] = f_4 * ish0_189[k]
                   - f_5 * ish1_189[k]
                   + f_3 * pc_y[k] * isi_253[k];

        t_328[k] = f_3 * pc_y[k] * isi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, hsi_257, ish0_190, ish0_191, \
                         ish0_194, ish1_190, ish1_191, ish1_194, isi_255, isi_256, \
                         isi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_15 * hsi_257[k]
                   + f_10 * ish0_194[k]
                   - f_11 * ish1_194[k]
                   + f_3 * pc_x[k] * isi_257[k];

        t_330[k] = f_6 * ish0_190[k]
                   - f_7 * ish1_190[k]
                   + f_3 * pc_y[k] * isi_255[k];

        t_331[k] = f_4 * ish0_191[k]
                   - f_5 * ish1_191[k]
                   + f_3 * pc_y[k] * isi_256[k];

        t_332[k] = f_3 * pc_y[k] * isi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, hsi_261, ish0_192, ish0_193, \
                         ish0_198, ish1_192, ish1_193, ish1_198, isi_258, isi_259, \
                         isi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_15 * hsi_261[k]
                   + f_8 * ish0_198[k]
                   - f_9 * ish1_198[k]
                   + f_3 * pc_x[k] * isi_261[k];

        t_334[k] = f_8 * ish0_192[k]
                   - f_9 * ish1_192[k]
                   + f_3 * pc_y[k] * isi_258[k];

        t_335[k] = f_6 * ish0_193[k]
                   - f_7 * ish1_193[k]
                   + f_3 * pc_y[k] * isi_259[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, hsi_266, ish0_194, ish0_203, \
                         ish1_194, ish1_203, isi_260, isi_261, \
                         isi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_4 * ish0_194[k]
                   - f_5 * ish1_194[k]
                   + f_3 * pc_y[k] * isi_260[k];

        t_337[k] = f_3 * pc_y[k] * isi_261[k];

        t_338[k] = f_15 * hsi_266[k]
                   + f_6 * ish0_203[k]
                   - f_7 * ish1_203[k]
                   + f_3 * pc_x[k] * isi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_y, ish0_195, ish0_196, ish0_197, ish1_195, \
                         ish1_196, ish1_197, isi_262, isi_263, \
                         isi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_10 * ish0_195[k]
                   - f_11 * ish1_195[k]
                   + f_3 * pc_y[k] * isi_262[k];

        t_340[k] = f_8 * ish0_196[k]
                   - f_9 * ish1_196[k]
                   + f_3 * pc_y[k] * isi_263[k];

        t_341[k] = f_6 * ish0_197[k]
                   - f_7 * ish1_197[k]
                   + f_3 * pc_y[k] * isi_264[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_y, hsi_272, hsi_273, ish0_198, \
                         ish0_209, ish1_198, ish1_209, isi_265, isi_266, isi_272, \
                         isi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_4 * ish0_198[k]
                   - f_5 * ish1_198[k]
                   + f_3 * pc_y[k] * isi_265[k];

        t_343[k] = f_3 * pc_y[k] * isi_266[k];

        t_344[k] = f_15 * hsi_272[k]
                   + f_4 * ish0_209[k]
                   - f_5 * ish1_209[k]
                   + f_3 * pc_x[k] * isi_272[k];

        t_345[k] = f_15 * hsi_273[k]
                   + f_3 * pc_x[k] * isi_273[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pc_x, pc_y, hsi_274, hsi_275, \
                         hsi_276, hsi_277, isi_272, isi_274, isi_275, isi_276, \
                         isi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_15 * hsi_274[k]
                   + f_3 * pc_x[k] * isi_274[k];

        t_347[k] = f_15 * hsi_275[k]
                   + f_3 * pc_x[k] * isi_275[k];

        t_348[k] = f_15 * hsi_276[k]
                   + f_3 * pc_x[k] * isi_276[k];

        t_349[k] = f_15 * hsi_277[k]
                   + f_3 * pc_x[k] * isi_277[k];

        t_350[k] = f_3 * pc_y[k] * isi_272[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_x, pc_y, hsi_279, ish0_204, ish0_205, \
                         ish1_204, ish1_205, isi_273, isi_274, \
                         isi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_15 * hsi_279[k]
                   + f_3 * pc_x[k] * isi_279[k];

        t_352[k] = f_1 * ish0_204[k]
                   - f_2 * ish1_204[k]
                   + f_3 * pc_y[k] * isi_273[k];

        t_353[k] = f_18 * ish0_205[k]
                   - f_19 * ish1_205[k]
                   + f_3 * pc_y[k] * isi_274[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_y, ish0_206, ish0_207, ish0_208, ish1_206, \
                         ish1_207, ish1_208, isi_275, isi_276, \
                         isi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_10 * ish0_206[k]
                   - f_11 * ish1_206[k]
                   + f_3 * pc_y[k] * isi_275[k];

        t_355[k] = f_8 * ish0_207[k]
                   - f_9 * ish1_207[k]
                   + f_3 * pc_y[k] * isi_276[k];

        t_356[k] = f_6 * ish0_208[k]
                   - f_7 * ish1_208[k]
                   + f_3 * pc_y[k] * isi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pc_x, pc_y, pc_z, hsi_167, hsi_280, \
                         ish0_209, ish0_210, ish1_209, ish1_210, isi_278, isi_279, \
                         isi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_4 * ish0_209[k]
                   - f_5 * ish1_209[k]
                   + f_3 * pc_y[k] * isi_278[k];

        t_358[k] = f_3 * pc_y[k] * isi_279[k];

        t_359[k] = f_15 * hsi_167[k]
                   + f_1 * ish0_209[k]
                   - f_2 * ish1_209[k]
                   + f_3 * pc_z[k] * isi_279[k];

        t_360[k] = f_14 * hsi_280[k]
                   + f_1 * ish0_210[k]
                   - f_2 * ish1_210[k]
                   + f_3 * pc_x[k] * isi_280[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, hsi_168, hsi_283, \
                         ish0_213, ish1_213, isi_280, isi_281, \
                         isi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_16 * hsi_168[k]
                   + f_3 * pc_y[k] * isi_280[k];

        t_362[k] = f_3 * pc_z[k] * isi_280[k];

        t_363[k] = f_14 * hsi_283[k]
                   + f_10 * ish0_213[k]
                   - f_11 * ish1_213[k]
                   + f_3 * pc_x[k] * isi_283[k];

        t_364[k] = f_3 * pc_z[k] * isi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pc_x, pc_z, hsi_286, ish0_210, ish0_216, \
                         ish1_210, ish1_216, isi_282, isi_283, \
                         isi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * ish0_210[k]
                   - f_5 * ish1_210[k]
                   + f_3 * pc_z[k] * isi_282[k];

        t_366[k] = f_14 * hsi_286[k]
                   + f_8 * ish0_216[k]
                   - f_9 * ish1_216[k]
                   + f_3 * pc_x[k] * isi_286[k];

        t_367[k] = f_3 * pc_z[k] * isi_283[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
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

    const auto *hsk0_216 = buffer.data(hsk0 + 216);
    const auto *hsk0_219 = buffer.data(hsk0 + 219);
    const auto *hsk0_222 = buffer.data(hsk0 + 222);
    const auto *hsk0_226 = buffer.data(hsk0 + 226);
    const auto *hsk0_228 = buffer.data(hsk0 + 228);
    const auto *hsk0_231 = buffer.data(hsk0 + 231);
    const auto *hsk0_233 = buffer.data(hsk0 + 233);
    const auto *hsk0_234 = buffer.data(hsk0 + 234);
    const auto *hsk0_244 = buffer.data(hsk0 + 244);
    const auto *hsk0_324 = buffer.data(hsk0 + 324);
    const auto *hsk0_327 = buffer.data(hsk0 + 327);
    const auto *hsk0_329 = buffer.data(hsk0 + 329);
    const auto *hsk0_330 = buffer.data(hsk0 + 330);
    const auto *hsk0_333 = buffer.data(hsk0 + 333);
    const auto *hsk0_334 = buffer.data(hsk0 + 334);
    const auto *hsk0_336 = buffer.data(hsk0 + 336);

    const auto *hsi_168 = buffer.data(hsi + 168);
    const auto *hsi_171 = buffer.data(hsi + 171);
    const auto *hsi_173 = buffer.data(hsi + 173);
    const auto *hsi_174 = buffer.data(hsi + 174);
    const auto *hsi_175 = buffer.data(hsi + 175);
    const auto *hsi_177 = buffer.data(hsi + 177);
    const auto *hsi_178 = buffer.data(hsi + 178);
    const auto *hsi_179 = buffer.data(hsi + 179);
    const auto *hsi_180 = buffer.data(hsi + 180);
    const auto *hsi_182 = buffer.data(hsi + 182);
    const auto *hsi_189 = buffer.data(hsi + 189);
    const auto *hsi_195 = buffer.data(hsi + 195);
    const auto *hsi_196 = buffer.data(hsi + 196);
    const auto *hsi_198 = buffer.data(hsi + 198);
    const auto *hsi_199 = buffer.data(hsi + 199);
    const auto *hsi_201 = buffer.data(hsi + 201);
    const auto *hsi_202 = buffer.data(hsi + 202);
    const auto *hsi_205 = buffer.data(hsi + 205);
    const auto *hsi_206 = buffer.data(hsi + 206);
    const auto *hsi_210 = buffer.data(hsi + 210);
    const auto *hsi_217 = buffer.data(hsi + 217);
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
    const auto *hsi_238 = buffer.data(hsi + 238);
    const auto *hsi_245 = buffer.data(hsi + 245);
    const auto *hsi_247 = buffer.data(hsi + 247);
    const auto *hsi_248 = buffer.data(hsi + 248);
    const auto *hsi_249 = buffer.data(hsi + 249);
    const auto *hsi_250 = buffer.data(hsi + 250);
    const auto *hsi_251 = buffer.data(hsi + 251);
    const auto *hsi_252 = buffer.data(hsi + 252);
    const auto *hsi_253 = buffer.data(hsi + 253);
    const auto *hsi_254 = buffer.data(hsi + 254);
    const auto *hsi_255 = buffer.data(hsi + 255);
    const auto *hsi_257 = buffer.data(hsi + 257);
    const auto *hsi_258 = buffer.data(hsi + 258);
    const auto *hsi_260 = buffer.data(hsi + 260);
    const auto *hsi_290 = buffer.data(hsi + 290);
    const auto *hsi_295 = buffer.data(hsi + 295);
    const auto *hsi_301 = buffer.data(hsi + 301);
    const auto *hsi_303 = buffer.data(hsi + 303);
    const auto *hsi_304 = buffer.data(hsi + 304);
    const auto *hsi_305 = buffer.data(hsi + 305);
    const auto *hsi_306 = buffer.data(hsi + 306);
    const auto *hsi_307 = buffer.data(hsi + 307);
    const auto *hsi_313 = buffer.data(hsi + 313);
    const auto *hsi_317 = buffer.data(hsi + 317);
    const auto *hsi_322 = buffer.data(hsi + 322);
    const auto *hsi_328 = buffer.data(hsi + 328);
    const auto *hsi_329 = buffer.data(hsi + 329);
    const auto *hsi_330 = buffer.data(hsi + 330);
    const auto *hsi_331 = buffer.data(hsi + 331);
    const auto *hsi_332 = buffer.data(hsi + 332);
    const auto *hsi_333 = buffer.data(hsi + 333);
    const auto *hsi_334 = buffer.data(hsi + 334);
    const auto *hsi_335 = buffer.data(hsi + 335);
    const auto *hsi_336 = buffer.data(hsi + 336);
    const auto *hsi_339 = buffer.data(hsi + 339);
    const auto *hsi_341 = buffer.data(hsi + 341);
    const auto *hsi_342 = buffer.data(hsi + 342);
    const auto *hsi_345 = buffer.data(hsi + 345);
    const auto *hsi_346 = buffer.data(hsi + 346);
    const auto *hsi_348 = buffer.data(hsi + 348);
    const auto *hsi_350 = buffer.data(hsi + 350);
    const auto *hsi_351 = buffer.data(hsi + 351);
    const auto *hsi_353 = buffer.data(hsi + 353);
    const auto *hsi_354 = buffer.data(hsi + 354);
    const auto *hsi_356 = buffer.data(hsi + 356);
    const auto *hsi_357 = buffer.data(hsi + 357);
    const auto *hsi_358 = buffer.data(hsi + 358);
    const auto *hsi_359 = buffer.data(hsi + 359);
    const auto *hsi_360 = buffer.data(hsi + 360);
    const auto *hsi_361 = buffer.data(hsi + 361);
    const auto *hsi_362 = buffer.data(hsi + 362);
    const auto *hsi_363 = buffer.data(hsi + 363);

    const auto *hsk1_216 = buffer.data(hsk1 + 216);
    const auto *hsk1_219 = buffer.data(hsk1 + 219);
    const auto *hsk1_222 = buffer.data(hsk1 + 222);
    const auto *hsk1_226 = buffer.data(hsk1 + 226);
    const auto *hsk1_228 = buffer.data(hsk1 + 228);
    const auto *hsk1_231 = buffer.data(hsk1 + 231);
    const auto *hsk1_233 = buffer.data(hsk1 + 233);
    const auto *hsk1_234 = buffer.data(hsk1 + 234);
    const auto *hsk1_244 = buffer.data(hsk1 + 244);
    const auto *hsk1_324 = buffer.data(hsk1 + 324);
    const auto *hsk1_327 = buffer.data(hsk1 + 327);
    const auto *hsk1_329 = buffer.data(hsk1 + 329);
    const auto *hsk1_330 = buffer.data(hsk1 + 330);
    const auto *hsk1_333 = buffer.data(hsk1 + 333);
    const auto *hsk1_334 = buffer.data(hsk1 + 334);
    const auto *hsk1_336 = buffer.data(hsk1 + 336);

    const auto *ish0_212 = buffer.data(ish0 + 212);
    const auto *ish0_213 = buffer.data(ish0 + 213);
    const auto *ish0_215 = buffer.data(ish0 + 215);
    const auto *ish0_216 = buffer.data(ish0 + 216);
    const auto *ish0_217 = buffer.data(ish0 + 217);
    const auto *ish0_219 = buffer.data(ish0 + 219);
    const auto *ish0_220 = buffer.data(ish0 + 220);
    const auto *ish0_225 = buffer.data(ish0 + 225);
    const auto *ish0_226 = buffer.data(ish0 + 226);
    const auto *ish0_227 = buffer.data(ish0 + 227);
    const auto *ish0_228 = buffer.data(ish0 + 228);
    const auto *ish0_230 = buffer.data(ish0 + 230);
    const auto *ish0_236 = buffer.data(ish0 + 236);
    const auto *ish0_240 = buffer.data(ish0 + 240);
    const auto *ish0_245 = buffer.data(ish0 + 245);
    const auto *ish0_248 = buffer.data(ish0 + 248);
    const auto *ish0_249 = buffer.data(ish0 + 249);
    const auto *ish0_250 = buffer.data(ish0 + 250);
    const auto *ish0_251 = buffer.data(ish0 + 251);
    const auto *ish0_252 = buffer.data(ish0 + 252);
    const auto *ish0_255 = buffer.data(ish0 + 255);
    const auto *ish0_257 = buffer.data(ish0 + 257);
    const auto *ish0_258 = buffer.data(ish0 + 258);
    const auto *ish0_261 = buffer.data(ish0 + 261);
    const auto *ish0_262 = buffer.data(ish0 + 262);
    const auto *ish0_264 = buffer.data(ish0 + 264);
    const auto *ish0_266 = buffer.data(ish0 + 266);
    const auto *ish0_267 = buffer.data(ish0 + 267);
    const auto *ish0_269 = buffer.data(ish0 + 269);
    const auto *ish0_270 = buffer.data(ish0 + 270);
    const auto *ish0_271 = buffer.data(ish0 + 271);
    const auto *ish0_272 = buffer.data(ish0 + 272);

    const auto *ish1_212 = buffer.data(ish1 + 212);
    const auto *ish1_213 = buffer.data(ish1 + 213);
    const auto *ish1_215 = buffer.data(ish1 + 215);
    const auto *ish1_216 = buffer.data(ish1 + 216);
    const auto *ish1_217 = buffer.data(ish1 + 217);
    const auto *ish1_219 = buffer.data(ish1 + 219);
    const auto *ish1_220 = buffer.data(ish1 + 220);
    const auto *ish1_225 = buffer.data(ish1 + 225);
    const auto *ish1_226 = buffer.data(ish1 + 226);
    const auto *ish1_227 = buffer.data(ish1 + 227);
    const auto *ish1_228 = buffer.data(ish1 + 228);
    const auto *ish1_230 = buffer.data(ish1 + 230);
    const auto *ish1_236 = buffer.data(ish1 + 236);
    const auto *ish1_240 = buffer.data(ish1 + 240);
    const auto *ish1_245 = buffer.data(ish1 + 245);
    const auto *ish1_248 = buffer.data(ish1 + 248);
    const auto *ish1_249 = buffer.data(ish1 + 249);
    const auto *ish1_250 = buffer.data(ish1 + 250);
    const auto *ish1_251 = buffer.data(ish1 + 251);
    const auto *ish1_252 = buffer.data(ish1 + 252);
    const auto *ish1_255 = buffer.data(ish1 + 255);
    const auto *ish1_257 = buffer.data(ish1 + 257);
    const auto *ish1_258 = buffer.data(ish1 + 258);
    const auto *ish1_261 = buffer.data(ish1 + 261);
    const auto *ish1_262 = buffer.data(ish1 + 262);
    const auto *ish1_264 = buffer.data(ish1 + 264);
    const auto *ish1_266 = buffer.data(ish1 + 266);
    const auto *ish1_267 = buffer.data(ish1 + 267);
    const auto *ish1_269 = buffer.data(ish1 + 269);
    const auto *ish1_270 = buffer.data(ish1 + 270);
    const auto *ish1_271 = buffer.data(ish1 + 271);
    const auto *ish1_272 = buffer.data(ish1 + 272);

    const auto *isi_285 = buffer.data(isi + 285);
    const auto *isi_286 = buffer.data(isi + 286);
    const auto *isi_287 = buffer.data(isi + 287);
    const auto *isi_289 = buffer.data(isi + 289);
    const auto *isi_290 = buffer.data(isi + 290);
    const auto *isi_291 = buffer.data(isi + 291);
    const auto *isi_292 = buffer.data(isi + 292);
    const auto *isi_294 = buffer.data(isi + 294);
    const auto *isi_295 = buffer.data(isi + 295);
    const auto *isi_301 = buffer.data(isi + 301);
    const auto *isi_302 = buffer.data(isi + 302);
    const auto *isi_303 = buffer.data(isi + 303);
    const auto *isi_304 = buffer.data(isi + 304);
    const auto *isi_305 = buffer.data(isi + 305);
    const auto *isi_306 = buffer.data(isi + 306);
    const auto *isi_307 = buffer.data(isi + 307);
    const auto *isi_308 = buffer.data(isi + 308);
    const auto *isi_310 = buffer.data(isi + 310);
    const auto *isi_311 = buffer.data(isi + 311);
    const auto *isi_313 = buffer.data(isi + 313);
    const auto *isi_314 = buffer.data(isi + 314);
    const auto *isi_317 = buffer.data(isi + 317);
    const auto *isi_318 = buffer.data(isi + 318);
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
    const auto *isi_338 = buffer.data(isi + 338);
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
    const auto *isi_364 = buffer.data(isi + 364);
    const auto *isi_366 = buffer.data(isi + 366);
    const auto *isi_367 = buffer.data(isi + 367);
    const auto *isi_369 = buffer.data(isi + 369);
    const auto *isi_370 = buffer.data(isi + 370);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pc_x, pc_y, pc_z, hsi_173, hsi_290, \
                         ish0_212, ish0_220, ish1_212, ish1_220, isi_285, isi_286, \
                         isi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_16 * hsi_173[k]
                   + f_3 * pc_y[k] * isi_285[k];

        t_369[k] = f_6 * ish0_212[k]
                   - f_7 * ish1_212[k]
                   + f_3 * pc_z[k] * isi_285[k];

        t_370[k] = f_14 * hsi_290[k]
                   + f_6 * ish0_220[k]
                   - f_7 * ish1_220[k]
                   + f_3 * pc_x[k] * isi_290[k];

        t_371[k] = f_3 * pc_z[k] * isi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, hsi_177, ish0_213, ish0_215, \
                         ish1_213, ish1_215, isi_287, isi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * ish0_213[k]
                   - f_5 * ish1_213[k]
                   + f_3 * pc_z[k] * isi_287[k];

        t_373[k] = f_16 * hsi_177[k]
                   + f_3 * pc_y[k] * isi_289[k];

        t_374[k] = f_8 * ish0_215[k]
                   - f_9 * ish1_215[k]
                   + f_3 * pc_z[k] * isi_289[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_z, hsi_295, ish0_216, ish0_225, \
                         ish1_216, ish1_225, isi_290, isi_291, \
                         isi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_14 * hsi_295[k]
                   + f_4 * ish0_225[k]
                   - f_5 * ish1_225[k]
                   + f_3 * pc_x[k] * isi_295[k];

        t_376[k] = f_3 * pc_z[k] * isi_290[k];

        t_377[k] = f_4 * ish0_216[k]
                   - f_5 * ish1_216[k]
                   + f_3 * pc_z[k] * isi_291[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, hsi_182, hsi_301, \
                         ish0_217, ish0_219, ish1_217, ish1_219, isi_292, isi_294, \
                         isi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * ish0_217[k]
                   - f_7 * ish1_217[k]
                   + f_3 * pc_z[k] * isi_292[k];

        t_379[k] = f_16 * hsi_182[k]
                   + f_3 * pc_y[k] * isi_294[k];

        t_380[k] = f_10 * ish0_219[k]
                   - f_11 * ish1_219[k]
                   + f_3 * pc_z[k] * isi_294[k];

        t_381[k] = f_14 * hsi_301[k]
                   + f_3 * pc_x[k] * isi_301[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_z, hsi_303, hsi_304, \
                         hsi_305, hsi_306, isi_295, isi_303, isi_304, isi_305, \
                         isi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_z[k] * isi_295[k];

        t_383[k] = f_14 * hsi_303[k]
                   + f_3 * pc_x[k] * isi_303[k];

        t_384[k] = f_14 * hsi_304[k]
                   + f_3 * pc_x[k] * isi_304[k];

        t_385[k] = f_14 * hsi_305[k]
                   + f_3 * pc_x[k] * isi_305[k];

        t_386[k] = f_14 * hsi_306[k]
                   + f_3 * pc_x[k] * isi_306[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pc_y, pc_z, hsi_189, hsi_307, \
                         ish0_225, ish1_225, isi_301, isi_302, \
                         isi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_14 * hsi_307[k]
                   + f_3 * pc_x[k] * isi_307[k];

        t_388[k] = f_16 * hsi_189[k]
                   + f_1 * ish0_225[k]
                   - f_2 * ish1_225[k]
                   + f_3 * pc_y[k] * isi_301[k];

        t_389[k] = f_3 * pc_z[k] * isi_301[k];

        t_390[k] = f_4 * ish0_225[k]
                   - f_5 * ish1_225[k]
                   + f_3 * pc_z[k] * isi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_z, ish0_226, ish0_227, ish0_228, ish1_226, \
                         ish1_227, ish1_228, isi_303, isi_304, \
                         isi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_6 * ish0_226[k]
                   - f_7 * ish1_226[k]
                   + f_3 * pc_z[k] * isi_303[k];

        t_392[k] = f_8 * ish0_227[k]
                   - f_9 * ish1_227[k]
                   + f_3 * pc_z[k] * isi_304[k];

        t_393[k] = f_10 * ish0_228[k]
                   - f_11 * ish1_228[k]
                   + f_3 * pc_z[k] * isi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_z, pc_y, pc_z, hsk0_216, hsi_195, \
                         hsi_196, hsk1_216, ish0_230, ish1_230, isi_307, \
                         isi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * hsi_195[k]
                   + f_3 * pc_y[k] * isi_307[k];

        t_395[k] = f_1 * ish0_230[k]
                   - f_2 * ish1_230[k]
                   + f_3 * pc_z[k] * isi_307[k];

        t_396[k] = pa_z[k] * hsk0_216[k]
                   - f_12 * pc_z[k] * hsk1_216[k];

        t_397[k] = f_15 * hsi_196[k]
                   + f_3 * pc_y[k] * isi_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_z, pc_y, pc_z, hsk0_219, hsi_168, hsi_198, \
                         hsk1_219, isi_308, isi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * hsi_168[k]
                   + f_3 * pc_z[k] * isi_308[k];

        t_399[k] = pa_z[k] * hsk0_219[k]
                   - f_12 * pc_z[k] * hsk1_219[k];

        t_400[k] = f_15 * hsi_198[k]
                   + f_3 * pc_y[k] * isi_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_z, pc_x, pc_z, hsk0_222, hsi_171, hsi_313, \
                         hsk1_222, ish0_236, ish1_236, isi_311, \
                         isi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_14 * hsi_313[k]
                   + f_10 * ish0_236[k]
                   - f_11 * ish1_236[k]
                   + f_3 * pc_x[k] * isi_313[k];

        t_402[k] = pa_z[k] * hsk0_222[k]
                   - f_12 * pc_z[k] * hsk1_222[k];

        t_403[k] = f_13 * hsi_171[k]
                   + f_3 * pc_z[k] * isi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_z, pc_x, pc_y, pc_z, hsk0_226, hsi_201, \
                         hsi_317, hsk1_226, ish0_240, ish1_240, isi_313, \
                         isi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * hsi_201[k]
                   + f_3 * pc_y[k] * isi_313[k];

        t_405[k] = f_14 * hsi_317[k]
                   + f_8 * ish0_240[k]
                   - f_9 * ish1_240[k]
                   + f_3 * pc_x[k] * isi_317[k];

        t_406[k] = pa_z[k] * hsk0_226[k]
                   - f_12 * pc_z[k] * hsk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_z, pc_y, pc_z, hsk0_228, hsi_174, hsi_175, \
                         hsi_205, hsk1_228, isi_314, isi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * hsi_174[k]
                   + f_3 * pc_z[k] * isi_314[k];

        t_408[k] = pa_z[k] * hsk0_228[k]
                   + f_14 * hsi_175[k]
                   - f_12 * pc_z[k] * hsk1_228[k];

        t_409[k] = f_15 * hsi_205[k]
                   + f_3 * pc_y[k] * isi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pa_z, pc_x, pc_z, hsk0_231, hsi_178, hsi_322, \
                         hsk1_231, ish0_245, ish1_245, isi_318, \
                         isi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_14 * hsi_322[k]
                   + f_6 * ish0_245[k]
                   - f_7 * ish1_245[k]
                   + f_3 * pc_x[k] * isi_322[k];

        t_411[k] = pa_z[k] * hsk0_231[k]
                   - f_12 * pc_z[k] * hsk1_231[k];

        t_412[k] = f_13 * hsi_178[k]
                   + f_3 * pc_z[k] * isi_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pa_z, pc_y, pc_z, hsk0_233, hsk0_234, hsi_179, \
                         hsi_180, hsi_210, hsk1_233, hsk1_234, \
                         isi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_z[k] * hsk0_233[k]
                   + f_14 * hsi_179[k]
                   - f_12 * pc_z[k] * hsk1_233[k];

        t_414[k] = pa_z[k] * hsk0_234[k]
                   + f_15 * hsi_180[k]
                   - f_12 * pc_z[k] * hsk1_234[k];

        t_415[k] = f_15 * hsi_210[k]
                   + f_3 * pc_y[k] * isi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, hsi_328, hsi_329, hsi_330, hsi_331, \
                         ish0_251, ish1_251, isi_328, isi_329, isi_330, \
                         isi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_14 * hsi_328[k]
                   + f_4 * ish0_251[k]
                   - f_5 * ish1_251[k]
                   + f_3 * pc_x[k] * isi_328[k];

        t_417[k] = f_14 * hsi_329[k]
                   + f_3 * pc_x[k] * isi_329[k];

        t_418[k] = f_14 * hsi_330[k]
                   + f_3 * pc_x[k] * isi_330[k];

        t_419[k] = f_14 * hsi_331[k]
                   + f_3 * pc_x[k] * isi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, hsi_332, hsi_333, hsi_334, hsi_335, \
                         isi_332, isi_333, isi_334, isi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_14 * hsi_332[k]
                   + f_3 * pc_x[k] * isi_332[k];

        t_421[k] = f_14 * hsi_333[k]
                   + f_3 * pc_x[k] * isi_333[k];

        t_422[k] = f_14 * hsi_334[k]
                   + f_3 * pc_x[k] * isi_334[k];

        t_423[k] = f_14 * hsi_335[k]
                   + f_3 * pc_x[k] * isi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_z, pc_y, pc_z, hsk0_244, hsi_189, hsi_219, \
                         hsk1_244, ish0_248, ish1_248, isi_329, \
                         isi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * hsk0_244[k]
                   - f_12 * pc_z[k] * hsk1_244[k];

        t_425[k] = f_13 * hsi_189[k]
                   + f_3 * pc_z[k] * isi_329[k];

        t_426[k] = f_15 * hsi_219[k]
                   + f_10 * ish0_248[k]
                   - f_11 * ish1_248[k]
                   + f_3 * pc_y[k] * isi_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, hsi_220, hsi_221, hsi_222, ish0_249, \
                         ish0_250, ish0_251, ish1_249, ish1_250, ish1_251, isi_332, isi_333, \
                         isi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * hsi_220[k]
                   + f_8 * ish0_249[k]
                   - f_9 * ish1_249[k]
                   + f_3 * pc_y[k] * isi_332[k];

        t_428[k] = f_15 * hsi_221[k]
                   + f_6 * ish0_250[k]
                   - f_7 * ish1_250[k]
                   + f_3 * pc_y[k] * isi_333[k];

        t_429[k] = f_15 * hsi_222[k]
                   + f_4 * ish0_251[k]
                   - f_5 * ish1_251[k]
                   + f_3 * pc_y[k] * isi_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, hsi_195, hsi_223, hsi_336, \
                         ish0_251, ish0_252, ish1_251, ish1_252, isi_335, \
                         isi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * hsi_223[k]
                   + f_3 * pc_y[k] * isi_335[k];

        t_431[k] = f_13 * hsi_195[k]
                   + f_1 * ish0_251[k]
                   - f_2 * ish1_251[k]
                   + f_3 * pc_z[k] * isi_335[k];

        t_432[k] = f_14 * hsi_336[k]
                   + f_1 * ish0_252[k]
                   - f_2 * ish1_252[k]
                   + f_3 * pc_x[k] * isi_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, hsi_196, hsi_224, \
                         hsi_226, hsi_339, ish0_255, ish1_255, isi_336, isi_338, \
                         isi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * hsi_224[k]
                   + f_3 * pc_y[k] * isi_336[k];

        t_434[k] = f_14 * hsi_196[k]
                   + f_3 * pc_z[k] * isi_336[k];

        t_435[k] = f_14 * hsi_339[k]
                   + f_10 * ish0_255[k]
                   - f_11 * ish1_255[k]
                   + f_3 * pc_x[k] * isi_339[k];

        t_436[k] = f_14 * hsi_226[k]
                   + f_3 * pc_y[k] * isi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, hsi_199, hsi_341, hsi_342, ish0_257, \
                         ish0_258, ish1_257, ish1_258, isi_339, isi_341, \
                         isi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_14 * hsi_341[k]
                   + f_10 * ish0_257[k]
                   - f_11 * ish1_257[k]
                   + f_3 * pc_x[k] * isi_341[k];

        t_438[k] = f_14 * hsi_342[k]
                   + f_8 * ish0_258[k]
                   - f_9 * ish1_258[k]
                   + f_3 * pc_x[k] * isi_342[k];

        t_439[k] = f_14 * hsi_199[k]
                   + f_3 * pc_z[k] * isi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, hsi_229, hsi_345, hsi_346, ish0_261, \
                         ish0_262, ish1_261, ish1_262, isi_341, isi_345, \
                         isi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * hsi_229[k]
                   + f_3 * pc_y[k] * isi_341[k];

        t_441[k] = f_14 * hsi_345[k]
                   + f_8 * ish0_261[k]
                   - f_9 * ish1_261[k]
                   + f_3 * pc_x[k] * isi_345[k];

        t_442[k] = f_14 * hsi_346[k]
                   + f_6 * ish0_262[k]
                   - f_7 * ish1_262[k]
                   + f_3 * pc_x[k] * isi_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, hsi_202, hsi_233, hsi_348, \
                         ish0_264, ish1_264, isi_342, isi_345, \
                         isi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * hsi_202[k]
                   + f_3 * pc_z[k] * isi_342[k];

        t_444[k] = f_14 * hsi_348[k]
                   + f_6 * ish0_264[k]
                   - f_7 * ish1_264[k]
                   + f_3 * pc_x[k] * isi_348[k];

        t_445[k] = f_14 * hsi_233[k]
                   + f_3 * pc_y[k] * isi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, hsi_206, hsi_350, hsi_351, ish0_266, \
                         ish0_267, ish1_266, ish1_267, isi_346, isi_350, \
                         isi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_14 * hsi_350[k]
                   + f_6 * ish0_266[k]
                   - f_7 * ish1_266[k]
                   + f_3 * pc_x[k] * isi_350[k];

        t_447[k] = f_14 * hsi_351[k]
                   + f_4 * ish0_267[k]
                   - f_5 * ish1_267[k]
                   + f_3 * pc_x[k] * isi_351[k];

        t_448[k] = f_14 * hsi_206[k]
                   + f_3 * pc_z[k] * isi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, hsi_238, hsi_353, hsi_354, ish0_269, \
                         ish0_270, ish1_269, ish1_270, isi_350, isi_353, \
                         isi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * hsi_353[k]
                   + f_4 * ish0_269[k]
                   - f_5 * ish1_269[k]
                   + f_3 * pc_x[k] * isi_353[k];

        t_450[k] = f_14 * hsi_354[k]
                   + f_4 * ish0_270[k]
                   - f_5 * ish1_270[k]
                   + f_3 * pc_x[k] * isi_354[k];

        t_451[k] = f_14 * hsi_238[k]
                   + f_3 * pc_y[k] * isi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, hsi_356, hsi_357, hsi_358, hsi_359, \
                         ish0_272, ish1_272, isi_356, isi_357, isi_358, \
                         isi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * hsi_356[k]
                   + f_4 * ish0_272[k]
                   - f_5 * ish1_272[k]
                   + f_3 * pc_x[k] * isi_356[k];

        t_453[k] = f_14 * hsi_357[k]
                   + f_3 * pc_x[k] * isi_357[k];

        t_454[k] = f_14 * hsi_358[k]
                   + f_3 * pc_x[k] * isi_358[k];

        t_455[k] = f_14 * hsi_359[k]
                   + f_3 * pc_x[k] * isi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, hsi_360, hsi_361, hsi_362, hsi_363, \
                         isi_360, isi_361, isi_362, isi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * hsi_360[k]
                   + f_3 * pc_x[k] * isi_360[k];

        t_457[k] = f_14 * hsi_361[k]
                   + f_3 * pc_x[k] * isi_361[k];

        t_458[k] = f_14 * hsi_362[k]
                   + f_3 * pc_x[k] * isi_362[k];

        t_459[k] = f_14 * hsi_363[k]
                   + f_3 * pc_x[k] * isi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, hsi_217, hsi_245, hsi_247, ish0_267, \
                         ish0_269, ish1_267, ish1_269, isi_357, \
                         isi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * hsi_245[k]
                   + f_1 * ish0_267[k]
                   - f_2 * ish1_267[k]
                   + f_3 * pc_y[k] * isi_357[k];

        t_461[k] = f_14 * hsi_217[k]
                   + f_3 * pc_z[k] * isi_357[k];

        t_462[k] = f_14 * hsi_247[k]
                   + f_10 * ish0_269[k]
                   - f_11 * ish1_269[k]
                   + f_3 * pc_y[k] * isi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, hsi_248, hsi_249, hsi_250, ish0_270, \
                         ish0_271, ish0_272, ish1_270, ish1_271, ish1_272, isi_360, isi_361, \
                         isi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * hsi_248[k]
                   + f_8 * ish0_270[k]
                   - f_9 * ish1_270[k]
                   + f_3 * pc_y[k] * isi_360[k];

        t_464[k] = f_14 * hsi_249[k]
                   + f_6 * ish0_271[k]
                   - f_7 * ish1_271[k]
                   + f_3 * pc_y[k] * isi_361[k];

        t_465[k] = f_14 * hsi_250[k]
                   + f_4 * ish0_272[k]
                   - f_5 * ish1_272[k]
                   + f_3 * pc_y[k] * isi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_y, pc_y, pc_z, hsk0_324, hsi_223, \
                         hsi_251, hsi_252, hsk1_324, ish0_272, ish1_272, isi_363, \
                         isi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * hsi_251[k]
                   + f_3 * pc_y[k] * isi_363[k];

        t_467[k] = f_14 * hsi_223[k]
                   + f_1 * ish0_272[k]
                   - f_2 * ish1_272[k]
                   + f_3 * pc_z[k] * isi_363[k];

        t_468[k] = pa_y[k] * hsk0_324[k]
                   - f_12 * pc_y[k] * hsk1_324[k];

        t_469[k] = f_13 * hsi_252[k]
                   + f_3 * pc_y[k] * isi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_y, pc_y, pc_z, hsk0_327, hsk0_329, \
                         hsi_224, hsi_253, hsi_254, hsk1_327, hsk1_329, isi_364, \
                         isi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * hsi_224[k]
                   + f_3 * pc_z[k] * isi_364[k];

        t_471[k] = pa_y[k] * hsk0_327[k]
                   + f_14 * hsi_253[k]
                   - f_12 * pc_y[k] * hsk1_327[k];

        t_472[k] = f_13 * hsi_254[k]
                   + f_3 * pc_y[k] * isi_366[k];

        t_473[k] = pa_y[k] * hsk0_329[k]
                   - f_12 * pc_y[k] * hsk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_y, pc_y, pc_z, hsk0_330, hsk0_333, \
                         hsi_227, hsi_255, hsi_257, hsk1_330, hsk1_333, isi_367, \
                         isi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_y[k] * hsk0_330[k]
                   + f_15 * hsi_255[k]
                   - f_12 * pc_y[k] * hsk1_330[k];

        t_475[k] = f_15 * hsi_227[k]
                   + f_3 * pc_z[k] * isi_367[k];

        t_476[k] = f_13 * hsi_257[k]
                   + f_3 * pc_y[k] * isi_369[k];

        t_477[k] = pa_y[k] * hsk0_333[k]
                   - f_12 * pc_y[k] * hsk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_y, pc_y, pc_z, hsk0_334, hsk0_336, hsi_230, \
                         hsi_258, hsi_260, hsk1_334, hsk1_336, \
                         isi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pa_y[k] * hsk0_334[k]
                   + f_16 * hsi_258[k]
                   - f_12 * pc_y[k] * hsk1_334[k];

        t_479[k] = f_15 * hsi_230[k]
                   + f_3 * pc_z[k] * isi_370[k];

        t_480[k] = pa_y[k] * hsk0_336[k]
                   + f_14 * hsi_260[k]
                   - f_12 * pc_y[k] * hsk1_336[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
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
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);
    const auto f_20 = 3.5 / q;

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
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsk0_338 = buffer.data(hsk0 + 338);
    const auto *hsk0_339 = buffer.data(hsk0 + 339);
    const auto *hsk0_341 = buffer.data(hsk0 + 341);
    const auto *hsk0_342 = buffer.data(hsk0 + 342);
    const auto *hsk0_344 = buffer.data(hsk0 + 344);
    const auto *hsk0_359 = buffer.data(hsk0 + 359);
    const auto *hsk0_360 = buffer.data(hsk0 + 360);
    const auto *hsk0_363 = buffer.data(hsk0 + 363);
    const auto *hsk0_366 = buffer.data(hsk0 + 366);
    const auto *hsk0_370 = buffer.data(hsk0 + 370);
    const auto *hsk0_375 = buffer.data(hsk0 + 375);
    const auto *hsk0_540 = buffer.data(hsk0 + 540);
    const auto *hsk0_543 = buffer.data(hsk0 + 543);
    const auto *hsk0_546 = buffer.data(hsk0 + 546);
    const auto *hsk0_550 = buffer.data(hsk0 + 550);
    const auto *hsk0_555 = buffer.data(hsk0 + 555);
    const auto *hsk0_568 = buffer.data(hsk0 + 568);
    const auto *hsk0_570 = buffer.data(hsk0 + 570);
    const auto *hsk0_571 = buffer.data(hsk0 + 571);
    const auto *hsk0_572 = buffer.data(hsk0 + 572);
    const auto *hsk0_573 = buffer.data(hsk0 + 573);
    const auto *hsk0_575 = buffer.data(hsk0 + 575);
    const auto *hsk0_581 = buffer.data(hsk0 + 581);
    const auto *hsk0_585 = buffer.data(hsk0 + 585);
    const auto *hsk0_588 = buffer.data(hsk0 + 588);
    const auto *hsk0_590 = buffer.data(hsk0 + 590);
    const auto *hsk0_593 = buffer.data(hsk0 + 593);
    const auto *hsk0_594 = buffer.data(hsk0 + 594);
    const auto *hsk0_596 = buffer.data(hsk0 + 596);

    const auto *hsi_234 = buffer.data(hsi + 234);
    const auto *hsi_245 = buffer.data(hsi + 245);
    const auto *hsi_252 = buffer.data(hsi + 252);
    const auto *hsi_261 = buffer.data(hsi + 261);
    const auto *hsi_262 = buffer.data(hsi + 262);
    const auto *hsi_264 = buffer.data(hsi + 264);
    const auto *hsi_265 = buffer.data(hsi + 265);
    const auto *hsi_266 = buffer.data(hsi + 266);
    const auto *hsi_273 = buffer.data(hsi + 273);
    const auto *hsi_275 = buffer.data(hsi + 275);
    const auto *hsi_276 = buffer.data(hsi + 276);
    const auto *hsi_277 = buffer.data(hsi + 277);
    const auto *hsi_278 = buffer.data(hsi + 278);
    const auto *hsi_279 = buffer.data(hsi + 279);
    const auto *hsi_280 = buffer.data(hsi + 280);
    const auto *hsi_283 = buffer.data(hsi + 283);
    const auto *hsi_285 = buffer.data(hsi + 285);
    const auto *hsi_286 = buffer.data(hsi + 286);
    const auto *hsi_289 = buffer.data(hsi + 289);
    const auto *hsi_290 = buffer.data(hsi + 290);
    const auto *hsi_294 = buffer.data(hsi + 294);
    const auto *hsi_307 = buffer.data(hsi + 307);
    const auto *hsi_308 = buffer.data(hsi + 308);
    const auto *hsi_310 = buffer.data(hsi + 310);
    const auto *hsi_313 = buffer.data(hsi + 313);
    const auto *hsi_317 = buffer.data(hsi + 317);
    const auto *hsi_322 = buffer.data(hsi + 322);
    const auto *hsi_385 = buffer.data(hsi + 385);
    const auto *hsi_386 = buffer.data(hsi + 386);
    const auto *hsi_387 = buffer.data(hsi + 387);
    const auto *hsi_388 = buffer.data(hsi + 388);
    const auto *hsi_389 = buffer.data(hsi + 389);
    const auto *hsi_390 = buffer.data(hsi + 390);
    const auto *hsi_391 = buffer.data(hsi + 391);
    const auto *hsi_392 = buffer.data(hsi + 392);
    const auto *hsi_397 = buffer.data(hsi + 397);
    const auto *hsi_401 = buffer.data(hsi + 401);
    const auto *hsi_406 = buffer.data(hsi + 406);
    const auto *hsi_412 = buffer.data(hsi + 412);
    const auto *hsi_413 = buffer.data(hsi + 413);
    const auto *hsi_414 = buffer.data(hsi + 414);
    const auto *hsi_415 = buffer.data(hsi + 415);
    const auto *hsi_416 = buffer.data(hsi + 416);
    const auto *hsi_417 = buffer.data(hsi + 417);
    const auto *hsi_419 = buffer.data(hsi + 419);
    const auto *hsi_420 = buffer.data(hsi + 420);
    const auto *hsi_423 = buffer.data(hsi + 423);
    const auto *hsi_426 = buffer.data(hsi + 426);
    const auto *hsi_430 = buffer.data(hsi + 430);
    const auto *hsi_435 = buffer.data(hsi + 435);
    const auto *hsi_441 = buffer.data(hsi + 441);
    const auto *hsi_443 = buffer.data(hsi + 443);
    const auto *hsi_444 = buffer.data(hsi + 444);
    const auto *hsi_445 = buffer.data(hsi + 445);
    const auto *hsi_446 = buffer.data(hsi + 446);
    const auto *hsi_447 = buffer.data(hsi + 447);
    const auto *hsi_453 = buffer.data(hsi + 453);
    const auto *hsi_457 = buffer.data(hsi + 457);
    const auto *hsi_460 = buffer.data(hsi + 460);
    const auto *hsi_462 = buffer.data(hsi + 462);
    const auto *hsi_465 = buffer.data(hsi + 465);
    const auto *hsi_466 = buffer.data(hsi + 466);
    const auto *hsi_468 = buffer.data(hsi + 468);
    const auto *hsi_469 = buffer.data(hsi + 469);
    const auto *hsi_470 = buffer.data(hsi + 470);
    const auto *hsi_471 = buffer.data(hsi + 471);
    const auto *hsi_472 = buffer.data(hsi + 472);
    const auto *hsi_473 = buffer.data(hsi + 473);

    const auto *hsk1_338 = buffer.data(hsk1 + 338);
    const auto *hsk1_339 = buffer.data(hsk1 + 339);
    const auto *hsk1_341 = buffer.data(hsk1 + 341);
    const auto *hsk1_342 = buffer.data(hsk1 + 342);
    const auto *hsk1_344 = buffer.data(hsk1 + 344);
    const auto *hsk1_359 = buffer.data(hsk1 + 359);
    const auto *hsk1_360 = buffer.data(hsk1 + 360);
    const auto *hsk1_363 = buffer.data(hsk1 + 363);
    const auto *hsk1_366 = buffer.data(hsk1 + 366);
    const auto *hsk1_370 = buffer.data(hsk1 + 370);
    const auto *hsk1_375 = buffer.data(hsk1 + 375);
    const auto *hsk1_540 = buffer.data(hsk1 + 540);
    const auto *hsk1_543 = buffer.data(hsk1 + 543);
    const auto *hsk1_546 = buffer.data(hsk1 + 546);
    const auto *hsk1_550 = buffer.data(hsk1 + 550);
    const auto *hsk1_555 = buffer.data(hsk1 + 555);
    const auto *hsk1_568 = buffer.data(hsk1 + 568);
    const auto *hsk1_570 = buffer.data(hsk1 + 570);
    const auto *hsk1_571 = buffer.data(hsk1 + 571);
    const auto *hsk1_572 = buffer.data(hsk1 + 572);
    const auto *hsk1_573 = buffer.data(hsk1 + 573);
    const auto *hsk1_575 = buffer.data(hsk1 + 575);
    const auto *hsk1_581 = buffer.data(hsk1 + 581);
    const auto *hsk1_585 = buffer.data(hsk1 + 585);
    const auto *hsk1_588 = buffer.data(hsk1 + 588);
    const auto *hsk1_590 = buffer.data(hsk1 + 590);
    const auto *hsk1_593 = buffer.data(hsk1 + 593);
    const auto *hsk1_594 = buffer.data(hsk1 + 594);
    const auto *hsk1_596 = buffer.data(hsk1 + 596);

    const auto *ish0_288 = buffer.data(ish0 + 288);
    const auto *ish0_290 = buffer.data(ish0 + 290);
    const auto *ish0_291 = buffer.data(ish0 + 291);
    const auto *ish0_292 = buffer.data(ish0 + 292);
    const auto *ish0_293 = buffer.data(ish0 + 293);
    const auto *ish0_294 = buffer.data(ish0 + 294);
    const auto *ish0_295 = buffer.data(ish0 + 295);
    const auto *ish0_296 = buffer.data(ish0 + 296);
    const auto *ish0_297 = buffer.data(ish0 + 297);
    const auto *ish0_298 = buffer.data(ish0 + 298);
    const auto *ish0_299 = buffer.data(ish0 + 299);
    const auto *ish0_300 = buffer.data(ish0 + 300);
    const auto *ish0_301 = buffer.data(ish0 + 301);
    const auto *ish0_302 = buffer.data(ish0 + 302);
    const auto *ish0_303 = buffer.data(ish0 + 303);
    const auto *ish0_308 = buffer.data(ish0 + 308);
    const auto *ish0_309 = buffer.data(ish0 + 309);
    const auto *ish0_310 = buffer.data(ish0 + 310);
    const auto *ish0_311 = buffer.data(ish0 + 311);
    const auto *ish0_312 = buffer.data(ish0 + 312);
    const auto *ish0_313 = buffer.data(ish0 + 313);
    const auto *ish0_314 = buffer.data(ish0 + 314);
    const auto *ish0_315 = buffer.data(ish0 + 315);
    const auto *ish0_317 = buffer.data(ish0 + 317);
    const auto *ish0_318 = buffer.data(ish0 + 318);
    const auto *ish0_320 = buffer.data(ish0 + 320);
    const auto *ish0_321 = buffer.data(ish0 + 321);
    const auto *ish0_322 = buffer.data(ish0 + 322);
    const auto *ish0_324 = buffer.data(ish0 + 324);

    const auto *ish1_288 = buffer.data(ish1 + 288);
    const auto *ish1_290 = buffer.data(ish1 + 290);
    const auto *ish1_291 = buffer.data(ish1 + 291);
    const auto *ish1_292 = buffer.data(ish1 + 292);
    const auto *ish1_293 = buffer.data(ish1 + 293);
    const auto *ish1_294 = buffer.data(ish1 + 294);
    const auto *ish1_295 = buffer.data(ish1 + 295);
    const auto *ish1_296 = buffer.data(ish1 + 296);
    const auto *ish1_297 = buffer.data(ish1 + 297);
    const auto *ish1_298 = buffer.data(ish1 + 298);
    const auto *ish1_299 = buffer.data(ish1 + 299);
    const auto *ish1_300 = buffer.data(ish1 + 300);
    const auto *ish1_301 = buffer.data(ish1 + 301);
    const auto *ish1_302 = buffer.data(ish1 + 302);
    const auto *ish1_303 = buffer.data(ish1 + 303);
    const auto *ish1_308 = buffer.data(ish1 + 308);
    const auto *ish1_309 = buffer.data(ish1 + 309);
    const auto *ish1_310 = buffer.data(ish1 + 310);
    const auto *ish1_311 = buffer.data(ish1 + 311);
    const auto *ish1_312 = buffer.data(ish1 + 312);
    const auto *ish1_313 = buffer.data(ish1 + 313);
    const auto *ish1_314 = buffer.data(ish1 + 314);
    const auto *ish1_315 = buffer.data(ish1 + 315);
    const auto *ish1_317 = buffer.data(ish1 + 317);
    const auto *ish1_318 = buffer.data(ish1 + 318);
    const auto *ish1_320 = buffer.data(ish1 + 320);
    const auto *ish1_321 = buffer.data(ish1 + 321);
    const auto *ish1_322 = buffer.data(ish1 + 322);
    const auto *ish1_324 = buffer.data(ish1 + 324);

    const auto *isi_373 = buffer.data(isi + 373);
    const auto *isi_374 = buffer.data(isi + 374);
    const auto *isi_378 = buffer.data(isi + 378);
    const auto *isi_385 = buffer.data(isi + 385);
    const auto *isi_386 = buffer.data(isi + 386);
    const auto *isi_387 = buffer.data(isi + 387);
    const auto *isi_388 = buffer.data(isi + 388);
    const auto *isi_389 = buffer.data(isi + 389);
    const auto *isi_390 = buffer.data(isi + 390);
    const auto *isi_391 = buffer.data(isi + 391);
    const auto *isi_392 = buffer.data(isi + 392);
    const auto *isi_393 = buffer.data(isi + 393);
    const auto *isi_394 = buffer.data(isi + 394);
    const auto *isi_395 = buffer.data(isi + 395);
    const auto *isi_396 = buffer.data(isi + 396);
    const auto *isi_397 = buffer.data(isi + 397);
    const auto *isi_398 = buffer.data(isi + 398);
    const auto *isi_399 = buffer.data(isi + 399);
    const auto *isi_400 = buffer.data(isi + 400);
    const auto *isi_401 = buffer.data(isi + 401);
    const auto *isi_402 = buffer.data(isi + 402);
    const auto *isi_403 = buffer.data(isi + 403);
    const auto *isi_404 = buffer.data(isi + 404);
    const auto *isi_405 = buffer.data(isi + 405);
    const auto *isi_406 = buffer.data(isi + 406);
    const auto *isi_412 = buffer.data(isi + 412);
    const auto *isi_413 = buffer.data(isi + 413);
    const auto *isi_414 = buffer.data(isi + 414);
    const auto *isi_415 = buffer.data(isi + 415);
    const auto *isi_416 = buffer.data(isi + 416);
    const auto *isi_417 = buffer.data(isi + 417);
    const auto *isi_418 = buffer.data(isi + 418);
    const auto *isi_419 = buffer.data(isi + 419);
    const auto *isi_420 = buffer.data(isi + 420);
    const auto *isi_421 = buffer.data(isi + 421);
    const auto *isi_422 = buffer.data(isi + 422);
    const auto *isi_423 = buffer.data(isi + 423);
    const auto *isi_425 = buffer.data(isi + 425);
    const auto *isi_426 = buffer.data(isi + 426);
    const auto *isi_427 = buffer.data(isi + 427);
    const auto *isi_429 = buffer.data(isi + 429);
    const auto *isi_430 = buffer.data(isi + 430);
    const auto *isi_431 = buffer.data(isi + 431);
    const auto *isi_432 = buffer.data(isi + 432);
    const auto *isi_434 = buffer.data(isi + 434);
    const auto *isi_435 = buffer.data(isi + 435);
    const auto *isi_441 = buffer.data(isi + 441);
    const auto *isi_443 = buffer.data(isi + 443);
    const auto *isi_444 = buffer.data(isi + 444);
    const auto *isi_445 = buffer.data(isi + 445);
    const auto *isi_446 = buffer.data(isi + 446);
    const auto *isi_447 = buffer.data(isi + 447);
    const auto *isi_448 = buffer.data(isi + 448);
    const auto *isi_450 = buffer.data(isi + 450);
    const auto *isi_451 = buffer.data(isi + 451);
    const auto *isi_453 = buffer.data(isi + 453);
    const auto *isi_454 = buffer.data(isi + 454);
    const auto *isi_457 = buffer.data(isi + 457);
    const auto *isi_458 = buffer.data(isi + 458);
    const auto *isi_462 = buffer.data(isi + 462);
    const auto *isi_469 = buffer.data(isi + 469);
    const auto *isi_470 = buffer.data(isi + 470);
    const auto *isi_471 = buffer.data(isi + 471);
    const auto *isi_472 = buffer.data(isi + 472);
    const auto *isi_473 = buffer.data(isi + 473);

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_y, pc_y, pc_z, hsk0_338, hsk0_339, \
                         hsi_234, hsi_261, hsi_262, hsk1_338, hsk1_339, isi_373, \
                         isi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * hsi_261[k]
                   + f_3 * pc_y[k] * isi_373[k];

        t_482[k] = pa_y[k] * hsk0_338[k]
                   - f_12 * pc_y[k] * hsk1_338[k];

        t_483[k] = pa_y[k] * hsk0_339[k]
                   + f_17 * hsi_262[k]
                   - f_12 * pc_y[k] * hsk1_339[k];

        t_484[k] = f_15 * hsi_234[k]
                   + f_3 * pc_z[k] * isi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_y, pc_y, hsk0_341, hsk0_342, hsk0_344, \
                         hsi_264, hsi_265, hsi_266, hsk1_341, hsk1_342, hsk1_344, \
                         isi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_y[k] * hsk0_341[k]
                   + f_15 * hsi_264[k]
                   - f_12 * pc_y[k] * hsk1_341[k];

        t_486[k] = pa_y[k] * hsk0_342[k]
                   + f_14 * hsi_265[k]
                   - f_12 * pc_y[k] * hsk1_342[k];

        t_487[k] = f_13 * hsi_266[k]
                   + f_3 * pc_y[k] * isi_378[k];

        t_488[k] = pa_y[k] * hsk0_344[k]
                   - f_12 * pc_y[k] * hsk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, hsi_385, hsi_386, hsi_387, \
                         hsi_388, hsi_389, isi_385, isi_386, isi_387, isi_388, \
                         isi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_14 * hsi_385[k]
                   + f_3 * pc_x[k] * isi_385[k];

        t_490[k] = f_14 * hsi_386[k]
                   + f_3 * pc_x[k] * isi_386[k];

        t_491[k] = f_14 * hsi_387[k]
                   + f_3 * pc_x[k] * isi_387[k];

        t_492[k] = f_14 * hsi_388[k]
                   + f_3 * pc_x[k] * isi_388[k];

        t_493[k] = f_14 * hsi_389[k]
                   + f_3 * pc_x[k] * isi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, hsi_245, hsi_273, \
                         hsi_390, hsi_391, ish0_288, ish1_288, isi_385, isi_390, \
                         isi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_14 * hsi_390[k]
                   + f_3 * pc_x[k] * isi_390[k];

        t_495[k] = f_14 * hsi_391[k]
                   + f_3 * pc_x[k] * isi_391[k];

        t_496[k] = f_13 * hsi_273[k]
                   + f_1 * ish0_288[k]
                   - f_2 * ish1_288[k]
                   + f_3 * pc_y[k] * isi_385[k];

        t_497[k] = f_15 * hsi_245[k]
                   + f_3 * pc_z[k] * isi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, hsi_275, hsi_276, hsi_277, ish0_290, \
                         ish0_291, ish0_292, ish1_290, ish1_291, ish1_292, isi_387, isi_388, \
                         isi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * hsi_275[k]
                   + f_10 * ish0_290[k]
                   - f_11 * ish1_290[k]
                   + f_3 * pc_y[k] * isi_387[k];

        t_499[k] = f_13 * hsi_276[k]
                   + f_8 * ish0_291[k]
                   - f_9 * ish1_291[k]
                   + f_3 * pc_y[k] * isi_388[k];

        t_500[k] = f_13 * hsi_277[k]
                   + f_6 * ish0_292[k]
                   - f_7 * ish1_292[k]
                   + f_3 * pc_y[k] * isi_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pc_y, hsk0_359, hsi_278, hsi_279, \
                         hsk1_359, ish0_293, ish1_293, isi_390, \
                         isi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * hsi_278[k]
                   + f_4 * ish0_293[k]
                   - f_5 * ish1_293[k]
                   + f_3 * pc_y[k] * isi_390[k];

        t_502[k] = f_13 * hsi_279[k]
                   + f_3 * pc_y[k] * isi_391[k];

        t_503[k] = pa_y[k] * hsk0_359[k]
                   - f_12 * pc_y[k] * hsk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, hsi_252, \
                         hsi_392, ish0_294, ish1_294, isi_392, isi_393, \
                         isi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_14 * hsi_392[k]
                   + f_1 * ish0_294[k]
                   - f_2 * ish1_294[k]
                   + f_3 * pc_x[k] * isi_392[k];

        t_505[k] = f_3 * pc_y[k] * isi_392[k];

        t_506[k] = f_16 * hsi_252[k]
                   + f_3 * pc_z[k] * isi_392[k];

        t_507[k] = f_4 * ish0_294[k]
                   - f_5 * ish1_294[k]
                   + f_3 * pc_y[k] * isi_393[k];

        t_508[k] = f_3 * pc_y[k] * isi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, hsi_397, ish0_295, ish0_296, \
                         ish0_299, ish1_295, ish1_296, ish1_299, isi_395, isi_396, \
                         isi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_14 * hsi_397[k]
                   + f_10 * ish0_299[k]
                   - f_11 * ish1_299[k]
                   + f_3 * pc_x[k] * isi_397[k];

        t_510[k] = f_6 * ish0_295[k]
                   - f_7 * ish1_295[k]
                   + f_3 * pc_y[k] * isi_395[k];

        t_511[k] = f_4 * ish0_296[k]
                   - f_5 * ish1_296[k]
                   + f_3 * pc_y[k] * isi_396[k];

        t_512[k] = f_3 * pc_y[k] * isi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_y, hsi_401, ish0_297, ish0_298, \
                         ish0_303, ish1_297, ish1_298, ish1_303, isi_398, isi_399, \
                         isi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_14 * hsi_401[k]
                   + f_8 * ish0_303[k]
                   - f_9 * ish1_303[k]
                   + f_3 * pc_x[k] * isi_401[k];

        t_514[k] = f_8 * ish0_297[k]
                   - f_9 * ish1_297[k]
                   + f_3 * pc_y[k] * isi_398[k];

        t_515[k] = f_6 * ish0_298[k]
                   - f_7 * ish1_298[k]
                   + f_3 * pc_y[k] * isi_399[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_x, pc_y, hsi_406, ish0_299, ish0_308, \
                         ish1_299, ish1_308, isi_400, isi_401, \
                         isi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_4 * ish0_299[k]
                   - f_5 * ish1_299[k]
                   + f_3 * pc_y[k] * isi_400[k];

        t_517[k] = f_3 * pc_y[k] * isi_401[k];

        t_518[k] = f_14 * hsi_406[k]
                   + f_6 * ish0_308[k]
                   - f_7 * ish1_308[k]
                   + f_3 * pc_x[k] * isi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, ish0_300, ish0_301, ish0_302, ish1_300, \
                         ish1_301, ish1_302, isi_402, isi_403, \
                         isi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_10 * ish0_300[k]
                   - f_11 * ish1_300[k]
                   + f_3 * pc_y[k] * isi_402[k];

        t_520[k] = f_8 * ish0_301[k]
                   - f_9 * ish1_301[k]
                   + f_3 * pc_y[k] * isi_403[k];

        t_521[k] = f_6 * ish0_302[k]
                   - f_7 * ish1_302[k]
                   + f_3 * pc_y[k] * isi_404[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, hsi_412, hsi_413, ish0_303, \
                         ish0_314, ish1_303, ish1_314, isi_405, isi_406, isi_412, \
                         isi_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_4 * ish0_303[k]
                   - f_5 * ish1_303[k]
                   + f_3 * pc_y[k] * isi_405[k];

        t_523[k] = f_3 * pc_y[k] * isi_406[k];

        t_524[k] = f_14 * hsi_412[k]
                   + f_4 * ish0_314[k]
                   - f_5 * ish1_314[k]
                   + f_3 * pc_x[k] * isi_412[k];

        t_525[k] = f_14 * hsi_413[k]
                   + f_3 * pc_x[k] * isi_413[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, pc_x, pc_y, hsi_414, hsi_415, \
                         hsi_416, hsi_417, isi_412, isi_414, isi_415, isi_416, \
                         isi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_14 * hsi_414[k]
                   + f_3 * pc_x[k] * isi_414[k];

        t_527[k] = f_14 * hsi_415[k]
                   + f_3 * pc_x[k] * isi_415[k];

        t_528[k] = f_14 * hsi_416[k]
                   + f_3 * pc_x[k] * isi_416[k];

        t_529[k] = f_14 * hsi_417[k]
                   + f_3 * pc_x[k] * isi_417[k];

        t_530[k] = f_3 * pc_y[k] * isi_412[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pc_x, pc_y, hsi_419, ish0_309, ish0_310, \
                         ish1_309, ish1_310, isi_413, isi_414, \
                         isi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_14 * hsi_419[k]
                   + f_3 * pc_x[k] * isi_419[k];

        t_532[k] = f_1 * ish0_309[k]
                   - f_2 * ish1_309[k]
                   + f_3 * pc_y[k] * isi_413[k];

        t_533[k] = f_18 * ish0_310[k]
                   - f_19 * ish1_310[k]
                   + f_3 * pc_y[k] * isi_414[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, ish0_311, ish0_312, ish0_313, ish1_311, \
                         ish1_312, ish1_313, isi_415, isi_416, \
                         isi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * ish0_311[k]
                   - f_11 * ish1_311[k]
                   + f_3 * pc_y[k] * isi_415[k];

        t_535[k] = f_8 * ish0_312[k]
                   - f_9 * ish1_312[k]
                   + f_3 * pc_y[k] * isi_416[k];

        t_536[k] = f_6 * ish0_313[k]
                   - f_7 * ish1_313[k]
                   + f_3 * pc_y[k] * isi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_x, pc_x, pc_y, pc_z, hsk0_540, \
                         hsi_279, hsi_420, hsk1_540, ish0_314, ish1_314, isi_418, \
                         isi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * ish0_314[k]
                   - f_5 * ish1_314[k]
                   + f_3 * pc_y[k] * isi_418[k];

        t_538[k] = f_3 * pc_y[k] * isi_419[k];

        t_539[k] = f_16 * hsi_279[k]
                   + f_1 * ish0_314[k]
                   - f_2 * ish1_314[k]
                   + f_3 * pc_z[k] * isi_419[k];

        t_540[k] = pa_x[k] * hsk0_540[k]
                   + f_20 * hsi_420[k]
                   - f_12 * pc_x[k] * hsk1_540[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pa_x, pc_x, pc_y, pc_z, hsk0_543, \
                         hsi_280, hsi_423, hsk1_543, isi_420, isi_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_17 * hsi_280[k]
                   + f_3 * pc_y[k] * isi_420[k];

        t_542[k] = f_3 * pc_z[k] * isi_420[k];

        t_543[k] = pa_x[k] * hsk0_543[k]
                   + f_17 * hsi_423[k]
                   - f_12 * pc_x[k] * hsk1_543[k];

        t_544[k] = f_3 * pc_z[k] * isi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pa_x, pc_x, pc_z, hsk0_546, hsi_426, hsk1_546, \
                         ish0_315, ish1_315, isi_422, isi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * ish0_315[k]
                   - f_5 * ish1_315[k]
                   + f_3 * pc_z[k] * isi_422[k];

        t_546[k] = pa_x[k] * hsk0_546[k]
                   + f_16 * hsi_426[k]
                   - f_12 * pc_x[k] * hsk1_546[k];

        t_547[k] = f_3 * pc_z[k] * isi_423[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pa_x, pc_x, pc_y, pc_z, hsk0_550, \
                         hsi_285, hsi_430, hsk1_550, ish0_317, ish1_317, isi_425, \
                         isi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * hsi_285[k]
                   + f_3 * pc_y[k] * isi_425[k];

        t_549[k] = f_6 * ish0_317[k]
                   - f_7 * ish1_317[k]
                   + f_3 * pc_z[k] * isi_425[k];

        t_550[k] = pa_x[k] * hsk0_550[k]
                   + f_15 * hsi_430[k]
                   - f_12 * pc_x[k] * hsk1_550[k];

        t_551[k] = f_3 * pc_z[k] * isi_426[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_y, pc_z, hsi_289, ish0_318, ish0_320, \
                         ish1_318, ish1_320, isi_427, isi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_4 * ish0_318[k]
                   - f_5 * ish1_318[k]
                   + f_3 * pc_z[k] * isi_427[k];

        t_553[k] = f_17 * hsi_289[k]
                   + f_3 * pc_y[k] * isi_429[k];

        t_554[k] = f_8 * ish0_320[k]
                   - f_9 * ish1_320[k]
                   + f_3 * pc_z[k] * isi_429[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pa_x, pc_x, pc_z, hsk0_555, hsi_435, hsk1_555, \
                         ish0_321, ish1_321, isi_430, isi_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_x[k] * hsk0_555[k]
                   + f_14 * hsi_435[k]
                   - f_12 * pc_x[k] * hsk1_555[k];

        t_556[k] = f_3 * pc_z[k] * isi_430[k];

        t_557[k] = f_4 * ish0_321[k]
                   - f_5 * ish1_321[k]
                   + f_3 * pc_z[k] * isi_431[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pc_x, pc_y, pc_z, hsi_294, hsi_441, \
                         ish0_322, ish0_324, ish1_322, ish1_324, isi_432, isi_434, \
                         isi_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_6 * ish0_322[k]
                   - f_7 * ish1_322[k]
                   + f_3 * pc_z[k] * isi_432[k];

        t_559[k] = f_17 * hsi_294[k]
                   + f_3 * pc_y[k] * isi_434[k];

        t_560[k] = f_10 * ish0_324[k]
                   - f_11 * ish1_324[k]
                   + f_3 * pc_z[k] * isi_434[k];

        t_561[k] = f_13 * hsi_441[k]
                   + f_3 * pc_x[k] * isi_441[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pc_x, pc_z, hsi_443, hsi_444, \
                         hsi_445, hsi_446, isi_435, isi_443, isi_444, isi_445, \
                         isi_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_3 * pc_z[k] * isi_435[k];

        t_563[k] = f_13 * hsi_443[k]
                   + f_3 * pc_x[k] * isi_443[k];

        t_564[k] = f_13 * hsi_444[k]
                   + f_3 * pc_x[k] * isi_444[k];

        t_565[k] = f_13 * hsi_445[k]
                   + f_3 * pc_x[k] * isi_445[k];

        t_566[k] = f_13 * hsi_446[k]
                   + f_3 * pc_x[k] * isi_446[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pa_x, pc_x, pc_z, hsk0_568, hsk0_570, \
                         hsi_447, hsk1_568, hsk1_570, isi_441, \
                         isi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_13 * hsi_447[k]
                   + f_3 * pc_x[k] * isi_447[k];

        t_568[k] = pa_x[k] * hsk0_568[k]
                   - f_12 * pc_x[k] * hsk1_568[k];

        t_569[k] = f_3 * pc_z[k] * isi_441[k];

        t_570[k] = pa_x[k] * hsk0_570[k]
                   - f_12 * pc_x[k] * hsk1_570[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, pa_x, pc_x, pc_y, hsk0_571, hsk0_572, \
                         hsk0_573, hsi_307, hsk1_571, hsk1_572, hsk1_573, \
                         isi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = pa_x[k] * hsk0_571[k]
                   - f_12 * pc_x[k] * hsk1_571[k];

        t_572[k] = pa_x[k] * hsk0_572[k]
                   - f_12 * pc_x[k] * hsk1_572[k];

        t_573[k] = pa_x[k] * hsk0_573[k]
                   - f_12 * pc_x[k] * hsk1_573[k];

        t_574[k] = f_17 * hsi_307[k]
                   + f_3 * pc_y[k] * isi_447[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pa_x, pa_z, pc_x, pc_y, pc_z, hsk0_360, \
                         hsk0_575, hsi_280, hsi_308, hsk1_360, hsk1_575, \
                         isi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = pa_x[k] * hsk0_575[k]
                   - f_12 * pc_x[k] * hsk1_575[k];

        t_576[k] = pa_z[k] * hsk0_360[k]
                   - f_12 * pc_z[k] * hsk1_360[k];

        t_577[k] = f_16 * hsi_308[k]
                   + f_3 * pc_y[k] * isi_448[k];

        t_578[k] = f_13 * hsi_280[k]
                   + f_3 * pc_z[k] * isi_448[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pa_x, pa_z, pc_x, pc_y, pc_z, hsk0_363, \
                         hsk0_581, hsi_310, hsi_453, hsk1_363, hsk1_581, \
                         isi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = pa_z[k] * hsk0_363[k]
                   - f_12 * pc_z[k] * hsk1_363[k];

        t_580[k] = f_16 * hsi_310[k]
                   + f_3 * pc_y[k] * isi_450[k];

        t_581[k] = pa_x[k] * hsk0_581[k]
                   + f_17 * hsi_453[k]
                   - f_12 * pc_x[k] * hsk1_581[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pa_z, pc_y, pc_z, hsk0_366, hsi_283, hsi_313, \
                         hsk1_366, isi_451, isi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = pa_z[k] * hsk0_366[k]
                   - f_12 * pc_z[k] * hsk1_366[k];

        t_583[k] = f_13 * hsi_283[k]
                   + f_3 * pc_z[k] * isi_451[k];

        t_584[k] = f_16 * hsi_313[k]
                   + f_3 * pc_y[k] * isi_453[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pa_x, pa_z, pc_x, pc_z, hsk0_370, hsk0_585, \
                         hsi_286, hsi_457, hsk1_370, hsk1_585, \
                         isi_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = pa_x[k] * hsk0_585[k]
                   + f_16 * hsi_457[k]
                   - f_12 * pc_x[k] * hsk1_585[k];

        t_586[k] = pa_z[k] * hsk0_370[k]
                   - f_12 * pc_z[k] * hsk1_370[k];

        t_587[k] = f_13 * hsi_286[k]
                   + f_3 * pc_z[k] * isi_454[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pa_x, pc_x, pc_y, hsk0_588, hsk0_590, hsi_317, \
                         hsi_460, hsi_462, hsk1_588, hsk1_590, \
                         isi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_x[k] * hsk0_588[k]
                   + f_15 * hsi_460[k]
                   - f_12 * pc_x[k] * hsk1_588[k];

        t_589[k] = f_16 * hsi_317[k]
                   + f_3 * pc_y[k] * isi_457[k];

        t_590[k] = pa_x[k] * hsk0_590[k]
                   + f_15 * hsi_462[k]
                   - f_12 * pc_x[k] * hsk1_590[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, pa_x, pa_z, pc_x, pc_z, hsk0_375, hsk0_593, \
                         hsi_290, hsi_465, hsk1_375, hsk1_593, \
                         isi_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = pa_z[k] * hsk0_375[k]
                   - f_12 * pc_z[k] * hsk1_375[k];

        t_592[k] = f_13 * hsi_290[k]
                   + f_3 * pc_z[k] * isi_458[k];

        t_593[k] = pa_x[k] * hsk0_593[k]
                   + f_14 * hsi_465[k]
                   - f_12 * pc_x[k] * hsk1_593[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pa_x, pc_x, pc_y, hsk0_594, hsk0_596, hsi_322, \
                         hsi_466, hsi_468, hsk1_594, hsk1_596, \
                         isi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = pa_x[k] * hsk0_594[k]
                   + f_14 * hsi_466[k]
                   - f_12 * pc_x[k] * hsk1_594[k];

        t_595[k] = f_16 * hsi_322[k]
                   + f_3 * pc_y[k] * isi_462[k];

        t_596[k] = pa_x[k] * hsk0_596[k]
                   + f_14 * hsi_468[k]
                   - f_12 * pc_x[k] * hsk1_596[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, pc_x, hsi_469, hsi_470, hsi_471, \
                         hsi_472, hsi_473, isi_469, isi_470, isi_471, isi_472, \
                         isi_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_13 * hsi_469[k]
                   + f_3 * pc_x[k] * isi_469[k];

        t_598[k] = f_13 * hsi_470[k]
                   + f_3 * pc_x[k] * isi_470[k];

        t_599[k] = f_13 * hsi_471[k]
                   + f_3 * pc_x[k] * isi_471[k];

        t_600[k] = f_13 * hsi_472[k]
                   + f_3 * pc_x[k] * isi_472[k];

        t_601[k] = f_13 * hsi_473[k]
                   + f_3 * pc_x[k] * isi_473[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t isi, const size_t ncols,
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
    const auto f_20 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsk0_504 = buffer.data(hsk0 + 504);
    const auto *hsk0_509 = buffer.data(hsk0 + 509);
    const auto *hsk0_513 = buffer.data(hsk0 + 513);
    const auto *hsk0_518 = buffer.data(hsk0 + 518);
    const auto *hsk0_524 = buffer.data(hsk0 + 524);
    const auto *hsk0_604 = buffer.data(hsk0 + 604);
    const auto *hsk0_606 = buffer.data(hsk0 + 606);
    const auto *hsk0_607 = buffer.data(hsk0 + 607);
    const auto *hsk0_608 = buffer.data(hsk0 + 608);
    const auto *hsk0_609 = buffer.data(hsk0 + 609);
    const auto *hsk0_611 = buffer.data(hsk0 + 611);
    const auto *hsk0_612 = buffer.data(hsk0 + 612);
    const auto *hsk0_615 = buffer.data(hsk0 + 615);
    const auto *hsk0_617 = buffer.data(hsk0 + 617);
    const auto *hsk0_618 = buffer.data(hsk0 + 618);
    const auto *hsk0_621 = buffer.data(hsk0 + 621);
    const auto *hsk0_622 = buffer.data(hsk0 + 622);
    const auto *hsk0_624 = buffer.data(hsk0 + 624);
    const auto *hsk0_626 = buffer.data(hsk0 + 626);
    const auto *hsk0_627 = buffer.data(hsk0 + 627);
    const auto *hsk0_629 = buffer.data(hsk0 + 629);
    const auto *hsk0_630 = buffer.data(hsk0 + 630);
    const auto *hsk0_632 = buffer.data(hsk0 + 632);
    const auto *hsk0_640 = buffer.data(hsk0 + 640);
    const auto *hsk0_642 = buffer.data(hsk0 + 642);
    const auto *hsk0_643 = buffer.data(hsk0 + 643);
    const auto *hsk0_644 = buffer.data(hsk0 + 644);
    const auto *hsk0_645 = buffer.data(hsk0 + 645);
    const auto *hsk0_647 = buffer.data(hsk0 + 647);
    const auto *hsk0_648 = buffer.data(hsk0 + 648);
    const auto *hsk0_651 = buffer.data(hsk0 + 651);
    const auto *hsk0_653 = buffer.data(hsk0 + 653);
    const auto *hsk0_654 = buffer.data(hsk0 + 654);
    const auto *hsk0_657 = buffer.data(hsk0 + 657);
    const auto *hsk0_658 = buffer.data(hsk0 + 658);
    const auto *hsk0_660 = buffer.data(hsk0 + 660);
    const auto *hsk0_662 = buffer.data(hsk0 + 662);
    const auto *hsk0_663 = buffer.data(hsk0 + 663);
    const auto *hsk0_665 = buffer.data(hsk0 + 665);
    const auto *hsk0_666 = buffer.data(hsk0 + 666);
    const auto *hsk0_668 = buffer.data(hsk0 + 668);
    const auto *hsk0_676 = buffer.data(hsk0 + 676);
    const auto *hsk0_678 = buffer.data(hsk0 + 678);
    const auto *hsk0_679 = buffer.data(hsk0 + 679);
    const auto *hsk0_680 = buffer.data(hsk0 + 680);
    const auto *hsk0_681 = buffer.data(hsk0 + 681);
    const auto *hsk0_683 = buffer.data(hsk0 + 683);
    const auto *hsk0_687 = buffer.data(hsk0 + 687);
    const auto *hsk0_690 = buffer.data(hsk0 + 690);
    const auto *hsk0_694 = buffer.data(hsk0 + 694);
    const auto *hsk0_696 = buffer.data(hsk0 + 696);
    const auto *hsk0_699 = buffer.data(hsk0 + 699);
    const auto *hsk0_701 = buffer.data(hsk0 + 701);
    const auto *hsk0_702 = buffer.data(hsk0 + 702);
    const auto *hsk0_712 = buffer.data(hsk0 + 712);
    const auto *hsk0_714 = buffer.data(hsk0 + 714);
    const auto *hsk0_715 = buffer.data(hsk0 + 715);
    const auto *hsk0_716 = buffer.data(hsk0 + 716);
    const auto *hsk0_717 = buffer.data(hsk0 + 717);
    const auto *hsk0_719 = buffer.data(hsk0 + 719);

    const auto *hsi_301 = buffer.data(hsi + 301);
    const auto *hsi_308 = buffer.data(hsi + 308);
    const auto *hsi_311 = buffer.data(hsi + 311);
    const auto *hsi_314 = buffer.data(hsi + 314);
    const auto *hsi_318 = buffer.data(hsi + 318);
    const auto *hsi_329 = buffer.data(hsi + 329);
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
    const auto *hsi_363 = buffer.data(hsi + 363);
    const auto *hsi_364 = buffer.data(hsi + 364);
    const auto *hsi_366 = buffer.data(hsi + 366);
    const auto *hsi_367 = buffer.data(hsi + 367);
    const auto *hsi_369 = buffer.data(hsi + 369);
    const auto *hsi_370 = buffer.data(hsi + 370);
    const auto *hsi_373 = buffer.data(hsi + 373);
    const auto *hsi_374 = buffer.data(hsi + 374);
    const auto *hsi_378 = buffer.data(hsi + 378);
    const auto *hsi_385 = buffer.data(hsi + 385);
    const auto *hsi_391 = buffer.data(hsi + 391);
    const auto *hsi_392 = buffer.data(hsi + 392);
    const auto *hsi_394 = buffer.data(hsi + 394);
    const auto *hsi_397 = buffer.data(hsi + 397);
    const auto *hsi_401 = buffer.data(hsi + 401);
    const auto *hsi_406 = buffer.data(hsi + 406);
    const auto *hsi_419 = buffer.data(hsi + 419);
    const auto *hsi_474 = buffer.data(hsi + 474);
    const auto *hsi_475 = buffer.data(hsi + 475);
    const auto *hsi_476 = buffer.data(hsi + 476);
    const auto *hsi_479 = buffer.data(hsi + 479);
    const auto *hsi_481 = buffer.data(hsi + 481);
    const auto *hsi_482 = buffer.data(hsi + 482);
    const auto *hsi_485 = buffer.data(hsi + 485);
    const auto *hsi_486 = buffer.data(hsi + 486);
    const auto *hsi_488 = buffer.data(hsi + 488);
    const auto *hsi_490 = buffer.data(hsi + 490);
    const auto *hsi_491 = buffer.data(hsi + 491);
    const auto *hsi_493 = buffer.data(hsi + 493);
    const auto *hsi_494 = buffer.data(hsi + 494);
    const auto *hsi_496 = buffer.data(hsi + 496);
    const auto *hsi_497 = buffer.data(hsi + 497);
    const auto *hsi_498 = buffer.data(hsi + 498);
    const auto *hsi_499 = buffer.data(hsi + 499);
    const auto *hsi_500 = buffer.data(hsi + 500);
    const auto *hsi_501 = buffer.data(hsi + 501);
    const auto *hsi_502 = buffer.data(hsi + 502);
    const auto *hsi_503 = buffer.data(hsi + 503);
    const auto *hsi_504 = buffer.data(hsi + 504);
    const auto *hsi_507 = buffer.data(hsi + 507);
    const auto *hsi_509 = buffer.data(hsi + 509);
    const auto *hsi_510 = buffer.data(hsi + 510);
    const auto *hsi_513 = buffer.data(hsi + 513);
    const auto *hsi_514 = buffer.data(hsi + 514);
    const auto *hsi_516 = buffer.data(hsi + 516);
    const auto *hsi_518 = buffer.data(hsi + 518);
    const auto *hsi_519 = buffer.data(hsi + 519);
    const auto *hsi_521 = buffer.data(hsi + 521);
    const auto *hsi_522 = buffer.data(hsi + 522);
    const auto *hsi_524 = buffer.data(hsi + 524);
    const auto *hsi_525 = buffer.data(hsi + 525);
    const auto *hsi_526 = buffer.data(hsi + 526);
    const auto *hsi_527 = buffer.data(hsi + 527);
    const auto *hsi_528 = buffer.data(hsi + 528);
    const auto *hsi_529 = buffer.data(hsi + 529);
    const auto *hsi_530 = buffer.data(hsi + 530);
    const auto *hsi_531 = buffer.data(hsi + 531);
    const auto *hsi_535 = buffer.data(hsi + 535);
    const auto *hsi_538 = buffer.data(hsi + 538);
    const auto *hsi_542 = buffer.data(hsi + 542);
    const auto *hsi_544 = buffer.data(hsi + 544);
    const auto *hsi_547 = buffer.data(hsi + 547);
    const auto *hsi_549 = buffer.data(hsi + 549);
    const auto *hsi_550 = buffer.data(hsi + 550);
    const auto *hsi_553 = buffer.data(hsi + 553);
    const auto *hsi_554 = buffer.data(hsi + 554);
    const auto *hsi_555 = buffer.data(hsi + 555);
    const auto *hsi_556 = buffer.data(hsi + 556);
    const auto *hsi_557 = buffer.data(hsi + 557);
    const auto *hsi_558 = buffer.data(hsi + 558);
    const auto *hsi_559 = buffer.data(hsi + 559);

    const auto *hsk1_504 = buffer.data(hsk1 + 504);
    const auto *hsk1_509 = buffer.data(hsk1 + 509);
    const auto *hsk1_513 = buffer.data(hsk1 + 513);
    const auto *hsk1_518 = buffer.data(hsk1 + 518);
    const auto *hsk1_524 = buffer.data(hsk1 + 524);
    const auto *hsk1_604 = buffer.data(hsk1 + 604);
    const auto *hsk1_606 = buffer.data(hsk1 + 606);
    const auto *hsk1_607 = buffer.data(hsk1 + 607);
    const auto *hsk1_608 = buffer.data(hsk1 + 608);
    const auto *hsk1_609 = buffer.data(hsk1 + 609);
    const auto *hsk1_611 = buffer.data(hsk1 + 611);
    const auto *hsk1_612 = buffer.data(hsk1 + 612);
    const auto *hsk1_615 = buffer.data(hsk1 + 615);
    const auto *hsk1_617 = buffer.data(hsk1 + 617);
    const auto *hsk1_618 = buffer.data(hsk1 + 618);
    const auto *hsk1_621 = buffer.data(hsk1 + 621);
    const auto *hsk1_622 = buffer.data(hsk1 + 622);
    const auto *hsk1_624 = buffer.data(hsk1 + 624);
    const auto *hsk1_626 = buffer.data(hsk1 + 626);
    const auto *hsk1_627 = buffer.data(hsk1 + 627);
    const auto *hsk1_629 = buffer.data(hsk1 + 629);
    const auto *hsk1_630 = buffer.data(hsk1 + 630);
    const auto *hsk1_632 = buffer.data(hsk1 + 632);
    const auto *hsk1_640 = buffer.data(hsk1 + 640);
    const auto *hsk1_642 = buffer.data(hsk1 + 642);
    const auto *hsk1_643 = buffer.data(hsk1 + 643);
    const auto *hsk1_644 = buffer.data(hsk1 + 644);
    const auto *hsk1_645 = buffer.data(hsk1 + 645);
    const auto *hsk1_647 = buffer.data(hsk1 + 647);
    const auto *hsk1_648 = buffer.data(hsk1 + 648);
    const auto *hsk1_651 = buffer.data(hsk1 + 651);
    const auto *hsk1_653 = buffer.data(hsk1 + 653);
    const auto *hsk1_654 = buffer.data(hsk1 + 654);
    const auto *hsk1_657 = buffer.data(hsk1 + 657);
    const auto *hsk1_658 = buffer.data(hsk1 + 658);
    const auto *hsk1_660 = buffer.data(hsk1 + 660);
    const auto *hsk1_662 = buffer.data(hsk1 + 662);
    const auto *hsk1_663 = buffer.data(hsk1 + 663);
    const auto *hsk1_665 = buffer.data(hsk1 + 665);
    const auto *hsk1_666 = buffer.data(hsk1 + 666);
    const auto *hsk1_668 = buffer.data(hsk1 + 668);
    const auto *hsk1_676 = buffer.data(hsk1 + 676);
    const auto *hsk1_678 = buffer.data(hsk1 + 678);
    const auto *hsk1_679 = buffer.data(hsk1 + 679);
    const auto *hsk1_680 = buffer.data(hsk1 + 680);
    const auto *hsk1_681 = buffer.data(hsk1 + 681);
    const auto *hsk1_683 = buffer.data(hsk1 + 683);
    const auto *hsk1_687 = buffer.data(hsk1 + 687);
    const auto *hsk1_690 = buffer.data(hsk1 + 690);
    const auto *hsk1_694 = buffer.data(hsk1 + 694);
    const auto *hsk1_696 = buffer.data(hsk1 + 696);
    const auto *hsk1_699 = buffer.data(hsk1 + 699);
    const auto *hsk1_701 = buffer.data(hsk1 + 701);
    const auto *hsk1_702 = buffer.data(hsk1 + 702);
    const auto *hsk1_712 = buffer.data(hsk1 + 712);
    const auto *hsk1_714 = buffer.data(hsk1 + 714);
    const auto *hsk1_715 = buffer.data(hsk1 + 715);
    const auto *hsk1_716 = buffer.data(hsk1 + 716);
    const auto *hsk1_717 = buffer.data(hsk1 + 717);
    const auto *hsk1_719 = buffer.data(hsk1 + 719);

    const auto *isi_469 = buffer.data(isi + 469);
    const auto *isi_474 = buffer.data(isi + 474);
    const auto *isi_475 = buffer.data(isi + 475);
    const auto *isi_476 = buffer.data(isi + 476);
    const auto *isi_478 = buffer.data(isi + 478);
    const auto *isi_479 = buffer.data(isi + 479);
    const auto *isi_481 = buffer.data(isi + 481);
    const auto *isi_482 = buffer.data(isi + 482);
    const auto *isi_485 = buffer.data(isi + 485);
    const auto *isi_486 = buffer.data(isi + 486);
    const auto *isi_490 = buffer.data(isi + 490);
    const auto *isi_497 = buffer.data(isi + 497);
    const auto *isi_498 = buffer.data(isi + 498);
    const auto *isi_499 = buffer.data(isi + 499);
    const auto *isi_500 = buffer.data(isi + 500);
    const auto *isi_501 = buffer.data(isi + 501);
    const auto *isi_502 = buffer.data(isi + 502);
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
    const auto *isi_526 = buffer.data(isi + 526);
    const auto *isi_527 = buffer.data(isi + 527);
    const auto *isi_528 = buffer.data(isi + 528);
    const auto *isi_529 = buffer.data(isi + 529);
    const auto *isi_530 = buffer.data(isi + 530);
    const auto *isi_531 = buffer.data(isi + 531);
    const auto *isi_532 = buffer.data(isi + 532);
    const auto *isi_534 = buffer.data(isi + 534);
    const auto *isi_535 = buffer.data(isi + 535);
    const auto *isi_537 = buffer.data(isi + 537);
    const auto *isi_538 = buffer.data(isi + 538);
    const auto *isi_541 = buffer.data(isi + 541);
    const auto *isi_542 = buffer.data(isi + 542);
    const auto *isi_546 = buffer.data(isi + 546);
    const auto *isi_553 = buffer.data(isi + 553);
    const auto *isi_554 = buffer.data(isi + 554);
    const auto *isi_555 = buffer.data(isi + 555);
    const auto *isi_556 = buffer.data(isi + 556);
    const auto *isi_557 = buffer.data(isi + 557);
    const auto *isi_558 = buffer.data(isi + 558);
    const auto *isi_559 = buffer.data(isi + 559);

#pragma omp simd aligned(t_602, t_603, t_604, t_605, pa_x, pc_x, pc_z, hsk0_604, hsi_301, \
                         hsi_474, hsi_475, hsk1_604, isi_469, isi_474, \
                         isi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_13 * hsi_474[k]
                   + f_3 * pc_x[k] * isi_474[k];

        t_603[k] = f_13 * hsi_475[k]
                   + f_3 * pc_x[k] * isi_475[k];

        t_604[k] = pa_x[k] * hsk0_604[k]
                   - f_12 * pc_x[k] * hsk1_604[k];

        t_605[k] = f_13 * hsi_301[k]
                   + f_3 * pc_z[k] * isi_469[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, pa_x, pc_x, hsk0_606, hsk0_607, hsk0_608, \
                         hsk0_609, hsk1_606, hsk1_607, hsk1_608, \
                         hsk1_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = pa_x[k] * hsk0_606[k]
                   - f_12 * pc_x[k] * hsk1_606[k];

        t_607[k] = pa_x[k] * hsk0_607[k]
                   - f_12 * pc_x[k] * hsk1_607[k];

        t_608[k] = pa_x[k] * hsk0_608[k]
                   - f_12 * pc_x[k] * hsk1_608[k];

        t_609[k] = pa_x[k] * hsk0_609[k]
                   - f_12 * pc_x[k] * hsk1_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_x, pc_x, pc_y, hsk0_611, hsk0_612, \
                         hsi_335, hsi_336, hsi_476, hsk1_611, hsk1_612, isi_475, \
                         isi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * hsi_335[k]
                   + f_3 * pc_y[k] * isi_475[k];

        t_611[k] = pa_x[k] * hsk0_611[k]
                   - f_12 * pc_x[k] * hsk1_611[k];

        t_612[k] = pa_x[k] * hsk0_612[k]
                   + f_20 * hsi_476[k]
                   - f_12 * pc_x[k] * hsk1_612[k];

        t_613[k] = f_15 * hsi_336[k]
                   + f_3 * pc_y[k] * isi_476[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, pa_x, pc_x, pc_y, pc_z, hsk0_615, hsi_308, \
                         hsi_338, hsi_479, hsk1_615, isi_476, isi_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_14 * hsi_308[k]
                   + f_3 * pc_z[k] * isi_476[k];

        t_615[k] = pa_x[k] * hsk0_615[k]
                   + f_17 * hsi_479[k]
                   - f_12 * pc_x[k] * hsk1_615[k];

        t_616[k] = f_15 * hsi_338[k]
                   + f_3 * pc_y[k] * isi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pa_x, pc_x, pc_z, hsk0_617, hsk0_618, hsi_311, \
                         hsi_481, hsi_482, hsk1_617, hsk1_618, \
                         isi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = pa_x[k] * hsk0_617[k]
                   + f_17 * hsi_481[k]
                   - f_12 * pc_x[k] * hsk1_617[k];

        t_618[k] = pa_x[k] * hsk0_618[k]
                   + f_16 * hsi_482[k]
                   - f_12 * pc_x[k] * hsk1_618[k];

        t_619[k] = f_14 * hsi_311[k]
                   + f_3 * pc_z[k] * isi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pa_x, pc_x, pc_y, hsk0_621, hsk0_622, hsi_341, \
                         hsi_485, hsi_486, hsk1_621, hsk1_622, \
                         isi_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * hsi_341[k]
                   + f_3 * pc_y[k] * isi_481[k];

        t_621[k] = pa_x[k] * hsk0_621[k]
                   + f_16 * hsi_485[k]
                   - f_12 * pc_x[k] * hsk1_621[k];

        t_622[k] = pa_x[k] * hsk0_622[k]
                   + f_15 * hsi_486[k]
                   - f_12 * pc_x[k] * hsk1_622[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pa_x, pc_x, pc_y, pc_z, hsk0_624, hsi_314, \
                         hsi_345, hsi_488, hsk1_624, isi_482, isi_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * hsi_314[k]
                   + f_3 * pc_z[k] * isi_482[k];

        t_624[k] = pa_x[k] * hsk0_624[k]
                   + f_15 * hsi_488[k]
                   - f_12 * pc_x[k] * hsk1_624[k];

        t_625[k] = f_15 * hsi_345[k]
                   + f_3 * pc_y[k] * isi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_x, pc_x, pc_z, hsk0_626, hsk0_627, hsi_318, \
                         hsi_490, hsi_491, hsk1_626, hsk1_627, \
                         isi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_x[k] * hsk0_626[k]
                   + f_15 * hsi_490[k]
                   - f_12 * pc_x[k] * hsk1_626[k];

        t_627[k] = pa_x[k] * hsk0_627[k]
                   + f_14 * hsi_491[k]
                   - f_12 * pc_x[k] * hsk1_627[k];

        t_628[k] = f_14 * hsi_318[k]
                   + f_3 * pc_z[k] * isi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pa_x, pc_x, pc_y, hsk0_629, hsk0_630, hsi_350, \
                         hsi_493, hsi_494, hsk1_629, hsk1_630, \
                         isi_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_x[k] * hsk0_629[k]
                   + f_14 * hsi_493[k]
                   - f_12 * pc_x[k] * hsk1_629[k];

        t_630[k] = pa_x[k] * hsk0_630[k]
                   + f_14 * hsi_494[k]
                   - f_12 * pc_x[k] * hsk1_630[k];

        t_631[k] = f_15 * hsi_350[k]
                   + f_3 * pc_y[k] * isi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pa_x, pc_x, hsk0_632, hsi_496, hsi_497, \
                         hsi_498, hsi_499, hsk1_632, isi_497, isi_498, \
                         isi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = pa_x[k] * hsk0_632[k]
                   + f_14 * hsi_496[k]
                   - f_12 * pc_x[k] * hsk1_632[k];

        t_633[k] = f_13 * hsi_497[k]
                   + f_3 * pc_x[k] * isi_497[k];

        t_634[k] = f_13 * hsi_498[k]
                   + f_3 * pc_x[k] * isi_498[k];

        t_635[k] = f_13 * hsi_499[k]
                   + f_3 * pc_x[k] * isi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, hsi_500, hsi_501, hsi_502, hsi_503, \
                         isi_500, isi_501, isi_502, isi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_13 * hsi_500[k]
                   + f_3 * pc_x[k] * isi_500[k];

        t_637[k] = f_13 * hsi_501[k]
                   + f_3 * pc_x[k] * isi_501[k];

        t_638[k] = f_13 * hsi_502[k]
                   + f_3 * pc_x[k] * isi_502[k];

        t_639[k] = f_13 * hsi_503[k]
                   + f_3 * pc_x[k] * isi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, pa_x, pc_x, pc_z, hsk0_640, hsk0_642, \
                         hsk0_643, hsi_329, hsk1_640, hsk1_642, hsk1_643, \
                         isi_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = pa_x[k] * hsk0_640[k]
                   - f_12 * pc_x[k] * hsk1_640[k];

        t_641[k] = f_14 * hsi_329[k]
                   + f_3 * pc_z[k] * isi_497[k];

        t_642[k] = pa_x[k] * hsk0_642[k]
                   - f_12 * pc_x[k] * hsk1_642[k];

        t_643[k] = pa_x[k] * hsk0_643[k]
                   - f_12 * pc_x[k] * hsk1_643[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, pa_x, pc_x, pc_y, hsk0_644, hsk0_645, \
                         hsk0_647, hsi_363, hsk1_644, hsk1_645, hsk1_647, \
                         isi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = pa_x[k] * hsk0_644[k]
                   - f_12 * pc_x[k] * hsk1_644[k];

        t_645[k] = pa_x[k] * hsk0_645[k]
                   - f_12 * pc_x[k] * hsk1_645[k];

        t_646[k] = f_15 * hsi_363[k]
                   + f_3 * pc_y[k] * isi_503[k];

        t_647[k] = pa_x[k] * hsk0_647[k]
                   - f_12 * pc_x[k] * hsk1_647[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pa_x, pc_x, pc_y, pc_z, hsk0_648, hsi_336, \
                         hsi_364, hsi_504, hsk1_648, isi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = pa_x[k] * hsk0_648[k]
                   + f_20 * hsi_504[k]
                   - f_12 * pc_x[k] * hsk1_648[k];

        t_649[k] = f_14 * hsi_364[k]
                   + f_3 * pc_y[k] * isi_504[k];

        t_650[k] = f_15 * hsi_336[k]
                   + f_3 * pc_z[k] * isi_504[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pa_x, pc_x, pc_y, hsk0_651, hsk0_653, hsi_366, \
                         hsi_507, hsi_509, hsk1_651, hsk1_653, \
                         isi_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = pa_x[k] * hsk0_651[k]
                   + f_17 * hsi_507[k]
                   - f_12 * pc_x[k] * hsk1_651[k];

        t_652[k] = f_14 * hsi_366[k]
                   + f_3 * pc_y[k] * isi_506[k];

        t_653[k] = pa_x[k] * hsk0_653[k]
                   + f_17 * hsi_509[k]
                   - f_12 * pc_x[k] * hsk1_653[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pa_x, pc_x, pc_y, pc_z, hsk0_654, hsi_339, \
                         hsi_369, hsi_510, hsk1_654, isi_507, isi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = pa_x[k] * hsk0_654[k]
                   + f_16 * hsi_510[k]
                   - f_12 * pc_x[k] * hsk1_654[k];

        t_655[k] = f_15 * hsi_339[k]
                   + f_3 * pc_z[k] * isi_507[k];

        t_656[k] = f_14 * hsi_369[k]
                   + f_3 * pc_y[k] * isi_509[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_x, pc_x, pc_z, hsk0_657, hsk0_658, hsi_342, \
                         hsi_513, hsi_514, hsk1_657, hsk1_658, \
                         isi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pa_x[k] * hsk0_657[k]
                   + f_16 * hsi_513[k]
                   - f_12 * pc_x[k] * hsk1_657[k];

        t_658[k] = pa_x[k] * hsk0_658[k]
                   + f_15 * hsi_514[k]
                   - f_12 * pc_x[k] * hsk1_658[k];

        t_659[k] = f_15 * hsi_342[k]
                   + f_3 * pc_z[k] * isi_510[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, pa_x, pc_x, pc_y, hsk0_660, hsk0_662, hsi_373, \
                         hsi_516, hsi_518, hsk1_660, hsk1_662, \
                         isi_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = pa_x[k] * hsk0_660[k]
                   + f_15 * hsi_516[k]
                   - f_12 * pc_x[k] * hsk1_660[k];

        t_661[k] = f_14 * hsi_373[k]
                   + f_3 * pc_y[k] * isi_513[k];

        t_662[k] = pa_x[k] * hsk0_662[k]
                   + f_15 * hsi_518[k]
                   - f_12 * pc_x[k] * hsk1_662[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, pa_x, pc_x, pc_z, hsk0_663, hsk0_665, hsi_346, \
                         hsi_519, hsi_521, hsk1_663, hsk1_665, \
                         isi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = pa_x[k] * hsk0_663[k]
                   + f_14 * hsi_519[k]
                   - f_12 * pc_x[k] * hsk1_663[k];

        t_664[k] = f_15 * hsi_346[k]
                   + f_3 * pc_z[k] * isi_514[k];

        t_665[k] = pa_x[k] * hsk0_665[k]
                   + f_14 * hsi_521[k]
                   - f_12 * pc_x[k] * hsk1_665[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, pa_x, pc_x, pc_y, hsk0_666, hsk0_668, hsi_378, \
                         hsi_522, hsi_524, hsk1_666, hsk1_668, \
                         isi_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = pa_x[k] * hsk0_666[k]
                   + f_14 * hsi_522[k]
                   - f_12 * pc_x[k] * hsk1_666[k];

        t_667[k] = f_14 * hsi_378[k]
                   + f_3 * pc_y[k] * isi_518[k];

        t_668[k] = pa_x[k] * hsk0_668[k]
                   + f_14 * hsi_524[k]
                   - f_12 * pc_x[k] * hsk1_668[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, pc_x, hsi_525, hsi_526, hsi_527, \
                         hsi_528, hsi_529, isi_525, isi_526, isi_527, isi_528, \
                         isi_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_13 * hsi_525[k]
                   + f_3 * pc_x[k] * isi_525[k];

        t_670[k] = f_13 * hsi_526[k]
                   + f_3 * pc_x[k] * isi_526[k];

        t_671[k] = f_13 * hsi_527[k]
                   + f_3 * pc_x[k] * isi_527[k];

        t_672[k] = f_13 * hsi_528[k]
                   + f_3 * pc_x[k] * isi_528[k];

        t_673[k] = f_13 * hsi_529[k]
                   + f_3 * pc_x[k] * isi_529[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pa_x, pc_x, pc_z, hsk0_676, hsi_357, \
                         hsi_530, hsi_531, hsk1_676, isi_525, isi_530, \
                         isi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_13 * hsi_530[k]
                   + f_3 * pc_x[k] * isi_530[k];

        t_675[k] = f_13 * hsi_531[k]
                   + f_3 * pc_x[k] * isi_531[k];

        t_676[k] = pa_x[k] * hsk0_676[k]
                   - f_12 * pc_x[k] * hsk1_676[k];

        t_677[k] = f_15 * hsi_357[k]
                   + f_3 * pc_z[k] * isi_525[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, t_681, pa_x, pc_x, hsk0_678, hsk0_679, hsk0_680, \
                         hsk0_681, hsk1_678, hsk1_679, hsk1_680, \
                         hsk1_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = pa_x[k] * hsk0_678[k]
                   - f_12 * pc_x[k] * hsk1_678[k];

        t_679[k] = pa_x[k] * hsk0_679[k]
                   - f_12 * pc_x[k] * hsk1_679[k];

        t_680[k] = pa_x[k] * hsk0_680[k]
                   - f_12 * pc_x[k] * hsk1_680[k];

        t_681[k] = pa_x[k] * hsk0_681[k]
                   - f_12 * pc_x[k] * hsk1_681[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_x, pa_y, pc_x, pc_y, hsk0_504, \
                         hsk0_683, hsi_391, hsi_392, hsk1_504, hsk1_683, isi_531, \
                         isi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * hsi_391[k]
                   + f_3 * pc_y[k] * isi_531[k];

        t_683[k] = pa_x[k] * hsk0_683[k]
                   - f_12 * pc_x[k] * hsk1_683[k];

        t_684[k] = pa_y[k] * hsk0_504[k]
                   - f_12 * pc_y[k] * hsk1_504[k];

        t_685[k] = f_13 * hsi_392[k]
                   + f_3 * pc_y[k] * isi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, pa_x, pc_x, pc_y, pc_z, hsk0_687, hsi_364, \
                         hsi_394, hsi_535, hsk1_687, isi_532, isi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * hsi_364[k]
                   + f_3 * pc_z[k] * isi_532[k];

        t_687[k] = pa_x[k] * hsk0_687[k]
                   + f_17 * hsi_535[k]
                   - f_12 * pc_x[k] * hsk1_687[k];

        t_688[k] = f_13 * hsi_394[k]
                   + f_3 * pc_y[k] * isi_534[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pa_x, pa_y, pc_x, pc_y, pc_z, hsk0_509, \
                         hsk0_690, hsi_367, hsi_538, hsk1_509, hsk1_690, \
                         isi_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = pa_y[k] * hsk0_509[k]
                   - f_12 * pc_y[k] * hsk1_509[k];

        t_690[k] = pa_x[k] * hsk0_690[k]
                   + f_16 * hsi_538[k]
                   - f_12 * pc_x[k] * hsk1_690[k];

        t_691[k] = f_16 * hsi_367[k]
                   + f_3 * pc_z[k] * isi_535[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pa_x, pa_y, pc_x, pc_y, hsk0_513, hsk0_694, \
                         hsi_397, hsi_542, hsk1_513, hsk1_694, \
                         isi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_13 * hsi_397[k]
                   + f_3 * pc_y[k] * isi_537[k];

        t_693[k] = pa_y[k] * hsk0_513[k]
                   - f_12 * pc_y[k] * hsk1_513[k];

        t_694[k] = pa_x[k] * hsk0_694[k]
                   + f_15 * hsi_542[k]
                   - f_12 * pc_x[k] * hsk1_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pa_x, pc_x, pc_y, pc_z, hsk0_696, hsi_370, \
                         hsi_401, hsi_544, hsk1_696, isi_538, isi_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_16 * hsi_370[k]
                   + f_3 * pc_z[k] * isi_538[k];

        t_696[k] = pa_x[k] * hsk0_696[k]
                   + f_15 * hsi_544[k]
                   - f_12 * pc_x[k] * hsk1_696[k];

        t_697[k] = f_13 * hsi_401[k]
                   + f_3 * pc_y[k] * isi_541[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pa_x, pa_y, pc_x, pc_y, pc_z, hsk0_518, \
                         hsk0_699, hsi_374, hsi_547, hsk1_518, hsk1_699, \
                         isi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = pa_y[k] * hsk0_518[k]
                   - f_12 * pc_y[k] * hsk1_518[k];

        t_699[k] = pa_x[k] * hsk0_699[k]
                   + f_14 * hsi_547[k]
                   - f_12 * pc_x[k] * hsk1_699[k];

        t_700[k] = f_16 * hsi_374[k]
                   + f_3 * pc_z[k] * isi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, pa_x, pc_x, pc_y, hsk0_701, hsk0_702, hsi_406, \
                         hsi_549, hsi_550, hsk1_701, hsk1_702, \
                         isi_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pa_x[k] * hsk0_701[k]
                   + f_14 * hsi_549[k]
                   - f_12 * pc_x[k] * hsk1_701[k];

        t_702[k] = pa_x[k] * hsk0_702[k]
                   + f_14 * hsi_550[k]
                   - f_12 * pc_x[k] * hsk1_702[k];

        t_703[k] = f_13 * hsi_406[k]
                   + f_3 * pc_y[k] * isi_546[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, pa_y, pc_x, pc_y, hsk0_524, hsi_553, \
                         hsi_554, hsi_555, hsk1_524, isi_553, isi_554, \
                         isi_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = pa_y[k] * hsk0_524[k]
                   - f_12 * pc_y[k] * hsk1_524[k];

        t_705[k] = f_13 * hsi_553[k]
                   + f_3 * pc_x[k] * isi_553[k];

        t_706[k] = f_13 * hsi_554[k]
                   + f_3 * pc_x[k] * isi_554[k];

        t_707[k] = f_13 * hsi_555[k]
                   + f_3 * pc_x[k] * isi_555[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pc_x, hsi_556, hsi_557, hsi_558, hsi_559, \
                         isi_556, isi_557, isi_558, isi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_13 * hsi_556[k]
                   + f_3 * pc_x[k] * isi_556[k];

        t_709[k] = f_13 * hsi_557[k]
                   + f_3 * pc_x[k] * isi_557[k];

        t_710[k] = f_13 * hsi_558[k]
                   + f_3 * pc_x[k] * isi_558[k];

        t_711[k] = f_13 * hsi_559[k]
                   + f_3 * pc_x[k] * isi_559[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_x, pc_x, pc_z, hsk0_712, hsk0_714, \
                         hsk0_715, hsi_385, hsk1_712, hsk1_714, hsk1_715, \
                         isi_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = pa_x[k] * hsk0_712[k]
                   - f_12 * pc_x[k] * hsk1_712[k];

        t_713[k] = f_16 * hsi_385[k]
                   + f_3 * pc_z[k] * isi_553[k];

        t_714[k] = pa_x[k] * hsk0_714[k]
                   - f_12 * pc_x[k] * hsk1_714[k];

        t_715[k] = pa_x[k] * hsk0_715[k]
                   - f_12 * pc_x[k] * hsk1_715[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_x, pc_x, pc_y, hsk0_716, hsk0_717, \
                         hsk0_719, hsi_419, hsk1_716, hsk1_717, hsk1_719, \
                         isi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = pa_x[k] * hsk0_716[k]
                   - f_12 * pc_x[k] * hsk1_716[k];

        t_717[k] = pa_x[k] * hsk0_717[k]
                   - f_12 * pc_x[k] * hsk1_717[k];

        t_718[k] = f_13 * hsi_419[k]
                   + f_3 * pc_y[k] * isi_559[k];

        t_719[k] = pa_x[k] * hsk0_719[k]
                   - f_12 * pc_x[k] * hsk1_719[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);
    const auto f_20 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsk0_540 = buffer.data(hsk0 + 540);
    const auto *hsk0_541 = buffer.data(hsk0 + 541);
    const auto *hsk0_543 = buffer.data(hsk0 + 543);
    const auto *hsk0_546 = buffer.data(hsk0 + 546);
    const auto *hsk0_550 = buffer.data(hsk0 + 550);
    const auto *hsk0_555 = buffer.data(hsk0 + 555);
    const auto *hsk0_568 = buffer.data(hsk0 + 568);
    const auto *hsk0_570 = buffer.data(hsk0 + 570);
    const auto *hsk0_571 = buffer.data(hsk0 + 571);
    const auto *hsk0_572 = buffer.data(hsk0 + 572);
    const auto *hsk0_573 = buffer.data(hsk0 + 573);
    const auto *hsk0_720 = buffer.data(hsk0 + 720);
    const auto *hsk0_725 = buffer.data(hsk0 + 725);
    const auto *hsk0_729 = buffer.data(hsk0 + 729);
    const auto *hsk0_734 = buffer.data(hsk0 + 734);
    const auto *hsk0_740 = buffer.data(hsk0 + 740);
    const auto *hsk0_748 = buffer.data(hsk0 + 748);
    const auto *hsk0_749 = buffer.data(hsk0 + 749);
    const auto *hsk0_750 = buffer.data(hsk0 + 750);
    const auto *hsk0_751 = buffer.data(hsk0 + 751);
    const auto *hsk0_752 = buffer.data(hsk0 + 752);
    const auto *hsk0_753 = buffer.data(hsk0 + 753);
    const auto *hsk0_755 = buffer.data(hsk0 + 755);

    const auto *hsi_392 = buffer.data(hsi + 392);
    const auto *hsi_441 = buffer.data(hsi + 441);
    const auto *hsi_442 = buffer.data(hsi + 442);
    const auto *hsi_443 = buffer.data(hsi + 443);
    const auto *hsi_444 = buffer.data(hsi + 444);
    const auto *hsi_445 = buffer.data(hsi + 445);
    const auto *hsi_447 = buffer.data(hsi + 447);
    const auto *hsi_475 = buffer.data(hsi + 475);
    const auto *hsi_560 = buffer.data(hsi + 560);
    const auto *hsi_565 = buffer.data(hsi + 565);
    const auto *hsi_569 = buffer.data(hsi + 569);
    const auto *hsi_574 = buffer.data(hsi + 574);
    const auto *hsi_580 = buffer.data(hsi + 580);
    const auto *hsi_581 = buffer.data(hsi + 581);
    const auto *hsi_582 = buffer.data(hsi + 582);
    const auto *hsi_583 = buffer.data(hsi + 583);
    const auto *hsi_584 = buffer.data(hsi + 584);
    const auto *hsi_585 = buffer.data(hsi + 585);
    const auto *hsi_587 = buffer.data(hsi + 587);

    const auto *hsk1_540 = buffer.data(hsk1 + 540);
    const auto *hsk1_541 = buffer.data(hsk1 + 541);
    const auto *hsk1_543 = buffer.data(hsk1 + 543);
    const auto *hsk1_546 = buffer.data(hsk1 + 546);
    const auto *hsk1_550 = buffer.data(hsk1 + 550);
    const auto *hsk1_555 = buffer.data(hsk1 + 555);
    const auto *hsk1_568 = buffer.data(hsk1 + 568);
    const auto *hsk1_570 = buffer.data(hsk1 + 570);
    const auto *hsk1_571 = buffer.data(hsk1 + 571);
    const auto *hsk1_572 = buffer.data(hsk1 + 572);
    const auto *hsk1_573 = buffer.data(hsk1 + 573);
    const auto *hsk1_720 = buffer.data(hsk1 + 720);
    const auto *hsk1_725 = buffer.data(hsk1 + 725);
    const auto *hsk1_729 = buffer.data(hsk1 + 729);
    const auto *hsk1_734 = buffer.data(hsk1 + 734);
    const auto *hsk1_740 = buffer.data(hsk1 + 740);
    const auto *hsk1_748 = buffer.data(hsk1 + 748);
    const auto *hsk1_749 = buffer.data(hsk1 + 749);
    const auto *hsk1_750 = buffer.data(hsk1 + 750);
    const auto *hsk1_751 = buffer.data(hsk1 + 751);
    const auto *hsk1_752 = buffer.data(hsk1 + 752);
    const auto *hsk1_753 = buffer.data(hsk1 + 753);
    const auto *hsk1_755 = buffer.data(hsk1 + 755);

    const auto *ish0_420 = buffer.data(ish0 + 420);
    const auto *ish0_421 = buffer.data(ish0 + 421);
    const auto *ish0_422 = buffer.data(ish0 + 422);
    const auto *ish0_423 = buffer.data(ish0 + 423);
    const auto *ish0_424 = buffer.data(ish0 + 424);
    const auto *ish0_425 = buffer.data(ish0 + 425);
    const auto *ish0_426 = buffer.data(ish0 + 426);
    const auto *ish0_427 = buffer.data(ish0 + 427);
    const auto *ish0_428 = buffer.data(ish0 + 428);
    const auto *ish0_429 = buffer.data(ish0 + 429);
    const auto *ish0_441 = buffer.data(ish0 + 441);
    const auto *ish0_442 = buffer.data(ish0 + 442);
    const auto *ish0_444 = buffer.data(ish0 + 444);
    const auto *ish0_446 = buffer.data(ish0 + 446);
    const auto *ish0_447 = buffer.data(ish0 + 447);
    const auto *ish0_449 = buffer.data(ish0 + 449);
    const auto *ish0_450 = buffer.data(ish0 + 450);
    const auto *ish0_451 = buffer.data(ish0 + 451);
    const auto *ish0_453 = buffer.data(ish0 + 453);
    const auto *ish0_454 = buffer.data(ish0 + 454);
    const auto *ish0_455 = buffer.data(ish0 + 455);
    const auto *ish0_456 = buffer.data(ish0 + 456);
    const auto *ish0_457 = buffer.data(ish0 + 457);
    const auto *ish0_458 = buffer.data(ish0 + 458);
    const auto *ish0_459 = buffer.data(ish0 + 459);
    const auto *ish0_460 = buffer.data(ish0 + 460);
    const auto *ish0_461 = buffer.data(ish0 + 461);
    const auto *ish0_464 = buffer.data(ish0 + 464);
    const auto *ish0_466 = buffer.data(ish0 + 466);
    const auto *ish0_467 = buffer.data(ish0 + 467);
    const auto *ish0_469 = buffer.data(ish0 + 469);
    const auto *ish0_470 = buffer.data(ish0 + 470);
    const auto *ish0_471 = buffer.data(ish0 + 471);
    const auto *ish0_473 = buffer.data(ish0 + 473);
    const auto *ish0_474 = buffer.data(ish0 + 474);
    const auto *ish0_475 = buffer.data(ish0 + 475);
    const auto *ish0_476 = buffer.data(ish0 + 476);
    const auto *ish0_478 = buffer.data(ish0 + 478);
    const auto *ish0_479 = buffer.data(ish0 + 479);
    const auto *ish0_480 = buffer.data(ish0 + 480);
    const auto *ish0_481 = buffer.data(ish0 + 481);
    const auto *ish0_482 = buffer.data(ish0 + 482);
    const auto *ish0_483 = buffer.data(ish0 + 483);
    const auto *ish0_484 = buffer.data(ish0 + 484);
    const auto *ish0_485 = buffer.data(ish0 + 485);
    const auto *ish0_486 = buffer.data(ish0 + 486);
    const auto *ish0_487 = buffer.data(ish0 + 487);
    const auto *ish0_488 = buffer.data(ish0 + 488);
    const auto *ish0_489 = buffer.data(ish0 + 489);
    const auto *ish0_490 = buffer.data(ish0 + 490);
    const auto *ish0_491 = buffer.data(ish0 + 491);
    const auto *ish0_492 = buffer.data(ish0 + 492);
    const auto *ish0_493 = buffer.data(ish0 + 493);
    const auto *ish0_494 = buffer.data(ish0 + 494);
    const auto *ish0_495 = buffer.data(ish0 + 495);
    const auto *ish0_496 = buffer.data(ish0 + 496);
    const auto *ish0_497 = buffer.data(ish0 + 497);

    const auto *ish1_420 = buffer.data(ish1 + 420);
    const auto *ish1_421 = buffer.data(ish1 + 421);
    const auto *ish1_422 = buffer.data(ish1 + 422);
    const auto *ish1_423 = buffer.data(ish1 + 423);
    const auto *ish1_424 = buffer.data(ish1 + 424);
    const auto *ish1_425 = buffer.data(ish1 + 425);
    const auto *ish1_426 = buffer.data(ish1 + 426);
    const auto *ish1_427 = buffer.data(ish1 + 427);
    const auto *ish1_428 = buffer.data(ish1 + 428);
    const auto *ish1_429 = buffer.data(ish1 + 429);
    const auto *ish1_441 = buffer.data(ish1 + 441);
    const auto *ish1_442 = buffer.data(ish1 + 442);
    const auto *ish1_444 = buffer.data(ish1 + 444);
    const auto *ish1_446 = buffer.data(ish1 + 446);
    const auto *ish1_447 = buffer.data(ish1 + 447);
    const auto *ish1_449 = buffer.data(ish1 + 449);
    const auto *ish1_450 = buffer.data(ish1 + 450);
    const auto *ish1_451 = buffer.data(ish1 + 451);
    const auto *ish1_453 = buffer.data(ish1 + 453);
    const auto *ish1_454 = buffer.data(ish1 + 454);
    const auto *ish1_455 = buffer.data(ish1 + 455);
    const auto *ish1_456 = buffer.data(ish1 + 456);
    const auto *ish1_457 = buffer.data(ish1 + 457);
    const auto *ish1_458 = buffer.data(ish1 + 458);
    const auto *ish1_459 = buffer.data(ish1 + 459);
    const auto *ish1_460 = buffer.data(ish1 + 460);
    const auto *ish1_461 = buffer.data(ish1 + 461);
    const auto *ish1_464 = buffer.data(ish1 + 464);
    const auto *ish1_466 = buffer.data(ish1 + 466);
    const auto *ish1_467 = buffer.data(ish1 + 467);
    const auto *ish1_469 = buffer.data(ish1 + 469);
    const auto *ish1_470 = buffer.data(ish1 + 470);
    const auto *ish1_471 = buffer.data(ish1 + 471);
    const auto *ish1_473 = buffer.data(ish1 + 473);
    const auto *ish1_474 = buffer.data(ish1 + 474);
    const auto *ish1_475 = buffer.data(ish1 + 475);
    const auto *ish1_476 = buffer.data(ish1 + 476);
    const auto *ish1_478 = buffer.data(ish1 + 478);
    const auto *ish1_479 = buffer.data(ish1 + 479);
    const auto *ish1_480 = buffer.data(ish1 + 480);
    const auto *ish1_481 = buffer.data(ish1 + 481);
    const auto *ish1_482 = buffer.data(ish1 + 482);
    const auto *ish1_483 = buffer.data(ish1 + 483);
    const auto *ish1_484 = buffer.data(ish1 + 484);
    const auto *ish1_485 = buffer.data(ish1 + 485);
    const auto *ish1_486 = buffer.data(ish1 + 486);
    const auto *ish1_487 = buffer.data(ish1 + 487);
    const auto *ish1_488 = buffer.data(ish1 + 488);
    const auto *ish1_489 = buffer.data(ish1 + 489);
    const auto *ish1_490 = buffer.data(ish1 + 490);
    const auto *ish1_491 = buffer.data(ish1 + 491);
    const auto *ish1_492 = buffer.data(ish1 + 492);
    const auto *ish1_493 = buffer.data(ish1 + 493);
    const auto *ish1_494 = buffer.data(ish1 + 494);
    const auto *ish1_495 = buffer.data(ish1 + 495);
    const auto *ish1_496 = buffer.data(ish1 + 496);
    const auto *ish1_497 = buffer.data(ish1 + 497);

    const auto *isi_560 = buffer.data(isi + 560);
    const auto *isi_561 = buffer.data(isi + 561);
    const auto *isi_562 = buffer.data(isi + 562);
    const auto *isi_563 = buffer.data(isi + 563);
    const auto *isi_564 = buffer.data(isi + 564);
    const auto *isi_565 = buffer.data(isi + 565);
    const auto *isi_566 = buffer.data(isi + 566);
    const auto *isi_567 = buffer.data(isi + 567);
    const auto *isi_568 = buffer.data(isi + 568);
    const auto *isi_569 = buffer.data(isi + 569);
    const auto *isi_570 = buffer.data(isi + 570);
    const auto *isi_571 = buffer.data(isi + 571);
    const auto *isi_572 = buffer.data(isi + 572);
    const auto *isi_573 = buffer.data(isi + 573);
    const auto *isi_574 = buffer.data(isi + 574);
    const auto *isi_580 = buffer.data(isi + 580);
    const auto *isi_581 = buffer.data(isi + 581);
    const auto *isi_582 = buffer.data(isi + 582);
    const auto *isi_583 = buffer.data(isi + 583);
    const auto *isi_584 = buffer.data(isi + 584);
    const auto *isi_585 = buffer.data(isi + 585);
    const auto *isi_587 = buffer.data(isi + 587);
    const auto *isi_588 = buffer.data(isi + 588);
    const auto *isi_589 = buffer.data(isi + 589);
    const auto *isi_591 = buffer.data(isi + 591);
    const auto *isi_593 = buffer.data(isi + 593);
    const auto *isi_594 = buffer.data(isi + 594);
    const auto *isi_596 = buffer.data(isi + 596);
    const auto *isi_597 = buffer.data(isi + 597);
    const auto *isi_598 = buffer.data(isi + 598);
    const auto *isi_600 = buffer.data(isi + 600);
    const auto *isi_601 = buffer.data(isi + 601);
    const auto *isi_602 = buffer.data(isi + 602);
    const auto *isi_603 = buffer.data(isi + 603);
    const auto *isi_605 = buffer.data(isi + 605);
    const auto *isi_606 = buffer.data(isi + 606);
    const auto *isi_607 = buffer.data(isi + 607);
    const auto *isi_608 = buffer.data(isi + 608);
    const auto *isi_609 = buffer.data(isi + 609);
    const auto *isi_610 = buffer.data(isi + 610);
    const auto *isi_611 = buffer.data(isi + 611);
    const auto *isi_612 = buffer.data(isi + 612);
    const auto *isi_613 = buffer.data(isi + 613);
    const auto *isi_614 = buffer.data(isi + 614);
    const auto *isi_615 = buffer.data(isi + 615);
    const auto *isi_618 = buffer.data(isi + 618);
    const auto *isi_620 = buffer.data(isi + 620);
    const auto *isi_621 = buffer.data(isi + 621);
    const auto *isi_623 = buffer.data(isi + 623);
    const auto *isi_624 = buffer.data(isi + 624);
    const auto *isi_625 = buffer.data(isi + 625);
    const auto *isi_627 = buffer.data(isi + 627);
    const auto *isi_628 = buffer.data(isi + 628);
    const auto *isi_629 = buffer.data(isi + 629);
    const auto *isi_630 = buffer.data(isi + 630);
    const auto *isi_632 = buffer.data(isi + 632);
    const auto *isi_633 = buffer.data(isi + 633);
    const auto *isi_634 = buffer.data(isi + 634);
    const auto *isi_635 = buffer.data(isi + 635);
    const auto *isi_636 = buffer.data(isi + 636);
    const auto *isi_637 = buffer.data(isi + 637);
    const auto *isi_638 = buffer.data(isi + 638);
    const auto *isi_639 = buffer.data(isi + 639);
    const auto *isi_640 = buffer.data(isi + 640);
    const auto *isi_641 = buffer.data(isi + 641);
    const auto *isi_642 = buffer.data(isi + 642);
    const auto *isi_643 = buffer.data(isi + 643);
    const auto *isi_644 = buffer.data(isi + 644);
    const auto *isi_645 = buffer.data(isi + 645);
    const auto *isi_646 = buffer.data(isi + 646);
    const auto *isi_647 = buffer.data(isi + 647);
    const auto *isi_648 = buffer.data(isi + 648);
    const auto *isi_649 = buffer.data(isi + 649);
    const auto *isi_650 = buffer.data(isi + 650);
    const auto *isi_651 = buffer.data(isi + 651);
    const auto *isi_652 = buffer.data(isi + 652);
    const auto *isi_653 = buffer.data(isi + 653);
    const auto *isi_654 = buffer.data(isi + 654);
    const auto *isi_655 = buffer.data(isi + 655);
    const auto *isi_656 = buffer.data(isi + 656);
    const auto *isi_657 = buffer.data(isi + 657);
    const auto *isi_658 = buffer.data(isi + 658);

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_x, pc_x, pc_y, pc_z, hsk0_720, \
                         hsi_392, hsi_560, hsk1_720, ish0_420, ish1_420, isi_560, \
                         isi_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_x[k] * hsk0_720[k]
                   + f_20 * hsi_560[k]
                   - f_12 * pc_x[k] * hsk1_720[k];

        t_721[k] = f_3 * pc_y[k] * isi_560[k];

        t_722[k] = f_17 * hsi_392[k]
                   + f_3 * pc_z[k] * isi_560[k];

        t_723[k] = f_4 * ish0_420[k]
                   - f_5 * ish1_420[k]
                   + f_3 * pc_y[k] * isi_561[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pa_x, pc_x, pc_y, hsk0_725, hsi_565, hsk1_725, \
                         ish0_421, ish1_421, isi_562, isi_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_3 * pc_y[k] * isi_562[k];

        t_725[k] = pa_x[k] * hsk0_725[k]
                   + f_17 * hsi_565[k]
                   - f_12 * pc_x[k] * hsk1_725[k];

        t_726[k] = f_6 * ish0_421[k]
                   - f_7 * ish1_421[k]
                   + f_3 * pc_y[k] * isi_563[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pa_x, pc_x, pc_y, hsk0_729, hsi_569, hsk1_729, \
                         ish0_422, ish1_422, isi_564, isi_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_4 * ish0_422[k]
                   - f_5 * ish1_422[k]
                   + f_3 * pc_y[k] * isi_564[k];

        t_728[k] = f_3 * pc_y[k] * isi_565[k];

        t_729[k] = pa_x[k] * hsk0_729[k]
                   + f_16 * hsi_569[k]
                   - f_12 * pc_x[k] * hsk1_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pc_y, ish0_423, ish0_424, ish0_425, \
                         ish1_423, ish1_424, ish1_425, isi_566, isi_567, isi_568, \
                         isi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_8 * ish0_423[k]
                   - f_9 * ish1_423[k]
                   + f_3 * pc_y[k] * isi_566[k];

        t_731[k] = f_6 * ish0_424[k]
                   - f_7 * ish1_424[k]
                   + f_3 * pc_y[k] * isi_567[k];

        t_732[k] = f_4 * ish0_425[k]
                   - f_5 * ish1_425[k]
                   + f_3 * pc_y[k] * isi_568[k];

        t_733[k] = f_3 * pc_y[k] * isi_569[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pa_x, pc_x, pc_y, hsk0_734, hsi_574, hsk1_734, \
                         ish0_426, ish0_427, ish1_426, ish1_427, isi_570, \
                         isi_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_x[k] * hsk0_734[k]
                   + f_15 * hsi_574[k]
                   - f_12 * pc_x[k] * hsk1_734[k];

        t_735[k] = f_10 * ish0_426[k]
                   - f_11 * ish1_426[k]
                   + f_3 * pc_y[k] * isi_570[k];

        t_736[k] = f_8 * ish0_427[k]
                   - f_9 * ish1_427[k]
                   + f_3 * pc_y[k] * isi_571[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pc_y, ish0_428, ish0_429, ish1_428, ish1_429, \
                         isi_572, isi_573, isi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_6 * ish0_428[k]
                   - f_7 * ish1_428[k]
                   + f_3 * pc_y[k] * isi_572[k];

        t_738[k] = f_4 * ish0_429[k]
                   - f_5 * ish1_429[k]
                   + f_3 * pc_y[k] * isi_573[k];

        t_739[k] = f_3 * pc_y[k] * isi_574[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pa_x, pc_x, hsk0_740, hsi_580, hsi_581, \
                         hsi_582, hsi_583, hsk1_740, isi_581, isi_582, \
                         isi_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = pa_x[k] * hsk0_740[k]
                   + f_14 * hsi_580[k]
                   - f_12 * pc_x[k] * hsk1_740[k];

        t_741[k] = f_13 * hsi_581[k]
                   + f_3 * pc_x[k] * isi_581[k];

        t_742[k] = f_13 * hsi_582[k]
                   + f_3 * pc_x[k] * isi_582[k];

        t_743[k] = f_13 * hsi_583[k]
                   + f_3 * pc_x[k] * isi_583[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pc_x, pc_y, hsi_584, hsi_585, hsi_587, \
                         isi_580, isi_584, isi_585, isi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_13 * hsi_584[k]
                   + f_3 * pc_x[k] * isi_584[k];

        t_745[k] = f_13 * hsi_585[k]
                   + f_3 * pc_x[k] * isi_585[k];

        t_746[k] = f_3 * pc_y[k] * isi_580[k];

        t_747[k] = f_13 * hsi_587[k]
                   + f_3 * pc_x[k] * isi_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pa_x, pc_x, hsk0_748, hsk0_749, hsk0_750, \
                         hsk0_751, hsk1_748, hsk1_749, hsk1_750, \
                         hsk1_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = pa_x[k] * hsk0_748[k]
                   - f_12 * pc_x[k] * hsk1_748[k];

        t_749[k] = pa_x[k] * hsk0_749[k]
                   - f_12 * pc_x[k] * hsk1_749[k];

        t_750[k] = pa_x[k] * hsk0_750[k]
                   - f_12 * pc_x[k] * hsk1_750[k];

        t_751[k] = pa_x[k] * hsk0_751[k]
                   - f_12 * pc_x[k] * hsk1_751[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pa_x, pc_x, pc_y, hsk0_752, hsk0_753, \
                         hsk0_755, hsk1_752, hsk1_753, hsk1_755, \
                         isi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = pa_x[k] * hsk0_752[k]
                   - f_12 * pc_x[k] * hsk1_752[k];

        t_753[k] = pa_x[k] * hsk0_753[k]
                   - f_12 * pc_x[k] * hsk1_753[k];

        t_754[k] = f_3 * pc_y[k] * isi_587[k];

        t_755[k] = pa_x[k] * hsk0_755[k]
                   - f_12 * pc_x[k] * hsk1_755[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, pc_x, pc_z, ish0_441, ish0_442, \
                         ish0_444, ish1_441, ish1_442, ish1_444, isi_588, isi_589, \
                         isi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_1 * ish0_441[k]
                   - f_2 * ish1_441[k]
                   + f_3 * pc_x[k] * isi_588[k];

        t_757[k] = f_18 * ish0_442[k]
                   - f_19 * ish1_442[k]
                   + f_3 * pc_x[k] * isi_589[k];

        t_758[k] = f_3 * pc_z[k] * isi_588[k];

        t_759[k] = f_10 * ish0_444[k]
                   - f_11 * ish1_444[k]
                   + f_3 * pc_x[k] * isi_591[k];

        t_760[k] = f_3 * pc_z[k] * isi_589[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pc_x, pc_z, ish0_446, ish0_447, ish0_449, \
                         ish1_446, ish1_447, ish1_449, isi_591, isi_593, isi_594, \
                         isi_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_10 * ish0_446[k]
                   - f_11 * ish1_446[k]
                   + f_3 * pc_x[k] * isi_593[k];

        t_762[k] = f_8 * ish0_447[k]
                   - f_9 * ish1_447[k]
                   + f_3 * pc_x[k] * isi_594[k];

        t_763[k] = f_3 * pc_z[k] * isi_591[k];

        t_764[k] = f_8 * ish0_449[k]
                   - f_9 * ish1_449[k]
                   + f_3 * pc_x[k] * isi_596[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, pc_x, pc_z, ish0_450, ish0_451, ish0_453, \
                         ish1_450, ish1_451, ish1_453, isi_594, isi_597, isi_598, \
                         isi_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_8 * ish0_450[k]
                   - f_9 * ish1_450[k]
                   + f_3 * pc_x[k] * isi_597[k];

        t_766[k] = f_6 * ish0_451[k]
                   - f_7 * ish1_451[k]
                   + f_3 * pc_x[k] * isi_598[k];

        t_767[k] = f_3 * pc_z[k] * isi_594[k];

        t_768[k] = f_6 * ish0_453[k]
                   - f_7 * ish1_453[k]
                   + f_3 * pc_x[k] * isi_600[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, pc_x, pc_z, ish0_454, ish0_455, ish0_456, \
                         ish1_454, ish1_455, ish1_456, isi_598, isi_601, isi_602, \
                         isi_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_6 * ish0_454[k]
                   - f_7 * ish1_454[k]
                   + f_3 * pc_x[k] * isi_601[k];

        t_770[k] = f_6 * ish0_455[k]
                   - f_7 * ish1_455[k]
                   + f_3 * pc_x[k] * isi_602[k];

        t_771[k] = f_4 * ish0_456[k]
                   - f_5 * ish1_456[k]
                   + f_3 * pc_x[k] * isi_603[k];

        t_772[k] = f_3 * pc_z[k] * isi_598[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, ish0_458, ish0_459, ish0_460, ish1_458, \
                         ish1_459, ish1_460, isi_605, isi_606, \
                         isi_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_4 * ish0_458[k]
                   - f_5 * ish1_458[k]
                   + f_3 * pc_x[k] * isi_605[k];

        t_774[k] = f_4 * ish0_459[k]
                   - f_5 * ish1_459[k]
                   + f_3 * pc_x[k] * isi_606[k];

        t_775[k] = f_4 * ish0_460[k]
                   - f_5 * ish1_460[k]
                   + f_3 * pc_x[k] * isi_607[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, t_780, t_781, pc_x, ish0_461, ish1_461, \
                         isi_608, isi_609, isi_610, isi_611, isi_612, \
                         isi_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_4 * ish0_461[k]
                   - f_5 * ish1_461[k]
                   + f_3 * pc_x[k] * isi_608[k];

        t_777[k] = f_3 * pc_x[k] * isi_609[k];

        t_778[k] = f_3 * pc_x[k] * isi_610[k];

        t_779[k] = f_3 * pc_x[k] * isi_611[k];

        t_780[k] = f_3 * pc_x[k] * isi_612[k];

        t_781[k] = f_3 * pc_x[k] * isi_613[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, t_786, pc_x, pc_y, pc_z, hsi_441, \
                         ish0_456, ish1_456, isi_609, isi_610, isi_614, \
                         isi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_3 * pc_x[k] * isi_614[k];

        t_783[k] = f_3 * pc_x[k] * isi_615[k];

        t_784[k] = f_0 * hsi_441[k]
                   + f_1 * ish0_456[k]
                   - f_2 * ish1_456[k]
                   + f_3 * pc_y[k] * isi_609[k];

        t_785[k] = f_3 * pc_z[k] * isi_609[k];

        t_786[k] = f_4 * ish0_456[k]
                   - f_5 * ish1_456[k]
                   + f_3 * pc_z[k] * isi_610[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_z, ish0_457, ish0_458, ish0_459, ish1_457, \
                         ish1_458, ish1_459, isi_611, isi_612, \
                         isi_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_6 * ish0_457[k]
                   - f_7 * ish1_457[k]
                   + f_3 * pc_z[k] * isi_611[k];

        t_788[k] = f_8 * ish0_458[k]
                   - f_9 * ish1_458[k]
                   + f_3 * pc_z[k] * isi_612[k];

        t_789[k] = f_10 * ish0_459[k]
                   - f_11 * ish1_459[k]
                   + f_3 * pc_z[k] * isi_613[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_y, pc_z, hsk0_540, hsk0_541, \
                         hsi_447, hsk1_540, hsk1_541, ish0_461, ish1_461, \
                         isi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_0 * hsi_447[k]
                   + f_3 * pc_y[k] * isi_615[k];

        t_791[k] = f_1 * ish0_461[k]
                   - f_2 * ish1_461[k]
                   + f_3 * pc_z[k] * isi_615[k];

        t_792[k] = pa_z[k] * hsk0_540[k]
                   - f_12 * pc_z[k] * hsk1_540[k];

        t_793[k] = pa_z[k] * hsk0_541[k]
                   - f_12 * pc_z[k] * hsk1_541[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pa_z, pc_x, pc_z, hsk0_543, hsk1_543, ish0_464, \
                         ish0_466, ish1_464, ish1_466, isi_618, \
                         isi_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_18 * ish0_464[k]
                   - f_19 * ish1_464[k]
                   + f_3 * pc_x[k] * isi_618[k];

        t_795[k] = pa_z[k] * hsk0_543[k]
                   - f_12 * pc_z[k] * hsk1_543[k];

        t_796[k] = f_10 * ish0_466[k]
                   - f_11 * ish1_466[k]
                   + f_3 * pc_x[k] * isi_620[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pa_z, pc_x, pc_z, hsk0_546, hsk1_546, ish0_467, \
                         ish0_469, ish1_467, ish1_469, isi_621, \
                         isi_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_10 * ish0_467[k]
                   - f_11 * ish1_467[k]
                   + f_3 * pc_x[k] * isi_621[k];

        t_798[k] = pa_z[k] * hsk0_546[k]
                   - f_12 * pc_z[k] * hsk1_546[k];

        t_799[k] = f_8 * ish0_469[k]
                   - f_9 * ish1_469[k]
                   + f_3 * pc_x[k] * isi_623[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pa_z, pc_x, pc_z, hsk0_550, hsk1_550, ish0_470, \
                         ish0_471, ish1_470, ish1_471, isi_624, \
                         isi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_8 * ish0_470[k]
                   - f_9 * ish1_470[k]
                   + f_3 * pc_x[k] * isi_624[k];

        t_801[k] = f_8 * ish0_471[k]
                   - f_9 * ish1_471[k]
                   + f_3 * pc_x[k] * isi_625[k];

        t_802[k] = pa_z[k] * hsk0_550[k]
                   - f_12 * pc_z[k] * hsk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, ish0_473, ish0_474, ish0_475, ish1_473, \
                         ish1_474, ish1_475, isi_627, isi_628, \
                         isi_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_6 * ish0_473[k]
                   - f_7 * ish1_473[k]
                   + f_3 * pc_x[k] * isi_627[k];

        t_804[k] = f_6 * ish0_474[k]
                   - f_7 * ish1_474[k]
                   + f_3 * pc_x[k] * isi_628[k];

        t_805[k] = f_6 * ish0_475[k]
                   - f_7 * ish1_475[k]
                   + f_3 * pc_x[k] * isi_629[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pa_z, pc_x, pc_z, hsk0_555, hsk1_555, ish0_476, \
                         ish0_478, ish1_476, ish1_478, isi_630, \
                         isi_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_6 * ish0_476[k]
                   - f_7 * ish1_476[k]
                   + f_3 * pc_x[k] * isi_630[k];

        t_807[k] = pa_z[k] * hsk0_555[k]
                   - f_12 * pc_z[k] * hsk1_555[k];

        t_808[k] = f_4 * ish0_478[k]
                   - f_5 * ish1_478[k]
                   + f_3 * pc_x[k] * isi_632[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, ish0_479, ish0_480, ish0_481, ish1_479, \
                         ish1_480, ish1_481, isi_633, isi_634, \
                         isi_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_4 * ish0_479[k]
                   - f_5 * ish1_479[k]
                   + f_3 * pc_x[k] * isi_633[k];

        t_810[k] = f_4 * ish0_480[k]
                   - f_5 * ish1_480[k]
                   + f_3 * pc_x[k] * isi_634[k];

        t_811[k] = f_4 * ish0_481[k]
                   - f_5 * ish1_481[k]
                   + f_3 * pc_x[k] * isi_635[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, t_816, t_817, pc_x, ish0_482, ish1_482, \
                         isi_636, isi_637, isi_638, isi_639, isi_640, \
                         isi_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_4 * ish0_482[k]
                   - f_5 * ish1_482[k]
                   + f_3 * pc_x[k] * isi_636[k];

        t_813[k] = f_3 * pc_x[k] * isi_637[k];

        t_814[k] = f_3 * pc_x[k] * isi_638[k];

        t_815[k] = f_3 * pc_x[k] * isi_639[k];

        t_816[k] = f_3 * pc_x[k] * isi_640[k];

        t_817[k] = f_3 * pc_x[k] * isi_641[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, pa_z, pc_x, pc_z, hsk0_568, hsi_441, \
                         hsk1_568, isi_637, isi_642, isi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_3 * pc_x[k] * isi_642[k];

        t_819[k] = f_3 * pc_x[k] * isi_643[k];

        t_820[k] = pa_z[k] * hsk0_568[k]
                   - f_12 * pc_z[k] * hsk1_568[k];

        t_821[k] = f_13 * hsi_441[k]
                   + f_3 * pc_z[k] * isi_637[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pa_z, pc_z, hsk0_570, hsk0_571, hsk0_572, \
                         hsi_442, hsi_443, hsi_444, hsk1_570, hsk1_571, \
                         hsk1_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pa_z[k] * hsk0_570[k]
                   + f_14 * hsi_442[k]
                   - f_12 * pc_z[k] * hsk1_570[k];

        t_823[k] = pa_z[k] * hsk0_571[k]
                   + f_15 * hsi_443[k]
                   - f_12 * pc_z[k] * hsk1_571[k];

        t_824[k] = pa_z[k] * hsk0_572[k]
                   + f_16 * hsi_444[k]
                   - f_12 * pc_z[k] * hsk1_572[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pa_z, pc_y, pc_z, hsk0_573, hsi_445, hsi_447, \
                         hsi_475, hsk1_573, ish0_482, ish1_482, \
                         isi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = pa_z[k] * hsk0_573[k]
                   + f_17 * hsi_445[k]
                   - f_12 * pc_z[k] * hsk1_573[k];

        t_826[k] = f_17 * hsi_475[k]
                   + f_3 * pc_y[k] * isi_643[k];

        t_827[k] = f_13 * hsi_447[k]
                   + f_1 * ish0_482[k]
                   - f_2 * ish1_482[k]
                   + f_3 * pc_z[k] * isi_643[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, pc_x, ish0_483, ish0_484, ish0_485, ish1_483, \
                         ish1_484, ish1_485, isi_644, isi_645, \
                         isi_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = f_1 * ish0_483[k]
                   - f_2 * ish1_483[k]
                   + f_3 * pc_x[k] * isi_644[k];

        t_829[k] = f_18 * ish0_484[k]
                   - f_19 * ish1_484[k]
                   + f_3 * pc_x[k] * isi_645[k];

        t_830[k] = f_18 * ish0_485[k]
                   - f_19 * ish1_485[k]
                   + f_3 * pc_x[k] * isi_646[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, pc_x, ish0_486, ish0_487, ish0_488, ish1_486, \
                         ish1_487, ish1_488, isi_647, isi_648, \
                         isi_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_10 * ish0_486[k]
                   - f_11 * ish1_486[k]
                   + f_3 * pc_x[k] * isi_647[k];

        t_832[k] = f_10 * ish0_487[k]
                   - f_11 * ish1_487[k]
                   + f_3 * pc_x[k] * isi_648[k];

        t_833[k] = f_10 * ish0_488[k]
                   - f_11 * ish1_488[k]
                   + f_3 * pc_x[k] * isi_649[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_x, ish0_489, ish0_490, ish0_491, ish1_489, \
                         ish1_490, ish1_491, isi_650, isi_651, \
                         isi_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_8 * ish0_489[k]
                   - f_9 * ish1_489[k]
                   + f_3 * pc_x[k] * isi_650[k];

        t_835[k] = f_8 * ish0_490[k]
                   - f_9 * ish1_490[k]
                   + f_3 * pc_x[k] * isi_651[k];

        t_836[k] = f_8 * ish0_491[k]
                   - f_9 * ish1_491[k]
                   + f_3 * pc_x[k] * isi_652[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_x, ish0_492, ish0_493, ish0_494, ish1_492, \
                         ish1_493, ish1_494, isi_653, isi_654, \
                         isi_655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_8 * ish0_492[k]
                   - f_9 * ish1_492[k]
                   + f_3 * pc_x[k] * isi_653[k];

        t_838[k] = f_6 * ish0_493[k]
                   - f_7 * ish1_493[k]
                   + f_3 * pc_x[k] * isi_654[k];

        t_839[k] = f_6 * ish0_494[k]
                   - f_7 * ish1_494[k]
                   + f_3 * pc_x[k] * isi_655[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, ish0_495, ish0_496, ish0_497, ish1_495, \
                         ish1_496, ish1_497, isi_656, isi_657, \
                         isi_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_6 * ish0_495[k]
                   - f_7 * ish1_495[k]
                   + f_3 * pc_x[k] * isi_656[k];

        t_841[k] = f_6 * ish0_496[k]
                   - f_7 * ish1_496[k]
                   + f_3 * pc_x[k] * isi_657[k];

        t_842[k] = f_6 * ish0_497[k]
                   - f_7 * ish1_497[k]
                   + f_3 * pc_x[k] * isi_658[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
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
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);

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
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsk0_720 = buffer.data(hsk0 + 720);
    const auto *hsk0_722 = buffer.data(hsk0 + 722);
    const auto *hsk0_725 = buffer.data(hsk0 + 725);
    const auto *hsk0_729 = buffer.data(hsk0 + 729);
    const auto *hsk0_734 = buffer.data(hsk0 + 734);

    const auto *hsi_469 = buffer.data(hsi + 469);
    const auto *hsi_475 = buffer.data(hsi + 475);
    const auto *hsi_497 = buffer.data(hsi + 497);
    const auto *hsi_499 = buffer.data(hsi + 499);
    const auto *hsi_500 = buffer.data(hsi + 500);
    const auto *hsi_501 = buffer.data(hsi + 501);
    const auto *hsi_502 = buffer.data(hsi + 502);
    const auto *hsi_503 = buffer.data(hsi + 503);
    const auto *hsi_525 = buffer.data(hsi + 525);
    const auto *hsi_527 = buffer.data(hsi + 527);
    const auto *hsi_528 = buffer.data(hsi + 528);
    const auto *hsi_529 = buffer.data(hsi + 529);
    const auto *hsi_530 = buffer.data(hsi + 530);
    const auto *hsi_531 = buffer.data(hsi + 531);
    const auto *hsi_553 = buffer.data(hsi + 553);
    const auto *hsi_555 = buffer.data(hsi + 555);
    const auto *hsi_556 = buffer.data(hsi + 556);
    const auto *hsi_557 = buffer.data(hsi + 557);
    const auto *hsi_558 = buffer.data(hsi + 558);
    const auto *hsi_559 = buffer.data(hsi + 559);

    const auto *hsk1_720 = buffer.data(hsk1 + 720);
    const auto *hsk1_722 = buffer.data(hsk1 + 722);
    const auto *hsk1_725 = buffer.data(hsk1 + 725);
    const auto *hsk1_729 = buffer.data(hsk1 + 729);
    const auto *hsk1_734 = buffer.data(hsk1 + 734);

    const auto *ish0_498 = buffer.data(ish0 + 498);
    const auto *ish0_499 = buffer.data(ish0 + 499);
    const auto *ish0_500 = buffer.data(ish0 + 500);
    const auto *ish0_501 = buffer.data(ish0 + 501);
    const auto *ish0_502 = buffer.data(ish0 + 502);
    const auto *ish0_503 = buffer.data(ish0 + 503);
    const auto *ish0_504 = buffer.data(ish0 + 504);
    const auto *ish0_505 = buffer.data(ish0 + 505);
    const auto *ish0_506 = buffer.data(ish0 + 506);
    const auto *ish0_507 = buffer.data(ish0 + 507);
    const auto *ish0_508 = buffer.data(ish0 + 508);
    const auto *ish0_509 = buffer.data(ish0 + 509);
    const auto *ish0_510 = buffer.data(ish0 + 510);
    const auto *ish0_511 = buffer.data(ish0 + 511);
    const auto *ish0_512 = buffer.data(ish0 + 512);
    const auto *ish0_513 = buffer.data(ish0 + 513);
    const auto *ish0_514 = buffer.data(ish0 + 514);
    const auto *ish0_515 = buffer.data(ish0 + 515);
    const auto *ish0_516 = buffer.data(ish0 + 516);
    const auto *ish0_517 = buffer.data(ish0 + 517);
    const auto *ish0_518 = buffer.data(ish0 + 518);
    const auto *ish0_519 = buffer.data(ish0 + 519);
    const auto *ish0_520 = buffer.data(ish0 + 520);
    const auto *ish0_521 = buffer.data(ish0 + 521);
    const auto *ish0_522 = buffer.data(ish0 + 522);
    const auto *ish0_523 = buffer.data(ish0 + 523);
    const auto *ish0_524 = buffer.data(ish0 + 524);
    const auto *ish0_525 = buffer.data(ish0 + 525);
    const auto *ish0_526 = buffer.data(ish0 + 526);
    const auto *ish0_527 = buffer.data(ish0 + 527);
    const auto *ish0_528 = buffer.data(ish0 + 528);
    const auto *ish0_529 = buffer.data(ish0 + 529);
    const auto *ish0_530 = buffer.data(ish0 + 530);
    const auto *ish0_531 = buffer.data(ish0 + 531);
    const auto *ish0_532 = buffer.data(ish0 + 532);
    const auto *ish0_533 = buffer.data(ish0 + 533);
    const auto *ish0_534 = buffer.data(ish0 + 534);
    const auto *ish0_535 = buffer.data(ish0 + 535);
    const auto *ish0_536 = buffer.data(ish0 + 536);
    const auto *ish0_537 = buffer.data(ish0 + 537);
    const auto *ish0_538 = buffer.data(ish0 + 538);
    const auto *ish0_539 = buffer.data(ish0 + 539);
    const auto *ish0_540 = buffer.data(ish0 + 540);
    const auto *ish0_541 = buffer.data(ish0 + 541);
    const auto *ish0_542 = buffer.data(ish0 + 542);
    const auto *ish0_543 = buffer.data(ish0 + 543);
    const auto *ish0_544 = buffer.data(ish0 + 544);
    const auto *ish0_545 = buffer.data(ish0 + 545);
    const auto *ish0_547 = buffer.data(ish0 + 547);
    const auto *ish0_549 = buffer.data(ish0 + 549);
    const auto *ish0_550 = buffer.data(ish0 + 550);
    const auto *ish0_552 = buffer.data(ish0 + 552);
    const auto *ish0_553 = buffer.data(ish0 + 553);
    const auto *ish0_554 = buffer.data(ish0 + 554);
    const auto *ish0_556 = buffer.data(ish0 + 556);
    const auto *ish0_557 = buffer.data(ish0 + 557);
    const auto *ish0_558 = buffer.data(ish0 + 558);
    const auto *ish0_559 = buffer.data(ish0 + 559);
    const auto *ish0_561 = buffer.data(ish0 + 561);
    const auto *ish0_562 = buffer.data(ish0 + 562);
    const auto *ish0_563 = buffer.data(ish0 + 563);
    const auto *ish0_564 = buffer.data(ish0 + 564);

    const auto *ish1_498 = buffer.data(ish1 + 498);
    const auto *ish1_499 = buffer.data(ish1 + 499);
    const auto *ish1_500 = buffer.data(ish1 + 500);
    const auto *ish1_501 = buffer.data(ish1 + 501);
    const auto *ish1_502 = buffer.data(ish1 + 502);
    const auto *ish1_503 = buffer.data(ish1 + 503);
    const auto *ish1_504 = buffer.data(ish1 + 504);
    const auto *ish1_505 = buffer.data(ish1 + 505);
    const auto *ish1_506 = buffer.data(ish1 + 506);
    const auto *ish1_507 = buffer.data(ish1 + 507);
    const auto *ish1_508 = buffer.data(ish1 + 508);
    const auto *ish1_509 = buffer.data(ish1 + 509);
    const auto *ish1_510 = buffer.data(ish1 + 510);
    const auto *ish1_511 = buffer.data(ish1 + 511);
    const auto *ish1_512 = buffer.data(ish1 + 512);
    const auto *ish1_513 = buffer.data(ish1 + 513);
    const auto *ish1_514 = buffer.data(ish1 + 514);
    const auto *ish1_515 = buffer.data(ish1 + 515);
    const auto *ish1_516 = buffer.data(ish1 + 516);
    const auto *ish1_517 = buffer.data(ish1 + 517);
    const auto *ish1_518 = buffer.data(ish1 + 518);
    const auto *ish1_519 = buffer.data(ish1 + 519);
    const auto *ish1_520 = buffer.data(ish1 + 520);
    const auto *ish1_521 = buffer.data(ish1 + 521);
    const auto *ish1_522 = buffer.data(ish1 + 522);
    const auto *ish1_523 = buffer.data(ish1 + 523);
    const auto *ish1_524 = buffer.data(ish1 + 524);
    const auto *ish1_525 = buffer.data(ish1 + 525);
    const auto *ish1_526 = buffer.data(ish1 + 526);
    const auto *ish1_527 = buffer.data(ish1 + 527);
    const auto *ish1_528 = buffer.data(ish1 + 528);
    const auto *ish1_529 = buffer.data(ish1 + 529);
    const auto *ish1_530 = buffer.data(ish1 + 530);
    const auto *ish1_531 = buffer.data(ish1 + 531);
    const auto *ish1_532 = buffer.data(ish1 + 532);
    const auto *ish1_533 = buffer.data(ish1 + 533);
    const auto *ish1_534 = buffer.data(ish1 + 534);
    const auto *ish1_535 = buffer.data(ish1 + 535);
    const auto *ish1_536 = buffer.data(ish1 + 536);
    const auto *ish1_537 = buffer.data(ish1 + 537);
    const auto *ish1_538 = buffer.data(ish1 + 538);
    const auto *ish1_539 = buffer.data(ish1 + 539);
    const auto *ish1_540 = buffer.data(ish1 + 540);
    const auto *ish1_541 = buffer.data(ish1 + 541);
    const auto *ish1_542 = buffer.data(ish1 + 542);
    const auto *ish1_543 = buffer.data(ish1 + 543);
    const auto *ish1_544 = buffer.data(ish1 + 544);
    const auto *ish1_545 = buffer.data(ish1 + 545);
    const auto *ish1_547 = buffer.data(ish1 + 547);
    const auto *ish1_549 = buffer.data(ish1 + 549);
    const auto *ish1_550 = buffer.data(ish1 + 550);
    const auto *ish1_552 = buffer.data(ish1 + 552);
    const auto *ish1_553 = buffer.data(ish1 + 553);
    const auto *ish1_554 = buffer.data(ish1 + 554);
    const auto *ish1_556 = buffer.data(ish1 + 556);
    const auto *ish1_557 = buffer.data(ish1 + 557);
    const auto *ish1_558 = buffer.data(ish1 + 558);
    const auto *ish1_559 = buffer.data(ish1 + 559);
    const auto *ish1_561 = buffer.data(ish1 + 561);
    const auto *ish1_562 = buffer.data(ish1 + 562);
    const auto *ish1_563 = buffer.data(ish1 + 563);
    const auto *ish1_564 = buffer.data(ish1 + 564);

    const auto *isi_659 = buffer.data(isi + 659);
    const auto *isi_660 = buffer.data(isi + 660);
    const auto *isi_661 = buffer.data(isi + 661);
    const auto *isi_662 = buffer.data(isi + 662);
    const auto *isi_663 = buffer.data(isi + 663);
    const auto *isi_664 = buffer.data(isi + 664);
    const auto *isi_665 = buffer.data(isi + 665);
    const auto *isi_666 = buffer.data(isi + 666);
    const auto *isi_667 = buffer.data(isi + 667);
    const auto *isi_668 = buffer.data(isi + 668);
    const auto *isi_669 = buffer.data(isi + 669);
    const auto *isi_670 = buffer.data(isi + 670);
    const auto *isi_671 = buffer.data(isi + 671);
    const auto *isi_672 = buffer.data(isi + 672);
    const auto *isi_673 = buffer.data(isi + 673);
    const auto *isi_674 = buffer.data(isi + 674);
    const auto *isi_675 = buffer.data(isi + 675);
    const auto *isi_676 = buffer.data(isi + 676);
    const auto *isi_677 = buffer.data(isi + 677);
    const auto *isi_678 = buffer.data(isi + 678);
    const auto *isi_679 = buffer.data(isi + 679);
    const auto *isi_680 = buffer.data(isi + 680);
    const auto *isi_681 = buffer.data(isi + 681);
    const auto *isi_682 = buffer.data(isi + 682);
    const auto *isi_683 = buffer.data(isi + 683);
    const auto *isi_684 = buffer.data(isi + 684);
    const auto *isi_685 = buffer.data(isi + 685);
    const auto *isi_686 = buffer.data(isi + 686);
    const auto *isi_687 = buffer.data(isi + 687);
    const auto *isi_688 = buffer.data(isi + 688);
    const auto *isi_689 = buffer.data(isi + 689);
    const auto *isi_690 = buffer.data(isi + 690);
    const auto *isi_691 = buffer.data(isi + 691);
    const auto *isi_692 = buffer.data(isi + 692);
    const auto *isi_693 = buffer.data(isi + 693);
    const auto *isi_694 = buffer.data(isi + 694);
    const auto *isi_695 = buffer.data(isi + 695);
    const auto *isi_696 = buffer.data(isi + 696);
    const auto *isi_697 = buffer.data(isi + 697);
    const auto *isi_698 = buffer.data(isi + 698);
    const auto *isi_699 = buffer.data(isi + 699);
    const auto *isi_700 = buffer.data(isi + 700);
    const auto *isi_701 = buffer.data(isi + 701);
    const auto *isi_702 = buffer.data(isi + 702);
    const auto *isi_703 = buffer.data(isi + 703);
    const auto *isi_704 = buffer.data(isi + 704);
    const auto *isi_705 = buffer.data(isi + 705);
    const auto *isi_706 = buffer.data(isi + 706);
    const auto *isi_707 = buffer.data(isi + 707);
    const auto *isi_708 = buffer.data(isi + 708);
    const auto *isi_709 = buffer.data(isi + 709);
    const auto *isi_710 = buffer.data(isi + 710);
    const auto *isi_711 = buffer.data(isi + 711);
    const auto *isi_712 = buffer.data(isi + 712);
    const auto *isi_713 = buffer.data(isi + 713);
    const auto *isi_714 = buffer.data(isi + 714);
    const auto *isi_715 = buffer.data(isi + 715);
    const auto *isi_716 = buffer.data(isi + 716);
    const auto *isi_717 = buffer.data(isi + 717);
    const auto *isi_718 = buffer.data(isi + 718);
    const auto *isi_719 = buffer.data(isi + 719);
    const auto *isi_720 = buffer.data(isi + 720);
    const auto *isi_721 = buffer.data(isi + 721);
    const auto *isi_722 = buffer.data(isi + 722);
    const auto *isi_723 = buffer.data(isi + 723);
    const auto *isi_724 = buffer.data(isi + 724);
    const auto *isi_725 = buffer.data(isi + 725);
    const auto *isi_726 = buffer.data(isi + 726);
    const auto *isi_727 = buffer.data(isi + 727);
    const auto *isi_729 = buffer.data(isi + 729);
    const auto *isi_731 = buffer.data(isi + 731);
    const auto *isi_732 = buffer.data(isi + 732);
    const auto *isi_734 = buffer.data(isi + 734);
    const auto *isi_735 = buffer.data(isi + 735);
    const auto *isi_736 = buffer.data(isi + 736);
    const auto *isi_738 = buffer.data(isi + 738);
    const auto *isi_739 = buffer.data(isi + 739);
    const auto *isi_740 = buffer.data(isi + 740);
    const auto *isi_741 = buffer.data(isi + 741);
    const auto *isi_743 = buffer.data(isi + 743);
    const auto *isi_744 = buffer.data(isi + 744);
    const auto *isi_745 = buffer.data(isi + 745);
    const auto *isi_746 = buffer.data(isi + 746);

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, ish0_498, ish0_499, ish0_500, ish1_498, \
                         ish1_499, ish1_500, isi_659, isi_660, \
                         isi_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_4 * ish0_498[k]
                   - f_5 * ish1_498[k]
                   + f_3 * pc_x[k] * isi_659[k];

        t_844[k] = f_4 * ish0_499[k]
                   - f_5 * ish1_499[k]
                   + f_3 * pc_x[k] * isi_660[k];

        t_845[k] = f_4 * ish0_500[k]
                   - f_5 * ish1_500[k]
                   + f_3 * pc_x[k] * isi_661[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, t_849, pc_x, ish0_501, ish0_502, ish0_503, \
                         ish1_501, ish1_502, ish1_503, isi_662, isi_663, isi_664, \
                         isi_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_4 * ish0_501[k]
                   - f_5 * ish1_501[k]
                   + f_3 * pc_x[k] * isi_662[k];

        t_847[k] = f_4 * ish0_502[k]
                   - f_5 * ish1_502[k]
                   + f_3 * pc_x[k] * isi_663[k];

        t_848[k] = f_4 * ish0_503[k]
                   - f_5 * ish1_503[k]
                   + f_3 * pc_x[k] * isi_664[k];

        t_849[k] = f_3 * pc_x[k] * isi_665[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, t_855, pc_x, isi_666, isi_667, \
                         isi_668, isi_669, isi_670, isi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_3 * pc_x[k] * isi_666[k];

        t_851[k] = f_3 * pc_x[k] * isi_667[k];

        t_852[k] = f_3 * pc_x[k] * isi_668[k];

        t_853[k] = f_3 * pc_x[k] * isi_669[k];

        t_854[k] = f_3 * pc_x[k] * isi_670[k];

        t_855[k] = f_3 * pc_x[k] * isi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, hsi_469, hsi_497, hsi_499, ish0_498, \
                         ish0_500, ish1_498, ish1_500, isi_665, \
                         isi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * hsi_497[k]
                   + f_1 * ish0_498[k]
                   - f_2 * ish1_498[k]
                   + f_3 * pc_y[k] * isi_665[k];

        t_857[k] = f_14 * hsi_469[k]
                   + f_3 * pc_z[k] * isi_665[k];

        t_858[k] = f_16 * hsi_499[k]
                   + f_10 * ish0_500[k]
                   - f_11 * ish1_500[k]
                   + f_3 * pc_y[k] * isi_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, hsi_500, hsi_501, hsi_502, ish0_501, \
                         ish0_502, ish0_503, ish1_501, ish1_502, ish1_503, isi_668, isi_669, \
                         isi_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * hsi_500[k]
                   + f_8 * ish0_501[k]
                   - f_9 * ish1_501[k]
                   + f_3 * pc_y[k] * isi_668[k];

        t_860[k] = f_16 * hsi_501[k]
                   + f_6 * ish0_502[k]
                   - f_7 * ish1_502[k]
                   + f_3 * pc_y[k] * isi_669[k];

        t_861[k] = f_16 * hsi_502[k]
                   + f_4 * ish0_503[k]
                   - f_5 * ish1_503[k]
                   + f_3 * pc_y[k] * isi_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, hsi_475, hsi_503, ish0_503, \
                         ish0_504, ish1_503, ish1_504, isi_671, \
                         isi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * hsi_503[k]
                   + f_3 * pc_y[k] * isi_671[k];

        t_863[k] = f_14 * hsi_475[k]
                   + f_1 * ish0_503[k]
                   - f_2 * ish1_503[k]
                   + f_3 * pc_z[k] * isi_671[k];

        t_864[k] = f_1 * ish0_504[k]
                   - f_2 * ish1_504[k]
                   + f_3 * pc_x[k] * isi_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, pc_x, ish0_505, ish0_506, ish0_507, ish1_505, \
                         ish1_506, ish1_507, isi_673, isi_674, \
                         isi_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_18 * ish0_505[k]
                   - f_19 * ish1_505[k]
                   + f_3 * pc_x[k] * isi_673[k];

        t_866[k] = f_18 * ish0_506[k]
                   - f_19 * ish1_506[k]
                   + f_3 * pc_x[k] * isi_674[k];

        t_867[k] = f_10 * ish0_507[k]
                   - f_11 * ish1_507[k]
                   + f_3 * pc_x[k] * isi_675[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, pc_x, ish0_508, ish0_509, ish0_510, ish1_508, \
                         ish1_509, ish1_510, isi_676, isi_677, \
                         isi_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_10 * ish0_508[k]
                   - f_11 * ish1_508[k]
                   + f_3 * pc_x[k] * isi_676[k];

        t_869[k] = f_10 * ish0_509[k]
                   - f_11 * ish1_509[k]
                   + f_3 * pc_x[k] * isi_677[k];

        t_870[k] = f_8 * ish0_510[k]
                   - f_9 * ish1_510[k]
                   + f_3 * pc_x[k] * isi_678[k];
    }

#pragma omp simd aligned(t_871, t_872, t_873, pc_x, ish0_511, ish0_512, ish0_513, ish1_511, \
                         ish1_512, ish1_513, isi_679, isi_680, \
                         isi_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_871[k] = f_8 * ish0_511[k]
                   - f_9 * ish1_511[k]
                   + f_3 * pc_x[k] * isi_679[k];

        t_872[k] = f_8 * ish0_512[k]
                   - f_9 * ish1_512[k]
                   + f_3 * pc_x[k] * isi_680[k];

        t_873[k] = f_8 * ish0_513[k]
                   - f_9 * ish1_513[k]
                   + f_3 * pc_x[k] * isi_681[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, pc_x, ish0_514, ish0_515, ish0_516, ish1_514, \
                         ish1_515, ish1_516, isi_682, isi_683, \
                         isi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_6 * ish0_514[k]
                   - f_7 * ish1_514[k]
                   + f_3 * pc_x[k] * isi_682[k];

        t_875[k] = f_6 * ish0_515[k]
                   - f_7 * ish1_515[k]
                   + f_3 * pc_x[k] * isi_683[k];

        t_876[k] = f_6 * ish0_516[k]
                   - f_7 * ish1_516[k]
                   + f_3 * pc_x[k] * isi_684[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_x, ish0_517, ish0_518, ish0_519, ish1_517, \
                         ish1_518, ish1_519, isi_685, isi_686, \
                         isi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_6 * ish0_517[k]
                   - f_7 * ish1_517[k]
                   + f_3 * pc_x[k] * isi_685[k];

        t_878[k] = f_6 * ish0_518[k]
                   - f_7 * ish1_518[k]
                   + f_3 * pc_x[k] * isi_686[k];

        t_879[k] = f_4 * ish0_519[k]
                   - f_5 * ish1_519[k]
                   + f_3 * pc_x[k] * isi_687[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_x, ish0_520, ish0_521, ish0_522, ish1_520, \
                         ish1_521, ish1_522, isi_688, isi_689, \
                         isi_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_4 * ish0_520[k]
                   - f_5 * ish1_520[k]
                   + f_3 * pc_x[k] * isi_688[k];

        t_881[k] = f_4 * ish0_521[k]
                   - f_5 * ish1_521[k]
                   + f_3 * pc_x[k] * isi_689[k];

        t_882[k] = f_4 * ish0_522[k]
                   - f_5 * ish1_522[k]
                   + f_3 * pc_x[k] * isi_690[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, t_887, pc_x, ish0_523, ish0_524, \
                         ish1_523, ish1_524, isi_691, isi_692, isi_693, isi_694, \
                         isi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_4 * ish0_523[k]
                   - f_5 * ish1_523[k]
                   + f_3 * pc_x[k] * isi_691[k];

        t_884[k] = f_4 * ish0_524[k]
                   - f_5 * ish1_524[k]
                   + f_3 * pc_x[k] * isi_692[k];

        t_885[k] = f_3 * pc_x[k] * isi_693[k];

        t_886[k] = f_3 * pc_x[k] * isi_694[k];

        t_887[k] = f_3 * pc_x[k] * isi_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, t_892, pc_x, pc_y, hsi_525, ish0_519, \
                         ish1_519, isi_693, isi_696, isi_697, isi_698, \
                         isi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_3 * pc_x[k] * isi_696[k];

        t_889[k] = f_3 * pc_x[k] * isi_697[k];

        t_890[k] = f_3 * pc_x[k] * isi_698[k];

        t_891[k] = f_3 * pc_x[k] * isi_699[k];

        t_892[k] = f_15 * hsi_525[k]
                   + f_1 * ish0_519[k]
                   - f_2 * ish1_519[k]
                   + f_3 * pc_y[k] * isi_693[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, hsi_497, hsi_527, hsi_528, ish0_521, \
                         ish0_522, ish1_521, ish1_522, isi_693, isi_695, \
                         isi_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * hsi_497[k]
                   + f_3 * pc_z[k] * isi_693[k];

        t_894[k] = f_15 * hsi_527[k]
                   + f_10 * ish0_521[k]
                   - f_11 * ish1_521[k]
                   + f_3 * pc_y[k] * isi_695[k];

        t_895[k] = f_15 * hsi_528[k]
                   + f_8 * ish0_522[k]
                   - f_9 * ish1_522[k]
                   + f_3 * pc_y[k] * isi_696[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, hsi_529, hsi_530, hsi_531, ish0_523, \
                         ish0_524, ish1_523, ish1_524, isi_697, isi_698, \
                         isi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * hsi_529[k]
                   + f_6 * ish0_523[k]
                   - f_7 * ish1_523[k]
                   + f_3 * pc_y[k] * isi_697[k];

        t_897[k] = f_15 * hsi_530[k]
                   + f_4 * ish0_524[k]
                   - f_5 * ish1_524[k]
                   + f_3 * pc_y[k] * isi_698[k];

        t_898[k] = f_15 * hsi_531[k]
                   + f_3 * pc_y[k] * isi_699[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_z, hsi_503, ish0_524, ish0_525, \
                         ish0_526, ish1_524, ish1_525, ish1_526, isi_699, isi_700, \
                         isi_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_15 * hsi_503[k]
                   + f_1 * ish0_524[k]
                   - f_2 * ish1_524[k]
                   + f_3 * pc_z[k] * isi_699[k];

        t_900[k] = f_1 * ish0_525[k]
                   - f_2 * ish1_525[k]
                   + f_3 * pc_x[k] * isi_700[k];

        t_901[k] = f_18 * ish0_526[k]
                   - f_19 * ish1_526[k]
                   + f_3 * pc_x[k] * isi_701[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, ish0_527, ish0_528, ish0_529, ish1_527, \
                         ish1_528, ish1_529, isi_702, isi_703, \
                         isi_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_18 * ish0_527[k]
                   - f_19 * ish1_527[k]
                   + f_3 * pc_x[k] * isi_702[k];

        t_903[k] = f_10 * ish0_528[k]
                   - f_11 * ish1_528[k]
                   + f_3 * pc_x[k] * isi_703[k];

        t_904[k] = f_10 * ish0_529[k]
                   - f_11 * ish1_529[k]
                   + f_3 * pc_x[k] * isi_704[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, ish0_530, ish0_531, ish0_532, ish1_530, \
                         ish1_531, ish1_532, isi_705, isi_706, \
                         isi_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_10 * ish0_530[k]
                   - f_11 * ish1_530[k]
                   + f_3 * pc_x[k] * isi_705[k];

        t_906[k] = f_8 * ish0_531[k]
                   - f_9 * ish1_531[k]
                   + f_3 * pc_x[k] * isi_706[k];

        t_907[k] = f_8 * ish0_532[k]
                   - f_9 * ish1_532[k]
                   + f_3 * pc_x[k] * isi_707[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, ish0_533, ish0_534, ish0_535, ish1_533, \
                         ish1_534, ish1_535, isi_708, isi_709, \
                         isi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_8 * ish0_533[k]
                   - f_9 * ish1_533[k]
                   + f_3 * pc_x[k] * isi_708[k];

        t_909[k] = f_8 * ish0_534[k]
                   - f_9 * ish1_534[k]
                   + f_3 * pc_x[k] * isi_709[k];

        t_910[k] = f_6 * ish0_535[k]
                   - f_7 * ish1_535[k]
                   + f_3 * pc_x[k] * isi_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, ish0_536, ish0_537, ish0_538, ish1_536, \
                         ish1_537, ish1_538, isi_711, isi_712, \
                         isi_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_6 * ish0_536[k]
                   - f_7 * ish1_536[k]
                   + f_3 * pc_x[k] * isi_711[k];

        t_912[k] = f_6 * ish0_537[k]
                   - f_7 * ish1_537[k]
                   + f_3 * pc_x[k] * isi_712[k];

        t_913[k] = f_6 * ish0_538[k]
                   - f_7 * ish1_538[k]
                   + f_3 * pc_x[k] * isi_713[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, ish0_539, ish0_540, ish0_541, ish1_539, \
                         ish1_540, ish1_541, isi_714, isi_715, \
                         isi_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_6 * ish0_539[k]
                   - f_7 * ish1_539[k]
                   + f_3 * pc_x[k] * isi_714[k];

        t_915[k] = f_4 * ish0_540[k]
                   - f_5 * ish1_540[k]
                   + f_3 * pc_x[k] * isi_715[k];

        t_916[k] = f_4 * ish0_541[k]
                   - f_5 * ish1_541[k]
                   + f_3 * pc_x[k] * isi_716[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, ish0_542, ish0_543, ish0_544, ish1_542, \
                         ish1_543, ish1_544, isi_717, isi_718, \
                         isi_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_4 * ish0_542[k]
                   - f_5 * ish1_542[k]
                   + f_3 * pc_x[k] * isi_717[k];

        t_918[k] = f_4 * ish0_543[k]
                   - f_5 * ish1_543[k]
                   + f_3 * pc_x[k] * isi_718[k];

        t_919[k] = f_4 * ish0_544[k]
                   - f_5 * ish1_544[k]
                   + f_3 * pc_x[k] * isi_719[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, t_925, pc_x, ish0_545, ish1_545, \
                         isi_720, isi_721, isi_722, isi_723, isi_724, \
                         isi_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_4 * ish0_545[k]
                   - f_5 * ish1_545[k]
                   + f_3 * pc_x[k] * isi_720[k];

        t_921[k] = f_3 * pc_x[k] * isi_721[k];

        t_922[k] = f_3 * pc_x[k] * isi_722[k];

        t_923[k] = f_3 * pc_x[k] * isi_723[k];

        t_924[k] = f_3 * pc_x[k] * isi_724[k];

        t_925[k] = f_3 * pc_x[k] * isi_725[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, t_929, pc_x, pc_y, pc_z, hsi_525, hsi_553, \
                         ish0_540, ish1_540, isi_721, isi_726, \
                         isi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_3 * pc_x[k] * isi_726[k];

        t_927[k] = f_3 * pc_x[k] * isi_727[k];

        t_928[k] = f_14 * hsi_553[k]
                   + f_1 * ish0_540[k]
                   - f_2 * ish1_540[k]
                   + f_3 * pc_y[k] * isi_721[k];

        t_929[k] = f_16 * hsi_525[k]
                   + f_3 * pc_z[k] * isi_721[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, pc_y, hsi_555, hsi_556, hsi_557, ish0_542, \
                         ish0_543, ish0_544, ish1_542, ish1_543, ish1_544, isi_723, isi_724, \
                         isi_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_14 * hsi_555[k]
                   + f_10 * ish0_542[k]
                   - f_11 * ish1_542[k]
                   + f_3 * pc_y[k] * isi_723[k];

        t_931[k] = f_14 * hsi_556[k]
                   + f_8 * ish0_543[k]
                   - f_9 * ish1_543[k]
                   + f_3 * pc_y[k] * isi_724[k];

        t_932[k] = f_14 * hsi_557[k]
                   + f_6 * ish0_544[k]
                   - f_7 * ish1_544[k]
                   + f_3 * pc_y[k] * isi_725[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, pa_y, pc_y, pc_z, hsk0_720, hsi_531, \
                         hsi_558, hsi_559, hsk1_720, ish0_545, ish1_545, isi_726, \
                         isi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_14 * hsi_558[k]
                   + f_4 * ish0_545[k]
                   - f_5 * ish1_545[k]
                   + f_3 * pc_y[k] * isi_726[k];

        t_934[k] = f_14 * hsi_559[k]
                   + f_3 * pc_y[k] * isi_727[k];

        t_935[k] = f_16 * hsi_531[k]
                   + f_1 * ish0_545[k]
                   - f_2 * ish1_545[k]
                   + f_3 * pc_z[k] * isi_727[k];

        t_936[k] = pa_y[k] * hsk0_720[k]
                   - f_12 * pc_y[k] * hsk1_720[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, pa_y, pc_x, pc_y, hsk0_722, hsk1_722, ish0_547, \
                         ish0_549, ish1_547, ish1_549, isi_729, \
                         isi_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = f_18 * ish0_547[k]
                   - f_19 * ish1_547[k]
                   + f_3 * pc_x[k] * isi_729[k];

        t_938[k] = pa_y[k] * hsk0_722[k]
                   - f_12 * pc_y[k] * hsk1_722[k];

        t_939[k] = f_10 * ish0_549[k]
                   - f_11 * ish1_549[k]
                   + f_3 * pc_x[k] * isi_731[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, pa_y, pc_x, pc_y, hsk0_725, hsk1_725, ish0_550, \
                         ish0_552, ish1_550, ish1_552, isi_732, \
                         isi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_10 * ish0_550[k]
                   - f_11 * ish1_550[k]
                   + f_3 * pc_x[k] * isi_732[k];

        t_941[k] = pa_y[k] * hsk0_725[k]
                   - f_12 * pc_y[k] * hsk1_725[k];

        t_942[k] = f_8 * ish0_552[k]
                   - f_9 * ish1_552[k]
                   + f_3 * pc_x[k] * isi_734[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, pa_y, pc_x, pc_y, hsk0_729, hsk1_729, ish0_553, \
                         ish0_554, ish1_553, ish1_554, isi_735, \
                         isi_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_8 * ish0_553[k]
                   - f_9 * ish1_553[k]
                   + f_3 * pc_x[k] * isi_735[k];

        t_944[k] = f_8 * ish0_554[k]
                   - f_9 * ish1_554[k]
                   + f_3 * pc_x[k] * isi_736[k];

        t_945[k] = pa_y[k] * hsk0_729[k]
                   - f_12 * pc_y[k] * hsk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_x, ish0_556, ish0_557, ish0_558, ish1_556, \
                         ish1_557, ish1_558, isi_738, isi_739, \
                         isi_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_6 * ish0_556[k]
                   - f_7 * ish1_556[k]
                   + f_3 * pc_x[k] * isi_738[k];

        t_947[k] = f_6 * ish0_557[k]
                   - f_7 * ish1_557[k]
                   + f_3 * pc_x[k] * isi_739[k];

        t_948[k] = f_6 * ish0_558[k]
                   - f_7 * ish1_558[k]
                   + f_3 * pc_x[k] * isi_740[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, pa_y, pc_x, pc_y, hsk0_734, hsk1_734, ish0_559, \
                         ish0_561, ish1_559, ish1_561, isi_741, \
                         isi_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_6 * ish0_559[k]
                   - f_7 * ish1_559[k]
                   + f_3 * pc_x[k] * isi_741[k];

        t_950[k] = pa_y[k] * hsk0_734[k]
                   - f_12 * pc_y[k] * hsk1_734[k];

        t_951[k] = f_4 * ish0_561[k]
                   - f_5 * ish1_561[k]
                   + f_3 * pc_x[k] * isi_743[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, pc_x, ish0_562, ish0_563, ish0_564, ish1_562, \
                         ish1_563, ish1_564, isi_744, isi_745, \
                         isi_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = f_4 * ish0_562[k]
                   - f_5 * ish1_562[k]
                   + f_3 * pc_x[k] * isi_744[k];

        t_953[k] = f_4 * ish0_563[k]
                   - f_5 * ish1_563[k]
                   + f_3 * pc_x[k] * isi_745[k];

        t_954[k] = f_4 * ish0_564[k]
                   - f_5 * ish1_564[k]
                   + f_3 * pc_x[k] * isi_746[k];
    }
}

static auto
compute_prim_isk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsk0,
                                                          const size_t hsi, const size_t hsk1,
                                                          const size_t ish0, const size_t ish1,
                                                          const size_t isi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_18 = 2.5 / gamma;
    const auto f_19 = 2.5 * p / (gamma * q);
    const auto f_20 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsk0_740 = buffer.data(hsk0 + 740);
    const auto *hsk0_748 = buffer.data(hsk0 + 748);
    const auto *hsk0_750 = buffer.data(hsk0 + 750);
    const auto *hsk0_751 = buffer.data(hsk0 + 751);
    const auto *hsk0_752 = buffer.data(hsk0 + 752);
    const auto *hsk0_753 = buffer.data(hsk0 + 753);
    const auto *hsk0_755 = buffer.data(hsk0 + 755);

    const auto *hsi_553 = buffer.data(hsi + 553);
    const auto *hsi_581 = buffer.data(hsi + 581);
    const auto *hsi_583 = buffer.data(hsi + 583);
    const auto *hsi_584 = buffer.data(hsi + 584);
    const auto *hsi_585 = buffer.data(hsi + 585);
    const auto *hsi_586 = buffer.data(hsi + 586);
    const auto *hsi_587 = buffer.data(hsi + 587);

    const auto *hsk1_740 = buffer.data(hsk1 + 740);
    const auto *hsk1_748 = buffer.data(hsk1 + 748);
    const auto *hsk1_750 = buffer.data(hsk1 + 750);
    const auto *hsk1_751 = buffer.data(hsk1 + 751);
    const auto *hsk1_752 = buffer.data(hsk1 + 752);
    const auto *hsk1_753 = buffer.data(hsk1 + 753);
    const auto *hsk1_755 = buffer.data(hsk1 + 755);

    const auto *ish0_565 = buffer.data(ish0 + 565);
    const auto *ish0_567 = buffer.data(ish0 + 567);
    const auto *ish0_569 = buffer.data(ish0 + 569);
    const auto *ish0_570 = buffer.data(ish0 + 570);
    const auto *ish0_572 = buffer.data(ish0 + 572);
    const auto *ish0_573 = buffer.data(ish0 + 573);
    const auto *ish0_574 = buffer.data(ish0 + 574);
    const auto *ish0_576 = buffer.data(ish0 + 576);
    const auto *ish0_577 = buffer.data(ish0 + 577);
    const auto *ish0_578 = buffer.data(ish0 + 578);
    const auto *ish0_579 = buffer.data(ish0 + 579);
    const auto *ish0_581 = buffer.data(ish0 + 581);
    const auto *ish0_582 = buffer.data(ish0 + 582);
    const auto *ish0_583 = buffer.data(ish0 + 583);
    const auto *ish0_584 = buffer.data(ish0 + 584);
    const auto *ish0_585 = buffer.data(ish0 + 585);
    const auto *ish0_586 = buffer.data(ish0 + 586);
    const auto *ish0_587 = buffer.data(ish0 + 587);

    const auto *ish1_565 = buffer.data(ish1 + 565);
    const auto *ish1_567 = buffer.data(ish1 + 567);
    const auto *ish1_569 = buffer.data(ish1 + 569);
    const auto *ish1_570 = buffer.data(ish1 + 570);
    const auto *ish1_572 = buffer.data(ish1 + 572);
    const auto *ish1_573 = buffer.data(ish1 + 573);
    const auto *ish1_574 = buffer.data(ish1 + 574);
    const auto *ish1_576 = buffer.data(ish1 + 576);
    const auto *ish1_577 = buffer.data(ish1 + 577);
    const auto *ish1_578 = buffer.data(ish1 + 578);
    const auto *ish1_579 = buffer.data(ish1 + 579);
    const auto *ish1_581 = buffer.data(ish1 + 581);
    const auto *ish1_582 = buffer.data(ish1 + 582);
    const auto *ish1_583 = buffer.data(ish1 + 583);
    const auto *ish1_584 = buffer.data(ish1 + 584);
    const auto *ish1_585 = buffer.data(ish1 + 585);
    const auto *ish1_586 = buffer.data(ish1 + 586);
    const auto *ish1_587 = buffer.data(ish1 + 587);

    const auto *isi_747 = buffer.data(isi + 747);
    const auto *isi_749 = buffer.data(isi + 749);
    const auto *isi_750 = buffer.data(isi + 750);
    const auto *isi_751 = buffer.data(isi + 751);
    const auto *isi_752 = buffer.data(isi + 752);
    const auto *isi_753 = buffer.data(isi + 753);
    const auto *isi_754 = buffer.data(isi + 754);
    const auto *isi_755 = buffer.data(isi + 755);
    const auto *isi_756 = buffer.data(isi + 756);
    const auto *isi_758 = buffer.data(isi + 758);
    const auto *isi_759 = buffer.data(isi + 759);
    const auto *isi_761 = buffer.data(isi + 761);
    const auto *isi_762 = buffer.data(isi + 762);
    const auto *isi_763 = buffer.data(isi + 763);
    const auto *isi_765 = buffer.data(isi + 765);
    const auto *isi_766 = buffer.data(isi + 766);
    const auto *isi_767 = buffer.data(isi + 767);
    const auto *isi_768 = buffer.data(isi + 768);
    const auto *isi_770 = buffer.data(isi + 770);
    const auto *isi_771 = buffer.data(isi + 771);
    const auto *isi_772 = buffer.data(isi + 772);
    const auto *isi_773 = buffer.data(isi + 773);
    const auto *isi_774 = buffer.data(isi + 774);
    const auto *isi_776 = buffer.data(isi + 776);
    const auto *isi_777 = buffer.data(isi + 777);
    const auto *isi_778 = buffer.data(isi + 778);
    const auto *isi_779 = buffer.data(isi + 779);
    const auto *isi_780 = buffer.data(isi + 780);
    const auto *isi_781 = buffer.data(isi + 781);
    const auto *isi_782 = buffer.data(isi + 782);
    const auto *isi_783 = buffer.data(isi + 783);

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, pa_y, pc_x, pc_y, hsk0_740, \
                         hsk1_740, ish0_565, ish1_565, isi_747, isi_749, isi_750, \
                         isi_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_4 * ish0_565[k]
                   - f_5 * ish1_565[k]
                   + f_3 * pc_x[k] * isi_747[k];

        t_956[k] = pa_y[k] * hsk0_740[k]
                   - f_12 * pc_y[k] * hsk1_740[k];

        t_957[k] = f_3 * pc_x[k] * isi_749[k];

        t_958[k] = f_3 * pc_x[k] * isi_750[k];

        t_959[k] = f_3 * pc_x[k] * isi_751[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, pa_y, pc_x, pc_y, hsk0_748, \
                         hsi_581, hsk1_748, isi_752, isi_753, isi_754, \
                         isi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_3 * pc_x[k] * isi_752[k];

        t_961[k] = f_3 * pc_x[k] * isi_753[k];

        t_962[k] = f_3 * pc_x[k] * isi_754[k];

        t_963[k] = f_3 * pc_x[k] * isi_755[k];

        t_964[k] = pa_y[k] * hsk0_748[k]
                   + f_20 * hsi_581[k]
                   - f_12 * pc_y[k] * hsk1_748[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, pa_y, pc_y, pc_z, hsk0_750, hsk0_751, hsi_553, \
                         hsi_583, hsi_584, hsk1_750, hsk1_751, \
                         isi_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_17 * hsi_553[k]
                   + f_3 * pc_z[k] * isi_749[k];

        t_966[k] = pa_y[k] * hsk0_750[k]
                   + f_17 * hsi_583[k]
                   - f_12 * pc_y[k] * hsk1_750[k];

        t_967[k] = pa_y[k] * hsk0_751[k]
                   + f_16 * hsi_584[k]
                   - f_12 * pc_y[k] * hsk1_751[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, pa_y, pc_y, hsk0_752, hsk0_753, hsk0_755, \
                         hsi_585, hsi_586, hsi_587, hsk1_752, hsk1_753, hsk1_755, \
                         isi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = pa_y[k] * hsk0_752[k]
                   + f_15 * hsi_585[k]
                   - f_12 * pc_y[k] * hsk1_752[k];

        t_969[k] = pa_y[k] * hsk0_753[k]
                   + f_14 * hsi_586[k]
                   - f_12 * pc_y[k] * hsk1_753[k];

        t_970[k] = f_13 * hsi_587[k]
                   + f_3 * pc_y[k] * isi_755[k];

        t_971[k] = pa_y[k] * hsk0_755[k]
                   - f_12 * pc_y[k] * hsk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, pc_x, pc_y, ish0_567, ish0_569, \
                         ish0_570, ish1_567, ish1_569, ish1_570, isi_756, isi_758, \
                         isi_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_1 * ish0_567[k]
                   - f_2 * ish1_567[k]
                   + f_3 * pc_x[k] * isi_756[k];

        t_973[k] = f_3 * pc_y[k] * isi_756[k];

        t_974[k] = f_18 * ish0_569[k]
                   - f_19 * ish1_569[k]
                   + f_3 * pc_x[k] * isi_758[k];

        t_975[k] = f_10 * ish0_570[k]
                   - f_11 * ish1_570[k]
                   + f_3 * pc_x[k] * isi_759[k];

        t_976[k] = f_3 * pc_y[k] * isi_758[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, pc_x, pc_y, ish0_572, ish0_573, ish0_574, \
                         ish1_572, ish1_573, ish1_574, isi_761, isi_762, \
                         isi_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_10 * ish0_572[k]
                   - f_11 * ish1_572[k]
                   + f_3 * pc_x[k] * isi_761[k];

        t_978[k] = f_8 * ish0_573[k]
                   - f_9 * ish1_573[k]
                   + f_3 * pc_x[k] * isi_762[k];

        t_979[k] = f_8 * ish0_574[k]
                   - f_9 * ish1_574[k]
                   + f_3 * pc_x[k] * isi_763[k];

        t_980[k] = f_3 * pc_y[k] * isi_761[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_x, ish0_576, ish0_577, ish0_578, ish1_576, \
                         ish1_577, ish1_578, isi_765, isi_766, \
                         isi_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_8 * ish0_576[k]
                   - f_9 * ish1_576[k]
                   + f_3 * pc_x[k] * isi_765[k];

        t_982[k] = f_6 * ish0_577[k]
                   - f_7 * ish1_577[k]
                   + f_3 * pc_x[k] * isi_766[k];

        t_983[k] = f_6 * ish0_578[k]
                   - f_7 * ish1_578[k]
                   + f_3 * pc_x[k] * isi_767[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, pc_x, pc_y, ish0_579, ish0_581, ish0_582, \
                         ish1_579, ish1_581, ish1_582, isi_765, isi_768, isi_770, \
                         isi_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_6 * ish0_579[k]
                   - f_7 * ish1_579[k]
                   + f_3 * pc_x[k] * isi_768[k];

        t_985[k] = f_3 * pc_y[k] * isi_765[k];

        t_986[k] = f_6 * ish0_581[k]
                   - f_7 * ish1_581[k]
                   + f_3 * pc_x[k] * isi_770[k];

        t_987[k] = f_4 * ish0_582[k]
                   - f_5 * ish1_582[k]
                   + f_3 * pc_x[k] * isi_771[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, pc_x, pc_y, ish0_583, ish0_584, ish0_585, \
                         ish1_583, ish1_584, ish1_585, isi_770, isi_772, isi_773, \
                         isi_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_4 * ish0_583[k]
                   - f_5 * ish1_583[k]
                   + f_3 * pc_x[k] * isi_772[k];

        t_989[k] = f_4 * ish0_584[k]
                   - f_5 * ish1_584[k]
                   + f_3 * pc_x[k] * isi_773[k];

        t_990[k] = f_4 * ish0_585[k]
                   - f_5 * ish1_585[k]
                   + f_3 * pc_x[k] * isi_774[k];

        t_991[k] = f_3 * pc_y[k] * isi_770[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, t_997, pc_x, ish0_587, ish1_587, \
                         isi_776, isi_777, isi_778, isi_779, isi_780, \
                         isi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * ish0_587[k]
                   - f_5 * ish1_587[k]
                   + f_3 * pc_x[k] * isi_776[k];

        t_993[k] = f_3 * pc_x[k] * isi_777[k];

        t_994[k] = f_3 * pc_x[k] * isi_778[k];

        t_995[k] = f_3 * pc_x[k] * isi_779[k];

        t_996[k] = f_3 * pc_x[k] * isi_780[k];

        t_997[k] = f_3 * pc_x[k] * isi_781[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, pc_x, pc_y, ish0_582, ish0_583, \
                         ish1_582, ish1_583, isi_777, isi_778, isi_782, \
                         isi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_3 * pc_x[k] * isi_782[k];

        t_999[k] = f_3 * pc_x[k] * isi_783[k];

        t_1000[k] = f_1 * ish0_582[k]
                    - f_2 * ish1_582[k]
                    + f_3 * pc_y[k] * isi_777[k];

        t_1001[k] = f_18 * ish0_583[k]
                    - f_19 * ish1_583[k]
                    + f_3 * pc_y[k] * isi_778[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, ish0_584, ish0_585, ish0_586, ish1_584, \
                         ish1_585, ish1_586, isi_779, isi_780, \
                         isi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_10 * ish0_584[k]
                    - f_11 * ish1_584[k]
                    + f_3 * pc_y[k] * isi_779[k];

        t_1003[k] = f_8 * ish0_585[k]
                    - f_9 * ish1_585[k]
                    + f_3 * pc_y[k] * isi_780[k];

        t_1004[k] = f_6 * ish0_586[k]
                    - f_7 * ish1_586[k]
                    + f_3 * pc_y[k] * isi_781[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pc_y, pc_z, hsi_587, ish0_587, ish1_587, \
                         isi_782, isi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_4 * ish0_587[k]
                    - f_5 * ish1_587[k]
                    + f_3 * pc_y[k] * isi_782[k];

        t_1006[k] = f_3 * pc_y[k] * isi_783[k];

        t_1007[k] = f_0 * hsi_587[k]
                    + f_1 * ish0_587[k]
                    - f_2 * ish1_587[k]
                    + f_3 * pc_z[k] * isi_783[k];
    }
}

auto
compute_prim_isk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsk0, const size_t hsi,
                                                   const size_t hsk1, const size_t ish0,
                                                   const size_t ish1, const size_t isi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_isk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, isi, ncols, gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);

    compute_prim_isk_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, hsk0, hsi,
                                                              hsk1, ish0, ish1, isi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
