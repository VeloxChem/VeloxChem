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


#include "SimdThreeCenterElectronRepulsionVrrRecDSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dsl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t psl0,
                                                          const size_t psk, const size_t psl1,
                                                          const size_t dsi0, const size_t dsi1,
                                                          const size_t dsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 1.5 / q;

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
    auto *t_129 = buffer.data(target + 129);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psl0_0 = buffer.data(psl0 + 0);
    const auto *psl0_3 = buffer.data(psl0 + 3);
    const auto *psl0_5 = buffer.data(psl0 + 5);
    const auto *psl0_6 = buffer.data(psl0 + 6);
    const auto *psl0_9 = buffer.data(psl0 + 9);
    const auto *psl0_10 = buffer.data(psl0 + 10);
    const auto *psl0_14 = buffer.data(psl0 + 14);
    const auto *psl0_15 = buffer.data(psl0 + 15);
    const auto *psl0_20 = buffer.data(psl0 + 20);
    const auto *psl0_21 = buffer.data(psl0 + 21);
    const auto *psl0_27 = buffer.data(psl0 + 27);
    const auto *psl0_48 = buffer.data(psl0 + 48);
    const auto *psl0_51 = buffer.data(psl0 + 51);
    const auto *psl0_55 = buffer.data(psl0 + 55);
    const auto *psl0_60 = buffer.data(psl0 + 60);
    const auto *psl0_66 = buffer.data(psl0 + 66);
    const auto *psl0_81 = buffer.data(psl0 + 81);
    const auto *psl0_83 = buffer.data(psl0 + 83);
    const auto *psl0_84 = buffer.data(psl0 + 84);
    const auto *psl0_85 = buffer.data(psl0 + 85);
    const auto *psl0_86 = buffer.data(psl0 + 86);
    const auto *psl0_87 = buffer.data(psl0 + 87);
    const auto *psl0_89 = buffer.data(psl0 + 89);
    const auto *psl0_95 = buffer.data(psl0 + 95);
    const auto *psl0_99 = buffer.data(psl0 + 99);
    const auto *psl0_104 = buffer.data(psl0 + 104);
    const auto *psl0_110 = buffer.data(psl0 + 110);
    const auto *psl0_117 = buffer.data(psl0 + 117);
    const auto *psl0_126 = buffer.data(psl0 + 126);
    const auto *psl0_127 = buffer.data(psl0 + 127);
    const auto *psl0_128 = buffer.data(psl0 + 128);
    const auto *psl0_129 = buffer.data(psl0 + 129);

    const auto *psk_0 = buffer.data(psk + 0);
    const auto *psk_5 = buffer.data(psk + 5);
    const auto *psk_9 = buffer.data(psk + 9);
    const auto *psk_14 = buffer.data(psk + 14);
    const auto *psk_20 = buffer.data(psk + 20);
    const auto *psk_28 = buffer.data(psk + 28);
    const auto *psk_30 = buffer.data(psk + 30);
    const auto *psk_31 = buffer.data(psk + 31);
    const auto *psk_32 = buffer.data(psk + 32);
    const auto *psk_33 = buffer.data(psk + 33);
    const auto *psk_35 = buffer.data(psk + 35);
    const auto *psk_39 = buffer.data(psk + 39);
    const auto *psk_42 = buffer.data(psk + 42);
    const auto *psk_46 = buffer.data(psk + 46);
    const auto *psk_51 = buffer.data(psk + 51);
    const auto *psk_57 = buffer.data(psk + 57);
    const auto *psk_64 = buffer.data(psk + 64);
    const auto *psk_66 = buffer.data(psk + 66);
    const auto *psk_67 = buffer.data(psk + 67);
    const auto *psk_68 = buffer.data(psk + 68);
    const auto *psk_69 = buffer.data(psk + 69);
    const auto *psk_70 = buffer.data(psk + 70);
    const auto *psk_71 = buffer.data(psk + 71);
    const auto *psk_77 = buffer.data(psk + 77);
    const auto *psk_81 = buffer.data(psk + 81);
    const auto *psk_86 = buffer.data(psk + 86);
    const auto *psk_92 = buffer.data(psk + 92);
    const auto *psk_99 = buffer.data(psk + 99);
    const auto *psk_100 = buffer.data(psk + 100);
    const auto *psk_101 = buffer.data(psk + 101);
    const auto *psk_102 = buffer.data(psk + 102);
    const auto *psk_103 = buffer.data(psk + 103);
    const auto *psk_104 = buffer.data(psk + 104);
    const auto *psk_105 = buffer.data(psk + 105);
    const auto *psk_107 = buffer.data(psk + 107);

    const auto *psl1_0 = buffer.data(psl1 + 0);
    const auto *psl1_3 = buffer.data(psl1 + 3);
    const auto *psl1_5 = buffer.data(psl1 + 5);
    const auto *psl1_6 = buffer.data(psl1 + 6);
    const auto *psl1_9 = buffer.data(psl1 + 9);
    const auto *psl1_10 = buffer.data(psl1 + 10);
    const auto *psl1_14 = buffer.data(psl1 + 14);
    const auto *psl1_15 = buffer.data(psl1 + 15);
    const auto *psl1_20 = buffer.data(psl1 + 20);
    const auto *psl1_21 = buffer.data(psl1 + 21);
    const auto *psl1_27 = buffer.data(psl1 + 27);
    const auto *psl1_48 = buffer.data(psl1 + 48);
    const auto *psl1_51 = buffer.data(psl1 + 51);
    const auto *psl1_55 = buffer.data(psl1 + 55);
    const auto *psl1_60 = buffer.data(psl1 + 60);
    const auto *psl1_66 = buffer.data(psl1 + 66);
    const auto *psl1_81 = buffer.data(psl1 + 81);
    const auto *psl1_83 = buffer.data(psl1 + 83);
    const auto *psl1_84 = buffer.data(psl1 + 84);
    const auto *psl1_85 = buffer.data(psl1 + 85);
    const auto *psl1_86 = buffer.data(psl1 + 86);
    const auto *psl1_87 = buffer.data(psl1 + 87);
    const auto *psl1_89 = buffer.data(psl1 + 89);
    const auto *psl1_95 = buffer.data(psl1 + 95);
    const auto *psl1_99 = buffer.data(psl1 + 99);
    const auto *psl1_104 = buffer.data(psl1 + 104);
    const auto *psl1_110 = buffer.data(psl1 + 110);
    const auto *psl1_117 = buffer.data(psl1 + 117);
    const auto *psl1_126 = buffer.data(psl1 + 126);
    const auto *psl1_127 = buffer.data(psl1 + 127);
    const auto *psl1_128 = buffer.data(psl1 + 128);
    const auto *psl1_129 = buffer.data(psl1 + 129);

    const auto *dsi0_0 = buffer.data(dsi0 + 0);
    const auto *dsi0_1 = buffer.data(dsi0 + 1);
    const auto *dsi0_2 = buffer.data(dsi0 + 2);
    const auto *dsi0_3 = buffer.data(dsi0 + 3);
    const auto *dsi0_5 = buffer.data(dsi0 + 5);
    const auto *dsi0_6 = buffer.data(dsi0 + 6);
    const auto *dsi0_8 = buffer.data(dsi0 + 8);
    const auto *dsi0_9 = buffer.data(dsi0 + 9);
    const auto *dsi0_10 = buffer.data(dsi0 + 10);
    const auto *dsi0_12 = buffer.data(dsi0 + 12);
    const auto *dsi0_13 = buffer.data(dsi0 + 13);
    const auto *dsi0_14 = buffer.data(dsi0 + 14);
    const auto *dsi0_21 = buffer.data(dsi0 + 21);
    const auto *dsi0_23 = buffer.data(dsi0 + 23);
    const auto *dsi0_24 = buffer.data(dsi0 + 24);
    const auto *dsi0_25 = buffer.data(dsi0 + 25);
    const auto *dsi0_26 = buffer.data(dsi0 + 26);
    const auto *dsi0_27 = buffer.data(dsi0 + 27);
    const auto *dsi0_31 = buffer.data(dsi0 + 31);
    const auto *dsi0_34 = buffer.data(dsi0 + 34);
    const auto *dsi0_35 = buffer.data(dsi0 + 35);
    const auto *dsi0_38 = buffer.data(dsi0 + 38);
    const auto *dsi0_39 = buffer.data(dsi0 + 39);
    const auto *dsi0_40 = buffer.data(dsi0 + 40);
    const auto *dsi0_58 = buffer.data(dsi0 + 58);
    const auto *dsi0_60 = buffer.data(dsi0 + 60);
    const auto *dsi0_61 = buffer.data(dsi0 + 61);
    const auto *dsi0_63 = buffer.data(dsi0 + 63);
    const auto *dsi0_64 = buffer.data(dsi0 + 64);
    const auto *dsi0_65 = buffer.data(dsi0 + 65);
    const auto *dsi0_67 = buffer.data(dsi0 + 67);
    const auto *dsi0_68 = buffer.data(dsi0 + 68);
    const auto *dsi0_69 = buffer.data(dsi0 + 69);
    const auto *dsi0_70 = buffer.data(dsi0 + 70);

    const auto *dsi1_0 = buffer.data(dsi1 + 0);
    const auto *dsi1_1 = buffer.data(dsi1 + 1);
    const auto *dsi1_2 = buffer.data(dsi1 + 2);
    const auto *dsi1_3 = buffer.data(dsi1 + 3);
    const auto *dsi1_5 = buffer.data(dsi1 + 5);
    const auto *dsi1_6 = buffer.data(dsi1 + 6);
    const auto *dsi1_8 = buffer.data(dsi1 + 8);
    const auto *dsi1_9 = buffer.data(dsi1 + 9);
    const auto *dsi1_10 = buffer.data(dsi1 + 10);
    const auto *dsi1_12 = buffer.data(dsi1 + 12);
    const auto *dsi1_13 = buffer.data(dsi1 + 13);
    const auto *dsi1_14 = buffer.data(dsi1 + 14);
    const auto *dsi1_21 = buffer.data(dsi1 + 21);
    const auto *dsi1_23 = buffer.data(dsi1 + 23);
    const auto *dsi1_24 = buffer.data(dsi1 + 24);
    const auto *dsi1_25 = buffer.data(dsi1 + 25);
    const auto *dsi1_26 = buffer.data(dsi1 + 26);
    const auto *dsi1_27 = buffer.data(dsi1 + 27);
    const auto *dsi1_31 = buffer.data(dsi1 + 31);
    const auto *dsi1_34 = buffer.data(dsi1 + 34);
    const auto *dsi1_35 = buffer.data(dsi1 + 35);
    const auto *dsi1_38 = buffer.data(dsi1 + 38);
    const auto *dsi1_39 = buffer.data(dsi1 + 39);
    const auto *dsi1_40 = buffer.data(dsi1 + 40);
    const auto *dsi1_58 = buffer.data(dsi1 + 58);
    const auto *dsi1_60 = buffer.data(dsi1 + 60);
    const auto *dsi1_61 = buffer.data(dsi1 + 61);
    const auto *dsi1_63 = buffer.data(dsi1 + 63);
    const auto *dsi1_64 = buffer.data(dsi1 + 64);
    const auto *dsi1_65 = buffer.data(dsi1 + 65);
    const auto *dsi1_67 = buffer.data(dsi1 + 67);
    const auto *dsi1_68 = buffer.data(dsi1 + 68);
    const auto *dsi1_69 = buffer.data(dsi1 + 69);
    const auto *dsi1_70 = buffer.data(dsi1 + 70);

    const auto *dsk_0 = buffer.data(dsk + 0);
    const auto *dsk_1 = buffer.data(dsk + 1);
    const auto *dsk_2 = buffer.data(dsk + 2);
    const auto *dsk_3 = buffer.data(dsk + 3);
    const auto *dsk_5 = buffer.data(dsk + 5);
    const auto *dsk_6 = buffer.data(dsk + 6);
    const auto *dsk_8 = buffer.data(dsk + 8);
    const auto *dsk_9 = buffer.data(dsk + 9);
    const auto *dsk_10 = buffer.data(dsk + 10);
    const auto *dsk_12 = buffer.data(dsk + 12);
    const auto *dsk_13 = buffer.data(dsk + 13);
    const auto *dsk_14 = buffer.data(dsk + 14);
    const auto *dsk_15 = buffer.data(dsk + 15);
    const auto *dsk_17 = buffer.data(dsk + 17);
    const auto *dsk_18 = buffer.data(dsk + 18);
    const auto *dsk_19 = buffer.data(dsk + 19);
    const auto *dsk_20 = buffer.data(dsk + 20);
    const auto *dsk_21 = buffer.data(dsk + 21);
    const auto *dsk_27 = buffer.data(dsk + 27);
    const auto *dsk_28 = buffer.data(dsk + 28);
    const auto *dsk_30 = buffer.data(dsk + 30);
    const auto *dsk_31 = buffer.data(dsk + 31);
    const auto *dsk_32 = buffer.data(dsk + 32);
    const auto *dsk_33 = buffer.data(dsk + 33);
    const auto *dsk_34 = buffer.data(dsk + 34);
    const auto *dsk_35 = buffer.data(dsk + 35);
    const auto *dsk_36 = buffer.data(dsk + 36);
    const auto *dsk_37 = buffer.data(dsk + 37);
    const auto *dsk_39 = buffer.data(dsk + 39);
    const auto *dsk_41 = buffer.data(dsk + 41);
    const auto *dsk_42 = buffer.data(dsk + 42);
    const auto *dsk_43 = buffer.data(dsk + 43);
    const auto *dsk_45 = buffer.data(dsk + 45);
    const auto *dsk_46 = buffer.data(dsk + 46);
    const auto *dsk_47 = buffer.data(dsk + 47);
    const auto *dsk_48 = buffer.data(dsk + 48);
    const auto *dsk_50 = buffer.data(dsk + 50);
    const auto *dsk_51 = buffer.data(dsk + 51);
    const auto *dsk_52 = buffer.data(dsk + 52);
    const auto *dsk_53 = buffer.data(dsk + 53);
    const auto *dsk_54 = buffer.data(dsk + 54);
    const auto *dsk_56 = buffer.data(dsk + 56);
    const auto *dsk_57 = buffer.data(dsk + 57);
    const auto *dsk_64 = buffer.data(dsk + 64);
    const auto *dsk_66 = buffer.data(dsk + 66);
    const auto *dsk_67 = buffer.data(dsk + 67);
    const auto *dsk_68 = buffer.data(dsk + 68);
    const auto *dsk_69 = buffer.data(dsk + 69);
    const auto *dsk_70 = buffer.data(dsk + 70);
    const auto *dsk_71 = buffer.data(dsk + 71);
    const auto *dsk_72 = buffer.data(dsk + 72);
    const auto *dsk_74 = buffer.data(dsk + 74);
    const auto *dsk_76 = buffer.data(dsk + 76);
    const auto *dsk_77 = buffer.data(dsk + 77);
    const auto *dsk_79 = buffer.data(dsk + 79);
    const auto *dsk_80 = buffer.data(dsk + 80);
    const auto *dsk_81 = buffer.data(dsk + 81);
    const auto *dsk_83 = buffer.data(dsk + 83);
    const auto *dsk_84 = buffer.data(dsk + 84);
    const auto *dsk_85 = buffer.data(dsk + 85);
    const auto *dsk_86 = buffer.data(dsk + 86);
    const auto *dsk_88 = buffer.data(dsk + 88);
    const auto *dsk_89 = buffer.data(dsk + 89);
    const auto *dsk_90 = buffer.data(dsk + 90);
    const auto *dsk_91 = buffer.data(dsk + 91);
    const auto *dsk_92 = buffer.data(dsk + 92);
    const auto *dsk_99 = buffer.data(dsk + 99);
    const auto *dsk_100 = buffer.data(dsk + 100);
    const auto *dsk_101 = buffer.data(dsk + 101);
    const auto *dsk_102 = buffer.data(dsk + 102);
    const auto *dsk_103 = buffer.data(dsk + 103);
    const auto *dsk_104 = buffer.data(dsk + 104);
    const auto *dsk_105 = buffer.data(dsk + 105);
    const auto *dsk_107 = buffer.data(dsk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, psk_0, dsi0_0, \
                         dsi1_0, dsk_0, dsk_1, dsk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psk_0[k]
                 + f_1 * dsi0_0[k]
                 - f_2 * dsi1_0[k]
                 + f_3 * pc_x[k] * dsk_0[k];

        t_1[k] = f_3 * pc_y[k] * dsk_0[k];

        t_2[k] = f_3 * pc_z[k] * dsk_0[k];

        t_3[k] = f_4 * dsi0_0[k]
                 - f_5 * dsi1_0[k]
                 + f_3 * pc_y[k] * dsk_1[k];

        t_4[k] = f_3 * pc_y[k] * dsk_2[k];

        t_5[k] = f_4 * dsi0_0[k]
                 - f_5 * dsi1_0[k]
                 + f_3 * pc_z[k] * dsk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, dsi0_1, dsi0_2, dsi0_3, dsi1_1, \
                         dsi1_2, dsi1_3, dsk_3, dsk_5, dsk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * dsi0_1[k]
                 - f_7 * dsi1_1[k]
                 + f_3 * pc_y[k] * dsk_3[k];

        t_7[k] = f_3 * pc_z[k] * dsk_3[k];

        t_8[k] = f_3 * pc_y[k] * dsk_5[k];

        t_9[k] = f_6 * dsi0_2[k]
                 - f_7 * dsi1_2[k]
                 + f_3 * pc_z[k] * dsk_5[k];

        t_10[k] = f_8 * dsi0_3[k]
                  - f_9 * dsi1_3[k]
                  + f_3 * pc_y[k] * dsk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, dsi0_5, dsi0_6, \
                         dsi1_5, dsi1_6, dsk_6, dsk_8, dsk_9, dsk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * dsk_6[k];

        t_12[k] = f_4 * dsi0_5[k]
                  - f_5 * dsi1_5[k]
                  + f_3 * pc_y[k] * dsk_8[k];

        t_13[k] = f_3 * pc_y[k] * dsk_9[k];

        t_14[k] = f_8 * dsi0_5[k]
                  - f_9 * dsi1_5[k]
                  + f_3 * pc_z[k] * dsk_9[k];

        t_15[k] = f_10 * dsi0_6[k]
                  - f_11 * dsi1_6[k]
                  + f_3 * pc_y[k] * dsk_10[k];

        t_16[k] = f_3 * pc_z[k] * dsk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, dsi0_8, dsi0_9, dsi1_8, dsi1_9, \
                         dsk_12, dsk_13, dsk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * dsi0_8[k]
                  - f_7 * dsi1_8[k]
                  + f_3 * pc_y[k] * dsk_12[k];

        t_18[k] = f_4 * dsi0_9[k]
                  - f_5 * dsi1_9[k]
                  + f_3 * pc_y[k] * dsk_13[k];

        t_19[k] = f_3 * pc_y[k] * dsk_14[k];

        t_20[k] = f_10 * dsi0_9[k]
                  - f_11 * dsi1_9[k]
                  + f_3 * pc_z[k] * dsk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, dsi0_10, dsi0_12, dsi0_13, \
                         dsi1_10, dsi1_12, dsi1_13, dsk_15, dsk_17, \
                         dsk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * dsi0_10[k]
                  - f_13 * dsi1_10[k]
                  + f_3 * pc_y[k] * dsk_15[k];

        t_22[k] = f_3 * pc_z[k] * dsk_15[k];

        t_23[k] = f_8 * dsi0_12[k]
                  - f_9 * dsi1_12[k]
                  + f_3 * pc_y[k] * dsk_17[k];

        t_24[k] = f_6 * dsi0_13[k]
                  - f_7 * dsi1_13[k]
                  + f_3 * pc_y[k] * dsk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, psk_28, dsi0_14, \
                         dsi1_14, dsk_19, dsk_20, dsk_21, dsk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * dsi0_14[k]
                  - f_5 * dsi1_14[k]
                  + f_3 * pc_y[k] * dsk_19[k];

        t_26[k] = f_3 * pc_y[k] * dsk_20[k];

        t_27[k] = f_12 * dsi0_14[k]
                  - f_13 * dsi1_14[k]
                  + f_3 * pc_z[k] * dsk_20[k];

        t_28[k] = f_0 * psk_28[k]
                  + f_3 * pc_x[k] * dsk_28[k];

        t_29[k] = f_3 * pc_z[k] * dsk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, psk_30, psk_31, psk_32, \
                         psk_33, dsk_27, dsk_30, dsk_31, dsk_32, \
                         dsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * psk_30[k]
                  + f_3 * pc_x[k] * dsk_30[k];

        t_31[k] = f_0 * psk_31[k]
                  + f_3 * pc_x[k] * dsk_31[k];

        t_32[k] = f_0 * psk_32[k]
                  + f_3 * pc_x[k] * dsk_32[k];

        t_33[k] = f_0 * psk_33[k]
                  + f_3 * pc_x[k] * dsk_33[k];

        t_34[k] = f_3 * pc_y[k] * dsk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, psk_35, dsi0_21, dsi0_23, \
                         dsi1_21, dsi1_23, dsk_28, dsk_30, dsk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * psk_35[k]
                  + f_3 * pc_x[k] * dsk_35[k];

        t_36[k] = f_1 * dsi0_21[k]
                  - f_2 * dsi1_21[k]
                  + f_3 * pc_y[k] * dsk_28[k];

        t_37[k] = f_3 * pc_z[k] * dsk_28[k];

        t_38[k] = f_12 * dsi0_23[k]
                  - f_13 * dsi1_23[k]
                  + f_3 * pc_y[k] * dsk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, dsi0_24, dsi0_25, dsi0_26, dsi1_24, dsi1_25, \
                         dsi1_26, dsk_31, dsk_32, dsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * dsi0_24[k]
                  - f_11 * dsi1_24[k]
                  + f_3 * pc_y[k] * dsk_31[k];

        t_40[k] = f_8 * dsi0_25[k]
                  - f_9 * dsi1_25[k]
                  + f_3 * pc_y[k] * dsk_32[k];

        t_41[k] = f_6 * dsi0_26[k]
                  - f_7 * dsi1_26[k]
                  + f_3 * pc_y[k] * dsk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, psl0_0, psk_0, \
                         psl1_0, dsi0_27, dsi1_27, dsk_34, dsk_35, \
                         dsk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * dsi0_27[k]
                  - f_5 * dsi1_27[k]
                  + f_3 * pc_y[k] * dsk_34[k];

        t_43[k] = f_3 * pc_y[k] * dsk_35[k];

        t_44[k] = f_1 * dsi0_27[k]
                  - f_2 * dsi1_27[k]
                  + f_3 * pc_z[k] * dsk_35[k];

        t_45[k] = pa_y[k] * psl0_0[k]
                  - f_14 * pc_y[k] * psl1_0[k];

        t_46[k] = f_15 * psk_0[k]
                  + f_3 * pc_y[k] * dsk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pa_y, pc_x, pc_y, pc_z, psl0_5, \
                         psl0_48, psk_39, psl1_5, psl1_48, dsk_36, \
                         dsk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * dsk_36[k];

        t_48[k] = pa_x[k] * psl0_48[k]
                  + f_16 * psk_39[k]
                  - f_14 * pc_x[k] * psl1_48[k];

        t_49[k] = f_3 * pc_z[k] * dsk_37[k];

        t_50[k] = pa_y[k] * psl0_5[k]
                  - f_14 * pc_y[k] * psl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_x, pc_x, pc_y, pc_z, psl0_51, psk_5, psk_42, \
                         psl1_51, dsk_39, dsk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_x[k] * psl0_51[k]
                  + f_17 * psk_42[k]
                  - f_14 * pc_x[k] * psl1_51[k];

        t_52[k] = f_3 * pc_z[k] * dsk_39[k];

        t_53[k] = f_15 * psk_5[k]
                  + f_3 * pc_y[k] * dsk_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_x, pa_y, pc_x, pc_y, pc_z, psl0_9, psl0_55, \
                         psk_46, psl1_9, psl1_55, dsk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_y[k] * psl0_9[k]
                  - f_14 * pc_y[k] * psl1_9[k];

        t_55[k] = pa_x[k] * psl0_55[k]
                  + f_18 * psk_46[k]
                  - f_14 * pc_x[k] * psl1_55[k];

        t_56[k] = f_3 * pc_z[k] * dsk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pc_y, pc_z, psl0_14, psk_9, psl1_14, dsi0_31, \
                         dsi1_31, dsk_43, dsk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_4 * dsi0_31[k]
                  - f_5 * dsi1_31[k]
                  + f_3 * pc_z[k] * dsk_43[k];

        t_58[k] = f_15 * psk_9[k]
                  + f_3 * pc_y[k] * dsk_45[k];

        t_59[k] = pa_y[k] * psl0_14[k]
                  - f_14 * pc_y[k] * psl1_14[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pc_x, pc_z, psl0_60, psk_51, psl1_60, \
                         dsi0_34, dsi1_34, dsk_46, dsk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * psl0_60[k]
                  + f_19 * psk_51[k]
                  - f_14 * pc_x[k] * psl1_60[k];

        t_61[k] = f_3 * pc_z[k] * dsk_46[k];

        t_62[k] = f_4 * dsi0_34[k]
                  - f_5 * dsi1_34[k]
                  + f_3 * pc_z[k] * dsk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, psl0_20, psk_14, psl1_20, \
                         dsi0_35, dsi1_35, dsk_48, dsk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * dsi0_35[k]
                  - f_7 * dsi1_35[k]
                  + f_3 * pc_z[k] * dsk_48[k];

        t_64[k] = f_15 * psk_14[k]
                  + f_3 * pc_y[k] * dsk_50[k];

        t_65[k] = pa_y[k] * psl0_20[k]
                  - f_14 * pc_y[k] * psl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_x, pc_x, pc_z, psl0_66, psk_57, psl1_66, \
                         dsi0_38, dsi1_38, dsk_51, dsk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * psl0_66[k]
                  + f_0 * psk_57[k]
                  - f_14 * pc_x[k] * psl1_66[k];

        t_67[k] = f_3 * pc_z[k] * dsk_51[k];

        t_68[k] = f_4 * dsi0_38[k]
                  - f_5 * dsi1_38[k]
                  + f_3 * pc_z[k] * dsk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, psk_20, dsi0_39, dsi0_40, dsi1_39, \
                         dsi1_40, dsk_53, dsk_54, dsk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * dsi0_39[k]
                  - f_7 * dsi1_39[k]
                  + f_3 * pc_z[k] * dsk_53[k];

        t_70[k] = f_8 * dsi0_40[k]
                  - f_9 * dsi1_40[k]
                  + f_3 * pc_z[k] * dsk_54[k];

        t_71[k] = f_15 * psk_20[k]
                  + f_3 * pc_y[k] * dsk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, psl0_27, psk_64, \
                         psk_66, psl1_27, dsk_57, dsk_64, dsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * psl0_27[k]
                  - f_14 * pc_y[k] * psl1_27[k];

        t_73[k] = f_15 * psk_64[k]
                  + f_3 * pc_x[k] * dsk_64[k];

        t_74[k] = f_3 * pc_z[k] * dsk_57[k];

        t_75[k] = f_15 * psk_66[k]
                  + f_3 * pc_x[k] * dsk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, psk_67, psk_68, psk_69, psk_70, \
                         psk_71, dsk_67, dsk_68, dsk_69, dsk_70, \
                         dsk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_15 * psk_67[k]
                  + f_3 * pc_x[k] * dsk_67[k];

        t_77[k] = f_15 * psk_68[k]
                  + f_3 * pc_x[k] * dsk_68[k];

        t_78[k] = f_15 * psk_69[k]
                  + f_3 * pc_x[k] * dsk_69[k];

        t_79[k] = f_15 * psk_70[k]
                  + f_3 * pc_x[k] * dsk_70[k];

        t_80[k] = f_15 * psk_71[k]
                  + f_3 * pc_x[k] * dsk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_x, pc_x, pc_z, psl0_81, psl0_83, psl0_84, \
                         psl1_81, psl1_83, psl1_84, dsk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_x[k] * psl0_81[k]
                  - f_14 * pc_x[k] * psl1_81[k];

        t_82[k] = f_3 * pc_z[k] * dsk_64[k];

        t_83[k] = pa_x[k] * psl0_83[k]
                  - f_14 * pc_x[k] * psl1_83[k];

        t_84[k] = pa_x[k] * psl0_84[k]
                  - f_14 * pc_x[k] * psl1_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_x, pc_x, pc_y, psl0_85, psl0_86, psl0_87, \
                         psk_35, psl1_85, psl1_86, psl1_87, dsk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_x[k] * psl0_85[k]
                  - f_14 * pc_x[k] * psl1_85[k];

        t_86[k] = pa_x[k] * psl0_86[k]
                  - f_14 * pc_x[k] * psl1_86[k];

        t_87[k] = pa_x[k] * psl0_87[k]
                  - f_14 * pc_x[k] * psl1_87[k];

        t_88[k] = f_15 * psk_35[k]
                  + f_3 * pc_y[k] * dsk_71[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pa_z, pc_x, pc_y, pc_z, psl0_0, \
                         psl0_89, psk_0, psl1_0, psl1_89, dsk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_x[k] * psl0_89[k]
                  - f_14 * pc_x[k] * psl1_89[k];

        t_90[k] = pa_z[k] * psl0_0[k]
                  - f_14 * pc_z[k] * psl1_0[k];

        t_91[k] = f_3 * pc_y[k] * dsk_72[k];

        t_92[k] = f_15 * psk_0[k]
                  + f_3 * pc_z[k] * dsk_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_x, pa_z, pc_x, pc_y, pc_z, psl0_3, psl0_95, \
                         psk_77, psl1_3, psl1_95, dsk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * psl0_3[k]
                  - f_14 * pc_z[k] * psl1_3[k];

        t_94[k] = f_3 * pc_y[k] * dsk_74[k];

        t_95[k] = pa_x[k] * psl0_95[k]
                  + f_16 * psk_77[k]
                  - f_14 * pc_x[k] * psl1_95[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_z, pc_y, pc_z, psl0_6, psl1_6, dsi0_58, dsi1_58, \
                         dsk_76, dsk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * psl0_6[k]
                  - f_14 * pc_z[k] * psl1_6[k];

        t_97[k] = f_4 * dsi0_58[k]
                  - f_5 * dsi1_58[k]
                  + f_3 * pc_y[k] * dsk_76[k];

        t_98[k] = f_3 * pc_y[k] * dsk_77[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_x, pa_z, pc_x, pc_y, pc_z, psl0_10, psl0_99, \
                         psk_81, psl1_10, psl1_99, dsi0_60, dsi1_60, \
                         dsk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * psl0_99[k]
                  + f_17 * psk_81[k]
                  - f_14 * pc_x[k] * psl1_99[k];

        t_100[k] = pa_z[k] * psl0_10[k]
                   - f_14 * pc_z[k] * psl1_10[k];

        t_101[k] = f_6 * dsi0_60[k]
                   - f_7 * dsi1_60[k]
                   + f_3 * pc_y[k] * dsk_79[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_x, pc_x, pc_y, psl0_104, psk_86, psl1_104, \
                         dsi0_61, dsi1_61, dsk_80, dsk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * dsi0_61[k]
                   - f_5 * dsi1_61[k]
                   + f_3 * pc_y[k] * dsk_80[k];

        t_103[k] = f_3 * pc_y[k] * dsk_81[k];

        t_104[k] = pa_x[k] * psl0_104[k]
                   + f_18 * psk_86[k]
                   - f_14 * pc_x[k] * psl1_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_z, pc_y, pc_z, psl0_15, psl1_15, dsi0_63, \
                         dsi0_64, dsi1_63, dsi1_64, dsk_83, dsk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * psl0_15[k]
                   - f_14 * pc_z[k] * psl1_15[k];

        t_106[k] = f_8 * dsi0_63[k]
                   - f_9 * dsi1_63[k]
                   + f_3 * pc_y[k] * dsk_83[k];

        t_107[k] = f_6 * dsi0_64[k]
                   - f_7 * dsi1_64[k]
                   + f_3 * pc_y[k] * dsk_84[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_x, pc_x, pc_y, psl0_110, psk_92, psl1_110, \
                         dsi0_65, dsi1_65, dsk_85, dsk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_4 * dsi0_65[k]
                   - f_5 * dsi1_65[k]
                   + f_3 * pc_y[k] * dsk_85[k];

        t_109[k] = f_3 * pc_y[k] * dsk_86[k];

        t_110[k] = pa_x[k] * psl0_110[k]
                   + f_19 * psk_92[k]
                   - f_14 * pc_x[k] * psl1_110[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_z, pc_y, pc_z, psl0_21, psl1_21, dsi0_67, \
                         dsi0_68, dsi1_67, dsi1_68, dsk_88, dsk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_z[k] * psl0_21[k]
                   - f_14 * pc_z[k] * psl1_21[k];

        t_112[k] = f_10 * dsi0_67[k]
                   - f_11 * dsi1_67[k]
                   + f_3 * pc_y[k] * dsk_88[k];

        t_113[k] = f_8 * dsi0_68[k]
                   - f_9 * dsi1_68[k]
                   + f_3 * pc_y[k] * dsk_89[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pc_y, dsi0_69, dsi0_70, dsi1_69, dsi1_70, \
                         dsk_90, dsk_91, dsk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_6 * dsi0_69[k]
                   - f_7 * dsi1_69[k]
                   + f_3 * pc_y[k] * dsk_90[k];

        t_115[k] = f_4 * dsi0_70[k]
                   - f_5 * dsi1_70[k]
                   + f_3 * pc_y[k] * dsk_91[k];

        t_116[k] = f_3 * pc_y[k] * dsk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pc_x, psl0_117, psk_99, psk_100, \
                         psk_101, psk_102, psl1_117, dsk_100, dsk_101, \
                         dsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_x[k] * psl0_117[k]
                   + f_0 * psk_99[k]
                   - f_14 * pc_x[k] * psl1_117[k];

        t_118[k] = f_15 * psk_100[k]
                   + f_3 * pc_x[k] * dsk_100[k];

        t_119[k] = f_15 * psk_101[k]
                   + f_3 * pc_x[k] * dsk_101[k];

        t_120[k] = f_15 * psk_102[k]
                   + f_3 * pc_x[k] * dsk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, psk_103, psk_104, \
                         psk_105, psk_107, dsk_99, dsk_103, dsk_104, dsk_105, \
                         dsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_15 * psk_103[k]
                   + f_3 * pc_x[k] * dsk_103[k];

        t_122[k] = f_15 * psk_104[k]
                   + f_3 * pc_x[k] * dsk_104[k];

        t_123[k] = f_15 * psk_105[k]
                   + f_3 * pc_x[k] * dsk_105[k];

        t_124[k] = f_3 * pc_y[k] * dsk_99[k];

        t_125[k] = f_15 * psk_107[k]
                   + f_3 * pc_x[k] * dsk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_x, pc_x, psl0_126, psl0_127, psl0_128, \
                         psl0_129, psl1_126, psl1_127, psl1_128, \
                         psl1_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_x[k] * psl0_126[k]
                   - f_14 * pc_x[k] * psl1_126[k];

        t_127[k] = pa_x[k] * psl0_127[k]
                   - f_14 * pc_x[k] * psl1_127[k];

        t_128[k] = pa_x[k] * psl0_128[k]
                   - f_14 * pc_x[k] * psl1_128[k];

        t_129[k] = pa_x[k] * psl0_129[k]
                   - f_14 * pc_x[k] * psl1_129[k];
    }
}

static auto
compute_prim_dsl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t psl0,
                                                          const size_t psk, const size_t psl1,
                                                          const size_t dsi0, const size_t dsi1,
                                                          const size_t dsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 1.5 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psl0_46 = buffer.data(psl0 + 46);
    const auto *psl0_48 = buffer.data(psl0 + 48);
    const auto *psl0_51 = buffer.data(psl0 + 51);
    const auto *psl0_55 = buffer.data(psl0 + 55);
    const auto *psl0_60 = buffer.data(psl0 + 60);
    const auto *psl0_66 = buffer.data(psl0 + 66);
    const auto *psl0_81 = buffer.data(psl0 + 81);
    const auto *psl0_90 = buffer.data(psl0 + 90);
    const auto *psl0_92 = buffer.data(psl0 + 92);
    const auto *psl0_95 = buffer.data(psl0 + 95);
    const auto *psl0_99 = buffer.data(psl0 + 99);
    const auto *psl0_104 = buffer.data(psl0 + 104);
    const auto *psl0_110 = buffer.data(psl0 + 110);
    const auto *psl0_117 = buffer.data(psl0 + 117);
    const auto *psl0_128 = buffer.data(psl0 + 128);
    const auto *psl0_129 = buffer.data(psl0 + 129);
    const auto *psl0_130 = buffer.data(psl0 + 130);
    const auto *psl0_131 = buffer.data(psl0 + 131);
    const auto *psl0_132 = buffer.data(psl0 + 132);
    const auto *psl0_134 = buffer.data(psl0 + 134);

    const auto *psk_64 = buffer.data(psk + 64);
    const auto *psk_71 = buffer.data(psk + 71);
    const auto *psk_102 = buffer.data(psk + 102);
    const auto *psk_103 = buffer.data(psk + 103);
    const auto *psk_104 = buffer.data(psk + 104);
    const auto *psk_105 = buffer.data(psk + 105);
    const auto *psk_106 = buffer.data(psk + 106);
    const auto *psk_107 = buffer.data(psk + 107);

    const auto *psl1_46 = buffer.data(psl1 + 46);
    const auto *psl1_48 = buffer.data(psl1 + 48);
    const auto *psl1_51 = buffer.data(psl1 + 51);
    const auto *psl1_55 = buffer.data(psl1 + 55);
    const auto *psl1_60 = buffer.data(psl1 + 60);
    const auto *psl1_66 = buffer.data(psl1 + 66);
    const auto *psl1_81 = buffer.data(psl1 + 81);
    const auto *psl1_90 = buffer.data(psl1 + 90);
    const auto *psl1_92 = buffer.data(psl1 + 92);
    const auto *psl1_95 = buffer.data(psl1 + 95);
    const auto *psl1_99 = buffer.data(psl1 + 99);
    const auto *psl1_104 = buffer.data(psl1 + 104);
    const auto *psl1_110 = buffer.data(psl1 + 110);
    const auto *psl1_117 = buffer.data(psl1 + 117);
    const auto *psl1_128 = buffer.data(psl1 + 128);
    const auto *psl1_129 = buffer.data(psl1 + 129);
    const auto *psl1_130 = buffer.data(psl1 + 130);
    const auto *psl1_131 = buffer.data(psl1 + 131);
    const auto *psl1_132 = buffer.data(psl1 + 132);
    const auto *psl1_134 = buffer.data(psl1 + 134);

    const auto *dsi0_84 = buffer.data(dsi0 + 84);
    const auto *dsi0_85 = buffer.data(dsi0 + 85);
    const auto *dsi0_87 = buffer.data(dsi0 + 87);
    const auto *dsi0_89 = buffer.data(dsi0 + 89);
    const auto *dsi0_90 = buffer.data(dsi0 + 90);
    const auto *dsi0_92 = buffer.data(dsi0 + 92);
    const auto *dsi0_93 = buffer.data(dsi0 + 93);
    const auto *dsi0_94 = buffer.data(dsi0 + 94);
    const auto *dsi0_96 = buffer.data(dsi0 + 96);
    const auto *dsi0_97 = buffer.data(dsi0 + 97);
    const auto *dsi0_98 = buffer.data(dsi0 + 98);
    const auto *dsi0_99 = buffer.data(dsi0 + 99);
    const auto *dsi0_101 = buffer.data(dsi0 + 101);
    const auto *dsi0_102 = buffer.data(dsi0 + 102);
    const auto *dsi0_103 = buffer.data(dsi0 + 103);
    const auto *dsi0_104 = buffer.data(dsi0 + 104);
    const auto *dsi0_105 = buffer.data(dsi0 + 105);
    const auto *dsi0_106 = buffer.data(dsi0 + 106);
    const auto *dsi0_107 = buffer.data(dsi0 + 107);
    const auto *dsi0_108 = buffer.data(dsi0 + 108);
    const auto *dsi0_109 = buffer.data(dsi0 + 109);
    const auto *dsi0_110 = buffer.data(dsi0 + 110);
    const auto *dsi0_111 = buffer.data(dsi0 + 111);
    const auto *dsi0_116 = buffer.data(dsi0 + 116);
    const auto *dsi0_119 = buffer.data(dsi0 + 119);
    const auto *dsi0_120 = buffer.data(dsi0 + 120);
    const auto *dsi0_123 = buffer.data(dsi0 + 123);
    const auto *dsi0_124 = buffer.data(dsi0 + 124);
    const auto *dsi0_125 = buffer.data(dsi0 + 125);
    const auto *dsi0_128 = buffer.data(dsi0 + 128);
    const auto *dsi0_129 = buffer.data(dsi0 + 129);
    const auto *dsi0_130 = buffer.data(dsi0 + 130);
    const auto *dsi0_131 = buffer.data(dsi0 + 131);
    const auto *dsi0_134 = buffer.data(dsi0 + 134);
    const auto *dsi0_135 = buffer.data(dsi0 + 135);
    const auto *dsi0_136 = buffer.data(dsi0 + 136);
    const auto *dsi0_137 = buffer.data(dsi0 + 137);
    const auto *dsi0_138 = buffer.data(dsi0 + 138);
    const auto *dsi0_140 = buffer.data(dsi0 + 140);
    const auto *dsi0_142 = buffer.data(dsi0 + 142);
    const auto *dsi0_143 = buffer.data(dsi0 + 143);
    const auto *dsi0_145 = buffer.data(dsi0 + 145);
    const auto *dsi0_146 = buffer.data(dsi0 + 146);
    const auto *dsi0_147 = buffer.data(dsi0 + 147);
    const auto *dsi0_149 = buffer.data(dsi0 + 149);
    const auto *dsi0_150 = buffer.data(dsi0 + 150);
    const auto *dsi0_151 = buffer.data(dsi0 + 151);
    const auto *dsi0_152 = buffer.data(dsi0 + 152);
    const auto *dsi0_154 = buffer.data(dsi0 + 154);
    const auto *dsi0_155 = buffer.data(dsi0 + 155);
    const auto *dsi0_156 = buffer.data(dsi0 + 156);
    const auto *dsi0_157 = buffer.data(dsi0 + 157);
    const auto *dsi0_158 = buffer.data(dsi0 + 158);
    const auto *dsi0_160 = buffer.data(dsi0 + 160);
    const auto *dsi0_161 = buffer.data(dsi0 + 161);
    const auto *dsi0_162 = buffer.data(dsi0 + 162);
    const auto *dsi0_163 = buffer.data(dsi0 + 163);
    const auto *dsi0_164 = buffer.data(dsi0 + 164);
    const auto *dsi0_165 = buffer.data(dsi0 + 165);
    const auto *dsi0_167 = buffer.data(dsi0 + 167);

    const auto *dsi1_84 = buffer.data(dsi1 + 84);
    const auto *dsi1_85 = buffer.data(dsi1 + 85);
    const auto *dsi1_87 = buffer.data(dsi1 + 87);
    const auto *dsi1_89 = buffer.data(dsi1 + 89);
    const auto *dsi1_90 = buffer.data(dsi1 + 90);
    const auto *dsi1_92 = buffer.data(dsi1 + 92);
    const auto *dsi1_93 = buffer.data(dsi1 + 93);
    const auto *dsi1_94 = buffer.data(dsi1 + 94);
    const auto *dsi1_96 = buffer.data(dsi1 + 96);
    const auto *dsi1_97 = buffer.data(dsi1 + 97);
    const auto *dsi1_98 = buffer.data(dsi1 + 98);
    const auto *dsi1_99 = buffer.data(dsi1 + 99);
    const auto *dsi1_101 = buffer.data(dsi1 + 101);
    const auto *dsi1_102 = buffer.data(dsi1 + 102);
    const auto *dsi1_103 = buffer.data(dsi1 + 103);
    const auto *dsi1_104 = buffer.data(dsi1 + 104);
    const auto *dsi1_105 = buffer.data(dsi1 + 105);
    const auto *dsi1_106 = buffer.data(dsi1 + 106);
    const auto *dsi1_107 = buffer.data(dsi1 + 107);
    const auto *dsi1_108 = buffer.data(dsi1 + 108);
    const auto *dsi1_109 = buffer.data(dsi1 + 109);
    const auto *dsi1_110 = buffer.data(dsi1 + 110);
    const auto *dsi1_111 = buffer.data(dsi1 + 111);
    const auto *dsi1_116 = buffer.data(dsi1 + 116);
    const auto *dsi1_119 = buffer.data(dsi1 + 119);
    const auto *dsi1_120 = buffer.data(dsi1 + 120);
    const auto *dsi1_123 = buffer.data(dsi1 + 123);
    const auto *dsi1_124 = buffer.data(dsi1 + 124);
    const auto *dsi1_125 = buffer.data(dsi1 + 125);
    const auto *dsi1_128 = buffer.data(dsi1 + 128);
    const auto *dsi1_129 = buffer.data(dsi1 + 129);
    const auto *dsi1_130 = buffer.data(dsi1 + 130);
    const auto *dsi1_131 = buffer.data(dsi1 + 131);
    const auto *dsi1_134 = buffer.data(dsi1 + 134);
    const auto *dsi1_135 = buffer.data(dsi1 + 135);
    const auto *dsi1_136 = buffer.data(dsi1 + 136);
    const auto *dsi1_137 = buffer.data(dsi1 + 137);
    const auto *dsi1_138 = buffer.data(dsi1 + 138);
    const auto *dsi1_140 = buffer.data(dsi1 + 140);
    const auto *dsi1_142 = buffer.data(dsi1 + 142);
    const auto *dsi1_143 = buffer.data(dsi1 + 143);
    const auto *dsi1_145 = buffer.data(dsi1 + 145);
    const auto *dsi1_146 = buffer.data(dsi1 + 146);
    const auto *dsi1_147 = buffer.data(dsi1 + 147);
    const auto *dsi1_149 = buffer.data(dsi1 + 149);
    const auto *dsi1_150 = buffer.data(dsi1 + 150);
    const auto *dsi1_151 = buffer.data(dsi1 + 151);
    const auto *dsi1_152 = buffer.data(dsi1 + 152);
    const auto *dsi1_154 = buffer.data(dsi1 + 154);
    const auto *dsi1_155 = buffer.data(dsi1 + 155);
    const auto *dsi1_156 = buffer.data(dsi1 + 156);
    const auto *dsi1_157 = buffer.data(dsi1 + 157);
    const auto *dsi1_158 = buffer.data(dsi1 + 158);
    const auto *dsi1_160 = buffer.data(dsi1 + 160);
    const auto *dsi1_161 = buffer.data(dsi1 + 161);
    const auto *dsi1_162 = buffer.data(dsi1 + 162);
    const auto *dsi1_163 = buffer.data(dsi1 + 163);
    const auto *dsi1_164 = buffer.data(dsi1 + 164);
    const auto *dsi1_165 = buffer.data(dsi1 + 165);
    const auto *dsi1_167 = buffer.data(dsi1 + 167);

    const auto *dsk_107 = buffer.data(dsk + 107);
    const auto *dsk_108 = buffer.data(dsk + 108);
    const auto *dsk_109 = buffer.data(dsk + 109);
    const auto *dsk_111 = buffer.data(dsk + 111);
    const auto *dsk_113 = buffer.data(dsk + 113);
    const auto *dsk_114 = buffer.data(dsk + 114);
    const auto *dsk_116 = buffer.data(dsk + 116);
    const auto *dsk_117 = buffer.data(dsk + 117);
    const auto *dsk_118 = buffer.data(dsk + 118);
    const auto *dsk_120 = buffer.data(dsk + 120);
    const auto *dsk_121 = buffer.data(dsk + 121);
    const auto *dsk_122 = buffer.data(dsk + 122);
    const auto *dsk_123 = buffer.data(dsk + 123);
    const auto *dsk_125 = buffer.data(dsk + 125);
    const auto *dsk_126 = buffer.data(dsk + 126);
    const auto *dsk_127 = buffer.data(dsk + 127);
    const auto *dsk_128 = buffer.data(dsk + 128);
    const auto *dsk_129 = buffer.data(dsk + 129);
    const auto *dsk_131 = buffer.data(dsk + 131);
    const auto *dsk_132 = buffer.data(dsk + 132);
    const auto *dsk_133 = buffer.data(dsk + 133);
    const auto *dsk_134 = buffer.data(dsk + 134);
    const auto *dsk_135 = buffer.data(dsk + 135);
    const auto *dsk_136 = buffer.data(dsk + 136);
    const auto *dsk_137 = buffer.data(dsk + 137);
    const auto *dsk_138 = buffer.data(dsk + 138);
    const auto *dsk_139 = buffer.data(dsk + 139);
    const auto *dsk_140 = buffer.data(dsk + 140);
    const auto *dsk_141 = buffer.data(dsk + 141);
    const auto *dsk_142 = buffer.data(dsk + 142);
    const auto *dsk_143 = buffer.data(dsk + 143);
    const auto *dsk_148 = buffer.data(dsk + 148);
    const auto *dsk_151 = buffer.data(dsk + 151);
    const auto *dsk_152 = buffer.data(dsk + 152);
    const auto *dsk_155 = buffer.data(dsk + 155);
    const auto *dsk_156 = buffer.data(dsk + 156);
    const auto *dsk_157 = buffer.data(dsk + 157);
    const auto *dsk_160 = buffer.data(dsk + 160);
    const auto *dsk_161 = buffer.data(dsk + 161);
    const auto *dsk_162 = buffer.data(dsk + 162);
    const auto *dsk_163 = buffer.data(dsk + 163);
    const auto *dsk_166 = buffer.data(dsk + 166);
    const auto *dsk_167 = buffer.data(dsk + 167);
    const auto *dsk_168 = buffer.data(dsk + 168);
    const auto *dsk_169 = buffer.data(dsk + 169);
    const auto *dsk_170 = buffer.data(dsk + 170);
    const auto *dsk_172 = buffer.data(dsk + 172);
    const auto *dsk_173 = buffer.data(dsk + 173);
    const auto *dsk_174 = buffer.data(dsk + 174);
    const auto *dsk_175 = buffer.data(dsk + 175);
    const auto *dsk_176 = buffer.data(dsk + 176);
    const auto *dsk_177 = buffer.data(dsk + 177);
    const auto *dsk_178 = buffer.data(dsk + 178);
    const auto *dsk_179 = buffer.data(dsk + 179);
    const auto *dsk_180 = buffer.data(dsk + 180);
    const auto *dsk_182 = buffer.data(dsk + 182);
    const auto *dsk_183 = buffer.data(dsk + 183);
    const auto *dsk_185 = buffer.data(dsk + 185);
    const auto *dsk_186 = buffer.data(dsk + 186);
    const auto *dsk_187 = buffer.data(dsk + 187);
    const auto *dsk_189 = buffer.data(dsk + 189);
    const auto *dsk_190 = buffer.data(dsk + 190);
    const auto *dsk_191 = buffer.data(dsk + 191);
    const auto *dsk_192 = buffer.data(dsk + 192);
    const auto *dsk_194 = buffer.data(dsk + 194);
    const auto *dsk_195 = buffer.data(dsk + 195);
    const auto *dsk_196 = buffer.data(dsk + 196);
    const auto *dsk_197 = buffer.data(dsk + 197);
    const auto *dsk_198 = buffer.data(dsk + 198);
    const auto *dsk_200 = buffer.data(dsk + 200);
    const auto *dsk_201 = buffer.data(dsk + 201);
    const auto *dsk_202 = buffer.data(dsk + 202);
    const auto *dsk_203 = buffer.data(dsk + 203);
    const auto *dsk_204 = buffer.data(dsk + 204);
    const auto *dsk_205 = buffer.data(dsk + 205);
    const auto *dsk_207 = buffer.data(dsk + 207);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_x, pc_x, pc_y, psl0_130, psl0_131, \
                         psl0_132, psl1_130, psl1_131, psl1_132, \
                         dsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_x[k] * psl0_130[k]
                   - f_14 * pc_x[k] * psl1_130[k];

        t_131[k] = pa_x[k] * psl0_131[k]
                   - f_14 * pc_x[k] * psl1_131[k];

        t_132[k] = pa_x[k] * psl0_132[k]
                   - f_14 * pc_x[k] * psl1_132[k];

        t_133[k] = f_3 * pc_y[k] * dsk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pc_x, pc_z, psl0_134, psl1_134, \
                         dsi0_84, dsi0_85, dsi1_84, dsi1_85, dsk_108, \
                         dsk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_x[k] * psl0_134[k]
                   - f_14 * pc_x[k] * psl1_134[k];

        t_135[k] = f_1 * dsi0_84[k]
                   - f_2 * dsi1_84[k]
                   + f_3 * pc_x[k] * dsk_108[k];

        t_136[k] = f_20 * dsi0_85[k]
                   - f_21 * dsi1_85[k]
                   + f_3 * pc_x[k] * dsk_109[k];

        t_137[k] = f_3 * pc_z[k] * dsk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_z, dsi0_87, dsi0_89, dsi0_90, \
                         dsi1_87, dsi1_89, dsi1_90, dsk_109, dsk_111, dsk_113, \
                         dsk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_12 * dsi0_87[k]
                   - f_13 * dsi1_87[k]
                   + f_3 * pc_x[k] * dsk_111[k];

        t_139[k] = f_3 * pc_z[k] * dsk_109[k];

        t_140[k] = f_12 * dsi0_89[k]
                   - f_13 * dsi1_89[k]
                   + f_3 * pc_x[k] * dsk_113[k];

        t_141[k] = f_10 * dsi0_90[k]
                   - f_11 * dsi1_90[k]
                   + f_3 * pc_x[k] * dsk_114[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, dsi0_92, dsi0_93, dsi0_94, \
                         dsi1_92, dsi1_93, dsi1_94, dsk_111, dsk_116, dsk_117, \
                         dsk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * dsk_111[k];

        t_143[k] = f_10 * dsi0_92[k]
                   - f_11 * dsi1_92[k]
                   + f_3 * pc_x[k] * dsk_116[k];

        t_144[k] = f_10 * dsi0_93[k]
                   - f_11 * dsi1_93[k]
                   + f_3 * pc_x[k] * dsk_117[k];

        t_145[k] = f_8 * dsi0_94[k]
                   - f_9 * dsi1_94[k]
                   + f_3 * pc_x[k] * dsk_118[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_x, pc_z, dsi0_96, dsi0_97, dsi0_98, \
                         dsi1_96, dsi1_97, dsi1_98, dsk_114, dsk_120, dsk_121, \
                         dsk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_3 * pc_z[k] * dsk_114[k];

        t_147[k] = f_8 * dsi0_96[k]
                   - f_9 * dsi1_96[k]
                   + f_3 * pc_x[k] * dsk_120[k];

        t_148[k] = f_8 * dsi0_97[k]
                   - f_9 * dsi1_97[k]
                   + f_3 * pc_x[k] * dsk_121[k];

        t_149[k] = f_8 * dsi0_98[k]
                   - f_9 * dsi1_98[k]
                   + f_3 * pc_x[k] * dsk_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_z, dsi0_99, dsi0_101, dsi0_102, \
                         dsi1_99, dsi1_101, dsi1_102, dsk_118, dsk_123, dsk_125, \
                         dsk_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_6 * dsi0_99[k]
                   - f_7 * dsi1_99[k]
                   + f_3 * pc_x[k] * dsk_123[k];

        t_151[k] = f_3 * pc_z[k] * dsk_118[k];

        t_152[k] = f_6 * dsi0_101[k]
                   - f_7 * dsi1_101[k]
                   + f_3 * pc_x[k] * dsk_125[k];

        t_153[k] = f_6 * dsi0_102[k]
                   - f_7 * dsi1_102[k]
                   + f_3 * pc_x[k] * dsk_126[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_z, dsi0_103, dsi0_104, dsi0_105, \
                         dsi1_103, dsi1_104, dsi1_105, dsk_123, dsk_127, dsk_128, \
                         dsk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_6 * dsi0_103[k]
                   - f_7 * dsi1_103[k]
                   + f_3 * pc_x[k] * dsk_127[k];

        t_155[k] = f_6 * dsi0_104[k]
                   - f_7 * dsi1_104[k]
                   + f_3 * pc_x[k] * dsk_128[k];

        t_156[k] = f_4 * dsi0_105[k]
                   - f_5 * dsi1_105[k]
                   + f_3 * pc_x[k] * dsk_129[k];

        t_157[k] = f_3 * pc_z[k] * dsk_123[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pc_x, dsi0_107, dsi0_108, dsi0_109, dsi1_107, \
                         dsi1_108, dsi1_109, dsk_131, dsk_132, \
                         dsk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_4 * dsi0_107[k]
                   - f_5 * dsi1_107[k]
                   + f_3 * pc_x[k] * dsk_131[k];

        t_159[k] = f_4 * dsi0_108[k]
                   - f_5 * dsi1_108[k]
                   + f_3 * pc_x[k] * dsk_132[k];

        t_160[k] = f_4 * dsi0_109[k]
                   - f_5 * dsi1_109[k]
                   + f_3 * pc_x[k] * dsk_133[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, pc_x, dsi0_110, dsi0_111, \
                         dsi1_110, dsi1_111, dsk_134, dsk_135, dsk_136, dsk_137, \
                         dsk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_4 * dsi0_110[k]
                   - f_5 * dsi1_110[k]
                   + f_3 * pc_x[k] * dsk_134[k];

        t_162[k] = f_4 * dsi0_111[k]
                   - f_5 * dsi1_111[k]
                   + f_3 * pc_x[k] * dsk_135[k];

        t_163[k] = f_3 * pc_x[k] * dsk_136[k];

        t_164[k] = f_3 * pc_x[k] * dsk_137[k];

        t_165[k] = f_3 * pc_x[k] * dsk_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, dsk_139, dsk_140, dsk_141, \
                         dsk_142, dsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_3 * pc_x[k] * dsk_139[k];

        t_167[k] = f_3 * pc_x[k] * dsk_140[k];

        t_168[k] = f_3 * pc_x[k] * dsk_141[k];

        t_169[k] = f_3 * pc_x[k] * dsk_142[k];

        t_170[k] = f_3 * pc_x[k] * dsk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_y, pc_z, psk_64, dsi0_105, dsi0_106, \
                         dsi1_105, dsi1_106, dsk_136, dsk_137, \
                         dsk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_0 * psk_64[k]
                   + f_1 * dsi0_105[k]
                   - f_2 * dsi1_105[k]
                   + f_3 * pc_y[k] * dsk_136[k];

        t_172[k] = f_3 * pc_z[k] * dsk_136[k];

        t_173[k] = f_4 * dsi0_105[k]
                   - f_5 * dsi1_105[k]
                   + f_3 * pc_z[k] * dsk_137[k];

        t_174[k] = f_6 * dsi0_106[k]
                   - f_7 * dsi1_106[k]
                   + f_3 * pc_z[k] * dsk_138[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_z, dsi0_107, dsi0_108, dsi0_109, dsi1_107, \
                         dsi1_108, dsi1_109, dsk_139, dsk_140, \
                         dsk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_8 * dsi0_107[k]
                   - f_9 * dsi1_107[k]
                   + f_3 * pc_z[k] * dsk_139[k];

        t_176[k] = f_10 * dsi0_108[k]
                   - f_11 * dsi1_108[k]
                   + f_3 * pc_z[k] * dsk_140[k];

        t_177[k] = f_12 * dsi0_109[k]
                   - f_13 * dsi1_109[k]
                   + f_3 * pc_z[k] * dsk_141[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_y, pa_z, pc_y, pc_z, psl0_46, psl0_90, \
                         psk_71, psl1_46, psl1_90, dsi0_111, dsi1_111, \
                         dsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_0 * psk_71[k]
                   + f_3 * pc_y[k] * dsk_143[k];

        t_179[k] = f_1 * dsi0_111[k]
                   - f_2 * dsi1_111[k]
                   + f_3 * pc_z[k] * dsk_143[k];

        t_180[k] = pa_y[k] * psl0_90[k]
                   - f_14 * pc_y[k] * psl1_90[k];

        t_181[k] = pa_z[k] * psl0_46[k]
                   - f_14 * pc_z[k] * psl1_46[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pa_y, pa_z, pc_x, pc_y, pc_z, psl0_48, psl0_92, \
                         psl1_48, psl1_92, dsi0_116, dsi1_116, \
                         dsk_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = pa_y[k] * psl0_92[k]
                   - f_14 * pc_y[k] * psl1_92[k];

        t_183[k] = pa_z[k] * psl0_48[k]
                   - f_14 * pc_z[k] * psl1_48[k];

        t_184[k] = f_12 * dsi0_116[k]
                   - f_13 * dsi1_116[k]
                   + f_3 * pc_x[k] * dsk_148[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pa_y, pa_z, pc_x, pc_y, pc_z, psl0_51, psl0_95, \
                         psl1_51, psl1_95, dsi0_119, dsi1_119, \
                         dsk_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pa_y[k] * psl0_95[k]
                   - f_14 * pc_y[k] * psl1_95[k];

        t_186[k] = pa_z[k] * psl0_51[k]
                   - f_14 * pc_z[k] * psl1_51[k];

        t_187[k] = f_10 * dsi0_119[k]
                   - f_11 * dsi1_119[k]
                   + f_3 * pc_x[k] * dsk_151[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_y, pa_z, pc_x, pc_y, pc_z, psl0_55, psl0_99, \
                         psl1_55, psl1_99, dsi0_120, dsi1_120, \
                         dsk_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_10 * dsi0_120[k]
                   - f_11 * dsi1_120[k]
                   + f_3 * pc_x[k] * dsk_152[k];

        t_189[k] = pa_y[k] * psl0_99[k]
                   - f_14 * pc_y[k] * psl1_99[k];

        t_190[k] = pa_z[k] * psl0_55[k]
                   - f_14 * pc_z[k] * psl1_55[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_x, dsi0_123, dsi0_124, dsi0_125, dsi1_123, \
                         dsi1_124, dsi1_125, dsk_155, dsk_156, \
                         dsk_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_8 * dsi0_123[k]
                   - f_9 * dsi1_123[k]
                   + f_3 * pc_x[k] * dsk_155[k];

        t_192[k] = f_8 * dsi0_124[k]
                   - f_9 * dsi1_124[k]
                   + f_3 * pc_x[k] * dsk_156[k];

        t_193[k] = f_8 * dsi0_125[k]
                   - f_9 * dsi1_125[k]
                   + f_3 * pc_x[k] * dsk_157[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_y, pa_z, pc_x, pc_y, pc_z, psl0_60, psl0_104, \
                         psl1_60, psl1_104, dsi0_128, dsi1_128, \
                         dsk_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_y[k] * psl0_104[k]
                   - f_14 * pc_y[k] * psl1_104[k];

        t_195[k] = pa_z[k] * psl0_60[k]
                   - f_14 * pc_z[k] * psl1_60[k];

        t_196[k] = f_6 * dsi0_128[k]
                   - f_7 * dsi1_128[k]
                   + f_3 * pc_x[k] * dsk_160[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pc_x, dsi0_129, dsi0_130, dsi0_131, dsi1_129, \
                         dsi1_130, dsi1_131, dsk_161, dsk_162, \
                         dsk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_6 * dsi0_129[k]
                   - f_7 * dsi1_129[k]
                   + f_3 * pc_x[k] * dsk_161[k];

        t_198[k] = f_6 * dsi0_130[k]
                   - f_7 * dsi1_130[k]
                   + f_3 * pc_x[k] * dsk_162[k];

        t_199[k] = f_6 * dsi0_131[k]
                   - f_7 * dsi1_131[k]
                   + f_3 * pc_x[k] * dsk_163[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pa_y, pa_z, pc_x, pc_y, pc_z, psl0_66, psl0_110, \
                         psl1_66, psl1_110, dsi0_134, dsi1_134, \
                         dsk_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * psl0_110[k]
                   - f_14 * pc_y[k] * psl1_110[k];

        t_201[k] = pa_z[k] * psl0_66[k]
                   - f_14 * pc_z[k] * psl1_66[k];

        t_202[k] = f_4 * dsi0_134[k]
                   - f_5 * dsi1_134[k]
                   + f_3 * pc_x[k] * dsk_166[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, dsi0_135, dsi0_136, dsi0_137, dsi1_135, \
                         dsi1_136, dsi1_137, dsk_167, dsk_168, \
                         dsk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_4 * dsi0_135[k]
                   - f_5 * dsi1_135[k]
                   + f_3 * pc_x[k] * dsk_167[k];

        t_204[k] = f_4 * dsi0_136[k]
                   - f_5 * dsi1_136[k]
                   + f_3 * pc_x[k] * dsk_168[k];

        t_205[k] = f_4 * dsi0_137[k]
                   - f_5 * dsi1_137[k]
                   + f_3 * pc_x[k] * dsk_169[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, pa_y, pc_x, pc_y, psl0_117, \
                         psl1_117, dsi0_138, dsi1_138, dsk_170, dsk_172, dsk_173, \
                         dsk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_4 * dsi0_138[k]
                   - f_5 * dsi1_138[k]
                   + f_3 * pc_x[k] * dsk_170[k];

        t_207[k] = pa_y[k] * psl0_117[k]
                   - f_14 * pc_y[k] * psl1_117[k];

        t_208[k] = f_3 * pc_x[k] * dsk_172[k];

        t_209[k] = f_3 * pc_x[k] * dsk_173[k];

        t_210[k] = f_3 * pc_x[k] * dsk_174[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, psl0_81, \
                         psl1_81, dsk_175, dsk_176, dsk_177, dsk_178, \
                         dsk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_3 * pc_x[k] * dsk_175[k];

        t_212[k] = f_3 * pc_x[k] * dsk_176[k];

        t_213[k] = f_3 * pc_x[k] * dsk_177[k];

        t_214[k] = f_3 * pc_x[k] * dsk_178[k];

        t_215[k] = f_3 * pc_x[k] * dsk_179[k];

        t_216[k] = pa_z[k] * psl0_81[k]
                   - f_14 * pc_z[k] * psl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pa_y, pc_y, pc_z, psl0_128, psl0_129, psk_64, \
                         psk_102, psk_103, psl1_128, psl1_129, \
                         dsk_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * psk_64[k]
                   + f_3 * pc_z[k] * dsk_172[k];

        t_218[k] = pa_y[k] * psl0_128[k]
                   + f_16 * psk_102[k]
                   - f_14 * pc_y[k] * psl1_128[k];

        t_219[k] = pa_y[k] * psl0_129[k]
                   + f_17 * psk_103[k]
                   - f_14 * pc_y[k] * psl1_129[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pa_y, pc_y, psl0_130, psl0_131, psl0_132, \
                         psk_104, psk_105, psk_106, psl1_130, psl1_131, \
                         psl1_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_y[k] * psl0_130[k]
                   + f_18 * psk_104[k]
                   - f_14 * pc_y[k] * psl1_130[k];

        t_221[k] = pa_y[k] * psl0_131[k]
                   + f_19 * psk_105[k]
                   - f_14 * pc_y[k] * psl1_131[k];

        t_222[k] = pa_y[k] * psl0_132[k]
                   + f_0 * psk_106[k]
                   - f_14 * pc_y[k] * psl1_132[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, psl0_134, psk_107, \
                         psl1_134, dsi0_140, dsi1_140, dsk_179, \
                         dsk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * psk_107[k]
                   + f_3 * pc_y[k] * dsk_179[k];

        t_224[k] = pa_y[k] * psl0_134[k]
                   - f_14 * pc_y[k] * psl1_134[k];

        t_225[k] = f_1 * dsi0_140[k]
                   - f_2 * dsi1_140[k]
                   + f_3 * pc_x[k] * dsk_180[k];

        t_226[k] = f_3 * pc_y[k] * dsk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pc_x, pc_y, dsi0_142, dsi0_143, dsi0_145, \
                         dsi1_142, dsi1_143, dsi1_145, dsk_182, dsk_183, \
                         dsk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_20 * dsi0_142[k]
                   - f_21 * dsi1_142[k]
                   + f_3 * pc_x[k] * dsk_182[k];

        t_228[k] = f_12 * dsi0_143[k]
                   - f_13 * dsi1_143[k]
                   + f_3 * pc_x[k] * dsk_183[k];

        t_229[k] = f_3 * pc_y[k] * dsk_182[k];

        t_230[k] = f_12 * dsi0_145[k]
                   - f_13 * dsi1_145[k]
                   + f_3 * pc_x[k] * dsk_185[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pc_x, pc_y, dsi0_146, dsi0_147, dsi0_149, \
                         dsi1_146, dsi1_147, dsi1_149, dsk_185, dsk_186, dsk_187, \
                         dsk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_10 * dsi0_146[k]
                   - f_11 * dsi1_146[k]
                   + f_3 * pc_x[k] * dsk_186[k];

        t_232[k] = f_10 * dsi0_147[k]
                   - f_11 * dsi1_147[k]
                   + f_3 * pc_x[k] * dsk_187[k];

        t_233[k] = f_3 * pc_y[k] * dsk_185[k];

        t_234[k] = f_10 * dsi0_149[k]
                   - f_11 * dsi1_149[k]
                   + f_3 * pc_x[k] * dsk_189[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, dsi0_150, dsi0_151, dsi0_152, \
                         dsi1_150, dsi1_151, dsi1_152, dsk_189, dsk_190, dsk_191, \
                         dsk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_8 * dsi0_150[k]
                   - f_9 * dsi1_150[k]
                   + f_3 * pc_x[k] * dsk_190[k];

        t_236[k] = f_8 * dsi0_151[k]
                   - f_9 * dsi1_151[k]
                   + f_3 * pc_x[k] * dsk_191[k];

        t_237[k] = f_8 * dsi0_152[k]
                   - f_9 * dsi1_152[k]
                   + f_3 * pc_x[k] * dsk_192[k];

        t_238[k] = f_3 * pc_y[k] * dsk_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, dsi0_154, dsi0_155, dsi0_156, dsi1_154, \
                         dsi1_155, dsi1_156, dsk_194, dsk_195, \
                         dsk_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_8 * dsi0_154[k]
                   - f_9 * dsi1_154[k]
                   + f_3 * pc_x[k] * dsk_194[k];

        t_240[k] = f_6 * dsi0_155[k]
                   - f_7 * dsi1_155[k]
                   + f_3 * pc_x[k] * dsk_195[k];

        t_241[k] = f_6 * dsi0_156[k]
                   - f_7 * dsi1_156[k]
                   + f_3 * pc_x[k] * dsk_196[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pc_x, pc_y, dsi0_157, dsi0_158, dsi0_160, \
                         dsi1_157, dsi1_158, dsi1_160, dsk_194, dsk_197, dsk_198, \
                         dsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_6 * dsi0_157[k]
                   - f_7 * dsi1_157[k]
                   + f_3 * pc_x[k] * dsk_197[k];

        t_243[k] = f_6 * dsi0_158[k]
                   - f_7 * dsi1_158[k]
                   + f_3 * pc_x[k] * dsk_198[k];

        t_244[k] = f_3 * pc_y[k] * dsk_194[k];

        t_245[k] = f_6 * dsi0_160[k]
                   - f_7 * dsi1_160[k]
                   + f_3 * pc_x[k] * dsk_200[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pc_x, dsi0_161, dsi0_162, dsi0_163, dsi1_161, \
                         dsi1_162, dsi1_163, dsk_201, dsk_202, \
                         dsk_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_4 * dsi0_161[k]
                   - f_5 * dsi1_161[k]
                   + f_3 * pc_x[k] * dsk_201[k];

        t_247[k] = f_4 * dsi0_162[k]
                   - f_5 * dsi1_162[k]
                   + f_3 * pc_x[k] * dsk_202[k];

        t_248[k] = f_4 * dsi0_163[k]
                   - f_5 * dsi1_163[k]
                   + f_3 * pc_x[k] * dsk_203[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pc_x, pc_y, dsi0_164, dsi0_165, dsi0_167, \
                         dsi1_164, dsi1_165, dsi1_167, dsk_200, dsk_204, dsk_205, \
                         dsk_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_4 * dsi0_164[k]
                   - f_5 * dsi1_164[k]
                   + f_3 * pc_x[k] * dsk_204[k];

        t_250[k] = f_4 * dsi0_165[k]
                   - f_5 * dsi1_165[k]
                   + f_3 * pc_x[k] * dsk_205[k];

        t_251[k] = f_3 * pc_y[k] * dsk_200[k];

        t_252[k] = f_4 * dsi0_167[k]
                   - f_5 * dsi1_167[k]
                   + f_3 * pc_x[k] * dsk_207[k];
    }
}

static auto
compute_prim_dsl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t psk, const size_t dsi0,
                                                          const size_t dsi1, const size_t dsk,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psk_107 = buffer.data(psk + 107);

    const auto *dsi0_161 = buffer.data(dsi0 + 161);
    const auto *dsi0_162 = buffer.data(dsi0 + 162);
    const auto *dsi0_163 = buffer.data(dsi0 + 163);
    const auto *dsi0_164 = buffer.data(dsi0 + 164);
    const auto *dsi0_165 = buffer.data(dsi0 + 165);
    const auto *dsi0_166 = buffer.data(dsi0 + 166);
    const auto *dsi0_167 = buffer.data(dsi0 + 167);

    const auto *dsi1_161 = buffer.data(dsi1 + 161);
    const auto *dsi1_162 = buffer.data(dsi1 + 162);
    const auto *dsi1_163 = buffer.data(dsi1 + 163);
    const auto *dsi1_164 = buffer.data(dsi1 + 164);
    const auto *dsi1_165 = buffer.data(dsi1 + 165);
    const auto *dsi1_166 = buffer.data(dsi1 + 166);
    const auto *dsi1_167 = buffer.data(dsi1 + 167);

    const auto *dsk_208 = buffer.data(dsk + 208);
    const auto *dsk_209 = buffer.data(dsk + 209);
    const auto *dsk_210 = buffer.data(dsk + 210);
    const auto *dsk_211 = buffer.data(dsk + 211);
    const auto *dsk_212 = buffer.data(dsk + 212);
    const auto *dsk_213 = buffer.data(dsk + 213);
    const auto *dsk_214 = buffer.data(dsk + 214);
    const auto *dsk_215 = buffer.data(dsk + 215);

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, t_258, t_259, pc_x, dsk_208, \
                         dsk_209, dsk_210, dsk_211, dsk_212, dsk_213, \
                         dsk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_3 * pc_x[k] * dsk_208[k];

        t_254[k] = f_3 * pc_x[k] * dsk_209[k];

        t_255[k] = f_3 * pc_x[k] * dsk_210[k];

        t_256[k] = f_3 * pc_x[k] * dsk_211[k];

        t_257[k] = f_3 * pc_x[k] * dsk_212[k];

        t_258[k] = f_3 * pc_x[k] * dsk_213[k];

        t_259[k] = f_3 * pc_x[k] * dsk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, dsi0_161, dsi0_162, dsi0_163, \
                         dsi1_161, dsi1_162, dsi1_163, dsk_208, dsk_209, dsk_210, \
                         dsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_3 * pc_x[k] * dsk_215[k];

        t_261[k] = f_1 * dsi0_161[k]
                   - f_2 * dsi1_161[k]
                   + f_3 * pc_y[k] * dsk_208[k];

        t_262[k] = f_20 * dsi0_162[k]
                   - f_21 * dsi1_162[k]
                   + f_3 * pc_y[k] * dsk_209[k];

        t_263[k] = f_12 * dsi0_163[k]
                   - f_13 * dsi1_163[k]
                   + f_3 * pc_y[k] * dsk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, dsi0_164, dsi0_165, dsi0_166, dsi1_164, \
                         dsi1_165, dsi1_166, dsk_211, dsk_212, \
                         dsk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * dsi0_164[k]
                   - f_11 * dsi1_164[k]
                   + f_3 * pc_y[k] * dsk_211[k];

        t_265[k] = f_8 * dsi0_165[k]
                   - f_9 * dsi1_165[k]
                   + f_3 * pc_y[k] * dsk_212[k];

        t_266[k] = f_6 * dsi0_166[k]
                   - f_7 * dsi1_166[k]
                   + f_3 * pc_y[k] * dsk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, psk_107, dsi0_167, dsi1_167, \
                         dsk_214, dsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * dsi0_167[k]
                   - f_5 * dsi1_167[k]
                   + f_3 * pc_y[k] * dsk_214[k];

        t_268[k] = f_3 * pc_y[k] * dsk_215[k];

        t_269[k] = f_0 * psk_107[k]
                   + f_1 * dsi0_167[k]
                   - f_2 * dsi1_167[k]
                   + f_3 * pc_z[k] * dsk_215[k];
    }
}

auto
compute_prim_dsl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psl0, const size_t psk,
                                                   const size_t psl1, const size_t dsi0,
                                                   const size_t dsi1, const size_t dsk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dsl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, psl0, psk,
                                                              psl1, dsi0, dsi1, dsk, ncols,
                                                              gamma, p, q);

    compute_prim_dsl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, psl0, psk,
                                                              psl1, dsi0, dsi1, dsk, ncols,
                                                              gamma, p, q);

    compute_prim_dsl_three_center_electron_repulsion_0_piece2(buffer, target, pc, psk, dsi0,
                                                              dsi1, dsk, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
