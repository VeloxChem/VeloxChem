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


#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fsk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsk0,
                                                          const size_t dsi, const size_t dsk1,
                                                          const size_t fsh0, const size_t fsh1,
                                                          const size_t fsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_15 = 2.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsk0_0 = buffer.data(dsk0 + 0);
    const auto *dsk0_3 = buffer.data(dsk0 + 3);
    const auto *dsk0_5 = buffer.data(dsk0 + 5);
    const auto *dsk0_6 = buffer.data(dsk0 + 6);
    const auto *dsk0_9 = buffer.data(dsk0 + 9);
    const auto *dsk0_10 = buffer.data(dsk0 + 10);
    const auto *dsk0_14 = buffer.data(dsk0 + 14);
    const auto *dsk0_15 = buffer.data(dsk0 + 15);
    const auto *dsk0_20 = buffer.data(dsk0 + 20);
    const auto *dsk0_28 = buffer.data(dsk0 + 28);
    const auto *dsk0_35 = buffer.data(dsk0 + 35);
    const auto *dsk0_108 = buffer.data(dsk0 + 108);
    const auto *dsk0_111 = buffer.data(dsk0 + 111);
    const auto *dsk0_114 = buffer.data(dsk0 + 114);
    const auto *dsk0_118 = buffer.data(dsk0 + 118);
    const auto *dsk0_123 = buffer.data(dsk0 + 123);

    const auto *dsi_0 = buffer.data(dsi + 0);
    const auto *dsi_1 = buffer.data(dsi + 1);
    const auto *dsi_2 = buffer.data(dsi + 2);
    const auto *dsi_3 = buffer.data(dsi + 3);
    const auto *dsi_5 = buffer.data(dsi + 5);
    const auto *dsi_6 = buffer.data(dsi + 6);
    const auto *dsi_9 = buffer.data(dsi + 9);
    const auto *dsi_10 = buffer.data(dsi + 10);
    const auto *dsi_14 = buffer.data(dsi + 14);
    const auto *dsi_21 = buffer.data(dsi + 21);
    const auto *dsi_23 = buffer.data(dsi + 23);
    const auto *dsi_24 = buffer.data(dsi + 24);
    const auto *dsi_25 = buffer.data(dsi + 25);
    const auto *dsi_27 = buffer.data(dsi + 27);
    const auto *dsi_28 = buffer.data(dsi + 28);
    const auto *dsi_33 = buffer.data(dsi + 33);
    const auto *dsi_37 = buffer.data(dsi + 37);
    const auto *dsi_42 = buffer.data(dsi + 42);
    const auto *dsi_49 = buffer.data(dsi + 49);
    const auto *dsi_51 = buffer.data(dsi + 51);
    const auto *dsi_52 = buffer.data(dsi + 52);
    const auto *dsi_53 = buffer.data(dsi + 53);
    const auto *dsi_54 = buffer.data(dsi + 54);
    const auto *dsi_55 = buffer.data(dsi + 55);
    const auto *dsi_77 = buffer.data(dsi + 77);
    const auto *dsi_78 = buffer.data(dsi + 78);
    const auto *dsi_79 = buffer.data(dsi + 79);
    const auto *dsi_80 = buffer.data(dsi + 80);
    const auto *dsi_81 = buffer.data(dsi + 81);
    const auto *dsi_83 = buffer.data(dsi + 83);
    const auto *dsi_84 = buffer.data(dsi + 84);
    const auto *dsi_87 = buffer.data(dsi + 87);
    const auto *dsi_90 = buffer.data(dsi + 90);
    const auto *dsi_94 = buffer.data(dsi + 94);
    const auto *dsi_99 = buffer.data(dsi + 99);

    const auto *dsk1_0 = buffer.data(dsk1 + 0);
    const auto *dsk1_3 = buffer.data(dsk1 + 3);
    const auto *dsk1_5 = buffer.data(dsk1 + 5);
    const auto *dsk1_6 = buffer.data(dsk1 + 6);
    const auto *dsk1_9 = buffer.data(dsk1 + 9);
    const auto *dsk1_10 = buffer.data(dsk1 + 10);
    const auto *dsk1_14 = buffer.data(dsk1 + 14);
    const auto *dsk1_15 = buffer.data(dsk1 + 15);
    const auto *dsk1_20 = buffer.data(dsk1 + 20);
    const auto *dsk1_28 = buffer.data(dsk1 + 28);
    const auto *dsk1_35 = buffer.data(dsk1 + 35);
    const auto *dsk1_108 = buffer.data(dsk1 + 108);
    const auto *dsk1_111 = buffer.data(dsk1 + 111);
    const auto *dsk1_114 = buffer.data(dsk1 + 114);
    const auto *dsk1_118 = buffer.data(dsk1 + 118);
    const auto *dsk1_123 = buffer.data(dsk1 + 123);

    const auto *fsh0_0 = buffer.data(fsh0 + 0);
    const auto *fsh0_1 = buffer.data(fsh0 + 1);
    const auto *fsh0_2 = buffer.data(fsh0 + 2);
    const auto *fsh0_3 = buffer.data(fsh0 + 3);
    const auto *fsh0_5 = buffer.data(fsh0 + 5);
    const auto *fsh0_6 = buffer.data(fsh0 + 6);
    const auto *fsh0_8 = buffer.data(fsh0 + 8);
    const auto *fsh0_9 = buffer.data(fsh0 + 9);
    const auto *fsh0_15 = buffer.data(fsh0 + 15);
    const auto *fsh0_17 = buffer.data(fsh0 + 17);
    const auto *fsh0_18 = buffer.data(fsh0 + 18);
    const auto *fsh0_19 = buffer.data(fsh0 + 19);
    const auto *fsh0_20 = buffer.data(fsh0 + 20);
    const auto *fsh0_24 = buffer.data(fsh0 + 24);
    const auto *fsh0_27 = buffer.data(fsh0 + 27);
    const auto *fsh0_28 = buffer.data(fsh0 + 28);
    const auto *fsh0_36 = buffer.data(fsh0 + 36);
    const auto *fsh0_37 = buffer.data(fsh0 + 37);
    const auto *fsh0_38 = buffer.data(fsh0 + 38);
    const auto *fsh0_39 = buffer.data(fsh0 + 39);
    const auto *fsh0_44 = buffer.data(fsh0 + 44);
    const auto *fsh0_46 = buffer.data(fsh0 + 46);
    const auto *fsh0_47 = buffer.data(fsh0 + 47);
    const auto *fsh0_49 = buffer.data(fsh0 + 49);
    const auto *fsh0_50 = buffer.data(fsh0 + 50);
    const auto *fsh0_51 = buffer.data(fsh0 + 51);
    const auto *fsh0_58 = buffer.data(fsh0 + 58);
    const auto *fsh0_59 = buffer.data(fsh0 + 59);
    const auto *fsh0_60 = buffer.data(fsh0 + 60);
    const auto *fsh0_61 = buffer.data(fsh0 + 61);
    const auto *fsh0_62 = buffer.data(fsh0 + 62);
    const auto *fsh0_63 = buffer.data(fsh0 + 63);
    const auto *fsh0_65 = buffer.data(fsh0 + 65);
    const auto *fsh0_66 = buffer.data(fsh0 + 66);
    const auto *fsh0_68 = buffer.data(fsh0 + 68);
    const auto *fsh0_69 = buffer.data(fsh0 + 69);
    const auto *fsh0_70 = buffer.data(fsh0 + 70);
    const auto *fsh0_72 = buffer.data(fsh0 + 72);

    const auto *fsh1_0 = buffer.data(fsh1 + 0);
    const auto *fsh1_1 = buffer.data(fsh1 + 1);
    const auto *fsh1_2 = buffer.data(fsh1 + 2);
    const auto *fsh1_3 = buffer.data(fsh1 + 3);
    const auto *fsh1_5 = buffer.data(fsh1 + 5);
    const auto *fsh1_6 = buffer.data(fsh1 + 6);
    const auto *fsh1_8 = buffer.data(fsh1 + 8);
    const auto *fsh1_9 = buffer.data(fsh1 + 9);
    const auto *fsh1_15 = buffer.data(fsh1 + 15);
    const auto *fsh1_17 = buffer.data(fsh1 + 17);
    const auto *fsh1_18 = buffer.data(fsh1 + 18);
    const auto *fsh1_19 = buffer.data(fsh1 + 19);
    const auto *fsh1_20 = buffer.data(fsh1 + 20);
    const auto *fsh1_24 = buffer.data(fsh1 + 24);
    const auto *fsh1_27 = buffer.data(fsh1 + 27);
    const auto *fsh1_28 = buffer.data(fsh1 + 28);
    const auto *fsh1_36 = buffer.data(fsh1 + 36);
    const auto *fsh1_37 = buffer.data(fsh1 + 37);
    const auto *fsh1_38 = buffer.data(fsh1 + 38);
    const auto *fsh1_39 = buffer.data(fsh1 + 39);
    const auto *fsh1_44 = buffer.data(fsh1 + 44);
    const auto *fsh1_46 = buffer.data(fsh1 + 46);
    const auto *fsh1_47 = buffer.data(fsh1 + 47);
    const auto *fsh1_49 = buffer.data(fsh1 + 49);
    const auto *fsh1_50 = buffer.data(fsh1 + 50);
    const auto *fsh1_51 = buffer.data(fsh1 + 51);
    const auto *fsh1_58 = buffer.data(fsh1 + 58);
    const auto *fsh1_59 = buffer.data(fsh1 + 59);
    const auto *fsh1_60 = buffer.data(fsh1 + 60);
    const auto *fsh1_61 = buffer.data(fsh1 + 61);
    const auto *fsh1_62 = buffer.data(fsh1 + 62);
    const auto *fsh1_63 = buffer.data(fsh1 + 63);
    const auto *fsh1_65 = buffer.data(fsh1 + 65);
    const auto *fsh1_66 = buffer.data(fsh1 + 66);
    const auto *fsh1_68 = buffer.data(fsh1 + 68);
    const auto *fsh1_69 = buffer.data(fsh1 + 69);
    const auto *fsh1_70 = buffer.data(fsh1 + 70);
    const auto *fsh1_72 = buffer.data(fsh1 + 72);

    const auto *fsi_0 = buffer.data(fsi + 0);
    const auto *fsi_1 = buffer.data(fsi + 1);
    const auto *fsi_2 = buffer.data(fsi + 2);
    const auto *fsi_3 = buffer.data(fsi + 3);
    const auto *fsi_5 = buffer.data(fsi + 5);
    const auto *fsi_6 = buffer.data(fsi + 6);
    const auto *fsi_8 = buffer.data(fsi + 8);
    const auto *fsi_9 = buffer.data(fsi + 9);
    const auto *fsi_10 = buffer.data(fsi + 10);
    const auto *fsi_12 = buffer.data(fsi + 12);
    const auto *fsi_13 = buffer.data(fsi + 13);
    const auto *fsi_14 = buffer.data(fsi + 14);
    const auto *fsi_15 = buffer.data(fsi + 15);
    const auto *fsi_20 = buffer.data(fsi + 20);
    const auto *fsi_21 = buffer.data(fsi + 21);
    const auto *fsi_23 = buffer.data(fsi + 23);
    const auto *fsi_24 = buffer.data(fsi + 24);
    const auto *fsi_25 = buffer.data(fsi + 25);
    const auto *fsi_26 = buffer.data(fsi + 26);
    const auto *fsi_27 = buffer.data(fsi + 27);
    const auto *fsi_28 = buffer.data(fsi + 28);
    const auto *fsi_29 = buffer.data(fsi + 29);
    const auto *fsi_31 = buffer.data(fsi + 31);
    const auto *fsi_33 = buffer.data(fsi + 33);
    const auto *fsi_34 = buffer.data(fsi + 34);
    const auto *fsi_35 = buffer.data(fsi + 35);
    const auto *fsi_37 = buffer.data(fsi + 37);
    const auto *fsi_38 = buffer.data(fsi + 38);
    const auto *fsi_39 = buffer.data(fsi + 39);
    const auto *fsi_40 = buffer.data(fsi + 40);
    const auto *fsi_42 = buffer.data(fsi + 42);
    const auto *fsi_43 = buffer.data(fsi + 43);
    const auto *fsi_49 = buffer.data(fsi + 49);
    const auto *fsi_50 = buffer.data(fsi + 50);
    const auto *fsi_51 = buffer.data(fsi + 51);
    const auto *fsi_52 = buffer.data(fsi + 52);
    const auto *fsi_53 = buffer.data(fsi + 53);
    const auto *fsi_54 = buffer.data(fsi + 54);
    const auto *fsi_55 = buffer.data(fsi + 55);
    const auto *fsi_56 = buffer.data(fsi + 56);
    const auto *fsi_58 = buffer.data(fsi + 58);
    const auto *fsi_60 = buffer.data(fsi + 60);
    const auto *fsi_61 = buffer.data(fsi + 61);
    const auto *fsi_63 = buffer.data(fsi + 63);
    const auto *fsi_64 = buffer.data(fsi + 64);
    const auto *fsi_65 = buffer.data(fsi + 65);
    const auto *fsi_67 = buffer.data(fsi + 67);
    const auto *fsi_68 = buffer.data(fsi + 68);
    const auto *fsi_69 = buffer.data(fsi + 69);
    const auto *fsi_70 = buffer.data(fsi + 70);
    const auto *fsi_76 = buffer.data(fsi + 76);
    const auto *fsi_77 = buffer.data(fsi + 77);
    const auto *fsi_78 = buffer.data(fsi + 78);
    const auto *fsi_79 = buffer.data(fsi + 79);
    const auto *fsi_80 = buffer.data(fsi + 80);
    const auto *fsi_81 = buffer.data(fsi + 81);
    const auto *fsi_82 = buffer.data(fsi + 82);
    const auto *fsi_83 = buffer.data(fsi + 83);
    const auto *fsi_84 = buffer.data(fsi + 84);
    const auto *fsi_85 = buffer.data(fsi + 85);
    const auto *fsi_86 = buffer.data(fsi + 86);
    const auto *fsi_87 = buffer.data(fsi + 87);
    const auto *fsi_89 = buffer.data(fsi + 89);
    const auto *fsi_90 = buffer.data(fsi + 90);
    const auto *fsi_91 = buffer.data(fsi + 91);
    const auto *fsi_93 = buffer.data(fsi + 93);
    const auto *fsi_94 = buffer.data(fsi + 94);
    const auto *fsi_95 = buffer.data(fsi + 95);
    const auto *fsi_96 = buffer.data(fsi + 96);
    const auto *fsi_98 = buffer.data(fsi + 98);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dsi_0, fsh0_0, \
                         fsh1_0, fsi_0, fsi_1, fsi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsi_0[k]
                 + f_1 * fsh0_0[k]
                 - f_2 * fsh1_0[k]
                 + f_3 * pc_x[k] * fsi_0[k];

        t_1[k] = f_3 * pc_y[k] * fsi_0[k];

        t_2[k] = f_3 * pc_z[k] * fsi_0[k];

        t_3[k] = f_4 * fsh0_0[k]
                 - f_5 * fsh1_0[k]
                 + f_3 * pc_y[k] * fsi_1[k];

        t_4[k] = f_3 * pc_y[k] * fsi_2[k];

        t_5[k] = f_4 * fsh0_0[k]
                 - f_5 * fsh1_0[k]
                 + f_3 * pc_z[k] * fsi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, fsh0_1, fsh0_2, fsh0_3, fsh1_1, \
                         fsh1_2, fsh1_3, fsi_3, fsi_5, fsi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * fsh0_1[k]
                 - f_7 * fsh1_1[k]
                 + f_3 * pc_y[k] * fsi_3[k];

        t_7[k] = f_3 * pc_z[k] * fsi_3[k];

        t_8[k] = f_3 * pc_y[k] * fsi_5[k];

        t_9[k] = f_6 * fsh0_2[k]
                 - f_7 * fsh1_2[k]
                 + f_3 * pc_z[k] * fsi_5[k];

        t_10[k] = f_8 * fsh0_3[k]
                  - f_9 * fsh1_3[k]
                  + f_3 * pc_y[k] * fsi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, fsh0_5, fsh0_6, \
                         fsh1_5, fsh1_6, fsi_6, fsi_8, fsi_9, fsi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * fsi_6[k];

        t_12[k] = f_4 * fsh0_5[k]
                  - f_5 * fsh1_5[k]
                  + f_3 * pc_y[k] * fsi_8[k];

        t_13[k] = f_3 * pc_y[k] * fsi_9[k];

        t_14[k] = f_8 * fsh0_5[k]
                  - f_9 * fsh1_5[k]
                  + f_3 * pc_z[k] * fsi_9[k];

        t_15[k] = f_10 * fsh0_6[k]
                  - f_11 * fsh1_6[k]
                  + f_3 * pc_y[k] * fsi_10[k];

        t_16[k] = f_3 * pc_z[k] * fsi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, fsh0_8, fsh0_9, fsh1_8, fsh1_9, \
                         fsi_12, fsi_13, fsi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fsh0_8[k]
                  - f_7 * fsh1_8[k]
                  + f_3 * pc_y[k] * fsi_12[k];

        t_18[k] = f_4 * fsh0_9[k]
                  - f_5 * fsh1_9[k]
                  + f_3 * pc_y[k] * fsi_13[k];

        t_19[k] = f_3 * pc_y[k] * fsi_14[k];

        t_20[k] = f_10 * fsh0_9[k]
                  - f_11 * fsh1_9[k]
                  + f_3 * pc_z[k] * fsi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, dsi_21, dsi_23, dsi_24, \
                         dsi_25, fsi_15, fsi_21, fsi_23, fsi_24, \
                         fsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * dsi_21[k]
                  + f_3 * pc_x[k] * fsi_21[k];

        t_22[k] = f_3 * pc_z[k] * fsi_15[k];

        t_23[k] = f_0 * dsi_23[k]
                  + f_3 * pc_x[k] * fsi_23[k];

        t_24[k] = f_0 * dsi_24[k]
                  + f_3 * pc_x[k] * fsi_24[k];

        t_25[k] = f_0 * dsi_25[k]
                  + f_3 * pc_x[k] * fsi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, dsi_27, fsh0_15, fsh1_15, \
                         fsi_20, fsi_21, fsi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * fsi_20[k];

        t_27[k] = f_0 * dsi_27[k]
                  + f_3 * pc_x[k] * fsi_27[k];

        t_28[k] = f_1 * fsh0_15[k]
                  - f_2 * fsh1_15[k]
                  + f_3 * pc_y[k] * fsi_21[k];

        t_29[k] = f_3 * pc_z[k] * fsi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, fsh0_17, fsh0_18, fsh0_19, fsh1_17, fsh1_18, \
                         fsh1_19, fsi_23, fsi_24, fsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * fsh0_17[k]
                  - f_11 * fsh1_17[k]
                  + f_3 * pc_y[k] * fsi_23[k];

        t_31[k] = f_8 * fsh0_18[k]
                  - f_9 * fsh1_18[k]
                  + f_3 * pc_y[k] * fsi_24[k];

        t_32[k] = f_6 * fsh0_19[k]
                  - f_7 * fsh1_19[k]
                  + f_3 * pc_y[k] * fsi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, dsk0_0, dsi_0, \
                         dsk1_0, fsh0_20, fsh1_20, fsi_26, fsi_27, \
                         fsi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * fsh0_20[k]
                  - f_5 * fsh1_20[k]
                  + f_3 * pc_y[k] * fsi_26[k];

        t_34[k] = f_3 * pc_y[k] * fsi_27[k];

        t_35[k] = f_1 * fsh0_20[k]
                  - f_2 * fsh1_20[k]
                  + f_3 * pc_z[k] * fsi_27[k];

        t_36[k] = pa_y[k] * dsk0_0[k]
                  - f_12 * pc_y[k] * dsk1_0[k];

        t_37[k] = f_13 * dsi_0[k]
                  + f_3 * pc_y[k] * fsi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, dsk0_3, dsk0_5, dsi_1, \
                         dsk1_3, dsk1_5, fsi_28, fsi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * fsi_28[k];

        t_39[k] = pa_y[k] * dsk0_3[k]
                  + f_14 * dsi_1[k]
                  - f_12 * pc_y[k] * dsk1_3[k];

        t_40[k] = f_3 * pc_z[k] * fsi_29[k];

        t_41[k] = pa_y[k] * dsk0_5[k]
                  - f_12 * pc_y[k] * dsk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, dsk0_6, dsk0_9, dsi_3, \
                         dsi_5, dsk1_6, dsk1_9, fsi_31, fsi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * dsk0_6[k]
                  + f_0 * dsi_3[k]
                  - f_12 * pc_y[k] * dsk1_6[k];

        t_43[k] = f_3 * pc_z[k] * fsi_31[k];

        t_44[k] = f_13 * dsi_5[k]
                  + f_3 * pc_y[k] * fsi_33[k];

        t_45[k] = pa_y[k] * dsk0_9[k]
                  - f_12 * pc_y[k] * dsk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, dsk0_10, dsi_6, dsi_9, \
                         dsk1_10, fsh0_24, fsh1_24, fsi_34, fsi_35, \
                         fsi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * dsk0_10[k]
                  + f_15 * dsi_6[k]
                  - f_12 * pc_y[k] * dsk1_10[k];

        t_47[k] = f_3 * pc_z[k] * fsi_34[k];

        t_48[k] = f_4 * fsh0_24[k]
                  - f_5 * fsh1_24[k]
                  + f_3 * pc_z[k] * fsi_35[k];

        t_49[k] = f_13 * dsi_9[k]
                  + f_3 * pc_y[k] * fsi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, dsk0_14, dsk0_15, dsi_10, \
                         dsk1_14, dsk1_15, fsh0_27, fsh1_27, fsi_38, \
                         fsi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * dsk0_14[k]
                  - f_12 * pc_y[k] * dsk1_14[k];

        t_51[k] = pa_y[k] * dsk0_15[k]
                  + f_16 * dsi_10[k]
                  - f_12 * pc_y[k] * dsk1_15[k];

        t_52[k] = f_3 * pc_z[k] * fsi_38[k];

        t_53[k] = f_4 * fsh0_27[k]
                  - f_5 * fsh1_27[k]
                  + f_3 * pc_z[k] * fsi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, dsk0_20, dsi_14, dsk1_20, \
                         fsh0_28, fsh1_28, fsi_40, fsi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * fsh0_28[k]
                  - f_7 * fsh1_28[k]
                  + f_3 * pc_z[k] * fsi_40[k];

        t_55[k] = f_13 * dsi_14[k]
                  + f_3 * pc_y[k] * fsi_42[k];

        t_56[k] = pa_y[k] * dsk0_20[k]
                  - f_12 * pc_y[k] * dsk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, dsi_49, dsi_51, dsi_52, \
                         dsi_53, fsi_43, fsi_49, fsi_51, fsi_52, \
                         fsi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_14 * dsi_49[k]
                  + f_3 * pc_x[k] * fsi_49[k];

        t_58[k] = f_3 * pc_z[k] * fsi_43[k];

        t_59[k] = f_14 * dsi_51[k]
                  + f_3 * pc_x[k] * fsi_51[k];

        t_60[k] = f_14 * dsi_52[k]
                  + f_3 * pc_x[k] * fsi_52[k];

        t_61[k] = f_14 * dsi_53[k]
                  + f_3 * pc_x[k] * fsi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, dsi_21, dsi_54, dsi_55, \
                         fsh0_36, fsh1_36, fsi_49, fsi_54, fsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_14 * dsi_54[k]
                  + f_3 * pc_x[k] * fsi_54[k];

        t_63[k] = f_14 * dsi_55[k]
                  + f_3 * pc_x[k] * fsi_55[k];

        t_64[k] = f_13 * dsi_21[k]
                  + f_1 * fsh0_36[k]
                  - f_2 * fsh1_36[k]
                  + f_3 * pc_y[k] * fsi_49[k];

        t_65[k] = f_3 * pc_z[k] * fsi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, fsh0_36, fsh0_37, fsh0_38, fsh1_36, fsh1_37, \
                         fsh1_38, fsi_50, fsi_51, fsi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * fsh0_36[k]
                  - f_5 * fsh1_36[k]
                  + f_3 * pc_z[k] * fsi_50[k];

        t_67[k] = f_6 * fsh0_37[k]
                  - f_7 * fsh1_37[k]
                  + f_3 * pc_z[k] * fsi_51[k];

        t_68[k] = f_8 * fsh0_38[k]
                  - f_9 * fsh1_38[k]
                  + f_3 * pc_z[k] * fsi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, dsk0_35, dsi_27, dsk1_35, \
                         fsh0_39, fsh1_39, fsi_53, fsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * fsh0_39[k]
                  - f_11 * fsh1_39[k]
                  + f_3 * pc_z[k] * fsi_53[k];

        t_70[k] = f_13 * dsi_27[k]
                  + f_3 * pc_y[k] * fsi_55[k];

        t_71[k] = pa_y[k] * dsk0_35[k]
                  - f_12 * pc_y[k] * dsk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, dsk0_0, dsk0_3, \
                         dsi_0, dsk1_0, dsk1_3, fsi_56, fsi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * dsk0_0[k]
                  - f_12 * pc_z[k] * dsk1_0[k];

        t_73[k] = f_3 * pc_y[k] * fsi_56[k];

        t_74[k] = f_13 * dsi_0[k]
                  + f_3 * pc_z[k] * fsi_56[k];

        t_75[k] = pa_z[k] * dsk0_3[k]
                  - f_12 * pc_z[k] * dsk1_3[k];

        t_76[k] = f_3 * pc_y[k] * fsi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, dsk0_5, dsk0_6, dsi_2, \
                         dsk1_5, dsk1_6, fsh0_44, fsh1_44, fsi_60, \
                         fsi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * dsk0_5[k]
                  + f_14 * dsi_2[k]
                  - f_12 * pc_z[k] * dsk1_5[k];

        t_78[k] = pa_z[k] * dsk0_6[k]
                  - f_12 * pc_z[k] * dsk1_6[k];

        t_79[k] = f_4 * fsh0_44[k]
                  - f_5 * fsh1_44[k]
                  + f_3 * pc_y[k] * fsi_60[k];

        t_80[k] = f_3 * pc_y[k] * fsi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, dsk0_9, dsk0_10, dsi_5, dsk1_9, \
                         dsk1_10, fsh0_46, fsh1_46, fsi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * dsk0_9[k]
                  + f_0 * dsi_5[k]
                  - f_12 * pc_z[k] * dsk1_9[k];

        t_82[k] = pa_z[k] * dsk0_10[k]
                  - f_12 * pc_z[k] * dsk1_10[k];

        t_83[k] = f_6 * fsh0_46[k]
                  - f_7 * fsh1_46[k]
                  + f_3 * pc_y[k] * fsi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, dsk0_14, dsk0_15, dsi_9, \
                         dsk1_14, dsk1_15, fsh0_47, fsh1_47, fsi_64, \
                         fsi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * fsh0_47[k]
                  - f_5 * fsh1_47[k]
                  + f_3 * pc_y[k] * fsi_64[k];

        t_85[k] = f_3 * pc_y[k] * fsi_65[k];

        t_86[k] = pa_z[k] * dsk0_14[k]
                  + f_15 * dsi_9[k]
                  - f_12 * pc_z[k] * dsk1_14[k];

        t_87[k] = pa_z[k] * dsk0_15[k]
                  - f_12 * pc_z[k] * dsk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, fsh0_49, fsh0_50, fsh0_51, fsh1_49, \
                         fsh1_50, fsh1_51, fsi_67, fsi_68, fsi_69, \
                         fsi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * fsh0_49[k]
                  - f_9 * fsh1_49[k]
                  + f_3 * pc_y[k] * fsi_67[k];

        t_89[k] = f_6 * fsh0_50[k]
                  - f_7 * fsh1_50[k]
                  + f_3 * pc_y[k] * fsi_68[k];

        t_90[k] = f_4 * fsh0_51[k]
                  - f_5 * fsh1_51[k]
                  + f_3 * pc_y[k] * fsi_69[k];

        t_91[k] = f_3 * pc_y[k] * fsi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, dsk0_20, dsi_14, dsi_77, \
                         dsi_78, dsi_79, dsk1_20, fsi_77, fsi_78, \
                         fsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * dsk0_20[k]
                  + f_16 * dsi_14[k]
                  - f_12 * pc_z[k] * dsk1_20[k];

        t_93[k] = f_14 * dsi_77[k]
                  + f_3 * pc_x[k] * fsi_77[k];

        t_94[k] = f_14 * dsi_78[k]
                  + f_3 * pc_x[k] * fsi_78[k];

        t_95[k] = f_14 * dsi_79[k]
                  + f_3 * pc_x[k] * fsi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, dsi_80, dsi_81, dsi_83, fsi_76, \
                         fsi_80, fsi_81, fsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_14 * dsi_80[k]
                  + f_3 * pc_x[k] * fsi_80[k];

        t_97[k] = f_14 * dsi_81[k]
                  + f_3 * pc_x[k] * fsi_81[k];

        t_98[k] = f_3 * pc_y[k] * fsi_76[k];

        t_99[k] = f_14 * dsi_83[k]
                  + f_3 * pc_x[k] * fsi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, dsk0_28, dsk1_28, fsh0_58, \
                         fsh0_59, fsh1_58, fsh1_59, fsi_78, fsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * dsk0_28[k]
                   - f_12 * pc_z[k] * dsk1_28[k];

        t_101[k] = f_17 * fsh0_58[k]
                   - f_18 * fsh1_58[k]
                   + f_3 * pc_y[k] * fsi_78[k];

        t_102[k] = f_10 * fsh0_59[k]
                   - f_11 * fsh1_59[k]
                   + f_3 * pc_y[k] * fsi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, fsh0_60, fsh0_61, fsh0_62, fsh1_60, \
                         fsh1_61, fsh1_62, fsi_80, fsi_81, fsi_82, \
                         fsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * fsh0_60[k]
                   - f_9 * fsh1_60[k]
                   + f_3 * pc_y[k] * fsi_80[k];

        t_104[k] = f_6 * fsh0_61[k]
                   - f_7 * fsh1_61[k]
                   + f_3 * pc_y[k] * fsi_81[k];

        t_105[k] = f_4 * fsh0_62[k]
                   - f_5 * fsh1_62[k]
                   + f_3 * pc_y[k] * fsi_82[k];

        t_106[k] = f_3 * pc_y[k] * fsi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_x, pc_x, pc_y, pc_z, dsk0_108, dsi_27, \
                         dsi_28, dsi_84, dsk1_108, fsh0_62, fsh1_62, fsi_83, \
                         fsi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * dsi_27[k]
                   + f_1 * fsh0_62[k]
                   - f_2 * fsh1_62[k]
                   + f_3 * pc_z[k] * fsi_83[k];

        t_108[k] = pa_x[k] * dsk0_108[k]
                   + f_19 * dsi_84[k]
                   - f_12 * pc_x[k] * dsk1_108[k];

        t_109[k] = f_14 * dsi_28[k]
                   + f_3 * pc_y[k] * fsi_84[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_x, pc_x, pc_z, dsk0_111, dsi_87, \
                         dsk1_111, fsh0_63, fsh1_63, fsi_84, fsi_85, \
                         fsi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * pc_z[k] * fsi_84[k];

        t_111[k] = pa_x[k] * dsk0_111[k]
                   + f_16 * dsi_87[k]
                   - f_12 * pc_x[k] * dsk1_111[k];

        t_112[k] = f_3 * pc_z[k] * fsi_85[k];

        t_113[k] = f_4 * fsh0_63[k]
                   - f_5 * fsh1_63[k]
                   + f_3 * pc_z[k] * fsi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_x, pc_x, pc_y, pc_z, dsk0_114, dsi_33, \
                         dsi_90, dsk1_114, fsh0_65, fsh1_65, fsi_87, \
                         fsi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_x[k] * dsk0_114[k]
                   + f_15 * dsi_90[k]
                   - f_12 * pc_x[k] * dsk1_114[k];

        t_115[k] = f_3 * pc_z[k] * fsi_87[k];

        t_116[k] = f_14 * dsi_33[k]
                   + f_3 * pc_y[k] * fsi_89[k];

        t_117[k] = f_6 * fsh0_65[k]
                   - f_7 * fsh1_65[k]
                   + f_3 * pc_z[k] * fsi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_x, pc_x, pc_z, dsk0_118, dsi_94, dsk1_118, \
                         fsh0_66, fsh1_66, fsi_90, fsi_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_x[k] * dsk0_118[k]
                   + f_0 * dsi_94[k]
                   - f_12 * pc_x[k] * dsk1_118[k];

        t_119[k] = f_3 * pc_z[k] * fsi_90[k];

        t_120[k] = f_4 * fsh0_66[k]
                   - f_5 * fsh1_66[k]
                   + f_3 * pc_z[k] * fsi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pc_x, pc_y, pc_z, dsk0_123, dsi_37, \
                         dsi_99, dsk1_123, fsh0_68, fsh1_68, fsi_93, \
                         fsi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * dsi_37[k]
                   + f_3 * pc_y[k] * fsi_93[k];

        t_122[k] = f_8 * fsh0_68[k]
                   - f_9 * fsh1_68[k]
                   + f_3 * pc_z[k] * fsi_93[k];

        t_123[k] = pa_x[k] * dsk0_123[k]
                   + f_14 * dsi_99[k]
                   - f_12 * pc_x[k] * dsk1_123[k];

        t_124[k] = f_3 * pc_z[k] * fsi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, dsi_42, fsh0_69, fsh0_70, \
                         fsh0_72, fsh1_69, fsh1_70, fsh1_72, fsi_95, fsi_96, \
                         fsi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * fsh0_69[k]
                   - f_5 * fsh1_69[k]
                   + f_3 * pc_z[k] * fsi_95[k];

        t_126[k] = f_6 * fsh0_70[k]
                   - f_7 * fsh1_70[k]
                   + f_3 * pc_z[k] * fsi_96[k];

        t_127[k] = f_14 * dsi_42[k]
                   + f_3 * pc_y[k] * fsi_98[k];

        t_128[k] = f_10 * fsh0_72[k]
                   - f_11 * fsh1_72[k]
                   + f_3 * pc_z[k] * fsi_98[k];
    }
}

static auto
compute_prim_fsk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsk0,
                                                          const size_t dsi, const size_t dsk1,
                                                          const size_t fsh0, const size_t fsh1,
                                                          const size_t fsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_15 = 2.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsk0_39 = buffer.data(dsk0 + 39);
    const auto *dsk0_42 = buffer.data(dsk0 + 42);
    const auto *dsk0_46 = buffer.data(dsk0 + 46);
    const auto *dsk0_51 = buffer.data(dsk0 + 51);
    const auto *dsk0_72 = buffer.data(dsk0 + 72);
    const auto *dsk0_77 = buffer.data(dsk0 + 77);
    const auto *dsk0_81 = buffer.data(dsk0 + 81);
    const auto *dsk0_86 = buffer.data(dsk0 + 86);
    const auto *dsk0_92 = buffer.data(dsk0 + 92);
    const auto *dsk0_108 = buffer.data(dsk0 + 108);
    const auto *dsk0_109 = buffer.data(dsk0 + 109);
    const auto *dsk0_111 = buffer.data(dsk0 + 111);
    const auto *dsk0_136 = buffer.data(dsk0 + 136);
    const auto *dsk0_138 = buffer.data(dsk0 + 138);
    const auto *dsk0_139 = buffer.data(dsk0 + 139);
    const auto *dsk0_140 = buffer.data(dsk0 + 140);
    const auto *dsk0_141 = buffer.data(dsk0 + 141);
    const auto *dsk0_143 = buffer.data(dsk0 + 143);
    const auto *dsk0_156 = buffer.data(dsk0 + 156);
    const auto *dsk0_161 = buffer.data(dsk0 + 161);
    const auto *dsk0_162 = buffer.data(dsk0 + 162);
    const auto *dsk0_172 = buffer.data(dsk0 + 172);
    const auto *dsk0_174 = buffer.data(dsk0 + 174);
    const auto *dsk0_175 = buffer.data(dsk0 + 175);
    const auto *dsk0_176 = buffer.data(dsk0 + 176);
    const auto *dsk0_177 = buffer.data(dsk0 + 177);
    const auto *dsk0_179 = buffer.data(dsk0 + 179);
    const auto *dsk0_180 = buffer.data(dsk0 + 180);
    const auto *dsk0_185 = buffer.data(dsk0 + 185);
    const auto *dsk0_189 = buffer.data(dsk0 + 189);
    const auto *dsk0_194 = buffer.data(dsk0 + 194);
    const auto *dsk0_200 = buffer.data(dsk0 + 200);
    const auto *dsk0_208 = buffer.data(dsk0 + 208);
    const auto *dsk0_209 = buffer.data(dsk0 + 209);
    const auto *dsk0_210 = buffer.data(dsk0 + 210);
    const auto *dsk0_211 = buffer.data(dsk0 + 211);
    const auto *dsk0_212 = buffer.data(dsk0 + 212);
    const auto *dsk0_213 = buffer.data(dsk0 + 213);
    const auto *dsk0_215 = buffer.data(dsk0 + 215);

    const auto *dsi_28 = buffer.data(dsi + 28);
    const auto *dsi_31 = buffer.data(dsi + 31);
    const auto *dsi_34 = buffer.data(dsi + 34);
    const auto *dsi_38 = buffer.data(dsi + 38);
    const auto *dsi_49 = buffer.data(dsi + 49);
    const auto *dsi_55 = buffer.data(dsi + 55);
    const auto *dsi_56 = buffer.data(dsi + 56);
    const auto *dsi_58 = buffer.data(dsi + 58);
    const auto *dsi_61 = buffer.data(dsi + 61);
    const auto *dsi_65 = buffer.data(dsi + 65);
    const auto *dsi_70 = buffer.data(dsi + 70);
    const auto *dsi_83 = buffer.data(dsi + 83);
    const auto *dsi_105 = buffer.data(dsi + 105);
    const auto *dsi_107 = buffer.data(dsi + 107);
    const auto *dsi_108 = buffer.data(dsi + 108);
    const auto *dsi_109 = buffer.data(dsi + 109);
    const auto *dsi_110 = buffer.data(dsi + 110);
    const auto *dsi_111 = buffer.data(dsi + 111);
    const auto *dsi_124 = buffer.data(dsi + 124);
    const auto *dsi_129 = buffer.data(dsi + 129);
    const auto *dsi_130 = buffer.data(dsi + 130);
    const auto *dsi_133 = buffer.data(dsi + 133);
    const auto *dsi_134 = buffer.data(dsi + 134);
    const auto *dsi_135 = buffer.data(dsi + 135);
    const auto *dsi_136 = buffer.data(dsi + 136);
    const auto *dsi_137 = buffer.data(dsi + 137);
    const auto *dsi_138 = buffer.data(dsi + 138);
    const auto *dsi_139 = buffer.data(dsi + 139);
    const auto *dsi_140 = buffer.data(dsi + 140);
    const auto *dsi_145 = buffer.data(dsi + 145);
    const auto *dsi_149 = buffer.data(dsi + 149);
    const auto *dsi_154 = buffer.data(dsi + 154);
    const auto *dsi_160 = buffer.data(dsi + 160);
    const auto *dsi_161 = buffer.data(dsi + 161);
    const auto *dsi_162 = buffer.data(dsi + 162);
    const auto *dsi_163 = buffer.data(dsi + 163);
    const auto *dsi_164 = buffer.data(dsi + 164);
    const auto *dsi_165 = buffer.data(dsi + 165);
    const auto *dsi_167 = buffer.data(dsi + 167);

    const auto *dsk1_39 = buffer.data(dsk1 + 39);
    const auto *dsk1_42 = buffer.data(dsk1 + 42);
    const auto *dsk1_46 = buffer.data(dsk1 + 46);
    const auto *dsk1_51 = buffer.data(dsk1 + 51);
    const auto *dsk1_72 = buffer.data(dsk1 + 72);
    const auto *dsk1_77 = buffer.data(dsk1 + 77);
    const auto *dsk1_81 = buffer.data(dsk1 + 81);
    const auto *dsk1_86 = buffer.data(dsk1 + 86);
    const auto *dsk1_92 = buffer.data(dsk1 + 92);
    const auto *dsk1_108 = buffer.data(dsk1 + 108);
    const auto *dsk1_109 = buffer.data(dsk1 + 109);
    const auto *dsk1_111 = buffer.data(dsk1 + 111);
    const auto *dsk1_136 = buffer.data(dsk1 + 136);
    const auto *dsk1_138 = buffer.data(dsk1 + 138);
    const auto *dsk1_139 = buffer.data(dsk1 + 139);
    const auto *dsk1_140 = buffer.data(dsk1 + 140);
    const auto *dsk1_141 = buffer.data(dsk1 + 141);
    const auto *dsk1_143 = buffer.data(dsk1 + 143);
    const auto *dsk1_156 = buffer.data(dsk1 + 156);
    const auto *dsk1_161 = buffer.data(dsk1 + 161);
    const auto *dsk1_162 = buffer.data(dsk1 + 162);
    const auto *dsk1_172 = buffer.data(dsk1 + 172);
    const auto *dsk1_174 = buffer.data(dsk1 + 174);
    const auto *dsk1_175 = buffer.data(dsk1 + 175);
    const auto *dsk1_176 = buffer.data(dsk1 + 176);
    const auto *dsk1_177 = buffer.data(dsk1 + 177);
    const auto *dsk1_179 = buffer.data(dsk1 + 179);
    const auto *dsk1_180 = buffer.data(dsk1 + 180);
    const auto *dsk1_185 = buffer.data(dsk1 + 185);
    const auto *dsk1_189 = buffer.data(dsk1 + 189);
    const auto *dsk1_194 = buffer.data(dsk1 + 194);
    const auto *dsk1_200 = buffer.data(dsk1 + 200);
    const auto *dsk1_208 = buffer.data(dsk1 + 208);
    const auto *dsk1_209 = buffer.data(dsk1 + 209);
    const auto *dsk1_210 = buffer.data(dsk1 + 210);
    const auto *dsk1_211 = buffer.data(dsk1 + 211);
    const auto *dsk1_212 = buffer.data(dsk1 + 212);
    const auto *dsk1_213 = buffer.data(dsk1 + 213);
    const auto *dsk1_215 = buffer.data(dsk1 + 215);

    const auto *fsh0_105 = buffer.data(fsh0 + 105);
    const auto *fsh0_106 = buffer.data(fsh0 + 106);
    const auto *fsh0_107 = buffer.data(fsh0 + 107);
    const auto *fsh0_108 = buffer.data(fsh0 + 108);
    const auto *fsh0_109 = buffer.data(fsh0 + 109);
    const auto *fsh0_110 = buffer.data(fsh0 + 110);
    const auto *fsh0_111 = buffer.data(fsh0 + 111);
    const auto *fsh0_112 = buffer.data(fsh0 + 112);
    const auto *fsh0_113 = buffer.data(fsh0 + 113);
    const auto *fsh0_114 = buffer.data(fsh0 + 114);
    const auto *fsh0_126 = buffer.data(fsh0 + 126);
    const auto *fsh0_127 = buffer.data(fsh0 + 127);
    const auto *fsh0_129 = buffer.data(fsh0 + 129);
    const auto *fsh0_131 = buffer.data(fsh0 + 131);
    const auto *fsh0_132 = buffer.data(fsh0 + 132);
    const auto *fsh0_134 = buffer.data(fsh0 + 134);
    const auto *fsh0_135 = buffer.data(fsh0 + 135);
    const auto *fsh0_136 = buffer.data(fsh0 + 136);
    const auto *fsh0_138 = buffer.data(fsh0 + 138);
    const auto *fsh0_139 = buffer.data(fsh0 + 139);
    const auto *fsh0_140 = buffer.data(fsh0 + 140);
    const auto *fsh0_141 = buffer.data(fsh0 + 141);
    const auto *fsh0_142 = buffer.data(fsh0 + 142);
    const auto *fsh0_143 = buffer.data(fsh0 + 143);
    const auto *fsh0_144 = buffer.data(fsh0 + 144);
    const auto *fsh0_145 = buffer.data(fsh0 + 145);
    const auto *fsh0_146 = buffer.data(fsh0 + 146);
    const auto *fsh0_149 = buffer.data(fsh0 + 149);
    const auto *fsh0_151 = buffer.data(fsh0 + 151);

    const auto *fsh1_105 = buffer.data(fsh1 + 105);
    const auto *fsh1_106 = buffer.data(fsh1 + 106);
    const auto *fsh1_107 = buffer.data(fsh1 + 107);
    const auto *fsh1_108 = buffer.data(fsh1 + 108);
    const auto *fsh1_109 = buffer.data(fsh1 + 109);
    const auto *fsh1_110 = buffer.data(fsh1 + 110);
    const auto *fsh1_111 = buffer.data(fsh1 + 111);
    const auto *fsh1_112 = buffer.data(fsh1 + 112);
    const auto *fsh1_113 = buffer.data(fsh1 + 113);
    const auto *fsh1_114 = buffer.data(fsh1 + 114);
    const auto *fsh1_126 = buffer.data(fsh1 + 126);
    const auto *fsh1_127 = buffer.data(fsh1 + 127);
    const auto *fsh1_129 = buffer.data(fsh1 + 129);
    const auto *fsh1_131 = buffer.data(fsh1 + 131);
    const auto *fsh1_132 = buffer.data(fsh1 + 132);
    const auto *fsh1_134 = buffer.data(fsh1 + 134);
    const auto *fsh1_135 = buffer.data(fsh1 + 135);
    const auto *fsh1_136 = buffer.data(fsh1 + 136);
    const auto *fsh1_138 = buffer.data(fsh1 + 138);
    const auto *fsh1_139 = buffer.data(fsh1 + 139);
    const auto *fsh1_140 = buffer.data(fsh1 + 140);
    const auto *fsh1_141 = buffer.data(fsh1 + 141);
    const auto *fsh1_142 = buffer.data(fsh1 + 142);
    const auto *fsh1_143 = buffer.data(fsh1 + 143);
    const auto *fsh1_144 = buffer.data(fsh1 + 144);
    const auto *fsh1_145 = buffer.data(fsh1 + 145);
    const auto *fsh1_146 = buffer.data(fsh1 + 146);
    const auto *fsh1_149 = buffer.data(fsh1 + 149);
    const auto *fsh1_151 = buffer.data(fsh1 + 151);

    const auto *fsi_99 = buffer.data(fsi + 99);
    const auto *fsi_105 = buffer.data(fsi + 105);
    const auto *fsi_107 = buffer.data(fsi + 107);
    const auto *fsi_108 = buffer.data(fsi + 108);
    const auto *fsi_109 = buffer.data(fsi + 109);
    const auto *fsi_110 = buffer.data(fsi + 110);
    const auto *fsi_111 = buffer.data(fsi + 111);
    const auto *fsi_112 = buffer.data(fsi + 112);
    const auto *fsi_114 = buffer.data(fsi + 114);
    const auto *fsi_115 = buffer.data(fsi + 115);
    const auto *fsi_117 = buffer.data(fsi + 117);
    const auto *fsi_118 = buffer.data(fsi + 118);
    const auto *fsi_121 = buffer.data(fsi + 121);
    const auto *fsi_122 = buffer.data(fsi + 122);
    const auto *fsi_126 = buffer.data(fsi + 126);
    const auto *fsi_133 = buffer.data(fsi + 133);
    const auto *fsi_134 = buffer.data(fsi + 134);
    const auto *fsi_135 = buffer.data(fsi + 135);
    const auto *fsi_136 = buffer.data(fsi + 136);
    const auto *fsi_137 = buffer.data(fsi + 137);
    const auto *fsi_138 = buffer.data(fsi + 138);
    const auto *fsi_139 = buffer.data(fsi + 139);
    const auto *fsi_140 = buffer.data(fsi + 140);
    const auto *fsi_141 = buffer.data(fsi + 141);
    const auto *fsi_142 = buffer.data(fsi + 142);
    const auto *fsi_143 = buffer.data(fsi + 143);
    const auto *fsi_144 = buffer.data(fsi + 144);
    const auto *fsi_145 = buffer.data(fsi + 145);
    const auto *fsi_146 = buffer.data(fsi + 146);
    const auto *fsi_147 = buffer.data(fsi + 147);
    const auto *fsi_148 = buffer.data(fsi + 148);
    const auto *fsi_149 = buffer.data(fsi + 149);
    const auto *fsi_150 = buffer.data(fsi + 150);
    const auto *fsi_151 = buffer.data(fsi + 151);
    const auto *fsi_152 = buffer.data(fsi + 152);
    const auto *fsi_153 = buffer.data(fsi + 153);
    const auto *fsi_154 = buffer.data(fsi + 154);
    const auto *fsi_160 = buffer.data(fsi + 160);
    const auto *fsi_161 = buffer.data(fsi + 161);
    const auto *fsi_162 = buffer.data(fsi + 162);
    const auto *fsi_163 = buffer.data(fsi + 163);
    const auto *fsi_164 = buffer.data(fsi + 164);
    const auto *fsi_165 = buffer.data(fsi + 165);
    const auto *fsi_167 = buffer.data(fsi + 167);
    const auto *fsi_168 = buffer.data(fsi + 168);
    const auto *fsi_169 = buffer.data(fsi + 169);
    const auto *fsi_171 = buffer.data(fsi + 171);
    const auto *fsi_173 = buffer.data(fsi + 173);
    const auto *fsi_174 = buffer.data(fsi + 174);
    const auto *fsi_176 = buffer.data(fsi + 176);
    const auto *fsi_177 = buffer.data(fsi + 177);
    const auto *fsi_178 = buffer.data(fsi + 178);
    const auto *fsi_180 = buffer.data(fsi + 180);
    const auto *fsi_181 = buffer.data(fsi + 181);
    const auto *fsi_182 = buffer.data(fsi + 182);
    const auto *fsi_183 = buffer.data(fsi + 183);
    const auto *fsi_185 = buffer.data(fsi + 185);
    const auto *fsi_186 = buffer.data(fsi + 186);
    const auto *fsi_187 = buffer.data(fsi + 187);
    const auto *fsi_188 = buffer.data(fsi + 188);
    const auto *fsi_189 = buffer.data(fsi + 189);
    const auto *fsi_190 = buffer.data(fsi + 190);
    const auto *fsi_191 = buffer.data(fsi + 191);
    const auto *fsi_192 = buffer.data(fsi + 192);
    const auto *fsi_193 = buffer.data(fsi + 193);
    const auto *fsi_194 = buffer.data(fsi + 194);
    const auto *fsi_195 = buffer.data(fsi + 195);
    const auto *fsi_198 = buffer.data(fsi + 198);
    const auto *fsi_200 = buffer.data(fsi + 200);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, dsi_105, dsi_107, \
                         dsi_108, dsi_109, fsi_99, fsi_105, fsi_107, fsi_108, \
                         fsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_13 * dsi_105[k]
                   + f_3 * pc_x[k] * fsi_105[k];

        t_130[k] = f_3 * pc_z[k] * fsi_99[k];

        t_131[k] = f_13 * dsi_107[k]
                   + f_3 * pc_x[k] * fsi_107[k];

        t_132[k] = f_13 * dsi_108[k]
                   + f_3 * pc_x[k] * fsi_108[k];

        t_133[k] = f_13 * dsi_109[k]
                   + f_3 * pc_x[k] * fsi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pc_x, pc_z, dsk0_136, dsi_110, \
                         dsi_111, dsk1_136, fsi_105, fsi_110, fsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_13 * dsi_110[k]
                   + f_3 * pc_x[k] * fsi_110[k];

        t_135[k] = f_13 * dsi_111[k]
                   + f_3 * pc_x[k] * fsi_111[k];

        t_136[k] = pa_x[k] * dsk0_136[k]
                   - f_12 * pc_x[k] * dsk1_136[k];

        t_137[k] = f_3 * pc_z[k] * fsi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_x, pc_x, dsk0_138, dsk0_139, dsk0_140, \
                         dsk0_141, dsk1_138, dsk1_139, dsk1_140, \
                         dsk1_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pa_x[k] * dsk0_138[k]
                   - f_12 * pc_x[k] * dsk1_138[k];

        t_139[k] = pa_x[k] * dsk0_139[k]
                   - f_12 * pc_x[k] * dsk1_139[k];

        t_140[k] = pa_x[k] * dsk0_140[k]
                   - f_12 * pc_x[k] * dsk1_140[k];

        t_141[k] = pa_x[k] * dsk0_141[k]
                   - f_12 * pc_x[k] * dsk1_141[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_x, pa_y, pc_x, pc_y, dsk0_72, \
                         dsk0_143, dsi_55, dsi_56, dsk1_72, dsk1_143, fsi_111, \
                         fsi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * dsi_55[k]
                   + f_3 * pc_y[k] * fsi_111[k];

        t_143[k] = pa_x[k] * dsk0_143[k]
                   - f_12 * pc_x[k] * dsk1_143[k];

        t_144[k] = pa_y[k] * dsk0_72[k]
                   - f_12 * pc_y[k] * dsk1_72[k];

        t_145[k] = f_13 * dsi_56[k]
                   + f_3 * pc_y[k] * fsi_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_y, pa_z, pc_y, pc_z, dsk0_39, dsk0_77, \
                         dsi_28, dsi_58, dsk1_39, dsk1_77, fsi_112, \
                         fsi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * dsi_28[k]
                   + f_3 * pc_z[k] * fsi_112[k];

        t_147[k] = pa_z[k] * dsk0_39[k]
                   - f_12 * pc_z[k] * dsk1_39[k];

        t_148[k] = f_13 * dsi_58[k]
                   + f_3 * pc_y[k] * fsi_114[k];

        t_149[k] = pa_y[k] * dsk0_77[k]
                   - f_12 * pc_y[k] * dsk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pa_z, pc_y, pc_z, dsk0_42, dsk0_81, \
                         dsi_31, dsi_61, dsk1_42, dsk1_81, fsi_115, \
                         fsi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * dsk0_42[k]
                   - f_12 * pc_z[k] * dsk1_42[k];

        t_151[k] = f_13 * dsi_31[k]
                   + f_3 * pc_z[k] * fsi_115[k];

        t_152[k] = f_13 * dsi_61[k]
                   + f_3 * pc_y[k] * fsi_117[k];

        t_153[k] = pa_y[k] * dsk0_81[k]
                   - f_12 * pc_y[k] * dsk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_x, pa_z, pc_x, pc_z, dsk0_46, dsk0_156, \
                         dsi_34, dsi_124, dsk1_46, dsk1_156, fsi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_z[k] * dsk0_46[k]
                   - f_12 * pc_z[k] * dsk1_46[k];

        t_155[k] = f_13 * dsi_34[k]
                   + f_3 * pc_z[k] * fsi_118[k];

        t_156[k] = pa_x[k] * dsk0_156[k]
                   + f_0 * dsi_124[k]
                   - f_12 * pc_x[k] * dsk1_156[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pa_y, pa_z, pc_y, pc_z, dsk0_51, dsk0_86, \
                         dsi_38, dsi_65, dsk1_51, dsk1_86, fsi_121, \
                         fsi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * dsi_65[k]
                   + f_3 * pc_y[k] * fsi_121[k];

        t_158[k] = pa_y[k] * dsk0_86[k]
                   - f_12 * pc_y[k] * dsk1_86[k];

        t_159[k] = pa_z[k] * dsk0_51[k]
                   - f_12 * pc_z[k] * dsk1_51[k];

        t_160[k] = f_13 * dsi_38[k]
                   + f_3 * pc_z[k] * fsi_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pa_x, pc_x, pc_y, dsk0_161, dsk0_162, dsi_70, \
                         dsi_129, dsi_130, dsk1_161, dsk1_162, \
                         fsi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pa_x[k] * dsk0_161[k]
                   + f_14 * dsi_129[k]
                   - f_12 * pc_x[k] * dsk1_161[k];

        t_162[k] = pa_x[k] * dsk0_162[k]
                   + f_14 * dsi_130[k]
                   - f_12 * pc_x[k] * dsk1_162[k];

        t_163[k] = f_13 * dsi_70[k]
                   + f_3 * pc_y[k] * fsi_126[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_y, pc_x, pc_y, dsk0_92, dsi_133, \
                         dsi_134, dsi_135, dsk1_92, fsi_133, fsi_134, \
                         fsi_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_y[k] * dsk0_92[k]
                   - f_12 * pc_y[k] * dsk1_92[k];

        t_165[k] = f_13 * dsi_133[k]
                   + f_3 * pc_x[k] * fsi_133[k];

        t_166[k] = f_13 * dsi_134[k]
                   + f_3 * pc_x[k] * fsi_134[k];

        t_167[k] = f_13 * dsi_135[k]
                   + f_3 * pc_x[k] * fsi_135[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, dsi_136, dsi_137, dsi_138, dsi_139, \
                         fsi_136, fsi_137, fsi_138, fsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_13 * dsi_136[k]
                   + f_3 * pc_x[k] * fsi_136[k];

        t_169[k] = f_13 * dsi_137[k]
                   + f_3 * pc_x[k] * fsi_137[k];

        t_170[k] = f_13 * dsi_138[k]
                   + f_3 * pc_x[k] * fsi_138[k];

        t_171[k] = f_13 * dsi_139[k]
                   + f_3 * pc_x[k] * fsi_139[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_x, pc_x, pc_z, dsk0_172, dsk0_174, \
                         dsk0_175, dsi_49, dsk1_172, dsk1_174, dsk1_175, \
                         fsi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_x[k] * dsk0_172[k]
                   - f_12 * pc_x[k] * dsk1_172[k];

        t_173[k] = f_13 * dsi_49[k]
                   + f_3 * pc_z[k] * fsi_133[k];

        t_174[k] = pa_x[k] * dsk0_174[k]
                   - f_12 * pc_x[k] * dsk1_174[k];

        t_175[k] = pa_x[k] * dsk0_175[k]
                   - f_12 * pc_x[k] * dsk1_175[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pc_x, pc_y, dsk0_176, dsk0_177, \
                         dsk0_179, dsi_83, dsk1_176, dsk1_177, dsk1_179, \
                         fsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_x[k] * dsk0_176[k]
                   - f_12 * pc_x[k] * dsk1_176[k];

        t_177[k] = pa_x[k] * dsk0_177[k]
                   - f_12 * pc_x[k] * dsk1_177[k];

        t_178[k] = f_13 * dsi_83[k]
                   + f_3 * pc_y[k] * fsi_139[k];

        t_179[k] = pa_x[k] * dsk0_179[k]
                   - f_12 * pc_x[k] * dsk1_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pc_x, pc_y, pc_z, dsk0_180, dsi_56, \
                         dsi_140, dsk1_180, fsh0_105, fsh1_105, fsi_140, \
                         fsi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_x[k] * dsk0_180[k]
                   + f_19 * dsi_140[k]
                   - f_12 * pc_x[k] * dsk1_180[k];

        t_181[k] = f_3 * pc_y[k] * fsi_140[k];

        t_182[k] = f_14 * dsi_56[k]
                   + f_3 * pc_z[k] * fsi_140[k];

        t_183[k] = f_4 * fsh0_105[k]
                   - f_5 * fsh1_105[k]
                   + f_3 * pc_y[k] * fsi_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_x, pc_x, pc_y, dsk0_185, dsi_145, dsk1_185, \
                         fsh0_106, fsh1_106, fsi_142, fsi_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * fsi_142[k];

        t_185[k] = pa_x[k] * dsk0_185[k]
                   + f_16 * dsi_145[k]
                   - f_12 * pc_x[k] * dsk1_185[k];

        t_186[k] = f_6 * fsh0_106[k]
                   - f_7 * fsh1_106[k]
                   + f_3 * pc_y[k] * fsi_143[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_x, pc_x, pc_y, dsk0_189, dsi_149, dsk1_189, \
                         fsh0_107, fsh1_107, fsi_144, fsi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_4 * fsh0_107[k]
                   - f_5 * fsh1_107[k]
                   + f_3 * pc_y[k] * fsi_144[k];

        t_188[k] = f_3 * pc_y[k] * fsi_145[k];

        t_189[k] = pa_x[k] * dsk0_189[k]
                   + f_15 * dsi_149[k]
                   - f_12 * pc_x[k] * dsk1_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pc_y, fsh0_108, fsh0_109, fsh0_110, \
                         fsh1_108, fsh1_109, fsh1_110, fsi_146, fsi_147, fsi_148, \
                         fsi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * fsh0_108[k]
                   - f_9 * fsh1_108[k]
                   + f_3 * pc_y[k] * fsi_146[k];

        t_191[k] = f_6 * fsh0_109[k]
                   - f_7 * fsh1_109[k]
                   + f_3 * pc_y[k] * fsi_147[k];

        t_192[k] = f_4 * fsh0_110[k]
                   - f_5 * fsh1_110[k]
                   + f_3 * pc_y[k] * fsi_148[k];

        t_193[k] = f_3 * pc_y[k] * fsi_149[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_x, pc_x, pc_y, dsk0_194, dsi_154, dsk1_194, \
                         fsh0_111, fsh0_112, fsh1_111, fsh1_112, fsi_150, \
                         fsi_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * dsk0_194[k]
                   + f_0 * dsi_154[k]
                   - f_12 * pc_x[k] * dsk1_194[k];

        t_195[k] = f_10 * fsh0_111[k]
                   - f_11 * fsh1_111[k]
                   + f_3 * pc_y[k] * fsi_150[k];

        t_196[k] = f_8 * fsh0_112[k]
                   - f_9 * fsh1_112[k]
                   + f_3 * pc_y[k] * fsi_151[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pc_y, fsh0_113, fsh0_114, fsh1_113, fsh1_114, \
                         fsi_152, fsi_153, fsi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_6 * fsh0_113[k]
                   - f_7 * fsh1_113[k]
                   + f_3 * pc_y[k] * fsi_152[k];

        t_198[k] = f_4 * fsh0_114[k]
                   - f_5 * fsh1_114[k]
                   + f_3 * pc_y[k] * fsi_153[k];

        t_199[k] = f_3 * pc_y[k] * fsi_154[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_x, pc_x, dsk0_200, dsi_160, dsi_161, \
                         dsi_162, dsi_163, dsk1_200, fsi_161, fsi_162, \
                         fsi_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_x[k] * dsk0_200[k]
                   + f_14 * dsi_160[k]
                   - f_12 * pc_x[k] * dsk1_200[k];

        t_201[k] = f_13 * dsi_161[k]
                   + f_3 * pc_x[k] * fsi_161[k];

        t_202[k] = f_13 * dsi_162[k]
                   + f_3 * pc_x[k] * fsi_162[k];

        t_203[k] = f_13 * dsi_163[k]
                   + f_3 * pc_x[k] * fsi_163[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_x, pc_y, dsi_164, dsi_165, dsi_167, \
                         fsi_160, fsi_164, fsi_165, fsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_13 * dsi_164[k]
                   + f_3 * pc_x[k] * fsi_164[k];

        t_205[k] = f_13 * dsi_165[k]
                   + f_3 * pc_x[k] * fsi_165[k];

        t_206[k] = f_3 * pc_y[k] * fsi_160[k];

        t_207[k] = f_13 * dsi_167[k]
                   + f_3 * pc_x[k] * fsi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pc_x, dsk0_208, dsk0_209, dsk0_210, \
                         dsk0_211, dsk1_208, dsk1_209, dsk1_210, \
                         dsk1_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_x[k] * dsk0_208[k]
                   - f_12 * pc_x[k] * dsk1_208[k];

        t_209[k] = pa_x[k] * dsk0_209[k]
                   - f_12 * pc_x[k] * dsk1_209[k];

        t_210[k] = pa_x[k] * dsk0_210[k]
                   - f_12 * pc_x[k] * dsk1_210[k];

        t_211[k] = pa_x[k] * dsk0_211[k]
                   - f_12 * pc_x[k] * dsk1_211[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pc_x, pc_y, dsk0_212, dsk0_213, \
                         dsk0_215, dsk1_212, dsk1_213, dsk1_215, \
                         fsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pa_x[k] * dsk0_212[k]
                   - f_12 * pc_x[k] * dsk1_212[k];

        t_213[k] = pa_x[k] * dsk0_213[k]
                   - f_12 * pc_x[k] * dsk1_213[k];

        t_214[k] = f_3 * pc_y[k] * fsi_167[k];

        t_215[k] = pa_x[k] * dsk0_215[k]
                   - f_12 * pc_x[k] * dsk1_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pc_x, pc_z, fsh0_126, fsh0_127, \
                         fsh0_129, fsh1_126, fsh1_127, fsh1_129, fsi_168, fsi_169, \
                         fsi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_1 * fsh0_126[k]
                   - f_2 * fsh1_126[k]
                   + f_3 * pc_x[k] * fsi_168[k];

        t_217[k] = f_17 * fsh0_127[k]
                   - f_18 * fsh1_127[k]
                   + f_3 * pc_x[k] * fsi_169[k];

        t_218[k] = f_3 * pc_z[k] * fsi_168[k];

        t_219[k] = f_10 * fsh0_129[k]
                   - f_11 * fsh1_129[k]
                   + f_3 * pc_x[k] * fsi_171[k];

        t_220[k] = f_3 * pc_z[k] * fsi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pc_x, pc_z, fsh0_131, fsh0_132, fsh0_134, \
                         fsh1_131, fsh1_132, fsh1_134, fsi_171, fsi_173, fsi_174, \
                         fsi_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_10 * fsh0_131[k]
                   - f_11 * fsh1_131[k]
                   + f_3 * pc_x[k] * fsi_173[k];

        t_222[k] = f_8 * fsh0_132[k]
                   - f_9 * fsh1_132[k]
                   + f_3 * pc_x[k] * fsi_174[k];

        t_223[k] = f_3 * pc_z[k] * fsi_171[k];

        t_224[k] = f_8 * fsh0_134[k]
                   - f_9 * fsh1_134[k]
                   + f_3 * pc_x[k] * fsi_176[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pc_x, pc_z, fsh0_135, fsh0_136, fsh0_138, \
                         fsh1_135, fsh1_136, fsh1_138, fsi_174, fsi_177, fsi_178, \
                         fsi_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_8 * fsh0_135[k]
                   - f_9 * fsh1_135[k]
                   + f_3 * pc_x[k] * fsi_177[k];

        t_226[k] = f_6 * fsh0_136[k]
                   - f_7 * fsh1_136[k]
                   + f_3 * pc_x[k] * fsi_178[k];

        t_227[k] = f_3 * pc_z[k] * fsi_174[k];

        t_228[k] = f_6 * fsh0_138[k]
                   - f_7 * fsh1_138[k]
                   + f_3 * pc_x[k] * fsi_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pc_x, pc_z, fsh0_139, fsh0_140, fsh0_141, \
                         fsh1_139, fsh1_140, fsh1_141, fsi_178, fsi_181, fsi_182, \
                         fsi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_6 * fsh0_139[k]
                   - f_7 * fsh1_139[k]
                   + f_3 * pc_x[k] * fsi_181[k];

        t_230[k] = f_6 * fsh0_140[k]
                   - f_7 * fsh1_140[k]
                   + f_3 * pc_x[k] * fsi_182[k];

        t_231[k] = f_4 * fsh0_141[k]
                   - f_5 * fsh1_141[k]
                   + f_3 * pc_x[k] * fsi_183[k];

        t_232[k] = f_3 * pc_z[k] * fsi_178[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, fsh0_143, fsh0_144, fsh0_145, fsh1_143, \
                         fsh1_144, fsh1_145, fsi_185, fsi_186, \
                         fsi_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_4 * fsh0_143[k]
                   - f_5 * fsh1_143[k]
                   + f_3 * pc_x[k] * fsi_185[k];

        t_234[k] = f_4 * fsh0_144[k]
                   - f_5 * fsh1_144[k]
                   + f_3 * pc_x[k] * fsi_186[k];

        t_235[k] = f_4 * fsh0_145[k]
                   - f_5 * fsh1_145[k]
                   + f_3 * pc_x[k] * fsi_187[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, pc_x, fsh0_146, fsh1_146, \
                         fsi_188, fsi_189, fsi_190, fsi_191, fsi_192, \
                         fsi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_4 * fsh0_146[k]
                   - f_5 * fsh1_146[k]
                   + f_3 * pc_x[k] * fsi_188[k];

        t_237[k] = f_3 * pc_x[k] * fsi_189[k];

        t_238[k] = f_3 * pc_x[k] * fsi_190[k];

        t_239[k] = f_3 * pc_x[k] * fsi_191[k];

        t_240[k] = f_3 * pc_x[k] * fsi_192[k];

        t_241[k] = f_3 * pc_x[k] * fsi_193[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, dsi_105, \
                         fsh0_141, fsh1_141, fsi_189, fsi_190, fsi_194, \
                         fsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * pc_x[k] * fsi_194[k];

        t_243[k] = f_3 * pc_x[k] * fsi_195[k];

        t_244[k] = f_0 * dsi_105[k]
                   + f_1 * fsh0_141[k]
                   - f_2 * fsh1_141[k]
                   + f_3 * pc_y[k] * fsi_189[k];

        t_245[k] = f_3 * pc_z[k] * fsi_189[k];

        t_246[k] = f_4 * fsh0_141[k]
                   - f_5 * fsh1_141[k]
                   + f_3 * pc_z[k] * fsi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, fsh0_142, fsh0_143, fsh0_144, fsh1_142, \
                         fsh1_143, fsh1_144, fsi_191, fsi_192, \
                         fsi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * fsh0_142[k]
                   - f_7 * fsh1_142[k]
                   + f_3 * pc_z[k] * fsi_191[k];

        t_248[k] = f_8 * fsh0_143[k]
                   - f_9 * fsh1_143[k]
                   + f_3 * pc_z[k] * fsi_192[k];

        t_249[k] = f_10 * fsh0_144[k]
                   - f_11 * fsh1_144[k]
                   + f_3 * pc_z[k] * fsi_193[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, dsk0_108, dsk0_109, \
                         dsi_111, dsk1_108, dsk1_109, fsh0_146, fsh1_146, \
                         fsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * dsi_111[k]
                   + f_3 * pc_y[k] * fsi_195[k];

        t_251[k] = f_1 * fsh0_146[k]
                   - f_2 * fsh1_146[k]
                   + f_3 * pc_z[k] * fsi_195[k];

        t_252[k] = pa_z[k] * dsk0_108[k]
                   - f_12 * pc_z[k] * dsk1_108[k];

        t_253[k] = pa_z[k] * dsk0_109[k]
                   - f_12 * pc_z[k] * dsk1_109[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_x, pc_z, dsk0_111, dsk1_111, fsh0_149, \
                         fsh0_151, fsh1_149, fsh1_151, fsi_198, \
                         fsi_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_17 * fsh0_149[k]
                   - f_18 * fsh1_149[k]
                   + f_3 * pc_x[k] * fsi_198[k];

        t_255[k] = pa_z[k] * dsk0_111[k]
                   - f_12 * pc_z[k] * dsk1_111[k];

        t_256[k] = f_10 * fsh0_151[k]
                   - f_11 * fsh1_151[k]
                   + f_3 * pc_x[k] * fsi_200[k];
    }
}

static auto
compute_prim_fsk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsk0,
                                                          const size_t dsi, const size_t dsk1,
                                                          const size_t fsh0, const size_t fsh1,
                                                          const size_t fsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_15 = 2.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsk0_114 = buffer.data(dsk0 + 114);
    const auto *dsk0_118 = buffer.data(dsk0 + 118);
    const auto *dsk0_123 = buffer.data(dsk0 + 123);
    const auto *dsk0_136 = buffer.data(dsk0 + 136);
    const auto *dsk0_138 = buffer.data(dsk0 + 138);
    const auto *dsk0_139 = buffer.data(dsk0 + 139);
    const auto *dsk0_140 = buffer.data(dsk0 + 140);
    const auto *dsk0_141 = buffer.data(dsk0 + 141);
    const auto *dsk0_180 = buffer.data(dsk0 + 180);
    const auto *dsk0_182 = buffer.data(dsk0 + 182);
    const auto *dsk0_185 = buffer.data(dsk0 + 185);
    const auto *dsk0_189 = buffer.data(dsk0 + 189);
    const auto *dsk0_194 = buffer.data(dsk0 + 194);
    const auto *dsk0_200 = buffer.data(dsk0 + 200);
    const auto *dsk0_208 = buffer.data(dsk0 + 208);
    const auto *dsk0_210 = buffer.data(dsk0 + 210);
    const auto *dsk0_211 = buffer.data(dsk0 + 211);
    const auto *dsk0_212 = buffer.data(dsk0 + 212);
    const auto *dsk0_213 = buffer.data(dsk0 + 213);
    const auto *dsk0_215 = buffer.data(dsk0 + 215);

    const auto *dsi_105 = buffer.data(dsi + 105);
    const auto *dsi_106 = buffer.data(dsi + 106);
    const auto *dsi_107 = buffer.data(dsi + 107);
    const auto *dsi_108 = buffer.data(dsi + 108);
    const auto *dsi_109 = buffer.data(dsi + 109);
    const auto *dsi_111 = buffer.data(dsi + 111);
    const auto *dsi_133 = buffer.data(dsi + 133);
    const auto *dsi_139 = buffer.data(dsi + 139);
    const auto *dsi_161 = buffer.data(dsi + 161);
    const auto *dsi_163 = buffer.data(dsi + 163);
    const auto *dsi_164 = buffer.data(dsi + 164);
    const auto *dsi_165 = buffer.data(dsi + 165);
    const auto *dsi_166 = buffer.data(dsi + 166);
    const auto *dsi_167 = buffer.data(dsi + 167);

    const auto *dsk1_114 = buffer.data(dsk1 + 114);
    const auto *dsk1_118 = buffer.data(dsk1 + 118);
    const auto *dsk1_123 = buffer.data(dsk1 + 123);
    const auto *dsk1_136 = buffer.data(dsk1 + 136);
    const auto *dsk1_138 = buffer.data(dsk1 + 138);
    const auto *dsk1_139 = buffer.data(dsk1 + 139);
    const auto *dsk1_140 = buffer.data(dsk1 + 140);
    const auto *dsk1_141 = buffer.data(dsk1 + 141);
    const auto *dsk1_180 = buffer.data(dsk1 + 180);
    const auto *dsk1_182 = buffer.data(dsk1 + 182);
    const auto *dsk1_185 = buffer.data(dsk1 + 185);
    const auto *dsk1_189 = buffer.data(dsk1 + 189);
    const auto *dsk1_194 = buffer.data(dsk1 + 194);
    const auto *dsk1_200 = buffer.data(dsk1 + 200);
    const auto *dsk1_208 = buffer.data(dsk1 + 208);
    const auto *dsk1_210 = buffer.data(dsk1 + 210);
    const auto *dsk1_211 = buffer.data(dsk1 + 211);
    const auto *dsk1_212 = buffer.data(dsk1 + 212);
    const auto *dsk1_213 = buffer.data(dsk1 + 213);
    const auto *dsk1_215 = buffer.data(dsk1 + 215);

    const auto *fsh0_152 = buffer.data(fsh0 + 152);
    const auto *fsh0_154 = buffer.data(fsh0 + 154);
    const auto *fsh0_155 = buffer.data(fsh0 + 155);
    const auto *fsh0_156 = buffer.data(fsh0 + 156);
    const auto *fsh0_158 = buffer.data(fsh0 + 158);
    const auto *fsh0_159 = buffer.data(fsh0 + 159);
    const auto *fsh0_160 = buffer.data(fsh0 + 160);
    const auto *fsh0_161 = buffer.data(fsh0 + 161);
    const auto *fsh0_163 = buffer.data(fsh0 + 163);
    const auto *fsh0_164 = buffer.data(fsh0 + 164);
    const auto *fsh0_165 = buffer.data(fsh0 + 165);
    const auto *fsh0_166 = buffer.data(fsh0 + 166);
    const auto *fsh0_167 = buffer.data(fsh0 + 167);
    const auto *fsh0_169 = buffer.data(fsh0 + 169);
    const auto *fsh0_171 = buffer.data(fsh0 + 171);
    const auto *fsh0_172 = buffer.data(fsh0 + 172);
    const auto *fsh0_174 = buffer.data(fsh0 + 174);
    const auto *fsh0_175 = buffer.data(fsh0 + 175);
    const auto *fsh0_176 = buffer.data(fsh0 + 176);
    const auto *fsh0_178 = buffer.data(fsh0 + 178);
    const auto *fsh0_179 = buffer.data(fsh0 + 179);
    const auto *fsh0_180 = buffer.data(fsh0 + 180);
    const auto *fsh0_181 = buffer.data(fsh0 + 181);
    const auto *fsh0_183 = buffer.data(fsh0 + 183);
    const auto *fsh0_184 = buffer.data(fsh0 + 184);
    const auto *fsh0_185 = buffer.data(fsh0 + 185);
    const auto *fsh0_186 = buffer.data(fsh0 + 186);
    const auto *fsh0_187 = buffer.data(fsh0 + 187);
    const auto *fsh0_189 = buffer.data(fsh0 + 189);
    const auto *fsh0_191 = buffer.data(fsh0 + 191);
    const auto *fsh0_192 = buffer.data(fsh0 + 192);
    const auto *fsh0_194 = buffer.data(fsh0 + 194);
    const auto *fsh0_195 = buffer.data(fsh0 + 195);
    const auto *fsh0_196 = buffer.data(fsh0 + 196);
    const auto *fsh0_198 = buffer.data(fsh0 + 198);
    const auto *fsh0_199 = buffer.data(fsh0 + 199);
    const auto *fsh0_200 = buffer.data(fsh0 + 200);
    const auto *fsh0_201 = buffer.data(fsh0 + 201);
    const auto *fsh0_203 = buffer.data(fsh0 + 203);
    const auto *fsh0_204 = buffer.data(fsh0 + 204);
    const auto *fsh0_205 = buffer.data(fsh0 + 205);
    const auto *fsh0_206 = buffer.data(fsh0 + 206);
    const auto *fsh0_207 = buffer.data(fsh0 + 207);
    const auto *fsh0_208 = buffer.data(fsh0 + 208);
    const auto *fsh0_209 = buffer.data(fsh0 + 209);

    const auto *fsh1_152 = buffer.data(fsh1 + 152);
    const auto *fsh1_154 = buffer.data(fsh1 + 154);
    const auto *fsh1_155 = buffer.data(fsh1 + 155);
    const auto *fsh1_156 = buffer.data(fsh1 + 156);
    const auto *fsh1_158 = buffer.data(fsh1 + 158);
    const auto *fsh1_159 = buffer.data(fsh1 + 159);
    const auto *fsh1_160 = buffer.data(fsh1 + 160);
    const auto *fsh1_161 = buffer.data(fsh1 + 161);
    const auto *fsh1_163 = buffer.data(fsh1 + 163);
    const auto *fsh1_164 = buffer.data(fsh1 + 164);
    const auto *fsh1_165 = buffer.data(fsh1 + 165);
    const auto *fsh1_166 = buffer.data(fsh1 + 166);
    const auto *fsh1_167 = buffer.data(fsh1 + 167);
    const auto *fsh1_169 = buffer.data(fsh1 + 169);
    const auto *fsh1_171 = buffer.data(fsh1 + 171);
    const auto *fsh1_172 = buffer.data(fsh1 + 172);
    const auto *fsh1_174 = buffer.data(fsh1 + 174);
    const auto *fsh1_175 = buffer.data(fsh1 + 175);
    const auto *fsh1_176 = buffer.data(fsh1 + 176);
    const auto *fsh1_178 = buffer.data(fsh1 + 178);
    const auto *fsh1_179 = buffer.data(fsh1 + 179);
    const auto *fsh1_180 = buffer.data(fsh1 + 180);
    const auto *fsh1_181 = buffer.data(fsh1 + 181);
    const auto *fsh1_183 = buffer.data(fsh1 + 183);
    const auto *fsh1_184 = buffer.data(fsh1 + 184);
    const auto *fsh1_185 = buffer.data(fsh1 + 185);
    const auto *fsh1_186 = buffer.data(fsh1 + 186);
    const auto *fsh1_187 = buffer.data(fsh1 + 187);
    const auto *fsh1_189 = buffer.data(fsh1 + 189);
    const auto *fsh1_191 = buffer.data(fsh1 + 191);
    const auto *fsh1_192 = buffer.data(fsh1 + 192);
    const auto *fsh1_194 = buffer.data(fsh1 + 194);
    const auto *fsh1_195 = buffer.data(fsh1 + 195);
    const auto *fsh1_196 = buffer.data(fsh1 + 196);
    const auto *fsh1_198 = buffer.data(fsh1 + 198);
    const auto *fsh1_199 = buffer.data(fsh1 + 199);
    const auto *fsh1_200 = buffer.data(fsh1 + 200);
    const auto *fsh1_201 = buffer.data(fsh1 + 201);
    const auto *fsh1_203 = buffer.data(fsh1 + 203);
    const auto *fsh1_204 = buffer.data(fsh1 + 204);
    const auto *fsh1_205 = buffer.data(fsh1 + 205);
    const auto *fsh1_206 = buffer.data(fsh1 + 206);
    const auto *fsh1_207 = buffer.data(fsh1 + 207);
    const auto *fsh1_208 = buffer.data(fsh1 + 208);
    const auto *fsh1_209 = buffer.data(fsh1 + 209);

    const auto *fsi_201 = buffer.data(fsi + 201);
    const auto *fsi_203 = buffer.data(fsi + 203);
    const auto *fsi_204 = buffer.data(fsi + 204);
    const auto *fsi_205 = buffer.data(fsi + 205);
    const auto *fsi_207 = buffer.data(fsi + 207);
    const auto *fsi_208 = buffer.data(fsi + 208);
    const auto *fsi_209 = buffer.data(fsi + 209);
    const auto *fsi_210 = buffer.data(fsi + 210);
    const auto *fsi_212 = buffer.data(fsi + 212);
    const auto *fsi_213 = buffer.data(fsi + 213);
    const auto *fsi_214 = buffer.data(fsi + 214);
    const auto *fsi_215 = buffer.data(fsi + 215);
    const auto *fsi_216 = buffer.data(fsi + 216);
    const auto *fsi_217 = buffer.data(fsi + 217);
    const auto *fsi_218 = buffer.data(fsi + 218);
    const auto *fsi_219 = buffer.data(fsi + 219);
    const auto *fsi_220 = buffer.data(fsi + 220);
    const auto *fsi_221 = buffer.data(fsi + 221);
    const auto *fsi_222 = buffer.data(fsi + 222);
    const auto *fsi_223 = buffer.data(fsi + 223);
    const auto *fsi_225 = buffer.data(fsi + 225);
    const auto *fsi_227 = buffer.data(fsi + 227);
    const auto *fsi_228 = buffer.data(fsi + 228);
    const auto *fsi_230 = buffer.data(fsi + 230);
    const auto *fsi_231 = buffer.data(fsi + 231);
    const auto *fsi_232 = buffer.data(fsi + 232);
    const auto *fsi_234 = buffer.data(fsi + 234);
    const auto *fsi_235 = buffer.data(fsi + 235);
    const auto *fsi_236 = buffer.data(fsi + 236);
    const auto *fsi_237 = buffer.data(fsi + 237);
    const auto *fsi_239 = buffer.data(fsi + 239);
    const auto *fsi_240 = buffer.data(fsi + 240);
    const auto *fsi_241 = buffer.data(fsi + 241);
    const auto *fsi_242 = buffer.data(fsi + 242);
    const auto *fsi_243 = buffer.data(fsi + 243);
    const auto *fsi_245 = buffer.data(fsi + 245);
    const auto *fsi_246 = buffer.data(fsi + 246);
    const auto *fsi_247 = buffer.data(fsi + 247);
    const auto *fsi_248 = buffer.data(fsi + 248);
    const auto *fsi_249 = buffer.data(fsi + 249);
    const auto *fsi_250 = buffer.data(fsi + 250);
    const auto *fsi_251 = buffer.data(fsi + 251);
    const auto *fsi_252 = buffer.data(fsi + 252);
    const auto *fsi_254 = buffer.data(fsi + 254);
    const auto *fsi_255 = buffer.data(fsi + 255);
    const auto *fsi_257 = buffer.data(fsi + 257);
    const auto *fsi_258 = buffer.data(fsi + 258);
    const auto *fsi_259 = buffer.data(fsi + 259);
    const auto *fsi_261 = buffer.data(fsi + 261);
    const auto *fsi_262 = buffer.data(fsi + 262);
    const auto *fsi_263 = buffer.data(fsi + 263);
    const auto *fsi_264 = buffer.data(fsi + 264);
    const auto *fsi_266 = buffer.data(fsi + 266);
    const auto *fsi_267 = buffer.data(fsi + 267);
    const auto *fsi_268 = buffer.data(fsi + 268);
    const auto *fsi_269 = buffer.data(fsi + 269);
    const auto *fsi_270 = buffer.data(fsi + 270);
    const auto *fsi_272 = buffer.data(fsi + 272);
    const auto *fsi_273 = buffer.data(fsi + 273);
    const auto *fsi_274 = buffer.data(fsi + 274);
    const auto *fsi_275 = buffer.data(fsi + 275);
    const auto *fsi_276 = buffer.data(fsi + 276);
    const auto *fsi_277 = buffer.data(fsi + 277);
    const auto *fsi_278 = buffer.data(fsi + 278);
    const auto *fsi_279 = buffer.data(fsi + 279);

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, dsk0_114, dsk1_114, fsh0_152, \
                         fsh0_154, fsh1_152, fsh1_154, fsi_201, \
                         fsi_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_10 * fsh0_152[k]
                   - f_11 * fsh1_152[k]
                   + f_3 * pc_x[k] * fsi_201[k];

        t_258[k] = pa_z[k] * dsk0_114[k]
                   - f_12 * pc_z[k] * dsk1_114[k];

        t_259[k] = f_8 * fsh0_154[k]
                   - f_9 * fsh1_154[k]
                   + f_3 * pc_x[k] * fsi_203[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_z, dsk0_118, dsk1_118, fsh0_155, \
                         fsh0_156, fsh1_155, fsh1_156, fsi_204, \
                         fsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_8 * fsh0_155[k]
                   - f_9 * fsh1_155[k]
                   + f_3 * pc_x[k] * fsi_204[k];

        t_261[k] = f_8 * fsh0_156[k]
                   - f_9 * fsh1_156[k]
                   + f_3 * pc_x[k] * fsi_205[k];

        t_262[k] = pa_z[k] * dsk0_118[k]
                   - f_12 * pc_z[k] * dsk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pc_x, fsh0_158, fsh0_159, fsh0_160, fsh1_158, \
                         fsh1_159, fsh1_160, fsi_207, fsi_208, \
                         fsi_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_6 * fsh0_158[k]
                   - f_7 * fsh1_158[k]
                   + f_3 * pc_x[k] * fsi_207[k];

        t_264[k] = f_6 * fsh0_159[k]
                   - f_7 * fsh1_159[k]
                   + f_3 * pc_x[k] * fsi_208[k];

        t_265[k] = f_6 * fsh0_160[k]
                   - f_7 * fsh1_160[k]
                   + f_3 * pc_x[k] * fsi_209[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, dsk0_123, dsk1_123, fsh0_161, \
                         fsh0_163, fsh1_161, fsh1_163, fsi_210, \
                         fsi_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_6 * fsh0_161[k]
                   - f_7 * fsh1_161[k]
                   + f_3 * pc_x[k] * fsi_210[k];

        t_267[k] = pa_z[k] * dsk0_123[k]
                   - f_12 * pc_z[k] * dsk1_123[k];

        t_268[k] = f_4 * fsh0_163[k]
                   - f_5 * fsh1_163[k]
                   + f_3 * pc_x[k] * fsi_212[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pc_x, fsh0_164, fsh0_165, fsh0_166, fsh1_164, \
                         fsh1_165, fsh1_166, fsi_213, fsi_214, \
                         fsi_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_4 * fsh0_164[k]
                   - f_5 * fsh1_164[k]
                   + f_3 * pc_x[k] * fsi_213[k];

        t_270[k] = f_4 * fsh0_165[k]
                   - f_5 * fsh1_165[k]
                   + f_3 * pc_x[k] * fsi_214[k];

        t_271[k] = f_4 * fsh0_166[k]
                   - f_5 * fsh1_166[k]
                   + f_3 * pc_x[k] * fsi_215[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, t_277, pc_x, fsh0_167, fsh1_167, \
                         fsi_216, fsi_217, fsi_218, fsi_219, fsi_220, \
                         fsi_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_4 * fsh0_167[k]
                   - f_5 * fsh1_167[k]
                   + f_3 * pc_x[k] * fsi_216[k];

        t_273[k] = f_3 * pc_x[k] * fsi_217[k];

        t_274[k] = f_3 * pc_x[k] * fsi_218[k];

        t_275[k] = f_3 * pc_x[k] * fsi_219[k];

        t_276[k] = f_3 * pc_x[k] * fsi_220[k];

        t_277[k] = f_3 * pc_x[k] * fsi_221[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_z, pc_x, pc_z, dsk0_136, dsi_105, \
                         dsk1_136, fsi_217, fsi_222, fsi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_3 * pc_x[k] * fsi_222[k];

        t_279[k] = f_3 * pc_x[k] * fsi_223[k];

        t_280[k] = pa_z[k] * dsk0_136[k]
                   - f_12 * pc_z[k] * dsk1_136[k];

        t_281[k] = f_13 * dsi_105[k]
                   + f_3 * pc_z[k] * fsi_217[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pa_z, pc_z, dsk0_138, dsk0_139, dsk0_140, \
                         dsi_106, dsi_107, dsi_108, dsk1_138, dsk1_139, \
                         dsk1_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pa_z[k] * dsk0_138[k]
                   + f_14 * dsi_106[k]
                   - f_12 * pc_z[k] * dsk1_138[k];

        t_283[k] = pa_z[k] * dsk0_139[k]
                   + f_0 * dsi_107[k]
                   - f_12 * pc_z[k] * dsk1_139[k];

        t_284[k] = pa_z[k] * dsk0_140[k]
                   + f_15 * dsi_108[k]
                   - f_12 * pc_z[k] * dsk1_140[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_z, pc_y, pc_z, dsk0_141, dsi_109, dsi_111, \
                         dsi_139, dsk1_141, fsh0_167, fsh1_167, \
                         fsi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = pa_z[k] * dsk0_141[k]
                   + f_16 * dsi_109[k]
                   - f_12 * pc_z[k] * dsk1_141[k];

        t_286[k] = f_14 * dsi_139[k]
                   + f_3 * pc_y[k] * fsi_223[k];

        t_287[k] = f_13 * dsi_111[k]
                   + f_1 * fsh0_167[k]
                   - f_2 * fsh1_167[k]
                   + f_3 * pc_z[k] * fsi_223[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_y, pc_x, pc_y, dsk0_180, dsk0_182, dsk1_180, \
                         dsk1_182, fsh0_169, fsh1_169, fsi_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pa_y[k] * dsk0_180[k]
                   - f_12 * pc_y[k] * dsk1_180[k];

        t_289[k] = f_17 * fsh0_169[k]
                   - f_18 * fsh1_169[k]
                   + f_3 * pc_x[k] * fsi_225[k];

        t_290[k] = pa_y[k] * dsk0_182[k]
                   - f_12 * pc_y[k] * dsk1_182[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_x, pc_y, dsk0_185, dsk1_185, fsh0_171, \
                         fsh0_172, fsh1_171, fsh1_172, fsi_227, \
                         fsi_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_10 * fsh0_171[k]
                   - f_11 * fsh1_171[k]
                   + f_3 * pc_x[k] * fsi_227[k];

        t_292[k] = f_10 * fsh0_172[k]
                   - f_11 * fsh1_172[k]
                   + f_3 * pc_x[k] * fsi_228[k];

        t_293[k] = pa_y[k] * dsk0_185[k]
                   - f_12 * pc_y[k] * dsk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, fsh0_174, fsh0_175, fsh0_176, fsh1_174, \
                         fsh1_175, fsh1_176, fsi_230, fsi_231, \
                         fsi_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_8 * fsh0_174[k]
                   - f_9 * fsh1_174[k]
                   + f_3 * pc_x[k] * fsi_230[k];

        t_295[k] = f_8 * fsh0_175[k]
                   - f_9 * fsh1_175[k]
                   + f_3 * pc_x[k] * fsi_231[k];

        t_296[k] = f_8 * fsh0_176[k]
                   - f_9 * fsh1_176[k]
                   + f_3 * pc_x[k] * fsi_232[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pa_y, pc_x, pc_y, dsk0_189, dsk1_189, fsh0_178, \
                         fsh0_179, fsh1_178, fsh1_179, fsi_234, \
                         fsi_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pa_y[k] * dsk0_189[k]
                   - f_12 * pc_y[k] * dsk1_189[k];

        t_298[k] = f_6 * fsh0_178[k]
                   - f_7 * fsh1_178[k]
                   + f_3 * pc_x[k] * fsi_234[k];

        t_299[k] = f_6 * fsh0_179[k]
                   - f_7 * fsh1_179[k]
                   + f_3 * pc_x[k] * fsi_235[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pa_y, pc_x, pc_y, dsk0_194, dsk1_194, fsh0_180, \
                         fsh0_181, fsh1_180, fsh1_181, fsi_236, \
                         fsi_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_6 * fsh0_180[k]
                   - f_7 * fsh1_180[k]
                   + f_3 * pc_x[k] * fsi_236[k];

        t_301[k] = f_6 * fsh0_181[k]
                   - f_7 * fsh1_181[k]
                   + f_3 * pc_x[k] * fsi_237[k];

        t_302[k] = pa_y[k] * dsk0_194[k]
                   - f_12 * pc_y[k] * dsk1_194[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pc_x, fsh0_183, fsh0_184, fsh0_185, fsh1_183, \
                         fsh1_184, fsh1_185, fsi_239, fsi_240, \
                         fsi_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_4 * fsh0_183[k]
                   - f_5 * fsh1_183[k]
                   + f_3 * pc_x[k] * fsi_239[k];

        t_304[k] = f_4 * fsh0_184[k]
                   - f_5 * fsh1_184[k]
                   + f_3 * pc_x[k] * fsi_240[k];

        t_305[k] = f_4 * fsh0_185[k]
                   - f_5 * fsh1_185[k]
                   + f_3 * pc_x[k] * fsi_241[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pa_y, pc_x, pc_y, dsk0_200, dsk1_200, \
                         fsh0_186, fsh0_187, fsh1_186, fsh1_187, fsi_242, fsi_243, \
                         fsi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_4 * fsh0_186[k]
                   - f_5 * fsh1_186[k]
                   + f_3 * pc_x[k] * fsi_242[k];

        t_307[k] = f_4 * fsh0_187[k]
                   - f_5 * fsh1_187[k]
                   + f_3 * pc_x[k] * fsi_243[k];

        t_308[k] = pa_y[k] * dsk0_200[k]
                   - f_12 * pc_y[k] * dsk1_200[k];

        t_309[k] = f_3 * pc_x[k] * fsi_245[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, pc_x, fsi_246, fsi_247, \
                         fsi_248, fsi_249, fsi_250, fsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_3 * pc_x[k] * fsi_246[k];

        t_311[k] = f_3 * pc_x[k] * fsi_247[k];

        t_312[k] = f_3 * pc_x[k] * fsi_248[k];

        t_313[k] = f_3 * pc_x[k] * fsi_249[k];

        t_314[k] = f_3 * pc_x[k] * fsi_250[k];

        t_315[k] = f_3 * pc_x[k] * fsi_251[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_y, pc_y, pc_z, dsk0_208, dsk0_210, dsi_133, \
                         dsi_161, dsi_163, dsk1_208, dsk1_210, \
                         fsi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_y[k] * dsk0_208[k]
                   + f_19 * dsi_161[k]
                   - f_12 * pc_y[k] * dsk1_208[k];

        t_317[k] = f_14 * dsi_133[k]
                   + f_3 * pc_z[k] * fsi_245[k];

        t_318[k] = pa_y[k] * dsk0_210[k]
                   + f_16 * dsi_163[k]
                   - f_12 * pc_y[k] * dsk1_210[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_y, pc_y, dsk0_211, dsk0_212, dsk0_213, \
                         dsi_164, dsi_165, dsi_166, dsk1_211, dsk1_212, \
                         dsk1_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = pa_y[k] * dsk0_211[k]
                   + f_15 * dsi_164[k]
                   - f_12 * pc_y[k] * dsk1_211[k];

        t_320[k] = pa_y[k] * dsk0_212[k]
                   + f_0 * dsi_165[k]
                   - f_12 * pc_y[k] * dsk1_212[k];

        t_321[k] = pa_y[k] * dsk0_213[k]
                   + f_14 * dsi_166[k]
                   - f_12 * pc_y[k] * dsk1_213[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pa_y, pc_x, pc_y, dsk0_215, dsi_167, \
                         dsk1_215, fsh0_189, fsh1_189, fsi_251, \
                         fsi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_13 * dsi_167[k]
                   + f_3 * pc_y[k] * fsi_251[k];

        t_323[k] = pa_y[k] * dsk0_215[k]
                   - f_12 * pc_y[k] * dsk1_215[k];

        t_324[k] = f_1 * fsh0_189[k]
                   - f_2 * fsh1_189[k]
                   + f_3 * pc_x[k] * fsi_252[k];

        t_325[k] = f_3 * pc_y[k] * fsi_252[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_x, pc_y, fsh0_191, fsh0_192, fsh0_194, \
                         fsh1_191, fsh1_192, fsh1_194, fsi_254, fsi_255, \
                         fsi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_17 * fsh0_191[k]
                   - f_18 * fsh1_191[k]
                   + f_3 * pc_x[k] * fsi_254[k];

        t_327[k] = f_10 * fsh0_192[k]
                   - f_11 * fsh1_192[k]
                   + f_3 * pc_x[k] * fsi_255[k];

        t_328[k] = f_3 * pc_y[k] * fsi_254[k];

        t_329[k] = f_10 * fsh0_194[k]
                   - f_11 * fsh1_194[k]
                   + f_3 * pc_x[k] * fsi_257[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, pc_y, fsh0_195, fsh0_196, fsh0_198, \
                         fsh1_195, fsh1_196, fsh1_198, fsi_257, fsi_258, fsi_259, \
                         fsi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_8 * fsh0_195[k]
                   - f_9 * fsh1_195[k]
                   + f_3 * pc_x[k] * fsi_258[k];

        t_331[k] = f_8 * fsh0_196[k]
                   - f_9 * fsh1_196[k]
                   + f_3 * pc_x[k] * fsi_259[k];

        t_332[k] = f_3 * pc_y[k] * fsi_257[k];

        t_333[k] = f_8 * fsh0_198[k]
                   - f_9 * fsh1_198[k]
                   + f_3 * pc_x[k] * fsi_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, fsh0_199, fsh0_200, fsh0_201, \
                         fsh1_199, fsh1_200, fsh1_201, fsi_261, fsi_262, fsi_263, \
                         fsi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_6 * fsh0_199[k]
                   - f_7 * fsh1_199[k]
                   + f_3 * pc_x[k] * fsi_262[k];

        t_335[k] = f_6 * fsh0_200[k]
                   - f_7 * fsh1_200[k]
                   + f_3 * pc_x[k] * fsi_263[k];

        t_336[k] = f_6 * fsh0_201[k]
                   - f_7 * fsh1_201[k]
                   + f_3 * pc_x[k] * fsi_264[k];

        t_337[k] = f_3 * pc_y[k] * fsi_261[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pc_x, fsh0_203, fsh0_204, fsh0_205, fsh1_203, \
                         fsh1_204, fsh1_205, fsi_266, fsi_267, \
                         fsi_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_6 * fsh0_203[k]
                   - f_7 * fsh1_203[k]
                   + f_3 * pc_x[k] * fsi_266[k];

        t_339[k] = f_4 * fsh0_204[k]
                   - f_5 * fsh1_204[k]
                   + f_3 * pc_x[k] * fsi_267[k];

        t_340[k] = f_4 * fsh0_205[k]
                   - f_5 * fsh1_205[k]
                   + f_3 * pc_x[k] * fsi_268[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, fsh0_206, fsh0_207, fsh0_209, \
                         fsh1_206, fsh1_207, fsh1_209, fsi_266, fsi_269, fsi_270, \
                         fsi_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_4 * fsh0_206[k]
                   - f_5 * fsh1_206[k]
                   + f_3 * pc_x[k] * fsi_269[k];

        t_342[k] = f_4 * fsh0_207[k]
                   - f_5 * fsh1_207[k]
                   + f_3 * pc_x[k] * fsi_270[k];

        t_343[k] = f_3 * pc_y[k] * fsi_266[k];

        t_344[k] = f_4 * fsh0_209[k]
                   - f_5 * fsh1_209[k]
                   + f_3 * pc_x[k] * fsi_272[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, t_350, t_351, pc_x, fsi_273, \
                         fsi_274, fsi_275, fsi_276, fsi_277, fsi_278, \
                         fsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_3 * pc_x[k] * fsi_273[k];

        t_346[k] = f_3 * pc_x[k] * fsi_274[k];

        t_347[k] = f_3 * pc_x[k] * fsi_275[k];

        t_348[k] = f_3 * pc_x[k] * fsi_276[k];

        t_349[k] = f_3 * pc_x[k] * fsi_277[k];

        t_350[k] = f_3 * pc_x[k] * fsi_278[k];

        t_351[k] = f_3 * pc_x[k] * fsi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_y, fsh0_204, fsh0_205, fsh0_206, fsh1_204, \
                         fsh1_205, fsh1_206, fsi_273, fsi_274, \
                         fsi_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * fsh0_204[k]
                   - f_2 * fsh1_204[k]
                   + f_3 * pc_y[k] * fsi_273[k];

        t_353[k] = f_17 * fsh0_205[k]
                   - f_18 * fsh1_205[k]
                   + f_3 * pc_y[k] * fsi_274[k];

        t_354[k] = f_10 * fsh0_206[k]
                   - f_11 * fsh1_206[k]
                   + f_3 * pc_y[k] * fsi_275[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_y, fsh0_207, fsh0_208, fsh0_209, \
                         fsh1_207, fsh1_208, fsh1_209, fsi_276, fsi_277, fsi_278, \
                         fsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_8 * fsh0_207[k]
                   - f_9 * fsh1_207[k]
                   + f_3 * pc_y[k] * fsi_276[k];

        t_356[k] = f_6 * fsh0_208[k]
                   - f_7 * fsh1_208[k]
                   + f_3 * pc_y[k] * fsi_277[k];

        t_357[k] = f_4 * fsh0_209[k]
                   - f_5 * fsh1_209[k]
                   + f_3 * pc_y[k] * fsi_278[k];

        t_358[k] = f_3 * pc_y[k] * fsi_279[k];
    }

#pragma omp simd aligned(t_359, pc_z, dsi_167, fsh0_209, fsh1_209, \
                         fsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_0 * dsi_167[k]
                   + f_1 * fsh0_209[k]
                   - f_2 * fsh1_209[k]
                   + f_3 * pc_z[k] * fsi_279[k];
    }
}

auto
compute_prim_fsk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsk0, const size_t dsi,
                                                   const size_t dsk1, const size_t fsh0,
                                                   const size_t fsh1, const size_t fsi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fsk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, dsk0, dsi,
                                                              dsk1, fsh0, fsh1, fsi, ncols,
                                                              gamma, p, q);

    compute_prim_fsk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, dsk0, dsi,
                                                              dsk1, fsh0, fsh1, fsi, ncols,
                                                              gamma, p, q);

    compute_prim_fsk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, dsk0, dsi,
                                                              dsk1, fsh0, fsh1, fsi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
