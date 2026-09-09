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


#include "SimdThreeCenterElectronRepulsionVrrRecDSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dsk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t psk0,
                                                          const size_t psi, const size_t psk1,
                                                          const size_t dsh0, const size_t dsh1,
                                                          const size_t dsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 1.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psk0_0 = buffer.data(psk0 + 0);
    const auto *psk0_3 = buffer.data(psk0 + 3);
    const auto *psk0_5 = buffer.data(psk0 + 5);
    const auto *psk0_6 = buffer.data(psk0 + 6);
    const auto *psk0_9 = buffer.data(psk0 + 9);
    const auto *psk0_10 = buffer.data(psk0 + 10);
    const auto *psk0_14 = buffer.data(psk0 + 14);
    const auto *psk0_15 = buffer.data(psk0 + 15);
    const auto *psk0_20 = buffer.data(psk0 + 20);
    const auto *psk0_39 = buffer.data(psk0 + 39);
    const auto *psk0_42 = buffer.data(psk0 + 42);
    const auto *psk0_46 = buffer.data(psk0 + 46);
    const auto *psk0_51 = buffer.data(psk0 + 51);
    const auto *psk0_64 = buffer.data(psk0 + 64);
    const auto *psk0_66 = buffer.data(psk0 + 66);
    const auto *psk0_67 = buffer.data(psk0 + 67);
    const auto *psk0_68 = buffer.data(psk0 + 68);
    const auto *psk0_69 = buffer.data(psk0 + 69);
    const auto *psk0_71 = buffer.data(psk0 + 71);
    const auto *psk0_77 = buffer.data(psk0 + 77);
    const auto *psk0_81 = buffer.data(psk0 + 81);
    const auto *psk0_86 = buffer.data(psk0 + 86);
    const auto *psk0_92 = buffer.data(psk0 + 92);
    const auto *psk0_100 = buffer.data(psk0 + 100);
    const auto *psk0_101 = buffer.data(psk0 + 101);
    const auto *psk0_102 = buffer.data(psk0 + 102);
    const auto *psk0_103 = buffer.data(psk0 + 103);
    const auto *psk0_104 = buffer.data(psk0 + 104);
    const auto *psk0_105 = buffer.data(psk0 + 105);
    const auto *psk0_107 = buffer.data(psk0 + 107);

    const auto *psi_0 = buffer.data(psi + 0);
    const auto *psi_5 = buffer.data(psi + 5);
    const auto *psi_9 = buffer.data(psi + 9);
    const auto *psi_14 = buffer.data(psi + 14);
    const auto *psi_21 = buffer.data(psi + 21);
    const auto *psi_23 = buffer.data(psi + 23);
    const auto *psi_24 = buffer.data(psi + 24);
    const auto *psi_25 = buffer.data(psi + 25);
    const auto *psi_27 = buffer.data(psi + 27);
    const auto *psi_31 = buffer.data(psi + 31);
    const auto *psi_34 = buffer.data(psi + 34);
    const auto *psi_38 = buffer.data(psi + 38);
    const auto *psi_43 = buffer.data(psi + 43);
    const auto *psi_49 = buffer.data(psi + 49);
    const auto *psi_51 = buffer.data(psi + 51);
    const auto *psi_52 = buffer.data(psi + 52);
    const auto *psi_53 = buffer.data(psi + 53);
    const auto *psi_54 = buffer.data(psi + 54);
    const auto *psi_55 = buffer.data(psi + 55);
    const auto *psi_61 = buffer.data(psi + 61);
    const auto *psi_65 = buffer.data(psi + 65);
    const auto *psi_70 = buffer.data(psi + 70);
    const auto *psi_76 = buffer.data(psi + 76);
    const auto *psi_77 = buffer.data(psi + 77);
    const auto *psi_78 = buffer.data(psi + 78);
    const auto *psi_79 = buffer.data(psi + 79);
    const auto *psi_80 = buffer.data(psi + 80);
    const auto *psi_81 = buffer.data(psi + 81);
    const auto *psi_83 = buffer.data(psi + 83);

    const auto *psk1_0 = buffer.data(psk1 + 0);
    const auto *psk1_3 = buffer.data(psk1 + 3);
    const auto *psk1_5 = buffer.data(psk1 + 5);
    const auto *psk1_6 = buffer.data(psk1 + 6);
    const auto *psk1_9 = buffer.data(psk1 + 9);
    const auto *psk1_10 = buffer.data(psk1 + 10);
    const auto *psk1_14 = buffer.data(psk1 + 14);
    const auto *psk1_15 = buffer.data(psk1 + 15);
    const auto *psk1_20 = buffer.data(psk1 + 20);
    const auto *psk1_39 = buffer.data(psk1 + 39);
    const auto *psk1_42 = buffer.data(psk1 + 42);
    const auto *psk1_46 = buffer.data(psk1 + 46);
    const auto *psk1_51 = buffer.data(psk1 + 51);
    const auto *psk1_64 = buffer.data(psk1 + 64);
    const auto *psk1_66 = buffer.data(psk1 + 66);
    const auto *psk1_67 = buffer.data(psk1 + 67);
    const auto *psk1_68 = buffer.data(psk1 + 68);
    const auto *psk1_69 = buffer.data(psk1 + 69);
    const auto *psk1_71 = buffer.data(psk1 + 71);
    const auto *psk1_77 = buffer.data(psk1 + 77);
    const auto *psk1_81 = buffer.data(psk1 + 81);
    const auto *psk1_86 = buffer.data(psk1 + 86);
    const auto *psk1_92 = buffer.data(psk1 + 92);
    const auto *psk1_100 = buffer.data(psk1 + 100);
    const auto *psk1_101 = buffer.data(psk1 + 101);
    const auto *psk1_102 = buffer.data(psk1 + 102);
    const auto *psk1_103 = buffer.data(psk1 + 103);
    const auto *psk1_104 = buffer.data(psk1 + 104);
    const auto *psk1_105 = buffer.data(psk1 + 105);
    const auto *psk1_107 = buffer.data(psk1 + 107);

    const auto *dsh0_0 = buffer.data(dsh0 + 0);
    const auto *dsh0_1 = buffer.data(dsh0 + 1);
    const auto *dsh0_2 = buffer.data(dsh0 + 2);
    const auto *dsh0_3 = buffer.data(dsh0 + 3);
    const auto *dsh0_5 = buffer.data(dsh0 + 5);
    const auto *dsh0_6 = buffer.data(dsh0 + 6);
    const auto *dsh0_8 = buffer.data(dsh0 + 8);
    const auto *dsh0_9 = buffer.data(dsh0 + 9);
    const auto *dsh0_15 = buffer.data(dsh0 + 15);
    const auto *dsh0_17 = buffer.data(dsh0 + 17);
    const auto *dsh0_18 = buffer.data(dsh0 + 18);
    const auto *dsh0_19 = buffer.data(dsh0 + 19);
    const auto *dsh0_20 = buffer.data(dsh0 + 20);
    const auto *dsh0_24 = buffer.data(dsh0 + 24);
    const auto *dsh0_27 = buffer.data(dsh0 + 27);
    const auto *dsh0_28 = buffer.data(dsh0 + 28);
    const auto *dsh0_44 = buffer.data(dsh0 + 44);
    const auto *dsh0_46 = buffer.data(dsh0 + 46);
    const auto *dsh0_47 = buffer.data(dsh0 + 47);
    const auto *dsh0_49 = buffer.data(dsh0 + 49);
    const auto *dsh0_50 = buffer.data(dsh0 + 50);
    const auto *dsh0_51 = buffer.data(dsh0 + 51);
    const auto *dsh0_63 = buffer.data(dsh0 + 63);
    const auto *dsh0_64 = buffer.data(dsh0 + 64);
    const auto *dsh0_66 = buffer.data(dsh0 + 66);
    const auto *dsh0_68 = buffer.data(dsh0 + 68);
    const auto *dsh0_69 = buffer.data(dsh0 + 69);
    const auto *dsh0_71 = buffer.data(dsh0 + 71);
    const auto *dsh0_72 = buffer.data(dsh0 + 72);
    const auto *dsh0_73 = buffer.data(dsh0 + 73);
    const auto *dsh0_75 = buffer.data(dsh0 + 75);
    const auto *dsh0_76 = buffer.data(dsh0 + 76);
    const auto *dsh0_77 = buffer.data(dsh0 + 77);
    const auto *dsh0_78 = buffer.data(dsh0 + 78);
    const auto *dsh0_80 = buffer.data(dsh0 + 80);
    const auto *dsh0_81 = buffer.data(dsh0 + 81);
    const auto *dsh0_82 = buffer.data(dsh0 + 82);

    const auto *dsh1_0 = buffer.data(dsh1 + 0);
    const auto *dsh1_1 = buffer.data(dsh1 + 1);
    const auto *dsh1_2 = buffer.data(dsh1 + 2);
    const auto *dsh1_3 = buffer.data(dsh1 + 3);
    const auto *dsh1_5 = buffer.data(dsh1 + 5);
    const auto *dsh1_6 = buffer.data(dsh1 + 6);
    const auto *dsh1_8 = buffer.data(dsh1 + 8);
    const auto *dsh1_9 = buffer.data(dsh1 + 9);
    const auto *dsh1_15 = buffer.data(dsh1 + 15);
    const auto *dsh1_17 = buffer.data(dsh1 + 17);
    const auto *dsh1_18 = buffer.data(dsh1 + 18);
    const auto *dsh1_19 = buffer.data(dsh1 + 19);
    const auto *dsh1_20 = buffer.data(dsh1 + 20);
    const auto *dsh1_24 = buffer.data(dsh1 + 24);
    const auto *dsh1_27 = buffer.data(dsh1 + 27);
    const auto *dsh1_28 = buffer.data(dsh1 + 28);
    const auto *dsh1_44 = buffer.data(dsh1 + 44);
    const auto *dsh1_46 = buffer.data(dsh1 + 46);
    const auto *dsh1_47 = buffer.data(dsh1 + 47);
    const auto *dsh1_49 = buffer.data(dsh1 + 49);
    const auto *dsh1_50 = buffer.data(dsh1 + 50);
    const auto *dsh1_51 = buffer.data(dsh1 + 51);
    const auto *dsh1_63 = buffer.data(dsh1 + 63);
    const auto *dsh1_64 = buffer.data(dsh1 + 64);
    const auto *dsh1_66 = buffer.data(dsh1 + 66);
    const auto *dsh1_68 = buffer.data(dsh1 + 68);
    const auto *dsh1_69 = buffer.data(dsh1 + 69);
    const auto *dsh1_71 = buffer.data(dsh1 + 71);
    const auto *dsh1_72 = buffer.data(dsh1 + 72);
    const auto *dsh1_73 = buffer.data(dsh1 + 73);
    const auto *dsh1_75 = buffer.data(dsh1 + 75);
    const auto *dsh1_76 = buffer.data(dsh1 + 76);
    const auto *dsh1_77 = buffer.data(dsh1 + 77);
    const auto *dsh1_78 = buffer.data(dsh1 + 78);
    const auto *dsh1_80 = buffer.data(dsh1 + 80);
    const auto *dsh1_81 = buffer.data(dsh1 + 81);
    const auto *dsh1_82 = buffer.data(dsh1 + 82);

    const auto *dsi_0 = buffer.data(dsi + 0);
    const auto *dsi_1 = buffer.data(dsi + 1);
    const auto *dsi_2 = buffer.data(dsi + 2);
    const auto *dsi_3 = buffer.data(dsi + 3);
    const auto *dsi_5 = buffer.data(dsi + 5);
    const auto *dsi_6 = buffer.data(dsi + 6);
    const auto *dsi_8 = buffer.data(dsi + 8);
    const auto *dsi_9 = buffer.data(dsi + 9);
    const auto *dsi_10 = buffer.data(dsi + 10);
    const auto *dsi_12 = buffer.data(dsi + 12);
    const auto *dsi_13 = buffer.data(dsi + 13);
    const auto *dsi_14 = buffer.data(dsi + 14);
    const auto *dsi_15 = buffer.data(dsi + 15);
    const auto *dsi_20 = buffer.data(dsi + 20);
    const auto *dsi_21 = buffer.data(dsi + 21);
    const auto *dsi_23 = buffer.data(dsi + 23);
    const auto *dsi_24 = buffer.data(dsi + 24);
    const auto *dsi_25 = buffer.data(dsi + 25);
    const auto *dsi_26 = buffer.data(dsi + 26);
    const auto *dsi_27 = buffer.data(dsi + 27);
    const auto *dsi_28 = buffer.data(dsi + 28);
    const auto *dsi_29 = buffer.data(dsi + 29);
    const auto *dsi_31 = buffer.data(dsi + 31);
    const auto *dsi_33 = buffer.data(dsi + 33);
    const auto *dsi_34 = buffer.data(dsi + 34);
    const auto *dsi_35 = buffer.data(dsi + 35);
    const auto *dsi_37 = buffer.data(dsi + 37);
    const auto *dsi_38 = buffer.data(dsi + 38);
    const auto *dsi_39 = buffer.data(dsi + 39);
    const auto *dsi_40 = buffer.data(dsi + 40);
    const auto *dsi_42 = buffer.data(dsi + 42);
    const auto *dsi_43 = buffer.data(dsi + 43);
    const auto *dsi_49 = buffer.data(dsi + 49);
    const auto *dsi_51 = buffer.data(dsi + 51);
    const auto *dsi_52 = buffer.data(dsi + 52);
    const auto *dsi_53 = buffer.data(dsi + 53);
    const auto *dsi_54 = buffer.data(dsi + 54);
    const auto *dsi_55 = buffer.data(dsi + 55);
    const auto *dsi_56 = buffer.data(dsi + 56);
    const auto *dsi_58 = buffer.data(dsi + 58);
    const auto *dsi_60 = buffer.data(dsi + 60);
    const auto *dsi_61 = buffer.data(dsi + 61);
    const auto *dsi_63 = buffer.data(dsi + 63);
    const auto *dsi_64 = buffer.data(dsi + 64);
    const auto *dsi_65 = buffer.data(dsi + 65);
    const auto *dsi_67 = buffer.data(dsi + 67);
    const auto *dsi_68 = buffer.data(dsi + 68);
    const auto *dsi_69 = buffer.data(dsi + 69);
    const auto *dsi_70 = buffer.data(dsi + 70);
    const auto *dsi_76 = buffer.data(dsi + 76);
    const auto *dsi_77 = buffer.data(dsi + 77);
    const auto *dsi_78 = buffer.data(dsi + 78);
    const auto *dsi_79 = buffer.data(dsi + 79);
    const auto *dsi_80 = buffer.data(dsi + 80);
    const auto *dsi_81 = buffer.data(dsi + 81);
    const auto *dsi_83 = buffer.data(dsi + 83);
    const auto *dsi_84 = buffer.data(dsi + 84);
    const auto *dsi_85 = buffer.data(dsi + 85);
    const auto *dsi_87 = buffer.data(dsi + 87);
    const auto *dsi_89 = buffer.data(dsi + 89);
    const auto *dsi_90 = buffer.data(dsi + 90);
    const auto *dsi_92 = buffer.data(dsi + 92);
    const auto *dsi_93 = buffer.data(dsi + 93);
    const auto *dsi_94 = buffer.data(dsi + 94);
    const auto *dsi_96 = buffer.data(dsi + 96);
    const auto *dsi_97 = buffer.data(dsi + 97);
    const auto *dsi_98 = buffer.data(dsi + 98);
    const auto *dsi_99 = buffer.data(dsi + 99);
    const auto *dsi_101 = buffer.data(dsi + 101);
    const auto *dsi_102 = buffer.data(dsi + 102);
    const auto *dsi_103 = buffer.data(dsi + 103);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, psi_0, dsh0_0, \
                         dsh1_0, dsi_0, dsi_1, dsi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psi_0[k]
                 + f_1 * dsh0_0[k]
                 - f_2 * dsh1_0[k]
                 + f_3 * pc_x[k] * dsi_0[k];

        t_1[k] = f_3 * pc_y[k] * dsi_0[k];

        t_2[k] = f_3 * pc_z[k] * dsi_0[k];

        t_3[k] = f_4 * dsh0_0[k]
                 - f_5 * dsh1_0[k]
                 + f_3 * pc_y[k] * dsi_1[k];

        t_4[k] = f_3 * pc_y[k] * dsi_2[k];

        t_5[k] = f_4 * dsh0_0[k]
                 - f_5 * dsh1_0[k]
                 + f_3 * pc_z[k] * dsi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, dsh0_1, dsh0_2, dsh0_3, dsh1_1, \
                         dsh1_2, dsh1_3, dsi_3, dsi_5, dsi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * dsh0_1[k]
                 - f_7 * dsh1_1[k]
                 + f_3 * pc_y[k] * dsi_3[k];

        t_7[k] = f_3 * pc_z[k] * dsi_3[k];

        t_8[k] = f_3 * pc_y[k] * dsi_5[k];

        t_9[k] = f_6 * dsh0_2[k]
                 - f_7 * dsh1_2[k]
                 + f_3 * pc_z[k] * dsi_5[k];

        t_10[k] = f_8 * dsh0_3[k]
                  - f_9 * dsh1_3[k]
                  + f_3 * pc_y[k] * dsi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, dsh0_5, dsh0_6, \
                         dsh1_5, dsh1_6, dsi_6, dsi_8, dsi_9, dsi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * dsi_6[k];

        t_12[k] = f_4 * dsh0_5[k]
                  - f_5 * dsh1_5[k]
                  + f_3 * pc_y[k] * dsi_8[k];

        t_13[k] = f_3 * pc_y[k] * dsi_9[k];

        t_14[k] = f_8 * dsh0_5[k]
                  - f_9 * dsh1_5[k]
                  + f_3 * pc_z[k] * dsi_9[k];

        t_15[k] = f_10 * dsh0_6[k]
                  - f_11 * dsh1_6[k]
                  + f_3 * pc_y[k] * dsi_10[k];

        t_16[k] = f_3 * pc_z[k] * dsi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, dsh0_8, dsh0_9, dsh1_8, dsh1_9, \
                         dsi_12, dsi_13, dsi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * dsh0_8[k]
                  - f_7 * dsh1_8[k]
                  + f_3 * pc_y[k] * dsi_12[k];

        t_18[k] = f_4 * dsh0_9[k]
                  - f_5 * dsh1_9[k]
                  + f_3 * pc_y[k] * dsi_13[k];

        t_19[k] = f_3 * pc_y[k] * dsi_14[k];

        t_20[k] = f_10 * dsh0_9[k]
                  - f_11 * dsh1_9[k]
                  + f_3 * pc_z[k] * dsi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, psi_21, psi_23, psi_24, \
                         psi_25, dsi_15, dsi_21, dsi_23, dsi_24, \
                         dsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * psi_21[k]
                  + f_3 * pc_x[k] * dsi_21[k];

        t_22[k] = f_3 * pc_z[k] * dsi_15[k];

        t_23[k] = f_0 * psi_23[k]
                  + f_3 * pc_x[k] * dsi_23[k];

        t_24[k] = f_0 * psi_24[k]
                  + f_3 * pc_x[k] * dsi_24[k];

        t_25[k] = f_0 * psi_25[k]
                  + f_3 * pc_x[k] * dsi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, psi_27, dsh0_15, dsh1_15, \
                         dsi_20, dsi_21, dsi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * dsi_20[k];

        t_27[k] = f_0 * psi_27[k]
                  + f_3 * pc_x[k] * dsi_27[k];

        t_28[k] = f_1 * dsh0_15[k]
                  - f_2 * dsh1_15[k]
                  + f_3 * pc_y[k] * dsi_21[k];

        t_29[k] = f_3 * pc_z[k] * dsi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, dsh0_17, dsh0_18, dsh0_19, dsh1_17, dsh1_18, \
                         dsh1_19, dsi_23, dsi_24, dsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * dsh0_17[k]
                  - f_11 * dsh1_17[k]
                  + f_3 * pc_y[k] * dsi_23[k];

        t_31[k] = f_8 * dsh0_18[k]
                  - f_9 * dsh1_18[k]
                  + f_3 * pc_y[k] * dsi_24[k];

        t_32[k] = f_6 * dsh0_19[k]
                  - f_7 * dsh1_19[k]
                  + f_3 * pc_y[k] * dsi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, psk0_0, psi_0, \
                         psk1_0, dsh0_20, dsh1_20, dsi_26, dsi_27, \
                         dsi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * dsh0_20[k]
                  - f_5 * dsh1_20[k]
                  + f_3 * pc_y[k] * dsi_26[k];

        t_34[k] = f_3 * pc_y[k] * dsi_27[k];

        t_35[k] = f_1 * dsh0_20[k]
                  - f_2 * dsh1_20[k]
                  + f_3 * pc_z[k] * dsi_27[k];

        t_36[k] = pa_y[k] * psk0_0[k]
                  - f_12 * pc_y[k] * psk1_0[k];

        t_37[k] = f_13 * psi_0[k]
                  + f_3 * pc_y[k] * dsi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pa_y, pc_x, pc_y, pc_z, psk0_5, \
                         psk0_39, psi_31, psk1_5, psk1_39, dsi_28, \
                         dsi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * dsi_28[k];

        t_39[k] = pa_x[k] * psk0_39[k]
                  + f_14 * psi_31[k]
                  - f_12 * pc_x[k] * psk1_39[k];

        t_40[k] = f_3 * pc_z[k] * dsi_29[k];

        t_41[k] = pa_y[k] * psk0_5[k]
                  - f_12 * pc_y[k] * psk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pc_x, pc_y, pc_z, psk0_42, psi_5, psi_34, \
                         psk1_42, dsi_31, dsi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * psk0_42[k]
                  + f_15 * psi_34[k]
                  - f_12 * pc_x[k] * psk1_42[k];

        t_43[k] = f_3 * pc_z[k] * dsi_31[k];

        t_44[k] = f_13 * psi_5[k]
                  + f_3 * pc_y[k] * dsi_33[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_x, pa_y, pc_x, pc_y, pc_z, psk0_9, psk0_46, \
                         psi_38, psk1_9, psk1_46, dsi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * psk0_9[k]
                  - f_12 * pc_y[k] * psk1_9[k];

        t_46[k] = pa_x[k] * psk0_46[k]
                  + f_16 * psi_38[k]
                  - f_12 * pc_x[k] * psk1_46[k];

        t_47[k] = f_3 * pc_z[k] * dsi_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_y, pc_y, pc_z, psk0_14, psi_9, psk1_14, dsh0_24, \
                         dsh1_24, dsi_35, dsi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_4 * dsh0_24[k]
                  - f_5 * dsh1_24[k]
                  + f_3 * pc_z[k] * dsi_35[k];

        t_49[k] = f_13 * psi_9[k]
                  + f_3 * pc_y[k] * dsi_37[k];

        t_50[k] = pa_y[k] * psk0_14[k]
                  - f_12 * pc_y[k] * psk1_14[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_x, pc_x, pc_z, psk0_51, psi_43, psk1_51, \
                         dsh0_27, dsh1_27, dsi_38, dsi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_x[k] * psk0_51[k]
                  + f_0 * psi_43[k]
                  - f_12 * pc_x[k] * psk1_51[k];

        t_52[k] = f_3 * pc_z[k] * dsi_38[k];

        t_53[k] = f_4 * dsh0_27[k]
                  - f_5 * dsh1_27[k]
                  + f_3 * pc_z[k] * dsi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, psk0_20, psi_14, psk1_20, \
                         dsh0_28, dsh1_28, dsi_40, dsi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * dsh0_28[k]
                  - f_7 * dsh1_28[k]
                  + f_3 * pc_z[k] * dsi_40[k];

        t_55[k] = f_13 * psi_14[k]
                  + f_3 * pc_y[k] * dsi_42[k];

        t_56[k] = pa_y[k] * psk0_20[k]
                  - f_12 * pc_y[k] * psk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, psi_49, psi_51, psi_52, \
                         psi_53, dsi_43, dsi_49, dsi_51, dsi_52, \
                         dsi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_13 * psi_49[k]
                  + f_3 * pc_x[k] * dsi_49[k];

        t_58[k] = f_3 * pc_z[k] * dsi_43[k];

        t_59[k] = f_13 * psi_51[k]
                  + f_3 * pc_x[k] * dsi_51[k];

        t_60[k] = f_13 * psi_52[k]
                  + f_3 * pc_x[k] * dsi_52[k];

        t_61[k] = f_13 * psi_53[k]
                  + f_3 * pc_x[k] * dsi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pc_x, pc_z, psk0_64, psi_54, psi_55, \
                         psk1_64, dsi_49, dsi_54, dsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_13 * psi_54[k]
                  + f_3 * pc_x[k] * dsi_54[k];

        t_63[k] = f_13 * psi_55[k]
                  + f_3 * pc_x[k] * dsi_55[k];

        t_64[k] = pa_x[k] * psk0_64[k]
                  - f_12 * pc_x[k] * psk1_64[k];

        t_65[k] = f_3 * pc_z[k] * dsi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pc_x, psk0_66, psk0_67, psk0_68, \
                         psk0_69, psk1_66, psk1_67, psk1_68, psk1_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * psk0_66[k]
                  - f_12 * pc_x[k] * psk1_66[k];

        t_67[k] = pa_x[k] * psk0_67[k]
                  - f_12 * pc_x[k] * psk1_67[k];

        t_68[k] = pa_x[k] * psk0_68[k]
                  - f_12 * pc_x[k] * psk1_68[k];

        t_69[k] = pa_x[k] * psk0_69[k]
                  - f_12 * pc_x[k] * psk1_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pa_z, pc_x, pc_y, pc_z, psk0_0, \
                         psk0_71, psi_27, psk1_0, psk1_71, dsi_55, \
                         dsi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * psi_27[k]
                  + f_3 * pc_y[k] * dsi_55[k];

        t_71[k] = pa_x[k] * psk0_71[k]
                  - f_12 * pc_x[k] * psk1_71[k];

        t_72[k] = pa_z[k] * psk0_0[k]
                  - f_12 * pc_z[k] * psk1_0[k];

        t_73[k] = f_3 * pc_y[k] * dsi_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_z, pc_y, pc_z, psk0_3, psi_0, psk1_3, dsi_56, \
                         dsi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * psi_0[k]
                  + f_3 * pc_z[k] * dsi_56[k];

        t_75[k] = pa_z[k] * psk0_3[k]
                  - f_12 * pc_z[k] * psk1_3[k];

        t_76[k] = f_3 * pc_y[k] * dsi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_x, pa_z, pc_x, pc_y, pc_z, psk0_6, psk0_77, \
                         psi_61, psk1_6, psk1_77, dsh0_44, dsh1_44, \
                         dsi_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_x[k] * psk0_77[k]
                  + f_14 * psi_61[k]
                  - f_12 * pc_x[k] * psk1_77[k];

        t_78[k] = pa_z[k] * psk0_6[k]
                  - f_12 * pc_z[k] * psk1_6[k];

        t_79[k] = f_4 * dsh0_44[k]
                  - f_5 * dsh1_44[k]
                  + f_3 * pc_y[k] * dsi_60[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_x, pa_z, pc_x, pc_y, pc_z, psk0_10, psk0_81, \
                         psi_65, psk1_10, psk1_81, dsi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_3 * pc_y[k] * dsi_61[k];

        t_81[k] = pa_x[k] * psk0_81[k]
                  + f_15 * psi_65[k]
                  - f_12 * pc_x[k] * psk1_81[k];

        t_82[k] = pa_z[k] * psk0_10[k]
                  - f_12 * pc_z[k] * psk1_10[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pc_y, dsh0_46, dsh0_47, dsh1_46, dsh1_47, dsi_63, \
                         dsi_64, dsi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_6 * dsh0_46[k]
                  - f_7 * dsh1_46[k]
                  + f_3 * pc_y[k] * dsi_63[k];

        t_84[k] = f_4 * dsh0_47[k]
                  - f_5 * dsh1_47[k]
                  + f_3 * pc_y[k] * dsi_64[k];

        t_85[k] = f_3 * pc_y[k] * dsi_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_x, pa_z, pc_x, pc_y, pc_z, psk0_15, psk0_86, \
                         psi_70, psk1_15, psk1_86, dsh0_49, dsh1_49, \
                         dsi_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * psk0_86[k]
                  + f_16 * psi_70[k]
                  - f_12 * pc_x[k] * psk1_86[k];

        t_87[k] = pa_z[k] * psk0_15[k]
                  - f_12 * pc_z[k] * psk1_15[k];

        t_88[k] = f_8 * dsh0_49[k]
                  - f_9 * dsh1_49[k]
                  + f_3 * pc_y[k] * dsi_67[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_y, dsh0_50, dsh0_51, dsh1_50, dsh1_51, dsi_68, \
                         dsi_69, dsi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_6 * dsh0_50[k]
                  - f_7 * dsh1_50[k]
                  + f_3 * pc_y[k] * dsi_68[k];

        t_90[k] = f_4 * dsh0_51[k]
                  - f_5 * dsh1_51[k]
                  + f_3 * pc_y[k] * dsi_69[k];

        t_91[k] = f_3 * pc_y[k] * dsi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pc_x, psk0_92, psi_76, psi_77, psi_78, \
                         psi_79, psk1_92, dsi_77, dsi_78, dsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_x[k] * psk0_92[k]
                  + f_0 * psi_76[k]
                  - f_12 * pc_x[k] * psk1_92[k];

        t_93[k] = f_13 * psi_77[k]
                  + f_3 * pc_x[k] * dsi_77[k];

        t_94[k] = f_13 * psi_78[k]
                  + f_3 * pc_x[k] * dsi_78[k];

        t_95[k] = f_13 * psi_79[k]
                  + f_3 * pc_x[k] * dsi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, psi_80, psi_81, psi_83, dsi_76, \
                         dsi_80, dsi_81, dsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_13 * psi_80[k]
                  + f_3 * pc_x[k] * dsi_80[k];

        t_97[k] = f_13 * psi_81[k]
                  + f_3 * pc_x[k] * dsi_81[k];

        t_98[k] = f_3 * pc_y[k] * dsi_76[k];

        t_99[k] = f_13 * psi_83[k]
                  + f_3 * pc_x[k] * dsi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pc_x, psk0_100, psk0_101, psk0_102, \
                         psk0_103, psk1_100, psk1_101, psk1_102, \
                         psk1_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * psk0_100[k]
                   - f_12 * pc_x[k] * psk1_100[k];

        t_101[k] = pa_x[k] * psk0_101[k]
                   - f_12 * pc_x[k] * psk1_101[k];

        t_102[k] = pa_x[k] * psk0_102[k]
                   - f_12 * pc_x[k] * psk1_102[k];

        t_103[k] = pa_x[k] * psk0_103[k]
                   - f_12 * pc_x[k] * psk1_103[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pc_x, pc_y, psk0_104, psk0_105, \
                         psk0_107, psk1_104, psk1_105, psk1_107, \
                         dsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_x[k] * psk0_104[k]
                   - f_12 * pc_x[k] * psk1_104[k];

        t_105[k] = pa_x[k] * psk0_105[k]
                   - f_12 * pc_x[k] * psk1_105[k];

        t_106[k] = f_3 * pc_y[k] * dsi_83[k];

        t_107[k] = pa_x[k] * psk0_107[k]
                   - f_12 * pc_x[k] * psk1_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pc_x, pc_z, dsh0_63, dsh0_64, \
                         dsh0_66, dsh1_63, dsh1_64, dsh1_66, dsi_84, dsi_85, \
                         dsi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * dsh0_63[k]
                   - f_2 * dsh1_63[k]
                   + f_3 * pc_x[k] * dsi_84[k];

        t_109[k] = f_17 * dsh0_64[k]
                   - f_18 * dsh1_64[k]
                   + f_3 * pc_x[k] * dsi_85[k];

        t_110[k] = f_3 * pc_z[k] * dsi_84[k];

        t_111[k] = f_10 * dsh0_66[k]
                   - f_11 * dsh1_66[k]
                   + f_3 * pc_x[k] * dsi_87[k];

        t_112[k] = f_3 * pc_z[k] * dsi_85[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_x, pc_z, dsh0_68, dsh0_69, dsh0_71, \
                         dsh1_68, dsh1_69, dsh1_71, dsi_87, dsi_89, dsi_90, \
                         dsi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * dsh0_68[k]
                   - f_11 * dsh1_68[k]
                   + f_3 * pc_x[k] * dsi_89[k];

        t_114[k] = f_8 * dsh0_69[k]
                   - f_9 * dsh1_69[k]
                   + f_3 * pc_x[k] * dsi_90[k];

        t_115[k] = f_3 * pc_z[k] * dsi_87[k];

        t_116[k] = f_8 * dsh0_71[k]
                   - f_9 * dsh1_71[k]
                   + f_3 * pc_x[k] * dsi_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pc_x, pc_z, dsh0_72, dsh0_73, dsh0_75, \
                         dsh1_72, dsh1_73, dsh1_75, dsi_90, dsi_93, dsi_94, \
                         dsi_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * dsh0_72[k]
                   - f_9 * dsh1_72[k]
                   + f_3 * pc_x[k] * dsi_93[k];

        t_118[k] = f_6 * dsh0_73[k]
                   - f_7 * dsh1_73[k]
                   + f_3 * pc_x[k] * dsi_94[k];

        t_119[k] = f_3 * pc_z[k] * dsi_90[k];

        t_120[k] = f_6 * dsh0_75[k]
                   - f_7 * dsh1_75[k]
                   + f_3 * pc_x[k] * dsi_96[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_z, dsh0_76, dsh0_77, dsh0_78, \
                         dsh1_76, dsh1_77, dsh1_78, dsi_94, dsi_97, dsi_98, \
                         dsi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_6 * dsh0_76[k]
                   - f_7 * dsh1_76[k]
                   + f_3 * pc_x[k] * dsi_97[k];

        t_122[k] = f_6 * dsh0_77[k]
                   - f_7 * dsh1_77[k]
                   + f_3 * pc_x[k] * dsi_98[k];

        t_123[k] = f_4 * dsh0_78[k]
                   - f_5 * dsh1_78[k]
                   + f_3 * pc_x[k] * dsi_99[k];

        t_124[k] = f_3 * pc_z[k] * dsi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, dsh0_80, dsh0_81, dsh0_82, dsh1_80, \
                         dsh1_81, dsh1_82, dsi_101, dsi_102, dsi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * dsh0_80[k]
                   - f_5 * dsh1_80[k]
                   + f_3 * pc_x[k] * dsi_101[k];

        t_126[k] = f_4 * dsh0_81[k]
                   - f_5 * dsh1_81[k]
                   + f_3 * pc_x[k] * dsi_102[k];

        t_127[k] = f_4 * dsh0_82[k]
                   - f_5 * dsh1_82[k]
                   + f_3 * pc_x[k] * dsi_103[k];
    }
}

static auto
compute_prim_dsk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t psk0,
                                                          const size_t psi, const size_t psk1,
                                                          const size_t dsh0, const size_t dsh1,
                                                          const size_t dsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 1.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);

    auto *t_128 = buffer.data(target + 128);
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psk0_37 = buffer.data(psk0 + 37);
    const auto *psk0_39 = buffer.data(psk0 + 39);
    const auto *psk0_42 = buffer.data(psk0 + 42);
    const auto *psk0_46 = buffer.data(psk0 + 46);
    const auto *psk0_51 = buffer.data(psk0 + 51);
    const auto *psk0_64 = buffer.data(psk0 + 64);
    const auto *psk0_72 = buffer.data(psk0 + 72);
    const auto *psk0_74 = buffer.data(psk0 + 74);
    const auto *psk0_77 = buffer.data(psk0 + 77);
    const auto *psk0_81 = buffer.data(psk0 + 81);
    const auto *psk0_86 = buffer.data(psk0 + 86);
    const auto *psk0_92 = buffer.data(psk0 + 92);
    const auto *psk0_102 = buffer.data(psk0 + 102);
    const auto *psk0_103 = buffer.data(psk0 + 103);
    const auto *psk0_104 = buffer.data(psk0 + 104);
    const auto *psk0_105 = buffer.data(psk0 + 105);
    const auto *psk0_107 = buffer.data(psk0 + 107);

    const auto *psi_49 = buffer.data(psi + 49);
    const auto *psi_55 = buffer.data(psi + 55);
    const auto *psi_79 = buffer.data(psi + 79);
    const auto *psi_80 = buffer.data(psi + 80);
    const auto *psi_81 = buffer.data(psi + 81);
    const auto *psi_82 = buffer.data(psi + 82);
    const auto *psi_83 = buffer.data(psi + 83);

    const auto *psk1_37 = buffer.data(psk1 + 37);
    const auto *psk1_39 = buffer.data(psk1 + 39);
    const auto *psk1_42 = buffer.data(psk1 + 42);
    const auto *psk1_46 = buffer.data(psk1 + 46);
    const auto *psk1_51 = buffer.data(psk1 + 51);
    const auto *psk1_64 = buffer.data(psk1 + 64);
    const auto *psk1_72 = buffer.data(psk1 + 72);
    const auto *psk1_74 = buffer.data(psk1 + 74);
    const auto *psk1_77 = buffer.data(psk1 + 77);
    const auto *psk1_81 = buffer.data(psk1 + 81);
    const auto *psk1_86 = buffer.data(psk1 + 86);
    const auto *psk1_92 = buffer.data(psk1 + 92);
    const auto *psk1_102 = buffer.data(psk1 + 102);
    const auto *psk1_103 = buffer.data(psk1 + 103);
    const auto *psk1_104 = buffer.data(psk1 + 104);
    const auto *psk1_105 = buffer.data(psk1 + 105);
    const auto *psk1_107 = buffer.data(psk1 + 107);

    const auto *dsh0_78 = buffer.data(dsh0 + 78);
    const auto *dsh0_79 = buffer.data(dsh0 + 79);
    const auto *dsh0_80 = buffer.data(dsh0 + 80);
    const auto *dsh0_81 = buffer.data(dsh0 + 81);
    const auto *dsh0_83 = buffer.data(dsh0 + 83);
    const auto *dsh0_88 = buffer.data(dsh0 + 88);
    const auto *dsh0_91 = buffer.data(dsh0 + 91);
    const auto *dsh0_92 = buffer.data(dsh0 + 92);
    const auto *dsh0_95 = buffer.data(dsh0 + 95);
    const auto *dsh0_96 = buffer.data(dsh0 + 96);
    const auto *dsh0_97 = buffer.data(dsh0 + 97);
    const auto *dsh0_100 = buffer.data(dsh0 + 100);
    const auto *dsh0_101 = buffer.data(dsh0 + 101);
    const auto *dsh0_102 = buffer.data(dsh0 + 102);
    const auto *dsh0_103 = buffer.data(dsh0 + 103);
    const auto *dsh0_105 = buffer.data(dsh0 + 105);
    const auto *dsh0_107 = buffer.data(dsh0 + 107);
    const auto *dsh0_108 = buffer.data(dsh0 + 108);
    const auto *dsh0_110 = buffer.data(dsh0 + 110);
    const auto *dsh0_111 = buffer.data(dsh0 + 111);
    const auto *dsh0_112 = buffer.data(dsh0 + 112);
    const auto *dsh0_114 = buffer.data(dsh0 + 114);
    const auto *dsh0_115 = buffer.data(dsh0 + 115);
    const auto *dsh0_116 = buffer.data(dsh0 + 116);
    const auto *dsh0_117 = buffer.data(dsh0 + 117);
    const auto *dsh0_119 = buffer.data(dsh0 + 119);
    const auto *dsh0_120 = buffer.data(dsh0 + 120);
    const auto *dsh0_121 = buffer.data(dsh0 + 121);
    const auto *dsh0_122 = buffer.data(dsh0 + 122);
    const auto *dsh0_123 = buffer.data(dsh0 + 123);
    const auto *dsh0_124 = buffer.data(dsh0 + 124);
    const auto *dsh0_125 = buffer.data(dsh0 + 125);

    const auto *dsh1_78 = buffer.data(dsh1 + 78);
    const auto *dsh1_79 = buffer.data(dsh1 + 79);
    const auto *dsh1_80 = buffer.data(dsh1 + 80);
    const auto *dsh1_81 = buffer.data(dsh1 + 81);
    const auto *dsh1_83 = buffer.data(dsh1 + 83);
    const auto *dsh1_88 = buffer.data(dsh1 + 88);
    const auto *dsh1_91 = buffer.data(dsh1 + 91);
    const auto *dsh1_92 = buffer.data(dsh1 + 92);
    const auto *dsh1_95 = buffer.data(dsh1 + 95);
    const auto *dsh1_96 = buffer.data(dsh1 + 96);
    const auto *dsh1_97 = buffer.data(dsh1 + 97);
    const auto *dsh1_100 = buffer.data(dsh1 + 100);
    const auto *dsh1_101 = buffer.data(dsh1 + 101);
    const auto *dsh1_102 = buffer.data(dsh1 + 102);
    const auto *dsh1_103 = buffer.data(dsh1 + 103);
    const auto *dsh1_105 = buffer.data(dsh1 + 105);
    const auto *dsh1_107 = buffer.data(dsh1 + 107);
    const auto *dsh1_108 = buffer.data(dsh1 + 108);
    const auto *dsh1_110 = buffer.data(dsh1 + 110);
    const auto *dsh1_111 = buffer.data(dsh1 + 111);
    const auto *dsh1_112 = buffer.data(dsh1 + 112);
    const auto *dsh1_114 = buffer.data(dsh1 + 114);
    const auto *dsh1_115 = buffer.data(dsh1 + 115);
    const auto *dsh1_116 = buffer.data(dsh1 + 116);
    const auto *dsh1_117 = buffer.data(dsh1 + 117);
    const auto *dsh1_119 = buffer.data(dsh1 + 119);
    const auto *dsh1_120 = buffer.data(dsh1 + 120);
    const auto *dsh1_121 = buffer.data(dsh1 + 121);
    const auto *dsh1_122 = buffer.data(dsh1 + 122);
    const auto *dsh1_123 = buffer.data(dsh1 + 123);
    const auto *dsh1_124 = buffer.data(dsh1 + 124);
    const auto *dsh1_125 = buffer.data(dsh1 + 125);

    const auto *dsi_104 = buffer.data(dsi + 104);
    const auto *dsi_105 = buffer.data(dsi + 105);
    const auto *dsi_106 = buffer.data(dsi + 106);
    const auto *dsi_107 = buffer.data(dsi + 107);
    const auto *dsi_108 = buffer.data(dsi + 108);
    const auto *dsi_109 = buffer.data(dsi + 109);
    const auto *dsi_110 = buffer.data(dsi + 110);
    const auto *dsi_111 = buffer.data(dsi + 111);
    const auto *dsi_116 = buffer.data(dsi + 116);
    const auto *dsi_119 = buffer.data(dsi + 119);
    const auto *dsi_120 = buffer.data(dsi + 120);
    const auto *dsi_123 = buffer.data(dsi + 123);
    const auto *dsi_124 = buffer.data(dsi + 124);
    const auto *dsi_125 = buffer.data(dsi + 125);
    const auto *dsi_128 = buffer.data(dsi + 128);
    const auto *dsi_129 = buffer.data(dsi + 129);
    const auto *dsi_130 = buffer.data(dsi + 130);
    const auto *dsi_131 = buffer.data(dsi + 131);
    const auto *dsi_133 = buffer.data(dsi + 133);
    const auto *dsi_134 = buffer.data(dsi + 134);
    const auto *dsi_135 = buffer.data(dsi + 135);
    const auto *dsi_136 = buffer.data(dsi + 136);
    const auto *dsi_137 = buffer.data(dsi + 137);
    const auto *dsi_138 = buffer.data(dsi + 138);
    const auto *dsi_139 = buffer.data(dsi + 139);
    const auto *dsi_140 = buffer.data(dsi + 140);
    const auto *dsi_142 = buffer.data(dsi + 142);
    const auto *dsi_143 = buffer.data(dsi + 143);
    const auto *dsi_145 = buffer.data(dsi + 145);
    const auto *dsi_146 = buffer.data(dsi + 146);
    const auto *dsi_147 = buffer.data(dsi + 147);
    const auto *dsi_149 = buffer.data(dsi + 149);
    const auto *dsi_150 = buffer.data(dsi + 150);
    const auto *dsi_151 = buffer.data(dsi + 151);
    const auto *dsi_152 = buffer.data(dsi + 152);
    const auto *dsi_154 = buffer.data(dsi + 154);
    const auto *dsi_155 = buffer.data(dsi + 155);
    const auto *dsi_156 = buffer.data(dsi + 156);
    const auto *dsi_157 = buffer.data(dsi + 157);
    const auto *dsi_158 = buffer.data(dsi + 158);
    const auto *dsi_160 = buffer.data(dsi + 160);
    const auto *dsi_161 = buffer.data(dsi + 161);
    const auto *dsi_162 = buffer.data(dsi + 162);
    const auto *dsi_163 = buffer.data(dsi + 163);
    const auto *dsi_164 = buffer.data(dsi + 164);
    const auto *dsi_165 = buffer.data(dsi + 165);
    const auto *dsi_166 = buffer.data(dsi + 166);
    const auto *dsi_167 = buffer.data(dsi + 167);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pc_x, dsh0_83, dsh1_83, \
                         dsi_104, dsi_105, dsi_106, dsi_107, dsi_108, \
                         dsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_4 * dsh0_83[k]
                   - f_5 * dsh1_83[k]
                   + f_3 * pc_x[k] * dsi_104[k];

        t_129[k] = f_3 * pc_x[k] * dsi_105[k];

        t_130[k] = f_3 * pc_x[k] * dsi_106[k];

        t_131[k] = f_3 * pc_x[k] * dsi_107[k];

        t_132[k] = f_3 * pc_x[k] * dsi_108[k];

        t_133[k] = f_3 * pc_x[k] * dsi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, psi_49, dsh0_78, \
                         dsh1_78, dsi_105, dsi_106, dsi_110, dsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_3 * pc_x[k] * dsi_110[k];

        t_135[k] = f_3 * pc_x[k] * dsi_111[k];

        t_136[k] = f_0 * psi_49[k]
                   + f_1 * dsh0_78[k]
                   - f_2 * dsh1_78[k]
                   + f_3 * pc_y[k] * dsi_105[k];

        t_137[k] = f_3 * pc_z[k] * dsi_105[k];

        t_138[k] = f_4 * dsh0_78[k]
                   - f_5 * dsh1_78[k]
                   + f_3 * pc_z[k] * dsi_106[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_z, dsh0_79, dsh0_80, dsh0_81, dsh1_79, \
                         dsh1_80, dsh1_81, dsi_107, dsi_108, dsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_6 * dsh0_79[k]
                   - f_7 * dsh1_79[k]
                   + f_3 * pc_z[k] * dsi_107[k];

        t_140[k] = f_8 * dsh0_80[k]
                   - f_9 * dsh1_80[k]
                   + f_3 * pc_z[k] * dsi_108[k];

        t_141[k] = f_10 * dsh0_81[k]
                   - f_11 * dsh1_81[k]
                   + f_3 * pc_z[k] * dsi_109[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_y, pa_z, pc_y, pc_z, psk0_37, psk0_72, \
                         psi_55, psk1_37, psk1_72, dsh0_83, dsh1_83, \
                         dsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_0 * psi_55[k]
                   + f_3 * pc_y[k] * dsi_111[k];

        t_143[k] = f_1 * dsh0_83[k]
                   - f_2 * dsh1_83[k]
                   + f_3 * pc_z[k] * dsi_111[k];

        t_144[k] = pa_y[k] * psk0_72[k]
                   - f_12 * pc_y[k] * psk1_72[k];

        t_145[k] = pa_z[k] * psk0_37[k]
                   - f_12 * pc_z[k] * psk1_37[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_y, pa_z, pc_x, pc_y, pc_z, psk0_39, psk0_74, \
                         psk1_39, psk1_74, dsh0_88, dsh1_88, dsi_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_y[k] * psk0_74[k]
                   - f_12 * pc_y[k] * psk1_74[k];

        t_147[k] = pa_z[k] * psk0_39[k]
                   - f_12 * pc_z[k] * psk1_39[k];

        t_148[k] = f_10 * dsh0_88[k]
                   - f_11 * dsh1_88[k]
                   + f_3 * pc_x[k] * dsi_116[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_y, pa_z, pc_x, pc_y, pc_z, psk0_42, psk0_77, \
                         psk1_42, psk1_77, dsh0_91, dsh1_91, dsi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * psk0_77[k]
                   - f_12 * pc_y[k] * psk1_77[k];

        t_150[k] = pa_z[k] * psk0_42[k]
                   - f_12 * pc_z[k] * psk1_42[k];

        t_151[k] = f_8 * dsh0_91[k]
                   - f_9 * dsh1_91[k]
                   + f_3 * pc_x[k] * dsi_119[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pa_y, pa_z, pc_x, pc_y, pc_z, psk0_46, psk0_81, \
                         psk1_46, psk1_81, dsh0_92, dsh1_92, dsi_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_8 * dsh0_92[k]
                   - f_9 * dsh1_92[k]
                   + f_3 * pc_x[k] * dsi_120[k];

        t_153[k] = pa_y[k] * psk0_81[k]
                   - f_12 * pc_y[k] * psk1_81[k];

        t_154[k] = pa_z[k] * psk0_46[k]
                   - f_12 * pc_z[k] * psk1_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pc_x, dsh0_95, dsh0_96, dsh0_97, dsh1_95, \
                         dsh1_96, dsh1_97, dsi_123, dsi_124, dsi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_6 * dsh0_95[k]
                   - f_7 * dsh1_95[k]
                   + f_3 * pc_x[k] * dsi_123[k];

        t_156[k] = f_6 * dsh0_96[k]
                   - f_7 * dsh1_96[k]
                   + f_3 * pc_x[k] * dsi_124[k];

        t_157[k] = f_6 * dsh0_97[k]
                   - f_7 * dsh1_97[k]
                   + f_3 * pc_x[k] * dsi_125[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pa_z, pc_x, pc_y, pc_z, psk0_51, psk0_86, \
                         psk1_51, psk1_86, dsh0_100, dsh1_100, \
                         dsi_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = pa_y[k] * psk0_86[k]
                   - f_12 * pc_y[k] * psk1_86[k];

        t_159[k] = pa_z[k] * psk0_51[k]
                   - f_12 * pc_z[k] * psk1_51[k];

        t_160[k] = f_4 * dsh0_100[k]
                   - f_5 * dsh1_100[k]
                   + f_3 * pc_x[k] * dsi_128[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_x, dsh0_101, dsh0_102, dsh0_103, dsh1_101, \
                         dsh1_102, dsh1_103, dsi_129, dsi_130, \
                         dsi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_4 * dsh0_101[k]
                   - f_5 * dsh1_101[k]
                   + f_3 * pc_x[k] * dsi_129[k];

        t_162[k] = f_4 * dsh0_102[k]
                   - f_5 * dsh1_102[k]
                   + f_3 * pc_x[k] * dsi_130[k];

        t_163[k] = f_4 * dsh0_103[k]
                   - f_5 * dsh1_103[k]
                   + f_3 * pc_x[k] * dsi_131[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, pa_y, pc_x, pc_y, psk0_92, \
                         psk1_92, dsi_133, dsi_134, dsi_135, dsi_136, \
                         dsi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_y[k] * psk0_92[k]
                   - f_12 * pc_y[k] * psk1_92[k];

        t_165[k] = f_3 * pc_x[k] * dsi_133[k];

        t_166[k] = f_3 * pc_x[k] * dsi_134[k];

        t_167[k] = f_3 * pc_x[k] * dsi_135[k];

        t_168[k] = f_3 * pc_x[k] * dsi_136[k];

        t_169[k] = f_3 * pc_x[k] * dsi_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_z, pc_x, pc_z, psk0_64, psi_49, \
                         psk1_64, dsi_133, dsi_138, dsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_3 * pc_x[k] * dsi_138[k];

        t_171[k] = f_3 * pc_x[k] * dsi_139[k];

        t_172[k] = pa_z[k] * psk0_64[k]
                   - f_12 * pc_z[k] * psk1_64[k];

        t_173[k] = f_13 * psi_49[k]
                   + f_3 * pc_z[k] * dsi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pa_y, pc_y, psk0_102, psk0_103, psk0_104, \
                         psi_79, psi_80, psi_81, psk1_102, psk1_103, \
                         psk1_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_y[k] * psk0_102[k]
                   + f_14 * psi_79[k]
                   - f_12 * pc_y[k] * psk1_102[k];

        t_175[k] = pa_y[k] * psk0_103[k]
                   + f_15 * psi_80[k]
                   - f_12 * pc_y[k] * psk1_103[k];

        t_176[k] = pa_y[k] * psk0_104[k]
                   + f_16 * psi_81[k]
                   - f_12 * pc_y[k] * psk1_104[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, psk0_105, psk0_107, psi_82, psi_83, \
                         psk1_105, psk1_107, dsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_y[k] * psk0_105[k]
                   + f_0 * psi_82[k]
                   - f_12 * pc_y[k] * psk1_105[k];

        t_178[k] = f_13 * psi_83[k]
                   + f_3 * pc_y[k] * dsi_139[k];

        t_179[k] = pa_y[k] * psk0_107[k]
                   - f_12 * pc_y[k] * psk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, dsh0_105, dsh0_107, \
                         dsh0_108, dsh1_105, dsh1_107, dsh1_108, dsi_140, dsi_142, \
                         dsi_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * dsh0_105[k]
                   - f_2 * dsh1_105[k]
                   + f_3 * pc_x[k] * dsi_140[k];

        t_181[k] = f_3 * pc_y[k] * dsi_140[k];

        t_182[k] = f_17 * dsh0_107[k]
                   - f_18 * dsh1_107[k]
                   + f_3 * pc_x[k] * dsi_142[k];

        t_183[k] = f_10 * dsh0_108[k]
                   - f_11 * dsh1_108[k]
                   + f_3 * pc_x[k] * dsi_143[k];

        t_184[k] = f_3 * pc_y[k] * dsi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, dsh0_110, dsh0_111, dsh0_112, \
                         dsh1_110, dsh1_111, dsh1_112, dsi_145, dsi_146, \
                         dsi_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_10 * dsh0_110[k]
                   - f_11 * dsh1_110[k]
                   + f_3 * pc_x[k] * dsi_145[k];

        t_186[k] = f_8 * dsh0_111[k]
                   - f_9 * dsh1_111[k]
                   + f_3 * pc_x[k] * dsi_146[k];

        t_187[k] = f_8 * dsh0_112[k]
                   - f_9 * dsh1_112[k]
                   + f_3 * pc_x[k] * dsi_147[k];

        t_188[k] = f_3 * pc_y[k] * dsi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, dsh0_114, dsh0_115, dsh0_116, dsh1_114, \
                         dsh1_115, dsh1_116, dsi_149, dsi_150, \
                         dsi_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_8 * dsh0_114[k]
                   - f_9 * dsh1_114[k]
                   + f_3 * pc_x[k] * dsi_149[k];

        t_190[k] = f_6 * dsh0_115[k]
                   - f_7 * dsh1_115[k]
                   + f_3 * pc_x[k] * dsi_150[k];

        t_191[k] = f_6 * dsh0_116[k]
                   - f_7 * dsh1_116[k]
                   + f_3 * pc_x[k] * dsi_151[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_x, pc_y, dsh0_117, dsh0_119, dsh0_120, \
                         dsh1_117, dsh1_119, dsh1_120, dsi_149, dsi_152, dsi_154, \
                         dsi_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_6 * dsh0_117[k]
                   - f_7 * dsh1_117[k]
                   + f_3 * pc_x[k] * dsi_152[k];

        t_193[k] = f_3 * pc_y[k] * dsi_149[k];

        t_194[k] = f_6 * dsh0_119[k]
                   - f_7 * dsh1_119[k]
                   + f_3 * pc_x[k] * dsi_154[k];

        t_195[k] = f_4 * dsh0_120[k]
                   - f_5 * dsh1_120[k]
                   + f_3 * pc_x[k] * dsi_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, dsh0_121, dsh0_122, dsh0_123, \
                         dsh1_121, dsh1_122, dsh1_123, dsi_154, dsi_156, dsi_157, \
                         dsi_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_4 * dsh0_121[k]
                   - f_5 * dsh1_121[k]
                   + f_3 * pc_x[k] * dsi_156[k];

        t_197[k] = f_4 * dsh0_122[k]
                   - f_5 * dsh1_122[k]
                   + f_3 * pc_x[k] * dsi_157[k];

        t_198[k] = f_4 * dsh0_123[k]
                   - f_5 * dsh1_123[k]
                   + f_3 * pc_x[k] * dsi_158[k];

        t_199[k] = f_3 * pc_y[k] * dsi_154[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, t_205, pc_x, dsh0_125, dsh1_125, \
                         dsi_160, dsi_161, dsi_162, dsi_163, dsi_164, \
                         dsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_4 * dsh0_125[k]
                   - f_5 * dsh1_125[k]
                   + f_3 * pc_x[k] * dsi_160[k];

        t_201[k] = f_3 * pc_x[k] * dsi_161[k];

        t_202[k] = f_3 * pc_x[k] * dsi_162[k];

        t_203[k] = f_3 * pc_x[k] * dsi_163[k];

        t_204[k] = f_3 * pc_x[k] * dsi_164[k];

        t_205[k] = f_3 * pc_x[k] * dsi_165[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_x, pc_y, dsh0_120, dsh0_121, dsh1_120, \
                         dsh1_121, dsi_161, dsi_162, dsi_166, dsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_3 * pc_x[k] * dsi_166[k];

        t_207[k] = f_3 * pc_x[k] * dsi_167[k];

        t_208[k] = f_1 * dsh0_120[k]
                   - f_2 * dsh1_120[k]
                   + f_3 * pc_y[k] * dsi_161[k];

        t_209[k] = f_17 * dsh0_121[k]
                   - f_18 * dsh1_121[k]
                   + f_3 * pc_y[k] * dsi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, dsh0_122, dsh0_123, dsh0_124, dsh1_122, \
                         dsh1_123, dsh1_124, dsi_163, dsi_164, \
                         dsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * dsh0_122[k]
                   - f_11 * dsh1_122[k]
                   + f_3 * pc_y[k] * dsi_163[k];

        t_211[k] = f_8 * dsh0_123[k]
                   - f_9 * dsh1_123[k]
                   + f_3 * pc_y[k] * dsi_164[k];

        t_212[k] = f_6 * dsh0_124[k]
                   - f_7 * dsh1_124[k]
                   + f_3 * pc_y[k] * dsi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pc_y, pc_z, psi_83, dsh0_125, dsh1_125, dsi_166, \
                         dsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * dsh0_125[k]
                   - f_5 * dsh1_125[k]
                   + f_3 * pc_y[k] * dsi_166[k];

        t_214[k] = f_3 * pc_y[k] * dsi_167[k];

        t_215[k] = f_0 * psi_83[k]
                   + f_1 * dsh0_125[k]
                   - f_2 * dsh1_125[k]
                   + f_3 * pc_z[k] * dsi_167[k];
    }
}

auto
compute_prim_dsk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psk0, const size_t psi,
                                                   const size_t psk1, const size_t dsh0,
                                                   const size_t dsh1, const size_t dsi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dsk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, psk0, psi,
                                                              psk1, dsh0, dsh1, dsi, ncols,
                                                              gamma, p, q);

    compute_prim_dsk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, psk0, psi,
                                                              psk1, dsh0, dsh1, dsi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
