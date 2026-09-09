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


#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t psi0,
                                                          const size_t psh, const size_t psi1,
                                                          const size_t dsg0, const size_t dsg1,
                                                          const size_t dsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);

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
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psi0_0 = buffer.data(psi0 + 0);
    const auto *psi0_3 = buffer.data(psi0 + 3);
    const auto *psi0_5 = buffer.data(psi0 + 5);
    const auto *psi0_6 = buffer.data(psi0 + 6);
    const auto *psi0_9 = buffer.data(psi0 + 9);
    const auto *psi0_10 = buffer.data(psi0 + 10);
    const auto *psi0_14 = buffer.data(psi0 + 14);
    const auto *psi0_29 = buffer.data(psi0 + 29);
    const auto *psi0_31 = buffer.data(psi0 + 31);
    const auto *psi0_34 = buffer.data(psi0 + 34);
    const auto *psi0_38 = buffer.data(psi0 + 38);
    const auto *psi0_49 = buffer.data(psi0 + 49);
    const auto *psi0_51 = buffer.data(psi0 + 51);
    const auto *psi0_52 = buffer.data(psi0 + 52);
    const auto *psi0_53 = buffer.data(psi0 + 53);
    const auto *psi0_55 = buffer.data(psi0 + 55);
    const auto *psi0_56 = buffer.data(psi0 + 56);
    const auto *psi0_58 = buffer.data(psi0 + 58);
    const auto *psi0_61 = buffer.data(psi0 + 61);
    const auto *psi0_65 = buffer.data(psi0 + 65);
    const auto *psi0_70 = buffer.data(psi0 + 70);
    const auto *psi0_77 = buffer.data(psi0 + 77);
    const auto *psi0_78 = buffer.data(psi0 + 78);
    const auto *psi0_79 = buffer.data(psi0 + 79);
    const auto *psi0_80 = buffer.data(psi0 + 80);
    const auto *psi0_81 = buffer.data(psi0 + 81);
    const auto *psi0_83 = buffer.data(psi0 + 83);

    const auto *psh_0 = buffer.data(psh + 0);
    const auto *psh_5 = buffer.data(psh + 5);
    const auto *psh_9 = buffer.data(psh + 9);
    const auto *psh_15 = buffer.data(psh + 15);
    const auto *psh_17 = buffer.data(psh + 17);
    const auto *psh_18 = buffer.data(psh + 18);
    const auto *psh_20 = buffer.data(psh + 20);
    const auto *psh_24 = buffer.data(psh + 24);
    const auto *psh_27 = buffer.data(psh + 27);
    const auto *psh_31 = buffer.data(psh + 31);
    const auto *psh_36 = buffer.data(psh + 36);
    const auto *psh_38 = buffer.data(psh + 38);
    const auto *psh_39 = buffer.data(psh + 39);
    const auto *psh_40 = buffer.data(psh + 40);
    const auto *psh_41 = buffer.data(psh + 41);
    const auto *psh_47 = buffer.data(psh + 47);
    const auto *psh_51 = buffer.data(psh + 51);
    const auto *psh_56 = buffer.data(psh + 56);
    const auto *psh_57 = buffer.data(psh + 57);
    const auto *psh_58 = buffer.data(psh + 58);
    const auto *psh_59 = buffer.data(psh + 59);
    const auto *psh_60 = buffer.data(psh + 60);
    const auto *psh_62 = buffer.data(psh + 62);

    const auto *psi1_0 = buffer.data(psi1 + 0);
    const auto *psi1_3 = buffer.data(psi1 + 3);
    const auto *psi1_5 = buffer.data(psi1 + 5);
    const auto *psi1_6 = buffer.data(psi1 + 6);
    const auto *psi1_9 = buffer.data(psi1 + 9);
    const auto *psi1_10 = buffer.data(psi1 + 10);
    const auto *psi1_14 = buffer.data(psi1 + 14);
    const auto *psi1_29 = buffer.data(psi1 + 29);
    const auto *psi1_31 = buffer.data(psi1 + 31);
    const auto *psi1_34 = buffer.data(psi1 + 34);
    const auto *psi1_38 = buffer.data(psi1 + 38);
    const auto *psi1_49 = buffer.data(psi1 + 49);
    const auto *psi1_51 = buffer.data(psi1 + 51);
    const auto *psi1_52 = buffer.data(psi1 + 52);
    const auto *psi1_53 = buffer.data(psi1 + 53);
    const auto *psi1_55 = buffer.data(psi1 + 55);
    const auto *psi1_56 = buffer.data(psi1 + 56);
    const auto *psi1_58 = buffer.data(psi1 + 58);
    const auto *psi1_61 = buffer.data(psi1 + 61);
    const auto *psi1_65 = buffer.data(psi1 + 65);
    const auto *psi1_70 = buffer.data(psi1 + 70);
    const auto *psi1_77 = buffer.data(psi1 + 77);
    const auto *psi1_78 = buffer.data(psi1 + 78);
    const auto *psi1_79 = buffer.data(psi1 + 79);
    const auto *psi1_80 = buffer.data(psi1 + 80);
    const auto *psi1_81 = buffer.data(psi1 + 81);
    const auto *psi1_83 = buffer.data(psi1 + 83);

    const auto *dsg0_0 = buffer.data(dsg0 + 0);
    const auto *dsg0_1 = buffer.data(dsg0 + 1);
    const auto *dsg0_2 = buffer.data(dsg0 + 2);
    const auto *dsg0_3 = buffer.data(dsg0 + 3);
    const auto *dsg0_5 = buffer.data(dsg0 + 5);
    const auto *dsg0_10 = buffer.data(dsg0 + 10);
    const auto *dsg0_12 = buffer.data(dsg0 + 12);
    const auto *dsg0_13 = buffer.data(dsg0 + 13);
    const auto *dsg0_14 = buffer.data(dsg0 + 14);
    const auto *dsg0_18 = buffer.data(dsg0 + 18);
    const auto *dsg0_32 = buffer.data(dsg0 + 32);
    const auto *dsg0_34 = buffer.data(dsg0 + 34);
    const auto *dsg0_35 = buffer.data(dsg0 + 35);
    const auto *dsg0_45 = buffer.data(dsg0 + 45);
    const auto *dsg0_46 = buffer.data(dsg0 + 46);
    const auto *dsg0_48 = buffer.data(dsg0 + 48);
    const auto *dsg0_50 = buffer.data(dsg0 + 50);
    const auto *dsg0_51 = buffer.data(dsg0 + 51);
    const auto *dsg0_53 = buffer.data(dsg0 + 53);
    const auto *dsg0_54 = buffer.data(dsg0 + 54);
    const auto *dsg0_55 = buffer.data(dsg0 + 55);
    const auto *dsg0_56 = buffer.data(dsg0 + 56);
    const auto *dsg0_57 = buffer.data(dsg0 + 57);
    const auto *dsg0_58 = buffer.data(dsg0 + 58);
    const auto *dsg0_59 = buffer.data(dsg0 + 59);
    const auto *dsg0_64 = buffer.data(dsg0 + 64);
    const auto *dsg0_67 = buffer.data(dsg0 + 67);
    const auto *dsg0_68 = buffer.data(dsg0 + 68);
    const auto *dsg0_71 = buffer.data(dsg0 + 71);
    const auto *dsg0_72 = buffer.data(dsg0 + 72);
    const auto *dsg0_73 = buffer.data(dsg0 + 73);

    const auto *dsg1_0 = buffer.data(dsg1 + 0);
    const auto *dsg1_1 = buffer.data(dsg1 + 1);
    const auto *dsg1_2 = buffer.data(dsg1 + 2);
    const auto *dsg1_3 = buffer.data(dsg1 + 3);
    const auto *dsg1_5 = buffer.data(dsg1 + 5);
    const auto *dsg1_10 = buffer.data(dsg1 + 10);
    const auto *dsg1_12 = buffer.data(dsg1 + 12);
    const auto *dsg1_13 = buffer.data(dsg1 + 13);
    const auto *dsg1_14 = buffer.data(dsg1 + 14);
    const auto *dsg1_18 = buffer.data(dsg1 + 18);
    const auto *dsg1_32 = buffer.data(dsg1 + 32);
    const auto *dsg1_34 = buffer.data(dsg1 + 34);
    const auto *dsg1_35 = buffer.data(dsg1 + 35);
    const auto *dsg1_45 = buffer.data(dsg1 + 45);
    const auto *dsg1_46 = buffer.data(dsg1 + 46);
    const auto *dsg1_48 = buffer.data(dsg1 + 48);
    const auto *dsg1_50 = buffer.data(dsg1 + 50);
    const auto *dsg1_51 = buffer.data(dsg1 + 51);
    const auto *dsg1_53 = buffer.data(dsg1 + 53);
    const auto *dsg1_54 = buffer.data(dsg1 + 54);
    const auto *dsg1_55 = buffer.data(dsg1 + 55);
    const auto *dsg1_56 = buffer.data(dsg1 + 56);
    const auto *dsg1_57 = buffer.data(dsg1 + 57);
    const auto *dsg1_58 = buffer.data(dsg1 + 58);
    const auto *dsg1_59 = buffer.data(dsg1 + 59);
    const auto *dsg1_64 = buffer.data(dsg1 + 64);
    const auto *dsg1_67 = buffer.data(dsg1 + 67);
    const auto *dsg1_68 = buffer.data(dsg1 + 68);
    const auto *dsg1_71 = buffer.data(dsg1 + 71);
    const auto *dsg1_72 = buffer.data(dsg1 + 72);
    const auto *dsg1_73 = buffer.data(dsg1 + 73);

    const auto *dsh_0 = buffer.data(dsh + 0);
    const auto *dsh_1 = buffer.data(dsh + 1);
    const auto *dsh_2 = buffer.data(dsh + 2);
    const auto *dsh_3 = buffer.data(dsh + 3);
    const auto *dsh_5 = buffer.data(dsh + 5);
    const auto *dsh_6 = buffer.data(dsh + 6);
    const auto *dsh_8 = buffer.data(dsh + 8);
    const auto *dsh_9 = buffer.data(dsh + 9);
    const auto *dsh_10 = buffer.data(dsh + 10);
    const auto *dsh_14 = buffer.data(dsh + 14);
    const auto *dsh_15 = buffer.data(dsh + 15);
    const auto *dsh_17 = buffer.data(dsh + 17);
    const auto *dsh_18 = buffer.data(dsh + 18);
    const auto *dsh_19 = buffer.data(dsh + 19);
    const auto *dsh_20 = buffer.data(dsh + 20);
    const auto *dsh_21 = buffer.data(dsh + 21);
    const auto *dsh_22 = buffer.data(dsh + 22);
    const auto *dsh_24 = buffer.data(dsh + 24);
    const auto *dsh_26 = buffer.data(dsh + 26);
    const auto *dsh_27 = buffer.data(dsh + 27);
    const auto *dsh_28 = buffer.data(dsh + 28);
    const auto *dsh_30 = buffer.data(dsh + 30);
    const auto *dsh_31 = buffer.data(dsh + 31);
    const auto *dsh_36 = buffer.data(dsh + 36);
    const auto *dsh_38 = buffer.data(dsh + 38);
    const auto *dsh_39 = buffer.data(dsh + 39);
    const auto *dsh_40 = buffer.data(dsh + 40);
    const auto *dsh_41 = buffer.data(dsh + 41);
    const auto *dsh_42 = buffer.data(dsh + 42);
    const auto *dsh_44 = buffer.data(dsh + 44);
    const auto *dsh_46 = buffer.data(dsh + 46);
    const auto *dsh_47 = buffer.data(dsh + 47);
    const auto *dsh_49 = buffer.data(dsh + 49);
    const auto *dsh_50 = buffer.data(dsh + 50);
    const auto *dsh_51 = buffer.data(dsh + 51);
    const auto *dsh_56 = buffer.data(dsh + 56);
    const auto *dsh_57 = buffer.data(dsh + 57);
    const auto *dsh_58 = buffer.data(dsh + 58);
    const auto *dsh_59 = buffer.data(dsh + 59);
    const auto *dsh_60 = buffer.data(dsh + 60);
    const auto *dsh_62 = buffer.data(dsh + 62);
    const auto *dsh_63 = buffer.data(dsh + 63);
    const auto *dsh_64 = buffer.data(dsh + 64);
    const auto *dsh_66 = buffer.data(dsh + 66);
    const auto *dsh_68 = buffer.data(dsh + 68);
    const auto *dsh_69 = buffer.data(dsh + 69);
    const auto *dsh_71 = buffer.data(dsh + 71);
    const auto *dsh_72 = buffer.data(dsh + 72);
    const auto *dsh_73 = buffer.data(dsh + 73);
    const auto *dsh_75 = buffer.data(dsh + 75);
    const auto *dsh_76 = buffer.data(dsh + 76);
    const auto *dsh_77 = buffer.data(dsh + 77);
    const auto *dsh_78 = buffer.data(dsh + 78);
    const auto *dsh_79 = buffer.data(dsh + 79);
    const auto *dsh_80 = buffer.data(dsh + 80);
    const auto *dsh_81 = buffer.data(dsh + 81);
    const auto *dsh_82 = buffer.data(dsh + 82);
    const auto *dsh_83 = buffer.data(dsh + 83);
    const auto *dsh_88 = buffer.data(dsh + 88);
    const auto *dsh_91 = buffer.data(dsh + 91);
    const auto *dsh_92 = buffer.data(dsh + 92);
    const auto *dsh_95 = buffer.data(dsh + 95);
    const auto *dsh_96 = buffer.data(dsh + 96);
    const auto *dsh_97 = buffer.data(dsh + 97);
    const auto *dsh_99 = buffer.data(dsh + 99);
    const auto *dsh_100 = buffer.data(dsh + 100);
    const auto *dsh_101 = buffer.data(dsh + 101);
    const auto *dsh_102 = buffer.data(dsh + 102);
    const auto *dsh_103 = buffer.data(dsh + 103);
    const auto *dsh_104 = buffer.data(dsh + 104);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, psh_0, dsg0_0, \
                         dsg1_0, dsh_0, dsh_1, dsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psh_0[k]
                 + f_1 * dsg0_0[k]
                 - f_2 * dsg1_0[k]
                 + f_3 * pc_x[k] * dsh_0[k];

        t_1[k] = f_3 * pc_y[k] * dsh_0[k];

        t_2[k] = f_3 * pc_z[k] * dsh_0[k];

        t_3[k] = f_4 * dsg0_0[k]
                 - f_5 * dsg1_0[k]
                 + f_3 * pc_y[k] * dsh_1[k];

        t_4[k] = f_3 * pc_y[k] * dsh_2[k];

        t_5[k] = f_4 * dsg0_0[k]
                 - f_5 * dsg1_0[k]
                 + f_3 * pc_z[k] * dsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, dsg0_1, dsg0_2, dsg0_3, dsg1_1, \
                         dsg1_2, dsg1_3, dsh_3, dsh_5, dsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * dsg0_1[k]
                 - f_7 * dsg1_1[k]
                 + f_3 * pc_y[k] * dsh_3[k];

        t_7[k] = f_3 * pc_z[k] * dsh_3[k];

        t_8[k] = f_3 * pc_y[k] * dsh_5[k];

        t_9[k] = f_6 * dsg0_2[k]
                 - f_7 * dsg1_2[k]
                 + f_3 * pc_z[k] * dsh_5[k];

        t_10[k] = f_8 * dsg0_3[k]
                  - f_9 * dsg1_3[k]
                  + f_3 * pc_y[k] * dsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, psh_15, dsg0_5, \
                         dsg1_5, dsh_6, dsh_8, dsh_9, dsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * dsh_6[k];

        t_12[k] = f_4 * dsg0_5[k]
                  - f_5 * dsg1_5[k]
                  + f_3 * pc_y[k] * dsh_8[k];

        t_13[k] = f_3 * pc_y[k] * dsh_9[k];

        t_14[k] = f_8 * dsg0_5[k]
                  - f_9 * dsg1_5[k]
                  + f_3 * pc_z[k] * dsh_9[k];

        t_15[k] = f_0 * psh_15[k]
                  + f_3 * pc_x[k] * dsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, psh_17, psh_18, \
                         psh_20, dsh_10, dsh_14, dsh_17, dsh_18, \
                         dsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * dsh_10[k];

        t_17[k] = f_0 * psh_17[k]
                  + f_3 * pc_x[k] * dsh_17[k];

        t_18[k] = f_0 * psh_18[k]
                  + f_3 * pc_x[k] * dsh_18[k];

        t_19[k] = f_3 * pc_y[k] * dsh_14[k];

        t_20[k] = f_0 * psh_20[k]
                  + f_3 * pc_x[k] * dsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, dsg0_10, dsg0_12, dsg0_13, \
                         dsg1_10, dsg1_12, dsg1_13, dsh_15, dsh_17, \
                         dsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * dsg0_10[k]
                  - f_2 * dsg1_10[k]
                  + f_3 * pc_y[k] * dsh_15[k];

        t_22[k] = f_3 * pc_z[k] * dsh_15[k];

        t_23[k] = f_8 * dsg0_12[k]
                  - f_9 * dsg1_12[k]
                  + f_3 * pc_y[k] * dsh_17[k];

        t_24[k] = f_6 * dsg0_13[k]
                  - f_7 * dsg1_13[k]
                  + f_3 * pc_y[k] * dsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, psi0_0, psh_0, \
                         psi1_0, dsg0_14, dsg1_14, dsh_19, dsh_20, \
                         dsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * dsg0_14[k]
                  - f_5 * dsg1_14[k]
                  + f_3 * pc_y[k] * dsh_19[k];

        t_26[k] = f_3 * pc_y[k] * dsh_20[k];

        t_27[k] = f_1 * dsg0_14[k]
                  - f_2 * dsg1_14[k]
                  + f_3 * pc_z[k] * dsh_20[k];

        t_28[k] = pa_y[k] * psi0_0[k]
                  - f_10 * pc_y[k] * psi1_0[k];

        t_29[k] = f_11 * psh_0[k]
                  + f_3 * pc_y[k] * dsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pc_x, pc_y, pc_z, psi0_5, \
                         psi0_31, psh_24, psi1_5, psi1_31, dsh_21, \
                         dsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * dsh_21[k];

        t_31[k] = pa_x[k] * psi0_31[k]
                  + f_12 * psh_24[k]
                  - f_10 * pc_x[k] * psi1_31[k];

        t_32[k] = f_3 * pc_z[k] * dsh_22[k];

        t_33[k] = pa_y[k] * psi0_5[k]
                  - f_10 * pc_y[k] * psi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pc_x, pc_y, pc_z, psi0_34, psh_5, psh_27, \
                         psi1_34, dsh_24, dsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * psi0_34[k]
                  + f_13 * psh_27[k]
                  - f_10 * pc_x[k] * psi1_34[k];

        t_35[k] = f_3 * pc_z[k] * dsh_24[k];

        t_36[k] = f_11 * psh_5[k]
                  + f_3 * pc_y[k] * dsh_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pa_y, pc_x, pc_y, pc_z, psi0_9, psi0_38, \
                         psh_31, psi1_9, psi1_38, dsh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * psi0_9[k]
                  - f_10 * pc_y[k] * psi1_9[k];

        t_38[k] = pa_x[k] * psi0_38[k]
                  + f_0 * psh_31[k]
                  - f_10 * pc_x[k] * psi1_38[k];

        t_39[k] = f_3 * pc_z[k] * dsh_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pc_y, pc_z, psi0_14, psh_9, psi1_14, dsg0_18, \
                         dsg1_18, dsh_28, dsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_4 * dsg0_18[k]
                  - f_5 * dsg1_18[k]
                  + f_3 * pc_z[k] * dsh_28[k];

        t_41[k] = f_11 * psh_9[k]
                  + f_3 * pc_y[k] * dsh_30[k];

        t_42[k] = pa_y[k] * psi0_14[k]
                  - f_10 * pc_y[k] * psi1_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pc_x, pc_z, psh_36, psh_38, psh_39, \
                         psh_40, dsh_31, dsh_36, dsh_38, dsh_39, \
                         dsh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * psh_36[k]
                  + f_3 * pc_x[k] * dsh_36[k];

        t_44[k] = f_3 * pc_z[k] * dsh_31[k];

        t_45[k] = f_11 * psh_38[k]
                  + f_3 * pc_x[k] * dsh_38[k];

        t_46[k] = f_11 * psh_39[k]
                  + f_3 * pc_x[k] * dsh_39[k];

        t_47[k] = f_11 * psh_40[k]
                  + f_3 * pc_x[k] * dsh_40[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pc_x, pc_z, psi0_49, psi0_51, psh_41, \
                         psi1_49, psi1_51, dsh_36, dsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * psh_41[k]
                  + f_3 * pc_x[k] * dsh_41[k];

        t_49[k] = pa_x[k] * psi0_49[k]
                  - f_10 * pc_x[k] * psi1_49[k];

        t_50[k] = f_3 * pc_z[k] * dsh_36[k];

        t_51[k] = pa_x[k] * psi0_51[k]
                  - f_10 * pc_x[k] * psi1_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pc_x, pc_y, psi0_52, psi0_53, psi0_55, \
                         psh_20, psi1_52, psi1_53, psi1_55, dsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_x[k] * psi0_52[k]
                  - f_10 * pc_x[k] * psi1_52[k];

        t_53[k] = pa_x[k] * psi0_53[k]
                  - f_10 * pc_x[k] * psi1_53[k];

        t_54[k] = f_11 * psh_20[k]
                  + f_3 * pc_y[k] * dsh_41[k];

        t_55[k] = pa_x[k] * psi0_55[k]
                  - f_10 * pc_x[k] * psi1_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_z, pc_y, pc_z, psi0_0, psi0_3, \
                         psh_0, psi1_0, psi1_3, dsh_42, dsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_z[k] * psi0_0[k]
                  - f_10 * pc_z[k] * psi1_0[k];

        t_57[k] = f_3 * pc_y[k] * dsh_42[k];

        t_58[k] = f_11 * psh_0[k]
                  + f_3 * pc_z[k] * dsh_42[k];

        t_59[k] = pa_z[k] * psi0_3[k]
                  - f_10 * pc_z[k] * psi1_3[k];

        t_60[k] = f_3 * pc_y[k] * dsh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_x, pa_z, pc_x, pc_y, pc_z, psi0_6, psi0_61, \
                         psh_47, psi1_6, psi1_61, dsg0_32, dsg1_32, \
                         dsh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_x[k] * psi0_61[k]
                  + f_12 * psh_47[k]
                  - f_10 * pc_x[k] * psi1_61[k];

        t_62[k] = pa_z[k] * psi0_6[k]
                  - f_10 * pc_z[k] * psi1_6[k];

        t_63[k] = f_4 * dsg0_32[k]
                  - f_5 * dsg1_32[k]
                  + f_3 * pc_y[k] * dsh_46[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_x, pa_z, pc_x, pc_y, pc_z, psi0_10, psi0_65, \
                         psh_51, psi1_10, psi1_65, dsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * pc_y[k] * dsh_47[k];

        t_65[k] = pa_x[k] * psi0_65[k]
                  + f_13 * psh_51[k]
                  - f_10 * pc_x[k] * psi1_65[k];

        t_66[k] = pa_z[k] * psi0_10[k]
                  - f_10 * pc_z[k] * psi1_10[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, dsg0_34, dsg0_35, dsg1_34, dsg1_35, dsh_49, \
                         dsh_50, dsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_6 * dsg0_34[k]
                  - f_7 * dsg1_34[k]
                  + f_3 * pc_y[k] * dsh_49[k];

        t_68[k] = f_4 * dsg0_35[k]
                  - f_5 * dsg1_35[k]
                  + f_3 * pc_y[k] * dsh_50[k];

        t_69[k] = f_3 * pc_y[k] * dsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pc_x, psi0_70, psh_56, psh_57, psh_58, \
                         psh_59, psi1_70, dsh_57, dsh_58, dsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * psi0_70[k]
                  + f_0 * psh_56[k]
                  - f_10 * pc_x[k] * psi1_70[k];

        t_71[k] = f_11 * psh_57[k]
                  + f_3 * pc_x[k] * dsh_57[k];

        t_72[k] = f_11 * psh_58[k]
                  + f_3 * pc_x[k] * dsh_58[k];

        t_73[k] = f_11 * psh_59[k]
                  + f_3 * pc_x[k] * dsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_y, psi0_77, psh_60, psh_62, \
                         psi1_77, dsh_56, dsh_60, dsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * psh_60[k]
                  + f_3 * pc_x[k] * dsh_60[k];

        t_75[k] = f_3 * pc_y[k] * dsh_56[k];

        t_76[k] = f_11 * psh_62[k]
                  + f_3 * pc_x[k] * dsh_62[k];

        t_77[k] = pa_x[k] * psi0_77[k]
                  - f_10 * pc_x[k] * psi1_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pc_x, psi0_78, psi0_79, psi0_80, \
                         psi0_81, psi1_78, psi1_79, psi1_80, psi1_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * psi0_78[k]
                  - f_10 * pc_x[k] * psi1_78[k];

        t_79[k] = pa_x[k] * psi0_79[k]
                  - f_10 * pc_x[k] * psi1_79[k];

        t_80[k] = pa_x[k] * psi0_80[k]
                  - f_10 * pc_x[k] * psi1_80[k];

        t_81[k] = pa_x[k] * psi0_81[k]
                  - f_10 * pc_x[k] * psi1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pc_x, pc_y, psi0_83, psi1_83, dsg0_45, \
                         dsg0_46, dsg1_45, dsg1_46, dsh_62, dsh_63, \
                         dsh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_y[k] * dsh_62[k];

        t_83[k] = pa_x[k] * psi0_83[k]
                  - f_10 * pc_x[k] * psi1_83[k];

        t_84[k] = f_1 * dsg0_45[k]
                  - f_2 * dsg1_45[k]
                  + f_3 * pc_x[k] * dsh_63[k];

        t_85[k] = f_14 * dsg0_46[k]
                  - f_15 * dsg1_46[k]
                  + f_3 * pc_x[k] * dsh_64[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_z, dsg0_48, dsg0_50, dsg1_48, \
                         dsg1_50, dsh_63, dsh_64, dsh_66, dsh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * pc_z[k] * dsh_63[k];

        t_87[k] = f_8 * dsg0_48[k]
                  - f_9 * dsg1_48[k]
                  + f_3 * pc_x[k] * dsh_66[k];

        t_88[k] = f_3 * pc_z[k] * dsh_64[k];

        t_89[k] = f_8 * dsg0_50[k]
                  - f_9 * dsg1_50[k]
                  + f_3 * pc_x[k] * dsh_68[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_z, dsg0_51, dsg0_53, dsg0_54, \
                         dsg1_51, dsg1_53, dsg1_54, dsh_66, dsh_69, dsh_71, \
                         dsh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_6 * dsg0_51[k]
                  - f_7 * dsg1_51[k]
                  + f_3 * pc_x[k] * dsh_69[k];

        t_91[k] = f_3 * pc_z[k] * dsh_66[k];

        t_92[k] = f_6 * dsg0_53[k]
                  - f_7 * dsg1_53[k]
                  + f_3 * pc_x[k] * dsh_71[k];

        t_93[k] = f_6 * dsg0_54[k]
                  - f_7 * dsg1_54[k]
                  + f_3 * pc_x[k] * dsh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pc_x, pc_z, dsg0_55, dsg0_57, dsg0_58, \
                         dsg1_55, dsg1_57, dsg1_58, dsh_69, dsh_73, dsh_75, \
                         dsh_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_4 * dsg0_55[k]
                  - f_5 * dsg1_55[k]
                  + f_3 * pc_x[k] * dsh_73[k];

        t_95[k] = f_3 * pc_z[k] * dsh_69[k];

        t_96[k] = f_4 * dsg0_57[k]
                  - f_5 * dsg1_57[k]
                  + f_3 * pc_x[k] * dsh_75[k];

        t_97[k] = f_4 * dsg0_58[k]
                  - f_5 * dsg1_58[k]
                  + f_3 * pc_x[k] * dsh_76[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, pc_x, dsg0_59, dsg1_59, \
                         dsh_77, dsh_78, dsh_79, dsh_80, dsh_81, \
                         dsh_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_4 * dsg0_59[k]
                  - f_5 * dsg1_59[k]
                  + f_3 * pc_x[k] * dsh_77[k];

        t_99[k] = f_3 * pc_x[k] * dsh_78[k];

        t_100[k] = f_3 * pc_x[k] * dsh_79[k];

        t_101[k] = f_3 * pc_x[k] * dsh_80[k];

        t_102[k] = f_3 * pc_x[k] * dsh_81[k];

        t_103[k] = f_3 * pc_x[k] * dsh_82[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, psh_36, dsg0_55, \
                         dsg1_55, dsh_78, dsh_79, dsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * pc_x[k] * dsh_83[k];

        t_105[k] = f_0 * psh_36[k]
                   + f_1 * dsg0_55[k]
                   - f_2 * dsg1_55[k]
                   + f_3 * pc_y[k] * dsh_78[k];

        t_106[k] = f_3 * pc_z[k] * dsh_78[k];

        t_107[k] = f_4 * dsg0_55[k]
                   - f_5 * dsg1_55[k]
                   + f_3 * pc_z[k] * dsh_79[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, psh_41, dsg0_56, dsg0_57, \
                         dsg0_59, dsg1_56, dsg1_57, dsg1_59, dsh_80, dsh_81, \
                         dsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_6 * dsg0_56[k]
                   - f_7 * dsg1_56[k]
                   + f_3 * pc_z[k] * dsh_80[k];

        t_109[k] = f_8 * dsg0_57[k]
                   - f_9 * dsg1_57[k]
                   + f_3 * pc_z[k] * dsh_81[k];

        t_110[k] = f_0 * psh_41[k]
                   + f_3 * pc_y[k] * dsh_83[k];

        t_111[k] = f_1 * dsg0_59[k]
                   - f_2 * dsg1_59[k]
                   + f_3 * pc_z[k] * dsh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_y, pa_z, pc_y, pc_z, psi0_29, psi0_31, \
                         psi0_56, psi0_58, psi1_29, psi1_31, psi1_56, \
                         psi1_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pa_y[k] * psi0_56[k]
                   - f_10 * pc_y[k] * psi1_56[k];

        t_113[k] = pa_z[k] * psi0_29[k]
                   - f_10 * pc_z[k] * psi1_29[k];

        t_114[k] = pa_y[k] * psi0_58[k]
                   - f_10 * pc_y[k] * psi1_58[k];

        t_115[k] = pa_z[k] * psi0_31[k]
                   - f_10 * pc_z[k] * psi1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pa_y, pa_z, pc_x, pc_y, pc_z, psi0_34, psi0_61, \
                         psi1_34, psi1_61, dsg0_64, dsg1_64, dsh_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_8 * dsg0_64[k]
                   - f_9 * dsg1_64[k]
                   + f_3 * pc_x[k] * dsh_88[k];

        t_117[k] = pa_y[k] * psi0_61[k]
                   - f_10 * pc_y[k] * psi1_61[k];

        t_118[k] = pa_z[k] * psi0_34[k]
                   - f_10 * pc_z[k] * psi1_34[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pa_y, pc_x, pc_y, psi0_65, psi1_65, dsg0_67, \
                         dsg0_68, dsg1_67, dsg1_68, dsh_91, dsh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_6 * dsg0_67[k]
                   - f_7 * dsg1_67[k]
                   + f_3 * pc_x[k] * dsh_91[k];

        t_120[k] = f_6 * dsg0_68[k]
                   - f_7 * dsg1_68[k]
                   + f_3 * pc_x[k] * dsh_92[k];

        t_121[k] = pa_y[k] * psi0_65[k]
                   - f_10 * pc_y[k] * psi1_65[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_z, pc_x, pc_z, psi0_38, psi1_38, dsg0_71, \
                         dsg0_72, dsg1_71, dsg1_72, dsh_95, dsh_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = pa_z[k] * psi0_38[k]
                   - f_10 * pc_z[k] * psi1_38[k];

        t_123[k] = f_4 * dsg0_71[k]
                   - f_5 * dsg1_71[k]
                   + f_3 * pc_x[k] * dsh_95[k];

        t_124[k] = f_4 * dsg0_72[k]
                   - f_5 * dsg1_72[k]
                   + f_3 * pc_x[k] * dsh_96[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_y, pc_x, pc_y, psi0_70, \
                         psi1_70, dsg0_73, dsg1_73, dsh_97, dsh_99, dsh_100, \
                         dsh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * dsg0_73[k]
                   - f_5 * dsg1_73[k]
                   + f_3 * pc_x[k] * dsh_97[k];

        t_126[k] = pa_y[k] * psi0_70[k]
                   - f_10 * pc_y[k] * psi1_70[k];

        t_127[k] = f_3 * pc_x[k] * dsh_99[k];

        t_128[k] = f_3 * pc_x[k] * dsh_100[k];

        t_129[k] = f_3 * pc_x[k] * dsh_101[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_z, pc_x, pc_z, psi0_49, psh_36, \
                         psi1_49, dsh_99, dsh_102, dsh_103, dsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_x[k] * dsh_102[k];

        t_131[k] = f_3 * pc_x[k] * dsh_103[k];

        t_132[k] = f_3 * pc_x[k] * dsh_104[k];

        t_133[k] = pa_z[k] * psi0_49[k]
                   - f_10 * pc_z[k] * psi1_49[k];

        t_134[k] = f_11 * psh_36[k]
                   + f_3 * pc_z[k] * dsh_99[k];
    }
}

static auto
compute_prim_dsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t psi0,
                                                          const size_t psh, const size_t psi1,
                                                          const size_t dsg0, const size_t dsg1,
                                                          const size_t dsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psi0_79 = buffer.data(psi0 + 79);
    const auto *psi0_80 = buffer.data(psi0 + 80);
    const auto *psi0_81 = buffer.data(psi0 + 81);
    const auto *psi0_83 = buffer.data(psi0 + 83);

    const auto *psh_59 = buffer.data(psh + 59);
    const auto *psh_60 = buffer.data(psh + 60);
    const auto *psh_61 = buffer.data(psh + 61);
    const auto *psh_62 = buffer.data(psh + 62);

    const auto *psi1_79 = buffer.data(psi1 + 79);
    const auto *psi1_80 = buffer.data(psi1 + 80);
    const auto *psi1_81 = buffer.data(psi1 + 81);
    const auto *psi1_83 = buffer.data(psi1 + 83);

    const auto *dsg0_75 = buffer.data(dsg0 + 75);
    const auto *dsg0_77 = buffer.data(dsg0 + 77);
    const auto *dsg0_78 = buffer.data(dsg0 + 78);
    const auto *dsg0_80 = buffer.data(dsg0 + 80);
    const auto *dsg0_81 = buffer.data(dsg0 + 81);
    const auto *dsg0_82 = buffer.data(dsg0 + 82);
    const auto *dsg0_84 = buffer.data(dsg0 + 84);
    const auto *dsg0_85 = buffer.data(dsg0 + 85);
    const auto *dsg0_86 = buffer.data(dsg0 + 86);
    const auto *dsg0_87 = buffer.data(dsg0 + 87);
    const auto *dsg0_88 = buffer.data(dsg0 + 88);
    const auto *dsg0_89 = buffer.data(dsg0 + 89);

    const auto *dsg1_75 = buffer.data(dsg1 + 75);
    const auto *dsg1_77 = buffer.data(dsg1 + 77);
    const auto *dsg1_78 = buffer.data(dsg1 + 78);
    const auto *dsg1_80 = buffer.data(dsg1 + 80);
    const auto *dsg1_81 = buffer.data(dsg1 + 81);
    const auto *dsg1_82 = buffer.data(dsg1 + 82);
    const auto *dsg1_84 = buffer.data(dsg1 + 84);
    const auto *dsg1_85 = buffer.data(dsg1 + 85);
    const auto *dsg1_86 = buffer.data(dsg1 + 86);
    const auto *dsg1_87 = buffer.data(dsg1 + 87);
    const auto *dsg1_88 = buffer.data(dsg1 + 88);
    const auto *dsg1_89 = buffer.data(dsg1 + 89);

    const auto *dsh_104 = buffer.data(dsh + 104);
    const auto *dsh_105 = buffer.data(dsh + 105);
    const auto *dsh_107 = buffer.data(dsh + 107);
    const auto *dsh_108 = buffer.data(dsh + 108);
    const auto *dsh_110 = buffer.data(dsh + 110);
    const auto *dsh_111 = buffer.data(dsh + 111);
    const auto *dsh_112 = buffer.data(dsh + 112);
    const auto *dsh_114 = buffer.data(dsh + 114);
    const auto *dsh_115 = buffer.data(dsh + 115);
    const auto *dsh_116 = buffer.data(dsh + 116);
    const auto *dsh_117 = buffer.data(dsh + 117);
    const auto *dsh_119 = buffer.data(dsh + 119);
    const auto *dsh_120 = buffer.data(dsh + 120);
    const auto *dsh_121 = buffer.data(dsh + 121);
    const auto *dsh_122 = buffer.data(dsh + 122);
    const auto *dsh_123 = buffer.data(dsh + 123);
    const auto *dsh_124 = buffer.data(dsh + 124);
    const auto *dsh_125 = buffer.data(dsh + 125);

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pc_y, psi0_79, psi0_80, psi0_81, psh_59, \
                         psh_60, psh_61, psi1_79, psi1_80, psi1_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pa_y[k] * psi0_79[k]
                   + f_12 * psh_59[k]
                   - f_10 * pc_y[k] * psi1_79[k];

        t_136[k] = pa_y[k] * psi0_80[k]
                   + f_13 * psh_60[k]
                   - f_10 * pc_y[k] * psi1_80[k];

        t_137[k] = pa_y[k] * psi0_81[k]
                   + f_0 * psh_61[k]
                   - f_10 * pc_y[k] * psi1_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pc_x, pc_y, psi0_83, psh_62, \
                         psi1_83, dsg0_75, dsg1_75, dsh_104, dsh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_11 * psh_62[k]
                   + f_3 * pc_y[k] * dsh_104[k];

        t_139[k] = pa_y[k] * psi0_83[k]
                   - f_10 * pc_y[k] * psi1_83[k];

        t_140[k] = f_1 * dsg0_75[k]
                   - f_2 * dsg1_75[k]
                   + f_3 * pc_x[k] * dsh_105[k];

        t_141[k] = f_3 * pc_y[k] * dsh_105[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_y, dsg0_77, dsg0_78, dsg0_80, \
                         dsg1_77, dsg1_78, dsg1_80, dsh_107, dsh_108, \
                         dsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * dsg0_77[k]
                   - f_15 * dsg1_77[k]
                   + f_3 * pc_x[k] * dsh_107[k];

        t_143[k] = f_8 * dsg0_78[k]
                   - f_9 * dsg1_78[k]
                   + f_3 * pc_x[k] * dsh_108[k];

        t_144[k] = f_3 * pc_y[k] * dsh_107[k];

        t_145[k] = f_8 * dsg0_80[k]
                   - f_9 * dsg1_80[k]
                   + f_3 * pc_x[k] * dsh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_x, pc_y, dsg0_81, dsg0_82, dsg0_84, \
                         dsg1_81, dsg1_82, dsg1_84, dsh_110, dsh_111, dsh_112, \
                         dsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * dsg0_81[k]
                   - f_7 * dsg1_81[k]
                   + f_3 * pc_x[k] * dsh_111[k];

        t_147[k] = f_6 * dsg0_82[k]
                   - f_7 * dsg1_82[k]
                   + f_3 * pc_x[k] * dsh_112[k];

        t_148[k] = f_3 * pc_y[k] * dsh_110[k];

        t_149[k] = f_6 * dsg0_84[k]
                   - f_7 * dsg1_84[k]
                   + f_3 * pc_x[k] * dsh_114[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, dsg0_85, dsg0_86, dsg0_87, \
                         dsg1_85, dsg1_86, dsg1_87, dsh_114, dsh_115, dsh_116, \
                         dsh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_4 * dsg0_85[k]
                   - f_5 * dsg1_85[k]
                   + f_3 * pc_x[k] * dsh_115[k];

        t_151[k] = f_4 * dsg0_86[k]
                   - f_5 * dsg1_86[k]
                   + f_3 * pc_x[k] * dsh_116[k];

        t_152[k] = f_4 * dsg0_87[k]
                   - f_5 * dsg1_87[k]
                   + f_3 * pc_x[k] * dsh_117[k];

        t_153[k] = f_3 * pc_y[k] * dsh_114[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, t_159, pc_x, dsg0_89, dsg1_89, \
                         dsh_119, dsh_120, dsh_121, dsh_122, dsh_123, \
                         dsh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_4 * dsg0_89[k]
                   - f_5 * dsg1_89[k]
                   + f_3 * pc_x[k] * dsh_119[k];

        t_155[k] = f_3 * pc_x[k] * dsh_120[k];

        t_156[k] = f_3 * pc_x[k] * dsh_121[k];

        t_157[k] = f_3 * pc_x[k] * dsh_122[k];

        t_158[k] = f_3 * pc_x[k] * dsh_123[k];

        t_159[k] = f_3 * pc_x[k] * dsh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, dsg0_85, dsg0_86, dsg0_87, \
                         dsg1_85, dsg1_86, dsg1_87, dsh_120, dsh_121, dsh_122, \
                         dsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_3 * pc_x[k] * dsh_125[k];

        t_161[k] = f_1 * dsg0_85[k]
                   - f_2 * dsg1_85[k]
                   + f_3 * pc_y[k] * dsh_120[k];

        t_162[k] = f_14 * dsg0_86[k]
                   - f_15 * dsg1_86[k]
                   + f_3 * pc_y[k] * dsh_121[k];

        t_163[k] = f_8 * dsg0_87[k]
                   - f_9 * dsg1_87[k]
                   + f_3 * pc_y[k] * dsh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, psh_62, dsg0_88, dsg0_89, \
                         dsg1_88, dsg1_89, dsh_123, dsh_124, dsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * dsg0_88[k]
                   - f_7 * dsg1_88[k]
                   + f_3 * pc_y[k] * dsh_123[k];

        t_165[k] = f_4 * dsg0_89[k]
                   - f_5 * dsg1_89[k]
                   + f_3 * pc_y[k] * dsh_124[k];

        t_166[k] = f_3 * pc_y[k] * dsh_125[k];

        t_167[k] = f_0 * psh_62[k]
                   + f_1 * dsg0_89[k]
                   - f_2 * dsg1_89[k]
                   + f_3 * pc_z[k] * dsh_125[k];
    }
}

auto
compute_prim_dsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psi0, const size_t psh,
                                                   const size_t psi1, const size_t dsg0,
                                                   const size_t dsg1, const size_t dsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, psi0, psh,
                                                              psi1, dsg0, dsg1, dsh, ncols,
                                                              gamma, p, q);

    compute_prim_dsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, psi0, psh,
                                                              psi1, dsg0, dsg1, dsh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
