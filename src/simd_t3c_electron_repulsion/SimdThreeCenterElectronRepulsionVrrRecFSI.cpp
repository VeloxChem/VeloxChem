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


#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsi0,
                                                          const size_t dsh, const size_t dsi1,
                                                          const size_t fsg0, const size_t fsg1,
                                                          const size_t fsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);
    const auto f_16 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsi0_0 = buffer.data(dsi0 + 0);
    const auto *dsi0_3 = buffer.data(dsi0 + 3);
    const auto *dsi0_5 = buffer.data(dsi0 + 5);
    const auto *dsi0_6 = buffer.data(dsi0 + 6);
    const auto *dsi0_9 = buffer.data(dsi0 + 9);
    const auto *dsi0_10 = buffer.data(dsi0 + 10);
    const auto *dsi0_14 = buffer.data(dsi0 + 14);
    const auto *dsi0_21 = buffer.data(dsi0 + 21);
    const auto *dsi0_27 = buffer.data(dsi0 + 27);
    const auto *dsi0_31 = buffer.data(dsi0 + 31);
    const auto *dsi0_34 = buffer.data(dsi0 + 34);
    const auto *dsi0_38 = buffer.data(dsi0 + 38);
    const auto *dsi0_56 = buffer.data(dsi0 + 56);
    const auto *dsi0_61 = buffer.data(dsi0 + 61);
    const auto *dsi0_65 = buffer.data(dsi0 + 65);
    const auto *dsi0_70 = buffer.data(dsi0 + 70);
    const auto *dsi0_84 = buffer.data(dsi0 + 84);
    const auto *dsi0_87 = buffer.data(dsi0 + 87);
    const auto *dsi0_90 = buffer.data(dsi0 + 90);
    const auto *dsi0_94 = buffer.data(dsi0 + 94);
    const auto *dsi0_105 = buffer.data(dsi0 + 105);
    const auto *dsi0_107 = buffer.data(dsi0 + 107);
    const auto *dsi0_108 = buffer.data(dsi0 + 108);
    const auto *dsi0_109 = buffer.data(dsi0 + 109);
    const auto *dsi0_111 = buffer.data(dsi0 + 111);
    const auto *dsi0_124 = buffer.data(dsi0 + 124);

    const auto *dsh_0 = buffer.data(dsh + 0);
    const auto *dsh_1 = buffer.data(dsh + 1);
    const auto *dsh_2 = buffer.data(dsh + 2);
    const auto *dsh_3 = buffer.data(dsh + 3);
    const auto *dsh_5 = buffer.data(dsh + 5);
    const auto *dsh_6 = buffer.data(dsh + 6);
    const auto *dsh_9 = buffer.data(dsh + 9);
    const auto *dsh_15 = buffer.data(dsh + 15);
    const auto *dsh_17 = buffer.data(dsh + 17);
    const auto *dsh_18 = buffer.data(dsh + 18);
    const auto *dsh_20 = buffer.data(dsh + 20);
    const auto *dsh_21 = buffer.data(dsh + 21);
    const auto *dsh_24 = buffer.data(dsh + 24);
    const auto *dsh_26 = buffer.data(dsh + 26);
    const auto *dsh_27 = buffer.data(dsh + 27);
    const auto *dsh_30 = buffer.data(dsh + 30);
    const auto *dsh_36 = buffer.data(dsh + 36);
    const auto *dsh_38 = buffer.data(dsh + 38);
    const auto *dsh_39 = buffer.data(dsh + 39);
    const auto *dsh_40 = buffer.data(dsh + 40);
    const auto *dsh_41 = buffer.data(dsh + 41);
    const auto *dsh_42 = buffer.data(dsh + 42);
    const auto *dsh_44 = buffer.data(dsh + 44);
    const auto *dsh_47 = buffer.data(dsh + 47);
    const auto *dsh_51 = buffer.data(dsh + 51);
    const auto *dsh_57 = buffer.data(dsh + 57);
    const auto *dsh_58 = buffer.data(dsh + 58);
    const auto *dsh_59 = buffer.data(dsh + 59);
    const auto *dsh_60 = buffer.data(dsh + 60);
    const auto *dsh_62 = buffer.data(dsh + 62);
    const auto *dsh_63 = buffer.data(dsh + 63);
    const auto *dsh_66 = buffer.data(dsh + 66);
    const auto *dsh_69 = buffer.data(dsh + 69);
    const auto *dsh_73 = buffer.data(dsh + 73);
    const auto *dsh_78 = buffer.data(dsh + 78);
    const auto *dsh_80 = buffer.data(dsh + 80);
    const auto *dsh_81 = buffer.data(dsh + 81);
    const auto *dsh_82 = buffer.data(dsh + 82);
    const auto *dsh_83 = buffer.data(dsh + 83);
    const auto *dsh_96 = buffer.data(dsh + 96);
    const auto *dsh_99 = buffer.data(dsh + 99);
    const auto *dsh_100 = buffer.data(dsh + 100);
    const auto *dsh_101 = buffer.data(dsh + 101);
    const auto *dsh_102 = buffer.data(dsh + 102);
    const auto *dsh_103 = buffer.data(dsh + 103);

    const auto *dsi1_0 = buffer.data(dsi1 + 0);
    const auto *dsi1_3 = buffer.data(dsi1 + 3);
    const auto *dsi1_5 = buffer.data(dsi1 + 5);
    const auto *dsi1_6 = buffer.data(dsi1 + 6);
    const auto *dsi1_9 = buffer.data(dsi1 + 9);
    const auto *dsi1_10 = buffer.data(dsi1 + 10);
    const auto *dsi1_14 = buffer.data(dsi1 + 14);
    const auto *dsi1_21 = buffer.data(dsi1 + 21);
    const auto *dsi1_27 = buffer.data(dsi1 + 27);
    const auto *dsi1_31 = buffer.data(dsi1 + 31);
    const auto *dsi1_34 = buffer.data(dsi1 + 34);
    const auto *dsi1_38 = buffer.data(dsi1 + 38);
    const auto *dsi1_56 = buffer.data(dsi1 + 56);
    const auto *dsi1_61 = buffer.data(dsi1 + 61);
    const auto *dsi1_65 = buffer.data(dsi1 + 65);
    const auto *dsi1_70 = buffer.data(dsi1 + 70);
    const auto *dsi1_84 = buffer.data(dsi1 + 84);
    const auto *dsi1_87 = buffer.data(dsi1 + 87);
    const auto *dsi1_90 = buffer.data(dsi1 + 90);
    const auto *dsi1_94 = buffer.data(dsi1 + 94);
    const auto *dsi1_105 = buffer.data(dsi1 + 105);
    const auto *dsi1_107 = buffer.data(dsi1 + 107);
    const auto *dsi1_108 = buffer.data(dsi1 + 108);
    const auto *dsi1_109 = buffer.data(dsi1 + 109);
    const auto *dsi1_111 = buffer.data(dsi1 + 111);
    const auto *dsi1_124 = buffer.data(dsi1 + 124);

    const auto *fsg0_0 = buffer.data(fsg0 + 0);
    const auto *fsg0_1 = buffer.data(fsg0 + 1);
    const auto *fsg0_2 = buffer.data(fsg0 + 2);
    const auto *fsg0_3 = buffer.data(fsg0 + 3);
    const auto *fsg0_5 = buffer.data(fsg0 + 5);
    const auto *fsg0_10 = buffer.data(fsg0 + 10);
    const auto *fsg0_12 = buffer.data(fsg0 + 12);
    const auto *fsg0_13 = buffer.data(fsg0 + 13);
    const auto *fsg0_14 = buffer.data(fsg0 + 14);
    const auto *fsg0_18 = buffer.data(fsg0 + 18);
    const auto *fsg0_25 = buffer.data(fsg0 + 25);
    const auto *fsg0_26 = buffer.data(fsg0 + 26);
    const auto *fsg0_27 = buffer.data(fsg0 + 27);
    const auto *fsg0_32 = buffer.data(fsg0 + 32);
    const auto *fsg0_34 = buffer.data(fsg0 + 34);
    const auto *fsg0_35 = buffer.data(fsg0 + 35);
    const auto *fsg0_41 = buffer.data(fsg0 + 41);
    const auto *fsg0_42 = buffer.data(fsg0 + 42);
    const auto *fsg0_43 = buffer.data(fsg0 + 43);
    const auto *fsg0_44 = buffer.data(fsg0 + 44);
    const auto *fsg0_45 = buffer.data(fsg0 + 45);
    const auto *fsg0_47 = buffer.data(fsg0 + 47);
    const auto *fsg0_48 = buffer.data(fsg0 + 48);
    const auto *fsg0_50 = buffer.data(fsg0 + 50);

    const auto *fsg1_0 = buffer.data(fsg1 + 0);
    const auto *fsg1_1 = buffer.data(fsg1 + 1);
    const auto *fsg1_2 = buffer.data(fsg1 + 2);
    const auto *fsg1_3 = buffer.data(fsg1 + 3);
    const auto *fsg1_5 = buffer.data(fsg1 + 5);
    const auto *fsg1_10 = buffer.data(fsg1 + 10);
    const auto *fsg1_12 = buffer.data(fsg1 + 12);
    const auto *fsg1_13 = buffer.data(fsg1 + 13);
    const auto *fsg1_14 = buffer.data(fsg1 + 14);
    const auto *fsg1_18 = buffer.data(fsg1 + 18);
    const auto *fsg1_25 = buffer.data(fsg1 + 25);
    const auto *fsg1_26 = buffer.data(fsg1 + 26);
    const auto *fsg1_27 = buffer.data(fsg1 + 27);
    const auto *fsg1_32 = buffer.data(fsg1 + 32);
    const auto *fsg1_34 = buffer.data(fsg1 + 34);
    const auto *fsg1_35 = buffer.data(fsg1 + 35);
    const auto *fsg1_41 = buffer.data(fsg1 + 41);
    const auto *fsg1_42 = buffer.data(fsg1 + 42);
    const auto *fsg1_43 = buffer.data(fsg1 + 43);
    const auto *fsg1_44 = buffer.data(fsg1 + 44);
    const auto *fsg1_45 = buffer.data(fsg1 + 45);
    const auto *fsg1_47 = buffer.data(fsg1 + 47);
    const auto *fsg1_48 = buffer.data(fsg1 + 48);
    const auto *fsg1_50 = buffer.data(fsg1 + 50);

    const auto *fsh_0 = buffer.data(fsh + 0);
    const auto *fsh_1 = buffer.data(fsh + 1);
    const auto *fsh_2 = buffer.data(fsh + 2);
    const auto *fsh_3 = buffer.data(fsh + 3);
    const auto *fsh_5 = buffer.data(fsh + 5);
    const auto *fsh_6 = buffer.data(fsh + 6);
    const auto *fsh_8 = buffer.data(fsh + 8);
    const auto *fsh_9 = buffer.data(fsh + 9);
    const auto *fsh_10 = buffer.data(fsh + 10);
    const auto *fsh_14 = buffer.data(fsh + 14);
    const auto *fsh_15 = buffer.data(fsh + 15);
    const auto *fsh_17 = buffer.data(fsh + 17);
    const auto *fsh_18 = buffer.data(fsh + 18);
    const auto *fsh_19 = buffer.data(fsh + 19);
    const auto *fsh_20 = buffer.data(fsh + 20);
    const auto *fsh_21 = buffer.data(fsh + 21);
    const auto *fsh_22 = buffer.data(fsh + 22);
    const auto *fsh_24 = buffer.data(fsh + 24);
    const auto *fsh_26 = buffer.data(fsh + 26);
    const auto *fsh_27 = buffer.data(fsh + 27);
    const auto *fsh_28 = buffer.data(fsh + 28);
    const auto *fsh_30 = buffer.data(fsh + 30);
    const auto *fsh_31 = buffer.data(fsh + 31);
    const auto *fsh_36 = buffer.data(fsh + 36);
    const auto *fsh_37 = buffer.data(fsh + 37);
    const auto *fsh_38 = buffer.data(fsh + 38);
    const auto *fsh_39 = buffer.data(fsh + 39);
    const auto *fsh_40 = buffer.data(fsh + 40);
    const auto *fsh_41 = buffer.data(fsh + 41);
    const auto *fsh_42 = buffer.data(fsh + 42);
    const auto *fsh_44 = buffer.data(fsh + 44);
    const auto *fsh_46 = buffer.data(fsh + 46);
    const auto *fsh_47 = buffer.data(fsh + 47);
    const auto *fsh_49 = buffer.data(fsh + 49);
    const auto *fsh_50 = buffer.data(fsh + 50);
    const auto *fsh_51 = buffer.data(fsh + 51);
    const auto *fsh_56 = buffer.data(fsh + 56);
    const auto *fsh_57 = buffer.data(fsh + 57);
    const auto *fsh_58 = buffer.data(fsh + 58);
    const auto *fsh_59 = buffer.data(fsh + 59);
    const auto *fsh_60 = buffer.data(fsh + 60);
    const auto *fsh_61 = buffer.data(fsh + 61);
    const auto *fsh_62 = buffer.data(fsh + 62);
    const auto *fsh_63 = buffer.data(fsh + 63);
    const auto *fsh_64 = buffer.data(fsh + 64);
    const auto *fsh_65 = buffer.data(fsh + 65);
    const auto *fsh_66 = buffer.data(fsh + 66);
    const auto *fsh_68 = buffer.data(fsh + 68);
    const auto *fsh_69 = buffer.data(fsh + 69);
    const auto *fsh_70 = buffer.data(fsh + 70);
    const auto *fsh_72 = buffer.data(fsh + 72);
    const auto *fsh_73 = buffer.data(fsh + 73);
    const auto *fsh_78 = buffer.data(fsh + 78);
    const auto *fsh_80 = buffer.data(fsh + 80);
    const auto *fsh_81 = buffer.data(fsh + 81);
    const auto *fsh_82 = buffer.data(fsh + 82);
    const auto *fsh_83 = buffer.data(fsh + 83);
    const auto *fsh_84 = buffer.data(fsh + 84);
    const auto *fsh_86 = buffer.data(fsh + 86);
    const auto *fsh_87 = buffer.data(fsh + 87);
    const auto *fsh_89 = buffer.data(fsh + 89);
    const auto *fsh_90 = buffer.data(fsh + 90);
    const auto *fsh_93 = buffer.data(fsh + 93);
    const auto *fsh_99 = buffer.data(fsh + 99);
    const auto *fsh_100 = buffer.data(fsh + 100);
    const auto *fsh_101 = buffer.data(fsh + 101);
    const auto *fsh_102 = buffer.data(fsh + 102);
    const auto *fsh_103 = buffer.data(fsh + 103);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dsh_0, fsg0_0, \
                         fsg1_0, fsh_0, fsh_1, fsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsh_0[k]
                 + f_1 * fsg0_0[k]
                 - f_2 * fsg1_0[k]
                 + f_3 * pc_x[k] * fsh_0[k];

        t_1[k] = f_3 * pc_y[k] * fsh_0[k];

        t_2[k] = f_3 * pc_z[k] * fsh_0[k];

        t_3[k] = f_4 * fsg0_0[k]
                 - f_5 * fsg1_0[k]
                 + f_3 * pc_y[k] * fsh_1[k];

        t_4[k] = f_3 * pc_y[k] * fsh_2[k];

        t_5[k] = f_4 * fsg0_0[k]
                 - f_5 * fsg1_0[k]
                 + f_3 * pc_z[k] * fsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, fsg0_1, fsg0_2, fsg0_3, fsg1_1, \
                         fsg1_2, fsg1_3, fsh_3, fsh_5, fsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * fsg0_1[k]
                 - f_7 * fsg1_1[k]
                 + f_3 * pc_y[k] * fsh_3[k];

        t_7[k] = f_3 * pc_z[k] * fsh_3[k];

        t_8[k] = f_3 * pc_y[k] * fsh_5[k];

        t_9[k] = f_6 * fsg0_2[k]
                 - f_7 * fsg1_2[k]
                 + f_3 * pc_z[k] * fsh_5[k];

        t_10[k] = f_8 * fsg0_3[k]
                  - f_9 * fsg1_3[k]
                  + f_3 * pc_y[k] * fsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, dsh_15, fsg0_5, \
                         fsg1_5, fsh_6, fsh_8, fsh_9, fsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * fsh_6[k];

        t_12[k] = f_4 * fsg0_5[k]
                  - f_5 * fsg1_5[k]
                  + f_3 * pc_y[k] * fsh_8[k];

        t_13[k] = f_3 * pc_y[k] * fsh_9[k];

        t_14[k] = f_8 * fsg0_5[k]
                  - f_9 * fsg1_5[k]
                  + f_3 * pc_z[k] * fsh_9[k];

        t_15[k] = f_0 * dsh_15[k]
                  + f_3 * pc_x[k] * fsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, dsh_17, dsh_18, \
                         dsh_20, fsh_10, fsh_14, fsh_17, fsh_18, \
                         fsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * fsh_10[k];

        t_17[k] = f_0 * dsh_17[k]
                  + f_3 * pc_x[k] * fsh_17[k];

        t_18[k] = f_0 * dsh_18[k]
                  + f_3 * pc_x[k] * fsh_18[k];

        t_19[k] = f_3 * pc_y[k] * fsh_14[k];

        t_20[k] = f_0 * dsh_20[k]
                  + f_3 * pc_x[k] * fsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, fsg0_10, fsg0_12, fsg0_13, \
                         fsg1_10, fsg1_12, fsg1_13, fsh_15, fsh_17, \
                         fsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * fsg0_10[k]
                  - f_2 * fsg1_10[k]
                  + f_3 * pc_y[k] * fsh_15[k];

        t_22[k] = f_3 * pc_z[k] * fsh_15[k];

        t_23[k] = f_8 * fsg0_12[k]
                  - f_9 * fsg1_12[k]
                  + f_3 * pc_y[k] * fsh_17[k];

        t_24[k] = f_6 * fsg0_13[k]
                  - f_7 * fsg1_13[k]
                  + f_3 * pc_y[k] * fsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, dsi0_0, dsh_0, \
                         dsi1_0, fsg0_14, fsg1_14, fsh_19, fsh_20, \
                         fsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * fsg0_14[k]
                  - f_5 * fsg1_14[k]
                  + f_3 * pc_y[k] * fsh_19[k];

        t_26[k] = f_3 * pc_y[k] * fsh_20[k];

        t_27[k] = f_1 * fsg0_14[k]
                  - f_2 * fsg1_14[k]
                  + f_3 * pc_z[k] * fsh_20[k];

        t_28[k] = pa_y[k] * dsi0_0[k]
                  - f_10 * pc_y[k] * dsi1_0[k];

        t_29[k] = f_11 * dsh_0[k]
                  + f_3 * pc_y[k] * fsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, dsi0_3, dsi0_5, dsh_1, \
                         dsi1_3, dsi1_5, fsh_21, fsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * fsh_21[k];

        t_31[k] = pa_y[k] * dsi0_3[k]
                  + f_12 * dsh_1[k]
                  - f_10 * pc_y[k] * dsi1_3[k];

        t_32[k] = f_3 * pc_z[k] * fsh_22[k];

        t_33[k] = pa_y[k] * dsi0_5[k]
                  - f_10 * pc_y[k] * dsi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, dsi0_6, dsi0_9, dsh_3, \
                         dsh_5, dsi1_6, dsi1_9, fsh_24, fsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * dsi0_6[k]
                  + f_0 * dsh_3[k]
                  - f_10 * pc_y[k] * dsi1_6[k];

        t_35[k] = f_3 * pc_z[k] * fsh_24[k];

        t_36[k] = f_11 * dsh_5[k]
                  + f_3 * pc_y[k] * fsh_26[k];

        t_37[k] = pa_y[k] * dsi0_9[k]
                  - f_10 * pc_y[k] * dsi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, dsi0_10, dsh_6, dsh_9, \
                         dsi1_10, fsg0_18, fsg1_18, fsh_27, fsh_28, \
                         fsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * dsi0_10[k]
                  + f_13 * dsh_6[k]
                  - f_10 * pc_y[k] * dsi1_10[k];

        t_39[k] = f_3 * pc_z[k] * fsh_27[k];

        t_40[k] = f_4 * fsg0_18[k]
                  - f_5 * fsg1_18[k]
                  + f_3 * pc_z[k] * fsh_28[k];

        t_41[k] = f_11 * dsh_9[k]
                  + f_3 * pc_y[k] * fsh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, dsi0_14, dsh_36, \
                         dsh_38, dsi1_14, fsh_31, fsh_36, fsh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * dsi0_14[k]
                  - f_10 * pc_y[k] * dsi1_14[k];

        t_43[k] = f_12 * dsh_36[k]
                  + f_3 * pc_x[k] * fsh_36[k];

        t_44[k] = f_3 * pc_z[k] * fsh_31[k];

        t_45[k] = f_12 * dsh_38[k]
                  + f_3 * pc_x[k] * fsh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, dsh_15, dsh_39, dsh_40, dsh_41, \
                         fsg0_25, fsg1_25, fsh_36, fsh_39, fsh_40, \
                         fsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * dsh_39[k]
                  + f_3 * pc_x[k] * fsh_39[k];

        t_47[k] = f_12 * dsh_40[k]
                  + f_3 * pc_x[k] * fsh_40[k];

        t_48[k] = f_12 * dsh_41[k]
                  + f_3 * pc_x[k] * fsh_41[k];

        t_49[k] = f_11 * dsh_15[k]
                  + f_1 * fsg0_25[k]
                  - f_2 * fsg1_25[k]
                  + f_3 * pc_y[k] * fsh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, fsg0_25, fsg0_26, fsg0_27, fsg1_25, \
                         fsg1_26, fsg1_27, fsh_36, fsh_37, fsh_38, \
                         fsh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * fsh_36[k];

        t_51[k] = f_4 * fsg0_25[k]
                  - f_5 * fsg1_25[k]
                  + f_3 * pc_z[k] * fsh_37[k];

        t_52[k] = f_6 * fsg0_26[k]
                  - f_7 * fsg1_26[k]
                  + f_3 * pc_z[k] * fsh_38[k];

        t_53[k] = f_8 * fsg0_27[k]
                  - f_9 * fsg1_27[k]
                  + f_3 * pc_z[k] * fsh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, dsi0_0, dsi0_27, \
                         dsh_20, dsi1_0, dsi1_27, fsh_41, fsh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * dsh_20[k]
                  + f_3 * pc_y[k] * fsh_41[k];

        t_55[k] = pa_y[k] * dsi0_27[k]
                  - f_10 * pc_y[k] * dsi1_27[k];

        t_56[k] = pa_z[k] * dsi0_0[k]
                  - f_10 * pc_z[k] * dsi1_0[k];

        t_57[k] = f_3 * pc_y[k] * fsh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, dsi0_3, dsi0_5, dsh_0, \
                         dsh_2, dsi1_3, dsi1_5, fsh_42, fsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * dsh_0[k]
                  + f_3 * pc_z[k] * fsh_42[k];

        t_59[k] = pa_z[k] * dsi0_3[k]
                  - f_10 * pc_z[k] * dsi1_3[k];

        t_60[k] = f_3 * pc_y[k] * fsh_44[k];

        t_61[k] = pa_z[k] * dsi0_5[k]
                  + f_12 * dsh_2[k]
                  - f_10 * pc_z[k] * dsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, dsi0_6, dsi0_9, dsh_5, \
                         dsi1_6, dsi1_9, fsg0_32, fsg1_32, fsh_46, \
                         fsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * dsi0_6[k]
                  - f_10 * pc_z[k] * dsi1_6[k];

        t_63[k] = f_4 * fsg0_32[k]
                  - f_5 * fsg1_32[k]
                  + f_3 * pc_y[k] * fsh_46[k];

        t_64[k] = f_3 * pc_y[k] * fsh_47[k];

        t_65[k] = pa_z[k] * dsi0_9[k]
                  + f_0 * dsh_5[k]
                  - f_10 * pc_z[k] * dsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, dsi0_10, dsi1_10, fsg0_34, \
                         fsg0_35, fsg1_34, fsg1_35, fsh_49, fsh_50, \
                         fsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * dsi0_10[k]
                  - f_10 * pc_z[k] * dsi1_10[k];

        t_67[k] = f_6 * fsg0_34[k]
                  - f_7 * fsg1_34[k]
                  + f_3 * pc_y[k] * fsh_49[k];

        t_68[k] = f_4 * fsg0_35[k]
                  - f_5 * fsg1_35[k]
                  + f_3 * pc_y[k] * fsh_50[k];

        t_69[k] = f_3 * pc_y[k] * fsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, dsi0_14, dsh_9, dsh_57, \
                         dsh_58, dsh_59, dsi1_14, fsh_57, fsh_58, \
                         fsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * dsi0_14[k]
                  + f_13 * dsh_9[k]
                  - f_10 * pc_z[k] * dsi1_14[k];

        t_71[k] = f_12 * dsh_57[k]
                  + f_3 * pc_x[k] * fsh_57[k];

        t_72[k] = f_12 * dsh_58[k]
                  + f_3 * pc_x[k] * fsh_58[k];

        t_73[k] = f_12 * dsh_59[k]
                  + f_3 * pc_x[k] * fsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, dsi0_21, dsh_60, \
                         dsh_62, dsi1_21, fsh_56, fsh_60, fsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_12 * dsh_60[k]
                  + f_3 * pc_x[k] * fsh_60[k];

        t_75[k] = f_3 * pc_y[k] * fsh_56[k];

        t_76[k] = f_12 * dsh_62[k]
                  + f_3 * pc_x[k] * fsh_62[k];

        t_77[k] = pa_z[k] * dsi0_21[k]
                  - f_10 * pc_z[k] * dsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, fsg0_41, fsg0_42, fsg0_43, fsg1_41, fsg1_42, \
                         fsg1_43, fsh_58, fsh_59, fsh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_14 * fsg0_41[k]
                  - f_15 * fsg1_41[k]
                  + f_3 * pc_y[k] * fsh_58[k];

        t_79[k] = f_8 * fsg0_42[k]
                  - f_9 * fsg1_42[k]
                  + f_3 * pc_y[k] * fsh_59[k];

        t_80[k] = f_6 * fsg0_43[k]
                  - f_7 * fsg1_43[k]
                  + f_3 * pc_y[k] * fsh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_x, pc_x, pc_y, pc_z, dsi0_84, dsh_20, \
                         dsh_63, dsi1_84, fsg0_44, fsg1_44, fsh_61, \
                         fsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * fsg0_44[k]
                  - f_5 * fsg1_44[k]
                  + f_3 * pc_y[k] * fsh_61[k];

        t_82[k] = f_3 * pc_y[k] * fsh_62[k];

        t_83[k] = f_11 * dsh_20[k]
                  + f_1 * fsg0_44[k]
                  - f_2 * fsg1_44[k]
                  + f_3 * pc_z[k] * fsh_62[k];

        t_84[k] = pa_x[k] * dsi0_84[k]
                  + f_16 * dsh_63[k]
                  - f_10 * pc_x[k] * dsi1_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_x, pc_x, pc_y, pc_z, dsi0_87, dsh_21, \
                         dsh_66, dsi1_87, fsh_63, fsh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * dsh_21[k]
                  + f_3 * pc_y[k] * fsh_63[k];

        t_86[k] = f_3 * pc_z[k] * fsh_63[k];

        t_87[k] = pa_x[k] * dsi0_87[k]
                  + f_13 * dsh_66[k]
                  - f_10 * pc_x[k] * dsi1_87[k];

        t_88[k] = f_3 * pc_z[k] * fsh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pc_x, pc_z, dsi0_90, dsh_69, dsi1_90, \
                         fsg0_45, fsg1_45, fsh_65, fsh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * fsg0_45[k]
                  - f_5 * fsg1_45[k]
                  + f_3 * pc_z[k] * fsh_65[k];

        t_90[k] = pa_x[k] * dsi0_90[k]
                  + f_0 * dsh_69[k]
                  - f_10 * pc_x[k] * dsi1_90[k];

        t_91[k] = f_3 * pc_z[k] * fsh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pc_x, pc_y, pc_z, dsi0_94, dsh_26, \
                         dsh_73, dsi1_94, fsg0_47, fsg1_47, fsh_68, \
                         fsh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * dsh_26[k]
                  + f_3 * pc_y[k] * fsh_68[k];

        t_93[k] = f_6 * fsg0_47[k]
                  - f_7 * fsg1_47[k]
                  + f_3 * pc_z[k] * fsh_68[k];

        t_94[k] = pa_x[k] * dsi0_94[k]
                  + f_12 * dsh_73[k]
                  - f_10 * pc_x[k] * dsi1_94[k];

        t_95[k] = f_3 * pc_z[k] * fsh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, dsh_30, dsh_78, fsg0_48, \
                         fsg0_50, fsg1_48, fsg1_50, fsh_70, fsh_72, \
                         fsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * fsg0_48[k]
                  - f_5 * fsg1_48[k]
                  + f_3 * pc_z[k] * fsh_70[k];

        t_97[k] = f_12 * dsh_30[k]
                  + f_3 * pc_y[k] * fsh_72[k];

        t_98[k] = f_8 * fsg0_50[k]
                  - f_9 * fsg1_50[k]
                  + f_3 * pc_z[k] * fsh_72[k];

        t_99[k] = f_11 * dsh_78[k]
                  + f_3 * pc_x[k] * fsh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, dsh_80, dsh_81, \
                         dsh_82, dsh_83, fsh_73, fsh_80, fsh_81, fsh_82, \
                         fsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * fsh_73[k];

        t_101[k] = f_11 * dsh_80[k]
                   + f_3 * pc_x[k] * fsh_80[k];

        t_102[k] = f_11 * dsh_81[k]
                   + f_3 * pc_x[k] * fsh_81[k];

        t_103[k] = f_11 * dsh_82[k]
                   + f_3 * pc_x[k] * fsh_82[k];

        t_104[k] = f_11 * dsh_83[k]
                   + f_3 * pc_x[k] * fsh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pc_x, pc_z, dsi0_105, dsi0_107, \
                         dsi0_108, dsi1_105, dsi1_107, dsi1_108, \
                         fsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_x[k] * dsi0_105[k]
                   - f_10 * pc_x[k] * dsi1_105[k];

        t_106[k] = f_3 * pc_z[k] * fsh_78[k];

        t_107[k] = pa_x[k] * dsi0_107[k]
                   - f_10 * pc_x[k] * dsi1_107[k];

        t_108[k] = pa_x[k] * dsi0_108[k]
                   - f_10 * pc_x[k] * dsi1_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_y, pc_x, pc_y, dsi0_56, \
                         dsi0_109, dsi0_111, dsh_41, dsi1_56, dsi1_109, dsi1_111, \
                         fsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * dsi0_109[k]
                   - f_10 * pc_x[k] * dsi1_109[k];

        t_110[k] = f_12 * dsh_41[k]
                   + f_3 * pc_y[k] * fsh_83[k];

        t_111[k] = pa_x[k] * dsi0_111[k]
                   - f_10 * pc_x[k] * dsi1_111[k];

        t_112[k] = pa_y[k] * dsi0_56[k]
                   - f_10 * pc_y[k] * dsi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, dsi0_31, dsh_21, \
                         dsh_42, dsh_44, dsi1_31, fsh_84, fsh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * dsh_42[k]
                   + f_3 * pc_y[k] * fsh_84[k];

        t_114[k] = f_11 * dsh_21[k]
                   + f_3 * pc_z[k] * fsh_84[k];

        t_115[k] = pa_z[k] * dsi0_31[k]
                   - f_10 * pc_z[k] * dsi1_31[k];

        t_116[k] = f_11 * dsh_44[k]
                   + f_3 * pc_y[k] * fsh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, dsi0_34, dsi0_61, \
                         dsh_24, dsh_47, dsi1_34, dsi1_61, fsh_87, \
                         fsh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * dsi0_61[k]
                   - f_10 * pc_y[k] * dsi1_61[k];

        t_118[k] = pa_z[k] * dsi0_34[k]
                   - f_10 * pc_z[k] * dsi1_34[k];

        t_119[k] = f_11 * dsh_24[k]
                   + f_3 * pc_z[k] * fsh_87[k];

        t_120[k] = f_11 * dsh_47[k]
                   + f_3 * pc_y[k] * fsh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, dsi0_38, dsi0_65, \
                         dsh_27, dsi1_38, dsi1_65, fsh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * dsi0_65[k]
                   - f_10 * pc_y[k] * dsi1_65[k];

        t_122[k] = pa_z[k] * dsi0_38[k]
                   - f_10 * pc_z[k] * dsi1_38[k];

        t_123[k] = f_11 * dsh_27[k]
                   + f_3 * pc_z[k] * fsh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pa_y, pc_x, pc_y, dsi0_70, dsi0_124, \
                         dsh_51, dsh_96, dsi1_70, dsi1_124, fsh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_x[k] * dsi0_124[k]
                   + f_12 * dsh_96[k]
                   - f_10 * pc_x[k] * dsi1_124[k];

        t_125[k] = f_11 * dsh_51[k]
                   + f_3 * pc_y[k] * fsh_93[k];

        t_126[k] = pa_y[k] * dsi0_70[k]
                   - f_10 * pc_y[k] * dsi1_70[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, dsh_99, dsh_100, dsh_101, \
                         dsh_102, dsh_103, fsh_99, fsh_100, fsh_101, fsh_102, \
                         fsh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_11 * dsh_99[k]
                   + f_3 * pc_x[k] * fsh_99[k];

        t_128[k] = f_11 * dsh_100[k]
                   + f_3 * pc_x[k] * fsh_100[k];

        t_129[k] = f_11 * dsh_101[k]
                   + f_3 * pc_x[k] * fsh_101[k];

        t_130[k] = f_11 * dsh_102[k]
                   + f_3 * pc_x[k] * fsh_102[k];

        t_131[k] = f_11 * dsh_103[k]
                   + f_3 * pc_x[k] * fsh_103[k];
    }
}

static auto
compute_prim_fsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsi0,
                                                          const size_t dsh, const size_t dsi1,
                                                          const size_t fsg0, const size_t fsg1,
                                                          const size_t fsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);
    const auto f_16 = 3.0 / q;

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

    const auto *dsi0_84 = buffer.data(dsi0 + 84);
    const auto *dsi0_85 = buffer.data(dsi0 + 85);
    const auto *dsi0_87 = buffer.data(dsi0 + 87);
    const auto *dsi0_90 = buffer.data(dsi0 + 90);
    const auto *dsi0_94 = buffer.data(dsi0 + 94);
    const auto *dsi0_105 = buffer.data(dsi0 + 105);
    const auto *dsi0_107 = buffer.data(dsi0 + 107);
    const auto *dsi0_108 = buffer.data(dsi0 + 108);
    const auto *dsi0_109 = buffer.data(dsi0 + 109);
    const auto *dsi0_133 = buffer.data(dsi0 + 133);
    const auto *dsi0_135 = buffer.data(dsi0 + 135);
    const auto *dsi0_136 = buffer.data(dsi0 + 136);
    const auto *dsi0_137 = buffer.data(dsi0 + 137);
    const auto *dsi0_139 = buffer.data(dsi0 + 139);
    const auto *dsi0_140 = buffer.data(dsi0 + 140);
    const auto *dsi0_142 = buffer.data(dsi0 + 142);
    const auto *dsi0_145 = buffer.data(dsi0 + 145);
    const auto *dsi0_149 = buffer.data(dsi0 + 149);
    const auto *dsi0_154 = buffer.data(dsi0 + 154);
    const auto *dsi0_161 = buffer.data(dsi0 + 161);
    const auto *dsi0_162 = buffer.data(dsi0 + 162);
    const auto *dsi0_163 = buffer.data(dsi0 + 163);
    const auto *dsi0_164 = buffer.data(dsi0 + 164);
    const auto *dsi0_165 = buffer.data(dsi0 + 165);
    const auto *dsi0_167 = buffer.data(dsi0 + 167);

    const auto *dsh_36 = buffer.data(dsh + 36);
    const auto *dsh_42 = buffer.data(dsh + 42);
    const auto *dsh_62 = buffer.data(dsh + 62);
    const auto *dsh_78 = buffer.data(dsh + 78);
    const auto *dsh_79 = buffer.data(dsh + 79);
    const auto *dsh_80 = buffer.data(dsh + 80);
    const auto *dsh_81 = buffer.data(dsh + 81);
    const auto *dsh_83 = buffer.data(dsh + 83);
    const auto *dsh_99 = buffer.data(dsh + 99);
    const auto *dsh_104 = buffer.data(dsh + 104);
    const auto *dsh_105 = buffer.data(dsh + 105);
    const auto *dsh_110 = buffer.data(dsh + 110);
    const auto *dsh_114 = buffer.data(dsh + 114);
    const auto *dsh_119 = buffer.data(dsh + 119);
    const auto *dsh_120 = buffer.data(dsh + 120);
    const auto *dsh_121 = buffer.data(dsh + 121);
    const auto *dsh_122 = buffer.data(dsh + 122);
    const auto *dsh_123 = buffer.data(dsh + 123);
    const auto *dsh_124 = buffer.data(dsh + 124);
    const auto *dsh_125 = buffer.data(dsh + 125);

    const auto *dsi1_84 = buffer.data(dsi1 + 84);
    const auto *dsi1_85 = buffer.data(dsi1 + 85);
    const auto *dsi1_87 = buffer.data(dsi1 + 87);
    const auto *dsi1_90 = buffer.data(dsi1 + 90);
    const auto *dsi1_94 = buffer.data(dsi1 + 94);
    const auto *dsi1_105 = buffer.data(dsi1 + 105);
    const auto *dsi1_107 = buffer.data(dsi1 + 107);
    const auto *dsi1_108 = buffer.data(dsi1 + 108);
    const auto *dsi1_109 = buffer.data(dsi1 + 109);
    const auto *dsi1_133 = buffer.data(dsi1 + 133);
    const auto *dsi1_135 = buffer.data(dsi1 + 135);
    const auto *dsi1_136 = buffer.data(dsi1 + 136);
    const auto *dsi1_137 = buffer.data(dsi1 + 137);
    const auto *dsi1_139 = buffer.data(dsi1 + 139);
    const auto *dsi1_140 = buffer.data(dsi1 + 140);
    const auto *dsi1_142 = buffer.data(dsi1 + 142);
    const auto *dsi1_145 = buffer.data(dsi1 + 145);
    const auto *dsi1_149 = buffer.data(dsi1 + 149);
    const auto *dsi1_154 = buffer.data(dsi1 + 154);
    const auto *dsi1_161 = buffer.data(dsi1 + 161);
    const auto *dsi1_162 = buffer.data(dsi1 + 162);
    const auto *dsi1_163 = buffer.data(dsi1 + 163);
    const auto *dsi1_164 = buffer.data(dsi1 + 164);
    const auto *dsi1_165 = buffer.data(dsi1 + 165);
    const auto *dsi1_167 = buffer.data(dsi1 + 167);

    const auto *fsg0_75 = buffer.data(fsg0 + 75);
    const auto *fsg0_76 = buffer.data(fsg0 + 76);
    const auto *fsg0_77 = buffer.data(fsg0 + 77);
    const auto *fsg0_78 = buffer.data(fsg0 + 78);
    const auto *fsg0_79 = buffer.data(fsg0 + 79);
    const auto *fsg0_80 = buffer.data(fsg0 + 80);
    const auto *fsg0_90 = buffer.data(fsg0 + 90);
    const auto *fsg0_91 = buffer.data(fsg0 + 91);
    const auto *fsg0_93 = buffer.data(fsg0 + 93);
    const auto *fsg0_95 = buffer.data(fsg0 + 95);
    const auto *fsg0_96 = buffer.data(fsg0 + 96);
    const auto *fsg0_98 = buffer.data(fsg0 + 98);
    const auto *fsg0_99 = buffer.data(fsg0 + 99);
    const auto *fsg0_100 = buffer.data(fsg0 + 100);
    const auto *fsg0_101 = buffer.data(fsg0 + 101);
    const auto *fsg0_102 = buffer.data(fsg0 + 102);
    const auto *fsg0_103 = buffer.data(fsg0 + 103);
    const auto *fsg0_104 = buffer.data(fsg0 + 104);
    const auto *fsg0_107 = buffer.data(fsg0 + 107);
    const auto *fsg0_109 = buffer.data(fsg0 + 109);
    const auto *fsg0_110 = buffer.data(fsg0 + 110);
    const auto *fsg0_112 = buffer.data(fsg0 + 112);
    const auto *fsg0_113 = buffer.data(fsg0 + 113);
    const auto *fsg0_114 = buffer.data(fsg0 + 114);
    const auto *fsg0_116 = buffer.data(fsg0 + 116);
    const auto *fsg0_117 = buffer.data(fsg0 + 117);
    const auto *fsg0_118 = buffer.data(fsg0 + 118);
    const auto *fsg0_119 = buffer.data(fsg0 + 119);
    const auto *fsg0_121 = buffer.data(fsg0 + 121);
    const auto *fsg0_123 = buffer.data(fsg0 + 123);
    const auto *fsg0_124 = buffer.data(fsg0 + 124);
    const auto *fsg0_126 = buffer.data(fsg0 + 126);
    const auto *fsg0_127 = buffer.data(fsg0 + 127);
    const auto *fsg0_128 = buffer.data(fsg0 + 128);
    const auto *fsg0_130 = buffer.data(fsg0 + 130);
    const auto *fsg0_131 = buffer.data(fsg0 + 131);
    const auto *fsg0_132 = buffer.data(fsg0 + 132);
    const auto *fsg0_133 = buffer.data(fsg0 + 133);
    const auto *fsg0_135 = buffer.data(fsg0 + 135);
    const auto *fsg0_137 = buffer.data(fsg0 + 137);
    const auto *fsg0_138 = buffer.data(fsg0 + 138);

    const auto *fsg1_75 = buffer.data(fsg1 + 75);
    const auto *fsg1_76 = buffer.data(fsg1 + 76);
    const auto *fsg1_77 = buffer.data(fsg1 + 77);
    const auto *fsg1_78 = buffer.data(fsg1 + 78);
    const auto *fsg1_79 = buffer.data(fsg1 + 79);
    const auto *fsg1_80 = buffer.data(fsg1 + 80);
    const auto *fsg1_90 = buffer.data(fsg1 + 90);
    const auto *fsg1_91 = buffer.data(fsg1 + 91);
    const auto *fsg1_93 = buffer.data(fsg1 + 93);
    const auto *fsg1_95 = buffer.data(fsg1 + 95);
    const auto *fsg1_96 = buffer.data(fsg1 + 96);
    const auto *fsg1_98 = buffer.data(fsg1 + 98);
    const auto *fsg1_99 = buffer.data(fsg1 + 99);
    const auto *fsg1_100 = buffer.data(fsg1 + 100);
    const auto *fsg1_101 = buffer.data(fsg1 + 101);
    const auto *fsg1_102 = buffer.data(fsg1 + 102);
    const auto *fsg1_103 = buffer.data(fsg1 + 103);
    const auto *fsg1_104 = buffer.data(fsg1 + 104);
    const auto *fsg1_107 = buffer.data(fsg1 + 107);
    const auto *fsg1_109 = buffer.data(fsg1 + 109);
    const auto *fsg1_110 = buffer.data(fsg1 + 110);
    const auto *fsg1_112 = buffer.data(fsg1 + 112);
    const auto *fsg1_113 = buffer.data(fsg1 + 113);
    const auto *fsg1_114 = buffer.data(fsg1 + 114);
    const auto *fsg1_116 = buffer.data(fsg1 + 116);
    const auto *fsg1_117 = buffer.data(fsg1 + 117);
    const auto *fsg1_118 = buffer.data(fsg1 + 118);
    const auto *fsg1_119 = buffer.data(fsg1 + 119);
    const auto *fsg1_121 = buffer.data(fsg1 + 121);
    const auto *fsg1_123 = buffer.data(fsg1 + 123);
    const auto *fsg1_124 = buffer.data(fsg1 + 124);
    const auto *fsg1_126 = buffer.data(fsg1 + 126);
    const auto *fsg1_127 = buffer.data(fsg1 + 127);
    const auto *fsg1_128 = buffer.data(fsg1 + 128);
    const auto *fsg1_130 = buffer.data(fsg1 + 130);
    const auto *fsg1_131 = buffer.data(fsg1 + 131);
    const auto *fsg1_132 = buffer.data(fsg1 + 132);
    const auto *fsg1_133 = buffer.data(fsg1 + 133);
    const auto *fsg1_135 = buffer.data(fsg1 + 135);
    const auto *fsg1_137 = buffer.data(fsg1 + 137);
    const auto *fsg1_138 = buffer.data(fsg1 + 138);

    const auto *fsh_99 = buffer.data(fsh + 99);
    const auto *fsh_104 = buffer.data(fsh + 104);
    const auto *fsh_105 = buffer.data(fsh + 105);
    const auto *fsh_106 = buffer.data(fsh + 106);
    const auto *fsh_107 = buffer.data(fsh + 107);
    const auto *fsh_108 = buffer.data(fsh + 108);
    const auto *fsh_109 = buffer.data(fsh + 109);
    const auto *fsh_110 = buffer.data(fsh + 110);
    const auto *fsh_111 = buffer.data(fsh + 111);
    const auto *fsh_112 = buffer.data(fsh + 112);
    const auto *fsh_113 = buffer.data(fsh + 113);
    const auto *fsh_114 = buffer.data(fsh + 114);
    const auto *fsh_119 = buffer.data(fsh + 119);
    const auto *fsh_120 = buffer.data(fsh + 120);
    const auto *fsh_121 = buffer.data(fsh + 121);
    const auto *fsh_122 = buffer.data(fsh + 122);
    const auto *fsh_123 = buffer.data(fsh + 123);
    const auto *fsh_125 = buffer.data(fsh + 125);
    const auto *fsh_126 = buffer.data(fsh + 126);
    const auto *fsh_127 = buffer.data(fsh + 127);
    const auto *fsh_129 = buffer.data(fsh + 129);
    const auto *fsh_131 = buffer.data(fsh + 131);
    const auto *fsh_132 = buffer.data(fsh + 132);
    const auto *fsh_134 = buffer.data(fsh + 134);
    const auto *fsh_135 = buffer.data(fsh + 135);
    const auto *fsh_136 = buffer.data(fsh + 136);
    const auto *fsh_138 = buffer.data(fsh + 138);
    const auto *fsh_139 = buffer.data(fsh + 139);
    const auto *fsh_140 = buffer.data(fsh + 140);
    const auto *fsh_141 = buffer.data(fsh + 141);
    const auto *fsh_142 = buffer.data(fsh + 142);
    const auto *fsh_143 = buffer.data(fsh + 143);
    const auto *fsh_144 = buffer.data(fsh + 144);
    const auto *fsh_145 = buffer.data(fsh + 145);
    const auto *fsh_146 = buffer.data(fsh + 146);
    const auto *fsh_149 = buffer.data(fsh + 149);
    const auto *fsh_151 = buffer.data(fsh + 151);
    const auto *fsh_152 = buffer.data(fsh + 152);
    const auto *fsh_154 = buffer.data(fsh + 154);
    const auto *fsh_155 = buffer.data(fsh + 155);
    const auto *fsh_156 = buffer.data(fsh + 156);
    const auto *fsh_158 = buffer.data(fsh + 158);
    const auto *fsh_159 = buffer.data(fsh + 159);
    const auto *fsh_160 = buffer.data(fsh + 160);
    const auto *fsh_161 = buffer.data(fsh + 161);
    const auto *fsh_162 = buffer.data(fsh + 162);
    const auto *fsh_163 = buffer.data(fsh + 163);
    const auto *fsh_164 = buffer.data(fsh + 164);
    const auto *fsh_165 = buffer.data(fsh + 165);
    const auto *fsh_166 = buffer.data(fsh + 166);
    const auto *fsh_167 = buffer.data(fsh + 167);
    const auto *fsh_169 = buffer.data(fsh + 169);
    const auto *fsh_171 = buffer.data(fsh + 171);
    const auto *fsh_172 = buffer.data(fsh + 172);
    const auto *fsh_174 = buffer.data(fsh + 174);
    const auto *fsh_175 = buffer.data(fsh + 175);
    const auto *fsh_176 = buffer.data(fsh + 176);
    const auto *fsh_178 = buffer.data(fsh + 178);
    const auto *fsh_179 = buffer.data(fsh + 179);
    const auto *fsh_180 = buffer.data(fsh + 180);
    const auto *fsh_181 = buffer.data(fsh + 181);
    const auto *fsh_183 = buffer.data(fsh + 183);
    const auto *fsh_184 = buffer.data(fsh + 184);
    const auto *fsh_185 = buffer.data(fsh + 185);
    const auto *fsh_186 = buffer.data(fsh + 186);
    const auto *fsh_187 = buffer.data(fsh + 187);
    const auto *fsh_188 = buffer.data(fsh + 188);
    const auto *fsh_189 = buffer.data(fsh + 189);
    const auto *fsh_191 = buffer.data(fsh + 191);
    const auto *fsh_192 = buffer.data(fsh + 192);

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pc_x, pc_z, dsi0_133, dsi0_135, \
                         dsh_36, dsh_104, dsi1_133, dsi1_135, fsh_99, \
                         fsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_11 * dsh_104[k]
                   + f_3 * pc_x[k] * fsh_104[k];

        t_133[k] = pa_x[k] * dsi0_133[k]
                   - f_10 * pc_x[k] * dsi1_133[k];

        t_134[k] = f_11 * dsh_36[k]
                   + f_3 * pc_z[k] * fsh_99[k];

        t_135[k] = pa_x[k] * dsi0_135[k]
                   - f_10 * pc_x[k] * dsi1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, pc_x, pc_y, dsi0_136, dsi0_137, \
                         dsi0_139, dsh_62, dsi1_136, dsi1_137, dsi1_139, \
                         fsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * dsi0_136[k]
                   - f_10 * pc_x[k] * dsi1_136[k];

        t_137[k] = pa_x[k] * dsi0_137[k]
                   - f_10 * pc_x[k] * dsi1_137[k];

        t_138[k] = f_11 * dsh_62[k]
                   + f_3 * pc_y[k] * fsh_104[k];

        t_139[k] = pa_x[k] * dsi0_139[k]
                   - f_10 * pc_x[k] * dsi1_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pc_x, pc_y, pc_z, dsi0_140, dsh_42, \
                         dsh_105, dsi1_140, fsg0_75, fsg1_75, fsh_105, \
                         fsh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pa_x[k] * dsi0_140[k]
                   + f_16 * dsh_105[k]
                   - f_10 * pc_x[k] * dsi1_140[k];

        t_141[k] = f_3 * pc_y[k] * fsh_105[k];

        t_142[k] = f_12 * dsh_42[k]
                   + f_3 * pc_z[k] * fsh_105[k];

        t_143[k] = f_4 * fsg0_75[k]
                   - f_5 * fsg1_75[k]
                   + f_3 * pc_y[k] * fsh_106[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pa_x, pc_x, pc_y, dsi0_145, dsh_110, dsi1_145, \
                         fsg0_76, fsg1_76, fsh_107, fsh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_3 * pc_y[k] * fsh_107[k];

        t_145[k] = pa_x[k] * dsi0_145[k]
                   + f_13 * dsh_110[k]
                   - f_10 * pc_x[k] * dsi1_145[k];

        t_146[k] = f_6 * fsg0_76[k]
                   - f_7 * fsg1_76[k]
                   + f_3 * pc_y[k] * fsh_108[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_x, pc_x, pc_y, dsi0_149, dsh_114, dsi1_149, \
                         fsg0_77, fsg1_77, fsh_109, fsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * fsg0_77[k]
                   - f_5 * fsg1_77[k]
                   + f_3 * pc_y[k] * fsh_109[k];

        t_148[k] = f_3 * pc_y[k] * fsh_110[k];

        t_149[k] = pa_x[k] * dsi0_149[k]
                   + f_0 * dsh_114[k]
                   - f_10 * pc_x[k] * dsi1_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_y, fsg0_78, fsg0_79, fsg0_80, fsg1_78, \
                         fsg1_79, fsg1_80, fsh_111, fsh_112, fsh_113, \
                         fsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_8 * fsg0_78[k]
                   - f_9 * fsg1_78[k]
                   + f_3 * pc_y[k] * fsh_111[k];

        t_151[k] = f_6 * fsg0_79[k]
                   - f_7 * fsg1_79[k]
                   + f_3 * pc_y[k] * fsh_112[k];

        t_152[k] = f_4 * fsg0_80[k]
                   - f_5 * fsg1_80[k]
                   + f_3 * pc_y[k] * fsh_113[k];

        t_153[k] = f_3 * pc_y[k] * fsh_114[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pc_x, dsi0_154, dsh_119, dsh_120, \
                         dsh_121, dsh_122, dsi1_154, fsh_120, fsh_121, \
                         fsh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_x[k] * dsi0_154[k]
                   + f_12 * dsh_119[k]
                   - f_10 * pc_x[k] * dsi1_154[k];

        t_155[k] = f_11 * dsh_120[k]
                   + f_3 * pc_x[k] * fsh_120[k];

        t_156[k] = f_11 * dsh_121[k]
                   + f_3 * pc_x[k] * fsh_121[k];

        t_157[k] = f_11 * dsh_122[k]
                   + f_3 * pc_x[k] * fsh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_x, pc_x, pc_y, dsi0_161, dsh_123, \
                         dsh_125, dsi1_161, fsh_119, fsh_123, fsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_11 * dsh_123[k]
                   + f_3 * pc_x[k] * fsh_123[k];

        t_159[k] = f_3 * pc_y[k] * fsh_119[k];

        t_160[k] = f_11 * dsh_125[k]
                   + f_3 * pc_x[k] * fsh_125[k];

        t_161[k] = pa_x[k] * dsi0_161[k]
                   - f_10 * pc_x[k] * dsi1_161[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_x, pc_x, dsi0_162, dsi0_163, dsi0_164, \
                         dsi0_165, dsi1_162, dsi1_163, dsi1_164, \
                         dsi1_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_x[k] * dsi0_162[k]
                   - f_10 * pc_x[k] * dsi1_162[k];

        t_163[k] = pa_x[k] * dsi0_163[k]
                   - f_10 * pc_x[k] * dsi1_163[k];

        t_164[k] = pa_x[k] * dsi0_164[k]
                   - f_10 * pc_x[k] * dsi1_164[k];

        t_165[k] = pa_x[k] * dsi0_165[k]
                   - f_10 * pc_x[k] * dsi1_165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_x, pc_x, pc_y, dsi0_167, dsi1_167, \
                         fsg0_90, fsg0_91, fsg1_90, fsg1_91, fsh_125, fsh_126, \
                         fsh_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_3 * pc_y[k] * fsh_125[k];

        t_167[k] = pa_x[k] * dsi0_167[k]
                   - f_10 * pc_x[k] * dsi1_167[k];

        t_168[k] = f_1 * fsg0_90[k]
                   - f_2 * fsg1_90[k]
                   + f_3 * pc_x[k] * fsh_126[k];

        t_169[k] = f_14 * fsg0_91[k]
                   - f_15 * fsg1_91[k]
                   + f_3 * pc_x[k] * fsh_127[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_z, fsg0_93, fsg0_95, fsg1_93, \
                         fsg1_95, fsh_126, fsh_127, fsh_129, fsh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_3 * pc_z[k] * fsh_126[k];

        t_171[k] = f_8 * fsg0_93[k]
                   - f_9 * fsg1_93[k]
                   + f_3 * pc_x[k] * fsh_129[k];

        t_172[k] = f_3 * pc_z[k] * fsh_127[k];

        t_173[k] = f_8 * fsg0_95[k]
                   - f_9 * fsg1_95[k]
                   + f_3 * pc_x[k] * fsh_131[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_z, fsg0_96, fsg0_98, fsg0_99, \
                         fsg1_96, fsg1_98, fsg1_99, fsh_129, fsh_132, fsh_134, \
                         fsh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * fsg0_96[k]
                   - f_7 * fsg1_96[k]
                   + f_3 * pc_x[k] * fsh_132[k];

        t_175[k] = f_3 * pc_z[k] * fsh_129[k];

        t_176[k] = f_6 * fsg0_98[k]
                   - f_7 * fsg1_98[k]
                   + f_3 * pc_x[k] * fsh_134[k];

        t_177[k] = f_6 * fsg0_99[k]
                   - f_7 * fsg1_99[k]
                   + f_3 * pc_x[k] * fsh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pc_x, pc_z, fsg0_100, fsg0_102, fsg0_103, \
                         fsg1_100, fsg1_102, fsg1_103, fsh_132, fsh_136, fsh_138, \
                         fsh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_4 * fsg0_100[k]
                   - f_5 * fsg1_100[k]
                   + f_3 * pc_x[k] * fsh_136[k];

        t_179[k] = f_3 * pc_z[k] * fsh_132[k];

        t_180[k] = f_4 * fsg0_102[k]
                   - f_5 * fsg1_102[k]
                   + f_3 * pc_x[k] * fsh_138[k];

        t_181[k] = f_4 * fsg0_103[k]
                   - f_5 * fsg1_103[k]
                   + f_3 * pc_x[k] * fsh_139[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, t_187, pc_x, fsg0_104, fsg1_104, \
                         fsh_140, fsh_141, fsh_142, fsh_143, fsh_144, \
                         fsh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_4 * fsg0_104[k]
                   - f_5 * fsg1_104[k]
                   + f_3 * pc_x[k] * fsh_140[k];

        t_183[k] = f_3 * pc_x[k] * fsh_141[k];

        t_184[k] = f_3 * pc_x[k] * fsh_142[k];

        t_185[k] = f_3 * pc_x[k] * fsh_143[k];

        t_186[k] = f_3 * pc_x[k] * fsh_144[k];

        t_187[k] = f_3 * pc_x[k] * fsh_145[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pc_x, pc_y, pc_z, dsh_78, fsg0_100, \
                         fsg1_100, fsh_141, fsh_142, fsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_3 * pc_x[k] * fsh_146[k];

        t_189[k] = f_0 * dsh_78[k]
                   + f_1 * fsg0_100[k]
                   - f_2 * fsg1_100[k]
                   + f_3 * pc_y[k] * fsh_141[k];

        t_190[k] = f_3 * pc_z[k] * fsh_141[k];

        t_191[k] = f_4 * fsg0_100[k]
                   - f_5 * fsg1_100[k]
                   + f_3 * pc_z[k] * fsh_142[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, dsh_83, fsg0_101, fsg0_102, \
                         fsg0_104, fsg1_101, fsg1_102, fsg1_104, fsh_143, fsh_144, \
                         fsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_6 * fsg0_101[k]
                   - f_7 * fsg1_101[k]
                   + f_3 * pc_z[k] * fsh_143[k];

        t_193[k] = f_8 * fsg0_102[k]
                   - f_9 * fsg1_102[k]
                   + f_3 * pc_z[k] * fsh_144[k];

        t_194[k] = f_0 * dsh_83[k]
                   + f_3 * pc_y[k] * fsh_146[k];

        t_195[k] = f_1 * fsg0_104[k]
                   - f_2 * fsg1_104[k]
                   + f_3 * pc_z[k] * fsh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_z, pc_x, pc_z, dsi0_84, dsi0_85, \
                         dsi0_87, dsi1_84, dsi1_85, dsi1_87, fsg0_107, fsg1_107, \
                         fsh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_z[k] * dsi0_84[k]
                   - f_10 * pc_z[k] * dsi1_84[k];

        t_197[k] = pa_z[k] * dsi0_85[k]
                   - f_10 * pc_z[k] * dsi1_85[k];

        t_198[k] = f_14 * fsg0_107[k]
                   - f_15 * fsg1_107[k]
                   + f_3 * pc_x[k] * fsh_149[k];

        t_199[k] = pa_z[k] * dsi0_87[k]
                   - f_10 * pc_z[k] * dsi1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pa_z, pc_x, pc_z, dsi0_90, dsi1_90, fsg0_109, \
                         fsg0_110, fsg1_109, fsg1_110, fsh_151, \
                         fsh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_8 * fsg0_109[k]
                   - f_9 * fsg1_109[k]
                   + f_3 * pc_x[k] * fsh_151[k];

        t_201[k] = f_8 * fsg0_110[k]
                   - f_9 * fsg1_110[k]
                   + f_3 * pc_x[k] * fsh_152[k];

        t_202[k] = pa_z[k] * dsi0_90[k]
                   - f_10 * pc_z[k] * dsi1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, fsg0_112, fsg0_113, fsg0_114, fsg1_112, \
                         fsg1_113, fsg1_114, fsh_154, fsh_155, \
                         fsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_6 * fsg0_112[k]
                   - f_7 * fsg1_112[k]
                   + f_3 * pc_x[k] * fsh_154[k];

        t_204[k] = f_6 * fsg0_113[k]
                   - f_7 * fsg1_113[k]
                   + f_3 * pc_x[k] * fsh_155[k];

        t_205[k] = f_6 * fsg0_114[k]
                   - f_7 * fsg1_114[k]
                   + f_3 * pc_x[k] * fsh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pa_z, pc_x, pc_z, dsi0_94, dsi1_94, fsg0_116, \
                         fsg0_117, fsg1_116, fsg1_117, fsh_158, \
                         fsh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_z[k] * dsi0_94[k]
                   - f_10 * pc_z[k] * dsi1_94[k];

        t_207[k] = f_4 * fsg0_116[k]
                   - f_5 * fsg1_116[k]
                   + f_3 * pc_x[k] * fsh_158[k];

        t_208[k] = f_4 * fsg0_117[k]
                   - f_5 * fsg1_117[k]
                   + f_3 * pc_x[k] * fsh_159[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pc_x, fsg0_118, fsg0_119, \
                         fsg1_118, fsg1_119, fsh_160, fsh_161, fsh_162, fsh_163, \
                         fsh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_4 * fsg0_118[k]
                   - f_5 * fsg1_118[k]
                   + f_3 * pc_x[k] * fsh_160[k];

        t_210[k] = f_4 * fsg0_119[k]
                   - f_5 * fsg1_119[k]
                   + f_3 * pc_x[k] * fsh_161[k];

        t_211[k] = f_3 * pc_x[k] * fsh_162[k];

        t_212[k] = f_3 * pc_x[k] * fsh_163[k];

        t_213[k] = f_3 * pc_x[k] * fsh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pa_z, pc_x, pc_z, dsi0_105, \
                         dsh_78, dsi1_105, fsh_162, fsh_165, fsh_166, \
                         fsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_x[k] * fsh_165[k];

        t_215[k] = f_3 * pc_x[k] * fsh_166[k];

        t_216[k] = f_3 * pc_x[k] * fsh_167[k];

        t_217[k] = pa_z[k] * dsi0_105[k]
                   - f_10 * pc_z[k] * dsi1_105[k];

        t_218[k] = f_11 * dsh_78[k]
                   + f_3 * pc_z[k] * fsh_162[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pa_z, pc_z, dsi0_107, dsi0_108, dsi0_109, \
                         dsh_79, dsh_80, dsh_81, dsi1_107, dsi1_108, \
                         dsi1_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = pa_z[k] * dsi0_107[k]
                   + f_12 * dsh_79[k]
                   - f_10 * pc_z[k] * dsi1_107[k];

        t_220[k] = pa_z[k] * dsi0_108[k]
                   + f_0 * dsh_80[k]
                   - f_10 * pc_z[k] * dsi1_108[k];

        t_221[k] = pa_z[k] * dsi0_109[k]
                   + f_13 * dsh_81[k]
                   - f_10 * pc_z[k] * dsi1_109[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_y, pc_y, pc_z, dsi0_140, dsh_83, dsh_104, \
                         dsi1_140, fsg0_119, fsg1_119, fsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_12 * dsh_104[k]
                   + f_3 * pc_y[k] * fsh_167[k];

        t_223[k] = f_11 * dsh_83[k]
                   + f_1 * fsg0_119[k]
                   - f_2 * fsg1_119[k]
                   + f_3 * pc_z[k] * fsh_167[k];

        t_224[k] = pa_y[k] * dsi0_140[k]
                   - f_10 * pc_y[k] * dsi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pc_x, pc_y, dsi0_142, dsi1_142, fsg0_121, \
                         fsg0_123, fsg1_121, fsg1_123, fsh_169, \
                         fsh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_14 * fsg0_121[k]
                   - f_15 * fsg1_121[k]
                   + f_3 * pc_x[k] * fsh_169[k];

        t_226[k] = pa_y[k] * dsi0_142[k]
                   - f_10 * pc_y[k] * dsi1_142[k];

        t_227[k] = f_8 * fsg0_123[k]
                   - f_9 * fsg1_123[k]
                   + f_3 * pc_x[k] * fsh_171[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_y, pc_x, pc_y, dsi0_145, dsi1_145, fsg0_124, \
                         fsg0_126, fsg1_124, fsg1_126, fsh_172, \
                         fsh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_8 * fsg0_124[k]
                   - f_9 * fsg1_124[k]
                   + f_3 * pc_x[k] * fsh_172[k];

        t_229[k] = pa_y[k] * dsi0_145[k]
                   - f_10 * pc_y[k] * dsi1_145[k];

        t_230[k] = f_6 * fsg0_126[k]
                   - f_7 * fsg1_126[k]
                   + f_3 * pc_x[k] * fsh_174[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pa_y, pc_x, pc_y, dsi0_149, dsi1_149, fsg0_127, \
                         fsg0_128, fsg1_127, fsg1_128, fsh_175, \
                         fsh_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_6 * fsg0_127[k]
                   - f_7 * fsg1_127[k]
                   + f_3 * pc_x[k] * fsh_175[k];

        t_232[k] = f_6 * fsg0_128[k]
                   - f_7 * fsg1_128[k]
                   + f_3 * pc_x[k] * fsh_176[k];

        t_233[k] = pa_y[k] * dsi0_149[k]
                   - f_10 * pc_y[k] * dsi1_149[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, fsg0_130, fsg0_131, fsg0_132, fsg1_130, \
                         fsg1_131, fsg1_132, fsh_178, fsh_179, \
                         fsh_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_4 * fsg0_130[k]
                   - f_5 * fsg1_130[k]
                   + f_3 * pc_x[k] * fsh_178[k];

        t_235[k] = f_4 * fsg0_131[k]
                   - f_5 * fsg1_131[k]
                   + f_3 * pc_x[k] * fsh_179[k];

        t_236[k] = f_4 * fsg0_132[k]
                   - f_5 * fsg1_132[k]
                   + f_3 * pc_x[k] * fsh_180[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pa_y, pc_x, pc_y, dsi0_154, \
                         dsi1_154, fsg0_133, fsg1_133, fsh_181, fsh_183, fsh_184, \
                         fsh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * fsg0_133[k]
                   - f_5 * fsg1_133[k]
                   + f_3 * pc_x[k] * fsh_181[k];

        t_238[k] = pa_y[k] * dsi0_154[k]
                   - f_10 * pc_y[k] * dsi1_154[k];

        t_239[k] = f_3 * pc_x[k] * fsh_183[k];

        t_240[k] = f_3 * pc_x[k] * fsh_184[k];

        t_241[k] = f_3 * pc_x[k] * fsh_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_y, pc_x, pc_y, dsi0_161, dsh_120, \
                         dsi1_161, fsh_186, fsh_187, fsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * pc_x[k] * fsh_186[k];

        t_243[k] = f_3 * pc_x[k] * fsh_187[k];

        t_244[k] = f_3 * pc_x[k] * fsh_188[k];

        t_245[k] = pa_y[k] * dsi0_161[k]
                   + f_16 * dsh_120[k]
                   - f_10 * pc_y[k] * dsi1_161[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_y, pc_y, pc_z, dsi0_163, dsi0_164, dsh_99, \
                         dsh_122, dsh_123, dsi1_163, dsi1_164, \
                         fsh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * dsh_99[k]
                   + f_3 * pc_z[k] * fsh_183[k];

        t_247[k] = pa_y[k] * dsi0_163[k]
                   + f_13 * dsh_122[k]
                   - f_10 * pc_y[k] * dsi1_163[k];

        t_248[k] = pa_y[k] * dsi0_164[k]
                   + f_0 * dsh_123[k]
                   - f_10 * pc_y[k] * dsi1_164[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_y, pc_y, dsi0_165, dsi0_167, dsh_124, \
                         dsh_125, dsi1_165, dsi1_167, fsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pa_y[k] * dsi0_165[k]
                   + f_12 * dsh_124[k]
                   - f_10 * pc_y[k] * dsi1_165[k];

        t_250[k] = f_11 * dsh_125[k]
                   + f_3 * pc_y[k] * fsh_188[k];

        t_251[k] = pa_y[k] * dsi0_167[k]
                   - f_10 * pc_y[k] * dsi1_167[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, pc_x, pc_y, fsg0_135, fsg0_137, \
                         fsg0_138, fsg1_135, fsg1_137, fsg1_138, fsh_189, fsh_191, \
                         fsh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * fsg0_135[k]
                   - f_2 * fsg1_135[k]
                   + f_3 * pc_x[k] * fsh_189[k];

        t_253[k] = f_3 * pc_y[k] * fsh_189[k];

        t_254[k] = f_14 * fsg0_137[k]
                   - f_15 * fsg1_137[k]
                   + f_3 * pc_x[k] * fsh_191[k];

        t_255[k] = f_8 * fsg0_138[k]
                   - f_9 * fsg1_138[k]
                   + f_3 * pc_x[k] * fsh_192[k];

        t_256[k] = f_3 * pc_y[k] * fsh_191[k];
    }
}

static auto
compute_prim_fsi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t dsh, const size_t fsg0,
                                                          const size_t fsg1, const size_t fsh,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsh_125 = buffer.data(dsh + 125);

    const auto *fsg0_140 = buffer.data(fsg0 + 140);
    const auto *fsg0_141 = buffer.data(fsg0 + 141);
    const auto *fsg0_142 = buffer.data(fsg0 + 142);
    const auto *fsg0_144 = buffer.data(fsg0 + 144);
    const auto *fsg0_145 = buffer.data(fsg0 + 145);
    const auto *fsg0_146 = buffer.data(fsg0 + 146);
    const auto *fsg0_147 = buffer.data(fsg0 + 147);
    const auto *fsg0_148 = buffer.data(fsg0 + 148);
    const auto *fsg0_149 = buffer.data(fsg0 + 149);

    const auto *fsg1_140 = buffer.data(fsg1 + 140);
    const auto *fsg1_141 = buffer.data(fsg1 + 141);
    const auto *fsg1_142 = buffer.data(fsg1 + 142);
    const auto *fsg1_144 = buffer.data(fsg1 + 144);
    const auto *fsg1_145 = buffer.data(fsg1 + 145);
    const auto *fsg1_146 = buffer.data(fsg1 + 146);
    const auto *fsg1_147 = buffer.data(fsg1 + 147);
    const auto *fsg1_148 = buffer.data(fsg1 + 148);
    const auto *fsg1_149 = buffer.data(fsg1 + 149);

    const auto *fsh_194 = buffer.data(fsh + 194);
    const auto *fsh_195 = buffer.data(fsh + 195);
    const auto *fsh_196 = buffer.data(fsh + 196);
    const auto *fsh_198 = buffer.data(fsh + 198);
    const auto *fsh_199 = buffer.data(fsh + 199);
    const auto *fsh_200 = buffer.data(fsh + 200);
    const auto *fsh_201 = buffer.data(fsh + 201);
    const auto *fsh_203 = buffer.data(fsh + 203);
    const auto *fsh_204 = buffer.data(fsh + 204);
    const auto *fsh_205 = buffer.data(fsh + 205);
    const auto *fsh_206 = buffer.data(fsh + 206);
    const auto *fsh_207 = buffer.data(fsh + 207);
    const auto *fsh_208 = buffer.data(fsh + 208);
    const auto *fsh_209 = buffer.data(fsh + 209);

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pc_x, pc_y, fsg0_140, fsg0_141, fsg0_142, \
                         fsg1_140, fsg1_141, fsg1_142, fsh_194, fsh_195, \
                         fsh_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_8 * fsg0_140[k]
                   - f_9 * fsg1_140[k]
                   + f_3 * pc_x[k] * fsh_194[k];

        t_258[k] = f_6 * fsg0_141[k]
                   - f_7 * fsg1_141[k]
                   + f_3 * pc_x[k] * fsh_195[k];

        t_259[k] = f_6 * fsg0_142[k]
                   - f_7 * fsg1_142[k]
                   + f_3 * pc_x[k] * fsh_196[k];

        t_260[k] = f_3 * pc_y[k] * fsh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, fsg0_144, fsg0_145, fsg0_146, fsg1_144, \
                         fsg1_145, fsg1_146, fsh_198, fsh_199, \
                         fsh_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_6 * fsg0_144[k]
                   - f_7 * fsg1_144[k]
                   + f_3 * pc_x[k] * fsh_198[k];

        t_262[k] = f_4 * fsg0_145[k]
                   - f_5 * fsg1_145[k]
                   + f_3 * pc_x[k] * fsh_199[k];

        t_263[k] = f_4 * fsg0_146[k]
                   - f_5 * fsg1_146[k]
                   + f_3 * pc_x[k] * fsh_200[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, fsg0_147, fsg0_149, \
                         fsg1_147, fsg1_149, fsh_198, fsh_201, fsh_203, fsh_204, \
                         fsh_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * fsg0_147[k]
                   - f_5 * fsg1_147[k]
                   + f_3 * pc_x[k] * fsh_201[k];

        t_265[k] = f_3 * pc_y[k] * fsh_198[k];

        t_266[k] = f_4 * fsg0_149[k]
                   - f_5 * fsg1_149[k]
                   + f_3 * pc_x[k] * fsh_203[k];

        t_267[k] = f_3 * pc_x[k] * fsh_204[k];

        t_268[k] = f_3 * pc_x[k] * fsh_205[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, pc_x, pc_y, fsg0_145, fsg1_145, \
                         fsh_204, fsh_206, fsh_207, fsh_208, fsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_3 * pc_x[k] * fsh_206[k];

        t_270[k] = f_3 * pc_x[k] * fsh_207[k];

        t_271[k] = f_3 * pc_x[k] * fsh_208[k];

        t_272[k] = f_3 * pc_x[k] * fsh_209[k];

        t_273[k] = f_1 * fsg0_145[k]
                   - f_2 * fsg1_145[k]
                   + f_3 * pc_y[k] * fsh_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pc_y, fsg0_146, fsg0_147, fsg0_148, fsg1_146, \
                         fsg1_147, fsg1_148, fsh_205, fsh_206, \
                         fsh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_14 * fsg0_146[k]
                   - f_15 * fsg1_146[k]
                   + f_3 * pc_y[k] * fsh_205[k];

        t_275[k] = f_8 * fsg0_147[k]
                   - f_9 * fsg1_147[k]
                   + f_3 * pc_y[k] * fsh_206[k];

        t_276[k] = f_6 * fsg0_148[k]
                   - f_7 * fsg1_148[k]
                   + f_3 * pc_y[k] * fsh_207[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pc_y, pc_z, dsh_125, fsg0_149, fsg1_149, \
                         fsh_208, fsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_4 * fsg0_149[k]
                   - f_5 * fsg1_149[k]
                   + f_3 * pc_y[k] * fsh_208[k];

        t_278[k] = f_3 * pc_y[k] * fsh_209[k];

        t_279[k] = f_0 * dsh_125[k]
                   + f_1 * fsg0_149[k]
                   - f_2 * fsg1_149[k]
                   + f_3 * pc_z[k] * fsh_209[k];
    }
}

auto
compute_prim_fsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsi0, const size_t dsh,
                                                   const size_t dsi1, const size_t fsg0,
                                                   const size_t fsg1, const size_t fsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, dsi0, dsh,
                                                              dsi1, fsg0, fsg1, fsh, ncols,
                                                              gamma, p, q);

    compute_prim_fsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, dsi0, dsh,
                                                              dsi1, fsg0, fsg1, fsh, ncols,
                                                              gamma, p, q);

    compute_prim_fsi_three_center_electron_repulsion_0_piece2(buffer, target, pc, dsh, fsg0,
                                                              fsg1, fsh, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
