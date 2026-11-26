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


#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fsh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsh0,
                                                          const size_t dsg, const size_t dsh1,
                                                          const size_t fsf0, const size_t fsf1,
                                                          const size_t fsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / gamma;
    const auto f_12 = 1.5 * p / (gamma * q);
    const auto f_13 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsh0_0 = buffer.data(dsh0 + 0);
    const auto *dsh0_3 = buffer.data(dsh0 + 3);
    const auto *dsh0_5 = buffer.data(dsh0 + 5);
    const auto *dsh0_6 = buffer.data(dsh0 + 6);
    const auto *dsh0_9 = buffer.data(dsh0 + 9);
    const auto *dsh0_15 = buffer.data(dsh0 + 15);
    const auto *dsh0_20 = buffer.data(dsh0 + 20);
    const auto *dsh0_24 = buffer.data(dsh0 + 24);
    const auto *dsh0_27 = buffer.data(dsh0 + 27);
    const auto *dsh0_42 = buffer.data(dsh0 + 42);
    const auto *dsh0_47 = buffer.data(dsh0 + 47);
    const auto *dsh0_51 = buffer.data(dsh0 + 51);
    const auto *dsh0_63 = buffer.data(dsh0 + 63);
    const auto *dsh0_66 = buffer.data(dsh0 + 66);
    const auto *dsh0_69 = buffer.data(dsh0 + 69);
    const auto *dsh0_78 = buffer.data(dsh0 + 78);
    const auto *dsh0_80 = buffer.data(dsh0 + 80);
    const auto *dsh0_81 = buffer.data(dsh0 + 81);
    const auto *dsh0_83 = buffer.data(dsh0 + 83);
    const auto *dsh0_99 = buffer.data(dsh0 + 99);
    const auto *dsh0_101 = buffer.data(dsh0 + 101);
    const auto *dsh0_102 = buffer.data(dsh0 + 102);
    const auto *dsh0_104 = buffer.data(dsh0 + 104);
    const auto *dsh0_105 = buffer.data(dsh0 + 105);
    const auto *dsh0_110 = buffer.data(dsh0 + 110);
    const auto *dsh0_114 = buffer.data(dsh0 + 114);
    const auto *dsh0_120 = buffer.data(dsh0 + 120);
    const auto *dsh0_121 = buffer.data(dsh0 + 121);
    const auto *dsh0_122 = buffer.data(dsh0 + 122);
    const auto *dsh0_123 = buffer.data(dsh0 + 123);
    const auto *dsh0_125 = buffer.data(dsh0 + 125);

    const auto *dsg_0 = buffer.data(dsg + 0);
    const auto *dsg_1 = buffer.data(dsg + 1);
    const auto *dsg_2 = buffer.data(dsg + 2);
    const auto *dsg_3 = buffer.data(dsg + 3);
    const auto *dsg_5 = buffer.data(dsg + 5);
    const auto *dsg_10 = buffer.data(dsg + 10);
    const auto *dsg_12 = buffer.data(dsg + 12);
    const auto *dsg_14 = buffer.data(dsg + 14);
    const auto *dsg_15 = buffer.data(dsg + 15);
    const auto *dsg_18 = buffer.data(dsg + 18);
    const auto *dsg_20 = buffer.data(dsg + 20);
    const auto *dsg_25 = buffer.data(dsg + 25);
    const auto *dsg_27 = buffer.data(dsg + 27);
    const auto *dsg_28 = buffer.data(dsg + 28);
    const auto *dsg_29 = buffer.data(dsg + 29);
    const auto *dsg_30 = buffer.data(dsg + 30);
    const auto *dsg_32 = buffer.data(dsg + 32);
    const auto *dsg_35 = buffer.data(dsg + 35);
    const auto *dsg_40 = buffer.data(dsg + 40);
    const auto *dsg_41 = buffer.data(dsg + 41);
    const auto *dsg_42 = buffer.data(dsg + 42);
    const auto *dsg_44 = buffer.data(dsg + 44);
    const auto *dsg_45 = buffer.data(dsg + 45);
    const auto *dsg_48 = buffer.data(dsg + 48);
    const auto *dsg_51 = buffer.data(dsg + 51);
    const auto *dsg_55 = buffer.data(dsg + 55);
    const auto *dsg_57 = buffer.data(dsg + 57);
    const auto *dsg_58 = buffer.data(dsg + 58);
    const auto *dsg_59 = buffer.data(dsg + 59);
    const auto *dsg_70 = buffer.data(dsg + 70);
    const auto *dsg_71 = buffer.data(dsg + 71);
    const auto *dsg_72 = buffer.data(dsg + 72);
    const auto *dsg_73 = buffer.data(dsg + 73);
    const auto *dsg_74 = buffer.data(dsg + 74);
    const auto *dsg_75 = buffer.data(dsg + 75);
    const auto *dsg_80 = buffer.data(dsg + 80);
    const auto *dsg_84 = buffer.data(dsg + 84);
    const auto *dsg_85 = buffer.data(dsg + 85);
    const auto *dsg_86 = buffer.data(dsg + 86);
    const auto *dsg_87 = buffer.data(dsg + 87);
    const auto *dsg_89 = buffer.data(dsg + 89);

    const auto *dsh1_0 = buffer.data(dsh1 + 0);
    const auto *dsh1_3 = buffer.data(dsh1 + 3);
    const auto *dsh1_5 = buffer.data(dsh1 + 5);
    const auto *dsh1_6 = buffer.data(dsh1 + 6);
    const auto *dsh1_9 = buffer.data(dsh1 + 9);
    const auto *dsh1_15 = buffer.data(dsh1 + 15);
    const auto *dsh1_20 = buffer.data(dsh1 + 20);
    const auto *dsh1_24 = buffer.data(dsh1 + 24);
    const auto *dsh1_27 = buffer.data(dsh1 + 27);
    const auto *dsh1_42 = buffer.data(dsh1 + 42);
    const auto *dsh1_47 = buffer.data(dsh1 + 47);
    const auto *dsh1_51 = buffer.data(dsh1 + 51);
    const auto *dsh1_63 = buffer.data(dsh1 + 63);
    const auto *dsh1_66 = buffer.data(dsh1 + 66);
    const auto *dsh1_69 = buffer.data(dsh1 + 69);
    const auto *dsh1_78 = buffer.data(dsh1 + 78);
    const auto *dsh1_80 = buffer.data(dsh1 + 80);
    const auto *dsh1_81 = buffer.data(dsh1 + 81);
    const auto *dsh1_83 = buffer.data(dsh1 + 83);
    const auto *dsh1_99 = buffer.data(dsh1 + 99);
    const auto *dsh1_101 = buffer.data(dsh1 + 101);
    const auto *dsh1_102 = buffer.data(dsh1 + 102);
    const auto *dsh1_104 = buffer.data(dsh1 + 104);
    const auto *dsh1_105 = buffer.data(dsh1 + 105);
    const auto *dsh1_110 = buffer.data(dsh1 + 110);
    const auto *dsh1_114 = buffer.data(dsh1 + 114);
    const auto *dsh1_120 = buffer.data(dsh1 + 120);
    const auto *dsh1_121 = buffer.data(dsh1 + 121);
    const auto *dsh1_122 = buffer.data(dsh1 + 122);
    const auto *dsh1_123 = buffer.data(dsh1 + 123);
    const auto *dsh1_125 = buffer.data(dsh1 + 125);

    const auto *fsf0_0 = buffer.data(fsf0 + 0);
    const auto *fsf0_1 = buffer.data(fsf0 + 1);
    const auto *fsf0_2 = buffer.data(fsf0 + 2);
    const auto *fsf0_6 = buffer.data(fsf0 + 6);
    const auto *fsf0_8 = buffer.data(fsf0 + 8);
    const auto *fsf0_9 = buffer.data(fsf0 + 9);
    const auto *fsf0_16 = buffer.data(fsf0 + 16);
    const auto *fsf0_17 = buffer.data(fsf0 + 17);
    const auto *fsf0_22 = buffer.data(fsf0 + 22);
    const auto *fsf0_27 = buffer.data(fsf0 + 27);
    const auto *fsf0_28 = buffer.data(fsf0 + 28);
    const auto *fsf0_29 = buffer.data(fsf0 + 29);
    const auto *fsf0_30 = buffer.data(fsf0 + 30);
    const auto *fsf0_32 = buffer.data(fsf0 + 32);
    const auto *fsf0_50 = buffer.data(fsf0 + 50);
    const auto *fsf0_51 = buffer.data(fsf0 + 51);
    const auto *fsf0_52 = buffer.data(fsf0 + 52);
    const auto *fsf0_60 = buffer.data(fsf0 + 60);
    const auto *fsf0_61 = buffer.data(fsf0 + 61);
    const auto *fsf0_63 = buffer.data(fsf0 + 63);
    const auto *fsf0_65 = buffer.data(fsf0 + 65);
    const auto *fsf0_66 = buffer.data(fsf0 + 66);

    const auto *fsf1_0 = buffer.data(fsf1 + 0);
    const auto *fsf1_1 = buffer.data(fsf1 + 1);
    const auto *fsf1_2 = buffer.data(fsf1 + 2);
    const auto *fsf1_6 = buffer.data(fsf1 + 6);
    const auto *fsf1_8 = buffer.data(fsf1 + 8);
    const auto *fsf1_9 = buffer.data(fsf1 + 9);
    const auto *fsf1_16 = buffer.data(fsf1 + 16);
    const auto *fsf1_17 = buffer.data(fsf1 + 17);
    const auto *fsf1_22 = buffer.data(fsf1 + 22);
    const auto *fsf1_27 = buffer.data(fsf1 + 27);
    const auto *fsf1_28 = buffer.data(fsf1 + 28);
    const auto *fsf1_29 = buffer.data(fsf1 + 29);
    const auto *fsf1_30 = buffer.data(fsf1 + 30);
    const auto *fsf1_32 = buffer.data(fsf1 + 32);
    const auto *fsf1_50 = buffer.data(fsf1 + 50);
    const auto *fsf1_51 = buffer.data(fsf1 + 51);
    const auto *fsf1_52 = buffer.data(fsf1 + 52);
    const auto *fsf1_60 = buffer.data(fsf1 + 60);
    const auto *fsf1_61 = buffer.data(fsf1 + 61);
    const auto *fsf1_63 = buffer.data(fsf1 + 63);
    const auto *fsf1_65 = buffer.data(fsf1 + 65);
    const auto *fsf1_66 = buffer.data(fsf1 + 66);

    const auto *fsg_0 = buffer.data(fsg + 0);
    const auto *fsg_1 = buffer.data(fsg + 1);
    const auto *fsg_2 = buffer.data(fsg + 2);
    const auto *fsg_3 = buffer.data(fsg + 3);
    const auto *fsg_5 = buffer.data(fsg + 5);
    const auto *fsg_6 = buffer.data(fsg + 6);
    const auto *fsg_9 = buffer.data(fsg + 9);
    const auto *fsg_10 = buffer.data(fsg + 10);
    const auto *fsg_12 = buffer.data(fsg + 12);
    const auto *fsg_13 = buffer.data(fsg + 13);
    const auto *fsg_14 = buffer.data(fsg + 14);
    const auto *fsg_15 = buffer.data(fsg + 15);
    const auto *fsg_16 = buffer.data(fsg + 16);
    const auto *fsg_18 = buffer.data(fsg + 18);
    const auto *fsg_20 = buffer.data(fsg + 20);
    const auto *fsg_21 = buffer.data(fsg + 21);
    const auto *fsg_25 = buffer.data(fsg + 25);
    const auto *fsg_26 = buffer.data(fsg + 26);
    const auto *fsg_27 = buffer.data(fsg + 27);
    const auto *fsg_28 = buffer.data(fsg + 28);
    const auto *fsg_29 = buffer.data(fsg + 29);
    const auto *fsg_30 = buffer.data(fsg + 30);
    const auto *fsg_32 = buffer.data(fsg + 32);
    const auto *fsg_34 = buffer.data(fsg + 34);
    const auto *fsg_35 = buffer.data(fsg + 35);
    const auto *fsg_39 = buffer.data(fsg + 39);
    const auto *fsg_40 = buffer.data(fsg + 40);
    const auto *fsg_41 = buffer.data(fsg + 41);
    const auto *fsg_42 = buffer.data(fsg + 42);
    const auto *fsg_43 = buffer.data(fsg + 43);
    const auto *fsg_44 = buffer.data(fsg + 44);
    const auto *fsg_45 = buffer.data(fsg + 45);
    const auto *fsg_46 = buffer.data(fsg + 46);
    const auto *fsg_47 = buffer.data(fsg + 47);
    const auto *fsg_48 = buffer.data(fsg + 48);
    const auto *fsg_50 = buffer.data(fsg + 50);
    const auto *fsg_51 = buffer.data(fsg + 51);
    const auto *fsg_55 = buffer.data(fsg + 55);
    const auto *fsg_57 = buffer.data(fsg + 57);
    const auto *fsg_58 = buffer.data(fsg + 58);
    const auto *fsg_59 = buffer.data(fsg + 59);
    const auto *fsg_60 = buffer.data(fsg + 60);
    const auto *fsg_62 = buffer.data(fsg + 62);
    const auto *fsg_63 = buffer.data(fsg + 63);
    const auto *fsg_65 = buffer.data(fsg + 65);
    const auto *fsg_70 = buffer.data(fsg + 70);
    const auto *fsg_71 = buffer.data(fsg + 71);
    const auto *fsg_72 = buffer.data(fsg + 72);
    const auto *fsg_73 = buffer.data(fsg + 73);
    const auto *fsg_74 = buffer.data(fsg + 74);
    const auto *fsg_75 = buffer.data(fsg + 75);
    const auto *fsg_76 = buffer.data(fsg + 76);
    const auto *fsg_77 = buffer.data(fsg + 77);
    const auto *fsg_78 = buffer.data(fsg + 78);
    const auto *fsg_79 = buffer.data(fsg + 79);
    const auto *fsg_80 = buffer.data(fsg + 80);
    const auto *fsg_84 = buffer.data(fsg + 84);
    const auto *fsg_85 = buffer.data(fsg + 85);
    const auto *fsg_86 = buffer.data(fsg + 86);
    const auto *fsg_87 = buffer.data(fsg + 87);
    const auto *fsg_89 = buffer.data(fsg + 89);
    const auto *fsg_90 = buffer.data(fsg + 90);
    const auto *fsg_91 = buffer.data(fsg + 91);
    const auto *fsg_93 = buffer.data(fsg + 93);
    const auto *fsg_95 = buffer.data(fsg + 95);
    const auto *fsg_96 = buffer.data(fsg + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dsg_0, fsf0_0, \
                         fsf1_0, fsg_0, fsg_1, fsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsg_0[k]
                 + f_1 * fsf0_0[k]
                 - f_2 * fsf1_0[k]
                 + f_3 * pc_x[k] * fsg_0[k];

        t_1[k] = f_3 * pc_y[k] * fsg_0[k];

        t_2[k] = f_3 * pc_z[k] * fsg_0[k];

        t_3[k] = f_4 * fsf0_0[k]
                 - f_5 * fsf1_0[k]
                 + f_3 * pc_y[k] * fsg_1[k];

        t_4[k] = f_3 * pc_y[k] * fsg_2[k];

        t_5[k] = f_4 * fsf0_0[k]
                 - f_5 * fsf1_0[k]
                 + f_3 * pc_z[k] * fsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, dsg_10, fsf0_1, fsf0_2, \
                         fsf1_1, fsf1_2, fsg_3, fsg_5, fsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * fsf0_1[k]
                 - f_7 * fsf1_1[k]
                 + f_3 * pc_y[k] * fsg_3[k];

        t_7[k] = f_3 * pc_z[k] * fsg_3[k];

        t_8[k] = f_3 * pc_y[k] * fsg_5[k];

        t_9[k] = f_6 * fsf0_2[k]
                 - f_7 * fsf1_2[k]
                 + f_3 * pc_z[k] * fsg_5[k];

        t_10[k] = f_0 * dsg_10[k]
                  + f_3 * pc_x[k] * fsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, dsg_12, dsg_14, fsg_6, \
                         fsg_9, fsg_12, fsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * fsg_6[k];

        t_12[k] = f_0 * dsg_12[k]
                  + f_3 * pc_x[k] * fsg_12[k];

        t_13[k] = f_3 * pc_y[k] * fsg_9[k];

        t_14[k] = f_0 * dsg_14[k]
                  + f_3 * pc_x[k] * fsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, fsf0_6, fsf0_8, fsf0_9, fsf1_6, \
                         fsf1_8, fsf1_9, fsg_10, fsg_12, fsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * fsf0_6[k]
                  - f_2 * fsf1_6[k]
                  + f_3 * pc_y[k] * fsg_10[k];

        t_16[k] = f_3 * pc_z[k] * fsg_10[k];

        t_17[k] = f_6 * fsf0_8[k]
                  - f_7 * fsf1_8[k]
                  + f_3 * pc_y[k] * fsg_12[k];

        t_18[k] = f_4 * fsf0_9[k]
                  - f_5 * fsf1_9[k]
                  + f_3 * pc_y[k] * fsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, dsh0_0, dsg_0, \
                         dsh1_0, fsf0_9, fsf1_9, fsg_14, fsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * fsg_14[k];

        t_20[k] = f_1 * fsf0_9[k]
                  - f_2 * fsf1_9[k]
                  + f_3 * pc_z[k] * fsg_14[k];

        t_21[k] = pa_y[k] * dsh0_0[k]
                  - f_8 * pc_y[k] * dsh1_0[k];

        t_22[k] = f_9 * dsg_0[k]
                  + f_3 * pc_y[k] * fsg_15[k];

        t_23[k] = f_3 * pc_z[k] * fsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, dsh0_3, dsh0_5, dsh0_6, \
                         dsg_1, dsg_3, dsh1_3, dsh1_5, dsh1_6, fsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * dsh0_3[k]
                  + f_10 * dsg_1[k]
                  - f_8 * pc_y[k] * dsh1_3[k];

        t_25[k] = f_3 * pc_z[k] * fsg_16[k];

        t_26[k] = pa_y[k] * dsh0_5[k]
                  - f_8 * pc_y[k] * dsh1_5[k];

        t_27[k] = pa_y[k] * dsh0_6[k]
                  + f_0 * dsg_3[k]
                  - f_8 * pc_y[k] * dsh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, dsh0_9, dsg_5, \
                         dsg_25, dsh1_9, fsg_18, fsg_20, fsg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * fsg_18[k];

        t_29[k] = f_9 * dsg_5[k]
                  + f_3 * pc_y[k] * fsg_20[k];

        t_30[k] = pa_y[k] * dsh0_9[k]
                  - f_8 * pc_y[k] * dsh1_9[k];

        t_31[k] = f_10 * dsg_25[k]
                  + f_3 * pc_x[k] * fsg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, dsg_27, dsg_28, dsg_29, fsg_21, \
                         fsg_27, fsg_28, fsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * fsg_21[k];

        t_33[k] = f_10 * dsg_27[k]
                  + f_3 * pc_x[k] * fsg_27[k];

        t_34[k] = f_10 * dsg_28[k]
                  + f_3 * pc_x[k] * fsg_28[k];

        t_35[k] = f_10 * dsg_29[k]
                  + f_3 * pc_x[k] * fsg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, dsg_10, fsf0_16, fsf0_17, \
                         fsf1_16, fsf1_17, fsg_25, fsg_26, fsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * dsg_10[k]
                  + f_1 * fsf0_16[k]
                  - f_2 * fsf1_16[k]
                  + f_3 * pc_y[k] * fsg_25[k];

        t_37[k] = f_3 * pc_z[k] * fsg_25[k];

        t_38[k] = f_4 * fsf0_16[k]
                  - f_5 * fsf1_16[k]
                  + f_3 * pc_z[k] * fsg_26[k];

        t_39[k] = f_6 * fsf0_17[k]
                  - f_7 * fsf1_17[k]
                  + f_3 * pc_z[k] * fsg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, dsh0_0, dsh0_20, \
                         dsg_14, dsh1_0, dsh1_20, fsg_29, fsg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * dsg_14[k]
                  + f_3 * pc_y[k] * fsg_29[k];

        t_41[k] = pa_y[k] * dsh0_20[k]
                  - f_8 * pc_y[k] * dsh1_20[k];

        t_42[k] = pa_z[k] * dsh0_0[k]
                  - f_8 * pc_z[k] * dsh1_0[k];

        t_43[k] = f_3 * pc_y[k] * fsg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, dsh0_3, dsh0_5, dsg_0, \
                         dsg_2, dsh1_3, dsh1_5, fsg_30, fsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * dsg_0[k]
                  + f_3 * pc_z[k] * fsg_30[k];

        t_45[k] = pa_z[k] * dsh0_3[k]
                  - f_8 * pc_z[k] * dsh1_3[k];

        t_46[k] = f_3 * pc_y[k] * fsg_32[k];

        t_47[k] = pa_z[k] * dsh0_5[k]
                  + f_10 * dsg_2[k]
                  - f_8 * pc_z[k] * dsh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, dsh0_6, dsh0_9, dsg_5, \
                         dsh1_6, dsh1_9, fsf0_22, fsf1_22, fsg_34, \
                         fsg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * dsh0_6[k]
                  - f_8 * pc_z[k] * dsh1_6[k];

        t_49[k] = f_4 * fsf0_22[k]
                  - f_5 * fsf1_22[k]
                  + f_3 * pc_y[k] * fsg_34[k];

        t_50[k] = f_3 * pc_y[k] * fsg_35[k];

        t_51[k] = pa_z[k] * dsh0_9[k]
                  + f_0 * dsg_5[k]
                  - f_8 * pc_z[k] * dsh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, dsg_40, dsg_41, dsg_42, \
                         dsg_44, fsg_39, fsg_40, fsg_41, fsg_42, \
                         fsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * dsg_40[k]
                  + f_3 * pc_x[k] * fsg_40[k];

        t_53[k] = f_10 * dsg_41[k]
                  + f_3 * pc_x[k] * fsg_41[k];

        t_54[k] = f_10 * dsg_42[k]
                  + f_3 * pc_x[k] * fsg_42[k];

        t_55[k] = f_3 * pc_y[k] * fsg_39[k];

        t_56[k] = f_10 * dsg_44[k]
                  + f_3 * pc_x[k] * fsg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, dsh0_15, dsh1_15, fsf0_27, \
                         fsf0_28, fsf1_27, fsf1_28, fsg_41, fsg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * dsh0_15[k]
                  - f_8 * pc_z[k] * dsh1_15[k];

        t_58[k] = f_11 * fsf0_27[k]
                  - f_12 * fsf1_27[k]
                  + f_3 * pc_y[k] * fsg_41[k];

        t_59[k] = f_6 * fsf0_28[k]
                  - f_7 * fsf1_28[k]
                  + f_3 * pc_y[k] * fsg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pc_x, pc_y, pc_z, dsh0_63, dsg_14, \
                         dsg_45, dsh1_63, fsf0_29, fsf1_29, fsg_43, \
                         fsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * fsf0_29[k]
                  - f_5 * fsf1_29[k]
                  + f_3 * pc_y[k] * fsg_43[k];

        t_61[k] = f_3 * pc_y[k] * fsg_44[k];

        t_62[k] = f_9 * dsg_14[k]
                  + f_1 * fsf0_29[k]
                  - f_2 * fsf1_29[k]
                  + f_3 * pc_z[k] * fsg_44[k];

        t_63[k] = pa_x[k] * dsh0_63[k]
                  + f_13 * dsg_45[k]
                  - f_8 * pc_x[k] * dsh1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pc_x, pc_y, pc_z, dsh0_66, dsg_15, \
                         dsg_48, dsh1_66, fsg_45, fsg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * dsg_15[k]
                  + f_3 * pc_y[k] * fsg_45[k];

        t_65[k] = f_3 * pc_z[k] * fsg_45[k];

        t_66[k] = pa_x[k] * dsh0_66[k]
                  + f_0 * dsg_48[k]
                  - f_8 * pc_x[k] * dsh1_66[k];

        t_67[k] = f_3 * pc_z[k] * fsg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_x, pc_x, pc_z, dsh0_69, dsg_51, dsh1_69, \
                         fsf0_30, fsf1_30, fsg_47, fsg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * fsf0_30[k]
                  - f_5 * fsf1_30[k]
                  + f_3 * pc_z[k] * fsg_47[k];

        t_69[k] = pa_x[k] * dsh0_69[k]
                  + f_10 * dsg_51[k]
                  - f_8 * pc_x[k] * dsh1_69[k];

        t_70[k] = f_3 * pc_z[k] * fsg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, dsg_20, dsg_55, fsf0_32, \
                         fsf1_32, fsg_50, fsg_51, fsg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * dsg_20[k]
                  + f_3 * pc_y[k] * fsg_50[k];

        t_72[k] = f_6 * fsf0_32[k]
                  - f_7 * fsf1_32[k]
                  + f_3 * pc_z[k] * fsg_50[k];

        t_73[k] = f_9 * dsg_55[k]
                  + f_3 * pc_x[k] * fsg_55[k];

        t_74[k] = f_3 * pc_z[k] * fsg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_x, pc_x, dsh0_78, dsg_57, dsg_58, dsg_59, \
                         dsh1_78, fsg_57, fsg_58, fsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * dsg_57[k]
                  + f_3 * pc_x[k] * fsg_57[k];

        t_76[k] = f_9 * dsg_58[k]
                  + f_3 * pc_x[k] * fsg_58[k];

        t_77[k] = f_9 * dsg_59[k]
                  + f_3 * pc_x[k] * fsg_59[k];

        t_78[k] = pa_x[k] * dsh0_78[k]
                  - f_8 * pc_x[k] * dsh1_78[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_x, pc_x, pc_y, pc_z, dsh0_80, dsh0_81, \
                         dsg_29, dsh1_80, dsh1_81, fsg_55, fsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * fsg_55[k];

        t_80[k] = pa_x[k] * dsh0_80[k]
                  - f_8 * pc_x[k] * dsh1_80[k];

        t_81[k] = pa_x[k] * dsh0_81[k]
                  - f_8 * pc_x[k] * dsh1_81[k];

        t_82[k] = f_10 * dsg_29[k]
                  + f_3 * pc_y[k] * fsg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pa_y, pc_x, pc_y, pc_z, dsh0_42, \
                         dsh0_83, dsg_15, dsg_30, dsh1_42, dsh1_83, \
                         fsg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pa_x[k] * dsh0_83[k]
                  - f_8 * pc_x[k] * dsh1_83[k];

        t_84[k] = pa_y[k] * dsh0_42[k]
                  - f_8 * pc_y[k] * dsh1_42[k];

        t_85[k] = f_9 * dsg_30[k]
                  + f_3 * pc_y[k] * fsg_60[k];

        t_86[k] = f_9 * dsg_15[k]
                  + f_3 * pc_z[k] * fsg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, dsh0_24, dsh0_27, \
                         dsh0_47, dsg_32, dsh1_24, dsh1_27, dsh1_47, \
                         fsg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * dsh0_24[k]
                  - f_8 * pc_z[k] * dsh1_24[k];

        t_88[k] = f_9 * dsg_32[k]
                  + f_3 * pc_y[k] * fsg_62[k];

        t_89[k] = pa_y[k] * dsh0_47[k]
                  - f_8 * pc_y[k] * dsh1_47[k];

        t_90[k] = pa_z[k] * dsh0_27[k]
                  - f_8 * pc_z[k] * dsh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, dsh0_51, dsg_18, \
                         dsg_35, dsg_70, dsh1_51, fsg_63, fsg_65, \
                         fsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * dsg_18[k]
                  + f_3 * pc_z[k] * fsg_63[k];

        t_92[k] = f_9 * dsg_35[k]
                  + f_3 * pc_y[k] * fsg_65[k];

        t_93[k] = pa_y[k] * dsh0_51[k]
                  - f_8 * pc_y[k] * dsh1_51[k];

        t_94[k] = f_9 * dsg_70[k]
                  + f_3 * pc_x[k] * fsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, dsg_71, dsg_72, dsg_73, dsg_74, fsg_71, \
                         fsg_72, fsg_73, fsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_9 * dsg_71[k]
                  + f_3 * pc_x[k] * fsg_71[k];

        t_96[k] = f_9 * dsg_72[k]
                  + f_3 * pc_x[k] * fsg_72[k];

        t_97[k] = f_9 * dsg_73[k]
                  + f_3 * pc_x[k] * fsg_73[k];

        t_98[k] = f_9 * dsg_74[k]
                  + f_3 * pc_x[k] * fsg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pc_x, pc_z, dsh0_99, dsh0_101, \
                         dsh0_102, dsg_25, dsh1_99, dsh1_101, dsh1_102, \
                         fsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * dsh0_99[k]
                  - f_8 * pc_x[k] * dsh1_99[k];

        t_100[k] = f_9 * dsg_25[k]
                   + f_3 * pc_z[k] * fsg_70[k];

        t_101[k] = pa_x[k] * dsh0_101[k]
                   - f_8 * pc_x[k] * dsh1_101[k];

        t_102[k] = pa_x[k] * dsh0_102[k]
                   - f_8 * pc_x[k] * dsh1_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_x, pc_x, pc_y, dsh0_104, dsh0_105, \
                         dsg_44, dsg_75, dsh1_104, dsh1_105, fsg_74, \
                         fsg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_9 * dsg_44[k]
                   + f_3 * pc_y[k] * fsg_74[k];

        t_104[k] = pa_x[k] * dsh0_104[k]
                   - f_8 * pc_x[k] * dsh1_104[k];

        t_105[k] = pa_x[k] * dsh0_105[k]
                   + f_13 * dsg_75[k]
                   - f_8 * pc_x[k] * dsh1_105[k];

        t_106[k] = f_3 * pc_y[k] * fsg_75[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, pc_z, dsg_30, fsf0_50, fsf1_50, fsg_75, \
                         fsg_76, fsg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_10 * dsg_30[k]
                   + f_3 * pc_z[k] * fsg_75[k];

        t_108[k] = f_4 * fsf0_50[k]
                   - f_5 * fsf1_50[k]
                   + f_3 * pc_y[k] * fsg_76[k];

        t_109[k] = f_3 * pc_y[k] * fsg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_x, pc_x, pc_y, dsh0_110, dsg_80, dsh1_110, \
                         fsf0_51, fsf0_52, fsf1_51, fsf1_52, fsg_78, \
                         fsg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_x[k] * dsh0_110[k]
                   + f_0 * dsg_80[k]
                   - f_8 * pc_x[k] * dsh1_110[k];

        t_111[k] = f_6 * fsf0_51[k]
                   - f_7 * fsf1_51[k]
                   + f_3 * pc_y[k] * fsg_78[k];

        t_112[k] = f_4 * fsf0_52[k]
                   - f_5 * fsf1_52[k]
                   + f_3 * pc_y[k] * fsg_79[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pc_x, pc_y, dsh0_114, dsg_84, \
                         dsg_85, dsg_86, dsh1_114, fsg_80, fsg_85, \
                         fsg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * pc_y[k] * fsg_80[k];

        t_114[k] = pa_x[k] * dsh0_114[k]
                   + f_10 * dsg_84[k]
                   - f_8 * pc_x[k] * dsh1_114[k];

        t_115[k] = f_9 * dsg_85[k]
                   + f_3 * pc_x[k] * fsg_85[k];

        t_116[k] = f_9 * dsg_86[k]
                   + f_3 * pc_x[k] * fsg_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pc_x, pc_y, dsh0_120, dsg_87, \
                         dsg_89, dsh1_120, fsg_84, fsg_87, fsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_9 * dsg_87[k]
                   + f_3 * pc_x[k] * fsg_87[k];

        t_118[k] = f_3 * pc_y[k] * fsg_84[k];

        t_119[k] = f_9 * dsg_89[k]
                   + f_3 * pc_x[k] * fsg_89[k];

        t_120[k] = pa_x[k] * dsh0_120[k]
                   - f_8 * pc_x[k] * dsh1_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pc_x, pc_y, dsh0_121, dsh0_122, \
                         dsh0_123, dsh1_121, dsh1_122, dsh1_123, \
                         fsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_x[k] * dsh0_121[k]
                   - f_8 * pc_x[k] * dsh1_121[k];

        t_122[k] = pa_x[k] * dsh0_122[k]
                   - f_8 * pc_x[k] * dsh1_122[k];

        t_123[k] = pa_x[k] * dsh0_123[k]
                   - f_8 * pc_x[k] * dsh1_123[k];

        t_124[k] = f_3 * pc_y[k] * fsg_89[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_x, pc_x, pc_z, dsh0_125, dsh1_125, \
                         fsf0_60, fsf0_61, fsf1_60, fsf1_61, fsg_90, \
                         fsg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_x[k] * dsh0_125[k]
                   - f_8 * pc_x[k] * dsh1_125[k];

        t_126[k] = f_1 * fsf0_60[k]
                   - f_2 * fsf1_60[k]
                   + f_3 * pc_x[k] * fsg_90[k];

        t_127[k] = f_11 * fsf0_61[k]
                   - f_12 * fsf1_61[k]
                   + f_3 * pc_x[k] * fsg_91[k];

        t_128[k] = f_3 * pc_z[k] * fsg_90[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pc_x, pc_z, fsf0_63, fsf0_65, fsf0_66, \
                         fsf1_63, fsf1_65, fsf1_66, fsg_91, fsg_93, fsg_95, \
                         fsg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_6 * fsf0_63[k]
                   - f_7 * fsf1_63[k]
                   + f_3 * pc_x[k] * fsg_93[k];

        t_130[k] = f_3 * pc_z[k] * fsg_91[k];

        t_131[k] = f_6 * fsf0_65[k]
                   - f_7 * fsf1_65[k]
                   + f_3 * pc_x[k] * fsg_95[k];

        t_132[k] = f_4 * fsf0_66[k]
                   - f_5 * fsf1_66[k]
                   + f_3 * pc_x[k] * fsg_96[k];
    }
}

static auto
compute_prim_fsh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsh0,
                                                          const size_t dsg, const size_t dsh1,
                                                          const size_t fsf0, const size_t fsf1,
                                                          const size_t fsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / gamma;
    const auto f_12 = 1.5 * p / (gamma * q);
    const auto f_13 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsh0_63 = buffer.data(dsh0 + 63);
    const auto *dsh0_64 = buffer.data(dsh0 + 64);
    const auto *dsh0_66 = buffer.data(dsh0 + 66);
    const auto *dsh0_69 = buffer.data(dsh0 + 69);
    const auto *dsh0_78 = buffer.data(dsh0 + 78);
    const auto *dsh0_80 = buffer.data(dsh0 + 80);
    const auto *dsh0_81 = buffer.data(dsh0 + 81);
    const auto *dsh0_105 = buffer.data(dsh0 + 105);
    const auto *dsh0_107 = buffer.data(dsh0 + 107);
    const auto *dsh0_110 = buffer.data(dsh0 + 110);
    const auto *dsh0_114 = buffer.data(dsh0 + 114);
    const auto *dsh0_120 = buffer.data(dsh0 + 120);
    const auto *dsh0_122 = buffer.data(dsh0 + 122);
    const auto *dsh0_123 = buffer.data(dsh0 + 123);
    const auto *dsh0_125 = buffer.data(dsh0 + 125);

    const auto *dsg_55 = buffer.data(dsg + 55);
    const auto *dsg_56 = buffer.data(dsg + 56);
    const auto *dsg_57 = buffer.data(dsg + 57);
    const auto *dsg_59 = buffer.data(dsg + 59);
    const auto *dsg_70 = buffer.data(dsg + 70);
    const auto *dsg_74 = buffer.data(dsg + 74);
    const auto *dsg_85 = buffer.data(dsg + 85);
    const auto *dsg_87 = buffer.data(dsg + 87);
    const auto *dsg_88 = buffer.data(dsg + 88);
    const auto *dsg_89 = buffer.data(dsg + 89);

    const auto *dsh1_63 = buffer.data(dsh1 + 63);
    const auto *dsh1_64 = buffer.data(dsh1 + 64);
    const auto *dsh1_66 = buffer.data(dsh1 + 66);
    const auto *dsh1_69 = buffer.data(dsh1 + 69);
    const auto *dsh1_78 = buffer.data(dsh1 + 78);
    const auto *dsh1_80 = buffer.data(dsh1 + 80);
    const auto *dsh1_81 = buffer.data(dsh1 + 81);
    const auto *dsh1_105 = buffer.data(dsh1 + 105);
    const auto *dsh1_107 = buffer.data(dsh1 + 107);
    const auto *dsh1_110 = buffer.data(dsh1 + 110);
    const auto *dsh1_114 = buffer.data(dsh1 + 114);
    const auto *dsh1_120 = buffer.data(dsh1 + 120);
    const auto *dsh1_122 = buffer.data(dsh1 + 122);
    const auto *dsh1_123 = buffer.data(dsh1 + 123);
    const auto *dsh1_125 = buffer.data(dsh1 + 125);

    const auto *fsf0_66 = buffer.data(fsf0 + 66);
    const auto *fsf0_67 = buffer.data(fsf0 + 67);
    const auto *fsf0_68 = buffer.data(fsf0 + 68);
    const auto *fsf0_69 = buffer.data(fsf0 + 69);
    const auto *fsf0_72 = buffer.data(fsf0 + 72);
    const auto *fsf0_74 = buffer.data(fsf0 + 74);
    const auto *fsf0_75 = buffer.data(fsf0 + 75);
    const auto *fsf0_77 = buffer.data(fsf0 + 77);
    const auto *fsf0_78 = buffer.data(fsf0 + 78);
    const auto *fsf0_79 = buffer.data(fsf0 + 79);
    const auto *fsf0_81 = buffer.data(fsf0 + 81);
    const auto *fsf0_83 = buffer.data(fsf0 + 83);
    const auto *fsf0_84 = buffer.data(fsf0 + 84);
    const auto *fsf0_86 = buffer.data(fsf0 + 86);
    const auto *fsf0_87 = buffer.data(fsf0 + 87);
    const auto *fsf0_88 = buffer.data(fsf0 + 88);
    const auto *fsf0_90 = buffer.data(fsf0 + 90);
    const auto *fsf0_92 = buffer.data(fsf0 + 92);
    const auto *fsf0_93 = buffer.data(fsf0 + 93);
    const auto *fsf0_95 = buffer.data(fsf0 + 95);
    const auto *fsf0_96 = buffer.data(fsf0 + 96);
    const auto *fsf0_97 = buffer.data(fsf0 + 97);
    const auto *fsf0_98 = buffer.data(fsf0 + 98);
    const auto *fsf0_99 = buffer.data(fsf0 + 99);

    const auto *fsf1_66 = buffer.data(fsf1 + 66);
    const auto *fsf1_67 = buffer.data(fsf1 + 67);
    const auto *fsf1_68 = buffer.data(fsf1 + 68);
    const auto *fsf1_69 = buffer.data(fsf1 + 69);
    const auto *fsf1_72 = buffer.data(fsf1 + 72);
    const auto *fsf1_74 = buffer.data(fsf1 + 74);
    const auto *fsf1_75 = buffer.data(fsf1 + 75);
    const auto *fsf1_77 = buffer.data(fsf1 + 77);
    const auto *fsf1_78 = buffer.data(fsf1 + 78);
    const auto *fsf1_79 = buffer.data(fsf1 + 79);
    const auto *fsf1_81 = buffer.data(fsf1 + 81);
    const auto *fsf1_83 = buffer.data(fsf1 + 83);
    const auto *fsf1_84 = buffer.data(fsf1 + 84);
    const auto *fsf1_86 = buffer.data(fsf1 + 86);
    const auto *fsf1_87 = buffer.data(fsf1 + 87);
    const auto *fsf1_88 = buffer.data(fsf1 + 88);
    const auto *fsf1_90 = buffer.data(fsf1 + 90);
    const auto *fsf1_92 = buffer.data(fsf1 + 92);
    const auto *fsf1_93 = buffer.data(fsf1 + 93);
    const auto *fsf1_95 = buffer.data(fsf1 + 95);
    const auto *fsf1_96 = buffer.data(fsf1 + 96);
    const auto *fsf1_97 = buffer.data(fsf1 + 97);
    const auto *fsf1_98 = buffer.data(fsf1 + 98);
    const auto *fsf1_99 = buffer.data(fsf1 + 99);

    const auto *fsg_93 = buffer.data(fsg + 93);
    const auto *fsg_98 = buffer.data(fsg + 98);
    const auto *fsg_99 = buffer.data(fsg + 99);
    const auto *fsg_100 = buffer.data(fsg + 100);
    const auto *fsg_101 = buffer.data(fsg + 101);
    const auto *fsg_102 = buffer.data(fsg + 102);
    const auto *fsg_103 = buffer.data(fsg + 103);
    const auto *fsg_104 = buffer.data(fsg + 104);
    const auto *fsg_107 = buffer.data(fsg + 107);
    const auto *fsg_109 = buffer.data(fsg + 109);
    const auto *fsg_110 = buffer.data(fsg + 110);
    const auto *fsg_112 = buffer.data(fsg + 112);
    const auto *fsg_113 = buffer.data(fsg + 113);
    const auto *fsg_114 = buffer.data(fsg + 114);
    const auto *fsg_115 = buffer.data(fsg + 115);
    const auto *fsg_116 = buffer.data(fsg + 116);
    const auto *fsg_117 = buffer.data(fsg + 117);
    const auto *fsg_118 = buffer.data(fsg + 118);
    const auto *fsg_119 = buffer.data(fsg + 119);
    const auto *fsg_121 = buffer.data(fsg + 121);
    const auto *fsg_123 = buffer.data(fsg + 123);
    const auto *fsg_124 = buffer.data(fsg + 124);
    const auto *fsg_126 = buffer.data(fsg + 126);
    const auto *fsg_127 = buffer.data(fsg + 127);
    const auto *fsg_128 = buffer.data(fsg + 128);
    const auto *fsg_130 = buffer.data(fsg + 130);
    const auto *fsg_131 = buffer.data(fsg + 131);
    const auto *fsg_132 = buffer.data(fsg + 132);
    const auto *fsg_133 = buffer.data(fsg + 133);
    const auto *fsg_134 = buffer.data(fsg + 134);
    const auto *fsg_135 = buffer.data(fsg + 135);
    const auto *fsg_137 = buffer.data(fsg + 137);
    const auto *fsg_138 = buffer.data(fsg + 138);
    const auto *fsg_140 = buffer.data(fsg + 140);
    const auto *fsg_141 = buffer.data(fsg + 141);
    const auto *fsg_142 = buffer.data(fsg + 142);
    const auto *fsg_144 = buffer.data(fsg + 144);
    const auto *fsg_145 = buffer.data(fsg + 145);
    const auto *fsg_146 = buffer.data(fsg + 146);
    const auto *fsg_147 = buffer.data(fsg + 147);
    const auto *fsg_148 = buffer.data(fsg + 148);
    const auto *fsg_149 = buffer.data(fsg + 149);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pc_x, pc_z, fsf0_68, fsf0_69, \
                         fsf1_68, fsf1_69, fsg_93, fsg_98, fsg_99, fsg_100, \
                         fsg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * fsg_93[k];

        t_134[k] = f_4 * fsf0_68[k]
                   - f_5 * fsf1_68[k]
                   + f_3 * pc_x[k] * fsg_98[k];

        t_135[k] = f_4 * fsf0_69[k]
                   - f_5 * fsf1_69[k]
                   + f_3 * pc_x[k] * fsg_99[k];

        t_136[k] = f_3 * pc_x[k] * fsg_100[k];

        t_137[k] = f_3 * pc_x[k] * fsg_101[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pc_x, pc_y, pc_z, dsg_55, fsf0_66, \
                         fsf1_66, fsg_100, fsg_102, fsg_103, fsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * pc_x[k] * fsg_102[k];

        t_139[k] = f_3 * pc_x[k] * fsg_103[k];

        t_140[k] = f_3 * pc_x[k] * fsg_104[k];

        t_141[k] = f_0 * dsg_55[k]
                   + f_1 * fsf0_66[k]
                   - f_2 * fsf1_66[k]
                   + f_3 * pc_y[k] * fsg_100[k];

        t_142[k] = f_3 * pc_z[k] * fsg_100[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_y, pc_z, dsg_59, fsf0_66, fsf0_67, \
                         fsf0_69, fsf1_66, fsf1_67, fsf1_69, fsg_101, fsg_102, \
                         fsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * fsf0_66[k]
                   - f_5 * fsf1_66[k]
                   + f_3 * pc_z[k] * fsg_101[k];

        t_144[k] = f_6 * fsf0_67[k]
                   - f_7 * fsf1_67[k]
                   + f_3 * pc_z[k] * fsg_102[k];

        t_145[k] = f_0 * dsg_59[k]
                   + f_3 * pc_y[k] * fsg_104[k];

        t_146[k] = f_1 * fsf0_69[k]
                   - f_2 * fsf1_69[k]
                   + f_3 * pc_z[k] * fsg_104[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pc_x, pc_z, dsh0_63, dsh0_64, \
                         dsh0_66, dsh1_63, dsh1_64, dsh1_66, fsf0_72, fsf1_72, \
                         fsg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pa_z[k] * dsh0_63[k]
                   - f_8 * pc_z[k] * dsh1_63[k];

        t_148[k] = pa_z[k] * dsh0_64[k]
                   - f_8 * pc_z[k] * dsh1_64[k];

        t_149[k] = f_11 * fsf0_72[k]
                   - f_12 * fsf1_72[k]
                   + f_3 * pc_x[k] * fsg_107[k];

        t_150[k] = pa_z[k] * dsh0_66[k]
                   - f_8 * pc_z[k] * dsh1_66[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_z, pc_x, pc_z, dsh0_69, dsh1_69, fsf0_74, \
                         fsf0_75, fsf1_74, fsf1_75, fsg_109, fsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_6 * fsf0_74[k]
                   - f_7 * fsf1_74[k]
                   + f_3 * pc_x[k] * fsg_109[k];

        t_152[k] = f_6 * fsf0_75[k]
                   - f_7 * fsf1_75[k]
                   + f_3 * pc_x[k] * fsg_110[k];

        t_153[k] = pa_z[k] * dsh0_69[k]
                   - f_8 * pc_z[k] * dsh1_69[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, fsf0_77, fsf0_78, fsf0_79, fsf1_77, \
                         fsf1_78, fsf1_79, fsg_112, fsg_113, fsg_114, \
                         fsg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_4 * fsf0_77[k]
                   - f_5 * fsf1_77[k]
                   + f_3 * pc_x[k] * fsg_112[k];

        t_155[k] = f_4 * fsf0_78[k]
                   - f_5 * fsf1_78[k]
                   + f_3 * pc_x[k] * fsg_113[k];

        t_156[k] = f_4 * fsf0_79[k]
                   - f_5 * fsf1_79[k]
                   + f_3 * pc_x[k] * fsg_114[k];

        t_157[k] = f_3 * pc_x[k] * fsg_115[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pc_x, pc_z, dsh0_78, \
                         dsh1_78, fsg_116, fsg_117, fsg_118, fsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * pc_x[k] * fsg_116[k];

        t_159[k] = f_3 * pc_x[k] * fsg_117[k];

        t_160[k] = f_3 * pc_x[k] * fsg_118[k];

        t_161[k] = f_3 * pc_x[k] * fsg_119[k];

        t_162[k] = pa_z[k] * dsh0_78[k]
                   - f_8 * pc_z[k] * dsh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_z, pc_z, dsh0_80, dsh0_81, dsg_55, dsg_56, \
                         dsg_57, dsh1_80, dsh1_81, fsg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * dsg_55[k]
                   + f_3 * pc_z[k] * fsg_115[k];

        t_164[k] = pa_z[k] * dsh0_80[k]
                   + f_10 * dsg_56[k]
                   - f_8 * pc_z[k] * dsh1_80[k];

        t_165[k] = pa_z[k] * dsh0_81[k]
                   + f_0 * dsg_57[k]
                   - f_8 * pc_z[k] * dsh1_81[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_y, pc_y, pc_z, dsh0_105, dsg_59, dsg_74, \
                         dsh1_105, fsf0_79, fsf1_79, fsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * dsg_74[k]
                   + f_3 * pc_y[k] * fsg_119[k];

        t_167[k] = f_9 * dsg_59[k]
                   + f_1 * fsf0_79[k]
                   - f_2 * fsf1_79[k]
                   + f_3 * pc_z[k] * fsg_119[k];

        t_168[k] = pa_y[k] * dsh0_105[k]
                   - f_8 * pc_y[k] * dsh1_105[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pa_y, pc_x, pc_y, dsh0_107, dsh1_107, fsf0_81, \
                         fsf0_83, fsf1_81, fsf1_83, fsg_121, fsg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_11 * fsf0_81[k]
                   - f_12 * fsf1_81[k]
                   + f_3 * pc_x[k] * fsg_121[k];

        t_170[k] = pa_y[k] * dsh0_107[k]
                   - f_8 * pc_y[k] * dsh1_107[k];

        t_171[k] = f_6 * fsf0_83[k]
                   - f_7 * fsf1_83[k]
                   + f_3 * pc_x[k] * fsg_123[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pa_y, pc_x, pc_y, dsh0_110, dsh1_110, fsf0_84, \
                         fsf0_86, fsf1_84, fsf1_86, fsg_124, fsg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_6 * fsf0_84[k]
                   - f_7 * fsf1_84[k]
                   + f_3 * pc_x[k] * fsg_124[k];

        t_173[k] = pa_y[k] * dsh0_110[k]
                   - f_8 * pc_y[k] * dsh1_110[k];

        t_174[k] = f_4 * fsf0_86[k]
                   - f_5 * fsf1_86[k]
                   + f_3 * pc_x[k] * fsg_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, dsh0_114, dsh1_114, \
                         fsf0_87, fsf0_88, fsf1_87, fsf1_88, fsg_127, fsg_128, \
                         fsg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_4 * fsf0_87[k]
                   - f_5 * fsf1_87[k]
                   + f_3 * pc_x[k] * fsg_127[k];

        t_176[k] = f_4 * fsf0_88[k]
                   - f_5 * fsf1_88[k]
                   + f_3 * pc_x[k] * fsg_128[k];

        t_177[k] = pa_y[k] * dsh0_114[k]
                   - f_8 * pc_y[k] * dsh1_114[k];

        t_178[k] = f_3 * pc_x[k] * fsg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pa_y, pc_x, pc_y, dsh0_120, \
                         dsg_85, dsh1_120, fsg_131, fsg_132, fsg_133, \
                         fsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_3 * pc_x[k] * fsg_131[k];

        t_180[k] = f_3 * pc_x[k] * fsg_132[k];

        t_181[k] = f_3 * pc_x[k] * fsg_133[k];

        t_182[k] = f_3 * pc_x[k] * fsg_134[k];

        t_183[k] = pa_y[k] * dsh0_120[k]
                   + f_13 * dsg_85[k]
                   - f_8 * pc_y[k] * dsh1_120[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_y, pc_y, pc_z, dsh0_122, dsh0_123, dsg_70, \
                         dsg_87, dsg_88, dsh1_122, dsh1_123, fsg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_10 * dsg_70[k]
                   + f_3 * pc_z[k] * fsg_130[k];

        t_185[k] = pa_y[k] * dsh0_122[k]
                   + f_0 * dsg_87[k]
                   - f_8 * pc_y[k] * dsh1_122[k];

        t_186[k] = pa_y[k] * dsh0_123[k]
                   + f_10 * dsg_88[k]
                   - f_8 * pc_y[k] * dsh1_123[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pa_y, pc_x, pc_y, dsh0_125, dsg_89, \
                         dsh1_125, fsf0_90, fsf1_90, fsg_134, fsg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_9 * dsg_89[k]
                   + f_3 * pc_y[k] * fsg_134[k];

        t_188[k] = pa_y[k] * dsh0_125[k]
                   - f_8 * pc_y[k] * dsh1_125[k];

        t_189[k] = f_1 * fsf0_90[k]
                   - f_2 * fsf1_90[k]
                   + f_3 * pc_x[k] * fsg_135[k];

        t_190[k] = f_3 * pc_y[k] * fsg_135[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pc_x, pc_y, fsf0_92, fsf0_93, fsf0_95, \
                         fsf1_92, fsf1_93, fsf1_95, fsg_137, fsg_138, \
                         fsg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_11 * fsf0_92[k]
                   - f_12 * fsf1_92[k]
                   + f_3 * pc_x[k] * fsg_137[k];

        t_192[k] = f_6 * fsf0_93[k]
                   - f_7 * fsf1_93[k]
                   + f_3 * pc_x[k] * fsg_138[k];

        t_193[k] = f_3 * pc_y[k] * fsg_137[k];

        t_194[k] = f_6 * fsf0_95[k]
                   - f_7 * fsf1_95[k]
                   + f_3 * pc_x[k] * fsg_140[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pc_x, pc_y, fsf0_96, fsf0_97, fsf0_99, \
                         fsf1_96, fsf1_97, fsf1_99, fsg_140, fsg_141, fsg_142, \
                         fsg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_4 * fsf0_96[k]
                   - f_5 * fsf1_96[k]
                   + f_3 * pc_x[k] * fsg_141[k];

        t_196[k] = f_4 * fsf0_97[k]
                   - f_5 * fsf1_97[k]
                   + f_3 * pc_x[k] * fsg_142[k];

        t_197[k] = f_3 * pc_y[k] * fsg_140[k];

        t_198[k] = f_4 * fsf0_99[k]
                   - f_5 * fsf1_99[k]
                   + f_3 * pc_x[k] * fsg_144[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, fsf0_96, \
                         fsf1_96, fsg_145, fsg_146, fsg_147, fsg_148, \
                         fsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_x[k] * fsg_145[k];

        t_200[k] = f_3 * pc_x[k] * fsg_146[k];

        t_201[k] = f_3 * pc_x[k] * fsg_147[k];

        t_202[k] = f_3 * pc_x[k] * fsg_148[k];

        t_203[k] = f_3 * pc_x[k] * fsg_149[k];

        t_204[k] = f_1 * fsf0_96[k]
                   - f_2 * fsf1_96[k]
                   + f_3 * pc_y[k] * fsg_145[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_y, fsf0_97, fsf0_98, fsf0_99, fsf1_97, \
                         fsf1_98, fsf1_99, fsg_146, fsg_147, fsg_148, \
                         fsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_11 * fsf0_97[k]
                   - f_12 * fsf1_97[k]
                   + f_3 * pc_y[k] * fsg_146[k];

        t_206[k] = f_6 * fsf0_98[k]
                   - f_7 * fsf1_98[k]
                   + f_3 * pc_y[k] * fsg_147[k];

        t_207[k] = f_4 * fsf0_99[k]
                   - f_5 * fsf1_99[k]
                   + f_3 * pc_y[k] * fsg_148[k];

        t_208[k] = f_3 * pc_y[k] * fsg_149[k];
    }

#pragma omp simd aligned(t_209, pc_z, dsg_89, fsf0_99, fsf1_99, \
                         fsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_0 * dsg_89[k]
                   + f_1 * fsf0_99[k]
                   - f_2 * fsf1_99[k]
                   + f_3 * pc_z[k] * fsg_149[k];
    }
}

auto
compute_prim_fsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsh0, const size_t dsg,
                                                   const size_t dsh1, const size_t fsf0,
                                                   const size_t fsf1, const size_t fsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fsh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, dsh0, dsg,
                                                              dsh1, fsf0, fsf1, fsg, ncols,
                                                              gamma, p, q);

    compute_prim_fsh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, dsh0, dsg,
                                                              dsh1, fsf0, fsf1, fsg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
