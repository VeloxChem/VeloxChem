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


#include "SimdThreeCenterElectronRepulsionVrrRecGSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gsl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsl0,
                                                          const size_t fsk, const size_t fsl1,
                                                          const size_t gsi0, const size_t gsi1,
                                                          const size_t gsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *fsl0_0 = buffer.data(fsl0 + 0);
    const auto *fsl0_3 = buffer.data(fsl0 + 3);
    const auto *fsl0_5 = buffer.data(fsl0 + 5);
    const auto *fsl0_6 = buffer.data(fsl0 + 6);
    const auto *fsl0_9 = buffer.data(fsl0 + 9);
    const auto *fsl0_10 = buffer.data(fsl0 + 10);
    const auto *fsl0_14 = buffer.data(fsl0 + 14);
    const auto *fsl0_15 = buffer.data(fsl0 + 15);
    const auto *fsl0_20 = buffer.data(fsl0 + 20);
    const auto *fsl0_21 = buffer.data(fsl0 + 21);
    const auto *fsl0_27 = buffer.data(fsl0 + 27);
    const auto *fsl0_36 = buffer.data(fsl0 + 36);
    const auto *fsl0_44 = buffer.data(fsl0 + 44);

    const auto *fsk_0 = buffer.data(fsk + 0);
    const auto *fsk_1 = buffer.data(fsk + 1);
    const auto *fsk_2 = buffer.data(fsk + 2);
    const auto *fsk_3 = buffer.data(fsk + 3);
    const auto *fsk_5 = buffer.data(fsk + 5);
    const auto *fsk_6 = buffer.data(fsk + 6);
    const auto *fsk_9 = buffer.data(fsk + 9);
    const auto *fsk_10 = buffer.data(fsk + 10);
    const auto *fsk_14 = buffer.data(fsk + 14);
    const auto *fsk_15 = buffer.data(fsk + 15);
    const auto *fsk_20 = buffer.data(fsk + 20);
    const auto *fsk_28 = buffer.data(fsk + 28);
    const auto *fsk_30 = buffer.data(fsk + 30);
    const auto *fsk_31 = buffer.data(fsk + 31);
    const auto *fsk_32 = buffer.data(fsk + 32);
    const auto *fsk_33 = buffer.data(fsk + 33);
    const auto *fsk_35 = buffer.data(fsk + 35);
    const auto *fsk_64 = buffer.data(fsk + 64);
    const auto *fsk_66 = buffer.data(fsk + 66);
    const auto *fsk_67 = buffer.data(fsk + 67);
    const auto *fsk_68 = buffer.data(fsk + 68);
    const auto *fsk_69 = buffer.data(fsk + 69);
    const auto *fsk_70 = buffer.data(fsk + 70);
    const auto *fsk_71 = buffer.data(fsk + 71);
    const auto *fsk_100 = buffer.data(fsk + 100);
    const auto *fsk_101 = buffer.data(fsk + 101);
    const auto *fsk_102 = buffer.data(fsk + 102);
    const auto *fsk_103 = buffer.data(fsk + 103);
    const auto *fsk_104 = buffer.data(fsk + 104);
    const auto *fsk_105 = buffer.data(fsk + 105);
    const auto *fsk_107 = buffer.data(fsk + 107);

    const auto *fsl1_0 = buffer.data(fsl1 + 0);
    const auto *fsl1_3 = buffer.data(fsl1 + 3);
    const auto *fsl1_5 = buffer.data(fsl1 + 5);
    const auto *fsl1_6 = buffer.data(fsl1 + 6);
    const auto *fsl1_9 = buffer.data(fsl1 + 9);
    const auto *fsl1_10 = buffer.data(fsl1 + 10);
    const auto *fsl1_14 = buffer.data(fsl1 + 14);
    const auto *fsl1_15 = buffer.data(fsl1 + 15);
    const auto *fsl1_20 = buffer.data(fsl1 + 20);
    const auto *fsl1_21 = buffer.data(fsl1 + 21);
    const auto *fsl1_27 = buffer.data(fsl1 + 27);
    const auto *fsl1_36 = buffer.data(fsl1 + 36);
    const auto *fsl1_44 = buffer.data(fsl1 + 44);

    const auto *gsi0_0 = buffer.data(gsi0 + 0);
    const auto *gsi0_1 = buffer.data(gsi0 + 1);
    const auto *gsi0_2 = buffer.data(gsi0 + 2);
    const auto *gsi0_3 = buffer.data(gsi0 + 3);
    const auto *gsi0_5 = buffer.data(gsi0 + 5);
    const auto *gsi0_6 = buffer.data(gsi0 + 6);
    const auto *gsi0_8 = buffer.data(gsi0 + 8);
    const auto *gsi0_9 = buffer.data(gsi0 + 9);
    const auto *gsi0_10 = buffer.data(gsi0 + 10);
    const auto *gsi0_12 = buffer.data(gsi0 + 12);
    const auto *gsi0_13 = buffer.data(gsi0 + 13);
    const auto *gsi0_14 = buffer.data(gsi0 + 14);
    const auto *gsi0_21 = buffer.data(gsi0 + 21);
    const auto *gsi0_23 = buffer.data(gsi0 + 23);
    const auto *gsi0_24 = buffer.data(gsi0 + 24);
    const auto *gsi0_25 = buffer.data(gsi0 + 25);
    const auto *gsi0_26 = buffer.data(gsi0 + 26);
    const auto *gsi0_27 = buffer.data(gsi0 + 27);
    const auto *gsi0_31 = buffer.data(gsi0 + 31);
    const auto *gsi0_34 = buffer.data(gsi0 + 34);
    const auto *gsi0_35 = buffer.data(gsi0 + 35);
    const auto *gsi0_38 = buffer.data(gsi0 + 38);
    const auto *gsi0_39 = buffer.data(gsi0 + 39);
    const auto *gsi0_40 = buffer.data(gsi0 + 40);
    const auto *gsi0_49 = buffer.data(gsi0 + 49);
    const auto *gsi0_50 = buffer.data(gsi0 + 50);
    const auto *gsi0_51 = buffer.data(gsi0 + 51);
    const auto *gsi0_52 = buffer.data(gsi0 + 52);
    const auto *gsi0_53 = buffer.data(gsi0 + 53);
    const auto *gsi0_58 = buffer.data(gsi0 + 58);
    const auto *gsi0_60 = buffer.data(gsi0 + 60);
    const auto *gsi0_61 = buffer.data(gsi0 + 61);
    const auto *gsi0_63 = buffer.data(gsi0 + 63);
    const auto *gsi0_64 = buffer.data(gsi0 + 64);
    const auto *gsi0_65 = buffer.data(gsi0 + 65);
    const auto *gsi0_67 = buffer.data(gsi0 + 67);
    const auto *gsi0_68 = buffer.data(gsi0 + 68);
    const auto *gsi0_69 = buffer.data(gsi0 + 69);
    const auto *gsi0_70 = buffer.data(gsi0 + 70);
    const auto *gsi0_78 = buffer.data(gsi0 + 78);
    const auto *gsi0_79 = buffer.data(gsi0 + 79);

    const auto *gsi1_0 = buffer.data(gsi1 + 0);
    const auto *gsi1_1 = buffer.data(gsi1 + 1);
    const auto *gsi1_2 = buffer.data(gsi1 + 2);
    const auto *gsi1_3 = buffer.data(gsi1 + 3);
    const auto *gsi1_5 = buffer.data(gsi1 + 5);
    const auto *gsi1_6 = buffer.data(gsi1 + 6);
    const auto *gsi1_8 = buffer.data(gsi1 + 8);
    const auto *gsi1_9 = buffer.data(gsi1 + 9);
    const auto *gsi1_10 = buffer.data(gsi1 + 10);
    const auto *gsi1_12 = buffer.data(gsi1 + 12);
    const auto *gsi1_13 = buffer.data(gsi1 + 13);
    const auto *gsi1_14 = buffer.data(gsi1 + 14);
    const auto *gsi1_21 = buffer.data(gsi1 + 21);
    const auto *gsi1_23 = buffer.data(gsi1 + 23);
    const auto *gsi1_24 = buffer.data(gsi1 + 24);
    const auto *gsi1_25 = buffer.data(gsi1 + 25);
    const auto *gsi1_26 = buffer.data(gsi1 + 26);
    const auto *gsi1_27 = buffer.data(gsi1 + 27);
    const auto *gsi1_31 = buffer.data(gsi1 + 31);
    const auto *gsi1_34 = buffer.data(gsi1 + 34);
    const auto *gsi1_35 = buffer.data(gsi1 + 35);
    const auto *gsi1_38 = buffer.data(gsi1 + 38);
    const auto *gsi1_39 = buffer.data(gsi1 + 39);
    const auto *gsi1_40 = buffer.data(gsi1 + 40);
    const auto *gsi1_49 = buffer.data(gsi1 + 49);
    const auto *gsi1_50 = buffer.data(gsi1 + 50);
    const auto *gsi1_51 = buffer.data(gsi1 + 51);
    const auto *gsi1_52 = buffer.data(gsi1 + 52);
    const auto *gsi1_53 = buffer.data(gsi1 + 53);
    const auto *gsi1_58 = buffer.data(gsi1 + 58);
    const auto *gsi1_60 = buffer.data(gsi1 + 60);
    const auto *gsi1_61 = buffer.data(gsi1 + 61);
    const auto *gsi1_63 = buffer.data(gsi1 + 63);
    const auto *gsi1_64 = buffer.data(gsi1 + 64);
    const auto *gsi1_65 = buffer.data(gsi1 + 65);
    const auto *gsi1_67 = buffer.data(gsi1 + 67);
    const auto *gsi1_68 = buffer.data(gsi1 + 68);
    const auto *gsi1_69 = buffer.data(gsi1 + 69);
    const auto *gsi1_70 = buffer.data(gsi1 + 70);
    const auto *gsi1_78 = buffer.data(gsi1 + 78);
    const auto *gsi1_79 = buffer.data(gsi1 + 79);

    const auto *gsk_0 = buffer.data(gsk + 0);
    const auto *gsk_1 = buffer.data(gsk + 1);
    const auto *gsk_2 = buffer.data(gsk + 2);
    const auto *gsk_3 = buffer.data(gsk + 3);
    const auto *gsk_5 = buffer.data(gsk + 5);
    const auto *gsk_6 = buffer.data(gsk + 6);
    const auto *gsk_8 = buffer.data(gsk + 8);
    const auto *gsk_9 = buffer.data(gsk + 9);
    const auto *gsk_10 = buffer.data(gsk + 10);
    const auto *gsk_12 = buffer.data(gsk + 12);
    const auto *gsk_13 = buffer.data(gsk + 13);
    const auto *gsk_14 = buffer.data(gsk + 14);
    const auto *gsk_15 = buffer.data(gsk + 15);
    const auto *gsk_17 = buffer.data(gsk + 17);
    const auto *gsk_18 = buffer.data(gsk + 18);
    const auto *gsk_19 = buffer.data(gsk + 19);
    const auto *gsk_20 = buffer.data(gsk + 20);
    const auto *gsk_21 = buffer.data(gsk + 21);
    const auto *gsk_27 = buffer.data(gsk + 27);
    const auto *gsk_28 = buffer.data(gsk + 28);
    const auto *gsk_30 = buffer.data(gsk + 30);
    const auto *gsk_31 = buffer.data(gsk + 31);
    const auto *gsk_32 = buffer.data(gsk + 32);
    const auto *gsk_33 = buffer.data(gsk + 33);
    const auto *gsk_34 = buffer.data(gsk + 34);
    const auto *gsk_35 = buffer.data(gsk + 35);
    const auto *gsk_36 = buffer.data(gsk + 36);
    const auto *gsk_37 = buffer.data(gsk + 37);
    const auto *gsk_39 = buffer.data(gsk + 39);
    const auto *gsk_41 = buffer.data(gsk + 41);
    const auto *gsk_42 = buffer.data(gsk + 42);
    const auto *gsk_43 = buffer.data(gsk + 43);
    const auto *gsk_45 = buffer.data(gsk + 45);
    const auto *gsk_46 = buffer.data(gsk + 46);
    const auto *gsk_47 = buffer.data(gsk + 47);
    const auto *gsk_48 = buffer.data(gsk + 48);
    const auto *gsk_50 = buffer.data(gsk + 50);
    const auto *gsk_51 = buffer.data(gsk + 51);
    const auto *gsk_52 = buffer.data(gsk + 52);
    const auto *gsk_53 = buffer.data(gsk + 53);
    const auto *gsk_54 = buffer.data(gsk + 54);
    const auto *gsk_56 = buffer.data(gsk + 56);
    const auto *gsk_57 = buffer.data(gsk + 57);
    const auto *gsk_64 = buffer.data(gsk + 64);
    const auto *gsk_65 = buffer.data(gsk + 65);
    const auto *gsk_66 = buffer.data(gsk + 66);
    const auto *gsk_67 = buffer.data(gsk + 67);
    const auto *gsk_68 = buffer.data(gsk + 68);
    const auto *gsk_69 = buffer.data(gsk + 69);
    const auto *gsk_70 = buffer.data(gsk + 70);
    const auto *gsk_71 = buffer.data(gsk + 71);
    const auto *gsk_72 = buffer.data(gsk + 72);
    const auto *gsk_74 = buffer.data(gsk + 74);
    const auto *gsk_76 = buffer.data(gsk + 76);
    const auto *gsk_77 = buffer.data(gsk + 77);
    const auto *gsk_79 = buffer.data(gsk + 79);
    const auto *gsk_80 = buffer.data(gsk + 80);
    const auto *gsk_81 = buffer.data(gsk + 81);
    const auto *gsk_83 = buffer.data(gsk + 83);
    const auto *gsk_84 = buffer.data(gsk + 84);
    const auto *gsk_85 = buffer.data(gsk + 85);
    const auto *gsk_86 = buffer.data(gsk + 86);
    const auto *gsk_88 = buffer.data(gsk + 88);
    const auto *gsk_89 = buffer.data(gsk + 89);
    const auto *gsk_90 = buffer.data(gsk + 90);
    const auto *gsk_91 = buffer.data(gsk + 91);
    const auto *gsk_92 = buffer.data(gsk + 92);
    const auto *gsk_99 = buffer.data(gsk + 99);
    const auto *gsk_100 = buffer.data(gsk + 100);
    const auto *gsk_101 = buffer.data(gsk + 101);
    const auto *gsk_102 = buffer.data(gsk + 102);
    const auto *gsk_103 = buffer.data(gsk + 103);
    const auto *gsk_104 = buffer.data(gsk + 104);
    const auto *gsk_105 = buffer.data(gsk + 105);
    const auto *gsk_107 = buffer.data(gsk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fsk_0, gsi0_0, \
                         gsi1_0, gsk_0, gsk_1, gsk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsk_0[k]
                 + f_1 * gsi0_0[k]
                 - f_2 * gsi1_0[k]
                 + f_3 * pc_x[k] * gsk_0[k];

        t_1[k] = f_3 * pc_y[k] * gsk_0[k];

        t_2[k] = f_3 * pc_z[k] * gsk_0[k];

        t_3[k] = f_4 * gsi0_0[k]
                 - f_5 * gsi1_0[k]
                 + f_3 * pc_y[k] * gsk_1[k];

        t_4[k] = f_3 * pc_y[k] * gsk_2[k];

        t_5[k] = f_4 * gsi0_0[k]
                 - f_5 * gsi1_0[k]
                 + f_3 * pc_z[k] * gsk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, gsi0_1, gsi0_2, gsi0_3, gsi1_1, \
                         gsi1_2, gsi1_3, gsk_3, gsk_5, gsk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gsi0_1[k]
                 - f_7 * gsi1_1[k]
                 + f_3 * pc_y[k] * gsk_3[k];

        t_7[k] = f_3 * pc_z[k] * gsk_3[k];

        t_8[k] = f_3 * pc_y[k] * gsk_5[k];

        t_9[k] = f_6 * gsi0_2[k]
                 - f_7 * gsi1_2[k]
                 + f_3 * pc_z[k] * gsk_5[k];

        t_10[k] = f_8 * gsi0_3[k]
                  - f_9 * gsi1_3[k]
                  + f_3 * pc_y[k] * gsk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, gsi0_5, gsi0_6, \
                         gsi1_5, gsi1_6, gsk_6, gsk_8, gsk_9, gsk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * gsk_6[k];

        t_12[k] = f_4 * gsi0_5[k]
                  - f_5 * gsi1_5[k]
                  + f_3 * pc_y[k] * gsk_8[k];

        t_13[k] = f_3 * pc_y[k] * gsk_9[k];

        t_14[k] = f_8 * gsi0_5[k]
                  - f_9 * gsi1_5[k]
                  + f_3 * pc_z[k] * gsk_9[k];

        t_15[k] = f_10 * gsi0_6[k]
                  - f_11 * gsi1_6[k]
                  + f_3 * pc_y[k] * gsk_10[k];

        t_16[k] = f_3 * pc_z[k] * gsk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, gsi0_8, gsi0_9, gsi1_8, gsi1_9, \
                         gsk_12, gsk_13, gsk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * gsi0_8[k]
                  - f_7 * gsi1_8[k]
                  + f_3 * pc_y[k] * gsk_12[k];

        t_18[k] = f_4 * gsi0_9[k]
                  - f_5 * gsi1_9[k]
                  + f_3 * pc_y[k] * gsk_13[k];

        t_19[k] = f_3 * pc_y[k] * gsk_14[k];

        t_20[k] = f_10 * gsi0_9[k]
                  - f_11 * gsi1_9[k]
                  + f_3 * pc_z[k] * gsk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, gsi0_10, gsi0_12, gsi0_13, \
                         gsi1_10, gsi1_12, gsi1_13, gsk_15, gsk_17, \
                         gsk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * gsi0_10[k]
                  - f_13 * gsi1_10[k]
                  + f_3 * pc_y[k] * gsk_15[k];

        t_22[k] = f_3 * pc_z[k] * gsk_15[k];

        t_23[k] = f_8 * gsi0_12[k]
                  - f_9 * gsi1_12[k]
                  + f_3 * pc_y[k] * gsk_17[k];

        t_24[k] = f_6 * gsi0_13[k]
                  - f_7 * gsi1_13[k]
                  + f_3 * pc_y[k] * gsk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, fsk_28, gsi0_14, \
                         gsi1_14, gsk_19, gsk_20, gsk_21, gsk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * gsi0_14[k]
                  - f_5 * gsi1_14[k]
                  + f_3 * pc_y[k] * gsk_19[k];

        t_26[k] = f_3 * pc_y[k] * gsk_20[k];

        t_27[k] = f_12 * gsi0_14[k]
                  - f_13 * gsi1_14[k]
                  + f_3 * pc_z[k] * gsk_20[k];

        t_28[k] = f_0 * fsk_28[k]
                  + f_3 * pc_x[k] * gsk_28[k];

        t_29[k] = f_3 * pc_z[k] * gsk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, fsk_30, fsk_31, fsk_32, \
                         fsk_33, gsk_27, gsk_30, gsk_31, gsk_32, \
                         gsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * fsk_30[k]
                  + f_3 * pc_x[k] * gsk_30[k];

        t_31[k] = f_0 * fsk_31[k]
                  + f_3 * pc_x[k] * gsk_31[k];

        t_32[k] = f_0 * fsk_32[k]
                  + f_3 * pc_x[k] * gsk_32[k];

        t_33[k] = f_0 * fsk_33[k]
                  + f_3 * pc_x[k] * gsk_33[k];

        t_34[k] = f_3 * pc_y[k] * gsk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, fsk_35, gsi0_21, gsi0_23, \
                         gsi1_21, gsi1_23, gsk_28, gsk_30, gsk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * fsk_35[k]
                  + f_3 * pc_x[k] * gsk_35[k];

        t_36[k] = f_1 * gsi0_21[k]
                  - f_2 * gsi1_21[k]
                  + f_3 * pc_y[k] * gsk_28[k];

        t_37[k] = f_3 * pc_z[k] * gsk_28[k];

        t_38[k] = f_12 * gsi0_23[k]
                  - f_13 * gsi1_23[k]
                  + f_3 * pc_y[k] * gsk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, gsi0_24, gsi0_25, gsi0_26, gsi1_24, gsi1_25, \
                         gsi1_26, gsk_31, gsk_32, gsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * gsi0_24[k]
                  - f_11 * gsi1_24[k]
                  + f_3 * pc_y[k] * gsk_31[k];

        t_40[k] = f_8 * gsi0_25[k]
                  - f_9 * gsi1_25[k]
                  + f_3 * pc_y[k] * gsk_32[k];

        t_41[k] = f_6 * gsi0_26[k]
                  - f_7 * gsi1_26[k]
                  + f_3 * pc_y[k] * gsk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, fsl0_0, fsk_0, \
                         fsl1_0, gsi0_27, gsi1_27, gsk_34, gsk_35, \
                         gsk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * gsi0_27[k]
                  - f_5 * gsi1_27[k]
                  + f_3 * pc_y[k] * gsk_34[k];

        t_43[k] = f_3 * pc_y[k] * gsk_35[k];

        t_44[k] = f_1 * gsi0_27[k]
                  - f_2 * gsi1_27[k]
                  + f_3 * pc_z[k] * gsk_35[k];

        t_45[k] = pa_y[k] * fsl0_0[k]
                  - f_14 * pc_y[k] * fsl1_0[k];

        t_46[k] = f_15 * fsk_0[k]
                  + f_3 * pc_y[k] * gsk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, fsl0_3, fsl0_5, fsk_1, \
                         fsl1_3, fsl1_5, gsk_36, gsk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * gsk_36[k];

        t_48[k] = pa_y[k] * fsl0_3[k]
                  + f_16 * fsk_1[k]
                  - f_14 * pc_y[k] * fsl1_3[k];

        t_49[k] = f_3 * pc_z[k] * gsk_37[k];

        t_50[k] = pa_y[k] * fsl0_5[k]
                  - f_14 * pc_y[k] * fsl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, fsl0_6, fsl0_9, fsk_3, \
                         fsk_5, fsl1_6, fsl1_9, gsk_39, gsk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * fsl0_6[k]
                  + f_17 * fsk_3[k]
                  - f_14 * pc_y[k] * fsl1_6[k];

        t_52[k] = f_3 * pc_z[k] * gsk_39[k];

        t_53[k] = f_15 * fsk_5[k]
                  + f_3 * pc_y[k] * gsk_41[k];

        t_54[k] = pa_y[k] * fsl0_9[k]
                  - f_14 * pc_y[k] * fsl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, fsl0_10, fsk_6, fsk_9, \
                         fsl1_10, gsi0_31, gsi1_31, gsk_42, gsk_43, \
                         gsk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * fsl0_10[k]
                  + f_0 * fsk_6[k]
                  - f_14 * pc_y[k] * fsl1_10[k];

        t_56[k] = f_3 * pc_z[k] * gsk_42[k];

        t_57[k] = f_4 * gsi0_31[k]
                  - f_5 * gsi1_31[k]
                  + f_3 * pc_z[k] * gsk_43[k];

        t_58[k] = f_15 * fsk_9[k]
                  + f_3 * pc_y[k] * gsk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, fsl0_14, fsl0_15, fsk_10, \
                         fsl1_14, fsl1_15, gsi0_34, gsi1_34, gsk_46, \
                         gsk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * fsl0_14[k]
                  - f_14 * pc_y[k] * fsl1_14[k];

        t_60[k] = pa_y[k] * fsl0_15[k]
                  + f_18 * fsk_10[k]
                  - f_14 * pc_y[k] * fsl1_15[k];

        t_61[k] = f_3 * pc_z[k] * gsk_46[k];

        t_62[k] = f_4 * gsi0_34[k]
                  - f_5 * gsi1_34[k]
                  + f_3 * pc_z[k] * gsk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, fsl0_20, fsk_14, fsl1_20, \
                         gsi0_35, gsi1_35, gsk_48, gsk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * gsi0_35[k]
                  - f_7 * gsi1_35[k]
                  + f_3 * pc_z[k] * gsk_48[k];

        t_64[k] = f_15 * fsk_14[k]
                  + f_3 * pc_y[k] * gsk_50[k];

        t_65[k] = pa_y[k] * fsl0_20[k]
                  - f_14 * pc_y[k] * fsl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, fsl0_21, fsk_15, fsl1_21, \
                         gsi0_38, gsi1_38, gsk_51, gsk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * fsl0_21[k]
                  + f_19 * fsk_15[k]
                  - f_14 * pc_y[k] * fsl1_21[k];

        t_67[k] = f_3 * pc_z[k] * gsk_51[k];

        t_68[k] = f_4 * gsi0_38[k]
                  - f_5 * gsi1_38[k]
                  + f_3 * pc_z[k] * gsk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, fsk_20, gsi0_39, gsi0_40, gsi1_39, \
                         gsi1_40, gsk_53, gsk_54, gsk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * gsi0_39[k]
                  - f_7 * gsi1_39[k]
                  + f_3 * pc_z[k] * gsk_53[k];

        t_70[k] = f_8 * gsi0_40[k]
                  - f_9 * gsi1_40[k]
                  + f_3 * pc_z[k] * gsk_54[k];

        t_71[k] = f_15 * fsk_20[k]
                  + f_3 * pc_y[k] * gsk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, fsl0_27, fsk_64, \
                         fsk_66, fsl1_27, gsk_57, gsk_64, gsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * fsl0_27[k]
                  - f_14 * pc_y[k] * fsl1_27[k];

        t_73[k] = f_17 * fsk_64[k]
                  + f_3 * pc_x[k] * gsk_64[k];

        t_74[k] = f_3 * pc_z[k] * gsk_57[k];

        t_75[k] = f_17 * fsk_66[k]
                  + f_3 * pc_x[k] * gsk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, fsk_67, fsk_68, fsk_69, fsk_70, \
                         fsk_71, gsk_67, gsk_68, gsk_69, gsk_70, \
                         gsk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_17 * fsk_67[k]
                  + f_3 * pc_x[k] * gsk_67[k];

        t_77[k] = f_17 * fsk_68[k]
                  + f_3 * pc_x[k] * gsk_68[k];

        t_78[k] = f_17 * fsk_69[k]
                  + f_3 * pc_x[k] * gsk_69[k];

        t_79[k] = f_17 * fsk_70[k]
                  + f_3 * pc_x[k] * gsk_70[k];

        t_80[k] = f_17 * fsk_71[k]
                  + f_3 * pc_x[k] * gsk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, fsk_28, gsi0_49, gsi0_50, \
                         gsi1_49, gsi1_50, gsk_64, gsk_65, gsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * fsk_28[k]
                  + f_1 * gsi0_49[k]
                  - f_2 * gsi1_49[k]
                  + f_3 * pc_y[k] * gsk_64[k];

        t_82[k] = f_3 * pc_z[k] * gsk_64[k];

        t_83[k] = f_4 * gsi0_49[k]
                  - f_5 * gsi1_49[k]
                  + f_3 * pc_z[k] * gsk_65[k];

        t_84[k] = f_6 * gsi0_50[k]
                  - f_7 * gsi1_50[k]
                  + f_3 * pc_z[k] * gsk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, gsi0_51, gsi0_52, gsi0_53, gsi1_51, gsi1_52, \
                         gsi1_53, gsk_67, gsk_68, gsk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * gsi0_51[k]
                  - f_9 * gsi1_51[k]
                  + f_3 * pc_z[k] * gsk_67[k];

        t_86[k] = f_10 * gsi0_52[k]
                  - f_11 * gsi1_52[k]
                  + f_3 * pc_z[k] * gsk_68[k];

        t_87[k] = f_12 * gsi0_53[k]
                  - f_13 * gsi1_53[k]
                  + f_3 * pc_z[k] * gsk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, fsl0_0, fsl0_44, \
                         fsk_35, fsl1_0, fsl1_44, gsk_71, gsk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * fsk_35[k]
                  + f_3 * pc_y[k] * gsk_71[k];

        t_89[k] = pa_y[k] * fsl0_44[k]
                  - f_14 * pc_y[k] * fsl1_44[k];

        t_90[k] = pa_z[k] * fsl0_0[k]
                  - f_14 * pc_z[k] * fsl1_0[k];

        t_91[k] = f_3 * pc_y[k] * gsk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, fsl0_3, fsl0_5, fsk_0, \
                         fsk_2, fsl1_3, fsl1_5, gsk_72, gsk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * fsk_0[k]
                  + f_3 * pc_z[k] * gsk_72[k];

        t_93[k] = pa_z[k] * fsl0_3[k]
                  - f_14 * pc_z[k] * fsl1_3[k];

        t_94[k] = f_3 * pc_y[k] * gsk_74[k];

        t_95[k] = pa_z[k] * fsl0_5[k]
                  + f_16 * fsk_2[k]
                  - f_14 * pc_z[k] * fsl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, fsl0_6, fsl0_9, fsk_5, \
                         fsl1_6, fsl1_9, gsi0_58, gsi1_58, gsk_76, \
                         gsk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * fsl0_6[k]
                  - f_14 * pc_z[k] * fsl1_6[k];

        t_97[k] = f_4 * gsi0_58[k]
                  - f_5 * gsi1_58[k]
                  + f_3 * pc_y[k] * gsk_76[k];

        t_98[k] = f_3 * pc_y[k] * gsk_77[k];

        t_99[k] = pa_z[k] * fsl0_9[k]
                  + f_17 * fsk_5[k]
                  - f_14 * pc_z[k] * fsl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, fsl0_10, fsl1_10, \
                         gsi0_60, gsi0_61, gsi1_60, gsi1_61, gsk_79, gsk_80, \
                         gsk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * fsl0_10[k]
                   - f_14 * pc_z[k] * fsl1_10[k];

        t_101[k] = f_6 * gsi0_60[k]
                   - f_7 * gsi1_60[k]
                   + f_3 * pc_y[k] * gsk_79[k];

        t_102[k] = f_4 * gsi0_61[k]
                   - f_5 * gsi1_61[k]
                   + f_3 * pc_y[k] * gsk_80[k];

        t_103[k] = f_3 * pc_y[k] * gsk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, fsl0_14, fsl0_15, fsk_9, \
                         fsl1_14, fsl1_15, gsi0_63, gsi1_63, gsk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * fsl0_14[k]
                   + f_0 * fsk_9[k]
                   - f_14 * pc_z[k] * fsl1_14[k];

        t_105[k] = pa_z[k] * fsl0_15[k]
                   - f_14 * pc_z[k] * fsl1_15[k];

        t_106[k] = f_8 * gsi0_63[k]
                   - f_9 * gsi1_63[k]
                   + f_3 * pc_y[k] * gsk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, gsi0_64, gsi0_65, gsi1_64, gsi1_65, \
                         gsk_84, gsk_85, gsk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * gsi0_64[k]
                   - f_7 * gsi1_64[k]
                   + f_3 * pc_y[k] * gsk_84[k];

        t_108[k] = f_4 * gsi0_65[k]
                   - f_5 * gsi1_65[k]
                   + f_3 * pc_y[k] * gsk_85[k];

        t_109[k] = f_3 * pc_y[k] * gsk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, fsl0_20, fsl0_21, fsk_14, \
                         fsl1_20, fsl1_21, gsi0_67, gsi1_67, gsk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * fsl0_20[k]
                   + f_18 * fsk_14[k]
                   - f_14 * pc_z[k] * fsl1_20[k];

        t_111[k] = pa_z[k] * fsl0_21[k]
                   - f_14 * pc_z[k] * fsl1_21[k];

        t_112[k] = f_10 * gsi0_67[k]
                   - f_11 * gsi1_67[k]
                   + f_3 * pc_y[k] * gsk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, gsi0_68, gsi0_69, gsi0_70, gsi1_68, \
                         gsi1_69, gsi1_70, gsk_89, gsk_90, gsk_91, \
                         gsk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * gsi0_68[k]
                   - f_9 * gsi1_68[k]
                   + f_3 * pc_y[k] * gsk_89[k];

        t_114[k] = f_6 * gsi0_69[k]
                   - f_7 * gsi1_69[k]
                   + f_3 * pc_y[k] * gsk_90[k];

        t_115[k] = f_4 * gsi0_70[k]
                   - f_5 * gsi1_70[k]
                   + f_3 * pc_y[k] * gsk_91[k];

        t_116[k] = f_3 * pc_y[k] * gsk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, fsl0_27, fsk_20, \
                         fsk_100, fsk_101, fsk_102, fsl1_27, gsk_100, gsk_101, \
                         gsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * fsl0_27[k]
                   + f_19 * fsk_20[k]
                   - f_14 * pc_z[k] * fsl1_27[k];

        t_118[k] = f_17 * fsk_100[k]
                   + f_3 * pc_x[k] * gsk_100[k];

        t_119[k] = f_17 * fsk_101[k]
                   + f_3 * pc_x[k] * gsk_101[k];

        t_120[k] = f_17 * fsk_102[k]
                   + f_3 * pc_x[k] * gsk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, fsk_103, fsk_104, \
                         fsk_105, fsk_107, gsk_99, gsk_103, gsk_104, gsk_105, \
                         gsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_17 * fsk_103[k]
                   + f_3 * pc_x[k] * gsk_103[k];

        t_122[k] = f_17 * fsk_104[k]
                   + f_3 * pc_x[k] * gsk_104[k];

        t_123[k] = f_17 * fsk_105[k]
                   + f_3 * pc_x[k] * gsk_105[k];

        t_124[k] = f_3 * pc_y[k] * gsk_99[k];

        t_125[k] = f_17 * fsk_107[k]
                   + f_3 * pc_x[k] * gsk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, fsl0_36, fsl1_36, gsi0_78, \
                         gsi0_79, gsi1_78, gsi1_79, gsk_101, gsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * fsl0_36[k]
                   - f_14 * pc_z[k] * fsl1_36[k];

        t_127[k] = f_20 * gsi0_78[k]
                   - f_21 * gsi1_78[k]
                   + f_3 * pc_y[k] * gsk_101[k];

        t_128[k] = f_12 * gsi0_79[k]
                   - f_13 * gsi1_79[k]
                   + f_3 * pc_y[k] * gsk_102[k];
    }
}

static auto
compute_prim_gsl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsl0,
                                                          const size_t fsk, const size_t fsl1,
                                                          const size_t gsi0, const size_t gsi1,
                                                          const size_t gsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsl0_48 = buffer.data(fsl0 + 48);
    const auto *fsl0_51 = buffer.data(fsl0 + 51);
    const auto *fsl0_55 = buffer.data(fsl0 + 55);
    const auto *fsl0_60 = buffer.data(fsl0 + 60);
    const auto *fsl0_66 = buffer.data(fsl0 + 66);
    const auto *fsl0_81 = buffer.data(fsl0 + 81);
    const auto *fsl0_90 = buffer.data(fsl0 + 90);
    const auto *fsl0_95 = buffer.data(fsl0 + 95);
    const auto *fsl0_99 = buffer.data(fsl0 + 99);
    const auto *fsl0_102 = buffer.data(fsl0 + 102);
    const auto *fsl0_104 = buffer.data(fsl0 + 104);
    const auto *fsl0_107 = buffer.data(fsl0 + 107);
    const auto *fsl0_108 = buffer.data(fsl0 + 108);
    const auto *fsl0_110 = buffer.data(fsl0 + 110);
    const auto *fsl0_113 = buffer.data(fsl0 + 113);
    const auto *fsl0_114 = buffer.data(fsl0 + 114);
    const auto *fsl0_115 = buffer.data(fsl0 + 115);
    const auto *fsl0_117 = buffer.data(fsl0 + 117);
    const auto *fsl0_134 = buffer.data(fsl0 + 134);

    const auto *fsk_35 = buffer.data(fsk + 35);
    const auto *fsk_36 = buffer.data(fsk + 36);
    const auto *fsk_39 = buffer.data(fsk + 39);
    const auto *fsk_41 = buffer.data(fsk + 41);
    const auto *fsk_42 = buffer.data(fsk + 42);
    const auto *fsk_45 = buffer.data(fsk + 45);
    const auto *fsk_46 = buffer.data(fsk + 46);
    const auto *fsk_50 = buffer.data(fsk + 50);
    const auto *fsk_51 = buffer.data(fsk + 51);
    const auto *fsk_56 = buffer.data(fsk + 56);
    const auto *fsk_64 = buffer.data(fsk + 64);
    const auto *fsk_71 = buffer.data(fsk + 71);
    const auto *fsk_72 = buffer.data(fsk + 72);
    const auto *fsk_74 = buffer.data(fsk + 74);
    const auto *fsk_77 = buffer.data(fsk + 77);
    const auto *fsk_80 = buffer.data(fsk + 80);
    const auto *fsk_81 = buffer.data(fsk + 81);
    const auto *fsk_84 = buffer.data(fsk + 84);
    const auto *fsk_85 = buffer.data(fsk + 85);
    const auto *fsk_86 = buffer.data(fsk + 86);
    const auto *fsk_89 = buffer.data(fsk + 89);
    const auto *fsk_90 = buffer.data(fsk + 90);
    const auto *fsk_91 = buffer.data(fsk + 91);
    const auto *fsk_92 = buffer.data(fsk + 92);
    const auto *fsk_102 = buffer.data(fsk + 102);
    const auto *fsk_103 = buffer.data(fsk + 103);
    const auto *fsk_104 = buffer.data(fsk + 104);
    const auto *fsk_105 = buffer.data(fsk + 105);
    const auto *fsk_106 = buffer.data(fsk + 106);
    const auto *fsk_107 = buffer.data(fsk + 107);
    const auto *fsk_108 = buffer.data(fsk + 108);
    const auto *fsk_111 = buffer.data(fsk + 111);
    const auto *fsk_114 = buffer.data(fsk + 114);
    const auto *fsk_118 = buffer.data(fsk + 118);
    const auto *fsk_123 = buffer.data(fsk + 123);
    const auto *fsk_129 = buffer.data(fsk + 129);
    const auto *fsk_136 = buffer.data(fsk + 136);
    const auto *fsk_138 = buffer.data(fsk + 138);
    const auto *fsk_139 = buffer.data(fsk + 139);
    const auto *fsk_140 = buffer.data(fsk + 140);
    const auto *fsk_141 = buffer.data(fsk + 141);
    const auto *fsk_142 = buffer.data(fsk + 142);
    const auto *fsk_143 = buffer.data(fsk + 143);
    const auto *fsk_172 = buffer.data(fsk + 172);
    const auto *fsk_173 = buffer.data(fsk + 173);
    const auto *fsk_174 = buffer.data(fsk + 174);
    const auto *fsk_175 = buffer.data(fsk + 175);
    const auto *fsk_176 = buffer.data(fsk + 176);
    const auto *fsk_177 = buffer.data(fsk + 177);
    const auto *fsk_178 = buffer.data(fsk + 178);
    const auto *fsk_179 = buffer.data(fsk + 179);
    const auto *fsk_180 = buffer.data(fsk + 180);
    const auto *fsk_185 = buffer.data(fsk + 185);
    const auto *fsk_189 = buffer.data(fsk + 189);
    const auto *fsk_194 = buffer.data(fsk + 194);
    const auto *fsk_200 = buffer.data(fsk + 200);

    const auto *fsl1_48 = buffer.data(fsl1 + 48);
    const auto *fsl1_51 = buffer.data(fsl1 + 51);
    const auto *fsl1_55 = buffer.data(fsl1 + 55);
    const auto *fsl1_60 = buffer.data(fsl1 + 60);
    const auto *fsl1_66 = buffer.data(fsl1 + 66);
    const auto *fsl1_81 = buffer.data(fsl1 + 81);
    const auto *fsl1_90 = buffer.data(fsl1 + 90);
    const auto *fsl1_95 = buffer.data(fsl1 + 95);
    const auto *fsl1_99 = buffer.data(fsl1 + 99);
    const auto *fsl1_102 = buffer.data(fsl1 + 102);
    const auto *fsl1_104 = buffer.data(fsl1 + 104);
    const auto *fsl1_107 = buffer.data(fsl1 + 107);
    const auto *fsl1_108 = buffer.data(fsl1 + 108);
    const auto *fsl1_110 = buffer.data(fsl1 + 110);
    const auto *fsl1_113 = buffer.data(fsl1 + 113);
    const auto *fsl1_114 = buffer.data(fsl1 + 114);
    const auto *fsl1_115 = buffer.data(fsl1 + 115);
    const auto *fsl1_117 = buffer.data(fsl1 + 117);
    const auto *fsl1_134 = buffer.data(fsl1 + 134);

    const auto *gsi0_80 = buffer.data(gsi0 + 80);
    const auto *gsi0_81 = buffer.data(gsi0 + 81);
    const auto *gsi0_82 = buffer.data(gsi0 + 82);
    const auto *gsi0_83 = buffer.data(gsi0 + 83);
    const auto *gsi0_84 = buffer.data(gsi0 + 84);
    const auto *gsi0_86 = buffer.data(gsi0 + 86);
    const auto *gsi0_87 = buffer.data(gsi0 + 87);
    const auto *gsi0_89 = buffer.data(gsi0 + 89);
    const auto *gsi0_90 = buffer.data(gsi0 + 90);
    const auto *gsi0_91 = buffer.data(gsi0 + 91);
    const auto *gsi0_93 = buffer.data(gsi0 + 93);
    const auto *gsi0_94 = buffer.data(gsi0 + 94);
    const auto *gsi0_95 = buffer.data(gsi0 + 95);
    const auto *gsi0_96 = buffer.data(gsi0 + 96);
    const auto *gsi0_98 = buffer.data(gsi0 + 98);
    const auto *gsi0_99 = buffer.data(gsi0 + 99);
    const auto *gsi0_105 = buffer.data(gsi0 + 105);
    const auto *gsi0_106 = buffer.data(gsi0 + 106);
    const auto *gsi0_107 = buffer.data(gsi0 + 107);
    const auto *gsi0_108 = buffer.data(gsi0 + 108);
    const auto *gsi0_109 = buffer.data(gsi0 + 109);
    const auto *gsi0_111 = buffer.data(gsi0 + 111);
    const auto *gsi0_135 = buffer.data(gsi0 + 135);
    const auto *gsi0_136 = buffer.data(gsi0 + 136);
    const auto *gsi0_137 = buffer.data(gsi0 + 137);
    const auto *gsi0_138 = buffer.data(gsi0 + 138);
    const auto *gsi0_139 = buffer.data(gsi0 + 139);
    const auto *gsi0_140 = buffer.data(gsi0 + 140);
    const auto *gsi0_141 = buffer.data(gsi0 + 141);
    const auto *gsi0_142 = buffer.data(gsi0 + 142);
    const auto *gsi0_143 = buffer.data(gsi0 + 143);
    const auto *gsi0_144 = buffer.data(gsi0 + 144);
    const auto *gsi0_145 = buffer.data(gsi0 + 145);
    const auto *gsi0_146 = buffer.data(gsi0 + 146);
    const auto *gsi0_147 = buffer.data(gsi0 + 147);
    const auto *gsi0_148 = buffer.data(gsi0 + 148);
    const auto *gsi0_149 = buffer.data(gsi0 + 149);
    const auto *gsi0_154 = buffer.data(gsi0 + 154);
    const auto *gsi0_160 = buffer.data(gsi0 + 160);

    const auto *gsi1_80 = buffer.data(gsi1 + 80);
    const auto *gsi1_81 = buffer.data(gsi1 + 81);
    const auto *gsi1_82 = buffer.data(gsi1 + 82);
    const auto *gsi1_83 = buffer.data(gsi1 + 83);
    const auto *gsi1_84 = buffer.data(gsi1 + 84);
    const auto *gsi1_86 = buffer.data(gsi1 + 86);
    const auto *gsi1_87 = buffer.data(gsi1 + 87);
    const auto *gsi1_89 = buffer.data(gsi1 + 89);
    const auto *gsi1_90 = buffer.data(gsi1 + 90);
    const auto *gsi1_91 = buffer.data(gsi1 + 91);
    const auto *gsi1_93 = buffer.data(gsi1 + 93);
    const auto *gsi1_94 = buffer.data(gsi1 + 94);
    const auto *gsi1_95 = buffer.data(gsi1 + 95);
    const auto *gsi1_96 = buffer.data(gsi1 + 96);
    const auto *gsi1_98 = buffer.data(gsi1 + 98);
    const auto *gsi1_99 = buffer.data(gsi1 + 99);
    const auto *gsi1_105 = buffer.data(gsi1 + 105);
    const auto *gsi1_106 = buffer.data(gsi1 + 106);
    const auto *gsi1_107 = buffer.data(gsi1 + 107);
    const auto *gsi1_108 = buffer.data(gsi1 + 108);
    const auto *gsi1_109 = buffer.data(gsi1 + 109);
    const auto *gsi1_111 = buffer.data(gsi1 + 111);
    const auto *gsi1_135 = buffer.data(gsi1 + 135);
    const auto *gsi1_136 = buffer.data(gsi1 + 136);
    const auto *gsi1_137 = buffer.data(gsi1 + 137);
    const auto *gsi1_138 = buffer.data(gsi1 + 138);
    const auto *gsi1_139 = buffer.data(gsi1 + 139);
    const auto *gsi1_140 = buffer.data(gsi1 + 140);
    const auto *gsi1_141 = buffer.data(gsi1 + 141);
    const auto *gsi1_142 = buffer.data(gsi1 + 142);
    const auto *gsi1_143 = buffer.data(gsi1 + 143);
    const auto *gsi1_144 = buffer.data(gsi1 + 144);
    const auto *gsi1_145 = buffer.data(gsi1 + 145);
    const auto *gsi1_146 = buffer.data(gsi1 + 146);
    const auto *gsi1_147 = buffer.data(gsi1 + 147);
    const auto *gsi1_148 = buffer.data(gsi1 + 148);
    const auto *gsi1_149 = buffer.data(gsi1 + 149);
    const auto *gsi1_154 = buffer.data(gsi1 + 154);
    const auto *gsi1_160 = buffer.data(gsi1 + 160);

    const auto *gsk_103 = buffer.data(gsk + 103);
    const auto *gsk_104 = buffer.data(gsk + 104);
    const auto *gsk_105 = buffer.data(gsk + 105);
    const auto *gsk_106 = buffer.data(gsk + 106);
    const auto *gsk_107 = buffer.data(gsk + 107);
    const auto *gsk_108 = buffer.data(gsk + 108);
    const auto *gsk_109 = buffer.data(gsk + 109);
    const auto *gsk_110 = buffer.data(gsk + 110);
    const auto *gsk_111 = buffer.data(gsk + 111);
    const auto *gsk_113 = buffer.data(gsk + 113);
    const auto *gsk_114 = buffer.data(gsk + 114);
    const auto *gsk_115 = buffer.data(gsk + 115);
    const auto *gsk_117 = buffer.data(gsk + 117);
    const auto *gsk_118 = buffer.data(gsk + 118);
    const auto *gsk_119 = buffer.data(gsk + 119);
    const auto *gsk_120 = buffer.data(gsk + 120);
    const auto *gsk_122 = buffer.data(gsk + 122);
    const auto *gsk_123 = buffer.data(gsk + 123);
    const auto *gsk_124 = buffer.data(gsk + 124);
    const auto *gsk_125 = buffer.data(gsk + 125);
    const auto *gsk_126 = buffer.data(gsk + 126);
    const auto *gsk_128 = buffer.data(gsk + 128);
    const auto *gsk_129 = buffer.data(gsk + 129);
    const auto *gsk_136 = buffer.data(gsk + 136);
    const auto *gsk_137 = buffer.data(gsk + 137);
    const auto *gsk_138 = buffer.data(gsk + 138);
    const auto *gsk_139 = buffer.data(gsk + 139);
    const auto *gsk_140 = buffer.data(gsk + 140);
    const auto *gsk_141 = buffer.data(gsk + 141);
    const auto *gsk_142 = buffer.data(gsk + 142);
    const auto *gsk_143 = buffer.data(gsk + 143);
    const auto *gsk_144 = buffer.data(gsk + 144);
    const auto *gsk_146 = buffer.data(gsk + 146);
    const auto *gsk_147 = buffer.data(gsk + 147);
    const auto *gsk_149 = buffer.data(gsk + 149);
    const auto *gsk_150 = buffer.data(gsk + 150);
    const auto *gsk_153 = buffer.data(gsk + 153);
    const auto *gsk_154 = buffer.data(gsk + 154);
    const auto *gsk_158 = buffer.data(gsk + 158);
    const auto *gsk_159 = buffer.data(gsk + 159);
    const auto *gsk_164 = buffer.data(gsk + 164);
    const auto *gsk_172 = buffer.data(gsk + 172);
    const auto *gsk_173 = buffer.data(gsk + 173);
    const auto *gsk_174 = buffer.data(gsk + 174);
    const auto *gsk_175 = buffer.data(gsk + 175);
    const auto *gsk_176 = buffer.data(gsk + 176);
    const auto *gsk_177 = buffer.data(gsk + 177);
    const auto *gsk_178 = buffer.data(gsk + 178);
    const auto *gsk_179 = buffer.data(gsk + 179);
    const auto *gsk_180 = buffer.data(gsk + 180);
    const auto *gsk_181 = buffer.data(gsk + 181);
    const auto *gsk_182 = buffer.data(gsk + 182);
    const auto *gsk_183 = buffer.data(gsk + 183);
    const auto *gsk_184 = buffer.data(gsk + 184);
    const auto *gsk_185 = buffer.data(gsk + 185);
    const auto *gsk_186 = buffer.data(gsk + 186);
    const auto *gsk_187 = buffer.data(gsk + 187);
    const auto *gsk_188 = buffer.data(gsk + 188);
    const auto *gsk_189 = buffer.data(gsk + 189);
    const auto *gsk_190 = buffer.data(gsk + 190);
    const auto *gsk_191 = buffer.data(gsk + 191);
    const auto *gsk_192 = buffer.data(gsk + 192);
    const auto *gsk_193 = buffer.data(gsk + 193);
    const auto *gsk_194 = buffer.data(gsk + 194);
    const auto *gsk_200 = buffer.data(gsk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, gsi0_80, gsi0_81, gsi0_82, gsi1_80, \
                         gsi1_81, gsi1_82, gsk_103, gsk_104, gsk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * gsi0_80[k]
                   - f_11 * gsi1_80[k]
                   + f_3 * pc_y[k] * gsk_103[k];

        t_130[k] = f_8 * gsi0_81[k]
                   - f_9 * gsi1_81[k]
                   + f_3 * pc_y[k] * gsk_104[k];

        t_131[k] = f_6 * gsi0_82[k]
                   - f_7 * gsi1_82[k]
                   + f_3 * pc_y[k] * gsk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, fsk_35, fsk_108, \
                         gsi0_83, gsi0_84, gsi1_83, gsi1_84, gsk_106, gsk_107, \
                         gsk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * gsi0_83[k]
                   - f_5 * gsi1_83[k]
                   + f_3 * pc_y[k] * gsk_106[k];

        t_133[k] = f_3 * pc_y[k] * gsk_107[k];

        t_134[k] = f_15 * fsk_35[k]
                   + f_1 * gsi0_83[k]
                   - f_2 * gsi1_83[k]
                   + f_3 * pc_z[k] * gsk_107[k];

        t_135[k] = f_16 * fsk_108[k]
                   + f_1 * gsi0_84[k]
                   - f_2 * gsi1_84[k]
                   + f_3 * pc_x[k] * gsk_108[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, fsk_36, fsk_111, \
                         gsi0_87, gsi1_87, gsk_108, gsk_109, gsk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * fsk_36[k]
                   + f_3 * pc_y[k] * gsk_108[k];

        t_137[k] = f_3 * pc_z[k] * gsk_108[k];

        t_138[k] = f_16 * fsk_111[k]
                   + f_12 * gsi0_87[k]
                   - f_13 * gsi1_87[k]
                   + f_3 * pc_x[k] * gsk_111[k];

        t_139[k] = f_3 * pc_z[k] * gsk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_z, fsk_114, gsi0_84, gsi0_90, gsi1_84, \
                         gsi1_90, gsk_110, gsk_111, gsk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * gsi0_84[k]
                   - f_5 * gsi1_84[k]
                   + f_3 * pc_z[k] * gsk_110[k];

        t_141[k] = f_16 * fsk_114[k]
                   + f_10 * gsi0_90[k]
                   - f_11 * gsi1_90[k]
                   + f_3 * pc_x[k] * gsk_114[k];

        t_142[k] = f_3 * pc_z[k] * gsk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, fsk_41, fsk_118, \
                         gsi0_86, gsi0_94, gsi1_86, gsi1_94, gsk_113, gsk_114, \
                         gsk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * fsk_41[k]
                   + f_3 * pc_y[k] * gsk_113[k];

        t_144[k] = f_6 * gsi0_86[k]
                   - f_7 * gsi1_86[k]
                   + f_3 * pc_z[k] * gsk_113[k];

        t_145[k] = f_16 * fsk_118[k]
                   + f_8 * gsi0_94[k]
                   - f_9 * gsi1_94[k]
                   + f_3 * pc_x[k] * gsk_118[k];

        t_146[k] = f_3 * pc_z[k] * gsk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, fsk_45, gsi0_87, gsi0_89, gsi1_87, \
                         gsi1_89, gsk_115, gsk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * gsi0_87[k]
                   - f_5 * gsi1_87[k]
                   + f_3 * pc_z[k] * gsk_115[k];

        t_148[k] = f_16 * fsk_45[k]
                   + f_3 * pc_y[k] * gsk_117[k];

        t_149[k] = f_8 * gsi0_89[k]
                   - f_9 * gsi1_89[k]
                   + f_3 * pc_z[k] * gsk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, fsk_123, gsi0_90, gsi0_99, gsi1_90, \
                         gsi1_99, gsk_118, gsk_119, gsk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_16 * fsk_123[k]
                   + f_6 * gsi0_99[k]
                   - f_7 * gsi1_99[k]
                   + f_3 * pc_x[k] * gsk_123[k];

        t_151[k] = f_3 * pc_z[k] * gsk_118[k];

        t_152[k] = f_4 * gsi0_90[k]
                   - f_5 * gsi1_90[k]
                   + f_3 * pc_z[k] * gsk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, fsk_50, gsi0_91, gsi0_93, gsi1_91, \
                         gsi1_93, gsk_120, gsk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * gsi0_91[k]
                   - f_7 * gsi1_91[k]
                   + f_3 * pc_z[k] * gsk_120[k];

        t_154[k] = f_16 * fsk_50[k]
                   + f_3 * pc_y[k] * gsk_122[k];

        t_155[k] = f_10 * gsi0_93[k]
                   - f_11 * gsi1_93[k]
                   + f_3 * pc_z[k] * gsk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, fsk_129, gsi0_94, gsi0_105, gsi1_94, \
                         gsi1_105, gsk_123, gsk_124, gsk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * fsk_129[k]
                   + f_4 * gsi0_105[k]
                   - f_5 * gsi1_105[k]
                   + f_3 * pc_x[k] * gsk_129[k];

        t_157[k] = f_3 * pc_z[k] * gsk_123[k];

        t_158[k] = f_4 * gsi0_94[k]
                   - f_5 * gsi1_94[k]
                   + f_3 * pc_z[k] * gsk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, fsk_56, gsi0_95, gsi0_96, \
                         gsi0_98, gsi1_95, gsi1_96, gsi1_98, gsk_125, gsk_126, \
                         gsk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * gsi0_95[k]
                   - f_7 * gsi1_95[k]
                   + f_3 * pc_z[k] * gsk_125[k];

        t_160[k] = f_8 * gsi0_96[k]
                   - f_9 * gsi1_96[k]
                   + f_3 * pc_z[k] * gsk_126[k];

        t_161[k] = f_16 * fsk_56[k]
                   + f_3 * pc_y[k] * gsk_128[k];

        t_162[k] = f_12 * gsi0_98[k]
                   - f_13 * gsi1_98[k]
                   + f_3 * pc_z[k] * gsk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, fsk_136, fsk_138, \
                         fsk_139, fsk_140, gsk_129, gsk_136, gsk_138, gsk_139, \
                         gsk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_16 * fsk_136[k]
                   + f_3 * pc_x[k] * gsk_136[k];

        t_164[k] = f_3 * pc_z[k] * gsk_129[k];

        t_165[k] = f_16 * fsk_138[k]
                   + f_3 * pc_x[k] * gsk_138[k];

        t_166[k] = f_16 * fsk_139[k]
                   + f_3 * pc_x[k] * gsk_139[k];

        t_167[k] = f_16 * fsk_140[k]
                   + f_3 * pc_x[k] * gsk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, fsk_64, fsk_141, fsk_142, \
                         fsk_143, gsi0_105, gsi1_105, gsk_136, gsk_141, gsk_142, \
                         gsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_16 * fsk_141[k]
                   + f_3 * pc_x[k] * gsk_141[k];

        t_169[k] = f_16 * fsk_142[k]
                   + f_3 * pc_x[k] * gsk_142[k];

        t_170[k] = f_16 * fsk_143[k]
                   + f_3 * pc_x[k] * gsk_143[k];

        t_171[k] = f_16 * fsk_64[k]
                   + f_1 * gsi0_105[k]
                   - f_2 * gsi1_105[k]
                   + f_3 * pc_y[k] * gsk_136[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_z, gsi0_105, gsi0_106, gsi0_107, \
                         gsi1_105, gsi1_106, gsi1_107, gsk_136, gsk_137, gsk_138, \
                         gsk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * gsk_136[k];

        t_173[k] = f_4 * gsi0_105[k]
                   - f_5 * gsi1_105[k]
                   + f_3 * pc_z[k] * gsk_137[k];

        t_174[k] = f_6 * gsi0_106[k]
                   - f_7 * gsi1_106[k]
                   + f_3 * pc_z[k] * gsk_138[k];

        t_175[k] = f_8 * gsi0_107[k]
                   - f_9 * gsi1_107[k]
                   + f_3 * pc_z[k] * gsk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, fsk_71, gsi0_108, gsi0_109, \
                         gsi0_111, gsi1_108, gsi1_109, gsi1_111, gsk_140, gsk_141, \
                         gsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * gsi0_108[k]
                   - f_11 * gsi1_108[k]
                   + f_3 * pc_z[k] * gsk_140[k];

        t_177[k] = f_12 * gsi0_109[k]
                   - f_13 * gsi1_109[k]
                   + f_3 * pc_z[k] * gsk_141[k];

        t_178[k] = f_16 * fsk_71[k]
                   + f_3 * pc_y[k] * gsk_143[k];

        t_179[k] = f_1 * gsi0_111[k]
                   - f_2 * gsi1_111[k]
                   + f_3 * pc_z[k] * gsk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, fsl0_48, fsl0_90, \
                         fsk_36, fsk_72, fsl1_48, fsl1_90, gsk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * fsl0_90[k]
                   - f_14 * pc_y[k] * fsl1_90[k];

        t_181[k] = f_15 * fsk_72[k]
                   + f_3 * pc_y[k] * gsk_144[k];

        t_182[k] = f_15 * fsk_36[k]
                   + f_3 * pc_z[k] * gsk_144[k];

        t_183[k] = pa_z[k] * fsl0_48[k]
                   - f_14 * pc_z[k] * fsl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, fsl0_51, fsl0_95, \
                         fsk_39, fsk_74, fsl1_51, fsl1_95, gsk_146, \
                         gsk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * fsk_74[k]
                   + f_3 * pc_y[k] * gsk_146[k];

        t_185[k] = pa_y[k] * fsl0_95[k]
                   - f_14 * pc_y[k] * fsl1_95[k];

        t_186[k] = pa_z[k] * fsl0_51[k]
                   - f_14 * pc_z[k] * fsl1_51[k];

        t_187[k] = f_15 * fsk_39[k]
                   + f_3 * pc_z[k] * gsk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, fsl0_55, fsl0_99, \
                         fsk_42, fsk_77, fsl1_55, fsl1_99, gsk_149, \
                         gsk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * fsk_77[k]
                   + f_3 * pc_y[k] * gsk_149[k];

        t_189[k] = pa_y[k] * fsl0_99[k]
                   - f_14 * pc_y[k] * fsl1_99[k];

        t_190[k] = pa_z[k] * fsl0_55[k]
                   - f_14 * pc_z[k] * fsl1_55[k];

        t_191[k] = f_15 * fsk_42[k]
                   + f_3 * pc_z[k] * gsk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, fsl0_102, fsl0_104, fsk_80, fsk_81, \
                         fsl1_102, fsl1_104, gsk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * fsl0_102[k]
                   + f_16 * fsk_80[k]
                   - f_14 * pc_y[k] * fsl1_102[k];

        t_193[k] = f_15 * fsk_81[k]
                   + f_3 * pc_y[k] * gsk_153[k];

        t_194[k] = pa_y[k] * fsl0_104[k]
                   - f_14 * pc_y[k] * fsl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pa_z, pc_y, pc_z, fsl0_60, fsl0_107, \
                         fsk_46, fsk_84, fsl1_60, fsl1_107, gsk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * fsl0_60[k]
                   - f_14 * pc_z[k] * fsl1_60[k];

        t_196[k] = f_15 * fsk_46[k]
                   + f_3 * pc_z[k] * gsk_154[k];

        t_197[k] = pa_y[k] * fsl0_107[k]
                   + f_17 * fsk_84[k]
                   - f_14 * pc_y[k] * fsl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pc_y, fsl0_108, fsl0_110, fsk_85, fsk_86, \
                         fsl1_108, fsl1_110, gsk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_y[k] * fsl0_108[k]
                   + f_16 * fsk_85[k]
                   - f_14 * pc_y[k] * fsl1_108[k];

        t_199[k] = f_15 * fsk_86[k]
                   + f_3 * pc_y[k] * gsk_158[k];

        t_200[k] = pa_y[k] * fsl0_110[k]
                   - f_14 * pc_y[k] * fsl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pc_y, pc_z, fsl0_66, fsl0_113, \
                         fsk_51, fsk_89, fsl1_66, fsl1_113, gsk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * fsl0_66[k]
                   - f_14 * pc_z[k] * fsl1_66[k];

        t_202[k] = f_15 * fsk_51[k]
                   + f_3 * pc_z[k] * gsk_159[k];

        t_203[k] = pa_y[k] * fsl0_113[k]
                   + f_0 * fsk_89[k]
                   - f_14 * pc_y[k] * fsl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_y, fsl0_114, fsl0_115, fsl0_117, \
                         fsk_90, fsk_91, fsk_92, fsl1_114, fsl1_115, fsl1_117, \
                         gsk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * fsl0_114[k]
                   + f_17 * fsk_90[k]
                   - f_14 * pc_y[k] * fsl1_114[k];

        t_205[k] = pa_y[k] * fsl0_115[k]
                   + f_16 * fsk_91[k]
                   - f_14 * pc_y[k] * fsl1_115[k];

        t_206[k] = f_15 * fsk_92[k]
                   + f_3 * pc_y[k] * gsk_164[k];

        t_207[k] = pa_y[k] * fsl0_117[k]
                   - f_14 * pc_y[k] * fsl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, fsk_172, fsk_173, fsk_174, \
                         fsk_175, fsk_176, gsk_172, gsk_173, gsk_174, gsk_175, \
                         gsk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_16 * fsk_172[k]
                   + f_3 * pc_x[k] * gsk_172[k];

        t_209[k] = f_16 * fsk_173[k]
                   + f_3 * pc_x[k] * gsk_173[k];

        t_210[k] = f_16 * fsk_174[k]
                   + f_3 * pc_x[k] * gsk_174[k];

        t_211[k] = f_16 * fsk_175[k]
                   + f_3 * pc_x[k] * gsk_175[k];

        t_212[k] = f_16 * fsk_176[k]
                   + f_3 * pc_x[k] * gsk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, fsl0_81, fsk_177, \
                         fsk_178, fsk_179, fsl1_81, gsk_177, gsk_178, \
                         gsk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_16 * fsk_177[k]
                   + f_3 * pc_x[k] * gsk_177[k];

        t_214[k] = f_16 * fsk_178[k]
                   + f_3 * pc_x[k] * gsk_178[k];

        t_215[k] = f_16 * fsk_179[k]
                   + f_3 * pc_x[k] * gsk_179[k];

        t_216[k] = pa_z[k] * fsl0_81[k]
                   - f_14 * pc_z[k] * fsl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, fsk_64, fsk_102, fsk_103, gsi0_135, \
                         gsi0_136, gsi1_135, gsi1_136, gsk_172, gsk_174, \
                         gsk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * fsk_64[k]
                   + f_3 * pc_z[k] * gsk_172[k];

        t_218[k] = f_15 * fsk_102[k]
                   + f_12 * gsi0_135[k]
                   - f_13 * gsi1_135[k]
                   + f_3 * pc_y[k] * gsk_174[k];

        t_219[k] = f_15 * fsk_103[k]
                   + f_10 * gsi0_136[k]
                   - f_11 * gsi1_136[k]
                   + f_3 * pc_y[k] * gsk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, fsk_104, fsk_105, fsk_106, gsi0_137, \
                         gsi0_138, gsi0_139, gsi1_137, gsi1_138, gsi1_139, gsk_176, gsk_177, \
                         gsk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * fsk_104[k]
                   + f_8 * gsi0_137[k]
                   - f_9 * gsi1_137[k]
                   + f_3 * pc_y[k] * gsk_176[k];

        t_221[k] = f_15 * fsk_105[k]
                   + f_6 * gsi0_138[k]
                   - f_7 * gsi1_138[k]
                   + f_3 * pc_y[k] * gsk_177[k];

        t_222[k] = f_15 * fsk_106[k]
                   + f_4 * gsi0_139[k]
                   - f_5 * gsi1_139[k]
                   + f_3 * pc_y[k] * gsk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, fsl0_134, fsk_107, \
                         fsk_180, fsl1_134, gsi0_140, gsi1_140, gsk_179, \
                         gsk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * fsk_107[k]
                   + f_3 * pc_y[k] * gsk_179[k];

        t_224[k] = pa_y[k] * fsl0_134[k]
                   - f_14 * pc_y[k] * fsl1_134[k];

        t_225[k] = f_16 * fsk_180[k]
                   + f_1 * gsi0_140[k]
                   - f_2 * gsi1_140[k]
                   + f_3 * pc_x[k] * gsk_180[k];

        t_226[k] = f_3 * pc_y[k] * gsk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, fsk_72, gsi0_140, gsi1_140, gsk_180, \
                         gsk_181, gsk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * fsk_72[k]
                   + f_3 * pc_z[k] * gsk_180[k];

        t_228[k] = f_4 * gsi0_140[k]
                   - f_5 * gsi1_140[k]
                   + f_3 * pc_y[k] * gsk_181[k];

        t_229[k] = f_3 * pc_y[k] * gsk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, fsk_185, gsi0_141, gsi0_142, \
                         gsi0_145, gsi1_141, gsi1_142, gsi1_145, gsk_183, gsk_184, \
                         gsk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_16 * fsk_185[k]
                   + f_12 * gsi0_145[k]
                   - f_13 * gsi1_145[k]
                   + f_3 * pc_x[k] * gsk_185[k];

        t_231[k] = f_6 * gsi0_141[k]
                   - f_7 * gsi1_141[k]
                   + f_3 * pc_y[k] * gsk_183[k];

        t_232[k] = f_4 * gsi0_142[k]
                   - f_5 * gsi1_142[k]
                   + f_3 * pc_y[k] * gsk_184[k];

        t_233[k] = f_3 * pc_y[k] * gsk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, fsk_189, gsi0_143, gsi0_144, \
                         gsi0_149, gsi1_143, gsi1_144, gsi1_149, gsk_186, gsk_187, \
                         gsk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_16 * fsk_189[k]
                   + f_10 * gsi0_149[k]
                   - f_11 * gsi1_149[k]
                   + f_3 * pc_x[k] * gsk_189[k];

        t_235[k] = f_8 * gsi0_143[k]
                   - f_9 * gsi1_143[k]
                   + f_3 * pc_y[k] * gsk_186[k];

        t_236[k] = f_6 * gsi0_144[k]
                   - f_7 * gsi1_144[k]
                   + f_3 * pc_y[k] * gsk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_y, fsk_194, gsi0_145, gsi0_154, \
                         gsi1_145, gsi1_154, gsk_188, gsk_189, \
                         gsk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * gsi0_145[k]
                   - f_5 * gsi1_145[k]
                   + f_3 * pc_y[k] * gsk_188[k];

        t_238[k] = f_3 * pc_y[k] * gsk_189[k];

        t_239[k] = f_16 * fsk_194[k]
                   + f_8 * gsi0_154[k]
                   - f_9 * gsi1_154[k]
                   + f_3 * pc_x[k] * gsk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, gsi0_146, gsi0_147, gsi0_148, gsi1_146, \
                         gsi1_147, gsi1_148, gsk_190, gsk_191, \
                         gsk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * gsi0_146[k]
                   - f_11 * gsi1_146[k]
                   + f_3 * pc_y[k] * gsk_190[k];

        t_241[k] = f_8 * gsi0_147[k]
                   - f_9 * gsi1_147[k]
                   + f_3 * pc_y[k] * gsk_191[k];

        t_242[k] = f_6 * gsi0_148[k]
                   - f_7 * gsi1_148[k]
                   + f_3 * pc_y[k] * gsk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_y, fsk_200, gsi0_149, gsi0_160, \
                         gsi1_149, gsi1_160, gsk_193, gsk_194, \
                         gsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * gsi0_149[k]
                   - f_5 * gsi1_149[k]
                   + f_3 * pc_y[k] * gsk_193[k];

        t_244[k] = f_3 * pc_y[k] * gsk_194[k];

        t_245[k] = f_16 * fsk_200[k]
                   + f_6 * gsi0_160[k]
                   - f_7 * gsi1_160[k]
                   + f_3 * pc_x[k] * gsk_200[k];
    }
}

static auto
compute_prim_gsl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsl0,
                                                          const size_t fsk, const size_t fsl1,
                                                          const size_t gsi0, const size_t gsi1,
                                                          const size_t gsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsl0_135 = buffer.data(fsl0 + 135);
    const auto *fsl0_138 = buffer.data(fsl0 + 138);
    const auto *fsl0_141 = buffer.data(fsl0 + 141);
    const auto *fsl0_145 = buffer.data(fsl0 + 145);
    const auto *fsl0_150 = buffer.data(fsl0 + 150);
    const auto *fsl0_156 = buffer.data(fsl0 + 156);
    const auto *fsl0_225 = buffer.data(fsl0 + 225);
    const auto *fsl0_230 = buffer.data(fsl0 + 230);
    const auto *fsl0_270 = buffer.data(fsl0 + 270);
    const auto *fsl0_273 = buffer.data(fsl0 + 273);
    const auto *fsl0_276 = buffer.data(fsl0 + 276);
    const auto *fsl0_280 = buffer.data(fsl0 + 280);
    const auto *fsl0_285 = buffer.data(fsl0 + 285);
    const auto *fsl0_291 = buffer.data(fsl0 + 291);
    const auto *fsl0_306 = buffer.data(fsl0 + 306);
    const auto *fsl0_308 = buffer.data(fsl0 + 308);
    const auto *fsl0_309 = buffer.data(fsl0 + 309);
    const auto *fsl0_310 = buffer.data(fsl0 + 310);
    const auto *fsl0_311 = buffer.data(fsl0 + 311);
    const auto *fsl0_312 = buffer.data(fsl0 + 312);
    const auto *fsl0_314 = buffer.data(fsl0 + 314);
    const auto *fsl0_320 = buffer.data(fsl0 + 320);
    const auto *fsl0_324 = buffer.data(fsl0 + 324);
    const auto *fsl0_327 = buffer.data(fsl0 + 327);
    const auto *fsl0_329 = buffer.data(fsl0 + 329);
    const auto *fsl0_332 = buffer.data(fsl0 + 332);
    const auto *fsl0_333 = buffer.data(fsl0 + 333);
    const auto *fsl0_335 = buffer.data(fsl0 + 335);
    const auto *fsl0_338 = buffer.data(fsl0 + 338);
    const auto *fsl0_339 = buffer.data(fsl0 + 339);
    const auto *fsl0_340 = buffer.data(fsl0 + 340);
    const auto *fsl0_342 = buffer.data(fsl0 + 342);
    const auto *fsl0_351 = buffer.data(fsl0 + 351);
    const auto *fsl0_353 = buffer.data(fsl0 + 353);
    const auto *fsl0_354 = buffer.data(fsl0 + 354);
    const auto *fsl0_355 = buffer.data(fsl0 + 355);
    const auto *fsl0_356 = buffer.data(fsl0 + 356);
    const auto *fsl0_357 = buffer.data(fsl0 + 357);
    const auto *fsl0_359 = buffer.data(fsl0 + 359);
    const auto *fsl0_363 = buffer.data(fsl0 + 363);

    const auto *fsk_107 = buffer.data(fsk + 107);
    const auto *fsk_108 = buffer.data(fsk + 108);
    const auto *fsk_111 = buffer.data(fsk + 111);
    const auto *fsk_113 = buffer.data(fsk + 113);
    const auto *fsk_114 = buffer.data(fsk + 114);
    const auto *fsk_117 = buffer.data(fsk + 117);
    const auto *fsk_118 = buffer.data(fsk + 118);
    const auto *fsk_122 = buffer.data(fsk + 122);
    const auto *fsk_123 = buffer.data(fsk + 123);
    const auto *fsk_128 = buffer.data(fsk + 128);
    const auto *fsk_136 = buffer.data(fsk + 136);
    const auto *fsk_143 = buffer.data(fsk + 143);
    const auto *fsk_144 = buffer.data(fsk + 144);
    const auto *fsk_146 = buffer.data(fsk + 146);
    const auto *fsk_149 = buffer.data(fsk + 149);
    const auto *fsk_153 = buffer.data(fsk + 153);
    const auto *fsk_158 = buffer.data(fsk + 158);
    const auto *fsk_164 = buffer.data(fsk + 164);
    const auto *fsk_179 = buffer.data(fsk + 179);
    const auto *fsk_180 = buffer.data(fsk + 180);
    const auto *fsk_182 = buffer.data(fsk + 182);
    const auto *fsk_207 = buffer.data(fsk + 207);
    const auto *fsk_208 = buffer.data(fsk + 208);
    const auto *fsk_209 = buffer.data(fsk + 209);
    const auto *fsk_210 = buffer.data(fsk + 210);
    const auto *fsk_211 = buffer.data(fsk + 211);
    const auto *fsk_212 = buffer.data(fsk + 212);
    const auto *fsk_213 = buffer.data(fsk + 213);
    const auto *fsk_215 = buffer.data(fsk + 215);
    const auto *fsk_216 = buffer.data(fsk + 216);
    const auto *fsk_219 = buffer.data(fsk + 219);
    const auto *fsk_222 = buffer.data(fsk + 222);
    const auto *fsk_226 = buffer.data(fsk + 226);
    const auto *fsk_231 = buffer.data(fsk + 231);
    const auto *fsk_237 = buffer.data(fsk + 237);
    const auto *fsk_244 = buffer.data(fsk + 244);
    const auto *fsk_246 = buffer.data(fsk + 246);
    const auto *fsk_247 = buffer.data(fsk + 247);
    const auto *fsk_248 = buffer.data(fsk + 248);
    const auto *fsk_249 = buffer.data(fsk + 249);
    const auto *fsk_250 = buffer.data(fsk + 250);
    const auto *fsk_251 = buffer.data(fsk + 251);
    const auto *fsk_257 = buffer.data(fsk + 257);
    const auto *fsk_261 = buffer.data(fsk + 261);
    const auto *fsk_264 = buffer.data(fsk + 264);
    const auto *fsk_266 = buffer.data(fsk + 266);
    const auto *fsk_269 = buffer.data(fsk + 269);
    const auto *fsk_270 = buffer.data(fsk + 270);
    const auto *fsk_272 = buffer.data(fsk + 272);
    const auto *fsk_275 = buffer.data(fsk + 275);
    const auto *fsk_276 = buffer.data(fsk + 276);
    const auto *fsk_277 = buffer.data(fsk + 277);
    const auto *fsk_279 = buffer.data(fsk + 279);
    const auto *fsk_280 = buffer.data(fsk + 280);
    const auto *fsk_281 = buffer.data(fsk + 281);
    const auto *fsk_282 = buffer.data(fsk + 282);
    const auto *fsk_283 = buffer.data(fsk + 283);
    const auto *fsk_284 = buffer.data(fsk + 284);
    const auto *fsk_285 = buffer.data(fsk + 285);
    const auto *fsk_286 = buffer.data(fsk + 286);
    const auto *fsk_287 = buffer.data(fsk + 287);
    const auto *fsk_291 = buffer.data(fsk + 291);

    const auto *fsl1_135 = buffer.data(fsl1 + 135);
    const auto *fsl1_138 = buffer.data(fsl1 + 138);
    const auto *fsl1_141 = buffer.data(fsl1 + 141);
    const auto *fsl1_145 = buffer.data(fsl1 + 145);
    const auto *fsl1_150 = buffer.data(fsl1 + 150);
    const auto *fsl1_156 = buffer.data(fsl1 + 156);
    const auto *fsl1_225 = buffer.data(fsl1 + 225);
    const auto *fsl1_230 = buffer.data(fsl1 + 230);
    const auto *fsl1_270 = buffer.data(fsl1 + 270);
    const auto *fsl1_273 = buffer.data(fsl1 + 273);
    const auto *fsl1_276 = buffer.data(fsl1 + 276);
    const auto *fsl1_280 = buffer.data(fsl1 + 280);
    const auto *fsl1_285 = buffer.data(fsl1 + 285);
    const auto *fsl1_291 = buffer.data(fsl1 + 291);
    const auto *fsl1_306 = buffer.data(fsl1 + 306);
    const auto *fsl1_308 = buffer.data(fsl1 + 308);
    const auto *fsl1_309 = buffer.data(fsl1 + 309);
    const auto *fsl1_310 = buffer.data(fsl1 + 310);
    const auto *fsl1_311 = buffer.data(fsl1 + 311);
    const auto *fsl1_312 = buffer.data(fsl1 + 312);
    const auto *fsl1_314 = buffer.data(fsl1 + 314);
    const auto *fsl1_320 = buffer.data(fsl1 + 320);
    const auto *fsl1_324 = buffer.data(fsl1 + 324);
    const auto *fsl1_327 = buffer.data(fsl1 + 327);
    const auto *fsl1_329 = buffer.data(fsl1 + 329);
    const auto *fsl1_332 = buffer.data(fsl1 + 332);
    const auto *fsl1_333 = buffer.data(fsl1 + 333);
    const auto *fsl1_335 = buffer.data(fsl1 + 335);
    const auto *fsl1_338 = buffer.data(fsl1 + 338);
    const auto *fsl1_339 = buffer.data(fsl1 + 339);
    const auto *fsl1_340 = buffer.data(fsl1 + 340);
    const auto *fsl1_342 = buffer.data(fsl1 + 342);
    const auto *fsl1_351 = buffer.data(fsl1 + 351);
    const auto *fsl1_353 = buffer.data(fsl1 + 353);
    const auto *fsl1_354 = buffer.data(fsl1 + 354);
    const auto *fsl1_355 = buffer.data(fsl1 + 355);
    const auto *fsl1_356 = buffer.data(fsl1 + 356);
    const auto *fsl1_357 = buffer.data(fsl1 + 357);
    const auto *fsl1_359 = buffer.data(fsl1 + 359);
    const auto *fsl1_363 = buffer.data(fsl1 + 363);

    const auto *gsi0_150 = buffer.data(gsi0 + 150);
    const auto *gsi0_151 = buffer.data(gsi0 + 151);
    const auto *gsi0_152 = buffer.data(gsi0 + 152);
    const auto *gsi0_153 = buffer.data(gsi0 + 153);
    const auto *gsi0_154 = buffer.data(gsi0 + 154);
    const auto *gsi0_161 = buffer.data(gsi0 + 161);
    const auto *gsi0_162 = buffer.data(gsi0 + 162);
    const auto *gsi0_163 = buffer.data(gsi0 + 163);
    const auto *gsi0_164 = buffer.data(gsi0 + 164);
    const auto *gsi0_165 = buffer.data(gsi0 + 165);
    const auto *gsi0_166 = buffer.data(gsi0 + 166);
    const auto *gsi0_167 = buffer.data(gsi0 + 167);
    const auto *gsi0_168 = buffer.data(gsi0 + 168);
    const auto *gsi0_170 = buffer.data(gsi0 + 170);
    const auto *gsi0_171 = buffer.data(gsi0 + 171);
    const auto *gsi0_173 = buffer.data(gsi0 + 173);
    const auto *gsi0_174 = buffer.data(gsi0 + 174);
    const auto *gsi0_175 = buffer.data(gsi0 + 175);
    const auto *gsi0_177 = buffer.data(gsi0 + 177);
    const auto *gsi0_178 = buffer.data(gsi0 + 178);
    const auto *gsi0_179 = buffer.data(gsi0 + 179);
    const auto *gsi0_180 = buffer.data(gsi0 + 180);
    const auto *gsi0_182 = buffer.data(gsi0 + 182);

    const auto *gsi1_150 = buffer.data(gsi1 + 150);
    const auto *gsi1_151 = buffer.data(gsi1 + 151);
    const auto *gsi1_152 = buffer.data(gsi1 + 152);
    const auto *gsi1_153 = buffer.data(gsi1 + 153);
    const auto *gsi1_154 = buffer.data(gsi1 + 154);
    const auto *gsi1_161 = buffer.data(gsi1 + 161);
    const auto *gsi1_162 = buffer.data(gsi1 + 162);
    const auto *gsi1_163 = buffer.data(gsi1 + 163);
    const auto *gsi1_164 = buffer.data(gsi1 + 164);
    const auto *gsi1_165 = buffer.data(gsi1 + 165);
    const auto *gsi1_166 = buffer.data(gsi1 + 166);
    const auto *gsi1_167 = buffer.data(gsi1 + 167);
    const auto *gsi1_168 = buffer.data(gsi1 + 168);
    const auto *gsi1_170 = buffer.data(gsi1 + 170);
    const auto *gsi1_171 = buffer.data(gsi1 + 171);
    const auto *gsi1_173 = buffer.data(gsi1 + 173);
    const auto *gsi1_174 = buffer.data(gsi1 + 174);
    const auto *gsi1_175 = buffer.data(gsi1 + 175);
    const auto *gsi1_177 = buffer.data(gsi1 + 177);
    const auto *gsi1_178 = buffer.data(gsi1 + 178);
    const auto *gsi1_179 = buffer.data(gsi1 + 179);
    const auto *gsi1_180 = buffer.data(gsi1 + 180);
    const auto *gsi1_182 = buffer.data(gsi1 + 182);

    const auto *gsk_195 = buffer.data(gsk + 195);
    const auto *gsk_196 = buffer.data(gsk + 196);
    const auto *gsk_197 = buffer.data(gsk + 197);
    const auto *gsk_198 = buffer.data(gsk + 198);
    const auto *gsk_199 = buffer.data(gsk + 199);
    const auto *gsk_200 = buffer.data(gsk + 200);
    const auto *gsk_207 = buffer.data(gsk + 207);
    const auto *gsk_208 = buffer.data(gsk + 208);
    const auto *gsk_209 = buffer.data(gsk + 209);
    const auto *gsk_210 = buffer.data(gsk + 210);
    const auto *gsk_211 = buffer.data(gsk + 211);
    const auto *gsk_212 = buffer.data(gsk + 212);
    const auto *gsk_213 = buffer.data(gsk + 213);
    const auto *gsk_214 = buffer.data(gsk + 214);
    const auto *gsk_215 = buffer.data(gsk + 215);
    const auto *gsk_216 = buffer.data(gsk + 216);
    const auto *gsk_217 = buffer.data(gsk + 217);
    const auto *gsk_218 = buffer.data(gsk + 218);
    const auto *gsk_219 = buffer.data(gsk + 219);
    const auto *gsk_221 = buffer.data(gsk + 221);
    const auto *gsk_222 = buffer.data(gsk + 222);
    const auto *gsk_223 = buffer.data(gsk + 223);
    const auto *gsk_225 = buffer.data(gsk + 225);
    const auto *gsk_226 = buffer.data(gsk + 226);
    const auto *gsk_227 = buffer.data(gsk + 227);
    const auto *gsk_228 = buffer.data(gsk + 228);
    const auto *gsk_230 = buffer.data(gsk + 230);
    const auto *gsk_231 = buffer.data(gsk + 231);
    const auto *gsk_232 = buffer.data(gsk + 232);
    const auto *gsk_233 = buffer.data(gsk + 233);
    const auto *gsk_234 = buffer.data(gsk + 234);
    const auto *gsk_236 = buffer.data(gsk + 236);
    const auto *gsk_237 = buffer.data(gsk + 237);
    const auto *gsk_244 = buffer.data(gsk + 244);
    const auto *gsk_246 = buffer.data(gsk + 246);
    const auto *gsk_247 = buffer.data(gsk + 247);
    const auto *gsk_248 = buffer.data(gsk + 248);
    const auto *gsk_249 = buffer.data(gsk + 249);
    const auto *gsk_250 = buffer.data(gsk + 250);
    const auto *gsk_251 = buffer.data(gsk + 251);
    const auto *gsk_252 = buffer.data(gsk + 252);
    const auto *gsk_254 = buffer.data(gsk + 254);
    const auto *gsk_255 = buffer.data(gsk + 255);
    const auto *gsk_257 = buffer.data(gsk + 257);
    const auto *gsk_258 = buffer.data(gsk + 258);
    const auto *gsk_261 = buffer.data(gsk + 261);
    const auto *gsk_262 = buffer.data(gsk + 262);
    const auto *gsk_266 = buffer.data(gsk + 266);
    const auto *gsk_267 = buffer.data(gsk + 267);
    const auto *gsk_272 = buffer.data(gsk + 272);
    const auto *gsk_280 = buffer.data(gsk + 280);
    const auto *gsk_281 = buffer.data(gsk + 281);
    const auto *gsk_282 = buffer.data(gsk + 282);
    const auto *gsk_283 = buffer.data(gsk + 283);
    const auto *gsk_284 = buffer.data(gsk + 284);
    const auto *gsk_285 = buffer.data(gsk + 285);
    const auto *gsk_286 = buffer.data(gsk + 286);
    const auto *gsk_287 = buffer.data(gsk + 287);
    const auto *gsk_288 = buffer.data(gsk + 288);
    const auto *gsk_290 = buffer.data(gsk + 290);

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, gsi0_150, gsi0_151, gsi0_152, gsi1_150, \
                         gsi1_151, gsi1_152, gsk_195, gsk_196, \
                         gsk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * gsi0_150[k]
                   - f_13 * gsi1_150[k]
                   + f_3 * pc_y[k] * gsk_195[k];

        t_247[k] = f_10 * gsi0_151[k]
                   - f_11 * gsi1_151[k]
                   + f_3 * pc_y[k] * gsk_196[k];

        t_248[k] = f_8 * gsi0_152[k]
                   - f_9 * gsi1_152[k]
                   + f_3 * pc_y[k] * gsk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, gsi0_153, gsi0_154, gsi1_153, gsi1_154, \
                         gsk_198, gsk_199, gsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * gsi0_153[k]
                   - f_7 * gsi1_153[k]
                   + f_3 * pc_y[k] * gsk_198[k];

        t_250[k] = f_4 * gsi0_154[k]
                   - f_5 * gsi1_154[k]
                   + f_3 * pc_y[k] * gsk_199[k];

        t_251[k] = f_3 * pc_y[k] * gsk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, fsk_207, fsk_208, fsk_209, fsk_210, \
                         gsi0_167, gsi1_167, gsk_207, gsk_208, gsk_209, \
                         gsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_16 * fsk_207[k]
                   + f_4 * gsi0_167[k]
                   - f_5 * gsi1_167[k]
                   + f_3 * pc_x[k] * gsk_207[k];

        t_253[k] = f_16 * fsk_208[k]
                   + f_3 * pc_x[k] * gsk_208[k];

        t_254[k] = f_16 * fsk_209[k]
                   + f_3 * pc_x[k] * gsk_209[k];

        t_255[k] = f_16 * fsk_210[k]
                   + f_3 * pc_x[k] * gsk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, fsk_211, fsk_212, \
                         fsk_213, fsk_215, gsk_207, gsk_211, gsk_212, gsk_213, \
                         gsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_16 * fsk_211[k]
                   + f_3 * pc_x[k] * gsk_211[k];

        t_257[k] = f_16 * fsk_212[k]
                   + f_3 * pc_x[k] * gsk_212[k];

        t_258[k] = f_16 * fsk_213[k]
                   + f_3 * pc_x[k] * gsk_213[k];

        t_259[k] = f_3 * pc_y[k] * gsk_207[k];

        t_260[k] = f_16 * fsk_215[k]
                   + f_3 * pc_x[k] * gsk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_y, gsi0_161, gsi0_162, gsi0_163, gsi1_161, \
                         gsi1_162, gsi1_163, gsk_208, gsk_209, \
                         gsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * gsi0_161[k]
                   - f_2 * gsi1_161[k]
                   + f_3 * pc_y[k] * gsk_208[k];

        t_262[k] = f_20 * gsi0_162[k]
                   - f_21 * gsi1_162[k]
                   + f_3 * pc_y[k] * gsk_209[k];

        t_263[k] = f_12 * gsi0_163[k]
                   - f_13 * gsi1_163[k]
                   + f_3 * pc_y[k] * gsk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, gsi0_164, gsi0_165, gsi0_166, gsi1_164, \
                         gsi1_165, gsi1_166, gsk_211, gsk_212, \
                         gsk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * gsi0_164[k]
                   - f_11 * gsi1_164[k]
                   + f_3 * pc_y[k] * gsk_211[k];

        t_265[k] = f_8 * gsi0_165[k]
                   - f_9 * gsi1_165[k]
                   + f_3 * pc_y[k] * gsk_212[k];

        t_266[k] = f_6 * gsi0_166[k]
                   - f_7 * gsi1_166[k]
                   + f_3 * pc_y[k] * gsk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pa_x, pc_x, pc_y, pc_z, fsl0_270, \
                         fsk_107, fsk_216, fsl1_270, gsi0_167, gsi1_167, gsk_214, \
                         gsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * gsi0_167[k]
                   - f_5 * gsi1_167[k]
                   + f_3 * pc_y[k] * gsk_214[k];

        t_268[k] = f_3 * pc_y[k] * gsk_215[k];

        t_269[k] = f_16 * fsk_107[k]
                   + f_1 * gsi0_167[k]
                   - f_2 * gsi1_167[k]
                   + f_3 * pc_z[k] * gsk_215[k];

        t_270[k] = pa_x[k] * fsl0_270[k]
                   + f_22 * fsk_216[k]
                   - f_14 * pc_x[k] * fsl1_270[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_x, pc_x, pc_y, pc_z, fsl0_273, \
                         fsk_108, fsk_219, fsl1_273, gsk_216, gsk_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * fsk_108[k]
                   + f_3 * pc_y[k] * gsk_216[k];

        t_272[k] = f_3 * pc_z[k] * gsk_216[k];

        t_273[k] = pa_x[k] * fsl0_273[k]
                   + f_19 * fsk_219[k]
                   - f_14 * pc_x[k] * fsl1_273[k];

        t_274[k] = f_3 * pc_z[k] * gsk_217[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pa_x, pc_x, pc_z, fsl0_276, fsk_222, fsl1_276, \
                         gsi0_168, gsi1_168, gsk_218, gsk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_4 * gsi0_168[k]
                   - f_5 * gsi1_168[k]
                   + f_3 * pc_z[k] * gsk_218[k];

        t_276[k] = pa_x[k] * fsl0_276[k]
                   + f_18 * fsk_222[k]
                   - f_14 * pc_x[k] * fsl1_276[k];

        t_277[k] = f_3 * pc_z[k] * gsk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_x, pc_x, pc_y, pc_z, fsl0_280, \
                         fsk_113, fsk_226, fsl1_280, gsi0_170, gsi1_170, gsk_221, \
                         gsk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * fsk_113[k]
                   + f_3 * pc_y[k] * gsk_221[k];

        t_279[k] = f_6 * gsi0_170[k]
                   - f_7 * gsi1_170[k]
                   + f_3 * pc_z[k] * gsk_221[k];

        t_280[k] = pa_x[k] * fsl0_280[k]
                   + f_0 * fsk_226[k]
                   - f_14 * pc_x[k] * fsl1_280[k];

        t_281[k] = f_3 * pc_z[k] * gsk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pc_y, pc_z, fsk_117, gsi0_171, gsi0_173, \
                         gsi1_171, gsi1_173, gsk_223, gsk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * gsi0_171[k]
                   - f_5 * gsi1_171[k]
                   + f_3 * pc_z[k] * gsk_223[k];

        t_283[k] = f_17 * fsk_117[k]
                   + f_3 * pc_y[k] * gsk_225[k];

        t_284[k] = f_8 * gsi0_173[k]
                   - f_9 * gsi1_173[k]
                   + f_3 * pc_z[k] * gsk_225[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_x, pc_x, pc_z, fsl0_285, fsk_231, fsl1_285, \
                         gsi0_174, gsi1_174, gsk_226, gsk_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = pa_x[k] * fsl0_285[k]
                   + f_17 * fsk_231[k]
                   - f_14 * pc_x[k] * fsl1_285[k];

        t_286[k] = f_3 * pc_z[k] * gsk_226[k];

        t_287[k] = f_4 * gsi0_174[k]
                   - f_5 * gsi1_174[k]
                   + f_3 * pc_z[k] * gsk_227[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, fsk_122, gsi0_175, gsi0_177, \
                         gsi1_175, gsi1_177, gsk_228, gsk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * gsi0_175[k]
                   - f_7 * gsi1_175[k]
                   + f_3 * pc_z[k] * gsk_228[k];

        t_289[k] = f_17 * fsk_122[k]
                   + f_3 * pc_y[k] * gsk_230[k];

        t_290[k] = f_10 * gsi0_177[k]
                   - f_11 * gsi1_177[k]
                   + f_3 * pc_z[k] * gsk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_x, pc_x, pc_z, fsl0_291, fsk_237, fsl1_291, \
                         gsi0_178, gsi1_178, gsk_231, gsk_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = pa_x[k] * fsl0_291[k]
                   + f_16 * fsk_237[k]
                   - f_14 * pc_x[k] * fsl1_291[k];

        t_292[k] = f_3 * pc_z[k] * gsk_231[k];

        t_293[k] = f_4 * gsi0_178[k]
                   - f_5 * gsi1_178[k]
                   + f_3 * pc_z[k] * gsk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_y, pc_z, fsk_128, gsi0_179, gsi0_180, \
                         gsi0_182, gsi1_179, gsi1_180, gsi1_182, gsk_233, gsk_234, \
                         gsk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * gsi0_179[k]
                   - f_7 * gsi1_179[k]
                   + f_3 * pc_z[k] * gsk_233[k];

        t_295[k] = f_8 * gsi0_180[k]
                   - f_9 * gsi1_180[k]
                   + f_3 * pc_z[k] * gsk_234[k];

        t_296[k] = f_17 * fsk_128[k]
                   + f_3 * pc_y[k] * gsk_236[k];

        t_297[k] = f_12 * gsi0_182[k]
                   - f_13 * gsi1_182[k]
                   + f_3 * pc_z[k] * gsk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_z, fsk_244, fsk_246, \
                         fsk_247, fsk_248, gsk_237, gsk_244, gsk_246, gsk_247, \
                         gsk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_15 * fsk_244[k]
                   + f_3 * pc_x[k] * gsk_244[k];

        t_299[k] = f_3 * pc_z[k] * gsk_237[k];

        t_300[k] = f_15 * fsk_246[k]
                   + f_3 * pc_x[k] * gsk_246[k];

        t_301[k] = f_15 * fsk_247[k]
                   + f_3 * pc_x[k] * gsk_247[k];

        t_302[k] = f_15 * fsk_248[k]
                   + f_3 * pc_x[k] * gsk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pc_x, fsl0_306, fsk_249, fsk_250, \
                         fsk_251, fsl1_306, gsk_249, gsk_250, gsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_15 * fsk_249[k]
                   + f_3 * pc_x[k] * gsk_249[k];

        t_304[k] = f_15 * fsk_250[k]
                   + f_3 * pc_x[k] * gsk_250[k];

        t_305[k] = f_15 * fsk_251[k]
                   + f_3 * pc_x[k] * gsk_251[k];

        t_306[k] = pa_x[k] * fsl0_306[k]
                   - f_14 * pc_x[k] * fsl1_306[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_x, pc_x, pc_z, fsl0_308, fsl0_309, \
                         fsl0_310, fsl1_308, fsl1_309, fsl1_310, \
                         gsk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_z[k] * gsk_244[k];

        t_308[k] = pa_x[k] * fsl0_308[k]
                   - f_14 * pc_x[k] * fsl1_308[k];

        t_309[k] = pa_x[k] * fsl0_309[k]
                   - f_14 * pc_x[k] * fsl1_309[k];

        t_310[k] = pa_x[k] * fsl0_310[k]
                   - f_14 * pc_x[k] * fsl1_310[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pa_x, pc_x, pc_y, fsl0_311, fsl0_312, \
                         fsl0_314, fsk_143, fsl1_311, fsl1_312, fsl1_314, \
                         gsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = pa_x[k] * fsl0_311[k]
                   - f_14 * pc_x[k] * fsl1_311[k];

        t_312[k] = pa_x[k] * fsl0_312[k]
                   - f_14 * pc_x[k] * fsl1_312[k];

        t_313[k] = f_17 * fsk_143[k]
                   + f_3 * pc_y[k] * gsk_251[k];

        t_314[k] = pa_x[k] * fsl0_314[k]
                   - f_14 * pc_x[k] * fsl1_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_z, pc_y, pc_z, fsl0_135, fsl0_138, \
                         fsk_108, fsk_144, fsl1_135, fsl1_138, \
                         gsk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * fsl0_135[k]
                   - f_14 * pc_z[k] * fsl1_135[k];

        t_316[k] = f_16 * fsk_144[k]
                   + f_3 * pc_y[k] * gsk_252[k];

        t_317[k] = f_15 * fsk_108[k]
                   + f_3 * pc_z[k] * gsk_252[k];

        t_318[k] = pa_z[k] * fsl0_138[k]
                   - f_14 * pc_z[k] * fsl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_x, pa_z, pc_x, pc_y, pc_z, fsl0_141, \
                         fsl0_320, fsk_146, fsk_257, fsl1_141, fsl1_320, \
                         gsk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * fsk_146[k]
                   + f_3 * pc_y[k] * gsk_254[k];

        t_320[k] = pa_x[k] * fsl0_320[k]
                   + f_19 * fsk_257[k]
                   - f_14 * pc_x[k] * fsl1_320[k];

        t_321[k] = pa_z[k] * fsl0_141[k]
                   - f_14 * pc_z[k] * fsl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pa_x, pc_x, pc_y, pc_z, fsl0_324, fsk_111, \
                         fsk_149, fsk_261, fsl1_324, gsk_255, gsk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * fsk_111[k]
                   + f_3 * pc_z[k] * gsk_255[k];

        t_323[k] = f_16 * fsk_149[k]
                   + f_3 * pc_y[k] * gsk_257[k];

        t_324[k] = pa_x[k] * fsl0_324[k]
                   + f_18 * fsk_261[k]
                   - f_14 * pc_x[k] * fsl1_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, pa_x, pa_z, pc_x, pc_z, fsl0_145, fsl0_327, \
                         fsk_114, fsk_264, fsl1_145, fsl1_327, \
                         gsk_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * fsl0_145[k]
                   - f_14 * pc_z[k] * fsl1_145[k];

        t_326[k] = f_15 * fsk_114[k]
                   + f_3 * pc_z[k] * gsk_258[k];

        t_327[k] = pa_x[k] * fsl0_327[k]
                   + f_0 * fsk_264[k]
                   - f_14 * pc_x[k] * fsl1_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_x, pa_z, pc_x, pc_y, pc_z, fsl0_150, \
                         fsl0_329, fsk_153, fsk_266, fsl1_150, fsl1_329, \
                         gsk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_16 * fsk_153[k]
                   + f_3 * pc_y[k] * gsk_261[k];

        t_329[k] = pa_x[k] * fsl0_329[k]
                   + f_0 * fsk_266[k]
                   - f_14 * pc_x[k] * fsl1_329[k];

        t_330[k] = pa_z[k] * fsl0_150[k]
                   - f_14 * pc_z[k] * fsl1_150[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_x, pc_x, pc_z, fsl0_332, fsl0_333, fsk_118, \
                         fsk_269, fsk_270, fsl1_332, fsl1_333, \
                         gsk_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * fsk_118[k]
                   + f_3 * pc_z[k] * gsk_262[k];

        t_332[k] = pa_x[k] * fsl0_332[k]
                   + f_17 * fsk_269[k]
                   - f_14 * pc_x[k] * fsl1_332[k];

        t_333[k] = pa_x[k] * fsl0_333[k]
                   + f_17 * fsk_270[k]
                   - f_14 * pc_x[k] * fsl1_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pa_x, pa_z, pc_x, pc_y, pc_z, fsl0_156, \
                         fsl0_335, fsk_158, fsk_272, fsl1_156, fsl1_335, \
                         gsk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_16 * fsk_158[k]
                   + f_3 * pc_y[k] * gsk_266[k];

        t_335[k] = pa_x[k] * fsl0_335[k]
                   + f_17 * fsk_272[k]
                   - f_14 * pc_x[k] * fsl1_335[k];

        t_336[k] = pa_z[k] * fsl0_156[k]
                   - f_14 * pc_z[k] * fsl1_156[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pa_x, pc_x, pc_z, fsl0_338, fsl0_339, fsk_123, \
                         fsk_275, fsk_276, fsl1_338, fsl1_339, \
                         gsk_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_15 * fsk_123[k]
                   + f_3 * pc_z[k] * gsk_267[k];

        t_338[k] = pa_x[k] * fsl0_338[k]
                   + f_16 * fsk_275[k]
                   - f_14 * pc_x[k] * fsl1_338[k];

        t_339[k] = pa_x[k] * fsl0_339[k]
                   + f_16 * fsk_276[k]
                   - f_14 * pc_x[k] * fsl1_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_x, pc_x, pc_y, fsl0_340, fsl0_342, fsk_164, \
                         fsk_277, fsk_279, fsl1_340, fsl1_342, \
                         gsk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_x[k] * fsl0_340[k]
                   + f_16 * fsk_277[k]
                   - f_14 * pc_x[k] * fsl1_340[k];

        t_341[k] = f_16 * fsk_164[k]
                   + f_3 * pc_y[k] * gsk_272[k];

        t_342[k] = pa_x[k] * fsl0_342[k]
                   + f_16 * fsk_279[k]
                   - f_14 * pc_x[k] * fsl1_342[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, pc_x, fsk_280, fsk_281, fsk_282, \
                         fsk_283, fsk_284, gsk_280, gsk_281, gsk_282, gsk_283, \
                         gsk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_15 * fsk_280[k]
                   + f_3 * pc_x[k] * gsk_280[k];

        t_344[k] = f_15 * fsk_281[k]
                   + f_3 * pc_x[k] * gsk_281[k];

        t_345[k] = f_15 * fsk_282[k]
                   + f_3 * pc_x[k] * gsk_282[k];

        t_346[k] = f_15 * fsk_283[k]
                   + f_3 * pc_x[k] * gsk_283[k];

        t_347[k] = f_15 * fsk_284[k]
                   + f_3 * pc_x[k] * gsk_284[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_x, pc_x, fsl0_351, fsk_285, fsk_286, \
                         fsk_287, fsl1_351, gsk_285, gsk_286, gsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_15 * fsk_285[k]
                   + f_3 * pc_x[k] * gsk_285[k];

        t_349[k] = f_15 * fsk_286[k]
                   + f_3 * pc_x[k] * gsk_286[k];

        t_350[k] = f_15 * fsk_287[k]
                   + f_3 * pc_x[k] * gsk_287[k];

        t_351[k] = pa_x[k] * fsl0_351[k]
                   - f_14 * pc_x[k] * fsl1_351[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_x, pc_x, pc_z, fsl0_353, fsl0_354, \
                         fsl0_355, fsk_136, fsl1_353, fsl1_354, fsl1_355, \
                         gsk_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_15 * fsk_136[k]
                   + f_3 * pc_z[k] * gsk_280[k];

        t_353[k] = pa_x[k] * fsl0_353[k]
                   - f_14 * pc_x[k] * fsl1_353[k];

        t_354[k] = pa_x[k] * fsl0_354[k]
                   - f_14 * pc_x[k] * fsl1_354[k];

        t_355[k] = pa_x[k] * fsl0_355[k]
                   - f_14 * pc_x[k] * fsl1_355[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pc_x, pc_y, fsl0_356, fsl0_357, \
                         fsl0_359, fsk_179, fsl1_356, fsl1_357, fsl1_359, \
                         gsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * fsl0_356[k]
                   - f_14 * pc_x[k] * fsl1_356[k];

        t_357[k] = pa_x[k] * fsl0_357[k]
                   - f_14 * pc_x[k] * fsl1_357[k];

        t_358[k] = f_16 * fsk_179[k]
                   + f_3 * pc_y[k] * gsk_287[k];

        t_359[k] = pa_x[k] * fsl0_359[k]
                   - f_14 * pc_x[k] * fsl1_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pc_y, pc_z, fsl0_225, fsk_144, fsk_180, \
                         fsl1_225, gsk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pa_y[k] * fsl0_225[k]
                   - f_14 * pc_y[k] * fsl1_225[k];

        t_361[k] = f_15 * fsk_180[k]
                   + f_3 * pc_y[k] * gsk_288[k];

        t_362[k] = f_16 * fsk_144[k]
                   + f_3 * pc_z[k] * gsk_288[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_x, pa_y, pc_x, pc_y, fsl0_230, fsl0_363, \
                         fsk_182, fsk_291, fsl1_230, fsl1_363, \
                         gsk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_x[k] * fsl0_363[k]
                   + f_19 * fsk_291[k]
                   - f_14 * pc_x[k] * fsl1_363[k];

        t_364[k] = f_15 * fsk_182[k]
                   + f_3 * pc_y[k] * gsk_290[k];

        t_365[k] = pa_y[k] * fsl0_230[k]
                   - f_14 * pc_y[k] * fsl1_230[k];
    }
}

static auto
compute_prim_gsl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsl0,
                                                          const size_t fsk, const size_t fsl1,
                                                          const size_t gsi0, const size_t gsi1,
                                                          const size_t gsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
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
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsl0_234 = buffer.data(fsl0 + 234);
    const auto *fsl0_239 = buffer.data(fsl0 + 239);
    const auto *fsl0_245 = buffer.data(fsl0 + 245);
    const auto *fsl0_252 = buffer.data(fsl0 + 252);
    const auto *fsl0_366 = buffer.data(fsl0 + 366);
    const auto *fsl0_370 = buffer.data(fsl0 + 370);
    const auto *fsl0_372 = buffer.data(fsl0 + 372);
    const auto *fsl0_375 = buffer.data(fsl0 + 375);
    const auto *fsl0_377 = buffer.data(fsl0 + 377);
    const auto *fsl0_378 = buffer.data(fsl0 + 378);
    const auto *fsl0_381 = buffer.data(fsl0 + 381);
    const auto *fsl0_383 = buffer.data(fsl0 + 383);
    const auto *fsl0_384 = buffer.data(fsl0 + 384);
    const auto *fsl0_385 = buffer.data(fsl0 + 385);
    const auto *fsl0_396 = buffer.data(fsl0 + 396);
    const auto *fsl0_398 = buffer.data(fsl0 + 398);
    const auto *fsl0_399 = buffer.data(fsl0 + 399);
    const auto *fsl0_400 = buffer.data(fsl0 + 400);
    const auto *fsl0_401 = buffer.data(fsl0 + 401);
    const auto *fsl0_402 = buffer.data(fsl0 + 402);
    const auto *fsl0_404 = buffer.data(fsl0 + 404);
    const auto *fsl0_405 = buffer.data(fsl0 + 405);
    const auto *fsl0_410 = buffer.data(fsl0 + 410);
    const auto *fsl0_414 = buffer.data(fsl0 + 414);
    const auto *fsl0_419 = buffer.data(fsl0 + 419);
    const auto *fsl0_425 = buffer.data(fsl0 + 425);
    const auto *fsl0_432 = buffer.data(fsl0 + 432);
    const auto *fsl0_441 = buffer.data(fsl0 + 441);
    const auto *fsl0_442 = buffer.data(fsl0 + 442);
    const auto *fsl0_443 = buffer.data(fsl0 + 443);
    const auto *fsl0_444 = buffer.data(fsl0 + 444);
    const auto *fsl0_445 = buffer.data(fsl0 + 445);
    const auto *fsl0_446 = buffer.data(fsl0 + 446);
    const auto *fsl0_447 = buffer.data(fsl0 + 447);
    const auto *fsl0_449 = buffer.data(fsl0 + 449);

    const auto *fsk_147 = buffer.data(fsk + 147);
    const auto *fsk_150 = buffer.data(fsk + 150);
    const auto *fsk_154 = buffer.data(fsk + 154);
    const auto *fsk_159 = buffer.data(fsk + 159);
    const auto *fsk_172 = buffer.data(fsk + 172);
    const auto *fsk_180 = buffer.data(fsk + 180);
    const auto *fsk_185 = buffer.data(fsk + 185);
    const auto *fsk_189 = buffer.data(fsk + 189);
    const auto *fsk_194 = buffer.data(fsk + 194);
    const auto *fsk_200 = buffer.data(fsk + 200);
    const auto *fsk_215 = buffer.data(fsk + 215);
    const auto *fsk_244 = buffer.data(fsk + 244);
    const auto *fsk_294 = buffer.data(fsk + 294);
    const auto *fsk_298 = buffer.data(fsk + 298);
    const auto *fsk_300 = buffer.data(fsk + 300);
    const auto *fsk_303 = buffer.data(fsk + 303);
    const auto *fsk_305 = buffer.data(fsk + 305);
    const auto *fsk_306 = buffer.data(fsk + 306);
    const auto *fsk_309 = buffer.data(fsk + 309);
    const auto *fsk_311 = buffer.data(fsk + 311);
    const auto *fsk_312 = buffer.data(fsk + 312);
    const auto *fsk_313 = buffer.data(fsk + 313);
    const auto *fsk_316 = buffer.data(fsk + 316);
    const auto *fsk_317 = buffer.data(fsk + 317);
    const auto *fsk_318 = buffer.data(fsk + 318);
    const auto *fsk_319 = buffer.data(fsk + 319);
    const auto *fsk_320 = buffer.data(fsk + 320);
    const auto *fsk_321 = buffer.data(fsk + 321);
    const auto *fsk_322 = buffer.data(fsk + 322);
    const auto *fsk_323 = buffer.data(fsk + 323);
    const auto *fsk_324 = buffer.data(fsk + 324);
    const auto *fsk_329 = buffer.data(fsk + 329);
    const auto *fsk_333 = buffer.data(fsk + 333);
    const auto *fsk_338 = buffer.data(fsk + 338);
    const auto *fsk_344 = buffer.data(fsk + 344);
    const auto *fsk_351 = buffer.data(fsk + 351);
    const auto *fsk_352 = buffer.data(fsk + 352);
    const auto *fsk_353 = buffer.data(fsk + 353);
    const auto *fsk_354 = buffer.data(fsk + 354);
    const auto *fsk_355 = buffer.data(fsk + 355);
    const auto *fsk_356 = buffer.data(fsk + 356);
    const auto *fsk_357 = buffer.data(fsk + 357);
    const auto *fsk_359 = buffer.data(fsk + 359);

    const auto *fsl1_234 = buffer.data(fsl1 + 234);
    const auto *fsl1_239 = buffer.data(fsl1 + 239);
    const auto *fsl1_245 = buffer.data(fsl1 + 245);
    const auto *fsl1_252 = buffer.data(fsl1 + 252);
    const auto *fsl1_366 = buffer.data(fsl1 + 366);
    const auto *fsl1_370 = buffer.data(fsl1 + 370);
    const auto *fsl1_372 = buffer.data(fsl1 + 372);
    const auto *fsl1_375 = buffer.data(fsl1 + 375);
    const auto *fsl1_377 = buffer.data(fsl1 + 377);
    const auto *fsl1_378 = buffer.data(fsl1 + 378);
    const auto *fsl1_381 = buffer.data(fsl1 + 381);
    const auto *fsl1_383 = buffer.data(fsl1 + 383);
    const auto *fsl1_384 = buffer.data(fsl1 + 384);
    const auto *fsl1_385 = buffer.data(fsl1 + 385);
    const auto *fsl1_396 = buffer.data(fsl1 + 396);
    const auto *fsl1_398 = buffer.data(fsl1 + 398);
    const auto *fsl1_399 = buffer.data(fsl1 + 399);
    const auto *fsl1_400 = buffer.data(fsl1 + 400);
    const auto *fsl1_401 = buffer.data(fsl1 + 401);
    const auto *fsl1_402 = buffer.data(fsl1 + 402);
    const auto *fsl1_404 = buffer.data(fsl1 + 404);
    const auto *fsl1_405 = buffer.data(fsl1 + 405);
    const auto *fsl1_410 = buffer.data(fsl1 + 410);
    const auto *fsl1_414 = buffer.data(fsl1 + 414);
    const auto *fsl1_419 = buffer.data(fsl1 + 419);
    const auto *fsl1_425 = buffer.data(fsl1 + 425);
    const auto *fsl1_432 = buffer.data(fsl1 + 432);
    const auto *fsl1_441 = buffer.data(fsl1 + 441);
    const auto *fsl1_442 = buffer.data(fsl1 + 442);
    const auto *fsl1_443 = buffer.data(fsl1 + 443);
    const auto *fsl1_444 = buffer.data(fsl1 + 444);
    const auto *fsl1_445 = buffer.data(fsl1 + 445);
    const auto *fsl1_446 = buffer.data(fsl1 + 446);
    const auto *fsl1_447 = buffer.data(fsl1 + 447);
    const auto *fsl1_449 = buffer.data(fsl1 + 449);

    const auto *gsi0_252 = buffer.data(gsi0 + 252);
    const auto *gsi0_253 = buffer.data(gsi0 + 253);
    const auto *gsi0_254 = buffer.data(gsi0 + 254);
    const auto *gsi0_255 = buffer.data(gsi0 + 255);
    const auto *gsi0_256 = buffer.data(gsi0 + 256);
    const auto *gsi0_257 = buffer.data(gsi0 + 257);
    const auto *gsi0_258 = buffer.data(gsi0 + 258);
    const auto *gsi0_259 = buffer.data(gsi0 + 259);
    const auto *gsi0_260 = buffer.data(gsi0 + 260);
    const auto *gsi0_261 = buffer.data(gsi0 + 261);
    const auto *gsi0_262 = buffer.data(gsi0 + 262);
    const auto *gsi0_263 = buffer.data(gsi0 + 263);
    const auto *gsi0_264 = buffer.data(gsi0 + 264);
    const auto *gsi0_265 = buffer.data(gsi0 + 265);
    const auto *gsi0_266 = buffer.data(gsi0 + 266);
    const auto *gsi0_280 = buffer.data(gsi0 + 280);
    const auto *gsi0_281 = buffer.data(gsi0 + 281);
    const auto *gsi0_283 = buffer.data(gsi0 + 283);
    const auto *gsi0_285 = buffer.data(gsi0 + 285);
    const auto *gsi0_286 = buffer.data(gsi0 + 286);
    const auto *gsi0_288 = buffer.data(gsi0 + 288);
    const auto *gsi0_289 = buffer.data(gsi0 + 289);
    const auto *gsi0_290 = buffer.data(gsi0 + 290);
    const auto *gsi0_292 = buffer.data(gsi0 + 292);
    const auto *gsi0_293 = buffer.data(gsi0 + 293);
    const auto *gsi0_294 = buffer.data(gsi0 + 294);
    const auto *gsi0_295 = buffer.data(gsi0 + 295);
    const auto *gsi0_297 = buffer.data(gsi0 + 297);
    const auto *gsi0_298 = buffer.data(gsi0 + 298);
    const auto *gsi0_299 = buffer.data(gsi0 + 299);
    const auto *gsi0_300 = buffer.data(gsi0 + 300);
    const auto *gsi0_301 = buffer.data(gsi0 + 301);
    const auto *gsi0_302 = buffer.data(gsi0 + 302);
    const auto *gsi0_303 = buffer.data(gsi0 + 303);
    const auto *gsi0_304 = buffer.data(gsi0 + 304);
    const auto *gsi0_305 = buffer.data(gsi0 + 305);
    const auto *gsi0_306 = buffer.data(gsi0 + 306);
    const auto *gsi0_307 = buffer.data(gsi0 + 307);

    const auto *gsi1_252 = buffer.data(gsi1 + 252);
    const auto *gsi1_253 = buffer.data(gsi1 + 253);
    const auto *gsi1_254 = buffer.data(gsi1 + 254);
    const auto *gsi1_255 = buffer.data(gsi1 + 255);
    const auto *gsi1_256 = buffer.data(gsi1 + 256);
    const auto *gsi1_257 = buffer.data(gsi1 + 257);
    const auto *gsi1_258 = buffer.data(gsi1 + 258);
    const auto *gsi1_259 = buffer.data(gsi1 + 259);
    const auto *gsi1_260 = buffer.data(gsi1 + 260);
    const auto *gsi1_261 = buffer.data(gsi1 + 261);
    const auto *gsi1_262 = buffer.data(gsi1 + 262);
    const auto *gsi1_263 = buffer.data(gsi1 + 263);
    const auto *gsi1_264 = buffer.data(gsi1 + 264);
    const auto *gsi1_265 = buffer.data(gsi1 + 265);
    const auto *gsi1_266 = buffer.data(gsi1 + 266);
    const auto *gsi1_280 = buffer.data(gsi1 + 280);
    const auto *gsi1_281 = buffer.data(gsi1 + 281);
    const auto *gsi1_283 = buffer.data(gsi1 + 283);
    const auto *gsi1_285 = buffer.data(gsi1 + 285);
    const auto *gsi1_286 = buffer.data(gsi1 + 286);
    const auto *gsi1_288 = buffer.data(gsi1 + 288);
    const auto *gsi1_289 = buffer.data(gsi1 + 289);
    const auto *gsi1_290 = buffer.data(gsi1 + 290);
    const auto *gsi1_292 = buffer.data(gsi1 + 292);
    const auto *gsi1_293 = buffer.data(gsi1 + 293);
    const auto *gsi1_294 = buffer.data(gsi1 + 294);
    const auto *gsi1_295 = buffer.data(gsi1 + 295);
    const auto *gsi1_297 = buffer.data(gsi1 + 297);
    const auto *gsi1_298 = buffer.data(gsi1 + 298);
    const auto *gsi1_299 = buffer.data(gsi1 + 299);
    const auto *gsi1_300 = buffer.data(gsi1 + 300);
    const auto *gsi1_301 = buffer.data(gsi1 + 301);
    const auto *gsi1_302 = buffer.data(gsi1 + 302);
    const auto *gsi1_303 = buffer.data(gsi1 + 303);
    const auto *gsi1_304 = buffer.data(gsi1 + 304);
    const auto *gsi1_305 = buffer.data(gsi1 + 305);
    const auto *gsi1_306 = buffer.data(gsi1 + 306);
    const auto *gsi1_307 = buffer.data(gsi1 + 307);

    const auto *gsk_291 = buffer.data(gsk + 291);
    const auto *gsk_293 = buffer.data(gsk + 293);
    const auto *gsk_294 = buffer.data(gsk + 294);
    const auto *gsk_297 = buffer.data(gsk + 297);
    const auto *gsk_298 = buffer.data(gsk + 298);
    const auto *gsk_302 = buffer.data(gsk + 302);
    const auto *gsk_303 = buffer.data(gsk + 303);
    const auto *gsk_308 = buffer.data(gsk + 308);
    const auto *gsk_316 = buffer.data(gsk + 316);
    const auto *gsk_317 = buffer.data(gsk + 317);
    const auto *gsk_318 = buffer.data(gsk + 318);
    const auto *gsk_319 = buffer.data(gsk + 319);
    const auto *gsk_320 = buffer.data(gsk + 320);
    const auto *gsk_321 = buffer.data(gsk + 321);
    const auto *gsk_322 = buffer.data(gsk + 322);
    const auto *gsk_323 = buffer.data(gsk + 323);
    const auto *gsk_324 = buffer.data(gsk + 324);
    const auto *gsk_325 = buffer.data(gsk + 325);
    const auto *gsk_326 = buffer.data(gsk + 326);
    const auto *gsk_327 = buffer.data(gsk + 327);
    const auto *gsk_328 = buffer.data(gsk + 328);
    const auto *gsk_329 = buffer.data(gsk + 329);
    const auto *gsk_330 = buffer.data(gsk + 330);
    const auto *gsk_331 = buffer.data(gsk + 331);
    const auto *gsk_332 = buffer.data(gsk + 332);
    const auto *gsk_333 = buffer.data(gsk + 333);
    const auto *gsk_334 = buffer.data(gsk + 334);
    const auto *gsk_335 = buffer.data(gsk + 335);
    const auto *gsk_336 = buffer.data(gsk + 336);
    const auto *gsk_337 = buffer.data(gsk + 337);
    const auto *gsk_338 = buffer.data(gsk + 338);
    const auto *gsk_339 = buffer.data(gsk + 339);
    const auto *gsk_340 = buffer.data(gsk + 340);
    const auto *gsk_341 = buffer.data(gsk + 341);
    const auto *gsk_342 = buffer.data(gsk + 342);
    const auto *gsk_343 = buffer.data(gsk + 343);
    const auto *gsk_344 = buffer.data(gsk + 344);
    const auto *gsk_351 = buffer.data(gsk + 351);
    const auto *gsk_352 = buffer.data(gsk + 352);
    const auto *gsk_353 = buffer.data(gsk + 353);
    const auto *gsk_354 = buffer.data(gsk + 354);
    const auto *gsk_355 = buffer.data(gsk + 355);
    const auto *gsk_356 = buffer.data(gsk + 356);
    const auto *gsk_357 = buffer.data(gsk + 357);
    const auto *gsk_359 = buffer.data(gsk + 359);
    const auto *gsk_360 = buffer.data(gsk + 360);
    const auto *gsk_361 = buffer.data(gsk + 361);
    const auto *gsk_363 = buffer.data(gsk + 363);
    const auto *gsk_365 = buffer.data(gsk + 365);
    const auto *gsk_366 = buffer.data(gsk + 366);
    const auto *gsk_368 = buffer.data(gsk + 368);
    const auto *gsk_369 = buffer.data(gsk + 369);
    const auto *gsk_370 = buffer.data(gsk + 370);
    const auto *gsk_372 = buffer.data(gsk + 372);
    const auto *gsk_373 = buffer.data(gsk + 373);
    const auto *gsk_374 = buffer.data(gsk + 374);
    const auto *gsk_375 = buffer.data(gsk + 375);
    const auto *gsk_377 = buffer.data(gsk + 377);
    const auto *gsk_378 = buffer.data(gsk + 378);
    const auto *gsk_379 = buffer.data(gsk + 379);
    const auto *gsk_380 = buffer.data(gsk + 380);
    const auto *gsk_381 = buffer.data(gsk + 381);
    const auto *gsk_383 = buffer.data(gsk + 383);
    const auto *gsk_384 = buffer.data(gsk + 384);
    const auto *gsk_385 = buffer.data(gsk + 385);
    const auto *gsk_386 = buffer.data(gsk + 386);
    const auto *gsk_387 = buffer.data(gsk + 387);
    const auto *gsk_388 = buffer.data(gsk + 388);
    const auto *gsk_389 = buffer.data(gsk + 389);
    const auto *gsk_390 = buffer.data(gsk + 390);
    const auto *gsk_391 = buffer.data(gsk + 391);
    const auto *gsk_392 = buffer.data(gsk + 392);
    const auto *gsk_393 = buffer.data(gsk + 393);
    const auto *gsk_394 = buffer.data(gsk + 394);
    const auto *gsk_395 = buffer.data(gsk + 395);

#pragma omp simd aligned(t_366, t_367, t_368, pa_x, pc_x, pc_y, pc_z, fsl0_366, fsk_147, \
                         fsk_185, fsk_294, fsl1_366, gsk_291, gsk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_x[k] * fsl0_366[k]
                   + f_18 * fsk_294[k]
                   - f_14 * pc_x[k] * fsl1_366[k];

        t_367[k] = f_16 * fsk_147[k]
                   + f_3 * pc_z[k] * gsk_291[k];

        t_368[k] = f_15 * fsk_185[k]
                   + f_3 * pc_y[k] * gsk_293[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_x, pa_y, pc_x, pc_y, pc_z, fsl0_234, \
                         fsl0_370, fsk_150, fsk_298, fsl1_234, fsl1_370, \
                         gsk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pa_y[k] * fsl0_234[k]
                   - f_14 * pc_y[k] * fsl1_234[k];

        t_370[k] = pa_x[k] * fsl0_370[k]
                   + f_0 * fsk_298[k]
                   - f_14 * pc_x[k] * fsl1_370[k];

        t_371[k] = f_16 * fsk_150[k]
                   + f_3 * pc_z[k] * gsk_294[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_x, pa_y, pc_x, pc_y, fsl0_239, fsl0_372, \
                         fsk_189, fsk_300, fsl1_239, fsl1_372, \
                         gsk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = pa_x[k] * fsl0_372[k]
                   + f_0 * fsk_300[k]
                   - f_14 * pc_x[k] * fsl1_372[k];

        t_373[k] = f_15 * fsk_189[k]
                   + f_3 * pc_y[k] * gsk_297[k];

        t_374[k] = pa_y[k] * fsl0_239[k]
                   - f_14 * pc_y[k] * fsl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_x, pc_x, pc_z, fsl0_375, fsl0_377, fsk_154, \
                         fsk_303, fsk_305, fsl1_375, fsl1_377, \
                         gsk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_x[k] * fsl0_375[k]
                   + f_17 * fsk_303[k]
                   - f_14 * pc_x[k] * fsl1_375[k];

        t_376[k] = f_16 * fsk_154[k]
                   + f_3 * pc_z[k] * gsk_298[k];

        t_377[k] = pa_x[k] * fsl0_377[k]
                   + f_17 * fsk_305[k]
                   - f_14 * pc_x[k] * fsl1_377[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pa_x, pa_y, pc_x, pc_y, fsl0_245, fsl0_378, \
                         fsk_194, fsk_306, fsl1_245, fsl1_378, \
                         gsk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_x[k] * fsl0_378[k]
                   + f_17 * fsk_306[k]
                   - f_14 * pc_x[k] * fsl1_378[k];

        t_379[k] = f_15 * fsk_194[k]
                   + f_3 * pc_y[k] * gsk_302[k];

        t_380[k] = pa_y[k] * fsl0_245[k]
                   - f_14 * pc_y[k] * fsl1_245[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pa_x, pc_x, pc_z, fsl0_381, fsl0_383, fsk_159, \
                         fsk_309, fsk_311, fsl1_381, fsl1_383, \
                         gsk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = pa_x[k] * fsl0_381[k]
                   + f_16 * fsk_309[k]
                   - f_14 * pc_x[k] * fsl1_381[k];

        t_382[k] = f_16 * fsk_159[k]
                   + f_3 * pc_z[k] * gsk_303[k];

        t_383[k] = pa_x[k] * fsl0_383[k]
                   + f_16 * fsk_311[k]
                   - f_14 * pc_x[k] * fsl1_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pa_x, pc_x, pc_y, fsl0_384, fsl0_385, fsk_200, \
                         fsk_312, fsk_313, fsl1_384, fsl1_385, \
                         gsk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_x[k] * fsl0_384[k]
                   + f_16 * fsk_312[k]
                   - f_14 * pc_x[k] * fsl1_384[k];

        t_385[k] = pa_x[k] * fsl0_385[k]
                   + f_16 * fsk_313[k]
                   - f_14 * pc_x[k] * fsl1_385[k];

        t_386[k] = f_15 * fsk_200[k]
                   + f_3 * pc_y[k] * gsk_308[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_x, pc_y, fsl0_252, fsk_316, \
                         fsk_317, fsk_318, fsl1_252, gsk_316, gsk_317, \
                         gsk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pa_y[k] * fsl0_252[k]
                   - f_14 * pc_y[k] * fsl1_252[k];

        t_388[k] = f_15 * fsk_316[k]
                   + f_3 * pc_x[k] * gsk_316[k];

        t_389[k] = f_15 * fsk_317[k]
                   + f_3 * pc_x[k] * gsk_317[k];

        t_390[k] = f_15 * fsk_318[k]
                   + f_3 * pc_x[k] * gsk_318[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, pc_x, fsk_319, fsk_320, fsk_321, \
                         fsk_322, fsk_323, gsk_319, gsk_320, gsk_321, gsk_322, \
                         gsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_15 * fsk_319[k]
                   + f_3 * pc_x[k] * gsk_319[k];

        t_392[k] = f_15 * fsk_320[k]
                   + f_3 * pc_x[k] * gsk_320[k];

        t_393[k] = f_15 * fsk_321[k]
                   + f_3 * pc_x[k] * gsk_321[k];

        t_394[k] = f_15 * fsk_322[k]
                   + f_3 * pc_x[k] * gsk_322[k];

        t_395[k] = f_15 * fsk_323[k]
                   + f_3 * pc_x[k] * gsk_323[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pa_x, pc_x, pc_z, fsl0_396, fsl0_398, \
                         fsl0_399, fsk_172, fsl1_396, fsl1_398, fsl1_399, \
                         gsk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pa_x[k] * fsl0_396[k]
                   - f_14 * pc_x[k] * fsl1_396[k];

        t_397[k] = f_16 * fsk_172[k]
                   + f_3 * pc_z[k] * gsk_316[k];

        t_398[k] = pa_x[k] * fsl0_398[k]
                   - f_14 * pc_x[k] * fsl1_398[k];

        t_399[k] = pa_x[k] * fsl0_399[k]
                   - f_14 * pc_x[k] * fsl1_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, pa_x, pc_x, pc_y, fsl0_400, fsl0_401, \
                         fsl0_402, fsk_215, fsl1_400, fsl1_401, fsl1_402, \
                         gsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pa_x[k] * fsl0_400[k]
                   - f_14 * pc_x[k] * fsl1_400[k];

        t_401[k] = pa_x[k] * fsl0_401[k]
                   - f_14 * pc_x[k] * fsl1_401[k];

        t_402[k] = pa_x[k] * fsl0_402[k]
                   - f_14 * pc_x[k] * fsl1_402[k];

        t_403[k] = f_15 * fsk_215[k]
                   + f_3 * pc_y[k] * gsk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_x, pc_x, pc_y, pc_z, fsl0_404, \
                         fsl0_405, fsk_180, fsk_324, fsl1_404, fsl1_405, \
                         gsk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_x[k] * fsl0_404[k]
                   - f_14 * pc_x[k] * fsl1_404[k];

        t_405[k] = pa_x[k] * fsl0_405[k]
                   + f_22 * fsk_324[k]
                   - f_14 * pc_x[k] * fsl1_405[k];

        t_406[k] = f_3 * pc_y[k] * gsk_324[k];

        t_407[k] = f_17 * fsk_180[k]
                   + f_3 * pc_z[k] * gsk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pa_x, pc_x, pc_y, fsl0_410, fsk_329, fsl1_410, \
                         gsi0_252, gsi1_252, gsk_325, gsk_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_4 * gsi0_252[k]
                   - f_5 * gsi1_252[k]
                   + f_3 * pc_y[k] * gsk_325[k];

        t_409[k] = f_3 * pc_y[k] * gsk_326[k];

        t_410[k] = pa_x[k] * fsl0_410[k]
                   + f_19 * fsk_329[k]
                   - f_14 * pc_x[k] * fsl1_410[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, gsi0_253, gsi0_254, gsi1_253, gsi1_254, \
                         gsk_327, gsk_328, gsk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * gsi0_253[k]
                   - f_7 * gsi1_253[k]
                   + f_3 * pc_y[k] * gsk_327[k];

        t_412[k] = f_4 * gsi0_254[k]
                   - f_5 * gsi1_254[k]
                   + f_3 * pc_y[k] * gsk_328[k];

        t_413[k] = f_3 * pc_y[k] * gsk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pa_x, pc_x, pc_y, fsl0_414, fsk_333, fsl1_414, \
                         gsi0_255, gsi0_256, gsi1_255, gsi1_256, gsk_330, \
                         gsk_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pa_x[k] * fsl0_414[k]
                   + f_18 * fsk_333[k]
                   - f_14 * pc_x[k] * fsl1_414[k];

        t_415[k] = f_8 * gsi0_255[k]
                   - f_9 * gsi1_255[k]
                   + f_3 * pc_y[k] * gsk_330[k];

        t_416[k] = f_6 * gsi0_256[k]
                   - f_7 * gsi1_256[k]
                   + f_3 * pc_y[k] * gsk_331[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_x, pc_x, pc_y, fsl0_419, fsk_338, fsl1_419, \
                         gsi0_257, gsi1_257, gsk_332, gsk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * gsi0_257[k]
                   - f_5 * gsi1_257[k]
                   + f_3 * pc_y[k] * gsk_332[k];

        t_418[k] = f_3 * pc_y[k] * gsk_333[k];

        t_419[k] = pa_x[k] * fsl0_419[k]
                   + f_0 * fsk_338[k]
                   - f_14 * pc_x[k] * fsl1_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_y, gsi0_258, gsi0_259, gsi0_260, gsi1_258, \
                         gsi1_259, gsi1_260, gsk_334, gsk_335, \
                         gsk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * gsi0_258[k]
                   - f_11 * gsi1_258[k]
                   + f_3 * pc_y[k] * gsk_334[k];

        t_421[k] = f_8 * gsi0_259[k]
                   - f_9 * gsi1_259[k]
                   + f_3 * pc_y[k] * gsk_335[k];

        t_422[k] = f_6 * gsi0_260[k]
                   - f_7 * gsi1_260[k]
                   + f_3 * pc_y[k] * gsk_336[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pa_x, pc_x, pc_y, fsl0_425, fsk_344, fsl1_425, \
                         gsi0_261, gsi1_261, gsk_337, gsk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * gsi0_261[k]
                   - f_5 * gsi1_261[k]
                   + f_3 * pc_y[k] * gsk_337[k];

        t_424[k] = f_3 * pc_y[k] * gsk_338[k];

        t_425[k] = pa_x[k] * fsl0_425[k]
                   + f_17 * fsk_344[k]
                   - f_14 * pc_x[k] * fsl1_425[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_y, gsi0_262, gsi0_263, gsi0_264, gsi1_262, \
                         gsi1_263, gsi1_264, gsk_339, gsk_340, \
                         gsk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * gsi0_262[k]
                   - f_13 * gsi1_262[k]
                   + f_3 * pc_y[k] * gsk_339[k];

        t_427[k] = f_10 * gsi0_263[k]
                   - f_11 * gsi1_263[k]
                   + f_3 * pc_y[k] * gsk_340[k];

        t_428[k] = f_8 * gsi0_264[k]
                   - f_9 * gsi1_264[k]
                   + f_3 * pc_y[k] * gsk_341[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_y, gsi0_265, gsi0_266, gsi1_265, gsi1_266, \
                         gsk_342, gsk_343, gsk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * gsi0_265[k]
                   - f_7 * gsi1_265[k]
                   + f_3 * pc_y[k] * gsk_342[k];

        t_430[k] = f_4 * gsi0_266[k]
                   - f_5 * gsi1_266[k]
                   + f_3 * pc_y[k] * gsk_343[k];

        t_431[k] = f_3 * pc_y[k] * gsk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pa_x, pc_x, fsl0_432, fsk_351, fsk_352, \
                         fsk_353, fsk_354, fsl1_432, gsk_352, gsk_353, \
                         gsk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = pa_x[k] * fsl0_432[k]
                   + f_16 * fsk_351[k]
                   - f_14 * pc_x[k] * fsl1_432[k];

        t_433[k] = f_15 * fsk_352[k]
                   + f_3 * pc_x[k] * gsk_352[k];

        t_434[k] = f_15 * fsk_353[k]
                   + f_3 * pc_x[k] * gsk_353[k];

        t_435[k] = f_15 * fsk_354[k]
                   + f_3 * pc_x[k] * gsk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_y, fsk_355, fsk_356, \
                         fsk_357, fsk_359, gsk_351, gsk_355, gsk_356, gsk_357, \
                         gsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_15 * fsk_355[k]
                   + f_3 * pc_x[k] * gsk_355[k];

        t_437[k] = f_15 * fsk_356[k]
                   + f_3 * pc_x[k] * gsk_356[k];

        t_438[k] = f_15 * fsk_357[k]
                   + f_3 * pc_x[k] * gsk_357[k];

        t_439[k] = f_3 * pc_y[k] * gsk_351[k];

        t_440[k] = f_15 * fsk_359[k]
                   + f_3 * pc_x[k] * gsk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pa_x, pc_x, fsl0_441, fsl0_442, fsl0_443, \
                         fsl0_444, fsl1_441, fsl1_442, fsl1_443, \
                         fsl1_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pa_x[k] * fsl0_441[k]
                   - f_14 * pc_x[k] * fsl1_441[k];

        t_442[k] = pa_x[k] * fsl0_442[k]
                   - f_14 * pc_x[k] * fsl1_442[k];

        t_443[k] = pa_x[k] * fsl0_443[k]
                   - f_14 * pc_x[k] * fsl1_443[k];

        t_444[k] = pa_x[k] * fsl0_444[k]
                   - f_14 * pc_x[k] * fsl1_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_x, pc_x, pc_y, fsl0_445, fsl0_446, \
                         fsl0_447, fsl1_445, fsl1_446, fsl1_447, \
                         gsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_x[k] * fsl0_445[k]
                   - f_14 * pc_x[k] * fsl1_445[k];

        t_446[k] = pa_x[k] * fsl0_446[k]
                   - f_14 * pc_x[k] * fsl1_446[k];

        t_447[k] = pa_x[k] * fsl0_447[k]
                   - f_14 * pc_x[k] * fsl1_447[k];

        t_448[k] = f_3 * pc_y[k] * gsk_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_x, pc_x, pc_z, fsl0_449, fsl1_449, \
                         gsi0_280, gsi0_281, gsi1_280, gsi1_281, gsk_360, \
                         gsk_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_x[k] * fsl0_449[k]
                   - f_14 * pc_x[k] * fsl1_449[k];

        t_450[k] = f_1 * gsi0_280[k]
                   - f_2 * gsi1_280[k]
                   + f_3 * pc_x[k] * gsk_360[k];

        t_451[k] = f_20 * gsi0_281[k]
                   - f_21 * gsi1_281[k]
                   + f_3 * pc_x[k] * gsk_361[k];

        t_452[k] = f_3 * pc_z[k] * gsk_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_z, gsi0_283, gsi0_285, gsi0_286, \
                         gsi1_283, gsi1_285, gsi1_286, gsk_361, gsk_363, gsk_365, \
                         gsk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_12 * gsi0_283[k]
                   - f_13 * gsi1_283[k]
                   + f_3 * pc_x[k] * gsk_363[k];

        t_454[k] = f_3 * pc_z[k] * gsk_361[k];

        t_455[k] = f_12 * gsi0_285[k]
                   - f_13 * gsi1_285[k]
                   + f_3 * pc_x[k] * gsk_365[k];

        t_456[k] = f_10 * gsi0_286[k]
                   - f_11 * gsi1_286[k]
                   + f_3 * pc_x[k] * gsk_366[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_x, pc_z, gsi0_288, gsi0_289, gsi0_290, \
                         gsi1_288, gsi1_289, gsi1_290, gsk_363, gsk_368, gsk_369, \
                         gsk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * gsk_363[k];

        t_458[k] = f_10 * gsi0_288[k]
                   - f_11 * gsi1_288[k]
                   + f_3 * pc_x[k] * gsk_368[k];

        t_459[k] = f_10 * gsi0_289[k]
                   - f_11 * gsi1_289[k]
                   + f_3 * pc_x[k] * gsk_369[k];

        t_460[k] = f_8 * gsi0_290[k]
                   - f_9 * gsi1_290[k]
                   + f_3 * pc_x[k] * gsk_370[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pc_x, pc_z, gsi0_292, gsi0_293, gsi0_294, \
                         gsi1_292, gsi1_293, gsi1_294, gsk_366, gsk_372, gsk_373, \
                         gsk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_3 * pc_z[k] * gsk_366[k];

        t_462[k] = f_8 * gsi0_292[k]
                   - f_9 * gsi1_292[k]
                   + f_3 * pc_x[k] * gsk_372[k];

        t_463[k] = f_8 * gsi0_293[k]
                   - f_9 * gsi1_293[k]
                   + f_3 * pc_x[k] * gsk_373[k];

        t_464[k] = f_8 * gsi0_294[k]
                   - f_9 * gsi1_294[k]
                   + f_3 * pc_x[k] * gsk_374[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pc_x, pc_z, gsi0_295, gsi0_297, gsi0_298, \
                         gsi1_295, gsi1_297, gsi1_298, gsk_370, gsk_375, gsk_377, \
                         gsk_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_6 * gsi0_295[k]
                   - f_7 * gsi1_295[k]
                   + f_3 * pc_x[k] * gsk_375[k];

        t_466[k] = f_3 * pc_z[k] * gsk_370[k];

        t_467[k] = f_6 * gsi0_297[k]
                   - f_7 * gsi1_297[k]
                   + f_3 * pc_x[k] * gsk_377[k];

        t_468[k] = f_6 * gsi0_298[k]
                   - f_7 * gsi1_298[k]
                   + f_3 * pc_x[k] * gsk_378[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pc_x, pc_z, gsi0_299, gsi0_300, gsi0_301, \
                         gsi1_299, gsi1_300, gsi1_301, gsk_375, gsk_379, gsk_380, \
                         gsk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_6 * gsi0_299[k]
                   - f_7 * gsi1_299[k]
                   + f_3 * pc_x[k] * gsk_379[k];

        t_470[k] = f_6 * gsi0_300[k]
                   - f_7 * gsi1_300[k]
                   + f_3 * pc_x[k] * gsk_380[k];

        t_471[k] = f_4 * gsi0_301[k]
                   - f_5 * gsi1_301[k]
                   + f_3 * pc_x[k] * gsk_381[k];

        t_472[k] = f_3 * pc_z[k] * gsk_375[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_x, gsi0_303, gsi0_304, gsi0_305, gsi1_303, \
                         gsi1_304, gsi1_305, gsk_383, gsk_384, \
                         gsk_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_4 * gsi0_303[k]
                   - f_5 * gsi1_303[k]
                   + f_3 * pc_x[k] * gsk_383[k];

        t_474[k] = f_4 * gsi0_304[k]
                   - f_5 * gsi1_304[k]
                   + f_3 * pc_x[k] * gsk_384[k];

        t_475[k] = f_4 * gsi0_305[k]
                   - f_5 * gsi1_305[k]
                   + f_3 * pc_x[k] * gsk_385[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, pc_x, gsi0_306, gsi0_307, \
                         gsi1_306, gsi1_307, gsk_386, gsk_387, gsk_388, gsk_389, \
                         gsk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_4 * gsi0_306[k]
                   - f_5 * gsi1_306[k]
                   + f_3 * pc_x[k] * gsk_386[k];

        t_477[k] = f_4 * gsi0_307[k]
                   - f_5 * gsi1_307[k]
                   + f_3 * pc_x[k] * gsk_387[k];

        t_478[k] = f_3 * pc_x[k] * gsk_388[k];

        t_479[k] = f_3 * pc_x[k] * gsk_389[k];

        t_480[k] = f_3 * pc_x[k] * gsk_390[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, pc_x, gsk_391, gsk_392, gsk_393, \
                         gsk_394, gsk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_3 * pc_x[k] * gsk_391[k];

        t_482[k] = f_3 * pc_x[k] * gsk_392[k];

        t_483[k] = f_3 * pc_x[k] * gsk_393[k];

        t_484[k] = f_3 * pc_x[k] * gsk_394[k];

        t_485[k] = f_3 * pc_x[k] * gsk_395[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, pc_y, pc_z, fsk_244, gsi0_301, gsi0_302, \
                         gsi1_301, gsi1_302, gsk_388, gsk_389, \
                         gsk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_0 * fsk_244[k]
                   + f_1 * gsi0_301[k]
                   - f_2 * gsi1_301[k]
                   + f_3 * pc_y[k] * gsk_388[k];

        t_487[k] = f_3 * pc_z[k] * gsk_388[k];

        t_488[k] = f_4 * gsi0_301[k]
                   - f_5 * gsi1_301[k]
                   + f_3 * pc_z[k] * gsk_389[k];

        t_489[k] = f_6 * gsi0_302[k]
                   - f_7 * gsi1_302[k]
                   + f_3 * pc_z[k] * gsk_390[k];
    }
}

static auto
compute_prim_gsl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsl0,
                                                          const size_t fsk, const size_t fsl1,
                                                          const size_t gsi0, const size_t gsi1,
                                                          const size_t gsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsl0_270 = buffer.data(fsl0 + 270);
    const auto *fsl0_271 = buffer.data(fsl0 + 271);
    const auto *fsl0_273 = buffer.data(fsl0 + 273);
    const auto *fsl0_276 = buffer.data(fsl0 + 276);
    const auto *fsl0_280 = buffer.data(fsl0 + 280);
    const auto *fsl0_285 = buffer.data(fsl0 + 285);
    const auto *fsl0_291 = buffer.data(fsl0 + 291);
    const auto *fsl0_306 = buffer.data(fsl0 + 306);
    const auto *fsl0_308 = buffer.data(fsl0 + 308);
    const auto *fsl0_309 = buffer.data(fsl0 + 309);
    const auto *fsl0_310 = buffer.data(fsl0 + 310);
    const auto *fsl0_311 = buffer.data(fsl0 + 311);
    const auto *fsl0_312 = buffer.data(fsl0 + 312);
    const auto *fsl0_405 = buffer.data(fsl0 + 405);
    const auto *fsl0_407 = buffer.data(fsl0 + 407);
    const auto *fsl0_410 = buffer.data(fsl0 + 410);
    const auto *fsl0_414 = buffer.data(fsl0 + 414);
    const auto *fsl0_419 = buffer.data(fsl0 + 419);

    const auto *fsk_244 = buffer.data(fsk + 244);
    const auto *fsk_245 = buffer.data(fsk + 245);
    const auto *fsk_246 = buffer.data(fsk + 246);
    const auto *fsk_247 = buffer.data(fsk + 247);
    const auto *fsk_248 = buffer.data(fsk + 248);
    const auto *fsk_249 = buffer.data(fsk + 249);
    const auto *fsk_251 = buffer.data(fsk + 251);
    const auto *fsk_280 = buffer.data(fsk + 280);
    const auto *fsk_287 = buffer.data(fsk + 287);
    const auto *fsk_316 = buffer.data(fsk + 316);
    const auto *fsk_318 = buffer.data(fsk + 318);
    const auto *fsk_319 = buffer.data(fsk + 319);
    const auto *fsk_320 = buffer.data(fsk + 320);
    const auto *fsk_321 = buffer.data(fsk + 321);
    const auto *fsk_322 = buffer.data(fsk + 322);
    const auto *fsk_323 = buffer.data(fsk + 323);

    const auto *fsl1_270 = buffer.data(fsl1 + 270);
    const auto *fsl1_271 = buffer.data(fsl1 + 271);
    const auto *fsl1_273 = buffer.data(fsl1 + 273);
    const auto *fsl1_276 = buffer.data(fsl1 + 276);
    const auto *fsl1_280 = buffer.data(fsl1 + 280);
    const auto *fsl1_285 = buffer.data(fsl1 + 285);
    const auto *fsl1_291 = buffer.data(fsl1 + 291);
    const auto *fsl1_306 = buffer.data(fsl1 + 306);
    const auto *fsl1_308 = buffer.data(fsl1 + 308);
    const auto *fsl1_309 = buffer.data(fsl1 + 309);
    const auto *fsl1_310 = buffer.data(fsl1 + 310);
    const auto *fsl1_311 = buffer.data(fsl1 + 311);
    const auto *fsl1_312 = buffer.data(fsl1 + 312);
    const auto *fsl1_405 = buffer.data(fsl1 + 405);
    const auto *fsl1_407 = buffer.data(fsl1 + 407);
    const auto *fsl1_410 = buffer.data(fsl1 + 410);
    const auto *fsl1_414 = buffer.data(fsl1 + 414);
    const auto *fsl1_419 = buffer.data(fsl1 + 419);

    const auto *gsi0_303 = buffer.data(gsi0 + 303);
    const auto *gsi0_304 = buffer.data(gsi0 + 304);
    const auto *gsi0_305 = buffer.data(gsi0 + 305);
    const auto *gsi0_307 = buffer.data(gsi0 + 307);
    const auto *gsi0_310 = buffer.data(gsi0 + 310);
    const auto *gsi0_312 = buffer.data(gsi0 + 312);
    const auto *gsi0_313 = buffer.data(gsi0 + 313);
    const auto *gsi0_315 = buffer.data(gsi0 + 315);
    const auto *gsi0_316 = buffer.data(gsi0 + 316);
    const auto *gsi0_317 = buffer.data(gsi0 + 317);
    const auto *gsi0_319 = buffer.data(gsi0 + 319);
    const auto *gsi0_320 = buffer.data(gsi0 + 320);
    const auto *gsi0_321 = buffer.data(gsi0 + 321);
    const auto *gsi0_322 = buffer.data(gsi0 + 322);
    const auto *gsi0_324 = buffer.data(gsi0 + 324);
    const auto *gsi0_325 = buffer.data(gsi0 + 325);
    const auto *gsi0_326 = buffer.data(gsi0 + 326);
    const auto *gsi0_327 = buffer.data(gsi0 + 327);
    const auto *gsi0_328 = buffer.data(gsi0 + 328);
    const auto *gsi0_330 = buffer.data(gsi0 + 330);
    const auto *gsi0_331 = buffer.data(gsi0 + 331);
    const auto *gsi0_332 = buffer.data(gsi0 + 332);
    const auto *gsi0_333 = buffer.data(gsi0 + 333);
    const auto *gsi0_334 = buffer.data(gsi0 + 334);
    const auto *gsi0_335 = buffer.data(gsi0 + 335);
    const auto *gsi0_336 = buffer.data(gsi0 + 336);
    const auto *gsi0_337 = buffer.data(gsi0 + 337);
    const auto *gsi0_338 = buffer.data(gsi0 + 338);
    const auto *gsi0_339 = buffer.data(gsi0 + 339);
    const auto *gsi0_340 = buffer.data(gsi0 + 340);
    const auto *gsi0_341 = buffer.data(gsi0 + 341);
    const auto *gsi0_342 = buffer.data(gsi0 + 342);
    const auto *gsi0_343 = buffer.data(gsi0 + 343);
    const auto *gsi0_344 = buffer.data(gsi0 + 344);
    const auto *gsi0_345 = buffer.data(gsi0 + 345);
    const auto *gsi0_346 = buffer.data(gsi0 + 346);
    const auto *gsi0_347 = buffer.data(gsi0 + 347);
    const auto *gsi0_348 = buffer.data(gsi0 + 348);
    const auto *gsi0_349 = buffer.data(gsi0 + 349);
    const auto *gsi0_350 = buffer.data(gsi0 + 350);
    const auto *gsi0_351 = buffer.data(gsi0 + 351);
    const auto *gsi0_352 = buffer.data(gsi0 + 352);
    const auto *gsi0_353 = buffer.data(gsi0 + 353);
    const auto *gsi0_354 = buffer.data(gsi0 + 354);
    const auto *gsi0_355 = buffer.data(gsi0 + 355);
    const auto *gsi0_356 = buffer.data(gsi0 + 356);
    const auto *gsi0_357 = buffer.data(gsi0 + 357);
    const auto *gsi0_358 = buffer.data(gsi0 + 358);
    const auto *gsi0_359 = buffer.data(gsi0 + 359);
    const auto *gsi0_360 = buffer.data(gsi0 + 360);
    const auto *gsi0_361 = buffer.data(gsi0 + 361);
    const auto *gsi0_362 = buffer.data(gsi0 + 362);
    const auto *gsi0_363 = buffer.data(gsi0 + 363);
    const auto *gsi0_365 = buffer.data(gsi0 + 365);
    const auto *gsi0_367 = buffer.data(gsi0 + 367);
    const auto *gsi0_368 = buffer.data(gsi0 + 368);
    const auto *gsi0_370 = buffer.data(gsi0 + 370);
    const auto *gsi0_371 = buffer.data(gsi0 + 371);
    const auto *gsi0_372 = buffer.data(gsi0 + 372);
    const auto *gsi0_374 = buffer.data(gsi0 + 374);
    const auto *gsi0_375 = buffer.data(gsi0 + 375);
    const auto *gsi0_376 = buffer.data(gsi0 + 376);
    const auto *gsi0_377 = buffer.data(gsi0 + 377);
    const auto *gsi0_379 = buffer.data(gsi0 + 379);
    const auto *gsi0_380 = buffer.data(gsi0 + 380);
    const auto *gsi0_381 = buffer.data(gsi0 + 381);
    const auto *gsi0_382 = buffer.data(gsi0 + 382);

    const auto *gsi1_303 = buffer.data(gsi1 + 303);
    const auto *gsi1_304 = buffer.data(gsi1 + 304);
    const auto *gsi1_305 = buffer.data(gsi1 + 305);
    const auto *gsi1_307 = buffer.data(gsi1 + 307);
    const auto *gsi1_310 = buffer.data(gsi1 + 310);
    const auto *gsi1_312 = buffer.data(gsi1 + 312);
    const auto *gsi1_313 = buffer.data(gsi1 + 313);
    const auto *gsi1_315 = buffer.data(gsi1 + 315);
    const auto *gsi1_316 = buffer.data(gsi1 + 316);
    const auto *gsi1_317 = buffer.data(gsi1 + 317);
    const auto *gsi1_319 = buffer.data(gsi1 + 319);
    const auto *gsi1_320 = buffer.data(gsi1 + 320);
    const auto *gsi1_321 = buffer.data(gsi1 + 321);
    const auto *gsi1_322 = buffer.data(gsi1 + 322);
    const auto *gsi1_324 = buffer.data(gsi1 + 324);
    const auto *gsi1_325 = buffer.data(gsi1 + 325);
    const auto *gsi1_326 = buffer.data(gsi1 + 326);
    const auto *gsi1_327 = buffer.data(gsi1 + 327);
    const auto *gsi1_328 = buffer.data(gsi1 + 328);
    const auto *gsi1_330 = buffer.data(gsi1 + 330);
    const auto *gsi1_331 = buffer.data(gsi1 + 331);
    const auto *gsi1_332 = buffer.data(gsi1 + 332);
    const auto *gsi1_333 = buffer.data(gsi1 + 333);
    const auto *gsi1_334 = buffer.data(gsi1 + 334);
    const auto *gsi1_335 = buffer.data(gsi1 + 335);
    const auto *gsi1_336 = buffer.data(gsi1 + 336);
    const auto *gsi1_337 = buffer.data(gsi1 + 337);
    const auto *gsi1_338 = buffer.data(gsi1 + 338);
    const auto *gsi1_339 = buffer.data(gsi1 + 339);
    const auto *gsi1_340 = buffer.data(gsi1 + 340);
    const auto *gsi1_341 = buffer.data(gsi1 + 341);
    const auto *gsi1_342 = buffer.data(gsi1 + 342);
    const auto *gsi1_343 = buffer.data(gsi1 + 343);
    const auto *gsi1_344 = buffer.data(gsi1 + 344);
    const auto *gsi1_345 = buffer.data(gsi1 + 345);
    const auto *gsi1_346 = buffer.data(gsi1 + 346);
    const auto *gsi1_347 = buffer.data(gsi1 + 347);
    const auto *gsi1_348 = buffer.data(gsi1 + 348);
    const auto *gsi1_349 = buffer.data(gsi1 + 349);
    const auto *gsi1_350 = buffer.data(gsi1 + 350);
    const auto *gsi1_351 = buffer.data(gsi1 + 351);
    const auto *gsi1_352 = buffer.data(gsi1 + 352);
    const auto *gsi1_353 = buffer.data(gsi1 + 353);
    const auto *gsi1_354 = buffer.data(gsi1 + 354);
    const auto *gsi1_355 = buffer.data(gsi1 + 355);
    const auto *gsi1_356 = buffer.data(gsi1 + 356);
    const auto *gsi1_357 = buffer.data(gsi1 + 357);
    const auto *gsi1_358 = buffer.data(gsi1 + 358);
    const auto *gsi1_359 = buffer.data(gsi1 + 359);
    const auto *gsi1_360 = buffer.data(gsi1 + 360);
    const auto *gsi1_361 = buffer.data(gsi1 + 361);
    const auto *gsi1_362 = buffer.data(gsi1 + 362);
    const auto *gsi1_363 = buffer.data(gsi1 + 363);
    const auto *gsi1_365 = buffer.data(gsi1 + 365);
    const auto *gsi1_367 = buffer.data(gsi1 + 367);
    const auto *gsi1_368 = buffer.data(gsi1 + 368);
    const auto *gsi1_370 = buffer.data(gsi1 + 370);
    const auto *gsi1_371 = buffer.data(gsi1 + 371);
    const auto *gsi1_372 = buffer.data(gsi1 + 372);
    const auto *gsi1_374 = buffer.data(gsi1 + 374);
    const auto *gsi1_375 = buffer.data(gsi1 + 375);
    const auto *gsi1_376 = buffer.data(gsi1 + 376);
    const auto *gsi1_377 = buffer.data(gsi1 + 377);
    const auto *gsi1_379 = buffer.data(gsi1 + 379);
    const auto *gsi1_380 = buffer.data(gsi1 + 380);
    const auto *gsi1_381 = buffer.data(gsi1 + 381);
    const auto *gsi1_382 = buffer.data(gsi1 + 382);

    const auto *gsk_391 = buffer.data(gsk + 391);
    const auto *gsk_392 = buffer.data(gsk + 392);
    const auto *gsk_393 = buffer.data(gsk + 393);
    const auto *gsk_395 = buffer.data(gsk + 395);
    const auto *gsk_398 = buffer.data(gsk + 398);
    const auto *gsk_400 = buffer.data(gsk + 400);
    const auto *gsk_401 = buffer.data(gsk + 401);
    const auto *gsk_403 = buffer.data(gsk + 403);
    const auto *gsk_404 = buffer.data(gsk + 404);
    const auto *gsk_405 = buffer.data(gsk + 405);
    const auto *gsk_407 = buffer.data(gsk + 407);
    const auto *gsk_408 = buffer.data(gsk + 408);
    const auto *gsk_409 = buffer.data(gsk + 409);
    const auto *gsk_410 = buffer.data(gsk + 410);
    const auto *gsk_412 = buffer.data(gsk + 412);
    const auto *gsk_413 = buffer.data(gsk + 413);
    const auto *gsk_414 = buffer.data(gsk + 414);
    const auto *gsk_415 = buffer.data(gsk + 415);
    const auto *gsk_416 = buffer.data(gsk + 416);
    const auto *gsk_418 = buffer.data(gsk + 418);
    const auto *gsk_419 = buffer.data(gsk + 419);
    const auto *gsk_420 = buffer.data(gsk + 420);
    const auto *gsk_421 = buffer.data(gsk + 421);
    const auto *gsk_422 = buffer.data(gsk + 422);
    const auto *gsk_423 = buffer.data(gsk + 423);
    const auto *gsk_424 = buffer.data(gsk + 424);
    const auto *gsk_425 = buffer.data(gsk + 425);
    const auto *gsk_426 = buffer.data(gsk + 426);
    const auto *gsk_427 = buffer.data(gsk + 427);
    const auto *gsk_428 = buffer.data(gsk + 428);
    const auto *gsk_429 = buffer.data(gsk + 429);
    const auto *gsk_430 = buffer.data(gsk + 430);
    const auto *gsk_431 = buffer.data(gsk + 431);
    const auto *gsk_432 = buffer.data(gsk + 432);
    const auto *gsk_433 = buffer.data(gsk + 433);
    const auto *gsk_434 = buffer.data(gsk + 434);
    const auto *gsk_435 = buffer.data(gsk + 435);
    const auto *gsk_436 = buffer.data(gsk + 436);
    const auto *gsk_437 = buffer.data(gsk + 437);
    const auto *gsk_438 = buffer.data(gsk + 438);
    const auto *gsk_439 = buffer.data(gsk + 439);
    const auto *gsk_440 = buffer.data(gsk + 440);
    const auto *gsk_441 = buffer.data(gsk + 441);
    const auto *gsk_442 = buffer.data(gsk + 442);
    const auto *gsk_443 = buffer.data(gsk + 443);
    const auto *gsk_444 = buffer.data(gsk + 444);
    const auto *gsk_445 = buffer.data(gsk + 445);
    const auto *gsk_446 = buffer.data(gsk + 446);
    const auto *gsk_447 = buffer.data(gsk + 447);
    const auto *gsk_448 = buffer.data(gsk + 448);
    const auto *gsk_449 = buffer.data(gsk + 449);
    const auto *gsk_450 = buffer.data(gsk + 450);
    const auto *gsk_451 = buffer.data(gsk + 451);
    const auto *gsk_452 = buffer.data(gsk + 452);
    const auto *gsk_453 = buffer.data(gsk + 453);
    const auto *gsk_454 = buffer.data(gsk + 454);
    const auto *gsk_455 = buffer.data(gsk + 455);
    const auto *gsk_456 = buffer.data(gsk + 456);
    const auto *gsk_457 = buffer.data(gsk + 457);
    const auto *gsk_458 = buffer.data(gsk + 458);
    const auto *gsk_459 = buffer.data(gsk + 459);
    const auto *gsk_460 = buffer.data(gsk + 460);
    const auto *gsk_461 = buffer.data(gsk + 461);
    const auto *gsk_462 = buffer.data(gsk + 462);
    const auto *gsk_463 = buffer.data(gsk + 463);
    const auto *gsk_464 = buffer.data(gsk + 464);
    const auto *gsk_465 = buffer.data(gsk + 465);
    const auto *gsk_466 = buffer.data(gsk + 466);
    const auto *gsk_467 = buffer.data(gsk + 467);
    const auto *gsk_469 = buffer.data(gsk + 469);
    const auto *gsk_471 = buffer.data(gsk + 471);
    const auto *gsk_472 = buffer.data(gsk + 472);
    const auto *gsk_474 = buffer.data(gsk + 474);
    const auto *gsk_475 = buffer.data(gsk + 475);
    const auto *gsk_476 = buffer.data(gsk + 476);
    const auto *gsk_478 = buffer.data(gsk + 478);
    const auto *gsk_479 = buffer.data(gsk + 479);
    const auto *gsk_480 = buffer.data(gsk + 480);
    const auto *gsk_481 = buffer.data(gsk + 481);
    const auto *gsk_483 = buffer.data(gsk + 483);
    const auto *gsk_484 = buffer.data(gsk + 484);
    const auto *gsk_485 = buffer.data(gsk + 485);
    const auto *gsk_486 = buffer.data(gsk + 486);

#pragma omp simd aligned(t_490, t_491, t_492, pc_z, gsi0_303, gsi0_304, gsi0_305, gsi1_303, \
                         gsi1_304, gsi1_305, gsk_391, gsk_392, \
                         gsk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_8 * gsi0_303[k]
                   - f_9 * gsi1_303[k]
                   + f_3 * pc_z[k] * gsk_391[k];

        t_491[k] = f_10 * gsi0_304[k]
                   - f_11 * gsi1_304[k]
                   + f_3 * pc_z[k] * gsk_392[k];

        t_492[k] = f_12 * gsi0_305[k]
                   - f_13 * gsi1_305[k]
                   + f_3 * pc_z[k] * gsk_393[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, pa_z, pc_y, pc_z, fsl0_270, fsl0_271, \
                         fsk_251, fsl1_270, fsl1_271, gsi0_307, gsi1_307, \
                         gsk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_0 * fsk_251[k]
                   + f_3 * pc_y[k] * gsk_395[k];

        t_494[k] = f_1 * gsi0_307[k]
                   - f_2 * gsi1_307[k]
                   + f_3 * pc_z[k] * gsk_395[k];

        t_495[k] = pa_z[k] * fsl0_270[k]
                   - f_14 * pc_z[k] * fsl1_270[k];

        t_496[k] = pa_z[k] * fsl0_271[k]
                   - f_14 * pc_z[k] * fsl1_271[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pa_z, pc_x, pc_z, fsl0_273, fsl1_273, gsi0_310, \
                         gsi0_312, gsi1_310, gsi1_312, gsk_398, \
                         gsk_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_20 * gsi0_310[k]
                   - f_21 * gsi1_310[k]
                   + f_3 * pc_x[k] * gsk_398[k];

        t_498[k] = pa_z[k] * fsl0_273[k]
                   - f_14 * pc_z[k] * fsl1_273[k];

        t_499[k] = f_12 * gsi0_312[k]
                   - f_13 * gsi1_312[k]
                   + f_3 * pc_x[k] * gsk_400[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pa_z, pc_x, pc_z, fsl0_276, fsl1_276, gsi0_313, \
                         gsi0_315, gsi1_313, gsi1_315, gsk_401, \
                         gsk_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_12 * gsi0_313[k]
                   - f_13 * gsi1_313[k]
                   + f_3 * pc_x[k] * gsk_401[k];

        t_501[k] = pa_z[k] * fsl0_276[k]
                   - f_14 * pc_z[k] * fsl1_276[k];

        t_502[k] = f_10 * gsi0_315[k]
                   - f_11 * gsi1_315[k]
                   + f_3 * pc_x[k] * gsk_403[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pa_z, pc_x, pc_z, fsl0_280, fsl1_280, gsi0_316, \
                         gsi0_317, gsi1_316, gsi1_317, gsk_404, \
                         gsk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_10 * gsi0_316[k]
                   - f_11 * gsi1_316[k]
                   + f_3 * pc_x[k] * gsk_404[k];

        t_504[k] = f_10 * gsi0_317[k]
                   - f_11 * gsi1_317[k]
                   + f_3 * pc_x[k] * gsk_405[k];

        t_505[k] = pa_z[k] * fsl0_280[k]
                   - f_14 * pc_z[k] * fsl1_280[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pc_x, gsi0_319, gsi0_320, gsi0_321, gsi1_319, \
                         gsi1_320, gsi1_321, gsk_407, gsk_408, \
                         gsk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_8 * gsi0_319[k]
                   - f_9 * gsi1_319[k]
                   + f_3 * pc_x[k] * gsk_407[k];

        t_507[k] = f_8 * gsi0_320[k]
                   - f_9 * gsi1_320[k]
                   + f_3 * pc_x[k] * gsk_408[k];

        t_508[k] = f_8 * gsi0_321[k]
                   - f_9 * gsi1_321[k]
                   + f_3 * pc_x[k] * gsk_409[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_z, pc_x, pc_z, fsl0_285, fsl1_285, gsi0_322, \
                         gsi0_324, gsi1_322, gsi1_324, gsk_410, \
                         gsk_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_8 * gsi0_322[k]
                   - f_9 * gsi1_322[k]
                   + f_3 * pc_x[k] * gsk_410[k];

        t_510[k] = pa_z[k] * fsl0_285[k]
                   - f_14 * pc_z[k] * fsl1_285[k];

        t_511[k] = f_6 * gsi0_324[k]
                   - f_7 * gsi1_324[k]
                   + f_3 * pc_x[k] * gsk_412[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, gsi0_325, gsi0_326, gsi0_327, gsi1_325, \
                         gsi1_326, gsi1_327, gsk_413, gsk_414, \
                         gsk_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_6 * gsi0_325[k]
                   - f_7 * gsi1_325[k]
                   + f_3 * pc_x[k] * gsk_413[k];

        t_513[k] = f_6 * gsi0_326[k]
                   - f_7 * gsi1_326[k]
                   + f_3 * pc_x[k] * gsk_414[k];

        t_514[k] = f_6 * gsi0_327[k]
                   - f_7 * gsi1_327[k]
                   + f_3 * pc_x[k] * gsk_415[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pa_z, pc_x, pc_z, fsl0_291, fsl1_291, gsi0_328, \
                         gsi0_330, gsi1_328, gsi1_330, gsk_416, \
                         gsk_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_6 * gsi0_328[k]
                   - f_7 * gsi1_328[k]
                   + f_3 * pc_x[k] * gsk_416[k];

        t_516[k] = pa_z[k] * fsl0_291[k]
                   - f_14 * pc_z[k] * fsl1_291[k];

        t_517[k] = f_4 * gsi0_330[k]
                   - f_5 * gsi1_330[k]
                   + f_3 * pc_x[k] * gsk_418[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pc_x, gsi0_331, gsi0_332, gsi0_333, gsi1_331, \
                         gsi1_332, gsi1_333, gsk_419, gsk_420, \
                         gsk_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_4 * gsi0_331[k]
                   - f_5 * gsi1_331[k]
                   + f_3 * pc_x[k] * gsk_419[k];

        t_519[k] = f_4 * gsi0_332[k]
                   - f_5 * gsi1_332[k]
                   + f_3 * pc_x[k] * gsk_420[k];

        t_520[k] = f_4 * gsi0_333[k]
                   - f_5 * gsi1_333[k]
                   + f_3 * pc_x[k] * gsk_421[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, pc_x, gsi0_334, gsi0_335, \
                         gsi1_334, gsi1_335, gsk_422, gsk_423, gsk_424, gsk_425, \
                         gsk_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_4 * gsi0_334[k]
                   - f_5 * gsi1_334[k]
                   + f_3 * pc_x[k] * gsk_422[k];

        t_522[k] = f_4 * gsi0_335[k]
                   - f_5 * gsi1_335[k]
                   + f_3 * pc_x[k] * gsk_423[k];

        t_523[k] = f_3 * pc_x[k] * gsk_424[k];

        t_524[k] = f_3 * pc_x[k] * gsk_425[k];

        t_525[k] = f_3 * pc_x[k] * gsk_426[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, t_531, pa_z, pc_x, pc_z, fsl0_306, \
                         fsl1_306, gsk_427, gsk_428, gsk_429, gsk_430, \
                         gsk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_3 * pc_x[k] * gsk_427[k];

        t_527[k] = f_3 * pc_x[k] * gsk_428[k];

        t_528[k] = f_3 * pc_x[k] * gsk_429[k];

        t_529[k] = f_3 * pc_x[k] * gsk_430[k];

        t_530[k] = f_3 * pc_x[k] * gsk_431[k];

        t_531[k] = pa_z[k] * fsl0_306[k]
                   - f_14 * pc_z[k] * fsl1_306[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, pa_z, pc_z, fsl0_308, fsl0_309, fsk_244, \
                         fsk_245, fsk_246, fsl1_308, fsl1_309, \
                         gsk_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_15 * fsk_244[k]
                   + f_3 * pc_z[k] * gsk_424[k];

        t_533[k] = pa_z[k] * fsl0_308[k]
                   + f_16 * fsk_245[k]
                   - f_14 * pc_z[k] * fsl1_308[k];

        t_534[k] = pa_z[k] * fsl0_309[k]
                   + f_17 * fsk_246[k]
                   - f_14 * pc_z[k] * fsl1_309[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, pa_z, pc_z, fsl0_310, fsl0_311, fsl0_312, \
                         fsk_247, fsk_248, fsk_249, fsl1_310, fsl1_311, \
                         fsl1_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = pa_z[k] * fsl0_310[k]
                   + f_0 * fsk_247[k]
                   - f_14 * pc_z[k] * fsl1_310[k];

        t_536[k] = pa_z[k] * fsl0_311[k]
                   + f_18 * fsk_248[k]
                   - f_14 * pc_z[k] * fsl1_311[k];

        t_537[k] = pa_z[k] * fsl0_312[k]
                   + f_19 * fsk_249[k]
                   - f_14 * pc_z[k] * fsl1_312[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, pc_x, pc_y, pc_z, fsk_251, fsk_287, gsi0_335, \
                         gsi0_336, gsi1_335, gsi1_336, gsk_431, \
                         gsk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_17 * fsk_287[k]
                   + f_3 * pc_y[k] * gsk_431[k];

        t_539[k] = f_15 * fsk_251[k]
                   + f_1 * gsi0_335[k]
                   - f_2 * gsi1_335[k]
                   + f_3 * pc_z[k] * gsk_431[k];

        t_540[k] = f_1 * gsi0_336[k]
                   - f_2 * gsi1_336[k]
                   + f_3 * pc_x[k] * gsk_432[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_x, gsi0_337, gsi0_338, gsi0_339, gsi1_337, \
                         gsi1_338, gsi1_339, gsk_433, gsk_434, \
                         gsk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_20 * gsi0_337[k]
                   - f_21 * gsi1_337[k]
                   + f_3 * pc_x[k] * gsk_433[k];

        t_542[k] = f_20 * gsi0_338[k]
                   - f_21 * gsi1_338[k]
                   + f_3 * pc_x[k] * gsk_434[k];

        t_543[k] = f_12 * gsi0_339[k]
                   - f_13 * gsi1_339[k]
                   + f_3 * pc_x[k] * gsk_435[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pc_x, gsi0_340, gsi0_341, gsi0_342, gsi1_340, \
                         gsi1_341, gsi1_342, gsk_436, gsk_437, \
                         gsk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_12 * gsi0_340[k]
                   - f_13 * gsi1_340[k]
                   + f_3 * pc_x[k] * gsk_436[k];

        t_545[k] = f_12 * gsi0_341[k]
                   - f_13 * gsi1_341[k]
                   + f_3 * pc_x[k] * gsk_437[k];

        t_546[k] = f_10 * gsi0_342[k]
                   - f_11 * gsi1_342[k]
                   + f_3 * pc_x[k] * gsk_438[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pc_x, gsi0_343, gsi0_344, gsi0_345, gsi1_343, \
                         gsi1_344, gsi1_345, gsk_439, gsk_440, \
                         gsk_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_10 * gsi0_343[k]
                   - f_11 * gsi1_343[k]
                   + f_3 * pc_x[k] * gsk_439[k];

        t_548[k] = f_10 * gsi0_344[k]
                   - f_11 * gsi1_344[k]
                   + f_3 * pc_x[k] * gsk_440[k];

        t_549[k] = f_10 * gsi0_345[k]
                   - f_11 * gsi1_345[k]
                   + f_3 * pc_x[k] * gsk_441[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pc_x, gsi0_346, gsi0_347, gsi0_348, gsi1_346, \
                         gsi1_347, gsi1_348, gsk_442, gsk_443, \
                         gsk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_8 * gsi0_346[k]
                   - f_9 * gsi1_346[k]
                   + f_3 * pc_x[k] * gsk_442[k];

        t_551[k] = f_8 * gsi0_347[k]
                   - f_9 * gsi1_347[k]
                   + f_3 * pc_x[k] * gsk_443[k];

        t_552[k] = f_8 * gsi0_348[k]
                   - f_9 * gsi1_348[k]
                   + f_3 * pc_x[k] * gsk_444[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_x, gsi0_349, gsi0_350, gsi0_351, gsi1_349, \
                         gsi1_350, gsi1_351, gsk_445, gsk_446, \
                         gsk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_8 * gsi0_349[k]
                   - f_9 * gsi1_349[k]
                   + f_3 * pc_x[k] * gsk_445[k];

        t_554[k] = f_8 * gsi0_350[k]
                   - f_9 * gsi1_350[k]
                   + f_3 * pc_x[k] * gsk_446[k];

        t_555[k] = f_6 * gsi0_351[k]
                   - f_7 * gsi1_351[k]
                   + f_3 * pc_x[k] * gsk_447[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_x, gsi0_352, gsi0_353, gsi0_354, gsi1_352, \
                         gsi1_353, gsi1_354, gsk_448, gsk_449, \
                         gsk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_6 * gsi0_352[k]
                   - f_7 * gsi1_352[k]
                   + f_3 * pc_x[k] * gsk_448[k];

        t_557[k] = f_6 * gsi0_353[k]
                   - f_7 * gsi1_353[k]
                   + f_3 * pc_x[k] * gsk_449[k];

        t_558[k] = f_6 * gsi0_354[k]
                   - f_7 * gsi1_354[k]
                   + f_3 * pc_x[k] * gsk_450[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, pc_x, gsi0_355, gsi0_356, gsi0_357, gsi1_355, \
                         gsi1_356, gsi1_357, gsk_451, gsk_452, \
                         gsk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_6 * gsi0_355[k]
                   - f_7 * gsi1_355[k]
                   + f_3 * pc_x[k] * gsk_451[k];

        t_560[k] = f_6 * gsi0_356[k]
                   - f_7 * gsi1_356[k]
                   + f_3 * pc_x[k] * gsk_452[k];

        t_561[k] = f_4 * gsi0_357[k]
                   - f_5 * gsi1_357[k]
                   + f_3 * pc_x[k] * gsk_453[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, pc_x, gsi0_358, gsi0_359, gsi0_360, gsi1_358, \
                         gsi1_359, gsi1_360, gsk_454, gsk_455, \
                         gsk_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_4 * gsi0_358[k]
                   - f_5 * gsi1_358[k]
                   + f_3 * pc_x[k] * gsk_454[k];

        t_563[k] = f_4 * gsi0_359[k]
                   - f_5 * gsi1_359[k]
                   + f_3 * pc_x[k] * gsk_455[k];

        t_564[k] = f_4 * gsi0_360[k]
                   - f_5 * gsi1_360[k]
                   + f_3 * pc_x[k] * gsk_456[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pc_x, gsi0_361, gsi0_362, gsi0_363, \
                         gsi1_361, gsi1_362, gsi1_363, gsk_457, gsk_458, gsk_459, \
                         gsk_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_4 * gsi0_361[k]
                   - f_5 * gsi1_361[k]
                   + f_3 * pc_x[k] * gsk_457[k];

        t_566[k] = f_4 * gsi0_362[k]
                   - f_5 * gsi1_362[k]
                   + f_3 * pc_x[k] * gsk_458[k];

        t_567[k] = f_4 * gsi0_363[k]
                   - f_5 * gsi1_363[k]
                   + f_3 * pc_x[k] * gsk_459[k];

        t_568[k] = f_3 * pc_x[k] * gsk_460[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, t_574, t_575, pc_x, gsk_461, \
                         gsk_462, gsk_463, gsk_464, gsk_465, gsk_466, \
                         gsk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_3 * pc_x[k] * gsk_461[k];

        t_570[k] = f_3 * pc_x[k] * gsk_462[k];

        t_571[k] = f_3 * pc_x[k] * gsk_463[k];

        t_572[k] = f_3 * pc_x[k] * gsk_464[k];

        t_573[k] = f_3 * pc_x[k] * gsk_465[k];

        t_574[k] = f_3 * pc_x[k] * gsk_466[k];

        t_575[k] = f_3 * pc_x[k] * gsk_467[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, pc_y, pc_z, fsk_280, fsk_316, fsk_318, gsi0_357, \
                         gsi0_359, gsi1_357, gsi1_359, gsk_460, \
                         gsk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_16 * fsk_316[k]
                   + f_1 * gsi0_357[k]
                   - f_2 * gsi1_357[k]
                   + f_3 * pc_y[k] * gsk_460[k];

        t_577[k] = f_16 * fsk_280[k]
                   + f_3 * pc_z[k] * gsk_460[k];

        t_578[k] = f_16 * fsk_318[k]
                   + f_12 * gsi0_359[k]
                   - f_13 * gsi1_359[k]
                   + f_3 * pc_y[k] * gsk_462[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pc_y, fsk_319, fsk_320, fsk_321, gsi0_360, \
                         gsi0_361, gsi0_362, gsi1_360, gsi1_361, gsi1_362, gsk_463, gsk_464, \
                         gsk_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_16 * fsk_319[k]
                   + f_10 * gsi0_360[k]
                   - f_11 * gsi1_360[k]
                   + f_3 * pc_y[k] * gsk_463[k];

        t_580[k] = f_16 * fsk_320[k]
                   + f_8 * gsi0_361[k]
                   - f_9 * gsi1_361[k]
                   + f_3 * pc_y[k] * gsk_464[k];

        t_581[k] = f_16 * fsk_321[k]
                   + f_6 * gsi0_362[k]
                   - f_7 * gsi1_362[k]
                   + f_3 * pc_y[k] * gsk_465[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pa_y, pc_y, pc_z, fsl0_405, fsk_287, \
                         fsk_322, fsk_323, fsl1_405, gsi0_363, gsi1_363, gsk_466, \
                         gsk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_16 * fsk_322[k]
                   + f_4 * gsi0_363[k]
                   - f_5 * gsi1_363[k]
                   + f_3 * pc_y[k] * gsk_466[k];

        t_583[k] = f_16 * fsk_323[k]
                   + f_3 * pc_y[k] * gsk_467[k];

        t_584[k] = f_16 * fsk_287[k]
                   + f_1 * gsi0_363[k]
                   - f_2 * gsi1_363[k]
                   + f_3 * pc_z[k] * gsk_467[k];

        t_585[k] = pa_y[k] * fsl0_405[k]
                   - f_14 * pc_y[k] * fsl1_405[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, pa_y, pc_x, pc_y, fsl0_407, fsl1_407, gsi0_365, \
                         gsi0_367, gsi1_365, gsi1_367, gsk_469, \
                         gsk_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_20 * gsi0_365[k]
                   - f_21 * gsi1_365[k]
                   + f_3 * pc_x[k] * gsk_469[k];

        t_587[k] = pa_y[k] * fsl0_407[k]
                   - f_14 * pc_y[k] * fsl1_407[k];

        t_588[k] = f_12 * gsi0_367[k]
                   - f_13 * gsi1_367[k]
                   + f_3 * pc_x[k] * gsk_471[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, pa_y, pc_x, pc_y, fsl0_410, fsl1_410, gsi0_368, \
                         gsi0_370, gsi1_368, gsi1_370, gsk_472, \
                         gsk_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_12 * gsi0_368[k]
                   - f_13 * gsi1_368[k]
                   + f_3 * pc_x[k] * gsk_472[k];

        t_590[k] = pa_y[k] * fsl0_410[k]
                   - f_14 * pc_y[k] * fsl1_410[k];

        t_591[k] = f_10 * gsi0_370[k]
                   - f_11 * gsi1_370[k]
                   + f_3 * pc_x[k] * gsk_474[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pa_y, pc_x, pc_y, fsl0_414, fsl1_414, gsi0_371, \
                         gsi0_372, gsi1_371, gsi1_372, gsk_475, \
                         gsk_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_10 * gsi0_371[k]
                   - f_11 * gsi1_371[k]
                   + f_3 * pc_x[k] * gsk_475[k];

        t_593[k] = f_10 * gsi0_372[k]
                   - f_11 * gsi1_372[k]
                   + f_3 * pc_x[k] * gsk_476[k];

        t_594[k] = pa_y[k] * fsl0_414[k]
                   - f_14 * pc_y[k] * fsl1_414[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_x, gsi0_374, gsi0_375, gsi0_376, gsi1_374, \
                         gsi1_375, gsi1_376, gsk_478, gsk_479, \
                         gsk_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_8 * gsi0_374[k]
                   - f_9 * gsi1_374[k]
                   + f_3 * pc_x[k] * gsk_478[k];

        t_596[k] = f_8 * gsi0_375[k]
                   - f_9 * gsi1_375[k]
                   + f_3 * pc_x[k] * gsk_479[k];

        t_597[k] = f_8 * gsi0_376[k]
                   - f_9 * gsi1_376[k]
                   + f_3 * pc_x[k] * gsk_480[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pa_y, pc_x, pc_y, fsl0_419, fsl1_419, gsi0_377, \
                         gsi0_379, gsi1_377, gsi1_379, gsk_481, \
                         gsk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_8 * gsi0_377[k]
                   - f_9 * gsi1_377[k]
                   + f_3 * pc_x[k] * gsk_481[k];

        t_599[k] = pa_y[k] * fsl0_419[k]
                   - f_14 * pc_y[k] * fsl1_419[k];

        t_600[k] = f_6 * gsi0_379[k]
                   - f_7 * gsi1_379[k]
                   + f_3 * pc_x[k] * gsk_483[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, pc_x, gsi0_380, gsi0_381, gsi0_382, gsi1_380, \
                         gsi1_381, gsi1_382, gsk_484, gsk_485, \
                         gsk_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_6 * gsi0_380[k]
                   - f_7 * gsi1_380[k]
                   + f_3 * pc_x[k] * gsk_484[k];

        t_602[k] = f_6 * gsi0_381[k]
                   - f_7 * gsi1_381[k]
                   + f_3 * pc_x[k] * gsk_485[k];

        t_603[k] = f_6 * gsi0_382[k]
                   - f_7 * gsi1_382[k]
                   + f_3 * pc_x[k] * gsk_486[k];
    }
}

static auto
compute_prim_gsl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsl0,
                                                          const size_t fsk, const size_t fsl1,
                                                          const size_t gsi0, const size_t gsi1,
                                                          const size_t gsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsl0_425 = buffer.data(fsl0 + 425);
    const auto *fsl0_432 = buffer.data(fsl0 + 432);
    const auto *fsl0_441 = buffer.data(fsl0 + 441);
    const auto *fsl0_443 = buffer.data(fsl0 + 443);
    const auto *fsl0_444 = buffer.data(fsl0 + 444);
    const auto *fsl0_445 = buffer.data(fsl0 + 445);
    const auto *fsl0_446 = buffer.data(fsl0 + 446);
    const auto *fsl0_447 = buffer.data(fsl0 + 447);
    const auto *fsl0_449 = buffer.data(fsl0 + 449);

    const auto *fsk_316 = buffer.data(fsk + 316);
    const auto *fsk_352 = buffer.data(fsk + 352);
    const auto *fsk_354 = buffer.data(fsk + 354);
    const auto *fsk_355 = buffer.data(fsk + 355);
    const auto *fsk_356 = buffer.data(fsk + 356);
    const auto *fsk_357 = buffer.data(fsk + 357);
    const auto *fsk_358 = buffer.data(fsk + 358);
    const auto *fsk_359 = buffer.data(fsk + 359);

    const auto *fsl1_425 = buffer.data(fsl1 + 425);
    const auto *fsl1_432 = buffer.data(fsl1 + 432);
    const auto *fsl1_441 = buffer.data(fsl1 + 441);
    const auto *fsl1_443 = buffer.data(fsl1 + 443);
    const auto *fsl1_444 = buffer.data(fsl1 + 444);
    const auto *fsl1_445 = buffer.data(fsl1 + 445);
    const auto *fsl1_446 = buffer.data(fsl1 + 446);
    const auto *fsl1_447 = buffer.data(fsl1 + 447);
    const auto *fsl1_449 = buffer.data(fsl1 + 449);

    const auto *gsi0_383 = buffer.data(gsi0 + 383);
    const auto *gsi0_385 = buffer.data(gsi0 + 385);
    const auto *gsi0_386 = buffer.data(gsi0 + 386);
    const auto *gsi0_387 = buffer.data(gsi0 + 387);
    const auto *gsi0_388 = buffer.data(gsi0 + 388);
    const auto *gsi0_389 = buffer.data(gsi0 + 389);
    const auto *gsi0_390 = buffer.data(gsi0 + 390);
    const auto *gsi0_392 = buffer.data(gsi0 + 392);
    const auto *gsi0_394 = buffer.data(gsi0 + 394);
    const auto *gsi0_395 = buffer.data(gsi0 + 395);
    const auto *gsi0_397 = buffer.data(gsi0 + 397);
    const auto *gsi0_398 = buffer.data(gsi0 + 398);
    const auto *gsi0_399 = buffer.data(gsi0 + 399);
    const auto *gsi0_401 = buffer.data(gsi0 + 401);
    const auto *gsi0_402 = buffer.data(gsi0 + 402);
    const auto *gsi0_403 = buffer.data(gsi0 + 403);
    const auto *gsi0_404 = buffer.data(gsi0 + 404);
    const auto *gsi0_406 = buffer.data(gsi0 + 406);
    const auto *gsi0_407 = buffer.data(gsi0 + 407);
    const auto *gsi0_408 = buffer.data(gsi0 + 408);
    const auto *gsi0_409 = buffer.data(gsi0 + 409);
    const auto *gsi0_410 = buffer.data(gsi0 + 410);
    const auto *gsi0_412 = buffer.data(gsi0 + 412);
    const auto *gsi0_413 = buffer.data(gsi0 + 413);
    const auto *gsi0_414 = buffer.data(gsi0 + 414);
    const auto *gsi0_415 = buffer.data(gsi0 + 415);
    const auto *gsi0_416 = buffer.data(gsi0 + 416);
    const auto *gsi0_417 = buffer.data(gsi0 + 417);
    const auto *gsi0_418 = buffer.data(gsi0 + 418);
    const auto *gsi0_419 = buffer.data(gsi0 + 419);

    const auto *gsi1_383 = buffer.data(gsi1 + 383);
    const auto *gsi1_385 = buffer.data(gsi1 + 385);
    const auto *gsi1_386 = buffer.data(gsi1 + 386);
    const auto *gsi1_387 = buffer.data(gsi1 + 387);
    const auto *gsi1_388 = buffer.data(gsi1 + 388);
    const auto *gsi1_389 = buffer.data(gsi1 + 389);
    const auto *gsi1_390 = buffer.data(gsi1 + 390);
    const auto *gsi1_392 = buffer.data(gsi1 + 392);
    const auto *gsi1_394 = buffer.data(gsi1 + 394);
    const auto *gsi1_395 = buffer.data(gsi1 + 395);
    const auto *gsi1_397 = buffer.data(gsi1 + 397);
    const auto *gsi1_398 = buffer.data(gsi1 + 398);
    const auto *gsi1_399 = buffer.data(gsi1 + 399);
    const auto *gsi1_401 = buffer.data(gsi1 + 401);
    const auto *gsi1_402 = buffer.data(gsi1 + 402);
    const auto *gsi1_403 = buffer.data(gsi1 + 403);
    const auto *gsi1_404 = buffer.data(gsi1 + 404);
    const auto *gsi1_406 = buffer.data(gsi1 + 406);
    const auto *gsi1_407 = buffer.data(gsi1 + 407);
    const auto *gsi1_408 = buffer.data(gsi1 + 408);
    const auto *gsi1_409 = buffer.data(gsi1 + 409);
    const auto *gsi1_410 = buffer.data(gsi1 + 410);
    const auto *gsi1_412 = buffer.data(gsi1 + 412);
    const auto *gsi1_413 = buffer.data(gsi1 + 413);
    const auto *gsi1_414 = buffer.data(gsi1 + 414);
    const auto *gsi1_415 = buffer.data(gsi1 + 415);
    const auto *gsi1_416 = buffer.data(gsi1 + 416);
    const auto *gsi1_417 = buffer.data(gsi1 + 417);
    const auto *gsi1_418 = buffer.data(gsi1 + 418);
    const auto *gsi1_419 = buffer.data(gsi1 + 419);

    const auto *gsk_487 = buffer.data(gsk + 487);
    const auto *gsk_489 = buffer.data(gsk + 489);
    const auto *gsk_490 = buffer.data(gsk + 490);
    const auto *gsk_491 = buffer.data(gsk + 491);
    const auto *gsk_492 = buffer.data(gsk + 492);
    const auto *gsk_493 = buffer.data(gsk + 493);
    const auto *gsk_494 = buffer.data(gsk + 494);
    const auto *gsk_496 = buffer.data(gsk + 496);
    const auto *gsk_497 = buffer.data(gsk + 497);
    const auto *gsk_498 = buffer.data(gsk + 498);
    const auto *gsk_499 = buffer.data(gsk + 499);
    const auto *gsk_500 = buffer.data(gsk + 500);
    const auto *gsk_501 = buffer.data(gsk + 501);
    const auto *gsk_502 = buffer.data(gsk + 502);
    const auto *gsk_503 = buffer.data(gsk + 503);
    const auto *gsk_504 = buffer.data(gsk + 504);
    const auto *gsk_506 = buffer.data(gsk + 506);
    const auto *gsk_507 = buffer.data(gsk + 507);
    const auto *gsk_509 = buffer.data(gsk + 509);
    const auto *gsk_510 = buffer.data(gsk + 510);
    const auto *gsk_511 = buffer.data(gsk + 511);
    const auto *gsk_513 = buffer.data(gsk + 513);
    const auto *gsk_514 = buffer.data(gsk + 514);
    const auto *gsk_515 = buffer.data(gsk + 515);
    const auto *gsk_516 = buffer.data(gsk + 516);
    const auto *gsk_518 = buffer.data(gsk + 518);
    const auto *gsk_519 = buffer.data(gsk + 519);
    const auto *gsk_520 = buffer.data(gsk + 520);
    const auto *gsk_521 = buffer.data(gsk + 521);
    const auto *gsk_522 = buffer.data(gsk + 522);
    const auto *gsk_524 = buffer.data(gsk + 524);
    const auto *gsk_525 = buffer.data(gsk + 525);
    const auto *gsk_526 = buffer.data(gsk + 526);
    const auto *gsk_527 = buffer.data(gsk + 527);
    const auto *gsk_528 = buffer.data(gsk + 528);
    const auto *gsk_529 = buffer.data(gsk + 529);
    const auto *gsk_531 = buffer.data(gsk + 531);
    const auto *gsk_532 = buffer.data(gsk + 532);
    const auto *gsk_533 = buffer.data(gsk + 533);
    const auto *gsk_534 = buffer.data(gsk + 534);
    const auto *gsk_535 = buffer.data(gsk + 535);
    const auto *gsk_536 = buffer.data(gsk + 536);
    const auto *gsk_537 = buffer.data(gsk + 537);
    const auto *gsk_538 = buffer.data(gsk + 538);
    const auto *gsk_539 = buffer.data(gsk + 539);

#pragma omp simd aligned(t_604, t_605, t_606, pa_y, pc_x, pc_y, fsl0_425, fsl1_425, gsi0_383, \
                         gsi0_385, gsi1_383, gsi1_385, gsk_487, \
                         gsk_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_6 * gsi0_383[k]
                   - f_7 * gsi1_383[k]
                   + f_3 * pc_x[k] * gsk_487[k];

        t_605[k] = pa_y[k] * fsl0_425[k]
                   - f_14 * pc_y[k] * fsl1_425[k];

        t_606[k] = f_4 * gsi0_385[k]
                   - f_5 * gsi1_385[k]
                   + f_3 * pc_x[k] * gsk_489[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_x, gsi0_386, gsi0_387, gsi0_388, gsi1_386, \
                         gsi1_387, gsi1_388, gsk_490, gsk_491, \
                         gsk_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_4 * gsi0_386[k]
                   - f_5 * gsi1_386[k]
                   + f_3 * pc_x[k] * gsk_490[k];

        t_608[k] = f_4 * gsi0_387[k]
                   - f_5 * gsi1_387[k]
                   + f_3 * pc_x[k] * gsk_491[k];

        t_609[k] = f_4 * gsi0_388[k]
                   - f_5 * gsi1_388[k]
                   + f_3 * pc_x[k] * gsk_492[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_y, pc_x, pc_y, fsl0_432, fsl1_432, \
                         gsi0_389, gsi0_390, gsi1_389, gsi1_390, gsk_493, gsk_494, \
                         gsk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_4 * gsi0_389[k]
                   - f_5 * gsi1_389[k]
                   + f_3 * pc_x[k] * gsk_493[k];

        t_611[k] = f_4 * gsi0_390[k]
                   - f_5 * gsi1_390[k]
                   + f_3 * pc_x[k] * gsk_494[k];

        t_612[k] = pa_y[k] * fsl0_432[k]
                   - f_14 * pc_y[k] * fsl1_432[k];

        t_613[k] = f_3 * pc_x[k] * gsk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, t_619, t_620, pc_x, gsk_497, \
                         gsk_498, gsk_499, gsk_500, gsk_501, gsk_502, \
                         gsk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_3 * pc_x[k] * gsk_497[k];

        t_615[k] = f_3 * pc_x[k] * gsk_498[k];

        t_616[k] = f_3 * pc_x[k] * gsk_499[k];

        t_617[k] = f_3 * pc_x[k] * gsk_500[k];

        t_618[k] = f_3 * pc_x[k] * gsk_501[k];

        t_619[k] = f_3 * pc_x[k] * gsk_502[k];

        t_620[k] = f_3 * pc_x[k] * gsk_503[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_y, pc_y, pc_z, fsl0_441, fsl0_443, fsk_316, \
                         fsk_352, fsk_354, fsl1_441, fsl1_443, \
                         gsk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = pa_y[k] * fsl0_441[k]
                   + f_22 * fsk_352[k]
                   - f_14 * pc_y[k] * fsl1_441[k];

        t_622[k] = f_17 * fsk_316[k]
                   + f_3 * pc_z[k] * gsk_496[k];

        t_623[k] = pa_y[k] * fsl0_443[k]
                   + f_19 * fsk_354[k]
                   - f_14 * pc_y[k] * fsl1_443[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_y, pc_y, fsl0_444, fsl0_445, fsl0_446, \
                         fsk_355, fsk_356, fsk_357, fsl1_444, fsl1_445, \
                         fsl1_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pa_y[k] * fsl0_444[k]
                   + f_18 * fsk_355[k]
                   - f_14 * pc_y[k] * fsl1_444[k];

        t_625[k] = pa_y[k] * fsl0_445[k]
                   + f_0 * fsk_356[k]
                   - f_14 * pc_y[k] * fsl1_445[k];

        t_626[k] = pa_y[k] * fsl0_446[k]
                   + f_17 * fsk_357[k]
                   - f_14 * pc_y[k] * fsl1_446[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_y, pc_y, fsl0_447, fsl0_449, fsk_358, \
                         fsk_359, fsl1_447, fsl1_449, gsk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = pa_y[k] * fsl0_447[k]
                   + f_16 * fsk_358[k]
                   - f_14 * pc_y[k] * fsl1_447[k];

        t_628[k] = f_15 * fsk_359[k]
                   + f_3 * pc_y[k] * gsk_503[k];

        t_629[k] = pa_y[k] * fsl0_449[k]
                   - f_14 * pc_y[k] * fsl1_449[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, pc_x, pc_y, gsi0_392, gsi0_394, \
                         gsi0_395, gsi1_392, gsi1_394, gsi1_395, gsk_504, gsk_506, \
                         gsk_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_1 * gsi0_392[k]
                   - f_2 * gsi1_392[k]
                   + f_3 * pc_x[k] * gsk_504[k];

        t_631[k] = f_3 * pc_y[k] * gsk_504[k];

        t_632[k] = f_20 * gsi0_394[k]
                   - f_21 * gsi1_394[k]
                   + f_3 * pc_x[k] * gsk_506[k];

        t_633[k] = f_12 * gsi0_395[k]
                   - f_13 * gsi1_395[k]
                   + f_3 * pc_x[k] * gsk_507[k];

        t_634[k] = f_3 * pc_y[k] * gsk_506[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pc_x, pc_y, gsi0_397, gsi0_398, gsi0_399, \
                         gsi1_397, gsi1_398, gsi1_399, gsk_509, gsk_510, \
                         gsk_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_12 * gsi0_397[k]
                   - f_13 * gsi1_397[k]
                   + f_3 * pc_x[k] * gsk_509[k];

        t_636[k] = f_10 * gsi0_398[k]
                   - f_11 * gsi1_398[k]
                   + f_3 * pc_x[k] * gsk_510[k];

        t_637[k] = f_10 * gsi0_399[k]
                   - f_11 * gsi1_399[k]
                   + f_3 * pc_x[k] * gsk_511[k];

        t_638[k] = f_3 * pc_y[k] * gsk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, gsi0_401, gsi0_402, gsi0_403, gsi1_401, \
                         gsi1_402, gsi1_403, gsk_513, gsk_514, \
                         gsk_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_10 * gsi0_401[k]
                   - f_11 * gsi1_401[k]
                   + f_3 * pc_x[k] * gsk_513[k];

        t_640[k] = f_8 * gsi0_402[k]
                   - f_9 * gsi1_402[k]
                   + f_3 * pc_x[k] * gsk_514[k];

        t_641[k] = f_8 * gsi0_403[k]
                   - f_9 * gsi1_403[k]
                   + f_3 * pc_x[k] * gsk_515[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, gsi0_404, gsi0_406, gsi0_407, \
                         gsi1_404, gsi1_406, gsi1_407, gsk_513, gsk_516, gsk_518, \
                         gsk_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_8 * gsi0_404[k]
                   - f_9 * gsi1_404[k]
                   + f_3 * pc_x[k] * gsk_516[k];

        t_643[k] = f_3 * pc_y[k] * gsk_513[k];

        t_644[k] = f_8 * gsi0_406[k]
                   - f_9 * gsi1_406[k]
                   + f_3 * pc_x[k] * gsk_518[k];

        t_645[k] = f_6 * gsi0_407[k]
                   - f_7 * gsi1_407[k]
                   + f_3 * pc_x[k] * gsk_519[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pc_x, pc_y, gsi0_408, gsi0_409, gsi0_410, \
                         gsi1_408, gsi1_409, gsi1_410, gsk_518, gsk_520, gsk_521, \
                         gsk_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_6 * gsi0_408[k]
                   - f_7 * gsi1_408[k]
                   + f_3 * pc_x[k] * gsk_520[k];

        t_647[k] = f_6 * gsi0_409[k]
                   - f_7 * gsi1_409[k]
                   + f_3 * pc_x[k] * gsk_521[k];

        t_648[k] = f_6 * gsi0_410[k]
                   - f_7 * gsi1_410[k]
                   + f_3 * pc_x[k] * gsk_522[k];

        t_649[k] = f_3 * pc_y[k] * gsk_518[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, gsi0_412, gsi0_413, gsi0_414, gsi1_412, \
                         gsi1_413, gsi1_414, gsk_524, gsk_525, \
                         gsk_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_6 * gsi0_412[k]
                   - f_7 * gsi1_412[k]
                   + f_3 * pc_x[k] * gsk_524[k];

        t_651[k] = f_4 * gsi0_413[k]
                   - f_5 * gsi1_413[k]
                   + f_3 * pc_x[k] * gsk_525[k];

        t_652[k] = f_4 * gsi0_414[k]
                   - f_5 * gsi1_414[k]
                   + f_3 * pc_x[k] * gsk_526[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, pc_x, pc_y, gsi0_415, gsi0_416, gsi0_417, \
                         gsi1_415, gsi1_416, gsi1_417, gsk_524, gsk_527, gsk_528, \
                         gsk_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_4 * gsi0_415[k]
                   - f_5 * gsi1_415[k]
                   + f_3 * pc_x[k] * gsk_527[k];

        t_654[k] = f_4 * gsi0_416[k]
                   - f_5 * gsi1_416[k]
                   + f_3 * pc_x[k] * gsk_528[k];

        t_655[k] = f_4 * gsi0_417[k]
                   - f_5 * gsi1_417[k]
                   + f_3 * pc_x[k] * gsk_529[k];

        t_656[k] = f_3 * pc_y[k] * gsk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, t_662, pc_x, gsi0_419, gsi1_419, \
                         gsk_531, gsk_532, gsk_533, gsk_534, gsk_535, \
                         gsk_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_4 * gsi0_419[k]
                   - f_5 * gsi1_419[k]
                   + f_3 * pc_x[k] * gsk_531[k];

        t_658[k] = f_3 * pc_x[k] * gsk_532[k];

        t_659[k] = f_3 * pc_x[k] * gsk_533[k];

        t_660[k] = f_3 * pc_x[k] * gsk_534[k];

        t_661[k] = f_3 * pc_x[k] * gsk_535[k];

        t_662[k] = f_3 * pc_x[k] * gsk_536[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, t_667, pc_x, pc_y, gsi0_413, gsi0_414, \
                         gsi1_413, gsi1_414, gsk_532, gsk_533, gsk_537, gsk_538, \
                         gsk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_3 * pc_x[k] * gsk_537[k];

        t_664[k] = f_3 * pc_x[k] * gsk_538[k];

        t_665[k] = f_3 * pc_x[k] * gsk_539[k];

        t_666[k] = f_1 * gsi0_413[k]
                   - f_2 * gsi1_413[k]
                   + f_3 * pc_y[k] * gsk_532[k];

        t_667[k] = f_20 * gsi0_414[k]
                   - f_21 * gsi1_414[k]
                   + f_3 * pc_y[k] * gsk_533[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, pc_y, gsi0_415, gsi0_416, gsi0_417, gsi1_415, \
                         gsi1_416, gsi1_417, gsk_534, gsk_535, \
                         gsk_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_12 * gsi0_415[k]
                   - f_13 * gsi1_415[k]
                   + f_3 * pc_y[k] * gsk_534[k];

        t_669[k] = f_10 * gsi0_416[k]
                   - f_11 * gsi1_416[k]
                   + f_3 * pc_y[k] * gsk_535[k];

        t_670[k] = f_8 * gsi0_417[k]
                   - f_9 * gsi1_417[k]
                   + f_3 * pc_y[k] * gsk_536[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, pc_y, pc_z, fsk_359, gsi0_418, gsi0_419, \
                         gsi1_418, gsi1_419, gsk_537, gsk_538, \
                         gsk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = f_6 * gsi0_418[k]
                   - f_7 * gsi1_418[k]
                   + f_3 * pc_y[k] * gsk_537[k];

        t_672[k] = f_4 * gsi0_419[k]
                   - f_5 * gsi1_419[k]
                   + f_3 * pc_y[k] * gsk_538[k];

        t_673[k] = f_3 * pc_y[k] * gsk_539[k];

        t_674[k] = f_0 * fsk_359[k]
                   + f_1 * gsi0_419[k]
                   - f_2 * gsi1_419[k]
                   + f_3 * pc_z[k] * gsk_539[k];
    }
}

auto
compute_prim_gsl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsl0, const size_t fsk,
                                                   const size_t fsl1, const size_t gsi0,
                                                   const size_t gsi1, const size_t gsk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gsl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, fsl0, fsk,
                                                              fsl1, gsi0, gsi1, gsk, ncols,
                                                              gamma, p, q);

    compute_prim_gsl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, fsl0, fsk,
                                                              fsl1, gsi0, gsi1, gsk, ncols,
                                                              gamma, p, q);

    compute_prim_gsl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, fsl0, fsk,
                                                              fsl1, gsi0, gsi1, gsk, ncols,
                                                              gamma, p, q);

    compute_prim_gsl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, fsl0, fsk,
                                                              fsl1, gsi0, gsi1, gsk, ncols,
                                                              gamma, p, q);

    compute_prim_gsl_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, fsl0, fsk,
                                                              fsl1, gsi0, gsi1, gsk, ncols,
                                                              gamma, p, q);

    compute_prim_gsl_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, fsl0, fsk,
                                                              fsl1, gsi0, gsi1, gsk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
