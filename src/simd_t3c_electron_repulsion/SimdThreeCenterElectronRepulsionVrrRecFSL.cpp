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


#include "SimdThreeCenterElectronRepulsionVrrRecFSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fsl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsl0,
                                                          const size_t dsk, const size_t dsl1,
                                                          const size_t fsi0, const size_t fsi1,
                                                          const size_t fsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
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

    const auto *dsl0_0 = buffer.data(dsl0 + 0);
    const auto *dsl0_3 = buffer.data(dsl0 + 3);
    const auto *dsl0_5 = buffer.data(dsl0 + 5);
    const auto *dsl0_6 = buffer.data(dsl0 + 6);
    const auto *dsl0_9 = buffer.data(dsl0 + 9);
    const auto *dsl0_10 = buffer.data(dsl0 + 10);
    const auto *dsl0_14 = buffer.data(dsl0 + 14);
    const auto *dsl0_15 = buffer.data(dsl0 + 15);
    const auto *dsl0_20 = buffer.data(dsl0 + 20);
    const auto *dsl0_21 = buffer.data(dsl0 + 21);
    const auto *dsl0_27 = buffer.data(dsl0 + 27);
    const auto *dsl0_36 = buffer.data(dsl0 + 36);
    const auto *dsl0_44 = buffer.data(dsl0 + 44);

    const auto *dsk_0 = buffer.data(dsk + 0);
    const auto *dsk_1 = buffer.data(dsk + 1);
    const auto *dsk_2 = buffer.data(dsk + 2);
    const auto *dsk_3 = buffer.data(dsk + 3);
    const auto *dsk_5 = buffer.data(dsk + 5);
    const auto *dsk_6 = buffer.data(dsk + 6);
    const auto *dsk_9 = buffer.data(dsk + 9);
    const auto *dsk_10 = buffer.data(dsk + 10);
    const auto *dsk_14 = buffer.data(dsk + 14);
    const auto *dsk_15 = buffer.data(dsk + 15);
    const auto *dsk_20 = buffer.data(dsk + 20);
    const auto *dsk_28 = buffer.data(dsk + 28);
    const auto *dsk_30 = buffer.data(dsk + 30);
    const auto *dsk_31 = buffer.data(dsk + 31);
    const auto *dsk_32 = buffer.data(dsk + 32);
    const auto *dsk_33 = buffer.data(dsk + 33);
    const auto *dsk_35 = buffer.data(dsk + 35);
    const auto *dsk_64 = buffer.data(dsk + 64);
    const auto *dsk_66 = buffer.data(dsk + 66);
    const auto *dsk_67 = buffer.data(dsk + 67);
    const auto *dsk_68 = buffer.data(dsk + 68);
    const auto *dsk_69 = buffer.data(dsk + 69);
    const auto *dsk_70 = buffer.data(dsk + 70);
    const auto *dsk_71 = buffer.data(dsk + 71);
    const auto *dsk_100 = buffer.data(dsk + 100);
    const auto *dsk_101 = buffer.data(dsk + 101);
    const auto *dsk_102 = buffer.data(dsk + 102);
    const auto *dsk_103 = buffer.data(dsk + 103);
    const auto *dsk_104 = buffer.data(dsk + 104);
    const auto *dsk_105 = buffer.data(dsk + 105);
    const auto *dsk_107 = buffer.data(dsk + 107);

    const auto *dsl1_0 = buffer.data(dsl1 + 0);
    const auto *dsl1_3 = buffer.data(dsl1 + 3);
    const auto *dsl1_5 = buffer.data(dsl1 + 5);
    const auto *dsl1_6 = buffer.data(dsl1 + 6);
    const auto *dsl1_9 = buffer.data(dsl1 + 9);
    const auto *dsl1_10 = buffer.data(dsl1 + 10);
    const auto *dsl1_14 = buffer.data(dsl1 + 14);
    const auto *dsl1_15 = buffer.data(dsl1 + 15);
    const auto *dsl1_20 = buffer.data(dsl1 + 20);
    const auto *dsl1_21 = buffer.data(dsl1 + 21);
    const auto *dsl1_27 = buffer.data(dsl1 + 27);
    const auto *dsl1_36 = buffer.data(dsl1 + 36);
    const auto *dsl1_44 = buffer.data(dsl1 + 44);

    const auto *fsi0_0 = buffer.data(fsi0 + 0);
    const auto *fsi0_1 = buffer.data(fsi0 + 1);
    const auto *fsi0_2 = buffer.data(fsi0 + 2);
    const auto *fsi0_3 = buffer.data(fsi0 + 3);
    const auto *fsi0_5 = buffer.data(fsi0 + 5);
    const auto *fsi0_6 = buffer.data(fsi0 + 6);
    const auto *fsi0_8 = buffer.data(fsi0 + 8);
    const auto *fsi0_9 = buffer.data(fsi0 + 9);
    const auto *fsi0_10 = buffer.data(fsi0 + 10);
    const auto *fsi0_12 = buffer.data(fsi0 + 12);
    const auto *fsi0_13 = buffer.data(fsi0 + 13);
    const auto *fsi0_14 = buffer.data(fsi0 + 14);
    const auto *fsi0_21 = buffer.data(fsi0 + 21);
    const auto *fsi0_23 = buffer.data(fsi0 + 23);
    const auto *fsi0_24 = buffer.data(fsi0 + 24);
    const auto *fsi0_25 = buffer.data(fsi0 + 25);
    const auto *fsi0_26 = buffer.data(fsi0 + 26);
    const auto *fsi0_27 = buffer.data(fsi0 + 27);
    const auto *fsi0_31 = buffer.data(fsi0 + 31);
    const auto *fsi0_34 = buffer.data(fsi0 + 34);
    const auto *fsi0_35 = buffer.data(fsi0 + 35);
    const auto *fsi0_38 = buffer.data(fsi0 + 38);
    const auto *fsi0_39 = buffer.data(fsi0 + 39);
    const auto *fsi0_40 = buffer.data(fsi0 + 40);
    const auto *fsi0_49 = buffer.data(fsi0 + 49);
    const auto *fsi0_50 = buffer.data(fsi0 + 50);
    const auto *fsi0_51 = buffer.data(fsi0 + 51);
    const auto *fsi0_52 = buffer.data(fsi0 + 52);
    const auto *fsi0_53 = buffer.data(fsi0 + 53);
    const auto *fsi0_58 = buffer.data(fsi0 + 58);
    const auto *fsi0_60 = buffer.data(fsi0 + 60);
    const auto *fsi0_61 = buffer.data(fsi0 + 61);
    const auto *fsi0_63 = buffer.data(fsi0 + 63);
    const auto *fsi0_64 = buffer.data(fsi0 + 64);
    const auto *fsi0_65 = buffer.data(fsi0 + 65);
    const auto *fsi0_67 = buffer.data(fsi0 + 67);
    const auto *fsi0_68 = buffer.data(fsi0 + 68);
    const auto *fsi0_69 = buffer.data(fsi0 + 69);
    const auto *fsi0_70 = buffer.data(fsi0 + 70);
    const auto *fsi0_78 = buffer.data(fsi0 + 78);
    const auto *fsi0_79 = buffer.data(fsi0 + 79);

    const auto *fsi1_0 = buffer.data(fsi1 + 0);
    const auto *fsi1_1 = buffer.data(fsi1 + 1);
    const auto *fsi1_2 = buffer.data(fsi1 + 2);
    const auto *fsi1_3 = buffer.data(fsi1 + 3);
    const auto *fsi1_5 = buffer.data(fsi1 + 5);
    const auto *fsi1_6 = buffer.data(fsi1 + 6);
    const auto *fsi1_8 = buffer.data(fsi1 + 8);
    const auto *fsi1_9 = buffer.data(fsi1 + 9);
    const auto *fsi1_10 = buffer.data(fsi1 + 10);
    const auto *fsi1_12 = buffer.data(fsi1 + 12);
    const auto *fsi1_13 = buffer.data(fsi1 + 13);
    const auto *fsi1_14 = buffer.data(fsi1 + 14);
    const auto *fsi1_21 = buffer.data(fsi1 + 21);
    const auto *fsi1_23 = buffer.data(fsi1 + 23);
    const auto *fsi1_24 = buffer.data(fsi1 + 24);
    const auto *fsi1_25 = buffer.data(fsi1 + 25);
    const auto *fsi1_26 = buffer.data(fsi1 + 26);
    const auto *fsi1_27 = buffer.data(fsi1 + 27);
    const auto *fsi1_31 = buffer.data(fsi1 + 31);
    const auto *fsi1_34 = buffer.data(fsi1 + 34);
    const auto *fsi1_35 = buffer.data(fsi1 + 35);
    const auto *fsi1_38 = buffer.data(fsi1 + 38);
    const auto *fsi1_39 = buffer.data(fsi1 + 39);
    const auto *fsi1_40 = buffer.data(fsi1 + 40);
    const auto *fsi1_49 = buffer.data(fsi1 + 49);
    const auto *fsi1_50 = buffer.data(fsi1 + 50);
    const auto *fsi1_51 = buffer.data(fsi1 + 51);
    const auto *fsi1_52 = buffer.data(fsi1 + 52);
    const auto *fsi1_53 = buffer.data(fsi1 + 53);
    const auto *fsi1_58 = buffer.data(fsi1 + 58);
    const auto *fsi1_60 = buffer.data(fsi1 + 60);
    const auto *fsi1_61 = buffer.data(fsi1 + 61);
    const auto *fsi1_63 = buffer.data(fsi1 + 63);
    const auto *fsi1_64 = buffer.data(fsi1 + 64);
    const auto *fsi1_65 = buffer.data(fsi1 + 65);
    const auto *fsi1_67 = buffer.data(fsi1 + 67);
    const auto *fsi1_68 = buffer.data(fsi1 + 68);
    const auto *fsi1_69 = buffer.data(fsi1 + 69);
    const auto *fsi1_70 = buffer.data(fsi1 + 70);
    const auto *fsi1_78 = buffer.data(fsi1 + 78);
    const auto *fsi1_79 = buffer.data(fsi1 + 79);

    const auto *fsk_0 = buffer.data(fsk + 0);
    const auto *fsk_1 = buffer.data(fsk + 1);
    const auto *fsk_2 = buffer.data(fsk + 2);
    const auto *fsk_3 = buffer.data(fsk + 3);
    const auto *fsk_5 = buffer.data(fsk + 5);
    const auto *fsk_6 = buffer.data(fsk + 6);
    const auto *fsk_8 = buffer.data(fsk + 8);
    const auto *fsk_9 = buffer.data(fsk + 9);
    const auto *fsk_10 = buffer.data(fsk + 10);
    const auto *fsk_12 = buffer.data(fsk + 12);
    const auto *fsk_13 = buffer.data(fsk + 13);
    const auto *fsk_14 = buffer.data(fsk + 14);
    const auto *fsk_15 = buffer.data(fsk + 15);
    const auto *fsk_17 = buffer.data(fsk + 17);
    const auto *fsk_18 = buffer.data(fsk + 18);
    const auto *fsk_19 = buffer.data(fsk + 19);
    const auto *fsk_20 = buffer.data(fsk + 20);
    const auto *fsk_21 = buffer.data(fsk + 21);
    const auto *fsk_27 = buffer.data(fsk + 27);
    const auto *fsk_28 = buffer.data(fsk + 28);
    const auto *fsk_30 = buffer.data(fsk + 30);
    const auto *fsk_31 = buffer.data(fsk + 31);
    const auto *fsk_32 = buffer.data(fsk + 32);
    const auto *fsk_33 = buffer.data(fsk + 33);
    const auto *fsk_34 = buffer.data(fsk + 34);
    const auto *fsk_35 = buffer.data(fsk + 35);
    const auto *fsk_36 = buffer.data(fsk + 36);
    const auto *fsk_37 = buffer.data(fsk + 37);
    const auto *fsk_39 = buffer.data(fsk + 39);
    const auto *fsk_41 = buffer.data(fsk + 41);
    const auto *fsk_42 = buffer.data(fsk + 42);
    const auto *fsk_43 = buffer.data(fsk + 43);
    const auto *fsk_45 = buffer.data(fsk + 45);
    const auto *fsk_46 = buffer.data(fsk + 46);
    const auto *fsk_47 = buffer.data(fsk + 47);
    const auto *fsk_48 = buffer.data(fsk + 48);
    const auto *fsk_50 = buffer.data(fsk + 50);
    const auto *fsk_51 = buffer.data(fsk + 51);
    const auto *fsk_52 = buffer.data(fsk + 52);
    const auto *fsk_53 = buffer.data(fsk + 53);
    const auto *fsk_54 = buffer.data(fsk + 54);
    const auto *fsk_56 = buffer.data(fsk + 56);
    const auto *fsk_57 = buffer.data(fsk + 57);
    const auto *fsk_64 = buffer.data(fsk + 64);
    const auto *fsk_65 = buffer.data(fsk + 65);
    const auto *fsk_66 = buffer.data(fsk + 66);
    const auto *fsk_67 = buffer.data(fsk + 67);
    const auto *fsk_68 = buffer.data(fsk + 68);
    const auto *fsk_69 = buffer.data(fsk + 69);
    const auto *fsk_70 = buffer.data(fsk + 70);
    const auto *fsk_71 = buffer.data(fsk + 71);
    const auto *fsk_72 = buffer.data(fsk + 72);
    const auto *fsk_74 = buffer.data(fsk + 74);
    const auto *fsk_76 = buffer.data(fsk + 76);
    const auto *fsk_77 = buffer.data(fsk + 77);
    const auto *fsk_79 = buffer.data(fsk + 79);
    const auto *fsk_80 = buffer.data(fsk + 80);
    const auto *fsk_81 = buffer.data(fsk + 81);
    const auto *fsk_83 = buffer.data(fsk + 83);
    const auto *fsk_84 = buffer.data(fsk + 84);
    const auto *fsk_85 = buffer.data(fsk + 85);
    const auto *fsk_86 = buffer.data(fsk + 86);
    const auto *fsk_88 = buffer.data(fsk + 88);
    const auto *fsk_89 = buffer.data(fsk + 89);
    const auto *fsk_90 = buffer.data(fsk + 90);
    const auto *fsk_91 = buffer.data(fsk + 91);
    const auto *fsk_92 = buffer.data(fsk + 92);
    const auto *fsk_99 = buffer.data(fsk + 99);
    const auto *fsk_100 = buffer.data(fsk + 100);
    const auto *fsk_101 = buffer.data(fsk + 101);
    const auto *fsk_102 = buffer.data(fsk + 102);
    const auto *fsk_103 = buffer.data(fsk + 103);
    const auto *fsk_104 = buffer.data(fsk + 104);
    const auto *fsk_105 = buffer.data(fsk + 105);
    const auto *fsk_107 = buffer.data(fsk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dsk_0, fsi0_0, \
                         fsi1_0, fsk_0, fsk_1, fsk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsk_0[k]
                 + f_1 * fsi0_0[k]
                 - f_2 * fsi1_0[k]
                 + f_3 * pc_x[k] * fsk_0[k];

        t_1[k] = f_3 * pc_y[k] * fsk_0[k];

        t_2[k] = f_3 * pc_z[k] * fsk_0[k];

        t_3[k] = f_4 * fsi0_0[k]
                 - f_5 * fsi1_0[k]
                 + f_3 * pc_y[k] * fsk_1[k];

        t_4[k] = f_3 * pc_y[k] * fsk_2[k];

        t_5[k] = f_4 * fsi0_0[k]
                 - f_5 * fsi1_0[k]
                 + f_3 * pc_z[k] * fsk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, fsi0_1, fsi0_2, fsi0_3, fsi1_1, \
                         fsi1_2, fsi1_3, fsk_3, fsk_5, fsk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * fsi0_1[k]
                 - f_7 * fsi1_1[k]
                 + f_3 * pc_y[k] * fsk_3[k];

        t_7[k] = f_3 * pc_z[k] * fsk_3[k];

        t_8[k] = f_3 * pc_y[k] * fsk_5[k];

        t_9[k] = f_6 * fsi0_2[k]
                 - f_7 * fsi1_2[k]
                 + f_3 * pc_z[k] * fsk_5[k];

        t_10[k] = f_8 * fsi0_3[k]
                  - f_9 * fsi1_3[k]
                  + f_3 * pc_y[k] * fsk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, fsi0_5, fsi0_6, \
                         fsi1_5, fsi1_6, fsk_6, fsk_8, fsk_9, fsk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * fsk_6[k];

        t_12[k] = f_4 * fsi0_5[k]
                  - f_5 * fsi1_5[k]
                  + f_3 * pc_y[k] * fsk_8[k];

        t_13[k] = f_3 * pc_y[k] * fsk_9[k];

        t_14[k] = f_8 * fsi0_5[k]
                  - f_9 * fsi1_5[k]
                  + f_3 * pc_z[k] * fsk_9[k];

        t_15[k] = f_10 * fsi0_6[k]
                  - f_11 * fsi1_6[k]
                  + f_3 * pc_y[k] * fsk_10[k];

        t_16[k] = f_3 * pc_z[k] * fsk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, fsi0_8, fsi0_9, fsi1_8, fsi1_9, \
                         fsk_12, fsk_13, fsk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fsi0_8[k]
                  - f_7 * fsi1_8[k]
                  + f_3 * pc_y[k] * fsk_12[k];

        t_18[k] = f_4 * fsi0_9[k]
                  - f_5 * fsi1_9[k]
                  + f_3 * pc_y[k] * fsk_13[k];

        t_19[k] = f_3 * pc_y[k] * fsk_14[k];

        t_20[k] = f_10 * fsi0_9[k]
                  - f_11 * fsi1_9[k]
                  + f_3 * pc_z[k] * fsk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, fsi0_10, fsi0_12, fsi0_13, \
                         fsi1_10, fsi1_12, fsi1_13, fsk_15, fsk_17, \
                         fsk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * fsi0_10[k]
                  - f_13 * fsi1_10[k]
                  + f_3 * pc_y[k] * fsk_15[k];

        t_22[k] = f_3 * pc_z[k] * fsk_15[k];

        t_23[k] = f_8 * fsi0_12[k]
                  - f_9 * fsi1_12[k]
                  + f_3 * pc_y[k] * fsk_17[k];

        t_24[k] = f_6 * fsi0_13[k]
                  - f_7 * fsi1_13[k]
                  + f_3 * pc_y[k] * fsk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, dsk_28, fsi0_14, \
                         fsi1_14, fsk_19, fsk_20, fsk_21, fsk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * fsi0_14[k]
                  - f_5 * fsi1_14[k]
                  + f_3 * pc_y[k] * fsk_19[k];

        t_26[k] = f_3 * pc_y[k] * fsk_20[k];

        t_27[k] = f_12 * fsi0_14[k]
                  - f_13 * fsi1_14[k]
                  + f_3 * pc_z[k] * fsk_20[k];

        t_28[k] = f_0 * dsk_28[k]
                  + f_3 * pc_x[k] * fsk_28[k];

        t_29[k] = f_3 * pc_z[k] * fsk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, dsk_30, dsk_31, dsk_32, \
                         dsk_33, fsk_27, fsk_30, fsk_31, fsk_32, \
                         fsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dsk_30[k]
                  + f_3 * pc_x[k] * fsk_30[k];

        t_31[k] = f_0 * dsk_31[k]
                  + f_3 * pc_x[k] * fsk_31[k];

        t_32[k] = f_0 * dsk_32[k]
                  + f_3 * pc_x[k] * fsk_32[k];

        t_33[k] = f_0 * dsk_33[k]
                  + f_3 * pc_x[k] * fsk_33[k];

        t_34[k] = f_3 * pc_y[k] * fsk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, dsk_35, fsi0_21, fsi0_23, \
                         fsi1_21, fsi1_23, fsk_28, fsk_30, fsk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * dsk_35[k]
                  + f_3 * pc_x[k] * fsk_35[k];

        t_36[k] = f_1 * fsi0_21[k]
                  - f_2 * fsi1_21[k]
                  + f_3 * pc_y[k] * fsk_28[k];

        t_37[k] = f_3 * pc_z[k] * fsk_28[k];

        t_38[k] = f_12 * fsi0_23[k]
                  - f_13 * fsi1_23[k]
                  + f_3 * pc_y[k] * fsk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, fsi0_24, fsi0_25, fsi0_26, fsi1_24, fsi1_25, \
                         fsi1_26, fsk_31, fsk_32, fsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * fsi0_24[k]
                  - f_11 * fsi1_24[k]
                  + f_3 * pc_y[k] * fsk_31[k];

        t_40[k] = f_8 * fsi0_25[k]
                  - f_9 * fsi1_25[k]
                  + f_3 * pc_y[k] * fsk_32[k];

        t_41[k] = f_6 * fsi0_26[k]
                  - f_7 * fsi1_26[k]
                  + f_3 * pc_y[k] * fsk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, dsl0_0, dsk_0, \
                         dsl1_0, fsi0_27, fsi1_27, fsk_34, fsk_35, \
                         fsk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * fsi0_27[k]
                  - f_5 * fsi1_27[k]
                  + f_3 * pc_y[k] * fsk_34[k];

        t_43[k] = f_3 * pc_y[k] * fsk_35[k];

        t_44[k] = f_1 * fsi0_27[k]
                  - f_2 * fsi1_27[k]
                  + f_3 * pc_z[k] * fsk_35[k];

        t_45[k] = pa_y[k] * dsl0_0[k]
                  - f_14 * pc_y[k] * dsl1_0[k];

        t_46[k] = f_15 * dsk_0[k]
                  + f_3 * pc_y[k] * fsk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, dsl0_3, dsl0_5, dsk_1, \
                         dsl1_3, dsl1_5, fsk_36, fsk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * fsk_36[k];

        t_48[k] = pa_y[k] * dsl0_3[k]
                  + f_16 * dsk_1[k]
                  - f_14 * pc_y[k] * dsl1_3[k];

        t_49[k] = f_3 * pc_z[k] * fsk_37[k];

        t_50[k] = pa_y[k] * dsl0_5[k]
                  - f_14 * pc_y[k] * dsl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, dsl0_6, dsl0_9, dsk_3, \
                         dsk_5, dsl1_6, dsl1_9, fsk_39, fsk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * dsl0_6[k]
                  + f_0 * dsk_3[k]
                  - f_14 * pc_y[k] * dsl1_6[k];

        t_52[k] = f_3 * pc_z[k] * fsk_39[k];

        t_53[k] = f_15 * dsk_5[k]
                  + f_3 * pc_y[k] * fsk_41[k];

        t_54[k] = pa_y[k] * dsl0_9[k]
                  - f_14 * pc_y[k] * dsl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, dsl0_10, dsk_6, dsk_9, \
                         dsl1_10, fsi0_31, fsi1_31, fsk_42, fsk_43, \
                         fsk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * dsl0_10[k]
                  + f_17 * dsk_6[k]
                  - f_14 * pc_y[k] * dsl1_10[k];

        t_56[k] = f_3 * pc_z[k] * fsk_42[k];

        t_57[k] = f_4 * fsi0_31[k]
                  - f_5 * fsi1_31[k]
                  + f_3 * pc_z[k] * fsk_43[k];

        t_58[k] = f_15 * dsk_9[k]
                  + f_3 * pc_y[k] * fsk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, dsl0_14, dsl0_15, dsk_10, \
                         dsl1_14, dsl1_15, fsi0_34, fsi1_34, fsk_46, \
                         fsk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * dsl0_14[k]
                  - f_14 * pc_y[k] * dsl1_14[k];

        t_60[k] = pa_y[k] * dsl0_15[k]
                  + f_18 * dsk_10[k]
                  - f_14 * pc_y[k] * dsl1_15[k];

        t_61[k] = f_3 * pc_z[k] * fsk_46[k];

        t_62[k] = f_4 * fsi0_34[k]
                  - f_5 * fsi1_34[k]
                  + f_3 * pc_z[k] * fsk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, dsl0_20, dsk_14, dsl1_20, \
                         fsi0_35, fsi1_35, fsk_48, fsk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * fsi0_35[k]
                  - f_7 * fsi1_35[k]
                  + f_3 * pc_z[k] * fsk_48[k];

        t_64[k] = f_15 * dsk_14[k]
                  + f_3 * pc_y[k] * fsk_50[k];

        t_65[k] = pa_y[k] * dsl0_20[k]
                  - f_14 * pc_y[k] * dsl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, dsl0_21, dsk_15, dsl1_21, \
                         fsi0_38, fsi1_38, fsk_51, fsk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * dsl0_21[k]
                  + f_19 * dsk_15[k]
                  - f_14 * pc_y[k] * dsl1_21[k];

        t_67[k] = f_3 * pc_z[k] * fsk_51[k];

        t_68[k] = f_4 * fsi0_38[k]
                  - f_5 * fsi1_38[k]
                  + f_3 * pc_z[k] * fsk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, dsk_20, fsi0_39, fsi0_40, fsi1_39, \
                         fsi1_40, fsk_53, fsk_54, fsk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * fsi0_39[k]
                  - f_7 * fsi1_39[k]
                  + f_3 * pc_z[k] * fsk_53[k];

        t_70[k] = f_8 * fsi0_40[k]
                  - f_9 * fsi1_40[k]
                  + f_3 * pc_z[k] * fsk_54[k];

        t_71[k] = f_15 * dsk_20[k]
                  + f_3 * pc_y[k] * fsk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, dsl0_27, dsk_64, \
                         dsk_66, dsl1_27, fsk_57, fsk_64, fsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * dsl0_27[k]
                  - f_14 * pc_y[k] * dsl1_27[k];

        t_73[k] = f_16 * dsk_64[k]
                  + f_3 * pc_x[k] * fsk_64[k];

        t_74[k] = f_3 * pc_z[k] * fsk_57[k];

        t_75[k] = f_16 * dsk_66[k]
                  + f_3 * pc_x[k] * fsk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, dsk_67, dsk_68, dsk_69, dsk_70, \
                         dsk_71, fsk_67, fsk_68, fsk_69, fsk_70, \
                         fsk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_16 * dsk_67[k]
                  + f_3 * pc_x[k] * fsk_67[k];

        t_77[k] = f_16 * dsk_68[k]
                  + f_3 * pc_x[k] * fsk_68[k];

        t_78[k] = f_16 * dsk_69[k]
                  + f_3 * pc_x[k] * fsk_69[k];

        t_79[k] = f_16 * dsk_70[k]
                  + f_3 * pc_x[k] * fsk_70[k];

        t_80[k] = f_16 * dsk_71[k]
                  + f_3 * pc_x[k] * fsk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, dsk_28, fsi0_49, fsi0_50, \
                         fsi1_49, fsi1_50, fsk_64, fsk_65, fsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * dsk_28[k]
                  + f_1 * fsi0_49[k]
                  - f_2 * fsi1_49[k]
                  + f_3 * pc_y[k] * fsk_64[k];

        t_82[k] = f_3 * pc_z[k] * fsk_64[k];

        t_83[k] = f_4 * fsi0_49[k]
                  - f_5 * fsi1_49[k]
                  + f_3 * pc_z[k] * fsk_65[k];

        t_84[k] = f_6 * fsi0_50[k]
                  - f_7 * fsi1_50[k]
                  + f_3 * pc_z[k] * fsk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, fsi0_51, fsi0_52, fsi0_53, fsi1_51, fsi1_52, \
                         fsi1_53, fsk_67, fsk_68, fsk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * fsi0_51[k]
                  - f_9 * fsi1_51[k]
                  + f_3 * pc_z[k] * fsk_67[k];

        t_86[k] = f_10 * fsi0_52[k]
                  - f_11 * fsi1_52[k]
                  + f_3 * pc_z[k] * fsk_68[k];

        t_87[k] = f_12 * fsi0_53[k]
                  - f_13 * fsi1_53[k]
                  + f_3 * pc_z[k] * fsk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, dsl0_0, dsl0_44, \
                         dsk_35, dsl1_0, dsl1_44, fsk_71, fsk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * dsk_35[k]
                  + f_3 * pc_y[k] * fsk_71[k];

        t_89[k] = pa_y[k] * dsl0_44[k]
                  - f_14 * pc_y[k] * dsl1_44[k];

        t_90[k] = pa_z[k] * dsl0_0[k]
                  - f_14 * pc_z[k] * dsl1_0[k];

        t_91[k] = f_3 * pc_y[k] * fsk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, dsl0_3, dsl0_5, dsk_0, \
                         dsk_2, dsl1_3, dsl1_5, fsk_72, fsk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * dsk_0[k]
                  + f_3 * pc_z[k] * fsk_72[k];

        t_93[k] = pa_z[k] * dsl0_3[k]
                  - f_14 * pc_z[k] * dsl1_3[k];

        t_94[k] = f_3 * pc_y[k] * fsk_74[k];

        t_95[k] = pa_z[k] * dsl0_5[k]
                  + f_16 * dsk_2[k]
                  - f_14 * pc_z[k] * dsl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, dsl0_6, dsl0_9, dsk_5, \
                         dsl1_6, dsl1_9, fsi0_58, fsi1_58, fsk_76, \
                         fsk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * dsl0_6[k]
                  - f_14 * pc_z[k] * dsl1_6[k];

        t_97[k] = f_4 * fsi0_58[k]
                  - f_5 * fsi1_58[k]
                  + f_3 * pc_y[k] * fsk_76[k];

        t_98[k] = f_3 * pc_y[k] * fsk_77[k];

        t_99[k] = pa_z[k] * dsl0_9[k]
                  + f_0 * dsk_5[k]
                  - f_14 * pc_z[k] * dsl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, dsl0_10, dsl1_10, \
                         fsi0_60, fsi0_61, fsi1_60, fsi1_61, fsk_79, fsk_80, \
                         fsk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * dsl0_10[k]
                   - f_14 * pc_z[k] * dsl1_10[k];

        t_101[k] = f_6 * fsi0_60[k]
                   - f_7 * fsi1_60[k]
                   + f_3 * pc_y[k] * fsk_79[k];

        t_102[k] = f_4 * fsi0_61[k]
                   - f_5 * fsi1_61[k]
                   + f_3 * pc_y[k] * fsk_80[k];

        t_103[k] = f_3 * pc_y[k] * fsk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, dsl0_14, dsl0_15, dsk_9, \
                         dsl1_14, dsl1_15, fsi0_63, fsi1_63, fsk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * dsl0_14[k]
                   + f_17 * dsk_9[k]
                   - f_14 * pc_z[k] * dsl1_14[k];

        t_105[k] = pa_z[k] * dsl0_15[k]
                   - f_14 * pc_z[k] * dsl1_15[k];

        t_106[k] = f_8 * fsi0_63[k]
                   - f_9 * fsi1_63[k]
                   + f_3 * pc_y[k] * fsk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, fsi0_64, fsi0_65, fsi1_64, fsi1_65, \
                         fsk_84, fsk_85, fsk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * fsi0_64[k]
                   - f_7 * fsi1_64[k]
                   + f_3 * pc_y[k] * fsk_84[k];

        t_108[k] = f_4 * fsi0_65[k]
                   - f_5 * fsi1_65[k]
                   + f_3 * pc_y[k] * fsk_85[k];

        t_109[k] = f_3 * pc_y[k] * fsk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, dsl0_20, dsl0_21, dsk_14, \
                         dsl1_20, dsl1_21, fsi0_67, fsi1_67, fsk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * dsl0_20[k]
                   + f_18 * dsk_14[k]
                   - f_14 * pc_z[k] * dsl1_20[k];

        t_111[k] = pa_z[k] * dsl0_21[k]
                   - f_14 * pc_z[k] * dsl1_21[k];

        t_112[k] = f_10 * fsi0_67[k]
                   - f_11 * fsi1_67[k]
                   + f_3 * pc_y[k] * fsk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, fsi0_68, fsi0_69, fsi0_70, fsi1_68, \
                         fsi1_69, fsi1_70, fsk_89, fsk_90, fsk_91, \
                         fsk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * fsi0_68[k]
                   - f_9 * fsi1_68[k]
                   + f_3 * pc_y[k] * fsk_89[k];

        t_114[k] = f_6 * fsi0_69[k]
                   - f_7 * fsi1_69[k]
                   + f_3 * pc_y[k] * fsk_90[k];

        t_115[k] = f_4 * fsi0_70[k]
                   - f_5 * fsi1_70[k]
                   + f_3 * pc_y[k] * fsk_91[k];

        t_116[k] = f_3 * pc_y[k] * fsk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, dsl0_27, dsk_20, \
                         dsk_100, dsk_101, dsk_102, dsl1_27, fsk_100, fsk_101, \
                         fsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * dsl0_27[k]
                   + f_19 * dsk_20[k]
                   - f_14 * pc_z[k] * dsl1_27[k];

        t_118[k] = f_16 * dsk_100[k]
                   + f_3 * pc_x[k] * fsk_100[k];

        t_119[k] = f_16 * dsk_101[k]
                   + f_3 * pc_x[k] * fsk_101[k];

        t_120[k] = f_16 * dsk_102[k]
                   + f_3 * pc_x[k] * fsk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, dsk_103, dsk_104, \
                         dsk_105, dsk_107, fsk_99, fsk_103, fsk_104, fsk_105, \
                         fsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_16 * dsk_103[k]
                   + f_3 * pc_x[k] * fsk_103[k];

        t_122[k] = f_16 * dsk_104[k]
                   + f_3 * pc_x[k] * fsk_104[k];

        t_123[k] = f_16 * dsk_105[k]
                   + f_3 * pc_x[k] * fsk_105[k];

        t_124[k] = f_3 * pc_y[k] * fsk_99[k];

        t_125[k] = f_16 * dsk_107[k]
                   + f_3 * pc_x[k] * fsk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, dsl0_36, dsl1_36, fsi0_78, \
                         fsi0_79, fsi1_78, fsi1_79, fsk_101, fsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * dsl0_36[k]
                   - f_14 * pc_z[k] * dsl1_36[k];

        t_127[k] = f_20 * fsi0_78[k]
                   - f_21 * fsi1_78[k]
                   + f_3 * pc_y[k] * fsk_101[k];

        t_128[k] = f_12 * fsi0_79[k]
                   - f_13 * fsi1_79[k]
                   + f_3 * pc_y[k] * fsk_102[k];
    }
}

static auto
compute_prim_fsl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsl0,
                                                          const size_t dsk, const size_t dsl1,
                                                          const size_t fsi0, const size_t fsi1,
                                                          const size_t fsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsl0_48 = buffer.data(dsl0 + 48);
    const auto *dsl0_51 = buffer.data(dsl0 + 51);
    const auto *dsl0_55 = buffer.data(dsl0 + 55);
    const auto *dsl0_60 = buffer.data(dsl0 + 60);
    const auto *dsl0_66 = buffer.data(dsl0 + 66);
    const auto *dsl0_90 = buffer.data(dsl0 + 90);
    const auto *dsl0_95 = buffer.data(dsl0 + 95);
    const auto *dsl0_99 = buffer.data(dsl0 + 99);
    const auto *dsl0_104 = buffer.data(dsl0 + 104);
    const auto *dsl0_110 = buffer.data(dsl0 + 110);
    const auto *dsl0_117 = buffer.data(dsl0 + 117);
    const auto *dsl0_135 = buffer.data(dsl0 + 135);
    const auto *dsl0_138 = buffer.data(dsl0 + 138);
    const auto *dsl0_141 = buffer.data(dsl0 + 141);
    const auto *dsl0_145 = buffer.data(dsl0 + 145);
    const auto *dsl0_150 = buffer.data(dsl0 + 150);
    const auto *dsl0_156 = buffer.data(dsl0 + 156);
    const auto *dsl0_171 = buffer.data(dsl0 + 171);
    const auto *dsl0_173 = buffer.data(dsl0 + 173);
    const auto *dsl0_174 = buffer.data(dsl0 + 174);
    const auto *dsl0_175 = buffer.data(dsl0 + 175);
    const auto *dsl0_176 = buffer.data(dsl0 + 176);
    const auto *dsl0_177 = buffer.data(dsl0 + 177);
    const auto *dsl0_179 = buffer.data(dsl0 + 179);
    const auto *dsl0_192 = buffer.data(dsl0 + 192);
    const auto *dsl0_197 = buffer.data(dsl0 + 197);
    const auto *dsl0_198 = buffer.data(dsl0 + 198);
    const auto *dsl0_203 = buffer.data(dsl0 + 203);
    const auto *dsl0_204 = buffer.data(dsl0 + 204);
    const auto *dsl0_205 = buffer.data(dsl0 + 205);
    const auto *dsl0_216 = buffer.data(dsl0 + 216);
    const auto *dsl0_218 = buffer.data(dsl0 + 218);
    const auto *dsl0_219 = buffer.data(dsl0 + 219);
    const auto *dsl0_220 = buffer.data(dsl0 + 220);
    const auto *dsl0_221 = buffer.data(dsl0 + 221);
    const auto *dsl0_222 = buffer.data(dsl0 + 222);
    const auto *dsl0_224 = buffer.data(dsl0 + 224);
    const auto *dsl0_225 = buffer.data(dsl0 + 225);
    const auto *dsl0_230 = buffer.data(dsl0 + 230);
    const auto *dsl0_234 = buffer.data(dsl0 + 234);
    const auto *dsl0_239 = buffer.data(dsl0 + 239);
    const auto *dsl0_245 = buffer.data(dsl0 + 245);

    const auto *dsk_35 = buffer.data(dsk + 35);
    const auto *dsk_36 = buffer.data(dsk + 36);
    const auto *dsk_39 = buffer.data(dsk + 39);
    const auto *dsk_41 = buffer.data(dsk + 41);
    const auto *dsk_42 = buffer.data(dsk + 42);
    const auto *dsk_45 = buffer.data(dsk + 45);
    const auto *dsk_46 = buffer.data(dsk + 46);
    const auto *dsk_50 = buffer.data(dsk + 50);
    const auto *dsk_51 = buffer.data(dsk + 51);
    const auto *dsk_56 = buffer.data(dsk + 56);
    const auto *dsk_64 = buffer.data(dsk + 64);
    const auto *dsk_71 = buffer.data(dsk + 71);
    const auto *dsk_72 = buffer.data(dsk + 72);
    const auto *dsk_74 = buffer.data(dsk + 74);
    const auto *dsk_77 = buffer.data(dsk + 77);
    const auto *dsk_81 = buffer.data(dsk + 81);
    const auto *dsk_86 = buffer.data(dsk + 86);
    const auto *dsk_92 = buffer.data(dsk + 92);
    const auto *dsk_107 = buffer.data(dsk + 107);
    const auto *dsk_108 = buffer.data(dsk + 108);
    const auto *dsk_111 = buffer.data(dsk + 111);
    const auto *dsk_114 = buffer.data(dsk + 114);
    const auto *dsk_118 = buffer.data(dsk + 118);
    const auto *dsk_123 = buffer.data(dsk + 123);
    const auto *dsk_129 = buffer.data(dsk + 129);
    const auto *dsk_136 = buffer.data(dsk + 136);
    const auto *dsk_138 = buffer.data(dsk + 138);
    const auto *dsk_139 = buffer.data(dsk + 139);
    const auto *dsk_140 = buffer.data(dsk + 140);
    const auto *dsk_141 = buffer.data(dsk + 141);
    const auto *dsk_142 = buffer.data(dsk + 142);
    const auto *dsk_143 = buffer.data(dsk + 143);
    const auto *dsk_156 = buffer.data(dsk + 156);
    const auto *dsk_161 = buffer.data(dsk + 161);
    const auto *dsk_162 = buffer.data(dsk + 162);
    const auto *dsk_167 = buffer.data(dsk + 167);
    const auto *dsk_168 = buffer.data(dsk + 168);
    const auto *dsk_169 = buffer.data(dsk + 169);
    const auto *dsk_172 = buffer.data(dsk + 172);
    const auto *dsk_173 = buffer.data(dsk + 173);
    const auto *dsk_174 = buffer.data(dsk + 174);
    const auto *dsk_175 = buffer.data(dsk + 175);
    const auto *dsk_176 = buffer.data(dsk + 176);
    const auto *dsk_177 = buffer.data(dsk + 177);
    const auto *dsk_178 = buffer.data(dsk + 178);
    const auto *dsk_179 = buffer.data(dsk + 179);
    const auto *dsk_180 = buffer.data(dsk + 180);
    const auto *dsk_185 = buffer.data(dsk + 185);
    const auto *dsk_189 = buffer.data(dsk + 189);
    const auto *dsk_194 = buffer.data(dsk + 194);
    const auto *dsk_200 = buffer.data(dsk + 200);

    const auto *dsl1_48 = buffer.data(dsl1 + 48);
    const auto *dsl1_51 = buffer.data(dsl1 + 51);
    const auto *dsl1_55 = buffer.data(dsl1 + 55);
    const auto *dsl1_60 = buffer.data(dsl1 + 60);
    const auto *dsl1_66 = buffer.data(dsl1 + 66);
    const auto *dsl1_90 = buffer.data(dsl1 + 90);
    const auto *dsl1_95 = buffer.data(dsl1 + 95);
    const auto *dsl1_99 = buffer.data(dsl1 + 99);
    const auto *dsl1_104 = buffer.data(dsl1 + 104);
    const auto *dsl1_110 = buffer.data(dsl1 + 110);
    const auto *dsl1_117 = buffer.data(dsl1 + 117);
    const auto *dsl1_135 = buffer.data(dsl1 + 135);
    const auto *dsl1_138 = buffer.data(dsl1 + 138);
    const auto *dsl1_141 = buffer.data(dsl1 + 141);
    const auto *dsl1_145 = buffer.data(dsl1 + 145);
    const auto *dsl1_150 = buffer.data(dsl1 + 150);
    const auto *dsl1_156 = buffer.data(dsl1 + 156);
    const auto *dsl1_171 = buffer.data(dsl1 + 171);
    const auto *dsl1_173 = buffer.data(dsl1 + 173);
    const auto *dsl1_174 = buffer.data(dsl1 + 174);
    const auto *dsl1_175 = buffer.data(dsl1 + 175);
    const auto *dsl1_176 = buffer.data(dsl1 + 176);
    const auto *dsl1_177 = buffer.data(dsl1 + 177);
    const auto *dsl1_179 = buffer.data(dsl1 + 179);
    const auto *dsl1_192 = buffer.data(dsl1 + 192);
    const auto *dsl1_197 = buffer.data(dsl1 + 197);
    const auto *dsl1_198 = buffer.data(dsl1 + 198);
    const auto *dsl1_203 = buffer.data(dsl1 + 203);
    const auto *dsl1_204 = buffer.data(dsl1 + 204);
    const auto *dsl1_205 = buffer.data(dsl1 + 205);
    const auto *dsl1_216 = buffer.data(dsl1 + 216);
    const auto *dsl1_218 = buffer.data(dsl1 + 218);
    const auto *dsl1_219 = buffer.data(dsl1 + 219);
    const auto *dsl1_220 = buffer.data(dsl1 + 220);
    const auto *dsl1_221 = buffer.data(dsl1 + 221);
    const auto *dsl1_222 = buffer.data(dsl1 + 222);
    const auto *dsl1_224 = buffer.data(dsl1 + 224);
    const auto *dsl1_225 = buffer.data(dsl1 + 225);
    const auto *dsl1_230 = buffer.data(dsl1 + 230);
    const auto *dsl1_234 = buffer.data(dsl1 + 234);
    const auto *dsl1_239 = buffer.data(dsl1 + 239);
    const auto *dsl1_245 = buffer.data(dsl1 + 245);

    const auto *fsi0_80 = buffer.data(fsi0 + 80);
    const auto *fsi0_81 = buffer.data(fsi0 + 81);
    const auto *fsi0_82 = buffer.data(fsi0 + 82);
    const auto *fsi0_83 = buffer.data(fsi0 + 83);
    const auto *fsi0_84 = buffer.data(fsi0 + 84);
    const auto *fsi0_86 = buffer.data(fsi0 + 86);
    const auto *fsi0_87 = buffer.data(fsi0 + 87);
    const auto *fsi0_89 = buffer.data(fsi0 + 89);
    const auto *fsi0_90 = buffer.data(fsi0 + 90);
    const auto *fsi0_91 = buffer.data(fsi0 + 91);
    const auto *fsi0_93 = buffer.data(fsi0 + 93);
    const auto *fsi0_94 = buffer.data(fsi0 + 94);
    const auto *fsi0_95 = buffer.data(fsi0 + 95);
    const auto *fsi0_96 = buffer.data(fsi0 + 96);
    const auto *fsi0_98 = buffer.data(fsi0 + 98);
    const auto *fsi0_140 = buffer.data(fsi0 + 140);
    const auto *fsi0_141 = buffer.data(fsi0 + 141);
    const auto *fsi0_142 = buffer.data(fsi0 + 142);
    const auto *fsi0_143 = buffer.data(fsi0 + 143);
    const auto *fsi0_144 = buffer.data(fsi0 + 144);
    const auto *fsi0_145 = buffer.data(fsi0 + 145);
    const auto *fsi0_146 = buffer.data(fsi0 + 146);
    const auto *fsi0_147 = buffer.data(fsi0 + 147);
    const auto *fsi0_148 = buffer.data(fsi0 + 148);
    const auto *fsi0_149 = buffer.data(fsi0 + 149);
    const auto *fsi0_150 = buffer.data(fsi0 + 150);
    const auto *fsi0_151 = buffer.data(fsi0 + 151);
    const auto *fsi0_152 = buffer.data(fsi0 + 152);
    const auto *fsi0_153 = buffer.data(fsi0 + 153);
    const auto *fsi0_154 = buffer.data(fsi0 + 154);

    const auto *fsi1_80 = buffer.data(fsi1 + 80);
    const auto *fsi1_81 = buffer.data(fsi1 + 81);
    const auto *fsi1_82 = buffer.data(fsi1 + 82);
    const auto *fsi1_83 = buffer.data(fsi1 + 83);
    const auto *fsi1_84 = buffer.data(fsi1 + 84);
    const auto *fsi1_86 = buffer.data(fsi1 + 86);
    const auto *fsi1_87 = buffer.data(fsi1 + 87);
    const auto *fsi1_89 = buffer.data(fsi1 + 89);
    const auto *fsi1_90 = buffer.data(fsi1 + 90);
    const auto *fsi1_91 = buffer.data(fsi1 + 91);
    const auto *fsi1_93 = buffer.data(fsi1 + 93);
    const auto *fsi1_94 = buffer.data(fsi1 + 94);
    const auto *fsi1_95 = buffer.data(fsi1 + 95);
    const auto *fsi1_96 = buffer.data(fsi1 + 96);
    const auto *fsi1_98 = buffer.data(fsi1 + 98);
    const auto *fsi1_140 = buffer.data(fsi1 + 140);
    const auto *fsi1_141 = buffer.data(fsi1 + 141);
    const auto *fsi1_142 = buffer.data(fsi1 + 142);
    const auto *fsi1_143 = buffer.data(fsi1 + 143);
    const auto *fsi1_144 = buffer.data(fsi1 + 144);
    const auto *fsi1_145 = buffer.data(fsi1 + 145);
    const auto *fsi1_146 = buffer.data(fsi1 + 146);
    const auto *fsi1_147 = buffer.data(fsi1 + 147);
    const auto *fsi1_148 = buffer.data(fsi1 + 148);
    const auto *fsi1_149 = buffer.data(fsi1 + 149);
    const auto *fsi1_150 = buffer.data(fsi1 + 150);
    const auto *fsi1_151 = buffer.data(fsi1 + 151);
    const auto *fsi1_152 = buffer.data(fsi1 + 152);
    const auto *fsi1_153 = buffer.data(fsi1 + 153);
    const auto *fsi1_154 = buffer.data(fsi1 + 154);

    const auto *fsk_103 = buffer.data(fsk + 103);
    const auto *fsk_104 = buffer.data(fsk + 104);
    const auto *fsk_105 = buffer.data(fsk + 105);
    const auto *fsk_106 = buffer.data(fsk + 106);
    const auto *fsk_107 = buffer.data(fsk + 107);
    const auto *fsk_108 = buffer.data(fsk + 108);
    const auto *fsk_109 = buffer.data(fsk + 109);
    const auto *fsk_110 = buffer.data(fsk + 110);
    const auto *fsk_111 = buffer.data(fsk + 111);
    const auto *fsk_113 = buffer.data(fsk + 113);
    const auto *fsk_114 = buffer.data(fsk + 114);
    const auto *fsk_115 = buffer.data(fsk + 115);
    const auto *fsk_117 = buffer.data(fsk + 117);
    const auto *fsk_118 = buffer.data(fsk + 118);
    const auto *fsk_119 = buffer.data(fsk + 119);
    const auto *fsk_120 = buffer.data(fsk + 120);
    const auto *fsk_122 = buffer.data(fsk + 122);
    const auto *fsk_123 = buffer.data(fsk + 123);
    const auto *fsk_124 = buffer.data(fsk + 124);
    const auto *fsk_125 = buffer.data(fsk + 125);
    const auto *fsk_126 = buffer.data(fsk + 126);
    const auto *fsk_128 = buffer.data(fsk + 128);
    const auto *fsk_129 = buffer.data(fsk + 129);
    const auto *fsk_136 = buffer.data(fsk + 136);
    const auto *fsk_138 = buffer.data(fsk + 138);
    const auto *fsk_139 = buffer.data(fsk + 139);
    const auto *fsk_140 = buffer.data(fsk + 140);
    const auto *fsk_141 = buffer.data(fsk + 141);
    const auto *fsk_142 = buffer.data(fsk + 142);
    const auto *fsk_143 = buffer.data(fsk + 143);
    const auto *fsk_144 = buffer.data(fsk + 144);
    const auto *fsk_146 = buffer.data(fsk + 146);
    const auto *fsk_147 = buffer.data(fsk + 147);
    const auto *fsk_149 = buffer.data(fsk + 149);
    const auto *fsk_150 = buffer.data(fsk + 150);
    const auto *fsk_153 = buffer.data(fsk + 153);
    const auto *fsk_154 = buffer.data(fsk + 154);
    const auto *fsk_158 = buffer.data(fsk + 158);
    const auto *fsk_159 = buffer.data(fsk + 159);
    const auto *fsk_164 = buffer.data(fsk + 164);
    const auto *fsk_172 = buffer.data(fsk + 172);
    const auto *fsk_173 = buffer.data(fsk + 173);
    const auto *fsk_174 = buffer.data(fsk + 174);
    const auto *fsk_175 = buffer.data(fsk + 175);
    const auto *fsk_176 = buffer.data(fsk + 176);
    const auto *fsk_177 = buffer.data(fsk + 177);
    const auto *fsk_178 = buffer.data(fsk + 178);
    const auto *fsk_179 = buffer.data(fsk + 179);
    const auto *fsk_180 = buffer.data(fsk + 180);
    const auto *fsk_181 = buffer.data(fsk + 181);
    const auto *fsk_182 = buffer.data(fsk + 182);
    const auto *fsk_183 = buffer.data(fsk + 183);
    const auto *fsk_184 = buffer.data(fsk + 184);
    const auto *fsk_185 = buffer.data(fsk + 185);
    const auto *fsk_186 = buffer.data(fsk + 186);
    const auto *fsk_187 = buffer.data(fsk + 187);
    const auto *fsk_188 = buffer.data(fsk + 188);
    const auto *fsk_189 = buffer.data(fsk + 189);
    const auto *fsk_190 = buffer.data(fsk + 190);
    const auto *fsk_191 = buffer.data(fsk + 191);
    const auto *fsk_192 = buffer.data(fsk + 192);
    const auto *fsk_193 = buffer.data(fsk + 193);
    const auto *fsk_194 = buffer.data(fsk + 194);
    const auto *fsk_195 = buffer.data(fsk + 195);
    const auto *fsk_196 = buffer.data(fsk + 196);
    const auto *fsk_197 = buffer.data(fsk + 197);
    const auto *fsk_198 = buffer.data(fsk + 198);
    const auto *fsk_199 = buffer.data(fsk + 199);
    const auto *fsk_200 = buffer.data(fsk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, fsi0_80, fsi0_81, fsi0_82, fsi1_80, \
                         fsi1_81, fsi1_82, fsk_103, fsk_104, fsk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * fsi0_80[k]
                   - f_11 * fsi1_80[k]
                   + f_3 * pc_y[k] * fsk_103[k];

        t_130[k] = f_8 * fsi0_81[k]
                   - f_9 * fsi1_81[k]
                   + f_3 * pc_y[k] * fsk_104[k];

        t_131[k] = f_6 * fsi0_82[k]
                   - f_7 * fsi1_82[k]
                   + f_3 * pc_y[k] * fsk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pc_x, pc_y, pc_z, dsl0_135, dsk_35, \
                         dsk_108, dsl1_135, fsi0_83, fsi1_83, fsk_106, \
                         fsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * fsi0_83[k]
                   - f_5 * fsi1_83[k]
                   + f_3 * pc_y[k] * fsk_106[k];

        t_133[k] = f_3 * pc_y[k] * fsk_107[k];

        t_134[k] = f_15 * dsk_35[k]
                   + f_1 * fsi0_83[k]
                   - f_2 * fsi1_83[k]
                   + f_3 * pc_z[k] * fsk_107[k];

        t_135[k] = pa_x[k] * dsl0_135[k]
                   + f_22 * dsk_108[k]
                   - f_14 * pc_x[k] * dsl1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, pc_x, pc_y, pc_z, dsl0_138, dsk_36, \
                         dsk_111, dsl1_138, fsk_108, fsk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * dsk_36[k]
                   + f_3 * pc_y[k] * fsk_108[k];

        t_137[k] = f_3 * pc_z[k] * fsk_108[k];

        t_138[k] = pa_x[k] * dsl0_138[k]
                   + f_19 * dsk_111[k]
                   - f_14 * pc_x[k] * dsl1_138[k];

        t_139[k] = f_3 * pc_z[k] * fsk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_x, pc_x, pc_z, dsl0_141, dsk_114, dsl1_141, \
                         fsi0_84, fsi1_84, fsk_110, fsk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * fsi0_84[k]
                   - f_5 * fsi1_84[k]
                   + f_3 * pc_z[k] * fsk_110[k];

        t_141[k] = pa_x[k] * dsl0_141[k]
                   + f_18 * dsk_114[k]
                   - f_14 * pc_x[k] * dsl1_141[k];

        t_142[k] = f_3 * pc_z[k] * fsk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_x, pc_x, pc_y, pc_z, dsl0_145, dsk_41, \
                         dsk_118, dsl1_145, fsi0_86, fsi1_86, fsk_113, \
                         fsk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * dsk_41[k]
                   + f_3 * pc_y[k] * fsk_113[k];

        t_144[k] = f_6 * fsi0_86[k]
                   - f_7 * fsi1_86[k]
                   + f_3 * pc_z[k] * fsk_113[k];

        t_145[k] = pa_x[k] * dsl0_145[k]
                   + f_17 * dsk_118[k]
                   - f_14 * pc_x[k] * dsl1_145[k];

        t_146[k] = f_3 * pc_z[k] * fsk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, dsk_45, fsi0_87, fsi0_89, fsi1_87, \
                         fsi1_89, fsk_115, fsk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * fsi0_87[k]
                   - f_5 * fsi1_87[k]
                   + f_3 * pc_z[k] * fsk_115[k];

        t_148[k] = f_16 * dsk_45[k]
                   + f_3 * pc_y[k] * fsk_117[k];

        t_149[k] = f_8 * fsi0_89[k]
                   - f_9 * fsi1_89[k]
                   + f_3 * pc_z[k] * fsk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_x, pc_x, pc_z, dsl0_150, dsk_123, dsl1_150, \
                         fsi0_90, fsi1_90, fsk_118, fsk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_x[k] * dsl0_150[k]
                   + f_0 * dsk_123[k]
                   - f_14 * pc_x[k] * dsl1_150[k];

        t_151[k] = f_3 * pc_z[k] * fsk_118[k];

        t_152[k] = f_4 * fsi0_90[k]
                   - f_5 * fsi1_90[k]
                   + f_3 * pc_z[k] * fsk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, dsk_50, fsi0_91, fsi0_93, fsi1_91, \
                         fsi1_93, fsk_120, fsk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * fsi0_91[k]
                   - f_7 * fsi1_91[k]
                   + f_3 * pc_z[k] * fsk_120[k];

        t_154[k] = f_16 * dsk_50[k]
                   + f_3 * pc_y[k] * fsk_122[k];

        t_155[k] = f_10 * fsi0_93[k]
                   - f_11 * fsi1_93[k]
                   + f_3 * pc_z[k] * fsk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_x, pc_x, pc_z, dsl0_156, dsk_129, dsl1_156, \
                         fsi0_94, fsi1_94, fsk_123, fsk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_x[k] * dsl0_156[k]
                   + f_16 * dsk_129[k]
                   - f_14 * pc_x[k] * dsl1_156[k];

        t_157[k] = f_3 * pc_z[k] * fsk_123[k];

        t_158[k] = f_4 * fsi0_94[k]
                   - f_5 * fsi1_94[k]
                   + f_3 * pc_z[k] * fsk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, dsk_56, fsi0_95, fsi0_96, \
                         fsi0_98, fsi1_95, fsi1_96, fsi1_98, fsk_125, fsk_126, \
                         fsk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * fsi0_95[k]
                   - f_7 * fsi1_95[k]
                   + f_3 * pc_z[k] * fsk_125[k];

        t_160[k] = f_8 * fsi0_96[k]
                   - f_9 * fsi1_96[k]
                   + f_3 * pc_z[k] * fsk_126[k];

        t_161[k] = f_16 * dsk_56[k]
                   + f_3 * pc_y[k] * fsk_128[k];

        t_162[k] = f_12 * fsi0_98[k]
                   - f_13 * fsi1_98[k]
                   + f_3 * pc_z[k] * fsk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, dsk_136, dsk_138, \
                         dsk_139, dsk_140, fsk_129, fsk_136, fsk_138, fsk_139, \
                         fsk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_15 * dsk_136[k]
                   + f_3 * pc_x[k] * fsk_136[k];

        t_164[k] = f_3 * pc_z[k] * fsk_129[k];

        t_165[k] = f_15 * dsk_138[k]
                   + f_3 * pc_x[k] * fsk_138[k];

        t_166[k] = f_15 * dsk_139[k]
                   + f_3 * pc_x[k] * fsk_139[k];

        t_167[k] = f_15 * dsk_140[k]
                   + f_3 * pc_x[k] * fsk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pc_x, dsl0_171, dsk_141, dsk_142, \
                         dsk_143, dsl1_171, fsk_141, fsk_142, fsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_15 * dsk_141[k]
                   + f_3 * pc_x[k] * fsk_141[k];

        t_169[k] = f_15 * dsk_142[k]
                   + f_3 * pc_x[k] * fsk_142[k];

        t_170[k] = f_15 * dsk_143[k]
                   + f_3 * pc_x[k] * fsk_143[k];

        t_171[k] = pa_x[k] * dsl0_171[k]
                   - f_14 * pc_x[k] * dsl1_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_x, pc_x, pc_z, dsl0_173, dsl0_174, \
                         dsl0_175, dsl1_173, dsl1_174, dsl1_175, \
                         fsk_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * fsk_136[k];

        t_173[k] = pa_x[k] * dsl0_173[k]
                   - f_14 * pc_x[k] * dsl1_173[k];

        t_174[k] = pa_x[k] * dsl0_174[k]
                   - f_14 * pc_x[k] * dsl1_174[k];

        t_175[k] = pa_x[k] * dsl0_175[k]
                   - f_14 * pc_x[k] * dsl1_175[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pc_x, pc_y, dsl0_176, dsl0_177, \
                         dsl0_179, dsk_71, dsl1_176, dsl1_177, dsl1_179, \
                         fsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_x[k] * dsl0_176[k]
                   - f_14 * pc_x[k] * dsl1_176[k];

        t_177[k] = pa_x[k] * dsl0_177[k]
                   - f_14 * pc_x[k] * dsl1_177[k];

        t_178[k] = f_16 * dsk_71[k]
                   + f_3 * pc_y[k] * fsk_143[k];

        t_179[k] = pa_x[k] * dsl0_179[k]
                   - f_14 * pc_x[k] * dsl1_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, dsl0_48, dsl0_90, \
                         dsk_36, dsk_72, dsl1_48, dsl1_90, fsk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * dsl0_90[k]
                   - f_14 * pc_y[k] * dsl1_90[k];

        t_181[k] = f_15 * dsk_72[k]
                   + f_3 * pc_y[k] * fsk_144[k];

        t_182[k] = f_15 * dsk_36[k]
                   + f_3 * pc_z[k] * fsk_144[k];

        t_183[k] = pa_z[k] * dsl0_48[k]
                   - f_14 * pc_z[k] * dsl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, dsl0_51, dsl0_95, \
                         dsk_39, dsk_74, dsl1_51, dsl1_95, fsk_146, \
                         fsk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * dsk_74[k]
                   + f_3 * pc_y[k] * fsk_146[k];

        t_185[k] = pa_y[k] * dsl0_95[k]
                   - f_14 * pc_y[k] * dsl1_95[k];

        t_186[k] = pa_z[k] * dsl0_51[k]
                   - f_14 * pc_z[k] * dsl1_51[k];

        t_187[k] = f_15 * dsk_39[k]
                   + f_3 * pc_z[k] * fsk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, dsl0_55, dsl0_99, \
                         dsk_42, dsk_77, dsl1_55, dsl1_99, fsk_149, \
                         fsk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * dsk_77[k]
                   + f_3 * pc_y[k] * fsk_149[k];

        t_189[k] = pa_y[k] * dsl0_99[k]
                   - f_14 * pc_y[k] * dsl1_99[k];

        t_190[k] = pa_z[k] * dsl0_55[k]
                   - f_14 * pc_z[k] * dsl1_55[k];

        t_191[k] = f_15 * dsk_42[k]
                   + f_3 * pc_z[k] * fsk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_x, pa_y, pc_x, pc_y, dsl0_104, dsl0_192, \
                         dsk_81, dsk_156, dsl1_104, dsl1_192, fsk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_x[k] * dsl0_192[k]
                   + f_17 * dsk_156[k]
                   - f_14 * pc_x[k] * dsl1_192[k];

        t_193[k] = f_15 * dsk_81[k]
                   + f_3 * pc_y[k] * fsk_153[k];

        t_194[k] = pa_y[k] * dsl0_104[k]
                   - f_14 * pc_y[k] * dsl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_x, pa_z, pc_x, pc_z, dsl0_60, dsl0_197, \
                         dsk_46, dsk_161, dsl1_60, dsl1_197, fsk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * dsl0_60[k]
                   - f_14 * pc_z[k] * dsl1_60[k];

        t_196[k] = f_15 * dsk_46[k]
                   + f_3 * pc_z[k] * fsk_154[k];

        t_197[k] = pa_x[k] * dsl0_197[k]
                   + f_0 * dsk_161[k]
                   - f_14 * pc_x[k] * dsl1_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_x, pa_y, pc_x, pc_y, dsl0_110, dsl0_198, \
                         dsk_86, dsk_162, dsl1_110, dsl1_198, fsk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_x[k] * dsl0_198[k]
                   + f_0 * dsk_162[k]
                   - f_14 * pc_x[k] * dsl1_198[k];

        t_199[k] = f_15 * dsk_86[k]
                   + f_3 * pc_y[k] * fsk_158[k];

        t_200[k] = pa_y[k] * dsl0_110[k]
                   - f_14 * pc_y[k] * dsl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_x, pa_z, pc_x, pc_z, dsl0_66, dsl0_203, \
                         dsk_51, dsk_167, dsl1_66, dsl1_203, fsk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * dsl0_66[k]
                   - f_14 * pc_z[k] * dsl1_66[k];

        t_202[k] = f_15 * dsk_51[k]
                   + f_3 * pc_z[k] * fsk_159[k];

        t_203[k] = pa_x[k] * dsl0_203[k]
                   + f_16 * dsk_167[k]
                   - f_14 * pc_x[k] * dsl1_203[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_x, pc_x, pc_y, dsl0_204, dsl0_205, dsk_92, \
                         dsk_168, dsk_169, dsl1_204, dsl1_205, \
                         fsk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_x[k] * dsl0_204[k]
                   + f_16 * dsk_168[k]
                   - f_14 * pc_x[k] * dsl1_204[k];

        t_205[k] = pa_x[k] * dsl0_205[k]
                   + f_16 * dsk_169[k]
                   - f_14 * pc_x[k] * dsl1_205[k];

        t_206[k] = f_15 * dsk_92[k]
                   + f_3 * pc_y[k] * fsk_164[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pa_y, pc_x, pc_y, dsl0_117, dsk_172, \
                         dsk_173, dsk_174, dsl1_117, fsk_172, fsk_173, \
                         fsk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_y[k] * dsl0_117[k]
                   - f_14 * pc_y[k] * dsl1_117[k];

        t_208[k] = f_15 * dsk_172[k]
                   + f_3 * pc_x[k] * fsk_172[k];

        t_209[k] = f_15 * dsk_173[k]
                   + f_3 * pc_x[k] * fsk_173[k];

        t_210[k] = f_15 * dsk_174[k]
                   + f_3 * pc_x[k] * fsk_174[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pc_x, dsk_175, dsk_176, dsk_177, \
                         dsk_178, dsk_179, fsk_175, fsk_176, fsk_177, fsk_178, \
                         fsk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_15 * dsk_175[k]
                   + f_3 * pc_x[k] * fsk_175[k];

        t_212[k] = f_15 * dsk_176[k]
                   + f_3 * pc_x[k] * fsk_176[k];

        t_213[k] = f_15 * dsk_177[k]
                   + f_3 * pc_x[k] * fsk_177[k];

        t_214[k] = f_15 * dsk_178[k]
                   + f_3 * pc_x[k] * fsk_178[k];

        t_215[k] = f_15 * dsk_179[k]
                   + f_3 * pc_x[k] * fsk_179[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_x, pc_x, pc_z, dsl0_216, dsl0_218, \
                         dsl0_219, dsk_64, dsl1_216, dsl1_218, dsl1_219, \
                         fsk_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pa_x[k] * dsl0_216[k]
                   - f_14 * pc_x[k] * dsl1_216[k];

        t_217[k] = f_15 * dsk_64[k]
                   + f_3 * pc_z[k] * fsk_172[k];

        t_218[k] = pa_x[k] * dsl0_218[k]
                   - f_14 * pc_x[k] * dsl1_218[k];

        t_219[k] = pa_x[k] * dsl0_219[k]
                   - f_14 * pc_x[k] * dsl1_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_x, pc_x, pc_y, dsl0_220, dsl0_221, \
                         dsl0_222, dsk_107, dsl1_220, dsl1_221, dsl1_222, \
                         fsk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_x[k] * dsl0_220[k]
                   - f_14 * pc_x[k] * dsl1_220[k];

        t_221[k] = pa_x[k] * dsl0_221[k]
                   - f_14 * pc_x[k] * dsl1_221[k];

        t_222[k] = pa_x[k] * dsl0_222[k]
                   - f_14 * pc_x[k] * dsl1_222[k];

        t_223[k] = f_15 * dsk_107[k]
                   + f_3 * pc_y[k] * fsk_179[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pc_x, pc_y, pc_z, dsl0_224, \
                         dsl0_225, dsk_72, dsk_180, dsl1_224, dsl1_225, \
                         fsk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_x[k] * dsl0_224[k]
                   - f_14 * pc_x[k] * dsl1_224[k];

        t_225[k] = pa_x[k] * dsl0_225[k]
                   + f_22 * dsk_180[k]
                   - f_14 * pc_x[k] * dsl1_225[k];

        t_226[k] = f_3 * pc_y[k] * fsk_180[k];

        t_227[k] = f_16 * dsk_72[k]
                   + f_3 * pc_z[k] * fsk_180[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_x, pc_x, pc_y, dsl0_230, dsk_185, dsl1_230, \
                         fsi0_140, fsi1_140, fsk_181, fsk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * fsi0_140[k]
                   - f_5 * fsi1_140[k]
                   + f_3 * pc_y[k] * fsk_181[k];

        t_229[k] = f_3 * pc_y[k] * fsk_182[k];

        t_230[k] = pa_x[k] * dsl0_230[k]
                   + f_19 * dsk_185[k]
                   - f_14 * pc_x[k] * dsl1_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_y, fsi0_141, fsi0_142, fsi1_141, fsi1_142, \
                         fsk_183, fsk_184, fsk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_6 * fsi0_141[k]
                   - f_7 * fsi1_141[k]
                   + f_3 * pc_y[k] * fsk_183[k];

        t_232[k] = f_4 * fsi0_142[k]
                   - f_5 * fsi1_142[k]
                   + f_3 * pc_y[k] * fsk_184[k];

        t_233[k] = f_3 * pc_y[k] * fsk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_x, pc_x, pc_y, dsl0_234, dsk_189, dsl1_234, \
                         fsi0_143, fsi0_144, fsi1_143, fsi1_144, fsk_186, \
                         fsk_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_x[k] * dsl0_234[k]
                   + f_18 * dsk_189[k]
                   - f_14 * pc_x[k] * dsl1_234[k];

        t_235[k] = f_8 * fsi0_143[k]
                   - f_9 * fsi1_143[k]
                   + f_3 * pc_y[k] * fsk_186[k];

        t_236[k] = f_6 * fsi0_144[k]
                   - f_7 * fsi1_144[k]
                   + f_3 * pc_y[k] * fsk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_x, pc_x, pc_y, dsl0_239, dsk_194, dsl1_239, \
                         fsi0_145, fsi1_145, fsk_188, fsk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * fsi0_145[k]
                   - f_5 * fsi1_145[k]
                   + f_3 * pc_y[k] * fsk_188[k];

        t_238[k] = f_3 * pc_y[k] * fsk_189[k];

        t_239[k] = pa_x[k] * dsl0_239[k]
                   + f_17 * dsk_194[k]
                   - f_14 * pc_x[k] * dsl1_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, fsi0_146, fsi0_147, fsi0_148, fsi1_146, \
                         fsi1_147, fsi1_148, fsk_190, fsk_191, \
                         fsk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * fsi0_146[k]
                   - f_11 * fsi1_146[k]
                   + f_3 * pc_y[k] * fsk_190[k];

        t_241[k] = f_8 * fsi0_147[k]
                   - f_9 * fsi1_147[k]
                   + f_3 * pc_y[k] * fsk_191[k];

        t_242[k] = f_6 * fsi0_148[k]
                   - f_7 * fsi1_148[k]
                   + f_3 * pc_y[k] * fsk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_x, pc_x, pc_y, dsl0_245, dsk_200, dsl1_245, \
                         fsi0_149, fsi1_149, fsk_193, fsk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * fsi0_149[k]
                   - f_5 * fsi1_149[k]
                   + f_3 * pc_y[k] * fsk_193[k];

        t_244[k] = f_3 * pc_y[k] * fsk_194[k];

        t_245[k] = pa_x[k] * dsl0_245[k]
                   + f_0 * dsk_200[k]
                   - f_14 * pc_x[k] * dsl1_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, fsi0_150, fsi0_151, fsi0_152, fsi1_150, \
                         fsi1_151, fsi1_152, fsk_195, fsk_196, \
                         fsk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * fsi0_150[k]
                   - f_13 * fsi1_150[k]
                   + f_3 * pc_y[k] * fsk_195[k];

        t_247[k] = f_10 * fsi0_151[k]
                   - f_11 * fsi1_151[k]
                   + f_3 * pc_y[k] * fsk_196[k];

        t_248[k] = f_8 * fsi0_152[k]
                   - f_9 * fsi1_152[k]
                   + f_3 * pc_y[k] * fsk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, fsi0_153, fsi0_154, fsi1_153, fsi1_154, \
                         fsk_198, fsk_199, fsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * fsi0_153[k]
                   - f_7 * fsi1_153[k]
                   + f_3 * pc_y[k] * fsk_198[k];

        t_250[k] = f_4 * fsi0_154[k]
                   - f_5 * fsi1_154[k]
                   + f_3 * pc_y[k] * fsk_199[k];

        t_251[k] = f_3 * pc_y[k] * fsk_200[k];
    }
}

static auto
compute_prim_fsl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsl0,
                                                          const size_t dsk, const size_t dsl1,
                                                          const size_t fsi0, const size_t fsi1,
                                                          const size_t fsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsl0_135 = buffer.data(dsl0 + 135);
    const auto *dsl0_136 = buffer.data(dsl0 + 136);
    const auto *dsl0_138 = buffer.data(dsl0 + 138);
    const auto *dsl0_141 = buffer.data(dsl0 + 141);
    const auto *dsl0_145 = buffer.data(dsl0 + 145);
    const auto *dsl0_150 = buffer.data(dsl0 + 150);
    const auto *dsl0_156 = buffer.data(dsl0 + 156);
    const auto *dsl0_171 = buffer.data(dsl0 + 171);
    const auto *dsl0_173 = buffer.data(dsl0 + 173);
    const auto *dsl0_174 = buffer.data(dsl0 + 174);
    const auto *dsl0_175 = buffer.data(dsl0 + 175);
    const auto *dsl0_176 = buffer.data(dsl0 + 176);
    const auto *dsl0_177 = buffer.data(dsl0 + 177);
    const auto *dsl0_225 = buffer.data(dsl0 + 225);
    const auto *dsl0_227 = buffer.data(dsl0 + 227);
    const auto *dsl0_230 = buffer.data(dsl0 + 230);
    const auto *dsl0_234 = buffer.data(dsl0 + 234);
    const auto *dsl0_252 = buffer.data(dsl0 + 252);
    const auto *dsl0_261 = buffer.data(dsl0 + 261);
    const auto *dsl0_262 = buffer.data(dsl0 + 262);
    const auto *dsl0_263 = buffer.data(dsl0 + 263);
    const auto *dsl0_264 = buffer.data(dsl0 + 264);
    const auto *dsl0_265 = buffer.data(dsl0 + 265);
    const auto *dsl0_266 = buffer.data(dsl0 + 266);
    const auto *dsl0_267 = buffer.data(dsl0 + 267);
    const auto *dsl0_269 = buffer.data(dsl0 + 269);

    const auto *dsk_136 = buffer.data(dsk + 136);
    const auto *dsk_137 = buffer.data(dsk + 137);
    const auto *dsk_138 = buffer.data(dsk + 138);
    const auto *dsk_139 = buffer.data(dsk + 139);
    const auto *dsk_140 = buffer.data(dsk + 140);
    const auto *dsk_141 = buffer.data(dsk + 141);
    const auto *dsk_143 = buffer.data(dsk + 143);
    const auto *dsk_179 = buffer.data(dsk + 179);
    const auto *dsk_207 = buffer.data(dsk + 207);
    const auto *dsk_208 = buffer.data(dsk + 208);
    const auto *dsk_209 = buffer.data(dsk + 209);
    const auto *dsk_210 = buffer.data(dsk + 210);
    const auto *dsk_211 = buffer.data(dsk + 211);
    const auto *dsk_212 = buffer.data(dsk + 212);
    const auto *dsk_213 = buffer.data(dsk + 213);
    const auto *dsk_215 = buffer.data(dsk + 215);

    const auto *dsl1_135 = buffer.data(dsl1 + 135);
    const auto *dsl1_136 = buffer.data(dsl1 + 136);
    const auto *dsl1_138 = buffer.data(dsl1 + 138);
    const auto *dsl1_141 = buffer.data(dsl1 + 141);
    const auto *dsl1_145 = buffer.data(dsl1 + 145);
    const auto *dsl1_150 = buffer.data(dsl1 + 150);
    const auto *dsl1_156 = buffer.data(dsl1 + 156);
    const auto *dsl1_171 = buffer.data(dsl1 + 171);
    const auto *dsl1_173 = buffer.data(dsl1 + 173);
    const auto *dsl1_174 = buffer.data(dsl1 + 174);
    const auto *dsl1_175 = buffer.data(dsl1 + 175);
    const auto *dsl1_176 = buffer.data(dsl1 + 176);
    const auto *dsl1_177 = buffer.data(dsl1 + 177);
    const auto *dsl1_225 = buffer.data(dsl1 + 225);
    const auto *dsl1_227 = buffer.data(dsl1 + 227);
    const auto *dsl1_230 = buffer.data(dsl1 + 230);
    const auto *dsl1_234 = buffer.data(dsl1 + 234);
    const auto *dsl1_252 = buffer.data(dsl1 + 252);
    const auto *dsl1_261 = buffer.data(dsl1 + 261);
    const auto *dsl1_262 = buffer.data(dsl1 + 262);
    const auto *dsl1_263 = buffer.data(dsl1 + 263);
    const auto *dsl1_264 = buffer.data(dsl1 + 264);
    const auto *dsl1_265 = buffer.data(dsl1 + 265);
    const auto *dsl1_266 = buffer.data(dsl1 + 266);
    const auto *dsl1_267 = buffer.data(dsl1 + 267);
    const auto *dsl1_269 = buffer.data(dsl1 + 269);

    const auto *fsi0_168 = buffer.data(fsi0 + 168);
    const auto *fsi0_169 = buffer.data(fsi0 + 169);
    const auto *fsi0_171 = buffer.data(fsi0 + 171);
    const auto *fsi0_173 = buffer.data(fsi0 + 173);
    const auto *fsi0_174 = buffer.data(fsi0 + 174);
    const auto *fsi0_176 = buffer.data(fsi0 + 176);
    const auto *fsi0_177 = buffer.data(fsi0 + 177);
    const auto *fsi0_178 = buffer.data(fsi0 + 178);
    const auto *fsi0_180 = buffer.data(fsi0 + 180);
    const auto *fsi0_181 = buffer.data(fsi0 + 181);
    const auto *fsi0_182 = buffer.data(fsi0 + 182);
    const auto *fsi0_183 = buffer.data(fsi0 + 183);
    const auto *fsi0_185 = buffer.data(fsi0 + 185);
    const auto *fsi0_186 = buffer.data(fsi0 + 186);
    const auto *fsi0_187 = buffer.data(fsi0 + 187);
    const auto *fsi0_188 = buffer.data(fsi0 + 188);
    const auto *fsi0_189 = buffer.data(fsi0 + 189);
    const auto *fsi0_190 = buffer.data(fsi0 + 190);
    const auto *fsi0_191 = buffer.data(fsi0 + 191);
    const auto *fsi0_192 = buffer.data(fsi0 + 192);
    const auto *fsi0_193 = buffer.data(fsi0 + 193);
    const auto *fsi0_194 = buffer.data(fsi0 + 194);
    const auto *fsi0_195 = buffer.data(fsi0 + 195);
    const auto *fsi0_198 = buffer.data(fsi0 + 198);
    const auto *fsi0_200 = buffer.data(fsi0 + 200);
    const auto *fsi0_201 = buffer.data(fsi0 + 201);
    const auto *fsi0_203 = buffer.data(fsi0 + 203);
    const auto *fsi0_204 = buffer.data(fsi0 + 204);
    const auto *fsi0_205 = buffer.data(fsi0 + 205);
    const auto *fsi0_207 = buffer.data(fsi0 + 207);
    const auto *fsi0_208 = buffer.data(fsi0 + 208);
    const auto *fsi0_209 = buffer.data(fsi0 + 209);
    const auto *fsi0_210 = buffer.data(fsi0 + 210);
    const auto *fsi0_212 = buffer.data(fsi0 + 212);
    const auto *fsi0_213 = buffer.data(fsi0 + 213);
    const auto *fsi0_214 = buffer.data(fsi0 + 214);
    const auto *fsi0_215 = buffer.data(fsi0 + 215);
    const auto *fsi0_216 = buffer.data(fsi0 + 216);
    const auto *fsi0_218 = buffer.data(fsi0 + 218);
    const auto *fsi0_219 = buffer.data(fsi0 + 219);
    const auto *fsi0_220 = buffer.data(fsi0 + 220);
    const auto *fsi0_221 = buffer.data(fsi0 + 221);
    const auto *fsi0_222 = buffer.data(fsi0 + 222);
    const auto *fsi0_223 = buffer.data(fsi0 + 223);
    const auto *fsi0_225 = buffer.data(fsi0 + 225);
    const auto *fsi0_227 = buffer.data(fsi0 + 227);
    const auto *fsi0_228 = buffer.data(fsi0 + 228);
    const auto *fsi0_230 = buffer.data(fsi0 + 230);
    const auto *fsi0_231 = buffer.data(fsi0 + 231);
    const auto *fsi0_232 = buffer.data(fsi0 + 232);
    const auto *fsi0_234 = buffer.data(fsi0 + 234);
    const auto *fsi0_235 = buffer.data(fsi0 + 235);
    const auto *fsi0_236 = buffer.data(fsi0 + 236);

    const auto *fsi1_168 = buffer.data(fsi1 + 168);
    const auto *fsi1_169 = buffer.data(fsi1 + 169);
    const auto *fsi1_171 = buffer.data(fsi1 + 171);
    const auto *fsi1_173 = buffer.data(fsi1 + 173);
    const auto *fsi1_174 = buffer.data(fsi1 + 174);
    const auto *fsi1_176 = buffer.data(fsi1 + 176);
    const auto *fsi1_177 = buffer.data(fsi1 + 177);
    const auto *fsi1_178 = buffer.data(fsi1 + 178);
    const auto *fsi1_180 = buffer.data(fsi1 + 180);
    const auto *fsi1_181 = buffer.data(fsi1 + 181);
    const auto *fsi1_182 = buffer.data(fsi1 + 182);
    const auto *fsi1_183 = buffer.data(fsi1 + 183);
    const auto *fsi1_185 = buffer.data(fsi1 + 185);
    const auto *fsi1_186 = buffer.data(fsi1 + 186);
    const auto *fsi1_187 = buffer.data(fsi1 + 187);
    const auto *fsi1_188 = buffer.data(fsi1 + 188);
    const auto *fsi1_189 = buffer.data(fsi1 + 189);
    const auto *fsi1_190 = buffer.data(fsi1 + 190);
    const auto *fsi1_191 = buffer.data(fsi1 + 191);
    const auto *fsi1_192 = buffer.data(fsi1 + 192);
    const auto *fsi1_193 = buffer.data(fsi1 + 193);
    const auto *fsi1_194 = buffer.data(fsi1 + 194);
    const auto *fsi1_195 = buffer.data(fsi1 + 195);
    const auto *fsi1_198 = buffer.data(fsi1 + 198);
    const auto *fsi1_200 = buffer.data(fsi1 + 200);
    const auto *fsi1_201 = buffer.data(fsi1 + 201);
    const auto *fsi1_203 = buffer.data(fsi1 + 203);
    const auto *fsi1_204 = buffer.data(fsi1 + 204);
    const auto *fsi1_205 = buffer.data(fsi1 + 205);
    const auto *fsi1_207 = buffer.data(fsi1 + 207);
    const auto *fsi1_208 = buffer.data(fsi1 + 208);
    const auto *fsi1_209 = buffer.data(fsi1 + 209);
    const auto *fsi1_210 = buffer.data(fsi1 + 210);
    const auto *fsi1_212 = buffer.data(fsi1 + 212);
    const auto *fsi1_213 = buffer.data(fsi1 + 213);
    const auto *fsi1_214 = buffer.data(fsi1 + 214);
    const auto *fsi1_215 = buffer.data(fsi1 + 215);
    const auto *fsi1_216 = buffer.data(fsi1 + 216);
    const auto *fsi1_218 = buffer.data(fsi1 + 218);
    const auto *fsi1_219 = buffer.data(fsi1 + 219);
    const auto *fsi1_220 = buffer.data(fsi1 + 220);
    const auto *fsi1_221 = buffer.data(fsi1 + 221);
    const auto *fsi1_222 = buffer.data(fsi1 + 222);
    const auto *fsi1_223 = buffer.data(fsi1 + 223);
    const auto *fsi1_225 = buffer.data(fsi1 + 225);
    const auto *fsi1_227 = buffer.data(fsi1 + 227);
    const auto *fsi1_228 = buffer.data(fsi1 + 228);
    const auto *fsi1_230 = buffer.data(fsi1 + 230);
    const auto *fsi1_231 = buffer.data(fsi1 + 231);
    const auto *fsi1_232 = buffer.data(fsi1 + 232);
    const auto *fsi1_234 = buffer.data(fsi1 + 234);
    const auto *fsi1_235 = buffer.data(fsi1 + 235);
    const auto *fsi1_236 = buffer.data(fsi1 + 236);

    const auto *fsk_207 = buffer.data(fsk + 207);
    const auto *fsk_208 = buffer.data(fsk + 208);
    const auto *fsk_209 = buffer.data(fsk + 209);
    const auto *fsk_210 = buffer.data(fsk + 210);
    const auto *fsk_211 = buffer.data(fsk + 211);
    const auto *fsk_212 = buffer.data(fsk + 212);
    const auto *fsk_213 = buffer.data(fsk + 213);
    const auto *fsk_215 = buffer.data(fsk + 215);
    const auto *fsk_216 = buffer.data(fsk + 216);
    const auto *fsk_217 = buffer.data(fsk + 217);
    const auto *fsk_219 = buffer.data(fsk + 219);
    const auto *fsk_221 = buffer.data(fsk + 221);
    const auto *fsk_222 = buffer.data(fsk + 222);
    const auto *fsk_224 = buffer.data(fsk + 224);
    const auto *fsk_225 = buffer.data(fsk + 225);
    const auto *fsk_226 = buffer.data(fsk + 226);
    const auto *fsk_228 = buffer.data(fsk + 228);
    const auto *fsk_229 = buffer.data(fsk + 229);
    const auto *fsk_230 = buffer.data(fsk + 230);
    const auto *fsk_231 = buffer.data(fsk + 231);
    const auto *fsk_233 = buffer.data(fsk + 233);
    const auto *fsk_234 = buffer.data(fsk + 234);
    const auto *fsk_235 = buffer.data(fsk + 235);
    const auto *fsk_236 = buffer.data(fsk + 236);
    const auto *fsk_237 = buffer.data(fsk + 237);
    const auto *fsk_239 = buffer.data(fsk + 239);
    const auto *fsk_240 = buffer.data(fsk + 240);
    const auto *fsk_241 = buffer.data(fsk + 241);
    const auto *fsk_242 = buffer.data(fsk + 242);
    const auto *fsk_243 = buffer.data(fsk + 243);
    const auto *fsk_244 = buffer.data(fsk + 244);
    const auto *fsk_245 = buffer.data(fsk + 245);
    const auto *fsk_246 = buffer.data(fsk + 246);
    const auto *fsk_247 = buffer.data(fsk + 247);
    const auto *fsk_248 = buffer.data(fsk + 248);
    const auto *fsk_249 = buffer.data(fsk + 249);
    const auto *fsk_250 = buffer.data(fsk + 250);
    const auto *fsk_251 = buffer.data(fsk + 251);
    const auto *fsk_254 = buffer.data(fsk + 254);
    const auto *fsk_256 = buffer.data(fsk + 256);
    const auto *fsk_257 = buffer.data(fsk + 257);
    const auto *fsk_259 = buffer.data(fsk + 259);
    const auto *fsk_260 = buffer.data(fsk + 260);
    const auto *fsk_261 = buffer.data(fsk + 261);
    const auto *fsk_263 = buffer.data(fsk + 263);
    const auto *fsk_264 = buffer.data(fsk + 264);
    const auto *fsk_265 = buffer.data(fsk + 265);
    const auto *fsk_266 = buffer.data(fsk + 266);
    const auto *fsk_268 = buffer.data(fsk + 268);
    const auto *fsk_269 = buffer.data(fsk + 269);
    const auto *fsk_270 = buffer.data(fsk + 270);
    const auto *fsk_271 = buffer.data(fsk + 271);
    const auto *fsk_272 = buffer.data(fsk + 272);
    const auto *fsk_274 = buffer.data(fsk + 274);
    const auto *fsk_275 = buffer.data(fsk + 275);
    const auto *fsk_276 = buffer.data(fsk + 276);
    const auto *fsk_277 = buffer.data(fsk + 277);
    const auto *fsk_278 = buffer.data(fsk + 278);
    const auto *fsk_279 = buffer.data(fsk + 279);
    const auto *fsk_280 = buffer.data(fsk + 280);
    const auto *fsk_281 = buffer.data(fsk + 281);
    const auto *fsk_282 = buffer.data(fsk + 282);
    const auto *fsk_283 = buffer.data(fsk + 283);
    const auto *fsk_284 = buffer.data(fsk + 284);
    const auto *fsk_285 = buffer.data(fsk + 285);
    const auto *fsk_286 = buffer.data(fsk + 286);
    const auto *fsk_287 = buffer.data(fsk + 287);
    const auto *fsk_289 = buffer.data(fsk + 289);
    const auto *fsk_291 = buffer.data(fsk + 291);
    const auto *fsk_292 = buffer.data(fsk + 292);
    const auto *fsk_294 = buffer.data(fsk + 294);
    const auto *fsk_295 = buffer.data(fsk + 295);
    const auto *fsk_296 = buffer.data(fsk + 296);
    const auto *fsk_298 = buffer.data(fsk + 298);
    const auto *fsk_299 = buffer.data(fsk + 299);
    const auto *fsk_300 = buffer.data(fsk + 300);

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_x, pc_x, dsl0_252, dsk_207, dsk_208, \
                         dsk_209, dsk_210, dsl1_252, fsk_208, fsk_209, \
                         fsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pa_x[k] * dsl0_252[k]
                   + f_16 * dsk_207[k]
                   - f_14 * pc_x[k] * dsl1_252[k];

        t_253[k] = f_15 * dsk_208[k]
                   + f_3 * pc_x[k] * fsk_208[k];

        t_254[k] = f_15 * dsk_209[k]
                   + f_3 * pc_x[k] * fsk_209[k];

        t_255[k] = f_15 * dsk_210[k]
                   + f_3 * pc_x[k] * fsk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, dsk_211, dsk_212, \
                         dsk_213, dsk_215, fsk_207, fsk_211, fsk_212, fsk_213, \
                         fsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_15 * dsk_211[k]
                   + f_3 * pc_x[k] * fsk_211[k];

        t_257[k] = f_15 * dsk_212[k]
                   + f_3 * pc_x[k] * fsk_212[k];

        t_258[k] = f_15 * dsk_213[k]
                   + f_3 * pc_x[k] * fsk_213[k];

        t_259[k] = f_3 * pc_y[k] * fsk_207[k];

        t_260[k] = f_15 * dsk_215[k]
                   + f_3 * pc_x[k] * fsk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_x, pc_x, dsl0_261, dsl0_262, dsl0_263, \
                         dsl0_264, dsl1_261, dsl1_262, dsl1_263, \
                         dsl1_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_x[k] * dsl0_261[k]
                   - f_14 * pc_x[k] * dsl1_261[k];

        t_262[k] = pa_x[k] * dsl0_262[k]
                   - f_14 * pc_x[k] * dsl1_262[k];

        t_263[k] = pa_x[k] * dsl0_263[k]
                   - f_14 * pc_x[k] * dsl1_263[k];

        t_264[k] = pa_x[k] * dsl0_264[k]
                   - f_14 * pc_x[k] * dsl1_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pa_x, pc_x, pc_y, dsl0_265, dsl0_266, \
                         dsl0_267, dsl1_265, dsl1_266, dsl1_267, \
                         fsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = pa_x[k] * dsl0_265[k]
                   - f_14 * pc_x[k] * dsl1_265[k];

        t_266[k] = pa_x[k] * dsl0_266[k]
                   - f_14 * pc_x[k] * dsl1_266[k];

        t_267[k] = pa_x[k] * dsl0_267[k]
                   - f_14 * pc_x[k] * dsl1_267[k];

        t_268[k] = f_3 * pc_y[k] * fsk_215[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pc_x, pc_z, dsl0_269, dsl1_269, \
                         fsi0_168, fsi0_169, fsi1_168, fsi1_169, fsk_216, \
                         fsk_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_x[k] * dsl0_269[k]
                   - f_14 * pc_x[k] * dsl1_269[k];

        t_270[k] = f_1 * fsi0_168[k]
                   - f_2 * fsi1_168[k]
                   + f_3 * pc_x[k] * fsk_216[k];

        t_271[k] = f_20 * fsi0_169[k]
                   - f_21 * fsi1_169[k]
                   + f_3 * pc_x[k] * fsk_217[k];

        t_272[k] = f_3 * pc_z[k] * fsk_216[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pc_x, pc_z, fsi0_171, fsi0_173, fsi0_174, \
                         fsi1_171, fsi1_173, fsi1_174, fsk_217, fsk_219, fsk_221, \
                         fsk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_12 * fsi0_171[k]
                   - f_13 * fsi1_171[k]
                   + f_3 * pc_x[k] * fsk_219[k];

        t_274[k] = f_3 * pc_z[k] * fsk_217[k];

        t_275[k] = f_12 * fsi0_173[k]
                   - f_13 * fsi1_173[k]
                   + f_3 * pc_x[k] * fsk_221[k];

        t_276[k] = f_10 * fsi0_174[k]
                   - f_11 * fsi1_174[k]
                   + f_3 * pc_x[k] * fsk_222[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pc_x, pc_z, fsi0_176, fsi0_177, fsi0_178, \
                         fsi1_176, fsi1_177, fsi1_178, fsk_219, fsk_224, fsk_225, \
                         fsk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_3 * pc_z[k] * fsk_219[k];

        t_278[k] = f_10 * fsi0_176[k]
                   - f_11 * fsi1_176[k]
                   + f_3 * pc_x[k] * fsk_224[k];

        t_279[k] = f_10 * fsi0_177[k]
                   - f_11 * fsi1_177[k]
                   + f_3 * pc_x[k] * fsk_225[k];

        t_280[k] = f_8 * fsi0_178[k]
                   - f_9 * fsi1_178[k]
                   + f_3 * pc_x[k] * fsk_226[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pc_x, pc_z, fsi0_180, fsi0_181, fsi0_182, \
                         fsi1_180, fsi1_181, fsi1_182, fsk_222, fsk_228, fsk_229, \
                         fsk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * fsk_222[k];

        t_282[k] = f_8 * fsi0_180[k]
                   - f_9 * fsi1_180[k]
                   + f_3 * pc_x[k] * fsk_228[k];

        t_283[k] = f_8 * fsi0_181[k]
                   - f_9 * fsi1_181[k]
                   + f_3 * pc_x[k] * fsk_229[k];

        t_284[k] = f_8 * fsi0_182[k]
                   - f_9 * fsi1_182[k]
                   + f_3 * pc_x[k] * fsk_230[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pc_x, pc_z, fsi0_183, fsi0_185, fsi0_186, \
                         fsi1_183, fsi1_185, fsi1_186, fsk_226, fsk_231, fsk_233, \
                         fsk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_6 * fsi0_183[k]
                   - f_7 * fsi1_183[k]
                   + f_3 * pc_x[k] * fsk_231[k];

        t_286[k] = f_3 * pc_z[k] * fsk_226[k];

        t_287[k] = f_6 * fsi0_185[k]
                   - f_7 * fsi1_185[k]
                   + f_3 * pc_x[k] * fsk_233[k];

        t_288[k] = f_6 * fsi0_186[k]
                   - f_7 * fsi1_186[k]
                   + f_3 * pc_x[k] * fsk_234[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pc_x, pc_z, fsi0_187, fsi0_188, fsi0_189, \
                         fsi1_187, fsi1_188, fsi1_189, fsk_231, fsk_235, fsk_236, \
                         fsk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_6 * fsi0_187[k]
                   - f_7 * fsi1_187[k]
                   + f_3 * pc_x[k] * fsk_235[k];

        t_290[k] = f_6 * fsi0_188[k]
                   - f_7 * fsi1_188[k]
                   + f_3 * pc_x[k] * fsk_236[k];

        t_291[k] = f_4 * fsi0_189[k]
                   - f_5 * fsi1_189[k]
                   + f_3 * pc_x[k] * fsk_237[k];

        t_292[k] = f_3 * pc_z[k] * fsk_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pc_x, fsi0_191, fsi0_192, fsi0_193, fsi1_191, \
                         fsi1_192, fsi1_193, fsk_239, fsk_240, \
                         fsk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_4 * fsi0_191[k]
                   - f_5 * fsi1_191[k]
                   + f_3 * pc_x[k] * fsk_239[k];

        t_294[k] = f_4 * fsi0_192[k]
                   - f_5 * fsi1_192[k]
                   + f_3 * pc_x[k] * fsk_240[k];

        t_295[k] = f_4 * fsi0_193[k]
                   - f_5 * fsi1_193[k]
                   + f_3 * pc_x[k] * fsk_241[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, fsi0_194, fsi0_195, \
                         fsi1_194, fsi1_195, fsk_242, fsk_243, fsk_244, fsk_245, \
                         fsk_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_4 * fsi0_194[k]
                   - f_5 * fsi1_194[k]
                   + f_3 * pc_x[k] * fsk_242[k];

        t_297[k] = f_4 * fsi0_195[k]
                   - f_5 * fsi1_195[k]
                   + f_3 * pc_x[k] * fsk_243[k];

        t_298[k] = f_3 * pc_x[k] * fsk_244[k];

        t_299[k] = f_3 * pc_x[k] * fsk_245[k];

        t_300[k] = f_3 * pc_x[k] * fsk_246[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, t_305, pc_x, fsk_247, fsk_248, fsk_249, \
                         fsk_250, fsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_3 * pc_x[k] * fsk_247[k];

        t_302[k] = f_3 * pc_x[k] * fsk_248[k];

        t_303[k] = f_3 * pc_x[k] * fsk_249[k];

        t_304[k] = f_3 * pc_x[k] * fsk_250[k];

        t_305[k] = f_3 * pc_x[k] * fsk_251[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_y, pc_z, dsk_136, fsi0_189, fsi0_190, \
                         fsi1_189, fsi1_190, fsk_244, fsk_245, \
                         fsk_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * dsk_136[k]
                   + f_1 * fsi0_189[k]
                   - f_2 * fsi1_189[k]
                   + f_3 * pc_y[k] * fsk_244[k];

        t_307[k] = f_3 * pc_z[k] * fsk_244[k];

        t_308[k] = f_4 * fsi0_189[k]
                   - f_5 * fsi1_189[k]
                   + f_3 * pc_z[k] * fsk_245[k];

        t_309[k] = f_6 * fsi0_190[k]
                   - f_7 * fsi1_190[k]
                   + f_3 * pc_z[k] * fsk_246[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pc_z, fsi0_191, fsi0_192, fsi0_193, fsi1_191, \
                         fsi1_192, fsi1_193, fsk_247, fsk_248, \
                         fsk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_8 * fsi0_191[k]
                   - f_9 * fsi1_191[k]
                   + f_3 * pc_z[k] * fsk_247[k];

        t_311[k] = f_10 * fsi0_192[k]
                   - f_11 * fsi1_192[k]
                   + f_3 * pc_z[k] * fsk_248[k];

        t_312[k] = f_12 * fsi0_193[k]
                   - f_13 * fsi1_193[k]
                   + f_3 * pc_z[k] * fsk_249[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_z, pc_y, pc_z, dsl0_135, dsl0_136, \
                         dsk_143, dsl1_135, dsl1_136, fsi0_195, fsi1_195, \
                         fsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_0 * dsk_143[k]
                   + f_3 * pc_y[k] * fsk_251[k];

        t_314[k] = f_1 * fsi0_195[k]
                   - f_2 * fsi1_195[k]
                   + f_3 * pc_z[k] * fsk_251[k];

        t_315[k] = pa_z[k] * dsl0_135[k]
                   - f_14 * pc_z[k] * dsl1_135[k];

        t_316[k] = pa_z[k] * dsl0_136[k]
                   - f_14 * pc_z[k] * dsl1_136[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pa_z, pc_x, pc_z, dsl0_138, dsl1_138, fsi0_198, \
                         fsi0_200, fsi1_198, fsi1_200, fsk_254, \
                         fsk_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_20 * fsi0_198[k]
                   - f_21 * fsi1_198[k]
                   + f_3 * pc_x[k] * fsk_254[k];

        t_318[k] = pa_z[k] * dsl0_138[k]
                   - f_14 * pc_z[k] * dsl1_138[k];

        t_319[k] = f_12 * fsi0_200[k]
                   - f_13 * fsi1_200[k]
                   + f_3 * pc_x[k] * fsk_256[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pa_z, pc_x, pc_z, dsl0_141, dsl1_141, fsi0_201, \
                         fsi0_203, fsi1_201, fsi1_203, fsk_257, \
                         fsk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_12 * fsi0_201[k]
                   - f_13 * fsi1_201[k]
                   + f_3 * pc_x[k] * fsk_257[k];

        t_321[k] = pa_z[k] * dsl0_141[k]
                   - f_14 * pc_z[k] * dsl1_141[k];

        t_322[k] = f_10 * fsi0_203[k]
                   - f_11 * fsi1_203[k]
                   + f_3 * pc_x[k] * fsk_259[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, pa_z, pc_x, pc_z, dsl0_145, dsl1_145, fsi0_204, \
                         fsi0_205, fsi1_204, fsi1_205, fsk_260, \
                         fsk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_10 * fsi0_204[k]
                   - f_11 * fsi1_204[k]
                   + f_3 * pc_x[k] * fsk_260[k];

        t_324[k] = f_10 * fsi0_205[k]
                   - f_11 * fsi1_205[k]
                   + f_3 * pc_x[k] * fsk_261[k];

        t_325[k] = pa_z[k] * dsl0_145[k]
                   - f_14 * pc_z[k] * dsl1_145[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, pc_x, fsi0_207, fsi0_208, fsi0_209, fsi1_207, \
                         fsi1_208, fsi1_209, fsk_263, fsk_264, \
                         fsk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_8 * fsi0_207[k]
                   - f_9 * fsi1_207[k]
                   + f_3 * pc_x[k] * fsk_263[k];

        t_327[k] = f_8 * fsi0_208[k]
                   - f_9 * fsi1_208[k]
                   + f_3 * pc_x[k] * fsk_264[k];

        t_328[k] = f_8 * fsi0_209[k]
                   - f_9 * fsi1_209[k]
                   + f_3 * pc_x[k] * fsk_265[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pc_x, pc_z, dsl0_150, dsl1_150, fsi0_210, \
                         fsi0_212, fsi1_210, fsi1_212, fsk_266, \
                         fsk_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_8 * fsi0_210[k]
                   - f_9 * fsi1_210[k]
                   + f_3 * pc_x[k] * fsk_266[k];

        t_330[k] = pa_z[k] * dsl0_150[k]
                   - f_14 * pc_z[k] * dsl1_150[k];

        t_331[k] = f_6 * fsi0_212[k]
                   - f_7 * fsi1_212[k]
                   + f_3 * pc_x[k] * fsk_268[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pc_x, fsi0_213, fsi0_214, fsi0_215, fsi1_213, \
                         fsi1_214, fsi1_215, fsk_269, fsk_270, \
                         fsk_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_6 * fsi0_213[k]
                   - f_7 * fsi1_213[k]
                   + f_3 * pc_x[k] * fsk_269[k];

        t_333[k] = f_6 * fsi0_214[k]
                   - f_7 * fsi1_214[k]
                   + f_3 * pc_x[k] * fsk_270[k];

        t_334[k] = f_6 * fsi0_215[k]
                   - f_7 * fsi1_215[k]
                   + f_3 * pc_x[k] * fsk_271[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_x, pc_z, dsl0_156, dsl1_156, fsi0_216, \
                         fsi0_218, fsi1_216, fsi1_218, fsk_272, \
                         fsk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_6 * fsi0_216[k]
                   - f_7 * fsi1_216[k]
                   + f_3 * pc_x[k] * fsk_272[k];

        t_336[k] = pa_z[k] * dsl0_156[k]
                   - f_14 * pc_z[k] * dsl1_156[k];

        t_337[k] = f_4 * fsi0_218[k]
                   - f_5 * fsi1_218[k]
                   + f_3 * pc_x[k] * fsk_274[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pc_x, fsi0_219, fsi0_220, fsi0_221, fsi1_219, \
                         fsi1_220, fsi1_221, fsk_275, fsk_276, \
                         fsk_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_4 * fsi0_219[k]
                   - f_5 * fsi1_219[k]
                   + f_3 * pc_x[k] * fsk_275[k];

        t_339[k] = f_4 * fsi0_220[k]
                   - f_5 * fsi1_220[k]
                   + f_3 * pc_x[k] * fsk_276[k];

        t_340[k] = f_4 * fsi0_221[k]
                   - f_5 * fsi1_221[k]
                   + f_3 * pc_x[k] * fsk_277[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, pc_x, fsi0_222, fsi0_223, \
                         fsi1_222, fsi1_223, fsk_278, fsk_279, fsk_280, fsk_281, \
                         fsk_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_4 * fsi0_222[k]
                   - f_5 * fsi1_222[k]
                   + f_3 * pc_x[k] * fsk_278[k];

        t_342[k] = f_4 * fsi0_223[k]
                   - f_5 * fsi1_223[k]
                   + f_3 * pc_x[k] * fsk_279[k];

        t_343[k] = f_3 * pc_x[k] * fsk_280[k];

        t_344[k] = f_3 * pc_x[k] * fsk_281[k];

        t_345[k] = f_3 * pc_x[k] * fsk_282[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, t_351, pa_z, pc_x, pc_z, dsl0_171, \
                         dsl1_171, fsk_283, fsk_284, fsk_285, fsk_286, \
                         fsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_3 * pc_x[k] * fsk_283[k];

        t_347[k] = f_3 * pc_x[k] * fsk_284[k];

        t_348[k] = f_3 * pc_x[k] * fsk_285[k];

        t_349[k] = f_3 * pc_x[k] * fsk_286[k];

        t_350[k] = f_3 * pc_x[k] * fsk_287[k];

        t_351[k] = pa_z[k] * dsl0_171[k]
                   - f_14 * pc_z[k] * dsl1_171[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pa_z, pc_z, dsl0_173, dsl0_174, dsk_136, \
                         dsk_137, dsk_138, dsl1_173, dsl1_174, \
                         fsk_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_15 * dsk_136[k]
                   + f_3 * pc_z[k] * fsk_280[k];

        t_353[k] = pa_z[k] * dsl0_173[k]
                   + f_16 * dsk_137[k]
                   - f_14 * pc_z[k] * dsl1_173[k];

        t_354[k] = pa_z[k] * dsl0_174[k]
                   + f_0 * dsk_138[k]
                   - f_14 * pc_z[k] * dsl1_174[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pa_z, pc_z, dsl0_175, dsl0_176, dsl0_177, \
                         dsk_139, dsk_140, dsk_141, dsl1_175, dsl1_176, \
                         dsl1_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = pa_z[k] * dsl0_175[k]
                   + f_17 * dsk_139[k]
                   - f_14 * pc_z[k] * dsl1_175[k];

        t_356[k] = pa_z[k] * dsl0_176[k]
                   + f_18 * dsk_140[k]
                   - f_14 * pc_z[k] * dsl1_176[k];

        t_357[k] = pa_z[k] * dsl0_177[k]
                   + f_19 * dsk_141[k]
                   - f_14 * pc_z[k] * dsl1_177[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pa_y, pc_y, pc_z, dsl0_225, dsk_143, dsk_179, \
                         dsl1_225, fsi0_223, fsi1_223, fsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_16 * dsk_179[k]
                   + f_3 * pc_y[k] * fsk_287[k];

        t_359[k] = f_15 * dsk_143[k]
                   + f_1 * fsi0_223[k]
                   - f_2 * fsi1_223[k]
                   + f_3 * pc_z[k] * fsk_287[k];

        t_360[k] = pa_y[k] * dsl0_225[k]
                   - f_14 * pc_y[k] * dsl1_225[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pa_y, pc_x, pc_y, dsl0_227, dsl1_227, fsi0_225, \
                         fsi0_227, fsi1_225, fsi1_227, fsk_289, \
                         fsk_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_20 * fsi0_225[k]
                   - f_21 * fsi1_225[k]
                   + f_3 * pc_x[k] * fsk_289[k];

        t_362[k] = pa_y[k] * dsl0_227[k]
                   - f_14 * pc_y[k] * dsl1_227[k];

        t_363[k] = f_12 * fsi0_227[k]
                   - f_13 * fsi1_227[k]
                   + f_3 * pc_x[k] * fsk_291[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pa_y, pc_x, pc_y, dsl0_230, dsl1_230, fsi0_228, \
                         fsi0_230, fsi1_228, fsi1_230, fsk_292, \
                         fsk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_12 * fsi0_228[k]
                   - f_13 * fsi1_228[k]
                   + f_3 * pc_x[k] * fsk_292[k];

        t_365[k] = pa_y[k] * dsl0_230[k]
                   - f_14 * pc_y[k] * dsl1_230[k];

        t_366[k] = f_10 * fsi0_230[k]
                   - f_11 * fsi1_230[k]
                   + f_3 * pc_x[k] * fsk_294[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pa_y, pc_x, pc_y, dsl0_234, dsl1_234, fsi0_231, \
                         fsi0_232, fsi1_231, fsi1_232, fsk_295, \
                         fsk_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_10 * fsi0_231[k]
                   - f_11 * fsi1_231[k]
                   + f_3 * pc_x[k] * fsk_295[k];

        t_368[k] = f_10 * fsi0_232[k]
                   - f_11 * fsi1_232[k]
                   + f_3 * pc_x[k] * fsk_296[k];

        t_369[k] = pa_y[k] * dsl0_234[k]
                   - f_14 * pc_y[k] * dsl1_234[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, fsi0_234, fsi0_235, fsi0_236, fsi1_234, \
                         fsi1_235, fsi1_236, fsk_298, fsk_299, \
                         fsk_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_8 * fsi0_234[k]
                   - f_9 * fsi1_234[k]
                   + f_3 * pc_x[k] * fsk_298[k];

        t_371[k] = f_8 * fsi0_235[k]
                   - f_9 * fsi1_235[k]
                   + f_3 * pc_x[k] * fsk_299[k];

        t_372[k] = f_8 * fsi0_236[k]
                   - f_9 * fsi1_236[k]
                   + f_3 * pc_x[k] * fsk_300[k];
    }
}

static auto
compute_prim_fsl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsl0,
                                                          const size_t dsk, const size_t dsl1,
                                                          const size_t fsi0, const size_t fsi1,
                                                          const size_t fsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsl0_239 = buffer.data(dsl0 + 239);
    const auto *dsl0_245 = buffer.data(dsl0 + 245);
    const auto *dsl0_252 = buffer.data(dsl0 + 252);
    const auto *dsl0_261 = buffer.data(dsl0 + 261);
    const auto *dsl0_263 = buffer.data(dsl0 + 263);
    const auto *dsl0_264 = buffer.data(dsl0 + 264);
    const auto *dsl0_265 = buffer.data(dsl0 + 265);
    const auto *dsl0_266 = buffer.data(dsl0 + 266);
    const auto *dsl0_267 = buffer.data(dsl0 + 267);
    const auto *dsl0_269 = buffer.data(dsl0 + 269);

    const auto *dsk_172 = buffer.data(dsk + 172);
    const auto *dsk_208 = buffer.data(dsk + 208);
    const auto *dsk_210 = buffer.data(dsk + 210);
    const auto *dsk_211 = buffer.data(dsk + 211);
    const auto *dsk_212 = buffer.data(dsk + 212);
    const auto *dsk_213 = buffer.data(dsk + 213);
    const auto *dsk_214 = buffer.data(dsk + 214);
    const auto *dsk_215 = buffer.data(dsk + 215);

    const auto *dsl1_239 = buffer.data(dsl1 + 239);
    const auto *dsl1_245 = buffer.data(dsl1 + 245);
    const auto *dsl1_252 = buffer.data(dsl1 + 252);
    const auto *dsl1_261 = buffer.data(dsl1 + 261);
    const auto *dsl1_263 = buffer.data(dsl1 + 263);
    const auto *dsl1_264 = buffer.data(dsl1 + 264);
    const auto *dsl1_265 = buffer.data(dsl1 + 265);
    const auto *dsl1_266 = buffer.data(dsl1 + 266);
    const auto *dsl1_267 = buffer.data(dsl1 + 267);
    const auto *dsl1_269 = buffer.data(dsl1 + 269);

    const auto *fsi0_237 = buffer.data(fsi0 + 237);
    const auto *fsi0_239 = buffer.data(fsi0 + 239);
    const auto *fsi0_240 = buffer.data(fsi0 + 240);
    const auto *fsi0_241 = buffer.data(fsi0 + 241);
    const auto *fsi0_242 = buffer.data(fsi0 + 242);
    const auto *fsi0_243 = buffer.data(fsi0 + 243);
    const auto *fsi0_245 = buffer.data(fsi0 + 245);
    const auto *fsi0_246 = buffer.data(fsi0 + 246);
    const auto *fsi0_247 = buffer.data(fsi0 + 247);
    const auto *fsi0_248 = buffer.data(fsi0 + 248);
    const auto *fsi0_249 = buffer.data(fsi0 + 249);
    const auto *fsi0_250 = buffer.data(fsi0 + 250);
    const auto *fsi0_252 = buffer.data(fsi0 + 252);
    const auto *fsi0_254 = buffer.data(fsi0 + 254);
    const auto *fsi0_255 = buffer.data(fsi0 + 255);
    const auto *fsi0_257 = buffer.data(fsi0 + 257);
    const auto *fsi0_258 = buffer.data(fsi0 + 258);
    const auto *fsi0_259 = buffer.data(fsi0 + 259);
    const auto *fsi0_261 = buffer.data(fsi0 + 261);
    const auto *fsi0_262 = buffer.data(fsi0 + 262);
    const auto *fsi0_263 = buffer.data(fsi0 + 263);
    const auto *fsi0_264 = buffer.data(fsi0 + 264);
    const auto *fsi0_266 = buffer.data(fsi0 + 266);
    const auto *fsi0_267 = buffer.data(fsi0 + 267);
    const auto *fsi0_268 = buffer.data(fsi0 + 268);
    const auto *fsi0_269 = buffer.data(fsi0 + 269);
    const auto *fsi0_270 = buffer.data(fsi0 + 270);
    const auto *fsi0_272 = buffer.data(fsi0 + 272);
    const auto *fsi0_273 = buffer.data(fsi0 + 273);
    const auto *fsi0_274 = buffer.data(fsi0 + 274);
    const auto *fsi0_275 = buffer.data(fsi0 + 275);
    const auto *fsi0_276 = buffer.data(fsi0 + 276);
    const auto *fsi0_277 = buffer.data(fsi0 + 277);
    const auto *fsi0_278 = buffer.data(fsi0 + 278);
    const auto *fsi0_279 = buffer.data(fsi0 + 279);

    const auto *fsi1_237 = buffer.data(fsi1 + 237);
    const auto *fsi1_239 = buffer.data(fsi1 + 239);
    const auto *fsi1_240 = buffer.data(fsi1 + 240);
    const auto *fsi1_241 = buffer.data(fsi1 + 241);
    const auto *fsi1_242 = buffer.data(fsi1 + 242);
    const auto *fsi1_243 = buffer.data(fsi1 + 243);
    const auto *fsi1_245 = buffer.data(fsi1 + 245);
    const auto *fsi1_246 = buffer.data(fsi1 + 246);
    const auto *fsi1_247 = buffer.data(fsi1 + 247);
    const auto *fsi1_248 = buffer.data(fsi1 + 248);
    const auto *fsi1_249 = buffer.data(fsi1 + 249);
    const auto *fsi1_250 = buffer.data(fsi1 + 250);
    const auto *fsi1_252 = buffer.data(fsi1 + 252);
    const auto *fsi1_254 = buffer.data(fsi1 + 254);
    const auto *fsi1_255 = buffer.data(fsi1 + 255);
    const auto *fsi1_257 = buffer.data(fsi1 + 257);
    const auto *fsi1_258 = buffer.data(fsi1 + 258);
    const auto *fsi1_259 = buffer.data(fsi1 + 259);
    const auto *fsi1_261 = buffer.data(fsi1 + 261);
    const auto *fsi1_262 = buffer.data(fsi1 + 262);
    const auto *fsi1_263 = buffer.data(fsi1 + 263);
    const auto *fsi1_264 = buffer.data(fsi1 + 264);
    const auto *fsi1_266 = buffer.data(fsi1 + 266);
    const auto *fsi1_267 = buffer.data(fsi1 + 267);
    const auto *fsi1_268 = buffer.data(fsi1 + 268);
    const auto *fsi1_269 = buffer.data(fsi1 + 269);
    const auto *fsi1_270 = buffer.data(fsi1 + 270);
    const auto *fsi1_272 = buffer.data(fsi1 + 272);
    const auto *fsi1_273 = buffer.data(fsi1 + 273);
    const auto *fsi1_274 = buffer.data(fsi1 + 274);
    const auto *fsi1_275 = buffer.data(fsi1 + 275);
    const auto *fsi1_276 = buffer.data(fsi1 + 276);
    const auto *fsi1_277 = buffer.data(fsi1 + 277);
    const auto *fsi1_278 = buffer.data(fsi1 + 278);
    const auto *fsi1_279 = buffer.data(fsi1 + 279);

    const auto *fsk_301 = buffer.data(fsk + 301);
    const auto *fsk_303 = buffer.data(fsk + 303);
    const auto *fsk_304 = buffer.data(fsk + 304);
    const auto *fsk_305 = buffer.data(fsk + 305);
    const auto *fsk_306 = buffer.data(fsk + 306);
    const auto *fsk_307 = buffer.data(fsk + 307);
    const auto *fsk_309 = buffer.data(fsk + 309);
    const auto *fsk_310 = buffer.data(fsk + 310);
    const auto *fsk_311 = buffer.data(fsk + 311);
    const auto *fsk_312 = buffer.data(fsk + 312);
    const auto *fsk_313 = buffer.data(fsk + 313);
    const auto *fsk_314 = buffer.data(fsk + 314);
    const auto *fsk_316 = buffer.data(fsk + 316);
    const auto *fsk_317 = buffer.data(fsk + 317);
    const auto *fsk_318 = buffer.data(fsk + 318);
    const auto *fsk_319 = buffer.data(fsk + 319);
    const auto *fsk_320 = buffer.data(fsk + 320);
    const auto *fsk_321 = buffer.data(fsk + 321);
    const auto *fsk_322 = buffer.data(fsk + 322);
    const auto *fsk_323 = buffer.data(fsk + 323);
    const auto *fsk_324 = buffer.data(fsk + 324);
    const auto *fsk_326 = buffer.data(fsk + 326);
    const auto *fsk_327 = buffer.data(fsk + 327);
    const auto *fsk_329 = buffer.data(fsk + 329);
    const auto *fsk_330 = buffer.data(fsk + 330);
    const auto *fsk_331 = buffer.data(fsk + 331);
    const auto *fsk_333 = buffer.data(fsk + 333);
    const auto *fsk_334 = buffer.data(fsk + 334);
    const auto *fsk_335 = buffer.data(fsk + 335);
    const auto *fsk_336 = buffer.data(fsk + 336);
    const auto *fsk_338 = buffer.data(fsk + 338);
    const auto *fsk_339 = buffer.data(fsk + 339);
    const auto *fsk_340 = buffer.data(fsk + 340);
    const auto *fsk_341 = buffer.data(fsk + 341);
    const auto *fsk_342 = buffer.data(fsk + 342);
    const auto *fsk_344 = buffer.data(fsk + 344);
    const auto *fsk_345 = buffer.data(fsk + 345);
    const auto *fsk_346 = buffer.data(fsk + 346);
    const auto *fsk_347 = buffer.data(fsk + 347);
    const auto *fsk_348 = buffer.data(fsk + 348);
    const auto *fsk_349 = buffer.data(fsk + 349);
    const auto *fsk_351 = buffer.data(fsk + 351);
    const auto *fsk_352 = buffer.data(fsk + 352);
    const auto *fsk_353 = buffer.data(fsk + 353);
    const auto *fsk_354 = buffer.data(fsk + 354);
    const auto *fsk_355 = buffer.data(fsk + 355);
    const auto *fsk_356 = buffer.data(fsk + 356);
    const auto *fsk_357 = buffer.data(fsk + 357);
    const auto *fsk_358 = buffer.data(fsk + 358);
    const auto *fsk_359 = buffer.data(fsk + 359);

#pragma omp simd aligned(t_373, t_374, t_375, pa_y, pc_x, pc_y, dsl0_239, dsl1_239, fsi0_237, \
                         fsi0_239, fsi1_237, fsi1_239, fsk_301, \
                         fsk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_8 * fsi0_237[k]
                   - f_9 * fsi1_237[k]
                   + f_3 * pc_x[k] * fsk_301[k];

        t_374[k] = pa_y[k] * dsl0_239[k]
                   - f_14 * pc_y[k] * dsl1_239[k];

        t_375[k] = f_6 * fsi0_239[k]
                   - f_7 * fsi1_239[k]
                   + f_3 * pc_x[k] * fsk_303[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, fsi0_240, fsi0_241, fsi0_242, fsi1_240, \
                         fsi1_241, fsi1_242, fsk_304, fsk_305, \
                         fsk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_6 * fsi0_240[k]
                   - f_7 * fsi1_240[k]
                   + f_3 * pc_x[k] * fsk_304[k];

        t_377[k] = f_6 * fsi0_241[k]
                   - f_7 * fsi1_241[k]
                   + f_3 * pc_x[k] * fsk_305[k];

        t_378[k] = f_6 * fsi0_242[k]
                   - f_7 * fsi1_242[k]
                   + f_3 * pc_x[k] * fsk_306[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, pa_y, pc_x, pc_y, dsl0_245, dsl1_245, fsi0_243, \
                         fsi0_245, fsi1_243, fsi1_245, fsk_307, \
                         fsk_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_6 * fsi0_243[k]
                   - f_7 * fsi1_243[k]
                   + f_3 * pc_x[k] * fsk_307[k];

        t_380[k] = pa_y[k] * dsl0_245[k]
                   - f_14 * pc_y[k] * dsl1_245[k];

        t_381[k] = f_4 * fsi0_245[k]
                   - f_5 * fsi1_245[k]
                   + f_3 * pc_x[k] * fsk_309[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pc_x, fsi0_246, fsi0_247, fsi0_248, fsi1_246, \
                         fsi1_247, fsi1_248, fsk_310, fsk_311, \
                         fsk_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_4 * fsi0_246[k]
                   - f_5 * fsi1_246[k]
                   + f_3 * pc_x[k] * fsk_310[k];

        t_383[k] = f_4 * fsi0_247[k]
                   - f_5 * fsi1_247[k]
                   + f_3 * pc_x[k] * fsk_311[k];

        t_384[k] = f_4 * fsi0_248[k]
                   - f_5 * fsi1_248[k]
                   + f_3 * pc_x[k] * fsk_312[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pc_x, pc_y, dsl0_252, dsl1_252, \
                         fsi0_249, fsi0_250, fsi1_249, fsi1_250, fsk_313, fsk_314, \
                         fsk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_4 * fsi0_249[k]
                   - f_5 * fsi1_249[k]
                   + f_3 * pc_x[k] * fsk_313[k];

        t_386[k] = f_4 * fsi0_250[k]
                   - f_5 * fsi1_250[k]
                   + f_3 * pc_x[k] * fsk_314[k];

        t_387[k] = pa_y[k] * dsl0_252[k]
                   - f_14 * pc_y[k] * dsl1_252[k];

        t_388[k] = f_3 * pc_x[k] * fsk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, t_394, t_395, pc_x, fsk_317, \
                         fsk_318, fsk_319, fsk_320, fsk_321, fsk_322, \
                         fsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_3 * pc_x[k] * fsk_317[k];

        t_390[k] = f_3 * pc_x[k] * fsk_318[k];

        t_391[k] = f_3 * pc_x[k] * fsk_319[k];

        t_392[k] = f_3 * pc_x[k] * fsk_320[k];

        t_393[k] = f_3 * pc_x[k] * fsk_321[k];

        t_394[k] = f_3 * pc_x[k] * fsk_322[k];

        t_395[k] = f_3 * pc_x[k] * fsk_323[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, pa_y, pc_y, pc_z, dsl0_261, dsl0_263, dsk_172, \
                         dsk_208, dsk_210, dsl1_261, dsl1_263, \
                         fsk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pa_y[k] * dsl0_261[k]
                   + f_22 * dsk_208[k]
                   - f_14 * pc_y[k] * dsl1_261[k];

        t_397[k] = f_16 * dsk_172[k]
                   + f_3 * pc_z[k] * fsk_316[k];

        t_398[k] = pa_y[k] * dsl0_263[k]
                   + f_19 * dsk_210[k]
                   - f_14 * pc_y[k] * dsl1_263[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pa_y, pc_y, dsl0_264, dsl0_265, dsl0_266, \
                         dsk_211, dsk_212, dsk_213, dsl1_264, dsl1_265, \
                         dsl1_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pa_y[k] * dsl0_264[k]
                   + f_18 * dsk_211[k]
                   - f_14 * pc_y[k] * dsl1_264[k];

        t_400[k] = pa_y[k] * dsl0_265[k]
                   + f_17 * dsk_212[k]
                   - f_14 * pc_y[k] * dsl1_265[k];

        t_401[k] = pa_y[k] * dsl0_266[k]
                   + f_0 * dsk_213[k]
                   - f_14 * pc_y[k] * dsl1_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, dsl0_267, dsl0_269, dsk_214, \
                         dsk_215, dsl1_267, dsl1_269, fsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_y[k] * dsl0_267[k]
                   + f_16 * dsk_214[k]
                   - f_14 * pc_y[k] * dsl1_267[k];

        t_403[k] = f_15 * dsk_215[k]
                   + f_3 * pc_y[k] * fsk_323[k];

        t_404[k] = pa_y[k] * dsl0_269[k]
                   - f_14 * pc_y[k] * dsl1_269[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, fsi0_252, fsi0_254, \
                         fsi0_255, fsi1_252, fsi1_254, fsi1_255, fsk_324, fsk_326, \
                         fsk_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_1 * fsi0_252[k]
                   - f_2 * fsi1_252[k]
                   + f_3 * pc_x[k] * fsk_324[k];

        t_406[k] = f_3 * pc_y[k] * fsk_324[k];

        t_407[k] = f_20 * fsi0_254[k]
                   - f_21 * fsi1_254[k]
                   + f_3 * pc_x[k] * fsk_326[k];

        t_408[k] = f_12 * fsi0_255[k]
                   - f_13 * fsi1_255[k]
                   + f_3 * pc_x[k] * fsk_327[k];

        t_409[k] = f_3 * pc_y[k] * fsk_326[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, fsi0_257, fsi0_258, fsi0_259, \
                         fsi1_257, fsi1_258, fsi1_259, fsk_329, fsk_330, \
                         fsk_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_12 * fsi0_257[k]
                   - f_13 * fsi1_257[k]
                   + f_3 * pc_x[k] * fsk_329[k];

        t_411[k] = f_10 * fsi0_258[k]
                   - f_11 * fsi1_258[k]
                   + f_3 * pc_x[k] * fsk_330[k];

        t_412[k] = f_10 * fsi0_259[k]
                   - f_11 * fsi1_259[k]
                   + f_3 * pc_x[k] * fsk_331[k];

        t_413[k] = f_3 * pc_y[k] * fsk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, fsi0_261, fsi0_262, fsi0_263, fsi1_261, \
                         fsi1_262, fsi1_263, fsk_333, fsk_334, \
                         fsk_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_10 * fsi0_261[k]
                   - f_11 * fsi1_261[k]
                   + f_3 * pc_x[k] * fsk_333[k];

        t_415[k] = f_8 * fsi0_262[k]
                   - f_9 * fsi1_262[k]
                   + f_3 * pc_x[k] * fsk_334[k];

        t_416[k] = f_8 * fsi0_263[k]
                   - f_9 * fsi1_263[k]
                   + f_3 * pc_x[k] * fsk_335[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pc_x, pc_y, fsi0_264, fsi0_266, fsi0_267, \
                         fsi1_264, fsi1_266, fsi1_267, fsk_333, fsk_336, fsk_338, \
                         fsk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_8 * fsi0_264[k]
                   - f_9 * fsi1_264[k]
                   + f_3 * pc_x[k] * fsk_336[k];

        t_418[k] = f_3 * pc_y[k] * fsk_333[k];

        t_419[k] = f_8 * fsi0_266[k]
                   - f_9 * fsi1_266[k]
                   + f_3 * pc_x[k] * fsk_338[k];

        t_420[k] = f_6 * fsi0_267[k]
                   - f_7 * fsi1_267[k]
                   + f_3 * pc_x[k] * fsk_339[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, fsi0_268, fsi0_269, fsi0_270, \
                         fsi1_268, fsi1_269, fsi1_270, fsk_338, fsk_340, fsk_341, \
                         fsk_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_6 * fsi0_268[k]
                   - f_7 * fsi1_268[k]
                   + f_3 * pc_x[k] * fsk_340[k];

        t_422[k] = f_6 * fsi0_269[k]
                   - f_7 * fsi1_269[k]
                   + f_3 * pc_x[k] * fsk_341[k];

        t_423[k] = f_6 * fsi0_270[k]
                   - f_7 * fsi1_270[k]
                   + f_3 * pc_x[k] * fsk_342[k];

        t_424[k] = f_3 * pc_y[k] * fsk_338[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, pc_x, fsi0_272, fsi0_273, fsi0_274, fsi1_272, \
                         fsi1_273, fsi1_274, fsk_344, fsk_345, \
                         fsk_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_6 * fsi0_272[k]
                   - f_7 * fsi1_272[k]
                   + f_3 * pc_x[k] * fsk_344[k];

        t_426[k] = f_4 * fsi0_273[k]
                   - f_5 * fsi1_273[k]
                   + f_3 * pc_x[k] * fsk_345[k];

        t_427[k] = f_4 * fsi0_274[k]
                   - f_5 * fsi1_274[k]
                   + f_3 * pc_x[k] * fsk_346[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, fsi0_275, fsi0_276, fsi0_277, \
                         fsi1_275, fsi1_276, fsi1_277, fsk_344, fsk_347, fsk_348, \
                         fsk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_4 * fsi0_275[k]
                   - f_5 * fsi1_275[k]
                   + f_3 * pc_x[k] * fsk_347[k];

        t_429[k] = f_4 * fsi0_276[k]
                   - f_5 * fsi1_276[k]
                   + f_3 * pc_x[k] * fsk_348[k];

        t_430[k] = f_4 * fsi0_277[k]
                   - f_5 * fsi1_277[k]
                   + f_3 * pc_x[k] * fsk_349[k];

        t_431[k] = f_3 * pc_y[k] * fsk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, t_437, pc_x, fsi0_279, fsi1_279, \
                         fsk_351, fsk_352, fsk_353, fsk_354, fsk_355, \
                         fsk_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * fsi0_279[k]
                   - f_5 * fsi1_279[k]
                   + f_3 * pc_x[k] * fsk_351[k];

        t_433[k] = f_3 * pc_x[k] * fsk_352[k];

        t_434[k] = f_3 * pc_x[k] * fsk_353[k];

        t_435[k] = f_3 * pc_x[k] * fsk_354[k];

        t_436[k] = f_3 * pc_x[k] * fsk_355[k];

        t_437[k] = f_3 * pc_x[k] * fsk_356[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pc_x, pc_y, fsi0_273, fsi0_274, \
                         fsi1_273, fsi1_274, fsk_352, fsk_353, fsk_357, fsk_358, \
                         fsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_3 * pc_x[k] * fsk_357[k];

        t_439[k] = f_3 * pc_x[k] * fsk_358[k];

        t_440[k] = f_3 * pc_x[k] * fsk_359[k];

        t_441[k] = f_1 * fsi0_273[k]
                   - f_2 * fsi1_273[k]
                   + f_3 * pc_y[k] * fsk_352[k];

        t_442[k] = f_20 * fsi0_274[k]
                   - f_21 * fsi1_274[k]
                   + f_3 * pc_y[k] * fsk_353[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_y, fsi0_275, fsi0_276, fsi0_277, fsi1_275, \
                         fsi1_276, fsi1_277, fsk_354, fsk_355, \
                         fsk_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_12 * fsi0_275[k]
                   - f_13 * fsi1_275[k]
                   + f_3 * pc_y[k] * fsk_354[k];

        t_444[k] = f_10 * fsi0_276[k]
                   - f_11 * fsi1_276[k]
                   + f_3 * pc_y[k] * fsk_355[k];

        t_445[k] = f_8 * fsi0_277[k]
                   - f_9 * fsi1_277[k]
                   + f_3 * pc_y[k] * fsk_356[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pc_y, pc_z, dsk_215, fsi0_278, fsi0_279, \
                         fsi1_278, fsi1_279, fsk_357, fsk_358, \
                         fsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_6 * fsi0_278[k]
                   - f_7 * fsi1_278[k]
                   + f_3 * pc_y[k] * fsk_357[k];

        t_447[k] = f_4 * fsi0_279[k]
                   - f_5 * fsi1_279[k]
                   + f_3 * pc_y[k] * fsk_358[k];

        t_448[k] = f_3 * pc_y[k] * fsk_359[k];

        t_449[k] = f_0 * dsk_215[k]
                   + f_1 * fsi0_279[k]
                   - f_2 * fsi1_279[k]
                   + f_3 * pc_z[k] * fsk_359[k];
    }
}

auto
compute_prim_fsl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsl0, const size_t dsk,
                                                   const size_t dsl1, const size_t fsi0,
                                                   const size_t fsi1, const size_t fsk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fsl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, dsl0, dsk,
                                                              dsl1, fsi0, fsi1, fsk, ncols,
                                                              gamma, p, q);

    compute_prim_fsl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, dsl0, dsk,
                                                              dsl1, fsi0, fsi1, fsk, ncols,
                                                              gamma, p, q);

    compute_prim_fsl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, dsl0, dsk,
                                                              dsl1, fsi0, fsi1, fsk, ncols,
                                                              gamma, p, q);

    compute_prim_fsl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, dsl0, dsk,
                                                              dsl1, fsi0, fsi1, fsk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
