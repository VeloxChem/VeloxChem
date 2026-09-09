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


#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_0 = buffer.data(nsi0 + 0);
    const auto *nsi0_3 = buffer.data(nsi0 + 3);
    const auto *nsi0_5 = buffer.data(nsi0 + 5);
    const auto *nsi0_6 = buffer.data(nsi0 + 6);
    const auto *nsi0_9 = buffer.data(nsi0 + 9);
    const auto *nsi0_10 = buffer.data(nsi0 + 10);
    const auto *nsi0_14 = buffer.data(nsi0 + 14);
    const auto *nsi0_21 = buffer.data(nsi0 + 21);
    const auto *nsi0_27 = buffer.data(nsi0 + 27);
    const auto *nsi0_31 = buffer.data(nsi0 + 31);
    const auto *nsi0_34 = buffer.data(nsi0 + 34);
    const auto *nsi0_38 = buffer.data(nsi0 + 38);
    const auto *nsi0_56 = buffer.data(nsi0 + 56);
    const auto *nsi0_61 = buffer.data(nsi0 + 61);
    const auto *nsi0_65 = buffer.data(nsi0 + 65);
    const auto *nsi0_68 = buffer.data(nsi0 + 68);
    const auto *nsi0_70 = buffer.data(nsi0 + 70);

    const auto *nsh_0 = buffer.data(nsh + 0);
    const auto *nsh_1 = buffer.data(nsh + 1);
    const auto *nsh_2 = buffer.data(nsh + 2);
    const auto *nsh_3 = buffer.data(nsh + 3);
    const auto *nsh_5 = buffer.data(nsh + 5);
    const auto *nsh_6 = buffer.data(nsh + 6);
    const auto *nsh_9 = buffer.data(nsh + 9);
    const auto *nsh_15 = buffer.data(nsh + 15);
    const auto *nsh_17 = buffer.data(nsh + 17);
    const auto *nsh_18 = buffer.data(nsh + 18);
    const auto *nsh_20 = buffer.data(nsh + 20);
    const auto *nsh_21 = buffer.data(nsh + 21);
    const auto *nsh_24 = buffer.data(nsh + 24);
    const auto *nsh_26 = buffer.data(nsh + 26);
    const auto *nsh_27 = buffer.data(nsh + 27);
    const auto *nsh_30 = buffer.data(nsh + 30);
    const auto *nsh_36 = buffer.data(nsh + 36);
    const auto *nsh_38 = buffer.data(nsh + 38);
    const auto *nsh_39 = buffer.data(nsh + 39);
    const auto *nsh_40 = buffer.data(nsh + 40);
    const auto *nsh_41 = buffer.data(nsh + 41);
    const auto *nsh_42 = buffer.data(nsh + 42);
    const auto *nsh_44 = buffer.data(nsh + 44);
    const auto *nsh_47 = buffer.data(nsh + 47);
    const auto *nsh_50 = buffer.data(nsh + 50);
    const auto *nsh_51 = buffer.data(nsh + 51);
    const auto *nsh_57 = buffer.data(nsh + 57);
    const auto *nsh_58 = buffer.data(nsh + 58);
    const auto *nsh_59 = buffer.data(nsh + 59);
    const auto *nsh_60 = buffer.data(nsh + 60);
    const auto *nsh_62 = buffer.data(nsh + 62);
    const auto *nsh_63 = buffer.data(nsh + 63);
    const auto *nsh_66 = buffer.data(nsh + 66);
    const auto *nsh_69 = buffer.data(nsh + 69);
    const auto *nsh_73 = buffer.data(nsh + 73);
    const auto *nsh_78 = buffer.data(nsh + 78);
    const auto *nsh_80 = buffer.data(nsh + 80);
    const auto *nsh_81 = buffer.data(nsh + 81);
    const auto *nsh_82 = buffer.data(nsh + 82);
    const auto *nsh_83 = buffer.data(nsh + 83);
    const auto *nsh_99 = buffer.data(nsh + 99);

    const auto *nsi1_0 = buffer.data(nsi1 + 0);
    const auto *nsi1_3 = buffer.data(nsi1 + 3);
    const auto *nsi1_5 = buffer.data(nsi1 + 5);
    const auto *nsi1_6 = buffer.data(nsi1 + 6);
    const auto *nsi1_9 = buffer.data(nsi1 + 9);
    const auto *nsi1_10 = buffer.data(nsi1 + 10);
    const auto *nsi1_14 = buffer.data(nsi1 + 14);
    const auto *nsi1_21 = buffer.data(nsi1 + 21);
    const auto *nsi1_27 = buffer.data(nsi1 + 27);
    const auto *nsi1_31 = buffer.data(nsi1 + 31);
    const auto *nsi1_34 = buffer.data(nsi1 + 34);
    const auto *nsi1_38 = buffer.data(nsi1 + 38);
    const auto *nsi1_56 = buffer.data(nsi1 + 56);
    const auto *nsi1_61 = buffer.data(nsi1 + 61);
    const auto *nsi1_65 = buffer.data(nsi1 + 65);
    const auto *nsi1_68 = buffer.data(nsi1 + 68);
    const auto *nsi1_70 = buffer.data(nsi1 + 70);

    const auto *osg0_0 = buffer.data(osg0 + 0);
    const auto *osg0_1 = buffer.data(osg0 + 1);
    const auto *osg0_2 = buffer.data(osg0 + 2);
    const auto *osg0_3 = buffer.data(osg0 + 3);
    const auto *osg0_5 = buffer.data(osg0 + 5);
    const auto *osg0_10 = buffer.data(osg0 + 10);
    const auto *osg0_12 = buffer.data(osg0 + 12);
    const auto *osg0_13 = buffer.data(osg0 + 13);
    const auto *osg0_14 = buffer.data(osg0 + 14);
    const auto *osg0_18 = buffer.data(osg0 + 18);
    const auto *osg0_25 = buffer.data(osg0 + 25);
    const auto *osg0_26 = buffer.data(osg0 + 26);
    const auto *osg0_27 = buffer.data(osg0 + 27);
    const auto *osg0_32 = buffer.data(osg0 + 32);
    const auto *osg0_34 = buffer.data(osg0 + 34);
    const auto *osg0_35 = buffer.data(osg0 + 35);
    const auto *osg0_41 = buffer.data(osg0 + 41);
    const auto *osg0_42 = buffer.data(osg0 + 42);
    const auto *osg0_43 = buffer.data(osg0 + 43);
    const auto *osg0_44 = buffer.data(osg0 + 44);
    const auto *osg0_45 = buffer.data(osg0 + 45);
    const auto *osg0_47 = buffer.data(osg0 + 47);
    const auto *osg0_48 = buffer.data(osg0 + 48);
    const auto *osg0_50 = buffer.data(osg0 + 50);
    const auto *osg0_51 = buffer.data(osg0 + 51);
    const auto *osg0_55 = buffer.data(osg0 + 55);
    const auto *osg0_56 = buffer.data(osg0 + 56);
    const auto *osg0_57 = buffer.data(osg0 + 57);
    const auto *osg0_59 = buffer.data(osg0 + 59);

    const auto *osg1_0 = buffer.data(osg1 + 0);
    const auto *osg1_1 = buffer.data(osg1 + 1);
    const auto *osg1_2 = buffer.data(osg1 + 2);
    const auto *osg1_3 = buffer.data(osg1 + 3);
    const auto *osg1_5 = buffer.data(osg1 + 5);
    const auto *osg1_10 = buffer.data(osg1 + 10);
    const auto *osg1_12 = buffer.data(osg1 + 12);
    const auto *osg1_13 = buffer.data(osg1 + 13);
    const auto *osg1_14 = buffer.data(osg1 + 14);
    const auto *osg1_18 = buffer.data(osg1 + 18);
    const auto *osg1_25 = buffer.data(osg1 + 25);
    const auto *osg1_26 = buffer.data(osg1 + 26);
    const auto *osg1_27 = buffer.data(osg1 + 27);
    const auto *osg1_32 = buffer.data(osg1 + 32);
    const auto *osg1_34 = buffer.data(osg1 + 34);
    const auto *osg1_35 = buffer.data(osg1 + 35);
    const auto *osg1_41 = buffer.data(osg1 + 41);
    const auto *osg1_42 = buffer.data(osg1 + 42);
    const auto *osg1_43 = buffer.data(osg1 + 43);
    const auto *osg1_44 = buffer.data(osg1 + 44);
    const auto *osg1_45 = buffer.data(osg1 + 45);
    const auto *osg1_47 = buffer.data(osg1 + 47);
    const auto *osg1_48 = buffer.data(osg1 + 48);
    const auto *osg1_50 = buffer.data(osg1 + 50);
    const auto *osg1_51 = buffer.data(osg1 + 51);
    const auto *osg1_55 = buffer.data(osg1 + 55);
    const auto *osg1_56 = buffer.data(osg1 + 56);
    const auto *osg1_57 = buffer.data(osg1 + 57);
    const auto *osg1_59 = buffer.data(osg1 + 59);

    const auto *osh_0 = buffer.data(osh + 0);
    const auto *osh_1 = buffer.data(osh + 1);
    const auto *osh_2 = buffer.data(osh + 2);
    const auto *osh_3 = buffer.data(osh + 3);
    const auto *osh_5 = buffer.data(osh + 5);
    const auto *osh_6 = buffer.data(osh + 6);
    const auto *osh_8 = buffer.data(osh + 8);
    const auto *osh_9 = buffer.data(osh + 9);
    const auto *osh_10 = buffer.data(osh + 10);
    const auto *osh_14 = buffer.data(osh + 14);
    const auto *osh_15 = buffer.data(osh + 15);
    const auto *osh_17 = buffer.data(osh + 17);
    const auto *osh_18 = buffer.data(osh + 18);
    const auto *osh_19 = buffer.data(osh + 19);
    const auto *osh_20 = buffer.data(osh + 20);
    const auto *osh_21 = buffer.data(osh + 21);
    const auto *osh_22 = buffer.data(osh + 22);
    const auto *osh_24 = buffer.data(osh + 24);
    const auto *osh_26 = buffer.data(osh + 26);
    const auto *osh_27 = buffer.data(osh + 27);
    const auto *osh_28 = buffer.data(osh + 28);
    const auto *osh_30 = buffer.data(osh + 30);
    const auto *osh_31 = buffer.data(osh + 31);
    const auto *osh_36 = buffer.data(osh + 36);
    const auto *osh_37 = buffer.data(osh + 37);
    const auto *osh_38 = buffer.data(osh + 38);
    const auto *osh_39 = buffer.data(osh + 39);
    const auto *osh_40 = buffer.data(osh + 40);
    const auto *osh_41 = buffer.data(osh + 41);
    const auto *osh_42 = buffer.data(osh + 42);
    const auto *osh_44 = buffer.data(osh + 44);
    const auto *osh_46 = buffer.data(osh + 46);
    const auto *osh_47 = buffer.data(osh + 47);
    const auto *osh_49 = buffer.data(osh + 49);
    const auto *osh_50 = buffer.data(osh + 50);
    const auto *osh_51 = buffer.data(osh + 51);
    const auto *osh_56 = buffer.data(osh + 56);
    const auto *osh_57 = buffer.data(osh + 57);
    const auto *osh_58 = buffer.data(osh + 58);
    const auto *osh_59 = buffer.data(osh + 59);
    const auto *osh_60 = buffer.data(osh + 60);
    const auto *osh_61 = buffer.data(osh + 61);
    const auto *osh_62 = buffer.data(osh + 62);
    const auto *osh_63 = buffer.data(osh + 63);
    const auto *osh_64 = buffer.data(osh + 64);
    const auto *osh_65 = buffer.data(osh + 65);
    const auto *osh_66 = buffer.data(osh + 66);
    const auto *osh_68 = buffer.data(osh + 68);
    const auto *osh_69 = buffer.data(osh + 69);
    const auto *osh_70 = buffer.data(osh + 70);
    const auto *osh_72 = buffer.data(osh + 72);
    const auto *osh_73 = buffer.data(osh + 73);
    const auto *osh_78 = buffer.data(osh + 78);
    const auto *osh_79 = buffer.data(osh + 79);
    const auto *osh_80 = buffer.data(osh + 80);
    const auto *osh_81 = buffer.data(osh + 81);
    const auto *osh_82 = buffer.data(osh + 82);
    const auto *osh_83 = buffer.data(osh + 83);
    const auto *osh_84 = buffer.data(osh + 84);
    const auto *osh_86 = buffer.data(osh + 86);
    const auto *osh_87 = buffer.data(osh + 87);
    const auto *osh_89 = buffer.data(osh + 89);
    const auto *osh_90 = buffer.data(osh + 90);
    const auto *osh_93 = buffer.data(osh + 93);
    const auto *osh_99 = buffer.data(osh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, nsh_0, osg0_0, \
                         osg1_0, osh_0, osh_1, osh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * nsh_0[k]
                 + f_1 * osg0_0[k]
                 - f_2 * osg1_0[k]
                 + f_3 * pc_x[k] * osh_0[k];

        t_1[k] = f_3 * pc_y[k] * osh_0[k];

        t_2[k] = f_3 * pc_z[k] * osh_0[k];

        t_3[k] = f_4 * osg0_0[k]
                 - f_5 * osg1_0[k]
                 + f_3 * pc_y[k] * osh_1[k];

        t_4[k] = f_3 * pc_y[k] * osh_2[k];

        t_5[k] = f_4 * osg0_0[k]
                 - f_5 * osg1_0[k]
                 + f_3 * pc_z[k] * osh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, osg0_1, osg0_2, osg0_3, osg1_1, \
                         osg1_2, osg1_3, osh_3, osh_5, osh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * osg0_1[k]
                 - f_7 * osg1_1[k]
                 + f_3 * pc_y[k] * osh_3[k];

        t_7[k] = f_3 * pc_z[k] * osh_3[k];

        t_8[k] = f_3 * pc_y[k] * osh_5[k];

        t_9[k] = f_6 * osg0_2[k]
                 - f_7 * osg1_2[k]
                 + f_3 * pc_z[k] * osh_5[k];

        t_10[k] = f_8 * osg0_3[k]
                  - f_9 * osg1_3[k]
                  + f_3 * pc_y[k] * osh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, nsh_15, osg0_5, \
                         osg1_5, osh_6, osh_8, osh_9, osh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * osh_6[k];

        t_12[k] = f_4 * osg0_5[k]
                  - f_5 * osg1_5[k]
                  + f_3 * pc_y[k] * osh_8[k];

        t_13[k] = f_3 * pc_y[k] * osh_9[k];

        t_14[k] = f_8 * osg0_5[k]
                  - f_9 * osg1_5[k]
                  + f_3 * pc_z[k] * osh_9[k];

        t_15[k] = f_0 * nsh_15[k]
                  + f_3 * pc_x[k] * osh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, nsh_17, nsh_18, \
                         nsh_20, osh_10, osh_14, osh_17, osh_18, \
                         osh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * osh_10[k];

        t_17[k] = f_0 * nsh_17[k]
                  + f_3 * pc_x[k] * osh_17[k];

        t_18[k] = f_0 * nsh_18[k]
                  + f_3 * pc_x[k] * osh_18[k];

        t_19[k] = f_3 * pc_y[k] * osh_14[k];

        t_20[k] = f_0 * nsh_20[k]
                  + f_3 * pc_x[k] * osh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, osg0_10, osg0_12, osg0_13, \
                         osg1_10, osg1_12, osg1_13, osh_15, osh_17, \
                         osh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * osg0_10[k]
                  - f_2 * osg1_10[k]
                  + f_3 * pc_y[k] * osh_15[k];

        t_22[k] = f_3 * pc_z[k] * osh_15[k];

        t_23[k] = f_8 * osg0_12[k]
                  - f_9 * osg1_12[k]
                  + f_3 * pc_y[k] * osh_17[k];

        t_24[k] = f_6 * osg0_13[k]
                  - f_7 * osg1_13[k]
                  + f_3 * pc_y[k] * osh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, nsi0_0, nsh_0, \
                         nsi1_0, osg0_14, osg1_14, osh_19, osh_20, \
                         osh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * osg0_14[k]
                  - f_5 * osg1_14[k]
                  + f_3 * pc_y[k] * osh_19[k];

        t_26[k] = f_3 * pc_y[k] * osh_20[k];

        t_27[k] = f_1 * osg0_14[k]
                  - f_2 * osg1_14[k]
                  + f_3 * pc_z[k] * osh_20[k];

        t_28[k] = pa_y[k] * nsi0_0[k]
                  - f_10 * pc_y[k] * nsi1_0[k];

        t_29[k] = f_11 * nsh_0[k]
                  + f_3 * pc_y[k] * osh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, nsi0_3, nsi0_5, nsh_1, \
                         nsi1_3, nsi1_5, osh_21, osh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * osh_21[k];

        t_31[k] = pa_y[k] * nsi0_3[k]
                  + f_12 * nsh_1[k]
                  - f_10 * pc_y[k] * nsi1_3[k];

        t_32[k] = f_3 * pc_z[k] * osh_22[k];

        t_33[k] = pa_y[k] * nsi0_5[k]
                  - f_10 * pc_y[k] * nsi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, nsi0_6, nsi0_9, nsh_3, \
                         nsh_5, nsi1_6, nsi1_9, osh_24, osh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * nsi0_6[k]
                  + f_13 * nsh_3[k]
                  - f_10 * pc_y[k] * nsi1_6[k];

        t_35[k] = f_3 * pc_z[k] * osh_24[k];

        t_36[k] = f_11 * nsh_5[k]
                  + f_3 * pc_y[k] * osh_26[k];

        t_37[k] = pa_y[k] * nsi0_9[k]
                  - f_10 * pc_y[k] * nsi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, nsi0_10, nsh_6, nsh_9, \
                         nsi1_10, osg0_18, osg1_18, osh_27, osh_28, \
                         osh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * nsi0_10[k]
                  + f_14 * nsh_6[k]
                  - f_10 * pc_y[k] * nsi1_10[k];

        t_39[k] = f_3 * pc_z[k] * osh_27[k];

        t_40[k] = f_4 * osg0_18[k]
                  - f_5 * osg1_18[k]
                  + f_3 * pc_z[k] * osh_28[k];

        t_41[k] = f_11 * nsh_9[k]
                  + f_3 * pc_y[k] * osh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, nsi0_14, nsh_36, \
                         nsh_38, nsi1_14, osh_31, osh_36, osh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * nsi0_14[k]
                  - f_10 * pc_y[k] * nsi1_14[k];

        t_43[k] = f_15 * nsh_36[k]
                  + f_3 * pc_x[k] * osh_36[k];

        t_44[k] = f_3 * pc_z[k] * osh_31[k];

        t_45[k] = f_15 * nsh_38[k]
                  + f_3 * pc_x[k] * osh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, nsh_15, nsh_39, nsh_40, nsh_41, \
                         osg0_25, osg1_25, osh_36, osh_39, osh_40, \
                         osh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * nsh_39[k]
                  + f_3 * pc_x[k] * osh_39[k];

        t_47[k] = f_15 * nsh_40[k]
                  + f_3 * pc_x[k] * osh_40[k];

        t_48[k] = f_15 * nsh_41[k]
                  + f_3 * pc_x[k] * osh_41[k];

        t_49[k] = f_11 * nsh_15[k]
                  + f_1 * osg0_25[k]
                  - f_2 * osg1_25[k]
                  + f_3 * pc_y[k] * osh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, osg0_25, osg0_26, osg0_27, osg1_25, \
                         osg1_26, osg1_27, osh_36, osh_37, osh_38, \
                         osh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * osh_36[k];

        t_51[k] = f_4 * osg0_25[k]
                  - f_5 * osg1_25[k]
                  + f_3 * pc_z[k] * osh_37[k];

        t_52[k] = f_6 * osg0_26[k]
                  - f_7 * osg1_26[k]
                  + f_3 * pc_z[k] * osh_38[k];

        t_53[k] = f_8 * osg0_27[k]
                  - f_9 * osg1_27[k]
                  + f_3 * pc_z[k] * osh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, nsi0_0, nsi0_27, \
                         nsh_20, nsi1_0, nsi1_27, osh_41, osh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * nsh_20[k]
                  + f_3 * pc_y[k] * osh_41[k];

        t_55[k] = pa_y[k] * nsi0_27[k]
                  - f_10 * pc_y[k] * nsi1_27[k];

        t_56[k] = pa_z[k] * nsi0_0[k]
                  - f_10 * pc_z[k] * nsi1_0[k];

        t_57[k] = f_3 * pc_y[k] * osh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, nsi0_3, nsi0_5, nsh_0, \
                         nsh_2, nsi1_3, nsi1_5, osh_42, osh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * nsh_0[k]
                  + f_3 * pc_z[k] * osh_42[k];

        t_59[k] = pa_z[k] * nsi0_3[k]
                  - f_10 * pc_z[k] * nsi1_3[k];

        t_60[k] = f_3 * pc_y[k] * osh_44[k];

        t_61[k] = pa_z[k] * nsi0_5[k]
                  + f_12 * nsh_2[k]
                  - f_10 * pc_z[k] * nsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, nsi0_6, nsi0_9, nsh_5, \
                         nsi1_6, nsi1_9, osg0_32, osg1_32, osh_46, \
                         osh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * nsi0_6[k]
                  - f_10 * pc_z[k] * nsi1_6[k];

        t_63[k] = f_4 * osg0_32[k]
                  - f_5 * osg1_32[k]
                  + f_3 * pc_y[k] * osh_46[k];

        t_64[k] = f_3 * pc_y[k] * osh_47[k];

        t_65[k] = pa_z[k] * nsi0_9[k]
                  + f_13 * nsh_5[k]
                  - f_10 * pc_z[k] * nsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, nsi0_10, nsi1_10, osg0_34, \
                         osg0_35, osg1_34, osg1_35, osh_49, osh_50, \
                         osh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * nsi0_10[k]
                  - f_10 * pc_z[k] * nsi1_10[k];

        t_67[k] = f_6 * osg0_34[k]
                  - f_7 * osg1_34[k]
                  + f_3 * pc_y[k] * osh_49[k];

        t_68[k] = f_4 * osg0_35[k]
                  - f_5 * osg1_35[k]
                  + f_3 * pc_y[k] * osh_50[k];

        t_69[k] = f_3 * pc_y[k] * osh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, nsi0_14, nsh_9, nsh_57, \
                         nsh_58, nsh_59, nsi1_14, osh_57, osh_58, \
                         osh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * nsi0_14[k]
                  + f_14 * nsh_9[k]
                  - f_10 * pc_z[k] * nsi1_14[k];

        t_71[k] = f_15 * nsh_57[k]
                  + f_3 * pc_x[k] * osh_57[k];

        t_72[k] = f_15 * nsh_58[k]
                  + f_3 * pc_x[k] * osh_58[k];

        t_73[k] = f_15 * nsh_59[k]
                  + f_3 * pc_x[k] * osh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, nsi0_21, nsh_60, \
                         nsh_62, nsi1_21, osh_56, osh_60, osh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * nsh_60[k]
                  + f_3 * pc_x[k] * osh_60[k];

        t_75[k] = f_3 * pc_y[k] * osh_56[k];

        t_76[k] = f_15 * nsh_62[k]
                  + f_3 * pc_x[k] * osh_62[k];

        t_77[k] = pa_z[k] * nsi0_21[k]
                  - f_10 * pc_z[k] * nsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, osg0_41, osg0_42, osg0_43, osg1_41, osg1_42, \
                         osg1_43, osh_58, osh_59, osh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * osg0_41[k]
                  - f_17 * osg1_41[k]
                  + f_3 * pc_y[k] * osh_58[k];

        t_79[k] = f_8 * osg0_42[k]
                  - f_9 * osg1_42[k]
                  + f_3 * pc_y[k] * osh_59[k];

        t_80[k] = f_6 * osg0_43[k]
                  - f_7 * osg1_43[k]
                  + f_3 * pc_y[k] * osh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, nsh_20, nsh_63, osg0_44, \
                         osg0_45, osg1_44, osg1_45, osh_61, osh_62, \
                         osh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * osg0_44[k]
                  - f_5 * osg1_44[k]
                  + f_3 * pc_y[k] * osh_61[k];

        t_82[k] = f_3 * pc_y[k] * osh_62[k];

        t_83[k] = f_11 * nsh_20[k]
                  + f_1 * osg0_44[k]
                  - f_2 * osg1_44[k]
                  + f_3 * pc_z[k] * osh_62[k];

        t_84[k] = f_18 * nsh_63[k]
                  + f_1 * osg0_45[k]
                  - f_2 * osg1_45[k]
                  + f_3 * pc_x[k] * osh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, nsh_21, nsh_66, osg0_48, \
                         osg1_48, osh_63, osh_64, osh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * nsh_21[k]
                  + f_3 * pc_y[k] * osh_63[k];

        t_86[k] = f_3 * pc_z[k] * osh_63[k];

        t_87[k] = f_18 * nsh_66[k]
                  + f_8 * osg0_48[k]
                  - f_9 * osg1_48[k]
                  + f_3 * pc_x[k] * osh_66[k];

        t_88[k] = f_3 * pc_z[k] * osh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, nsh_69, osg0_45, osg0_51, osg1_45, \
                         osg1_51, osh_65, osh_66, osh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * osg0_45[k]
                  - f_5 * osg1_45[k]
                  + f_3 * pc_z[k] * osh_65[k];

        t_90[k] = f_18 * nsh_69[k]
                  + f_6 * osg0_51[k]
                  - f_7 * osg1_51[k]
                  + f_3 * pc_x[k] * osh_69[k];

        t_91[k] = f_3 * pc_z[k] * osh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, nsh_26, nsh_73, osg0_47, \
                         osg0_55, osg1_47, osg1_55, osh_68, osh_69, \
                         osh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * nsh_26[k]
                  + f_3 * pc_y[k] * osh_68[k];

        t_93[k] = f_6 * osg0_47[k]
                  - f_7 * osg1_47[k]
                  + f_3 * pc_z[k] * osh_68[k];

        t_94[k] = f_18 * nsh_73[k]
                  + f_4 * osg0_55[k]
                  - f_5 * osg1_55[k]
                  + f_3 * pc_x[k] * osh_73[k];

        t_95[k] = f_3 * pc_z[k] * osh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, nsh_30, nsh_78, osg0_48, \
                         osg0_50, osg1_48, osg1_50, osh_70, osh_72, \
                         osh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * osg0_48[k]
                  - f_5 * osg1_48[k]
                  + f_3 * pc_z[k] * osh_70[k];

        t_97[k] = f_12 * nsh_30[k]
                  + f_3 * pc_y[k] * osh_72[k];

        t_98[k] = f_8 * osg0_50[k]
                  - f_9 * osg1_50[k]
                  + f_3 * pc_z[k] * osh_72[k];

        t_99[k] = f_18 * nsh_78[k]
                  + f_3 * pc_x[k] * osh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, nsh_80, nsh_81, \
                         nsh_82, nsh_83, osh_73, osh_80, osh_81, osh_82, \
                         osh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * osh_73[k];

        t_101[k] = f_18 * nsh_80[k]
                   + f_3 * pc_x[k] * osh_80[k];

        t_102[k] = f_18 * nsh_81[k]
                   + f_3 * pc_x[k] * osh_81[k];

        t_103[k] = f_18 * nsh_82[k]
                   + f_3 * pc_x[k] * osh_82[k];

        t_104[k] = f_18 * nsh_83[k]
                   + f_3 * pc_x[k] * osh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, nsh_36, osg0_55, osg0_56, \
                         osg1_55, osg1_56, osh_78, osh_79, osh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * nsh_36[k]
                   + f_1 * osg0_55[k]
                   - f_2 * osg1_55[k]
                   + f_3 * pc_y[k] * osh_78[k];

        t_106[k] = f_3 * pc_z[k] * osh_78[k];

        t_107[k] = f_4 * osg0_55[k]
                   - f_5 * osg1_55[k]
                   + f_3 * pc_z[k] * osh_79[k];

        t_108[k] = f_6 * osg0_56[k]
                   - f_7 * osg1_56[k]
                   + f_3 * pc_z[k] * osh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, nsi0_56, nsh_41, \
                         nsi1_56, osg0_57, osg0_59, osg1_57, osg1_59, osh_81, \
                         osh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * osg0_57[k]
                   - f_9 * osg1_57[k]
                   + f_3 * pc_z[k] * osh_81[k];

        t_110[k] = f_12 * nsh_41[k]
                   + f_3 * pc_y[k] * osh_83[k];

        t_111[k] = f_1 * osg0_59[k]
                   - f_2 * osg1_59[k]
                   + f_3 * pc_z[k] * osh_83[k];

        t_112[k] = pa_y[k] * nsi0_56[k]
                   - f_10 * pc_y[k] * nsi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, nsi0_31, nsh_21, \
                         nsh_42, nsh_44, nsi1_31, osh_84, osh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * nsh_42[k]
                   + f_3 * pc_y[k] * osh_84[k];

        t_114[k] = f_11 * nsh_21[k]
                   + f_3 * pc_z[k] * osh_84[k];

        t_115[k] = pa_z[k] * nsi0_31[k]
                   - f_10 * pc_z[k] * nsi1_31[k];

        t_116[k] = f_11 * nsh_44[k]
                   + f_3 * pc_y[k] * osh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, nsi0_34, nsi0_61, \
                         nsh_24, nsh_47, nsi1_34, nsi1_61, osh_87, \
                         osh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * nsi0_61[k]
                   - f_10 * pc_y[k] * nsi1_61[k];

        t_118[k] = pa_z[k] * nsi0_34[k]
                   - f_10 * pc_z[k] * nsi1_34[k];

        t_119[k] = f_11 * nsh_24[k]
                   + f_3 * pc_z[k] * osh_87[k];

        t_120[k] = f_11 * nsh_47[k]
                   + f_3 * pc_y[k] * osh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, nsi0_38, nsi0_65, \
                         nsh_27, nsi1_38, nsi1_65, osh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * nsi0_65[k]
                   - f_10 * pc_y[k] * nsi1_65[k];

        t_122[k] = pa_z[k] * nsi0_38[k]
                   - f_10 * pc_z[k] * nsi1_38[k];

        t_123[k] = f_11 * nsh_27[k]
                   + f_3 * pc_z[k] * osh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, nsi0_68, nsi0_70, \
                         nsh_50, nsh_51, nsh_99, nsi1_68, nsi1_70, osh_93, \
                         osh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * nsi0_68[k]
                   + f_12 * nsh_50[k]
                   - f_10 * pc_y[k] * nsi1_68[k];

        t_125[k] = f_11 * nsh_51[k]
                   + f_3 * pc_y[k] * osh_93[k];

        t_126[k] = pa_y[k] * nsi0_70[k]
                   - f_10 * pc_y[k] * nsi1_70[k];

        t_127[k] = f_18 * nsh_99[k]
                   + f_3 * pc_x[k] * osh_99[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_49 = buffer.data(nsi0 + 49);
    const auto *nsi0_83 = buffer.data(nsi0 + 83);
    const auto *nsi0_84 = buffer.data(nsi0 + 84);
    const auto *nsi0_87 = buffer.data(nsi0 + 87);
    const auto *nsi0_90 = buffer.data(nsi0 + 90);
    const auto *nsi0_94 = buffer.data(nsi0 + 94);
    const auto *nsi0_96 = buffer.data(nsi0 + 96);
    const auto *nsi0_105 = buffer.data(nsi0 + 105);
    const auto *nsi0_140 = buffer.data(nsi0 + 140);
    const auto *nsi0_143 = buffer.data(nsi0 + 143);
    const auto *nsi0_145 = buffer.data(nsi0 + 145);
    const auto *nsi0_146 = buffer.data(nsi0 + 146);
    const auto *nsi0_149 = buffer.data(nsi0 + 149);
    const auto *nsi0_150 = buffer.data(nsi0 + 150);
    const auto *nsi0_152 = buffer.data(nsi0 + 152);
    const auto *nsi0_154 = buffer.data(nsi0 + 154);

    const auto *nsh_36 = buffer.data(nsh + 36);
    const auto *nsh_42 = buffer.data(nsh + 42);
    const auto *nsh_59 = buffer.data(nsh + 59);
    const auto *nsh_60 = buffer.data(nsh + 60);
    const auto *nsh_61 = buffer.data(nsh + 61);
    const auto *nsh_62 = buffer.data(nsh + 62);
    const auto *nsh_63 = buffer.data(nsh + 63);
    const auto *nsh_66 = buffer.data(nsh + 66);
    const auto *nsh_68 = buffer.data(nsh + 68);
    const auto *nsh_69 = buffer.data(nsh + 69);
    const auto *nsh_70 = buffer.data(nsh + 70);
    const auto *nsh_72 = buffer.data(nsh + 72);
    const auto *nsh_78 = buffer.data(nsh + 78);
    const auto *nsh_83 = buffer.data(nsh + 83);
    const auto *nsh_84 = buffer.data(nsh + 84);
    const auto *nsh_86 = buffer.data(nsh + 86);
    const auto *nsh_87 = buffer.data(nsh + 87);
    const auto *nsh_89 = buffer.data(nsh + 89);
    const auto *nsh_90 = buffer.data(nsh + 90);
    const auto *nsh_93 = buffer.data(nsh + 93);
    const auto *nsh_99 = buffer.data(nsh + 99);
    const auto *nsh_100 = buffer.data(nsh + 100);
    const auto *nsh_101 = buffer.data(nsh + 101);
    const auto *nsh_102 = buffer.data(nsh + 102);
    const auto *nsh_103 = buffer.data(nsh + 103);
    const auto *nsh_104 = buffer.data(nsh + 104);
    const auto *nsh_105 = buffer.data(nsh + 105);
    const auto *nsh_106 = buffer.data(nsh + 106);
    const auto *nsh_107 = buffer.data(nsh + 107);
    const auto *nsh_108 = buffer.data(nsh + 108);
    const auto *nsh_110 = buffer.data(nsh + 110);
    const auto *nsh_111 = buffer.data(nsh + 111);
    const auto *nsh_113 = buffer.data(nsh + 113);
    const auto *nsh_114 = buffer.data(nsh + 114);
    const auto *nsh_119 = buffer.data(nsh + 119);
    const auto *nsh_120 = buffer.data(nsh + 120);
    const auto *nsh_121 = buffer.data(nsh + 121);
    const auto *nsh_122 = buffer.data(nsh + 122);
    const auto *nsh_123 = buffer.data(nsh + 123);
    const auto *nsh_125 = buffer.data(nsh + 125);
    const auto *nsh_126 = buffer.data(nsh + 126);
    const auto *nsh_129 = buffer.data(nsh + 129);
    const auto *nsh_132 = buffer.data(nsh + 132);
    const auto *nsh_136 = buffer.data(nsh + 136);
    const auto *nsh_141 = buffer.data(nsh + 141);
    const auto *nsh_143 = buffer.data(nsh + 143);
    const auto *nsh_144 = buffer.data(nsh + 144);
    const auto *nsh_145 = buffer.data(nsh + 145);
    const auto *nsh_146 = buffer.data(nsh + 146);
    const auto *nsh_152 = buffer.data(nsh + 152);
    const auto *nsh_156 = buffer.data(nsh + 156);
    const auto *nsh_161 = buffer.data(nsh + 161);
    const auto *nsh_162 = buffer.data(nsh + 162);
    const auto *nsh_163 = buffer.data(nsh + 163);
    const auto *nsh_164 = buffer.data(nsh + 164);
    const auto *nsh_165 = buffer.data(nsh + 165);
    const auto *nsh_166 = buffer.data(nsh + 166);
    const auto *nsh_167 = buffer.data(nsh + 167);
    const auto *nsh_183 = buffer.data(nsh + 183);
    const auto *nsh_184 = buffer.data(nsh + 184);
    const auto *nsh_185 = buffer.data(nsh + 185);
    const auto *nsh_186 = buffer.data(nsh + 186);
    const auto *nsh_187 = buffer.data(nsh + 187);
    const auto *nsh_188 = buffer.data(nsh + 188);

    const auto *nsi1_49 = buffer.data(nsi1 + 49);
    const auto *nsi1_83 = buffer.data(nsi1 + 83);
    const auto *nsi1_84 = buffer.data(nsi1 + 84);
    const auto *nsi1_87 = buffer.data(nsi1 + 87);
    const auto *nsi1_90 = buffer.data(nsi1 + 90);
    const auto *nsi1_94 = buffer.data(nsi1 + 94);
    const auto *nsi1_96 = buffer.data(nsi1 + 96);
    const auto *nsi1_105 = buffer.data(nsi1 + 105);
    const auto *nsi1_140 = buffer.data(nsi1 + 140);
    const auto *nsi1_143 = buffer.data(nsi1 + 143);
    const auto *nsi1_145 = buffer.data(nsi1 + 145);
    const auto *nsi1_146 = buffer.data(nsi1 + 146);
    const auto *nsi1_149 = buffer.data(nsi1 + 149);
    const auto *nsi1_150 = buffer.data(nsi1 + 150);
    const auto *nsi1_152 = buffer.data(nsi1 + 152);
    const auto *nsi1_154 = buffer.data(nsi1 + 154);

    const auto *osg0_72 = buffer.data(osg0 + 72);
    const auto *osg0_73 = buffer.data(osg0 + 73);
    const auto *osg0_74 = buffer.data(osg0 + 74);
    const auto *osg0_75 = buffer.data(osg0 + 75);
    const auto *osg0_76 = buffer.data(osg0 + 76);
    const auto *osg0_77 = buffer.data(osg0 + 77);
    const auto *osg0_78 = buffer.data(osg0 + 78);
    const auto *osg0_79 = buffer.data(osg0 + 79);
    const auto *osg0_80 = buffer.data(osg0 + 80);
    const auto *osg0_84 = buffer.data(osg0 + 84);
    const auto *osg0_85 = buffer.data(osg0 + 85);
    const auto *osg0_86 = buffer.data(osg0 + 86);
    const auto *osg0_87 = buffer.data(osg0 + 87);
    const auto *osg0_88 = buffer.data(osg0 + 88);
    const auto *osg0_89 = buffer.data(osg0 + 89);
    const auto *osg0_90 = buffer.data(osg0 + 90);
    const auto *osg0_92 = buffer.data(osg0 + 92);
    const auto *osg0_93 = buffer.data(osg0 + 93);
    const auto *osg0_95 = buffer.data(osg0 + 95);
    const auto *osg0_96 = buffer.data(osg0 + 96);
    const auto *osg0_100 = buffer.data(osg0 + 100);
    const auto *osg0_101 = buffer.data(osg0 + 101);
    const auto *osg0_102 = buffer.data(osg0 + 102);
    const auto *osg0_104 = buffer.data(osg0 + 104);
    const auto *osg0_110 = buffer.data(osg0 + 110);
    const auto *osg0_114 = buffer.data(osg0 + 114);
    const auto *osg0_117 = buffer.data(osg0 + 117);
    const auto *osg0_118 = buffer.data(osg0 + 118);
    const auto *osg0_119 = buffer.data(osg0 + 119);
    const auto *osg0_130 = buffer.data(osg0 + 130);
    const auto *osg0_132 = buffer.data(osg0 + 132);

    const auto *osg1_72 = buffer.data(osg1 + 72);
    const auto *osg1_73 = buffer.data(osg1 + 73);
    const auto *osg1_74 = buffer.data(osg1 + 74);
    const auto *osg1_75 = buffer.data(osg1 + 75);
    const auto *osg1_76 = buffer.data(osg1 + 76);
    const auto *osg1_77 = buffer.data(osg1 + 77);
    const auto *osg1_78 = buffer.data(osg1 + 78);
    const auto *osg1_79 = buffer.data(osg1 + 79);
    const auto *osg1_80 = buffer.data(osg1 + 80);
    const auto *osg1_84 = buffer.data(osg1 + 84);
    const auto *osg1_85 = buffer.data(osg1 + 85);
    const auto *osg1_86 = buffer.data(osg1 + 86);
    const auto *osg1_87 = buffer.data(osg1 + 87);
    const auto *osg1_88 = buffer.data(osg1 + 88);
    const auto *osg1_89 = buffer.data(osg1 + 89);
    const auto *osg1_90 = buffer.data(osg1 + 90);
    const auto *osg1_92 = buffer.data(osg1 + 92);
    const auto *osg1_93 = buffer.data(osg1 + 93);
    const auto *osg1_95 = buffer.data(osg1 + 95);
    const auto *osg1_96 = buffer.data(osg1 + 96);
    const auto *osg1_100 = buffer.data(osg1 + 100);
    const auto *osg1_101 = buffer.data(osg1 + 101);
    const auto *osg1_102 = buffer.data(osg1 + 102);
    const auto *osg1_104 = buffer.data(osg1 + 104);
    const auto *osg1_110 = buffer.data(osg1 + 110);
    const auto *osg1_114 = buffer.data(osg1 + 114);
    const auto *osg1_117 = buffer.data(osg1 + 117);
    const auto *osg1_118 = buffer.data(osg1 + 118);
    const auto *osg1_119 = buffer.data(osg1 + 119);
    const auto *osg1_130 = buffer.data(osg1 + 130);
    const auto *osg1_132 = buffer.data(osg1 + 132);

    const auto *osh_99 = buffer.data(osh + 99);
    const auto *osh_100 = buffer.data(osh + 100);
    const auto *osh_101 = buffer.data(osh + 101);
    const auto *osh_102 = buffer.data(osh + 102);
    const auto *osh_103 = buffer.data(osh + 103);
    const auto *osh_104 = buffer.data(osh + 104);
    const auto *osh_105 = buffer.data(osh + 105);
    const auto *osh_106 = buffer.data(osh + 106);
    const auto *osh_107 = buffer.data(osh + 107);
    const auto *osh_108 = buffer.data(osh + 108);
    const auto *osh_109 = buffer.data(osh + 109);
    const auto *osh_110 = buffer.data(osh + 110);
    const auto *osh_111 = buffer.data(osh + 111);
    const auto *osh_112 = buffer.data(osh + 112);
    const auto *osh_113 = buffer.data(osh + 113);
    const auto *osh_114 = buffer.data(osh + 114);
    const auto *osh_119 = buffer.data(osh + 119);
    const auto *osh_120 = buffer.data(osh + 120);
    const auto *osh_121 = buffer.data(osh + 121);
    const auto *osh_122 = buffer.data(osh + 122);
    const auto *osh_123 = buffer.data(osh + 123);
    const auto *osh_124 = buffer.data(osh + 124);
    const auto *osh_125 = buffer.data(osh + 125);
    const auto *osh_126 = buffer.data(osh + 126);
    const auto *osh_127 = buffer.data(osh + 127);
    const auto *osh_128 = buffer.data(osh + 128);
    const auto *osh_129 = buffer.data(osh + 129);
    const auto *osh_131 = buffer.data(osh + 131);
    const auto *osh_132 = buffer.data(osh + 132);
    const auto *osh_133 = buffer.data(osh + 133);
    const auto *osh_135 = buffer.data(osh + 135);
    const auto *osh_136 = buffer.data(osh + 136);
    const auto *osh_141 = buffer.data(osh + 141);
    const auto *osh_142 = buffer.data(osh + 142);
    const auto *osh_143 = buffer.data(osh + 143);
    const auto *osh_144 = buffer.data(osh + 144);
    const auto *osh_145 = buffer.data(osh + 145);
    const auto *osh_146 = buffer.data(osh + 146);
    const auto *osh_147 = buffer.data(osh + 147);
    const auto *osh_149 = buffer.data(osh + 149);
    const auto *osh_150 = buffer.data(osh + 150);
    const auto *osh_152 = buffer.data(osh + 152);
    const auto *osh_153 = buffer.data(osh + 153);
    const auto *osh_156 = buffer.data(osh + 156);
    const auto *osh_161 = buffer.data(osh + 161);
    const auto *osh_162 = buffer.data(osh + 162);
    const auto *osh_163 = buffer.data(osh + 163);
    const auto *osh_164 = buffer.data(osh + 164);
    const auto *osh_165 = buffer.data(osh + 165);
    const auto *osh_166 = buffer.data(osh + 166);
    const auto *osh_167 = buffer.data(osh + 167);
    const auto *osh_168 = buffer.data(osh + 168);
    const auto *osh_170 = buffer.data(osh + 170);
    const auto *osh_171 = buffer.data(osh + 171);
    const auto *osh_173 = buffer.data(osh + 173);
    const auto *osh_174 = buffer.data(osh + 174);
    const auto *osh_177 = buffer.data(osh + 177);
    const auto *osh_183 = buffer.data(osh + 183);
    const auto *osh_184 = buffer.data(osh + 184);
    const auto *osh_185 = buffer.data(osh + 185);
    const auto *osh_186 = buffer.data(osh + 186);
    const auto *osh_187 = buffer.data(osh + 187);
    const auto *osh_188 = buffer.data(osh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, nsh_100, nsh_101, nsh_102, \
                         nsh_103, nsh_104, osh_100, osh_101, osh_102, osh_103, \
                         osh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_18 * nsh_100[k]
                   + f_3 * pc_x[k] * osh_100[k];

        t_129[k] = f_18 * nsh_101[k]
                   + f_3 * pc_x[k] * osh_101[k];

        t_130[k] = f_18 * nsh_102[k]
                   + f_3 * pc_x[k] * osh_102[k];

        t_131[k] = f_18 * nsh_103[k]
                   + f_3 * pc_x[k] * osh_103[k];

        t_132[k] = f_18 * nsh_104[k]
                   + f_3 * pc_x[k] * osh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, nsi0_49, nsh_36, nsh_59, \
                         nsi1_49, osg0_72, osg1_72, osh_99, osh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * nsi0_49[k]
                   - f_10 * pc_z[k] * nsi1_49[k];

        t_134[k] = f_11 * nsh_36[k]
                   + f_3 * pc_z[k] * osh_99[k];

        t_135[k] = f_11 * nsh_59[k]
                   + f_8 * osg0_72[k]
                   - f_9 * osg1_72[k]
                   + f_3 * pc_y[k] * osh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, nsh_60, nsh_61, nsh_62, osg0_73, osg0_74, \
                         osg1_73, osg1_74, osh_102, osh_103, osh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * nsh_60[k]
                   + f_6 * osg0_73[k]
                   - f_7 * osg1_73[k]
                   + f_3 * pc_y[k] * osh_102[k];

        t_137[k] = f_11 * nsh_61[k]
                   + f_4 * osg0_74[k]
                   - f_5 * osg1_74[k]
                   + f_3 * pc_y[k] * osh_103[k];

        t_138[k] = f_11 * nsh_62[k]
                   + f_3 * pc_y[k] * osh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, nsi0_83, nsh_42, \
                         nsh_105, nsi1_83, osg0_75, osg1_75, osh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * nsi0_83[k]
                   - f_10 * pc_y[k] * nsi1_83[k];

        t_140[k] = f_18 * nsh_105[k]
                   + f_1 * osg0_75[k]
                   - f_2 * osg1_75[k]
                   + f_3 * pc_x[k] * osh_105[k];

        t_141[k] = f_3 * pc_y[k] * osh_105[k];

        t_142[k] = f_12 * nsh_42[k]
                   + f_3 * pc_z[k] * osh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, nsh_110, osg0_75, osg0_80, osg1_75, \
                         osg1_80, osh_106, osh_107, osh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * osg0_75[k]
                   - f_5 * osg1_75[k]
                   + f_3 * pc_y[k] * osh_106[k];

        t_144[k] = f_3 * pc_y[k] * osh_107[k];

        t_145[k] = f_18 * nsh_110[k]
                   + f_8 * osg0_80[k]
                   - f_9 * osg1_80[k]
                   + f_3 * pc_x[k] * osh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, osg0_76, osg0_77, osg1_76, osg1_77, \
                         osh_108, osh_109, osh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * osg0_76[k]
                   - f_7 * osg1_76[k]
                   + f_3 * pc_y[k] * osh_108[k];

        t_147[k] = f_4 * osg0_77[k]
                   - f_5 * osg1_77[k]
                   + f_3 * pc_y[k] * osh_109[k];

        t_148[k] = f_3 * pc_y[k] * osh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, nsh_114, osg0_78, osg0_79, osg0_84, \
                         osg1_78, osg1_79, osg1_84, osh_111, osh_112, \
                         osh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_18 * nsh_114[k]
                   + f_6 * osg0_84[k]
                   - f_7 * osg1_84[k]
                   + f_3 * pc_x[k] * osh_114[k];

        t_150[k] = f_8 * osg0_78[k]
                   - f_9 * osg1_78[k]
                   + f_3 * pc_y[k] * osh_111[k];

        t_151[k] = f_6 * osg0_79[k]
                   - f_7 * osg1_79[k]
                   + f_3 * pc_y[k] * osh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, nsh_119, nsh_120, osg0_80, \
                         osg0_89, osg1_80, osg1_89, osh_113, osh_114, osh_119, \
                         osh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * osg0_80[k]
                   - f_5 * osg1_80[k]
                   + f_3 * pc_y[k] * osh_113[k];

        t_153[k] = f_3 * pc_y[k] * osh_114[k];

        t_154[k] = f_18 * nsh_119[k]
                   + f_4 * osg0_89[k]
                   - f_5 * osg1_89[k]
                   + f_3 * pc_x[k] * osh_119[k];

        t_155[k] = f_18 * nsh_120[k]
                   + f_3 * pc_x[k] * osh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, nsh_121, nsh_122, \
                         nsh_123, nsh_125, osh_119, osh_121, osh_122, osh_123, \
                         osh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * nsh_121[k]
                   + f_3 * pc_x[k] * osh_121[k];

        t_157[k] = f_18 * nsh_122[k]
                   + f_3 * pc_x[k] * osh_122[k];

        t_158[k] = f_18 * nsh_123[k]
                   + f_3 * pc_x[k] * osh_123[k];

        t_159[k] = f_3 * pc_y[k] * osh_119[k];

        t_160[k] = f_18 * nsh_125[k]
                   + f_3 * pc_x[k] * osh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, osg0_85, osg0_86, osg0_87, osg1_85, \
                         osg1_86, osg1_87, osh_120, osh_121, osh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * osg0_85[k]
                   - f_2 * osg1_85[k]
                   + f_3 * pc_y[k] * osh_120[k];

        t_162[k] = f_16 * osg0_86[k]
                   - f_17 * osg1_86[k]
                   + f_3 * pc_y[k] * osh_121[k];

        t_163[k] = f_8 * osg0_87[k]
                   - f_9 * osg1_87[k]
                   + f_3 * pc_y[k] * osh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, nsh_62, osg0_88, osg0_89, \
                         osg1_88, osg1_89, osh_123, osh_124, osh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * osg0_88[k]
                   - f_7 * osg1_88[k]
                   + f_3 * pc_y[k] * osh_123[k];

        t_165[k] = f_4 * osg0_89[k]
                   - f_5 * osg1_89[k]
                   + f_3 * pc_y[k] * osh_124[k];

        t_166[k] = f_3 * pc_y[k] * osh_125[k];

        t_167[k] = f_12 * nsh_62[k]
                   + f_1 * osg0_89[k]
                   - f_2 * osg1_89[k]
                   + f_3 * pc_z[k] * osh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, nsh_63, nsh_126, \
                         nsh_129, osg0_90, osg0_93, osg1_90, osg1_93, osh_126, \
                         osh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * nsh_126[k]
                   + f_1 * osg0_90[k]
                   - f_2 * osg1_90[k]
                   + f_3 * pc_x[k] * osh_126[k];

        t_169[k] = f_13 * nsh_63[k]
                   + f_3 * pc_y[k] * osh_126[k];

        t_170[k] = f_3 * pc_z[k] * osh_126[k];

        t_171[k] = f_19 * nsh_129[k]
                   + f_8 * osg0_93[k]
                   - f_9 * osg1_93[k]
                   + f_3 * pc_x[k] * osh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, nsh_132, osg0_90, osg0_96, \
                         osg1_90, osg1_96, osh_127, osh_128, osh_129, \
                         osh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * osh_127[k];

        t_173[k] = f_4 * osg0_90[k]
                   - f_5 * osg1_90[k]
                   + f_3 * pc_z[k] * osh_128[k];

        t_174[k] = f_19 * nsh_132[k]
                   + f_6 * osg0_96[k]
                   - f_7 * osg1_96[k]
                   + f_3 * pc_x[k] * osh_132[k];

        t_175[k] = f_3 * pc_z[k] * osh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, nsh_68, nsh_136, \
                         osg0_92, osg0_100, osg1_92, osg1_100, osh_131, osh_132, \
                         osh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * nsh_68[k]
                   + f_3 * pc_y[k] * osh_131[k];

        t_177[k] = f_6 * osg0_92[k]
                   - f_7 * osg1_92[k]
                   + f_3 * pc_z[k] * osh_131[k];

        t_178[k] = f_19 * nsh_136[k]
                   + f_4 * osg0_100[k]
                   - f_5 * osg1_100[k]
                   + f_3 * pc_x[k] * osh_136[k];

        t_179[k] = f_3 * pc_z[k] * osh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, nsh_72, nsh_141, \
                         osg0_93, osg0_95, osg1_93, osg1_95, osh_133, osh_135, \
                         osh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * osg0_93[k]
                   - f_5 * osg1_93[k]
                   + f_3 * pc_z[k] * osh_133[k];

        t_181[k] = f_13 * nsh_72[k]
                   + f_3 * pc_y[k] * osh_135[k];

        t_182[k] = f_8 * osg0_95[k]
                   - f_9 * osg1_95[k]
                   + f_3 * pc_z[k] * osh_135[k];

        t_183[k] = f_19 * nsh_141[k]
                   + f_3 * pc_x[k] * osh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, nsh_143, nsh_144, \
                         nsh_145, nsh_146, osh_136, osh_143, osh_144, osh_145, \
                         osh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * osh_136[k];

        t_185[k] = f_19 * nsh_143[k]
                   + f_3 * pc_x[k] * osh_143[k];

        t_186[k] = f_19 * nsh_144[k]
                   + f_3 * pc_x[k] * osh_144[k];

        t_187[k] = f_19 * nsh_145[k]
                   + f_3 * pc_x[k] * osh_145[k];

        t_188[k] = f_19 * nsh_146[k]
                   + f_3 * pc_x[k] * osh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, nsh_78, osg0_100, osg0_101, \
                         osg1_100, osg1_101, osh_141, osh_142, \
                         osh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * nsh_78[k]
                   + f_1 * osg0_100[k]
                   - f_2 * osg1_100[k]
                   + f_3 * pc_y[k] * osh_141[k];

        t_190[k] = f_3 * pc_z[k] * osh_141[k];

        t_191[k] = f_4 * osg0_100[k]
                   - f_5 * osg1_100[k]
                   + f_3 * pc_z[k] * osh_142[k];

        t_192[k] = f_6 * osg0_101[k]
                   - f_7 * osg1_101[k]
                   + f_3 * pc_z[k] * osh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, nsi0_84, nsh_83, \
                         nsi1_84, osg0_102, osg0_104, osg1_102, osg1_104, osh_144, \
                         osh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * osg0_102[k]
                   - f_9 * osg1_102[k]
                   + f_3 * pc_z[k] * osh_144[k];

        t_194[k] = f_13 * nsh_83[k]
                   + f_3 * pc_y[k] * osh_146[k];

        t_195[k] = f_1 * osg0_104[k]
                   - f_2 * osg1_104[k]
                   + f_3 * pc_z[k] * osh_146[k];

        t_196[k] = pa_z[k] * nsi0_84[k]
                   - f_10 * pc_z[k] * nsi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, nsi0_87, nsh_63, \
                         nsh_84, nsh_86, nsi1_87, osh_147, osh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * nsh_84[k]
                   + f_3 * pc_y[k] * osh_147[k];

        t_198[k] = f_11 * nsh_63[k]
                   + f_3 * pc_z[k] * osh_147[k];

        t_199[k] = pa_z[k] * nsi0_87[k]
                   - f_10 * pc_z[k] * nsi1_87[k];

        t_200[k] = f_12 * nsh_86[k]
                   + f_3 * pc_y[k] * osh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, nsi0_90, nsh_66, nsh_152, \
                         nsi1_90, osg0_110, osg1_110, osh_150, \
                         osh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_19 * nsh_152[k]
                   + f_8 * osg0_110[k]
                   - f_9 * osg1_110[k]
                   + f_3 * pc_x[k] * osh_152[k];

        t_202[k] = pa_z[k] * nsi0_90[k]
                   - f_10 * pc_z[k] * nsi1_90[k];

        t_203[k] = f_11 * nsh_66[k]
                   + f_3 * pc_z[k] * osh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, nsi0_94, nsh_89, \
                         nsh_156, nsi1_94, osg0_114, osg1_114, osh_152, \
                         osh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * nsh_89[k]
                   + f_3 * pc_y[k] * osh_152[k];

        t_205[k] = f_19 * nsh_156[k]
                   + f_6 * osg0_114[k]
                   - f_7 * osg1_114[k]
                   + f_3 * pc_x[k] * osh_156[k];

        t_206[k] = pa_z[k] * nsi0_94[k]
                   - f_10 * pc_z[k] * nsi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, nsi0_96, nsh_69, nsh_70, \
                         nsh_93, nsi1_96, osh_153, osh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * nsh_69[k]
                   + f_3 * pc_z[k] * osh_153[k];

        t_208[k] = pa_z[k] * nsi0_96[k]
                   + f_12 * nsh_70[k]
                   - f_10 * pc_z[k] * nsi1_96[k];

        t_209[k] = f_12 * nsh_93[k]
                   + f_3 * pc_y[k] * osh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, nsh_161, nsh_162, nsh_163, nsh_164, \
                         osg0_119, osg1_119, osh_161, osh_162, osh_163, \
                         osh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_19 * nsh_161[k]
                   + f_4 * osg0_119[k]
                   - f_5 * osg1_119[k]
                   + f_3 * pc_x[k] * osh_161[k];

        t_211[k] = f_19 * nsh_162[k]
                   + f_3 * pc_x[k] * osh_162[k];

        t_212[k] = f_19 * nsh_163[k]
                   + f_3 * pc_x[k] * osh_163[k];

        t_213[k] = f_19 * nsh_164[k]
                   + f_3 * pc_x[k] * osh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, nsi0_105, nsh_165, \
                         nsh_166, nsh_167, nsi1_105, osh_165, osh_166, \
                         osh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_19 * nsh_165[k]
                   + f_3 * pc_x[k] * osh_165[k];

        t_215[k] = f_19 * nsh_166[k]
                   + f_3 * pc_x[k] * osh_166[k];

        t_216[k] = f_19 * nsh_167[k]
                   + f_3 * pc_x[k] * osh_167[k];

        t_217[k] = pa_z[k] * nsi0_105[k]
                   - f_10 * pc_z[k] * nsi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, nsh_78, nsh_101, nsh_102, osg0_117, \
                         osg0_118, osg1_117, osg1_118, osh_162, osh_164, \
                         osh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * nsh_78[k]
                   + f_3 * pc_z[k] * osh_162[k];

        t_219[k] = f_12 * nsh_101[k]
                   + f_8 * osg0_117[k]
                   - f_9 * osg1_117[k]
                   + f_3 * pc_y[k] * osh_164[k];

        t_220[k] = f_12 * nsh_102[k]
                   + f_6 * osg0_118[k]
                   - f_7 * osg1_118[k]
                   + f_3 * pc_y[k] * osh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, nsi0_140, nsh_83, \
                         nsh_103, nsh_104, nsi1_140, osg0_119, osg1_119, osh_166, \
                         osh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * nsh_103[k]
                   + f_4 * osg0_119[k]
                   - f_5 * osg1_119[k]
                   + f_3 * pc_y[k] * osh_166[k];

        t_222[k] = f_12 * nsh_104[k]
                   + f_3 * pc_y[k] * osh_167[k];

        t_223[k] = f_11 * nsh_83[k]
                   + f_1 * osg0_119[k]
                   - f_2 * osg1_119[k]
                   + f_3 * pc_z[k] * osh_167[k];

        t_224[k] = pa_y[k] * nsi0_140[k]
                   - f_10 * pc_y[k] * nsi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, nsi0_143, nsh_84, \
                         nsh_105, nsh_106, nsh_107, nsi1_143, osh_168, \
                         osh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * nsh_105[k]
                   + f_3 * pc_y[k] * osh_168[k];

        t_226[k] = f_12 * nsh_84[k]
                   + f_3 * pc_z[k] * osh_168[k];

        t_227[k] = pa_y[k] * nsi0_143[k]
                   + f_12 * nsh_106[k]
                   - f_10 * pc_y[k] * nsi1_143[k];

        t_228[k] = f_11 * nsh_107[k]
                   + f_3 * pc_y[k] * osh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, nsi0_145, nsi0_146, \
                         nsh_87, nsh_108, nsh_110, nsi1_145, nsi1_146, osh_171, \
                         osh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * nsi0_145[k]
                   - f_10 * pc_y[k] * nsi1_145[k];

        t_230[k] = pa_y[k] * nsi0_146[k]
                   + f_13 * nsh_108[k]
                   - f_10 * pc_y[k] * nsi1_146[k];

        t_231[k] = f_12 * nsh_87[k]
                   + f_3 * pc_z[k] * osh_171[k];

        t_232[k] = f_11 * nsh_110[k]
                   + f_3 * pc_y[k] * osh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, nsi0_149, nsi0_150, nsh_90, \
                         nsh_111, nsi1_149, nsi1_150, osh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * nsi0_149[k]
                   - f_10 * pc_y[k] * nsi1_149[k];

        t_234[k] = pa_y[k] * nsi0_150[k]
                   + f_14 * nsh_111[k]
                   - f_10 * pc_y[k] * nsi1_150[k];

        t_235[k] = f_12 * nsh_90[k]
                   + f_3 * pc_z[k] * osh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, nsi0_152, nsi0_154, \
                         nsh_113, nsh_114, nsh_183, nsi1_152, nsi1_154, osh_177, \
                         osh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * nsi0_152[k]
                   + f_12 * nsh_113[k]
                   - f_10 * pc_y[k] * nsi1_152[k];

        t_237[k] = f_11 * nsh_114[k]
                   + f_3 * pc_y[k] * osh_177[k];

        t_238[k] = pa_y[k] * nsi0_154[k]
                   - f_10 * pc_y[k] * nsi1_154[k];

        t_239[k] = f_19 * nsh_183[k]
                   + f_3 * pc_x[k] * osh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, nsh_184, nsh_185, nsh_186, \
                         nsh_187, nsh_188, osh_184, osh_185, osh_186, osh_187, \
                         osh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_19 * nsh_184[k]
                   + f_3 * pc_x[k] * osh_184[k];

        t_241[k] = f_19 * nsh_185[k]
                   + f_3 * pc_x[k] * osh_185[k];

        t_242[k] = f_19 * nsh_186[k]
                   + f_3 * pc_x[k] * osh_186[k];

        t_243[k] = f_19 * nsh_187[k]
                   + f_3 * pc_x[k] * osh_187[k];

        t_244[k] = f_19 * nsh_188[k]
                   + f_3 * pc_x[k] * osh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, nsh_99, nsh_120, nsh_122, osg0_130, \
                         osg0_132, osg1_130, osg1_132, osh_183, \
                         osh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * nsh_120[k]
                   + f_1 * osg0_130[k]
                   - f_2 * osg1_130[k]
                   + f_3 * pc_y[k] * osh_183[k];

        t_246[k] = f_12 * nsh_99[k]
                   + f_3 * pc_z[k] * osh_183[k];

        t_247[k] = f_11 * nsh_122[k]
                   + f_8 * osg0_132[k]
                   - f_9 * osg1_132[k]
                   + f_3 * pc_y[k] * osh_185[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_167 = buffer.data(nsi0 + 167);
    const auto *nsi0_168 = buffer.data(nsi0 + 168);
    const auto *nsi0_171 = buffer.data(nsi0 + 171);
    const auto *nsi0_174 = buffer.data(nsi0 + 174);
    const auto *nsi0_178 = buffer.data(nsi0 + 178);
    const auto *nsi0_180 = buffer.data(nsi0 + 180);
    const auto *nsi0_189 = buffer.data(nsi0 + 189);

    const auto *nsh_105 = buffer.data(nsh + 105);
    const auto *nsh_123 = buffer.data(nsh + 123);
    const auto *nsh_124 = buffer.data(nsh + 124);
    const auto *nsh_125 = buffer.data(nsh + 125);
    const auto *nsh_126 = buffer.data(nsh + 126);
    const auto *nsh_129 = buffer.data(nsh + 129);
    const auto *nsh_131 = buffer.data(nsh + 131);
    const auto *nsh_132 = buffer.data(nsh + 132);
    const auto *nsh_133 = buffer.data(nsh + 133);
    const auto *nsh_135 = buffer.data(nsh + 135);
    const auto *nsh_141 = buffer.data(nsh + 141);
    const auto *nsh_146 = buffer.data(nsh + 146);
    const auto *nsh_147 = buffer.data(nsh + 147);
    const auto *nsh_149 = buffer.data(nsh + 149);
    const auto *nsh_150 = buffer.data(nsh + 150);
    const auto *nsh_152 = buffer.data(nsh + 152);
    const auto *nsh_153 = buffer.data(nsh + 153);
    const auto *nsh_156 = buffer.data(nsh + 156);
    const auto *nsh_162 = buffer.data(nsh + 162);
    const auto *nsh_164 = buffer.data(nsh + 164);
    const auto *nsh_165 = buffer.data(nsh + 165);
    const auto *nsh_166 = buffer.data(nsh + 166);
    const auto *nsh_167 = buffer.data(nsh + 167);
    const auto *nsh_168 = buffer.data(nsh + 168);
    const auto *nsh_170 = buffer.data(nsh + 170);
    const auto *nsh_173 = buffer.data(nsh + 173);
    const auto *nsh_177 = buffer.data(nsh + 177);
    const auto *nsh_183 = buffer.data(nsh + 183);
    const auto *nsh_185 = buffer.data(nsh + 185);
    const auto *nsh_186 = buffer.data(nsh + 186);
    const auto *nsh_187 = buffer.data(nsh + 187);
    const auto *nsh_189 = buffer.data(nsh + 189);
    const auto *nsh_194 = buffer.data(nsh + 194);
    const auto *nsh_198 = buffer.data(nsh + 198);
    const auto *nsh_203 = buffer.data(nsh + 203);
    const auto *nsh_204 = buffer.data(nsh + 204);
    const auto *nsh_205 = buffer.data(nsh + 205);
    const auto *nsh_206 = buffer.data(nsh + 206);
    const auto *nsh_207 = buffer.data(nsh + 207);
    const auto *nsh_209 = buffer.data(nsh + 209);
    const auto *nsh_210 = buffer.data(nsh + 210);
    const auto *nsh_213 = buffer.data(nsh + 213);
    const auto *nsh_216 = buffer.data(nsh + 216);
    const auto *nsh_220 = buffer.data(nsh + 220);
    const auto *nsh_225 = buffer.data(nsh + 225);
    const auto *nsh_227 = buffer.data(nsh + 227);
    const auto *nsh_228 = buffer.data(nsh + 228);
    const auto *nsh_229 = buffer.data(nsh + 229);
    const auto *nsh_230 = buffer.data(nsh + 230);
    const auto *nsh_236 = buffer.data(nsh + 236);
    const auto *nsh_240 = buffer.data(nsh + 240);
    const auto *nsh_245 = buffer.data(nsh + 245);
    const auto *nsh_246 = buffer.data(nsh + 246);
    const auto *nsh_247 = buffer.data(nsh + 247);
    const auto *nsh_248 = buffer.data(nsh + 248);
    const auto *nsh_249 = buffer.data(nsh + 249);
    const auto *nsh_250 = buffer.data(nsh + 250);
    const auto *nsh_251 = buffer.data(nsh + 251);
    const auto *nsh_252 = buffer.data(nsh + 252);
    const auto *nsh_255 = buffer.data(nsh + 255);
    const auto *nsh_257 = buffer.data(nsh + 257);
    const auto *nsh_258 = buffer.data(nsh + 258);
    const auto *nsh_261 = buffer.data(nsh + 261);
    const auto *nsh_262 = buffer.data(nsh + 262);
    const auto *nsh_264 = buffer.data(nsh + 264);
    const auto *nsh_266 = buffer.data(nsh + 266);
    const auto *nsh_267 = buffer.data(nsh + 267);
    const auto *nsh_268 = buffer.data(nsh + 268);
    const auto *nsh_269 = buffer.data(nsh + 269);
    const auto *nsh_270 = buffer.data(nsh + 270);
    const auto *nsh_271 = buffer.data(nsh + 271);
    const auto *nsh_272 = buffer.data(nsh + 272);

    const auto *nsi1_167 = buffer.data(nsi1 + 167);
    const auto *nsi1_168 = buffer.data(nsi1 + 168);
    const auto *nsi1_171 = buffer.data(nsi1 + 171);
    const auto *nsi1_174 = buffer.data(nsi1 + 174);
    const auto *nsi1_178 = buffer.data(nsi1 + 178);
    const auto *nsi1_180 = buffer.data(nsi1 + 180);
    const auto *nsi1_189 = buffer.data(nsi1 + 189);

    const auto *osg0_133 = buffer.data(osg0 + 133);
    const auto *osg0_134 = buffer.data(osg0 + 134);
    const auto *osg0_135 = buffer.data(osg0 + 135);
    const auto *osg0_136 = buffer.data(osg0 + 136);
    const auto *osg0_137 = buffer.data(osg0 + 137);
    const auto *osg0_138 = buffer.data(osg0 + 138);
    const auto *osg0_139 = buffer.data(osg0 + 139);
    const auto *osg0_140 = buffer.data(osg0 + 140);
    const auto *osg0_144 = buffer.data(osg0 + 144);
    const auto *osg0_145 = buffer.data(osg0 + 145);
    const auto *osg0_146 = buffer.data(osg0 + 146);
    const auto *osg0_147 = buffer.data(osg0 + 147);
    const auto *osg0_148 = buffer.data(osg0 + 148);
    const auto *osg0_149 = buffer.data(osg0 + 149);
    const auto *osg0_150 = buffer.data(osg0 + 150);
    const auto *osg0_152 = buffer.data(osg0 + 152);
    const auto *osg0_153 = buffer.data(osg0 + 153);
    const auto *osg0_155 = buffer.data(osg0 + 155);
    const auto *osg0_156 = buffer.data(osg0 + 156);
    const auto *osg0_160 = buffer.data(osg0 + 160);
    const auto *osg0_161 = buffer.data(osg0 + 161);
    const auto *osg0_162 = buffer.data(osg0 + 162);
    const auto *osg0_164 = buffer.data(osg0 + 164);
    const auto *osg0_170 = buffer.data(osg0 + 170);
    const auto *osg0_174 = buffer.data(osg0 + 174);
    const auto *osg0_177 = buffer.data(osg0 + 177);
    const auto *osg0_178 = buffer.data(osg0 + 178);
    const auto *osg0_179 = buffer.data(osg0 + 179);
    const auto *osg0_180 = buffer.data(osg0 + 180);
    const auto *osg0_183 = buffer.data(osg0 + 183);
    const auto *osg0_185 = buffer.data(osg0 + 185);
    const auto *osg0_186 = buffer.data(osg0 + 186);
    const auto *osg0_189 = buffer.data(osg0 + 189);
    const auto *osg0_190 = buffer.data(osg0 + 190);
    const auto *osg0_192 = buffer.data(osg0 + 192);
    const auto *osg0_193 = buffer.data(osg0 + 193);
    const auto *osg0_194 = buffer.data(osg0 + 194);

    const auto *osg1_133 = buffer.data(osg1 + 133);
    const auto *osg1_134 = buffer.data(osg1 + 134);
    const auto *osg1_135 = buffer.data(osg1 + 135);
    const auto *osg1_136 = buffer.data(osg1 + 136);
    const auto *osg1_137 = buffer.data(osg1 + 137);
    const auto *osg1_138 = buffer.data(osg1 + 138);
    const auto *osg1_139 = buffer.data(osg1 + 139);
    const auto *osg1_140 = buffer.data(osg1 + 140);
    const auto *osg1_144 = buffer.data(osg1 + 144);
    const auto *osg1_145 = buffer.data(osg1 + 145);
    const auto *osg1_146 = buffer.data(osg1 + 146);
    const auto *osg1_147 = buffer.data(osg1 + 147);
    const auto *osg1_148 = buffer.data(osg1 + 148);
    const auto *osg1_149 = buffer.data(osg1 + 149);
    const auto *osg1_150 = buffer.data(osg1 + 150);
    const auto *osg1_152 = buffer.data(osg1 + 152);
    const auto *osg1_153 = buffer.data(osg1 + 153);
    const auto *osg1_155 = buffer.data(osg1 + 155);
    const auto *osg1_156 = buffer.data(osg1 + 156);
    const auto *osg1_160 = buffer.data(osg1 + 160);
    const auto *osg1_161 = buffer.data(osg1 + 161);
    const auto *osg1_162 = buffer.data(osg1 + 162);
    const auto *osg1_164 = buffer.data(osg1 + 164);
    const auto *osg1_170 = buffer.data(osg1 + 170);
    const auto *osg1_174 = buffer.data(osg1 + 174);
    const auto *osg1_177 = buffer.data(osg1 + 177);
    const auto *osg1_178 = buffer.data(osg1 + 178);
    const auto *osg1_179 = buffer.data(osg1 + 179);
    const auto *osg1_180 = buffer.data(osg1 + 180);
    const auto *osg1_183 = buffer.data(osg1 + 183);
    const auto *osg1_185 = buffer.data(osg1 + 185);
    const auto *osg1_186 = buffer.data(osg1 + 186);
    const auto *osg1_189 = buffer.data(osg1 + 189);
    const auto *osg1_190 = buffer.data(osg1 + 190);
    const auto *osg1_192 = buffer.data(osg1 + 192);
    const auto *osg1_193 = buffer.data(osg1 + 193);
    const auto *osg1_194 = buffer.data(osg1 + 194);

    const auto *osh_186 = buffer.data(osh + 186);
    const auto *osh_187 = buffer.data(osh + 187);
    const auto *osh_188 = buffer.data(osh + 188);
    const auto *osh_189 = buffer.data(osh + 189);
    const auto *osh_190 = buffer.data(osh + 190);
    const auto *osh_191 = buffer.data(osh + 191);
    const auto *osh_192 = buffer.data(osh + 192);
    const auto *osh_193 = buffer.data(osh + 193);
    const auto *osh_194 = buffer.data(osh + 194);
    const auto *osh_195 = buffer.data(osh + 195);
    const auto *osh_196 = buffer.data(osh + 196);
    const auto *osh_197 = buffer.data(osh + 197);
    const auto *osh_198 = buffer.data(osh + 198);
    const auto *osh_203 = buffer.data(osh + 203);
    const auto *osh_204 = buffer.data(osh + 204);
    const auto *osh_205 = buffer.data(osh + 205);
    const auto *osh_206 = buffer.data(osh + 206);
    const auto *osh_207 = buffer.data(osh + 207);
    const auto *osh_208 = buffer.data(osh + 208);
    const auto *osh_209 = buffer.data(osh + 209);
    const auto *osh_210 = buffer.data(osh + 210);
    const auto *osh_211 = buffer.data(osh + 211);
    const auto *osh_212 = buffer.data(osh + 212);
    const auto *osh_213 = buffer.data(osh + 213);
    const auto *osh_215 = buffer.data(osh + 215);
    const auto *osh_216 = buffer.data(osh + 216);
    const auto *osh_217 = buffer.data(osh + 217);
    const auto *osh_219 = buffer.data(osh + 219);
    const auto *osh_220 = buffer.data(osh + 220);
    const auto *osh_225 = buffer.data(osh + 225);
    const auto *osh_226 = buffer.data(osh + 226);
    const auto *osh_227 = buffer.data(osh + 227);
    const auto *osh_228 = buffer.data(osh + 228);
    const auto *osh_229 = buffer.data(osh + 229);
    const auto *osh_230 = buffer.data(osh + 230);
    const auto *osh_231 = buffer.data(osh + 231);
    const auto *osh_233 = buffer.data(osh + 233);
    const auto *osh_234 = buffer.data(osh + 234);
    const auto *osh_236 = buffer.data(osh + 236);
    const auto *osh_237 = buffer.data(osh + 237);
    const auto *osh_240 = buffer.data(osh + 240);
    const auto *osh_245 = buffer.data(osh + 245);
    const auto *osh_246 = buffer.data(osh + 246);
    const auto *osh_247 = buffer.data(osh + 247);
    const auto *osh_248 = buffer.data(osh + 248);
    const auto *osh_249 = buffer.data(osh + 249);
    const auto *osh_250 = buffer.data(osh + 250);
    const auto *osh_251 = buffer.data(osh + 251);
    const auto *osh_252 = buffer.data(osh + 252);
    const auto *osh_254 = buffer.data(osh + 254);
    const auto *osh_255 = buffer.data(osh + 255);
    const auto *osh_257 = buffer.data(osh + 257);
    const auto *osh_258 = buffer.data(osh + 258);
    const auto *osh_261 = buffer.data(osh + 261);
    const auto *osh_262 = buffer.data(osh + 262);
    const auto *osh_264 = buffer.data(osh + 264);
    const auto *osh_266 = buffer.data(osh + 266);
    const auto *osh_267 = buffer.data(osh + 267);
    const auto *osh_268 = buffer.data(osh + 268);
    const auto *osh_269 = buffer.data(osh + 269);
    const auto *osh_270 = buffer.data(osh + 270);
    const auto *osh_271 = buffer.data(osh + 271);
    const auto *osh_272 = buffer.data(osh + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, nsh_123, nsh_124, nsh_125, osg0_133, \
                         osg0_134, osg1_133, osg1_134, osh_186, osh_187, \
                         osh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * nsh_123[k]
                   + f_6 * osg0_133[k]
                   - f_7 * osg1_133[k]
                   + f_3 * pc_y[k] * osh_186[k];

        t_249[k] = f_11 * nsh_124[k]
                   + f_4 * osg0_134[k]
                   - f_5 * osg1_134[k]
                   + f_3 * pc_y[k] * osh_187[k];

        t_250[k] = f_11 * nsh_125[k]
                   + f_3 * pc_y[k] * osh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, nsi0_167, \
                         nsh_105, nsh_189, nsi1_167, osg0_135, osg1_135, \
                         osh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * nsi0_167[k]
                   - f_10 * pc_y[k] * nsi1_167[k];

        t_252[k] = f_19 * nsh_189[k]
                   + f_1 * osg0_135[k]
                   - f_2 * osg1_135[k]
                   + f_3 * pc_x[k] * osh_189[k];

        t_253[k] = f_3 * pc_y[k] * osh_189[k];

        t_254[k] = f_13 * nsh_105[k]
                   + f_3 * pc_z[k] * osh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, nsh_194, osg0_135, osg0_140, \
                         osg1_135, osg1_140, osh_190, osh_191, \
                         osh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * osg0_135[k]
                   - f_5 * osg1_135[k]
                   + f_3 * pc_y[k] * osh_190[k];

        t_256[k] = f_3 * pc_y[k] * osh_191[k];

        t_257[k] = f_19 * nsh_194[k]
                   + f_8 * osg0_140[k]
                   - f_9 * osg1_140[k]
                   + f_3 * pc_x[k] * osh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, osg0_136, osg0_137, osg1_136, osg1_137, \
                         osh_192, osh_193, osh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * osg0_136[k]
                   - f_7 * osg1_136[k]
                   + f_3 * pc_y[k] * osh_192[k];

        t_259[k] = f_4 * osg0_137[k]
                   - f_5 * osg1_137[k]
                   + f_3 * pc_y[k] * osh_193[k];

        t_260[k] = f_3 * pc_y[k] * osh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, nsh_198, osg0_138, osg0_139, \
                         osg0_144, osg1_138, osg1_139, osg1_144, osh_195, osh_196, \
                         osh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_19 * nsh_198[k]
                   + f_6 * osg0_144[k]
                   - f_7 * osg1_144[k]
                   + f_3 * pc_x[k] * osh_198[k];

        t_262[k] = f_8 * osg0_138[k]
                   - f_9 * osg1_138[k]
                   + f_3 * pc_y[k] * osh_195[k];

        t_263[k] = f_6 * osg0_139[k]
                   - f_7 * osg1_139[k]
                   + f_3 * pc_y[k] * osh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, nsh_203, nsh_204, osg0_140, \
                         osg0_149, osg1_140, osg1_149, osh_197, osh_198, osh_203, \
                         osh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * osg0_140[k]
                   - f_5 * osg1_140[k]
                   + f_3 * pc_y[k] * osh_197[k];

        t_265[k] = f_3 * pc_y[k] * osh_198[k];

        t_266[k] = f_19 * nsh_203[k]
                   + f_4 * osg0_149[k]
                   - f_5 * osg1_149[k]
                   + f_3 * pc_x[k] * osh_203[k];

        t_267[k] = f_19 * nsh_204[k]
                   + f_3 * pc_x[k] * osh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, nsh_205, nsh_206, \
                         nsh_207, nsh_209, osh_203, osh_205, osh_206, osh_207, \
                         osh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_19 * nsh_205[k]
                   + f_3 * pc_x[k] * osh_205[k];

        t_269[k] = f_19 * nsh_206[k]
                   + f_3 * pc_x[k] * osh_206[k];

        t_270[k] = f_19 * nsh_207[k]
                   + f_3 * pc_x[k] * osh_207[k];

        t_271[k] = f_3 * pc_y[k] * osh_203[k];

        t_272[k] = f_19 * nsh_209[k]
                   + f_3 * pc_x[k] * osh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, osg0_145, osg0_146, osg0_147, osg1_145, \
                         osg1_146, osg1_147, osh_204, osh_205, \
                         osh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * osg0_145[k]
                   - f_2 * osg1_145[k]
                   + f_3 * pc_y[k] * osh_204[k];

        t_274[k] = f_16 * osg0_146[k]
                   - f_17 * osg1_146[k]
                   + f_3 * pc_y[k] * osh_205[k];

        t_275[k] = f_8 * osg0_147[k]
                   - f_9 * osg1_147[k]
                   + f_3 * pc_y[k] * osh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, nsh_125, osg0_148, osg0_149, \
                         osg1_148, osg1_149, osh_207, osh_208, \
                         osh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * osg0_148[k]
                   - f_7 * osg1_148[k]
                   + f_3 * pc_y[k] * osh_207[k];

        t_277[k] = f_4 * osg0_149[k]
                   - f_5 * osg1_149[k]
                   + f_3 * pc_y[k] * osh_208[k];

        t_278[k] = f_3 * pc_y[k] * osh_209[k];

        t_279[k] = f_13 * nsh_125[k]
                   + f_1 * osg0_149[k]
                   - f_2 * osg1_149[k]
                   + f_3 * pc_z[k] * osh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, nsh_126, nsh_210, \
                         nsh_213, osg0_150, osg0_153, osg1_150, osg1_153, osh_210, \
                         osh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_20 * nsh_210[k]
                   + f_1 * osg0_150[k]
                   - f_2 * osg1_150[k]
                   + f_3 * pc_x[k] * osh_210[k];

        t_281[k] = f_14 * nsh_126[k]
                   + f_3 * pc_y[k] * osh_210[k];

        t_282[k] = f_3 * pc_z[k] * osh_210[k];

        t_283[k] = f_20 * nsh_213[k]
                   + f_8 * osg0_153[k]
                   - f_9 * osg1_153[k]
                   + f_3 * pc_x[k] * osh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, nsh_216, osg0_150, osg0_156, \
                         osg1_150, osg1_156, osh_211, osh_212, osh_213, \
                         osh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * osh_211[k];

        t_285[k] = f_4 * osg0_150[k]
                   - f_5 * osg1_150[k]
                   + f_3 * pc_z[k] * osh_212[k];

        t_286[k] = f_20 * nsh_216[k]
                   + f_6 * osg0_156[k]
                   - f_7 * osg1_156[k]
                   + f_3 * pc_x[k] * osh_216[k];

        t_287[k] = f_3 * pc_z[k] * osh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, nsh_131, nsh_220, \
                         osg0_152, osg0_160, osg1_152, osg1_160, osh_215, osh_216, \
                         osh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * nsh_131[k]
                   + f_3 * pc_y[k] * osh_215[k];

        t_289[k] = f_6 * osg0_152[k]
                   - f_7 * osg1_152[k]
                   + f_3 * pc_z[k] * osh_215[k];

        t_290[k] = f_20 * nsh_220[k]
                   + f_4 * osg0_160[k]
                   - f_5 * osg1_160[k]
                   + f_3 * pc_x[k] * osh_220[k];

        t_291[k] = f_3 * pc_z[k] * osh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, nsh_135, nsh_225, \
                         osg0_153, osg0_155, osg1_153, osg1_155, osh_217, osh_219, \
                         osh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * osg0_153[k]
                   - f_5 * osg1_153[k]
                   + f_3 * pc_z[k] * osh_217[k];

        t_293[k] = f_14 * nsh_135[k]
                   + f_3 * pc_y[k] * osh_219[k];

        t_294[k] = f_8 * osg0_155[k]
                   - f_9 * osg1_155[k]
                   + f_3 * pc_z[k] * osh_219[k];

        t_295[k] = f_20 * nsh_225[k]
                   + f_3 * pc_x[k] * osh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, nsh_227, nsh_228, \
                         nsh_229, nsh_230, osh_220, osh_227, osh_228, osh_229, \
                         osh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * osh_220[k];

        t_297[k] = f_20 * nsh_227[k]
                   + f_3 * pc_x[k] * osh_227[k];

        t_298[k] = f_20 * nsh_228[k]
                   + f_3 * pc_x[k] * osh_228[k];

        t_299[k] = f_20 * nsh_229[k]
                   + f_3 * pc_x[k] * osh_229[k];

        t_300[k] = f_20 * nsh_230[k]
                   + f_3 * pc_x[k] * osh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, nsh_141, osg0_160, osg0_161, \
                         osg1_160, osg1_161, osh_225, osh_226, \
                         osh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * nsh_141[k]
                   + f_1 * osg0_160[k]
                   - f_2 * osg1_160[k]
                   + f_3 * pc_y[k] * osh_225[k];

        t_302[k] = f_3 * pc_z[k] * osh_225[k];

        t_303[k] = f_4 * osg0_160[k]
                   - f_5 * osg1_160[k]
                   + f_3 * pc_z[k] * osh_226[k];

        t_304[k] = f_6 * osg0_161[k]
                   - f_7 * osg1_161[k]
                   + f_3 * pc_z[k] * osh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, nsi0_168, nsh_146, \
                         nsi1_168, osg0_162, osg0_164, osg1_162, osg1_164, osh_228, \
                         osh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * osg0_162[k]
                   - f_9 * osg1_162[k]
                   + f_3 * pc_z[k] * osh_228[k];

        t_306[k] = f_14 * nsh_146[k]
                   + f_3 * pc_y[k] * osh_230[k];

        t_307[k] = f_1 * osg0_164[k]
                   - f_2 * osg1_164[k]
                   + f_3 * pc_z[k] * osh_230[k];

        t_308[k] = pa_z[k] * nsi0_168[k]
                   - f_10 * pc_z[k] * nsi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, nsi0_171, nsh_126, \
                         nsh_147, nsh_149, nsi1_171, osh_231, osh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * nsh_147[k]
                   + f_3 * pc_y[k] * osh_231[k];

        t_310[k] = f_11 * nsh_126[k]
                   + f_3 * pc_z[k] * osh_231[k];

        t_311[k] = pa_z[k] * nsi0_171[k]
                   - f_10 * pc_z[k] * nsi1_171[k];

        t_312[k] = f_13 * nsh_149[k]
                   + f_3 * pc_y[k] * osh_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, nsi0_174, nsh_129, nsh_236, \
                         nsi1_174, osg0_170, osg1_170, osh_234, \
                         osh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_20 * nsh_236[k]
                   + f_8 * osg0_170[k]
                   - f_9 * osg1_170[k]
                   + f_3 * pc_x[k] * osh_236[k];

        t_314[k] = pa_z[k] * nsi0_174[k]
                   - f_10 * pc_z[k] * nsi1_174[k];

        t_315[k] = f_11 * nsh_129[k]
                   + f_3 * pc_z[k] * osh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, nsi0_178, nsh_152, \
                         nsh_240, nsi1_178, osg0_174, osg1_174, osh_236, \
                         osh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * nsh_152[k]
                   + f_3 * pc_y[k] * osh_236[k];

        t_317[k] = f_20 * nsh_240[k]
                   + f_6 * osg0_174[k]
                   - f_7 * osg1_174[k]
                   + f_3 * pc_x[k] * osh_240[k];

        t_318[k] = pa_z[k] * nsi0_178[k]
                   - f_10 * pc_z[k] * nsi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, nsi0_180, nsh_132, nsh_133, \
                         nsh_156, nsi1_180, osh_237, osh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * nsh_132[k]
                   + f_3 * pc_z[k] * osh_237[k];

        t_320[k] = pa_z[k] * nsi0_180[k]
                   + f_12 * nsh_133[k]
                   - f_10 * pc_z[k] * nsi1_180[k];

        t_321[k] = f_13 * nsh_156[k]
                   + f_3 * pc_y[k] * osh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, nsh_245, nsh_246, nsh_247, nsh_248, \
                         osg0_179, osg1_179, osh_245, osh_246, osh_247, \
                         osh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_20 * nsh_245[k]
                   + f_4 * osg0_179[k]
                   - f_5 * osg1_179[k]
                   + f_3 * pc_x[k] * osh_245[k];

        t_323[k] = f_20 * nsh_246[k]
                   + f_3 * pc_x[k] * osh_246[k];

        t_324[k] = f_20 * nsh_247[k]
                   + f_3 * pc_x[k] * osh_247[k];

        t_325[k] = f_20 * nsh_248[k]
                   + f_3 * pc_x[k] * osh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, nsi0_189, nsh_249, \
                         nsh_250, nsh_251, nsi1_189, osh_249, osh_250, \
                         osh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_20 * nsh_249[k]
                   + f_3 * pc_x[k] * osh_249[k];

        t_327[k] = f_20 * nsh_250[k]
                   + f_3 * pc_x[k] * osh_250[k];

        t_328[k] = f_20 * nsh_251[k]
                   + f_3 * pc_x[k] * osh_251[k];

        t_329[k] = pa_z[k] * nsi0_189[k]
                   - f_10 * pc_z[k] * nsi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, nsh_141, nsh_164, nsh_165, osg0_177, \
                         osg0_178, osg1_177, osg1_178, osh_246, osh_248, \
                         osh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * nsh_141[k]
                   + f_3 * pc_z[k] * osh_246[k];

        t_331[k] = f_13 * nsh_164[k]
                   + f_8 * osg0_177[k]
                   - f_9 * osg1_177[k]
                   + f_3 * pc_y[k] * osh_248[k];

        t_332[k] = f_13 * nsh_165[k]
                   + f_6 * osg0_178[k]
                   - f_7 * osg1_178[k]
                   + f_3 * pc_y[k] * osh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, nsh_146, nsh_166, nsh_167, osg0_179, \
                         osg1_179, osh_250, osh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * nsh_166[k]
                   + f_4 * osg0_179[k]
                   - f_5 * osg1_179[k]
                   + f_3 * pc_y[k] * osh_250[k];

        t_334[k] = f_13 * nsh_167[k]
                   + f_3 * pc_y[k] * osh_251[k];

        t_335[k] = f_11 * nsh_146[k]
                   + f_1 * osg0_179[k]
                   - f_2 * osg1_179[k]
                   + f_3 * pc_z[k] * osh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, nsh_147, nsh_168, nsh_252, \
                         osg0_180, osg1_180, osh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_20 * nsh_252[k]
                   + f_1 * osg0_180[k]
                   - f_2 * osg1_180[k]
                   + f_3 * pc_x[k] * osh_252[k];

        t_337[k] = f_12 * nsh_168[k]
                   + f_3 * pc_y[k] * osh_252[k];

        t_338[k] = f_12 * nsh_147[k]
                   + f_3 * pc_z[k] * osh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, nsh_170, nsh_255, nsh_257, osg0_183, \
                         osg0_185, osg1_183, osg1_185, osh_254, osh_255, \
                         osh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_20 * nsh_255[k]
                   + f_8 * osg0_183[k]
                   - f_9 * osg1_183[k]
                   + f_3 * pc_x[k] * osh_255[k];

        t_340[k] = f_12 * nsh_170[k]
                   + f_3 * pc_y[k] * osh_254[k];

        t_341[k] = f_20 * nsh_257[k]
                   + f_8 * osg0_185[k]
                   - f_9 * osg1_185[k]
                   + f_3 * pc_x[k] * osh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, nsh_150, nsh_173, nsh_258, \
                         osg0_186, osg1_186, osh_255, osh_257, \
                         osh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_20 * nsh_258[k]
                   + f_6 * osg0_186[k]
                   - f_7 * osg1_186[k]
                   + f_3 * pc_x[k] * osh_258[k];

        t_343[k] = f_12 * nsh_150[k]
                   + f_3 * pc_z[k] * osh_255[k];

        t_344[k] = f_12 * nsh_173[k]
                   + f_3 * pc_y[k] * osh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, nsh_153, nsh_261, nsh_262, osg0_189, \
                         osg0_190, osg1_189, osg1_190, osh_258, osh_261, \
                         osh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_20 * nsh_261[k]
                   + f_6 * osg0_189[k]
                   - f_7 * osg1_189[k]
                   + f_3 * pc_x[k] * osh_261[k];

        t_346[k] = f_20 * nsh_262[k]
                   + f_4 * osg0_190[k]
                   - f_5 * osg1_190[k]
                   + f_3 * pc_x[k] * osh_262[k];

        t_347[k] = f_12 * nsh_153[k]
                   + f_3 * pc_z[k] * osh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, nsh_177, nsh_264, nsh_266, osg0_192, \
                         osg0_194, osg1_192, osg1_194, osh_261, osh_264, \
                         osh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_20 * nsh_264[k]
                   + f_4 * osg0_192[k]
                   - f_5 * osg1_192[k]
                   + f_3 * pc_x[k] * osh_264[k];

        t_349[k] = f_12 * nsh_177[k]
                   + f_3 * pc_y[k] * osh_261[k];

        t_350[k] = f_20 * nsh_266[k]
                   + f_4 * osg0_194[k]
                   - f_5 * osg1_194[k]
                   + f_3 * pc_x[k] * osh_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, nsh_267, nsh_268, nsh_269, \
                         nsh_270, nsh_271, osh_267, osh_268, osh_269, osh_270, \
                         osh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_20 * nsh_267[k]
                   + f_3 * pc_x[k] * osh_267[k];

        t_352[k] = f_20 * nsh_268[k]
                   + f_3 * pc_x[k] * osh_268[k];

        t_353[k] = f_20 * nsh_269[k]
                   + f_3 * pc_x[k] * osh_269[k];

        t_354[k] = f_20 * nsh_270[k]
                   + f_3 * pc_x[k] * osh_270[k];

        t_355[k] = f_20 * nsh_271[k]
                   + f_3 * pc_x[k] * osh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, nsh_162, nsh_183, nsh_272, \
                         osg0_190, osg1_190, osh_267, osh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_20 * nsh_272[k]
                   + f_3 * pc_x[k] * osh_272[k];

        t_357[k] = f_12 * nsh_183[k]
                   + f_1 * osg0_190[k]
                   - f_2 * osg1_190[k]
                   + f_3 * pc_y[k] * osh_267[k];

        t_358[k] = f_12 * nsh_162[k]
                   + f_3 * pc_z[k] * osh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, nsh_185, nsh_186, nsh_187, osg0_192, \
                         osg0_193, osg0_194, osg1_192, osg1_193, osg1_194, osh_269, osh_270, \
                         osh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * nsh_185[k]
                   + f_8 * osg0_192[k]
                   - f_9 * osg1_192[k]
                   + f_3 * pc_y[k] * osh_269[k];

        t_360[k] = f_12 * nsh_186[k]
                   + f_6 * osg0_193[k]
                   - f_7 * osg1_193[k]
                   + f_3 * pc_y[k] * osh_270[k];

        t_361[k] = f_12 * nsh_187[k]
                   + f_4 * osg0_194[k]
                   - f_5 * osg1_194[k]
                   + f_3 * pc_y[k] * osh_271[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_252 = buffer.data(nsi0 + 252);
    const auto *nsi0_255 = buffer.data(nsi0 + 255);
    const auto *nsi0_257 = buffer.data(nsi0 + 257);
    const auto *nsi0_258 = buffer.data(nsi0 + 258);
    const auto *nsi0_261 = buffer.data(nsi0 + 261);
    const auto *nsi0_262 = buffer.data(nsi0 + 262);
    const auto *nsi0_264 = buffer.data(nsi0 + 264);
    const auto *nsi0_266 = buffer.data(nsi0 + 266);
    const auto *nsi0_279 = buffer.data(nsi0 + 279);
    const auto *nsi0_280 = buffer.data(nsi0 + 280);
    const auto *nsi0_283 = buffer.data(nsi0 + 283);
    const auto *nsi0_286 = buffer.data(nsi0 + 286);
    const auto *nsi0_290 = buffer.data(nsi0 + 290);
    const auto *nsi0_292 = buffer.data(nsi0 + 292);
    const auto *nsi0_301 = buffer.data(nsi0 + 301);

    const auto *nsh_167 = buffer.data(nsh + 167);
    const auto *nsh_168 = buffer.data(nsh + 168);
    const auto *nsh_171 = buffer.data(nsh + 171);
    const auto *nsh_174 = buffer.data(nsh + 174);
    const auto *nsh_183 = buffer.data(nsh + 183);
    const auto *nsh_188 = buffer.data(nsh + 188);
    const auto *nsh_189 = buffer.data(nsh + 189);
    const auto *nsh_190 = buffer.data(nsh + 190);
    const auto *nsh_191 = buffer.data(nsh + 191);
    const auto *nsh_192 = buffer.data(nsh + 192);
    const auto *nsh_194 = buffer.data(nsh + 194);
    const auto *nsh_195 = buffer.data(nsh + 195);
    const auto *nsh_197 = buffer.data(nsh + 197);
    const auto *nsh_198 = buffer.data(nsh + 198);
    const auto *nsh_204 = buffer.data(nsh + 204);
    const auto *nsh_206 = buffer.data(nsh + 206);
    const auto *nsh_207 = buffer.data(nsh + 207);
    const auto *nsh_208 = buffer.data(nsh + 208);
    const auto *nsh_209 = buffer.data(nsh + 209);
    const auto *nsh_210 = buffer.data(nsh + 210);
    const auto *nsh_213 = buffer.data(nsh + 213);
    const auto *nsh_215 = buffer.data(nsh + 215);
    const auto *nsh_216 = buffer.data(nsh + 216);
    const auto *nsh_217 = buffer.data(nsh + 217);
    const auto *nsh_219 = buffer.data(nsh + 219);
    const auto *nsh_225 = buffer.data(nsh + 225);
    const auto *nsh_230 = buffer.data(nsh + 230);
    const auto *nsh_231 = buffer.data(nsh + 231);
    const auto *nsh_233 = buffer.data(nsh + 233);
    const auto *nsh_236 = buffer.data(nsh + 236);
    const auto *nsh_240 = buffer.data(nsh + 240);
    const auto *nsh_248 = buffer.data(nsh + 248);
    const auto *nsh_249 = buffer.data(nsh + 249);
    const auto *nsh_250 = buffer.data(nsh + 250);
    const auto *nsh_251 = buffer.data(nsh + 251);
    const auto *nsh_252 = buffer.data(nsh + 252);
    const auto *nsh_288 = buffer.data(nsh + 288);
    const auto *nsh_289 = buffer.data(nsh + 289);
    const auto *nsh_290 = buffer.data(nsh + 290);
    const auto *nsh_291 = buffer.data(nsh + 291);
    const auto *nsh_292 = buffer.data(nsh + 292);
    const auto *nsh_293 = buffer.data(nsh + 293);
    const auto *nsh_294 = buffer.data(nsh + 294);
    const auto *nsh_299 = buffer.data(nsh + 299);
    const auto *nsh_303 = buffer.data(nsh + 303);
    const auto *nsh_308 = buffer.data(nsh + 308);
    const auto *nsh_309 = buffer.data(nsh + 309);
    const auto *nsh_310 = buffer.data(nsh + 310);
    const auto *nsh_311 = buffer.data(nsh + 311);
    const auto *nsh_312 = buffer.data(nsh + 312);
    const auto *nsh_314 = buffer.data(nsh + 314);
    const auto *nsh_315 = buffer.data(nsh + 315);
    const auto *nsh_318 = buffer.data(nsh + 318);
    const auto *nsh_321 = buffer.data(nsh + 321);
    const auto *nsh_325 = buffer.data(nsh + 325);
    const auto *nsh_330 = buffer.data(nsh + 330);
    const auto *nsh_332 = buffer.data(nsh + 332);
    const auto *nsh_333 = buffer.data(nsh + 333);
    const auto *nsh_334 = buffer.data(nsh + 334);
    const auto *nsh_335 = buffer.data(nsh + 335);
    const auto *nsh_341 = buffer.data(nsh + 341);
    const auto *nsh_345 = buffer.data(nsh + 345);
    const auto *nsh_350 = buffer.data(nsh + 350);
    const auto *nsh_351 = buffer.data(nsh + 351);
    const auto *nsh_352 = buffer.data(nsh + 352);
    const auto *nsh_353 = buffer.data(nsh + 353);
    const auto *nsh_354 = buffer.data(nsh + 354);
    const auto *nsh_355 = buffer.data(nsh + 355);
    const auto *nsh_356 = buffer.data(nsh + 356);
    const auto *nsh_357 = buffer.data(nsh + 357);

    const auto *nsi1_252 = buffer.data(nsi1 + 252);
    const auto *nsi1_255 = buffer.data(nsi1 + 255);
    const auto *nsi1_257 = buffer.data(nsi1 + 257);
    const auto *nsi1_258 = buffer.data(nsi1 + 258);
    const auto *nsi1_261 = buffer.data(nsi1 + 261);
    const auto *nsi1_262 = buffer.data(nsi1 + 262);
    const auto *nsi1_264 = buffer.data(nsi1 + 264);
    const auto *nsi1_266 = buffer.data(nsi1 + 266);
    const auto *nsi1_279 = buffer.data(nsi1 + 279);
    const auto *nsi1_280 = buffer.data(nsi1 + 280);
    const auto *nsi1_283 = buffer.data(nsi1 + 283);
    const auto *nsi1_286 = buffer.data(nsi1 + 286);
    const auto *nsi1_290 = buffer.data(nsi1 + 290);
    const auto *nsi1_292 = buffer.data(nsi1 + 292);
    const auto *nsi1_301 = buffer.data(nsi1 + 301);

    const auto *osg0_194 = buffer.data(osg0 + 194);
    const auto *osg0_205 = buffer.data(osg0 + 205);
    const auto *osg0_207 = buffer.data(osg0 + 207);
    const auto *osg0_208 = buffer.data(osg0 + 208);
    const auto *osg0_209 = buffer.data(osg0 + 209);
    const auto *osg0_210 = buffer.data(osg0 + 210);
    const auto *osg0_211 = buffer.data(osg0 + 211);
    const auto *osg0_212 = buffer.data(osg0 + 212);
    const auto *osg0_213 = buffer.data(osg0 + 213);
    const auto *osg0_214 = buffer.data(osg0 + 214);
    const auto *osg0_215 = buffer.data(osg0 + 215);
    const auto *osg0_219 = buffer.data(osg0 + 219);
    const auto *osg0_220 = buffer.data(osg0 + 220);
    const auto *osg0_221 = buffer.data(osg0 + 221);
    const auto *osg0_222 = buffer.data(osg0 + 222);
    const auto *osg0_223 = buffer.data(osg0 + 223);
    const auto *osg0_224 = buffer.data(osg0 + 224);
    const auto *osg0_225 = buffer.data(osg0 + 225);
    const auto *osg0_227 = buffer.data(osg0 + 227);
    const auto *osg0_228 = buffer.data(osg0 + 228);
    const auto *osg0_230 = buffer.data(osg0 + 230);
    const auto *osg0_231 = buffer.data(osg0 + 231);
    const auto *osg0_235 = buffer.data(osg0 + 235);
    const auto *osg0_236 = buffer.data(osg0 + 236);
    const auto *osg0_237 = buffer.data(osg0 + 237);
    const auto *osg0_239 = buffer.data(osg0 + 239);
    const auto *osg0_245 = buffer.data(osg0 + 245);
    const auto *osg0_249 = buffer.data(osg0 + 249);
    const auto *osg0_252 = buffer.data(osg0 + 252);
    const auto *osg0_253 = buffer.data(osg0 + 253);
    const auto *osg0_254 = buffer.data(osg0 + 254);
    const auto *osg0_255 = buffer.data(osg0 + 255);

    const auto *osg1_194 = buffer.data(osg1 + 194);
    const auto *osg1_205 = buffer.data(osg1 + 205);
    const auto *osg1_207 = buffer.data(osg1 + 207);
    const auto *osg1_208 = buffer.data(osg1 + 208);
    const auto *osg1_209 = buffer.data(osg1 + 209);
    const auto *osg1_210 = buffer.data(osg1 + 210);
    const auto *osg1_211 = buffer.data(osg1 + 211);
    const auto *osg1_212 = buffer.data(osg1 + 212);
    const auto *osg1_213 = buffer.data(osg1 + 213);
    const auto *osg1_214 = buffer.data(osg1 + 214);
    const auto *osg1_215 = buffer.data(osg1 + 215);
    const auto *osg1_219 = buffer.data(osg1 + 219);
    const auto *osg1_220 = buffer.data(osg1 + 220);
    const auto *osg1_221 = buffer.data(osg1 + 221);
    const auto *osg1_222 = buffer.data(osg1 + 222);
    const auto *osg1_223 = buffer.data(osg1 + 223);
    const auto *osg1_224 = buffer.data(osg1 + 224);
    const auto *osg1_225 = buffer.data(osg1 + 225);
    const auto *osg1_227 = buffer.data(osg1 + 227);
    const auto *osg1_228 = buffer.data(osg1 + 228);
    const auto *osg1_230 = buffer.data(osg1 + 230);
    const auto *osg1_231 = buffer.data(osg1 + 231);
    const auto *osg1_235 = buffer.data(osg1 + 235);
    const auto *osg1_236 = buffer.data(osg1 + 236);
    const auto *osg1_237 = buffer.data(osg1 + 237);
    const auto *osg1_239 = buffer.data(osg1 + 239);
    const auto *osg1_245 = buffer.data(osg1 + 245);
    const auto *osg1_249 = buffer.data(osg1 + 249);
    const auto *osg1_252 = buffer.data(osg1 + 252);
    const auto *osg1_253 = buffer.data(osg1 + 253);
    const auto *osg1_254 = buffer.data(osg1 + 254);
    const auto *osg1_255 = buffer.data(osg1 + 255);

    const auto *osh_272 = buffer.data(osh + 272);
    const auto *osh_273 = buffer.data(osh + 273);
    const auto *osh_275 = buffer.data(osh + 275);
    const auto *osh_276 = buffer.data(osh + 276);
    const auto *osh_278 = buffer.data(osh + 278);
    const auto *osh_279 = buffer.data(osh + 279);
    const auto *osh_282 = buffer.data(osh + 282);
    const auto *osh_288 = buffer.data(osh + 288);
    const auto *osh_289 = buffer.data(osh + 289);
    const auto *osh_290 = buffer.data(osh + 290);
    const auto *osh_291 = buffer.data(osh + 291);
    const auto *osh_292 = buffer.data(osh + 292);
    const auto *osh_293 = buffer.data(osh + 293);
    const auto *osh_294 = buffer.data(osh + 294);
    const auto *osh_295 = buffer.data(osh + 295);
    const auto *osh_296 = buffer.data(osh + 296);
    const auto *osh_297 = buffer.data(osh + 297);
    const auto *osh_298 = buffer.data(osh + 298);
    const auto *osh_299 = buffer.data(osh + 299);
    const auto *osh_300 = buffer.data(osh + 300);
    const auto *osh_301 = buffer.data(osh + 301);
    const auto *osh_302 = buffer.data(osh + 302);
    const auto *osh_303 = buffer.data(osh + 303);
    const auto *osh_308 = buffer.data(osh + 308);
    const auto *osh_309 = buffer.data(osh + 309);
    const auto *osh_310 = buffer.data(osh + 310);
    const auto *osh_311 = buffer.data(osh + 311);
    const auto *osh_312 = buffer.data(osh + 312);
    const auto *osh_313 = buffer.data(osh + 313);
    const auto *osh_314 = buffer.data(osh + 314);
    const auto *osh_315 = buffer.data(osh + 315);
    const auto *osh_316 = buffer.data(osh + 316);
    const auto *osh_317 = buffer.data(osh + 317);
    const auto *osh_318 = buffer.data(osh + 318);
    const auto *osh_320 = buffer.data(osh + 320);
    const auto *osh_321 = buffer.data(osh + 321);
    const auto *osh_322 = buffer.data(osh + 322);
    const auto *osh_324 = buffer.data(osh + 324);
    const auto *osh_325 = buffer.data(osh + 325);
    const auto *osh_330 = buffer.data(osh + 330);
    const auto *osh_331 = buffer.data(osh + 331);
    const auto *osh_332 = buffer.data(osh + 332);
    const auto *osh_333 = buffer.data(osh + 333);
    const auto *osh_334 = buffer.data(osh + 334);
    const auto *osh_335 = buffer.data(osh + 335);
    const auto *osh_336 = buffer.data(osh + 336);
    const auto *osh_338 = buffer.data(osh + 338);
    const auto *osh_339 = buffer.data(osh + 339);
    const auto *osh_341 = buffer.data(osh + 341);
    const auto *osh_342 = buffer.data(osh + 342);
    const auto *osh_345 = buffer.data(osh + 345);
    const auto *osh_350 = buffer.data(osh + 350);
    const auto *osh_351 = buffer.data(osh + 351);
    const auto *osh_352 = buffer.data(osh + 352);
    const auto *osh_353 = buffer.data(osh + 353);
    const auto *osh_354 = buffer.data(osh + 354);
    const auto *osh_355 = buffer.data(osh + 355);
    const auto *osh_356 = buffer.data(osh + 356);
    const auto *osh_357 = buffer.data(osh + 357);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, nsi0_252, nsh_167, \
                         nsh_188, nsh_189, nsi1_252, osg0_194, osg1_194, osh_272, \
                         osh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * nsh_188[k]
                   + f_3 * pc_y[k] * osh_272[k];

        t_363[k] = f_12 * nsh_167[k]
                   + f_1 * osg0_194[k]
                   - f_2 * osg1_194[k]
                   + f_3 * pc_z[k] * osh_272[k];

        t_364[k] = pa_y[k] * nsi0_252[k]
                   - f_10 * pc_y[k] * nsi1_252[k];

        t_365[k] = f_11 * nsh_189[k]
                   + f_3 * pc_y[k] * osh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, nsi0_255, nsi0_257, \
                         nsh_168, nsh_190, nsh_191, nsi1_255, nsi1_257, osh_273, \
                         osh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * nsh_168[k]
                   + f_3 * pc_z[k] * osh_273[k];

        t_367[k] = pa_y[k] * nsi0_255[k]
                   + f_12 * nsh_190[k]
                   - f_10 * pc_y[k] * nsi1_255[k];

        t_368[k] = f_11 * nsh_191[k]
                   + f_3 * pc_y[k] * osh_275[k];

        t_369[k] = pa_y[k] * nsi0_257[k]
                   - f_10 * pc_y[k] * nsi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, nsi0_258, nsi0_261, \
                         nsh_171, nsh_192, nsh_194, nsi1_258, nsi1_261, osh_276, \
                         osh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * nsi0_258[k]
                   + f_13 * nsh_192[k]
                   - f_10 * pc_y[k] * nsi1_258[k];

        t_371[k] = f_13 * nsh_171[k]
                   + f_3 * pc_z[k] * osh_276[k];

        t_372[k] = f_11 * nsh_194[k]
                   + f_3 * pc_y[k] * osh_278[k];

        t_373[k] = pa_y[k] * nsi0_261[k]
                   - f_10 * pc_y[k] * nsi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, nsi0_262, nsi0_264, nsh_174, \
                         nsh_195, nsh_197, nsi1_262, nsi1_264, \
                         osh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * nsi0_262[k]
                   + f_14 * nsh_195[k]
                   - f_10 * pc_y[k] * nsi1_262[k];

        t_375[k] = f_13 * nsh_174[k]
                   + f_3 * pc_z[k] * osh_279[k];

        t_376[k] = pa_y[k] * nsi0_264[k]
                   + f_12 * nsh_197[k]
                   - f_10 * pc_y[k] * nsi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, nsi0_266, nsh_198, \
                         nsh_288, nsh_289, nsi1_266, osh_282, osh_288, \
                         osh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * nsh_198[k]
                   + f_3 * pc_y[k] * osh_282[k];

        t_378[k] = pa_y[k] * nsi0_266[k]
                   - f_10 * pc_y[k] * nsi1_266[k];

        t_379[k] = f_20 * nsh_288[k]
                   + f_3 * pc_x[k] * osh_288[k];

        t_380[k] = f_20 * nsh_289[k]
                   + f_3 * pc_x[k] * osh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, nsh_290, nsh_291, nsh_292, nsh_293, \
                         osh_290, osh_291, osh_292, osh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_20 * nsh_290[k]
                   + f_3 * pc_x[k] * osh_290[k];

        t_382[k] = f_20 * nsh_291[k]
                   + f_3 * pc_x[k] * osh_291[k];

        t_383[k] = f_20 * nsh_292[k]
                   + f_3 * pc_x[k] * osh_292[k];

        t_384[k] = f_20 * nsh_293[k]
                   + f_3 * pc_x[k] * osh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, nsh_183, nsh_204, nsh_206, osg0_205, \
                         osg0_207, osg1_205, osg1_207, osh_288, \
                         osh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * nsh_204[k]
                   + f_1 * osg0_205[k]
                   - f_2 * osg1_205[k]
                   + f_3 * pc_y[k] * osh_288[k];

        t_386[k] = f_13 * nsh_183[k]
                   + f_3 * pc_z[k] * osh_288[k];

        t_387[k] = f_11 * nsh_206[k]
                   + f_8 * osg0_207[k]
                   - f_9 * osg1_207[k]
                   + f_3 * pc_y[k] * osh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, nsh_207, nsh_208, nsh_209, osg0_208, \
                         osg0_209, osg1_208, osg1_209, osh_291, osh_292, \
                         osh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * nsh_207[k]
                   + f_6 * osg0_208[k]
                   - f_7 * osg1_208[k]
                   + f_3 * pc_y[k] * osh_291[k];

        t_389[k] = f_11 * nsh_208[k]
                   + f_4 * osg0_209[k]
                   - f_5 * osg1_209[k]
                   + f_3 * pc_y[k] * osh_292[k];

        t_390[k] = f_11 * nsh_209[k]
                   + f_3 * pc_y[k] * osh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, nsi0_279, \
                         nsh_189, nsh_294, nsi1_279, osg0_210, osg1_210, \
                         osh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * nsi0_279[k]
                   - f_10 * pc_y[k] * nsi1_279[k];

        t_392[k] = f_20 * nsh_294[k]
                   + f_1 * osg0_210[k]
                   - f_2 * osg1_210[k]
                   + f_3 * pc_x[k] * osh_294[k];

        t_393[k] = f_3 * pc_y[k] * osh_294[k];

        t_394[k] = f_14 * nsh_189[k]
                   + f_3 * pc_z[k] * osh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, nsh_299, osg0_210, osg0_215, \
                         osg1_210, osg1_215, osh_295, osh_296, \
                         osh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * osg0_210[k]
                   - f_5 * osg1_210[k]
                   + f_3 * pc_y[k] * osh_295[k];

        t_396[k] = f_3 * pc_y[k] * osh_296[k];

        t_397[k] = f_20 * nsh_299[k]
                   + f_8 * osg0_215[k]
                   - f_9 * osg1_215[k]
                   + f_3 * pc_x[k] * osh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, osg0_211, osg0_212, osg1_211, osg1_212, \
                         osh_297, osh_298, osh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * osg0_211[k]
                   - f_7 * osg1_211[k]
                   + f_3 * pc_y[k] * osh_297[k];

        t_399[k] = f_4 * osg0_212[k]
                   - f_5 * osg1_212[k]
                   + f_3 * pc_y[k] * osh_298[k];

        t_400[k] = f_3 * pc_y[k] * osh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, nsh_303, osg0_213, osg0_214, \
                         osg0_219, osg1_213, osg1_214, osg1_219, osh_300, osh_301, \
                         osh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_20 * nsh_303[k]
                   + f_6 * osg0_219[k]
                   - f_7 * osg1_219[k]
                   + f_3 * pc_x[k] * osh_303[k];

        t_402[k] = f_8 * osg0_213[k]
                   - f_9 * osg1_213[k]
                   + f_3 * pc_y[k] * osh_300[k];

        t_403[k] = f_6 * osg0_214[k]
                   - f_7 * osg1_214[k]
                   + f_3 * pc_y[k] * osh_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, nsh_308, nsh_309, osg0_215, \
                         osg0_224, osg1_215, osg1_224, osh_302, osh_303, osh_308, \
                         osh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * osg0_215[k]
                   - f_5 * osg1_215[k]
                   + f_3 * pc_y[k] * osh_302[k];

        t_405[k] = f_3 * pc_y[k] * osh_303[k];

        t_406[k] = f_20 * nsh_308[k]
                   + f_4 * osg0_224[k]
                   - f_5 * osg1_224[k]
                   + f_3 * pc_x[k] * osh_308[k];

        t_407[k] = f_20 * nsh_309[k]
                   + f_3 * pc_x[k] * osh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, nsh_310, nsh_311, \
                         nsh_312, nsh_314, osh_308, osh_310, osh_311, osh_312, \
                         osh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_20 * nsh_310[k]
                   + f_3 * pc_x[k] * osh_310[k];

        t_409[k] = f_20 * nsh_311[k]
                   + f_3 * pc_x[k] * osh_311[k];

        t_410[k] = f_20 * nsh_312[k]
                   + f_3 * pc_x[k] * osh_312[k];

        t_411[k] = f_3 * pc_y[k] * osh_308[k];

        t_412[k] = f_20 * nsh_314[k]
                   + f_3 * pc_x[k] * osh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, osg0_220, osg0_221, osg0_222, osg1_220, \
                         osg1_221, osg1_222, osh_309, osh_310, \
                         osh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * osg0_220[k]
                   - f_2 * osg1_220[k]
                   + f_3 * pc_y[k] * osh_309[k];

        t_414[k] = f_16 * osg0_221[k]
                   - f_17 * osg1_221[k]
                   + f_3 * pc_y[k] * osh_310[k];

        t_415[k] = f_8 * osg0_222[k]
                   - f_9 * osg1_222[k]
                   + f_3 * pc_y[k] * osh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, nsh_209, osg0_223, osg0_224, \
                         osg1_223, osg1_224, osh_312, osh_313, \
                         osh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * osg0_223[k]
                   - f_7 * osg1_223[k]
                   + f_3 * pc_y[k] * osh_312[k];

        t_417[k] = f_4 * osg0_224[k]
                   - f_5 * osg1_224[k]
                   + f_3 * pc_y[k] * osh_313[k];

        t_418[k] = f_3 * pc_y[k] * osh_314[k];

        t_419[k] = f_14 * nsh_209[k]
                   + f_1 * osg0_224[k]
                   - f_2 * osg1_224[k]
                   + f_3 * pc_z[k] * osh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, nsh_210, nsh_315, \
                         nsh_318, osg0_225, osg0_228, osg1_225, osg1_228, osh_315, \
                         osh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_21 * nsh_315[k]
                   + f_1 * osg0_225[k]
                   - f_2 * osg1_225[k]
                   + f_3 * pc_x[k] * osh_315[k];

        t_421[k] = f_22 * nsh_210[k]
                   + f_3 * pc_y[k] * osh_315[k];

        t_422[k] = f_3 * pc_z[k] * osh_315[k];

        t_423[k] = f_21 * nsh_318[k]
                   + f_8 * osg0_228[k]
                   - f_9 * osg1_228[k]
                   + f_3 * pc_x[k] * osh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pc_x, pc_z, nsh_321, osg0_225, osg0_231, \
                         osg1_225, osg1_231, osh_316, osh_317, osh_318, \
                         osh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * osh_316[k];

        t_425[k] = f_4 * osg0_225[k]
                   - f_5 * osg1_225[k]
                   + f_3 * pc_z[k] * osh_317[k];

        t_426[k] = f_21 * nsh_321[k]
                   + f_6 * osg0_231[k]
                   - f_7 * osg1_231[k]
                   + f_3 * pc_x[k] * osh_321[k];

        t_427[k] = f_3 * pc_z[k] * osh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, nsh_215, nsh_325, \
                         osg0_227, osg0_235, osg1_227, osg1_235, osh_320, osh_321, \
                         osh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_22 * nsh_215[k]
                   + f_3 * pc_y[k] * osh_320[k];

        t_429[k] = f_6 * osg0_227[k]
                   - f_7 * osg1_227[k]
                   + f_3 * pc_z[k] * osh_320[k];

        t_430[k] = f_21 * nsh_325[k]
                   + f_4 * osg0_235[k]
                   - f_5 * osg1_235[k]
                   + f_3 * pc_x[k] * osh_325[k];

        t_431[k] = f_3 * pc_z[k] * osh_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, nsh_219, nsh_330, \
                         osg0_228, osg0_230, osg1_228, osg1_230, osh_322, osh_324, \
                         osh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * osg0_228[k]
                   - f_5 * osg1_228[k]
                   + f_3 * pc_z[k] * osh_322[k];

        t_433[k] = f_22 * nsh_219[k]
                   + f_3 * pc_y[k] * osh_324[k];

        t_434[k] = f_8 * osg0_230[k]
                   - f_9 * osg1_230[k]
                   + f_3 * pc_z[k] * osh_324[k];

        t_435[k] = f_21 * nsh_330[k]
                   + f_3 * pc_x[k] * osh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, nsh_332, nsh_333, \
                         nsh_334, nsh_335, osh_325, osh_332, osh_333, osh_334, \
                         osh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * osh_325[k];

        t_437[k] = f_21 * nsh_332[k]
                   + f_3 * pc_x[k] * osh_332[k];

        t_438[k] = f_21 * nsh_333[k]
                   + f_3 * pc_x[k] * osh_333[k];

        t_439[k] = f_21 * nsh_334[k]
                   + f_3 * pc_x[k] * osh_334[k];

        t_440[k] = f_21 * nsh_335[k]
                   + f_3 * pc_x[k] * osh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, nsh_225, osg0_235, osg0_236, \
                         osg1_235, osg1_236, osh_330, osh_331, \
                         osh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_22 * nsh_225[k]
                   + f_1 * osg0_235[k]
                   - f_2 * osg1_235[k]
                   + f_3 * pc_y[k] * osh_330[k];

        t_442[k] = f_3 * pc_z[k] * osh_330[k];

        t_443[k] = f_4 * osg0_235[k]
                   - f_5 * osg1_235[k]
                   + f_3 * pc_z[k] * osh_331[k];

        t_444[k] = f_6 * osg0_236[k]
                   - f_7 * osg1_236[k]
                   + f_3 * pc_z[k] * osh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, nsi0_280, nsh_230, \
                         nsi1_280, osg0_237, osg0_239, osg1_237, osg1_239, osh_333, \
                         osh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * osg0_237[k]
                   - f_9 * osg1_237[k]
                   + f_3 * pc_z[k] * osh_333[k];

        t_446[k] = f_22 * nsh_230[k]
                   + f_3 * pc_y[k] * osh_335[k];

        t_447[k] = f_1 * osg0_239[k]
                   - f_2 * osg1_239[k]
                   + f_3 * pc_z[k] * osh_335[k];

        t_448[k] = pa_z[k] * nsi0_280[k]
                   - f_10 * pc_z[k] * nsi1_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_z, pc_y, pc_z, nsi0_283, nsh_210, \
                         nsh_231, nsh_233, nsi1_283, osh_336, osh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * nsh_231[k]
                   + f_3 * pc_y[k] * osh_336[k];

        t_450[k] = f_11 * nsh_210[k]
                   + f_3 * pc_z[k] * osh_336[k];

        t_451[k] = pa_z[k] * nsi0_283[k]
                   - f_10 * pc_z[k] * nsi1_283[k];

        t_452[k] = f_14 * nsh_233[k]
                   + f_3 * pc_y[k] * osh_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pc_x, pc_z, nsi0_286, nsh_213, nsh_341, \
                         nsi1_286, osg0_245, osg1_245, osh_339, \
                         osh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_21 * nsh_341[k]
                   + f_8 * osg0_245[k]
                   - f_9 * osg1_245[k]
                   + f_3 * pc_x[k] * osh_341[k];

        t_454[k] = pa_z[k] * nsi0_286[k]
                   - f_10 * pc_z[k] * nsi1_286[k];

        t_455[k] = f_11 * nsh_213[k]
                   + f_3 * pc_z[k] * osh_339[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pc_x, pc_y, pc_z, nsi0_290, nsh_236, \
                         nsh_345, nsi1_290, osg0_249, osg1_249, osh_341, \
                         osh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * nsh_236[k]
                   + f_3 * pc_y[k] * osh_341[k];

        t_457[k] = f_21 * nsh_345[k]
                   + f_6 * osg0_249[k]
                   - f_7 * osg1_249[k]
                   + f_3 * pc_x[k] * osh_345[k];

        t_458[k] = pa_z[k] * nsi0_290[k]
                   - f_10 * pc_z[k] * nsi1_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pc_y, pc_z, nsi0_292, nsh_216, nsh_217, \
                         nsh_240, nsi1_292, osh_342, osh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * nsh_216[k]
                   + f_3 * pc_z[k] * osh_342[k];

        t_460[k] = pa_z[k] * nsi0_292[k]
                   + f_12 * nsh_217[k]
                   - f_10 * pc_z[k] * nsi1_292[k];

        t_461[k] = f_14 * nsh_240[k]
                   + f_3 * pc_y[k] * osh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, nsh_350, nsh_351, nsh_352, nsh_353, \
                         osg0_254, osg1_254, osh_350, osh_351, osh_352, \
                         osh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_21 * nsh_350[k]
                   + f_4 * osg0_254[k]
                   - f_5 * osg1_254[k]
                   + f_3 * pc_x[k] * osh_350[k];

        t_463[k] = f_21 * nsh_351[k]
                   + f_3 * pc_x[k] * osh_351[k];

        t_464[k] = f_21 * nsh_352[k]
                   + f_3 * pc_x[k] * osh_352[k];

        t_465[k] = f_21 * nsh_353[k]
                   + f_3 * pc_x[k] * osh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_x, pc_z, nsi0_301, nsh_354, \
                         nsh_355, nsh_356, nsi1_301, osh_354, osh_355, \
                         osh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_21 * nsh_354[k]
                   + f_3 * pc_x[k] * osh_354[k];

        t_467[k] = f_21 * nsh_355[k]
                   + f_3 * pc_x[k] * osh_355[k];

        t_468[k] = f_21 * nsh_356[k]
                   + f_3 * pc_x[k] * osh_356[k];

        t_469[k] = pa_z[k] * nsi0_301[k]
                   - f_10 * pc_z[k] * nsi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, nsh_225, nsh_248, nsh_249, osg0_252, \
                         osg0_253, osg1_252, osg1_253, osh_351, osh_353, \
                         osh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * nsh_225[k]
                   + f_3 * pc_z[k] * osh_351[k];

        t_471[k] = f_14 * nsh_248[k]
                   + f_8 * osg0_252[k]
                   - f_9 * osg1_252[k]
                   + f_3 * pc_y[k] * osh_353[k];

        t_472[k] = f_14 * nsh_249[k]
                   + f_6 * osg0_253[k]
                   - f_7 * osg1_253[k]
                   + f_3 * pc_y[k] * osh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, nsh_230, nsh_250, nsh_251, osg0_254, \
                         osg1_254, osh_355, osh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * nsh_250[k]
                   + f_4 * osg0_254[k]
                   - f_5 * osg1_254[k]
                   + f_3 * pc_y[k] * osh_355[k];

        t_474[k] = f_14 * nsh_251[k]
                   + f_3 * pc_y[k] * osh_356[k];

        t_475[k] = f_11 * nsh_230[k]
                   + f_1 * osg0_254[k]
                   - f_2 * osg1_254[k]
                   + f_3 * pc_z[k] * osh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, nsh_231, nsh_252, nsh_357, \
                         osg0_255, osg1_255, osh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_21 * nsh_357[k]
                   + f_1 * osg0_255[k]
                   - f_2 * osg1_255[k]
                   + f_3 * pc_x[k] * osh_357[k];

        t_477[k] = f_13 * nsh_252[k]
                   + f_3 * pc_y[k] * osh_357[k];

        t_478[k] = f_12 * nsh_231[k]
                   + f_3 * pc_z[k] * osh_357[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_392 = buffer.data(nsi0 + 392);
    const auto *nsi0_395 = buffer.data(nsi0 + 395);
    const auto *nsi0_397 = buffer.data(nsi0 + 397);
    const auto *nsi0_398 = buffer.data(nsi0 + 398);
    const auto *nsi0_401 = buffer.data(nsi0 + 401);
    const auto *nsi0_402 = buffer.data(nsi0 + 402);
    const auto *nsi0_404 = buffer.data(nsi0 + 404);
    const auto *nsi0_406 = buffer.data(nsi0 + 406);
    const auto *nsi0_419 = buffer.data(nsi0 + 419);

    const auto *nsh_234 = buffer.data(nsh + 234);
    const auto *nsh_237 = buffer.data(nsh + 237);
    const auto *nsh_246 = buffer.data(nsh + 246);
    const auto *nsh_251 = buffer.data(nsh + 251);
    const auto *nsh_252 = buffer.data(nsh + 252);
    const auto *nsh_254 = buffer.data(nsh + 254);
    const auto *nsh_255 = buffer.data(nsh + 255);
    const auto *nsh_257 = buffer.data(nsh + 257);
    const auto *nsh_258 = buffer.data(nsh + 258);
    const auto *nsh_261 = buffer.data(nsh + 261);
    const auto *nsh_267 = buffer.data(nsh + 267);
    const auto *nsh_269 = buffer.data(nsh + 269);
    const auto *nsh_270 = buffer.data(nsh + 270);
    const auto *nsh_271 = buffer.data(nsh + 271);
    const auto *nsh_272 = buffer.data(nsh + 272);
    const auto *nsh_273 = buffer.data(nsh + 273);
    const auto *nsh_275 = buffer.data(nsh + 275);
    const auto *nsh_276 = buffer.data(nsh + 276);
    const auto *nsh_278 = buffer.data(nsh + 278);
    const auto *nsh_279 = buffer.data(nsh + 279);
    const auto *nsh_282 = buffer.data(nsh + 282);
    const auto *nsh_288 = buffer.data(nsh + 288);
    const auto *nsh_290 = buffer.data(nsh + 290);
    const auto *nsh_291 = buffer.data(nsh + 291);
    const auto *nsh_292 = buffer.data(nsh + 292);
    const auto *nsh_293 = buffer.data(nsh + 293);
    const auto *nsh_294 = buffer.data(nsh + 294);
    const auto *nsh_295 = buffer.data(nsh + 295);
    const auto *nsh_296 = buffer.data(nsh + 296);
    const auto *nsh_297 = buffer.data(nsh + 297);
    const auto *nsh_299 = buffer.data(nsh + 299);
    const auto *nsh_300 = buffer.data(nsh + 300);
    const auto *nsh_302 = buffer.data(nsh + 302);
    const auto *nsh_303 = buffer.data(nsh + 303);
    const auto *nsh_309 = buffer.data(nsh + 309);
    const auto *nsh_311 = buffer.data(nsh + 311);
    const auto *nsh_312 = buffer.data(nsh + 312);
    const auto *nsh_313 = buffer.data(nsh + 313);
    const auto *nsh_314 = buffer.data(nsh + 314);
    const auto *nsh_315 = buffer.data(nsh + 315);
    const auto *nsh_360 = buffer.data(nsh + 360);
    const auto *nsh_362 = buffer.data(nsh + 362);
    const auto *nsh_363 = buffer.data(nsh + 363);
    const auto *nsh_366 = buffer.data(nsh + 366);
    const auto *nsh_367 = buffer.data(nsh + 367);
    const auto *nsh_369 = buffer.data(nsh + 369);
    const auto *nsh_371 = buffer.data(nsh + 371);
    const auto *nsh_372 = buffer.data(nsh + 372);
    const auto *nsh_373 = buffer.data(nsh + 373);
    const auto *nsh_374 = buffer.data(nsh + 374);
    const auto *nsh_375 = buffer.data(nsh + 375);
    const auto *nsh_376 = buffer.data(nsh + 376);
    const auto *nsh_377 = buffer.data(nsh + 377);
    const auto *nsh_378 = buffer.data(nsh + 378);
    const auto *nsh_381 = buffer.data(nsh + 381);
    const auto *nsh_383 = buffer.data(nsh + 383);
    const auto *nsh_384 = buffer.data(nsh + 384);
    const auto *nsh_387 = buffer.data(nsh + 387);
    const auto *nsh_388 = buffer.data(nsh + 388);
    const auto *nsh_390 = buffer.data(nsh + 390);
    const auto *nsh_392 = buffer.data(nsh + 392);
    const auto *nsh_393 = buffer.data(nsh + 393);
    const auto *nsh_394 = buffer.data(nsh + 394);
    const auto *nsh_395 = buffer.data(nsh + 395);
    const auto *nsh_396 = buffer.data(nsh + 396);
    const auto *nsh_397 = buffer.data(nsh + 397);
    const auto *nsh_398 = buffer.data(nsh + 398);
    const auto *nsh_414 = buffer.data(nsh + 414);
    const auto *nsh_415 = buffer.data(nsh + 415);
    const auto *nsh_416 = buffer.data(nsh + 416);
    const auto *nsh_417 = buffer.data(nsh + 417);
    const auto *nsh_418 = buffer.data(nsh + 418);
    const auto *nsh_419 = buffer.data(nsh + 419);
    const auto *nsh_420 = buffer.data(nsh + 420);
    const auto *nsh_425 = buffer.data(nsh + 425);
    const auto *nsh_429 = buffer.data(nsh + 429);
    const auto *nsh_434 = buffer.data(nsh + 434);
    const auto *nsh_435 = buffer.data(nsh + 435);
    const auto *nsh_436 = buffer.data(nsh + 436);
    const auto *nsh_437 = buffer.data(nsh + 437);
    const auto *nsh_438 = buffer.data(nsh + 438);
    const auto *nsh_440 = buffer.data(nsh + 440);
    const auto *nsh_441 = buffer.data(nsh + 441);
    const auto *nsh_444 = buffer.data(nsh + 444);

    const auto *nsi1_392 = buffer.data(nsi1 + 392);
    const auto *nsi1_395 = buffer.data(nsi1 + 395);
    const auto *nsi1_397 = buffer.data(nsi1 + 397);
    const auto *nsi1_398 = buffer.data(nsi1 + 398);
    const auto *nsi1_401 = buffer.data(nsi1 + 401);
    const auto *nsi1_402 = buffer.data(nsi1 + 402);
    const auto *nsi1_404 = buffer.data(nsi1 + 404);
    const auto *nsi1_406 = buffer.data(nsi1 + 406);
    const auto *nsi1_419 = buffer.data(nsi1 + 419);

    const auto *osg0_258 = buffer.data(osg0 + 258);
    const auto *osg0_260 = buffer.data(osg0 + 260);
    const auto *osg0_261 = buffer.data(osg0 + 261);
    const auto *osg0_264 = buffer.data(osg0 + 264);
    const auto *osg0_265 = buffer.data(osg0 + 265);
    const auto *osg0_267 = buffer.data(osg0 + 267);
    const auto *osg0_268 = buffer.data(osg0 + 268);
    const auto *osg0_269 = buffer.data(osg0 + 269);
    const auto *osg0_270 = buffer.data(osg0 + 270);
    const auto *osg0_273 = buffer.data(osg0 + 273);
    const auto *osg0_275 = buffer.data(osg0 + 275);
    const auto *osg0_276 = buffer.data(osg0 + 276);
    const auto *osg0_279 = buffer.data(osg0 + 279);
    const auto *osg0_280 = buffer.data(osg0 + 280);
    const auto *osg0_282 = buffer.data(osg0 + 282);
    const auto *osg0_283 = buffer.data(osg0 + 283);
    const auto *osg0_284 = buffer.data(osg0 + 284);
    const auto *osg0_295 = buffer.data(osg0 + 295);
    const auto *osg0_297 = buffer.data(osg0 + 297);
    const auto *osg0_298 = buffer.data(osg0 + 298);
    const auto *osg0_299 = buffer.data(osg0 + 299);
    const auto *osg0_300 = buffer.data(osg0 + 300);
    const auto *osg0_301 = buffer.data(osg0 + 301);
    const auto *osg0_302 = buffer.data(osg0 + 302);
    const auto *osg0_303 = buffer.data(osg0 + 303);
    const auto *osg0_304 = buffer.data(osg0 + 304);
    const auto *osg0_305 = buffer.data(osg0 + 305);
    const auto *osg0_309 = buffer.data(osg0 + 309);
    const auto *osg0_310 = buffer.data(osg0 + 310);
    const auto *osg0_311 = buffer.data(osg0 + 311);
    const auto *osg0_312 = buffer.data(osg0 + 312);
    const auto *osg0_313 = buffer.data(osg0 + 313);
    const auto *osg0_314 = buffer.data(osg0 + 314);
    const auto *osg0_315 = buffer.data(osg0 + 315);
    const auto *osg0_318 = buffer.data(osg0 + 318);

    const auto *osg1_258 = buffer.data(osg1 + 258);
    const auto *osg1_260 = buffer.data(osg1 + 260);
    const auto *osg1_261 = buffer.data(osg1 + 261);
    const auto *osg1_264 = buffer.data(osg1 + 264);
    const auto *osg1_265 = buffer.data(osg1 + 265);
    const auto *osg1_267 = buffer.data(osg1 + 267);
    const auto *osg1_268 = buffer.data(osg1 + 268);
    const auto *osg1_269 = buffer.data(osg1 + 269);
    const auto *osg1_270 = buffer.data(osg1 + 270);
    const auto *osg1_273 = buffer.data(osg1 + 273);
    const auto *osg1_275 = buffer.data(osg1 + 275);
    const auto *osg1_276 = buffer.data(osg1 + 276);
    const auto *osg1_279 = buffer.data(osg1 + 279);
    const auto *osg1_280 = buffer.data(osg1 + 280);
    const auto *osg1_282 = buffer.data(osg1 + 282);
    const auto *osg1_283 = buffer.data(osg1 + 283);
    const auto *osg1_284 = buffer.data(osg1 + 284);
    const auto *osg1_295 = buffer.data(osg1 + 295);
    const auto *osg1_297 = buffer.data(osg1 + 297);
    const auto *osg1_298 = buffer.data(osg1 + 298);
    const auto *osg1_299 = buffer.data(osg1 + 299);
    const auto *osg1_300 = buffer.data(osg1 + 300);
    const auto *osg1_301 = buffer.data(osg1 + 301);
    const auto *osg1_302 = buffer.data(osg1 + 302);
    const auto *osg1_303 = buffer.data(osg1 + 303);
    const auto *osg1_304 = buffer.data(osg1 + 304);
    const auto *osg1_305 = buffer.data(osg1 + 305);
    const auto *osg1_309 = buffer.data(osg1 + 309);
    const auto *osg1_310 = buffer.data(osg1 + 310);
    const auto *osg1_311 = buffer.data(osg1 + 311);
    const auto *osg1_312 = buffer.data(osg1 + 312);
    const auto *osg1_313 = buffer.data(osg1 + 313);
    const auto *osg1_314 = buffer.data(osg1 + 314);
    const auto *osg1_315 = buffer.data(osg1 + 315);
    const auto *osg1_318 = buffer.data(osg1 + 318);

    const auto *osh_359 = buffer.data(osh + 359);
    const auto *osh_360 = buffer.data(osh + 360);
    const auto *osh_362 = buffer.data(osh + 362);
    const auto *osh_363 = buffer.data(osh + 363);
    const auto *osh_366 = buffer.data(osh + 366);
    const auto *osh_367 = buffer.data(osh + 367);
    const auto *osh_369 = buffer.data(osh + 369);
    const auto *osh_371 = buffer.data(osh + 371);
    const auto *osh_372 = buffer.data(osh + 372);
    const auto *osh_373 = buffer.data(osh + 373);
    const auto *osh_374 = buffer.data(osh + 374);
    const auto *osh_375 = buffer.data(osh + 375);
    const auto *osh_376 = buffer.data(osh + 376);
    const auto *osh_377 = buffer.data(osh + 377);
    const auto *osh_378 = buffer.data(osh + 378);
    const auto *osh_380 = buffer.data(osh + 380);
    const auto *osh_381 = buffer.data(osh + 381);
    const auto *osh_383 = buffer.data(osh + 383);
    const auto *osh_384 = buffer.data(osh + 384);
    const auto *osh_387 = buffer.data(osh + 387);
    const auto *osh_388 = buffer.data(osh + 388);
    const auto *osh_390 = buffer.data(osh + 390);
    const auto *osh_392 = buffer.data(osh + 392);
    const auto *osh_393 = buffer.data(osh + 393);
    const auto *osh_394 = buffer.data(osh + 394);
    const auto *osh_395 = buffer.data(osh + 395);
    const auto *osh_396 = buffer.data(osh + 396);
    const auto *osh_397 = buffer.data(osh + 397);
    const auto *osh_398 = buffer.data(osh + 398);
    const auto *osh_399 = buffer.data(osh + 399);
    const auto *osh_401 = buffer.data(osh + 401);
    const auto *osh_402 = buffer.data(osh + 402);
    const auto *osh_404 = buffer.data(osh + 404);
    const auto *osh_405 = buffer.data(osh + 405);
    const auto *osh_408 = buffer.data(osh + 408);
    const auto *osh_414 = buffer.data(osh + 414);
    const auto *osh_415 = buffer.data(osh + 415);
    const auto *osh_416 = buffer.data(osh + 416);
    const auto *osh_417 = buffer.data(osh + 417);
    const auto *osh_418 = buffer.data(osh + 418);
    const auto *osh_419 = buffer.data(osh + 419);
    const auto *osh_420 = buffer.data(osh + 420);
    const auto *osh_421 = buffer.data(osh + 421);
    const auto *osh_422 = buffer.data(osh + 422);
    const auto *osh_423 = buffer.data(osh + 423);
    const auto *osh_424 = buffer.data(osh + 424);
    const auto *osh_425 = buffer.data(osh + 425);
    const auto *osh_426 = buffer.data(osh + 426);
    const auto *osh_427 = buffer.data(osh + 427);
    const auto *osh_428 = buffer.data(osh + 428);
    const auto *osh_429 = buffer.data(osh + 429);
    const auto *osh_434 = buffer.data(osh + 434);
    const auto *osh_435 = buffer.data(osh + 435);
    const auto *osh_436 = buffer.data(osh + 436);
    const auto *osh_437 = buffer.data(osh + 437);
    const auto *osh_438 = buffer.data(osh + 438);
    const auto *osh_439 = buffer.data(osh + 439);
    const auto *osh_440 = buffer.data(osh + 440);
    const auto *osh_441 = buffer.data(osh + 441);
    const auto *osh_444 = buffer.data(osh + 444);

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, nsh_254, nsh_360, nsh_362, osg0_258, \
                         osg0_260, osg1_258, osg1_260, osh_359, osh_360, \
                         osh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_21 * nsh_360[k]
                   + f_8 * osg0_258[k]
                   - f_9 * osg1_258[k]
                   + f_3 * pc_x[k] * osh_360[k];

        t_480[k] = f_13 * nsh_254[k]
                   + f_3 * pc_y[k] * osh_359[k];

        t_481[k] = f_21 * nsh_362[k]
                   + f_8 * osg0_260[k]
                   - f_9 * osg1_260[k]
                   + f_3 * pc_x[k] * osh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, nsh_234, nsh_257, nsh_363, \
                         osg0_261, osg1_261, osh_360, osh_362, \
                         osh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_21 * nsh_363[k]
                   + f_6 * osg0_261[k]
                   - f_7 * osg1_261[k]
                   + f_3 * pc_x[k] * osh_363[k];

        t_483[k] = f_12 * nsh_234[k]
                   + f_3 * pc_z[k] * osh_360[k];

        t_484[k] = f_13 * nsh_257[k]
                   + f_3 * pc_y[k] * osh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, nsh_237, nsh_366, nsh_367, osg0_264, \
                         osg0_265, osg1_264, osg1_265, osh_363, osh_366, \
                         osh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_21 * nsh_366[k]
                   + f_6 * osg0_264[k]
                   - f_7 * osg1_264[k]
                   + f_3 * pc_x[k] * osh_366[k];

        t_486[k] = f_21 * nsh_367[k]
                   + f_4 * osg0_265[k]
                   - f_5 * osg1_265[k]
                   + f_3 * pc_x[k] * osh_367[k];

        t_487[k] = f_12 * nsh_237[k]
                   + f_3 * pc_z[k] * osh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, nsh_261, nsh_369, nsh_371, osg0_267, \
                         osg0_269, osg1_267, osg1_269, osh_366, osh_369, \
                         osh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_21 * nsh_369[k]
                   + f_4 * osg0_267[k]
                   - f_5 * osg1_267[k]
                   + f_3 * pc_x[k] * osh_369[k];

        t_489[k] = f_13 * nsh_261[k]
                   + f_3 * pc_y[k] * osh_366[k];

        t_490[k] = f_21 * nsh_371[k]
                   + f_4 * osg0_269[k]
                   - f_5 * osg1_269[k]
                   + f_3 * pc_x[k] * osh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, nsh_372, nsh_373, nsh_374, \
                         nsh_375, nsh_376, osh_372, osh_373, osh_374, osh_375, \
                         osh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_21 * nsh_372[k]
                   + f_3 * pc_x[k] * osh_372[k];

        t_492[k] = f_21 * nsh_373[k]
                   + f_3 * pc_x[k] * osh_373[k];

        t_493[k] = f_21 * nsh_374[k]
                   + f_3 * pc_x[k] * osh_374[k];

        t_494[k] = f_21 * nsh_375[k]
                   + f_3 * pc_x[k] * osh_375[k];

        t_495[k] = f_21 * nsh_376[k]
                   + f_3 * pc_x[k] * osh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, nsh_246, nsh_267, nsh_377, \
                         osg0_265, osg1_265, osh_372, osh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_21 * nsh_377[k]
                   + f_3 * pc_x[k] * osh_377[k];

        t_497[k] = f_13 * nsh_267[k]
                   + f_1 * osg0_265[k]
                   - f_2 * osg1_265[k]
                   + f_3 * pc_y[k] * osh_372[k];

        t_498[k] = f_12 * nsh_246[k]
                   + f_3 * pc_z[k] * osh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, nsh_269, nsh_270, nsh_271, osg0_267, \
                         osg0_268, osg0_269, osg1_267, osg1_268, osg1_269, osh_374, osh_375, \
                         osh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * nsh_269[k]
                   + f_8 * osg0_267[k]
                   - f_9 * osg1_267[k]
                   + f_3 * pc_y[k] * osh_374[k];

        t_500[k] = f_13 * nsh_270[k]
                   + f_6 * osg0_268[k]
                   - f_7 * osg1_268[k]
                   + f_3 * pc_y[k] * osh_375[k];

        t_501[k] = f_13 * nsh_271[k]
                   + f_4 * osg0_269[k]
                   - f_5 * osg1_269[k]
                   + f_3 * pc_y[k] * osh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, nsh_251, nsh_272, nsh_378, \
                         osg0_269, osg0_270, osg1_269, osg1_270, osh_377, \
                         osh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * nsh_272[k]
                   + f_3 * pc_y[k] * osh_377[k];

        t_503[k] = f_12 * nsh_251[k]
                   + f_1 * osg0_269[k]
                   - f_2 * osg1_269[k]
                   + f_3 * pc_z[k] * osh_377[k];

        t_504[k] = f_21 * nsh_378[k]
                   + f_1 * osg0_270[k]
                   - f_2 * osg1_270[k]
                   + f_3 * pc_x[k] * osh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, nsh_252, nsh_273, \
                         nsh_275, nsh_381, osg0_273, osg1_273, osh_378, osh_380, \
                         osh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * nsh_273[k]
                   + f_3 * pc_y[k] * osh_378[k];

        t_506[k] = f_13 * nsh_252[k]
                   + f_3 * pc_z[k] * osh_378[k];

        t_507[k] = f_21 * nsh_381[k]
                   + f_8 * osg0_273[k]
                   - f_9 * osg1_273[k]
                   + f_3 * pc_x[k] * osh_381[k];

        t_508[k] = f_12 * nsh_275[k]
                   + f_3 * pc_y[k] * osh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, nsh_255, nsh_383, nsh_384, osg0_275, \
                         osg0_276, osg1_275, osg1_276, osh_381, osh_383, \
                         osh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_21 * nsh_383[k]
                   + f_8 * osg0_275[k]
                   - f_9 * osg1_275[k]
                   + f_3 * pc_x[k] * osh_383[k];

        t_510[k] = f_21 * nsh_384[k]
                   + f_6 * osg0_276[k]
                   - f_7 * osg1_276[k]
                   + f_3 * pc_x[k] * osh_384[k];

        t_511[k] = f_13 * nsh_255[k]
                   + f_3 * pc_z[k] * osh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, nsh_278, nsh_387, nsh_388, osg0_279, \
                         osg0_280, osg1_279, osg1_280, osh_383, osh_387, \
                         osh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * nsh_278[k]
                   + f_3 * pc_y[k] * osh_383[k];

        t_513[k] = f_21 * nsh_387[k]
                   + f_6 * osg0_279[k]
                   - f_7 * osg1_279[k]
                   + f_3 * pc_x[k] * osh_387[k];

        t_514[k] = f_21 * nsh_388[k]
                   + f_4 * osg0_280[k]
                   - f_5 * osg1_280[k]
                   + f_3 * pc_x[k] * osh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, nsh_258, nsh_282, nsh_390, \
                         osg0_282, osg1_282, osh_384, osh_387, \
                         osh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * nsh_258[k]
                   + f_3 * pc_z[k] * osh_384[k];

        t_516[k] = f_21 * nsh_390[k]
                   + f_4 * osg0_282[k]
                   - f_5 * osg1_282[k]
                   + f_3 * pc_x[k] * osh_390[k];

        t_517[k] = f_12 * nsh_282[k]
                   + f_3 * pc_y[k] * osh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, nsh_392, nsh_393, nsh_394, nsh_395, \
                         osg0_284, osg1_284, osh_392, osh_393, osh_394, \
                         osh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_21 * nsh_392[k]
                   + f_4 * osg0_284[k]
                   - f_5 * osg1_284[k]
                   + f_3 * pc_x[k] * osh_392[k];

        t_519[k] = f_21 * nsh_393[k]
                   + f_3 * pc_x[k] * osh_393[k];

        t_520[k] = f_21 * nsh_394[k]
                   + f_3 * pc_x[k] * osh_394[k];

        t_521[k] = f_21 * nsh_395[k]
                   + f_3 * pc_x[k] * osh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, nsh_288, nsh_396, nsh_397, \
                         nsh_398, osg0_280, osg1_280, osh_393, osh_396, osh_397, \
                         osh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_21 * nsh_396[k]
                   + f_3 * pc_x[k] * osh_396[k];

        t_523[k] = f_21 * nsh_397[k]
                   + f_3 * pc_x[k] * osh_397[k];

        t_524[k] = f_21 * nsh_398[k]
                   + f_3 * pc_x[k] * osh_398[k];

        t_525[k] = f_12 * nsh_288[k]
                   + f_1 * osg0_280[k]
                   - f_2 * osg1_280[k]
                   + f_3 * pc_y[k] * osh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, nsh_267, nsh_290, nsh_291, osg0_282, \
                         osg0_283, osg1_282, osg1_283, osh_393, osh_395, \
                         osh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * nsh_267[k]
                   + f_3 * pc_z[k] * osh_393[k];

        t_527[k] = f_12 * nsh_290[k]
                   + f_8 * osg0_282[k]
                   - f_9 * osg1_282[k]
                   + f_3 * pc_y[k] * osh_395[k];

        t_528[k] = f_12 * nsh_291[k]
                   + f_6 * osg0_283[k]
                   - f_7 * osg1_283[k]
                   + f_3 * pc_y[k] * osh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_y, pc_y, pc_z, nsi0_392, nsh_272, \
                         nsh_292, nsh_293, nsi1_392, osg0_284, osg1_284, osh_397, \
                         osh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * nsh_292[k]
                   + f_4 * osg0_284[k]
                   - f_5 * osg1_284[k]
                   + f_3 * pc_y[k] * osh_397[k];

        t_530[k] = f_12 * nsh_293[k]
                   + f_3 * pc_y[k] * osh_398[k];

        t_531[k] = f_13 * nsh_272[k]
                   + f_1 * osg0_284[k]
                   - f_2 * osg1_284[k]
                   + f_3 * pc_z[k] * osh_398[k];

        t_532[k] = pa_y[k] * nsi0_392[k]
                   - f_10 * pc_y[k] * nsi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pc_y, pc_z, nsi0_395, nsh_273, \
                         nsh_294, nsh_295, nsh_296, nsi1_395, osh_399, \
                         osh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * nsh_294[k]
                   + f_3 * pc_y[k] * osh_399[k];

        t_534[k] = f_14 * nsh_273[k]
                   + f_3 * pc_z[k] * osh_399[k];

        t_535[k] = pa_y[k] * nsi0_395[k]
                   + f_12 * nsh_295[k]
                   - f_10 * pc_y[k] * nsi1_395[k];

        t_536[k] = f_11 * nsh_296[k]
                   + f_3 * pc_y[k] * osh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_y, pc_y, pc_z, nsi0_397, nsi0_398, \
                         nsh_276, nsh_297, nsh_299, nsi1_397, nsi1_398, osh_402, \
                         osh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * nsi0_397[k]
                   - f_10 * pc_y[k] * nsi1_397[k];

        t_538[k] = pa_y[k] * nsi0_398[k]
                   + f_13 * nsh_297[k]
                   - f_10 * pc_y[k] * nsi1_398[k];

        t_539[k] = f_14 * nsh_276[k]
                   + f_3 * pc_z[k] * osh_402[k];

        t_540[k] = f_11 * nsh_299[k]
                   + f_3 * pc_y[k] * osh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_y, pc_y, pc_z, nsi0_401, nsi0_402, nsh_279, \
                         nsh_300, nsi1_401, nsi1_402, osh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * nsi0_401[k]
                   - f_10 * pc_y[k] * nsi1_401[k];

        t_542[k] = pa_y[k] * nsi0_402[k]
                   + f_14 * nsh_300[k]
                   - f_10 * pc_y[k] * nsi1_402[k];

        t_543[k] = f_14 * nsh_279[k]
                   + f_3 * pc_z[k] * osh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_x, pc_y, nsi0_404, nsi0_406, \
                         nsh_302, nsh_303, nsh_414, nsi1_404, nsi1_406, osh_408, \
                         osh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pa_y[k] * nsi0_404[k]
                   + f_12 * nsh_302[k]
                   - f_10 * pc_y[k] * nsi1_404[k];

        t_545[k] = f_11 * nsh_303[k]
                   + f_3 * pc_y[k] * osh_408[k];

        t_546[k] = pa_y[k] * nsi0_406[k]
                   - f_10 * pc_y[k] * nsi1_406[k];

        t_547[k] = f_21 * nsh_414[k]
                   + f_3 * pc_x[k] * osh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, nsh_415, nsh_416, nsh_417, \
                         nsh_418, nsh_419, osh_415, osh_416, osh_417, osh_418, \
                         osh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_21 * nsh_415[k]
                   + f_3 * pc_x[k] * osh_415[k];

        t_549[k] = f_21 * nsh_416[k]
                   + f_3 * pc_x[k] * osh_416[k];

        t_550[k] = f_21 * nsh_417[k]
                   + f_3 * pc_x[k] * osh_417[k];

        t_551[k] = f_21 * nsh_418[k]
                   + f_3 * pc_x[k] * osh_418[k];

        t_552[k] = f_21 * nsh_419[k]
                   + f_3 * pc_x[k] * osh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, nsh_288, nsh_309, nsh_311, osg0_295, \
                         osg0_297, osg1_295, osg1_297, osh_414, \
                         osh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * nsh_309[k]
                   + f_1 * osg0_295[k]
                   - f_2 * osg1_295[k]
                   + f_3 * pc_y[k] * osh_414[k];

        t_554[k] = f_14 * nsh_288[k]
                   + f_3 * pc_z[k] * osh_414[k];

        t_555[k] = f_11 * nsh_311[k]
                   + f_8 * osg0_297[k]
                   - f_9 * osg1_297[k]
                   + f_3 * pc_y[k] * osh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, nsh_312, nsh_313, nsh_314, osg0_298, \
                         osg0_299, osg1_298, osg1_299, osh_417, osh_418, \
                         osh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * nsh_312[k]
                   + f_6 * osg0_298[k]
                   - f_7 * osg1_298[k]
                   + f_3 * pc_y[k] * osh_417[k];

        t_557[k] = f_11 * nsh_313[k]
                   + f_4 * osg0_299[k]
                   - f_5 * osg1_299[k]
                   + f_3 * pc_y[k] * osh_418[k];

        t_558[k] = f_11 * nsh_314[k]
                   + f_3 * pc_y[k] * osh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, nsi0_419, \
                         nsh_294, nsh_420, nsi1_419, osg0_300, osg1_300, \
                         osh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_y[k] * nsi0_419[k]
                   - f_10 * pc_y[k] * nsi1_419[k];

        t_560[k] = f_21 * nsh_420[k]
                   + f_1 * osg0_300[k]
                   - f_2 * osg1_300[k]
                   + f_3 * pc_x[k] * osh_420[k];

        t_561[k] = f_3 * pc_y[k] * osh_420[k];

        t_562[k] = f_22 * nsh_294[k]
                   + f_3 * pc_z[k] * osh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, nsh_425, osg0_300, osg0_305, \
                         osg1_300, osg1_305, osh_421, osh_422, \
                         osh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_4 * osg0_300[k]
                   - f_5 * osg1_300[k]
                   + f_3 * pc_y[k] * osh_421[k];

        t_564[k] = f_3 * pc_y[k] * osh_422[k];

        t_565[k] = f_21 * nsh_425[k]
                   + f_8 * osg0_305[k]
                   - f_9 * osg1_305[k]
                   + f_3 * pc_x[k] * osh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_y, osg0_301, osg0_302, osg1_301, osg1_302, \
                         osh_423, osh_424, osh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_6 * osg0_301[k]
                   - f_7 * osg1_301[k]
                   + f_3 * pc_y[k] * osh_423[k];

        t_567[k] = f_4 * osg0_302[k]
                   - f_5 * osg1_302[k]
                   + f_3 * pc_y[k] * osh_424[k];

        t_568[k] = f_3 * pc_y[k] * osh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_y, nsh_429, osg0_303, osg0_304, \
                         osg0_309, osg1_303, osg1_304, osg1_309, osh_426, osh_427, \
                         osh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_21 * nsh_429[k]
                   + f_6 * osg0_309[k]
                   - f_7 * osg1_309[k]
                   + f_3 * pc_x[k] * osh_429[k];

        t_570[k] = f_8 * osg0_303[k]
                   - f_9 * osg1_303[k]
                   + f_3 * pc_y[k] * osh_426[k];

        t_571[k] = f_6 * osg0_304[k]
                   - f_7 * osg1_304[k]
                   + f_3 * pc_y[k] * osh_427[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, nsh_434, nsh_435, osg0_305, \
                         osg0_314, osg1_305, osg1_314, osh_428, osh_429, osh_434, \
                         osh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * osg0_305[k]
                   - f_5 * osg1_305[k]
                   + f_3 * pc_y[k] * osh_428[k];

        t_573[k] = f_3 * pc_y[k] * osh_429[k];

        t_574[k] = f_21 * nsh_434[k]
                   + f_4 * osg0_314[k]
                   - f_5 * osg1_314[k]
                   + f_3 * pc_x[k] * osh_434[k];

        t_575[k] = f_21 * nsh_435[k]
                   + f_3 * pc_x[k] * osh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pc_x, pc_y, nsh_436, nsh_437, \
                         nsh_438, nsh_440, osh_434, osh_436, osh_437, osh_438, \
                         osh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_21 * nsh_436[k]
                   + f_3 * pc_x[k] * osh_436[k];

        t_577[k] = f_21 * nsh_437[k]
                   + f_3 * pc_x[k] * osh_437[k];

        t_578[k] = f_21 * nsh_438[k]
                   + f_3 * pc_x[k] * osh_438[k];

        t_579[k] = f_3 * pc_y[k] * osh_434[k];

        t_580[k] = f_21 * nsh_440[k]
                   + f_3 * pc_x[k] * osh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, osg0_310, osg0_311, osg0_312, osg1_310, \
                         osg1_311, osg1_312, osh_435, osh_436, \
                         osh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * osg0_310[k]
                   - f_2 * osg1_310[k]
                   + f_3 * pc_y[k] * osh_435[k];

        t_582[k] = f_16 * osg0_311[k]
                   - f_17 * osg1_311[k]
                   + f_3 * pc_y[k] * osh_436[k];

        t_583[k] = f_8 * osg0_312[k]
                   - f_9 * osg1_312[k]
                   + f_3 * pc_y[k] * osh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, nsh_314, osg0_313, osg0_314, \
                         osg1_313, osg1_314, osh_438, osh_439, \
                         osh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * osg0_313[k]
                   - f_7 * osg1_313[k]
                   + f_3 * pc_y[k] * osh_438[k];

        t_585[k] = f_4 * osg0_314[k]
                   - f_5 * osg1_314[k]
                   + f_3 * pc_y[k] * osh_439[k];

        t_586[k] = f_3 * pc_y[k] * osh_440[k];

        t_587[k] = f_22 * nsh_314[k]
                   + f_1 * osg0_314[k]
                   - f_2 * osg1_314[k]
                   + f_3 * pc_z[k] * osh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, nsh_315, nsh_441, \
                         nsh_444, osg0_315, osg0_318, osg1_315, osg1_318, osh_441, \
                         osh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_22 * nsh_441[k]
                   + f_1 * osg0_315[k]
                   - f_2 * osg1_315[k]
                   + f_3 * pc_x[k] * osh_441[k];

        t_589[k] = f_21 * nsh_315[k]
                   + f_3 * pc_y[k] * osh_441[k];

        t_590[k] = f_3 * pc_z[k] * osh_441[k];

        t_591[k] = f_22 * nsh_444[k]
                   + f_8 * osg0_318[k]
                   - f_9 * osg1_318[k]
                   + f_3 * pc_x[k] * osh_444[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_420 = buffer.data(nsi0 + 420);
    const auto *nsi0_423 = buffer.data(nsi0 + 423);
    const auto *nsi0_426 = buffer.data(nsi0 + 426);
    const auto *nsi0_430 = buffer.data(nsi0 + 430);
    const auto *nsi0_432 = buffer.data(nsi0 + 432);
    const auto *nsi0_441 = buffer.data(nsi0 + 441);

    const auto *nsh_315 = buffer.data(nsh + 315);
    const auto *nsh_318 = buffer.data(nsh + 318);
    const auto *nsh_320 = buffer.data(nsh + 320);
    const auto *nsh_321 = buffer.data(nsh + 321);
    const auto *nsh_322 = buffer.data(nsh + 322);
    const auto *nsh_324 = buffer.data(nsh + 324);
    const auto *nsh_330 = buffer.data(nsh + 330);
    const auto *nsh_335 = buffer.data(nsh + 335);
    const auto *nsh_336 = buffer.data(nsh + 336);
    const auto *nsh_338 = buffer.data(nsh + 338);
    const auto *nsh_339 = buffer.data(nsh + 339);
    const auto *nsh_341 = buffer.data(nsh + 341);
    const auto *nsh_342 = buffer.data(nsh + 342);
    const auto *nsh_345 = buffer.data(nsh + 345);
    const auto *nsh_351 = buffer.data(nsh + 351);
    const auto *nsh_353 = buffer.data(nsh + 353);
    const auto *nsh_354 = buffer.data(nsh + 354);
    const auto *nsh_355 = buffer.data(nsh + 355);
    const auto *nsh_356 = buffer.data(nsh + 356);
    const auto *nsh_357 = buffer.data(nsh + 357);
    const auto *nsh_359 = buffer.data(nsh + 359);
    const auto *nsh_360 = buffer.data(nsh + 360);
    const auto *nsh_362 = buffer.data(nsh + 362);
    const auto *nsh_363 = buffer.data(nsh + 363);
    const auto *nsh_366 = buffer.data(nsh + 366);
    const auto *nsh_372 = buffer.data(nsh + 372);
    const auto *nsh_374 = buffer.data(nsh + 374);
    const auto *nsh_375 = buffer.data(nsh + 375);
    const auto *nsh_376 = buffer.data(nsh + 376);
    const auto *nsh_377 = buffer.data(nsh + 377);
    const auto *nsh_378 = buffer.data(nsh + 378);
    const auto *nsh_380 = buffer.data(nsh + 380);
    const auto *nsh_383 = buffer.data(nsh + 383);
    const auto *nsh_387 = buffer.data(nsh + 387);
    const auto *nsh_393 = buffer.data(nsh + 393);
    const auto *nsh_395 = buffer.data(nsh + 395);
    const auto *nsh_396 = buffer.data(nsh + 396);
    const auto *nsh_397 = buffer.data(nsh + 397);
    const auto *nsh_398 = buffer.data(nsh + 398);
    const auto *nsh_399 = buffer.data(nsh + 399);
    const auto *nsh_447 = buffer.data(nsh + 447);
    const auto *nsh_451 = buffer.data(nsh + 451);
    const auto *nsh_456 = buffer.data(nsh + 456);
    const auto *nsh_458 = buffer.data(nsh + 458);
    const auto *nsh_459 = buffer.data(nsh + 459);
    const auto *nsh_460 = buffer.data(nsh + 460);
    const auto *nsh_461 = buffer.data(nsh + 461);
    const auto *nsh_467 = buffer.data(nsh + 467);
    const auto *nsh_471 = buffer.data(nsh + 471);
    const auto *nsh_476 = buffer.data(nsh + 476);
    const auto *nsh_477 = buffer.data(nsh + 477);
    const auto *nsh_478 = buffer.data(nsh + 478);
    const auto *nsh_479 = buffer.data(nsh + 479);
    const auto *nsh_480 = buffer.data(nsh + 480);
    const auto *nsh_481 = buffer.data(nsh + 481);
    const auto *nsh_482 = buffer.data(nsh + 482);
    const auto *nsh_483 = buffer.data(nsh + 483);
    const auto *nsh_486 = buffer.data(nsh + 486);
    const auto *nsh_488 = buffer.data(nsh + 488);
    const auto *nsh_489 = buffer.data(nsh + 489);
    const auto *nsh_492 = buffer.data(nsh + 492);
    const auto *nsh_493 = buffer.data(nsh + 493);
    const auto *nsh_495 = buffer.data(nsh + 495);
    const auto *nsh_497 = buffer.data(nsh + 497);
    const auto *nsh_498 = buffer.data(nsh + 498);
    const auto *nsh_499 = buffer.data(nsh + 499);
    const auto *nsh_500 = buffer.data(nsh + 500);
    const auto *nsh_501 = buffer.data(nsh + 501);
    const auto *nsh_502 = buffer.data(nsh + 502);
    const auto *nsh_503 = buffer.data(nsh + 503);
    const auto *nsh_504 = buffer.data(nsh + 504);
    const auto *nsh_507 = buffer.data(nsh + 507);
    const auto *nsh_509 = buffer.data(nsh + 509);
    const auto *nsh_510 = buffer.data(nsh + 510);
    const auto *nsh_513 = buffer.data(nsh + 513);
    const auto *nsh_514 = buffer.data(nsh + 514);
    const auto *nsh_516 = buffer.data(nsh + 516);
    const auto *nsh_518 = buffer.data(nsh + 518);
    const auto *nsh_519 = buffer.data(nsh + 519);
    const auto *nsh_520 = buffer.data(nsh + 520);
    const auto *nsh_521 = buffer.data(nsh + 521);
    const auto *nsh_522 = buffer.data(nsh + 522);
    const auto *nsh_523 = buffer.data(nsh + 523);
    const auto *nsh_524 = buffer.data(nsh + 524);
    const auto *nsh_525 = buffer.data(nsh + 525);

    const auto *nsi1_420 = buffer.data(nsi1 + 420);
    const auto *nsi1_423 = buffer.data(nsi1 + 423);
    const auto *nsi1_426 = buffer.data(nsi1 + 426);
    const auto *nsi1_430 = buffer.data(nsi1 + 430);
    const auto *nsi1_432 = buffer.data(nsi1 + 432);
    const auto *nsi1_441 = buffer.data(nsi1 + 441);

    const auto *osg0_315 = buffer.data(osg0 + 315);
    const auto *osg0_317 = buffer.data(osg0 + 317);
    const auto *osg0_318 = buffer.data(osg0 + 318);
    const auto *osg0_320 = buffer.data(osg0 + 320);
    const auto *osg0_321 = buffer.data(osg0 + 321);
    const auto *osg0_325 = buffer.data(osg0 + 325);
    const auto *osg0_326 = buffer.data(osg0 + 326);
    const auto *osg0_327 = buffer.data(osg0 + 327);
    const auto *osg0_329 = buffer.data(osg0 + 329);
    const auto *osg0_335 = buffer.data(osg0 + 335);
    const auto *osg0_339 = buffer.data(osg0 + 339);
    const auto *osg0_342 = buffer.data(osg0 + 342);
    const auto *osg0_343 = buffer.data(osg0 + 343);
    const auto *osg0_344 = buffer.data(osg0 + 344);
    const auto *osg0_345 = buffer.data(osg0 + 345);
    const auto *osg0_348 = buffer.data(osg0 + 348);
    const auto *osg0_350 = buffer.data(osg0 + 350);
    const auto *osg0_351 = buffer.data(osg0 + 351);
    const auto *osg0_354 = buffer.data(osg0 + 354);
    const auto *osg0_355 = buffer.data(osg0 + 355);
    const auto *osg0_357 = buffer.data(osg0 + 357);
    const auto *osg0_358 = buffer.data(osg0 + 358);
    const auto *osg0_359 = buffer.data(osg0 + 359);
    const auto *osg0_360 = buffer.data(osg0 + 360);
    const auto *osg0_363 = buffer.data(osg0 + 363);
    const auto *osg0_365 = buffer.data(osg0 + 365);
    const auto *osg0_366 = buffer.data(osg0 + 366);
    const auto *osg0_369 = buffer.data(osg0 + 369);
    const auto *osg0_370 = buffer.data(osg0 + 370);
    const auto *osg0_372 = buffer.data(osg0 + 372);
    const auto *osg0_373 = buffer.data(osg0 + 373);
    const auto *osg0_374 = buffer.data(osg0 + 374);
    const auto *osg0_375 = buffer.data(osg0 + 375);

    const auto *osg1_315 = buffer.data(osg1 + 315);
    const auto *osg1_317 = buffer.data(osg1 + 317);
    const auto *osg1_318 = buffer.data(osg1 + 318);
    const auto *osg1_320 = buffer.data(osg1 + 320);
    const auto *osg1_321 = buffer.data(osg1 + 321);
    const auto *osg1_325 = buffer.data(osg1 + 325);
    const auto *osg1_326 = buffer.data(osg1 + 326);
    const auto *osg1_327 = buffer.data(osg1 + 327);
    const auto *osg1_329 = buffer.data(osg1 + 329);
    const auto *osg1_335 = buffer.data(osg1 + 335);
    const auto *osg1_339 = buffer.data(osg1 + 339);
    const auto *osg1_342 = buffer.data(osg1 + 342);
    const auto *osg1_343 = buffer.data(osg1 + 343);
    const auto *osg1_344 = buffer.data(osg1 + 344);
    const auto *osg1_345 = buffer.data(osg1 + 345);
    const auto *osg1_348 = buffer.data(osg1 + 348);
    const auto *osg1_350 = buffer.data(osg1 + 350);
    const auto *osg1_351 = buffer.data(osg1 + 351);
    const auto *osg1_354 = buffer.data(osg1 + 354);
    const auto *osg1_355 = buffer.data(osg1 + 355);
    const auto *osg1_357 = buffer.data(osg1 + 357);
    const auto *osg1_358 = buffer.data(osg1 + 358);
    const auto *osg1_359 = buffer.data(osg1 + 359);
    const auto *osg1_360 = buffer.data(osg1 + 360);
    const auto *osg1_363 = buffer.data(osg1 + 363);
    const auto *osg1_365 = buffer.data(osg1 + 365);
    const auto *osg1_366 = buffer.data(osg1 + 366);
    const auto *osg1_369 = buffer.data(osg1 + 369);
    const auto *osg1_370 = buffer.data(osg1 + 370);
    const auto *osg1_372 = buffer.data(osg1 + 372);
    const auto *osg1_373 = buffer.data(osg1 + 373);
    const auto *osg1_374 = buffer.data(osg1 + 374);
    const auto *osg1_375 = buffer.data(osg1 + 375);

    const auto *osh_442 = buffer.data(osh + 442);
    const auto *osh_443 = buffer.data(osh + 443);
    const auto *osh_444 = buffer.data(osh + 444);
    const auto *osh_446 = buffer.data(osh + 446);
    const auto *osh_447 = buffer.data(osh + 447);
    const auto *osh_448 = buffer.data(osh + 448);
    const auto *osh_450 = buffer.data(osh + 450);
    const auto *osh_451 = buffer.data(osh + 451);
    const auto *osh_456 = buffer.data(osh + 456);
    const auto *osh_457 = buffer.data(osh + 457);
    const auto *osh_458 = buffer.data(osh + 458);
    const auto *osh_459 = buffer.data(osh + 459);
    const auto *osh_460 = buffer.data(osh + 460);
    const auto *osh_461 = buffer.data(osh + 461);
    const auto *osh_462 = buffer.data(osh + 462);
    const auto *osh_464 = buffer.data(osh + 464);
    const auto *osh_465 = buffer.data(osh + 465);
    const auto *osh_467 = buffer.data(osh + 467);
    const auto *osh_468 = buffer.data(osh + 468);
    const auto *osh_471 = buffer.data(osh + 471);
    const auto *osh_476 = buffer.data(osh + 476);
    const auto *osh_477 = buffer.data(osh + 477);
    const auto *osh_478 = buffer.data(osh + 478);
    const auto *osh_479 = buffer.data(osh + 479);
    const auto *osh_480 = buffer.data(osh + 480);
    const auto *osh_481 = buffer.data(osh + 481);
    const auto *osh_482 = buffer.data(osh + 482);
    const auto *osh_483 = buffer.data(osh + 483);
    const auto *osh_485 = buffer.data(osh + 485);
    const auto *osh_486 = buffer.data(osh + 486);
    const auto *osh_488 = buffer.data(osh + 488);
    const auto *osh_489 = buffer.data(osh + 489);
    const auto *osh_492 = buffer.data(osh + 492);
    const auto *osh_493 = buffer.data(osh + 493);
    const auto *osh_495 = buffer.data(osh + 495);
    const auto *osh_497 = buffer.data(osh + 497);
    const auto *osh_498 = buffer.data(osh + 498);
    const auto *osh_499 = buffer.data(osh + 499);
    const auto *osh_500 = buffer.data(osh + 500);
    const auto *osh_501 = buffer.data(osh + 501);
    const auto *osh_502 = buffer.data(osh + 502);
    const auto *osh_503 = buffer.data(osh + 503);
    const auto *osh_504 = buffer.data(osh + 504);
    const auto *osh_506 = buffer.data(osh + 506);
    const auto *osh_507 = buffer.data(osh + 507);
    const auto *osh_509 = buffer.data(osh + 509);
    const auto *osh_510 = buffer.data(osh + 510);
    const auto *osh_513 = buffer.data(osh + 513);
    const auto *osh_514 = buffer.data(osh + 514);
    const auto *osh_516 = buffer.data(osh + 516);
    const auto *osh_518 = buffer.data(osh + 518);
    const auto *osh_519 = buffer.data(osh + 519);
    const auto *osh_520 = buffer.data(osh + 520);
    const auto *osh_521 = buffer.data(osh + 521);
    const auto *osh_522 = buffer.data(osh + 522);
    const auto *osh_523 = buffer.data(osh + 523);
    const auto *osh_524 = buffer.data(osh + 524);
    const auto *osh_525 = buffer.data(osh + 525);

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, nsh_447, osg0_315, osg0_321, \
                         osg1_315, osg1_321, osh_442, osh_443, osh_444, \
                         osh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * osh_442[k];

        t_593[k] = f_4 * osg0_315[k]
                   - f_5 * osg1_315[k]
                   + f_3 * pc_z[k] * osh_443[k];

        t_594[k] = f_22 * nsh_447[k]
                   + f_6 * osg0_321[k]
                   - f_7 * osg1_321[k]
                   + f_3 * pc_x[k] * osh_447[k];

        t_595[k] = f_3 * pc_z[k] * osh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, nsh_320, nsh_451, \
                         osg0_317, osg0_325, osg1_317, osg1_325, osh_446, osh_447, \
                         osh_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_21 * nsh_320[k]
                   + f_3 * pc_y[k] * osh_446[k];

        t_597[k] = f_6 * osg0_317[k]
                   - f_7 * osg1_317[k]
                   + f_3 * pc_z[k] * osh_446[k];

        t_598[k] = f_22 * nsh_451[k]
                   + f_4 * osg0_325[k]
                   - f_5 * osg1_325[k]
                   + f_3 * pc_x[k] * osh_451[k];

        t_599[k] = f_3 * pc_z[k] * osh_447[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, nsh_324, nsh_456, \
                         osg0_318, osg0_320, osg1_318, osg1_320, osh_448, osh_450, \
                         osh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_4 * osg0_318[k]
                   - f_5 * osg1_318[k]
                   + f_3 * pc_z[k] * osh_448[k];

        t_601[k] = f_21 * nsh_324[k]
                   + f_3 * pc_y[k] * osh_450[k];

        t_602[k] = f_8 * osg0_320[k]
                   - f_9 * osg1_320[k]
                   + f_3 * pc_z[k] * osh_450[k];

        t_603[k] = f_22 * nsh_456[k]
                   + f_3 * pc_x[k] * osh_456[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_z, nsh_458, nsh_459, \
                         nsh_460, nsh_461, osh_451, osh_458, osh_459, osh_460, \
                         osh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * osh_451[k];

        t_605[k] = f_22 * nsh_458[k]
                   + f_3 * pc_x[k] * osh_458[k];

        t_606[k] = f_22 * nsh_459[k]
                   + f_3 * pc_x[k] * osh_459[k];

        t_607[k] = f_22 * nsh_460[k]
                   + f_3 * pc_x[k] * osh_460[k];

        t_608[k] = f_22 * nsh_461[k]
                   + f_3 * pc_x[k] * osh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pc_y, pc_z, nsh_330, osg0_325, osg0_326, \
                         osg1_325, osg1_326, osh_456, osh_457, \
                         osh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_21 * nsh_330[k]
                   + f_1 * osg0_325[k]
                   - f_2 * osg1_325[k]
                   + f_3 * pc_y[k] * osh_456[k];

        t_610[k] = f_3 * pc_z[k] * osh_456[k];

        t_611[k] = f_4 * osg0_325[k]
                   - f_5 * osg1_325[k]
                   + f_3 * pc_z[k] * osh_457[k];

        t_612[k] = f_6 * osg0_326[k]
                   - f_7 * osg1_326[k]
                   + f_3 * pc_z[k] * osh_458[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pc_y, pc_z, nsi0_420, nsh_335, \
                         nsi1_420, osg0_327, osg0_329, osg1_327, osg1_329, osh_459, \
                         osh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_8 * osg0_327[k]
                   - f_9 * osg1_327[k]
                   + f_3 * pc_z[k] * osh_459[k];

        t_614[k] = f_21 * nsh_335[k]
                   + f_3 * pc_y[k] * osh_461[k];

        t_615[k] = f_1 * osg0_329[k]
                   - f_2 * osg1_329[k]
                   + f_3 * pc_z[k] * osh_461[k];

        t_616[k] = pa_z[k] * nsi0_420[k]
                   - f_10 * pc_z[k] * nsi1_420[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pc_y, pc_z, nsi0_423, nsh_315, \
                         nsh_336, nsh_338, nsi1_423, osh_462, osh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_22 * nsh_336[k]
                   + f_3 * pc_y[k] * osh_462[k];

        t_618[k] = f_11 * nsh_315[k]
                   + f_3 * pc_z[k] * osh_462[k];

        t_619[k] = pa_z[k] * nsi0_423[k]
                   - f_10 * pc_z[k] * nsi1_423[k];

        t_620[k] = f_22 * nsh_338[k]
                   + f_3 * pc_y[k] * osh_464[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_z, pc_x, pc_z, nsi0_426, nsh_318, nsh_467, \
                         nsi1_426, osg0_335, osg1_335, osh_465, \
                         osh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_22 * nsh_467[k]
                   + f_8 * osg0_335[k]
                   - f_9 * osg1_335[k]
                   + f_3 * pc_x[k] * osh_467[k];

        t_622[k] = pa_z[k] * nsi0_426[k]
                   - f_10 * pc_z[k] * nsi1_426[k];

        t_623[k] = f_11 * nsh_318[k]
                   + f_3 * pc_z[k] * osh_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pc_x, pc_y, pc_z, nsi0_430, nsh_341, \
                         nsh_471, nsi1_430, osg0_339, osg1_339, osh_467, \
                         osh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_22 * nsh_341[k]
                   + f_3 * pc_y[k] * osh_467[k];

        t_625[k] = f_22 * nsh_471[k]
                   + f_6 * osg0_339[k]
                   - f_7 * osg1_339[k]
                   + f_3 * pc_x[k] * osh_471[k];

        t_626[k] = pa_z[k] * nsi0_430[k]
                   - f_10 * pc_z[k] * nsi1_430[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_z, pc_y, pc_z, nsi0_432, nsh_321, nsh_322, \
                         nsh_345, nsi1_432, osh_468, osh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_11 * nsh_321[k]
                   + f_3 * pc_z[k] * osh_468[k];

        t_628[k] = pa_z[k] * nsi0_432[k]
                   + f_12 * nsh_322[k]
                   - f_10 * pc_z[k] * nsi1_432[k];

        t_629[k] = f_22 * nsh_345[k]
                   + f_3 * pc_y[k] * osh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, nsh_476, nsh_477, nsh_478, nsh_479, \
                         osg0_344, osg1_344, osh_476, osh_477, osh_478, \
                         osh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_22 * nsh_476[k]
                   + f_4 * osg0_344[k]
                   - f_5 * osg1_344[k]
                   + f_3 * pc_x[k] * osh_476[k];

        t_631[k] = f_22 * nsh_477[k]
                   + f_3 * pc_x[k] * osh_477[k];

        t_632[k] = f_22 * nsh_478[k]
                   + f_3 * pc_x[k] * osh_478[k];

        t_633[k] = f_22 * nsh_479[k]
                   + f_3 * pc_x[k] * osh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_z, pc_x, pc_z, nsi0_441, nsh_480, \
                         nsh_481, nsh_482, nsi1_441, osh_480, osh_481, \
                         osh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_22 * nsh_480[k]
                   + f_3 * pc_x[k] * osh_480[k];

        t_635[k] = f_22 * nsh_481[k]
                   + f_3 * pc_x[k] * osh_481[k];

        t_636[k] = f_22 * nsh_482[k]
                   + f_3 * pc_x[k] * osh_482[k];

        t_637[k] = pa_z[k] * nsi0_441[k]
                   - f_10 * pc_z[k] * nsi1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, nsh_330, nsh_353, nsh_354, osg0_342, \
                         osg0_343, osg1_342, osg1_343, osh_477, osh_479, \
                         osh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * nsh_330[k]
                   + f_3 * pc_z[k] * osh_477[k];

        t_639[k] = f_22 * nsh_353[k]
                   + f_8 * osg0_342[k]
                   - f_9 * osg1_342[k]
                   + f_3 * pc_y[k] * osh_479[k];

        t_640[k] = f_22 * nsh_354[k]
                   + f_6 * osg0_343[k]
                   - f_7 * osg1_343[k]
                   + f_3 * pc_y[k] * osh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, nsh_335, nsh_355, nsh_356, osg0_344, \
                         osg1_344, osh_481, osh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_22 * nsh_355[k]
                   + f_4 * osg0_344[k]
                   - f_5 * osg1_344[k]
                   + f_3 * pc_y[k] * osh_481[k];

        t_642[k] = f_22 * nsh_356[k]
                   + f_3 * pc_y[k] * osh_482[k];

        t_643[k] = f_11 * nsh_335[k]
                   + f_1 * osg0_344[k]
                   - f_2 * osg1_344[k]
                   + f_3 * pc_z[k] * osh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, nsh_336, nsh_357, nsh_483, \
                         osg0_345, osg1_345, osh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_22 * nsh_483[k]
                   + f_1 * osg0_345[k]
                   - f_2 * osg1_345[k]
                   + f_3 * pc_x[k] * osh_483[k];

        t_645[k] = f_14 * nsh_357[k]
                   + f_3 * pc_y[k] * osh_483[k];

        t_646[k] = f_12 * nsh_336[k]
                   + f_3 * pc_z[k] * osh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, nsh_359, nsh_486, nsh_488, osg0_348, \
                         osg0_350, osg1_348, osg1_350, osh_485, osh_486, \
                         osh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_22 * nsh_486[k]
                   + f_8 * osg0_348[k]
                   - f_9 * osg1_348[k]
                   + f_3 * pc_x[k] * osh_486[k];

        t_648[k] = f_14 * nsh_359[k]
                   + f_3 * pc_y[k] * osh_485[k];

        t_649[k] = f_22 * nsh_488[k]
                   + f_8 * osg0_350[k]
                   - f_9 * osg1_350[k]
                   + f_3 * pc_x[k] * osh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, nsh_339, nsh_362, nsh_489, \
                         osg0_351, osg1_351, osh_486, osh_488, \
                         osh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_22 * nsh_489[k]
                   + f_6 * osg0_351[k]
                   - f_7 * osg1_351[k]
                   + f_3 * pc_x[k] * osh_489[k];

        t_651[k] = f_12 * nsh_339[k]
                   + f_3 * pc_z[k] * osh_486[k];

        t_652[k] = f_14 * nsh_362[k]
                   + f_3 * pc_y[k] * osh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, nsh_342, nsh_492, nsh_493, osg0_354, \
                         osg0_355, osg1_354, osg1_355, osh_489, osh_492, \
                         osh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_22 * nsh_492[k]
                   + f_6 * osg0_354[k]
                   - f_7 * osg1_354[k]
                   + f_3 * pc_x[k] * osh_492[k];

        t_654[k] = f_22 * nsh_493[k]
                   + f_4 * osg0_355[k]
                   - f_5 * osg1_355[k]
                   + f_3 * pc_x[k] * osh_493[k];

        t_655[k] = f_12 * nsh_342[k]
                   + f_3 * pc_z[k] * osh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, nsh_366, nsh_495, nsh_497, osg0_357, \
                         osg0_359, osg1_357, osg1_359, osh_492, osh_495, \
                         osh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_22 * nsh_495[k]
                   + f_4 * osg0_357[k]
                   - f_5 * osg1_357[k]
                   + f_3 * pc_x[k] * osh_495[k];

        t_657[k] = f_14 * nsh_366[k]
                   + f_3 * pc_y[k] * osh_492[k];

        t_658[k] = f_22 * nsh_497[k]
                   + f_4 * osg0_359[k]
                   - f_5 * osg1_359[k]
                   + f_3 * pc_x[k] * osh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, nsh_498, nsh_499, nsh_500, \
                         nsh_501, nsh_502, osh_498, osh_499, osh_500, osh_501, \
                         osh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_22 * nsh_498[k]
                   + f_3 * pc_x[k] * osh_498[k];

        t_660[k] = f_22 * nsh_499[k]
                   + f_3 * pc_x[k] * osh_499[k];

        t_661[k] = f_22 * nsh_500[k]
                   + f_3 * pc_x[k] * osh_500[k];

        t_662[k] = f_22 * nsh_501[k]
                   + f_3 * pc_x[k] * osh_501[k];

        t_663[k] = f_22 * nsh_502[k]
                   + f_3 * pc_x[k] * osh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, nsh_351, nsh_372, nsh_503, \
                         osg0_355, osg1_355, osh_498, osh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_22 * nsh_503[k]
                   + f_3 * pc_x[k] * osh_503[k];

        t_665[k] = f_14 * nsh_372[k]
                   + f_1 * osg0_355[k]
                   - f_2 * osg1_355[k]
                   + f_3 * pc_y[k] * osh_498[k];

        t_666[k] = f_12 * nsh_351[k]
                   + f_3 * pc_z[k] * osh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, nsh_374, nsh_375, nsh_376, osg0_357, \
                         osg0_358, osg0_359, osg1_357, osg1_358, osg1_359, osh_500, osh_501, \
                         osh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * nsh_374[k]
                   + f_8 * osg0_357[k]
                   - f_9 * osg1_357[k]
                   + f_3 * pc_y[k] * osh_500[k];

        t_668[k] = f_14 * nsh_375[k]
                   + f_6 * osg0_358[k]
                   - f_7 * osg1_358[k]
                   + f_3 * pc_y[k] * osh_501[k];

        t_669[k] = f_14 * nsh_376[k]
                   + f_4 * osg0_359[k]
                   - f_5 * osg1_359[k]
                   + f_3 * pc_y[k] * osh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, nsh_356, nsh_377, nsh_504, \
                         osg0_359, osg0_360, osg1_359, osg1_360, osh_503, \
                         osh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * nsh_377[k]
                   + f_3 * pc_y[k] * osh_503[k];

        t_671[k] = f_12 * nsh_356[k]
                   + f_1 * osg0_359[k]
                   - f_2 * osg1_359[k]
                   + f_3 * pc_z[k] * osh_503[k];

        t_672[k] = f_22 * nsh_504[k]
                   + f_1 * osg0_360[k]
                   - f_2 * osg1_360[k]
                   + f_3 * pc_x[k] * osh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, nsh_357, nsh_378, \
                         nsh_380, nsh_507, osg0_363, osg1_363, osh_504, osh_506, \
                         osh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * nsh_378[k]
                   + f_3 * pc_y[k] * osh_504[k];

        t_674[k] = f_13 * nsh_357[k]
                   + f_3 * pc_z[k] * osh_504[k];

        t_675[k] = f_22 * nsh_507[k]
                   + f_8 * osg0_363[k]
                   - f_9 * osg1_363[k]
                   + f_3 * pc_x[k] * osh_507[k];

        t_676[k] = f_13 * nsh_380[k]
                   + f_3 * pc_y[k] * osh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, nsh_360, nsh_509, nsh_510, osg0_365, \
                         osg0_366, osg1_365, osg1_366, osh_507, osh_509, \
                         osh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_22 * nsh_509[k]
                   + f_8 * osg0_365[k]
                   - f_9 * osg1_365[k]
                   + f_3 * pc_x[k] * osh_509[k];

        t_678[k] = f_22 * nsh_510[k]
                   + f_6 * osg0_366[k]
                   - f_7 * osg1_366[k]
                   + f_3 * pc_x[k] * osh_510[k];

        t_679[k] = f_13 * nsh_360[k]
                   + f_3 * pc_z[k] * osh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, nsh_383, nsh_513, nsh_514, osg0_369, \
                         osg0_370, osg1_369, osg1_370, osh_509, osh_513, \
                         osh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * nsh_383[k]
                   + f_3 * pc_y[k] * osh_509[k];

        t_681[k] = f_22 * nsh_513[k]
                   + f_6 * osg0_369[k]
                   - f_7 * osg1_369[k]
                   + f_3 * pc_x[k] * osh_513[k];

        t_682[k] = f_22 * nsh_514[k]
                   + f_4 * osg0_370[k]
                   - f_5 * osg1_370[k]
                   + f_3 * pc_x[k] * osh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, nsh_363, nsh_387, nsh_516, \
                         osg0_372, osg1_372, osh_510, osh_513, \
                         osh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * nsh_363[k]
                   + f_3 * pc_z[k] * osh_510[k];

        t_684[k] = f_22 * nsh_516[k]
                   + f_4 * osg0_372[k]
                   - f_5 * osg1_372[k]
                   + f_3 * pc_x[k] * osh_516[k];

        t_685[k] = f_13 * nsh_387[k]
                   + f_3 * pc_y[k] * osh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, nsh_518, nsh_519, nsh_520, nsh_521, \
                         osg0_374, osg1_374, osh_518, osh_519, osh_520, \
                         osh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_22 * nsh_518[k]
                   + f_4 * osg0_374[k]
                   - f_5 * osg1_374[k]
                   + f_3 * pc_x[k] * osh_518[k];

        t_687[k] = f_22 * nsh_519[k]
                   + f_3 * pc_x[k] * osh_519[k];

        t_688[k] = f_22 * nsh_520[k]
                   + f_3 * pc_x[k] * osh_520[k];

        t_689[k] = f_22 * nsh_521[k]
                   + f_3 * pc_x[k] * osh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, nsh_393, nsh_522, nsh_523, \
                         nsh_524, osg0_370, osg1_370, osh_519, osh_522, osh_523, \
                         osh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_22 * nsh_522[k]
                   + f_3 * pc_x[k] * osh_522[k];

        t_691[k] = f_22 * nsh_523[k]
                   + f_3 * pc_x[k] * osh_523[k];

        t_692[k] = f_22 * nsh_524[k]
                   + f_3 * pc_x[k] * osh_524[k];

        t_693[k] = f_13 * nsh_393[k]
                   + f_1 * osg0_370[k]
                   - f_2 * osg1_370[k]
                   + f_3 * pc_y[k] * osh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, nsh_372, nsh_395, nsh_396, osg0_372, \
                         osg0_373, osg1_372, osg1_373, osh_519, osh_521, \
                         osh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * nsh_372[k]
                   + f_3 * pc_z[k] * osh_519[k];

        t_695[k] = f_13 * nsh_395[k]
                   + f_8 * osg0_372[k]
                   - f_9 * osg1_372[k]
                   + f_3 * pc_y[k] * osh_521[k];

        t_696[k] = f_13 * nsh_396[k]
                   + f_6 * osg0_373[k]
                   - f_7 * osg1_373[k]
                   + f_3 * pc_y[k] * osh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, nsh_377, nsh_397, nsh_398, osg0_374, \
                         osg1_374, osh_523, osh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * nsh_397[k]
                   + f_4 * osg0_374[k]
                   - f_5 * osg1_374[k]
                   + f_3 * pc_y[k] * osh_523[k];

        t_698[k] = f_13 * nsh_398[k]
                   + f_3 * pc_y[k] * osh_524[k];

        t_699[k] = f_13 * nsh_377[k]
                   + f_1 * osg0_374[k]
                   - f_2 * osg1_374[k]
                   + f_3 * pc_z[k] * osh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, nsh_378, nsh_399, nsh_525, \
                         osg0_375, osg1_375, osh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_22 * nsh_525[k]
                   + f_1 * osg0_375[k]
                   - f_2 * osg1_375[k]
                   + f_3 * pc_x[k] * osh_525[k];

        t_701[k] = f_12 * nsh_399[k]
                   + f_3 * pc_y[k] * osh_525[k];

        t_702[k] = f_14 * nsh_378[k]
                   + f_3 * pc_z[k] * osh_525[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_560 = buffer.data(nsi0 + 560);
    const auto *nsi0_563 = buffer.data(nsi0 + 563);
    const auto *nsi0_565 = buffer.data(nsi0 + 565);
    const auto *nsi0_566 = buffer.data(nsi0 + 566);
    const auto *nsi0_569 = buffer.data(nsi0 + 569);
    const auto *nsi0_570 = buffer.data(nsi0 + 570);
    const auto *nsi0_572 = buffer.data(nsi0 + 572);
    const auto *nsi0_574 = buffer.data(nsi0 + 574);
    const auto *nsi0_587 = buffer.data(nsi0 + 587);
    const auto *nsi0_588 = buffer.data(nsi0 + 588);
    const auto *nsi0_591 = buffer.data(nsi0 + 591);
    const auto *nsi0_594 = buffer.data(nsi0 + 594);

    const auto *nsh_381 = buffer.data(nsh + 381);
    const auto *nsh_384 = buffer.data(nsh + 384);
    const auto *nsh_393 = buffer.data(nsh + 393);
    const auto *nsh_398 = buffer.data(nsh + 398);
    const auto *nsh_399 = buffer.data(nsh + 399);
    const auto *nsh_401 = buffer.data(nsh + 401);
    const auto *nsh_402 = buffer.data(nsh + 402);
    const auto *nsh_404 = buffer.data(nsh + 404);
    const auto *nsh_405 = buffer.data(nsh + 405);
    const auto *nsh_408 = buffer.data(nsh + 408);
    const auto *nsh_414 = buffer.data(nsh + 414);
    const auto *nsh_416 = buffer.data(nsh + 416);
    const auto *nsh_417 = buffer.data(nsh + 417);
    const auto *nsh_418 = buffer.data(nsh + 418);
    const auto *nsh_419 = buffer.data(nsh + 419);
    const auto *nsh_420 = buffer.data(nsh + 420);
    const auto *nsh_421 = buffer.data(nsh + 421);
    const auto *nsh_422 = buffer.data(nsh + 422);
    const auto *nsh_423 = buffer.data(nsh + 423);
    const auto *nsh_425 = buffer.data(nsh + 425);
    const auto *nsh_426 = buffer.data(nsh + 426);
    const auto *nsh_428 = buffer.data(nsh + 428);
    const auto *nsh_429 = buffer.data(nsh + 429);
    const auto *nsh_435 = buffer.data(nsh + 435);
    const auto *nsh_437 = buffer.data(nsh + 437);
    const auto *nsh_438 = buffer.data(nsh + 438);
    const auto *nsh_439 = buffer.data(nsh + 439);
    const auto *nsh_440 = buffer.data(nsh + 440);
    const auto *nsh_441 = buffer.data(nsh + 441);
    const auto *nsh_444 = buffer.data(nsh + 444);
    const auto *nsh_446 = buffer.data(nsh + 446);
    const auto *nsh_450 = buffer.data(nsh + 450);
    const auto *nsh_456 = buffer.data(nsh + 456);
    const auto *nsh_461 = buffer.data(nsh + 461);
    const auto *nsh_462 = buffer.data(nsh + 462);
    const auto *nsh_464 = buffer.data(nsh + 464);
    const auto *nsh_528 = buffer.data(nsh + 528);
    const auto *nsh_530 = buffer.data(nsh + 530);
    const auto *nsh_531 = buffer.data(nsh + 531);
    const auto *nsh_534 = buffer.data(nsh + 534);
    const auto *nsh_535 = buffer.data(nsh + 535);
    const auto *nsh_537 = buffer.data(nsh + 537);
    const auto *nsh_539 = buffer.data(nsh + 539);
    const auto *nsh_540 = buffer.data(nsh + 540);
    const auto *nsh_541 = buffer.data(nsh + 541);
    const auto *nsh_542 = buffer.data(nsh + 542);
    const auto *nsh_543 = buffer.data(nsh + 543);
    const auto *nsh_544 = buffer.data(nsh + 544);
    const auto *nsh_545 = buffer.data(nsh + 545);
    const auto *nsh_561 = buffer.data(nsh + 561);
    const auto *nsh_562 = buffer.data(nsh + 562);
    const auto *nsh_563 = buffer.data(nsh + 563);
    const auto *nsh_564 = buffer.data(nsh + 564);
    const auto *nsh_565 = buffer.data(nsh + 565);
    const auto *nsh_566 = buffer.data(nsh + 566);
    const auto *nsh_567 = buffer.data(nsh + 567);
    const auto *nsh_572 = buffer.data(nsh + 572);
    const auto *nsh_576 = buffer.data(nsh + 576);
    const auto *nsh_581 = buffer.data(nsh + 581);
    const auto *nsh_582 = buffer.data(nsh + 582);
    const auto *nsh_583 = buffer.data(nsh + 583);
    const auto *nsh_584 = buffer.data(nsh + 584);
    const auto *nsh_585 = buffer.data(nsh + 585);
    const auto *nsh_587 = buffer.data(nsh + 587);
    const auto *nsh_588 = buffer.data(nsh + 588);
    const auto *nsh_591 = buffer.data(nsh + 591);
    const auto *nsh_594 = buffer.data(nsh + 594);
    const auto *nsh_598 = buffer.data(nsh + 598);
    const auto *nsh_603 = buffer.data(nsh + 603);
    const auto *nsh_605 = buffer.data(nsh + 605);
    const auto *nsh_606 = buffer.data(nsh + 606);
    const auto *nsh_607 = buffer.data(nsh + 607);
    const auto *nsh_608 = buffer.data(nsh + 608);
    const auto *nsh_614 = buffer.data(nsh + 614);

    const auto *nsi1_560 = buffer.data(nsi1 + 560);
    const auto *nsi1_563 = buffer.data(nsi1 + 563);
    const auto *nsi1_565 = buffer.data(nsi1 + 565);
    const auto *nsi1_566 = buffer.data(nsi1 + 566);
    const auto *nsi1_569 = buffer.data(nsi1 + 569);
    const auto *nsi1_570 = buffer.data(nsi1 + 570);
    const auto *nsi1_572 = buffer.data(nsi1 + 572);
    const auto *nsi1_574 = buffer.data(nsi1 + 574);
    const auto *nsi1_587 = buffer.data(nsi1 + 587);
    const auto *nsi1_588 = buffer.data(nsi1 + 588);
    const auto *nsi1_591 = buffer.data(nsi1 + 591);
    const auto *nsi1_594 = buffer.data(nsi1 + 594);

    const auto *osg0_378 = buffer.data(osg0 + 378);
    const auto *osg0_380 = buffer.data(osg0 + 380);
    const auto *osg0_381 = buffer.data(osg0 + 381);
    const auto *osg0_384 = buffer.data(osg0 + 384);
    const auto *osg0_385 = buffer.data(osg0 + 385);
    const auto *osg0_387 = buffer.data(osg0 + 387);
    const auto *osg0_388 = buffer.data(osg0 + 388);
    const auto *osg0_389 = buffer.data(osg0 + 389);
    const auto *osg0_400 = buffer.data(osg0 + 400);
    const auto *osg0_402 = buffer.data(osg0 + 402);
    const auto *osg0_403 = buffer.data(osg0 + 403);
    const auto *osg0_404 = buffer.data(osg0 + 404);
    const auto *osg0_405 = buffer.data(osg0 + 405);
    const auto *osg0_406 = buffer.data(osg0 + 406);
    const auto *osg0_407 = buffer.data(osg0 + 407);
    const auto *osg0_408 = buffer.data(osg0 + 408);
    const auto *osg0_409 = buffer.data(osg0 + 409);
    const auto *osg0_410 = buffer.data(osg0 + 410);
    const auto *osg0_414 = buffer.data(osg0 + 414);
    const auto *osg0_415 = buffer.data(osg0 + 415);
    const auto *osg0_416 = buffer.data(osg0 + 416);
    const auto *osg0_417 = buffer.data(osg0 + 417);
    const auto *osg0_418 = buffer.data(osg0 + 418);
    const auto *osg0_419 = buffer.data(osg0 + 419);
    const auto *osg0_420 = buffer.data(osg0 + 420);
    const auto *osg0_422 = buffer.data(osg0 + 422);
    const auto *osg0_423 = buffer.data(osg0 + 423);
    const auto *osg0_425 = buffer.data(osg0 + 425);
    const auto *osg0_426 = buffer.data(osg0 + 426);
    const auto *osg0_430 = buffer.data(osg0 + 430);
    const auto *osg0_431 = buffer.data(osg0 + 431);
    const auto *osg0_432 = buffer.data(osg0 + 432);
    const auto *osg0_434 = buffer.data(osg0 + 434);
    const auto *osg0_440 = buffer.data(osg0 + 440);

    const auto *osg1_378 = buffer.data(osg1 + 378);
    const auto *osg1_380 = buffer.data(osg1 + 380);
    const auto *osg1_381 = buffer.data(osg1 + 381);
    const auto *osg1_384 = buffer.data(osg1 + 384);
    const auto *osg1_385 = buffer.data(osg1 + 385);
    const auto *osg1_387 = buffer.data(osg1 + 387);
    const auto *osg1_388 = buffer.data(osg1 + 388);
    const auto *osg1_389 = buffer.data(osg1 + 389);
    const auto *osg1_400 = buffer.data(osg1 + 400);
    const auto *osg1_402 = buffer.data(osg1 + 402);
    const auto *osg1_403 = buffer.data(osg1 + 403);
    const auto *osg1_404 = buffer.data(osg1 + 404);
    const auto *osg1_405 = buffer.data(osg1 + 405);
    const auto *osg1_406 = buffer.data(osg1 + 406);
    const auto *osg1_407 = buffer.data(osg1 + 407);
    const auto *osg1_408 = buffer.data(osg1 + 408);
    const auto *osg1_409 = buffer.data(osg1 + 409);
    const auto *osg1_410 = buffer.data(osg1 + 410);
    const auto *osg1_414 = buffer.data(osg1 + 414);
    const auto *osg1_415 = buffer.data(osg1 + 415);
    const auto *osg1_416 = buffer.data(osg1 + 416);
    const auto *osg1_417 = buffer.data(osg1 + 417);
    const auto *osg1_418 = buffer.data(osg1 + 418);
    const auto *osg1_419 = buffer.data(osg1 + 419);
    const auto *osg1_420 = buffer.data(osg1 + 420);
    const auto *osg1_422 = buffer.data(osg1 + 422);
    const auto *osg1_423 = buffer.data(osg1 + 423);
    const auto *osg1_425 = buffer.data(osg1 + 425);
    const auto *osg1_426 = buffer.data(osg1 + 426);
    const auto *osg1_430 = buffer.data(osg1 + 430);
    const auto *osg1_431 = buffer.data(osg1 + 431);
    const auto *osg1_432 = buffer.data(osg1 + 432);
    const auto *osg1_434 = buffer.data(osg1 + 434);
    const auto *osg1_440 = buffer.data(osg1 + 440);

    const auto *osh_527 = buffer.data(osh + 527);
    const auto *osh_528 = buffer.data(osh + 528);
    const auto *osh_530 = buffer.data(osh + 530);
    const auto *osh_531 = buffer.data(osh + 531);
    const auto *osh_534 = buffer.data(osh + 534);
    const auto *osh_535 = buffer.data(osh + 535);
    const auto *osh_537 = buffer.data(osh + 537);
    const auto *osh_539 = buffer.data(osh + 539);
    const auto *osh_540 = buffer.data(osh + 540);
    const auto *osh_541 = buffer.data(osh + 541);
    const auto *osh_542 = buffer.data(osh + 542);
    const auto *osh_543 = buffer.data(osh + 543);
    const auto *osh_544 = buffer.data(osh + 544);
    const auto *osh_545 = buffer.data(osh + 545);
    const auto *osh_546 = buffer.data(osh + 546);
    const auto *osh_548 = buffer.data(osh + 548);
    const auto *osh_549 = buffer.data(osh + 549);
    const auto *osh_551 = buffer.data(osh + 551);
    const auto *osh_552 = buffer.data(osh + 552);
    const auto *osh_555 = buffer.data(osh + 555);
    const auto *osh_561 = buffer.data(osh + 561);
    const auto *osh_562 = buffer.data(osh + 562);
    const auto *osh_563 = buffer.data(osh + 563);
    const auto *osh_564 = buffer.data(osh + 564);
    const auto *osh_565 = buffer.data(osh + 565);
    const auto *osh_566 = buffer.data(osh + 566);
    const auto *osh_567 = buffer.data(osh + 567);
    const auto *osh_568 = buffer.data(osh + 568);
    const auto *osh_569 = buffer.data(osh + 569);
    const auto *osh_570 = buffer.data(osh + 570);
    const auto *osh_571 = buffer.data(osh + 571);
    const auto *osh_572 = buffer.data(osh + 572);
    const auto *osh_573 = buffer.data(osh + 573);
    const auto *osh_574 = buffer.data(osh + 574);
    const auto *osh_575 = buffer.data(osh + 575);
    const auto *osh_576 = buffer.data(osh + 576);
    const auto *osh_581 = buffer.data(osh + 581);
    const auto *osh_582 = buffer.data(osh + 582);
    const auto *osh_583 = buffer.data(osh + 583);
    const auto *osh_584 = buffer.data(osh + 584);
    const auto *osh_585 = buffer.data(osh + 585);
    const auto *osh_586 = buffer.data(osh + 586);
    const auto *osh_587 = buffer.data(osh + 587);
    const auto *osh_588 = buffer.data(osh + 588);
    const auto *osh_589 = buffer.data(osh + 589);
    const auto *osh_590 = buffer.data(osh + 590);
    const auto *osh_591 = buffer.data(osh + 591);
    const auto *osh_593 = buffer.data(osh + 593);
    const auto *osh_594 = buffer.data(osh + 594);
    const auto *osh_595 = buffer.data(osh + 595);
    const auto *osh_597 = buffer.data(osh + 597);
    const auto *osh_598 = buffer.data(osh + 598);
    const auto *osh_603 = buffer.data(osh + 603);
    const auto *osh_604 = buffer.data(osh + 604);
    const auto *osh_605 = buffer.data(osh + 605);
    const auto *osh_606 = buffer.data(osh + 606);
    const auto *osh_607 = buffer.data(osh + 607);
    const auto *osh_608 = buffer.data(osh + 608);
    const auto *osh_609 = buffer.data(osh + 609);
    const auto *osh_611 = buffer.data(osh + 611);
    const auto *osh_612 = buffer.data(osh + 612);
    const auto *osh_614 = buffer.data(osh + 614);

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, nsh_401, nsh_528, nsh_530, osg0_378, \
                         osg0_380, osg1_378, osg1_380, osh_527, osh_528, \
                         osh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_22 * nsh_528[k]
                   + f_8 * osg0_378[k]
                   - f_9 * osg1_378[k]
                   + f_3 * pc_x[k] * osh_528[k];

        t_704[k] = f_12 * nsh_401[k]
                   + f_3 * pc_y[k] * osh_527[k];

        t_705[k] = f_22 * nsh_530[k]
                   + f_8 * osg0_380[k]
                   - f_9 * osg1_380[k]
                   + f_3 * pc_x[k] * osh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, nsh_381, nsh_404, nsh_531, \
                         osg0_381, osg1_381, osh_528, osh_530, \
                         osh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_22 * nsh_531[k]
                   + f_6 * osg0_381[k]
                   - f_7 * osg1_381[k]
                   + f_3 * pc_x[k] * osh_531[k];

        t_707[k] = f_14 * nsh_381[k]
                   + f_3 * pc_z[k] * osh_528[k];

        t_708[k] = f_12 * nsh_404[k]
                   + f_3 * pc_y[k] * osh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, nsh_384, nsh_534, nsh_535, osg0_384, \
                         osg0_385, osg1_384, osg1_385, osh_531, osh_534, \
                         osh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_22 * nsh_534[k]
                   + f_6 * osg0_384[k]
                   - f_7 * osg1_384[k]
                   + f_3 * pc_x[k] * osh_534[k];

        t_710[k] = f_22 * nsh_535[k]
                   + f_4 * osg0_385[k]
                   - f_5 * osg1_385[k]
                   + f_3 * pc_x[k] * osh_535[k];

        t_711[k] = f_14 * nsh_384[k]
                   + f_3 * pc_z[k] * osh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, nsh_408, nsh_537, nsh_539, osg0_387, \
                         osg0_389, osg1_387, osg1_389, osh_534, osh_537, \
                         osh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_22 * nsh_537[k]
                   + f_4 * osg0_387[k]
                   - f_5 * osg1_387[k]
                   + f_3 * pc_x[k] * osh_537[k];

        t_713[k] = f_12 * nsh_408[k]
                   + f_3 * pc_y[k] * osh_534[k];

        t_714[k] = f_22 * nsh_539[k]
                   + f_4 * osg0_389[k]
                   - f_5 * osg1_389[k]
                   + f_3 * pc_x[k] * osh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, nsh_540, nsh_541, nsh_542, \
                         nsh_543, nsh_544, osh_540, osh_541, osh_542, osh_543, \
                         osh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_22 * nsh_540[k]
                   + f_3 * pc_x[k] * osh_540[k];

        t_716[k] = f_22 * nsh_541[k]
                   + f_3 * pc_x[k] * osh_541[k];

        t_717[k] = f_22 * nsh_542[k]
                   + f_3 * pc_x[k] * osh_542[k];

        t_718[k] = f_22 * nsh_543[k]
                   + f_3 * pc_x[k] * osh_543[k];

        t_719[k] = f_22 * nsh_544[k]
                   + f_3 * pc_x[k] * osh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, nsh_393, nsh_414, nsh_545, \
                         osg0_385, osg1_385, osh_540, osh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_22 * nsh_545[k]
                   + f_3 * pc_x[k] * osh_545[k];

        t_721[k] = f_12 * nsh_414[k]
                   + f_1 * osg0_385[k]
                   - f_2 * osg1_385[k]
                   + f_3 * pc_y[k] * osh_540[k];

        t_722[k] = f_14 * nsh_393[k]
                   + f_3 * pc_z[k] * osh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, nsh_416, nsh_417, nsh_418, osg0_387, \
                         osg0_388, osg0_389, osg1_387, osg1_388, osg1_389, osh_542, osh_543, \
                         osh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * nsh_416[k]
                   + f_8 * osg0_387[k]
                   - f_9 * osg1_387[k]
                   + f_3 * pc_y[k] * osh_542[k];

        t_724[k] = f_12 * nsh_417[k]
                   + f_6 * osg0_388[k]
                   - f_7 * osg1_388[k]
                   + f_3 * pc_y[k] * osh_543[k];

        t_725[k] = f_12 * nsh_418[k]
                   + f_4 * osg0_389[k]
                   - f_5 * osg1_389[k]
                   + f_3 * pc_y[k] * osh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_y, pc_y, pc_z, nsi0_560, nsh_398, \
                         nsh_419, nsh_420, nsi1_560, osg0_389, osg1_389, osh_545, \
                         osh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * nsh_419[k]
                   + f_3 * pc_y[k] * osh_545[k];

        t_727[k] = f_14 * nsh_398[k]
                   + f_1 * osg0_389[k]
                   - f_2 * osg1_389[k]
                   + f_3 * pc_z[k] * osh_545[k];

        t_728[k] = pa_y[k] * nsi0_560[k]
                   - f_10 * pc_y[k] * nsi1_560[k];

        t_729[k] = f_11 * nsh_420[k]
                   + f_3 * pc_y[k] * osh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_y, pc_y, pc_z, nsi0_563, nsi0_565, \
                         nsh_399, nsh_421, nsh_422, nsi1_563, nsi1_565, osh_546, \
                         osh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_22 * nsh_399[k]
                   + f_3 * pc_z[k] * osh_546[k];

        t_731[k] = pa_y[k] * nsi0_563[k]
                   + f_12 * nsh_421[k]
                   - f_10 * pc_y[k] * nsi1_563[k];

        t_732[k] = f_11 * nsh_422[k]
                   + f_3 * pc_y[k] * osh_548[k];

        t_733[k] = pa_y[k] * nsi0_565[k]
                   - f_10 * pc_y[k] * nsi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_y, pc_y, pc_z, nsi0_566, nsi0_569, \
                         nsh_402, nsh_423, nsh_425, nsi1_566, nsi1_569, osh_549, \
                         osh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_y[k] * nsi0_566[k]
                   + f_13 * nsh_423[k]
                   - f_10 * pc_y[k] * nsi1_566[k];

        t_735[k] = f_22 * nsh_402[k]
                   + f_3 * pc_z[k] * osh_549[k];

        t_736[k] = f_11 * nsh_425[k]
                   + f_3 * pc_y[k] * osh_551[k];

        t_737[k] = pa_y[k] * nsi0_569[k]
                   - f_10 * pc_y[k] * nsi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pa_y, pc_y, pc_z, nsi0_570, nsi0_572, nsh_405, \
                         nsh_426, nsh_428, nsi1_570, nsi1_572, \
                         osh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pa_y[k] * nsi0_570[k]
                   + f_14 * nsh_426[k]
                   - f_10 * pc_y[k] * nsi1_570[k];

        t_739[k] = f_22 * nsh_405[k]
                   + f_3 * pc_z[k] * osh_552[k];

        t_740[k] = pa_y[k] * nsi0_572[k]
                   + f_12 * nsh_428[k]
                   - f_10 * pc_y[k] * nsi1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_y, pc_x, pc_y, nsi0_574, nsh_429, \
                         nsh_561, nsh_562, nsi1_574, osh_555, osh_561, \
                         osh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * nsh_429[k]
                   + f_3 * pc_y[k] * osh_555[k];

        t_742[k] = pa_y[k] * nsi0_574[k]
                   - f_10 * pc_y[k] * nsi1_574[k];

        t_743[k] = f_22 * nsh_561[k]
                   + f_3 * pc_x[k] * osh_561[k];

        t_744[k] = f_22 * nsh_562[k]
                   + f_3 * pc_x[k] * osh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, nsh_563, nsh_564, nsh_565, nsh_566, \
                         osh_563, osh_564, osh_565, osh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_22 * nsh_563[k]
                   + f_3 * pc_x[k] * osh_563[k];

        t_746[k] = f_22 * nsh_564[k]
                   + f_3 * pc_x[k] * osh_564[k];

        t_747[k] = f_22 * nsh_565[k]
                   + f_3 * pc_x[k] * osh_565[k];

        t_748[k] = f_22 * nsh_566[k]
                   + f_3 * pc_x[k] * osh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, nsh_414, nsh_435, nsh_437, osg0_400, \
                         osg0_402, osg1_400, osg1_402, osh_561, \
                         osh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * nsh_435[k]
                   + f_1 * osg0_400[k]
                   - f_2 * osg1_400[k]
                   + f_3 * pc_y[k] * osh_561[k];

        t_750[k] = f_22 * nsh_414[k]
                   + f_3 * pc_z[k] * osh_561[k];

        t_751[k] = f_11 * nsh_437[k]
                   + f_8 * osg0_402[k]
                   - f_9 * osg1_402[k]
                   + f_3 * pc_y[k] * osh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, nsh_438, nsh_439, nsh_440, osg0_403, \
                         osg0_404, osg1_403, osg1_404, osh_564, osh_565, \
                         osh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * nsh_438[k]
                   + f_6 * osg0_403[k]
                   - f_7 * osg1_403[k]
                   + f_3 * pc_y[k] * osh_564[k];

        t_753[k] = f_11 * nsh_439[k]
                   + f_4 * osg0_404[k]
                   - f_5 * osg1_404[k]
                   + f_3 * pc_y[k] * osh_565[k];

        t_754[k] = f_11 * nsh_440[k]
                   + f_3 * pc_y[k] * osh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_y, pc_x, pc_y, pc_z, nsi0_587, \
                         nsh_420, nsh_567, nsi1_587, osg0_405, osg1_405, \
                         osh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_y[k] * nsi0_587[k]
                   - f_10 * pc_y[k] * nsi1_587[k];

        t_756[k] = f_22 * nsh_567[k]
                   + f_1 * osg0_405[k]
                   - f_2 * osg1_405[k]
                   + f_3 * pc_x[k] * osh_567[k];

        t_757[k] = f_3 * pc_y[k] * osh_567[k];

        t_758[k] = f_21 * nsh_420[k]
                   + f_3 * pc_z[k] * osh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, nsh_572, osg0_405, osg0_410, \
                         osg1_405, osg1_410, osh_568, osh_569, \
                         osh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_4 * osg0_405[k]
                   - f_5 * osg1_405[k]
                   + f_3 * pc_y[k] * osh_568[k];

        t_760[k] = f_3 * pc_y[k] * osh_569[k];

        t_761[k] = f_22 * nsh_572[k]
                   + f_8 * osg0_410[k]
                   - f_9 * osg1_410[k]
                   + f_3 * pc_x[k] * osh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, osg0_406, osg0_407, osg1_406, osg1_407, \
                         osh_570, osh_571, osh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_6 * osg0_406[k]
                   - f_7 * osg1_406[k]
                   + f_3 * pc_y[k] * osh_570[k];

        t_763[k] = f_4 * osg0_407[k]
                   - f_5 * osg1_407[k]
                   + f_3 * pc_y[k] * osh_571[k];

        t_764[k] = f_3 * pc_y[k] * osh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_y, nsh_576, osg0_408, osg0_409, \
                         osg0_414, osg1_408, osg1_409, osg1_414, osh_573, osh_574, \
                         osh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_22 * nsh_576[k]
                   + f_6 * osg0_414[k]
                   - f_7 * osg1_414[k]
                   + f_3 * pc_x[k] * osh_576[k];

        t_766[k] = f_8 * osg0_408[k]
                   - f_9 * osg1_408[k]
                   + f_3 * pc_y[k] * osh_573[k];

        t_767[k] = f_6 * osg0_409[k]
                   - f_7 * osg1_409[k]
                   + f_3 * pc_y[k] * osh_574[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, nsh_581, nsh_582, osg0_410, \
                         osg0_419, osg1_410, osg1_419, osh_575, osh_576, osh_581, \
                         osh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * osg0_410[k]
                   - f_5 * osg1_410[k]
                   + f_3 * pc_y[k] * osh_575[k];

        t_769[k] = f_3 * pc_y[k] * osh_576[k];

        t_770[k] = f_22 * nsh_581[k]
                   + f_4 * osg0_419[k]
                   - f_5 * osg1_419[k]
                   + f_3 * pc_x[k] * osh_581[k];

        t_771[k] = f_22 * nsh_582[k]
                   + f_3 * pc_x[k] * osh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pc_y, nsh_583, nsh_584, \
                         nsh_585, nsh_587, osh_581, osh_583, osh_584, osh_585, \
                         osh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_22 * nsh_583[k]
                   + f_3 * pc_x[k] * osh_583[k];

        t_773[k] = f_22 * nsh_584[k]
                   + f_3 * pc_x[k] * osh_584[k];

        t_774[k] = f_22 * nsh_585[k]
                   + f_3 * pc_x[k] * osh_585[k];

        t_775[k] = f_3 * pc_y[k] * osh_581[k];

        t_776[k] = f_22 * nsh_587[k]
                   + f_3 * pc_x[k] * osh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, osg0_415, osg0_416, osg0_417, osg1_415, \
                         osg1_416, osg1_417, osh_582, osh_583, \
                         osh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_1 * osg0_415[k]
                   - f_2 * osg1_415[k]
                   + f_3 * pc_y[k] * osh_582[k];

        t_778[k] = f_16 * osg0_416[k]
                   - f_17 * osg1_416[k]
                   + f_3 * pc_y[k] * osh_583[k];

        t_779[k] = f_8 * osg0_417[k]
                   - f_9 * osg1_417[k]
                   + f_3 * pc_y[k] * osh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, nsh_440, osg0_418, osg0_419, \
                         osg1_418, osg1_419, osh_585, osh_586, \
                         osh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * osg0_418[k]
                   - f_7 * osg1_418[k]
                   + f_3 * pc_y[k] * osh_585[k];

        t_781[k] = f_4 * osg0_419[k]
                   - f_5 * osg1_419[k]
                   + f_3 * pc_y[k] * osh_586[k];

        t_782[k] = f_3 * pc_y[k] * osh_587[k];

        t_783[k] = f_21 * nsh_440[k]
                   + f_1 * osg0_419[k]
                   - f_2 * osg1_419[k]
                   + f_3 * pc_z[k] * osh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, nsh_441, nsh_588, \
                         nsh_591, osg0_420, osg0_423, osg1_420, osg1_423, osh_588, \
                         osh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_14 * nsh_588[k]
                   + f_1 * osg0_420[k]
                   - f_2 * osg1_420[k]
                   + f_3 * pc_x[k] * osh_588[k];

        t_785[k] = f_20 * nsh_441[k]
                   + f_3 * pc_y[k] * osh_588[k];

        t_786[k] = f_3 * pc_z[k] * osh_588[k];

        t_787[k] = f_14 * nsh_591[k]
                   + f_8 * osg0_423[k]
                   - f_9 * osg1_423[k]
                   + f_3 * pc_x[k] * osh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pc_x, pc_z, nsh_594, osg0_420, osg0_426, \
                         osg1_420, osg1_426, osh_589, osh_590, osh_591, \
                         osh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_3 * pc_z[k] * osh_589[k];

        t_789[k] = f_4 * osg0_420[k]
                   - f_5 * osg1_420[k]
                   + f_3 * pc_z[k] * osh_590[k];

        t_790[k] = f_14 * nsh_594[k]
                   + f_6 * osg0_426[k]
                   - f_7 * osg1_426[k]
                   + f_3 * pc_x[k] * osh_594[k];

        t_791[k] = f_3 * pc_z[k] * osh_591[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pc_x, pc_y, pc_z, nsh_446, nsh_598, \
                         osg0_422, osg0_430, osg1_422, osg1_430, osh_593, osh_594, \
                         osh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_20 * nsh_446[k]
                   + f_3 * pc_y[k] * osh_593[k];

        t_793[k] = f_6 * osg0_422[k]
                   - f_7 * osg1_422[k]
                   + f_3 * pc_z[k] * osh_593[k];

        t_794[k] = f_14 * nsh_598[k]
                   + f_4 * osg0_430[k]
                   - f_5 * osg1_430[k]
                   + f_3 * pc_x[k] * osh_598[k];

        t_795[k] = f_3 * pc_z[k] * osh_594[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, nsh_450, nsh_603, \
                         osg0_423, osg0_425, osg1_423, osg1_425, osh_595, osh_597, \
                         osh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_4 * osg0_423[k]
                   - f_5 * osg1_423[k]
                   + f_3 * pc_z[k] * osh_595[k];

        t_797[k] = f_20 * nsh_450[k]
                   + f_3 * pc_y[k] * osh_597[k];

        t_798[k] = f_8 * osg0_425[k]
                   - f_9 * osg1_425[k]
                   + f_3 * pc_z[k] * osh_597[k];

        t_799[k] = f_14 * nsh_603[k]
                   + f_3 * pc_x[k] * osh_603[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pc_x, pc_z, nsh_605, nsh_606, \
                         nsh_607, nsh_608, osh_598, osh_605, osh_606, osh_607, \
                         osh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_3 * pc_z[k] * osh_598[k];

        t_801[k] = f_14 * nsh_605[k]
                   + f_3 * pc_x[k] * osh_605[k];

        t_802[k] = f_14 * nsh_606[k]
                   + f_3 * pc_x[k] * osh_606[k];

        t_803[k] = f_14 * nsh_607[k]
                   + f_3 * pc_x[k] * osh_607[k];

        t_804[k] = f_14 * nsh_608[k]
                   + f_3 * pc_x[k] * osh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pc_y, pc_z, nsh_456, osg0_430, osg0_431, \
                         osg1_430, osg1_431, osh_603, osh_604, \
                         osh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_20 * nsh_456[k]
                   + f_1 * osg0_430[k]
                   - f_2 * osg1_430[k]
                   + f_3 * pc_y[k] * osh_603[k];

        t_806[k] = f_3 * pc_z[k] * osh_603[k];

        t_807[k] = f_4 * osg0_430[k]
                   - f_5 * osg1_430[k]
                   + f_3 * pc_z[k] * osh_604[k];

        t_808[k] = f_6 * osg0_431[k]
                   - f_7 * osg1_431[k]
                   + f_3 * pc_z[k] * osh_605[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pa_z, pc_y, pc_z, nsi0_588, nsh_461, \
                         nsi1_588, osg0_432, osg0_434, osg1_432, osg1_434, osh_606, \
                         osh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_8 * osg0_432[k]
                   - f_9 * osg1_432[k]
                   + f_3 * pc_z[k] * osh_606[k];

        t_810[k] = f_20 * nsh_461[k]
                   + f_3 * pc_y[k] * osh_608[k];

        t_811[k] = f_1 * osg0_434[k]
                   - f_2 * osg1_434[k]
                   + f_3 * pc_z[k] * osh_608[k];

        t_812[k] = pa_z[k] * nsi0_588[k]
                   - f_10 * pc_z[k] * nsi1_588[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, pa_z, pc_y, pc_z, nsi0_591, nsh_441, \
                         nsh_462, nsh_464, nsi1_591, osh_609, osh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_21 * nsh_462[k]
                   + f_3 * pc_y[k] * osh_609[k];

        t_814[k] = f_11 * nsh_441[k]
                   + f_3 * pc_z[k] * osh_609[k];

        t_815[k] = pa_z[k] * nsi0_591[k]
                   - f_10 * pc_z[k] * nsi1_591[k];

        t_816[k] = f_21 * nsh_464[k]
                   + f_3 * pc_y[k] * osh_611[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pa_z, pc_x, pc_z, nsi0_594, nsh_444, nsh_614, \
                         nsi1_594, osg0_440, osg1_440, osh_612, \
                         osh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_14 * nsh_614[k]
                   + f_8 * osg0_440[k]
                   - f_9 * osg1_440[k]
                   + f_3 * pc_x[k] * osh_614[k];

        t_818[k] = pa_z[k] * nsi0_594[k]
                   - f_10 * pc_z[k] * nsi1_594[k];

        t_819[k] = f_11 * nsh_444[k]
                   + f_3 * pc_z[k] * osh_612[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_598 = buffer.data(nsi0 + 598);
    const auto *nsi0_600 = buffer.data(nsi0 + 600);
    const auto *nsi0_609 = buffer.data(nsi0 + 609);

    const auto *nsh_447 = buffer.data(nsh + 447);
    const auto *nsh_448 = buffer.data(nsh + 448);
    const auto *nsh_456 = buffer.data(nsh + 456);
    const auto *nsh_461 = buffer.data(nsh + 461);
    const auto *nsh_462 = buffer.data(nsh + 462);
    const auto *nsh_465 = buffer.data(nsh + 465);
    const auto *nsh_467 = buffer.data(nsh + 467);
    const auto *nsh_468 = buffer.data(nsh + 468);
    const auto *nsh_471 = buffer.data(nsh + 471);
    const auto *nsh_477 = buffer.data(nsh + 477);
    const auto *nsh_479 = buffer.data(nsh + 479);
    const auto *nsh_480 = buffer.data(nsh + 480);
    const auto *nsh_481 = buffer.data(nsh + 481);
    const auto *nsh_482 = buffer.data(nsh + 482);
    const auto *nsh_483 = buffer.data(nsh + 483);
    const auto *nsh_485 = buffer.data(nsh + 485);
    const auto *nsh_486 = buffer.data(nsh + 486);
    const auto *nsh_488 = buffer.data(nsh + 488);
    const auto *nsh_489 = buffer.data(nsh + 489);
    const auto *nsh_492 = buffer.data(nsh + 492);
    const auto *nsh_498 = buffer.data(nsh + 498);
    const auto *nsh_500 = buffer.data(nsh + 500);
    const auto *nsh_501 = buffer.data(nsh + 501);
    const auto *nsh_502 = buffer.data(nsh + 502);
    const auto *nsh_503 = buffer.data(nsh + 503);
    const auto *nsh_504 = buffer.data(nsh + 504);
    const auto *nsh_506 = buffer.data(nsh + 506);
    const auto *nsh_507 = buffer.data(nsh + 507);
    const auto *nsh_509 = buffer.data(nsh + 509);
    const auto *nsh_510 = buffer.data(nsh + 510);
    const auto *nsh_513 = buffer.data(nsh + 513);
    const auto *nsh_519 = buffer.data(nsh + 519);
    const auto *nsh_521 = buffer.data(nsh + 521);
    const auto *nsh_522 = buffer.data(nsh + 522);
    const auto *nsh_523 = buffer.data(nsh + 523);
    const auto *nsh_524 = buffer.data(nsh + 524);
    const auto *nsh_525 = buffer.data(nsh + 525);
    const auto *nsh_527 = buffer.data(nsh + 527);
    const auto *nsh_530 = buffer.data(nsh + 530);
    const auto *nsh_534 = buffer.data(nsh + 534);
    const auto *nsh_540 = buffer.data(nsh + 540);
    const auto *nsh_542 = buffer.data(nsh + 542);
    const auto *nsh_543 = buffer.data(nsh + 543);
    const auto *nsh_544 = buffer.data(nsh + 544);
    const auto *nsh_545 = buffer.data(nsh + 545);
    const auto *nsh_618 = buffer.data(nsh + 618);
    const auto *nsh_623 = buffer.data(nsh + 623);
    const auto *nsh_624 = buffer.data(nsh + 624);
    const auto *nsh_625 = buffer.data(nsh + 625);
    const auto *nsh_626 = buffer.data(nsh + 626);
    const auto *nsh_627 = buffer.data(nsh + 627);
    const auto *nsh_628 = buffer.data(nsh + 628);
    const auto *nsh_629 = buffer.data(nsh + 629);
    const auto *nsh_630 = buffer.data(nsh + 630);
    const auto *nsh_633 = buffer.data(nsh + 633);
    const auto *nsh_635 = buffer.data(nsh + 635);
    const auto *nsh_636 = buffer.data(nsh + 636);
    const auto *nsh_639 = buffer.data(nsh + 639);
    const auto *nsh_640 = buffer.data(nsh + 640);
    const auto *nsh_642 = buffer.data(nsh + 642);
    const auto *nsh_644 = buffer.data(nsh + 644);
    const auto *nsh_645 = buffer.data(nsh + 645);
    const auto *nsh_646 = buffer.data(nsh + 646);
    const auto *nsh_647 = buffer.data(nsh + 647);
    const auto *nsh_648 = buffer.data(nsh + 648);
    const auto *nsh_649 = buffer.data(nsh + 649);
    const auto *nsh_650 = buffer.data(nsh + 650);
    const auto *nsh_651 = buffer.data(nsh + 651);
    const auto *nsh_654 = buffer.data(nsh + 654);
    const auto *nsh_656 = buffer.data(nsh + 656);
    const auto *nsh_657 = buffer.data(nsh + 657);
    const auto *nsh_660 = buffer.data(nsh + 660);
    const auto *nsh_661 = buffer.data(nsh + 661);
    const auto *nsh_663 = buffer.data(nsh + 663);
    const auto *nsh_665 = buffer.data(nsh + 665);
    const auto *nsh_666 = buffer.data(nsh + 666);
    const auto *nsh_667 = buffer.data(nsh + 667);
    const auto *nsh_668 = buffer.data(nsh + 668);
    const auto *nsh_669 = buffer.data(nsh + 669);
    const auto *nsh_670 = buffer.data(nsh + 670);
    const auto *nsh_671 = buffer.data(nsh + 671);
    const auto *nsh_672 = buffer.data(nsh + 672);
    const auto *nsh_675 = buffer.data(nsh + 675);
    const auto *nsh_677 = buffer.data(nsh + 677);
    const auto *nsh_678 = buffer.data(nsh + 678);
    const auto *nsh_681 = buffer.data(nsh + 681);
    const auto *nsh_682 = buffer.data(nsh + 682);
    const auto *nsh_684 = buffer.data(nsh + 684);
    const auto *nsh_686 = buffer.data(nsh + 686);
    const auto *nsh_687 = buffer.data(nsh + 687);
    const auto *nsh_688 = buffer.data(nsh + 688);
    const auto *nsh_689 = buffer.data(nsh + 689);
    const auto *nsh_690 = buffer.data(nsh + 690);
    const auto *nsh_691 = buffer.data(nsh + 691);
    const auto *nsh_692 = buffer.data(nsh + 692);
    const auto *nsh_693 = buffer.data(nsh + 693);

    const auto *nsi1_598 = buffer.data(nsi1 + 598);
    const auto *nsi1_600 = buffer.data(nsi1 + 600);
    const auto *nsi1_609 = buffer.data(nsi1 + 609);

    const auto *osg0_444 = buffer.data(osg0 + 444);
    const auto *osg0_447 = buffer.data(osg0 + 447);
    const auto *osg0_448 = buffer.data(osg0 + 448);
    const auto *osg0_449 = buffer.data(osg0 + 449);
    const auto *osg0_450 = buffer.data(osg0 + 450);
    const auto *osg0_453 = buffer.data(osg0 + 453);
    const auto *osg0_455 = buffer.data(osg0 + 455);
    const auto *osg0_456 = buffer.data(osg0 + 456);
    const auto *osg0_459 = buffer.data(osg0 + 459);
    const auto *osg0_460 = buffer.data(osg0 + 460);
    const auto *osg0_462 = buffer.data(osg0 + 462);
    const auto *osg0_463 = buffer.data(osg0 + 463);
    const auto *osg0_464 = buffer.data(osg0 + 464);
    const auto *osg0_465 = buffer.data(osg0 + 465);
    const auto *osg0_468 = buffer.data(osg0 + 468);
    const auto *osg0_470 = buffer.data(osg0 + 470);
    const auto *osg0_471 = buffer.data(osg0 + 471);
    const auto *osg0_474 = buffer.data(osg0 + 474);
    const auto *osg0_475 = buffer.data(osg0 + 475);
    const auto *osg0_477 = buffer.data(osg0 + 477);
    const auto *osg0_478 = buffer.data(osg0 + 478);
    const auto *osg0_479 = buffer.data(osg0 + 479);
    const auto *osg0_480 = buffer.data(osg0 + 480);
    const auto *osg0_483 = buffer.data(osg0 + 483);
    const auto *osg0_485 = buffer.data(osg0 + 485);
    const auto *osg0_486 = buffer.data(osg0 + 486);
    const auto *osg0_489 = buffer.data(osg0 + 489);
    const auto *osg0_490 = buffer.data(osg0 + 490);
    const auto *osg0_492 = buffer.data(osg0 + 492);
    const auto *osg0_493 = buffer.data(osg0 + 493);
    const auto *osg0_494 = buffer.data(osg0 + 494);
    const auto *osg0_495 = buffer.data(osg0 + 495);

    const auto *osg1_444 = buffer.data(osg1 + 444);
    const auto *osg1_447 = buffer.data(osg1 + 447);
    const auto *osg1_448 = buffer.data(osg1 + 448);
    const auto *osg1_449 = buffer.data(osg1 + 449);
    const auto *osg1_450 = buffer.data(osg1 + 450);
    const auto *osg1_453 = buffer.data(osg1 + 453);
    const auto *osg1_455 = buffer.data(osg1 + 455);
    const auto *osg1_456 = buffer.data(osg1 + 456);
    const auto *osg1_459 = buffer.data(osg1 + 459);
    const auto *osg1_460 = buffer.data(osg1 + 460);
    const auto *osg1_462 = buffer.data(osg1 + 462);
    const auto *osg1_463 = buffer.data(osg1 + 463);
    const auto *osg1_464 = buffer.data(osg1 + 464);
    const auto *osg1_465 = buffer.data(osg1 + 465);
    const auto *osg1_468 = buffer.data(osg1 + 468);
    const auto *osg1_470 = buffer.data(osg1 + 470);
    const auto *osg1_471 = buffer.data(osg1 + 471);
    const auto *osg1_474 = buffer.data(osg1 + 474);
    const auto *osg1_475 = buffer.data(osg1 + 475);
    const auto *osg1_477 = buffer.data(osg1 + 477);
    const auto *osg1_478 = buffer.data(osg1 + 478);
    const auto *osg1_479 = buffer.data(osg1 + 479);
    const auto *osg1_480 = buffer.data(osg1 + 480);
    const auto *osg1_483 = buffer.data(osg1 + 483);
    const auto *osg1_485 = buffer.data(osg1 + 485);
    const auto *osg1_486 = buffer.data(osg1 + 486);
    const auto *osg1_489 = buffer.data(osg1 + 489);
    const auto *osg1_490 = buffer.data(osg1 + 490);
    const auto *osg1_492 = buffer.data(osg1 + 492);
    const auto *osg1_493 = buffer.data(osg1 + 493);
    const auto *osg1_494 = buffer.data(osg1 + 494);
    const auto *osg1_495 = buffer.data(osg1 + 495);

    const auto *osh_614 = buffer.data(osh + 614);
    const auto *osh_615 = buffer.data(osh + 615);
    const auto *osh_618 = buffer.data(osh + 618);
    const auto *osh_623 = buffer.data(osh + 623);
    const auto *osh_624 = buffer.data(osh + 624);
    const auto *osh_625 = buffer.data(osh + 625);
    const auto *osh_626 = buffer.data(osh + 626);
    const auto *osh_627 = buffer.data(osh + 627);
    const auto *osh_628 = buffer.data(osh + 628);
    const auto *osh_629 = buffer.data(osh + 629);
    const auto *osh_630 = buffer.data(osh + 630);
    const auto *osh_632 = buffer.data(osh + 632);
    const auto *osh_633 = buffer.data(osh + 633);
    const auto *osh_635 = buffer.data(osh + 635);
    const auto *osh_636 = buffer.data(osh + 636);
    const auto *osh_639 = buffer.data(osh + 639);
    const auto *osh_640 = buffer.data(osh + 640);
    const auto *osh_642 = buffer.data(osh + 642);
    const auto *osh_644 = buffer.data(osh + 644);
    const auto *osh_645 = buffer.data(osh + 645);
    const auto *osh_646 = buffer.data(osh + 646);
    const auto *osh_647 = buffer.data(osh + 647);
    const auto *osh_648 = buffer.data(osh + 648);
    const auto *osh_649 = buffer.data(osh + 649);
    const auto *osh_650 = buffer.data(osh + 650);
    const auto *osh_651 = buffer.data(osh + 651);
    const auto *osh_653 = buffer.data(osh + 653);
    const auto *osh_654 = buffer.data(osh + 654);
    const auto *osh_656 = buffer.data(osh + 656);
    const auto *osh_657 = buffer.data(osh + 657);
    const auto *osh_660 = buffer.data(osh + 660);
    const auto *osh_661 = buffer.data(osh + 661);
    const auto *osh_663 = buffer.data(osh + 663);
    const auto *osh_665 = buffer.data(osh + 665);
    const auto *osh_666 = buffer.data(osh + 666);
    const auto *osh_667 = buffer.data(osh + 667);
    const auto *osh_668 = buffer.data(osh + 668);
    const auto *osh_669 = buffer.data(osh + 669);
    const auto *osh_670 = buffer.data(osh + 670);
    const auto *osh_671 = buffer.data(osh + 671);
    const auto *osh_672 = buffer.data(osh + 672);
    const auto *osh_674 = buffer.data(osh + 674);
    const auto *osh_675 = buffer.data(osh + 675);
    const auto *osh_677 = buffer.data(osh + 677);
    const auto *osh_678 = buffer.data(osh + 678);
    const auto *osh_681 = buffer.data(osh + 681);
    const auto *osh_682 = buffer.data(osh + 682);
    const auto *osh_684 = buffer.data(osh + 684);
    const auto *osh_686 = buffer.data(osh + 686);
    const auto *osh_687 = buffer.data(osh + 687);
    const auto *osh_688 = buffer.data(osh + 688);
    const auto *osh_689 = buffer.data(osh + 689);
    const auto *osh_690 = buffer.data(osh + 690);
    const auto *osh_691 = buffer.data(osh + 691);
    const auto *osh_692 = buffer.data(osh + 692);
    const auto *osh_693 = buffer.data(osh + 693);

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_x, pc_y, pc_z, nsi0_598, nsh_467, \
                         nsh_618, nsi1_598, osg0_444, osg1_444, osh_614, \
                         osh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_21 * nsh_467[k]
                   + f_3 * pc_y[k] * osh_614[k];

        t_821[k] = f_14 * nsh_618[k]
                   + f_6 * osg0_444[k]
                   - f_7 * osg1_444[k]
                   + f_3 * pc_x[k] * osh_618[k];

        t_822[k] = pa_z[k] * nsi0_598[k]
                   - f_10 * pc_z[k] * nsi1_598[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pa_z, pc_y, pc_z, nsi0_600, nsh_447, nsh_448, \
                         nsh_471, nsi1_600, osh_615, osh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_11 * nsh_447[k]
                   + f_3 * pc_z[k] * osh_615[k];

        t_824[k] = pa_z[k] * nsi0_600[k]
                   + f_12 * nsh_448[k]
                   - f_10 * pc_z[k] * nsi1_600[k];

        t_825[k] = f_21 * nsh_471[k]
                   + f_3 * pc_y[k] * osh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, nsh_623, nsh_624, nsh_625, nsh_626, \
                         osg0_449, osg1_449, osh_623, osh_624, osh_625, \
                         osh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_14 * nsh_623[k]
                   + f_4 * osg0_449[k]
                   - f_5 * osg1_449[k]
                   + f_3 * pc_x[k] * osh_623[k];

        t_827[k] = f_14 * nsh_624[k]
                   + f_3 * pc_x[k] * osh_624[k];

        t_828[k] = f_14 * nsh_625[k]
                   + f_3 * pc_x[k] * osh_625[k];

        t_829[k] = f_14 * nsh_626[k]
                   + f_3 * pc_x[k] * osh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pa_z, pc_x, pc_z, nsi0_609, nsh_627, \
                         nsh_628, nsh_629, nsi1_609, osh_627, osh_628, \
                         osh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_14 * nsh_627[k]
                   + f_3 * pc_x[k] * osh_627[k];

        t_831[k] = f_14 * nsh_628[k]
                   + f_3 * pc_x[k] * osh_628[k];

        t_832[k] = f_14 * nsh_629[k]
                   + f_3 * pc_x[k] * osh_629[k];

        t_833[k] = pa_z[k] * nsi0_609[k]
                   - f_10 * pc_z[k] * nsi1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, nsh_456, nsh_479, nsh_480, osg0_447, \
                         osg0_448, osg1_447, osg1_448, osh_624, osh_626, \
                         osh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * nsh_456[k]
                   + f_3 * pc_z[k] * osh_624[k];

        t_835[k] = f_21 * nsh_479[k]
                   + f_8 * osg0_447[k]
                   - f_9 * osg1_447[k]
                   + f_3 * pc_y[k] * osh_626[k];

        t_836[k] = f_21 * nsh_480[k]
                   + f_6 * osg0_448[k]
                   - f_7 * osg1_448[k]
                   + f_3 * pc_y[k] * osh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, nsh_461, nsh_481, nsh_482, osg0_449, \
                         osg1_449, osh_628, osh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_21 * nsh_481[k]
                   + f_4 * osg0_449[k]
                   - f_5 * osg1_449[k]
                   + f_3 * pc_y[k] * osh_628[k];

        t_838[k] = f_21 * nsh_482[k]
                   + f_3 * pc_y[k] * osh_629[k];

        t_839[k] = f_11 * nsh_461[k]
                   + f_1 * osg0_449[k]
                   - f_2 * osg1_449[k]
                   + f_3 * pc_z[k] * osh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, nsh_462, nsh_483, nsh_630, \
                         osg0_450, osg1_450, osh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_14 * nsh_630[k]
                   + f_1 * osg0_450[k]
                   - f_2 * osg1_450[k]
                   + f_3 * pc_x[k] * osh_630[k];

        t_841[k] = f_22 * nsh_483[k]
                   + f_3 * pc_y[k] * osh_630[k];

        t_842[k] = f_12 * nsh_462[k]
                   + f_3 * pc_z[k] * osh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, nsh_485, nsh_633, nsh_635, osg0_453, \
                         osg0_455, osg1_453, osg1_455, osh_632, osh_633, \
                         osh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_14 * nsh_633[k]
                   + f_8 * osg0_453[k]
                   - f_9 * osg1_453[k]
                   + f_3 * pc_x[k] * osh_633[k];

        t_844[k] = f_22 * nsh_485[k]
                   + f_3 * pc_y[k] * osh_632[k];

        t_845[k] = f_14 * nsh_635[k]
                   + f_8 * osg0_455[k]
                   - f_9 * osg1_455[k]
                   + f_3 * pc_x[k] * osh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, nsh_465, nsh_488, nsh_636, \
                         osg0_456, osg1_456, osh_633, osh_635, \
                         osh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_14 * nsh_636[k]
                   + f_6 * osg0_456[k]
                   - f_7 * osg1_456[k]
                   + f_3 * pc_x[k] * osh_636[k];

        t_847[k] = f_12 * nsh_465[k]
                   + f_3 * pc_z[k] * osh_633[k];

        t_848[k] = f_22 * nsh_488[k]
                   + f_3 * pc_y[k] * osh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, nsh_468, nsh_639, nsh_640, osg0_459, \
                         osg0_460, osg1_459, osg1_460, osh_636, osh_639, \
                         osh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_14 * nsh_639[k]
                   + f_6 * osg0_459[k]
                   - f_7 * osg1_459[k]
                   + f_3 * pc_x[k] * osh_639[k];

        t_850[k] = f_14 * nsh_640[k]
                   + f_4 * osg0_460[k]
                   - f_5 * osg1_460[k]
                   + f_3 * pc_x[k] * osh_640[k];

        t_851[k] = f_12 * nsh_468[k]
                   + f_3 * pc_z[k] * osh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, nsh_492, nsh_642, nsh_644, osg0_462, \
                         osg0_464, osg1_462, osg1_464, osh_639, osh_642, \
                         osh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_14 * nsh_642[k]
                   + f_4 * osg0_462[k]
                   - f_5 * osg1_462[k]
                   + f_3 * pc_x[k] * osh_642[k];

        t_853[k] = f_22 * nsh_492[k]
                   + f_3 * pc_y[k] * osh_639[k];

        t_854[k] = f_14 * nsh_644[k]
                   + f_4 * osg0_464[k]
                   - f_5 * osg1_464[k]
                   + f_3 * pc_x[k] * osh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, nsh_645, nsh_646, nsh_647, \
                         nsh_648, nsh_649, osh_645, osh_646, osh_647, osh_648, \
                         osh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_14 * nsh_645[k]
                   + f_3 * pc_x[k] * osh_645[k];

        t_856[k] = f_14 * nsh_646[k]
                   + f_3 * pc_x[k] * osh_646[k];

        t_857[k] = f_14 * nsh_647[k]
                   + f_3 * pc_x[k] * osh_647[k];

        t_858[k] = f_14 * nsh_648[k]
                   + f_3 * pc_x[k] * osh_648[k];

        t_859[k] = f_14 * nsh_649[k]
                   + f_3 * pc_x[k] * osh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, nsh_477, nsh_498, nsh_650, \
                         osg0_460, osg1_460, osh_645, osh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_14 * nsh_650[k]
                   + f_3 * pc_x[k] * osh_650[k];

        t_861[k] = f_22 * nsh_498[k]
                   + f_1 * osg0_460[k]
                   - f_2 * osg1_460[k]
                   + f_3 * pc_y[k] * osh_645[k];

        t_862[k] = f_12 * nsh_477[k]
                   + f_3 * pc_z[k] * osh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, nsh_500, nsh_501, nsh_502, osg0_462, \
                         osg0_463, osg0_464, osg1_462, osg1_463, osg1_464, osh_647, osh_648, \
                         osh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_22 * nsh_500[k]
                   + f_8 * osg0_462[k]
                   - f_9 * osg1_462[k]
                   + f_3 * pc_y[k] * osh_647[k];

        t_864[k] = f_22 * nsh_501[k]
                   + f_6 * osg0_463[k]
                   - f_7 * osg1_463[k]
                   + f_3 * pc_y[k] * osh_648[k];

        t_865[k] = f_22 * nsh_502[k]
                   + f_4 * osg0_464[k]
                   - f_5 * osg1_464[k]
                   + f_3 * pc_y[k] * osh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, nsh_482, nsh_503, nsh_651, \
                         osg0_464, osg0_465, osg1_464, osg1_465, osh_650, \
                         osh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_22 * nsh_503[k]
                   + f_3 * pc_y[k] * osh_650[k];

        t_867[k] = f_12 * nsh_482[k]
                   + f_1 * osg0_464[k]
                   - f_2 * osg1_464[k]
                   + f_3 * pc_z[k] * osh_650[k];

        t_868[k] = f_14 * nsh_651[k]
                   + f_1 * osg0_465[k]
                   - f_2 * osg1_465[k]
                   + f_3 * pc_x[k] * osh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, nsh_483, nsh_504, \
                         nsh_506, nsh_654, osg0_468, osg1_468, osh_651, osh_653, \
                         osh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * nsh_504[k]
                   + f_3 * pc_y[k] * osh_651[k];

        t_870[k] = f_13 * nsh_483[k]
                   + f_3 * pc_z[k] * osh_651[k];

        t_871[k] = f_14 * nsh_654[k]
                   + f_8 * osg0_468[k]
                   - f_9 * osg1_468[k]
                   + f_3 * pc_x[k] * osh_654[k];

        t_872[k] = f_14 * nsh_506[k]
                   + f_3 * pc_y[k] * osh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, nsh_486, nsh_656, nsh_657, osg0_470, \
                         osg0_471, osg1_470, osg1_471, osh_654, osh_656, \
                         osh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_14 * nsh_656[k]
                   + f_8 * osg0_470[k]
                   - f_9 * osg1_470[k]
                   + f_3 * pc_x[k] * osh_656[k];

        t_874[k] = f_14 * nsh_657[k]
                   + f_6 * osg0_471[k]
                   - f_7 * osg1_471[k]
                   + f_3 * pc_x[k] * osh_657[k];

        t_875[k] = f_13 * nsh_486[k]
                   + f_3 * pc_z[k] * osh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, nsh_509, nsh_660, nsh_661, osg0_474, \
                         osg0_475, osg1_474, osg1_475, osh_656, osh_660, \
                         osh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * nsh_509[k]
                   + f_3 * pc_y[k] * osh_656[k];

        t_877[k] = f_14 * nsh_660[k]
                   + f_6 * osg0_474[k]
                   - f_7 * osg1_474[k]
                   + f_3 * pc_x[k] * osh_660[k];

        t_878[k] = f_14 * nsh_661[k]
                   + f_4 * osg0_475[k]
                   - f_5 * osg1_475[k]
                   + f_3 * pc_x[k] * osh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, nsh_489, nsh_513, nsh_663, \
                         osg0_477, osg1_477, osh_657, osh_660, \
                         osh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * nsh_489[k]
                   + f_3 * pc_z[k] * osh_657[k];

        t_880[k] = f_14 * nsh_663[k]
                   + f_4 * osg0_477[k]
                   - f_5 * osg1_477[k]
                   + f_3 * pc_x[k] * osh_663[k];

        t_881[k] = f_14 * nsh_513[k]
                   + f_3 * pc_y[k] * osh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, nsh_665, nsh_666, nsh_667, nsh_668, \
                         osg0_479, osg1_479, osh_665, osh_666, osh_667, \
                         osh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_14 * nsh_665[k]
                   + f_4 * osg0_479[k]
                   - f_5 * osg1_479[k]
                   + f_3 * pc_x[k] * osh_665[k];

        t_883[k] = f_14 * nsh_666[k]
                   + f_3 * pc_x[k] * osh_666[k];

        t_884[k] = f_14 * nsh_667[k]
                   + f_3 * pc_x[k] * osh_667[k];

        t_885[k] = f_14 * nsh_668[k]
                   + f_3 * pc_x[k] * osh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, nsh_519, nsh_669, nsh_670, \
                         nsh_671, osg0_475, osg1_475, osh_666, osh_669, osh_670, \
                         osh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_14 * nsh_669[k]
                   + f_3 * pc_x[k] * osh_669[k];

        t_887[k] = f_14 * nsh_670[k]
                   + f_3 * pc_x[k] * osh_670[k];

        t_888[k] = f_14 * nsh_671[k]
                   + f_3 * pc_x[k] * osh_671[k];

        t_889[k] = f_14 * nsh_519[k]
                   + f_1 * osg0_475[k]
                   - f_2 * osg1_475[k]
                   + f_3 * pc_y[k] * osh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, nsh_498, nsh_521, nsh_522, osg0_477, \
                         osg0_478, osg1_477, osg1_478, osh_666, osh_668, \
                         osh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * nsh_498[k]
                   + f_3 * pc_z[k] * osh_666[k];

        t_891[k] = f_14 * nsh_521[k]
                   + f_8 * osg0_477[k]
                   - f_9 * osg1_477[k]
                   + f_3 * pc_y[k] * osh_668[k];

        t_892[k] = f_14 * nsh_522[k]
                   + f_6 * osg0_478[k]
                   - f_7 * osg1_478[k]
                   + f_3 * pc_y[k] * osh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, nsh_503, nsh_523, nsh_524, osg0_479, \
                         osg1_479, osh_670, osh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * nsh_523[k]
                   + f_4 * osg0_479[k]
                   - f_5 * osg1_479[k]
                   + f_3 * pc_y[k] * osh_670[k];

        t_894[k] = f_14 * nsh_524[k]
                   + f_3 * pc_y[k] * osh_671[k];

        t_895[k] = f_13 * nsh_503[k]
                   + f_1 * osg0_479[k]
                   - f_2 * osg1_479[k]
                   + f_3 * pc_z[k] * osh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, nsh_504, nsh_525, nsh_672, \
                         osg0_480, osg1_480, osh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_14 * nsh_672[k]
                   + f_1 * osg0_480[k]
                   - f_2 * osg1_480[k]
                   + f_3 * pc_x[k] * osh_672[k];

        t_897[k] = f_13 * nsh_525[k]
                   + f_3 * pc_y[k] * osh_672[k];

        t_898[k] = f_14 * nsh_504[k]
                   + f_3 * pc_z[k] * osh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, nsh_527, nsh_675, nsh_677, osg0_483, \
                         osg0_485, osg1_483, osg1_485, osh_674, osh_675, \
                         osh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_14 * nsh_675[k]
                   + f_8 * osg0_483[k]
                   - f_9 * osg1_483[k]
                   + f_3 * pc_x[k] * osh_675[k];

        t_900[k] = f_13 * nsh_527[k]
                   + f_3 * pc_y[k] * osh_674[k];

        t_901[k] = f_14 * nsh_677[k]
                   + f_8 * osg0_485[k]
                   - f_9 * osg1_485[k]
                   + f_3 * pc_x[k] * osh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, nsh_507, nsh_530, nsh_678, \
                         osg0_486, osg1_486, osh_675, osh_677, \
                         osh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_14 * nsh_678[k]
                   + f_6 * osg0_486[k]
                   - f_7 * osg1_486[k]
                   + f_3 * pc_x[k] * osh_678[k];

        t_903[k] = f_14 * nsh_507[k]
                   + f_3 * pc_z[k] * osh_675[k];

        t_904[k] = f_13 * nsh_530[k]
                   + f_3 * pc_y[k] * osh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, nsh_510, nsh_681, nsh_682, osg0_489, \
                         osg0_490, osg1_489, osg1_490, osh_678, osh_681, \
                         osh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_14 * nsh_681[k]
                   + f_6 * osg0_489[k]
                   - f_7 * osg1_489[k]
                   + f_3 * pc_x[k] * osh_681[k];

        t_906[k] = f_14 * nsh_682[k]
                   + f_4 * osg0_490[k]
                   - f_5 * osg1_490[k]
                   + f_3 * pc_x[k] * osh_682[k];

        t_907[k] = f_14 * nsh_510[k]
                   + f_3 * pc_z[k] * osh_678[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, nsh_534, nsh_684, nsh_686, osg0_492, \
                         osg0_494, osg1_492, osg1_494, osh_681, osh_684, \
                         osh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * nsh_684[k]
                   + f_4 * osg0_492[k]
                   - f_5 * osg1_492[k]
                   + f_3 * pc_x[k] * osh_684[k];

        t_909[k] = f_13 * nsh_534[k]
                   + f_3 * pc_y[k] * osh_681[k];

        t_910[k] = f_14 * nsh_686[k]
                   + f_4 * osg0_494[k]
                   - f_5 * osg1_494[k]
                   + f_3 * pc_x[k] * osh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, nsh_687, nsh_688, nsh_689, \
                         nsh_690, nsh_691, osh_687, osh_688, osh_689, osh_690, \
                         osh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_14 * nsh_687[k]
                   + f_3 * pc_x[k] * osh_687[k];

        t_912[k] = f_14 * nsh_688[k]
                   + f_3 * pc_x[k] * osh_688[k];

        t_913[k] = f_14 * nsh_689[k]
                   + f_3 * pc_x[k] * osh_689[k];

        t_914[k] = f_14 * nsh_690[k]
                   + f_3 * pc_x[k] * osh_690[k];

        t_915[k] = f_14 * nsh_691[k]
                   + f_3 * pc_x[k] * osh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, nsh_519, nsh_540, nsh_692, \
                         osg0_490, osg1_490, osh_687, osh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_14 * nsh_692[k]
                   + f_3 * pc_x[k] * osh_692[k];

        t_917[k] = f_13 * nsh_540[k]
                   + f_1 * osg0_490[k]
                   - f_2 * osg1_490[k]
                   + f_3 * pc_y[k] * osh_687[k];

        t_918[k] = f_14 * nsh_519[k]
                   + f_3 * pc_z[k] * osh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, nsh_542, nsh_543, nsh_544, osg0_492, \
                         osg0_493, osg0_494, osg1_492, osg1_493, osg1_494, osh_689, osh_690, \
                         osh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * nsh_542[k]
                   + f_8 * osg0_492[k]
                   - f_9 * osg1_492[k]
                   + f_3 * pc_y[k] * osh_689[k];

        t_920[k] = f_13 * nsh_543[k]
                   + f_6 * osg0_493[k]
                   - f_7 * osg1_493[k]
                   + f_3 * pc_y[k] * osh_690[k];

        t_921[k] = f_13 * nsh_544[k]
                   + f_4 * osg0_494[k]
                   - f_5 * osg1_494[k]
                   + f_3 * pc_y[k] * osh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, nsh_524, nsh_545, nsh_693, \
                         osg0_494, osg0_495, osg1_494, osg1_495, osh_692, \
                         osh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * nsh_545[k]
                   + f_3 * pc_y[k] * osh_692[k];

        t_923[k] = f_14 * nsh_524[k]
                   + f_1 * osg0_494[k]
                   - f_2 * osg1_494[k]
                   + f_3 * pc_z[k] * osh_692[k];

        t_924[k] = f_14 * nsh_693[k]
                   + f_1 * osg0_495[k]
                   - f_2 * osg1_495[k]
                   + f_3 * pc_x[k] * osh_693[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_756 = buffer.data(nsi0 + 756);
    const auto *nsi0_759 = buffer.data(nsi0 + 759);
    const auto *nsi0_761 = buffer.data(nsi0 + 761);
    const auto *nsi0_762 = buffer.data(nsi0 + 762);
    const auto *nsi0_765 = buffer.data(nsi0 + 765);
    const auto *nsi0_766 = buffer.data(nsi0 + 766);
    const auto *nsi0_768 = buffer.data(nsi0 + 768);
    const auto *nsi0_770 = buffer.data(nsi0 + 770);
    const auto *nsi0_783 = buffer.data(nsi0 + 783);
    const auto *nsi0_784 = buffer.data(nsi0 + 784);
    const auto *nsi0_787 = buffer.data(nsi0 + 787);
    const auto *nsi0_790 = buffer.data(nsi0 + 790);

    const auto *nsh_525 = buffer.data(nsh + 525);
    const auto *nsh_528 = buffer.data(nsh + 528);
    const auto *nsh_531 = buffer.data(nsh + 531);
    const auto *nsh_540 = buffer.data(nsh + 540);
    const auto *nsh_545 = buffer.data(nsh + 545);
    const auto *nsh_546 = buffer.data(nsh + 546);
    const auto *nsh_548 = buffer.data(nsh + 548);
    const auto *nsh_549 = buffer.data(nsh + 549);
    const auto *nsh_551 = buffer.data(nsh + 551);
    const auto *nsh_552 = buffer.data(nsh + 552);
    const auto *nsh_555 = buffer.data(nsh + 555);
    const auto *nsh_561 = buffer.data(nsh + 561);
    const auto *nsh_563 = buffer.data(nsh + 563);
    const auto *nsh_564 = buffer.data(nsh + 564);
    const auto *nsh_565 = buffer.data(nsh + 565);
    const auto *nsh_566 = buffer.data(nsh + 566);
    const auto *nsh_567 = buffer.data(nsh + 567);
    const auto *nsh_568 = buffer.data(nsh + 568);
    const auto *nsh_569 = buffer.data(nsh + 569);
    const auto *nsh_570 = buffer.data(nsh + 570);
    const auto *nsh_572 = buffer.data(nsh + 572);
    const auto *nsh_573 = buffer.data(nsh + 573);
    const auto *nsh_575 = buffer.data(nsh + 575);
    const auto *nsh_576 = buffer.data(nsh + 576);
    const auto *nsh_582 = buffer.data(nsh + 582);
    const auto *nsh_584 = buffer.data(nsh + 584);
    const auto *nsh_585 = buffer.data(nsh + 585);
    const auto *nsh_586 = buffer.data(nsh + 586);
    const auto *nsh_587 = buffer.data(nsh + 587);
    const auto *nsh_588 = buffer.data(nsh + 588);
    const auto *nsh_591 = buffer.data(nsh + 591);
    const auto *nsh_593 = buffer.data(nsh + 593);
    const auto *nsh_597 = buffer.data(nsh + 597);
    const auto *nsh_603 = buffer.data(nsh + 603);
    const auto *nsh_608 = buffer.data(nsh + 608);
    const auto *nsh_609 = buffer.data(nsh + 609);
    const auto *nsh_611 = buffer.data(nsh + 611);
    const auto *nsh_696 = buffer.data(nsh + 696);
    const auto *nsh_698 = buffer.data(nsh + 698);
    const auto *nsh_699 = buffer.data(nsh + 699);
    const auto *nsh_702 = buffer.data(nsh + 702);
    const auto *nsh_703 = buffer.data(nsh + 703);
    const auto *nsh_705 = buffer.data(nsh + 705);
    const auto *nsh_707 = buffer.data(nsh + 707);
    const auto *nsh_708 = buffer.data(nsh + 708);
    const auto *nsh_709 = buffer.data(nsh + 709);
    const auto *nsh_710 = buffer.data(nsh + 710);
    const auto *nsh_711 = buffer.data(nsh + 711);
    const auto *nsh_712 = buffer.data(nsh + 712);
    const auto *nsh_713 = buffer.data(nsh + 713);
    const auto *nsh_729 = buffer.data(nsh + 729);
    const auto *nsh_730 = buffer.data(nsh + 730);
    const auto *nsh_731 = buffer.data(nsh + 731);
    const auto *nsh_732 = buffer.data(nsh + 732);
    const auto *nsh_733 = buffer.data(nsh + 733);
    const auto *nsh_734 = buffer.data(nsh + 734);
    const auto *nsh_735 = buffer.data(nsh + 735);
    const auto *nsh_740 = buffer.data(nsh + 740);
    const auto *nsh_744 = buffer.data(nsh + 744);
    const auto *nsh_749 = buffer.data(nsh + 749);
    const auto *nsh_750 = buffer.data(nsh + 750);
    const auto *nsh_751 = buffer.data(nsh + 751);
    const auto *nsh_752 = buffer.data(nsh + 752);
    const auto *nsh_753 = buffer.data(nsh + 753);
    const auto *nsh_755 = buffer.data(nsh + 755);
    const auto *nsh_756 = buffer.data(nsh + 756);
    const auto *nsh_759 = buffer.data(nsh + 759);
    const auto *nsh_762 = buffer.data(nsh + 762);
    const auto *nsh_766 = buffer.data(nsh + 766);
    const auto *nsh_771 = buffer.data(nsh + 771);
    const auto *nsh_773 = buffer.data(nsh + 773);
    const auto *nsh_774 = buffer.data(nsh + 774);
    const auto *nsh_775 = buffer.data(nsh + 775);
    const auto *nsh_776 = buffer.data(nsh + 776);
    const auto *nsh_782 = buffer.data(nsh + 782);

    const auto *nsi1_756 = buffer.data(nsi1 + 756);
    const auto *nsi1_759 = buffer.data(nsi1 + 759);
    const auto *nsi1_761 = buffer.data(nsi1 + 761);
    const auto *nsi1_762 = buffer.data(nsi1 + 762);
    const auto *nsi1_765 = buffer.data(nsi1 + 765);
    const auto *nsi1_766 = buffer.data(nsi1 + 766);
    const auto *nsi1_768 = buffer.data(nsi1 + 768);
    const auto *nsi1_770 = buffer.data(nsi1 + 770);
    const auto *nsi1_783 = buffer.data(nsi1 + 783);
    const auto *nsi1_784 = buffer.data(nsi1 + 784);
    const auto *nsi1_787 = buffer.data(nsi1 + 787);
    const auto *nsi1_790 = buffer.data(nsi1 + 790);

    const auto *osg0_498 = buffer.data(osg0 + 498);
    const auto *osg0_500 = buffer.data(osg0 + 500);
    const auto *osg0_501 = buffer.data(osg0 + 501);
    const auto *osg0_504 = buffer.data(osg0 + 504);
    const auto *osg0_505 = buffer.data(osg0 + 505);
    const auto *osg0_507 = buffer.data(osg0 + 507);
    const auto *osg0_508 = buffer.data(osg0 + 508);
    const auto *osg0_509 = buffer.data(osg0 + 509);
    const auto *osg0_520 = buffer.data(osg0 + 520);
    const auto *osg0_522 = buffer.data(osg0 + 522);
    const auto *osg0_523 = buffer.data(osg0 + 523);
    const auto *osg0_524 = buffer.data(osg0 + 524);
    const auto *osg0_525 = buffer.data(osg0 + 525);
    const auto *osg0_526 = buffer.data(osg0 + 526);
    const auto *osg0_527 = buffer.data(osg0 + 527);
    const auto *osg0_528 = buffer.data(osg0 + 528);
    const auto *osg0_529 = buffer.data(osg0 + 529);
    const auto *osg0_530 = buffer.data(osg0 + 530);
    const auto *osg0_534 = buffer.data(osg0 + 534);
    const auto *osg0_535 = buffer.data(osg0 + 535);
    const auto *osg0_536 = buffer.data(osg0 + 536);
    const auto *osg0_537 = buffer.data(osg0 + 537);
    const auto *osg0_538 = buffer.data(osg0 + 538);
    const auto *osg0_539 = buffer.data(osg0 + 539);
    const auto *osg0_540 = buffer.data(osg0 + 540);
    const auto *osg0_542 = buffer.data(osg0 + 542);
    const auto *osg0_543 = buffer.data(osg0 + 543);
    const auto *osg0_545 = buffer.data(osg0 + 545);
    const auto *osg0_546 = buffer.data(osg0 + 546);
    const auto *osg0_550 = buffer.data(osg0 + 550);
    const auto *osg0_551 = buffer.data(osg0 + 551);
    const auto *osg0_552 = buffer.data(osg0 + 552);
    const auto *osg0_554 = buffer.data(osg0 + 554);
    const auto *osg0_560 = buffer.data(osg0 + 560);

    const auto *osg1_498 = buffer.data(osg1 + 498);
    const auto *osg1_500 = buffer.data(osg1 + 500);
    const auto *osg1_501 = buffer.data(osg1 + 501);
    const auto *osg1_504 = buffer.data(osg1 + 504);
    const auto *osg1_505 = buffer.data(osg1 + 505);
    const auto *osg1_507 = buffer.data(osg1 + 507);
    const auto *osg1_508 = buffer.data(osg1 + 508);
    const auto *osg1_509 = buffer.data(osg1 + 509);
    const auto *osg1_520 = buffer.data(osg1 + 520);
    const auto *osg1_522 = buffer.data(osg1 + 522);
    const auto *osg1_523 = buffer.data(osg1 + 523);
    const auto *osg1_524 = buffer.data(osg1 + 524);
    const auto *osg1_525 = buffer.data(osg1 + 525);
    const auto *osg1_526 = buffer.data(osg1 + 526);
    const auto *osg1_527 = buffer.data(osg1 + 527);
    const auto *osg1_528 = buffer.data(osg1 + 528);
    const auto *osg1_529 = buffer.data(osg1 + 529);
    const auto *osg1_530 = buffer.data(osg1 + 530);
    const auto *osg1_534 = buffer.data(osg1 + 534);
    const auto *osg1_535 = buffer.data(osg1 + 535);
    const auto *osg1_536 = buffer.data(osg1 + 536);
    const auto *osg1_537 = buffer.data(osg1 + 537);
    const auto *osg1_538 = buffer.data(osg1 + 538);
    const auto *osg1_539 = buffer.data(osg1 + 539);
    const auto *osg1_540 = buffer.data(osg1 + 540);
    const auto *osg1_542 = buffer.data(osg1 + 542);
    const auto *osg1_543 = buffer.data(osg1 + 543);
    const auto *osg1_545 = buffer.data(osg1 + 545);
    const auto *osg1_546 = buffer.data(osg1 + 546);
    const auto *osg1_550 = buffer.data(osg1 + 550);
    const auto *osg1_551 = buffer.data(osg1 + 551);
    const auto *osg1_552 = buffer.data(osg1 + 552);
    const auto *osg1_554 = buffer.data(osg1 + 554);
    const auto *osg1_560 = buffer.data(osg1 + 560);

    const auto *osh_693 = buffer.data(osh + 693);
    const auto *osh_695 = buffer.data(osh + 695);
    const auto *osh_696 = buffer.data(osh + 696);
    const auto *osh_698 = buffer.data(osh + 698);
    const auto *osh_699 = buffer.data(osh + 699);
    const auto *osh_702 = buffer.data(osh + 702);
    const auto *osh_703 = buffer.data(osh + 703);
    const auto *osh_705 = buffer.data(osh + 705);
    const auto *osh_707 = buffer.data(osh + 707);
    const auto *osh_708 = buffer.data(osh + 708);
    const auto *osh_709 = buffer.data(osh + 709);
    const auto *osh_710 = buffer.data(osh + 710);
    const auto *osh_711 = buffer.data(osh + 711);
    const auto *osh_712 = buffer.data(osh + 712);
    const auto *osh_713 = buffer.data(osh + 713);
    const auto *osh_714 = buffer.data(osh + 714);
    const auto *osh_716 = buffer.data(osh + 716);
    const auto *osh_717 = buffer.data(osh + 717);
    const auto *osh_719 = buffer.data(osh + 719);
    const auto *osh_720 = buffer.data(osh + 720);
    const auto *osh_723 = buffer.data(osh + 723);
    const auto *osh_729 = buffer.data(osh + 729);
    const auto *osh_730 = buffer.data(osh + 730);
    const auto *osh_731 = buffer.data(osh + 731);
    const auto *osh_732 = buffer.data(osh + 732);
    const auto *osh_733 = buffer.data(osh + 733);
    const auto *osh_734 = buffer.data(osh + 734);
    const auto *osh_735 = buffer.data(osh + 735);
    const auto *osh_736 = buffer.data(osh + 736);
    const auto *osh_737 = buffer.data(osh + 737);
    const auto *osh_738 = buffer.data(osh + 738);
    const auto *osh_739 = buffer.data(osh + 739);
    const auto *osh_740 = buffer.data(osh + 740);
    const auto *osh_741 = buffer.data(osh + 741);
    const auto *osh_742 = buffer.data(osh + 742);
    const auto *osh_743 = buffer.data(osh + 743);
    const auto *osh_744 = buffer.data(osh + 744);
    const auto *osh_749 = buffer.data(osh + 749);
    const auto *osh_750 = buffer.data(osh + 750);
    const auto *osh_751 = buffer.data(osh + 751);
    const auto *osh_752 = buffer.data(osh + 752);
    const auto *osh_753 = buffer.data(osh + 753);
    const auto *osh_754 = buffer.data(osh + 754);
    const auto *osh_755 = buffer.data(osh + 755);
    const auto *osh_756 = buffer.data(osh + 756);
    const auto *osh_757 = buffer.data(osh + 757);
    const auto *osh_758 = buffer.data(osh + 758);
    const auto *osh_759 = buffer.data(osh + 759);
    const auto *osh_761 = buffer.data(osh + 761);
    const auto *osh_762 = buffer.data(osh + 762);
    const auto *osh_763 = buffer.data(osh + 763);
    const auto *osh_765 = buffer.data(osh + 765);
    const auto *osh_766 = buffer.data(osh + 766);
    const auto *osh_771 = buffer.data(osh + 771);
    const auto *osh_772 = buffer.data(osh + 772);
    const auto *osh_773 = buffer.data(osh + 773);
    const auto *osh_774 = buffer.data(osh + 774);
    const auto *osh_775 = buffer.data(osh + 775);
    const auto *osh_776 = buffer.data(osh + 776);
    const auto *osh_777 = buffer.data(osh + 777);
    const auto *osh_779 = buffer.data(osh + 779);
    const auto *osh_780 = buffer.data(osh + 780);
    const auto *osh_782 = buffer.data(osh + 782);

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, nsh_525, nsh_546, \
                         nsh_548, nsh_696, osg0_498, osg1_498, osh_693, osh_695, \
                         osh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * nsh_546[k]
                   + f_3 * pc_y[k] * osh_693[k];

        t_926[k] = f_22 * nsh_525[k]
                   + f_3 * pc_z[k] * osh_693[k];

        t_927[k] = f_14 * nsh_696[k]
                   + f_8 * osg0_498[k]
                   - f_9 * osg1_498[k]
                   + f_3 * pc_x[k] * osh_696[k];

        t_928[k] = f_12 * nsh_548[k]
                   + f_3 * pc_y[k] * osh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, nsh_528, nsh_698, nsh_699, osg0_500, \
                         osg0_501, osg1_500, osg1_501, osh_696, osh_698, \
                         osh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_14 * nsh_698[k]
                   + f_8 * osg0_500[k]
                   - f_9 * osg1_500[k]
                   + f_3 * pc_x[k] * osh_698[k];

        t_930[k] = f_14 * nsh_699[k]
                   + f_6 * osg0_501[k]
                   - f_7 * osg1_501[k]
                   + f_3 * pc_x[k] * osh_699[k];

        t_931[k] = f_22 * nsh_528[k]
                   + f_3 * pc_z[k] * osh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, nsh_551, nsh_702, nsh_703, osg0_504, \
                         osg0_505, osg1_504, osg1_505, osh_698, osh_702, \
                         osh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * nsh_551[k]
                   + f_3 * pc_y[k] * osh_698[k];

        t_933[k] = f_14 * nsh_702[k]
                   + f_6 * osg0_504[k]
                   - f_7 * osg1_504[k]
                   + f_3 * pc_x[k] * osh_702[k];

        t_934[k] = f_14 * nsh_703[k]
                   + f_4 * osg0_505[k]
                   - f_5 * osg1_505[k]
                   + f_3 * pc_x[k] * osh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, nsh_531, nsh_555, nsh_705, \
                         osg0_507, osg1_507, osh_699, osh_702, \
                         osh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_22 * nsh_531[k]
                   + f_3 * pc_z[k] * osh_699[k];

        t_936[k] = f_14 * nsh_705[k]
                   + f_4 * osg0_507[k]
                   - f_5 * osg1_507[k]
                   + f_3 * pc_x[k] * osh_705[k];

        t_937[k] = f_12 * nsh_555[k]
                   + f_3 * pc_y[k] * osh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, nsh_707, nsh_708, nsh_709, nsh_710, \
                         osg0_509, osg1_509, osh_707, osh_708, osh_709, \
                         osh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_14 * nsh_707[k]
                   + f_4 * osg0_509[k]
                   - f_5 * osg1_509[k]
                   + f_3 * pc_x[k] * osh_707[k];

        t_939[k] = f_14 * nsh_708[k]
                   + f_3 * pc_x[k] * osh_708[k];

        t_940[k] = f_14 * nsh_709[k]
                   + f_3 * pc_x[k] * osh_709[k];

        t_941[k] = f_14 * nsh_710[k]
                   + f_3 * pc_x[k] * osh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, nsh_561, nsh_711, nsh_712, \
                         nsh_713, osg0_505, osg1_505, osh_708, osh_711, osh_712, \
                         osh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_14 * nsh_711[k]
                   + f_3 * pc_x[k] * osh_711[k];

        t_943[k] = f_14 * nsh_712[k]
                   + f_3 * pc_x[k] * osh_712[k];

        t_944[k] = f_14 * nsh_713[k]
                   + f_3 * pc_x[k] * osh_713[k];

        t_945[k] = f_12 * nsh_561[k]
                   + f_1 * osg0_505[k]
                   - f_2 * osg1_505[k]
                   + f_3 * pc_y[k] * osh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, nsh_540, nsh_563, nsh_564, osg0_507, \
                         osg0_508, osg1_507, osg1_508, osh_708, osh_710, \
                         osh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_22 * nsh_540[k]
                   + f_3 * pc_z[k] * osh_708[k];

        t_947[k] = f_12 * nsh_563[k]
                   + f_8 * osg0_507[k]
                   - f_9 * osg1_507[k]
                   + f_3 * pc_y[k] * osh_710[k];

        t_948[k] = f_12 * nsh_564[k]
                   + f_6 * osg0_508[k]
                   - f_7 * osg1_508[k]
                   + f_3 * pc_y[k] * osh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, nsi0_756, nsh_545, \
                         nsh_565, nsh_566, nsi1_756, osg0_509, osg1_509, osh_712, \
                         osh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * nsh_565[k]
                   + f_4 * osg0_509[k]
                   - f_5 * osg1_509[k]
                   + f_3 * pc_y[k] * osh_712[k];

        t_950[k] = f_12 * nsh_566[k]
                   + f_3 * pc_y[k] * osh_713[k];

        t_951[k] = f_22 * nsh_545[k]
                   + f_1 * osg0_509[k]
                   - f_2 * osg1_509[k]
                   + f_3 * pc_z[k] * osh_713[k];

        t_952[k] = pa_y[k] * nsi0_756[k]
                   - f_10 * pc_y[k] * nsi1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, pc_z, nsi0_759, nsh_546, \
                         nsh_567, nsh_568, nsh_569, nsi1_759, osh_714, \
                         osh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * nsh_567[k]
                   + f_3 * pc_y[k] * osh_714[k];

        t_954[k] = f_21 * nsh_546[k]
                   + f_3 * pc_z[k] * osh_714[k];

        t_955[k] = pa_y[k] * nsi0_759[k]
                   + f_12 * nsh_568[k]
                   - f_10 * pc_y[k] * nsi1_759[k];

        t_956[k] = f_11 * nsh_569[k]
                   + f_3 * pc_y[k] * osh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_y, pc_y, pc_z, nsi0_761, nsi0_762, \
                         nsh_549, nsh_570, nsh_572, nsi1_761, nsi1_762, osh_717, \
                         osh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pa_y[k] * nsi0_761[k]
                   - f_10 * pc_y[k] * nsi1_761[k];

        t_958[k] = pa_y[k] * nsi0_762[k]
                   + f_13 * nsh_570[k]
                   - f_10 * pc_y[k] * nsi1_762[k];

        t_959[k] = f_21 * nsh_549[k]
                   + f_3 * pc_z[k] * osh_717[k];

        t_960[k] = f_11 * nsh_572[k]
                   + f_3 * pc_y[k] * osh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pa_y, pc_y, pc_z, nsi0_765, nsi0_766, nsh_552, \
                         nsh_573, nsi1_765, nsi1_766, osh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pa_y[k] * nsi0_765[k]
                   - f_10 * pc_y[k] * nsi1_765[k];

        t_962[k] = pa_y[k] * nsi0_766[k]
                   + f_14 * nsh_573[k]
                   - f_10 * pc_y[k] * nsi1_766[k];

        t_963[k] = f_21 * nsh_552[k]
                   + f_3 * pc_z[k] * osh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pa_y, pc_x, pc_y, nsi0_768, nsi0_770, \
                         nsh_575, nsh_576, nsh_729, nsi1_768, nsi1_770, osh_723, \
                         osh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pa_y[k] * nsi0_768[k]
                   + f_12 * nsh_575[k]
                   - f_10 * pc_y[k] * nsi1_768[k];

        t_965[k] = f_11 * nsh_576[k]
                   + f_3 * pc_y[k] * osh_723[k];

        t_966[k] = pa_y[k] * nsi0_770[k]
                   - f_10 * pc_y[k] * nsi1_770[k];

        t_967[k] = f_14 * nsh_729[k]
                   + f_3 * pc_x[k] * osh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, nsh_730, nsh_731, nsh_732, \
                         nsh_733, nsh_734, osh_730, osh_731, osh_732, osh_733, \
                         osh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_14 * nsh_730[k]
                   + f_3 * pc_x[k] * osh_730[k];

        t_969[k] = f_14 * nsh_731[k]
                   + f_3 * pc_x[k] * osh_731[k];

        t_970[k] = f_14 * nsh_732[k]
                   + f_3 * pc_x[k] * osh_732[k];

        t_971[k] = f_14 * nsh_733[k]
                   + f_3 * pc_x[k] * osh_733[k];

        t_972[k] = f_14 * nsh_734[k]
                   + f_3 * pc_x[k] * osh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, nsh_561, nsh_582, nsh_584, osg0_520, \
                         osg0_522, osg1_520, osg1_522, osh_729, \
                         osh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * nsh_582[k]
                   + f_1 * osg0_520[k]
                   - f_2 * osg1_520[k]
                   + f_3 * pc_y[k] * osh_729[k];

        t_974[k] = f_21 * nsh_561[k]
                   + f_3 * pc_z[k] * osh_729[k];

        t_975[k] = f_11 * nsh_584[k]
                   + f_8 * osg0_522[k]
                   - f_9 * osg1_522[k]
                   + f_3 * pc_y[k] * osh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, nsh_585, nsh_586, nsh_587, osg0_523, \
                         osg0_524, osg1_523, osg1_524, osh_732, osh_733, \
                         osh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * nsh_585[k]
                   + f_6 * osg0_523[k]
                   - f_7 * osg1_523[k]
                   + f_3 * pc_y[k] * osh_732[k];

        t_977[k] = f_11 * nsh_586[k]
                   + f_4 * osg0_524[k]
                   - f_5 * osg1_524[k]
                   + f_3 * pc_y[k] * osh_733[k];

        t_978[k] = f_11 * nsh_587[k]
                   + f_3 * pc_y[k] * osh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_y, pc_x, pc_y, pc_z, nsi0_783, \
                         nsh_567, nsh_735, nsi1_783, osg0_525, osg1_525, \
                         osh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pa_y[k] * nsi0_783[k]
                   - f_10 * pc_y[k] * nsi1_783[k];

        t_980[k] = f_14 * nsh_735[k]
                   + f_1 * osg0_525[k]
                   - f_2 * osg1_525[k]
                   + f_3 * pc_x[k] * osh_735[k];

        t_981[k] = f_3 * pc_y[k] * osh_735[k];

        t_982[k] = f_20 * nsh_567[k]
                   + f_3 * pc_z[k] * osh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, nsh_740, osg0_525, osg0_530, \
                         osg1_525, osg1_530, osh_736, osh_737, \
                         osh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_4 * osg0_525[k]
                   - f_5 * osg1_525[k]
                   + f_3 * pc_y[k] * osh_736[k];

        t_984[k] = f_3 * pc_y[k] * osh_737[k];

        t_985[k] = f_14 * nsh_740[k]
                   + f_8 * osg0_530[k]
                   - f_9 * osg1_530[k]
                   + f_3 * pc_x[k] * osh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_y, osg0_526, osg0_527, osg1_526, osg1_527, \
                         osh_738, osh_739, osh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_6 * osg0_526[k]
                   - f_7 * osg1_526[k]
                   + f_3 * pc_y[k] * osh_738[k];

        t_987[k] = f_4 * osg0_527[k]
                   - f_5 * osg1_527[k]
                   + f_3 * pc_y[k] * osh_739[k];

        t_988[k] = f_3 * pc_y[k] * osh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, nsh_744, osg0_528, osg0_529, \
                         osg0_534, osg1_528, osg1_529, osg1_534, osh_741, osh_742, \
                         osh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_14 * nsh_744[k]
                   + f_6 * osg0_534[k]
                   - f_7 * osg1_534[k]
                   + f_3 * pc_x[k] * osh_744[k];

        t_990[k] = f_8 * osg0_528[k]
                   - f_9 * osg1_528[k]
                   + f_3 * pc_y[k] * osh_741[k];

        t_991[k] = f_6 * osg0_529[k]
                   - f_7 * osg1_529[k]
                   + f_3 * pc_y[k] * osh_742[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pc_x, pc_y, nsh_749, nsh_750, osg0_530, \
                         osg0_539, osg1_530, osg1_539, osh_743, osh_744, osh_749, \
                         osh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * osg0_530[k]
                   - f_5 * osg1_530[k]
                   + f_3 * pc_y[k] * osh_743[k];

        t_993[k] = f_3 * pc_y[k] * osh_744[k];

        t_994[k] = f_14 * nsh_749[k]
                   + f_4 * osg0_539[k]
                   - f_5 * osg1_539[k]
                   + f_3 * pc_x[k] * osh_749[k];

        t_995[k] = f_14 * nsh_750[k]
                   + f_3 * pc_x[k] * osh_750[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, pc_x, pc_y, nsh_751, nsh_752, \
                         nsh_753, nsh_755, osh_749, osh_751, osh_752, osh_753, \
                         osh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_14 * nsh_751[k]
                   + f_3 * pc_x[k] * osh_751[k];

        t_997[k] = f_14 * nsh_752[k]
                   + f_3 * pc_x[k] * osh_752[k];

        t_998[k] = f_14 * nsh_753[k]
                   + f_3 * pc_x[k] * osh_753[k];

        t_999[k] = f_3 * pc_y[k] * osh_749[k];

        t_1000[k] = f_14 * nsh_755[k]
                    + f_3 * pc_x[k] * osh_755[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pc_y, osg0_535, osg0_536, osg0_537, osg1_535, \
                         osg1_536, osg1_537, osh_750, osh_751, \
                         osh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_1 * osg0_535[k]
                    - f_2 * osg1_535[k]
                    + f_3 * pc_y[k] * osh_750[k];

        t_1002[k] = f_16 * osg0_536[k]
                    - f_17 * osg1_536[k]
                    + f_3 * pc_y[k] * osh_751[k];

        t_1003[k] = f_8 * osg0_537[k]
                    - f_9 * osg1_537[k]
                    + f_3 * pc_y[k] * osh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, nsh_587, osg0_538, \
                         osg0_539, osg1_538, osg1_539, osh_753, osh_754, \
                         osh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * osg0_538[k]
                    - f_7 * osg1_538[k]
                    + f_3 * pc_y[k] * osh_753[k];

        t_1005[k] = f_4 * osg0_539[k]
                    - f_5 * osg1_539[k]
                    + f_3 * pc_y[k] * osh_754[k];

        t_1006[k] = f_3 * pc_y[k] * osh_755[k];

        t_1007[k] = f_20 * nsh_587[k]
                    + f_1 * osg0_539[k]
                    - f_2 * osg1_539[k]
                    + f_3 * pc_z[k] * osh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, nsh_588, nsh_756, \
                         nsh_759, osg0_540, osg0_543, osg1_540, osg1_543, osh_756, \
                         osh_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_13 * nsh_756[k]
                    + f_1 * osg0_540[k]
                    - f_2 * osg1_540[k]
                    + f_3 * pc_x[k] * osh_756[k];

        t_1009[k] = f_19 * nsh_588[k]
                    + f_3 * pc_y[k] * osh_756[k];

        t_1010[k] = f_3 * pc_z[k] * osh_756[k];

        t_1011[k] = f_13 * nsh_759[k]
                    + f_8 * osg0_543[k]
                    - f_9 * osg1_543[k]
                    + f_3 * pc_x[k] * osh_759[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pc_x, pc_z, nsh_762, osg0_540, \
                         osg0_546, osg1_540, osg1_546, osh_757, osh_758, osh_759, \
                         osh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_3 * pc_z[k] * osh_757[k];

        t_1013[k] = f_4 * osg0_540[k]
                    - f_5 * osg1_540[k]
                    + f_3 * pc_z[k] * osh_758[k];

        t_1014[k] = f_13 * nsh_762[k]
                    + f_6 * osg0_546[k]
                    - f_7 * osg1_546[k]
                    + f_3 * pc_x[k] * osh_762[k];

        t_1015[k] = f_3 * pc_z[k] * osh_759[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, pc_z, nsh_593, nsh_766, \
                         osg0_542, osg0_550, osg1_542, osg1_550, osh_761, osh_762, \
                         osh_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * nsh_593[k]
                    + f_3 * pc_y[k] * osh_761[k];

        t_1017[k] = f_6 * osg0_542[k]
                    - f_7 * osg1_542[k]
                    + f_3 * pc_z[k] * osh_761[k];

        t_1018[k] = f_13 * nsh_766[k]
                    + f_4 * osg0_550[k]
                    - f_5 * osg1_550[k]
                    + f_3 * pc_x[k] * osh_766[k];

        t_1019[k] = f_3 * pc_z[k] * osh_762[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, pc_z, nsh_597, nsh_771, \
                         osg0_543, osg0_545, osg1_543, osg1_545, osh_763, osh_765, \
                         osh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * osg0_543[k]
                    - f_5 * osg1_543[k]
                    + f_3 * pc_z[k] * osh_763[k];

        t_1021[k] = f_19 * nsh_597[k]
                    + f_3 * pc_y[k] * osh_765[k];

        t_1022[k] = f_8 * osg0_545[k]
                    - f_9 * osg1_545[k]
                    + f_3 * pc_z[k] * osh_765[k];

        t_1023[k] = f_13 * nsh_771[k]
                    + f_3 * pc_x[k] * osh_771[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, t_1027, t_1028, pc_x, pc_z, nsh_773, nsh_774, \
                         nsh_775, nsh_776, osh_766, osh_773, osh_774, osh_775, \
                         osh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * osh_766[k];

        t_1025[k] = f_13 * nsh_773[k]
                    + f_3 * pc_x[k] * osh_773[k];

        t_1026[k] = f_13 * nsh_774[k]
                    + f_3 * pc_x[k] * osh_774[k];

        t_1027[k] = f_13 * nsh_775[k]
                    + f_3 * pc_x[k] * osh_775[k];

        t_1028[k] = f_13 * nsh_776[k]
                    + f_3 * pc_x[k] * osh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pc_y, pc_z, nsh_603, osg0_550, \
                         osg0_551, osg1_550, osg1_551, osh_771, osh_772, \
                         osh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_19 * nsh_603[k]
                    + f_1 * osg0_550[k]
                    - f_2 * osg1_550[k]
                    + f_3 * pc_y[k] * osh_771[k];

        t_1030[k] = f_3 * pc_z[k] * osh_771[k];

        t_1031[k] = f_4 * osg0_550[k]
                    - f_5 * osg1_550[k]
                    + f_3 * pc_z[k] * osh_772[k];

        t_1032[k] = f_6 * osg0_551[k]
                    - f_7 * osg1_551[k]
                    + f_3 * pc_z[k] * osh_773[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, pa_z, pc_y, pc_z, nsi0_784, nsh_608, \
                         nsi1_784, osg0_552, osg0_554, osg1_552, osg1_554, osh_774, \
                         osh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_8 * osg0_552[k]
                    - f_9 * osg1_552[k]
                    + f_3 * pc_z[k] * osh_774[k];

        t_1034[k] = f_19 * nsh_608[k]
                    + f_3 * pc_y[k] * osh_776[k];

        t_1035[k] = f_1 * osg0_554[k]
                    - f_2 * osg1_554[k]
                    + f_3 * pc_z[k] * osh_776[k];

        t_1036[k] = pa_z[k] * nsi0_784[k]
                    - f_10 * pc_z[k] * nsi1_784[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, pa_z, pc_y, pc_z, nsi0_787, nsh_588, \
                         nsh_609, nsh_611, nsi1_787, osh_777, osh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_20 * nsh_609[k]
                    + f_3 * pc_y[k] * osh_777[k];

        t_1038[k] = f_11 * nsh_588[k]
                    + f_3 * pc_z[k] * osh_777[k];

        t_1039[k] = pa_z[k] * nsi0_787[k]
                    - f_10 * pc_z[k] * nsi1_787[k];

        t_1040[k] = f_20 * nsh_611[k]
                    + f_3 * pc_y[k] * osh_779[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, pa_z, pc_x, pc_z, nsi0_790, nsh_591, nsh_782, \
                         nsi1_790, osg0_560, osg1_560, osh_780, \
                         osh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_13 * nsh_782[k]
                    + f_8 * osg0_560[k]
                    - f_9 * osg1_560[k]
                    + f_3 * pc_x[k] * osh_782[k];

        t_1042[k] = pa_z[k] * nsi0_790[k]
                    - f_10 * pc_z[k] * nsi1_790[k];

        t_1043[k] = f_11 * nsh_591[k]
                    + f_3 * pc_z[k] * osh_780[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsi0,
                                                          const size_t nsh, const size_t nsi1,
                                                          const size_t osg0, const size_t osg1,
                                                          const size_t osh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_794 = buffer.data(nsi0 + 794);
    const auto *nsi0_796 = buffer.data(nsi0 + 796);
    const auto *nsi0_805 = buffer.data(nsi0 + 805);

    const auto *nsh_594 = buffer.data(nsh + 594);
    const auto *nsh_595 = buffer.data(nsh + 595);
    const auto *nsh_603 = buffer.data(nsh + 603);
    const auto *nsh_608 = buffer.data(nsh + 608);
    const auto *nsh_609 = buffer.data(nsh + 609);
    const auto *nsh_612 = buffer.data(nsh + 612);
    const auto *nsh_614 = buffer.data(nsh + 614);
    const auto *nsh_615 = buffer.data(nsh + 615);
    const auto *nsh_618 = buffer.data(nsh + 618);
    const auto *nsh_624 = buffer.data(nsh + 624);
    const auto *nsh_626 = buffer.data(nsh + 626);
    const auto *nsh_627 = buffer.data(nsh + 627);
    const auto *nsh_628 = buffer.data(nsh + 628);
    const auto *nsh_629 = buffer.data(nsh + 629);
    const auto *nsh_630 = buffer.data(nsh + 630);
    const auto *nsh_632 = buffer.data(nsh + 632);
    const auto *nsh_633 = buffer.data(nsh + 633);
    const auto *nsh_635 = buffer.data(nsh + 635);
    const auto *nsh_636 = buffer.data(nsh + 636);
    const auto *nsh_639 = buffer.data(nsh + 639);
    const auto *nsh_645 = buffer.data(nsh + 645);
    const auto *nsh_647 = buffer.data(nsh + 647);
    const auto *nsh_648 = buffer.data(nsh + 648);
    const auto *nsh_649 = buffer.data(nsh + 649);
    const auto *nsh_650 = buffer.data(nsh + 650);
    const auto *nsh_651 = buffer.data(nsh + 651);
    const auto *nsh_653 = buffer.data(nsh + 653);
    const auto *nsh_654 = buffer.data(nsh + 654);
    const auto *nsh_656 = buffer.data(nsh + 656);
    const auto *nsh_657 = buffer.data(nsh + 657);
    const auto *nsh_660 = buffer.data(nsh + 660);
    const auto *nsh_666 = buffer.data(nsh + 666);
    const auto *nsh_668 = buffer.data(nsh + 668);
    const auto *nsh_669 = buffer.data(nsh + 669);
    const auto *nsh_670 = buffer.data(nsh + 670);
    const auto *nsh_671 = buffer.data(nsh + 671);
    const auto *nsh_672 = buffer.data(nsh + 672);
    const auto *nsh_674 = buffer.data(nsh + 674);
    const auto *nsh_677 = buffer.data(nsh + 677);
    const auto *nsh_681 = buffer.data(nsh + 681);
    const auto *nsh_687 = buffer.data(nsh + 687);
    const auto *nsh_689 = buffer.data(nsh + 689);
    const auto *nsh_690 = buffer.data(nsh + 690);
    const auto *nsh_691 = buffer.data(nsh + 691);
    const auto *nsh_692 = buffer.data(nsh + 692);
    const auto *nsh_786 = buffer.data(nsh + 786);
    const auto *nsh_791 = buffer.data(nsh + 791);
    const auto *nsh_792 = buffer.data(nsh + 792);
    const auto *nsh_793 = buffer.data(nsh + 793);
    const auto *nsh_794 = buffer.data(nsh + 794);
    const auto *nsh_795 = buffer.data(nsh + 795);
    const auto *nsh_796 = buffer.data(nsh + 796);
    const auto *nsh_797 = buffer.data(nsh + 797);
    const auto *nsh_798 = buffer.data(nsh + 798);
    const auto *nsh_801 = buffer.data(nsh + 801);
    const auto *nsh_803 = buffer.data(nsh + 803);
    const auto *nsh_804 = buffer.data(nsh + 804);
    const auto *nsh_807 = buffer.data(nsh + 807);
    const auto *nsh_808 = buffer.data(nsh + 808);
    const auto *nsh_810 = buffer.data(nsh + 810);
    const auto *nsh_812 = buffer.data(nsh + 812);
    const auto *nsh_813 = buffer.data(nsh + 813);
    const auto *nsh_814 = buffer.data(nsh + 814);
    const auto *nsh_815 = buffer.data(nsh + 815);
    const auto *nsh_816 = buffer.data(nsh + 816);
    const auto *nsh_817 = buffer.data(nsh + 817);
    const auto *nsh_818 = buffer.data(nsh + 818);
    const auto *nsh_819 = buffer.data(nsh + 819);
    const auto *nsh_822 = buffer.data(nsh + 822);
    const auto *nsh_824 = buffer.data(nsh + 824);
    const auto *nsh_825 = buffer.data(nsh + 825);
    const auto *nsh_828 = buffer.data(nsh + 828);
    const auto *nsh_829 = buffer.data(nsh + 829);
    const auto *nsh_831 = buffer.data(nsh + 831);
    const auto *nsh_833 = buffer.data(nsh + 833);
    const auto *nsh_834 = buffer.data(nsh + 834);
    const auto *nsh_835 = buffer.data(nsh + 835);
    const auto *nsh_836 = buffer.data(nsh + 836);
    const auto *nsh_837 = buffer.data(nsh + 837);
    const auto *nsh_838 = buffer.data(nsh + 838);
    const auto *nsh_839 = buffer.data(nsh + 839);
    const auto *nsh_840 = buffer.data(nsh + 840);
    const auto *nsh_843 = buffer.data(nsh + 843);
    const auto *nsh_845 = buffer.data(nsh + 845);
    const auto *nsh_846 = buffer.data(nsh + 846);
    const auto *nsh_849 = buffer.data(nsh + 849);
    const auto *nsh_850 = buffer.data(nsh + 850);
    const auto *nsh_852 = buffer.data(nsh + 852);
    const auto *nsh_854 = buffer.data(nsh + 854);
    const auto *nsh_855 = buffer.data(nsh + 855);
    const auto *nsh_856 = buffer.data(nsh + 856);
    const auto *nsh_857 = buffer.data(nsh + 857);
    const auto *nsh_858 = buffer.data(nsh + 858);
    const auto *nsh_859 = buffer.data(nsh + 859);
    const auto *nsh_860 = buffer.data(nsh + 860);
    const auto *nsh_861 = buffer.data(nsh + 861);

    const auto *nsi1_794 = buffer.data(nsi1 + 794);
    const auto *nsi1_796 = buffer.data(nsi1 + 796);
    const auto *nsi1_805 = buffer.data(nsi1 + 805);

    const auto *osg0_564 = buffer.data(osg0 + 564);
    const auto *osg0_567 = buffer.data(osg0 + 567);
    const auto *osg0_568 = buffer.data(osg0 + 568);
    const auto *osg0_569 = buffer.data(osg0 + 569);
    const auto *osg0_570 = buffer.data(osg0 + 570);
    const auto *osg0_573 = buffer.data(osg0 + 573);
    const auto *osg0_575 = buffer.data(osg0 + 575);
    const auto *osg0_576 = buffer.data(osg0 + 576);
    const auto *osg0_579 = buffer.data(osg0 + 579);
    const auto *osg0_580 = buffer.data(osg0 + 580);
    const auto *osg0_582 = buffer.data(osg0 + 582);
    const auto *osg0_583 = buffer.data(osg0 + 583);
    const auto *osg0_584 = buffer.data(osg0 + 584);
    const auto *osg0_585 = buffer.data(osg0 + 585);
    const auto *osg0_588 = buffer.data(osg0 + 588);
    const auto *osg0_590 = buffer.data(osg0 + 590);
    const auto *osg0_591 = buffer.data(osg0 + 591);
    const auto *osg0_594 = buffer.data(osg0 + 594);
    const auto *osg0_595 = buffer.data(osg0 + 595);
    const auto *osg0_597 = buffer.data(osg0 + 597);
    const auto *osg0_598 = buffer.data(osg0 + 598);
    const auto *osg0_599 = buffer.data(osg0 + 599);
    const auto *osg0_600 = buffer.data(osg0 + 600);
    const auto *osg0_603 = buffer.data(osg0 + 603);
    const auto *osg0_605 = buffer.data(osg0 + 605);
    const auto *osg0_606 = buffer.data(osg0 + 606);
    const auto *osg0_609 = buffer.data(osg0 + 609);
    const auto *osg0_610 = buffer.data(osg0 + 610);
    const auto *osg0_612 = buffer.data(osg0 + 612);
    const auto *osg0_613 = buffer.data(osg0 + 613);
    const auto *osg0_614 = buffer.data(osg0 + 614);
    const auto *osg0_615 = buffer.data(osg0 + 615);

    const auto *osg1_564 = buffer.data(osg1 + 564);
    const auto *osg1_567 = buffer.data(osg1 + 567);
    const auto *osg1_568 = buffer.data(osg1 + 568);
    const auto *osg1_569 = buffer.data(osg1 + 569);
    const auto *osg1_570 = buffer.data(osg1 + 570);
    const auto *osg1_573 = buffer.data(osg1 + 573);
    const auto *osg1_575 = buffer.data(osg1 + 575);
    const auto *osg1_576 = buffer.data(osg1 + 576);
    const auto *osg1_579 = buffer.data(osg1 + 579);
    const auto *osg1_580 = buffer.data(osg1 + 580);
    const auto *osg1_582 = buffer.data(osg1 + 582);
    const auto *osg1_583 = buffer.data(osg1 + 583);
    const auto *osg1_584 = buffer.data(osg1 + 584);
    const auto *osg1_585 = buffer.data(osg1 + 585);
    const auto *osg1_588 = buffer.data(osg1 + 588);
    const auto *osg1_590 = buffer.data(osg1 + 590);
    const auto *osg1_591 = buffer.data(osg1 + 591);
    const auto *osg1_594 = buffer.data(osg1 + 594);
    const auto *osg1_595 = buffer.data(osg1 + 595);
    const auto *osg1_597 = buffer.data(osg1 + 597);
    const auto *osg1_598 = buffer.data(osg1 + 598);
    const auto *osg1_599 = buffer.data(osg1 + 599);
    const auto *osg1_600 = buffer.data(osg1 + 600);
    const auto *osg1_603 = buffer.data(osg1 + 603);
    const auto *osg1_605 = buffer.data(osg1 + 605);
    const auto *osg1_606 = buffer.data(osg1 + 606);
    const auto *osg1_609 = buffer.data(osg1 + 609);
    const auto *osg1_610 = buffer.data(osg1 + 610);
    const auto *osg1_612 = buffer.data(osg1 + 612);
    const auto *osg1_613 = buffer.data(osg1 + 613);
    const auto *osg1_614 = buffer.data(osg1 + 614);
    const auto *osg1_615 = buffer.data(osg1 + 615);

    const auto *osh_782 = buffer.data(osh + 782);
    const auto *osh_783 = buffer.data(osh + 783);
    const auto *osh_786 = buffer.data(osh + 786);
    const auto *osh_791 = buffer.data(osh + 791);
    const auto *osh_792 = buffer.data(osh + 792);
    const auto *osh_793 = buffer.data(osh + 793);
    const auto *osh_794 = buffer.data(osh + 794);
    const auto *osh_795 = buffer.data(osh + 795);
    const auto *osh_796 = buffer.data(osh + 796);
    const auto *osh_797 = buffer.data(osh + 797);
    const auto *osh_798 = buffer.data(osh + 798);
    const auto *osh_800 = buffer.data(osh + 800);
    const auto *osh_801 = buffer.data(osh + 801);
    const auto *osh_803 = buffer.data(osh + 803);
    const auto *osh_804 = buffer.data(osh + 804);
    const auto *osh_807 = buffer.data(osh + 807);
    const auto *osh_808 = buffer.data(osh + 808);
    const auto *osh_810 = buffer.data(osh + 810);
    const auto *osh_812 = buffer.data(osh + 812);
    const auto *osh_813 = buffer.data(osh + 813);
    const auto *osh_814 = buffer.data(osh + 814);
    const auto *osh_815 = buffer.data(osh + 815);
    const auto *osh_816 = buffer.data(osh + 816);
    const auto *osh_817 = buffer.data(osh + 817);
    const auto *osh_818 = buffer.data(osh + 818);
    const auto *osh_819 = buffer.data(osh + 819);
    const auto *osh_821 = buffer.data(osh + 821);
    const auto *osh_822 = buffer.data(osh + 822);
    const auto *osh_824 = buffer.data(osh + 824);
    const auto *osh_825 = buffer.data(osh + 825);
    const auto *osh_828 = buffer.data(osh + 828);
    const auto *osh_829 = buffer.data(osh + 829);
    const auto *osh_831 = buffer.data(osh + 831);
    const auto *osh_833 = buffer.data(osh + 833);
    const auto *osh_834 = buffer.data(osh + 834);
    const auto *osh_835 = buffer.data(osh + 835);
    const auto *osh_836 = buffer.data(osh + 836);
    const auto *osh_837 = buffer.data(osh + 837);
    const auto *osh_838 = buffer.data(osh + 838);
    const auto *osh_839 = buffer.data(osh + 839);
    const auto *osh_840 = buffer.data(osh + 840);
    const auto *osh_842 = buffer.data(osh + 842);
    const auto *osh_843 = buffer.data(osh + 843);
    const auto *osh_845 = buffer.data(osh + 845);
    const auto *osh_846 = buffer.data(osh + 846);
    const auto *osh_849 = buffer.data(osh + 849);
    const auto *osh_850 = buffer.data(osh + 850);
    const auto *osh_852 = buffer.data(osh + 852);
    const auto *osh_854 = buffer.data(osh + 854);
    const auto *osh_855 = buffer.data(osh + 855);
    const auto *osh_856 = buffer.data(osh + 856);
    const auto *osh_857 = buffer.data(osh + 857);
    const auto *osh_858 = buffer.data(osh + 858);
    const auto *osh_859 = buffer.data(osh + 859);
    const auto *osh_860 = buffer.data(osh + 860);
    const auto *osh_861 = buffer.data(osh + 861);

#pragma omp simd aligned(t_1044, t_1045, t_1046, pa_z, pc_x, pc_y, pc_z, nsi0_794, nsh_614, \
                         nsh_786, nsi1_794, osg0_564, osg1_564, osh_782, \
                         osh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_20 * nsh_614[k]
                    + f_3 * pc_y[k] * osh_782[k];

        t_1045[k] = f_13 * nsh_786[k]
                    + f_6 * osg0_564[k]
                    - f_7 * osg1_564[k]
                    + f_3 * pc_x[k] * osh_786[k];

        t_1046[k] = pa_z[k] * nsi0_794[k]
                    - f_10 * pc_z[k] * nsi1_794[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pa_z, pc_y, pc_z, nsi0_796, nsh_594, nsh_595, \
                         nsh_618, nsi1_796, osh_783, osh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_11 * nsh_594[k]
                    + f_3 * pc_z[k] * osh_783[k];

        t_1048[k] = pa_z[k] * nsi0_796[k]
                    + f_12 * nsh_595[k]
                    - f_10 * pc_z[k] * nsi1_796[k];

        t_1049[k] = f_20 * nsh_618[k]
                    + f_3 * pc_y[k] * osh_786[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, nsh_791, nsh_792, nsh_793, \
                         nsh_794, osg0_569, osg1_569, osh_791, osh_792, osh_793, \
                         osh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_13 * nsh_791[k]
                    + f_4 * osg0_569[k]
                    - f_5 * osg1_569[k]
                    + f_3 * pc_x[k] * osh_791[k];

        t_1051[k] = f_13 * nsh_792[k]
                    + f_3 * pc_x[k] * osh_792[k];

        t_1052[k] = f_13 * nsh_793[k]
                    + f_3 * pc_x[k] * osh_793[k];

        t_1053[k] = f_13 * nsh_794[k]
                    + f_3 * pc_x[k] * osh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, pa_z, pc_x, pc_z, nsi0_805, nsh_795, \
                         nsh_796, nsh_797, nsi1_805, osh_795, osh_796, \
                         osh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_13 * nsh_795[k]
                    + f_3 * pc_x[k] * osh_795[k];

        t_1055[k] = f_13 * nsh_796[k]
                    + f_3 * pc_x[k] * osh_796[k];

        t_1056[k] = f_13 * nsh_797[k]
                    + f_3 * pc_x[k] * osh_797[k];

        t_1057[k] = pa_z[k] * nsi0_805[k]
                    - f_10 * pc_z[k] * nsi1_805[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_y, pc_z, nsh_603, nsh_626, nsh_627, \
                         osg0_567, osg0_568, osg1_567, osg1_568, osh_792, osh_794, \
                         osh_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_11 * nsh_603[k]
                    + f_3 * pc_z[k] * osh_792[k];

        t_1059[k] = f_20 * nsh_626[k]
                    + f_8 * osg0_567[k]
                    - f_9 * osg1_567[k]
                    + f_3 * pc_y[k] * osh_794[k];

        t_1060[k] = f_20 * nsh_627[k]
                    + f_6 * osg0_568[k]
                    - f_7 * osg1_568[k]
                    + f_3 * pc_y[k] * osh_795[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_y, pc_z, nsh_608, nsh_628, nsh_629, \
                         osg0_569, osg1_569, osh_796, osh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_20 * nsh_628[k]
                    + f_4 * osg0_569[k]
                    - f_5 * osg1_569[k]
                    + f_3 * pc_y[k] * osh_796[k];

        t_1062[k] = f_20 * nsh_629[k]
                    + f_3 * pc_y[k] * osh_797[k];

        t_1063[k] = f_11 * nsh_608[k]
                    + f_1 * osg0_569[k]
                    - f_2 * osg1_569[k]
                    + f_3 * pc_z[k] * osh_797[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pc_x, pc_y, pc_z, nsh_609, nsh_630, nsh_798, \
                         osg0_570, osg1_570, osh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_13 * nsh_798[k]
                    + f_1 * osg0_570[k]
                    - f_2 * osg1_570[k]
                    + f_3 * pc_x[k] * osh_798[k];

        t_1065[k] = f_21 * nsh_630[k]
                    + f_3 * pc_y[k] * osh_798[k];

        t_1066[k] = f_12 * nsh_609[k]
                    + f_3 * pc_z[k] * osh_798[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_x, pc_y, nsh_632, nsh_801, nsh_803, \
                         osg0_573, osg0_575, osg1_573, osg1_575, osh_800, osh_801, \
                         osh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_13 * nsh_801[k]
                    + f_8 * osg0_573[k]
                    - f_9 * osg1_573[k]
                    + f_3 * pc_x[k] * osh_801[k];

        t_1068[k] = f_21 * nsh_632[k]
                    + f_3 * pc_y[k] * osh_800[k];

        t_1069[k] = f_13 * nsh_803[k]
                    + f_8 * osg0_575[k]
                    - f_9 * osg1_575[k]
                    + f_3 * pc_x[k] * osh_803[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, nsh_612, nsh_635, nsh_804, \
                         osg0_576, osg1_576, osh_801, osh_803, \
                         osh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_13 * nsh_804[k]
                    + f_6 * osg0_576[k]
                    - f_7 * osg1_576[k]
                    + f_3 * pc_x[k] * osh_804[k];

        t_1071[k] = f_12 * nsh_612[k]
                    + f_3 * pc_z[k] * osh_801[k];

        t_1072[k] = f_21 * nsh_635[k]
                    + f_3 * pc_y[k] * osh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_z, nsh_615, nsh_807, nsh_808, \
                         osg0_579, osg0_580, osg1_579, osg1_580, osh_804, osh_807, \
                         osh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_13 * nsh_807[k]
                    + f_6 * osg0_579[k]
                    - f_7 * osg1_579[k]
                    + f_3 * pc_x[k] * osh_807[k];

        t_1074[k] = f_13 * nsh_808[k]
                    + f_4 * osg0_580[k]
                    - f_5 * osg1_580[k]
                    + f_3 * pc_x[k] * osh_808[k];

        t_1075[k] = f_12 * nsh_615[k]
                    + f_3 * pc_z[k] * osh_804[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_y, nsh_639, nsh_810, nsh_812, \
                         osg0_582, osg0_584, osg1_582, osg1_584, osh_807, osh_810, \
                         osh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_13 * nsh_810[k]
                    + f_4 * osg0_582[k]
                    - f_5 * osg1_582[k]
                    + f_3 * pc_x[k] * osh_810[k];

        t_1077[k] = f_21 * nsh_639[k]
                    + f_3 * pc_y[k] * osh_807[k];

        t_1078[k] = f_13 * nsh_812[k]
                    + f_4 * osg0_584[k]
                    - f_5 * osg1_584[k]
                    + f_3 * pc_x[k] * osh_812[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, nsh_813, nsh_814, \
                         nsh_815, nsh_816, nsh_817, osh_813, osh_814, osh_815, osh_816, \
                         osh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_13 * nsh_813[k]
                    + f_3 * pc_x[k] * osh_813[k];

        t_1080[k] = f_13 * nsh_814[k]
                    + f_3 * pc_x[k] * osh_814[k];

        t_1081[k] = f_13 * nsh_815[k]
                    + f_3 * pc_x[k] * osh_815[k];

        t_1082[k] = f_13 * nsh_816[k]
                    + f_3 * pc_x[k] * osh_816[k];

        t_1083[k] = f_13 * nsh_817[k]
                    + f_3 * pc_x[k] * osh_817[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, pc_x, pc_y, pc_z, nsh_624, nsh_645, nsh_818, \
                         osg0_580, osg1_580, osh_813, osh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_13 * nsh_818[k]
                    + f_3 * pc_x[k] * osh_818[k];

        t_1085[k] = f_21 * nsh_645[k]
                    + f_1 * osg0_580[k]
                    - f_2 * osg1_580[k]
                    + f_3 * pc_y[k] * osh_813[k];

        t_1086[k] = f_12 * nsh_624[k]
                    + f_3 * pc_z[k] * osh_813[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, nsh_647, nsh_648, nsh_649, osg0_582, \
                         osg0_583, osg0_584, osg1_582, osg1_583, osg1_584, osh_815, osh_816, \
                         osh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_21 * nsh_647[k]
                    + f_8 * osg0_582[k]
                    - f_9 * osg1_582[k]
                    + f_3 * pc_y[k] * osh_815[k];

        t_1088[k] = f_21 * nsh_648[k]
                    + f_6 * osg0_583[k]
                    - f_7 * osg1_583[k]
                    + f_3 * pc_y[k] * osh_816[k];

        t_1089[k] = f_21 * nsh_649[k]
                    + f_4 * osg0_584[k]
                    - f_5 * osg1_584[k]
                    + f_3 * pc_y[k] * osh_817[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, nsh_629, nsh_650, nsh_819, \
                         osg0_584, osg0_585, osg1_584, osg1_585, osh_818, \
                         osh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_21 * nsh_650[k]
                    + f_3 * pc_y[k] * osh_818[k];

        t_1091[k] = f_12 * nsh_629[k]
                    + f_1 * osg0_584[k]
                    - f_2 * osg1_584[k]
                    + f_3 * pc_z[k] * osh_818[k];

        t_1092[k] = f_13 * nsh_819[k]
                    + f_1 * osg0_585[k]
                    - f_2 * osg1_585[k]
                    + f_3 * pc_x[k] * osh_819[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, nsh_630, nsh_651, \
                         nsh_653, nsh_822, osg0_588, osg1_588, osh_819, osh_821, \
                         osh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_22 * nsh_651[k]
                    + f_3 * pc_y[k] * osh_819[k];

        t_1094[k] = f_13 * nsh_630[k]
                    + f_3 * pc_z[k] * osh_819[k];

        t_1095[k] = f_13 * nsh_822[k]
                    + f_8 * osg0_588[k]
                    - f_9 * osg1_588[k]
                    + f_3 * pc_x[k] * osh_822[k];

        t_1096[k] = f_22 * nsh_653[k]
                    + f_3 * pc_y[k] * osh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, nsh_633, nsh_824, nsh_825, \
                         osg0_590, osg0_591, osg1_590, osg1_591, osh_822, osh_824, \
                         osh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_13 * nsh_824[k]
                    + f_8 * osg0_590[k]
                    - f_9 * osg1_590[k]
                    + f_3 * pc_x[k] * osh_824[k];

        t_1098[k] = f_13 * nsh_825[k]
                    + f_6 * osg0_591[k]
                    - f_7 * osg1_591[k]
                    + f_3 * pc_x[k] * osh_825[k];

        t_1099[k] = f_13 * nsh_633[k]
                    + f_3 * pc_z[k] * osh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_y, nsh_656, nsh_828, nsh_829, \
                         osg0_594, osg0_595, osg1_594, osg1_595, osh_824, osh_828, \
                         osh_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_22 * nsh_656[k]
                    + f_3 * pc_y[k] * osh_824[k];

        t_1101[k] = f_13 * nsh_828[k]
                    + f_6 * osg0_594[k]
                    - f_7 * osg1_594[k]
                    + f_3 * pc_x[k] * osh_828[k];

        t_1102[k] = f_13 * nsh_829[k]
                    + f_4 * osg0_595[k]
                    - f_5 * osg1_595[k]
                    + f_3 * pc_x[k] * osh_829[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, pc_y, pc_z, nsh_636, nsh_660, nsh_831, \
                         osg0_597, osg1_597, osh_825, osh_828, \
                         osh_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * nsh_636[k]
                    + f_3 * pc_z[k] * osh_825[k];

        t_1104[k] = f_13 * nsh_831[k]
                    + f_4 * osg0_597[k]
                    - f_5 * osg1_597[k]
                    + f_3 * pc_x[k] * osh_831[k];

        t_1105[k] = f_22 * nsh_660[k]
                    + f_3 * pc_y[k] * osh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, nsh_833, nsh_834, nsh_835, \
                         nsh_836, osg0_599, osg1_599, osh_833, osh_834, osh_835, \
                         osh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_13 * nsh_833[k]
                    + f_4 * osg0_599[k]
                    - f_5 * osg1_599[k]
                    + f_3 * pc_x[k] * osh_833[k];

        t_1107[k] = f_13 * nsh_834[k]
                    + f_3 * pc_x[k] * osh_834[k];

        t_1108[k] = f_13 * nsh_835[k]
                    + f_3 * pc_x[k] * osh_835[k];

        t_1109[k] = f_13 * nsh_836[k]
                    + f_3 * pc_x[k] * osh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pc_x, pc_y, nsh_666, nsh_837, \
                         nsh_838, nsh_839, osg0_595, osg1_595, osh_834, osh_837, osh_838, \
                         osh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_13 * nsh_837[k]
                    + f_3 * pc_x[k] * osh_837[k];

        t_1111[k] = f_13 * nsh_838[k]
                    + f_3 * pc_x[k] * osh_838[k];

        t_1112[k] = f_13 * nsh_839[k]
                    + f_3 * pc_x[k] * osh_839[k];

        t_1113[k] = f_22 * nsh_666[k]
                    + f_1 * osg0_595[k]
                    - f_2 * osg1_595[k]
                    + f_3 * pc_y[k] * osh_834[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_y, pc_z, nsh_645, nsh_668, nsh_669, \
                         osg0_597, osg0_598, osg1_597, osg1_598, osh_834, osh_836, \
                         osh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * nsh_645[k]
                    + f_3 * pc_z[k] * osh_834[k];

        t_1115[k] = f_22 * nsh_668[k]
                    + f_8 * osg0_597[k]
                    - f_9 * osg1_597[k]
                    + f_3 * pc_y[k] * osh_836[k];

        t_1116[k] = f_22 * nsh_669[k]
                    + f_6 * osg0_598[k]
                    - f_7 * osg1_598[k]
                    + f_3 * pc_y[k] * osh_837[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, pc_y, pc_z, nsh_650, nsh_670, nsh_671, \
                         osg0_599, osg1_599, osh_838, osh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_22 * nsh_670[k]
                    + f_4 * osg0_599[k]
                    - f_5 * osg1_599[k]
                    + f_3 * pc_y[k] * osh_838[k];

        t_1118[k] = f_22 * nsh_671[k]
                    + f_3 * pc_y[k] * osh_839[k];

        t_1119[k] = f_13 * nsh_650[k]
                    + f_1 * osg0_599[k]
                    - f_2 * osg1_599[k]
                    + f_3 * pc_z[k] * osh_839[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, pc_y, pc_z, nsh_651, nsh_672, nsh_840, \
                         osg0_600, osg1_600, osh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_13 * nsh_840[k]
                    + f_1 * osg0_600[k]
                    - f_2 * osg1_600[k]
                    + f_3 * pc_x[k] * osh_840[k];

        t_1121[k] = f_14 * nsh_672[k]
                    + f_3 * pc_y[k] * osh_840[k];

        t_1122[k] = f_14 * nsh_651[k]
                    + f_3 * pc_z[k] * osh_840[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, pc_y, nsh_674, nsh_843, nsh_845, \
                         osg0_603, osg0_605, osg1_603, osg1_605, osh_842, osh_843, \
                         osh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_13 * nsh_843[k]
                    + f_8 * osg0_603[k]
                    - f_9 * osg1_603[k]
                    + f_3 * pc_x[k] * osh_843[k];

        t_1124[k] = f_14 * nsh_674[k]
                    + f_3 * pc_y[k] * osh_842[k];

        t_1125[k] = f_13 * nsh_845[k]
                    + f_8 * osg0_605[k]
                    - f_9 * osg1_605[k]
                    + f_3 * pc_x[k] * osh_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, nsh_654, nsh_677, nsh_846, \
                         osg0_606, osg1_606, osh_843, osh_845, \
                         osh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_13 * nsh_846[k]
                    + f_6 * osg0_606[k]
                    - f_7 * osg1_606[k]
                    + f_3 * pc_x[k] * osh_846[k];

        t_1127[k] = f_14 * nsh_654[k]
                    + f_3 * pc_z[k] * osh_843[k];

        t_1128[k] = f_14 * nsh_677[k]
                    + f_3 * pc_y[k] * osh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, pc_z, nsh_657, nsh_849, nsh_850, \
                         osg0_609, osg0_610, osg1_609, osg1_610, osh_846, osh_849, \
                         osh_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_13 * nsh_849[k]
                    + f_6 * osg0_609[k]
                    - f_7 * osg1_609[k]
                    + f_3 * pc_x[k] * osh_849[k];

        t_1130[k] = f_13 * nsh_850[k]
                    + f_4 * osg0_610[k]
                    - f_5 * osg1_610[k]
                    + f_3 * pc_x[k] * osh_850[k];

        t_1131[k] = f_14 * nsh_657[k]
                    + f_3 * pc_z[k] * osh_846[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, pc_x, pc_y, nsh_681, nsh_852, nsh_854, \
                         osg0_612, osg0_614, osg1_612, osg1_614, osh_849, osh_852, \
                         osh_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_13 * nsh_852[k]
                    + f_4 * osg0_612[k]
                    - f_5 * osg1_612[k]
                    + f_3 * pc_x[k] * osh_852[k];

        t_1133[k] = f_14 * nsh_681[k]
                    + f_3 * pc_y[k] * osh_849[k];

        t_1134[k] = f_13 * nsh_854[k]
                    + f_4 * osg0_614[k]
                    - f_5 * osg1_614[k]
                    + f_3 * pc_x[k] * osh_854[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, nsh_855, nsh_856, \
                         nsh_857, nsh_858, nsh_859, osh_855, osh_856, osh_857, osh_858, \
                         osh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_13 * nsh_855[k]
                    + f_3 * pc_x[k] * osh_855[k];

        t_1136[k] = f_13 * nsh_856[k]
                    + f_3 * pc_x[k] * osh_856[k];

        t_1137[k] = f_13 * nsh_857[k]
                    + f_3 * pc_x[k] * osh_857[k];

        t_1138[k] = f_13 * nsh_858[k]
                    + f_3 * pc_x[k] * osh_858[k];

        t_1139[k] = f_13 * nsh_859[k]
                    + f_3 * pc_x[k] * osh_859[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, nsh_666, nsh_687, nsh_860, \
                         osg0_610, osg1_610, osh_855, osh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_13 * nsh_860[k]
                    + f_3 * pc_x[k] * osh_860[k];

        t_1141[k] = f_14 * nsh_687[k]
                    + f_1 * osg0_610[k]
                    - f_2 * osg1_610[k]
                    + f_3 * pc_y[k] * osh_855[k];

        t_1142[k] = f_14 * nsh_666[k]
                    + f_3 * pc_z[k] * osh_855[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_y, nsh_689, nsh_690, nsh_691, osg0_612, \
                         osg0_613, osg0_614, osg1_612, osg1_613, osg1_614, osh_857, osh_858, \
                         osh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_14 * nsh_689[k]
                    + f_8 * osg0_612[k]
                    - f_9 * osg1_612[k]
                    + f_3 * pc_y[k] * osh_857[k];

        t_1144[k] = f_14 * nsh_690[k]
                    + f_6 * osg0_613[k]
                    - f_7 * osg1_613[k]
                    + f_3 * pc_y[k] * osh_858[k];

        t_1145[k] = f_14 * nsh_691[k]
                    + f_4 * osg0_614[k]
                    - f_5 * osg1_614[k]
                    + f_3 * pc_y[k] * osh_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_y, pc_z, nsh_671, nsh_692, nsh_861, \
                         osg0_614, osg0_615, osg1_614, osg1_615, osh_860, \
                         osh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * nsh_692[k]
                    + f_3 * pc_y[k] * osh_860[k];

        t_1147[k] = f_14 * nsh_671[k]
                    + f_1 * osg0_614[k]
                    - f_2 * osg1_614[k]
                    + f_3 * pc_z[k] * osh_860[k];

        t_1148[k] = f_13 * nsh_861[k]
                    + f_1 * osg0_615[k]
                    - f_2 * osg1_615[k]
                    + f_3 * pc_x[k] * osh_861[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osg0, const size_t osg1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_980 = buffer.data(nsi0 + 980);
    const auto *nsi0_983 = buffer.data(nsi0 + 983);
    const auto *nsi0_985 = buffer.data(nsi0 + 985);
    const auto *nsi0_986 = buffer.data(nsi0 + 986);
    const auto *nsi0_989 = buffer.data(nsi0 + 989);
    const auto *nsi0_990 = buffer.data(nsi0 + 990);
    const auto *nsi0_992 = buffer.data(nsi0 + 992);
    const auto *nsi0_994 = buffer.data(nsi0 + 994);
    const auto *nsi0_1007 = buffer.data(nsi0 + 1007);

    const auto *nsh_672 = buffer.data(nsh + 672);
    const auto *nsh_675 = buffer.data(nsh + 675);
    const auto *nsh_678 = buffer.data(nsh + 678);
    const auto *nsh_687 = buffer.data(nsh + 687);
    const auto *nsh_692 = buffer.data(nsh + 692);
    const auto *nsh_693 = buffer.data(nsh + 693);
    const auto *nsh_695 = buffer.data(nsh + 695);
    const auto *nsh_696 = buffer.data(nsh + 696);
    const auto *nsh_698 = buffer.data(nsh + 698);
    const auto *nsh_699 = buffer.data(nsh + 699);
    const auto *nsh_702 = buffer.data(nsh + 702);
    const auto *nsh_708 = buffer.data(nsh + 708);
    const auto *nsh_710 = buffer.data(nsh + 710);
    const auto *nsh_711 = buffer.data(nsh + 711);
    const auto *nsh_712 = buffer.data(nsh + 712);
    const auto *nsh_713 = buffer.data(nsh + 713);
    const auto *nsh_714 = buffer.data(nsh + 714);
    const auto *nsh_716 = buffer.data(nsh + 716);
    const auto *nsh_717 = buffer.data(nsh + 717);
    const auto *nsh_719 = buffer.data(nsh + 719);
    const auto *nsh_720 = buffer.data(nsh + 720);
    const auto *nsh_723 = buffer.data(nsh + 723);
    const auto *nsh_729 = buffer.data(nsh + 729);
    const auto *nsh_731 = buffer.data(nsh + 731);
    const auto *nsh_732 = buffer.data(nsh + 732);
    const auto *nsh_733 = buffer.data(nsh + 733);
    const auto *nsh_734 = buffer.data(nsh + 734);
    const auto *nsh_735 = buffer.data(nsh + 735);
    const auto *nsh_736 = buffer.data(nsh + 736);
    const auto *nsh_737 = buffer.data(nsh + 737);
    const auto *nsh_738 = buffer.data(nsh + 738);
    const auto *nsh_740 = buffer.data(nsh + 740);
    const auto *nsh_741 = buffer.data(nsh + 741);
    const auto *nsh_743 = buffer.data(nsh + 743);
    const auto *nsh_744 = buffer.data(nsh + 744);
    const auto *nsh_750 = buffer.data(nsh + 750);
    const auto *nsh_752 = buffer.data(nsh + 752);
    const auto *nsh_753 = buffer.data(nsh + 753);
    const auto *nsh_754 = buffer.data(nsh + 754);
    const auto *nsh_755 = buffer.data(nsh + 755);
    const auto *nsh_864 = buffer.data(nsh + 864);
    const auto *nsh_866 = buffer.data(nsh + 866);
    const auto *nsh_867 = buffer.data(nsh + 867);
    const auto *nsh_870 = buffer.data(nsh + 870);
    const auto *nsh_871 = buffer.data(nsh + 871);
    const auto *nsh_873 = buffer.data(nsh + 873);
    const auto *nsh_875 = buffer.data(nsh + 875);
    const auto *nsh_876 = buffer.data(nsh + 876);
    const auto *nsh_877 = buffer.data(nsh + 877);
    const auto *nsh_878 = buffer.data(nsh + 878);
    const auto *nsh_879 = buffer.data(nsh + 879);
    const auto *nsh_880 = buffer.data(nsh + 880);
    const auto *nsh_881 = buffer.data(nsh + 881);
    const auto *nsh_882 = buffer.data(nsh + 882);
    const auto *nsh_885 = buffer.data(nsh + 885);
    const auto *nsh_887 = buffer.data(nsh + 887);
    const auto *nsh_888 = buffer.data(nsh + 888);
    const auto *nsh_891 = buffer.data(nsh + 891);
    const auto *nsh_892 = buffer.data(nsh + 892);
    const auto *nsh_894 = buffer.data(nsh + 894);
    const auto *nsh_896 = buffer.data(nsh + 896);
    const auto *nsh_897 = buffer.data(nsh + 897);
    const auto *nsh_898 = buffer.data(nsh + 898);
    const auto *nsh_899 = buffer.data(nsh + 899);
    const auto *nsh_900 = buffer.data(nsh + 900);
    const auto *nsh_901 = buffer.data(nsh + 901);
    const auto *nsh_902 = buffer.data(nsh + 902);
    const auto *nsh_918 = buffer.data(nsh + 918);
    const auto *nsh_919 = buffer.data(nsh + 919);
    const auto *nsh_920 = buffer.data(nsh + 920);
    const auto *nsh_921 = buffer.data(nsh + 921);
    const auto *nsh_922 = buffer.data(nsh + 922);
    const auto *nsh_923 = buffer.data(nsh + 923);
    const auto *nsh_924 = buffer.data(nsh + 924);
    const auto *nsh_929 = buffer.data(nsh + 929);
    const auto *nsh_933 = buffer.data(nsh + 933);
    const auto *nsh_938 = buffer.data(nsh + 938);
    const auto *nsh_939 = buffer.data(nsh + 939);
    const auto *nsh_940 = buffer.data(nsh + 940);
    const auto *nsh_941 = buffer.data(nsh + 941);
    const auto *nsh_942 = buffer.data(nsh + 942);
    const auto *nsh_944 = buffer.data(nsh + 944);

    const auto *nsi1_980 = buffer.data(nsi1 + 980);
    const auto *nsi1_983 = buffer.data(nsi1 + 983);
    const auto *nsi1_985 = buffer.data(nsi1 + 985);
    const auto *nsi1_986 = buffer.data(nsi1 + 986);
    const auto *nsi1_989 = buffer.data(nsi1 + 989);
    const auto *nsi1_990 = buffer.data(nsi1 + 990);
    const auto *nsi1_992 = buffer.data(nsi1 + 992);
    const auto *nsi1_994 = buffer.data(nsi1 + 994);
    const auto *nsi1_1007 = buffer.data(nsi1 + 1007);

    const auto *osg0_618 = buffer.data(osg0 + 618);
    const auto *osg0_620 = buffer.data(osg0 + 620);
    const auto *osg0_621 = buffer.data(osg0 + 621);
    const auto *osg0_624 = buffer.data(osg0 + 624);
    const auto *osg0_625 = buffer.data(osg0 + 625);
    const auto *osg0_627 = buffer.data(osg0 + 627);
    const auto *osg0_628 = buffer.data(osg0 + 628);
    const auto *osg0_629 = buffer.data(osg0 + 629);
    const auto *osg0_630 = buffer.data(osg0 + 630);
    const auto *osg0_633 = buffer.data(osg0 + 633);
    const auto *osg0_635 = buffer.data(osg0 + 635);
    const auto *osg0_636 = buffer.data(osg0 + 636);
    const auto *osg0_639 = buffer.data(osg0 + 639);
    const auto *osg0_640 = buffer.data(osg0 + 640);
    const auto *osg0_642 = buffer.data(osg0 + 642);
    const auto *osg0_643 = buffer.data(osg0 + 643);
    const auto *osg0_644 = buffer.data(osg0 + 644);
    const auto *osg0_655 = buffer.data(osg0 + 655);
    const auto *osg0_657 = buffer.data(osg0 + 657);
    const auto *osg0_658 = buffer.data(osg0 + 658);
    const auto *osg0_659 = buffer.data(osg0 + 659);
    const auto *osg0_660 = buffer.data(osg0 + 660);
    const auto *osg0_661 = buffer.data(osg0 + 661);
    const auto *osg0_662 = buffer.data(osg0 + 662);
    const auto *osg0_663 = buffer.data(osg0 + 663);
    const auto *osg0_664 = buffer.data(osg0 + 664);
    const auto *osg0_665 = buffer.data(osg0 + 665);
    const auto *osg0_669 = buffer.data(osg0 + 669);
    const auto *osg0_670 = buffer.data(osg0 + 670);
    const auto *osg0_671 = buffer.data(osg0 + 671);
    const auto *osg0_672 = buffer.data(osg0 + 672);
    const auto *osg0_673 = buffer.data(osg0 + 673);
    const auto *osg0_674 = buffer.data(osg0 + 674);

    const auto *osg1_618 = buffer.data(osg1 + 618);
    const auto *osg1_620 = buffer.data(osg1 + 620);
    const auto *osg1_621 = buffer.data(osg1 + 621);
    const auto *osg1_624 = buffer.data(osg1 + 624);
    const auto *osg1_625 = buffer.data(osg1 + 625);
    const auto *osg1_627 = buffer.data(osg1 + 627);
    const auto *osg1_628 = buffer.data(osg1 + 628);
    const auto *osg1_629 = buffer.data(osg1 + 629);
    const auto *osg1_630 = buffer.data(osg1 + 630);
    const auto *osg1_633 = buffer.data(osg1 + 633);
    const auto *osg1_635 = buffer.data(osg1 + 635);
    const auto *osg1_636 = buffer.data(osg1 + 636);
    const auto *osg1_639 = buffer.data(osg1 + 639);
    const auto *osg1_640 = buffer.data(osg1 + 640);
    const auto *osg1_642 = buffer.data(osg1 + 642);
    const auto *osg1_643 = buffer.data(osg1 + 643);
    const auto *osg1_644 = buffer.data(osg1 + 644);
    const auto *osg1_655 = buffer.data(osg1 + 655);
    const auto *osg1_657 = buffer.data(osg1 + 657);
    const auto *osg1_658 = buffer.data(osg1 + 658);
    const auto *osg1_659 = buffer.data(osg1 + 659);
    const auto *osg1_660 = buffer.data(osg1 + 660);
    const auto *osg1_661 = buffer.data(osg1 + 661);
    const auto *osg1_662 = buffer.data(osg1 + 662);
    const auto *osg1_663 = buffer.data(osg1 + 663);
    const auto *osg1_664 = buffer.data(osg1 + 664);
    const auto *osg1_665 = buffer.data(osg1 + 665);
    const auto *osg1_669 = buffer.data(osg1 + 669);
    const auto *osg1_670 = buffer.data(osg1 + 670);
    const auto *osg1_671 = buffer.data(osg1 + 671);
    const auto *osg1_672 = buffer.data(osg1 + 672);
    const auto *osg1_673 = buffer.data(osg1 + 673);
    const auto *osg1_674 = buffer.data(osg1 + 674);

    const auto *osh_861 = buffer.data(osh + 861);
    const auto *osh_863 = buffer.data(osh + 863);
    const auto *osh_864 = buffer.data(osh + 864);
    const auto *osh_866 = buffer.data(osh + 866);
    const auto *osh_867 = buffer.data(osh + 867);
    const auto *osh_870 = buffer.data(osh + 870);
    const auto *osh_871 = buffer.data(osh + 871);
    const auto *osh_873 = buffer.data(osh + 873);
    const auto *osh_875 = buffer.data(osh + 875);
    const auto *osh_876 = buffer.data(osh + 876);
    const auto *osh_877 = buffer.data(osh + 877);
    const auto *osh_878 = buffer.data(osh + 878);
    const auto *osh_879 = buffer.data(osh + 879);
    const auto *osh_880 = buffer.data(osh + 880);
    const auto *osh_881 = buffer.data(osh + 881);
    const auto *osh_882 = buffer.data(osh + 882);
    const auto *osh_884 = buffer.data(osh + 884);
    const auto *osh_885 = buffer.data(osh + 885);
    const auto *osh_887 = buffer.data(osh + 887);
    const auto *osh_888 = buffer.data(osh + 888);
    const auto *osh_891 = buffer.data(osh + 891);
    const auto *osh_892 = buffer.data(osh + 892);
    const auto *osh_894 = buffer.data(osh + 894);
    const auto *osh_896 = buffer.data(osh + 896);
    const auto *osh_897 = buffer.data(osh + 897);
    const auto *osh_898 = buffer.data(osh + 898);
    const auto *osh_899 = buffer.data(osh + 899);
    const auto *osh_900 = buffer.data(osh + 900);
    const auto *osh_901 = buffer.data(osh + 901);
    const auto *osh_902 = buffer.data(osh + 902);
    const auto *osh_903 = buffer.data(osh + 903);
    const auto *osh_905 = buffer.data(osh + 905);
    const auto *osh_906 = buffer.data(osh + 906);
    const auto *osh_908 = buffer.data(osh + 908);
    const auto *osh_909 = buffer.data(osh + 909);
    const auto *osh_912 = buffer.data(osh + 912);
    const auto *osh_918 = buffer.data(osh + 918);
    const auto *osh_919 = buffer.data(osh + 919);
    const auto *osh_920 = buffer.data(osh + 920);
    const auto *osh_921 = buffer.data(osh + 921);
    const auto *osh_922 = buffer.data(osh + 922);
    const auto *osh_923 = buffer.data(osh + 923);
    const auto *osh_924 = buffer.data(osh + 924);
    const auto *osh_925 = buffer.data(osh + 925);
    const auto *osh_926 = buffer.data(osh + 926);
    const auto *osh_927 = buffer.data(osh + 927);
    const auto *osh_928 = buffer.data(osh + 928);
    const auto *osh_929 = buffer.data(osh + 929);
    const auto *osh_930 = buffer.data(osh + 930);
    const auto *osh_931 = buffer.data(osh + 931);
    const auto *osh_932 = buffer.data(osh + 932);
    const auto *osh_933 = buffer.data(osh + 933);
    const auto *osh_938 = buffer.data(osh + 938);
    const auto *osh_939 = buffer.data(osh + 939);
    const auto *osh_940 = buffer.data(osh + 940);
    const auto *osh_941 = buffer.data(osh + 941);
    const auto *osh_942 = buffer.data(osh + 942);
    const auto *osh_943 = buffer.data(osh + 943);
    const auto *osh_944 = buffer.data(osh + 944);

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, nsh_672, nsh_693, \
                         nsh_695, nsh_864, osg0_618, osg1_618, osh_861, osh_863, \
                         osh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_13 * nsh_693[k]
                    + f_3 * pc_y[k] * osh_861[k];

        t_1150[k] = f_22 * nsh_672[k]
                    + f_3 * pc_z[k] * osh_861[k];

        t_1151[k] = f_13 * nsh_864[k]
                    + f_8 * osg0_618[k]
                    - f_9 * osg1_618[k]
                    + f_3 * pc_x[k] * osh_864[k];

        t_1152[k] = f_13 * nsh_695[k]
                    + f_3 * pc_y[k] * osh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pc_x, pc_z, nsh_675, nsh_866, nsh_867, \
                         osg0_620, osg0_621, osg1_620, osg1_621, osh_864, osh_866, \
                         osh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_13 * nsh_866[k]
                    + f_8 * osg0_620[k]
                    - f_9 * osg1_620[k]
                    + f_3 * pc_x[k] * osh_866[k];

        t_1154[k] = f_13 * nsh_867[k]
                    + f_6 * osg0_621[k]
                    - f_7 * osg1_621[k]
                    + f_3 * pc_x[k] * osh_867[k];

        t_1155[k] = f_22 * nsh_675[k]
                    + f_3 * pc_z[k] * osh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pc_x, pc_y, nsh_698, nsh_870, nsh_871, \
                         osg0_624, osg0_625, osg1_624, osg1_625, osh_866, osh_870, \
                         osh_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * nsh_698[k]
                    + f_3 * pc_y[k] * osh_866[k];

        t_1157[k] = f_13 * nsh_870[k]
                    + f_6 * osg0_624[k]
                    - f_7 * osg1_624[k]
                    + f_3 * pc_x[k] * osh_870[k];

        t_1158[k] = f_13 * nsh_871[k]
                    + f_4 * osg0_625[k]
                    - f_5 * osg1_625[k]
                    + f_3 * pc_x[k] * osh_871[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pc_x, pc_y, pc_z, nsh_678, nsh_702, nsh_873, \
                         osg0_627, osg1_627, osh_867, osh_870, \
                         osh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_22 * nsh_678[k]
                    + f_3 * pc_z[k] * osh_867[k];

        t_1160[k] = f_13 * nsh_873[k]
                    + f_4 * osg0_627[k]
                    - f_5 * osg1_627[k]
                    + f_3 * pc_x[k] * osh_873[k];

        t_1161[k] = f_13 * nsh_702[k]
                    + f_3 * pc_y[k] * osh_870[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pc_x, nsh_875, nsh_876, nsh_877, \
                         nsh_878, osg0_629, osg1_629, osh_875, osh_876, osh_877, \
                         osh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = f_13 * nsh_875[k]
                    + f_4 * osg0_629[k]
                    - f_5 * osg1_629[k]
                    + f_3 * pc_x[k] * osh_875[k];

        t_1163[k] = f_13 * nsh_876[k]
                    + f_3 * pc_x[k] * osh_876[k];

        t_1164[k] = f_13 * nsh_877[k]
                    + f_3 * pc_x[k] * osh_877[k];

        t_1165[k] = f_13 * nsh_878[k]
                    + f_3 * pc_x[k] * osh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pc_x, pc_y, nsh_708, nsh_879, \
                         nsh_880, nsh_881, osg0_625, osg1_625, osh_876, osh_879, osh_880, \
                         osh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_13 * nsh_879[k]
                    + f_3 * pc_x[k] * osh_879[k];

        t_1167[k] = f_13 * nsh_880[k]
                    + f_3 * pc_x[k] * osh_880[k];

        t_1168[k] = f_13 * nsh_881[k]
                    + f_3 * pc_x[k] * osh_881[k];

        t_1169[k] = f_13 * nsh_708[k]
                    + f_1 * osg0_625[k]
                    - f_2 * osg1_625[k]
                    + f_3 * pc_y[k] * osh_876[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pc_y, pc_z, nsh_687, nsh_710, nsh_711, \
                         osg0_627, osg0_628, osg1_627, osg1_628, osh_876, osh_878, \
                         osh_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_22 * nsh_687[k]
                    + f_3 * pc_z[k] * osh_876[k];

        t_1171[k] = f_13 * nsh_710[k]
                    + f_8 * osg0_627[k]
                    - f_9 * osg1_627[k]
                    + f_3 * pc_y[k] * osh_878[k];

        t_1172[k] = f_13 * nsh_711[k]
                    + f_6 * osg0_628[k]
                    - f_7 * osg1_628[k]
                    + f_3 * pc_y[k] * osh_879[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pc_y, pc_z, nsh_692, nsh_712, nsh_713, \
                         osg0_629, osg1_629, osh_880, osh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_13 * nsh_712[k]
                    + f_4 * osg0_629[k]
                    - f_5 * osg1_629[k]
                    + f_3 * pc_y[k] * osh_880[k];

        t_1174[k] = f_13 * nsh_713[k]
                    + f_3 * pc_y[k] * osh_881[k];

        t_1175[k] = f_22 * nsh_692[k]
                    + f_1 * osg0_629[k]
                    - f_2 * osg1_629[k]
                    + f_3 * pc_z[k] * osh_881[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, nsh_693, nsh_714, nsh_882, \
                         osg0_630, osg1_630, osh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_13 * nsh_882[k]
                    + f_1 * osg0_630[k]
                    - f_2 * osg1_630[k]
                    + f_3 * pc_x[k] * osh_882[k];

        t_1177[k] = f_12 * nsh_714[k]
                    + f_3 * pc_y[k] * osh_882[k];

        t_1178[k] = f_21 * nsh_693[k]
                    + f_3 * pc_z[k] * osh_882[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pc_x, pc_y, nsh_716, nsh_885, nsh_887, \
                         osg0_633, osg0_635, osg1_633, osg1_635, osh_884, osh_885, \
                         osh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_13 * nsh_885[k]
                    + f_8 * osg0_633[k]
                    - f_9 * osg1_633[k]
                    + f_3 * pc_x[k] * osh_885[k];

        t_1180[k] = f_12 * nsh_716[k]
                    + f_3 * pc_y[k] * osh_884[k];

        t_1181[k] = f_13 * nsh_887[k]
                    + f_8 * osg0_635[k]
                    - f_9 * osg1_635[k]
                    + f_3 * pc_x[k] * osh_887[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pc_z, nsh_696, nsh_719, nsh_888, \
                         osg0_636, osg1_636, osh_885, osh_887, \
                         osh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_13 * nsh_888[k]
                    + f_6 * osg0_636[k]
                    - f_7 * osg1_636[k]
                    + f_3 * pc_x[k] * osh_888[k];

        t_1183[k] = f_21 * nsh_696[k]
                    + f_3 * pc_z[k] * osh_885[k];

        t_1184[k] = f_12 * nsh_719[k]
                    + f_3 * pc_y[k] * osh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, nsh_699, nsh_891, nsh_892, \
                         osg0_639, osg0_640, osg1_639, osg1_640, osh_888, osh_891, \
                         osh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_13 * nsh_891[k]
                    + f_6 * osg0_639[k]
                    - f_7 * osg1_639[k]
                    + f_3 * pc_x[k] * osh_891[k];

        t_1186[k] = f_13 * nsh_892[k]
                    + f_4 * osg0_640[k]
                    - f_5 * osg1_640[k]
                    + f_3 * pc_x[k] * osh_892[k];

        t_1187[k] = f_21 * nsh_699[k]
                    + f_3 * pc_z[k] * osh_888[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pc_x, pc_y, nsh_723, nsh_894, nsh_896, \
                         osg0_642, osg0_644, osg1_642, osg1_644, osh_891, osh_894, \
                         osh_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_13 * nsh_894[k]
                    + f_4 * osg0_642[k]
                    - f_5 * osg1_642[k]
                    + f_3 * pc_x[k] * osh_894[k];

        t_1189[k] = f_12 * nsh_723[k]
                    + f_3 * pc_y[k] * osh_891[k];

        t_1190[k] = f_13 * nsh_896[k]
                    + f_4 * osg0_644[k]
                    - f_5 * osg1_644[k]
                    + f_3 * pc_x[k] * osh_896[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, t_1195, pc_x, nsh_897, nsh_898, \
                         nsh_899, nsh_900, nsh_901, osh_897, osh_898, osh_899, osh_900, \
                         osh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_13 * nsh_897[k]
                    + f_3 * pc_x[k] * osh_897[k];

        t_1192[k] = f_13 * nsh_898[k]
                    + f_3 * pc_x[k] * osh_898[k];

        t_1193[k] = f_13 * nsh_899[k]
                    + f_3 * pc_x[k] * osh_899[k];

        t_1194[k] = f_13 * nsh_900[k]
                    + f_3 * pc_x[k] * osh_900[k];

        t_1195[k] = f_13 * nsh_901[k]
                    + f_3 * pc_x[k] * osh_901[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, nsh_708, nsh_729, nsh_902, \
                         osg0_640, osg1_640, osh_897, osh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_13 * nsh_902[k]
                    + f_3 * pc_x[k] * osh_902[k];

        t_1197[k] = f_12 * nsh_729[k]
                    + f_1 * osg0_640[k]
                    - f_2 * osg1_640[k]
                    + f_3 * pc_y[k] * osh_897[k];

        t_1198[k] = f_21 * nsh_708[k]
                    + f_3 * pc_z[k] * osh_897[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_y, nsh_731, nsh_732, nsh_733, osg0_642, \
                         osg0_643, osg0_644, osg1_642, osg1_643, osg1_644, osh_899, osh_900, \
                         osh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_12 * nsh_731[k]
                    + f_8 * osg0_642[k]
                    - f_9 * osg1_642[k]
                    + f_3 * pc_y[k] * osh_899[k];

        t_1200[k] = f_12 * nsh_732[k]
                    + f_6 * osg0_643[k]
                    - f_7 * osg1_643[k]
                    + f_3 * pc_y[k] * osh_900[k];

        t_1201[k] = f_12 * nsh_733[k]
                    + f_4 * osg0_644[k]
                    - f_5 * osg1_644[k]
                    + f_3 * pc_y[k] * osh_901[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pa_y, pc_y, pc_z, nsi0_980, nsh_713, \
                         nsh_734, nsh_735, nsi1_980, osg0_644, osg1_644, osh_902, \
                         osh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * nsh_734[k]
                    + f_3 * pc_y[k] * osh_902[k];

        t_1203[k] = f_21 * nsh_713[k]
                    + f_1 * osg0_644[k]
                    - f_2 * osg1_644[k]
                    + f_3 * pc_z[k] * osh_902[k];

        t_1204[k] = pa_y[k] * nsi0_980[k]
                    - f_10 * pc_y[k] * nsi1_980[k];

        t_1205[k] = f_11 * nsh_735[k]
                    + f_3 * pc_y[k] * osh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pa_y, pc_y, pc_z, nsi0_983, nsi0_985, \
                         nsh_714, nsh_736, nsh_737, nsi1_983, nsi1_985, osh_903, \
                         osh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_20 * nsh_714[k]
                    + f_3 * pc_z[k] * osh_903[k];

        t_1207[k] = pa_y[k] * nsi0_983[k]
                    + f_12 * nsh_736[k]
                    - f_10 * pc_y[k] * nsi1_983[k];

        t_1208[k] = f_11 * nsh_737[k]
                    + f_3 * pc_y[k] * osh_905[k];

        t_1209[k] = pa_y[k] * nsi0_985[k]
                    - f_10 * pc_y[k] * nsi1_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pa_y, pc_y, pc_z, nsi0_986, nsi0_989, \
                         nsh_717, nsh_738, nsh_740, nsi1_986, nsi1_989, osh_906, \
                         osh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pa_y[k] * nsi0_986[k]
                    + f_13 * nsh_738[k]
                    - f_10 * pc_y[k] * nsi1_986[k];

        t_1211[k] = f_20 * nsh_717[k]
                    + f_3 * pc_z[k] * osh_906[k];

        t_1212[k] = f_11 * nsh_740[k]
                    + f_3 * pc_y[k] * osh_908[k];

        t_1213[k] = pa_y[k] * nsi0_989[k]
                    - f_10 * pc_y[k] * nsi1_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pa_y, pc_y, pc_z, nsi0_990, nsi0_992, \
                         nsh_720, nsh_741, nsh_743, nsi1_990, nsi1_992, \
                         osh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_y[k] * nsi0_990[k]
                    + f_14 * nsh_741[k]
                    - f_10 * pc_y[k] * nsi1_990[k];

        t_1215[k] = f_20 * nsh_720[k]
                    + f_3 * pc_z[k] * osh_909[k];

        t_1216[k] = pa_y[k] * nsi0_992[k]
                    + f_12 * nsh_743[k]
                    - f_10 * pc_y[k] * nsi1_992[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, pa_y, pc_x, pc_y, nsi0_994, nsh_744, \
                         nsh_918, nsh_919, nsi1_994, osh_912, osh_918, \
                         osh_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_11 * nsh_744[k]
                    + f_3 * pc_y[k] * osh_912[k];

        t_1218[k] = pa_y[k] * nsi0_994[k]
                    - f_10 * pc_y[k] * nsi1_994[k];

        t_1219[k] = f_13 * nsh_918[k]
                    + f_3 * pc_x[k] * osh_918[k];

        t_1220[k] = f_13 * nsh_919[k]
                    + f_3 * pc_x[k] * osh_919[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, pc_x, nsh_920, nsh_921, nsh_922, \
                         nsh_923, osh_920, osh_921, osh_922, osh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_13 * nsh_920[k]
                    + f_3 * pc_x[k] * osh_920[k];

        t_1222[k] = f_13 * nsh_921[k]
                    + f_3 * pc_x[k] * osh_921[k];

        t_1223[k] = f_13 * nsh_922[k]
                    + f_3 * pc_x[k] * osh_922[k];

        t_1224[k] = f_13 * nsh_923[k]
                    + f_3 * pc_x[k] * osh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pc_y, pc_z, nsh_729, nsh_750, nsh_752, \
                         osg0_655, osg0_657, osg1_655, osg1_657, osh_918, \
                         osh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_11 * nsh_750[k]
                    + f_1 * osg0_655[k]
                    - f_2 * osg1_655[k]
                    + f_3 * pc_y[k] * osh_918[k];

        t_1226[k] = f_20 * nsh_729[k]
                    + f_3 * pc_z[k] * osh_918[k];

        t_1227[k] = f_11 * nsh_752[k]
                    + f_8 * osg0_657[k]
                    - f_9 * osg1_657[k]
                    + f_3 * pc_y[k] * osh_920[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, pc_y, nsh_753, nsh_754, nsh_755, osg0_658, \
                         osg0_659, osg1_658, osg1_659, osh_921, osh_922, \
                         osh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = f_11 * nsh_753[k]
                    + f_6 * osg0_658[k]
                    - f_7 * osg1_658[k]
                    + f_3 * pc_y[k] * osh_921[k];

        t_1229[k] = f_11 * nsh_754[k]
                    + f_4 * osg0_659[k]
                    - f_5 * osg1_659[k]
                    + f_3 * pc_y[k] * osh_922[k];

        t_1230[k] = f_11 * nsh_755[k]
                    + f_3 * pc_y[k] * osh_923[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pa_y, pc_x, pc_y, pc_z, nsi0_1007, \
                         nsh_735, nsh_924, nsi1_1007, osg0_660, osg1_660, \
                         osh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = pa_y[k] * nsi0_1007[k]
                    - f_10 * pc_y[k] * nsi1_1007[k];

        t_1232[k] = f_13 * nsh_924[k]
                    + f_1 * osg0_660[k]
                    - f_2 * osg1_660[k]
                    + f_3 * pc_x[k] * osh_924[k];

        t_1233[k] = f_3 * pc_y[k] * osh_924[k];

        t_1234[k] = f_19 * nsh_735[k]
                    + f_3 * pc_z[k] * osh_924[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, pc_y, nsh_929, osg0_660, osg0_665, \
                         osg1_660, osg1_665, osh_925, osh_926, \
                         osh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_4 * osg0_660[k]
                    - f_5 * osg1_660[k]
                    + f_3 * pc_y[k] * osh_925[k];

        t_1236[k] = f_3 * pc_y[k] * osh_926[k];

        t_1237[k] = f_13 * nsh_929[k]
                    + f_8 * osg0_665[k]
                    - f_9 * osg1_665[k]
                    + f_3 * pc_x[k] * osh_929[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pc_y, osg0_661, osg0_662, osg1_661, osg1_662, \
                         osh_927, osh_928, osh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_6 * osg0_661[k]
                    - f_7 * osg1_661[k]
                    + f_3 * pc_y[k] * osh_927[k];

        t_1239[k] = f_4 * osg0_662[k]
                    - f_5 * osg1_662[k]
                    + f_3 * pc_y[k] * osh_928[k];

        t_1240[k] = f_3 * pc_y[k] * osh_929[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_y, nsh_933, osg0_663, osg0_664, \
                         osg0_669, osg1_663, osg1_664, osg1_669, osh_930, osh_931, \
                         osh_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_13 * nsh_933[k]
                    + f_6 * osg0_669[k]
                    - f_7 * osg1_669[k]
                    + f_3 * pc_x[k] * osh_933[k];

        t_1242[k] = f_8 * osg0_663[k]
                    - f_9 * osg1_663[k]
                    + f_3 * pc_y[k] * osh_930[k];

        t_1243[k] = f_6 * osg0_664[k]
                    - f_7 * osg1_664[k]
                    + f_3 * pc_y[k] * osh_931[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, pc_x, pc_y, nsh_938, nsh_939, \
                         osg0_665, osg0_674, osg1_665, osg1_674, osh_932, osh_933, osh_938, \
                         osh_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_4 * osg0_665[k]
                    - f_5 * osg1_665[k]
                    + f_3 * pc_y[k] * osh_932[k];

        t_1245[k] = f_3 * pc_y[k] * osh_933[k];

        t_1246[k] = f_13 * nsh_938[k]
                    + f_4 * osg0_674[k]
                    - f_5 * osg1_674[k]
                    + f_3 * pc_x[k] * osh_938[k];

        t_1247[k] = f_13 * nsh_939[k]
                    + f_3 * pc_x[k] * osh_939[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, pc_x, pc_y, nsh_940, nsh_941, \
                         nsh_942, nsh_944, osh_938, osh_940, osh_941, osh_942, \
                         osh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_13 * nsh_940[k]
                    + f_3 * pc_x[k] * osh_940[k];

        t_1249[k] = f_13 * nsh_941[k]
                    + f_3 * pc_x[k] * osh_941[k];

        t_1250[k] = f_13 * nsh_942[k]
                    + f_3 * pc_x[k] * osh_942[k];

        t_1251[k] = f_3 * pc_y[k] * osh_938[k];

        t_1252[k] = f_13 * nsh_944[k]
                    + f_3 * pc_x[k] * osh_944[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, pc_y, osg0_670, osg0_671, osg0_672, osg1_670, \
                         osg1_671, osg1_672, osh_939, osh_940, \
                         osh_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = f_1 * osg0_670[k]
                    - f_2 * osg1_670[k]
                    + f_3 * pc_y[k] * osh_939[k];

        t_1254[k] = f_16 * osg0_671[k]
                    - f_17 * osg1_671[k]
                    + f_3 * pc_y[k] * osh_940[k];

        t_1255[k] = f_8 * osg0_672[k]
                    - f_9 * osg1_672[k]
                    + f_3 * pc_y[k] * osh_941[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, nsh_755, osg0_673, \
                         osg0_674, osg1_673, osg1_674, osh_942, osh_943, \
                         osh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_6 * osg0_673[k]
                    - f_7 * osg1_673[k]
                    + f_3 * pc_y[k] * osh_942[k];

        t_1257[k] = f_4 * osg0_674[k]
                    - f_5 * osg1_674[k]
                    + f_3 * pc_y[k] * osh_943[k];

        t_1258[k] = f_3 * pc_y[k] * osh_944[k];

        t_1259[k] = f_19 * nsh_755[k]
                    + f_1 * osg0_674[k]
                    - f_2 * osg1_674[k]
                    + f_3 * pc_z[k] * osh_944[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osg0, const size_t osg1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_1008 = buffer.data(nsi0 + 1008);
    const auto *nsi0_1011 = buffer.data(nsi0 + 1011);
    const auto *nsi0_1014 = buffer.data(nsi0 + 1014);
    const auto *nsi0_1018 = buffer.data(nsi0 + 1018);
    const auto *nsi0_1020 = buffer.data(nsi0 + 1020);
    const auto *nsi0_1029 = buffer.data(nsi0 + 1029);

    const auto *nsh_756 = buffer.data(nsh + 756);
    const auto *nsh_759 = buffer.data(nsh + 759);
    const auto *nsh_761 = buffer.data(nsh + 761);
    const auto *nsh_762 = buffer.data(nsh + 762);
    const auto *nsh_763 = buffer.data(nsh + 763);
    const auto *nsh_765 = buffer.data(nsh + 765);
    const auto *nsh_771 = buffer.data(nsh + 771);
    const auto *nsh_776 = buffer.data(nsh + 776);
    const auto *nsh_777 = buffer.data(nsh + 777);
    const auto *nsh_779 = buffer.data(nsh + 779);
    const auto *nsh_780 = buffer.data(nsh + 780);
    const auto *nsh_782 = buffer.data(nsh + 782);
    const auto *nsh_783 = buffer.data(nsh + 783);
    const auto *nsh_786 = buffer.data(nsh + 786);
    const auto *nsh_792 = buffer.data(nsh + 792);
    const auto *nsh_794 = buffer.data(nsh + 794);
    const auto *nsh_795 = buffer.data(nsh + 795);
    const auto *nsh_796 = buffer.data(nsh + 796);
    const auto *nsh_797 = buffer.data(nsh + 797);
    const auto *nsh_798 = buffer.data(nsh + 798);
    const auto *nsh_800 = buffer.data(nsh + 800);
    const auto *nsh_801 = buffer.data(nsh + 801);
    const auto *nsh_803 = buffer.data(nsh + 803);
    const auto *nsh_804 = buffer.data(nsh + 804);
    const auto *nsh_807 = buffer.data(nsh + 807);
    const auto *nsh_813 = buffer.data(nsh + 813);
    const auto *nsh_815 = buffer.data(nsh + 815);
    const auto *nsh_816 = buffer.data(nsh + 816);
    const auto *nsh_817 = buffer.data(nsh + 817);
    const auto *nsh_818 = buffer.data(nsh + 818);
    const auto *nsh_819 = buffer.data(nsh + 819);
    const auto *nsh_821 = buffer.data(nsh + 821);
    const auto *nsh_824 = buffer.data(nsh + 824);
    const auto *nsh_828 = buffer.data(nsh + 828);
    const auto *nsh_834 = buffer.data(nsh + 834);
    const auto *nsh_836 = buffer.data(nsh + 836);
    const auto *nsh_837 = buffer.data(nsh + 837);
    const auto *nsh_838 = buffer.data(nsh + 838);
    const auto *nsh_839 = buffer.data(nsh + 839);
    const auto *nsh_945 = buffer.data(nsh + 945);
    const auto *nsh_948 = buffer.data(nsh + 948);
    const auto *nsh_951 = buffer.data(nsh + 951);
    const auto *nsh_955 = buffer.data(nsh + 955);
    const auto *nsh_960 = buffer.data(nsh + 960);
    const auto *nsh_962 = buffer.data(nsh + 962);
    const auto *nsh_963 = buffer.data(nsh + 963);
    const auto *nsh_964 = buffer.data(nsh + 964);
    const auto *nsh_965 = buffer.data(nsh + 965);
    const auto *nsh_971 = buffer.data(nsh + 971);
    const auto *nsh_975 = buffer.data(nsh + 975);
    const auto *nsh_980 = buffer.data(nsh + 980);
    const auto *nsh_981 = buffer.data(nsh + 981);
    const auto *nsh_982 = buffer.data(nsh + 982);
    const auto *nsh_983 = buffer.data(nsh + 983);
    const auto *nsh_984 = buffer.data(nsh + 984);
    const auto *nsh_985 = buffer.data(nsh + 985);
    const auto *nsh_986 = buffer.data(nsh + 986);
    const auto *nsh_987 = buffer.data(nsh + 987);
    const auto *nsh_990 = buffer.data(nsh + 990);
    const auto *nsh_992 = buffer.data(nsh + 992);
    const auto *nsh_993 = buffer.data(nsh + 993);
    const auto *nsh_996 = buffer.data(nsh + 996);
    const auto *nsh_997 = buffer.data(nsh + 997);
    const auto *nsh_999 = buffer.data(nsh + 999);
    const auto *nsh_1001 = buffer.data(nsh + 1001);
    const auto *nsh_1002 = buffer.data(nsh + 1002);
    const auto *nsh_1003 = buffer.data(nsh + 1003);
    const auto *nsh_1004 = buffer.data(nsh + 1004);
    const auto *nsh_1005 = buffer.data(nsh + 1005);
    const auto *nsh_1006 = buffer.data(nsh + 1006);
    const auto *nsh_1007 = buffer.data(nsh + 1007);
    const auto *nsh_1008 = buffer.data(nsh + 1008);
    const auto *nsh_1011 = buffer.data(nsh + 1011);
    const auto *nsh_1013 = buffer.data(nsh + 1013);
    const auto *nsh_1014 = buffer.data(nsh + 1014);
    const auto *nsh_1017 = buffer.data(nsh + 1017);
    const auto *nsh_1018 = buffer.data(nsh + 1018);
    const auto *nsh_1020 = buffer.data(nsh + 1020);
    const auto *nsh_1022 = buffer.data(nsh + 1022);
    const auto *nsh_1023 = buffer.data(nsh + 1023);
    const auto *nsh_1024 = buffer.data(nsh + 1024);
    const auto *nsh_1025 = buffer.data(nsh + 1025);
    const auto *nsh_1026 = buffer.data(nsh + 1026);
    const auto *nsh_1027 = buffer.data(nsh + 1027);
    const auto *nsh_1028 = buffer.data(nsh + 1028);

    const auto *nsi1_1008 = buffer.data(nsi1 + 1008);
    const auto *nsi1_1011 = buffer.data(nsi1 + 1011);
    const auto *nsi1_1014 = buffer.data(nsi1 + 1014);
    const auto *nsi1_1018 = buffer.data(nsi1 + 1018);
    const auto *nsi1_1020 = buffer.data(nsi1 + 1020);
    const auto *nsi1_1029 = buffer.data(nsi1 + 1029);

    const auto *osg0_675 = buffer.data(osg0 + 675);
    const auto *osg0_677 = buffer.data(osg0 + 677);
    const auto *osg0_678 = buffer.data(osg0 + 678);
    const auto *osg0_680 = buffer.data(osg0 + 680);
    const auto *osg0_681 = buffer.data(osg0 + 681);
    const auto *osg0_685 = buffer.data(osg0 + 685);
    const auto *osg0_686 = buffer.data(osg0 + 686);
    const auto *osg0_687 = buffer.data(osg0 + 687);
    const auto *osg0_689 = buffer.data(osg0 + 689);
    const auto *osg0_695 = buffer.data(osg0 + 695);
    const auto *osg0_699 = buffer.data(osg0 + 699);
    const auto *osg0_702 = buffer.data(osg0 + 702);
    const auto *osg0_703 = buffer.data(osg0 + 703);
    const auto *osg0_704 = buffer.data(osg0 + 704);
    const auto *osg0_705 = buffer.data(osg0 + 705);
    const auto *osg0_708 = buffer.data(osg0 + 708);
    const auto *osg0_710 = buffer.data(osg0 + 710);
    const auto *osg0_711 = buffer.data(osg0 + 711);
    const auto *osg0_714 = buffer.data(osg0 + 714);
    const auto *osg0_715 = buffer.data(osg0 + 715);
    const auto *osg0_717 = buffer.data(osg0 + 717);
    const auto *osg0_718 = buffer.data(osg0 + 718);
    const auto *osg0_719 = buffer.data(osg0 + 719);
    const auto *osg0_720 = buffer.data(osg0 + 720);
    const auto *osg0_723 = buffer.data(osg0 + 723);
    const auto *osg0_725 = buffer.data(osg0 + 725);
    const auto *osg0_726 = buffer.data(osg0 + 726);
    const auto *osg0_729 = buffer.data(osg0 + 729);
    const auto *osg0_730 = buffer.data(osg0 + 730);
    const auto *osg0_732 = buffer.data(osg0 + 732);
    const auto *osg0_733 = buffer.data(osg0 + 733);
    const auto *osg0_734 = buffer.data(osg0 + 734);

    const auto *osg1_675 = buffer.data(osg1 + 675);
    const auto *osg1_677 = buffer.data(osg1 + 677);
    const auto *osg1_678 = buffer.data(osg1 + 678);
    const auto *osg1_680 = buffer.data(osg1 + 680);
    const auto *osg1_681 = buffer.data(osg1 + 681);
    const auto *osg1_685 = buffer.data(osg1 + 685);
    const auto *osg1_686 = buffer.data(osg1 + 686);
    const auto *osg1_687 = buffer.data(osg1 + 687);
    const auto *osg1_689 = buffer.data(osg1 + 689);
    const auto *osg1_695 = buffer.data(osg1 + 695);
    const auto *osg1_699 = buffer.data(osg1 + 699);
    const auto *osg1_702 = buffer.data(osg1 + 702);
    const auto *osg1_703 = buffer.data(osg1 + 703);
    const auto *osg1_704 = buffer.data(osg1 + 704);
    const auto *osg1_705 = buffer.data(osg1 + 705);
    const auto *osg1_708 = buffer.data(osg1 + 708);
    const auto *osg1_710 = buffer.data(osg1 + 710);
    const auto *osg1_711 = buffer.data(osg1 + 711);
    const auto *osg1_714 = buffer.data(osg1 + 714);
    const auto *osg1_715 = buffer.data(osg1 + 715);
    const auto *osg1_717 = buffer.data(osg1 + 717);
    const auto *osg1_718 = buffer.data(osg1 + 718);
    const auto *osg1_719 = buffer.data(osg1 + 719);
    const auto *osg1_720 = buffer.data(osg1 + 720);
    const auto *osg1_723 = buffer.data(osg1 + 723);
    const auto *osg1_725 = buffer.data(osg1 + 725);
    const auto *osg1_726 = buffer.data(osg1 + 726);
    const auto *osg1_729 = buffer.data(osg1 + 729);
    const auto *osg1_730 = buffer.data(osg1 + 730);
    const auto *osg1_732 = buffer.data(osg1 + 732);
    const auto *osg1_733 = buffer.data(osg1 + 733);
    const auto *osg1_734 = buffer.data(osg1 + 734);

    const auto *osh_945 = buffer.data(osh + 945);
    const auto *osh_946 = buffer.data(osh + 946);
    const auto *osh_947 = buffer.data(osh + 947);
    const auto *osh_948 = buffer.data(osh + 948);
    const auto *osh_950 = buffer.data(osh + 950);
    const auto *osh_951 = buffer.data(osh + 951);
    const auto *osh_952 = buffer.data(osh + 952);
    const auto *osh_954 = buffer.data(osh + 954);
    const auto *osh_955 = buffer.data(osh + 955);
    const auto *osh_960 = buffer.data(osh + 960);
    const auto *osh_961 = buffer.data(osh + 961);
    const auto *osh_962 = buffer.data(osh + 962);
    const auto *osh_963 = buffer.data(osh + 963);
    const auto *osh_964 = buffer.data(osh + 964);
    const auto *osh_965 = buffer.data(osh + 965);
    const auto *osh_966 = buffer.data(osh + 966);
    const auto *osh_968 = buffer.data(osh + 968);
    const auto *osh_969 = buffer.data(osh + 969);
    const auto *osh_971 = buffer.data(osh + 971);
    const auto *osh_972 = buffer.data(osh + 972);
    const auto *osh_975 = buffer.data(osh + 975);
    const auto *osh_980 = buffer.data(osh + 980);
    const auto *osh_981 = buffer.data(osh + 981);
    const auto *osh_982 = buffer.data(osh + 982);
    const auto *osh_983 = buffer.data(osh + 983);
    const auto *osh_984 = buffer.data(osh + 984);
    const auto *osh_985 = buffer.data(osh + 985);
    const auto *osh_986 = buffer.data(osh + 986);
    const auto *osh_987 = buffer.data(osh + 987);
    const auto *osh_989 = buffer.data(osh + 989);
    const auto *osh_990 = buffer.data(osh + 990);
    const auto *osh_992 = buffer.data(osh + 992);
    const auto *osh_993 = buffer.data(osh + 993);
    const auto *osh_996 = buffer.data(osh + 996);
    const auto *osh_997 = buffer.data(osh + 997);
    const auto *osh_999 = buffer.data(osh + 999);
    const auto *osh_1001 = buffer.data(osh + 1001);
    const auto *osh_1002 = buffer.data(osh + 1002);
    const auto *osh_1003 = buffer.data(osh + 1003);
    const auto *osh_1004 = buffer.data(osh + 1004);
    const auto *osh_1005 = buffer.data(osh + 1005);
    const auto *osh_1006 = buffer.data(osh + 1006);
    const auto *osh_1007 = buffer.data(osh + 1007);
    const auto *osh_1008 = buffer.data(osh + 1008);
    const auto *osh_1010 = buffer.data(osh + 1010);
    const auto *osh_1011 = buffer.data(osh + 1011);
    const auto *osh_1013 = buffer.data(osh + 1013);
    const auto *osh_1014 = buffer.data(osh + 1014);
    const auto *osh_1017 = buffer.data(osh + 1017);
    const auto *osh_1018 = buffer.data(osh + 1018);
    const auto *osh_1020 = buffer.data(osh + 1020);
    const auto *osh_1022 = buffer.data(osh + 1022);
    const auto *osh_1023 = buffer.data(osh + 1023);
    const auto *osh_1024 = buffer.data(osh + 1024);
    const auto *osh_1025 = buffer.data(osh + 1025);
    const auto *osh_1026 = buffer.data(osh + 1026);
    const auto *osh_1027 = buffer.data(osh + 1027);
    const auto *osh_1028 = buffer.data(osh + 1028);

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pc_x, pc_y, pc_z, nsh_756, nsh_945, \
                         nsh_948, osg0_675, osg0_678, osg1_675, osg1_678, osh_945, \
                         osh_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_12 * nsh_945[k]
                    + f_1 * osg0_675[k]
                    - f_2 * osg1_675[k]
                    + f_3 * pc_x[k] * osh_945[k];

        t_1261[k] = f_18 * nsh_756[k]
                    + f_3 * pc_y[k] * osh_945[k];

        t_1262[k] = f_3 * pc_z[k] * osh_945[k];

        t_1263[k] = f_12 * nsh_948[k]
                    + f_8 * osg0_678[k]
                    - f_9 * osg1_678[k]
                    + f_3 * pc_x[k] * osh_948[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, t_1267, pc_x, pc_z, nsh_951, osg0_675, \
                         osg0_681, osg1_675, osg1_681, osh_946, osh_947, osh_948, \
                         osh_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_3 * pc_z[k] * osh_946[k];

        t_1265[k] = f_4 * osg0_675[k]
                    - f_5 * osg1_675[k]
                    + f_3 * pc_z[k] * osh_947[k];

        t_1266[k] = f_12 * nsh_951[k]
                    + f_6 * osg0_681[k]
                    - f_7 * osg1_681[k]
                    + f_3 * pc_x[k] * osh_951[k];

        t_1267[k] = f_3 * pc_z[k] * osh_948[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pc_x, pc_y, pc_z, nsh_761, nsh_955, \
                         osg0_677, osg0_685, osg1_677, osg1_685, osh_950, osh_951, \
                         osh_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_18 * nsh_761[k]
                    + f_3 * pc_y[k] * osh_950[k];

        t_1269[k] = f_6 * osg0_677[k]
                    - f_7 * osg1_677[k]
                    + f_3 * pc_z[k] * osh_950[k];

        t_1270[k] = f_12 * nsh_955[k]
                    + f_4 * osg0_685[k]
                    - f_5 * osg1_685[k]
                    + f_3 * pc_x[k] * osh_955[k];

        t_1271[k] = f_3 * pc_z[k] * osh_951[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, pc_z, nsh_765, nsh_960, \
                         osg0_678, osg0_680, osg1_678, osg1_680, osh_952, osh_954, \
                         osh_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * osg0_678[k]
                    - f_5 * osg1_678[k]
                    + f_3 * pc_z[k] * osh_952[k];

        t_1273[k] = f_18 * nsh_765[k]
                    + f_3 * pc_y[k] * osh_954[k];

        t_1274[k] = f_8 * osg0_680[k]
                    - f_9 * osg1_680[k]
                    + f_3 * pc_z[k] * osh_954[k];

        t_1275[k] = f_12 * nsh_960[k]
                    + f_3 * pc_x[k] * osh_960[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, t_1279, t_1280, pc_x, pc_z, nsh_962, nsh_963, \
                         nsh_964, nsh_965, osh_955, osh_962, osh_963, osh_964, \
                         osh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_3 * pc_z[k] * osh_955[k];

        t_1277[k] = f_12 * nsh_962[k]
                    + f_3 * pc_x[k] * osh_962[k];

        t_1278[k] = f_12 * nsh_963[k]
                    + f_3 * pc_x[k] * osh_963[k];

        t_1279[k] = f_12 * nsh_964[k]
                    + f_3 * pc_x[k] * osh_964[k];

        t_1280[k] = f_12 * nsh_965[k]
                    + f_3 * pc_x[k] * osh_965[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, pc_y, pc_z, nsh_771, osg0_685, \
                         osg0_686, osg1_685, osg1_686, osh_960, osh_961, \
                         osh_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_18 * nsh_771[k]
                    + f_1 * osg0_685[k]
                    - f_2 * osg1_685[k]
                    + f_3 * pc_y[k] * osh_960[k];

        t_1282[k] = f_3 * pc_z[k] * osh_960[k];

        t_1283[k] = f_4 * osg0_685[k]
                    - f_5 * osg1_685[k]
                    + f_3 * pc_z[k] * osh_961[k];

        t_1284[k] = f_6 * osg0_686[k]
                    - f_7 * osg1_686[k]
                    + f_3 * pc_z[k] * osh_962[k];
    }

#pragma omp simd aligned(t_1285, t_1286, t_1287, t_1288, pa_z, pc_y, pc_z, nsi0_1008, nsh_776, \
                         nsi1_1008, osg0_687, osg0_689, osg1_687, osg1_689, osh_963, \
                         osh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1285[k] = f_8 * osg0_687[k]
                    - f_9 * osg1_687[k]
                    + f_3 * pc_z[k] * osh_963[k];

        t_1286[k] = f_18 * nsh_776[k]
                    + f_3 * pc_y[k] * osh_965[k];

        t_1287[k] = f_1 * osg0_689[k]
                    - f_2 * osg1_689[k]
                    + f_3 * pc_z[k] * osh_965[k];

        t_1288[k] = pa_z[k] * nsi0_1008[k]
                    - f_10 * pc_z[k] * nsi1_1008[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, pa_z, pc_y, pc_z, nsi0_1011, nsh_756, \
                         nsh_777, nsh_779, nsi1_1011, osh_966, \
                         osh_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_19 * nsh_777[k]
                    + f_3 * pc_y[k] * osh_966[k];

        t_1290[k] = f_11 * nsh_756[k]
                    + f_3 * pc_z[k] * osh_966[k];

        t_1291[k] = pa_z[k] * nsi0_1011[k]
                    - f_10 * pc_z[k] * nsi1_1011[k];

        t_1292[k] = f_19 * nsh_779[k]
                    + f_3 * pc_y[k] * osh_968[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, pa_z, pc_x, pc_z, nsi0_1014, nsh_759, \
                         nsh_971, nsi1_1014, osg0_695, osg1_695, osh_969, \
                         osh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_12 * nsh_971[k]
                    + f_8 * osg0_695[k]
                    - f_9 * osg1_695[k]
                    + f_3 * pc_x[k] * osh_971[k];

        t_1294[k] = pa_z[k] * nsi0_1014[k]
                    - f_10 * pc_z[k] * nsi1_1014[k];

        t_1295[k] = f_11 * nsh_759[k]
                    + f_3 * pc_z[k] * osh_969[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, pa_z, pc_x, pc_y, pc_z, nsi0_1018, nsh_782, \
                         nsh_975, nsi1_1018, osg0_699, osg1_699, osh_971, \
                         osh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_19 * nsh_782[k]
                    + f_3 * pc_y[k] * osh_971[k];

        t_1297[k] = f_12 * nsh_975[k]
                    + f_6 * osg0_699[k]
                    - f_7 * osg1_699[k]
                    + f_3 * pc_x[k] * osh_975[k];

        t_1298[k] = pa_z[k] * nsi0_1018[k]
                    - f_10 * pc_z[k] * nsi1_1018[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, pa_z, pc_y, pc_z, nsi0_1020, nsh_762, \
                         nsh_763, nsh_786, nsi1_1020, osh_972, \
                         osh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = f_11 * nsh_762[k]
                    + f_3 * pc_z[k] * osh_972[k];

        t_1300[k] = pa_z[k] * nsi0_1020[k]
                    + f_12 * nsh_763[k]
                    - f_10 * pc_z[k] * nsi1_1020[k];

        t_1301[k] = f_19 * nsh_786[k]
                    + f_3 * pc_y[k] * osh_975[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, t_1305, pc_x, nsh_980, nsh_981, nsh_982, \
                         nsh_983, osg0_704, osg1_704, osh_980, osh_981, osh_982, \
                         osh_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_12 * nsh_980[k]
                    + f_4 * osg0_704[k]
                    - f_5 * osg1_704[k]
                    + f_3 * pc_x[k] * osh_980[k];

        t_1303[k] = f_12 * nsh_981[k]
                    + f_3 * pc_x[k] * osh_981[k];

        t_1304[k] = f_12 * nsh_982[k]
                    + f_3 * pc_x[k] * osh_982[k];

        t_1305[k] = f_12 * nsh_983[k]
                    + f_3 * pc_x[k] * osh_983[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, pa_z, pc_x, pc_z, nsi0_1029, nsh_984, \
                         nsh_985, nsh_986, nsi1_1029, osh_984, osh_985, \
                         osh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_12 * nsh_984[k]
                    + f_3 * pc_x[k] * osh_984[k];

        t_1307[k] = f_12 * nsh_985[k]
                    + f_3 * pc_x[k] * osh_985[k];

        t_1308[k] = f_12 * nsh_986[k]
                    + f_3 * pc_x[k] * osh_986[k];

        t_1309[k] = pa_z[k] * nsi0_1029[k]
                    - f_10 * pc_z[k] * nsi1_1029[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, pc_y, pc_z, nsh_771, nsh_794, nsh_795, \
                         osg0_702, osg0_703, osg1_702, osg1_703, osh_981, osh_983, \
                         osh_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_11 * nsh_771[k]
                    + f_3 * pc_z[k] * osh_981[k];

        t_1311[k] = f_19 * nsh_794[k]
                    + f_8 * osg0_702[k]
                    - f_9 * osg1_702[k]
                    + f_3 * pc_y[k] * osh_983[k];

        t_1312[k] = f_19 * nsh_795[k]
                    + f_6 * osg0_703[k]
                    - f_7 * osg1_703[k]
                    + f_3 * pc_y[k] * osh_984[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pc_y, pc_z, nsh_776, nsh_796, nsh_797, \
                         osg0_704, osg1_704, osh_985, osh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_19 * nsh_796[k]
                    + f_4 * osg0_704[k]
                    - f_5 * osg1_704[k]
                    + f_3 * pc_y[k] * osh_985[k];

        t_1314[k] = f_19 * nsh_797[k]
                    + f_3 * pc_y[k] * osh_986[k];

        t_1315[k] = f_11 * nsh_776[k]
                    + f_1 * osg0_704[k]
                    - f_2 * osg1_704[k]
                    + f_3 * pc_z[k] * osh_986[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, pc_x, pc_y, pc_z, nsh_777, nsh_798, nsh_987, \
                         osg0_705, osg1_705, osh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_12 * nsh_987[k]
                    + f_1 * osg0_705[k]
                    - f_2 * osg1_705[k]
                    + f_3 * pc_x[k] * osh_987[k];

        t_1317[k] = f_20 * nsh_798[k]
                    + f_3 * pc_y[k] * osh_987[k];

        t_1318[k] = f_12 * nsh_777[k]
                    + f_3 * pc_z[k] * osh_987[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pc_x, pc_y, nsh_800, nsh_990, nsh_992, \
                         osg0_708, osg0_710, osg1_708, osg1_710, osh_989, osh_990, \
                         osh_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_12 * nsh_990[k]
                    + f_8 * osg0_708[k]
                    - f_9 * osg1_708[k]
                    + f_3 * pc_x[k] * osh_990[k];

        t_1320[k] = f_20 * nsh_800[k]
                    + f_3 * pc_y[k] * osh_989[k];

        t_1321[k] = f_12 * nsh_992[k]
                    + f_8 * osg0_710[k]
                    - f_9 * osg1_710[k]
                    + f_3 * pc_x[k] * osh_992[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pc_x, pc_y, pc_z, nsh_780, nsh_803, nsh_993, \
                         osg0_711, osg1_711, osh_990, osh_992, \
                         osh_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_12 * nsh_993[k]
                    + f_6 * osg0_711[k]
                    - f_7 * osg1_711[k]
                    + f_3 * pc_x[k] * osh_993[k];

        t_1323[k] = f_12 * nsh_780[k]
                    + f_3 * pc_z[k] * osh_990[k];

        t_1324[k] = f_20 * nsh_803[k]
                    + f_3 * pc_y[k] * osh_992[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_x, pc_z, nsh_783, nsh_996, nsh_997, \
                         osg0_714, osg0_715, osg1_714, osg1_715, osh_993, osh_996, \
                         osh_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_12 * nsh_996[k]
                    + f_6 * osg0_714[k]
                    - f_7 * osg1_714[k]
                    + f_3 * pc_x[k] * osh_996[k];

        t_1326[k] = f_12 * nsh_997[k]
                    + f_4 * osg0_715[k]
                    - f_5 * osg1_715[k]
                    + f_3 * pc_x[k] * osh_997[k];

        t_1327[k] = f_12 * nsh_783[k]
                    + f_3 * pc_z[k] * osh_993[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pc_x, pc_y, nsh_807, nsh_999, nsh_1001, \
                         osg0_717, osg0_719, osg1_717, osg1_719, osh_996, osh_999, \
                         osh_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_12 * nsh_999[k]
                    + f_4 * osg0_717[k]
                    - f_5 * osg1_717[k]
                    + f_3 * pc_x[k] * osh_999[k];

        t_1329[k] = f_20 * nsh_807[k]
                    + f_3 * pc_y[k] * osh_996[k];

        t_1330[k] = f_12 * nsh_1001[k]
                    + f_4 * osg0_719[k]
                    - f_5 * osg1_719[k]
                    + f_3 * pc_x[k] * osh_1001[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, t_1335, pc_x, nsh_1002, nsh_1003, \
                         nsh_1004, nsh_1005, nsh_1006, osh_1002, osh_1003, osh_1004, osh_1005, \
                         osh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_12 * nsh_1002[k]
                    + f_3 * pc_x[k] * osh_1002[k];

        t_1332[k] = f_12 * nsh_1003[k]
                    + f_3 * pc_x[k] * osh_1003[k];

        t_1333[k] = f_12 * nsh_1004[k]
                    + f_3 * pc_x[k] * osh_1004[k];

        t_1334[k] = f_12 * nsh_1005[k]
                    + f_3 * pc_x[k] * osh_1005[k];

        t_1335[k] = f_12 * nsh_1006[k]
                    + f_3 * pc_x[k] * osh_1006[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, pc_x, pc_y, pc_z, nsh_792, nsh_813, nsh_1007, \
                         osg0_715, osg1_715, osh_1002, osh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_12 * nsh_1007[k]
                    + f_3 * pc_x[k] * osh_1007[k];

        t_1337[k] = f_20 * nsh_813[k]
                    + f_1 * osg0_715[k]
                    - f_2 * osg1_715[k]
                    + f_3 * pc_y[k] * osh_1002[k];

        t_1338[k] = f_12 * nsh_792[k]
                    + f_3 * pc_z[k] * osh_1002[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, pc_y, nsh_815, nsh_816, nsh_817, osg0_717, \
                         osg0_718, osg0_719, osg1_717, osg1_718, osg1_719, osh_1004, osh_1005, \
                         osh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_20 * nsh_815[k]
                    + f_8 * osg0_717[k]
                    - f_9 * osg1_717[k]
                    + f_3 * pc_y[k] * osh_1004[k];

        t_1340[k] = f_20 * nsh_816[k]
                    + f_6 * osg0_718[k]
                    - f_7 * osg1_718[k]
                    + f_3 * pc_y[k] * osh_1005[k];

        t_1341[k] = f_20 * nsh_817[k]
                    + f_4 * osg0_719[k]
                    - f_5 * osg1_719[k]
                    + f_3 * pc_y[k] * osh_1006[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pc_x, pc_y, pc_z, nsh_797, nsh_818, nsh_1008, \
                         osg0_719, osg0_720, osg1_719, osg1_720, osh_1007, \
                         osh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_20 * nsh_818[k]
                    + f_3 * pc_y[k] * osh_1007[k];

        t_1343[k] = f_12 * nsh_797[k]
                    + f_1 * osg0_719[k]
                    - f_2 * osg1_719[k]
                    + f_3 * pc_z[k] * osh_1007[k];

        t_1344[k] = f_12 * nsh_1008[k]
                    + f_1 * osg0_720[k]
                    - f_2 * osg1_720[k]
                    + f_3 * pc_x[k] * osh_1008[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, t_1348, pc_x, pc_y, pc_z, nsh_798, nsh_819, \
                         nsh_821, nsh_1011, osg0_723, osg1_723, osh_1008, osh_1010, \
                         osh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_21 * nsh_819[k]
                    + f_3 * pc_y[k] * osh_1008[k];

        t_1346[k] = f_13 * nsh_798[k]
                    + f_3 * pc_z[k] * osh_1008[k];

        t_1347[k] = f_12 * nsh_1011[k]
                    + f_8 * osg0_723[k]
                    - f_9 * osg1_723[k]
                    + f_3 * pc_x[k] * osh_1011[k];

        t_1348[k] = f_21 * nsh_821[k]
                    + f_3 * pc_y[k] * osh_1010[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pc_x, pc_z, nsh_801, nsh_1013, nsh_1014, \
                         osg0_725, osg0_726, osg1_725, osg1_726, osh_1011, osh_1013, \
                         osh_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_12 * nsh_1013[k]
                    + f_8 * osg0_725[k]
                    - f_9 * osg1_725[k]
                    + f_3 * pc_x[k] * osh_1013[k];

        t_1350[k] = f_12 * nsh_1014[k]
                    + f_6 * osg0_726[k]
                    - f_7 * osg1_726[k]
                    + f_3 * pc_x[k] * osh_1014[k];

        t_1351[k] = f_13 * nsh_801[k]
                    + f_3 * pc_z[k] * osh_1011[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pc_x, pc_y, nsh_824, nsh_1017, nsh_1018, \
                         osg0_729, osg0_730, osg1_729, osg1_730, osh_1013, osh_1017, \
                         osh_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_21 * nsh_824[k]
                    + f_3 * pc_y[k] * osh_1013[k];

        t_1353[k] = f_12 * nsh_1017[k]
                    + f_6 * osg0_729[k]
                    - f_7 * osg1_729[k]
                    + f_3 * pc_x[k] * osh_1017[k];

        t_1354[k] = f_12 * nsh_1018[k]
                    + f_4 * osg0_730[k]
                    - f_5 * osg1_730[k]
                    + f_3 * pc_x[k] * osh_1018[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pc_x, pc_y, pc_z, nsh_804, nsh_828, nsh_1020, \
                         osg0_732, osg1_732, osh_1014, osh_1017, \
                         osh_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_13 * nsh_804[k]
                    + f_3 * pc_z[k] * osh_1014[k];

        t_1356[k] = f_12 * nsh_1020[k]
                    + f_4 * osg0_732[k]
                    - f_5 * osg1_732[k]
                    + f_3 * pc_x[k] * osh_1020[k];

        t_1357[k] = f_21 * nsh_828[k]
                    + f_3 * pc_y[k] * osh_1017[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pc_x, nsh_1022, nsh_1023, nsh_1024, \
                         nsh_1025, osg0_734, osg1_734, osh_1022, osh_1023, osh_1024, \
                         osh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_12 * nsh_1022[k]
                    + f_4 * osg0_734[k]
                    - f_5 * osg1_734[k]
                    + f_3 * pc_x[k] * osh_1022[k];

        t_1359[k] = f_12 * nsh_1023[k]
                    + f_3 * pc_x[k] * osh_1023[k];

        t_1360[k] = f_12 * nsh_1024[k]
                    + f_3 * pc_x[k] * osh_1024[k];

        t_1361[k] = f_12 * nsh_1025[k]
                    + f_3 * pc_x[k] * osh_1025[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, pc_x, pc_y, nsh_834, nsh_1026, \
                         nsh_1027, nsh_1028, osg0_730, osg1_730, osh_1023, osh_1026, osh_1027, \
                         osh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_12 * nsh_1026[k]
                    + f_3 * pc_x[k] * osh_1026[k];

        t_1363[k] = f_12 * nsh_1027[k]
                    + f_3 * pc_x[k] * osh_1027[k];

        t_1364[k] = f_12 * nsh_1028[k]
                    + f_3 * pc_x[k] * osh_1028[k];

        t_1365[k] = f_21 * nsh_834[k]
                    + f_1 * osg0_730[k]
                    - f_2 * osg1_730[k]
                    + f_3 * pc_y[k] * osh_1023[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, pc_y, pc_z, nsh_813, nsh_836, nsh_837, \
                         osg0_732, osg0_733, osg1_732, osg1_733, osh_1023, osh_1025, \
                         osh_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_13 * nsh_813[k]
                    + f_3 * pc_z[k] * osh_1023[k];

        t_1367[k] = f_21 * nsh_836[k]
                    + f_8 * osg0_732[k]
                    - f_9 * osg1_732[k]
                    + f_3 * pc_y[k] * osh_1025[k];

        t_1368[k] = f_21 * nsh_837[k]
                    + f_6 * osg0_733[k]
                    - f_7 * osg1_733[k]
                    + f_3 * pc_y[k] * osh_1026[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, pc_y, pc_z, nsh_818, nsh_838, nsh_839, \
                         osg0_734, osg1_734, osh_1027, osh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_21 * nsh_838[k]
                    + f_4 * osg0_734[k]
                    - f_5 * osg1_734[k]
                    + f_3 * pc_y[k] * osh_1027[k];

        t_1370[k] = f_21 * nsh_839[k]
                    + f_3 * pc_y[k] * osh_1028[k];

        t_1371[k] = f_13 * nsh_818[k]
                    + f_1 * osg0_734[k]
                    - f_2 * osg1_734[k]
                    + f_3 * pc_z[k] * osh_1028[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t nsh, const size_t osg0,
                                                           const size_t osg1, const size_t osh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);
    auto *t_1388 = buffer.data(target + 1388);
    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);
    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh_819 = buffer.data(nsh + 819);
    const auto *nsh_822 = buffer.data(nsh + 822);
    const auto *nsh_825 = buffer.data(nsh + 825);
    const auto *nsh_834 = buffer.data(nsh + 834);
    const auto *nsh_839 = buffer.data(nsh + 839);
    const auto *nsh_840 = buffer.data(nsh + 840);
    const auto *nsh_842 = buffer.data(nsh + 842);
    const auto *nsh_843 = buffer.data(nsh + 843);
    const auto *nsh_845 = buffer.data(nsh + 845);
    const auto *nsh_846 = buffer.data(nsh + 846);
    const auto *nsh_849 = buffer.data(nsh + 849);
    const auto *nsh_855 = buffer.data(nsh + 855);
    const auto *nsh_857 = buffer.data(nsh + 857);
    const auto *nsh_858 = buffer.data(nsh + 858);
    const auto *nsh_859 = buffer.data(nsh + 859);
    const auto *nsh_860 = buffer.data(nsh + 860);
    const auto *nsh_861 = buffer.data(nsh + 861);
    const auto *nsh_863 = buffer.data(nsh + 863);
    const auto *nsh_864 = buffer.data(nsh + 864);
    const auto *nsh_866 = buffer.data(nsh + 866);
    const auto *nsh_867 = buffer.data(nsh + 867);
    const auto *nsh_870 = buffer.data(nsh + 870);
    const auto *nsh_876 = buffer.data(nsh + 876);
    const auto *nsh_878 = buffer.data(nsh + 878);
    const auto *nsh_879 = buffer.data(nsh + 879);
    const auto *nsh_880 = buffer.data(nsh + 880);
    const auto *nsh_881 = buffer.data(nsh + 881);
    const auto *nsh_882 = buffer.data(nsh + 882);
    const auto *nsh_884 = buffer.data(nsh + 884);
    const auto *nsh_885 = buffer.data(nsh + 885);
    const auto *nsh_887 = buffer.data(nsh + 887);
    const auto *nsh_888 = buffer.data(nsh + 888);
    const auto *nsh_891 = buffer.data(nsh + 891);
    const auto *nsh_897 = buffer.data(nsh + 897);
    const auto *nsh_899 = buffer.data(nsh + 899);
    const auto *nsh_900 = buffer.data(nsh + 900);
    const auto *nsh_901 = buffer.data(nsh + 901);
    const auto *nsh_902 = buffer.data(nsh + 902);
    const auto *nsh_903 = buffer.data(nsh + 903);
    const auto *nsh_905 = buffer.data(nsh + 905);
    const auto *nsh_908 = buffer.data(nsh + 908);
    const auto *nsh_912 = buffer.data(nsh + 912);
    const auto *nsh_918 = buffer.data(nsh + 918);
    const auto *nsh_1029 = buffer.data(nsh + 1029);
    const auto *nsh_1032 = buffer.data(nsh + 1032);
    const auto *nsh_1034 = buffer.data(nsh + 1034);
    const auto *nsh_1035 = buffer.data(nsh + 1035);
    const auto *nsh_1038 = buffer.data(nsh + 1038);
    const auto *nsh_1039 = buffer.data(nsh + 1039);
    const auto *nsh_1041 = buffer.data(nsh + 1041);
    const auto *nsh_1043 = buffer.data(nsh + 1043);
    const auto *nsh_1044 = buffer.data(nsh + 1044);
    const auto *nsh_1045 = buffer.data(nsh + 1045);
    const auto *nsh_1046 = buffer.data(nsh + 1046);
    const auto *nsh_1047 = buffer.data(nsh + 1047);
    const auto *nsh_1048 = buffer.data(nsh + 1048);
    const auto *nsh_1049 = buffer.data(nsh + 1049);
    const auto *nsh_1050 = buffer.data(nsh + 1050);
    const auto *nsh_1053 = buffer.data(nsh + 1053);
    const auto *nsh_1055 = buffer.data(nsh + 1055);
    const auto *nsh_1056 = buffer.data(nsh + 1056);
    const auto *nsh_1059 = buffer.data(nsh + 1059);
    const auto *nsh_1060 = buffer.data(nsh + 1060);
    const auto *nsh_1062 = buffer.data(nsh + 1062);
    const auto *nsh_1064 = buffer.data(nsh + 1064);
    const auto *nsh_1065 = buffer.data(nsh + 1065);
    const auto *nsh_1066 = buffer.data(nsh + 1066);
    const auto *nsh_1067 = buffer.data(nsh + 1067);
    const auto *nsh_1068 = buffer.data(nsh + 1068);
    const auto *nsh_1069 = buffer.data(nsh + 1069);
    const auto *nsh_1070 = buffer.data(nsh + 1070);
    const auto *nsh_1071 = buffer.data(nsh + 1071);
    const auto *nsh_1074 = buffer.data(nsh + 1074);
    const auto *nsh_1076 = buffer.data(nsh + 1076);
    const auto *nsh_1077 = buffer.data(nsh + 1077);
    const auto *nsh_1080 = buffer.data(nsh + 1080);
    const auto *nsh_1081 = buffer.data(nsh + 1081);
    const auto *nsh_1083 = buffer.data(nsh + 1083);
    const auto *nsh_1085 = buffer.data(nsh + 1085);
    const auto *nsh_1086 = buffer.data(nsh + 1086);
    const auto *nsh_1087 = buffer.data(nsh + 1087);
    const auto *nsh_1088 = buffer.data(nsh + 1088);
    const auto *nsh_1089 = buffer.data(nsh + 1089);
    const auto *nsh_1090 = buffer.data(nsh + 1090);
    const auto *nsh_1091 = buffer.data(nsh + 1091);
    const auto *nsh_1092 = buffer.data(nsh + 1092);
    const auto *nsh_1095 = buffer.data(nsh + 1095);
    const auto *nsh_1097 = buffer.data(nsh + 1097);
    const auto *nsh_1098 = buffer.data(nsh + 1098);
    const auto *nsh_1101 = buffer.data(nsh + 1101);
    const auto *nsh_1102 = buffer.data(nsh + 1102);
    const auto *nsh_1104 = buffer.data(nsh + 1104);
    const auto *nsh_1106 = buffer.data(nsh + 1106);
    const auto *nsh_1107 = buffer.data(nsh + 1107);
    const auto *nsh_1108 = buffer.data(nsh + 1108);
    const auto *nsh_1109 = buffer.data(nsh + 1109);
    const auto *nsh_1110 = buffer.data(nsh + 1110);
    const auto *nsh_1111 = buffer.data(nsh + 1111);
    const auto *nsh_1112 = buffer.data(nsh + 1112);

    const auto *osg0_735 = buffer.data(osg0 + 735);
    const auto *osg0_738 = buffer.data(osg0 + 738);
    const auto *osg0_740 = buffer.data(osg0 + 740);
    const auto *osg0_741 = buffer.data(osg0 + 741);
    const auto *osg0_744 = buffer.data(osg0 + 744);
    const auto *osg0_745 = buffer.data(osg0 + 745);
    const auto *osg0_747 = buffer.data(osg0 + 747);
    const auto *osg0_748 = buffer.data(osg0 + 748);
    const auto *osg0_749 = buffer.data(osg0 + 749);
    const auto *osg0_750 = buffer.data(osg0 + 750);
    const auto *osg0_753 = buffer.data(osg0 + 753);
    const auto *osg0_755 = buffer.data(osg0 + 755);
    const auto *osg0_756 = buffer.data(osg0 + 756);
    const auto *osg0_759 = buffer.data(osg0 + 759);
    const auto *osg0_760 = buffer.data(osg0 + 760);
    const auto *osg0_762 = buffer.data(osg0 + 762);
    const auto *osg0_763 = buffer.data(osg0 + 763);
    const auto *osg0_764 = buffer.data(osg0 + 764);
    const auto *osg0_765 = buffer.data(osg0 + 765);
    const auto *osg0_768 = buffer.data(osg0 + 768);
    const auto *osg0_770 = buffer.data(osg0 + 770);
    const auto *osg0_771 = buffer.data(osg0 + 771);
    const auto *osg0_774 = buffer.data(osg0 + 774);
    const auto *osg0_775 = buffer.data(osg0 + 775);
    const auto *osg0_777 = buffer.data(osg0 + 777);
    const auto *osg0_778 = buffer.data(osg0 + 778);
    const auto *osg0_779 = buffer.data(osg0 + 779);
    const auto *osg0_780 = buffer.data(osg0 + 780);
    const auto *osg0_783 = buffer.data(osg0 + 783);
    const auto *osg0_785 = buffer.data(osg0 + 785);
    const auto *osg0_786 = buffer.data(osg0 + 786);
    const auto *osg0_789 = buffer.data(osg0 + 789);
    const auto *osg0_790 = buffer.data(osg0 + 790);
    const auto *osg0_792 = buffer.data(osg0 + 792);
    const auto *osg0_794 = buffer.data(osg0 + 794);

    const auto *osg1_735 = buffer.data(osg1 + 735);
    const auto *osg1_738 = buffer.data(osg1 + 738);
    const auto *osg1_740 = buffer.data(osg1 + 740);
    const auto *osg1_741 = buffer.data(osg1 + 741);
    const auto *osg1_744 = buffer.data(osg1 + 744);
    const auto *osg1_745 = buffer.data(osg1 + 745);
    const auto *osg1_747 = buffer.data(osg1 + 747);
    const auto *osg1_748 = buffer.data(osg1 + 748);
    const auto *osg1_749 = buffer.data(osg1 + 749);
    const auto *osg1_750 = buffer.data(osg1 + 750);
    const auto *osg1_753 = buffer.data(osg1 + 753);
    const auto *osg1_755 = buffer.data(osg1 + 755);
    const auto *osg1_756 = buffer.data(osg1 + 756);
    const auto *osg1_759 = buffer.data(osg1 + 759);
    const auto *osg1_760 = buffer.data(osg1 + 760);
    const auto *osg1_762 = buffer.data(osg1 + 762);
    const auto *osg1_763 = buffer.data(osg1 + 763);
    const auto *osg1_764 = buffer.data(osg1 + 764);
    const auto *osg1_765 = buffer.data(osg1 + 765);
    const auto *osg1_768 = buffer.data(osg1 + 768);
    const auto *osg1_770 = buffer.data(osg1 + 770);
    const auto *osg1_771 = buffer.data(osg1 + 771);
    const auto *osg1_774 = buffer.data(osg1 + 774);
    const auto *osg1_775 = buffer.data(osg1 + 775);
    const auto *osg1_777 = buffer.data(osg1 + 777);
    const auto *osg1_778 = buffer.data(osg1 + 778);
    const auto *osg1_779 = buffer.data(osg1 + 779);
    const auto *osg1_780 = buffer.data(osg1 + 780);
    const auto *osg1_783 = buffer.data(osg1 + 783);
    const auto *osg1_785 = buffer.data(osg1 + 785);
    const auto *osg1_786 = buffer.data(osg1 + 786);
    const auto *osg1_789 = buffer.data(osg1 + 789);
    const auto *osg1_790 = buffer.data(osg1 + 790);
    const auto *osg1_792 = buffer.data(osg1 + 792);
    const auto *osg1_794 = buffer.data(osg1 + 794);

    const auto *osh_1029 = buffer.data(osh + 1029);
    const auto *osh_1031 = buffer.data(osh + 1031);
    const auto *osh_1032 = buffer.data(osh + 1032);
    const auto *osh_1034 = buffer.data(osh + 1034);
    const auto *osh_1035 = buffer.data(osh + 1035);
    const auto *osh_1038 = buffer.data(osh + 1038);
    const auto *osh_1039 = buffer.data(osh + 1039);
    const auto *osh_1041 = buffer.data(osh + 1041);
    const auto *osh_1043 = buffer.data(osh + 1043);
    const auto *osh_1044 = buffer.data(osh + 1044);
    const auto *osh_1045 = buffer.data(osh + 1045);
    const auto *osh_1046 = buffer.data(osh + 1046);
    const auto *osh_1047 = buffer.data(osh + 1047);
    const auto *osh_1048 = buffer.data(osh + 1048);
    const auto *osh_1049 = buffer.data(osh + 1049);
    const auto *osh_1050 = buffer.data(osh + 1050);
    const auto *osh_1052 = buffer.data(osh + 1052);
    const auto *osh_1053 = buffer.data(osh + 1053);
    const auto *osh_1055 = buffer.data(osh + 1055);
    const auto *osh_1056 = buffer.data(osh + 1056);
    const auto *osh_1059 = buffer.data(osh + 1059);
    const auto *osh_1060 = buffer.data(osh + 1060);
    const auto *osh_1062 = buffer.data(osh + 1062);
    const auto *osh_1064 = buffer.data(osh + 1064);
    const auto *osh_1065 = buffer.data(osh + 1065);
    const auto *osh_1066 = buffer.data(osh + 1066);
    const auto *osh_1067 = buffer.data(osh + 1067);
    const auto *osh_1068 = buffer.data(osh + 1068);
    const auto *osh_1069 = buffer.data(osh + 1069);
    const auto *osh_1070 = buffer.data(osh + 1070);
    const auto *osh_1071 = buffer.data(osh + 1071);
    const auto *osh_1073 = buffer.data(osh + 1073);
    const auto *osh_1074 = buffer.data(osh + 1074);
    const auto *osh_1076 = buffer.data(osh + 1076);
    const auto *osh_1077 = buffer.data(osh + 1077);
    const auto *osh_1080 = buffer.data(osh + 1080);
    const auto *osh_1081 = buffer.data(osh + 1081);
    const auto *osh_1083 = buffer.data(osh + 1083);
    const auto *osh_1085 = buffer.data(osh + 1085);
    const auto *osh_1086 = buffer.data(osh + 1086);
    const auto *osh_1087 = buffer.data(osh + 1087);
    const auto *osh_1088 = buffer.data(osh + 1088);
    const auto *osh_1089 = buffer.data(osh + 1089);
    const auto *osh_1090 = buffer.data(osh + 1090);
    const auto *osh_1091 = buffer.data(osh + 1091);
    const auto *osh_1092 = buffer.data(osh + 1092);
    const auto *osh_1094 = buffer.data(osh + 1094);
    const auto *osh_1095 = buffer.data(osh + 1095);
    const auto *osh_1097 = buffer.data(osh + 1097);
    const auto *osh_1098 = buffer.data(osh + 1098);
    const auto *osh_1101 = buffer.data(osh + 1101);
    const auto *osh_1102 = buffer.data(osh + 1102);
    const auto *osh_1104 = buffer.data(osh + 1104);
    const auto *osh_1106 = buffer.data(osh + 1106);
    const auto *osh_1107 = buffer.data(osh + 1107);
    const auto *osh_1108 = buffer.data(osh + 1108);
    const auto *osh_1109 = buffer.data(osh + 1109);
    const auto *osh_1110 = buffer.data(osh + 1110);
    const auto *osh_1111 = buffer.data(osh + 1111);
    const auto *osh_1112 = buffer.data(osh + 1112);

#pragma omp simd aligned(t_1372, t_1373, t_1374, pc_x, pc_y, pc_z, nsh_819, nsh_840, nsh_1029, \
                         osg0_735, osg1_735, osh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1372[k] = f_12 * nsh_1029[k]
                    + f_1 * osg0_735[k]
                    - f_2 * osg1_735[k]
                    + f_3 * pc_x[k] * osh_1029[k];

        t_1373[k] = f_22 * nsh_840[k]
                    + f_3 * pc_y[k] * osh_1029[k];

        t_1374[k] = f_14 * nsh_819[k]
                    + f_3 * pc_z[k] * osh_1029[k];
    }

#pragma omp simd aligned(t_1375, t_1376, t_1377, pc_x, pc_y, nsh_842, nsh_1032, nsh_1034, \
                         osg0_738, osg0_740, osg1_738, osg1_740, osh_1031, osh_1032, \
                         osh_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1375[k] = f_12 * nsh_1032[k]
                    + f_8 * osg0_738[k]
                    - f_9 * osg1_738[k]
                    + f_3 * pc_x[k] * osh_1032[k];

        t_1376[k] = f_22 * nsh_842[k]
                    + f_3 * pc_y[k] * osh_1031[k];

        t_1377[k] = f_12 * nsh_1034[k]
                    + f_8 * osg0_740[k]
                    - f_9 * osg1_740[k]
                    + f_3 * pc_x[k] * osh_1034[k];
    }

#pragma omp simd aligned(t_1378, t_1379, t_1380, pc_x, pc_y, pc_z, nsh_822, nsh_845, nsh_1035, \
                         osg0_741, osg1_741, osh_1032, osh_1034, \
                         osh_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1378[k] = f_12 * nsh_1035[k]
                    + f_6 * osg0_741[k]
                    - f_7 * osg1_741[k]
                    + f_3 * pc_x[k] * osh_1035[k];

        t_1379[k] = f_14 * nsh_822[k]
                    + f_3 * pc_z[k] * osh_1032[k];

        t_1380[k] = f_22 * nsh_845[k]
                    + f_3 * pc_y[k] * osh_1034[k];
    }

#pragma omp simd aligned(t_1381, t_1382, t_1383, pc_x, pc_z, nsh_825, nsh_1038, nsh_1039, \
                         osg0_744, osg0_745, osg1_744, osg1_745, osh_1035, osh_1038, \
                         osh_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1381[k] = f_12 * nsh_1038[k]
                    + f_6 * osg0_744[k]
                    - f_7 * osg1_744[k]
                    + f_3 * pc_x[k] * osh_1038[k];

        t_1382[k] = f_12 * nsh_1039[k]
                    + f_4 * osg0_745[k]
                    - f_5 * osg1_745[k]
                    + f_3 * pc_x[k] * osh_1039[k];

        t_1383[k] = f_14 * nsh_825[k]
                    + f_3 * pc_z[k] * osh_1035[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, pc_x, pc_y, nsh_849, nsh_1041, nsh_1043, \
                         osg0_747, osg0_749, osg1_747, osg1_749, osh_1038, osh_1041, \
                         osh_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_12 * nsh_1041[k]
                    + f_4 * osg0_747[k]
                    - f_5 * osg1_747[k]
                    + f_3 * pc_x[k] * osh_1041[k];

        t_1385[k] = f_22 * nsh_849[k]
                    + f_3 * pc_y[k] * osh_1038[k];

        t_1386[k] = f_12 * nsh_1043[k]
                    + f_4 * osg0_749[k]
                    - f_5 * osg1_749[k]
                    + f_3 * pc_x[k] * osh_1043[k];
    }

#pragma omp simd aligned(t_1387, t_1388, t_1389, t_1390, t_1391, pc_x, nsh_1044, nsh_1045, \
                         nsh_1046, nsh_1047, nsh_1048, osh_1044, osh_1045, osh_1046, osh_1047, \
                         osh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1387[k] = f_12 * nsh_1044[k]
                    + f_3 * pc_x[k] * osh_1044[k];

        t_1388[k] = f_12 * nsh_1045[k]
                    + f_3 * pc_x[k] * osh_1045[k];

        t_1389[k] = f_12 * nsh_1046[k]
                    + f_3 * pc_x[k] * osh_1046[k];

        t_1390[k] = f_12 * nsh_1047[k]
                    + f_3 * pc_x[k] * osh_1047[k];

        t_1391[k] = f_12 * nsh_1048[k]
                    + f_3 * pc_x[k] * osh_1048[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, pc_x, pc_y, pc_z, nsh_834, nsh_855, nsh_1049, \
                         osg0_745, osg1_745, osh_1044, osh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_12 * nsh_1049[k]
                    + f_3 * pc_x[k] * osh_1049[k];

        t_1393[k] = f_22 * nsh_855[k]
                    + f_1 * osg0_745[k]
                    - f_2 * osg1_745[k]
                    + f_3 * pc_y[k] * osh_1044[k];

        t_1394[k] = f_14 * nsh_834[k]
                    + f_3 * pc_z[k] * osh_1044[k];
    }

#pragma omp simd aligned(t_1395, t_1396, t_1397, pc_y, nsh_857, nsh_858, nsh_859, osg0_747, \
                         osg0_748, osg0_749, osg1_747, osg1_748, osg1_749, osh_1046, osh_1047, \
                         osh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1395[k] = f_22 * nsh_857[k]
                    + f_8 * osg0_747[k]
                    - f_9 * osg1_747[k]
                    + f_3 * pc_y[k] * osh_1046[k];

        t_1396[k] = f_22 * nsh_858[k]
                    + f_6 * osg0_748[k]
                    - f_7 * osg1_748[k]
                    + f_3 * pc_y[k] * osh_1047[k];

        t_1397[k] = f_22 * nsh_859[k]
                    + f_4 * osg0_749[k]
                    - f_5 * osg1_749[k]
                    + f_3 * pc_y[k] * osh_1048[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, pc_x, pc_y, pc_z, nsh_839, nsh_860, nsh_1050, \
                         osg0_749, osg0_750, osg1_749, osg1_750, osh_1049, \
                         osh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_22 * nsh_860[k]
                    + f_3 * pc_y[k] * osh_1049[k];

        t_1399[k] = f_14 * nsh_839[k]
                    + f_1 * osg0_749[k]
                    - f_2 * osg1_749[k]
                    + f_3 * pc_z[k] * osh_1049[k];

        t_1400[k] = f_12 * nsh_1050[k]
                    + f_1 * osg0_750[k]
                    - f_2 * osg1_750[k]
                    + f_3 * pc_x[k] * osh_1050[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, t_1404, pc_x, pc_y, pc_z, nsh_840, nsh_861, \
                         nsh_863, nsh_1053, osg0_753, osg1_753, osh_1050, osh_1052, \
                         osh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = f_14 * nsh_861[k]
                    + f_3 * pc_y[k] * osh_1050[k];

        t_1402[k] = f_22 * nsh_840[k]
                    + f_3 * pc_z[k] * osh_1050[k];

        t_1403[k] = f_12 * nsh_1053[k]
                    + f_8 * osg0_753[k]
                    - f_9 * osg1_753[k]
                    + f_3 * pc_x[k] * osh_1053[k];

        t_1404[k] = f_14 * nsh_863[k]
                    + f_3 * pc_y[k] * osh_1052[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, pc_x, pc_z, nsh_843, nsh_1055, nsh_1056, \
                         osg0_755, osg0_756, osg1_755, osg1_756, osh_1053, osh_1055, \
                         osh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_12 * nsh_1055[k]
                    + f_8 * osg0_755[k]
                    - f_9 * osg1_755[k]
                    + f_3 * pc_x[k] * osh_1055[k];

        t_1406[k] = f_12 * nsh_1056[k]
                    + f_6 * osg0_756[k]
                    - f_7 * osg1_756[k]
                    + f_3 * pc_x[k] * osh_1056[k];

        t_1407[k] = f_22 * nsh_843[k]
                    + f_3 * pc_z[k] * osh_1053[k];
    }

#pragma omp simd aligned(t_1408, t_1409, t_1410, pc_x, pc_y, nsh_866, nsh_1059, nsh_1060, \
                         osg0_759, osg0_760, osg1_759, osg1_760, osh_1055, osh_1059, \
                         osh_1060 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1408[k] = f_14 * nsh_866[k]
                    + f_3 * pc_y[k] * osh_1055[k];

        t_1409[k] = f_12 * nsh_1059[k]
                    + f_6 * osg0_759[k]
                    - f_7 * osg1_759[k]
                    + f_3 * pc_x[k] * osh_1059[k];

        t_1410[k] = f_12 * nsh_1060[k]
                    + f_4 * osg0_760[k]
                    - f_5 * osg1_760[k]
                    + f_3 * pc_x[k] * osh_1060[k];
    }

#pragma omp simd aligned(t_1411, t_1412, t_1413, pc_x, pc_y, pc_z, nsh_846, nsh_870, nsh_1062, \
                         osg0_762, osg1_762, osh_1056, osh_1059, \
                         osh_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1411[k] = f_22 * nsh_846[k]
                    + f_3 * pc_z[k] * osh_1056[k];

        t_1412[k] = f_12 * nsh_1062[k]
                    + f_4 * osg0_762[k]
                    - f_5 * osg1_762[k]
                    + f_3 * pc_x[k] * osh_1062[k];

        t_1413[k] = f_14 * nsh_870[k]
                    + f_3 * pc_y[k] * osh_1059[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, t_1417, pc_x, nsh_1064, nsh_1065, nsh_1066, \
                         nsh_1067, osg0_764, osg1_764, osh_1064, osh_1065, osh_1066, \
                         osh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = f_12 * nsh_1064[k]
                    + f_4 * osg0_764[k]
                    - f_5 * osg1_764[k]
                    + f_3 * pc_x[k] * osh_1064[k];

        t_1415[k] = f_12 * nsh_1065[k]
                    + f_3 * pc_x[k] * osh_1065[k];

        t_1416[k] = f_12 * nsh_1066[k]
                    + f_3 * pc_x[k] * osh_1066[k];

        t_1417[k] = f_12 * nsh_1067[k]
                    + f_3 * pc_x[k] * osh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, pc_x, pc_y, nsh_876, nsh_1068, \
                         nsh_1069, nsh_1070, osg0_760, osg1_760, osh_1065, osh_1068, osh_1069, \
                         osh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_12 * nsh_1068[k]
                    + f_3 * pc_x[k] * osh_1068[k];

        t_1419[k] = f_12 * nsh_1069[k]
                    + f_3 * pc_x[k] * osh_1069[k];

        t_1420[k] = f_12 * nsh_1070[k]
                    + f_3 * pc_x[k] * osh_1070[k];

        t_1421[k] = f_14 * nsh_876[k]
                    + f_1 * osg0_760[k]
                    - f_2 * osg1_760[k]
                    + f_3 * pc_y[k] * osh_1065[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, pc_y, pc_z, nsh_855, nsh_878, nsh_879, \
                         osg0_762, osg0_763, osg1_762, osg1_763, osh_1065, osh_1067, \
                         osh_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_22 * nsh_855[k]
                    + f_3 * pc_z[k] * osh_1065[k];

        t_1423[k] = f_14 * nsh_878[k]
                    + f_8 * osg0_762[k]
                    - f_9 * osg1_762[k]
                    + f_3 * pc_y[k] * osh_1067[k];

        t_1424[k] = f_14 * nsh_879[k]
                    + f_6 * osg0_763[k]
                    - f_7 * osg1_763[k]
                    + f_3 * pc_y[k] * osh_1068[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, pc_y, pc_z, nsh_860, nsh_880, nsh_881, \
                         osg0_764, osg1_764, osh_1069, osh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_14 * nsh_880[k]
                    + f_4 * osg0_764[k]
                    - f_5 * osg1_764[k]
                    + f_3 * pc_y[k] * osh_1069[k];

        t_1426[k] = f_14 * nsh_881[k]
                    + f_3 * pc_y[k] * osh_1070[k];

        t_1427[k] = f_22 * nsh_860[k]
                    + f_1 * osg0_764[k]
                    - f_2 * osg1_764[k]
                    + f_3 * pc_z[k] * osh_1070[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, pc_x, pc_y, pc_z, nsh_861, nsh_882, nsh_1071, \
                         osg0_765, osg1_765, osh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_12 * nsh_1071[k]
                    + f_1 * osg0_765[k]
                    - f_2 * osg1_765[k]
                    + f_3 * pc_x[k] * osh_1071[k];

        t_1429[k] = f_13 * nsh_882[k]
                    + f_3 * pc_y[k] * osh_1071[k];

        t_1430[k] = f_21 * nsh_861[k]
                    + f_3 * pc_z[k] * osh_1071[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, pc_x, pc_y, nsh_884, nsh_1074, nsh_1076, \
                         osg0_768, osg0_770, osg1_768, osg1_770, osh_1073, osh_1074, \
                         osh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = f_12 * nsh_1074[k]
                    + f_8 * osg0_768[k]
                    - f_9 * osg1_768[k]
                    + f_3 * pc_x[k] * osh_1074[k];

        t_1432[k] = f_13 * nsh_884[k]
                    + f_3 * pc_y[k] * osh_1073[k];

        t_1433[k] = f_12 * nsh_1076[k]
                    + f_8 * osg0_770[k]
                    - f_9 * osg1_770[k]
                    + f_3 * pc_x[k] * osh_1076[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, pc_x, pc_y, pc_z, nsh_864, nsh_887, nsh_1077, \
                         osg0_771, osg1_771, osh_1074, osh_1076, \
                         osh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = f_12 * nsh_1077[k]
                    + f_6 * osg0_771[k]
                    - f_7 * osg1_771[k]
                    + f_3 * pc_x[k] * osh_1077[k];

        t_1435[k] = f_21 * nsh_864[k]
                    + f_3 * pc_z[k] * osh_1074[k];

        t_1436[k] = f_13 * nsh_887[k]
                    + f_3 * pc_y[k] * osh_1076[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, pc_x, pc_z, nsh_867, nsh_1080, nsh_1081, \
                         osg0_774, osg0_775, osg1_774, osg1_775, osh_1077, osh_1080, \
                         osh_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = f_12 * nsh_1080[k]
                    + f_6 * osg0_774[k]
                    - f_7 * osg1_774[k]
                    + f_3 * pc_x[k] * osh_1080[k];

        t_1438[k] = f_12 * nsh_1081[k]
                    + f_4 * osg0_775[k]
                    - f_5 * osg1_775[k]
                    + f_3 * pc_x[k] * osh_1081[k];

        t_1439[k] = f_21 * nsh_867[k]
                    + f_3 * pc_z[k] * osh_1077[k];
    }

#pragma omp simd aligned(t_1440, t_1441, t_1442, pc_x, pc_y, nsh_891, nsh_1083, nsh_1085, \
                         osg0_777, osg0_779, osg1_777, osg1_779, osh_1080, osh_1083, \
                         osh_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1440[k] = f_12 * nsh_1083[k]
                    + f_4 * osg0_777[k]
                    - f_5 * osg1_777[k]
                    + f_3 * pc_x[k] * osh_1083[k];

        t_1441[k] = f_13 * nsh_891[k]
                    + f_3 * pc_y[k] * osh_1080[k];

        t_1442[k] = f_12 * nsh_1085[k]
                    + f_4 * osg0_779[k]
                    - f_5 * osg1_779[k]
                    + f_3 * pc_x[k] * osh_1085[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, t_1446, t_1447, pc_x, nsh_1086, nsh_1087, \
                         nsh_1088, nsh_1089, nsh_1090, osh_1086, osh_1087, osh_1088, osh_1089, \
                         osh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = f_12 * nsh_1086[k]
                    + f_3 * pc_x[k] * osh_1086[k];

        t_1444[k] = f_12 * nsh_1087[k]
                    + f_3 * pc_x[k] * osh_1087[k];

        t_1445[k] = f_12 * nsh_1088[k]
                    + f_3 * pc_x[k] * osh_1088[k];

        t_1446[k] = f_12 * nsh_1089[k]
                    + f_3 * pc_x[k] * osh_1089[k];

        t_1447[k] = f_12 * nsh_1090[k]
                    + f_3 * pc_x[k] * osh_1090[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, pc_z, nsh_876, nsh_897, nsh_1091, \
                         osg0_775, osg1_775, osh_1086, osh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_12 * nsh_1091[k]
                    + f_3 * pc_x[k] * osh_1091[k];

        t_1449[k] = f_13 * nsh_897[k]
                    + f_1 * osg0_775[k]
                    - f_2 * osg1_775[k]
                    + f_3 * pc_y[k] * osh_1086[k];

        t_1450[k] = f_21 * nsh_876[k]
                    + f_3 * pc_z[k] * osh_1086[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_y, nsh_899, nsh_900, nsh_901, osg0_777, \
                         osg0_778, osg0_779, osg1_777, osg1_778, osg1_779, osh_1088, osh_1089, \
                         osh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_13 * nsh_899[k]
                    + f_8 * osg0_777[k]
                    - f_9 * osg1_777[k]
                    + f_3 * pc_y[k] * osh_1088[k];

        t_1452[k] = f_13 * nsh_900[k]
                    + f_6 * osg0_778[k]
                    - f_7 * osg1_778[k]
                    + f_3 * pc_y[k] * osh_1089[k];

        t_1453[k] = f_13 * nsh_901[k]
                    + f_4 * osg0_779[k]
                    - f_5 * osg1_779[k]
                    + f_3 * pc_y[k] * osh_1090[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_y, pc_z, nsh_881, nsh_902, nsh_1092, \
                         osg0_779, osg0_780, osg1_779, osg1_780, osh_1091, \
                         osh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * nsh_902[k]
                    + f_3 * pc_y[k] * osh_1091[k];

        t_1455[k] = f_21 * nsh_881[k]
                    + f_1 * osg0_779[k]
                    - f_2 * osg1_779[k]
                    + f_3 * pc_z[k] * osh_1091[k];

        t_1456[k] = f_12 * nsh_1092[k]
                    + f_1 * osg0_780[k]
                    - f_2 * osg1_780[k]
                    + f_3 * pc_x[k] * osh_1092[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, t_1460, pc_x, pc_y, pc_z, nsh_882, nsh_903, \
                         nsh_905, nsh_1095, osg0_783, osg1_783, osh_1092, osh_1094, \
                         osh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_12 * nsh_903[k]
                    + f_3 * pc_y[k] * osh_1092[k];

        t_1458[k] = f_20 * nsh_882[k]
                    + f_3 * pc_z[k] * osh_1092[k];

        t_1459[k] = f_12 * nsh_1095[k]
                    + f_8 * osg0_783[k]
                    - f_9 * osg1_783[k]
                    + f_3 * pc_x[k] * osh_1095[k];

        t_1460[k] = f_12 * nsh_905[k]
                    + f_3 * pc_y[k] * osh_1094[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pc_x, pc_z, nsh_885, nsh_1097, nsh_1098, \
                         osg0_785, osg0_786, osg1_785, osg1_786, osh_1095, osh_1097, \
                         osh_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = f_12 * nsh_1097[k]
                    + f_8 * osg0_785[k]
                    - f_9 * osg1_785[k]
                    + f_3 * pc_x[k] * osh_1097[k];

        t_1462[k] = f_12 * nsh_1098[k]
                    + f_6 * osg0_786[k]
                    - f_7 * osg1_786[k]
                    + f_3 * pc_x[k] * osh_1098[k];

        t_1463[k] = f_20 * nsh_885[k]
                    + f_3 * pc_z[k] * osh_1095[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pc_x, pc_y, nsh_908, nsh_1101, nsh_1102, \
                         osg0_789, osg0_790, osg1_789, osg1_790, osh_1097, osh_1101, \
                         osh_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_12 * nsh_908[k]
                    + f_3 * pc_y[k] * osh_1097[k];

        t_1465[k] = f_12 * nsh_1101[k]
                    + f_6 * osg0_789[k]
                    - f_7 * osg1_789[k]
                    + f_3 * pc_x[k] * osh_1101[k];

        t_1466[k] = f_12 * nsh_1102[k]
                    + f_4 * osg0_790[k]
                    - f_5 * osg1_790[k]
                    + f_3 * pc_x[k] * osh_1102[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, pc_x, pc_y, pc_z, nsh_888, nsh_912, nsh_1104, \
                         osg0_792, osg1_792, osh_1098, osh_1101, \
                         osh_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_20 * nsh_888[k]
                    + f_3 * pc_z[k] * osh_1098[k];

        t_1468[k] = f_12 * nsh_1104[k]
                    + f_4 * osg0_792[k]
                    - f_5 * osg1_792[k]
                    + f_3 * pc_x[k] * osh_1104[k];

        t_1469[k] = f_12 * nsh_912[k]
                    + f_3 * pc_y[k] * osh_1101[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, pc_x, nsh_1106, nsh_1107, nsh_1108, \
                         nsh_1109, osg0_794, osg1_794, osh_1106, osh_1107, osh_1108, \
                         osh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_12 * nsh_1106[k]
                    + f_4 * osg0_794[k]
                    - f_5 * osg1_794[k]
                    + f_3 * pc_x[k] * osh_1106[k];

        t_1471[k] = f_12 * nsh_1107[k]
                    + f_3 * pc_x[k] * osh_1107[k];

        t_1472[k] = f_12 * nsh_1108[k]
                    + f_3 * pc_x[k] * osh_1108[k];

        t_1473[k] = f_12 * nsh_1109[k]
                    + f_3 * pc_x[k] * osh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pc_x, pc_y, nsh_918, nsh_1110, \
                         nsh_1111, nsh_1112, osg0_790, osg1_790, osh_1107, osh_1110, osh_1111, \
                         osh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_12 * nsh_1110[k]
                    + f_3 * pc_x[k] * osh_1110[k];

        t_1475[k] = f_12 * nsh_1111[k]
                    + f_3 * pc_x[k] * osh_1111[k];

        t_1476[k] = f_12 * nsh_1112[k]
                    + f_3 * pc_x[k] * osh_1112[k];

        t_1477[k] = f_12 * nsh_918[k]
                    + f_1 * osg0_790[k]
                    - f_2 * osg1_790[k]
                    + f_3 * pc_y[k] * osh_1107[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osg0, const size_t osg1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);
    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);
    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);
    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);
    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);
    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_1232 = buffer.data(nsi0 + 1232);
    const auto *nsi0_1235 = buffer.data(nsi0 + 1235);
    const auto *nsi0_1237 = buffer.data(nsi0 + 1237);
    const auto *nsi0_1238 = buffer.data(nsi0 + 1238);
    const auto *nsi0_1241 = buffer.data(nsi0 + 1241);
    const auto *nsi0_1242 = buffer.data(nsi0 + 1242);
    const auto *nsi0_1244 = buffer.data(nsi0 + 1244);
    const auto *nsi0_1246 = buffer.data(nsi0 + 1246);
    const auto *nsi0_1259 = buffer.data(nsi0 + 1259);
    const auto *nsi0_1260 = buffer.data(nsi0 + 1260);
    const auto *nsi0_1263 = buffer.data(nsi0 + 1263);
    const auto *nsi0_1266 = buffer.data(nsi0 + 1266);
    const auto *nsi0_1270 = buffer.data(nsi0 + 1270);
    const auto *nsi0_1540 = buffer.data(nsi0 + 1540);
    const auto *nsi0_1543 = buffer.data(nsi0 + 1543);
    const auto *nsi0_1546 = buffer.data(nsi0 + 1546);
    const auto *nsi0_1550 = buffer.data(nsi0 + 1550);
    const auto *nsi0_1561 = buffer.data(nsi0 + 1561);
    const auto *nsi0_1563 = buffer.data(nsi0 + 1563);
    const auto *nsi0_1564 = buffer.data(nsi0 + 1564);
    const auto *nsi0_1565 = buffer.data(nsi0 + 1565);
    const auto *nsi0_1567 = buffer.data(nsi0 + 1567);
    const auto *nsi0_1573 = buffer.data(nsi0 + 1573);
    const auto *nsi0_1577 = buffer.data(nsi0 + 1577);
    const auto *nsi0_1580 = buffer.data(nsi0 + 1580);
    const auto *nsi0_1582 = buffer.data(nsi0 + 1582);
    const auto *nsi0_1589 = buffer.data(nsi0 + 1589);
    const auto *nsi0_1591 = buffer.data(nsi0 + 1591);
    const auto *nsi0_1592 = buffer.data(nsi0 + 1592);
    const auto *nsi0_1593 = buffer.data(nsi0 + 1593);
    const auto *nsi0_1595 = buffer.data(nsi0 + 1595);
    const auto *nsi0_1596 = buffer.data(nsi0 + 1596);
    const auto *nsi0_1599 = buffer.data(nsi0 + 1599);

    const auto *nsh_897 = buffer.data(nsh + 897);
    const auto *nsh_902 = buffer.data(nsh + 902);
    const auto *nsh_903 = buffer.data(nsh + 903);
    const auto *nsh_906 = buffer.data(nsh + 906);
    const auto *nsh_909 = buffer.data(nsh + 909);
    const auto *nsh_918 = buffer.data(nsh + 918);
    const auto *nsh_920 = buffer.data(nsh + 920);
    const auto *nsh_921 = buffer.data(nsh + 921);
    const auto *nsh_922 = buffer.data(nsh + 922);
    const auto *nsh_923 = buffer.data(nsh + 923);
    const auto *nsh_924 = buffer.data(nsh + 924);
    const auto *nsh_925 = buffer.data(nsh + 925);
    const auto *nsh_926 = buffer.data(nsh + 926);
    const auto *nsh_927 = buffer.data(nsh + 927);
    const auto *nsh_929 = buffer.data(nsh + 929);
    const auto *nsh_930 = buffer.data(nsh + 930);
    const auto *nsh_932 = buffer.data(nsh + 932);
    const auto *nsh_933 = buffer.data(nsh + 933);
    const auto *nsh_939 = buffer.data(nsh + 939);
    const auto *nsh_941 = buffer.data(nsh + 941);
    const auto *nsh_942 = buffer.data(nsh + 942);
    const auto *nsh_943 = buffer.data(nsh + 943);
    const auto *nsh_944 = buffer.data(nsh + 944);
    const auto *nsh_945 = buffer.data(nsh + 945);
    const auto *nsh_948 = buffer.data(nsh + 948);
    const auto *nsh_950 = buffer.data(nsh + 950);
    const auto *nsh_951 = buffer.data(nsh + 951);
    const auto *nsh_954 = buffer.data(nsh + 954);
    const auto *nsh_960 = buffer.data(nsh + 960);
    const auto *nsh_965 = buffer.data(nsh + 965);
    const auto *nsh_966 = buffer.data(nsh + 966);
    const auto *nsh_968 = buffer.data(nsh + 968);
    const auto *nsh_971 = buffer.data(nsh + 971);
    const auto *nsh_975 = buffer.data(nsh + 975);
    const auto *nsh_986 = buffer.data(nsh + 986);
    const auto *nsh_987 = buffer.data(nsh + 987);
    const auto *nsh_989 = buffer.data(nsh + 989);
    const auto *nsh_1128 = buffer.data(nsh + 1128);
    const auto *nsh_1129 = buffer.data(nsh + 1129);
    const auto *nsh_1130 = buffer.data(nsh + 1130);
    const auto *nsh_1131 = buffer.data(nsh + 1131);
    const auto *nsh_1132 = buffer.data(nsh + 1132);
    const auto *nsh_1133 = buffer.data(nsh + 1133);
    const auto *nsh_1134 = buffer.data(nsh + 1134);
    const auto *nsh_1139 = buffer.data(nsh + 1139);
    const auto *nsh_1143 = buffer.data(nsh + 1143);
    const auto *nsh_1148 = buffer.data(nsh + 1148);
    const auto *nsh_1149 = buffer.data(nsh + 1149);
    const auto *nsh_1150 = buffer.data(nsh + 1150);
    const auto *nsh_1151 = buffer.data(nsh + 1151);
    const auto *nsh_1152 = buffer.data(nsh + 1152);
    const auto *nsh_1154 = buffer.data(nsh + 1154);
    const auto *nsh_1155 = buffer.data(nsh + 1155);
    const auto *nsh_1158 = buffer.data(nsh + 1158);
    const auto *nsh_1161 = buffer.data(nsh + 1161);
    const auto *nsh_1165 = buffer.data(nsh + 1165);
    const auto *nsh_1170 = buffer.data(nsh + 1170);
    const auto *nsh_1172 = buffer.data(nsh + 1172);
    const auto *nsh_1173 = buffer.data(nsh + 1173);
    const auto *nsh_1174 = buffer.data(nsh + 1174);
    const auto *nsh_1175 = buffer.data(nsh + 1175);
    const auto *nsh_1181 = buffer.data(nsh + 1181);
    const auto *nsh_1185 = buffer.data(nsh + 1185);
    const auto *nsh_1188 = buffer.data(nsh + 1188);
    const auto *nsh_1190 = buffer.data(nsh + 1190);
    const auto *nsh_1191 = buffer.data(nsh + 1191);
    const auto *nsh_1192 = buffer.data(nsh + 1192);
    const auto *nsh_1193 = buffer.data(nsh + 1193);
    const auto *nsh_1194 = buffer.data(nsh + 1194);
    const auto *nsh_1195 = buffer.data(nsh + 1195);
    const auto *nsh_1196 = buffer.data(nsh + 1196);
    const auto *nsh_1197 = buffer.data(nsh + 1197);
    const auto *nsh_1200 = buffer.data(nsh + 1200);

    const auto *nsi1_1232 = buffer.data(nsi1 + 1232);
    const auto *nsi1_1235 = buffer.data(nsi1 + 1235);
    const auto *nsi1_1237 = buffer.data(nsi1 + 1237);
    const auto *nsi1_1238 = buffer.data(nsi1 + 1238);
    const auto *nsi1_1241 = buffer.data(nsi1 + 1241);
    const auto *nsi1_1242 = buffer.data(nsi1 + 1242);
    const auto *nsi1_1244 = buffer.data(nsi1 + 1244);
    const auto *nsi1_1246 = buffer.data(nsi1 + 1246);
    const auto *nsi1_1259 = buffer.data(nsi1 + 1259);
    const auto *nsi1_1260 = buffer.data(nsi1 + 1260);
    const auto *nsi1_1263 = buffer.data(nsi1 + 1263);
    const auto *nsi1_1266 = buffer.data(nsi1 + 1266);
    const auto *nsi1_1270 = buffer.data(nsi1 + 1270);
    const auto *nsi1_1540 = buffer.data(nsi1 + 1540);
    const auto *nsi1_1543 = buffer.data(nsi1 + 1543);
    const auto *nsi1_1546 = buffer.data(nsi1 + 1546);
    const auto *nsi1_1550 = buffer.data(nsi1 + 1550);
    const auto *nsi1_1561 = buffer.data(nsi1 + 1561);
    const auto *nsi1_1563 = buffer.data(nsi1 + 1563);
    const auto *nsi1_1564 = buffer.data(nsi1 + 1564);
    const auto *nsi1_1565 = buffer.data(nsi1 + 1565);
    const auto *nsi1_1567 = buffer.data(nsi1 + 1567);
    const auto *nsi1_1573 = buffer.data(nsi1 + 1573);
    const auto *nsi1_1577 = buffer.data(nsi1 + 1577);
    const auto *nsi1_1580 = buffer.data(nsi1 + 1580);
    const auto *nsi1_1582 = buffer.data(nsi1 + 1582);
    const auto *nsi1_1589 = buffer.data(nsi1 + 1589);
    const auto *nsi1_1591 = buffer.data(nsi1 + 1591);
    const auto *nsi1_1592 = buffer.data(nsi1 + 1592);
    const auto *nsi1_1593 = buffer.data(nsi1 + 1593);
    const auto *nsi1_1595 = buffer.data(nsi1 + 1595);
    const auto *nsi1_1596 = buffer.data(nsi1 + 1596);
    const auto *nsi1_1599 = buffer.data(nsi1 + 1599);

    const auto *osg0_792 = buffer.data(osg0 + 792);
    const auto *osg0_793 = buffer.data(osg0 + 793);
    const auto *osg0_794 = buffer.data(osg0 + 794);
    const auto *osg0_805 = buffer.data(osg0 + 805);
    const auto *osg0_807 = buffer.data(osg0 + 807);
    const auto *osg0_808 = buffer.data(osg0 + 808);
    const auto *osg0_809 = buffer.data(osg0 + 809);
    const auto *osg0_810 = buffer.data(osg0 + 810);
    const auto *osg0_811 = buffer.data(osg0 + 811);
    const auto *osg0_812 = buffer.data(osg0 + 812);
    const auto *osg0_813 = buffer.data(osg0 + 813);
    const auto *osg0_814 = buffer.data(osg0 + 814);
    const auto *osg0_815 = buffer.data(osg0 + 815);
    const auto *osg0_819 = buffer.data(osg0 + 819);
    const auto *osg0_820 = buffer.data(osg0 + 820);
    const auto *osg0_821 = buffer.data(osg0 + 821);
    const auto *osg0_822 = buffer.data(osg0 + 822);
    const auto *osg0_823 = buffer.data(osg0 + 823);
    const auto *osg0_824 = buffer.data(osg0 + 824);
    const auto *osg0_825 = buffer.data(osg0 + 825);
    const auto *osg0_827 = buffer.data(osg0 + 827);
    const auto *osg0_828 = buffer.data(osg0 + 828);
    const auto *osg0_830 = buffer.data(osg0 + 830);

    const auto *osg1_792 = buffer.data(osg1 + 792);
    const auto *osg1_793 = buffer.data(osg1 + 793);
    const auto *osg1_794 = buffer.data(osg1 + 794);
    const auto *osg1_805 = buffer.data(osg1 + 805);
    const auto *osg1_807 = buffer.data(osg1 + 807);
    const auto *osg1_808 = buffer.data(osg1 + 808);
    const auto *osg1_809 = buffer.data(osg1 + 809);
    const auto *osg1_810 = buffer.data(osg1 + 810);
    const auto *osg1_811 = buffer.data(osg1 + 811);
    const auto *osg1_812 = buffer.data(osg1 + 812);
    const auto *osg1_813 = buffer.data(osg1 + 813);
    const auto *osg1_814 = buffer.data(osg1 + 814);
    const auto *osg1_815 = buffer.data(osg1 + 815);
    const auto *osg1_819 = buffer.data(osg1 + 819);
    const auto *osg1_820 = buffer.data(osg1 + 820);
    const auto *osg1_821 = buffer.data(osg1 + 821);
    const auto *osg1_822 = buffer.data(osg1 + 822);
    const auto *osg1_823 = buffer.data(osg1 + 823);
    const auto *osg1_824 = buffer.data(osg1 + 824);
    const auto *osg1_825 = buffer.data(osg1 + 825);
    const auto *osg1_827 = buffer.data(osg1 + 827);
    const auto *osg1_828 = buffer.data(osg1 + 828);
    const auto *osg1_830 = buffer.data(osg1 + 830);

    const auto *osh_1107 = buffer.data(osh + 1107);
    const auto *osh_1109 = buffer.data(osh + 1109);
    const auto *osh_1110 = buffer.data(osh + 1110);
    const auto *osh_1111 = buffer.data(osh + 1111);
    const auto *osh_1112 = buffer.data(osh + 1112);
    const auto *osh_1113 = buffer.data(osh + 1113);
    const auto *osh_1115 = buffer.data(osh + 1115);
    const auto *osh_1116 = buffer.data(osh + 1116);
    const auto *osh_1118 = buffer.data(osh + 1118);
    const auto *osh_1119 = buffer.data(osh + 1119);
    const auto *osh_1122 = buffer.data(osh + 1122);
    const auto *osh_1128 = buffer.data(osh + 1128);
    const auto *osh_1129 = buffer.data(osh + 1129);
    const auto *osh_1130 = buffer.data(osh + 1130);
    const auto *osh_1131 = buffer.data(osh + 1131);
    const auto *osh_1132 = buffer.data(osh + 1132);
    const auto *osh_1133 = buffer.data(osh + 1133);
    const auto *osh_1134 = buffer.data(osh + 1134);
    const auto *osh_1135 = buffer.data(osh + 1135);
    const auto *osh_1136 = buffer.data(osh + 1136);
    const auto *osh_1137 = buffer.data(osh + 1137);
    const auto *osh_1138 = buffer.data(osh + 1138);
    const auto *osh_1139 = buffer.data(osh + 1139);
    const auto *osh_1140 = buffer.data(osh + 1140);
    const auto *osh_1141 = buffer.data(osh + 1141);
    const auto *osh_1142 = buffer.data(osh + 1142);
    const auto *osh_1143 = buffer.data(osh + 1143);
    const auto *osh_1148 = buffer.data(osh + 1148);
    const auto *osh_1149 = buffer.data(osh + 1149);
    const auto *osh_1150 = buffer.data(osh + 1150);
    const auto *osh_1151 = buffer.data(osh + 1151);
    const auto *osh_1152 = buffer.data(osh + 1152);
    const auto *osh_1153 = buffer.data(osh + 1153);
    const auto *osh_1154 = buffer.data(osh + 1154);
    const auto *osh_1155 = buffer.data(osh + 1155);
    const auto *osh_1156 = buffer.data(osh + 1156);
    const auto *osh_1157 = buffer.data(osh + 1157);
    const auto *osh_1158 = buffer.data(osh + 1158);
    const auto *osh_1160 = buffer.data(osh + 1160);
    const auto *osh_1161 = buffer.data(osh + 1161);
    const auto *osh_1162 = buffer.data(osh + 1162);
    const auto *osh_1164 = buffer.data(osh + 1164);
    const auto *osh_1165 = buffer.data(osh + 1165);
    const auto *osh_1170 = buffer.data(osh + 1170);
    const auto *osh_1172 = buffer.data(osh + 1172);
    const auto *osh_1173 = buffer.data(osh + 1173);
    const auto *osh_1174 = buffer.data(osh + 1174);
    const auto *osh_1175 = buffer.data(osh + 1175);
    const auto *osh_1176 = buffer.data(osh + 1176);
    const auto *osh_1178 = buffer.data(osh + 1178);
    const auto *osh_1179 = buffer.data(osh + 1179);
    const auto *osh_1181 = buffer.data(osh + 1181);
    const auto *osh_1182 = buffer.data(osh + 1182);
    const auto *osh_1185 = buffer.data(osh + 1185);
    const auto *osh_1191 = buffer.data(osh + 1191);
    const auto *osh_1192 = buffer.data(osh + 1192);
    const auto *osh_1193 = buffer.data(osh + 1193);
    const auto *osh_1194 = buffer.data(osh + 1194);
    const auto *osh_1195 = buffer.data(osh + 1195);
    const auto *osh_1196 = buffer.data(osh + 1196);
    const auto *osh_1197 = buffer.data(osh + 1197);
    const auto *osh_1199 = buffer.data(osh + 1199);

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_y, pc_z, nsh_897, nsh_920, nsh_921, \
                         osg0_792, osg0_793, osg1_792, osg1_793, osh_1107, osh_1109, \
                         osh_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_20 * nsh_897[k]
                    + f_3 * pc_z[k] * osh_1107[k];

        t_1479[k] = f_12 * nsh_920[k]
                    + f_8 * osg0_792[k]
                    - f_9 * osg1_792[k]
                    + f_3 * pc_y[k] * osh_1109[k];

        t_1480[k] = f_12 * nsh_921[k]
                    + f_6 * osg0_793[k]
                    - f_7 * osg1_793[k]
                    + f_3 * pc_y[k] * osh_1110[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, t_1484, pa_y, pc_y, pc_z, nsi0_1232, nsh_902, \
                         nsh_922, nsh_923, nsi1_1232, osg0_794, osg1_794, osh_1111, \
                         osh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_12 * nsh_922[k]
                    + f_4 * osg0_794[k]
                    - f_5 * osg1_794[k]
                    + f_3 * pc_y[k] * osh_1111[k];

        t_1482[k] = f_12 * nsh_923[k]
                    + f_3 * pc_y[k] * osh_1112[k];

        t_1483[k] = f_20 * nsh_902[k]
                    + f_1 * osg0_794[k]
                    - f_2 * osg1_794[k]
                    + f_3 * pc_z[k] * osh_1112[k];

        t_1484[k] = pa_y[k] * nsi0_1232[k]
                    - f_10 * pc_y[k] * nsi1_1232[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, t_1488, pa_y, pc_y, pc_z, nsi0_1235, nsh_903, \
                         nsh_924, nsh_925, nsh_926, nsi1_1235, osh_1113, \
                         osh_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_11 * nsh_924[k]
                    + f_3 * pc_y[k] * osh_1113[k];

        t_1486[k] = f_19 * nsh_903[k]
                    + f_3 * pc_z[k] * osh_1113[k];

        t_1487[k] = pa_y[k] * nsi0_1235[k]
                    + f_12 * nsh_925[k]
                    - f_10 * pc_y[k] * nsi1_1235[k];

        t_1488[k] = f_11 * nsh_926[k]
                    + f_3 * pc_y[k] * osh_1115[k];
    }

#pragma omp simd aligned(t_1489, t_1490, t_1491, t_1492, pa_y, pc_y, pc_z, nsi0_1237, \
                         nsi0_1238, nsh_906, nsh_927, nsh_929, nsi1_1237, nsi1_1238, osh_1116, \
                         osh_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1489[k] = pa_y[k] * nsi0_1237[k]
                    - f_10 * pc_y[k] * nsi1_1237[k];

        t_1490[k] = pa_y[k] * nsi0_1238[k]
                    + f_13 * nsh_927[k]
                    - f_10 * pc_y[k] * nsi1_1238[k];

        t_1491[k] = f_19 * nsh_906[k]
                    + f_3 * pc_z[k] * osh_1116[k];

        t_1492[k] = f_11 * nsh_929[k]
                    + f_3 * pc_y[k] * osh_1118[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pa_y, pc_y, pc_z, nsi0_1241, nsi0_1242, \
                         nsh_909, nsh_930, nsi1_1241, nsi1_1242, \
                         osh_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = pa_y[k] * nsi0_1241[k]
                    - f_10 * pc_y[k] * nsi1_1241[k];

        t_1494[k] = pa_y[k] * nsi0_1242[k]
                    + f_14 * nsh_930[k]
                    - f_10 * pc_y[k] * nsi1_1242[k];

        t_1495[k] = f_19 * nsh_909[k]
                    + f_3 * pc_z[k] * osh_1119[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pa_y, pc_x, pc_y, nsi0_1244, \
                         nsi0_1246, nsh_932, nsh_933, nsh_1128, nsi1_1244, nsi1_1246, \
                         osh_1122, osh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = pa_y[k] * nsi0_1244[k]
                    + f_12 * nsh_932[k]
                    - f_10 * pc_y[k] * nsi1_1244[k];

        t_1497[k] = f_11 * nsh_933[k]
                    + f_3 * pc_y[k] * osh_1122[k];

        t_1498[k] = pa_y[k] * nsi0_1246[k]
                    - f_10 * pc_y[k] * nsi1_1246[k];

        t_1499[k] = f_12 * nsh_1128[k]
                    + f_3 * pc_x[k] * osh_1128[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, t_1504, pc_x, nsh_1129, nsh_1130, \
                         nsh_1131, nsh_1132, nsh_1133, osh_1129, osh_1130, osh_1131, osh_1132, \
                         osh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_12 * nsh_1129[k]
                    + f_3 * pc_x[k] * osh_1129[k];

        t_1501[k] = f_12 * nsh_1130[k]
                    + f_3 * pc_x[k] * osh_1130[k];

        t_1502[k] = f_12 * nsh_1131[k]
                    + f_3 * pc_x[k] * osh_1131[k];

        t_1503[k] = f_12 * nsh_1132[k]
                    + f_3 * pc_x[k] * osh_1132[k];

        t_1504[k] = f_12 * nsh_1133[k]
                    + f_3 * pc_x[k] * osh_1133[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pc_y, pc_z, nsh_918, nsh_939, nsh_941, \
                         osg0_805, osg0_807, osg1_805, osg1_807, osh_1128, \
                         osh_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_11 * nsh_939[k]
                    + f_1 * osg0_805[k]
                    - f_2 * osg1_805[k]
                    + f_3 * pc_y[k] * osh_1128[k];

        t_1506[k] = f_19 * nsh_918[k]
                    + f_3 * pc_z[k] * osh_1128[k];

        t_1507[k] = f_11 * nsh_941[k]
                    + f_8 * osg0_807[k]
                    - f_9 * osg1_807[k]
                    + f_3 * pc_y[k] * osh_1130[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pc_y, nsh_942, nsh_943, nsh_944, osg0_808, \
                         osg0_809, osg1_808, osg1_809, osh_1131, osh_1132, \
                         osh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_11 * nsh_942[k]
                    + f_6 * osg0_808[k]
                    - f_7 * osg1_808[k]
                    + f_3 * pc_y[k] * osh_1131[k];

        t_1509[k] = f_11 * nsh_943[k]
                    + f_4 * osg0_809[k]
                    - f_5 * osg1_809[k]
                    + f_3 * pc_y[k] * osh_1132[k];

        t_1510[k] = f_11 * nsh_944[k]
                    + f_3 * pc_y[k] * osh_1133[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pa_y, pc_x, pc_y, pc_z, nsi0_1259, \
                         nsh_924, nsh_1134, nsi1_1259, osg0_810, osg1_810, \
                         osh_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = pa_y[k] * nsi0_1259[k]
                    - f_10 * pc_y[k] * nsi1_1259[k];

        t_1512[k] = f_12 * nsh_1134[k]
                    + f_1 * osg0_810[k]
                    - f_2 * osg1_810[k]
                    + f_3 * pc_x[k] * osh_1134[k];

        t_1513[k] = f_3 * pc_y[k] * osh_1134[k];

        t_1514[k] = f_18 * nsh_924[k]
                    + f_3 * pc_z[k] * osh_1134[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, pc_x, pc_y, nsh_1139, osg0_810, osg0_815, \
                         osg1_810, osg1_815, osh_1135, osh_1136, \
                         osh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_4 * osg0_810[k]
                    - f_5 * osg1_810[k]
                    + f_3 * pc_y[k] * osh_1135[k];

        t_1516[k] = f_3 * pc_y[k] * osh_1136[k];

        t_1517[k] = f_12 * nsh_1139[k]
                    + f_8 * osg0_815[k]
                    - f_9 * osg1_815[k]
                    + f_3 * pc_x[k] * osh_1139[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, pc_y, osg0_811, osg0_812, osg1_811, osg1_812, \
                         osh_1137, osh_1138, osh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_6 * osg0_811[k]
                    - f_7 * osg1_811[k]
                    + f_3 * pc_y[k] * osh_1137[k];

        t_1519[k] = f_4 * osg0_812[k]
                    - f_5 * osg1_812[k]
                    + f_3 * pc_y[k] * osh_1138[k];

        t_1520[k] = f_3 * pc_y[k] * osh_1139[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, pc_x, pc_y, nsh_1143, osg0_813, osg0_814, \
                         osg0_819, osg1_813, osg1_814, osg1_819, osh_1140, osh_1141, \
                         osh_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_12 * nsh_1143[k]
                    + f_6 * osg0_819[k]
                    - f_7 * osg1_819[k]
                    + f_3 * pc_x[k] * osh_1143[k];

        t_1522[k] = f_8 * osg0_813[k]
                    - f_9 * osg1_813[k]
                    + f_3 * pc_y[k] * osh_1140[k];

        t_1523[k] = f_6 * osg0_814[k]
                    - f_7 * osg1_814[k]
                    + f_3 * pc_y[k] * osh_1141[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, t_1527, pc_x, pc_y, nsh_1148, nsh_1149, \
                         osg0_815, osg0_824, osg1_815, osg1_824, osh_1142, osh_1143, osh_1148, \
                         osh_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = f_4 * osg0_815[k]
                    - f_5 * osg1_815[k]
                    + f_3 * pc_y[k] * osh_1142[k];

        t_1525[k] = f_3 * pc_y[k] * osh_1143[k];

        t_1526[k] = f_12 * nsh_1148[k]
                    + f_4 * osg0_824[k]
                    - f_5 * osg1_824[k]
                    + f_3 * pc_x[k] * osh_1148[k];

        t_1527[k] = f_12 * nsh_1149[k]
                    + f_3 * pc_x[k] * osh_1149[k];
    }

#pragma omp simd aligned(t_1528, t_1529, t_1530, t_1531, t_1532, pc_x, pc_y, nsh_1150, \
                         nsh_1151, nsh_1152, nsh_1154, osh_1148, osh_1150, osh_1151, osh_1152, \
                         osh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1528[k] = f_12 * nsh_1150[k]
                    + f_3 * pc_x[k] * osh_1150[k];

        t_1529[k] = f_12 * nsh_1151[k]
                    + f_3 * pc_x[k] * osh_1151[k];

        t_1530[k] = f_12 * nsh_1152[k]
                    + f_3 * pc_x[k] * osh_1152[k];

        t_1531[k] = f_3 * pc_y[k] * osh_1148[k];

        t_1532[k] = f_12 * nsh_1154[k]
                    + f_3 * pc_x[k] * osh_1154[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, pc_y, osg0_820, osg0_821, osg0_822, osg1_820, \
                         osg1_821, osg1_822, osh_1149, osh_1150, \
                         osh_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = f_1 * osg0_820[k]
                    - f_2 * osg1_820[k]
                    + f_3 * pc_y[k] * osh_1149[k];

        t_1534[k] = f_16 * osg0_821[k]
                    - f_17 * osg1_821[k]
                    + f_3 * pc_y[k] * osh_1150[k];

        t_1535[k] = f_8 * osg0_822[k]
                    - f_9 * osg1_822[k]
                    + f_3 * pc_y[k] * osh_1151[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pc_y, pc_z, nsh_944, osg0_823, \
                         osg0_824, osg1_823, osg1_824, osh_1152, osh_1153, \
                         osh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_6 * osg0_823[k]
                    - f_7 * osg1_823[k]
                    + f_3 * pc_y[k] * osh_1152[k];

        t_1537[k] = f_4 * osg0_824[k]
                    - f_5 * osg1_824[k]
                    + f_3 * pc_y[k] * osh_1153[k];

        t_1538[k] = f_3 * pc_y[k] * osh_1154[k];

        t_1539[k] = f_18 * nsh_944[k]
                    + f_1 * osg0_824[k]
                    - f_2 * osg1_824[k]
                    + f_3 * pc_z[k] * osh_1154[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, t_1543, pa_x, pc_x, pc_y, pc_z, nsi0_1540, \
                         nsi0_1543, nsh_945, nsh_1155, nsh_1158, nsi1_1540, nsi1_1543, \
                         osh_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = pa_x[k] * nsi0_1540[k]
                    + f_21 * nsh_1155[k]
                    - f_10 * pc_x[k] * nsi1_1540[k];

        t_1541[k] = f_15 * nsh_945[k]
                    + f_3 * pc_y[k] * osh_1155[k];

        t_1542[k] = f_3 * pc_z[k] * osh_1155[k];

        t_1543[k] = pa_x[k] * nsi0_1543[k]
                    + f_14 * nsh_1158[k]
                    - f_10 * pc_x[k] * nsi1_1543[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, t_1547, pa_x, pc_x, pc_z, nsi0_1546, \
                         nsh_1161, nsi1_1546, osg0_825, osg1_825, osh_1156, osh_1157, \
                         osh_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = f_3 * pc_z[k] * osh_1156[k];

        t_1545[k] = f_4 * osg0_825[k]
                    - f_5 * osg1_825[k]
                    + f_3 * pc_z[k] * osh_1157[k];

        t_1546[k] = pa_x[k] * nsi0_1546[k]
                    + f_13 * nsh_1161[k]
                    - f_10 * pc_x[k] * nsi1_1546[k];

        t_1547[k] = f_3 * pc_z[k] * osh_1158[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, t_1551, pa_x, pc_x, pc_y, pc_z, nsi0_1550, \
                         nsh_950, nsh_1165, nsi1_1550, osg0_827, osg1_827, osh_1160, \
                         osh_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = f_15 * nsh_950[k]
                    + f_3 * pc_y[k] * osh_1160[k];

        t_1549[k] = f_6 * osg0_827[k]
                    - f_7 * osg1_827[k]
                    + f_3 * pc_z[k] * osh_1160[k];

        t_1550[k] = pa_x[k] * nsi0_1550[k]
                    + f_12 * nsh_1165[k]
                    - f_10 * pc_x[k] * nsi1_1550[k];

        t_1551[k] = f_3 * pc_z[k] * osh_1161[k];
    }

#pragma omp simd aligned(t_1552, t_1553, t_1554, t_1555, pc_x, pc_y, pc_z, nsh_954, nsh_1170, \
                         osg0_828, osg0_830, osg1_828, osg1_830, osh_1162, osh_1164, \
                         osh_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1552[k] = f_4 * osg0_828[k]
                    - f_5 * osg1_828[k]
                    + f_3 * pc_z[k] * osh_1162[k];

        t_1553[k] = f_15 * nsh_954[k]
                    + f_3 * pc_y[k] * osh_1164[k];

        t_1554[k] = f_8 * osg0_830[k]
                    - f_9 * osg1_830[k]
                    + f_3 * pc_z[k] * osh_1164[k];

        t_1555[k] = f_11 * nsh_1170[k]
                    + f_3 * pc_x[k] * osh_1170[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, t_1559, t_1560, pc_x, pc_z, nsh_1172, \
                         nsh_1173, nsh_1174, nsh_1175, osh_1165, osh_1172, osh_1173, osh_1174, \
                         osh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_3 * pc_z[k] * osh_1165[k];

        t_1557[k] = f_11 * nsh_1172[k]
                    + f_3 * pc_x[k] * osh_1172[k];

        t_1558[k] = f_11 * nsh_1173[k]
                    + f_3 * pc_x[k] * osh_1173[k];

        t_1559[k] = f_11 * nsh_1174[k]
                    + f_3 * pc_x[k] * osh_1174[k];

        t_1560[k] = f_11 * nsh_1175[k]
                    + f_3 * pc_x[k] * osh_1175[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, pa_x, pc_x, pc_z, nsi0_1561, \
                         nsi0_1563, nsi0_1564, nsi1_1561, nsi1_1563, nsi1_1564, \
                         osh_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = pa_x[k] * nsi0_1561[k]
                    - f_10 * pc_x[k] * nsi1_1561[k];

        t_1562[k] = f_3 * pc_z[k] * osh_1170[k];

        t_1563[k] = pa_x[k] * nsi0_1563[k]
                    - f_10 * pc_x[k] * nsi1_1563[k];

        t_1564[k] = pa_x[k] * nsi0_1564[k]
                    - f_10 * pc_x[k] * nsi1_1564[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, pa_x, pc_x, pc_y, nsi0_1565, nsi0_1567, \
                         nsh_965, nsi1_1565, nsi1_1567, osh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = pa_x[k] * nsi0_1565[k]
                    - f_10 * pc_x[k] * nsi1_1565[k];

        t_1566[k] = f_15 * nsh_965[k]
                    + f_3 * pc_y[k] * osh_1175[k];

        t_1567[k] = pa_x[k] * nsi0_1567[k]
                    - f_10 * pc_x[k] * nsi1_1567[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, pa_z, pc_y, pc_z, nsi0_1260, \
                         nsi0_1263, nsh_945, nsh_966, nsi1_1260, nsi1_1263, \
                         osh_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pa_z[k] * nsi0_1260[k]
                    - f_10 * pc_z[k] * nsi1_1260[k];

        t_1569[k] = f_18 * nsh_966[k]
                    + f_3 * pc_y[k] * osh_1176[k];

        t_1570[k] = f_11 * nsh_945[k]
                    + f_3 * pc_z[k] * osh_1176[k];

        t_1571[k] = pa_z[k] * nsi0_1263[k]
                    - f_10 * pc_z[k] * nsi1_1263[k];
    }

#pragma omp simd aligned(t_1572, t_1573, t_1574, pa_x, pa_z, pc_x, pc_y, pc_z, nsi0_1266, \
                         nsi0_1573, nsh_968, nsh_1181, nsi1_1266, nsi1_1573, \
                         osh_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = f_18 * nsh_968[k]
                    + f_3 * pc_y[k] * osh_1178[k];

        t_1573[k] = pa_x[k] * nsi0_1573[k]
                    + f_14 * nsh_1181[k]
                    - f_10 * pc_x[k] * nsi1_1573[k];

        t_1574[k] = pa_z[k] * nsi0_1266[k]
                    - f_10 * pc_z[k] * nsi1_1266[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, pa_x, pc_x, pc_y, pc_z, nsi0_1577, nsh_948, \
                         nsh_971, nsh_1185, nsi1_1577, osh_1179, \
                         osh_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_11 * nsh_948[k]
                    + f_3 * pc_z[k] * osh_1179[k];

        t_1576[k] = f_18 * nsh_971[k]
                    + f_3 * pc_y[k] * osh_1181[k];

        t_1577[k] = pa_x[k] * nsi0_1577[k]
                    + f_13 * nsh_1185[k]
                    - f_10 * pc_x[k] * nsi1_1577[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pa_x, pa_z, pc_x, pc_z, nsi0_1270, nsi0_1580, \
                         nsh_951, nsh_1188, nsi1_1270, nsi1_1580, \
                         osh_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pa_z[k] * nsi0_1270[k]
                    - f_10 * pc_z[k] * nsi1_1270[k];

        t_1579[k] = f_11 * nsh_951[k]
                    + f_3 * pc_z[k] * osh_1182[k];

        t_1580[k] = pa_x[k] * nsi0_1580[k]
                    + f_12 * nsh_1188[k]
                    - f_10 * pc_x[k] * nsi1_1580[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, t_1584, pa_x, pc_x, pc_y, nsi0_1582, nsh_975, \
                         nsh_1190, nsh_1191, nsh_1192, nsi1_1582, osh_1185, osh_1191, \
                         osh_1192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_18 * nsh_975[k]
                    + f_3 * pc_y[k] * osh_1185[k];

        t_1582[k] = pa_x[k] * nsi0_1582[k]
                    + f_12 * nsh_1190[k]
                    - f_10 * pc_x[k] * nsi1_1582[k];

        t_1583[k] = f_11 * nsh_1191[k]
                    + f_3 * pc_x[k] * osh_1191[k];

        t_1584[k] = f_11 * nsh_1192[k]
                    + f_3 * pc_x[k] * osh_1192[k];
    }

#pragma omp simd aligned(t_1585, t_1586, t_1587, t_1588, pc_x, nsh_1193, nsh_1194, nsh_1195, \
                         nsh_1196, osh_1193, osh_1194, osh_1195, \
                         osh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1585[k] = f_11 * nsh_1193[k]
                    + f_3 * pc_x[k] * osh_1193[k];

        t_1586[k] = f_11 * nsh_1194[k]
                    + f_3 * pc_x[k] * osh_1194[k];

        t_1587[k] = f_11 * nsh_1195[k]
                    + f_3 * pc_x[k] * osh_1195[k];

        t_1588[k] = f_11 * nsh_1196[k]
                    + f_3 * pc_x[k] * osh_1196[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, t_1592, pa_x, pc_x, pc_z, nsi0_1589, \
                         nsi0_1591, nsi0_1592, nsh_960, nsi1_1589, nsi1_1591, nsi1_1592, \
                         osh_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = pa_x[k] * nsi0_1589[k]
                    - f_10 * pc_x[k] * nsi1_1589[k];

        t_1590[k] = f_11 * nsh_960[k]
                    + f_3 * pc_z[k] * osh_1191[k];

        t_1591[k] = pa_x[k] * nsi0_1591[k]
                    - f_10 * pc_x[k] * nsi1_1591[k];

        t_1592[k] = pa_x[k] * nsi0_1592[k]
                    - f_10 * pc_x[k] * nsi1_1592[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, t_1596, pa_x, pc_x, pc_y, nsi0_1593, \
                         nsi0_1595, nsi0_1596, nsh_986, nsh_1197, nsi1_1593, nsi1_1595, \
                         nsi1_1596, osh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = pa_x[k] * nsi0_1593[k]
                    - f_10 * pc_x[k] * nsi1_1593[k];

        t_1594[k] = f_18 * nsh_986[k]
                    + f_3 * pc_y[k] * osh_1196[k];

        t_1595[k] = pa_x[k] * nsi0_1595[k]
                    - f_10 * pc_x[k] * nsi1_1595[k];

        t_1596[k] = pa_x[k] * nsi0_1596[k]
                    + f_21 * nsh_1197[k]
                    - f_10 * pc_x[k] * nsi1_1596[k];
    }

#pragma omp simd aligned(t_1597, t_1598, t_1599, t_1600, pa_x, pc_x, pc_y, pc_z, nsi0_1599, \
                         nsh_966, nsh_987, nsh_989, nsh_1200, nsi1_1599, osh_1197, \
                         osh_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1597[k] = f_19 * nsh_987[k]
                    + f_3 * pc_y[k] * osh_1197[k];

        t_1598[k] = f_12 * nsh_966[k]
                    + f_3 * pc_z[k] * osh_1197[k];

        t_1599[k] = pa_x[k] * nsi0_1599[k]
                    + f_14 * nsh_1200[k]
                    - f_10 * pc_x[k] * nsi1_1599[k];

        t_1600[k] = f_19 * nsh_989[k]
                    + f_3 * pc_y[k] * osh_1199[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);
    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);
    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);
    auto *t_1620 = buffer.data(target + 1620);
    auto *t_1621 = buffer.data(target + 1621);
    auto *t_1622 = buffer.data(target + 1622);
    auto *t_1623 = buffer.data(target + 1623);
    auto *t_1624 = buffer.data(target + 1624);
    auto *t_1625 = buffer.data(target + 1625);
    auto *t_1626 = buffer.data(target + 1626);
    auto *t_1627 = buffer.data(target + 1627);
    auto *t_1628 = buffer.data(target + 1628);
    auto *t_1629 = buffer.data(target + 1629);
    auto *t_1630 = buffer.data(target + 1630);
    auto *t_1631 = buffer.data(target + 1631);
    auto *t_1632 = buffer.data(target + 1632);
    auto *t_1633 = buffer.data(target + 1633);
    auto *t_1634 = buffer.data(target + 1634);
    auto *t_1635 = buffer.data(target + 1635);
    auto *t_1636 = buffer.data(target + 1636);
    auto *t_1637 = buffer.data(target + 1637);
    auto *t_1638 = buffer.data(target + 1638);
    auto *t_1639 = buffer.data(target + 1639);
    auto *t_1640 = buffer.data(target + 1640);
    auto *t_1641 = buffer.data(target + 1641);
    auto *t_1642 = buffer.data(target + 1642);
    auto *t_1643 = buffer.data(target + 1643);
    auto *t_1644 = buffer.data(target + 1644);
    auto *t_1645 = buffer.data(target + 1645);
    auto *t_1646 = buffer.data(target + 1646);
    auto *t_1647 = buffer.data(target + 1647);
    auto *t_1648 = buffer.data(target + 1648);
    auto *t_1649 = buffer.data(target + 1649);
    auto *t_1650 = buffer.data(target + 1650);
    auto *t_1651 = buffer.data(target + 1651);
    auto *t_1652 = buffer.data(target + 1652);
    auto *t_1653 = buffer.data(target + 1653);
    auto *t_1654 = buffer.data(target + 1654);
    auto *t_1655 = buffer.data(target + 1655);
    auto *t_1656 = buffer.data(target + 1656);
    auto *t_1657 = buffer.data(target + 1657);
    auto *t_1658 = buffer.data(target + 1658);
    auto *t_1659 = buffer.data(target + 1659);
    auto *t_1660 = buffer.data(target + 1660);
    auto *t_1661 = buffer.data(target + 1661);
    auto *t_1662 = buffer.data(target + 1662);
    auto *t_1663 = buffer.data(target + 1663);
    auto *t_1664 = buffer.data(target + 1664);
    auto *t_1665 = buffer.data(target + 1665);
    auto *t_1666 = buffer.data(target + 1666);
    auto *t_1667 = buffer.data(target + 1667);
    auto *t_1668 = buffer.data(target + 1668);
    auto *t_1669 = buffer.data(target + 1669);
    auto *t_1670 = buffer.data(target + 1670);
    auto *t_1671 = buffer.data(target + 1671);
    auto *t_1672 = buffer.data(target + 1672);
    auto *t_1673 = buffer.data(target + 1673);
    auto *t_1674 = buffer.data(target + 1674);
    auto *t_1675 = buffer.data(target + 1675);
    auto *t_1676 = buffer.data(target + 1676);
    auto *t_1677 = buffer.data(target + 1677);
    auto *t_1678 = buffer.data(target + 1678);
    auto *t_1679 = buffer.data(target + 1679);
    auto *t_1680 = buffer.data(target + 1680);
    auto *t_1681 = buffer.data(target + 1681);
    auto *t_1682 = buffer.data(target + 1682);
    auto *t_1683 = buffer.data(target + 1683);
    auto *t_1684 = buffer.data(target + 1684);
    auto *t_1685 = buffer.data(target + 1685);
    auto *t_1686 = buffer.data(target + 1686);
    auto *t_1687 = buffer.data(target + 1687);
    auto *t_1688 = buffer.data(target + 1688);
    auto *t_1689 = buffer.data(target + 1689);
    auto *t_1690 = buffer.data(target + 1690);
    auto *t_1691 = buffer.data(target + 1691);
    auto *t_1692 = buffer.data(target + 1692);
    auto *t_1693 = buffer.data(target + 1693);
    auto *t_1694 = buffer.data(target + 1694);
    auto *t_1695 = buffer.data(target + 1695);
    auto *t_1696 = buffer.data(target + 1696);
    auto *t_1697 = buffer.data(target + 1697);
    auto *t_1698 = buffer.data(target + 1698);
    auto *t_1699 = buffer.data(target + 1699);
    auto *t_1700 = buffer.data(target + 1700);
    auto *t_1701 = buffer.data(target + 1701);
    auto *t_1702 = buffer.data(target + 1702);
    auto *t_1703 = buffer.data(target + 1703);
    auto *t_1704 = buffer.data(target + 1704);
    auto *t_1705 = buffer.data(target + 1705);
    auto *t_1706 = buffer.data(target + 1706);
    auto *t_1707 = buffer.data(target + 1707);
    auto *t_1708 = buffer.data(target + 1708);
    auto *t_1709 = buffer.data(target + 1709);
    auto *t_1710 = buffer.data(target + 1710);
    auto *t_1711 = buffer.data(target + 1711);
    auto *t_1712 = buffer.data(target + 1712);
    auto *t_1713 = buffer.data(target + 1713);
    auto *t_1714 = buffer.data(target + 1714);
    auto *t_1715 = buffer.data(target + 1715);
    auto *t_1716 = buffer.data(target + 1716);
    auto *t_1717 = buffer.data(target + 1717);
    auto *t_1718 = buffer.data(target + 1718);

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_1601 = buffer.data(nsi0 + 1601);
    const auto *nsi0_1602 = buffer.data(nsi0 + 1602);
    const auto *nsi0_1605 = buffer.data(nsi0 + 1605);
    const auto *nsi0_1606 = buffer.data(nsi0 + 1606);
    const auto *nsi0_1608 = buffer.data(nsi0 + 1608);
    const auto *nsi0_1610 = buffer.data(nsi0 + 1610);
    const auto *nsi0_1617 = buffer.data(nsi0 + 1617);
    const auto *nsi0_1619 = buffer.data(nsi0 + 1619);
    const auto *nsi0_1620 = buffer.data(nsi0 + 1620);
    const auto *nsi0_1621 = buffer.data(nsi0 + 1621);
    const auto *nsi0_1623 = buffer.data(nsi0 + 1623);
    const auto *nsi0_1624 = buffer.data(nsi0 + 1624);
    const auto *nsi0_1627 = buffer.data(nsi0 + 1627);
    const auto *nsi0_1629 = buffer.data(nsi0 + 1629);
    const auto *nsi0_1630 = buffer.data(nsi0 + 1630);
    const auto *nsi0_1633 = buffer.data(nsi0 + 1633);
    const auto *nsi0_1634 = buffer.data(nsi0 + 1634);
    const auto *nsi0_1636 = buffer.data(nsi0 + 1636);
    const auto *nsi0_1638 = buffer.data(nsi0 + 1638);
    const auto *nsi0_1645 = buffer.data(nsi0 + 1645);
    const auto *nsi0_1647 = buffer.data(nsi0 + 1647);
    const auto *nsi0_1648 = buffer.data(nsi0 + 1648);
    const auto *nsi0_1649 = buffer.data(nsi0 + 1649);
    const auto *nsi0_1651 = buffer.data(nsi0 + 1651);
    const auto *nsi0_1652 = buffer.data(nsi0 + 1652);
    const auto *nsi0_1655 = buffer.data(nsi0 + 1655);
    const auto *nsi0_1657 = buffer.data(nsi0 + 1657);
    const auto *nsi0_1658 = buffer.data(nsi0 + 1658);
    const auto *nsi0_1661 = buffer.data(nsi0 + 1661);
    const auto *nsi0_1662 = buffer.data(nsi0 + 1662);
    const auto *nsi0_1664 = buffer.data(nsi0 + 1664);
    const auto *nsi0_1666 = buffer.data(nsi0 + 1666);
    const auto *nsi0_1673 = buffer.data(nsi0 + 1673);
    const auto *nsi0_1675 = buffer.data(nsi0 + 1675);
    const auto *nsi0_1676 = buffer.data(nsi0 + 1676);
    const auto *nsi0_1677 = buffer.data(nsi0 + 1677);
    const auto *nsi0_1679 = buffer.data(nsi0 + 1679);
    const auto *nsi0_1680 = buffer.data(nsi0 + 1680);
    const auto *nsi0_1683 = buffer.data(nsi0 + 1683);
    const auto *nsi0_1685 = buffer.data(nsi0 + 1685);
    const auto *nsi0_1686 = buffer.data(nsi0 + 1686);
    const auto *nsi0_1689 = buffer.data(nsi0 + 1689);
    const auto *nsi0_1690 = buffer.data(nsi0 + 1690);
    const auto *nsi0_1692 = buffer.data(nsi0 + 1692);
    const auto *nsi0_1694 = buffer.data(nsi0 + 1694);
    const auto *nsi0_1701 = buffer.data(nsi0 + 1701);
    const auto *nsi0_1703 = buffer.data(nsi0 + 1703);
    const auto *nsi0_1704 = buffer.data(nsi0 + 1704);
    const auto *nsi0_1705 = buffer.data(nsi0 + 1705);
    const auto *nsi0_1707 = buffer.data(nsi0 + 1707);
    const auto *nsi0_1708 = buffer.data(nsi0 + 1708);
    const auto *nsi0_1711 = buffer.data(nsi0 + 1711);
    const auto *nsi0_1713 = buffer.data(nsi0 + 1713);
    const auto *nsi0_1714 = buffer.data(nsi0 + 1714);
    const auto *nsi0_1717 = buffer.data(nsi0 + 1717);
    const auto *nsi0_1718 = buffer.data(nsi0 + 1718);

    const auto *nsh_969 = buffer.data(nsh + 969);
    const auto *nsh_972 = buffer.data(nsh + 972);
    const auto *nsh_981 = buffer.data(nsh + 981);
    const auto *nsh_987 = buffer.data(nsh + 987);
    const auto *nsh_990 = buffer.data(nsh + 990);
    const auto *nsh_992 = buffer.data(nsh + 992);
    const auto *nsh_993 = buffer.data(nsh + 993);
    const auto *nsh_996 = buffer.data(nsh + 996);
    const auto *nsh_1002 = buffer.data(nsh + 1002);
    const auto *nsh_1007 = buffer.data(nsh + 1007);
    const auto *nsh_1008 = buffer.data(nsh + 1008);
    const auto *nsh_1010 = buffer.data(nsh + 1010);
    const auto *nsh_1011 = buffer.data(nsh + 1011);
    const auto *nsh_1013 = buffer.data(nsh + 1013);
    const auto *nsh_1014 = buffer.data(nsh + 1014);
    const auto *nsh_1017 = buffer.data(nsh + 1017);
    const auto *nsh_1023 = buffer.data(nsh + 1023);
    const auto *nsh_1028 = buffer.data(nsh + 1028);
    const auto *nsh_1029 = buffer.data(nsh + 1029);
    const auto *nsh_1031 = buffer.data(nsh + 1031);
    const auto *nsh_1032 = buffer.data(nsh + 1032);
    const auto *nsh_1034 = buffer.data(nsh + 1034);
    const auto *nsh_1035 = buffer.data(nsh + 1035);
    const auto *nsh_1038 = buffer.data(nsh + 1038);
    const auto *nsh_1044 = buffer.data(nsh + 1044);
    const auto *nsh_1049 = buffer.data(nsh + 1049);
    const auto *nsh_1050 = buffer.data(nsh + 1050);
    const auto *nsh_1052 = buffer.data(nsh + 1052);
    const auto *nsh_1053 = buffer.data(nsh + 1053);
    const auto *nsh_1055 = buffer.data(nsh + 1055);
    const auto *nsh_1059 = buffer.data(nsh + 1059);
    const auto *nsh_1070 = buffer.data(nsh + 1070);
    const auto *nsh_1071 = buffer.data(nsh + 1071);
    const auto *nsh_1073 = buffer.data(nsh + 1073);
    const auto *nsh_1076 = buffer.data(nsh + 1076);
    const auto *nsh_1202 = buffer.data(nsh + 1202);
    const auto *nsh_1203 = buffer.data(nsh + 1203);
    const auto *nsh_1206 = buffer.data(nsh + 1206);
    const auto *nsh_1207 = buffer.data(nsh + 1207);
    const auto *nsh_1209 = buffer.data(nsh + 1209);
    const auto *nsh_1211 = buffer.data(nsh + 1211);
    const auto *nsh_1212 = buffer.data(nsh + 1212);
    const auto *nsh_1213 = buffer.data(nsh + 1213);
    const auto *nsh_1214 = buffer.data(nsh + 1214);
    const auto *nsh_1215 = buffer.data(nsh + 1215);
    const auto *nsh_1216 = buffer.data(nsh + 1216);
    const auto *nsh_1217 = buffer.data(nsh + 1217);
    const auto *nsh_1218 = buffer.data(nsh + 1218);
    const auto *nsh_1221 = buffer.data(nsh + 1221);
    const auto *nsh_1223 = buffer.data(nsh + 1223);
    const auto *nsh_1224 = buffer.data(nsh + 1224);
    const auto *nsh_1227 = buffer.data(nsh + 1227);
    const auto *nsh_1228 = buffer.data(nsh + 1228);
    const auto *nsh_1230 = buffer.data(nsh + 1230);
    const auto *nsh_1232 = buffer.data(nsh + 1232);
    const auto *nsh_1233 = buffer.data(nsh + 1233);
    const auto *nsh_1234 = buffer.data(nsh + 1234);
    const auto *nsh_1235 = buffer.data(nsh + 1235);
    const auto *nsh_1236 = buffer.data(nsh + 1236);
    const auto *nsh_1237 = buffer.data(nsh + 1237);
    const auto *nsh_1238 = buffer.data(nsh + 1238);
    const auto *nsh_1239 = buffer.data(nsh + 1239);
    const auto *nsh_1242 = buffer.data(nsh + 1242);
    const auto *nsh_1244 = buffer.data(nsh + 1244);
    const auto *nsh_1245 = buffer.data(nsh + 1245);
    const auto *nsh_1248 = buffer.data(nsh + 1248);
    const auto *nsh_1249 = buffer.data(nsh + 1249);
    const auto *nsh_1251 = buffer.data(nsh + 1251);
    const auto *nsh_1253 = buffer.data(nsh + 1253);
    const auto *nsh_1254 = buffer.data(nsh + 1254);
    const auto *nsh_1255 = buffer.data(nsh + 1255);
    const auto *nsh_1256 = buffer.data(nsh + 1256);
    const auto *nsh_1257 = buffer.data(nsh + 1257);
    const auto *nsh_1258 = buffer.data(nsh + 1258);
    const auto *nsh_1259 = buffer.data(nsh + 1259);
    const auto *nsh_1260 = buffer.data(nsh + 1260);
    const auto *nsh_1263 = buffer.data(nsh + 1263);
    const auto *nsh_1265 = buffer.data(nsh + 1265);
    const auto *nsh_1266 = buffer.data(nsh + 1266);
    const auto *nsh_1269 = buffer.data(nsh + 1269);
    const auto *nsh_1270 = buffer.data(nsh + 1270);
    const auto *nsh_1272 = buffer.data(nsh + 1272);
    const auto *nsh_1274 = buffer.data(nsh + 1274);
    const auto *nsh_1275 = buffer.data(nsh + 1275);
    const auto *nsh_1276 = buffer.data(nsh + 1276);
    const auto *nsh_1277 = buffer.data(nsh + 1277);
    const auto *nsh_1278 = buffer.data(nsh + 1278);
    const auto *nsh_1279 = buffer.data(nsh + 1279);
    const auto *nsh_1280 = buffer.data(nsh + 1280);
    const auto *nsh_1281 = buffer.data(nsh + 1281);
    const auto *nsh_1284 = buffer.data(nsh + 1284);
    const auto *nsh_1286 = buffer.data(nsh + 1286);
    const auto *nsh_1287 = buffer.data(nsh + 1287);
    const auto *nsh_1290 = buffer.data(nsh + 1290);
    const auto *nsh_1291 = buffer.data(nsh + 1291);

    const auto *nsi1_1601 = buffer.data(nsi1 + 1601);
    const auto *nsi1_1602 = buffer.data(nsi1 + 1602);
    const auto *nsi1_1605 = buffer.data(nsi1 + 1605);
    const auto *nsi1_1606 = buffer.data(nsi1 + 1606);
    const auto *nsi1_1608 = buffer.data(nsi1 + 1608);
    const auto *nsi1_1610 = buffer.data(nsi1 + 1610);
    const auto *nsi1_1617 = buffer.data(nsi1 + 1617);
    const auto *nsi1_1619 = buffer.data(nsi1 + 1619);
    const auto *nsi1_1620 = buffer.data(nsi1 + 1620);
    const auto *nsi1_1621 = buffer.data(nsi1 + 1621);
    const auto *nsi1_1623 = buffer.data(nsi1 + 1623);
    const auto *nsi1_1624 = buffer.data(nsi1 + 1624);
    const auto *nsi1_1627 = buffer.data(nsi1 + 1627);
    const auto *nsi1_1629 = buffer.data(nsi1 + 1629);
    const auto *nsi1_1630 = buffer.data(nsi1 + 1630);
    const auto *nsi1_1633 = buffer.data(nsi1 + 1633);
    const auto *nsi1_1634 = buffer.data(nsi1 + 1634);
    const auto *nsi1_1636 = buffer.data(nsi1 + 1636);
    const auto *nsi1_1638 = buffer.data(nsi1 + 1638);
    const auto *nsi1_1645 = buffer.data(nsi1 + 1645);
    const auto *nsi1_1647 = buffer.data(nsi1 + 1647);
    const auto *nsi1_1648 = buffer.data(nsi1 + 1648);
    const auto *nsi1_1649 = buffer.data(nsi1 + 1649);
    const auto *nsi1_1651 = buffer.data(nsi1 + 1651);
    const auto *nsi1_1652 = buffer.data(nsi1 + 1652);
    const auto *nsi1_1655 = buffer.data(nsi1 + 1655);
    const auto *nsi1_1657 = buffer.data(nsi1 + 1657);
    const auto *nsi1_1658 = buffer.data(nsi1 + 1658);
    const auto *nsi1_1661 = buffer.data(nsi1 + 1661);
    const auto *nsi1_1662 = buffer.data(nsi1 + 1662);
    const auto *nsi1_1664 = buffer.data(nsi1 + 1664);
    const auto *nsi1_1666 = buffer.data(nsi1 + 1666);
    const auto *nsi1_1673 = buffer.data(nsi1 + 1673);
    const auto *nsi1_1675 = buffer.data(nsi1 + 1675);
    const auto *nsi1_1676 = buffer.data(nsi1 + 1676);
    const auto *nsi1_1677 = buffer.data(nsi1 + 1677);
    const auto *nsi1_1679 = buffer.data(nsi1 + 1679);
    const auto *nsi1_1680 = buffer.data(nsi1 + 1680);
    const auto *nsi1_1683 = buffer.data(nsi1 + 1683);
    const auto *nsi1_1685 = buffer.data(nsi1 + 1685);
    const auto *nsi1_1686 = buffer.data(nsi1 + 1686);
    const auto *nsi1_1689 = buffer.data(nsi1 + 1689);
    const auto *nsi1_1690 = buffer.data(nsi1 + 1690);
    const auto *nsi1_1692 = buffer.data(nsi1 + 1692);
    const auto *nsi1_1694 = buffer.data(nsi1 + 1694);
    const auto *nsi1_1701 = buffer.data(nsi1 + 1701);
    const auto *nsi1_1703 = buffer.data(nsi1 + 1703);
    const auto *nsi1_1704 = buffer.data(nsi1 + 1704);
    const auto *nsi1_1705 = buffer.data(nsi1 + 1705);
    const auto *nsi1_1707 = buffer.data(nsi1 + 1707);
    const auto *nsi1_1708 = buffer.data(nsi1 + 1708);
    const auto *nsi1_1711 = buffer.data(nsi1 + 1711);
    const auto *nsi1_1713 = buffer.data(nsi1 + 1713);
    const auto *nsi1_1714 = buffer.data(nsi1 + 1714);
    const auto *nsi1_1717 = buffer.data(nsi1 + 1717);
    const auto *nsi1_1718 = buffer.data(nsi1 + 1718);

    const auto *osh_1200 = buffer.data(osh + 1200);
    const auto *osh_1202 = buffer.data(osh + 1202);
    const auto *osh_1203 = buffer.data(osh + 1203);
    const auto *osh_1206 = buffer.data(osh + 1206);
    const auto *osh_1212 = buffer.data(osh + 1212);
    const auto *osh_1213 = buffer.data(osh + 1213);
    const auto *osh_1214 = buffer.data(osh + 1214);
    const auto *osh_1215 = buffer.data(osh + 1215);
    const auto *osh_1216 = buffer.data(osh + 1216);
    const auto *osh_1217 = buffer.data(osh + 1217);
    const auto *osh_1218 = buffer.data(osh + 1218);
    const auto *osh_1220 = buffer.data(osh + 1220);
    const auto *osh_1221 = buffer.data(osh + 1221);
    const auto *osh_1223 = buffer.data(osh + 1223);
    const auto *osh_1224 = buffer.data(osh + 1224);
    const auto *osh_1227 = buffer.data(osh + 1227);
    const auto *osh_1233 = buffer.data(osh + 1233);
    const auto *osh_1234 = buffer.data(osh + 1234);
    const auto *osh_1235 = buffer.data(osh + 1235);
    const auto *osh_1236 = buffer.data(osh + 1236);
    const auto *osh_1237 = buffer.data(osh + 1237);
    const auto *osh_1238 = buffer.data(osh + 1238);
    const auto *osh_1239 = buffer.data(osh + 1239);
    const auto *osh_1241 = buffer.data(osh + 1241);
    const auto *osh_1242 = buffer.data(osh + 1242);
    const auto *osh_1244 = buffer.data(osh + 1244);
    const auto *osh_1245 = buffer.data(osh + 1245);
    const auto *osh_1248 = buffer.data(osh + 1248);
    const auto *osh_1254 = buffer.data(osh + 1254);
    const auto *osh_1255 = buffer.data(osh + 1255);
    const auto *osh_1256 = buffer.data(osh + 1256);
    const auto *osh_1257 = buffer.data(osh + 1257);
    const auto *osh_1258 = buffer.data(osh + 1258);
    const auto *osh_1259 = buffer.data(osh + 1259);
    const auto *osh_1260 = buffer.data(osh + 1260);
    const auto *osh_1262 = buffer.data(osh + 1262);
    const auto *osh_1263 = buffer.data(osh + 1263);
    const auto *osh_1265 = buffer.data(osh + 1265);
    const auto *osh_1266 = buffer.data(osh + 1266);
    const auto *osh_1269 = buffer.data(osh + 1269);
    const auto *osh_1275 = buffer.data(osh + 1275);
    const auto *osh_1276 = buffer.data(osh + 1276);
    const auto *osh_1277 = buffer.data(osh + 1277);
    const auto *osh_1278 = buffer.data(osh + 1278);
    const auto *osh_1279 = buffer.data(osh + 1279);
    const auto *osh_1280 = buffer.data(osh + 1280);
    const auto *osh_1281 = buffer.data(osh + 1281);
    const auto *osh_1283 = buffer.data(osh + 1283);
    const auto *osh_1284 = buffer.data(osh + 1284);
    const auto *osh_1286 = buffer.data(osh + 1286);

#pragma omp simd aligned(t_1601, t_1602, t_1603, pa_x, pc_x, pc_z, nsi0_1601, nsi0_1602, \
                         nsh_969, nsh_1202, nsh_1203, nsi1_1601, nsi1_1602, \
                         osh_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = pa_x[k] * nsi0_1601[k]
                    + f_14 * nsh_1202[k]
                    - f_10 * pc_x[k] * nsi1_1601[k];

        t_1602[k] = pa_x[k] * nsi0_1602[k]
                    + f_13 * nsh_1203[k]
                    - f_10 * pc_x[k] * nsi1_1602[k];

        t_1603[k] = f_12 * nsh_969[k]
                    + f_3 * pc_z[k] * osh_1200[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, pa_x, pc_x, pc_y, nsi0_1605, nsi0_1606, \
                         nsh_992, nsh_1206, nsh_1207, nsi1_1605, nsi1_1606, \
                         osh_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = f_19 * nsh_992[k]
                    + f_3 * pc_y[k] * osh_1202[k];

        t_1605[k] = pa_x[k] * nsi0_1605[k]
                    + f_13 * nsh_1206[k]
                    - f_10 * pc_x[k] * nsi1_1605[k];

        t_1606[k] = pa_x[k] * nsi0_1606[k]
                    + f_12 * nsh_1207[k]
                    - f_10 * pc_x[k] * nsi1_1606[k];
    }

#pragma omp simd aligned(t_1607, t_1608, t_1609, pa_x, pc_x, pc_y, pc_z, nsi0_1608, nsh_972, \
                         nsh_996, nsh_1209, nsi1_1608, osh_1203, \
                         osh_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1607[k] = f_12 * nsh_972[k]
                    + f_3 * pc_z[k] * osh_1203[k];

        t_1608[k] = pa_x[k] * nsi0_1608[k]
                    + f_12 * nsh_1209[k]
                    - f_10 * pc_x[k] * nsi1_1608[k];

        t_1609[k] = f_19 * nsh_996[k]
                    + f_3 * pc_y[k] * osh_1206[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, pa_x, pc_x, nsi0_1610, nsh_1211, \
                         nsh_1212, nsh_1213, nsh_1214, nsi1_1610, osh_1212, osh_1213, \
                         osh_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = pa_x[k] * nsi0_1610[k]
                    + f_12 * nsh_1211[k]
                    - f_10 * pc_x[k] * nsi1_1610[k];

        t_1611[k] = f_11 * nsh_1212[k]
                    + f_3 * pc_x[k] * osh_1212[k];

        t_1612[k] = f_11 * nsh_1213[k]
                    + f_3 * pc_x[k] * osh_1213[k];

        t_1613[k] = f_11 * nsh_1214[k]
                    + f_3 * pc_x[k] * osh_1214[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, t_1617, pa_x, pc_x, nsi0_1617, nsh_1215, \
                         nsh_1216, nsh_1217, nsi1_1617, osh_1215, osh_1216, \
                         osh_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_11 * nsh_1215[k]
                    + f_3 * pc_x[k] * osh_1215[k];

        t_1615[k] = f_11 * nsh_1216[k]
                    + f_3 * pc_x[k] * osh_1216[k];

        t_1616[k] = f_11 * nsh_1217[k]
                    + f_3 * pc_x[k] * osh_1217[k];

        t_1617[k] = pa_x[k] * nsi0_1617[k]
                    - f_10 * pc_x[k] * nsi1_1617[k];
    }

#pragma omp simd aligned(t_1618, t_1619, t_1620, t_1621, pa_x, pc_x, pc_z, nsi0_1619, \
                         nsi0_1620, nsi0_1621, nsh_981, nsi1_1619, nsi1_1620, nsi1_1621, \
                         osh_1212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1618[k] = f_12 * nsh_981[k]
                    + f_3 * pc_z[k] * osh_1212[k];

        t_1619[k] = pa_x[k] * nsi0_1619[k]
                    - f_10 * pc_x[k] * nsi1_1619[k];

        t_1620[k] = pa_x[k] * nsi0_1620[k]
                    - f_10 * pc_x[k] * nsi1_1620[k];

        t_1621[k] = pa_x[k] * nsi0_1621[k]
                    - f_10 * pc_x[k] * nsi1_1621[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, t_1625, pa_x, pc_x, pc_y, nsi0_1623, \
                         nsi0_1624, nsh_1007, nsh_1008, nsh_1218, nsi1_1623, nsi1_1624, \
                         osh_1217, osh_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_19 * nsh_1007[k]
                    + f_3 * pc_y[k] * osh_1217[k];

        t_1623[k] = pa_x[k] * nsi0_1623[k]
                    - f_10 * pc_x[k] * nsi1_1623[k];

        t_1624[k] = pa_x[k] * nsi0_1624[k]
                    + f_21 * nsh_1218[k]
                    - f_10 * pc_x[k] * nsi1_1624[k];

        t_1625[k] = f_20 * nsh_1008[k]
                    + f_3 * pc_y[k] * osh_1218[k];
    }

#pragma omp simd aligned(t_1626, t_1627, t_1628, pa_x, pc_x, pc_y, pc_z, nsi0_1627, nsh_987, \
                         nsh_1010, nsh_1221, nsi1_1627, osh_1218, \
                         osh_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1626[k] = f_13 * nsh_987[k]
                    + f_3 * pc_z[k] * osh_1218[k];

        t_1627[k] = pa_x[k] * nsi0_1627[k]
                    + f_14 * nsh_1221[k]
                    - f_10 * pc_x[k] * nsi1_1627[k];

        t_1628[k] = f_20 * nsh_1010[k]
                    + f_3 * pc_y[k] * osh_1220[k];
    }

#pragma omp simd aligned(t_1629, t_1630, t_1631, pa_x, pc_x, pc_z, nsi0_1629, nsi0_1630, \
                         nsh_990, nsh_1223, nsh_1224, nsi1_1629, nsi1_1630, \
                         osh_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1629[k] = pa_x[k] * nsi0_1629[k]
                    + f_14 * nsh_1223[k]
                    - f_10 * pc_x[k] * nsi1_1629[k];

        t_1630[k] = pa_x[k] * nsi0_1630[k]
                    + f_13 * nsh_1224[k]
                    - f_10 * pc_x[k] * nsi1_1630[k];

        t_1631[k] = f_13 * nsh_990[k]
                    + f_3 * pc_z[k] * osh_1221[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pa_x, pc_x, pc_y, nsi0_1633, nsi0_1634, \
                         nsh_1013, nsh_1227, nsh_1228, nsi1_1633, nsi1_1634, \
                         osh_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_20 * nsh_1013[k]
                    + f_3 * pc_y[k] * osh_1223[k];

        t_1633[k] = pa_x[k] * nsi0_1633[k]
                    + f_13 * nsh_1227[k]
                    - f_10 * pc_x[k] * nsi1_1633[k];

        t_1634[k] = pa_x[k] * nsi0_1634[k]
                    + f_12 * nsh_1228[k]
                    - f_10 * pc_x[k] * nsi1_1634[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, pa_x, pc_x, pc_y, pc_z, nsi0_1636, nsh_993, \
                         nsh_1017, nsh_1230, nsi1_1636, osh_1224, \
                         osh_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = f_13 * nsh_993[k]
                    + f_3 * pc_z[k] * osh_1224[k];

        t_1636[k] = pa_x[k] * nsi0_1636[k]
                    + f_12 * nsh_1230[k]
                    - f_10 * pc_x[k] * nsi1_1636[k];

        t_1637[k] = f_20 * nsh_1017[k]
                    + f_3 * pc_y[k] * osh_1227[k];
    }

#pragma omp simd aligned(t_1638, t_1639, t_1640, t_1641, pa_x, pc_x, nsi0_1638, nsh_1232, \
                         nsh_1233, nsh_1234, nsh_1235, nsi1_1638, osh_1233, osh_1234, \
                         osh_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1638[k] = pa_x[k] * nsi0_1638[k]
                    + f_12 * nsh_1232[k]
                    - f_10 * pc_x[k] * nsi1_1638[k];

        t_1639[k] = f_11 * nsh_1233[k]
                    + f_3 * pc_x[k] * osh_1233[k];

        t_1640[k] = f_11 * nsh_1234[k]
                    + f_3 * pc_x[k] * osh_1234[k];

        t_1641[k] = f_11 * nsh_1235[k]
                    + f_3 * pc_x[k] * osh_1235[k];
    }

#pragma omp simd aligned(t_1642, t_1643, t_1644, t_1645, pa_x, pc_x, nsi0_1645, nsh_1236, \
                         nsh_1237, nsh_1238, nsi1_1645, osh_1236, osh_1237, \
                         osh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1642[k] = f_11 * nsh_1236[k]
                    + f_3 * pc_x[k] * osh_1236[k];

        t_1643[k] = f_11 * nsh_1237[k]
                    + f_3 * pc_x[k] * osh_1237[k];

        t_1644[k] = f_11 * nsh_1238[k]
                    + f_3 * pc_x[k] * osh_1238[k];

        t_1645[k] = pa_x[k] * nsi0_1645[k]
                    - f_10 * pc_x[k] * nsi1_1645[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, t_1649, pa_x, pc_x, pc_z, nsi0_1647, \
                         nsi0_1648, nsi0_1649, nsh_1002, nsi1_1647, nsi1_1648, nsi1_1649, \
                         osh_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = f_13 * nsh_1002[k]
                    + f_3 * pc_z[k] * osh_1233[k];

        t_1647[k] = pa_x[k] * nsi0_1647[k]
                    - f_10 * pc_x[k] * nsi1_1647[k];

        t_1648[k] = pa_x[k] * nsi0_1648[k]
                    - f_10 * pc_x[k] * nsi1_1648[k];

        t_1649[k] = pa_x[k] * nsi0_1649[k]
                    - f_10 * pc_x[k] * nsi1_1649[k];
    }

#pragma omp simd aligned(t_1650, t_1651, t_1652, t_1653, pa_x, pc_x, pc_y, nsi0_1651, \
                         nsi0_1652, nsh_1028, nsh_1029, nsh_1239, nsi1_1651, nsi1_1652, \
                         osh_1238, osh_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1650[k] = f_20 * nsh_1028[k]
                    + f_3 * pc_y[k] * osh_1238[k];

        t_1651[k] = pa_x[k] * nsi0_1651[k]
                    - f_10 * pc_x[k] * nsi1_1651[k];

        t_1652[k] = pa_x[k] * nsi0_1652[k]
                    + f_21 * nsh_1239[k]
                    - f_10 * pc_x[k] * nsi1_1652[k];

        t_1653[k] = f_21 * nsh_1029[k]
                    + f_3 * pc_y[k] * osh_1239[k];
    }

#pragma omp simd aligned(t_1654, t_1655, t_1656, pa_x, pc_x, pc_y, pc_z, nsi0_1655, nsh_1008, \
                         nsh_1031, nsh_1242, nsi1_1655, osh_1239, \
                         osh_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1654[k] = f_14 * nsh_1008[k]
                    + f_3 * pc_z[k] * osh_1239[k];

        t_1655[k] = pa_x[k] * nsi0_1655[k]
                    + f_14 * nsh_1242[k]
                    - f_10 * pc_x[k] * nsi1_1655[k];

        t_1656[k] = f_21 * nsh_1031[k]
                    + f_3 * pc_y[k] * osh_1241[k];
    }

#pragma omp simd aligned(t_1657, t_1658, t_1659, pa_x, pc_x, pc_z, nsi0_1657, nsi0_1658, \
                         nsh_1011, nsh_1244, nsh_1245, nsi1_1657, nsi1_1658, \
                         osh_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1657[k] = pa_x[k] * nsi0_1657[k]
                    + f_14 * nsh_1244[k]
                    - f_10 * pc_x[k] * nsi1_1657[k];

        t_1658[k] = pa_x[k] * nsi0_1658[k]
                    + f_13 * nsh_1245[k]
                    - f_10 * pc_x[k] * nsi1_1658[k];

        t_1659[k] = f_14 * nsh_1011[k]
                    + f_3 * pc_z[k] * osh_1242[k];
    }

#pragma omp simd aligned(t_1660, t_1661, t_1662, pa_x, pc_x, pc_y, nsi0_1661, nsi0_1662, \
                         nsh_1034, nsh_1248, nsh_1249, nsi1_1661, nsi1_1662, \
                         osh_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1660[k] = f_21 * nsh_1034[k]
                    + f_3 * pc_y[k] * osh_1244[k];

        t_1661[k] = pa_x[k] * nsi0_1661[k]
                    + f_13 * nsh_1248[k]
                    - f_10 * pc_x[k] * nsi1_1661[k];

        t_1662[k] = pa_x[k] * nsi0_1662[k]
                    + f_12 * nsh_1249[k]
                    - f_10 * pc_x[k] * nsi1_1662[k];
    }

#pragma omp simd aligned(t_1663, t_1664, t_1665, pa_x, pc_x, pc_y, pc_z, nsi0_1664, nsh_1014, \
                         nsh_1038, nsh_1251, nsi1_1664, osh_1245, \
                         osh_1248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1663[k] = f_14 * nsh_1014[k]
                    + f_3 * pc_z[k] * osh_1245[k];

        t_1664[k] = pa_x[k] * nsi0_1664[k]
                    + f_12 * nsh_1251[k]
                    - f_10 * pc_x[k] * nsi1_1664[k];

        t_1665[k] = f_21 * nsh_1038[k]
                    + f_3 * pc_y[k] * osh_1248[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, t_1669, pa_x, pc_x, nsi0_1666, nsh_1253, \
                         nsh_1254, nsh_1255, nsh_1256, nsi1_1666, osh_1254, osh_1255, \
                         osh_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = pa_x[k] * nsi0_1666[k]
                    + f_12 * nsh_1253[k]
                    - f_10 * pc_x[k] * nsi1_1666[k];

        t_1667[k] = f_11 * nsh_1254[k]
                    + f_3 * pc_x[k] * osh_1254[k];

        t_1668[k] = f_11 * nsh_1255[k]
                    + f_3 * pc_x[k] * osh_1255[k];

        t_1669[k] = f_11 * nsh_1256[k]
                    + f_3 * pc_x[k] * osh_1256[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, t_1673, pa_x, pc_x, nsi0_1673, nsh_1257, \
                         nsh_1258, nsh_1259, nsi1_1673, osh_1257, osh_1258, \
                         osh_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_11 * nsh_1257[k]
                    + f_3 * pc_x[k] * osh_1257[k];

        t_1671[k] = f_11 * nsh_1258[k]
                    + f_3 * pc_x[k] * osh_1258[k];

        t_1672[k] = f_11 * nsh_1259[k]
                    + f_3 * pc_x[k] * osh_1259[k];

        t_1673[k] = pa_x[k] * nsi0_1673[k]
                    - f_10 * pc_x[k] * nsi1_1673[k];
    }

#pragma omp simd aligned(t_1674, t_1675, t_1676, t_1677, pa_x, pc_x, pc_z, nsi0_1675, \
                         nsi0_1676, nsi0_1677, nsh_1023, nsi1_1675, nsi1_1676, nsi1_1677, \
                         osh_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1674[k] = f_14 * nsh_1023[k]
                    + f_3 * pc_z[k] * osh_1254[k];

        t_1675[k] = pa_x[k] * nsi0_1675[k]
                    - f_10 * pc_x[k] * nsi1_1675[k];

        t_1676[k] = pa_x[k] * nsi0_1676[k]
                    - f_10 * pc_x[k] * nsi1_1676[k];

        t_1677[k] = pa_x[k] * nsi0_1677[k]
                    - f_10 * pc_x[k] * nsi1_1677[k];
    }

#pragma omp simd aligned(t_1678, t_1679, t_1680, t_1681, pa_x, pc_x, pc_y, nsi0_1679, \
                         nsi0_1680, nsh_1049, nsh_1050, nsh_1260, nsi1_1679, nsi1_1680, \
                         osh_1259, osh_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1678[k] = f_21 * nsh_1049[k]
                    + f_3 * pc_y[k] * osh_1259[k];

        t_1679[k] = pa_x[k] * nsi0_1679[k]
                    - f_10 * pc_x[k] * nsi1_1679[k];

        t_1680[k] = pa_x[k] * nsi0_1680[k]
                    + f_21 * nsh_1260[k]
                    - f_10 * pc_x[k] * nsi1_1680[k];

        t_1681[k] = f_22 * nsh_1050[k]
                    + f_3 * pc_y[k] * osh_1260[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, pa_x, pc_x, pc_y, pc_z, nsi0_1683, nsh_1029, \
                         nsh_1052, nsh_1263, nsi1_1683, osh_1260, \
                         osh_1262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_22 * nsh_1029[k]
                    + f_3 * pc_z[k] * osh_1260[k];

        t_1683[k] = pa_x[k] * nsi0_1683[k]
                    + f_14 * nsh_1263[k]
                    - f_10 * pc_x[k] * nsi1_1683[k];

        t_1684[k] = f_22 * nsh_1052[k]
                    + f_3 * pc_y[k] * osh_1262[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pa_x, pc_x, pc_z, nsi0_1685, nsi0_1686, \
                         nsh_1032, nsh_1265, nsh_1266, nsi1_1685, nsi1_1686, \
                         osh_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = pa_x[k] * nsi0_1685[k]
                    + f_14 * nsh_1265[k]
                    - f_10 * pc_x[k] * nsi1_1685[k];

        t_1686[k] = pa_x[k] * nsi0_1686[k]
                    + f_13 * nsh_1266[k]
                    - f_10 * pc_x[k] * nsi1_1686[k];

        t_1687[k] = f_22 * nsh_1032[k]
                    + f_3 * pc_z[k] * osh_1263[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, pa_x, pc_x, pc_y, nsi0_1689, nsi0_1690, \
                         nsh_1055, nsh_1269, nsh_1270, nsi1_1689, nsi1_1690, \
                         osh_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = f_22 * nsh_1055[k]
                    + f_3 * pc_y[k] * osh_1265[k];

        t_1689[k] = pa_x[k] * nsi0_1689[k]
                    + f_13 * nsh_1269[k]
                    - f_10 * pc_x[k] * nsi1_1689[k];

        t_1690[k] = pa_x[k] * nsi0_1690[k]
                    + f_12 * nsh_1270[k]
                    - f_10 * pc_x[k] * nsi1_1690[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, pa_x, pc_x, pc_y, pc_z, nsi0_1692, nsh_1035, \
                         nsh_1059, nsh_1272, nsi1_1692, osh_1266, \
                         osh_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_22 * nsh_1035[k]
                    + f_3 * pc_z[k] * osh_1266[k];

        t_1692[k] = pa_x[k] * nsi0_1692[k]
                    + f_12 * nsh_1272[k]
                    - f_10 * pc_x[k] * nsi1_1692[k];

        t_1693[k] = f_22 * nsh_1059[k]
                    + f_3 * pc_y[k] * osh_1269[k];
    }

#pragma omp simd aligned(t_1694, t_1695, t_1696, t_1697, pa_x, pc_x, nsi0_1694, nsh_1274, \
                         nsh_1275, nsh_1276, nsh_1277, nsi1_1694, osh_1275, osh_1276, \
                         osh_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1694[k] = pa_x[k] * nsi0_1694[k]
                    + f_12 * nsh_1274[k]
                    - f_10 * pc_x[k] * nsi1_1694[k];

        t_1695[k] = f_11 * nsh_1275[k]
                    + f_3 * pc_x[k] * osh_1275[k];

        t_1696[k] = f_11 * nsh_1276[k]
                    + f_3 * pc_x[k] * osh_1276[k];

        t_1697[k] = f_11 * nsh_1277[k]
                    + f_3 * pc_x[k] * osh_1277[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, t_1701, pa_x, pc_x, nsi0_1701, nsh_1278, \
                         nsh_1279, nsh_1280, nsi1_1701, osh_1278, osh_1279, \
                         osh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_11 * nsh_1278[k]
                    + f_3 * pc_x[k] * osh_1278[k];

        t_1699[k] = f_11 * nsh_1279[k]
                    + f_3 * pc_x[k] * osh_1279[k];

        t_1700[k] = f_11 * nsh_1280[k]
                    + f_3 * pc_x[k] * osh_1280[k];

        t_1701[k] = pa_x[k] * nsi0_1701[k]
                    - f_10 * pc_x[k] * nsi1_1701[k];
    }

#pragma omp simd aligned(t_1702, t_1703, t_1704, t_1705, pa_x, pc_x, pc_z, nsi0_1703, \
                         nsi0_1704, nsi0_1705, nsh_1044, nsi1_1703, nsi1_1704, nsi1_1705, \
                         osh_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1702[k] = f_22 * nsh_1044[k]
                    + f_3 * pc_z[k] * osh_1275[k];

        t_1703[k] = pa_x[k] * nsi0_1703[k]
                    - f_10 * pc_x[k] * nsi1_1703[k];

        t_1704[k] = pa_x[k] * nsi0_1704[k]
                    - f_10 * pc_x[k] * nsi1_1704[k];

        t_1705[k] = pa_x[k] * nsi0_1705[k]
                    - f_10 * pc_x[k] * nsi1_1705[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, t_1709, pa_x, pc_x, pc_y, nsi0_1707, \
                         nsi0_1708, nsh_1070, nsh_1071, nsh_1281, nsi1_1707, nsi1_1708, \
                         osh_1280, osh_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = f_22 * nsh_1070[k]
                    + f_3 * pc_y[k] * osh_1280[k];

        t_1707[k] = pa_x[k] * nsi0_1707[k]
                    - f_10 * pc_x[k] * nsi1_1707[k];

        t_1708[k] = pa_x[k] * nsi0_1708[k]
                    + f_21 * nsh_1281[k]
                    - f_10 * pc_x[k] * nsi1_1708[k];

        t_1709[k] = f_14 * nsh_1071[k]
                    + f_3 * pc_y[k] * osh_1281[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, pa_x, pc_x, pc_y, pc_z, nsi0_1711, nsh_1050, \
                         nsh_1073, nsh_1284, nsi1_1711, osh_1281, \
                         osh_1283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = f_21 * nsh_1050[k]
                    + f_3 * pc_z[k] * osh_1281[k];

        t_1711[k] = pa_x[k] * nsi0_1711[k]
                    + f_14 * nsh_1284[k]
                    - f_10 * pc_x[k] * nsi1_1711[k];

        t_1712[k] = f_14 * nsh_1073[k]
                    + f_3 * pc_y[k] * osh_1283[k];
    }

#pragma omp simd aligned(t_1713, t_1714, t_1715, pa_x, pc_x, pc_z, nsi0_1713, nsi0_1714, \
                         nsh_1053, nsh_1286, nsh_1287, nsi1_1713, nsi1_1714, \
                         osh_1284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1713[k] = pa_x[k] * nsi0_1713[k]
                    + f_14 * nsh_1286[k]
                    - f_10 * pc_x[k] * nsi1_1713[k];

        t_1714[k] = pa_x[k] * nsi0_1714[k]
                    + f_13 * nsh_1287[k]
                    - f_10 * pc_x[k] * nsi1_1714[k];

        t_1715[k] = f_21 * nsh_1053[k]
                    + f_3 * pc_z[k] * osh_1284[k];
    }

#pragma omp simd aligned(t_1716, t_1717, t_1718, pa_x, pc_x, pc_y, nsi0_1717, nsi0_1718, \
                         nsh_1076, nsh_1290, nsh_1291, nsi1_1717, nsi1_1718, \
                         osh_1286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1716[k] = f_14 * nsh_1076[k]
                    + f_3 * pc_y[k] * osh_1286[k];

        t_1717[k] = pa_x[k] * nsi0_1717[k]
                    + f_13 * nsh_1290[k]
                    - f_10 * pc_x[k] * nsi1_1717[k];

        t_1718[k] = pa_x[k] * nsi0_1718[k]
                    + f_12 * nsh_1291[k]
                    - f_10 * pc_x[k] * nsi1_1718[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osg0, const size_t osg1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.0 / q;
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

    auto *t_1719 = buffer.data(target + 1719);
    auto *t_1720 = buffer.data(target + 1720);
    auto *t_1721 = buffer.data(target + 1721);
    auto *t_1722 = buffer.data(target + 1722);
    auto *t_1723 = buffer.data(target + 1723);
    auto *t_1724 = buffer.data(target + 1724);
    auto *t_1725 = buffer.data(target + 1725);
    auto *t_1726 = buffer.data(target + 1726);
    auto *t_1727 = buffer.data(target + 1727);
    auto *t_1728 = buffer.data(target + 1728);
    auto *t_1729 = buffer.data(target + 1729);
    auto *t_1730 = buffer.data(target + 1730);
    auto *t_1731 = buffer.data(target + 1731);
    auto *t_1732 = buffer.data(target + 1732);
    auto *t_1733 = buffer.data(target + 1733);
    auto *t_1734 = buffer.data(target + 1734);
    auto *t_1735 = buffer.data(target + 1735);
    auto *t_1736 = buffer.data(target + 1736);
    auto *t_1737 = buffer.data(target + 1737);
    auto *t_1738 = buffer.data(target + 1738);
    auto *t_1739 = buffer.data(target + 1739);
    auto *t_1740 = buffer.data(target + 1740);
    auto *t_1741 = buffer.data(target + 1741);
    auto *t_1742 = buffer.data(target + 1742);
    auto *t_1743 = buffer.data(target + 1743);
    auto *t_1744 = buffer.data(target + 1744);
    auto *t_1745 = buffer.data(target + 1745);
    auto *t_1746 = buffer.data(target + 1746);
    auto *t_1747 = buffer.data(target + 1747);
    auto *t_1748 = buffer.data(target + 1748);
    auto *t_1749 = buffer.data(target + 1749);
    auto *t_1750 = buffer.data(target + 1750);
    auto *t_1751 = buffer.data(target + 1751);
    auto *t_1752 = buffer.data(target + 1752);
    auto *t_1753 = buffer.data(target + 1753);
    auto *t_1754 = buffer.data(target + 1754);
    auto *t_1755 = buffer.data(target + 1755);
    auto *t_1756 = buffer.data(target + 1756);
    auto *t_1757 = buffer.data(target + 1757);
    auto *t_1758 = buffer.data(target + 1758);
    auto *t_1759 = buffer.data(target + 1759);
    auto *t_1760 = buffer.data(target + 1760);
    auto *t_1761 = buffer.data(target + 1761);
    auto *t_1762 = buffer.data(target + 1762);
    auto *t_1763 = buffer.data(target + 1763);
    auto *t_1764 = buffer.data(target + 1764);
    auto *t_1765 = buffer.data(target + 1765);
    auto *t_1766 = buffer.data(target + 1766);
    auto *t_1767 = buffer.data(target + 1767);
    auto *t_1768 = buffer.data(target + 1768);
    auto *t_1769 = buffer.data(target + 1769);
    auto *t_1770 = buffer.data(target + 1770);
    auto *t_1771 = buffer.data(target + 1771);
    auto *t_1772 = buffer.data(target + 1772);
    auto *t_1773 = buffer.data(target + 1773);
    auto *t_1774 = buffer.data(target + 1774);
    auto *t_1775 = buffer.data(target + 1775);
    auto *t_1776 = buffer.data(target + 1776);
    auto *t_1777 = buffer.data(target + 1777);
    auto *t_1778 = buffer.data(target + 1778);
    auto *t_1779 = buffer.data(target + 1779);
    auto *t_1780 = buffer.data(target + 1780);
    auto *t_1781 = buffer.data(target + 1781);
    auto *t_1782 = buffer.data(target + 1782);
    auto *t_1783 = buffer.data(target + 1783);
    auto *t_1784 = buffer.data(target + 1784);
    auto *t_1785 = buffer.data(target + 1785);
    auto *t_1786 = buffer.data(target + 1786);
    auto *t_1787 = buffer.data(target + 1787);
    auto *t_1788 = buffer.data(target + 1788);
    auto *t_1789 = buffer.data(target + 1789);
    auto *t_1790 = buffer.data(target + 1790);
    auto *t_1791 = buffer.data(target + 1791);
    auto *t_1792 = buffer.data(target + 1792);
    auto *t_1793 = buffer.data(target + 1793);
    auto *t_1794 = buffer.data(target + 1794);
    auto *t_1795 = buffer.data(target + 1795);
    auto *t_1796 = buffer.data(target + 1796);
    auto *t_1797 = buffer.data(target + 1797);
    auto *t_1798 = buffer.data(target + 1798);
    auto *t_1799 = buffer.data(target + 1799);
    auto *t_1800 = buffer.data(target + 1800);
    auto *t_1801 = buffer.data(target + 1801);
    auto *t_1802 = buffer.data(target + 1802);
    auto *t_1803 = buffer.data(target + 1803);
    auto *t_1804 = buffer.data(target + 1804);
    auto *t_1805 = buffer.data(target + 1805);
    auto *t_1806 = buffer.data(target + 1806);
    auto *t_1807 = buffer.data(target + 1807);
    auto *t_1808 = buffer.data(target + 1808);
    auto *t_1809 = buffer.data(target + 1809);
    auto *t_1810 = buffer.data(target + 1810);
    auto *t_1811 = buffer.data(target + 1811);
    auto *t_1812 = buffer.data(target + 1812);
    auto *t_1813 = buffer.data(target + 1813);
    auto *t_1814 = buffer.data(target + 1814);
    auto *t_1815 = buffer.data(target + 1815);
    auto *t_1816 = buffer.data(target + 1816);
    auto *t_1817 = buffer.data(target + 1817);
    auto *t_1818 = buffer.data(target + 1818);
    auto *t_1819 = buffer.data(target + 1819);
    auto *t_1820 = buffer.data(target + 1820);
    auto *t_1821 = buffer.data(target + 1821);
    auto *t_1822 = buffer.data(target + 1822);
    auto *t_1823 = buffer.data(target + 1823);
    auto *t_1824 = buffer.data(target + 1824);
    auto *t_1825 = buffer.data(target + 1825);
    auto *t_1826 = buffer.data(target + 1826);
    auto *t_1827 = buffer.data(target + 1827);
    auto *t_1828 = buffer.data(target + 1828);
    auto *t_1829 = buffer.data(target + 1829);
    auto *t_1830 = buffer.data(target + 1830);
    auto *t_1831 = buffer.data(target + 1831);
    auto *t_1832 = buffer.data(target + 1832);
    auto *t_1833 = buffer.data(target + 1833);
    auto *t_1834 = buffer.data(target + 1834);
    auto *t_1835 = buffer.data(target + 1835);
    auto *t_1836 = buffer.data(target + 1836);
    auto *t_1837 = buffer.data(target + 1837);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_1512 = buffer.data(nsi0 + 1512);
    const auto *nsi0_1517 = buffer.data(nsi0 + 1517);
    const auto *nsi0_1521 = buffer.data(nsi0 + 1521);
    const auto *nsi0_1526 = buffer.data(nsi0 + 1526);
    const auto *nsi0_1720 = buffer.data(nsi0 + 1720);
    const auto *nsi0_1722 = buffer.data(nsi0 + 1722);
    const auto *nsi0_1729 = buffer.data(nsi0 + 1729);
    const auto *nsi0_1731 = buffer.data(nsi0 + 1731);
    const auto *nsi0_1732 = buffer.data(nsi0 + 1732);
    const auto *nsi0_1733 = buffer.data(nsi0 + 1733);
    const auto *nsi0_1735 = buffer.data(nsi0 + 1735);
    const auto *nsi0_1736 = buffer.data(nsi0 + 1736);
    const auto *nsi0_1739 = buffer.data(nsi0 + 1739);
    const auto *nsi0_1741 = buffer.data(nsi0 + 1741);
    const auto *nsi0_1742 = buffer.data(nsi0 + 1742);
    const auto *nsi0_1745 = buffer.data(nsi0 + 1745);
    const auto *nsi0_1746 = buffer.data(nsi0 + 1746);
    const auto *nsi0_1748 = buffer.data(nsi0 + 1748);
    const auto *nsi0_1750 = buffer.data(nsi0 + 1750);
    const auto *nsi0_1757 = buffer.data(nsi0 + 1757);
    const auto *nsi0_1759 = buffer.data(nsi0 + 1759);
    const auto *nsi0_1760 = buffer.data(nsi0 + 1760);
    const auto *nsi0_1761 = buffer.data(nsi0 + 1761);
    const auto *nsi0_1763 = buffer.data(nsi0 + 1763);
    const auto *nsi0_1764 = buffer.data(nsi0 + 1764);
    const auto *nsi0_1767 = buffer.data(nsi0 + 1767);
    const auto *nsi0_1769 = buffer.data(nsi0 + 1769);
    const auto *nsi0_1770 = buffer.data(nsi0 + 1770);
    const auto *nsi0_1773 = buffer.data(nsi0 + 1773);
    const auto *nsi0_1774 = buffer.data(nsi0 + 1774);
    const auto *nsi0_1776 = buffer.data(nsi0 + 1776);
    const auto *nsi0_1778 = buffer.data(nsi0 + 1778);
    const auto *nsi0_1785 = buffer.data(nsi0 + 1785);
    const auto *nsi0_1787 = buffer.data(nsi0 + 1787);
    const auto *nsi0_1788 = buffer.data(nsi0 + 1788);
    const auto *nsi0_1789 = buffer.data(nsi0 + 1789);
    const auto *nsi0_1791 = buffer.data(nsi0 + 1791);
    const auto *nsi0_1795 = buffer.data(nsi0 + 1795);
    const auto *nsi0_1798 = buffer.data(nsi0 + 1798);
    const auto *nsi0_1802 = buffer.data(nsi0 + 1802);
    const auto *nsi0_1804 = buffer.data(nsi0 + 1804);
    const auto *nsi0_1813 = buffer.data(nsi0 + 1813);
    const auto *nsi0_1815 = buffer.data(nsi0 + 1815);
    const auto *nsi0_1816 = buffer.data(nsi0 + 1816);
    const auto *nsi0_1817 = buffer.data(nsi0 + 1817);
    const auto *nsi0_1819 = buffer.data(nsi0 + 1819);
    const auto *nsi0_1820 = buffer.data(nsi0 + 1820);
    const auto *nsi0_1825 = buffer.data(nsi0 + 1825);
    const auto *nsi0_1829 = buffer.data(nsi0 + 1829);
    const auto *nsi0_1834 = buffer.data(nsi0 + 1834);

    const auto *nsh_1056 = buffer.data(nsh + 1056);
    const auto *nsh_1065 = buffer.data(nsh + 1065);
    const auto *nsh_1071 = buffer.data(nsh + 1071);
    const auto *nsh_1074 = buffer.data(nsh + 1074);
    const auto *nsh_1077 = buffer.data(nsh + 1077);
    const auto *nsh_1080 = buffer.data(nsh + 1080);
    const auto *nsh_1086 = buffer.data(nsh + 1086);
    const auto *nsh_1091 = buffer.data(nsh + 1091);
    const auto *nsh_1092 = buffer.data(nsh + 1092);
    const auto *nsh_1094 = buffer.data(nsh + 1094);
    const auto *nsh_1095 = buffer.data(nsh + 1095);
    const auto *nsh_1097 = buffer.data(nsh + 1097);
    const auto *nsh_1098 = buffer.data(nsh + 1098);
    const auto *nsh_1101 = buffer.data(nsh + 1101);
    const auto *nsh_1107 = buffer.data(nsh + 1107);
    const auto *nsh_1112 = buffer.data(nsh + 1112);
    const auto *nsh_1113 = buffer.data(nsh + 1113);
    const auto *nsh_1115 = buffer.data(nsh + 1115);
    const auto *nsh_1116 = buffer.data(nsh + 1116);
    const auto *nsh_1118 = buffer.data(nsh + 1118);
    const auto *nsh_1119 = buffer.data(nsh + 1119);
    const auto *nsh_1122 = buffer.data(nsh + 1122);
    const auto *nsh_1128 = buffer.data(nsh + 1128);
    const auto *nsh_1133 = buffer.data(nsh + 1133);
    const auto *nsh_1134 = buffer.data(nsh + 1134);
    const auto *nsh_1136 = buffer.data(nsh + 1136);
    const auto *nsh_1139 = buffer.data(nsh + 1139);
    const auto *nsh_1143 = buffer.data(nsh + 1143);
    const auto *nsh_1154 = buffer.data(nsh + 1154);
    const auto *nsh_1293 = buffer.data(nsh + 1293);
    const auto *nsh_1295 = buffer.data(nsh + 1295);
    const auto *nsh_1296 = buffer.data(nsh + 1296);
    const auto *nsh_1297 = buffer.data(nsh + 1297);
    const auto *nsh_1298 = buffer.data(nsh + 1298);
    const auto *nsh_1299 = buffer.data(nsh + 1299);
    const auto *nsh_1300 = buffer.data(nsh + 1300);
    const auto *nsh_1301 = buffer.data(nsh + 1301);
    const auto *nsh_1302 = buffer.data(nsh + 1302);
    const auto *nsh_1305 = buffer.data(nsh + 1305);
    const auto *nsh_1307 = buffer.data(nsh + 1307);
    const auto *nsh_1308 = buffer.data(nsh + 1308);
    const auto *nsh_1311 = buffer.data(nsh + 1311);
    const auto *nsh_1312 = buffer.data(nsh + 1312);
    const auto *nsh_1314 = buffer.data(nsh + 1314);
    const auto *nsh_1316 = buffer.data(nsh + 1316);
    const auto *nsh_1317 = buffer.data(nsh + 1317);
    const auto *nsh_1318 = buffer.data(nsh + 1318);
    const auto *nsh_1319 = buffer.data(nsh + 1319);
    const auto *nsh_1320 = buffer.data(nsh + 1320);
    const auto *nsh_1321 = buffer.data(nsh + 1321);
    const auto *nsh_1322 = buffer.data(nsh + 1322);
    const auto *nsh_1323 = buffer.data(nsh + 1323);
    const auto *nsh_1326 = buffer.data(nsh + 1326);
    const auto *nsh_1328 = buffer.data(nsh + 1328);
    const auto *nsh_1329 = buffer.data(nsh + 1329);
    const auto *nsh_1332 = buffer.data(nsh + 1332);
    const auto *nsh_1333 = buffer.data(nsh + 1333);
    const auto *nsh_1335 = buffer.data(nsh + 1335);
    const auto *nsh_1337 = buffer.data(nsh + 1337);
    const auto *nsh_1338 = buffer.data(nsh + 1338);
    const auto *nsh_1339 = buffer.data(nsh + 1339);
    const auto *nsh_1340 = buffer.data(nsh + 1340);
    const auto *nsh_1341 = buffer.data(nsh + 1341);
    const auto *nsh_1342 = buffer.data(nsh + 1342);
    const auto *nsh_1343 = buffer.data(nsh + 1343);
    const auto *nsh_1347 = buffer.data(nsh + 1347);
    const auto *nsh_1350 = buffer.data(nsh + 1350);
    const auto *nsh_1354 = buffer.data(nsh + 1354);
    const auto *nsh_1356 = buffer.data(nsh + 1356);
    const auto *nsh_1359 = buffer.data(nsh + 1359);
    const auto *nsh_1360 = buffer.data(nsh + 1360);
    const auto *nsh_1361 = buffer.data(nsh + 1361);
    const auto *nsh_1362 = buffer.data(nsh + 1362);
    const auto *nsh_1363 = buffer.data(nsh + 1363);
    const auto *nsh_1364 = buffer.data(nsh + 1364);
    const auto *nsh_1365 = buffer.data(nsh + 1365);
    const auto *nsh_1370 = buffer.data(nsh + 1370);
    const auto *nsh_1374 = buffer.data(nsh + 1374);
    const auto *nsh_1379 = buffer.data(nsh + 1379);
    const auto *nsh_1380 = buffer.data(nsh + 1380);
    const auto *nsh_1381 = buffer.data(nsh + 1381);
    const auto *nsh_1382 = buffer.data(nsh + 1382);

    const auto *nsi1_1512 = buffer.data(nsi1 + 1512);
    const auto *nsi1_1517 = buffer.data(nsi1 + 1517);
    const auto *nsi1_1521 = buffer.data(nsi1 + 1521);
    const auto *nsi1_1526 = buffer.data(nsi1 + 1526);
    const auto *nsi1_1720 = buffer.data(nsi1 + 1720);
    const auto *nsi1_1722 = buffer.data(nsi1 + 1722);
    const auto *nsi1_1729 = buffer.data(nsi1 + 1729);
    const auto *nsi1_1731 = buffer.data(nsi1 + 1731);
    const auto *nsi1_1732 = buffer.data(nsi1 + 1732);
    const auto *nsi1_1733 = buffer.data(nsi1 + 1733);
    const auto *nsi1_1735 = buffer.data(nsi1 + 1735);
    const auto *nsi1_1736 = buffer.data(nsi1 + 1736);
    const auto *nsi1_1739 = buffer.data(nsi1 + 1739);
    const auto *nsi1_1741 = buffer.data(nsi1 + 1741);
    const auto *nsi1_1742 = buffer.data(nsi1 + 1742);
    const auto *nsi1_1745 = buffer.data(nsi1 + 1745);
    const auto *nsi1_1746 = buffer.data(nsi1 + 1746);
    const auto *nsi1_1748 = buffer.data(nsi1 + 1748);
    const auto *nsi1_1750 = buffer.data(nsi1 + 1750);
    const auto *nsi1_1757 = buffer.data(nsi1 + 1757);
    const auto *nsi1_1759 = buffer.data(nsi1 + 1759);
    const auto *nsi1_1760 = buffer.data(nsi1 + 1760);
    const auto *nsi1_1761 = buffer.data(nsi1 + 1761);
    const auto *nsi1_1763 = buffer.data(nsi1 + 1763);
    const auto *nsi1_1764 = buffer.data(nsi1 + 1764);
    const auto *nsi1_1767 = buffer.data(nsi1 + 1767);
    const auto *nsi1_1769 = buffer.data(nsi1 + 1769);
    const auto *nsi1_1770 = buffer.data(nsi1 + 1770);
    const auto *nsi1_1773 = buffer.data(nsi1 + 1773);
    const auto *nsi1_1774 = buffer.data(nsi1 + 1774);
    const auto *nsi1_1776 = buffer.data(nsi1 + 1776);
    const auto *nsi1_1778 = buffer.data(nsi1 + 1778);
    const auto *nsi1_1785 = buffer.data(nsi1 + 1785);
    const auto *nsi1_1787 = buffer.data(nsi1 + 1787);
    const auto *nsi1_1788 = buffer.data(nsi1 + 1788);
    const auto *nsi1_1789 = buffer.data(nsi1 + 1789);
    const auto *nsi1_1791 = buffer.data(nsi1 + 1791);
    const auto *nsi1_1795 = buffer.data(nsi1 + 1795);
    const auto *nsi1_1798 = buffer.data(nsi1 + 1798);
    const auto *nsi1_1802 = buffer.data(nsi1 + 1802);
    const auto *nsi1_1804 = buffer.data(nsi1 + 1804);
    const auto *nsi1_1813 = buffer.data(nsi1 + 1813);
    const auto *nsi1_1815 = buffer.data(nsi1 + 1815);
    const auto *nsi1_1816 = buffer.data(nsi1 + 1816);
    const auto *nsi1_1817 = buffer.data(nsi1 + 1817);
    const auto *nsi1_1819 = buffer.data(nsi1 + 1819);
    const auto *nsi1_1820 = buffer.data(nsi1 + 1820);
    const auto *nsi1_1825 = buffer.data(nsi1 + 1825);
    const auto *nsi1_1829 = buffer.data(nsi1 + 1829);
    const auto *nsi1_1834 = buffer.data(nsi1 + 1834);

    const auto *osg0_975 = buffer.data(osg0 + 975);
    const auto *osg0_976 = buffer.data(osg0 + 976);
    const auto *osg0_977 = buffer.data(osg0 + 977);
    const auto *osg0_978 = buffer.data(osg0 + 978);
    const auto *osg0_979 = buffer.data(osg0 + 979);
    const auto *osg0_980 = buffer.data(osg0 + 980);

    const auto *osg1_975 = buffer.data(osg1 + 975);
    const auto *osg1_976 = buffer.data(osg1 + 976);
    const auto *osg1_977 = buffer.data(osg1 + 977);
    const auto *osg1_978 = buffer.data(osg1 + 978);
    const auto *osg1_979 = buffer.data(osg1 + 979);
    const auto *osg1_980 = buffer.data(osg1 + 980);

    const auto *osh_1287 = buffer.data(osh + 1287);
    const auto *osh_1290 = buffer.data(osh + 1290);
    const auto *osh_1296 = buffer.data(osh + 1296);
    const auto *osh_1297 = buffer.data(osh + 1297);
    const auto *osh_1298 = buffer.data(osh + 1298);
    const auto *osh_1299 = buffer.data(osh + 1299);
    const auto *osh_1300 = buffer.data(osh + 1300);
    const auto *osh_1301 = buffer.data(osh + 1301);
    const auto *osh_1302 = buffer.data(osh + 1302);
    const auto *osh_1304 = buffer.data(osh + 1304);
    const auto *osh_1305 = buffer.data(osh + 1305);
    const auto *osh_1307 = buffer.data(osh + 1307);
    const auto *osh_1308 = buffer.data(osh + 1308);
    const auto *osh_1311 = buffer.data(osh + 1311);
    const auto *osh_1317 = buffer.data(osh + 1317);
    const auto *osh_1318 = buffer.data(osh + 1318);
    const auto *osh_1319 = buffer.data(osh + 1319);
    const auto *osh_1320 = buffer.data(osh + 1320);
    const auto *osh_1321 = buffer.data(osh + 1321);
    const auto *osh_1322 = buffer.data(osh + 1322);
    const auto *osh_1323 = buffer.data(osh + 1323);
    const auto *osh_1325 = buffer.data(osh + 1325);
    const auto *osh_1326 = buffer.data(osh + 1326);
    const auto *osh_1328 = buffer.data(osh + 1328);
    const auto *osh_1329 = buffer.data(osh + 1329);
    const auto *osh_1332 = buffer.data(osh + 1332);
    const auto *osh_1338 = buffer.data(osh + 1338);
    const auto *osh_1339 = buffer.data(osh + 1339);
    const auto *osh_1340 = buffer.data(osh + 1340);
    const auto *osh_1341 = buffer.data(osh + 1341);
    const auto *osh_1342 = buffer.data(osh + 1342);
    const auto *osh_1343 = buffer.data(osh + 1343);
    const auto *osh_1344 = buffer.data(osh + 1344);
    const auto *osh_1346 = buffer.data(osh + 1346);
    const auto *osh_1347 = buffer.data(osh + 1347);
    const auto *osh_1349 = buffer.data(osh + 1349);
    const auto *osh_1350 = buffer.data(osh + 1350);
    const auto *osh_1353 = buffer.data(osh + 1353);
    const auto *osh_1359 = buffer.data(osh + 1359);
    const auto *osh_1360 = buffer.data(osh + 1360);
    const auto *osh_1361 = buffer.data(osh + 1361);
    const auto *osh_1362 = buffer.data(osh + 1362);
    const auto *osh_1363 = buffer.data(osh + 1363);
    const auto *osh_1364 = buffer.data(osh + 1364);
    const auto *osh_1365 = buffer.data(osh + 1365);
    const auto *osh_1366 = buffer.data(osh + 1366);
    const auto *osh_1367 = buffer.data(osh + 1367);
    const auto *osh_1368 = buffer.data(osh + 1368);
    const auto *osh_1369 = buffer.data(osh + 1369);
    const auto *osh_1370 = buffer.data(osh + 1370);
    const auto *osh_1371 = buffer.data(osh + 1371);
    const auto *osh_1372 = buffer.data(osh + 1372);
    const auto *osh_1373 = buffer.data(osh + 1373);
    const auto *osh_1374 = buffer.data(osh + 1374);
    const auto *osh_1380 = buffer.data(osh + 1380);
    const auto *osh_1381 = buffer.data(osh + 1381);
    const auto *osh_1382 = buffer.data(osh + 1382);

#pragma omp simd aligned(t_1719, t_1720, t_1721, pa_x, pc_x, pc_y, pc_z, nsi0_1720, nsh_1056, \
                         nsh_1080, nsh_1293, nsi1_1720, osh_1287, \
                         osh_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1719[k] = f_21 * nsh_1056[k]
                    + f_3 * pc_z[k] * osh_1287[k];

        t_1720[k] = pa_x[k] * nsi0_1720[k]
                    + f_12 * nsh_1293[k]
                    - f_10 * pc_x[k] * nsi1_1720[k];

        t_1721[k] = f_14 * nsh_1080[k]
                    + f_3 * pc_y[k] * osh_1290[k];
    }

#pragma omp simd aligned(t_1722, t_1723, t_1724, t_1725, pa_x, pc_x, nsi0_1722, nsh_1295, \
                         nsh_1296, nsh_1297, nsh_1298, nsi1_1722, osh_1296, osh_1297, \
                         osh_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1722[k] = pa_x[k] * nsi0_1722[k]
                    + f_12 * nsh_1295[k]
                    - f_10 * pc_x[k] * nsi1_1722[k];

        t_1723[k] = f_11 * nsh_1296[k]
                    + f_3 * pc_x[k] * osh_1296[k];

        t_1724[k] = f_11 * nsh_1297[k]
                    + f_3 * pc_x[k] * osh_1297[k];

        t_1725[k] = f_11 * nsh_1298[k]
                    + f_3 * pc_x[k] * osh_1298[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, t_1729, pa_x, pc_x, nsi0_1729, nsh_1299, \
                         nsh_1300, nsh_1301, nsi1_1729, osh_1299, osh_1300, \
                         osh_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_11 * nsh_1299[k]
                    + f_3 * pc_x[k] * osh_1299[k];

        t_1727[k] = f_11 * nsh_1300[k]
                    + f_3 * pc_x[k] * osh_1300[k];

        t_1728[k] = f_11 * nsh_1301[k]
                    + f_3 * pc_x[k] * osh_1301[k];

        t_1729[k] = pa_x[k] * nsi0_1729[k]
                    - f_10 * pc_x[k] * nsi1_1729[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, t_1733, pa_x, pc_x, pc_z, nsi0_1731, \
                         nsi0_1732, nsi0_1733, nsh_1065, nsi1_1731, nsi1_1732, nsi1_1733, \
                         osh_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = f_21 * nsh_1065[k]
                    + f_3 * pc_z[k] * osh_1296[k];

        t_1731[k] = pa_x[k] * nsi0_1731[k]
                    - f_10 * pc_x[k] * nsi1_1731[k];

        t_1732[k] = pa_x[k] * nsi0_1732[k]
                    - f_10 * pc_x[k] * nsi1_1732[k];

        t_1733[k] = pa_x[k] * nsi0_1733[k]
                    - f_10 * pc_x[k] * nsi1_1733[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, t_1737, pa_x, pc_x, pc_y, nsi0_1735, \
                         nsi0_1736, nsh_1091, nsh_1092, nsh_1302, nsi1_1735, nsi1_1736, \
                         osh_1301, osh_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = f_14 * nsh_1091[k]
                    + f_3 * pc_y[k] * osh_1301[k];

        t_1735[k] = pa_x[k] * nsi0_1735[k]
                    - f_10 * pc_x[k] * nsi1_1735[k];

        t_1736[k] = pa_x[k] * nsi0_1736[k]
                    + f_21 * nsh_1302[k]
                    - f_10 * pc_x[k] * nsi1_1736[k];

        t_1737[k] = f_13 * nsh_1092[k]
                    + f_3 * pc_y[k] * osh_1302[k];
    }

#pragma omp simd aligned(t_1738, t_1739, t_1740, pa_x, pc_x, pc_y, pc_z, nsi0_1739, nsh_1071, \
                         nsh_1094, nsh_1305, nsi1_1739, osh_1302, \
                         osh_1304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1738[k] = f_20 * nsh_1071[k]
                    + f_3 * pc_z[k] * osh_1302[k];

        t_1739[k] = pa_x[k] * nsi0_1739[k]
                    + f_14 * nsh_1305[k]
                    - f_10 * pc_x[k] * nsi1_1739[k];

        t_1740[k] = f_13 * nsh_1094[k]
                    + f_3 * pc_y[k] * osh_1304[k];
    }

#pragma omp simd aligned(t_1741, t_1742, t_1743, pa_x, pc_x, pc_z, nsi0_1741, nsi0_1742, \
                         nsh_1074, nsh_1307, nsh_1308, nsi1_1741, nsi1_1742, \
                         osh_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1741[k] = pa_x[k] * nsi0_1741[k]
                    + f_14 * nsh_1307[k]
                    - f_10 * pc_x[k] * nsi1_1741[k];

        t_1742[k] = pa_x[k] * nsi0_1742[k]
                    + f_13 * nsh_1308[k]
                    - f_10 * pc_x[k] * nsi1_1742[k];

        t_1743[k] = f_20 * nsh_1074[k]
                    + f_3 * pc_z[k] * osh_1305[k];
    }

#pragma omp simd aligned(t_1744, t_1745, t_1746, pa_x, pc_x, pc_y, nsi0_1745, nsi0_1746, \
                         nsh_1097, nsh_1311, nsh_1312, nsi1_1745, nsi1_1746, \
                         osh_1307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1744[k] = f_13 * nsh_1097[k]
                    + f_3 * pc_y[k] * osh_1307[k];

        t_1745[k] = pa_x[k] * nsi0_1745[k]
                    + f_13 * nsh_1311[k]
                    - f_10 * pc_x[k] * nsi1_1745[k];

        t_1746[k] = pa_x[k] * nsi0_1746[k]
                    + f_12 * nsh_1312[k]
                    - f_10 * pc_x[k] * nsi1_1746[k];
    }

#pragma omp simd aligned(t_1747, t_1748, t_1749, pa_x, pc_x, pc_y, pc_z, nsi0_1748, nsh_1077, \
                         nsh_1101, nsh_1314, nsi1_1748, osh_1308, \
                         osh_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1747[k] = f_20 * nsh_1077[k]
                    + f_3 * pc_z[k] * osh_1308[k];

        t_1748[k] = pa_x[k] * nsi0_1748[k]
                    + f_12 * nsh_1314[k]
                    - f_10 * pc_x[k] * nsi1_1748[k];

        t_1749[k] = f_13 * nsh_1101[k]
                    + f_3 * pc_y[k] * osh_1311[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, pa_x, pc_x, nsi0_1750, nsh_1316, \
                         nsh_1317, nsh_1318, nsh_1319, nsi1_1750, osh_1317, osh_1318, \
                         osh_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = pa_x[k] * nsi0_1750[k]
                    + f_12 * nsh_1316[k]
                    - f_10 * pc_x[k] * nsi1_1750[k];

        t_1751[k] = f_11 * nsh_1317[k]
                    + f_3 * pc_x[k] * osh_1317[k];

        t_1752[k] = f_11 * nsh_1318[k]
                    + f_3 * pc_x[k] * osh_1318[k];

        t_1753[k] = f_11 * nsh_1319[k]
                    + f_3 * pc_x[k] * osh_1319[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, pa_x, pc_x, nsi0_1757, nsh_1320, \
                         nsh_1321, nsh_1322, nsi1_1757, osh_1320, osh_1321, \
                         osh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = f_11 * nsh_1320[k]
                    + f_3 * pc_x[k] * osh_1320[k];

        t_1755[k] = f_11 * nsh_1321[k]
                    + f_3 * pc_x[k] * osh_1321[k];

        t_1756[k] = f_11 * nsh_1322[k]
                    + f_3 * pc_x[k] * osh_1322[k];

        t_1757[k] = pa_x[k] * nsi0_1757[k]
                    - f_10 * pc_x[k] * nsi1_1757[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, t_1761, pa_x, pc_x, pc_z, nsi0_1759, \
                         nsi0_1760, nsi0_1761, nsh_1086, nsi1_1759, nsi1_1760, nsi1_1761, \
                         osh_1317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = f_20 * nsh_1086[k]
                    + f_3 * pc_z[k] * osh_1317[k];

        t_1759[k] = pa_x[k] * nsi0_1759[k]
                    - f_10 * pc_x[k] * nsi1_1759[k];

        t_1760[k] = pa_x[k] * nsi0_1760[k]
                    - f_10 * pc_x[k] * nsi1_1760[k];

        t_1761[k] = pa_x[k] * nsi0_1761[k]
                    - f_10 * pc_x[k] * nsi1_1761[k];
    }

#pragma omp simd aligned(t_1762, t_1763, t_1764, t_1765, pa_x, pc_x, pc_y, nsi0_1763, \
                         nsi0_1764, nsh_1112, nsh_1113, nsh_1323, nsi1_1763, nsi1_1764, \
                         osh_1322, osh_1323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1762[k] = f_13 * nsh_1112[k]
                    + f_3 * pc_y[k] * osh_1322[k];

        t_1763[k] = pa_x[k] * nsi0_1763[k]
                    - f_10 * pc_x[k] * nsi1_1763[k];

        t_1764[k] = pa_x[k] * nsi0_1764[k]
                    + f_21 * nsh_1323[k]
                    - f_10 * pc_x[k] * nsi1_1764[k];

        t_1765[k] = f_12 * nsh_1113[k]
                    + f_3 * pc_y[k] * osh_1323[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pa_x, pc_x, pc_y, pc_z, nsi0_1767, nsh_1092, \
                         nsh_1115, nsh_1326, nsi1_1767, osh_1323, \
                         osh_1325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_19 * nsh_1092[k]
                    + f_3 * pc_z[k] * osh_1323[k];

        t_1767[k] = pa_x[k] * nsi0_1767[k]
                    + f_14 * nsh_1326[k]
                    - f_10 * pc_x[k] * nsi1_1767[k];

        t_1768[k] = f_12 * nsh_1115[k]
                    + f_3 * pc_y[k] * osh_1325[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, pa_x, pc_x, pc_z, nsi0_1769, nsi0_1770, \
                         nsh_1095, nsh_1328, nsh_1329, nsi1_1769, nsi1_1770, \
                         osh_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = pa_x[k] * nsi0_1769[k]
                    + f_14 * nsh_1328[k]
                    - f_10 * pc_x[k] * nsi1_1769[k];

        t_1770[k] = pa_x[k] * nsi0_1770[k]
                    + f_13 * nsh_1329[k]
                    - f_10 * pc_x[k] * nsi1_1770[k];

        t_1771[k] = f_19 * nsh_1095[k]
                    + f_3 * pc_z[k] * osh_1326[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, pa_x, pc_x, pc_y, nsi0_1773, nsi0_1774, \
                         nsh_1118, nsh_1332, nsh_1333, nsi1_1773, nsi1_1774, \
                         osh_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = f_12 * nsh_1118[k]
                    + f_3 * pc_y[k] * osh_1328[k];

        t_1773[k] = pa_x[k] * nsi0_1773[k]
                    + f_13 * nsh_1332[k]
                    - f_10 * pc_x[k] * nsi1_1773[k];

        t_1774[k] = pa_x[k] * nsi0_1774[k]
                    + f_12 * nsh_1333[k]
                    - f_10 * pc_x[k] * nsi1_1774[k];
    }

#pragma omp simd aligned(t_1775, t_1776, t_1777, pa_x, pc_x, pc_y, pc_z, nsi0_1776, nsh_1098, \
                         nsh_1122, nsh_1335, nsi1_1776, osh_1329, \
                         osh_1332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1775[k] = f_19 * nsh_1098[k]
                    + f_3 * pc_z[k] * osh_1329[k];

        t_1776[k] = pa_x[k] * nsi0_1776[k]
                    + f_12 * nsh_1335[k]
                    - f_10 * pc_x[k] * nsi1_1776[k];

        t_1777[k] = f_12 * nsh_1122[k]
                    + f_3 * pc_y[k] * osh_1332[k];
    }

#pragma omp simd aligned(t_1778, t_1779, t_1780, t_1781, pa_x, pc_x, nsi0_1778, nsh_1337, \
                         nsh_1338, nsh_1339, nsh_1340, nsi1_1778, osh_1338, osh_1339, \
                         osh_1340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1778[k] = pa_x[k] * nsi0_1778[k]
                    + f_12 * nsh_1337[k]
                    - f_10 * pc_x[k] * nsi1_1778[k];

        t_1779[k] = f_11 * nsh_1338[k]
                    + f_3 * pc_x[k] * osh_1338[k];

        t_1780[k] = f_11 * nsh_1339[k]
                    + f_3 * pc_x[k] * osh_1339[k];

        t_1781[k] = f_11 * nsh_1340[k]
                    + f_3 * pc_x[k] * osh_1340[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, pa_x, pc_x, nsi0_1785, nsh_1341, \
                         nsh_1342, nsh_1343, nsi1_1785, osh_1341, osh_1342, \
                         osh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_11 * nsh_1341[k]
                    + f_3 * pc_x[k] * osh_1341[k];

        t_1783[k] = f_11 * nsh_1342[k]
                    + f_3 * pc_x[k] * osh_1342[k];

        t_1784[k] = f_11 * nsh_1343[k]
                    + f_3 * pc_x[k] * osh_1343[k];

        t_1785[k] = pa_x[k] * nsi0_1785[k]
                    - f_10 * pc_x[k] * nsi1_1785[k];
    }

#pragma omp simd aligned(t_1786, t_1787, t_1788, t_1789, pa_x, pc_x, pc_z, nsi0_1787, \
                         nsi0_1788, nsi0_1789, nsh_1107, nsi1_1787, nsi1_1788, nsi1_1789, \
                         osh_1338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1786[k] = f_19 * nsh_1107[k]
                    + f_3 * pc_z[k] * osh_1338[k];

        t_1787[k] = pa_x[k] * nsi0_1787[k]
                    - f_10 * pc_x[k] * nsi1_1787[k];

        t_1788[k] = pa_x[k] * nsi0_1788[k]
                    - f_10 * pc_x[k] * nsi1_1788[k];

        t_1789[k] = pa_x[k] * nsi0_1789[k]
                    - f_10 * pc_x[k] * nsi1_1789[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, t_1793, pa_x, pa_y, pc_x, pc_y, nsi0_1512, \
                         nsi0_1791, nsh_1133, nsh_1134, nsi1_1512, nsi1_1791, osh_1343, \
                         osh_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_12 * nsh_1133[k]
                    + f_3 * pc_y[k] * osh_1343[k];

        t_1791[k] = pa_x[k] * nsi0_1791[k]
                    - f_10 * pc_x[k] * nsi1_1791[k];

        t_1792[k] = pa_y[k] * nsi0_1512[k]
                    - f_10 * pc_y[k] * nsi1_1512[k];

        t_1793[k] = f_11 * nsh_1134[k]
                    + f_3 * pc_y[k] * osh_1344[k];
    }

#pragma omp simd aligned(t_1794, t_1795, t_1796, pa_x, pc_x, pc_y, pc_z, nsi0_1795, nsh_1113, \
                         nsh_1136, nsh_1347, nsi1_1795, osh_1344, \
                         osh_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1794[k] = f_18 * nsh_1113[k]
                    + f_3 * pc_z[k] * osh_1344[k];

        t_1795[k] = pa_x[k] * nsi0_1795[k]
                    + f_14 * nsh_1347[k]
                    - f_10 * pc_x[k] * nsi1_1795[k];

        t_1796[k] = f_11 * nsh_1136[k]
                    + f_3 * pc_y[k] * osh_1346[k];
    }

#pragma omp simd aligned(t_1797, t_1798, t_1799, pa_x, pa_y, pc_x, pc_y, pc_z, nsi0_1517, \
                         nsi0_1798, nsh_1116, nsh_1350, nsi1_1517, nsi1_1798, \
                         osh_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1797[k] = pa_y[k] * nsi0_1517[k]
                    - f_10 * pc_y[k] * nsi1_1517[k];

        t_1798[k] = pa_x[k] * nsi0_1798[k]
                    + f_13 * nsh_1350[k]
                    - f_10 * pc_x[k] * nsi1_1798[k];

        t_1799[k] = f_18 * nsh_1116[k]
                    + f_3 * pc_z[k] * osh_1347[k];
    }

#pragma omp simd aligned(t_1800, t_1801, t_1802, pa_x, pa_y, pc_x, pc_y, nsi0_1521, nsi0_1802, \
                         nsh_1139, nsh_1354, nsi1_1521, nsi1_1802, \
                         osh_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1800[k] = f_11 * nsh_1139[k]
                    + f_3 * pc_y[k] * osh_1349[k];

        t_1801[k] = pa_y[k] * nsi0_1521[k]
                    - f_10 * pc_y[k] * nsi1_1521[k];

        t_1802[k] = pa_x[k] * nsi0_1802[k]
                    + f_12 * nsh_1354[k]
                    - f_10 * pc_x[k] * nsi1_1802[k];
    }

#pragma omp simd aligned(t_1803, t_1804, t_1805, pa_x, pc_x, pc_y, pc_z, nsi0_1804, nsh_1119, \
                         nsh_1143, nsh_1356, nsi1_1804, osh_1350, \
                         osh_1353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1803[k] = f_18 * nsh_1119[k]
                    + f_3 * pc_z[k] * osh_1350[k];

        t_1804[k] = pa_x[k] * nsi0_1804[k]
                    + f_12 * nsh_1356[k]
                    - f_10 * pc_x[k] * nsi1_1804[k];

        t_1805[k] = f_11 * nsh_1143[k]
                    + f_3 * pc_y[k] * osh_1353[k];
    }

#pragma omp simd aligned(t_1806, t_1807, t_1808, t_1809, pa_y, pc_x, pc_y, nsi0_1526, \
                         nsh_1359, nsh_1360, nsh_1361, nsi1_1526, osh_1359, osh_1360, \
                         osh_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1806[k] = pa_y[k] * nsi0_1526[k]
                    - f_10 * pc_y[k] * nsi1_1526[k];

        t_1807[k] = f_11 * nsh_1359[k]
                    + f_3 * pc_x[k] * osh_1359[k];

        t_1808[k] = f_11 * nsh_1360[k]
                    + f_3 * pc_x[k] * osh_1360[k];

        t_1809[k] = f_11 * nsh_1361[k]
                    + f_3 * pc_x[k] * osh_1361[k];
    }

#pragma omp simd aligned(t_1810, t_1811, t_1812, t_1813, pa_x, pc_x, nsi0_1813, nsh_1362, \
                         nsh_1363, nsh_1364, nsi1_1813, osh_1362, osh_1363, \
                         osh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1810[k] = f_11 * nsh_1362[k]
                    + f_3 * pc_x[k] * osh_1362[k];

        t_1811[k] = f_11 * nsh_1363[k]
                    + f_3 * pc_x[k] * osh_1363[k];

        t_1812[k] = f_11 * nsh_1364[k]
                    + f_3 * pc_x[k] * osh_1364[k];

        t_1813[k] = pa_x[k] * nsi0_1813[k]
                    - f_10 * pc_x[k] * nsi1_1813[k];
    }

#pragma omp simd aligned(t_1814, t_1815, t_1816, t_1817, pa_x, pc_x, pc_z, nsi0_1815, \
                         nsi0_1816, nsi0_1817, nsh_1128, nsi1_1815, nsi1_1816, nsi1_1817, \
                         osh_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = f_18 * nsh_1128[k]
                    + f_3 * pc_z[k] * osh_1359[k];

        t_1815[k] = pa_x[k] * nsi0_1815[k]
                    - f_10 * pc_x[k] * nsi1_1815[k];

        t_1816[k] = pa_x[k] * nsi0_1816[k]
                    - f_10 * pc_x[k] * nsi1_1816[k];

        t_1817[k] = pa_x[k] * nsi0_1817[k]
                    - f_10 * pc_x[k] * nsi1_1817[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, t_1821, pa_x, pc_x, pc_y, nsi0_1819, \
                         nsi0_1820, nsh_1154, nsh_1365, nsi1_1819, nsi1_1820, osh_1364, \
                         osh_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = f_11 * nsh_1154[k]
                    + f_3 * pc_y[k] * osh_1364[k];

        t_1819[k] = pa_x[k] * nsi0_1819[k]
                    - f_10 * pc_x[k] * nsi1_1819[k];

        t_1820[k] = pa_x[k] * nsi0_1820[k]
                    + f_21 * nsh_1365[k]
                    - f_10 * pc_x[k] * nsi1_1820[k];

        t_1821[k] = f_3 * pc_y[k] * osh_1365[k];
    }

#pragma omp simd aligned(t_1822, t_1823, t_1824, pc_y, pc_z, nsh_1134, osg0_975, osg1_975, \
                         osh_1365, osh_1366, osh_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1822[k] = f_15 * nsh_1134[k]
                    + f_3 * pc_z[k] * osh_1365[k];

        t_1823[k] = f_4 * osg0_975[k]
                    - f_5 * osg1_975[k]
                    + f_3 * pc_y[k] * osh_1366[k];

        t_1824[k] = f_3 * pc_y[k] * osh_1367[k];
    }

#pragma omp simd aligned(t_1825, t_1826, t_1827, pa_x, pc_x, pc_y, nsi0_1825, nsh_1370, \
                         nsi1_1825, osg0_976, osg0_977, osg1_976, osg1_977, osh_1368, \
                         osh_1369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1825[k] = pa_x[k] * nsi0_1825[k]
                    + f_14 * nsh_1370[k]
                    - f_10 * pc_x[k] * nsi1_1825[k];

        t_1826[k] = f_6 * osg0_976[k]
                    - f_7 * osg1_976[k]
                    + f_3 * pc_y[k] * osh_1368[k];

        t_1827[k] = f_4 * osg0_977[k]
                    - f_5 * osg1_977[k]
                    + f_3 * pc_y[k] * osh_1369[k];
    }

#pragma omp simd aligned(t_1828, t_1829, t_1830, pa_x, pc_x, pc_y, nsi0_1829, nsh_1374, \
                         nsi1_1829, osg0_978, osg1_978, osh_1370, \
                         osh_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1828[k] = f_3 * pc_y[k] * osh_1370[k];

        t_1829[k] = pa_x[k] * nsi0_1829[k]
                    + f_13 * nsh_1374[k]
                    - f_10 * pc_x[k] * nsi1_1829[k];

        t_1830[k] = f_8 * osg0_978[k]
                    - f_9 * osg1_978[k]
                    + f_3 * pc_y[k] * osh_1371[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, pc_y, osg0_979, osg0_980, osg1_979, osg1_980, \
                         osh_1372, osh_1373, osh_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_6 * osg0_979[k]
                    - f_7 * osg1_979[k]
                    + f_3 * pc_y[k] * osh_1372[k];

        t_1832[k] = f_4 * osg0_980[k]
                    - f_5 * osg1_980[k]
                    + f_3 * pc_y[k] * osh_1373[k];

        t_1833[k] = f_3 * pc_y[k] * osh_1374[k];
    }

#pragma omp simd aligned(t_1834, t_1835, t_1836, t_1837, pa_x, pc_x, nsi0_1834, nsh_1379, \
                         nsh_1380, nsh_1381, nsh_1382, nsi1_1834, osh_1380, osh_1381, \
                         osh_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = pa_x[k] * nsi0_1834[k]
                    + f_12 * nsh_1379[k]
                    - f_10 * pc_x[k] * nsi1_1834[k];

        t_1835[k] = f_11 * nsh_1380[k]
                    + f_3 * pc_x[k] * osh_1380[k];

        t_1836[k] = f_11 * nsh_1381[k]
                    + f_3 * pc_x[k] * osh_1381[k];

        t_1837[k] = f_11 * nsh_1382[k]
                    + f_3 * pc_x[k] * osh_1382[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osg0, const size_t osg1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;

    auto *t_1838 = buffer.data(target + 1838);
    auto *t_1839 = buffer.data(target + 1839);
    auto *t_1840 = buffer.data(target + 1840);
    auto *t_1841 = buffer.data(target + 1841);
    auto *t_1842 = buffer.data(target + 1842);
    auto *t_1843 = buffer.data(target + 1843);
    auto *t_1844 = buffer.data(target + 1844);
    auto *t_1845 = buffer.data(target + 1845);
    auto *t_1846 = buffer.data(target + 1846);
    auto *t_1847 = buffer.data(target + 1847);
    auto *t_1848 = buffer.data(target + 1848);
    auto *t_1849 = buffer.data(target + 1849);
    auto *t_1850 = buffer.data(target + 1850);
    auto *t_1851 = buffer.data(target + 1851);
    auto *t_1852 = buffer.data(target + 1852);
    auto *t_1853 = buffer.data(target + 1853);
    auto *t_1854 = buffer.data(target + 1854);
    auto *t_1855 = buffer.data(target + 1855);
    auto *t_1856 = buffer.data(target + 1856);
    auto *t_1857 = buffer.data(target + 1857);
    auto *t_1858 = buffer.data(target + 1858);
    auto *t_1859 = buffer.data(target + 1859);
    auto *t_1860 = buffer.data(target + 1860);
    auto *t_1861 = buffer.data(target + 1861);
    auto *t_1862 = buffer.data(target + 1862);
    auto *t_1863 = buffer.data(target + 1863);
    auto *t_1864 = buffer.data(target + 1864);
    auto *t_1865 = buffer.data(target + 1865);
    auto *t_1866 = buffer.data(target + 1866);
    auto *t_1867 = buffer.data(target + 1867);
    auto *t_1868 = buffer.data(target + 1868);
    auto *t_1869 = buffer.data(target + 1869);
    auto *t_1870 = buffer.data(target + 1870);
    auto *t_1871 = buffer.data(target + 1871);
    auto *t_1872 = buffer.data(target + 1872);
    auto *t_1873 = buffer.data(target + 1873);
    auto *t_1874 = buffer.data(target + 1874);
    auto *t_1875 = buffer.data(target + 1875);
    auto *t_1876 = buffer.data(target + 1876);
    auto *t_1877 = buffer.data(target + 1877);
    auto *t_1878 = buffer.data(target + 1878);
    auto *t_1879 = buffer.data(target + 1879);
    auto *t_1880 = buffer.data(target + 1880);
    auto *t_1881 = buffer.data(target + 1881);
    auto *t_1882 = buffer.data(target + 1882);
    auto *t_1883 = buffer.data(target + 1883);
    auto *t_1884 = buffer.data(target + 1884);
    auto *t_1885 = buffer.data(target + 1885);
    auto *t_1886 = buffer.data(target + 1886);
    auto *t_1887 = buffer.data(target + 1887);
    auto *t_1888 = buffer.data(target + 1888);
    auto *t_1889 = buffer.data(target + 1889);
    auto *t_1890 = buffer.data(target + 1890);
    auto *t_1891 = buffer.data(target + 1891);
    auto *t_1892 = buffer.data(target + 1892);
    auto *t_1893 = buffer.data(target + 1893);
    auto *t_1894 = buffer.data(target + 1894);
    auto *t_1895 = buffer.data(target + 1895);
    auto *t_1896 = buffer.data(target + 1896);
    auto *t_1897 = buffer.data(target + 1897);
    auto *t_1898 = buffer.data(target + 1898);
    auto *t_1899 = buffer.data(target + 1899);
    auto *t_1900 = buffer.data(target + 1900);
    auto *t_1901 = buffer.data(target + 1901);
    auto *t_1902 = buffer.data(target + 1902);
    auto *t_1903 = buffer.data(target + 1903);
    auto *t_1904 = buffer.data(target + 1904);
    auto *t_1905 = buffer.data(target + 1905);
    auto *t_1906 = buffer.data(target + 1906);
    auto *t_1907 = buffer.data(target + 1907);
    auto *t_1908 = buffer.data(target + 1908);
    auto *t_1909 = buffer.data(target + 1909);
    auto *t_1910 = buffer.data(target + 1910);
    auto *t_1911 = buffer.data(target + 1911);
    auto *t_1912 = buffer.data(target + 1912);
    auto *t_1913 = buffer.data(target + 1913);
    auto *t_1914 = buffer.data(target + 1914);
    auto *t_1915 = buffer.data(target + 1915);
    auto *t_1916 = buffer.data(target + 1916);
    auto *t_1917 = buffer.data(target + 1917);
    auto *t_1918 = buffer.data(target + 1918);
    auto *t_1919 = buffer.data(target + 1919);
    auto *t_1920 = buffer.data(target + 1920);
    auto *t_1921 = buffer.data(target + 1921);
    auto *t_1922 = buffer.data(target + 1922);
    auto *t_1923 = buffer.data(target + 1923);
    auto *t_1924 = buffer.data(target + 1924);
    auto *t_1925 = buffer.data(target + 1925);
    auto *t_1926 = buffer.data(target + 1926);
    auto *t_1927 = buffer.data(target + 1927);
    auto *t_1928 = buffer.data(target + 1928);
    auto *t_1929 = buffer.data(target + 1929);
    auto *t_1930 = buffer.data(target + 1930);
    auto *t_1931 = buffer.data(target + 1931);
    auto *t_1932 = buffer.data(target + 1932);
    auto *t_1933 = buffer.data(target + 1933);
    auto *t_1934 = buffer.data(target + 1934);
    auto *t_1935 = buffer.data(target + 1935);
    auto *t_1936 = buffer.data(target + 1936);
    auto *t_1937 = buffer.data(target + 1937);
    auto *t_1938 = buffer.data(target + 1938);
    auto *t_1939 = buffer.data(target + 1939);
    auto *t_1940 = buffer.data(target + 1940);
    auto *t_1941 = buffer.data(target + 1941);
    auto *t_1942 = buffer.data(target + 1942);
    auto *t_1943 = buffer.data(target + 1943);
    auto *t_1944 = buffer.data(target + 1944);
    auto *t_1945 = buffer.data(target + 1945);
    auto *t_1946 = buffer.data(target + 1946);
    auto *t_1947 = buffer.data(target + 1947);
    auto *t_1948 = buffer.data(target + 1948);
    auto *t_1949 = buffer.data(target + 1949);
    auto *t_1950 = buffer.data(target + 1950);
    auto *t_1951 = buffer.data(target + 1951);
    auto *t_1952 = buffer.data(target + 1952);
    auto *t_1953 = buffer.data(target + 1953);
    auto *t_1954 = buffer.data(target + 1954);
    auto *t_1955 = buffer.data(target + 1955);
    auto *t_1956 = buffer.data(target + 1956);
    auto *t_1957 = buffer.data(target + 1957);
    auto *t_1958 = buffer.data(target + 1958);
    auto *t_1959 = buffer.data(target + 1959);
    auto *t_1960 = buffer.data(target + 1960);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_1540 = buffer.data(nsi0 + 1540);
    const auto *nsi0_1541 = buffer.data(nsi0 + 1541);
    const auto *nsi0_1543 = buffer.data(nsi0 + 1543);
    const auto *nsi0_1546 = buffer.data(nsi0 + 1546);
    const auto *nsi0_1550 = buffer.data(nsi0 + 1550);
    const auto *nsi0_1561 = buffer.data(nsi0 + 1561);
    const auto *nsi0_1563 = buffer.data(nsi0 + 1563);
    const auto *nsi0_1564 = buffer.data(nsi0 + 1564);
    const auto *nsi0_1565 = buffer.data(nsi0 + 1565);
    const auto *nsi0_1841 = buffer.data(nsi0 + 1841);
    const auto *nsi0_1842 = buffer.data(nsi0 + 1842);
    const auto *nsi0_1843 = buffer.data(nsi0 + 1843);
    const auto *nsi0_1844 = buffer.data(nsi0 + 1844);
    const auto *nsi0_1845 = buffer.data(nsi0 + 1845);
    const auto *nsi0_1847 = buffer.data(nsi0 + 1847);

    const auto *nsh_1170 = buffer.data(nsh + 1170);
    const auto *nsh_1171 = buffer.data(nsh + 1171);
    const auto *nsh_1172 = buffer.data(nsh + 1172);
    const auto *nsh_1173 = buffer.data(nsh + 1173);
    const auto *nsh_1175 = buffer.data(nsh + 1175);
    const auto *nsh_1191 = buffer.data(nsh + 1191);
    const auto *nsh_1196 = buffer.data(nsh + 1196);
    const auto *nsh_1212 = buffer.data(nsh + 1212);
    const auto *nsh_1214 = buffer.data(nsh + 1214);
    const auto *nsh_1215 = buffer.data(nsh + 1215);
    const auto *nsh_1216 = buffer.data(nsh + 1216);
    const auto *nsh_1217 = buffer.data(nsh + 1217);
    const auto *nsh_1233 = buffer.data(nsh + 1233);
    const auto *nsh_1235 = buffer.data(nsh + 1235);
    const auto *nsh_1236 = buffer.data(nsh + 1236);
    const auto *nsh_1237 = buffer.data(nsh + 1237);
    const auto *nsh_1238 = buffer.data(nsh + 1238);
    const auto *nsh_1383 = buffer.data(nsh + 1383);
    const auto *nsh_1385 = buffer.data(nsh + 1385);

    const auto *nsi1_1540 = buffer.data(nsi1 + 1540);
    const auto *nsi1_1541 = buffer.data(nsi1 + 1541);
    const auto *nsi1_1543 = buffer.data(nsi1 + 1543);
    const auto *nsi1_1546 = buffer.data(nsi1 + 1546);
    const auto *nsi1_1550 = buffer.data(nsi1 + 1550);
    const auto *nsi1_1561 = buffer.data(nsi1 + 1561);
    const auto *nsi1_1563 = buffer.data(nsi1 + 1563);
    const auto *nsi1_1564 = buffer.data(nsi1 + 1564);
    const auto *nsi1_1565 = buffer.data(nsi1 + 1565);
    const auto *nsi1_1841 = buffer.data(nsi1 + 1841);
    const auto *nsi1_1842 = buffer.data(nsi1 + 1842);
    const auto *nsi1_1843 = buffer.data(nsi1 + 1843);
    const auto *nsi1_1844 = buffer.data(nsi1 + 1844);
    const auto *nsi1_1845 = buffer.data(nsi1 + 1845);
    const auto *nsi1_1847 = buffer.data(nsi1 + 1847);

    const auto *osg0_990 = buffer.data(osg0 + 990);
    const auto *osg0_991 = buffer.data(osg0 + 991);
    const auto *osg0_993 = buffer.data(osg0 + 993);
    const auto *osg0_995 = buffer.data(osg0 + 995);
    const auto *osg0_996 = buffer.data(osg0 + 996);
    const auto *osg0_998 = buffer.data(osg0 + 998);
    const auto *osg0_999 = buffer.data(osg0 + 999);
    const auto *osg0_1000 = buffer.data(osg0 + 1000);
    const auto *osg0_1001 = buffer.data(osg0 + 1001);
    const auto *osg0_1002 = buffer.data(osg0 + 1002);
    const auto *osg0_1003 = buffer.data(osg0 + 1003);
    const auto *osg0_1004 = buffer.data(osg0 + 1004);
    const auto *osg0_1007 = buffer.data(osg0 + 1007);
    const auto *osg0_1009 = buffer.data(osg0 + 1009);
    const auto *osg0_1010 = buffer.data(osg0 + 1010);
    const auto *osg0_1012 = buffer.data(osg0 + 1012);
    const auto *osg0_1013 = buffer.data(osg0 + 1013);
    const auto *osg0_1014 = buffer.data(osg0 + 1014);
    const auto *osg0_1016 = buffer.data(osg0 + 1016);
    const auto *osg0_1017 = buffer.data(osg0 + 1017);
    const auto *osg0_1018 = buffer.data(osg0 + 1018);
    const auto *osg0_1019 = buffer.data(osg0 + 1019);
    const auto *osg0_1020 = buffer.data(osg0 + 1020);
    const auto *osg0_1021 = buffer.data(osg0 + 1021);
    const auto *osg0_1022 = buffer.data(osg0 + 1022);
    const auto *osg0_1023 = buffer.data(osg0 + 1023);
    const auto *osg0_1024 = buffer.data(osg0 + 1024);
    const auto *osg0_1025 = buffer.data(osg0 + 1025);
    const auto *osg0_1026 = buffer.data(osg0 + 1026);
    const auto *osg0_1027 = buffer.data(osg0 + 1027);
    const auto *osg0_1028 = buffer.data(osg0 + 1028);
    const auto *osg0_1029 = buffer.data(osg0 + 1029);
    const auto *osg0_1030 = buffer.data(osg0 + 1030);
    const auto *osg0_1031 = buffer.data(osg0 + 1031);
    const auto *osg0_1032 = buffer.data(osg0 + 1032);
    const auto *osg0_1033 = buffer.data(osg0 + 1033);
    const auto *osg0_1034 = buffer.data(osg0 + 1034);
    const auto *osg0_1035 = buffer.data(osg0 + 1035);
    const auto *osg0_1036 = buffer.data(osg0 + 1036);
    const auto *osg0_1037 = buffer.data(osg0 + 1037);
    const auto *osg0_1038 = buffer.data(osg0 + 1038);
    const auto *osg0_1039 = buffer.data(osg0 + 1039);
    const auto *osg0_1040 = buffer.data(osg0 + 1040);
    const auto *osg0_1041 = buffer.data(osg0 + 1041);
    const auto *osg0_1042 = buffer.data(osg0 + 1042);
    const auto *osg0_1043 = buffer.data(osg0 + 1043);
    const auto *osg0_1044 = buffer.data(osg0 + 1044);
    const auto *osg0_1045 = buffer.data(osg0 + 1045);
    const auto *osg0_1046 = buffer.data(osg0 + 1046);
    const auto *osg0_1047 = buffer.data(osg0 + 1047);
    const auto *osg0_1048 = buffer.data(osg0 + 1048);
    const auto *osg0_1049 = buffer.data(osg0 + 1049);
    const auto *osg0_1050 = buffer.data(osg0 + 1050);

    const auto *osg1_990 = buffer.data(osg1 + 990);
    const auto *osg1_991 = buffer.data(osg1 + 991);
    const auto *osg1_993 = buffer.data(osg1 + 993);
    const auto *osg1_995 = buffer.data(osg1 + 995);
    const auto *osg1_996 = buffer.data(osg1 + 996);
    const auto *osg1_998 = buffer.data(osg1 + 998);
    const auto *osg1_999 = buffer.data(osg1 + 999);
    const auto *osg1_1000 = buffer.data(osg1 + 1000);
    const auto *osg1_1001 = buffer.data(osg1 + 1001);
    const auto *osg1_1002 = buffer.data(osg1 + 1002);
    const auto *osg1_1003 = buffer.data(osg1 + 1003);
    const auto *osg1_1004 = buffer.data(osg1 + 1004);
    const auto *osg1_1007 = buffer.data(osg1 + 1007);
    const auto *osg1_1009 = buffer.data(osg1 + 1009);
    const auto *osg1_1010 = buffer.data(osg1 + 1010);
    const auto *osg1_1012 = buffer.data(osg1 + 1012);
    const auto *osg1_1013 = buffer.data(osg1 + 1013);
    const auto *osg1_1014 = buffer.data(osg1 + 1014);
    const auto *osg1_1016 = buffer.data(osg1 + 1016);
    const auto *osg1_1017 = buffer.data(osg1 + 1017);
    const auto *osg1_1018 = buffer.data(osg1 + 1018);
    const auto *osg1_1019 = buffer.data(osg1 + 1019);
    const auto *osg1_1020 = buffer.data(osg1 + 1020);
    const auto *osg1_1021 = buffer.data(osg1 + 1021);
    const auto *osg1_1022 = buffer.data(osg1 + 1022);
    const auto *osg1_1023 = buffer.data(osg1 + 1023);
    const auto *osg1_1024 = buffer.data(osg1 + 1024);
    const auto *osg1_1025 = buffer.data(osg1 + 1025);
    const auto *osg1_1026 = buffer.data(osg1 + 1026);
    const auto *osg1_1027 = buffer.data(osg1 + 1027);
    const auto *osg1_1028 = buffer.data(osg1 + 1028);
    const auto *osg1_1029 = buffer.data(osg1 + 1029);
    const auto *osg1_1030 = buffer.data(osg1 + 1030);
    const auto *osg1_1031 = buffer.data(osg1 + 1031);
    const auto *osg1_1032 = buffer.data(osg1 + 1032);
    const auto *osg1_1033 = buffer.data(osg1 + 1033);
    const auto *osg1_1034 = buffer.data(osg1 + 1034);
    const auto *osg1_1035 = buffer.data(osg1 + 1035);
    const auto *osg1_1036 = buffer.data(osg1 + 1036);
    const auto *osg1_1037 = buffer.data(osg1 + 1037);
    const auto *osg1_1038 = buffer.data(osg1 + 1038);
    const auto *osg1_1039 = buffer.data(osg1 + 1039);
    const auto *osg1_1040 = buffer.data(osg1 + 1040);
    const auto *osg1_1041 = buffer.data(osg1 + 1041);
    const auto *osg1_1042 = buffer.data(osg1 + 1042);
    const auto *osg1_1043 = buffer.data(osg1 + 1043);
    const auto *osg1_1044 = buffer.data(osg1 + 1044);
    const auto *osg1_1045 = buffer.data(osg1 + 1045);
    const auto *osg1_1046 = buffer.data(osg1 + 1046);
    const auto *osg1_1047 = buffer.data(osg1 + 1047);
    const auto *osg1_1048 = buffer.data(osg1 + 1048);
    const auto *osg1_1049 = buffer.data(osg1 + 1049);
    const auto *osg1_1050 = buffer.data(osg1 + 1050);

    const auto *osh_1379 = buffer.data(osh + 1379);
    const auto *osh_1383 = buffer.data(osh + 1383);
    const auto *osh_1385 = buffer.data(osh + 1385);
    const auto *osh_1386 = buffer.data(osh + 1386);
    const auto *osh_1387 = buffer.data(osh + 1387);
    const auto *osh_1389 = buffer.data(osh + 1389);
    const auto *osh_1391 = buffer.data(osh + 1391);
    const auto *osh_1392 = buffer.data(osh + 1392);
    const auto *osh_1394 = buffer.data(osh + 1394);
    const auto *osh_1395 = buffer.data(osh + 1395);
    const auto *osh_1396 = buffer.data(osh + 1396);
    const auto *osh_1398 = buffer.data(osh + 1398);
    const auto *osh_1399 = buffer.data(osh + 1399);
    const auto *osh_1400 = buffer.data(osh + 1400);
    const auto *osh_1401 = buffer.data(osh + 1401);
    const auto *osh_1402 = buffer.data(osh + 1402);
    const auto *osh_1403 = buffer.data(osh + 1403);
    const auto *osh_1404 = buffer.data(osh + 1404);
    const auto *osh_1405 = buffer.data(osh + 1405);
    const auto *osh_1406 = buffer.data(osh + 1406);
    const auto *osh_1409 = buffer.data(osh + 1409);
    const auto *osh_1411 = buffer.data(osh + 1411);
    const auto *osh_1412 = buffer.data(osh + 1412);
    const auto *osh_1414 = buffer.data(osh + 1414);
    const auto *osh_1415 = buffer.data(osh + 1415);
    const auto *osh_1416 = buffer.data(osh + 1416);
    const auto *osh_1418 = buffer.data(osh + 1418);
    const auto *osh_1419 = buffer.data(osh + 1419);
    const auto *osh_1420 = buffer.data(osh + 1420);
    const auto *osh_1421 = buffer.data(osh + 1421);
    const auto *osh_1422 = buffer.data(osh + 1422);
    const auto *osh_1423 = buffer.data(osh + 1423);
    const auto *osh_1424 = buffer.data(osh + 1424);
    const auto *osh_1425 = buffer.data(osh + 1425);
    const auto *osh_1426 = buffer.data(osh + 1426);
    const auto *osh_1427 = buffer.data(osh + 1427);
    const auto *osh_1428 = buffer.data(osh + 1428);
    const auto *osh_1429 = buffer.data(osh + 1429);
    const auto *osh_1430 = buffer.data(osh + 1430);
    const auto *osh_1431 = buffer.data(osh + 1431);
    const auto *osh_1432 = buffer.data(osh + 1432);
    const auto *osh_1433 = buffer.data(osh + 1433);
    const auto *osh_1434 = buffer.data(osh + 1434);
    const auto *osh_1435 = buffer.data(osh + 1435);
    const auto *osh_1436 = buffer.data(osh + 1436);
    const auto *osh_1437 = buffer.data(osh + 1437);
    const auto *osh_1438 = buffer.data(osh + 1438);
    const auto *osh_1439 = buffer.data(osh + 1439);
    const auto *osh_1440 = buffer.data(osh + 1440);
    const auto *osh_1441 = buffer.data(osh + 1441);
    const auto *osh_1442 = buffer.data(osh + 1442);
    const auto *osh_1443 = buffer.data(osh + 1443);
    const auto *osh_1444 = buffer.data(osh + 1444);
    const auto *osh_1445 = buffer.data(osh + 1445);
    const auto *osh_1446 = buffer.data(osh + 1446);
    const auto *osh_1447 = buffer.data(osh + 1447);
    const auto *osh_1448 = buffer.data(osh + 1448);
    const auto *osh_1449 = buffer.data(osh + 1449);
    const auto *osh_1450 = buffer.data(osh + 1450);
    const auto *osh_1451 = buffer.data(osh + 1451);
    const auto *osh_1452 = buffer.data(osh + 1452);
    const auto *osh_1453 = buffer.data(osh + 1453);
    const auto *osh_1454 = buffer.data(osh + 1454);
    const auto *osh_1455 = buffer.data(osh + 1455);
    const auto *osh_1456 = buffer.data(osh + 1456);
    const auto *osh_1457 = buffer.data(osh + 1457);
    const auto *osh_1458 = buffer.data(osh + 1458);
    const auto *osh_1459 = buffer.data(osh + 1459);
    const auto *osh_1460 = buffer.data(osh + 1460);
    const auto *osh_1461 = buffer.data(osh + 1461);
    const auto *osh_1462 = buffer.data(osh + 1462);
    const auto *osh_1463 = buffer.data(osh + 1463);
    const auto *osh_1464 = buffer.data(osh + 1464);
    const auto *osh_1465 = buffer.data(osh + 1465);
    const auto *osh_1466 = buffer.data(osh + 1466);
    const auto *osh_1467 = buffer.data(osh + 1467);
    const auto *osh_1468 = buffer.data(osh + 1468);
    const auto *osh_1469 = buffer.data(osh + 1469);
    const auto *osh_1470 = buffer.data(osh + 1470);

#pragma omp simd aligned(t_1838, t_1839, t_1840, t_1841, pa_x, pc_x, pc_y, nsi0_1841, \
                         nsh_1383, nsh_1385, nsi1_1841, osh_1379, osh_1383, \
                         osh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1838[k] = f_11 * nsh_1383[k]
                    + f_3 * pc_x[k] * osh_1383[k];

        t_1839[k] = f_3 * pc_y[k] * osh_1379[k];

        t_1840[k] = f_11 * nsh_1385[k]
                    + f_3 * pc_x[k] * osh_1385[k];

        t_1841[k] = pa_x[k] * nsi0_1841[k]
                    - f_10 * pc_x[k] * nsi1_1841[k];
    }

#pragma omp simd aligned(t_1842, t_1843, t_1844, t_1845, pa_x, pc_x, nsi0_1842, nsi0_1843, \
                         nsi0_1844, nsi0_1845, nsi1_1842, nsi1_1843, nsi1_1844, \
                         nsi1_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1842[k] = pa_x[k] * nsi0_1842[k]
                    - f_10 * pc_x[k] * nsi1_1842[k];

        t_1843[k] = pa_x[k] * nsi0_1843[k]
                    - f_10 * pc_x[k] * nsi1_1843[k];

        t_1844[k] = pa_x[k] * nsi0_1844[k]
                    - f_10 * pc_x[k] * nsi1_1844[k];

        t_1845[k] = pa_x[k] * nsi0_1845[k]
                    - f_10 * pc_x[k] * nsi1_1845[k];
    }

#pragma omp simd aligned(t_1846, t_1847, t_1848, t_1849, pa_x, pc_x, pc_y, nsi0_1847, \
                         nsi1_1847, osg0_990, osg0_991, osg1_990, osg1_991, osh_1385, \
                         osh_1386, osh_1387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1846[k] = f_3 * pc_y[k] * osh_1385[k];

        t_1847[k] = pa_x[k] * nsi0_1847[k]
                    - f_10 * pc_x[k] * nsi1_1847[k];

        t_1848[k] = f_1 * osg0_990[k]
                    - f_2 * osg1_990[k]
                    + f_3 * pc_x[k] * osh_1386[k];

        t_1849[k] = f_16 * osg0_991[k]
                    - f_17 * osg1_991[k]
                    + f_3 * pc_x[k] * osh_1387[k];
    }

#pragma omp simd aligned(t_1850, t_1851, t_1852, t_1853, pc_x, pc_z, osg0_993, osg0_995, \
                         osg1_993, osg1_995, osh_1386, osh_1387, osh_1389, \
                         osh_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1850[k] = f_3 * pc_z[k] * osh_1386[k];

        t_1851[k] = f_8 * osg0_993[k]
                    - f_9 * osg1_993[k]
                    + f_3 * pc_x[k] * osh_1389[k];

        t_1852[k] = f_3 * pc_z[k] * osh_1387[k];

        t_1853[k] = f_8 * osg0_995[k]
                    - f_9 * osg1_995[k]
                    + f_3 * pc_x[k] * osh_1391[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, t_1857, pc_x, pc_z, osg0_996, osg0_998, \
                         osg0_999, osg1_996, osg1_998, osg1_999, osh_1389, osh_1392, osh_1394, \
                         osh_1395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_6 * osg0_996[k]
                    - f_7 * osg1_996[k]
                    + f_3 * pc_x[k] * osh_1392[k];

        t_1855[k] = f_3 * pc_z[k] * osh_1389[k];

        t_1856[k] = f_6 * osg0_998[k]
                    - f_7 * osg1_998[k]
                    + f_3 * pc_x[k] * osh_1394[k];

        t_1857[k] = f_6 * osg0_999[k]
                    - f_7 * osg1_999[k]
                    + f_3 * pc_x[k] * osh_1395[k];
    }

#pragma omp simd aligned(t_1858, t_1859, t_1860, t_1861, pc_x, pc_z, osg0_1000, osg0_1002, \
                         osg0_1003, osg1_1000, osg1_1002, osg1_1003, osh_1392, osh_1396, \
                         osh_1398, osh_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1858[k] = f_4 * osg0_1000[k]
                    - f_5 * osg1_1000[k]
                    + f_3 * pc_x[k] * osh_1396[k];

        t_1859[k] = f_3 * pc_z[k] * osh_1392[k];

        t_1860[k] = f_4 * osg0_1002[k]
                    - f_5 * osg1_1002[k]
                    + f_3 * pc_x[k] * osh_1398[k];

        t_1861[k] = f_4 * osg0_1003[k]
                    - f_5 * osg1_1003[k]
                    + f_3 * pc_x[k] * osh_1399[k];
    }

#pragma omp simd aligned(t_1862, t_1863, t_1864, t_1865, t_1866, t_1867, pc_x, osg0_1004, \
                         osg1_1004, osh_1400, osh_1401, osh_1402, osh_1403, osh_1404, \
                         osh_1405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1862[k] = f_4 * osg0_1004[k]
                    - f_5 * osg1_1004[k]
                    + f_3 * pc_x[k] * osh_1400[k];

        t_1863[k] = f_3 * pc_x[k] * osh_1401[k];

        t_1864[k] = f_3 * pc_x[k] * osh_1402[k];

        t_1865[k] = f_3 * pc_x[k] * osh_1403[k];

        t_1866[k] = f_3 * pc_x[k] * osh_1404[k];

        t_1867[k] = f_3 * pc_x[k] * osh_1405[k];
    }

#pragma omp simd aligned(t_1868, t_1869, t_1870, t_1871, pc_x, pc_y, pc_z, nsh_1170, \
                         osg0_1000, osg1_1000, osh_1401, osh_1402, \
                         osh_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1868[k] = f_3 * pc_x[k] * osh_1406[k];

        t_1869[k] = f_0 * nsh_1170[k]
                    + f_1 * osg0_1000[k]
                    - f_2 * osg1_1000[k]
                    + f_3 * pc_y[k] * osh_1401[k];

        t_1870[k] = f_3 * pc_z[k] * osh_1401[k];

        t_1871[k] = f_4 * osg0_1000[k]
                    - f_5 * osg1_1000[k]
                    + f_3 * pc_z[k] * osh_1402[k];
    }

#pragma omp simd aligned(t_1872, t_1873, t_1874, t_1875, pc_y, pc_z, nsh_1175, osg0_1001, \
                         osg0_1002, osg0_1004, osg1_1001, osg1_1002, osg1_1004, osh_1403, \
                         osh_1404, osh_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1872[k] = f_6 * osg0_1001[k]
                    - f_7 * osg1_1001[k]
                    + f_3 * pc_z[k] * osh_1403[k];

        t_1873[k] = f_8 * osg0_1002[k]
                    - f_9 * osg1_1002[k]
                    + f_3 * pc_z[k] * osh_1404[k];

        t_1874[k] = f_0 * nsh_1175[k]
                    + f_3 * pc_y[k] * osh_1406[k];

        t_1875[k] = f_1 * osg0_1004[k]
                    - f_2 * osg1_1004[k]
                    + f_3 * pc_z[k] * osh_1406[k];
    }

#pragma omp simd aligned(t_1876, t_1877, t_1878, t_1879, pa_z, pc_x, pc_z, nsi0_1540, \
                         nsi0_1541, nsi0_1543, nsi1_1540, nsi1_1541, nsi1_1543, osg0_1007, \
                         osg1_1007, osh_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1876[k] = pa_z[k] * nsi0_1540[k]
                    - f_10 * pc_z[k] * nsi1_1540[k];

        t_1877[k] = pa_z[k] * nsi0_1541[k]
                    - f_10 * pc_z[k] * nsi1_1541[k];

        t_1878[k] = f_16 * osg0_1007[k]
                    - f_17 * osg1_1007[k]
                    + f_3 * pc_x[k] * osh_1409[k];

        t_1879[k] = pa_z[k] * nsi0_1543[k]
                    - f_10 * pc_z[k] * nsi1_1543[k];
    }

#pragma omp simd aligned(t_1880, t_1881, t_1882, pa_z, pc_x, pc_z, nsi0_1546, nsi1_1546, \
                         osg0_1009, osg0_1010, osg1_1009, osg1_1010, osh_1411, \
                         osh_1412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1880[k] = f_8 * osg0_1009[k]
                    - f_9 * osg1_1009[k]
                    + f_3 * pc_x[k] * osh_1411[k];

        t_1881[k] = f_8 * osg0_1010[k]
                    - f_9 * osg1_1010[k]
                    + f_3 * pc_x[k] * osh_1412[k];

        t_1882[k] = pa_z[k] * nsi0_1546[k]
                    - f_10 * pc_z[k] * nsi1_1546[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pc_x, osg0_1012, osg0_1013, osg0_1014, \
                         osg1_1012, osg1_1013, osg1_1014, osh_1414, osh_1415, \
                         osh_1416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_6 * osg0_1012[k]
                    - f_7 * osg1_1012[k]
                    + f_3 * pc_x[k] * osh_1414[k];

        t_1884[k] = f_6 * osg0_1013[k]
                    - f_7 * osg1_1013[k]
                    + f_3 * pc_x[k] * osh_1415[k];

        t_1885[k] = f_6 * osg0_1014[k]
                    - f_7 * osg1_1014[k]
                    + f_3 * pc_x[k] * osh_1416[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pa_z, pc_x, pc_z, nsi0_1550, nsi1_1550, \
                         osg0_1016, osg0_1017, osg1_1016, osg1_1017, osh_1418, \
                         osh_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = pa_z[k] * nsi0_1550[k]
                    - f_10 * pc_z[k] * nsi1_1550[k];

        t_1887[k] = f_4 * osg0_1016[k]
                    - f_5 * osg1_1016[k]
                    + f_3 * pc_x[k] * osh_1418[k];

        t_1888[k] = f_4 * osg0_1017[k]
                    - f_5 * osg1_1017[k]
                    + f_3 * pc_x[k] * osh_1419[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, t_1892, t_1893, pc_x, osg0_1018, osg0_1019, \
                         osg1_1018, osg1_1019, osh_1420, osh_1421, osh_1422, osh_1423, \
                         osh_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = f_4 * osg0_1018[k]
                    - f_5 * osg1_1018[k]
                    + f_3 * pc_x[k] * osh_1420[k];

        t_1890[k] = f_4 * osg0_1019[k]
                    - f_5 * osg1_1019[k]
                    + f_3 * pc_x[k] * osh_1421[k];

        t_1891[k] = f_3 * pc_x[k] * osh_1422[k];

        t_1892[k] = f_3 * pc_x[k] * osh_1423[k];

        t_1893[k] = f_3 * pc_x[k] * osh_1424[k];
    }

#pragma omp simd aligned(t_1894, t_1895, t_1896, t_1897, t_1898, pa_z, pc_x, pc_z, nsi0_1561, \
                         nsh_1170, nsi1_1561, osh_1422, osh_1425, osh_1426, \
                         osh_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1894[k] = f_3 * pc_x[k] * osh_1425[k];

        t_1895[k] = f_3 * pc_x[k] * osh_1426[k];

        t_1896[k] = f_3 * pc_x[k] * osh_1427[k];

        t_1897[k] = pa_z[k] * nsi0_1561[k]
                    - f_10 * pc_z[k] * nsi1_1561[k];

        t_1898[k] = f_11 * nsh_1170[k]
                    + f_3 * pc_z[k] * osh_1422[k];
    }

#pragma omp simd aligned(t_1899, t_1900, t_1901, pa_z, pc_z, nsi0_1563, nsi0_1564, nsi0_1565, \
                         nsh_1171, nsh_1172, nsh_1173, nsi1_1563, nsi1_1564, \
                         nsi1_1565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1899[k] = pa_z[k] * nsi0_1563[k]
                    + f_12 * nsh_1171[k]
                    - f_10 * pc_z[k] * nsi1_1563[k];

        t_1900[k] = pa_z[k] * nsi0_1564[k]
                    + f_13 * nsh_1172[k]
                    - f_10 * pc_z[k] * nsi1_1564[k];

        t_1901[k] = pa_z[k] * nsi0_1565[k]
                    + f_14 * nsh_1173[k]
                    - f_10 * pc_z[k] * nsi1_1565[k];
    }

#pragma omp simd aligned(t_1902, t_1903, t_1904, pc_x, pc_y, pc_z, nsh_1175, nsh_1196, \
                         osg0_1019, osg0_1020, osg1_1019, osg1_1020, osh_1427, \
                         osh_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1902[k] = f_15 * nsh_1196[k]
                    + f_3 * pc_y[k] * osh_1427[k];

        t_1903[k] = f_11 * nsh_1175[k]
                    + f_1 * osg0_1019[k]
                    - f_2 * osg1_1019[k]
                    + f_3 * pc_z[k] * osh_1427[k];

        t_1904[k] = f_1 * osg0_1020[k]
                    - f_2 * osg1_1020[k]
                    + f_3 * pc_x[k] * osh_1428[k];
    }

#pragma omp simd aligned(t_1905, t_1906, t_1907, pc_x, osg0_1021, osg0_1022, osg0_1023, \
                         osg1_1021, osg1_1022, osg1_1023, osh_1429, osh_1430, \
                         osh_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1905[k] = f_16 * osg0_1021[k]
                    - f_17 * osg1_1021[k]
                    + f_3 * pc_x[k] * osh_1429[k];

        t_1906[k] = f_16 * osg0_1022[k]
                    - f_17 * osg1_1022[k]
                    + f_3 * pc_x[k] * osh_1430[k];

        t_1907[k] = f_8 * osg0_1023[k]
                    - f_9 * osg1_1023[k]
                    + f_3 * pc_x[k] * osh_1431[k];
    }

#pragma omp simd aligned(t_1908, t_1909, t_1910, pc_x, osg0_1024, osg0_1025, osg0_1026, \
                         osg1_1024, osg1_1025, osg1_1026, osh_1432, osh_1433, \
                         osh_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1908[k] = f_8 * osg0_1024[k]
                    - f_9 * osg1_1024[k]
                    + f_3 * pc_x[k] * osh_1432[k];

        t_1909[k] = f_8 * osg0_1025[k]
                    - f_9 * osg1_1025[k]
                    + f_3 * pc_x[k] * osh_1433[k];

        t_1910[k] = f_6 * osg0_1026[k]
                    - f_7 * osg1_1026[k]
                    + f_3 * pc_x[k] * osh_1434[k];
    }

#pragma omp simd aligned(t_1911, t_1912, t_1913, pc_x, osg0_1027, osg0_1028, osg0_1029, \
                         osg1_1027, osg1_1028, osg1_1029, osh_1435, osh_1436, \
                         osh_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1911[k] = f_6 * osg0_1027[k]
                    - f_7 * osg1_1027[k]
                    + f_3 * pc_x[k] * osh_1435[k];

        t_1912[k] = f_6 * osg0_1028[k]
                    - f_7 * osg1_1028[k]
                    + f_3 * pc_x[k] * osh_1436[k];

        t_1913[k] = f_6 * osg0_1029[k]
                    - f_7 * osg1_1029[k]
                    + f_3 * pc_x[k] * osh_1437[k];
    }

#pragma omp simd aligned(t_1914, t_1915, t_1916, pc_x, osg0_1030, osg0_1031, osg0_1032, \
                         osg1_1030, osg1_1031, osg1_1032, osh_1438, osh_1439, \
                         osh_1440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1914[k] = f_4 * osg0_1030[k]
                    - f_5 * osg1_1030[k]
                    + f_3 * pc_x[k] * osh_1438[k];

        t_1915[k] = f_4 * osg0_1031[k]
                    - f_5 * osg1_1031[k]
                    + f_3 * pc_x[k] * osh_1439[k];

        t_1916[k] = f_4 * osg0_1032[k]
                    - f_5 * osg1_1032[k]
                    + f_3 * pc_x[k] * osh_1440[k];
    }

#pragma omp simd aligned(t_1917, t_1918, t_1919, t_1920, t_1921, pc_x, osg0_1033, osg0_1034, \
                         osg1_1033, osg1_1034, osh_1441, osh_1442, osh_1443, osh_1444, \
                         osh_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1917[k] = f_4 * osg0_1033[k]
                    - f_5 * osg1_1033[k]
                    + f_3 * pc_x[k] * osh_1441[k];

        t_1918[k] = f_4 * osg0_1034[k]
                    - f_5 * osg1_1034[k]
                    + f_3 * pc_x[k] * osh_1442[k];

        t_1919[k] = f_3 * pc_x[k] * osh_1443[k];

        t_1920[k] = f_3 * pc_x[k] * osh_1444[k];

        t_1921[k] = f_3 * pc_x[k] * osh_1445[k];
    }

#pragma omp simd aligned(t_1922, t_1923, t_1924, t_1925, t_1926, pc_x, pc_y, pc_z, nsh_1191, \
                         nsh_1212, osg0_1030, osg1_1030, osh_1443, osh_1446, osh_1447, \
                         osh_1448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1922[k] = f_3 * pc_x[k] * osh_1446[k];

        t_1923[k] = f_3 * pc_x[k] * osh_1447[k];

        t_1924[k] = f_3 * pc_x[k] * osh_1448[k];

        t_1925[k] = f_18 * nsh_1212[k]
                    + f_1 * osg0_1030[k]
                    - f_2 * osg1_1030[k]
                    + f_3 * pc_y[k] * osh_1443[k];

        t_1926[k] = f_12 * nsh_1191[k]
                    + f_3 * pc_z[k] * osh_1443[k];
    }

#pragma omp simd aligned(t_1927, t_1928, t_1929, pc_y, nsh_1214, nsh_1215, nsh_1216, \
                         osg0_1032, osg0_1033, osg0_1034, osg1_1032, osg1_1033, osg1_1034, \
                         osh_1445, osh_1446, osh_1447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1927[k] = f_18 * nsh_1214[k]
                    + f_8 * osg0_1032[k]
                    - f_9 * osg1_1032[k]
                    + f_3 * pc_y[k] * osh_1445[k];

        t_1928[k] = f_18 * nsh_1215[k]
                    + f_6 * osg0_1033[k]
                    - f_7 * osg1_1033[k]
                    + f_3 * pc_y[k] * osh_1446[k];

        t_1929[k] = f_18 * nsh_1216[k]
                    + f_4 * osg0_1034[k]
                    - f_5 * osg1_1034[k]
                    + f_3 * pc_y[k] * osh_1447[k];
    }

#pragma omp simd aligned(t_1930, t_1931, t_1932, pc_x, pc_y, pc_z, nsh_1196, nsh_1217, \
                         osg0_1034, osg0_1035, osg1_1034, osg1_1035, osh_1448, \
                         osh_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1930[k] = f_18 * nsh_1217[k]
                    + f_3 * pc_y[k] * osh_1448[k];

        t_1931[k] = f_12 * nsh_1196[k]
                    + f_1 * osg0_1034[k]
                    - f_2 * osg1_1034[k]
                    + f_3 * pc_z[k] * osh_1448[k];

        t_1932[k] = f_1 * osg0_1035[k]
                    - f_2 * osg1_1035[k]
                    + f_3 * pc_x[k] * osh_1449[k];
    }

#pragma omp simd aligned(t_1933, t_1934, t_1935, pc_x, osg0_1036, osg0_1037, osg0_1038, \
                         osg1_1036, osg1_1037, osg1_1038, osh_1450, osh_1451, \
                         osh_1452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1933[k] = f_16 * osg0_1036[k]
                    - f_17 * osg1_1036[k]
                    + f_3 * pc_x[k] * osh_1450[k];

        t_1934[k] = f_16 * osg0_1037[k]
                    - f_17 * osg1_1037[k]
                    + f_3 * pc_x[k] * osh_1451[k];

        t_1935[k] = f_8 * osg0_1038[k]
                    - f_9 * osg1_1038[k]
                    + f_3 * pc_x[k] * osh_1452[k];
    }

#pragma omp simd aligned(t_1936, t_1937, t_1938, pc_x, osg0_1039, osg0_1040, osg0_1041, \
                         osg1_1039, osg1_1040, osg1_1041, osh_1453, osh_1454, \
                         osh_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1936[k] = f_8 * osg0_1039[k]
                    - f_9 * osg1_1039[k]
                    + f_3 * pc_x[k] * osh_1453[k];

        t_1937[k] = f_8 * osg0_1040[k]
                    - f_9 * osg1_1040[k]
                    + f_3 * pc_x[k] * osh_1454[k];

        t_1938[k] = f_6 * osg0_1041[k]
                    - f_7 * osg1_1041[k]
                    + f_3 * pc_x[k] * osh_1455[k];
    }

#pragma omp simd aligned(t_1939, t_1940, t_1941, pc_x, osg0_1042, osg0_1043, osg0_1044, \
                         osg1_1042, osg1_1043, osg1_1044, osh_1456, osh_1457, \
                         osh_1458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1939[k] = f_6 * osg0_1042[k]
                    - f_7 * osg1_1042[k]
                    + f_3 * pc_x[k] * osh_1456[k];

        t_1940[k] = f_6 * osg0_1043[k]
                    - f_7 * osg1_1043[k]
                    + f_3 * pc_x[k] * osh_1457[k];

        t_1941[k] = f_6 * osg0_1044[k]
                    - f_7 * osg1_1044[k]
                    + f_3 * pc_x[k] * osh_1458[k];
    }

#pragma omp simd aligned(t_1942, t_1943, t_1944, pc_x, osg0_1045, osg0_1046, osg0_1047, \
                         osg1_1045, osg1_1046, osg1_1047, osh_1459, osh_1460, \
                         osh_1461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1942[k] = f_4 * osg0_1045[k]
                    - f_5 * osg1_1045[k]
                    + f_3 * pc_x[k] * osh_1459[k];

        t_1943[k] = f_4 * osg0_1046[k]
                    - f_5 * osg1_1046[k]
                    + f_3 * pc_x[k] * osh_1460[k];

        t_1944[k] = f_4 * osg0_1047[k]
                    - f_5 * osg1_1047[k]
                    + f_3 * pc_x[k] * osh_1461[k];
    }

#pragma omp simd aligned(t_1945, t_1946, t_1947, t_1948, t_1949, pc_x, osg0_1048, osg0_1049, \
                         osg1_1048, osg1_1049, osh_1462, osh_1463, osh_1464, osh_1465, \
                         osh_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1945[k] = f_4 * osg0_1048[k]
                    - f_5 * osg1_1048[k]
                    + f_3 * pc_x[k] * osh_1462[k];

        t_1946[k] = f_4 * osg0_1049[k]
                    - f_5 * osg1_1049[k]
                    + f_3 * pc_x[k] * osh_1463[k];

        t_1947[k] = f_3 * pc_x[k] * osh_1464[k];

        t_1948[k] = f_3 * pc_x[k] * osh_1465[k];

        t_1949[k] = f_3 * pc_x[k] * osh_1466[k];
    }

#pragma omp simd aligned(t_1950, t_1951, t_1952, t_1953, t_1954, pc_x, pc_y, pc_z, nsh_1212, \
                         nsh_1233, osg0_1045, osg1_1045, osh_1464, osh_1467, osh_1468, \
                         osh_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1950[k] = f_3 * pc_x[k] * osh_1467[k];

        t_1951[k] = f_3 * pc_x[k] * osh_1468[k];

        t_1952[k] = f_3 * pc_x[k] * osh_1469[k];

        t_1953[k] = f_19 * nsh_1233[k]
                    + f_1 * osg0_1045[k]
                    - f_2 * osg1_1045[k]
                    + f_3 * pc_y[k] * osh_1464[k];

        t_1954[k] = f_13 * nsh_1212[k]
                    + f_3 * pc_z[k] * osh_1464[k];
    }

#pragma omp simd aligned(t_1955, t_1956, t_1957, pc_y, nsh_1235, nsh_1236, nsh_1237, \
                         osg0_1047, osg0_1048, osg0_1049, osg1_1047, osg1_1048, osg1_1049, \
                         osh_1466, osh_1467, osh_1468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1955[k] = f_19 * nsh_1235[k]
                    + f_8 * osg0_1047[k]
                    - f_9 * osg1_1047[k]
                    + f_3 * pc_y[k] * osh_1466[k];

        t_1956[k] = f_19 * nsh_1236[k]
                    + f_6 * osg0_1048[k]
                    - f_7 * osg1_1048[k]
                    + f_3 * pc_y[k] * osh_1467[k];

        t_1957[k] = f_19 * nsh_1237[k]
                    + f_4 * osg0_1049[k]
                    - f_5 * osg1_1049[k]
                    + f_3 * pc_y[k] * osh_1468[k];
    }

#pragma omp simd aligned(t_1958, t_1959, t_1960, pc_x, pc_y, pc_z, nsh_1217, nsh_1238, \
                         osg0_1049, osg0_1050, osg1_1049, osg1_1050, osh_1469, \
                         osh_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1958[k] = f_19 * nsh_1238[k]
                    + f_3 * pc_y[k] * osh_1469[k];

        t_1959[k] = f_13 * nsh_1217[k]
                    + f_1 * osg0_1049[k]
                    - f_2 * osg1_1049[k]
                    + f_3 * pc_z[k] * osh_1469[k];

        t_1960[k] = f_1 * osg0_1050[k]
                    - f_2 * osg1_1050[k]
                    + f_3 * pc_x[k] * osh_1470[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t nsh, const size_t osg0,
                                                           const size_t osg1, const size_t osh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;
    const auto f_22 = 2.5 / q;

    auto *t_1961 = buffer.data(target + 1961);
    auto *t_1962 = buffer.data(target + 1962);
    auto *t_1963 = buffer.data(target + 1963);
    auto *t_1964 = buffer.data(target + 1964);
    auto *t_1965 = buffer.data(target + 1965);
    auto *t_1966 = buffer.data(target + 1966);
    auto *t_1967 = buffer.data(target + 1967);
    auto *t_1968 = buffer.data(target + 1968);
    auto *t_1969 = buffer.data(target + 1969);
    auto *t_1970 = buffer.data(target + 1970);
    auto *t_1971 = buffer.data(target + 1971);
    auto *t_1972 = buffer.data(target + 1972);
    auto *t_1973 = buffer.data(target + 1973);
    auto *t_1974 = buffer.data(target + 1974);
    auto *t_1975 = buffer.data(target + 1975);
    auto *t_1976 = buffer.data(target + 1976);
    auto *t_1977 = buffer.data(target + 1977);
    auto *t_1978 = buffer.data(target + 1978);
    auto *t_1979 = buffer.data(target + 1979);
    auto *t_1980 = buffer.data(target + 1980);
    auto *t_1981 = buffer.data(target + 1981);
    auto *t_1982 = buffer.data(target + 1982);
    auto *t_1983 = buffer.data(target + 1983);
    auto *t_1984 = buffer.data(target + 1984);
    auto *t_1985 = buffer.data(target + 1985);
    auto *t_1986 = buffer.data(target + 1986);
    auto *t_1987 = buffer.data(target + 1987);
    auto *t_1988 = buffer.data(target + 1988);
    auto *t_1989 = buffer.data(target + 1989);
    auto *t_1990 = buffer.data(target + 1990);
    auto *t_1991 = buffer.data(target + 1991);
    auto *t_1992 = buffer.data(target + 1992);
    auto *t_1993 = buffer.data(target + 1993);
    auto *t_1994 = buffer.data(target + 1994);
    auto *t_1995 = buffer.data(target + 1995);
    auto *t_1996 = buffer.data(target + 1996);
    auto *t_1997 = buffer.data(target + 1997);
    auto *t_1998 = buffer.data(target + 1998);
    auto *t_1999 = buffer.data(target + 1999);
    auto *t_2000 = buffer.data(target + 2000);
    auto *t_2001 = buffer.data(target + 2001);
    auto *t_2002 = buffer.data(target + 2002);
    auto *t_2003 = buffer.data(target + 2003);
    auto *t_2004 = buffer.data(target + 2004);
    auto *t_2005 = buffer.data(target + 2005);
    auto *t_2006 = buffer.data(target + 2006);
    auto *t_2007 = buffer.data(target + 2007);
    auto *t_2008 = buffer.data(target + 2008);
    auto *t_2009 = buffer.data(target + 2009);
    auto *t_2010 = buffer.data(target + 2010);
    auto *t_2011 = buffer.data(target + 2011);
    auto *t_2012 = buffer.data(target + 2012);
    auto *t_2013 = buffer.data(target + 2013);
    auto *t_2014 = buffer.data(target + 2014);
    auto *t_2015 = buffer.data(target + 2015);
    auto *t_2016 = buffer.data(target + 2016);
    auto *t_2017 = buffer.data(target + 2017);
    auto *t_2018 = buffer.data(target + 2018);
    auto *t_2019 = buffer.data(target + 2019);
    auto *t_2020 = buffer.data(target + 2020);
    auto *t_2021 = buffer.data(target + 2021);
    auto *t_2022 = buffer.data(target + 2022);
    auto *t_2023 = buffer.data(target + 2023);
    auto *t_2024 = buffer.data(target + 2024);
    auto *t_2025 = buffer.data(target + 2025);
    auto *t_2026 = buffer.data(target + 2026);
    auto *t_2027 = buffer.data(target + 2027);
    auto *t_2028 = buffer.data(target + 2028);
    auto *t_2029 = buffer.data(target + 2029);
    auto *t_2030 = buffer.data(target + 2030);
    auto *t_2031 = buffer.data(target + 2031);
    auto *t_2032 = buffer.data(target + 2032);
    auto *t_2033 = buffer.data(target + 2033);
    auto *t_2034 = buffer.data(target + 2034);
    auto *t_2035 = buffer.data(target + 2035);
    auto *t_2036 = buffer.data(target + 2036);
    auto *t_2037 = buffer.data(target + 2037);
    auto *t_2038 = buffer.data(target + 2038);
    auto *t_2039 = buffer.data(target + 2039);
    auto *t_2040 = buffer.data(target + 2040);
    auto *t_2041 = buffer.data(target + 2041);
    auto *t_2042 = buffer.data(target + 2042);
    auto *t_2043 = buffer.data(target + 2043);
    auto *t_2044 = buffer.data(target + 2044);
    auto *t_2045 = buffer.data(target + 2045);
    auto *t_2046 = buffer.data(target + 2046);
    auto *t_2047 = buffer.data(target + 2047);
    auto *t_2048 = buffer.data(target + 2048);
    auto *t_2049 = buffer.data(target + 2049);
    auto *t_2050 = buffer.data(target + 2050);
    auto *t_2051 = buffer.data(target + 2051);
    auto *t_2052 = buffer.data(target + 2052);
    auto *t_2053 = buffer.data(target + 2053);
    auto *t_2054 = buffer.data(target + 2054);
    auto *t_2055 = buffer.data(target + 2055);
    auto *t_2056 = buffer.data(target + 2056);
    auto *t_2057 = buffer.data(target + 2057);
    auto *t_2058 = buffer.data(target + 2058);
    auto *t_2059 = buffer.data(target + 2059);
    auto *t_2060 = buffer.data(target + 2060);
    auto *t_2061 = buffer.data(target + 2061);
    auto *t_2062 = buffer.data(target + 2062);
    auto *t_2063 = buffer.data(target + 2063);
    auto *t_2064 = buffer.data(target + 2064);
    auto *t_2065 = buffer.data(target + 2065);
    auto *t_2066 = buffer.data(target + 2066);
    auto *t_2067 = buffer.data(target + 2067);
    auto *t_2068 = buffer.data(target + 2068);
    auto *t_2069 = buffer.data(target + 2069);
    auto *t_2070 = buffer.data(target + 2070);
    auto *t_2071 = buffer.data(target + 2071);
    auto *t_2072 = buffer.data(target + 2072);
    auto *t_2073 = buffer.data(target + 2073);
    auto *t_2074 = buffer.data(target + 2074);
    auto *t_2075 = buffer.data(target + 2075);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh_1233 = buffer.data(nsh + 1233);
    const auto *nsh_1238 = buffer.data(nsh + 1238);
    const auto *nsh_1254 = buffer.data(nsh + 1254);
    const auto *nsh_1256 = buffer.data(nsh + 1256);
    const auto *nsh_1257 = buffer.data(nsh + 1257);
    const auto *nsh_1258 = buffer.data(nsh + 1258);
    const auto *nsh_1259 = buffer.data(nsh + 1259);
    const auto *nsh_1275 = buffer.data(nsh + 1275);
    const auto *nsh_1277 = buffer.data(nsh + 1277);
    const auto *nsh_1278 = buffer.data(nsh + 1278);
    const auto *nsh_1279 = buffer.data(nsh + 1279);
    const auto *nsh_1280 = buffer.data(nsh + 1280);
    const auto *nsh_1296 = buffer.data(nsh + 1296);
    const auto *nsh_1298 = buffer.data(nsh + 1298);
    const auto *nsh_1299 = buffer.data(nsh + 1299);
    const auto *nsh_1300 = buffer.data(nsh + 1300);
    const auto *nsh_1301 = buffer.data(nsh + 1301);
    const auto *nsh_1317 = buffer.data(nsh + 1317);
    const auto *nsh_1319 = buffer.data(nsh + 1319);
    const auto *nsh_1320 = buffer.data(nsh + 1320);
    const auto *nsh_1321 = buffer.data(nsh + 1321);
    const auto *nsh_1322 = buffer.data(nsh + 1322);

    const auto *osg0_1051 = buffer.data(osg0 + 1051);
    const auto *osg0_1052 = buffer.data(osg0 + 1052);
    const auto *osg0_1053 = buffer.data(osg0 + 1053);
    const auto *osg0_1054 = buffer.data(osg0 + 1054);
    const auto *osg0_1055 = buffer.data(osg0 + 1055);
    const auto *osg0_1056 = buffer.data(osg0 + 1056);
    const auto *osg0_1057 = buffer.data(osg0 + 1057);
    const auto *osg0_1058 = buffer.data(osg0 + 1058);
    const auto *osg0_1059 = buffer.data(osg0 + 1059);
    const auto *osg0_1060 = buffer.data(osg0 + 1060);
    const auto *osg0_1061 = buffer.data(osg0 + 1061);
    const auto *osg0_1062 = buffer.data(osg0 + 1062);
    const auto *osg0_1063 = buffer.data(osg0 + 1063);
    const auto *osg0_1064 = buffer.data(osg0 + 1064);
    const auto *osg0_1065 = buffer.data(osg0 + 1065);
    const auto *osg0_1066 = buffer.data(osg0 + 1066);
    const auto *osg0_1067 = buffer.data(osg0 + 1067);
    const auto *osg0_1068 = buffer.data(osg0 + 1068);
    const auto *osg0_1069 = buffer.data(osg0 + 1069);
    const auto *osg0_1070 = buffer.data(osg0 + 1070);
    const auto *osg0_1071 = buffer.data(osg0 + 1071);
    const auto *osg0_1072 = buffer.data(osg0 + 1072);
    const auto *osg0_1073 = buffer.data(osg0 + 1073);
    const auto *osg0_1074 = buffer.data(osg0 + 1074);
    const auto *osg0_1075 = buffer.data(osg0 + 1075);
    const auto *osg0_1076 = buffer.data(osg0 + 1076);
    const auto *osg0_1077 = buffer.data(osg0 + 1077);
    const auto *osg0_1078 = buffer.data(osg0 + 1078);
    const auto *osg0_1079 = buffer.data(osg0 + 1079);
    const auto *osg0_1080 = buffer.data(osg0 + 1080);
    const auto *osg0_1081 = buffer.data(osg0 + 1081);
    const auto *osg0_1082 = buffer.data(osg0 + 1082);
    const auto *osg0_1083 = buffer.data(osg0 + 1083);
    const auto *osg0_1084 = buffer.data(osg0 + 1084);
    const auto *osg0_1085 = buffer.data(osg0 + 1085);
    const auto *osg0_1086 = buffer.data(osg0 + 1086);
    const auto *osg0_1087 = buffer.data(osg0 + 1087);
    const auto *osg0_1088 = buffer.data(osg0 + 1088);
    const auto *osg0_1089 = buffer.data(osg0 + 1089);
    const auto *osg0_1090 = buffer.data(osg0 + 1090);
    const auto *osg0_1091 = buffer.data(osg0 + 1091);
    const auto *osg0_1092 = buffer.data(osg0 + 1092);
    const auto *osg0_1093 = buffer.data(osg0 + 1093);
    const auto *osg0_1094 = buffer.data(osg0 + 1094);
    const auto *osg0_1095 = buffer.data(osg0 + 1095);
    const auto *osg0_1096 = buffer.data(osg0 + 1096);
    const auto *osg0_1097 = buffer.data(osg0 + 1097);
    const auto *osg0_1098 = buffer.data(osg0 + 1098);
    const auto *osg0_1099 = buffer.data(osg0 + 1099);
    const auto *osg0_1100 = buffer.data(osg0 + 1100);
    const auto *osg0_1101 = buffer.data(osg0 + 1101);
    const auto *osg0_1102 = buffer.data(osg0 + 1102);
    const auto *osg0_1103 = buffer.data(osg0 + 1103);
    const auto *osg0_1104 = buffer.data(osg0 + 1104);
    const auto *osg0_1105 = buffer.data(osg0 + 1105);
    const auto *osg0_1106 = buffer.data(osg0 + 1106);
    const auto *osg0_1107 = buffer.data(osg0 + 1107);
    const auto *osg0_1108 = buffer.data(osg0 + 1108);
    const auto *osg0_1109 = buffer.data(osg0 + 1109);
    const auto *osg0_1110 = buffer.data(osg0 + 1110);
    const auto *osg0_1111 = buffer.data(osg0 + 1111);
    const auto *osg0_1112 = buffer.data(osg0 + 1112);
    const auto *osg0_1113 = buffer.data(osg0 + 1113);

    const auto *osg1_1051 = buffer.data(osg1 + 1051);
    const auto *osg1_1052 = buffer.data(osg1 + 1052);
    const auto *osg1_1053 = buffer.data(osg1 + 1053);
    const auto *osg1_1054 = buffer.data(osg1 + 1054);
    const auto *osg1_1055 = buffer.data(osg1 + 1055);
    const auto *osg1_1056 = buffer.data(osg1 + 1056);
    const auto *osg1_1057 = buffer.data(osg1 + 1057);
    const auto *osg1_1058 = buffer.data(osg1 + 1058);
    const auto *osg1_1059 = buffer.data(osg1 + 1059);
    const auto *osg1_1060 = buffer.data(osg1 + 1060);
    const auto *osg1_1061 = buffer.data(osg1 + 1061);
    const auto *osg1_1062 = buffer.data(osg1 + 1062);
    const auto *osg1_1063 = buffer.data(osg1 + 1063);
    const auto *osg1_1064 = buffer.data(osg1 + 1064);
    const auto *osg1_1065 = buffer.data(osg1 + 1065);
    const auto *osg1_1066 = buffer.data(osg1 + 1066);
    const auto *osg1_1067 = buffer.data(osg1 + 1067);
    const auto *osg1_1068 = buffer.data(osg1 + 1068);
    const auto *osg1_1069 = buffer.data(osg1 + 1069);
    const auto *osg1_1070 = buffer.data(osg1 + 1070);
    const auto *osg1_1071 = buffer.data(osg1 + 1071);
    const auto *osg1_1072 = buffer.data(osg1 + 1072);
    const auto *osg1_1073 = buffer.data(osg1 + 1073);
    const auto *osg1_1074 = buffer.data(osg1 + 1074);
    const auto *osg1_1075 = buffer.data(osg1 + 1075);
    const auto *osg1_1076 = buffer.data(osg1 + 1076);
    const auto *osg1_1077 = buffer.data(osg1 + 1077);
    const auto *osg1_1078 = buffer.data(osg1 + 1078);
    const auto *osg1_1079 = buffer.data(osg1 + 1079);
    const auto *osg1_1080 = buffer.data(osg1 + 1080);
    const auto *osg1_1081 = buffer.data(osg1 + 1081);
    const auto *osg1_1082 = buffer.data(osg1 + 1082);
    const auto *osg1_1083 = buffer.data(osg1 + 1083);
    const auto *osg1_1084 = buffer.data(osg1 + 1084);
    const auto *osg1_1085 = buffer.data(osg1 + 1085);
    const auto *osg1_1086 = buffer.data(osg1 + 1086);
    const auto *osg1_1087 = buffer.data(osg1 + 1087);
    const auto *osg1_1088 = buffer.data(osg1 + 1088);
    const auto *osg1_1089 = buffer.data(osg1 + 1089);
    const auto *osg1_1090 = buffer.data(osg1 + 1090);
    const auto *osg1_1091 = buffer.data(osg1 + 1091);
    const auto *osg1_1092 = buffer.data(osg1 + 1092);
    const auto *osg1_1093 = buffer.data(osg1 + 1093);
    const auto *osg1_1094 = buffer.data(osg1 + 1094);
    const auto *osg1_1095 = buffer.data(osg1 + 1095);
    const auto *osg1_1096 = buffer.data(osg1 + 1096);
    const auto *osg1_1097 = buffer.data(osg1 + 1097);
    const auto *osg1_1098 = buffer.data(osg1 + 1098);
    const auto *osg1_1099 = buffer.data(osg1 + 1099);
    const auto *osg1_1100 = buffer.data(osg1 + 1100);
    const auto *osg1_1101 = buffer.data(osg1 + 1101);
    const auto *osg1_1102 = buffer.data(osg1 + 1102);
    const auto *osg1_1103 = buffer.data(osg1 + 1103);
    const auto *osg1_1104 = buffer.data(osg1 + 1104);
    const auto *osg1_1105 = buffer.data(osg1 + 1105);
    const auto *osg1_1106 = buffer.data(osg1 + 1106);
    const auto *osg1_1107 = buffer.data(osg1 + 1107);
    const auto *osg1_1108 = buffer.data(osg1 + 1108);
    const auto *osg1_1109 = buffer.data(osg1 + 1109);
    const auto *osg1_1110 = buffer.data(osg1 + 1110);
    const auto *osg1_1111 = buffer.data(osg1 + 1111);
    const auto *osg1_1112 = buffer.data(osg1 + 1112);
    const auto *osg1_1113 = buffer.data(osg1 + 1113);

    const auto *osh_1471 = buffer.data(osh + 1471);
    const auto *osh_1472 = buffer.data(osh + 1472);
    const auto *osh_1473 = buffer.data(osh + 1473);
    const auto *osh_1474 = buffer.data(osh + 1474);
    const auto *osh_1475 = buffer.data(osh + 1475);
    const auto *osh_1476 = buffer.data(osh + 1476);
    const auto *osh_1477 = buffer.data(osh + 1477);
    const auto *osh_1478 = buffer.data(osh + 1478);
    const auto *osh_1479 = buffer.data(osh + 1479);
    const auto *osh_1480 = buffer.data(osh + 1480);
    const auto *osh_1481 = buffer.data(osh + 1481);
    const auto *osh_1482 = buffer.data(osh + 1482);
    const auto *osh_1483 = buffer.data(osh + 1483);
    const auto *osh_1484 = buffer.data(osh + 1484);
    const auto *osh_1485 = buffer.data(osh + 1485);
    const auto *osh_1486 = buffer.data(osh + 1486);
    const auto *osh_1487 = buffer.data(osh + 1487);
    const auto *osh_1488 = buffer.data(osh + 1488);
    const auto *osh_1489 = buffer.data(osh + 1489);
    const auto *osh_1490 = buffer.data(osh + 1490);
    const auto *osh_1491 = buffer.data(osh + 1491);
    const auto *osh_1492 = buffer.data(osh + 1492);
    const auto *osh_1493 = buffer.data(osh + 1493);
    const auto *osh_1494 = buffer.data(osh + 1494);
    const auto *osh_1495 = buffer.data(osh + 1495);
    const auto *osh_1496 = buffer.data(osh + 1496);
    const auto *osh_1497 = buffer.data(osh + 1497);
    const auto *osh_1498 = buffer.data(osh + 1498);
    const auto *osh_1499 = buffer.data(osh + 1499);
    const auto *osh_1500 = buffer.data(osh + 1500);
    const auto *osh_1501 = buffer.data(osh + 1501);
    const auto *osh_1502 = buffer.data(osh + 1502);
    const auto *osh_1503 = buffer.data(osh + 1503);
    const auto *osh_1504 = buffer.data(osh + 1504);
    const auto *osh_1505 = buffer.data(osh + 1505);
    const auto *osh_1506 = buffer.data(osh + 1506);
    const auto *osh_1507 = buffer.data(osh + 1507);
    const auto *osh_1508 = buffer.data(osh + 1508);
    const auto *osh_1509 = buffer.data(osh + 1509);
    const auto *osh_1510 = buffer.data(osh + 1510);
    const auto *osh_1511 = buffer.data(osh + 1511);
    const auto *osh_1512 = buffer.data(osh + 1512);
    const auto *osh_1513 = buffer.data(osh + 1513);
    const auto *osh_1514 = buffer.data(osh + 1514);
    const auto *osh_1515 = buffer.data(osh + 1515);
    const auto *osh_1516 = buffer.data(osh + 1516);
    const auto *osh_1517 = buffer.data(osh + 1517);
    const auto *osh_1518 = buffer.data(osh + 1518);
    const auto *osh_1519 = buffer.data(osh + 1519);
    const auto *osh_1520 = buffer.data(osh + 1520);
    const auto *osh_1521 = buffer.data(osh + 1521);
    const auto *osh_1522 = buffer.data(osh + 1522);
    const auto *osh_1523 = buffer.data(osh + 1523);
    const auto *osh_1524 = buffer.data(osh + 1524);
    const auto *osh_1525 = buffer.data(osh + 1525);
    const auto *osh_1526 = buffer.data(osh + 1526);
    const auto *osh_1527 = buffer.data(osh + 1527);
    const auto *osh_1528 = buffer.data(osh + 1528);
    const auto *osh_1529 = buffer.data(osh + 1529);
    const auto *osh_1530 = buffer.data(osh + 1530);
    const auto *osh_1531 = buffer.data(osh + 1531);
    const auto *osh_1532 = buffer.data(osh + 1532);
    const auto *osh_1533 = buffer.data(osh + 1533);
    const auto *osh_1534 = buffer.data(osh + 1534);
    const auto *osh_1535 = buffer.data(osh + 1535);
    const auto *osh_1536 = buffer.data(osh + 1536);
    const auto *osh_1537 = buffer.data(osh + 1537);
    const auto *osh_1538 = buffer.data(osh + 1538);
    const auto *osh_1539 = buffer.data(osh + 1539);
    const auto *osh_1540 = buffer.data(osh + 1540);
    const auto *osh_1541 = buffer.data(osh + 1541);
    const auto *osh_1542 = buffer.data(osh + 1542);
    const auto *osh_1543 = buffer.data(osh + 1543);
    const auto *osh_1544 = buffer.data(osh + 1544);
    const auto *osh_1545 = buffer.data(osh + 1545);
    const auto *osh_1546 = buffer.data(osh + 1546);
    const auto *osh_1547 = buffer.data(osh + 1547);
    const auto *osh_1548 = buffer.data(osh + 1548);
    const auto *osh_1549 = buffer.data(osh + 1549);
    const auto *osh_1550 = buffer.data(osh + 1550);
    const auto *osh_1551 = buffer.data(osh + 1551);
    const auto *osh_1552 = buffer.data(osh + 1552);
    const auto *osh_1553 = buffer.data(osh + 1553);
    const auto *osh_1554 = buffer.data(osh + 1554);
    const auto *osh_1555 = buffer.data(osh + 1555);
    const auto *osh_1556 = buffer.data(osh + 1556);
    const auto *osh_1557 = buffer.data(osh + 1557);

#pragma omp simd aligned(t_1961, t_1962, t_1963, pc_x, osg0_1051, osg0_1052, osg0_1053, \
                         osg1_1051, osg1_1052, osg1_1053, osh_1471, osh_1472, \
                         osh_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1961[k] = f_16 * osg0_1051[k]
                    - f_17 * osg1_1051[k]
                    + f_3 * pc_x[k] * osh_1471[k];

        t_1962[k] = f_16 * osg0_1052[k]
                    - f_17 * osg1_1052[k]
                    + f_3 * pc_x[k] * osh_1472[k];

        t_1963[k] = f_8 * osg0_1053[k]
                    - f_9 * osg1_1053[k]
                    + f_3 * pc_x[k] * osh_1473[k];
    }

#pragma omp simd aligned(t_1964, t_1965, t_1966, pc_x, osg0_1054, osg0_1055, osg0_1056, \
                         osg1_1054, osg1_1055, osg1_1056, osh_1474, osh_1475, \
                         osh_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1964[k] = f_8 * osg0_1054[k]
                    - f_9 * osg1_1054[k]
                    + f_3 * pc_x[k] * osh_1474[k];

        t_1965[k] = f_8 * osg0_1055[k]
                    - f_9 * osg1_1055[k]
                    + f_3 * pc_x[k] * osh_1475[k];

        t_1966[k] = f_6 * osg0_1056[k]
                    - f_7 * osg1_1056[k]
                    + f_3 * pc_x[k] * osh_1476[k];
    }

#pragma omp simd aligned(t_1967, t_1968, t_1969, pc_x, osg0_1057, osg0_1058, osg0_1059, \
                         osg1_1057, osg1_1058, osg1_1059, osh_1477, osh_1478, \
                         osh_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1967[k] = f_6 * osg0_1057[k]
                    - f_7 * osg1_1057[k]
                    + f_3 * pc_x[k] * osh_1477[k];

        t_1968[k] = f_6 * osg0_1058[k]
                    - f_7 * osg1_1058[k]
                    + f_3 * pc_x[k] * osh_1478[k];

        t_1969[k] = f_6 * osg0_1059[k]
                    - f_7 * osg1_1059[k]
                    + f_3 * pc_x[k] * osh_1479[k];
    }

#pragma omp simd aligned(t_1970, t_1971, t_1972, pc_x, osg0_1060, osg0_1061, osg0_1062, \
                         osg1_1060, osg1_1061, osg1_1062, osh_1480, osh_1481, \
                         osh_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1970[k] = f_4 * osg0_1060[k]
                    - f_5 * osg1_1060[k]
                    + f_3 * pc_x[k] * osh_1480[k];

        t_1971[k] = f_4 * osg0_1061[k]
                    - f_5 * osg1_1061[k]
                    + f_3 * pc_x[k] * osh_1481[k];

        t_1972[k] = f_4 * osg0_1062[k]
                    - f_5 * osg1_1062[k]
                    + f_3 * pc_x[k] * osh_1482[k];
    }

#pragma omp simd aligned(t_1973, t_1974, t_1975, t_1976, t_1977, pc_x, osg0_1063, osg0_1064, \
                         osg1_1063, osg1_1064, osh_1483, osh_1484, osh_1485, osh_1486, \
                         osh_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1973[k] = f_4 * osg0_1063[k]
                    - f_5 * osg1_1063[k]
                    + f_3 * pc_x[k] * osh_1483[k];

        t_1974[k] = f_4 * osg0_1064[k]
                    - f_5 * osg1_1064[k]
                    + f_3 * pc_x[k] * osh_1484[k];

        t_1975[k] = f_3 * pc_x[k] * osh_1485[k];

        t_1976[k] = f_3 * pc_x[k] * osh_1486[k];

        t_1977[k] = f_3 * pc_x[k] * osh_1487[k];
    }

#pragma omp simd aligned(t_1978, t_1979, t_1980, t_1981, t_1982, pc_x, pc_y, pc_z, nsh_1233, \
                         nsh_1254, osg0_1060, osg1_1060, osh_1485, osh_1488, osh_1489, \
                         osh_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1978[k] = f_3 * pc_x[k] * osh_1488[k];

        t_1979[k] = f_3 * pc_x[k] * osh_1489[k];

        t_1980[k] = f_3 * pc_x[k] * osh_1490[k];

        t_1981[k] = f_20 * nsh_1254[k]
                    + f_1 * osg0_1060[k]
                    - f_2 * osg1_1060[k]
                    + f_3 * pc_y[k] * osh_1485[k];

        t_1982[k] = f_14 * nsh_1233[k]
                    + f_3 * pc_z[k] * osh_1485[k];
    }

#pragma omp simd aligned(t_1983, t_1984, t_1985, pc_y, nsh_1256, nsh_1257, nsh_1258, \
                         osg0_1062, osg0_1063, osg0_1064, osg1_1062, osg1_1063, osg1_1064, \
                         osh_1487, osh_1488, osh_1489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1983[k] = f_20 * nsh_1256[k]
                    + f_8 * osg0_1062[k]
                    - f_9 * osg1_1062[k]
                    + f_3 * pc_y[k] * osh_1487[k];

        t_1984[k] = f_20 * nsh_1257[k]
                    + f_6 * osg0_1063[k]
                    - f_7 * osg1_1063[k]
                    + f_3 * pc_y[k] * osh_1488[k];

        t_1985[k] = f_20 * nsh_1258[k]
                    + f_4 * osg0_1064[k]
                    - f_5 * osg1_1064[k]
                    + f_3 * pc_y[k] * osh_1489[k];
    }

#pragma omp simd aligned(t_1986, t_1987, t_1988, pc_x, pc_y, pc_z, nsh_1238, nsh_1259, \
                         osg0_1064, osg0_1065, osg1_1064, osg1_1065, osh_1490, \
                         osh_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1986[k] = f_20 * nsh_1259[k]
                    + f_3 * pc_y[k] * osh_1490[k];

        t_1987[k] = f_14 * nsh_1238[k]
                    + f_1 * osg0_1064[k]
                    - f_2 * osg1_1064[k]
                    + f_3 * pc_z[k] * osh_1490[k];

        t_1988[k] = f_1 * osg0_1065[k]
                    - f_2 * osg1_1065[k]
                    + f_3 * pc_x[k] * osh_1491[k];
    }

#pragma omp simd aligned(t_1989, t_1990, t_1991, pc_x, osg0_1066, osg0_1067, osg0_1068, \
                         osg1_1066, osg1_1067, osg1_1068, osh_1492, osh_1493, \
                         osh_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1989[k] = f_16 * osg0_1066[k]
                    - f_17 * osg1_1066[k]
                    + f_3 * pc_x[k] * osh_1492[k];

        t_1990[k] = f_16 * osg0_1067[k]
                    - f_17 * osg1_1067[k]
                    + f_3 * pc_x[k] * osh_1493[k];

        t_1991[k] = f_8 * osg0_1068[k]
                    - f_9 * osg1_1068[k]
                    + f_3 * pc_x[k] * osh_1494[k];
    }

#pragma omp simd aligned(t_1992, t_1993, t_1994, pc_x, osg0_1069, osg0_1070, osg0_1071, \
                         osg1_1069, osg1_1070, osg1_1071, osh_1495, osh_1496, \
                         osh_1497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1992[k] = f_8 * osg0_1069[k]
                    - f_9 * osg1_1069[k]
                    + f_3 * pc_x[k] * osh_1495[k];

        t_1993[k] = f_8 * osg0_1070[k]
                    - f_9 * osg1_1070[k]
                    + f_3 * pc_x[k] * osh_1496[k];

        t_1994[k] = f_6 * osg0_1071[k]
                    - f_7 * osg1_1071[k]
                    + f_3 * pc_x[k] * osh_1497[k];
    }

#pragma omp simd aligned(t_1995, t_1996, t_1997, pc_x, osg0_1072, osg0_1073, osg0_1074, \
                         osg1_1072, osg1_1073, osg1_1074, osh_1498, osh_1499, \
                         osh_1500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1995[k] = f_6 * osg0_1072[k]
                    - f_7 * osg1_1072[k]
                    + f_3 * pc_x[k] * osh_1498[k];

        t_1996[k] = f_6 * osg0_1073[k]
                    - f_7 * osg1_1073[k]
                    + f_3 * pc_x[k] * osh_1499[k];

        t_1997[k] = f_6 * osg0_1074[k]
                    - f_7 * osg1_1074[k]
                    + f_3 * pc_x[k] * osh_1500[k];
    }

#pragma omp simd aligned(t_1998, t_1999, t_2000, pc_x, osg0_1075, osg0_1076, osg0_1077, \
                         osg1_1075, osg1_1076, osg1_1077, osh_1501, osh_1502, \
                         osh_1503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1998[k] = f_4 * osg0_1075[k]
                    - f_5 * osg1_1075[k]
                    + f_3 * pc_x[k] * osh_1501[k];

        t_1999[k] = f_4 * osg0_1076[k]
                    - f_5 * osg1_1076[k]
                    + f_3 * pc_x[k] * osh_1502[k];

        t_2000[k] = f_4 * osg0_1077[k]
                    - f_5 * osg1_1077[k]
                    + f_3 * pc_x[k] * osh_1503[k];
    }

#pragma omp simd aligned(t_2001, t_2002, t_2003, t_2004, t_2005, pc_x, osg0_1078, osg0_1079, \
                         osg1_1078, osg1_1079, osh_1504, osh_1505, osh_1506, osh_1507, \
                         osh_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2001[k] = f_4 * osg0_1078[k]
                    - f_5 * osg1_1078[k]
                    + f_3 * pc_x[k] * osh_1504[k];

        t_2002[k] = f_4 * osg0_1079[k]
                    - f_5 * osg1_1079[k]
                    + f_3 * pc_x[k] * osh_1505[k];

        t_2003[k] = f_3 * pc_x[k] * osh_1506[k];

        t_2004[k] = f_3 * pc_x[k] * osh_1507[k];

        t_2005[k] = f_3 * pc_x[k] * osh_1508[k];
    }

#pragma omp simd aligned(t_2006, t_2007, t_2008, t_2009, t_2010, pc_x, pc_y, pc_z, nsh_1254, \
                         nsh_1275, osg0_1075, osg1_1075, osh_1506, osh_1509, osh_1510, \
                         osh_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2006[k] = f_3 * pc_x[k] * osh_1509[k];

        t_2007[k] = f_3 * pc_x[k] * osh_1510[k];

        t_2008[k] = f_3 * pc_x[k] * osh_1511[k];

        t_2009[k] = f_21 * nsh_1275[k]
                    + f_1 * osg0_1075[k]
                    - f_2 * osg1_1075[k]
                    + f_3 * pc_y[k] * osh_1506[k];

        t_2010[k] = f_22 * nsh_1254[k]
                    + f_3 * pc_z[k] * osh_1506[k];
    }

#pragma omp simd aligned(t_2011, t_2012, t_2013, pc_y, nsh_1277, nsh_1278, nsh_1279, \
                         osg0_1077, osg0_1078, osg0_1079, osg1_1077, osg1_1078, osg1_1079, \
                         osh_1508, osh_1509, osh_1510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2011[k] = f_21 * nsh_1277[k]
                    + f_8 * osg0_1077[k]
                    - f_9 * osg1_1077[k]
                    + f_3 * pc_y[k] * osh_1508[k];

        t_2012[k] = f_21 * nsh_1278[k]
                    + f_6 * osg0_1078[k]
                    - f_7 * osg1_1078[k]
                    + f_3 * pc_y[k] * osh_1509[k];

        t_2013[k] = f_21 * nsh_1279[k]
                    + f_4 * osg0_1079[k]
                    - f_5 * osg1_1079[k]
                    + f_3 * pc_y[k] * osh_1510[k];
    }

#pragma omp simd aligned(t_2014, t_2015, t_2016, pc_x, pc_y, pc_z, nsh_1259, nsh_1280, \
                         osg0_1079, osg0_1080, osg1_1079, osg1_1080, osh_1511, \
                         osh_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2014[k] = f_21 * nsh_1280[k]
                    + f_3 * pc_y[k] * osh_1511[k];

        t_2015[k] = f_22 * nsh_1259[k]
                    + f_1 * osg0_1079[k]
                    - f_2 * osg1_1079[k]
                    + f_3 * pc_z[k] * osh_1511[k];

        t_2016[k] = f_1 * osg0_1080[k]
                    - f_2 * osg1_1080[k]
                    + f_3 * pc_x[k] * osh_1512[k];
    }

#pragma omp simd aligned(t_2017, t_2018, t_2019, pc_x, osg0_1081, osg0_1082, osg0_1083, \
                         osg1_1081, osg1_1082, osg1_1083, osh_1513, osh_1514, \
                         osh_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2017[k] = f_16 * osg0_1081[k]
                    - f_17 * osg1_1081[k]
                    + f_3 * pc_x[k] * osh_1513[k];

        t_2018[k] = f_16 * osg0_1082[k]
                    - f_17 * osg1_1082[k]
                    + f_3 * pc_x[k] * osh_1514[k];

        t_2019[k] = f_8 * osg0_1083[k]
                    - f_9 * osg1_1083[k]
                    + f_3 * pc_x[k] * osh_1515[k];
    }

#pragma omp simd aligned(t_2020, t_2021, t_2022, pc_x, osg0_1084, osg0_1085, osg0_1086, \
                         osg1_1084, osg1_1085, osg1_1086, osh_1516, osh_1517, \
                         osh_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2020[k] = f_8 * osg0_1084[k]
                    - f_9 * osg1_1084[k]
                    + f_3 * pc_x[k] * osh_1516[k];

        t_2021[k] = f_8 * osg0_1085[k]
                    - f_9 * osg1_1085[k]
                    + f_3 * pc_x[k] * osh_1517[k];

        t_2022[k] = f_6 * osg0_1086[k]
                    - f_7 * osg1_1086[k]
                    + f_3 * pc_x[k] * osh_1518[k];
    }

#pragma omp simd aligned(t_2023, t_2024, t_2025, pc_x, osg0_1087, osg0_1088, osg0_1089, \
                         osg1_1087, osg1_1088, osg1_1089, osh_1519, osh_1520, \
                         osh_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2023[k] = f_6 * osg0_1087[k]
                    - f_7 * osg1_1087[k]
                    + f_3 * pc_x[k] * osh_1519[k];

        t_2024[k] = f_6 * osg0_1088[k]
                    - f_7 * osg1_1088[k]
                    + f_3 * pc_x[k] * osh_1520[k];

        t_2025[k] = f_6 * osg0_1089[k]
                    - f_7 * osg1_1089[k]
                    + f_3 * pc_x[k] * osh_1521[k];
    }

#pragma omp simd aligned(t_2026, t_2027, t_2028, pc_x, osg0_1090, osg0_1091, osg0_1092, \
                         osg1_1090, osg1_1091, osg1_1092, osh_1522, osh_1523, \
                         osh_1524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2026[k] = f_4 * osg0_1090[k]
                    - f_5 * osg1_1090[k]
                    + f_3 * pc_x[k] * osh_1522[k];

        t_2027[k] = f_4 * osg0_1091[k]
                    - f_5 * osg1_1091[k]
                    + f_3 * pc_x[k] * osh_1523[k];

        t_2028[k] = f_4 * osg0_1092[k]
                    - f_5 * osg1_1092[k]
                    + f_3 * pc_x[k] * osh_1524[k];
    }

#pragma omp simd aligned(t_2029, t_2030, t_2031, t_2032, t_2033, pc_x, osg0_1093, osg0_1094, \
                         osg1_1093, osg1_1094, osh_1525, osh_1526, osh_1527, osh_1528, \
                         osh_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2029[k] = f_4 * osg0_1093[k]
                    - f_5 * osg1_1093[k]
                    + f_3 * pc_x[k] * osh_1525[k];

        t_2030[k] = f_4 * osg0_1094[k]
                    - f_5 * osg1_1094[k]
                    + f_3 * pc_x[k] * osh_1526[k];

        t_2031[k] = f_3 * pc_x[k] * osh_1527[k];

        t_2032[k] = f_3 * pc_x[k] * osh_1528[k];

        t_2033[k] = f_3 * pc_x[k] * osh_1529[k];
    }

#pragma omp simd aligned(t_2034, t_2035, t_2036, t_2037, t_2038, pc_x, pc_y, pc_z, nsh_1275, \
                         nsh_1296, osg0_1090, osg1_1090, osh_1527, osh_1530, osh_1531, \
                         osh_1532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2034[k] = f_3 * pc_x[k] * osh_1530[k];

        t_2035[k] = f_3 * pc_x[k] * osh_1531[k];

        t_2036[k] = f_3 * pc_x[k] * osh_1532[k];

        t_2037[k] = f_22 * nsh_1296[k]
                    + f_1 * osg0_1090[k]
                    - f_2 * osg1_1090[k]
                    + f_3 * pc_y[k] * osh_1527[k];

        t_2038[k] = f_21 * nsh_1275[k]
                    + f_3 * pc_z[k] * osh_1527[k];
    }

#pragma omp simd aligned(t_2039, t_2040, t_2041, pc_y, nsh_1298, nsh_1299, nsh_1300, \
                         osg0_1092, osg0_1093, osg0_1094, osg1_1092, osg1_1093, osg1_1094, \
                         osh_1529, osh_1530, osh_1531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2039[k] = f_22 * nsh_1298[k]
                    + f_8 * osg0_1092[k]
                    - f_9 * osg1_1092[k]
                    + f_3 * pc_y[k] * osh_1529[k];

        t_2040[k] = f_22 * nsh_1299[k]
                    + f_6 * osg0_1093[k]
                    - f_7 * osg1_1093[k]
                    + f_3 * pc_y[k] * osh_1530[k];

        t_2041[k] = f_22 * nsh_1300[k]
                    + f_4 * osg0_1094[k]
                    - f_5 * osg1_1094[k]
                    + f_3 * pc_y[k] * osh_1531[k];
    }

#pragma omp simd aligned(t_2042, t_2043, t_2044, pc_x, pc_y, pc_z, nsh_1280, nsh_1301, \
                         osg0_1094, osg0_1095, osg1_1094, osg1_1095, osh_1532, \
                         osh_1533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2042[k] = f_22 * nsh_1301[k]
                    + f_3 * pc_y[k] * osh_1532[k];

        t_2043[k] = f_21 * nsh_1280[k]
                    + f_1 * osg0_1094[k]
                    - f_2 * osg1_1094[k]
                    + f_3 * pc_z[k] * osh_1532[k];

        t_2044[k] = f_1 * osg0_1095[k]
                    - f_2 * osg1_1095[k]
                    + f_3 * pc_x[k] * osh_1533[k];
    }

#pragma omp simd aligned(t_2045, t_2046, t_2047, pc_x, osg0_1096, osg0_1097, osg0_1098, \
                         osg1_1096, osg1_1097, osg1_1098, osh_1534, osh_1535, \
                         osh_1536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2045[k] = f_16 * osg0_1096[k]
                    - f_17 * osg1_1096[k]
                    + f_3 * pc_x[k] * osh_1534[k];

        t_2046[k] = f_16 * osg0_1097[k]
                    - f_17 * osg1_1097[k]
                    + f_3 * pc_x[k] * osh_1535[k];

        t_2047[k] = f_8 * osg0_1098[k]
                    - f_9 * osg1_1098[k]
                    + f_3 * pc_x[k] * osh_1536[k];
    }

#pragma omp simd aligned(t_2048, t_2049, t_2050, pc_x, osg0_1099, osg0_1100, osg0_1101, \
                         osg1_1099, osg1_1100, osg1_1101, osh_1537, osh_1538, \
                         osh_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2048[k] = f_8 * osg0_1099[k]
                    - f_9 * osg1_1099[k]
                    + f_3 * pc_x[k] * osh_1537[k];

        t_2049[k] = f_8 * osg0_1100[k]
                    - f_9 * osg1_1100[k]
                    + f_3 * pc_x[k] * osh_1538[k];

        t_2050[k] = f_6 * osg0_1101[k]
                    - f_7 * osg1_1101[k]
                    + f_3 * pc_x[k] * osh_1539[k];
    }

#pragma omp simd aligned(t_2051, t_2052, t_2053, pc_x, osg0_1102, osg0_1103, osg0_1104, \
                         osg1_1102, osg1_1103, osg1_1104, osh_1540, osh_1541, \
                         osh_1542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2051[k] = f_6 * osg0_1102[k]
                    - f_7 * osg1_1102[k]
                    + f_3 * pc_x[k] * osh_1540[k];

        t_2052[k] = f_6 * osg0_1103[k]
                    - f_7 * osg1_1103[k]
                    + f_3 * pc_x[k] * osh_1541[k];

        t_2053[k] = f_6 * osg0_1104[k]
                    - f_7 * osg1_1104[k]
                    + f_3 * pc_x[k] * osh_1542[k];
    }

#pragma omp simd aligned(t_2054, t_2055, t_2056, pc_x, osg0_1105, osg0_1106, osg0_1107, \
                         osg1_1105, osg1_1106, osg1_1107, osh_1543, osh_1544, \
                         osh_1545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2054[k] = f_4 * osg0_1105[k]
                    - f_5 * osg1_1105[k]
                    + f_3 * pc_x[k] * osh_1543[k];

        t_2055[k] = f_4 * osg0_1106[k]
                    - f_5 * osg1_1106[k]
                    + f_3 * pc_x[k] * osh_1544[k];

        t_2056[k] = f_4 * osg0_1107[k]
                    - f_5 * osg1_1107[k]
                    + f_3 * pc_x[k] * osh_1545[k];
    }

#pragma omp simd aligned(t_2057, t_2058, t_2059, t_2060, t_2061, pc_x, osg0_1108, osg0_1109, \
                         osg1_1108, osg1_1109, osh_1546, osh_1547, osh_1548, osh_1549, \
                         osh_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2057[k] = f_4 * osg0_1108[k]
                    - f_5 * osg1_1108[k]
                    + f_3 * pc_x[k] * osh_1546[k];

        t_2058[k] = f_4 * osg0_1109[k]
                    - f_5 * osg1_1109[k]
                    + f_3 * pc_x[k] * osh_1547[k];

        t_2059[k] = f_3 * pc_x[k] * osh_1548[k];

        t_2060[k] = f_3 * pc_x[k] * osh_1549[k];

        t_2061[k] = f_3 * pc_x[k] * osh_1550[k];
    }

#pragma omp simd aligned(t_2062, t_2063, t_2064, t_2065, t_2066, pc_x, pc_y, pc_z, nsh_1296, \
                         nsh_1317, osg0_1105, osg1_1105, osh_1548, osh_1551, osh_1552, \
                         osh_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2062[k] = f_3 * pc_x[k] * osh_1551[k];

        t_2063[k] = f_3 * pc_x[k] * osh_1552[k];

        t_2064[k] = f_3 * pc_x[k] * osh_1553[k];

        t_2065[k] = f_14 * nsh_1317[k]
                    + f_1 * osg0_1105[k]
                    - f_2 * osg1_1105[k]
                    + f_3 * pc_y[k] * osh_1548[k];

        t_2066[k] = f_20 * nsh_1296[k]
                    + f_3 * pc_z[k] * osh_1548[k];
    }

#pragma omp simd aligned(t_2067, t_2068, t_2069, pc_y, nsh_1319, nsh_1320, nsh_1321, \
                         osg0_1107, osg0_1108, osg0_1109, osg1_1107, osg1_1108, osg1_1109, \
                         osh_1550, osh_1551, osh_1552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2067[k] = f_14 * nsh_1319[k]
                    + f_8 * osg0_1107[k]
                    - f_9 * osg1_1107[k]
                    + f_3 * pc_y[k] * osh_1550[k];

        t_2068[k] = f_14 * nsh_1320[k]
                    + f_6 * osg0_1108[k]
                    - f_7 * osg1_1108[k]
                    + f_3 * pc_y[k] * osh_1551[k];

        t_2069[k] = f_14 * nsh_1321[k]
                    + f_4 * osg0_1109[k]
                    - f_5 * osg1_1109[k]
                    + f_3 * pc_y[k] * osh_1552[k];
    }

#pragma omp simd aligned(t_2070, t_2071, t_2072, pc_x, pc_y, pc_z, nsh_1301, nsh_1322, \
                         osg0_1109, osg0_1110, osg1_1109, osg1_1110, osh_1553, \
                         osh_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2070[k] = f_14 * nsh_1322[k]
                    + f_3 * pc_y[k] * osh_1553[k];

        t_2071[k] = f_20 * nsh_1301[k]
                    + f_1 * osg0_1109[k]
                    - f_2 * osg1_1109[k]
                    + f_3 * pc_z[k] * osh_1553[k];

        t_2072[k] = f_1 * osg0_1110[k]
                    - f_2 * osg1_1110[k]
                    + f_3 * pc_x[k] * osh_1554[k];
    }

#pragma omp simd aligned(t_2073, t_2074, t_2075, pc_x, osg0_1111, osg0_1112, osg0_1113, \
                         osg1_1111, osg1_1112, osg1_1113, osh_1555, osh_1556, \
                         osh_1557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2073[k] = f_16 * osg0_1111[k]
                    - f_17 * osg1_1111[k]
                    + f_3 * pc_x[k] * osh_1555[k];

        t_2074[k] = f_16 * osg0_1112[k]
                    - f_17 * osg1_1112[k]
                    + f_3 * pc_x[k] * osh_1556[k];

        t_2075[k] = f_8 * osg0_1113[k]
                    - f_9 * osg1_1113[k]
                    + f_3 * pc_x[k] * osh_1557[k];
    }
}

static auto
compute_prim_osi_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsi0,
                                                           const size_t nsh, const size_t nsi1,
                                                           const size_t osg0, const size_t osg1,
                                                           const size_t osh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_21 = 3.0 / q;

    auto *t_2076 = buffer.data(target + 2076);
    auto *t_2077 = buffer.data(target + 2077);
    auto *t_2078 = buffer.data(target + 2078);
    auto *t_2079 = buffer.data(target + 2079);
    auto *t_2080 = buffer.data(target + 2080);
    auto *t_2081 = buffer.data(target + 2081);
    auto *t_2082 = buffer.data(target + 2082);
    auto *t_2083 = buffer.data(target + 2083);
    auto *t_2084 = buffer.data(target + 2084);
    auto *t_2085 = buffer.data(target + 2085);
    auto *t_2086 = buffer.data(target + 2086);
    auto *t_2087 = buffer.data(target + 2087);
    auto *t_2088 = buffer.data(target + 2088);
    auto *t_2089 = buffer.data(target + 2089);
    auto *t_2090 = buffer.data(target + 2090);
    auto *t_2091 = buffer.data(target + 2091);
    auto *t_2092 = buffer.data(target + 2092);
    auto *t_2093 = buffer.data(target + 2093);
    auto *t_2094 = buffer.data(target + 2094);
    auto *t_2095 = buffer.data(target + 2095);
    auto *t_2096 = buffer.data(target + 2096);
    auto *t_2097 = buffer.data(target + 2097);
    auto *t_2098 = buffer.data(target + 2098);
    auto *t_2099 = buffer.data(target + 2099);
    auto *t_2100 = buffer.data(target + 2100);
    auto *t_2101 = buffer.data(target + 2101);
    auto *t_2102 = buffer.data(target + 2102);
    auto *t_2103 = buffer.data(target + 2103);
    auto *t_2104 = buffer.data(target + 2104);
    auto *t_2105 = buffer.data(target + 2105);
    auto *t_2106 = buffer.data(target + 2106);
    auto *t_2107 = buffer.data(target + 2107);
    auto *t_2108 = buffer.data(target + 2108);
    auto *t_2109 = buffer.data(target + 2109);
    auto *t_2110 = buffer.data(target + 2110);
    auto *t_2111 = buffer.data(target + 2111);
    auto *t_2112 = buffer.data(target + 2112);
    auto *t_2113 = buffer.data(target + 2113);
    auto *t_2114 = buffer.data(target + 2114);
    auto *t_2115 = buffer.data(target + 2115);
    auto *t_2116 = buffer.data(target + 2116);
    auto *t_2117 = buffer.data(target + 2117);
    auto *t_2118 = buffer.data(target + 2118);
    auto *t_2119 = buffer.data(target + 2119);
    auto *t_2120 = buffer.data(target + 2120);
    auto *t_2121 = buffer.data(target + 2121);
    auto *t_2122 = buffer.data(target + 2122);
    auto *t_2123 = buffer.data(target + 2123);
    auto *t_2124 = buffer.data(target + 2124);
    auto *t_2125 = buffer.data(target + 2125);
    auto *t_2126 = buffer.data(target + 2126);
    auto *t_2127 = buffer.data(target + 2127);
    auto *t_2128 = buffer.data(target + 2128);
    auto *t_2129 = buffer.data(target + 2129);
    auto *t_2130 = buffer.data(target + 2130);
    auto *t_2131 = buffer.data(target + 2131);
    auto *t_2132 = buffer.data(target + 2132);
    auto *t_2133 = buffer.data(target + 2133);
    auto *t_2134 = buffer.data(target + 2134);
    auto *t_2135 = buffer.data(target + 2135);
    auto *t_2136 = buffer.data(target + 2136);
    auto *t_2137 = buffer.data(target + 2137);
    auto *t_2138 = buffer.data(target + 2138);
    auto *t_2139 = buffer.data(target + 2139);
    auto *t_2140 = buffer.data(target + 2140);
    auto *t_2141 = buffer.data(target + 2141);
    auto *t_2142 = buffer.data(target + 2142);
    auto *t_2143 = buffer.data(target + 2143);
    auto *t_2144 = buffer.data(target + 2144);
    auto *t_2145 = buffer.data(target + 2145);
    auto *t_2146 = buffer.data(target + 2146);
    auto *t_2147 = buffer.data(target + 2147);
    auto *t_2148 = buffer.data(target + 2148);
    auto *t_2149 = buffer.data(target + 2149);
    auto *t_2150 = buffer.data(target + 2150);
    auto *t_2151 = buffer.data(target + 2151);
    auto *t_2152 = buffer.data(target + 2152);
    auto *t_2153 = buffer.data(target + 2153);
    auto *t_2154 = buffer.data(target + 2154);
    auto *t_2155 = buffer.data(target + 2155);
    auto *t_2156 = buffer.data(target + 2156);
    auto *t_2157 = buffer.data(target + 2157);
    auto *t_2158 = buffer.data(target + 2158);
    auto *t_2159 = buffer.data(target + 2159);
    auto *t_2160 = buffer.data(target + 2160);
    auto *t_2161 = buffer.data(target + 2161);
    auto *t_2162 = buffer.data(target + 2162);
    auto *t_2163 = buffer.data(target + 2163);
    auto *t_2164 = buffer.data(target + 2164);
    auto *t_2165 = buffer.data(target + 2165);
    auto *t_2166 = buffer.data(target + 2166);
    auto *t_2167 = buffer.data(target + 2167);
    auto *t_2168 = buffer.data(target + 2168);
    auto *t_2169 = buffer.data(target + 2169);
    auto *t_2170 = buffer.data(target + 2170);
    auto *t_2171 = buffer.data(target + 2171);
    auto *t_2172 = buffer.data(target + 2172);
    auto *t_2173 = buffer.data(target + 2173);
    auto *t_2174 = buffer.data(target + 2174);
    auto *t_2175 = buffer.data(target + 2175);
    auto *t_2176 = buffer.data(target + 2176);
    auto *t_2177 = buffer.data(target + 2177);
    auto *t_2178 = buffer.data(target + 2178);
    auto *t_2179 = buffer.data(target + 2179);
    auto *t_2180 = buffer.data(target + 2180);
    auto *t_2181 = buffer.data(target + 2181);
    auto *t_2182 = buffer.data(target + 2182);
    auto *t_2183 = buffer.data(target + 2183);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsi0_1820 = buffer.data(nsi0 + 1820);
    const auto *nsi0_1822 = buffer.data(nsi0 + 1822);
    const auto *nsi0_1825 = buffer.data(nsi0 + 1825);
    const auto *nsi0_1829 = buffer.data(nsi0 + 1829);
    const auto *nsi0_1834 = buffer.data(nsi0 + 1834);
    const auto *nsi0_1841 = buffer.data(nsi0 + 1841);
    const auto *nsi0_1843 = buffer.data(nsi0 + 1843);
    const auto *nsi0_1844 = buffer.data(nsi0 + 1844);
    const auto *nsi0_1845 = buffer.data(nsi0 + 1845);
    const auto *nsi0_1847 = buffer.data(nsi0 + 1847);

    const auto *nsh_1317 = buffer.data(nsh + 1317);
    const auto *nsh_1322 = buffer.data(nsh + 1322);
    const auto *nsh_1338 = buffer.data(nsh + 1338);
    const auto *nsh_1340 = buffer.data(nsh + 1340);
    const auto *nsh_1341 = buffer.data(nsh + 1341);
    const auto *nsh_1342 = buffer.data(nsh + 1342);
    const auto *nsh_1343 = buffer.data(nsh + 1343);
    const auto *nsh_1359 = buffer.data(nsh + 1359);
    const auto *nsh_1361 = buffer.data(nsh + 1361);
    const auto *nsh_1362 = buffer.data(nsh + 1362);
    const auto *nsh_1363 = buffer.data(nsh + 1363);
    const auto *nsh_1364 = buffer.data(nsh + 1364);
    const auto *nsh_1380 = buffer.data(nsh + 1380);
    const auto *nsh_1382 = buffer.data(nsh + 1382);
    const auto *nsh_1383 = buffer.data(nsh + 1383);
    const auto *nsh_1384 = buffer.data(nsh + 1384);
    const auto *nsh_1385 = buffer.data(nsh + 1385);

    const auto *nsi1_1820 = buffer.data(nsi1 + 1820);
    const auto *nsi1_1822 = buffer.data(nsi1 + 1822);
    const auto *nsi1_1825 = buffer.data(nsi1 + 1825);
    const auto *nsi1_1829 = buffer.data(nsi1 + 1829);
    const auto *nsi1_1834 = buffer.data(nsi1 + 1834);
    const auto *nsi1_1841 = buffer.data(nsi1 + 1841);
    const auto *nsi1_1843 = buffer.data(nsi1 + 1843);
    const auto *nsi1_1844 = buffer.data(nsi1 + 1844);
    const auto *nsi1_1845 = buffer.data(nsi1 + 1845);
    const auto *nsi1_1847 = buffer.data(nsi1 + 1847);

    const auto *osg0_1114 = buffer.data(osg0 + 1114);
    const auto *osg0_1115 = buffer.data(osg0 + 1115);
    const auto *osg0_1116 = buffer.data(osg0 + 1116);
    const auto *osg0_1117 = buffer.data(osg0 + 1117);
    const auto *osg0_1118 = buffer.data(osg0 + 1118);
    const auto *osg0_1119 = buffer.data(osg0 + 1119);
    const auto *osg0_1120 = buffer.data(osg0 + 1120);
    const auto *osg0_1121 = buffer.data(osg0 + 1121);
    const auto *osg0_1122 = buffer.data(osg0 + 1122);
    const auto *osg0_1123 = buffer.data(osg0 + 1123);
    const auto *osg0_1124 = buffer.data(osg0 + 1124);
    const auto *osg0_1125 = buffer.data(osg0 + 1125);
    const auto *osg0_1126 = buffer.data(osg0 + 1126);
    const auto *osg0_1127 = buffer.data(osg0 + 1127);
    const auto *osg0_1128 = buffer.data(osg0 + 1128);
    const auto *osg0_1129 = buffer.data(osg0 + 1129);
    const auto *osg0_1130 = buffer.data(osg0 + 1130);
    const auto *osg0_1131 = buffer.data(osg0 + 1131);
    const auto *osg0_1132 = buffer.data(osg0 + 1132);
    const auto *osg0_1133 = buffer.data(osg0 + 1133);
    const auto *osg0_1134 = buffer.data(osg0 + 1134);
    const auto *osg0_1135 = buffer.data(osg0 + 1135);
    const auto *osg0_1136 = buffer.data(osg0 + 1136);
    const auto *osg0_1137 = buffer.data(osg0 + 1137);
    const auto *osg0_1138 = buffer.data(osg0 + 1138);
    const auto *osg0_1139 = buffer.data(osg0 + 1139);
    const auto *osg0_1141 = buffer.data(osg0 + 1141);
    const auto *osg0_1143 = buffer.data(osg0 + 1143);
    const auto *osg0_1144 = buffer.data(osg0 + 1144);
    const auto *osg0_1146 = buffer.data(osg0 + 1146);
    const auto *osg0_1147 = buffer.data(osg0 + 1147);
    const auto *osg0_1148 = buffer.data(osg0 + 1148);
    const auto *osg0_1150 = buffer.data(osg0 + 1150);
    const auto *osg0_1151 = buffer.data(osg0 + 1151);
    const auto *osg0_1152 = buffer.data(osg0 + 1152);
    const auto *osg0_1153 = buffer.data(osg0 + 1153);
    const auto *osg0_1155 = buffer.data(osg0 + 1155);
    const auto *osg0_1157 = buffer.data(osg0 + 1157);
    const auto *osg0_1158 = buffer.data(osg0 + 1158);
    const auto *osg0_1160 = buffer.data(osg0 + 1160);
    const auto *osg0_1161 = buffer.data(osg0 + 1161);
    const auto *osg0_1162 = buffer.data(osg0 + 1162);
    const auto *osg0_1164 = buffer.data(osg0 + 1164);
    const auto *osg0_1165 = buffer.data(osg0 + 1165);
    const auto *osg0_1166 = buffer.data(osg0 + 1166);
    const auto *osg0_1167 = buffer.data(osg0 + 1167);
    const auto *osg0_1168 = buffer.data(osg0 + 1168);
    const auto *osg0_1169 = buffer.data(osg0 + 1169);

    const auto *osg1_1114 = buffer.data(osg1 + 1114);
    const auto *osg1_1115 = buffer.data(osg1 + 1115);
    const auto *osg1_1116 = buffer.data(osg1 + 1116);
    const auto *osg1_1117 = buffer.data(osg1 + 1117);
    const auto *osg1_1118 = buffer.data(osg1 + 1118);
    const auto *osg1_1119 = buffer.data(osg1 + 1119);
    const auto *osg1_1120 = buffer.data(osg1 + 1120);
    const auto *osg1_1121 = buffer.data(osg1 + 1121);
    const auto *osg1_1122 = buffer.data(osg1 + 1122);
    const auto *osg1_1123 = buffer.data(osg1 + 1123);
    const auto *osg1_1124 = buffer.data(osg1 + 1124);
    const auto *osg1_1125 = buffer.data(osg1 + 1125);
    const auto *osg1_1126 = buffer.data(osg1 + 1126);
    const auto *osg1_1127 = buffer.data(osg1 + 1127);
    const auto *osg1_1128 = buffer.data(osg1 + 1128);
    const auto *osg1_1129 = buffer.data(osg1 + 1129);
    const auto *osg1_1130 = buffer.data(osg1 + 1130);
    const auto *osg1_1131 = buffer.data(osg1 + 1131);
    const auto *osg1_1132 = buffer.data(osg1 + 1132);
    const auto *osg1_1133 = buffer.data(osg1 + 1133);
    const auto *osg1_1134 = buffer.data(osg1 + 1134);
    const auto *osg1_1135 = buffer.data(osg1 + 1135);
    const auto *osg1_1136 = buffer.data(osg1 + 1136);
    const auto *osg1_1137 = buffer.data(osg1 + 1137);
    const auto *osg1_1138 = buffer.data(osg1 + 1138);
    const auto *osg1_1139 = buffer.data(osg1 + 1139);
    const auto *osg1_1141 = buffer.data(osg1 + 1141);
    const auto *osg1_1143 = buffer.data(osg1 + 1143);
    const auto *osg1_1144 = buffer.data(osg1 + 1144);
    const auto *osg1_1146 = buffer.data(osg1 + 1146);
    const auto *osg1_1147 = buffer.data(osg1 + 1147);
    const auto *osg1_1148 = buffer.data(osg1 + 1148);
    const auto *osg1_1150 = buffer.data(osg1 + 1150);
    const auto *osg1_1151 = buffer.data(osg1 + 1151);
    const auto *osg1_1152 = buffer.data(osg1 + 1152);
    const auto *osg1_1153 = buffer.data(osg1 + 1153);
    const auto *osg1_1155 = buffer.data(osg1 + 1155);
    const auto *osg1_1157 = buffer.data(osg1 + 1157);
    const auto *osg1_1158 = buffer.data(osg1 + 1158);
    const auto *osg1_1160 = buffer.data(osg1 + 1160);
    const auto *osg1_1161 = buffer.data(osg1 + 1161);
    const auto *osg1_1162 = buffer.data(osg1 + 1162);
    const auto *osg1_1164 = buffer.data(osg1 + 1164);
    const auto *osg1_1165 = buffer.data(osg1 + 1165);
    const auto *osg1_1166 = buffer.data(osg1 + 1166);
    const auto *osg1_1167 = buffer.data(osg1 + 1167);
    const auto *osg1_1168 = buffer.data(osg1 + 1168);
    const auto *osg1_1169 = buffer.data(osg1 + 1169);

    const auto *osh_1558 = buffer.data(osh + 1558);
    const auto *osh_1559 = buffer.data(osh + 1559);
    const auto *osh_1560 = buffer.data(osh + 1560);
    const auto *osh_1561 = buffer.data(osh + 1561);
    const auto *osh_1562 = buffer.data(osh + 1562);
    const auto *osh_1563 = buffer.data(osh + 1563);
    const auto *osh_1564 = buffer.data(osh + 1564);
    const auto *osh_1565 = buffer.data(osh + 1565);
    const auto *osh_1566 = buffer.data(osh + 1566);
    const auto *osh_1567 = buffer.data(osh + 1567);
    const auto *osh_1568 = buffer.data(osh + 1568);
    const auto *osh_1569 = buffer.data(osh + 1569);
    const auto *osh_1570 = buffer.data(osh + 1570);
    const auto *osh_1571 = buffer.data(osh + 1571);
    const auto *osh_1572 = buffer.data(osh + 1572);
    const auto *osh_1573 = buffer.data(osh + 1573);
    const auto *osh_1574 = buffer.data(osh + 1574);
    const auto *osh_1575 = buffer.data(osh + 1575);
    const auto *osh_1576 = buffer.data(osh + 1576);
    const auto *osh_1577 = buffer.data(osh + 1577);
    const auto *osh_1578 = buffer.data(osh + 1578);
    const auto *osh_1579 = buffer.data(osh + 1579);
    const auto *osh_1580 = buffer.data(osh + 1580);
    const auto *osh_1581 = buffer.data(osh + 1581);
    const auto *osh_1582 = buffer.data(osh + 1582);
    const auto *osh_1583 = buffer.data(osh + 1583);
    const auto *osh_1584 = buffer.data(osh + 1584);
    const auto *osh_1585 = buffer.data(osh + 1585);
    const auto *osh_1586 = buffer.data(osh + 1586);
    const auto *osh_1587 = buffer.data(osh + 1587);
    const auto *osh_1588 = buffer.data(osh + 1588);
    const auto *osh_1589 = buffer.data(osh + 1589);
    const auto *osh_1590 = buffer.data(osh + 1590);
    const auto *osh_1591 = buffer.data(osh + 1591);
    const auto *osh_1592 = buffer.data(osh + 1592);
    const auto *osh_1593 = buffer.data(osh + 1593);
    const auto *osh_1594 = buffer.data(osh + 1594);
    const auto *osh_1595 = buffer.data(osh + 1595);
    const auto *osh_1597 = buffer.data(osh + 1597);
    const auto *osh_1599 = buffer.data(osh + 1599);
    const auto *osh_1600 = buffer.data(osh + 1600);
    const auto *osh_1602 = buffer.data(osh + 1602);
    const auto *osh_1603 = buffer.data(osh + 1603);
    const auto *osh_1604 = buffer.data(osh + 1604);
    const auto *osh_1606 = buffer.data(osh + 1606);
    const auto *osh_1607 = buffer.data(osh + 1607);
    const auto *osh_1608 = buffer.data(osh + 1608);
    const auto *osh_1609 = buffer.data(osh + 1609);
    const auto *osh_1611 = buffer.data(osh + 1611);
    const auto *osh_1612 = buffer.data(osh + 1612);
    const auto *osh_1613 = buffer.data(osh + 1613);
    const auto *osh_1614 = buffer.data(osh + 1614);
    const auto *osh_1615 = buffer.data(osh + 1615);
    const auto *osh_1616 = buffer.data(osh + 1616);
    const auto *osh_1617 = buffer.data(osh + 1617);
    const auto *osh_1619 = buffer.data(osh + 1619);
    const auto *osh_1620 = buffer.data(osh + 1620);
    const auto *osh_1622 = buffer.data(osh + 1622);
    const auto *osh_1623 = buffer.data(osh + 1623);
    const auto *osh_1624 = buffer.data(osh + 1624);
    const auto *osh_1626 = buffer.data(osh + 1626);
    const auto *osh_1627 = buffer.data(osh + 1627);
    const auto *osh_1628 = buffer.data(osh + 1628);
    const auto *osh_1629 = buffer.data(osh + 1629);
    const auto *osh_1631 = buffer.data(osh + 1631);
    const auto *osh_1632 = buffer.data(osh + 1632);
    const auto *osh_1633 = buffer.data(osh + 1633);
    const auto *osh_1634 = buffer.data(osh + 1634);
    const auto *osh_1635 = buffer.data(osh + 1635);
    const auto *osh_1636 = buffer.data(osh + 1636);
    const auto *osh_1637 = buffer.data(osh + 1637);

#pragma omp simd aligned(t_2076, t_2077, t_2078, pc_x, osg0_1114, osg0_1115, osg0_1116, \
                         osg1_1114, osg1_1115, osg1_1116, osh_1558, osh_1559, \
                         osh_1560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2076[k] = f_8 * osg0_1114[k]
                    - f_9 * osg1_1114[k]
                    + f_3 * pc_x[k] * osh_1558[k];

        t_2077[k] = f_8 * osg0_1115[k]
                    - f_9 * osg1_1115[k]
                    + f_3 * pc_x[k] * osh_1559[k];

        t_2078[k] = f_6 * osg0_1116[k]
                    - f_7 * osg1_1116[k]
                    + f_3 * pc_x[k] * osh_1560[k];
    }

#pragma omp simd aligned(t_2079, t_2080, t_2081, pc_x, osg0_1117, osg0_1118, osg0_1119, \
                         osg1_1117, osg1_1118, osg1_1119, osh_1561, osh_1562, \
                         osh_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2079[k] = f_6 * osg0_1117[k]
                    - f_7 * osg1_1117[k]
                    + f_3 * pc_x[k] * osh_1561[k];

        t_2080[k] = f_6 * osg0_1118[k]
                    - f_7 * osg1_1118[k]
                    + f_3 * pc_x[k] * osh_1562[k];

        t_2081[k] = f_6 * osg0_1119[k]
                    - f_7 * osg1_1119[k]
                    + f_3 * pc_x[k] * osh_1563[k];
    }

#pragma omp simd aligned(t_2082, t_2083, t_2084, pc_x, osg0_1120, osg0_1121, osg0_1122, \
                         osg1_1120, osg1_1121, osg1_1122, osh_1564, osh_1565, \
                         osh_1566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2082[k] = f_4 * osg0_1120[k]
                    - f_5 * osg1_1120[k]
                    + f_3 * pc_x[k] * osh_1564[k];

        t_2083[k] = f_4 * osg0_1121[k]
                    - f_5 * osg1_1121[k]
                    + f_3 * pc_x[k] * osh_1565[k];

        t_2084[k] = f_4 * osg0_1122[k]
                    - f_5 * osg1_1122[k]
                    + f_3 * pc_x[k] * osh_1566[k];
    }

#pragma omp simd aligned(t_2085, t_2086, t_2087, t_2088, t_2089, pc_x, osg0_1123, osg0_1124, \
                         osg1_1123, osg1_1124, osh_1567, osh_1568, osh_1569, osh_1570, \
                         osh_1571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2085[k] = f_4 * osg0_1123[k]
                    - f_5 * osg1_1123[k]
                    + f_3 * pc_x[k] * osh_1567[k];

        t_2086[k] = f_4 * osg0_1124[k]
                    - f_5 * osg1_1124[k]
                    + f_3 * pc_x[k] * osh_1568[k];

        t_2087[k] = f_3 * pc_x[k] * osh_1569[k];

        t_2088[k] = f_3 * pc_x[k] * osh_1570[k];

        t_2089[k] = f_3 * pc_x[k] * osh_1571[k];
    }

#pragma omp simd aligned(t_2090, t_2091, t_2092, t_2093, t_2094, pc_x, pc_y, pc_z, nsh_1317, \
                         nsh_1338, osg0_1120, osg1_1120, osh_1569, osh_1572, osh_1573, \
                         osh_1574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2090[k] = f_3 * pc_x[k] * osh_1572[k];

        t_2091[k] = f_3 * pc_x[k] * osh_1573[k];

        t_2092[k] = f_3 * pc_x[k] * osh_1574[k];

        t_2093[k] = f_13 * nsh_1338[k]
                    + f_1 * osg0_1120[k]
                    - f_2 * osg1_1120[k]
                    + f_3 * pc_y[k] * osh_1569[k];

        t_2094[k] = f_19 * nsh_1317[k]
                    + f_3 * pc_z[k] * osh_1569[k];
    }

#pragma omp simd aligned(t_2095, t_2096, t_2097, pc_y, nsh_1340, nsh_1341, nsh_1342, \
                         osg0_1122, osg0_1123, osg0_1124, osg1_1122, osg1_1123, osg1_1124, \
                         osh_1571, osh_1572, osh_1573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2095[k] = f_13 * nsh_1340[k]
                    + f_8 * osg0_1122[k]
                    - f_9 * osg1_1122[k]
                    + f_3 * pc_y[k] * osh_1571[k];

        t_2096[k] = f_13 * nsh_1341[k]
                    + f_6 * osg0_1123[k]
                    - f_7 * osg1_1123[k]
                    + f_3 * pc_y[k] * osh_1572[k];

        t_2097[k] = f_13 * nsh_1342[k]
                    + f_4 * osg0_1124[k]
                    - f_5 * osg1_1124[k]
                    + f_3 * pc_y[k] * osh_1573[k];
    }

#pragma omp simd aligned(t_2098, t_2099, t_2100, pc_x, pc_y, pc_z, nsh_1322, nsh_1343, \
                         osg0_1124, osg0_1125, osg1_1124, osg1_1125, osh_1574, \
                         osh_1575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2098[k] = f_13 * nsh_1343[k]
                    + f_3 * pc_y[k] * osh_1574[k];

        t_2099[k] = f_19 * nsh_1322[k]
                    + f_1 * osg0_1124[k]
                    - f_2 * osg1_1124[k]
                    + f_3 * pc_z[k] * osh_1574[k];

        t_2100[k] = f_1 * osg0_1125[k]
                    - f_2 * osg1_1125[k]
                    + f_3 * pc_x[k] * osh_1575[k];
    }

#pragma omp simd aligned(t_2101, t_2102, t_2103, pc_x, osg0_1126, osg0_1127, osg0_1128, \
                         osg1_1126, osg1_1127, osg1_1128, osh_1576, osh_1577, \
                         osh_1578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2101[k] = f_16 * osg0_1126[k]
                    - f_17 * osg1_1126[k]
                    + f_3 * pc_x[k] * osh_1576[k];

        t_2102[k] = f_16 * osg0_1127[k]
                    - f_17 * osg1_1127[k]
                    + f_3 * pc_x[k] * osh_1577[k];

        t_2103[k] = f_8 * osg0_1128[k]
                    - f_9 * osg1_1128[k]
                    + f_3 * pc_x[k] * osh_1578[k];
    }

#pragma omp simd aligned(t_2104, t_2105, t_2106, pc_x, osg0_1129, osg0_1130, osg0_1131, \
                         osg1_1129, osg1_1130, osg1_1131, osh_1579, osh_1580, \
                         osh_1581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2104[k] = f_8 * osg0_1129[k]
                    - f_9 * osg1_1129[k]
                    + f_3 * pc_x[k] * osh_1579[k];

        t_2105[k] = f_8 * osg0_1130[k]
                    - f_9 * osg1_1130[k]
                    + f_3 * pc_x[k] * osh_1580[k];

        t_2106[k] = f_6 * osg0_1131[k]
                    - f_7 * osg1_1131[k]
                    + f_3 * pc_x[k] * osh_1581[k];
    }

#pragma omp simd aligned(t_2107, t_2108, t_2109, pc_x, osg0_1132, osg0_1133, osg0_1134, \
                         osg1_1132, osg1_1133, osg1_1134, osh_1582, osh_1583, \
                         osh_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2107[k] = f_6 * osg0_1132[k]
                    - f_7 * osg1_1132[k]
                    + f_3 * pc_x[k] * osh_1582[k];

        t_2108[k] = f_6 * osg0_1133[k]
                    - f_7 * osg1_1133[k]
                    + f_3 * pc_x[k] * osh_1583[k];

        t_2109[k] = f_6 * osg0_1134[k]
                    - f_7 * osg1_1134[k]
                    + f_3 * pc_x[k] * osh_1584[k];
    }

#pragma omp simd aligned(t_2110, t_2111, t_2112, pc_x, osg0_1135, osg0_1136, osg0_1137, \
                         osg1_1135, osg1_1136, osg1_1137, osh_1585, osh_1586, \
                         osh_1587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2110[k] = f_4 * osg0_1135[k]
                    - f_5 * osg1_1135[k]
                    + f_3 * pc_x[k] * osh_1585[k];

        t_2111[k] = f_4 * osg0_1136[k]
                    - f_5 * osg1_1136[k]
                    + f_3 * pc_x[k] * osh_1586[k];

        t_2112[k] = f_4 * osg0_1137[k]
                    - f_5 * osg1_1137[k]
                    + f_3 * pc_x[k] * osh_1587[k];
    }

#pragma omp simd aligned(t_2113, t_2114, t_2115, t_2116, t_2117, pc_x, osg0_1138, osg0_1139, \
                         osg1_1138, osg1_1139, osh_1588, osh_1589, osh_1590, osh_1591, \
                         osh_1592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2113[k] = f_4 * osg0_1138[k]
                    - f_5 * osg1_1138[k]
                    + f_3 * pc_x[k] * osh_1588[k];

        t_2114[k] = f_4 * osg0_1139[k]
                    - f_5 * osg1_1139[k]
                    + f_3 * pc_x[k] * osh_1589[k];

        t_2115[k] = f_3 * pc_x[k] * osh_1590[k];

        t_2116[k] = f_3 * pc_x[k] * osh_1591[k];

        t_2117[k] = f_3 * pc_x[k] * osh_1592[k];
    }

#pragma omp simd aligned(t_2118, t_2119, t_2120, t_2121, t_2122, pc_x, pc_y, pc_z, nsh_1338, \
                         nsh_1359, osg0_1135, osg1_1135, osh_1590, osh_1593, osh_1594, \
                         osh_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2118[k] = f_3 * pc_x[k] * osh_1593[k];

        t_2119[k] = f_3 * pc_x[k] * osh_1594[k];

        t_2120[k] = f_3 * pc_x[k] * osh_1595[k];

        t_2121[k] = f_12 * nsh_1359[k]
                    + f_1 * osg0_1135[k]
                    - f_2 * osg1_1135[k]
                    + f_3 * pc_y[k] * osh_1590[k];

        t_2122[k] = f_18 * nsh_1338[k]
                    + f_3 * pc_z[k] * osh_1590[k];
    }

#pragma omp simd aligned(t_2123, t_2124, t_2125, pc_y, nsh_1361, nsh_1362, nsh_1363, \
                         osg0_1137, osg0_1138, osg0_1139, osg1_1137, osg1_1138, osg1_1139, \
                         osh_1592, osh_1593, osh_1594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2123[k] = f_12 * nsh_1361[k]
                    + f_8 * osg0_1137[k]
                    - f_9 * osg1_1137[k]
                    + f_3 * pc_y[k] * osh_1592[k];

        t_2124[k] = f_12 * nsh_1362[k]
                    + f_6 * osg0_1138[k]
                    - f_7 * osg1_1138[k]
                    + f_3 * pc_y[k] * osh_1593[k];

        t_2125[k] = f_12 * nsh_1363[k]
                    + f_4 * osg0_1139[k]
                    - f_5 * osg1_1139[k]
                    + f_3 * pc_y[k] * osh_1594[k];
    }

#pragma omp simd aligned(t_2126, t_2127, t_2128, pa_y, pc_y, pc_z, nsi0_1820, nsh_1343, \
                         nsh_1364, nsi1_1820, osg0_1139, osg1_1139, \
                         osh_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2126[k] = f_12 * nsh_1364[k]
                    + f_3 * pc_y[k] * osh_1595[k];

        t_2127[k] = f_18 * nsh_1343[k]
                    + f_1 * osg0_1139[k]
                    - f_2 * osg1_1139[k]
                    + f_3 * pc_z[k] * osh_1595[k];

        t_2128[k] = pa_y[k] * nsi0_1820[k]
                    - f_10 * pc_y[k] * nsi1_1820[k];
    }

#pragma omp simd aligned(t_2129, t_2130, t_2131, pa_y, pc_x, pc_y, nsi0_1822, nsi1_1822, \
                         osg0_1141, osg0_1143, osg1_1141, osg1_1143, osh_1597, \
                         osh_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2129[k] = f_16 * osg0_1141[k]
                    - f_17 * osg1_1141[k]
                    + f_3 * pc_x[k] * osh_1597[k];

        t_2130[k] = pa_y[k] * nsi0_1822[k]
                    - f_10 * pc_y[k] * nsi1_1822[k];

        t_2131[k] = f_8 * osg0_1143[k]
                    - f_9 * osg1_1143[k]
                    + f_3 * pc_x[k] * osh_1599[k];
    }

#pragma omp simd aligned(t_2132, t_2133, t_2134, pa_y, pc_x, pc_y, nsi0_1825, nsi1_1825, \
                         osg0_1144, osg0_1146, osg1_1144, osg1_1146, osh_1600, \
                         osh_1602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2132[k] = f_8 * osg0_1144[k]
                    - f_9 * osg1_1144[k]
                    + f_3 * pc_x[k] * osh_1600[k];

        t_2133[k] = pa_y[k] * nsi0_1825[k]
                    - f_10 * pc_y[k] * nsi1_1825[k];

        t_2134[k] = f_6 * osg0_1146[k]
                    - f_7 * osg1_1146[k]
                    + f_3 * pc_x[k] * osh_1602[k];
    }

#pragma omp simd aligned(t_2135, t_2136, t_2137, pa_y, pc_x, pc_y, nsi0_1829, nsi1_1829, \
                         osg0_1147, osg0_1148, osg1_1147, osg1_1148, osh_1603, \
                         osh_1604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2135[k] = f_6 * osg0_1147[k]
                    - f_7 * osg1_1147[k]
                    + f_3 * pc_x[k] * osh_1603[k];

        t_2136[k] = f_6 * osg0_1148[k]
                    - f_7 * osg1_1148[k]
                    + f_3 * pc_x[k] * osh_1604[k];

        t_2137[k] = pa_y[k] * nsi0_1829[k]
                    - f_10 * pc_y[k] * nsi1_1829[k];
    }

#pragma omp simd aligned(t_2138, t_2139, t_2140, pc_x, osg0_1150, osg0_1151, osg0_1152, \
                         osg1_1150, osg1_1151, osg1_1152, osh_1606, osh_1607, \
                         osh_1608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2138[k] = f_4 * osg0_1150[k]
                    - f_5 * osg1_1150[k]
                    + f_3 * pc_x[k] * osh_1606[k];

        t_2139[k] = f_4 * osg0_1151[k]
                    - f_5 * osg1_1151[k]
                    + f_3 * pc_x[k] * osh_1607[k];

        t_2140[k] = f_4 * osg0_1152[k]
                    - f_5 * osg1_1152[k]
                    + f_3 * pc_x[k] * osh_1608[k];
    }

#pragma omp simd aligned(t_2141, t_2142, t_2143, t_2144, t_2145, pa_y, pc_x, pc_y, nsi0_1834, \
                         nsi1_1834, osg0_1153, osg1_1153, osh_1609, osh_1611, osh_1612, \
                         osh_1613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2141[k] = f_4 * osg0_1153[k]
                    - f_5 * osg1_1153[k]
                    + f_3 * pc_x[k] * osh_1609[k];

        t_2142[k] = pa_y[k] * nsi0_1834[k]
                    - f_10 * pc_y[k] * nsi1_1834[k];

        t_2143[k] = f_3 * pc_x[k] * osh_1611[k];

        t_2144[k] = f_3 * pc_x[k] * osh_1612[k];

        t_2145[k] = f_3 * pc_x[k] * osh_1613[k];
    }

#pragma omp simd aligned(t_2146, t_2147, t_2148, t_2149, pa_y, pc_x, pc_y, nsi0_1841, \
                         nsh_1380, nsi1_1841, osh_1614, osh_1615, \
                         osh_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2146[k] = f_3 * pc_x[k] * osh_1614[k];

        t_2147[k] = f_3 * pc_x[k] * osh_1615[k];

        t_2148[k] = f_3 * pc_x[k] * osh_1616[k];

        t_2149[k] = pa_y[k] * nsi0_1841[k]
                    + f_21 * nsh_1380[k]
                    - f_10 * pc_y[k] * nsi1_1841[k];
    }

#pragma omp simd aligned(t_2150, t_2151, t_2152, pa_y, pc_y, pc_z, nsi0_1843, nsi0_1844, \
                         nsh_1359, nsh_1382, nsh_1383, nsi1_1843, nsi1_1844, \
                         osh_1611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2150[k] = f_15 * nsh_1359[k]
                    + f_3 * pc_z[k] * osh_1611[k];

        t_2151[k] = pa_y[k] * nsi0_1843[k]
                    + f_14 * nsh_1382[k]
                    - f_10 * pc_y[k] * nsi1_1843[k];

        t_2152[k] = pa_y[k] * nsi0_1844[k]
                    + f_13 * nsh_1383[k]
                    - f_10 * pc_y[k] * nsi1_1844[k];
    }

#pragma omp simd aligned(t_2153, t_2154, t_2155, pa_y, pc_y, nsi0_1845, nsi0_1847, nsh_1384, \
                         nsh_1385, nsi1_1845, nsi1_1847, osh_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2153[k] = pa_y[k] * nsi0_1845[k]
                    + f_12 * nsh_1384[k]
                    - f_10 * pc_y[k] * nsi1_1845[k];

        t_2154[k] = f_11 * nsh_1385[k]
                    + f_3 * pc_y[k] * osh_1616[k];

        t_2155[k] = pa_y[k] * nsi0_1847[k]
                    - f_10 * pc_y[k] * nsi1_1847[k];
    }

#pragma omp simd aligned(t_2156, t_2157, t_2158, t_2159, t_2160, pc_x, pc_y, osg0_1155, \
                         osg0_1157, osg0_1158, osg1_1155, osg1_1157, osg1_1158, osh_1617, \
                         osh_1619, osh_1620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2156[k] = f_1 * osg0_1155[k]
                    - f_2 * osg1_1155[k]
                    + f_3 * pc_x[k] * osh_1617[k];

        t_2157[k] = f_3 * pc_y[k] * osh_1617[k];

        t_2158[k] = f_16 * osg0_1157[k]
                    - f_17 * osg1_1157[k]
                    + f_3 * pc_x[k] * osh_1619[k];

        t_2159[k] = f_8 * osg0_1158[k]
                    - f_9 * osg1_1158[k]
                    + f_3 * pc_x[k] * osh_1620[k];

        t_2160[k] = f_3 * pc_y[k] * osh_1619[k];
    }

#pragma omp simd aligned(t_2161, t_2162, t_2163, t_2164, pc_x, pc_y, osg0_1160, osg0_1161, \
                         osg0_1162, osg1_1160, osg1_1161, osg1_1162, osh_1622, osh_1623, \
                         osh_1624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2161[k] = f_8 * osg0_1160[k]
                    - f_9 * osg1_1160[k]
                    + f_3 * pc_x[k] * osh_1622[k];

        t_2162[k] = f_6 * osg0_1161[k]
                    - f_7 * osg1_1161[k]
                    + f_3 * pc_x[k] * osh_1623[k];

        t_2163[k] = f_6 * osg0_1162[k]
                    - f_7 * osg1_1162[k]
                    + f_3 * pc_x[k] * osh_1624[k];

        t_2164[k] = f_3 * pc_y[k] * osh_1622[k];
    }

#pragma omp simd aligned(t_2165, t_2166, t_2167, pc_x, osg0_1164, osg0_1165, osg0_1166, \
                         osg1_1164, osg1_1165, osg1_1166, osh_1626, osh_1627, \
                         osh_1628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2165[k] = f_6 * osg0_1164[k]
                    - f_7 * osg1_1164[k]
                    + f_3 * pc_x[k] * osh_1626[k];

        t_2166[k] = f_4 * osg0_1165[k]
                    - f_5 * osg1_1165[k]
                    + f_3 * pc_x[k] * osh_1627[k];

        t_2167[k] = f_4 * osg0_1166[k]
                    - f_5 * osg1_1166[k]
                    + f_3 * pc_x[k] * osh_1628[k];
    }

#pragma omp simd aligned(t_2168, t_2169, t_2170, t_2171, t_2172, pc_x, pc_y, osg0_1167, \
                         osg0_1169, osg1_1167, osg1_1169, osh_1626, osh_1629, osh_1631, \
                         osh_1632, osh_1633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2168[k] = f_4 * osg0_1167[k]
                    - f_5 * osg1_1167[k]
                    + f_3 * pc_x[k] * osh_1629[k];

        t_2169[k] = f_3 * pc_y[k] * osh_1626[k];

        t_2170[k] = f_4 * osg0_1169[k]
                    - f_5 * osg1_1169[k]
                    + f_3 * pc_x[k] * osh_1631[k];

        t_2171[k] = f_3 * pc_x[k] * osh_1632[k];

        t_2172[k] = f_3 * pc_x[k] * osh_1633[k];
    }

#pragma omp simd aligned(t_2173, t_2174, t_2175, t_2176, t_2177, pc_x, pc_y, osg0_1165, \
                         osg1_1165, osh_1632, osh_1634, osh_1635, osh_1636, \
                         osh_1637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2173[k] = f_3 * pc_x[k] * osh_1634[k];

        t_2174[k] = f_3 * pc_x[k] * osh_1635[k];

        t_2175[k] = f_3 * pc_x[k] * osh_1636[k];

        t_2176[k] = f_3 * pc_x[k] * osh_1637[k];

        t_2177[k] = f_1 * osg0_1165[k]
                    - f_2 * osg1_1165[k]
                    + f_3 * pc_y[k] * osh_1632[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, pc_y, osg0_1166, osg0_1167, osg0_1168, \
                         osg1_1166, osg1_1167, osg1_1168, osh_1633, osh_1634, \
                         osh_1635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = f_16 * osg0_1166[k]
                    - f_17 * osg1_1166[k]
                    + f_3 * pc_y[k] * osh_1633[k];

        t_2179[k] = f_8 * osg0_1167[k]
                    - f_9 * osg1_1167[k]
                    + f_3 * pc_y[k] * osh_1634[k];

        t_2180[k] = f_6 * osg0_1168[k]
                    - f_7 * osg1_1168[k]
                    + f_3 * pc_y[k] * osh_1635[k];
    }

#pragma omp simd aligned(t_2181, t_2182, t_2183, pc_y, pc_z, nsh_1385, osg0_1169, osg1_1169, \
                         osh_1636, osh_1637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2181[k] = f_4 * osg0_1169[k]
                    - f_5 * osg1_1169[k]
                    + f_3 * pc_y[k] * osh_1636[k];

        t_2182[k] = f_3 * pc_y[k] * osh_1637[k];

        t_2183[k] = f_0 * nsh_1385[k]
                    + f_1 * osg0_1169[k]
                    - f_2 * osg1_1169[k]
                    + f_3 * pc_z[k] * osh_1637[k];
    }
}

auto
compute_prim_osi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t nsi0, const size_t nsh,
                                                   const size_t nsi1, const size_t osg0,
                                                   const size_t osg1, const size_t osh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_osi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, nsi0, nsh,
                                                              nsi1, osg0, osg1, osh, ncols,
                                                              gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osg0, osg1, osh,
                                                               ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osg0, osg1, osh,
                                                               ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece12(buffer, target, pc, nsh, osg0,
                                                               osg1, osh, ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osg0, osg1, osh,
                                                               ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece14(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osh, ncols, gamma, p,
                                                               q);

    compute_prim_osi_three_center_electron_repulsion_0_piece15(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osg0, osg1, osh,
                                                               ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osg0, osg1, osh,
                                                               ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece17(buffer, target, pc, nsh, osg0,
                                                               osg1, osh, ncols, gamma, p, q);

    compute_prim_osi_three_center_electron_repulsion_0_piece18(buffer, target, pa, pc, nsi0,
                                                               nsh, nsi1, osg0, osg1, osh,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
