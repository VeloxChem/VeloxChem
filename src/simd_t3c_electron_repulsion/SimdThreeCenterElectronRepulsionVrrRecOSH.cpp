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


#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_0 = buffer.data(nsh0 + 0);
    const auto *nsh0_3 = buffer.data(nsh0 + 3);
    const auto *nsh0_5 = buffer.data(nsh0 + 5);
    const auto *nsh0_6 = buffer.data(nsh0 + 6);
    const auto *nsh0_9 = buffer.data(nsh0 + 9);
    const auto *nsh0_15 = buffer.data(nsh0 + 15);
    const auto *nsh0_20 = buffer.data(nsh0 + 20);
    const auto *nsh0_24 = buffer.data(nsh0 + 24);
    const auto *nsh0_27 = buffer.data(nsh0 + 27);
    const auto *nsh0_36 = buffer.data(nsh0 + 36);
    const auto *nsh0_42 = buffer.data(nsh0 + 42);
    const auto *nsh0_47 = buffer.data(nsh0 + 47);
    const auto *nsh0_51 = buffer.data(nsh0 + 51);
    const auto *nsh0_62 = buffer.data(nsh0 + 62);

    const auto *nsg_0 = buffer.data(nsg + 0);
    const auto *nsg_1 = buffer.data(nsg + 1);
    const auto *nsg_2 = buffer.data(nsg + 2);
    const auto *nsg_3 = buffer.data(nsg + 3);
    const auto *nsg_5 = buffer.data(nsg + 5);
    const auto *nsg_10 = buffer.data(nsg + 10);
    const auto *nsg_12 = buffer.data(nsg + 12);
    const auto *nsg_14 = buffer.data(nsg + 14);
    const auto *nsg_15 = buffer.data(nsg + 15);
    const auto *nsg_18 = buffer.data(nsg + 18);
    const auto *nsg_20 = buffer.data(nsg + 20);
    const auto *nsg_25 = buffer.data(nsg + 25);
    const auto *nsg_27 = buffer.data(nsg + 27);
    const auto *nsg_28 = buffer.data(nsg + 28);
    const auto *nsg_29 = buffer.data(nsg + 29);
    const auto *nsg_30 = buffer.data(nsg + 30);
    const auto *nsg_32 = buffer.data(nsg + 32);
    const auto *nsg_35 = buffer.data(nsg + 35);
    const auto *nsg_40 = buffer.data(nsg + 40);
    const auto *nsg_41 = buffer.data(nsg + 41);
    const auto *nsg_42 = buffer.data(nsg + 42);
    const auto *nsg_43 = buffer.data(nsg + 43);
    const auto *nsg_44 = buffer.data(nsg + 44);
    const auto *nsg_45 = buffer.data(nsg + 45);
    const auto *nsg_48 = buffer.data(nsg + 48);
    const auto *nsg_51 = buffer.data(nsg + 51);
    const auto *nsg_55 = buffer.data(nsg + 55);
    const auto *nsg_57 = buffer.data(nsg + 57);
    const auto *nsg_58 = buffer.data(nsg + 58);
    const auto *nsg_59 = buffer.data(nsg + 59);
    const auto *nsg_70 = buffer.data(nsg + 70);
    const auto *nsg_71 = buffer.data(nsg + 71);
    const auto *nsg_72 = buffer.data(nsg + 72);
    const auto *nsg_73 = buffer.data(nsg + 73);
    const auto *nsg_74 = buffer.data(nsg + 74);
    const auto *nsg_75 = buffer.data(nsg + 75);
    const auto *nsg_80 = buffer.data(nsg + 80);
    const auto *nsg_84 = buffer.data(nsg + 84);
    const auto *nsg_85 = buffer.data(nsg + 85);
    const auto *nsg_86 = buffer.data(nsg + 86);
    const auto *nsg_87 = buffer.data(nsg + 87);
    const auto *nsg_89 = buffer.data(nsg + 89);
    const auto *nsg_90 = buffer.data(nsg + 90);
    const auto *nsg_93 = buffer.data(nsg + 93);

    const auto *nsh1_0 = buffer.data(nsh1 + 0);
    const auto *nsh1_3 = buffer.data(nsh1 + 3);
    const auto *nsh1_5 = buffer.data(nsh1 + 5);
    const auto *nsh1_6 = buffer.data(nsh1 + 6);
    const auto *nsh1_9 = buffer.data(nsh1 + 9);
    const auto *nsh1_15 = buffer.data(nsh1 + 15);
    const auto *nsh1_20 = buffer.data(nsh1 + 20);
    const auto *nsh1_24 = buffer.data(nsh1 + 24);
    const auto *nsh1_27 = buffer.data(nsh1 + 27);
    const auto *nsh1_36 = buffer.data(nsh1 + 36);
    const auto *nsh1_42 = buffer.data(nsh1 + 42);
    const auto *nsh1_47 = buffer.data(nsh1 + 47);
    const auto *nsh1_51 = buffer.data(nsh1 + 51);
    const auto *nsh1_62 = buffer.data(nsh1 + 62);

    const auto *osf0_0 = buffer.data(osf0 + 0);
    const auto *osf0_1 = buffer.data(osf0 + 1);
    const auto *osf0_2 = buffer.data(osf0 + 2);
    const auto *osf0_6 = buffer.data(osf0 + 6);
    const auto *osf0_8 = buffer.data(osf0 + 8);
    const auto *osf0_9 = buffer.data(osf0 + 9);
    const auto *osf0_16 = buffer.data(osf0 + 16);
    const auto *osf0_17 = buffer.data(osf0 + 17);
    const auto *osf0_22 = buffer.data(osf0 + 22);
    const auto *osf0_27 = buffer.data(osf0 + 27);
    const auto *osf0_28 = buffer.data(osf0 + 28);
    const auto *osf0_29 = buffer.data(osf0 + 29);
    const auto *osf0_30 = buffer.data(osf0 + 30);
    const auto *osf0_32 = buffer.data(osf0 + 32);
    const auto *osf0_33 = buffer.data(osf0 + 33);
    const auto *osf0_36 = buffer.data(osf0 + 36);
    const auto *osf0_37 = buffer.data(osf0 + 37);
    const auto *osf0_39 = buffer.data(osf0 + 39);
    const auto *osf0_48 = buffer.data(osf0 + 48);
    const auto *osf0_49 = buffer.data(osf0 + 49);
    const auto *osf0_50 = buffer.data(osf0 + 50);
    const auto *osf0_51 = buffer.data(osf0 + 51);
    const auto *osf0_52 = buffer.data(osf0 + 52);
    const auto *osf0_55 = buffer.data(osf0 + 55);
    const auto *osf0_56 = buffer.data(osf0 + 56);
    const auto *osf0_57 = buffer.data(osf0 + 57);
    const auto *osf0_58 = buffer.data(osf0 + 58);
    const auto *osf0_59 = buffer.data(osf0 + 59);
    const auto *osf0_60 = buffer.data(osf0 + 60);
    const auto *osf0_63 = buffer.data(osf0 + 63);

    const auto *osf1_0 = buffer.data(osf1 + 0);
    const auto *osf1_1 = buffer.data(osf1 + 1);
    const auto *osf1_2 = buffer.data(osf1 + 2);
    const auto *osf1_6 = buffer.data(osf1 + 6);
    const auto *osf1_8 = buffer.data(osf1 + 8);
    const auto *osf1_9 = buffer.data(osf1 + 9);
    const auto *osf1_16 = buffer.data(osf1 + 16);
    const auto *osf1_17 = buffer.data(osf1 + 17);
    const auto *osf1_22 = buffer.data(osf1 + 22);
    const auto *osf1_27 = buffer.data(osf1 + 27);
    const auto *osf1_28 = buffer.data(osf1 + 28);
    const auto *osf1_29 = buffer.data(osf1 + 29);
    const auto *osf1_30 = buffer.data(osf1 + 30);
    const auto *osf1_32 = buffer.data(osf1 + 32);
    const auto *osf1_33 = buffer.data(osf1 + 33);
    const auto *osf1_36 = buffer.data(osf1 + 36);
    const auto *osf1_37 = buffer.data(osf1 + 37);
    const auto *osf1_39 = buffer.data(osf1 + 39);
    const auto *osf1_48 = buffer.data(osf1 + 48);
    const auto *osf1_49 = buffer.data(osf1 + 49);
    const auto *osf1_50 = buffer.data(osf1 + 50);
    const auto *osf1_51 = buffer.data(osf1 + 51);
    const auto *osf1_52 = buffer.data(osf1 + 52);
    const auto *osf1_55 = buffer.data(osf1 + 55);
    const auto *osf1_56 = buffer.data(osf1 + 56);
    const auto *osf1_57 = buffer.data(osf1 + 57);
    const auto *osf1_58 = buffer.data(osf1 + 58);
    const auto *osf1_59 = buffer.data(osf1 + 59);
    const auto *osf1_60 = buffer.data(osf1 + 60);
    const auto *osf1_63 = buffer.data(osf1 + 63);

    const auto *osg_0 = buffer.data(osg + 0);
    const auto *osg_1 = buffer.data(osg + 1);
    const auto *osg_2 = buffer.data(osg + 2);
    const auto *osg_3 = buffer.data(osg + 3);
    const auto *osg_5 = buffer.data(osg + 5);
    const auto *osg_6 = buffer.data(osg + 6);
    const auto *osg_9 = buffer.data(osg + 9);
    const auto *osg_10 = buffer.data(osg + 10);
    const auto *osg_12 = buffer.data(osg + 12);
    const auto *osg_13 = buffer.data(osg + 13);
    const auto *osg_14 = buffer.data(osg + 14);
    const auto *osg_15 = buffer.data(osg + 15);
    const auto *osg_16 = buffer.data(osg + 16);
    const auto *osg_18 = buffer.data(osg + 18);
    const auto *osg_20 = buffer.data(osg + 20);
    const auto *osg_21 = buffer.data(osg + 21);
    const auto *osg_25 = buffer.data(osg + 25);
    const auto *osg_26 = buffer.data(osg + 26);
    const auto *osg_27 = buffer.data(osg + 27);
    const auto *osg_28 = buffer.data(osg + 28);
    const auto *osg_29 = buffer.data(osg + 29);
    const auto *osg_30 = buffer.data(osg + 30);
    const auto *osg_32 = buffer.data(osg + 32);
    const auto *osg_34 = buffer.data(osg + 34);
    const auto *osg_35 = buffer.data(osg + 35);
    const auto *osg_39 = buffer.data(osg + 39);
    const auto *osg_40 = buffer.data(osg + 40);
    const auto *osg_41 = buffer.data(osg + 41);
    const auto *osg_42 = buffer.data(osg + 42);
    const auto *osg_43 = buffer.data(osg + 43);
    const auto *osg_44 = buffer.data(osg + 44);
    const auto *osg_45 = buffer.data(osg + 45);
    const auto *osg_46 = buffer.data(osg + 46);
    const auto *osg_47 = buffer.data(osg + 47);
    const auto *osg_48 = buffer.data(osg + 48);
    const auto *osg_50 = buffer.data(osg + 50);
    const auto *osg_51 = buffer.data(osg + 51);
    const auto *osg_55 = buffer.data(osg + 55);
    const auto *osg_56 = buffer.data(osg + 56);
    const auto *osg_57 = buffer.data(osg + 57);
    const auto *osg_58 = buffer.data(osg + 58);
    const auto *osg_59 = buffer.data(osg + 59);
    const auto *osg_60 = buffer.data(osg + 60);
    const auto *osg_62 = buffer.data(osg + 62);
    const auto *osg_63 = buffer.data(osg + 63);
    const auto *osg_65 = buffer.data(osg + 65);
    const auto *osg_70 = buffer.data(osg + 70);
    const auto *osg_71 = buffer.data(osg + 71);
    const auto *osg_72 = buffer.data(osg + 72);
    const auto *osg_73 = buffer.data(osg + 73);
    const auto *osg_74 = buffer.data(osg + 74);
    const auto *osg_75 = buffer.data(osg + 75);
    const auto *osg_76 = buffer.data(osg + 76);
    const auto *osg_77 = buffer.data(osg + 77);
    const auto *osg_78 = buffer.data(osg + 78);
    const auto *osg_79 = buffer.data(osg + 79);
    const auto *osg_80 = buffer.data(osg + 80);
    const auto *osg_84 = buffer.data(osg + 84);
    const auto *osg_85 = buffer.data(osg + 85);
    const auto *osg_86 = buffer.data(osg + 86);
    const auto *osg_87 = buffer.data(osg + 87);
    const auto *osg_88 = buffer.data(osg + 88);
    const auto *osg_89 = buffer.data(osg + 89);
    const auto *osg_90 = buffer.data(osg + 90);
    const auto *osg_93 = buffer.data(osg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, nsg_0, osf0_0, \
                         osf1_0, osg_0, osg_1, osg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * nsg_0[k]
                 + f_1 * osf0_0[k]
                 - f_2 * osf1_0[k]
                 + f_3 * pc_x[k] * osg_0[k];

        t_1[k] = f_3 * pc_y[k] * osg_0[k];

        t_2[k] = f_3 * pc_z[k] * osg_0[k];

        t_3[k] = f_4 * osf0_0[k]
                 - f_5 * osf1_0[k]
                 + f_3 * pc_y[k] * osg_1[k];

        t_4[k] = f_3 * pc_y[k] * osg_2[k];

        t_5[k] = f_4 * osf0_0[k]
                 - f_5 * osf1_0[k]
                 + f_3 * pc_z[k] * osg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, nsg_10, osf0_1, osf0_2, \
                         osf1_1, osf1_2, osg_3, osg_5, osg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * osf0_1[k]
                 - f_7 * osf1_1[k]
                 + f_3 * pc_y[k] * osg_3[k];

        t_7[k] = f_3 * pc_z[k] * osg_3[k];

        t_8[k] = f_3 * pc_y[k] * osg_5[k];

        t_9[k] = f_6 * osf0_2[k]
                 - f_7 * osf1_2[k]
                 + f_3 * pc_z[k] * osg_5[k];

        t_10[k] = f_0 * nsg_10[k]
                  + f_3 * pc_x[k] * osg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, nsg_12, nsg_14, osg_6, \
                         osg_9, osg_12, osg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * osg_6[k];

        t_12[k] = f_0 * nsg_12[k]
                  + f_3 * pc_x[k] * osg_12[k];

        t_13[k] = f_3 * pc_y[k] * osg_9[k];

        t_14[k] = f_0 * nsg_14[k]
                  + f_3 * pc_x[k] * osg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, osf0_6, osf0_8, osf0_9, osf1_6, \
                         osf1_8, osf1_9, osg_10, osg_12, osg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * osf0_6[k]
                  - f_2 * osf1_6[k]
                  + f_3 * pc_y[k] * osg_10[k];

        t_16[k] = f_3 * pc_z[k] * osg_10[k];

        t_17[k] = f_6 * osf0_8[k]
                  - f_7 * osf1_8[k]
                  + f_3 * pc_y[k] * osg_12[k];

        t_18[k] = f_4 * osf0_9[k]
                  - f_5 * osf1_9[k]
                  + f_3 * pc_y[k] * osg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, nsh0_0, nsg_0, \
                         nsh1_0, osf0_9, osf1_9, osg_14, osg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * osg_14[k];

        t_20[k] = f_1 * osf0_9[k]
                  - f_2 * osf1_9[k]
                  + f_3 * pc_z[k] * osg_14[k];

        t_21[k] = pa_y[k] * nsh0_0[k]
                  - f_8 * pc_y[k] * nsh1_0[k];

        t_22[k] = f_9 * nsg_0[k]
                  + f_3 * pc_y[k] * osg_15[k];

        t_23[k] = f_3 * pc_z[k] * osg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, nsh0_3, nsh0_5, nsh0_6, \
                         nsg_1, nsg_3, nsh1_3, nsh1_5, nsh1_6, osg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * nsh0_3[k]
                  + f_10 * nsg_1[k]
                  - f_8 * pc_y[k] * nsh1_3[k];

        t_25[k] = f_3 * pc_z[k] * osg_16[k];

        t_26[k] = pa_y[k] * nsh0_5[k]
                  - f_8 * pc_y[k] * nsh1_5[k];

        t_27[k] = pa_y[k] * nsh0_6[k]
                  + f_11 * nsg_3[k]
                  - f_8 * pc_y[k] * nsh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, nsh0_9, nsg_5, \
                         nsg_25, nsh1_9, osg_18, osg_20, osg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * osg_18[k];

        t_29[k] = f_9 * nsg_5[k]
                  + f_3 * pc_y[k] * osg_20[k];

        t_30[k] = pa_y[k] * nsh0_9[k]
                  - f_8 * pc_y[k] * nsh1_9[k];

        t_31[k] = f_12 * nsg_25[k]
                  + f_3 * pc_x[k] * osg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, nsg_27, nsg_28, nsg_29, osg_21, \
                         osg_27, osg_28, osg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * osg_21[k];

        t_33[k] = f_12 * nsg_27[k]
                  + f_3 * pc_x[k] * osg_27[k];

        t_34[k] = f_12 * nsg_28[k]
                  + f_3 * pc_x[k] * osg_28[k];

        t_35[k] = f_12 * nsg_29[k]
                  + f_3 * pc_x[k] * osg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, nsg_10, osf0_16, osf0_17, \
                         osf1_16, osf1_17, osg_25, osg_26, osg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * nsg_10[k]
                  + f_1 * osf0_16[k]
                  - f_2 * osf1_16[k]
                  + f_3 * pc_y[k] * osg_25[k];

        t_37[k] = f_3 * pc_z[k] * osg_25[k];

        t_38[k] = f_4 * osf0_16[k]
                  - f_5 * osf1_16[k]
                  + f_3 * pc_z[k] * osg_26[k];

        t_39[k] = f_6 * osf0_17[k]
                  - f_7 * osf1_17[k]
                  + f_3 * pc_z[k] * osg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, nsh0_0, nsh0_20, \
                         nsg_14, nsh1_0, nsh1_20, osg_29, osg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * nsg_14[k]
                  + f_3 * pc_y[k] * osg_29[k];

        t_41[k] = pa_y[k] * nsh0_20[k]
                  - f_8 * pc_y[k] * nsh1_20[k];

        t_42[k] = pa_z[k] * nsh0_0[k]
                  - f_8 * pc_z[k] * nsh1_0[k];

        t_43[k] = f_3 * pc_y[k] * osg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, nsh0_3, nsh0_5, nsg_0, \
                         nsg_2, nsh1_3, nsh1_5, osg_30, osg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * nsg_0[k]
                  + f_3 * pc_z[k] * osg_30[k];

        t_45[k] = pa_z[k] * nsh0_3[k]
                  - f_8 * pc_z[k] * nsh1_3[k];

        t_46[k] = f_3 * pc_y[k] * osg_32[k];

        t_47[k] = pa_z[k] * nsh0_5[k]
                  + f_10 * nsg_2[k]
                  - f_8 * pc_z[k] * nsh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, nsh0_6, nsh0_9, nsg_5, \
                         nsh1_6, nsh1_9, osf0_22, osf1_22, osg_34, \
                         osg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * nsh0_6[k]
                  - f_8 * pc_z[k] * nsh1_6[k];

        t_49[k] = f_4 * osf0_22[k]
                  - f_5 * osf1_22[k]
                  + f_3 * pc_y[k] * osg_34[k];

        t_50[k] = f_3 * pc_y[k] * osg_35[k];

        t_51[k] = pa_z[k] * nsh0_9[k]
                  + f_11 * nsg_5[k]
                  - f_8 * pc_z[k] * nsh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, nsg_40, nsg_41, nsg_42, \
                         nsg_44, osg_39, osg_40, osg_41, osg_42, \
                         osg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * nsg_40[k]
                  + f_3 * pc_x[k] * osg_40[k];

        t_53[k] = f_12 * nsg_41[k]
                  + f_3 * pc_x[k] * osg_41[k];

        t_54[k] = f_12 * nsg_42[k]
                  + f_3 * pc_x[k] * osg_42[k];

        t_55[k] = f_3 * pc_y[k] * osg_39[k];

        t_56[k] = f_12 * nsg_44[k]
                  + f_3 * pc_x[k] * osg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, nsh0_15, nsh1_15, osf0_27, \
                         osf0_28, osf1_27, osf1_28, osg_41, osg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * nsh0_15[k]
                  - f_8 * pc_z[k] * nsh1_15[k];

        t_58[k] = f_13 * osf0_27[k]
                  - f_14 * osf1_27[k]
                  + f_3 * pc_y[k] * osg_41[k];

        t_59[k] = f_6 * osf0_28[k]
                  - f_7 * osf1_28[k]
                  + f_3 * pc_y[k] * osg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, nsg_14, nsg_45, osf0_29, \
                         osf0_30, osf1_29, osf1_30, osg_43, osg_44, \
                         osg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * osf0_29[k]
                  - f_5 * osf1_29[k]
                  + f_3 * pc_y[k] * osg_43[k];

        t_61[k] = f_3 * pc_y[k] * osg_44[k];

        t_62[k] = f_9 * nsg_14[k]
                  + f_1 * osf0_29[k]
                  - f_2 * osf1_29[k]
                  + f_3 * pc_z[k] * osg_44[k];

        t_63[k] = f_15 * nsg_45[k]
                  + f_1 * osf0_30[k]
                  - f_2 * osf1_30[k]
                  + f_3 * pc_x[k] * osg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, nsg_15, nsg_48, osf0_33, \
                         osf1_33, osg_45, osg_46, osg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * nsg_15[k]
                  + f_3 * pc_y[k] * osg_45[k];

        t_65[k] = f_3 * pc_z[k] * osg_45[k];

        t_66[k] = f_15 * nsg_48[k]
                  + f_6 * osf0_33[k]
                  - f_7 * osf1_33[k]
                  + f_3 * pc_x[k] * osg_48[k];

        t_67[k] = f_3 * pc_z[k] * osg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, nsg_51, osf0_30, osf0_36, osf1_30, \
                         osf1_36, osg_47, osg_48, osg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * osf0_30[k]
                  - f_5 * osf1_30[k]
                  + f_3 * pc_z[k] * osg_47[k];

        t_69[k] = f_15 * nsg_51[k]
                  + f_4 * osf0_36[k]
                  - f_5 * osf1_36[k]
                  + f_3 * pc_x[k] * osg_51[k];

        t_70[k] = f_3 * pc_z[k] * osg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, nsg_20, nsg_55, osf0_32, \
                         osf1_32, osg_50, osg_51, osg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * nsg_20[k]
                  + f_3 * pc_y[k] * osg_50[k];

        t_72[k] = f_6 * osf0_32[k]
                  - f_7 * osf1_32[k]
                  + f_3 * pc_z[k] * osg_50[k];

        t_73[k] = f_15 * nsg_55[k]
                  + f_3 * pc_x[k] * osg_55[k];

        t_74[k] = f_3 * pc_z[k] * osg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, nsg_25, nsg_57, nsg_58, nsg_59, \
                         osf0_36, osf1_36, osg_55, osg_57, osg_58, \
                         osg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * nsg_57[k]
                  + f_3 * pc_x[k] * osg_57[k];

        t_76[k] = f_15 * nsg_58[k]
                  + f_3 * pc_x[k] * osg_58[k];

        t_77[k] = f_15 * nsg_59[k]
                  + f_3 * pc_x[k] * osg_59[k];

        t_78[k] = f_10 * nsg_25[k]
                  + f_1 * osf0_36[k]
                  - f_2 * osf1_36[k]
                  + f_3 * pc_y[k] * osg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, nsg_29, osf0_36, osf0_37, \
                         osf1_36, osf1_37, osg_55, osg_56, osg_57, \
                         osg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * osg_55[k];

        t_80[k] = f_4 * osf0_36[k]
                  - f_5 * osf1_36[k]
                  + f_3 * pc_z[k] * osg_56[k];

        t_81[k] = f_6 * osf0_37[k]
                  - f_7 * osf1_37[k]
                  + f_3 * pc_z[k] * osg_57[k];

        t_82[k] = f_10 * nsg_29[k]
                  + f_3 * pc_y[k] * osg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, nsh0_42, nsg_15, nsg_30, \
                         nsh1_42, osf0_39, osf1_39, osg_59, osg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * osf0_39[k]
                  - f_2 * osf1_39[k]
                  + f_3 * pc_z[k] * osg_59[k];

        t_84[k] = pa_y[k] * nsh0_42[k]
                  - f_8 * pc_y[k] * nsh1_42[k];

        t_85[k] = f_9 * nsg_30[k]
                  + f_3 * pc_y[k] * osg_60[k];

        t_86[k] = f_9 * nsg_15[k]
                  + f_3 * pc_z[k] * osg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, nsh0_24, nsh0_27, \
                         nsh0_47, nsg_32, nsh1_24, nsh1_27, nsh1_47, \
                         osg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * nsh0_24[k]
                  - f_8 * pc_z[k] * nsh1_24[k];

        t_88[k] = f_9 * nsg_32[k]
                  + f_3 * pc_y[k] * osg_62[k];

        t_89[k] = pa_y[k] * nsh0_47[k]
                  - f_8 * pc_y[k] * nsh1_47[k];

        t_90[k] = pa_z[k] * nsh0_27[k]
                  - f_8 * pc_z[k] * nsh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, nsh0_51, nsg_18, \
                         nsg_35, nsg_70, nsh1_51, osg_63, osg_65, \
                         osg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * nsg_18[k]
                  + f_3 * pc_z[k] * osg_63[k];

        t_92[k] = f_9 * nsg_35[k]
                  + f_3 * pc_y[k] * osg_65[k];

        t_93[k] = pa_y[k] * nsh0_51[k]
                  - f_8 * pc_y[k] * nsh1_51[k];

        t_94[k] = f_15 * nsg_70[k]
                  + f_3 * pc_x[k] * osg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, nsg_71, nsg_72, nsg_73, nsg_74, osg_71, \
                         osg_72, osg_73, osg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * nsg_71[k]
                  + f_3 * pc_x[k] * osg_71[k];

        t_96[k] = f_15 * nsg_72[k]
                  + f_3 * pc_x[k] * osg_72[k];

        t_97[k] = f_15 * nsg_73[k]
                  + f_3 * pc_x[k] * osg_73[k];

        t_98[k] = f_15 * nsg_74[k]
                  + f_3 * pc_x[k] * osg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, nsh0_36, nsg_25, nsg_42, \
                         nsh1_36, osf0_48, osf1_48, osg_70, osg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * nsh0_36[k]
                  - f_8 * pc_z[k] * nsh1_36[k];

        t_100[k] = f_9 * nsg_25[k]
                   + f_3 * pc_z[k] * osg_70[k];

        t_101[k] = f_9 * nsg_42[k]
                   + f_6 * osf0_48[k]
                   - f_7 * osf1_48[k]
                   + f_3 * pc_y[k] * osg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, nsh0_62, nsg_43, nsg_44, nsh1_62, \
                         osf0_49, osf1_49, osg_73, osg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * nsg_43[k]
                   + f_4 * osf0_49[k]
                   - f_5 * osf1_49[k]
                   + f_3 * pc_y[k] * osg_73[k];

        t_103[k] = f_9 * nsg_44[k]
                   + f_3 * pc_y[k] * osg_74[k];

        t_104[k] = pa_y[k] * nsh0_62[k]
                   - f_8 * pc_y[k] * nsh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, nsg_30, nsg_75, \
                         osf0_50, osf1_50, osg_75, osg_76, osg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * nsg_75[k]
                   + f_1 * osf0_50[k]
                   - f_2 * osf1_50[k]
                   + f_3 * pc_x[k] * osg_75[k];

        t_106[k] = f_3 * pc_y[k] * osg_75[k];

        t_107[k] = f_10 * nsg_30[k]
                   + f_3 * pc_z[k] * osg_75[k];

        t_108[k] = f_4 * osf0_50[k]
                   - f_5 * osf1_50[k]
                   + f_3 * pc_y[k] * osg_76[k];

        t_109[k] = f_3 * pc_y[k] * osg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, nsg_80, osf0_51, osf0_52, \
                         osf0_55, osf1_51, osf1_52, osf1_55, osg_78, osg_79, \
                         osg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * nsg_80[k]
                   + f_6 * osf0_55[k]
                   - f_7 * osf1_55[k]
                   + f_3 * pc_x[k] * osg_80[k];

        t_111[k] = f_6 * osf0_51[k]
                   - f_7 * osf1_51[k]
                   + f_3 * pc_y[k] * osg_78[k];

        t_112[k] = f_4 * osf0_52[k]
                   - f_5 * osf1_52[k]
                   + f_3 * pc_y[k] * osg_79[k];

        t_113[k] = f_3 * pc_y[k] * osg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, nsg_84, nsg_85, nsg_86, nsg_87, \
                         osf0_59, osf1_59, osg_84, osg_85, osg_86, \
                         osg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * nsg_84[k]
                   + f_4 * osf0_59[k]
                   - f_5 * osf1_59[k]
                   + f_3 * pc_x[k] * osg_84[k];

        t_115[k] = f_15 * nsg_85[k]
                   + f_3 * pc_x[k] * osg_85[k];

        t_116[k] = f_15 * nsg_86[k]
                   + f_3 * pc_x[k] * osg_86[k];

        t_117[k] = f_15 * nsg_87[k]
                   + f_3 * pc_x[k] * osg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, nsg_89, osf0_56, osf0_57, \
                         osf1_56, osf1_57, osg_84, osg_85, osg_86, \
                         osg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * osg_84[k];

        t_119[k] = f_15 * nsg_89[k]
                   + f_3 * pc_x[k] * osg_89[k];

        t_120[k] = f_1 * osf0_56[k]
                   - f_2 * osf1_56[k]
                   + f_3 * pc_y[k] * osg_85[k];

        t_121[k] = f_13 * osf0_57[k]
                   - f_14 * osf1_57[k]
                   + f_3 * pc_y[k] * osg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, nsg_44, osf0_58, osf0_59, \
                         osf1_58, osf1_59, osg_87, osg_88, osg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * osf0_58[k]
                   - f_7 * osf1_58[k]
                   + f_3 * pc_y[k] * osg_87[k];

        t_123[k] = f_4 * osf0_59[k]
                   - f_5 * osf1_59[k]
                   + f_3 * pc_y[k] * osg_88[k];

        t_124[k] = f_3 * pc_y[k] * osg_89[k];

        t_125[k] = f_10 * nsg_44[k]
                   + f_1 * osf0_59[k]
                   - f_2 * osf1_59[k]
                   + f_3 * pc_z[k] * osg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, nsg_45, nsg_90, nsg_93, \
                         osf0_60, osf0_63, osf1_60, osf1_63, osg_90, \
                         osg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_16 * nsg_90[k]
                   + f_1 * osf0_60[k]
                   - f_2 * osf1_60[k]
                   + f_3 * pc_x[k] * osg_90[k];

        t_127[k] = f_11 * nsg_45[k]
                   + f_3 * pc_y[k] * osg_90[k];

        t_128[k] = f_3 * pc_z[k] * osg_90[k];

        t_129[k] = f_16 * nsg_93[k]
                   + f_6 * osf0_63[k]
                   - f_7 * osf1_63[k]
                   + f_3 * pc_x[k] * osg_93[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_63 = buffer.data(nsh0 + 63);
    const auto *nsh0_66 = buffer.data(nsh0 + 66);
    const auto *nsh0_69 = buffer.data(nsh0 + 69);
    const auto *nsh0_78 = buffer.data(nsh0 + 78);
    const auto *nsh0_105 = buffer.data(nsh0 + 105);
    const auto *nsh0_108 = buffer.data(nsh0 + 108);
    const auto *nsh0_110 = buffer.data(nsh0 + 110);
    const auto *nsh0_111 = buffer.data(nsh0 + 111);
    const auto *nsh0_114 = buffer.data(nsh0 + 114);
    const auto *nsh0_125 = buffer.data(nsh0 + 125);
    const auto *nsh0_126 = buffer.data(nsh0 + 126);
    const auto *nsh0_129 = buffer.data(nsh0 + 129);
    const auto *nsh0_132 = buffer.data(nsh0 + 132);
    const auto *nsh0_141 = buffer.data(nsh0 + 141);

    const auto *nsg_45 = buffer.data(nsg + 45);
    const auto *nsg_48 = buffer.data(nsg + 48);
    const auto *nsg_50 = buffer.data(nsg + 50);
    const auto *nsg_55 = buffer.data(nsg + 55);
    const auto *nsg_59 = buffer.data(nsg + 59);
    const auto *nsg_60 = buffer.data(nsg + 60);
    const auto *nsg_62 = buffer.data(nsg + 62);
    const auto *nsg_63 = buffer.data(nsg + 63);
    const auto *nsg_65 = buffer.data(nsg + 65);
    const auto *nsg_70 = buffer.data(nsg + 70);
    const auto *nsg_72 = buffer.data(nsg + 72);
    const auto *nsg_73 = buffer.data(nsg + 73);
    const auto *nsg_74 = buffer.data(nsg + 74);
    const auto *nsg_75 = buffer.data(nsg + 75);
    const auto *nsg_76 = buffer.data(nsg + 76);
    const auto *nsg_77 = buffer.data(nsg + 77);
    const auto *nsg_78 = buffer.data(nsg + 78);
    const auto *nsg_80 = buffer.data(nsg + 80);
    const auto *nsg_85 = buffer.data(nsg + 85);
    const auto *nsg_87 = buffer.data(nsg + 87);
    const auto *nsg_88 = buffer.data(nsg + 88);
    const auto *nsg_89 = buffer.data(nsg + 89);
    const auto *nsg_90 = buffer.data(nsg + 90);
    const auto *nsg_93 = buffer.data(nsg + 93);
    const auto *nsg_95 = buffer.data(nsg + 95);
    const auto *nsg_96 = buffer.data(nsg + 96);
    const auto *nsg_100 = buffer.data(nsg + 100);
    const auto *nsg_102 = buffer.data(nsg + 102);
    const auto *nsg_103 = buffer.data(nsg + 103);
    const auto *nsg_104 = buffer.data(nsg + 104);
    const auto *nsg_105 = buffer.data(nsg + 105);
    const auto *nsg_107 = buffer.data(nsg + 107);
    const auto *nsg_110 = buffer.data(nsg + 110);
    const auto *nsg_114 = buffer.data(nsg + 114);
    const auto *nsg_115 = buffer.data(nsg + 115);
    const auto *nsg_116 = buffer.data(nsg + 116);
    const auto *nsg_117 = buffer.data(nsg + 117);
    const auto *nsg_118 = buffer.data(nsg + 118);
    const auto *nsg_119 = buffer.data(nsg + 119);
    const auto *nsg_130 = buffer.data(nsg + 130);
    const auto *nsg_131 = buffer.data(nsg + 131);
    const auto *nsg_132 = buffer.data(nsg + 132);
    const auto *nsg_133 = buffer.data(nsg + 133);
    const auto *nsg_134 = buffer.data(nsg + 134);
    const auto *nsg_135 = buffer.data(nsg + 135);
    const auto *nsg_140 = buffer.data(nsg + 140);
    const auto *nsg_144 = buffer.data(nsg + 144);
    const auto *nsg_145 = buffer.data(nsg + 145);
    const auto *nsg_146 = buffer.data(nsg + 146);
    const auto *nsg_147 = buffer.data(nsg + 147);
    const auto *nsg_149 = buffer.data(nsg + 149);
    const auto *nsg_150 = buffer.data(nsg + 150);
    const auto *nsg_153 = buffer.data(nsg + 153);
    const auto *nsg_156 = buffer.data(nsg + 156);
    const auto *nsg_160 = buffer.data(nsg + 160);
    const auto *nsg_162 = buffer.data(nsg + 162);
    const auto *nsg_163 = buffer.data(nsg + 163);
    const auto *nsg_164 = buffer.data(nsg + 164);
    const auto *nsg_170 = buffer.data(nsg + 170);
    const auto *nsg_174 = buffer.data(nsg + 174);
    const auto *nsg_175 = buffer.data(nsg + 175);
    const auto *nsg_176 = buffer.data(nsg + 176);
    const auto *nsg_177 = buffer.data(nsg + 177);
    const auto *nsg_178 = buffer.data(nsg + 178);
    const auto *nsg_179 = buffer.data(nsg + 179);

    const auto *nsh1_63 = buffer.data(nsh1 + 63);
    const auto *nsh1_66 = buffer.data(nsh1 + 66);
    const auto *nsh1_69 = buffer.data(nsh1 + 69);
    const auto *nsh1_78 = buffer.data(nsh1 + 78);
    const auto *nsh1_105 = buffer.data(nsh1 + 105);
    const auto *nsh1_108 = buffer.data(nsh1 + 108);
    const auto *nsh1_110 = buffer.data(nsh1 + 110);
    const auto *nsh1_111 = buffer.data(nsh1 + 111);
    const auto *nsh1_114 = buffer.data(nsh1 + 114);
    const auto *nsh1_125 = buffer.data(nsh1 + 125);
    const auto *nsh1_126 = buffer.data(nsh1 + 126);
    const auto *nsh1_129 = buffer.data(nsh1 + 129);
    const auto *nsh1_132 = buffer.data(nsh1 + 132);
    const auto *nsh1_141 = buffer.data(nsh1 + 141);

    const auto *osf0_60 = buffer.data(osf0 + 60);
    const auto *osf0_62 = buffer.data(osf0 + 62);
    const auto *osf0_66 = buffer.data(osf0 + 66);
    const auto *osf0_67 = buffer.data(osf0 + 67);
    const auto *osf0_69 = buffer.data(osf0 + 69);
    const auto *osf0_75 = buffer.data(osf0 + 75);
    const auto *osf0_78 = buffer.data(osf0 + 78);
    const auto *osf0_79 = buffer.data(osf0 + 79);
    const auto *osf0_86 = buffer.data(osf0 + 86);
    const auto *osf0_88 = buffer.data(osf0 + 88);
    const auto *osf0_89 = buffer.data(osf0 + 89);
    const auto *osf0_90 = buffer.data(osf0 + 90);
    const auto *osf0_91 = buffer.data(osf0 + 91);
    const auto *osf0_92 = buffer.data(osf0 + 92);
    const auto *osf0_95 = buffer.data(osf0 + 95);
    const auto *osf0_96 = buffer.data(osf0 + 96);
    const auto *osf0_97 = buffer.data(osf0 + 97);
    const auto *osf0_98 = buffer.data(osf0 + 98);
    const auto *osf0_99 = buffer.data(osf0 + 99);
    const auto *osf0_100 = buffer.data(osf0 + 100);
    const auto *osf0_102 = buffer.data(osf0 + 102);
    const auto *osf0_103 = buffer.data(osf0 + 103);
    const auto *osf0_106 = buffer.data(osf0 + 106);
    const auto *osf0_107 = buffer.data(osf0 + 107);
    const auto *osf0_109 = buffer.data(osf0 + 109);
    const auto *osf0_115 = buffer.data(osf0 + 115);
    const auto *osf0_118 = buffer.data(osf0 + 118);
    const auto *osf0_119 = buffer.data(osf0 + 119);

    const auto *osf1_60 = buffer.data(osf1 + 60);
    const auto *osf1_62 = buffer.data(osf1 + 62);
    const auto *osf1_66 = buffer.data(osf1 + 66);
    const auto *osf1_67 = buffer.data(osf1 + 67);
    const auto *osf1_69 = buffer.data(osf1 + 69);
    const auto *osf1_75 = buffer.data(osf1 + 75);
    const auto *osf1_78 = buffer.data(osf1 + 78);
    const auto *osf1_79 = buffer.data(osf1 + 79);
    const auto *osf1_86 = buffer.data(osf1 + 86);
    const auto *osf1_88 = buffer.data(osf1 + 88);
    const auto *osf1_89 = buffer.data(osf1 + 89);
    const auto *osf1_90 = buffer.data(osf1 + 90);
    const auto *osf1_91 = buffer.data(osf1 + 91);
    const auto *osf1_92 = buffer.data(osf1 + 92);
    const auto *osf1_95 = buffer.data(osf1 + 95);
    const auto *osf1_96 = buffer.data(osf1 + 96);
    const auto *osf1_97 = buffer.data(osf1 + 97);
    const auto *osf1_98 = buffer.data(osf1 + 98);
    const auto *osf1_99 = buffer.data(osf1 + 99);
    const auto *osf1_100 = buffer.data(osf1 + 100);
    const auto *osf1_102 = buffer.data(osf1 + 102);
    const auto *osf1_103 = buffer.data(osf1 + 103);
    const auto *osf1_106 = buffer.data(osf1 + 106);
    const auto *osf1_107 = buffer.data(osf1 + 107);
    const auto *osf1_109 = buffer.data(osf1 + 109);
    const auto *osf1_115 = buffer.data(osf1 + 115);
    const auto *osf1_118 = buffer.data(osf1 + 118);
    const auto *osf1_119 = buffer.data(osf1 + 119);

    const auto *osg_91 = buffer.data(osg + 91);
    const auto *osg_92 = buffer.data(osg + 92);
    const auto *osg_93 = buffer.data(osg + 93);
    const auto *osg_95 = buffer.data(osg + 95);
    const auto *osg_96 = buffer.data(osg + 96);
    const auto *osg_100 = buffer.data(osg + 100);
    const auto *osg_101 = buffer.data(osg + 101);
    const auto *osg_102 = buffer.data(osg + 102);
    const auto *osg_103 = buffer.data(osg + 103);
    const auto *osg_104 = buffer.data(osg + 104);
    const auto *osg_105 = buffer.data(osg + 105);
    const auto *osg_107 = buffer.data(osg + 107);
    const auto *osg_108 = buffer.data(osg + 108);
    const auto *osg_110 = buffer.data(osg + 110);
    const auto *osg_114 = buffer.data(osg + 114);
    const auto *osg_115 = buffer.data(osg + 115);
    const auto *osg_116 = buffer.data(osg + 116);
    const auto *osg_117 = buffer.data(osg + 117);
    const auto *osg_118 = buffer.data(osg + 118);
    const auto *osg_119 = buffer.data(osg + 119);
    const auto *osg_120 = buffer.data(osg + 120);
    const auto *osg_122 = buffer.data(osg + 122);
    const auto *osg_123 = buffer.data(osg + 123);
    const auto *osg_125 = buffer.data(osg + 125);
    const auto *osg_130 = buffer.data(osg + 130);
    const auto *osg_131 = buffer.data(osg + 131);
    const auto *osg_132 = buffer.data(osg + 132);
    const auto *osg_133 = buffer.data(osg + 133);
    const auto *osg_134 = buffer.data(osg + 134);
    const auto *osg_135 = buffer.data(osg + 135);
    const auto *osg_136 = buffer.data(osg + 136);
    const auto *osg_137 = buffer.data(osg + 137);
    const auto *osg_138 = buffer.data(osg + 138);
    const auto *osg_139 = buffer.data(osg + 139);
    const auto *osg_140 = buffer.data(osg + 140);
    const auto *osg_144 = buffer.data(osg + 144);
    const auto *osg_145 = buffer.data(osg + 145);
    const auto *osg_146 = buffer.data(osg + 146);
    const auto *osg_147 = buffer.data(osg + 147);
    const auto *osg_148 = buffer.data(osg + 148);
    const auto *osg_149 = buffer.data(osg + 149);
    const auto *osg_150 = buffer.data(osg + 150);
    const auto *osg_151 = buffer.data(osg + 151);
    const auto *osg_152 = buffer.data(osg + 152);
    const auto *osg_153 = buffer.data(osg + 153);
    const auto *osg_155 = buffer.data(osg + 155);
    const auto *osg_156 = buffer.data(osg + 156);
    const auto *osg_160 = buffer.data(osg + 160);
    const auto *osg_161 = buffer.data(osg + 161);
    const auto *osg_162 = buffer.data(osg + 162);
    const auto *osg_163 = buffer.data(osg + 163);
    const auto *osg_164 = buffer.data(osg + 164);
    const auto *osg_165 = buffer.data(osg + 165);
    const auto *osg_167 = buffer.data(osg + 167);
    const auto *osg_168 = buffer.data(osg + 168);
    const auto *osg_170 = buffer.data(osg + 170);
    const auto *osg_174 = buffer.data(osg + 174);
    const auto *osg_175 = buffer.data(osg + 175);
    const auto *osg_176 = buffer.data(osg + 176);
    const auto *osg_177 = buffer.data(osg + 177);
    const auto *osg_178 = buffer.data(osg + 178);
    const auto *osg_179 = buffer.data(osg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, nsg_96, osf0_60, osf0_66, \
                         osf1_60, osf1_66, osg_91, osg_92, osg_93, \
                         osg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * osg_91[k];

        t_131[k] = f_4 * osf0_60[k]
                   - f_5 * osf1_60[k]
                   + f_3 * pc_z[k] * osg_92[k];

        t_132[k] = f_16 * nsg_96[k]
                   + f_4 * osf0_66[k]
                   - f_5 * osf1_66[k]
                   + f_3 * pc_x[k] * osg_96[k];

        t_133[k] = f_3 * pc_z[k] * osg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, nsg_50, nsg_100, \
                         osf0_62, osf1_62, osg_95, osg_96, osg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * nsg_50[k]
                   + f_3 * pc_y[k] * osg_95[k];

        t_135[k] = f_6 * osf0_62[k]
                   - f_7 * osf1_62[k]
                   + f_3 * pc_z[k] * osg_95[k];

        t_136[k] = f_16 * nsg_100[k]
                   + f_3 * pc_x[k] * osg_100[k];

        t_137[k] = f_3 * pc_z[k] * osg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, nsg_55, nsg_102, nsg_103, \
                         nsg_104, osf0_66, osf1_66, osg_100, osg_102, osg_103, \
                         osg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * nsg_102[k]
                   + f_3 * pc_x[k] * osg_102[k];

        t_139[k] = f_16 * nsg_103[k]
                   + f_3 * pc_x[k] * osg_103[k];

        t_140[k] = f_16 * nsg_104[k]
                   + f_3 * pc_x[k] * osg_104[k];

        t_141[k] = f_11 * nsg_55[k]
                   + f_1 * osf0_66[k]
                   - f_2 * osf1_66[k]
                   + f_3 * pc_y[k] * osg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, nsg_59, osf0_66, osf0_67, \
                         osf1_66, osf1_67, osg_100, osg_101, osg_102, \
                         osg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * osg_100[k];

        t_143[k] = f_4 * osf0_66[k]
                   - f_5 * osf1_66[k]
                   + f_3 * pc_z[k] * osg_101[k];

        t_144[k] = f_6 * osf0_67[k]
                   - f_7 * osf1_67[k]
                   + f_3 * pc_z[k] * osg_102[k];

        t_145[k] = f_11 * nsg_59[k]
                   + f_3 * pc_y[k] * osg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, nsh0_63, nsg_45, \
                         nsg_60, nsh1_63, osf0_69, osf1_69, osg_104, \
                         osg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * osf0_69[k]
                   - f_2 * osf1_69[k]
                   + f_3 * pc_z[k] * osg_104[k];

        t_147[k] = pa_z[k] * nsh0_63[k]
                   - f_8 * pc_z[k] * nsh1_63[k];

        t_148[k] = f_10 * nsg_60[k]
                   + f_3 * pc_y[k] * osg_105[k];

        t_149[k] = f_9 * nsg_45[k]
                   + f_3 * pc_z[k] * osg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, nsh0_66, nsg_62, \
                         nsg_110, nsh1_66, osf0_75, osf1_75, osg_107, \
                         osg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * nsh0_66[k]
                   - f_8 * pc_z[k] * nsh1_66[k];

        t_151[k] = f_10 * nsg_62[k]
                   + f_3 * pc_y[k] * osg_107[k];

        t_152[k] = f_16 * nsg_110[k]
                   + f_6 * osf0_75[k]
                   - f_7 * osf1_75[k]
                   + f_3 * pc_x[k] * osg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, nsh0_69, nsg_48, nsg_65, \
                         nsh1_69, osg_108, osg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * nsh0_69[k]
                   - f_8 * pc_z[k] * nsh1_69[k];

        t_154[k] = f_9 * nsg_48[k]
                   + f_3 * pc_z[k] * osg_108[k];

        t_155[k] = f_10 * nsg_65[k]
                   + f_3 * pc_y[k] * osg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, nsg_114, nsg_115, nsg_116, nsg_117, \
                         osf0_79, osf1_79, osg_114, osg_115, osg_116, \
                         osg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * nsg_114[k]
                   + f_4 * osf0_79[k]
                   - f_5 * osf1_79[k]
                   + f_3 * pc_x[k] * osg_114[k];

        t_157[k] = f_16 * nsg_115[k]
                   + f_3 * pc_x[k] * osg_115[k];

        t_158[k] = f_16 * nsg_116[k]
                   + f_3 * pc_x[k] * osg_116[k];

        t_159[k] = f_16 * nsg_117[k]
                   + f_3 * pc_x[k] * osg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, nsh0_78, nsg_55, \
                         nsg_118, nsg_119, nsh1_78, osg_115, osg_118, \
                         osg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * nsg_118[k]
                   + f_3 * pc_x[k] * osg_118[k];

        t_161[k] = f_16 * nsg_119[k]
                   + f_3 * pc_x[k] * osg_119[k];

        t_162[k] = pa_z[k] * nsh0_78[k]
                   - f_8 * pc_z[k] * nsh1_78[k];

        t_163[k] = f_9 * nsg_55[k]
                   + f_3 * pc_z[k] * osg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, nsg_72, nsg_73, nsg_74, osf0_78, osf0_79, \
                         osf1_78, osf1_79, osg_117, osg_118, osg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * nsg_72[k]
                   + f_6 * osf0_78[k]
                   - f_7 * osf1_78[k]
                   + f_3 * pc_y[k] * osg_117[k];

        t_165[k] = f_10 * nsg_73[k]
                   + f_4 * osf0_79[k]
                   - f_5 * osf1_79[k]
                   + f_3 * pc_y[k] * osg_118[k];

        t_166[k] = f_10 * nsg_74[k]
                   + f_3 * pc_y[k] * osg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, nsh0_105, nsg_59, \
                         nsg_60, nsg_75, nsh1_105, osf0_79, osf1_79, osg_119, \
                         osg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * nsg_59[k]
                   + f_1 * osf0_79[k]
                   - f_2 * osf1_79[k]
                   + f_3 * pc_z[k] * osg_119[k];

        t_168[k] = pa_y[k] * nsh0_105[k]
                   - f_8 * pc_y[k] * nsh1_105[k];

        t_169[k] = f_9 * nsg_75[k]
                   + f_3 * pc_y[k] * osg_120[k];

        t_170[k] = f_10 * nsg_60[k]
                   + f_3 * pc_z[k] * osg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, nsh0_108, nsh0_110, nsh0_111, \
                         nsg_76, nsg_77, nsg_78, nsh1_108, nsh1_110, nsh1_111, \
                         osg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * nsh0_108[k]
                   + f_10 * nsg_76[k]
                   - f_8 * pc_y[k] * nsh1_108[k];

        t_172[k] = f_9 * nsg_77[k]
                   + f_3 * pc_y[k] * osg_122[k];

        t_173[k] = pa_y[k] * nsh0_110[k]
                   - f_8 * pc_y[k] * nsh1_110[k];

        t_174[k] = pa_y[k] * nsh0_111[k]
                   + f_11 * nsg_78[k]
                   - f_8 * pc_y[k] * nsh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, nsh0_114, nsg_63, \
                         nsg_80, nsg_130, nsh1_114, osg_123, osg_125, \
                         osg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * nsg_63[k]
                   + f_3 * pc_z[k] * osg_123[k];

        t_176[k] = f_9 * nsg_80[k]
                   + f_3 * pc_y[k] * osg_125[k];

        t_177[k] = pa_y[k] * nsh0_114[k]
                   - f_8 * pc_y[k] * nsh1_114[k];

        t_178[k] = f_16 * nsg_130[k]
                   + f_3 * pc_x[k] * osg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, nsg_131, nsg_132, nsg_133, nsg_134, \
                         osg_131, osg_132, osg_133, osg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_16 * nsg_131[k]
                   + f_3 * pc_x[k] * osg_131[k];

        t_180[k] = f_16 * nsg_132[k]
                   + f_3 * pc_x[k] * osg_132[k];

        t_181[k] = f_16 * nsg_133[k]
                   + f_3 * pc_x[k] * osg_133[k];

        t_182[k] = f_16 * nsg_134[k]
                   + f_3 * pc_x[k] * osg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, nsg_70, nsg_85, nsg_87, osf0_86, \
                         osf0_88, osf1_86, osf1_88, osg_130, osg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * nsg_85[k]
                   + f_1 * osf0_86[k]
                   - f_2 * osf1_86[k]
                   + f_3 * pc_y[k] * osg_130[k];

        t_184[k] = f_10 * nsg_70[k]
                   + f_3 * pc_z[k] * osg_130[k];

        t_185[k] = f_9 * nsg_87[k]
                   + f_6 * osf0_88[k]
                   - f_7 * osf1_88[k]
                   + f_3 * pc_y[k] * osg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, nsh0_125, nsg_88, nsg_89, nsh1_125, \
                         osf0_89, osf1_89, osg_133, osg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * nsg_88[k]
                   + f_4 * osf0_89[k]
                   - f_5 * osf1_89[k]
                   + f_3 * pc_y[k] * osg_133[k];

        t_187[k] = f_9 * nsg_89[k]
                   + f_3 * pc_y[k] * osg_134[k];

        t_188[k] = pa_y[k] * nsh0_125[k]
                   - f_8 * pc_y[k] * nsh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, nsg_75, nsg_135, \
                         osf0_90, osf1_90, osg_135, osg_136, osg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * nsg_135[k]
                   + f_1 * osf0_90[k]
                   - f_2 * osf1_90[k]
                   + f_3 * pc_x[k] * osg_135[k];

        t_190[k] = f_3 * pc_y[k] * osg_135[k];

        t_191[k] = f_11 * nsg_75[k]
                   + f_3 * pc_z[k] * osg_135[k];

        t_192[k] = f_4 * osf0_90[k]
                   - f_5 * osf1_90[k]
                   + f_3 * pc_y[k] * osg_136[k];

        t_193[k] = f_3 * pc_y[k] * osg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, nsg_140, osf0_91, osf0_92, \
                         osf0_95, osf1_91, osf1_92, osf1_95, osg_138, osg_139, \
                         osg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_16 * nsg_140[k]
                   + f_6 * osf0_95[k]
                   - f_7 * osf1_95[k]
                   + f_3 * pc_x[k] * osg_140[k];

        t_195[k] = f_6 * osf0_91[k]
                   - f_7 * osf1_91[k]
                   + f_3 * pc_y[k] * osg_138[k];

        t_196[k] = f_4 * osf0_92[k]
                   - f_5 * osf1_92[k]
                   + f_3 * pc_y[k] * osg_139[k];

        t_197[k] = f_3 * pc_y[k] * osg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, nsg_144, nsg_145, nsg_146, nsg_147, \
                         osf0_99, osf1_99, osg_144, osg_145, osg_146, \
                         osg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * nsg_144[k]
                   + f_4 * osf0_99[k]
                   - f_5 * osf1_99[k]
                   + f_3 * pc_x[k] * osg_144[k];

        t_199[k] = f_16 * nsg_145[k]
                   + f_3 * pc_x[k] * osg_145[k];

        t_200[k] = f_16 * nsg_146[k]
                   + f_3 * pc_x[k] * osg_146[k];

        t_201[k] = f_16 * nsg_147[k]
                   + f_3 * pc_x[k] * osg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, nsg_149, osf0_96, osf0_97, \
                         osf1_96, osf1_97, osg_144, osg_145, osg_146, \
                         osg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * osg_144[k];

        t_203[k] = f_16 * nsg_149[k]
                   + f_3 * pc_x[k] * osg_149[k];

        t_204[k] = f_1 * osf0_96[k]
                   - f_2 * osf1_96[k]
                   + f_3 * pc_y[k] * osg_145[k];

        t_205[k] = f_13 * osf0_97[k]
                   - f_14 * osf1_97[k]
                   + f_3 * pc_y[k] * osg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, nsg_89, osf0_98, osf0_99, \
                         osf1_98, osf1_99, osg_147, osg_148, osg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * osf0_98[k]
                   - f_7 * osf1_98[k]
                   + f_3 * pc_y[k] * osg_147[k];

        t_207[k] = f_4 * osf0_99[k]
                   - f_5 * osf1_99[k]
                   + f_3 * pc_y[k] * osg_148[k];

        t_208[k] = f_3 * pc_y[k] * osg_149[k];

        t_209[k] = f_11 * nsg_89[k]
                   + f_1 * osf0_99[k]
                   - f_2 * osf1_99[k]
                   + f_3 * pc_z[k] * osg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, nsg_90, nsg_150, \
                         nsg_153, osf0_100, osf0_103, osf1_100, osf1_103, osg_150, \
                         osg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * nsg_150[k]
                   + f_1 * osf0_100[k]
                   - f_2 * osf1_100[k]
                   + f_3 * pc_x[k] * osg_150[k];

        t_211[k] = f_18 * nsg_90[k]
                   + f_3 * pc_y[k] * osg_150[k];

        t_212[k] = f_3 * pc_z[k] * osg_150[k];

        t_213[k] = f_17 * nsg_153[k]
                   + f_6 * osf0_103[k]
                   - f_7 * osf1_103[k]
                   + f_3 * pc_x[k] * osg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, nsg_156, osf0_100, osf0_106, \
                         osf1_100, osf1_106, osg_151, osg_152, osg_153, \
                         osg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * osg_151[k];

        t_215[k] = f_4 * osf0_100[k]
                   - f_5 * osf1_100[k]
                   + f_3 * pc_z[k] * osg_152[k];

        t_216[k] = f_17 * nsg_156[k]
                   + f_4 * osf0_106[k]
                   - f_5 * osf1_106[k]
                   + f_3 * pc_x[k] * osg_156[k];

        t_217[k] = f_3 * pc_z[k] * osg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, nsg_95, nsg_160, \
                         osf0_102, osf1_102, osg_155, osg_156, \
                         osg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_18 * nsg_95[k]
                   + f_3 * pc_y[k] * osg_155[k];

        t_219[k] = f_6 * osf0_102[k]
                   - f_7 * osf1_102[k]
                   + f_3 * pc_z[k] * osg_155[k];

        t_220[k] = f_17 * nsg_160[k]
                   + f_3 * pc_x[k] * osg_160[k];

        t_221[k] = f_3 * pc_z[k] * osg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, nsg_100, nsg_162, nsg_163, \
                         nsg_164, osf0_106, osf1_106, osg_160, osg_162, osg_163, \
                         osg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_17 * nsg_162[k]
                   + f_3 * pc_x[k] * osg_162[k];

        t_223[k] = f_17 * nsg_163[k]
                   + f_3 * pc_x[k] * osg_163[k];

        t_224[k] = f_17 * nsg_164[k]
                   + f_3 * pc_x[k] * osg_164[k];

        t_225[k] = f_18 * nsg_100[k]
                   + f_1 * osf0_106[k]
                   - f_2 * osf1_106[k]
                   + f_3 * pc_y[k] * osg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, nsg_104, osf0_106, osf0_107, \
                         osf1_106, osf1_107, osg_160, osg_161, osg_162, \
                         osg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * osg_160[k];

        t_227[k] = f_4 * osf0_106[k]
                   - f_5 * osf1_106[k]
                   + f_3 * pc_z[k] * osg_161[k];

        t_228[k] = f_6 * osf0_107[k]
                   - f_7 * osf1_107[k]
                   + f_3 * pc_z[k] * osg_162[k];

        t_229[k] = f_18 * nsg_104[k]
                   + f_3 * pc_y[k] * osg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, nsh0_126, nsg_90, \
                         nsg_105, nsh1_126, osf0_109, osf1_109, osg_164, \
                         osg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * osf0_109[k]
                   - f_2 * osf1_109[k]
                   + f_3 * pc_z[k] * osg_164[k];

        t_231[k] = pa_z[k] * nsh0_126[k]
                   - f_8 * pc_z[k] * nsh1_126[k];

        t_232[k] = f_11 * nsg_105[k]
                   + f_3 * pc_y[k] * osg_165[k];

        t_233[k] = f_9 * nsg_90[k]
                   + f_3 * pc_z[k] * osg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, nsh0_129, nsg_107, \
                         nsg_170, nsh1_129, osf0_115, osf1_115, osg_167, \
                         osg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * nsh0_129[k]
                   - f_8 * pc_z[k] * nsh1_129[k];

        t_235[k] = f_11 * nsg_107[k]
                   + f_3 * pc_y[k] * osg_167[k];

        t_236[k] = f_17 * nsg_170[k]
                   + f_6 * osf0_115[k]
                   - f_7 * osf1_115[k]
                   + f_3 * pc_x[k] * osg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, nsh0_132, nsg_93, nsg_110, \
                         nsh1_132, osg_168, osg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * nsh0_132[k]
                   - f_8 * pc_z[k] * nsh1_132[k];

        t_238[k] = f_9 * nsg_93[k]
                   + f_3 * pc_z[k] * osg_168[k];

        t_239[k] = f_11 * nsg_110[k]
                   + f_3 * pc_y[k] * osg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, nsg_174, nsg_175, nsg_176, nsg_177, \
                         osf0_119, osf1_119, osg_174, osg_175, osg_176, \
                         osg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * nsg_174[k]
                   + f_4 * osf0_119[k]
                   - f_5 * osf1_119[k]
                   + f_3 * pc_x[k] * osg_174[k];

        t_241[k] = f_17 * nsg_175[k]
                   + f_3 * pc_x[k] * osg_175[k];

        t_242[k] = f_17 * nsg_176[k]
                   + f_3 * pc_x[k] * osg_176[k];

        t_243[k] = f_17 * nsg_177[k]
                   + f_3 * pc_x[k] * osg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, nsh0_141, nsg_100, \
                         nsg_178, nsg_179, nsh1_141, osg_175, osg_178, \
                         osg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_17 * nsg_178[k]
                   + f_3 * pc_x[k] * osg_178[k];

        t_245[k] = f_17 * nsg_179[k]
                   + f_3 * pc_x[k] * osg_179[k];

        t_246[k] = pa_z[k] * nsh0_141[k]
                   - f_8 * pc_z[k] * nsh1_141[k];

        t_247[k] = f_9 * nsg_100[k]
                   + f_3 * pc_z[k] * osg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, nsg_117, nsg_118, nsg_119, osf0_118, \
                         osf0_119, osf1_118, osf1_119, osg_177, osg_178, \
                         osg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * nsg_117[k]
                   + f_6 * osf0_118[k]
                   - f_7 * osf1_118[k]
                   + f_3 * pc_y[k] * osg_177[k];

        t_249[k] = f_11 * nsg_118[k]
                   + f_4 * osf0_119[k]
                   - f_5 * osf1_119[k]
                   + f_3 * pc_y[k] * osg_178[k];

        t_250[k] = f_11 * nsg_119[k]
                   + f_3 * pc_y[k] * osg_179[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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
    auto *t_368 = buffer.data(target + 368);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_189 = buffer.data(nsh0 + 189);
    const auto *nsh0_192 = buffer.data(nsh0 + 192);
    const auto *nsh0_194 = buffer.data(nsh0 + 194);
    const auto *nsh0_195 = buffer.data(nsh0 + 195);
    const auto *nsh0_198 = buffer.data(nsh0 + 198);
    const auto *nsh0_209 = buffer.data(nsh0 + 209);
    const auto *nsh0_210 = buffer.data(nsh0 + 210);
    const auto *nsh0_213 = buffer.data(nsh0 + 213);
    const auto *nsh0_216 = buffer.data(nsh0 + 216);
    const auto *nsh0_225 = buffer.data(nsh0 + 225);

    const auto *nsg_104 = buffer.data(nsg + 104);
    const auto *nsg_105 = buffer.data(nsg + 105);
    const auto *nsg_108 = buffer.data(nsg + 108);
    const auto *nsg_115 = buffer.data(nsg + 115);
    const auto *nsg_119 = buffer.data(nsg + 119);
    const auto *nsg_120 = buffer.data(nsg + 120);
    const auto *nsg_122 = buffer.data(nsg + 122);
    const auto *nsg_123 = buffer.data(nsg + 123);
    const auto *nsg_125 = buffer.data(nsg + 125);
    const auto *nsg_130 = buffer.data(nsg + 130);
    const auto *nsg_132 = buffer.data(nsg + 132);
    const auto *nsg_133 = buffer.data(nsg + 133);
    const auto *nsg_134 = buffer.data(nsg + 134);
    const auto *nsg_135 = buffer.data(nsg + 135);
    const auto *nsg_136 = buffer.data(nsg + 136);
    const auto *nsg_137 = buffer.data(nsg + 137);
    const auto *nsg_138 = buffer.data(nsg + 138);
    const auto *nsg_140 = buffer.data(nsg + 140);
    const auto *nsg_145 = buffer.data(nsg + 145);
    const auto *nsg_147 = buffer.data(nsg + 147);
    const auto *nsg_148 = buffer.data(nsg + 148);
    const auto *nsg_149 = buffer.data(nsg + 149);
    const auto *nsg_150 = buffer.data(nsg + 150);
    const auto *nsg_153 = buffer.data(nsg + 153);
    const auto *nsg_155 = buffer.data(nsg + 155);
    const auto *nsg_160 = buffer.data(nsg + 160);
    const auto *nsg_164 = buffer.data(nsg + 164);
    const auto *nsg_165 = buffer.data(nsg + 165);
    const auto *nsg_167 = buffer.data(nsg + 167);
    const auto *nsg_168 = buffer.data(nsg + 168);
    const auto *nsg_170 = buffer.data(nsg + 170);
    const auto *nsg_177 = buffer.data(nsg + 177);
    const auto *nsg_178 = buffer.data(nsg + 178);
    const auto *nsg_179 = buffer.data(nsg + 179);
    const auto *nsg_180 = buffer.data(nsg + 180);
    const auto *nsg_182 = buffer.data(nsg + 182);
    const auto *nsg_183 = buffer.data(nsg + 183);
    const auto *nsg_185 = buffer.data(nsg + 185);
    const auto *nsg_186 = buffer.data(nsg + 186);
    const auto *nsg_189 = buffer.data(nsg + 189);
    const auto *nsg_190 = buffer.data(nsg + 190);
    const auto *nsg_191 = buffer.data(nsg + 191);
    const auto *nsg_192 = buffer.data(nsg + 192);
    const auto *nsg_193 = buffer.data(nsg + 193);
    const auto *nsg_194 = buffer.data(nsg + 194);
    const auto *nsg_205 = buffer.data(nsg + 205);
    const auto *nsg_206 = buffer.data(nsg + 206);
    const auto *nsg_207 = buffer.data(nsg + 207);
    const auto *nsg_208 = buffer.data(nsg + 208);
    const auto *nsg_209 = buffer.data(nsg + 209);
    const auto *nsg_210 = buffer.data(nsg + 210);
    const auto *nsg_215 = buffer.data(nsg + 215);
    const auto *nsg_219 = buffer.data(nsg + 219);
    const auto *nsg_220 = buffer.data(nsg + 220);
    const auto *nsg_221 = buffer.data(nsg + 221);
    const auto *nsg_222 = buffer.data(nsg + 222);
    const auto *nsg_224 = buffer.data(nsg + 224);
    const auto *nsg_225 = buffer.data(nsg + 225);
    const auto *nsg_228 = buffer.data(nsg + 228);
    const auto *nsg_231 = buffer.data(nsg + 231);
    const auto *nsg_235 = buffer.data(nsg + 235);
    const auto *nsg_237 = buffer.data(nsg + 237);
    const auto *nsg_238 = buffer.data(nsg + 238);
    const auto *nsg_239 = buffer.data(nsg + 239);
    const auto *nsg_245 = buffer.data(nsg + 245);
    const auto *nsg_249 = buffer.data(nsg + 249);
    const auto *nsg_250 = buffer.data(nsg + 250);
    const auto *nsg_251 = buffer.data(nsg + 251);
    const auto *nsg_252 = buffer.data(nsg + 252);
    const auto *nsg_253 = buffer.data(nsg + 253);
    const auto *nsg_254 = buffer.data(nsg + 254);
    const auto *nsg_255 = buffer.data(nsg + 255);
    const auto *nsg_258 = buffer.data(nsg + 258);
    const auto *nsg_260 = buffer.data(nsg + 260);
    const auto *nsg_261 = buffer.data(nsg + 261);
    const auto *nsg_264 = buffer.data(nsg + 264);
    const auto *nsg_265 = buffer.data(nsg + 265);
    const auto *nsg_266 = buffer.data(nsg + 266);

    const auto *nsh1_189 = buffer.data(nsh1 + 189);
    const auto *nsh1_192 = buffer.data(nsh1 + 192);
    const auto *nsh1_194 = buffer.data(nsh1 + 194);
    const auto *nsh1_195 = buffer.data(nsh1 + 195);
    const auto *nsh1_198 = buffer.data(nsh1 + 198);
    const auto *nsh1_209 = buffer.data(nsh1 + 209);
    const auto *nsh1_210 = buffer.data(nsh1 + 210);
    const auto *nsh1_213 = buffer.data(nsh1 + 213);
    const auto *nsh1_216 = buffer.data(nsh1 + 216);
    const auto *nsh1_225 = buffer.data(nsh1 + 225);

    const auto *osf0_119 = buffer.data(osf0 + 119);
    const auto *osf0_120 = buffer.data(osf0 + 120);
    const auto *osf0_123 = buffer.data(osf0 + 123);
    const auto *osf0_125 = buffer.data(osf0 + 125);
    const auto *osf0_126 = buffer.data(osf0 + 126);
    const auto *osf0_128 = buffer.data(osf0 + 128);
    const auto *osf0_129 = buffer.data(osf0 + 129);
    const auto *osf0_136 = buffer.data(osf0 + 136);
    const auto *osf0_138 = buffer.data(osf0 + 138);
    const auto *osf0_139 = buffer.data(osf0 + 139);
    const auto *osf0_140 = buffer.data(osf0 + 140);
    const auto *osf0_141 = buffer.data(osf0 + 141);
    const auto *osf0_142 = buffer.data(osf0 + 142);
    const auto *osf0_145 = buffer.data(osf0 + 145);
    const auto *osf0_146 = buffer.data(osf0 + 146);
    const auto *osf0_147 = buffer.data(osf0 + 147);
    const auto *osf0_148 = buffer.data(osf0 + 148);
    const auto *osf0_149 = buffer.data(osf0 + 149);
    const auto *osf0_150 = buffer.data(osf0 + 150);
    const auto *osf0_152 = buffer.data(osf0 + 152);
    const auto *osf0_153 = buffer.data(osf0 + 153);
    const auto *osf0_156 = buffer.data(osf0 + 156);
    const auto *osf0_157 = buffer.data(osf0 + 157);
    const auto *osf0_159 = buffer.data(osf0 + 159);
    const auto *osf0_165 = buffer.data(osf0 + 165);
    const auto *osf0_168 = buffer.data(osf0 + 168);
    const auto *osf0_169 = buffer.data(osf0 + 169);
    const auto *osf0_170 = buffer.data(osf0 + 170);
    const auto *osf0_173 = buffer.data(osf0 + 173);
    const auto *osf0_175 = buffer.data(osf0 + 175);
    const auto *osf0_176 = buffer.data(osf0 + 176);
    const auto *osf0_179 = buffer.data(osf0 + 179);

    const auto *osf1_119 = buffer.data(osf1 + 119);
    const auto *osf1_120 = buffer.data(osf1 + 120);
    const auto *osf1_123 = buffer.data(osf1 + 123);
    const auto *osf1_125 = buffer.data(osf1 + 125);
    const auto *osf1_126 = buffer.data(osf1 + 126);
    const auto *osf1_128 = buffer.data(osf1 + 128);
    const auto *osf1_129 = buffer.data(osf1 + 129);
    const auto *osf1_136 = buffer.data(osf1 + 136);
    const auto *osf1_138 = buffer.data(osf1 + 138);
    const auto *osf1_139 = buffer.data(osf1 + 139);
    const auto *osf1_140 = buffer.data(osf1 + 140);
    const auto *osf1_141 = buffer.data(osf1 + 141);
    const auto *osf1_142 = buffer.data(osf1 + 142);
    const auto *osf1_145 = buffer.data(osf1 + 145);
    const auto *osf1_146 = buffer.data(osf1 + 146);
    const auto *osf1_147 = buffer.data(osf1 + 147);
    const auto *osf1_148 = buffer.data(osf1 + 148);
    const auto *osf1_149 = buffer.data(osf1 + 149);
    const auto *osf1_150 = buffer.data(osf1 + 150);
    const auto *osf1_152 = buffer.data(osf1 + 152);
    const auto *osf1_153 = buffer.data(osf1 + 153);
    const auto *osf1_156 = buffer.data(osf1 + 156);
    const auto *osf1_157 = buffer.data(osf1 + 157);
    const auto *osf1_159 = buffer.data(osf1 + 159);
    const auto *osf1_165 = buffer.data(osf1 + 165);
    const auto *osf1_168 = buffer.data(osf1 + 168);
    const auto *osf1_169 = buffer.data(osf1 + 169);
    const auto *osf1_170 = buffer.data(osf1 + 170);
    const auto *osf1_173 = buffer.data(osf1 + 173);
    const auto *osf1_175 = buffer.data(osf1 + 175);
    const auto *osf1_176 = buffer.data(osf1 + 176);
    const auto *osf1_179 = buffer.data(osf1 + 179);

    const auto *osg_179 = buffer.data(osg + 179);
    const auto *osg_180 = buffer.data(osg + 180);
    const auto *osg_182 = buffer.data(osg + 182);
    const auto *osg_183 = buffer.data(osg + 183);
    const auto *osg_185 = buffer.data(osg + 185);
    const auto *osg_186 = buffer.data(osg + 186);
    const auto *osg_189 = buffer.data(osg + 189);
    const auto *osg_190 = buffer.data(osg + 190);
    const auto *osg_191 = buffer.data(osg + 191);
    const auto *osg_192 = buffer.data(osg + 192);
    const auto *osg_193 = buffer.data(osg + 193);
    const auto *osg_194 = buffer.data(osg + 194);
    const auto *osg_195 = buffer.data(osg + 195);
    const auto *osg_197 = buffer.data(osg + 197);
    const auto *osg_198 = buffer.data(osg + 198);
    const auto *osg_200 = buffer.data(osg + 200);
    const auto *osg_205 = buffer.data(osg + 205);
    const auto *osg_206 = buffer.data(osg + 206);
    const auto *osg_207 = buffer.data(osg + 207);
    const auto *osg_208 = buffer.data(osg + 208);
    const auto *osg_209 = buffer.data(osg + 209);
    const auto *osg_210 = buffer.data(osg + 210);
    const auto *osg_211 = buffer.data(osg + 211);
    const auto *osg_212 = buffer.data(osg + 212);
    const auto *osg_213 = buffer.data(osg + 213);
    const auto *osg_214 = buffer.data(osg + 214);
    const auto *osg_215 = buffer.data(osg + 215);
    const auto *osg_219 = buffer.data(osg + 219);
    const auto *osg_220 = buffer.data(osg + 220);
    const auto *osg_221 = buffer.data(osg + 221);
    const auto *osg_222 = buffer.data(osg + 222);
    const auto *osg_223 = buffer.data(osg + 223);
    const auto *osg_224 = buffer.data(osg + 224);
    const auto *osg_225 = buffer.data(osg + 225);
    const auto *osg_226 = buffer.data(osg + 226);
    const auto *osg_227 = buffer.data(osg + 227);
    const auto *osg_228 = buffer.data(osg + 228);
    const auto *osg_230 = buffer.data(osg + 230);
    const auto *osg_231 = buffer.data(osg + 231);
    const auto *osg_235 = buffer.data(osg + 235);
    const auto *osg_236 = buffer.data(osg + 236);
    const auto *osg_237 = buffer.data(osg + 237);
    const auto *osg_238 = buffer.data(osg + 238);
    const auto *osg_239 = buffer.data(osg + 239);
    const auto *osg_240 = buffer.data(osg + 240);
    const auto *osg_242 = buffer.data(osg + 242);
    const auto *osg_243 = buffer.data(osg + 243);
    const auto *osg_245 = buffer.data(osg + 245);
    const auto *osg_249 = buffer.data(osg + 249);
    const auto *osg_250 = buffer.data(osg + 250);
    const auto *osg_251 = buffer.data(osg + 251);
    const auto *osg_252 = buffer.data(osg + 252);
    const auto *osg_253 = buffer.data(osg + 253);
    const auto *osg_254 = buffer.data(osg + 254);
    const auto *osg_255 = buffer.data(osg + 255);
    const auto *osg_257 = buffer.data(osg + 257);
    const auto *osg_258 = buffer.data(osg + 258);
    const auto *osg_260 = buffer.data(osg + 260);
    const auto *osg_261 = buffer.data(osg + 261);
    const auto *osg_264 = buffer.data(osg + 264);
    const auto *osg_265 = buffer.data(osg + 265);
    const auto *osg_266 = buffer.data(osg + 266);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, nsg_104, nsg_120, nsg_180, \
                         osf0_119, osf0_120, osf1_119, osf1_120, osg_179, \
                         osg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * nsg_104[k]
                   + f_1 * osf0_119[k]
                   - f_2 * osf1_119[k]
                   + f_3 * pc_z[k] * osg_179[k];

        t_252[k] = f_17 * nsg_180[k]
                   + f_1 * osf0_120[k]
                   - f_2 * osf1_120[k]
                   + f_3 * pc_x[k] * osg_180[k];

        t_253[k] = f_10 * nsg_120[k]
                   + f_3 * pc_y[k] * osg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, nsg_105, nsg_122, nsg_183, \
                         osf0_123, osf1_123, osg_180, osg_182, \
                         osg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * nsg_105[k]
                   + f_3 * pc_z[k] * osg_180[k];

        t_255[k] = f_17 * nsg_183[k]
                   + f_6 * osf0_123[k]
                   - f_7 * osf1_123[k]
                   + f_3 * pc_x[k] * osg_183[k];

        t_256[k] = f_10 * nsg_122[k]
                   + f_3 * pc_y[k] * osg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, nsg_108, nsg_185, nsg_186, osf0_125, \
                         osf0_126, osf1_125, osf1_126, osg_183, osg_185, \
                         osg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * nsg_185[k]
                   + f_6 * osf0_125[k]
                   - f_7 * osf1_125[k]
                   + f_3 * pc_x[k] * osg_185[k];

        t_258[k] = f_17 * nsg_186[k]
                   + f_4 * osf0_126[k]
                   - f_5 * osf1_126[k]
                   + f_3 * pc_x[k] * osg_186[k];

        t_259[k] = f_10 * nsg_108[k]
                   + f_3 * pc_z[k] * osg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, nsg_125, nsg_189, nsg_190, \
                         nsg_191, osf0_129, osf1_129, osg_185, osg_189, osg_190, \
                         osg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * nsg_125[k]
                   + f_3 * pc_y[k] * osg_185[k];

        t_261[k] = f_17 * nsg_189[k]
                   + f_4 * osf0_129[k]
                   - f_5 * osf1_129[k]
                   + f_3 * pc_x[k] * osg_189[k];

        t_262[k] = f_17 * nsg_190[k]
                   + f_3 * pc_x[k] * osg_190[k];

        t_263[k] = f_17 * nsg_191[k]
                   + f_3 * pc_x[k] * osg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, nsg_130, nsg_192, nsg_193, \
                         nsg_194, osf0_126, osf1_126, osg_190, osg_192, osg_193, \
                         osg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * nsg_192[k]
                   + f_3 * pc_x[k] * osg_192[k];

        t_265[k] = f_17 * nsg_193[k]
                   + f_3 * pc_x[k] * osg_193[k];

        t_266[k] = f_17 * nsg_194[k]
                   + f_3 * pc_x[k] * osg_194[k];

        t_267[k] = f_10 * nsg_130[k]
                   + f_1 * osf0_126[k]
                   - f_2 * osf1_126[k]
                   + f_3 * pc_y[k] * osg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, nsg_115, nsg_132, nsg_133, osf0_128, \
                         osf0_129, osf1_128, osf1_129, osg_190, osg_192, \
                         osg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * nsg_115[k]
                   + f_3 * pc_z[k] * osg_190[k];

        t_269[k] = f_10 * nsg_132[k]
                   + f_6 * osf0_128[k]
                   - f_7 * osf1_128[k]
                   + f_3 * pc_y[k] * osg_192[k];

        t_270[k] = f_10 * nsg_133[k]
                   + f_4 * osf0_129[k]
                   - f_5 * osf1_129[k]
                   + f_3 * pc_y[k] * osg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, nsh0_189, nsg_119, \
                         nsg_134, nsg_135, nsh1_189, osf0_129, osf1_129, osg_194, \
                         osg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * nsg_134[k]
                   + f_3 * pc_y[k] * osg_194[k];

        t_272[k] = f_10 * nsg_119[k]
                   + f_1 * osf0_129[k]
                   - f_2 * osf1_129[k]
                   + f_3 * pc_z[k] * osg_194[k];

        t_273[k] = pa_y[k] * nsh0_189[k]
                   - f_8 * pc_y[k] * nsh1_189[k];

        t_274[k] = f_9 * nsg_135[k]
                   + f_3 * pc_y[k] * osg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, nsh0_192, nsh0_194, \
                         nsg_120, nsg_136, nsg_137, nsh1_192, nsh1_194, osg_195, \
                         osg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * nsg_120[k]
                   + f_3 * pc_z[k] * osg_195[k];

        t_276[k] = pa_y[k] * nsh0_192[k]
                   + f_10 * nsg_136[k]
                   - f_8 * pc_y[k] * nsh1_192[k];

        t_277[k] = f_9 * nsg_137[k]
                   + f_3 * pc_y[k] * osg_197[k];

        t_278[k] = pa_y[k] * nsh0_194[k]
                   - f_8 * pc_y[k] * nsh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, nsh0_195, nsh0_198, \
                         nsg_123, nsg_138, nsg_140, nsh1_195, nsh1_198, osg_198, \
                         osg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * nsh0_195[k]
                   + f_11 * nsg_138[k]
                   - f_8 * pc_y[k] * nsh1_195[k];

        t_280[k] = f_11 * nsg_123[k]
                   + f_3 * pc_z[k] * osg_198[k];

        t_281[k] = f_9 * nsg_140[k]
                   + f_3 * pc_y[k] * osg_200[k];

        t_282[k] = pa_y[k] * nsh0_198[k]
                   - f_8 * pc_y[k] * nsh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, nsg_205, nsg_206, nsg_207, \
                         nsg_208, nsg_209, osg_205, osg_206, osg_207, osg_208, \
                         osg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_17 * nsg_205[k]
                   + f_3 * pc_x[k] * osg_205[k];

        t_284[k] = f_17 * nsg_206[k]
                   + f_3 * pc_x[k] * osg_206[k];

        t_285[k] = f_17 * nsg_207[k]
                   + f_3 * pc_x[k] * osg_207[k];

        t_286[k] = f_17 * nsg_208[k]
                   + f_3 * pc_x[k] * osg_208[k];

        t_287[k] = f_17 * nsg_209[k]
                   + f_3 * pc_x[k] * osg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, nsg_130, nsg_145, nsg_147, osf0_136, \
                         osf0_138, osf1_136, osf1_138, osg_205, \
                         osg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * nsg_145[k]
                   + f_1 * osf0_136[k]
                   - f_2 * osf1_136[k]
                   + f_3 * pc_y[k] * osg_205[k];

        t_289[k] = f_11 * nsg_130[k]
                   + f_3 * pc_z[k] * osg_205[k];

        t_290[k] = f_9 * nsg_147[k]
                   + f_6 * osf0_138[k]
                   - f_7 * osf1_138[k]
                   + f_3 * pc_y[k] * osg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, nsh0_209, nsg_148, nsg_149, \
                         nsh1_209, osf0_139, osf1_139, osg_208, \
                         osg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * nsg_148[k]
                   + f_4 * osf0_139[k]
                   - f_5 * osf1_139[k]
                   + f_3 * pc_y[k] * osg_208[k];

        t_292[k] = f_9 * nsg_149[k]
                   + f_3 * pc_y[k] * osg_209[k];

        t_293[k] = pa_y[k] * nsh0_209[k]
                   - f_8 * pc_y[k] * nsh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, nsg_135, \
                         nsg_210, osf0_140, osf1_140, osg_210, osg_211, \
                         osg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_17 * nsg_210[k]
                   + f_1 * osf0_140[k]
                   - f_2 * osf1_140[k]
                   + f_3 * pc_x[k] * osg_210[k];

        t_295[k] = f_3 * pc_y[k] * osg_210[k];

        t_296[k] = f_18 * nsg_135[k]
                   + f_3 * pc_z[k] * osg_210[k];

        t_297[k] = f_4 * osf0_140[k]
                   - f_5 * osf1_140[k]
                   + f_3 * pc_y[k] * osg_211[k];

        t_298[k] = f_3 * pc_y[k] * osg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, nsg_215, osf0_141, osf0_142, \
                         osf0_145, osf1_141, osf1_142, osf1_145, osg_213, osg_214, \
                         osg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_17 * nsg_215[k]
                   + f_6 * osf0_145[k]
                   - f_7 * osf1_145[k]
                   + f_3 * pc_x[k] * osg_215[k];

        t_300[k] = f_6 * osf0_141[k]
                   - f_7 * osf1_141[k]
                   + f_3 * pc_y[k] * osg_213[k];

        t_301[k] = f_4 * osf0_142[k]
                   - f_5 * osf1_142[k]
                   + f_3 * pc_y[k] * osg_214[k];

        t_302[k] = f_3 * pc_y[k] * osg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, nsg_219, nsg_220, nsg_221, nsg_222, \
                         osf0_149, osf1_149, osg_219, osg_220, osg_221, \
                         osg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_17 * nsg_219[k]
                   + f_4 * osf0_149[k]
                   - f_5 * osf1_149[k]
                   + f_3 * pc_x[k] * osg_219[k];

        t_304[k] = f_17 * nsg_220[k]
                   + f_3 * pc_x[k] * osg_220[k];

        t_305[k] = f_17 * nsg_221[k]
                   + f_3 * pc_x[k] * osg_221[k];

        t_306[k] = f_17 * nsg_222[k]
                   + f_3 * pc_x[k] * osg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, nsg_224, osf0_146, osf0_147, \
                         osf1_146, osf1_147, osg_219, osg_220, osg_221, \
                         osg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * osg_219[k];

        t_308[k] = f_17 * nsg_224[k]
                   + f_3 * pc_x[k] * osg_224[k];

        t_309[k] = f_1 * osf0_146[k]
                   - f_2 * osf1_146[k]
                   + f_3 * pc_y[k] * osg_220[k];

        t_310[k] = f_13 * osf0_147[k]
                   - f_14 * osf1_147[k]
                   + f_3 * pc_y[k] * osg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, nsg_149, osf0_148, osf0_149, \
                         osf1_148, osf1_149, osg_222, osg_223, \
                         osg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * osf0_148[k]
                   - f_7 * osf1_148[k]
                   + f_3 * pc_y[k] * osg_222[k];

        t_312[k] = f_4 * osf0_149[k]
                   - f_5 * osf1_149[k]
                   + f_3 * pc_y[k] * osg_223[k];

        t_313[k] = f_3 * pc_y[k] * osg_224[k];

        t_314[k] = f_18 * nsg_149[k]
                   + f_1 * osf0_149[k]
                   - f_2 * osf1_149[k]
                   + f_3 * pc_z[k] * osg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, nsg_150, nsg_225, \
                         nsg_228, osf0_150, osf0_153, osf1_150, osf1_153, osg_225, \
                         osg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_19 * nsg_225[k]
                   + f_1 * osf0_150[k]
                   - f_2 * osf1_150[k]
                   + f_3 * pc_x[k] * osg_225[k];

        t_316[k] = f_20 * nsg_150[k]
                   + f_3 * pc_y[k] * osg_225[k];

        t_317[k] = f_3 * pc_z[k] * osg_225[k];

        t_318[k] = f_19 * nsg_228[k]
                   + f_6 * osf0_153[k]
                   - f_7 * osf1_153[k]
                   + f_3 * pc_x[k] * osg_228[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_z, nsg_231, osf0_150, osf0_156, \
                         osf1_150, osf1_156, osg_226, osg_227, osg_228, \
                         osg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * osg_226[k];

        t_320[k] = f_4 * osf0_150[k]
                   - f_5 * osf1_150[k]
                   + f_3 * pc_z[k] * osg_227[k];

        t_321[k] = f_19 * nsg_231[k]
                   + f_4 * osf0_156[k]
                   - f_5 * osf1_156[k]
                   + f_3 * pc_x[k] * osg_231[k];

        t_322[k] = f_3 * pc_z[k] * osg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, nsg_155, nsg_235, \
                         osf0_152, osf1_152, osg_230, osg_231, \
                         osg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_20 * nsg_155[k]
                   + f_3 * pc_y[k] * osg_230[k];

        t_324[k] = f_6 * osf0_152[k]
                   - f_7 * osf1_152[k]
                   + f_3 * pc_z[k] * osg_230[k];

        t_325[k] = f_19 * nsg_235[k]
                   + f_3 * pc_x[k] * osg_235[k];

        t_326[k] = f_3 * pc_z[k] * osg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, nsg_160, nsg_237, nsg_238, \
                         nsg_239, osf0_156, osf1_156, osg_235, osg_237, osg_238, \
                         osg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_19 * nsg_237[k]
                   + f_3 * pc_x[k] * osg_237[k];

        t_328[k] = f_19 * nsg_238[k]
                   + f_3 * pc_x[k] * osg_238[k];

        t_329[k] = f_19 * nsg_239[k]
                   + f_3 * pc_x[k] * osg_239[k];

        t_330[k] = f_20 * nsg_160[k]
                   + f_1 * osf0_156[k]
                   - f_2 * osf1_156[k]
                   + f_3 * pc_y[k] * osg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, nsg_164, osf0_156, osf0_157, \
                         osf1_156, osf1_157, osg_235, osg_236, osg_237, \
                         osg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * osg_235[k];

        t_332[k] = f_4 * osf0_156[k]
                   - f_5 * osf1_156[k]
                   + f_3 * pc_z[k] * osg_236[k];

        t_333[k] = f_6 * osf0_157[k]
                   - f_7 * osf1_157[k]
                   + f_3 * pc_z[k] * osg_237[k];

        t_334[k] = f_20 * nsg_164[k]
                   + f_3 * pc_y[k] * osg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pc_y, pc_z, nsh0_210, nsg_150, \
                         nsg_165, nsh1_210, osf0_159, osf1_159, osg_239, \
                         osg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * osf0_159[k]
                   - f_2 * osf1_159[k]
                   + f_3 * pc_z[k] * osg_239[k];

        t_336[k] = pa_z[k] * nsh0_210[k]
                   - f_8 * pc_z[k] * nsh1_210[k];

        t_337[k] = f_18 * nsg_165[k]
                   + f_3 * pc_y[k] * osg_240[k];

        t_338[k] = f_9 * nsg_150[k]
                   + f_3 * pc_z[k] * osg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_x, pc_y, pc_z, nsh0_213, nsg_167, \
                         nsg_245, nsh1_213, osf0_165, osf1_165, osg_242, \
                         osg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * nsh0_213[k]
                   - f_8 * pc_z[k] * nsh1_213[k];

        t_340[k] = f_18 * nsg_167[k]
                   + f_3 * pc_y[k] * osg_242[k];

        t_341[k] = f_19 * nsg_245[k]
                   + f_6 * osf0_165[k]
                   - f_7 * osf1_165[k]
                   + f_3 * pc_x[k] * osg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, nsh0_216, nsg_153, nsg_170, \
                         nsh1_216, osg_243, osg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * nsh0_216[k]
                   - f_8 * pc_z[k] * nsh1_216[k];

        t_343[k] = f_9 * nsg_153[k]
                   + f_3 * pc_z[k] * osg_243[k];

        t_344[k] = f_18 * nsg_170[k]
                   + f_3 * pc_y[k] * osg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, nsg_249, nsg_250, nsg_251, nsg_252, \
                         osf0_169, osf1_169, osg_249, osg_250, osg_251, \
                         osg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_19 * nsg_249[k]
                   + f_4 * osf0_169[k]
                   - f_5 * osf1_169[k]
                   + f_3 * pc_x[k] * osg_249[k];

        t_346[k] = f_19 * nsg_250[k]
                   + f_3 * pc_x[k] * osg_250[k];

        t_347[k] = f_19 * nsg_251[k]
                   + f_3 * pc_x[k] * osg_251[k];

        t_348[k] = f_19 * nsg_252[k]
                   + f_3 * pc_x[k] * osg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, nsh0_225, nsg_160, \
                         nsg_253, nsg_254, nsh1_225, osg_250, osg_253, \
                         osg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_19 * nsg_253[k]
                   + f_3 * pc_x[k] * osg_253[k];

        t_350[k] = f_19 * nsg_254[k]
                   + f_3 * pc_x[k] * osg_254[k];

        t_351[k] = pa_z[k] * nsh0_225[k]
                   - f_8 * pc_z[k] * nsh1_225[k];

        t_352[k] = f_9 * nsg_160[k]
                   + f_3 * pc_z[k] * osg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, nsg_177, nsg_178, nsg_179, osf0_168, \
                         osf0_169, osf1_168, osf1_169, osg_252, osg_253, \
                         osg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_18 * nsg_177[k]
                   + f_6 * osf0_168[k]
                   - f_7 * osf1_168[k]
                   + f_3 * pc_y[k] * osg_252[k];

        t_354[k] = f_18 * nsg_178[k]
                   + f_4 * osf0_169[k]
                   - f_5 * osf1_169[k]
                   + f_3 * pc_y[k] * osg_253[k];

        t_355[k] = f_18 * nsg_179[k]
                   + f_3 * pc_y[k] * osg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, nsg_164, nsg_180, nsg_255, \
                         osf0_169, osf0_170, osf1_169, osf1_170, osg_254, \
                         osg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * nsg_164[k]
                   + f_1 * osf0_169[k]
                   - f_2 * osf1_169[k]
                   + f_3 * pc_z[k] * osg_254[k];

        t_357[k] = f_19 * nsg_255[k]
                   + f_1 * osf0_170[k]
                   - f_2 * osf1_170[k]
                   + f_3 * pc_x[k] * osg_255[k];

        t_358[k] = f_11 * nsg_180[k]
                   + f_3 * pc_y[k] * osg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, nsg_165, nsg_182, nsg_258, \
                         osf0_173, osf1_173, osg_255, osg_257, \
                         osg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * nsg_165[k]
                   + f_3 * pc_z[k] * osg_255[k];

        t_360[k] = f_19 * nsg_258[k]
                   + f_6 * osf0_173[k]
                   - f_7 * osf1_173[k]
                   + f_3 * pc_x[k] * osg_258[k];

        t_361[k] = f_11 * nsg_182[k]
                   + f_3 * pc_y[k] * osg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, nsg_168, nsg_260, nsg_261, osf0_175, \
                         osf0_176, osf1_175, osf1_176, osg_258, osg_260, \
                         osg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_19 * nsg_260[k]
                   + f_6 * osf0_175[k]
                   - f_7 * osf1_175[k]
                   + f_3 * pc_x[k] * osg_260[k];

        t_363[k] = f_19 * nsg_261[k]
                   + f_4 * osf0_176[k]
                   - f_5 * osf1_176[k]
                   + f_3 * pc_x[k] * osg_261[k];

        t_364[k] = f_10 * nsg_168[k]
                   + f_3 * pc_z[k] * osg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, nsg_185, nsg_264, nsg_265, \
                         nsg_266, osf0_179, osf1_179, osg_260, osg_264, osg_265, \
                         osg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * nsg_185[k]
                   + f_3 * pc_y[k] * osg_260[k];

        t_366[k] = f_19 * nsg_264[k]
                   + f_4 * osf0_179[k]
                   - f_5 * osf1_179[k]
                   + f_3 * pc_x[k] * osg_264[k];

        t_367[k] = f_19 * nsg_265[k]
                   + f_3 * pc_x[k] * osg_265[k];

        t_368[k] = f_19 * nsg_266[k]
                   + f_3 * pc_x[k] * osg_266[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_294 = buffer.data(nsh0 + 294);
    const auto *nsh0_297 = buffer.data(nsh0 + 297);
    const auto *nsh0_299 = buffer.data(nsh0 + 299);
    const auto *nsh0_300 = buffer.data(nsh0 + 300);
    const auto *nsh0_303 = buffer.data(nsh0 + 303);
    const auto *nsh0_314 = buffer.data(nsh0 + 314);
    const auto *nsh0_315 = buffer.data(nsh0 + 315);
    const auto *nsh0_318 = buffer.data(nsh0 + 318);
    const auto *nsh0_321 = buffer.data(nsh0 + 321);
    const auto *nsh0_330 = buffer.data(nsh0 + 330);

    const auto *nsg_175 = buffer.data(nsg + 175);
    const auto *nsg_179 = buffer.data(nsg + 179);
    const auto *nsg_180 = buffer.data(nsg + 180);
    const auto *nsg_183 = buffer.data(nsg + 183);
    const auto *nsg_190 = buffer.data(nsg + 190);
    const auto *nsg_192 = buffer.data(nsg + 192);
    const auto *nsg_193 = buffer.data(nsg + 193);
    const auto *nsg_194 = buffer.data(nsg + 194);
    const auto *nsg_195 = buffer.data(nsg + 195);
    const auto *nsg_197 = buffer.data(nsg + 197);
    const auto *nsg_198 = buffer.data(nsg + 198);
    const auto *nsg_200 = buffer.data(nsg + 200);
    const auto *nsg_205 = buffer.data(nsg + 205);
    const auto *nsg_207 = buffer.data(nsg + 207);
    const auto *nsg_208 = buffer.data(nsg + 208);
    const auto *nsg_209 = buffer.data(nsg + 209);
    const auto *nsg_210 = buffer.data(nsg + 210);
    const auto *nsg_211 = buffer.data(nsg + 211);
    const auto *nsg_212 = buffer.data(nsg + 212);
    const auto *nsg_213 = buffer.data(nsg + 213);
    const auto *nsg_215 = buffer.data(nsg + 215);
    const auto *nsg_220 = buffer.data(nsg + 220);
    const auto *nsg_222 = buffer.data(nsg + 222);
    const auto *nsg_223 = buffer.data(nsg + 223);
    const auto *nsg_224 = buffer.data(nsg + 224);
    const auto *nsg_225 = buffer.data(nsg + 225);
    const auto *nsg_228 = buffer.data(nsg + 228);
    const auto *nsg_230 = buffer.data(nsg + 230);
    const auto *nsg_235 = buffer.data(nsg + 235);
    const auto *nsg_239 = buffer.data(nsg + 239);
    const auto *nsg_240 = buffer.data(nsg + 240);
    const auto *nsg_242 = buffer.data(nsg + 242);
    const auto *nsg_245 = buffer.data(nsg + 245);
    const auto *nsg_252 = buffer.data(nsg + 252);
    const auto *nsg_253 = buffer.data(nsg + 253);
    const auto *nsg_254 = buffer.data(nsg + 254);
    const auto *nsg_255 = buffer.data(nsg + 255);
    const auto *nsg_257 = buffer.data(nsg + 257);
    const auto *nsg_267 = buffer.data(nsg + 267);
    const auto *nsg_268 = buffer.data(nsg + 268);
    const auto *nsg_269 = buffer.data(nsg + 269);
    const auto *nsg_270 = buffer.data(nsg + 270);
    const auto *nsg_273 = buffer.data(nsg + 273);
    const auto *nsg_275 = buffer.data(nsg + 275);
    const auto *nsg_276 = buffer.data(nsg + 276);
    const auto *nsg_279 = buffer.data(nsg + 279);
    const auto *nsg_280 = buffer.data(nsg + 280);
    const auto *nsg_281 = buffer.data(nsg + 281);
    const auto *nsg_282 = buffer.data(nsg + 282);
    const auto *nsg_283 = buffer.data(nsg + 283);
    const auto *nsg_284 = buffer.data(nsg + 284);
    const auto *nsg_295 = buffer.data(nsg + 295);
    const auto *nsg_296 = buffer.data(nsg + 296);
    const auto *nsg_297 = buffer.data(nsg + 297);
    const auto *nsg_298 = buffer.data(nsg + 298);
    const auto *nsg_299 = buffer.data(nsg + 299);
    const auto *nsg_300 = buffer.data(nsg + 300);
    const auto *nsg_305 = buffer.data(nsg + 305);
    const auto *nsg_309 = buffer.data(nsg + 309);
    const auto *nsg_310 = buffer.data(nsg + 310);
    const auto *nsg_311 = buffer.data(nsg + 311);
    const auto *nsg_312 = buffer.data(nsg + 312);
    const auto *nsg_314 = buffer.data(nsg + 314);
    const auto *nsg_315 = buffer.data(nsg + 315);
    const auto *nsg_318 = buffer.data(nsg + 318);
    const auto *nsg_321 = buffer.data(nsg + 321);
    const auto *nsg_325 = buffer.data(nsg + 325);
    const auto *nsg_327 = buffer.data(nsg + 327);
    const auto *nsg_328 = buffer.data(nsg + 328);
    const auto *nsg_329 = buffer.data(nsg + 329);
    const auto *nsg_335 = buffer.data(nsg + 335);
    const auto *nsg_339 = buffer.data(nsg + 339);
    const auto *nsg_340 = buffer.data(nsg + 340);
    const auto *nsg_341 = buffer.data(nsg + 341);
    const auto *nsg_342 = buffer.data(nsg + 342);
    const auto *nsg_343 = buffer.data(nsg + 343);
    const auto *nsg_344 = buffer.data(nsg + 344);
    const auto *nsg_345 = buffer.data(nsg + 345);
    const auto *nsg_348 = buffer.data(nsg + 348);

    const auto *nsh1_294 = buffer.data(nsh1 + 294);
    const auto *nsh1_297 = buffer.data(nsh1 + 297);
    const auto *nsh1_299 = buffer.data(nsh1 + 299);
    const auto *nsh1_300 = buffer.data(nsh1 + 300);
    const auto *nsh1_303 = buffer.data(nsh1 + 303);
    const auto *nsh1_314 = buffer.data(nsh1 + 314);
    const auto *nsh1_315 = buffer.data(nsh1 + 315);
    const auto *nsh1_318 = buffer.data(nsh1 + 318);
    const auto *nsh1_321 = buffer.data(nsh1 + 321);
    const auto *nsh1_330 = buffer.data(nsh1 + 330);

    const auto *osf0_176 = buffer.data(osf0 + 176);
    const auto *osf0_178 = buffer.data(osf0 + 178);
    const auto *osf0_179 = buffer.data(osf0 + 179);
    const auto *osf0_180 = buffer.data(osf0 + 180);
    const auto *osf0_183 = buffer.data(osf0 + 183);
    const auto *osf0_185 = buffer.data(osf0 + 185);
    const auto *osf0_186 = buffer.data(osf0 + 186);
    const auto *osf0_188 = buffer.data(osf0 + 188);
    const auto *osf0_189 = buffer.data(osf0 + 189);
    const auto *osf0_196 = buffer.data(osf0 + 196);
    const auto *osf0_198 = buffer.data(osf0 + 198);
    const auto *osf0_199 = buffer.data(osf0 + 199);
    const auto *osf0_200 = buffer.data(osf0 + 200);
    const auto *osf0_201 = buffer.data(osf0 + 201);
    const auto *osf0_202 = buffer.data(osf0 + 202);
    const auto *osf0_205 = buffer.data(osf0 + 205);
    const auto *osf0_206 = buffer.data(osf0 + 206);
    const auto *osf0_207 = buffer.data(osf0 + 207);
    const auto *osf0_208 = buffer.data(osf0 + 208);
    const auto *osf0_209 = buffer.data(osf0 + 209);
    const auto *osf0_210 = buffer.data(osf0 + 210);
    const auto *osf0_212 = buffer.data(osf0 + 212);
    const auto *osf0_213 = buffer.data(osf0 + 213);
    const auto *osf0_216 = buffer.data(osf0 + 216);
    const auto *osf0_217 = buffer.data(osf0 + 217);
    const auto *osf0_219 = buffer.data(osf0 + 219);
    const auto *osf0_225 = buffer.data(osf0 + 225);
    const auto *osf0_228 = buffer.data(osf0 + 228);
    const auto *osf0_229 = buffer.data(osf0 + 229);
    const auto *osf0_230 = buffer.data(osf0 + 230);
    const auto *osf0_233 = buffer.data(osf0 + 233);

    const auto *osf1_176 = buffer.data(osf1 + 176);
    const auto *osf1_178 = buffer.data(osf1 + 178);
    const auto *osf1_179 = buffer.data(osf1 + 179);
    const auto *osf1_180 = buffer.data(osf1 + 180);
    const auto *osf1_183 = buffer.data(osf1 + 183);
    const auto *osf1_185 = buffer.data(osf1 + 185);
    const auto *osf1_186 = buffer.data(osf1 + 186);
    const auto *osf1_188 = buffer.data(osf1 + 188);
    const auto *osf1_189 = buffer.data(osf1 + 189);
    const auto *osf1_196 = buffer.data(osf1 + 196);
    const auto *osf1_198 = buffer.data(osf1 + 198);
    const auto *osf1_199 = buffer.data(osf1 + 199);
    const auto *osf1_200 = buffer.data(osf1 + 200);
    const auto *osf1_201 = buffer.data(osf1 + 201);
    const auto *osf1_202 = buffer.data(osf1 + 202);
    const auto *osf1_205 = buffer.data(osf1 + 205);
    const auto *osf1_206 = buffer.data(osf1 + 206);
    const auto *osf1_207 = buffer.data(osf1 + 207);
    const auto *osf1_208 = buffer.data(osf1 + 208);
    const auto *osf1_209 = buffer.data(osf1 + 209);
    const auto *osf1_210 = buffer.data(osf1 + 210);
    const auto *osf1_212 = buffer.data(osf1 + 212);
    const auto *osf1_213 = buffer.data(osf1 + 213);
    const auto *osf1_216 = buffer.data(osf1 + 216);
    const auto *osf1_217 = buffer.data(osf1 + 217);
    const auto *osf1_219 = buffer.data(osf1 + 219);
    const auto *osf1_225 = buffer.data(osf1 + 225);
    const auto *osf1_228 = buffer.data(osf1 + 228);
    const auto *osf1_229 = buffer.data(osf1 + 229);
    const auto *osf1_230 = buffer.data(osf1 + 230);
    const auto *osf1_233 = buffer.data(osf1 + 233);

    const auto *osg_265 = buffer.data(osg + 265);
    const auto *osg_267 = buffer.data(osg + 267);
    const auto *osg_268 = buffer.data(osg + 268);
    const auto *osg_269 = buffer.data(osg + 269);
    const auto *osg_270 = buffer.data(osg + 270);
    const auto *osg_272 = buffer.data(osg + 272);
    const auto *osg_273 = buffer.data(osg + 273);
    const auto *osg_275 = buffer.data(osg + 275);
    const auto *osg_276 = buffer.data(osg + 276);
    const auto *osg_279 = buffer.data(osg + 279);
    const auto *osg_280 = buffer.data(osg + 280);
    const auto *osg_281 = buffer.data(osg + 281);
    const auto *osg_282 = buffer.data(osg + 282);
    const auto *osg_283 = buffer.data(osg + 283);
    const auto *osg_284 = buffer.data(osg + 284);
    const auto *osg_285 = buffer.data(osg + 285);
    const auto *osg_287 = buffer.data(osg + 287);
    const auto *osg_288 = buffer.data(osg + 288);
    const auto *osg_290 = buffer.data(osg + 290);
    const auto *osg_295 = buffer.data(osg + 295);
    const auto *osg_296 = buffer.data(osg + 296);
    const auto *osg_297 = buffer.data(osg + 297);
    const auto *osg_298 = buffer.data(osg + 298);
    const auto *osg_299 = buffer.data(osg + 299);
    const auto *osg_300 = buffer.data(osg + 300);
    const auto *osg_301 = buffer.data(osg + 301);
    const auto *osg_302 = buffer.data(osg + 302);
    const auto *osg_303 = buffer.data(osg + 303);
    const auto *osg_304 = buffer.data(osg + 304);
    const auto *osg_305 = buffer.data(osg + 305);
    const auto *osg_309 = buffer.data(osg + 309);
    const auto *osg_310 = buffer.data(osg + 310);
    const auto *osg_311 = buffer.data(osg + 311);
    const auto *osg_312 = buffer.data(osg + 312);
    const auto *osg_313 = buffer.data(osg + 313);
    const auto *osg_314 = buffer.data(osg + 314);
    const auto *osg_315 = buffer.data(osg + 315);
    const auto *osg_316 = buffer.data(osg + 316);
    const auto *osg_317 = buffer.data(osg + 317);
    const auto *osg_318 = buffer.data(osg + 318);
    const auto *osg_320 = buffer.data(osg + 320);
    const auto *osg_321 = buffer.data(osg + 321);
    const auto *osg_325 = buffer.data(osg + 325);
    const auto *osg_326 = buffer.data(osg + 326);
    const auto *osg_327 = buffer.data(osg + 327);
    const auto *osg_328 = buffer.data(osg + 328);
    const auto *osg_329 = buffer.data(osg + 329);
    const auto *osg_330 = buffer.data(osg + 330);
    const auto *osg_332 = buffer.data(osg + 332);
    const auto *osg_333 = buffer.data(osg + 333);
    const auto *osg_335 = buffer.data(osg + 335);
    const auto *osg_339 = buffer.data(osg + 339);
    const auto *osg_340 = buffer.data(osg + 340);
    const auto *osg_341 = buffer.data(osg + 341);
    const auto *osg_342 = buffer.data(osg + 342);
    const auto *osg_343 = buffer.data(osg + 343);
    const auto *osg_344 = buffer.data(osg + 344);
    const auto *osg_345 = buffer.data(osg + 345);
    const auto *osg_347 = buffer.data(osg + 347);
    const auto *osg_348 = buffer.data(osg + 348);

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, nsg_190, nsg_267, nsg_268, \
                         nsg_269, osf0_176, osf1_176, osg_265, osg_267, osg_268, \
                         osg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_19 * nsg_267[k]
                   + f_3 * pc_x[k] * osg_267[k];

        t_370[k] = f_19 * nsg_268[k]
                   + f_3 * pc_x[k] * osg_268[k];

        t_371[k] = f_19 * nsg_269[k]
                   + f_3 * pc_x[k] * osg_269[k];

        t_372[k] = f_11 * nsg_190[k]
                   + f_1 * osf0_176[k]
                   - f_2 * osf1_176[k]
                   + f_3 * pc_y[k] * osg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, nsg_175, nsg_192, nsg_193, osf0_178, \
                         osf0_179, osf1_178, osf1_179, osg_265, osg_267, \
                         osg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * nsg_175[k]
                   + f_3 * pc_z[k] * osg_265[k];

        t_374[k] = f_11 * nsg_192[k]
                   + f_6 * osf0_178[k]
                   - f_7 * osf1_178[k]
                   + f_3 * pc_y[k] * osg_267[k];

        t_375[k] = f_11 * nsg_193[k]
                   + f_4 * osf0_179[k]
                   - f_5 * osf1_179[k]
                   + f_3 * pc_y[k] * osg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, nsg_179, nsg_194, nsg_270, \
                         osf0_179, osf0_180, osf1_179, osf1_180, osg_269, \
                         osg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * nsg_194[k]
                   + f_3 * pc_y[k] * osg_269[k];

        t_377[k] = f_10 * nsg_179[k]
                   + f_1 * osf0_179[k]
                   - f_2 * osf1_179[k]
                   + f_3 * pc_z[k] * osg_269[k];

        t_378[k] = f_19 * nsg_270[k]
                   + f_1 * osf0_180[k]
                   - f_2 * osf1_180[k]
                   + f_3 * pc_x[k] * osg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, nsg_180, nsg_195, \
                         nsg_197, nsg_273, osf0_183, osf1_183, osg_270, osg_272, \
                         osg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * nsg_195[k]
                   + f_3 * pc_y[k] * osg_270[k];

        t_380[k] = f_11 * nsg_180[k]
                   + f_3 * pc_z[k] * osg_270[k];

        t_381[k] = f_19 * nsg_273[k]
                   + f_6 * osf0_183[k]
                   - f_7 * osf1_183[k]
                   + f_3 * pc_x[k] * osg_273[k];

        t_382[k] = f_10 * nsg_197[k]
                   + f_3 * pc_y[k] * osg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, nsg_183, nsg_275, nsg_276, osf0_185, \
                         osf0_186, osf1_185, osf1_186, osg_273, osg_275, \
                         osg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_19 * nsg_275[k]
                   + f_6 * osf0_185[k]
                   - f_7 * osf1_185[k]
                   + f_3 * pc_x[k] * osg_275[k];

        t_384[k] = f_19 * nsg_276[k]
                   + f_4 * osf0_186[k]
                   - f_5 * osf1_186[k]
                   + f_3 * pc_x[k] * osg_276[k];

        t_385[k] = f_11 * nsg_183[k]
                   + f_3 * pc_z[k] * osg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, nsg_200, nsg_279, nsg_280, \
                         nsg_281, osf0_189, osf1_189, osg_275, osg_279, osg_280, \
                         osg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * nsg_200[k]
                   + f_3 * pc_y[k] * osg_275[k];

        t_387[k] = f_19 * nsg_279[k]
                   + f_4 * osf0_189[k]
                   - f_5 * osf1_189[k]
                   + f_3 * pc_x[k] * osg_279[k];

        t_388[k] = f_19 * nsg_280[k]
                   + f_3 * pc_x[k] * osg_280[k];

        t_389[k] = f_19 * nsg_281[k]
                   + f_3 * pc_x[k] * osg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, nsg_205, nsg_282, nsg_283, \
                         nsg_284, osf0_186, osf1_186, osg_280, osg_282, osg_283, \
                         osg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_19 * nsg_282[k]
                   + f_3 * pc_x[k] * osg_282[k];

        t_391[k] = f_19 * nsg_283[k]
                   + f_3 * pc_x[k] * osg_283[k];

        t_392[k] = f_19 * nsg_284[k]
                   + f_3 * pc_x[k] * osg_284[k];

        t_393[k] = f_10 * nsg_205[k]
                   + f_1 * osf0_186[k]
                   - f_2 * osf1_186[k]
                   + f_3 * pc_y[k] * osg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, nsg_190, nsg_207, nsg_208, osf0_188, \
                         osf0_189, osf1_188, osf1_189, osg_280, osg_282, \
                         osg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * nsg_190[k]
                   + f_3 * pc_z[k] * osg_280[k];

        t_395[k] = f_10 * nsg_207[k]
                   + f_6 * osf0_188[k]
                   - f_7 * osf1_188[k]
                   + f_3 * pc_y[k] * osg_282[k];

        t_396[k] = f_10 * nsg_208[k]
                   + f_4 * osf0_189[k]
                   - f_5 * osf1_189[k]
                   + f_3 * pc_y[k] * osg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pc_y, pc_z, nsh0_294, nsg_194, \
                         nsg_209, nsg_210, nsh1_294, osf0_189, osf1_189, osg_284, \
                         osg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * nsg_209[k]
                   + f_3 * pc_y[k] * osg_284[k];

        t_398[k] = f_11 * nsg_194[k]
                   + f_1 * osf0_189[k]
                   - f_2 * osf1_189[k]
                   + f_3 * pc_z[k] * osg_284[k];

        t_399[k] = pa_y[k] * nsh0_294[k]
                   - f_8 * pc_y[k] * nsh1_294[k];

        t_400[k] = f_9 * nsg_210[k]
                   + f_3 * pc_y[k] * osg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pc_y, pc_z, nsh0_297, nsh0_299, \
                         nsg_195, nsg_211, nsg_212, nsh1_297, nsh1_299, osg_285, \
                         osg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * nsg_195[k]
                   + f_3 * pc_z[k] * osg_285[k];

        t_402[k] = pa_y[k] * nsh0_297[k]
                   + f_10 * nsg_211[k]
                   - f_8 * pc_y[k] * nsh1_297[k];

        t_403[k] = f_9 * nsg_212[k]
                   + f_3 * pc_y[k] * osg_287[k];

        t_404[k] = pa_y[k] * nsh0_299[k]
                   - f_8 * pc_y[k] * nsh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_y, pc_y, pc_z, nsh0_300, nsh0_303, \
                         nsg_198, nsg_213, nsg_215, nsh1_300, nsh1_303, osg_288, \
                         osg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_y[k] * nsh0_300[k]
                   + f_11 * nsg_213[k]
                   - f_8 * pc_y[k] * nsh1_300[k];

        t_406[k] = f_18 * nsg_198[k]
                   + f_3 * pc_z[k] * osg_288[k];

        t_407[k] = f_9 * nsg_215[k]
                   + f_3 * pc_y[k] * osg_290[k];

        t_408[k] = pa_y[k] * nsh0_303[k]
                   - f_8 * pc_y[k] * nsh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, nsg_295, nsg_296, nsg_297, \
                         nsg_298, nsg_299, osg_295, osg_296, osg_297, osg_298, \
                         osg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_19 * nsg_295[k]
                   + f_3 * pc_x[k] * osg_295[k];

        t_410[k] = f_19 * nsg_296[k]
                   + f_3 * pc_x[k] * osg_296[k];

        t_411[k] = f_19 * nsg_297[k]
                   + f_3 * pc_x[k] * osg_297[k];

        t_412[k] = f_19 * nsg_298[k]
                   + f_3 * pc_x[k] * osg_298[k];

        t_413[k] = f_19 * nsg_299[k]
                   + f_3 * pc_x[k] * osg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, nsg_205, nsg_220, nsg_222, osf0_196, \
                         osf0_198, osf1_196, osf1_198, osg_295, \
                         osg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * nsg_220[k]
                   + f_1 * osf0_196[k]
                   - f_2 * osf1_196[k]
                   + f_3 * pc_y[k] * osg_295[k];

        t_415[k] = f_18 * nsg_205[k]
                   + f_3 * pc_z[k] * osg_295[k];

        t_416[k] = f_9 * nsg_222[k]
                   + f_6 * osf0_198[k]
                   - f_7 * osf1_198[k]
                   + f_3 * pc_y[k] * osg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, nsh0_314, nsg_223, nsg_224, \
                         nsh1_314, osf0_199, osf1_199, osg_298, \
                         osg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * nsg_223[k]
                   + f_4 * osf0_199[k]
                   - f_5 * osf1_199[k]
                   + f_3 * pc_y[k] * osg_298[k];

        t_418[k] = f_9 * nsg_224[k]
                   + f_3 * pc_y[k] * osg_299[k];

        t_419[k] = pa_y[k] * nsh0_314[k]
                   - f_8 * pc_y[k] * nsh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, nsg_210, \
                         nsg_300, osf0_200, osf1_200, osg_300, osg_301, \
                         osg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_19 * nsg_300[k]
                   + f_1 * osf0_200[k]
                   - f_2 * osf1_200[k]
                   + f_3 * pc_x[k] * osg_300[k];

        t_421[k] = f_3 * pc_y[k] * osg_300[k];

        t_422[k] = f_20 * nsg_210[k]
                   + f_3 * pc_z[k] * osg_300[k];

        t_423[k] = f_4 * osf0_200[k]
                   - f_5 * osf1_200[k]
                   + f_3 * pc_y[k] * osg_301[k];

        t_424[k] = f_3 * pc_y[k] * osg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, nsg_305, osf0_201, osf0_202, \
                         osf0_205, osf1_201, osf1_202, osf1_205, osg_303, osg_304, \
                         osg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_19 * nsg_305[k]
                   + f_6 * osf0_205[k]
                   - f_7 * osf1_205[k]
                   + f_3 * pc_x[k] * osg_305[k];

        t_426[k] = f_6 * osf0_201[k]
                   - f_7 * osf1_201[k]
                   + f_3 * pc_y[k] * osg_303[k];

        t_427[k] = f_4 * osf0_202[k]
                   - f_5 * osf1_202[k]
                   + f_3 * pc_y[k] * osg_304[k];

        t_428[k] = f_3 * pc_y[k] * osg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, nsg_309, nsg_310, nsg_311, nsg_312, \
                         osf0_209, osf1_209, osg_309, osg_310, osg_311, \
                         osg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_19 * nsg_309[k]
                   + f_4 * osf0_209[k]
                   - f_5 * osf1_209[k]
                   + f_3 * pc_x[k] * osg_309[k];

        t_430[k] = f_19 * nsg_310[k]
                   + f_3 * pc_x[k] * osg_310[k];

        t_431[k] = f_19 * nsg_311[k]
                   + f_3 * pc_x[k] * osg_311[k];

        t_432[k] = f_19 * nsg_312[k]
                   + f_3 * pc_x[k] * osg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, nsg_314, osf0_206, osf0_207, \
                         osf1_206, osf1_207, osg_309, osg_310, osg_311, \
                         osg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_y[k] * osg_309[k];

        t_434[k] = f_19 * nsg_314[k]
                   + f_3 * pc_x[k] * osg_314[k];

        t_435[k] = f_1 * osf0_206[k]
                   - f_2 * osf1_206[k]
                   + f_3 * pc_y[k] * osg_310[k];

        t_436[k] = f_13 * osf0_207[k]
                   - f_14 * osf1_207[k]
                   + f_3 * pc_y[k] * osg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_y, pc_z, nsg_224, osf0_208, osf0_209, \
                         osf1_208, osf1_209, osg_312, osg_313, \
                         osg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_6 * osf0_208[k]
                   - f_7 * osf1_208[k]
                   + f_3 * pc_y[k] * osg_312[k];

        t_438[k] = f_4 * osf0_209[k]
                   - f_5 * osf1_209[k]
                   + f_3 * pc_y[k] * osg_313[k];

        t_439[k] = f_3 * pc_y[k] * osg_314[k];

        t_440[k] = f_20 * nsg_224[k]
                   + f_1 * osf0_209[k]
                   - f_2 * osf1_209[k]
                   + f_3 * pc_z[k] * osg_314[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, nsg_225, nsg_315, \
                         nsg_318, osf0_210, osf0_213, osf1_210, osf1_213, osg_315, \
                         osg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_20 * nsg_315[k]
                   + f_1 * osf0_210[k]
                   - f_2 * osf1_210[k]
                   + f_3 * pc_x[k] * osg_315[k];

        t_442[k] = f_19 * nsg_225[k]
                   + f_3 * pc_y[k] * osg_315[k];

        t_443[k] = f_3 * pc_z[k] * osg_315[k];

        t_444[k] = f_20 * nsg_318[k]
                   + f_6 * osf0_213[k]
                   - f_7 * osf1_213[k]
                   + f_3 * pc_x[k] * osg_318[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_z, nsg_321, osf0_210, osf0_216, \
                         osf1_210, osf1_216, osg_316, osg_317, osg_318, \
                         osg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_z[k] * osg_316[k];

        t_446[k] = f_4 * osf0_210[k]
                   - f_5 * osf1_210[k]
                   + f_3 * pc_z[k] * osg_317[k];

        t_447[k] = f_20 * nsg_321[k]
                   + f_4 * osf0_216[k]
                   - f_5 * osf1_216[k]
                   + f_3 * pc_x[k] * osg_321[k];

        t_448[k] = f_3 * pc_z[k] * osg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, nsg_230, nsg_325, \
                         osf0_212, osf1_212, osg_320, osg_321, \
                         osg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_19 * nsg_230[k]
                   + f_3 * pc_y[k] * osg_320[k];

        t_450[k] = f_6 * osf0_212[k]
                   - f_7 * osf1_212[k]
                   + f_3 * pc_z[k] * osg_320[k];

        t_451[k] = f_20 * nsg_325[k]
                   + f_3 * pc_x[k] * osg_325[k];

        t_452[k] = f_3 * pc_z[k] * osg_321[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, nsg_235, nsg_327, nsg_328, \
                         nsg_329, osf0_216, osf1_216, osg_325, osg_327, osg_328, \
                         osg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_20 * nsg_327[k]
                   + f_3 * pc_x[k] * osg_327[k];

        t_454[k] = f_20 * nsg_328[k]
                   + f_3 * pc_x[k] * osg_328[k];

        t_455[k] = f_20 * nsg_329[k]
                   + f_3 * pc_x[k] * osg_329[k];

        t_456[k] = f_19 * nsg_235[k]
                   + f_1 * osf0_216[k]
                   - f_2 * osf1_216[k]
                   + f_3 * pc_y[k] * osg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_y, pc_z, nsg_239, osf0_216, osf0_217, \
                         osf1_216, osf1_217, osg_325, osg_326, osg_327, \
                         osg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * osg_325[k];

        t_458[k] = f_4 * osf0_216[k]
                   - f_5 * osf1_216[k]
                   + f_3 * pc_z[k] * osg_326[k];

        t_459[k] = f_6 * osf0_217[k]
                   - f_7 * osf1_217[k]
                   + f_3 * pc_z[k] * osg_327[k];

        t_460[k] = f_19 * nsg_239[k]
                   + f_3 * pc_y[k] * osg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_z, pc_y, pc_z, nsh0_315, nsg_225, \
                         nsg_240, nsh1_315, osf0_219, osf1_219, osg_329, \
                         osg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_1 * osf0_219[k]
                   - f_2 * osf1_219[k]
                   + f_3 * pc_z[k] * osg_329[k];

        t_462[k] = pa_z[k] * nsh0_315[k]
                   - f_8 * pc_z[k] * nsh1_315[k];

        t_463[k] = f_20 * nsg_240[k]
                   + f_3 * pc_y[k] * osg_330[k];

        t_464[k] = f_9 * nsg_225[k]
                   + f_3 * pc_z[k] * osg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_z, pc_x, pc_y, pc_z, nsh0_318, nsg_242, \
                         nsg_335, nsh1_318, osf0_225, osf1_225, osg_332, \
                         osg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * nsh0_318[k]
                   - f_8 * pc_z[k] * nsh1_318[k];

        t_466[k] = f_20 * nsg_242[k]
                   + f_3 * pc_y[k] * osg_332[k];

        t_467[k] = f_20 * nsg_335[k]
                   + f_6 * osf0_225[k]
                   - f_7 * osf1_225[k]
                   + f_3 * pc_x[k] * osg_335[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_z, pc_y, pc_z, nsh0_321, nsg_228, nsg_245, \
                         nsh1_321, osg_333, osg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_z[k] * nsh0_321[k]
                   - f_8 * pc_z[k] * nsh1_321[k];

        t_469[k] = f_9 * nsg_228[k]
                   + f_3 * pc_z[k] * osg_333[k];

        t_470[k] = f_20 * nsg_245[k]
                   + f_3 * pc_y[k] * osg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, nsg_339, nsg_340, nsg_341, nsg_342, \
                         osf0_229, osf1_229, osg_339, osg_340, osg_341, \
                         osg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_20 * nsg_339[k]
                   + f_4 * osf0_229[k]
                   - f_5 * osf1_229[k]
                   + f_3 * pc_x[k] * osg_339[k];

        t_472[k] = f_20 * nsg_340[k]
                   + f_3 * pc_x[k] * osg_340[k];

        t_473[k] = f_20 * nsg_341[k]
                   + f_3 * pc_x[k] * osg_341[k];

        t_474[k] = f_20 * nsg_342[k]
                   + f_3 * pc_x[k] * osg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pc_x, pc_z, nsh0_330, nsg_235, \
                         nsg_343, nsg_344, nsh1_330, osg_340, osg_343, \
                         osg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_20 * nsg_343[k]
                   + f_3 * pc_x[k] * osg_343[k];

        t_476[k] = f_20 * nsg_344[k]
                   + f_3 * pc_x[k] * osg_344[k];

        t_477[k] = pa_z[k] * nsh0_330[k]
                   - f_8 * pc_z[k] * nsh1_330[k];

        t_478[k] = f_9 * nsg_235[k]
                   + f_3 * pc_z[k] * osg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_y, nsg_252, nsg_253, nsg_254, osf0_228, \
                         osf0_229, osf1_228, osf1_229, osg_342, osg_343, \
                         osg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_20 * nsg_252[k]
                   + f_6 * osf0_228[k]
                   - f_7 * osf1_228[k]
                   + f_3 * pc_y[k] * osg_342[k];

        t_480[k] = f_20 * nsg_253[k]
                   + f_4 * osf0_229[k]
                   - f_5 * osf1_229[k]
                   + f_3 * pc_y[k] * osg_343[k];

        t_481[k] = f_20 * nsg_254[k]
                   + f_3 * pc_y[k] * osg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, nsg_239, nsg_255, nsg_345, \
                         osf0_229, osf0_230, osf1_229, osf1_230, osg_344, \
                         osg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * nsg_239[k]
                   + f_1 * osf0_229[k]
                   - f_2 * osf1_229[k]
                   + f_3 * pc_z[k] * osg_344[k];

        t_483[k] = f_20 * nsg_345[k]
                   + f_1 * osf0_230[k]
                   - f_2 * osf1_230[k]
                   + f_3 * pc_x[k] * osg_345[k];

        t_484[k] = f_18 * nsg_255[k]
                   + f_3 * pc_y[k] * osg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, nsg_240, nsg_257, nsg_348, \
                         osf0_233, osf1_233, osg_345, osg_347, \
                         osg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * nsg_240[k]
                   + f_3 * pc_z[k] * osg_345[k];

        t_486[k] = f_20 * nsg_348[k]
                   + f_6 * osf0_233[k]
                   - f_7 * osf1_233[k]
                   + f_3 * pc_x[k] * osg_348[k];

        t_487[k] = f_18 * nsg_257[k]
                   + f_3 * pc_y[k] * osg_347[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_420 = buffer.data(nsh0 + 420);
    const auto *nsh0_423 = buffer.data(nsh0 + 423);
    const auto *nsh0_425 = buffer.data(nsh0 + 425);
    const auto *nsh0_426 = buffer.data(nsh0 + 426);
    const auto *nsh0_429 = buffer.data(nsh0 + 429);
    const auto *nsh0_440 = buffer.data(nsh0 + 440);

    const auto *nsg_243 = buffer.data(nsg + 243);
    const auto *nsg_250 = buffer.data(nsg + 250);
    const auto *nsg_254 = buffer.data(nsg + 254);
    const auto *nsg_255 = buffer.data(nsg + 255);
    const auto *nsg_258 = buffer.data(nsg + 258);
    const auto *nsg_260 = buffer.data(nsg + 260);
    const auto *nsg_265 = buffer.data(nsg + 265);
    const auto *nsg_267 = buffer.data(nsg + 267);
    const auto *nsg_268 = buffer.data(nsg + 268);
    const auto *nsg_269 = buffer.data(nsg + 269);
    const auto *nsg_270 = buffer.data(nsg + 270);
    const auto *nsg_272 = buffer.data(nsg + 272);
    const auto *nsg_273 = buffer.data(nsg + 273);
    const auto *nsg_275 = buffer.data(nsg + 275);
    const auto *nsg_280 = buffer.data(nsg + 280);
    const auto *nsg_282 = buffer.data(nsg + 282);
    const auto *nsg_283 = buffer.data(nsg + 283);
    const auto *nsg_284 = buffer.data(nsg + 284);
    const auto *nsg_285 = buffer.data(nsg + 285);
    const auto *nsg_287 = buffer.data(nsg + 287);
    const auto *nsg_288 = buffer.data(nsg + 288);
    const auto *nsg_290 = buffer.data(nsg + 290);
    const auto *nsg_295 = buffer.data(nsg + 295);
    const auto *nsg_297 = buffer.data(nsg + 297);
    const auto *nsg_298 = buffer.data(nsg + 298);
    const auto *nsg_299 = buffer.data(nsg + 299);
    const auto *nsg_300 = buffer.data(nsg + 300);
    const auto *nsg_301 = buffer.data(nsg + 301);
    const auto *nsg_302 = buffer.data(nsg + 302);
    const auto *nsg_303 = buffer.data(nsg + 303);
    const auto *nsg_305 = buffer.data(nsg + 305);
    const auto *nsg_310 = buffer.data(nsg + 310);
    const auto *nsg_312 = buffer.data(nsg + 312);
    const auto *nsg_313 = buffer.data(nsg + 313);
    const auto *nsg_314 = buffer.data(nsg + 314);
    const auto *nsg_315 = buffer.data(nsg + 315);
    const auto *nsg_320 = buffer.data(nsg + 320);
    const auto *nsg_325 = buffer.data(nsg + 325);
    const auto *nsg_350 = buffer.data(nsg + 350);
    const auto *nsg_351 = buffer.data(nsg + 351);
    const auto *nsg_354 = buffer.data(nsg + 354);
    const auto *nsg_355 = buffer.data(nsg + 355);
    const auto *nsg_356 = buffer.data(nsg + 356);
    const auto *nsg_357 = buffer.data(nsg + 357);
    const auto *nsg_358 = buffer.data(nsg + 358);
    const auto *nsg_359 = buffer.data(nsg + 359);
    const auto *nsg_360 = buffer.data(nsg + 360);
    const auto *nsg_363 = buffer.data(nsg + 363);
    const auto *nsg_365 = buffer.data(nsg + 365);
    const auto *nsg_366 = buffer.data(nsg + 366);
    const auto *nsg_369 = buffer.data(nsg + 369);
    const auto *nsg_370 = buffer.data(nsg + 370);
    const auto *nsg_371 = buffer.data(nsg + 371);
    const auto *nsg_372 = buffer.data(nsg + 372);
    const auto *nsg_373 = buffer.data(nsg + 373);
    const auto *nsg_374 = buffer.data(nsg + 374);
    const auto *nsg_375 = buffer.data(nsg + 375);
    const auto *nsg_378 = buffer.data(nsg + 378);
    const auto *nsg_380 = buffer.data(nsg + 380);
    const auto *nsg_381 = buffer.data(nsg + 381);
    const auto *nsg_384 = buffer.data(nsg + 384);
    const auto *nsg_385 = buffer.data(nsg + 385);
    const auto *nsg_386 = buffer.data(nsg + 386);
    const auto *nsg_387 = buffer.data(nsg + 387);
    const auto *nsg_388 = buffer.data(nsg + 388);
    const auto *nsg_389 = buffer.data(nsg + 389);
    const auto *nsg_400 = buffer.data(nsg + 400);
    const auto *nsg_401 = buffer.data(nsg + 401);
    const auto *nsg_402 = buffer.data(nsg + 402);
    const auto *nsg_403 = buffer.data(nsg + 403);
    const auto *nsg_404 = buffer.data(nsg + 404);
    const auto *nsg_405 = buffer.data(nsg + 405);
    const auto *nsg_410 = buffer.data(nsg + 410);
    const auto *nsg_414 = buffer.data(nsg + 414);
    const auto *nsg_415 = buffer.data(nsg + 415);
    const auto *nsg_416 = buffer.data(nsg + 416);
    const auto *nsg_417 = buffer.data(nsg + 417);
    const auto *nsg_419 = buffer.data(nsg + 419);
    const auto *nsg_420 = buffer.data(nsg + 420);
    const auto *nsg_423 = buffer.data(nsg + 423);
    const auto *nsg_426 = buffer.data(nsg + 426);
    const auto *nsg_430 = buffer.data(nsg + 430);
    const auto *nsg_432 = buffer.data(nsg + 432);
    const auto *nsg_433 = buffer.data(nsg + 433);
    const auto *nsg_434 = buffer.data(nsg + 434);

    const auto *nsh1_420 = buffer.data(nsh1 + 420);
    const auto *nsh1_423 = buffer.data(nsh1 + 423);
    const auto *nsh1_425 = buffer.data(nsh1 + 425);
    const auto *nsh1_426 = buffer.data(nsh1 + 426);
    const auto *nsh1_429 = buffer.data(nsh1 + 429);
    const auto *nsh1_440 = buffer.data(nsh1 + 440);

    const auto *osf0_235 = buffer.data(osf0 + 235);
    const auto *osf0_236 = buffer.data(osf0 + 236);
    const auto *osf0_238 = buffer.data(osf0 + 238);
    const auto *osf0_239 = buffer.data(osf0 + 239);
    const auto *osf0_240 = buffer.data(osf0 + 240);
    const auto *osf0_243 = buffer.data(osf0 + 243);
    const auto *osf0_245 = buffer.data(osf0 + 245);
    const auto *osf0_246 = buffer.data(osf0 + 246);
    const auto *osf0_248 = buffer.data(osf0 + 248);
    const auto *osf0_249 = buffer.data(osf0 + 249);
    const auto *osf0_250 = buffer.data(osf0 + 250);
    const auto *osf0_253 = buffer.data(osf0 + 253);
    const auto *osf0_255 = buffer.data(osf0 + 255);
    const auto *osf0_256 = buffer.data(osf0 + 256);
    const auto *osf0_258 = buffer.data(osf0 + 258);
    const auto *osf0_259 = buffer.data(osf0 + 259);
    const auto *osf0_266 = buffer.data(osf0 + 266);
    const auto *osf0_268 = buffer.data(osf0 + 268);
    const auto *osf0_269 = buffer.data(osf0 + 269);
    const auto *osf0_270 = buffer.data(osf0 + 270);
    const auto *osf0_271 = buffer.data(osf0 + 271);
    const auto *osf0_272 = buffer.data(osf0 + 272);
    const auto *osf0_275 = buffer.data(osf0 + 275);
    const auto *osf0_276 = buffer.data(osf0 + 276);
    const auto *osf0_277 = buffer.data(osf0 + 277);
    const auto *osf0_278 = buffer.data(osf0 + 278);
    const auto *osf0_279 = buffer.data(osf0 + 279);
    const auto *osf0_280 = buffer.data(osf0 + 280);
    const auto *osf0_282 = buffer.data(osf0 + 282);
    const auto *osf0_283 = buffer.data(osf0 + 283);
    const auto *osf0_286 = buffer.data(osf0 + 286);

    const auto *osf1_235 = buffer.data(osf1 + 235);
    const auto *osf1_236 = buffer.data(osf1 + 236);
    const auto *osf1_238 = buffer.data(osf1 + 238);
    const auto *osf1_239 = buffer.data(osf1 + 239);
    const auto *osf1_240 = buffer.data(osf1 + 240);
    const auto *osf1_243 = buffer.data(osf1 + 243);
    const auto *osf1_245 = buffer.data(osf1 + 245);
    const auto *osf1_246 = buffer.data(osf1 + 246);
    const auto *osf1_248 = buffer.data(osf1 + 248);
    const auto *osf1_249 = buffer.data(osf1 + 249);
    const auto *osf1_250 = buffer.data(osf1 + 250);
    const auto *osf1_253 = buffer.data(osf1 + 253);
    const auto *osf1_255 = buffer.data(osf1 + 255);
    const auto *osf1_256 = buffer.data(osf1 + 256);
    const auto *osf1_258 = buffer.data(osf1 + 258);
    const auto *osf1_259 = buffer.data(osf1 + 259);
    const auto *osf1_266 = buffer.data(osf1 + 266);
    const auto *osf1_268 = buffer.data(osf1 + 268);
    const auto *osf1_269 = buffer.data(osf1 + 269);
    const auto *osf1_270 = buffer.data(osf1 + 270);
    const auto *osf1_271 = buffer.data(osf1 + 271);
    const auto *osf1_272 = buffer.data(osf1 + 272);
    const auto *osf1_275 = buffer.data(osf1 + 275);
    const auto *osf1_276 = buffer.data(osf1 + 276);
    const auto *osf1_277 = buffer.data(osf1 + 277);
    const auto *osf1_278 = buffer.data(osf1 + 278);
    const auto *osf1_279 = buffer.data(osf1 + 279);
    const auto *osf1_280 = buffer.data(osf1 + 280);
    const auto *osf1_282 = buffer.data(osf1 + 282);
    const auto *osf1_283 = buffer.data(osf1 + 283);
    const auto *osf1_286 = buffer.data(osf1 + 286);

    const auto *osg_348 = buffer.data(osg + 348);
    const auto *osg_350 = buffer.data(osg + 350);
    const auto *osg_351 = buffer.data(osg + 351);
    const auto *osg_354 = buffer.data(osg + 354);
    const auto *osg_355 = buffer.data(osg + 355);
    const auto *osg_356 = buffer.data(osg + 356);
    const auto *osg_357 = buffer.data(osg + 357);
    const auto *osg_358 = buffer.data(osg + 358);
    const auto *osg_359 = buffer.data(osg + 359);
    const auto *osg_360 = buffer.data(osg + 360);
    const auto *osg_362 = buffer.data(osg + 362);
    const auto *osg_363 = buffer.data(osg + 363);
    const auto *osg_365 = buffer.data(osg + 365);
    const auto *osg_366 = buffer.data(osg + 366);
    const auto *osg_369 = buffer.data(osg + 369);
    const auto *osg_370 = buffer.data(osg + 370);
    const auto *osg_371 = buffer.data(osg + 371);
    const auto *osg_372 = buffer.data(osg + 372);
    const auto *osg_373 = buffer.data(osg + 373);
    const auto *osg_374 = buffer.data(osg + 374);
    const auto *osg_375 = buffer.data(osg + 375);
    const auto *osg_377 = buffer.data(osg + 377);
    const auto *osg_378 = buffer.data(osg + 378);
    const auto *osg_380 = buffer.data(osg + 380);
    const auto *osg_381 = buffer.data(osg + 381);
    const auto *osg_384 = buffer.data(osg + 384);
    const auto *osg_385 = buffer.data(osg + 385);
    const auto *osg_386 = buffer.data(osg + 386);
    const auto *osg_387 = buffer.data(osg + 387);
    const auto *osg_388 = buffer.data(osg + 388);
    const auto *osg_389 = buffer.data(osg + 389);
    const auto *osg_390 = buffer.data(osg + 390);
    const auto *osg_392 = buffer.data(osg + 392);
    const auto *osg_393 = buffer.data(osg + 393);
    const auto *osg_395 = buffer.data(osg + 395);
    const auto *osg_400 = buffer.data(osg + 400);
    const auto *osg_401 = buffer.data(osg + 401);
    const auto *osg_402 = buffer.data(osg + 402);
    const auto *osg_403 = buffer.data(osg + 403);
    const auto *osg_404 = buffer.data(osg + 404);
    const auto *osg_405 = buffer.data(osg + 405);
    const auto *osg_406 = buffer.data(osg + 406);
    const auto *osg_407 = buffer.data(osg + 407);
    const auto *osg_408 = buffer.data(osg + 408);
    const auto *osg_409 = buffer.data(osg + 409);
    const auto *osg_410 = buffer.data(osg + 410);
    const auto *osg_414 = buffer.data(osg + 414);
    const auto *osg_415 = buffer.data(osg + 415);
    const auto *osg_416 = buffer.data(osg + 416);
    const auto *osg_417 = buffer.data(osg + 417);
    const auto *osg_418 = buffer.data(osg + 418);
    const auto *osg_419 = buffer.data(osg + 419);
    const auto *osg_420 = buffer.data(osg + 420);
    const auto *osg_421 = buffer.data(osg + 421);
    const auto *osg_422 = buffer.data(osg + 422);
    const auto *osg_423 = buffer.data(osg + 423);
    const auto *osg_425 = buffer.data(osg + 425);
    const auto *osg_426 = buffer.data(osg + 426);
    const auto *osg_430 = buffer.data(osg + 430);
    const auto *osg_432 = buffer.data(osg + 432);
    const auto *osg_433 = buffer.data(osg + 433);
    const auto *osg_434 = buffer.data(osg + 434);

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, nsg_243, nsg_350, nsg_351, osf0_235, \
                         osf0_236, osf1_235, osf1_236, osg_348, osg_350, \
                         osg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_20 * nsg_350[k]
                   + f_6 * osf0_235[k]
                   - f_7 * osf1_235[k]
                   + f_3 * pc_x[k] * osg_350[k];

        t_489[k] = f_20 * nsg_351[k]
                   + f_4 * osf0_236[k]
                   - f_5 * osf1_236[k]
                   + f_3 * pc_x[k] * osg_351[k];

        t_490[k] = f_10 * nsg_243[k]
                   + f_3 * pc_z[k] * osg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, nsg_260, nsg_354, nsg_355, \
                         nsg_356, osf0_239, osf1_239, osg_350, osg_354, osg_355, \
                         osg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_18 * nsg_260[k]
                   + f_3 * pc_y[k] * osg_350[k];

        t_492[k] = f_20 * nsg_354[k]
                   + f_4 * osf0_239[k]
                   - f_5 * osf1_239[k]
                   + f_3 * pc_x[k] * osg_354[k];

        t_493[k] = f_20 * nsg_355[k]
                   + f_3 * pc_x[k] * osg_355[k];

        t_494[k] = f_20 * nsg_356[k]
                   + f_3 * pc_x[k] * osg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, nsg_265, nsg_357, nsg_358, \
                         nsg_359, osf0_236, osf1_236, osg_355, osg_357, osg_358, \
                         osg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_20 * nsg_357[k]
                   + f_3 * pc_x[k] * osg_357[k];

        t_496[k] = f_20 * nsg_358[k]
                   + f_3 * pc_x[k] * osg_358[k];

        t_497[k] = f_20 * nsg_359[k]
                   + f_3 * pc_x[k] * osg_359[k];

        t_498[k] = f_18 * nsg_265[k]
                   + f_1 * osf0_236[k]
                   - f_2 * osf1_236[k]
                   + f_3 * pc_y[k] * osg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, nsg_250, nsg_267, nsg_268, osf0_238, \
                         osf0_239, osf1_238, osf1_239, osg_355, osg_357, \
                         osg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * nsg_250[k]
                   + f_3 * pc_z[k] * osg_355[k];

        t_500[k] = f_18 * nsg_267[k]
                   + f_6 * osf0_238[k]
                   - f_7 * osf1_238[k]
                   + f_3 * pc_y[k] * osg_357[k];

        t_501[k] = f_18 * nsg_268[k]
                   + f_4 * osf0_239[k]
                   - f_5 * osf1_239[k]
                   + f_3 * pc_y[k] * osg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, nsg_254, nsg_269, nsg_360, \
                         osf0_239, osf0_240, osf1_239, osf1_240, osg_359, \
                         osg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_18 * nsg_269[k]
                   + f_3 * pc_y[k] * osg_359[k];

        t_503[k] = f_10 * nsg_254[k]
                   + f_1 * osf0_239[k]
                   - f_2 * osf1_239[k]
                   + f_3 * pc_z[k] * osg_359[k];

        t_504[k] = f_20 * nsg_360[k]
                   + f_1 * osf0_240[k]
                   - f_2 * osf1_240[k]
                   + f_3 * pc_x[k] * osg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, nsg_255, nsg_270, \
                         nsg_272, nsg_363, osf0_243, osf1_243, osg_360, osg_362, \
                         osg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * nsg_270[k]
                   + f_3 * pc_y[k] * osg_360[k];

        t_506[k] = f_11 * nsg_255[k]
                   + f_3 * pc_z[k] * osg_360[k];

        t_507[k] = f_20 * nsg_363[k]
                   + f_6 * osf0_243[k]
                   - f_7 * osf1_243[k]
                   + f_3 * pc_x[k] * osg_363[k];

        t_508[k] = f_11 * nsg_272[k]
                   + f_3 * pc_y[k] * osg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, nsg_258, nsg_365, nsg_366, osf0_245, \
                         osf0_246, osf1_245, osf1_246, osg_363, osg_365, \
                         osg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_20 * nsg_365[k]
                   + f_6 * osf0_245[k]
                   - f_7 * osf1_245[k]
                   + f_3 * pc_x[k] * osg_365[k];

        t_510[k] = f_20 * nsg_366[k]
                   + f_4 * osf0_246[k]
                   - f_5 * osf1_246[k]
                   + f_3 * pc_x[k] * osg_366[k];

        t_511[k] = f_11 * nsg_258[k]
                   + f_3 * pc_z[k] * osg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, nsg_275, nsg_369, nsg_370, \
                         nsg_371, osf0_249, osf1_249, osg_365, osg_369, osg_370, \
                         osg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * nsg_275[k]
                   + f_3 * pc_y[k] * osg_365[k];

        t_513[k] = f_20 * nsg_369[k]
                   + f_4 * osf0_249[k]
                   - f_5 * osf1_249[k]
                   + f_3 * pc_x[k] * osg_369[k];

        t_514[k] = f_20 * nsg_370[k]
                   + f_3 * pc_x[k] * osg_370[k];

        t_515[k] = f_20 * nsg_371[k]
                   + f_3 * pc_x[k] * osg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, nsg_280, nsg_372, nsg_373, \
                         nsg_374, osf0_246, osf1_246, osg_370, osg_372, osg_373, \
                         osg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_20 * nsg_372[k]
                   + f_3 * pc_x[k] * osg_372[k];

        t_517[k] = f_20 * nsg_373[k]
                   + f_3 * pc_x[k] * osg_373[k];

        t_518[k] = f_20 * nsg_374[k]
                   + f_3 * pc_x[k] * osg_374[k];

        t_519[k] = f_11 * nsg_280[k]
                   + f_1 * osf0_246[k]
                   - f_2 * osf1_246[k]
                   + f_3 * pc_y[k] * osg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, nsg_265, nsg_282, nsg_283, osf0_248, \
                         osf0_249, osf1_248, osf1_249, osg_370, osg_372, \
                         osg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * nsg_265[k]
                   + f_3 * pc_z[k] * osg_370[k];

        t_521[k] = f_11 * nsg_282[k]
                   + f_6 * osf0_248[k]
                   - f_7 * osf1_248[k]
                   + f_3 * pc_y[k] * osg_372[k];

        t_522[k] = f_11 * nsg_283[k]
                   + f_4 * osf0_249[k]
                   - f_5 * osf1_249[k]
                   + f_3 * pc_y[k] * osg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, nsg_269, nsg_284, nsg_375, \
                         osf0_249, osf0_250, osf1_249, osf1_250, osg_374, \
                         osg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * nsg_284[k]
                   + f_3 * pc_y[k] * osg_374[k];

        t_524[k] = f_11 * nsg_269[k]
                   + f_1 * osf0_249[k]
                   - f_2 * osf1_249[k]
                   + f_3 * pc_z[k] * osg_374[k];

        t_525[k] = f_20 * nsg_375[k]
                   + f_1 * osf0_250[k]
                   - f_2 * osf1_250[k]
                   + f_3 * pc_x[k] * osg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, nsg_270, nsg_285, \
                         nsg_287, nsg_378, osf0_253, osf1_253, osg_375, osg_377, \
                         osg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * nsg_285[k]
                   + f_3 * pc_y[k] * osg_375[k];

        t_527[k] = f_18 * nsg_270[k]
                   + f_3 * pc_z[k] * osg_375[k];

        t_528[k] = f_20 * nsg_378[k]
                   + f_6 * osf0_253[k]
                   - f_7 * osf1_253[k]
                   + f_3 * pc_x[k] * osg_378[k];

        t_529[k] = f_10 * nsg_287[k]
                   + f_3 * pc_y[k] * osg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, nsg_273, nsg_380, nsg_381, osf0_255, \
                         osf0_256, osf1_255, osf1_256, osg_378, osg_380, \
                         osg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_20 * nsg_380[k]
                   + f_6 * osf0_255[k]
                   - f_7 * osf1_255[k]
                   + f_3 * pc_x[k] * osg_380[k];

        t_531[k] = f_20 * nsg_381[k]
                   + f_4 * osf0_256[k]
                   - f_5 * osf1_256[k]
                   + f_3 * pc_x[k] * osg_381[k];

        t_532[k] = f_18 * nsg_273[k]
                   + f_3 * pc_z[k] * osg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, nsg_290, nsg_384, nsg_385, \
                         nsg_386, osf0_259, osf1_259, osg_380, osg_384, osg_385, \
                         osg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * nsg_290[k]
                   + f_3 * pc_y[k] * osg_380[k];

        t_534[k] = f_20 * nsg_384[k]
                   + f_4 * osf0_259[k]
                   - f_5 * osf1_259[k]
                   + f_3 * pc_x[k] * osg_384[k];

        t_535[k] = f_20 * nsg_385[k]
                   + f_3 * pc_x[k] * osg_385[k];

        t_536[k] = f_20 * nsg_386[k]
                   + f_3 * pc_x[k] * osg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, nsg_295, nsg_387, nsg_388, \
                         nsg_389, osf0_256, osf1_256, osg_385, osg_387, osg_388, \
                         osg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_20 * nsg_387[k]
                   + f_3 * pc_x[k] * osg_387[k];

        t_538[k] = f_20 * nsg_388[k]
                   + f_3 * pc_x[k] * osg_388[k];

        t_539[k] = f_20 * nsg_389[k]
                   + f_3 * pc_x[k] * osg_389[k];

        t_540[k] = f_10 * nsg_295[k]
                   + f_1 * osf0_256[k]
                   - f_2 * osf1_256[k]
                   + f_3 * pc_y[k] * osg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, nsg_280, nsg_297, nsg_298, osf0_258, \
                         osf0_259, osf1_258, osf1_259, osg_385, osg_387, \
                         osg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_18 * nsg_280[k]
                   + f_3 * pc_z[k] * osg_385[k];

        t_542[k] = f_10 * nsg_297[k]
                   + f_6 * osf0_258[k]
                   - f_7 * osf1_258[k]
                   + f_3 * pc_y[k] * osg_387[k];

        t_543[k] = f_10 * nsg_298[k]
                   + f_4 * osf0_259[k]
                   - f_5 * osf1_259[k]
                   + f_3 * pc_y[k] * osg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_y, pc_z, nsh0_420, nsg_284, \
                         nsg_299, nsg_300, nsh1_420, osf0_259, osf1_259, osg_389, \
                         osg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * nsg_299[k]
                   + f_3 * pc_y[k] * osg_389[k];

        t_545[k] = f_18 * nsg_284[k]
                   + f_1 * osf0_259[k]
                   - f_2 * osf1_259[k]
                   + f_3 * pc_z[k] * osg_389[k];

        t_546[k] = pa_y[k] * nsh0_420[k]
                   - f_8 * pc_y[k] * nsh1_420[k];

        t_547[k] = f_9 * nsg_300[k]
                   + f_3 * pc_y[k] * osg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pa_y, pc_y, pc_z, nsh0_423, nsh0_425, \
                         nsg_285, nsg_301, nsg_302, nsh1_423, nsh1_425, osg_390, \
                         osg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_20 * nsg_285[k]
                   + f_3 * pc_z[k] * osg_390[k];

        t_549[k] = pa_y[k] * nsh0_423[k]
                   + f_10 * nsg_301[k]
                   - f_8 * pc_y[k] * nsh1_423[k];

        t_550[k] = f_9 * nsg_302[k]
                   + f_3 * pc_y[k] * osg_392[k];

        t_551[k] = pa_y[k] * nsh0_425[k]
                   - f_8 * pc_y[k] * nsh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_y, pc_y, pc_z, nsh0_426, nsh0_429, \
                         nsg_288, nsg_303, nsg_305, nsh1_426, nsh1_429, osg_393, \
                         osg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_y[k] * nsh0_426[k]
                   + f_11 * nsg_303[k]
                   - f_8 * pc_y[k] * nsh1_426[k];

        t_553[k] = f_20 * nsg_288[k]
                   + f_3 * pc_z[k] * osg_393[k];

        t_554[k] = f_9 * nsg_305[k]
                   + f_3 * pc_y[k] * osg_395[k];

        t_555[k] = pa_y[k] * nsh0_429[k]
                   - f_8 * pc_y[k] * nsh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, nsg_400, nsg_401, nsg_402, \
                         nsg_403, nsg_404, osg_400, osg_401, osg_402, osg_403, \
                         osg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_20 * nsg_400[k]
                   + f_3 * pc_x[k] * osg_400[k];

        t_557[k] = f_20 * nsg_401[k]
                   + f_3 * pc_x[k] * osg_401[k];

        t_558[k] = f_20 * nsg_402[k]
                   + f_3 * pc_x[k] * osg_402[k];

        t_559[k] = f_20 * nsg_403[k]
                   + f_3 * pc_x[k] * osg_403[k];

        t_560[k] = f_20 * nsg_404[k]
                   + f_3 * pc_x[k] * osg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, nsg_295, nsg_310, nsg_312, osf0_266, \
                         osf0_268, osf1_266, osf1_268, osg_400, \
                         osg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * nsg_310[k]
                   + f_1 * osf0_266[k]
                   - f_2 * osf1_266[k]
                   + f_3 * pc_y[k] * osg_400[k];

        t_562[k] = f_20 * nsg_295[k]
                   + f_3 * pc_z[k] * osg_400[k];

        t_563[k] = f_9 * nsg_312[k]
                   + f_6 * osf0_268[k]
                   - f_7 * osf1_268[k]
                   + f_3 * pc_y[k] * osg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_y, pc_y, nsh0_440, nsg_313, nsg_314, \
                         nsh1_440, osf0_269, osf1_269, osg_403, \
                         osg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * nsg_313[k]
                   + f_4 * osf0_269[k]
                   - f_5 * osf1_269[k]
                   + f_3 * pc_y[k] * osg_403[k];

        t_565[k] = f_9 * nsg_314[k]
                   + f_3 * pc_y[k] * osg_404[k];

        t_566[k] = pa_y[k] * nsh0_440[k]
                   - f_8 * pc_y[k] * nsh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, nsg_300, \
                         nsg_405, osf0_270, osf1_270, osg_405, osg_406, \
                         osg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_20 * nsg_405[k]
                   + f_1 * osf0_270[k]
                   - f_2 * osf1_270[k]
                   + f_3 * pc_x[k] * osg_405[k];

        t_568[k] = f_3 * pc_y[k] * osg_405[k];

        t_569[k] = f_19 * nsg_300[k]
                   + f_3 * pc_z[k] * osg_405[k];

        t_570[k] = f_4 * osf0_270[k]
                   - f_5 * osf1_270[k]
                   + f_3 * pc_y[k] * osg_406[k];

        t_571[k] = f_3 * pc_y[k] * osg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, nsg_410, osf0_271, osf0_272, \
                         osf0_275, osf1_271, osf1_272, osf1_275, osg_408, osg_409, \
                         osg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_20 * nsg_410[k]
                   + f_6 * osf0_275[k]
                   - f_7 * osf1_275[k]
                   + f_3 * pc_x[k] * osg_410[k];

        t_573[k] = f_6 * osf0_271[k]
                   - f_7 * osf1_271[k]
                   + f_3 * pc_y[k] * osg_408[k];

        t_574[k] = f_4 * osf0_272[k]
                   - f_5 * osf1_272[k]
                   + f_3 * pc_y[k] * osg_409[k];

        t_575[k] = f_3 * pc_y[k] * osg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, nsg_414, nsg_415, nsg_416, nsg_417, \
                         osf0_279, osf1_279, osg_414, osg_415, osg_416, \
                         osg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_20 * nsg_414[k]
                   + f_4 * osf0_279[k]
                   - f_5 * osf1_279[k]
                   + f_3 * pc_x[k] * osg_414[k];

        t_577[k] = f_20 * nsg_415[k]
                   + f_3 * pc_x[k] * osg_415[k];

        t_578[k] = f_20 * nsg_416[k]
                   + f_3 * pc_x[k] * osg_416[k];

        t_579[k] = f_20 * nsg_417[k]
                   + f_3 * pc_x[k] * osg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, nsg_419, osf0_276, osf0_277, \
                         osf1_276, osf1_277, osg_414, osg_415, osg_416, \
                         osg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_3 * pc_y[k] * osg_414[k];

        t_581[k] = f_20 * nsg_419[k]
                   + f_3 * pc_x[k] * osg_419[k];

        t_582[k] = f_1 * osf0_276[k]
                   - f_2 * osf1_276[k]
                   + f_3 * pc_y[k] * osg_415[k];

        t_583[k] = f_13 * osf0_277[k]
                   - f_14 * osf1_277[k]
                   + f_3 * pc_y[k] * osg_416[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, nsg_314, osf0_278, osf0_279, \
                         osf1_278, osf1_279, osg_417, osg_418, \
                         osg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * osf0_278[k]
                   - f_7 * osf1_278[k]
                   + f_3 * pc_y[k] * osg_417[k];

        t_585[k] = f_4 * osf0_279[k]
                   - f_5 * osf1_279[k]
                   + f_3 * pc_y[k] * osg_418[k];

        t_586[k] = f_3 * pc_y[k] * osg_419[k];

        t_587[k] = f_19 * nsg_314[k]
                   + f_1 * osf0_279[k]
                   - f_2 * osf1_279[k]
                   + f_3 * pc_z[k] * osg_419[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, nsg_315, nsg_420, \
                         nsg_423, osf0_280, osf0_283, osf1_280, osf1_283, osg_420, \
                         osg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_18 * nsg_420[k]
                   + f_1 * osf0_280[k]
                   - f_2 * osf1_280[k]
                   + f_3 * pc_x[k] * osg_420[k];

        t_589[k] = f_17 * nsg_315[k]
                   + f_3 * pc_y[k] * osg_420[k];

        t_590[k] = f_3 * pc_z[k] * osg_420[k];

        t_591[k] = f_18 * nsg_423[k]
                   + f_6 * osf0_283[k]
                   - f_7 * osf1_283[k]
                   + f_3 * pc_x[k] * osg_423[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, nsg_426, osf0_280, osf0_286, \
                         osf1_280, osf1_286, osg_421, osg_422, osg_423, \
                         osg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * osg_421[k];

        t_593[k] = f_4 * osf0_280[k]
                   - f_5 * osf1_280[k]
                   + f_3 * pc_z[k] * osg_422[k];

        t_594[k] = f_18 * nsg_426[k]
                   + f_4 * osf0_286[k]
                   - f_5 * osf1_286[k]
                   + f_3 * pc_x[k] * osg_426[k];

        t_595[k] = f_3 * pc_z[k] * osg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, nsg_320, nsg_430, \
                         osf0_282, osf1_282, osg_425, osg_426, \
                         osg_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * nsg_320[k]
                   + f_3 * pc_y[k] * osg_425[k];

        t_597[k] = f_6 * osf0_282[k]
                   - f_7 * osf1_282[k]
                   + f_3 * pc_z[k] * osg_425[k];

        t_598[k] = f_18 * nsg_430[k]
                   + f_3 * pc_x[k] * osg_430[k];

        t_599[k] = f_3 * pc_z[k] * osg_426[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, nsg_325, nsg_432, nsg_433, \
                         nsg_434, osf0_286, osf1_286, osg_430, osg_432, osg_433, \
                         osg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_18 * nsg_432[k]
                   + f_3 * pc_x[k] * osg_432[k];

        t_601[k] = f_18 * nsg_433[k]
                   + f_3 * pc_x[k] * osg_433[k];

        t_602[k] = f_18 * nsg_434[k]
                   + f_3 * pc_x[k] * osg_434[k];

        t_603[k] = f_17 * nsg_325[k]
                   + f_1 * osf0_286[k]
                   - f_2 * osf1_286[k]
                   + f_3 * pc_y[k] * osg_430[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_441 = buffer.data(nsh0 + 441);
    const auto *nsh0_444 = buffer.data(nsh0 + 444);
    const auto *nsh0_447 = buffer.data(nsh0 + 447);
    const auto *nsh0_456 = buffer.data(nsh0 + 456);

    const auto *nsg_315 = buffer.data(nsg + 315);
    const auto *nsg_318 = buffer.data(nsg + 318);
    const auto *nsg_325 = buffer.data(nsg + 325);
    const auto *nsg_329 = buffer.data(nsg + 329);
    const auto *nsg_330 = buffer.data(nsg + 330);
    const auto *nsg_332 = buffer.data(nsg + 332);
    const auto *nsg_333 = buffer.data(nsg + 333);
    const auto *nsg_335 = buffer.data(nsg + 335);
    const auto *nsg_340 = buffer.data(nsg + 340);
    const auto *nsg_342 = buffer.data(nsg + 342);
    const auto *nsg_343 = buffer.data(nsg + 343);
    const auto *nsg_344 = buffer.data(nsg + 344);
    const auto *nsg_345 = buffer.data(nsg + 345);
    const auto *nsg_347 = buffer.data(nsg + 347);
    const auto *nsg_348 = buffer.data(nsg + 348);
    const auto *nsg_350 = buffer.data(nsg + 350);
    const auto *nsg_355 = buffer.data(nsg + 355);
    const auto *nsg_357 = buffer.data(nsg + 357);
    const auto *nsg_358 = buffer.data(nsg + 358);
    const auto *nsg_359 = buffer.data(nsg + 359);
    const auto *nsg_360 = buffer.data(nsg + 360);
    const auto *nsg_362 = buffer.data(nsg + 362);
    const auto *nsg_363 = buffer.data(nsg + 363);
    const auto *nsg_365 = buffer.data(nsg + 365);
    const auto *nsg_370 = buffer.data(nsg + 370);
    const auto *nsg_372 = buffer.data(nsg + 372);
    const auto *nsg_373 = buffer.data(nsg + 373);
    const auto *nsg_374 = buffer.data(nsg + 374);
    const auto *nsg_375 = buffer.data(nsg + 375);
    const auto *nsg_377 = buffer.data(nsg + 377);
    const auto *nsg_378 = buffer.data(nsg + 378);
    const auto *nsg_380 = buffer.data(nsg + 380);
    const auto *nsg_385 = buffer.data(nsg + 385);
    const auto *nsg_387 = buffer.data(nsg + 387);
    const auto *nsg_388 = buffer.data(nsg + 388);
    const auto *nsg_389 = buffer.data(nsg + 389);
    const auto *nsg_390 = buffer.data(nsg + 390);
    const auto *nsg_392 = buffer.data(nsg + 392);
    const auto *nsg_395 = buffer.data(nsg + 395);
    const auto *nsg_400 = buffer.data(nsg + 400);
    const auto *nsg_402 = buffer.data(nsg + 402);
    const auto *nsg_403 = buffer.data(nsg + 403);
    const auto *nsg_440 = buffer.data(nsg + 440);
    const auto *nsg_444 = buffer.data(nsg + 444);
    const auto *nsg_445 = buffer.data(nsg + 445);
    const auto *nsg_446 = buffer.data(nsg + 446);
    const auto *nsg_447 = buffer.data(nsg + 447);
    const auto *nsg_448 = buffer.data(nsg + 448);
    const auto *nsg_449 = buffer.data(nsg + 449);
    const auto *nsg_450 = buffer.data(nsg + 450);
    const auto *nsg_453 = buffer.data(nsg + 453);
    const auto *nsg_455 = buffer.data(nsg + 455);
    const auto *nsg_456 = buffer.data(nsg + 456);
    const auto *nsg_459 = buffer.data(nsg + 459);
    const auto *nsg_460 = buffer.data(nsg + 460);
    const auto *nsg_461 = buffer.data(nsg + 461);
    const auto *nsg_462 = buffer.data(nsg + 462);
    const auto *nsg_463 = buffer.data(nsg + 463);
    const auto *nsg_464 = buffer.data(nsg + 464);
    const auto *nsg_465 = buffer.data(nsg + 465);
    const auto *nsg_468 = buffer.data(nsg + 468);
    const auto *nsg_470 = buffer.data(nsg + 470);
    const auto *nsg_471 = buffer.data(nsg + 471);
    const auto *nsg_474 = buffer.data(nsg + 474);
    const auto *nsg_475 = buffer.data(nsg + 475);
    const auto *nsg_476 = buffer.data(nsg + 476);
    const auto *nsg_477 = buffer.data(nsg + 477);
    const auto *nsg_478 = buffer.data(nsg + 478);
    const auto *nsg_479 = buffer.data(nsg + 479);
    const auto *nsg_480 = buffer.data(nsg + 480);
    const auto *nsg_483 = buffer.data(nsg + 483);
    const auto *nsg_485 = buffer.data(nsg + 485);
    const auto *nsg_486 = buffer.data(nsg + 486);
    const auto *nsg_489 = buffer.data(nsg + 489);
    const auto *nsg_490 = buffer.data(nsg + 490);
    const auto *nsg_491 = buffer.data(nsg + 491);
    const auto *nsg_492 = buffer.data(nsg + 492);
    const auto *nsg_493 = buffer.data(nsg + 493);
    const auto *nsg_494 = buffer.data(nsg + 494);
    const auto *nsg_495 = buffer.data(nsg + 495);
    const auto *nsg_498 = buffer.data(nsg + 498);
    const auto *nsg_500 = buffer.data(nsg + 500);
    const auto *nsg_501 = buffer.data(nsg + 501);
    const auto *nsg_504 = buffer.data(nsg + 504);
    const auto *nsg_505 = buffer.data(nsg + 505);
    const auto *nsg_506 = buffer.data(nsg + 506);
    const auto *nsg_507 = buffer.data(nsg + 507);
    const auto *nsg_508 = buffer.data(nsg + 508);
    const auto *nsg_509 = buffer.data(nsg + 509);

    const auto *nsh1_441 = buffer.data(nsh1 + 441);
    const auto *nsh1_444 = buffer.data(nsh1 + 444);
    const auto *nsh1_447 = buffer.data(nsh1 + 447);
    const auto *nsh1_456 = buffer.data(nsh1 + 456);

    const auto *osf0_286 = buffer.data(osf0 + 286);
    const auto *osf0_287 = buffer.data(osf0 + 287);
    const auto *osf0_289 = buffer.data(osf0 + 289);
    const auto *osf0_295 = buffer.data(osf0 + 295);
    const auto *osf0_298 = buffer.data(osf0 + 298);
    const auto *osf0_299 = buffer.data(osf0 + 299);
    const auto *osf0_300 = buffer.data(osf0 + 300);
    const auto *osf0_303 = buffer.data(osf0 + 303);
    const auto *osf0_305 = buffer.data(osf0 + 305);
    const auto *osf0_306 = buffer.data(osf0 + 306);
    const auto *osf0_308 = buffer.data(osf0 + 308);
    const auto *osf0_309 = buffer.data(osf0 + 309);
    const auto *osf0_310 = buffer.data(osf0 + 310);
    const auto *osf0_313 = buffer.data(osf0 + 313);
    const auto *osf0_315 = buffer.data(osf0 + 315);
    const auto *osf0_316 = buffer.data(osf0 + 316);
    const auto *osf0_318 = buffer.data(osf0 + 318);
    const auto *osf0_319 = buffer.data(osf0 + 319);
    const auto *osf0_320 = buffer.data(osf0 + 320);
    const auto *osf0_323 = buffer.data(osf0 + 323);
    const auto *osf0_325 = buffer.data(osf0 + 325);
    const auto *osf0_326 = buffer.data(osf0 + 326);
    const auto *osf0_328 = buffer.data(osf0 + 328);
    const auto *osf0_329 = buffer.data(osf0 + 329);
    const auto *osf0_330 = buffer.data(osf0 + 330);
    const auto *osf0_333 = buffer.data(osf0 + 333);
    const auto *osf0_335 = buffer.data(osf0 + 335);
    const auto *osf0_336 = buffer.data(osf0 + 336);
    const auto *osf0_338 = buffer.data(osf0 + 338);
    const auto *osf0_339 = buffer.data(osf0 + 339);

    const auto *osf1_286 = buffer.data(osf1 + 286);
    const auto *osf1_287 = buffer.data(osf1 + 287);
    const auto *osf1_289 = buffer.data(osf1 + 289);
    const auto *osf1_295 = buffer.data(osf1 + 295);
    const auto *osf1_298 = buffer.data(osf1 + 298);
    const auto *osf1_299 = buffer.data(osf1 + 299);
    const auto *osf1_300 = buffer.data(osf1 + 300);
    const auto *osf1_303 = buffer.data(osf1 + 303);
    const auto *osf1_305 = buffer.data(osf1 + 305);
    const auto *osf1_306 = buffer.data(osf1 + 306);
    const auto *osf1_308 = buffer.data(osf1 + 308);
    const auto *osf1_309 = buffer.data(osf1 + 309);
    const auto *osf1_310 = buffer.data(osf1 + 310);
    const auto *osf1_313 = buffer.data(osf1 + 313);
    const auto *osf1_315 = buffer.data(osf1 + 315);
    const auto *osf1_316 = buffer.data(osf1 + 316);
    const auto *osf1_318 = buffer.data(osf1 + 318);
    const auto *osf1_319 = buffer.data(osf1 + 319);
    const auto *osf1_320 = buffer.data(osf1 + 320);
    const auto *osf1_323 = buffer.data(osf1 + 323);
    const auto *osf1_325 = buffer.data(osf1 + 325);
    const auto *osf1_326 = buffer.data(osf1 + 326);
    const auto *osf1_328 = buffer.data(osf1 + 328);
    const auto *osf1_329 = buffer.data(osf1 + 329);
    const auto *osf1_330 = buffer.data(osf1 + 330);
    const auto *osf1_333 = buffer.data(osf1 + 333);
    const auto *osf1_335 = buffer.data(osf1 + 335);
    const auto *osf1_336 = buffer.data(osf1 + 336);
    const auto *osf1_338 = buffer.data(osf1 + 338);
    const auto *osf1_339 = buffer.data(osf1 + 339);

    const auto *osg_430 = buffer.data(osg + 430);
    const auto *osg_431 = buffer.data(osg + 431);
    const auto *osg_432 = buffer.data(osg + 432);
    const auto *osg_434 = buffer.data(osg + 434);
    const auto *osg_435 = buffer.data(osg + 435);
    const auto *osg_437 = buffer.data(osg + 437);
    const auto *osg_438 = buffer.data(osg + 438);
    const auto *osg_440 = buffer.data(osg + 440);
    const auto *osg_444 = buffer.data(osg + 444);
    const auto *osg_445 = buffer.data(osg + 445);
    const auto *osg_446 = buffer.data(osg + 446);
    const auto *osg_447 = buffer.data(osg + 447);
    const auto *osg_448 = buffer.data(osg + 448);
    const auto *osg_449 = buffer.data(osg + 449);
    const auto *osg_450 = buffer.data(osg + 450);
    const auto *osg_452 = buffer.data(osg + 452);
    const auto *osg_453 = buffer.data(osg + 453);
    const auto *osg_455 = buffer.data(osg + 455);
    const auto *osg_456 = buffer.data(osg + 456);
    const auto *osg_459 = buffer.data(osg + 459);
    const auto *osg_460 = buffer.data(osg + 460);
    const auto *osg_461 = buffer.data(osg + 461);
    const auto *osg_462 = buffer.data(osg + 462);
    const auto *osg_463 = buffer.data(osg + 463);
    const auto *osg_464 = buffer.data(osg + 464);
    const auto *osg_465 = buffer.data(osg + 465);
    const auto *osg_467 = buffer.data(osg + 467);
    const auto *osg_468 = buffer.data(osg + 468);
    const auto *osg_470 = buffer.data(osg + 470);
    const auto *osg_471 = buffer.data(osg + 471);
    const auto *osg_474 = buffer.data(osg + 474);
    const auto *osg_475 = buffer.data(osg + 475);
    const auto *osg_476 = buffer.data(osg + 476);
    const auto *osg_477 = buffer.data(osg + 477);
    const auto *osg_478 = buffer.data(osg + 478);
    const auto *osg_479 = buffer.data(osg + 479);
    const auto *osg_480 = buffer.data(osg + 480);
    const auto *osg_482 = buffer.data(osg + 482);
    const auto *osg_483 = buffer.data(osg + 483);
    const auto *osg_485 = buffer.data(osg + 485);
    const auto *osg_486 = buffer.data(osg + 486);
    const auto *osg_489 = buffer.data(osg + 489);
    const auto *osg_490 = buffer.data(osg + 490);
    const auto *osg_491 = buffer.data(osg + 491);
    const auto *osg_492 = buffer.data(osg + 492);
    const auto *osg_493 = buffer.data(osg + 493);
    const auto *osg_494 = buffer.data(osg + 494);
    const auto *osg_495 = buffer.data(osg + 495);
    const auto *osg_497 = buffer.data(osg + 497);
    const auto *osg_498 = buffer.data(osg + 498);
    const auto *osg_500 = buffer.data(osg + 500);
    const auto *osg_501 = buffer.data(osg + 501);
    const auto *osg_504 = buffer.data(osg + 504);
    const auto *osg_505 = buffer.data(osg + 505);
    const auto *osg_506 = buffer.data(osg + 506);
    const auto *osg_507 = buffer.data(osg + 507);
    const auto *osg_508 = buffer.data(osg + 508);
    const auto *osg_509 = buffer.data(osg + 509);

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pc_y, pc_z, nsg_329, osf0_286, osf0_287, \
                         osf1_286, osf1_287, osg_430, osg_431, osg_432, \
                         osg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * osg_430[k];

        t_605[k] = f_4 * osf0_286[k]
                   - f_5 * osf1_286[k]
                   + f_3 * pc_z[k] * osg_431[k];

        t_606[k] = f_6 * osf0_287[k]
                   - f_7 * osf1_287[k]
                   + f_3 * pc_z[k] * osg_432[k];

        t_607[k] = f_17 * nsg_329[k]
                   + f_3 * pc_y[k] * osg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_z, pc_y, pc_z, nsh0_441, nsg_315, \
                         nsg_330, nsh1_441, osf0_289, osf1_289, osg_434, \
                         osg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_1 * osf0_289[k]
                   - f_2 * osf1_289[k]
                   + f_3 * pc_z[k] * osg_434[k];

        t_609[k] = pa_z[k] * nsh0_441[k]
                   - f_8 * pc_z[k] * nsh1_441[k];

        t_610[k] = f_19 * nsg_330[k]
                   + f_3 * pc_y[k] * osg_435[k];

        t_611[k] = f_9 * nsg_315[k]
                   + f_3 * pc_z[k] * osg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_z, pc_x, pc_y, pc_z, nsh0_444, nsg_332, \
                         nsg_440, nsh1_444, osf0_295, osf1_295, osg_437, \
                         osg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_z[k] * nsh0_444[k]
                   - f_8 * pc_z[k] * nsh1_444[k];

        t_613[k] = f_19 * nsg_332[k]
                   + f_3 * pc_y[k] * osg_437[k];

        t_614[k] = f_18 * nsg_440[k]
                   + f_6 * osf0_295[k]
                   - f_7 * osf1_295[k]
                   + f_3 * pc_x[k] * osg_440[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pa_z, pc_y, pc_z, nsh0_447, nsg_318, nsg_335, \
                         nsh1_447, osg_438, osg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pa_z[k] * nsh0_447[k]
                   - f_8 * pc_z[k] * nsh1_447[k];

        t_616[k] = f_9 * nsg_318[k]
                   + f_3 * pc_z[k] * osg_438[k];

        t_617[k] = f_19 * nsg_335[k]
                   + f_3 * pc_y[k] * osg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pc_x, nsg_444, nsg_445, nsg_446, nsg_447, \
                         osf0_299, osf1_299, osg_444, osg_445, osg_446, \
                         osg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_18 * nsg_444[k]
                   + f_4 * osf0_299[k]
                   - f_5 * osf1_299[k]
                   + f_3 * pc_x[k] * osg_444[k];

        t_619[k] = f_18 * nsg_445[k]
                   + f_3 * pc_x[k] * osg_445[k];

        t_620[k] = f_18 * nsg_446[k]
                   + f_3 * pc_x[k] * osg_446[k];

        t_621[k] = f_18 * nsg_447[k]
                   + f_3 * pc_x[k] * osg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_z, pc_x, pc_z, nsh0_456, nsg_325, \
                         nsg_448, nsg_449, nsh1_456, osg_445, osg_448, \
                         osg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_18 * nsg_448[k]
                   + f_3 * pc_x[k] * osg_448[k];

        t_623[k] = f_18 * nsg_449[k]
                   + f_3 * pc_x[k] * osg_449[k];

        t_624[k] = pa_z[k] * nsh0_456[k]
                   - f_8 * pc_z[k] * nsh1_456[k];

        t_625[k] = f_9 * nsg_325[k]
                   + f_3 * pc_z[k] * osg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, nsg_342, nsg_343, nsg_344, osf0_298, \
                         osf0_299, osf1_298, osf1_299, osg_447, osg_448, \
                         osg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_19 * nsg_342[k]
                   + f_6 * osf0_298[k]
                   - f_7 * osf1_298[k]
                   + f_3 * pc_y[k] * osg_447[k];

        t_627[k] = f_19 * nsg_343[k]
                   + f_4 * osf0_299[k]
                   - f_5 * osf1_299[k]
                   + f_3 * pc_y[k] * osg_448[k];

        t_628[k] = f_19 * nsg_344[k]
                   + f_3 * pc_y[k] * osg_449[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, pc_z, nsg_329, nsg_345, nsg_450, \
                         osf0_299, osf0_300, osf1_299, osf1_300, osg_449, \
                         osg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_9 * nsg_329[k]
                   + f_1 * osf0_299[k]
                   - f_2 * osf1_299[k]
                   + f_3 * pc_z[k] * osg_449[k];

        t_630[k] = f_18 * nsg_450[k]
                   + f_1 * osf0_300[k]
                   - f_2 * osf1_300[k]
                   + f_3 * pc_x[k] * osg_450[k];

        t_631[k] = f_20 * nsg_345[k]
                   + f_3 * pc_y[k] * osg_450[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pc_x, pc_y, pc_z, nsg_330, nsg_347, nsg_453, \
                         osf0_303, osf1_303, osg_450, osg_452, \
                         osg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_10 * nsg_330[k]
                   + f_3 * pc_z[k] * osg_450[k];

        t_633[k] = f_18 * nsg_453[k]
                   + f_6 * osf0_303[k]
                   - f_7 * osf1_303[k]
                   + f_3 * pc_x[k] * osg_453[k];

        t_634[k] = f_20 * nsg_347[k]
                   + f_3 * pc_y[k] * osg_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, nsg_333, nsg_455, nsg_456, osf0_305, \
                         osf0_306, osf1_305, osf1_306, osg_453, osg_455, \
                         osg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_18 * nsg_455[k]
                   + f_6 * osf0_305[k]
                   - f_7 * osf1_305[k]
                   + f_3 * pc_x[k] * osg_455[k];

        t_636[k] = f_18 * nsg_456[k]
                   + f_4 * osf0_306[k]
                   - f_5 * osf1_306[k]
                   + f_3 * pc_x[k] * osg_456[k];

        t_637[k] = f_10 * nsg_333[k]
                   + f_3 * pc_z[k] * osg_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, nsg_350, nsg_459, nsg_460, \
                         nsg_461, osf0_309, osf1_309, osg_455, osg_459, osg_460, \
                         osg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_20 * nsg_350[k]
                   + f_3 * pc_y[k] * osg_455[k];

        t_639[k] = f_18 * nsg_459[k]
                   + f_4 * osf0_309[k]
                   - f_5 * osf1_309[k]
                   + f_3 * pc_x[k] * osg_459[k];

        t_640[k] = f_18 * nsg_460[k]
                   + f_3 * pc_x[k] * osg_460[k];

        t_641[k] = f_18 * nsg_461[k]
                   + f_3 * pc_x[k] * osg_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, nsg_355, nsg_462, nsg_463, \
                         nsg_464, osf0_306, osf1_306, osg_460, osg_462, osg_463, \
                         osg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_18 * nsg_462[k]
                   + f_3 * pc_x[k] * osg_462[k];

        t_643[k] = f_18 * nsg_463[k]
                   + f_3 * pc_x[k] * osg_463[k];

        t_644[k] = f_18 * nsg_464[k]
                   + f_3 * pc_x[k] * osg_464[k];

        t_645[k] = f_20 * nsg_355[k]
                   + f_1 * osf0_306[k]
                   - f_2 * osf1_306[k]
                   + f_3 * pc_y[k] * osg_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, nsg_340, nsg_357, nsg_358, osf0_308, \
                         osf0_309, osf1_308, osf1_309, osg_460, osg_462, \
                         osg_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * nsg_340[k]
                   + f_3 * pc_z[k] * osg_460[k];

        t_647[k] = f_20 * nsg_357[k]
                   + f_6 * osf0_308[k]
                   - f_7 * osf1_308[k]
                   + f_3 * pc_y[k] * osg_462[k];

        t_648[k] = f_20 * nsg_358[k]
                   + f_4 * osf0_309[k]
                   - f_5 * osf1_309[k]
                   + f_3 * pc_y[k] * osg_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, nsg_344, nsg_359, nsg_465, \
                         osf0_309, osf0_310, osf1_309, osf1_310, osg_464, \
                         osg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_20 * nsg_359[k]
                   + f_3 * pc_y[k] * osg_464[k];

        t_650[k] = f_10 * nsg_344[k]
                   + f_1 * osf0_309[k]
                   - f_2 * osf1_309[k]
                   + f_3 * pc_z[k] * osg_464[k];

        t_651[k] = f_18 * nsg_465[k]
                   + f_1 * osf0_310[k]
                   - f_2 * osf1_310[k]
                   + f_3 * pc_x[k] * osg_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, nsg_345, nsg_360, \
                         nsg_362, nsg_468, osf0_313, osf1_313, osg_465, osg_467, \
                         osg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_18 * nsg_360[k]
                   + f_3 * pc_y[k] * osg_465[k];

        t_653[k] = f_11 * nsg_345[k]
                   + f_3 * pc_z[k] * osg_465[k];

        t_654[k] = f_18 * nsg_468[k]
                   + f_6 * osf0_313[k]
                   - f_7 * osf1_313[k]
                   + f_3 * pc_x[k] * osg_468[k];

        t_655[k] = f_18 * nsg_362[k]
                   + f_3 * pc_y[k] * osg_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, nsg_348, nsg_470, nsg_471, osf0_315, \
                         osf0_316, osf1_315, osf1_316, osg_468, osg_470, \
                         osg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_18 * nsg_470[k]
                   + f_6 * osf0_315[k]
                   - f_7 * osf1_315[k]
                   + f_3 * pc_x[k] * osg_470[k];

        t_657[k] = f_18 * nsg_471[k]
                   + f_4 * osf0_316[k]
                   - f_5 * osf1_316[k]
                   + f_3 * pc_x[k] * osg_471[k];

        t_658[k] = f_11 * nsg_348[k]
                   + f_3 * pc_z[k] * osg_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, nsg_365, nsg_474, nsg_475, \
                         nsg_476, osf0_319, osf1_319, osg_470, osg_474, osg_475, \
                         osg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_18 * nsg_365[k]
                   + f_3 * pc_y[k] * osg_470[k];

        t_660[k] = f_18 * nsg_474[k]
                   + f_4 * osf0_319[k]
                   - f_5 * osf1_319[k]
                   + f_3 * pc_x[k] * osg_474[k];

        t_661[k] = f_18 * nsg_475[k]
                   + f_3 * pc_x[k] * osg_475[k];

        t_662[k] = f_18 * nsg_476[k]
                   + f_3 * pc_x[k] * osg_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, nsg_370, nsg_477, nsg_478, \
                         nsg_479, osf0_316, osf1_316, osg_475, osg_477, osg_478, \
                         osg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_18 * nsg_477[k]
                   + f_3 * pc_x[k] * osg_477[k];

        t_664[k] = f_18 * nsg_478[k]
                   + f_3 * pc_x[k] * osg_478[k];

        t_665[k] = f_18 * nsg_479[k]
                   + f_3 * pc_x[k] * osg_479[k];

        t_666[k] = f_18 * nsg_370[k]
                   + f_1 * osf0_316[k]
                   - f_2 * osf1_316[k]
                   + f_3 * pc_y[k] * osg_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, nsg_355, nsg_372, nsg_373, osf0_318, \
                         osf0_319, osf1_318, osf1_319, osg_475, osg_477, \
                         osg_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * nsg_355[k]
                   + f_3 * pc_z[k] * osg_475[k];

        t_668[k] = f_18 * nsg_372[k]
                   + f_6 * osf0_318[k]
                   - f_7 * osf1_318[k]
                   + f_3 * pc_y[k] * osg_477[k];

        t_669[k] = f_18 * nsg_373[k]
                   + f_4 * osf0_319[k]
                   - f_5 * osf1_319[k]
                   + f_3 * pc_y[k] * osg_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, nsg_359, nsg_374, nsg_480, \
                         osf0_319, osf0_320, osf1_319, osf1_320, osg_479, \
                         osg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_18 * nsg_374[k]
                   + f_3 * pc_y[k] * osg_479[k];

        t_671[k] = f_11 * nsg_359[k]
                   + f_1 * osf0_319[k]
                   - f_2 * osf1_319[k]
                   + f_3 * pc_z[k] * osg_479[k];

        t_672[k] = f_18 * nsg_480[k]
                   + f_1 * osf0_320[k]
                   - f_2 * osf1_320[k]
                   + f_3 * pc_x[k] * osg_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, nsg_360, nsg_375, \
                         nsg_377, nsg_483, osf0_323, osf1_323, osg_480, osg_482, \
                         osg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * nsg_375[k]
                   + f_3 * pc_y[k] * osg_480[k];

        t_674[k] = f_18 * nsg_360[k]
                   + f_3 * pc_z[k] * osg_480[k];

        t_675[k] = f_18 * nsg_483[k]
                   + f_6 * osf0_323[k]
                   - f_7 * osf1_323[k]
                   + f_3 * pc_x[k] * osg_483[k];

        t_676[k] = f_11 * nsg_377[k]
                   + f_3 * pc_y[k] * osg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, nsg_363, nsg_485, nsg_486, osf0_325, \
                         osf0_326, osf1_325, osf1_326, osg_483, osg_485, \
                         osg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_18 * nsg_485[k]
                   + f_6 * osf0_325[k]
                   - f_7 * osf1_325[k]
                   + f_3 * pc_x[k] * osg_485[k];

        t_678[k] = f_18 * nsg_486[k]
                   + f_4 * osf0_326[k]
                   - f_5 * osf1_326[k]
                   + f_3 * pc_x[k] * osg_486[k];

        t_679[k] = f_18 * nsg_363[k]
                   + f_3 * pc_z[k] * osg_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, nsg_380, nsg_489, nsg_490, \
                         nsg_491, osf0_329, osf1_329, osg_485, osg_489, osg_490, \
                         osg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * nsg_380[k]
                   + f_3 * pc_y[k] * osg_485[k];

        t_681[k] = f_18 * nsg_489[k]
                   + f_4 * osf0_329[k]
                   - f_5 * osf1_329[k]
                   + f_3 * pc_x[k] * osg_489[k];

        t_682[k] = f_18 * nsg_490[k]
                   + f_3 * pc_x[k] * osg_490[k];

        t_683[k] = f_18 * nsg_491[k]
                   + f_3 * pc_x[k] * osg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, nsg_385, nsg_492, nsg_493, \
                         nsg_494, osf0_326, osf1_326, osg_490, osg_492, osg_493, \
                         osg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_18 * nsg_492[k]
                   + f_3 * pc_x[k] * osg_492[k];

        t_685[k] = f_18 * nsg_493[k]
                   + f_3 * pc_x[k] * osg_493[k];

        t_686[k] = f_18 * nsg_494[k]
                   + f_3 * pc_x[k] * osg_494[k];

        t_687[k] = f_11 * nsg_385[k]
                   + f_1 * osf0_326[k]
                   - f_2 * osf1_326[k]
                   + f_3 * pc_y[k] * osg_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, nsg_370, nsg_387, nsg_388, osf0_328, \
                         osf0_329, osf1_328, osf1_329, osg_490, osg_492, \
                         osg_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_18 * nsg_370[k]
                   + f_3 * pc_z[k] * osg_490[k];

        t_689[k] = f_11 * nsg_387[k]
                   + f_6 * osf0_328[k]
                   - f_7 * osf1_328[k]
                   + f_3 * pc_y[k] * osg_492[k];

        t_690[k] = f_11 * nsg_388[k]
                   + f_4 * osf0_329[k]
                   - f_5 * osf1_329[k]
                   + f_3 * pc_y[k] * osg_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, nsg_374, nsg_389, nsg_495, \
                         osf0_329, osf0_330, osf1_329, osf1_330, osg_494, \
                         osg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * nsg_389[k]
                   + f_3 * pc_y[k] * osg_494[k];

        t_692[k] = f_18 * nsg_374[k]
                   + f_1 * osf0_329[k]
                   - f_2 * osf1_329[k]
                   + f_3 * pc_z[k] * osg_494[k];

        t_693[k] = f_18 * nsg_495[k]
                   + f_1 * osf0_330[k]
                   - f_2 * osf1_330[k]
                   + f_3 * pc_x[k] * osg_495[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, nsg_375, nsg_390, \
                         nsg_392, nsg_498, osf0_333, osf1_333, osg_495, osg_497, \
                         osg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * nsg_390[k]
                   + f_3 * pc_y[k] * osg_495[k];

        t_695[k] = f_20 * nsg_375[k]
                   + f_3 * pc_z[k] * osg_495[k];

        t_696[k] = f_18 * nsg_498[k]
                   + f_6 * osf0_333[k]
                   - f_7 * osf1_333[k]
                   + f_3 * pc_x[k] * osg_498[k];

        t_697[k] = f_10 * nsg_392[k]
                   + f_3 * pc_y[k] * osg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, nsg_378, nsg_500, nsg_501, osf0_335, \
                         osf0_336, osf1_335, osf1_336, osg_498, osg_500, \
                         osg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_18 * nsg_500[k]
                   + f_6 * osf0_335[k]
                   - f_7 * osf1_335[k]
                   + f_3 * pc_x[k] * osg_500[k];

        t_699[k] = f_18 * nsg_501[k]
                   + f_4 * osf0_336[k]
                   - f_5 * osf1_336[k]
                   + f_3 * pc_x[k] * osg_501[k];

        t_700[k] = f_20 * nsg_378[k]
                   + f_3 * pc_z[k] * osg_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, nsg_395, nsg_504, nsg_505, \
                         nsg_506, osf0_339, osf1_339, osg_500, osg_504, osg_505, \
                         osg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * nsg_395[k]
                   + f_3 * pc_y[k] * osg_500[k];

        t_702[k] = f_18 * nsg_504[k]
                   + f_4 * osf0_339[k]
                   - f_5 * osf1_339[k]
                   + f_3 * pc_x[k] * osg_504[k];

        t_703[k] = f_18 * nsg_505[k]
                   + f_3 * pc_x[k] * osg_505[k];

        t_704[k] = f_18 * nsg_506[k]
                   + f_3 * pc_x[k] * osg_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, nsg_400, nsg_507, nsg_508, \
                         nsg_509, osf0_336, osf1_336, osg_505, osg_507, osg_508, \
                         osg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_18 * nsg_507[k]
                   + f_3 * pc_x[k] * osg_507[k];

        t_706[k] = f_18 * nsg_508[k]
                   + f_3 * pc_x[k] * osg_508[k];

        t_707[k] = f_18 * nsg_509[k]
                   + f_3 * pc_x[k] * osg_509[k];

        t_708[k] = f_10 * nsg_400[k]
                   + f_1 * osf0_336[k]
                   - f_2 * osf1_336[k]
                   + f_3 * pc_y[k] * osg_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, nsg_385, nsg_402, nsg_403, osf0_338, \
                         osf0_339, osf1_338, osf1_339, osg_505, osg_507, \
                         osg_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_20 * nsg_385[k]
                   + f_3 * pc_z[k] * osg_505[k];

        t_710[k] = f_10 * nsg_402[k]
                   + f_6 * osf0_338[k]
                   - f_7 * osf1_338[k]
                   + f_3 * pc_y[k] * osg_507[k];

        t_711[k] = f_10 * nsg_403[k]
                   + f_4 * osf0_339[k]
                   - f_5 * osf1_339[k]
                   + f_3 * pc_y[k] * osg_508[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_567 = buffer.data(nsh0 + 567);
    const auto *nsh0_570 = buffer.data(nsh0 + 570);
    const auto *nsh0_572 = buffer.data(nsh0 + 572);
    const auto *nsh0_573 = buffer.data(nsh0 + 573);
    const auto *nsh0_576 = buffer.data(nsh0 + 576);
    const auto *nsh0_587 = buffer.data(nsh0 + 587);
    const auto *nsh0_588 = buffer.data(nsh0 + 588);
    const auto *nsh0_591 = buffer.data(nsh0 + 591);
    const auto *nsh0_594 = buffer.data(nsh0 + 594);
    const auto *nsh0_603 = buffer.data(nsh0 + 603);

    const auto *nsg_389 = buffer.data(nsg + 389);
    const auto *nsg_390 = buffer.data(nsg + 390);
    const auto *nsg_393 = buffer.data(nsg + 393);
    const auto *nsg_400 = buffer.data(nsg + 400);
    const auto *nsg_404 = buffer.data(nsg + 404);
    const auto *nsg_405 = buffer.data(nsg + 405);
    const auto *nsg_406 = buffer.data(nsg + 406);
    const auto *nsg_407 = buffer.data(nsg + 407);
    const auto *nsg_408 = buffer.data(nsg + 408);
    const auto *nsg_410 = buffer.data(nsg + 410);
    const auto *nsg_415 = buffer.data(nsg + 415);
    const auto *nsg_417 = buffer.data(nsg + 417);
    const auto *nsg_418 = buffer.data(nsg + 418);
    const auto *nsg_419 = buffer.data(nsg + 419);
    const auto *nsg_420 = buffer.data(nsg + 420);
    const auto *nsg_423 = buffer.data(nsg + 423);
    const auto *nsg_425 = buffer.data(nsg + 425);
    const auto *nsg_430 = buffer.data(nsg + 430);
    const auto *nsg_434 = buffer.data(nsg + 434);
    const auto *nsg_435 = buffer.data(nsg + 435);
    const auto *nsg_437 = buffer.data(nsg + 437);
    const auto *nsg_438 = buffer.data(nsg + 438);
    const auto *nsg_440 = buffer.data(nsg + 440);
    const auto *nsg_445 = buffer.data(nsg + 445);
    const auto *nsg_447 = buffer.data(nsg + 447);
    const auto *nsg_448 = buffer.data(nsg + 448);
    const auto *nsg_449 = buffer.data(nsg + 449);
    const auto *nsg_450 = buffer.data(nsg + 450);
    const auto *nsg_452 = buffer.data(nsg + 452);
    const auto *nsg_453 = buffer.data(nsg + 453);
    const auto *nsg_455 = buffer.data(nsg + 455);
    const auto *nsg_460 = buffer.data(nsg + 460);
    const auto *nsg_462 = buffer.data(nsg + 462);
    const auto *nsg_463 = buffer.data(nsg + 463);
    const auto *nsg_464 = buffer.data(nsg + 464);
    const auto *nsg_465 = buffer.data(nsg + 465);
    const auto *nsg_467 = buffer.data(nsg + 467);
    const auto *nsg_470 = buffer.data(nsg + 470);
    const auto *nsg_520 = buffer.data(nsg + 520);
    const auto *nsg_521 = buffer.data(nsg + 521);
    const auto *nsg_522 = buffer.data(nsg + 522);
    const auto *nsg_523 = buffer.data(nsg + 523);
    const auto *nsg_524 = buffer.data(nsg + 524);
    const auto *nsg_525 = buffer.data(nsg + 525);
    const auto *nsg_530 = buffer.data(nsg + 530);
    const auto *nsg_534 = buffer.data(nsg + 534);
    const auto *nsg_535 = buffer.data(nsg + 535);
    const auto *nsg_536 = buffer.data(nsg + 536);
    const auto *nsg_537 = buffer.data(nsg + 537);
    const auto *nsg_539 = buffer.data(nsg + 539);
    const auto *nsg_540 = buffer.data(nsg + 540);
    const auto *nsg_543 = buffer.data(nsg + 543);
    const auto *nsg_546 = buffer.data(nsg + 546);
    const auto *nsg_550 = buffer.data(nsg + 550);
    const auto *nsg_552 = buffer.data(nsg + 552);
    const auto *nsg_553 = buffer.data(nsg + 553);
    const auto *nsg_554 = buffer.data(nsg + 554);
    const auto *nsg_560 = buffer.data(nsg + 560);
    const auto *nsg_564 = buffer.data(nsg + 564);
    const auto *nsg_565 = buffer.data(nsg + 565);
    const auto *nsg_566 = buffer.data(nsg + 566);
    const auto *nsg_567 = buffer.data(nsg + 567);
    const auto *nsg_568 = buffer.data(nsg + 568);
    const auto *nsg_569 = buffer.data(nsg + 569);
    const auto *nsg_570 = buffer.data(nsg + 570);
    const auto *nsg_573 = buffer.data(nsg + 573);
    const auto *nsg_575 = buffer.data(nsg + 575);
    const auto *nsg_576 = buffer.data(nsg + 576);
    const auto *nsg_579 = buffer.data(nsg + 579);
    const auto *nsg_580 = buffer.data(nsg + 580);
    const auto *nsg_581 = buffer.data(nsg + 581);
    const auto *nsg_582 = buffer.data(nsg + 582);
    const auto *nsg_583 = buffer.data(nsg + 583);
    const auto *nsg_584 = buffer.data(nsg + 584);
    const auto *nsg_585 = buffer.data(nsg + 585);
    const auto *nsg_588 = buffer.data(nsg + 588);
    const auto *nsg_590 = buffer.data(nsg + 590);
    const auto *nsg_591 = buffer.data(nsg + 591);
    const auto *nsg_594 = buffer.data(nsg + 594);
    const auto *nsg_595 = buffer.data(nsg + 595);
    const auto *nsg_596 = buffer.data(nsg + 596);

    const auto *nsh1_567 = buffer.data(nsh1 + 567);
    const auto *nsh1_570 = buffer.data(nsh1 + 570);
    const auto *nsh1_572 = buffer.data(nsh1 + 572);
    const auto *nsh1_573 = buffer.data(nsh1 + 573);
    const auto *nsh1_576 = buffer.data(nsh1 + 576);
    const auto *nsh1_587 = buffer.data(nsh1 + 587);
    const auto *nsh1_588 = buffer.data(nsh1 + 588);
    const auto *nsh1_591 = buffer.data(nsh1 + 591);
    const auto *nsh1_594 = buffer.data(nsh1 + 594);
    const auto *nsh1_603 = buffer.data(nsh1 + 603);

    const auto *osf0_339 = buffer.data(osf0 + 339);
    const auto *osf0_346 = buffer.data(osf0 + 346);
    const auto *osf0_348 = buffer.data(osf0 + 348);
    const auto *osf0_349 = buffer.data(osf0 + 349);
    const auto *osf0_350 = buffer.data(osf0 + 350);
    const auto *osf0_351 = buffer.data(osf0 + 351);
    const auto *osf0_352 = buffer.data(osf0 + 352);
    const auto *osf0_355 = buffer.data(osf0 + 355);
    const auto *osf0_356 = buffer.data(osf0 + 356);
    const auto *osf0_357 = buffer.data(osf0 + 357);
    const auto *osf0_358 = buffer.data(osf0 + 358);
    const auto *osf0_359 = buffer.data(osf0 + 359);
    const auto *osf0_360 = buffer.data(osf0 + 360);
    const auto *osf0_362 = buffer.data(osf0 + 362);
    const auto *osf0_363 = buffer.data(osf0 + 363);
    const auto *osf0_366 = buffer.data(osf0 + 366);
    const auto *osf0_367 = buffer.data(osf0 + 367);
    const auto *osf0_369 = buffer.data(osf0 + 369);
    const auto *osf0_375 = buffer.data(osf0 + 375);
    const auto *osf0_378 = buffer.data(osf0 + 378);
    const auto *osf0_379 = buffer.data(osf0 + 379);
    const auto *osf0_380 = buffer.data(osf0 + 380);
    const auto *osf0_383 = buffer.data(osf0 + 383);
    const auto *osf0_385 = buffer.data(osf0 + 385);
    const auto *osf0_386 = buffer.data(osf0 + 386);
    const auto *osf0_388 = buffer.data(osf0 + 388);
    const auto *osf0_389 = buffer.data(osf0 + 389);
    const auto *osf0_390 = buffer.data(osf0 + 390);
    const auto *osf0_393 = buffer.data(osf0 + 393);
    const auto *osf0_395 = buffer.data(osf0 + 395);
    const auto *osf0_396 = buffer.data(osf0 + 396);
    const auto *osf0_399 = buffer.data(osf0 + 399);

    const auto *osf1_339 = buffer.data(osf1 + 339);
    const auto *osf1_346 = buffer.data(osf1 + 346);
    const auto *osf1_348 = buffer.data(osf1 + 348);
    const auto *osf1_349 = buffer.data(osf1 + 349);
    const auto *osf1_350 = buffer.data(osf1 + 350);
    const auto *osf1_351 = buffer.data(osf1 + 351);
    const auto *osf1_352 = buffer.data(osf1 + 352);
    const auto *osf1_355 = buffer.data(osf1 + 355);
    const auto *osf1_356 = buffer.data(osf1 + 356);
    const auto *osf1_357 = buffer.data(osf1 + 357);
    const auto *osf1_358 = buffer.data(osf1 + 358);
    const auto *osf1_359 = buffer.data(osf1 + 359);
    const auto *osf1_360 = buffer.data(osf1 + 360);
    const auto *osf1_362 = buffer.data(osf1 + 362);
    const auto *osf1_363 = buffer.data(osf1 + 363);
    const auto *osf1_366 = buffer.data(osf1 + 366);
    const auto *osf1_367 = buffer.data(osf1 + 367);
    const auto *osf1_369 = buffer.data(osf1 + 369);
    const auto *osf1_375 = buffer.data(osf1 + 375);
    const auto *osf1_378 = buffer.data(osf1 + 378);
    const auto *osf1_379 = buffer.data(osf1 + 379);
    const auto *osf1_380 = buffer.data(osf1 + 380);
    const auto *osf1_383 = buffer.data(osf1 + 383);
    const auto *osf1_385 = buffer.data(osf1 + 385);
    const auto *osf1_386 = buffer.data(osf1 + 386);
    const auto *osf1_388 = buffer.data(osf1 + 388);
    const auto *osf1_389 = buffer.data(osf1 + 389);
    const auto *osf1_390 = buffer.data(osf1 + 390);
    const auto *osf1_393 = buffer.data(osf1 + 393);
    const auto *osf1_395 = buffer.data(osf1 + 395);
    const auto *osf1_396 = buffer.data(osf1 + 396);
    const auto *osf1_399 = buffer.data(osf1 + 399);

    const auto *osg_509 = buffer.data(osg + 509);
    const auto *osg_510 = buffer.data(osg + 510);
    const auto *osg_512 = buffer.data(osg + 512);
    const auto *osg_513 = buffer.data(osg + 513);
    const auto *osg_515 = buffer.data(osg + 515);
    const auto *osg_520 = buffer.data(osg + 520);
    const auto *osg_521 = buffer.data(osg + 521);
    const auto *osg_522 = buffer.data(osg + 522);
    const auto *osg_523 = buffer.data(osg + 523);
    const auto *osg_524 = buffer.data(osg + 524);
    const auto *osg_525 = buffer.data(osg + 525);
    const auto *osg_526 = buffer.data(osg + 526);
    const auto *osg_527 = buffer.data(osg + 527);
    const auto *osg_528 = buffer.data(osg + 528);
    const auto *osg_529 = buffer.data(osg + 529);
    const auto *osg_530 = buffer.data(osg + 530);
    const auto *osg_534 = buffer.data(osg + 534);
    const auto *osg_535 = buffer.data(osg + 535);
    const auto *osg_536 = buffer.data(osg + 536);
    const auto *osg_537 = buffer.data(osg + 537);
    const auto *osg_538 = buffer.data(osg + 538);
    const auto *osg_539 = buffer.data(osg + 539);
    const auto *osg_540 = buffer.data(osg + 540);
    const auto *osg_541 = buffer.data(osg + 541);
    const auto *osg_542 = buffer.data(osg + 542);
    const auto *osg_543 = buffer.data(osg + 543);
    const auto *osg_545 = buffer.data(osg + 545);
    const auto *osg_546 = buffer.data(osg + 546);
    const auto *osg_550 = buffer.data(osg + 550);
    const auto *osg_551 = buffer.data(osg + 551);
    const auto *osg_552 = buffer.data(osg + 552);
    const auto *osg_553 = buffer.data(osg + 553);
    const auto *osg_554 = buffer.data(osg + 554);
    const auto *osg_555 = buffer.data(osg + 555);
    const auto *osg_557 = buffer.data(osg + 557);
    const auto *osg_558 = buffer.data(osg + 558);
    const auto *osg_560 = buffer.data(osg + 560);
    const auto *osg_564 = buffer.data(osg + 564);
    const auto *osg_565 = buffer.data(osg + 565);
    const auto *osg_566 = buffer.data(osg + 566);
    const auto *osg_567 = buffer.data(osg + 567);
    const auto *osg_568 = buffer.data(osg + 568);
    const auto *osg_569 = buffer.data(osg + 569);
    const auto *osg_570 = buffer.data(osg + 570);
    const auto *osg_572 = buffer.data(osg + 572);
    const auto *osg_573 = buffer.data(osg + 573);
    const auto *osg_575 = buffer.data(osg + 575);
    const auto *osg_576 = buffer.data(osg + 576);
    const auto *osg_579 = buffer.data(osg + 579);
    const auto *osg_580 = buffer.data(osg + 580);
    const auto *osg_581 = buffer.data(osg + 581);
    const auto *osg_582 = buffer.data(osg + 582);
    const auto *osg_583 = buffer.data(osg + 583);
    const auto *osg_584 = buffer.data(osg + 584);
    const auto *osg_585 = buffer.data(osg + 585);
    const auto *osg_587 = buffer.data(osg + 587);
    const auto *osg_588 = buffer.data(osg + 588);
    const auto *osg_590 = buffer.data(osg + 590);
    const auto *osg_591 = buffer.data(osg + 591);
    const auto *osg_594 = buffer.data(osg + 594);
    const auto *osg_595 = buffer.data(osg + 595);
    const auto *osg_596 = buffer.data(osg + 596);

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_y, pc_y, pc_z, nsh0_567, nsg_389, \
                         nsg_404, nsg_405, nsh1_567, osf0_339, osf1_339, osg_509, \
                         osg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * nsg_404[k]
                   + f_3 * pc_y[k] * osg_509[k];

        t_713[k] = f_20 * nsg_389[k]
                   + f_1 * osf0_339[k]
                   - f_2 * osf1_339[k]
                   + f_3 * pc_z[k] * osg_509[k];

        t_714[k] = pa_y[k] * nsh0_567[k]
                   - f_8 * pc_y[k] * nsh1_567[k];

        t_715[k] = f_9 * nsg_405[k]
                   + f_3 * pc_y[k] * osg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_y, pc_y, pc_z, nsh0_570, nsh0_572, \
                         nsg_390, nsg_406, nsg_407, nsh1_570, nsh1_572, osg_510, \
                         osg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_19 * nsg_390[k]
                   + f_3 * pc_z[k] * osg_510[k];

        t_717[k] = pa_y[k] * nsh0_570[k]
                   + f_10 * nsg_406[k]
                   - f_8 * pc_y[k] * nsh1_570[k];

        t_718[k] = f_9 * nsg_407[k]
                   + f_3 * pc_y[k] * osg_512[k];

        t_719[k] = pa_y[k] * nsh0_572[k]
                   - f_8 * pc_y[k] * nsh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_y, pc_y, pc_z, nsh0_573, nsh0_576, \
                         nsg_393, nsg_408, nsg_410, nsh1_573, nsh1_576, osg_513, \
                         osg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_y[k] * nsh0_573[k]
                   + f_11 * nsg_408[k]
                   - f_8 * pc_y[k] * nsh1_573[k];

        t_721[k] = f_19 * nsg_393[k]
                   + f_3 * pc_z[k] * osg_513[k];

        t_722[k] = f_9 * nsg_410[k]
                   + f_3 * pc_y[k] * osg_515[k];

        t_723[k] = pa_y[k] * nsh0_576[k]
                   - f_8 * pc_y[k] * nsh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, nsg_520, nsg_521, nsg_522, \
                         nsg_523, nsg_524, osg_520, osg_521, osg_522, osg_523, \
                         osg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_18 * nsg_520[k]
                   + f_3 * pc_x[k] * osg_520[k];

        t_725[k] = f_18 * nsg_521[k]
                   + f_3 * pc_x[k] * osg_521[k];

        t_726[k] = f_18 * nsg_522[k]
                   + f_3 * pc_x[k] * osg_522[k];

        t_727[k] = f_18 * nsg_523[k]
                   + f_3 * pc_x[k] * osg_523[k];

        t_728[k] = f_18 * nsg_524[k]
                   + f_3 * pc_x[k] * osg_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, nsg_400, nsg_415, nsg_417, osf0_346, \
                         osf0_348, osf1_346, osf1_348, osg_520, \
                         osg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * nsg_415[k]
                   + f_1 * osf0_346[k]
                   - f_2 * osf1_346[k]
                   + f_3 * pc_y[k] * osg_520[k];

        t_730[k] = f_19 * nsg_400[k]
                   + f_3 * pc_z[k] * osg_520[k];

        t_731[k] = f_9 * nsg_417[k]
                   + f_6 * osf0_348[k]
                   - f_7 * osf1_348[k]
                   + f_3 * pc_y[k] * osg_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pc_y, nsh0_587, nsg_418, nsg_419, \
                         nsh1_587, osf0_349, osf1_349, osg_523, \
                         osg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * nsg_418[k]
                   + f_4 * osf0_349[k]
                   - f_5 * osf1_349[k]
                   + f_3 * pc_y[k] * osg_523[k];

        t_733[k] = f_9 * nsg_419[k]
                   + f_3 * pc_y[k] * osg_524[k];

        t_734[k] = pa_y[k] * nsh0_587[k]
                   - f_8 * pc_y[k] * nsh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, nsg_405, \
                         nsg_525, osf0_350, osf1_350, osg_525, osg_526, \
                         osg_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_18 * nsg_525[k]
                   + f_1 * osf0_350[k]
                   - f_2 * osf1_350[k]
                   + f_3 * pc_x[k] * osg_525[k];

        t_736[k] = f_3 * pc_y[k] * osg_525[k];

        t_737[k] = f_17 * nsg_405[k]
                   + f_3 * pc_z[k] * osg_525[k];

        t_738[k] = f_4 * osf0_350[k]
                   - f_5 * osf1_350[k]
                   + f_3 * pc_y[k] * osg_526[k];

        t_739[k] = f_3 * pc_y[k] * osg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pc_y, nsg_530, osf0_351, osf0_352, \
                         osf0_355, osf1_351, osf1_352, osf1_355, osg_528, osg_529, \
                         osg_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_18 * nsg_530[k]
                   + f_6 * osf0_355[k]
                   - f_7 * osf1_355[k]
                   + f_3 * pc_x[k] * osg_530[k];

        t_741[k] = f_6 * osf0_351[k]
                   - f_7 * osf1_351[k]
                   + f_3 * pc_y[k] * osg_528[k];

        t_742[k] = f_4 * osf0_352[k]
                   - f_5 * osf1_352[k]
                   + f_3 * pc_y[k] * osg_529[k];

        t_743[k] = f_3 * pc_y[k] * osg_530[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pc_x, nsg_534, nsg_535, nsg_536, nsg_537, \
                         osf0_359, osf1_359, osg_534, osg_535, osg_536, \
                         osg_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_18 * nsg_534[k]
                   + f_4 * osf0_359[k]
                   - f_5 * osf1_359[k]
                   + f_3 * pc_x[k] * osg_534[k];

        t_745[k] = f_18 * nsg_535[k]
                   + f_3 * pc_x[k] * osg_535[k];

        t_746[k] = f_18 * nsg_536[k]
                   + f_3 * pc_x[k] * osg_536[k];

        t_747[k] = f_18 * nsg_537[k]
                   + f_3 * pc_x[k] * osg_537[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_x, pc_y, nsg_539, osf0_356, osf0_357, \
                         osf1_356, osf1_357, osg_534, osg_535, osg_536, \
                         osg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_3 * pc_y[k] * osg_534[k];

        t_749[k] = f_18 * nsg_539[k]
                   + f_3 * pc_x[k] * osg_539[k];

        t_750[k] = f_1 * osf0_356[k]
                   - f_2 * osf1_356[k]
                   + f_3 * pc_y[k] * osg_535[k];

        t_751[k] = f_13 * osf0_357[k]
                   - f_14 * osf1_357[k]
                   + f_3 * pc_y[k] * osg_536[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, nsg_419, osf0_358, osf0_359, \
                         osf1_358, osf1_359, osg_537, osg_538, \
                         osg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_6 * osf0_358[k]
                   - f_7 * osf1_358[k]
                   + f_3 * pc_y[k] * osg_537[k];

        t_753[k] = f_4 * osf0_359[k]
                   - f_5 * osf1_359[k]
                   + f_3 * pc_y[k] * osg_538[k];

        t_754[k] = f_3 * pc_y[k] * osg_539[k];

        t_755[k] = f_17 * nsg_419[k]
                   + f_1 * osf0_359[k]
                   - f_2 * osf1_359[k]
                   + f_3 * pc_z[k] * osg_539[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, nsg_420, nsg_540, \
                         nsg_543, osf0_360, osf0_363, osf1_360, osf1_363, osg_540, \
                         osg_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_11 * nsg_540[k]
                   + f_1 * osf0_360[k]
                   - f_2 * osf1_360[k]
                   + f_3 * pc_x[k] * osg_540[k];

        t_757[k] = f_16 * nsg_420[k]
                   + f_3 * pc_y[k] * osg_540[k];

        t_758[k] = f_3 * pc_z[k] * osg_540[k];

        t_759[k] = f_11 * nsg_543[k]
                   + f_6 * osf0_363[k]
                   - f_7 * osf1_363[k]
                   + f_3 * pc_x[k] * osg_543[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, pc_x, pc_z, nsg_546, osf0_360, osf0_366, \
                         osf1_360, osf1_366, osg_541, osg_542, osg_543, \
                         osg_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_3 * pc_z[k] * osg_541[k];

        t_761[k] = f_4 * osf0_360[k]
                   - f_5 * osf1_360[k]
                   + f_3 * pc_z[k] * osg_542[k];

        t_762[k] = f_11 * nsg_546[k]
                   + f_4 * osf0_366[k]
                   - f_5 * osf1_366[k]
                   + f_3 * pc_x[k] * osg_546[k];

        t_763[k] = f_3 * pc_z[k] * osg_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, nsg_425, nsg_550, \
                         osf0_362, osf1_362, osg_545, osg_546, \
                         osg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_16 * nsg_425[k]
                   + f_3 * pc_y[k] * osg_545[k];

        t_765[k] = f_6 * osf0_362[k]
                   - f_7 * osf1_362[k]
                   + f_3 * pc_z[k] * osg_545[k];

        t_766[k] = f_11 * nsg_550[k]
                   + f_3 * pc_x[k] * osg_550[k];

        t_767[k] = f_3 * pc_z[k] * osg_546[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, nsg_430, nsg_552, nsg_553, \
                         nsg_554, osf0_366, osf1_366, osg_550, osg_552, osg_553, \
                         osg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_11 * nsg_552[k]
                   + f_3 * pc_x[k] * osg_552[k];

        t_769[k] = f_11 * nsg_553[k]
                   + f_3 * pc_x[k] * osg_553[k];

        t_770[k] = f_11 * nsg_554[k]
                   + f_3 * pc_x[k] * osg_554[k];

        t_771[k] = f_16 * nsg_430[k]
                   + f_1 * osf0_366[k]
                   - f_2 * osf1_366[k]
                   + f_3 * pc_y[k] * osg_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, pc_y, pc_z, nsg_434, osf0_366, osf0_367, \
                         osf1_366, osf1_367, osg_550, osg_551, osg_552, \
                         osg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * osg_550[k];

        t_773[k] = f_4 * osf0_366[k]
                   - f_5 * osf1_366[k]
                   + f_3 * pc_z[k] * osg_551[k];

        t_774[k] = f_6 * osf0_367[k]
                   - f_7 * osf1_367[k]
                   + f_3 * pc_z[k] * osg_552[k];

        t_775[k] = f_16 * nsg_434[k]
                   + f_3 * pc_y[k] * osg_554[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pa_z, pc_y, pc_z, nsh0_588, nsg_420, \
                         nsg_435, nsh1_588, osf0_369, osf1_369, osg_554, \
                         osg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_1 * osf0_369[k]
                   - f_2 * osf1_369[k]
                   + f_3 * pc_z[k] * osg_554[k];

        t_777[k] = pa_z[k] * nsh0_588[k]
                   - f_8 * pc_z[k] * nsh1_588[k];

        t_778[k] = f_17 * nsg_435[k]
                   + f_3 * pc_y[k] * osg_555[k];

        t_779[k] = f_9 * nsg_420[k]
                   + f_3 * pc_z[k] * osg_555[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pa_z, pc_x, pc_y, pc_z, nsh0_591, nsg_437, \
                         nsg_560, nsh1_591, osf0_375, osf1_375, osg_557, \
                         osg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pa_z[k] * nsh0_591[k]
                   - f_8 * pc_z[k] * nsh1_591[k];

        t_781[k] = f_17 * nsg_437[k]
                   + f_3 * pc_y[k] * osg_557[k];

        t_782[k] = f_11 * nsg_560[k]
                   + f_6 * osf0_375[k]
                   - f_7 * osf1_375[k]
                   + f_3 * pc_x[k] * osg_560[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pa_z, pc_y, pc_z, nsh0_594, nsg_423, nsg_440, \
                         nsh1_594, osg_558, osg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pa_z[k] * nsh0_594[k]
                   - f_8 * pc_z[k] * nsh1_594[k];

        t_784[k] = f_9 * nsg_423[k]
                   + f_3 * pc_z[k] * osg_558[k];

        t_785[k] = f_17 * nsg_440[k]
                   + f_3 * pc_y[k] * osg_560[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pc_x, nsg_564, nsg_565, nsg_566, nsg_567, \
                         osf0_379, osf1_379, osg_564, osg_565, osg_566, \
                         osg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_11 * nsg_564[k]
                   + f_4 * osf0_379[k]
                   - f_5 * osf1_379[k]
                   + f_3 * pc_x[k] * osg_564[k];

        t_787[k] = f_11 * nsg_565[k]
                   + f_3 * pc_x[k] * osg_565[k];

        t_788[k] = f_11 * nsg_566[k]
                   + f_3 * pc_x[k] * osg_566[k];

        t_789[k] = f_11 * nsg_567[k]
                   + f_3 * pc_x[k] * osg_567[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_x, pc_z, nsh0_603, nsg_430, \
                         nsg_568, nsg_569, nsh1_603, osg_565, osg_568, \
                         osg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_11 * nsg_568[k]
                   + f_3 * pc_x[k] * osg_568[k];

        t_791[k] = f_11 * nsg_569[k]
                   + f_3 * pc_x[k] * osg_569[k];

        t_792[k] = pa_z[k] * nsh0_603[k]
                   - f_8 * pc_z[k] * nsh1_603[k];

        t_793[k] = f_9 * nsg_430[k]
                   + f_3 * pc_z[k] * osg_565[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pc_y, nsg_447, nsg_448, nsg_449, osf0_378, \
                         osf0_379, osf1_378, osf1_379, osg_567, osg_568, \
                         osg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_17 * nsg_447[k]
                   + f_6 * osf0_378[k]
                   - f_7 * osf1_378[k]
                   + f_3 * pc_y[k] * osg_567[k];

        t_795[k] = f_17 * nsg_448[k]
                   + f_4 * osf0_379[k]
                   - f_5 * osf1_379[k]
                   + f_3 * pc_y[k] * osg_568[k];

        t_796[k] = f_17 * nsg_449[k]
                   + f_3 * pc_y[k] * osg_569[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pc_x, pc_y, pc_z, nsg_434, nsg_450, nsg_570, \
                         osf0_379, osf0_380, osf1_379, osf1_380, osg_569, \
                         osg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_9 * nsg_434[k]
                   + f_1 * osf0_379[k]
                   - f_2 * osf1_379[k]
                   + f_3 * pc_z[k] * osg_569[k];

        t_798[k] = f_11 * nsg_570[k]
                   + f_1 * osf0_380[k]
                   - f_2 * osf1_380[k]
                   + f_3 * pc_x[k] * osg_570[k];

        t_799[k] = f_19 * nsg_450[k]
                   + f_3 * pc_y[k] * osg_570[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, nsg_435, nsg_452, nsg_573, \
                         osf0_383, osf1_383, osg_570, osg_572, \
                         osg_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_10 * nsg_435[k]
                   + f_3 * pc_z[k] * osg_570[k];

        t_801[k] = f_11 * nsg_573[k]
                   + f_6 * osf0_383[k]
                   - f_7 * osf1_383[k]
                   + f_3 * pc_x[k] * osg_573[k];

        t_802[k] = f_19 * nsg_452[k]
                   + f_3 * pc_y[k] * osg_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, pc_z, nsg_438, nsg_575, nsg_576, osf0_385, \
                         osf0_386, osf1_385, osf1_386, osg_573, osg_575, \
                         osg_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_11 * nsg_575[k]
                   + f_6 * osf0_385[k]
                   - f_7 * osf1_385[k]
                   + f_3 * pc_x[k] * osg_575[k];

        t_804[k] = f_11 * nsg_576[k]
                   + f_4 * osf0_386[k]
                   - f_5 * osf1_386[k]
                   + f_3 * pc_x[k] * osg_576[k];

        t_805[k] = f_10 * nsg_438[k]
                   + f_3 * pc_z[k] * osg_573[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, pc_x, pc_y, nsg_455, nsg_579, nsg_580, \
                         nsg_581, osf0_389, osf1_389, osg_575, osg_579, osg_580, \
                         osg_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_19 * nsg_455[k]
                   + f_3 * pc_y[k] * osg_575[k];

        t_807[k] = f_11 * nsg_579[k]
                   + f_4 * osf0_389[k]
                   - f_5 * osf1_389[k]
                   + f_3 * pc_x[k] * osg_579[k];

        t_808[k] = f_11 * nsg_580[k]
                   + f_3 * pc_x[k] * osg_580[k];

        t_809[k] = f_11 * nsg_581[k]
                   + f_3 * pc_x[k] * osg_581[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pc_x, pc_y, nsg_460, nsg_582, nsg_583, \
                         nsg_584, osf0_386, osf1_386, osg_580, osg_582, osg_583, \
                         osg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_11 * nsg_582[k]
                   + f_3 * pc_x[k] * osg_582[k];

        t_811[k] = f_11 * nsg_583[k]
                   + f_3 * pc_x[k] * osg_583[k];

        t_812[k] = f_11 * nsg_584[k]
                   + f_3 * pc_x[k] * osg_584[k];

        t_813[k] = f_19 * nsg_460[k]
                   + f_1 * osf0_386[k]
                   - f_2 * osf1_386[k]
                   + f_3 * pc_y[k] * osg_580[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pc_y, pc_z, nsg_445, nsg_462, nsg_463, osf0_388, \
                         osf0_389, osf1_388, osf1_389, osg_580, osg_582, \
                         osg_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_10 * nsg_445[k]
                   + f_3 * pc_z[k] * osg_580[k];

        t_815[k] = f_19 * nsg_462[k]
                   + f_6 * osf0_388[k]
                   - f_7 * osf1_388[k]
                   + f_3 * pc_y[k] * osg_582[k];

        t_816[k] = f_19 * nsg_463[k]
                   + f_4 * osf0_389[k]
                   - f_5 * osf1_389[k]
                   + f_3 * pc_y[k] * osg_583[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pc_x, pc_y, pc_z, nsg_449, nsg_464, nsg_585, \
                         osf0_389, osf0_390, osf1_389, osf1_390, osg_584, \
                         osg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_19 * nsg_464[k]
                   + f_3 * pc_y[k] * osg_584[k];

        t_818[k] = f_10 * nsg_449[k]
                   + f_1 * osf0_389[k]
                   - f_2 * osf1_389[k]
                   + f_3 * pc_z[k] * osg_584[k];

        t_819[k] = f_11 * nsg_585[k]
                   + f_1 * osf0_390[k]
                   - f_2 * osf1_390[k]
                   + f_3 * pc_x[k] * osg_585[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_x, pc_y, pc_z, nsg_450, nsg_465, \
                         nsg_467, nsg_588, osf0_393, osf1_393, osg_585, osg_587, \
                         osg_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_20 * nsg_465[k]
                   + f_3 * pc_y[k] * osg_585[k];

        t_821[k] = f_11 * nsg_450[k]
                   + f_3 * pc_z[k] * osg_585[k];

        t_822[k] = f_11 * nsg_588[k]
                   + f_6 * osf0_393[k]
                   - f_7 * osf1_393[k]
                   + f_3 * pc_x[k] * osg_588[k];

        t_823[k] = f_20 * nsg_467[k]
                   + f_3 * pc_y[k] * osg_587[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, nsg_453, nsg_590, nsg_591, osf0_395, \
                         osf0_396, osf1_395, osf1_396, osg_588, osg_590, \
                         osg_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_11 * nsg_590[k]
                   + f_6 * osf0_395[k]
                   - f_7 * osf1_395[k]
                   + f_3 * pc_x[k] * osg_590[k];

        t_825[k] = f_11 * nsg_591[k]
                   + f_4 * osf0_396[k]
                   - f_5 * osf1_396[k]
                   + f_3 * pc_x[k] * osg_591[k];

        t_826[k] = f_11 * nsg_453[k]
                   + f_3 * pc_z[k] * osg_588[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pc_x, pc_y, nsg_470, nsg_594, nsg_595, \
                         nsg_596, osf0_399, osf1_399, osg_590, osg_594, osg_595, \
                         osg_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_20 * nsg_470[k]
                   + f_3 * pc_y[k] * osg_590[k];

        t_828[k] = f_11 * nsg_594[k]
                   + f_4 * osf0_399[k]
                   - f_5 * osf1_399[k]
                   + f_3 * pc_x[k] * osg_594[k];

        t_829[k] = f_11 * nsg_595[k]
                   + f_3 * pc_x[k] * osg_595[k];

        t_830[k] = f_11 * nsg_596[k]
                   + f_3 * pc_x[k] * osg_596[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_735 = buffer.data(nsh0 + 735);
    const auto *nsh0_738 = buffer.data(nsh0 + 738);
    const auto *nsh0_740 = buffer.data(nsh0 + 740);
    const auto *nsh0_741 = buffer.data(nsh0 + 741);
    const auto *nsh0_744 = buffer.data(nsh0 + 744);
    const auto *nsh0_755 = buffer.data(nsh0 + 755);

    const auto *nsg_460 = buffer.data(nsg + 460);
    const auto *nsg_464 = buffer.data(nsg + 464);
    const auto *nsg_465 = buffer.data(nsg + 465);
    const auto *nsg_468 = buffer.data(nsg + 468);
    const auto *nsg_475 = buffer.data(nsg + 475);
    const auto *nsg_477 = buffer.data(nsg + 477);
    const auto *nsg_478 = buffer.data(nsg + 478);
    const auto *nsg_479 = buffer.data(nsg + 479);
    const auto *nsg_480 = buffer.data(nsg + 480);
    const auto *nsg_482 = buffer.data(nsg + 482);
    const auto *nsg_483 = buffer.data(nsg + 483);
    const auto *nsg_485 = buffer.data(nsg + 485);
    const auto *nsg_490 = buffer.data(nsg + 490);
    const auto *nsg_492 = buffer.data(nsg + 492);
    const auto *nsg_493 = buffer.data(nsg + 493);
    const auto *nsg_494 = buffer.data(nsg + 494);
    const auto *nsg_495 = buffer.data(nsg + 495);
    const auto *nsg_497 = buffer.data(nsg + 497);
    const auto *nsg_498 = buffer.data(nsg + 498);
    const auto *nsg_500 = buffer.data(nsg + 500);
    const auto *nsg_505 = buffer.data(nsg + 505);
    const auto *nsg_507 = buffer.data(nsg + 507);
    const auto *nsg_508 = buffer.data(nsg + 508);
    const auto *nsg_509 = buffer.data(nsg + 509);
    const auto *nsg_510 = buffer.data(nsg + 510);
    const auto *nsg_512 = buffer.data(nsg + 512);
    const auto *nsg_513 = buffer.data(nsg + 513);
    const auto *nsg_515 = buffer.data(nsg + 515);
    const auto *nsg_520 = buffer.data(nsg + 520);
    const auto *nsg_522 = buffer.data(nsg + 522);
    const auto *nsg_523 = buffer.data(nsg + 523);
    const auto *nsg_524 = buffer.data(nsg + 524);
    const auto *nsg_525 = buffer.data(nsg + 525);
    const auto *nsg_526 = buffer.data(nsg + 526);
    const auto *nsg_527 = buffer.data(nsg + 527);
    const auto *nsg_528 = buffer.data(nsg + 528);
    const auto *nsg_530 = buffer.data(nsg + 530);
    const auto *nsg_535 = buffer.data(nsg + 535);
    const auto *nsg_537 = buffer.data(nsg + 537);
    const auto *nsg_538 = buffer.data(nsg + 538);
    const auto *nsg_539 = buffer.data(nsg + 539);
    const auto *nsg_597 = buffer.data(nsg + 597);
    const auto *nsg_598 = buffer.data(nsg + 598);
    const auto *nsg_599 = buffer.data(nsg + 599);
    const auto *nsg_600 = buffer.data(nsg + 600);
    const auto *nsg_603 = buffer.data(nsg + 603);
    const auto *nsg_605 = buffer.data(nsg + 605);
    const auto *nsg_606 = buffer.data(nsg + 606);
    const auto *nsg_609 = buffer.data(nsg + 609);
    const auto *nsg_610 = buffer.data(nsg + 610);
    const auto *nsg_611 = buffer.data(nsg + 611);
    const auto *nsg_612 = buffer.data(nsg + 612);
    const auto *nsg_613 = buffer.data(nsg + 613);
    const auto *nsg_614 = buffer.data(nsg + 614);
    const auto *nsg_615 = buffer.data(nsg + 615);
    const auto *nsg_618 = buffer.data(nsg + 618);
    const auto *nsg_620 = buffer.data(nsg + 620);
    const auto *nsg_621 = buffer.data(nsg + 621);
    const auto *nsg_624 = buffer.data(nsg + 624);
    const auto *nsg_625 = buffer.data(nsg + 625);
    const auto *nsg_626 = buffer.data(nsg + 626);
    const auto *nsg_627 = buffer.data(nsg + 627);
    const auto *nsg_628 = buffer.data(nsg + 628);
    const auto *nsg_629 = buffer.data(nsg + 629);
    const auto *nsg_630 = buffer.data(nsg + 630);
    const auto *nsg_633 = buffer.data(nsg + 633);
    const auto *nsg_635 = buffer.data(nsg + 635);
    const auto *nsg_636 = buffer.data(nsg + 636);
    const auto *nsg_639 = buffer.data(nsg + 639);
    const auto *nsg_640 = buffer.data(nsg + 640);
    const auto *nsg_641 = buffer.data(nsg + 641);
    const auto *nsg_642 = buffer.data(nsg + 642);
    const auto *nsg_643 = buffer.data(nsg + 643);
    const auto *nsg_644 = buffer.data(nsg + 644);
    const auto *nsg_655 = buffer.data(nsg + 655);
    const auto *nsg_656 = buffer.data(nsg + 656);
    const auto *nsg_657 = buffer.data(nsg + 657);
    const auto *nsg_658 = buffer.data(nsg + 658);
    const auto *nsg_659 = buffer.data(nsg + 659);
    const auto *nsg_660 = buffer.data(nsg + 660);
    const auto *nsg_665 = buffer.data(nsg + 665);
    const auto *nsg_669 = buffer.data(nsg + 669);
    const auto *nsg_670 = buffer.data(nsg + 670);
    const auto *nsg_671 = buffer.data(nsg + 671);
    const auto *nsg_672 = buffer.data(nsg + 672);
    const auto *nsg_674 = buffer.data(nsg + 674);

    const auto *nsh1_735 = buffer.data(nsh1 + 735);
    const auto *nsh1_738 = buffer.data(nsh1 + 738);
    const auto *nsh1_740 = buffer.data(nsh1 + 740);
    const auto *nsh1_741 = buffer.data(nsh1 + 741);
    const auto *nsh1_744 = buffer.data(nsh1 + 744);
    const auto *nsh1_755 = buffer.data(nsh1 + 755);

    const auto *osf0_396 = buffer.data(osf0 + 396);
    const auto *osf0_398 = buffer.data(osf0 + 398);
    const auto *osf0_399 = buffer.data(osf0 + 399);
    const auto *osf0_400 = buffer.data(osf0 + 400);
    const auto *osf0_403 = buffer.data(osf0 + 403);
    const auto *osf0_405 = buffer.data(osf0 + 405);
    const auto *osf0_406 = buffer.data(osf0 + 406);
    const auto *osf0_408 = buffer.data(osf0 + 408);
    const auto *osf0_409 = buffer.data(osf0 + 409);
    const auto *osf0_410 = buffer.data(osf0 + 410);
    const auto *osf0_413 = buffer.data(osf0 + 413);
    const auto *osf0_415 = buffer.data(osf0 + 415);
    const auto *osf0_416 = buffer.data(osf0 + 416);
    const auto *osf0_418 = buffer.data(osf0 + 418);
    const auto *osf0_419 = buffer.data(osf0 + 419);
    const auto *osf0_420 = buffer.data(osf0 + 420);
    const auto *osf0_423 = buffer.data(osf0 + 423);
    const auto *osf0_425 = buffer.data(osf0 + 425);
    const auto *osf0_426 = buffer.data(osf0 + 426);
    const auto *osf0_428 = buffer.data(osf0 + 428);
    const auto *osf0_429 = buffer.data(osf0 + 429);
    const auto *osf0_436 = buffer.data(osf0 + 436);
    const auto *osf0_438 = buffer.data(osf0 + 438);
    const auto *osf0_439 = buffer.data(osf0 + 439);
    const auto *osf0_440 = buffer.data(osf0 + 440);
    const auto *osf0_441 = buffer.data(osf0 + 441);
    const auto *osf0_442 = buffer.data(osf0 + 442);
    const auto *osf0_445 = buffer.data(osf0 + 445);
    const auto *osf0_446 = buffer.data(osf0 + 446);
    const auto *osf0_447 = buffer.data(osf0 + 447);
    const auto *osf0_448 = buffer.data(osf0 + 448);
    const auto *osf0_449 = buffer.data(osf0 + 449);

    const auto *osf1_396 = buffer.data(osf1 + 396);
    const auto *osf1_398 = buffer.data(osf1 + 398);
    const auto *osf1_399 = buffer.data(osf1 + 399);
    const auto *osf1_400 = buffer.data(osf1 + 400);
    const auto *osf1_403 = buffer.data(osf1 + 403);
    const auto *osf1_405 = buffer.data(osf1 + 405);
    const auto *osf1_406 = buffer.data(osf1 + 406);
    const auto *osf1_408 = buffer.data(osf1 + 408);
    const auto *osf1_409 = buffer.data(osf1 + 409);
    const auto *osf1_410 = buffer.data(osf1 + 410);
    const auto *osf1_413 = buffer.data(osf1 + 413);
    const auto *osf1_415 = buffer.data(osf1 + 415);
    const auto *osf1_416 = buffer.data(osf1 + 416);
    const auto *osf1_418 = buffer.data(osf1 + 418);
    const auto *osf1_419 = buffer.data(osf1 + 419);
    const auto *osf1_420 = buffer.data(osf1 + 420);
    const auto *osf1_423 = buffer.data(osf1 + 423);
    const auto *osf1_425 = buffer.data(osf1 + 425);
    const auto *osf1_426 = buffer.data(osf1 + 426);
    const auto *osf1_428 = buffer.data(osf1 + 428);
    const auto *osf1_429 = buffer.data(osf1 + 429);
    const auto *osf1_436 = buffer.data(osf1 + 436);
    const auto *osf1_438 = buffer.data(osf1 + 438);
    const auto *osf1_439 = buffer.data(osf1 + 439);
    const auto *osf1_440 = buffer.data(osf1 + 440);
    const auto *osf1_441 = buffer.data(osf1 + 441);
    const auto *osf1_442 = buffer.data(osf1 + 442);
    const auto *osf1_445 = buffer.data(osf1 + 445);
    const auto *osf1_446 = buffer.data(osf1 + 446);
    const auto *osf1_447 = buffer.data(osf1 + 447);
    const auto *osf1_448 = buffer.data(osf1 + 448);
    const auto *osf1_449 = buffer.data(osf1 + 449);

    const auto *osg_595 = buffer.data(osg + 595);
    const auto *osg_597 = buffer.data(osg + 597);
    const auto *osg_598 = buffer.data(osg + 598);
    const auto *osg_599 = buffer.data(osg + 599);
    const auto *osg_600 = buffer.data(osg + 600);
    const auto *osg_602 = buffer.data(osg + 602);
    const auto *osg_603 = buffer.data(osg + 603);
    const auto *osg_605 = buffer.data(osg + 605);
    const auto *osg_606 = buffer.data(osg + 606);
    const auto *osg_609 = buffer.data(osg + 609);
    const auto *osg_610 = buffer.data(osg + 610);
    const auto *osg_611 = buffer.data(osg + 611);
    const auto *osg_612 = buffer.data(osg + 612);
    const auto *osg_613 = buffer.data(osg + 613);
    const auto *osg_614 = buffer.data(osg + 614);
    const auto *osg_615 = buffer.data(osg + 615);
    const auto *osg_617 = buffer.data(osg + 617);
    const auto *osg_618 = buffer.data(osg + 618);
    const auto *osg_620 = buffer.data(osg + 620);
    const auto *osg_621 = buffer.data(osg + 621);
    const auto *osg_624 = buffer.data(osg + 624);
    const auto *osg_625 = buffer.data(osg + 625);
    const auto *osg_626 = buffer.data(osg + 626);
    const auto *osg_627 = buffer.data(osg + 627);
    const auto *osg_628 = buffer.data(osg + 628);
    const auto *osg_629 = buffer.data(osg + 629);
    const auto *osg_630 = buffer.data(osg + 630);
    const auto *osg_632 = buffer.data(osg + 632);
    const auto *osg_633 = buffer.data(osg + 633);
    const auto *osg_635 = buffer.data(osg + 635);
    const auto *osg_636 = buffer.data(osg + 636);
    const auto *osg_639 = buffer.data(osg + 639);
    const auto *osg_640 = buffer.data(osg + 640);
    const auto *osg_641 = buffer.data(osg + 641);
    const auto *osg_642 = buffer.data(osg + 642);
    const auto *osg_643 = buffer.data(osg + 643);
    const auto *osg_644 = buffer.data(osg + 644);
    const auto *osg_645 = buffer.data(osg + 645);
    const auto *osg_647 = buffer.data(osg + 647);
    const auto *osg_648 = buffer.data(osg + 648);
    const auto *osg_650 = buffer.data(osg + 650);
    const auto *osg_655 = buffer.data(osg + 655);
    const auto *osg_656 = buffer.data(osg + 656);
    const auto *osg_657 = buffer.data(osg + 657);
    const auto *osg_658 = buffer.data(osg + 658);
    const auto *osg_659 = buffer.data(osg + 659);
    const auto *osg_660 = buffer.data(osg + 660);
    const auto *osg_661 = buffer.data(osg + 661);
    const auto *osg_662 = buffer.data(osg + 662);
    const auto *osg_663 = buffer.data(osg + 663);
    const auto *osg_664 = buffer.data(osg + 664);
    const auto *osg_665 = buffer.data(osg + 665);
    const auto *osg_669 = buffer.data(osg + 669);
    const auto *osg_670 = buffer.data(osg + 670);
    const auto *osg_671 = buffer.data(osg + 671);
    const auto *osg_672 = buffer.data(osg + 672);
    const auto *osg_673 = buffer.data(osg + 673);
    const auto *osg_674 = buffer.data(osg + 674);

#pragma omp simd aligned(t_831, t_832, t_833, t_834, pc_x, pc_y, nsg_475, nsg_597, nsg_598, \
                         nsg_599, osf0_396, osf1_396, osg_595, osg_597, osg_598, \
                         osg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_11 * nsg_597[k]
                   + f_3 * pc_x[k] * osg_597[k];

        t_832[k] = f_11 * nsg_598[k]
                   + f_3 * pc_x[k] * osg_598[k];

        t_833[k] = f_11 * nsg_599[k]
                   + f_3 * pc_x[k] * osg_599[k];

        t_834[k] = f_20 * nsg_475[k]
                   + f_1 * osf0_396[k]
                   - f_2 * osf1_396[k]
                   + f_3 * pc_y[k] * osg_595[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pc_y, pc_z, nsg_460, nsg_477, nsg_478, osf0_398, \
                         osf0_399, osf1_398, osf1_399, osg_595, osg_597, \
                         osg_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_11 * nsg_460[k]
                   + f_3 * pc_z[k] * osg_595[k];

        t_836[k] = f_20 * nsg_477[k]
                   + f_6 * osf0_398[k]
                   - f_7 * osf1_398[k]
                   + f_3 * pc_y[k] * osg_597[k];

        t_837[k] = f_20 * nsg_478[k]
                   + f_4 * osf0_399[k]
                   - f_5 * osf1_399[k]
                   + f_3 * pc_y[k] * osg_598[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pc_x, pc_y, pc_z, nsg_464, nsg_479, nsg_600, \
                         osf0_399, osf0_400, osf1_399, osf1_400, osg_599, \
                         osg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_20 * nsg_479[k]
                   + f_3 * pc_y[k] * osg_599[k];

        t_839[k] = f_11 * nsg_464[k]
                   + f_1 * osf0_399[k]
                   - f_2 * osf1_399[k]
                   + f_3 * pc_z[k] * osg_599[k];

        t_840[k] = f_11 * nsg_600[k]
                   + f_1 * osf0_400[k]
                   - f_2 * osf1_400[k]
                   + f_3 * pc_x[k] * osg_600[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pc_x, pc_y, pc_z, nsg_465, nsg_480, \
                         nsg_482, nsg_603, osf0_403, osf1_403, osg_600, osg_602, \
                         osg_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_18 * nsg_480[k]
                   + f_3 * pc_y[k] * osg_600[k];

        t_842[k] = f_18 * nsg_465[k]
                   + f_3 * pc_z[k] * osg_600[k];

        t_843[k] = f_11 * nsg_603[k]
                   + f_6 * osf0_403[k]
                   - f_7 * osf1_403[k]
                   + f_3 * pc_x[k] * osg_603[k];

        t_844[k] = f_18 * nsg_482[k]
                   + f_3 * pc_y[k] * osg_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_z, nsg_468, nsg_605, nsg_606, osf0_405, \
                         osf0_406, osf1_405, osf1_406, osg_603, osg_605, \
                         osg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_11 * nsg_605[k]
                   + f_6 * osf0_405[k]
                   - f_7 * osf1_405[k]
                   + f_3 * pc_x[k] * osg_605[k];

        t_846[k] = f_11 * nsg_606[k]
                   + f_4 * osf0_406[k]
                   - f_5 * osf1_406[k]
                   + f_3 * pc_x[k] * osg_606[k];

        t_847[k] = f_18 * nsg_468[k]
                   + f_3 * pc_z[k] * osg_603[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, pc_y, nsg_485, nsg_609, nsg_610, \
                         nsg_611, osf0_409, osf1_409, osg_605, osg_609, osg_610, \
                         osg_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_18 * nsg_485[k]
                   + f_3 * pc_y[k] * osg_605[k];

        t_849[k] = f_11 * nsg_609[k]
                   + f_4 * osf0_409[k]
                   - f_5 * osf1_409[k]
                   + f_3 * pc_x[k] * osg_609[k];

        t_850[k] = f_11 * nsg_610[k]
                   + f_3 * pc_x[k] * osg_610[k];

        t_851[k] = f_11 * nsg_611[k]
                   + f_3 * pc_x[k] * osg_611[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, pc_y, nsg_490, nsg_612, nsg_613, \
                         nsg_614, osf0_406, osf1_406, osg_610, osg_612, osg_613, \
                         osg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_11 * nsg_612[k]
                   + f_3 * pc_x[k] * osg_612[k];

        t_853[k] = f_11 * nsg_613[k]
                   + f_3 * pc_x[k] * osg_613[k];

        t_854[k] = f_11 * nsg_614[k]
                   + f_3 * pc_x[k] * osg_614[k];

        t_855[k] = f_18 * nsg_490[k]
                   + f_1 * osf0_406[k]
                   - f_2 * osf1_406[k]
                   + f_3 * pc_y[k] * osg_610[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, nsg_475, nsg_492, nsg_493, osf0_408, \
                         osf0_409, osf1_408, osf1_409, osg_610, osg_612, \
                         osg_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_18 * nsg_475[k]
                   + f_3 * pc_z[k] * osg_610[k];

        t_857[k] = f_18 * nsg_492[k]
                   + f_6 * osf0_408[k]
                   - f_7 * osf1_408[k]
                   + f_3 * pc_y[k] * osg_612[k];

        t_858[k] = f_18 * nsg_493[k]
                   + f_4 * osf0_409[k]
                   - f_5 * osf1_409[k]
                   + f_3 * pc_y[k] * osg_613[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_x, pc_y, pc_z, nsg_479, nsg_494, nsg_615, \
                         osf0_409, osf0_410, osf1_409, osf1_410, osg_614, \
                         osg_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_18 * nsg_494[k]
                   + f_3 * pc_y[k] * osg_614[k];

        t_860[k] = f_18 * nsg_479[k]
                   + f_1 * osf0_409[k]
                   - f_2 * osf1_409[k]
                   + f_3 * pc_z[k] * osg_614[k];

        t_861[k] = f_11 * nsg_615[k]
                   + f_1 * osf0_410[k]
                   - f_2 * osf1_410[k]
                   + f_3 * pc_x[k] * osg_615[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pc_x, pc_y, pc_z, nsg_480, nsg_495, \
                         nsg_497, nsg_618, osf0_413, osf1_413, osg_615, osg_617, \
                         osg_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_11 * nsg_495[k]
                   + f_3 * pc_y[k] * osg_615[k];

        t_863[k] = f_20 * nsg_480[k]
                   + f_3 * pc_z[k] * osg_615[k];

        t_864[k] = f_11 * nsg_618[k]
                   + f_6 * osf0_413[k]
                   - f_7 * osf1_413[k]
                   + f_3 * pc_x[k] * osg_618[k];

        t_865[k] = f_11 * nsg_497[k]
                   + f_3 * pc_y[k] * osg_617[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_z, nsg_483, nsg_620, nsg_621, osf0_415, \
                         osf0_416, osf1_415, osf1_416, osg_618, osg_620, \
                         osg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_11 * nsg_620[k]
                   + f_6 * osf0_415[k]
                   - f_7 * osf1_415[k]
                   + f_3 * pc_x[k] * osg_620[k];

        t_867[k] = f_11 * nsg_621[k]
                   + f_4 * osf0_416[k]
                   - f_5 * osf1_416[k]
                   + f_3 * pc_x[k] * osg_621[k];

        t_868[k] = f_20 * nsg_483[k]
                   + f_3 * pc_z[k] * osg_618[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, nsg_500, nsg_624, nsg_625, \
                         nsg_626, osf0_419, osf1_419, osg_620, osg_624, osg_625, \
                         osg_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_11 * nsg_500[k]
                   + f_3 * pc_y[k] * osg_620[k];

        t_870[k] = f_11 * nsg_624[k]
                   + f_4 * osf0_419[k]
                   - f_5 * osf1_419[k]
                   + f_3 * pc_x[k] * osg_624[k];

        t_871[k] = f_11 * nsg_625[k]
                   + f_3 * pc_x[k] * osg_625[k];

        t_872[k] = f_11 * nsg_626[k]
                   + f_3 * pc_x[k] * osg_626[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pc_x, pc_y, nsg_505, nsg_627, nsg_628, \
                         nsg_629, osf0_416, osf1_416, osg_625, osg_627, osg_628, \
                         osg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_11 * nsg_627[k]
                   + f_3 * pc_x[k] * osg_627[k];

        t_874[k] = f_11 * nsg_628[k]
                   + f_3 * pc_x[k] * osg_628[k];

        t_875[k] = f_11 * nsg_629[k]
                   + f_3 * pc_x[k] * osg_629[k];

        t_876[k] = f_11 * nsg_505[k]
                   + f_1 * osf0_416[k]
                   - f_2 * osf1_416[k]
                   + f_3 * pc_y[k] * osg_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, nsg_490, nsg_507, nsg_508, osf0_418, \
                         osf0_419, osf1_418, osf1_419, osg_625, osg_627, \
                         osg_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_20 * nsg_490[k]
                   + f_3 * pc_z[k] * osg_625[k];

        t_878[k] = f_11 * nsg_507[k]
                   + f_6 * osf0_418[k]
                   - f_7 * osf1_418[k]
                   + f_3 * pc_y[k] * osg_627[k];

        t_879[k] = f_11 * nsg_508[k]
                   + f_4 * osf0_419[k]
                   - f_5 * osf1_419[k]
                   + f_3 * pc_y[k] * osg_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_x, pc_y, pc_z, nsg_494, nsg_509, nsg_630, \
                         osf0_419, osf0_420, osf1_419, osf1_420, osg_629, \
                         osg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * nsg_509[k]
                   + f_3 * pc_y[k] * osg_629[k];

        t_881[k] = f_20 * nsg_494[k]
                   + f_1 * osf0_419[k]
                   - f_2 * osf1_419[k]
                   + f_3 * pc_z[k] * osg_629[k];

        t_882[k] = f_11 * nsg_630[k]
                   + f_1 * osf0_420[k]
                   - f_2 * osf1_420[k]
                   + f_3 * pc_x[k] * osg_630[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pc_x, pc_y, pc_z, nsg_495, nsg_510, \
                         nsg_512, nsg_633, osf0_423, osf1_423, osg_630, osg_632, \
                         osg_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_10 * nsg_510[k]
                   + f_3 * pc_y[k] * osg_630[k];

        t_884[k] = f_19 * nsg_495[k]
                   + f_3 * pc_z[k] * osg_630[k];

        t_885[k] = f_11 * nsg_633[k]
                   + f_6 * osf0_423[k]
                   - f_7 * osf1_423[k]
                   + f_3 * pc_x[k] * osg_633[k];

        t_886[k] = f_10 * nsg_512[k]
                   + f_3 * pc_y[k] * osg_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pc_x, pc_z, nsg_498, nsg_635, nsg_636, osf0_425, \
                         osf0_426, osf1_425, osf1_426, osg_633, osg_635, \
                         osg_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_11 * nsg_635[k]
                   + f_6 * osf0_425[k]
                   - f_7 * osf1_425[k]
                   + f_3 * pc_x[k] * osg_635[k];

        t_888[k] = f_11 * nsg_636[k]
                   + f_4 * osf0_426[k]
                   - f_5 * osf1_426[k]
                   + f_3 * pc_x[k] * osg_636[k];

        t_889[k] = f_19 * nsg_498[k]
                   + f_3 * pc_z[k] * osg_633[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pc_x, pc_y, nsg_515, nsg_639, nsg_640, \
                         nsg_641, osf0_429, osf1_429, osg_635, osg_639, osg_640, \
                         osg_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_10 * nsg_515[k]
                   + f_3 * pc_y[k] * osg_635[k];

        t_891[k] = f_11 * nsg_639[k]
                   + f_4 * osf0_429[k]
                   - f_5 * osf1_429[k]
                   + f_3 * pc_x[k] * osg_639[k];

        t_892[k] = f_11 * nsg_640[k]
                   + f_3 * pc_x[k] * osg_640[k];

        t_893[k] = f_11 * nsg_641[k]
                   + f_3 * pc_x[k] * osg_641[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pc_x, pc_y, nsg_520, nsg_642, nsg_643, \
                         nsg_644, osf0_426, osf1_426, osg_640, osg_642, osg_643, \
                         osg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_11 * nsg_642[k]
                   + f_3 * pc_x[k] * osg_642[k];

        t_895[k] = f_11 * nsg_643[k]
                   + f_3 * pc_x[k] * osg_643[k];

        t_896[k] = f_11 * nsg_644[k]
                   + f_3 * pc_x[k] * osg_644[k];

        t_897[k] = f_10 * nsg_520[k]
                   + f_1 * osf0_426[k]
                   - f_2 * osf1_426[k]
                   + f_3 * pc_y[k] * osg_640[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_y, pc_z, nsg_505, nsg_522, nsg_523, osf0_428, \
                         osf0_429, osf1_428, osf1_429, osg_640, osg_642, \
                         osg_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_19 * nsg_505[k]
                   + f_3 * pc_z[k] * osg_640[k];

        t_899[k] = f_10 * nsg_522[k]
                   + f_6 * osf0_428[k]
                   - f_7 * osf1_428[k]
                   + f_3 * pc_y[k] * osg_642[k];

        t_900[k] = f_10 * nsg_523[k]
                   + f_4 * osf0_429[k]
                   - f_5 * osf1_429[k]
                   + f_3 * pc_y[k] * osg_643[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pa_y, pc_y, pc_z, nsh0_735, nsg_509, \
                         nsg_524, nsg_525, nsh1_735, osf0_429, osf1_429, osg_644, \
                         osg_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_10 * nsg_524[k]
                   + f_3 * pc_y[k] * osg_644[k];

        t_902[k] = f_19 * nsg_509[k]
                   + f_1 * osf0_429[k]
                   - f_2 * osf1_429[k]
                   + f_3 * pc_z[k] * osg_644[k];

        t_903[k] = pa_y[k] * nsh0_735[k]
                   - f_8 * pc_y[k] * nsh1_735[k];

        t_904[k] = f_9 * nsg_525[k]
                   + f_3 * pc_y[k] * osg_645[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pa_y, pc_y, pc_z, nsh0_738, nsh0_740, \
                         nsg_510, nsg_526, nsg_527, nsh1_738, nsh1_740, osg_645, \
                         osg_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_17 * nsg_510[k]
                   + f_3 * pc_z[k] * osg_645[k];

        t_906[k] = pa_y[k] * nsh0_738[k]
                   + f_10 * nsg_526[k]
                   - f_8 * pc_y[k] * nsh1_738[k];

        t_907[k] = f_9 * nsg_527[k]
                   + f_3 * pc_y[k] * osg_647[k];

        t_908[k] = pa_y[k] * nsh0_740[k]
                   - f_8 * pc_y[k] * nsh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pa_y, pc_y, pc_z, nsh0_741, nsh0_744, \
                         nsg_513, nsg_528, nsg_530, nsh1_741, nsh1_744, osg_648, \
                         osg_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pa_y[k] * nsh0_741[k]
                   + f_11 * nsg_528[k]
                   - f_8 * pc_y[k] * nsh1_741[k];

        t_910[k] = f_17 * nsg_513[k]
                   + f_3 * pc_z[k] * osg_648[k];

        t_911[k] = f_9 * nsg_530[k]
                   + f_3 * pc_y[k] * osg_650[k];

        t_912[k] = pa_y[k] * nsh0_744[k]
                   - f_8 * pc_y[k] * nsh1_744[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, pc_x, nsg_655, nsg_656, nsg_657, \
                         nsg_658, nsg_659, osg_655, osg_656, osg_657, osg_658, \
                         osg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_11 * nsg_655[k]
                   + f_3 * pc_x[k] * osg_655[k];

        t_914[k] = f_11 * nsg_656[k]
                   + f_3 * pc_x[k] * osg_656[k];

        t_915[k] = f_11 * nsg_657[k]
                   + f_3 * pc_x[k] * osg_657[k];

        t_916[k] = f_11 * nsg_658[k]
                   + f_3 * pc_x[k] * osg_658[k];

        t_917[k] = f_11 * nsg_659[k]
                   + f_3 * pc_x[k] * osg_659[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_y, pc_z, nsg_520, nsg_535, nsg_537, osf0_436, \
                         osf0_438, osf1_436, osf1_438, osg_655, \
                         osg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_9 * nsg_535[k]
                   + f_1 * osf0_436[k]
                   - f_2 * osf1_436[k]
                   + f_3 * pc_y[k] * osg_655[k];

        t_919[k] = f_17 * nsg_520[k]
                   + f_3 * pc_z[k] * osg_655[k];

        t_920[k] = f_9 * nsg_537[k]
                   + f_6 * osf0_438[k]
                   - f_7 * osf1_438[k]
                   + f_3 * pc_y[k] * osg_657[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pa_y, pc_y, nsh0_755, nsg_538, nsg_539, \
                         nsh1_755, osf0_439, osf1_439, osg_658, \
                         osg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_9 * nsg_538[k]
                   + f_4 * osf0_439[k]
                   - f_5 * osf1_439[k]
                   + f_3 * pc_y[k] * osg_658[k];

        t_922[k] = f_9 * nsg_539[k]
                   + f_3 * pc_y[k] * osg_659[k];

        t_923[k] = pa_y[k] * nsh0_755[k]
                   - f_8 * pc_y[k] * nsh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, nsg_525, \
                         nsg_660, osf0_440, osf1_440, osg_660, osg_661, \
                         osg_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_11 * nsg_660[k]
                   + f_1 * osf0_440[k]
                   - f_2 * osf1_440[k]
                   + f_3 * pc_x[k] * osg_660[k];

        t_925[k] = f_3 * pc_y[k] * osg_660[k];

        t_926[k] = f_16 * nsg_525[k]
                   + f_3 * pc_z[k] * osg_660[k];

        t_927[k] = f_4 * osf0_440[k]
                   - f_5 * osf1_440[k]
                   + f_3 * pc_y[k] * osg_661[k];

        t_928[k] = f_3 * pc_y[k] * osg_662[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, t_932, pc_x, pc_y, nsg_665, osf0_441, osf0_442, \
                         osf0_445, osf1_441, osf1_442, osf1_445, osg_663, osg_664, \
                         osg_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_11 * nsg_665[k]
                   + f_6 * osf0_445[k]
                   - f_7 * osf1_445[k]
                   + f_3 * pc_x[k] * osg_665[k];

        t_930[k] = f_6 * osf0_441[k]
                   - f_7 * osf1_441[k]
                   + f_3 * pc_y[k] * osg_663[k];

        t_931[k] = f_4 * osf0_442[k]
                   - f_5 * osf1_442[k]
                   + f_3 * pc_y[k] * osg_664[k];

        t_932[k] = f_3 * pc_y[k] * osg_665[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, pc_x, nsg_669, nsg_670, nsg_671, nsg_672, \
                         osf0_449, osf1_449, osg_669, osg_670, osg_671, \
                         osg_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_11 * nsg_669[k]
                   + f_4 * osf0_449[k]
                   - f_5 * osf1_449[k]
                   + f_3 * pc_x[k] * osg_669[k];

        t_934[k] = f_11 * nsg_670[k]
                   + f_3 * pc_x[k] * osg_670[k];

        t_935[k] = f_11 * nsg_671[k]
                   + f_3 * pc_x[k] * osg_671[k];

        t_936[k] = f_11 * nsg_672[k]
                   + f_3 * pc_x[k] * osg_672[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, pc_x, pc_y, nsg_674, osf0_446, osf0_447, \
                         osf1_446, osf1_447, osg_669, osg_670, osg_671, \
                         osg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = f_3 * pc_y[k] * osg_669[k];

        t_938[k] = f_11 * nsg_674[k]
                   + f_3 * pc_x[k] * osg_674[k];

        t_939[k] = f_1 * osf0_446[k]
                   - f_2 * osf1_446[k]
                   + f_3 * pc_y[k] * osg_670[k];

        t_940[k] = f_13 * osf0_447[k]
                   - f_14 * osf1_447[k]
                   + f_3 * pc_y[k] * osg_671[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, t_944, pc_y, pc_z, nsg_539, osf0_448, osf0_449, \
                         osf1_448, osf1_449, osg_672, osg_673, \
                         osg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = f_6 * osf0_448[k]
                   - f_7 * osf1_448[k]
                   + f_3 * pc_y[k] * osg_672[k];

        t_942[k] = f_4 * osf0_449[k]
                   - f_5 * osf1_449[k]
                   + f_3 * pc_y[k] * osg_673[k];

        t_943[k] = f_3 * pc_y[k] * osg_674[k];

        t_944[k] = f_16 * nsg_539[k]
                   + f_1 * osf0_449[k]
                   - f_2 * osf1_449[k]
                   + f_3 * pc_z[k] * osg_674[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_756 = buffer.data(nsh0 + 756);
    const auto *nsh0_759 = buffer.data(nsh0 + 759);
    const auto *nsh0_762 = buffer.data(nsh0 + 762);
    const auto *nsh0_771 = buffer.data(nsh0 + 771);

    const auto *nsg_540 = buffer.data(nsg + 540);
    const auto *nsg_543 = buffer.data(nsg + 543);
    const auto *nsg_545 = buffer.data(nsg + 545);
    const auto *nsg_550 = buffer.data(nsg + 550);
    const auto *nsg_554 = buffer.data(nsg + 554);
    const auto *nsg_555 = buffer.data(nsg + 555);
    const auto *nsg_557 = buffer.data(nsg + 557);
    const auto *nsg_558 = buffer.data(nsg + 558);
    const auto *nsg_560 = buffer.data(nsg + 560);
    const auto *nsg_565 = buffer.data(nsg + 565);
    const auto *nsg_567 = buffer.data(nsg + 567);
    const auto *nsg_568 = buffer.data(nsg + 568);
    const auto *nsg_569 = buffer.data(nsg + 569);
    const auto *nsg_570 = buffer.data(nsg + 570);
    const auto *nsg_572 = buffer.data(nsg + 572);
    const auto *nsg_573 = buffer.data(nsg + 573);
    const auto *nsg_575 = buffer.data(nsg + 575);
    const auto *nsg_580 = buffer.data(nsg + 580);
    const auto *nsg_582 = buffer.data(nsg + 582);
    const auto *nsg_583 = buffer.data(nsg + 583);
    const auto *nsg_584 = buffer.data(nsg + 584);
    const auto *nsg_585 = buffer.data(nsg + 585);
    const auto *nsg_587 = buffer.data(nsg + 587);
    const auto *nsg_588 = buffer.data(nsg + 588);
    const auto *nsg_590 = buffer.data(nsg + 590);
    const auto *nsg_595 = buffer.data(nsg + 595);
    const auto *nsg_597 = buffer.data(nsg + 597);
    const auto *nsg_598 = buffer.data(nsg + 598);
    const auto *nsg_599 = buffer.data(nsg + 599);
    const auto *nsg_600 = buffer.data(nsg + 600);
    const auto *nsg_602 = buffer.data(nsg + 602);
    const auto *nsg_603 = buffer.data(nsg + 603);
    const auto *nsg_605 = buffer.data(nsg + 605);
    const auto *nsg_610 = buffer.data(nsg + 610);
    const auto *nsg_612 = buffer.data(nsg + 612);
    const auto *nsg_613 = buffer.data(nsg + 613);
    const auto *nsg_614 = buffer.data(nsg + 614);
    const auto *nsg_615 = buffer.data(nsg + 615);
    const auto *nsg_617 = buffer.data(nsg + 617);
    const auto *nsg_675 = buffer.data(nsg + 675);
    const auto *nsg_678 = buffer.data(nsg + 678);
    const auto *nsg_681 = buffer.data(nsg + 681);
    const auto *nsg_685 = buffer.data(nsg + 685);
    const auto *nsg_687 = buffer.data(nsg + 687);
    const auto *nsg_688 = buffer.data(nsg + 688);
    const auto *nsg_689 = buffer.data(nsg + 689);
    const auto *nsg_695 = buffer.data(nsg + 695);
    const auto *nsg_699 = buffer.data(nsg + 699);
    const auto *nsg_700 = buffer.data(nsg + 700);
    const auto *nsg_701 = buffer.data(nsg + 701);
    const auto *nsg_702 = buffer.data(nsg + 702);
    const auto *nsg_703 = buffer.data(nsg + 703);
    const auto *nsg_704 = buffer.data(nsg + 704);
    const auto *nsg_705 = buffer.data(nsg + 705);
    const auto *nsg_708 = buffer.data(nsg + 708);
    const auto *nsg_710 = buffer.data(nsg + 710);
    const auto *nsg_711 = buffer.data(nsg + 711);
    const auto *nsg_714 = buffer.data(nsg + 714);
    const auto *nsg_715 = buffer.data(nsg + 715);
    const auto *nsg_716 = buffer.data(nsg + 716);
    const auto *nsg_717 = buffer.data(nsg + 717);
    const auto *nsg_718 = buffer.data(nsg + 718);
    const auto *nsg_719 = buffer.data(nsg + 719);
    const auto *nsg_720 = buffer.data(nsg + 720);
    const auto *nsg_723 = buffer.data(nsg + 723);
    const auto *nsg_725 = buffer.data(nsg + 725);
    const auto *nsg_726 = buffer.data(nsg + 726);
    const auto *nsg_729 = buffer.data(nsg + 729);
    const auto *nsg_730 = buffer.data(nsg + 730);
    const auto *nsg_731 = buffer.data(nsg + 731);
    const auto *nsg_732 = buffer.data(nsg + 732);
    const auto *nsg_733 = buffer.data(nsg + 733);
    const auto *nsg_734 = buffer.data(nsg + 734);
    const auto *nsg_735 = buffer.data(nsg + 735);
    const auto *nsg_738 = buffer.data(nsg + 738);
    const auto *nsg_740 = buffer.data(nsg + 740);
    const auto *nsg_741 = buffer.data(nsg + 741);
    const auto *nsg_744 = buffer.data(nsg + 744);
    const auto *nsg_745 = buffer.data(nsg + 745);
    const auto *nsg_746 = buffer.data(nsg + 746);
    const auto *nsg_747 = buffer.data(nsg + 747);
    const auto *nsg_748 = buffer.data(nsg + 748);
    const auto *nsg_749 = buffer.data(nsg + 749);
    const auto *nsg_750 = buffer.data(nsg + 750);
    const auto *nsg_753 = buffer.data(nsg + 753);
    const auto *nsg_755 = buffer.data(nsg + 755);
    const auto *nsg_756 = buffer.data(nsg + 756);

    const auto *nsh1_756 = buffer.data(nsh1 + 756);
    const auto *nsh1_759 = buffer.data(nsh1 + 759);
    const auto *nsh1_762 = buffer.data(nsh1 + 762);
    const auto *nsh1_771 = buffer.data(nsh1 + 771);

    const auto *osf0_450 = buffer.data(osf0 + 450);
    const auto *osf0_452 = buffer.data(osf0 + 452);
    const auto *osf0_453 = buffer.data(osf0 + 453);
    const auto *osf0_456 = buffer.data(osf0 + 456);
    const auto *osf0_457 = buffer.data(osf0 + 457);
    const auto *osf0_459 = buffer.data(osf0 + 459);
    const auto *osf0_465 = buffer.data(osf0 + 465);
    const auto *osf0_468 = buffer.data(osf0 + 468);
    const auto *osf0_469 = buffer.data(osf0 + 469);
    const auto *osf0_470 = buffer.data(osf0 + 470);
    const auto *osf0_473 = buffer.data(osf0 + 473);
    const auto *osf0_475 = buffer.data(osf0 + 475);
    const auto *osf0_476 = buffer.data(osf0 + 476);
    const auto *osf0_478 = buffer.data(osf0 + 478);
    const auto *osf0_479 = buffer.data(osf0 + 479);
    const auto *osf0_480 = buffer.data(osf0 + 480);
    const auto *osf0_483 = buffer.data(osf0 + 483);
    const auto *osf0_485 = buffer.data(osf0 + 485);
    const auto *osf0_486 = buffer.data(osf0 + 486);
    const auto *osf0_488 = buffer.data(osf0 + 488);
    const auto *osf0_489 = buffer.data(osf0 + 489);
    const auto *osf0_490 = buffer.data(osf0 + 490);
    const auto *osf0_493 = buffer.data(osf0 + 493);
    const auto *osf0_495 = buffer.data(osf0 + 495);
    const auto *osf0_496 = buffer.data(osf0 + 496);
    const auto *osf0_498 = buffer.data(osf0 + 498);
    const auto *osf0_499 = buffer.data(osf0 + 499);
    const auto *osf0_500 = buffer.data(osf0 + 500);
    const auto *osf0_503 = buffer.data(osf0 + 503);
    const auto *osf0_505 = buffer.data(osf0 + 505);
    const auto *osf0_506 = buffer.data(osf0 + 506);

    const auto *osf1_450 = buffer.data(osf1 + 450);
    const auto *osf1_452 = buffer.data(osf1 + 452);
    const auto *osf1_453 = buffer.data(osf1 + 453);
    const auto *osf1_456 = buffer.data(osf1 + 456);
    const auto *osf1_457 = buffer.data(osf1 + 457);
    const auto *osf1_459 = buffer.data(osf1 + 459);
    const auto *osf1_465 = buffer.data(osf1 + 465);
    const auto *osf1_468 = buffer.data(osf1 + 468);
    const auto *osf1_469 = buffer.data(osf1 + 469);
    const auto *osf1_470 = buffer.data(osf1 + 470);
    const auto *osf1_473 = buffer.data(osf1 + 473);
    const auto *osf1_475 = buffer.data(osf1 + 475);
    const auto *osf1_476 = buffer.data(osf1 + 476);
    const auto *osf1_478 = buffer.data(osf1 + 478);
    const auto *osf1_479 = buffer.data(osf1 + 479);
    const auto *osf1_480 = buffer.data(osf1 + 480);
    const auto *osf1_483 = buffer.data(osf1 + 483);
    const auto *osf1_485 = buffer.data(osf1 + 485);
    const auto *osf1_486 = buffer.data(osf1 + 486);
    const auto *osf1_488 = buffer.data(osf1 + 488);
    const auto *osf1_489 = buffer.data(osf1 + 489);
    const auto *osf1_490 = buffer.data(osf1 + 490);
    const auto *osf1_493 = buffer.data(osf1 + 493);
    const auto *osf1_495 = buffer.data(osf1 + 495);
    const auto *osf1_496 = buffer.data(osf1 + 496);
    const auto *osf1_498 = buffer.data(osf1 + 498);
    const auto *osf1_499 = buffer.data(osf1 + 499);
    const auto *osf1_500 = buffer.data(osf1 + 500);
    const auto *osf1_503 = buffer.data(osf1 + 503);
    const auto *osf1_505 = buffer.data(osf1 + 505);
    const auto *osf1_506 = buffer.data(osf1 + 506);

    const auto *osg_675 = buffer.data(osg + 675);
    const auto *osg_676 = buffer.data(osg + 676);
    const auto *osg_677 = buffer.data(osg + 677);
    const auto *osg_678 = buffer.data(osg + 678);
    const auto *osg_680 = buffer.data(osg + 680);
    const auto *osg_681 = buffer.data(osg + 681);
    const auto *osg_685 = buffer.data(osg + 685);
    const auto *osg_686 = buffer.data(osg + 686);
    const auto *osg_687 = buffer.data(osg + 687);
    const auto *osg_688 = buffer.data(osg + 688);
    const auto *osg_689 = buffer.data(osg + 689);
    const auto *osg_690 = buffer.data(osg + 690);
    const auto *osg_692 = buffer.data(osg + 692);
    const auto *osg_693 = buffer.data(osg + 693);
    const auto *osg_695 = buffer.data(osg + 695);
    const auto *osg_699 = buffer.data(osg + 699);
    const auto *osg_700 = buffer.data(osg + 700);
    const auto *osg_701 = buffer.data(osg + 701);
    const auto *osg_702 = buffer.data(osg + 702);
    const auto *osg_703 = buffer.data(osg + 703);
    const auto *osg_704 = buffer.data(osg + 704);
    const auto *osg_705 = buffer.data(osg + 705);
    const auto *osg_707 = buffer.data(osg + 707);
    const auto *osg_708 = buffer.data(osg + 708);
    const auto *osg_710 = buffer.data(osg + 710);
    const auto *osg_711 = buffer.data(osg + 711);
    const auto *osg_714 = buffer.data(osg + 714);
    const auto *osg_715 = buffer.data(osg + 715);
    const auto *osg_716 = buffer.data(osg + 716);
    const auto *osg_717 = buffer.data(osg + 717);
    const auto *osg_718 = buffer.data(osg + 718);
    const auto *osg_719 = buffer.data(osg + 719);
    const auto *osg_720 = buffer.data(osg + 720);
    const auto *osg_722 = buffer.data(osg + 722);
    const auto *osg_723 = buffer.data(osg + 723);
    const auto *osg_725 = buffer.data(osg + 725);
    const auto *osg_726 = buffer.data(osg + 726);
    const auto *osg_729 = buffer.data(osg + 729);
    const auto *osg_730 = buffer.data(osg + 730);
    const auto *osg_731 = buffer.data(osg + 731);
    const auto *osg_732 = buffer.data(osg + 732);
    const auto *osg_733 = buffer.data(osg + 733);
    const auto *osg_734 = buffer.data(osg + 734);
    const auto *osg_735 = buffer.data(osg + 735);
    const auto *osg_737 = buffer.data(osg + 737);
    const auto *osg_738 = buffer.data(osg + 738);
    const auto *osg_740 = buffer.data(osg + 740);
    const auto *osg_741 = buffer.data(osg + 741);
    const auto *osg_744 = buffer.data(osg + 744);
    const auto *osg_745 = buffer.data(osg + 745);
    const auto *osg_746 = buffer.data(osg + 746);
    const auto *osg_747 = buffer.data(osg + 747);
    const auto *osg_748 = buffer.data(osg + 748);
    const auto *osg_749 = buffer.data(osg + 749);
    const auto *osg_750 = buffer.data(osg + 750);
    const auto *osg_752 = buffer.data(osg + 752);
    const auto *osg_753 = buffer.data(osg + 753);
    const auto *osg_755 = buffer.data(osg + 755);
    const auto *osg_756 = buffer.data(osg + 756);

#pragma omp simd aligned(t_945, t_946, t_947, t_948, pc_x, pc_y, pc_z, nsg_540, nsg_675, \
                         nsg_678, osf0_450, osf0_453, osf1_450, osf1_453, osg_675, \
                         osg_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_10 * nsg_675[k]
                   + f_1 * osf0_450[k]
                   - f_2 * osf1_450[k]
                   + f_3 * pc_x[k] * osg_675[k];

        t_946[k] = f_15 * nsg_540[k]
                   + f_3 * pc_y[k] * osg_675[k];

        t_947[k] = f_3 * pc_z[k] * osg_675[k];

        t_948[k] = f_10 * nsg_678[k]
                   + f_6 * osf0_453[k]
                   - f_7 * osf1_453[k]
                   + f_3 * pc_x[k] * osg_678[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pc_x, pc_z, nsg_681, osf0_450, osf0_456, \
                         osf1_450, osf1_456, osg_676, osg_677, osg_678, \
                         osg_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_3 * pc_z[k] * osg_676[k];

        t_950[k] = f_4 * osf0_450[k]
                   - f_5 * osf1_450[k]
                   + f_3 * pc_z[k] * osg_677[k];

        t_951[k] = f_10 * nsg_681[k]
                   + f_4 * osf0_456[k]
                   - f_5 * osf1_456[k]
                   + f_3 * pc_x[k] * osg_681[k];

        t_952[k] = f_3 * pc_z[k] * osg_678[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pc_x, pc_y, pc_z, nsg_545, nsg_685, \
                         osf0_452, osf1_452, osg_680, osg_681, \
                         osg_685 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_15 * nsg_545[k]
                   + f_3 * pc_y[k] * osg_680[k];

        t_954[k] = f_6 * osf0_452[k]
                   - f_7 * osf1_452[k]
                   + f_3 * pc_z[k] * osg_680[k];

        t_955[k] = f_10 * nsg_685[k]
                   + f_3 * pc_x[k] * osg_685[k];

        t_956[k] = f_3 * pc_z[k] * osg_681[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pc_x, pc_y, nsg_550, nsg_687, nsg_688, \
                         nsg_689, osf0_456, osf1_456, osg_685, osg_687, osg_688, \
                         osg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_10 * nsg_687[k]
                   + f_3 * pc_x[k] * osg_687[k];

        t_958[k] = f_10 * nsg_688[k]
                   + f_3 * pc_x[k] * osg_688[k];

        t_959[k] = f_10 * nsg_689[k]
                   + f_3 * pc_x[k] * osg_689[k];

        t_960[k] = f_15 * nsg_550[k]
                   + f_1 * osf0_456[k]
                   - f_2 * osf1_456[k]
                   + f_3 * pc_y[k] * osg_685[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pc_y, pc_z, nsg_554, osf0_456, osf0_457, \
                         osf1_456, osf1_457, osg_685, osg_686, osg_687, \
                         osg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * osg_685[k];

        t_962[k] = f_4 * osf0_456[k]
                   - f_5 * osf1_456[k]
                   + f_3 * pc_z[k] * osg_686[k];

        t_963[k] = f_6 * osf0_457[k]
                   - f_7 * osf1_457[k]
                   + f_3 * pc_z[k] * osg_687[k];

        t_964[k] = f_15 * nsg_554[k]
                   + f_3 * pc_y[k] * osg_689[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pa_z, pc_y, pc_z, nsh0_756, nsg_540, \
                         nsg_555, nsh1_756, osf0_459, osf1_459, osg_689, \
                         osg_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_1 * osf0_459[k]
                   - f_2 * osf1_459[k]
                   + f_3 * pc_z[k] * osg_689[k];

        t_966[k] = pa_z[k] * nsh0_756[k]
                   - f_8 * pc_z[k] * nsh1_756[k];

        t_967[k] = f_16 * nsg_555[k]
                   + f_3 * pc_y[k] * osg_690[k];

        t_968[k] = f_9 * nsg_540[k]
                   + f_3 * pc_z[k] * osg_690[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_z, pc_x, pc_y, pc_z, nsh0_759, nsg_557, \
                         nsg_695, nsh1_759, osf0_465, osf1_465, osg_692, \
                         osg_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = pa_z[k] * nsh0_759[k]
                   - f_8 * pc_z[k] * nsh1_759[k];

        t_970[k] = f_16 * nsg_557[k]
                   + f_3 * pc_y[k] * osg_692[k];

        t_971[k] = f_10 * nsg_695[k]
                   + f_6 * osf0_465[k]
                   - f_7 * osf1_465[k]
                   + f_3 * pc_x[k] * osg_695[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, pa_z, pc_y, pc_z, nsh0_762, nsg_543, nsg_560, \
                         nsh1_762, osg_693, osg_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = pa_z[k] * nsh0_762[k]
                   - f_8 * pc_z[k] * nsh1_762[k];

        t_973[k] = f_9 * nsg_543[k]
                   + f_3 * pc_z[k] * osg_693[k];

        t_974[k] = f_16 * nsg_560[k]
                   + f_3 * pc_y[k] * osg_695[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, pc_x, nsg_699, nsg_700, nsg_701, nsg_702, \
                         osf0_469, osf1_469, osg_699, osg_700, osg_701, \
                         osg_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_10 * nsg_699[k]
                   + f_4 * osf0_469[k]
                   - f_5 * osf1_469[k]
                   + f_3 * pc_x[k] * osg_699[k];

        t_976[k] = f_10 * nsg_700[k]
                   + f_3 * pc_x[k] * osg_700[k];

        t_977[k] = f_10 * nsg_701[k]
                   + f_3 * pc_x[k] * osg_701[k];

        t_978[k] = f_10 * nsg_702[k]
                   + f_3 * pc_x[k] * osg_702[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_z, pc_x, pc_z, nsh0_771, nsg_550, \
                         nsg_703, nsg_704, nsh1_771, osg_700, osg_703, \
                         osg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_10 * nsg_703[k]
                   + f_3 * pc_x[k] * osg_703[k];

        t_980[k] = f_10 * nsg_704[k]
                   + f_3 * pc_x[k] * osg_704[k];

        t_981[k] = pa_z[k] * nsh0_771[k]
                   - f_8 * pc_z[k] * nsh1_771[k];

        t_982[k] = f_9 * nsg_550[k]
                   + f_3 * pc_z[k] * osg_700[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_y, nsg_567, nsg_568, nsg_569, osf0_468, \
                         osf0_469, osf1_468, osf1_469, osg_702, osg_703, \
                         osg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_16 * nsg_567[k]
                   + f_6 * osf0_468[k]
                   - f_7 * osf1_468[k]
                   + f_3 * pc_y[k] * osg_702[k];

        t_984[k] = f_16 * nsg_568[k]
                   + f_4 * osf0_469[k]
                   - f_5 * osf1_469[k]
                   + f_3 * pc_y[k] * osg_703[k];

        t_985[k] = f_16 * nsg_569[k]
                   + f_3 * pc_y[k] * osg_704[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_x, pc_y, pc_z, nsg_554, nsg_570, nsg_705, \
                         osf0_469, osf0_470, osf1_469, osf1_470, osg_704, \
                         osg_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_9 * nsg_554[k]
                   + f_1 * osf0_469[k]
                   - f_2 * osf1_469[k]
                   + f_3 * pc_z[k] * osg_704[k];

        t_987[k] = f_10 * nsg_705[k]
                   + f_1 * osf0_470[k]
                   - f_2 * osf1_470[k]
                   + f_3 * pc_x[k] * osg_705[k];

        t_988[k] = f_17 * nsg_570[k]
                   + f_3 * pc_y[k] * osg_705[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, pc_z, nsg_555, nsg_572, nsg_708, \
                         osf0_473, osf1_473, osg_705, osg_707, \
                         osg_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_10 * nsg_555[k]
                   + f_3 * pc_z[k] * osg_705[k];

        t_990[k] = f_10 * nsg_708[k]
                   + f_6 * osf0_473[k]
                   - f_7 * osf1_473[k]
                   + f_3 * pc_x[k] * osg_708[k];

        t_991[k] = f_17 * nsg_572[k]
                   + f_3 * pc_y[k] * osg_707[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, pc_z, nsg_558, nsg_710, nsg_711, osf0_475, \
                         osf0_476, osf1_475, osf1_476, osg_708, osg_710, \
                         osg_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_10 * nsg_710[k]
                   + f_6 * osf0_475[k]
                   - f_7 * osf1_475[k]
                   + f_3 * pc_x[k] * osg_710[k];

        t_993[k] = f_10 * nsg_711[k]
                   + f_4 * osf0_476[k]
                   - f_5 * osf1_476[k]
                   + f_3 * pc_x[k] * osg_711[k];

        t_994[k] = f_10 * nsg_558[k]
                   + f_3 * pc_z[k] * osg_708[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pc_x, pc_y, nsg_575, nsg_714, nsg_715, \
                         nsg_716, osf0_479, osf1_479, osg_710, osg_714, osg_715, \
                         osg_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_17 * nsg_575[k]
                   + f_3 * pc_y[k] * osg_710[k];

        t_996[k] = f_10 * nsg_714[k]
                   + f_4 * osf0_479[k]
                   - f_5 * osf1_479[k]
                   + f_3 * pc_x[k] * osg_714[k];

        t_997[k] = f_10 * nsg_715[k]
                   + f_3 * pc_x[k] * osg_715[k];

        t_998[k] = f_10 * nsg_716[k]
                   + f_3 * pc_x[k] * osg_716[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, t_1002, pc_x, pc_y, nsg_580, nsg_717, nsg_718, \
                         nsg_719, osf0_476, osf1_476, osg_715, osg_717, osg_718, \
                         osg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_10 * nsg_717[k]
                   + f_3 * pc_x[k] * osg_717[k];

        t_1000[k] = f_10 * nsg_718[k]
                    + f_3 * pc_x[k] * osg_718[k];

        t_1001[k] = f_10 * nsg_719[k]
                    + f_3 * pc_x[k] * osg_719[k];

        t_1002[k] = f_17 * nsg_580[k]
                    + f_1 * osf0_476[k]
                    - f_2 * osf1_476[k]
                    + f_3 * pc_y[k] * osg_715[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, pc_y, pc_z, nsg_565, nsg_582, nsg_583, \
                         osf0_478, osf0_479, osf1_478, osf1_479, osg_715, osg_717, \
                         osg_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_10 * nsg_565[k]
                    + f_3 * pc_z[k] * osg_715[k];

        t_1004[k] = f_17 * nsg_582[k]
                    + f_6 * osf0_478[k]
                    - f_7 * osf1_478[k]
                    + f_3 * pc_y[k] * osg_717[k];

        t_1005[k] = f_17 * nsg_583[k]
                    + f_4 * osf0_479[k]
                    - f_5 * osf1_479[k]
                    + f_3 * pc_y[k] * osg_718[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, pc_x, pc_y, pc_z, nsg_569, nsg_584, nsg_720, \
                         osf0_479, osf0_480, osf1_479, osf1_480, osg_719, \
                         osg_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_17 * nsg_584[k]
                    + f_3 * pc_y[k] * osg_719[k];

        t_1007[k] = f_10 * nsg_569[k]
                    + f_1 * osf0_479[k]
                    - f_2 * osf1_479[k]
                    + f_3 * pc_z[k] * osg_719[k];

        t_1008[k] = f_10 * nsg_720[k]
                    + f_1 * osf0_480[k]
                    - f_2 * osf1_480[k]
                    + f_3 * pc_x[k] * osg_720[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pc_x, pc_y, pc_z, nsg_570, nsg_585, \
                         nsg_587, nsg_723, osf0_483, osf1_483, osg_720, osg_722, \
                         osg_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_19 * nsg_585[k]
                    + f_3 * pc_y[k] * osg_720[k];

        t_1010[k] = f_11 * nsg_570[k]
                    + f_3 * pc_z[k] * osg_720[k];

        t_1011[k] = f_10 * nsg_723[k]
                    + f_6 * osf0_483[k]
                    - f_7 * osf1_483[k]
                    + f_3 * pc_x[k] * osg_723[k];

        t_1012[k] = f_19 * nsg_587[k]
                    + f_3 * pc_y[k] * osg_722[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, pc_z, nsg_573, nsg_725, nsg_726, \
                         osf0_485, osf0_486, osf1_485, osf1_486, osg_723, osg_725, \
                         osg_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_10 * nsg_725[k]
                    + f_6 * osf0_485[k]
                    - f_7 * osf1_485[k]
                    + f_3 * pc_x[k] * osg_725[k];

        t_1014[k] = f_10 * nsg_726[k]
                    + f_4 * osf0_486[k]
                    - f_5 * osf1_486[k]
                    + f_3 * pc_x[k] * osg_726[k];

        t_1015[k] = f_11 * nsg_573[k]
                    + f_3 * pc_z[k] * osg_723[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, nsg_590, nsg_729, \
                         nsg_730, nsg_731, osf0_489, osf1_489, osg_725, osg_729, osg_730, \
                         osg_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * nsg_590[k]
                    + f_3 * pc_y[k] * osg_725[k];

        t_1017[k] = f_10 * nsg_729[k]
                    + f_4 * osf0_489[k]
                    - f_5 * osf1_489[k]
                    + f_3 * pc_x[k] * osg_729[k];

        t_1018[k] = f_10 * nsg_730[k]
                    + f_3 * pc_x[k] * osg_730[k];

        t_1019[k] = f_10 * nsg_731[k]
                    + f_3 * pc_x[k] * osg_731[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, nsg_595, nsg_732, \
                         nsg_733, nsg_734, osf0_486, osf1_486, osg_730, osg_732, osg_733, \
                         osg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_10 * nsg_732[k]
                    + f_3 * pc_x[k] * osg_732[k];

        t_1021[k] = f_10 * nsg_733[k]
                    + f_3 * pc_x[k] * osg_733[k];

        t_1022[k] = f_10 * nsg_734[k]
                    + f_3 * pc_x[k] * osg_734[k];

        t_1023[k] = f_19 * nsg_595[k]
                    + f_1 * osf0_486[k]
                    - f_2 * osf1_486[k]
                    + f_3 * pc_y[k] * osg_730[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pc_y, pc_z, nsg_580, nsg_597, nsg_598, \
                         osf0_488, osf0_489, osf1_488, osf1_489, osg_730, osg_732, \
                         osg_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_11 * nsg_580[k]
                    + f_3 * pc_z[k] * osg_730[k];

        t_1025[k] = f_19 * nsg_597[k]
                    + f_6 * osf0_488[k]
                    - f_7 * osf1_488[k]
                    + f_3 * pc_y[k] * osg_732[k];

        t_1026[k] = f_19 * nsg_598[k]
                    + f_4 * osf0_489[k]
                    - f_5 * osf1_489[k]
                    + f_3 * pc_y[k] * osg_733[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, pc_x, pc_y, pc_z, nsg_584, nsg_599, nsg_735, \
                         osf0_489, osf0_490, osf1_489, osf1_490, osg_734, \
                         osg_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_19 * nsg_599[k]
                    + f_3 * pc_y[k] * osg_734[k];

        t_1028[k] = f_11 * nsg_584[k]
                    + f_1 * osf0_489[k]
                    - f_2 * osf1_489[k]
                    + f_3 * pc_z[k] * osg_734[k];

        t_1029[k] = f_10 * nsg_735[k]
                    + f_1 * osf0_490[k]
                    - f_2 * osf1_490[k]
                    + f_3 * pc_x[k] * osg_735[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, pc_x, pc_y, pc_z, nsg_585, nsg_600, \
                         nsg_602, nsg_738, osf0_493, osf1_493, osg_735, osg_737, \
                         osg_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_20 * nsg_600[k]
                    + f_3 * pc_y[k] * osg_735[k];

        t_1031[k] = f_18 * nsg_585[k]
                    + f_3 * pc_z[k] * osg_735[k];

        t_1032[k] = f_10 * nsg_738[k]
                    + f_6 * osf0_493[k]
                    - f_7 * osf1_493[k]
                    + f_3 * pc_x[k] * osg_738[k];

        t_1033[k] = f_20 * nsg_602[k]
                    + f_3 * pc_y[k] * osg_737[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_z, nsg_588, nsg_740, nsg_741, \
                         osf0_495, osf0_496, osf1_495, osf1_496, osg_738, osg_740, \
                         osg_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_10 * nsg_740[k]
                    + f_6 * osf0_495[k]
                    - f_7 * osf1_495[k]
                    + f_3 * pc_x[k] * osg_740[k];

        t_1035[k] = f_10 * nsg_741[k]
                    + f_4 * osf0_496[k]
                    - f_5 * osf1_496[k]
                    + f_3 * pc_x[k] * osg_741[k];

        t_1036[k] = f_18 * nsg_588[k]
                    + f_3 * pc_z[k] * osg_738[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, pc_x, pc_y, nsg_605, nsg_744, \
                         nsg_745, nsg_746, osf0_499, osf1_499, osg_740, osg_744, osg_745, \
                         osg_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_20 * nsg_605[k]
                    + f_3 * pc_y[k] * osg_740[k];

        t_1038[k] = f_10 * nsg_744[k]
                    + f_4 * osf0_499[k]
                    - f_5 * osf1_499[k]
                    + f_3 * pc_x[k] * osg_744[k];

        t_1039[k] = f_10 * nsg_745[k]
                    + f_3 * pc_x[k] * osg_745[k];

        t_1040[k] = f_10 * nsg_746[k]
                    + f_3 * pc_x[k] * osg_746[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, t_1044, pc_x, pc_y, nsg_610, nsg_747, \
                         nsg_748, nsg_749, osf0_496, osf1_496, osg_745, osg_747, osg_748, \
                         osg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_10 * nsg_747[k]
                    + f_3 * pc_x[k] * osg_747[k];

        t_1042[k] = f_10 * nsg_748[k]
                    + f_3 * pc_x[k] * osg_748[k];

        t_1043[k] = f_10 * nsg_749[k]
                    + f_3 * pc_x[k] * osg_749[k];

        t_1044[k] = f_20 * nsg_610[k]
                    + f_1 * osf0_496[k]
                    - f_2 * osf1_496[k]
                    + f_3 * pc_y[k] * osg_745[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pc_y, pc_z, nsg_595, nsg_612, nsg_613, \
                         osf0_498, osf0_499, osf1_498, osf1_499, osg_745, osg_747, \
                         osg_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = f_18 * nsg_595[k]
                    + f_3 * pc_z[k] * osg_745[k];

        t_1046[k] = f_20 * nsg_612[k]
                    + f_6 * osf0_498[k]
                    - f_7 * osf1_498[k]
                    + f_3 * pc_y[k] * osg_747[k];

        t_1047[k] = f_20 * nsg_613[k]
                    + f_4 * osf0_499[k]
                    - f_5 * osf1_499[k]
                    + f_3 * pc_y[k] * osg_748[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pc_x, pc_y, pc_z, nsg_599, nsg_614, nsg_750, \
                         osf0_499, osf0_500, osf1_499, osf1_500, osg_749, \
                         osg_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_20 * nsg_614[k]
                    + f_3 * pc_y[k] * osg_749[k];

        t_1049[k] = f_18 * nsg_599[k]
                    + f_1 * osf0_499[k]
                    - f_2 * osf1_499[k]
                    + f_3 * pc_z[k] * osg_749[k];

        t_1050[k] = f_10 * nsg_750[k]
                    + f_1 * osf0_500[k]
                    - f_2 * osf1_500[k]
                    + f_3 * pc_x[k] * osg_750[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, pc_x, pc_y, pc_z, nsg_600, nsg_615, \
                         nsg_617, nsg_753, osf0_503, osf1_503, osg_750, osg_752, \
                         osg_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_18 * nsg_615[k]
                    + f_3 * pc_y[k] * osg_750[k];

        t_1052[k] = f_20 * nsg_600[k]
                    + f_3 * pc_z[k] * osg_750[k];

        t_1053[k] = f_10 * nsg_753[k]
                    + f_6 * osf0_503[k]
                    - f_7 * osf1_503[k]
                    + f_3 * pc_x[k] * osg_753[k];

        t_1054[k] = f_18 * nsg_617[k]
                    + f_3 * pc_y[k] * osg_752[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, nsg_603, nsg_755, nsg_756, \
                         osf0_505, osf0_506, osf1_505, osf1_506, osg_753, osg_755, \
                         osg_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_10 * nsg_755[k]
                    + f_6 * osf0_505[k]
                    - f_7 * osf1_505[k]
                    + f_3 * pc_x[k] * osg_755[k];

        t_1056[k] = f_10 * nsg_756[k]
                    + f_4 * osf0_506[k]
                    - f_5 * osf1_506[k]
                    + f_3 * pc_x[k] * osg_756[k];

        t_1057[k] = f_20 * nsg_603[k]
                    + f_3 * pc_z[k] * osg_753[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsh0,
                                                          const size_t nsg, const size_t nsh1,
                                                          const size_t osf0, const size_t osf1,
                                                          const size_t osg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_924 = buffer.data(nsh0 + 924);
    const auto *nsh0_927 = buffer.data(nsh0 + 927);
    const auto *nsh0_929 = buffer.data(nsh0 + 929);
    const auto *nsh0_930 = buffer.data(nsh0 + 930);
    const auto *nsh0_933 = buffer.data(nsh0 + 933);
    const auto *nsh0_944 = buffer.data(nsh0 + 944);
    const auto *nsh0_1155 = buffer.data(nsh0 + 1155);
    const auto *nsh0_1158 = buffer.data(nsh0 + 1158);
    const auto *nsh0_1161 = buffer.data(nsh0 + 1161);
    const auto *nsh0_1170 = buffer.data(nsh0 + 1170);
    const auto *nsh0_1172 = buffer.data(nsh0 + 1172);
    const auto *nsh0_1173 = buffer.data(nsh0 + 1173);

    const auto *nsg_610 = buffer.data(nsg + 610);
    const auto *nsg_614 = buffer.data(nsg + 614);
    const auto *nsg_615 = buffer.data(nsg + 615);
    const auto *nsg_618 = buffer.data(nsg + 618);
    const auto *nsg_620 = buffer.data(nsg + 620);
    const auto *nsg_625 = buffer.data(nsg + 625);
    const auto *nsg_627 = buffer.data(nsg + 627);
    const auto *nsg_628 = buffer.data(nsg + 628);
    const auto *nsg_629 = buffer.data(nsg + 629);
    const auto *nsg_630 = buffer.data(nsg + 630);
    const auto *nsg_632 = buffer.data(nsg + 632);
    const auto *nsg_633 = buffer.data(nsg + 633);
    const auto *nsg_635 = buffer.data(nsg + 635);
    const auto *nsg_640 = buffer.data(nsg + 640);
    const auto *nsg_642 = buffer.data(nsg + 642);
    const auto *nsg_643 = buffer.data(nsg + 643);
    const auto *nsg_644 = buffer.data(nsg + 644);
    const auto *nsg_645 = buffer.data(nsg + 645);
    const auto *nsg_647 = buffer.data(nsg + 647);
    const auto *nsg_648 = buffer.data(nsg + 648);
    const auto *nsg_650 = buffer.data(nsg + 650);
    const auto *nsg_655 = buffer.data(nsg + 655);
    const auto *nsg_657 = buffer.data(nsg + 657);
    const auto *nsg_658 = buffer.data(nsg + 658);
    const auto *nsg_659 = buffer.data(nsg + 659);
    const auto *nsg_660 = buffer.data(nsg + 660);
    const auto *nsg_661 = buffer.data(nsg + 661);
    const auto *nsg_662 = buffer.data(nsg + 662);
    const auto *nsg_663 = buffer.data(nsg + 663);
    const auto *nsg_665 = buffer.data(nsg + 665);
    const auto *nsg_670 = buffer.data(nsg + 670);
    const auto *nsg_672 = buffer.data(nsg + 672);
    const auto *nsg_673 = buffer.data(nsg + 673);
    const auto *nsg_674 = buffer.data(nsg + 674);
    const auto *nsg_675 = buffer.data(nsg + 675);
    const auto *nsg_680 = buffer.data(nsg + 680);
    const auto *nsg_689 = buffer.data(nsg + 689);
    const auto *nsg_759 = buffer.data(nsg + 759);
    const auto *nsg_760 = buffer.data(nsg + 760);
    const auto *nsg_761 = buffer.data(nsg + 761);
    const auto *nsg_762 = buffer.data(nsg + 762);
    const auto *nsg_763 = buffer.data(nsg + 763);
    const auto *nsg_764 = buffer.data(nsg + 764);
    const auto *nsg_765 = buffer.data(nsg + 765);
    const auto *nsg_768 = buffer.data(nsg + 768);
    const auto *nsg_770 = buffer.data(nsg + 770);
    const auto *nsg_771 = buffer.data(nsg + 771);
    const auto *nsg_774 = buffer.data(nsg + 774);
    const auto *nsg_775 = buffer.data(nsg + 775);
    const auto *nsg_776 = buffer.data(nsg + 776);
    const auto *nsg_777 = buffer.data(nsg + 777);
    const auto *nsg_778 = buffer.data(nsg + 778);
    const auto *nsg_779 = buffer.data(nsg + 779);
    const auto *nsg_780 = buffer.data(nsg + 780);
    const auto *nsg_783 = buffer.data(nsg + 783);
    const auto *nsg_785 = buffer.data(nsg + 785);
    const auto *nsg_786 = buffer.data(nsg + 786);
    const auto *nsg_789 = buffer.data(nsg + 789);
    const auto *nsg_790 = buffer.data(nsg + 790);
    const auto *nsg_791 = buffer.data(nsg + 791);
    const auto *nsg_792 = buffer.data(nsg + 792);
    const auto *nsg_793 = buffer.data(nsg + 793);
    const auto *nsg_794 = buffer.data(nsg + 794);
    const auto *nsg_805 = buffer.data(nsg + 805);
    const auto *nsg_806 = buffer.data(nsg + 806);
    const auto *nsg_807 = buffer.data(nsg + 807);
    const auto *nsg_808 = buffer.data(nsg + 808);
    const auto *nsg_809 = buffer.data(nsg + 809);
    const auto *nsg_810 = buffer.data(nsg + 810);
    const auto *nsg_815 = buffer.data(nsg + 815);
    const auto *nsg_819 = buffer.data(nsg + 819);
    const auto *nsg_820 = buffer.data(nsg + 820);
    const auto *nsg_821 = buffer.data(nsg + 821);
    const auto *nsg_822 = buffer.data(nsg + 822);
    const auto *nsg_824 = buffer.data(nsg + 824);
    const auto *nsg_825 = buffer.data(nsg + 825);
    const auto *nsg_828 = buffer.data(nsg + 828);
    const auto *nsg_831 = buffer.data(nsg + 831);
    const auto *nsg_835 = buffer.data(nsg + 835);
    const auto *nsg_837 = buffer.data(nsg + 837);
    const auto *nsg_838 = buffer.data(nsg + 838);
    const auto *nsg_839 = buffer.data(nsg + 839);

    const auto *nsh1_924 = buffer.data(nsh1 + 924);
    const auto *nsh1_927 = buffer.data(nsh1 + 927);
    const auto *nsh1_929 = buffer.data(nsh1 + 929);
    const auto *nsh1_930 = buffer.data(nsh1 + 930);
    const auto *nsh1_933 = buffer.data(nsh1 + 933);
    const auto *nsh1_944 = buffer.data(nsh1 + 944);
    const auto *nsh1_1155 = buffer.data(nsh1 + 1155);
    const auto *nsh1_1158 = buffer.data(nsh1 + 1158);
    const auto *nsh1_1161 = buffer.data(nsh1 + 1161);
    const auto *nsh1_1170 = buffer.data(nsh1 + 1170);
    const auto *nsh1_1172 = buffer.data(nsh1 + 1172);
    const auto *nsh1_1173 = buffer.data(nsh1 + 1173);

    const auto *osf0_506 = buffer.data(osf0 + 506);
    const auto *osf0_508 = buffer.data(osf0 + 508);
    const auto *osf0_509 = buffer.data(osf0 + 509);
    const auto *osf0_510 = buffer.data(osf0 + 510);
    const auto *osf0_513 = buffer.data(osf0 + 513);
    const auto *osf0_515 = buffer.data(osf0 + 515);
    const auto *osf0_516 = buffer.data(osf0 + 516);
    const auto *osf0_518 = buffer.data(osf0 + 518);
    const auto *osf0_519 = buffer.data(osf0 + 519);
    const auto *osf0_520 = buffer.data(osf0 + 520);
    const auto *osf0_523 = buffer.data(osf0 + 523);
    const auto *osf0_525 = buffer.data(osf0 + 525);
    const auto *osf0_526 = buffer.data(osf0 + 526);
    const auto *osf0_528 = buffer.data(osf0 + 528);
    const auto *osf0_529 = buffer.data(osf0 + 529);
    const auto *osf0_536 = buffer.data(osf0 + 536);
    const auto *osf0_538 = buffer.data(osf0 + 538);
    const auto *osf0_539 = buffer.data(osf0 + 539);
    const auto *osf0_540 = buffer.data(osf0 + 540);
    const auto *osf0_541 = buffer.data(osf0 + 541);
    const auto *osf0_542 = buffer.data(osf0 + 542);
    const auto *osf0_545 = buffer.data(osf0 + 545);
    const auto *osf0_546 = buffer.data(osf0 + 546);
    const auto *osf0_547 = buffer.data(osf0 + 547);
    const auto *osf0_548 = buffer.data(osf0 + 548);
    const auto *osf0_549 = buffer.data(osf0 + 549);
    const auto *osf0_550 = buffer.data(osf0 + 550);
    const auto *osf0_552 = buffer.data(osf0 + 552);

    const auto *osf1_506 = buffer.data(osf1 + 506);
    const auto *osf1_508 = buffer.data(osf1 + 508);
    const auto *osf1_509 = buffer.data(osf1 + 509);
    const auto *osf1_510 = buffer.data(osf1 + 510);
    const auto *osf1_513 = buffer.data(osf1 + 513);
    const auto *osf1_515 = buffer.data(osf1 + 515);
    const auto *osf1_516 = buffer.data(osf1 + 516);
    const auto *osf1_518 = buffer.data(osf1 + 518);
    const auto *osf1_519 = buffer.data(osf1 + 519);
    const auto *osf1_520 = buffer.data(osf1 + 520);
    const auto *osf1_523 = buffer.data(osf1 + 523);
    const auto *osf1_525 = buffer.data(osf1 + 525);
    const auto *osf1_526 = buffer.data(osf1 + 526);
    const auto *osf1_528 = buffer.data(osf1 + 528);
    const auto *osf1_529 = buffer.data(osf1 + 529);
    const auto *osf1_536 = buffer.data(osf1 + 536);
    const auto *osf1_538 = buffer.data(osf1 + 538);
    const auto *osf1_539 = buffer.data(osf1 + 539);
    const auto *osf1_540 = buffer.data(osf1 + 540);
    const auto *osf1_541 = buffer.data(osf1 + 541);
    const auto *osf1_542 = buffer.data(osf1 + 542);
    const auto *osf1_545 = buffer.data(osf1 + 545);
    const auto *osf1_546 = buffer.data(osf1 + 546);
    const auto *osf1_547 = buffer.data(osf1 + 547);
    const auto *osf1_548 = buffer.data(osf1 + 548);
    const auto *osf1_549 = buffer.data(osf1 + 549);
    const auto *osf1_550 = buffer.data(osf1 + 550);
    const auto *osf1_552 = buffer.data(osf1 + 552);

    const auto *osg_755 = buffer.data(osg + 755);
    const auto *osg_759 = buffer.data(osg + 759);
    const auto *osg_760 = buffer.data(osg + 760);
    const auto *osg_761 = buffer.data(osg + 761);
    const auto *osg_762 = buffer.data(osg + 762);
    const auto *osg_763 = buffer.data(osg + 763);
    const auto *osg_764 = buffer.data(osg + 764);
    const auto *osg_765 = buffer.data(osg + 765);
    const auto *osg_767 = buffer.data(osg + 767);
    const auto *osg_768 = buffer.data(osg + 768);
    const auto *osg_770 = buffer.data(osg + 770);
    const auto *osg_771 = buffer.data(osg + 771);
    const auto *osg_774 = buffer.data(osg + 774);
    const auto *osg_775 = buffer.data(osg + 775);
    const auto *osg_776 = buffer.data(osg + 776);
    const auto *osg_777 = buffer.data(osg + 777);
    const auto *osg_778 = buffer.data(osg + 778);
    const auto *osg_779 = buffer.data(osg + 779);
    const auto *osg_780 = buffer.data(osg + 780);
    const auto *osg_782 = buffer.data(osg + 782);
    const auto *osg_783 = buffer.data(osg + 783);
    const auto *osg_785 = buffer.data(osg + 785);
    const auto *osg_786 = buffer.data(osg + 786);
    const auto *osg_789 = buffer.data(osg + 789);
    const auto *osg_790 = buffer.data(osg + 790);
    const auto *osg_791 = buffer.data(osg + 791);
    const auto *osg_792 = buffer.data(osg + 792);
    const auto *osg_793 = buffer.data(osg + 793);
    const auto *osg_794 = buffer.data(osg + 794);
    const auto *osg_795 = buffer.data(osg + 795);
    const auto *osg_797 = buffer.data(osg + 797);
    const auto *osg_798 = buffer.data(osg + 798);
    const auto *osg_800 = buffer.data(osg + 800);
    const auto *osg_805 = buffer.data(osg + 805);
    const auto *osg_806 = buffer.data(osg + 806);
    const auto *osg_807 = buffer.data(osg + 807);
    const auto *osg_808 = buffer.data(osg + 808);
    const auto *osg_809 = buffer.data(osg + 809);
    const auto *osg_810 = buffer.data(osg + 810);
    const auto *osg_811 = buffer.data(osg + 811);
    const auto *osg_812 = buffer.data(osg + 812);
    const auto *osg_813 = buffer.data(osg + 813);
    const auto *osg_814 = buffer.data(osg + 814);
    const auto *osg_815 = buffer.data(osg + 815);
    const auto *osg_819 = buffer.data(osg + 819);
    const auto *osg_820 = buffer.data(osg + 820);
    const auto *osg_821 = buffer.data(osg + 821);
    const auto *osg_822 = buffer.data(osg + 822);
    const auto *osg_823 = buffer.data(osg + 823);
    const auto *osg_824 = buffer.data(osg + 824);
    const auto *osg_825 = buffer.data(osg + 825);
    const auto *osg_826 = buffer.data(osg + 826);
    const auto *osg_827 = buffer.data(osg + 827);
    const auto *osg_828 = buffer.data(osg + 828);
    const auto *osg_830 = buffer.data(osg + 830);
    const auto *osg_831 = buffer.data(osg + 831);
    const auto *osg_835 = buffer.data(osg + 835);
    const auto *osg_837 = buffer.data(osg + 837);
    const auto *osg_838 = buffer.data(osg + 838);
    const auto *osg_839 = buffer.data(osg + 839);

#pragma omp simd aligned(t_1058, t_1059, t_1060, t_1061, pc_x, pc_y, nsg_620, nsg_759, \
                         nsg_760, nsg_761, osf0_509, osf1_509, osg_755, osg_759, osg_760, \
                         osg_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_18 * nsg_620[k]
                    + f_3 * pc_y[k] * osg_755[k];

        t_1059[k] = f_10 * nsg_759[k]
                    + f_4 * osf0_509[k]
                    - f_5 * osf1_509[k]
                    + f_3 * pc_x[k] * osg_759[k];

        t_1060[k] = f_10 * nsg_760[k]
                    + f_3 * pc_x[k] * osg_760[k];

        t_1061[k] = f_10 * nsg_761[k]
                    + f_3 * pc_x[k] * osg_761[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, pc_x, pc_y, nsg_625, nsg_762, \
                         nsg_763, nsg_764, osf0_506, osf1_506, osg_760, osg_762, osg_763, \
                         osg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_10 * nsg_762[k]
                    + f_3 * pc_x[k] * osg_762[k];

        t_1063[k] = f_10 * nsg_763[k]
                    + f_3 * pc_x[k] * osg_763[k];

        t_1064[k] = f_10 * nsg_764[k]
                    + f_3 * pc_x[k] * osg_764[k];

        t_1065[k] = f_18 * nsg_625[k]
                    + f_1 * osf0_506[k]
                    - f_2 * osf1_506[k]
                    + f_3 * pc_y[k] * osg_760[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, pc_y, pc_z, nsg_610, nsg_627, nsg_628, \
                         osf0_508, osf0_509, osf1_508, osf1_509, osg_760, osg_762, \
                         osg_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_20 * nsg_610[k]
                    + f_3 * pc_z[k] * osg_760[k];

        t_1067[k] = f_18 * nsg_627[k]
                    + f_6 * osf0_508[k]
                    - f_7 * osf1_508[k]
                    + f_3 * pc_y[k] * osg_762[k];

        t_1068[k] = f_18 * nsg_628[k]
                    + f_4 * osf0_509[k]
                    - f_5 * osf1_509[k]
                    + f_3 * pc_y[k] * osg_763[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, pc_x, pc_y, pc_z, nsg_614, nsg_629, nsg_765, \
                         osf0_509, osf0_510, osf1_509, osf1_510, osg_764, \
                         osg_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_18 * nsg_629[k]
                    + f_3 * pc_y[k] * osg_764[k];

        t_1070[k] = f_20 * nsg_614[k]
                    + f_1 * osf0_509[k]
                    - f_2 * osf1_509[k]
                    + f_3 * pc_z[k] * osg_764[k];

        t_1071[k] = f_10 * nsg_765[k]
                    + f_1 * osf0_510[k]
                    - f_2 * osf1_510[k]
                    + f_3 * pc_x[k] * osg_765[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, pc_x, pc_y, pc_z, nsg_615, nsg_630, \
                         nsg_632, nsg_768, osf0_513, osf1_513, osg_765, osg_767, \
                         osg_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = f_11 * nsg_630[k]
                    + f_3 * pc_y[k] * osg_765[k];

        t_1073[k] = f_19 * nsg_615[k]
                    + f_3 * pc_z[k] * osg_765[k];

        t_1074[k] = f_10 * nsg_768[k]
                    + f_6 * osf0_513[k]
                    - f_7 * osf1_513[k]
                    + f_3 * pc_x[k] * osg_768[k];

        t_1075[k] = f_11 * nsg_632[k]
                    + f_3 * pc_y[k] * osg_767[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_z, nsg_618, nsg_770, nsg_771, \
                         osf0_515, osf0_516, osf1_515, osf1_516, osg_768, osg_770, \
                         osg_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_10 * nsg_770[k]
                    + f_6 * osf0_515[k]
                    - f_7 * osf1_515[k]
                    + f_3 * pc_x[k] * osg_770[k];

        t_1077[k] = f_10 * nsg_771[k]
                    + f_4 * osf0_516[k]
                    - f_5 * osf1_516[k]
                    + f_3 * pc_x[k] * osg_771[k];

        t_1078[k] = f_19 * nsg_618[k]
                    + f_3 * pc_z[k] * osg_768[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pc_x, pc_y, nsg_635, nsg_774, \
                         nsg_775, nsg_776, osf0_519, osf1_519, osg_770, osg_774, osg_775, \
                         osg_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_11 * nsg_635[k]
                    + f_3 * pc_y[k] * osg_770[k];

        t_1080[k] = f_10 * nsg_774[k]
                    + f_4 * osf0_519[k]
                    - f_5 * osf1_519[k]
                    + f_3 * pc_x[k] * osg_774[k];

        t_1081[k] = f_10 * nsg_775[k]
                    + f_3 * pc_x[k] * osg_775[k];

        t_1082[k] = f_10 * nsg_776[k]
                    + f_3 * pc_x[k] * osg_776[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, pc_x, pc_y, nsg_640, nsg_777, \
                         nsg_778, nsg_779, osf0_516, osf1_516, osg_775, osg_777, osg_778, \
                         osg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_10 * nsg_777[k]
                    + f_3 * pc_x[k] * osg_777[k];

        t_1084[k] = f_10 * nsg_778[k]
                    + f_3 * pc_x[k] * osg_778[k];

        t_1085[k] = f_10 * nsg_779[k]
                    + f_3 * pc_x[k] * osg_779[k];

        t_1086[k] = f_11 * nsg_640[k]
                    + f_1 * osf0_516[k]
                    - f_2 * osf1_516[k]
                    + f_3 * pc_y[k] * osg_775[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, pc_z, nsg_625, nsg_642, nsg_643, \
                         osf0_518, osf0_519, osf1_518, osf1_519, osg_775, osg_777, \
                         osg_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_19 * nsg_625[k]
                    + f_3 * pc_z[k] * osg_775[k];

        t_1088[k] = f_11 * nsg_642[k]
                    + f_6 * osf0_518[k]
                    - f_7 * osf1_518[k]
                    + f_3 * pc_y[k] * osg_777[k];

        t_1089[k] = f_11 * nsg_643[k]
                    + f_4 * osf0_519[k]
                    - f_5 * osf1_519[k]
                    + f_3 * pc_y[k] * osg_778[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, nsg_629, nsg_644, nsg_780, \
                         osf0_519, osf0_520, osf1_519, osf1_520, osg_779, \
                         osg_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_11 * nsg_644[k]
                    + f_3 * pc_y[k] * osg_779[k];

        t_1091[k] = f_19 * nsg_629[k]
                    + f_1 * osf0_519[k]
                    - f_2 * osf1_519[k]
                    + f_3 * pc_z[k] * osg_779[k];

        t_1092[k] = f_10 * nsg_780[k]
                    + f_1 * osf0_520[k]
                    - f_2 * osf1_520[k]
                    + f_3 * pc_x[k] * osg_780[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, nsg_630, nsg_645, \
                         nsg_647, nsg_783, osf0_523, osf1_523, osg_780, osg_782, \
                         osg_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_10 * nsg_645[k]
                    + f_3 * pc_y[k] * osg_780[k];

        t_1094[k] = f_17 * nsg_630[k]
                    + f_3 * pc_z[k] * osg_780[k];

        t_1095[k] = f_10 * nsg_783[k]
                    + f_6 * osf0_523[k]
                    - f_7 * osf1_523[k]
                    + f_3 * pc_x[k] * osg_783[k];

        t_1096[k] = f_10 * nsg_647[k]
                    + f_3 * pc_y[k] * osg_782[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, nsg_633, nsg_785, nsg_786, \
                         osf0_525, osf0_526, osf1_525, osf1_526, osg_783, osg_785, \
                         osg_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_10 * nsg_785[k]
                    + f_6 * osf0_525[k]
                    - f_7 * osf1_525[k]
                    + f_3 * pc_x[k] * osg_785[k];

        t_1098[k] = f_10 * nsg_786[k]
                    + f_4 * osf0_526[k]
                    - f_5 * osf1_526[k]
                    + f_3 * pc_x[k] * osg_786[k];

        t_1099[k] = f_17 * nsg_633[k]
                    + f_3 * pc_z[k] * osg_783[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, pc_y, nsg_650, nsg_789, \
                         nsg_790, nsg_791, osf0_529, osf1_529, osg_785, osg_789, osg_790, \
                         osg_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_10 * nsg_650[k]
                    + f_3 * pc_y[k] * osg_785[k];

        t_1101[k] = f_10 * nsg_789[k]
                    + f_4 * osf0_529[k]
                    - f_5 * osf1_529[k]
                    + f_3 * pc_x[k] * osg_789[k];

        t_1102[k] = f_10 * nsg_790[k]
                    + f_3 * pc_x[k] * osg_790[k];

        t_1103[k] = f_10 * nsg_791[k]
                    + f_3 * pc_x[k] * osg_791[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, pc_y, nsg_655, nsg_792, \
                         nsg_793, nsg_794, osf0_526, osf1_526, osg_790, osg_792, osg_793, \
                         osg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_10 * nsg_792[k]
                    + f_3 * pc_x[k] * osg_792[k];

        t_1105[k] = f_10 * nsg_793[k]
                    + f_3 * pc_x[k] * osg_793[k];

        t_1106[k] = f_10 * nsg_794[k]
                    + f_3 * pc_x[k] * osg_794[k];

        t_1107[k] = f_10 * nsg_655[k]
                    + f_1 * osf0_526[k]
                    - f_2 * osf1_526[k]
                    + f_3 * pc_y[k] * osg_790[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, nsg_640, nsg_657, nsg_658, \
                         osf0_528, osf0_529, osf1_528, osf1_529, osg_790, osg_792, \
                         osg_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * nsg_640[k]
                    + f_3 * pc_z[k] * osg_790[k];

        t_1109[k] = f_10 * nsg_657[k]
                    + f_6 * osf0_528[k]
                    - f_7 * osf1_528[k]
                    + f_3 * pc_y[k] * osg_792[k];

        t_1110[k] = f_10 * nsg_658[k]
                    + f_4 * osf0_529[k]
                    - f_5 * osf1_529[k]
                    + f_3 * pc_y[k] * osg_793[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, pa_y, pc_y, pc_z, nsh0_924, nsg_644, \
                         nsg_659, nsg_660, nsh1_924, osf0_529, osf1_529, osg_794, \
                         osg_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_10 * nsg_659[k]
                    + f_3 * pc_y[k] * osg_794[k];

        t_1112[k] = f_17 * nsg_644[k]
                    + f_1 * osf0_529[k]
                    - f_2 * osf1_529[k]
                    + f_3 * pc_z[k] * osg_794[k];

        t_1113[k] = pa_y[k] * nsh0_924[k]
                    - f_8 * pc_y[k] * nsh1_924[k];

        t_1114[k] = f_9 * nsg_660[k]
                    + f_3 * pc_y[k] * osg_795[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, pa_y, pc_y, pc_z, nsh0_927, nsh0_929, \
                         nsg_645, nsg_661, nsg_662, nsh1_927, nsh1_929, osg_795, \
                         osg_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_16 * nsg_645[k]
                    + f_3 * pc_z[k] * osg_795[k];

        t_1116[k] = pa_y[k] * nsh0_927[k]
                    + f_10 * nsg_661[k]
                    - f_8 * pc_y[k] * nsh1_927[k];

        t_1117[k] = f_9 * nsg_662[k]
                    + f_3 * pc_y[k] * osg_797[k];

        t_1118[k] = pa_y[k] * nsh0_929[k]
                    - f_8 * pc_y[k] * nsh1_929[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, t_1122, pa_y, pc_y, pc_z, nsh0_930, nsh0_933, \
                         nsg_648, nsg_663, nsg_665, nsh1_930, nsh1_933, osg_798, \
                         osg_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pa_y[k] * nsh0_930[k]
                    + f_11 * nsg_663[k]
                    - f_8 * pc_y[k] * nsh1_930[k];

        t_1120[k] = f_16 * nsg_648[k]
                    + f_3 * pc_z[k] * osg_798[k];

        t_1121[k] = f_9 * nsg_665[k]
                    + f_3 * pc_y[k] * osg_800[k];

        t_1122[k] = pa_y[k] * nsh0_933[k]
                    - f_8 * pc_y[k] * nsh1_933[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, t_1127, pc_x, nsg_805, nsg_806, \
                         nsg_807, nsg_808, nsg_809, osg_805, osg_806, osg_807, osg_808, \
                         osg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_10 * nsg_805[k]
                    + f_3 * pc_x[k] * osg_805[k];

        t_1124[k] = f_10 * nsg_806[k]
                    + f_3 * pc_x[k] * osg_806[k];

        t_1125[k] = f_10 * nsg_807[k]
                    + f_3 * pc_x[k] * osg_807[k];

        t_1126[k] = f_10 * nsg_808[k]
                    + f_3 * pc_x[k] * osg_808[k];

        t_1127[k] = f_10 * nsg_809[k]
                    + f_3 * pc_x[k] * osg_809[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pc_y, pc_z, nsg_655, nsg_670, nsg_672, \
                         osf0_536, osf0_538, osf1_536, osf1_538, osg_805, \
                         osg_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_9 * nsg_670[k]
                    + f_1 * osf0_536[k]
                    - f_2 * osf1_536[k]
                    + f_3 * pc_y[k] * osg_805[k];

        t_1129[k] = f_16 * nsg_655[k]
                    + f_3 * pc_z[k] * osg_805[k];

        t_1130[k] = f_9 * nsg_672[k]
                    + f_6 * osf0_538[k]
                    - f_7 * osf1_538[k]
                    + f_3 * pc_y[k] * osg_807[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pa_y, pc_y, nsh0_944, nsg_673, nsg_674, \
                         nsh1_944, osf0_539, osf1_539, osg_808, \
                         osg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_9 * nsg_673[k]
                    + f_4 * osf0_539[k]
                    - f_5 * osf1_539[k]
                    + f_3 * pc_y[k] * osg_808[k];

        t_1132[k] = f_9 * nsg_674[k]
                    + f_3 * pc_y[k] * osg_809[k];

        t_1133[k] = pa_y[k] * nsh0_944[k]
                    - f_8 * pc_y[k] * nsh1_944[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, nsg_660, \
                         nsg_810, osf0_540, osf1_540, osg_810, osg_811, \
                         osg_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_10 * nsg_810[k]
                    + f_1 * osf0_540[k]
                    - f_2 * osf1_540[k]
                    + f_3 * pc_x[k] * osg_810[k];

        t_1135[k] = f_3 * pc_y[k] * osg_810[k];

        t_1136[k] = f_15 * nsg_660[k]
                    + f_3 * pc_z[k] * osg_810[k];

        t_1137[k] = f_4 * osf0_540[k]
                    - f_5 * osf1_540[k]
                    + f_3 * pc_y[k] * osg_811[k];

        t_1138[k] = f_3 * pc_y[k] * osg_812[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, pc_x, pc_y, nsg_815, osf0_541, \
                         osf0_542, osf0_545, osf1_541, osf1_542, osf1_545, osg_813, osg_814, \
                         osg_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_10 * nsg_815[k]
                    + f_6 * osf0_545[k]
                    - f_7 * osf1_545[k]
                    + f_3 * pc_x[k] * osg_815[k];

        t_1140[k] = f_6 * osf0_541[k]
                    - f_7 * osf1_541[k]
                    + f_3 * pc_y[k] * osg_813[k];

        t_1141[k] = f_4 * osf0_542[k]
                    - f_5 * osf1_542[k]
                    + f_3 * pc_y[k] * osg_814[k];

        t_1142[k] = f_3 * pc_y[k] * osg_815[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, t_1146, pc_x, nsg_819, nsg_820, nsg_821, \
                         nsg_822, osf0_549, osf1_549, osg_819, osg_820, osg_821, \
                         osg_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_10 * nsg_819[k]
                    + f_4 * osf0_549[k]
                    - f_5 * osf1_549[k]
                    + f_3 * pc_x[k] * osg_819[k];

        t_1144[k] = f_10 * nsg_820[k]
                    + f_3 * pc_x[k] * osg_820[k];

        t_1145[k] = f_10 * nsg_821[k]
                    + f_3 * pc_x[k] * osg_821[k];

        t_1146[k] = f_10 * nsg_822[k]
                    + f_3 * pc_x[k] * osg_822[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, pc_x, pc_y, nsg_824, osf0_546, \
                         osf0_547, osf1_546, osf1_547, osg_819, osg_820, osg_821, \
                         osg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_3 * pc_y[k] * osg_819[k];

        t_1148[k] = f_10 * nsg_824[k]
                    + f_3 * pc_x[k] * osg_824[k];

        t_1149[k] = f_1 * osf0_546[k]
                    - f_2 * osf1_546[k]
                    + f_3 * pc_y[k] * osg_820[k];

        t_1150[k] = f_13 * osf0_547[k]
                    - f_14 * osf1_547[k]
                    + f_3 * pc_y[k] * osg_821[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pc_y, pc_z, nsg_674, osf0_548, \
                         osf0_549, osf1_548, osf1_549, osg_822, osg_823, \
                         osg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_6 * osf0_548[k]
                    - f_7 * osf1_548[k]
                    + f_3 * pc_y[k] * osg_822[k];

        t_1152[k] = f_4 * osf0_549[k]
                    - f_5 * osf1_549[k]
                    + f_3 * pc_y[k] * osg_823[k];

        t_1153[k] = f_3 * pc_y[k] * osg_824[k];

        t_1154[k] = f_15 * nsg_674[k]
                    + f_1 * osf0_549[k]
                    - f_2 * osf1_549[k]
                    + f_3 * pc_z[k] * osg_824[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, pa_x, pc_x, pc_y, pc_z, nsh0_1155, \
                         nsh0_1158, nsg_675, nsg_825, nsg_828, nsh1_1155, nsh1_1158, \
                         osg_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = pa_x[k] * nsh0_1155[k]
                    + f_20 * nsg_825[k]
                    - f_8 * pc_x[k] * nsh1_1155[k];

        t_1156[k] = f_12 * nsg_675[k]
                    + f_3 * pc_y[k] * osg_825[k];

        t_1157[k] = f_3 * pc_z[k] * osg_825[k];

        t_1158[k] = pa_x[k] * nsh0_1158[k]
                    + f_11 * nsg_828[k]
                    - f_8 * pc_x[k] * nsh1_1158[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, t_1162, pa_x, pc_x, pc_z, nsh0_1161, nsg_831, \
                         nsh1_1161, osf0_550, osf1_550, osg_826, osg_827, \
                         osg_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_3 * pc_z[k] * osg_826[k];

        t_1160[k] = f_4 * osf0_550[k]
                    - f_5 * osf1_550[k]
                    + f_3 * pc_z[k] * osg_827[k];

        t_1161[k] = pa_x[k] * nsh0_1161[k]
                    + f_10 * nsg_831[k]
                    - f_8 * pc_x[k] * nsh1_1161[k];

        t_1162[k] = f_3 * pc_z[k] * osg_828[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, t_1166, pc_x, pc_y, pc_z, nsg_680, nsg_835, \
                         osf0_552, osf1_552, osg_830, osg_831, \
                         osg_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_12 * nsg_680[k]
                    + f_3 * pc_y[k] * osg_830[k];

        t_1164[k] = f_6 * osf0_552[k]
                    - f_7 * osf1_552[k]
                    + f_3 * pc_z[k] * osg_830[k];

        t_1165[k] = f_9 * nsg_835[k]
                    + f_3 * pc_x[k] * osg_835[k];

        t_1166[k] = f_3 * pc_z[k] * osg_831[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, pa_x, pc_x, nsh0_1170, nsg_837, \
                         nsg_838, nsg_839, nsh1_1170, osg_837, osg_838, \
                         osg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_9 * nsg_837[k]
                    + f_3 * pc_x[k] * osg_837[k];

        t_1168[k] = f_9 * nsg_838[k]
                    + f_3 * pc_x[k] * osg_838[k];

        t_1169[k] = f_9 * nsg_839[k]
                    + f_3 * pc_x[k] * osg_839[k];

        t_1170[k] = pa_x[k] * nsh0_1170[k]
                    - f_8 * pc_x[k] * nsh1_1170[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, pa_x, pc_x, pc_y, pc_z, nsh0_1172, \
                         nsh0_1173, nsg_689, nsh1_1172, nsh1_1173, osg_835, \
                         osg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_3 * pc_z[k] * osg_835[k];

        t_1172[k] = pa_x[k] * nsh0_1172[k]
                    - f_8 * pc_x[k] * nsh1_1172[k];

        t_1173[k] = pa_x[k] * nsh0_1173[k]
                    - f_8 * pc_x[k] * nsh1_1173[k];

        t_1174[k] = f_12 * nsg_689[k]
                    + f_3 * pc_y[k] * osg_839[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsh0,
                                                           const size_t nsg, const size_t nsh1,
                                                           const size_t osg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_945 = buffer.data(nsh0 + 945);
    const auto *nsh0_948 = buffer.data(nsh0 + 948);
    const auto *nsh0_951 = buffer.data(nsh0 + 951);
    const auto *nsh0_1175 = buffer.data(nsh0 + 1175);
    const auto *nsh0_1181 = buffer.data(nsh0 + 1181);
    const auto *nsh0_1185 = buffer.data(nsh0 + 1185);
    const auto *nsh0_1191 = buffer.data(nsh0 + 1191);
    const auto *nsh0_1193 = buffer.data(nsh0 + 1193);
    const auto *nsh0_1194 = buffer.data(nsh0 + 1194);
    const auto *nsh0_1196 = buffer.data(nsh0 + 1196);
    const auto *nsh0_1197 = buffer.data(nsh0 + 1197);
    const auto *nsh0_1200 = buffer.data(nsh0 + 1200);
    const auto *nsh0_1202 = buffer.data(nsh0 + 1202);
    const auto *nsh0_1203 = buffer.data(nsh0 + 1203);
    const auto *nsh0_1206 = buffer.data(nsh0 + 1206);
    const auto *nsh0_1212 = buffer.data(nsh0 + 1212);
    const auto *nsh0_1214 = buffer.data(nsh0 + 1214);
    const auto *nsh0_1215 = buffer.data(nsh0 + 1215);
    const auto *nsh0_1217 = buffer.data(nsh0 + 1217);
    const auto *nsh0_1218 = buffer.data(nsh0 + 1218);
    const auto *nsh0_1221 = buffer.data(nsh0 + 1221);
    const auto *nsh0_1223 = buffer.data(nsh0 + 1223);
    const auto *nsh0_1224 = buffer.data(nsh0 + 1224);
    const auto *nsh0_1227 = buffer.data(nsh0 + 1227);
    const auto *nsh0_1233 = buffer.data(nsh0 + 1233);
    const auto *nsh0_1235 = buffer.data(nsh0 + 1235);
    const auto *nsh0_1236 = buffer.data(nsh0 + 1236);
    const auto *nsh0_1238 = buffer.data(nsh0 + 1238);
    const auto *nsh0_1239 = buffer.data(nsh0 + 1239);
    const auto *nsh0_1242 = buffer.data(nsh0 + 1242);
    const auto *nsh0_1244 = buffer.data(nsh0 + 1244);
    const auto *nsh0_1245 = buffer.data(nsh0 + 1245);
    const auto *nsh0_1248 = buffer.data(nsh0 + 1248);
    const auto *nsh0_1254 = buffer.data(nsh0 + 1254);
    const auto *nsh0_1256 = buffer.data(nsh0 + 1256);
    const auto *nsh0_1257 = buffer.data(nsh0 + 1257);
    const auto *nsh0_1259 = buffer.data(nsh0 + 1259);
    const auto *nsh0_1260 = buffer.data(nsh0 + 1260);
    const auto *nsh0_1263 = buffer.data(nsh0 + 1263);
    const auto *nsh0_1265 = buffer.data(nsh0 + 1265);
    const auto *nsh0_1266 = buffer.data(nsh0 + 1266);
    const auto *nsh0_1269 = buffer.data(nsh0 + 1269);
    const auto *nsh0_1275 = buffer.data(nsh0 + 1275);
    const auto *nsh0_1277 = buffer.data(nsh0 + 1277);
    const auto *nsh0_1278 = buffer.data(nsh0 + 1278);
    const auto *nsh0_1280 = buffer.data(nsh0 + 1280);
    const auto *nsh0_1281 = buffer.data(nsh0 + 1281);
    const auto *nsh0_1284 = buffer.data(nsh0 + 1284);
    const auto *nsh0_1286 = buffer.data(nsh0 + 1286);
    const auto *nsh0_1287 = buffer.data(nsh0 + 1287);
    const auto *nsh0_1290 = buffer.data(nsh0 + 1290);

    const auto *nsg_675 = buffer.data(nsg + 675);
    const auto *nsg_678 = buffer.data(nsg + 678);
    const auto *nsg_685 = buffer.data(nsg + 685);
    const auto *nsg_690 = buffer.data(nsg + 690);
    const auto *nsg_692 = buffer.data(nsg + 692);
    const auto *nsg_693 = buffer.data(nsg + 693);
    const auto *nsg_695 = buffer.data(nsg + 695);
    const auto *nsg_700 = buffer.data(nsg + 700);
    const auto *nsg_704 = buffer.data(nsg + 704);
    const auto *nsg_705 = buffer.data(nsg + 705);
    const auto *nsg_707 = buffer.data(nsg + 707);
    const auto *nsg_708 = buffer.data(nsg + 708);
    const auto *nsg_710 = buffer.data(nsg + 710);
    const auto *nsg_715 = buffer.data(nsg + 715);
    const auto *nsg_719 = buffer.data(nsg + 719);
    const auto *nsg_720 = buffer.data(nsg + 720);
    const auto *nsg_722 = buffer.data(nsg + 722);
    const auto *nsg_723 = buffer.data(nsg + 723);
    const auto *nsg_725 = buffer.data(nsg + 725);
    const auto *nsg_730 = buffer.data(nsg + 730);
    const auto *nsg_734 = buffer.data(nsg + 734);
    const auto *nsg_735 = buffer.data(nsg + 735);
    const auto *nsg_737 = buffer.data(nsg + 737);
    const auto *nsg_738 = buffer.data(nsg + 738);
    const auto *nsg_740 = buffer.data(nsg + 740);
    const auto *nsg_745 = buffer.data(nsg + 745);
    const auto *nsg_749 = buffer.data(nsg + 749);
    const auto *nsg_750 = buffer.data(nsg + 750);
    const auto *nsg_752 = buffer.data(nsg + 752);
    const auto *nsg_753 = buffer.data(nsg + 753);
    const auto *nsg_755 = buffer.data(nsg + 755);
    const auto *nsg_764 = buffer.data(nsg + 764);
    const auto *nsg_765 = buffer.data(nsg + 765);
    const auto *nsg_767 = buffer.data(nsg + 767);
    const auto *nsg_770 = buffer.data(nsg + 770);
    const auto *nsg_845 = buffer.data(nsg + 845);
    const auto *nsg_849 = buffer.data(nsg + 849);
    const auto *nsg_850 = buffer.data(nsg + 850);
    const auto *nsg_851 = buffer.data(nsg + 851);
    const auto *nsg_852 = buffer.data(nsg + 852);
    const auto *nsg_853 = buffer.data(nsg + 853);
    const auto *nsg_854 = buffer.data(nsg + 854);
    const auto *nsg_855 = buffer.data(nsg + 855);
    const auto *nsg_858 = buffer.data(nsg + 858);
    const auto *nsg_860 = buffer.data(nsg + 860);
    const auto *nsg_861 = buffer.data(nsg + 861);
    const auto *nsg_864 = buffer.data(nsg + 864);
    const auto *nsg_865 = buffer.data(nsg + 865);
    const auto *nsg_866 = buffer.data(nsg + 866);
    const auto *nsg_867 = buffer.data(nsg + 867);
    const auto *nsg_868 = buffer.data(nsg + 868);
    const auto *nsg_869 = buffer.data(nsg + 869);
    const auto *nsg_870 = buffer.data(nsg + 870);
    const auto *nsg_873 = buffer.data(nsg + 873);
    const auto *nsg_875 = buffer.data(nsg + 875);
    const auto *nsg_876 = buffer.data(nsg + 876);
    const auto *nsg_879 = buffer.data(nsg + 879);
    const auto *nsg_880 = buffer.data(nsg + 880);
    const auto *nsg_881 = buffer.data(nsg + 881);
    const auto *nsg_882 = buffer.data(nsg + 882);
    const auto *nsg_883 = buffer.data(nsg + 883);
    const auto *nsg_884 = buffer.data(nsg + 884);
    const auto *nsg_885 = buffer.data(nsg + 885);
    const auto *nsg_888 = buffer.data(nsg + 888);
    const auto *nsg_890 = buffer.data(nsg + 890);
    const auto *nsg_891 = buffer.data(nsg + 891);
    const auto *nsg_894 = buffer.data(nsg + 894);
    const auto *nsg_895 = buffer.data(nsg + 895);
    const auto *nsg_896 = buffer.data(nsg + 896);
    const auto *nsg_897 = buffer.data(nsg + 897);
    const auto *nsg_898 = buffer.data(nsg + 898);
    const auto *nsg_899 = buffer.data(nsg + 899);
    const auto *nsg_900 = buffer.data(nsg + 900);
    const auto *nsg_903 = buffer.data(nsg + 903);
    const auto *nsg_905 = buffer.data(nsg + 905);
    const auto *nsg_906 = buffer.data(nsg + 906);
    const auto *nsg_909 = buffer.data(nsg + 909);
    const auto *nsg_910 = buffer.data(nsg + 910);
    const auto *nsg_911 = buffer.data(nsg + 911);
    const auto *nsg_912 = buffer.data(nsg + 912);
    const auto *nsg_913 = buffer.data(nsg + 913);
    const auto *nsg_914 = buffer.data(nsg + 914);
    const auto *nsg_915 = buffer.data(nsg + 915);
    const auto *nsg_918 = buffer.data(nsg + 918);
    const auto *nsg_920 = buffer.data(nsg + 920);
    const auto *nsg_921 = buffer.data(nsg + 921);
    const auto *nsg_924 = buffer.data(nsg + 924);
    const auto *nsg_925 = buffer.data(nsg + 925);
    const auto *nsg_926 = buffer.data(nsg + 926);
    const auto *nsg_927 = buffer.data(nsg + 927);

    const auto *nsh1_945 = buffer.data(nsh1 + 945);
    const auto *nsh1_948 = buffer.data(nsh1 + 948);
    const auto *nsh1_951 = buffer.data(nsh1 + 951);
    const auto *nsh1_1175 = buffer.data(nsh1 + 1175);
    const auto *nsh1_1181 = buffer.data(nsh1 + 1181);
    const auto *nsh1_1185 = buffer.data(nsh1 + 1185);
    const auto *nsh1_1191 = buffer.data(nsh1 + 1191);
    const auto *nsh1_1193 = buffer.data(nsh1 + 1193);
    const auto *nsh1_1194 = buffer.data(nsh1 + 1194);
    const auto *nsh1_1196 = buffer.data(nsh1 + 1196);
    const auto *nsh1_1197 = buffer.data(nsh1 + 1197);
    const auto *nsh1_1200 = buffer.data(nsh1 + 1200);
    const auto *nsh1_1202 = buffer.data(nsh1 + 1202);
    const auto *nsh1_1203 = buffer.data(nsh1 + 1203);
    const auto *nsh1_1206 = buffer.data(nsh1 + 1206);
    const auto *nsh1_1212 = buffer.data(nsh1 + 1212);
    const auto *nsh1_1214 = buffer.data(nsh1 + 1214);
    const auto *nsh1_1215 = buffer.data(nsh1 + 1215);
    const auto *nsh1_1217 = buffer.data(nsh1 + 1217);
    const auto *nsh1_1218 = buffer.data(nsh1 + 1218);
    const auto *nsh1_1221 = buffer.data(nsh1 + 1221);
    const auto *nsh1_1223 = buffer.data(nsh1 + 1223);
    const auto *nsh1_1224 = buffer.data(nsh1 + 1224);
    const auto *nsh1_1227 = buffer.data(nsh1 + 1227);
    const auto *nsh1_1233 = buffer.data(nsh1 + 1233);
    const auto *nsh1_1235 = buffer.data(nsh1 + 1235);
    const auto *nsh1_1236 = buffer.data(nsh1 + 1236);
    const auto *nsh1_1238 = buffer.data(nsh1 + 1238);
    const auto *nsh1_1239 = buffer.data(nsh1 + 1239);
    const auto *nsh1_1242 = buffer.data(nsh1 + 1242);
    const auto *nsh1_1244 = buffer.data(nsh1 + 1244);
    const auto *nsh1_1245 = buffer.data(nsh1 + 1245);
    const auto *nsh1_1248 = buffer.data(nsh1 + 1248);
    const auto *nsh1_1254 = buffer.data(nsh1 + 1254);
    const auto *nsh1_1256 = buffer.data(nsh1 + 1256);
    const auto *nsh1_1257 = buffer.data(nsh1 + 1257);
    const auto *nsh1_1259 = buffer.data(nsh1 + 1259);
    const auto *nsh1_1260 = buffer.data(nsh1 + 1260);
    const auto *nsh1_1263 = buffer.data(nsh1 + 1263);
    const auto *nsh1_1265 = buffer.data(nsh1 + 1265);
    const auto *nsh1_1266 = buffer.data(nsh1 + 1266);
    const auto *nsh1_1269 = buffer.data(nsh1 + 1269);
    const auto *nsh1_1275 = buffer.data(nsh1 + 1275);
    const auto *nsh1_1277 = buffer.data(nsh1 + 1277);
    const auto *nsh1_1278 = buffer.data(nsh1 + 1278);
    const auto *nsh1_1280 = buffer.data(nsh1 + 1280);
    const auto *nsh1_1281 = buffer.data(nsh1 + 1281);
    const auto *nsh1_1284 = buffer.data(nsh1 + 1284);
    const auto *nsh1_1286 = buffer.data(nsh1 + 1286);
    const auto *nsh1_1287 = buffer.data(nsh1 + 1287);
    const auto *nsh1_1290 = buffer.data(nsh1 + 1290);

    const auto *osg_840 = buffer.data(osg + 840);
    const auto *osg_842 = buffer.data(osg + 842);
    const auto *osg_843 = buffer.data(osg + 843);
    const auto *osg_845 = buffer.data(osg + 845);
    const auto *osg_850 = buffer.data(osg + 850);
    const auto *osg_851 = buffer.data(osg + 851);
    const auto *osg_852 = buffer.data(osg + 852);
    const auto *osg_853 = buffer.data(osg + 853);
    const auto *osg_854 = buffer.data(osg + 854);
    const auto *osg_855 = buffer.data(osg + 855);
    const auto *osg_857 = buffer.data(osg + 857);
    const auto *osg_858 = buffer.data(osg + 858);
    const auto *osg_860 = buffer.data(osg + 860);
    const auto *osg_865 = buffer.data(osg + 865);
    const auto *osg_866 = buffer.data(osg + 866);
    const auto *osg_867 = buffer.data(osg + 867);
    const auto *osg_868 = buffer.data(osg + 868);
    const auto *osg_869 = buffer.data(osg + 869);
    const auto *osg_870 = buffer.data(osg + 870);
    const auto *osg_872 = buffer.data(osg + 872);
    const auto *osg_873 = buffer.data(osg + 873);
    const auto *osg_875 = buffer.data(osg + 875);
    const auto *osg_880 = buffer.data(osg + 880);
    const auto *osg_881 = buffer.data(osg + 881);
    const auto *osg_882 = buffer.data(osg + 882);
    const auto *osg_883 = buffer.data(osg + 883);
    const auto *osg_884 = buffer.data(osg + 884);
    const auto *osg_885 = buffer.data(osg + 885);
    const auto *osg_887 = buffer.data(osg + 887);
    const auto *osg_888 = buffer.data(osg + 888);
    const auto *osg_890 = buffer.data(osg + 890);
    const auto *osg_895 = buffer.data(osg + 895);
    const auto *osg_896 = buffer.data(osg + 896);
    const auto *osg_897 = buffer.data(osg + 897);
    const auto *osg_898 = buffer.data(osg + 898);
    const auto *osg_899 = buffer.data(osg + 899);
    const auto *osg_900 = buffer.data(osg + 900);
    const auto *osg_902 = buffer.data(osg + 902);
    const auto *osg_903 = buffer.data(osg + 903);
    const auto *osg_905 = buffer.data(osg + 905);
    const auto *osg_910 = buffer.data(osg + 910);
    const auto *osg_911 = buffer.data(osg + 911);
    const auto *osg_912 = buffer.data(osg + 912);
    const auto *osg_913 = buffer.data(osg + 913);
    const auto *osg_914 = buffer.data(osg + 914);
    const auto *osg_915 = buffer.data(osg + 915);
    const auto *osg_917 = buffer.data(osg + 917);
    const auto *osg_918 = buffer.data(osg + 918);
    const auto *osg_920 = buffer.data(osg + 920);
    const auto *osg_925 = buffer.data(osg + 925);
    const auto *osg_926 = buffer.data(osg + 926);
    const auto *osg_927 = buffer.data(osg + 927);

#pragma omp simd aligned(t_1175, t_1176, t_1177, t_1178, pa_x, pa_z, pc_x, pc_y, pc_z, \
                         nsh0_945, nsh0_1175, nsg_675, nsg_690, nsh1_945, nsh1_1175, \
                         osg_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = pa_x[k] * nsh0_1175[k]
                    - f_8 * pc_x[k] * nsh1_1175[k];

        t_1176[k] = pa_z[k] * nsh0_945[k]
                    - f_8 * pc_z[k] * nsh1_945[k];

        t_1177[k] = f_15 * nsg_690[k]
                    + f_3 * pc_y[k] * osg_840[k];

        t_1178[k] = f_9 * nsg_675[k]
                    + f_3 * pc_z[k] * osg_840[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pa_x, pa_z, pc_x, pc_y, pc_z, nsh0_948, \
                         nsh0_1181, nsg_692, nsg_845, nsh1_948, nsh1_1181, \
                         osg_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = pa_z[k] * nsh0_948[k]
                    - f_8 * pc_z[k] * nsh1_948[k];

        t_1180[k] = f_15 * nsg_692[k]
                    + f_3 * pc_y[k] * osg_842[k];

        t_1181[k] = pa_x[k] * nsh0_1181[k]
                    + f_11 * nsg_845[k]
                    - f_8 * pc_x[k] * nsh1_1181[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pa_z, pc_y, pc_z, nsh0_951, nsg_678, nsg_695, \
                         nsh1_951, osg_843, osg_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = pa_z[k] * nsh0_951[k]
                    - f_8 * pc_z[k] * nsh1_951[k];

        t_1183[k] = f_9 * nsg_678[k]
                    + f_3 * pc_z[k] * osg_843[k];

        t_1184[k] = f_15 * nsg_695[k]
                    + f_3 * pc_y[k] * osg_845[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, pa_x, pc_x, nsh0_1185, nsg_849, \
                         nsg_850, nsg_851, nsg_852, nsh1_1185, osg_850, osg_851, \
                         osg_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pa_x[k] * nsh0_1185[k]
                    + f_10 * nsg_849[k]
                    - f_8 * pc_x[k] * nsh1_1185[k];

        t_1186[k] = f_9 * nsg_850[k]
                    + f_3 * pc_x[k] * osg_850[k];

        t_1187[k] = f_9 * nsg_851[k]
                    + f_3 * pc_x[k] * osg_851[k];

        t_1188[k] = f_9 * nsg_852[k]
                    + f_3 * pc_x[k] * osg_852[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_x, pc_x, pc_z, nsh0_1191, nsg_685, \
                         nsg_853, nsg_854, nsh1_1191, osg_850, osg_853, \
                         osg_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_9 * nsg_853[k]
                    + f_3 * pc_x[k] * osg_853[k];

        t_1190[k] = f_9 * nsg_854[k]
                    + f_3 * pc_x[k] * osg_854[k];

        t_1191[k] = pa_x[k] * nsh0_1191[k]
                    - f_8 * pc_x[k] * nsh1_1191[k];

        t_1192[k] = f_9 * nsg_685[k]
                    + f_3 * pc_z[k] * osg_850[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, pa_x, pc_x, pc_y, nsh0_1193, \
                         nsh0_1194, nsh0_1196, nsg_704, nsh1_1193, nsh1_1194, nsh1_1196, \
                         osg_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = pa_x[k] * nsh0_1193[k]
                    - f_8 * pc_x[k] * nsh1_1193[k];

        t_1194[k] = pa_x[k] * nsh0_1194[k]
                    - f_8 * pc_x[k] * nsh1_1194[k];

        t_1195[k] = f_15 * nsg_704[k]
                    + f_3 * pc_y[k] * osg_854[k];

        t_1196[k] = pa_x[k] * nsh0_1196[k]
                    - f_8 * pc_x[k] * nsh1_1196[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pa_x, pc_x, pc_y, pc_z, nsh0_1197, nsg_690, \
                         nsg_705, nsg_855, nsh1_1197, osg_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = pa_x[k] * nsh0_1197[k]
                    + f_20 * nsg_855[k]
                    - f_8 * pc_x[k] * nsh1_1197[k];

        t_1198[k] = f_16 * nsg_705[k]
                    + f_3 * pc_y[k] * osg_855[k];

        t_1199[k] = f_10 * nsg_690[k]
                    + f_3 * pc_z[k] * osg_855[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pa_x, pc_x, pc_y, nsh0_1200, nsh0_1202, \
                         nsg_707, nsg_858, nsg_860, nsh1_1200, nsh1_1202, \
                         osg_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = pa_x[k] * nsh0_1200[k]
                    + f_11 * nsg_858[k]
                    - f_8 * pc_x[k] * nsh1_1200[k];

        t_1201[k] = f_16 * nsg_707[k]
                    + f_3 * pc_y[k] * osg_857[k];

        t_1202[k] = pa_x[k] * nsh0_1202[k]
                    + f_11 * nsg_860[k]
                    - f_8 * pc_x[k] * nsh1_1202[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, pa_x, pc_x, pc_y, pc_z, nsh0_1203, nsg_693, \
                         nsg_710, nsg_861, nsh1_1203, osg_858, \
                         osg_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = pa_x[k] * nsh0_1203[k]
                    + f_10 * nsg_861[k]
                    - f_8 * pc_x[k] * nsh1_1203[k];

        t_1204[k] = f_10 * nsg_693[k]
                    + f_3 * pc_z[k] * osg_858[k];

        t_1205[k] = f_16 * nsg_710[k]
                    + f_3 * pc_y[k] * osg_860[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pa_x, pc_x, nsh0_1206, nsg_864, \
                         nsg_865, nsg_866, nsg_867, nsh1_1206, osg_865, osg_866, \
                         osg_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = pa_x[k] * nsh0_1206[k]
                    + f_10 * nsg_864[k]
                    - f_8 * pc_x[k] * nsh1_1206[k];

        t_1207[k] = f_9 * nsg_865[k]
                    + f_3 * pc_x[k] * osg_865[k];

        t_1208[k] = f_9 * nsg_866[k]
                    + f_3 * pc_x[k] * osg_866[k];

        t_1209[k] = f_9 * nsg_867[k]
                    + f_3 * pc_x[k] * osg_867[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pa_x, pc_x, pc_z, nsh0_1212, nsg_700, \
                         nsg_868, nsg_869, nsh1_1212, osg_865, osg_868, \
                         osg_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_9 * nsg_868[k]
                    + f_3 * pc_x[k] * osg_868[k];

        t_1211[k] = f_9 * nsg_869[k]
                    + f_3 * pc_x[k] * osg_869[k];

        t_1212[k] = pa_x[k] * nsh0_1212[k]
                    - f_8 * pc_x[k] * nsh1_1212[k];

        t_1213[k] = f_10 * nsg_700[k]
                    + f_3 * pc_z[k] * osg_865[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pa_x, pc_x, pc_y, nsh0_1214, \
                         nsh0_1215, nsh0_1217, nsg_719, nsh1_1214, nsh1_1215, nsh1_1217, \
                         osg_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_x[k] * nsh0_1214[k]
                    - f_8 * pc_x[k] * nsh1_1214[k];

        t_1215[k] = pa_x[k] * nsh0_1215[k]
                    - f_8 * pc_x[k] * nsh1_1215[k];

        t_1216[k] = f_16 * nsg_719[k]
                    + f_3 * pc_y[k] * osg_869[k];

        t_1217[k] = pa_x[k] * nsh0_1217[k]
                    - f_8 * pc_x[k] * nsh1_1217[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pa_x, pc_x, pc_y, pc_z, nsh0_1218, nsg_705, \
                         nsg_720, nsg_870, nsh1_1218, osg_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pa_x[k] * nsh0_1218[k]
                    + f_20 * nsg_870[k]
                    - f_8 * pc_x[k] * nsh1_1218[k];

        t_1219[k] = f_17 * nsg_720[k]
                    + f_3 * pc_y[k] * osg_870[k];

        t_1220[k] = f_11 * nsg_705[k]
                    + f_3 * pc_z[k] * osg_870[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pa_x, pc_x, pc_y, nsh0_1221, nsh0_1223, \
                         nsg_722, nsg_873, nsg_875, nsh1_1221, nsh1_1223, \
                         osg_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = pa_x[k] * nsh0_1221[k]
                    + f_11 * nsg_873[k]
                    - f_8 * pc_x[k] * nsh1_1221[k];

        t_1222[k] = f_17 * nsg_722[k]
                    + f_3 * pc_y[k] * osg_872[k];

        t_1223[k] = pa_x[k] * nsh0_1223[k]
                    + f_11 * nsg_875[k]
                    - f_8 * pc_x[k] * nsh1_1223[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pa_x, pc_x, pc_y, pc_z, nsh0_1224, nsg_708, \
                         nsg_725, nsg_876, nsh1_1224, osg_873, \
                         osg_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = pa_x[k] * nsh0_1224[k]
                    + f_10 * nsg_876[k]
                    - f_8 * pc_x[k] * nsh1_1224[k];

        t_1225[k] = f_11 * nsg_708[k]
                    + f_3 * pc_z[k] * osg_873[k];

        t_1226[k] = f_17 * nsg_725[k]
                    + f_3 * pc_y[k] * osg_875[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, pa_x, pc_x, nsh0_1227, nsg_879, \
                         nsg_880, nsg_881, nsg_882, nsh1_1227, osg_880, osg_881, \
                         osg_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = pa_x[k] * nsh0_1227[k]
                    + f_10 * nsg_879[k]
                    - f_8 * pc_x[k] * nsh1_1227[k];

        t_1228[k] = f_9 * nsg_880[k]
                    + f_3 * pc_x[k] * osg_880[k];

        t_1229[k] = f_9 * nsg_881[k]
                    + f_3 * pc_x[k] * osg_881[k];

        t_1230[k] = f_9 * nsg_882[k]
                    + f_3 * pc_x[k] * osg_882[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pa_x, pc_x, pc_z, nsh0_1233, nsg_715, \
                         nsg_883, nsg_884, nsh1_1233, osg_880, osg_883, \
                         osg_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_9 * nsg_883[k]
                    + f_3 * pc_x[k] * osg_883[k];

        t_1232[k] = f_9 * nsg_884[k]
                    + f_3 * pc_x[k] * osg_884[k];

        t_1233[k] = pa_x[k] * nsh0_1233[k]
                    - f_8 * pc_x[k] * nsh1_1233[k];

        t_1234[k] = f_11 * nsg_715[k]
                    + f_3 * pc_z[k] * osg_880[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, pa_x, pc_x, pc_y, nsh0_1235, \
                         nsh0_1236, nsh0_1238, nsg_734, nsh1_1235, nsh1_1236, nsh1_1238, \
                         osg_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = pa_x[k] * nsh0_1235[k]
                    - f_8 * pc_x[k] * nsh1_1235[k];

        t_1236[k] = pa_x[k] * nsh0_1236[k]
                    - f_8 * pc_x[k] * nsh1_1236[k];

        t_1237[k] = f_17 * nsg_734[k]
                    + f_3 * pc_y[k] * osg_884[k];

        t_1238[k] = pa_x[k] * nsh0_1238[k]
                    - f_8 * pc_x[k] * nsh1_1238[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pa_x, pc_x, pc_y, pc_z, nsh0_1239, nsg_720, \
                         nsg_735, nsg_885, nsh1_1239, osg_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = pa_x[k] * nsh0_1239[k]
                    + f_20 * nsg_885[k]
                    - f_8 * pc_x[k] * nsh1_1239[k];

        t_1240[k] = f_19 * nsg_735[k]
                    + f_3 * pc_y[k] * osg_885[k];

        t_1241[k] = f_18 * nsg_720[k]
                    + f_3 * pc_z[k] * osg_885[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, pa_x, pc_x, pc_y, nsh0_1242, nsh0_1244, \
                         nsg_737, nsg_888, nsg_890, nsh1_1242, nsh1_1244, \
                         osg_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = pa_x[k] * nsh0_1242[k]
                    + f_11 * nsg_888[k]
                    - f_8 * pc_x[k] * nsh1_1242[k];

        t_1243[k] = f_19 * nsg_737[k]
                    + f_3 * pc_y[k] * osg_887[k];

        t_1244[k] = pa_x[k] * nsh0_1244[k]
                    + f_11 * nsg_890[k]
                    - f_8 * pc_x[k] * nsh1_1244[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, pa_x, pc_x, pc_y, pc_z, nsh0_1245, nsg_723, \
                         nsg_740, nsg_891, nsh1_1245, osg_888, \
                         osg_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = pa_x[k] * nsh0_1245[k]
                    + f_10 * nsg_891[k]
                    - f_8 * pc_x[k] * nsh1_1245[k];

        t_1246[k] = f_18 * nsg_723[k]
                    + f_3 * pc_z[k] * osg_888[k];

        t_1247[k] = f_19 * nsg_740[k]
                    + f_3 * pc_y[k] * osg_890[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, pa_x, pc_x, nsh0_1248, nsg_894, \
                         nsg_895, nsg_896, nsg_897, nsh1_1248, osg_895, osg_896, \
                         osg_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = pa_x[k] * nsh0_1248[k]
                    + f_10 * nsg_894[k]
                    - f_8 * pc_x[k] * nsh1_1248[k];

        t_1249[k] = f_9 * nsg_895[k]
                    + f_3 * pc_x[k] * osg_895[k];

        t_1250[k] = f_9 * nsg_896[k]
                    + f_3 * pc_x[k] * osg_896[k];

        t_1251[k] = f_9 * nsg_897[k]
                    + f_3 * pc_x[k] * osg_897[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pa_x, pc_x, pc_z, nsh0_1254, nsg_730, \
                         nsg_898, nsg_899, nsh1_1254, osg_895, osg_898, \
                         osg_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_9 * nsg_898[k]
                    + f_3 * pc_x[k] * osg_898[k];

        t_1253[k] = f_9 * nsg_899[k]
                    + f_3 * pc_x[k] * osg_899[k];

        t_1254[k] = pa_x[k] * nsh0_1254[k]
                    - f_8 * pc_x[k] * nsh1_1254[k];

        t_1255[k] = f_18 * nsg_730[k]
                    + f_3 * pc_z[k] * osg_895[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pa_x, pc_x, pc_y, nsh0_1256, \
                         nsh0_1257, nsh0_1259, nsg_749, nsh1_1256, nsh1_1257, nsh1_1259, \
                         osg_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = pa_x[k] * nsh0_1256[k]
                    - f_8 * pc_x[k] * nsh1_1256[k];

        t_1257[k] = pa_x[k] * nsh0_1257[k]
                    - f_8 * pc_x[k] * nsh1_1257[k];

        t_1258[k] = f_19 * nsg_749[k]
                    + f_3 * pc_y[k] * osg_899[k];

        t_1259[k] = pa_x[k] * nsh0_1259[k]
                    - f_8 * pc_x[k] * nsh1_1259[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, pa_x, pc_x, pc_y, pc_z, nsh0_1260, nsg_735, \
                         nsg_750, nsg_900, nsh1_1260, osg_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = pa_x[k] * nsh0_1260[k]
                    + f_20 * nsg_900[k]
                    - f_8 * pc_x[k] * nsh1_1260[k];

        t_1261[k] = f_20 * nsg_750[k]
                    + f_3 * pc_y[k] * osg_900[k];

        t_1262[k] = f_20 * nsg_735[k]
                    + f_3 * pc_z[k] * osg_900[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, pa_x, pc_x, pc_y, nsh0_1263, nsh0_1265, \
                         nsg_752, nsg_903, nsg_905, nsh1_1263, nsh1_1265, \
                         osg_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = pa_x[k] * nsh0_1263[k]
                    + f_11 * nsg_903[k]
                    - f_8 * pc_x[k] * nsh1_1263[k];

        t_1264[k] = f_20 * nsg_752[k]
                    + f_3 * pc_y[k] * osg_902[k];

        t_1265[k] = pa_x[k] * nsh0_1265[k]
                    + f_11 * nsg_905[k]
                    - f_8 * pc_x[k] * nsh1_1265[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, pa_x, pc_x, pc_y, pc_z, nsh0_1266, nsg_738, \
                         nsg_755, nsg_906, nsh1_1266, osg_903, \
                         osg_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = pa_x[k] * nsh0_1266[k]
                    + f_10 * nsg_906[k]
                    - f_8 * pc_x[k] * nsh1_1266[k];

        t_1267[k] = f_20 * nsg_738[k]
                    + f_3 * pc_z[k] * osg_903[k];

        t_1268[k] = f_20 * nsg_755[k]
                    + f_3 * pc_y[k] * osg_905[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, t_1272, pa_x, pc_x, nsh0_1269, nsg_909, \
                         nsg_910, nsg_911, nsg_912, nsh1_1269, osg_910, osg_911, \
                         osg_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = pa_x[k] * nsh0_1269[k]
                    + f_10 * nsg_909[k]
                    - f_8 * pc_x[k] * nsh1_1269[k];

        t_1270[k] = f_9 * nsg_910[k]
                    + f_3 * pc_x[k] * osg_910[k];

        t_1271[k] = f_9 * nsg_911[k]
                    + f_3 * pc_x[k] * osg_911[k];

        t_1272[k] = f_9 * nsg_912[k]
                    + f_3 * pc_x[k] * osg_912[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, pa_x, pc_x, pc_z, nsh0_1275, nsg_745, \
                         nsg_913, nsg_914, nsh1_1275, osg_910, osg_913, \
                         osg_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_9 * nsg_913[k]
                    + f_3 * pc_x[k] * osg_913[k];

        t_1274[k] = f_9 * nsg_914[k]
                    + f_3 * pc_x[k] * osg_914[k];

        t_1275[k] = pa_x[k] * nsh0_1275[k]
                    - f_8 * pc_x[k] * nsh1_1275[k];

        t_1276[k] = f_20 * nsg_745[k]
                    + f_3 * pc_z[k] * osg_910[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, t_1280, pa_x, pc_x, pc_y, nsh0_1277, \
                         nsh0_1278, nsh0_1280, nsg_764, nsh1_1277, nsh1_1278, nsh1_1280, \
                         osg_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = pa_x[k] * nsh0_1277[k]
                    - f_8 * pc_x[k] * nsh1_1277[k];

        t_1278[k] = pa_x[k] * nsh0_1278[k]
                    - f_8 * pc_x[k] * nsh1_1278[k];

        t_1279[k] = f_20 * nsg_764[k]
                    + f_3 * pc_y[k] * osg_914[k];

        t_1280[k] = pa_x[k] * nsh0_1280[k]
                    - f_8 * pc_x[k] * nsh1_1280[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pa_x, pc_x, pc_y, pc_z, nsh0_1281, nsg_750, \
                         nsg_765, nsg_915, nsh1_1281, osg_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = pa_x[k] * nsh0_1281[k]
                    + f_20 * nsg_915[k]
                    - f_8 * pc_x[k] * nsh1_1281[k];

        t_1282[k] = f_18 * nsg_765[k]
                    + f_3 * pc_y[k] * osg_915[k];

        t_1283[k] = f_19 * nsg_750[k]
                    + f_3 * pc_z[k] * osg_915[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, pa_x, pc_x, pc_y, nsh0_1284, nsh0_1286, \
                         nsg_767, nsg_918, nsg_920, nsh1_1284, nsh1_1286, \
                         osg_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = pa_x[k] * nsh0_1284[k]
                    + f_11 * nsg_918[k]
                    - f_8 * pc_x[k] * nsh1_1284[k];

        t_1285[k] = f_18 * nsg_767[k]
                    + f_3 * pc_y[k] * osg_917[k];

        t_1286[k] = pa_x[k] * nsh0_1286[k]
                    + f_11 * nsg_920[k]
                    - f_8 * pc_x[k] * nsh1_1286[k];
    }

#pragma omp simd aligned(t_1287, t_1288, t_1289, pa_x, pc_x, pc_y, pc_z, nsh0_1287, nsg_753, \
                         nsg_770, nsg_921, nsh1_1287, osg_918, \
                         osg_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1287[k] = pa_x[k] * nsh0_1287[k]
                    + f_10 * nsg_921[k]
                    - f_8 * pc_x[k] * nsh1_1287[k];

        t_1288[k] = f_19 * nsg_753[k]
                    + f_3 * pc_z[k] * osg_918[k];

        t_1289[k] = f_18 * nsg_770[k]
                    + f_3 * pc_y[k] * osg_920[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, t_1293, pa_x, pc_x, nsh0_1290, nsg_924, \
                         nsg_925, nsg_926, nsg_927, nsh1_1290, osg_925, osg_926, \
                         osg_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = pa_x[k] * nsh0_1290[k]
                    + f_10 * nsg_924[k]
                    - f_8 * pc_x[k] * nsh1_1290[k];

        t_1291[k] = f_9 * nsg_925[k]
                    + f_3 * pc_x[k] * osg_925[k];

        t_1292[k] = f_9 * nsg_926[k]
                    + f_3 * pc_x[k] * osg_926[k];

        t_1293[k] = f_9 * nsg_927[k]
                    + f_3 * pc_x[k] * osg_927[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsh0,
                                                           const size_t nsg, const size_t nsh1,
                                                           const size_t osf0, const size_t osf1,
                                                           const size_t osg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_1134 = buffer.data(nsh0 + 1134);
    const auto *nsh0_1139 = buffer.data(nsh0 + 1139);
    const auto *nsh0_1143 = buffer.data(nsh0 + 1143);
    const auto *nsh0_1155 = buffer.data(nsh0 + 1155);
    const auto *nsh0_1156 = buffer.data(nsh0 + 1156);
    const auto *nsh0_1158 = buffer.data(nsh0 + 1158);
    const auto *nsh0_1161 = buffer.data(nsh0 + 1161);
    const auto *nsh0_1296 = buffer.data(nsh0 + 1296);
    const auto *nsh0_1298 = buffer.data(nsh0 + 1298);
    const auto *nsh0_1299 = buffer.data(nsh0 + 1299);
    const auto *nsh0_1301 = buffer.data(nsh0 + 1301);
    const auto *nsh0_1302 = buffer.data(nsh0 + 1302);
    const auto *nsh0_1305 = buffer.data(nsh0 + 1305);
    const auto *nsh0_1307 = buffer.data(nsh0 + 1307);
    const auto *nsh0_1308 = buffer.data(nsh0 + 1308);
    const auto *nsh0_1311 = buffer.data(nsh0 + 1311);
    const auto *nsh0_1317 = buffer.data(nsh0 + 1317);
    const auto *nsh0_1319 = buffer.data(nsh0 + 1319);
    const auto *nsh0_1320 = buffer.data(nsh0 + 1320);
    const auto *nsh0_1322 = buffer.data(nsh0 + 1322);
    const auto *nsh0_1323 = buffer.data(nsh0 + 1323);
    const auto *nsh0_1326 = buffer.data(nsh0 + 1326);
    const auto *nsh0_1328 = buffer.data(nsh0 + 1328);
    const auto *nsh0_1329 = buffer.data(nsh0 + 1329);
    const auto *nsh0_1332 = buffer.data(nsh0 + 1332);
    const auto *nsh0_1338 = buffer.data(nsh0 + 1338);
    const auto *nsh0_1340 = buffer.data(nsh0 + 1340);
    const auto *nsh0_1341 = buffer.data(nsh0 + 1341);
    const auto *nsh0_1343 = buffer.data(nsh0 + 1343);
    const auto *nsh0_1347 = buffer.data(nsh0 + 1347);
    const auto *nsh0_1350 = buffer.data(nsh0 + 1350);
    const auto *nsh0_1359 = buffer.data(nsh0 + 1359);
    const auto *nsh0_1361 = buffer.data(nsh0 + 1361);
    const auto *nsh0_1362 = buffer.data(nsh0 + 1362);
    const auto *nsh0_1364 = buffer.data(nsh0 + 1364);
    const auto *nsh0_1365 = buffer.data(nsh0 + 1365);
    const auto *nsh0_1370 = buffer.data(nsh0 + 1370);
    const auto *nsh0_1374 = buffer.data(nsh0 + 1374);
    const auto *nsh0_1380 = buffer.data(nsh0 + 1380);
    const auto *nsh0_1381 = buffer.data(nsh0 + 1381);
    const auto *nsh0_1382 = buffer.data(nsh0 + 1382);
    const auto *nsh0_1383 = buffer.data(nsh0 + 1383);
    const auto *nsh0_1385 = buffer.data(nsh0 + 1385);

    const auto *nsg_760 = buffer.data(nsg + 760);
    const auto *nsg_765 = buffer.data(nsg + 765);
    const auto *nsg_768 = buffer.data(nsg + 768);
    const auto *nsg_775 = buffer.data(nsg + 775);
    const auto *nsg_779 = buffer.data(nsg + 779);
    const auto *nsg_780 = buffer.data(nsg + 780);
    const auto *nsg_782 = buffer.data(nsg + 782);
    const auto *nsg_783 = buffer.data(nsg + 783);
    const auto *nsg_785 = buffer.data(nsg + 785);
    const auto *nsg_790 = buffer.data(nsg + 790);
    const auto *nsg_794 = buffer.data(nsg + 794);
    const auto *nsg_795 = buffer.data(nsg + 795);
    const auto *nsg_797 = buffer.data(nsg + 797);
    const auto *nsg_798 = buffer.data(nsg + 798);
    const auto *nsg_800 = buffer.data(nsg + 800);
    const auto *nsg_805 = buffer.data(nsg + 805);
    const auto *nsg_809 = buffer.data(nsg + 809);
    const auto *nsg_810 = buffer.data(nsg + 810);
    const auto *nsg_812 = buffer.data(nsg + 812);
    const auto *nsg_815 = buffer.data(nsg + 815);
    const auto *nsg_824 = buffer.data(nsg + 824);
    const auto *nsg_835 = buffer.data(nsg + 835);
    const auto *nsg_839 = buffer.data(nsg + 839);
    const auto *nsg_928 = buffer.data(nsg + 928);
    const auto *nsg_929 = buffer.data(nsg + 929);
    const auto *nsg_930 = buffer.data(nsg + 930);
    const auto *nsg_933 = buffer.data(nsg + 933);
    const auto *nsg_935 = buffer.data(nsg + 935);
    const auto *nsg_936 = buffer.data(nsg + 936);
    const auto *nsg_939 = buffer.data(nsg + 939);
    const auto *nsg_940 = buffer.data(nsg + 940);
    const auto *nsg_941 = buffer.data(nsg + 941);
    const auto *nsg_942 = buffer.data(nsg + 942);
    const auto *nsg_943 = buffer.data(nsg + 943);
    const auto *nsg_944 = buffer.data(nsg + 944);
    const auto *nsg_945 = buffer.data(nsg + 945);
    const auto *nsg_948 = buffer.data(nsg + 948);
    const auto *nsg_950 = buffer.data(nsg + 950);
    const auto *nsg_951 = buffer.data(nsg + 951);
    const auto *nsg_954 = buffer.data(nsg + 954);
    const auto *nsg_955 = buffer.data(nsg + 955);
    const auto *nsg_956 = buffer.data(nsg + 956);
    const auto *nsg_957 = buffer.data(nsg + 957);
    const auto *nsg_958 = buffer.data(nsg + 958);
    const auto *nsg_959 = buffer.data(nsg + 959);
    const auto *nsg_963 = buffer.data(nsg + 963);
    const auto *nsg_966 = buffer.data(nsg + 966);
    const auto *nsg_970 = buffer.data(nsg + 970);
    const auto *nsg_971 = buffer.data(nsg + 971);
    const auto *nsg_972 = buffer.data(nsg + 972);
    const auto *nsg_973 = buffer.data(nsg + 973);
    const auto *nsg_974 = buffer.data(nsg + 974);
    const auto *nsg_975 = buffer.data(nsg + 975);
    const auto *nsg_980 = buffer.data(nsg + 980);
    const auto *nsg_984 = buffer.data(nsg + 984);
    const auto *nsg_985 = buffer.data(nsg + 985);
    const auto *nsg_986 = buffer.data(nsg + 986);
    const auto *nsg_987 = buffer.data(nsg + 987);
    const auto *nsg_989 = buffer.data(nsg + 989);

    const auto *nsh1_1134 = buffer.data(nsh1 + 1134);
    const auto *nsh1_1139 = buffer.data(nsh1 + 1139);
    const auto *nsh1_1143 = buffer.data(nsh1 + 1143);
    const auto *nsh1_1155 = buffer.data(nsh1 + 1155);
    const auto *nsh1_1156 = buffer.data(nsh1 + 1156);
    const auto *nsh1_1158 = buffer.data(nsh1 + 1158);
    const auto *nsh1_1161 = buffer.data(nsh1 + 1161);
    const auto *nsh1_1296 = buffer.data(nsh1 + 1296);
    const auto *nsh1_1298 = buffer.data(nsh1 + 1298);
    const auto *nsh1_1299 = buffer.data(nsh1 + 1299);
    const auto *nsh1_1301 = buffer.data(nsh1 + 1301);
    const auto *nsh1_1302 = buffer.data(nsh1 + 1302);
    const auto *nsh1_1305 = buffer.data(nsh1 + 1305);
    const auto *nsh1_1307 = buffer.data(nsh1 + 1307);
    const auto *nsh1_1308 = buffer.data(nsh1 + 1308);
    const auto *nsh1_1311 = buffer.data(nsh1 + 1311);
    const auto *nsh1_1317 = buffer.data(nsh1 + 1317);
    const auto *nsh1_1319 = buffer.data(nsh1 + 1319);
    const auto *nsh1_1320 = buffer.data(nsh1 + 1320);
    const auto *nsh1_1322 = buffer.data(nsh1 + 1322);
    const auto *nsh1_1323 = buffer.data(nsh1 + 1323);
    const auto *nsh1_1326 = buffer.data(nsh1 + 1326);
    const auto *nsh1_1328 = buffer.data(nsh1 + 1328);
    const auto *nsh1_1329 = buffer.data(nsh1 + 1329);
    const auto *nsh1_1332 = buffer.data(nsh1 + 1332);
    const auto *nsh1_1338 = buffer.data(nsh1 + 1338);
    const auto *nsh1_1340 = buffer.data(nsh1 + 1340);
    const auto *nsh1_1341 = buffer.data(nsh1 + 1341);
    const auto *nsh1_1343 = buffer.data(nsh1 + 1343);
    const auto *nsh1_1347 = buffer.data(nsh1 + 1347);
    const auto *nsh1_1350 = buffer.data(nsh1 + 1350);
    const auto *nsh1_1359 = buffer.data(nsh1 + 1359);
    const auto *nsh1_1361 = buffer.data(nsh1 + 1361);
    const auto *nsh1_1362 = buffer.data(nsh1 + 1362);
    const auto *nsh1_1364 = buffer.data(nsh1 + 1364);
    const auto *nsh1_1365 = buffer.data(nsh1 + 1365);
    const auto *nsh1_1370 = buffer.data(nsh1 + 1370);
    const auto *nsh1_1374 = buffer.data(nsh1 + 1374);
    const auto *nsh1_1380 = buffer.data(nsh1 + 1380);
    const auto *nsh1_1381 = buffer.data(nsh1 + 1381);
    const auto *nsh1_1382 = buffer.data(nsh1 + 1382);
    const auto *nsh1_1383 = buffer.data(nsh1 + 1383);
    const auto *nsh1_1385 = buffer.data(nsh1 + 1385);

    const auto *osf0_650 = buffer.data(osf0 + 650);
    const auto *osf0_651 = buffer.data(osf0 + 651);
    const auto *osf0_652 = buffer.data(osf0 + 652);
    const auto *osf0_660 = buffer.data(osf0 + 660);
    const auto *osf0_661 = buffer.data(osf0 + 661);
    const auto *osf0_663 = buffer.data(osf0 + 663);
    const auto *osf0_665 = buffer.data(osf0 + 665);
    const auto *osf0_666 = buffer.data(osf0 + 666);
    const auto *osf0_667 = buffer.data(osf0 + 667);
    const auto *osf0_668 = buffer.data(osf0 + 668);
    const auto *osf0_669 = buffer.data(osf0 + 669);
    const auto *osf0_672 = buffer.data(osf0 + 672);
    const auto *osf0_674 = buffer.data(osf0 + 674);
    const auto *osf0_675 = buffer.data(osf0 + 675);
    const auto *osf0_677 = buffer.data(osf0 + 677);
    const auto *osf0_678 = buffer.data(osf0 + 678);
    const auto *osf0_679 = buffer.data(osf0 + 679);

    const auto *osf1_650 = buffer.data(osf1 + 650);
    const auto *osf1_651 = buffer.data(osf1 + 651);
    const auto *osf1_652 = buffer.data(osf1 + 652);
    const auto *osf1_660 = buffer.data(osf1 + 660);
    const auto *osf1_661 = buffer.data(osf1 + 661);
    const auto *osf1_663 = buffer.data(osf1 + 663);
    const auto *osf1_665 = buffer.data(osf1 + 665);
    const auto *osf1_666 = buffer.data(osf1 + 666);
    const auto *osf1_667 = buffer.data(osf1 + 667);
    const auto *osf1_668 = buffer.data(osf1 + 668);
    const auto *osf1_669 = buffer.data(osf1 + 669);
    const auto *osf1_672 = buffer.data(osf1 + 672);
    const auto *osf1_674 = buffer.data(osf1 + 674);
    const auto *osf1_675 = buffer.data(osf1 + 675);
    const auto *osf1_677 = buffer.data(osf1 + 677);
    const auto *osf1_678 = buffer.data(osf1 + 678);
    const auto *osf1_679 = buffer.data(osf1 + 679);

    const auto *osg_925 = buffer.data(osg + 925);
    const auto *osg_928 = buffer.data(osg + 928);
    const auto *osg_929 = buffer.data(osg + 929);
    const auto *osg_930 = buffer.data(osg + 930);
    const auto *osg_932 = buffer.data(osg + 932);
    const auto *osg_933 = buffer.data(osg + 933);
    const auto *osg_935 = buffer.data(osg + 935);
    const auto *osg_940 = buffer.data(osg + 940);
    const auto *osg_941 = buffer.data(osg + 941);
    const auto *osg_942 = buffer.data(osg + 942);
    const auto *osg_943 = buffer.data(osg + 943);
    const auto *osg_944 = buffer.data(osg + 944);
    const auto *osg_945 = buffer.data(osg + 945);
    const auto *osg_947 = buffer.data(osg + 947);
    const auto *osg_948 = buffer.data(osg + 948);
    const auto *osg_950 = buffer.data(osg + 950);
    const auto *osg_955 = buffer.data(osg + 955);
    const auto *osg_956 = buffer.data(osg + 956);
    const auto *osg_957 = buffer.data(osg + 957);
    const auto *osg_958 = buffer.data(osg + 958);
    const auto *osg_959 = buffer.data(osg + 959);
    const auto *osg_960 = buffer.data(osg + 960);
    const auto *osg_962 = buffer.data(osg + 962);
    const auto *osg_963 = buffer.data(osg + 963);
    const auto *osg_965 = buffer.data(osg + 965);
    const auto *osg_970 = buffer.data(osg + 970);
    const auto *osg_971 = buffer.data(osg + 971);
    const auto *osg_972 = buffer.data(osg + 972);
    const auto *osg_973 = buffer.data(osg + 973);
    const auto *osg_974 = buffer.data(osg + 974);
    const auto *osg_975 = buffer.data(osg + 975);
    const auto *osg_976 = buffer.data(osg + 976);
    const auto *osg_977 = buffer.data(osg + 977);
    const auto *osg_978 = buffer.data(osg + 978);
    const auto *osg_979 = buffer.data(osg + 979);
    const auto *osg_980 = buffer.data(osg + 980);
    const auto *osg_984 = buffer.data(osg + 984);
    const auto *osg_985 = buffer.data(osg + 985);
    const auto *osg_986 = buffer.data(osg + 986);
    const auto *osg_987 = buffer.data(osg + 987);
    const auto *osg_989 = buffer.data(osg + 989);
    const auto *osg_990 = buffer.data(osg + 990);
    const auto *osg_991 = buffer.data(osg + 991);
    const auto *osg_993 = buffer.data(osg + 993);
    const auto *osg_995 = buffer.data(osg + 995);
    const auto *osg_996 = buffer.data(osg + 996);
    const auto *osg_998 = buffer.data(osg + 998);
    const auto *osg_999 = buffer.data(osg + 999);
    const auto *osg_1000 = buffer.data(osg + 1000);
    const auto *osg_1001 = buffer.data(osg + 1001);
    const auto *osg_1002 = buffer.data(osg + 1002);
    const auto *osg_1003 = buffer.data(osg + 1003);
    const auto *osg_1004 = buffer.data(osg + 1004);
    const auto *osg_1007 = buffer.data(osg + 1007);
    const auto *osg_1009 = buffer.data(osg + 1009);
    const auto *osg_1010 = buffer.data(osg + 1010);
    const auto *osg_1012 = buffer.data(osg + 1012);
    const auto *osg_1013 = buffer.data(osg + 1013);
    const auto *osg_1014 = buffer.data(osg + 1014);
    const auto *osg_1015 = buffer.data(osg + 1015);
    const auto *osg_1016 = buffer.data(osg + 1016);
    const auto *osg_1017 = buffer.data(osg + 1017);

#pragma omp simd aligned(t_1294, t_1295, t_1296, t_1297, pa_x, pc_x, pc_z, nsh0_1296, nsg_760, \
                         nsg_928, nsg_929, nsh1_1296, osg_925, osg_928, \
                         osg_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1294[k] = f_9 * nsg_928[k]
                    + f_3 * pc_x[k] * osg_928[k];

        t_1295[k] = f_9 * nsg_929[k]
                    + f_3 * pc_x[k] * osg_929[k];

        t_1296[k] = pa_x[k] * nsh0_1296[k]
                    - f_8 * pc_x[k] * nsh1_1296[k];

        t_1297[k] = f_19 * nsg_760[k]
                    + f_3 * pc_z[k] * osg_925[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, t_1301, pa_x, pc_x, pc_y, nsh0_1298, \
                         nsh0_1299, nsh0_1301, nsg_779, nsh1_1298, nsh1_1299, nsh1_1301, \
                         osg_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pa_x[k] * nsh0_1298[k]
                    - f_8 * pc_x[k] * nsh1_1298[k];

        t_1299[k] = pa_x[k] * nsh0_1299[k]
                    - f_8 * pc_x[k] * nsh1_1299[k];

        t_1300[k] = f_18 * nsg_779[k]
                    + f_3 * pc_y[k] * osg_929[k];

        t_1301[k] = pa_x[k] * nsh0_1301[k]
                    - f_8 * pc_x[k] * nsh1_1301[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, pa_x, pc_x, pc_y, pc_z, nsh0_1302, nsg_765, \
                         nsg_780, nsg_930, nsh1_1302, osg_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = pa_x[k] * nsh0_1302[k]
                    + f_20 * nsg_930[k]
                    - f_8 * pc_x[k] * nsh1_1302[k];

        t_1303[k] = f_11 * nsg_780[k]
                    + f_3 * pc_y[k] * osg_930[k];

        t_1304[k] = f_17 * nsg_765[k]
                    + f_3 * pc_z[k] * osg_930[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, pa_x, pc_x, pc_y, nsh0_1305, nsh0_1307, \
                         nsg_782, nsg_933, nsg_935, nsh1_1305, nsh1_1307, \
                         osg_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = pa_x[k] * nsh0_1305[k]
                    + f_11 * nsg_933[k]
                    - f_8 * pc_x[k] * nsh1_1305[k];

        t_1306[k] = f_11 * nsg_782[k]
                    + f_3 * pc_y[k] * osg_932[k];

        t_1307[k] = pa_x[k] * nsh0_1307[k]
                    + f_11 * nsg_935[k]
                    - f_8 * pc_x[k] * nsh1_1307[k];
    }

#pragma omp simd aligned(t_1308, t_1309, t_1310, pa_x, pc_x, pc_y, pc_z, nsh0_1308, nsg_768, \
                         nsg_785, nsg_936, nsh1_1308, osg_933, \
                         osg_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1308[k] = pa_x[k] * nsh0_1308[k]
                    + f_10 * nsg_936[k]
                    - f_8 * pc_x[k] * nsh1_1308[k];

        t_1309[k] = f_17 * nsg_768[k]
                    + f_3 * pc_z[k] * osg_933[k];

        t_1310[k] = f_11 * nsg_785[k]
                    + f_3 * pc_y[k] * osg_935[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, t_1314, pa_x, pc_x, nsh0_1311, nsg_939, \
                         nsg_940, nsg_941, nsg_942, nsh1_1311, osg_940, osg_941, \
                         osg_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = pa_x[k] * nsh0_1311[k]
                    + f_10 * nsg_939[k]
                    - f_8 * pc_x[k] * nsh1_1311[k];

        t_1312[k] = f_9 * nsg_940[k]
                    + f_3 * pc_x[k] * osg_940[k];

        t_1313[k] = f_9 * nsg_941[k]
                    + f_3 * pc_x[k] * osg_941[k];

        t_1314[k] = f_9 * nsg_942[k]
                    + f_3 * pc_x[k] * osg_942[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, pa_x, pc_x, pc_z, nsh0_1317, nsg_775, \
                         nsg_943, nsg_944, nsh1_1317, osg_940, osg_943, \
                         osg_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = f_9 * nsg_943[k]
                    + f_3 * pc_x[k] * osg_943[k];

        t_1316[k] = f_9 * nsg_944[k]
                    + f_3 * pc_x[k] * osg_944[k];

        t_1317[k] = pa_x[k] * nsh0_1317[k]
                    - f_8 * pc_x[k] * nsh1_1317[k];

        t_1318[k] = f_17 * nsg_775[k]
                    + f_3 * pc_z[k] * osg_940[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, t_1322, pa_x, pc_x, pc_y, nsh0_1319, \
                         nsh0_1320, nsh0_1322, nsg_794, nsh1_1319, nsh1_1320, nsh1_1322, \
                         osg_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = pa_x[k] * nsh0_1319[k]
                    - f_8 * pc_x[k] * nsh1_1319[k];

        t_1320[k] = pa_x[k] * nsh0_1320[k]
                    - f_8 * pc_x[k] * nsh1_1320[k];

        t_1321[k] = f_11 * nsg_794[k]
                    + f_3 * pc_y[k] * osg_944[k];

        t_1322[k] = pa_x[k] * nsh0_1322[k]
                    - f_8 * pc_x[k] * nsh1_1322[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, pa_x, pc_x, pc_y, pc_z, nsh0_1323, nsg_780, \
                         nsg_795, nsg_945, nsh1_1323, osg_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = pa_x[k] * nsh0_1323[k]
                    + f_20 * nsg_945[k]
                    - f_8 * pc_x[k] * nsh1_1323[k];

        t_1324[k] = f_10 * nsg_795[k]
                    + f_3 * pc_y[k] * osg_945[k];

        t_1325[k] = f_16 * nsg_780[k]
                    + f_3 * pc_z[k] * osg_945[k];
    }

#pragma omp simd aligned(t_1326, t_1327, t_1328, pa_x, pc_x, pc_y, nsh0_1326, nsh0_1328, \
                         nsg_797, nsg_948, nsg_950, nsh1_1326, nsh1_1328, \
                         osg_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1326[k] = pa_x[k] * nsh0_1326[k]
                    + f_11 * nsg_948[k]
                    - f_8 * pc_x[k] * nsh1_1326[k];

        t_1327[k] = f_10 * nsg_797[k]
                    + f_3 * pc_y[k] * osg_947[k];

        t_1328[k] = pa_x[k] * nsh0_1328[k]
                    + f_11 * nsg_950[k]
                    - f_8 * pc_x[k] * nsh1_1328[k];
    }

#pragma omp simd aligned(t_1329, t_1330, t_1331, pa_x, pc_x, pc_y, pc_z, nsh0_1329, nsg_783, \
                         nsg_800, nsg_951, nsh1_1329, osg_948, \
                         osg_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1329[k] = pa_x[k] * nsh0_1329[k]
                    + f_10 * nsg_951[k]
                    - f_8 * pc_x[k] * nsh1_1329[k];

        t_1330[k] = f_16 * nsg_783[k]
                    + f_3 * pc_z[k] * osg_948[k];

        t_1331[k] = f_10 * nsg_800[k]
                    + f_3 * pc_y[k] * osg_950[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, pa_x, pc_x, nsh0_1332, nsg_954, \
                         nsg_955, nsg_956, nsg_957, nsh1_1332, osg_955, osg_956, \
                         osg_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = pa_x[k] * nsh0_1332[k]
                    + f_10 * nsg_954[k]
                    - f_8 * pc_x[k] * nsh1_1332[k];

        t_1333[k] = f_9 * nsg_955[k]
                    + f_3 * pc_x[k] * osg_955[k];

        t_1334[k] = f_9 * nsg_956[k]
                    + f_3 * pc_x[k] * osg_956[k];

        t_1335[k] = f_9 * nsg_957[k]
                    + f_3 * pc_x[k] * osg_957[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, t_1339, pa_x, pc_x, pc_z, nsh0_1338, nsg_790, \
                         nsg_958, nsg_959, nsh1_1338, osg_955, osg_958, \
                         osg_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_9 * nsg_958[k]
                    + f_3 * pc_x[k] * osg_958[k];

        t_1337[k] = f_9 * nsg_959[k]
                    + f_3 * pc_x[k] * osg_959[k];

        t_1338[k] = pa_x[k] * nsh0_1338[k]
                    - f_8 * pc_x[k] * nsh1_1338[k];

        t_1339[k] = f_16 * nsg_790[k]
                    + f_3 * pc_z[k] * osg_955[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, t_1343, pa_x, pc_x, pc_y, nsh0_1340, \
                         nsh0_1341, nsh0_1343, nsg_809, nsh1_1340, nsh1_1341, nsh1_1343, \
                         osg_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = pa_x[k] * nsh0_1340[k]
                    - f_8 * pc_x[k] * nsh1_1340[k];

        t_1341[k] = pa_x[k] * nsh0_1341[k]
                    - f_8 * pc_x[k] * nsh1_1341[k];

        t_1342[k] = f_10 * nsg_809[k]
                    + f_3 * pc_y[k] * osg_959[k];

        t_1343[k] = pa_x[k] * nsh0_1343[k]
                    - f_8 * pc_x[k] * nsh1_1343[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, pa_y, pc_y, pc_z, nsh0_1134, nsg_795, \
                         nsg_810, nsh1_1134, osg_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = pa_y[k] * nsh0_1134[k]
                    - f_8 * pc_y[k] * nsh1_1134[k];

        t_1345[k] = f_9 * nsg_810[k]
                    + f_3 * pc_y[k] * osg_960[k];

        t_1346[k] = f_15 * nsg_795[k]
                    + f_3 * pc_z[k] * osg_960[k];
    }

#pragma omp simd aligned(t_1347, t_1348, t_1349, pa_x, pa_y, pc_x, pc_y, nsh0_1139, nsh0_1347, \
                         nsg_812, nsg_963, nsh1_1139, nsh1_1347, \
                         osg_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1347[k] = pa_x[k] * nsh0_1347[k]
                    + f_11 * nsg_963[k]
                    - f_8 * pc_x[k] * nsh1_1347[k];

        t_1348[k] = f_9 * nsg_812[k]
                    + f_3 * pc_y[k] * osg_962[k];

        t_1349[k] = pa_y[k] * nsh0_1139[k]
                    - f_8 * pc_y[k] * nsh1_1139[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, pa_x, pc_x, pc_y, pc_z, nsh0_1350, nsg_798, \
                         nsg_815, nsg_966, nsh1_1350, osg_963, \
                         osg_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = pa_x[k] * nsh0_1350[k]
                    + f_10 * nsg_966[k]
                    - f_8 * pc_x[k] * nsh1_1350[k];

        t_1351[k] = f_15 * nsg_798[k]
                    + f_3 * pc_z[k] * osg_963[k];

        t_1352[k] = f_9 * nsg_815[k]
                    + f_3 * pc_y[k] * osg_965[k];
    }

#pragma omp simd aligned(t_1353, t_1354, t_1355, t_1356, pa_y, pc_x, pc_y, nsh0_1143, nsg_970, \
                         nsg_971, nsg_972, nsh1_1143, osg_970, osg_971, \
                         osg_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1353[k] = pa_y[k] * nsh0_1143[k]
                    - f_8 * pc_y[k] * nsh1_1143[k];

        t_1354[k] = f_9 * nsg_970[k]
                    + f_3 * pc_x[k] * osg_970[k];

        t_1355[k] = f_9 * nsg_971[k]
                    + f_3 * pc_x[k] * osg_971[k];

        t_1356[k] = f_9 * nsg_972[k]
                    + f_3 * pc_x[k] * osg_972[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, t_1360, pa_x, pc_x, pc_z, nsh0_1359, nsg_805, \
                         nsg_973, nsg_974, nsh1_1359, osg_970, osg_973, \
                         osg_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_9 * nsg_973[k]
                    + f_3 * pc_x[k] * osg_973[k];

        t_1358[k] = f_9 * nsg_974[k]
                    + f_3 * pc_x[k] * osg_974[k];

        t_1359[k] = pa_x[k] * nsh0_1359[k]
                    - f_8 * pc_x[k] * nsh1_1359[k];

        t_1360[k] = f_15 * nsg_805[k]
                    + f_3 * pc_z[k] * osg_970[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, t_1364, pa_x, pc_x, pc_y, nsh0_1361, \
                         nsh0_1362, nsh0_1364, nsg_824, nsh1_1361, nsh1_1362, nsh1_1364, \
                         osg_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = pa_x[k] * nsh0_1361[k]
                    - f_8 * pc_x[k] * nsh1_1361[k];

        t_1362[k] = pa_x[k] * nsh0_1362[k]
                    - f_8 * pc_x[k] * nsh1_1362[k];

        t_1363[k] = f_9 * nsg_824[k]
                    + f_3 * pc_y[k] * osg_974[k];

        t_1364[k] = pa_x[k] * nsh0_1364[k]
                    - f_8 * pc_x[k] * nsh1_1364[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, t_1368, pa_x, pc_x, pc_y, pc_z, nsh0_1365, \
                         nsg_810, nsg_975, nsh1_1365, osf0_650, osf1_650, osg_975, \
                         osg_976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = pa_x[k] * nsh0_1365[k]
                    + f_20 * nsg_975[k]
                    - f_8 * pc_x[k] * nsh1_1365[k];

        t_1366[k] = f_3 * pc_y[k] * osg_975[k];

        t_1367[k] = f_12 * nsg_810[k]
                    + f_3 * pc_z[k] * osg_975[k];

        t_1368[k] = f_4 * osf0_650[k]
                    - f_5 * osf1_650[k]
                    + f_3 * pc_y[k] * osg_976[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, pa_x, pc_x, pc_y, nsh0_1370, nsg_980, \
                         nsh1_1370, osf0_651, osf1_651, osg_977, \
                         osg_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_3 * pc_y[k] * osg_977[k];

        t_1370[k] = pa_x[k] * nsh0_1370[k]
                    + f_11 * nsg_980[k]
                    - f_8 * pc_x[k] * nsh1_1370[k];

        t_1371[k] = f_6 * osf0_651[k]
                    - f_7 * osf1_651[k]
                    + f_3 * pc_y[k] * osg_978[k];
    }

#pragma omp simd aligned(t_1372, t_1373, t_1374, t_1375, pa_x, pc_x, pc_y, nsh0_1374, nsg_984, \
                         nsg_985, nsh1_1374, osf0_652, osf1_652, osg_979, osg_980, \
                         osg_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1372[k] = f_4 * osf0_652[k]
                    - f_5 * osf1_652[k]
                    + f_3 * pc_y[k] * osg_979[k];

        t_1373[k] = f_3 * pc_y[k] * osg_980[k];

        t_1374[k] = pa_x[k] * nsh0_1374[k]
                    + f_10 * nsg_984[k]
                    - f_8 * pc_x[k] * nsh1_1374[k];

        t_1375[k] = f_9 * nsg_985[k]
                    + f_3 * pc_x[k] * osg_985[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, t_1379, pc_x, pc_y, nsg_986, nsg_987, \
                         nsg_989, osg_984, osg_986, osg_987, osg_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_9 * nsg_986[k]
                    + f_3 * pc_x[k] * osg_986[k];

        t_1377[k] = f_9 * nsg_987[k]
                    + f_3 * pc_x[k] * osg_987[k];

        t_1378[k] = f_3 * pc_y[k] * osg_984[k];

        t_1379[k] = f_9 * nsg_989[k]
                    + f_3 * pc_x[k] * osg_989[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, pa_x, pc_x, nsh0_1380, nsh0_1381, \
                         nsh0_1382, nsh0_1383, nsh1_1380, nsh1_1381, nsh1_1382, \
                         nsh1_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = pa_x[k] * nsh0_1380[k]
                    - f_8 * pc_x[k] * nsh1_1380[k];

        t_1381[k] = pa_x[k] * nsh0_1381[k]
                    - f_8 * pc_x[k] * nsh1_1381[k];

        t_1382[k] = pa_x[k] * nsh0_1382[k]
                    - f_8 * pc_x[k] * nsh1_1382[k];

        t_1383[k] = pa_x[k] * nsh0_1383[k]
                    - f_8 * pc_x[k] * nsh1_1383[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, t_1387, pa_x, pc_x, pc_y, nsh0_1385, \
                         nsh1_1385, osf0_660, osf0_661, osf1_660, osf1_661, osg_989, osg_990, \
                         osg_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_3 * pc_y[k] * osg_989[k];

        t_1385[k] = pa_x[k] * nsh0_1385[k]
                    - f_8 * pc_x[k] * nsh1_1385[k];

        t_1386[k] = f_1 * osf0_660[k]
                    - f_2 * osf1_660[k]
                    + f_3 * pc_x[k] * osg_990[k];

        t_1387[k] = f_13 * osf0_661[k]
                    - f_14 * osf1_661[k]
                    + f_3 * pc_x[k] * osg_991[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, t_1391, pc_x, pc_z, osf0_663, osf0_665, \
                         osf1_663, osf1_665, osg_990, osg_991, osg_993, \
                         osg_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_3 * pc_z[k] * osg_990[k];

        t_1389[k] = f_6 * osf0_663[k]
                    - f_7 * osf1_663[k]
                    + f_3 * pc_x[k] * osg_993[k];

        t_1390[k] = f_3 * pc_z[k] * osg_991[k];

        t_1391[k] = f_6 * osf0_665[k]
                    - f_7 * osf1_665[k]
                    + f_3 * pc_x[k] * osg_995[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, t_1395, pc_x, pc_z, osf0_666, osf0_668, \
                         osf0_669, osf1_666, osf1_668, osf1_669, osg_993, osg_996, osg_998, \
                         osg_999 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_4 * osf0_666[k]
                    - f_5 * osf1_666[k]
                    + f_3 * pc_x[k] * osg_996[k];

        t_1393[k] = f_3 * pc_z[k] * osg_993[k];

        t_1394[k] = f_4 * osf0_668[k]
                    - f_5 * osf1_668[k]
                    + f_3 * pc_x[k] * osg_998[k];

        t_1395[k] = f_4 * osf0_669[k]
                    - f_5 * osf1_669[k]
                    + f_3 * pc_x[k] * osg_999[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, t_1399, t_1400, t_1401, pc_x, pc_y, nsg_835, \
                         osf0_666, osf1_666, osg_1000, osg_1001, osg_1002, osg_1003, \
                         osg_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_3 * pc_x[k] * osg_1000[k];

        t_1397[k] = f_3 * pc_x[k] * osg_1001[k];

        t_1398[k] = f_3 * pc_x[k] * osg_1002[k];

        t_1399[k] = f_3 * pc_x[k] * osg_1003[k];

        t_1400[k] = f_3 * pc_x[k] * osg_1004[k];

        t_1401[k] = f_0 * nsg_835[k]
                    + f_1 * osf0_666[k]
                    - f_2 * osf1_666[k]
                    + f_3 * pc_y[k] * osg_1000[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, t_1405, pc_y, pc_z, nsg_839, osf0_666, \
                         osf0_667, osf1_666, osf1_667, osg_1000, osg_1001, osg_1002, \
                         osg_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_3 * pc_z[k] * osg_1000[k];

        t_1403[k] = f_4 * osf0_666[k]
                    - f_5 * osf1_666[k]
                    + f_3 * pc_z[k] * osg_1001[k];

        t_1404[k] = f_6 * osf0_667[k]
                    - f_7 * osf1_667[k]
                    + f_3 * pc_z[k] * osg_1002[k];

        t_1405[k] = f_0 * nsg_839[k]
                    + f_3 * pc_y[k] * osg_1004[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pa_z, pc_z, nsh0_1155, nsh0_1156, nsh1_1155, \
                         nsh1_1156, osf0_669, osf1_669, osg_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_1 * osf0_669[k]
                    - f_2 * osf1_669[k]
                    + f_3 * pc_z[k] * osg_1004[k];

        t_1407[k] = pa_z[k] * nsh0_1155[k]
                    - f_8 * pc_z[k] * nsh1_1155[k];

        t_1408[k] = pa_z[k] * nsh0_1156[k]
                    - f_8 * pc_z[k] * nsh1_1156[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pa_z, pc_x, pc_z, nsh0_1158, nsh1_1158, \
                         osf0_672, osf0_674, osf1_672, osf1_674, osg_1007, \
                         osg_1009 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_13 * osf0_672[k]
                    - f_14 * osf1_672[k]
                    + f_3 * pc_x[k] * osg_1007[k];

        t_1410[k] = pa_z[k] * nsh0_1158[k]
                    - f_8 * pc_z[k] * nsh1_1158[k];

        t_1411[k] = f_6 * osf0_674[k]
                    - f_7 * osf1_674[k]
                    + f_3 * pc_x[k] * osg_1009[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pa_z, pc_x, pc_z, nsh0_1161, nsh1_1161, \
                         osf0_675, osf0_677, osf1_675, osf1_677, osg_1010, \
                         osg_1012 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_6 * osf0_675[k]
                    - f_7 * osf1_675[k]
                    + f_3 * pc_x[k] * osg_1010[k];

        t_1413[k] = pa_z[k] * nsh0_1161[k]
                    - f_8 * pc_z[k] * nsh1_1161[k];

        t_1414[k] = f_4 * osf0_677[k]
                    - f_5 * osf1_677[k]
                    + f_3 * pc_x[k] * osg_1012[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, t_1418, t_1419, pc_x, osf0_678, osf0_679, \
                         osf1_678, osf1_679, osg_1013, osg_1014, osg_1015, osg_1016, \
                         osg_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_4 * osf0_678[k]
                    - f_5 * osf1_678[k]
                    + f_3 * pc_x[k] * osg_1013[k];

        t_1416[k] = f_4 * osf0_679[k]
                    - f_5 * osf1_679[k]
                    + f_3 * pc_x[k] * osg_1014[k];

        t_1417[k] = f_3 * pc_x[k] * osg_1015[k];

        t_1418[k] = f_3 * pc_x[k] * osg_1016[k];

        t_1419[k] = f_3 * pc_x[k] * osg_1017[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsh0,
                                                           const size_t nsg, const size_t nsh1,
                                                           const size_t osf0, const size_t osf1,
                                                           const size_t osg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_1170 = buffer.data(nsh0 + 1170);
    const auto *nsh0_1172 = buffer.data(nsh0 + 1172);
    const auto *nsh0_1173 = buffer.data(nsh0 + 1173);

    const auto *nsg_835 = buffer.data(nsg + 835);
    const auto *nsg_836 = buffer.data(nsg + 836);
    const auto *nsg_837 = buffer.data(nsg + 837);
    const auto *nsg_839 = buffer.data(nsg + 839);
    const auto *nsg_850 = buffer.data(nsg + 850);
    const auto *nsg_854 = buffer.data(nsg + 854);
    const auto *nsg_865 = buffer.data(nsg + 865);
    const auto *nsg_867 = buffer.data(nsg + 867);
    const auto *nsg_868 = buffer.data(nsg + 868);
    const auto *nsg_869 = buffer.data(nsg + 869);
    const auto *nsg_880 = buffer.data(nsg + 880);
    const auto *nsg_882 = buffer.data(nsg + 882);
    const auto *nsg_883 = buffer.data(nsg + 883);
    const auto *nsg_884 = buffer.data(nsg + 884);
    const auto *nsg_895 = buffer.data(nsg + 895);
    const auto *nsg_897 = buffer.data(nsg + 897);
    const auto *nsg_898 = buffer.data(nsg + 898);
    const auto *nsg_899 = buffer.data(nsg + 899);
    const auto *nsg_910 = buffer.data(nsg + 910);
    const auto *nsg_912 = buffer.data(nsg + 912);
    const auto *nsg_913 = buffer.data(nsg + 913);
    const auto *nsg_914 = buffer.data(nsg + 914);
    const auto *nsg_925 = buffer.data(nsg + 925);
    const auto *nsg_927 = buffer.data(nsg + 927);
    const auto *nsg_928 = buffer.data(nsg + 928);
    const auto *nsg_929 = buffer.data(nsg + 929);

    const auto *nsh1_1170 = buffer.data(nsh1 + 1170);
    const auto *nsh1_1172 = buffer.data(nsh1 + 1172);
    const auto *nsh1_1173 = buffer.data(nsh1 + 1173);

    const auto *osf0_679 = buffer.data(osf0 + 679);
    const auto *osf0_680 = buffer.data(osf0 + 680);
    const auto *osf0_681 = buffer.data(osf0 + 681);
    const auto *osf0_682 = buffer.data(osf0 + 682);
    const auto *osf0_683 = buffer.data(osf0 + 683);
    const auto *osf0_684 = buffer.data(osf0 + 684);
    const auto *osf0_685 = buffer.data(osf0 + 685);
    const auto *osf0_686 = buffer.data(osf0 + 686);
    const auto *osf0_687 = buffer.data(osf0 + 687);
    const auto *osf0_688 = buffer.data(osf0 + 688);
    const auto *osf0_689 = buffer.data(osf0 + 689);
    const auto *osf0_690 = buffer.data(osf0 + 690);
    const auto *osf0_691 = buffer.data(osf0 + 691);
    const auto *osf0_692 = buffer.data(osf0 + 692);
    const auto *osf0_693 = buffer.data(osf0 + 693);
    const auto *osf0_694 = buffer.data(osf0 + 694);
    const auto *osf0_695 = buffer.data(osf0 + 695);
    const auto *osf0_696 = buffer.data(osf0 + 696);
    const auto *osf0_697 = buffer.data(osf0 + 697);
    const auto *osf0_698 = buffer.data(osf0 + 698);
    const auto *osf0_699 = buffer.data(osf0 + 699);
    const auto *osf0_700 = buffer.data(osf0 + 700);
    const auto *osf0_701 = buffer.data(osf0 + 701);
    const auto *osf0_702 = buffer.data(osf0 + 702);
    const auto *osf0_703 = buffer.data(osf0 + 703);
    const auto *osf0_704 = buffer.data(osf0 + 704);
    const auto *osf0_705 = buffer.data(osf0 + 705);
    const auto *osf0_706 = buffer.data(osf0 + 706);
    const auto *osf0_707 = buffer.data(osf0 + 707);
    const auto *osf0_708 = buffer.data(osf0 + 708);
    const auto *osf0_709 = buffer.data(osf0 + 709);
    const auto *osf0_710 = buffer.data(osf0 + 710);
    const auto *osf0_711 = buffer.data(osf0 + 711);
    const auto *osf0_712 = buffer.data(osf0 + 712);
    const auto *osf0_713 = buffer.data(osf0 + 713);
    const auto *osf0_714 = buffer.data(osf0 + 714);
    const auto *osf0_715 = buffer.data(osf0 + 715);
    const auto *osf0_716 = buffer.data(osf0 + 716);
    const auto *osf0_717 = buffer.data(osf0 + 717);
    const auto *osf0_718 = buffer.data(osf0 + 718);
    const auto *osf0_719 = buffer.data(osf0 + 719);
    const auto *osf0_720 = buffer.data(osf0 + 720);
    const auto *osf0_721 = buffer.data(osf0 + 721);
    const auto *osf0_722 = buffer.data(osf0 + 722);
    const auto *osf0_723 = buffer.data(osf0 + 723);
    const auto *osf0_724 = buffer.data(osf0 + 724);
    const auto *osf0_725 = buffer.data(osf0 + 725);
    const auto *osf0_726 = buffer.data(osf0 + 726);
    const auto *osf0_727 = buffer.data(osf0 + 727);
    const auto *osf0_728 = buffer.data(osf0 + 728);
    const auto *osf0_729 = buffer.data(osf0 + 729);
    const auto *osf0_730 = buffer.data(osf0 + 730);
    const auto *osf0_731 = buffer.data(osf0 + 731);

    const auto *osf1_679 = buffer.data(osf1 + 679);
    const auto *osf1_680 = buffer.data(osf1 + 680);
    const auto *osf1_681 = buffer.data(osf1 + 681);
    const auto *osf1_682 = buffer.data(osf1 + 682);
    const auto *osf1_683 = buffer.data(osf1 + 683);
    const auto *osf1_684 = buffer.data(osf1 + 684);
    const auto *osf1_685 = buffer.data(osf1 + 685);
    const auto *osf1_686 = buffer.data(osf1 + 686);
    const auto *osf1_687 = buffer.data(osf1 + 687);
    const auto *osf1_688 = buffer.data(osf1 + 688);
    const auto *osf1_689 = buffer.data(osf1 + 689);
    const auto *osf1_690 = buffer.data(osf1 + 690);
    const auto *osf1_691 = buffer.data(osf1 + 691);
    const auto *osf1_692 = buffer.data(osf1 + 692);
    const auto *osf1_693 = buffer.data(osf1 + 693);
    const auto *osf1_694 = buffer.data(osf1 + 694);
    const auto *osf1_695 = buffer.data(osf1 + 695);
    const auto *osf1_696 = buffer.data(osf1 + 696);
    const auto *osf1_697 = buffer.data(osf1 + 697);
    const auto *osf1_698 = buffer.data(osf1 + 698);
    const auto *osf1_699 = buffer.data(osf1 + 699);
    const auto *osf1_700 = buffer.data(osf1 + 700);
    const auto *osf1_701 = buffer.data(osf1 + 701);
    const auto *osf1_702 = buffer.data(osf1 + 702);
    const auto *osf1_703 = buffer.data(osf1 + 703);
    const auto *osf1_704 = buffer.data(osf1 + 704);
    const auto *osf1_705 = buffer.data(osf1 + 705);
    const auto *osf1_706 = buffer.data(osf1 + 706);
    const auto *osf1_707 = buffer.data(osf1 + 707);
    const auto *osf1_708 = buffer.data(osf1 + 708);
    const auto *osf1_709 = buffer.data(osf1 + 709);
    const auto *osf1_710 = buffer.data(osf1 + 710);
    const auto *osf1_711 = buffer.data(osf1 + 711);
    const auto *osf1_712 = buffer.data(osf1 + 712);
    const auto *osf1_713 = buffer.data(osf1 + 713);
    const auto *osf1_714 = buffer.data(osf1 + 714);
    const auto *osf1_715 = buffer.data(osf1 + 715);
    const auto *osf1_716 = buffer.data(osf1 + 716);
    const auto *osf1_717 = buffer.data(osf1 + 717);
    const auto *osf1_718 = buffer.data(osf1 + 718);
    const auto *osf1_719 = buffer.data(osf1 + 719);
    const auto *osf1_720 = buffer.data(osf1 + 720);
    const auto *osf1_721 = buffer.data(osf1 + 721);
    const auto *osf1_722 = buffer.data(osf1 + 722);
    const auto *osf1_723 = buffer.data(osf1 + 723);
    const auto *osf1_724 = buffer.data(osf1 + 724);
    const auto *osf1_725 = buffer.data(osf1 + 725);
    const auto *osf1_726 = buffer.data(osf1 + 726);
    const auto *osf1_727 = buffer.data(osf1 + 727);
    const auto *osf1_728 = buffer.data(osf1 + 728);
    const auto *osf1_729 = buffer.data(osf1 + 729);
    const auto *osf1_730 = buffer.data(osf1 + 730);
    const auto *osf1_731 = buffer.data(osf1 + 731);

    const auto *osg_1015 = buffer.data(osg + 1015);
    const auto *osg_1018 = buffer.data(osg + 1018);
    const auto *osg_1019 = buffer.data(osg + 1019);
    const auto *osg_1020 = buffer.data(osg + 1020);
    const auto *osg_1021 = buffer.data(osg + 1021);
    const auto *osg_1022 = buffer.data(osg + 1022);
    const auto *osg_1023 = buffer.data(osg + 1023);
    const auto *osg_1024 = buffer.data(osg + 1024);
    const auto *osg_1025 = buffer.data(osg + 1025);
    const auto *osg_1026 = buffer.data(osg + 1026);
    const auto *osg_1027 = buffer.data(osg + 1027);
    const auto *osg_1028 = buffer.data(osg + 1028);
    const auto *osg_1029 = buffer.data(osg + 1029);
    const auto *osg_1030 = buffer.data(osg + 1030);
    const auto *osg_1031 = buffer.data(osg + 1031);
    const auto *osg_1032 = buffer.data(osg + 1032);
    const auto *osg_1033 = buffer.data(osg + 1033);
    const auto *osg_1034 = buffer.data(osg + 1034);
    const auto *osg_1035 = buffer.data(osg + 1035);
    const auto *osg_1036 = buffer.data(osg + 1036);
    const auto *osg_1037 = buffer.data(osg + 1037);
    const auto *osg_1038 = buffer.data(osg + 1038);
    const auto *osg_1039 = buffer.data(osg + 1039);
    const auto *osg_1040 = buffer.data(osg + 1040);
    const auto *osg_1041 = buffer.data(osg + 1041);
    const auto *osg_1042 = buffer.data(osg + 1042);
    const auto *osg_1043 = buffer.data(osg + 1043);
    const auto *osg_1044 = buffer.data(osg + 1044);
    const auto *osg_1045 = buffer.data(osg + 1045);
    const auto *osg_1046 = buffer.data(osg + 1046);
    const auto *osg_1047 = buffer.data(osg + 1047);
    const auto *osg_1048 = buffer.data(osg + 1048);
    const auto *osg_1049 = buffer.data(osg + 1049);
    const auto *osg_1050 = buffer.data(osg + 1050);
    const auto *osg_1051 = buffer.data(osg + 1051);
    const auto *osg_1052 = buffer.data(osg + 1052);
    const auto *osg_1053 = buffer.data(osg + 1053);
    const auto *osg_1054 = buffer.data(osg + 1054);
    const auto *osg_1055 = buffer.data(osg + 1055);
    const auto *osg_1056 = buffer.data(osg + 1056);
    const auto *osg_1057 = buffer.data(osg + 1057);
    const auto *osg_1058 = buffer.data(osg + 1058);
    const auto *osg_1059 = buffer.data(osg + 1059);
    const auto *osg_1060 = buffer.data(osg + 1060);
    const auto *osg_1061 = buffer.data(osg + 1061);
    const auto *osg_1062 = buffer.data(osg + 1062);
    const auto *osg_1063 = buffer.data(osg + 1063);
    const auto *osg_1064 = buffer.data(osg + 1064);
    const auto *osg_1065 = buffer.data(osg + 1065);
    const auto *osg_1066 = buffer.data(osg + 1066);
    const auto *osg_1067 = buffer.data(osg + 1067);
    const auto *osg_1068 = buffer.data(osg + 1068);
    const auto *osg_1069 = buffer.data(osg + 1069);
    const auto *osg_1070 = buffer.data(osg + 1070);
    const auto *osg_1071 = buffer.data(osg + 1071);
    const auto *osg_1072 = buffer.data(osg + 1072);
    const auto *osg_1073 = buffer.data(osg + 1073);
    const auto *osg_1074 = buffer.data(osg + 1074);
    const auto *osg_1075 = buffer.data(osg + 1075);
    const auto *osg_1076 = buffer.data(osg + 1076);
    const auto *osg_1077 = buffer.data(osg + 1077);
    const auto *osg_1078 = buffer.data(osg + 1078);
    const auto *osg_1079 = buffer.data(osg + 1079);
    const auto *osg_1080 = buffer.data(osg + 1080);
    const auto *osg_1081 = buffer.data(osg + 1081);
    const auto *osg_1082 = buffer.data(osg + 1082);
    const auto *osg_1083 = buffer.data(osg + 1083);
    const auto *osg_1084 = buffer.data(osg + 1084);
    const auto *osg_1085 = buffer.data(osg + 1085);
    const auto *osg_1086 = buffer.data(osg + 1086);
    const auto *osg_1087 = buffer.data(osg + 1087);
    const auto *osg_1088 = buffer.data(osg + 1088);
    const auto *osg_1089 = buffer.data(osg + 1089);
    const auto *osg_1090 = buffer.data(osg + 1090);
    const auto *osg_1091 = buffer.data(osg + 1091);
    const auto *osg_1092 = buffer.data(osg + 1092);
    const auto *osg_1093 = buffer.data(osg + 1093);
    const auto *osg_1094 = buffer.data(osg + 1094);
    const auto *osg_1095 = buffer.data(osg + 1095);
    const auto *osg_1096 = buffer.data(osg + 1096);

#pragma omp simd aligned(t_1420, t_1421, t_1422, t_1423, pa_z, pc_x, pc_z, nsh0_1170, nsg_835, \
                         nsh1_1170, osg_1015, osg_1018, osg_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1420[k] = f_3 * pc_x[k] * osg_1018[k];

        t_1421[k] = f_3 * pc_x[k] * osg_1019[k];

        t_1422[k] = pa_z[k] * nsh0_1170[k]
                    - f_8 * pc_z[k] * nsh1_1170[k];

        t_1423[k] = f_9 * nsg_835[k]
                    + f_3 * pc_z[k] * osg_1015[k];
    }

#pragma omp simd aligned(t_1424, t_1425, t_1426, pa_z, pc_y, pc_z, nsh0_1172, nsh0_1173, \
                         nsg_836, nsg_837, nsg_854, nsh1_1172, nsh1_1173, \
                         osg_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1424[k] = pa_z[k] * nsh0_1172[k]
                    + f_10 * nsg_836[k]
                    - f_8 * pc_z[k] * nsh1_1172[k];

        t_1425[k] = pa_z[k] * nsh0_1173[k]
                    + f_11 * nsg_837[k]
                    - f_8 * pc_z[k] * nsh1_1173[k];

        t_1426[k] = f_12 * nsg_854[k]
                    + f_3 * pc_y[k] * osg_1019[k];
    }

#pragma omp simd aligned(t_1427, t_1428, t_1429, pc_x, pc_z, nsg_839, osf0_679, osf0_680, \
                         osf0_681, osf1_679, osf1_680, osf1_681, osg_1019, osg_1020, \
                         osg_1021 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1427[k] = f_9 * nsg_839[k]
                    + f_1 * osf0_679[k]
                    - f_2 * osf1_679[k]
                    + f_3 * pc_z[k] * osg_1019[k];

        t_1428[k] = f_1 * osf0_680[k]
                    - f_2 * osf1_680[k]
                    + f_3 * pc_x[k] * osg_1020[k];

        t_1429[k] = f_13 * osf0_681[k]
                    - f_14 * osf1_681[k]
                    + f_3 * pc_x[k] * osg_1021[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pc_x, osf0_682, osf0_683, osf0_684, osf1_682, \
                         osf1_683, osf1_684, osg_1022, osg_1023, \
                         osg_1024 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_13 * osf0_682[k]
                    - f_14 * osf1_682[k]
                    + f_3 * pc_x[k] * osg_1022[k];

        t_1431[k] = f_6 * osf0_683[k]
                    - f_7 * osf1_683[k]
                    + f_3 * pc_x[k] * osg_1023[k];

        t_1432[k] = f_6 * osf0_684[k]
                    - f_7 * osf1_684[k]
                    + f_3 * pc_x[k] * osg_1024[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pc_x, osf0_685, osf0_686, osf0_687, osf1_685, \
                         osf1_686, osf1_687, osg_1025, osg_1026, \
                         osg_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_6 * osf0_685[k]
                    - f_7 * osf1_685[k]
                    + f_3 * pc_x[k] * osg_1025[k];

        t_1434[k] = f_4 * osf0_686[k]
                    - f_5 * osf1_686[k]
                    + f_3 * pc_x[k] * osg_1026[k];

        t_1435[k] = f_4 * osf0_687[k]
                    - f_5 * osf1_687[k]
                    + f_3 * pc_x[k] * osg_1027[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, t_1439, t_1440, pc_x, osf0_688, osf0_689, \
                         osf1_688, osf1_689, osg_1028, osg_1029, osg_1030, osg_1031, \
                         osg_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_4 * osf0_688[k]
                    - f_5 * osf1_688[k]
                    + f_3 * pc_x[k] * osg_1028[k];

        t_1437[k] = f_4 * osf0_689[k]
                    - f_5 * osf1_689[k]
                    + f_3 * pc_x[k] * osg_1029[k];

        t_1438[k] = f_3 * pc_x[k] * osg_1030[k];

        t_1439[k] = f_3 * pc_x[k] * osg_1031[k];

        t_1440[k] = f_3 * pc_x[k] * osg_1032[k];
    }

#pragma omp simd aligned(t_1441, t_1442, t_1443, t_1444, pc_x, pc_y, pc_z, nsg_850, nsg_865, \
                         osf0_686, osf1_686, osg_1030, osg_1033, \
                         osg_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1441[k] = f_3 * pc_x[k] * osg_1033[k];

        t_1442[k] = f_3 * pc_x[k] * osg_1034[k];

        t_1443[k] = f_15 * nsg_865[k]
                    + f_1 * osf0_686[k]
                    - f_2 * osf1_686[k]
                    + f_3 * pc_y[k] * osg_1030[k];

        t_1444[k] = f_10 * nsg_850[k]
                    + f_3 * pc_z[k] * osg_1030[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_y, nsg_867, nsg_868, nsg_869, osf0_688, \
                         osf0_689, osf1_688, osf1_689, osg_1032, osg_1033, \
                         osg_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_15 * nsg_867[k]
                    + f_6 * osf0_688[k]
                    - f_7 * osf1_688[k]
                    + f_3 * pc_y[k] * osg_1032[k];

        t_1446[k] = f_15 * nsg_868[k]
                    + f_4 * osf0_689[k]
                    - f_5 * osf1_689[k]
                    + f_3 * pc_y[k] * osg_1033[k];

        t_1447[k] = f_15 * nsg_869[k]
                    + f_3 * pc_y[k] * osg_1034[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_z, nsg_854, osf0_689, osf0_690, \
                         osf0_691, osf1_689, osf1_690, osf1_691, osg_1034, osg_1035, \
                         osg_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_10 * nsg_854[k]
                    + f_1 * osf0_689[k]
                    - f_2 * osf1_689[k]
                    + f_3 * pc_z[k] * osg_1034[k];

        t_1449[k] = f_1 * osf0_690[k]
                    - f_2 * osf1_690[k]
                    + f_3 * pc_x[k] * osg_1035[k];

        t_1450[k] = f_13 * osf0_691[k]
                    - f_14 * osf1_691[k]
                    + f_3 * pc_x[k] * osg_1036[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, osf0_692, osf0_693, osf0_694, osf1_692, \
                         osf1_693, osf1_694, osg_1037, osg_1038, \
                         osg_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_13 * osf0_692[k]
                    - f_14 * osf1_692[k]
                    + f_3 * pc_x[k] * osg_1037[k];

        t_1452[k] = f_6 * osf0_693[k]
                    - f_7 * osf1_693[k]
                    + f_3 * pc_x[k] * osg_1038[k];

        t_1453[k] = f_6 * osf0_694[k]
                    - f_7 * osf1_694[k]
                    + f_3 * pc_x[k] * osg_1039[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, osf0_695, osf0_696, osf0_697, osf1_695, \
                         osf1_696, osf1_697, osg_1040, osg_1041, \
                         osg_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_6 * osf0_695[k]
                    - f_7 * osf1_695[k]
                    + f_3 * pc_x[k] * osg_1040[k];

        t_1455[k] = f_4 * osf0_696[k]
                    - f_5 * osf1_696[k]
                    + f_3 * pc_x[k] * osg_1041[k];

        t_1456[k] = f_4 * osf0_697[k]
                    - f_5 * osf1_697[k]
                    + f_3 * pc_x[k] * osg_1042[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, t_1460, t_1461, pc_x, osf0_698, osf0_699, \
                         osf1_698, osf1_699, osg_1043, osg_1044, osg_1045, osg_1046, \
                         osg_1047 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_4 * osf0_698[k]
                    - f_5 * osf1_698[k]
                    + f_3 * pc_x[k] * osg_1043[k];

        t_1458[k] = f_4 * osf0_699[k]
                    - f_5 * osf1_699[k]
                    + f_3 * pc_x[k] * osg_1044[k];

        t_1459[k] = f_3 * pc_x[k] * osg_1045[k];

        t_1460[k] = f_3 * pc_x[k] * osg_1046[k];

        t_1461[k] = f_3 * pc_x[k] * osg_1047[k];
    }

#pragma omp simd aligned(t_1462, t_1463, t_1464, t_1465, pc_x, pc_y, pc_z, nsg_865, nsg_880, \
                         osf0_696, osf1_696, osg_1045, osg_1048, \
                         osg_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1462[k] = f_3 * pc_x[k] * osg_1048[k];

        t_1463[k] = f_3 * pc_x[k] * osg_1049[k];

        t_1464[k] = f_16 * nsg_880[k]
                    + f_1 * osf0_696[k]
                    - f_2 * osf1_696[k]
                    + f_3 * pc_y[k] * osg_1045[k];

        t_1465[k] = f_11 * nsg_865[k]
                    + f_3 * pc_z[k] * osg_1045[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, pc_y, nsg_882, nsg_883, nsg_884, osf0_698, \
                         osf0_699, osf1_698, osf1_699, osg_1047, osg_1048, \
                         osg_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_16 * nsg_882[k]
                    + f_6 * osf0_698[k]
                    - f_7 * osf1_698[k]
                    + f_3 * pc_y[k] * osg_1047[k];

        t_1467[k] = f_16 * nsg_883[k]
                    + f_4 * osf0_699[k]
                    - f_5 * osf1_699[k]
                    + f_3 * pc_y[k] * osg_1048[k];

        t_1468[k] = f_16 * nsg_884[k]
                    + f_3 * pc_y[k] * osg_1049[k];
    }

#pragma omp simd aligned(t_1469, t_1470, t_1471, pc_x, pc_z, nsg_869, osf0_699, osf0_700, \
                         osf0_701, osf1_699, osf1_700, osf1_701, osg_1049, osg_1050, \
                         osg_1051 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1469[k] = f_11 * nsg_869[k]
                    + f_1 * osf0_699[k]
                    - f_2 * osf1_699[k]
                    + f_3 * pc_z[k] * osg_1049[k];

        t_1470[k] = f_1 * osf0_700[k]
                    - f_2 * osf1_700[k]
                    + f_3 * pc_x[k] * osg_1050[k];

        t_1471[k] = f_13 * osf0_701[k]
                    - f_14 * osf1_701[k]
                    + f_3 * pc_x[k] * osg_1051[k];
    }

#pragma omp simd aligned(t_1472, t_1473, t_1474, pc_x, osf0_702, osf0_703, osf0_704, osf1_702, \
                         osf1_703, osf1_704, osg_1052, osg_1053, \
                         osg_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1472[k] = f_13 * osf0_702[k]
                    - f_14 * osf1_702[k]
                    + f_3 * pc_x[k] * osg_1052[k];

        t_1473[k] = f_6 * osf0_703[k]
                    - f_7 * osf1_703[k]
                    + f_3 * pc_x[k] * osg_1053[k];

        t_1474[k] = f_6 * osf0_704[k]
                    - f_7 * osf1_704[k]
                    + f_3 * pc_x[k] * osg_1054[k];
    }

#pragma omp simd aligned(t_1475, t_1476, t_1477, pc_x, osf0_705, osf0_706, osf0_707, osf1_705, \
                         osf1_706, osf1_707, osg_1055, osg_1056, \
                         osg_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1475[k] = f_6 * osf0_705[k]
                    - f_7 * osf1_705[k]
                    + f_3 * pc_x[k] * osg_1055[k];

        t_1476[k] = f_4 * osf0_706[k]
                    - f_5 * osf1_706[k]
                    + f_3 * pc_x[k] * osg_1056[k];

        t_1477[k] = f_4 * osf0_707[k]
                    - f_5 * osf1_707[k]
                    + f_3 * pc_x[k] * osg_1057[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, t_1481, t_1482, pc_x, osf0_708, osf0_709, \
                         osf1_708, osf1_709, osg_1058, osg_1059, osg_1060, osg_1061, \
                         osg_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_4 * osf0_708[k]
                    - f_5 * osf1_708[k]
                    + f_3 * pc_x[k] * osg_1058[k];

        t_1479[k] = f_4 * osf0_709[k]
                    - f_5 * osf1_709[k]
                    + f_3 * pc_x[k] * osg_1059[k];

        t_1480[k] = f_3 * pc_x[k] * osg_1060[k];

        t_1481[k] = f_3 * pc_x[k] * osg_1061[k];

        t_1482[k] = f_3 * pc_x[k] * osg_1062[k];
    }

#pragma omp simd aligned(t_1483, t_1484, t_1485, t_1486, pc_x, pc_y, pc_z, nsg_880, nsg_895, \
                         osf0_706, osf1_706, osg_1060, osg_1063, \
                         osg_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1483[k] = f_3 * pc_x[k] * osg_1063[k];

        t_1484[k] = f_3 * pc_x[k] * osg_1064[k];

        t_1485[k] = f_17 * nsg_895[k]
                    + f_1 * osf0_706[k]
                    - f_2 * osf1_706[k]
                    + f_3 * pc_y[k] * osg_1060[k];

        t_1486[k] = f_18 * nsg_880[k]
                    + f_3 * pc_z[k] * osg_1060[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pc_y, nsg_897, nsg_898, nsg_899, osf0_708, \
                         osf0_709, osf1_708, osf1_709, osg_1062, osg_1063, \
                         osg_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_17 * nsg_897[k]
                    + f_6 * osf0_708[k]
                    - f_7 * osf1_708[k]
                    + f_3 * pc_y[k] * osg_1062[k];

        t_1488[k] = f_17 * nsg_898[k]
                    + f_4 * osf0_709[k]
                    - f_5 * osf1_709[k]
                    + f_3 * pc_y[k] * osg_1063[k];

        t_1489[k] = f_17 * nsg_899[k]
                    + f_3 * pc_y[k] * osg_1064[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pc_x, pc_z, nsg_884, osf0_709, osf0_710, \
                         osf0_711, osf1_709, osf1_710, osf1_711, osg_1064, osg_1065, \
                         osg_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_18 * nsg_884[k]
                    + f_1 * osf0_709[k]
                    - f_2 * osf1_709[k]
                    + f_3 * pc_z[k] * osg_1064[k];

        t_1491[k] = f_1 * osf0_710[k]
                    - f_2 * osf1_710[k]
                    + f_3 * pc_x[k] * osg_1065[k];

        t_1492[k] = f_13 * osf0_711[k]
                    - f_14 * osf1_711[k]
                    + f_3 * pc_x[k] * osg_1066[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pc_x, osf0_712, osf0_713, osf0_714, osf1_712, \
                         osf1_713, osf1_714, osg_1067, osg_1068, \
                         osg_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_13 * osf0_712[k]
                    - f_14 * osf1_712[k]
                    + f_3 * pc_x[k] * osg_1067[k];

        t_1494[k] = f_6 * osf0_713[k]
                    - f_7 * osf1_713[k]
                    + f_3 * pc_x[k] * osg_1068[k];

        t_1495[k] = f_6 * osf0_714[k]
                    - f_7 * osf1_714[k]
                    + f_3 * pc_x[k] * osg_1069[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, pc_x, osf0_715, osf0_716, osf0_717, osf1_715, \
                         osf1_716, osf1_717, osg_1070, osg_1071, \
                         osg_1072 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_6 * osf0_715[k]
                    - f_7 * osf1_715[k]
                    + f_3 * pc_x[k] * osg_1070[k];

        t_1497[k] = f_4 * osf0_716[k]
                    - f_5 * osf1_716[k]
                    + f_3 * pc_x[k] * osg_1071[k];

        t_1498[k] = f_4 * osf0_717[k]
                    - f_5 * osf1_717[k]
                    + f_3 * pc_x[k] * osg_1072[k];
    }

#pragma omp simd aligned(t_1499, t_1500, t_1501, t_1502, t_1503, pc_x, osf0_718, osf0_719, \
                         osf1_718, osf1_719, osg_1073, osg_1074, osg_1075, osg_1076, \
                         osg_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1499[k] = f_4 * osf0_718[k]
                    - f_5 * osf1_718[k]
                    + f_3 * pc_x[k] * osg_1073[k];

        t_1500[k] = f_4 * osf0_719[k]
                    - f_5 * osf1_719[k]
                    + f_3 * pc_x[k] * osg_1074[k];

        t_1501[k] = f_3 * pc_x[k] * osg_1075[k];

        t_1502[k] = f_3 * pc_x[k] * osg_1076[k];

        t_1503[k] = f_3 * pc_x[k] * osg_1077[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, t_1507, pc_x, pc_y, pc_z, nsg_895, nsg_910, \
                         osf0_716, osf1_716, osg_1075, osg_1078, \
                         osg_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = f_3 * pc_x[k] * osg_1078[k];

        t_1505[k] = f_3 * pc_x[k] * osg_1079[k];

        t_1506[k] = f_19 * nsg_910[k]
                    + f_1 * osf0_716[k]
                    - f_2 * osf1_716[k]
                    + f_3 * pc_y[k] * osg_1075[k];

        t_1507[k] = f_20 * nsg_895[k]
                    + f_3 * pc_z[k] * osg_1075[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pc_y, nsg_912, nsg_913, nsg_914, osf0_718, \
                         osf0_719, osf1_718, osf1_719, osg_1077, osg_1078, \
                         osg_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_19 * nsg_912[k]
                    + f_6 * osf0_718[k]
                    - f_7 * osf1_718[k]
                    + f_3 * pc_y[k] * osg_1077[k];

        t_1509[k] = f_19 * nsg_913[k]
                    + f_4 * osf0_719[k]
                    - f_5 * osf1_719[k]
                    + f_3 * pc_y[k] * osg_1078[k];

        t_1510[k] = f_19 * nsg_914[k]
                    + f_3 * pc_y[k] * osg_1079[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, pc_x, pc_z, nsg_899, osf0_719, osf0_720, \
                         osf0_721, osf1_719, osf1_720, osf1_721, osg_1079, osg_1080, \
                         osg_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = f_20 * nsg_899[k]
                    + f_1 * osf0_719[k]
                    - f_2 * osf1_719[k]
                    + f_3 * pc_z[k] * osg_1079[k];

        t_1512[k] = f_1 * osf0_720[k]
                    - f_2 * osf1_720[k]
                    + f_3 * pc_x[k] * osg_1080[k];

        t_1513[k] = f_13 * osf0_721[k]
                    - f_14 * osf1_721[k]
                    + f_3 * pc_x[k] * osg_1081[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, pc_x, osf0_722, osf0_723, osf0_724, osf1_722, \
                         osf1_723, osf1_724, osg_1082, osg_1083, \
                         osg_1084 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_13 * osf0_722[k]
                    - f_14 * osf1_722[k]
                    + f_3 * pc_x[k] * osg_1082[k];

        t_1515[k] = f_6 * osf0_723[k]
                    - f_7 * osf1_723[k]
                    + f_3 * pc_x[k] * osg_1083[k];

        t_1516[k] = f_6 * osf0_724[k]
                    - f_7 * osf1_724[k]
                    + f_3 * pc_x[k] * osg_1084[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, pc_x, osf0_725, osf0_726, osf0_727, osf1_725, \
                         osf1_726, osf1_727, osg_1085, osg_1086, \
                         osg_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_6 * osf0_725[k]
                    - f_7 * osf1_725[k]
                    + f_3 * pc_x[k] * osg_1085[k];

        t_1518[k] = f_4 * osf0_726[k]
                    - f_5 * osf1_726[k]
                    + f_3 * pc_x[k] * osg_1086[k];

        t_1519[k] = f_4 * osf0_727[k]
                    - f_5 * osf1_727[k]
                    + f_3 * pc_x[k] * osg_1087[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, t_1523, t_1524, pc_x, osf0_728, osf0_729, \
                         osf1_728, osf1_729, osg_1088, osg_1089, osg_1090, osg_1091, \
                         osg_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_4 * osf0_728[k]
                    - f_5 * osf1_728[k]
                    + f_3 * pc_x[k] * osg_1088[k];

        t_1521[k] = f_4 * osf0_729[k]
                    - f_5 * osf1_729[k]
                    + f_3 * pc_x[k] * osg_1089[k];

        t_1522[k] = f_3 * pc_x[k] * osg_1090[k];

        t_1523[k] = f_3 * pc_x[k] * osg_1091[k];

        t_1524[k] = f_3 * pc_x[k] * osg_1092[k];
    }

#pragma omp simd aligned(t_1525, t_1526, t_1527, t_1528, pc_x, pc_y, pc_z, nsg_910, nsg_925, \
                         osf0_726, osf1_726, osg_1090, osg_1093, \
                         osg_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1525[k] = f_3 * pc_x[k] * osg_1093[k];

        t_1526[k] = f_3 * pc_x[k] * osg_1094[k];

        t_1527[k] = f_20 * nsg_925[k]
                    + f_1 * osf0_726[k]
                    - f_2 * osf1_726[k]
                    + f_3 * pc_y[k] * osg_1090[k];

        t_1528[k] = f_19 * nsg_910[k]
                    + f_3 * pc_z[k] * osg_1090[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, pc_y, nsg_927, nsg_928, nsg_929, osf0_728, \
                         osf0_729, osf1_728, osf1_729, osg_1092, osg_1093, \
                         osg_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_20 * nsg_927[k]
                    + f_6 * osf0_728[k]
                    - f_7 * osf1_728[k]
                    + f_3 * pc_y[k] * osg_1092[k];

        t_1530[k] = f_20 * nsg_928[k]
                    + f_4 * osf0_729[k]
                    - f_5 * osf1_729[k]
                    + f_3 * pc_y[k] * osg_1093[k];

        t_1531[k] = f_20 * nsg_929[k]
                    + f_3 * pc_y[k] * osg_1094[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, pc_x, pc_z, nsg_914, osf0_729, osf0_730, \
                         osf0_731, osf1_729, osf1_730, osf1_731, osg_1094, osg_1095, \
                         osg_1096 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = f_19 * nsg_914[k]
                    + f_1 * osf0_729[k]
                    - f_2 * osf1_729[k]
                    + f_3 * pc_z[k] * osg_1094[k];

        t_1533[k] = f_1 * osf0_730[k]
                    - f_2 * osf1_730[k]
                    + f_3 * pc_x[k] * osg_1095[k];

        t_1534[k] = f_13 * osf0_731[k]
                    - f_14 * osf1_731[k]
                    + f_3 * pc_x[k] * osg_1096[k];
    }
}

static auto
compute_prim_osh_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t nsh0,
                                                           const size_t nsg, const size_t nsh1,
                                                           const size_t osf0, const size_t osf1,
                                                           const size_t osg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsh0_1365 = buffer.data(nsh0 + 1365);
    const auto *nsh0_1367 = buffer.data(nsh0 + 1367);
    const auto *nsh0_1370 = buffer.data(nsh0 + 1370);
    const auto *nsh0_1374 = buffer.data(nsh0 + 1374);
    const auto *nsh0_1380 = buffer.data(nsh0 + 1380);
    const auto *nsh0_1382 = buffer.data(nsh0 + 1382);
    const auto *nsh0_1383 = buffer.data(nsh0 + 1383);
    const auto *nsh0_1385 = buffer.data(nsh0 + 1385);

    const auto *nsg_925 = buffer.data(nsg + 925);
    const auto *nsg_929 = buffer.data(nsg + 929);
    const auto *nsg_940 = buffer.data(nsg + 940);
    const auto *nsg_942 = buffer.data(nsg + 942);
    const auto *nsg_943 = buffer.data(nsg + 943);
    const auto *nsg_944 = buffer.data(nsg + 944);
    const auto *nsg_955 = buffer.data(nsg + 955);
    const auto *nsg_957 = buffer.data(nsg + 957);
    const auto *nsg_958 = buffer.data(nsg + 958);
    const auto *nsg_959 = buffer.data(nsg + 959);
    const auto *nsg_970 = buffer.data(nsg + 970);
    const auto *nsg_972 = buffer.data(nsg + 972);
    const auto *nsg_973 = buffer.data(nsg + 973);
    const auto *nsg_974 = buffer.data(nsg + 974);
    const auto *nsg_985 = buffer.data(nsg + 985);
    const auto *nsg_987 = buffer.data(nsg + 987);
    const auto *nsg_988 = buffer.data(nsg + 988);
    const auto *nsg_989 = buffer.data(nsg + 989);

    const auto *nsh1_1365 = buffer.data(nsh1 + 1365);
    const auto *nsh1_1367 = buffer.data(nsh1 + 1367);
    const auto *nsh1_1370 = buffer.data(nsh1 + 1370);
    const auto *nsh1_1374 = buffer.data(nsh1 + 1374);
    const auto *nsh1_1380 = buffer.data(nsh1 + 1380);
    const auto *nsh1_1382 = buffer.data(nsh1 + 1382);
    const auto *nsh1_1383 = buffer.data(nsh1 + 1383);
    const auto *nsh1_1385 = buffer.data(nsh1 + 1385);

    const auto *osf0_732 = buffer.data(osf0 + 732);
    const auto *osf0_733 = buffer.data(osf0 + 733);
    const auto *osf0_734 = buffer.data(osf0 + 734);
    const auto *osf0_735 = buffer.data(osf0 + 735);
    const auto *osf0_736 = buffer.data(osf0 + 736);
    const auto *osf0_737 = buffer.data(osf0 + 737);
    const auto *osf0_738 = buffer.data(osf0 + 738);
    const auto *osf0_739 = buffer.data(osf0 + 739);
    const auto *osf0_740 = buffer.data(osf0 + 740);
    const auto *osf0_741 = buffer.data(osf0 + 741);
    const auto *osf0_742 = buffer.data(osf0 + 742);
    const auto *osf0_743 = buffer.data(osf0 + 743);
    const auto *osf0_744 = buffer.data(osf0 + 744);
    const auto *osf0_745 = buffer.data(osf0 + 745);
    const auto *osf0_746 = buffer.data(osf0 + 746);
    const auto *osf0_747 = buffer.data(osf0 + 747);
    const auto *osf0_748 = buffer.data(osf0 + 748);
    const auto *osf0_749 = buffer.data(osf0 + 749);
    const auto *osf0_750 = buffer.data(osf0 + 750);
    const auto *osf0_751 = buffer.data(osf0 + 751);
    const auto *osf0_752 = buffer.data(osf0 + 752);
    const auto *osf0_753 = buffer.data(osf0 + 753);
    const auto *osf0_754 = buffer.data(osf0 + 754);
    const auto *osf0_755 = buffer.data(osf0 + 755);
    const auto *osf0_756 = buffer.data(osf0 + 756);
    const auto *osf0_757 = buffer.data(osf0 + 757);
    const auto *osf0_758 = buffer.data(osf0 + 758);
    const auto *osf0_759 = buffer.data(osf0 + 759);
    const auto *osf0_761 = buffer.data(osf0 + 761);
    const auto *osf0_763 = buffer.data(osf0 + 763);
    const auto *osf0_764 = buffer.data(osf0 + 764);
    const auto *osf0_766 = buffer.data(osf0 + 766);
    const auto *osf0_767 = buffer.data(osf0 + 767);
    const auto *osf0_768 = buffer.data(osf0 + 768);
    const auto *osf0_770 = buffer.data(osf0 + 770);
    const auto *osf0_772 = buffer.data(osf0 + 772);
    const auto *osf0_773 = buffer.data(osf0 + 773);
    const auto *osf0_775 = buffer.data(osf0 + 775);
    const auto *osf0_776 = buffer.data(osf0 + 776);
    const auto *osf0_777 = buffer.data(osf0 + 777);
    const auto *osf0_778 = buffer.data(osf0 + 778);
    const auto *osf0_779 = buffer.data(osf0 + 779);

    const auto *osf1_732 = buffer.data(osf1 + 732);
    const auto *osf1_733 = buffer.data(osf1 + 733);
    const auto *osf1_734 = buffer.data(osf1 + 734);
    const auto *osf1_735 = buffer.data(osf1 + 735);
    const auto *osf1_736 = buffer.data(osf1 + 736);
    const auto *osf1_737 = buffer.data(osf1 + 737);
    const auto *osf1_738 = buffer.data(osf1 + 738);
    const auto *osf1_739 = buffer.data(osf1 + 739);
    const auto *osf1_740 = buffer.data(osf1 + 740);
    const auto *osf1_741 = buffer.data(osf1 + 741);
    const auto *osf1_742 = buffer.data(osf1 + 742);
    const auto *osf1_743 = buffer.data(osf1 + 743);
    const auto *osf1_744 = buffer.data(osf1 + 744);
    const auto *osf1_745 = buffer.data(osf1 + 745);
    const auto *osf1_746 = buffer.data(osf1 + 746);
    const auto *osf1_747 = buffer.data(osf1 + 747);
    const auto *osf1_748 = buffer.data(osf1 + 748);
    const auto *osf1_749 = buffer.data(osf1 + 749);
    const auto *osf1_750 = buffer.data(osf1 + 750);
    const auto *osf1_751 = buffer.data(osf1 + 751);
    const auto *osf1_752 = buffer.data(osf1 + 752);
    const auto *osf1_753 = buffer.data(osf1 + 753);
    const auto *osf1_754 = buffer.data(osf1 + 754);
    const auto *osf1_755 = buffer.data(osf1 + 755);
    const auto *osf1_756 = buffer.data(osf1 + 756);
    const auto *osf1_757 = buffer.data(osf1 + 757);
    const auto *osf1_758 = buffer.data(osf1 + 758);
    const auto *osf1_759 = buffer.data(osf1 + 759);
    const auto *osf1_761 = buffer.data(osf1 + 761);
    const auto *osf1_763 = buffer.data(osf1 + 763);
    const auto *osf1_764 = buffer.data(osf1 + 764);
    const auto *osf1_766 = buffer.data(osf1 + 766);
    const auto *osf1_767 = buffer.data(osf1 + 767);
    const auto *osf1_768 = buffer.data(osf1 + 768);
    const auto *osf1_770 = buffer.data(osf1 + 770);
    const auto *osf1_772 = buffer.data(osf1 + 772);
    const auto *osf1_773 = buffer.data(osf1 + 773);
    const auto *osf1_775 = buffer.data(osf1 + 775);
    const auto *osf1_776 = buffer.data(osf1 + 776);
    const auto *osf1_777 = buffer.data(osf1 + 777);
    const auto *osf1_778 = buffer.data(osf1 + 778);
    const auto *osf1_779 = buffer.data(osf1 + 779);

    const auto *osg_1097 = buffer.data(osg + 1097);
    const auto *osg_1098 = buffer.data(osg + 1098);
    const auto *osg_1099 = buffer.data(osg + 1099);
    const auto *osg_1100 = buffer.data(osg + 1100);
    const auto *osg_1101 = buffer.data(osg + 1101);
    const auto *osg_1102 = buffer.data(osg + 1102);
    const auto *osg_1103 = buffer.data(osg + 1103);
    const auto *osg_1104 = buffer.data(osg + 1104);
    const auto *osg_1105 = buffer.data(osg + 1105);
    const auto *osg_1106 = buffer.data(osg + 1106);
    const auto *osg_1107 = buffer.data(osg + 1107);
    const auto *osg_1108 = buffer.data(osg + 1108);
    const auto *osg_1109 = buffer.data(osg + 1109);
    const auto *osg_1110 = buffer.data(osg + 1110);
    const auto *osg_1111 = buffer.data(osg + 1111);
    const auto *osg_1112 = buffer.data(osg + 1112);
    const auto *osg_1113 = buffer.data(osg + 1113);
    const auto *osg_1114 = buffer.data(osg + 1114);
    const auto *osg_1115 = buffer.data(osg + 1115);
    const auto *osg_1116 = buffer.data(osg + 1116);
    const auto *osg_1117 = buffer.data(osg + 1117);
    const auto *osg_1118 = buffer.data(osg + 1118);
    const auto *osg_1119 = buffer.data(osg + 1119);
    const auto *osg_1120 = buffer.data(osg + 1120);
    const auto *osg_1121 = buffer.data(osg + 1121);
    const auto *osg_1122 = buffer.data(osg + 1122);
    const auto *osg_1123 = buffer.data(osg + 1123);
    const auto *osg_1124 = buffer.data(osg + 1124);
    const auto *osg_1125 = buffer.data(osg + 1125);
    const auto *osg_1126 = buffer.data(osg + 1126);
    const auto *osg_1127 = buffer.data(osg + 1127);
    const auto *osg_1128 = buffer.data(osg + 1128);
    const auto *osg_1129 = buffer.data(osg + 1129);
    const auto *osg_1130 = buffer.data(osg + 1130);
    const auto *osg_1131 = buffer.data(osg + 1131);
    const auto *osg_1132 = buffer.data(osg + 1132);
    const auto *osg_1133 = buffer.data(osg + 1133);
    const auto *osg_1134 = buffer.data(osg + 1134);
    const auto *osg_1135 = buffer.data(osg + 1135);
    const auto *osg_1136 = buffer.data(osg + 1136);
    const auto *osg_1137 = buffer.data(osg + 1137);
    const auto *osg_1138 = buffer.data(osg + 1138);
    const auto *osg_1139 = buffer.data(osg + 1139);
    const auto *osg_1141 = buffer.data(osg + 1141);
    const auto *osg_1143 = buffer.data(osg + 1143);
    const auto *osg_1144 = buffer.data(osg + 1144);
    const auto *osg_1146 = buffer.data(osg + 1146);
    const auto *osg_1147 = buffer.data(osg + 1147);
    const auto *osg_1148 = buffer.data(osg + 1148);
    const auto *osg_1150 = buffer.data(osg + 1150);
    const auto *osg_1151 = buffer.data(osg + 1151);
    const auto *osg_1152 = buffer.data(osg + 1152);
    const auto *osg_1153 = buffer.data(osg + 1153);
    const auto *osg_1154 = buffer.data(osg + 1154);
    const auto *osg_1155 = buffer.data(osg + 1155);
    const auto *osg_1157 = buffer.data(osg + 1157);
    const auto *osg_1158 = buffer.data(osg + 1158);
    const auto *osg_1160 = buffer.data(osg + 1160);
    const auto *osg_1161 = buffer.data(osg + 1161);
    const auto *osg_1162 = buffer.data(osg + 1162);
    const auto *osg_1164 = buffer.data(osg + 1164);
    const auto *osg_1165 = buffer.data(osg + 1165);
    const auto *osg_1166 = buffer.data(osg + 1166);
    const auto *osg_1167 = buffer.data(osg + 1167);
    const auto *osg_1168 = buffer.data(osg + 1168);
    const auto *osg_1169 = buffer.data(osg + 1169);

#pragma omp simd aligned(t_1535, t_1536, t_1537, pc_x, osf0_732, osf0_733, osf0_734, osf1_732, \
                         osf1_733, osf1_734, osg_1097, osg_1098, \
                         osg_1099 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1535[k] = f_13 * osf0_732[k]
                    - f_14 * osf1_732[k]
                    + f_3 * pc_x[k] * osg_1097[k];

        t_1536[k] = f_6 * osf0_733[k]
                    - f_7 * osf1_733[k]
                    + f_3 * pc_x[k] * osg_1098[k];

        t_1537[k] = f_6 * osf0_734[k]
                    - f_7 * osf1_734[k]
                    + f_3 * pc_x[k] * osg_1099[k];
    }

#pragma omp simd aligned(t_1538, t_1539, t_1540, pc_x, osf0_735, osf0_736, osf0_737, osf1_735, \
                         osf1_736, osf1_737, osg_1100, osg_1101, \
                         osg_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1538[k] = f_6 * osf0_735[k]
                    - f_7 * osf1_735[k]
                    + f_3 * pc_x[k] * osg_1100[k];

        t_1539[k] = f_4 * osf0_736[k]
                    - f_5 * osf1_736[k]
                    + f_3 * pc_x[k] * osg_1101[k];

        t_1540[k] = f_4 * osf0_737[k]
                    - f_5 * osf1_737[k]
                    + f_3 * pc_x[k] * osg_1102[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, t_1544, t_1545, pc_x, osf0_738, osf0_739, \
                         osf1_738, osf1_739, osg_1103, osg_1104, osg_1105, osg_1106, \
                         osg_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_4 * osf0_738[k]
                    - f_5 * osf1_738[k]
                    + f_3 * pc_x[k] * osg_1103[k];

        t_1542[k] = f_4 * osf0_739[k]
                    - f_5 * osf1_739[k]
                    + f_3 * pc_x[k] * osg_1104[k];

        t_1543[k] = f_3 * pc_x[k] * osg_1105[k];

        t_1544[k] = f_3 * pc_x[k] * osg_1106[k];

        t_1545[k] = f_3 * pc_x[k] * osg_1107[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pc_x, pc_y, pc_z, nsg_925, nsg_940, \
                         osf0_736, osf1_736, osg_1105, osg_1108, \
                         osg_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_3 * pc_x[k] * osg_1108[k];

        t_1547[k] = f_3 * pc_x[k] * osg_1109[k];

        t_1548[k] = f_18 * nsg_940[k]
                    + f_1 * osf0_736[k]
                    - f_2 * osf1_736[k]
                    + f_3 * pc_y[k] * osg_1105[k];

        t_1549[k] = f_17 * nsg_925[k]
                    + f_3 * pc_z[k] * osg_1105[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pc_y, nsg_942, nsg_943, nsg_944, osf0_738, \
                         osf0_739, osf1_738, osf1_739, osg_1107, osg_1108, \
                         osg_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_18 * nsg_942[k]
                    + f_6 * osf0_738[k]
                    - f_7 * osf1_738[k]
                    + f_3 * pc_y[k] * osg_1107[k];

        t_1551[k] = f_18 * nsg_943[k]
                    + f_4 * osf0_739[k]
                    - f_5 * osf1_739[k]
                    + f_3 * pc_y[k] * osg_1108[k];

        t_1552[k] = f_18 * nsg_944[k]
                    + f_3 * pc_y[k] * osg_1109[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, pc_x, pc_z, nsg_929, osf0_739, osf0_740, \
                         osf0_741, osf1_739, osf1_740, osf1_741, osg_1109, osg_1110, \
                         osg_1111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = f_17 * nsg_929[k]
                    + f_1 * osf0_739[k]
                    - f_2 * osf1_739[k]
                    + f_3 * pc_z[k] * osg_1109[k];

        t_1554[k] = f_1 * osf0_740[k]
                    - f_2 * osf1_740[k]
                    + f_3 * pc_x[k] * osg_1110[k];

        t_1555[k] = f_13 * osf0_741[k]
                    - f_14 * osf1_741[k]
                    + f_3 * pc_x[k] * osg_1111[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, pc_x, osf0_742, osf0_743, osf0_744, osf1_742, \
                         osf1_743, osf1_744, osg_1112, osg_1113, \
                         osg_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_13 * osf0_742[k]
                    - f_14 * osf1_742[k]
                    + f_3 * pc_x[k] * osg_1112[k];

        t_1557[k] = f_6 * osf0_743[k]
                    - f_7 * osf1_743[k]
                    + f_3 * pc_x[k] * osg_1113[k];

        t_1558[k] = f_6 * osf0_744[k]
                    - f_7 * osf1_744[k]
                    + f_3 * pc_x[k] * osg_1114[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, pc_x, osf0_745, osf0_746, osf0_747, osf1_745, \
                         osf1_746, osf1_747, osg_1115, osg_1116, \
                         osg_1117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_6 * osf0_745[k]
                    - f_7 * osf1_745[k]
                    + f_3 * pc_x[k] * osg_1115[k];

        t_1560[k] = f_4 * osf0_746[k]
                    - f_5 * osf1_746[k]
                    + f_3 * pc_x[k] * osg_1116[k];

        t_1561[k] = f_4 * osf0_747[k]
                    - f_5 * osf1_747[k]
                    + f_3 * pc_x[k] * osg_1117[k];
    }

#pragma omp simd aligned(t_1562, t_1563, t_1564, t_1565, t_1566, pc_x, osf0_748, osf0_749, \
                         osf1_748, osf1_749, osg_1118, osg_1119, osg_1120, osg_1121, \
                         osg_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1562[k] = f_4 * osf0_748[k]
                    - f_5 * osf1_748[k]
                    + f_3 * pc_x[k] * osg_1118[k];

        t_1563[k] = f_4 * osf0_749[k]
                    - f_5 * osf1_749[k]
                    + f_3 * pc_x[k] * osg_1119[k];

        t_1564[k] = f_3 * pc_x[k] * osg_1120[k];

        t_1565[k] = f_3 * pc_x[k] * osg_1121[k];

        t_1566[k] = f_3 * pc_x[k] * osg_1122[k];
    }

#pragma omp simd aligned(t_1567, t_1568, t_1569, t_1570, pc_x, pc_y, pc_z, nsg_940, nsg_955, \
                         osf0_746, osf1_746, osg_1120, osg_1123, \
                         osg_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1567[k] = f_3 * pc_x[k] * osg_1123[k];

        t_1568[k] = f_3 * pc_x[k] * osg_1124[k];

        t_1569[k] = f_11 * nsg_955[k]
                    + f_1 * osf0_746[k]
                    - f_2 * osf1_746[k]
                    + f_3 * pc_y[k] * osg_1120[k];

        t_1570[k] = f_16 * nsg_940[k]
                    + f_3 * pc_z[k] * osg_1120[k];
    }

#pragma omp simd aligned(t_1571, t_1572, t_1573, pc_y, nsg_957, nsg_958, nsg_959, osf0_748, \
                         osf0_749, osf1_748, osf1_749, osg_1122, osg_1123, \
                         osg_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1571[k] = f_11 * nsg_957[k]
                    + f_6 * osf0_748[k]
                    - f_7 * osf1_748[k]
                    + f_3 * pc_y[k] * osg_1122[k];

        t_1572[k] = f_11 * nsg_958[k]
                    + f_4 * osf0_749[k]
                    - f_5 * osf1_749[k]
                    + f_3 * pc_y[k] * osg_1123[k];

        t_1573[k] = f_11 * nsg_959[k]
                    + f_3 * pc_y[k] * osg_1124[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, pc_x, pc_z, nsg_944, osf0_749, osf0_750, \
                         osf0_751, osf1_749, osf1_750, osf1_751, osg_1124, osg_1125, \
                         osg_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = f_16 * nsg_944[k]
                    + f_1 * osf0_749[k]
                    - f_2 * osf1_749[k]
                    + f_3 * pc_z[k] * osg_1124[k];

        t_1575[k] = f_1 * osf0_750[k]
                    - f_2 * osf1_750[k]
                    + f_3 * pc_x[k] * osg_1125[k];

        t_1576[k] = f_13 * osf0_751[k]
                    - f_14 * osf1_751[k]
                    + f_3 * pc_x[k] * osg_1126[k];
    }

#pragma omp simd aligned(t_1577, t_1578, t_1579, pc_x, osf0_752, osf0_753, osf0_754, osf1_752, \
                         osf1_753, osf1_754, osg_1127, osg_1128, \
                         osg_1129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1577[k] = f_13 * osf0_752[k]
                    - f_14 * osf1_752[k]
                    + f_3 * pc_x[k] * osg_1127[k];

        t_1578[k] = f_6 * osf0_753[k]
                    - f_7 * osf1_753[k]
                    + f_3 * pc_x[k] * osg_1128[k];

        t_1579[k] = f_6 * osf0_754[k]
                    - f_7 * osf1_754[k]
                    + f_3 * pc_x[k] * osg_1129[k];
    }

#pragma omp simd aligned(t_1580, t_1581, t_1582, pc_x, osf0_755, osf0_756, osf0_757, osf1_755, \
                         osf1_756, osf1_757, osg_1130, osg_1131, \
                         osg_1132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1580[k] = f_6 * osf0_755[k]
                    - f_7 * osf1_755[k]
                    + f_3 * pc_x[k] * osg_1130[k];

        t_1581[k] = f_4 * osf0_756[k]
                    - f_5 * osf1_756[k]
                    + f_3 * pc_x[k] * osg_1131[k];

        t_1582[k] = f_4 * osf0_757[k]
                    - f_5 * osf1_757[k]
                    + f_3 * pc_x[k] * osg_1132[k];
    }

#pragma omp simd aligned(t_1583, t_1584, t_1585, t_1586, t_1587, pc_x, osf0_758, osf0_759, \
                         osf1_758, osf1_759, osg_1133, osg_1134, osg_1135, osg_1136, \
                         osg_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1583[k] = f_4 * osf0_758[k]
                    - f_5 * osf1_758[k]
                    + f_3 * pc_x[k] * osg_1133[k];

        t_1584[k] = f_4 * osf0_759[k]
                    - f_5 * osf1_759[k]
                    + f_3 * pc_x[k] * osg_1134[k];

        t_1585[k] = f_3 * pc_x[k] * osg_1135[k];

        t_1586[k] = f_3 * pc_x[k] * osg_1136[k];

        t_1587[k] = f_3 * pc_x[k] * osg_1137[k];
    }

#pragma omp simd aligned(t_1588, t_1589, t_1590, t_1591, pc_x, pc_y, pc_z, nsg_955, nsg_970, \
                         osf0_756, osf1_756, osg_1135, osg_1138, \
                         osg_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1588[k] = f_3 * pc_x[k] * osg_1138[k];

        t_1589[k] = f_3 * pc_x[k] * osg_1139[k];

        t_1590[k] = f_10 * nsg_970[k]
                    + f_1 * osf0_756[k]
                    - f_2 * osf1_756[k]
                    + f_3 * pc_y[k] * osg_1135[k];

        t_1591[k] = f_15 * nsg_955[k]
                    + f_3 * pc_z[k] * osg_1135[k];
    }

#pragma omp simd aligned(t_1592, t_1593, t_1594, pc_y, nsg_972, nsg_973, nsg_974, osf0_758, \
                         osf0_759, osf1_758, osf1_759, osg_1137, osg_1138, \
                         osg_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1592[k] = f_10 * nsg_972[k]
                    + f_6 * osf0_758[k]
                    - f_7 * osf1_758[k]
                    + f_3 * pc_y[k] * osg_1137[k];

        t_1593[k] = f_10 * nsg_973[k]
                    + f_4 * osf0_759[k]
                    - f_5 * osf1_759[k]
                    + f_3 * pc_y[k] * osg_1138[k];

        t_1594[k] = f_10 * nsg_974[k]
                    + f_3 * pc_y[k] * osg_1139[k];
    }

#pragma omp simd aligned(t_1595, t_1596, t_1597, pa_y, pc_x, pc_y, pc_z, nsh0_1365, nsg_959, \
                         nsh1_1365, osf0_759, osf0_761, osf1_759, osf1_761, osg_1139, \
                         osg_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1595[k] = f_15 * nsg_959[k]
                    + f_1 * osf0_759[k]
                    - f_2 * osf1_759[k]
                    + f_3 * pc_z[k] * osg_1139[k];

        t_1596[k] = pa_y[k] * nsh0_1365[k]
                    - f_8 * pc_y[k] * nsh1_1365[k];

        t_1597[k] = f_13 * osf0_761[k]
                    - f_14 * osf1_761[k]
                    + f_3 * pc_x[k] * osg_1141[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, pa_y, pc_x, pc_y, nsh0_1367, nsh1_1367, \
                         osf0_763, osf0_764, osf1_763, osf1_764, osg_1143, \
                         osg_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = pa_y[k] * nsh0_1367[k]
                    - f_8 * pc_y[k] * nsh1_1367[k];

        t_1599[k] = f_6 * osf0_763[k]
                    - f_7 * osf1_763[k]
                    + f_3 * pc_x[k] * osg_1143[k];

        t_1600[k] = f_6 * osf0_764[k]
                    - f_7 * osf1_764[k]
                    + f_3 * pc_x[k] * osg_1144[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, pa_y, pc_x, pc_y, nsh0_1370, nsh1_1370, \
                         osf0_766, osf0_767, osf1_766, osf1_767, osg_1146, \
                         osg_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = pa_y[k] * nsh0_1370[k]
                    - f_8 * pc_y[k] * nsh1_1370[k];

        t_1602[k] = f_4 * osf0_766[k]
                    - f_5 * osf1_766[k]
                    + f_3 * pc_x[k] * osg_1146[k];

        t_1603[k] = f_4 * osf0_767[k]
                    - f_5 * osf1_767[k]
                    + f_3 * pc_x[k] * osg_1147[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, t_1607, t_1608, pa_y, pc_x, pc_y, nsh0_1374, \
                         nsh1_1374, osf0_768, osf1_768, osg_1148, osg_1150, osg_1151, \
                         osg_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = f_4 * osf0_768[k]
                    - f_5 * osf1_768[k]
                    + f_3 * pc_x[k] * osg_1148[k];

        t_1605[k] = pa_y[k] * nsh0_1374[k]
                    - f_8 * pc_y[k] * nsh1_1374[k];

        t_1606[k] = f_3 * pc_x[k] * osg_1150[k];

        t_1607[k] = f_3 * pc_x[k] * osg_1151[k];

        t_1608[k] = f_3 * pc_x[k] * osg_1152[k];
    }

#pragma omp simd aligned(t_1609, t_1610, t_1611, t_1612, pa_y, pc_x, pc_y, pc_z, nsh0_1380, \
                         nsg_970, nsg_985, nsh1_1380, osg_1150, osg_1153, \
                         osg_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1609[k] = f_3 * pc_x[k] * osg_1153[k];

        t_1610[k] = f_3 * pc_x[k] * osg_1154[k];

        t_1611[k] = pa_y[k] * nsh0_1380[k]
                    + f_20 * nsg_985[k]
                    - f_8 * pc_y[k] * nsh1_1380[k];

        t_1612[k] = f_12 * nsg_970[k]
                    + f_3 * pc_z[k] * osg_1150[k];
    }

#pragma omp simd aligned(t_1613, t_1614, t_1615, t_1616, pa_y, pc_y, nsh0_1382, nsh0_1383, \
                         nsh0_1385, nsg_987, nsg_988, nsg_989, nsh1_1382, nsh1_1383, \
                         nsh1_1385, osg_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1613[k] = pa_y[k] * nsh0_1382[k]
                    + f_11 * nsg_987[k]
                    - f_8 * pc_y[k] * nsh1_1382[k];

        t_1614[k] = pa_y[k] * nsh0_1383[k]
                    + f_10 * nsg_988[k]
                    - f_8 * pc_y[k] * nsh1_1383[k];

        t_1615[k] = f_9 * nsg_989[k]
                    + f_3 * pc_y[k] * osg_1154[k];

        t_1616[k] = pa_y[k] * nsh0_1385[k]
                    - f_8 * pc_y[k] * nsh1_1385[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, t_1620, t_1621, pc_x, pc_y, osf0_770, \
                         osf0_772, osf0_773, osf1_770, osf1_772, osf1_773, osg_1155, osg_1157, \
                         osg_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_1 * osf0_770[k]
                    - f_2 * osf1_770[k]
                    + f_3 * pc_x[k] * osg_1155[k];

        t_1618[k] = f_3 * pc_y[k] * osg_1155[k];

        t_1619[k] = f_13 * osf0_772[k]
                    - f_14 * osf1_772[k]
                    + f_3 * pc_x[k] * osg_1157[k];

        t_1620[k] = f_6 * osf0_773[k]
                    - f_7 * osf1_773[k]
                    + f_3 * pc_x[k] * osg_1158[k];

        t_1621[k] = f_3 * pc_y[k] * osg_1157[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, t_1625, pc_x, pc_y, osf0_775, osf0_776, \
                         osf0_777, osf1_775, osf1_776, osf1_777, osg_1160, osg_1161, \
                         osg_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_6 * osf0_775[k]
                    - f_7 * osf1_775[k]
                    + f_3 * pc_x[k] * osg_1160[k];

        t_1623[k] = f_4 * osf0_776[k]
                    - f_5 * osf1_776[k]
                    + f_3 * pc_x[k] * osg_1161[k];

        t_1624[k] = f_4 * osf0_777[k]
                    - f_5 * osf1_777[k]
                    + f_3 * pc_x[k] * osg_1162[k];

        t_1625[k] = f_3 * pc_y[k] * osg_1160[k];
    }

#pragma omp simd aligned(t_1626, t_1627, t_1628, t_1629, t_1630, t_1631, pc_x, osf0_779, \
                         osf1_779, osg_1164, osg_1165, osg_1166, osg_1167, osg_1168, \
                         osg_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1626[k] = f_4 * osf0_779[k]
                    - f_5 * osf1_779[k]
                    + f_3 * pc_x[k] * osg_1164[k];

        t_1627[k] = f_3 * pc_x[k] * osg_1165[k];

        t_1628[k] = f_3 * pc_x[k] * osg_1166[k];

        t_1629[k] = f_3 * pc_x[k] * osg_1167[k];

        t_1630[k] = f_3 * pc_x[k] * osg_1168[k];

        t_1631[k] = f_3 * pc_x[k] * osg_1169[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pc_y, osf0_776, osf0_777, osf0_778, osf1_776, \
                         osf1_777, osf1_778, osg_1165, osg_1166, \
                         osg_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_1 * osf0_776[k]
                    - f_2 * osf1_776[k]
                    + f_3 * pc_y[k] * osg_1165[k];

        t_1633[k] = f_13 * osf0_777[k]
                    - f_14 * osf1_777[k]
                    + f_3 * pc_y[k] * osg_1166[k];

        t_1634[k] = f_6 * osf0_778[k]
                    - f_7 * osf1_778[k]
                    + f_3 * pc_y[k] * osg_1167[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, pc_y, pc_z, nsg_989, osf0_779, osf1_779, \
                         osg_1168, osg_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = f_4 * osf0_779[k]
                    - f_5 * osf1_779[k]
                    + f_3 * pc_y[k] * osg_1168[k];

        t_1636[k] = f_3 * pc_y[k] * osg_1169[k];

        t_1637[k] = f_0 * nsg_989[k]
                    + f_1 * osf0_779[k]
                    - f_2 * osf1_779[k]
                    + f_3 * pc_z[k] * osg_1169[k];
    }
}

auto
compute_prim_osh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t nsh0, const size_t nsg,
                                                   const size_t nsh1, const size_t osf0,
                                                   const size_t osf1, const size_t osg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_osh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, nsh0, nsg,
                                                              nsh1, osf0, osf1, osg, ncols,
                                                              gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, nsh0,
                                                               nsg, nsh1, osg, ncols, gamma, p,
                                                               q);

    compute_prim_osh_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, nsh0,
                                                               nsg, nsh1, osf0, osf1, osg,
                                                               ncols, gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece12(buffer, target, pa, pc, nsh0,
                                                               nsg, nsh1, osf0, osf1, osg,
                                                               ncols, gamma, p, q);

    compute_prim_osh_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, nsh0,
                                                               nsg, nsh1, osf0, osf1, osg,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
