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


#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
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
    const auto f_12 = 5.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;

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

    const auto *osh0_0 = buffer.data(osh0 + 0);
    const auto *osh0_3 = buffer.data(osh0 + 3);
    const auto *osh0_5 = buffer.data(osh0 + 5);
    const auto *osh0_6 = buffer.data(osh0 + 6);
    const auto *osh0_9 = buffer.data(osh0 + 9);
    const auto *osh0_15 = buffer.data(osh0 + 15);
    const auto *osh0_20 = buffer.data(osh0 + 20);
    const auto *osh0_24 = buffer.data(osh0 + 24);
    const auto *osh0_27 = buffer.data(osh0 + 27);
    const auto *osh0_36 = buffer.data(osh0 + 36);
    const auto *osh0_42 = buffer.data(osh0 + 42);
    const auto *osh0_47 = buffer.data(osh0 + 47);
    const auto *osh0_51 = buffer.data(osh0 + 51);
    const auto *osh0_62 = buffer.data(osh0 + 62);

    const auto *osg_0 = buffer.data(osg + 0);
    const auto *osg_1 = buffer.data(osg + 1);
    const auto *osg_2 = buffer.data(osg + 2);
    const auto *osg_3 = buffer.data(osg + 3);
    const auto *osg_5 = buffer.data(osg + 5);
    const auto *osg_10 = buffer.data(osg + 10);
    const auto *osg_12 = buffer.data(osg + 12);
    const auto *osg_14 = buffer.data(osg + 14);
    const auto *osg_15 = buffer.data(osg + 15);
    const auto *osg_18 = buffer.data(osg + 18);
    const auto *osg_20 = buffer.data(osg + 20);
    const auto *osg_25 = buffer.data(osg + 25);
    const auto *osg_27 = buffer.data(osg + 27);
    const auto *osg_28 = buffer.data(osg + 28);
    const auto *osg_29 = buffer.data(osg + 29);
    const auto *osg_30 = buffer.data(osg + 30);
    const auto *osg_32 = buffer.data(osg + 32);
    const auto *osg_35 = buffer.data(osg + 35);
    const auto *osg_40 = buffer.data(osg + 40);
    const auto *osg_41 = buffer.data(osg + 41);
    const auto *osg_42 = buffer.data(osg + 42);
    const auto *osg_43 = buffer.data(osg + 43);
    const auto *osg_44 = buffer.data(osg + 44);
    const auto *osg_45 = buffer.data(osg + 45);
    const auto *osg_48 = buffer.data(osg + 48);
    const auto *osg_51 = buffer.data(osg + 51);
    const auto *osg_55 = buffer.data(osg + 55);
    const auto *osg_57 = buffer.data(osg + 57);
    const auto *osg_58 = buffer.data(osg + 58);
    const auto *osg_59 = buffer.data(osg + 59);
    const auto *osg_70 = buffer.data(osg + 70);
    const auto *osg_71 = buffer.data(osg + 71);
    const auto *osg_72 = buffer.data(osg + 72);
    const auto *osg_73 = buffer.data(osg + 73);
    const auto *osg_74 = buffer.data(osg + 74);
    const auto *osg_75 = buffer.data(osg + 75);
    const auto *osg_80 = buffer.data(osg + 80);
    const auto *osg_84 = buffer.data(osg + 84);
    const auto *osg_85 = buffer.data(osg + 85);
    const auto *osg_86 = buffer.data(osg + 86);
    const auto *osg_87 = buffer.data(osg + 87);
    const auto *osg_89 = buffer.data(osg + 89);
    const auto *osg_90 = buffer.data(osg + 90);
    const auto *osg_93 = buffer.data(osg + 93);

    const auto *osh1_0 = buffer.data(osh1 + 0);
    const auto *osh1_3 = buffer.data(osh1 + 3);
    const auto *osh1_5 = buffer.data(osh1 + 5);
    const auto *osh1_6 = buffer.data(osh1 + 6);
    const auto *osh1_9 = buffer.data(osh1 + 9);
    const auto *osh1_15 = buffer.data(osh1 + 15);
    const auto *osh1_20 = buffer.data(osh1 + 20);
    const auto *osh1_24 = buffer.data(osh1 + 24);
    const auto *osh1_27 = buffer.data(osh1 + 27);
    const auto *osh1_36 = buffer.data(osh1 + 36);
    const auto *osh1_42 = buffer.data(osh1 + 42);
    const auto *osh1_47 = buffer.data(osh1 + 47);
    const auto *osh1_51 = buffer.data(osh1 + 51);
    const auto *osh1_62 = buffer.data(osh1 + 62);

    const auto *qsf0_0 = buffer.data(qsf0 + 0);
    const auto *qsf0_1 = buffer.data(qsf0 + 1);
    const auto *qsf0_2 = buffer.data(qsf0 + 2);
    const auto *qsf0_6 = buffer.data(qsf0 + 6);
    const auto *qsf0_8 = buffer.data(qsf0 + 8);
    const auto *qsf0_9 = buffer.data(qsf0 + 9);
    const auto *qsf0_16 = buffer.data(qsf0 + 16);
    const auto *qsf0_17 = buffer.data(qsf0 + 17);
    const auto *qsf0_22 = buffer.data(qsf0 + 22);
    const auto *qsf0_27 = buffer.data(qsf0 + 27);
    const auto *qsf0_28 = buffer.data(qsf0 + 28);
    const auto *qsf0_29 = buffer.data(qsf0 + 29);
    const auto *qsf0_30 = buffer.data(qsf0 + 30);
    const auto *qsf0_32 = buffer.data(qsf0 + 32);
    const auto *qsf0_33 = buffer.data(qsf0 + 33);
    const auto *qsf0_36 = buffer.data(qsf0 + 36);
    const auto *qsf0_37 = buffer.data(qsf0 + 37);
    const auto *qsf0_39 = buffer.data(qsf0 + 39);
    const auto *qsf0_48 = buffer.data(qsf0 + 48);
    const auto *qsf0_49 = buffer.data(qsf0 + 49);
    const auto *qsf0_50 = buffer.data(qsf0 + 50);
    const auto *qsf0_51 = buffer.data(qsf0 + 51);
    const auto *qsf0_52 = buffer.data(qsf0 + 52);
    const auto *qsf0_55 = buffer.data(qsf0 + 55);
    const auto *qsf0_56 = buffer.data(qsf0 + 56);
    const auto *qsf0_57 = buffer.data(qsf0 + 57);
    const auto *qsf0_58 = buffer.data(qsf0 + 58);
    const auto *qsf0_59 = buffer.data(qsf0 + 59);
    const auto *qsf0_60 = buffer.data(qsf0 + 60);
    const auto *qsf0_63 = buffer.data(qsf0 + 63);

    const auto *qsf1_0 = buffer.data(qsf1 + 0);
    const auto *qsf1_1 = buffer.data(qsf1 + 1);
    const auto *qsf1_2 = buffer.data(qsf1 + 2);
    const auto *qsf1_6 = buffer.data(qsf1 + 6);
    const auto *qsf1_8 = buffer.data(qsf1 + 8);
    const auto *qsf1_9 = buffer.data(qsf1 + 9);
    const auto *qsf1_16 = buffer.data(qsf1 + 16);
    const auto *qsf1_17 = buffer.data(qsf1 + 17);
    const auto *qsf1_22 = buffer.data(qsf1 + 22);
    const auto *qsf1_27 = buffer.data(qsf1 + 27);
    const auto *qsf1_28 = buffer.data(qsf1 + 28);
    const auto *qsf1_29 = buffer.data(qsf1 + 29);
    const auto *qsf1_30 = buffer.data(qsf1 + 30);
    const auto *qsf1_32 = buffer.data(qsf1 + 32);
    const auto *qsf1_33 = buffer.data(qsf1 + 33);
    const auto *qsf1_36 = buffer.data(qsf1 + 36);
    const auto *qsf1_37 = buffer.data(qsf1 + 37);
    const auto *qsf1_39 = buffer.data(qsf1 + 39);
    const auto *qsf1_48 = buffer.data(qsf1 + 48);
    const auto *qsf1_49 = buffer.data(qsf1 + 49);
    const auto *qsf1_50 = buffer.data(qsf1 + 50);
    const auto *qsf1_51 = buffer.data(qsf1 + 51);
    const auto *qsf1_52 = buffer.data(qsf1 + 52);
    const auto *qsf1_55 = buffer.data(qsf1 + 55);
    const auto *qsf1_56 = buffer.data(qsf1 + 56);
    const auto *qsf1_57 = buffer.data(qsf1 + 57);
    const auto *qsf1_58 = buffer.data(qsf1 + 58);
    const auto *qsf1_59 = buffer.data(qsf1 + 59);
    const auto *qsf1_60 = buffer.data(qsf1 + 60);
    const auto *qsf1_63 = buffer.data(qsf1 + 63);

    const auto *qsg_0 = buffer.data(qsg + 0);
    const auto *qsg_1 = buffer.data(qsg + 1);
    const auto *qsg_2 = buffer.data(qsg + 2);
    const auto *qsg_3 = buffer.data(qsg + 3);
    const auto *qsg_5 = buffer.data(qsg + 5);
    const auto *qsg_6 = buffer.data(qsg + 6);
    const auto *qsg_9 = buffer.data(qsg + 9);
    const auto *qsg_10 = buffer.data(qsg + 10);
    const auto *qsg_12 = buffer.data(qsg + 12);
    const auto *qsg_13 = buffer.data(qsg + 13);
    const auto *qsg_14 = buffer.data(qsg + 14);
    const auto *qsg_15 = buffer.data(qsg + 15);
    const auto *qsg_16 = buffer.data(qsg + 16);
    const auto *qsg_18 = buffer.data(qsg + 18);
    const auto *qsg_20 = buffer.data(qsg + 20);
    const auto *qsg_21 = buffer.data(qsg + 21);
    const auto *qsg_25 = buffer.data(qsg + 25);
    const auto *qsg_26 = buffer.data(qsg + 26);
    const auto *qsg_27 = buffer.data(qsg + 27);
    const auto *qsg_28 = buffer.data(qsg + 28);
    const auto *qsg_29 = buffer.data(qsg + 29);
    const auto *qsg_30 = buffer.data(qsg + 30);
    const auto *qsg_32 = buffer.data(qsg + 32);
    const auto *qsg_34 = buffer.data(qsg + 34);
    const auto *qsg_35 = buffer.data(qsg + 35);
    const auto *qsg_39 = buffer.data(qsg + 39);
    const auto *qsg_40 = buffer.data(qsg + 40);
    const auto *qsg_41 = buffer.data(qsg + 41);
    const auto *qsg_42 = buffer.data(qsg + 42);
    const auto *qsg_43 = buffer.data(qsg + 43);
    const auto *qsg_44 = buffer.data(qsg + 44);
    const auto *qsg_45 = buffer.data(qsg + 45);
    const auto *qsg_46 = buffer.data(qsg + 46);
    const auto *qsg_47 = buffer.data(qsg + 47);
    const auto *qsg_48 = buffer.data(qsg + 48);
    const auto *qsg_50 = buffer.data(qsg + 50);
    const auto *qsg_51 = buffer.data(qsg + 51);
    const auto *qsg_55 = buffer.data(qsg + 55);
    const auto *qsg_56 = buffer.data(qsg + 56);
    const auto *qsg_57 = buffer.data(qsg + 57);
    const auto *qsg_58 = buffer.data(qsg + 58);
    const auto *qsg_59 = buffer.data(qsg + 59);
    const auto *qsg_60 = buffer.data(qsg + 60);
    const auto *qsg_62 = buffer.data(qsg + 62);
    const auto *qsg_63 = buffer.data(qsg + 63);
    const auto *qsg_65 = buffer.data(qsg + 65);
    const auto *qsg_70 = buffer.data(qsg + 70);
    const auto *qsg_71 = buffer.data(qsg + 71);
    const auto *qsg_72 = buffer.data(qsg + 72);
    const auto *qsg_73 = buffer.data(qsg + 73);
    const auto *qsg_74 = buffer.data(qsg + 74);
    const auto *qsg_75 = buffer.data(qsg + 75);
    const auto *qsg_76 = buffer.data(qsg + 76);
    const auto *qsg_77 = buffer.data(qsg + 77);
    const auto *qsg_78 = buffer.data(qsg + 78);
    const auto *qsg_79 = buffer.data(qsg + 79);
    const auto *qsg_80 = buffer.data(qsg + 80);
    const auto *qsg_84 = buffer.data(qsg + 84);
    const auto *qsg_85 = buffer.data(qsg + 85);
    const auto *qsg_86 = buffer.data(qsg + 86);
    const auto *qsg_87 = buffer.data(qsg + 87);
    const auto *qsg_88 = buffer.data(qsg + 88);
    const auto *qsg_89 = buffer.data(qsg + 89);
    const auto *qsg_90 = buffer.data(qsg + 90);
    const auto *qsg_93 = buffer.data(qsg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, osg_0, qsf0_0, \
                         qsf1_0, qsg_0, qsg_1, qsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * osg_0[k]
                 + f_1 * qsf0_0[k]
                 - f_2 * qsf1_0[k]
                 + f_3 * pc_x[k] * qsg_0[k];

        t_1[k] = f_3 * pc_y[k] * qsg_0[k];

        t_2[k] = f_3 * pc_z[k] * qsg_0[k];

        t_3[k] = f_4 * qsf0_0[k]
                 - f_5 * qsf1_0[k]
                 + f_3 * pc_y[k] * qsg_1[k];

        t_4[k] = f_3 * pc_y[k] * qsg_2[k];

        t_5[k] = f_4 * qsf0_0[k]
                 - f_5 * qsf1_0[k]
                 + f_3 * pc_z[k] * qsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, osg_10, qsf0_1, qsf0_2, \
                         qsf1_1, qsf1_2, qsg_3, qsg_5, qsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * qsf0_1[k]
                 - f_7 * qsf1_1[k]
                 + f_3 * pc_y[k] * qsg_3[k];

        t_7[k] = f_3 * pc_z[k] * qsg_3[k];

        t_8[k] = f_3 * pc_y[k] * qsg_5[k];

        t_9[k] = f_6 * qsf0_2[k]
                 - f_7 * qsf1_2[k]
                 + f_3 * pc_z[k] * qsg_5[k];

        t_10[k] = f_0 * osg_10[k]
                  + f_3 * pc_x[k] * qsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, osg_12, osg_14, qsg_6, \
                         qsg_9, qsg_12, qsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * qsg_6[k];

        t_12[k] = f_0 * osg_12[k]
                  + f_3 * pc_x[k] * qsg_12[k];

        t_13[k] = f_3 * pc_y[k] * qsg_9[k];

        t_14[k] = f_0 * osg_14[k]
                  + f_3 * pc_x[k] * qsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, qsf0_6, qsf0_8, qsf0_9, qsf1_6, \
                         qsf1_8, qsf1_9, qsg_10, qsg_12, qsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * qsf0_6[k]
                  - f_2 * qsf1_6[k]
                  + f_3 * pc_y[k] * qsg_10[k];

        t_16[k] = f_3 * pc_z[k] * qsg_10[k];

        t_17[k] = f_6 * qsf0_8[k]
                  - f_7 * qsf1_8[k]
                  + f_3 * pc_y[k] * qsg_12[k];

        t_18[k] = f_4 * qsf0_9[k]
                  - f_5 * qsf1_9[k]
                  + f_3 * pc_y[k] * qsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, osh0_0, osg_0, \
                         osh1_0, qsf0_9, qsf1_9, qsg_14, qsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * qsg_14[k];

        t_20[k] = f_1 * qsf0_9[k]
                  - f_2 * qsf1_9[k]
                  + f_3 * pc_z[k] * qsg_14[k];

        t_21[k] = pa_y[k] * osh0_0[k]
                  - f_8 * pc_y[k] * osh1_0[k];

        t_22[k] = f_9 * osg_0[k]
                  + f_3 * pc_y[k] * qsg_15[k];

        t_23[k] = f_3 * pc_z[k] * qsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, osh0_3, osh0_5, osh0_6, \
                         osg_1, osg_3, osh1_3, osh1_5, osh1_6, qsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * osh0_3[k]
                  + f_10 * osg_1[k]
                  - f_8 * pc_y[k] * osh1_3[k];

        t_25[k] = f_3 * pc_z[k] * qsg_16[k];

        t_26[k] = pa_y[k] * osh0_5[k]
                  - f_8 * pc_y[k] * osh1_5[k];

        t_27[k] = pa_y[k] * osh0_6[k]
                  + f_11 * osg_3[k]
                  - f_8 * pc_y[k] * osh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, osh0_9, osg_5, \
                         osg_25, osh1_9, qsg_18, qsg_20, qsg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * qsg_18[k];

        t_29[k] = f_9 * osg_5[k]
                  + f_3 * pc_y[k] * qsg_20[k];

        t_30[k] = pa_y[k] * osh0_9[k]
                  - f_8 * pc_y[k] * osh1_9[k];

        t_31[k] = f_12 * osg_25[k]
                  + f_3 * pc_x[k] * qsg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, osg_27, osg_28, osg_29, qsg_21, \
                         qsg_27, qsg_28, qsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * qsg_21[k];

        t_33[k] = f_12 * osg_27[k]
                  + f_3 * pc_x[k] * qsg_27[k];

        t_34[k] = f_12 * osg_28[k]
                  + f_3 * pc_x[k] * qsg_28[k];

        t_35[k] = f_12 * osg_29[k]
                  + f_3 * pc_x[k] * qsg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, osg_10, qsf0_16, qsf0_17, \
                         qsf1_16, qsf1_17, qsg_25, qsg_26, qsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * osg_10[k]
                  + f_1 * qsf0_16[k]
                  - f_2 * qsf1_16[k]
                  + f_3 * pc_y[k] * qsg_25[k];

        t_37[k] = f_3 * pc_z[k] * qsg_25[k];

        t_38[k] = f_4 * qsf0_16[k]
                  - f_5 * qsf1_16[k]
                  + f_3 * pc_z[k] * qsg_26[k];

        t_39[k] = f_6 * qsf0_17[k]
                  - f_7 * qsf1_17[k]
                  + f_3 * pc_z[k] * qsg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, osh0_0, osh0_20, \
                         osg_14, osh1_0, osh1_20, qsg_29, qsg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * osg_14[k]
                  + f_3 * pc_y[k] * qsg_29[k];

        t_41[k] = pa_y[k] * osh0_20[k]
                  - f_8 * pc_y[k] * osh1_20[k];

        t_42[k] = pa_z[k] * osh0_0[k]
                  - f_8 * pc_z[k] * osh1_0[k];

        t_43[k] = f_3 * pc_y[k] * qsg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, osh0_3, osh0_5, osg_0, \
                         osg_2, osh1_3, osh1_5, qsg_30, qsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * osg_0[k]
                  + f_3 * pc_z[k] * qsg_30[k];

        t_45[k] = pa_z[k] * osh0_3[k]
                  - f_8 * pc_z[k] * osh1_3[k];

        t_46[k] = f_3 * pc_y[k] * qsg_32[k];

        t_47[k] = pa_z[k] * osh0_5[k]
                  + f_10 * osg_2[k]
                  - f_8 * pc_z[k] * osh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, osh0_6, osh0_9, osg_5, \
                         osh1_6, osh1_9, qsf0_22, qsf1_22, qsg_34, \
                         qsg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * osh0_6[k]
                  - f_8 * pc_z[k] * osh1_6[k];

        t_49[k] = f_4 * qsf0_22[k]
                  - f_5 * qsf1_22[k]
                  + f_3 * pc_y[k] * qsg_34[k];

        t_50[k] = f_3 * pc_y[k] * qsg_35[k];

        t_51[k] = pa_z[k] * osh0_9[k]
                  + f_11 * osg_5[k]
                  - f_8 * pc_z[k] * osh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, osg_40, osg_41, osg_42, \
                         osg_44, qsg_39, qsg_40, qsg_41, qsg_42, \
                         qsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * osg_40[k]
                  + f_3 * pc_x[k] * qsg_40[k];

        t_53[k] = f_12 * osg_41[k]
                  + f_3 * pc_x[k] * qsg_41[k];

        t_54[k] = f_12 * osg_42[k]
                  + f_3 * pc_x[k] * qsg_42[k];

        t_55[k] = f_3 * pc_y[k] * qsg_39[k];

        t_56[k] = f_12 * osg_44[k]
                  + f_3 * pc_x[k] * qsg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, osh0_15, osh1_15, qsf0_27, \
                         qsf0_28, qsf1_27, qsf1_28, qsg_41, qsg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * osh0_15[k]
                  - f_8 * pc_z[k] * osh1_15[k];

        t_58[k] = f_13 * qsf0_27[k]
                  - f_14 * qsf1_27[k]
                  + f_3 * pc_y[k] * qsg_41[k];

        t_59[k] = f_6 * qsf0_28[k]
                  - f_7 * qsf1_28[k]
                  + f_3 * pc_y[k] * qsg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, osg_14, osg_45, qsf0_29, \
                         qsf0_30, qsf1_29, qsf1_30, qsg_43, qsg_44, \
                         qsg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * qsf0_29[k]
                  - f_5 * qsf1_29[k]
                  + f_3 * pc_y[k] * qsg_43[k];

        t_61[k] = f_3 * pc_y[k] * qsg_44[k];

        t_62[k] = f_9 * osg_14[k]
                  + f_1 * qsf0_29[k]
                  - f_2 * qsf1_29[k]
                  + f_3 * pc_z[k] * qsg_44[k];

        t_63[k] = f_15 * osg_45[k]
                  + f_1 * qsf0_30[k]
                  - f_2 * qsf1_30[k]
                  + f_3 * pc_x[k] * qsg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, osg_15, osg_48, qsf0_33, \
                         qsf1_33, qsg_45, qsg_46, qsg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * osg_15[k]
                  + f_3 * pc_y[k] * qsg_45[k];

        t_65[k] = f_3 * pc_z[k] * qsg_45[k];

        t_66[k] = f_15 * osg_48[k]
                  + f_6 * qsf0_33[k]
                  - f_7 * qsf1_33[k]
                  + f_3 * pc_x[k] * qsg_48[k];

        t_67[k] = f_3 * pc_z[k] * qsg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, osg_51, qsf0_30, qsf0_36, qsf1_30, \
                         qsf1_36, qsg_47, qsg_48, qsg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * qsf0_30[k]
                  - f_5 * qsf1_30[k]
                  + f_3 * pc_z[k] * qsg_47[k];

        t_69[k] = f_15 * osg_51[k]
                  + f_4 * qsf0_36[k]
                  - f_5 * qsf1_36[k]
                  + f_3 * pc_x[k] * qsg_51[k];

        t_70[k] = f_3 * pc_z[k] * qsg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, osg_20, osg_55, qsf0_32, \
                         qsf1_32, qsg_50, qsg_51, qsg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * osg_20[k]
                  + f_3 * pc_y[k] * qsg_50[k];

        t_72[k] = f_6 * qsf0_32[k]
                  - f_7 * qsf1_32[k]
                  + f_3 * pc_z[k] * qsg_50[k];

        t_73[k] = f_15 * osg_55[k]
                  + f_3 * pc_x[k] * qsg_55[k];

        t_74[k] = f_3 * pc_z[k] * qsg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, osg_25, osg_57, osg_58, osg_59, \
                         qsf0_36, qsf1_36, qsg_55, qsg_57, qsg_58, \
                         qsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * osg_57[k]
                  + f_3 * pc_x[k] * qsg_57[k];

        t_76[k] = f_15 * osg_58[k]
                  + f_3 * pc_x[k] * qsg_58[k];

        t_77[k] = f_15 * osg_59[k]
                  + f_3 * pc_x[k] * qsg_59[k];

        t_78[k] = f_10 * osg_25[k]
                  + f_1 * qsf0_36[k]
                  - f_2 * qsf1_36[k]
                  + f_3 * pc_y[k] * qsg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, osg_29, qsf0_36, qsf0_37, \
                         qsf1_36, qsf1_37, qsg_55, qsg_56, qsg_57, \
                         qsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * qsg_55[k];

        t_80[k] = f_4 * qsf0_36[k]
                  - f_5 * qsf1_36[k]
                  + f_3 * pc_z[k] * qsg_56[k];

        t_81[k] = f_6 * qsf0_37[k]
                  - f_7 * qsf1_37[k]
                  + f_3 * pc_z[k] * qsg_57[k];

        t_82[k] = f_10 * osg_29[k]
                  + f_3 * pc_y[k] * qsg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, osh0_42, osg_15, osg_30, \
                         osh1_42, qsf0_39, qsf1_39, qsg_59, qsg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * qsf0_39[k]
                  - f_2 * qsf1_39[k]
                  + f_3 * pc_z[k] * qsg_59[k];

        t_84[k] = pa_y[k] * osh0_42[k]
                  - f_8 * pc_y[k] * osh1_42[k];

        t_85[k] = f_9 * osg_30[k]
                  + f_3 * pc_y[k] * qsg_60[k];

        t_86[k] = f_9 * osg_15[k]
                  + f_3 * pc_z[k] * qsg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, osh0_24, osh0_27, \
                         osh0_47, osg_32, osh1_24, osh1_27, osh1_47, \
                         qsg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * osh0_24[k]
                  - f_8 * pc_z[k] * osh1_24[k];

        t_88[k] = f_9 * osg_32[k]
                  + f_3 * pc_y[k] * qsg_62[k];

        t_89[k] = pa_y[k] * osh0_47[k]
                  - f_8 * pc_y[k] * osh1_47[k];

        t_90[k] = pa_z[k] * osh0_27[k]
                  - f_8 * pc_z[k] * osh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, osh0_51, osg_18, \
                         osg_35, osg_70, osh1_51, qsg_63, qsg_65, \
                         qsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * osg_18[k]
                  + f_3 * pc_z[k] * qsg_63[k];

        t_92[k] = f_9 * osg_35[k]
                  + f_3 * pc_y[k] * qsg_65[k];

        t_93[k] = pa_y[k] * osh0_51[k]
                  - f_8 * pc_y[k] * osh1_51[k];

        t_94[k] = f_15 * osg_70[k]
                  + f_3 * pc_x[k] * qsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, osg_71, osg_72, osg_73, osg_74, qsg_71, \
                         qsg_72, qsg_73, qsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * osg_71[k]
                  + f_3 * pc_x[k] * qsg_71[k];

        t_96[k] = f_15 * osg_72[k]
                  + f_3 * pc_x[k] * qsg_72[k];

        t_97[k] = f_15 * osg_73[k]
                  + f_3 * pc_x[k] * qsg_73[k];

        t_98[k] = f_15 * osg_74[k]
                  + f_3 * pc_x[k] * qsg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, osh0_36, osg_25, osg_42, \
                         osh1_36, qsf0_48, qsf1_48, qsg_70, qsg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * osh0_36[k]
                  - f_8 * pc_z[k] * osh1_36[k];

        t_100[k] = f_9 * osg_25[k]
                   + f_3 * pc_z[k] * qsg_70[k];

        t_101[k] = f_9 * osg_42[k]
                   + f_6 * qsf0_48[k]
                   - f_7 * qsf1_48[k]
                   + f_3 * pc_y[k] * qsg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, osh0_62, osg_43, osg_44, osh1_62, \
                         qsf0_49, qsf1_49, qsg_73, qsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * osg_43[k]
                   + f_4 * qsf0_49[k]
                   - f_5 * qsf1_49[k]
                   + f_3 * pc_y[k] * qsg_73[k];

        t_103[k] = f_9 * osg_44[k]
                   + f_3 * pc_y[k] * qsg_74[k];

        t_104[k] = pa_y[k] * osh0_62[k]
                   - f_8 * pc_y[k] * osh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, osg_30, osg_75, \
                         qsf0_50, qsf1_50, qsg_75, qsg_76, qsg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * osg_75[k]
                   + f_1 * qsf0_50[k]
                   - f_2 * qsf1_50[k]
                   + f_3 * pc_x[k] * qsg_75[k];

        t_106[k] = f_3 * pc_y[k] * qsg_75[k];

        t_107[k] = f_10 * osg_30[k]
                   + f_3 * pc_z[k] * qsg_75[k];

        t_108[k] = f_4 * qsf0_50[k]
                   - f_5 * qsf1_50[k]
                   + f_3 * pc_y[k] * qsg_76[k];

        t_109[k] = f_3 * pc_y[k] * qsg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, osg_80, qsf0_51, qsf0_52, \
                         qsf0_55, qsf1_51, qsf1_52, qsf1_55, qsg_78, qsg_79, \
                         qsg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * osg_80[k]
                   + f_6 * qsf0_55[k]
                   - f_7 * qsf1_55[k]
                   + f_3 * pc_x[k] * qsg_80[k];

        t_111[k] = f_6 * qsf0_51[k]
                   - f_7 * qsf1_51[k]
                   + f_3 * pc_y[k] * qsg_78[k];

        t_112[k] = f_4 * qsf0_52[k]
                   - f_5 * qsf1_52[k]
                   + f_3 * pc_y[k] * qsg_79[k];

        t_113[k] = f_3 * pc_y[k] * qsg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, osg_84, osg_85, osg_86, osg_87, \
                         qsf0_59, qsf1_59, qsg_84, qsg_85, qsg_86, \
                         qsg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * osg_84[k]
                   + f_4 * qsf0_59[k]
                   - f_5 * qsf1_59[k]
                   + f_3 * pc_x[k] * qsg_84[k];

        t_115[k] = f_15 * osg_85[k]
                   + f_3 * pc_x[k] * qsg_85[k];

        t_116[k] = f_15 * osg_86[k]
                   + f_3 * pc_x[k] * qsg_86[k];

        t_117[k] = f_15 * osg_87[k]
                   + f_3 * pc_x[k] * qsg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, osg_89, qsf0_56, qsf0_57, \
                         qsf1_56, qsf1_57, qsg_84, qsg_85, qsg_86, \
                         qsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * qsg_84[k];

        t_119[k] = f_15 * osg_89[k]
                   + f_3 * pc_x[k] * qsg_89[k];

        t_120[k] = f_1 * qsf0_56[k]
                   - f_2 * qsf1_56[k]
                   + f_3 * pc_y[k] * qsg_85[k];

        t_121[k] = f_13 * qsf0_57[k]
                   - f_14 * qsf1_57[k]
                   + f_3 * pc_y[k] * qsg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, osg_44, qsf0_58, qsf0_59, \
                         qsf1_58, qsf1_59, qsg_87, qsg_88, qsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * qsf0_58[k]
                   - f_7 * qsf1_58[k]
                   + f_3 * pc_y[k] * qsg_87[k];

        t_123[k] = f_4 * qsf0_59[k]
                   - f_5 * qsf1_59[k]
                   + f_3 * pc_y[k] * qsg_88[k];

        t_124[k] = f_3 * pc_y[k] * qsg_89[k];

        t_125[k] = f_10 * osg_44[k]
                   + f_1 * qsf0_59[k]
                   - f_2 * qsf1_59[k]
                   + f_3 * pc_z[k] * qsg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, osg_45, osg_90, osg_93, \
                         qsf0_60, qsf0_63, qsf1_60, qsf1_63, qsg_90, \
                         qsg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_16 * osg_90[k]
                   + f_1 * qsf0_60[k]
                   - f_2 * qsf1_60[k]
                   + f_3 * pc_x[k] * qsg_90[k];

        t_127[k] = f_11 * osg_45[k]
                   + f_3 * pc_y[k] * qsg_90[k];

        t_128[k] = f_3 * pc_z[k] * qsg_90[k];

        t_129[k] = f_16 * osg_93[k]
                   + f_6 * qsf0_63[k]
                   - f_7 * qsf1_63[k]
                   + f_3 * pc_x[k] * qsg_93[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
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

    const auto *osh0_63 = buffer.data(osh0 + 63);
    const auto *osh0_66 = buffer.data(osh0 + 66);
    const auto *osh0_69 = buffer.data(osh0 + 69);
    const auto *osh0_78 = buffer.data(osh0 + 78);
    const auto *osh0_105 = buffer.data(osh0 + 105);
    const auto *osh0_108 = buffer.data(osh0 + 108);
    const auto *osh0_110 = buffer.data(osh0 + 110);
    const auto *osh0_111 = buffer.data(osh0 + 111);
    const auto *osh0_114 = buffer.data(osh0 + 114);
    const auto *osh0_125 = buffer.data(osh0 + 125);
    const auto *osh0_126 = buffer.data(osh0 + 126);
    const auto *osh0_129 = buffer.data(osh0 + 129);
    const auto *osh0_132 = buffer.data(osh0 + 132);
    const auto *osh0_141 = buffer.data(osh0 + 141);

    const auto *osg_45 = buffer.data(osg + 45);
    const auto *osg_48 = buffer.data(osg + 48);
    const auto *osg_50 = buffer.data(osg + 50);
    const auto *osg_55 = buffer.data(osg + 55);
    const auto *osg_59 = buffer.data(osg + 59);
    const auto *osg_60 = buffer.data(osg + 60);
    const auto *osg_62 = buffer.data(osg + 62);
    const auto *osg_63 = buffer.data(osg + 63);
    const auto *osg_65 = buffer.data(osg + 65);
    const auto *osg_70 = buffer.data(osg + 70);
    const auto *osg_72 = buffer.data(osg + 72);
    const auto *osg_73 = buffer.data(osg + 73);
    const auto *osg_74 = buffer.data(osg + 74);
    const auto *osg_75 = buffer.data(osg + 75);
    const auto *osg_76 = buffer.data(osg + 76);
    const auto *osg_77 = buffer.data(osg + 77);
    const auto *osg_78 = buffer.data(osg + 78);
    const auto *osg_80 = buffer.data(osg + 80);
    const auto *osg_85 = buffer.data(osg + 85);
    const auto *osg_87 = buffer.data(osg + 87);
    const auto *osg_88 = buffer.data(osg + 88);
    const auto *osg_89 = buffer.data(osg + 89);
    const auto *osg_90 = buffer.data(osg + 90);
    const auto *osg_93 = buffer.data(osg + 93);
    const auto *osg_95 = buffer.data(osg + 95);
    const auto *osg_96 = buffer.data(osg + 96);
    const auto *osg_100 = buffer.data(osg + 100);
    const auto *osg_102 = buffer.data(osg + 102);
    const auto *osg_103 = buffer.data(osg + 103);
    const auto *osg_104 = buffer.data(osg + 104);
    const auto *osg_105 = buffer.data(osg + 105);
    const auto *osg_107 = buffer.data(osg + 107);
    const auto *osg_110 = buffer.data(osg + 110);
    const auto *osg_114 = buffer.data(osg + 114);
    const auto *osg_115 = buffer.data(osg + 115);
    const auto *osg_116 = buffer.data(osg + 116);
    const auto *osg_117 = buffer.data(osg + 117);
    const auto *osg_118 = buffer.data(osg + 118);
    const auto *osg_119 = buffer.data(osg + 119);
    const auto *osg_130 = buffer.data(osg + 130);
    const auto *osg_131 = buffer.data(osg + 131);
    const auto *osg_132 = buffer.data(osg + 132);
    const auto *osg_133 = buffer.data(osg + 133);
    const auto *osg_134 = buffer.data(osg + 134);
    const auto *osg_135 = buffer.data(osg + 135);
    const auto *osg_140 = buffer.data(osg + 140);
    const auto *osg_144 = buffer.data(osg + 144);
    const auto *osg_145 = buffer.data(osg + 145);
    const auto *osg_146 = buffer.data(osg + 146);
    const auto *osg_147 = buffer.data(osg + 147);
    const auto *osg_149 = buffer.data(osg + 149);
    const auto *osg_150 = buffer.data(osg + 150);
    const auto *osg_153 = buffer.data(osg + 153);
    const auto *osg_156 = buffer.data(osg + 156);
    const auto *osg_160 = buffer.data(osg + 160);
    const auto *osg_162 = buffer.data(osg + 162);
    const auto *osg_163 = buffer.data(osg + 163);
    const auto *osg_164 = buffer.data(osg + 164);
    const auto *osg_170 = buffer.data(osg + 170);
    const auto *osg_174 = buffer.data(osg + 174);
    const auto *osg_175 = buffer.data(osg + 175);
    const auto *osg_176 = buffer.data(osg + 176);
    const auto *osg_177 = buffer.data(osg + 177);
    const auto *osg_178 = buffer.data(osg + 178);
    const auto *osg_179 = buffer.data(osg + 179);

    const auto *osh1_63 = buffer.data(osh1 + 63);
    const auto *osh1_66 = buffer.data(osh1 + 66);
    const auto *osh1_69 = buffer.data(osh1 + 69);
    const auto *osh1_78 = buffer.data(osh1 + 78);
    const auto *osh1_105 = buffer.data(osh1 + 105);
    const auto *osh1_108 = buffer.data(osh1 + 108);
    const auto *osh1_110 = buffer.data(osh1 + 110);
    const auto *osh1_111 = buffer.data(osh1 + 111);
    const auto *osh1_114 = buffer.data(osh1 + 114);
    const auto *osh1_125 = buffer.data(osh1 + 125);
    const auto *osh1_126 = buffer.data(osh1 + 126);
    const auto *osh1_129 = buffer.data(osh1 + 129);
    const auto *osh1_132 = buffer.data(osh1 + 132);
    const auto *osh1_141 = buffer.data(osh1 + 141);

    const auto *qsf0_60 = buffer.data(qsf0 + 60);
    const auto *qsf0_62 = buffer.data(qsf0 + 62);
    const auto *qsf0_66 = buffer.data(qsf0 + 66);
    const auto *qsf0_67 = buffer.data(qsf0 + 67);
    const auto *qsf0_69 = buffer.data(qsf0 + 69);
    const auto *qsf0_75 = buffer.data(qsf0 + 75);
    const auto *qsf0_78 = buffer.data(qsf0 + 78);
    const auto *qsf0_79 = buffer.data(qsf0 + 79);
    const auto *qsf0_86 = buffer.data(qsf0 + 86);
    const auto *qsf0_88 = buffer.data(qsf0 + 88);
    const auto *qsf0_89 = buffer.data(qsf0 + 89);
    const auto *qsf0_90 = buffer.data(qsf0 + 90);
    const auto *qsf0_91 = buffer.data(qsf0 + 91);
    const auto *qsf0_92 = buffer.data(qsf0 + 92);
    const auto *qsf0_95 = buffer.data(qsf0 + 95);
    const auto *qsf0_96 = buffer.data(qsf0 + 96);
    const auto *qsf0_97 = buffer.data(qsf0 + 97);
    const auto *qsf0_98 = buffer.data(qsf0 + 98);
    const auto *qsf0_99 = buffer.data(qsf0 + 99);
    const auto *qsf0_100 = buffer.data(qsf0 + 100);
    const auto *qsf0_102 = buffer.data(qsf0 + 102);
    const auto *qsf0_103 = buffer.data(qsf0 + 103);
    const auto *qsf0_106 = buffer.data(qsf0 + 106);
    const auto *qsf0_107 = buffer.data(qsf0 + 107);
    const auto *qsf0_109 = buffer.data(qsf0 + 109);
    const auto *qsf0_115 = buffer.data(qsf0 + 115);
    const auto *qsf0_118 = buffer.data(qsf0 + 118);
    const auto *qsf0_119 = buffer.data(qsf0 + 119);

    const auto *qsf1_60 = buffer.data(qsf1 + 60);
    const auto *qsf1_62 = buffer.data(qsf1 + 62);
    const auto *qsf1_66 = buffer.data(qsf1 + 66);
    const auto *qsf1_67 = buffer.data(qsf1 + 67);
    const auto *qsf1_69 = buffer.data(qsf1 + 69);
    const auto *qsf1_75 = buffer.data(qsf1 + 75);
    const auto *qsf1_78 = buffer.data(qsf1 + 78);
    const auto *qsf1_79 = buffer.data(qsf1 + 79);
    const auto *qsf1_86 = buffer.data(qsf1 + 86);
    const auto *qsf1_88 = buffer.data(qsf1 + 88);
    const auto *qsf1_89 = buffer.data(qsf1 + 89);
    const auto *qsf1_90 = buffer.data(qsf1 + 90);
    const auto *qsf1_91 = buffer.data(qsf1 + 91);
    const auto *qsf1_92 = buffer.data(qsf1 + 92);
    const auto *qsf1_95 = buffer.data(qsf1 + 95);
    const auto *qsf1_96 = buffer.data(qsf1 + 96);
    const auto *qsf1_97 = buffer.data(qsf1 + 97);
    const auto *qsf1_98 = buffer.data(qsf1 + 98);
    const auto *qsf1_99 = buffer.data(qsf1 + 99);
    const auto *qsf1_100 = buffer.data(qsf1 + 100);
    const auto *qsf1_102 = buffer.data(qsf1 + 102);
    const auto *qsf1_103 = buffer.data(qsf1 + 103);
    const auto *qsf1_106 = buffer.data(qsf1 + 106);
    const auto *qsf1_107 = buffer.data(qsf1 + 107);
    const auto *qsf1_109 = buffer.data(qsf1 + 109);
    const auto *qsf1_115 = buffer.data(qsf1 + 115);
    const auto *qsf1_118 = buffer.data(qsf1 + 118);
    const auto *qsf1_119 = buffer.data(qsf1 + 119);

    const auto *qsg_91 = buffer.data(qsg + 91);
    const auto *qsg_92 = buffer.data(qsg + 92);
    const auto *qsg_93 = buffer.data(qsg + 93);
    const auto *qsg_95 = buffer.data(qsg + 95);
    const auto *qsg_96 = buffer.data(qsg + 96);
    const auto *qsg_100 = buffer.data(qsg + 100);
    const auto *qsg_101 = buffer.data(qsg + 101);
    const auto *qsg_102 = buffer.data(qsg + 102);
    const auto *qsg_103 = buffer.data(qsg + 103);
    const auto *qsg_104 = buffer.data(qsg + 104);
    const auto *qsg_105 = buffer.data(qsg + 105);
    const auto *qsg_107 = buffer.data(qsg + 107);
    const auto *qsg_108 = buffer.data(qsg + 108);
    const auto *qsg_110 = buffer.data(qsg + 110);
    const auto *qsg_114 = buffer.data(qsg + 114);
    const auto *qsg_115 = buffer.data(qsg + 115);
    const auto *qsg_116 = buffer.data(qsg + 116);
    const auto *qsg_117 = buffer.data(qsg + 117);
    const auto *qsg_118 = buffer.data(qsg + 118);
    const auto *qsg_119 = buffer.data(qsg + 119);
    const auto *qsg_120 = buffer.data(qsg + 120);
    const auto *qsg_122 = buffer.data(qsg + 122);
    const auto *qsg_123 = buffer.data(qsg + 123);
    const auto *qsg_125 = buffer.data(qsg + 125);
    const auto *qsg_130 = buffer.data(qsg + 130);
    const auto *qsg_131 = buffer.data(qsg + 131);
    const auto *qsg_132 = buffer.data(qsg + 132);
    const auto *qsg_133 = buffer.data(qsg + 133);
    const auto *qsg_134 = buffer.data(qsg + 134);
    const auto *qsg_135 = buffer.data(qsg + 135);
    const auto *qsg_136 = buffer.data(qsg + 136);
    const auto *qsg_137 = buffer.data(qsg + 137);
    const auto *qsg_138 = buffer.data(qsg + 138);
    const auto *qsg_139 = buffer.data(qsg + 139);
    const auto *qsg_140 = buffer.data(qsg + 140);
    const auto *qsg_144 = buffer.data(qsg + 144);
    const auto *qsg_145 = buffer.data(qsg + 145);
    const auto *qsg_146 = buffer.data(qsg + 146);
    const auto *qsg_147 = buffer.data(qsg + 147);
    const auto *qsg_148 = buffer.data(qsg + 148);
    const auto *qsg_149 = buffer.data(qsg + 149);
    const auto *qsg_150 = buffer.data(qsg + 150);
    const auto *qsg_151 = buffer.data(qsg + 151);
    const auto *qsg_152 = buffer.data(qsg + 152);
    const auto *qsg_153 = buffer.data(qsg + 153);
    const auto *qsg_155 = buffer.data(qsg + 155);
    const auto *qsg_156 = buffer.data(qsg + 156);
    const auto *qsg_160 = buffer.data(qsg + 160);
    const auto *qsg_161 = buffer.data(qsg + 161);
    const auto *qsg_162 = buffer.data(qsg + 162);
    const auto *qsg_163 = buffer.data(qsg + 163);
    const auto *qsg_164 = buffer.data(qsg + 164);
    const auto *qsg_165 = buffer.data(qsg + 165);
    const auto *qsg_167 = buffer.data(qsg + 167);
    const auto *qsg_168 = buffer.data(qsg + 168);
    const auto *qsg_170 = buffer.data(qsg + 170);
    const auto *qsg_174 = buffer.data(qsg + 174);
    const auto *qsg_175 = buffer.data(qsg + 175);
    const auto *qsg_176 = buffer.data(qsg + 176);
    const auto *qsg_177 = buffer.data(qsg + 177);
    const auto *qsg_178 = buffer.data(qsg + 178);
    const auto *qsg_179 = buffer.data(qsg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, osg_96, qsf0_60, qsf0_66, \
                         qsf1_60, qsf1_66, qsg_91, qsg_92, qsg_93, \
                         qsg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * qsg_91[k];

        t_131[k] = f_4 * qsf0_60[k]
                   - f_5 * qsf1_60[k]
                   + f_3 * pc_z[k] * qsg_92[k];

        t_132[k] = f_16 * osg_96[k]
                   + f_4 * qsf0_66[k]
                   - f_5 * qsf1_66[k]
                   + f_3 * pc_x[k] * qsg_96[k];

        t_133[k] = f_3 * pc_z[k] * qsg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, osg_50, osg_100, \
                         qsf0_62, qsf1_62, qsg_95, qsg_96, qsg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * osg_50[k]
                   + f_3 * pc_y[k] * qsg_95[k];

        t_135[k] = f_6 * qsf0_62[k]
                   - f_7 * qsf1_62[k]
                   + f_3 * pc_z[k] * qsg_95[k];

        t_136[k] = f_16 * osg_100[k]
                   + f_3 * pc_x[k] * qsg_100[k];

        t_137[k] = f_3 * pc_z[k] * qsg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, osg_55, osg_102, osg_103, \
                         osg_104, qsf0_66, qsf1_66, qsg_100, qsg_102, qsg_103, \
                         qsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * osg_102[k]
                   + f_3 * pc_x[k] * qsg_102[k];

        t_139[k] = f_16 * osg_103[k]
                   + f_3 * pc_x[k] * qsg_103[k];

        t_140[k] = f_16 * osg_104[k]
                   + f_3 * pc_x[k] * qsg_104[k];

        t_141[k] = f_11 * osg_55[k]
                   + f_1 * qsf0_66[k]
                   - f_2 * qsf1_66[k]
                   + f_3 * pc_y[k] * qsg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, osg_59, qsf0_66, qsf0_67, \
                         qsf1_66, qsf1_67, qsg_100, qsg_101, qsg_102, \
                         qsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * qsg_100[k];

        t_143[k] = f_4 * qsf0_66[k]
                   - f_5 * qsf1_66[k]
                   + f_3 * pc_z[k] * qsg_101[k];

        t_144[k] = f_6 * qsf0_67[k]
                   - f_7 * qsf1_67[k]
                   + f_3 * pc_z[k] * qsg_102[k];

        t_145[k] = f_11 * osg_59[k]
                   + f_3 * pc_y[k] * qsg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, osh0_63, osg_45, \
                         osg_60, osh1_63, qsf0_69, qsf1_69, qsg_104, \
                         qsg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * qsf0_69[k]
                   - f_2 * qsf1_69[k]
                   + f_3 * pc_z[k] * qsg_104[k];

        t_147[k] = pa_z[k] * osh0_63[k]
                   - f_8 * pc_z[k] * osh1_63[k];

        t_148[k] = f_10 * osg_60[k]
                   + f_3 * pc_y[k] * qsg_105[k];

        t_149[k] = f_9 * osg_45[k]
                   + f_3 * pc_z[k] * qsg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, osh0_66, osg_62, \
                         osg_110, osh1_66, qsf0_75, qsf1_75, qsg_107, \
                         qsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * osh0_66[k]
                   - f_8 * pc_z[k] * osh1_66[k];

        t_151[k] = f_10 * osg_62[k]
                   + f_3 * pc_y[k] * qsg_107[k];

        t_152[k] = f_16 * osg_110[k]
                   + f_6 * qsf0_75[k]
                   - f_7 * qsf1_75[k]
                   + f_3 * pc_x[k] * qsg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, osh0_69, osg_48, osg_65, \
                         osh1_69, qsg_108, qsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * osh0_69[k]
                   - f_8 * pc_z[k] * osh1_69[k];

        t_154[k] = f_9 * osg_48[k]
                   + f_3 * pc_z[k] * qsg_108[k];

        t_155[k] = f_10 * osg_65[k]
                   + f_3 * pc_y[k] * qsg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, osg_114, osg_115, osg_116, osg_117, \
                         qsf0_79, qsf1_79, qsg_114, qsg_115, qsg_116, \
                         qsg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * osg_114[k]
                   + f_4 * qsf0_79[k]
                   - f_5 * qsf1_79[k]
                   + f_3 * pc_x[k] * qsg_114[k];

        t_157[k] = f_16 * osg_115[k]
                   + f_3 * pc_x[k] * qsg_115[k];

        t_158[k] = f_16 * osg_116[k]
                   + f_3 * pc_x[k] * qsg_116[k];

        t_159[k] = f_16 * osg_117[k]
                   + f_3 * pc_x[k] * qsg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, osh0_78, osg_55, \
                         osg_118, osg_119, osh1_78, qsg_115, qsg_118, \
                         qsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * osg_118[k]
                   + f_3 * pc_x[k] * qsg_118[k];

        t_161[k] = f_16 * osg_119[k]
                   + f_3 * pc_x[k] * qsg_119[k];

        t_162[k] = pa_z[k] * osh0_78[k]
                   - f_8 * pc_z[k] * osh1_78[k];

        t_163[k] = f_9 * osg_55[k]
                   + f_3 * pc_z[k] * qsg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, osg_72, osg_73, osg_74, qsf0_78, qsf0_79, \
                         qsf1_78, qsf1_79, qsg_117, qsg_118, qsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * osg_72[k]
                   + f_6 * qsf0_78[k]
                   - f_7 * qsf1_78[k]
                   + f_3 * pc_y[k] * qsg_117[k];

        t_165[k] = f_10 * osg_73[k]
                   + f_4 * qsf0_79[k]
                   - f_5 * qsf1_79[k]
                   + f_3 * pc_y[k] * qsg_118[k];

        t_166[k] = f_10 * osg_74[k]
                   + f_3 * pc_y[k] * qsg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, osh0_105, osg_59, \
                         osg_60, osg_75, osh1_105, qsf0_79, qsf1_79, qsg_119, \
                         qsg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * osg_59[k]
                   + f_1 * qsf0_79[k]
                   - f_2 * qsf1_79[k]
                   + f_3 * pc_z[k] * qsg_119[k];

        t_168[k] = pa_y[k] * osh0_105[k]
                   - f_8 * pc_y[k] * osh1_105[k];

        t_169[k] = f_9 * osg_75[k]
                   + f_3 * pc_y[k] * qsg_120[k];

        t_170[k] = f_10 * osg_60[k]
                   + f_3 * pc_z[k] * qsg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, osh0_108, osh0_110, osh0_111, \
                         osg_76, osg_77, osg_78, osh1_108, osh1_110, osh1_111, \
                         qsg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * osh0_108[k]
                   + f_10 * osg_76[k]
                   - f_8 * pc_y[k] * osh1_108[k];

        t_172[k] = f_9 * osg_77[k]
                   + f_3 * pc_y[k] * qsg_122[k];

        t_173[k] = pa_y[k] * osh0_110[k]
                   - f_8 * pc_y[k] * osh1_110[k];

        t_174[k] = pa_y[k] * osh0_111[k]
                   + f_11 * osg_78[k]
                   - f_8 * pc_y[k] * osh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, osh0_114, osg_63, \
                         osg_80, osg_130, osh1_114, qsg_123, qsg_125, \
                         qsg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * osg_63[k]
                   + f_3 * pc_z[k] * qsg_123[k];

        t_176[k] = f_9 * osg_80[k]
                   + f_3 * pc_y[k] * qsg_125[k];

        t_177[k] = pa_y[k] * osh0_114[k]
                   - f_8 * pc_y[k] * osh1_114[k];

        t_178[k] = f_16 * osg_130[k]
                   + f_3 * pc_x[k] * qsg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, osg_131, osg_132, osg_133, osg_134, \
                         qsg_131, qsg_132, qsg_133, qsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_16 * osg_131[k]
                   + f_3 * pc_x[k] * qsg_131[k];

        t_180[k] = f_16 * osg_132[k]
                   + f_3 * pc_x[k] * qsg_132[k];

        t_181[k] = f_16 * osg_133[k]
                   + f_3 * pc_x[k] * qsg_133[k];

        t_182[k] = f_16 * osg_134[k]
                   + f_3 * pc_x[k] * qsg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, osg_70, osg_85, osg_87, qsf0_86, \
                         qsf0_88, qsf1_86, qsf1_88, qsg_130, qsg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * osg_85[k]
                   + f_1 * qsf0_86[k]
                   - f_2 * qsf1_86[k]
                   + f_3 * pc_y[k] * qsg_130[k];

        t_184[k] = f_10 * osg_70[k]
                   + f_3 * pc_z[k] * qsg_130[k];

        t_185[k] = f_9 * osg_87[k]
                   + f_6 * qsf0_88[k]
                   - f_7 * qsf1_88[k]
                   + f_3 * pc_y[k] * qsg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, osh0_125, osg_88, osg_89, osh1_125, \
                         qsf0_89, qsf1_89, qsg_133, qsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * osg_88[k]
                   + f_4 * qsf0_89[k]
                   - f_5 * qsf1_89[k]
                   + f_3 * pc_y[k] * qsg_133[k];

        t_187[k] = f_9 * osg_89[k]
                   + f_3 * pc_y[k] * qsg_134[k];

        t_188[k] = pa_y[k] * osh0_125[k]
                   - f_8 * pc_y[k] * osh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, osg_75, osg_135, \
                         qsf0_90, qsf1_90, qsg_135, qsg_136, qsg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * osg_135[k]
                   + f_1 * qsf0_90[k]
                   - f_2 * qsf1_90[k]
                   + f_3 * pc_x[k] * qsg_135[k];

        t_190[k] = f_3 * pc_y[k] * qsg_135[k];

        t_191[k] = f_11 * osg_75[k]
                   + f_3 * pc_z[k] * qsg_135[k];

        t_192[k] = f_4 * qsf0_90[k]
                   - f_5 * qsf1_90[k]
                   + f_3 * pc_y[k] * qsg_136[k];

        t_193[k] = f_3 * pc_y[k] * qsg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, osg_140, qsf0_91, qsf0_92, \
                         qsf0_95, qsf1_91, qsf1_92, qsf1_95, qsg_138, qsg_139, \
                         qsg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_16 * osg_140[k]
                   + f_6 * qsf0_95[k]
                   - f_7 * qsf1_95[k]
                   + f_3 * pc_x[k] * qsg_140[k];

        t_195[k] = f_6 * qsf0_91[k]
                   - f_7 * qsf1_91[k]
                   + f_3 * pc_y[k] * qsg_138[k];

        t_196[k] = f_4 * qsf0_92[k]
                   - f_5 * qsf1_92[k]
                   + f_3 * pc_y[k] * qsg_139[k];

        t_197[k] = f_3 * pc_y[k] * qsg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, osg_144, osg_145, osg_146, osg_147, \
                         qsf0_99, qsf1_99, qsg_144, qsg_145, qsg_146, \
                         qsg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * osg_144[k]
                   + f_4 * qsf0_99[k]
                   - f_5 * qsf1_99[k]
                   + f_3 * pc_x[k] * qsg_144[k];

        t_199[k] = f_16 * osg_145[k]
                   + f_3 * pc_x[k] * qsg_145[k];

        t_200[k] = f_16 * osg_146[k]
                   + f_3 * pc_x[k] * qsg_146[k];

        t_201[k] = f_16 * osg_147[k]
                   + f_3 * pc_x[k] * qsg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, osg_149, qsf0_96, qsf0_97, \
                         qsf1_96, qsf1_97, qsg_144, qsg_145, qsg_146, \
                         qsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * qsg_144[k];

        t_203[k] = f_16 * osg_149[k]
                   + f_3 * pc_x[k] * qsg_149[k];

        t_204[k] = f_1 * qsf0_96[k]
                   - f_2 * qsf1_96[k]
                   + f_3 * pc_y[k] * qsg_145[k];

        t_205[k] = f_13 * qsf0_97[k]
                   - f_14 * qsf1_97[k]
                   + f_3 * pc_y[k] * qsg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, osg_89, qsf0_98, qsf0_99, \
                         qsf1_98, qsf1_99, qsg_147, qsg_148, qsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * qsf0_98[k]
                   - f_7 * qsf1_98[k]
                   + f_3 * pc_y[k] * qsg_147[k];

        t_207[k] = f_4 * qsf0_99[k]
                   - f_5 * qsf1_99[k]
                   + f_3 * pc_y[k] * qsg_148[k];

        t_208[k] = f_3 * pc_y[k] * qsg_149[k];

        t_209[k] = f_11 * osg_89[k]
                   + f_1 * qsf0_99[k]
                   - f_2 * qsf1_99[k]
                   + f_3 * pc_z[k] * qsg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, osg_90, osg_150, \
                         osg_153, qsf0_100, qsf0_103, qsf1_100, qsf1_103, qsg_150, \
                         qsg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * osg_150[k]
                   + f_1 * qsf0_100[k]
                   - f_2 * qsf1_100[k]
                   + f_3 * pc_x[k] * qsg_150[k];

        t_211[k] = f_18 * osg_90[k]
                   + f_3 * pc_y[k] * qsg_150[k];

        t_212[k] = f_3 * pc_z[k] * qsg_150[k];

        t_213[k] = f_17 * osg_153[k]
                   + f_6 * qsf0_103[k]
                   - f_7 * qsf1_103[k]
                   + f_3 * pc_x[k] * qsg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, osg_156, qsf0_100, qsf0_106, \
                         qsf1_100, qsf1_106, qsg_151, qsg_152, qsg_153, \
                         qsg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * qsg_151[k];

        t_215[k] = f_4 * qsf0_100[k]
                   - f_5 * qsf1_100[k]
                   + f_3 * pc_z[k] * qsg_152[k];

        t_216[k] = f_17 * osg_156[k]
                   + f_4 * qsf0_106[k]
                   - f_5 * qsf1_106[k]
                   + f_3 * pc_x[k] * qsg_156[k];

        t_217[k] = f_3 * pc_z[k] * qsg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, osg_95, osg_160, \
                         qsf0_102, qsf1_102, qsg_155, qsg_156, \
                         qsg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_18 * osg_95[k]
                   + f_3 * pc_y[k] * qsg_155[k];

        t_219[k] = f_6 * qsf0_102[k]
                   - f_7 * qsf1_102[k]
                   + f_3 * pc_z[k] * qsg_155[k];

        t_220[k] = f_17 * osg_160[k]
                   + f_3 * pc_x[k] * qsg_160[k];

        t_221[k] = f_3 * pc_z[k] * qsg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, osg_100, osg_162, osg_163, \
                         osg_164, qsf0_106, qsf1_106, qsg_160, qsg_162, qsg_163, \
                         qsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_17 * osg_162[k]
                   + f_3 * pc_x[k] * qsg_162[k];

        t_223[k] = f_17 * osg_163[k]
                   + f_3 * pc_x[k] * qsg_163[k];

        t_224[k] = f_17 * osg_164[k]
                   + f_3 * pc_x[k] * qsg_164[k];

        t_225[k] = f_18 * osg_100[k]
                   + f_1 * qsf0_106[k]
                   - f_2 * qsf1_106[k]
                   + f_3 * pc_y[k] * qsg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, osg_104, qsf0_106, qsf0_107, \
                         qsf1_106, qsf1_107, qsg_160, qsg_161, qsg_162, \
                         qsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * qsg_160[k];

        t_227[k] = f_4 * qsf0_106[k]
                   - f_5 * qsf1_106[k]
                   + f_3 * pc_z[k] * qsg_161[k];

        t_228[k] = f_6 * qsf0_107[k]
                   - f_7 * qsf1_107[k]
                   + f_3 * pc_z[k] * qsg_162[k];

        t_229[k] = f_18 * osg_104[k]
                   + f_3 * pc_y[k] * qsg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, osh0_126, osg_90, \
                         osg_105, osh1_126, qsf0_109, qsf1_109, qsg_164, \
                         qsg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * qsf0_109[k]
                   - f_2 * qsf1_109[k]
                   + f_3 * pc_z[k] * qsg_164[k];

        t_231[k] = pa_z[k] * osh0_126[k]
                   - f_8 * pc_z[k] * osh1_126[k];

        t_232[k] = f_11 * osg_105[k]
                   + f_3 * pc_y[k] * qsg_165[k];

        t_233[k] = f_9 * osg_90[k]
                   + f_3 * pc_z[k] * qsg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, osh0_129, osg_107, \
                         osg_170, osh1_129, qsf0_115, qsf1_115, qsg_167, \
                         qsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * osh0_129[k]
                   - f_8 * pc_z[k] * osh1_129[k];

        t_235[k] = f_11 * osg_107[k]
                   + f_3 * pc_y[k] * qsg_167[k];

        t_236[k] = f_17 * osg_170[k]
                   + f_6 * qsf0_115[k]
                   - f_7 * qsf1_115[k]
                   + f_3 * pc_x[k] * qsg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, osh0_132, osg_93, osg_110, \
                         osh1_132, qsg_168, qsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * osh0_132[k]
                   - f_8 * pc_z[k] * osh1_132[k];

        t_238[k] = f_9 * osg_93[k]
                   + f_3 * pc_z[k] * qsg_168[k];

        t_239[k] = f_11 * osg_110[k]
                   + f_3 * pc_y[k] * qsg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, osg_174, osg_175, osg_176, osg_177, \
                         qsf0_119, qsf1_119, qsg_174, qsg_175, qsg_176, \
                         qsg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * osg_174[k]
                   + f_4 * qsf0_119[k]
                   - f_5 * qsf1_119[k]
                   + f_3 * pc_x[k] * qsg_174[k];

        t_241[k] = f_17 * osg_175[k]
                   + f_3 * pc_x[k] * qsg_175[k];

        t_242[k] = f_17 * osg_176[k]
                   + f_3 * pc_x[k] * qsg_176[k];

        t_243[k] = f_17 * osg_177[k]
                   + f_3 * pc_x[k] * qsg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, osh0_141, osg_100, \
                         osg_178, osg_179, osh1_141, qsg_175, qsg_178, \
                         qsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_17 * osg_178[k]
                   + f_3 * pc_x[k] * qsg_178[k];

        t_245[k] = f_17 * osg_179[k]
                   + f_3 * pc_x[k] * qsg_179[k];

        t_246[k] = pa_z[k] * osh0_141[k]
                   - f_8 * pc_z[k] * osh1_141[k];

        t_247[k] = f_9 * osg_100[k]
                   + f_3 * pc_z[k] * qsg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, osg_117, osg_118, osg_119, qsf0_118, \
                         qsf0_119, qsf1_118, qsf1_119, qsg_177, qsg_178, \
                         qsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * osg_117[k]
                   + f_6 * qsf0_118[k]
                   - f_7 * qsf1_118[k]
                   + f_3 * pc_y[k] * qsg_177[k];

        t_249[k] = f_11 * osg_118[k]
                   + f_4 * qsf0_119[k]
                   - f_5 * qsf1_119[k]
                   + f_3 * pc_y[k] * qsg_178[k];

        t_250[k] = f_11 * osg_119[k]
                   + f_3 * pc_y[k] * qsg_179[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
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

    const auto *osh0_189 = buffer.data(osh0 + 189);
    const auto *osh0_192 = buffer.data(osh0 + 192);
    const auto *osh0_194 = buffer.data(osh0 + 194);
    const auto *osh0_195 = buffer.data(osh0 + 195);
    const auto *osh0_198 = buffer.data(osh0 + 198);
    const auto *osh0_209 = buffer.data(osh0 + 209);
    const auto *osh0_210 = buffer.data(osh0 + 210);
    const auto *osh0_213 = buffer.data(osh0 + 213);
    const auto *osh0_216 = buffer.data(osh0 + 216);
    const auto *osh0_225 = buffer.data(osh0 + 225);

    const auto *osg_104 = buffer.data(osg + 104);
    const auto *osg_105 = buffer.data(osg + 105);
    const auto *osg_108 = buffer.data(osg + 108);
    const auto *osg_115 = buffer.data(osg + 115);
    const auto *osg_119 = buffer.data(osg + 119);
    const auto *osg_120 = buffer.data(osg + 120);
    const auto *osg_122 = buffer.data(osg + 122);
    const auto *osg_123 = buffer.data(osg + 123);
    const auto *osg_125 = buffer.data(osg + 125);
    const auto *osg_130 = buffer.data(osg + 130);
    const auto *osg_132 = buffer.data(osg + 132);
    const auto *osg_133 = buffer.data(osg + 133);
    const auto *osg_134 = buffer.data(osg + 134);
    const auto *osg_135 = buffer.data(osg + 135);
    const auto *osg_136 = buffer.data(osg + 136);
    const auto *osg_137 = buffer.data(osg + 137);
    const auto *osg_138 = buffer.data(osg + 138);
    const auto *osg_140 = buffer.data(osg + 140);
    const auto *osg_145 = buffer.data(osg + 145);
    const auto *osg_147 = buffer.data(osg + 147);
    const auto *osg_148 = buffer.data(osg + 148);
    const auto *osg_149 = buffer.data(osg + 149);
    const auto *osg_150 = buffer.data(osg + 150);
    const auto *osg_153 = buffer.data(osg + 153);
    const auto *osg_155 = buffer.data(osg + 155);
    const auto *osg_160 = buffer.data(osg + 160);
    const auto *osg_164 = buffer.data(osg + 164);
    const auto *osg_165 = buffer.data(osg + 165);
    const auto *osg_167 = buffer.data(osg + 167);
    const auto *osg_168 = buffer.data(osg + 168);
    const auto *osg_170 = buffer.data(osg + 170);
    const auto *osg_177 = buffer.data(osg + 177);
    const auto *osg_178 = buffer.data(osg + 178);
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
    const auto *osg_205 = buffer.data(osg + 205);
    const auto *osg_206 = buffer.data(osg + 206);
    const auto *osg_207 = buffer.data(osg + 207);
    const auto *osg_208 = buffer.data(osg + 208);
    const auto *osg_209 = buffer.data(osg + 209);
    const auto *osg_210 = buffer.data(osg + 210);
    const auto *osg_215 = buffer.data(osg + 215);
    const auto *osg_219 = buffer.data(osg + 219);
    const auto *osg_220 = buffer.data(osg + 220);
    const auto *osg_221 = buffer.data(osg + 221);
    const auto *osg_222 = buffer.data(osg + 222);
    const auto *osg_224 = buffer.data(osg + 224);
    const auto *osg_225 = buffer.data(osg + 225);
    const auto *osg_228 = buffer.data(osg + 228);
    const auto *osg_231 = buffer.data(osg + 231);
    const auto *osg_235 = buffer.data(osg + 235);
    const auto *osg_237 = buffer.data(osg + 237);
    const auto *osg_238 = buffer.data(osg + 238);
    const auto *osg_239 = buffer.data(osg + 239);
    const auto *osg_245 = buffer.data(osg + 245);
    const auto *osg_249 = buffer.data(osg + 249);
    const auto *osg_250 = buffer.data(osg + 250);
    const auto *osg_251 = buffer.data(osg + 251);
    const auto *osg_252 = buffer.data(osg + 252);
    const auto *osg_253 = buffer.data(osg + 253);
    const auto *osg_254 = buffer.data(osg + 254);
    const auto *osg_255 = buffer.data(osg + 255);
    const auto *osg_258 = buffer.data(osg + 258);
    const auto *osg_260 = buffer.data(osg + 260);
    const auto *osg_261 = buffer.data(osg + 261);
    const auto *osg_264 = buffer.data(osg + 264);
    const auto *osg_265 = buffer.data(osg + 265);
    const auto *osg_266 = buffer.data(osg + 266);

    const auto *osh1_189 = buffer.data(osh1 + 189);
    const auto *osh1_192 = buffer.data(osh1 + 192);
    const auto *osh1_194 = buffer.data(osh1 + 194);
    const auto *osh1_195 = buffer.data(osh1 + 195);
    const auto *osh1_198 = buffer.data(osh1 + 198);
    const auto *osh1_209 = buffer.data(osh1 + 209);
    const auto *osh1_210 = buffer.data(osh1 + 210);
    const auto *osh1_213 = buffer.data(osh1 + 213);
    const auto *osh1_216 = buffer.data(osh1 + 216);
    const auto *osh1_225 = buffer.data(osh1 + 225);

    const auto *qsf0_119 = buffer.data(qsf0 + 119);
    const auto *qsf0_120 = buffer.data(qsf0 + 120);
    const auto *qsf0_123 = buffer.data(qsf0 + 123);
    const auto *qsf0_125 = buffer.data(qsf0 + 125);
    const auto *qsf0_126 = buffer.data(qsf0 + 126);
    const auto *qsf0_128 = buffer.data(qsf0 + 128);
    const auto *qsf0_129 = buffer.data(qsf0 + 129);
    const auto *qsf0_136 = buffer.data(qsf0 + 136);
    const auto *qsf0_138 = buffer.data(qsf0 + 138);
    const auto *qsf0_139 = buffer.data(qsf0 + 139);
    const auto *qsf0_140 = buffer.data(qsf0 + 140);
    const auto *qsf0_141 = buffer.data(qsf0 + 141);
    const auto *qsf0_142 = buffer.data(qsf0 + 142);
    const auto *qsf0_145 = buffer.data(qsf0 + 145);
    const auto *qsf0_146 = buffer.data(qsf0 + 146);
    const auto *qsf0_147 = buffer.data(qsf0 + 147);
    const auto *qsf0_148 = buffer.data(qsf0 + 148);
    const auto *qsf0_149 = buffer.data(qsf0 + 149);
    const auto *qsf0_150 = buffer.data(qsf0 + 150);
    const auto *qsf0_152 = buffer.data(qsf0 + 152);
    const auto *qsf0_153 = buffer.data(qsf0 + 153);
    const auto *qsf0_156 = buffer.data(qsf0 + 156);
    const auto *qsf0_157 = buffer.data(qsf0 + 157);
    const auto *qsf0_159 = buffer.data(qsf0 + 159);
    const auto *qsf0_165 = buffer.data(qsf0 + 165);
    const auto *qsf0_168 = buffer.data(qsf0 + 168);
    const auto *qsf0_169 = buffer.data(qsf0 + 169);
    const auto *qsf0_170 = buffer.data(qsf0 + 170);
    const auto *qsf0_173 = buffer.data(qsf0 + 173);
    const auto *qsf0_175 = buffer.data(qsf0 + 175);
    const auto *qsf0_176 = buffer.data(qsf0 + 176);
    const auto *qsf0_179 = buffer.data(qsf0 + 179);

    const auto *qsf1_119 = buffer.data(qsf1 + 119);
    const auto *qsf1_120 = buffer.data(qsf1 + 120);
    const auto *qsf1_123 = buffer.data(qsf1 + 123);
    const auto *qsf1_125 = buffer.data(qsf1 + 125);
    const auto *qsf1_126 = buffer.data(qsf1 + 126);
    const auto *qsf1_128 = buffer.data(qsf1 + 128);
    const auto *qsf1_129 = buffer.data(qsf1 + 129);
    const auto *qsf1_136 = buffer.data(qsf1 + 136);
    const auto *qsf1_138 = buffer.data(qsf1 + 138);
    const auto *qsf1_139 = buffer.data(qsf1 + 139);
    const auto *qsf1_140 = buffer.data(qsf1 + 140);
    const auto *qsf1_141 = buffer.data(qsf1 + 141);
    const auto *qsf1_142 = buffer.data(qsf1 + 142);
    const auto *qsf1_145 = buffer.data(qsf1 + 145);
    const auto *qsf1_146 = buffer.data(qsf1 + 146);
    const auto *qsf1_147 = buffer.data(qsf1 + 147);
    const auto *qsf1_148 = buffer.data(qsf1 + 148);
    const auto *qsf1_149 = buffer.data(qsf1 + 149);
    const auto *qsf1_150 = buffer.data(qsf1 + 150);
    const auto *qsf1_152 = buffer.data(qsf1 + 152);
    const auto *qsf1_153 = buffer.data(qsf1 + 153);
    const auto *qsf1_156 = buffer.data(qsf1 + 156);
    const auto *qsf1_157 = buffer.data(qsf1 + 157);
    const auto *qsf1_159 = buffer.data(qsf1 + 159);
    const auto *qsf1_165 = buffer.data(qsf1 + 165);
    const auto *qsf1_168 = buffer.data(qsf1 + 168);
    const auto *qsf1_169 = buffer.data(qsf1 + 169);
    const auto *qsf1_170 = buffer.data(qsf1 + 170);
    const auto *qsf1_173 = buffer.data(qsf1 + 173);
    const auto *qsf1_175 = buffer.data(qsf1 + 175);
    const auto *qsf1_176 = buffer.data(qsf1 + 176);
    const auto *qsf1_179 = buffer.data(qsf1 + 179);

    const auto *qsg_179 = buffer.data(qsg + 179);
    const auto *qsg_180 = buffer.data(qsg + 180);
    const auto *qsg_182 = buffer.data(qsg + 182);
    const auto *qsg_183 = buffer.data(qsg + 183);
    const auto *qsg_185 = buffer.data(qsg + 185);
    const auto *qsg_186 = buffer.data(qsg + 186);
    const auto *qsg_189 = buffer.data(qsg + 189);
    const auto *qsg_190 = buffer.data(qsg + 190);
    const auto *qsg_191 = buffer.data(qsg + 191);
    const auto *qsg_192 = buffer.data(qsg + 192);
    const auto *qsg_193 = buffer.data(qsg + 193);
    const auto *qsg_194 = buffer.data(qsg + 194);
    const auto *qsg_195 = buffer.data(qsg + 195);
    const auto *qsg_197 = buffer.data(qsg + 197);
    const auto *qsg_198 = buffer.data(qsg + 198);
    const auto *qsg_200 = buffer.data(qsg + 200);
    const auto *qsg_205 = buffer.data(qsg + 205);
    const auto *qsg_206 = buffer.data(qsg + 206);
    const auto *qsg_207 = buffer.data(qsg + 207);
    const auto *qsg_208 = buffer.data(qsg + 208);
    const auto *qsg_209 = buffer.data(qsg + 209);
    const auto *qsg_210 = buffer.data(qsg + 210);
    const auto *qsg_211 = buffer.data(qsg + 211);
    const auto *qsg_212 = buffer.data(qsg + 212);
    const auto *qsg_213 = buffer.data(qsg + 213);
    const auto *qsg_214 = buffer.data(qsg + 214);
    const auto *qsg_215 = buffer.data(qsg + 215);
    const auto *qsg_219 = buffer.data(qsg + 219);
    const auto *qsg_220 = buffer.data(qsg + 220);
    const auto *qsg_221 = buffer.data(qsg + 221);
    const auto *qsg_222 = buffer.data(qsg + 222);
    const auto *qsg_223 = buffer.data(qsg + 223);
    const auto *qsg_224 = buffer.data(qsg + 224);
    const auto *qsg_225 = buffer.data(qsg + 225);
    const auto *qsg_226 = buffer.data(qsg + 226);
    const auto *qsg_227 = buffer.data(qsg + 227);
    const auto *qsg_228 = buffer.data(qsg + 228);
    const auto *qsg_230 = buffer.data(qsg + 230);
    const auto *qsg_231 = buffer.data(qsg + 231);
    const auto *qsg_235 = buffer.data(qsg + 235);
    const auto *qsg_236 = buffer.data(qsg + 236);
    const auto *qsg_237 = buffer.data(qsg + 237);
    const auto *qsg_238 = buffer.data(qsg + 238);
    const auto *qsg_239 = buffer.data(qsg + 239);
    const auto *qsg_240 = buffer.data(qsg + 240);
    const auto *qsg_242 = buffer.data(qsg + 242);
    const auto *qsg_243 = buffer.data(qsg + 243);
    const auto *qsg_245 = buffer.data(qsg + 245);
    const auto *qsg_249 = buffer.data(qsg + 249);
    const auto *qsg_250 = buffer.data(qsg + 250);
    const auto *qsg_251 = buffer.data(qsg + 251);
    const auto *qsg_252 = buffer.data(qsg + 252);
    const auto *qsg_253 = buffer.data(qsg + 253);
    const auto *qsg_254 = buffer.data(qsg + 254);
    const auto *qsg_255 = buffer.data(qsg + 255);
    const auto *qsg_257 = buffer.data(qsg + 257);
    const auto *qsg_258 = buffer.data(qsg + 258);
    const auto *qsg_260 = buffer.data(qsg + 260);
    const auto *qsg_261 = buffer.data(qsg + 261);
    const auto *qsg_264 = buffer.data(qsg + 264);
    const auto *qsg_265 = buffer.data(qsg + 265);
    const auto *qsg_266 = buffer.data(qsg + 266);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, osg_104, osg_120, osg_180, \
                         qsf0_119, qsf0_120, qsf1_119, qsf1_120, qsg_179, \
                         qsg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * osg_104[k]
                   + f_1 * qsf0_119[k]
                   - f_2 * qsf1_119[k]
                   + f_3 * pc_z[k] * qsg_179[k];

        t_252[k] = f_17 * osg_180[k]
                   + f_1 * qsf0_120[k]
                   - f_2 * qsf1_120[k]
                   + f_3 * pc_x[k] * qsg_180[k];

        t_253[k] = f_10 * osg_120[k]
                   + f_3 * pc_y[k] * qsg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, osg_105, osg_122, osg_183, \
                         qsf0_123, qsf1_123, qsg_180, qsg_182, \
                         qsg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * osg_105[k]
                   + f_3 * pc_z[k] * qsg_180[k];

        t_255[k] = f_17 * osg_183[k]
                   + f_6 * qsf0_123[k]
                   - f_7 * qsf1_123[k]
                   + f_3 * pc_x[k] * qsg_183[k];

        t_256[k] = f_10 * osg_122[k]
                   + f_3 * pc_y[k] * qsg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, osg_108, osg_185, osg_186, qsf0_125, \
                         qsf0_126, qsf1_125, qsf1_126, qsg_183, qsg_185, \
                         qsg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * osg_185[k]
                   + f_6 * qsf0_125[k]
                   - f_7 * qsf1_125[k]
                   + f_3 * pc_x[k] * qsg_185[k];

        t_258[k] = f_17 * osg_186[k]
                   + f_4 * qsf0_126[k]
                   - f_5 * qsf1_126[k]
                   + f_3 * pc_x[k] * qsg_186[k];

        t_259[k] = f_10 * osg_108[k]
                   + f_3 * pc_z[k] * qsg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, osg_125, osg_189, osg_190, \
                         osg_191, qsf0_129, qsf1_129, qsg_185, qsg_189, qsg_190, \
                         qsg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * osg_125[k]
                   + f_3 * pc_y[k] * qsg_185[k];

        t_261[k] = f_17 * osg_189[k]
                   + f_4 * qsf0_129[k]
                   - f_5 * qsf1_129[k]
                   + f_3 * pc_x[k] * qsg_189[k];

        t_262[k] = f_17 * osg_190[k]
                   + f_3 * pc_x[k] * qsg_190[k];

        t_263[k] = f_17 * osg_191[k]
                   + f_3 * pc_x[k] * qsg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, osg_130, osg_192, osg_193, \
                         osg_194, qsf0_126, qsf1_126, qsg_190, qsg_192, qsg_193, \
                         qsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * osg_192[k]
                   + f_3 * pc_x[k] * qsg_192[k];

        t_265[k] = f_17 * osg_193[k]
                   + f_3 * pc_x[k] * qsg_193[k];

        t_266[k] = f_17 * osg_194[k]
                   + f_3 * pc_x[k] * qsg_194[k];

        t_267[k] = f_10 * osg_130[k]
                   + f_1 * qsf0_126[k]
                   - f_2 * qsf1_126[k]
                   + f_3 * pc_y[k] * qsg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, osg_115, osg_132, osg_133, qsf0_128, \
                         qsf0_129, qsf1_128, qsf1_129, qsg_190, qsg_192, \
                         qsg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * osg_115[k]
                   + f_3 * pc_z[k] * qsg_190[k];

        t_269[k] = f_10 * osg_132[k]
                   + f_6 * qsf0_128[k]
                   - f_7 * qsf1_128[k]
                   + f_3 * pc_y[k] * qsg_192[k];

        t_270[k] = f_10 * osg_133[k]
                   + f_4 * qsf0_129[k]
                   - f_5 * qsf1_129[k]
                   + f_3 * pc_y[k] * qsg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, osh0_189, osg_119, \
                         osg_134, osg_135, osh1_189, qsf0_129, qsf1_129, qsg_194, \
                         qsg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * osg_134[k]
                   + f_3 * pc_y[k] * qsg_194[k];

        t_272[k] = f_10 * osg_119[k]
                   + f_1 * qsf0_129[k]
                   - f_2 * qsf1_129[k]
                   + f_3 * pc_z[k] * qsg_194[k];

        t_273[k] = pa_y[k] * osh0_189[k]
                   - f_8 * pc_y[k] * osh1_189[k];

        t_274[k] = f_9 * osg_135[k]
                   + f_3 * pc_y[k] * qsg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, osh0_192, osh0_194, \
                         osg_120, osg_136, osg_137, osh1_192, osh1_194, qsg_195, \
                         qsg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * osg_120[k]
                   + f_3 * pc_z[k] * qsg_195[k];

        t_276[k] = pa_y[k] * osh0_192[k]
                   + f_10 * osg_136[k]
                   - f_8 * pc_y[k] * osh1_192[k];

        t_277[k] = f_9 * osg_137[k]
                   + f_3 * pc_y[k] * qsg_197[k];

        t_278[k] = pa_y[k] * osh0_194[k]
                   - f_8 * pc_y[k] * osh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, osh0_195, osh0_198, \
                         osg_123, osg_138, osg_140, osh1_195, osh1_198, qsg_198, \
                         qsg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * osh0_195[k]
                   + f_11 * osg_138[k]
                   - f_8 * pc_y[k] * osh1_195[k];

        t_280[k] = f_11 * osg_123[k]
                   + f_3 * pc_z[k] * qsg_198[k];

        t_281[k] = f_9 * osg_140[k]
                   + f_3 * pc_y[k] * qsg_200[k];

        t_282[k] = pa_y[k] * osh0_198[k]
                   - f_8 * pc_y[k] * osh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, osg_205, osg_206, osg_207, \
                         osg_208, osg_209, qsg_205, qsg_206, qsg_207, qsg_208, \
                         qsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_17 * osg_205[k]
                   + f_3 * pc_x[k] * qsg_205[k];

        t_284[k] = f_17 * osg_206[k]
                   + f_3 * pc_x[k] * qsg_206[k];

        t_285[k] = f_17 * osg_207[k]
                   + f_3 * pc_x[k] * qsg_207[k];

        t_286[k] = f_17 * osg_208[k]
                   + f_3 * pc_x[k] * qsg_208[k];

        t_287[k] = f_17 * osg_209[k]
                   + f_3 * pc_x[k] * qsg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, osg_130, osg_145, osg_147, qsf0_136, \
                         qsf0_138, qsf1_136, qsf1_138, qsg_205, \
                         qsg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * osg_145[k]
                   + f_1 * qsf0_136[k]
                   - f_2 * qsf1_136[k]
                   + f_3 * pc_y[k] * qsg_205[k];

        t_289[k] = f_11 * osg_130[k]
                   + f_3 * pc_z[k] * qsg_205[k];

        t_290[k] = f_9 * osg_147[k]
                   + f_6 * qsf0_138[k]
                   - f_7 * qsf1_138[k]
                   + f_3 * pc_y[k] * qsg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, osh0_209, osg_148, osg_149, \
                         osh1_209, qsf0_139, qsf1_139, qsg_208, \
                         qsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * osg_148[k]
                   + f_4 * qsf0_139[k]
                   - f_5 * qsf1_139[k]
                   + f_3 * pc_y[k] * qsg_208[k];

        t_292[k] = f_9 * osg_149[k]
                   + f_3 * pc_y[k] * qsg_209[k];

        t_293[k] = pa_y[k] * osh0_209[k]
                   - f_8 * pc_y[k] * osh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, osg_135, \
                         osg_210, qsf0_140, qsf1_140, qsg_210, qsg_211, \
                         qsg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_17 * osg_210[k]
                   + f_1 * qsf0_140[k]
                   - f_2 * qsf1_140[k]
                   + f_3 * pc_x[k] * qsg_210[k];

        t_295[k] = f_3 * pc_y[k] * qsg_210[k];

        t_296[k] = f_18 * osg_135[k]
                   + f_3 * pc_z[k] * qsg_210[k];

        t_297[k] = f_4 * qsf0_140[k]
                   - f_5 * qsf1_140[k]
                   + f_3 * pc_y[k] * qsg_211[k];

        t_298[k] = f_3 * pc_y[k] * qsg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, osg_215, qsf0_141, qsf0_142, \
                         qsf0_145, qsf1_141, qsf1_142, qsf1_145, qsg_213, qsg_214, \
                         qsg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_17 * osg_215[k]
                   + f_6 * qsf0_145[k]
                   - f_7 * qsf1_145[k]
                   + f_3 * pc_x[k] * qsg_215[k];

        t_300[k] = f_6 * qsf0_141[k]
                   - f_7 * qsf1_141[k]
                   + f_3 * pc_y[k] * qsg_213[k];

        t_301[k] = f_4 * qsf0_142[k]
                   - f_5 * qsf1_142[k]
                   + f_3 * pc_y[k] * qsg_214[k];

        t_302[k] = f_3 * pc_y[k] * qsg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, osg_219, osg_220, osg_221, osg_222, \
                         qsf0_149, qsf1_149, qsg_219, qsg_220, qsg_221, \
                         qsg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_17 * osg_219[k]
                   + f_4 * qsf0_149[k]
                   - f_5 * qsf1_149[k]
                   + f_3 * pc_x[k] * qsg_219[k];

        t_304[k] = f_17 * osg_220[k]
                   + f_3 * pc_x[k] * qsg_220[k];

        t_305[k] = f_17 * osg_221[k]
                   + f_3 * pc_x[k] * qsg_221[k];

        t_306[k] = f_17 * osg_222[k]
                   + f_3 * pc_x[k] * qsg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, osg_224, qsf0_146, qsf0_147, \
                         qsf1_146, qsf1_147, qsg_219, qsg_220, qsg_221, \
                         qsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * qsg_219[k];

        t_308[k] = f_17 * osg_224[k]
                   + f_3 * pc_x[k] * qsg_224[k];

        t_309[k] = f_1 * qsf0_146[k]
                   - f_2 * qsf1_146[k]
                   + f_3 * pc_y[k] * qsg_220[k];

        t_310[k] = f_13 * qsf0_147[k]
                   - f_14 * qsf1_147[k]
                   + f_3 * pc_y[k] * qsg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, osg_149, qsf0_148, qsf0_149, \
                         qsf1_148, qsf1_149, qsg_222, qsg_223, \
                         qsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * qsf0_148[k]
                   - f_7 * qsf1_148[k]
                   + f_3 * pc_y[k] * qsg_222[k];

        t_312[k] = f_4 * qsf0_149[k]
                   - f_5 * qsf1_149[k]
                   + f_3 * pc_y[k] * qsg_223[k];

        t_313[k] = f_3 * pc_y[k] * qsg_224[k];

        t_314[k] = f_18 * osg_149[k]
                   + f_1 * qsf0_149[k]
                   - f_2 * qsf1_149[k]
                   + f_3 * pc_z[k] * qsg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, osg_150, osg_225, \
                         osg_228, qsf0_150, qsf0_153, qsf1_150, qsf1_153, qsg_225, \
                         qsg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_19 * osg_225[k]
                   + f_1 * qsf0_150[k]
                   - f_2 * qsf1_150[k]
                   + f_3 * pc_x[k] * qsg_225[k];

        t_316[k] = f_20 * osg_150[k]
                   + f_3 * pc_y[k] * qsg_225[k];

        t_317[k] = f_3 * pc_z[k] * qsg_225[k];

        t_318[k] = f_19 * osg_228[k]
                   + f_6 * qsf0_153[k]
                   - f_7 * qsf1_153[k]
                   + f_3 * pc_x[k] * qsg_228[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_z, osg_231, qsf0_150, qsf0_156, \
                         qsf1_150, qsf1_156, qsg_226, qsg_227, qsg_228, \
                         qsg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * qsg_226[k];

        t_320[k] = f_4 * qsf0_150[k]
                   - f_5 * qsf1_150[k]
                   + f_3 * pc_z[k] * qsg_227[k];

        t_321[k] = f_19 * osg_231[k]
                   + f_4 * qsf0_156[k]
                   - f_5 * qsf1_156[k]
                   + f_3 * pc_x[k] * qsg_231[k];

        t_322[k] = f_3 * pc_z[k] * qsg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, osg_155, osg_235, \
                         qsf0_152, qsf1_152, qsg_230, qsg_231, \
                         qsg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_20 * osg_155[k]
                   + f_3 * pc_y[k] * qsg_230[k];

        t_324[k] = f_6 * qsf0_152[k]
                   - f_7 * qsf1_152[k]
                   + f_3 * pc_z[k] * qsg_230[k];

        t_325[k] = f_19 * osg_235[k]
                   + f_3 * pc_x[k] * qsg_235[k];

        t_326[k] = f_3 * pc_z[k] * qsg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, osg_160, osg_237, osg_238, \
                         osg_239, qsf0_156, qsf1_156, qsg_235, qsg_237, qsg_238, \
                         qsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_19 * osg_237[k]
                   + f_3 * pc_x[k] * qsg_237[k];

        t_328[k] = f_19 * osg_238[k]
                   + f_3 * pc_x[k] * qsg_238[k];

        t_329[k] = f_19 * osg_239[k]
                   + f_3 * pc_x[k] * qsg_239[k];

        t_330[k] = f_20 * osg_160[k]
                   + f_1 * qsf0_156[k]
                   - f_2 * qsf1_156[k]
                   + f_3 * pc_y[k] * qsg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, osg_164, qsf0_156, qsf0_157, \
                         qsf1_156, qsf1_157, qsg_235, qsg_236, qsg_237, \
                         qsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * qsg_235[k];

        t_332[k] = f_4 * qsf0_156[k]
                   - f_5 * qsf1_156[k]
                   + f_3 * pc_z[k] * qsg_236[k];

        t_333[k] = f_6 * qsf0_157[k]
                   - f_7 * qsf1_157[k]
                   + f_3 * pc_z[k] * qsg_237[k];

        t_334[k] = f_20 * osg_164[k]
                   + f_3 * pc_y[k] * qsg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pc_y, pc_z, osh0_210, osg_150, \
                         osg_165, osh1_210, qsf0_159, qsf1_159, qsg_239, \
                         qsg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * qsf0_159[k]
                   - f_2 * qsf1_159[k]
                   + f_3 * pc_z[k] * qsg_239[k];

        t_336[k] = pa_z[k] * osh0_210[k]
                   - f_8 * pc_z[k] * osh1_210[k];

        t_337[k] = f_18 * osg_165[k]
                   + f_3 * pc_y[k] * qsg_240[k];

        t_338[k] = f_9 * osg_150[k]
                   + f_3 * pc_z[k] * qsg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_x, pc_y, pc_z, osh0_213, osg_167, \
                         osg_245, osh1_213, qsf0_165, qsf1_165, qsg_242, \
                         qsg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * osh0_213[k]
                   - f_8 * pc_z[k] * osh1_213[k];

        t_340[k] = f_18 * osg_167[k]
                   + f_3 * pc_y[k] * qsg_242[k];

        t_341[k] = f_19 * osg_245[k]
                   + f_6 * qsf0_165[k]
                   - f_7 * qsf1_165[k]
                   + f_3 * pc_x[k] * qsg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, osh0_216, osg_153, osg_170, \
                         osh1_216, qsg_243, qsg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * osh0_216[k]
                   - f_8 * pc_z[k] * osh1_216[k];

        t_343[k] = f_9 * osg_153[k]
                   + f_3 * pc_z[k] * qsg_243[k];

        t_344[k] = f_18 * osg_170[k]
                   + f_3 * pc_y[k] * qsg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, osg_249, osg_250, osg_251, osg_252, \
                         qsf0_169, qsf1_169, qsg_249, qsg_250, qsg_251, \
                         qsg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_19 * osg_249[k]
                   + f_4 * qsf0_169[k]
                   - f_5 * qsf1_169[k]
                   + f_3 * pc_x[k] * qsg_249[k];

        t_346[k] = f_19 * osg_250[k]
                   + f_3 * pc_x[k] * qsg_250[k];

        t_347[k] = f_19 * osg_251[k]
                   + f_3 * pc_x[k] * qsg_251[k];

        t_348[k] = f_19 * osg_252[k]
                   + f_3 * pc_x[k] * qsg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, osh0_225, osg_160, \
                         osg_253, osg_254, osh1_225, qsg_250, qsg_253, \
                         qsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_19 * osg_253[k]
                   + f_3 * pc_x[k] * qsg_253[k];

        t_350[k] = f_19 * osg_254[k]
                   + f_3 * pc_x[k] * qsg_254[k];

        t_351[k] = pa_z[k] * osh0_225[k]
                   - f_8 * pc_z[k] * osh1_225[k];

        t_352[k] = f_9 * osg_160[k]
                   + f_3 * pc_z[k] * qsg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, osg_177, osg_178, osg_179, qsf0_168, \
                         qsf0_169, qsf1_168, qsf1_169, qsg_252, qsg_253, \
                         qsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_18 * osg_177[k]
                   + f_6 * qsf0_168[k]
                   - f_7 * qsf1_168[k]
                   + f_3 * pc_y[k] * qsg_252[k];

        t_354[k] = f_18 * osg_178[k]
                   + f_4 * qsf0_169[k]
                   - f_5 * qsf1_169[k]
                   + f_3 * pc_y[k] * qsg_253[k];

        t_355[k] = f_18 * osg_179[k]
                   + f_3 * pc_y[k] * qsg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, osg_164, osg_180, osg_255, \
                         qsf0_169, qsf0_170, qsf1_169, qsf1_170, qsg_254, \
                         qsg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * osg_164[k]
                   + f_1 * qsf0_169[k]
                   - f_2 * qsf1_169[k]
                   + f_3 * pc_z[k] * qsg_254[k];

        t_357[k] = f_19 * osg_255[k]
                   + f_1 * qsf0_170[k]
                   - f_2 * qsf1_170[k]
                   + f_3 * pc_x[k] * qsg_255[k];

        t_358[k] = f_11 * osg_180[k]
                   + f_3 * pc_y[k] * qsg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, osg_165, osg_182, osg_258, \
                         qsf0_173, qsf1_173, qsg_255, qsg_257, \
                         qsg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * osg_165[k]
                   + f_3 * pc_z[k] * qsg_255[k];

        t_360[k] = f_19 * osg_258[k]
                   + f_6 * qsf0_173[k]
                   - f_7 * qsf1_173[k]
                   + f_3 * pc_x[k] * qsg_258[k];

        t_361[k] = f_11 * osg_182[k]
                   + f_3 * pc_y[k] * qsg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, osg_168, osg_260, osg_261, qsf0_175, \
                         qsf0_176, qsf1_175, qsf1_176, qsg_258, qsg_260, \
                         qsg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_19 * osg_260[k]
                   + f_6 * qsf0_175[k]
                   - f_7 * qsf1_175[k]
                   + f_3 * pc_x[k] * qsg_260[k];

        t_363[k] = f_19 * osg_261[k]
                   + f_4 * qsf0_176[k]
                   - f_5 * qsf1_176[k]
                   + f_3 * pc_x[k] * qsg_261[k];

        t_364[k] = f_10 * osg_168[k]
                   + f_3 * pc_z[k] * qsg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, osg_185, osg_264, osg_265, \
                         osg_266, qsf0_179, qsf1_179, qsg_260, qsg_264, qsg_265, \
                         qsg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * osg_185[k]
                   + f_3 * pc_y[k] * qsg_260[k];

        t_366[k] = f_19 * osg_264[k]
                   + f_4 * qsf0_179[k]
                   - f_5 * qsf1_179[k]
                   + f_3 * pc_x[k] * qsg_264[k];

        t_367[k] = f_19 * osg_265[k]
                   + f_3 * pc_x[k] * qsg_265[k];

        t_368[k] = f_19 * osg_266[k]
                   + f_3 * pc_x[k] * qsg_266[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *osh0_294 = buffer.data(osh0 + 294);
    const auto *osh0_297 = buffer.data(osh0 + 297);
    const auto *osh0_299 = buffer.data(osh0 + 299);
    const auto *osh0_300 = buffer.data(osh0 + 300);
    const auto *osh0_303 = buffer.data(osh0 + 303);
    const auto *osh0_314 = buffer.data(osh0 + 314);
    const auto *osh0_315 = buffer.data(osh0 + 315);
    const auto *osh0_318 = buffer.data(osh0 + 318);
    const auto *osh0_321 = buffer.data(osh0 + 321);
    const auto *osh0_330 = buffer.data(osh0 + 330);

    const auto *osg_175 = buffer.data(osg + 175);
    const auto *osg_179 = buffer.data(osg + 179);
    const auto *osg_180 = buffer.data(osg + 180);
    const auto *osg_183 = buffer.data(osg + 183);
    const auto *osg_190 = buffer.data(osg + 190);
    const auto *osg_192 = buffer.data(osg + 192);
    const auto *osg_193 = buffer.data(osg + 193);
    const auto *osg_194 = buffer.data(osg + 194);
    const auto *osg_195 = buffer.data(osg + 195);
    const auto *osg_197 = buffer.data(osg + 197);
    const auto *osg_198 = buffer.data(osg + 198);
    const auto *osg_200 = buffer.data(osg + 200);
    const auto *osg_205 = buffer.data(osg + 205);
    const auto *osg_207 = buffer.data(osg + 207);
    const auto *osg_208 = buffer.data(osg + 208);
    const auto *osg_209 = buffer.data(osg + 209);
    const auto *osg_210 = buffer.data(osg + 210);
    const auto *osg_211 = buffer.data(osg + 211);
    const auto *osg_212 = buffer.data(osg + 212);
    const auto *osg_213 = buffer.data(osg + 213);
    const auto *osg_215 = buffer.data(osg + 215);
    const auto *osg_220 = buffer.data(osg + 220);
    const auto *osg_222 = buffer.data(osg + 222);
    const auto *osg_223 = buffer.data(osg + 223);
    const auto *osg_224 = buffer.data(osg + 224);
    const auto *osg_225 = buffer.data(osg + 225);
    const auto *osg_228 = buffer.data(osg + 228);
    const auto *osg_230 = buffer.data(osg + 230);
    const auto *osg_235 = buffer.data(osg + 235);
    const auto *osg_239 = buffer.data(osg + 239);
    const auto *osg_240 = buffer.data(osg + 240);
    const auto *osg_242 = buffer.data(osg + 242);
    const auto *osg_245 = buffer.data(osg + 245);
    const auto *osg_252 = buffer.data(osg + 252);
    const auto *osg_253 = buffer.data(osg + 253);
    const auto *osg_254 = buffer.data(osg + 254);
    const auto *osg_255 = buffer.data(osg + 255);
    const auto *osg_257 = buffer.data(osg + 257);
    const auto *osg_267 = buffer.data(osg + 267);
    const auto *osg_268 = buffer.data(osg + 268);
    const auto *osg_269 = buffer.data(osg + 269);
    const auto *osg_270 = buffer.data(osg + 270);
    const auto *osg_273 = buffer.data(osg + 273);
    const auto *osg_275 = buffer.data(osg + 275);
    const auto *osg_276 = buffer.data(osg + 276);
    const auto *osg_279 = buffer.data(osg + 279);
    const auto *osg_280 = buffer.data(osg + 280);
    const auto *osg_281 = buffer.data(osg + 281);
    const auto *osg_282 = buffer.data(osg + 282);
    const auto *osg_283 = buffer.data(osg + 283);
    const auto *osg_284 = buffer.data(osg + 284);
    const auto *osg_295 = buffer.data(osg + 295);
    const auto *osg_296 = buffer.data(osg + 296);
    const auto *osg_297 = buffer.data(osg + 297);
    const auto *osg_298 = buffer.data(osg + 298);
    const auto *osg_299 = buffer.data(osg + 299);
    const auto *osg_300 = buffer.data(osg + 300);
    const auto *osg_305 = buffer.data(osg + 305);
    const auto *osg_309 = buffer.data(osg + 309);
    const auto *osg_310 = buffer.data(osg + 310);
    const auto *osg_311 = buffer.data(osg + 311);
    const auto *osg_312 = buffer.data(osg + 312);
    const auto *osg_314 = buffer.data(osg + 314);
    const auto *osg_315 = buffer.data(osg + 315);
    const auto *osg_318 = buffer.data(osg + 318);
    const auto *osg_321 = buffer.data(osg + 321);
    const auto *osg_325 = buffer.data(osg + 325);
    const auto *osg_327 = buffer.data(osg + 327);
    const auto *osg_328 = buffer.data(osg + 328);
    const auto *osg_329 = buffer.data(osg + 329);
    const auto *osg_335 = buffer.data(osg + 335);
    const auto *osg_339 = buffer.data(osg + 339);
    const auto *osg_340 = buffer.data(osg + 340);
    const auto *osg_341 = buffer.data(osg + 341);
    const auto *osg_342 = buffer.data(osg + 342);
    const auto *osg_343 = buffer.data(osg + 343);
    const auto *osg_344 = buffer.data(osg + 344);
    const auto *osg_345 = buffer.data(osg + 345);
    const auto *osg_348 = buffer.data(osg + 348);

    const auto *osh1_294 = buffer.data(osh1 + 294);
    const auto *osh1_297 = buffer.data(osh1 + 297);
    const auto *osh1_299 = buffer.data(osh1 + 299);
    const auto *osh1_300 = buffer.data(osh1 + 300);
    const auto *osh1_303 = buffer.data(osh1 + 303);
    const auto *osh1_314 = buffer.data(osh1 + 314);
    const auto *osh1_315 = buffer.data(osh1 + 315);
    const auto *osh1_318 = buffer.data(osh1 + 318);
    const auto *osh1_321 = buffer.data(osh1 + 321);
    const auto *osh1_330 = buffer.data(osh1 + 330);

    const auto *qsf0_176 = buffer.data(qsf0 + 176);
    const auto *qsf0_178 = buffer.data(qsf0 + 178);
    const auto *qsf0_179 = buffer.data(qsf0 + 179);
    const auto *qsf0_180 = buffer.data(qsf0 + 180);
    const auto *qsf0_183 = buffer.data(qsf0 + 183);
    const auto *qsf0_185 = buffer.data(qsf0 + 185);
    const auto *qsf0_186 = buffer.data(qsf0 + 186);
    const auto *qsf0_188 = buffer.data(qsf0 + 188);
    const auto *qsf0_189 = buffer.data(qsf0 + 189);
    const auto *qsf0_196 = buffer.data(qsf0 + 196);
    const auto *qsf0_198 = buffer.data(qsf0 + 198);
    const auto *qsf0_199 = buffer.data(qsf0 + 199);
    const auto *qsf0_200 = buffer.data(qsf0 + 200);
    const auto *qsf0_201 = buffer.data(qsf0 + 201);
    const auto *qsf0_202 = buffer.data(qsf0 + 202);
    const auto *qsf0_205 = buffer.data(qsf0 + 205);
    const auto *qsf0_206 = buffer.data(qsf0 + 206);
    const auto *qsf0_207 = buffer.data(qsf0 + 207);
    const auto *qsf0_208 = buffer.data(qsf0 + 208);
    const auto *qsf0_209 = buffer.data(qsf0 + 209);
    const auto *qsf0_210 = buffer.data(qsf0 + 210);
    const auto *qsf0_212 = buffer.data(qsf0 + 212);
    const auto *qsf0_213 = buffer.data(qsf0 + 213);
    const auto *qsf0_216 = buffer.data(qsf0 + 216);
    const auto *qsf0_217 = buffer.data(qsf0 + 217);
    const auto *qsf0_219 = buffer.data(qsf0 + 219);
    const auto *qsf0_225 = buffer.data(qsf0 + 225);
    const auto *qsf0_228 = buffer.data(qsf0 + 228);
    const auto *qsf0_229 = buffer.data(qsf0 + 229);
    const auto *qsf0_230 = buffer.data(qsf0 + 230);
    const auto *qsf0_233 = buffer.data(qsf0 + 233);

    const auto *qsf1_176 = buffer.data(qsf1 + 176);
    const auto *qsf1_178 = buffer.data(qsf1 + 178);
    const auto *qsf1_179 = buffer.data(qsf1 + 179);
    const auto *qsf1_180 = buffer.data(qsf1 + 180);
    const auto *qsf1_183 = buffer.data(qsf1 + 183);
    const auto *qsf1_185 = buffer.data(qsf1 + 185);
    const auto *qsf1_186 = buffer.data(qsf1 + 186);
    const auto *qsf1_188 = buffer.data(qsf1 + 188);
    const auto *qsf1_189 = buffer.data(qsf1 + 189);
    const auto *qsf1_196 = buffer.data(qsf1 + 196);
    const auto *qsf1_198 = buffer.data(qsf1 + 198);
    const auto *qsf1_199 = buffer.data(qsf1 + 199);
    const auto *qsf1_200 = buffer.data(qsf1 + 200);
    const auto *qsf1_201 = buffer.data(qsf1 + 201);
    const auto *qsf1_202 = buffer.data(qsf1 + 202);
    const auto *qsf1_205 = buffer.data(qsf1 + 205);
    const auto *qsf1_206 = buffer.data(qsf1 + 206);
    const auto *qsf1_207 = buffer.data(qsf1 + 207);
    const auto *qsf1_208 = buffer.data(qsf1 + 208);
    const auto *qsf1_209 = buffer.data(qsf1 + 209);
    const auto *qsf1_210 = buffer.data(qsf1 + 210);
    const auto *qsf1_212 = buffer.data(qsf1 + 212);
    const auto *qsf1_213 = buffer.data(qsf1 + 213);
    const auto *qsf1_216 = buffer.data(qsf1 + 216);
    const auto *qsf1_217 = buffer.data(qsf1 + 217);
    const auto *qsf1_219 = buffer.data(qsf1 + 219);
    const auto *qsf1_225 = buffer.data(qsf1 + 225);
    const auto *qsf1_228 = buffer.data(qsf1 + 228);
    const auto *qsf1_229 = buffer.data(qsf1 + 229);
    const auto *qsf1_230 = buffer.data(qsf1 + 230);
    const auto *qsf1_233 = buffer.data(qsf1 + 233);

    const auto *qsg_265 = buffer.data(qsg + 265);
    const auto *qsg_267 = buffer.data(qsg + 267);
    const auto *qsg_268 = buffer.data(qsg + 268);
    const auto *qsg_269 = buffer.data(qsg + 269);
    const auto *qsg_270 = buffer.data(qsg + 270);
    const auto *qsg_272 = buffer.data(qsg + 272);
    const auto *qsg_273 = buffer.data(qsg + 273);
    const auto *qsg_275 = buffer.data(qsg + 275);
    const auto *qsg_276 = buffer.data(qsg + 276);
    const auto *qsg_279 = buffer.data(qsg + 279);
    const auto *qsg_280 = buffer.data(qsg + 280);
    const auto *qsg_281 = buffer.data(qsg + 281);
    const auto *qsg_282 = buffer.data(qsg + 282);
    const auto *qsg_283 = buffer.data(qsg + 283);
    const auto *qsg_284 = buffer.data(qsg + 284);
    const auto *qsg_285 = buffer.data(qsg + 285);
    const auto *qsg_287 = buffer.data(qsg + 287);
    const auto *qsg_288 = buffer.data(qsg + 288);
    const auto *qsg_290 = buffer.data(qsg + 290);
    const auto *qsg_295 = buffer.data(qsg + 295);
    const auto *qsg_296 = buffer.data(qsg + 296);
    const auto *qsg_297 = buffer.data(qsg + 297);
    const auto *qsg_298 = buffer.data(qsg + 298);
    const auto *qsg_299 = buffer.data(qsg + 299);
    const auto *qsg_300 = buffer.data(qsg + 300);
    const auto *qsg_301 = buffer.data(qsg + 301);
    const auto *qsg_302 = buffer.data(qsg + 302);
    const auto *qsg_303 = buffer.data(qsg + 303);
    const auto *qsg_304 = buffer.data(qsg + 304);
    const auto *qsg_305 = buffer.data(qsg + 305);
    const auto *qsg_309 = buffer.data(qsg + 309);
    const auto *qsg_310 = buffer.data(qsg + 310);
    const auto *qsg_311 = buffer.data(qsg + 311);
    const auto *qsg_312 = buffer.data(qsg + 312);
    const auto *qsg_313 = buffer.data(qsg + 313);
    const auto *qsg_314 = buffer.data(qsg + 314);
    const auto *qsg_315 = buffer.data(qsg + 315);
    const auto *qsg_316 = buffer.data(qsg + 316);
    const auto *qsg_317 = buffer.data(qsg + 317);
    const auto *qsg_318 = buffer.data(qsg + 318);
    const auto *qsg_320 = buffer.data(qsg + 320);
    const auto *qsg_321 = buffer.data(qsg + 321);
    const auto *qsg_325 = buffer.data(qsg + 325);
    const auto *qsg_326 = buffer.data(qsg + 326);
    const auto *qsg_327 = buffer.data(qsg + 327);
    const auto *qsg_328 = buffer.data(qsg + 328);
    const auto *qsg_329 = buffer.data(qsg + 329);
    const auto *qsg_330 = buffer.data(qsg + 330);
    const auto *qsg_332 = buffer.data(qsg + 332);
    const auto *qsg_333 = buffer.data(qsg + 333);
    const auto *qsg_335 = buffer.data(qsg + 335);
    const auto *qsg_339 = buffer.data(qsg + 339);
    const auto *qsg_340 = buffer.data(qsg + 340);
    const auto *qsg_341 = buffer.data(qsg + 341);
    const auto *qsg_342 = buffer.data(qsg + 342);
    const auto *qsg_343 = buffer.data(qsg + 343);
    const auto *qsg_344 = buffer.data(qsg + 344);
    const auto *qsg_345 = buffer.data(qsg + 345);
    const auto *qsg_347 = buffer.data(qsg + 347);
    const auto *qsg_348 = buffer.data(qsg + 348);

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, osg_190, osg_267, osg_268, \
                         osg_269, qsf0_176, qsf1_176, qsg_265, qsg_267, qsg_268, \
                         qsg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_19 * osg_267[k]
                   + f_3 * pc_x[k] * qsg_267[k];

        t_370[k] = f_19 * osg_268[k]
                   + f_3 * pc_x[k] * qsg_268[k];

        t_371[k] = f_19 * osg_269[k]
                   + f_3 * pc_x[k] * qsg_269[k];

        t_372[k] = f_11 * osg_190[k]
                   + f_1 * qsf0_176[k]
                   - f_2 * qsf1_176[k]
                   + f_3 * pc_y[k] * qsg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, osg_175, osg_192, osg_193, qsf0_178, \
                         qsf0_179, qsf1_178, qsf1_179, qsg_265, qsg_267, \
                         qsg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * osg_175[k]
                   + f_3 * pc_z[k] * qsg_265[k];

        t_374[k] = f_11 * osg_192[k]
                   + f_6 * qsf0_178[k]
                   - f_7 * qsf1_178[k]
                   + f_3 * pc_y[k] * qsg_267[k];

        t_375[k] = f_11 * osg_193[k]
                   + f_4 * qsf0_179[k]
                   - f_5 * qsf1_179[k]
                   + f_3 * pc_y[k] * qsg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, osg_179, osg_194, osg_270, \
                         qsf0_179, qsf0_180, qsf1_179, qsf1_180, qsg_269, \
                         qsg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * osg_194[k]
                   + f_3 * pc_y[k] * qsg_269[k];

        t_377[k] = f_10 * osg_179[k]
                   + f_1 * qsf0_179[k]
                   - f_2 * qsf1_179[k]
                   + f_3 * pc_z[k] * qsg_269[k];

        t_378[k] = f_19 * osg_270[k]
                   + f_1 * qsf0_180[k]
                   - f_2 * qsf1_180[k]
                   + f_3 * pc_x[k] * qsg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, osg_180, osg_195, \
                         osg_197, osg_273, qsf0_183, qsf1_183, qsg_270, qsg_272, \
                         qsg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * osg_195[k]
                   + f_3 * pc_y[k] * qsg_270[k];

        t_380[k] = f_11 * osg_180[k]
                   + f_3 * pc_z[k] * qsg_270[k];

        t_381[k] = f_19 * osg_273[k]
                   + f_6 * qsf0_183[k]
                   - f_7 * qsf1_183[k]
                   + f_3 * pc_x[k] * qsg_273[k];

        t_382[k] = f_10 * osg_197[k]
                   + f_3 * pc_y[k] * qsg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, osg_183, osg_275, osg_276, qsf0_185, \
                         qsf0_186, qsf1_185, qsf1_186, qsg_273, qsg_275, \
                         qsg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_19 * osg_275[k]
                   + f_6 * qsf0_185[k]
                   - f_7 * qsf1_185[k]
                   + f_3 * pc_x[k] * qsg_275[k];

        t_384[k] = f_19 * osg_276[k]
                   + f_4 * qsf0_186[k]
                   - f_5 * qsf1_186[k]
                   + f_3 * pc_x[k] * qsg_276[k];

        t_385[k] = f_11 * osg_183[k]
                   + f_3 * pc_z[k] * qsg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, osg_200, osg_279, osg_280, \
                         osg_281, qsf0_189, qsf1_189, qsg_275, qsg_279, qsg_280, \
                         qsg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * osg_200[k]
                   + f_3 * pc_y[k] * qsg_275[k];

        t_387[k] = f_19 * osg_279[k]
                   + f_4 * qsf0_189[k]
                   - f_5 * qsf1_189[k]
                   + f_3 * pc_x[k] * qsg_279[k];

        t_388[k] = f_19 * osg_280[k]
                   + f_3 * pc_x[k] * qsg_280[k];

        t_389[k] = f_19 * osg_281[k]
                   + f_3 * pc_x[k] * qsg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, osg_205, osg_282, osg_283, \
                         osg_284, qsf0_186, qsf1_186, qsg_280, qsg_282, qsg_283, \
                         qsg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_19 * osg_282[k]
                   + f_3 * pc_x[k] * qsg_282[k];

        t_391[k] = f_19 * osg_283[k]
                   + f_3 * pc_x[k] * qsg_283[k];

        t_392[k] = f_19 * osg_284[k]
                   + f_3 * pc_x[k] * qsg_284[k];

        t_393[k] = f_10 * osg_205[k]
                   + f_1 * qsf0_186[k]
                   - f_2 * qsf1_186[k]
                   + f_3 * pc_y[k] * qsg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, osg_190, osg_207, osg_208, qsf0_188, \
                         qsf0_189, qsf1_188, qsf1_189, qsg_280, qsg_282, \
                         qsg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * osg_190[k]
                   + f_3 * pc_z[k] * qsg_280[k];

        t_395[k] = f_10 * osg_207[k]
                   + f_6 * qsf0_188[k]
                   - f_7 * qsf1_188[k]
                   + f_3 * pc_y[k] * qsg_282[k];

        t_396[k] = f_10 * osg_208[k]
                   + f_4 * qsf0_189[k]
                   - f_5 * qsf1_189[k]
                   + f_3 * pc_y[k] * qsg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pc_y, pc_z, osh0_294, osg_194, \
                         osg_209, osg_210, osh1_294, qsf0_189, qsf1_189, qsg_284, \
                         qsg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * osg_209[k]
                   + f_3 * pc_y[k] * qsg_284[k];

        t_398[k] = f_11 * osg_194[k]
                   + f_1 * qsf0_189[k]
                   - f_2 * qsf1_189[k]
                   + f_3 * pc_z[k] * qsg_284[k];

        t_399[k] = pa_y[k] * osh0_294[k]
                   - f_8 * pc_y[k] * osh1_294[k];

        t_400[k] = f_9 * osg_210[k]
                   + f_3 * pc_y[k] * qsg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pc_y, pc_z, osh0_297, osh0_299, \
                         osg_195, osg_211, osg_212, osh1_297, osh1_299, qsg_285, \
                         qsg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * osg_195[k]
                   + f_3 * pc_z[k] * qsg_285[k];

        t_402[k] = pa_y[k] * osh0_297[k]
                   + f_10 * osg_211[k]
                   - f_8 * pc_y[k] * osh1_297[k];

        t_403[k] = f_9 * osg_212[k]
                   + f_3 * pc_y[k] * qsg_287[k];

        t_404[k] = pa_y[k] * osh0_299[k]
                   - f_8 * pc_y[k] * osh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_y, pc_y, pc_z, osh0_300, osh0_303, \
                         osg_198, osg_213, osg_215, osh1_300, osh1_303, qsg_288, \
                         qsg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_y[k] * osh0_300[k]
                   + f_11 * osg_213[k]
                   - f_8 * pc_y[k] * osh1_300[k];

        t_406[k] = f_18 * osg_198[k]
                   + f_3 * pc_z[k] * qsg_288[k];

        t_407[k] = f_9 * osg_215[k]
                   + f_3 * pc_y[k] * qsg_290[k];

        t_408[k] = pa_y[k] * osh0_303[k]
                   - f_8 * pc_y[k] * osh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, osg_295, osg_296, osg_297, \
                         osg_298, osg_299, qsg_295, qsg_296, qsg_297, qsg_298, \
                         qsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_19 * osg_295[k]
                   + f_3 * pc_x[k] * qsg_295[k];

        t_410[k] = f_19 * osg_296[k]
                   + f_3 * pc_x[k] * qsg_296[k];

        t_411[k] = f_19 * osg_297[k]
                   + f_3 * pc_x[k] * qsg_297[k];

        t_412[k] = f_19 * osg_298[k]
                   + f_3 * pc_x[k] * qsg_298[k];

        t_413[k] = f_19 * osg_299[k]
                   + f_3 * pc_x[k] * qsg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, osg_205, osg_220, osg_222, qsf0_196, \
                         qsf0_198, qsf1_196, qsf1_198, qsg_295, \
                         qsg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * osg_220[k]
                   + f_1 * qsf0_196[k]
                   - f_2 * qsf1_196[k]
                   + f_3 * pc_y[k] * qsg_295[k];

        t_415[k] = f_18 * osg_205[k]
                   + f_3 * pc_z[k] * qsg_295[k];

        t_416[k] = f_9 * osg_222[k]
                   + f_6 * qsf0_198[k]
                   - f_7 * qsf1_198[k]
                   + f_3 * pc_y[k] * qsg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, osh0_314, osg_223, osg_224, \
                         osh1_314, qsf0_199, qsf1_199, qsg_298, \
                         qsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * osg_223[k]
                   + f_4 * qsf0_199[k]
                   - f_5 * qsf1_199[k]
                   + f_3 * pc_y[k] * qsg_298[k];

        t_418[k] = f_9 * osg_224[k]
                   + f_3 * pc_y[k] * qsg_299[k];

        t_419[k] = pa_y[k] * osh0_314[k]
                   - f_8 * pc_y[k] * osh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, osg_210, \
                         osg_300, qsf0_200, qsf1_200, qsg_300, qsg_301, \
                         qsg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_19 * osg_300[k]
                   + f_1 * qsf0_200[k]
                   - f_2 * qsf1_200[k]
                   + f_3 * pc_x[k] * qsg_300[k];

        t_421[k] = f_3 * pc_y[k] * qsg_300[k];

        t_422[k] = f_20 * osg_210[k]
                   + f_3 * pc_z[k] * qsg_300[k];

        t_423[k] = f_4 * qsf0_200[k]
                   - f_5 * qsf1_200[k]
                   + f_3 * pc_y[k] * qsg_301[k];

        t_424[k] = f_3 * pc_y[k] * qsg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, osg_305, qsf0_201, qsf0_202, \
                         qsf0_205, qsf1_201, qsf1_202, qsf1_205, qsg_303, qsg_304, \
                         qsg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_19 * osg_305[k]
                   + f_6 * qsf0_205[k]
                   - f_7 * qsf1_205[k]
                   + f_3 * pc_x[k] * qsg_305[k];

        t_426[k] = f_6 * qsf0_201[k]
                   - f_7 * qsf1_201[k]
                   + f_3 * pc_y[k] * qsg_303[k];

        t_427[k] = f_4 * qsf0_202[k]
                   - f_5 * qsf1_202[k]
                   + f_3 * pc_y[k] * qsg_304[k];

        t_428[k] = f_3 * pc_y[k] * qsg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, osg_309, osg_310, osg_311, osg_312, \
                         qsf0_209, qsf1_209, qsg_309, qsg_310, qsg_311, \
                         qsg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_19 * osg_309[k]
                   + f_4 * qsf0_209[k]
                   - f_5 * qsf1_209[k]
                   + f_3 * pc_x[k] * qsg_309[k];

        t_430[k] = f_19 * osg_310[k]
                   + f_3 * pc_x[k] * qsg_310[k];

        t_431[k] = f_19 * osg_311[k]
                   + f_3 * pc_x[k] * qsg_311[k];

        t_432[k] = f_19 * osg_312[k]
                   + f_3 * pc_x[k] * qsg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, osg_314, qsf0_206, qsf0_207, \
                         qsf1_206, qsf1_207, qsg_309, qsg_310, qsg_311, \
                         qsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_y[k] * qsg_309[k];

        t_434[k] = f_19 * osg_314[k]
                   + f_3 * pc_x[k] * qsg_314[k];

        t_435[k] = f_1 * qsf0_206[k]
                   - f_2 * qsf1_206[k]
                   + f_3 * pc_y[k] * qsg_310[k];

        t_436[k] = f_13 * qsf0_207[k]
                   - f_14 * qsf1_207[k]
                   + f_3 * pc_y[k] * qsg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_y, pc_z, osg_224, qsf0_208, qsf0_209, \
                         qsf1_208, qsf1_209, qsg_312, qsg_313, \
                         qsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_6 * qsf0_208[k]
                   - f_7 * qsf1_208[k]
                   + f_3 * pc_y[k] * qsg_312[k];

        t_438[k] = f_4 * qsf0_209[k]
                   - f_5 * qsf1_209[k]
                   + f_3 * pc_y[k] * qsg_313[k];

        t_439[k] = f_3 * pc_y[k] * qsg_314[k];

        t_440[k] = f_20 * osg_224[k]
                   + f_1 * qsf0_209[k]
                   - f_2 * qsf1_209[k]
                   + f_3 * pc_z[k] * qsg_314[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, osg_225, osg_315, \
                         osg_318, qsf0_210, qsf0_213, qsf1_210, qsf1_213, qsg_315, \
                         qsg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_21 * osg_315[k]
                   + f_1 * qsf0_210[k]
                   - f_2 * qsf1_210[k]
                   + f_3 * pc_x[k] * qsg_315[k];

        t_442[k] = f_21 * osg_225[k]
                   + f_3 * pc_y[k] * qsg_315[k];

        t_443[k] = f_3 * pc_z[k] * qsg_315[k];

        t_444[k] = f_21 * osg_318[k]
                   + f_6 * qsf0_213[k]
                   - f_7 * qsf1_213[k]
                   + f_3 * pc_x[k] * qsg_318[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_z, osg_321, qsf0_210, qsf0_216, \
                         qsf1_210, qsf1_216, qsg_316, qsg_317, qsg_318, \
                         qsg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_z[k] * qsg_316[k];

        t_446[k] = f_4 * qsf0_210[k]
                   - f_5 * qsf1_210[k]
                   + f_3 * pc_z[k] * qsg_317[k];

        t_447[k] = f_21 * osg_321[k]
                   + f_4 * qsf0_216[k]
                   - f_5 * qsf1_216[k]
                   + f_3 * pc_x[k] * qsg_321[k];

        t_448[k] = f_3 * pc_z[k] * qsg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, osg_230, osg_325, \
                         qsf0_212, qsf1_212, qsg_320, qsg_321, \
                         qsg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_21 * osg_230[k]
                   + f_3 * pc_y[k] * qsg_320[k];

        t_450[k] = f_6 * qsf0_212[k]
                   - f_7 * qsf1_212[k]
                   + f_3 * pc_z[k] * qsg_320[k];

        t_451[k] = f_21 * osg_325[k]
                   + f_3 * pc_x[k] * qsg_325[k];

        t_452[k] = f_3 * pc_z[k] * qsg_321[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, osg_235, osg_327, osg_328, \
                         osg_329, qsf0_216, qsf1_216, qsg_325, qsg_327, qsg_328, \
                         qsg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_21 * osg_327[k]
                   + f_3 * pc_x[k] * qsg_327[k];

        t_454[k] = f_21 * osg_328[k]
                   + f_3 * pc_x[k] * qsg_328[k];

        t_455[k] = f_21 * osg_329[k]
                   + f_3 * pc_x[k] * qsg_329[k];

        t_456[k] = f_21 * osg_235[k]
                   + f_1 * qsf0_216[k]
                   - f_2 * qsf1_216[k]
                   + f_3 * pc_y[k] * qsg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_y, pc_z, osg_239, qsf0_216, qsf0_217, \
                         qsf1_216, qsf1_217, qsg_325, qsg_326, qsg_327, \
                         qsg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * qsg_325[k];

        t_458[k] = f_4 * qsf0_216[k]
                   - f_5 * qsf1_216[k]
                   + f_3 * pc_z[k] * qsg_326[k];

        t_459[k] = f_6 * qsf0_217[k]
                   - f_7 * qsf1_217[k]
                   + f_3 * pc_z[k] * qsg_327[k];

        t_460[k] = f_21 * osg_239[k]
                   + f_3 * pc_y[k] * qsg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_z, pc_y, pc_z, osh0_315, osg_225, \
                         osg_240, osh1_315, qsf0_219, qsf1_219, qsg_329, \
                         qsg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_1 * qsf0_219[k]
                   - f_2 * qsf1_219[k]
                   + f_3 * pc_z[k] * qsg_329[k];

        t_462[k] = pa_z[k] * osh0_315[k]
                   - f_8 * pc_z[k] * osh1_315[k];

        t_463[k] = f_20 * osg_240[k]
                   + f_3 * pc_y[k] * qsg_330[k];

        t_464[k] = f_9 * osg_225[k]
                   + f_3 * pc_z[k] * qsg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_z, pc_x, pc_y, pc_z, osh0_318, osg_242, \
                         osg_335, osh1_318, qsf0_225, qsf1_225, qsg_332, \
                         qsg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * osh0_318[k]
                   - f_8 * pc_z[k] * osh1_318[k];

        t_466[k] = f_20 * osg_242[k]
                   + f_3 * pc_y[k] * qsg_332[k];

        t_467[k] = f_21 * osg_335[k]
                   + f_6 * qsf0_225[k]
                   - f_7 * qsf1_225[k]
                   + f_3 * pc_x[k] * qsg_335[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_z, pc_y, pc_z, osh0_321, osg_228, osg_245, \
                         osh1_321, qsg_333, qsg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_z[k] * osh0_321[k]
                   - f_8 * pc_z[k] * osh1_321[k];

        t_469[k] = f_9 * osg_228[k]
                   + f_3 * pc_z[k] * qsg_333[k];

        t_470[k] = f_20 * osg_245[k]
                   + f_3 * pc_y[k] * qsg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, osg_339, osg_340, osg_341, osg_342, \
                         qsf0_229, qsf1_229, qsg_339, qsg_340, qsg_341, \
                         qsg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_21 * osg_339[k]
                   + f_4 * qsf0_229[k]
                   - f_5 * qsf1_229[k]
                   + f_3 * pc_x[k] * qsg_339[k];

        t_472[k] = f_21 * osg_340[k]
                   + f_3 * pc_x[k] * qsg_340[k];

        t_473[k] = f_21 * osg_341[k]
                   + f_3 * pc_x[k] * qsg_341[k];

        t_474[k] = f_21 * osg_342[k]
                   + f_3 * pc_x[k] * qsg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pc_x, pc_z, osh0_330, osg_235, \
                         osg_343, osg_344, osh1_330, qsg_340, qsg_343, \
                         qsg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_21 * osg_343[k]
                   + f_3 * pc_x[k] * qsg_343[k];

        t_476[k] = f_21 * osg_344[k]
                   + f_3 * pc_x[k] * qsg_344[k];

        t_477[k] = pa_z[k] * osh0_330[k]
                   - f_8 * pc_z[k] * osh1_330[k];

        t_478[k] = f_9 * osg_235[k]
                   + f_3 * pc_z[k] * qsg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_y, osg_252, osg_253, osg_254, qsf0_228, \
                         qsf0_229, qsf1_228, qsf1_229, qsg_342, qsg_343, \
                         qsg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_20 * osg_252[k]
                   + f_6 * qsf0_228[k]
                   - f_7 * qsf1_228[k]
                   + f_3 * pc_y[k] * qsg_342[k];

        t_480[k] = f_20 * osg_253[k]
                   + f_4 * qsf0_229[k]
                   - f_5 * qsf1_229[k]
                   + f_3 * pc_y[k] * qsg_343[k];

        t_481[k] = f_20 * osg_254[k]
                   + f_3 * pc_y[k] * qsg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, osg_239, osg_255, osg_345, \
                         qsf0_229, qsf0_230, qsf1_229, qsf1_230, qsg_344, \
                         qsg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * osg_239[k]
                   + f_1 * qsf0_229[k]
                   - f_2 * qsf1_229[k]
                   + f_3 * pc_z[k] * qsg_344[k];

        t_483[k] = f_21 * osg_345[k]
                   + f_1 * qsf0_230[k]
                   - f_2 * qsf1_230[k]
                   + f_3 * pc_x[k] * qsg_345[k];

        t_484[k] = f_18 * osg_255[k]
                   + f_3 * pc_y[k] * qsg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, osg_240, osg_257, osg_348, \
                         qsf0_233, qsf1_233, qsg_345, qsg_347, \
                         qsg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * osg_240[k]
                   + f_3 * pc_z[k] * qsg_345[k];

        t_486[k] = f_21 * osg_348[k]
                   + f_6 * qsf0_233[k]
                   - f_7 * qsf1_233[k]
                   + f_3 * pc_x[k] * qsg_348[k];

        t_487[k] = f_18 * osg_257[k]
                   + f_3 * pc_y[k] * qsg_347[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *osh0_420 = buffer.data(osh0 + 420);
    const auto *osh0_423 = buffer.data(osh0 + 423);
    const auto *osh0_425 = buffer.data(osh0 + 425);
    const auto *osh0_426 = buffer.data(osh0 + 426);
    const auto *osh0_429 = buffer.data(osh0 + 429);
    const auto *osh0_440 = buffer.data(osh0 + 440);

    const auto *osg_243 = buffer.data(osg + 243);
    const auto *osg_250 = buffer.data(osg + 250);
    const auto *osg_254 = buffer.data(osg + 254);
    const auto *osg_255 = buffer.data(osg + 255);
    const auto *osg_258 = buffer.data(osg + 258);
    const auto *osg_260 = buffer.data(osg + 260);
    const auto *osg_265 = buffer.data(osg + 265);
    const auto *osg_267 = buffer.data(osg + 267);
    const auto *osg_268 = buffer.data(osg + 268);
    const auto *osg_269 = buffer.data(osg + 269);
    const auto *osg_270 = buffer.data(osg + 270);
    const auto *osg_272 = buffer.data(osg + 272);
    const auto *osg_273 = buffer.data(osg + 273);
    const auto *osg_275 = buffer.data(osg + 275);
    const auto *osg_280 = buffer.data(osg + 280);
    const auto *osg_282 = buffer.data(osg + 282);
    const auto *osg_283 = buffer.data(osg + 283);
    const auto *osg_284 = buffer.data(osg + 284);
    const auto *osg_285 = buffer.data(osg + 285);
    const auto *osg_287 = buffer.data(osg + 287);
    const auto *osg_288 = buffer.data(osg + 288);
    const auto *osg_290 = buffer.data(osg + 290);
    const auto *osg_295 = buffer.data(osg + 295);
    const auto *osg_297 = buffer.data(osg + 297);
    const auto *osg_298 = buffer.data(osg + 298);
    const auto *osg_299 = buffer.data(osg + 299);
    const auto *osg_300 = buffer.data(osg + 300);
    const auto *osg_301 = buffer.data(osg + 301);
    const auto *osg_302 = buffer.data(osg + 302);
    const auto *osg_303 = buffer.data(osg + 303);
    const auto *osg_305 = buffer.data(osg + 305);
    const auto *osg_310 = buffer.data(osg + 310);
    const auto *osg_312 = buffer.data(osg + 312);
    const auto *osg_313 = buffer.data(osg + 313);
    const auto *osg_314 = buffer.data(osg + 314);
    const auto *osg_315 = buffer.data(osg + 315);
    const auto *osg_320 = buffer.data(osg + 320);
    const auto *osg_325 = buffer.data(osg + 325);
    const auto *osg_350 = buffer.data(osg + 350);
    const auto *osg_351 = buffer.data(osg + 351);
    const auto *osg_354 = buffer.data(osg + 354);
    const auto *osg_355 = buffer.data(osg + 355);
    const auto *osg_356 = buffer.data(osg + 356);
    const auto *osg_357 = buffer.data(osg + 357);
    const auto *osg_358 = buffer.data(osg + 358);
    const auto *osg_359 = buffer.data(osg + 359);
    const auto *osg_360 = buffer.data(osg + 360);
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
    const auto *osg_378 = buffer.data(osg + 378);
    const auto *osg_380 = buffer.data(osg + 380);
    const auto *osg_381 = buffer.data(osg + 381);
    const auto *osg_384 = buffer.data(osg + 384);
    const auto *osg_385 = buffer.data(osg + 385);
    const auto *osg_386 = buffer.data(osg + 386);
    const auto *osg_387 = buffer.data(osg + 387);
    const auto *osg_388 = buffer.data(osg + 388);
    const auto *osg_389 = buffer.data(osg + 389);
    const auto *osg_400 = buffer.data(osg + 400);
    const auto *osg_401 = buffer.data(osg + 401);
    const auto *osg_402 = buffer.data(osg + 402);
    const auto *osg_403 = buffer.data(osg + 403);
    const auto *osg_404 = buffer.data(osg + 404);
    const auto *osg_405 = buffer.data(osg + 405);
    const auto *osg_410 = buffer.data(osg + 410);
    const auto *osg_414 = buffer.data(osg + 414);
    const auto *osg_415 = buffer.data(osg + 415);
    const auto *osg_416 = buffer.data(osg + 416);
    const auto *osg_417 = buffer.data(osg + 417);
    const auto *osg_419 = buffer.data(osg + 419);
    const auto *osg_420 = buffer.data(osg + 420);
    const auto *osg_423 = buffer.data(osg + 423);
    const auto *osg_426 = buffer.data(osg + 426);
    const auto *osg_430 = buffer.data(osg + 430);
    const auto *osg_432 = buffer.data(osg + 432);
    const auto *osg_433 = buffer.data(osg + 433);
    const auto *osg_434 = buffer.data(osg + 434);

    const auto *osh1_420 = buffer.data(osh1 + 420);
    const auto *osh1_423 = buffer.data(osh1 + 423);
    const auto *osh1_425 = buffer.data(osh1 + 425);
    const auto *osh1_426 = buffer.data(osh1 + 426);
    const auto *osh1_429 = buffer.data(osh1 + 429);
    const auto *osh1_440 = buffer.data(osh1 + 440);

    const auto *qsf0_235 = buffer.data(qsf0 + 235);
    const auto *qsf0_236 = buffer.data(qsf0 + 236);
    const auto *qsf0_238 = buffer.data(qsf0 + 238);
    const auto *qsf0_239 = buffer.data(qsf0 + 239);
    const auto *qsf0_240 = buffer.data(qsf0 + 240);
    const auto *qsf0_243 = buffer.data(qsf0 + 243);
    const auto *qsf0_245 = buffer.data(qsf0 + 245);
    const auto *qsf0_246 = buffer.data(qsf0 + 246);
    const auto *qsf0_248 = buffer.data(qsf0 + 248);
    const auto *qsf0_249 = buffer.data(qsf0 + 249);
    const auto *qsf0_250 = buffer.data(qsf0 + 250);
    const auto *qsf0_253 = buffer.data(qsf0 + 253);
    const auto *qsf0_255 = buffer.data(qsf0 + 255);
    const auto *qsf0_256 = buffer.data(qsf0 + 256);
    const auto *qsf0_258 = buffer.data(qsf0 + 258);
    const auto *qsf0_259 = buffer.data(qsf0 + 259);
    const auto *qsf0_266 = buffer.data(qsf0 + 266);
    const auto *qsf0_268 = buffer.data(qsf0 + 268);
    const auto *qsf0_269 = buffer.data(qsf0 + 269);
    const auto *qsf0_270 = buffer.data(qsf0 + 270);
    const auto *qsf0_271 = buffer.data(qsf0 + 271);
    const auto *qsf0_272 = buffer.data(qsf0 + 272);
    const auto *qsf0_275 = buffer.data(qsf0 + 275);
    const auto *qsf0_276 = buffer.data(qsf0 + 276);
    const auto *qsf0_277 = buffer.data(qsf0 + 277);
    const auto *qsf0_278 = buffer.data(qsf0 + 278);
    const auto *qsf0_279 = buffer.data(qsf0 + 279);
    const auto *qsf0_280 = buffer.data(qsf0 + 280);
    const auto *qsf0_282 = buffer.data(qsf0 + 282);
    const auto *qsf0_283 = buffer.data(qsf0 + 283);
    const auto *qsf0_286 = buffer.data(qsf0 + 286);

    const auto *qsf1_235 = buffer.data(qsf1 + 235);
    const auto *qsf1_236 = buffer.data(qsf1 + 236);
    const auto *qsf1_238 = buffer.data(qsf1 + 238);
    const auto *qsf1_239 = buffer.data(qsf1 + 239);
    const auto *qsf1_240 = buffer.data(qsf1 + 240);
    const auto *qsf1_243 = buffer.data(qsf1 + 243);
    const auto *qsf1_245 = buffer.data(qsf1 + 245);
    const auto *qsf1_246 = buffer.data(qsf1 + 246);
    const auto *qsf1_248 = buffer.data(qsf1 + 248);
    const auto *qsf1_249 = buffer.data(qsf1 + 249);
    const auto *qsf1_250 = buffer.data(qsf1 + 250);
    const auto *qsf1_253 = buffer.data(qsf1 + 253);
    const auto *qsf1_255 = buffer.data(qsf1 + 255);
    const auto *qsf1_256 = buffer.data(qsf1 + 256);
    const auto *qsf1_258 = buffer.data(qsf1 + 258);
    const auto *qsf1_259 = buffer.data(qsf1 + 259);
    const auto *qsf1_266 = buffer.data(qsf1 + 266);
    const auto *qsf1_268 = buffer.data(qsf1 + 268);
    const auto *qsf1_269 = buffer.data(qsf1 + 269);
    const auto *qsf1_270 = buffer.data(qsf1 + 270);
    const auto *qsf1_271 = buffer.data(qsf1 + 271);
    const auto *qsf1_272 = buffer.data(qsf1 + 272);
    const auto *qsf1_275 = buffer.data(qsf1 + 275);
    const auto *qsf1_276 = buffer.data(qsf1 + 276);
    const auto *qsf1_277 = buffer.data(qsf1 + 277);
    const auto *qsf1_278 = buffer.data(qsf1 + 278);
    const auto *qsf1_279 = buffer.data(qsf1 + 279);
    const auto *qsf1_280 = buffer.data(qsf1 + 280);
    const auto *qsf1_282 = buffer.data(qsf1 + 282);
    const auto *qsf1_283 = buffer.data(qsf1 + 283);
    const auto *qsf1_286 = buffer.data(qsf1 + 286);

    const auto *qsg_348 = buffer.data(qsg + 348);
    const auto *qsg_350 = buffer.data(qsg + 350);
    const auto *qsg_351 = buffer.data(qsg + 351);
    const auto *qsg_354 = buffer.data(qsg + 354);
    const auto *qsg_355 = buffer.data(qsg + 355);
    const auto *qsg_356 = buffer.data(qsg + 356);
    const auto *qsg_357 = buffer.data(qsg + 357);
    const auto *qsg_358 = buffer.data(qsg + 358);
    const auto *qsg_359 = buffer.data(qsg + 359);
    const auto *qsg_360 = buffer.data(qsg + 360);
    const auto *qsg_362 = buffer.data(qsg + 362);
    const auto *qsg_363 = buffer.data(qsg + 363);
    const auto *qsg_365 = buffer.data(qsg + 365);
    const auto *qsg_366 = buffer.data(qsg + 366);
    const auto *qsg_369 = buffer.data(qsg + 369);
    const auto *qsg_370 = buffer.data(qsg + 370);
    const auto *qsg_371 = buffer.data(qsg + 371);
    const auto *qsg_372 = buffer.data(qsg + 372);
    const auto *qsg_373 = buffer.data(qsg + 373);
    const auto *qsg_374 = buffer.data(qsg + 374);
    const auto *qsg_375 = buffer.data(qsg + 375);
    const auto *qsg_377 = buffer.data(qsg + 377);
    const auto *qsg_378 = buffer.data(qsg + 378);
    const auto *qsg_380 = buffer.data(qsg + 380);
    const auto *qsg_381 = buffer.data(qsg + 381);
    const auto *qsg_384 = buffer.data(qsg + 384);
    const auto *qsg_385 = buffer.data(qsg + 385);
    const auto *qsg_386 = buffer.data(qsg + 386);
    const auto *qsg_387 = buffer.data(qsg + 387);
    const auto *qsg_388 = buffer.data(qsg + 388);
    const auto *qsg_389 = buffer.data(qsg + 389);
    const auto *qsg_390 = buffer.data(qsg + 390);
    const auto *qsg_392 = buffer.data(qsg + 392);
    const auto *qsg_393 = buffer.data(qsg + 393);
    const auto *qsg_395 = buffer.data(qsg + 395);
    const auto *qsg_400 = buffer.data(qsg + 400);
    const auto *qsg_401 = buffer.data(qsg + 401);
    const auto *qsg_402 = buffer.data(qsg + 402);
    const auto *qsg_403 = buffer.data(qsg + 403);
    const auto *qsg_404 = buffer.data(qsg + 404);
    const auto *qsg_405 = buffer.data(qsg + 405);
    const auto *qsg_406 = buffer.data(qsg + 406);
    const auto *qsg_407 = buffer.data(qsg + 407);
    const auto *qsg_408 = buffer.data(qsg + 408);
    const auto *qsg_409 = buffer.data(qsg + 409);
    const auto *qsg_410 = buffer.data(qsg + 410);
    const auto *qsg_414 = buffer.data(qsg + 414);
    const auto *qsg_415 = buffer.data(qsg + 415);
    const auto *qsg_416 = buffer.data(qsg + 416);
    const auto *qsg_417 = buffer.data(qsg + 417);
    const auto *qsg_418 = buffer.data(qsg + 418);
    const auto *qsg_419 = buffer.data(qsg + 419);
    const auto *qsg_420 = buffer.data(qsg + 420);
    const auto *qsg_421 = buffer.data(qsg + 421);
    const auto *qsg_422 = buffer.data(qsg + 422);
    const auto *qsg_423 = buffer.data(qsg + 423);
    const auto *qsg_425 = buffer.data(qsg + 425);
    const auto *qsg_426 = buffer.data(qsg + 426);
    const auto *qsg_430 = buffer.data(qsg + 430);
    const auto *qsg_432 = buffer.data(qsg + 432);
    const auto *qsg_433 = buffer.data(qsg + 433);
    const auto *qsg_434 = buffer.data(qsg + 434);

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, osg_243, osg_350, osg_351, qsf0_235, \
                         qsf0_236, qsf1_235, qsf1_236, qsg_348, qsg_350, \
                         qsg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_21 * osg_350[k]
                   + f_6 * qsf0_235[k]
                   - f_7 * qsf1_235[k]
                   + f_3 * pc_x[k] * qsg_350[k];

        t_489[k] = f_21 * osg_351[k]
                   + f_4 * qsf0_236[k]
                   - f_5 * qsf1_236[k]
                   + f_3 * pc_x[k] * qsg_351[k];

        t_490[k] = f_10 * osg_243[k]
                   + f_3 * pc_z[k] * qsg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, osg_260, osg_354, osg_355, \
                         osg_356, qsf0_239, qsf1_239, qsg_350, qsg_354, qsg_355, \
                         qsg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_18 * osg_260[k]
                   + f_3 * pc_y[k] * qsg_350[k];

        t_492[k] = f_21 * osg_354[k]
                   + f_4 * qsf0_239[k]
                   - f_5 * qsf1_239[k]
                   + f_3 * pc_x[k] * qsg_354[k];

        t_493[k] = f_21 * osg_355[k]
                   + f_3 * pc_x[k] * qsg_355[k];

        t_494[k] = f_21 * osg_356[k]
                   + f_3 * pc_x[k] * qsg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, osg_265, osg_357, osg_358, \
                         osg_359, qsf0_236, qsf1_236, qsg_355, qsg_357, qsg_358, \
                         qsg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_21 * osg_357[k]
                   + f_3 * pc_x[k] * qsg_357[k];

        t_496[k] = f_21 * osg_358[k]
                   + f_3 * pc_x[k] * qsg_358[k];

        t_497[k] = f_21 * osg_359[k]
                   + f_3 * pc_x[k] * qsg_359[k];

        t_498[k] = f_18 * osg_265[k]
                   + f_1 * qsf0_236[k]
                   - f_2 * qsf1_236[k]
                   + f_3 * pc_y[k] * qsg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, osg_250, osg_267, osg_268, qsf0_238, \
                         qsf0_239, qsf1_238, qsf1_239, qsg_355, qsg_357, \
                         qsg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * osg_250[k]
                   + f_3 * pc_z[k] * qsg_355[k];

        t_500[k] = f_18 * osg_267[k]
                   + f_6 * qsf0_238[k]
                   - f_7 * qsf1_238[k]
                   + f_3 * pc_y[k] * qsg_357[k];

        t_501[k] = f_18 * osg_268[k]
                   + f_4 * qsf0_239[k]
                   - f_5 * qsf1_239[k]
                   + f_3 * pc_y[k] * qsg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, osg_254, osg_269, osg_360, \
                         qsf0_239, qsf0_240, qsf1_239, qsf1_240, qsg_359, \
                         qsg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_18 * osg_269[k]
                   + f_3 * pc_y[k] * qsg_359[k];

        t_503[k] = f_10 * osg_254[k]
                   + f_1 * qsf0_239[k]
                   - f_2 * qsf1_239[k]
                   + f_3 * pc_z[k] * qsg_359[k];

        t_504[k] = f_21 * osg_360[k]
                   + f_1 * qsf0_240[k]
                   - f_2 * qsf1_240[k]
                   + f_3 * pc_x[k] * qsg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, osg_255, osg_270, \
                         osg_272, osg_363, qsf0_243, qsf1_243, qsg_360, qsg_362, \
                         qsg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * osg_270[k]
                   + f_3 * pc_y[k] * qsg_360[k];

        t_506[k] = f_11 * osg_255[k]
                   + f_3 * pc_z[k] * qsg_360[k];

        t_507[k] = f_21 * osg_363[k]
                   + f_6 * qsf0_243[k]
                   - f_7 * qsf1_243[k]
                   + f_3 * pc_x[k] * qsg_363[k];

        t_508[k] = f_11 * osg_272[k]
                   + f_3 * pc_y[k] * qsg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, osg_258, osg_365, osg_366, qsf0_245, \
                         qsf0_246, qsf1_245, qsf1_246, qsg_363, qsg_365, \
                         qsg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_21 * osg_365[k]
                   + f_6 * qsf0_245[k]
                   - f_7 * qsf1_245[k]
                   + f_3 * pc_x[k] * qsg_365[k];

        t_510[k] = f_21 * osg_366[k]
                   + f_4 * qsf0_246[k]
                   - f_5 * qsf1_246[k]
                   + f_3 * pc_x[k] * qsg_366[k];

        t_511[k] = f_11 * osg_258[k]
                   + f_3 * pc_z[k] * qsg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, osg_275, osg_369, osg_370, \
                         osg_371, qsf0_249, qsf1_249, qsg_365, qsg_369, qsg_370, \
                         qsg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * osg_275[k]
                   + f_3 * pc_y[k] * qsg_365[k];

        t_513[k] = f_21 * osg_369[k]
                   + f_4 * qsf0_249[k]
                   - f_5 * qsf1_249[k]
                   + f_3 * pc_x[k] * qsg_369[k];

        t_514[k] = f_21 * osg_370[k]
                   + f_3 * pc_x[k] * qsg_370[k];

        t_515[k] = f_21 * osg_371[k]
                   + f_3 * pc_x[k] * qsg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, osg_280, osg_372, osg_373, \
                         osg_374, qsf0_246, qsf1_246, qsg_370, qsg_372, qsg_373, \
                         qsg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_21 * osg_372[k]
                   + f_3 * pc_x[k] * qsg_372[k];

        t_517[k] = f_21 * osg_373[k]
                   + f_3 * pc_x[k] * qsg_373[k];

        t_518[k] = f_21 * osg_374[k]
                   + f_3 * pc_x[k] * qsg_374[k];

        t_519[k] = f_11 * osg_280[k]
                   + f_1 * qsf0_246[k]
                   - f_2 * qsf1_246[k]
                   + f_3 * pc_y[k] * qsg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, osg_265, osg_282, osg_283, qsf0_248, \
                         qsf0_249, qsf1_248, qsf1_249, qsg_370, qsg_372, \
                         qsg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * osg_265[k]
                   + f_3 * pc_z[k] * qsg_370[k];

        t_521[k] = f_11 * osg_282[k]
                   + f_6 * qsf0_248[k]
                   - f_7 * qsf1_248[k]
                   + f_3 * pc_y[k] * qsg_372[k];

        t_522[k] = f_11 * osg_283[k]
                   + f_4 * qsf0_249[k]
                   - f_5 * qsf1_249[k]
                   + f_3 * pc_y[k] * qsg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, osg_269, osg_284, osg_375, \
                         qsf0_249, qsf0_250, qsf1_249, qsf1_250, qsg_374, \
                         qsg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * osg_284[k]
                   + f_3 * pc_y[k] * qsg_374[k];

        t_524[k] = f_11 * osg_269[k]
                   + f_1 * qsf0_249[k]
                   - f_2 * qsf1_249[k]
                   + f_3 * pc_z[k] * qsg_374[k];

        t_525[k] = f_21 * osg_375[k]
                   + f_1 * qsf0_250[k]
                   - f_2 * qsf1_250[k]
                   + f_3 * pc_x[k] * qsg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, osg_270, osg_285, \
                         osg_287, osg_378, qsf0_253, qsf1_253, qsg_375, qsg_377, \
                         qsg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * osg_285[k]
                   + f_3 * pc_y[k] * qsg_375[k];

        t_527[k] = f_18 * osg_270[k]
                   + f_3 * pc_z[k] * qsg_375[k];

        t_528[k] = f_21 * osg_378[k]
                   + f_6 * qsf0_253[k]
                   - f_7 * qsf1_253[k]
                   + f_3 * pc_x[k] * qsg_378[k];

        t_529[k] = f_10 * osg_287[k]
                   + f_3 * pc_y[k] * qsg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, osg_273, osg_380, osg_381, qsf0_255, \
                         qsf0_256, qsf1_255, qsf1_256, qsg_378, qsg_380, \
                         qsg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_21 * osg_380[k]
                   + f_6 * qsf0_255[k]
                   - f_7 * qsf1_255[k]
                   + f_3 * pc_x[k] * qsg_380[k];

        t_531[k] = f_21 * osg_381[k]
                   + f_4 * qsf0_256[k]
                   - f_5 * qsf1_256[k]
                   + f_3 * pc_x[k] * qsg_381[k];

        t_532[k] = f_18 * osg_273[k]
                   + f_3 * pc_z[k] * qsg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, osg_290, osg_384, osg_385, \
                         osg_386, qsf0_259, qsf1_259, qsg_380, qsg_384, qsg_385, \
                         qsg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * osg_290[k]
                   + f_3 * pc_y[k] * qsg_380[k];

        t_534[k] = f_21 * osg_384[k]
                   + f_4 * qsf0_259[k]
                   - f_5 * qsf1_259[k]
                   + f_3 * pc_x[k] * qsg_384[k];

        t_535[k] = f_21 * osg_385[k]
                   + f_3 * pc_x[k] * qsg_385[k];

        t_536[k] = f_21 * osg_386[k]
                   + f_3 * pc_x[k] * qsg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, osg_295, osg_387, osg_388, \
                         osg_389, qsf0_256, qsf1_256, qsg_385, qsg_387, qsg_388, \
                         qsg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_21 * osg_387[k]
                   + f_3 * pc_x[k] * qsg_387[k];

        t_538[k] = f_21 * osg_388[k]
                   + f_3 * pc_x[k] * qsg_388[k];

        t_539[k] = f_21 * osg_389[k]
                   + f_3 * pc_x[k] * qsg_389[k];

        t_540[k] = f_10 * osg_295[k]
                   + f_1 * qsf0_256[k]
                   - f_2 * qsf1_256[k]
                   + f_3 * pc_y[k] * qsg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, osg_280, osg_297, osg_298, qsf0_258, \
                         qsf0_259, qsf1_258, qsf1_259, qsg_385, qsg_387, \
                         qsg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_18 * osg_280[k]
                   + f_3 * pc_z[k] * qsg_385[k];

        t_542[k] = f_10 * osg_297[k]
                   + f_6 * qsf0_258[k]
                   - f_7 * qsf1_258[k]
                   + f_3 * pc_y[k] * qsg_387[k];

        t_543[k] = f_10 * osg_298[k]
                   + f_4 * qsf0_259[k]
                   - f_5 * qsf1_259[k]
                   + f_3 * pc_y[k] * qsg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_y, pc_z, osh0_420, osg_284, \
                         osg_299, osg_300, osh1_420, qsf0_259, qsf1_259, qsg_389, \
                         qsg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * osg_299[k]
                   + f_3 * pc_y[k] * qsg_389[k];

        t_545[k] = f_18 * osg_284[k]
                   + f_1 * qsf0_259[k]
                   - f_2 * qsf1_259[k]
                   + f_3 * pc_z[k] * qsg_389[k];

        t_546[k] = pa_y[k] * osh0_420[k]
                   - f_8 * pc_y[k] * osh1_420[k];

        t_547[k] = f_9 * osg_300[k]
                   + f_3 * pc_y[k] * qsg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pa_y, pc_y, pc_z, osh0_423, osh0_425, \
                         osg_285, osg_301, osg_302, osh1_423, osh1_425, qsg_390, \
                         qsg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_20 * osg_285[k]
                   + f_3 * pc_z[k] * qsg_390[k];

        t_549[k] = pa_y[k] * osh0_423[k]
                   + f_10 * osg_301[k]
                   - f_8 * pc_y[k] * osh1_423[k];

        t_550[k] = f_9 * osg_302[k]
                   + f_3 * pc_y[k] * qsg_392[k];

        t_551[k] = pa_y[k] * osh0_425[k]
                   - f_8 * pc_y[k] * osh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_y, pc_y, pc_z, osh0_426, osh0_429, \
                         osg_288, osg_303, osg_305, osh1_426, osh1_429, qsg_393, \
                         qsg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_y[k] * osh0_426[k]
                   + f_11 * osg_303[k]
                   - f_8 * pc_y[k] * osh1_426[k];

        t_553[k] = f_20 * osg_288[k]
                   + f_3 * pc_z[k] * qsg_393[k];

        t_554[k] = f_9 * osg_305[k]
                   + f_3 * pc_y[k] * qsg_395[k];

        t_555[k] = pa_y[k] * osh0_429[k]
                   - f_8 * pc_y[k] * osh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, osg_400, osg_401, osg_402, \
                         osg_403, osg_404, qsg_400, qsg_401, qsg_402, qsg_403, \
                         qsg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_21 * osg_400[k]
                   + f_3 * pc_x[k] * qsg_400[k];

        t_557[k] = f_21 * osg_401[k]
                   + f_3 * pc_x[k] * qsg_401[k];

        t_558[k] = f_21 * osg_402[k]
                   + f_3 * pc_x[k] * qsg_402[k];

        t_559[k] = f_21 * osg_403[k]
                   + f_3 * pc_x[k] * qsg_403[k];

        t_560[k] = f_21 * osg_404[k]
                   + f_3 * pc_x[k] * qsg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, osg_295, osg_310, osg_312, qsf0_266, \
                         qsf0_268, qsf1_266, qsf1_268, qsg_400, \
                         qsg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * osg_310[k]
                   + f_1 * qsf0_266[k]
                   - f_2 * qsf1_266[k]
                   + f_3 * pc_y[k] * qsg_400[k];

        t_562[k] = f_20 * osg_295[k]
                   + f_3 * pc_z[k] * qsg_400[k];

        t_563[k] = f_9 * osg_312[k]
                   + f_6 * qsf0_268[k]
                   - f_7 * qsf1_268[k]
                   + f_3 * pc_y[k] * qsg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_y, pc_y, osh0_440, osg_313, osg_314, \
                         osh1_440, qsf0_269, qsf1_269, qsg_403, \
                         qsg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * osg_313[k]
                   + f_4 * qsf0_269[k]
                   - f_5 * qsf1_269[k]
                   + f_3 * pc_y[k] * qsg_403[k];

        t_565[k] = f_9 * osg_314[k]
                   + f_3 * pc_y[k] * qsg_404[k];

        t_566[k] = pa_y[k] * osh0_440[k]
                   - f_8 * pc_y[k] * osh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, osg_300, \
                         osg_405, qsf0_270, qsf1_270, qsg_405, qsg_406, \
                         qsg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_21 * osg_405[k]
                   + f_1 * qsf0_270[k]
                   - f_2 * qsf1_270[k]
                   + f_3 * pc_x[k] * qsg_405[k];

        t_568[k] = f_3 * pc_y[k] * qsg_405[k];

        t_569[k] = f_21 * osg_300[k]
                   + f_3 * pc_z[k] * qsg_405[k];

        t_570[k] = f_4 * qsf0_270[k]
                   - f_5 * qsf1_270[k]
                   + f_3 * pc_y[k] * qsg_406[k];

        t_571[k] = f_3 * pc_y[k] * qsg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, osg_410, qsf0_271, qsf0_272, \
                         qsf0_275, qsf1_271, qsf1_272, qsf1_275, qsg_408, qsg_409, \
                         qsg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_21 * osg_410[k]
                   + f_6 * qsf0_275[k]
                   - f_7 * qsf1_275[k]
                   + f_3 * pc_x[k] * qsg_410[k];

        t_573[k] = f_6 * qsf0_271[k]
                   - f_7 * qsf1_271[k]
                   + f_3 * pc_y[k] * qsg_408[k];

        t_574[k] = f_4 * qsf0_272[k]
                   - f_5 * qsf1_272[k]
                   + f_3 * pc_y[k] * qsg_409[k];

        t_575[k] = f_3 * pc_y[k] * qsg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, osg_414, osg_415, osg_416, osg_417, \
                         qsf0_279, qsf1_279, qsg_414, qsg_415, qsg_416, \
                         qsg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_21 * osg_414[k]
                   + f_4 * qsf0_279[k]
                   - f_5 * qsf1_279[k]
                   + f_3 * pc_x[k] * qsg_414[k];

        t_577[k] = f_21 * osg_415[k]
                   + f_3 * pc_x[k] * qsg_415[k];

        t_578[k] = f_21 * osg_416[k]
                   + f_3 * pc_x[k] * qsg_416[k];

        t_579[k] = f_21 * osg_417[k]
                   + f_3 * pc_x[k] * qsg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, osg_419, qsf0_276, qsf0_277, \
                         qsf1_276, qsf1_277, qsg_414, qsg_415, qsg_416, \
                         qsg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_3 * pc_y[k] * qsg_414[k];

        t_581[k] = f_21 * osg_419[k]
                   + f_3 * pc_x[k] * qsg_419[k];

        t_582[k] = f_1 * qsf0_276[k]
                   - f_2 * qsf1_276[k]
                   + f_3 * pc_y[k] * qsg_415[k];

        t_583[k] = f_13 * qsf0_277[k]
                   - f_14 * qsf1_277[k]
                   + f_3 * pc_y[k] * qsg_416[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, osg_314, qsf0_278, qsf0_279, \
                         qsf1_278, qsf1_279, qsg_417, qsg_418, \
                         qsg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * qsf0_278[k]
                   - f_7 * qsf1_278[k]
                   + f_3 * pc_y[k] * qsg_417[k];

        t_585[k] = f_4 * qsf0_279[k]
                   - f_5 * qsf1_279[k]
                   + f_3 * pc_y[k] * qsg_418[k];

        t_586[k] = f_3 * pc_y[k] * qsg_419[k];

        t_587[k] = f_21 * osg_314[k]
                   + f_1 * qsf0_279[k]
                   - f_2 * qsf1_279[k]
                   + f_3 * pc_z[k] * qsg_419[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, osg_315, osg_420, \
                         osg_423, qsf0_280, qsf0_283, qsf1_280, qsf1_283, qsg_420, \
                         qsg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_20 * osg_420[k]
                   + f_1 * qsf0_280[k]
                   - f_2 * qsf1_280[k]
                   + f_3 * pc_x[k] * qsg_420[k];

        t_589[k] = f_19 * osg_315[k]
                   + f_3 * pc_y[k] * qsg_420[k];

        t_590[k] = f_3 * pc_z[k] * qsg_420[k];

        t_591[k] = f_20 * osg_423[k]
                   + f_6 * qsf0_283[k]
                   - f_7 * qsf1_283[k]
                   + f_3 * pc_x[k] * qsg_423[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, osg_426, qsf0_280, qsf0_286, \
                         qsf1_280, qsf1_286, qsg_421, qsg_422, qsg_423, \
                         qsg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * qsg_421[k];

        t_593[k] = f_4 * qsf0_280[k]
                   - f_5 * qsf1_280[k]
                   + f_3 * pc_z[k] * qsg_422[k];

        t_594[k] = f_20 * osg_426[k]
                   + f_4 * qsf0_286[k]
                   - f_5 * qsf1_286[k]
                   + f_3 * pc_x[k] * qsg_426[k];

        t_595[k] = f_3 * pc_z[k] * qsg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, osg_320, osg_430, \
                         qsf0_282, qsf1_282, qsg_425, qsg_426, \
                         qsg_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_19 * osg_320[k]
                   + f_3 * pc_y[k] * qsg_425[k];

        t_597[k] = f_6 * qsf0_282[k]
                   - f_7 * qsf1_282[k]
                   + f_3 * pc_z[k] * qsg_425[k];

        t_598[k] = f_20 * osg_430[k]
                   + f_3 * pc_x[k] * qsg_430[k];

        t_599[k] = f_3 * pc_z[k] * qsg_426[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, osg_325, osg_432, osg_433, \
                         osg_434, qsf0_286, qsf1_286, qsg_430, qsg_432, qsg_433, \
                         qsg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_20 * osg_432[k]
                   + f_3 * pc_x[k] * qsg_432[k];

        t_601[k] = f_20 * osg_433[k]
                   + f_3 * pc_x[k] * qsg_433[k];

        t_602[k] = f_20 * osg_434[k]
                   + f_3 * pc_x[k] * qsg_434[k];

        t_603[k] = f_19 * osg_325[k]
                   + f_1 * qsf0_286[k]
                   - f_2 * qsf1_286[k]
                   + f_3 * pc_y[k] * qsg_430[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *osh0_441 = buffer.data(osh0 + 441);
    const auto *osh0_444 = buffer.data(osh0 + 444);
    const auto *osh0_447 = buffer.data(osh0 + 447);
    const auto *osh0_456 = buffer.data(osh0 + 456);

    const auto *osg_315 = buffer.data(osg + 315);
    const auto *osg_318 = buffer.data(osg + 318);
    const auto *osg_325 = buffer.data(osg + 325);
    const auto *osg_329 = buffer.data(osg + 329);
    const auto *osg_330 = buffer.data(osg + 330);
    const auto *osg_332 = buffer.data(osg + 332);
    const auto *osg_333 = buffer.data(osg + 333);
    const auto *osg_335 = buffer.data(osg + 335);
    const auto *osg_340 = buffer.data(osg + 340);
    const auto *osg_342 = buffer.data(osg + 342);
    const auto *osg_343 = buffer.data(osg + 343);
    const auto *osg_344 = buffer.data(osg + 344);
    const auto *osg_345 = buffer.data(osg + 345);
    const auto *osg_347 = buffer.data(osg + 347);
    const auto *osg_348 = buffer.data(osg + 348);
    const auto *osg_350 = buffer.data(osg + 350);
    const auto *osg_355 = buffer.data(osg + 355);
    const auto *osg_357 = buffer.data(osg + 357);
    const auto *osg_358 = buffer.data(osg + 358);
    const auto *osg_359 = buffer.data(osg + 359);
    const auto *osg_360 = buffer.data(osg + 360);
    const auto *osg_362 = buffer.data(osg + 362);
    const auto *osg_363 = buffer.data(osg + 363);
    const auto *osg_365 = buffer.data(osg + 365);
    const auto *osg_370 = buffer.data(osg + 370);
    const auto *osg_372 = buffer.data(osg + 372);
    const auto *osg_373 = buffer.data(osg + 373);
    const auto *osg_374 = buffer.data(osg + 374);
    const auto *osg_375 = buffer.data(osg + 375);
    const auto *osg_377 = buffer.data(osg + 377);
    const auto *osg_378 = buffer.data(osg + 378);
    const auto *osg_380 = buffer.data(osg + 380);
    const auto *osg_385 = buffer.data(osg + 385);
    const auto *osg_387 = buffer.data(osg + 387);
    const auto *osg_388 = buffer.data(osg + 388);
    const auto *osg_389 = buffer.data(osg + 389);
    const auto *osg_390 = buffer.data(osg + 390);
    const auto *osg_392 = buffer.data(osg + 392);
    const auto *osg_395 = buffer.data(osg + 395);
    const auto *osg_400 = buffer.data(osg + 400);
    const auto *osg_402 = buffer.data(osg + 402);
    const auto *osg_403 = buffer.data(osg + 403);
    const auto *osg_440 = buffer.data(osg + 440);
    const auto *osg_444 = buffer.data(osg + 444);
    const auto *osg_445 = buffer.data(osg + 445);
    const auto *osg_446 = buffer.data(osg + 446);
    const auto *osg_447 = buffer.data(osg + 447);
    const auto *osg_448 = buffer.data(osg + 448);
    const auto *osg_449 = buffer.data(osg + 449);
    const auto *osg_450 = buffer.data(osg + 450);
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
    const auto *osg_498 = buffer.data(osg + 498);
    const auto *osg_500 = buffer.data(osg + 500);
    const auto *osg_501 = buffer.data(osg + 501);
    const auto *osg_504 = buffer.data(osg + 504);
    const auto *osg_505 = buffer.data(osg + 505);
    const auto *osg_506 = buffer.data(osg + 506);
    const auto *osg_507 = buffer.data(osg + 507);
    const auto *osg_508 = buffer.data(osg + 508);
    const auto *osg_509 = buffer.data(osg + 509);

    const auto *osh1_441 = buffer.data(osh1 + 441);
    const auto *osh1_444 = buffer.data(osh1 + 444);
    const auto *osh1_447 = buffer.data(osh1 + 447);
    const auto *osh1_456 = buffer.data(osh1 + 456);

    const auto *qsf0_286 = buffer.data(qsf0 + 286);
    const auto *qsf0_287 = buffer.data(qsf0 + 287);
    const auto *qsf0_289 = buffer.data(qsf0 + 289);
    const auto *qsf0_295 = buffer.data(qsf0 + 295);
    const auto *qsf0_298 = buffer.data(qsf0 + 298);
    const auto *qsf0_299 = buffer.data(qsf0 + 299);
    const auto *qsf0_300 = buffer.data(qsf0 + 300);
    const auto *qsf0_303 = buffer.data(qsf0 + 303);
    const auto *qsf0_305 = buffer.data(qsf0 + 305);
    const auto *qsf0_306 = buffer.data(qsf0 + 306);
    const auto *qsf0_308 = buffer.data(qsf0 + 308);
    const auto *qsf0_309 = buffer.data(qsf0 + 309);
    const auto *qsf0_310 = buffer.data(qsf0 + 310);
    const auto *qsf0_313 = buffer.data(qsf0 + 313);
    const auto *qsf0_315 = buffer.data(qsf0 + 315);
    const auto *qsf0_316 = buffer.data(qsf0 + 316);
    const auto *qsf0_318 = buffer.data(qsf0 + 318);
    const auto *qsf0_319 = buffer.data(qsf0 + 319);
    const auto *qsf0_320 = buffer.data(qsf0 + 320);
    const auto *qsf0_323 = buffer.data(qsf0 + 323);
    const auto *qsf0_325 = buffer.data(qsf0 + 325);
    const auto *qsf0_326 = buffer.data(qsf0 + 326);
    const auto *qsf0_328 = buffer.data(qsf0 + 328);
    const auto *qsf0_329 = buffer.data(qsf0 + 329);
    const auto *qsf0_330 = buffer.data(qsf0 + 330);
    const auto *qsf0_333 = buffer.data(qsf0 + 333);
    const auto *qsf0_335 = buffer.data(qsf0 + 335);
    const auto *qsf0_336 = buffer.data(qsf0 + 336);
    const auto *qsf0_338 = buffer.data(qsf0 + 338);
    const auto *qsf0_339 = buffer.data(qsf0 + 339);

    const auto *qsf1_286 = buffer.data(qsf1 + 286);
    const auto *qsf1_287 = buffer.data(qsf1 + 287);
    const auto *qsf1_289 = buffer.data(qsf1 + 289);
    const auto *qsf1_295 = buffer.data(qsf1 + 295);
    const auto *qsf1_298 = buffer.data(qsf1 + 298);
    const auto *qsf1_299 = buffer.data(qsf1 + 299);
    const auto *qsf1_300 = buffer.data(qsf1 + 300);
    const auto *qsf1_303 = buffer.data(qsf1 + 303);
    const auto *qsf1_305 = buffer.data(qsf1 + 305);
    const auto *qsf1_306 = buffer.data(qsf1 + 306);
    const auto *qsf1_308 = buffer.data(qsf1 + 308);
    const auto *qsf1_309 = buffer.data(qsf1 + 309);
    const auto *qsf1_310 = buffer.data(qsf1 + 310);
    const auto *qsf1_313 = buffer.data(qsf1 + 313);
    const auto *qsf1_315 = buffer.data(qsf1 + 315);
    const auto *qsf1_316 = buffer.data(qsf1 + 316);
    const auto *qsf1_318 = buffer.data(qsf1 + 318);
    const auto *qsf1_319 = buffer.data(qsf1 + 319);
    const auto *qsf1_320 = buffer.data(qsf1 + 320);
    const auto *qsf1_323 = buffer.data(qsf1 + 323);
    const auto *qsf1_325 = buffer.data(qsf1 + 325);
    const auto *qsf1_326 = buffer.data(qsf1 + 326);
    const auto *qsf1_328 = buffer.data(qsf1 + 328);
    const auto *qsf1_329 = buffer.data(qsf1 + 329);
    const auto *qsf1_330 = buffer.data(qsf1 + 330);
    const auto *qsf1_333 = buffer.data(qsf1 + 333);
    const auto *qsf1_335 = buffer.data(qsf1 + 335);
    const auto *qsf1_336 = buffer.data(qsf1 + 336);
    const auto *qsf1_338 = buffer.data(qsf1 + 338);
    const auto *qsf1_339 = buffer.data(qsf1 + 339);

    const auto *qsg_430 = buffer.data(qsg + 430);
    const auto *qsg_431 = buffer.data(qsg + 431);
    const auto *qsg_432 = buffer.data(qsg + 432);
    const auto *qsg_434 = buffer.data(qsg + 434);
    const auto *qsg_435 = buffer.data(qsg + 435);
    const auto *qsg_437 = buffer.data(qsg + 437);
    const auto *qsg_438 = buffer.data(qsg + 438);
    const auto *qsg_440 = buffer.data(qsg + 440);
    const auto *qsg_444 = buffer.data(qsg + 444);
    const auto *qsg_445 = buffer.data(qsg + 445);
    const auto *qsg_446 = buffer.data(qsg + 446);
    const auto *qsg_447 = buffer.data(qsg + 447);
    const auto *qsg_448 = buffer.data(qsg + 448);
    const auto *qsg_449 = buffer.data(qsg + 449);
    const auto *qsg_450 = buffer.data(qsg + 450);
    const auto *qsg_452 = buffer.data(qsg + 452);
    const auto *qsg_453 = buffer.data(qsg + 453);
    const auto *qsg_455 = buffer.data(qsg + 455);
    const auto *qsg_456 = buffer.data(qsg + 456);
    const auto *qsg_459 = buffer.data(qsg + 459);
    const auto *qsg_460 = buffer.data(qsg + 460);
    const auto *qsg_461 = buffer.data(qsg + 461);
    const auto *qsg_462 = buffer.data(qsg + 462);
    const auto *qsg_463 = buffer.data(qsg + 463);
    const auto *qsg_464 = buffer.data(qsg + 464);
    const auto *qsg_465 = buffer.data(qsg + 465);
    const auto *qsg_467 = buffer.data(qsg + 467);
    const auto *qsg_468 = buffer.data(qsg + 468);
    const auto *qsg_470 = buffer.data(qsg + 470);
    const auto *qsg_471 = buffer.data(qsg + 471);
    const auto *qsg_474 = buffer.data(qsg + 474);
    const auto *qsg_475 = buffer.data(qsg + 475);
    const auto *qsg_476 = buffer.data(qsg + 476);
    const auto *qsg_477 = buffer.data(qsg + 477);
    const auto *qsg_478 = buffer.data(qsg + 478);
    const auto *qsg_479 = buffer.data(qsg + 479);
    const auto *qsg_480 = buffer.data(qsg + 480);
    const auto *qsg_482 = buffer.data(qsg + 482);
    const auto *qsg_483 = buffer.data(qsg + 483);
    const auto *qsg_485 = buffer.data(qsg + 485);
    const auto *qsg_486 = buffer.data(qsg + 486);
    const auto *qsg_489 = buffer.data(qsg + 489);
    const auto *qsg_490 = buffer.data(qsg + 490);
    const auto *qsg_491 = buffer.data(qsg + 491);
    const auto *qsg_492 = buffer.data(qsg + 492);
    const auto *qsg_493 = buffer.data(qsg + 493);
    const auto *qsg_494 = buffer.data(qsg + 494);
    const auto *qsg_495 = buffer.data(qsg + 495);
    const auto *qsg_497 = buffer.data(qsg + 497);
    const auto *qsg_498 = buffer.data(qsg + 498);
    const auto *qsg_500 = buffer.data(qsg + 500);
    const auto *qsg_501 = buffer.data(qsg + 501);
    const auto *qsg_504 = buffer.data(qsg + 504);
    const auto *qsg_505 = buffer.data(qsg + 505);
    const auto *qsg_506 = buffer.data(qsg + 506);
    const auto *qsg_507 = buffer.data(qsg + 507);
    const auto *qsg_508 = buffer.data(qsg + 508);
    const auto *qsg_509 = buffer.data(qsg + 509);

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pc_y, pc_z, osg_329, qsf0_286, qsf0_287, \
                         qsf1_286, qsf1_287, qsg_430, qsg_431, qsg_432, \
                         qsg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * qsg_430[k];

        t_605[k] = f_4 * qsf0_286[k]
                   - f_5 * qsf1_286[k]
                   + f_3 * pc_z[k] * qsg_431[k];

        t_606[k] = f_6 * qsf0_287[k]
                   - f_7 * qsf1_287[k]
                   + f_3 * pc_z[k] * qsg_432[k];

        t_607[k] = f_19 * osg_329[k]
                   + f_3 * pc_y[k] * qsg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_z, pc_y, pc_z, osh0_441, osg_315, \
                         osg_330, osh1_441, qsf0_289, qsf1_289, qsg_434, \
                         qsg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_1 * qsf0_289[k]
                   - f_2 * qsf1_289[k]
                   + f_3 * pc_z[k] * qsg_434[k];

        t_609[k] = pa_z[k] * osh0_441[k]
                   - f_8 * pc_z[k] * osh1_441[k];

        t_610[k] = f_21 * osg_330[k]
                   + f_3 * pc_y[k] * qsg_435[k];

        t_611[k] = f_9 * osg_315[k]
                   + f_3 * pc_z[k] * qsg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_z, pc_x, pc_y, pc_z, osh0_444, osg_332, \
                         osg_440, osh1_444, qsf0_295, qsf1_295, qsg_437, \
                         qsg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_z[k] * osh0_444[k]
                   - f_8 * pc_z[k] * osh1_444[k];

        t_613[k] = f_21 * osg_332[k]
                   + f_3 * pc_y[k] * qsg_437[k];

        t_614[k] = f_20 * osg_440[k]
                   + f_6 * qsf0_295[k]
                   - f_7 * qsf1_295[k]
                   + f_3 * pc_x[k] * qsg_440[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pa_z, pc_y, pc_z, osh0_447, osg_318, osg_335, \
                         osh1_447, qsg_438, qsg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pa_z[k] * osh0_447[k]
                   - f_8 * pc_z[k] * osh1_447[k];

        t_616[k] = f_9 * osg_318[k]
                   + f_3 * pc_z[k] * qsg_438[k];

        t_617[k] = f_21 * osg_335[k]
                   + f_3 * pc_y[k] * qsg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pc_x, osg_444, osg_445, osg_446, osg_447, \
                         qsf0_299, qsf1_299, qsg_444, qsg_445, qsg_446, \
                         qsg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_20 * osg_444[k]
                   + f_4 * qsf0_299[k]
                   - f_5 * qsf1_299[k]
                   + f_3 * pc_x[k] * qsg_444[k];

        t_619[k] = f_20 * osg_445[k]
                   + f_3 * pc_x[k] * qsg_445[k];

        t_620[k] = f_20 * osg_446[k]
                   + f_3 * pc_x[k] * qsg_446[k];

        t_621[k] = f_20 * osg_447[k]
                   + f_3 * pc_x[k] * qsg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_z, pc_x, pc_z, osh0_456, osg_325, \
                         osg_448, osg_449, osh1_456, qsg_445, qsg_448, \
                         qsg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_20 * osg_448[k]
                   + f_3 * pc_x[k] * qsg_448[k];

        t_623[k] = f_20 * osg_449[k]
                   + f_3 * pc_x[k] * qsg_449[k];

        t_624[k] = pa_z[k] * osh0_456[k]
                   - f_8 * pc_z[k] * osh1_456[k];

        t_625[k] = f_9 * osg_325[k]
                   + f_3 * pc_z[k] * qsg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, osg_342, osg_343, osg_344, qsf0_298, \
                         qsf0_299, qsf1_298, qsf1_299, qsg_447, qsg_448, \
                         qsg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_21 * osg_342[k]
                   + f_6 * qsf0_298[k]
                   - f_7 * qsf1_298[k]
                   + f_3 * pc_y[k] * qsg_447[k];

        t_627[k] = f_21 * osg_343[k]
                   + f_4 * qsf0_299[k]
                   - f_5 * qsf1_299[k]
                   + f_3 * pc_y[k] * qsg_448[k];

        t_628[k] = f_21 * osg_344[k]
                   + f_3 * pc_y[k] * qsg_449[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, pc_z, osg_329, osg_345, osg_450, \
                         qsf0_299, qsf0_300, qsf1_299, qsf1_300, qsg_449, \
                         qsg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_9 * osg_329[k]
                   + f_1 * qsf0_299[k]
                   - f_2 * qsf1_299[k]
                   + f_3 * pc_z[k] * qsg_449[k];

        t_630[k] = f_20 * osg_450[k]
                   + f_1 * qsf0_300[k]
                   - f_2 * qsf1_300[k]
                   + f_3 * pc_x[k] * qsg_450[k];

        t_631[k] = f_20 * osg_345[k]
                   + f_3 * pc_y[k] * qsg_450[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pc_x, pc_y, pc_z, osg_330, osg_347, osg_453, \
                         qsf0_303, qsf1_303, qsg_450, qsg_452, \
                         qsg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_10 * osg_330[k]
                   + f_3 * pc_z[k] * qsg_450[k];

        t_633[k] = f_20 * osg_453[k]
                   + f_6 * qsf0_303[k]
                   - f_7 * qsf1_303[k]
                   + f_3 * pc_x[k] * qsg_453[k];

        t_634[k] = f_20 * osg_347[k]
                   + f_3 * pc_y[k] * qsg_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, osg_333, osg_455, osg_456, qsf0_305, \
                         qsf0_306, qsf1_305, qsf1_306, qsg_453, qsg_455, \
                         qsg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_20 * osg_455[k]
                   + f_6 * qsf0_305[k]
                   - f_7 * qsf1_305[k]
                   + f_3 * pc_x[k] * qsg_455[k];

        t_636[k] = f_20 * osg_456[k]
                   + f_4 * qsf0_306[k]
                   - f_5 * qsf1_306[k]
                   + f_3 * pc_x[k] * qsg_456[k];

        t_637[k] = f_10 * osg_333[k]
                   + f_3 * pc_z[k] * qsg_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, osg_350, osg_459, osg_460, \
                         osg_461, qsf0_309, qsf1_309, qsg_455, qsg_459, qsg_460, \
                         qsg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_20 * osg_350[k]
                   + f_3 * pc_y[k] * qsg_455[k];

        t_639[k] = f_20 * osg_459[k]
                   + f_4 * qsf0_309[k]
                   - f_5 * qsf1_309[k]
                   + f_3 * pc_x[k] * qsg_459[k];

        t_640[k] = f_20 * osg_460[k]
                   + f_3 * pc_x[k] * qsg_460[k];

        t_641[k] = f_20 * osg_461[k]
                   + f_3 * pc_x[k] * qsg_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, osg_355, osg_462, osg_463, \
                         osg_464, qsf0_306, qsf1_306, qsg_460, qsg_462, qsg_463, \
                         qsg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_20 * osg_462[k]
                   + f_3 * pc_x[k] * qsg_462[k];

        t_643[k] = f_20 * osg_463[k]
                   + f_3 * pc_x[k] * qsg_463[k];

        t_644[k] = f_20 * osg_464[k]
                   + f_3 * pc_x[k] * qsg_464[k];

        t_645[k] = f_20 * osg_355[k]
                   + f_1 * qsf0_306[k]
                   - f_2 * qsf1_306[k]
                   + f_3 * pc_y[k] * qsg_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, osg_340, osg_357, osg_358, qsf0_308, \
                         qsf0_309, qsf1_308, qsf1_309, qsg_460, qsg_462, \
                         qsg_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * osg_340[k]
                   + f_3 * pc_z[k] * qsg_460[k];

        t_647[k] = f_20 * osg_357[k]
                   + f_6 * qsf0_308[k]
                   - f_7 * qsf1_308[k]
                   + f_3 * pc_y[k] * qsg_462[k];

        t_648[k] = f_20 * osg_358[k]
                   + f_4 * qsf0_309[k]
                   - f_5 * qsf1_309[k]
                   + f_3 * pc_y[k] * qsg_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, osg_344, osg_359, osg_465, \
                         qsf0_309, qsf0_310, qsf1_309, qsf1_310, qsg_464, \
                         qsg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_20 * osg_359[k]
                   + f_3 * pc_y[k] * qsg_464[k];

        t_650[k] = f_10 * osg_344[k]
                   + f_1 * qsf0_309[k]
                   - f_2 * qsf1_309[k]
                   + f_3 * pc_z[k] * qsg_464[k];

        t_651[k] = f_20 * osg_465[k]
                   + f_1 * qsf0_310[k]
                   - f_2 * qsf1_310[k]
                   + f_3 * pc_x[k] * qsg_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, osg_345, osg_360, \
                         osg_362, osg_468, qsf0_313, qsf1_313, qsg_465, qsg_467, \
                         qsg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_18 * osg_360[k]
                   + f_3 * pc_y[k] * qsg_465[k];

        t_653[k] = f_11 * osg_345[k]
                   + f_3 * pc_z[k] * qsg_465[k];

        t_654[k] = f_20 * osg_468[k]
                   + f_6 * qsf0_313[k]
                   - f_7 * qsf1_313[k]
                   + f_3 * pc_x[k] * qsg_468[k];

        t_655[k] = f_18 * osg_362[k]
                   + f_3 * pc_y[k] * qsg_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, osg_348, osg_470, osg_471, qsf0_315, \
                         qsf0_316, qsf1_315, qsf1_316, qsg_468, qsg_470, \
                         qsg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_20 * osg_470[k]
                   + f_6 * qsf0_315[k]
                   - f_7 * qsf1_315[k]
                   + f_3 * pc_x[k] * qsg_470[k];

        t_657[k] = f_20 * osg_471[k]
                   + f_4 * qsf0_316[k]
                   - f_5 * qsf1_316[k]
                   + f_3 * pc_x[k] * qsg_471[k];

        t_658[k] = f_11 * osg_348[k]
                   + f_3 * pc_z[k] * qsg_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, osg_365, osg_474, osg_475, \
                         osg_476, qsf0_319, qsf1_319, qsg_470, qsg_474, qsg_475, \
                         qsg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_18 * osg_365[k]
                   + f_3 * pc_y[k] * qsg_470[k];

        t_660[k] = f_20 * osg_474[k]
                   + f_4 * qsf0_319[k]
                   - f_5 * qsf1_319[k]
                   + f_3 * pc_x[k] * qsg_474[k];

        t_661[k] = f_20 * osg_475[k]
                   + f_3 * pc_x[k] * qsg_475[k];

        t_662[k] = f_20 * osg_476[k]
                   + f_3 * pc_x[k] * qsg_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, osg_370, osg_477, osg_478, \
                         osg_479, qsf0_316, qsf1_316, qsg_475, qsg_477, qsg_478, \
                         qsg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_20 * osg_477[k]
                   + f_3 * pc_x[k] * qsg_477[k];

        t_664[k] = f_20 * osg_478[k]
                   + f_3 * pc_x[k] * qsg_478[k];

        t_665[k] = f_20 * osg_479[k]
                   + f_3 * pc_x[k] * qsg_479[k];

        t_666[k] = f_18 * osg_370[k]
                   + f_1 * qsf0_316[k]
                   - f_2 * qsf1_316[k]
                   + f_3 * pc_y[k] * qsg_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, osg_355, osg_372, osg_373, qsf0_318, \
                         qsf0_319, qsf1_318, qsf1_319, qsg_475, qsg_477, \
                         qsg_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * osg_355[k]
                   + f_3 * pc_z[k] * qsg_475[k];

        t_668[k] = f_18 * osg_372[k]
                   + f_6 * qsf0_318[k]
                   - f_7 * qsf1_318[k]
                   + f_3 * pc_y[k] * qsg_477[k];

        t_669[k] = f_18 * osg_373[k]
                   + f_4 * qsf0_319[k]
                   - f_5 * qsf1_319[k]
                   + f_3 * pc_y[k] * qsg_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, osg_359, osg_374, osg_480, \
                         qsf0_319, qsf0_320, qsf1_319, qsf1_320, qsg_479, \
                         qsg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_18 * osg_374[k]
                   + f_3 * pc_y[k] * qsg_479[k];

        t_671[k] = f_11 * osg_359[k]
                   + f_1 * qsf0_319[k]
                   - f_2 * qsf1_319[k]
                   + f_3 * pc_z[k] * qsg_479[k];

        t_672[k] = f_20 * osg_480[k]
                   + f_1 * qsf0_320[k]
                   - f_2 * qsf1_320[k]
                   + f_3 * pc_x[k] * qsg_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, osg_360, osg_375, \
                         osg_377, osg_483, qsf0_323, qsf1_323, qsg_480, qsg_482, \
                         qsg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * osg_375[k]
                   + f_3 * pc_y[k] * qsg_480[k];

        t_674[k] = f_18 * osg_360[k]
                   + f_3 * pc_z[k] * qsg_480[k];

        t_675[k] = f_20 * osg_483[k]
                   + f_6 * qsf0_323[k]
                   - f_7 * qsf1_323[k]
                   + f_3 * pc_x[k] * qsg_483[k];

        t_676[k] = f_11 * osg_377[k]
                   + f_3 * pc_y[k] * qsg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, osg_363, osg_485, osg_486, qsf0_325, \
                         qsf0_326, qsf1_325, qsf1_326, qsg_483, qsg_485, \
                         qsg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_20 * osg_485[k]
                   + f_6 * qsf0_325[k]
                   - f_7 * qsf1_325[k]
                   + f_3 * pc_x[k] * qsg_485[k];

        t_678[k] = f_20 * osg_486[k]
                   + f_4 * qsf0_326[k]
                   - f_5 * qsf1_326[k]
                   + f_3 * pc_x[k] * qsg_486[k];

        t_679[k] = f_18 * osg_363[k]
                   + f_3 * pc_z[k] * qsg_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, osg_380, osg_489, osg_490, \
                         osg_491, qsf0_329, qsf1_329, qsg_485, qsg_489, qsg_490, \
                         qsg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * osg_380[k]
                   + f_3 * pc_y[k] * qsg_485[k];

        t_681[k] = f_20 * osg_489[k]
                   + f_4 * qsf0_329[k]
                   - f_5 * qsf1_329[k]
                   + f_3 * pc_x[k] * qsg_489[k];

        t_682[k] = f_20 * osg_490[k]
                   + f_3 * pc_x[k] * qsg_490[k];

        t_683[k] = f_20 * osg_491[k]
                   + f_3 * pc_x[k] * qsg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, osg_385, osg_492, osg_493, \
                         osg_494, qsf0_326, qsf1_326, qsg_490, qsg_492, qsg_493, \
                         qsg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_20 * osg_492[k]
                   + f_3 * pc_x[k] * qsg_492[k];

        t_685[k] = f_20 * osg_493[k]
                   + f_3 * pc_x[k] * qsg_493[k];

        t_686[k] = f_20 * osg_494[k]
                   + f_3 * pc_x[k] * qsg_494[k];

        t_687[k] = f_11 * osg_385[k]
                   + f_1 * qsf0_326[k]
                   - f_2 * qsf1_326[k]
                   + f_3 * pc_y[k] * qsg_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, osg_370, osg_387, osg_388, qsf0_328, \
                         qsf0_329, qsf1_328, qsf1_329, qsg_490, qsg_492, \
                         qsg_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_18 * osg_370[k]
                   + f_3 * pc_z[k] * qsg_490[k];

        t_689[k] = f_11 * osg_387[k]
                   + f_6 * qsf0_328[k]
                   - f_7 * qsf1_328[k]
                   + f_3 * pc_y[k] * qsg_492[k];

        t_690[k] = f_11 * osg_388[k]
                   + f_4 * qsf0_329[k]
                   - f_5 * qsf1_329[k]
                   + f_3 * pc_y[k] * qsg_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, osg_374, osg_389, osg_495, \
                         qsf0_329, qsf0_330, qsf1_329, qsf1_330, qsg_494, \
                         qsg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * osg_389[k]
                   + f_3 * pc_y[k] * qsg_494[k];

        t_692[k] = f_18 * osg_374[k]
                   + f_1 * qsf0_329[k]
                   - f_2 * qsf1_329[k]
                   + f_3 * pc_z[k] * qsg_494[k];

        t_693[k] = f_20 * osg_495[k]
                   + f_1 * qsf0_330[k]
                   - f_2 * qsf1_330[k]
                   + f_3 * pc_x[k] * qsg_495[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, osg_375, osg_390, \
                         osg_392, osg_498, qsf0_333, qsf1_333, qsg_495, qsg_497, \
                         qsg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * osg_390[k]
                   + f_3 * pc_y[k] * qsg_495[k];

        t_695[k] = f_20 * osg_375[k]
                   + f_3 * pc_z[k] * qsg_495[k];

        t_696[k] = f_20 * osg_498[k]
                   + f_6 * qsf0_333[k]
                   - f_7 * qsf1_333[k]
                   + f_3 * pc_x[k] * qsg_498[k];

        t_697[k] = f_10 * osg_392[k]
                   + f_3 * pc_y[k] * qsg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, osg_378, osg_500, osg_501, qsf0_335, \
                         qsf0_336, qsf1_335, qsf1_336, qsg_498, qsg_500, \
                         qsg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_20 * osg_500[k]
                   + f_6 * qsf0_335[k]
                   - f_7 * qsf1_335[k]
                   + f_3 * pc_x[k] * qsg_500[k];

        t_699[k] = f_20 * osg_501[k]
                   + f_4 * qsf0_336[k]
                   - f_5 * qsf1_336[k]
                   + f_3 * pc_x[k] * qsg_501[k];

        t_700[k] = f_20 * osg_378[k]
                   + f_3 * pc_z[k] * qsg_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, osg_395, osg_504, osg_505, \
                         osg_506, qsf0_339, qsf1_339, qsg_500, qsg_504, qsg_505, \
                         qsg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * osg_395[k]
                   + f_3 * pc_y[k] * qsg_500[k];

        t_702[k] = f_20 * osg_504[k]
                   + f_4 * qsf0_339[k]
                   - f_5 * qsf1_339[k]
                   + f_3 * pc_x[k] * qsg_504[k];

        t_703[k] = f_20 * osg_505[k]
                   + f_3 * pc_x[k] * qsg_505[k];

        t_704[k] = f_20 * osg_506[k]
                   + f_3 * pc_x[k] * qsg_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, osg_400, osg_507, osg_508, \
                         osg_509, qsf0_336, qsf1_336, qsg_505, qsg_507, qsg_508, \
                         qsg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_20 * osg_507[k]
                   + f_3 * pc_x[k] * qsg_507[k];

        t_706[k] = f_20 * osg_508[k]
                   + f_3 * pc_x[k] * qsg_508[k];

        t_707[k] = f_20 * osg_509[k]
                   + f_3 * pc_x[k] * qsg_509[k];

        t_708[k] = f_10 * osg_400[k]
                   + f_1 * qsf0_336[k]
                   - f_2 * qsf1_336[k]
                   + f_3 * pc_y[k] * qsg_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, osg_385, osg_402, osg_403, qsf0_338, \
                         qsf0_339, qsf1_338, qsf1_339, qsg_505, qsg_507, \
                         qsg_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_20 * osg_385[k]
                   + f_3 * pc_z[k] * qsg_505[k];

        t_710[k] = f_10 * osg_402[k]
                   + f_6 * qsf0_338[k]
                   - f_7 * qsf1_338[k]
                   + f_3 * pc_y[k] * qsg_507[k];

        t_711[k] = f_10 * osg_403[k]
                   + f_4 * qsf0_339[k]
                   - f_5 * qsf1_339[k]
                   + f_3 * pc_y[k] * qsg_508[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *osh0_567 = buffer.data(osh0 + 567);
    const auto *osh0_570 = buffer.data(osh0 + 570);
    const auto *osh0_572 = buffer.data(osh0 + 572);
    const auto *osh0_573 = buffer.data(osh0 + 573);
    const auto *osh0_576 = buffer.data(osh0 + 576);
    const auto *osh0_587 = buffer.data(osh0 + 587);
    const auto *osh0_588 = buffer.data(osh0 + 588);
    const auto *osh0_591 = buffer.data(osh0 + 591);
    const auto *osh0_594 = buffer.data(osh0 + 594);
    const auto *osh0_603 = buffer.data(osh0 + 603);

    const auto *osg_389 = buffer.data(osg + 389);
    const auto *osg_390 = buffer.data(osg + 390);
    const auto *osg_393 = buffer.data(osg + 393);
    const auto *osg_400 = buffer.data(osg + 400);
    const auto *osg_404 = buffer.data(osg + 404);
    const auto *osg_405 = buffer.data(osg + 405);
    const auto *osg_406 = buffer.data(osg + 406);
    const auto *osg_407 = buffer.data(osg + 407);
    const auto *osg_408 = buffer.data(osg + 408);
    const auto *osg_410 = buffer.data(osg + 410);
    const auto *osg_415 = buffer.data(osg + 415);
    const auto *osg_417 = buffer.data(osg + 417);
    const auto *osg_418 = buffer.data(osg + 418);
    const auto *osg_419 = buffer.data(osg + 419);
    const auto *osg_420 = buffer.data(osg + 420);
    const auto *osg_423 = buffer.data(osg + 423);
    const auto *osg_425 = buffer.data(osg + 425);
    const auto *osg_430 = buffer.data(osg + 430);
    const auto *osg_434 = buffer.data(osg + 434);
    const auto *osg_435 = buffer.data(osg + 435);
    const auto *osg_437 = buffer.data(osg + 437);
    const auto *osg_438 = buffer.data(osg + 438);
    const auto *osg_440 = buffer.data(osg + 440);
    const auto *osg_445 = buffer.data(osg + 445);
    const auto *osg_447 = buffer.data(osg + 447);
    const auto *osg_448 = buffer.data(osg + 448);
    const auto *osg_449 = buffer.data(osg + 449);
    const auto *osg_450 = buffer.data(osg + 450);
    const auto *osg_452 = buffer.data(osg + 452);
    const auto *osg_453 = buffer.data(osg + 453);
    const auto *osg_455 = buffer.data(osg + 455);
    const auto *osg_460 = buffer.data(osg + 460);
    const auto *osg_462 = buffer.data(osg + 462);
    const auto *osg_463 = buffer.data(osg + 463);
    const auto *osg_464 = buffer.data(osg + 464);
    const auto *osg_465 = buffer.data(osg + 465);
    const auto *osg_467 = buffer.data(osg + 467);
    const auto *osg_470 = buffer.data(osg + 470);
    const auto *osg_520 = buffer.data(osg + 520);
    const auto *osg_521 = buffer.data(osg + 521);
    const auto *osg_522 = buffer.data(osg + 522);
    const auto *osg_523 = buffer.data(osg + 523);
    const auto *osg_524 = buffer.data(osg + 524);
    const auto *osg_525 = buffer.data(osg + 525);
    const auto *osg_530 = buffer.data(osg + 530);
    const auto *osg_534 = buffer.data(osg + 534);
    const auto *osg_535 = buffer.data(osg + 535);
    const auto *osg_536 = buffer.data(osg + 536);
    const auto *osg_537 = buffer.data(osg + 537);
    const auto *osg_539 = buffer.data(osg + 539);
    const auto *osg_540 = buffer.data(osg + 540);
    const auto *osg_543 = buffer.data(osg + 543);
    const auto *osg_546 = buffer.data(osg + 546);
    const auto *osg_550 = buffer.data(osg + 550);
    const auto *osg_552 = buffer.data(osg + 552);
    const auto *osg_553 = buffer.data(osg + 553);
    const auto *osg_554 = buffer.data(osg + 554);
    const auto *osg_560 = buffer.data(osg + 560);
    const auto *osg_564 = buffer.data(osg + 564);
    const auto *osg_565 = buffer.data(osg + 565);
    const auto *osg_566 = buffer.data(osg + 566);
    const auto *osg_567 = buffer.data(osg + 567);
    const auto *osg_568 = buffer.data(osg + 568);
    const auto *osg_569 = buffer.data(osg + 569);
    const auto *osg_570 = buffer.data(osg + 570);
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
    const auto *osg_588 = buffer.data(osg + 588);
    const auto *osg_590 = buffer.data(osg + 590);
    const auto *osg_591 = buffer.data(osg + 591);
    const auto *osg_594 = buffer.data(osg + 594);
    const auto *osg_595 = buffer.data(osg + 595);
    const auto *osg_596 = buffer.data(osg + 596);

    const auto *osh1_567 = buffer.data(osh1 + 567);
    const auto *osh1_570 = buffer.data(osh1 + 570);
    const auto *osh1_572 = buffer.data(osh1 + 572);
    const auto *osh1_573 = buffer.data(osh1 + 573);
    const auto *osh1_576 = buffer.data(osh1 + 576);
    const auto *osh1_587 = buffer.data(osh1 + 587);
    const auto *osh1_588 = buffer.data(osh1 + 588);
    const auto *osh1_591 = buffer.data(osh1 + 591);
    const auto *osh1_594 = buffer.data(osh1 + 594);
    const auto *osh1_603 = buffer.data(osh1 + 603);

    const auto *qsf0_339 = buffer.data(qsf0 + 339);
    const auto *qsf0_346 = buffer.data(qsf0 + 346);
    const auto *qsf0_348 = buffer.data(qsf0 + 348);
    const auto *qsf0_349 = buffer.data(qsf0 + 349);
    const auto *qsf0_350 = buffer.data(qsf0 + 350);
    const auto *qsf0_351 = buffer.data(qsf0 + 351);
    const auto *qsf0_352 = buffer.data(qsf0 + 352);
    const auto *qsf0_355 = buffer.data(qsf0 + 355);
    const auto *qsf0_356 = buffer.data(qsf0 + 356);
    const auto *qsf0_357 = buffer.data(qsf0 + 357);
    const auto *qsf0_358 = buffer.data(qsf0 + 358);
    const auto *qsf0_359 = buffer.data(qsf0 + 359);
    const auto *qsf0_360 = buffer.data(qsf0 + 360);
    const auto *qsf0_362 = buffer.data(qsf0 + 362);
    const auto *qsf0_363 = buffer.data(qsf0 + 363);
    const auto *qsf0_366 = buffer.data(qsf0 + 366);
    const auto *qsf0_367 = buffer.data(qsf0 + 367);
    const auto *qsf0_369 = buffer.data(qsf0 + 369);
    const auto *qsf0_375 = buffer.data(qsf0 + 375);
    const auto *qsf0_378 = buffer.data(qsf0 + 378);
    const auto *qsf0_379 = buffer.data(qsf0 + 379);
    const auto *qsf0_380 = buffer.data(qsf0 + 380);
    const auto *qsf0_383 = buffer.data(qsf0 + 383);
    const auto *qsf0_385 = buffer.data(qsf0 + 385);
    const auto *qsf0_386 = buffer.data(qsf0 + 386);
    const auto *qsf0_388 = buffer.data(qsf0 + 388);
    const auto *qsf0_389 = buffer.data(qsf0 + 389);
    const auto *qsf0_390 = buffer.data(qsf0 + 390);
    const auto *qsf0_393 = buffer.data(qsf0 + 393);
    const auto *qsf0_395 = buffer.data(qsf0 + 395);
    const auto *qsf0_396 = buffer.data(qsf0 + 396);
    const auto *qsf0_399 = buffer.data(qsf0 + 399);

    const auto *qsf1_339 = buffer.data(qsf1 + 339);
    const auto *qsf1_346 = buffer.data(qsf1 + 346);
    const auto *qsf1_348 = buffer.data(qsf1 + 348);
    const auto *qsf1_349 = buffer.data(qsf1 + 349);
    const auto *qsf1_350 = buffer.data(qsf1 + 350);
    const auto *qsf1_351 = buffer.data(qsf1 + 351);
    const auto *qsf1_352 = buffer.data(qsf1 + 352);
    const auto *qsf1_355 = buffer.data(qsf1 + 355);
    const auto *qsf1_356 = buffer.data(qsf1 + 356);
    const auto *qsf1_357 = buffer.data(qsf1 + 357);
    const auto *qsf1_358 = buffer.data(qsf1 + 358);
    const auto *qsf1_359 = buffer.data(qsf1 + 359);
    const auto *qsf1_360 = buffer.data(qsf1 + 360);
    const auto *qsf1_362 = buffer.data(qsf1 + 362);
    const auto *qsf1_363 = buffer.data(qsf1 + 363);
    const auto *qsf1_366 = buffer.data(qsf1 + 366);
    const auto *qsf1_367 = buffer.data(qsf1 + 367);
    const auto *qsf1_369 = buffer.data(qsf1 + 369);
    const auto *qsf1_375 = buffer.data(qsf1 + 375);
    const auto *qsf1_378 = buffer.data(qsf1 + 378);
    const auto *qsf1_379 = buffer.data(qsf1 + 379);
    const auto *qsf1_380 = buffer.data(qsf1 + 380);
    const auto *qsf1_383 = buffer.data(qsf1 + 383);
    const auto *qsf1_385 = buffer.data(qsf1 + 385);
    const auto *qsf1_386 = buffer.data(qsf1 + 386);
    const auto *qsf1_388 = buffer.data(qsf1 + 388);
    const auto *qsf1_389 = buffer.data(qsf1 + 389);
    const auto *qsf1_390 = buffer.data(qsf1 + 390);
    const auto *qsf1_393 = buffer.data(qsf1 + 393);
    const auto *qsf1_395 = buffer.data(qsf1 + 395);
    const auto *qsf1_396 = buffer.data(qsf1 + 396);
    const auto *qsf1_399 = buffer.data(qsf1 + 399);

    const auto *qsg_509 = buffer.data(qsg + 509);
    const auto *qsg_510 = buffer.data(qsg + 510);
    const auto *qsg_512 = buffer.data(qsg + 512);
    const auto *qsg_513 = buffer.data(qsg + 513);
    const auto *qsg_515 = buffer.data(qsg + 515);
    const auto *qsg_520 = buffer.data(qsg + 520);
    const auto *qsg_521 = buffer.data(qsg + 521);
    const auto *qsg_522 = buffer.data(qsg + 522);
    const auto *qsg_523 = buffer.data(qsg + 523);
    const auto *qsg_524 = buffer.data(qsg + 524);
    const auto *qsg_525 = buffer.data(qsg + 525);
    const auto *qsg_526 = buffer.data(qsg + 526);
    const auto *qsg_527 = buffer.data(qsg + 527);
    const auto *qsg_528 = buffer.data(qsg + 528);
    const auto *qsg_529 = buffer.data(qsg + 529);
    const auto *qsg_530 = buffer.data(qsg + 530);
    const auto *qsg_534 = buffer.data(qsg + 534);
    const auto *qsg_535 = buffer.data(qsg + 535);
    const auto *qsg_536 = buffer.data(qsg + 536);
    const auto *qsg_537 = buffer.data(qsg + 537);
    const auto *qsg_538 = buffer.data(qsg + 538);
    const auto *qsg_539 = buffer.data(qsg + 539);
    const auto *qsg_540 = buffer.data(qsg + 540);
    const auto *qsg_541 = buffer.data(qsg + 541);
    const auto *qsg_542 = buffer.data(qsg + 542);
    const auto *qsg_543 = buffer.data(qsg + 543);
    const auto *qsg_545 = buffer.data(qsg + 545);
    const auto *qsg_546 = buffer.data(qsg + 546);
    const auto *qsg_550 = buffer.data(qsg + 550);
    const auto *qsg_551 = buffer.data(qsg + 551);
    const auto *qsg_552 = buffer.data(qsg + 552);
    const auto *qsg_553 = buffer.data(qsg + 553);
    const auto *qsg_554 = buffer.data(qsg + 554);
    const auto *qsg_555 = buffer.data(qsg + 555);
    const auto *qsg_557 = buffer.data(qsg + 557);
    const auto *qsg_558 = buffer.data(qsg + 558);
    const auto *qsg_560 = buffer.data(qsg + 560);
    const auto *qsg_564 = buffer.data(qsg + 564);
    const auto *qsg_565 = buffer.data(qsg + 565);
    const auto *qsg_566 = buffer.data(qsg + 566);
    const auto *qsg_567 = buffer.data(qsg + 567);
    const auto *qsg_568 = buffer.data(qsg + 568);
    const auto *qsg_569 = buffer.data(qsg + 569);
    const auto *qsg_570 = buffer.data(qsg + 570);
    const auto *qsg_572 = buffer.data(qsg + 572);
    const auto *qsg_573 = buffer.data(qsg + 573);
    const auto *qsg_575 = buffer.data(qsg + 575);
    const auto *qsg_576 = buffer.data(qsg + 576);
    const auto *qsg_579 = buffer.data(qsg + 579);
    const auto *qsg_580 = buffer.data(qsg + 580);
    const auto *qsg_581 = buffer.data(qsg + 581);
    const auto *qsg_582 = buffer.data(qsg + 582);
    const auto *qsg_583 = buffer.data(qsg + 583);
    const auto *qsg_584 = buffer.data(qsg + 584);
    const auto *qsg_585 = buffer.data(qsg + 585);
    const auto *qsg_587 = buffer.data(qsg + 587);
    const auto *qsg_588 = buffer.data(qsg + 588);
    const auto *qsg_590 = buffer.data(qsg + 590);
    const auto *qsg_591 = buffer.data(qsg + 591);
    const auto *qsg_594 = buffer.data(qsg + 594);
    const auto *qsg_595 = buffer.data(qsg + 595);
    const auto *qsg_596 = buffer.data(qsg + 596);

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_y, pc_y, pc_z, osh0_567, osg_389, \
                         osg_404, osg_405, osh1_567, qsf0_339, qsf1_339, qsg_509, \
                         qsg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * osg_404[k]
                   + f_3 * pc_y[k] * qsg_509[k];

        t_713[k] = f_20 * osg_389[k]
                   + f_1 * qsf0_339[k]
                   - f_2 * qsf1_339[k]
                   + f_3 * pc_z[k] * qsg_509[k];

        t_714[k] = pa_y[k] * osh0_567[k]
                   - f_8 * pc_y[k] * osh1_567[k];

        t_715[k] = f_9 * osg_405[k]
                   + f_3 * pc_y[k] * qsg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_y, pc_y, pc_z, osh0_570, osh0_572, \
                         osg_390, osg_406, osg_407, osh1_570, osh1_572, qsg_510, \
                         qsg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_21 * osg_390[k]
                   + f_3 * pc_z[k] * qsg_510[k];

        t_717[k] = pa_y[k] * osh0_570[k]
                   + f_10 * osg_406[k]
                   - f_8 * pc_y[k] * osh1_570[k];

        t_718[k] = f_9 * osg_407[k]
                   + f_3 * pc_y[k] * qsg_512[k];

        t_719[k] = pa_y[k] * osh0_572[k]
                   - f_8 * pc_y[k] * osh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_y, pc_y, pc_z, osh0_573, osh0_576, \
                         osg_393, osg_408, osg_410, osh1_573, osh1_576, qsg_513, \
                         qsg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_y[k] * osh0_573[k]
                   + f_11 * osg_408[k]
                   - f_8 * pc_y[k] * osh1_573[k];

        t_721[k] = f_21 * osg_393[k]
                   + f_3 * pc_z[k] * qsg_513[k];

        t_722[k] = f_9 * osg_410[k]
                   + f_3 * pc_y[k] * qsg_515[k];

        t_723[k] = pa_y[k] * osh0_576[k]
                   - f_8 * pc_y[k] * osh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, osg_520, osg_521, osg_522, \
                         osg_523, osg_524, qsg_520, qsg_521, qsg_522, qsg_523, \
                         qsg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_20 * osg_520[k]
                   + f_3 * pc_x[k] * qsg_520[k];

        t_725[k] = f_20 * osg_521[k]
                   + f_3 * pc_x[k] * qsg_521[k];

        t_726[k] = f_20 * osg_522[k]
                   + f_3 * pc_x[k] * qsg_522[k];

        t_727[k] = f_20 * osg_523[k]
                   + f_3 * pc_x[k] * qsg_523[k];

        t_728[k] = f_20 * osg_524[k]
                   + f_3 * pc_x[k] * qsg_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, osg_400, osg_415, osg_417, qsf0_346, \
                         qsf0_348, qsf1_346, qsf1_348, qsg_520, \
                         qsg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * osg_415[k]
                   + f_1 * qsf0_346[k]
                   - f_2 * qsf1_346[k]
                   + f_3 * pc_y[k] * qsg_520[k];

        t_730[k] = f_21 * osg_400[k]
                   + f_3 * pc_z[k] * qsg_520[k];

        t_731[k] = f_9 * osg_417[k]
                   + f_6 * qsf0_348[k]
                   - f_7 * qsf1_348[k]
                   + f_3 * pc_y[k] * qsg_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pc_y, osh0_587, osg_418, osg_419, \
                         osh1_587, qsf0_349, qsf1_349, qsg_523, \
                         qsg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * osg_418[k]
                   + f_4 * qsf0_349[k]
                   - f_5 * qsf1_349[k]
                   + f_3 * pc_y[k] * qsg_523[k];

        t_733[k] = f_9 * osg_419[k]
                   + f_3 * pc_y[k] * qsg_524[k];

        t_734[k] = pa_y[k] * osh0_587[k]
                   - f_8 * pc_y[k] * osh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, osg_405, \
                         osg_525, qsf0_350, qsf1_350, qsg_525, qsg_526, \
                         qsg_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_20 * osg_525[k]
                   + f_1 * qsf0_350[k]
                   - f_2 * qsf1_350[k]
                   + f_3 * pc_x[k] * qsg_525[k];

        t_736[k] = f_3 * pc_y[k] * qsg_525[k];

        t_737[k] = f_19 * osg_405[k]
                   + f_3 * pc_z[k] * qsg_525[k];

        t_738[k] = f_4 * qsf0_350[k]
                   - f_5 * qsf1_350[k]
                   + f_3 * pc_y[k] * qsg_526[k];

        t_739[k] = f_3 * pc_y[k] * qsg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pc_y, osg_530, qsf0_351, qsf0_352, \
                         qsf0_355, qsf1_351, qsf1_352, qsf1_355, qsg_528, qsg_529, \
                         qsg_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_20 * osg_530[k]
                   + f_6 * qsf0_355[k]
                   - f_7 * qsf1_355[k]
                   + f_3 * pc_x[k] * qsg_530[k];

        t_741[k] = f_6 * qsf0_351[k]
                   - f_7 * qsf1_351[k]
                   + f_3 * pc_y[k] * qsg_528[k];

        t_742[k] = f_4 * qsf0_352[k]
                   - f_5 * qsf1_352[k]
                   + f_3 * pc_y[k] * qsg_529[k];

        t_743[k] = f_3 * pc_y[k] * qsg_530[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pc_x, osg_534, osg_535, osg_536, osg_537, \
                         qsf0_359, qsf1_359, qsg_534, qsg_535, qsg_536, \
                         qsg_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_20 * osg_534[k]
                   + f_4 * qsf0_359[k]
                   - f_5 * qsf1_359[k]
                   + f_3 * pc_x[k] * qsg_534[k];

        t_745[k] = f_20 * osg_535[k]
                   + f_3 * pc_x[k] * qsg_535[k];

        t_746[k] = f_20 * osg_536[k]
                   + f_3 * pc_x[k] * qsg_536[k];

        t_747[k] = f_20 * osg_537[k]
                   + f_3 * pc_x[k] * qsg_537[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_x, pc_y, osg_539, qsf0_356, qsf0_357, \
                         qsf1_356, qsf1_357, qsg_534, qsg_535, qsg_536, \
                         qsg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_3 * pc_y[k] * qsg_534[k];

        t_749[k] = f_20 * osg_539[k]
                   + f_3 * pc_x[k] * qsg_539[k];

        t_750[k] = f_1 * qsf0_356[k]
                   - f_2 * qsf1_356[k]
                   + f_3 * pc_y[k] * qsg_535[k];

        t_751[k] = f_13 * qsf0_357[k]
                   - f_14 * qsf1_357[k]
                   + f_3 * pc_y[k] * qsg_536[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, osg_419, qsf0_358, qsf0_359, \
                         qsf1_358, qsf1_359, qsg_537, qsg_538, \
                         qsg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_6 * qsf0_358[k]
                   - f_7 * qsf1_358[k]
                   + f_3 * pc_y[k] * qsg_537[k];

        t_753[k] = f_4 * qsf0_359[k]
                   - f_5 * qsf1_359[k]
                   + f_3 * pc_y[k] * qsg_538[k];

        t_754[k] = f_3 * pc_y[k] * qsg_539[k];

        t_755[k] = f_19 * osg_419[k]
                   + f_1 * qsf0_359[k]
                   - f_2 * qsf1_359[k]
                   + f_3 * pc_z[k] * qsg_539[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, osg_420, osg_540, \
                         osg_543, qsf0_360, qsf0_363, qsf1_360, qsf1_363, qsg_540, \
                         qsg_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_18 * osg_540[k]
                   + f_1 * qsf0_360[k]
                   - f_2 * qsf1_360[k]
                   + f_3 * pc_x[k] * qsg_540[k];

        t_757[k] = f_17 * osg_420[k]
                   + f_3 * pc_y[k] * qsg_540[k];

        t_758[k] = f_3 * pc_z[k] * qsg_540[k];

        t_759[k] = f_18 * osg_543[k]
                   + f_6 * qsf0_363[k]
                   - f_7 * qsf1_363[k]
                   + f_3 * pc_x[k] * qsg_543[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, pc_x, pc_z, osg_546, qsf0_360, qsf0_366, \
                         qsf1_360, qsf1_366, qsg_541, qsg_542, qsg_543, \
                         qsg_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_3 * pc_z[k] * qsg_541[k];

        t_761[k] = f_4 * qsf0_360[k]
                   - f_5 * qsf1_360[k]
                   + f_3 * pc_z[k] * qsg_542[k];

        t_762[k] = f_18 * osg_546[k]
                   + f_4 * qsf0_366[k]
                   - f_5 * qsf1_366[k]
                   + f_3 * pc_x[k] * qsg_546[k];

        t_763[k] = f_3 * pc_z[k] * qsg_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, osg_425, osg_550, \
                         qsf0_362, qsf1_362, qsg_545, qsg_546, \
                         qsg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_17 * osg_425[k]
                   + f_3 * pc_y[k] * qsg_545[k];

        t_765[k] = f_6 * qsf0_362[k]
                   - f_7 * qsf1_362[k]
                   + f_3 * pc_z[k] * qsg_545[k];

        t_766[k] = f_18 * osg_550[k]
                   + f_3 * pc_x[k] * qsg_550[k];

        t_767[k] = f_3 * pc_z[k] * qsg_546[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, osg_430, osg_552, osg_553, \
                         osg_554, qsf0_366, qsf1_366, qsg_550, qsg_552, qsg_553, \
                         qsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_18 * osg_552[k]
                   + f_3 * pc_x[k] * qsg_552[k];

        t_769[k] = f_18 * osg_553[k]
                   + f_3 * pc_x[k] * qsg_553[k];

        t_770[k] = f_18 * osg_554[k]
                   + f_3 * pc_x[k] * qsg_554[k];

        t_771[k] = f_17 * osg_430[k]
                   + f_1 * qsf0_366[k]
                   - f_2 * qsf1_366[k]
                   + f_3 * pc_y[k] * qsg_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, pc_y, pc_z, osg_434, qsf0_366, qsf0_367, \
                         qsf1_366, qsf1_367, qsg_550, qsg_551, qsg_552, \
                         qsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * qsg_550[k];

        t_773[k] = f_4 * qsf0_366[k]
                   - f_5 * qsf1_366[k]
                   + f_3 * pc_z[k] * qsg_551[k];

        t_774[k] = f_6 * qsf0_367[k]
                   - f_7 * qsf1_367[k]
                   + f_3 * pc_z[k] * qsg_552[k];

        t_775[k] = f_17 * osg_434[k]
                   + f_3 * pc_y[k] * qsg_554[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pa_z, pc_y, pc_z, osh0_588, osg_420, \
                         osg_435, osh1_588, qsf0_369, qsf1_369, qsg_554, \
                         qsg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_1 * qsf0_369[k]
                   - f_2 * qsf1_369[k]
                   + f_3 * pc_z[k] * qsg_554[k];

        t_777[k] = pa_z[k] * osh0_588[k]
                   - f_8 * pc_z[k] * osh1_588[k];

        t_778[k] = f_19 * osg_435[k]
                   + f_3 * pc_y[k] * qsg_555[k];

        t_779[k] = f_9 * osg_420[k]
                   + f_3 * pc_z[k] * qsg_555[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pa_z, pc_x, pc_y, pc_z, osh0_591, osg_437, \
                         osg_560, osh1_591, qsf0_375, qsf1_375, qsg_557, \
                         qsg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pa_z[k] * osh0_591[k]
                   - f_8 * pc_z[k] * osh1_591[k];

        t_781[k] = f_19 * osg_437[k]
                   + f_3 * pc_y[k] * qsg_557[k];

        t_782[k] = f_18 * osg_560[k]
                   + f_6 * qsf0_375[k]
                   - f_7 * qsf1_375[k]
                   + f_3 * pc_x[k] * qsg_560[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pa_z, pc_y, pc_z, osh0_594, osg_423, osg_440, \
                         osh1_594, qsg_558, qsg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pa_z[k] * osh0_594[k]
                   - f_8 * pc_z[k] * osh1_594[k];

        t_784[k] = f_9 * osg_423[k]
                   + f_3 * pc_z[k] * qsg_558[k];

        t_785[k] = f_19 * osg_440[k]
                   + f_3 * pc_y[k] * qsg_560[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pc_x, osg_564, osg_565, osg_566, osg_567, \
                         qsf0_379, qsf1_379, qsg_564, qsg_565, qsg_566, \
                         qsg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_18 * osg_564[k]
                   + f_4 * qsf0_379[k]
                   - f_5 * qsf1_379[k]
                   + f_3 * pc_x[k] * qsg_564[k];

        t_787[k] = f_18 * osg_565[k]
                   + f_3 * pc_x[k] * qsg_565[k];

        t_788[k] = f_18 * osg_566[k]
                   + f_3 * pc_x[k] * qsg_566[k];

        t_789[k] = f_18 * osg_567[k]
                   + f_3 * pc_x[k] * qsg_567[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_x, pc_z, osh0_603, osg_430, \
                         osg_568, osg_569, osh1_603, qsg_565, qsg_568, \
                         qsg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_18 * osg_568[k]
                   + f_3 * pc_x[k] * qsg_568[k];

        t_791[k] = f_18 * osg_569[k]
                   + f_3 * pc_x[k] * qsg_569[k];

        t_792[k] = pa_z[k] * osh0_603[k]
                   - f_8 * pc_z[k] * osh1_603[k];

        t_793[k] = f_9 * osg_430[k]
                   + f_3 * pc_z[k] * qsg_565[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pc_y, osg_447, osg_448, osg_449, qsf0_378, \
                         qsf0_379, qsf1_378, qsf1_379, qsg_567, qsg_568, \
                         qsg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_19 * osg_447[k]
                   + f_6 * qsf0_378[k]
                   - f_7 * qsf1_378[k]
                   + f_3 * pc_y[k] * qsg_567[k];

        t_795[k] = f_19 * osg_448[k]
                   + f_4 * qsf0_379[k]
                   - f_5 * qsf1_379[k]
                   + f_3 * pc_y[k] * qsg_568[k];

        t_796[k] = f_19 * osg_449[k]
                   + f_3 * pc_y[k] * qsg_569[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pc_x, pc_y, pc_z, osg_434, osg_450, osg_570, \
                         qsf0_379, qsf0_380, qsf1_379, qsf1_380, qsg_569, \
                         qsg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_9 * osg_434[k]
                   + f_1 * qsf0_379[k]
                   - f_2 * qsf1_379[k]
                   + f_3 * pc_z[k] * qsg_569[k];

        t_798[k] = f_18 * osg_570[k]
                   + f_1 * qsf0_380[k]
                   - f_2 * qsf1_380[k]
                   + f_3 * pc_x[k] * qsg_570[k];

        t_799[k] = f_21 * osg_450[k]
                   + f_3 * pc_y[k] * qsg_570[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, osg_435, osg_452, osg_573, \
                         qsf0_383, qsf1_383, qsg_570, qsg_572, \
                         qsg_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_10 * osg_435[k]
                   + f_3 * pc_z[k] * qsg_570[k];

        t_801[k] = f_18 * osg_573[k]
                   + f_6 * qsf0_383[k]
                   - f_7 * qsf1_383[k]
                   + f_3 * pc_x[k] * qsg_573[k];

        t_802[k] = f_21 * osg_452[k]
                   + f_3 * pc_y[k] * qsg_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, pc_z, osg_438, osg_575, osg_576, qsf0_385, \
                         qsf0_386, qsf1_385, qsf1_386, qsg_573, qsg_575, \
                         qsg_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_18 * osg_575[k]
                   + f_6 * qsf0_385[k]
                   - f_7 * qsf1_385[k]
                   + f_3 * pc_x[k] * qsg_575[k];

        t_804[k] = f_18 * osg_576[k]
                   + f_4 * qsf0_386[k]
                   - f_5 * qsf1_386[k]
                   + f_3 * pc_x[k] * qsg_576[k];

        t_805[k] = f_10 * osg_438[k]
                   + f_3 * pc_z[k] * qsg_573[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, pc_x, pc_y, osg_455, osg_579, osg_580, \
                         osg_581, qsf0_389, qsf1_389, qsg_575, qsg_579, qsg_580, \
                         qsg_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_21 * osg_455[k]
                   + f_3 * pc_y[k] * qsg_575[k];

        t_807[k] = f_18 * osg_579[k]
                   + f_4 * qsf0_389[k]
                   - f_5 * qsf1_389[k]
                   + f_3 * pc_x[k] * qsg_579[k];

        t_808[k] = f_18 * osg_580[k]
                   + f_3 * pc_x[k] * qsg_580[k];

        t_809[k] = f_18 * osg_581[k]
                   + f_3 * pc_x[k] * qsg_581[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pc_x, pc_y, osg_460, osg_582, osg_583, \
                         osg_584, qsf0_386, qsf1_386, qsg_580, qsg_582, qsg_583, \
                         qsg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_18 * osg_582[k]
                   + f_3 * pc_x[k] * qsg_582[k];

        t_811[k] = f_18 * osg_583[k]
                   + f_3 * pc_x[k] * qsg_583[k];

        t_812[k] = f_18 * osg_584[k]
                   + f_3 * pc_x[k] * qsg_584[k];

        t_813[k] = f_21 * osg_460[k]
                   + f_1 * qsf0_386[k]
                   - f_2 * qsf1_386[k]
                   + f_3 * pc_y[k] * qsg_580[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pc_y, pc_z, osg_445, osg_462, osg_463, qsf0_388, \
                         qsf0_389, qsf1_388, qsf1_389, qsg_580, qsg_582, \
                         qsg_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_10 * osg_445[k]
                   + f_3 * pc_z[k] * qsg_580[k];

        t_815[k] = f_21 * osg_462[k]
                   + f_6 * qsf0_388[k]
                   - f_7 * qsf1_388[k]
                   + f_3 * pc_y[k] * qsg_582[k];

        t_816[k] = f_21 * osg_463[k]
                   + f_4 * qsf0_389[k]
                   - f_5 * qsf1_389[k]
                   + f_3 * pc_y[k] * qsg_583[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pc_x, pc_y, pc_z, osg_449, osg_464, osg_585, \
                         qsf0_389, qsf0_390, qsf1_389, qsf1_390, qsg_584, \
                         qsg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_21 * osg_464[k]
                   + f_3 * pc_y[k] * qsg_584[k];

        t_818[k] = f_10 * osg_449[k]
                   + f_1 * qsf0_389[k]
                   - f_2 * qsf1_389[k]
                   + f_3 * pc_z[k] * qsg_584[k];

        t_819[k] = f_18 * osg_585[k]
                   + f_1 * qsf0_390[k]
                   - f_2 * qsf1_390[k]
                   + f_3 * pc_x[k] * qsg_585[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_x, pc_y, pc_z, osg_450, osg_465, \
                         osg_467, osg_588, qsf0_393, qsf1_393, qsg_585, qsg_587, \
                         qsg_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_20 * osg_465[k]
                   + f_3 * pc_y[k] * qsg_585[k];

        t_821[k] = f_11 * osg_450[k]
                   + f_3 * pc_z[k] * qsg_585[k];

        t_822[k] = f_18 * osg_588[k]
                   + f_6 * qsf0_393[k]
                   - f_7 * qsf1_393[k]
                   + f_3 * pc_x[k] * qsg_588[k];

        t_823[k] = f_20 * osg_467[k]
                   + f_3 * pc_y[k] * qsg_587[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, osg_453, osg_590, osg_591, qsf0_395, \
                         qsf0_396, qsf1_395, qsf1_396, qsg_588, qsg_590, \
                         qsg_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_18 * osg_590[k]
                   + f_6 * qsf0_395[k]
                   - f_7 * qsf1_395[k]
                   + f_3 * pc_x[k] * qsg_590[k];

        t_825[k] = f_18 * osg_591[k]
                   + f_4 * qsf0_396[k]
                   - f_5 * qsf1_396[k]
                   + f_3 * pc_x[k] * qsg_591[k];

        t_826[k] = f_11 * osg_453[k]
                   + f_3 * pc_z[k] * qsg_588[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pc_x, pc_y, osg_470, osg_594, osg_595, \
                         osg_596, qsf0_399, qsf1_399, qsg_590, qsg_594, qsg_595, \
                         qsg_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_20 * osg_470[k]
                   + f_3 * pc_y[k] * qsg_590[k];

        t_828[k] = f_18 * osg_594[k]
                   + f_4 * qsf0_399[k]
                   - f_5 * qsf1_399[k]
                   + f_3 * pc_x[k] * qsg_594[k];

        t_829[k] = f_18 * osg_595[k]
                   + f_3 * pc_x[k] * qsg_595[k];

        t_830[k] = f_18 * osg_596[k]
                   + f_3 * pc_x[k] * qsg_596[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *osh0_735 = buffer.data(osh0 + 735);
    const auto *osh0_738 = buffer.data(osh0 + 738);
    const auto *osh0_740 = buffer.data(osh0 + 740);
    const auto *osh0_741 = buffer.data(osh0 + 741);
    const auto *osh0_744 = buffer.data(osh0 + 744);
    const auto *osh0_755 = buffer.data(osh0 + 755);

    const auto *osg_460 = buffer.data(osg + 460);
    const auto *osg_464 = buffer.data(osg + 464);
    const auto *osg_465 = buffer.data(osg + 465);
    const auto *osg_468 = buffer.data(osg + 468);
    const auto *osg_475 = buffer.data(osg + 475);
    const auto *osg_477 = buffer.data(osg + 477);
    const auto *osg_478 = buffer.data(osg + 478);
    const auto *osg_479 = buffer.data(osg + 479);
    const auto *osg_480 = buffer.data(osg + 480);
    const auto *osg_482 = buffer.data(osg + 482);
    const auto *osg_483 = buffer.data(osg + 483);
    const auto *osg_485 = buffer.data(osg + 485);
    const auto *osg_490 = buffer.data(osg + 490);
    const auto *osg_492 = buffer.data(osg + 492);
    const auto *osg_493 = buffer.data(osg + 493);
    const auto *osg_494 = buffer.data(osg + 494);
    const auto *osg_495 = buffer.data(osg + 495);
    const auto *osg_497 = buffer.data(osg + 497);
    const auto *osg_498 = buffer.data(osg + 498);
    const auto *osg_500 = buffer.data(osg + 500);
    const auto *osg_505 = buffer.data(osg + 505);
    const auto *osg_507 = buffer.data(osg + 507);
    const auto *osg_508 = buffer.data(osg + 508);
    const auto *osg_509 = buffer.data(osg + 509);
    const auto *osg_510 = buffer.data(osg + 510);
    const auto *osg_512 = buffer.data(osg + 512);
    const auto *osg_513 = buffer.data(osg + 513);
    const auto *osg_515 = buffer.data(osg + 515);
    const auto *osg_520 = buffer.data(osg + 520);
    const auto *osg_522 = buffer.data(osg + 522);
    const auto *osg_523 = buffer.data(osg + 523);
    const auto *osg_524 = buffer.data(osg + 524);
    const auto *osg_525 = buffer.data(osg + 525);
    const auto *osg_526 = buffer.data(osg + 526);
    const auto *osg_527 = buffer.data(osg + 527);
    const auto *osg_528 = buffer.data(osg + 528);
    const auto *osg_530 = buffer.data(osg + 530);
    const auto *osg_535 = buffer.data(osg + 535);
    const auto *osg_537 = buffer.data(osg + 537);
    const auto *osg_538 = buffer.data(osg + 538);
    const auto *osg_539 = buffer.data(osg + 539);
    const auto *osg_597 = buffer.data(osg + 597);
    const auto *osg_598 = buffer.data(osg + 598);
    const auto *osg_599 = buffer.data(osg + 599);
    const auto *osg_600 = buffer.data(osg + 600);
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
    const auto *osg_633 = buffer.data(osg + 633);
    const auto *osg_635 = buffer.data(osg + 635);
    const auto *osg_636 = buffer.data(osg + 636);
    const auto *osg_639 = buffer.data(osg + 639);
    const auto *osg_640 = buffer.data(osg + 640);
    const auto *osg_641 = buffer.data(osg + 641);
    const auto *osg_642 = buffer.data(osg + 642);
    const auto *osg_643 = buffer.data(osg + 643);
    const auto *osg_644 = buffer.data(osg + 644);
    const auto *osg_655 = buffer.data(osg + 655);
    const auto *osg_656 = buffer.data(osg + 656);
    const auto *osg_657 = buffer.data(osg + 657);
    const auto *osg_658 = buffer.data(osg + 658);
    const auto *osg_659 = buffer.data(osg + 659);
    const auto *osg_660 = buffer.data(osg + 660);
    const auto *osg_665 = buffer.data(osg + 665);
    const auto *osg_669 = buffer.data(osg + 669);
    const auto *osg_670 = buffer.data(osg + 670);
    const auto *osg_671 = buffer.data(osg + 671);
    const auto *osg_672 = buffer.data(osg + 672);
    const auto *osg_674 = buffer.data(osg + 674);

    const auto *osh1_735 = buffer.data(osh1 + 735);
    const auto *osh1_738 = buffer.data(osh1 + 738);
    const auto *osh1_740 = buffer.data(osh1 + 740);
    const auto *osh1_741 = buffer.data(osh1 + 741);
    const auto *osh1_744 = buffer.data(osh1 + 744);
    const auto *osh1_755 = buffer.data(osh1 + 755);

    const auto *qsf0_396 = buffer.data(qsf0 + 396);
    const auto *qsf0_398 = buffer.data(qsf0 + 398);
    const auto *qsf0_399 = buffer.data(qsf0 + 399);
    const auto *qsf0_400 = buffer.data(qsf0 + 400);
    const auto *qsf0_403 = buffer.data(qsf0 + 403);
    const auto *qsf0_405 = buffer.data(qsf0 + 405);
    const auto *qsf0_406 = buffer.data(qsf0 + 406);
    const auto *qsf0_408 = buffer.data(qsf0 + 408);
    const auto *qsf0_409 = buffer.data(qsf0 + 409);
    const auto *qsf0_410 = buffer.data(qsf0 + 410);
    const auto *qsf0_413 = buffer.data(qsf0 + 413);
    const auto *qsf0_415 = buffer.data(qsf0 + 415);
    const auto *qsf0_416 = buffer.data(qsf0 + 416);
    const auto *qsf0_418 = buffer.data(qsf0 + 418);
    const auto *qsf0_419 = buffer.data(qsf0 + 419);
    const auto *qsf0_420 = buffer.data(qsf0 + 420);
    const auto *qsf0_423 = buffer.data(qsf0 + 423);
    const auto *qsf0_425 = buffer.data(qsf0 + 425);
    const auto *qsf0_426 = buffer.data(qsf0 + 426);
    const auto *qsf0_428 = buffer.data(qsf0 + 428);
    const auto *qsf0_429 = buffer.data(qsf0 + 429);
    const auto *qsf0_436 = buffer.data(qsf0 + 436);
    const auto *qsf0_438 = buffer.data(qsf0 + 438);
    const auto *qsf0_439 = buffer.data(qsf0 + 439);
    const auto *qsf0_440 = buffer.data(qsf0 + 440);
    const auto *qsf0_441 = buffer.data(qsf0 + 441);
    const auto *qsf0_442 = buffer.data(qsf0 + 442);
    const auto *qsf0_445 = buffer.data(qsf0 + 445);
    const auto *qsf0_446 = buffer.data(qsf0 + 446);
    const auto *qsf0_447 = buffer.data(qsf0 + 447);
    const auto *qsf0_448 = buffer.data(qsf0 + 448);
    const auto *qsf0_449 = buffer.data(qsf0 + 449);

    const auto *qsf1_396 = buffer.data(qsf1 + 396);
    const auto *qsf1_398 = buffer.data(qsf1 + 398);
    const auto *qsf1_399 = buffer.data(qsf1 + 399);
    const auto *qsf1_400 = buffer.data(qsf1 + 400);
    const auto *qsf1_403 = buffer.data(qsf1 + 403);
    const auto *qsf1_405 = buffer.data(qsf1 + 405);
    const auto *qsf1_406 = buffer.data(qsf1 + 406);
    const auto *qsf1_408 = buffer.data(qsf1 + 408);
    const auto *qsf1_409 = buffer.data(qsf1 + 409);
    const auto *qsf1_410 = buffer.data(qsf1 + 410);
    const auto *qsf1_413 = buffer.data(qsf1 + 413);
    const auto *qsf1_415 = buffer.data(qsf1 + 415);
    const auto *qsf1_416 = buffer.data(qsf1 + 416);
    const auto *qsf1_418 = buffer.data(qsf1 + 418);
    const auto *qsf1_419 = buffer.data(qsf1 + 419);
    const auto *qsf1_420 = buffer.data(qsf1 + 420);
    const auto *qsf1_423 = buffer.data(qsf1 + 423);
    const auto *qsf1_425 = buffer.data(qsf1 + 425);
    const auto *qsf1_426 = buffer.data(qsf1 + 426);
    const auto *qsf1_428 = buffer.data(qsf1 + 428);
    const auto *qsf1_429 = buffer.data(qsf1 + 429);
    const auto *qsf1_436 = buffer.data(qsf1 + 436);
    const auto *qsf1_438 = buffer.data(qsf1 + 438);
    const auto *qsf1_439 = buffer.data(qsf1 + 439);
    const auto *qsf1_440 = buffer.data(qsf1 + 440);
    const auto *qsf1_441 = buffer.data(qsf1 + 441);
    const auto *qsf1_442 = buffer.data(qsf1 + 442);
    const auto *qsf1_445 = buffer.data(qsf1 + 445);
    const auto *qsf1_446 = buffer.data(qsf1 + 446);
    const auto *qsf1_447 = buffer.data(qsf1 + 447);
    const auto *qsf1_448 = buffer.data(qsf1 + 448);
    const auto *qsf1_449 = buffer.data(qsf1 + 449);

    const auto *qsg_595 = buffer.data(qsg + 595);
    const auto *qsg_597 = buffer.data(qsg + 597);
    const auto *qsg_598 = buffer.data(qsg + 598);
    const auto *qsg_599 = buffer.data(qsg + 599);
    const auto *qsg_600 = buffer.data(qsg + 600);
    const auto *qsg_602 = buffer.data(qsg + 602);
    const auto *qsg_603 = buffer.data(qsg + 603);
    const auto *qsg_605 = buffer.data(qsg + 605);
    const auto *qsg_606 = buffer.data(qsg + 606);
    const auto *qsg_609 = buffer.data(qsg + 609);
    const auto *qsg_610 = buffer.data(qsg + 610);
    const auto *qsg_611 = buffer.data(qsg + 611);
    const auto *qsg_612 = buffer.data(qsg + 612);
    const auto *qsg_613 = buffer.data(qsg + 613);
    const auto *qsg_614 = buffer.data(qsg + 614);
    const auto *qsg_615 = buffer.data(qsg + 615);
    const auto *qsg_617 = buffer.data(qsg + 617);
    const auto *qsg_618 = buffer.data(qsg + 618);
    const auto *qsg_620 = buffer.data(qsg + 620);
    const auto *qsg_621 = buffer.data(qsg + 621);
    const auto *qsg_624 = buffer.data(qsg + 624);
    const auto *qsg_625 = buffer.data(qsg + 625);
    const auto *qsg_626 = buffer.data(qsg + 626);
    const auto *qsg_627 = buffer.data(qsg + 627);
    const auto *qsg_628 = buffer.data(qsg + 628);
    const auto *qsg_629 = buffer.data(qsg + 629);
    const auto *qsg_630 = buffer.data(qsg + 630);
    const auto *qsg_632 = buffer.data(qsg + 632);
    const auto *qsg_633 = buffer.data(qsg + 633);
    const auto *qsg_635 = buffer.data(qsg + 635);
    const auto *qsg_636 = buffer.data(qsg + 636);
    const auto *qsg_639 = buffer.data(qsg + 639);
    const auto *qsg_640 = buffer.data(qsg + 640);
    const auto *qsg_641 = buffer.data(qsg + 641);
    const auto *qsg_642 = buffer.data(qsg + 642);
    const auto *qsg_643 = buffer.data(qsg + 643);
    const auto *qsg_644 = buffer.data(qsg + 644);
    const auto *qsg_645 = buffer.data(qsg + 645);
    const auto *qsg_647 = buffer.data(qsg + 647);
    const auto *qsg_648 = buffer.data(qsg + 648);
    const auto *qsg_650 = buffer.data(qsg + 650);
    const auto *qsg_655 = buffer.data(qsg + 655);
    const auto *qsg_656 = buffer.data(qsg + 656);
    const auto *qsg_657 = buffer.data(qsg + 657);
    const auto *qsg_658 = buffer.data(qsg + 658);
    const auto *qsg_659 = buffer.data(qsg + 659);
    const auto *qsg_660 = buffer.data(qsg + 660);
    const auto *qsg_661 = buffer.data(qsg + 661);
    const auto *qsg_662 = buffer.data(qsg + 662);
    const auto *qsg_663 = buffer.data(qsg + 663);
    const auto *qsg_664 = buffer.data(qsg + 664);
    const auto *qsg_665 = buffer.data(qsg + 665);
    const auto *qsg_669 = buffer.data(qsg + 669);
    const auto *qsg_670 = buffer.data(qsg + 670);
    const auto *qsg_671 = buffer.data(qsg + 671);
    const auto *qsg_672 = buffer.data(qsg + 672);
    const auto *qsg_673 = buffer.data(qsg + 673);
    const auto *qsg_674 = buffer.data(qsg + 674);

#pragma omp simd aligned(t_831, t_832, t_833, t_834, pc_x, pc_y, osg_475, osg_597, osg_598, \
                         osg_599, qsf0_396, qsf1_396, qsg_595, qsg_597, qsg_598, \
                         qsg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_18 * osg_597[k]
                   + f_3 * pc_x[k] * qsg_597[k];

        t_832[k] = f_18 * osg_598[k]
                   + f_3 * pc_x[k] * qsg_598[k];

        t_833[k] = f_18 * osg_599[k]
                   + f_3 * pc_x[k] * qsg_599[k];

        t_834[k] = f_20 * osg_475[k]
                   + f_1 * qsf0_396[k]
                   - f_2 * qsf1_396[k]
                   + f_3 * pc_y[k] * qsg_595[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pc_y, pc_z, osg_460, osg_477, osg_478, qsf0_398, \
                         qsf0_399, qsf1_398, qsf1_399, qsg_595, qsg_597, \
                         qsg_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_11 * osg_460[k]
                   + f_3 * pc_z[k] * qsg_595[k];

        t_836[k] = f_20 * osg_477[k]
                   + f_6 * qsf0_398[k]
                   - f_7 * qsf1_398[k]
                   + f_3 * pc_y[k] * qsg_597[k];

        t_837[k] = f_20 * osg_478[k]
                   + f_4 * qsf0_399[k]
                   - f_5 * qsf1_399[k]
                   + f_3 * pc_y[k] * qsg_598[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pc_x, pc_y, pc_z, osg_464, osg_479, osg_600, \
                         qsf0_399, qsf0_400, qsf1_399, qsf1_400, qsg_599, \
                         qsg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_20 * osg_479[k]
                   + f_3 * pc_y[k] * qsg_599[k];

        t_839[k] = f_11 * osg_464[k]
                   + f_1 * qsf0_399[k]
                   - f_2 * qsf1_399[k]
                   + f_3 * pc_z[k] * qsg_599[k];

        t_840[k] = f_18 * osg_600[k]
                   + f_1 * qsf0_400[k]
                   - f_2 * qsf1_400[k]
                   + f_3 * pc_x[k] * qsg_600[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pc_x, pc_y, pc_z, osg_465, osg_480, \
                         osg_482, osg_603, qsf0_403, qsf1_403, qsg_600, qsg_602, \
                         qsg_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_18 * osg_480[k]
                   + f_3 * pc_y[k] * qsg_600[k];

        t_842[k] = f_18 * osg_465[k]
                   + f_3 * pc_z[k] * qsg_600[k];

        t_843[k] = f_18 * osg_603[k]
                   + f_6 * qsf0_403[k]
                   - f_7 * qsf1_403[k]
                   + f_3 * pc_x[k] * qsg_603[k];

        t_844[k] = f_18 * osg_482[k]
                   + f_3 * pc_y[k] * qsg_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_z, osg_468, osg_605, osg_606, qsf0_405, \
                         qsf0_406, qsf1_405, qsf1_406, qsg_603, qsg_605, \
                         qsg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_18 * osg_605[k]
                   + f_6 * qsf0_405[k]
                   - f_7 * qsf1_405[k]
                   + f_3 * pc_x[k] * qsg_605[k];

        t_846[k] = f_18 * osg_606[k]
                   + f_4 * qsf0_406[k]
                   - f_5 * qsf1_406[k]
                   + f_3 * pc_x[k] * qsg_606[k];

        t_847[k] = f_18 * osg_468[k]
                   + f_3 * pc_z[k] * qsg_603[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, pc_y, osg_485, osg_609, osg_610, \
                         osg_611, qsf0_409, qsf1_409, qsg_605, qsg_609, qsg_610, \
                         qsg_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_18 * osg_485[k]
                   + f_3 * pc_y[k] * qsg_605[k];

        t_849[k] = f_18 * osg_609[k]
                   + f_4 * qsf0_409[k]
                   - f_5 * qsf1_409[k]
                   + f_3 * pc_x[k] * qsg_609[k];

        t_850[k] = f_18 * osg_610[k]
                   + f_3 * pc_x[k] * qsg_610[k];

        t_851[k] = f_18 * osg_611[k]
                   + f_3 * pc_x[k] * qsg_611[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, pc_y, osg_490, osg_612, osg_613, \
                         osg_614, qsf0_406, qsf1_406, qsg_610, qsg_612, qsg_613, \
                         qsg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_18 * osg_612[k]
                   + f_3 * pc_x[k] * qsg_612[k];

        t_853[k] = f_18 * osg_613[k]
                   + f_3 * pc_x[k] * qsg_613[k];

        t_854[k] = f_18 * osg_614[k]
                   + f_3 * pc_x[k] * qsg_614[k];

        t_855[k] = f_18 * osg_490[k]
                   + f_1 * qsf0_406[k]
                   - f_2 * qsf1_406[k]
                   + f_3 * pc_y[k] * qsg_610[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, osg_475, osg_492, osg_493, qsf0_408, \
                         qsf0_409, qsf1_408, qsf1_409, qsg_610, qsg_612, \
                         qsg_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_18 * osg_475[k]
                   + f_3 * pc_z[k] * qsg_610[k];

        t_857[k] = f_18 * osg_492[k]
                   + f_6 * qsf0_408[k]
                   - f_7 * qsf1_408[k]
                   + f_3 * pc_y[k] * qsg_612[k];

        t_858[k] = f_18 * osg_493[k]
                   + f_4 * qsf0_409[k]
                   - f_5 * qsf1_409[k]
                   + f_3 * pc_y[k] * qsg_613[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_x, pc_y, pc_z, osg_479, osg_494, osg_615, \
                         qsf0_409, qsf0_410, qsf1_409, qsf1_410, qsg_614, \
                         qsg_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_18 * osg_494[k]
                   + f_3 * pc_y[k] * qsg_614[k];

        t_860[k] = f_18 * osg_479[k]
                   + f_1 * qsf0_409[k]
                   - f_2 * qsf1_409[k]
                   + f_3 * pc_z[k] * qsg_614[k];

        t_861[k] = f_18 * osg_615[k]
                   + f_1 * qsf0_410[k]
                   - f_2 * qsf1_410[k]
                   + f_3 * pc_x[k] * qsg_615[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pc_x, pc_y, pc_z, osg_480, osg_495, \
                         osg_497, osg_618, qsf0_413, qsf1_413, qsg_615, qsg_617, \
                         qsg_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_11 * osg_495[k]
                   + f_3 * pc_y[k] * qsg_615[k];

        t_863[k] = f_20 * osg_480[k]
                   + f_3 * pc_z[k] * qsg_615[k];

        t_864[k] = f_18 * osg_618[k]
                   + f_6 * qsf0_413[k]
                   - f_7 * qsf1_413[k]
                   + f_3 * pc_x[k] * qsg_618[k];

        t_865[k] = f_11 * osg_497[k]
                   + f_3 * pc_y[k] * qsg_617[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_z, osg_483, osg_620, osg_621, qsf0_415, \
                         qsf0_416, qsf1_415, qsf1_416, qsg_618, qsg_620, \
                         qsg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * osg_620[k]
                   + f_6 * qsf0_415[k]
                   - f_7 * qsf1_415[k]
                   + f_3 * pc_x[k] * qsg_620[k];

        t_867[k] = f_18 * osg_621[k]
                   + f_4 * qsf0_416[k]
                   - f_5 * qsf1_416[k]
                   + f_3 * pc_x[k] * qsg_621[k];

        t_868[k] = f_20 * osg_483[k]
                   + f_3 * pc_z[k] * qsg_618[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, osg_500, osg_624, osg_625, \
                         osg_626, qsf0_419, qsf1_419, qsg_620, qsg_624, qsg_625, \
                         qsg_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_11 * osg_500[k]
                   + f_3 * pc_y[k] * qsg_620[k];

        t_870[k] = f_18 * osg_624[k]
                   + f_4 * qsf0_419[k]
                   - f_5 * qsf1_419[k]
                   + f_3 * pc_x[k] * qsg_624[k];

        t_871[k] = f_18 * osg_625[k]
                   + f_3 * pc_x[k] * qsg_625[k];

        t_872[k] = f_18 * osg_626[k]
                   + f_3 * pc_x[k] * qsg_626[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pc_x, pc_y, osg_505, osg_627, osg_628, \
                         osg_629, qsf0_416, qsf1_416, qsg_625, qsg_627, qsg_628, \
                         qsg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_18 * osg_627[k]
                   + f_3 * pc_x[k] * qsg_627[k];

        t_874[k] = f_18 * osg_628[k]
                   + f_3 * pc_x[k] * qsg_628[k];

        t_875[k] = f_18 * osg_629[k]
                   + f_3 * pc_x[k] * qsg_629[k];

        t_876[k] = f_11 * osg_505[k]
                   + f_1 * qsf0_416[k]
                   - f_2 * qsf1_416[k]
                   + f_3 * pc_y[k] * qsg_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, osg_490, osg_507, osg_508, qsf0_418, \
                         qsf0_419, qsf1_418, qsf1_419, qsg_625, qsg_627, \
                         qsg_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_20 * osg_490[k]
                   + f_3 * pc_z[k] * qsg_625[k];

        t_878[k] = f_11 * osg_507[k]
                   + f_6 * qsf0_418[k]
                   - f_7 * qsf1_418[k]
                   + f_3 * pc_y[k] * qsg_627[k];

        t_879[k] = f_11 * osg_508[k]
                   + f_4 * qsf0_419[k]
                   - f_5 * qsf1_419[k]
                   + f_3 * pc_y[k] * qsg_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_x, pc_y, pc_z, osg_494, osg_509, osg_630, \
                         qsf0_419, qsf0_420, qsf1_419, qsf1_420, qsg_629, \
                         qsg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * osg_509[k]
                   + f_3 * pc_y[k] * qsg_629[k];

        t_881[k] = f_20 * osg_494[k]
                   + f_1 * qsf0_419[k]
                   - f_2 * qsf1_419[k]
                   + f_3 * pc_z[k] * qsg_629[k];

        t_882[k] = f_18 * osg_630[k]
                   + f_1 * qsf0_420[k]
                   - f_2 * qsf1_420[k]
                   + f_3 * pc_x[k] * qsg_630[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pc_x, pc_y, pc_z, osg_495, osg_510, \
                         osg_512, osg_633, qsf0_423, qsf1_423, qsg_630, qsg_632, \
                         qsg_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_10 * osg_510[k]
                   + f_3 * pc_y[k] * qsg_630[k];

        t_884[k] = f_21 * osg_495[k]
                   + f_3 * pc_z[k] * qsg_630[k];

        t_885[k] = f_18 * osg_633[k]
                   + f_6 * qsf0_423[k]
                   - f_7 * qsf1_423[k]
                   + f_3 * pc_x[k] * qsg_633[k];

        t_886[k] = f_10 * osg_512[k]
                   + f_3 * pc_y[k] * qsg_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pc_x, pc_z, osg_498, osg_635, osg_636, qsf0_425, \
                         qsf0_426, qsf1_425, qsf1_426, qsg_633, qsg_635, \
                         qsg_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_18 * osg_635[k]
                   + f_6 * qsf0_425[k]
                   - f_7 * qsf1_425[k]
                   + f_3 * pc_x[k] * qsg_635[k];

        t_888[k] = f_18 * osg_636[k]
                   + f_4 * qsf0_426[k]
                   - f_5 * qsf1_426[k]
                   + f_3 * pc_x[k] * qsg_636[k];

        t_889[k] = f_21 * osg_498[k]
                   + f_3 * pc_z[k] * qsg_633[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pc_x, pc_y, osg_515, osg_639, osg_640, \
                         osg_641, qsf0_429, qsf1_429, qsg_635, qsg_639, qsg_640, \
                         qsg_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_10 * osg_515[k]
                   + f_3 * pc_y[k] * qsg_635[k];

        t_891[k] = f_18 * osg_639[k]
                   + f_4 * qsf0_429[k]
                   - f_5 * qsf1_429[k]
                   + f_3 * pc_x[k] * qsg_639[k];

        t_892[k] = f_18 * osg_640[k]
                   + f_3 * pc_x[k] * qsg_640[k];

        t_893[k] = f_18 * osg_641[k]
                   + f_3 * pc_x[k] * qsg_641[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pc_x, pc_y, osg_520, osg_642, osg_643, \
                         osg_644, qsf0_426, qsf1_426, qsg_640, qsg_642, qsg_643, \
                         qsg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_18 * osg_642[k]
                   + f_3 * pc_x[k] * qsg_642[k];

        t_895[k] = f_18 * osg_643[k]
                   + f_3 * pc_x[k] * qsg_643[k];

        t_896[k] = f_18 * osg_644[k]
                   + f_3 * pc_x[k] * qsg_644[k];

        t_897[k] = f_10 * osg_520[k]
                   + f_1 * qsf0_426[k]
                   - f_2 * qsf1_426[k]
                   + f_3 * pc_y[k] * qsg_640[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_y, pc_z, osg_505, osg_522, osg_523, qsf0_428, \
                         qsf0_429, qsf1_428, qsf1_429, qsg_640, qsg_642, \
                         qsg_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_21 * osg_505[k]
                   + f_3 * pc_z[k] * qsg_640[k];

        t_899[k] = f_10 * osg_522[k]
                   + f_6 * qsf0_428[k]
                   - f_7 * qsf1_428[k]
                   + f_3 * pc_y[k] * qsg_642[k];

        t_900[k] = f_10 * osg_523[k]
                   + f_4 * qsf0_429[k]
                   - f_5 * qsf1_429[k]
                   + f_3 * pc_y[k] * qsg_643[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pa_y, pc_y, pc_z, osh0_735, osg_509, \
                         osg_524, osg_525, osh1_735, qsf0_429, qsf1_429, qsg_644, \
                         qsg_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_10 * osg_524[k]
                   + f_3 * pc_y[k] * qsg_644[k];

        t_902[k] = f_21 * osg_509[k]
                   + f_1 * qsf0_429[k]
                   - f_2 * qsf1_429[k]
                   + f_3 * pc_z[k] * qsg_644[k];

        t_903[k] = pa_y[k] * osh0_735[k]
                   - f_8 * pc_y[k] * osh1_735[k];

        t_904[k] = f_9 * osg_525[k]
                   + f_3 * pc_y[k] * qsg_645[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pa_y, pc_y, pc_z, osh0_738, osh0_740, \
                         osg_510, osg_526, osg_527, osh1_738, osh1_740, qsg_645, \
                         qsg_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_19 * osg_510[k]
                   + f_3 * pc_z[k] * qsg_645[k];

        t_906[k] = pa_y[k] * osh0_738[k]
                   + f_10 * osg_526[k]
                   - f_8 * pc_y[k] * osh1_738[k];

        t_907[k] = f_9 * osg_527[k]
                   + f_3 * pc_y[k] * qsg_647[k];

        t_908[k] = pa_y[k] * osh0_740[k]
                   - f_8 * pc_y[k] * osh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pa_y, pc_y, pc_z, osh0_741, osh0_744, \
                         osg_513, osg_528, osg_530, osh1_741, osh1_744, qsg_648, \
                         qsg_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pa_y[k] * osh0_741[k]
                   + f_11 * osg_528[k]
                   - f_8 * pc_y[k] * osh1_741[k];

        t_910[k] = f_19 * osg_513[k]
                   + f_3 * pc_z[k] * qsg_648[k];

        t_911[k] = f_9 * osg_530[k]
                   + f_3 * pc_y[k] * qsg_650[k];

        t_912[k] = pa_y[k] * osh0_744[k]
                   - f_8 * pc_y[k] * osh1_744[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, pc_x, osg_655, osg_656, osg_657, \
                         osg_658, osg_659, qsg_655, qsg_656, qsg_657, qsg_658, \
                         qsg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_18 * osg_655[k]
                   + f_3 * pc_x[k] * qsg_655[k];

        t_914[k] = f_18 * osg_656[k]
                   + f_3 * pc_x[k] * qsg_656[k];

        t_915[k] = f_18 * osg_657[k]
                   + f_3 * pc_x[k] * qsg_657[k];

        t_916[k] = f_18 * osg_658[k]
                   + f_3 * pc_x[k] * qsg_658[k];

        t_917[k] = f_18 * osg_659[k]
                   + f_3 * pc_x[k] * qsg_659[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_y, pc_z, osg_520, osg_535, osg_537, qsf0_436, \
                         qsf0_438, qsf1_436, qsf1_438, qsg_655, \
                         qsg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_9 * osg_535[k]
                   + f_1 * qsf0_436[k]
                   - f_2 * qsf1_436[k]
                   + f_3 * pc_y[k] * qsg_655[k];

        t_919[k] = f_19 * osg_520[k]
                   + f_3 * pc_z[k] * qsg_655[k];

        t_920[k] = f_9 * osg_537[k]
                   + f_6 * qsf0_438[k]
                   - f_7 * qsf1_438[k]
                   + f_3 * pc_y[k] * qsg_657[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pa_y, pc_y, osh0_755, osg_538, osg_539, \
                         osh1_755, qsf0_439, qsf1_439, qsg_658, \
                         qsg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_9 * osg_538[k]
                   + f_4 * qsf0_439[k]
                   - f_5 * qsf1_439[k]
                   + f_3 * pc_y[k] * qsg_658[k];

        t_922[k] = f_9 * osg_539[k]
                   + f_3 * pc_y[k] * qsg_659[k];

        t_923[k] = pa_y[k] * osh0_755[k]
                   - f_8 * pc_y[k] * osh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, osg_525, \
                         osg_660, qsf0_440, qsf1_440, qsg_660, qsg_661, \
                         qsg_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_18 * osg_660[k]
                   + f_1 * qsf0_440[k]
                   - f_2 * qsf1_440[k]
                   + f_3 * pc_x[k] * qsg_660[k];

        t_925[k] = f_3 * pc_y[k] * qsg_660[k];

        t_926[k] = f_17 * osg_525[k]
                   + f_3 * pc_z[k] * qsg_660[k];

        t_927[k] = f_4 * qsf0_440[k]
                   - f_5 * qsf1_440[k]
                   + f_3 * pc_y[k] * qsg_661[k];

        t_928[k] = f_3 * pc_y[k] * qsg_662[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, t_932, pc_x, pc_y, osg_665, qsf0_441, qsf0_442, \
                         qsf0_445, qsf1_441, qsf1_442, qsf1_445, qsg_663, qsg_664, \
                         qsg_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_18 * osg_665[k]
                   + f_6 * qsf0_445[k]
                   - f_7 * qsf1_445[k]
                   + f_3 * pc_x[k] * qsg_665[k];

        t_930[k] = f_6 * qsf0_441[k]
                   - f_7 * qsf1_441[k]
                   + f_3 * pc_y[k] * qsg_663[k];

        t_931[k] = f_4 * qsf0_442[k]
                   - f_5 * qsf1_442[k]
                   + f_3 * pc_y[k] * qsg_664[k];

        t_932[k] = f_3 * pc_y[k] * qsg_665[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, pc_x, osg_669, osg_670, osg_671, osg_672, \
                         qsf0_449, qsf1_449, qsg_669, qsg_670, qsg_671, \
                         qsg_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_18 * osg_669[k]
                   + f_4 * qsf0_449[k]
                   - f_5 * qsf1_449[k]
                   + f_3 * pc_x[k] * qsg_669[k];

        t_934[k] = f_18 * osg_670[k]
                   + f_3 * pc_x[k] * qsg_670[k];

        t_935[k] = f_18 * osg_671[k]
                   + f_3 * pc_x[k] * qsg_671[k];

        t_936[k] = f_18 * osg_672[k]
                   + f_3 * pc_x[k] * qsg_672[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, pc_x, pc_y, osg_674, qsf0_446, qsf0_447, \
                         qsf1_446, qsf1_447, qsg_669, qsg_670, qsg_671, \
                         qsg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = f_3 * pc_y[k] * qsg_669[k];

        t_938[k] = f_18 * osg_674[k]
                   + f_3 * pc_x[k] * qsg_674[k];

        t_939[k] = f_1 * qsf0_446[k]
                   - f_2 * qsf1_446[k]
                   + f_3 * pc_y[k] * qsg_670[k];

        t_940[k] = f_13 * qsf0_447[k]
                   - f_14 * qsf1_447[k]
                   + f_3 * pc_y[k] * qsg_671[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, t_944, pc_y, pc_z, osg_539, qsf0_448, qsf0_449, \
                         qsf1_448, qsf1_449, qsg_672, qsg_673, \
                         qsg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = f_6 * qsf0_448[k]
                   - f_7 * qsf1_448[k]
                   + f_3 * pc_y[k] * qsg_672[k];

        t_942[k] = f_4 * qsf0_449[k]
                   - f_5 * qsf1_449[k]
                   + f_3 * pc_y[k] * qsg_673[k];

        t_943[k] = f_3 * pc_y[k] * qsg_674[k];

        t_944[k] = f_17 * osg_539[k]
                   + f_1 * qsf0_449[k]
                   - f_2 * qsf1_449[k]
                   + f_3 * pc_z[k] * qsg_674[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *osh0_756 = buffer.data(osh0 + 756);
    const auto *osh0_759 = buffer.data(osh0 + 759);
    const auto *osh0_762 = buffer.data(osh0 + 762);
    const auto *osh0_771 = buffer.data(osh0 + 771);

    const auto *osg_540 = buffer.data(osg + 540);
    const auto *osg_543 = buffer.data(osg + 543);
    const auto *osg_545 = buffer.data(osg + 545);
    const auto *osg_550 = buffer.data(osg + 550);
    const auto *osg_554 = buffer.data(osg + 554);
    const auto *osg_555 = buffer.data(osg + 555);
    const auto *osg_557 = buffer.data(osg + 557);
    const auto *osg_558 = buffer.data(osg + 558);
    const auto *osg_560 = buffer.data(osg + 560);
    const auto *osg_565 = buffer.data(osg + 565);
    const auto *osg_567 = buffer.data(osg + 567);
    const auto *osg_568 = buffer.data(osg + 568);
    const auto *osg_569 = buffer.data(osg + 569);
    const auto *osg_570 = buffer.data(osg + 570);
    const auto *osg_572 = buffer.data(osg + 572);
    const auto *osg_573 = buffer.data(osg + 573);
    const auto *osg_575 = buffer.data(osg + 575);
    const auto *osg_580 = buffer.data(osg + 580);
    const auto *osg_582 = buffer.data(osg + 582);
    const auto *osg_583 = buffer.data(osg + 583);
    const auto *osg_584 = buffer.data(osg + 584);
    const auto *osg_585 = buffer.data(osg + 585);
    const auto *osg_587 = buffer.data(osg + 587);
    const auto *osg_588 = buffer.data(osg + 588);
    const auto *osg_590 = buffer.data(osg + 590);
    const auto *osg_595 = buffer.data(osg + 595);
    const auto *osg_597 = buffer.data(osg + 597);
    const auto *osg_598 = buffer.data(osg + 598);
    const auto *osg_599 = buffer.data(osg + 599);
    const auto *osg_600 = buffer.data(osg + 600);
    const auto *osg_602 = buffer.data(osg + 602);
    const auto *osg_603 = buffer.data(osg + 603);
    const auto *osg_605 = buffer.data(osg + 605);
    const auto *osg_610 = buffer.data(osg + 610);
    const auto *osg_612 = buffer.data(osg + 612);
    const auto *osg_613 = buffer.data(osg + 613);
    const auto *osg_614 = buffer.data(osg + 614);
    const auto *osg_615 = buffer.data(osg + 615);
    const auto *osg_617 = buffer.data(osg + 617);
    const auto *osg_675 = buffer.data(osg + 675);
    const auto *osg_678 = buffer.data(osg + 678);
    const auto *osg_681 = buffer.data(osg + 681);
    const auto *osg_685 = buffer.data(osg + 685);
    const auto *osg_687 = buffer.data(osg + 687);
    const auto *osg_688 = buffer.data(osg + 688);
    const auto *osg_689 = buffer.data(osg + 689);
    const auto *osg_695 = buffer.data(osg + 695);
    const auto *osg_699 = buffer.data(osg + 699);
    const auto *osg_700 = buffer.data(osg + 700);
    const auto *osg_701 = buffer.data(osg + 701);
    const auto *osg_702 = buffer.data(osg + 702);
    const auto *osg_703 = buffer.data(osg + 703);
    const auto *osg_704 = buffer.data(osg + 704);
    const auto *osg_705 = buffer.data(osg + 705);
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
    const auto *osg_753 = buffer.data(osg + 753);
    const auto *osg_755 = buffer.data(osg + 755);
    const auto *osg_756 = buffer.data(osg + 756);

    const auto *osh1_756 = buffer.data(osh1 + 756);
    const auto *osh1_759 = buffer.data(osh1 + 759);
    const auto *osh1_762 = buffer.data(osh1 + 762);
    const auto *osh1_771 = buffer.data(osh1 + 771);

    const auto *qsf0_450 = buffer.data(qsf0 + 450);
    const auto *qsf0_452 = buffer.data(qsf0 + 452);
    const auto *qsf0_453 = buffer.data(qsf0 + 453);
    const auto *qsf0_456 = buffer.data(qsf0 + 456);
    const auto *qsf0_457 = buffer.data(qsf0 + 457);
    const auto *qsf0_459 = buffer.data(qsf0 + 459);
    const auto *qsf0_465 = buffer.data(qsf0 + 465);
    const auto *qsf0_468 = buffer.data(qsf0 + 468);
    const auto *qsf0_469 = buffer.data(qsf0 + 469);
    const auto *qsf0_470 = buffer.data(qsf0 + 470);
    const auto *qsf0_473 = buffer.data(qsf0 + 473);
    const auto *qsf0_475 = buffer.data(qsf0 + 475);
    const auto *qsf0_476 = buffer.data(qsf0 + 476);
    const auto *qsf0_478 = buffer.data(qsf0 + 478);
    const auto *qsf0_479 = buffer.data(qsf0 + 479);
    const auto *qsf0_480 = buffer.data(qsf0 + 480);
    const auto *qsf0_483 = buffer.data(qsf0 + 483);
    const auto *qsf0_485 = buffer.data(qsf0 + 485);
    const auto *qsf0_486 = buffer.data(qsf0 + 486);
    const auto *qsf0_488 = buffer.data(qsf0 + 488);
    const auto *qsf0_489 = buffer.data(qsf0 + 489);
    const auto *qsf0_490 = buffer.data(qsf0 + 490);
    const auto *qsf0_493 = buffer.data(qsf0 + 493);
    const auto *qsf0_495 = buffer.data(qsf0 + 495);
    const auto *qsf0_496 = buffer.data(qsf0 + 496);
    const auto *qsf0_498 = buffer.data(qsf0 + 498);
    const auto *qsf0_499 = buffer.data(qsf0 + 499);
    const auto *qsf0_500 = buffer.data(qsf0 + 500);
    const auto *qsf0_503 = buffer.data(qsf0 + 503);
    const auto *qsf0_505 = buffer.data(qsf0 + 505);
    const auto *qsf0_506 = buffer.data(qsf0 + 506);

    const auto *qsf1_450 = buffer.data(qsf1 + 450);
    const auto *qsf1_452 = buffer.data(qsf1 + 452);
    const auto *qsf1_453 = buffer.data(qsf1 + 453);
    const auto *qsf1_456 = buffer.data(qsf1 + 456);
    const auto *qsf1_457 = buffer.data(qsf1 + 457);
    const auto *qsf1_459 = buffer.data(qsf1 + 459);
    const auto *qsf1_465 = buffer.data(qsf1 + 465);
    const auto *qsf1_468 = buffer.data(qsf1 + 468);
    const auto *qsf1_469 = buffer.data(qsf1 + 469);
    const auto *qsf1_470 = buffer.data(qsf1 + 470);
    const auto *qsf1_473 = buffer.data(qsf1 + 473);
    const auto *qsf1_475 = buffer.data(qsf1 + 475);
    const auto *qsf1_476 = buffer.data(qsf1 + 476);
    const auto *qsf1_478 = buffer.data(qsf1 + 478);
    const auto *qsf1_479 = buffer.data(qsf1 + 479);
    const auto *qsf1_480 = buffer.data(qsf1 + 480);
    const auto *qsf1_483 = buffer.data(qsf1 + 483);
    const auto *qsf1_485 = buffer.data(qsf1 + 485);
    const auto *qsf1_486 = buffer.data(qsf1 + 486);
    const auto *qsf1_488 = buffer.data(qsf1 + 488);
    const auto *qsf1_489 = buffer.data(qsf1 + 489);
    const auto *qsf1_490 = buffer.data(qsf1 + 490);
    const auto *qsf1_493 = buffer.data(qsf1 + 493);
    const auto *qsf1_495 = buffer.data(qsf1 + 495);
    const auto *qsf1_496 = buffer.data(qsf1 + 496);
    const auto *qsf1_498 = buffer.data(qsf1 + 498);
    const auto *qsf1_499 = buffer.data(qsf1 + 499);
    const auto *qsf1_500 = buffer.data(qsf1 + 500);
    const auto *qsf1_503 = buffer.data(qsf1 + 503);
    const auto *qsf1_505 = buffer.data(qsf1 + 505);
    const auto *qsf1_506 = buffer.data(qsf1 + 506);

    const auto *qsg_675 = buffer.data(qsg + 675);
    const auto *qsg_676 = buffer.data(qsg + 676);
    const auto *qsg_677 = buffer.data(qsg + 677);
    const auto *qsg_678 = buffer.data(qsg + 678);
    const auto *qsg_680 = buffer.data(qsg + 680);
    const auto *qsg_681 = buffer.data(qsg + 681);
    const auto *qsg_685 = buffer.data(qsg + 685);
    const auto *qsg_686 = buffer.data(qsg + 686);
    const auto *qsg_687 = buffer.data(qsg + 687);
    const auto *qsg_688 = buffer.data(qsg + 688);
    const auto *qsg_689 = buffer.data(qsg + 689);
    const auto *qsg_690 = buffer.data(qsg + 690);
    const auto *qsg_692 = buffer.data(qsg + 692);
    const auto *qsg_693 = buffer.data(qsg + 693);
    const auto *qsg_695 = buffer.data(qsg + 695);
    const auto *qsg_699 = buffer.data(qsg + 699);
    const auto *qsg_700 = buffer.data(qsg + 700);
    const auto *qsg_701 = buffer.data(qsg + 701);
    const auto *qsg_702 = buffer.data(qsg + 702);
    const auto *qsg_703 = buffer.data(qsg + 703);
    const auto *qsg_704 = buffer.data(qsg + 704);
    const auto *qsg_705 = buffer.data(qsg + 705);
    const auto *qsg_707 = buffer.data(qsg + 707);
    const auto *qsg_708 = buffer.data(qsg + 708);
    const auto *qsg_710 = buffer.data(qsg + 710);
    const auto *qsg_711 = buffer.data(qsg + 711);
    const auto *qsg_714 = buffer.data(qsg + 714);
    const auto *qsg_715 = buffer.data(qsg + 715);
    const auto *qsg_716 = buffer.data(qsg + 716);
    const auto *qsg_717 = buffer.data(qsg + 717);
    const auto *qsg_718 = buffer.data(qsg + 718);
    const auto *qsg_719 = buffer.data(qsg + 719);
    const auto *qsg_720 = buffer.data(qsg + 720);
    const auto *qsg_722 = buffer.data(qsg + 722);
    const auto *qsg_723 = buffer.data(qsg + 723);
    const auto *qsg_725 = buffer.data(qsg + 725);
    const auto *qsg_726 = buffer.data(qsg + 726);
    const auto *qsg_729 = buffer.data(qsg + 729);
    const auto *qsg_730 = buffer.data(qsg + 730);
    const auto *qsg_731 = buffer.data(qsg + 731);
    const auto *qsg_732 = buffer.data(qsg + 732);
    const auto *qsg_733 = buffer.data(qsg + 733);
    const auto *qsg_734 = buffer.data(qsg + 734);
    const auto *qsg_735 = buffer.data(qsg + 735);
    const auto *qsg_737 = buffer.data(qsg + 737);
    const auto *qsg_738 = buffer.data(qsg + 738);
    const auto *qsg_740 = buffer.data(qsg + 740);
    const auto *qsg_741 = buffer.data(qsg + 741);
    const auto *qsg_744 = buffer.data(qsg + 744);
    const auto *qsg_745 = buffer.data(qsg + 745);
    const auto *qsg_746 = buffer.data(qsg + 746);
    const auto *qsg_747 = buffer.data(qsg + 747);
    const auto *qsg_748 = buffer.data(qsg + 748);
    const auto *qsg_749 = buffer.data(qsg + 749);
    const auto *qsg_750 = buffer.data(qsg + 750);
    const auto *qsg_752 = buffer.data(qsg + 752);
    const auto *qsg_753 = buffer.data(qsg + 753);
    const auto *qsg_755 = buffer.data(qsg + 755);
    const auto *qsg_756 = buffer.data(qsg + 756);

#pragma omp simd aligned(t_945, t_946, t_947, t_948, pc_x, pc_y, pc_z, osg_540, osg_675, \
                         osg_678, qsf0_450, qsf0_453, qsf1_450, qsf1_453, qsg_675, \
                         qsg_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_11 * osg_675[k]
                   + f_1 * qsf0_450[k]
                   - f_2 * qsf1_450[k]
                   + f_3 * pc_x[k] * qsg_675[k];

        t_946[k] = f_16 * osg_540[k]
                   + f_3 * pc_y[k] * qsg_675[k];

        t_947[k] = f_3 * pc_z[k] * qsg_675[k];

        t_948[k] = f_11 * osg_678[k]
                   + f_6 * qsf0_453[k]
                   - f_7 * qsf1_453[k]
                   + f_3 * pc_x[k] * qsg_678[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pc_x, pc_z, osg_681, qsf0_450, qsf0_456, \
                         qsf1_450, qsf1_456, qsg_676, qsg_677, qsg_678, \
                         qsg_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_3 * pc_z[k] * qsg_676[k];

        t_950[k] = f_4 * qsf0_450[k]
                   - f_5 * qsf1_450[k]
                   + f_3 * pc_z[k] * qsg_677[k];

        t_951[k] = f_11 * osg_681[k]
                   + f_4 * qsf0_456[k]
                   - f_5 * qsf1_456[k]
                   + f_3 * pc_x[k] * qsg_681[k];

        t_952[k] = f_3 * pc_z[k] * qsg_678[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pc_x, pc_y, pc_z, osg_545, osg_685, \
                         qsf0_452, qsf1_452, qsg_680, qsg_681, \
                         qsg_685 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_16 * osg_545[k]
                   + f_3 * pc_y[k] * qsg_680[k];

        t_954[k] = f_6 * qsf0_452[k]
                   - f_7 * qsf1_452[k]
                   + f_3 * pc_z[k] * qsg_680[k];

        t_955[k] = f_11 * osg_685[k]
                   + f_3 * pc_x[k] * qsg_685[k];

        t_956[k] = f_3 * pc_z[k] * qsg_681[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pc_x, pc_y, osg_550, osg_687, osg_688, \
                         osg_689, qsf0_456, qsf1_456, qsg_685, qsg_687, qsg_688, \
                         qsg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_11 * osg_687[k]
                   + f_3 * pc_x[k] * qsg_687[k];

        t_958[k] = f_11 * osg_688[k]
                   + f_3 * pc_x[k] * qsg_688[k];

        t_959[k] = f_11 * osg_689[k]
                   + f_3 * pc_x[k] * qsg_689[k];

        t_960[k] = f_16 * osg_550[k]
                   + f_1 * qsf0_456[k]
                   - f_2 * qsf1_456[k]
                   + f_3 * pc_y[k] * qsg_685[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pc_y, pc_z, osg_554, qsf0_456, qsf0_457, \
                         qsf1_456, qsf1_457, qsg_685, qsg_686, qsg_687, \
                         qsg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * qsg_685[k];

        t_962[k] = f_4 * qsf0_456[k]
                   - f_5 * qsf1_456[k]
                   + f_3 * pc_z[k] * qsg_686[k];

        t_963[k] = f_6 * qsf0_457[k]
                   - f_7 * qsf1_457[k]
                   + f_3 * pc_z[k] * qsg_687[k];

        t_964[k] = f_16 * osg_554[k]
                   + f_3 * pc_y[k] * qsg_689[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pa_z, pc_y, pc_z, osh0_756, osg_540, \
                         osg_555, osh1_756, qsf0_459, qsf1_459, qsg_689, \
                         qsg_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_1 * qsf0_459[k]
                   - f_2 * qsf1_459[k]
                   + f_3 * pc_z[k] * qsg_689[k];

        t_966[k] = pa_z[k] * osh0_756[k]
                   - f_8 * pc_z[k] * osh1_756[k];

        t_967[k] = f_17 * osg_555[k]
                   + f_3 * pc_y[k] * qsg_690[k];

        t_968[k] = f_9 * osg_540[k]
                   + f_3 * pc_z[k] * qsg_690[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_z, pc_x, pc_y, pc_z, osh0_759, osg_557, \
                         osg_695, osh1_759, qsf0_465, qsf1_465, qsg_692, \
                         qsg_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = pa_z[k] * osh0_759[k]
                   - f_8 * pc_z[k] * osh1_759[k];

        t_970[k] = f_17 * osg_557[k]
                   + f_3 * pc_y[k] * qsg_692[k];

        t_971[k] = f_11 * osg_695[k]
                   + f_6 * qsf0_465[k]
                   - f_7 * qsf1_465[k]
                   + f_3 * pc_x[k] * qsg_695[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, pa_z, pc_y, pc_z, osh0_762, osg_543, osg_560, \
                         osh1_762, qsg_693, qsg_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = pa_z[k] * osh0_762[k]
                   - f_8 * pc_z[k] * osh1_762[k];

        t_973[k] = f_9 * osg_543[k]
                   + f_3 * pc_z[k] * qsg_693[k];

        t_974[k] = f_17 * osg_560[k]
                   + f_3 * pc_y[k] * qsg_695[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, pc_x, osg_699, osg_700, osg_701, osg_702, \
                         qsf0_469, qsf1_469, qsg_699, qsg_700, qsg_701, \
                         qsg_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_11 * osg_699[k]
                   + f_4 * qsf0_469[k]
                   - f_5 * qsf1_469[k]
                   + f_3 * pc_x[k] * qsg_699[k];

        t_976[k] = f_11 * osg_700[k]
                   + f_3 * pc_x[k] * qsg_700[k];

        t_977[k] = f_11 * osg_701[k]
                   + f_3 * pc_x[k] * qsg_701[k];

        t_978[k] = f_11 * osg_702[k]
                   + f_3 * pc_x[k] * qsg_702[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_z, pc_x, pc_z, osh0_771, osg_550, \
                         osg_703, osg_704, osh1_771, qsg_700, qsg_703, \
                         qsg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_11 * osg_703[k]
                   + f_3 * pc_x[k] * qsg_703[k];

        t_980[k] = f_11 * osg_704[k]
                   + f_3 * pc_x[k] * qsg_704[k];

        t_981[k] = pa_z[k] * osh0_771[k]
                   - f_8 * pc_z[k] * osh1_771[k];

        t_982[k] = f_9 * osg_550[k]
                   + f_3 * pc_z[k] * qsg_700[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_y, osg_567, osg_568, osg_569, qsf0_468, \
                         qsf0_469, qsf1_468, qsf1_469, qsg_702, qsg_703, \
                         qsg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_17 * osg_567[k]
                   + f_6 * qsf0_468[k]
                   - f_7 * qsf1_468[k]
                   + f_3 * pc_y[k] * qsg_702[k];

        t_984[k] = f_17 * osg_568[k]
                   + f_4 * qsf0_469[k]
                   - f_5 * qsf1_469[k]
                   + f_3 * pc_y[k] * qsg_703[k];

        t_985[k] = f_17 * osg_569[k]
                   + f_3 * pc_y[k] * qsg_704[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_x, pc_y, pc_z, osg_554, osg_570, osg_705, \
                         qsf0_469, qsf0_470, qsf1_469, qsf1_470, qsg_704, \
                         qsg_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_9 * osg_554[k]
                   + f_1 * qsf0_469[k]
                   - f_2 * qsf1_469[k]
                   + f_3 * pc_z[k] * qsg_704[k];

        t_987[k] = f_11 * osg_705[k]
                   + f_1 * qsf0_470[k]
                   - f_2 * qsf1_470[k]
                   + f_3 * pc_x[k] * qsg_705[k];

        t_988[k] = f_19 * osg_570[k]
                   + f_3 * pc_y[k] * qsg_705[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, pc_z, osg_555, osg_572, osg_708, \
                         qsf0_473, qsf1_473, qsg_705, qsg_707, \
                         qsg_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_10 * osg_555[k]
                   + f_3 * pc_z[k] * qsg_705[k];

        t_990[k] = f_11 * osg_708[k]
                   + f_6 * qsf0_473[k]
                   - f_7 * qsf1_473[k]
                   + f_3 * pc_x[k] * qsg_708[k];

        t_991[k] = f_19 * osg_572[k]
                   + f_3 * pc_y[k] * qsg_707[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, pc_z, osg_558, osg_710, osg_711, qsf0_475, \
                         qsf0_476, qsf1_475, qsf1_476, qsg_708, qsg_710, \
                         qsg_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_11 * osg_710[k]
                   + f_6 * qsf0_475[k]
                   - f_7 * qsf1_475[k]
                   + f_3 * pc_x[k] * qsg_710[k];

        t_993[k] = f_11 * osg_711[k]
                   + f_4 * qsf0_476[k]
                   - f_5 * qsf1_476[k]
                   + f_3 * pc_x[k] * qsg_711[k];

        t_994[k] = f_10 * osg_558[k]
                   + f_3 * pc_z[k] * qsg_708[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pc_x, pc_y, osg_575, osg_714, osg_715, \
                         osg_716, qsf0_479, qsf1_479, qsg_710, qsg_714, qsg_715, \
                         qsg_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_19 * osg_575[k]
                   + f_3 * pc_y[k] * qsg_710[k];

        t_996[k] = f_11 * osg_714[k]
                   + f_4 * qsf0_479[k]
                   - f_5 * qsf1_479[k]
                   + f_3 * pc_x[k] * qsg_714[k];

        t_997[k] = f_11 * osg_715[k]
                   + f_3 * pc_x[k] * qsg_715[k];

        t_998[k] = f_11 * osg_716[k]
                   + f_3 * pc_x[k] * qsg_716[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, t_1002, pc_x, pc_y, osg_580, osg_717, osg_718, \
                         osg_719, qsf0_476, qsf1_476, qsg_715, qsg_717, qsg_718, \
                         qsg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_11 * osg_717[k]
                   + f_3 * pc_x[k] * qsg_717[k];

        t_1000[k] = f_11 * osg_718[k]
                    + f_3 * pc_x[k] * qsg_718[k];

        t_1001[k] = f_11 * osg_719[k]
                    + f_3 * pc_x[k] * qsg_719[k];

        t_1002[k] = f_19 * osg_580[k]
                    + f_1 * qsf0_476[k]
                    - f_2 * qsf1_476[k]
                    + f_3 * pc_y[k] * qsg_715[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, pc_y, pc_z, osg_565, osg_582, osg_583, \
                         qsf0_478, qsf0_479, qsf1_478, qsf1_479, qsg_715, qsg_717, \
                         qsg_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_10 * osg_565[k]
                    + f_3 * pc_z[k] * qsg_715[k];

        t_1004[k] = f_19 * osg_582[k]
                    + f_6 * qsf0_478[k]
                    - f_7 * qsf1_478[k]
                    + f_3 * pc_y[k] * qsg_717[k];

        t_1005[k] = f_19 * osg_583[k]
                    + f_4 * qsf0_479[k]
                    - f_5 * qsf1_479[k]
                    + f_3 * pc_y[k] * qsg_718[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, pc_x, pc_y, pc_z, osg_569, osg_584, osg_720, \
                         qsf0_479, qsf0_480, qsf1_479, qsf1_480, qsg_719, \
                         qsg_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_19 * osg_584[k]
                    + f_3 * pc_y[k] * qsg_719[k];

        t_1007[k] = f_10 * osg_569[k]
                    + f_1 * qsf0_479[k]
                    - f_2 * qsf1_479[k]
                    + f_3 * pc_z[k] * qsg_719[k];

        t_1008[k] = f_11 * osg_720[k]
                    + f_1 * qsf0_480[k]
                    - f_2 * qsf1_480[k]
                    + f_3 * pc_x[k] * qsg_720[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pc_x, pc_y, pc_z, osg_570, osg_585, \
                         osg_587, osg_723, qsf0_483, qsf1_483, qsg_720, qsg_722, \
                         qsg_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_21 * osg_585[k]
                    + f_3 * pc_y[k] * qsg_720[k];

        t_1010[k] = f_11 * osg_570[k]
                    + f_3 * pc_z[k] * qsg_720[k];

        t_1011[k] = f_11 * osg_723[k]
                    + f_6 * qsf0_483[k]
                    - f_7 * qsf1_483[k]
                    + f_3 * pc_x[k] * qsg_723[k];

        t_1012[k] = f_21 * osg_587[k]
                    + f_3 * pc_y[k] * qsg_722[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, pc_z, osg_573, osg_725, osg_726, \
                         qsf0_485, qsf0_486, qsf1_485, qsf1_486, qsg_723, qsg_725, \
                         qsg_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_11 * osg_725[k]
                    + f_6 * qsf0_485[k]
                    - f_7 * qsf1_485[k]
                    + f_3 * pc_x[k] * qsg_725[k];

        t_1014[k] = f_11 * osg_726[k]
                    + f_4 * qsf0_486[k]
                    - f_5 * qsf1_486[k]
                    + f_3 * pc_x[k] * qsg_726[k];

        t_1015[k] = f_11 * osg_573[k]
                    + f_3 * pc_z[k] * qsg_723[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, osg_590, osg_729, \
                         osg_730, osg_731, qsf0_489, qsf1_489, qsg_725, qsg_729, qsg_730, \
                         qsg_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_21 * osg_590[k]
                    + f_3 * pc_y[k] * qsg_725[k];

        t_1017[k] = f_11 * osg_729[k]
                    + f_4 * qsf0_489[k]
                    - f_5 * qsf1_489[k]
                    + f_3 * pc_x[k] * qsg_729[k];

        t_1018[k] = f_11 * osg_730[k]
                    + f_3 * pc_x[k] * qsg_730[k];

        t_1019[k] = f_11 * osg_731[k]
                    + f_3 * pc_x[k] * qsg_731[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, osg_595, osg_732, \
                         osg_733, osg_734, qsf0_486, qsf1_486, qsg_730, qsg_732, qsg_733, \
                         qsg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_11 * osg_732[k]
                    + f_3 * pc_x[k] * qsg_732[k];

        t_1021[k] = f_11 * osg_733[k]
                    + f_3 * pc_x[k] * qsg_733[k];

        t_1022[k] = f_11 * osg_734[k]
                    + f_3 * pc_x[k] * qsg_734[k];

        t_1023[k] = f_21 * osg_595[k]
                    + f_1 * qsf0_486[k]
                    - f_2 * qsf1_486[k]
                    + f_3 * pc_y[k] * qsg_730[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pc_y, pc_z, osg_580, osg_597, osg_598, \
                         qsf0_488, qsf0_489, qsf1_488, qsf1_489, qsg_730, qsg_732, \
                         qsg_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_11 * osg_580[k]
                    + f_3 * pc_z[k] * qsg_730[k];

        t_1025[k] = f_21 * osg_597[k]
                    + f_6 * qsf0_488[k]
                    - f_7 * qsf1_488[k]
                    + f_3 * pc_y[k] * qsg_732[k];

        t_1026[k] = f_21 * osg_598[k]
                    + f_4 * qsf0_489[k]
                    - f_5 * qsf1_489[k]
                    + f_3 * pc_y[k] * qsg_733[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, pc_x, pc_y, pc_z, osg_584, osg_599, osg_735, \
                         qsf0_489, qsf0_490, qsf1_489, qsf1_490, qsg_734, \
                         qsg_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_21 * osg_599[k]
                    + f_3 * pc_y[k] * qsg_734[k];

        t_1028[k] = f_11 * osg_584[k]
                    + f_1 * qsf0_489[k]
                    - f_2 * qsf1_489[k]
                    + f_3 * pc_z[k] * qsg_734[k];

        t_1029[k] = f_11 * osg_735[k]
                    + f_1 * qsf0_490[k]
                    - f_2 * qsf1_490[k]
                    + f_3 * pc_x[k] * qsg_735[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, pc_x, pc_y, pc_z, osg_585, osg_600, \
                         osg_602, osg_738, qsf0_493, qsf1_493, qsg_735, qsg_737, \
                         qsg_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_20 * osg_600[k]
                    + f_3 * pc_y[k] * qsg_735[k];

        t_1031[k] = f_18 * osg_585[k]
                    + f_3 * pc_z[k] * qsg_735[k];

        t_1032[k] = f_11 * osg_738[k]
                    + f_6 * qsf0_493[k]
                    - f_7 * qsf1_493[k]
                    + f_3 * pc_x[k] * qsg_738[k];

        t_1033[k] = f_20 * osg_602[k]
                    + f_3 * pc_y[k] * qsg_737[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_z, osg_588, osg_740, osg_741, \
                         qsf0_495, qsf0_496, qsf1_495, qsf1_496, qsg_738, qsg_740, \
                         qsg_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_11 * osg_740[k]
                    + f_6 * qsf0_495[k]
                    - f_7 * qsf1_495[k]
                    + f_3 * pc_x[k] * qsg_740[k];

        t_1035[k] = f_11 * osg_741[k]
                    + f_4 * qsf0_496[k]
                    - f_5 * qsf1_496[k]
                    + f_3 * pc_x[k] * qsg_741[k];

        t_1036[k] = f_18 * osg_588[k]
                    + f_3 * pc_z[k] * qsg_738[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, pc_x, pc_y, osg_605, osg_744, \
                         osg_745, osg_746, qsf0_499, qsf1_499, qsg_740, qsg_744, qsg_745, \
                         qsg_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_20 * osg_605[k]
                    + f_3 * pc_y[k] * qsg_740[k];

        t_1038[k] = f_11 * osg_744[k]
                    + f_4 * qsf0_499[k]
                    - f_5 * qsf1_499[k]
                    + f_3 * pc_x[k] * qsg_744[k];

        t_1039[k] = f_11 * osg_745[k]
                    + f_3 * pc_x[k] * qsg_745[k];

        t_1040[k] = f_11 * osg_746[k]
                    + f_3 * pc_x[k] * qsg_746[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, t_1044, pc_x, pc_y, osg_610, osg_747, \
                         osg_748, osg_749, qsf0_496, qsf1_496, qsg_745, qsg_747, qsg_748, \
                         qsg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_11 * osg_747[k]
                    + f_3 * pc_x[k] * qsg_747[k];

        t_1042[k] = f_11 * osg_748[k]
                    + f_3 * pc_x[k] * qsg_748[k];

        t_1043[k] = f_11 * osg_749[k]
                    + f_3 * pc_x[k] * qsg_749[k];

        t_1044[k] = f_20 * osg_610[k]
                    + f_1 * qsf0_496[k]
                    - f_2 * qsf1_496[k]
                    + f_3 * pc_y[k] * qsg_745[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pc_y, pc_z, osg_595, osg_612, osg_613, \
                         qsf0_498, qsf0_499, qsf1_498, qsf1_499, qsg_745, qsg_747, \
                         qsg_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = f_18 * osg_595[k]
                    + f_3 * pc_z[k] * qsg_745[k];

        t_1046[k] = f_20 * osg_612[k]
                    + f_6 * qsf0_498[k]
                    - f_7 * qsf1_498[k]
                    + f_3 * pc_y[k] * qsg_747[k];

        t_1047[k] = f_20 * osg_613[k]
                    + f_4 * qsf0_499[k]
                    - f_5 * qsf1_499[k]
                    + f_3 * pc_y[k] * qsg_748[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pc_x, pc_y, pc_z, osg_599, osg_614, osg_750, \
                         qsf0_499, qsf0_500, qsf1_499, qsf1_500, qsg_749, \
                         qsg_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_20 * osg_614[k]
                    + f_3 * pc_y[k] * qsg_749[k];

        t_1049[k] = f_18 * osg_599[k]
                    + f_1 * qsf0_499[k]
                    - f_2 * qsf1_499[k]
                    + f_3 * pc_z[k] * qsg_749[k];

        t_1050[k] = f_11 * osg_750[k]
                    + f_1 * qsf0_500[k]
                    - f_2 * qsf1_500[k]
                    + f_3 * pc_x[k] * qsg_750[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, pc_x, pc_y, pc_z, osg_600, osg_615, \
                         osg_617, osg_753, qsf0_503, qsf1_503, qsg_750, qsg_752, \
                         qsg_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_18 * osg_615[k]
                    + f_3 * pc_y[k] * qsg_750[k];

        t_1052[k] = f_20 * osg_600[k]
                    + f_3 * pc_z[k] * qsg_750[k];

        t_1053[k] = f_11 * osg_753[k]
                    + f_6 * qsf0_503[k]
                    - f_7 * qsf1_503[k]
                    + f_3 * pc_x[k] * qsg_753[k];

        t_1054[k] = f_18 * osg_617[k]
                    + f_3 * pc_y[k] * qsg_752[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, osg_603, osg_755, osg_756, \
                         qsf0_505, qsf0_506, qsf1_505, qsf1_506, qsg_753, qsg_755, \
                         qsg_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_11 * osg_755[k]
                    + f_6 * qsf0_505[k]
                    - f_7 * qsf1_505[k]
                    + f_3 * pc_x[k] * qsg_755[k];

        t_1056[k] = f_11 * osg_756[k]
                    + f_4 * qsf0_506[k]
                    - f_5 * qsf1_506[k]
                    + f_3 * pc_x[k] * qsg_756[k];

        t_1057[k] = f_20 * osg_603[k]
                    + f_3 * pc_z[k] * qsg_753[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osh0,
                                                          const size_t osg, const size_t osh1,
                                                          const size_t qsf0, const size_t qsf1,
                                                          const size_t qsg, const size_t ncols,
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
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_924 = buffer.data(osh0 + 924);
    const auto *osh0_927 = buffer.data(osh0 + 927);
    const auto *osh0_929 = buffer.data(osh0 + 929);
    const auto *osh0_930 = buffer.data(osh0 + 930);
    const auto *osh0_933 = buffer.data(osh0 + 933);
    const auto *osh0_944 = buffer.data(osh0 + 944);

    const auto *osg_610 = buffer.data(osg + 610);
    const auto *osg_614 = buffer.data(osg + 614);
    const auto *osg_615 = buffer.data(osg + 615);
    const auto *osg_618 = buffer.data(osg + 618);
    const auto *osg_620 = buffer.data(osg + 620);
    const auto *osg_625 = buffer.data(osg + 625);
    const auto *osg_627 = buffer.data(osg + 627);
    const auto *osg_628 = buffer.data(osg + 628);
    const auto *osg_629 = buffer.data(osg + 629);
    const auto *osg_630 = buffer.data(osg + 630);
    const auto *osg_632 = buffer.data(osg + 632);
    const auto *osg_633 = buffer.data(osg + 633);
    const auto *osg_635 = buffer.data(osg + 635);
    const auto *osg_640 = buffer.data(osg + 640);
    const auto *osg_642 = buffer.data(osg + 642);
    const auto *osg_643 = buffer.data(osg + 643);
    const auto *osg_644 = buffer.data(osg + 644);
    const auto *osg_645 = buffer.data(osg + 645);
    const auto *osg_647 = buffer.data(osg + 647);
    const auto *osg_648 = buffer.data(osg + 648);
    const auto *osg_650 = buffer.data(osg + 650);
    const auto *osg_655 = buffer.data(osg + 655);
    const auto *osg_657 = buffer.data(osg + 657);
    const auto *osg_658 = buffer.data(osg + 658);
    const auto *osg_659 = buffer.data(osg + 659);
    const auto *osg_660 = buffer.data(osg + 660);
    const auto *osg_661 = buffer.data(osg + 661);
    const auto *osg_662 = buffer.data(osg + 662);
    const auto *osg_663 = buffer.data(osg + 663);
    const auto *osg_665 = buffer.data(osg + 665);
    const auto *osg_670 = buffer.data(osg + 670);
    const auto *osg_672 = buffer.data(osg + 672);
    const auto *osg_673 = buffer.data(osg + 673);
    const auto *osg_674 = buffer.data(osg + 674);
    const auto *osg_675 = buffer.data(osg + 675);
    const auto *osg_680 = buffer.data(osg + 680);
    const auto *osg_685 = buffer.data(osg + 685);
    const auto *osg_689 = buffer.data(osg + 689);
    const auto *osg_759 = buffer.data(osg + 759);
    const auto *osg_760 = buffer.data(osg + 760);
    const auto *osg_761 = buffer.data(osg + 761);
    const auto *osg_762 = buffer.data(osg + 762);
    const auto *osg_763 = buffer.data(osg + 763);
    const auto *osg_764 = buffer.data(osg + 764);
    const auto *osg_765 = buffer.data(osg + 765);
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
    const auto *osg_783 = buffer.data(osg + 783);
    const auto *osg_785 = buffer.data(osg + 785);
    const auto *osg_786 = buffer.data(osg + 786);
    const auto *osg_789 = buffer.data(osg + 789);
    const auto *osg_790 = buffer.data(osg + 790);
    const auto *osg_791 = buffer.data(osg + 791);
    const auto *osg_792 = buffer.data(osg + 792);
    const auto *osg_793 = buffer.data(osg + 793);
    const auto *osg_794 = buffer.data(osg + 794);
    const auto *osg_805 = buffer.data(osg + 805);
    const auto *osg_806 = buffer.data(osg + 806);
    const auto *osg_807 = buffer.data(osg + 807);
    const auto *osg_808 = buffer.data(osg + 808);
    const auto *osg_809 = buffer.data(osg + 809);
    const auto *osg_810 = buffer.data(osg + 810);
    const auto *osg_815 = buffer.data(osg + 815);
    const auto *osg_819 = buffer.data(osg + 819);
    const auto *osg_820 = buffer.data(osg + 820);
    const auto *osg_821 = buffer.data(osg + 821);
    const auto *osg_822 = buffer.data(osg + 822);
    const auto *osg_824 = buffer.data(osg + 824);
    const auto *osg_825 = buffer.data(osg + 825);
    const auto *osg_828 = buffer.data(osg + 828);
    const auto *osg_831 = buffer.data(osg + 831);
    const auto *osg_835 = buffer.data(osg + 835);
    const auto *osg_837 = buffer.data(osg + 837);
    const auto *osg_838 = buffer.data(osg + 838);
    const auto *osg_839 = buffer.data(osg + 839);

    const auto *osh1_924 = buffer.data(osh1 + 924);
    const auto *osh1_927 = buffer.data(osh1 + 927);
    const auto *osh1_929 = buffer.data(osh1 + 929);
    const auto *osh1_930 = buffer.data(osh1 + 930);
    const auto *osh1_933 = buffer.data(osh1 + 933);
    const auto *osh1_944 = buffer.data(osh1 + 944);

    const auto *qsf0_506 = buffer.data(qsf0 + 506);
    const auto *qsf0_508 = buffer.data(qsf0 + 508);
    const auto *qsf0_509 = buffer.data(qsf0 + 509);
    const auto *qsf0_510 = buffer.data(qsf0 + 510);
    const auto *qsf0_513 = buffer.data(qsf0 + 513);
    const auto *qsf0_515 = buffer.data(qsf0 + 515);
    const auto *qsf0_516 = buffer.data(qsf0 + 516);
    const auto *qsf0_518 = buffer.data(qsf0 + 518);
    const auto *qsf0_519 = buffer.data(qsf0 + 519);
    const auto *qsf0_520 = buffer.data(qsf0 + 520);
    const auto *qsf0_523 = buffer.data(qsf0 + 523);
    const auto *qsf0_525 = buffer.data(qsf0 + 525);
    const auto *qsf0_526 = buffer.data(qsf0 + 526);
    const auto *qsf0_528 = buffer.data(qsf0 + 528);
    const auto *qsf0_529 = buffer.data(qsf0 + 529);
    const auto *qsf0_536 = buffer.data(qsf0 + 536);
    const auto *qsf0_538 = buffer.data(qsf0 + 538);
    const auto *qsf0_539 = buffer.data(qsf0 + 539);
    const auto *qsf0_540 = buffer.data(qsf0 + 540);
    const auto *qsf0_541 = buffer.data(qsf0 + 541);
    const auto *qsf0_542 = buffer.data(qsf0 + 542);
    const auto *qsf0_545 = buffer.data(qsf0 + 545);
    const auto *qsf0_546 = buffer.data(qsf0 + 546);
    const auto *qsf0_547 = buffer.data(qsf0 + 547);
    const auto *qsf0_548 = buffer.data(qsf0 + 548);
    const auto *qsf0_549 = buffer.data(qsf0 + 549);
    const auto *qsf0_550 = buffer.data(qsf0 + 550);
    const auto *qsf0_552 = buffer.data(qsf0 + 552);
    const auto *qsf0_553 = buffer.data(qsf0 + 553);
    const auto *qsf0_556 = buffer.data(qsf0 + 556);
    const auto *qsf0_557 = buffer.data(qsf0 + 557);

    const auto *qsf1_506 = buffer.data(qsf1 + 506);
    const auto *qsf1_508 = buffer.data(qsf1 + 508);
    const auto *qsf1_509 = buffer.data(qsf1 + 509);
    const auto *qsf1_510 = buffer.data(qsf1 + 510);
    const auto *qsf1_513 = buffer.data(qsf1 + 513);
    const auto *qsf1_515 = buffer.data(qsf1 + 515);
    const auto *qsf1_516 = buffer.data(qsf1 + 516);
    const auto *qsf1_518 = buffer.data(qsf1 + 518);
    const auto *qsf1_519 = buffer.data(qsf1 + 519);
    const auto *qsf1_520 = buffer.data(qsf1 + 520);
    const auto *qsf1_523 = buffer.data(qsf1 + 523);
    const auto *qsf1_525 = buffer.data(qsf1 + 525);
    const auto *qsf1_526 = buffer.data(qsf1 + 526);
    const auto *qsf1_528 = buffer.data(qsf1 + 528);
    const auto *qsf1_529 = buffer.data(qsf1 + 529);
    const auto *qsf1_536 = buffer.data(qsf1 + 536);
    const auto *qsf1_538 = buffer.data(qsf1 + 538);
    const auto *qsf1_539 = buffer.data(qsf1 + 539);
    const auto *qsf1_540 = buffer.data(qsf1 + 540);
    const auto *qsf1_541 = buffer.data(qsf1 + 541);
    const auto *qsf1_542 = buffer.data(qsf1 + 542);
    const auto *qsf1_545 = buffer.data(qsf1 + 545);
    const auto *qsf1_546 = buffer.data(qsf1 + 546);
    const auto *qsf1_547 = buffer.data(qsf1 + 547);
    const auto *qsf1_548 = buffer.data(qsf1 + 548);
    const auto *qsf1_549 = buffer.data(qsf1 + 549);
    const auto *qsf1_550 = buffer.data(qsf1 + 550);
    const auto *qsf1_552 = buffer.data(qsf1 + 552);
    const auto *qsf1_553 = buffer.data(qsf1 + 553);
    const auto *qsf1_556 = buffer.data(qsf1 + 556);
    const auto *qsf1_557 = buffer.data(qsf1 + 557);

    const auto *qsg_755 = buffer.data(qsg + 755);
    const auto *qsg_759 = buffer.data(qsg + 759);
    const auto *qsg_760 = buffer.data(qsg + 760);
    const auto *qsg_761 = buffer.data(qsg + 761);
    const auto *qsg_762 = buffer.data(qsg + 762);
    const auto *qsg_763 = buffer.data(qsg + 763);
    const auto *qsg_764 = buffer.data(qsg + 764);
    const auto *qsg_765 = buffer.data(qsg + 765);
    const auto *qsg_767 = buffer.data(qsg + 767);
    const auto *qsg_768 = buffer.data(qsg + 768);
    const auto *qsg_770 = buffer.data(qsg + 770);
    const auto *qsg_771 = buffer.data(qsg + 771);
    const auto *qsg_774 = buffer.data(qsg + 774);
    const auto *qsg_775 = buffer.data(qsg + 775);
    const auto *qsg_776 = buffer.data(qsg + 776);
    const auto *qsg_777 = buffer.data(qsg + 777);
    const auto *qsg_778 = buffer.data(qsg + 778);
    const auto *qsg_779 = buffer.data(qsg + 779);
    const auto *qsg_780 = buffer.data(qsg + 780);
    const auto *qsg_782 = buffer.data(qsg + 782);
    const auto *qsg_783 = buffer.data(qsg + 783);
    const auto *qsg_785 = buffer.data(qsg + 785);
    const auto *qsg_786 = buffer.data(qsg + 786);
    const auto *qsg_789 = buffer.data(qsg + 789);
    const auto *qsg_790 = buffer.data(qsg + 790);
    const auto *qsg_791 = buffer.data(qsg + 791);
    const auto *qsg_792 = buffer.data(qsg + 792);
    const auto *qsg_793 = buffer.data(qsg + 793);
    const auto *qsg_794 = buffer.data(qsg + 794);
    const auto *qsg_795 = buffer.data(qsg + 795);
    const auto *qsg_797 = buffer.data(qsg + 797);
    const auto *qsg_798 = buffer.data(qsg + 798);
    const auto *qsg_800 = buffer.data(qsg + 800);
    const auto *qsg_805 = buffer.data(qsg + 805);
    const auto *qsg_806 = buffer.data(qsg + 806);
    const auto *qsg_807 = buffer.data(qsg + 807);
    const auto *qsg_808 = buffer.data(qsg + 808);
    const auto *qsg_809 = buffer.data(qsg + 809);
    const auto *qsg_810 = buffer.data(qsg + 810);
    const auto *qsg_811 = buffer.data(qsg + 811);
    const auto *qsg_812 = buffer.data(qsg + 812);
    const auto *qsg_813 = buffer.data(qsg + 813);
    const auto *qsg_814 = buffer.data(qsg + 814);
    const auto *qsg_815 = buffer.data(qsg + 815);
    const auto *qsg_819 = buffer.data(qsg + 819);
    const auto *qsg_820 = buffer.data(qsg + 820);
    const auto *qsg_821 = buffer.data(qsg + 821);
    const auto *qsg_822 = buffer.data(qsg + 822);
    const auto *qsg_823 = buffer.data(qsg + 823);
    const auto *qsg_824 = buffer.data(qsg + 824);
    const auto *qsg_825 = buffer.data(qsg + 825);
    const auto *qsg_826 = buffer.data(qsg + 826);
    const auto *qsg_827 = buffer.data(qsg + 827);
    const auto *qsg_828 = buffer.data(qsg + 828);
    const auto *qsg_830 = buffer.data(qsg + 830);
    const auto *qsg_831 = buffer.data(qsg + 831);
    const auto *qsg_835 = buffer.data(qsg + 835);
    const auto *qsg_836 = buffer.data(qsg + 836);
    const auto *qsg_837 = buffer.data(qsg + 837);
    const auto *qsg_838 = buffer.data(qsg + 838);
    const auto *qsg_839 = buffer.data(qsg + 839);

#pragma omp simd aligned(t_1058, t_1059, t_1060, t_1061, pc_x, pc_y, osg_620, osg_759, \
                         osg_760, osg_761, qsf0_509, qsf1_509, qsg_755, qsg_759, qsg_760, \
                         qsg_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_18 * osg_620[k]
                    + f_3 * pc_y[k] * qsg_755[k];

        t_1059[k] = f_11 * osg_759[k]
                    + f_4 * qsf0_509[k]
                    - f_5 * qsf1_509[k]
                    + f_3 * pc_x[k] * qsg_759[k];

        t_1060[k] = f_11 * osg_760[k]
                    + f_3 * pc_x[k] * qsg_760[k];

        t_1061[k] = f_11 * osg_761[k]
                    + f_3 * pc_x[k] * qsg_761[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, pc_x, pc_y, osg_625, osg_762, \
                         osg_763, osg_764, qsf0_506, qsf1_506, qsg_760, qsg_762, qsg_763, \
                         qsg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_11 * osg_762[k]
                    + f_3 * pc_x[k] * qsg_762[k];

        t_1063[k] = f_11 * osg_763[k]
                    + f_3 * pc_x[k] * qsg_763[k];

        t_1064[k] = f_11 * osg_764[k]
                    + f_3 * pc_x[k] * qsg_764[k];

        t_1065[k] = f_18 * osg_625[k]
                    + f_1 * qsf0_506[k]
                    - f_2 * qsf1_506[k]
                    + f_3 * pc_y[k] * qsg_760[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, pc_y, pc_z, osg_610, osg_627, osg_628, \
                         qsf0_508, qsf0_509, qsf1_508, qsf1_509, qsg_760, qsg_762, \
                         qsg_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_20 * osg_610[k]
                    + f_3 * pc_z[k] * qsg_760[k];

        t_1067[k] = f_18 * osg_627[k]
                    + f_6 * qsf0_508[k]
                    - f_7 * qsf1_508[k]
                    + f_3 * pc_y[k] * qsg_762[k];

        t_1068[k] = f_18 * osg_628[k]
                    + f_4 * qsf0_509[k]
                    - f_5 * qsf1_509[k]
                    + f_3 * pc_y[k] * qsg_763[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, pc_x, pc_y, pc_z, osg_614, osg_629, osg_765, \
                         qsf0_509, qsf0_510, qsf1_509, qsf1_510, qsg_764, \
                         qsg_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_18 * osg_629[k]
                    + f_3 * pc_y[k] * qsg_764[k];

        t_1070[k] = f_20 * osg_614[k]
                    + f_1 * qsf0_509[k]
                    - f_2 * qsf1_509[k]
                    + f_3 * pc_z[k] * qsg_764[k];

        t_1071[k] = f_11 * osg_765[k]
                    + f_1 * qsf0_510[k]
                    - f_2 * qsf1_510[k]
                    + f_3 * pc_x[k] * qsg_765[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, pc_x, pc_y, pc_z, osg_615, osg_630, \
                         osg_632, osg_768, qsf0_513, qsf1_513, qsg_765, qsg_767, \
                         qsg_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = f_11 * osg_630[k]
                    + f_3 * pc_y[k] * qsg_765[k];

        t_1073[k] = f_21 * osg_615[k]
                    + f_3 * pc_z[k] * qsg_765[k];

        t_1074[k] = f_11 * osg_768[k]
                    + f_6 * qsf0_513[k]
                    - f_7 * qsf1_513[k]
                    + f_3 * pc_x[k] * qsg_768[k];

        t_1075[k] = f_11 * osg_632[k]
                    + f_3 * pc_y[k] * qsg_767[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_z, osg_618, osg_770, osg_771, \
                         qsf0_515, qsf0_516, qsf1_515, qsf1_516, qsg_768, qsg_770, \
                         qsg_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_11 * osg_770[k]
                    + f_6 * qsf0_515[k]
                    - f_7 * qsf1_515[k]
                    + f_3 * pc_x[k] * qsg_770[k];

        t_1077[k] = f_11 * osg_771[k]
                    + f_4 * qsf0_516[k]
                    - f_5 * qsf1_516[k]
                    + f_3 * pc_x[k] * qsg_771[k];

        t_1078[k] = f_21 * osg_618[k]
                    + f_3 * pc_z[k] * qsg_768[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pc_x, pc_y, osg_635, osg_774, \
                         osg_775, osg_776, qsf0_519, qsf1_519, qsg_770, qsg_774, qsg_775, \
                         qsg_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_11 * osg_635[k]
                    + f_3 * pc_y[k] * qsg_770[k];

        t_1080[k] = f_11 * osg_774[k]
                    + f_4 * qsf0_519[k]
                    - f_5 * qsf1_519[k]
                    + f_3 * pc_x[k] * qsg_774[k];

        t_1081[k] = f_11 * osg_775[k]
                    + f_3 * pc_x[k] * qsg_775[k];

        t_1082[k] = f_11 * osg_776[k]
                    + f_3 * pc_x[k] * qsg_776[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, pc_x, pc_y, osg_640, osg_777, \
                         osg_778, osg_779, qsf0_516, qsf1_516, qsg_775, qsg_777, qsg_778, \
                         qsg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_11 * osg_777[k]
                    + f_3 * pc_x[k] * qsg_777[k];

        t_1084[k] = f_11 * osg_778[k]
                    + f_3 * pc_x[k] * qsg_778[k];

        t_1085[k] = f_11 * osg_779[k]
                    + f_3 * pc_x[k] * qsg_779[k];

        t_1086[k] = f_11 * osg_640[k]
                    + f_1 * qsf0_516[k]
                    - f_2 * qsf1_516[k]
                    + f_3 * pc_y[k] * qsg_775[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, pc_z, osg_625, osg_642, osg_643, \
                         qsf0_518, qsf0_519, qsf1_518, qsf1_519, qsg_775, qsg_777, \
                         qsg_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_21 * osg_625[k]
                    + f_3 * pc_z[k] * qsg_775[k];

        t_1088[k] = f_11 * osg_642[k]
                    + f_6 * qsf0_518[k]
                    - f_7 * qsf1_518[k]
                    + f_3 * pc_y[k] * qsg_777[k];

        t_1089[k] = f_11 * osg_643[k]
                    + f_4 * qsf0_519[k]
                    - f_5 * qsf1_519[k]
                    + f_3 * pc_y[k] * qsg_778[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, osg_629, osg_644, osg_780, \
                         qsf0_519, qsf0_520, qsf1_519, qsf1_520, qsg_779, \
                         qsg_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_11 * osg_644[k]
                    + f_3 * pc_y[k] * qsg_779[k];

        t_1091[k] = f_21 * osg_629[k]
                    + f_1 * qsf0_519[k]
                    - f_2 * qsf1_519[k]
                    + f_3 * pc_z[k] * qsg_779[k];

        t_1092[k] = f_11 * osg_780[k]
                    + f_1 * qsf0_520[k]
                    - f_2 * qsf1_520[k]
                    + f_3 * pc_x[k] * qsg_780[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, osg_630, osg_645, \
                         osg_647, osg_783, qsf0_523, qsf1_523, qsg_780, qsg_782, \
                         qsg_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_10 * osg_645[k]
                    + f_3 * pc_y[k] * qsg_780[k];

        t_1094[k] = f_19 * osg_630[k]
                    + f_3 * pc_z[k] * qsg_780[k];

        t_1095[k] = f_11 * osg_783[k]
                    + f_6 * qsf0_523[k]
                    - f_7 * qsf1_523[k]
                    + f_3 * pc_x[k] * qsg_783[k];

        t_1096[k] = f_10 * osg_647[k]
                    + f_3 * pc_y[k] * qsg_782[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, osg_633, osg_785, osg_786, \
                         qsf0_525, qsf0_526, qsf1_525, qsf1_526, qsg_783, qsg_785, \
                         qsg_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_11 * osg_785[k]
                    + f_6 * qsf0_525[k]
                    - f_7 * qsf1_525[k]
                    + f_3 * pc_x[k] * qsg_785[k];

        t_1098[k] = f_11 * osg_786[k]
                    + f_4 * qsf0_526[k]
                    - f_5 * qsf1_526[k]
                    + f_3 * pc_x[k] * qsg_786[k];

        t_1099[k] = f_19 * osg_633[k]
                    + f_3 * pc_z[k] * qsg_783[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, pc_y, osg_650, osg_789, \
                         osg_790, osg_791, qsf0_529, qsf1_529, qsg_785, qsg_789, qsg_790, \
                         qsg_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_10 * osg_650[k]
                    + f_3 * pc_y[k] * qsg_785[k];

        t_1101[k] = f_11 * osg_789[k]
                    + f_4 * qsf0_529[k]
                    - f_5 * qsf1_529[k]
                    + f_3 * pc_x[k] * qsg_789[k];

        t_1102[k] = f_11 * osg_790[k]
                    + f_3 * pc_x[k] * qsg_790[k];

        t_1103[k] = f_11 * osg_791[k]
                    + f_3 * pc_x[k] * qsg_791[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, pc_y, osg_655, osg_792, \
                         osg_793, osg_794, qsf0_526, qsf1_526, qsg_790, qsg_792, qsg_793, \
                         qsg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_11 * osg_792[k]
                    + f_3 * pc_x[k] * qsg_792[k];

        t_1105[k] = f_11 * osg_793[k]
                    + f_3 * pc_x[k] * qsg_793[k];

        t_1106[k] = f_11 * osg_794[k]
                    + f_3 * pc_x[k] * qsg_794[k];

        t_1107[k] = f_10 * osg_655[k]
                    + f_1 * qsf0_526[k]
                    - f_2 * qsf1_526[k]
                    + f_3 * pc_y[k] * qsg_790[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, osg_640, osg_657, osg_658, \
                         qsf0_528, qsf0_529, qsf1_528, qsf1_529, qsg_790, qsg_792, \
                         qsg_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_19 * osg_640[k]
                    + f_3 * pc_z[k] * qsg_790[k];

        t_1109[k] = f_10 * osg_657[k]
                    + f_6 * qsf0_528[k]
                    - f_7 * qsf1_528[k]
                    + f_3 * pc_y[k] * qsg_792[k];

        t_1110[k] = f_10 * osg_658[k]
                    + f_4 * qsf0_529[k]
                    - f_5 * qsf1_529[k]
                    + f_3 * pc_y[k] * qsg_793[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, pa_y, pc_y, pc_z, osh0_924, osg_644, \
                         osg_659, osg_660, osh1_924, qsf0_529, qsf1_529, qsg_794, \
                         qsg_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_10 * osg_659[k]
                    + f_3 * pc_y[k] * qsg_794[k];

        t_1112[k] = f_19 * osg_644[k]
                    + f_1 * qsf0_529[k]
                    - f_2 * qsf1_529[k]
                    + f_3 * pc_z[k] * qsg_794[k];

        t_1113[k] = pa_y[k] * osh0_924[k]
                    - f_8 * pc_y[k] * osh1_924[k];

        t_1114[k] = f_9 * osg_660[k]
                    + f_3 * pc_y[k] * qsg_795[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, pa_y, pc_y, pc_z, osh0_927, osh0_929, \
                         osg_645, osg_661, osg_662, osh1_927, osh1_929, qsg_795, \
                         qsg_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_17 * osg_645[k]
                    + f_3 * pc_z[k] * qsg_795[k];

        t_1116[k] = pa_y[k] * osh0_927[k]
                    + f_10 * osg_661[k]
                    - f_8 * pc_y[k] * osh1_927[k];

        t_1117[k] = f_9 * osg_662[k]
                    + f_3 * pc_y[k] * qsg_797[k];

        t_1118[k] = pa_y[k] * osh0_929[k]
                    - f_8 * pc_y[k] * osh1_929[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, t_1122, pa_y, pc_y, pc_z, osh0_930, osh0_933, \
                         osg_648, osg_663, osg_665, osh1_930, osh1_933, qsg_798, \
                         qsg_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pa_y[k] * osh0_930[k]
                    + f_11 * osg_663[k]
                    - f_8 * pc_y[k] * osh1_930[k];

        t_1120[k] = f_17 * osg_648[k]
                    + f_3 * pc_z[k] * qsg_798[k];

        t_1121[k] = f_9 * osg_665[k]
                    + f_3 * pc_y[k] * qsg_800[k];

        t_1122[k] = pa_y[k] * osh0_933[k]
                    - f_8 * pc_y[k] * osh1_933[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, t_1127, pc_x, osg_805, osg_806, \
                         osg_807, osg_808, osg_809, qsg_805, qsg_806, qsg_807, qsg_808, \
                         qsg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_11 * osg_805[k]
                    + f_3 * pc_x[k] * qsg_805[k];

        t_1124[k] = f_11 * osg_806[k]
                    + f_3 * pc_x[k] * qsg_806[k];

        t_1125[k] = f_11 * osg_807[k]
                    + f_3 * pc_x[k] * qsg_807[k];

        t_1126[k] = f_11 * osg_808[k]
                    + f_3 * pc_x[k] * qsg_808[k];

        t_1127[k] = f_11 * osg_809[k]
                    + f_3 * pc_x[k] * qsg_809[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pc_y, pc_z, osg_655, osg_670, osg_672, \
                         qsf0_536, qsf0_538, qsf1_536, qsf1_538, qsg_805, \
                         qsg_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_9 * osg_670[k]
                    + f_1 * qsf0_536[k]
                    - f_2 * qsf1_536[k]
                    + f_3 * pc_y[k] * qsg_805[k];

        t_1129[k] = f_17 * osg_655[k]
                    + f_3 * pc_z[k] * qsg_805[k];

        t_1130[k] = f_9 * osg_672[k]
                    + f_6 * qsf0_538[k]
                    - f_7 * qsf1_538[k]
                    + f_3 * pc_y[k] * qsg_807[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pa_y, pc_y, osh0_944, osg_673, osg_674, \
                         osh1_944, qsf0_539, qsf1_539, qsg_808, \
                         qsg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_9 * osg_673[k]
                    + f_4 * qsf0_539[k]
                    - f_5 * qsf1_539[k]
                    + f_3 * pc_y[k] * qsg_808[k];

        t_1132[k] = f_9 * osg_674[k]
                    + f_3 * pc_y[k] * qsg_809[k];

        t_1133[k] = pa_y[k] * osh0_944[k]
                    - f_8 * pc_y[k] * osh1_944[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, osg_660, \
                         osg_810, qsf0_540, qsf1_540, qsg_810, qsg_811, \
                         qsg_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_11 * osg_810[k]
                    + f_1 * qsf0_540[k]
                    - f_2 * qsf1_540[k]
                    + f_3 * pc_x[k] * qsg_810[k];

        t_1135[k] = f_3 * pc_y[k] * qsg_810[k];

        t_1136[k] = f_16 * osg_660[k]
                    + f_3 * pc_z[k] * qsg_810[k];

        t_1137[k] = f_4 * qsf0_540[k]
                    - f_5 * qsf1_540[k]
                    + f_3 * pc_y[k] * qsg_811[k];

        t_1138[k] = f_3 * pc_y[k] * qsg_812[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, pc_x, pc_y, osg_815, qsf0_541, \
                         qsf0_542, qsf0_545, qsf1_541, qsf1_542, qsf1_545, qsg_813, qsg_814, \
                         qsg_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_11 * osg_815[k]
                    + f_6 * qsf0_545[k]
                    - f_7 * qsf1_545[k]
                    + f_3 * pc_x[k] * qsg_815[k];

        t_1140[k] = f_6 * qsf0_541[k]
                    - f_7 * qsf1_541[k]
                    + f_3 * pc_y[k] * qsg_813[k];

        t_1141[k] = f_4 * qsf0_542[k]
                    - f_5 * qsf1_542[k]
                    + f_3 * pc_y[k] * qsg_814[k];

        t_1142[k] = f_3 * pc_y[k] * qsg_815[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, t_1146, pc_x, osg_819, osg_820, osg_821, \
                         osg_822, qsf0_549, qsf1_549, qsg_819, qsg_820, qsg_821, \
                         qsg_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_11 * osg_819[k]
                    + f_4 * qsf0_549[k]
                    - f_5 * qsf1_549[k]
                    + f_3 * pc_x[k] * qsg_819[k];

        t_1144[k] = f_11 * osg_820[k]
                    + f_3 * pc_x[k] * qsg_820[k];

        t_1145[k] = f_11 * osg_821[k]
                    + f_3 * pc_x[k] * qsg_821[k];

        t_1146[k] = f_11 * osg_822[k]
                    + f_3 * pc_x[k] * qsg_822[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, pc_x, pc_y, osg_824, qsf0_546, \
                         qsf0_547, qsf1_546, qsf1_547, qsg_819, qsg_820, qsg_821, \
                         qsg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_3 * pc_y[k] * qsg_819[k];

        t_1148[k] = f_11 * osg_824[k]
                    + f_3 * pc_x[k] * qsg_824[k];

        t_1149[k] = f_1 * qsf0_546[k]
                    - f_2 * qsf1_546[k]
                    + f_3 * pc_y[k] * qsg_820[k];

        t_1150[k] = f_13 * qsf0_547[k]
                    - f_14 * qsf1_547[k]
                    + f_3 * pc_y[k] * qsg_821[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pc_y, pc_z, osg_674, qsf0_548, \
                         qsf0_549, qsf1_548, qsf1_549, qsg_822, qsg_823, \
                         qsg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_6 * qsf0_548[k]
                    - f_7 * qsf1_548[k]
                    + f_3 * pc_y[k] * qsg_822[k];

        t_1152[k] = f_4 * qsf0_549[k]
                    - f_5 * qsf1_549[k]
                    + f_3 * pc_y[k] * qsg_823[k];

        t_1153[k] = f_3 * pc_y[k] * qsg_824[k];

        t_1154[k] = f_16 * osg_674[k]
                    + f_1 * qsf0_549[k]
                    - f_2 * qsf1_549[k]
                    + f_3 * pc_z[k] * qsg_824[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, pc_x, pc_y, pc_z, osg_675, osg_825, \
                         osg_828, qsf0_550, qsf0_553, qsf1_550, qsf1_553, qsg_825, \
                         qsg_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_10 * osg_825[k]
                    + f_1 * qsf0_550[k]
                    - f_2 * qsf1_550[k]
                    + f_3 * pc_x[k] * qsg_825[k];

        t_1156[k] = f_15 * osg_675[k]
                    + f_3 * pc_y[k] * qsg_825[k];

        t_1157[k] = f_3 * pc_z[k] * qsg_825[k];

        t_1158[k] = f_10 * osg_828[k]
                    + f_6 * qsf0_553[k]
                    - f_7 * qsf1_553[k]
                    + f_3 * pc_x[k] * qsg_828[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, t_1162, pc_x, pc_z, osg_831, qsf0_550, \
                         qsf0_556, qsf1_550, qsf1_556, qsg_826, qsg_827, qsg_828, \
                         qsg_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_3 * pc_z[k] * qsg_826[k];

        t_1160[k] = f_4 * qsf0_550[k]
                    - f_5 * qsf1_550[k]
                    + f_3 * pc_z[k] * qsg_827[k];

        t_1161[k] = f_10 * osg_831[k]
                    + f_4 * qsf0_556[k]
                    - f_5 * qsf1_556[k]
                    + f_3 * pc_x[k] * qsg_831[k];

        t_1162[k] = f_3 * pc_z[k] * qsg_828[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, t_1166, pc_x, pc_y, pc_z, osg_680, osg_835, \
                         qsf0_552, qsf1_552, qsg_830, qsg_831, \
                         qsg_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_15 * osg_680[k]
                    + f_3 * pc_y[k] * qsg_830[k];

        t_1164[k] = f_6 * qsf0_552[k]
                    - f_7 * qsf1_552[k]
                    + f_3 * pc_z[k] * qsg_830[k];

        t_1165[k] = f_10 * osg_835[k]
                    + f_3 * pc_x[k] * qsg_835[k];

        t_1166[k] = f_3 * pc_z[k] * qsg_831[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, pc_x, pc_y, osg_685, osg_837, \
                         osg_838, osg_839, qsf0_556, qsf1_556, qsg_835, qsg_837, qsg_838, \
                         qsg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_10 * osg_837[k]
                    + f_3 * pc_x[k] * qsg_837[k];

        t_1168[k] = f_10 * osg_838[k]
                    + f_3 * pc_x[k] * qsg_838[k];

        t_1169[k] = f_10 * osg_839[k]
                    + f_3 * pc_x[k] * qsg_839[k];

        t_1170[k] = f_15 * osg_685[k]
                    + f_1 * qsf0_556[k]
                    - f_2 * qsf1_556[k]
                    + f_3 * pc_y[k] * qsg_835[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, pc_y, pc_z, osg_689, qsf0_556, \
                         qsf0_557, qsf1_556, qsf1_557, qsg_835, qsg_836, qsg_837, \
                         qsg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_3 * pc_z[k] * qsg_835[k];

        t_1172[k] = f_4 * qsf0_556[k]
                    - f_5 * qsf1_556[k]
                    + f_3 * pc_z[k] * qsg_836[k];

        t_1173[k] = f_6 * qsf0_557[k]
                    - f_7 * qsf1_557[k]
                    + f_3 * pc_z[k] * qsg_837[k];

        t_1174[k] = f_15 * osg_689[k]
                    + f_3 * pc_y[k] * qsg_839[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsf0, const size_t qsf1,
                                                           const size_t qsg, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_945 = buffer.data(osh0 + 945);
    const auto *osh0_948 = buffer.data(osh0 + 948);
    const auto *osh0_951 = buffer.data(osh0 + 951);
    const auto *osh0_960 = buffer.data(osh0 + 960);

    const auto *osg_675 = buffer.data(osg + 675);
    const auto *osg_678 = buffer.data(osg + 678);
    const auto *osg_685 = buffer.data(osg + 685);
    const auto *osg_689 = buffer.data(osg + 689);
    const auto *osg_690 = buffer.data(osg + 690);
    const auto *osg_692 = buffer.data(osg + 692);
    const auto *osg_693 = buffer.data(osg + 693);
    const auto *osg_695 = buffer.data(osg + 695);
    const auto *osg_700 = buffer.data(osg + 700);
    const auto *osg_702 = buffer.data(osg + 702);
    const auto *osg_703 = buffer.data(osg + 703);
    const auto *osg_704 = buffer.data(osg + 704);
    const auto *osg_705 = buffer.data(osg + 705);
    const auto *osg_707 = buffer.data(osg + 707);
    const auto *osg_708 = buffer.data(osg + 708);
    const auto *osg_710 = buffer.data(osg + 710);
    const auto *osg_715 = buffer.data(osg + 715);
    const auto *osg_717 = buffer.data(osg + 717);
    const auto *osg_718 = buffer.data(osg + 718);
    const auto *osg_719 = buffer.data(osg + 719);
    const auto *osg_720 = buffer.data(osg + 720);
    const auto *osg_722 = buffer.data(osg + 722);
    const auto *osg_723 = buffer.data(osg + 723);
    const auto *osg_725 = buffer.data(osg + 725);
    const auto *osg_730 = buffer.data(osg + 730);
    const auto *osg_732 = buffer.data(osg + 732);
    const auto *osg_733 = buffer.data(osg + 733);
    const auto *osg_734 = buffer.data(osg + 734);
    const auto *osg_735 = buffer.data(osg + 735);
    const auto *osg_737 = buffer.data(osg + 737);
    const auto *osg_738 = buffer.data(osg + 738);
    const auto *osg_740 = buffer.data(osg + 740);
    const auto *osg_745 = buffer.data(osg + 745);
    const auto *osg_747 = buffer.data(osg + 747);
    const auto *osg_748 = buffer.data(osg + 748);
    const auto *osg_749 = buffer.data(osg + 749);
    const auto *osg_750 = buffer.data(osg + 750);
    const auto *osg_752 = buffer.data(osg + 752);
    const auto *osg_755 = buffer.data(osg + 755);
    const auto *osg_760 = buffer.data(osg + 760);
    const auto *osg_762 = buffer.data(osg + 762);
    const auto *osg_763 = buffer.data(osg + 763);
    const auto *osg_764 = buffer.data(osg + 764);
    const auto *osg_845 = buffer.data(osg + 845);
    const auto *osg_849 = buffer.data(osg + 849);
    const auto *osg_850 = buffer.data(osg + 850);
    const auto *osg_851 = buffer.data(osg + 851);
    const auto *osg_852 = buffer.data(osg + 852);
    const auto *osg_853 = buffer.data(osg + 853);
    const auto *osg_854 = buffer.data(osg + 854);
    const auto *osg_855 = buffer.data(osg + 855);
    const auto *osg_858 = buffer.data(osg + 858);
    const auto *osg_860 = buffer.data(osg + 860);
    const auto *osg_861 = buffer.data(osg + 861);
    const auto *osg_864 = buffer.data(osg + 864);
    const auto *osg_865 = buffer.data(osg + 865);
    const auto *osg_866 = buffer.data(osg + 866);
    const auto *osg_867 = buffer.data(osg + 867);
    const auto *osg_868 = buffer.data(osg + 868);
    const auto *osg_869 = buffer.data(osg + 869);
    const auto *osg_870 = buffer.data(osg + 870);
    const auto *osg_873 = buffer.data(osg + 873);
    const auto *osg_875 = buffer.data(osg + 875);
    const auto *osg_876 = buffer.data(osg + 876);
    const auto *osg_879 = buffer.data(osg + 879);
    const auto *osg_880 = buffer.data(osg + 880);
    const auto *osg_881 = buffer.data(osg + 881);
    const auto *osg_882 = buffer.data(osg + 882);
    const auto *osg_883 = buffer.data(osg + 883);
    const auto *osg_884 = buffer.data(osg + 884);
    const auto *osg_885 = buffer.data(osg + 885);
    const auto *osg_888 = buffer.data(osg + 888);
    const auto *osg_890 = buffer.data(osg + 890);
    const auto *osg_891 = buffer.data(osg + 891);
    const auto *osg_894 = buffer.data(osg + 894);
    const auto *osg_895 = buffer.data(osg + 895);
    const auto *osg_896 = buffer.data(osg + 896);
    const auto *osg_897 = buffer.data(osg + 897);
    const auto *osg_898 = buffer.data(osg + 898);
    const auto *osg_899 = buffer.data(osg + 899);
    const auto *osg_900 = buffer.data(osg + 900);
    const auto *osg_903 = buffer.data(osg + 903);
    const auto *osg_905 = buffer.data(osg + 905);
    const auto *osg_906 = buffer.data(osg + 906);
    const auto *osg_909 = buffer.data(osg + 909);
    const auto *osg_910 = buffer.data(osg + 910);
    const auto *osg_911 = buffer.data(osg + 911);
    const auto *osg_912 = buffer.data(osg + 912);
    const auto *osg_913 = buffer.data(osg + 913);
    const auto *osg_914 = buffer.data(osg + 914);
    const auto *osg_915 = buffer.data(osg + 915);

    const auto *osh1_945 = buffer.data(osh1 + 945);
    const auto *osh1_948 = buffer.data(osh1 + 948);
    const auto *osh1_951 = buffer.data(osh1 + 951);
    const auto *osh1_960 = buffer.data(osh1 + 960);

    const auto *qsf0_559 = buffer.data(qsf0 + 559);
    const auto *qsf0_565 = buffer.data(qsf0 + 565);
    const auto *qsf0_568 = buffer.data(qsf0 + 568);
    const auto *qsf0_569 = buffer.data(qsf0 + 569);
    const auto *qsf0_570 = buffer.data(qsf0 + 570);
    const auto *qsf0_573 = buffer.data(qsf0 + 573);
    const auto *qsf0_575 = buffer.data(qsf0 + 575);
    const auto *qsf0_576 = buffer.data(qsf0 + 576);
    const auto *qsf0_578 = buffer.data(qsf0 + 578);
    const auto *qsf0_579 = buffer.data(qsf0 + 579);
    const auto *qsf0_580 = buffer.data(qsf0 + 580);
    const auto *qsf0_583 = buffer.data(qsf0 + 583);
    const auto *qsf0_585 = buffer.data(qsf0 + 585);
    const auto *qsf0_586 = buffer.data(qsf0 + 586);
    const auto *qsf0_588 = buffer.data(qsf0 + 588);
    const auto *qsf0_589 = buffer.data(qsf0 + 589);
    const auto *qsf0_590 = buffer.data(qsf0 + 590);
    const auto *qsf0_593 = buffer.data(qsf0 + 593);
    const auto *qsf0_595 = buffer.data(qsf0 + 595);
    const auto *qsf0_596 = buffer.data(qsf0 + 596);
    const auto *qsf0_598 = buffer.data(qsf0 + 598);
    const auto *qsf0_599 = buffer.data(qsf0 + 599);
    const auto *qsf0_600 = buffer.data(qsf0 + 600);
    const auto *qsf0_603 = buffer.data(qsf0 + 603);
    const auto *qsf0_605 = buffer.data(qsf0 + 605);
    const auto *qsf0_606 = buffer.data(qsf0 + 606);
    const auto *qsf0_608 = buffer.data(qsf0 + 608);
    const auto *qsf0_609 = buffer.data(qsf0 + 609);
    const auto *qsf0_610 = buffer.data(qsf0 + 610);

    const auto *qsf1_559 = buffer.data(qsf1 + 559);
    const auto *qsf1_565 = buffer.data(qsf1 + 565);
    const auto *qsf1_568 = buffer.data(qsf1 + 568);
    const auto *qsf1_569 = buffer.data(qsf1 + 569);
    const auto *qsf1_570 = buffer.data(qsf1 + 570);
    const auto *qsf1_573 = buffer.data(qsf1 + 573);
    const auto *qsf1_575 = buffer.data(qsf1 + 575);
    const auto *qsf1_576 = buffer.data(qsf1 + 576);
    const auto *qsf1_578 = buffer.data(qsf1 + 578);
    const auto *qsf1_579 = buffer.data(qsf1 + 579);
    const auto *qsf1_580 = buffer.data(qsf1 + 580);
    const auto *qsf1_583 = buffer.data(qsf1 + 583);
    const auto *qsf1_585 = buffer.data(qsf1 + 585);
    const auto *qsf1_586 = buffer.data(qsf1 + 586);
    const auto *qsf1_588 = buffer.data(qsf1 + 588);
    const auto *qsf1_589 = buffer.data(qsf1 + 589);
    const auto *qsf1_590 = buffer.data(qsf1 + 590);
    const auto *qsf1_593 = buffer.data(qsf1 + 593);
    const auto *qsf1_595 = buffer.data(qsf1 + 595);
    const auto *qsf1_596 = buffer.data(qsf1 + 596);
    const auto *qsf1_598 = buffer.data(qsf1 + 598);
    const auto *qsf1_599 = buffer.data(qsf1 + 599);
    const auto *qsf1_600 = buffer.data(qsf1 + 600);
    const auto *qsf1_603 = buffer.data(qsf1 + 603);
    const auto *qsf1_605 = buffer.data(qsf1 + 605);
    const auto *qsf1_606 = buffer.data(qsf1 + 606);
    const auto *qsf1_608 = buffer.data(qsf1 + 608);
    const auto *qsf1_609 = buffer.data(qsf1 + 609);
    const auto *qsf1_610 = buffer.data(qsf1 + 610);

    const auto *qsg_839 = buffer.data(qsg + 839);
    const auto *qsg_840 = buffer.data(qsg + 840);
    const auto *qsg_842 = buffer.data(qsg + 842);
    const auto *qsg_843 = buffer.data(qsg + 843);
    const auto *qsg_845 = buffer.data(qsg + 845);
    const auto *qsg_849 = buffer.data(qsg + 849);
    const auto *qsg_850 = buffer.data(qsg + 850);
    const auto *qsg_851 = buffer.data(qsg + 851);
    const auto *qsg_852 = buffer.data(qsg + 852);
    const auto *qsg_853 = buffer.data(qsg + 853);
    const auto *qsg_854 = buffer.data(qsg + 854);
    const auto *qsg_855 = buffer.data(qsg + 855);
    const auto *qsg_857 = buffer.data(qsg + 857);
    const auto *qsg_858 = buffer.data(qsg + 858);
    const auto *qsg_860 = buffer.data(qsg + 860);
    const auto *qsg_861 = buffer.data(qsg + 861);
    const auto *qsg_864 = buffer.data(qsg + 864);
    const auto *qsg_865 = buffer.data(qsg + 865);
    const auto *qsg_866 = buffer.data(qsg + 866);
    const auto *qsg_867 = buffer.data(qsg + 867);
    const auto *qsg_868 = buffer.data(qsg + 868);
    const auto *qsg_869 = buffer.data(qsg + 869);
    const auto *qsg_870 = buffer.data(qsg + 870);
    const auto *qsg_872 = buffer.data(qsg + 872);
    const auto *qsg_873 = buffer.data(qsg + 873);
    const auto *qsg_875 = buffer.data(qsg + 875);
    const auto *qsg_876 = buffer.data(qsg + 876);
    const auto *qsg_879 = buffer.data(qsg + 879);
    const auto *qsg_880 = buffer.data(qsg + 880);
    const auto *qsg_881 = buffer.data(qsg + 881);
    const auto *qsg_882 = buffer.data(qsg + 882);
    const auto *qsg_883 = buffer.data(qsg + 883);
    const auto *qsg_884 = buffer.data(qsg + 884);
    const auto *qsg_885 = buffer.data(qsg + 885);
    const auto *qsg_887 = buffer.data(qsg + 887);
    const auto *qsg_888 = buffer.data(qsg + 888);
    const auto *qsg_890 = buffer.data(qsg + 890);
    const auto *qsg_891 = buffer.data(qsg + 891);
    const auto *qsg_894 = buffer.data(qsg + 894);
    const auto *qsg_895 = buffer.data(qsg + 895);
    const auto *qsg_896 = buffer.data(qsg + 896);
    const auto *qsg_897 = buffer.data(qsg + 897);
    const auto *qsg_898 = buffer.data(qsg + 898);
    const auto *qsg_899 = buffer.data(qsg + 899);
    const auto *qsg_900 = buffer.data(qsg + 900);
    const auto *qsg_902 = buffer.data(qsg + 902);
    const auto *qsg_903 = buffer.data(qsg + 903);
    const auto *qsg_905 = buffer.data(qsg + 905);
    const auto *qsg_906 = buffer.data(qsg + 906);
    const auto *qsg_909 = buffer.data(qsg + 909);
    const auto *qsg_910 = buffer.data(qsg + 910);
    const auto *qsg_911 = buffer.data(qsg + 911);
    const auto *qsg_912 = buffer.data(qsg + 912);
    const auto *qsg_913 = buffer.data(qsg + 913);
    const auto *qsg_914 = buffer.data(qsg + 914);
    const auto *qsg_915 = buffer.data(qsg + 915);

#pragma omp simd aligned(t_1175, t_1176, t_1177, t_1178, pa_z, pc_y, pc_z, osh0_945, osg_675, \
                         osg_690, osh1_945, qsf0_559, qsf1_559, qsg_839, \
                         qsg_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = f_1 * qsf0_559[k]
                    - f_2 * qsf1_559[k]
                    + f_3 * pc_z[k] * qsg_839[k];

        t_1176[k] = pa_z[k] * osh0_945[k]
                    - f_8 * pc_z[k] * osh1_945[k];

        t_1177[k] = f_16 * osg_690[k]
                    + f_3 * pc_y[k] * qsg_840[k];

        t_1178[k] = f_9 * osg_675[k]
                    + f_3 * pc_z[k] * qsg_840[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pa_z, pc_x, pc_y, pc_z, osh0_948, osg_692, \
                         osg_845, osh1_948, qsf0_565, qsf1_565, qsg_842, \
                         qsg_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = pa_z[k] * osh0_948[k]
                    - f_8 * pc_z[k] * osh1_948[k];

        t_1180[k] = f_16 * osg_692[k]
                    + f_3 * pc_y[k] * qsg_842[k];

        t_1181[k] = f_10 * osg_845[k]
                    + f_6 * qsf0_565[k]
                    - f_7 * qsf1_565[k]
                    + f_3 * pc_x[k] * qsg_845[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pa_z, pc_y, pc_z, osh0_951, osg_678, osg_695, \
                         osh1_951, qsg_843, qsg_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = pa_z[k] * osh0_951[k]
                    - f_8 * pc_z[k] * osh1_951[k];

        t_1183[k] = f_9 * osg_678[k]
                    + f_3 * pc_z[k] * qsg_843[k];

        t_1184[k] = f_16 * osg_695[k]
                    + f_3 * pc_y[k] * qsg_845[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, pc_x, osg_849, osg_850, osg_851, \
                         osg_852, qsf0_569, qsf1_569, qsg_849, qsg_850, qsg_851, \
                         qsg_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_10 * osg_849[k]
                    + f_4 * qsf0_569[k]
                    - f_5 * qsf1_569[k]
                    + f_3 * pc_x[k] * qsg_849[k];

        t_1186[k] = f_10 * osg_850[k]
                    + f_3 * pc_x[k] * qsg_850[k];

        t_1187[k] = f_10 * osg_851[k]
                    + f_3 * pc_x[k] * qsg_851[k];

        t_1188[k] = f_10 * osg_852[k]
                    + f_3 * pc_x[k] * qsg_852[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_z, pc_x, pc_z, osh0_960, osg_685, \
                         osg_853, osg_854, osh1_960, qsg_850, qsg_853, \
                         qsg_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_10 * osg_853[k]
                    + f_3 * pc_x[k] * qsg_853[k];

        t_1190[k] = f_10 * osg_854[k]
                    + f_3 * pc_x[k] * qsg_854[k];

        t_1191[k] = pa_z[k] * osh0_960[k]
                    - f_8 * pc_z[k] * osh1_960[k];

        t_1192[k] = f_9 * osg_685[k]
                    + f_3 * pc_z[k] * qsg_850[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pc_y, osg_702, osg_703, osg_704, qsf0_568, \
                         qsf0_569, qsf1_568, qsf1_569, qsg_852, qsg_853, \
                         qsg_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_16 * osg_702[k]
                    + f_6 * qsf0_568[k]
                    - f_7 * qsf1_568[k]
                    + f_3 * pc_y[k] * qsg_852[k];

        t_1194[k] = f_16 * osg_703[k]
                    + f_4 * qsf0_569[k]
                    - f_5 * qsf1_569[k]
                    + f_3 * pc_y[k] * qsg_853[k];

        t_1195[k] = f_16 * osg_704[k]
                    + f_3 * pc_y[k] * qsg_854[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, osg_689, osg_705, osg_855, \
                         qsf0_569, qsf0_570, qsf1_569, qsf1_570, qsg_854, \
                         qsg_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_9 * osg_689[k]
                    + f_1 * qsf0_569[k]
                    - f_2 * qsf1_569[k]
                    + f_3 * pc_z[k] * qsg_854[k];

        t_1197[k] = f_10 * osg_855[k]
                    + f_1 * qsf0_570[k]
                    - f_2 * qsf1_570[k]
                    + f_3 * pc_x[k] * qsg_855[k];

        t_1198[k] = f_17 * osg_705[k]
                    + f_3 * pc_y[k] * qsg_855[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_y, pc_z, osg_690, osg_707, osg_858, \
                         qsf0_573, qsf1_573, qsg_855, qsg_857, \
                         qsg_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_10 * osg_690[k]
                    + f_3 * pc_z[k] * qsg_855[k];

        t_1200[k] = f_10 * osg_858[k]
                    + f_6 * qsf0_573[k]
                    - f_7 * qsf1_573[k]
                    + f_3 * pc_x[k] * qsg_858[k];

        t_1201[k] = f_17 * osg_707[k]
                    + f_3 * pc_y[k] * qsg_857[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, pc_z, osg_693, osg_860, osg_861, \
                         qsf0_575, qsf0_576, qsf1_575, qsf1_576, qsg_858, qsg_860, \
                         qsg_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_10 * osg_860[k]
                    + f_6 * qsf0_575[k]
                    - f_7 * qsf1_575[k]
                    + f_3 * pc_x[k] * qsg_860[k];

        t_1203[k] = f_10 * osg_861[k]
                    + f_4 * qsf0_576[k]
                    - f_5 * qsf1_576[k]
                    + f_3 * pc_x[k] * qsg_861[k];

        t_1204[k] = f_10 * osg_693[k]
                    + f_3 * pc_z[k] * qsg_858[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, pc_x, pc_y, osg_710, osg_864, \
                         osg_865, osg_866, qsf0_579, qsf1_579, qsg_860, qsg_864, qsg_865, \
                         qsg_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_17 * osg_710[k]
                    + f_3 * pc_y[k] * qsg_860[k];

        t_1206[k] = f_10 * osg_864[k]
                    + f_4 * qsf0_579[k]
                    - f_5 * qsf1_579[k]
                    + f_3 * pc_x[k] * qsg_864[k];

        t_1207[k] = f_10 * osg_865[k]
                    + f_3 * pc_x[k] * qsg_865[k];

        t_1208[k] = f_10 * osg_866[k]
                    + f_3 * pc_x[k] * qsg_866[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, t_1212, pc_x, pc_y, osg_715, osg_867, \
                         osg_868, osg_869, qsf0_576, qsf1_576, qsg_865, qsg_867, qsg_868, \
                         qsg_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = f_10 * osg_867[k]
                    + f_3 * pc_x[k] * qsg_867[k];

        t_1210[k] = f_10 * osg_868[k]
                    + f_3 * pc_x[k] * qsg_868[k];

        t_1211[k] = f_10 * osg_869[k]
                    + f_3 * pc_x[k] * qsg_869[k];

        t_1212[k] = f_17 * osg_715[k]
                    + f_1 * qsf0_576[k]
                    - f_2 * qsf1_576[k]
                    + f_3 * pc_y[k] * qsg_865[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, pc_y, pc_z, osg_700, osg_717, osg_718, \
                         qsf0_578, qsf0_579, qsf1_578, qsf1_579, qsg_865, qsg_867, \
                         qsg_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = f_10 * osg_700[k]
                    + f_3 * pc_z[k] * qsg_865[k];

        t_1214[k] = f_17 * osg_717[k]
                    + f_6 * qsf0_578[k]
                    - f_7 * qsf1_578[k]
                    + f_3 * pc_y[k] * qsg_867[k];

        t_1215[k] = f_17 * osg_718[k]
                    + f_4 * qsf0_579[k]
                    - f_5 * qsf1_579[k]
                    + f_3 * pc_y[k] * qsg_868[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_x, pc_y, pc_z, osg_704, osg_719, osg_870, \
                         qsf0_579, qsf0_580, qsf1_579, qsf1_580, qsg_869, \
                         qsg_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_17 * osg_719[k]
                    + f_3 * pc_y[k] * qsg_869[k];

        t_1217[k] = f_10 * osg_704[k]
                    + f_1 * qsf0_579[k]
                    - f_2 * qsf1_579[k]
                    + f_3 * pc_z[k] * qsg_869[k];

        t_1218[k] = f_10 * osg_870[k]
                    + f_1 * qsf0_580[k]
                    - f_2 * qsf1_580[k]
                    + f_3 * pc_x[k] * qsg_870[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, t_1222, pc_x, pc_y, pc_z, osg_705, osg_720, \
                         osg_722, osg_873, qsf0_583, qsf1_583, qsg_870, qsg_872, \
                         qsg_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_19 * osg_720[k]
                    + f_3 * pc_y[k] * qsg_870[k];

        t_1220[k] = f_11 * osg_705[k]
                    + f_3 * pc_z[k] * qsg_870[k];

        t_1221[k] = f_10 * osg_873[k]
                    + f_6 * qsf0_583[k]
                    - f_7 * qsf1_583[k]
                    + f_3 * pc_x[k] * qsg_873[k];

        t_1222[k] = f_19 * osg_722[k]
                    + f_3 * pc_y[k] * qsg_872[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, pc_x, pc_z, osg_708, osg_875, osg_876, \
                         qsf0_585, qsf0_586, qsf1_585, qsf1_586, qsg_873, qsg_875, \
                         qsg_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = f_10 * osg_875[k]
                    + f_6 * qsf0_585[k]
                    - f_7 * qsf1_585[k]
                    + f_3 * pc_x[k] * qsg_875[k];

        t_1224[k] = f_10 * osg_876[k]
                    + f_4 * qsf0_586[k]
                    - f_5 * qsf1_586[k]
                    + f_3 * pc_x[k] * qsg_876[k];

        t_1225[k] = f_11 * osg_708[k]
                    + f_3 * pc_z[k] * qsg_873[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pc_x, pc_y, osg_725, osg_879, \
                         osg_880, osg_881, qsf0_589, qsf1_589, qsg_875, qsg_879, qsg_880, \
                         qsg_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_19 * osg_725[k]
                    + f_3 * pc_y[k] * qsg_875[k];

        t_1227[k] = f_10 * osg_879[k]
                    + f_4 * qsf0_589[k]
                    - f_5 * qsf1_589[k]
                    + f_3 * pc_x[k] * qsg_879[k];

        t_1228[k] = f_10 * osg_880[k]
                    + f_3 * pc_x[k] * qsg_880[k];

        t_1229[k] = f_10 * osg_881[k]
                    + f_3 * pc_x[k] * qsg_881[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pc_x, pc_y, osg_730, osg_882, \
                         osg_883, osg_884, qsf0_586, qsf1_586, qsg_880, qsg_882, qsg_883, \
                         qsg_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_10 * osg_882[k]
                    + f_3 * pc_x[k] * qsg_882[k];

        t_1231[k] = f_10 * osg_883[k]
                    + f_3 * pc_x[k] * qsg_883[k];

        t_1232[k] = f_10 * osg_884[k]
                    + f_3 * pc_x[k] * qsg_884[k];

        t_1233[k] = f_19 * osg_730[k]
                    + f_1 * qsf0_586[k]
                    - f_2 * qsf1_586[k]
                    + f_3 * pc_y[k] * qsg_880[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pc_y, pc_z, osg_715, osg_732, osg_733, \
                         qsf0_588, qsf0_589, qsf1_588, qsf1_589, qsg_880, qsg_882, \
                         qsg_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = f_11 * osg_715[k]
                    + f_3 * pc_z[k] * qsg_880[k];

        t_1235[k] = f_19 * osg_732[k]
                    + f_6 * qsf0_588[k]
                    - f_7 * qsf1_588[k]
                    + f_3 * pc_y[k] * qsg_882[k];

        t_1236[k] = f_19 * osg_733[k]
                    + f_4 * qsf0_589[k]
                    - f_5 * qsf1_589[k]
                    + f_3 * pc_y[k] * qsg_883[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, pc_x, pc_y, pc_z, osg_719, osg_734, osg_885, \
                         qsf0_589, qsf0_590, qsf1_589, qsf1_590, qsg_884, \
                         qsg_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_19 * osg_734[k]
                    + f_3 * pc_y[k] * qsg_884[k];

        t_1238[k] = f_11 * osg_719[k]
                    + f_1 * qsf0_589[k]
                    - f_2 * qsf1_589[k]
                    + f_3 * pc_z[k] * qsg_884[k];

        t_1239[k] = f_10 * osg_885[k]
                    + f_1 * qsf0_590[k]
                    - f_2 * qsf1_590[k]
                    + f_3 * pc_x[k] * qsg_885[k];
    }

#pragma omp simd aligned(t_1240, t_1241, t_1242, t_1243, pc_x, pc_y, pc_z, osg_720, osg_735, \
                         osg_737, osg_888, qsf0_593, qsf1_593, qsg_885, qsg_887, \
                         qsg_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1240[k] = f_21 * osg_735[k]
                    + f_3 * pc_y[k] * qsg_885[k];

        t_1241[k] = f_18 * osg_720[k]
                    + f_3 * pc_z[k] * qsg_885[k];

        t_1242[k] = f_10 * osg_888[k]
                    + f_6 * qsf0_593[k]
                    - f_7 * qsf1_593[k]
                    + f_3 * pc_x[k] * qsg_888[k];

        t_1243[k] = f_21 * osg_737[k]
                    + f_3 * pc_y[k] * qsg_887[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, pc_x, pc_z, osg_723, osg_890, osg_891, \
                         qsf0_595, qsf0_596, qsf1_595, qsf1_596, qsg_888, qsg_890, \
                         qsg_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_10 * osg_890[k]
                    + f_6 * qsf0_595[k]
                    - f_7 * qsf1_595[k]
                    + f_3 * pc_x[k] * qsg_890[k];

        t_1245[k] = f_10 * osg_891[k]
                    + f_4 * qsf0_596[k]
                    - f_5 * qsf1_596[k]
                    + f_3 * pc_x[k] * qsg_891[k];

        t_1246[k] = f_18 * osg_723[k]
                    + f_3 * pc_z[k] * qsg_888[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, pc_x, pc_y, osg_740, osg_894, \
                         osg_895, osg_896, qsf0_599, qsf1_599, qsg_890, qsg_894, qsg_895, \
                         qsg_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_21 * osg_740[k]
                    + f_3 * pc_y[k] * qsg_890[k];

        t_1248[k] = f_10 * osg_894[k]
                    + f_4 * qsf0_599[k]
                    - f_5 * qsf1_599[k]
                    + f_3 * pc_x[k] * qsg_894[k];

        t_1249[k] = f_10 * osg_895[k]
                    + f_3 * pc_x[k] * qsg_895[k];

        t_1250[k] = f_10 * osg_896[k]
                    + f_3 * pc_x[k] * qsg_896[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, pc_x, pc_y, osg_745, osg_897, \
                         osg_898, osg_899, qsf0_596, qsf1_596, qsg_895, qsg_897, qsg_898, \
                         qsg_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_10 * osg_897[k]
                    + f_3 * pc_x[k] * qsg_897[k];

        t_1252[k] = f_10 * osg_898[k]
                    + f_3 * pc_x[k] * qsg_898[k];

        t_1253[k] = f_10 * osg_899[k]
                    + f_3 * pc_x[k] * qsg_899[k];

        t_1254[k] = f_21 * osg_745[k]
                    + f_1 * qsf0_596[k]
                    - f_2 * qsf1_596[k]
                    + f_3 * pc_y[k] * qsg_895[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, pc_y, pc_z, osg_730, osg_747, osg_748, \
                         qsf0_598, qsf0_599, qsf1_598, qsf1_599, qsg_895, qsg_897, \
                         qsg_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = f_18 * osg_730[k]
                    + f_3 * pc_z[k] * qsg_895[k];

        t_1256[k] = f_21 * osg_747[k]
                    + f_6 * qsf0_598[k]
                    - f_7 * qsf1_598[k]
                    + f_3 * pc_y[k] * qsg_897[k];

        t_1257[k] = f_21 * osg_748[k]
                    + f_4 * qsf0_599[k]
                    - f_5 * qsf1_599[k]
                    + f_3 * pc_y[k] * qsg_898[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, pc_x, pc_y, pc_z, osg_734, osg_749, osg_900, \
                         qsf0_599, qsf0_600, qsf1_599, qsf1_600, qsg_899, \
                         qsg_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = f_21 * osg_749[k]
                    + f_3 * pc_y[k] * qsg_899[k];

        t_1259[k] = f_18 * osg_734[k]
                    + f_1 * qsf0_599[k]
                    - f_2 * qsf1_599[k]
                    + f_3 * pc_z[k] * qsg_899[k];

        t_1260[k] = f_10 * osg_900[k]
                    + f_1 * qsf0_600[k]
                    - f_2 * qsf1_600[k]
                    + f_3 * pc_x[k] * qsg_900[k];
    }

#pragma omp simd aligned(t_1261, t_1262, t_1263, t_1264, pc_x, pc_y, pc_z, osg_735, osg_750, \
                         osg_752, osg_903, qsf0_603, qsf1_603, qsg_900, qsg_902, \
                         qsg_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1261[k] = f_20 * osg_750[k]
                    + f_3 * pc_y[k] * qsg_900[k];

        t_1262[k] = f_20 * osg_735[k]
                    + f_3 * pc_z[k] * qsg_900[k];

        t_1263[k] = f_10 * osg_903[k]
                    + f_6 * qsf0_603[k]
                    - f_7 * qsf1_603[k]
                    + f_3 * pc_x[k] * qsg_903[k];

        t_1264[k] = f_20 * osg_752[k]
                    + f_3 * pc_y[k] * qsg_902[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, pc_x, pc_z, osg_738, osg_905, osg_906, \
                         qsf0_605, qsf0_606, qsf1_605, qsf1_606, qsg_903, qsg_905, \
                         qsg_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_10 * osg_905[k]
                    + f_6 * qsf0_605[k]
                    - f_7 * qsf1_605[k]
                    + f_3 * pc_x[k] * qsg_905[k];

        t_1266[k] = f_10 * osg_906[k]
                    + f_4 * qsf0_606[k]
                    - f_5 * qsf1_606[k]
                    + f_3 * pc_x[k] * qsg_906[k];

        t_1267[k] = f_20 * osg_738[k]
                    + f_3 * pc_z[k] * qsg_903[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pc_x, pc_y, osg_755, osg_909, \
                         osg_910, osg_911, qsf0_609, qsf1_609, qsg_905, qsg_909, qsg_910, \
                         qsg_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_20 * osg_755[k]
                    + f_3 * pc_y[k] * qsg_905[k];

        t_1269[k] = f_10 * osg_909[k]
                    + f_4 * qsf0_609[k]
                    - f_5 * qsf1_609[k]
                    + f_3 * pc_x[k] * qsg_909[k];

        t_1270[k] = f_10 * osg_910[k]
                    + f_3 * pc_x[k] * qsg_910[k];

        t_1271[k] = f_10 * osg_911[k]
                    + f_3 * pc_x[k] * qsg_911[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, osg_760, osg_912, \
                         osg_913, osg_914, qsf0_606, qsf1_606, qsg_910, qsg_912, qsg_913, \
                         qsg_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_10 * osg_912[k]
                    + f_3 * pc_x[k] * qsg_912[k];

        t_1273[k] = f_10 * osg_913[k]
                    + f_3 * pc_x[k] * qsg_913[k];

        t_1274[k] = f_10 * osg_914[k]
                    + f_3 * pc_x[k] * qsg_914[k];

        t_1275[k] = f_20 * osg_760[k]
                    + f_1 * qsf0_606[k]
                    - f_2 * qsf1_606[k]
                    + f_3 * pc_y[k] * qsg_910[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, pc_y, pc_z, osg_745, osg_762, osg_763, \
                         qsf0_608, qsf0_609, qsf1_608, qsf1_609, qsg_910, qsg_912, \
                         qsg_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_20 * osg_745[k]
                    + f_3 * pc_z[k] * qsg_910[k];

        t_1277[k] = f_20 * osg_762[k]
                    + f_6 * qsf0_608[k]
                    - f_7 * qsf1_608[k]
                    + f_3 * pc_y[k] * qsg_912[k];

        t_1278[k] = f_20 * osg_763[k]
                    + f_4 * qsf0_609[k]
                    - f_5 * qsf1_609[k]
                    + f_3 * pc_y[k] * qsg_913[k];
    }

#pragma omp simd aligned(t_1279, t_1280, t_1281, pc_x, pc_y, pc_z, osg_749, osg_764, osg_915, \
                         qsf0_609, qsf0_610, qsf1_609, qsf1_610, qsg_914, \
                         qsg_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1279[k] = f_20 * osg_764[k]
                    + f_3 * pc_y[k] * qsg_914[k];

        t_1280[k] = f_20 * osg_749[k]
                    + f_1 * qsf0_609[k]
                    - f_2 * qsf1_609[k]
                    + f_3 * pc_z[k] * qsg_914[k];

        t_1281[k] = f_10 * osg_915[k]
                    + f_1 * qsf0_610[k]
                    - f_2 * qsf1_610[k]
                    + f_3 * pc_x[k] * qsg_915[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsf0, const size_t qsf1,
                                                           const size_t qsg, const size_t ncols,
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
    const auto f_12 = 5.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_1134 = buffer.data(osh0 + 1134);
    const auto *osh0_1137 = buffer.data(osh0 + 1137);
    const auto *osh0_1139 = buffer.data(osh0 + 1139);
    const auto *osh0_1140 = buffer.data(osh0 + 1140);
    const auto *osh0_1143 = buffer.data(osh0 + 1143);
    const auto *osh0_1154 = buffer.data(osh0 + 1154);
    const auto *osh0_1386 = buffer.data(osh0 + 1386);
    const auto *osh0_1389 = buffer.data(osh0 + 1389);
    const auto *osh0_1392 = buffer.data(osh0 + 1392);

    const auto *osg_750 = buffer.data(osg + 750);
    const auto *osg_753 = buffer.data(osg + 753);
    const auto *osg_760 = buffer.data(osg + 760);
    const auto *osg_764 = buffer.data(osg + 764);
    const auto *osg_765 = buffer.data(osg + 765);
    const auto *osg_767 = buffer.data(osg + 767);
    const auto *osg_768 = buffer.data(osg + 768);
    const auto *osg_770 = buffer.data(osg + 770);
    const auto *osg_775 = buffer.data(osg + 775);
    const auto *osg_777 = buffer.data(osg + 777);
    const auto *osg_778 = buffer.data(osg + 778);
    const auto *osg_779 = buffer.data(osg + 779);
    const auto *osg_780 = buffer.data(osg + 780);
    const auto *osg_782 = buffer.data(osg + 782);
    const auto *osg_783 = buffer.data(osg + 783);
    const auto *osg_785 = buffer.data(osg + 785);
    const auto *osg_790 = buffer.data(osg + 790);
    const auto *osg_792 = buffer.data(osg + 792);
    const auto *osg_793 = buffer.data(osg + 793);
    const auto *osg_794 = buffer.data(osg + 794);
    const auto *osg_795 = buffer.data(osg + 795);
    const auto *osg_797 = buffer.data(osg + 797);
    const auto *osg_798 = buffer.data(osg + 798);
    const auto *osg_800 = buffer.data(osg + 800);
    const auto *osg_805 = buffer.data(osg + 805);
    const auto *osg_807 = buffer.data(osg + 807);
    const auto *osg_808 = buffer.data(osg + 808);
    const auto *osg_809 = buffer.data(osg + 809);
    const auto *osg_810 = buffer.data(osg + 810);
    const auto *osg_811 = buffer.data(osg + 811);
    const auto *osg_812 = buffer.data(osg + 812);
    const auto *osg_813 = buffer.data(osg + 813);
    const auto *osg_815 = buffer.data(osg + 815);
    const auto *osg_820 = buffer.data(osg + 820);
    const auto *osg_822 = buffer.data(osg + 822);
    const auto *osg_823 = buffer.data(osg + 823);
    const auto *osg_824 = buffer.data(osg + 824);
    const auto *osg_825 = buffer.data(osg + 825);
    const auto *osg_830 = buffer.data(osg + 830);
    const auto *osg_918 = buffer.data(osg + 918);
    const auto *osg_920 = buffer.data(osg + 920);
    const auto *osg_921 = buffer.data(osg + 921);
    const auto *osg_924 = buffer.data(osg + 924);
    const auto *osg_925 = buffer.data(osg + 925);
    const auto *osg_926 = buffer.data(osg + 926);
    const auto *osg_927 = buffer.data(osg + 927);
    const auto *osg_928 = buffer.data(osg + 928);
    const auto *osg_929 = buffer.data(osg + 929);
    const auto *osg_930 = buffer.data(osg + 930);
    const auto *osg_933 = buffer.data(osg + 933);
    const auto *osg_935 = buffer.data(osg + 935);
    const auto *osg_936 = buffer.data(osg + 936);
    const auto *osg_939 = buffer.data(osg + 939);
    const auto *osg_940 = buffer.data(osg + 940);
    const auto *osg_941 = buffer.data(osg + 941);
    const auto *osg_942 = buffer.data(osg + 942);
    const auto *osg_943 = buffer.data(osg + 943);
    const auto *osg_944 = buffer.data(osg + 944);
    const auto *osg_945 = buffer.data(osg + 945);
    const auto *osg_948 = buffer.data(osg + 948);
    const auto *osg_950 = buffer.data(osg + 950);
    const auto *osg_951 = buffer.data(osg + 951);
    const auto *osg_954 = buffer.data(osg + 954);
    const auto *osg_955 = buffer.data(osg + 955);
    const auto *osg_956 = buffer.data(osg + 956);
    const auto *osg_957 = buffer.data(osg + 957);
    const auto *osg_958 = buffer.data(osg + 958);
    const auto *osg_959 = buffer.data(osg + 959);
    const auto *osg_970 = buffer.data(osg + 970);
    const auto *osg_971 = buffer.data(osg + 971);
    const auto *osg_972 = buffer.data(osg + 972);
    const auto *osg_973 = buffer.data(osg + 973);
    const auto *osg_974 = buffer.data(osg + 974);
    const auto *osg_975 = buffer.data(osg + 975);
    const auto *osg_980 = buffer.data(osg + 980);
    const auto *osg_984 = buffer.data(osg + 984);
    const auto *osg_985 = buffer.data(osg + 985);
    const auto *osg_986 = buffer.data(osg + 986);
    const auto *osg_987 = buffer.data(osg + 987);
    const auto *osg_989 = buffer.data(osg + 989);
    const auto *osg_990 = buffer.data(osg + 990);
    const auto *osg_993 = buffer.data(osg + 993);
    const auto *osg_996 = buffer.data(osg + 996);
    const auto *osg_1000 = buffer.data(osg + 1000);

    const auto *osh1_1134 = buffer.data(osh1 + 1134);
    const auto *osh1_1137 = buffer.data(osh1 + 1137);
    const auto *osh1_1139 = buffer.data(osh1 + 1139);
    const auto *osh1_1140 = buffer.data(osh1 + 1140);
    const auto *osh1_1143 = buffer.data(osh1 + 1143);
    const auto *osh1_1154 = buffer.data(osh1 + 1154);
    const auto *osh1_1386 = buffer.data(osh1 + 1386);
    const auto *osh1_1389 = buffer.data(osh1 + 1389);
    const auto *osh1_1392 = buffer.data(osh1 + 1392);

    const auto *qsf0_613 = buffer.data(qsf0 + 613);
    const auto *qsf0_615 = buffer.data(qsf0 + 615);
    const auto *qsf0_616 = buffer.data(qsf0 + 616);
    const auto *qsf0_618 = buffer.data(qsf0 + 618);
    const auto *qsf0_619 = buffer.data(qsf0 + 619);
    const auto *qsf0_620 = buffer.data(qsf0 + 620);
    const auto *qsf0_623 = buffer.data(qsf0 + 623);
    const auto *qsf0_625 = buffer.data(qsf0 + 625);
    const auto *qsf0_626 = buffer.data(qsf0 + 626);
    const auto *qsf0_628 = buffer.data(qsf0 + 628);
    const auto *qsf0_629 = buffer.data(qsf0 + 629);
    const auto *qsf0_630 = buffer.data(qsf0 + 630);
    const auto *qsf0_633 = buffer.data(qsf0 + 633);
    const auto *qsf0_635 = buffer.data(qsf0 + 635);
    const auto *qsf0_636 = buffer.data(qsf0 + 636);
    const auto *qsf0_638 = buffer.data(qsf0 + 638);
    const auto *qsf0_639 = buffer.data(qsf0 + 639);
    const auto *qsf0_646 = buffer.data(qsf0 + 646);
    const auto *qsf0_648 = buffer.data(qsf0 + 648);
    const auto *qsf0_649 = buffer.data(qsf0 + 649);
    const auto *qsf0_650 = buffer.data(qsf0 + 650);
    const auto *qsf0_651 = buffer.data(qsf0 + 651);
    const auto *qsf0_652 = buffer.data(qsf0 + 652);
    const auto *qsf0_655 = buffer.data(qsf0 + 655);
    const auto *qsf0_656 = buffer.data(qsf0 + 656);
    const auto *qsf0_657 = buffer.data(qsf0 + 657);
    const auto *qsf0_658 = buffer.data(qsf0 + 658);
    const auto *qsf0_659 = buffer.data(qsf0 + 659);
    const auto *qsf0_660 = buffer.data(qsf0 + 660);
    const auto *qsf0_662 = buffer.data(qsf0 + 662);

    const auto *qsf1_613 = buffer.data(qsf1 + 613);
    const auto *qsf1_615 = buffer.data(qsf1 + 615);
    const auto *qsf1_616 = buffer.data(qsf1 + 616);
    const auto *qsf1_618 = buffer.data(qsf1 + 618);
    const auto *qsf1_619 = buffer.data(qsf1 + 619);
    const auto *qsf1_620 = buffer.data(qsf1 + 620);
    const auto *qsf1_623 = buffer.data(qsf1 + 623);
    const auto *qsf1_625 = buffer.data(qsf1 + 625);
    const auto *qsf1_626 = buffer.data(qsf1 + 626);
    const auto *qsf1_628 = buffer.data(qsf1 + 628);
    const auto *qsf1_629 = buffer.data(qsf1 + 629);
    const auto *qsf1_630 = buffer.data(qsf1 + 630);
    const auto *qsf1_633 = buffer.data(qsf1 + 633);
    const auto *qsf1_635 = buffer.data(qsf1 + 635);
    const auto *qsf1_636 = buffer.data(qsf1 + 636);
    const auto *qsf1_638 = buffer.data(qsf1 + 638);
    const auto *qsf1_639 = buffer.data(qsf1 + 639);
    const auto *qsf1_646 = buffer.data(qsf1 + 646);
    const auto *qsf1_648 = buffer.data(qsf1 + 648);
    const auto *qsf1_649 = buffer.data(qsf1 + 649);
    const auto *qsf1_650 = buffer.data(qsf1 + 650);
    const auto *qsf1_651 = buffer.data(qsf1 + 651);
    const auto *qsf1_652 = buffer.data(qsf1 + 652);
    const auto *qsf1_655 = buffer.data(qsf1 + 655);
    const auto *qsf1_656 = buffer.data(qsf1 + 656);
    const auto *qsf1_657 = buffer.data(qsf1 + 657);
    const auto *qsf1_658 = buffer.data(qsf1 + 658);
    const auto *qsf1_659 = buffer.data(qsf1 + 659);
    const auto *qsf1_660 = buffer.data(qsf1 + 660);
    const auto *qsf1_662 = buffer.data(qsf1 + 662);

    const auto *qsg_915 = buffer.data(qsg + 915);
    const auto *qsg_917 = buffer.data(qsg + 917);
    const auto *qsg_918 = buffer.data(qsg + 918);
    const auto *qsg_920 = buffer.data(qsg + 920);
    const auto *qsg_921 = buffer.data(qsg + 921);
    const auto *qsg_924 = buffer.data(qsg + 924);
    const auto *qsg_925 = buffer.data(qsg + 925);
    const auto *qsg_926 = buffer.data(qsg + 926);
    const auto *qsg_927 = buffer.data(qsg + 927);
    const auto *qsg_928 = buffer.data(qsg + 928);
    const auto *qsg_929 = buffer.data(qsg + 929);
    const auto *qsg_930 = buffer.data(qsg + 930);
    const auto *qsg_932 = buffer.data(qsg + 932);
    const auto *qsg_933 = buffer.data(qsg + 933);
    const auto *qsg_935 = buffer.data(qsg + 935);
    const auto *qsg_936 = buffer.data(qsg + 936);
    const auto *qsg_939 = buffer.data(qsg + 939);
    const auto *qsg_940 = buffer.data(qsg + 940);
    const auto *qsg_941 = buffer.data(qsg + 941);
    const auto *qsg_942 = buffer.data(qsg + 942);
    const auto *qsg_943 = buffer.data(qsg + 943);
    const auto *qsg_944 = buffer.data(qsg + 944);
    const auto *qsg_945 = buffer.data(qsg + 945);
    const auto *qsg_947 = buffer.data(qsg + 947);
    const auto *qsg_948 = buffer.data(qsg + 948);
    const auto *qsg_950 = buffer.data(qsg + 950);
    const auto *qsg_951 = buffer.data(qsg + 951);
    const auto *qsg_954 = buffer.data(qsg + 954);
    const auto *qsg_955 = buffer.data(qsg + 955);
    const auto *qsg_956 = buffer.data(qsg + 956);
    const auto *qsg_957 = buffer.data(qsg + 957);
    const auto *qsg_958 = buffer.data(qsg + 958);
    const auto *qsg_959 = buffer.data(qsg + 959);
    const auto *qsg_960 = buffer.data(qsg + 960);
    const auto *qsg_962 = buffer.data(qsg + 962);
    const auto *qsg_963 = buffer.data(qsg + 963);
    const auto *qsg_965 = buffer.data(qsg + 965);
    const auto *qsg_970 = buffer.data(qsg + 970);
    const auto *qsg_971 = buffer.data(qsg + 971);
    const auto *qsg_972 = buffer.data(qsg + 972);
    const auto *qsg_973 = buffer.data(qsg + 973);
    const auto *qsg_974 = buffer.data(qsg + 974);
    const auto *qsg_975 = buffer.data(qsg + 975);
    const auto *qsg_976 = buffer.data(qsg + 976);
    const auto *qsg_977 = buffer.data(qsg + 977);
    const auto *qsg_978 = buffer.data(qsg + 978);
    const auto *qsg_979 = buffer.data(qsg + 979);
    const auto *qsg_980 = buffer.data(qsg + 980);
    const auto *qsg_984 = buffer.data(qsg + 984);
    const auto *qsg_985 = buffer.data(qsg + 985);
    const auto *qsg_986 = buffer.data(qsg + 986);
    const auto *qsg_987 = buffer.data(qsg + 987);
    const auto *qsg_988 = buffer.data(qsg + 988);
    const auto *qsg_989 = buffer.data(qsg + 989);
    const auto *qsg_990 = buffer.data(qsg + 990);
    const auto *qsg_991 = buffer.data(qsg + 991);
    const auto *qsg_992 = buffer.data(qsg + 992);
    const auto *qsg_993 = buffer.data(qsg + 993);
    const auto *qsg_995 = buffer.data(qsg + 995);
    const auto *qsg_996 = buffer.data(qsg + 996);
    const auto *qsg_1000 = buffer.data(qsg + 1000);

#pragma omp simd aligned(t_1282, t_1283, t_1284, t_1285, pc_x, pc_y, pc_z, osg_750, osg_765, \
                         osg_767, osg_918, qsf0_613, qsf1_613, qsg_915, qsg_917, \
                         qsg_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1282[k] = f_18 * osg_765[k]
                    + f_3 * pc_y[k] * qsg_915[k];

        t_1283[k] = f_21 * osg_750[k]
                    + f_3 * pc_z[k] * qsg_915[k];

        t_1284[k] = f_10 * osg_918[k]
                    + f_6 * qsf0_613[k]
                    - f_7 * qsf1_613[k]
                    + f_3 * pc_x[k] * qsg_918[k];

        t_1285[k] = f_18 * osg_767[k]
                    + f_3 * pc_y[k] * qsg_917[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, pc_x, pc_z, osg_753, osg_920, osg_921, \
                         qsf0_615, qsf0_616, qsf1_615, qsf1_616, qsg_918, qsg_920, \
                         qsg_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_10 * osg_920[k]
                    + f_6 * qsf0_615[k]
                    - f_7 * qsf1_615[k]
                    + f_3 * pc_x[k] * qsg_920[k];

        t_1287[k] = f_10 * osg_921[k]
                    + f_4 * qsf0_616[k]
                    - f_5 * qsf1_616[k]
                    + f_3 * pc_x[k] * qsg_921[k];

        t_1288[k] = f_21 * osg_753[k]
                    + f_3 * pc_z[k] * qsg_918[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, pc_x, pc_y, osg_770, osg_924, \
                         osg_925, osg_926, qsf0_619, qsf1_619, qsg_920, qsg_924, qsg_925, \
                         qsg_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_18 * osg_770[k]
                    + f_3 * pc_y[k] * qsg_920[k];

        t_1290[k] = f_10 * osg_924[k]
                    + f_4 * qsf0_619[k]
                    - f_5 * qsf1_619[k]
                    + f_3 * pc_x[k] * qsg_924[k];

        t_1291[k] = f_10 * osg_925[k]
                    + f_3 * pc_x[k] * qsg_925[k];

        t_1292[k] = f_10 * osg_926[k]
                    + f_3 * pc_x[k] * qsg_926[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pc_x, pc_y, osg_775, osg_927, \
                         osg_928, osg_929, qsf0_616, qsf1_616, qsg_925, qsg_927, qsg_928, \
                         qsg_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_10 * osg_927[k]
                    + f_3 * pc_x[k] * qsg_927[k];

        t_1294[k] = f_10 * osg_928[k]
                    + f_3 * pc_x[k] * qsg_928[k];

        t_1295[k] = f_10 * osg_929[k]
                    + f_3 * pc_x[k] * qsg_929[k];

        t_1296[k] = f_18 * osg_775[k]
                    + f_1 * qsf0_616[k]
                    - f_2 * qsf1_616[k]
                    + f_3 * pc_y[k] * qsg_925[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, pc_y, pc_z, osg_760, osg_777, osg_778, \
                         qsf0_618, qsf0_619, qsf1_618, qsf1_619, qsg_925, qsg_927, \
                         qsg_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = f_21 * osg_760[k]
                    + f_3 * pc_z[k] * qsg_925[k];

        t_1298[k] = f_18 * osg_777[k]
                    + f_6 * qsf0_618[k]
                    - f_7 * qsf1_618[k]
                    + f_3 * pc_y[k] * qsg_927[k];

        t_1299[k] = f_18 * osg_778[k]
                    + f_4 * qsf0_619[k]
                    - f_5 * qsf1_619[k]
                    + f_3 * pc_y[k] * qsg_928[k];
    }

#pragma omp simd aligned(t_1300, t_1301, t_1302, pc_x, pc_y, pc_z, osg_764, osg_779, osg_930, \
                         qsf0_619, qsf0_620, qsf1_619, qsf1_620, qsg_929, \
                         qsg_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1300[k] = f_18 * osg_779[k]
                    + f_3 * pc_y[k] * qsg_929[k];

        t_1301[k] = f_21 * osg_764[k]
                    + f_1 * qsf0_619[k]
                    - f_2 * qsf1_619[k]
                    + f_3 * pc_z[k] * qsg_929[k];

        t_1302[k] = f_10 * osg_930[k]
                    + f_1 * qsf0_620[k]
                    - f_2 * qsf1_620[k]
                    + f_3 * pc_x[k] * qsg_930[k];
    }

#pragma omp simd aligned(t_1303, t_1304, t_1305, t_1306, pc_x, pc_y, pc_z, osg_765, osg_780, \
                         osg_782, osg_933, qsf0_623, qsf1_623, qsg_930, qsg_932, \
                         qsg_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1303[k] = f_11 * osg_780[k]
                    + f_3 * pc_y[k] * qsg_930[k];

        t_1304[k] = f_19 * osg_765[k]
                    + f_3 * pc_z[k] * qsg_930[k];

        t_1305[k] = f_10 * osg_933[k]
                    + f_6 * qsf0_623[k]
                    - f_7 * qsf1_623[k]
                    + f_3 * pc_x[k] * qsg_933[k];

        t_1306[k] = f_11 * osg_782[k]
                    + f_3 * pc_y[k] * qsg_932[k];
    }

#pragma omp simd aligned(t_1307, t_1308, t_1309, pc_x, pc_z, osg_768, osg_935, osg_936, \
                         qsf0_625, qsf0_626, qsf1_625, qsf1_626, qsg_933, qsg_935, \
                         qsg_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1307[k] = f_10 * osg_935[k]
                    + f_6 * qsf0_625[k]
                    - f_7 * qsf1_625[k]
                    + f_3 * pc_x[k] * qsg_935[k];

        t_1308[k] = f_10 * osg_936[k]
                    + f_4 * qsf0_626[k]
                    - f_5 * qsf1_626[k]
                    + f_3 * pc_x[k] * qsg_936[k];

        t_1309[k] = f_19 * osg_768[k]
                    + f_3 * pc_z[k] * qsg_933[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, t_1313, pc_x, pc_y, osg_785, osg_939, \
                         osg_940, osg_941, qsf0_629, qsf1_629, qsg_935, qsg_939, qsg_940, \
                         qsg_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_11 * osg_785[k]
                    + f_3 * pc_y[k] * qsg_935[k];

        t_1311[k] = f_10 * osg_939[k]
                    + f_4 * qsf0_629[k]
                    - f_5 * qsf1_629[k]
                    + f_3 * pc_x[k] * qsg_939[k];

        t_1312[k] = f_10 * osg_940[k]
                    + f_3 * pc_x[k] * qsg_940[k];

        t_1313[k] = f_10 * osg_941[k]
                    + f_3 * pc_x[k] * qsg_941[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, t_1317, pc_x, pc_y, osg_790, osg_942, \
                         osg_943, osg_944, qsf0_626, qsf1_626, qsg_940, qsg_942, qsg_943, \
                         qsg_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_10 * osg_942[k]
                    + f_3 * pc_x[k] * qsg_942[k];

        t_1315[k] = f_10 * osg_943[k]
                    + f_3 * pc_x[k] * qsg_943[k];

        t_1316[k] = f_10 * osg_944[k]
                    + f_3 * pc_x[k] * qsg_944[k];

        t_1317[k] = f_11 * osg_790[k]
                    + f_1 * qsf0_626[k]
                    - f_2 * qsf1_626[k]
                    + f_3 * pc_y[k] * qsg_940[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, pc_y, pc_z, osg_775, osg_792, osg_793, \
                         qsf0_628, qsf0_629, qsf1_628, qsf1_629, qsg_940, qsg_942, \
                         qsg_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_19 * osg_775[k]
                    + f_3 * pc_z[k] * qsg_940[k];

        t_1319[k] = f_11 * osg_792[k]
                    + f_6 * qsf0_628[k]
                    - f_7 * qsf1_628[k]
                    + f_3 * pc_y[k] * qsg_942[k];

        t_1320[k] = f_11 * osg_793[k]
                    + f_4 * qsf0_629[k]
                    - f_5 * qsf1_629[k]
                    + f_3 * pc_y[k] * qsg_943[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, pc_x, pc_y, pc_z, osg_779, osg_794, osg_945, \
                         qsf0_629, qsf0_630, qsf1_629, qsf1_630, qsg_944, \
                         qsg_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = f_11 * osg_794[k]
                    + f_3 * pc_y[k] * qsg_944[k];

        t_1322[k] = f_19 * osg_779[k]
                    + f_1 * qsf0_629[k]
                    - f_2 * qsf1_629[k]
                    + f_3 * pc_z[k] * qsg_944[k];

        t_1323[k] = f_10 * osg_945[k]
                    + f_1 * qsf0_630[k]
                    - f_2 * qsf1_630[k]
                    + f_3 * pc_x[k] * qsg_945[k];
    }

#pragma omp simd aligned(t_1324, t_1325, t_1326, t_1327, pc_x, pc_y, pc_z, osg_780, osg_795, \
                         osg_797, osg_948, qsf0_633, qsf1_633, qsg_945, qsg_947, \
                         qsg_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1324[k] = f_10 * osg_795[k]
                    + f_3 * pc_y[k] * qsg_945[k];

        t_1325[k] = f_17 * osg_780[k]
                    + f_3 * pc_z[k] * qsg_945[k];

        t_1326[k] = f_10 * osg_948[k]
                    + f_6 * qsf0_633[k]
                    - f_7 * qsf1_633[k]
                    + f_3 * pc_x[k] * qsg_948[k];

        t_1327[k] = f_10 * osg_797[k]
                    + f_3 * pc_y[k] * qsg_947[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pc_x, pc_z, osg_783, osg_950, osg_951, \
                         qsf0_635, qsf0_636, qsf1_635, qsf1_636, qsg_948, qsg_950, \
                         qsg_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_10 * osg_950[k]
                    + f_6 * qsf0_635[k]
                    - f_7 * qsf1_635[k]
                    + f_3 * pc_x[k] * qsg_950[k];

        t_1329[k] = f_10 * osg_951[k]
                    + f_4 * qsf0_636[k]
                    - f_5 * qsf1_636[k]
                    + f_3 * pc_x[k] * qsg_951[k];

        t_1330[k] = f_17 * osg_783[k]
                    + f_3 * pc_z[k] * qsg_948[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, pc_x, pc_y, osg_800, osg_954, \
                         osg_955, osg_956, qsf0_639, qsf1_639, qsg_950, qsg_954, qsg_955, \
                         qsg_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_10 * osg_800[k]
                    + f_3 * pc_y[k] * qsg_950[k];

        t_1332[k] = f_10 * osg_954[k]
                    + f_4 * qsf0_639[k]
                    - f_5 * qsf1_639[k]
                    + f_3 * pc_x[k] * qsg_954[k];

        t_1333[k] = f_10 * osg_955[k]
                    + f_3 * pc_x[k] * qsg_955[k];

        t_1334[k] = f_10 * osg_956[k]
                    + f_3 * pc_x[k] * qsg_956[k];
    }

#pragma omp simd aligned(t_1335, t_1336, t_1337, t_1338, pc_x, pc_y, osg_805, osg_957, \
                         osg_958, osg_959, qsf0_636, qsf1_636, qsg_955, qsg_957, qsg_958, \
                         qsg_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1335[k] = f_10 * osg_957[k]
                    + f_3 * pc_x[k] * qsg_957[k];

        t_1336[k] = f_10 * osg_958[k]
                    + f_3 * pc_x[k] * qsg_958[k];

        t_1337[k] = f_10 * osg_959[k]
                    + f_3 * pc_x[k] * qsg_959[k];

        t_1338[k] = f_10 * osg_805[k]
                    + f_1 * qsf0_636[k]
                    - f_2 * qsf1_636[k]
                    + f_3 * pc_y[k] * qsg_955[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, pc_y, pc_z, osg_790, osg_807, osg_808, \
                         qsf0_638, qsf0_639, qsf1_638, qsf1_639, qsg_955, qsg_957, \
                         qsg_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_17 * osg_790[k]
                    + f_3 * pc_z[k] * qsg_955[k];

        t_1340[k] = f_10 * osg_807[k]
                    + f_6 * qsf0_638[k]
                    - f_7 * qsf1_638[k]
                    + f_3 * pc_y[k] * qsg_957[k];

        t_1341[k] = f_10 * osg_808[k]
                    + f_4 * qsf0_639[k]
                    - f_5 * qsf1_639[k]
                    + f_3 * pc_y[k] * qsg_958[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, t_1345, pa_y, pc_y, pc_z, osh0_1134, osg_794, \
                         osg_809, osg_810, osh1_1134, qsf0_639, qsf1_639, qsg_959, \
                         qsg_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_10 * osg_809[k]
                    + f_3 * pc_y[k] * qsg_959[k];

        t_1343[k] = f_17 * osg_794[k]
                    + f_1 * qsf0_639[k]
                    - f_2 * qsf1_639[k]
                    + f_3 * pc_z[k] * qsg_959[k];

        t_1344[k] = pa_y[k] * osh0_1134[k]
                    - f_8 * pc_y[k] * osh1_1134[k];

        t_1345[k] = f_9 * osg_810[k]
                    + f_3 * pc_y[k] * qsg_960[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, t_1349, pa_y, pc_y, pc_z, osh0_1137, \
                         osh0_1139, osg_795, osg_811, osg_812, osh1_1137, osh1_1139, qsg_960, \
                         qsg_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_16 * osg_795[k]
                    + f_3 * pc_z[k] * qsg_960[k];

        t_1347[k] = pa_y[k] * osh0_1137[k]
                    + f_10 * osg_811[k]
                    - f_8 * pc_y[k] * osh1_1137[k];

        t_1348[k] = f_9 * osg_812[k]
                    + f_3 * pc_y[k] * qsg_962[k];

        t_1349[k] = pa_y[k] * osh0_1139[k]
                    - f_8 * pc_y[k] * osh1_1139[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, t_1353, pa_y, pc_y, pc_z, osh0_1140, \
                         osh0_1143, osg_798, osg_813, osg_815, osh1_1140, osh1_1143, qsg_963, \
                         qsg_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = pa_y[k] * osh0_1140[k]
                    + f_11 * osg_813[k]
                    - f_8 * pc_y[k] * osh1_1140[k];

        t_1351[k] = f_16 * osg_798[k]
                    + f_3 * pc_z[k] * qsg_963[k];

        t_1352[k] = f_9 * osg_815[k]
                    + f_3 * pc_y[k] * qsg_965[k];

        t_1353[k] = pa_y[k] * osh0_1143[k]
                    - f_8 * pc_y[k] * osh1_1143[k];
    }

#pragma omp simd aligned(t_1354, t_1355, t_1356, t_1357, t_1358, pc_x, osg_970, osg_971, \
                         osg_972, osg_973, osg_974, qsg_970, qsg_971, qsg_972, qsg_973, \
                         qsg_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1354[k] = f_10 * osg_970[k]
                    + f_3 * pc_x[k] * qsg_970[k];

        t_1355[k] = f_10 * osg_971[k]
                    + f_3 * pc_x[k] * qsg_971[k];

        t_1356[k] = f_10 * osg_972[k]
                    + f_3 * pc_x[k] * qsg_972[k];

        t_1357[k] = f_10 * osg_973[k]
                    + f_3 * pc_x[k] * qsg_973[k];

        t_1358[k] = f_10 * osg_974[k]
                    + f_3 * pc_x[k] * qsg_974[k];
    }

#pragma omp simd aligned(t_1359, t_1360, t_1361, pc_y, pc_z, osg_805, osg_820, osg_822, \
                         qsf0_646, qsf0_648, qsf1_646, qsf1_648, qsg_970, \
                         qsg_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1359[k] = f_9 * osg_820[k]
                    + f_1 * qsf0_646[k]
                    - f_2 * qsf1_646[k]
                    + f_3 * pc_y[k] * qsg_970[k];

        t_1360[k] = f_16 * osg_805[k]
                    + f_3 * pc_z[k] * qsg_970[k];

        t_1361[k] = f_9 * osg_822[k]
                    + f_6 * qsf0_648[k]
                    - f_7 * qsf1_648[k]
                    + f_3 * pc_y[k] * qsg_972[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, pa_y, pc_y, osh0_1154, osg_823, osg_824, \
                         osh1_1154, qsf0_649, qsf1_649, qsg_973, \
                         qsg_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_9 * osg_823[k]
                    + f_4 * qsf0_649[k]
                    - f_5 * qsf1_649[k]
                    + f_3 * pc_y[k] * qsg_973[k];

        t_1363[k] = f_9 * osg_824[k]
                    + f_3 * pc_y[k] * qsg_974[k];

        t_1364[k] = pa_y[k] * osh0_1154[k]
                    - f_8 * pc_y[k] * osh1_1154[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, t_1368, t_1369, pc_x, pc_y, pc_z, osg_810, \
                         osg_975, qsf0_650, qsf1_650, qsg_975, qsg_976, \
                         qsg_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = f_10 * osg_975[k]
                    + f_1 * qsf0_650[k]
                    - f_2 * qsf1_650[k]
                    + f_3 * pc_x[k] * qsg_975[k];

        t_1366[k] = f_3 * pc_y[k] * qsg_975[k];

        t_1367[k] = f_15 * osg_810[k]
                    + f_3 * pc_z[k] * qsg_975[k];

        t_1368[k] = f_4 * qsf0_650[k]
                    - f_5 * qsf1_650[k]
                    + f_3 * pc_y[k] * qsg_976[k];

        t_1369[k] = f_3 * pc_y[k] * qsg_977[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pc_x, pc_y, osg_980, qsf0_651, \
                         qsf0_652, qsf0_655, qsf1_651, qsf1_652, qsf1_655, qsg_978, qsg_979, \
                         qsg_980 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_10 * osg_980[k]
                    + f_6 * qsf0_655[k]
                    - f_7 * qsf1_655[k]
                    + f_3 * pc_x[k] * qsg_980[k];

        t_1371[k] = f_6 * qsf0_651[k]
                    - f_7 * qsf1_651[k]
                    + f_3 * pc_y[k] * qsg_978[k];

        t_1372[k] = f_4 * qsf0_652[k]
                    - f_5 * qsf1_652[k]
                    + f_3 * pc_y[k] * qsg_979[k];

        t_1373[k] = f_3 * pc_y[k] * qsg_980[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, t_1377, pc_x, osg_984, osg_985, osg_986, \
                         osg_987, qsf0_659, qsf1_659, qsg_984, qsg_985, qsg_986, \
                         qsg_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_10 * osg_984[k]
                    + f_4 * qsf0_659[k]
                    - f_5 * qsf1_659[k]
                    + f_3 * pc_x[k] * qsg_984[k];

        t_1375[k] = f_10 * osg_985[k]
                    + f_3 * pc_x[k] * qsg_985[k];

        t_1376[k] = f_10 * osg_986[k]
                    + f_3 * pc_x[k] * qsg_986[k];

        t_1377[k] = f_10 * osg_987[k]
                    + f_3 * pc_x[k] * qsg_987[k];
    }

#pragma omp simd aligned(t_1378, t_1379, t_1380, t_1381, pc_x, pc_y, osg_989, qsf0_656, \
                         qsf0_657, qsf1_656, qsf1_657, qsg_984, qsg_985, qsg_986, \
                         qsg_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1378[k] = f_3 * pc_y[k] * qsg_984[k];

        t_1379[k] = f_10 * osg_989[k]
                    + f_3 * pc_x[k] * qsg_989[k];

        t_1380[k] = f_1 * qsf0_656[k]
                    - f_2 * qsf1_656[k]
                    + f_3 * pc_y[k] * qsg_985[k];

        t_1381[k] = f_13 * qsf0_657[k]
                    - f_14 * qsf1_657[k]
                    + f_3 * pc_y[k] * qsg_986[k];
    }

#pragma omp simd aligned(t_1382, t_1383, t_1384, t_1385, pc_y, pc_z, osg_824, qsf0_658, \
                         qsf0_659, qsf1_658, qsf1_659, qsg_987, qsg_988, \
                         qsg_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1382[k] = f_6 * qsf0_658[k]
                    - f_7 * qsf1_658[k]
                    + f_3 * pc_y[k] * qsg_987[k];

        t_1383[k] = f_4 * qsf0_659[k]
                    - f_5 * qsf1_659[k]
                    + f_3 * pc_y[k] * qsg_988[k];

        t_1384[k] = f_3 * pc_y[k] * qsg_989[k];

        t_1385[k] = f_15 * osg_824[k]
                    + f_1 * qsf0_659[k]
                    - f_2 * qsf1_659[k]
                    + f_3 * pc_z[k] * qsg_989[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pa_x, pc_x, pc_y, pc_z, osh0_1386, \
                         osh0_1389, osg_825, osg_990, osg_993, osh1_1386, osh1_1389, \
                         qsg_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = pa_x[k] * osh0_1386[k]
                    + f_20 * osg_990[k]
                    - f_8 * pc_x[k] * osh1_1386[k];

        t_1387[k] = f_12 * osg_825[k]
                    + f_3 * pc_y[k] * qsg_990[k];

        t_1388[k] = f_3 * pc_z[k] * qsg_990[k];

        t_1389[k] = pa_x[k] * osh0_1389[k]
                    + f_11 * osg_993[k]
                    - f_8 * pc_x[k] * osh1_1389[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, pa_x, pc_x, pc_z, osh0_1392, osg_996, \
                         osh1_1392, qsf0_660, qsf1_660, qsg_991, qsg_992, \
                         qsg_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_3 * pc_z[k] * qsg_991[k];

        t_1391[k] = f_4 * qsf0_660[k]
                    - f_5 * qsf1_660[k]
                    + f_3 * pc_z[k] * qsg_992[k];

        t_1392[k] = pa_x[k] * osh0_1392[k]
                    + f_10 * osg_996[k]
                    - f_8 * pc_x[k] * osh1_1392[k];

        t_1393[k] = f_3 * pc_z[k] * qsg_993[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, t_1397, pc_x, pc_y, pc_z, osg_830, osg_1000, \
                         qsf0_662, qsf1_662, qsg_995, qsg_996, \
                         qsg_1000 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_12 * osg_830[k]
                    + f_3 * pc_y[k] * qsg_995[k];

        t_1395[k] = f_6 * qsf0_662[k]
                    - f_7 * qsf1_662[k]
                    + f_3 * pc_z[k] * qsg_995[k];

        t_1396[k] = f_9 * osg_1000[k]
                    + f_3 * pc_x[k] * qsg_1000[k];

        t_1397[k] = f_3 * pc_z[k] * qsg_996[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsg, const size_t ncols,
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
    const auto f_12 = 5.5 / q;
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_1155 = buffer.data(osh0 + 1155);
    const auto *osh0_1158 = buffer.data(osh0 + 1158);
    const auto *osh0_1161 = buffer.data(osh0 + 1161);
    const auto *osh0_1401 = buffer.data(osh0 + 1401);
    const auto *osh0_1403 = buffer.data(osh0 + 1403);
    const auto *osh0_1404 = buffer.data(osh0 + 1404);
    const auto *osh0_1406 = buffer.data(osh0 + 1406);
    const auto *osh0_1412 = buffer.data(osh0 + 1412);
    const auto *osh0_1416 = buffer.data(osh0 + 1416);
    const auto *osh0_1422 = buffer.data(osh0 + 1422);
    const auto *osh0_1424 = buffer.data(osh0 + 1424);
    const auto *osh0_1425 = buffer.data(osh0 + 1425);
    const auto *osh0_1427 = buffer.data(osh0 + 1427);
    const auto *osh0_1428 = buffer.data(osh0 + 1428);
    const auto *osh0_1431 = buffer.data(osh0 + 1431);
    const auto *osh0_1433 = buffer.data(osh0 + 1433);
    const auto *osh0_1434 = buffer.data(osh0 + 1434);
    const auto *osh0_1437 = buffer.data(osh0 + 1437);
    const auto *osh0_1443 = buffer.data(osh0 + 1443);
    const auto *osh0_1445 = buffer.data(osh0 + 1445);
    const auto *osh0_1446 = buffer.data(osh0 + 1446);
    const auto *osh0_1448 = buffer.data(osh0 + 1448);
    const auto *osh0_1449 = buffer.data(osh0 + 1449);
    const auto *osh0_1452 = buffer.data(osh0 + 1452);
    const auto *osh0_1454 = buffer.data(osh0 + 1454);
    const auto *osh0_1455 = buffer.data(osh0 + 1455);
    const auto *osh0_1458 = buffer.data(osh0 + 1458);
    const auto *osh0_1464 = buffer.data(osh0 + 1464);
    const auto *osh0_1466 = buffer.data(osh0 + 1466);
    const auto *osh0_1467 = buffer.data(osh0 + 1467);
    const auto *osh0_1469 = buffer.data(osh0 + 1469);
    const auto *osh0_1470 = buffer.data(osh0 + 1470);
    const auto *osh0_1473 = buffer.data(osh0 + 1473);
    const auto *osh0_1475 = buffer.data(osh0 + 1475);
    const auto *osh0_1476 = buffer.data(osh0 + 1476);
    const auto *osh0_1479 = buffer.data(osh0 + 1479);
    const auto *osh0_1485 = buffer.data(osh0 + 1485);
    const auto *osh0_1487 = buffer.data(osh0 + 1487);
    const auto *osh0_1488 = buffer.data(osh0 + 1488);
    const auto *osh0_1490 = buffer.data(osh0 + 1490);
    const auto *osh0_1491 = buffer.data(osh0 + 1491);
    const auto *osh0_1494 = buffer.data(osh0 + 1494);
    const auto *osh0_1496 = buffer.data(osh0 + 1496);
    const auto *osh0_1497 = buffer.data(osh0 + 1497);
    const auto *osh0_1500 = buffer.data(osh0 + 1500);
    const auto *osh0_1506 = buffer.data(osh0 + 1506);
    const auto *osh0_1508 = buffer.data(osh0 + 1508);
    const auto *osh0_1509 = buffer.data(osh0 + 1509);
    const auto *osh0_1511 = buffer.data(osh0 + 1511);
    const auto *osh0_1512 = buffer.data(osh0 + 1512);
    const auto *osh0_1515 = buffer.data(osh0 + 1515);
    const auto *osh0_1517 = buffer.data(osh0 + 1517);

    const auto *osg_825 = buffer.data(osg + 825);
    const auto *osg_828 = buffer.data(osg + 828);
    const auto *osg_835 = buffer.data(osg + 835);
    const auto *osg_839 = buffer.data(osg + 839);
    const auto *osg_840 = buffer.data(osg + 840);
    const auto *osg_842 = buffer.data(osg + 842);
    const auto *osg_843 = buffer.data(osg + 843);
    const auto *osg_845 = buffer.data(osg + 845);
    const auto *osg_850 = buffer.data(osg + 850);
    const auto *osg_854 = buffer.data(osg + 854);
    const auto *osg_855 = buffer.data(osg + 855);
    const auto *osg_857 = buffer.data(osg + 857);
    const auto *osg_858 = buffer.data(osg + 858);
    const auto *osg_860 = buffer.data(osg + 860);
    const auto *osg_865 = buffer.data(osg + 865);
    const auto *osg_869 = buffer.data(osg + 869);
    const auto *osg_870 = buffer.data(osg + 870);
    const auto *osg_872 = buffer.data(osg + 872);
    const auto *osg_873 = buffer.data(osg + 873);
    const auto *osg_875 = buffer.data(osg + 875);
    const auto *osg_880 = buffer.data(osg + 880);
    const auto *osg_884 = buffer.data(osg + 884);
    const auto *osg_885 = buffer.data(osg + 885);
    const auto *osg_887 = buffer.data(osg + 887);
    const auto *osg_888 = buffer.data(osg + 888);
    const auto *osg_890 = buffer.data(osg + 890);
    const auto *osg_895 = buffer.data(osg + 895);
    const auto *osg_899 = buffer.data(osg + 899);
    const auto *osg_900 = buffer.data(osg + 900);
    const auto *osg_902 = buffer.data(osg + 902);
    const auto *osg_905 = buffer.data(osg + 905);
    const auto *osg_914 = buffer.data(osg + 914);
    const auto *osg_915 = buffer.data(osg + 915);
    const auto *osg_917 = buffer.data(osg + 917);
    const auto *osg_1002 = buffer.data(osg + 1002);
    const auto *osg_1003 = buffer.data(osg + 1003);
    const auto *osg_1004 = buffer.data(osg + 1004);
    const auto *osg_1010 = buffer.data(osg + 1010);
    const auto *osg_1014 = buffer.data(osg + 1014);
    const auto *osg_1015 = buffer.data(osg + 1015);
    const auto *osg_1016 = buffer.data(osg + 1016);
    const auto *osg_1017 = buffer.data(osg + 1017);
    const auto *osg_1018 = buffer.data(osg + 1018);
    const auto *osg_1019 = buffer.data(osg + 1019);
    const auto *osg_1020 = buffer.data(osg + 1020);
    const auto *osg_1023 = buffer.data(osg + 1023);
    const auto *osg_1025 = buffer.data(osg + 1025);
    const auto *osg_1026 = buffer.data(osg + 1026);
    const auto *osg_1029 = buffer.data(osg + 1029);
    const auto *osg_1030 = buffer.data(osg + 1030);
    const auto *osg_1031 = buffer.data(osg + 1031);
    const auto *osg_1032 = buffer.data(osg + 1032);
    const auto *osg_1033 = buffer.data(osg + 1033);
    const auto *osg_1034 = buffer.data(osg + 1034);
    const auto *osg_1035 = buffer.data(osg + 1035);
    const auto *osg_1038 = buffer.data(osg + 1038);
    const auto *osg_1040 = buffer.data(osg + 1040);
    const auto *osg_1041 = buffer.data(osg + 1041);
    const auto *osg_1044 = buffer.data(osg + 1044);
    const auto *osg_1045 = buffer.data(osg + 1045);
    const auto *osg_1046 = buffer.data(osg + 1046);
    const auto *osg_1047 = buffer.data(osg + 1047);
    const auto *osg_1048 = buffer.data(osg + 1048);
    const auto *osg_1049 = buffer.data(osg + 1049);
    const auto *osg_1050 = buffer.data(osg + 1050);
    const auto *osg_1053 = buffer.data(osg + 1053);
    const auto *osg_1055 = buffer.data(osg + 1055);
    const auto *osg_1056 = buffer.data(osg + 1056);
    const auto *osg_1059 = buffer.data(osg + 1059);
    const auto *osg_1060 = buffer.data(osg + 1060);
    const auto *osg_1061 = buffer.data(osg + 1061);
    const auto *osg_1062 = buffer.data(osg + 1062);
    const auto *osg_1063 = buffer.data(osg + 1063);
    const auto *osg_1064 = buffer.data(osg + 1064);
    const auto *osg_1065 = buffer.data(osg + 1065);
    const auto *osg_1068 = buffer.data(osg + 1068);
    const auto *osg_1070 = buffer.data(osg + 1070);
    const auto *osg_1071 = buffer.data(osg + 1071);
    const auto *osg_1074 = buffer.data(osg + 1074);
    const auto *osg_1075 = buffer.data(osg + 1075);
    const auto *osg_1076 = buffer.data(osg + 1076);
    const auto *osg_1077 = buffer.data(osg + 1077);
    const auto *osg_1078 = buffer.data(osg + 1078);
    const auto *osg_1079 = buffer.data(osg + 1079);
    const auto *osg_1080 = buffer.data(osg + 1080);
    const auto *osg_1083 = buffer.data(osg + 1083);
    const auto *osg_1085 = buffer.data(osg + 1085);

    const auto *osh1_1155 = buffer.data(osh1 + 1155);
    const auto *osh1_1158 = buffer.data(osh1 + 1158);
    const auto *osh1_1161 = buffer.data(osh1 + 1161);
    const auto *osh1_1401 = buffer.data(osh1 + 1401);
    const auto *osh1_1403 = buffer.data(osh1 + 1403);
    const auto *osh1_1404 = buffer.data(osh1 + 1404);
    const auto *osh1_1406 = buffer.data(osh1 + 1406);
    const auto *osh1_1412 = buffer.data(osh1 + 1412);
    const auto *osh1_1416 = buffer.data(osh1 + 1416);
    const auto *osh1_1422 = buffer.data(osh1 + 1422);
    const auto *osh1_1424 = buffer.data(osh1 + 1424);
    const auto *osh1_1425 = buffer.data(osh1 + 1425);
    const auto *osh1_1427 = buffer.data(osh1 + 1427);
    const auto *osh1_1428 = buffer.data(osh1 + 1428);
    const auto *osh1_1431 = buffer.data(osh1 + 1431);
    const auto *osh1_1433 = buffer.data(osh1 + 1433);
    const auto *osh1_1434 = buffer.data(osh1 + 1434);
    const auto *osh1_1437 = buffer.data(osh1 + 1437);
    const auto *osh1_1443 = buffer.data(osh1 + 1443);
    const auto *osh1_1445 = buffer.data(osh1 + 1445);
    const auto *osh1_1446 = buffer.data(osh1 + 1446);
    const auto *osh1_1448 = buffer.data(osh1 + 1448);
    const auto *osh1_1449 = buffer.data(osh1 + 1449);
    const auto *osh1_1452 = buffer.data(osh1 + 1452);
    const auto *osh1_1454 = buffer.data(osh1 + 1454);
    const auto *osh1_1455 = buffer.data(osh1 + 1455);
    const auto *osh1_1458 = buffer.data(osh1 + 1458);
    const auto *osh1_1464 = buffer.data(osh1 + 1464);
    const auto *osh1_1466 = buffer.data(osh1 + 1466);
    const auto *osh1_1467 = buffer.data(osh1 + 1467);
    const auto *osh1_1469 = buffer.data(osh1 + 1469);
    const auto *osh1_1470 = buffer.data(osh1 + 1470);
    const auto *osh1_1473 = buffer.data(osh1 + 1473);
    const auto *osh1_1475 = buffer.data(osh1 + 1475);
    const auto *osh1_1476 = buffer.data(osh1 + 1476);
    const auto *osh1_1479 = buffer.data(osh1 + 1479);
    const auto *osh1_1485 = buffer.data(osh1 + 1485);
    const auto *osh1_1487 = buffer.data(osh1 + 1487);
    const auto *osh1_1488 = buffer.data(osh1 + 1488);
    const auto *osh1_1490 = buffer.data(osh1 + 1490);
    const auto *osh1_1491 = buffer.data(osh1 + 1491);
    const auto *osh1_1494 = buffer.data(osh1 + 1494);
    const auto *osh1_1496 = buffer.data(osh1 + 1496);
    const auto *osh1_1497 = buffer.data(osh1 + 1497);
    const auto *osh1_1500 = buffer.data(osh1 + 1500);
    const auto *osh1_1506 = buffer.data(osh1 + 1506);
    const auto *osh1_1508 = buffer.data(osh1 + 1508);
    const auto *osh1_1509 = buffer.data(osh1 + 1509);
    const auto *osh1_1511 = buffer.data(osh1 + 1511);
    const auto *osh1_1512 = buffer.data(osh1 + 1512);
    const auto *osh1_1515 = buffer.data(osh1 + 1515);
    const auto *osh1_1517 = buffer.data(osh1 + 1517);

    const auto *qsg_1000 = buffer.data(qsg + 1000);
    const auto *qsg_1002 = buffer.data(qsg + 1002);
    const auto *qsg_1003 = buffer.data(qsg + 1003);
    const auto *qsg_1004 = buffer.data(qsg + 1004);
    const auto *qsg_1005 = buffer.data(qsg + 1005);
    const auto *qsg_1007 = buffer.data(qsg + 1007);
    const auto *qsg_1008 = buffer.data(qsg + 1008);
    const auto *qsg_1010 = buffer.data(qsg + 1010);
    const auto *qsg_1015 = buffer.data(qsg + 1015);
    const auto *qsg_1016 = buffer.data(qsg + 1016);
    const auto *qsg_1017 = buffer.data(qsg + 1017);
    const auto *qsg_1018 = buffer.data(qsg + 1018);
    const auto *qsg_1019 = buffer.data(qsg + 1019);
    const auto *qsg_1020 = buffer.data(qsg + 1020);
    const auto *qsg_1022 = buffer.data(qsg + 1022);
    const auto *qsg_1023 = buffer.data(qsg + 1023);
    const auto *qsg_1025 = buffer.data(qsg + 1025);
    const auto *qsg_1030 = buffer.data(qsg + 1030);
    const auto *qsg_1031 = buffer.data(qsg + 1031);
    const auto *qsg_1032 = buffer.data(qsg + 1032);
    const auto *qsg_1033 = buffer.data(qsg + 1033);
    const auto *qsg_1034 = buffer.data(qsg + 1034);
    const auto *qsg_1035 = buffer.data(qsg + 1035);
    const auto *qsg_1037 = buffer.data(qsg + 1037);
    const auto *qsg_1038 = buffer.data(qsg + 1038);
    const auto *qsg_1040 = buffer.data(qsg + 1040);
    const auto *qsg_1045 = buffer.data(qsg + 1045);
    const auto *qsg_1046 = buffer.data(qsg + 1046);
    const auto *qsg_1047 = buffer.data(qsg + 1047);
    const auto *qsg_1048 = buffer.data(qsg + 1048);
    const auto *qsg_1049 = buffer.data(qsg + 1049);
    const auto *qsg_1050 = buffer.data(qsg + 1050);
    const auto *qsg_1052 = buffer.data(qsg + 1052);
    const auto *qsg_1053 = buffer.data(qsg + 1053);
    const auto *qsg_1055 = buffer.data(qsg + 1055);
    const auto *qsg_1060 = buffer.data(qsg + 1060);
    const auto *qsg_1061 = buffer.data(qsg + 1061);
    const auto *qsg_1062 = buffer.data(qsg + 1062);
    const auto *qsg_1063 = buffer.data(qsg + 1063);
    const auto *qsg_1064 = buffer.data(qsg + 1064);
    const auto *qsg_1065 = buffer.data(qsg + 1065);
    const auto *qsg_1067 = buffer.data(qsg + 1067);
    const auto *qsg_1068 = buffer.data(qsg + 1068);
    const auto *qsg_1070 = buffer.data(qsg + 1070);
    const auto *qsg_1075 = buffer.data(qsg + 1075);
    const auto *qsg_1076 = buffer.data(qsg + 1076);
    const auto *qsg_1077 = buffer.data(qsg + 1077);
    const auto *qsg_1078 = buffer.data(qsg + 1078);
    const auto *qsg_1079 = buffer.data(qsg + 1079);
    const auto *qsg_1080 = buffer.data(qsg + 1080);
    const auto *qsg_1082 = buffer.data(qsg + 1082);

#pragma omp simd aligned(t_1398, t_1399, t_1400, t_1401, pa_x, pc_x, osh0_1401, osg_1002, \
                         osg_1003, osg_1004, osh1_1401, qsg_1002, qsg_1003, \
                         qsg_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_9 * osg_1002[k]
                    + f_3 * pc_x[k] * qsg_1002[k];

        t_1399[k] = f_9 * osg_1003[k]
                    + f_3 * pc_x[k] * qsg_1003[k];

        t_1400[k] = f_9 * osg_1004[k]
                    + f_3 * pc_x[k] * qsg_1004[k];

        t_1401[k] = pa_x[k] * osh0_1401[k]
                    - f_8 * pc_x[k] * osh1_1401[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, t_1405, pa_x, pc_x, pc_y, pc_z, osh0_1403, \
                         osh0_1404, osg_839, osh1_1403, osh1_1404, qsg_1000, \
                         qsg_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_3 * pc_z[k] * qsg_1000[k];

        t_1403[k] = pa_x[k] * osh0_1403[k]
                    - f_8 * pc_x[k] * osh1_1403[k];

        t_1404[k] = pa_x[k] * osh0_1404[k]
                    - f_8 * pc_x[k] * osh1_1404[k];

        t_1405[k] = f_12 * osg_839[k]
                    + f_3 * pc_y[k] * qsg_1004[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, t_1409, pa_x, pa_z, pc_x, pc_y, pc_z, \
                         osh0_1155, osh0_1406, osg_825, osg_840, osh1_1155, osh1_1406, \
                         qsg_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = pa_x[k] * osh0_1406[k]
                    - f_8 * pc_x[k] * osh1_1406[k];

        t_1407[k] = pa_z[k] * osh0_1155[k]
                    - f_8 * pc_z[k] * osh1_1155[k];

        t_1408[k] = f_15 * osg_840[k]
                    + f_3 * pc_y[k] * qsg_1005[k];

        t_1409[k] = f_9 * osg_825[k]
                    + f_3 * pc_z[k] * qsg_1005[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pa_x, pa_z, pc_x, pc_y, pc_z, osh0_1158, \
                         osh0_1412, osg_842, osg_1010, osh1_1158, osh1_1412, \
                         qsg_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = pa_z[k] * osh0_1158[k]
                    - f_8 * pc_z[k] * osh1_1158[k];

        t_1411[k] = f_15 * osg_842[k]
                    + f_3 * pc_y[k] * qsg_1007[k];

        t_1412[k] = pa_x[k] * osh0_1412[k]
                    + f_11 * osg_1010[k]
                    - f_8 * pc_x[k] * osh1_1412[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, pa_z, pc_y, pc_z, osh0_1161, osg_828, \
                         osg_845, osh1_1161, qsg_1008, qsg_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = pa_z[k] * osh0_1161[k]
                    - f_8 * pc_z[k] * osh1_1161[k];

        t_1414[k] = f_9 * osg_828[k]
                    + f_3 * pc_z[k] * qsg_1008[k];

        t_1415[k] = f_15 * osg_845[k]
                    + f_3 * pc_y[k] * qsg_1010[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, t_1419, pa_x, pc_x, osh0_1416, osg_1014, \
                         osg_1015, osg_1016, osg_1017, osh1_1416, qsg_1015, qsg_1016, \
                         qsg_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = pa_x[k] * osh0_1416[k]
                    + f_10 * osg_1014[k]
                    - f_8 * pc_x[k] * osh1_1416[k];

        t_1417[k] = f_9 * osg_1015[k]
                    + f_3 * pc_x[k] * qsg_1015[k];

        t_1418[k] = f_9 * osg_1016[k]
                    + f_3 * pc_x[k] * qsg_1016[k];

        t_1419[k] = f_9 * osg_1017[k]
                    + f_3 * pc_x[k] * qsg_1017[k];
    }

#pragma omp simd aligned(t_1420, t_1421, t_1422, t_1423, pa_x, pc_x, pc_z, osh0_1422, osg_835, \
                         osg_1018, osg_1019, osh1_1422, qsg_1015, qsg_1018, \
                         qsg_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1420[k] = f_9 * osg_1018[k]
                    + f_3 * pc_x[k] * qsg_1018[k];

        t_1421[k] = f_9 * osg_1019[k]
                    + f_3 * pc_x[k] * qsg_1019[k];

        t_1422[k] = pa_x[k] * osh0_1422[k]
                    - f_8 * pc_x[k] * osh1_1422[k];

        t_1423[k] = f_9 * osg_835[k]
                    + f_3 * pc_z[k] * qsg_1015[k];
    }

#pragma omp simd aligned(t_1424, t_1425, t_1426, t_1427, pa_x, pc_x, pc_y, osh0_1424, \
                         osh0_1425, osh0_1427, osg_854, osh1_1424, osh1_1425, osh1_1427, \
                         qsg_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1424[k] = pa_x[k] * osh0_1424[k]
                    - f_8 * pc_x[k] * osh1_1424[k];

        t_1425[k] = pa_x[k] * osh0_1425[k]
                    - f_8 * pc_x[k] * osh1_1425[k];

        t_1426[k] = f_15 * osg_854[k]
                    + f_3 * pc_y[k] * qsg_1019[k];

        t_1427[k] = pa_x[k] * osh0_1427[k]
                    - f_8 * pc_x[k] * osh1_1427[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, pa_x, pc_x, pc_y, pc_z, osh0_1428, osg_840, \
                         osg_855, osg_1020, osh1_1428, qsg_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = pa_x[k] * osh0_1428[k]
                    + f_20 * osg_1020[k]
                    - f_8 * pc_x[k] * osh1_1428[k];

        t_1429[k] = f_16 * osg_855[k]
                    + f_3 * pc_y[k] * qsg_1020[k];

        t_1430[k] = f_10 * osg_840[k]
                    + f_3 * pc_z[k] * qsg_1020[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, pa_x, pc_x, pc_y, osh0_1431, osh0_1433, \
                         osg_857, osg_1023, osg_1025, osh1_1431, osh1_1433, \
                         qsg_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = pa_x[k] * osh0_1431[k]
                    + f_11 * osg_1023[k]
                    - f_8 * pc_x[k] * osh1_1431[k];

        t_1432[k] = f_16 * osg_857[k]
                    + f_3 * pc_y[k] * qsg_1022[k];

        t_1433[k] = pa_x[k] * osh0_1433[k]
                    + f_11 * osg_1025[k]
                    - f_8 * pc_x[k] * osh1_1433[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, pa_x, pc_x, pc_y, pc_z, osh0_1434, osg_843, \
                         osg_860, osg_1026, osh1_1434, qsg_1023, \
                         qsg_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = pa_x[k] * osh0_1434[k]
                    + f_10 * osg_1026[k]
                    - f_8 * pc_x[k] * osh1_1434[k];

        t_1435[k] = f_10 * osg_843[k]
                    + f_3 * pc_z[k] * qsg_1023[k];

        t_1436[k] = f_16 * osg_860[k]
                    + f_3 * pc_y[k] * qsg_1025[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, t_1440, pa_x, pc_x, osh0_1437, osg_1029, \
                         osg_1030, osg_1031, osg_1032, osh1_1437, qsg_1030, qsg_1031, \
                         qsg_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = pa_x[k] * osh0_1437[k]
                    + f_10 * osg_1029[k]
                    - f_8 * pc_x[k] * osh1_1437[k];

        t_1438[k] = f_9 * osg_1030[k]
                    + f_3 * pc_x[k] * qsg_1030[k];

        t_1439[k] = f_9 * osg_1031[k]
                    + f_3 * pc_x[k] * qsg_1031[k];

        t_1440[k] = f_9 * osg_1032[k]
                    + f_3 * pc_x[k] * qsg_1032[k];
    }

#pragma omp simd aligned(t_1441, t_1442, t_1443, t_1444, pa_x, pc_x, pc_z, osh0_1443, osg_850, \
                         osg_1033, osg_1034, osh1_1443, qsg_1030, qsg_1033, \
                         qsg_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1441[k] = f_9 * osg_1033[k]
                    + f_3 * pc_x[k] * qsg_1033[k];

        t_1442[k] = f_9 * osg_1034[k]
                    + f_3 * pc_x[k] * qsg_1034[k];

        t_1443[k] = pa_x[k] * osh0_1443[k]
                    - f_8 * pc_x[k] * osh1_1443[k];

        t_1444[k] = f_10 * osg_850[k]
                    + f_3 * pc_z[k] * qsg_1030[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, t_1448, pa_x, pc_x, pc_y, osh0_1445, \
                         osh0_1446, osh0_1448, osg_869, osh1_1445, osh1_1446, osh1_1448, \
                         qsg_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = pa_x[k] * osh0_1445[k]
                    - f_8 * pc_x[k] * osh1_1445[k];

        t_1446[k] = pa_x[k] * osh0_1446[k]
                    - f_8 * pc_x[k] * osh1_1446[k];

        t_1447[k] = f_16 * osg_869[k]
                    + f_3 * pc_y[k] * qsg_1034[k];

        t_1448[k] = pa_x[k] * osh0_1448[k]
                    - f_8 * pc_x[k] * osh1_1448[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, pa_x, pc_x, pc_y, pc_z, osh0_1449, osg_855, \
                         osg_870, osg_1035, osh1_1449, qsg_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = pa_x[k] * osh0_1449[k]
                    + f_20 * osg_1035[k]
                    - f_8 * pc_x[k] * osh1_1449[k];

        t_1450[k] = f_17 * osg_870[k]
                    + f_3 * pc_y[k] * qsg_1035[k];

        t_1451[k] = f_11 * osg_855[k]
                    + f_3 * pc_z[k] * qsg_1035[k];
    }

#pragma omp simd aligned(t_1452, t_1453, t_1454, pa_x, pc_x, pc_y, osh0_1452, osh0_1454, \
                         osg_872, osg_1038, osg_1040, osh1_1452, osh1_1454, \
                         qsg_1037 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1452[k] = pa_x[k] * osh0_1452[k]
                    + f_11 * osg_1038[k]
                    - f_8 * pc_x[k] * osh1_1452[k];

        t_1453[k] = f_17 * osg_872[k]
                    + f_3 * pc_y[k] * qsg_1037[k];

        t_1454[k] = pa_x[k] * osh0_1454[k]
                    + f_11 * osg_1040[k]
                    - f_8 * pc_x[k] * osh1_1454[k];
    }

#pragma omp simd aligned(t_1455, t_1456, t_1457, pa_x, pc_x, pc_y, pc_z, osh0_1455, osg_858, \
                         osg_875, osg_1041, osh1_1455, qsg_1038, \
                         qsg_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1455[k] = pa_x[k] * osh0_1455[k]
                    + f_10 * osg_1041[k]
                    - f_8 * pc_x[k] * osh1_1455[k];

        t_1456[k] = f_11 * osg_858[k]
                    + f_3 * pc_z[k] * qsg_1038[k];

        t_1457[k] = f_17 * osg_875[k]
                    + f_3 * pc_y[k] * qsg_1040[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, t_1461, pa_x, pc_x, osh0_1458, osg_1044, \
                         osg_1045, osg_1046, osg_1047, osh1_1458, qsg_1045, qsg_1046, \
                         qsg_1047 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = pa_x[k] * osh0_1458[k]
                    + f_10 * osg_1044[k]
                    - f_8 * pc_x[k] * osh1_1458[k];

        t_1459[k] = f_9 * osg_1045[k]
                    + f_3 * pc_x[k] * qsg_1045[k];

        t_1460[k] = f_9 * osg_1046[k]
                    + f_3 * pc_x[k] * qsg_1046[k];

        t_1461[k] = f_9 * osg_1047[k]
                    + f_3 * pc_x[k] * qsg_1047[k];
    }

#pragma omp simd aligned(t_1462, t_1463, t_1464, t_1465, pa_x, pc_x, pc_z, osh0_1464, osg_865, \
                         osg_1048, osg_1049, osh1_1464, qsg_1045, qsg_1048, \
                         qsg_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1462[k] = f_9 * osg_1048[k]
                    + f_3 * pc_x[k] * qsg_1048[k];

        t_1463[k] = f_9 * osg_1049[k]
                    + f_3 * pc_x[k] * qsg_1049[k];

        t_1464[k] = pa_x[k] * osh0_1464[k]
                    - f_8 * pc_x[k] * osh1_1464[k];

        t_1465[k] = f_11 * osg_865[k]
                    + f_3 * pc_z[k] * qsg_1045[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pa_x, pc_x, pc_y, osh0_1466, \
                         osh0_1467, osh0_1469, osg_884, osh1_1466, osh1_1467, osh1_1469, \
                         qsg_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = pa_x[k] * osh0_1466[k]
                    - f_8 * pc_x[k] * osh1_1466[k];

        t_1467[k] = pa_x[k] * osh0_1467[k]
                    - f_8 * pc_x[k] * osh1_1467[k];

        t_1468[k] = f_17 * osg_884[k]
                    + f_3 * pc_y[k] * qsg_1049[k];

        t_1469[k] = pa_x[k] * osh0_1469[k]
                    - f_8 * pc_x[k] * osh1_1469[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, pa_x, pc_x, pc_y, pc_z, osh0_1470, osg_870, \
                         osg_885, osg_1050, osh1_1470, qsg_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = pa_x[k] * osh0_1470[k]
                    + f_20 * osg_1050[k]
                    - f_8 * pc_x[k] * osh1_1470[k];

        t_1471[k] = f_19 * osg_885[k]
                    + f_3 * pc_y[k] * qsg_1050[k];

        t_1472[k] = f_18 * osg_870[k]
                    + f_3 * pc_z[k] * qsg_1050[k];
    }

#pragma omp simd aligned(t_1473, t_1474, t_1475, pa_x, pc_x, pc_y, osh0_1473, osh0_1475, \
                         osg_887, osg_1053, osg_1055, osh1_1473, osh1_1475, \
                         qsg_1052 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1473[k] = pa_x[k] * osh0_1473[k]
                    + f_11 * osg_1053[k]
                    - f_8 * pc_x[k] * osh1_1473[k];

        t_1474[k] = f_19 * osg_887[k]
                    + f_3 * pc_y[k] * qsg_1052[k];

        t_1475[k] = pa_x[k] * osh0_1475[k]
                    + f_11 * osg_1055[k]
                    - f_8 * pc_x[k] * osh1_1475[k];
    }

#pragma omp simd aligned(t_1476, t_1477, t_1478, pa_x, pc_x, pc_y, pc_z, osh0_1476, osg_873, \
                         osg_890, osg_1056, osh1_1476, qsg_1053, \
                         qsg_1055 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1476[k] = pa_x[k] * osh0_1476[k]
                    + f_10 * osg_1056[k]
                    - f_8 * pc_x[k] * osh1_1476[k];

        t_1477[k] = f_18 * osg_873[k]
                    + f_3 * pc_z[k] * qsg_1053[k];

        t_1478[k] = f_19 * osg_890[k]
                    + f_3 * pc_y[k] * qsg_1055[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, t_1482, pa_x, pc_x, osh0_1479, osg_1059, \
                         osg_1060, osg_1061, osg_1062, osh1_1479, qsg_1060, qsg_1061, \
                         qsg_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = pa_x[k] * osh0_1479[k]
                    + f_10 * osg_1059[k]
                    - f_8 * pc_x[k] * osh1_1479[k];

        t_1480[k] = f_9 * osg_1060[k]
                    + f_3 * pc_x[k] * qsg_1060[k];

        t_1481[k] = f_9 * osg_1061[k]
                    + f_3 * pc_x[k] * qsg_1061[k];

        t_1482[k] = f_9 * osg_1062[k]
                    + f_3 * pc_x[k] * qsg_1062[k];
    }

#pragma omp simd aligned(t_1483, t_1484, t_1485, t_1486, pa_x, pc_x, pc_z, osh0_1485, osg_880, \
                         osg_1063, osg_1064, osh1_1485, qsg_1060, qsg_1063, \
                         qsg_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1483[k] = f_9 * osg_1063[k]
                    + f_3 * pc_x[k] * qsg_1063[k];

        t_1484[k] = f_9 * osg_1064[k]
                    + f_3 * pc_x[k] * qsg_1064[k];

        t_1485[k] = pa_x[k] * osh0_1485[k]
                    - f_8 * pc_x[k] * osh1_1485[k];

        t_1486[k] = f_18 * osg_880[k]
                    + f_3 * pc_z[k] * qsg_1060[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, t_1490, pa_x, pc_x, pc_y, osh0_1487, \
                         osh0_1488, osh0_1490, osg_899, osh1_1487, osh1_1488, osh1_1490, \
                         qsg_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = pa_x[k] * osh0_1487[k]
                    - f_8 * pc_x[k] * osh1_1487[k];

        t_1488[k] = pa_x[k] * osh0_1488[k]
                    - f_8 * pc_x[k] * osh1_1488[k];

        t_1489[k] = f_19 * osg_899[k]
                    + f_3 * pc_y[k] * qsg_1064[k];

        t_1490[k] = pa_x[k] * osh0_1490[k]
                    - f_8 * pc_x[k] * osh1_1490[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pa_x, pc_x, pc_y, pc_z, osh0_1491, osg_885, \
                         osg_900, osg_1065, osh1_1491, qsg_1065 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = pa_x[k] * osh0_1491[k]
                    + f_20 * osg_1065[k]
                    - f_8 * pc_x[k] * osh1_1491[k];

        t_1492[k] = f_21 * osg_900[k]
                    + f_3 * pc_y[k] * qsg_1065[k];

        t_1493[k] = f_20 * osg_885[k]
                    + f_3 * pc_z[k] * qsg_1065[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, pa_x, pc_x, pc_y, osh0_1494, osh0_1496, \
                         osg_902, osg_1068, osg_1070, osh1_1494, osh1_1496, \
                         qsg_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = pa_x[k] * osh0_1494[k]
                    + f_11 * osg_1068[k]
                    - f_8 * pc_x[k] * osh1_1494[k];

        t_1495[k] = f_21 * osg_902[k]
                    + f_3 * pc_y[k] * qsg_1067[k];

        t_1496[k] = pa_x[k] * osh0_1496[k]
                    + f_11 * osg_1070[k]
                    - f_8 * pc_x[k] * osh1_1496[k];
    }

#pragma omp simd aligned(t_1497, t_1498, t_1499, pa_x, pc_x, pc_y, pc_z, osh0_1497, osg_888, \
                         osg_905, osg_1071, osh1_1497, qsg_1068, \
                         qsg_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1497[k] = pa_x[k] * osh0_1497[k]
                    + f_10 * osg_1071[k]
                    - f_8 * pc_x[k] * osh1_1497[k];

        t_1498[k] = f_20 * osg_888[k]
                    + f_3 * pc_z[k] * qsg_1068[k];

        t_1499[k] = f_21 * osg_905[k]
                    + f_3 * pc_y[k] * qsg_1070[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, pa_x, pc_x, osh0_1500, osg_1074, \
                         osg_1075, osg_1076, osg_1077, osh1_1500, qsg_1075, qsg_1076, \
                         qsg_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = pa_x[k] * osh0_1500[k]
                    + f_10 * osg_1074[k]
                    - f_8 * pc_x[k] * osh1_1500[k];

        t_1501[k] = f_9 * osg_1075[k]
                    + f_3 * pc_x[k] * qsg_1075[k];

        t_1502[k] = f_9 * osg_1076[k]
                    + f_3 * pc_x[k] * qsg_1076[k];

        t_1503[k] = f_9 * osg_1077[k]
                    + f_3 * pc_x[k] * qsg_1077[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, t_1507, pa_x, pc_x, pc_z, osh0_1506, osg_895, \
                         osg_1078, osg_1079, osh1_1506, qsg_1075, qsg_1078, \
                         qsg_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = f_9 * osg_1078[k]
                    + f_3 * pc_x[k] * qsg_1078[k];

        t_1505[k] = f_9 * osg_1079[k]
                    + f_3 * pc_x[k] * qsg_1079[k];

        t_1506[k] = pa_x[k] * osh0_1506[k]
                    - f_8 * pc_x[k] * osh1_1506[k];

        t_1507[k] = f_20 * osg_895[k]
                    + f_3 * pc_z[k] * qsg_1075[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, t_1511, pa_x, pc_x, pc_y, osh0_1508, \
                         osh0_1509, osh0_1511, osg_914, osh1_1508, osh1_1509, osh1_1511, \
                         qsg_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = pa_x[k] * osh0_1508[k]
                    - f_8 * pc_x[k] * osh1_1508[k];

        t_1509[k] = pa_x[k] * osh0_1509[k]
                    - f_8 * pc_x[k] * osh1_1509[k];

        t_1510[k] = f_21 * osg_914[k]
                    + f_3 * pc_y[k] * qsg_1079[k];

        t_1511[k] = pa_x[k] * osh0_1511[k]
                    - f_8 * pc_x[k] * osh1_1511[k];
    }

#pragma omp simd aligned(t_1512, t_1513, t_1514, pa_x, pc_x, pc_y, pc_z, osh0_1512, osg_900, \
                         osg_915, osg_1080, osh1_1512, qsg_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = pa_x[k] * osh0_1512[k]
                    + f_20 * osg_1080[k]
                    - f_8 * pc_x[k] * osh1_1512[k];

        t_1513[k] = f_20 * osg_915[k]
                    + f_3 * pc_y[k] * qsg_1080[k];

        t_1514[k] = f_21 * osg_900[k]
                    + f_3 * pc_z[k] * qsg_1080[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, pa_x, pc_x, pc_y, osh0_1515, osh0_1517, \
                         osg_917, osg_1083, osg_1085, osh1_1515, osh1_1517, \
                         qsg_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = pa_x[k] * osh0_1515[k]
                    + f_11 * osg_1083[k]
                    - f_8 * pc_x[k] * osh1_1515[k];

        t_1516[k] = f_20 * osg_917[k]
                    + f_3 * pc_y[k] * qsg_1082[k];

        t_1517[k] = pa_x[k] * osh0_1517[k]
                    + f_11 * osg_1085[k]
                    - f_8 * pc_x[k] * osh1_1517[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsf0, const size_t qsf1,
                                                           const size_t qsg, const size_t ncols,
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
    const auto f_12 = 5.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_1365 = buffer.data(osh0 + 1365);
    const auto *osh0_1370 = buffer.data(osh0 + 1370);
    const auto *osh0_1374 = buffer.data(osh0 + 1374);
    const auto *osh0_1518 = buffer.data(osh0 + 1518);
    const auto *osh0_1521 = buffer.data(osh0 + 1521);
    const auto *osh0_1527 = buffer.data(osh0 + 1527);
    const auto *osh0_1529 = buffer.data(osh0 + 1529);
    const auto *osh0_1530 = buffer.data(osh0 + 1530);
    const auto *osh0_1532 = buffer.data(osh0 + 1532);
    const auto *osh0_1533 = buffer.data(osh0 + 1533);
    const auto *osh0_1536 = buffer.data(osh0 + 1536);
    const auto *osh0_1538 = buffer.data(osh0 + 1538);
    const auto *osh0_1539 = buffer.data(osh0 + 1539);
    const auto *osh0_1542 = buffer.data(osh0 + 1542);
    const auto *osh0_1548 = buffer.data(osh0 + 1548);
    const auto *osh0_1550 = buffer.data(osh0 + 1550);
    const auto *osh0_1551 = buffer.data(osh0 + 1551);
    const auto *osh0_1553 = buffer.data(osh0 + 1553);
    const auto *osh0_1554 = buffer.data(osh0 + 1554);
    const auto *osh0_1557 = buffer.data(osh0 + 1557);
    const auto *osh0_1559 = buffer.data(osh0 + 1559);
    const auto *osh0_1560 = buffer.data(osh0 + 1560);
    const auto *osh0_1563 = buffer.data(osh0 + 1563);
    const auto *osh0_1569 = buffer.data(osh0 + 1569);
    const auto *osh0_1571 = buffer.data(osh0 + 1571);
    const auto *osh0_1572 = buffer.data(osh0 + 1572);
    const auto *osh0_1574 = buffer.data(osh0 + 1574);
    const auto *osh0_1575 = buffer.data(osh0 + 1575);
    const auto *osh0_1578 = buffer.data(osh0 + 1578);
    const auto *osh0_1580 = buffer.data(osh0 + 1580);
    const auto *osh0_1581 = buffer.data(osh0 + 1581);
    const auto *osh0_1584 = buffer.data(osh0 + 1584);
    const auto *osh0_1590 = buffer.data(osh0 + 1590);
    const auto *osh0_1592 = buffer.data(osh0 + 1592);
    const auto *osh0_1593 = buffer.data(osh0 + 1593);
    const auto *osh0_1595 = buffer.data(osh0 + 1595);
    const auto *osh0_1599 = buffer.data(osh0 + 1599);
    const auto *osh0_1602 = buffer.data(osh0 + 1602);
    const auto *osh0_1611 = buffer.data(osh0 + 1611);
    const auto *osh0_1613 = buffer.data(osh0 + 1613);
    const auto *osh0_1614 = buffer.data(osh0 + 1614);
    const auto *osh0_1616 = buffer.data(osh0 + 1616);
    const auto *osh0_1617 = buffer.data(osh0 + 1617);
    const auto *osh0_1622 = buffer.data(osh0 + 1622);
    const auto *osh0_1626 = buffer.data(osh0 + 1626);
    const auto *osh0_1632 = buffer.data(osh0 + 1632);
    const auto *osh0_1633 = buffer.data(osh0 + 1633);
    const auto *osh0_1634 = buffer.data(osh0 + 1634);
    const auto *osh0_1635 = buffer.data(osh0 + 1635);
    const auto *osh0_1637 = buffer.data(osh0 + 1637);

    const auto *osg_903 = buffer.data(osg + 903);
    const auto *osg_910 = buffer.data(osg + 910);
    const auto *osg_915 = buffer.data(osg + 915);
    const auto *osg_918 = buffer.data(osg + 918);
    const auto *osg_920 = buffer.data(osg + 920);
    const auto *osg_925 = buffer.data(osg + 925);
    const auto *osg_929 = buffer.data(osg + 929);
    const auto *osg_930 = buffer.data(osg + 930);
    const auto *osg_932 = buffer.data(osg + 932);
    const auto *osg_933 = buffer.data(osg + 933);
    const auto *osg_935 = buffer.data(osg + 935);
    const auto *osg_940 = buffer.data(osg + 940);
    const auto *osg_944 = buffer.data(osg + 944);
    const auto *osg_945 = buffer.data(osg + 945);
    const auto *osg_947 = buffer.data(osg + 947);
    const auto *osg_948 = buffer.data(osg + 948);
    const auto *osg_950 = buffer.data(osg + 950);
    const auto *osg_955 = buffer.data(osg + 955);
    const auto *osg_959 = buffer.data(osg + 959);
    const auto *osg_960 = buffer.data(osg + 960);
    const auto *osg_962 = buffer.data(osg + 962);
    const auto *osg_963 = buffer.data(osg + 963);
    const auto *osg_965 = buffer.data(osg + 965);
    const auto *osg_970 = buffer.data(osg + 970);
    const auto *osg_974 = buffer.data(osg + 974);
    const auto *osg_975 = buffer.data(osg + 975);
    const auto *osg_977 = buffer.data(osg + 977);
    const auto *osg_980 = buffer.data(osg + 980);
    const auto *osg_989 = buffer.data(osg + 989);
    const auto *osg_1086 = buffer.data(osg + 1086);
    const auto *osg_1089 = buffer.data(osg + 1089);
    const auto *osg_1090 = buffer.data(osg + 1090);
    const auto *osg_1091 = buffer.data(osg + 1091);
    const auto *osg_1092 = buffer.data(osg + 1092);
    const auto *osg_1093 = buffer.data(osg + 1093);
    const auto *osg_1094 = buffer.data(osg + 1094);
    const auto *osg_1095 = buffer.data(osg + 1095);
    const auto *osg_1098 = buffer.data(osg + 1098);
    const auto *osg_1100 = buffer.data(osg + 1100);
    const auto *osg_1101 = buffer.data(osg + 1101);
    const auto *osg_1104 = buffer.data(osg + 1104);
    const auto *osg_1105 = buffer.data(osg + 1105);
    const auto *osg_1106 = buffer.data(osg + 1106);
    const auto *osg_1107 = buffer.data(osg + 1107);
    const auto *osg_1108 = buffer.data(osg + 1108);
    const auto *osg_1109 = buffer.data(osg + 1109);
    const auto *osg_1110 = buffer.data(osg + 1110);
    const auto *osg_1113 = buffer.data(osg + 1113);
    const auto *osg_1115 = buffer.data(osg + 1115);
    const auto *osg_1116 = buffer.data(osg + 1116);
    const auto *osg_1119 = buffer.data(osg + 1119);
    const auto *osg_1120 = buffer.data(osg + 1120);
    const auto *osg_1121 = buffer.data(osg + 1121);
    const auto *osg_1122 = buffer.data(osg + 1122);
    const auto *osg_1123 = buffer.data(osg + 1123);
    const auto *osg_1124 = buffer.data(osg + 1124);
    const auto *osg_1125 = buffer.data(osg + 1125);
    const auto *osg_1128 = buffer.data(osg + 1128);
    const auto *osg_1130 = buffer.data(osg + 1130);
    const auto *osg_1131 = buffer.data(osg + 1131);
    const auto *osg_1134 = buffer.data(osg + 1134);
    const auto *osg_1135 = buffer.data(osg + 1135);
    const auto *osg_1136 = buffer.data(osg + 1136);
    const auto *osg_1137 = buffer.data(osg + 1137);
    const auto *osg_1138 = buffer.data(osg + 1138);
    const auto *osg_1139 = buffer.data(osg + 1139);
    const auto *osg_1143 = buffer.data(osg + 1143);
    const auto *osg_1146 = buffer.data(osg + 1146);
    const auto *osg_1150 = buffer.data(osg + 1150);
    const auto *osg_1151 = buffer.data(osg + 1151);
    const auto *osg_1152 = buffer.data(osg + 1152);
    const auto *osg_1153 = buffer.data(osg + 1153);
    const auto *osg_1154 = buffer.data(osg + 1154);
    const auto *osg_1155 = buffer.data(osg + 1155);
    const auto *osg_1160 = buffer.data(osg + 1160);
    const auto *osg_1164 = buffer.data(osg + 1164);
    const auto *osg_1165 = buffer.data(osg + 1165);
    const auto *osg_1166 = buffer.data(osg + 1166);
    const auto *osg_1167 = buffer.data(osg + 1167);
    const auto *osg_1169 = buffer.data(osg + 1169);

    const auto *osh1_1365 = buffer.data(osh1 + 1365);
    const auto *osh1_1370 = buffer.data(osh1 + 1370);
    const auto *osh1_1374 = buffer.data(osh1 + 1374);
    const auto *osh1_1518 = buffer.data(osh1 + 1518);
    const auto *osh1_1521 = buffer.data(osh1 + 1521);
    const auto *osh1_1527 = buffer.data(osh1 + 1527);
    const auto *osh1_1529 = buffer.data(osh1 + 1529);
    const auto *osh1_1530 = buffer.data(osh1 + 1530);
    const auto *osh1_1532 = buffer.data(osh1 + 1532);
    const auto *osh1_1533 = buffer.data(osh1 + 1533);
    const auto *osh1_1536 = buffer.data(osh1 + 1536);
    const auto *osh1_1538 = buffer.data(osh1 + 1538);
    const auto *osh1_1539 = buffer.data(osh1 + 1539);
    const auto *osh1_1542 = buffer.data(osh1 + 1542);
    const auto *osh1_1548 = buffer.data(osh1 + 1548);
    const auto *osh1_1550 = buffer.data(osh1 + 1550);
    const auto *osh1_1551 = buffer.data(osh1 + 1551);
    const auto *osh1_1553 = buffer.data(osh1 + 1553);
    const auto *osh1_1554 = buffer.data(osh1 + 1554);
    const auto *osh1_1557 = buffer.data(osh1 + 1557);
    const auto *osh1_1559 = buffer.data(osh1 + 1559);
    const auto *osh1_1560 = buffer.data(osh1 + 1560);
    const auto *osh1_1563 = buffer.data(osh1 + 1563);
    const auto *osh1_1569 = buffer.data(osh1 + 1569);
    const auto *osh1_1571 = buffer.data(osh1 + 1571);
    const auto *osh1_1572 = buffer.data(osh1 + 1572);
    const auto *osh1_1574 = buffer.data(osh1 + 1574);
    const auto *osh1_1575 = buffer.data(osh1 + 1575);
    const auto *osh1_1578 = buffer.data(osh1 + 1578);
    const auto *osh1_1580 = buffer.data(osh1 + 1580);
    const auto *osh1_1581 = buffer.data(osh1 + 1581);
    const auto *osh1_1584 = buffer.data(osh1 + 1584);
    const auto *osh1_1590 = buffer.data(osh1 + 1590);
    const auto *osh1_1592 = buffer.data(osh1 + 1592);
    const auto *osh1_1593 = buffer.data(osh1 + 1593);
    const auto *osh1_1595 = buffer.data(osh1 + 1595);
    const auto *osh1_1599 = buffer.data(osh1 + 1599);
    const auto *osh1_1602 = buffer.data(osh1 + 1602);
    const auto *osh1_1611 = buffer.data(osh1 + 1611);
    const auto *osh1_1613 = buffer.data(osh1 + 1613);
    const auto *osh1_1614 = buffer.data(osh1 + 1614);
    const auto *osh1_1616 = buffer.data(osh1 + 1616);
    const auto *osh1_1617 = buffer.data(osh1 + 1617);
    const auto *osh1_1622 = buffer.data(osh1 + 1622);
    const auto *osh1_1626 = buffer.data(osh1 + 1626);
    const auto *osh1_1632 = buffer.data(osh1 + 1632);
    const auto *osh1_1633 = buffer.data(osh1 + 1633);
    const auto *osh1_1634 = buffer.data(osh1 + 1634);
    const auto *osh1_1635 = buffer.data(osh1 + 1635);
    const auto *osh1_1637 = buffer.data(osh1 + 1637);

    const auto *qsf0_770 = buffer.data(qsf0 + 770);
    const auto *qsf0_771 = buffer.data(qsf0 + 771);
    const auto *qsf0_772 = buffer.data(qsf0 + 772);
    const auto *qsf0_780 = buffer.data(qsf0 + 780);
    const auto *qsf0_781 = buffer.data(qsf0 + 781);

    const auto *qsf1_770 = buffer.data(qsf1 + 770);
    const auto *qsf1_771 = buffer.data(qsf1 + 771);
    const auto *qsf1_772 = buffer.data(qsf1 + 772);
    const auto *qsf1_780 = buffer.data(qsf1 + 780);
    const auto *qsf1_781 = buffer.data(qsf1 + 781);

    const auto *qsg_1083 = buffer.data(qsg + 1083);
    const auto *qsg_1085 = buffer.data(qsg + 1085);
    const auto *qsg_1090 = buffer.data(qsg + 1090);
    const auto *qsg_1091 = buffer.data(qsg + 1091);
    const auto *qsg_1092 = buffer.data(qsg + 1092);
    const auto *qsg_1093 = buffer.data(qsg + 1093);
    const auto *qsg_1094 = buffer.data(qsg + 1094);
    const auto *qsg_1095 = buffer.data(qsg + 1095);
    const auto *qsg_1097 = buffer.data(qsg + 1097);
    const auto *qsg_1098 = buffer.data(qsg + 1098);
    const auto *qsg_1100 = buffer.data(qsg + 1100);
    const auto *qsg_1105 = buffer.data(qsg + 1105);
    const auto *qsg_1106 = buffer.data(qsg + 1106);
    const auto *qsg_1107 = buffer.data(qsg + 1107);
    const auto *qsg_1108 = buffer.data(qsg + 1108);
    const auto *qsg_1109 = buffer.data(qsg + 1109);
    const auto *qsg_1110 = buffer.data(qsg + 1110);
    const auto *qsg_1112 = buffer.data(qsg + 1112);
    const auto *qsg_1113 = buffer.data(qsg + 1113);
    const auto *qsg_1115 = buffer.data(qsg + 1115);
    const auto *qsg_1120 = buffer.data(qsg + 1120);
    const auto *qsg_1121 = buffer.data(qsg + 1121);
    const auto *qsg_1122 = buffer.data(qsg + 1122);
    const auto *qsg_1123 = buffer.data(qsg + 1123);
    const auto *qsg_1124 = buffer.data(qsg + 1124);
    const auto *qsg_1125 = buffer.data(qsg + 1125);
    const auto *qsg_1127 = buffer.data(qsg + 1127);
    const auto *qsg_1128 = buffer.data(qsg + 1128);
    const auto *qsg_1130 = buffer.data(qsg + 1130);
    const auto *qsg_1135 = buffer.data(qsg + 1135);
    const auto *qsg_1136 = buffer.data(qsg + 1136);
    const auto *qsg_1137 = buffer.data(qsg + 1137);
    const auto *qsg_1138 = buffer.data(qsg + 1138);
    const auto *qsg_1139 = buffer.data(qsg + 1139);
    const auto *qsg_1140 = buffer.data(qsg + 1140);
    const auto *qsg_1142 = buffer.data(qsg + 1142);
    const auto *qsg_1143 = buffer.data(qsg + 1143);
    const auto *qsg_1145 = buffer.data(qsg + 1145);
    const auto *qsg_1150 = buffer.data(qsg + 1150);
    const auto *qsg_1151 = buffer.data(qsg + 1151);
    const auto *qsg_1152 = buffer.data(qsg + 1152);
    const auto *qsg_1153 = buffer.data(qsg + 1153);
    const auto *qsg_1154 = buffer.data(qsg + 1154);
    const auto *qsg_1155 = buffer.data(qsg + 1155);
    const auto *qsg_1156 = buffer.data(qsg + 1156);
    const auto *qsg_1157 = buffer.data(qsg + 1157);
    const auto *qsg_1158 = buffer.data(qsg + 1158);
    const auto *qsg_1159 = buffer.data(qsg + 1159);
    const auto *qsg_1160 = buffer.data(qsg + 1160);
    const auto *qsg_1164 = buffer.data(qsg + 1164);
    const auto *qsg_1165 = buffer.data(qsg + 1165);
    const auto *qsg_1166 = buffer.data(qsg + 1166);
    const auto *qsg_1167 = buffer.data(qsg + 1167);
    const auto *qsg_1169 = buffer.data(qsg + 1169);
    const auto *qsg_1170 = buffer.data(qsg + 1170);
    const auto *qsg_1171 = buffer.data(qsg + 1171);

#pragma omp simd aligned(t_1518, t_1519, t_1520, pa_x, pc_x, pc_y, pc_z, osh0_1518, osg_903, \
                         osg_920, osg_1086, osh1_1518, qsg_1083, \
                         qsg_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = pa_x[k] * osh0_1518[k]
                    + f_10 * osg_1086[k]
                    - f_8 * pc_x[k] * osh1_1518[k];

        t_1519[k] = f_21 * osg_903[k]
                    + f_3 * pc_z[k] * qsg_1083[k];

        t_1520[k] = f_20 * osg_920[k]
                    + f_3 * pc_y[k] * qsg_1085[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, t_1524, pa_x, pc_x, osh0_1521, osg_1089, \
                         osg_1090, osg_1091, osg_1092, osh1_1521, qsg_1090, qsg_1091, \
                         qsg_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = pa_x[k] * osh0_1521[k]
                    + f_10 * osg_1089[k]
                    - f_8 * pc_x[k] * osh1_1521[k];

        t_1522[k] = f_9 * osg_1090[k]
                    + f_3 * pc_x[k] * qsg_1090[k];

        t_1523[k] = f_9 * osg_1091[k]
                    + f_3 * pc_x[k] * qsg_1091[k];

        t_1524[k] = f_9 * osg_1092[k]
                    + f_3 * pc_x[k] * qsg_1092[k];
    }

#pragma omp simd aligned(t_1525, t_1526, t_1527, t_1528, pa_x, pc_x, pc_z, osh0_1527, osg_910, \
                         osg_1093, osg_1094, osh1_1527, qsg_1090, qsg_1093, \
                         qsg_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1525[k] = f_9 * osg_1093[k]
                    + f_3 * pc_x[k] * qsg_1093[k];

        t_1526[k] = f_9 * osg_1094[k]
                    + f_3 * pc_x[k] * qsg_1094[k];

        t_1527[k] = pa_x[k] * osh0_1527[k]
                    - f_8 * pc_x[k] * osh1_1527[k];

        t_1528[k] = f_21 * osg_910[k]
                    + f_3 * pc_z[k] * qsg_1090[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, pa_x, pc_x, pc_y, osh0_1529, \
                         osh0_1530, osh0_1532, osg_929, osh1_1529, osh1_1530, osh1_1532, \
                         qsg_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = pa_x[k] * osh0_1529[k]
                    - f_8 * pc_x[k] * osh1_1529[k];

        t_1530[k] = pa_x[k] * osh0_1530[k]
                    - f_8 * pc_x[k] * osh1_1530[k];

        t_1531[k] = f_20 * osg_929[k]
                    + f_3 * pc_y[k] * qsg_1094[k];

        t_1532[k] = pa_x[k] * osh0_1532[k]
                    - f_8 * pc_x[k] * osh1_1532[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, pa_x, pc_x, pc_y, pc_z, osh0_1533, osg_915, \
                         osg_930, osg_1095, osh1_1533, qsg_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = pa_x[k] * osh0_1533[k]
                    + f_20 * osg_1095[k]
                    - f_8 * pc_x[k] * osh1_1533[k];

        t_1534[k] = f_18 * osg_930[k]
                    + f_3 * pc_y[k] * qsg_1095[k];

        t_1535[k] = f_19 * osg_915[k]
                    + f_3 * pc_z[k] * qsg_1095[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, pa_x, pc_x, pc_y, osh0_1536, osh0_1538, \
                         osg_932, osg_1098, osg_1100, osh1_1536, osh1_1538, \
                         qsg_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = pa_x[k] * osh0_1536[k]
                    + f_11 * osg_1098[k]
                    - f_8 * pc_x[k] * osh1_1536[k];

        t_1537[k] = f_18 * osg_932[k]
                    + f_3 * pc_y[k] * qsg_1097[k];

        t_1538[k] = pa_x[k] * osh0_1538[k]
                    + f_11 * osg_1100[k]
                    - f_8 * pc_x[k] * osh1_1538[k];
    }

#pragma omp simd aligned(t_1539, t_1540, t_1541, pa_x, pc_x, pc_y, pc_z, osh0_1539, osg_918, \
                         osg_935, osg_1101, osh1_1539, qsg_1098, \
                         qsg_1100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1539[k] = pa_x[k] * osh0_1539[k]
                    + f_10 * osg_1101[k]
                    - f_8 * pc_x[k] * osh1_1539[k];

        t_1540[k] = f_19 * osg_918[k]
                    + f_3 * pc_z[k] * qsg_1098[k];

        t_1541[k] = f_18 * osg_935[k]
                    + f_3 * pc_y[k] * qsg_1100[k];
    }

#pragma omp simd aligned(t_1542, t_1543, t_1544, t_1545, pa_x, pc_x, osh0_1542, osg_1104, \
                         osg_1105, osg_1106, osg_1107, osh1_1542, qsg_1105, qsg_1106, \
                         qsg_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1542[k] = pa_x[k] * osh0_1542[k]
                    + f_10 * osg_1104[k]
                    - f_8 * pc_x[k] * osh1_1542[k];

        t_1543[k] = f_9 * osg_1105[k]
                    + f_3 * pc_x[k] * qsg_1105[k];

        t_1544[k] = f_9 * osg_1106[k]
                    + f_3 * pc_x[k] * qsg_1106[k];

        t_1545[k] = f_9 * osg_1107[k]
                    + f_3 * pc_x[k] * qsg_1107[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pa_x, pc_x, pc_z, osh0_1548, osg_925, \
                         osg_1108, osg_1109, osh1_1548, qsg_1105, qsg_1108, \
                         qsg_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_9 * osg_1108[k]
                    + f_3 * pc_x[k] * qsg_1108[k];

        t_1547[k] = f_9 * osg_1109[k]
                    + f_3 * pc_x[k] * qsg_1109[k];

        t_1548[k] = pa_x[k] * osh0_1548[k]
                    - f_8 * pc_x[k] * osh1_1548[k];

        t_1549[k] = f_19 * osg_925[k]
                    + f_3 * pc_z[k] * qsg_1105[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, t_1553, pa_x, pc_x, pc_y, osh0_1550, \
                         osh0_1551, osh0_1553, osg_944, osh1_1550, osh1_1551, osh1_1553, \
                         qsg_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = pa_x[k] * osh0_1550[k]
                    - f_8 * pc_x[k] * osh1_1550[k];

        t_1551[k] = pa_x[k] * osh0_1551[k]
                    - f_8 * pc_x[k] * osh1_1551[k];

        t_1552[k] = f_18 * osg_944[k]
                    + f_3 * pc_y[k] * qsg_1109[k];

        t_1553[k] = pa_x[k] * osh0_1553[k]
                    - f_8 * pc_x[k] * osh1_1553[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, pa_x, pc_x, pc_y, pc_z, osh0_1554, osg_930, \
                         osg_945, osg_1110, osh1_1554, qsg_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = pa_x[k] * osh0_1554[k]
                    + f_20 * osg_1110[k]
                    - f_8 * pc_x[k] * osh1_1554[k];

        t_1555[k] = f_11 * osg_945[k]
                    + f_3 * pc_y[k] * qsg_1110[k];

        t_1556[k] = f_17 * osg_930[k]
                    + f_3 * pc_z[k] * qsg_1110[k];
    }

#pragma omp simd aligned(t_1557, t_1558, t_1559, pa_x, pc_x, pc_y, osh0_1557, osh0_1559, \
                         osg_947, osg_1113, osg_1115, osh1_1557, osh1_1559, \
                         qsg_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1557[k] = pa_x[k] * osh0_1557[k]
                    + f_11 * osg_1113[k]
                    - f_8 * pc_x[k] * osh1_1557[k];

        t_1558[k] = f_11 * osg_947[k]
                    + f_3 * pc_y[k] * qsg_1112[k];

        t_1559[k] = pa_x[k] * osh0_1559[k]
                    + f_11 * osg_1115[k]
                    - f_8 * pc_x[k] * osh1_1559[k];
    }

#pragma omp simd aligned(t_1560, t_1561, t_1562, pa_x, pc_x, pc_y, pc_z, osh0_1560, osg_933, \
                         osg_950, osg_1116, osh1_1560, qsg_1113, \
                         qsg_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1560[k] = pa_x[k] * osh0_1560[k]
                    + f_10 * osg_1116[k]
                    - f_8 * pc_x[k] * osh1_1560[k];

        t_1561[k] = f_17 * osg_933[k]
                    + f_3 * pc_z[k] * qsg_1113[k];

        t_1562[k] = f_11 * osg_950[k]
                    + f_3 * pc_y[k] * qsg_1115[k];
    }

#pragma omp simd aligned(t_1563, t_1564, t_1565, t_1566, pa_x, pc_x, osh0_1563, osg_1119, \
                         osg_1120, osg_1121, osg_1122, osh1_1563, qsg_1120, qsg_1121, \
                         qsg_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1563[k] = pa_x[k] * osh0_1563[k]
                    + f_10 * osg_1119[k]
                    - f_8 * pc_x[k] * osh1_1563[k];

        t_1564[k] = f_9 * osg_1120[k]
                    + f_3 * pc_x[k] * qsg_1120[k];

        t_1565[k] = f_9 * osg_1121[k]
                    + f_3 * pc_x[k] * qsg_1121[k];

        t_1566[k] = f_9 * osg_1122[k]
                    + f_3 * pc_x[k] * qsg_1122[k];
    }

#pragma omp simd aligned(t_1567, t_1568, t_1569, t_1570, pa_x, pc_x, pc_z, osh0_1569, osg_940, \
                         osg_1123, osg_1124, osh1_1569, qsg_1120, qsg_1123, \
                         qsg_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1567[k] = f_9 * osg_1123[k]
                    + f_3 * pc_x[k] * qsg_1123[k];

        t_1568[k] = f_9 * osg_1124[k]
                    + f_3 * pc_x[k] * qsg_1124[k];

        t_1569[k] = pa_x[k] * osh0_1569[k]
                    - f_8 * pc_x[k] * osh1_1569[k];

        t_1570[k] = f_17 * osg_940[k]
                    + f_3 * pc_z[k] * qsg_1120[k];
    }

#pragma omp simd aligned(t_1571, t_1572, t_1573, t_1574, pa_x, pc_x, pc_y, osh0_1571, \
                         osh0_1572, osh0_1574, osg_959, osh1_1571, osh1_1572, osh1_1574, \
                         qsg_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1571[k] = pa_x[k] * osh0_1571[k]
                    - f_8 * pc_x[k] * osh1_1571[k];

        t_1572[k] = pa_x[k] * osh0_1572[k]
                    - f_8 * pc_x[k] * osh1_1572[k];

        t_1573[k] = f_11 * osg_959[k]
                    + f_3 * pc_y[k] * qsg_1124[k];

        t_1574[k] = pa_x[k] * osh0_1574[k]
                    - f_8 * pc_x[k] * osh1_1574[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, pa_x, pc_x, pc_y, pc_z, osh0_1575, osg_945, \
                         osg_960, osg_1125, osh1_1575, qsg_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = pa_x[k] * osh0_1575[k]
                    + f_20 * osg_1125[k]
                    - f_8 * pc_x[k] * osh1_1575[k];

        t_1576[k] = f_10 * osg_960[k]
                    + f_3 * pc_y[k] * qsg_1125[k];

        t_1577[k] = f_16 * osg_945[k]
                    + f_3 * pc_z[k] * qsg_1125[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pa_x, pc_x, pc_y, osh0_1578, osh0_1580, \
                         osg_962, osg_1128, osg_1130, osh1_1578, osh1_1580, \
                         qsg_1127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pa_x[k] * osh0_1578[k]
                    + f_11 * osg_1128[k]
                    - f_8 * pc_x[k] * osh1_1578[k];

        t_1579[k] = f_10 * osg_962[k]
                    + f_3 * pc_y[k] * qsg_1127[k];

        t_1580[k] = pa_x[k] * osh0_1580[k]
                    + f_11 * osg_1130[k]
                    - f_8 * pc_x[k] * osh1_1580[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pa_x, pc_x, pc_y, pc_z, osh0_1581, osg_948, \
                         osg_965, osg_1131, osh1_1581, qsg_1128, \
                         qsg_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = pa_x[k] * osh0_1581[k]
                    + f_10 * osg_1131[k]
                    - f_8 * pc_x[k] * osh1_1581[k];

        t_1582[k] = f_16 * osg_948[k]
                    + f_3 * pc_z[k] * qsg_1128[k];

        t_1583[k] = f_10 * osg_965[k]
                    + f_3 * pc_y[k] * qsg_1130[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, t_1587, pa_x, pc_x, osh0_1584, osg_1134, \
                         osg_1135, osg_1136, osg_1137, osh1_1584, qsg_1135, qsg_1136, \
                         qsg_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = pa_x[k] * osh0_1584[k]
                    + f_10 * osg_1134[k]
                    - f_8 * pc_x[k] * osh1_1584[k];

        t_1585[k] = f_9 * osg_1135[k]
                    + f_3 * pc_x[k] * qsg_1135[k];

        t_1586[k] = f_9 * osg_1136[k]
                    + f_3 * pc_x[k] * qsg_1136[k];

        t_1587[k] = f_9 * osg_1137[k]
                    + f_3 * pc_x[k] * qsg_1137[k];
    }

#pragma omp simd aligned(t_1588, t_1589, t_1590, t_1591, pa_x, pc_x, pc_z, osh0_1590, osg_955, \
                         osg_1138, osg_1139, osh1_1590, qsg_1135, qsg_1138, \
                         qsg_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1588[k] = f_9 * osg_1138[k]
                    + f_3 * pc_x[k] * qsg_1138[k];

        t_1589[k] = f_9 * osg_1139[k]
                    + f_3 * pc_x[k] * qsg_1139[k];

        t_1590[k] = pa_x[k] * osh0_1590[k]
                    - f_8 * pc_x[k] * osh1_1590[k];

        t_1591[k] = f_16 * osg_955[k]
                    + f_3 * pc_z[k] * qsg_1135[k];
    }

#pragma omp simd aligned(t_1592, t_1593, t_1594, t_1595, pa_x, pc_x, pc_y, osh0_1592, \
                         osh0_1593, osh0_1595, osg_974, osh1_1592, osh1_1593, osh1_1595, \
                         qsg_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1592[k] = pa_x[k] * osh0_1592[k]
                    - f_8 * pc_x[k] * osh1_1592[k];

        t_1593[k] = pa_x[k] * osh0_1593[k]
                    - f_8 * pc_x[k] * osh1_1593[k];

        t_1594[k] = f_10 * osg_974[k]
                    + f_3 * pc_y[k] * qsg_1139[k];

        t_1595[k] = pa_x[k] * osh0_1595[k]
                    - f_8 * pc_x[k] * osh1_1595[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, pa_y, pc_y, pc_z, osh0_1365, osg_960, \
                         osg_975, osh1_1365, qsg_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = pa_y[k] * osh0_1365[k]
                    - f_8 * pc_y[k] * osh1_1365[k];

        t_1597[k] = f_9 * osg_975[k]
                    + f_3 * pc_y[k] * qsg_1140[k];

        t_1598[k] = f_15 * osg_960[k]
                    + f_3 * pc_z[k] * qsg_1140[k];
    }

#pragma omp simd aligned(t_1599, t_1600, t_1601, pa_x, pa_y, pc_x, pc_y, osh0_1370, osh0_1599, \
                         osg_977, osg_1143, osh1_1370, osh1_1599, \
                         qsg_1142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = pa_x[k] * osh0_1599[k]
                    + f_11 * osg_1143[k]
                    - f_8 * pc_x[k] * osh1_1599[k];

        t_1600[k] = f_9 * osg_977[k]
                    + f_3 * pc_y[k] * qsg_1142[k];

        t_1601[k] = pa_y[k] * osh0_1370[k]
                    - f_8 * pc_y[k] * osh1_1370[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, pa_x, pc_x, pc_y, pc_z, osh0_1602, osg_963, \
                         osg_980, osg_1146, osh1_1602, qsg_1143, \
                         qsg_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = pa_x[k] * osh0_1602[k]
                    + f_10 * osg_1146[k]
                    - f_8 * pc_x[k] * osh1_1602[k];

        t_1603[k] = f_15 * osg_963[k]
                    + f_3 * pc_z[k] * qsg_1143[k];

        t_1604[k] = f_9 * osg_980[k]
                    + f_3 * pc_y[k] * qsg_1145[k];
    }

#pragma omp simd aligned(t_1605, t_1606, t_1607, t_1608, pa_y, pc_x, pc_y, osh0_1374, \
                         osg_1150, osg_1151, osg_1152, osh1_1374, qsg_1150, qsg_1151, \
                         qsg_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1605[k] = pa_y[k] * osh0_1374[k]
                    - f_8 * pc_y[k] * osh1_1374[k];

        t_1606[k] = f_9 * osg_1150[k]
                    + f_3 * pc_x[k] * qsg_1150[k];

        t_1607[k] = f_9 * osg_1151[k]
                    + f_3 * pc_x[k] * qsg_1151[k];

        t_1608[k] = f_9 * osg_1152[k]
                    + f_3 * pc_x[k] * qsg_1152[k];
    }

#pragma omp simd aligned(t_1609, t_1610, t_1611, t_1612, pa_x, pc_x, pc_z, osh0_1611, osg_970, \
                         osg_1153, osg_1154, osh1_1611, qsg_1150, qsg_1153, \
                         qsg_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1609[k] = f_9 * osg_1153[k]
                    + f_3 * pc_x[k] * qsg_1153[k];

        t_1610[k] = f_9 * osg_1154[k]
                    + f_3 * pc_x[k] * qsg_1154[k];

        t_1611[k] = pa_x[k] * osh0_1611[k]
                    - f_8 * pc_x[k] * osh1_1611[k];

        t_1612[k] = f_15 * osg_970[k]
                    + f_3 * pc_z[k] * qsg_1150[k];
    }

#pragma omp simd aligned(t_1613, t_1614, t_1615, t_1616, pa_x, pc_x, pc_y, osh0_1613, \
                         osh0_1614, osh0_1616, osg_989, osh1_1613, osh1_1614, osh1_1616, \
                         qsg_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1613[k] = pa_x[k] * osh0_1613[k]
                    - f_8 * pc_x[k] * osh1_1613[k];

        t_1614[k] = pa_x[k] * osh0_1614[k]
                    - f_8 * pc_x[k] * osh1_1614[k];

        t_1615[k] = f_9 * osg_989[k]
                    + f_3 * pc_y[k] * qsg_1154[k];

        t_1616[k] = pa_x[k] * osh0_1616[k]
                    - f_8 * pc_x[k] * osh1_1616[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, t_1620, pa_x, pc_x, pc_y, pc_z, osh0_1617, \
                         osg_975, osg_1155, osh1_1617, qsf0_770, qsf1_770, qsg_1155, \
                         qsg_1156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = pa_x[k] * osh0_1617[k]
                    + f_20 * osg_1155[k]
                    - f_8 * pc_x[k] * osh1_1617[k];

        t_1618[k] = f_3 * pc_y[k] * qsg_1155[k];

        t_1619[k] = f_12 * osg_975[k]
                    + f_3 * pc_z[k] * qsg_1155[k];

        t_1620[k] = f_4 * qsf0_770[k]
                    - f_5 * qsf1_770[k]
                    + f_3 * pc_y[k] * qsg_1156[k];
    }

#pragma omp simd aligned(t_1621, t_1622, t_1623, pa_x, pc_x, pc_y, osh0_1622, osg_1160, \
                         osh1_1622, qsf0_771, qsf1_771, qsg_1157, \
                         qsg_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1621[k] = f_3 * pc_y[k] * qsg_1157[k];

        t_1622[k] = pa_x[k] * osh0_1622[k]
                    + f_11 * osg_1160[k]
                    - f_8 * pc_x[k] * osh1_1622[k];

        t_1623[k] = f_6 * qsf0_771[k]
                    - f_7 * qsf1_771[k]
                    + f_3 * pc_y[k] * qsg_1158[k];
    }

#pragma omp simd aligned(t_1624, t_1625, t_1626, t_1627, pa_x, pc_x, pc_y, osh0_1626, \
                         osg_1164, osg_1165, osh1_1626, qsf0_772, qsf1_772, qsg_1159, \
                         qsg_1160, qsg_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1624[k] = f_4 * qsf0_772[k]
                    - f_5 * qsf1_772[k]
                    + f_3 * pc_y[k] * qsg_1159[k];

        t_1625[k] = f_3 * pc_y[k] * qsg_1160[k];

        t_1626[k] = pa_x[k] * osh0_1626[k]
                    + f_10 * osg_1164[k]
                    - f_8 * pc_x[k] * osh1_1626[k];

        t_1627[k] = f_9 * osg_1165[k]
                    + f_3 * pc_x[k] * qsg_1165[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, t_1631, pc_x, pc_y, osg_1166, osg_1167, \
                         osg_1169, qsg_1164, qsg_1166, qsg_1167, \
                         qsg_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_9 * osg_1166[k]
                    + f_3 * pc_x[k] * qsg_1166[k];

        t_1629[k] = f_9 * osg_1167[k]
                    + f_3 * pc_x[k] * qsg_1167[k];

        t_1630[k] = f_3 * pc_y[k] * qsg_1164[k];

        t_1631[k] = f_9 * osg_1169[k]
                    + f_3 * pc_x[k] * qsg_1169[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, t_1635, pa_x, pc_x, osh0_1632, osh0_1633, \
                         osh0_1634, osh0_1635, osh1_1632, osh1_1633, osh1_1634, \
                         osh1_1635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = pa_x[k] * osh0_1632[k]
                    - f_8 * pc_x[k] * osh1_1632[k];

        t_1633[k] = pa_x[k] * osh0_1633[k]
                    - f_8 * pc_x[k] * osh1_1633[k];

        t_1634[k] = pa_x[k] * osh0_1634[k]
                    - f_8 * pc_x[k] * osh1_1634[k];

        t_1635[k] = pa_x[k] * osh0_1635[k]
                    - f_8 * pc_x[k] * osh1_1635[k];
    }

#pragma omp simd aligned(t_1636, t_1637, t_1638, t_1639, pa_x, pc_x, pc_y, osh0_1637, \
                         osh1_1637, qsf0_780, qsf0_781, qsf1_780, qsf1_781, qsg_1169, \
                         qsg_1170, qsg_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1636[k] = f_3 * pc_y[k] * qsg_1169[k];

        t_1637[k] = pa_x[k] * osh0_1637[k]
                    - f_8 * pc_x[k] * osh1_1637[k];

        t_1638[k] = f_1 * qsf0_780[k]
                    - f_2 * qsf1_780[k]
                    + f_3 * pc_x[k] * qsg_1170[k];

        t_1639[k] = f_13 * qsf0_781[k]
                    - f_14 * qsf1_781[k]
                    + f_3 * pc_x[k] * qsg_1171[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsf0, const size_t qsf1,
                                                           const size_t qsg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
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
    const auto f_12 = 5.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_1386 = buffer.data(osh0 + 1386);
    const auto *osh0_1387 = buffer.data(osh0 + 1387);
    const auto *osh0_1389 = buffer.data(osh0 + 1389);
    const auto *osh0_1392 = buffer.data(osh0 + 1392);
    const auto *osh0_1401 = buffer.data(osh0 + 1401);
    const auto *osh0_1403 = buffer.data(osh0 + 1403);
    const auto *osh0_1404 = buffer.data(osh0 + 1404);

    const auto *osg_1000 = buffer.data(osg + 1000);
    const auto *osg_1001 = buffer.data(osg + 1001);
    const auto *osg_1002 = buffer.data(osg + 1002);
    const auto *osg_1004 = buffer.data(osg + 1004);
    const auto *osg_1015 = buffer.data(osg + 1015);
    const auto *osg_1019 = buffer.data(osg + 1019);
    const auto *osg_1030 = buffer.data(osg + 1030);
    const auto *osg_1032 = buffer.data(osg + 1032);
    const auto *osg_1033 = buffer.data(osg + 1033);
    const auto *osg_1034 = buffer.data(osg + 1034);
    const auto *osg_1045 = buffer.data(osg + 1045);
    const auto *osg_1047 = buffer.data(osg + 1047);
    const auto *osg_1048 = buffer.data(osg + 1048);
    const auto *osg_1049 = buffer.data(osg + 1049);
    const auto *osg_1060 = buffer.data(osg + 1060);
    const auto *osg_1062 = buffer.data(osg + 1062);
    const auto *osg_1063 = buffer.data(osg + 1063);
    const auto *osg_1064 = buffer.data(osg + 1064);
    const auto *osg_1075 = buffer.data(osg + 1075);

    const auto *osh1_1386 = buffer.data(osh1 + 1386);
    const auto *osh1_1387 = buffer.data(osh1 + 1387);
    const auto *osh1_1389 = buffer.data(osh1 + 1389);
    const auto *osh1_1392 = buffer.data(osh1 + 1392);
    const auto *osh1_1401 = buffer.data(osh1 + 1401);
    const auto *osh1_1403 = buffer.data(osh1 + 1403);
    const auto *osh1_1404 = buffer.data(osh1 + 1404);

    const auto *qsf0_783 = buffer.data(qsf0 + 783);
    const auto *qsf0_785 = buffer.data(qsf0 + 785);
    const auto *qsf0_786 = buffer.data(qsf0 + 786);
    const auto *qsf0_787 = buffer.data(qsf0 + 787);
    const auto *qsf0_788 = buffer.data(qsf0 + 788);
    const auto *qsf0_789 = buffer.data(qsf0 + 789);
    const auto *qsf0_792 = buffer.data(qsf0 + 792);
    const auto *qsf0_794 = buffer.data(qsf0 + 794);
    const auto *qsf0_795 = buffer.data(qsf0 + 795);
    const auto *qsf0_797 = buffer.data(qsf0 + 797);
    const auto *qsf0_798 = buffer.data(qsf0 + 798);
    const auto *qsf0_799 = buffer.data(qsf0 + 799);
    const auto *qsf0_800 = buffer.data(qsf0 + 800);
    const auto *qsf0_801 = buffer.data(qsf0 + 801);
    const auto *qsf0_802 = buffer.data(qsf0 + 802);
    const auto *qsf0_803 = buffer.data(qsf0 + 803);
    const auto *qsf0_804 = buffer.data(qsf0 + 804);
    const auto *qsf0_805 = buffer.data(qsf0 + 805);
    const auto *qsf0_806 = buffer.data(qsf0 + 806);
    const auto *qsf0_807 = buffer.data(qsf0 + 807);
    const auto *qsf0_808 = buffer.data(qsf0 + 808);
    const auto *qsf0_809 = buffer.data(qsf0 + 809);
    const auto *qsf0_810 = buffer.data(qsf0 + 810);
    const auto *qsf0_811 = buffer.data(qsf0 + 811);
    const auto *qsf0_812 = buffer.data(qsf0 + 812);
    const auto *qsf0_813 = buffer.data(qsf0 + 813);
    const auto *qsf0_814 = buffer.data(qsf0 + 814);
    const auto *qsf0_815 = buffer.data(qsf0 + 815);
    const auto *qsf0_816 = buffer.data(qsf0 + 816);
    const auto *qsf0_817 = buffer.data(qsf0 + 817);
    const auto *qsf0_818 = buffer.data(qsf0 + 818);
    const auto *qsf0_819 = buffer.data(qsf0 + 819);
    const auto *qsf0_820 = buffer.data(qsf0 + 820);
    const auto *qsf0_821 = buffer.data(qsf0 + 821);
    const auto *qsf0_822 = buffer.data(qsf0 + 822);
    const auto *qsf0_823 = buffer.data(qsf0 + 823);
    const auto *qsf0_824 = buffer.data(qsf0 + 824);
    const auto *qsf0_825 = buffer.data(qsf0 + 825);
    const auto *qsf0_826 = buffer.data(qsf0 + 826);
    const auto *qsf0_827 = buffer.data(qsf0 + 827);
    const auto *qsf0_828 = buffer.data(qsf0 + 828);
    const auto *qsf0_829 = buffer.data(qsf0 + 829);
    const auto *qsf0_830 = buffer.data(qsf0 + 830);
    const auto *qsf0_831 = buffer.data(qsf0 + 831);
    const auto *qsf0_832 = buffer.data(qsf0 + 832);
    const auto *qsf0_833 = buffer.data(qsf0 + 833);
    const auto *qsf0_834 = buffer.data(qsf0 + 834);
    const auto *qsf0_835 = buffer.data(qsf0 + 835);
    const auto *qsf0_836 = buffer.data(qsf0 + 836);
    const auto *qsf0_837 = buffer.data(qsf0 + 837);
    const auto *qsf0_838 = buffer.data(qsf0 + 838);
    const auto *qsf0_839 = buffer.data(qsf0 + 839);

    const auto *qsf1_783 = buffer.data(qsf1 + 783);
    const auto *qsf1_785 = buffer.data(qsf1 + 785);
    const auto *qsf1_786 = buffer.data(qsf1 + 786);
    const auto *qsf1_787 = buffer.data(qsf1 + 787);
    const auto *qsf1_788 = buffer.data(qsf1 + 788);
    const auto *qsf1_789 = buffer.data(qsf1 + 789);
    const auto *qsf1_792 = buffer.data(qsf1 + 792);
    const auto *qsf1_794 = buffer.data(qsf1 + 794);
    const auto *qsf1_795 = buffer.data(qsf1 + 795);
    const auto *qsf1_797 = buffer.data(qsf1 + 797);
    const auto *qsf1_798 = buffer.data(qsf1 + 798);
    const auto *qsf1_799 = buffer.data(qsf1 + 799);
    const auto *qsf1_800 = buffer.data(qsf1 + 800);
    const auto *qsf1_801 = buffer.data(qsf1 + 801);
    const auto *qsf1_802 = buffer.data(qsf1 + 802);
    const auto *qsf1_803 = buffer.data(qsf1 + 803);
    const auto *qsf1_804 = buffer.data(qsf1 + 804);
    const auto *qsf1_805 = buffer.data(qsf1 + 805);
    const auto *qsf1_806 = buffer.data(qsf1 + 806);
    const auto *qsf1_807 = buffer.data(qsf1 + 807);
    const auto *qsf1_808 = buffer.data(qsf1 + 808);
    const auto *qsf1_809 = buffer.data(qsf1 + 809);
    const auto *qsf1_810 = buffer.data(qsf1 + 810);
    const auto *qsf1_811 = buffer.data(qsf1 + 811);
    const auto *qsf1_812 = buffer.data(qsf1 + 812);
    const auto *qsf1_813 = buffer.data(qsf1 + 813);
    const auto *qsf1_814 = buffer.data(qsf1 + 814);
    const auto *qsf1_815 = buffer.data(qsf1 + 815);
    const auto *qsf1_816 = buffer.data(qsf1 + 816);
    const auto *qsf1_817 = buffer.data(qsf1 + 817);
    const auto *qsf1_818 = buffer.data(qsf1 + 818);
    const auto *qsf1_819 = buffer.data(qsf1 + 819);
    const auto *qsf1_820 = buffer.data(qsf1 + 820);
    const auto *qsf1_821 = buffer.data(qsf1 + 821);
    const auto *qsf1_822 = buffer.data(qsf1 + 822);
    const auto *qsf1_823 = buffer.data(qsf1 + 823);
    const auto *qsf1_824 = buffer.data(qsf1 + 824);
    const auto *qsf1_825 = buffer.data(qsf1 + 825);
    const auto *qsf1_826 = buffer.data(qsf1 + 826);
    const auto *qsf1_827 = buffer.data(qsf1 + 827);
    const auto *qsf1_828 = buffer.data(qsf1 + 828);
    const auto *qsf1_829 = buffer.data(qsf1 + 829);
    const auto *qsf1_830 = buffer.data(qsf1 + 830);
    const auto *qsf1_831 = buffer.data(qsf1 + 831);
    const auto *qsf1_832 = buffer.data(qsf1 + 832);
    const auto *qsf1_833 = buffer.data(qsf1 + 833);
    const auto *qsf1_834 = buffer.data(qsf1 + 834);
    const auto *qsf1_835 = buffer.data(qsf1 + 835);
    const auto *qsf1_836 = buffer.data(qsf1 + 836);
    const auto *qsf1_837 = buffer.data(qsf1 + 837);
    const auto *qsf1_838 = buffer.data(qsf1 + 838);
    const auto *qsf1_839 = buffer.data(qsf1 + 839);

    const auto *qsg_1170 = buffer.data(qsg + 1170);
    const auto *qsg_1171 = buffer.data(qsg + 1171);
    const auto *qsg_1173 = buffer.data(qsg + 1173);
    const auto *qsg_1175 = buffer.data(qsg + 1175);
    const auto *qsg_1176 = buffer.data(qsg + 1176);
    const auto *qsg_1178 = buffer.data(qsg + 1178);
    const auto *qsg_1179 = buffer.data(qsg + 1179);
    const auto *qsg_1180 = buffer.data(qsg + 1180);
    const auto *qsg_1181 = buffer.data(qsg + 1181);
    const auto *qsg_1182 = buffer.data(qsg + 1182);
    const auto *qsg_1183 = buffer.data(qsg + 1183);
    const auto *qsg_1184 = buffer.data(qsg + 1184);
    const auto *qsg_1187 = buffer.data(qsg + 1187);
    const auto *qsg_1189 = buffer.data(qsg + 1189);
    const auto *qsg_1190 = buffer.data(qsg + 1190);
    const auto *qsg_1192 = buffer.data(qsg + 1192);
    const auto *qsg_1193 = buffer.data(qsg + 1193);
    const auto *qsg_1194 = buffer.data(qsg + 1194);
    const auto *qsg_1195 = buffer.data(qsg + 1195);
    const auto *qsg_1196 = buffer.data(qsg + 1196);
    const auto *qsg_1197 = buffer.data(qsg + 1197);
    const auto *qsg_1198 = buffer.data(qsg + 1198);
    const auto *qsg_1199 = buffer.data(qsg + 1199);
    const auto *qsg_1200 = buffer.data(qsg + 1200);
    const auto *qsg_1201 = buffer.data(qsg + 1201);
    const auto *qsg_1202 = buffer.data(qsg + 1202);
    const auto *qsg_1203 = buffer.data(qsg + 1203);
    const auto *qsg_1204 = buffer.data(qsg + 1204);
    const auto *qsg_1205 = buffer.data(qsg + 1205);
    const auto *qsg_1206 = buffer.data(qsg + 1206);
    const auto *qsg_1207 = buffer.data(qsg + 1207);
    const auto *qsg_1208 = buffer.data(qsg + 1208);
    const auto *qsg_1209 = buffer.data(qsg + 1209);
    const auto *qsg_1210 = buffer.data(qsg + 1210);
    const auto *qsg_1211 = buffer.data(qsg + 1211);
    const auto *qsg_1212 = buffer.data(qsg + 1212);
    const auto *qsg_1213 = buffer.data(qsg + 1213);
    const auto *qsg_1214 = buffer.data(qsg + 1214);
    const auto *qsg_1215 = buffer.data(qsg + 1215);
    const auto *qsg_1216 = buffer.data(qsg + 1216);
    const auto *qsg_1217 = buffer.data(qsg + 1217);
    const auto *qsg_1218 = buffer.data(qsg + 1218);
    const auto *qsg_1219 = buffer.data(qsg + 1219);
    const auto *qsg_1220 = buffer.data(qsg + 1220);
    const auto *qsg_1221 = buffer.data(qsg + 1221);
    const auto *qsg_1222 = buffer.data(qsg + 1222);
    const auto *qsg_1223 = buffer.data(qsg + 1223);
    const auto *qsg_1224 = buffer.data(qsg + 1224);
    const auto *qsg_1225 = buffer.data(qsg + 1225);
    const auto *qsg_1226 = buffer.data(qsg + 1226);
    const auto *qsg_1227 = buffer.data(qsg + 1227);
    const auto *qsg_1228 = buffer.data(qsg + 1228);
    const auto *qsg_1229 = buffer.data(qsg + 1229);
    const auto *qsg_1230 = buffer.data(qsg + 1230);
    const auto *qsg_1231 = buffer.data(qsg + 1231);
    const auto *qsg_1232 = buffer.data(qsg + 1232);
    const auto *qsg_1233 = buffer.data(qsg + 1233);
    const auto *qsg_1234 = buffer.data(qsg + 1234);
    const auto *qsg_1235 = buffer.data(qsg + 1235);
    const auto *qsg_1236 = buffer.data(qsg + 1236);
    const auto *qsg_1237 = buffer.data(qsg + 1237);
    const auto *qsg_1238 = buffer.data(qsg + 1238);
    const auto *qsg_1239 = buffer.data(qsg + 1239);
    const auto *qsg_1240 = buffer.data(qsg + 1240);
    const auto *qsg_1241 = buffer.data(qsg + 1241);
    const auto *qsg_1242 = buffer.data(qsg + 1242);
    const auto *qsg_1243 = buffer.data(qsg + 1243);
    const auto *qsg_1244 = buffer.data(qsg + 1244);
    const auto *qsg_1245 = buffer.data(qsg + 1245);
    const auto *qsg_1246 = buffer.data(qsg + 1246);
    const auto *qsg_1247 = buffer.data(qsg + 1247);
    const auto *qsg_1248 = buffer.data(qsg + 1248);
    const auto *qsg_1249 = buffer.data(qsg + 1249);
    const auto *qsg_1250 = buffer.data(qsg + 1250);
    const auto *qsg_1251 = buffer.data(qsg + 1251);
    const auto *qsg_1252 = buffer.data(qsg + 1252);
    const auto *qsg_1253 = buffer.data(qsg + 1253);
    const auto *qsg_1254 = buffer.data(qsg + 1254);
    const auto *qsg_1255 = buffer.data(qsg + 1255);
    const auto *qsg_1256 = buffer.data(qsg + 1256);
    const auto *qsg_1257 = buffer.data(qsg + 1257);
    const auto *qsg_1258 = buffer.data(qsg + 1258);
    const auto *qsg_1259 = buffer.data(qsg + 1259);

#pragma omp simd aligned(t_1640, t_1641, t_1642, t_1643, pc_x, pc_z, qsf0_783, qsf0_785, \
                         qsf1_783, qsf1_785, qsg_1170, qsg_1171, qsg_1173, \
                         qsg_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1640[k] = f_3 * pc_z[k] * qsg_1170[k];

        t_1641[k] = f_6 * qsf0_783[k]
                    - f_7 * qsf1_783[k]
                    + f_3 * pc_x[k] * qsg_1173[k];

        t_1642[k] = f_3 * pc_z[k] * qsg_1171[k];

        t_1643[k] = f_6 * qsf0_785[k]
                    - f_7 * qsf1_785[k]
                    + f_3 * pc_x[k] * qsg_1175[k];
    }

#pragma omp simd aligned(t_1644, t_1645, t_1646, t_1647, pc_x, pc_z, qsf0_786, qsf0_788, \
                         qsf0_789, qsf1_786, qsf1_788, qsf1_789, qsg_1173, qsg_1176, qsg_1178, \
                         qsg_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1644[k] = f_4 * qsf0_786[k]
                    - f_5 * qsf1_786[k]
                    + f_3 * pc_x[k] * qsg_1176[k];

        t_1645[k] = f_3 * pc_z[k] * qsg_1173[k];

        t_1646[k] = f_4 * qsf0_788[k]
                    - f_5 * qsf1_788[k]
                    + f_3 * pc_x[k] * qsg_1178[k];

        t_1647[k] = f_4 * qsf0_789[k]
                    - f_5 * qsf1_789[k]
                    + f_3 * pc_x[k] * qsg_1179[k];
    }

#pragma omp simd aligned(t_1648, t_1649, t_1650, t_1651, t_1652, t_1653, pc_x, pc_y, osg_1000, \
                         qsf0_786, qsf1_786, qsg_1180, qsg_1181, qsg_1182, qsg_1183, \
                         qsg_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1648[k] = f_3 * pc_x[k] * qsg_1180[k];

        t_1649[k] = f_3 * pc_x[k] * qsg_1181[k];

        t_1650[k] = f_3 * pc_x[k] * qsg_1182[k];

        t_1651[k] = f_3 * pc_x[k] * qsg_1183[k];

        t_1652[k] = f_3 * pc_x[k] * qsg_1184[k];

        t_1653[k] = f_0 * osg_1000[k]
                    + f_1 * qsf0_786[k]
                    - f_2 * qsf1_786[k]
                    + f_3 * pc_y[k] * qsg_1180[k];
    }

#pragma omp simd aligned(t_1654, t_1655, t_1656, t_1657, pc_y, pc_z, osg_1004, qsf0_786, \
                         qsf0_787, qsf1_786, qsf1_787, qsg_1180, qsg_1181, qsg_1182, \
                         qsg_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1654[k] = f_3 * pc_z[k] * qsg_1180[k];

        t_1655[k] = f_4 * qsf0_786[k]
                    - f_5 * qsf1_786[k]
                    + f_3 * pc_z[k] * qsg_1181[k];

        t_1656[k] = f_6 * qsf0_787[k]
                    - f_7 * qsf1_787[k]
                    + f_3 * pc_z[k] * qsg_1182[k];

        t_1657[k] = f_0 * osg_1004[k]
                    + f_3 * pc_y[k] * qsg_1184[k];
    }

#pragma omp simd aligned(t_1658, t_1659, t_1660, pa_z, pc_z, osh0_1386, osh0_1387, osh1_1386, \
                         osh1_1387, qsf0_789, qsf1_789, qsg_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1658[k] = f_1 * qsf0_789[k]
                    - f_2 * qsf1_789[k]
                    + f_3 * pc_z[k] * qsg_1184[k];

        t_1659[k] = pa_z[k] * osh0_1386[k]
                    - f_8 * pc_z[k] * osh1_1386[k];

        t_1660[k] = pa_z[k] * osh0_1387[k]
                    - f_8 * pc_z[k] * osh1_1387[k];
    }

#pragma omp simd aligned(t_1661, t_1662, t_1663, pa_z, pc_x, pc_z, osh0_1389, osh1_1389, \
                         qsf0_792, qsf0_794, qsf1_792, qsf1_794, qsg_1187, \
                         qsg_1189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1661[k] = f_13 * qsf0_792[k]
                    - f_14 * qsf1_792[k]
                    + f_3 * pc_x[k] * qsg_1187[k];

        t_1662[k] = pa_z[k] * osh0_1389[k]
                    - f_8 * pc_z[k] * osh1_1389[k];

        t_1663[k] = f_6 * qsf0_794[k]
                    - f_7 * qsf1_794[k]
                    + f_3 * pc_x[k] * qsg_1189[k];
    }

#pragma omp simd aligned(t_1664, t_1665, t_1666, pa_z, pc_x, pc_z, osh0_1392, osh1_1392, \
                         qsf0_795, qsf0_797, qsf1_795, qsf1_797, qsg_1190, \
                         qsg_1192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1664[k] = f_6 * qsf0_795[k]
                    - f_7 * qsf1_795[k]
                    + f_3 * pc_x[k] * qsg_1190[k];

        t_1665[k] = pa_z[k] * osh0_1392[k]
                    - f_8 * pc_z[k] * osh1_1392[k];

        t_1666[k] = f_4 * qsf0_797[k]
                    - f_5 * qsf1_797[k]
                    + f_3 * pc_x[k] * qsg_1192[k];
    }

#pragma omp simd aligned(t_1667, t_1668, t_1669, t_1670, t_1671, pc_x, qsf0_798, qsf0_799, \
                         qsf1_798, qsf1_799, qsg_1193, qsg_1194, qsg_1195, qsg_1196, \
                         qsg_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1667[k] = f_4 * qsf0_798[k]
                    - f_5 * qsf1_798[k]
                    + f_3 * pc_x[k] * qsg_1193[k];

        t_1668[k] = f_4 * qsf0_799[k]
                    - f_5 * qsf1_799[k]
                    + f_3 * pc_x[k] * qsg_1194[k];

        t_1669[k] = f_3 * pc_x[k] * qsg_1195[k];

        t_1670[k] = f_3 * pc_x[k] * qsg_1196[k];

        t_1671[k] = f_3 * pc_x[k] * qsg_1197[k];
    }

#pragma omp simd aligned(t_1672, t_1673, t_1674, t_1675, pa_z, pc_x, pc_z, osh0_1401, \
                         osg_1000, osh1_1401, qsg_1195, qsg_1198, \
                         qsg_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1672[k] = f_3 * pc_x[k] * qsg_1198[k];

        t_1673[k] = f_3 * pc_x[k] * qsg_1199[k];

        t_1674[k] = pa_z[k] * osh0_1401[k]
                    - f_8 * pc_z[k] * osh1_1401[k];

        t_1675[k] = f_9 * osg_1000[k]
                    + f_3 * pc_z[k] * qsg_1195[k];
    }

#pragma omp simd aligned(t_1676, t_1677, t_1678, pa_z, pc_y, pc_z, osh0_1403, osh0_1404, \
                         osg_1001, osg_1002, osg_1019, osh1_1403, osh1_1404, \
                         qsg_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1676[k] = pa_z[k] * osh0_1403[k]
                    + f_10 * osg_1001[k]
                    - f_8 * pc_z[k] * osh1_1403[k];

        t_1677[k] = pa_z[k] * osh0_1404[k]
                    + f_11 * osg_1002[k]
                    - f_8 * pc_z[k] * osh1_1404[k];

        t_1678[k] = f_12 * osg_1019[k]
                    + f_3 * pc_y[k] * qsg_1199[k];
    }

#pragma omp simd aligned(t_1679, t_1680, t_1681, pc_x, pc_z, osg_1004, qsf0_799, qsf0_800, \
                         qsf0_801, qsf1_799, qsf1_800, qsf1_801, qsg_1199, qsg_1200, \
                         qsg_1201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1679[k] = f_9 * osg_1004[k]
                    + f_1 * qsf0_799[k]
                    - f_2 * qsf1_799[k]
                    + f_3 * pc_z[k] * qsg_1199[k];

        t_1680[k] = f_1 * qsf0_800[k]
                    - f_2 * qsf1_800[k]
                    + f_3 * pc_x[k] * qsg_1200[k];

        t_1681[k] = f_13 * qsf0_801[k]
                    - f_14 * qsf1_801[k]
                    + f_3 * pc_x[k] * qsg_1201[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, pc_x, qsf0_802, qsf0_803, qsf0_804, qsf1_802, \
                         qsf1_803, qsf1_804, qsg_1202, qsg_1203, \
                         qsg_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_13 * qsf0_802[k]
                    - f_14 * qsf1_802[k]
                    + f_3 * pc_x[k] * qsg_1202[k];

        t_1683[k] = f_6 * qsf0_803[k]
                    - f_7 * qsf1_803[k]
                    + f_3 * pc_x[k] * qsg_1203[k];

        t_1684[k] = f_6 * qsf0_804[k]
                    - f_7 * qsf1_804[k]
                    + f_3 * pc_x[k] * qsg_1204[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pc_x, qsf0_805, qsf0_806, qsf0_807, qsf1_805, \
                         qsf1_806, qsf1_807, qsg_1205, qsg_1206, \
                         qsg_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = f_6 * qsf0_805[k]
                    - f_7 * qsf1_805[k]
                    + f_3 * pc_x[k] * qsg_1205[k];

        t_1686[k] = f_4 * qsf0_806[k]
                    - f_5 * qsf1_806[k]
                    + f_3 * pc_x[k] * qsg_1206[k];

        t_1687[k] = f_4 * qsf0_807[k]
                    - f_5 * qsf1_807[k]
                    + f_3 * pc_x[k] * qsg_1207[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, t_1691, t_1692, pc_x, qsf0_808, qsf0_809, \
                         qsf1_808, qsf1_809, qsg_1208, qsg_1209, qsg_1210, qsg_1211, \
                         qsg_1212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = f_4 * qsf0_808[k]
                    - f_5 * qsf1_808[k]
                    + f_3 * pc_x[k] * qsg_1208[k];

        t_1689[k] = f_4 * qsf0_809[k]
                    - f_5 * qsf1_809[k]
                    + f_3 * pc_x[k] * qsg_1209[k];

        t_1690[k] = f_3 * pc_x[k] * qsg_1210[k];

        t_1691[k] = f_3 * pc_x[k] * qsg_1211[k];

        t_1692[k] = f_3 * pc_x[k] * qsg_1212[k];
    }

#pragma omp simd aligned(t_1693, t_1694, t_1695, t_1696, pc_x, pc_y, pc_z, osg_1015, osg_1030, \
                         qsf0_806, qsf1_806, qsg_1210, qsg_1213, \
                         qsg_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1693[k] = f_3 * pc_x[k] * qsg_1213[k];

        t_1694[k] = f_3 * pc_x[k] * qsg_1214[k];

        t_1695[k] = f_15 * osg_1030[k]
                    + f_1 * qsf0_806[k]
                    - f_2 * qsf1_806[k]
                    + f_3 * pc_y[k] * qsg_1210[k];

        t_1696[k] = f_10 * osg_1015[k]
                    + f_3 * pc_z[k] * qsg_1210[k];
    }

#pragma omp simd aligned(t_1697, t_1698, t_1699, pc_y, osg_1032, osg_1033, osg_1034, qsf0_808, \
                         qsf0_809, qsf1_808, qsf1_809, qsg_1212, qsg_1213, \
                         qsg_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1697[k] = f_15 * osg_1032[k]
                    + f_6 * qsf0_808[k]
                    - f_7 * qsf1_808[k]
                    + f_3 * pc_y[k] * qsg_1212[k];

        t_1698[k] = f_15 * osg_1033[k]
                    + f_4 * qsf0_809[k]
                    - f_5 * qsf1_809[k]
                    + f_3 * pc_y[k] * qsg_1213[k];

        t_1699[k] = f_15 * osg_1034[k]
                    + f_3 * pc_y[k] * qsg_1214[k];
    }

#pragma omp simd aligned(t_1700, t_1701, t_1702, pc_x, pc_z, osg_1019, qsf0_809, qsf0_810, \
                         qsf0_811, qsf1_809, qsf1_810, qsf1_811, qsg_1214, qsg_1215, \
                         qsg_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1700[k] = f_10 * osg_1019[k]
                    + f_1 * qsf0_809[k]
                    - f_2 * qsf1_809[k]
                    + f_3 * pc_z[k] * qsg_1214[k];

        t_1701[k] = f_1 * qsf0_810[k]
                    - f_2 * qsf1_810[k]
                    + f_3 * pc_x[k] * qsg_1215[k];

        t_1702[k] = f_13 * qsf0_811[k]
                    - f_14 * qsf1_811[k]
                    + f_3 * pc_x[k] * qsg_1216[k];
    }

#pragma omp simd aligned(t_1703, t_1704, t_1705, pc_x, qsf0_812, qsf0_813, qsf0_814, qsf1_812, \
                         qsf1_813, qsf1_814, qsg_1217, qsg_1218, \
                         qsg_1219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1703[k] = f_13 * qsf0_812[k]
                    - f_14 * qsf1_812[k]
                    + f_3 * pc_x[k] * qsg_1217[k];

        t_1704[k] = f_6 * qsf0_813[k]
                    - f_7 * qsf1_813[k]
                    + f_3 * pc_x[k] * qsg_1218[k];

        t_1705[k] = f_6 * qsf0_814[k]
                    - f_7 * qsf1_814[k]
                    + f_3 * pc_x[k] * qsg_1219[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, pc_x, qsf0_815, qsf0_816, qsf0_817, qsf1_815, \
                         qsf1_816, qsf1_817, qsg_1220, qsg_1221, \
                         qsg_1222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = f_6 * qsf0_815[k]
                    - f_7 * qsf1_815[k]
                    + f_3 * pc_x[k] * qsg_1220[k];

        t_1707[k] = f_4 * qsf0_816[k]
                    - f_5 * qsf1_816[k]
                    + f_3 * pc_x[k] * qsg_1221[k];

        t_1708[k] = f_4 * qsf0_817[k]
                    - f_5 * qsf1_817[k]
                    + f_3 * pc_x[k] * qsg_1222[k];
    }

#pragma omp simd aligned(t_1709, t_1710, t_1711, t_1712, t_1713, pc_x, qsf0_818, qsf0_819, \
                         qsf1_818, qsf1_819, qsg_1223, qsg_1224, qsg_1225, qsg_1226, \
                         qsg_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1709[k] = f_4 * qsf0_818[k]
                    - f_5 * qsf1_818[k]
                    + f_3 * pc_x[k] * qsg_1223[k];

        t_1710[k] = f_4 * qsf0_819[k]
                    - f_5 * qsf1_819[k]
                    + f_3 * pc_x[k] * qsg_1224[k];

        t_1711[k] = f_3 * pc_x[k] * qsg_1225[k];

        t_1712[k] = f_3 * pc_x[k] * qsg_1226[k];

        t_1713[k] = f_3 * pc_x[k] * qsg_1227[k];
    }

#pragma omp simd aligned(t_1714, t_1715, t_1716, t_1717, pc_x, pc_y, pc_z, osg_1030, osg_1045, \
                         qsf0_816, qsf1_816, qsg_1225, qsg_1228, \
                         qsg_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = f_3 * pc_x[k] * qsg_1228[k];

        t_1715[k] = f_3 * pc_x[k] * qsg_1229[k];

        t_1716[k] = f_16 * osg_1045[k]
                    + f_1 * qsf0_816[k]
                    - f_2 * qsf1_816[k]
                    + f_3 * pc_y[k] * qsg_1225[k];

        t_1717[k] = f_11 * osg_1030[k]
                    + f_3 * pc_z[k] * qsg_1225[k];
    }

#pragma omp simd aligned(t_1718, t_1719, t_1720, pc_y, osg_1047, osg_1048, osg_1049, qsf0_818, \
                         qsf0_819, qsf1_818, qsf1_819, qsg_1227, qsg_1228, \
                         qsg_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1718[k] = f_16 * osg_1047[k]
                    + f_6 * qsf0_818[k]
                    - f_7 * qsf1_818[k]
                    + f_3 * pc_y[k] * qsg_1227[k];

        t_1719[k] = f_16 * osg_1048[k]
                    + f_4 * qsf0_819[k]
                    - f_5 * qsf1_819[k]
                    + f_3 * pc_y[k] * qsg_1228[k];

        t_1720[k] = f_16 * osg_1049[k]
                    + f_3 * pc_y[k] * qsg_1229[k];
    }

#pragma omp simd aligned(t_1721, t_1722, t_1723, pc_x, pc_z, osg_1034, qsf0_819, qsf0_820, \
                         qsf0_821, qsf1_819, qsf1_820, qsf1_821, qsg_1229, qsg_1230, \
                         qsg_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1721[k] = f_11 * osg_1034[k]
                    + f_1 * qsf0_819[k]
                    - f_2 * qsf1_819[k]
                    + f_3 * pc_z[k] * qsg_1229[k];

        t_1722[k] = f_1 * qsf0_820[k]
                    - f_2 * qsf1_820[k]
                    + f_3 * pc_x[k] * qsg_1230[k];

        t_1723[k] = f_13 * qsf0_821[k]
                    - f_14 * qsf1_821[k]
                    + f_3 * pc_x[k] * qsg_1231[k];
    }

#pragma omp simd aligned(t_1724, t_1725, t_1726, pc_x, qsf0_822, qsf0_823, qsf0_824, qsf1_822, \
                         qsf1_823, qsf1_824, qsg_1232, qsg_1233, \
                         qsg_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1724[k] = f_13 * qsf0_822[k]
                    - f_14 * qsf1_822[k]
                    + f_3 * pc_x[k] * qsg_1232[k];

        t_1725[k] = f_6 * qsf0_823[k]
                    - f_7 * qsf1_823[k]
                    + f_3 * pc_x[k] * qsg_1233[k];

        t_1726[k] = f_6 * qsf0_824[k]
                    - f_7 * qsf1_824[k]
                    + f_3 * pc_x[k] * qsg_1234[k];
    }

#pragma omp simd aligned(t_1727, t_1728, t_1729, pc_x, qsf0_825, qsf0_826, qsf0_827, qsf1_825, \
                         qsf1_826, qsf1_827, qsg_1235, qsg_1236, \
                         qsg_1237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1727[k] = f_6 * qsf0_825[k]
                    - f_7 * qsf1_825[k]
                    + f_3 * pc_x[k] * qsg_1235[k];

        t_1728[k] = f_4 * qsf0_826[k]
                    - f_5 * qsf1_826[k]
                    + f_3 * pc_x[k] * qsg_1236[k];

        t_1729[k] = f_4 * qsf0_827[k]
                    - f_5 * qsf1_827[k]
                    + f_3 * pc_x[k] * qsg_1237[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, t_1733, t_1734, pc_x, qsf0_828, qsf0_829, \
                         qsf1_828, qsf1_829, qsg_1238, qsg_1239, qsg_1240, qsg_1241, \
                         qsg_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = f_4 * qsf0_828[k]
                    - f_5 * qsf1_828[k]
                    + f_3 * pc_x[k] * qsg_1238[k];

        t_1731[k] = f_4 * qsf0_829[k]
                    - f_5 * qsf1_829[k]
                    + f_3 * pc_x[k] * qsg_1239[k];

        t_1732[k] = f_3 * pc_x[k] * qsg_1240[k];

        t_1733[k] = f_3 * pc_x[k] * qsg_1241[k];

        t_1734[k] = f_3 * pc_x[k] * qsg_1242[k];
    }

#pragma omp simd aligned(t_1735, t_1736, t_1737, t_1738, pc_x, pc_y, pc_z, osg_1045, osg_1060, \
                         qsf0_826, qsf1_826, qsg_1240, qsg_1243, \
                         qsg_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1735[k] = f_3 * pc_x[k] * qsg_1243[k];

        t_1736[k] = f_3 * pc_x[k] * qsg_1244[k];

        t_1737[k] = f_17 * osg_1060[k]
                    + f_1 * qsf0_826[k]
                    - f_2 * qsf1_826[k]
                    + f_3 * pc_y[k] * qsg_1240[k];

        t_1738[k] = f_18 * osg_1045[k]
                    + f_3 * pc_z[k] * qsg_1240[k];
    }

#pragma omp simd aligned(t_1739, t_1740, t_1741, pc_y, osg_1062, osg_1063, osg_1064, qsf0_828, \
                         qsf0_829, qsf1_828, qsf1_829, qsg_1242, qsg_1243, \
                         qsg_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1739[k] = f_17 * osg_1062[k]
                    + f_6 * qsf0_828[k]
                    - f_7 * qsf1_828[k]
                    + f_3 * pc_y[k] * qsg_1242[k];

        t_1740[k] = f_17 * osg_1063[k]
                    + f_4 * qsf0_829[k]
                    - f_5 * qsf1_829[k]
                    + f_3 * pc_y[k] * qsg_1243[k];

        t_1741[k] = f_17 * osg_1064[k]
                    + f_3 * pc_y[k] * qsg_1244[k];
    }

#pragma omp simd aligned(t_1742, t_1743, t_1744, pc_x, pc_z, osg_1049, qsf0_829, qsf0_830, \
                         qsf0_831, qsf1_829, qsf1_830, qsf1_831, qsg_1244, qsg_1245, \
                         qsg_1246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1742[k] = f_18 * osg_1049[k]
                    + f_1 * qsf0_829[k]
                    - f_2 * qsf1_829[k]
                    + f_3 * pc_z[k] * qsg_1244[k];

        t_1743[k] = f_1 * qsf0_830[k]
                    - f_2 * qsf1_830[k]
                    + f_3 * pc_x[k] * qsg_1245[k];

        t_1744[k] = f_13 * qsf0_831[k]
                    - f_14 * qsf1_831[k]
                    + f_3 * pc_x[k] * qsg_1246[k];
    }

#pragma omp simd aligned(t_1745, t_1746, t_1747, pc_x, qsf0_832, qsf0_833, qsf0_834, qsf1_832, \
                         qsf1_833, qsf1_834, qsg_1247, qsg_1248, \
                         qsg_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1745[k] = f_13 * qsf0_832[k]
                    - f_14 * qsf1_832[k]
                    + f_3 * pc_x[k] * qsg_1247[k];

        t_1746[k] = f_6 * qsf0_833[k]
                    - f_7 * qsf1_833[k]
                    + f_3 * pc_x[k] * qsg_1248[k];

        t_1747[k] = f_6 * qsf0_834[k]
                    - f_7 * qsf1_834[k]
                    + f_3 * pc_x[k] * qsg_1249[k];
    }

#pragma omp simd aligned(t_1748, t_1749, t_1750, pc_x, qsf0_835, qsf0_836, qsf0_837, qsf1_835, \
                         qsf1_836, qsf1_837, qsg_1250, qsg_1251, \
                         qsg_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1748[k] = f_6 * qsf0_835[k]
                    - f_7 * qsf1_835[k]
                    + f_3 * pc_x[k] * qsg_1250[k];

        t_1749[k] = f_4 * qsf0_836[k]
                    - f_5 * qsf1_836[k]
                    + f_3 * pc_x[k] * qsg_1251[k];

        t_1750[k] = f_4 * qsf0_837[k]
                    - f_5 * qsf1_837[k]
                    + f_3 * pc_x[k] * qsg_1252[k];
    }

#pragma omp simd aligned(t_1751, t_1752, t_1753, t_1754, t_1755, pc_x, qsf0_838, qsf0_839, \
                         qsf1_838, qsf1_839, qsg_1253, qsg_1254, qsg_1255, qsg_1256, \
                         qsg_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1751[k] = f_4 * qsf0_838[k]
                    - f_5 * qsf1_838[k]
                    + f_3 * pc_x[k] * qsg_1253[k];

        t_1752[k] = f_4 * qsf0_839[k]
                    - f_5 * qsf1_839[k]
                    + f_3 * pc_x[k] * qsg_1254[k];

        t_1753[k] = f_3 * pc_x[k] * qsg_1255[k];

        t_1754[k] = f_3 * pc_x[k] * qsg_1256[k];

        t_1755[k] = f_3 * pc_x[k] * qsg_1257[k];
    }

#pragma omp simd aligned(t_1756, t_1757, t_1758, t_1759, pc_x, pc_y, pc_z, osg_1060, osg_1075, \
                         qsf0_836, qsf1_836, qsg_1255, qsg_1258, \
                         qsg_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1756[k] = f_3 * pc_x[k] * qsg_1258[k];

        t_1757[k] = f_3 * pc_x[k] * qsg_1259[k];

        t_1758[k] = f_19 * osg_1075[k]
                    + f_1 * qsf0_836[k]
                    - f_2 * qsf1_836[k]
                    + f_3 * pc_y[k] * qsg_1255[k];

        t_1759[k] = f_20 * osg_1060[k]
                    + f_3 * pc_z[k] * qsg_1255[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsf0, const size_t qsf1,
                                                           const size_t qsg, const size_t ncols,
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
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 2.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_1617 = buffer.data(osh0 + 1617);
    const auto *osh0_1619 = buffer.data(osh0 + 1619);

    const auto *osg_1064 = buffer.data(osg + 1064);
    const auto *osg_1075 = buffer.data(osg + 1075);
    const auto *osg_1077 = buffer.data(osg + 1077);
    const auto *osg_1078 = buffer.data(osg + 1078);
    const auto *osg_1079 = buffer.data(osg + 1079);
    const auto *osg_1090 = buffer.data(osg + 1090);
    const auto *osg_1092 = buffer.data(osg + 1092);
    const auto *osg_1093 = buffer.data(osg + 1093);
    const auto *osg_1094 = buffer.data(osg + 1094);
    const auto *osg_1105 = buffer.data(osg + 1105);
    const auto *osg_1107 = buffer.data(osg + 1107);
    const auto *osg_1108 = buffer.data(osg + 1108);
    const auto *osg_1109 = buffer.data(osg + 1109);
    const auto *osg_1120 = buffer.data(osg + 1120);
    const auto *osg_1122 = buffer.data(osg + 1122);
    const auto *osg_1123 = buffer.data(osg + 1123);
    const auto *osg_1124 = buffer.data(osg + 1124);
    const auto *osg_1135 = buffer.data(osg + 1135);
    const auto *osg_1137 = buffer.data(osg + 1137);
    const auto *osg_1138 = buffer.data(osg + 1138);
    const auto *osg_1139 = buffer.data(osg + 1139);
    const auto *osg_1150 = buffer.data(osg + 1150);
    const auto *osg_1152 = buffer.data(osg + 1152);
    const auto *osg_1153 = buffer.data(osg + 1153);
    const auto *osg_1154 = buffer.data(osg + 1154);

    const auto *osh1_1617 = buffer.data(osh1 + 1617);
    const auto *osh1_1619 = buffer.data(osh1 + 1619);

    const auto *qsf0_838 = buffer.data(qsf0 + 838);
    const auto *qsf0_839 = buffer.data(qsf0 + 839);
    const auto *qsf0_840 = buffer.data(qsf0 + 840);
    const auto *qsf0_841 = buffer.data(qsf0 + 841);
    const auto *qsf0_842 = buffer.data(qsf0 + 842);
    const auto *qsf0_843 = buffer.data(qsf0 + 843);
    const auto *qsf0_844 = buffer.data(qsf0 + 844);
    const auto *qsf0_845 = buffer.data(qsf0 + 845);
    const auto *qsf0_846 = buffer.data(qsf0 + 846);
    const auto *qsf0_847 = buffer.data(qsf0 + 847);
    const auto *qsf0_848 = buffer.data(qsf0 + 848);
    const auto *qsf0_849 = buffer.data(qsf0 + 849);
    const auto *qsf0_850 = buffer.data(qsf0 + 850);
    const auto *qsf0_851 = buffer.data(qsf0 + 851);
    const auto *qsf0_852 = buffer.data(qsf0 + 852);
    const auto *qsf0_853 = buffer.data(qsf0 + 853);
    const auto *qsf0_854 = buffer.data(qsf0 + 854);
    const auto *qsf0_855 = buffer.data(qsf0 + 855);
    const auto *qsf0_856 = buffer.data(qsf0 + 856);
    const auto *qsf0_857 = buffer.data(qsf0 + 857);
    const auto *qsf0_858 = buffer.data(qsf0 + 858);
    const auto *qsf0_859 = buffer.data(qsf0 + 859);
    const auto *qsf0_860 = buffer.data(qsf0 + 860);
    const auto *qsf0_861 = buffer.data(qsf0 + 861);
    const auto *qsf0_862 = buffer.data(qsf0 + 862);
    const auto *qsf0_863 = buffer.data(qsf0 + 863);
    const auto *qsf0_864 = buffer.data(qsf0 + 864);
    const auto *qsf0_865 = buffer.data(qsf0 + 865);
    const auto *qsf0_866 = buffer.data(qsf0 + 866);
    const auto *qsf0_867 = buffer.data(qsf0 + 867);
    const auto *qsf0_868 = buffer.data(qsf0 + 868);
    const auto *qsf0_869 = buffer.data(qsf0 + 869);
    const auto *qsf0_870 = buffer.data(qsf0 + 870);
    const auto *qsf0_871 = buffer.data(qsf0 + 871);
    const auto *qsf0_872 = buffer.data(qsf0 + 872);
    const auto *qsf0_873 = buffer.data(qsf0 + 873);
    const auto *qsf0_874 = buffer.data(qsf0 + 874);
    const auto *qsf0_875 = buffer.data(qsf0 + 875);
    const auto *qsf0_876 = buffer.data(qsf0 + 876);
    const auto *qsf0_877 = buffer.data(qsf0 + 877);
    const auto *qsf0_878 = buffer.data(qsf0 + 878);
    const auto *qsf0_879 = buffer.data(qsf0 + 879);
    const auto *qsf0_880 = buffer.data(qsf0 + 880);
    const auto *qsf0_881 = buffer.data(qsf0 + 881);
    const auto *qsf0_882 = buffer.data(qsf0 + 882);
    const auto *qsf0_883 = buffer.data(qsf0 + 883);
    const auto *qsf0_884 = buffer.data(qsf0 + 884);
    const auto *qsf0_885 = buffer.data(qsf0 + 885);
    const auto *qsf0_886 = buffer.data(qsf0 + 886);
    const auto *qsf0_887 = buffer.data(qsf0 + 887);
    const auto *qsf0_888 = buffer.data(qsf0 + 888);
    const auto *qsf0_889 = buffer.data(qsf0 + 889);
    const auto *qsf0_891 = buffer.data(qsf0 + 891);
    const auto *qsf0_893 = buffer.data(qsf0 + 893);
    const auto *qsf0_894 = buffer.data(qsf0 + 894);

    const auto *qsf1_838 = buffer.data(qsf1 + 838);
    const auto *qsf1_839 = buffer.data(qsf1 + 839);
    const auto *qsf1_840 = buffer.data(qsf1 + 840);
    const auto *qsf1_841 = buffer.data(qsf1 + 841);
    const auto *qsf1_842 = buffer.data(qsf1 + 842);
    const auto *qsf1_843 = buffer.data(qsf1 + 843);
    const auto *qsf1_844 = buffer.data(qsf1 + 844);
    const auto *qsf1_845 = buffer.data(qsf1 + 845);
    const auto *qsf1_846 = buffer.data(qsf1 + 846);
    const auto *qsf1_847 = buffer.data(qsf1 + 847);
    const auto *qsf1_848 = buffer.data(qsf1 + 848);
    const auto *qsf1_849 = buffer.data(qsf1 + 849);
    const auto *qsf1_850 = buffer.data(qsf1 + 850);
    const auto *qsf1_851 = buffer.data(qsf1 + 851);
    const auto *qsf1_852 = buffer.data(qsf1 + 852);
    const auto *qsf1_853 = buffer.data(qsf1 + 853);
    const auto *qsf1_854 = buffer.data(qsf1 + 854);
    const auto *qsf1_855 = buffer.data(qsf1 + 855);
    const auto *qsf1_856 = buffer.data(qsf1 + 856);
    const auto *qsf1_857 = buffer.data(qsf1 + 857);
    const auto *qsf1_858 = buffer.data(qsf1 + 858);
    const auto *qsf1_859 = buffer.data(qsf1 + 859);
    const auto *qsf1_860 = buffer.data(qsf1 + 860);
    const auto *qsf1_861 = buffer.data(qsf1 + 861);
    const auto *qsf1_862 = buffer.data(qsf1 + 862);
    const auto *qsf1_863 = buffer.data(qsf1 + 863);
    const auto *qsf1_864 = buffer.data(qsf1 + 864);
    const auto *qsf1_865 = buffer.data(qsf1 + 865);
    const auto *qsf1_866 = buffer.data(qsf1 + 866);
    const auto *qsf1_867 = buffer.data(qsf1 + 867);
    const auto *qsf1_868 = buffer.data(qsf1 + 868);
    const auto *qsf1_869 = buffer.data(qsf1 + 869);
    const auto *qsf1_870 = buffer.data(qsf1 + 870);
    const auto *qsf1_871 = buffer.data(qsf1 + 871);
    const auto *qsf1_872 = buffer.data(qsf1 + 872);
    const auto *qsf1_873 = buffer.data(qsf1 + 873);
    const auto *qsf1_874 = buffer.data(qsf1 + 874);
    const auto *qsf1_875 = buffer.data(qsf1 + 875);
    const auto *qsf1_876 = buffer.data(qsf1 + 876);
    const auto *qsf1_877 = buffer.data(qsf1 + 877);
    const auto *qsf1_878 = buffer.data(qsf1 + 878);
    const auto *qsf1_879 = buffer.data(qsf1 + 879);
    const auto *qsf1_880 = buffer.data(qsf1 + 880);
    const auto *qsf1_881 = buffer.data(qsf1 + 881);
    const auto *qsf1_882 = buffer.data(qsf1 + 882);
    const auto *qsf1_883 = buffer.data(qsf1 + 883);
    const auto *qsf1_884 = buffer.data(qsf1 + 884);
    const auto *qsf1_885 = buffer.data(qsf1 + 885);
    const auto *qsf1_886 = buffer.data(qsf1 + 886);
    const auto *qsf1_887 = buffer.data(qsf1 + 887);
    const auto *qsf1_888 = buffer.data(qsf1 + 888);
    const auto *qsf1_889 = buffer.data(qsf1 + 889);
    const auto *qsf1_891 = buffer.data(qsf1 + 891);
    const auto *qsf1_893 = buffer.data(qsf1 + 893);
    const auto *qsf1_894 = buffer.data(qsf1 + 894);

    const auto *qsg_1257 = buffer.data(qsg + 1257);
    const auto *qsg_1258 = buffer.data(qsg + 1258);
    const auto *qsg_1259 = buffer.data(qsg + 1259);
    const auto *qsg_1260 = buffer.data(qsg + 1260);
    const auto *qsg_1261 = buffer.data(qsg + 1261);
    const auto *qsg_1262 = buffer.data(qsg + 1262);
    const auto *qsg_1263 = buffer.data(qsg + 1263);
    const auto *qsg_1264 = buffer.data(qsg + 1264);
    const auto *qsg_1265 = buffer.data(qsg + 1265);
    const auto *qsg_1266 = buffer.data(qsg + 1266);
    const auto *qsg_1267 = buffer.data(qsg + 1267);
    const auto *qsg_1268 = buffer.data(qsg + 1268);
    const auto *qsg_1269 = buffer.data(qsg + 1269);
    const auto *qsg_1270 = buffer.data(qsg + 1270);
    const auto *qsg_1271 = buffer.data(qsg + 1271);
    const auto *qsg_1272 = buffer.data(qsg + 1272);
    const auto *qsg_1273 = buffer.data(qsg + 1273);
    const auto *qsg_1274 = buffer.data(qsg + 1274);
    const auto *qsg_1275 = buffer.data(qsg + 1275);
    const auto *qsg_1276 = buffer.data(qsg + 1276);
    const auto *qsg_1277 = buffer.data(qsg + 1277);
    const auto *qsg_1278 = buffer.data(qsg + 1278);
    const auto *qsg_1279 = buffer.data(qsg + 1279);
    const auto *qsg_1280 = buffer.data(qsg + 1280);
    const auto *qsg_1281 = buffer.data(qsg + 1281);
    const auto *qsg_1282 = buffer.data(qsg + 1282);
    const auto *qsg_1283 = buffer.data(qsg + 1283);
    const auto *qsg_1284 = buffer.data(qsg + 1284);
    const auto *qsg_1285 = buffer.data(qsg + 1285);
    const auto *qsg_1286 = buffer.data(qsg + 1286);
    const auto *qsg_1287 = buffer.data(qsg + 1287);
    const auto *qsg_1288 = buffer.data(qsg + 1288);
    const auto *qsg_1289 = buffer.data(qsg + 1289);
    const auto *qsg_1290 = buffer.data(qsg + 1290);
    const auto *qsg_1291 = buffer.data(qsg + 1291);
    const auto *qsg_1292 = buffer.data(qsg + 1292);
    const auto *qsg_1293 = buffer.data(qsg + 1293);
    const auto *qsg_1294 = buffer.data(qsg + 1294);
    const auto *qsg_1295 = buffer.data(qsg + 1295);
    const auto *qsg_1296 = buffer.data(qsg + 1296);
    const auto *qsg_1297 = buffer.data(qsg + 1297);
    const auto *qsg_1298 = buffer.data(qsg + 1298);
    const auto *qsg_1299 = buffer.data(qsg + 1299);
    const auto *qsg_1300 = buffer.data(qsg + 1300);
    const auto *qsg_1301 = buffer.data(qsg + 1301);
    const auto *qsg_1302 = buffer.data(qsg + 1302);
    const auto *qsg_1303 = buffer.data(qsg + 1303);
    const auto *qsg_1304 = buffer.data(qsg + 1304);
    const auto *qsg_1305 = buffer.data(qsg + 1305);
    const auto *qsg_1306 = buffer.data(qsg + 1306);
    const auto *qsg_1307 = buffer.data(qsg + 1307);
    const auto *qsg_1308 = buffer.data(qsg + 1308);
    const auto *qsg_1309 = buffer.data(qsg + 1309);
    const auto *qsg_1310 = buffer.data(qsg + 1310);
    const auto *qsg_1311 = buffer.data(qsg + 1311);
    const auto *qsg_1312 = buffer.data(qsg + 1312);
    const auto *qsg_1313 = buffer.data(qsg + 1313);
    const auto *qsg_1314 = buffer.data(qsg + 1314);
    const auto *qsg_1315 = buffer.data(qsg + 1315);
    const auto *qsg_1316 = buffer.data(qsg + 1316);
    const auto *qsg_1317 = buffer.data(qsg + 1317);
    const auto *qsg_1318 = buffer.data(qsg + 1318);
    const auto *qsg_1319 = buffer.data(qsg + 1319);
    const auto *qsg_1320 = buffer.data(qsg + 1320);
    const auto *qsg_1321 = buffer.data(qsg + 1321);
    const auto *qsg_1322 = buffer.data(qsg + 1322);
    const auto *qsg_1323 = buffer.data(qsg + 1323);
    const auto *qsg_1324 = buffer.data(qsg + 1324);
    const auto *qsg_1325 = buffer.data(qsg + 1325);
    const auto *qsg_1326 = buffer.data(qsg + 1326);
    const auto *qsg_1327 = buffer.data(qsg + 1327);
    const auto *qsg_1328 = buffer.data(qsg + 1328);
    const auto *qsg_1329 = buffer.data(qsg + 1329);
    const auto *qsg_1330 = buffer.data(qsg + 1330);
    const auto *qsg_1331 = buffer.data(qsg + 1331);
    const auto *qsg_1332 = buffer.data(qsg + 1332);
    const auto *qsg_1333 = buffer.data(qsg + 1333);
    const auto *qsg_1334 = buffer.data(qsg + 1334);
    const auto *qsg_1336 = buffer.data(qsg + 1336);
    const auto *qsg_1338 = buffer.data(qsg + 1338);
    const auto *qsg_1339 = buffer.data(qsg + 1339);

#pragma omp simd aligned(t_1760, t_1761, t_1762, pc_y, osg_1077, osg_1078, osg_1079, qsf0_838, \
                         qsf0_839, qsf1_838, qsf1_839, qsg_1257, qsg_1258, \
                         qsg_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1760[k] = f_19 * osg_1077[k]
                    + f_6 * qsf0_838[k]
                    - f_7 * qsf1_838[k]
                    + f_3 * pc_y[k] * qsg_1257[k];

        t_1761[k] = f_19 * osg_1078[k]
                    + f_4 * qsf0_839[k]
                    - f_5 * qsf1_839[k]
                    + f_3 * pc_y[k] * qsg_1258[k];

        t_1762[k] = f_19 * osg_1079[k]
                    + f_3 * pc_y[k] * qsg_1259[k];
    }

#pragma omp simd aligned(t_1763, t_1764, t_1765, pc_x, pc_z, osg_1064, qsf0_839, qsf0_840, \
                         qsf0_841, qsf1_839, qsf1_840, qsf1_841, qsg_1259, qsg_1260, \
                         qsg_1261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1763[k] = f_20 * osg_1064[k]
                    + f_1 * qsf0_839[k]
                    - f_2 * qsf1_839[k]
                    + f_3 * pc_z[k] * qsg_1259[k];

        t_1764[k] = f_1 * qsf0_840[k]
                    - f_2 * qsf1_840[k]
                    + f_3 * pc_x[k] * qsg_1260[k];

        t_1765[k] = f_13 * qsf0_841[k]
                    - f_14 * qsf1_841[k]
                    + f_3 * pc_x[k] * qsg_1261[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pc_x, qsf0_842, qsf0_843, qsf0_844, qsf1_842, \
                         qsf1_843, qsf1_844, qsg_1262, qsg_1263, \
                         qsg_1264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_13 * qsf0_842[k]
                    - f_14 * qsf1_842[k]
                    + f_3 * pc_x[k] * qsg_1262[k];

        t_1767[k] = f_6 * qsf0_843[k]
                    - f_7 * qsf1_843[k]
                    + f_3 * pc_x[k] * qsg_1263[k];

        t_1768[k] = f_6 * qsf0_844[k]
                    - f_7 * qsf1_844[k]
                    + f_3 * pc_x[k] * qsg_1264[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, pc_x, qsf0_845, qsf0_846, qsf0_847, qsf1_845, \
                         qsf1_846, qsf1_847, qsg_1265, qsg_1266, \
                         qsg_1267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = f_6 * qsf0_845[k]
                    - f_7 * qsf1_845[k]
                    + f_3 * pc_x[k] * qsg_1265[k];

        t_1770[k] = f_4 * qsf0_846[k]
                    - f_5 * qsf1_846[k]
                    + f_3 * pc_x[k] * qsg_1266[k];

        t_1771[k] = f_4 * qsf0_847[k]
                    - f_5 * qsf1_847[k]
                    + f_3 * pc_x[k] * qsg_1267[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, t_1775, t_1776, pc_x, qsf0_848, qsf0_849, \
                         qsf1_848, qsf1_849, qsg_1268, qsg_1269, qsg_1270, qsg_1271, \
                         qsg_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = f_4 * qsf0_848[k]
                    - f_5 * qsf1_848[k]
                    + f_3 * pc_x[k] * qsg_1268[k];

        t_1773[k] = f_4 * qsf0_849[k]
                    - f_5 * qsf1_849[k]
                    + f_3 * pc_x[k] * qsg_1269[k];

        t_1774[k] = f_3 * pc_x[k] * qsg_1270[k];

        t_1775[k] = f_3 * pc_x[k] * qsg_1271[k];

        t_1776[k] = f_3 * pc_x[k] * qsg_1272[k];
    }

#pragma omp simd aligned(t_1777, t_1778, t_1779, t_1780, pc_x, pc_y, pc_z, osg_1075, osg_1090, \
                         qsf0_846, qsf1_846, qsg_1270, qsg_1273, \
                         qsg_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1777[k] = f_3 * pc_x[k] * qsg_1273[k];

        t_1778[k] = f_3 * pc_x[k] * qsg_1274[k];

        t_1779[k] = f_21 * osg_1090[k]
                    + f_1 * qsf0_846[k]
                    - f_2 * qsf1_846[k]
                    + f_3 * pc_y[k] * qsg_1270[k];

        t_1780[k] = f_21 * osg_1075[k]
                    + f_3 * pc_z[k] * qsg_1270[k];
    }

#pragma omp simd aligned(t_1781, t_1782, t_1783, pc_y, osg_1092, osg_1093, osg_1094, qsf0_848, \
                         qsf0_849, qsf1_848, qsf1_849, qsg_1272, qsg_1273, \
                         qsg_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1781[k] = f_21 * osg_1092[k]
                    + f_6 * qsf0_848[k]
                    - f_7 * qsf1_848[k]
                    + f_3 * pc_y[k] * qsg_1272[k];

        t_1782[k] = f_21 * osg_1093[k]
                    + f_4 * qsf0_849[k]
                    - f_5 * qsf1_849[k]
                    + f_3 * pc_y[k] * qsg_1273[k];

        t_1783[k] = f_21 * osg_1094[k]
                    + f_3 * pc_y[k] * qsg_1274[k];
    }

#pragma omp simd aligned(t_1784, t_1785, t_1786, pc_x, pc_z, osg_1079, qsf0_849, qsf0_850, \
                         qsf0_851, qsf1_849, qsf1_850, qsf1_851, qsg_1274, qsg_1275, \
                         qsg_1276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1784[k] = f_21 * osg_1079[k]
                    + f_1 * qsf0_849[k]
                    - f_2 * qsf1_849[k]
                    + f_3 * pc_z[k] * qsg_1274[k];

        t_1785[k] = f_1 * qsf0_850[k]
                    - f_2 * qsf1_850[k]
                    + f_3 * pc_x[k] * qsg_1275[k];

        t_1786[k] = f_13 * qsf0_851[k]
                    - f_14 * qsf1_851[k]
                    + f_3 * pc_x[k] * qsg_1276[k];
    }

#pragma omp simd aligned(t_1787, t_1788, t_1789, pc_x, qsf0_852, qsf0_853, qsf0_854, qsf1_852, \
                         qsf1_853, qsf1_854, qsg_1277, qsg_1278, \
                         qsg_1279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1787[k] = f_13 * qsf0_852[k]
                    - f_14 * qsf1_852[k]
                    + f_3 * pc_x[k] * qsg_1277[k];

        t_1788[k] = f_6 * qsf0_853[k]
                    - f_7 * qsf1_853[k]
                    + f_3 * pc_x[k] * qsg_1278[k];

        t_1789[k] = f_6 * qsf0_854[k]
                    - f_7 * qsf1_854[k]
                    + f_3 * pc_x[k] * qsg_1279[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, pc_x, qsf0_855, qsf0_856, qsf0_857, qsf1_855, \
                         qsf1_856, qsf1_857, qsg_1280, qsg_1281, \
                         qsg_1282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_6 * qsf0_855[k]
                    - f_7 * qsf1_855[k]
                    + f_3 * pc_x[k] * qsg_1280[k];

        t_1791[k] = f_4 * qsf0_856[k]
                    - f_5 * qsf1_856[k]
                    + f_3 * pc_x[k] * qsg_1281[k];

        t_1792[k] = f_4 * qsf0_857[k]
                    - f_5 * qsf1_857[k]
                    + f_3 * pc_x[k] * qsg_1282[k];
    }

#pragma omp simd aligned(t_1793, t_1794, t_1795, t_1796, t_1797, pc_x, qsf0_858, qsf0_859, \
                         qsf1_858, qsf1_859, qsg_1283, qsg_1284, qsg_1285, qsg_1286, \
                         qsg_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1793[k] = f_4 * qsf0_858[k]
                    - f_5 * qsf1_858[k]
                    + f_3 * pc_x[k] * qsg_1283[k];

        t_1794[k] = f_4 * qsf0_859[k]
                    - f_5 * qsf1_859[k]
                    + f_3 * pc_x[k] * qsg_1284[k];

        t_1795[k] = f_3 * pc_x[k] * qsg_1285[k];

        t_1796[k] = f_3 * pc_x[k] * qsg_1286[k];

        t_1797[k] = f_3 * pc_x[k] * qsg_1287[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, t_1801, pc_x, pc_y, pc_z, osg_1090, osg_1105, \
                         qsf0_856, qsf1_856, qsg_1285, qsg_1288, \
                         qsg_1289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = f_3 * pc_x[k] * qsg_1288[k];

        t_1799[k] = f_3 * pc_x[k] * qsg_1289[k];

        t_1800[k] = f_20 * osg_1105[k]
                    + f_1 * qsf0_856[k]
                    - f_2 * qsf1_856[k]
                    + f_3 * pc_y[k] * qsg_1285[k];

        t_1801[k] = f_19 * osg_1090[k]
                    + f_3 * pc_z[k] * qsg_1285[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pc_y, osg_1107, osg_1108, osg_1109, qsf0_858, \
                         qsf0_859, qsf1_858, qsf1_859, qsg_1287, qsg_1288, \
                         qsg_1289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_20 * osg_1107[k]
                    + f_6 * qsf0_858[k]
                    - f_7 * qsf1_858[k]
                    + f_3 * pc_y[k] * qsg_1287[k];

        t_1803[k] = f_20 * osg_1108[k]
                    + f_4 * qsf0_859[k]
                    - f_5 * qsf1_859[k]
                    + f_3 * pc_y[k] * qsg_1288[k];

        t_1804[k] = f_20 * osg_1109[k]
                    + f_3 * pc_y[k] * qsg_1289[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, pc_x, pc_z, osg_1094, qsf0_859, qsf0_860, \
                         qsf0_861, qsf1_859, qsf1_860, qsf1_861, qsg_1289, qsg_1290, \
                         qsg_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_19 * osg_1094[k]
                    + f_1 * qsf0_859[k]
                    - f_2 * qsf1_859[k]
                    + f_3 * pc_z[k] * qsg_1289[k];

        t_1806[k] = f_1 * qsf0_860[k]
                    - f_2 * qsf1_860[k]
                    + f_3 * pc_x[k] * qsg_1290[k];

        t_1807[k] = f_13 * qsf0_861[k]
                    - f_14 * qsf1_861[k]
                    + f_3 * pc_x[k] * qsg_1291[k];
    }

#pragma omp simd aligned(t_1808, t_1809, t_1810, pc_x, qsf0_862, qsf0_863, qsf0_864, qsf1_862, \
                         qsf1_863, qsf1_864, qsg_1292, qsg_1293, \
                         qsg_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1808[k] = f_13 * qsf0_862[k]
                    - f_14 * qsf1_862[k]
                    + f_3 * pc_x[k] * qsg_1292[k];

        t_1809[k] = f_6 * qsf0_863[k]
                    - f_7 * qsf1_863[k]
                    + f_3 * pc_x[k] * qsg_1293[k];

        t_1810[k] = f_6 * qsf0_864[k]
                    - f_7 * qsf1_864[k]
                    + f_3 * pc_x[k] * qsg_1294[k];
    }

#pragma omp simd aligned(t_1811, t_1812, t_1813, pc_x, qsf0_865, qsf0_866, qsf0_867, qsf1_865, \
                         qsf1_866, qsf1_867, qsg_1295, qsg_1296, \
                         qsg_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1811[k] = f_6 * qsf0_865[k]
                    - f_7 * qsf1_865[k]
                    + f_3 * pc_x[k] * qsg_1295[k];

        t_1812[k] = f_4 * qsf0_866[k]
                    - f_5 * qsf1_866[k]
                    + f_3 * pc_x[k] * qsg_1296[k];

        t_1813[k] = f_4 * qsf0_867[k]
                    - f_5 * qsf1_867[k]
                    + f_3 * pc_x[k] * qsg_1297[k];
    }

#pragma omp simd aligned(t_1814, t_1815, t_1816, t_1817, t_1818, pc_x, qsf0_868, qsf0_869, \
                         qsf1_868, qsf1_869, qsg_1298, qsg_1299, qsg_1300, qsg_1301, \
                         qsg_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = f_4 * qsf0_868[k]
                    - f_5 * qsf1_868[k]
                    + f_3 * pc_x[k] * qsg_1298[k];

        t_1815[k] = f_4 * qsf0_869[k]
                    - f_5 * qsf1_869[k]
                    + f_3 * pc_x[k] * qsg_1299[k];

        t_1816[k] = f_3 * pc_x[k] * qsg_1300[k];

        t_1817[k] = f_3 * pc_x[k] * qsg_1301[k];

        t_1818[k] = f_3 * pc_x[k] * qsg_1302[k];
    }

#pragma omp simd aligned(t_1819, t_1820, t_1821, t_1822, pc_x, pc_y, pc_z, osg_1105, osg_1120, \
                         qsf0_866, qsf1_866, qsg_1300, qsg_1303, \
                         qsg_1304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1819[k] = f_3 * pc_x[k] * qsg_1303[k];

        t_1820[k] = f_3 * pc_x[k] * qsg_1304[k];

        t_1821[k] = f_18 * osg_1120[k]
                    + f_1 * qsf0_866[k]
                    - f_2 * qsf1_866[k]
                    + f_3 * pc_y[k] * qsg_1300[k];

        t_1822[k] = f_17 * osg_1105[k]
                    + f_3 * pc_z[k] * qsg_1300[k];
    }

#pragma omp simd aligned(t_1823, t_1824, t_1825, pc_y, osg_1122, osg_1123, osg_1124, qsf0_868, \
                         qsf0_869, qsf1_868, qsf1_869, qsg_1302, qsg_1303, \
                         qsg_1304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1823[k] = f_18 * osg_1122[k]
                    + f_6 * qsf0_868[k]
                    - f_7 * qsf1_868[k]
                    + f_3 * pc_y[k] * qsg_1302[k];

        t_1824[k] = f_18 * osg_1123[k]
                    + f_4 * qsf0_869[k]
                    - f_5 * qsf1_869[k]
                    + f_3 * pc_y[k] * qsg_1303[k];

        t_1825[k] = f_18 * osg_1124[k]
                    + f_3 * pc_y[k] * qsg_1304[k];
    }

#pragma omp simd aligned(t_1826, t_1827, t_1828, pc_x, pc_z, osg_1109, qsf0_869, qsf0_870, \
                         qsf0_871, qsf1_869, qsf1_870, qsf1_871, qsg_1304, qsg_1305, \
                         qsg_1306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1826[k] = f_17 * osg_1109[k]
                    + f_1 * qsf0_869[k]
                    - f_2 * qsf1_869[k]
                    + f_3 * pc_z[k] * qsg_1304[k];

        t_1827[k] = f_1 * qsf0_870[k]
                    - f_2 * qsf1_870[k]
                    + f_3 * pc_x[k] * qsg_1305[k];

        t_1828[k] = f_13 * qsf0_871[k]
                    - f_14 * qsf1_871[k]
                    + f_3 * pc_x[k] * qsg_1306[k];
    }

#pragma omp simd aligned(t_1829, t_1830, t_1831, pc_x, qsf0_872, qsf0_873, qsf0_874, qsf1_872, \
                         qsf1_873, qsf1_874, qsg_1307, qsg_1308, \
                         qsg_1309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1829[k] = f_13 * qsf0_872[k]
                    - f_14 * qsf1_872[k]
                    + f_3 * pc_x[k] * qsg_1307[k];

        t_1830[k] = f_6 * qsf0_873[k]
                    - f_7 * qsf1_873[k]
                    + f_3 * pc_x[k] * qsg_1308[k];

        t_1831[k] = f_6 * qsf0_874[k]
                    - f_7 * qsf1_874[k]
                    + f_3 * pc_x[k] * qsg_1309[k];
    }

#pragma omp simd aligned(t_1832, t_1833, t_1834, pc_x, qsf0_875, qsf0_876, qsf0_877, qsf1_875, \
                         qsf1_876, qsf1_877, qsg_1310, qsg_1311, \
                         qsg_1312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1832[k] = f_6 * qsf0_875[k]
                    - f_7 * qsf1_875[k]
                    + f_3 * pc_x[k] * qsg_1310[k];

        t_1833[k] = f_4 * qsf0_876[k]
                    - f_5 * qsf1_876[k]
                    + f_3 * pc_x[k] * qsg_1311[k];

        t_1834[k] = f_4 * qsf0_877[k]
                    - f_5 * qsf1_877[k]
                    + f_3 * pc_x[k] * qsg_1312[k];
    }

#pragma omp simd aligned(t_1835, t_1836, t_1837, t_1838, t_1839, pc_x, qsf0_878, qsf0_879, \
                         qsf1_878, qsf1_879, qsg_1313, qsg_1314, qsg_1315, qsg_1316, \
                         qsg_1317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1835[k] = f_4 * qsf0_878[k]
                    - f_5 * qsf1_878[k]
                    + f_3 * pc_x[k] * qsg_1313[k];

        t_1836[k] = f_4 * qsf0_879[k]
                    - f_5 * qsf1_879[k]
                    + f_3 * pc_x[k] * qsg_1314[k];

        t_1837[k] = f_3 * pc_x[k] * qsg_1315[k];

        t_1838[k] = f_3 * pc_x[k] * qsg_1316[k];

        t_1839[k] = f_3 * pc_x[k] * qsg_1317[k];
    }

#pragma omp simd aligned(t_1840, t_1841, t_1842, t_1843, pc_x, pc_y, pc_z, osg_1120, osg_1135, \
                         qsf0_876, qsf1_876, qsg_1315, qsg_1318, \
                         qsg_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1840[k] = f_3 * pc_x[k] * qsg_1318[k];

        t_1841[k] = f_3 * pc_x[k] * qsg_1319[k];

        t_1842[k] = f_11 * osg_1135[k]
                    + f_1 * qsf0_876[k]
                    - f_2 * qsf1_876[k]
                    + f_3 * pc_y[k] * qsg_1315[k];

        t_1843[k] = f_16 * osg_1120[k]
                    + f_3 * pc_z[k] * qsg_1315[k];
    }

#pragma omp simd aligned(t_1844, t_1845, t_1846, pc_y, osg_1137, osg_1138, osg_1139, qsf0_878, \
                         qsf0_879, qsf1_878, qsf1_879, qsg_1317, qsg_1318, \
                         qsg_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1844[k] = f_11 * osg_1137[k]
                    + f_6 * qsf0_878[k]
                    - f_7 * qsf1_878[k]
                    + f_3 * pc_y[k] * qsg_1317[k];

        t_1845[k] = f_11 * osg_1138[k]
                    + f_4 * qsf0_879[k]
                    - f_5 * qsf1_879[k]
                    + f_3 * pc_y[k] * qsg_1318[k];

        t_1846[k] = f_11 * osg_1139[k]
                    + f_3 * pc_y[k] * qsg_1319[k];
    }

#pragma omp simd aligned(t_1847, t_1848, t_1849, pc_x, pc_z, osg_1124, qsf0_879, qsf0_880, \
                         qsf0_881, qsf1_879, qsf1_880, qsf1_881, qsg_1319, qsg_1320, \
                         qsg_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1847[k] = f_16 * osg_1124[k]
                    + f_1 * qsf0_879[k]
                    - f_2 * qsf1_879[k]
                    + f_3 * pc_z[k] * qsg_1319[k];

        t_1848[k] = f_1 * qsf0_880[k]
                    - f_2 * qsf1_880[k]
                    + f_3 * pc_x[k] * qsg_1320[k];

        t_1849[k] = f_13 * qsf0_881[k]
                    - f_14 * qsf1_881[k]
                    + f_3 * pc_x[k] * qsg_1321[k];
    }

#pragma omp simd aligned(t_1850, t_1851, t_1852, pc_x, qsf0_882, qsf0_883, qsf0_884, qsf1_882, \
                         qsf1_883, qsf1_884, qsg_1322, qsg_1323, \
                         qsg_1324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1850[k] = f_13 * qsf0_882[k]
                    - f_14 * qsf1_882[k]
                    + f_3 * pc_x[k] * qsg_1322[k];

        t_1851[k] = f_6 * qsf0_883[k]
                    - f_7 * qsf1_883[k]
                    + f_3 * pc_x[k] * qsg_1323[k];

        t_1852[k] = f_6 * qsf0_884[k]
                    - f_7 * qsf1_884[k]
                    + f_3 * pc_x[k] * qsg_1324[k];
    }

#pragma omp simd aligned(t_1853, t_1854, t_1855, pc_x, qsf0_885, qsf0_886, qsf0_887, qsf1_885, \
                         qsf1_886, qsf1_887, qsg_1325, qsg_1326, \
                         qsg_1327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1853[k] = f_6 * qsf0_885[k]
                    - f_7 * qsf1_885[k]
                    + f_3 * pc_x[k] * qsg_1325[k];

        t_1854[k] = f_4 * qsf0_886[k]
                    - f_5 * qsf1_886[k]
                    + f_3 * pc_x[k] * qsg_1326[k];

        t_1855[k] = f_4 * qsf0_887[k]
                    - f_5 * qsf1_887[k]
                    + f_3 * pc_x[k] * qsg_1327[k];
    }

#pragma omp simd aligned(t_1856, t_1857, t_1858, t_1859, t_1860, pc_x, qsf0_888, qsf0_889, \
                         qsf1_888, qsf1_889, qsg_1328, qsg_1329, qsg_1330, qsg_1331, \
                         qsg_1332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1856[k] = f_4 * qsf0_888[k]
                    - f_5 * qsf1_888[k]
                    + f_3 * pc_x[k] * qsg_1328[k];

        t_1857[k] = f_4 * qsf0_889[k]
                    - f_5 * qsf1_889[k]
                    + f_3 * pc_x[k] * qsg_1329[k];

        t_1858[k] = f_3 * pc_x[k] * qsg_1330[k];

        t_1859[k] = f_3 * pc_x[k] * qsg_1331[k];

        t_1860[k] = f_3 * pc_x[k] * qsg_1332[k];
    }

#pragma omp simd aligned(t_1861, t_1862, t_1863, t_1864, pc_x, pc_y, pc_z, osg_1135, osg_1150, \
                         qsf0_886, qsf1_886, qsg_1330, qsg_1333, \
                         qsg_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1861[k] = f_3 * pc_x[k] * qsg_1333[k];

        t_1862[k] = f_3 * pc_x[k] * qsg_1334[k];

        t_1863[k] = f_10 * osg_1150[k]
                    + f_1 * qsf0_886[k]
                    - f_2 * qsf1_886[k]
                    + f_3 * pc_y[k] * qsg_1330[k];

        t_1864[k] = f_15 * osg_1135[k]
                    + f_3 * pc_z[k] * qsg_1330[k];
    }

#pragma omp simd aligned(t_1865, t_1866, t_1867, pc_y, osg_1152, osg_1153, osg_1154, qsf0_888, \
                         qsf0_889, qsf1_888, qsf1_889, qsg_1332, qsg_1333, \
                         qsg_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1865[k] = f_10 * osg_1152[k]
                    + f_6 * qsf0_888[k]
                    - f_7 * qsf1_888[k]
                    + f_3 * pc_y[k] * qsg_1332[k];

        t_1866[k] = f_10 * osg_1153[k]
                    + f_4 * qsf0_889[k]
                    - f_5 * qsf1_889[k]
                    + f_3 * pc_y[k] * qsg_1333[k];

        t_1867[k] = f_10 * osg_1154[k]
                    + f_3 * pc_y[k] * qsg_1334[k];
    }

#pragma omp simd aligned(t_1868, t_1869, t_1870, pa_y, pc_x, pc_y, pc_z, osh0_1617, osg_1139, \
                         osh1_1617, qsf0_889, qsf0_891, qsf1_889, qsf1_891, qsg_1334, \
                         qsg_1336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1868[k] = f_15 * osg_1139[k]
                    + f_1 * qsf0_889[k]
                    - f_2 * qsf1_889[k]
                    + f_3 * pc_z[k] * qsg_1334[k];

        t_1869[k] = pa_y[k] * osh0_1617[k]
                    - f_8 * pc_y[k] * osh1_1617[k];

        t_1870[k] = f_13 * qsf0_891[k]
                    - f_14 * qsf1_891[k]
                    + f_3 * pc_x[k] * qsg_1336[k];
    }

#pragma omp simd aligned(t_1871, t_1872, t_1873, pa_y, pc_x, pc_y, osh0_1619, osh1_1619, \
                         qsf0_893, qsf0_894, qsf1_893, qsf1_894, qsg_1338, \
                         qsg_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1871[k] = pa_y[k] * osh0_1619[k]
                    - f_8 * pc_y[k] * osh1_1619[k];

        t_1872[k] = f_6 * qsf0_893[k]
                    - f_7 * qsf1_893[k]
                    + f_3 * pc_x[k] * qsg_1338[k];

        t_1873[k] = f_6 * qsf0_894[k]
                    - f_7 * qsf1_894[k]
                    + f_3 * pc_x[k] * qsg_1339[k];
    }
}

static auto
compute_prim_qsh_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osh0,
                                                           const size_t osg, const size_t osh1,
                                                           const size_t qsf0, const size_t qsf1,
                                                           const size_t qsg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
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
    const auto f_12 = 5.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_20 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh0_1622 = buffer.data(osh0 + 1622);
    const auto *osh0_1626 = buffer.data(osh0 + 1626);
    const auto *osh0_1632 = buffer.data(osh0 + 1632);
    const auto *osh0_1634 = buffer.data(osh0 + 1634);
    const auto *osh0_1635 = buffer.data(osh0 + 1635);
    const auto *osh0_1637 = buffer.data(osh0 + 1637);

    const auto *osg_1150 = buffer.data(osg + 1150);
    const auto *osg_1165 = buffer.data(osg + 1165);
    const auto *osg_1167 = buffer.data(osg + 1167);
    const auto *osg_1168 = buffer.data(osg + 1168);
    const auto *osg_1169 = buffer.data(osg + 1169);

    const auto *osh1_1622 = buffer.data(osh1 + 1622);
    const auto *osh1_1626 = buffer.data(osh1 + 1626);
    const auto *osh1_1632 = buffer.data(osh1 + 1632);
    const auto *osh1_1634 = buffer.data(osh1 + 1634);
    const auto *osh1_1635 = buffer.data(osh1 + 1635);
    const auto *osh1_1637 = buffer.data(osh1 + 1637);

    const auto *qsf0_896 = buffer.data(qsf0 + 896);
    const auto *qsf0_897 = buffer.data(qsf0 + 897);
    const auto *qsf0_898 = buffer.data(qsf0 + 898);
    const auto *qsf0_900 = buffer.data(qsf0 + 900);
    const auto *qsf0_902 = buffer.data(qsf0 + 902);
    const auto *qsf0_903 = buffer.data(qsf0 + 903);
    const auto *qsf0_905 = buffer.data(qsf0 + 905);
    const auto *qsf0_906 = buffer.data(qsf0 + 906);
    const auto *qsf0_907 = buffer.data(qsf0 + 907);
    const auto *qsf0_908 = buffer.data(qsf0 + 908);
    const auto *qsf0_909 = buffer.data(qsf0 + 909);

    const auto *qsf1_896 = buffer.data(qsf1 + 896);
    const auto *qsf1_897 = buffer.data(qsf1 + 897);
    const auto *qsf1_898 = buffer.data(qsf1 + 898);
    const auto *qsf1_900 = buffer.data(qsf1 + 900);
    const auto *qsf1_902 = buffer.data(qsf1 + 902);
    const auto *qsf1_903 = buffer.data(qsf1 + 903);
    const auto *qsf1_905 = buffer.data(qsf1 + 905);
    const auto *qsf1_906 = buffer.data(qsf1 + 906);
    const auto *qsf1_907 = buffer.data(qsf1 + 907);
    const auto *qsf1_908 = buffer.data(qsf1 + 908);
    const auto *qsf1_909 = buffer.data(qsf1 + 909);

    const auto *qsg_1341 = buffer.data(qsg + 1341);
    const auto *qsg_1342 = buffer.data(qsg + 1342);
    const auto *qsg_1343 = buffer.data(qsg + 1343);
    const auto *qsg_1345 = buffer.data(qsg + 1345);
    const auto *qsg_1346 = buffer.data(qsg + 1346);
    const auto *qsg_1347 = buffer.data(qsg + 1347);
    const auto *qsg_1348 = buffer.data(qsg + 1348);
    const auto *qsg_1349 = buffer.data(qsg + 1349);
    const auto *qsg_1350 = buffer.data(qsg + 1350);
    const auto *qsg_1352 = buffer.data(qsg + 1352);
    const auto *qsg_1353 = buffer.data(qsg + 1353);
    const auto *qsg_1355 = buffer.data(qsg + 1355);
    const auto *qsg_1356 = buffer.data(qsg + 1356);
    const auto *qsg_1357 = buffer.data(qsg + 1357);
    const auto *qsg_1359 = buffer.data(qsg + 1359);
    const auto *qsg_1360 = buffer.data(qsg + 1360);
    const auto *qsg_1361 = buffer.data(qsg + 1361);
    const auto *qsg_1362 = buffer.data(qsg + 1362);
    const auto *qsg_1363 = buffer.data(qsg + 1363);
    const auto *qsg_1364 = buffer.data(qsg + 1364);

#pragma omp simd aligned(t_1874, t_1875, t_1876, pa_y, pc_x, pc_y, osh0_1622, osh1_1622, \
                         qsf0_896, qsf0_897, qsf1_896, qsf1_897, qsg_1341, \
                         qsg_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1874[k] = pa_y[k] * osh0_1622[k]
                    - f_8 * pc_y[k] * osh1_1622[k];

        t_1875[k] = f_4 * qsf0_896[k]
                    - f_5 * qsf1_896[k]
                    + f_3 * pc_x[k] * qsg_1341[k];

        t_1876[k] = f_4 * qsf0_897[k]
                    - f_5 * qsf1_897[k]
                    + f_3 * pc_x[k] * qsg_1342[k];
    }

#pragma omp simd aligned(t_1877, t_1878, t_1879, t_1880, t_1881, pa_y, pc_x, pc_y, osh0_1626, \
                         osh1_1626, qsf0_898, qsf1_898, qsg_1343, qsg_1345, qsg_1346, \
                         qsg_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1877[k] = f_4 * qsf0_898[k]
                    - f_5 * qsf1_898[k]
                    + f_3 * pc_x[k] * qsg_1343[k];

        t_1878[k] = pa_y[k] * osh0_1626[k]
                    - f_8 * pc_y[k] * osh1_1626[k];

        t_1879[k] = f_3 * pc_x[k] * qsg_1345[k];

        t_1880[k] = f_3 * pc_x[k] * qsg_1346[k];

        t_1881[k] = f_3 * pc_x[k] * qsg_1347[k];
    }

#pragma omp simd aligned(t_1882, t_1883, t_1884, t_1885, pa_y, pc_x, pc_y, pc_z, osh0_1632, \
                         osg_1150, osg_1165, osh1_1632, qsg_1345, qsg_1348, \
                         qsg_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1882[k] = f_3 * pc_x[k] * qsg_1348[k];

        t_1883[k] = f_3 * pc_x[k] * qsg_1349[k];

        t_1884[k] = pa_y[k] * osh0_1632[k]
                    + f_20 * osg_1165[k]
                    - f_8 * pc_y[k] * osh1_1632[k];

        t_1885[k] = f_12 * osg_1150[k]
                    + f_3 * pc_z[k] * qsg_1345[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, t_1889, pa_y, pc_y, osh0_1634, osh0_1635, \
                         osh0_1637, osg_1167, osg_1168, osg_1169, osh1_1634, osh1_1635, \
                         osh1_1637, qsg_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = pa_y[k] * osh0_1634[k]
                    + f_11 * osg_1167[k]
                    - f_8 * pc_y[k] * osh1_1634[k];

        t_1887[k] = pa_y[k] * osh0_1635[k]
                    + f_10 * osg_1168[k]
                    - f_8 * pc_y[k] * osh1_1635[k];

        t_1888[k] = f_9 * osg_1169[k]
                    + f_3 * pc_y[k] * qsg_1349[k];

        t_1889[k] = pa_y[k] * osh0_1637[k]
                    - f_8 * pc_y[k] * osh1_1637[k];
    }

#pragma omp simd aligned(t_1890, t_1891, t_1892, t_1893, t_1894, pc_x, pc_y, qsf0_900, \
                         qsf0_902, qsf0_903, qsf1_900, qsf1_902, qsf1_903, qsg_1350, qsg_1352, \
                         qsg_1353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1890[k] = f_1 * qsf0_900[k]
                    - f_2 * qsf1_900[k]
                    + f_3 * pc_x[k] * qsg_1350[k];

        t_1891[k] = f_3 * pc_y[k] * qsg_1350[k];

        t_1892[k] = f_13 * qsf0_902[k]
                    - f_14 * qsf1_902[k]
                    + f_3 * pc_x[k] * qsg_1352[k];

        t_1893[k] = f_6 * qsf0_903[k]
                    - f_7 * qsf1_903[k]
                    + f_3 * pc_x[k] * qsg_1353[k];

        t_1894[k] = f_3 * pc_y[k] * qsg_1352[k];
    }

#pragma omp simd aligned(t_1895, t_1896, t_1897, t_1898, pc_x, pc_y, qsf0_905, qsf0_906, \
                         qsf0_907, qsf1_905, qsf1_906, qsf1_907, qsg_1355, qsg_1356, \
                         qsg_1357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1895[k] = f_6 * qsf0_905[k]
                    - f_7 * qsf1_905[k]
                    + f_3 * pc_x[k] * qsg_1355[k];

        t_1896[k] = f_4 * qsf0_906[k]
                    - f_5 * qsf1_906[k]
                    + f_3 * pc_x[k] * qsg_1356[k];

        t_1897[k] = f_4 * qsf0_907[k]
                    - f_5 * qsf1_907[k]
                    + f_3 * pc_x[k] * qsg_1357[k];

        t_1898[k] = f_3 * pc_y[k] * qsg_1355[k];
    }

#pragma omp simd aligned(t_1899, t_1900, t_1901, t_1902, t_1903, t_1904, pc_x, qsf0_909, \
                         qsf1_909, qsg_1359, qsg_1360, qsg_1361, qsg_1362, qsg_1363, \
                         qsg_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1899[k] = f_4 * qsf0_909[k]
                    - f_5 * qsf1_909[k]
                    + f_3 * pc_x[k] * qsg_1359[k];

        t_1900[k] = f_3 * pc_x[k] * qsg_1360[k];

        t_1901[k] = f_3 * pc_x[k] * qsg_1361[k];

        t_1902[k] = f_3 * pc_x[k] * qsg_1362[k];

        t_1903[k] = f_3 * pc_x[k] * qsg_1363[k];

        t_1904[k] = f_3 * pc_x[k] * qsg_1364[k];
    }

#pragma omp simd aligned(t_1905, t_1906, t_1907, pc_y, qsf0_906, qsf0_907, qsf0_908, qsf1_906, \
                         qsf1_907, qsf1_908, qsg_1360, qsg_1361, \
                         qsg_1362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1905[k] = f_1 * qsf0_906[k]
                    - f_2 * qsf1_906[k]
                    + f_3 * pc_y[k] * qsg_1360[k];

        t_1906[k] = f_13 * qsf0_907[k]
                    - f_14 * qsf1_907[k]
                    + f_3 * pc_y[k] * qsg_1361[k];

        t_1907[k] = f_6 * qsf0_908[k]
                    - f_7 * qsf1_908[k]
                    + f_3 * pc_y[k] * qsg_1362[k];
    }

#pragma omp simd aligned(t_1908, t_1909, t_1910, pc_y, pc_z, osg_1169, qsf0_909, qsf1_909, \
                         qsg_1363, qsg_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1908[k] = f_4 * qsf0_909[k]
                    - f_5 * qsf1_909[k]
                    + f_3 * pc_y[k] * qsg_1363[k];

        t_1909[k] = f_3 * pc_y[k] * qsg_1364[k];

        t_1910[k] = f_0 * osg_1169[k]
                    + f_1 * qsf0_909[k]
                    - f_2 * qsf1_909[k]
                    + f_3 * pc_z[k] * qsg_1364[k];
    }
}

auto
compute_prim_qsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t osh0, const size_t osg,
                                                   const size_t osh1, const size_t qsf0,
                                                   const size_t qsf1, const size_t qsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_qsh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, osh0, osg,
                                                              osh1, qsf0, qsf1, qsg, ncols,
                                                              gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsf0, qsf1, qsg,
                                                               ncols, gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsf0, qsf1, qsg,
                                                               ncols, gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece12(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsg, ncols, gamma, p,
                                                               q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsf0, qsf1, qsg,
                                                               ncols, gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece14(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsf0, qsf1, qsg,
                                                               ncols, gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece15(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsf0, qsf1, qsg,
                                                               ncols, gamma, p, q);

    compute_prim_qsh_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, osh0,
                                                               osg, osh1, qsf0, qsf1, qsg,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
