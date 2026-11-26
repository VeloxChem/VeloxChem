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


#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *osg0_0 = buffer.data(osg0 + 0);
    const auto *osg0_3 = buffer.data(osg0 + 3);
    const auto *osg0_5 = buffer.data(osg0 + 5);
    const auto *osg0_10 = buffer.data(osg0 + 10);
    const auto *osg0_14 = buffer.data(osg0 + 14);
    const auto *osg0_18 = buffer.data(osg0 + 18);
    const auto *osg0_25 = buffer.data(osg0 + 25);
    const auto *osg0_30 = buffer.data(osg0 + 30);
    const auto *osg0_35 = buffer.data(osg0 + 35);
    const auto *osg0_44 = buffer.data(osg0 + 44);
    const auto *osg0_45 = buffer.data(osg0 + 45);
    const auto *osg0_48 = buffer.data(osg0 + 48);
    const auto *osg0_55 = buffer.data(osg0 + 55);
    const auto *osg0_75 = buffer.data(osg0 + 75);
    const auto *osg0_78 = buffer.data(osg0 + 78);
    const auto *osg0_80 = buffer.data(osg0 + 80);

    const auto *osf_0 = buffer.data(osf + 0);
    const auto *osf_1 = buffer.data(osf + 1);
    const auto *osf_2 = buffer.data(osf + 2);
    const auto *osf_6 = buffer.data(osf + 6);
    const auto *osf_9 = buffer.data(osf + 9);
    const auto *osf_10 = buffer.data(osf + 10);
    const auto *osf_16 = buffer.data(osf + 16);
    const auto *osf_18 = buffer.data(osf + 18);
    const auto *osf_19 = buffer.data(osf + 19);
    const auto *osf_20 = buffer.data(osf + 20);
    const auto *osf_22 = buffer.data(osf + 22);
    const auto *osf_26 = buffer.data(osf + 26);
    const auto *osf_27 = buffer.data(osf + 27);
    const auto *osf_28 = buffer.data(osf + 28);
    const auto *osf_29 = buffer.data(osf + 29);
    const auto *osf_30 = buffer.data(osf + 30);
    const auto *osf_33 = buffer.data(osf + 33);
    const auto *osf_36 = buffer.data(osf + 36);
    const auto *osf_38 = buffer.data(osf + 38);
    const auto *osf_39 = buffer.data(osf + 39);
    const auto *osf_40 = buffer.data(osf + 40);
    const auto *osf_42 = buffer.data(osf + 42);
    const auto *osf_46 = buffer.data(osf + 46);
    const auto *osf_47 = buffer.data(osf + 47);
    const auto *osf_48 = buffer.data(osf + 48);
    const auto *osf_49 = buffer.data(osf + 49);
    const auto *osf_50 = buffer.data(osf + 50);
    const auto *osf_51 = buffer.data(osf + 51);
    const auto *osf_52 = buffer.data(osf + 52);
    const auto *osf_55 = buffer.data(osf + 55);
    const auto *osf_56 = buffer.data(osf + 56);
    const auto *osf_57 = buffer.data(osf + 57);
    const auto *osf_59 = buffer.data(osf + 59);
    const auto *osf_60 = buffer.data(osf + 60);
    const auto *osf_63 = buffer.data(osf + 63);
    const auto *osf_66 = buffer.data(osf + 66);
    const auto *osf_68 = buffer.data(osf + 68);
    const auto *osf_69 = buffer.data(osf + 69);
    const auto *osf_75 = buffer.data(osf + 75);
    const auto *osf_76 = buffer.data(osf + 76);
    const auto *osf_77 = buffer.data(osf + 77);
    const auto *osf_78 = buffer.data(osf + 78);
    const auto *osf_79 = buffer.data(osf + 79);
    const auto *osf_86 = buffer.data(osf + 86);
    const auto *osf_87 = buffer.data(osf + 87);
    const auto *osf_88 = buffer.data(osf + 88);
    const auto *osf_89 = buffer.data(osf + 89);

    const auto *osg1_0 = buffer.data(osg1 + 0);
    const auto *osg1_3 = buffer.data(osg1 + 3);
    const auto *osg1_5 = buffer.data(osg1 + 5);
    const auto *osg1_10 = buffer.data(osg1 + 10);
    const auto *osg1_14 = buffer.data(osg1 + 14);
    const auto *osg1_18 = buffer.data(osg1 + 18);
    const auto *osg1_25 = buffer.data(osg1 + 25);
    const auto *osg1_30 = buffer.data(osg1 + 30);
    const auto *osg1_35 = buffer.data(osg1 + 35);
    const auto *osg1_44 = buffer.data(osg1 + 44);
    const auto *osg1_45 = buffer.data(osg1 + 45);
    const auto *osg1_48 = buffer.data(osg1 + 48);
    const auto *osg1_55 = buffer.data(osg1 + 55);
    const auto *osg1_75 = buffer.data(osg1 + 75);
    const auto *osg1_78 = buffer.data(osg1 + 78);
    const auto *osg1_80 = buffer.data(osg1 + 80);

    const auto *qsd0_0 = buffer.data(qsd0 + 0);
    const auto *qsd0_3 = buffer.data(qsd0 + 3);
    const auto *qsd0_5 = buffer.data(qsd0 + 5);
    const auto *qsd0_9 = buffer.data(qsd0 + 9);
    const auto *qsd0_16 = buffer.data(qsd0 + 16);
    const auto *qsd0_17 = buffer.data(qsd0 + 17);
    const auto *qsd0_18 = buffer.data(qsd0 + 18);
    const auto *qsd0_21 = buffer.data(qsd0 + 21);
    const auto *qsd0_23 = buffer.data(qsd0 + 23);
    const auto *qsd0_29 = buffer.data(qsd0 + 29);
    const auto *qsd0_30 = buffer.data(qsd0 + 30);
    const auto *qsd0_33 = buffer.data(qsd0 + 33);
    const auto *qsd0_34 = buffer.data(qsd0 + 34);
    const auto *qsd0_35 = buffer.data(qsd0 + 35);
    const auto *qsd0_36 = buffer.data(qsd0 + 36);
    const auto *qsd0_39 = buffer.data(qsd0 + 39);
    const auto *qsd0_41 = buffer.data(qsd0 + 41);
    const auto *qsd0_47 = buffer.data(qsd0 + 47);

    const auto *qsd1_0 = buffer.data(qsd1 + 0);
    const auto *qsd1_3 = buffer.data(qsd1 + 3);
    const auto *qsd1_5 = buffer.data(qsd1 + 5);
    const auto *qsd1_9 = buffer.data(qsd1 + 9);
    const auto *qsd1_16 = buffer.data(qsd1 + 16);
    const auto *qsd1_17 = buffer.data(qsd1 + 17);
    const auto *qsd1_18 = buffer.data(qsd1 + 18);
    const auto *qsd1_21 = buffer.data(qsd1 + 21);
    const auto *qsd1_23 = buffer.data(qsd1 + 23);
    const auto *qsd1_29 = buffer.data(qsd1 + 29);
    const auto *qsd1_30 = buffer.data(qsd1 + 30);
    const auto *qsd1_33 = buffer.data(qsd1 + 33);
    const auto *qsd1_34 = buffer.data(qsd1 + 34);
    const auto *qsd1_35 = buffer.data(qsd1 + 35);
    const auto *qsd1_36 = buffer.data(qsd1 + 36);
    const auto *qsd1_39 = buffer.data(qsd1 + 39);
    const auto *qsd1_41 = buffer.data(qsd1 + 41);
    const auto *qsd1_47 = buffer.data(qsd1 + 47);

    const auto *qsf_0 = buffer.data(qsf + 0);
    const auto *qsf_1 = buffer.data(qsf + 1);
    const auto *qsf_2 = buffer.data(qsf + 2);
    const auto *qsf_3 = buffer.data(qsf + 3);
    const auto *qsf_5 = buffer.data(qsf + 5);
    const auto *qsf_6 = buffer.data(qsf + 6);
    const auto *qsf_8 = buffer.data(qsf + 8);
    const auto *qsf_9 = buffer.data(qsf + 9);
    const auto *qsf_10 = buffer.data(qsf + 10);
    const auto *qsf_11 = buffer.data(qsf + 11);
    const auto *qsf_13 = buffer.data(qsf + 13);
    const auto *qsf_16 = buffer.data(qsf + 16);
    const auto *qsf_17 = buffer.data(qsf + 17);
    const auto *qsf_18 = buffer.data(qsf + 18);
    const auto *qsf_19 = buffer.data(qsf + 19);
    const auto *qsf_20 = buffer.data(qsf + 20);
    const auto *qsf_22 = buffer.data(qsf + 22);
    const auto *qsf_25 = buffer.data(qsf + 25);
    const auto *qsf_26 = buffer.data(qsf + 26);
    const auto *qsf_27 = buffer.data(qsf + 27);
    const auto *qsf_28 = buffer.data(qsf + 28);
    const auto *qsf_29 = buffer.data(qsf + 29);
    const auto *qsf_30 = buffer.data(qsf + 30);
    const auto *qsf_31 = buffer.data(qsf + 31);
    const auto *qsf_32 = buffer.data(qsf + 32);
    const auto *qsf_33 = buffer.data(qsf + 33);
    const auto *qsf_36 = buffer.data(qsf + 36);
    const auto *qsf_37 = buffer.data(qsf + 37);
    const auto *qsf_38 = buffer.data(qsf + 38);
    const auto *qsf_39 = buffer.data(qsf + 39);
    const auto *qsf_40 = buffer.data(qsf + 40);
    const auto *qsf_42 = buffer.data(qsf + 42);
    const auto *qsf_46 = buffer.data(qsf + 46);
    const auto *qsf_47 = buffer.data(qsf + 47);
    const auto *qsf_48 = buffer.data(qsf + 48);
    const auto *qsf_49 = buffer.data(qsf + 49);
    const auto *qsf_50 = buffer.data(qsf + 50);
    const auto *qsf_51 = buffer.data(qsf + 51);
    const auto *qsf_52 = buffer.data(qsf + 52);
    const auto *qsf_55 = buffer.data(qsf + 55);
    const auto *qsf_56 = buffer.data(qsf + 56);
    const auto *qsf_57 = buffer.data(qsf + 57);
    const auto *qsf_58 = buffer.data(qsf + 58);
    const auto *qsf_59 = buffer.data(qsf + 59);
    const auto *qsf_60 = buffer.data(qsf + 60);
    const auto *qsf_61 = buffer.data(qsf + 61);
    const auto *qsf_62 = buffer.data(qsf + 62);
    const auto *qsf_63 = buffer.data(qsf + 63);
    const auto *qsf_66 = buffer.data(qsf + 66);
    const auto *qsf_67 = buffer.data(qsf + 67);
    const auto *qsf_68 = buffer.data(qsf + 68);
    const auto *qsf_69 = buffer.data(qsf + 69);
    const auto *qsf_70 = buffer.data(qsf + 70);
    const auto *qsf_72 = buffer.data(qsf + 72);
    const auto *qsf_75 = buffer.data(qsf + 75);
    const auto *qsf_76 = buffer.data(qsf + 76);
    const auto *qsf_77 = buffer.data(qsf + 77);
    const auto *qsf_78 = buffer.data(qsf + 78);
    const auto *qsf_79 = buffer.data(qsf + 79);
    const auto *qsf_80 = buffer.data(qsf + 80);
    const auto *qsf_82 = buffer.data(qsf + 82);
    const auto *qsf_86 = buffer.data(qsf + 86);
    const auto *qsf_87 = buffer.data(qsf + 87);
    const auto *qsf_88 = buffer.data(qsf + 88);
    const auto *qsf_89 = buffer.data(qsf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, osf_0, qsd0_0, \
                         qsd1_0, qsf_0, qsf_1, qsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * osf_0[k]
                 + f_1 * qsd0_0[k]
                 - f_2 * qsd1_0[k]
                 + f_3 * pc_x[k] * qsf_0[k];

        t_1[k] = f_3 * pc_y[k] * qsf_0[k];

        t_2[k] = f_3 * pc_z[k] * qsf_0[k];

        t_3[k] = f_4 * qsd0_0[k]
                 - f_5 * qsd1_0[k]
                 + f_3 * pc_y[k] * qsf_1[k];

        t_4[k] = f_3 * pc_y[k] * qsf_2[k];

        t_5[k] = f_4 * qsd0_0[k]
                 - f_5 * qsd1_0[k]
                 + f_3 * pc_z[k] * qsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, osf_6, osf_9, qsd0_3, \
                         qsd1_3, qsf_3, qsf_5, qsf_6, qsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * osf_6[k]
                 + f_3 * pc_x[k] * qsf_6[k];

        t_7[k] = f_3 * pc_z[k] * qsf_3[k];

        t_8[k] = f_3 * pc_y[k] * qsf_5[k];

        t_9[k] = f_0 * osf_9[k]
                 + f_3 * pc_x[k] * qsf_9[k];

        t_10[k] = f_1 * qsd0_3[k]
                  - f_2 * qsd1_3[k]
                  + f_3 * pc_y[k] * qsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, osg0_0, osg1_0, \
                         qsd0_5, qsd1_5, qsf_6, qsf_8, qsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * qsf_6[k];

        t_12[k] = f_4 * qsd0_5[k]
                  - f_5 * qsd1_5[k]
                  + f_3 * pc_y[k] * qsf_8[k];

        t_13[k] = f_3 * pc_y[k] * qsf_9[k];

        t_14[k] = f_1 * qsd0_5[k]
                  - f_2 * qsd1_5[k]
                  + f_3 * pc_z[k] * qsf_9[k];

        t_15[k] = pa_y[k] * osg0_0[k]
                  - f_6 * pc_y[k] * osg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, osg0_3, osg0_5, \
                         osf_0, osf_1, osg1_3, osg1_5, qsf_10, qsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * osf_0[k]
                  + f_3 * pc_y[k] * qsf_10[k];

        t_17[k] = f_3 * pc_z[k] * qsf_10[k];

        t_18[k] = pa_y[k] * osg0_3[k]
                  + f_8 * osf_1[k]
                  - f_6 * pc_y[k] * osg1_3[k];

        t_19[k] = f_3 * pc_z[k] * qsf_11[k];

        t_20[k] = pa_y[k] * osg0_5[k]
                  - f_6 * pc_y[k] * osg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, osf_16, osf_18, osf_19, qsf_13, \
                         qsf_16, qsf_18, qsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * osf_16[k]
                  + f_3 * pc_x[k] * qsf_16[k];

        t_22[k] = f_3 * pc_z[k] * qsf_13[k];

        t_23[k] = f_9 * osf_18[k]
                  + f_3 * pc_x[k] * qsf_18[k];

        t_24[k] = f_9 * osf_19[k]
                  + f_3 * pc_x[k] * qsf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, osf_6, osf_9, qsd0_9, qsd1_9, \
                         qsf_16, qsf_17, qsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * osf_6[k]
                  + f_1 * qsd0_9[k]
                  - f_2 * qsd1_9[k]
                  + f_3 * pc_y[k] * qsf_16[k];

        t_26[k] = f_3 * pc_z[k] * qsf_16[k];

        t_27[k] = f_4 * qsd0_9[k]
                  - f_5 * qsd1_9[k]
                  + f_3 * pc_z[k] * qsf_17[k];

        t_28[k] = f_7 * osf_9[k]
                  + f_3 * pc_y[k] * qsf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, osg0_0, osg0_14, \
                         osf_0, osg1_0, osg1_14, qsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * osg0_14[k]
                  - f_6 * pc_y[k] * osg1_14[k];

        t_30[k] = pa_z[k] * osg0_0[k]
                  - f_6 * pc_z[k] * osg1_0[k];

        t_31[k] = f_3 * pc_y[k] * qsf_20[k];

        t_32[k] = f_7 * osf_0[k]
                  + f_3 * pc_z[k] * qsf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, osg0_3, osg0_5, \
                         osf_2, osf_26, osg1_3, osg1_5, qsf_22, \
                         qsf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * osg0_3[k]
                  - f_6 * pc_z[k] * osg1_3[k];

        t_34[k] = f_3 * pc_y[k] * qsf_22[k];

        t_35[k] = pa_z[k] * osg0_5[k]
                  + f_8 * osf_2[k]
                  - f_6 * pc_z[k] * osg1_5[k];

        t_36[k] = f_9 * osf_26[k]
                  + f_3 * pc_x[k] * qsf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, osg0_10, osf_27, \
                         osf_29, osg1_10, qsf_25, qsf_27, qsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * osf_27[k]
                  + f_3 * pc_x[k] * qsf_27[k];

        t_38[k] = f_3 * pc_y[k] * qsf_25[k];

        t_39[k] = f_9 * osf_29[k]
                  + f_3 * pc_x[k] * qsf_29[k];

        t_40[k] = pa_z[k] * osg0_10[k]
                  - f_6 * pc_z[k] * osg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, osf_9, qsd0_16, qsd0_17, qsd1_16, \
                         qsd1_17, qsf_27, qsf_28, qsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * qsd0_16[k]
                  - f_11 * qsd1_16[k]
                  + f_3 * pc_y[k] * qsf_27[k];

        t_42[k] = f_4 * qsd0_17[k]
                  - f_5 * qsd1_17[k]
                  + f_3 * pc_y[k] * qsf_28[k];

        t_43[k] = f_3 * pc_y[k] * qsf_29[k];

        t_44[k] = f_7 * osf_9[k]
                  + f_1 * qsd0_17[k]
                  - f_2 * qsd1_17[k]
                  + f_3 * pc_z[k] * qsf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, osf_10, osf_30, osf_33, \
                         qsd0_18, qsd0_21, qsd1_18, qsd1_21, qsf_30, \
                         qsf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * osf_30[k]
                  + f_1 * qsd0_18[k]
                  - f_2 * qsd1_18[k]
                  + f_3 * pc_x[k] * qsf_30[k];

        t_46[k] = f_8 * osf_10[k]
                  + f_3 * pc_y[k] * qsf_30[k];

        t_47[k] = f_3 * pc_z[k] * qsf_30[k];

        t_48[k] = f_12 * osf_33[k]
                  + f_4 * qsd0_21[k]
                  - f_5 * qsd1_21[k]
                  + f_3 * pc_x[k] * qsf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, osf_36, osf_38, qsd0_18, \
                         qsd1_18, qsf_31, qsf_32, qsf_33, qsf_36, \
                         qsf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * qsf_31[k];

        t_50[k] = f_4 * qsd0_18[k]
                  - f_5 * qsd1_18[k]
                  + f_3 * pc_z[k] * qsf_32[k];

        t_51[k] = f_12 * osf_36[k]
                  + f_3 * pc_x[k] * qsf_36[k];

        t_52[k] = f_3 * pc_z[k] * qsf_33[k];

        t_53[k] = f_12 * osf_38[k]
                  + f_3 * pc_x[k] * qsf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, osf_16, osf_19, \
                         osf_39, qsd0_21, qsd1_21, qsf_36, qsf_37, \
                         qsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * osf_39[k]
                  + f_3 * pc_x[k] * qsf_39[k];

        t_55[k] = f_8 * osf_16[k]
                  + f_1 * qsd0_21[k]
                  - f_2 * qsd1_21[k]
                  + f_3 * pc_y[k] * qsf_36[k];

        t_56[k] = f_3 * pc_z[k] * qsf_36[k];

        t_57[k] = f_4 * qsd0_21[k]
                  - f_5 * qsd1_21[k]
                  + f_3 * pc_z[k] * qsf_37[k];

        t_58[k] = f_8 * osf_19[k]
                  + f_3 * pc_y[k] * qsf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, osg0_30, osf_10, osf_20, \
                         osg1_30, qsd0_23, qsd1_23, qsf_39, qsf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * qsd0_23[k]
                  - f_2 * qsd1_23[k]
                  + f_3 * pc_z[k] * qsf_39[k];

        t_60[k] = pa_y[k] * osg0_30[k]
                  - f_6 * pc_y[k] * osg1_30[k];

        t_61[k] = f_7 * osf_20[k]
                  + f_3 * pc_y[k] * qsf_40[k];

        t_62[k] = f_7 * osf_10[k]
                  + f_3 * pc_z[k] * qsf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, osg0_18, osg0_35, osf_22, \
                         osg1_18, osg1_35, qsf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * osg0_18[k]
                  - f_6 * pc_z[k] * osg1_18[k];

        t_64[k] = f_7 * osf_22[k]
                  + f_3 * pc_y[k] * qsf_42[k];

        t_65[k] = pa_y[k] * osg0_35[k]
                  - f_6 * pc_y[k] * osg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, osf_46, osf_47, osf_48, osf_49, qsf_46, \
                         qsf_47, qsf_48, qsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * osf_46[k]
                  + f_3 * pc_x[k] * qsf_46[k];

        t_67[k] = f_12 * osf_47[k]
                  + f_3 * pc_x[k] * qsf_47[k];

        t_68[k] = f_12 * osf_48[k]
                  + f_3 * pc_x[k] * qsf_48[k];

        t_69[k] = f_12 * osf_49[k]
                  + f_3 * pc_x[k] * qsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, osg0_25, osf_16, osf_28, osg1_25, \
                         qsd0_29, qsd1_29, qsf_46, qsf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * osg0_25[k]
                  - f_6 * pc_z[k] * osg1_25[k];

        t_71[k] = f_7 * osf_16[k]
                  + f_3 * pc_z[k] * qsf_46[k];

        t_72[k] = f_7 * osf_28[k]
                  + f_4 * qsd0_29[k]
                  - f_5 * qsd1_29[k]
                  + f_3 * pc_y[k] * qsf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, osg0_44, osf_29, osf_50, \
                         osg1_44, qsd0_30, qsd1_30, qsf_49, qsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * osf_29[k]
                  + f_3 * pc_y[k] * qsf_49[k];

        t_74[k] = pa_y[k] * osg0_44[k]
                  - f_6 * pc_y[k] * osg1_44[k];

        t_75[k] = f_12 * osf_50[k]
                  + f_1 * qsd0_30[k]
                  - f_2 * qsd1_30[k]
                  + f_3 * pc_x[k] * qsf_50[k];

        t_76[k] = f_3 * pc_y[k] * qsf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, osf_20, qsd0_30, qsd1_30, qsf_50, \
                         qsf_51, qsf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * osf_20[k]
                  + f_3 * pc_z[k] * qsf_50[k];

        t_78[k] = f_4 * qsd0_30[k]
                  - f_5 * qsd1_30[k]
                  + f_3 * pc_y[k] * qsf_51[k];

        t_79[k] = f_3 * pc_y[k] * qsf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, osf_55, osf_56, osf_57, qsd0_35, \
                         qsd1_35, qsf_55, qsf_56, qsf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * osf_55[k]
                  + f_4 * qsd0_35[k]
                  - f_5 * qsd1_35[k]
                  + f_3 * pc_x[k] * qsf_55[k];

        t_81[k] = f_12 * osf_56[k]
                  + f_3 * pc_x[k] * qsf_56[k];

        t_82[k] = f_12 * osf_57[k]
                  + f_3 * pc_x[k] * qsf_57[k];

        t_83[k] = f_3 * pc_y[k] * qsf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, osf_59, qsd0_33, qsd0_34, qsd1_33, \
                         qsd1_34, qsf_56, qsf_57, qsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * osf_59[k]
                  + f_3 * pc_x[k] * qsf_59[k];

        t_85[k] = f_1 * qsd0_33[k]
                  - f_2 * qsd1_33[k]
                  + f_3 * pc_y[k] * qsf_56[k];

        t_86[k] = f_10 * qsd0_34[k]
                  - f_11 * qsd1_34[k]
                  + f_3 * pc_y[k] * qsf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, osf_29, osf_60, qsd0_35, \
                         qsd0_36, qsd1_35, qsd1_36, qsf_58, qsf_59, \
                         qsf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * qsd0_35[k]
                  - f_5 * qsd1_35[k]
                  + f_3 * pc_y[k] * qsf_58[k];

        t_88[k] = f_3 * pc_y[k] * qsf_59[k];

        t_89[k] = f_8 * osf_29[k]
                  + f_1 * qsd0_35[k]
                  - f_2 * qsd1_35[k]
                  + f_3 * pc_z[k] * qsf_59[k];

        t_90[k] = f_13 * osf_60[k]
                  + f_1 * qsd0_36[k]
                  - f_2 * qsd1_36[k]
                  + f_3 * pc_x[k] * qsf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, osf_30, osf_63, qsd0_39, \
                         qsd1_39, qsf_60, qsf_61, qsf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * osf_30[k]
                  + f_3 * pc_y[k] * qsf_60[k];

        t_92[k] = f_3 * pc_z[k] * qsf_60[k];

        t_93[k] = f_13 * osf_63[k]
                  + f_4 * qsd0_39[k]
                  - f_5 * qsd1_39[k]
                  + f_3 * pc_x[k] * qsf_63[k];

        t_94[k] = f_3 * pc_z[k] * qsf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, osf_66, osf_68, qsd0_36, qsd1_36, \
                         qsf_62, qsf_63, qsf_66, qsf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * qsd0_36[k]
                  - f_5 * qsd1_36[k]
                  + f_3 * pc_z[k] * qsf_62[k];

        t_96[k] = f_13 * osf_66[k]
                  + f_3 * pc_x[k] * qsf_66[k];

        t_97[k] = f_3 * pc_z[k] * qsf_63[k];

        t_98[k] = f_13 * osf_68[k]
                  + f_3 * pc_x[k] * qsf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, osf_36, osf_39, \
                         osf_69, qsd0_39, qsd1_39, qsf_66, qsf_67, \
                         qsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * osf_69[k]
                  + f_3 * pc_x[k] * qsf_69[k];

        t_100[k] = f_14 * osf_36[k]
                   + f_1 * qsd0_39[k]
                   - f_2 * qsd1_39[k]
                   + f_3 * pc_y[k] * qsf_66[k];

        t_101[k] = f_3 * pc_z[k] * qsf_66[k];

        t_102[k] = f_4 * qsd0_39[k]
                   - f_5 * qsd1_39[k]
                   + f_3 * pc_z[k] * qsf_67[k];

        t_103[k] = f_14 * osf_39[k]
                   + f_3 * pc_y[k] * qsf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, osg0_45, osf_30, \
                         osf_40, osg1_45, qsd0_41, qsd1_41, qsf_69, \
                         qsf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * qsd0_41[k]
                   - f_2 * qsd1_41[k]
                   + f_3 * pc_z[k] * qsf_69[k];

        t_105[k] = pa_z[k] * osg0_45[k]
                   - f_6 * pc_z[k] * osg1_45[k];

        t_106[k] = f_8 * osf_40[k]
                   + f_3 * pc_y[k] * qsf_70[k];

        t_107[k] = f_7 * osf_30[k]
                   + f_3 * pc_z[k] * qsf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, osg0_48, osf_42, osf_75, \
                         osg1_48, qsd0_47, qsd1_47, qsf_72, qsf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * osg0_48[k]
                   - f_6 * pc_z[k] * osg1_48[k];

        t_109[k] = f_8 * osf_42[k]
                   + f_3 * pc_y[k] * qsf_72[k];

        t_110[k] = f_13 * osf_75[k]
                   + f_4 * qsd0_47[k]
                   - f_5 * qsd1_47[k]
                   + f_3 * pc_x[k] * qsf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, osf_76, osf_77, osf_78, osf_79, \
                         qsf_76, qsf_77, qsf_78, qsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * osf_76[k]
                   + f_3 * pc_x[k] * qsf_76[k];

        t_112[k] = f_13 * osf_77[k]
                   + f_3 * pc_x[k] * qsf_77[k];

        t_113[k] = f_13 * osf_78[k]
                   + f_3 * pc_x[k] * qsf_78[k];

        t_114[k] = f_13 * osf_79[k]
                   + f_3 * pc_x[k] * qsf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, osg0_55, osf_36, osf_48, \
                         osg1_55, qsd0_47, qsd1_47, qsf_76, qsf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * osg0_55[k]
                   - f_6 * pc_z[k] * osg1_55[k];

        t_116[k] = f_7 * osf_36[k]
                   + f_3 * pc_z[k] * qsf_76[k];

        t_117[k] = f_8 * osf_48[k]
                   + f_4 * qsd0_47[k]
                   - f_5 * qsd1_47[k]
                   + f_3 * pc_y[k] * qsf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, osg0_75, osf_39, \
                         osf_49, osf_50, osg1_75, qsd0_47, qsd1_47, qsf_79, \
                         qsf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * osf_49[k]
                   + f_3 * pc_y[k] * qsf_79[k];

        t_119[k] = f_7 * osf_39[k]
                   + f_1 * qsd0_47[k]
                   - f_2 * qsd1_47[k]
                   + f_3 * pc_z[k] * qsf_79[k];

        t_120[k] = pa_y[k] * osg0_75[k]
                   - f_6 * pc_y[k] * osg1_75[k];

        t_121[k] = f_7 * osf_50[k]
                   + f_3 * pc_y[k] * qsf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, osg0_78, osg0_80, \
                         osf_40, osf_51, osf_52, osg1_78, osg1_80, qsf_80, \
                         qsf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * osf_40[k]
                   + f_3 * pc_z[k] * qsf_80[k];

        t_123[k] = pa_y[k] * osg0_78[k]
                   + f_8 * osf_51[k]
                   - f_6 * pc_y[k] * osg1_78[k];

        t_124[k] = f_7 * osf_52[k]
                   + f_3 * pc_y[k] * qsf_82[k];

        t_125[k] = pa_y[k] * osg0_80[k]
                   - f_6 * pc_y[k] * osg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, osf_86, osf_87, osf_88, osf_89, \
                         qsf_86, qsf_87, qsf_88, qsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * osf_86[k]
                   + f_3 * pc_x[k] * qsf_86[k];

        t_127[k] = f_13 * osf_87[k]
                   + f_3 * pc_x[k] * qsf_87[k];

        t_128[k] = f_13 * osf_88[k]
                   + f_3 * pc_x[k] * qsf_88[k];

        t_129[k] = f_13 * osf_89[k]
                   + f_3 * pc_x[k] * qsf_89[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_89 = buffer.data(osg0 + 89);
    const auto *osg0_90 = buffer.data(osg0 + 90);
    const auto *osg0_93 = buffer.data(osg0 + 93);
    const auto *osg0_100 = buffer.data(osg0 + 100);
    const auto *osg0_135 = buffer.data(osg0 + 135);
    const auto *osg0_138 = buffer.data(osg0 + 138);
    const auto *osg0_140 = buffer.data(osg0 + 140);
    const auto *osg0_149 = buffer.data(osg0 + 149);
    const auto *osg0_150 = buffer.data(osg0 + 150);
    const auto *osg0_153 = buffer.data(osg0 + 153);
    const auto *osg0_160 = buffer.data(osg0 + 160);

    const auto *osf_46 = buffer.data(osf + 46);
    const auto *osf_50 = buffer.data(osf + 50);
    const auto *osf_56 = buffer.data(osf + 56);
    const auto *osf_58 = buffer.data(osf + 58);
    const auto *osf_59 = buffer.data(osf + 59);
    const auto *osf_60 = buffer.data(osf + 60);
    const auto *osf_66 = buffer.data(osf + 66);
    const auto *osf_69 = buffer.data(osf + 69);
    const auto *osf_70 = buffer.data(osf + 70);
    const auto *osf_72 = buffer.data(osf + 72);
    const auto *osf_76 = buffer.data(osf + 76);
    const auto *osf_78 = buffer.data(osf + 78);
    const auto *osf_79 = buffer.data(osf + 79);
    const auto *osf_80 = buffer.data(osf + 80);
    const auto *osf_82 = buffer.data(osf + 82);
    const auto *osf_86 = buffer.data(osf + 86);
    const auto *osf_88 = buffer.data(osf + 88);
    const auto *osf_89 = buffer.data(osf + 89);
    const auto *osf_90 = buffer.data(osf + 90);
    const auto *osf_91 = buffer.data(osf + 91);
    const auto *osf_92 = buffer.data(osf + 92);
    const auto *osf_95 = buffer.data(osf + 95);
    const auto *osf_96 = buffer.data(osf + 96);
    const auto *osf_97 = buffer.data(osf + 97);
    const auto *osf_98 = buffer.data(osf + 98);
    const auto *osf_99 = buffer.data(osf + 99);
    const auto *osf_100 = buffer.data(osf + 100);
    const auto *osf_103 = buffer.data(osf + 103);
    const auto *osf_106 = buffer.data(osf + 106);
    const auto *osf_108 = buffer.data(osf + 108);
    const auto *osf_109 = buffer.data(osf + 109);
    const auto *osf_110 = buffer.data(osf + 110);
    const auto *osf_112 = buffer.data(osf + 112);
    const auto *osf_115 = buffer.data(osf + 115);
    const auto *osf_116 = buffer.data(osf + 116);
    const auto *osf_117 = buffer.data(osf + 117);
    const auto *osf_118 = buffer.data(osf + 118);
    const auto *osf_119 = buffer.data(osf + 119);
    const auto *osf_120 = buffer.data(osf + 120);
    const auto *osf_123 = buffer.data(osf + 123);
    const auto *osf_125 = buffer.data(osf + 125);
    const auto *osf_126 = buffer.data(osf + 126);
    const auto *osf_127 = buffer.data(osf + 127);
    const auto *osf_128 = buffer.data(osf + 128);
    const auto *osf_129 = buffer.data(osf + 129);
    const auto *osf_136 = buffer.data(osf + 136);
    const auto *osf_137 = buffer.data(osf + 137);
    const auto *osf_138 = buffer.data(osf + 138);
    const auto *osf_139 = buffer.data(osf + 139);
    const auto *osf_140 = buffer.data(osf + 140);
    const auto *osf_145 = buffer.data(osf + 145);
    const auto *osf_146 = buffer.data(osf + 146);
    const auto *osf_147 = buffer.data(osf + 147);
    const auto *osf_149 = buffer.data(osf + 149);
    const auto *osf_150 = buffer.data(osf + 150);
    const auto *osf_153 = buffer.data(osf + 153);
    const auto *osf_156 = buffer.data(osf + 156);
    const auto *osf_158 = buffer.data(osf + 158);
    const auto *osf_159 = buffer.data(osf + 159);
    const auto *osf_165 = buffer.data(osf + 165);
    const auto *osf_166 = buffer.data(osf + 166);
    const auto *osf_167 = buffer.data(osf + 167);
    const auto *osf_168 = buffer.data(osf + 168);
    const auto *osf_169 = buffer.data(osf + 169);

    const auto *osg1_89 = buffer.data(osg1 + 89);
    const auto *osg1_90 = buffer.data(osg1 + 90);
    const auto *osg1_93 = buffer.data(osg1 + 93);
    const auto *osg1_100 = buffer.data(osg1 + 100);
    const auto *osg1_135 = buffer.data(osg1 + 135);
    const auto *osg1_138 = buffer.data(osg1 + 138);
    const auto *osg1_140 = buffer.data(osg1 + 140);
    const auto *osg1_149 = buffer.data(osg1 + 149);
    const auto *osg1_150 = buffer.data(osg1 + 150);
    const auto *osg1_153 = buffer.data(osg1 + 153);
    const auto *osg1_160 = buffer.data(osg1 + 160);

    const auto *qsd0_51 = buffer.data(qsd0 + 51);
    const auto *qsd0_53 = buffer.data(qsd0 + 53);
    const auto *qsd0_54 = buffer.data(qsd0 + 54);
    const auto *qsd0_57 = buffer.data(qsd0 + 57);
    const auto *qsd0_58 = buffer.data(qsd0 + 58);
    const auto *qsd0_59 = buffer.data(qsd0 + 59);
    const auto *qsd0_60 = buffer.data(qsd0 + 60);
    const auto *qsd0_63 = buffer.data(qsd0 + 63);
    const auto *qsd0_65 = buffer.data(qsd0 + 65);
    const auto *qsd0_71 = buffer.data(qsd0 + 71);
    const auto *qsd0_72 = buffer.data(qsd0 + 72);
    const auto *qsd0_75 = buffer.data(qsd0 + 75);
    const auto *qsd0_77 = buffer.data(qsd0 + 77);
    const auto *qsd0_81 = buffer.data(qsd0 + 81);
    const auto *qsd0_83 = buffer.data(qsd0 + 83);
    const auto *qsd0_84 = buffer.data(qsd0 + 84);
    const auto *qsd0_87 = buffer.data(qsd0 + 87);
    const auto *qsd0_88 = buffer.data(qsd0 + 88);
    const auto *qsd0_89 = buffer.data(qsd0 + 89);
    const auto *qsd0_90 = buffer.data(qsd0 + 90);
    const auto *qsd0_93 = buffer.data(qsd0 + 93);
    const auto *qsd0_95 = buffer.data(qsd0 + 95);
    const auto *qsd0_101 = buffer.data(qsd0 + 101);

    const auto *qsd1_51 = buffer.data(qsd1 + 51);
    const auto *qsd1_53 = buffer.data(qsd1 + 53);
    const auto *qsd1_54 = buffer.data(qsd1 + 54);
    const auto *qsd1_57 = buffer.data(qsd1 + 57);
    const auto *qsd1_58 = buffer.data(qsd1 + 58);
    const auto *qsd1_59 = buffer.data(qsd1 + 59);
    const auto *qsd1_60 = buffer.data(qsd1 + 60);
    const auto *qsd1_63 = buffer.data(qsd1 + 63);
    const auto *qsd1_65 = buffer.data(qsd1 + 65);
    const auto *qsd1_71 = buffer.data(qsd1 + 71);
    const auto *qsd1_72 = buffer.data(qsd1 + 72);
    const auto *qsd1_75 = buffer.data(qsd1 + 75);
    const auto *qsd1_77 = buffer.data(qsd1 + 77);
    const auto *qsd1_81 = buffer.data(qsd1 + 81);
    const auto *qsd1_83 = buffer.data(qsd1 + 83);
    const auto *qsd1_84 = buffer.data(qsd1 + 84);
    const auto *qsd1_87 = buffer.data(qsd1 + 87);
    const auto *qsd1_88 = buffer.data(qsd1 + 88);
    const auto *qsd1_89 = buffer.data(qsd1 + 89);
    const auto *qsd1_90 = buffer.data(qsd1 + 90);
    const auto *qsd1_93 = buffer.data(qsd1 + 93);
    const auto *qsd1_95 = buffer.data(qsd1 + 95);
    const auto *qsd1_101 = buffer.data(qsd1 + 101);

    const auto *qsf_86 = buffer.data(qsf + 86);
    const auto *qsf_88 = buffer.data(qsf + 88);
    const auto *qsf_89 = buffer.data(qsf + 89);
    const auto *qsf_90 = buffer.data(qsf + 90);
    const auto *qsf_91 = buffer.data(qsf + 91);
    const auto *qsf_92 = buffer.data(qsf + 92);
    const auto *qsf_95 = buffer.data(qsf + 95);
    const auto *qsf_96 = buffer.data(qsf + 96);
    const auto *qsf_97 = buffer.data(qsf + 97);
    const auto *qsf_98 = buffer.data(qsf + 98);
    const auto *qsf_99 = buffer.data(qsf + 99);
    const auto *qsf_100 = buffer.data(qsf + 100);
    const auto *qsf_101 = buffer.data(qsf + 101);
    const auto *qsf_102 = buffer.data(qsf + 102);
    const auto *qsf_103 = buffer.data(qsf + 103);
    const auto *qsf_106 = buffer.data(qsf + 106);
    const auto *qsf_107 = buffer.data(qsf + 107);
    const auto *qsf_108 = buffer.data(qsf + 108);
    const auto *qsf_109 = buffer.data(qsf + 109);
    const auto *qsf_110 = buffer.data(qsf + 110);
    const auto *qsf_112 = buffer.data(qsf + 112);
    const auto *qsf_115 = buffer.data(qsf + 115);
    const auto *qsf_116 = buffer.data(qsf + 116);
    const auto *qsf_117 = buffer.data(qsf + 117);
    const auto *qsf_118 = buffer.data(qsf + 118);
    const auto *qsf_119 = buffer.data(qsf + 119);
    const auto *qsf_120 = buffer.data(qsf + 120);
    const auto *qsf_122 = buffer.data(qsf + 122);
    const auto *qsf_123 = buffer.data(qsf + 123);
    const auto *qsf_125 = buffer.data(qsf + 125);
    const auto *qsf_126 = buffer.data(qsf + 126);
    const auto *qsf_127 = buffer.data(qsf + 127);
    const auto *qsf_128 = buffer.data(qsf + 128);
    const auto *qsf_129 = buffer.data(qsf + 129);
    const auto *qsf_130 = buffer.data(qsf + 130);
    const auto *qsf_132 = buffer.data(qsf + 132);
    const auto *qsf_136 = buffer.data(qsf + 136);
    const auto *qsf_137 = buffer.data(qsf + 137);
    const auto *qsf_138 = buffer.data(qsf + 138);
    const auto *qsf_139 = buffer.data(qsf + 139);
    const auto *qsf_140 = buffer.data(qsf + 140);
    const auto *qsf_141 = buffer.data(qsf + 141);
    const auto *qsf_142 = buffer.data(qsf + 142);
    const auto *qsf_145 = buffer.data(qsf + 145);
    const auto *qsf_146 = buffer.data(qsf + 146);
    const auto *qsf_147 = buffer.data(qsf + 147);
    const auto *qsf_148 = buffer.data(qsf + 148);
    const auto *qsf_149 = buffer.data(qsf + 149);
    const auto *qsf_150 = buffer.data(qsf + 150);
    const auto *qsf_151 = buffer.data(qsf + 151);
    const auto *qsf_152 = buffer.data(qsf + 152);
    const auto *qsf_153 = buffer.data(qsf + 153);
    const auto *qsf_156 = buffer.data(qsf + 156);
    const auto *qsf_157 = buffer.data(qsf + 157);
    const auto *qsf_158 = buffer.data(qsf + 158);
    const auto *qsf_159 = buffer.data(qsf + 159);
    const auto *qsf_160 = buffer.data(qsf + 160);
    const auto *qsf_162 = buffer.data(qsf + 162);
    const auto *qsf_165 = buffer.data(qsf + 165);
    const auto *qsf_166 = buffer.data(qsf + 166);
    const auto *qsf_167 = buffer.data(qsf + 167);
    const auto *qsf_168 = buffer.data(qsf + 168);
    const auto *qsf_169 = buffer.data(qsf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, osf_46, osf_56, osf_58, qsd0_51, \
                         qsd0_53, qsd1_51, qsd1_53, qsf_86, qsf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * osf_56[k]
                   + f_1 * qsd0_51[k]
                   - f_2 * qsd1_51[k]
                   + f_3 * pc_y[k] * qsf_86[k];

        t_131[k] = f_8 * osf_46[k]
                   + f_3 * pc_z[k] * qsf_86[k];

        t_132[k] = f_7 * osf_58[k]
                   + f_4 * qsd0_53[k]
                   - f_5 * qsd1_53[k]
                   + f_3 * pc_y[k] * qsf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, osg0_89, osf_59, \
                         osf_90, osg1_89, qsd0_54, qsd1_54, qsf_89, \
                         qsf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * osf_59[k]
                   + f_3 * pc_y[k] * qsf_89[k];

        t_134[k] = pa_y[k] * osg0_89[k]
                   - f_6 * pc_y[k] * osg1_89[k];

        t_135[k] = f_13 * osf_90[k]
                   + f_1 * qsd0_54[k]
                   - f_2 * qsd1_54[k]
                   + f_3 * pc_x[k] * qsf_90[k];

        t_136[k] = f_3 * pc_y[k] * qsf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, osf_50, qsd0_54, qsd1_54, qsf_90, \
                         qsf_91, qsf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * osf_50[k]
                   + f_3 * pc_z[k] * qsf_90[k];

        t_138[k] = f_4 * qsd0_54[k]
                   - f_5 * qsd1_54[k]
                   + f_3 * pc_y[k] * qsf_91[k];

        t_139[k] = f_3 * pc_y[k] * qsf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, osf_95, osf_96, osf_97, \
                         qsd0_59, qsd1_59, qsf_95, qsf_96, qsf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * osf_95[k]
                   + f_4 * qsd0_59[k]
                   - f_5 * qsd1_59[k]
                   + f_3 * pc_x[k] * qsf_95[k];

        t_141[k] = f_13 * osf_96[k]
                   + f_3 * pc_x[k] * qsf_96[k];

        t_142[k] = f_13 * osf_97[k]
                   + f_3 * pc_x[k] * qsf_97[k];

        t_143[k] = f_3 * pc_y[k] * qsf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, osf_99, qsd0_57, qsd0_58, qsd1_57, \
                         qsd1_58, qsf_96, qsf_97, qsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * osf_99[k]
                   + f_3 * pc_x[k] * qsf_99[k];

        t_145[k] = f_1 * qsd0_57[k]
                   - f_2 * qsd1_57[k]
                   + f_3 * pc_y[k] * qsf_96[k];

        t_146[k] = f_10 * qsd0_58[k]
                   - f_11 * qsd1_58[k]
                   + f_3 * pc_y[k] * qsf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, osf_59, osf_100, \
                         qsd0_59, qsd0_60, qsd1_59, qsd1_60, qsf_98, qsf_99, \
                         qsf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * qsd0_59[k]
                   - f_5 * qsd1_59[k]
                   + f_3 * pc_y[k] * qsf_98[k];

        t_148[k] = f_3 * pc_y[k] * qsf_99[k];

        t_149[k] = f_14 * osf_59[k]
                   + f_1 * qsd0_59[k]
                   - f_2 * qsd1_59[k]
                   + f_3 * pc_z[k] * qsf_99[k];

        t_150[k] = f_15 * osf_100[k]
                   + f_1 * qsd0_60[k]
                   - f_2 * qsd1_60[k]
                   + f_3 * pc_x[k] * qsf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, osf_60, osf_103, \
                         qsd0_63, qsd1_63, qsf_100, qsf_101, qsf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_16 * osf_60[k]
                   + f_3 * pc_y[k] * qsf_100[k];

        t_152[k] = f_3 * pc_z[k] * qsf_100[k];

        t_153[k] = f_15 * osf_103[k]
                   + f_4 * qsd0_63[k]
                   - f_5 * qsd1_63[k]
                   + f_3 * pc_x[k] * qsf_103[k];

        t_154[k] = f_3 * pc_z[k] * qsf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, osf_106, osf_108, qsd0_60, \
                         qsd1_60, qsf_102, qsf_103, qsf_106, qsf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * qsd0_60[k]
                   - f_5 * qsd1_60[k]
                   + f_3 * pc_z[k] * qsf_102[k];

        t_156[k] = f_15 * osf_106[k]
                   + f_3 * pc_x[k] * qsf_106[k];

        t_157[k] = f_3 * pc_z[k] * qsf_103[k];

        t_158[k] = f_15 * osf_108[k]
                   + f_3 * pc_x[k] * qsf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, osf_66, osf_69, \
                         osf_109, qsd0_63, qsd1_63, qsf_106, qsf_107, \
                         qsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_15 * osf_109[k]
                   + f_3 * pc_x[k] * qsf_109[k];

        t_160[k] = f_16 * osf_66[k]
                   + f_1 * qsd0_63[k]
                   - f_2 * qsd1_63[k]
                   + f_3 * pc_y[k] * qsf_106[k];

        t_161[k] = f_3 * pc_z[k] * qsf_106[k];

        t_162[k] = f_4 * qsd0_63[k]
                   - f_5 * qsd1_63[k]
                   + f_3 * pc_z[k] * qsf_107[k];

        t_163[k] = f_16 * osf_69[k]
                   + f_3 * pc_y[k] * qsf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, osg0_90, osf_60, \
                         osf_70, osg1_90, qsd0_65, qsd1_65, qsf_109, \
                         qsf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * qsd0_65[k]
                   - f_2 * qsd1_65[k]
                   + f_3 * pc_z[k] * qsf_109[k];

        t_165[k] = pa_z[k] * osg0_90[k]
                   - f_6 * pc_z[k] * osg1_90[k];

        t_166[k] = f_14 * osf_70[k]
                   + f_3 * pc_y[k] * qsf_110[k];

        t_167[k] = f_7 * osf_60[k]
                   + f_3 * pc_z[k] * qsf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, osg0_93, osf_72, \
                         osf_115, osg1_93, qsd0_71, qsd1_71, qsf_112, \
                         qsf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * osg0_93[k]
                   - f_6 * pc_z[k] * osg1_93[k];

        t_169[k] = f_14 * osf_72[k]
                   + f_3 * pc_y[k] * qsf_112[k];

        t_170[k] = f_15 * osf_115[k]
                   + f_4 * qsd0_71[k]
                   - f_5 * qsd1_71[k]
                   + f_3 * pc_x[k] * qsf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, osf_116, osf_117, osf_118, osf_119, \
                         qsf_116, qsf_117, qsf_118, qsf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_15 * osf_116[k]
                   + f_3 * pc_x[k] * qsf_116[k];

        t_172[k] = f_15 * osf_117[k]
                   + f_3 * pc_x[k] * qsf_117[k];

        t_173[k] = f_15 * osf_118[k]
                   + f_3 * pc_x[k] * qsf_118[k];

        t_174[k] = f_15 * osf_119[k]
                   + f_3 * pc_x[k] * qsf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, osg0_100, osf_66, osf_78, \
                         osg1_100, qsd0_71, qsd1_71, qsf_116, qsf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * osg0_100[k]
                   - f_6 * pc_z[k] * osg1_100[k];

        t_176[k] = f_7 * osf_66[k]
                   + f_3 * pc_z[k] * qsf_116[k];

        t_177[k] = f_14 * osf_78[k]
                   + f_4 * qsd0_71[k]
                   - f_5 * qsd1_71[k]
                   + f_3 * pc_y[k] * qsf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, osf_69, osf_79, osf_120, \
                         qsd0_71, qsd0_72, qsd1_71, qsd1_72, qsf_119, \
                         qsf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * osf_79[k]
                   + f_3 * pc_y[k] * qsf_119[k];

        t_179[k] = f_7 * osf_69[k]
                   + f_1 * qsd0_71[k]
                   - f_2 * qsd1_71[k]
                   + f_3 * pc_z[k] * qsf_119[k];

        t_180[k] = f_15 * osf_120[k]
                   + f_1 * qsd0_72[k]
                   - f_2 * qsd1_72[k]
                   + f_3 * pc_x[k] * qsf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, osf_70, osf_80, osf_82, \
                         osf_123, qsd0_75, qsd1_75, qsf_120, qsf_122, \
                         qsf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * osf_80[k]
                   + f_3 * pc_y[k] * qsf_120[k];

        t_182[k] = f_8 * osf_70[k]
                   + f_3 * pc_z[k] * qsf_120[k];

        t_183[k] = f_15 * osf_123[k]
                   + f_4 * qsd0_75[k]
                   - f_5 * qsd1_75[k]
                   + f_3 * pc_x[k] * qsf_123[k];

        t_184[k] = f_8 * osf_82[k]
                   + f_3 * pc_y[k] * qsf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, osf_125, osf_126, osf_127, osf_128, \
                         qsd0_77, qsd1_77, qsf_125, qsf_126, qsf_127, \
                         qsf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * osf_125[k]
                   + f_4 * qsd0_77[k]
                   - f_5 * qsd1_77[k]
                   + f_3 * pc_x[k] * qsf_125[k];

        t_186[k] = f_15 * osf_126[k]
                   + f_3 * pc_x[k] * qsf_126[k];

        t_187[k] = f_15 * osf_127[k]
                   + f_3 * pc_x[k] * qsf_127[k];

        t_188[k] = f_15 * osf_128[k]
                   + f_3 * pc_x[k] * qsf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, osf_76, osf_86, osf_129, \
                         qsd0_75, qsd1_75, qsf_126, qsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_15 * osf_129[k]
                   + f_3 * pc_x[k] * qsf_129[k];

        t_190[k] = f_8 * osf_86[k]
                   + f_1 * qsd0_75[k]
                   - f_2 * qsd1_75[k]
                   + f_3 * pc_y[k] * qsf_126[k];

        t_191[k] = f_8 * osf_76[k]
                   + f_3 * pc_z[k] * qsf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, osg0_135, osf_79, \
                         osf_88, osf_89, osg1_135, qsd0_77, qsd1_77, qsf_128, \
                         qsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * osf_88[k]
                   + f_4 * qsd0_77[k]
                   - f_5 * qsd1_77[k]
                   + f_3 * pc_y[k] * qsf_128[k];

        t_193[k] = f_8 * osf_89[k]
                   + f_3 * pc_y[k] * qsf_129[k];

        t_194[k] = f_8 * osf_79[k]
                   + f_1 * qsd0_77[k]
                   - f_2 * qsd1_77[k]
                   + f_3 * pc_z[k] * qsf_129[k];

        t_195[k] = pa_y[k] * osg0_135[k]
                   - f_6 * pc_y[k] * osg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, osg0_138, osf_80, \
                         osf_90, osf_91, osf_92, osg1_138, qsf_130, \
                         qsf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * osf_90[k]
                   + f_3 * pc_y[k] * qsf_130[k];

        t_197[k] = f_14 * osf_80[k]
                   + f_3 * pc_z[k] * qsf_130[k];

        t_198[k] = pa_y[k] * osg0_138[k]
                   + f_8 * osf_91[k]
                   - f_6 * pc_y[k] * osg1_138[k];

        t_199[k] = f_7 * osf_92[k]
                   + f_3 * pc_y[k] * qsf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, osg0_140, osf_136, \
                         osf_137, osf_138, osg1_140, qsf_136, qsf_137, \
                         qsf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * osg0_140[k]
                   - f_6 * pc_y[k] * osg1_140[k];

        t_201[k] = f_15 * osf_136[k]
                   + f_3 * pc_x[k] * qsf_136[k];

        t_202[k] = f_15 * osf_137[k]
                   + f_3 * pc_x[k] * qsf_137[k];

        t_203[k] = f_15 * osf_138[k]
                   + f_3 * pc_x[k] * qsf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, osf_86, osf_96, osf_139, \
                         qsd0_81, qsd1_81, qsf_136, qsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_15 * osf_139[k]
                   + f_3 * pc_x[k] * qsf_139[k];

        t_205[k] = f_7 * osf_96[k]
                   + f_1 * qsd0_81[k]
                   - f_2 * qsd1_81[k]
                   + f_3 * pc_y[k] * qsf_136[k];

        t_206[k] = f_14 * osf_86[k]
                   + f_3 * pc_z[k] * qsf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, osg0_149, osf_98, osf_99, osg1_149, \
                         qsd0_83, qsd1_83, qsf_138, qsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * osf_98[k]
                   + f_4 * qsd0_83[k]
                   - f_5 * qsd1_83[k]
                   + f_3 * pc_y[k] * qsf_138[k];

        t_208[k] = f_7 * osf_99[k]
                   + f_3 * pc_y[k] * qsf_139[k];

        t_209[k] = pa_y[k] * osg0_149[k]
                   - f_6 * pc_y[k] * osg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, osf_90, osf_140, \
                         qsd0_84, qsd1_84, qsf_140, qsf_141, qsf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * osf_140[k]
                   + f_1 * qsd0_84[k]
                   - f_2 * qsd1_84[k]
                   + f_3 * pc_x[k] * qsf_140[k];

        t_211[k] = f_3 * pc_y[k] * qsf_140[k];

        t_212[k] = f_16 * osf_90[k]
                   + f_3 * pc_z[k] * qsf_140[k];

        t_213[k] = f_4 * qsd0_84[k]
                   - f_5 * qsd1_84[k]
                   + f_3 * pc_y[k] * qsf_141[k];

        t_214[k] = f_3 * pc_y[k] * qsf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, osf_145, osf_146, osf_147, \
                         qsd0_89, qsd1_89, qsf_145, qsf_146, qsf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * osf_145[k]
                   + f_4 * qsd0_89[k]
                   - f_5 * qsd1_89[k]
                   + f_3 * pc_x[k] * qsf_145[k];

        t_216[k] = f_15 * osf_146[k]
                   + f_3 * pc_x[k] * qsf_146[k];

        t_217[k] = f_15 * osf_147[k]
                   + f_3 * pc_x[k] * qsf_147[k];

        t_218[k] = f_3 * pc_y[k] * qsf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, osf_149, qsd0_87, qsd0_88, qsd1_87, \
                         qsd1_88, qsf_146, qsf_147, qsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_15 * osf_149[k]
                   + f_3 * pc_x[k] * qsf_149[k];

        t_220[k] = f_1 * qsd0_87[k]
                   - f_2 * qsd1_87[k]
                   + f_3 * pc_y[k] * qsf_146[k];

        t_221[k] = f_10 * qsd0_88[k]
                   - f_11 * qsd1_88[k]
                   + f_3 * pc_y[k] * qsf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, pc_z, osf_99, osf_150, \
                         qsd0_89, qsd0_90, qsd1_89, qsd1_90, qsf_148, qsf_149, \
                         qsf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * qsd0_89[k]
                   - f_5 * qsd1_89[k]
                   + f_3 * pc_y[k] * qsf_148[k];

        t_223[k] = f_3 * pc_y[k] * qsf_149[k];

        t_224[k] = f_16 * osf_99[k]
                   + f_1 * qsd0_89[k]
                   - f_2 * qsd1_89[k]
                   + f_3 * pc_z[k] * qsf_149[k];

        t_225[k] = f_17 * osf_150[k]
                   + f_1 * qsd0_90[k]
                   - f_2 * qsd1_90[k]
                   + f_3 * pc_x[k] * qsf_150[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, osf_100, osf_153, \
                         qsd0_93, qsd1_93, qsf_150, qsf_151, qsf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_18 * osf_100[k]
                   + f_3 * pc_y[k] * qsf_150[k];

        t_227[k] = f_3 * pc_z[k] * qsf_150[k];

        t_228[k] = f_17 * osf_153[k]
                   + f_4 * qsd0_93[k]
                   - f_5 * qsd1_93[k]
                   + f_3 * pc_x[k] * qsf_153[k];

        t_229[k] = f_3 * pc_z[k] * qsf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, osf_156, osf_158, qsd0_90, \
                         qsd1_90, qsf_152, qsf_153, qsf_156, qsf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * qsd0_90[k]
                   - f_5 * qsd1_90[k]
                   + f_3 * pc_z[k] * qsf_152[k];

        t_231[k] = f_17 * osf_156[k]
                   + f_3 * pc_x[k] * qsf_156[k];

        t_232[k] = f_3 * pc_z[k] * qsf_153[k];

        t_233[k] = f_17 * osf_158[k]
                   + f_3 * pc_x[k] * qsf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, osf_106, \
                         osf_109, osf_159, qsd0_93, qsd1_93, qsf_156, qsf_157, \
                         qsf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_17 * osf_159[k]
                   + f_3 * pc_x[k] * qsf_159[k];

        t_235[k] = f_18 * osf_106[k]
                   + f_1 * qsd0_93[k]
                   - f_2 * qsd1_93[k]
                   + f_3 * pc_y[k] * qsf_156[k];

        t_236[k] = f_3 * pc_z[k] * qsf_156[k];

        t_237[k] = f_4 * qsd0_93[k]
                   - f_5 * qsd1_93[k]
                   + f_3 * pc_z[k] * qsf_157[k];

        t_238[k] = f_18 * osf_109[k]
                   + f_3 * pc_y[k] * qsf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pc_y, pc_z, osg0_150, osf_100, \
                         osf_110, osg1_150, qsd0_95, qsd1_95, qsf_159, \
                         qsf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * qsd0_95[k]
                   - f_2 * qsd1_95[k]
                   + f_3 * pc_z[k] * qsf_159[k];

        t_240[k] = pa_z[k] * osg0_150[k]
                   - f_6 * pc_z[k] * osg1_150[k];

        t_241[k] = f_16 * osf_110[k]
                   + f_3 * pc_y[k] * qsf_160[k];

        t_242[k] = f_7 * osf_100[k]
                   + f_3 * pc_z[k] * qsf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, osg0_153, osf_112, \
                         osf_165, osg1_153, qsd0_101, qsd1_101, qsf_162, \
                         qsf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * osg0_153[k]
                   - f_6 * pc_z[k] * osg1_153[k];

        t_244[k] = f_16 * osf_112[k]
                   + f_3 * pc_y[k] * qsf_162[k];

        t_245[k] = f_17 * osf_165[k]
                   + f_4 * qsd0_101[k]
                   - f_5 * qsd1_101[k]
                   + f_3 * pc_x[k] * qsf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, osf_166, osf_167, osf_168, osf_169, \
                         qsf_166, qsf_167, qsf_168, qsf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_17 * osf_166[k]
                   + f_3 * pc_x[k] * qsf_166[k];

        t_247[k] = f_17 * osf_167[k]
                   + f_3 * pc_x[k] * qsf_167[k];

        t_248[k] = f_17 * osf_168[k]
                   + f_3 * pc_x[k] * qsf_168[k];

        t_249[k] = f_17 * osf_169[k]
                   + f_3 * pc_x[k] * qsf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_z, pc_y, pc_z, osg0_160, osf_106, osf_118, \
                         osg1_160, qsd0_101, qsd1_101, qsf_166, \
                         qsf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * osg0_160[k]
                   - f_6 * pc_z[k] * osg1_160[k];

        t_251[k] = f_7 * osf_106[k]
                   + f_3 * pc_z[k] * qsf_166[k];

        t_252[k] = f_16 * osf_118[k]
                   + f_4 * qsd0_101[k]
                   - f_5 * qsd1_101[k]
                   + f_3 * pc_y[k] * qsf_168[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_210 = buffer.data(osg0 + 210);
    const auto *osg0_213 = buffer.data(osg0 + 213);
    const auto *osg0_215 = buffer.data(osg0 + 215);
    const auto *osg0_224 = buffer.data(osg0 + 224);
    const auto *osg0_225 = buffer.data(osg0 + 225);
    const auto *osg0_228 = buffer.data(osg0 + 228);
    const auto *osg0_235 = buffer.data(osg0 + 235);

    const auto *osf_109 = buffer.data(osf + 109);
    const auto *osf_110 = buffer.data(osf + 110);
    const auto *osf_116 = buffer.data(osf + 116);
    const auto *osf_119 = buffer.data(osf + 119);
    const auto *osf_120 = buffer.data(osf + 120);
    const auto *osf_122 = buffer.data(osf + 122);
    const auto *osf_126 = buffer.data(osf + 126);
    const auto *osf_128 = buffer.data(osf + 128);
    const auto *osf_129 = buffer.data(osf + 129);
    const auto *osf_130 = buffer.data(osf + 130);
    const auto *osf_132 = buffer.data(osf + 132);
    const auto *osf_136 = buffer.data(osf + 136);
    const auto *osf_138 = buffer.data(osf + 138);
    const auto *osf_139 = buffer.data(osf + 139);
    const auto *osf_140 = buffer.data(osf + 140);
    const auto *osf_141 = buffer.data(osf + 141);
    const auto *osf_142 = buffer.data(osf + 142);
    const auto *osf_146 = buffer.data(osf + 146);
    const auto *osf_148 = buffer.data(osf + 148);
    const auto *osf_149 = buffer.data(osf + 149);
    const auto *osf_150 = buffer.data(osf + 150);
    const auto *osf_156 = buffer.data(osf + 156);
    const auto *osf_159 = buffer.data(osf + 159);
    const auto *osf_160 = buffer.data(osf + 160);
    const auto *osf_162 = buffer.data(osf + 162);
    const auto *osf_166 = buffer.data(osf + 166);
    const auto *osf_168 = buffer.data(osf + 168);
    const auto *osf_169 = buffer.data(osf + 169);
    const auto *osf_170 = buffer.data(osf + 170);
    const auto *osf_172 = buffer.data(osf + 172);
    const auto *osf_173 = buffer.data(osf + 173);
    const auto *osf_175 = buffer.data(osf + 175);
    const auto *osf_176 = buffer.data(osf + 176);
    const auto *osf_177 = buffer.data(osf + 177);
    const auto *osf_178 = buffer.data(osf + 178);
    const auto *osf_179 = buffer.data(osf + 179);
    const auto *osf_180 = buffer.data(osf + 180);
    const auto *osf_182 = buffer.data(osf + 182);
    const auto *osf_183 = buffer.data(osf + 183);
    const auto *osf_185 = buffer.data(osf + 185);
    const auto *osf_186 = buffer.data(osf + 186);
    const auto *osf_187 = buffer.data(osf + 187);
    const auto *osf_188 = buffer.data(osf + 188);
    const auto *osf_189 = buffer.data(osf + 189);
    const auto *osf_196 = buffer.data(osf + 196);
    const auto *osf_197 = buffer.data(osf + 197);
    const auto *osf_198 = buffer.data(osf + 198);
    const auto *osf_199 = buffer.data(osf + 199);
    const auto *osf_200 = buffer.data(osf + 200);
    const auto *osf_205 = buffer.data(osf + 205);
    const auto *osf_206 = buffer.data(osf + 206);
    const auto *osf_207 = buffer.data(osf + 207);
    const auto *osf_209 = buffer.data(osf + 209);
    const auto *osf_210 = buffer.data(osf + 210);
    const auto *osf_213 = buffer.data(osf + 213);
    const auto *osf_216 = buffer.data(osf + 216);
    const auto *osf_218 = buffer.data(osf + 218);
    const auto *osf_219 = buffer.data(osf + 219);
    const auto *osf_225 = buffer.data(osf + 225);
    const auto *osf_226 = buffer.data(osf + 226);
    const auto *osf_227 = buffer.data(osf + 227);
    const auto *osf_228 = buffer.data(osf + 228);
    const auto *osf_229 = buffer.data(osf + 229);
    const auto *osf_230 = buffer.data(osf + 230);
    const auto *osf_233 = buffer.data(osf + 233);
    const auto *osf_235 = buffer.data(osf + 235);
    const auto *osf_236 = buffer.data(osf + 236);
    const auto *osf_237 = buffer.data(osf + 237);
    const auto *osf_238 = buffer.data(osf + 238);
    const auto *osf_239 = buffer.data(osf + 239);
    const auto *osf_240 = buffer.data(osf + 240);
    const auto *osf_243 = buffer.data(osf + 243);
    const auto *osf_245 = buffer.data(osf + 245);
    const auto *osf_246 = buffer.data(osf + 246);
    const auto *osf_247 = buffer.data(osf + 247);
    const auto *osf_248 = buffer.data(osf + 248);
    const auto *osf_249 = buffer.data(osf + 249);

    const auto *osg1_210 = buffer.data(osg1 + 210);
    const auto *osg1_213 = buffer.data(osg1 + 213);
    const auto *osg1_215 = buffer.data(osg1 + 215);
    const auto *osg1_224 = buffer.data(osg1 + 224);
    const auto *osg1_225 = buffer.data(osg1 + 225);
    const auto *osg1_228 = buffer.data(osg1 + 228);
    const auto *osg1_235 = buffer.data(osg1 + 235);

    const auto *qsd0_101 = buffer.data(qsd0 + 101);
    const auto *qsd0_102 = buffer.data(qsd0 + 102);
    const auto *qsd0_105 = buffer.data(qsd0 + 105);
    const auto *qsd0_107 = buffer.data(qsd0 + 107);
    const auto *qsd0_108 = buffer.data(qsd0 + 108);
    const auto *qsd0_111 = buffer.data(qsd0 + 111);
    const auto *qsd0_113 = buffer.data(qsd0 + 113);
    const auto *qsd0_117 = buffer.data(qsd0 + 117);
    const auto *qsd0_119 = buffer.data(qsd0 + 119);
    const auto *qsd0_120 = buffer.data(qsd0 + 120);
    const auto *qsd0_123 = buffer.data(qsd0 + 123);
    const auto *qsd0_124 = buffer.data(qsd0 + 124);
    const auto *qsd0_125 = buffer.data(qsd0 + 125);
    const auto *qsd0_126 = buffer.data(qsd0 + 126);
    const auto *qsd0_129 = buffer.data(qsd0 + 129);
    const auto *qsd0_131 = buffer.data(qsd0 + 131);
    const auto *qsd0_137 = buffer.data(qsd0 + 137);
    const auto *qsd0_138 = buffer.data(qsd0 + 138);
    const auto *qsd0_141 = buffer.data(qsd0 + 141);
    const auto *qsd0_143 = buffer.data(qsd0 + 143);
    const auto *qsd0_144 = buffer.data(qsd0 + 144);
    const auto *qsd0_147 = buffer.data(qsd0 + 147);
    const auto *qsd0_149 = buffer.data(qsd0 + 149);

    const auto *qsd1_101 = buffer.data(qsd1 + 101);
    const auto *qsd1_102 = buffer.data(qsd1 + 102);
    const auto *qsd1_105 = buffer.data(qsd1 + 105);
    const auto *qsd1_107 = buffer.data(qsd1 + 107);
    const auto *qsd1_108 = buffer.data(qsd1 + 108);
    const auto *qsd1_111 = buffer.data(qsd1 + 111);
    const auto *qsd1_113 = buffer.data(qsd1 + 113);
    const auto *qsd1_117 = buffer.data(qsd1 + 117);
    const auto *qsd1_119 = buffer.data(qsd1 + 119);
    const auto *qsd1_120 = buffer.data(qsd1 + 120);
    const auto *qsd1_123 = buffer.data(qsd1 + 123);
    const auto *qsd1_124 = buffer.data(qsd1 + 124);
    const auto *qsd1_125 = buffer.data(qsd1 + 125);
    const auto *qsd1_126 = buffer.data(qsd1 + 126);
    const auto *qsd1_129 = buffer.data(qsd1 + 129);
    const auto *qsd1_131 = buffer.data(qsd1 + 131);
    const auto *qsd1_137 = buffer.data(qsd1 + 137);
    const auto *qsd1_138 = buffer.data(qsd1 + 138);
    const auto *qsd1_141 = buffer.data(qsd1 + 141);
    const auto *qsd1_143 = buffer.data(qsd1 + 143);
    const auto *qsd1_144 = buffer.data(qsd1 + 144);
    const auto *qsd1_147 = buffer.data(qsd1 + 147);
    const auto *qsd1_149 = buffer.data(qsd1 + 149);

    const auto *qsf_169 = buffer.data(qsf + 169);
    const auto *qsf_170 = buffer.data(qsf + 170);
    const auto *qsf_172 = buffer.data(qsf + 172);
    const auto *qsf_173 = buffer.data(qsf + 173);
    const auto *qsf_175 = buffer.data(qsf + 175);
    const auto *qsf_176 = buffer.data(qsf + 176);
    const auto *qsf_177 = buffer.data(qsf + 177);
    const auto *qsf_178 = buffer.data(qsf + 178);
    const auto *qsf_179 = buffer.data(qsf + 179);
    const auto *qsf_180 = buffer.data(qsf + 180);
    const auto *qsf_182 = buffer.data(qsf + 182);
    const auto *qsf_183 = buffer.data(qsf + 183);
    const auto *qsf_185 = buffer.data(qsf + 185);
    const auto *qsf_186 = buffer.data(qsf + 186);
    const auto *qsf_187 = buffer.data(qsf + 187);
    const auto *qsf_188 = buffer.data(qsf + 188);
    const auto *qsf_189 = buffer.data(qsf + 189);
    const auto *qsf_190 = buffer.data(qsf + 190);
    const auto *qsf_192 = buffer.data(qsf + 192);
    const auto *qsf_196 = buffer.data(qsf + 196);
    const auto *qsf_197 = buffer.data(qsf + 197);
    const auto *qsf_198 = buffer.data(qsf + 198);
    const auto *qsf_199 = buffer.data(qsf + 199);
    const auto *qsf_200 = buffer.data(qsf + 200);
    const auto *qsf_201 = buffer.data(qsf + 201);
    const auto *qsf_202 = buffer.data(qsf + 202);
    const auto *qsf_205 = buffer.data(qsf + 205);
    const auto *qsf_206 = buffer.data(qsf + 206);
    const auto *qsf_207 = buffer.data(qsf + 207);
    const auto *qsf_208 = buffer.data(qsf + 208);
    const auto *qsf_209 = buffer.data(qsf + 209);
    const auto *qsf_210 = buffer.data(qsf + 210);
    const auto *qsf_211 = buffer.data(qsf + 211);
    const auto *qsf_212 = buffer.data(qsf + 212);
    const auto *qsf_213 = buffer.data(qsf + 213);
    const auto *qsf_216 = buffer.data(qsf + 216);
    const auto *qsf_217 = buffer.data(qsf + 217);
    const auto *qsf_218 = buffer.data(qsf + 218);
    const auto *qsf_219 = buffer.data(qsf + 219);
    const auto *qsf_220 = buffer.data(qsf + 220);
    const auto *qsf_222 = buffer.data(qsf + 222);
    const auto *qsf_225 = buffer.data(qsf + 225);
    const auto *qsf_226 = buffer.data(qsf + 226);
    const auto *qsf_227 = buffer.data(qsf + 227);
    const auto *qsf_228 = buffer.data(qsf + 228);
    const auto *qsf_229 = buffer.data(qsf + 229);
    const auto *qsf_230 = buffer.data(qsf + 230);
    const auto *qsf_232 = buffer.data(qsf + 232);
    const auto *qsf_233 = buffer.data(qsf + 233);
    const auto *qsf_235 = buffer.data(qsf + 235);
    const auto *qsf_236 = buffer.data(qsf + 236);
    const auto *qsf_237 = buffer.data(qsf + 237);
    const auto *qsf_238 = buffer.data(qsf + 238);
    const auto *qsf_239 = buffer.data(qsf + 239);
    const auto *qsf_240 = buffer.data(qsf + 240);
    const auto *qsf_242 = buffer.data(qsf + 242);
    const auto *qsf_243 = buffer.data(qsf + 243);
    const auto *qsf_245 = buffer.data(qsf + 245);
    const auto *qsf_246 = buffer.data(qsf + 246);
    const auto *qsf_247 = buffer.data(qsf + 247);
    const auto *qsf_248 = buffer.data(qsf + 248);
    const auto *qsf_249 = buffer.data(qsf + 249);

#pragma omp simd aligned(t_253, t_254, t_255, pc_x, pc_y, pc_z, osf_109, osf_119, osf_170, \
                         qsd0_101, qsd0_102, qsd1_101, qsd1_102, qsf_169, \
                         qsf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_16 * osf_119[k]
                   + f_3 * pc_y[k] * qsf_169[k];

        t_254[k] = f_7 * osf_109[k]
                   + f_1 * qsd0_101[k]
                   - f_2 * qsd1_101[k]
                   + f_3 * pc_z[k] * qsf_169[k];

        t_255[k] = f_17 * osf_170[k]
                   + f_1 * qsd0_102[k]
                   - f_2 * qsd1_102[k]
                   + f_3 * pc_x[k] * qsf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, osf_110, osf_120, \
                         osf_122, osf_173, qsd0_105, qsd1_105, qsf_170, qsf_172, \
                         qsf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * osf_120[k]
                   + f_3 * pc_y[k] * qsf_170[k];

        t_257[k] = f_8 * osf_110[k]
                   + f_3 * pc_z[k] * qsf_170[k];

        t_258[k] = f_17 * osf_173[k]
                   + f_4 * qsd0_105[k]
                   - f_5 * qsd1_105[k]
                   + f_3 * pc_x[k] * qsf_173[k];

        t_259[k] = f_14 * osf_122[k]
                   + f_3 * pc_y[k] * qsf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, osf_175, osf_176, osf_177, osf_178, \
                         qsd0_107, qsd1_107, qsf_175, qsf_176, qsf_177, \
                         qsf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_17 * osf_175[k]
                   + f_4 * qsd0_107[k]
                   - f_5 * qsd1_107[k]
                   + f_3 * pc_x[k] * qsf_175[k];

        t_261[k] = f_17 * osf_176[k]
                   + f_3 * pc_x[k] * qsf_176[k];

        t_262[k] = f_17 * osf_177[k]
                   + f_3 * pc_x[k] * qsf_177[k];

        t_263[k] = f_17 * osf_178[k]
                   + f_3 * pc_x[k] * qsf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, pc_z, osf_116, osf_126, osf_179, \
                         qsd0_105, qsd1_105, qsf_176, qsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * osf_179[k]
                   + f_3 * pc_x[k] * qsf_179[k];

        t_265[k] = f_14 * osf_126[k]
                   + f_1 * qsd0_105[k]
                   - f_2 * qsd1_105[k]
                   + f_3 * pc_y[k] * qsf_176[k];

        t_266[k] = f_8 * osf_116[k]
                   + f_3 * pc_z[k] * qsf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, osf_119, osf_128, osf_129, qsd0_107, \
                         qsd1_107, qsf_178, qsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * osf_128[k]
                   + f_4 * qsd0_107[k]
                   - f_5 * qsd1_107[k]
                   + f_3 * pc_y[k] * qsf_178[k];

        t_268[k] = f_14 * osf_129[k]
                   + f_3 * pc_y[k] * qsf_179[k];

        t_269[k] = f_8 * osf_119[k]
                   + f_1 * qsd0_107[k]
                   - f_2 * qsd1_107[k]
                   + f_3 * pc_z[k] * qsf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, osf_120, osf_130, osf_180, \
                         qsd0_108, qsd1_108, qsf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_17 * osf_180[k]
                   + f_1 * qsd0_108[k]
                   - f_2 * qsd1_108[k]
                   + f_3 * pc_x[k] * qsf_180[k];

        t_271[k] = f_8 * osf_130[k]
                   + f_3 * pc_y[k] * qsf_180[k];

        t_272[k] = f_14 * osf_120[k]
                   + f_3 * pc_z[k] * qsf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, osf_132, osf_183, osf_185, qsd0_111, \
                         qsd0_113, qsd1_111, qsd1_113, qsf_182, qsf_183, \
                         qsf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_17 * osf_183[k]
                   + f_4 * qsd0_111[k]
                   - f_5 * qsd1_111[k]
                   + f_3 * pc_x[k] * qsf_183[k];

        t_274[k] = f_8 * osf_132[k]
                   + f_3 * pc_y[k] * qsf_182[k];

        t_275[k] = f_17 * osf_185[k]
                   + f_4 * qsd0_113[k]
                   - f_5 * qsd1_113[k]
                   + f_3 * pc_x[k] * qsf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, osf_186, osf_187, osf_188, osf_189, \
                         qsf_186, qsf_187, qsf_188, qsf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_17 * osf_186[k]
                   + f_3 * pc_x[k] * qsf_186[k];

        t_277[k] = f_17 * osf_187[k]
                   + f_3 * pc_x[k] * qsf_187[k];

        t_278[k] = f_17 * osf_188[k]
                   + f_3 * pc_x[k] * qsf_188[k];

        t_279[k] = f_17 * osf_189[k]
                   + f_3 * pc_x[k] * qsf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_y, pc_z, osf_126, osf_136, osf_138, qsd0_111, \
                         qsd0_113, qsd1_111, qsd1_113, qsf_186, \
                         qsf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * osf_136[k]
                   + f_1 * qsd0_111[k]
                   - f_2 * qsd1_111[k]
                   + f_3 * pc_y[k] * qsf_186[k];

        t_281[k] = f_14 * osf_126[k]
                   + f_3 * pc_z[k] * qsf_186[k];

        t_282[k] = f_8 * osf_138[k]
                   + f_4 * qsd0_113[k]
                   - f_5 * qsd1_113[k]
                   + f_3 * pc_y[k] * qsf_188[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_y, pc_y, pc_z, osg0_210, osf_129, \
                         osf_139, osf_140, osg1_210, qsd0_113, qsd1_113, qsf_189, \
                         qsf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * osf_139[k]
                   + f_3 * pc_y[k] * qsf_189[k];

        t_284[k] = f_14 * osf_129[k]
                   + f_1 * qsd0_113[k]
                   - f_2 * qsd1_113[k]
                   + f_3 * pc_z[k] * qsf_189[k];

        t_285[k] = pa_y[k] * osg0_210[k]
                   - f_6 * pc_y[k] * osg1_210[k];

        t_286[k] = f_7 * osf_140[k]
                   + f_3 * pc_y[k] * qsf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_y, pc_y, pc_z, osg0_213, osg0_215, \
                         osf_130, osf_141, osf_142, osg1_213, osg1_215, qsf_190, \
                         qsf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_16 * osf_130[k]
                   + f_3 * pc_z[k] * qsf_190[k];

        t_288[k] = pa_y[k] * osg0_213[k]
                   + f_8 * osf_141[k]
                   - f_6 * pc_y[k] * osg1_213[k];

        t_289[k] = f_7 * osf_142[k]
                   + f_3 * pc_y[k] * qsf_192[k];

        t_290[k] = pa_y[k] * osg0_215[k]
                   - f_6 * pc_y[k] * osg1_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, osf_196, osf_197, osf_198, osf_199, \
                         qsf_196, qsf_197, qsf_198, qsf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_17 * osf_196[k]
                   + f_3 * pc_x[k] * qsf_196[k];

        t_292[k] = f_17 * osf_197[k]
                   + f_3 * pc_x[k] * qsf_197[k];

        t_293[k] = f_17 * osf_198[k]
                   + f_3 * pc_x[k] * qsf_198[k];

        t_294[k] = f_17 * osf_199[k]
                   + f_3 * pc_x[k] * qsf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pc_y, pc_z, osf_136, osf_146, osf_148, qsd0_117, \
                         qsd0_119, qsd1_117, qsd1_119, qsf_196, \
                         qsf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_7 * osf_146[k]
                   + f_1 * qsd0_117[k]
                   - f_2 * qsd1_117[k]
                   + f_3 * pc_y[k] * qsf_196[k];

        t_296[k] = f_16 * osf_136[k]
                   + f_3 * pc_z[k] * qsf_196[k];

        t_297[k] = f_7 * osf_148[k]
                   + f_4 * qsd0_119[k]
                   - f_5 * qsd1_119[k]
                   + f_3 * pc_y[k] * qsf_198[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pc_x, pc_y, osg0_224, osf_149, \
                         osf_200, osg1_224, qsd0_120, qsd1_120, qsf_199, \
                         qsf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * osf_149[k]
                   + f_3 * pc_y[k] * qsf_199[k];

        t_299[k] = pa_y[k] * osg0_224[k]
                   - f_6 * pc_y[k] * osg1_224[k];

        t_300[k] = f_17 * osf_200[k]
                   + f_1 * qsd0_120[k]
                   - f_2 * qsd1_120[k]
                   + f_3 * pc_x[k] * qsf_200[k];

        t_301[k] = f_3 * pc_y[k] * qsf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, osf_140, qsd0_120, qsd1_120, \
                         qsf_200, qsf_201, qsf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_18 * osf_140[k]
                   + f_3 * pc_z[k] * qsf_200[k];

        t_303[k] = f_4 * qsd0_120[k]
                   - f_5 * qsd1_120[k]
                   + f_3 * pc_y[k] * qsf_201[k];

        t_304[k] = f_3 * pc_y[k] * qsf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, osf_205, osf_206, osf_207, \
                         qsd0_125, qsd1_125, qsf_205, qsf_206, \
                         qsf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_17 * osf_205[k]
                   + f_4 * qsd0_125[k]
                   - f_5 * qsd1_125[k]
                   + f_3 * pc_x[k] * qsf_205[k];

        t_306[k] = f_17 * osf_206[k]
                   + f_3 * pc_x[k] * qsf_206[k];

        t_307[k] = f_17 * osf_207[k]
                   + f_3 * pc_x[k] * qsf_207[k];

        t_308[k] = f_3 * pc_y[k] * qsf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_x, pc_y, osf_209, qsd0_123, qsd0_124, \
                         qsd1_123, qsd1_124, qsf_206, qsf_207, \
                         qsf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_17 * osf_209[k]
                   + f_3 * pc_x[k] * qsf_209[k];

        t_310[k] = f_1 * qsd0_123[k]
                   - f_2 * qsd1_123[k]
                   + f_3 * pc_y[k] * qsf_206[k];

        t_311[k] = f_10 * qsd0_124[k]
                   - f_11 * qsd1_124[k]
                   + f_3 * pc_y[k] * qsf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, osf_149, osf_210, \
                         qsd0_125, qsd0_126, qsd1_125, qsd1_126, qsf_208, qsf_209, \
                         qsf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * qsd0_125[k]
                   - f_5 * qsd1_125[k]
                   + f_3 * pc_y[k] * qsf_208[k];

        t_313[k] = f_3 * pc_y[k] * qsf_209[k];

        t_314[k] = f_18 * osf_149[k]
                   + f_1 * qsd0_125[k]
                   - f_2 * qsd1_125[k]
                   + f_3 * pc_z[k] * qsf_209[k];

        t_315[k] = f_19 * osf_210[k]
                   + f_1 * qsd0_126[k]
                   - f_2 * qsd1_126[k]
                   + f_3 * pc_x[k] * qsf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, osf_150, osf_213, \
                         qsd0_129, qsd1_129, qsf_210, qsf_211, \
                         qsf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_19 * osf_150[k]
                   + f_3 * pc_y[k] * qsf_210[k];

        t_317[k] = f_3 * pc_z[k] * qsf_210[k];

        t_318[k] = f_19 * osf_213[k]
                   + f_4 * qsd0_129[k]
                   - f_5 * qsd1_129[k]
                   + f_3 * pc_x[k] * qsf_213[k];

        t_319[k] = f_3 * pc_z[k] * qsf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, pc_z, osf_216, osf_218, qsd0_126, \
                         qsd1_126, qsf_212, qsf_213, qsf_216, qsf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * qsd0_126[k]
                   - f_5 * qsd1_126[k]
                   + f_3 * pc_z[k] * qsf_212[k];

        t_321[k] = f_19 * osf_216[k]
                   + f_3 * pc_x[k] * qsf_216[k];

        t_322[k] = f_3 * pc_z[k] * qsf_213[k];

        t_323[k] = f_19 * osf_218[k]
                   + f_3 * pc_x[k] * qsf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, osf_156, \
                         osf_159, osf_219, qsd0_129, qsd1_129, qsf_216, qsf_217, \
                         qsf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_19 * osf_219[k]
                   + f_3 * pc_x[k] * qsf_219[k];

        t_325[k] = f_19 * osf_156[k]
                   + f_1 * qsd0_129[k]
                   - f_2 * qsd1_129[k]
                   + f_3 * pc_y[k] * qsf_216[k];

        t_326[k] = f_3 * pc_z[k] * qsf_216[k];

        t_327[k] = f_4 * qsd0_129[k]
                   - f_5 * qsd1_129[k]
                   + f_3 * pc_z[k] * qsf_217[k];

        t_328[k] = f_19 * osf_159[k]
                   + f_3 * pc_y[k] * qsf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pa_z, pc_y, pc_z, osg0_225, osf_150, \
                         osf_160, osg1_225, qsd0_131, qsd1_131, qsf_219, \
                         qsf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * qsd0_131[k]
                   - f_2 * qsd1_131[k]
                   + f_3 * pc_z[k] * qsf_219[k];

        t_330[k] = pa_z[k] * osg0_225[k]
                   - f_6 * pc_z[k] * osg1_225[k];

        t_331[k] = f_18 * osf_160[k]
                   + f_3 * pc_y[k] * qsf_220[k];

        t_332[k] = f_7 * osf_150[k]
                   + f_3 * pc_z[k] * qsf_220[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_z, pc_x, pc_y, pc_z, osg0_228, osf_162, \
                         osf_225, osg1_228, qsd0_137, qsd1_137, qsf_222, \
                         qsf_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * osg0_228[k]
                   - f_6 * pc_z[k] * osg1_228[k];

        t_334[k] = f_18 * osf_162[k]
                   + f_3 * pc_y[k] * qsf_222[k];

        t_335[k] = f_19 * osf_225[k]
                   + f_4 * qsd0_137[k]
                   - f_5 * qsd1_137[k]
                   + f_3 * pc_x[k] * qsf_225[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_x, osf_226, osf_227, osf_228, osf_229, \
                         qsf_226, qsf_227, qsf_228, qsf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_19 * osf_226[k]
                   + f_3 * pc_x[k] * qsf_226[k];

        t_337[k] = f_19 * osf_227[k]
                   + f_3 * pc_x[k] * qsf_227[k];

        t_338[k] = f_19 * osf_228[k]
                   + f_3 * pc_x[k] * qsf_228[k];

        t_339[k] = f_19 * osf_229[k]
                   + f_3 * pc_x[k] * qsf_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_z, pc_y, pc_z, osg0_235, osf_156, osf_168, \
                         osg1_235, qsd0_137, qsd1_137, qsf_226, \
                         qsf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_z[k] * osg0_235[k]
                   - f_6 * pc_z[k] * osg1_235[k];

        t_341[k] = f_7 * osf_156[k]
                   + f_3 * pc_z[k] * qsf_226[k];

        t_342[k] = f_18 * osf_168[k]
                   + f_4 * qsd0_137[k]
                   - f_5 * qsd1_137[k]
                   + f_3 * pc_y[k] * qsf_228[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, pc_y, pc_z, osf_159, osf_169, osf_230, \
                         qsd0_137, qsd0_138, qsd1_137, qsd1_138, qsf_229, \
                         qsf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_18 * osf_169[k]
                   + f_3 * pc_y[k] * qsf_229[k];

        t_344[k] = f_7 * osf_159[k]
                   + f_1 * qsd0_137[k]
                   - f_2 * qsd1_137[k]
                   + f_3 * pc_z[k] * qsf_229[k];

        t_345[k] = f_19 * osf_230[k]
                   + f_1 * qsd0_138[k]
                   - f_2 * qsd1_138[k]
                   + f_3 * pc_x[k] * qsf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pc_x, pc_y, pc_z, osf_160, osf_170, \
                         osf_172, osf_233, qsd0_141, qsd1_141, qsf_230, qsf_232, \
                         qsf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_16 * osf_170[k]
                   + f_3 * pc_y[k] * qsf_230[k];

        t_347[k] = f_8 * osf_160[k]
                   + f_3 * pc_z[k] * qsf_230[k];

        t_348[k] = f_19 * osf_233[k]
                   + f_4 * qsd0_141[k]
                   - f_5 * qsd1_141[k]
                   + f_3 * pc_x[k] * qsf_233[k];

        t_349[k] = f_16 * osf_172[k]
                   + f_3 * pc_y[k] * qsf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pc_x, osf_235, osf_236, osf_237, osf_238, \
                         qsd0_143, qsd1_143, qsf_235, qsf_236, qsf_237, \
                         qsf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_19 * osf_235[k]
                   + f_4 * qsd0_143[k]
                   - f_5 * qsd1_143[k]
                   + f_3 * pc_x[k] * qsf_235[k];

        t_351[k] = f_19 * osf_236[k]
                   + f_3 * pc_x[k] * qsf_236[k];

        t_352[k] = f_19 * osf_237[k]
                   + f_3 * pc_x[k] * qsf_237[k];

        t_353[k] = f_19 * osf_238[k]
                   + f_3 * pc_x[k] * qsf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_x, pc_y, pc_z, osf_166, osf_176, osf_239, \
                         qsd0_141, qsd1_141, qsf_236, qsf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_19 * osf_239[k]
                   + f_3 * pc_x[k] * qsf_239[k];

        t_355[k] = f_16 * osf_176[k]
                   + f_1 * qsd0_141[k]
                   - f_2 * qsd1_141[k]
                   + f_3 * pc_y[k] * qsf_236[k];

        t_356[k] = f_8 * osf_166[k]
                   + f_3 * pc_z[k] * qsf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, osf_169, osf_178, osf_179, qsd0_143, \
                         qsd1_143, qsf_238, qsf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_16 * osf_178[k]
                   + f_4 * qsd0_143[k]
                   - f_5 * qsd1_143[k]
                   + f_3 * pc_y[k] * qsf_238[k];

        t_358[k] = f_16 * osf_179[k]
                   + f_3 * pc_y[k] * qsf_239[k];

        t_359[k] = f_8 * osf_169[k]
                   + f_1 * qsd0_143[k]
                   - f_2 * qsd1_143[k]
                   + f_3 * pc_z[k] * qsf_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_x, pc_y, pc_z, osf_170, osf_180, osf_240, \
                         qsd0_144, qsd1_144, qsf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_19 * osf_240[k]
                   + f_1 * qsd0_144[k]
                   - f_2 * qsd1_144[k]
                   + f_3 * pc_x[k] * qsf_240[k];

        t_361[k] = f_14 * osf_180[k]
                   + f_3 * pc_y[k] * qsf_240[k];

        t_362[k] = f_14 * osf_170[k]
                   + f_3 * pc_z[k] * qsf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_x, pc_y, osf_182, osf_243, osf_245, qsd0_147, \
                         qsd0_149, qsd1_147, qsd1_149, qsf_242, qsf_243, \
                         qsf_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_19 * osf_243[k]
                   + f_4 * qsd0_147[k]
                   - f_5 * qsd1_147[k]
                   + f_3 * pc_x[k] * qsf_243[k];

        t_364[k] = f_14 * osf_182[k]
                   + f_3 * pc_y[k] * qsf_242[k];

        t_365[k] = f_19 * osf_245[k]
                   + f_4 * qsd0_149[k]
                   - f_5 * qsd1_149[k]
                   + f_3 * pc_x[k] * qsf_245[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, osf_246, osf_247, osf_248, osf_249, \
                         qsf_246, qsf_247, qsf_248, qsf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_19 * osf_246[k]
                   + f_3 * pc_x[k] * qsf_246[k];

        t_367[k] = f_19 * osf_247[k]
                   + f_3 * pc_x[k] * qsf_247[k];

        t_368[k] = f_19 * osf_248[k]
                   + f_3 * pc_x[k] * qsf_248[k];

        t_369[k] = f_19 * osf_249[k]
                   + f_3 * pc_x[k] * qsf_249[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_300 = buffer.data(osg0 + 300);
    const auto *osg0_303 = buffer.data(osg0 + 303);
    const auto *osg0_305 = buffer.data(osg0 + 305);
    const auto *osg0_314 = buffer.data(osg0 + 314);
    const auto *osg0_315 = buffer.data(osg0 + 315);
    const auto *osg0_318 = buffer.data(osg0 + 318);
    const auto *osg0_325 = buffer.data(osg0 + 325);

    const auto *osf_176 = buffer.data(osf + 176);
    const auto *osf_179 = buffer.data(osf + 179);
    const auto *osf_180 = buffer.data(osf + 180);
    const auto *osf_186 = buffer.data(osf + 186);
    const auto *osf_188 = buffer.data(osf + 188);
    const auto *osf_189 = buffer.data(osf + 189);
    const auto *osf_190 = buffer.data(osf + 190);
    const auto *osf_192 = buffer.data(osf + 192);
    const auto *osf_196 = buffer.data(osf + 196);
    const auto *osf_198 = buffer.data(osf + 198);
    const auto *osf_199 = buffer.data(osf + 199);
    const auto *osf_200 = buffer.data(osf + 200);
    const auto *osf_201 = buffer.data(osf + 201);
    const auto *osf_202 = buffer.data(osf + 202);
    const auto *osf_206 = buffer.data(osf + 206);
    const auto *osf_208 = buffer.data(osf + 208);
    const auto *osf_209 = buffer.data(osf + 209);
    const auto *osf_210 = buffer.data(osf + 210);
    const auto *osf_216 = buffer.data(osf + 216);
    const auto *osf_219 = buffer.data(osf + 219);
    const auto *osf_220 = buffer.data(osf + 220);
    const auto *osf_222 = buffer.data(osf + 222);
    const auto *osf_226 = buffer.data(osf + 226);
    const auto *osf_228 = buffer.data(osf + 228);
    const auto *osf_229 = buffer.data(osf + 229);
    const auto *osf_230 = buffer.data(osf + 230);
    const auto *osf_232 = buffer.data(osf + 232);
    const auto *osf_236 = buffer.data(osf + 236);
    const auto *osf_238 = buffer.data(osf + 238);
    const auto *osf_239 = buffer.data(osf + 239);
    const auto *osf_240 = buffer.data(osf + 240);
    const auto *osf_242 = buffer.data(osf + 242);
    const auto *osf_246 = buffer.data(osf + 246);
    const auto *osf_248 = buffer.data(osf + 248);
    const auto *osf_249 = buffer.data(osf + 249);
    const auto *osf_250 = buffer.data(osf + 250);
    const auto *osf_252 = buffer.data(osf + 252);
    const auto *osf_253 = buffer.data(osf + 253);
    const auto *osf_255 = buffer.data(osf + 255);
    const auto *osf_256 = buffer.data(osf + 256);
    const auto *osf_257 = buffer.data(osf + 257);
    const auto *osf_258 = buffer.data(osf + 258);
    const auto *osf_259 = buffer.data(osf + 259);
    const auto *osf_266 = buffer.data(osf + 266);
    const auto *osf_267 = buffer.data(osf + 267);
    const auto *osf_268 = buffer.data(osf + 268);
    const auto *osf_269 = buffer.data(osf + 269);
    const auto *osf_270 = buffer.data(osf + 270);
    const auto *osf_275 = buffer.data(osf + 275);
    const auto *osf_276 = buffer.data(osf + 276);
    const auto *osf_277 = buffer.data(osf + 277);
    const auto *osf_279 = buffer.data(osf + 279);
    const auto *osf_280 = buffer.data(osf + 280);
    const auto *osf_283 = buffer.data(osf + 283);
    const auto *osf_286 = buffer.data(osf + 286);
    const auto *osf_288 = buffer.data(osf + 288);
    const auto *osf_289 = buffer.data(osf + 289);
    const auto *osf_295 = buffer.data(osf + 295);
    const auto *osf_296 = buffer.data(osf + 296);
    const auto *osf_297 = buffer.data(osf + 297);
    const auto *osf_298 = buffer.data(osf + 298);
    const auto *osf_299 = buffer.data(osf + 299);
    const auto *osf_300 = buffer.data(osf + 300);
    const auto *osf_303 = buffer.data(osf + 303);
    const auto *osf_305 = buffer.data(osf + 305);
    const auto *osf_306 = buffer.data(osf + 306);
    const auto *osf_307 = buffer.data(osf + 307);
    const auto *osf_308 = buffer.data(osf + 308);
    const auto *osf_309 = buffer.data(osf + 309);
    const auto *osf_310 = buffer.data(osf + 310);
    const auto *osf_313 = buffer.data(osf + 313);
    const auto *osf_315 = buffer.data(osf + 315);
    const auto *osf_316 = buffer.data(osf + 316);
    const auto *osf_317 = buffer.data(osf + 317);
    const auto *osf_318 = buffer.data(osf + 318);
    const auto *osf_319 = buffer.data(osf + 319);
    const auto *osf_320 = buffer.data(osf + 320);
    const auto *osf_323 = buffer.data(osf + 323);

    const auto *osg1_300 = buffer.data(osg1 + 300);
    const auto *osg1_303 = buffer.data(osg1 + 303);
    const auto *osg1_305 = buffer.data(osg1 + 305);
    const auto *osg1_314 = buffer.data(osg1 + 314);
    const auto *osg1_315 = buffer.data(osg1 + 315);
    const auto *osg1_318 = buffer.data(osg1 + 318);
    const auto *osg1_325 = buffer.data(osg1 + 325);

    const auto *qsd0_147 = buffer.data(qsd0 + 147);
    const auto *qsd0_149 = buffer.data(qsd0 + 149);
    const auto *qsd0_150 = buffer.data(qsd0 + 150);
    const auto *qsd0_153 = buffer.data(qsd0 + 153);
    const auto *qsd0_155 = buffer.data(qsd0 + 155);
    const auto *qsd0_159 = buffer.data(qsd0 + 159);
    const auto *qsd0_161 = buffer.data(qsd0 + 161);
    const auto *qsd0_162 = buffer.data(qsd0 + 162);
    const auto *qsd0_165 = buffer.data(qsd0 + 165);
    const auto *qsd0_166 = buffer.data(qsd0 + 166);
    const auto *qsd0_167 = buffer.data(qsd0 + 167);
    const auto *qsd0_168 = buffer.data(qsd0 + 168);
    const auto *qsd0_171 = buffer.data(qsd0 + 171);
    const auto *qsd0_173 = buffer.data(qsd0 + 173);
    const auto *qsd0_179 = buffer.data(qsd0 + 179);
    const auto *qsd0_180 = buffer.data(qsd0 + 180);
    const auto *qsd0_183 = buffer.data(qsd0 + 183);
    const auto *qsd0_185 = buffer.data(qsd0 + 185);
    const auto *qsd0_186 = buffer.data(qsd0 + 186);
    const auto *qsd0_189 = buffer.data(qsd0 + 189);
    const auto *qsd0_191 = buffer.data(qsd0 + 191);
    const auto *qsd0_192 = buffer.data(qsd0 + 192);
    const auto *qsd0_195 = buffer.data(qsd0 + 195);

    const auto *qsd1_147 = buffer.data(qsd1 + 147);
    const auto *qsd1_149 = buffer.data(qsd1 + 149);
    const auto *qsd1_150 = buffer.data(qsd1 + 150);
    const auto *qsd1_153 = buffer.data(qsd1 + 153);
    const auto *qsd1_155 = buffer.data(qsd1 + 155);
    const auto *qsd1_159 = buffer.data(qsd1 + 159);
    const auto *qsd1_161 = buffer.data(qsd1 + 161);
    const auto *qsd1_162 = buffer.data(qsd1 + 162);
    const auto *qsd1_165 = buffer.data(qsd1 + 165);
    const auto *qsd1_166 = buffer.data(qsd1 + 166);
    const auto *qsd1_167 = buffer.data(qsd1 + 167);
    const auto *qsd1_168 = buffer.data(qsd1 + 168);
    const auto *qsd1_171 = buffer.data(qsd1 + 171);
    const auto *qsd1_173 = buffer.data(qsd1 + 173);
    const auto *qsd1_179 = buffer.data(qsd1 + 179);
    const auto *qsd1_180 = buffer.data(qsd1 + 180);
    const auto *qsd1_183 = buffer.data(qsd1 + 183);
    const auto *qsd1_185 = buffer.data(qsd1 + 185);
    const auto *qsd1_186 = buffer.data(qsd1 + 186);
    const auto *qsd1_189 = buffer.data(qsd1 + 189);
    const auto *qsd1_191 = buffer.data(qsd1 + 191);
    const auto *qsd1_192 = buffer.data(qsd1 + 192);
    const auto *qsd1_195 = buffer.data(qsd1 + 195);

    const auto *qsf_246 = buffer.data(qsf + 246);
    const auto *qsf_248 = buffer.data(qsf + 248);
    const auto *qsf_249 = buffer.data(qsf + 249);
    const auto *qsf_250 = buffer.data(qsf + 250);
    const auto *qsf_252 = buffer.data(qsf + 252);
    const auto *qsf_253 = buffer.data(qsf + 253);
    const auto *qsf_255 = buffer.data(qsf + 255);
    const auto *qsf_256 = buffer.data(qsf + 256);
    const auto *qsf_257 = buffer.data(qsf + 257);
    const auto *qsf_258 = buffer.data(qsf + 258);
    const auto *qsf_259 = buffer.data(qsf + 259);
    const auto *qsf_260 = buffer.data(qsf + 260);
    const auto *qsf_262 = buffer.data(qsf + 262);
    const auto *qsf_266 = buffer.data(qsf + 266);
    const auto *qsf_267 = buffer.data(qsf + 267);
    const auto *qsf_268 = buffer.data(qsf + 268);
    const auto *qsf_269 = buffer.data(qsf + 269);
    const auto *qsf_270 = buffer.data(qsf + 270);
    const auto *qsf_271 = buffer.data(qsf + 271);
    const auto *qsf_272 = buffer.data(qsf + 272);
    const auto *qsf_275 = buffer.data(qsf + 275);
    const auto *qsf_276 = buffer.data(qsf + 276);
    const auto *qsf_277 = buffer.data(qsf + 277);
    const auto *qsf_278 = buffer.data(qsf + 278);
    const auto *qsf_279 = buffer.data(qsf + 279);
    const auto *qsf_280 = buffer.data(qsf + 280);
    const auto *qsf_281 = buffer.data(qsf + 281);
    const auto *qsf_282 = buffer.data(qsf + 282);
    const auto *qsf_283 = buffer.data(qsf + 283);
    const auto *qsf_286 = buffer.data(qsf + 286);
    const auto *qsf_287 = buffer.data(qsf + 287);
    const auto *qsf_288 = buffer.data(qsf + 288);
    const auto *qsf_289 = buffer.data(qsf + 289);
    const auto *qsf_290 = buffer.data(qsf + 290);
    const auto *qsf_292 = buffer.data(qsf + 292);
    const auto *qsf_295 = buffer.data(qsf + 295);
    const auto *qsf_296 = buffer.data(qsf + 296);
    const auto *qsf_297 = buffer.data(qsf + 297);
    const auto *qsf_298 = buffer.data(qsf + 298);
    const auto *qsf_299 = buffer.data(qsf + 299);
    const auto *qsf_300 = buffer.data(qsf + 300);
    const auto *qsf_302 = buffer.data(qsf + 302);
    const auto *qsf_303 = buffer.data(qsf + 303);
    const auto *qsf_305 = buffer.data(qsf + 305);
    const auto *qsf_306 = buffer.data(qsf + 306);
    const auto *qsf_307 = buffer.data(qsf + 307);
    const auto *qsf_308 = buffer.data(qsf + 308);
    const auto *qsf_309 = buffer.data(qsf + 309);
    const auto *qsf_310 = buffer.data(qsf + 310);
    const auto *qsf_312 = buffer.data(qsf + 312);
    const auto *qsf_313 = buffer.data(qsf + 313);
    const auto *qsf_315 = buffer.data(qsf + 315);
    const auto *qsf_316 = buffer.data(qsf + 316);
    const auto *qsf_317 = buffer.data(qsf + 317);
    const auto *qsf_318 = buffer.data(qsf + 318);
    const auto *qsf_319 = buffer.data(qsf + 319);
    const auto *qsf_320 = buffer.data(qsf + 320);
    const auto *qsf_322 = buffer.data(qsf + 322);
    const auto *qsf_323 = buffer.data(qsf + 323);

#pragma omp simd aligned(t_370, t_371, t_372, pc_y, pc_z, osf_176, osf_186, osf_188, qsd0_147, \
                         qsd0_149, qsd1_147, qsd1_149, qsf_246, \
                         qsf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * osf_186[k]
                   + f_1 * qsd0_147[k]
                   - f_2 * qsd1_147[k]
                   + f_3 * pc_y[k] * qsf_246[k];

        t_371[k] = f_14 * osf_176[k]
                   + f_3 * pc_z[k] * qsf_246[k];

        t_372[k] = f_14 * osf_188[k]
                   + f_4 * qsd0_149[k]
                   - f_5 * qsd1_149[k]
                   + f_3 * pc_y[k] * qsf_248[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, pc_z, osf_179, osf_189, osf_250, \
                         qsd0_149, qsd0_150, qsd1_149, qsd1_150, qsf_249, \
                         qsf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * osf_189[k]
                   + f_3 * pc_y[k] * qsf_249[k];

        t_374[k] = f_14 * osf_179[k]
                   + f_1 * qsd0_149[k]
                   - f_2 * qsd1_149[k]
                   + f_3 * pc_z[k] * qsf_249[k];

        t_375[k] = f_19 * osf_250[k]
                   + f_1 * qsd0_150[k]
                   - f_2 * qsd1_150[k]
                   + f_3 * pc_x[k] * qsf_250[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pc_x, pc_y, pc_z, osf_180, osf_190, \
                         osf_192, osf_253, qsd0_153, qsd1_153, qsf_250, qsf_252, \
                         qsf_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_8 * osf_190[k]
                   + f_3 * pc_y[k] * qsf_250[k];

        t_377[k] = f_16 * osf_180[k]
                   + f_3 * pc_z[k] * qsf_250[k];

        t_378[k] = f_19 * osf_253[k]
                   + f_4 * qsd0_153[k]
                   - f_5 * qsd1_153[k]
                   + f_3 * pc_x[k] * qsf_253[k];

        t_379[k] = f_8 * osf_192[k]
                   + f_3 * pc_y[k] * qsf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, osf_255, osf_256, osf_257, osf_258, \
                         qsd0_155, qsd1_155, qsf_255, qsf_256, qsf_257, \
                         qsf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_19 * osf_255[k]
                   + f_4 * qsd0_155[k]
                   - f_5 * qsd1_155[k]
                   + f_3 * pc_x[k] * qsf_255[k];

        t_381[k] = f_19 * osf_256[k]
                   + f_3 * pc_x[k] * qsf_256[k];

        t_382[k] = f_19 * osf_257[k]
                   + f_3 * pc_x[k] * qsf_257[k];

        t_383[k] = f_19 * osf_258[k]
                   + f_3 * pc_x[k] * qsf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, osf_186, osf_196, osf_259, \
                         qsd0_153, qsd1_153, qsf_256, qsf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_19 * osf_259[k]
                   + f_3 * pc_x[k] * qsf_259[k];

        t_385[k] = f_8 * osf_196[k]
                   + f_1 * qsd0_153[k]
                   - f_2 * qsd1_153[k]
                   + f_3 * pc_y[k] * qsf_256[k];

        t_386[k] = f_16 * osf_186[k]
                   + f_3 * pc_z[k] * qsf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_y, pc_z, osg0_300, osf_189, \
                         osf_198, osf_199, osg1_300, qsd0_155, qsd1_155, qsf_258, \
                         qsf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * osf_198[k]
                   + f_4 * qsd0_155[k]
                   - f_5 * qsd1_155[k]
                   + f_3 * pc_y[k] * qsf_258[k];

        t_388[k] = f_8 * osf_199[k]
                   + f_3 * pc_y[k] * qsf_259[k];

        t_389[k] = f_16 * osf_189[k]
                   + f_1 * qsd0_155[k]
                   - f_2 * qsd1_155[k]
                   + f_3 * pc_z[k] * qsf_259[k];

        t_390[k] = pa_y[k] * osg0_300[k]
                   - f_6 * pc_y[k] * osg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_y, pc_z, osg0_303, osf_190, \
                         osf_200, osf_201, osf_202, osg1_303, qsf_260, \
                         qsf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_7 * osf_200[k]
                   + f_3 * pc_y[k] * qsf_260[k];

        t_392[k] = f_18 * osf_190[k]
                   + f_3 * pc_z[k] * qsf_260[k];

        t_393[k] = pa_y[k] * osg0_303[k]
                   + f_8 * osf_201[k]
                   - f_6 * pc_y[k] * osg1_303[k];

        t_394[k] = f_7 * osf_202[k]
                   + f_3 * pc_y[k] * qsf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, osg0_305, osf_266, \
                         osf_267, osf_268, osg1_305, qsf_266, qsf_267, \
                         qsf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * osg0_305[k]
                   - f_6 * pc_y[k] * osg1_305[k];

        t_396[k] = f_19 * osf_266[k]
                   + f_3 * pc_x[k] * qsf_266[k];

        t_397[k] = f_19 * osf_267[k]
                   + f_3 * pc_x[k] * qsf_267[k];

        t_398[k] = f_19 * osf_268[k]
                   + f_3 * pc_x[k] * qsf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, osf_196, osf_206, osf_269, \
                         qsd0_159, qsd1_159, qsf_266, qsf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_19 * osf_269[k]
                   + f_3 * pc_x[k] * qsf_269[k];

        t_400[k] = f_7 * osf_206[k]
                   + f_1 * qsd0_159[k]
                   - f_2 * qsd1_159[k]
                   + f_3 * pc_y[k] * qsf_266[k];

        t_401[k] = f_18 * osf_196[k]
                   + f_3 * pc_z[k] * qsf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, osg0_314, osf_208, osf_209, \
                         osg1_314, qsd0_161, qsd1_161, qsf_268, \
                         qsf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_7 * osf_208[k]
                   + f_4 * qsd0_161[k]
                   - f_5 * qsd1_161[k]
                   + f_3 * pc_y[k] * qsf_268[k];

        t_403[k] = f_7 * osf_209[k]
                   + f_3 * pc_y[k] * qsf_269[k];

        t_404[k] = pa_y[k] * osg0_314[k]
                   - f_6 * pc_y[k] * osg1_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, pc_z, osf_200, \
                         osf_270, qsd0_162, qsd1_162, qsf_270, qsf_271, \
                         qsf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_19 * osf_270[k]
                   + f_1 * qsd0_162[k]
                   - f_2 * qsd1_162[k]
                   + f_3 * pc_x[k] * qsf_270[k];

        t_406[k] = f_3 * pc_y[k] * qsf_270[k];

        t_407[k] = f_19 * osf_200[k]
                   + f_3 * pc_z[k] * qsf_270[k];

        t_408[k] = f_4 * qsd0_162[k]
                   - f_5 * qsd1_162[k]
                   + f_3 * pc_y[k] * qsf_271[k];

        t_409[k] = f_3 * pc_y[k] * qsf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, osf_275, osf_276, osf_277, \
                         qsd0_167, qsd1_167, qsf_275, qsf_276, \
                         qsf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_19 * osf_275[k]
                   + f_4 * qsd0_167[k]
                   - f_5 * qsd1_167[k]
                   + f_3 * pc_x[k] * qsf_275[k];

        t_411[k] = f_19 * osf_276[k]
                   + f_3 * pc_x[k] * qsf_276[k];

        t_412[k] = f_19 * osf_277[k]
                   + f_3 * pc_x[k] * qsf_277[k];

        t_413[k] = f_3 * pc_y[k] * qsf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, osf_279, qsd0_165, qsd0_166, \
                         qsd1_165, qsd1_166, qsf_276, qsf_277, \
                         qsf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_19 * osf_279[k]
                   + f_3 * pc_x[k] * qsf_279[k];

        t_415[k] = f_1 * qsd0_165[k]
                   - f_2 * qsd1_165[k]
                   + f_3 * pc_y[k] * qsf_276[k];

        t_416[k] = f_10 * qsd0_166[k]
                   - f_11 * qsd1_166[k]
                   + f_3 * pc_y[k] * qsf_277[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pc_x, pc_y, pc_z, osf_209, osf_280, \
                         qsd0_167, qsd0_168, qsd1_167, qsd1_168, qsf_278, qsf_279, \
                         qsf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * qsd0_167[k]
                   - f_5 * qsd1_167[k]
                   + f_3 * pc_y[k] * qsf_278[k];

        t_418[k] = f_3 * pc_y[k] * qsf_279[k];

        t_419[k] = f_19 * osf_209[k]
                   + f_1 * qsd0_167[k]
                   - f_2 * qsd1_167[k]
                   + f_3 * pc_z[k] * qsf_279[k];

        t_420[k] = f_18 * osf_280[k]
                   + f_1 * qsd0_168[k]
                   - f_2 * qsd1_168[k]
                   + f_3 * pc_x[k] * qsf_280[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, osf_210, osf_283, \
                         qsd0_171, qsd1_171, qsf_280, qsf_281, \
                         qsf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_17 * osf_210[k]
                   + f_3 * pc_y[k] * qsf_280[k];

        t_422[k] = f_3 * pc_z[k] * qsf_280[k];

        t_423[k] = f_18 * osf_283[k]
                   + f_4 * qsd0_171[k]
                   - f_5 * qsd1_171[k]
                   + f_3 * pc_x[k] * qsf_283[k];

        t_424[k] = f_3 * pc_z[k] * qsf_281[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_z, osf_286, osf_288, qsd0_168, \
                         qsd1_168, qsf_282, qsf_283, qsf_286, qsf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_4 * qsd0_168[k]
                   - f_5 * qsd1_168[k]
                   + f_3 * pc_z[k] * qsf_282[k];

        t_426[k] = f_18 * osf_286[k]
                   + f_3 * pc_x[k] * qsf_286[k];

        t_427[k] = f_3 * pc_z[k] * qsf_283[k];

        t_428[k] = f_18 * osf_288[k]
                   + f_3 * pc_x[k] * qsf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, pc_x, pc_y, pc_z, osf_216, \
                         osf_219, osf_289, qsd0_171, qsd1_171, qsf_286, qsf_287, \
                         qsf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_18 * osf_289[k]
                   + f_3 * pc_x[k] * qsf_289[k];

        t_430[k] = f_17 * osf_216[k]
                   + f_1 * qsd0_171[k]
                   - f_2 * qsd1_171[k]
                   + f_3 * pc_y[k] * qsf_286[k];

        t_431[k] = f_3 * pc_z[k] * qsf_286[k];

        t_432[k] = f_4 * qsd0_171[k]
                   - f_5 * qsd1_171[k]
                   + f_3 * pc_z[k] * qsf_287[k];

        t_433[k] = f_17 * osf_219[k]
                   + f_3 * pc_y[k] * qsf_289[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_z, pc_y, pc_z, osg0_315, osf_210, \
                         osf_220, osg1_315, qsd0_173, qsd1_173, qsf_289, \
                         qsf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * qsd0_173[k]
                   - f_2 * qsd1_173[k]
                   + f_3 * pc_z[k] * qsf_289[k];

        t_435[k] = pa_z[k] * osg0_315[k]
                   - f_6 * pc_z[k] * osg1_315[k];

        t_436[k] = f_19 * osf_220[k]
                   + f_3 * pc_y[k] * qsf_290[k];

        t_437[k] = f_7 * osf_210[k]
                   + f_3 * pc_z[k] * qsf_290[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pa_z, pc_x, pc_y, pc_z, osg0_318, osf_222, \
                         osf_295, osg1_318, qsd0_179, qsd1_179, qsf_292, \
                         qsf_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = pa_z[k] * osg0_318[k]
                   - f_6 * pc_z[k] * osg1_318[k];

        t_439[k] = f_19 * osf_222[k]
                   + f_3 * pc_y[k] * qsf_292[k];

        t_440[k] = f_18 * osf_295[k]
                   + f_4 * qsd0_179[k]
                   - f_5 * qsd1_179[k]
                   + f_3 * pc_x[k] * qsf_295[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, osf_296, osf_297, osf_298, osf_299, \
                         qsf_296, qsf_297, qsf_298, qsf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_18 * osf_296[k]
                   + f_3 * pc_x[k] * qsf_296[k];

        t_442[k] = f_18 * osf_297[k]
                   + f_3 * pc_x[k] * qsf_297[k];

        t_443[k] = f_18 * osf_298[k]
                   + f_3 * pc_x[k] * qsf_298[k];

        t_444[k] = f_18 * osf_299[k]
                   + f_3 * pc_x[k] * qsf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_z, pc_y, pc_z, osg0_325, osf_216, osf_228, \
                         osg1_325, qsd0_179, qsd1_179, qsf_296, \
                         qsf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_z[k] * osg0_325[k]
                   - f_6 * pc_z[k] * osg1_325[k];

        t_446[k] = f_7 * osf_216[k]
                   + f_3 * pc_z[k] * qsf_296[k];

        t_447[k] = f_19 * osf_228[k]
                   + f_4 * qsd0_179[k]
                   - f_5 * qsd1_179[k]
                   + f_3 * pc_y[k] * qsf_298[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pc_x, pc_y, pc_z, osf_219, osf_229, osf_300, \
                         qsd0_179, qsd0_180, qsd1_179, qsd1_180, qsf_299, \
                         qsf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_19 * osf_229[k]
                   + f_3 * pc_y[k] * qsf_299[k];

        t_449[k] = f_7 * osf_219[k]
                   + f_1 * qsd0_179[k]
                   - f_2 * qsd1_179[k]
                   + f_3 * pc_z[k] * qsf_299[k];

        t_450[k] = f_18 * osf_300[k]
                   + f_1 * qsd0_180[k]
                   - f_2 * qsd1_180[k]
                   + f_3 * pc_x[k] * qsf_300[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, osf_220, osf_230, \
                         osf_232, osf_303, qsd0_183, qsd1_183, qsf_300, qsf_302, \
                         qsf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * osf_230[k]
                   + f_3 * pc_y[k] * qsf_300[k];

        t_452[k] = f_8 * osf_220[k]
                   + f_3 * pc_z[k] * qsf_300[k];

        t_453[k] = f_18 * osf_303[k]
                   + f_4 * qsd0_183[k]
                   - f_5 * qsd1_183[k]
                   + f_3 * pc_x[k] * qsf_303[k];

        t_454[k] = f_18 * osf_232[k]
                   + f_3 * pc_y[k] * qsf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pc_x, osf_305, osf_306, osf_307, osf_308, \
                         qsd0_185, qsd1_185, qsf_305, qsf_306, qsf_307, \
                         qsf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_18 * osf_305[k]
                   + f_4 * qsd0_185[k]
                   - f_5 * qsd1_185[k]
                   + f_3 * pc_x[k] * qsf_305[k];

        t_456[k] = f_18 * osf_306[k]
                   + f_3 * pc_x[k] * qsf_306[k];

        t_457[k] = f_18 * osf_307[k]
                   + f_3 * pc_x[k] * qsf_307[k];

        t_458[k] = f_18 * osf_308[k]
                   + f_3 * pc_x[k] * qsf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_y, pc_z, osf_226, osf_236, osf_309, \
                         qsd0_183, qsd1_183, qsf_306, qsf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_18 * osf_309[k]
                   + f_3 * pc_x[k] * qsf_309[k];

        t_460[k] = f_18 * osf_236[k]
                   + f_1 * qsd0_183[k]
                   - f_2 * qsd1_183[k]
                   + f_3 * pc_y[k] * qsf_306[k];

        t_461[k] = f_8 * osf_226[k]
                   + f_3 * pc_z[k] * qsf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, osf_229, osf_238, osf_239, qsd0_185, \
                         qsd1_185, qsf_308, qsf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_18 * osf_238[k]
                   + f_4 * qsd0_185[k]
                   - f_5 * qsd1_185[k]
                   + f_3 * pc_y[k] * qsf_308[k];

        t_463[k] = f_18 * osf_239[k]
                   + f_3 * pc_y[k] * qsf_309[k];

        t_464[k] = f_8 * osf_229[k]
                   + f_1 * qsd0_185[k]
                   - f_2 * qsd1_185[k]
                   + f_3 * pc_z[k] * qsf_309[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_y, pc_z, osf_230, osf_240, osf_310, \
                         qsd0_186, qsd1_186, qsf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_18 * osf_310[k]
                   + f_1 * qsd0_186[k]
                   - f_2 * qsd1_186[k]
                   + f_3 * pc_x[k] * qsf_310[k];

        t_466[k] = f_16 * osf_240[k]
                   + f_3 * pc_y[k] * qsf_310[k];

        t_467[k] = f_14 * osf_230[k]
                   + f_3 * pc_z[k] * qsf_310[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, osf_242, osf_313, osf_315, qsd0_189, \
                         qsd0_191, qsd1_189, qsd1_191, qsf_312, qsf_313, \
                         qsf_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_18 * osf_313[k]
                   + f_4 * qsd0_189[k]
                   - f_5 * qsd1_189[k]
                   + f_3 * pc_x[k] * qsf_313[k];

        t_469[k] = f_16 * osf_242[k]
                   + f_3 * pc_y[k] * qsf_312[k];

        t_470[k] = f_18 * osf_315[k]
                   + f_4 * qsd0_191[k]
                   - f_5 * qsd1_191[k]
                   + f_3 * pc_x[k] * qsf_315[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, osf_316, osf_317, osf_318, osf_319, \
                         qsf_316, qsf_317, qsf_318, qsf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_18 * osf_316[k]
                   + f_3 * pc_x[k] * qsf_316[k];

        t_472[k] = f_18 * osf_317[k]
                   + f_3 * pc_x[k] * qsf_317[k];

        t_473[k] = f_18 * osf_318[k]
                   + f_3 * pc_x[k] * qsf_318[k];

        t_474[k] = f_18 * osf_319[k]
                   + f_3 * pc_x[k] * qsf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pc_y, pc_z, osf_236, osf_246, osf_248, qsd0_189, \
                         qsd0_191, qsd1_189, qsd1_191, qsf_316, \
                         qsf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_16 * osf_246[k]
                   + f_1 * qsd0_189[k]
                   - f_2 * qsd1_189[k]
                   + f_3 * pc_y[k] * qsf_316[k];

        t_476[k] = f_14 * osf_236[k]
                   + f_3 * pc_z[k] * qsf_316[k];

        t_477[k] = f_16 * osf_248[k]
                   + f_4 * qsd0_191[k]
                   - f_5 * qsd1_191[k]
                   + f_3 * pc_y[k] * qsf_318[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, osf_239, osf_249, osf_320, \
                         qsd0_191, qsd0_192, qsd1_191, qsd1_192, qsf_319, \
                         qsf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_16 * osf_249[k]
                   + f_3 * pc_y[k] * qsf_319[k];

        t_479[k] = f_14 * osf_239[k]
                   + f_1 * qsd0_191[k]
                   - f_2 * qsd1_191[k]
                   + f_3 * pc_z[k] * qsf_319[k];

        t_480[k] = f_18 * osf_320[k]
                   + f_1 * qsd0_192[k]
                   - f_2 * qsd1_192[k]
                   + f_3 * pc_x[k] * qsf_320[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, osf_240, osf_250, \
                         osf_252, osf_323, qsd0_195, qsd1_195, qsf_320, qsf_322, \
                         qsf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_14 * osf_250[k]
                   + f_3 * pc_y[k] * qsf_320[k];

        t_482[k] = f_16 * osf_240[k]
                   + f_3 * pc_z[k] * qsf_320[k];

        t_483[k] = f_18 * osf_323[k]
                   + f_4 * qsd0_195[k]
                   - f_5 * qsd1_195[k]
                   + f_3 * pc_x[k] * qsf_323[k];

        t_484[k] = f_14 * osf_252[k]
                   + f_3 * pc_y[k] * qsf_322[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_405 = buffer.data(osg0 + 405);
    const auto *osg0_408 = buffer.data(osg0 + 408);
    const auto *osg0_410 = buffer.data(osg0 + 410);
    const auto *osg0_419 = buffer.data(osg0 + 419);
    const auto *osg0_420 = buffer.data(osg0 + 420);
    const auto *osg0_423 = buffer.data(osg0 + 423);
    const auto *osg0_430 = buffer.data(osg0 + 430);

    const auto *osf_246 = buffer.data(osf + 246);
    const auto *osf_249 = buffer.data(osf + 249);
    const auto *osf_250 = buffer.data(osf + 250);
    const auto *osf_256 = buffer.data(osf + 256);
    const auto *osf_258 = buffer.data(osf + 258);
    const auto *osf_259 = buffer.data(osf + 259);
    const auto *osf_260 = buffer.data(osf + 260);
    const auto *osf_262 = buffer.data(osf + 262);
    const auto *osf_266 = buffer.data(osf + 266);
    const auto *osf_268 = buffer.data(osf + 268);
    const auto *osf_269 = buffer.data(osf + 269);
    const auto *osf_270 = buffer.data(osf + 270);
    const auto *osf_271 = buffer.data(osf + 271);
    const auto *osf_272 = buffer.data(osf + 272);
    const auto *osf_276 = buffer.data(osf + 276);
    const auto *osf_278 = buffer.data(osf + 278);
    const auto *osf_279 = buffer.data(osf + 279);
    const auto *osf_280 = buffer.data(osf + 280);
    const auto *osf_286 = buffer.data(osf + 286);
    const auto *osf_289 = buffer.data(osf + 289);
    const auto *osf_290 = buffer.data(osf + 290);
    const auto *osf_292 = buffer.data(osf + 292);
    const auto *osf_296 = buffer.data(osf + 296);
    const auto *osf_298 = buffer.data(osf + 298);
    const auto *osf_299 = buffer.data(osf + 299);
    const auto *osf_300 = buffer.data(osf + 300);
    const auto *osf_302 = buffer.data(osf + 302);
    const auto *osf_306 = buffer.data(osf + 306);
    const auto *osf_308 = buffer.data(osf + 308);
    const auto *osf_309 = buffer.data(osf + 309);
    const auto *osf_310 = buffer.data(osf + 310);
    const auto *osf_312 = buffer.data(osf + 312);
    const auto *osf_316 = buffer.data(osf + 316);
    const auto *osf_318 = buffer.data(osf + 318);
    const auto *osf_319 = buffer.data(osf + 319);
    const auto *osf_325 = buffer.data(osf + 325);
    const auto *osf_326 = buffer.data(osf + 326);
    const auto *osf_327 = buffer.data(osf + 327);
    const auto *osf_328 = buffer.data(osf + 328);
    const auto *osf_329 = buffer.data(osf + 329);
    const auto *osf_330 = buffer.data(osf + 330);
    const auto *osf_333 = buffer.data(osf + 333);
    const auto *osf_335 = buffer.data(osf + 335);
    const auto *osf_336 = buffer.data(osf + 336);
    const auto *osf_337 = buffer.data(osf + 337);
    const auto *osf_338 = buffer.data(osf + 338);
    const auto *osf_339 = buffer.data(osf + 339);
    const auto *osf_346 = buffer.data(osf + 346);
    const auto *osf_347 = buffer.data(osf + 347);
    const auto *osf_348 = buffer.data(osf + 348);
    const auto *osf_349 = buffer.data(osf + 349);
    const auto *osf_350 = buffer.data(osf + 350);
    const auto *osf_355 = buffer.data(osf + 355);
    const auto *osf_356 = buffer.data(osf + 356);
    const auto *osf_357 = buffer.data(osf + 357);
    const auto *osf_359 = buffer.data(osf + 359);
    const auto *osf_360 = buffer.data(osf + 360);
    const auto *osf_363 = buffer.data(osf + 363);
    const auto *osf_366 = buffer.data(osf + 366);
    const auto *osf_368 = buffer.data(osf + 368);
    const auto *osf_369 = buffer.data(osf + 369);
    const auto *osf_375 = buffer.data(osf + 375);
    const auto *osf_376 = buffer.data(osf + 376);
    const auto *osf_377 = buffer.data(osf + 377);
    const auto *osf_378 = buffer.data(osf + 378);
    const auto *osf_379 = buffer.data(osf + 379);
    const auto *osf_380 = buffer.data(osf + 380);
    const auto *osf_383 = buffer.data(osf + 383);
    const auto *osf_385 = buffer.data(osf + 385);
    const auto *osf_386 = buffer.data(osf + 386);
    const auto *osf_387 = buffer.data(osf + 387);
    const auto *osf_388 = buffer.data(osf + 388);
    const auto *osf_389 = buffer.data(osf + 389);
    const auto *osf_390 = buffer.data(osf + 390);
    const auto *osf_393 = buffer.data(osf + 393);
    const auto *osf_395 = buffer.data(osf + 395);
    const auto *osf_396 = buffer.data(osf + 396);
    const auto *osf_397 = buffer.data(osf + 397);
    const auto *osf_398 = buffer.data(osf + 398);
    const auto *osf_399 = buffer.data(osf + 399);
    const auto *osf_400 = buffer.data(osf + 400);

    const auto *osg1_405 = buffer.data(osg1 + 405);
    const auto *osg1_408 = buffer.data(osg1 + 408);
    const auto *osg1_410 = buffer.data(osg1 + 410);
    const auto *osg1_419 = buffer.data(osg1 + 419);
    const auto *osg1_420 = buffer.data(osg1 + 420);
    const auto *osg1_423 = buffer.data(osg1 + 423);
    const auto *osg1_430 = buffer.data(osg1 + 430);

    const auto *qsd0_195 = buffer.data(qsd0 + 195);
    const auto *qsd0_197 = buffer.data(qsd0 + 197);
    const auto *qsd0_198 = buffer.data(qsd0 + 198);
    const auto *qsd0_201 = buffer.data(qsd0 + 201);
    const auto *qsd0_203 = buffer.data(qsd0 + 203);
    const auto *qsd0_207 = buffer.data(qsd0 + 207);
    const auto *qsd0_209 = buffer.data(qsd0 + 209);
    const auto *qsd0_210 = buffer.data(qsd0 + 210);
    const auto *qsd0_213 = buffer.data(qsd0 + 213);
    const auto *qsd0_214 = buffer.data(qsd0 + 214);
    const auto *qsd0_215 = buffer.data(qsd0 + 215);
    const auto *qsd0_216 = buffer.data(qsd0 + 216);
    const auto *qsd0_219 = buffer.data(qsd0 + 219);
    const auto *qsd0_221 = buffer.data(qsd0 + 221);
    const auto *qsd0_227 = buffer.data(qsd0 + 227);
    const auto *qsd0_228 = buffer.data(qsd0 + 228);
    const auto *qsd0_231 = buffer.data(qsd0 + 231);
    const auto *qsd0_233 = buffer.data(qsd0 + 233);
    const auto *qsd0_234 = buffer.data(qsd0 + 234);
    const auto *qsd0_237 = buffer.data(qsd0 + 237);
    const auto *qsd0_239 = buffer.data(qsd0 + 239);
    const auto *qsd0_240 = buffer.data(qsd0 + 240);

    const auto *qsd1_195 = buffer.data(qsd1 + 195);
    const auto *qsd1_197 = buffer.data(qsd1 + 197);
    const auto *qsd1_198 = buffer.data(qsd1 + 198);
    const auto *qsd1_201 = buffer.data(qsd1 + 201);
    const auto *qsd1_203 = buffer.data(qsd1 + 203);
    const auto *qsd1_207 = buffer.data(qsd1 + 207);
    const auto *qsd1_209 = buffer.data(qsd1 + 209);
    const auto *qsd1_210 = buffer.data(qsd1 + 210);
    const auto *qsd1_213 = buffer.data(qsd1 + 213);
    const auto *qsd1_214 = buffer.data(qsd1 + 214);
    const auto *qsd1_215 = buffer.data(qsd1 + 215);
    const auto *qsd1_216 = buffer.data(qsd1 + 216);
    const auto *qsd1_219 = buffer.data(qsd1 + 219);
    const auto *qsd1_221 = buffer.data(qsd1 + 221);
    const auto *qsd1_227 = buffer.data(qsd1 + 227);
    const auto *qsd1_228 = buffer.data(qsd1 + 228);
    const auto *qsd1_231 = buffer.data(qsd1 + 231);
    const auto *qsd1_233 = buffer.data(qsd1 + 233);
    const auto *qsd1_234 = buffer.data(qsd1 + 234);
    const auto *qsd1_237 = buffer.data(qsd1 + 237);
    const auto *qsd1_239 = buffer.data(qsd1 + 239);
    const auto *qsd1_240 = buffer.data(qsd1 + 240);

    const auto *qsf_325 = buffer.data(qsf + 325);
    const auto *qsf_326 = buffer.data(qsf + 326);
    const auto *qsf_327 = buffer.data(qsf + 327);
    const auto *qsf_328 = buffer.data(qsf + 328);
    const auto *qsf_329 = buffer.data(qsf + 329);
    const auto *qsf_330 = buffer.data(qsf + 330);
    const auto *qsf_332 = buffer.data(qsf + 332);
    const auto *qsf_333 = buffer.data(qsf + 333);
    const auto *qsf_335 = buffer.data(qsf + 335);
    const auto *qsf_336 = buffer.data(qsf + 336);
    const auto *qsf_337 = buffer.data(qsf + 337);
    const auto *qsf_338 = buffer.data(qsf + 338);
    const auto *qsf_339 = buffer.data(qsf + 339);
    const auto *qsf_340 = buffer.data(qsf + 340);
    const auto *qsf_342 = buffer.data(qsf + 342);
    const auto *qsf_346 = buffer.data(qsf + 346);
    const auto *qsf_347 = buffer.data(qsf + 347);
    const auto *qsf_348 = buffer.data(qsf + 348);
    const auto *qsf_349 = buffer.data(qsf + 349);
    const auto *qsf_350 = buffer.data(qsf + 350);
    const auto *qsf_351 = buffer.data(qsf + 351);
    const auto *qsf_352 = buffer.data(qsf + 352);
    const auto *qsf_355 = buffer.data(qsf + 355);
    const auto *qsf_356 = buffer.data(qsf + 356);
    const auto *qsf_357 = buffer.data(qsf + 357);
    const auto *qsf_358 = buffer.data(qsf + 358);
    const auto *qsf_359 = buffer.data(qsf + 359);
    const auto *qsf_360 = buffer.data(qsf + 360);
    const auto *qsf_361 = buffer.data(qsf + 361);
    const auto *qsf_362 = buffer.data(qsf + 362);
    const auto *qsf_363 = buffer.data(qsf + 363);
    const auto *qsf_366 = buffer.data(qsf + 366);
    const auto *qsf_367 = buffer.data(qsf + 367);
    const auto *qsf_368 = buffer.data(qsf + 368);
    const auto *qsf_369 = buffer.data(qsf + 369);
    const auto *qsf_370 = buffer.data(qsf + 370);
    const auto *qsf_372 = buffer.data(qsf + 372);
    const auto *qsf_375 = buffer.data(qsf + 375);
    const auto *qsf_376 = buffer.data(qsf + 376);
    const auto *qsf_377 = buffer.data(qsf + 377);
    const auto *qsf_378 = buffer.data(qsf + 378);
    const auto *qsf_379 = buffer.data(qsf + 379);
    const auto *qsf_380 = buffer.data(qsf + 380);
    const auto *qsf_382 = buffer.data(qsf + 382);
    const auto *qsf_383 = buffer.data(qsf + 383);
    const auto *qsf_385 = buffer.data(qsf + 385);
    const auto *qsf_386 = buffer.data(qsf + 386);
    const auto *qsf_387 = buffer.data(qsf + 387);
    const auto *qsf_388 = buffer.data(qsf + 388);
    const auto *qsf_389 = buffer.data(qsf + 389);
    const auto *qsf_390 = buffer.data(qsf + 390);
    const auto *qsf_392 = buffer.data(qsf + 392);
    const auto *qsf_393 = buffer.data(qsf + 393);
    const auto *qsf_395 = buffer.data(qsf + 395);
    const auto *qsf_396 = buffer.data(qsf + 396);
    const auto *qsf_397 = buffer.data(qsf + 397);
    const auto *qsf_398 = buffer.data(qsf + 398);
    const auto *qsf_399 = buffer.data(qsf + 399);
    const auto *qsf_400 = buffer.data(qsf + 400);

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, osf_325, osf_326, osf_327, osf_328, \
                         qsd0_197, qsd1_197, qsf_325, qsf_326, qsf_327, \
                         qsf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_18 * osf_325[k]
                   + f_4 * qsd0_197[k]
                   - f_5 * qsd1_197[k]
                   + f_3 * pc_x[k] * qsf_325[k];

        t_486[k] = f_18 * osf_326[k]
                   + f_3 * pc_x[k] * qsf_326[k];

        t_487[k] = f_18 * osf_327[k]
                   + f_3 * pc_x[k] * qsf_327[k];

        t_488[k] = f_18 * osf_328[k]
                   + f_3 * pc_x[k] * qsf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, osf_246, osf_256, osf_329, \
                         qsd0_195, qsd1_195, qsf_326, qsf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_18 * osf_329[k]
                   + f_3 * pc_x[k] * qsf_329[k];

        t_490[k] = f_14 * osf_256[k]
                   + f_1 * qsd0_195[k]
                   - f_2 * qsd1_195[k]
                   + f_3 * pc_y[k] * qsf_326[k];

        t_491[k] = f_16 * osf_246[k]
                   + f_3 * pc_z[k] * qsf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pc_y, pc_z, osf_249, osf_258, osf_259, qsd0_197, \
                         qsd1_197, qsf_328, qsf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * osf_258[k]
                   + f_4 * qsd0_197[k]
                   - f_5 * qsd1_197[k]
                   + f_3 * pc_y[k] * qsf_328[k];

        t_493[k] = f_14 * osf_259[k]
                   + f_3 * pc_y[k] * qsf_329[k];

        t_494[k] = f_16 * osf_249[k]
                   + f_1 * qsd0_197[k]
                   - f_2 * qsd1_197[k]
                   + f_3 * pc_z[k] * qsf_329[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pc_x, pc_y, pc_z, osf_250, osf_260, osf_330, \
                         qsd0_198, qsd1_198, qsf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_18 * osf_330[k]
                   + f_1 * qsd0_198[k]
                   - f_2 * qsd1_198[k]
                   + f_3 * pc_x[k] * qsf_330[k];

        t_496[k] = f_8 * osf_260[k]
                   + f_3 * pc_y[k] * qsf_330[k];

        t_497[k] = f_18 * osf_250[k]
                   + f_3 * pc_z[k] * qsf_330[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_x, pc_y, osf_262, osf_333, osf_335, qsd0_201, \
                         qsd0_203, qsd1_201, qsd1_203, qsf_332, qsf_333, \
                         qsf_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_18 * osf_333[k]
                   + f_4 * qsd0_201[k]
                   - f_5 * qsd1_201[k]
                   + f_3 * pc_x[k] * qsf_333[k];

        t_499[k] = f_8 * osf_262[k]
                   + f_3 * pc_y[k] * qsf_332[k];

        t_500[k] = f_18 * osf_335[k]
                   + f_4 * qsd0_203[k]
                   - f_5 * qsd1_203[k]
                   + f_3 * pc_x[k] * qsf_335[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pc_x, osf_336, osf_337, osf_338, osf_339, \
                         qsf_336, qsf_337, qsf_338, qsf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_18 * osf_336[k]
                   + f_3 * pc_x[k] * qsf_336[k];

        t_502[k] = f_18 * osf_337[k]
                   + f_3 * pc_x[k] * qsf_337[k];

        t_503[k] = f_18 * osf_338[k]
                   + f_3 * pc_x[k] * qsf_338[k];

        t_504[k] = f_18 * osf_339[k]
                   + f_3 * pc_x[k] * qsf_339[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pc_y, pc_z, osf_256, osf_266, osf_268, qsd0_201, \
                         qsd0_203, qsd1_201, qsd1_203, qsf_336, \
                         qsf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_8 * osf_266[k]
                   + f_1 * qsd0_201[k]
                   - f_2 * qsd1_201[k]
                   + f_3 * pc_y[k] * qsf_336[k];

        t_506[k] = f_18 * osf_256[k]
                   + f_3 * pc_z[k] * qsf_336[k];

        t_507[k] = f_8 * osf_268[k]
                   + f_4 * qsd0_203[k]
                   - f_5 * qsd1_203[k]
                   + f_3 * pc_y[k] * qsf_338[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_y, pc_y, pc_z, osg0_405, osf_259, \
                         osf_269, osf_270, osg1_405, qsd0_203, qsd1_203, qsf_339, \
                         qsf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * osf_269[k]
                   + f_3 * pc_y[k] * qsf_339[k];

        t_509[k] = f_18 * osf_259[k]
                   + f_1 * qsd0_203[k]
                   - f_2 * qsd1_203[k]
                   + f_3 * pc_z[k] * qsf_339[k];

        t_510[k] = pa_y[k] * osg0_405[k]
                   - f_6 * pc_y[k] * osg1_405[k];

        t_511[k] = f_7 * osf_270[k]
                   + f_3 * pc_y[k] * qsf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_y, pc_y, pc_z, osg0_408, osg0_410, \
                         osf_260, osf_271, osf_272, osg1_408, osg1_410, qsf_340, \
                         qsf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_19 * osf_260[k]
                   + f_3 * pc_z[k] * qsf_340[k];

        t_513[k] = pa_y[k] * osg0_408[k]
                   + f_8 * osf_271[k]
                   - f_6 * pc_y[k] * osg1_408[k];

        t_514[k] = f_7 * osf_272[k]
                   + f_3 * pc_y[k] * qsf_342[k];

        t_515[k] = pa_y[k] * osg0_410[k]
                   - f_6 * pc_y[k] * osg1_410[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, osf_346, osf_347, osf_348, osf_349, \
                         qsf_346, qsf_347, qsf_348, qsf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_18 * osf_346[k]
                   + f_3 * pc_x[k] * qsf_346[k];

        t_517[k] = f_18 * osf_347[k]
                   + f_3 * pc_x[k] * qsf_347[k];

        t_518[k] = f_18 * osf_348[k]
                   + f_3 * pc_x[k] * qsf_348[k];

        t_519[k] = f_18 * osf_349[k]
                   + f_3 * pc_x[k] * qsf_349[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, osf_266, osf_276, osf_278, qsd0_207, \
                         qsd0_209, qsd1_207, qsd1_209, qsf_346, \
                         qsf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_7 * osf_276[k]
                   + f_1 * qsd0_207[k]
                   - f_2 * qsd1_207[k]
                   + f_3 * pc_y[k] * qsf_346[k];

        t_521[k] = f_19 * osf_266[k]
                   + f_3 * pc_z[k] * qsf_346[k];

        t_522[k] = f_7 * osf_278[k]
                   + f_4 * qsd0_209[k]
                   - f_5 * qsd1_209[k]
                   + f_3 * pc_y[k] * qsf_348[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_y, pc_x, pc_y, osg0_419, osf_279, \
                         osf_350, osg1_419, qsd0_210, qsd1_210, qsf_349, \
                         qsf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * osf_279[k]
                   + f_3 * pc_y[k] * qsf_349[k];

        t_524[k] = pa_y[k] * osg0_419[k]
                   - f_6 * pc_y[k] * osg1_419[k];

        t_525[k] = f_18 * osf_350[k]
                   + f_1 * qsd0_210[k]
                   - f_2 * qsd1_210[k]
                   + f_3 * pc_x[k] * qsf_350[k];

        t_526[k] = f_3 * pc_y[k] * qsf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_y, pc_z, osf_270, qsd0_210, qsd1_210, \
                         qsf_350, qsf_351, qsf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_17 * osf_270[k]
                   + f_3 * pc_z[k] * qsf_350[k];

        t_528[k] = f_4 * qsd0_210[k]
                   - f_5 * qsd1_210[k]
                   + f_3 * pc_y[k] * qsf_351[k];

        t_529[k] = f_3 * pc_y[k] * qsf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, osf_355, osf_356, osf_357, \
                         qsd0_215, qsd1_215, qsf_355, qsf_356, \
                         qsf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_18 * osf_355[k]
                   + f_4 * qsd0_215[k]
                   - f_5 * qsd1_215[k]
                   + f_3 * pc_x[k] * qsf_355[k];

        t_531[k] = f_18 * osf_356[k]
                   + f_3 * pc_x[k] * qsf_356[k];

        t_532[k] = f_18 * osf_357[k]
                   + f_3 * pc_x[k] * qsf_357[k];

        t_533[k] = f_3 * pc_y[k] * qsf_355[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_x, pc_y, osf_359, qsd0_213, qsd0_214, \
                         qsd1_213, qsd1_214, qsf_356, qsf_357, \
                         qsf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_18 * osf_359[k]
                   + f_3 * pc_x[k] * qsf_359[k];

        t_535[k] = f_1 * qsd0_213[k]
                   - f_2 * qsd1_213[k]
                   + f_3 * pc_y[k] * qsf_356[k];

        t_536[k] = f_10 * qsd0_214[k]
                   - f_11 * qsd1_214[k]
                   + f_3 * pc_y[k] * qsf_357[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, osf_279, osf_360, \
                         qsd0_215, qsd0_216, qsd1_215, qsd1_216, qsf_358, qsf_359, \
                         qsf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * qsd0_215[k]
                   - f_5 * qsd1_215[k]
                   + f_3 * pc_y[k] * qsf_358[k];

        t_538[k] = f_3 * pc_y[k] * qsf_359[k];

        t_539[k] = f_17 * osf_279[k]
                   + f_1 * qsd0_215[k]
                   - f_2 * qsd1_215[k]
                   + f_3 * pc_z[k] * qsf_359[k];

        t_540[k] = f_16 * osf_360[k]
                   + f_1 * qsd0_216[k]
                   - f_2 * qsd1_216[k]
                   + f_3 * pc_x[k] * qsf_360[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, osf_280, osf_363, \
                         qsd0_219, qsd1_219, qsf_360, qsf_361, \
                         qsf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_15 * osf_280[k]
                   + f_3 * pc_y[k] * qsf_360[k];

        t_542[k] = f_3 * pc_z[k] * qsf_360[k];

        t_543[k] = f_16 * osf_363[k]
                   + f_4 * qsd0_219[k]
                   - f_5 * qsd1_219[k]
                   + f_3 * pc_x[k] * qsf_363[k];

        t_544[k] = f_3 * pc_z[k] * qsf_361[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, pc_z, osf_366, osf_368, qsd0_216, \
                         qsd1_216, qsf_362, qsf_363, qsf_366, qsf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * qsd0_216[k]
                   - f_5 * qsd1_216[k]
                   + f_3 * pc_z[k] * qsf_362[k];

        t_546[k] = f_16 * osf_366[k]
                   + f_3 * pc_x[k] * qsf_366[k];

        t_547[k] = f_3 * pc_z[k] * qsf_363[k];

        t_548[k] = f_16 * osf_368[k]
                   + f_3 * pc_x[k] * qsf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, t_553, pc_x, pc_y, pc_z, osf_286, \
                         osf_289, osf_369, qsd0_219, qsd1_219, qsf_366, qsf_367, \
                         qsf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_16 * osf_369[k]
                   + f_3 * pc_x[k] * qsf_369[k];

        t_550[k] = f_15 * osf_286[k]
                   + f_1 * qsd0_219[k]
                   - f_2 * qsd1_219[k]
                   + f_3 * pc_y[k] * qsf_366[k];

        t_551[k] = f_3 * pc_z[k] * qsf_366[k];

        t_552[k] = f_4 * qsd0_219[k]
                   - f_5 * qsd1_219[k]
                   + f_3 * pc_z[k] * qsf_367[k];

        t_553[k] = f_15 * osf_289[k]
                   + f_3 * pc_y[k] * qsf_369[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, pa_z, pc_y, pc_z, osg0_420, osf_280, \
                         osf_290, osg1_420, qsd0_221, qsd1_221, qsf_369, \
                         qsf_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_1 * qsd0_221[k]
                   - f_2 * qsd1_221[k]
                   + f_3 * pc_z[k] * qsf_369[k];

        t_555[k] = pa_z[k] * osg0_420[k]
                   - f_6 * pc_z[k] * osg1_420[k];

        t_556[k] = f_17 * osf_290[k]
                   + f_3 * pc_y[k] * qsf_370[k];

        t_557[k] = f_7 * osf_280[k]
                   + f_3 * pc_z[k] * qsf_370[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pa_z, pc_x, pc_y, pc_z, osg0_423, osf_292, \
                         osf_375, osg1_423, qsd0_227, qsd1_227, qsf_372, \
                         qsf_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = pa_z[k] * osg0_423[k]
                   - f_6 * pc_z[k] * osg1_423[k];

        t_559[k] = f_17 * osf_292[k]
                   + f_3 * pc_y[k] * qsf_372[k];

        t_560[k] = f_16 * osf_375[k]
                   + f_4 * qsd0_227[k]
                   - f_5 * qsd1_227[k]
                   + f_3 * pc_x[k] * qsf_375[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, pc_x, osf_376, osf_377, osf_378, osf_379, \
                         qsf_376, qsf_377, qsf_378, qsf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_16 * osf_376[k]
                   + f_3 * pc_x[k] * qsf_376[k];

        t_562[k] = f_16 * osf_377[k]
                   + f_3 * pc_x[k] * qsf_377[k];

        t_563[k] = f_16 * osf_378[k]
                   + f_3 * pc_x[k] * qsf_378[k];

        t_564[k] = f_16 * osf_379[k]
                   + f_3 * pc_x[k] * qsf_379[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, pa_z, pc_y, pc_z, osg0_430, osf_286, osf_298, \
                         osg1_430, qsd0_227, qsd1_227, qsf_376, \
                         qsf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pa_z[k] * osg0_430[k]
                   - f_6 * pc_z[k] * osg1_430[k];

        t_566[k] = f_7 * osf_286[k]
                   + f_3 * pc_z[k] * qsf_376[k];

        t_567[k] = f_17 * osf_298[k]
                   + f_4 * qsd0_227[k]
                   - f_5 * qsd1_227[k]
                   + f_3 * pc_y[k] * qsf_378[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_x, pc_y, pc_z, osf_289, osf_299, osf_380, \
                         qsd0_227, qsd0_228, qsd1_227, qsd1_228, qsf_379, \
                         qsf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_17 * osf_299[k]
                   + f_3 * pc_y[k] * qsf_379[k];

        t_569[k] = f_7 * osf_289[k]
                   + f_1 * qsd0_227[k]
                   - f_2 * qsd1_227[k]
                   + f_3 * pc_z[k] * qsf_379[k];

        t_570[k] = f_16 * osf_380[k]
                   + f_1 * qsd0_228[k]
                   - f_2 * qsd1_228[k]
                   + f_3 * pc_x[k] * qsf_380[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, pc_x, pc_y, pc_z, osf_290, osf_300, \
                         osf_302, osf_383, qsd0_231, qsd1_231, qsf_380, qsf_382, \
                         qsf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_19 * osf_300[k]
                   + f_3 * pc_y[k] * qsf_380[k];

        t_572[k] = f_8 * osf_290[k]
                   + f_3 * pc_z[k] * qsf_380[k];

        t_573[k] = f_16 * osf_383[k]
                   + f_4 * qsd0_231[k]
                   - f_5 * qsd1_231[k]
                   + f_3 * pc_x[k] * qsf_383[k];

        t_574[k] = f_19 * osf_302[k]
                   + f_3 * pc_y[k] * qsf_382[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pc_x, osf_385, osf_386, osf_387, osf_388, \
                         qsd0_233, qsd1_233, qsf_385, qsf_386, qsf_387, \
                         qsf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_16 * osf_385[k]
                   + f_4 * qsd0_233[k]
                   - f_5 * qsd1_233[k]
                   + f_3 * pc_x[k] * qsf_385[k];

        t_576[k] = f_16 * osf_386[k]
                   + f_3 * pc_x[k] * qsf_386[k];

        t_577[k] = f_16 * osf_387[k]
                   + f_3 * pc_x[k] * qsf_387[k];

        t_578[k] = f_16 * osf_388[k]
                   + f_3 * pc_x[k] * qsf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pc_x, pc_y, pc_z, osf_296, osf_306, osf_389, \
                         qsd0_231, qsd1_231, qsf_386, qsf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_16 * osf_389[k]
                   + f_3 * pc_x[k] * qsf_389[k];

        t_580[k] = f_19 * osf_306[k]
                   + f_1 * qsd0_231[k]
                   - f_2 * qsd1_231[k]
                   + f_3 * pc_y[k] * qsf_386[k];

        t_581[k] = f_8 * osf_296[k]
                   + f_3 * pc_z[k] * qsf_386[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, pc_z, osf_299, osf_308, osf_309, qsd0_233, \
                         qsd1_233, qsf_388, qsf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_19 * osf_308[k]
                   + f_4 * qsd0_233[k]
                   - f_5 * qsd1_233[k]
                   + f_3 * pc_y[k] * qsf_388[k];

        t_583[k] = f_19 * osf_309[k]
                   + f_3 * pc_y[k] * qsf_389[k];

        t_584[k] = f_8 * osf_299[k]
                   + f_1 * qsd0_233[k]
                   - f_2 * qsd1_233[k]
                   + f_3 * pc_z[k] * qsf_389[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pc_x, pc_y, pc_z, osf_300, osf_310, osf_390, \
                         qsd0_234, qsd1_234, qsf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_16 * osf_390[k]
                   + f_1 * qsd0_234[k]
                   - f_2 * qsd1_234[k]
                   + f_3 * pc_x[k] * qsf_390[k];

        t_586[k] = f_18 * osf_310[k]
                   + f_3 * pc_y[k] * qsf_390[k];

        t_587[k] = f_14 * osf_300[k]
                   + f_3 * pc_z[k] * qsf_390[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pc_x, pc_y, osf_312, osf_393, osf_395, qsd0_237, \
                         qsd0_239, qsd1_237, qsd1_239, qsf_392, qsf_393, \
                         qsf_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_16 * osf_393[k]
                   + f_4 * qsd0_237[k]
                   - f_5 * qsd1_237[k]
                   + f_3 * pc_x[k] * qsf_393[k];

        t_589[k] = f_18 * osf_312[k]
                   + f_3 * pc_y[k] * qsf_392[k];

        t_590[k] = f_16 * osf_395[k]
                   + f_4 * qsd0_239[k]
                   - f_5 * qsd1_239[k]
                   + f_3 * pc_x[k] * qsf_395[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pc_x, osf_396, osf_397, osf_398, osf_399, \
                         qsf_396, qsf_397, qsf_398, qsf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_16 * osf_396[k]
                   + f_3 * pc_x[k] * qsf_396[k];

        t_592[k] = f_16 * osf_397[k]
                   + f_3 * pc_x[k] * qsf_397[k];

        t_593[k] = f_16 * osf_398[k]
                   + f_3 * pc_x[k] * qsf_398[k];

        t_594[k] = f_16 * osf_399[k]
                   + f_3 * pc_x[k] * qsf_399[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_y, pc_z, osf_306, osf_316, osf_318, qsd0_237, \
                         qsd0_239, qsd1_237, qsd1_239, qsf_396, \
                         qsf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_18 * osf_316[k]
                   + f_1 * qsd0_237[k]
                   - f_2 * qsd1_237[k]
                   + f_3 * pc_y[k] * qsf_396[k];

        t_596[k] = f_14 * osf_306[k]
                   + f_3 * pc_z[k] * qsf_396[k];

        t_597[k] = f_18 * osf_318[k]
                   + f_4 * qsd0_239[k]
                   - f_5 * qsd1_239[k]
                   + f_3 * pc_y[k] * qsf_398[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_y, pc_z, osf_309, osf_319, osf_400, \
                         qsd0_239, qsd0_240, qsd1_239, qsd1_240, qsf_399, \
                         qsf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_18 * osf_319[k]
                   + f_3 * pc_y[k] * qsf_399[k];

        t_599[k] = f_14 * osf_309[k]
                   + f_1 * qsd0_239[k]
                   - f_2 * qsd1_239[k]
                   + f_3 * pc_z[k] * qsf_399[k];

        t_600[k] = f_16 * osf_400[k]
                   + f_1 * qsd0_240[k]
                   - f_2 * qsd1_240[k]
                   + f_3 * pc_x[k] * qsf_400[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_525 = buffer.data(osg0 + 525);
    const auto *osg0_528 = buffer.data(osg0 + 528);
    const auto *osg0_530 = buffer.data(osg0 + 530);
    const auto *osg0_539 = buffer.data(osg0 + 539);
    const auto *osg0_540 = buffer.data(osg0 + 540);
    const auto *osg0_543 = buffer.data(osg0 + 543);
    const auto *osg0_550 = buffer.data(osg0 + 550);

    const auto *osf_310 = buffer.data(osf + 310);
    const auto *osf_316 = buffer.data(osf + 316);
    const auto *osf_319 = buffer.data(osf + 319);
    const auto *osf_320 = buffer.data(osf + 320);
    const auto *osf_322 = buffer.data(osf + 322);
    const auto *osf_326 = buffer.data(osf + 326);
    const auto *osf_328 = buffer.data(osf + 328);
    const auto *osf_329 = buffer.data(osf + 329);
    const auto *osf_330 = buffer.data(osf + 330);
    const auto *osf_332 = buffer.data(osf + 332);
    const auto *osf_336 = buffer.data(osf + 336);
    const auto *osf_338 = buffer.data(osf + 338);
    const auto *osf_339 = buffer.data(osf + 339);
    const auto *osf_340 = buffer.data(osf + 340);
    const auto *osf_342 = buffer.data(osf + 342);
    const auto *osf_346 = buffer.data(osf + 346);
    const auto *osf_348 = buffer.data(osf + 348);
    const auto *osf_349 = buffer.data(osf + 349);
    const auto *osf_350 = buffer.data(osf + 350);
    const auto *osf_351 = buffer.data(osf + 351);
    const auto *osf_352 = buffer.data(osf + 352);
    const auto *osf_356 = buffer.data(osf + 356);
    const auto *osf_358 = buffer.data(osf + 358);
    const auto *osf_359 = buffer.data(osf + 359);
    const auto *osf_360 = buffer.data(osf + 360);
    const auto *osf_366 = buffer.data(osf + 366);
    const auto *osf_369 = buffer.data(osf + 369);
    const auto *osf_370 = buffer.data(osf + 370);
    const auto *osf_372 = buffer.data(osf + 372);
    const auto *osf_376 = buffer.data(osf + 376);
    const auto *osf_378 = buffer.data(osf + 378);
    const auto *osf_379 = buffer.data(osf + 379);
    const auto *osf_380 = buffer.data(osf + 380);
    const auto *osf_382 = buffer.data(osf + 382);
    const auto *osf_386 = buffer.data(osf + 386);
    const auto *osf_388 = buffer.data(osf + 388);
    const auto *osf_389 = buffer.data(osf + 389);
    const auto *osf_403 = buffer.data(osf + 403);
    const auto *osf_405 = buffer.data(osf + 405);
    const auto *osf_406 = buffer.data(osf + 406);
    const auto *osf_407 = buffer.data(osf + 407);
    const auto *osf_408 = buffer.data(osf + 408);
    const auto *osf_409 = buffer.data(osf + 409);
    const auto *osf_410 = buffer.data(osf + 410);
    const auto *osf_413 = buffer.data(osf + 413);
    const auto *osf_415 = buffer.data(osf + 415);
    const auto *osf_416 = buffer.data(osf + 416);
    const auto *osf_417 = buffer.data(osf + 417);
    const auto *osf_418 = buffer.data(osf + 418);
    const auto *osf_419 = buffer.data(osf + 419);
    const auto *osf_420 = buffer.data(osf + 420);
    const auto *osf_423 = buffer.data(osf + 423);
    const auto *osf_425 = buffer.data(osf + 425);
    const auto *osf_426 = buffer.data(osf + 426);
    const auto *osf_427 = buffer.data(osf + 427);
    const auto *osf_428 = buffer.data(osf + 428);
    const auto *osf_429 = buffer.data(osf + 429);
    const auto *osf_436 = buffer.data(osf + 436);
    const auto *osf_437 = buffer.data(osf + 437);
    const auto *osf_438 = buffer.data(osf + 438);
    const auto *osf_439 = buffer.data(osf + 439);
    const auto *osf_440 = buffer.data(osf + 440);
    const auto *osf_445 = buffer.data(osf + 445);
    const auto *osf_446 = buffer.data(osf + 446);
    const auto *osf_447 = buffer.data(osf + 447);
    const auto *osf_449 = buffer.data(osf + 449);
    const auto *osf_450 = buffer.data(osf + 450);
    const auto *osf_453 = buffer.data(osf + 453);
    const auto *osf_456 = buffer.data(osf + 456);
    const auto *osf_458 = buffer.data(osf + 458);
    const auto *osf_459 = buffer.data(osf + 459);
    const auto *osf_465 = buffer.data(osf + 465);
    const auto *osf_466 = buffer.data(osf + 466);
    const auto *osf_467 = buffer.data(osf + 467);
    const auto *osf_468 = buffer.data(osf + 468);
    const auto *osf_469 = buffer.data(osf + 469);
    const auto *osf_470 = buffer.data(osf + 470);
    const auto *osf_473 = buffer.data(osf + 473);
    const auto *osf_475 = buffer.data(osf + 475);
    const auto *osf_476 = buffer.data(osf + 476);
    const auto *osf_477 = buffer.data(osf + 477);
    const auto *osf_478 = buffer.data(osf + 478);
    const auto *osf_479 = buffer.data(osf + 479);

    const auto *osg1_525 = buffer.data(osg1 + 525);
    const auto *osg1_528 = buffer.data(osg1 + 528);
    const auto *osg1_530 = buffer.data(osg1 + 530);
    const auto *osg1_539 = buffer.data(osg1 + 539);
    const auto *osg1_540 = buffer.data(osg1 + 540);
    const auto *osg1_543 = buffer.data(osg1 + 543);
    const auto *osg1_550 = buffer.data(osg1 + 550);

    const auto *qsd0_243 = buffer.data(qsd0 + 243);
    const auto *qsd0_245 = buffer.data(qsd0 + 245);
    const auto *qsd0_246 = buffer.data(qsd0 + 246);
    const auto *qsd0_249 = buffer.data(qsd0 + 249);
    const auto *qsd0_251 = buffer.data(qsd0 + 251);
    const auto *qsd0_252 = buffer.data(qsd0 + 252);
    const auto *qsd0_255 = buffer.data(qsd0 + 255);
    const auto *qsd0_257 = buffer.data(qsd0 + 257);
    const auto *qsd0_261 = buffer.data(qsd0 + 261);
    const auto *qsd0_263 = buffer.data(qsd0 + 263);
    const auto *qsd0_264 = buffer.data(qsd0 + 264);
    const auto *qsd0_267 = buffer.data(qsd0 + 267);
    const auto *qsd0_268 = buffer.data(qsd0 + 268);
    const auto *qsd0_269 = buffer.data(qsd0 + 269);
    const auto *qsd0_270 = buffer.data(qsd0 + 270);
    const auto *qsd0_273 = buffer.data(qsd0 + 273);
    const auto *qsd0_275 = buffer.data(qsd0 + 275);
    const auto *qsd0_281 = buffer.data(qsd0 + 281);
    const auto *qsd0_282 = buffer.data(qsd0 + 282);
    const auto *qsd0_285 = buffer.data(qsd0 + 285);
    const auto *qsd0_287 = buffer.data(qsd0 + 287);

    const auto *qsd1_243 = buffer.data(qsd1 + 243);
    const auto *qsd1_245 = buffer.data(qsd1 + 245);
    const auto *qsd1_246 = buffer.data(qsd1 + 246);
    const auto *qsd1_249 = buffer.data(qsd1 + 249);
    const auto *qsd1_251 = buffer.data(qsd1 + 251);
    const auto *qsd1_252 = buffer.data(qsd1 + 252);
    const auto *qsd1_255 = buffer.data(qsd1 + 255);
    const auto *qsd1_257 = buffer.data(qsd1 + 257);
    const auto *qsd1_261 = buffer.data(qsd1 + 261);
    const auto *qsd1_263 = buffer.data(qsd1 + 263);
    const auto *qsd1_264 = buffer.data(qsd1 + 264);
    const auto *qsd1_267 = buffer.data(qsd1 + 267);
    const auto *qsd1_268 = buffer.data(qsd1 + 268);
    const auto *qsd1_269 = buffer.data(qsd1 + 269);
    const auto *qsd1_270 = buffer.data(qsd1 + 270);
    const auto *qsd1_273 = buffer.data(qsd1 + 273);
    const auto *qsd1_275 = buffer.data(qsd1 + 275);
    const auto *qsd1_281 = buffer.data(qsd1 + 281);
    const auto *qsd1_282 = buffer.data(qsd1 + 282);
    const auto *qsd1_285 = buffer.data(qsd1 + 285);
    const auto *qsd1_287 = buffer.data(qsd1 + 287);

    const auto *qsf_400 = buffer.data(qsf + 400);
    const auto *qsf_402 = buffer.data(qsf + 402);
    const auto *qsf_403 = buffer.data(qsf + 403);
    const auto *qsf_405 = buffer.data(qsf + 405);
    const auto *qsf_406 = buffer.data(qsf + 406);
    const auto *qsf_407 = buffer.data(qsf + 407);
    const auto *qsf_408 = buffer.data(qsf + 408);
    const auto *qsf_409 = buffer.data(qsf + 409);
    const auto *qsf_410 = buffer.data(qsf + 410);
    const auto *qsf_412 = buffer.data(qsf + 412);
    const auto *qsf_413 = buffer.data(qsf + 413);
    const auto *qsf_415 = buffer.data(qsf + 415);
    const auto *qsf_416 = buffer.data(qsf + 416);
    const auto *qsf_417 = buffer.data(qsf + 417);
    const auto *qsf_418 = buffer.data(qsf + 418);
    const auto *qsf_419 = buffer.data(qsf + 419);
    const auto *qsf_420 = buffer.data(qsf + 420);
    const auto *qsf_422 = buffer.data(qsf + 422);
    const auto *qsf_423 = buffer.data(qsf + 423);
    const auto *qsf_425 = buffer.data(qsf + 425);
    const auto *qsf_426 = buffer.data(qsf + 426);
    const auto *qsf_427 = buffer.data(qsf + 427);
    const auto *qsf_428 = buffer.data(qsf + 428);
    const auto *qsf_429 = buffer.data(qsf + 429);
    const auto *qsf_430 = buffer.data(qsf + 430);
    const auto *qsf_432 = buffer.data(qsf + 432);
    const auto *qsf_436 = buffer.data(qsf + 436);
    const auto *qsf_437 = buffer.data(qsf + 437);
    const auto *qsf_438 = buffer.data(qsf + 438);
    const auto *qsf_439 = buffer.data(qsf + 439);
    const auto *qsf_440 = buffer.data(qsf + 440);
    const auto *qsf_441 = buffer.data(qsf + 441);
    const auto *qsf_442 = buffer.data(qsf + 442);
    const auto *qsf_445 = buffer.data(qsf + 445);
    const auto *qsf_446 = buffer.data(qsf + 446);
    const auto *qsf_447 = buffer.data(qsf + 447);
    const auto *qsf_448 = buffer.data(qsf + 448);
    const auto *qsf_449 = buffer.data(qsf + 449);
    const auto *qsf_450 = buffer.data(qsf + 450);
    const auto *qsf_451 = buffer.data(qsf + 451);
    const auto *qsf_452 = buffer.data(qsf + 452);
    const auto *qsf_453 = buffer.data(qsf + 453);
    const auto *qsf_456 = buffer.data(qsf + 456);
    const auto *qsf_457 = buffer.data(qsf + 457);
    const auto *qsf_458 = buffer.data(qsf + 458);
    const auto *qsf_459 = buffer.data(qsf + 459);
    const auto *qsf_460 = buffer.data(qsf + 460);
    const auto *qsf_462 = buffer.data(qsf + 462);
    const auto *qsf_465 = buffer.data(qsf + 465);
    const auto *qsf_466 = buffer.data(qsf + 466);
    const auto *qsf_467 = buffer.data(qsf + 467);
    const auto *qsf_468 = buffer.data(qsf + 468);
    const auto *qsf_469 = buffer.data(qsf + 469);
    const auto *qsf_470 = buffer.data(qsf + 470);
    const auto *qsf_472 = buffer.data(qsf + 472);
    const auto *qsf_473 = buffer.data(qsf + 473);
    const auto *qsf_475 = buffer.data(qsf + 475);
    const auto *qsf_476 = buffer.data(qsf + 476);
    const auto *qsf_477 = buffer.data(qsf + 477);
    const auto *qsf_478 = buffer.data(qsf + 478);
    const auto *qsf_479 = buffer.data(qsf + 479);

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, pc_z, osf_310, osf_320, \
                         osf_322, osf_403, qsd0_243, qsd1_243, qsf_400, qsf_402, \
                         qsf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_16 * osf_320[k]
                   + f_3 * pc_y[k] * qsf_400[k];

        t_602[k] = f_16 * osf_310[k]
                   + f_3 * pc_z[k] * qsf_400[k];

        t_603[k] = f_16 * osf_403[k]
                   + f_4 * qsd0_243[k]
                   - f_5 * qsd1_243[k]
                   + f_3 * pc_x[k] * qsf_403[k];

        t_604[k] = f_16 * osf_322[k]
                   + f_3 * pc_y[k] * qsf_402[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, osf_405, osf_406, osf_407, osf_408, \
                         qsd0_245, qsd1_245, qsf_405, qsf_406, qsf_407, \
                         qsf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_16 * osf_405[k]
                   + f_4 * qsd0_245[k]
                   - f_5 * qsd1_245[k]
                   + f_3 * pc_x[k] * qsf_405[k];

        t_606[k] = f_16 * osf_406[k]
                   + f_3 * pc_x[k] * qsf_406[k];

        t_607[k] = f_16 * osf_407[k]
                   + f_3 * pc_x[k] * qsf_407[k];

        t_608[k] = f_16 * osf_408[k]
                   + f_3 * pc_x[k] * qsf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_x, pc_y, pc_z, osf_316, osf_326, osf_409, \
                         qsd0_243, qsd1_243, qsf_406, qsf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_16 * osf_409[k]
                   + f_3 * pc_x[k] * qsf_409[k];

        t_610[k] = f_16 * osf_326[k]
                   + f_1 * qsd0_243[k]
                   - f_2 * qsd1_243[k]
                   + f_3 * pc_y[k] * qsf_406[k];

        t_611[k] = f_16 * osf_316[k]
                   + f_3 * pc_z[k] * qsf_406[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pc_y, pc_z, osf_319, osf_328, osf_329, qsd0_245, \
                         qsd1_245, qsf_408, qsf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_16 * osf_328[k]
                   + f_4 * qsd0_245[k]
                   - f_5 * qsd1_245[k]
                   + f_3 * pc_y[k] * qsf_408[k];

        t_613[k] = f_16 * osf_329[k]
                   + f_3 * pc_y[k] * qsf_409[k];

        t_614[k] = f_16 * osf_319[k]
                   + f_1 * qsd0_245[k]
                   - f_2 * qsd1_245[k]
                   + f_3 * pc_z[k] * qsf_409[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, pc_y, pc_z, osf_320, osf_330, osf_410, \
                         qsd0_246, qsd1_246, qsf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_16 * osf_410[k]
                   + f_1 * qsd0_246[k]
                   - f_2 * qsd1_246[k]
                   + f_3 * pc_x[k] * qsf_410[k];

        t_616[k] = f_14 * osf_330[k]
                   + f_3 * pc_y[k] * qsf_410[k];

        t_617[k] = f_18 * osf_320[k]
                   + f_3 * pc_z[k] * qsf_410[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pc_x, pc_y, osf_332, osf_413, osf_415, qsd0_249, \
                         qsd0_251, qsd1_249, qsd1_251, qsf_412, qsf_413, \
                         qsf_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_16 * osf_413[k]
                   + f_4 * qsd0_249[k]
                   - f_5 * qsd1_249[k]
                   + f_3 * pc_x[k] * qsf_413[k];

        t_619[k] = f_14 * osf_332[k]
                   + f_3 * pc_y[k] * qsf_412[k];

        t_620[k] = f_16 * osf_415[k]
                   + f_4 * qsd0_251[k]
                   - f_5 * qsd1_251[k]
                   + f_3 * pc_x[k] * qsf_415[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pc_x, osf_416, osf_417, osf_418, osf_419, \
                         qsf_416, qsf_417, qsf_418, qsf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_16 * osf_416[k]
                   + f_3 * pc_x[k] * qsf_416[k];

        t_622[k] = f_16 * osf_417[k]
                   + f_3 * pc_x[k] * qsf_417[k];

        t_623[k] = f_16 * osf_418[k]
                   + f_3 * pc_x[k] * qsf_418[k];

        t_624[k] = f_16 * osf_419[k]
                   + f_3 * pc_x[k] * qsf_419[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pc_z, osf_326, osf_336, osf_338, qsd0_249, \
                         qsd0_251, qsd1_249, qsd1_251, qsf_416, \
                         qsf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_14 * osf_336[k]
                   + f_1 * qsd0_249[k]
                   - f_2 * qsd1_249[k]
                   + f_3 * pc_y[k] * qsf_416[k];

        t_626[k] = f_18 * osf_326[k]
                   + f_3 * pc_z[k] * qsf_416[k];

        t_627[k] = f_14 * osf_338[k]
                   + f_4 * qsd0_251[k]
                   - f_5 * qsd1_251[k]
                   + f_3 * pc_y[k] * qsf_418[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pc_z, osf_329, osf_339, osf_420, \
                         qsd0_251, qsd0_252, qsd1_251, qsd1_252, qsf_419, \
                         qsf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_14 * osf_339[k]
                   + f_3 * pc_y[k] * qsf_419[k];

        t_629[k] = f_18 * osf_329[k]
                   + f_1 * qsd0_251[k]
                   - f_2 * qsd1_251[k]
                   + f_3 * pc_z[k] * qsf_419[k];

        t_630[k] = f_16 * osf_420[k]
                   + f_1 * qsd0_252[k]
                   - f_2 * qsd1_252[k]
                   + f_3 * pc_x[k] * qsf_420[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pc_z, osf_330, osf_340, \
                         osf_342, osf_423, qsd0_255, qsd1_255, qsf_420, qsf_422, \
                         qsf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_8 * osf_340[k]
                   + f_3 * pc_y[k] * qsf_420[k];

        t_632[k] = f_19 * osf_330[k]
                   + f_3 * pc_z[k] * qsf_420[k];

        t_633[k] = f_16 * osf_423[k]
                   + f_4 * qsd0_255[k]
                   - f_5 * qsd1_255[k]
                   + f_3 * pc_x[k] * qsf_423[k];

        t_634[k] = f_8 * osf_342[k]
                   + f_3 * pc_y[k] * qsf_422[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pc_x, osf_425, osf_426, osf_427, osf_428, \
                         qsd0_257, qsd1_257, qsf_425, qsf_426, qsf_427, \
                         qsf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_16 * osf_425[k]
                   + f_4 * qsd0_257[k]
                   - f_5 * qsd1_257[k]
                   + f_3 * pc_x[k] * qsf_425[k];

        t_636[k] = f_16 * osf_426[k]
                   + f_3 * pc_x[k] * qsf_426[k];

        t_637[k] = f_16 * osf_427[k]
                   + f_3 * pc_x[k] * qsf_427[k];

        t_638[k] = f_16 * osf_428[k]
                   + f_3 * pc_x[k] * qsf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, pc_z, osf_336, osf_346, osf_429, \
                         qsd0_255, qsd1_255, qsf_426, qsf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_16 * osf_429[k]
                   + f_3 * pc_x[k] * qsf_429[k];

        t_640[k] = f_8 * osf_346[k]
                   + f_1 * qsd0_255[k]
                   - f_2 * qsd1_255[k]
                   + f_3 * pc_y[k] * qsf_426[k];

        t_641[k] = f_19 * osf_336[k]
                   + f_3 * pc_z[k] * qsf_426[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_y, pc_y, pc_z, osg0_525, osf_339, \
                         osf_348, osf_349, osg1_525, qsd0_257, qsd1_257, qsf_428, \
                         qsf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_8 * osf_348[k]
                   + f_4 * qsd0_257[k]
                   - f_5 * qsd1_257[k]
                   + f_3 * pc_y[k] * qsf_428[k];

        t_643[k] = f_8 * osf_349[k]
                   + f_3 * pc_y[k] * qsf_429[k];

        t_644[k] = f_19 * osf_339[k]
                   + f_1 * qsd0_257[k]
                   - f_2 * qsd1_257[k]
                   + f_3 * pc_z[k] * qsf_429[k];

        t_645[k] = pa_y[k] * osg0_525[k]
                   - f_6 * pc_y[k] * osg1_525[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pa_y, pc_y, pc_z, osg0_528, osf_340, \
                         osf_350, osf_351, osf_352, osg1_528, qsf_430, \
                         qsf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_7 * osf_350[k]
                   + f_3 * pc_y[k] * qsf_430[k];

        t_647[k] = f_17 * osf_340[k]
                   + f_3 * pc_z[k] * qsf_430[k];

        t_648[k] = pa_y[k] * osg0_528[k]
                   + f_8 * osf_351[k]
                   - f_6 * pc_y[k] * osg1_528[k];

        t_649[k] = f_7 * osf_352[k]
                   + f_3 * pc_y[k] * qsf_432[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, pa_y, pc_x, pc_y, osg0_530, osf_436, \
                         osf_437, osf_438, osg1_530, qsf_436, qsf_437, \
                         qsf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pa_y[k] * osg0_530[k]
                   - f_6 * pc_y[k] * osg1_530[k];

        t_651[k] = f_16 * osf_436[k]
                   + f_3 * pc_x[k] * qsf_436[k];

        t_652[k] = f_16 * osf_437[k]
                   + f_3 * pc_x[k] * qsf_437[k];

        t_653[k] = f_16 * osf_438[k]
                   + f_3 * pc_x[k] * qsf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, pc_z, osf_346, osf_356, osf_439, \
                         qsd0_261, qsd1_261, qsf_436, qsf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_16 * osf_439[k]
                   + f_3 * pc_x[k] * qsf_439[k];

        t_655[k] = f_7 * osf_356[k]
                   + f_1 * qsd0_261[k]
                   - f_2 * qsd1_261[k]
                   + f_3 * pc_y[k] * qsf_436[k];

        t_656[k] = f_17 * osf_346[k]
                   + f_3 * pc_z[k] * qsf_436[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_y, pc_y, osg0_539, osf_358, osf_359, \
                         osg1_539, qsd0_263, qsd1_263, qsf_438, \
                         qsf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_7 * osf_358[k]
                   + f_4 * qsd0_263[k]
                   - f_5 * qsd1_263[k]
                   + f_3 * pc_y[k] * qsf_438[k];

        t_658[k] = f_7 * osf_359[k]
                   + f_3 * pc_y[k] * qsf_439[k];

        t_659[k] = pa_y[k] * osg0_539[k]
                   - f_6 * pc_y[k] * osg1_539[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pc_y, pc_z, osf_350, \
                         osf_440, qsd0_264, qsd1_264, qsf_440, qsf_441, \
                         qsf_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_16 * osf_440[k]
                   + f_1 * qsd0_264[k]
                   - f_2 * qsd1_264[k]
                   + f_3 * pc_x[k] * qsf_440[k];

        t_661[k] = f_3 * pc_y[k] * qsf_440[k];

        t_662[k] = f_15 * osf_350[k]
                   + f_3 * pc_z[k] * qsf_440[k];

        t_663[k] = f_4 * qsd0_264[k]
                   - f_5 * qsd1_264[k]
                   + f_3 * pc_y[k] * qsf_441[k];

        t_664[k] = f_3 * pc_y[k] * qsf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pc_x, pc_y, osf_445, osf_446, osf_447, \
                         qsd0_269, qsd1_269, qsf_445, qsf_446, \
                         qsf_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_16 * osf_445[k]
                   + f_4 * qsd0_269[k]
                   - f_5 * qsd1_269[k]
                   + f_3 * pc_x[k] * qsf_445[k];

        t_666[k] = f_16 * osf_446[k]
                   + f_3 * pc_x[k] * qsf_446[k];

        t_667[k] = f_16 * osf_447[k]
                   + f_3 * pc_x[k] * qsf_447[k];

        t_668[k] = f_3 * pc_y[k] * qsf_445[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_x, pc_y, osf_449, qsd0_267, qsd0_268, \
                         qsd1_267, qsd1_268, qsf_446, qsf_447, \
                         qsf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_16 * osf_449[k]
                   + f_3 * pc_x[k] * qsf_449[k];

        t_670[k] = f_1 * qsd0_267[k]
                   - f_2 * qsd1_267[k]
                   + f_3 * pc_y[k] * qsf_446[k];

        t_671[k] = f_10 * qsd0_268[k]
                   - f_11 * qsd1_268[k]
                   + f_3 * pc_y[k] * qsf_447[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, pc_y, pc_z, osf_359, osf_450, \
                         qsd0_269, qsd0_270, qsd1_269, qsd1_270, qsf_448, qsf_449, \
                         qsf_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * qsd0_269[k]
                   - f_5 * qsd1_269[k]
                   + f_3 * pc_y[k] * qsf_448[k];

        t_673[k] = f_3 * pc_y[k] * qsf_449[k];

        t_674[k] = f_15 * osf_359[k]
                   + f_1 * qsd0_269[k]
                   - f_2 * qsd1_269[k]
                   + f_3 * pc_z[k] * qsf_449[k];

        t_675[k] = f_14 * osf_450[k]
                   + f_1 * qsd0_270[k]
                   - f_2 * qsd1_270[k]
                   + f_3 * pc_x[k] * qsf_450[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pc_x, pc_y, pc_z, osf_360, osf_453, \
                         qsd0_273, qsd1_273, qsf_450, qsf_451, \
                         qsf_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_13 * osf_360[k]
                   + f_3 * pc_y[k] * qsf_450[k];

        t_677[k] = f_3 * pc_z[k] * qsf_450[k];

        t_678[k] = f_14 * osf_453[k]
                   + f_4 * qsd0_273[k]
                   - f_5 * qsd1_273[k]
                   + f_3 * pc_x[k] * qsf_453[k];

        t_679[k] = f_3 * pc_z[k] * qsf_451[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_z, osf_456, osf_458, qsd0_270, \
                         qsd1_270, qsf_452, qsf_453, qsf_456, qsf_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * qsd0_270[k]
                   - f_5 * qsd1_270[k]
                   + f_3 * pc_z[k] * qsf_452[k];

        t_681[k] = f_14 * osf_456[k]
                   + f_3 * pc_x[k] * qsf_456[k];

        t_682[k] = f_3 * pc_z[k] * qsf_453[k];

        t_683[k] = f_14 * osf_458[k]
                   + f_3 * pc_x[k] * qsf_458[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, pc_x, pc_y, pc_z, osf_366, \
                         osf_369, osf_459, qsd0_273, qsd1_273, qsf_456, qsf_457, \
                         qsf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_14 * osf_459[k]
                   + f_3 * pc_x[k] * qsf_459[k];

        t_685[k] = f_13 * osf_366[k]
                   + f_1 * qsd0_273[k]
                   - f_2 * qsd1_273[k]
                   + f_3 * pc_y[k] * qsf_456[k];

        t_686[k] = f_3 * pc_z[k] * qsf_456[k];

        t_687[k] = f_4 * qsd0_273[k]
                   - f_5 * qsd1_273[k]
                   + f_3 * pc_z[k] * qsf_457[k];

        t_688[k] = f_13 * osf_369[k]
                   + f_3 * pc_y[k] * qsf_459[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, pa_z, pc_y, pc_z, osg0_540, osf_360, \
                         osf_370, osg1_540, qsd0_275, qsd1_275, qsf_459, \
                         qsf_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_1 * qsd0_275[k]
                   - f_2 * qsd1_275[k]
                   + f_3 * pc_z[k] * qsf_459[k];

        t_690[k] = pa_z[k] * osg0_540[k]
                   - f_6 * pc_z[k] * osg1_540[k];

        t_691[k] = f_15 * osf_370[k]
                   + f_3 * pc_y[k] * qsf_460[k];

        t_692[k] = f_7 * osf_360[k]
                   + f_3 * pc_z[k] * qsf_460[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pa_z, pc_x, pc_y, pc_z, osg0_543, osf_372, \
                         osf_465, osg1_543, qsd0_281, qsd1_281, qsf_462, \
                         qsf_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = pa_z[k] * osg0_543[k]
                   - f_6 * pc_z[k] * osg1_543[k];

        t_694[k] = f_15 * osf_372[k]
                   + f_3 * pc_y[k] * qsf_462[k];

        t_695[k] = f_14 * osf_465[k]
                   + f_4 * qsd0_281[k]
                   - f_5 * qsd1_281[k]
                   + f_3 * pc_x[k] * qsf_465[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, pc_x, osf_466, osf_467, osf_468, osf_469, \
                         qsf_466, qsf_467, qsf_468, qsf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_14 * osf_466[k]
                   + f_3 * pc_x[k] * qsf_466[k];

        t_697[k] = f_14 * osf_467[k]
                   + f_3 * pc_x[k] * qsf_467[k];

        t_698[k] = f_14 * osf_468[k]
                   + f_3 * pc_x[k] * qsf_468[k];

        t_699[k] = f_14 * osf_469[k]
                   + f_3 * pc_x[k] * qsf_469[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pa_z, pc_y, pc_z, osg0_550, osf_366, osf_378, \
                         osg1_550, qsd0_281, qsd1_281, qsf_466, \
                         qsf_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = pa_z[k] * osg0_550[k]
                   - f_6 * pc_z[k] * osg1_550[k];

        t_701[k] = f_7 * osf_366[k]
                   + f_3 * pc_z[k] * qsf_466[k];

        t_702[k] = f_15 * osf_378[k]
                   + f_4 * qsd0_281[k]
                   - f_5 * qsd1_281[k]
                   + f_3 * pc_y[k] * qsf_468[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, pc_z, osf_369, osf_379, osf_470, \
                         qsd0_281, qsd0_282, qsd1_281, qsd1_282, qsf_469, \
                         qsf_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_15 * osf_379[k]
                   + f_3 * pc_y[k] * qsf_469[k];

        t_704[k] = f_7 * osf_369[k]
                   + f_1 * qsd0_281[k]
                   - f_2 * qsd1_281[k]
                   + f_3 * pc_z[k] * qsf_469[k];

        t_705[k] = f_14 * osf_470[k]
                   + f_1 * qsd0_282[k]
                   - f_2 * qsd1_282[k]
                   + f_3 * pc_x[k] * qsf_470[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, pc_x, pc_y, pc_z, osf_370, osf_380, \
                         osf_382, osf_473, qsd0_285, qsd1_285, qsf_470, qsf_472, \
                         qsf_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_17 * osf_380[k]
                   + f_3 * pc_y[k] * qsf_470[k];

        t_707[k] = f_8 * osf_370[k]
                   + f_3 * pc_z[k] * qsf_470[k];

        t_708[k] = f_14 * osf_473[k]
                   + f_4 * qsd0_285[k]
                   - f_5 * qsd1_285[k]
                   + f_3 * pc_x[k] * qsf_473[k];

        t_709[k] = f_17 * osf_382[k]
                   + f_3 * pc_y[k] * qsf_472[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, osf_475, osf_476, osf_477, osf_478, \
                         qsd0_287, qsd1_287, qsf_475, qsf_476, qsf_477, \
                         qsf_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_14 * osf_475[k]
                   + f_4 * qsd0_287[k]
                   - f_5 * qsd1_287[k]
                   + f_3 * pc_x[k] * qsf_475[k];

        t_711[k] = f_14 * osf_476[k]
                   + f_3 * pc_x[k] * qsf_476[k];

        t_712[k] = f_14 * osf_477[k]
                   + f_3 * pc_x[k] * qsf_477[k];

        t_713[k] = f_14 * osf_478[k]
                   + f_3 * pc_x[k] * qsf_478[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_x, pc_y, pc_z, osf_376, osf_386, osf_479, \
                         qsd0_285, qsd1_285, qsf_476, qsf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_14 * osf_479[k]
                   + f_3 * pc_x[k] * qsf_479[k];

        t_715[k] = f_17 * osf_386[k]
                   + f_1 * qsd0_285[k]
                   - f_2 * qsd1_285[k]
                   + f_3 * pc_y[k] * qsf_476[k];

        t_716[k] = f_8 * osf_376[k]
                   + f_3 * pc_z[k] * qsf_476[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pc_y, pc_z, osf_379, osf_388, osf_389, qsd0_287, \
                         qsd1_287, qsf_478, qsf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_17 * osf_388[k]
                   + f_4 * qsd0_287[k]
                   - f_5 * qsd1_287[k]
                   + f_3 * pc_y[k] * qsf_478[k];

        t_718[k] = f_17 * osf_389[k]
                   + f_3 * pc_y[k] * qsf_479[k];

        t_719[k] = f_8 * osf_379[k]
                   + f_1 * qsd0_287[k]
                   - f_2 * qsd1_287[k]
                   + f_3 * pc_z[k] * qsf_479[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_660 = buffer.data(osg0 + 660);
    const auto *osg0_663 = buffer.data(osg0 + 663);
    const auto *osg0_665 = buffer.data(osg0 + 665);
    const auto *osg0_674 = buffer.data(osg0 + 674);

    const auto *osf_380 = buffer.data(osf + 380);
    const auto *osf_386 = buffer.data(osf + 386);
    const auto *osf_389 = buffer.data(osf + 389);
    const auto *osf_390 = buffer.data(osf + 390);
    const auto *osf_392 = buffer.data(osf + 392);
    const auto *osf_396 = buffer.data(osf + 396);
    const auto *osf_398 = buffer.data(osf + 398);
    const auto *osf_399 = buffer.data(osf + 399);
    const auto *osf_400 = buffer.data(osf + 400);
    const auto *osf_402 = buffer.data(osf + 402);
    const auto *osf_406 = buffer.data(osf + 406);
    const auto *osf_408 = buffer.data(osf + 408);
    const auto *osf_409 = buffer.data(osf + 409);
    const auto *osf_410 = buffer.data(osf + 410);
    const auto *osf_412 = buffer.data(osf + 412);
    const auto *osf_416 = buffer.data(osf + 416);
    const auto *osf_418 = buffer.data(osf + 418);
    const auto *osf_419 = buffer.data(osf + 419);
    const auto *osf_420 = buffer.data(osf + 420);
    const auto *osf_422 = buffer.data(osf + 422);
    const auto *osf_426 = buffer.data(osf + 426);
    const auto *osf_428 = buffer.data(osf + 428);
    const auto *osf_429 = buffer.data(osf + 429);
    const auto *osf_430 = buffer.data(osf + 430);
    const auto *osf_432 = buffer.data(osf + 432);
    const auto *osf_436 = buffer.data(osf + 436);
    const auto *osf_438 = buffer.data(osf + 438);
    const auto *osf_439 = buffer.data(osf + 439);
    const auto *osf_440 = buffer.data(osf + 440);
    const auto *osf_441 = buffer.data(osf + 441);
    const auto *osf_442 = buffer.data(osf + 442);
    const auto *osf_446 = buffer.data(osf + 446);
    const auto *osf_448 = buffer.data(osf + 448);
    const auto *osf_449 = buffer.data(osf + 449);
    const auto *osf_450 = buffer.data(osf + 450);
    const auto *osf_480 = buffer.data(osf + 480);
    const auto *osf_483 = buffer.data(osf + 483);
    const auto *osf_485 = buffer.data(osf + 485);
    const auto *osf_486 = buffer.data(osf + 486);
    const auto *osf_487 = buffer.data(osf + 487);
    const auto *osf_488 = buffer.data(osf + 488);
    const auto *osf_489 = buffer.data(osf + 489);
    const auto *osf_490 = buffer.data(osf + 490);
    const auto *osf_493 = buffer.data(osf + 493);
    const auto *osf_495 = buffer.data(osf + 495);
    const auto *osf_496 = buffer.data(osf + 496);
    const auto *osf_497 = buffer.data(osf + 497);
    const auto *osf_498 = buffer.data(osf + 498);
    const auto *osf_499 = buffer.data(osf + 499);
    const auto *osf_500 = buffer.data(osf + 500);
    const auto *osf_503 = buffer.data(osf + 503);
    const auto *osf_505 = buffer.data(osf + 505);
    const auto *osf_506 = buffer.data(osf + 506);
    const auto *osf_507 = buffer.data(osf + 507);
    const auto *osf_508 = buffer.data(osf + 508);
    const auto *osf_509 = buffer.data(osf + 509);
    const auto *osf_510 = buffer.data(osf + 510);
    const auto *osf_513 = buffer.data(osf + 513);
    const auto *osf_515 = buffer.data(osf + 515);
    const auto *osf_516 = buffer.data(osf + 516);
    const auto *osf_517 = buffer.data(osf + 517);
    const auto *osf_518 = buffer.data(osf + 518);
    const auto *osf_519 = buffer.data(osf + 519);
    const auto *osf_520 = buffer.data(osf + 520);
    const auto *osf_523 = buffer.data(osf + 523);
    const auto *osf_525 = buffer.data(osf + 525);
    const auto *osf_526 = buffer.data(osf + 526);
    const auto *osf_527 = buffer.data(osf + 527);
    const auto *osf_528 = buffer.data(osf + 528);
    const auto *osf_529 = buffer.data(osf + 529);
    const auto *osf_536 = buffer.data(osf + 536);
    const auto *osf_537 = buffer.data(osf + 537);
    const auto *osf_538 = buffer.data(osf + 538);
    const auto *osf_539 = buffer.data(osf + 539);
    const auto *osf_540 = buffer.data(osf + 540);
    const auto *osf_545 = buffer.data(osf + 545);
    const auto *osf_546 = buffer.data(osf + 546);
    const auto *osf_547 = buffer.data(osf + 547);
    const auto *osf_549 = buffer.data(osf + 549);
    const auto *osf_550 = buffer.data(osf + 550);
    const auto *osf_553 = buffer.data(osf + 553);
    const auto *osf_556 = buffer.data(osf + 556);
    const auto *osf_558 = buffer.data(osf + 558);

    const auto *osg1_660 = buffer.data(osg1 + 660);
    const auto *osg1_663 = buffer.data(osg1 + 663);
    const auto *osg1_665 = buffer.data(osg1 + 665);
    const auto *osg1_674 = buffer.data(osg1 + 674);

    const auto *qsd0_288 = buffer.data(qsd0 + 288);
    const auto *qsd0_291 = buffer.data(qsd0 + 291);
    const auto *qsd0_293 = buffer.data(qsd0 + 293);
    const auto *qsd0_294 = buffer.data(qsd0 + 294);
    const auto *qsd0_297 = buffer.data(qsd0 + 297);
    const auto *qsd0_299 = buffer.data(qsd0 + 299);
    const auto *qsd0_300 = buffer.data(qsd0 + 300);
    const auto *qsd0_303 = buffer.data(qsd0 + 303);
    const auto *qsd0_305 = buffer.data(qsd0 + 305);
    const auto *qsd0_306 = buffer.data(qsd0 + 306);
    const auto *qsd0_309 = buffer.data(qsd0 + 309);
    const auto *qsd0_311 = buffer.data(qsd0 + 311);
    const auto *qsd0_312 = buffer.data(qsd0 + 312);
    const auto *qsd0_315 = buffer.data(qsd0 + 315);
    const auto *qsd0_317 = buffer.data(qsd0 + 317);
    const auto *qsd0_321 = buffer.data(qsd0 + 321);
    const auto *qsd0_323 = buffer.data(qsd0 + 323);
    const auto *qsd0_324 = buffer.data(qsd0 + 324);
    const auto *qsd0_327 = buffer.data(qsd0 + 327);
    const auto *qsd0_328 = buffer.data(qsd0 + 328);
    const auto *qsd0_329 = buffer.data(qsd0 + 329);
    const auto *qsd0_330 = buffer.data(qsd0 + 330);
    const auto *qsd0_333 = buffer.data(qsd0 + 333);

    const auto *qsd1_288 = buffer.data(qsd1 + 288);
    const auto *qsd1_291 = buffer.data(qsd1 + 291);
    const auto *qsd1_293 = buffer.data(qsd1 + 293);
    const auto *qsd1_294 = buffer.data(qsd1 + 294);
    const auto *qsd1_297 = buffer.data(qsd1 + 297);
    const auto *qsd1_299 = buffer.data(qsd1 + 299);
    const auto *qsd1_300 = buffer.data(qsd1 + 300);
    const auto *qsd1_303 = buffer.data(qsd1 + 303);
    const auto *qsd1_305 = buffer.data(qsd1 + 305);
    const auto *qsd1_306 = buffer.data(qsd1 + 306);
    const auto *qsd1_309 = buffer.data(qsd1 + 309);
    const auto *qsd1_311 = buffer.data(qsd1 + 311);
    const auto *qsd1_312 = buffer.data(qsd1 + 312);
    const auto *qsd1_315 = buffer.data(qsd1 + 315);
    const auto *qsd1_317 = buffer.data(qsd1 + 317);
    const auto *qsd1_321 = buffer.data(qsd1 + 321);
    const auto *qsd1_323 = buffer.data(qsd1 + 323);
    const auto *qsd1_324 = buffer.data(qsd1 + 324);
    const auto *qsd1_327 = buffer.data(qsd1 + 327);
    const auto *qsd1_328 = buffer.data(qsd1 + 328);
    const auto *qsd1_329 = buffer.data(qsd1 + 329);
    const auto *qsd1_330 = buffer.data(qsd1 + 330);
    const auto *qsd1_333 = buffer.data(qsd1 + 333);

    const auto *qsf_480 = buffer.data(qsf + 480);
    const auto *qsf_482 = buffer.data(qsf + 482);
    const auto *qsf_483 = buffer.data(qsf + 483);
    const auto *qsf_485 = buffer.data(qsf + 485);
    const auto *qsf_486 = buffer.data(qsf + 486);
    const auto *qsf_487 = buffer.data(qsf + 487);
    const auto *qsf_488 = buffer.data(qsf + 488);
    const auto *qsf_489 = buffer.data(qsf + 489);
    const auto *qsf_490 = buffer.data(qsf + 490);
    const auto *qsf_492 = buffer.data(qsf + 492);
    const auto *qsf_493 = buffer.data(qsf + 493);
    const auto *qsf_495 = buffer.data(qsf + 495);
    const auto *qsf_496 = buffer.data(qsf + 496);
    const auto *qsf_497 = buffer.data(qsf + 497);
    const auto *qsf_498 = buffer.data(qsf + 498);
    const auto *qsf_499 = buffer.data(qsf + 499);
    const auto *qsf_500 = buffer.data(qsf + 500);
    const auto *qsf_502 = buffer.data(qsf + 502);
    const auto *qsf_503 = buffer.data(qsf + 503);
    const auto *qsf_505 = buffer.data(qsf + 505);
    const auto *qsf_506 = buffer.data(qsf + 506);
    const auto *qsf_507 = buffer.data(qsf + 507);
    const auto *qsf_508 = buffer.data(qsf + 508);
    const auto *qsf_509 = buffer.data(qsf + 509);
    const auto *qsf_510 = buffer.data(qsf + 510);
    const auto *qsf_512 = buffer.data(qsf + 512);
    const auto *qsf_513 = buffer.data(qsf + 513);
    const auto *qsf_515 = buffer.data(qsf + 515);
    const auto *qsf_516 = buffer.data(qsf + 516);
    const auto *qsf_517 = buffer.data(qsf + 517);
    const auto *qsf_518 = buffer.data(qsf + 518);
    const auto *qsf_519 = buffer.data(qsf + 519);
    const auto *qsf_520 = buffer.data(qsf + 520);
    const auto *qsf_522 = buffer.data(qsf + 522);
    const auto *qsf_523 = buffer.data(qsf + 523);
    const auto *qsf_525 = buffer.data(qsf + 525);
    const auto *qsf_526 = buffer.data(qsf + 526);
    const auto *qsf_527 = buffer.data(qsf + 527);
    const auto *qsf_528 = buffer.data(qsf + 528);
    const auto *qsf_529 = buffer.data(qsf + 529);
    const auto *qsf_530 = buffer.data(qsf + 530);
    const auto *qsf_532 = buffer.data(qsf + 532);
    const auto *qsf_536 = buffer.data(qsf + 536);
    const auto *qsf_537 = buffer.data(qsf + 537);
    const auto *qsf_538 = buffer.data(qsf + 538);
    const auto *qsf_539 = buffer.data(qsf + 539);
    const auto *qsf_540 = buffer.data(qsf + 540);
    const auto *qsf_541 = buffer.data(qsf + 541);
    const auto *qsf_542 = buffer.data(qsf + 542);
    const auto *qsf_545 = buffer.data(qsf + 545);
    const auto *qsf_546 = buffer.data(qsf + 546);
    const auto *qsf_547 = buffer.data(qsf + 547);
    const auto *qsf_548 = buffer.data(qsf + 548);
    const auto *qsf_549 = buffer.data(qsf + 549);
    const auto *qsf_550 = buffer.data(qsf + 550);
    const auto *qsf_551 = buffer.data(qsf + 551);
    const auto *qsf_552 = buffer.data(qsf + 552);
    const auto *qsf_553 = buffer.data(qsf + 553);
    const auto *qsf_556 = buffer.data(qsf + 556);
    const auto *qsf_558 = buffer.data(qsf + 558);

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, osf_380, osf_390, osf_480, \
                         qsd0_288, qsd1_288, qsf_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_14 * osf_480[k]
                   + f_1 * qsd0_288[k]
                   - f_2 * qsd1_288[k]
                   + f_3 * pc_x[k] * qsf_480[k];

        t_721[k] = f_19 * osf_390[k]
                   + f_3 * pc_y[k] * qsf_480[k];

        t_722[k] = f_14 * osf_380[k]
                   + f_3 * pc_z[k] * qsf_480[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_x, pc_y, osf_392, osf_483, osf_485, qsd0_291, \
                         qsd0_293, qsd1_291, qsd1_293, qsf_482, qsf_483, \
                         qsf_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_14 * osf_483[k]
                   + f_4 * qsd0_291[k]
                   - f_5 * qsd1_291[k]
                   + f_3 * pc_x[k] * qsf_483[k];

        t_724[k] = f_19 * osf_392[k]
                   + f_3 * pc_y[k] * qsf_482[k];

        t_725[k] = f_14 * osf_485[k]
                   + f_4 * qsd0_293[k]
                   - f_5 * qsd1_293[k]
                   + f_3 * pc_x[k] * qsf_485[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pc_x, osf_486, osf_487, osf_488, osf_489, \
                         qsf_486, qsf_487, qsf_488, qsf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_14 * osf_486[k]
                   + f_3 * pc_x[k] * qsf_486[k];

        t_727[k] = f_14 * osf_487[k]
                   + f_3 * pc_x[k] * qsf_487[k];

        t_728[k] = f_14 * osf_488[k]
                   + f_3 * pc_x[k] * qsf_488[k];

        t_729[k] = f_14 * osf_489[k]
                   + f_3 * pc_x[k] * qsf_489[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_y, pc_z, osf_386, osf_396, osf_398, qsd0_291, \
                         qsd0_293, qsd1_291, qsd1_293, qsf_486, \
                         qsf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_19 * osf_396[k]
                   + f_1 * qsd0_291[k]
                   - f_2 * qsd1_291[k]
                   + f_3 * pc_y[k] * qsf_486[k];

        t_731[k] = f_14 * osf_386[k]
                   + f_3 * pc_z[k] * qsf_486[k];

        t_732[k] = f_19 * osf_398[k]
                   + f_4 * qsd0_293[k]
                   - f_5 * qsd1_293[k]
                   + f_3 * pc_y[k] * qsf_488[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, pc_z, osf_389, osf_399, osf_490, \
                         qsd0_293, qsd0_294, qsd1_293, qsd1_294, qsf_489, \
                         qsf_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_19 * osf_399[k]
                   + f_3 * pc_y[k] * qsf_489[k];

        t_734[k] = f_14 * osf_389[k]
                   + f_1 * qsd0_293[k]
                   - f_2 * qsd1_293[k]
                   + f_3 * pc_z[k] * qsf_489[k];

        t_735[k] = f_14 * osf_490[k]
                   + f_1 * qsd0_294[k]
                   - f_2 * qsd1_294[k]
                   + f_3 * pc_x[k] * qsf_490[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, osf_390, osf_400, \
                         osf_402, osf_493, qsd0_297, qsd1_297, qsf_490, qsf_492, \
                         qsf_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_18 * osf_400[k]
                   + f_3 * pc_y[k] * qsf_490[k];

        t_737[k] = f_16 * osf_390[k]
                   + f_3 * pc_z[k] * qsf_490[k];

        t_738[k] = f_14 * osf_493[k]
                   + f_4 * qsd0_297[k]
                   - f_5 * qsd1_297[k]
                   + f_3 * pc_x[k] * qsf_493[k];

        t_739[k] = f_18 * osf_402[k]
                   + f_3 * pc_y[k] * qsf_492[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, osf_495, osf_496, osf_497, osf_498, \
                         qsd0_299, qsd1_299, qsf_495, qsf_496, qsf_497, \
                         qsf_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_14 * osf_495[k]
                   + f_4 * qsd0_299[k]
                   - f_5 * qsd1_299[k]
                   + f_3 * pc_x[k] * qsf_495[k];

        t_741[k] = f_14 * osf_496[k]
                   + f_3 * pc_x[k] * qsf_496[k];

        t_742[k] = f_14 * osf_497[k]
                   + f_3 * pc_x[k] * qsf_497[k];

        t_743[k] = f_14 * osf_498[k]
                   + f_3 * pc_x[k] * qsf_498[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, pc_x, pc_y, pc_z, osf_396, osf_406, osf_499, \
                         qsd0_297, qsd1_297, qsf_496, qsf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_14 * osf_499[k]
                   + f_3 * pc_x[k] * qsf_499[k];

        t_745[k] = f_18 * osf_406[k]
                   + f_1 * qsd0_297[k]
                   - f_2 * qsd1_297[k]
                   + f_3 * pc_y[k] * qsf_496[k];

        t_746[k] = f_16 * osf_396[k]
                   + f_3 * pc_z[k] * qsf_496[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_y, pc_z, osf_399, osf_408, osf_409, qsd0_299, \
                         qsd1_299, qsf_498, qsf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_18 * osf_408[k]
                   + f_4 * qsd0_299[k]
                   - f_5 * qsd1_299[k]
                   + f_3 * pc_y[k] * qsf_498[k];

        t_748[k] = f_18 * osf_409[k]
                   + f_3 * pc_y[k] * qsf_499[k];

        t_749[k] = f_16 * osf_399[k]
                   + f_1 * qsd0_299[k]
                   - f_2 * qsd1_299[k]
                   + f_3 * pc_z[k] * qsf_499[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_x, pc_y, pc_z, osf_400, osf_410, osf_500, \
                         qsd0_300, qsd1_300, qsf_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_14 * osf_500[k]
                   + f_1 * qsd0_300[k]
                   - f_2 * qsd1_300[k]
                   + f_3 * pc_x[k] * qsf_500[k];

        t_751[k] = f_16 * osf_410[k]
                   + f_3 * pc_y[k] * qsf_500[k];

        t_752[k] = f_18 * osf_400[k]
                   + f_3 * pc_z[k] * qsf_500[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pc_x, pc_y, osf_412, osf_503, osf_505, qsd0_303, \
                         qsd0_305, qsd1_303, qsd1_305, qsf_502, qsf_503, \
                         qsf_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_14 * osf_503[k]
                   + f_4 * qsd0_303[k]
                   - f_5 * qsd1_303[k]
                   + f_3 * pc_x[k] * qsf_503[k];

        t_754[k] = f_16 * osf_412[k]
                   + f_3 * pc_y[k] * qsf_502[k];

        t_755[k] = f_14 * osf_505[k]
                   + f_4 * qsd0_305[k]
                   - f_5 * qsd1_305[k]
                   + f_3 * pc_x[k] * qsf_505[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, osf_506, osf_507, osf_508, osf_509, \
                         qsf_506, qsf_507, qsf_508, qsf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_14 * osf_506[k]
                   + f_3 * pc_x[k] * qsf_506[k];

        t_757[k] = f_14 * osf_507[k]
                   + f_3 * pc_x[k] * qsf_507[k];

        t_758[k] = f_14 * osf_508[k]
                   + f_3 * pc_x[k] * qsf_508[k];

        t_759[k] = f_14 * osf_509[k]
                   + f_3 * pc_x[k] * qsf_509[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pc_y, pc_z, osf_406, osf_416, osf_418, qsd0_303, \
                         qsd0_305, qsd1_303, qsd1_305, qsf_506, \
                         qsf_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_16 * osf_416[k]
                   + f_1 * qsd0_303[k]
                   - f_2 * qsd1_303[k]
                   + f_3 * pc_y[k] * qsf_506[k];

        t_761[k] = f_18 * osf_406[k]
                   + f_3 * pc_z[k] * qsf_506[k];

        t_762[k] = f_16 * osf_418[k]
                   + f_4 * qsd0_305[k]
                   - f_5 * qsd1_305[k]
                   + f_3 * pc_y[k] * qsf_508[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pc_x, pc_y, pc_z, osf_409, osf_419, osf_510, \
                         qsd0_305, qsd0_306, qsd1_305, qsd1_306, qsf_509, \
                         qsf_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_16 * osf_419[k]
                   + f_3 * pc_y[k] * qsf_509[k];

        t_764[k] = f_18 * osf_409[k]
                   + f_1 * qsd0_305[k]
                   - f_2 * qsd1_305[k]
                   + f_3 * pc_z[k] * qsf_509[k];

        t_765[k] = f_14 * osf_510[k]
                   + f_1 * qsd0_306[k]
                   - f_2 * qsd1_306[k]
                   + f_3 * pc_x[k] * qsf_510[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, pc_x, pc_y, pc_z, osf_410, osf_420, \
                         osf_422, osf_513, qsd0_309, qsd1_309, qsf_510, qsf_512, \
                         qsf_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_14 * osf_420[k]
                   + f_3 * pc_y[k] * qsf_510[k];

        t_767[k] = f_19 * osf_410[k]
                   + f_3 * pc_z[k] * qsf_510[k];

        t_768[k] = f_14 * osf_513[k]
                   + f_4 * qsd0_309[k]
                   - f_5 * qsd1_309[k]
                   + f_3 * pc_x[k] * qsf_513[k];

        t_769[k] = f_14 * osf_422[k]
                   + f_3 * pc_y[k] * qsf_512[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pc_x, osf_515, osf_516, osf_517, osf_518, \
                         qsd0_311, qsd1_311, qsf_515, qsf_516, qsf_517, \
                         qsf_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_14 * osf_515[k]
                   + f_4 * qsd0_311[k]
                   - f_5 * qsd1_311[k]
                   + f_3 * pc_x[k] * qsf_515[k];

        t_771[k] = f_14 * osf_516[k]
                   + f_3 * pc_x[k] * qsf_516[k];

        t_772[k] = f_14 * osf_517[k]
                   + f_3 * pc_x[k] * qsf_517[k];

        t_773[k] = f_14 * osf_518[k]
                   + f_3 * pc_x[k] * qsf_518[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, pc_x, pc_y, pc_z, osf_416, osf_426, osf_519, \
                         qsd0_309, qsd1_309, qsf_516, qsf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_14 * osf_519[k]
                   + f_3 * pc_x[k] * qsf_519[k];

        t_775[k] = f_14 * osf_426[k]
                   + f_1 * qsd0_309[k]
                   - f_2 * qsd1_309[k]
                   + f_3 * pc_y[k] * qsf_516[k];

        t_776[k] = f_19 * osf_416[k]
                   + f_3 * pc_z[k] * qsf_516[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, pc_z, osf_419, osf_428, osf_429, qsd0_311, \
                         qsd1_311, qsf_518, qsf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_14 * osf_428[k]
                   + f_4 * qsd0_311[k]
                   - f_5 * qsd1_311[k]
                   + f_3 * pc_y[k] * qsf_518[k];

        t_778[k] = f_14 * osf_429[k]
                   + f_3 * pc_y[k] * qsf_519[k];

        t_779[k] = f_19 * osf_419[k]
                   + f_1 * qsd0_311[k]
                   - f_2 * qsd1_311[k]
                   + f_3 * pc_z[k] * qsf_519[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pc_x, pc_y, pc_z, osf_420, osf_430, osf_520, \
                         qsd0_312, qsd1_312, qsf_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_14 * osf_520[k]
                   + f_1 * qsd0_312[k]
                   - f_2 * qsd1_312[k]
                   + f_3 * pc_x[k] * qsf_520[k];

        t_781[k] = f_8 * osf_430[k]
                   + f_3 * pc_y[k] * qsf_520[k];

        t_782[k] = f_17 * osf_420[k]
                   + f_3 * pc_z[k] * qsf_520[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pc_x, pc_y, osf_432, osf_523, osf_525, qsd0_315, \
                         qsd0_317, qsd1_315, qsd1_317, qsf_522, qsf_523, \
                         qsf_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_14 * osf_523[k]
                   + f_4 * qsd0_315[k]
                   - f_5 * qsd1_315[k]
                   + f_3 * pc_x[k] * qsf_523[k];

        t_784[k] = f_8 * osf_432[k]
                   + f_3 * pc_y[k] * qsf_522[k];

        t_785[k] = f_14 * osf_525[k]
                   + f_4 * qsd0_317[k]
                   - f_5 * qsd1_317[k]
                   + f_3 * pc_x[k] * qsf_525[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pc_x, osf_526, osf_527, osf_528, osf_529, \
                         qsf_526, qsf_527, qsf_528, qsf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_14 * osf_526[k]
                   + f_3 * pc_x[k] * qsf_526[k];

        t_787[k] = f_14 * osf_527[k]
                   + f_3 * pc_x[k] * qsf_527[k];

        t_788[k] = f_14 * osf_528[k]
                   + f_3 * pc_x[k] * qsf_528[k];

        t_789[k] = f_14 * osf_529[k]
                   + f_3 * pc_x[k] * qsf_529[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, pc_y, pc_z, osf_426, osf_436, osf_438, qsd0_315, \
                         qsd0_317, qsd1_315, qsd1_317, qsf_526, \
                         qsf_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_8 * osf_436[k]
                   + f_1 * qsd0_315[k]
                   - f_2 * qsd1_315[k]
                   + f_3 * pc_y[k] * qsf_526[k];

        t_791[k] = f_17 * osf_426[k]
                   + f_3 * pc_z[k] * qsf_526[k];

        t_792[k] = f_8 * osf_438[k]
                   + f_4 * qsd0_317[k]
                   - f_5 * qsd1_317[k]
                   + f_3 * pc_y[k] * qsf_528[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, pa_y, pc_y, pc_z, osg0_660, osf_429, \
                         osf_439, osf_440, osg1_660, qsd0_317, qsd1_317, qsf_529, \
                         qsf_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_8 * osf_439[k]
                   + f_3 * pc_y[k] * qsf_529[k];

        t_794[k] = f_17 * osf_429[k]
                   + f_1 * qsd0_317[k]
                   - f_2 * qsd1_317[k]
                   + f_3 * pc_z[k] * qsf_529[k];

        t_795[k] = pa_y[k] * osg0_660[k]
                   - f_6 * pc_y[k] * osg1_660[k];

        t_796[k] = f_7 * osf_440[k]
                   + f_3 * pc_y[k] * qsf_530[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pa_y, pc_y, pc_z, osg0_663, osg0_665, \
                         osf_430, osf_441, osf_442, osg1_663, osg1_665, qsf_530, \
                         qsf_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_15 * osf_430[k]
                   + f_3 * pc_z[k] * qsf_530[k];

        t_798[k] = pa_y[k] * osg0_663[k]
                   + f_8 * osf_441[k]
                   - f_6 * pc_y[k] * osg1_663[k];

        t_799[k] = f_7 * osf_442[k]
                   + f_3 * pc_y[k] * qsf_532[k];

        t_800[k] = pa_y[k] * osg0_665[k]
                   - f_6 * pc_y[k] * osg1_665[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pc_x, osf_536, osf_537, osf_538, osf_539, \
                         qsf_536, qsf_537, qsf_538, qsf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_14 * osf_536[k]
                   + f_3 * pc_x[k] * qsf_536[k];

        t_802[k] = f_14 * osf_537[k]
                   + f_3 * pc_x[k] * qsf_537[k];

        t_803[k] = f_14 * osf_538[k]
                   + f_3 * pc_x[k] * qsf_538[k];

        t_804[k] = f_14 * osf_539[k]
                   + f_3 * pc_x[k] * qsf_539[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pc_y, pc_z, osf_436, osf_446, osf_448, qsd0_321, \
                         qsd0_323, qsd1_321, qsd1_323, qsf_536, \
                         qsf_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_7 * osf_446[k]
                   + f_1 * qsd0_321[k]
                   - f_2 * qsd1_321[k]
                   + f_3 * pc_y[k] * qsf_536[k];

        t_806[k] = f_15 * osf_436[k]
                   + f_3 * pc_z[k] * qsf_536[k];

        t_807[k] = f_7 * osf_448[k]
                   + f_4 * qsd0_323[k]
                   - f_5 * qsd1_323[k]
                   + f_3 * pc_y[k] * qsf_538[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pa_y, pc_x, pc_y, osg0_674, osf_449, \
                         osf_540, osg1_674, qsd0_324, qsd1_324, qsf_539, \
                         qsf_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_7 * osf_449[k]
                   + f_3 * pc_y[k] * qsf_539[k];

        t_809[k] = pa_y[k] * osg0_674[k]
                   - f_6 * pc_y[k] * osg1_674[k];

        t_810[k] = f_14 * osf_540[k]
                   + f_1 * qsd0_324[k]
                   - f_2 * qsd1_324[k]
                   + f_3 * pc_x[k] * qsf_540[k];

        t_811[k] = f_3 * pc_y[k] * qsf_540[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pc_y, pc_z, osf_440, qsd0_324, qsd1_324, \
                         qsf_540, qsf_541, qsf_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_13 * osf_440[k]
                   + f_3 * pc_z[k] * qsf_540[k];

        t_813[k] = f_4 * qsd0_324[k]
                   - f_5 * qsd1_324[k]
                   + f_3 * pc_y[k] * qsf_541[k];

        t_814[k] = f_3 * pc_y[k] * qsf_542[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pc_x, pc_y, osf_545, osf_546, osf_547, \
                         qsd0_329, qsd1_329, qsf_545, qsf_546, \
                         qsf_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_14 * osf_545[k]
                   + f_4 * qsd0_329[k]
                   - f_5 * qsd1_329[k]
                   + f_3 * pc_x[k] * qsf_545[k];

        t_816[k] = f_14 * osf_546[k]
                   + f_3 * pc_x[k] * qsf_546[k];

        t_817[k] = f_14 * osf_547[k]
                   + f_3 * pc_x[k] * qsf_547[k];

        t_818[k] = f_3 * pc_y[k] * qsf_545[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_y, osf_549, qsd0_327, qsd0_328, \
                         qsd1_327, qsd1_328, qsf_546, qsf_547, \
                         qsf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_14 * osf_549[k]
                   + f_3 * pc_x[k] * qsf_549[k];

        t_820[k] = f_1 * qsd0_327[k]
                   - f_2 * qsd1_327[k]
                   + f_3 * pc_y[k] * qsf_546[k];

        t_821[k] = f_10 * qsd0_328[k]
                   - f_11 * qsd1_328[k]
                   + f_3 * pc_y[k] * qsf_547[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pc_x, pc_y, pc_z, osf_449, osf_550, \
                         qsd0_329, qsd0_330, qsd1_329, qsd1_330, qsf_548, qsf_549, \
                         qsf_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_4 * qsd0_329[k]
                   - f_5 * qsd1_329[k]
                   + f_3 * pc_y[k] * qsf_548[k];

        t_823[k] = f_3 * pc_y[k] * qsf_549[k];

        t_824[k] = f_13 * osf_449[k]
                   + f_1 * qsd0_329[k]
                   - f_2 * qsd1_329[k]
                   + f_3 * pc_z[k] * qsf_549[k];

        t_825[k] = f_8 * osf_550[k]
                   + f_1 * qsd0_330[k]
                   - f_2 * qsd1_330[k]
                   + f_3 * pc_x[k] * qsf_550[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, pc_y, pc_z, osf_450, osf_553, \
                         qsd0_333, qsd1_333, qsf_550, qsf_551, \
                         qsf_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_12 * osf_450[k]
                   + f_3 * pc_y[k] * qsf_550[k];

        t_827[k] = f_3 * pc_z[k] * qsf_550[k];

        t_828[k] = f_8 * osf_553[k]
                   + f_4 * qsd0_333[k]
                   - f_5 * qsd1_333[k]
                   + f_3 * pc_x[k] * qsf_553[k];

        t_829[k] = f_3 * pc_z[k] * qsf_551[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pc_x, pc_z, osf_556, osf_558, qsd0_330, \
                         qsd1_330, qsf_552, qsf_553, qsf_556, qsf_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_4 * qsd0_330[k]
                   - f_5 * qsd1_330[k]
                   + f_3 * pc_z[k] * qsf_552[k];

        t_831[k] = f_8 * osf_556[k]
                   + f_3 * pc_x[k] * qsf_556[k];

        t_832[k] = f_3 * pc_z[k] * qsf_553[k];

        t_833[k] = f_8 * osf_558[k]
                   + f_3 * pc_x[k] * qsf_558[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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
    auto *t_945 = buffer.data(target + 945);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_675 = buffer.data(osg0 + 675);
    const auto *osg0_678 = buffer.data(osg0 + 678);
    const auto *osg0_685 = buffer.data(osg0 + 685);

    const auto *osf_450 = buffer.data(osf + 450);
    const auto *osf_456 = buffer.data(osf + 456);
    const auto *osf_459 = buffer.data(osf + 459);
    const auto *osf_460 = buffer.data(osf + 460);
    const auto *osf_462 = buffer.data(osf + 462);
    const auto *osf_466 = buffer.data(osf + 466);
    const auto *osf_468 = buffer.data(osf + 468);
    const auto *osf_469 = buffer.data(osf + 469);
    const auto *osf_470 = buffer.data(osf + 470);
    const auto *osf_472 = buffer.data(osf + 472);
    const auto *osf_476 = buffer.data(osf + 476);
    const auto *osf_478 = buffer.data(osf + 478);
    const auto *osf_479 = buffer.data(osf + 479);
    const auto *osf_480 = buffer.data(osf + 480);
    const auto *osf_482 = buffer.data(osf + 482);
    const auto *osf_486 = buffer.data(osf + 486);
    const auto *osf_488 = buffer.data(osf + 488);
    const auto *osf_489 = buffer.data(osf + 489);
    const auto *osf_490 = buffer.data(osf + 490);
    const auto *osf_492 = buffer.data(osf + 492);
    const auto *osf_496 = buffer.data(osf + 496);
    const auto *osf_498 = buffer.data(osf + 498);
    const auto *osf_499 = buffer.data(osf + 499);
    const auto *osf_500 = buffer.data(osf + 500);
    const auto *osf_502 = buffer.data(osf + 502);
    const auto *osf_506 = buffer.data(osf + 506);
    const auto *osf_508 = buffer.data(osf + 508);
    const auto *osf_509 = buffer.data(osf + 509);
    const auto *osf_510 = buffer.data(osf + 510);
    const auto *osf_512 = buffer.data(osf + 512);
    const auto *osf_516 = buffer.data(osf + 516);
    const auto *osf_518 = buffer.data(osf + 518);
    const auto *osf_519 = buffer.data(osf + 519);
    const auto *osf_520 = buffer.data(osf + 520);
    const auto *osf_522 = buffer.data(osf + 522);
    const auto *osf_526 = buffer.data(osf + 526);
    const auto *osf_528 = buffer.data(osf + 528);
    const auto *osf_529 = buffer.data(osf + 529);
    const auto *osf_559 = buffer.data(osf + 559);
    const auto *osf_565 = buffer.data(osf + 565);
    const auto *osf_566 = buffer.data(osf + 566);
    const auto *osf_567 = buffer.data(osf + 567);
    const auto *osf_568 = buffer.data(osf + 568);
    const auto *osf_569 = buffer.data(osf + 569);
    const auto *osf_570 = buffer.data(osf + 570);
    const auto *osf_573 = buffer.data(osf + 573);
    const auto *osf_575 = buffer.data(osf + 575);
    const auto *osf_576 = buffer.data(osf + 576);
    const auto *osf_577 = buffer.data(osf + 577);
    const auto *osf_578 = buffer.data(osf + 578);
    const auto *osf_579 = buffer.data(osf + 579);
    const auto *osf_580 = buffer.data(osf + 580);
    const auto *osf_583 = buffer.data(osf + 583);
    const auto *osf_585 = buffer.data(osf + 585);
    const auto *osf_586 = buffer.data(osf + 586);
    const auto *osf_587 = buffer.data(osf + 587);
    const auto *osf_588 = buffer.data(osf + 588);
    const auto *osf_589 = buffer.data(osf + 589);
    const auto *osf_590 = buffer.data(osf + 590);
    const auto *osf_593 = buffer.data(osf + 593);
    const auto *osf_595 = buffer.data(osf + 595);
    const auto *osf_596 = buffer.data(osf + 596);
    const auto *osf_597 = buffer.data(osf + 597);
    const auto *osf_598 = buffer.data(osf + 598);
    const auto *osf_599 = buffer.data(osf + 599);
    const auto *osf_600 = buffer.data(osf + 600);
    const auto *osf_603 = buffer.data(osf + 603);
    const auto *osf_605 = buffer.data(osf + 605);
    const auto *osf_606 = buffer.data(osf + 606);
    const auto *osf_607 = buffer.data(osf + 607);
    const auto *osf_608 = buffer.data(osf + 608);
    const auto *osf_609 = buffer.data(osf + 609);
    const auto *osf_610 = buffer.data(osf + 610);
    const auto *osf_613 = buffer.data(osf + 613);
    const auto *osf_615 = buffer.data(osf + 615);
    const auto *osf_616 = buffer.data(osf + 616);
    const auto *osf_617 = buffer.data(osf + 617);
    const auto *osf_618 = buffer.data(osf + 618);
    const auto *osf_619 = buffer.data(osf + 619);
    const auto *osf_620 = buffer.data(osf + 620);
    const auto *osf_623 = buffer.data(osf + 623);
    const auto *osf_625 = buffer.data(osf + 625);
    const auto *osf_626 = buffer.data(osf + 626);
    const auto *osf_627 = buffer.data(osf + 627);
    const auto *osf_628 = buffer.data(osf + 628);
    const auto *osf_629 = buffer.data(osf + 629);
    const auto *osf_630 = buffer.data(osf + 630);

    const auto *osg1_675 = buffer.data(osg1 + 675);
    const auto *osg1_678 = buffer.data(osg1 + 678);
    const auto *osg1_685 = buffer.data(osg1 + 685);

    const auto *qsd0_333 = buffer.data(qsd0 + 333);
    const auto *qsd0_335 = buffer.data(qsd0 + 335);
    const auto *qsd0_341 = buffer.data(qsd0 + 341);
    const auto *qsd0_342 = buffer.data(qsd0 + 342);
    const auto *qsd0_345 = buffer.data(qsd0 + 345);
    const auto *qsd0_347 = buffer.data(qsd0 + 347);
    const auto *qsd0_348 = buffer.data(qsd0 + 348);
    const auto *qsd0_351 = buffer.data(qsd0 + 351);
    const auto *qsd0_353 = buffer.data(qsd0 + 353);
    const auto *qsd0_354 = buffer.data(qsd0 + 354);
    const auto *qsd0_357 = buffer.data(qsd0 + 357);
    const auto *qsd0_359 = buffer.data(qsd0 + 359);
    const auto *qsd0_360 = buffer.data(qsd0 + 360);
    const auto *qsd0_363 = buffer.data(qsd0 + 363);
    const auto *qsd0_365 = buffer.data(qsd0 + 365);
    const auto *qsd0_366 = buffer.data(qsd0 + 366);
    const auto *qsd0_369 = buffer.data(qsd0 + 369);
    const auto *qsd0_371 = buffer.data(qsd0 + 371);
    const auto *qsd0_372 = buffer.data(qsd0 + 372);
    const auto *qsd0_375 = buffer.data(qsd0 + 375);
    const auto *qsd0_377 = buffer.data(qsd0 + 377);
    const auto *qsd0_378 = buffer.data(qsd0 + 378);

    const auto *qsd1_333 = buffer.data(qsd1 + 333);
    const auto *qsd1_335 = buffer.data(qsd1 + 335);
    const auto *qsd1_341 = buffer.data(qsd1 + 341);
    const auto *qsd1_342 = buffer.data(qsd1 + 342);
    const auto *qsd1_345 = buffer.data(qsd1 + 345);
    const auto *qsd1_347 = buffer.data(qsd1 + 347);
    const auto *qsd1_348 = buffer.data(qsd1 + 348);
    const auto *qsd1_351 = buffer.data(qsd1 + 351);
    const auto *qsd1_353 = buffer.data(qsd1 + 353);
    const auto *qsd1_354 = buffer.data(qsd1 + 354);
    const auto *qsd1_357 = buffer.data(qsd1 + 357);
    const auto *qsd1_359 = buffer.data(qsd1 + 359);
    const auto *qsd1_360 = buffer.data(qsd1 + 360);
    const auto *qsd1_363 = buffer.data(qsd1 + 363);
    const auto *qsd1_365 = buffer.data(qsd1 + 365);
    const auto *qsd1_366 = buffer.data(qsd1 + 366);
    const auto *qsd1_369 = buffer.data(qsd1 + 369);
    const auto *qsd1_371 = buffer.data(qsd1 + 371);
    const auto *qsd1_372 = buffer.data(qsd1 + 372);
    const auto *qsd1_375 = buffer.data(qsd1 + 375);
    const auto *qsd1_377 = buffer.data(qsd1 + 377);
    const auto *qsd1_378 = buffer.data(qsd1 + 378);

    const auto *qsf_556 = buffer.data(qsf + 556);
    const auto *qsf_557 = buffer.data(qsf + 557);
    const auto *qsf_559 = buffer.data(qsf + 559);
    const auto *qsf_560 = buffer.data(qsf + 560);
    const auto *qsf_562 = buffer.data(qsf + 562);
    const auto *qsf_565 = buffer.data(qsf + 565);
    const auto *qsf_566 = buffer.data(qsf + 566);
    const auto *qsf_567 = buffer.data(qsf + 567);
    const auto *qsf_568 = buffer.data(qsf + 568);
    const auto *qsf_569 = buffer.data(qsf + 569);
    const auto *qsf_570 = buffer.data(qsf + 570);
    const auto *qsf_572 = buffer.data(qsf + 572);
    const auto *qsf_573 = buffer.data(qsf + 573);
    const auto *qsf_575 = buffer.data(qsf + 575);
    const auto *qsf_576 = buffer.data(qsf + 576);
    const auto *qsf_577 = buffer.data(qsf + 577);
    const auto *qsf_578 = buffer.data(qsf + 578);
    const auto *qsf_579 = buffer.data(qsf + 579);
    const auto *qsf_580 = buffer.data(qsf + 580);
    const auto *qsf_582 = buffer.data(qsf + 582);
    const auto *qsf_583 = buffer.data(qsf + 583);
    const auto *qsf_585 = buffer.data(qsf + 585);
    const auto *qsf_586 = buffer.data(qsf + 586);
    const auto *qsf_587 = buffer.data(qsf + 587);
    const auto *qsf_588 = buffer.data(qsf + 588);
    const auto *qsf_589 = buffer.data(qsf + 589);
    const auto *qsf_590 = buffer.data(qsf + 590);
    const auto *qsf_592 = buffer.data(qsf + 592);
    const auto *qsf_593 = buffer.data(qsf + 593);
    const auto *qsf_595 = buffer.data(qsf + 595);
    const auto *qsf_596 = buffer.data(qsf + 596);
    const auto *qsf_597 = buffer.data(qsf + 597);
    const auto *qsf_598 = buffer.data(qsf + 598);
    const auto *qsf_599 = buffer.data(qsf + 599);
    const auto *qsf_600 = buffer.data(qsf + 600);
    const auto *qsf_602 = buffer.data(qsf + 602);
    const auto *qsf_603 = buffer.data(qsf + 603);
    const auto *qsf_605 = buffer.data(qsf + 605);
    const auto *qsf_606 = buffer.data(qsf + 606);
    const auto *qsf_607 = buffer.data(qsf + 607);
    const auto *qsf_608 = buffer.data(qsf + 608);
    const auto *qsf_609 = buffer.data(qsf + 609);
    const auto *qsf_610 = buffer.data(qsf + 610);
    const auto *qsf_612 = buffer.data(qsf + 612);
    const auto *qsf_613 = buffer.data(qsf + 613);
    const auto *qsf_615 = buffer.data(qsf + 615);
    const auto *qsf_616 = buffer.data(qsf + 616);
    const auto *qsf_617 = buffer.data(qsf + 617);
    const auto *qsf_618 = buffer.data(qsf + 618);
    const auto *qsf_619 = buffer.data(qsf + 619);
    const auto *qsf_620 = buffer.data(qsf + 620);
    const auto *qsf_622 = buffer.data(qsf + 622);
    const auto *qsf_623 = buffer.data(qsf + 623);
    const auto *qsf_625 = buffer.data(qsf + 625);
    const auto *qsf_626 = buffer.data(qsf + 626);
    const auto *qsf_627 = buffer.data(qsf + 627);
    const auto *qsf_628 = buffer.data(qsf + 628);
    const auto *qsf_629 = buffer.data(qsf + 629);
    const auto *qsf_630 = buffer.data(qsf + 630);

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, pc_x, pc_y, pc_z, osf_456, \
                         osf_459, osf_559, qsd0_333, qsd1_333, qsf_556, qsf_557, \
                         qsf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_8 * osf_559[k]
                   + f_3 * pc_x[k] * qsf_559[k];

        t_835[k] = f_12 * osf_456[k]
                   + f_1 * qsd0_333[k]
                   - f_2 * qsd1_333[k]
                   + f_3 * pc_y[k] * qsf_556[k];

        t_836[k] = f_3 * pc_z[k] * qsf_556[k];

        t_837[k] = f_4 * qsd0_333[k]
                   - f_5 * qsd1_333[k]
                   + f_3 * pc_z[k] * qsf_557[k];

        t_838[k] = f_12 * osf_459[k]
                   + f_3 * pc_y[k] * qsf_559[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, t_842, pa_z, pc_y, pc_z, osg0_675, osf_450, \
                         osf_460, osg1_675, qsd0_335, qsd1_335, qsf_559, \
                         qsf_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_1 * qsd0_335[k]
                   - f_2 * qsd1_335[k]
                   + f_3 * pc_z[k] * qsf_559[k];

        t_840[k] = pa_z[k] * osg0_675[k]
                   - f_6 * pc_z[k] * osg1_675[k];

        t_841[k] = f_13 * osf_460[k]
                   + f_3 * pc_y[k] * qsf_560[k];

        t_842[k] = f_7 * osf_450[k]
                   + f_3 * pc_z[k] * qsf_560[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pa_z, pc_x, pc_y, pc_z, osg0_678, osf_462, \
                         osf_565, osg1_678, qsd0_341, qsd1_341, qsf_562, \
                         qsf_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = pa_z[k] * osg0_678[k]
                   - f_6 * pc_z[k] * osg1_678[k];

        t_844[k] = f_13 * osf_462[k]
                   + f_3 * pc_y[k] * qsf_562[k];

        t_845[k] = f_8 * osf_565[k]
                   + f_4 * qsd0_341[k]
                   - f_5 * qsd1_341[k]
                   + f_3 * pc_x[k] * qsf_565[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, t_849, pc_x, osf_566, osf_567, osf_568, osf_569, \
                         qsf_566, qsf_567, qsf_568, qsf_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_8 * osf_566[k]
                   + f_3 * pc_x[k] * qsf_566[k];

        t_847[k] = f_8 * osf_567[k]
                   + f_3 * pc_x[k] * qsf_567[k];

        t_848[k] = f_8 * osf_568[k]
                   + f_3 * pc_x[k] * qsf_568[k];

        t_849[k] = f_8 * osf_569[k]
                   + f_3 * pc_x[k] * qsf_569[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, pa_z, pc_y, pc_z, osg0_685, osf_456, osf_468, \
                         osg1_685, qsd0_341, qsd1_341, qsf_566, \
                         qsf_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = pa_z[k] * osg0_685[k]
                   - f_6 * pc_z[k] * osg1_685[k];

        t_851[k] = f_7 * osf_456[k]
                   + f_3 * pc_z[k] * qsf_566[k];

        t_852[k] = f_13 * osf_468[k]
                   + f_4 * qsd0_341[k]
                   - f_5 * qsd1_341[k]
                   + f_3 * pc_y[k] * qsf_568[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, pc_x, pc_y, pc_z, osf_459, osf_469, osf_570, \
                         qsd0_341, qsd0_342, qsd1_341, qsd1_342, qsf_569, \
                         qsf_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_13 * osf_469[k]
                   + f_3 * pc_y[k] * qsf_569[k];

        t_854[k] = f_7 * osf_459[k]
                   + f_1 * qsd0_341[k]
                   - f_2 * qsd1_341[k]
                   + f_3 * pc_z[k] * qsf_569[k];

        t_855[k] = f_8 * osf_570[k]
                   + f_1 * qsd0_342[k]
                   - f_2 * qsd1_342[k]
                   + f_3 * pc_x[k] * qsf_570[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, t_859, pc_x, pc_y, pc_z, osf_460, osf_470, \
                         osf_472, osf_573, qsd0_345, qsd1_345, qsf_570, qsf_572, \
                         qsf_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_15 * osf_470[k]
                   + f_3 * pc_y[k] * qsf_570[k];

        t_857[k] = f_8 * osf_460[k]
                   + f_3 * pc_z[k] * qsf_570[k];

        t_858[k] = f_8 * osf_573[k]
                   + f_4 * qsd0_345[k]
                   - f_5 * qsd1_345[k]
                   + f_3 * pc_x[k] * qsf_573[k];

        t_859[k] = f_15 * osf_472[k]
                   + f_3 * pc_y[k] * qsf_572[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, pc_x, osf_575, osf_576, osf_577, osf_578, \
                         qsd0_347, qsd1_347, qsf_575, qsf_576, qsf_577, \
                         qsf_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_8 * osf_575[k]
                   + f_4 * qsd0_347[k]
                   - f_5 * qsd1_347[k]
                   + f_3 * pc_x[k] * qsf_575[k];

        t_861[k] = f_8 * osf_576[k]
                   + f_3 * pc_x[k] * qsf_576[k];

        t_862[k] = f_8 * osf_577[k]
                   + f_3 * pc_x[k] * qsf_577[k];

        t_863[k] = f_8 * osf_578[k]
                   + f_3 * pc_x[k] * qsf_578[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pc_x, pc_y, pc_z, osf_466, osf_476, osf_579, \
                         qsd0_345, qsd1_345, qsf_576, qsf_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_8 * osf_579[k]
                   + f_3 * pc_x[k] * qsf_579[k];

        t_865[k] = f_15 * osf_476[k]
                   + f_1 * qsd0_345[k]
                   - f_2 * qsd1_345[k]
                   + f_3 * pc_y[k] * qsf_576[k];

        t_866[k] = f_8 * osf_466[k]
                   + f_3 * pc_z[k] * qsf_576[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pc_y, pc_z, osf_469, osf_478, osf_479, qsd0_347, \
                         qsd1_347, qsf_578, qsf_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = f_15 * osf_478[k]
                   + f_4 * qsd0_347[k]
                   - f_5 * qsd1_347[k]
                   + f_3 * pc_y[k] * qsf_578[k];

        t_868[k] = f_15 * osf_479[k]
                   + f_3 * pc_y[k] * qsf_579[k];

        t_869[k] = f_8 * osf_469[k]
                   + f_1 * qsd0_347[k]
                   - f_2 * qsd1_347[k]
                   + f_3 * pc_z[k] * qsf_579[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pc_x, pc_y, pc_z, osf_470, osf_480, osf_580, \
                         qsd0_348, qsd1_348, qsf_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_8 * osf_580[k]
                   + f_1 * qsd0_348[k]
                   - f_2 * qsd1_348[k]
                   + f_3 * pc_x[k] * qsf_580[k];

        t_871[k] = f_17 * osf_480[k]
                   + f_3 * pc_y[k] * qsf_580[k];

        t_872[k] = f_14 * osf_470[k]
                   + f_3 * pc_z[k] * qsf_580[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_y, osf_482, osf_583, osf_585, qsd0_351, \
                         qsd0_353, qsd1_351, qsd1_353, qsf_582, qsf_583, \
                         qsf_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_8 * osf_583[k]
                   + f_4 * qsd0_351[k]
                   - f_5 * qsd1_351[k]
                   + f_3 * pc_x[k] * qsf_583[k];

        t_874[k] = f_17 * osf_482[k]
                   + f_3 * pc_y[k] * qsf_582[k];

        t_875[k] = f_8 * osf_585[k]
                   + f_4 * qsd0_353[k]
                   - f_5 * qsd1_353[k]
                   + f_3 * pc_x[k] * qsf_585[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, t_879, pc_x, osf_586, osf_587, osf_588, osf_589, \
                         qsf_586, qsf_587, qsf_588, qsf_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_8 * osf_586[k]
                   + f_3 * pc_x[k] * qsf_586[k];

        t_877[k] = f_8 * osf_587[k]
                   + f_3 * pc_x[k] * qsf_587[k];

        t_878[k] = f_8 * osf_588[k]
                   + f_3 * pc_x[k] * qsf_588[k];

        t_879[k] = f_8 * osf_589[k]
                   + f_3 * pc_x[k] * qsf_589[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_y, pc_z, osf_476, osf_486, osf_488, qsd0_351, \
                         qsd0_353, qsd1_351, qsd1_353, qsf_586, \
                         qsf_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_17 * osf_486[k]
                   + f_1 * qsd0_351[k]
                   - f_2 * qsd1_351[k]
                   + f_3 * pc_y[k] * qsf_586[k];

        t_881[k] = f_14 * osf_476[k]
                   + f_3 * pc_z[k] * qsf_586[k];

        t_882[k] = f_17 * osf_488[k]
                   + f_4 * qsd0_353[k]
                   - f_5 * qsd1_353[k]
                   + f_3 * pc_y[k] * qsf_588[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, pc_x, pc_y, pc_z, osf_479, osf_489, osf_590, \
                         qsd0_353, qsd0_354, qsd1_353, qsd1_354, qsf_589, \
                         qsf_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_17 * osf_489[k]
                   + f_3 * pc_y[k] * qsf_589[k];

        t_884[k] = f_14 * osf_479[k]
                   + f_1 * qsd0_353[k]
                   - f_2 * qsd1_353[k]
                   + f_3 * pc_z[k] * qsf_589[k];

        t_885[k] = f_8 * osf_590[k]
                   + f_1 * qsd0_354[k]
                   - f_2 * qsd1_354[k]
                   + f_3 * pc_x[k] * qsf_590[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, pc_z, osf_480, osf_490, \
                         osf_492, osf_593, qsd0_357, qsd1_357, qsf_590, qsf_592, \
                         qsf_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_19 * osf_490[k]
                   + f_3 * pc_y[k] * qsf_590[k];

        t_887[k] = f_16 * osf_480[k]
                   + f_3 * pc_z[k] * qsf_590[k];

        t_888[k] = f_8 * osf_593[k]
                   + f_4 * qsd0_357[k]
                   - f_5 * qsd1_357[k]
                   + f_3 * pc_x[k] * qsf_593[k];

        t_889[k] = f_19 * osf_492[k]
                   + f_3 * pc_y[k] * qsf_592[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pc_x, osf_595, osf_596, osf_597, osf_598, \
                         qsd0_359, qsd1_359, qsf_595, qsf_596, qsf_597, \
                         qsf_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_8 * osf_595[k]
                   + f_4 * qsd0_359[k]
                   - f_5 * qsd1_359[k]
                   + f_3 * pc_x[k] * qsf_595[k];

        t_891[k] = f_8 * osf_596[k]
                   + f_3 * pc_x[k] * qsf_596[k];

        t_892[k] = f_8 * osf_597[k]
                   + f_3 * pc_x[k] * qsf_597[k];

        t_893[k] = f_8 * osf_598[k]
                   + f_3 * pc_x[k] * qsf_598[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, pc_x, pc_y, pc_z, osf_486, osf_496, osf_599, \
                         qsd0_357, qsd1_357, qsf_596, qsf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_8 * osf_599[k]
                   + f_3 * pc_x[k] * qsf_599[k];

        t_895[k] = f_19 * osf_496[k]
                   + f_1 * qsd0_357[k]
                   - f_2 * qsd1_357[k]
                   + f_3 * pc_y[k] * qsf_596[k];

        t_896[k] = f_16 * osf_486[k]
                   + f_3 * pc_z[k] * qsf_596[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pc_y, pc_z, osf_489, osf_498, osf_499, qsd0_359, \
                         qsd1_359, qsf_598, qsf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_19 * osf_498[k]
                   + f_4 * qsd0_359[k]
                   - f_5 * qsd1_359[k]
                   + f_3 * pc_y[k] * qsf_598[k];

        t_898[k] = f_19 * osf_499[k]
                   + f_3 * pc_y[k] * qsf_599[k];

        t_899[k] = f_16 * osf_489[k]
                   + f_1 * qsd0_359[k]
                   - f_2 * qsd1_359[k]
                   + f_3 * pc_z[k] * qsf_599[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, pc_x, pc_y, pc_z, osf_490, osf_500, osf_600, \
                         qsd0_360, qsd1_360, qsf_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_8 * osf_600[k]
                   + f_1 * qsd0_360[k]
                   - f_2 * qsd1_360[k]
                   + f_3 * pc_x[k] * qsf_600[k];

        t_901[k] = f_18 * osf_500[k]
                   + f_3 * pc_y[k] * qsf_600[k];

        t_902[k] = f_18 * osf_490[k]
                   + f_3 * pc_z[k] * qsf_600[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, osf_502, osf_603, osf_605, qsd0_363, \
                         qsd0_365, qsd1_363, qsd1_365, qsf_602, qsf_603, \
                         qsf_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_8 * osf_603[k]
                   + f_4 * qsd0_363[k]
                   - f_5 * qsd1_363[k]
                   + f_3 * pc_x[k] * qsf_603[k];

        t_904[k] = f_18 * osf_502[k]
                   + f_3 * pc_y[k] * qsf_602[k];

        t_905[k] = f_8 * osf_605[k]
                   + f_4 * qsd0_365[k]
                   - f_5 * qsd1_365[k]
                   + f_3 * pc_x[k] * qsf_605[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, pc_x, osf_606, osf_607, osf_608, osf_609, \
                         qsf_606, qsf_607, qsf_608, qsf_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_8 * osf_606[k]
                   + f_3 * pc_x[k] * qsf_606[k];

        t_907[k] = f_8 * osf_607[k]
                   + f_3 * pc_x[k] * qsf_607[k];

        t_908[k] = f_8 * osf_608[k]
                   + f_3 * pc_x[k] * qsf_608[k];

        t_909[k] = f_8 * osf_609[k]
                   + f_3 * pc_x[k] * qsf_609[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, pc_y, pc_z, osf_496, osf_506, osf_508, qsd0_363, \
                         qsd0_365, qsd1_363, qsd1_365, qsf_606, \
                         qsf_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = f_18 * osf_506[k]
                   + f_1 * qsd0_363[k]
                   - f_2 * qsd1_363[k]
                   + f_3 * pc_y[k] * qsf_606[k];

        t_911[k] = f_18 * osf_496[k]
                   + f_3 * pc_z[k] * qsf_606[k];

        t_912[k] = f_18 * osf_508[k]
                   + f_4 * qsd0_365[k]
                   - f_5 * qsd1_365[k]
                   + f_3 * pc_y[k] * qsf_608[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, pc_x, pc_y, pc_z, osf_499, osf_509, osf_610, \
                         qsd0_365, qsd0_366, qsd1_365, qsd1_366, qsf_609, \
                         qsf_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_18 * osf_509[k]
                   + f_3 * pc_y[k] * qsf_609[k];

        t_914[k] = f_18 * osf_499[k]
                   + f_1 * qsd0_365[k]
                   - f_2 * qsd1_365[k]
                   + f_3 * pc_z[k] * qsf_609[k];

        t_915[k] = f_8 * osf_610[k]
                   + f_1 * qsd0_366[k]
                   - f_2 * qsd1_366[k]
                   + f_3 * pc_x[k] * qsf_610[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, t_919, pc_x, pc_y, pc_z, osf_500, osf_510, \
                         osf_512, osf_613, qsd0_369, qsd1_369, qsf_610, qsf_612, \
                         qsf_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_16 * osf_510[k]
                   + f_3 * pc_y[k] * qsf_610[k];

        t_917[k] = f_19 * osf_500[k]
                   + f_3 * pc_z[k] * qsf_610[k];

        t_918[k] = f_8 * osf_613[k]
                   + f_4 * qsd0_369[k]
                   - f_5 * qsd1_369[k]
                   + f_3 * pc_x[k] * qsf_613[k];

        t_919[k] = f_16 * osf_512[k]
                   + f_3 * pc_y[k] * qsf_612[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, osf_615, osf_616, osf_617, osf_618, \
                         qsd0_371, qsd1_371, qsf_615, qsf_616, qsf_617, \
                         qsf_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_8 * osf_615[k]
                   + f_4 * qsd0_371[k]
                   - f_5 * qsd1_371[k]
                   + f_3 * pc_x[k] * qsf_615[k];

        t_921[k] = f_8 * osf_616[k]
                   + f_3 * pc_x[k] * qsf_616[k];

        t_922[k] = f_8 * osf_617[k]
                   + f_3 * pc_x[k] * qsf_617[k];

        t_923[k] = f_8 * osf_618[k]
                   + f_3 * pc_x[k] * qsf_618[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_x, pc_y, pc_z, osf_506, osf_516, osf_619, \
                         qsd0_369, qsd1_369, qsf_616, qsf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_8 * osf_619[k]
                   + f_3 * pc_x[k] * qsf_619[k];

        t_925[k] = f_16 * osf_516[k]
                   + f_1 * qsd0_369[k]
                   - f_2 * qsd1_369[k]
                   + f_3 * pc_y[k] * qsf_616[k];

        t_926[k] = f_19 * osf_506[k]
                   + f_3 * pc_z[k] * qsf_616[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, pc_y, pc_z, osf_509, osf_518, osf_519, qsd0_371, \
                         qsd1_371, qsf_618, qsf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_16 * osf_518[k]
                   + f_4 * qsd0_371[k]
                   - f_5 * qsd1_371[k]
                   + f_3 * pc_y[k] * qsf_618[k];

        t_928[k] = f_16 * osf_519[k]
                   + f_3 * pc_y[k] * qsf_619[k];

        t_929[k] = f_19 * osf_509[k]
                   + f_1 * qsd0_371[k]
                   - f_2 * qsd1_371[k]
                   + f_3 * pc_z[k] * qsf_619[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, pc_x, pc_y, pc_z, osf_510, osf_520, osf_620, \
                         qsd0_372, qsd1_372, qsf_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_8 * osf_620[k]
                   + f_1 * qsd0_372[k]
                   - f_2 * qsd1_372[k]
                   + f_3 * pc_x[k] * qsf_620[k];

        t_931[k] = f_14 * osf_520[k]
                   + f_3 * pc_y[k] * qsf_620[k];

        t_932[k] = f_17 * osf_510[k]
                   + f_3 * pc_z[k] * qsf_620[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, pc_x, pc_y, osf_522, osf_623, osf_625, qsd0_375, \
                         qsd0_377, qsd1_375, qsd1_377, qsf_622, qsf_623, \
                         qsf_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_8 * osf_623[k]
                   + f_4 * qsd0_375[k]
                   - f_5 * qsd1_375[k]
                   + f_3 * pc_x[k] * qsf_623[k];

        t_934[k] = f_14 * osf_522[k]
                   + f_3 * pc_y[k] * qsf_622[k];

        t_935[k] = f_8 * osf_625[k]
                   + f_4 * qsd0_377[k]
                   - f_5 * qsd1_377[k]
                   + f_3 * pc_x[k] * qsf_625[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pc_x, osf_626, osf_627, osf_628, osf_629, \
                         qsf_626, qsf_627, qsf_628, qsf_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_8 * osf_626[k]
                   + f_3 * pc_x[k] * qsf_626[k];

        t_937[k] = f_8 * osf_627[k]
                   + f_3 * pc_x[k] * qsf_627[k];

        t_938[k] = f_8 * osf_628[k]
                   + f_3 * pc_x[k] * qsf_628[k];

        t_939[k] = f_8 * osf_629[k]
                   + f_3 * pc_x[k] * qsf_629[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, pc_y, pc_z, osf_516, osf_526, osf_528, qsd0_375, \
                         qsd0_377, qsd1_375, qsd1_377, qsf_626, \
                         qsf_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_14 * osf_526[k]
                   + f_1 * qsd0_375[k]
                   - f_2 * qsd1_375[k]
                   + f_3 * pc_y[k] * qsf_626[k];

        t_941[k] = f_17 * osf_516[k]
                   + f_3 * pc_z[k] * qsf_626[k];

        t_942[k] = f_14 * osf_528[k]
                   + f_4 * qsd0_377[k]
                   - f_5 * qsd1_377[k]
                   + f_3 * pc_y[k] * qsf_628[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, pc_x, pc_y, pc_z, osf_519, osf_529, osf_630, \
                         qsd0_377, qsd0_378, qsd1_377, qsd1_378, qsf_629, \
                         qsf_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_14 * osf_529[k]
                   + f_3 * pc_y[k] * qsf_629[k];

        t_944[k] = f_17 * osf_519[k]
                   + f_1 * qsd0_377[k]
                   - f_2 * qsd1_377[k]
                   + f_3 * pc_z[k] * qsf_629[k];

        t_945[k] = f_8 * osf_630[k]
                   + f_1 * qsd0_378[k]
                   - f_2 * qsd1_378[k]
                   + f_3 * pc_x[k] * qsf_630[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_810 = buffer.data(osg0 + 810);
    const auto *osg0_813 = buffer.data(osg0 + 813);
    const auto *osg0_815 = buffer.data(osg0 + 815);
    const auto *osg0_824 = buffer.data(osg0 + 824);
    const auto *osg0_825 = buffer.data(osg0 + 825);
    const auto *osg0_828 = buffer.data(osg0 + 828);
    const auto *osg0_990 = buffer.data(osg0 + 990);
    const auto *osg0_993 = buffer.data(osg0 + 993);
    const auto *osg0_1000 = buffer.data(osg0 + 1000);
    const auto *osg0_1002 = buffer.data(osg0 + 1002);
    const auto *osg0_1004 = buffer.data(osg0 + 1004);
    const auto *osg0_1010 = buffer.data(osg0 + 1010);
    const auto *osg0_1015 = buffer.data(osg0 + 1015);
    const auto *osg0_1017 = buffer.data(osg0 + 1017);
    const auto *osg0_1019 = buffer.data(osg0 + 1019);
    const auto *osg0_1020 = buffer.data(osg0 + 1020);
    const auto *osg0_1023 = buffer.data(osg0 + 1023);
    const auto *osg0_1025 = buffer.data(osg0 + 1025);
    const auto *osg0_1030 = buffer.data(osg0 + 1030);
    const auto *osg0_1032 = buffer.data(osg0 + 1032);
    const auto *osg0_1034 = buffer.data(osg0 + 1034);
    const auto *osg0_1035 = buffer.data(osg0 + 1035);
    const auto *osg0_1038 = buffer.data(osg0 + 1038);
    const auto *osg0_1040 = buffer.data(osg0 + 1040);
    const auto *osg0_1045 = buffer.data(osg0 + 1045);
    const auto *osg0_1047 = buffer.data(osg0 + 1047);
    const auto *osg0_1049 = buffer.data(osg0 + 1049);
    const auto *osg0_1050 = buffer.data(osg0 + 1050);
    const auto *osg0_1053 = buffer.data(osg0 + 1053);
    const auto *osg0_1055 = buffer.data(osg0 + 1055);
    const auto *osg0_1060 = buffer.data(osg0 + 1060);
    const auto *osg0_1062 = buffer.data(osg0 + 1062);
    const auto *osg0_1064 = buffer.data(osg0 + 1064);
    const auto *osg0_1065 = buffer.data(osg0 + 1065);
    const auto *osg0_1068 = buffer.data(osg0 + 1068);

    const auto *osf_520 = buffer.data(osf + 520);
    const auto *osf_526 = buffer.data(osf + 526);
    const auto *osf_529 = buffer.data(osf + 529);
    const auto *osf_530 = buffer.data(osf + 530);
    const auto *osf_532 = buffer.data(osf + 532);
    const auto *osf_536 = buffer.data(osf + 536);
    const auto *osf_538 = buffer.data(osf + 538);
    const auto *osf_539 = buffer.data(osf + 539);
    const auto *osf_540 = buffer.data(osf + 540);
    const auto *osf_541 = buffer.data(osf + 541);
    const auto *osf_542 = buffer.data(osf + 542);
    const auto *osf_546 = buffer.data(osf + 546);
    const auto *osf_548 = buffer.data(osf + 548);
    const auto *osf_549 = buffer.data(osf + 549);
    const auto *osf_550 = buffer.data(osf + 550);
    const auto *osf_556 = buffer.data(osf + 556);
    const auto *osf_559 = buffer.data(osf + 559);
    const auto *osf_560 = buffer.data(osf + 560);
    const auto *osf_562 = buffer.data(osf + 562);
    const auto *osf_566 = buffer.data(osf + 566);
    const auto *osf_569 = buffer.data(osf + 569);
    const auto *osf_570 = buffer.data(osf + 570);
    const auto *osf_572 = buffer.data(osf + 572);
    const auto *osf_576 = buffer.data(osf + 576);
    const auto *osf_579 = buffer.data(osf + 579);
    const auto *osf_580 = buffer.data(osf + 580);
    const auto *osf_582 = buffer.data(osf + 582);
    const auto *osf_586 = buffer.data(osf + 586);
    const auto *osf_589 = buffer.data(osf + 589);
    const auto *osf_590 = buffer.data(osf + 590);
    const auto *osf_592 = buffer.data(osf + 592);
    const auto *osf_599 = buffer.data(osf + 599);
    const auto *osf_600 = buffer.data(osf + 600);
    const auto *osf_602 = buffer.data(osf + 602);
    const auto *osf_633 = buffer.data(osf + 633);
    const auto *osf_635 = buffer.data(osf + 635);
    const auto *osf_636 = buffer.data(osf + 636);
    const auto *osf_637 = buffer.data(osf + 637);
    const auto *osf_638 = buffer.data(osf + 638);
    const auto *osf_639 = buffer.data(osf + 639);
    const auto *osf_646 = buffer.data(osf + 646);
    const auto *osf_647 = buffer.data(osf + 647);
    const auto *osf_648 = buffer.data(osf + 648);
    const auto *osf_649 = buffer.data(osf + 649);
    const auto *osf_650 = buffer.data(osf + 650);
    const auto *osf_655 = buffer.data(osf + 655);
    const auto *osf_656 = buffer.data(osf + 656);
    const auto *osf_657 = buffer.data(osf + 657);
    const auto *osf_659 = buffer.data(osf + 659);
    const auto *osf_660 = buffer.data(osf + 660);
    const auto *osf_663 = buffer.data(osf + 663);
    const auto *osf_666 = buffer.data(osf + 666);
    const auto *osf_668 = buffer.data(osf + 668);
    const auto *osf_669 = buffer.data(osf + 669);
    const auto *osf_675 = buffer.data(osf + 675);
    const auto *osf_676 = buffer.data(osf + 676);
    const auto *osf_677 = buffer.data(osf + 677);
    const auto *osf_678 = buffer.data(osf + 678);
    const auto *osf_679 = buffer.data(osf + 679);
    const auto *osf_680 = buffer.data(osf + 680);
    const auto *osf_683 = buffer.data(osf + 683);
    const auto *osf_685 = buffer.data(osf + 685);
    const auto *osf_686 = buffer.data(osf + 686);
    const auto *osf_687 = buffer.data(osf + 687);
    const auto *osf_688 = buffer.data(osf + 688);
    const auto *osf_689 = buffer.data(osf + 689);
    const auto *osf_690 = buffer.data(osf + 690);
    const auto *osf_693 = buffer.data(osf + 693);
    const auto *osf_695 = buffer.data(osf + 695);
    const auto *osf_696 = buffer.data(osf + 696);
    const auto *osf_697 = buffer.data(osf + 697);
    const auto *osf_698 = buffer.data(osf + 698);
    const auto *osf_699 = buffer.data(osf + 699);
    const auto *osf_700 = buffer.data(osf + 700);
    const auto *osf_703 = buffer.data(osf + 703);
    const auto *osf_705 = buffer.data(osf + 705);
    const auto *osf_706 = buffer.data(osf + 706);
    const auto *osf_707 = buffer.data(osf + 707);
    const auto *osf_708 = buffer.data(osf + 708);
    const auto *osf_709 = buffer.data(osf + 709);
    const auto *osf_710 = buffer.data(osf + 710);
    const auto *osf_713 = buffer.data(osf + 713);

    const auto *osg1_810 = buffer.data(osg1 + 810);
    const auto *osg1_813 = buffer.data(osg1 + 813);
    const auto *osg1_815 = buffer.data(osg1 + 815);
    const auto *osg1_824 = buffer.data(osg1 + 824);
    const auto *osg1_825 = buffer.data(osg1 + 825);
    const auto *osg1_828 = buffer.data(osg1 + 828);
    const auto *osg1_990 = buffer.data(osg1 + 990);
    const auto *osg1_993 = buffer.data(osg1 + 993);
    const auto *osg1_1000 = buffer.data(osg1 + 1000);
    const auto *osg1_1002 = buffer.data(osg1 + 1002);
    const auto *osg1_1004 = buffer.data(osg1 + 1004);
    const auto *osg1_1010 = buffer.data(osg1 + 1010);
    const auto *osg1_1015 = buffer.data(osg1 + 1015);
    const auto *osg1_1017 = buffer.data(osg1 + 1017);
    const auto *osg1_1019 = buffer.data(osg1 + 1019);
    const auto *osg1_1020 = buffer.data(osg1 + 1020);
    const auto *osg1_1023 = buffer.data(osg1 + 1023);
    const auto *osg1_1025 = buffer.data(osg1 + 1025);
    const auto *osg1_1030 = buffer.data(osg1 + 1030);
    const auto *osg1_1032 = buffer.data(osg1 + 1032);
    const auto *osg1_1034 = buffer.data(osg1 + 1034);
    const auto *osg1_1035 = buffer.data(osg1 + 1035);
    const auto *osg1_1038 = buffer.data(osg1 + 1038);
    const auto *osg1_1040 = buffer.data(osg1 + 1040);
    const auto *osg1_1045 = buffer.data(osg1 + 1045);
    const auto *osg1_1047 = buffer.data(osg1 + 1047);
    const auto *osg1_1049 = buffer.data(osg1 + 1049);
    const auto *osg1_1050 = buffer.data(osg1 + 1050);
    const auto *osg1_1053 = buffer.data(osg1 + 1053);
    const auto *osg1_1055 = buffer.data(osg1 + 1055);
    const auto *osg1_1060 = buffer.data(osg1 + 1060);
    const auto *osg1_1062 = buffer.data(osg1 + 1062);
    const auto *osg1_1064 = buffer.data(osg1 + 1064);
    const auto *osg1_1065 = buffer.data(osg1 + 1065);
    const auto *osg1_1068 = buffer.data(osg1 + 1068);

    const auto *qsd0_381 = buffer.data(qsd0 + 381);
    const auto *qsd0_383 = buffer.data(qsd0 + 383);
    const auto *qsd0_387 = buffer.data(qsd0 + 387);
    const auto *qsd0_389 = buffer.data(qsd0 + 389);
    const auto *qsd0_390 = buffer.data(qsd0 + 390);
    const auto *qsd0_393 = buffer.data(qsd0 + 393);
    const auto *qsd0_394 = buffer.data(qsd0 + 394);
    const auto *qsd0_395 = buffer.data(qsd0 + 395);
    const auto *qsd0_396 = buffer.data(qsd0 + 396);

    const auto *qsd1_381 = buffer.data(qsd1 + 381);
    const auto *qsd1_383 = buffer.data(qsd1 + 383);
    const auto *qsd1_387 = buffer.data(qsd1 + 387);
    const auto *qsd1_389 = buffer.data(qsd1 + 389);
    const auto *qsd1_390 = buffer.data(qsd1 + 390);
    const auto *qsd1_393 = buffer.data(qsd1 + 393);
    const auto *qsd1_394 = buffer.data(qsd1 + 394);
    const auto *qsd1_395 = buffer.data(qsd1 + 395);
    const auto *qsd1_396 = buffer.data(qsd1 + 396);

    const auto *qsf_630 = buffer.data(qsf + 630);
    const auto *qsf_632 = buffer.data(qsf + 632);
    const auto *qsf_633 = buffer.data(qsf + 633);
    const auto *qsf_635 = buffer.data(qsf + 635);
    const auto *qsf_636 = buffer.data(qsf + 636);
    const auto *qsf_637 = buffer.data(qsf + 637);
    const auto *qsf_638 = buffer.data(qsf + 638);
    const auto *qsf_639 = buffer.data(qsf + 639);
    const auto *qsf_640 = buffer.data(qsf + 640);
    const auto *qsf_642 = buffer.data(qsf + 642);
    const auto *qsf_646 = buffer.data(qsf + 646);
    const auto *qsf_647 = buffer.data(qsf + 647);
    const auto *qsf_648 = buffer.data(qsf + 648);
    const auto *qsf_649 = buffer.data(qsf + 649);
    const auto *qsf_650 = buffer.data(qsf + 650);
    const auto *qsf_651 = buffer.data(qsf + 651);
    const auto *qsf_652 = buffer.data(qsf + 652);
    const auto *qsf_655 = buffer.data(qsf + 655);
    const auto *qsf_656 = buffer.data(qsf + 656);
    const auto *qsf_657 = buffer.data(qsf + 657);
    const auto *qsf_658 = buffer.data(qsf + 658);
    const auto *qsf_659 = buffer.data(qsf + 659);
    const auto *qsf_660 = buffer.data(qsf + 660);
    const auto *qsf_661 = buffer.data(qsf + 661);
    const auto *qsf_662 = buffer.data(qsf + 662);
    const auto *qsf_663 = buffer.data(qsf + 663);
    const auto *qsf_666 = buffer.data(qsf + 666);
    const auto *qsf_668 = buffer.data(qsf + 668);
    const auto *qsf_669 = buffer.data(qsf + 669);
    const auto *qsf_670 = buffer.data(qsf + 670);
    const auto *qsf_672 = buffer.data(qsf + 672);
    const auto *qsf_676 = buffer.data(qsf + 676);
    const auto *qsf_677 = buffer.data(qsf + 677);
    const auto *qsf_678 = buffer.data(qsf + 678);
    const auto *qsf_679 = buffer.data(qsf + 679);
    const auto *qsf_680 = buffer.data(qsf + 680);
    const auto *qsf_682 = buffer.data(qsf + 682);
    const auto *qsf_686 = buffer.data(qsf + 686);
    const auto *qsf_687 = buffer.data(qsf + 687);
    const auto *qsf_688 = buffer.data(qsf + 688);
    const auto *qsf_689 = buffer.data(qsf + 689);
    const auto *qsf_690 = buffer.data(qsf + 690);
    const auto *qsf_692 = buffer.data(qsf + 692);
    const auto *qsf_696 = buffer.data(qsf + 696);
    const auto *qsf_697 = buffer.data(qsf + 697);
    const auto *qsf_698 = buffer.data(qsf + 698);
    const auto *qsf_699 = buffer.data(qsf + 699);
    const auto *qsf_700 = buffer.data(qsf + 700);
    const auto *qsf_702 = buffer.data(qsf + 702);
    const auto *qsf_706 = buffer.data(qsf + 706);
    const auto *qsf_707 = buffer.data(qsf + 707);
    const auto *qsf_708 = buffer.data(qsf + 708);
    const auto *qsf_709 = buffer.data(qsf + 709);
    const auto *qsf_710 = buffer.data(qsf + 710);
    const auto *qsf_712 = buffer.data(qsf + 712);

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pc_x, pc_y, pc_z, osf_520, osf_530, \
                         osf_532, osf_633, qsd0_381, qsd1_381, qsf_630, qsf_632, \
                         qsf_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_8 * osf_530[k]
                   + f_3 * pc_y[k] * qsf_630[k];

        t_947[k] = f_15 * osf_520[k]
                   + f_3 * pc_z[k] * qsf_630[k];

        t_948[k] = f_8 * osf_633[k]
                   + f_4 * qsd0_381[k]
                   - f_5 * qsd1_381[k]
                   + f_3 * pc_x[k] * qsf_633[k];

        t_949[k] = f_8 * osf_532[k]
                   + f_3 * pc_y[k] * qsf_632[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pc_x, osf_635, osf_636, osf_637, osf_638, \
                         qsd0_383, qsd1_383, qsf_635, qsf_636, qsf_637, \
                         qsf_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_8 * osf_635[k]
                   + f_4 * qsd0_383[k]
                   - f_5 * qsd1_383[k]
                   + f_3 * pc_x[k] * qsf_635[k];

        t_951[k] = f_8 * osf_636[k]
                   + f_3 * pc_x[k] * qsf_636[k];

        t_952[k] = f_8 * osf_637[k]
                   + f_3 * pc_x[k] * qsf_637[k];

        t_953[k] = f_8 * osf_638[k]
                   + f_3 * pc_x[k] * qsf_638[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pc_x, pc_y, pc_z, osf_526, osf_536, osf_639, \
                         qsd0_381, qsd1_381, qsf_636, qsf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_8 * osf_639[k]
                   + f_3 * pc_x[k] * qsf_639[k];

        t_955[k] = f_8 * osf_536[k]
                   + f_1 * qsd0_381[k]
                   - f_2 * qsd1_381[k]
                   + f_3 * pc_y[k] * qsf_636[k];

        t_956[k] = f_15 * osf_526[k]
                   + f_3 * pc_z[k] * qsf_636[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_y, pc_y, pc_z, osg0_810, osf_529, \
                         osf_538, osf_539, osg1_810, qsd0_383, qsd1_383, qsf_638, \
                         qsf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_8 * osf_538[k]
                   + f_4 * qsd0_383[k]
                   - f_5 * qsd1_383[k]
                   + f_3 * pc_y[k] * qsf_638[k];

        t_958[k] = f_8 * osf_539[k]
                   + f_3 * pc_y[k] * qsf_639[k];

        t_959[k] = f_15 * osf_529[k]
                   + f_1 * qsd0_383[k]
                   - f_2 * qsd1_383[k]
                   + f_3 * pc_z[k] * qsf_639[k];

        t_960[k] = pa_y[k] * osg0_810[k]
                   - f_6 * pc_y[k] * osg1_810[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pa_y, pc_y, pc_z, osg0_813, osf_530, \
                         osf_540, osf_541, osf_542, osg1_813, qsf_640, \
                         qsf_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_7 * osf_540[k]
                   + f_3 * pc_y[k] * qsf_640[k];

        t_962[k] = f_13 * osf_530[k]
                   + f_3 * pc_z[k] * qsf_640[k];

        t_963[k] = pa_y[k] * osg0_813[k]
                   + f_8 * osf_541[k]
                   - f_6 * pc_y[k] * osg1_813[k];

        t_964[k] = f_7 * osf_542[k]
                   + f_3 * pc_y[k] * qsf_642[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pa_y, pc_x, pc_y, osg0_815, osf_646, \
                         osf_647, osf_648, osg1_815, qsf_646, qsf_647, \
                         qsf_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pa_y[k] * osg0_815[k]
                   - f_6 * pc_y[k] * osg1_815[k];

        t_966[k] = f_8 * osf_646[k]
                   + f_3 * pc_x[k] * qsf_646[k];

        t_967[k] = f_8 * osf_647[k]
                   + f_3 * pc_x[k] * qsf_647[k];

        t_968[k] = f_8 * osf_648[k]
                   + f_3 * pc_x[k] * qsf_648[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pc_x, pc_y, pc_z, osf_536, osf_546, osf_649, \
                         qsd0_387, qsd1_387, qsf_646, qsf_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_8 * osf_649[k]
                   + f_3 * pc_x[k] * qsf_649[k];

        t_970[k] = f_7 * osf_546[k]
                   + f_1 * qsd0_387[k]
                   - f_2 * qsd1_387[k]
                   + f_3 * pc_y[k] * qsf_646[k];

        t_971[k] = f_13 * osf_536[k]
                   + f_3 * pc_z[k] * qsf_646[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, pa_y, pc_y, osg0_824, osf_548, osf_549, \
                         osg1_824, qsd0_389, qsd1_389, qsf_648, \
                         qsf_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_7 * osf_548[k]
                   + f_4 * qsd0_389[k]
                   - f_5 * qsd1_389[k]
                   + f_3 * pc_y[k] * qsf_648[k];

        t_973[k] = f_7 * osf_549[k]
                   + f_3 * pc_y[k] * qsf_649[k];

        t_974[k] = pa_y[k] * osg0_824[k]
                   - f_6 * pc_y[k] * osg1_824[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, pc_x, pc_y, pc_z, osf_540, \
                         osf_650, qsd0_390, qsd1_390, qsf_650, qsf_651, \
                         qsf_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_8 * osf_650[k]
                   + f_1 * qsd0_390[k]
                   - f_2 * qsd1_390[k]
                   + f_3 * pc_x[k] * qsf_650[k];

        t_976[k] = f_3 * pc_y[k] * qsf_650[k];

        t_977[k] = f_12 * osf_540[k]
                   + f_3 * pc_z[k] * qsf_650[k];

        t_978[k] = f_4 * qsd0_390[k]
                   - f_5 * qsd1_390[k]
                   + f_3 * pc_y[k] * qsf_651[k];

        t_979[k] = f_3 * pc_y[k] * qsf_652[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, pc_x, pc_y, osf_655, osf_656, osf_657, \
                         qsd0_395, qsd1_395, qsf_655, qsf_656, \
                         qsf_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_8 * osf_655[k]
                   + f_4 * qsd0_395[k]
                   - f_5 * qsd1_395[k]
                   + f_3 * pc_x[k] * qsf_655[k];

        t_981[k] = f_8 * osf_656[k]
                   + f_3 * pc_x[k] * qsf_656[k];

        t_982[k] = f_8 * osf_657[k]
                   + f_3 * pc_x[k] * qsf_657[k];

        t_983[k] = f_3 * pc_y[k] * qsf_655[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_x, pc_y, osf_659, qsd0_393, qsd0_394, \
                         qsd1_393, qsd1_394, qsf_656, qsf_657, \
                         qsf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_8 * osf_659[k]
                   + f_3 * pc_x[k] * qsf_659[k];

        t_985[k] = f_1 * qsd0_393[k]
                   - f_2 * qsd1_393[k]
                   + f_3 * pc_y[k] * qsf_656[k];

        t_986[k] = f_10 * qsd0_394[k]
                   - f_11 * qsd1_394[k]
                   + f_3 * pc_y[k] * qsf_657[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, pa_x, pc_x, pc_y, pc_z, osg0_990, \
                         osf_549, osf_660, osg1_990, qsd0_395, qsd1_395, qsf_658, \
                         qsf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_4 * qsd0_395[k]
                   - f_5 * qsd1_395[k]
                   + f_3 * pc_y[k] * qsf_658[k];

        t_988[k] = f_3 * pc_y[k] * qsf_659[k];

        t_989[k] = f_12 * osf_549[k]
                   + f_1 * qsd0_395[k]
                   - f_2 * qsd1_395[k]
                   + f_3 * pc_z[k] * qsf_659[k];

        t_990[k] = pa_x[k] * osg0_990[k]
                   + f_16 * osf_660[k]
                   - f_6 * pc_x[k] * osg1_990[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pa_x, pc_x, pc_y, pc_z, osg0_993, \
                         osf_550, osf_663, osg1_993, qsf_660, qsf_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_9 * osf_550[k]
                   + f_3 * pc_y[k] * qsf_660[k];

        t_992[k] = f_3 * pc_z[k] * qsf_660[k];

        t_993[k] = pa_x[k] * osg0_993[k]
                   + f_8 * osf_663[k]
                   - f_6 * pc_x[k] * osg1_993[k];

        t_994[k] = f_3 * pc_z[k] * qsf_661[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pc_x, pc_z, osf_666, osf_668, qsd0_396, \
                         qsd1_396, qsf_662, qsf_663, qsf_666, qsf_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_4 * qsd0_396[k]
                   - f_5 * qsd1_396[k]
                   + f_3 * pc_z[k] * qsf_662[k];

        t_996[k] = f_7 * osf_666[k]
                   + f_3 * pc_x[k] * qsf_666[k];

        t_997[k] = f_3 * pc_z[k] * qsf_663[k];

        t_998[k] = f_7 * osf_668[k]
                   + f_3 * pc_x[k] * qsf_668[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, t_1002, pa_x, pc_x, pc_z, osg0_1000, \
                         osg0_1002, osf_669, osg1_1000, osg1_1002, qsf_666, \
                         qsf_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_7 * osf_669[k]
                   + f_3 * pc_x[k] * qsf_669[k];

        t_1000[k] = pa_x[k] * osg0_1000[k]
                    - f_6 * pc_x[k] * osg1_1000[k];

        t_1001[k] = f_3 * pc_z[k] * qsf_666[k];

        t_1002[k] = pa_x[k] * osg0_1002[k]
                    - f_6 * pc_x[k] * osg1_1002[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, pa_x, pa_z, pc_x, pc_y, pc_z, osg0_825, \
                         osg0_1004, osf_559, osg1_825, osg1_1004, \
                         qsf_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_9 * osf_559[k]
                    + f_3 * pc_y[k] * qsf_669[k];

        t_1004[k] = pa_x[k] * osg0_1004[k]
                    - f_6 * pc_x[k] * osg1_1004[k];

        t_1005[k] = pa_z[k] * osg0_825[k]
                    - f_6 * pc_z[k] * osg1_825[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, t_1009, pa_z, pc_y, pc_z, osg0_828, osf_550, \
                         osf_560, osf_562, osg1_828, qsf_670, qsf_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_12 * osf_560[k]
                    + f_3 * pc_y[k] * qsf_670[k];

        t_1007[k] = f_7 * osf_550[k]
                    + f_3 * pc_z[k] * qsf_670[k];

        t_1008[k] = pa_z[k] * osg0_828[k]
                    - f_6 * pc_z[k] * osg1_828[k];

        t_1009[k] = f_12 * osf_562[k]
                    + f_3 * pc_y[k] * qsf_672[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, t_1013, pa_x, pc_x, osg0_1010, osf_675, \
                         osf_676, osf_677, osf_678, osg1_1010, qsf_676, qsf_677, \
                         qsf_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = pa_x[k] * osg0_1010[k]
                    + f_8 * osf_675[k]
                    - f_6 * pc_x[k] * osg1_1010[k];

        t_1011[k] = f_7 * osf_676[k]
                    + f_3 * pc_x[k] * qsf_676[k];

        t_1012[k] = f_7 * osf_677[k]
                    + f_3 * pc_x[k] * qsf_677[k];

        t_1013[k] = f_7 * osf_678[k]
                    + f_3 * pc_x[k] * qsf_678[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, t_1017, pa_x, pc_x, pc_z, osg0_1015, \
                         osg0_1017, osf_556, osf_679, osg1_1015, osg1_1017, qsf_676, \
                         qsf_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_7 * osf_679[k]
                    + f_3 * pc_x[k] * qsf_679[k];

        t_1015[k] = pa_x[k] * osg0_1015[k]
                    - f_6 * pc_x[k] * osg1_1015[k];

        t_1016[k] = f_7 * osf_556[k]
                    + f_3 * pc_z[k] * qsf_676[k];

        t_1017[k] = pa_x[k] * osg0_1017[k]
                    - f_6 * pc_x[k] * osg1_1017[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, t_1021, pa_x, pc_x, pc_y, osg0_1019, \
                         osg0_1020, osf_569, osf_570, osf_680, osg1_1019, osg1_1020, qsf_679, \
                         qsf_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_12 * osf_569[k]
                    + f_3 * pc_y[k] * qsf_679[k];

        t_1019[k] = pa_x[k] * osg0_1019[k]
                    - f_6 * pc_x[k] * osg1_1019[k];

        t_1020[k] = pa_x[k] * osg0_1020[k]
                    + f_16 * osf_680[k]
                    - f_6 * pc_x[k] * osg1_1020[k];

        t_1021[k] = f_13 * osf_570[k]
                    + f_3 * pc_y[k] * qsf_680[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, pa_x, pc_x, pc_y, pc_z, osg0_1023, osf_560, \
                         osf_572, osf_683, osg1_1023, qsf_680, \
                         qsf_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = f_8 * osf_560[k]
                    + f_3 * pc_z[k] * qsf_680[k];

        t_1023[k] = pa_x[k] * osg0_1023[k]
                    + f_8 * osf_683[k]
                    - f_6 * pc_x[k] * osg1_1023[k];

        t_1024[k] = f_13 * osf_572[k]
                    + f_3 * pc_y[k] * qsf_682[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pa_x, pc_x, osg0_1025, osf_685, \
                         osf_686, osf_687, osf_688, osg1_1025, qsf_686, qsf_687, \
                         qsf_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = pa_x[k] * osg0_1025[k]
                    + f_8 * osf_685[k]
                    - f_6 * pc_x[k] * osg1_1025[k];

        t_1026[k] = f_7 * osf_686[k]
                    + f_3 * pc_x[k] * qsf_686[k];

        t_1027[k] = f_7 * osf_687[k]
                    + f_3 * pc_x[k] * qsf_687[k];

        t_1028[k] = f_7 * osf_688[k]
                    + f_3 * pc_x[k] * qsf_688[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pa_x, pc_x, pc_z, osg0_1030, \
                         osg0_1032, osf_566, osf_689, osg1_1030, osg1_1032, qsf_686, \
                         qsf_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_7 * osf_689[k]
                    + f_3 * pc_x[k] * qsf_689[k];

        t_1030[k] = pa_x[k] * osg0_1030[k]
                    - f_6 * pc_x[k] * osg1_1030[k];

        t_1031[k] = f_8 * osf_566[k]
                    + f_3 * pc_z[k] * qsf_686[k];

        t_1032[k] = pa_x[k] * osg0_1032[k]
                    - f_6 * pc_x[k] * osg1_1032[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, pa_x, pc_x, pc_y, osg0_1034, \
                         osg0_1035, osf_579, osf_580, osf_690, osg1_1034, osg1_1035, qsf_689, \
                         qsf_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_13 * osf_579[k]
                    + f_3 * pc_y[k] * qsf_689[k];

        t_1034[k] = pa_x[k] * osg0_1034[k]
                    - f_6 * pc_x[k] * osg1_1034[k];

        t_1035[k] = pa_x[k] * osg0_1035[k]
                    + f_16 * osf_690[k]
                    - f_6 * pc_x[k] * osg1_1035[k];

        t_1036[k] = f_15 * osf_580[k]
                    + f_3 * pc_y[k] * qsf_690[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, pa_x, pc_x, pc_y, pc_z, osg0_1038, osf_570, \
                         osf_582, osf_693, osg1_1038, qsf_690, \
                         qsf_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_14 * osf_570[k]
                    + f_3 * pc_z[k] * qsf_690[k];

        t_1038[k] = pa_x[k] * osg0_1038[k]
                    + f_8 * osf_693[k]
                    - f_6 * pc_x[k] * osg1_1038[k];

        t_1039[k] = f_15 * osf_582[k]
                    + f_3 * pc_y[k] * qsf_692[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, t_1043, pa_x, pc_x, osg0_1040, osf_695, \
                         osf_696, osf_697, osf_698, osg1_1040, qsf_696, qsf_697, \
                         qsf_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = pa_x[k] * osg0_1040[k]
                    + f_8 * osf_695[k]
                    - f_6 * pc_x[k] * osg1_1040[k];

        t_1041[k] = f_7 * osf_696[k]
                    + f_3 * pc_x[k] * qsf_696[k];

        t_1042[k] = f_7 * osf_697[k]
                    + f_3 * pc_x[k] * qsf_697[k];

        t_1043[k] = f_7 * osf_698[k]
                    + f_3 * pc_x[k] * qsf_698[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, t_1047, pa_x, pc_x, pc_z, osg0_1045, \
                         osg0_1047, osf_576, osf_699, osg1_1045, osg1_1047, qsf_696, \
                         qsf_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_7 * osf_699[k]
                    + f_3 * pc_x[k] * qsf_699[k];

        t_1045[k] = pa_x[k] * osg0_1045[k]
                    - f_6 * pc_x[k] * osg1_1045[k];

        t_1046[k] = f_14 * osf_576[k]
                    + f_3 * pc_z[k] * qsf_696[k];

        t_1047[k] = pa_x[k] * osg0_1047[k]
                    - f_6 * pc_x[k] * osg1_1047[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, t_1051, pa_x, pc_x, pc_y, osg0_1049, \
                         osg0_1050, osf_589, osf_590, osf_700, osg1_1049, osg1_1050, qsf_699, \
                         qsf_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_15 * osf_589[k]
                    + f_3 * pc_y[k] * qsf_699[k];

        t_1049[k] = pa_x[k] * osg0_1049[k]
                    - f_6 * pc_x[k] * osg1_1049[k];

        t_1050[k] = pa_x[k] * osg0_1050[k]
                    + f_16 * osf_700[k]
                    - f_6 * pc_x[k] * osg1_1050[k];

        t_1051[k] = f_17 * osf_590[k]
                    + f_3 * pc_y[k] * qsf_700[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pa_x, pc_x, pc_y, pc_z, osg0_1053, osf_580, \
                         osf_592, osf_703, osg1_1053, qsf_700, \
                         qsf_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_16 * osf_580[k]
                    + f_3 * pc_z[k] * qsf_700[k];

        t_1053[k] = pa_x[k] * osg0_1053[k]
                    + f_8 * osf_703[k]
                    - f_6 * pc_x[k] * osg1_1053[k];

        t_1054[k] = f_17 * osf_592[k]
                    + f_3 * pc_y[k] * qsf_702[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, t_1058, pa_x, pc_x, osg0_1055, osf_705, \
                         osf_706, osf_707, osf_708, osg1_1055, qsf_706, qsf_707, \
                         qsf_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = pa_x[k] * osg0_1055[k]
                    + f_8 * osf_705[k]
                    - f_6 * pc_x[k] * osg1_1055[k];

        t_1056[k] = f_7 * osf_706[k]
                    + f_3 * pc_x[k] * qsf_706[k];

        t_1057[k] = f_7 * osf_707[k]
                    + f_3 * pc_x[k] * qsf_707[k];

        t_1058[k] = f_7 * osf_708[k]
                    + f_3 * pc_x[k] * qsf_708[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, pa_x, pc_x, pc_z, osg0_1060, \
                         osg0_1062, osf_586, osf_709, osg1_1060, osg1_1062, qsf_706, \
                         qsf_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = f_7 * osf_709[k]
                    + f_3 * pc_x[k] * qsf_709[k];

        t_1060[k] = pa_x[k] * osg0_1060[k]
                    - f_6 * pc_x[k] * osg1_1060[k];

        t_1061[k] = f_16 * osf_586[k]
                    + f_3 * pc_z[k] * qsf_706[k];

        t_1062[k] = pa_x[k] * osg0_1062[k]
                    - f_6 * pc_x[k] * osg1_1062[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, pa_x, pc_x, pc_y, osg0_1064, \
                         osg0_1065, osf_599, osf_600, osf_710, osg1_1064, osg1_1065, qsf_709, \
                         qsf_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_17 * osf_599[k]
                    + f_3 * pc_y[k] * qsf_709[k];

        t_1064[k] = pa_x[k] * osg0_1064[k]
                    - f_6 * pc_x[k] * osg1_1064[k];

        t_1065[k] = pa_x[k] * osg0_1065[k]
                    + f_16 * osf_710[k]
                    - f_6 * pc_x[k] * osg1_1065[k];

        t_1066[k] = f_19 * osf_600[k]
                    + f_3 * pc_y[k] * qsf_710[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pa_x, pc_x, pc_y, pc_z, osg0_1068, osf_590, \
                         osf_602, osf_713, osg1_1068, qsf_710, \
                         qsf_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_18 * osf_590[k]
                    + f_3 * pc_z[k] * qsf_710[k];

        t_1068[k] = pa_x[k] * osg0_1068[k]
                    + f_8 * osf_713[k]
                    - f_6 * pc_x[k] * osg1_1068[k];

        t_1069[k] = f_19 * osf_602[k]
                    + f_3 * pc_y[k] * qsf_712[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osg0,
                                                          const size_t osf, const size_t osg1,
                                                          const size_t qsd0, const size_t qsd1,
                                                          const size_t qsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_975 = buffer.data(osg0 + 975);
    const auto *osg0_980 = buffer.data(osg0 + 980);
    const auto *osg0_990 = buffer.data(osg0 + 990);
    const auto *osg0_991 = buffer.data(osg0 + 991);
    const auto *osg0_993 = buffer.data(osg0 + 993);
    const auto *osg0_1000 = buffer.data(osg0 + 1000);
    const auto *osg0_1002 = buffer.data(osg0 + 1002);
    const auto *osg0_1070 = buffer.data(osg0 + 1070);
    const auto *osg0_1075 = buffer.data(osg0 + 1075);
    const auto *osg0_1077 = buffer.data(osg0 + 1077);
    const auto *osg0_1079 = buffer.data(osg0 + 1079);
    const auto *osg0_1080 = buffer.data(osg0 + 1080);
    const auto *osg0_1083 = buffer.data(osg0 + 1083);
    const auto *osg0_1085 = buffer.data(osg0 + 1085);
    const auto *osg0_1090 = buffer.data(osg0 + 1090);
    const auto *osg0_1092 = buffer.data(osg0 + 1092);
    const auto *osg0_1094 = buffer.data(osg0 + 1094);
    const auto *osg0_1095 = buffer.data(osg0 + 1095);
    const auto *osg0_1098 = buffer.data(osg0 + 1098);
    const auto *osg0_1100 = buffer.data(osg0 + 1100);
    const auto *osg0_1105 = buffer.data(osg0 + 1105);
    const auto *osg0_1107 = buffer.data(osg0 + 1107);
    const auto *osg0_1109 = buffer.data(osg0 + 1109);
    const auto *osg0_1110 = buffer.data(osg0 + 1110);
    const auto *osg0_1113 = buffer.data(osg0 + 1113);
    const auto *osg0_1115 = buffer.data(osg0 + 1115);
    const auto *osg0_1120 = buffer.data(osg0 + 1120);
    const auto *osg0_1122 = buffer.data(osg0 + 1122);
    const auto *osg0_1124 = buffer.data(osg0 + 1124);
    const auto *osg0_1125 = buffer.data(osg0 + 1125);
    const auto *osg0_1128 = buffer.data(osg0 + 1128);
    const auto *osg0_1130 = buffer.data(osg0 + 1130);
    const auto *osg0_1135 = buffer.data(osg0 + 1135);
    const auto *osg0_1137 = buffer.data(osg0 + 1137);
    const auto *osg0_1139 = buffer.data(osg0 + 1139);
    const auto *osg0_1143 = buffer.data(osg0 + 1143);
    const auto *osg0_1150 = buffer.data(osg0 + 1150);
    const auto *osg0_1152 = buffer.data(osg0 + 1152);
    const auto *osg0_1154 = buffer.data(osg0 + 1154);
    const auto *osg0_1155 = buffer.data(osg0 + 1155);
    const auto *osg0_1160 = buffer.data(osg0 + 1160);
    const auto *osg0_1165 = buffer.data(osg0 + 1165);
    const auto *osg0_1166 = buffer.data(osg0 + 1166);
    const auto *osg0_1167 = buffer.data(osg0 + 1167);
    const auto *osg0_1169 = buffer.data(osg0 + 1169);

    const auto *osf_596 = buffer.data(osf + 596);
    const auto *osf_600 = buffer.data(osf + 600);
    const auto *osf_606 = buffer.data(osf + 606);
    const auto *osf_609 = buffer.data(osf + 609);
    const auto *osf_610 = buffer.data(osf + 610);
    const auto *osf_612 = buffer.data(osf + 612);
    const auto *osf_616 = buffer.data(osf + 616);
    const auto *osf_619 = buffer.data(osf + 619);
    const auto *osf_620 = buffer.data(osf + 620);
    const auto *osf_622 = buffer.data(osf + 622);
    const auto *osf_626 = buffer.data(osf + 626);
    const auto *osf_629 = buffer.data(osf + 629);
    const auto *osf_630 = buffer.data(osf + 630);
    const auto *osf_632 = buffer.data(osf + 632);
    const auto *osf_636 = buffer.data(osf + 636);
    const auto *osf_639 = buffer.data(osf + 639);
    const auto *osf_640 = buffer.data(osf + 640);
    const auto *osf_642 = buffer.data(osf + 642);
    const auto *osf_646 = buffer.data(osf + 646);
    const auto *osf_649 = buffer.data(osf + 649);
    const auto *osf_650 = buffer.data(osf + 650);
    const auto *osf_652 = buffer.data(osf + 652);
    const auto *osf_659 = buffer.data(osf + 659);
    const auto *osf_666 = buffer.data(osf + 666);
    const auto *osf_667 = buffer.data(osf + 667);
    const auto *osf_669 = buffer.data(osf + 669);
    const auto *osf_679 = buffer.data(osf + 679);
    const auto *osf_715 = buffer.data(osf + 715);
    const auto *osf_716 = buffer.data(osf + 716);
    const auto *osf_717 = buffer.data(osf + 717);
    const auto *osf_718 = buffer.data(osf + 718);
    const auto *osf_719 = buffer.data(osf + 719);
    const auto *osf_720 = buffer.data(osf + 720);
    const auto *osf_723 = buffer.data(osf + 723);
    const auto *osf_725 = buffer.data(osf + 725);
    const auto *osf_726 = buffer.data(osf + 726);
    const auto *osf_727 = buffer.data(osf + 727);
    const auto *osf_728 = buffer.data(osf + 728);
    const auto *osf_729 = buffer.data(osf + 729);
    const auto *osf_730 = buffer.data(osf + 730);
    const auto *osf_733 = buffer.data(osf + 733);
    const auto *osf_735 = buffer.data(osf + 735);
    const auto *osf_736 = buffer.data(osf + 736);
    const auto *osf_737 = buffer.data(osf + 737);
    const auto *osf_738 = buffer.data(osf + 738);
    const auto *osf_739 = buffer.data(osf + 739);
    const auto *osf_740 = buffer.data(osf + 740);
    const auto *osf_743 = buffer.data(osf + 743);
    const auto *osf_745 = buffer.data(osf + 745);
    const auto *osf_746 = buffer.data(osf + 746);
    const auto *osf_747 = buffer.data(osf + 747);
    const auto *osf_748 = buffer.data(osf + 748);
    const auto *osf_749 = buffer.data(osf + 749);
    const auto *osf_750 = buffer.data(osf + 750);
    const auto *osf_753 = buffer.data(osf + 753);
    const auto *osf_755 = buffer.data(osf + 755);
    const auto *osf_756 = buffer.data(osf + 756);
    const auto *osf_757 = buffer.data(osf + 757);
    const auto *osf_758 = buffer.data(osf + 758);
    const auto *osf_759 = buffer.data(osf + 759);
    const auto *osf_763 = buffer.data(osf + 763);
    const auto *osf_766 = buffer.data(osf + 766);
    const auto *osf_767 = buffer.data(osf + 767);
    const auto *osf_768 = buffer.data(osf + 768);
    const auto *osf_769 = buffer.data(osf + 769);
    const auto *osf_770 = buffer.data(osf + 770);
    const auto *osf_775 = buffer.data(osf + 775);
    const auto *osf_776 = buffer.data(osf + 776);
    const auto *osf_777 = buffer.data(osf + 777);
    const auto *osf_779 = buffer.data(osf + 779);

    const auto *osg1_975 = buffer.data(osg1 + 975);
    const auto *osg1_980 = buffer.data(osg1 + 980);
    const auto *osg1_990 = buffer.data(osg1 + 990);
    const auto *osg1_991 = buffer.data(osg1 + 991);
    const auto *osg1_993 = buffer.data(osg1 + 993);
    const auto *osg1_1000 = buffer.data(osg1 + 1000);
    const auto *osg1_1002 = buffer.data(osg1 + 1002);
    const auto *osg1_1070 = buffer.data(osg1 + 1070);
    const auto *osg1_1075 = buffer.data(osg1 + 1075);
    const auto *osg1_1077 = buffer.data(osg1 + 1077);
    const auto *osg1_1079 = buffer.data(osg1 + 1079);
    const auto *osg1_1080 = buffer.data(osg1 + 1080);
    const auto *osg1_1083 = buffer.data(osg1 + 1083);
    const auto *osg1_1085 = buffer.data(osg1 + 1085);
    const auto *osg1_1090 = buffer.data(osg1 + 1090);
    const auto *osg1_1092 = buffer.data(osg1 + 1092);
    const auto *osg1_1094 = buffer.data(osg1 + 1094);
    const auto *osg1_1095 = buffer.data(osg1 + 1095);
    const auto *osg1_1098 = buffer.data(osg1 + 1098);
    const auto *osg1_1100 = buffer.data(osg1 + 1100);
    const auto *osg1_1105 = buffer.data(osg1 + 1105);
    const auto *osg1_1107 = buffer.data(osg1 + 1107);
    const auto *osg1_1109 = buffer.data(osg1 + 1109);
    const auto *osg1_1110 = buffer.data(osg1 + 1110);
    const auto *osg1_1113 = buffer.data(osg1 + 1113);
    const auto *osg1_1115 = buffer.data(osg1 + 1115);
    const auto *osg1_1120 = buffer.data(osg1 + 1120);
    const auto *osg1_1122 = buffer.data(osg1 + 1122);
    const auto *osg1_1124 = buffer.data(osg1 + 1124);
    const auto *osg1_1125 = buffer.data(osg1 + 1125);
    const auto *osg1_1128 = buffer.data(osg1 + 1128);
    const auto *osg1_1130 = buffer.data(osg1 + 1130);
    const auto *osg1_1135 = buffer.data(osg1 + 1135);
    const auto *osg1_1137 = buffer.data(osg1 + 1137);
    const auto *osg1_1139 = buffer.data(osg1 + 1139);
    const auto *osg1_1143 = buffer.data(osg1 + 1143);
    const auto *osg1_1150 = buffer.data(osg1 + 1150);
    const auto *osg1_1152 = buffer.data(osg1 + 1152);
    const auto *osg1_1154 = buffer.data(osg1 + 1154);
    const auto *osg1_1155 = buffer.data(osg1 + 1155);
    const auto *osg1_1160 = buffer.data(osg1 + 1160);
    const auto *osg1_1165 = buffer.data(osg1 + 1165);
    const auto *osg1_1166 = buffer.data(osg1 + 1166);
    const auto *osg1_1167 = buffer.data(osg1 + 1167);
    const auto *osg1_1169 = buffer.data(osg1 + 1169);

    const auto *qsd0_462 = buffer.data(qsd0 + 462);
    const auto *qsd0_468 = buffer.data(qsd0 + 468);
    const auto *qsd0_469 = buffer.data(qsd0 + 469);
    const auto *qsd0_471 = buffer.data(qsd0 + 471);
    const auto *qsd0_473 = buffer.data(qsd0 + 473);
    const auto *qsd0_476 = buffer.data(qsd0 + 476);
    const auto *qsd0_478 = buffer.data(qsd0 + 478);
    const auto *qsd0_479 = buffer.data(qsd0 + 479);

    const auto *qsd1_462 = buffer.data(qsd1 + 462);
    const auto *qsd1_468 = buffer.data(qsd1 + 468);
    const auto *qsd1_469 = buffer.data(qsd1 + 469);
    const auto *qsd1_471 = buffer.data(qsd1 + 471);
    const auto *qsd1_473 = buffer.data(qsd1 + 473);
    const auto *qsd1_476 = buffer.data(qsd1 + 476);
    const auto *qsd1_478 = buffer.data(qsd1 + 478);
    const auto *qsd1_479 = buffer.data(qsd1 + 479);

    const auto *qsf_716 = buffer.data(qsf + 716);
    const auto *qsf_717 = buffer.data(qsf + 717);
    const auto *qsf_718 = buffer.data(qsf + 718);
    const auto *qsf_719 = buffer.data(qsf + 719);
    const auto *qsf_720 = buffer.data(qsf + 720);
    const auto *qsf_722 = buffer.data(qsf + 722);
    const auto *qsf_726 = buffer.data(qsf + 726);
    const auto *qsf_727 = buffer.data(qsf + 727);
    const auto *qsf_728 = buffer.data(qsf + 728);
    const auto *qsf_729 = buffer.data(qsf + 729);
    const auto *qsf_730 = buffer.data(qsf + 730);
    const auto *qsf_732 = buffer.data(qsf + 732);
    const auto *qsf_736 = buffer.data(qsf + 736);
    const auto *qsf_737 = buffer.data(qsf + 737);
    const auto *qsf_738 = buffer.data(qsf + 738);
    const auto *qsf_739 = buffer.data(qsf + 739);
    const auto *qsf_740 = buffer.data(qsf + 740);
    const auto *qsf_742 = buffer.data(qsf + 742);
    const auto *qsf_746 = buffer.data(qsf + 746);
    const auto *qsf_747 = buffer.data(qsf + 747);
    const auto *qsf_748 = buffer.data(qsf + 748);
    const auto *qsf_749 = buffer.data(qsf + 749);
    const auto *qsf_750 = buffer.data(qsf + 750);
    const auto *qsf_752 = buffer.data(qsf + 752);
    const auto *qsf_756 = buffer.data(qsf + 756);
    const auto *qsf_757 = buffer.data(qsf + 757);
    const auto *qsf_758 = buffer.data(qsf + 758);
    const auto *qsf_759 = buffer.data(qsf + 759);
    const auto *qsf_760 = buffer.data(qsf + 760);
    const auto *qsf_762 = buffer.data(qsf + 762);
    const auto *qsf_766 = buffer.data(qsf + 766);
    const auto *qsf_767 = buffer.data(qsf + 767);
    const auto *qsf_768 = buffer.data(qsf + 768);
    const auto *qsf_769 = buffer.data(qsf + 769);
    const auto *qsf_770 = buffer.data(qsf + 770);
    const auto *qsf_771 = buffer.data(qsf + 771);
    const auto *qsf_772 = buffer.data(qsf + 772);
    const auto *qsf_775 = buffer.data(qsf + 775);
    const auto *qsf_776 = buffer.data(qsf + 776);
    const auto *qsf_777 = buffer.data(qsf + 777);
    const auto *qsf_779 = buffer.data(qsf + 779);
    const auto *qsf_780 = buffer.data(qsf + 780);
    const auto *qsf_781 = buffer.data(qsf + 781);
    const auto *qsf_783 = buffer.data(qsf + 783);
    const auto *qsf_785 = buffer.data(qsf + 785);
    const auto *qsf_786 = buffer.data(qsf + 786);
    const auto *qsf_787 = buffer.data(qsf + 787);
    const auto *qsf_788 = buffer.data(qsf + 788);
    const auto *qsf_789 = buffer.data(qsf + 789);
    const auto *qsf_792 = buffer.data(qsf + 792);
    const auto *qsf_794 = buffer.data(qsf + 794);
    const auto *qsf_795 = buffer.data(qsf + 795);
    const auto *qsf_796 = buffer.data(qsf + 796);
    const auto *qsf_797 = buffer.data(qsf + 797);
    const auto *qsf_798 = buffer.data(qsf + 798);
    const auto *qsf_799 = buffer.data(qsf + 799);

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, pa_x, pc_x, osg0_1070, osf_715, \
                         osf_716, osf_717, osf_718, osg1_1070, qsf_716, qsf_717, \
                         qsf_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = pa_x[k] * osg0_1070[k]
                    + f_8 * osf_715[k]
                    - f_6 * pc_x[k] * osg1_1070[k];

        t_1071[k] = f_7 * osf_716[k]
                    + f_3 * pc_x[k] * qsf_716[k];

        t_1072[k] = f_7 * osf_717[k]
                    + f_3 * pc_x[k] * qsf_717[k];

        t_1073[k] = f_7 * osf_718[k]
                    + f_3 * pc_x[k] * qsf_718[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, t_1077, pa_x, pc_x, pc_z, osg0_1075, \
                         osg0_1077, osf_596, osf_719, osg1_1075, osg1_1077, qsf_716, \
                         qsf_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = f_7 * osf_719[k]
                    + f_3 * pc_x[k] * qsf_719[k];

        t_1075[k] = pa_x[k] * osg0_1075[k]
                    - f_6 * pc_x[k] * osg1_1075[k];

        t_1076[k] = f_18 * osf_596[k]
                    + f_3 * pc_z[k] * qsf_716[k];

        t_1077[k] = pa_x[k] * osg0_1077[k]
                    - f_6 * pc_x[k] * osg1_1077[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, t_1081, pa_x, pc_x, pc_y, osg0_1079, \
                         osg0_1080, osf_609, osf_610, osf_720, osg1_1079, osg1_1080, qsf_719, \
                         qsf_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_19 * osf_609[k]
                    + f_3 * pc_y[k] * qsf_719[k];

        t_1079[k] = pa_x[k] * osg0_1079[k]
                    - f_6 * pc_x[k] * osg1_1079[k];

        t_1080[k] = pa_x[k] * osg0_1080[k]
                    + f_16 * osf_720[k]
                    - f_6 * pc_x[k] * osg1_1080[k];

        t_1081[k] = f_18 * osf_610[k]
                    + f_3 * pc_y[k] * qsf_720[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pa_x, pc_x, pc_y, pc_z, osg0_1083, osf_600, \
                         osf_612, osf_723, osg1_1083, qsf_720, \
                         qsf_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_19 * osf_600[k]
                    + f_3 * pc_z[k] * qsf_720[k];

        t_1083[k] = pa_x[k] * osg0_1083[k]
                    + f_8 * osf_723[k]
                    - f_6 * pc_x[k] * osg1_1083[k];

        t_1084[k] = f_18 * osf_612[k]
                    + f_3 * pc_y[k] * qsf_722[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, t_1088, pa_x, pc_x, osg0_1085, osf_725, \
                         osf_726, osf_727, osf_728, osg1_1085, qsf_726, qsf_727, \
                         qsf_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = pa_x[k] * osg0_1085[k]
                    + f_8 * osf_725[k]
                    - f_6 * pc_x[k] * osg1_1085[k];

        t_1086[k] = f_7 * osf_726[k]
                    + f_3 * pc_x[k] * qsf_726[k];

        t_1087[k] = f_7 * osf_727[k]
                    + f_3 * pc_x[k] * qsf_727[k];

        t_1088[k] = f_7 * osf_728[k]
                    + f_3 * pc_x[k] * qsf_728[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, t_1092, pa_x, pc_x, pc_z, osg0_1090, \
                         osg0_1092, osf_606, osf_729, osg1_1090, osg1_1092, qsf_726, \
                         qsf_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_7 * osf_729[k]
                    + f_3 * pc_x[k] * qsf_729[k];

        t_1090[k] = pa_x[k] * osg0_1090[k]
                    - f_6 * pc_x[k] * osg1_1090[k];

        t_1091[k] = f_19 * osf_606[k]
                    + f_3 * pc_z[k] * qsf_726[k];

        t_1092[k] = pa_x[k] * osg0_1092[k]
                    - f_6 * pc_x[k] * osg1_1092[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pa_x, pc_x, pc_y, osg0_1094, \
                         osg0_1095, osf_619, osf_620, osf_730, osg1_1094, osg1_1095, qsf_729, \
                         qsf_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_18 * osf_619[k]
                    + f_3 * pc_y[k] * qsf_729[k];

        t_1094[k] = pa_x[k] * osg0_1094[k]
                    - f_6 * pc_x[k] * osg1_1094[k];

        t_1095[k] = pa_x[k] * osg0_1095[k]
                    + f_16 * osf_730[k]
                    - f_6 * pc_x[k] * osg1_1095[k];

        t_1096[k] = f_16 * osf_620[k]
                    + f_3 * pc_y[k] * qsf_730[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pa_x, pc_x, pc_y, pc_z, osg0_1098, osf_610, \
                         osf_622, osf_733, osg1_1098, qsf_730, \
                         qsf_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_17 * osf_610[k]
                    + f_3 * pc_z[k] * qsf_730[k];

        t_1098[k] = pa_x[k] * osg0_1098[k]
                    + f_8 * osf_733[k]
                    - f_6 * pc_x[k] * osg1_1098[k];

        t_1099[k] = f_16 * osf_622[k]
                    + f_3 * pc_y[k] * qsf_732[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pa_x, pc_x, osg0_1100, osf_735, \
                         osf_736, osf_737, osf_738, osg1_1100, qsf_736, qsf_737, \
                         qsf_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = pa_x[k] * osg0_1100[k]
                    + f_8 * osf_735[k]
                    - f_6 * pc_x[k] * osg1_1100[k];

        t_1101[k] = f_7 * osf_736[k]
                    + f_3 * pc_x[k] * qsf_736[k];

        t_1102[k] = f_7 * osf_737[k]
                    + f_3 * pc_x[k] * qsf_737[k];

        t_1103[k] = f_7 * osf_738[k]
                    + f_3 * pc_x[k] * qsf_738[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pa_x, pc_x, pc_z, osg0_1105, \
                         osg0_1107, osf_616, osf_739, osg1_1105, osg1_1107, qsf_736, \
                         qsf_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_7 * osf_739[k]
                    + f_3 * pc_x[k] * qsf_739[k];

        t_1105[k] = pa_x[k] * osg0_1105[k]
                    - f_6 * pc_x[k] * osg1_1105[k];

        t_1106[k] = f_17 * osf_616[k]
                    + f_3 * pc_z[k] * qsf_736[k];

        t_1107[k] = pa_x[k] * osg0_1107[k]
                    - f_6 * pc_x[k] * osg1_1107[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, t_1111, pa_x, pc_x, pc_y, osg0_1109, \
                         osg0_1110, osf_629, osf_630, osf_740, osg1_1109, osg1_1110, qsf_739, \
                         qsf_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_16 * osf_629[k]
                    + f_3 * pc_y[k] * qsf_739[k];

        t_1109[k] = pa_x[k] * osg0_1109[k]
                    - f_6 * pc_x[k] * osg1_1109[k];

        t_1110[k] = pa_x[k] * osg0_1110[k]
                    + f_16 * osf_740[k]
                    - f_6 * pc_x[k] * osg1_1110[k];

        t_1111[k] = f_14 * osf_630[k]
                    + f_3 * pc_y[k] * qsf_740[k];
    }

#pragma omp simd aligned(t_1112, t_1113, t_1114, pa_x, pc_x, pc_y, pc_z, osg0_1113, osf_620, \
                         osf_632, osf_743, osg1_1113, qsf_740, \
                         qsf_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1112[k] = f_15 * osf_620[k]
                    + f_3 * pc_z[k] * qsf_740[k];

        t_1113[k] = pa_x[k] * osg0_1113[k]
                    + f_8 * osf_743[k]
                    - f_6 * pc_x[k] * osg1_1113[k];

        t_1114[k] = f_14 * osf_632[k]
                    + f_3 * pc_y[k] * qsf_742[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, pa_x, pc_x, osg0_1115, osf_745, \
                         osf_746, osf_747, osf_748, osg1_1115, qsf_746, qsf_747, \
                         qsf_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = pa_x[k] * osg0_1115[k]
                    + f_8 * osf_745[k]
                    - f_6 * pc_x[k] * osg1_1115[k];

        t_1116[k] = f_7 * osf_746[k]
                    + f_3 * pc_x[k] * qsf_746[k];

        t_1117[k] = f_7 * osf_747[k]
                    + f_3 * pc_x[k] * qsf_747[k];

        t_1118[k] = f_7 * osf_748[k]
                    + f_3 * pc_x[k] * qsf_748[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, t_1122, pa_x, pc_x, pc_z, osg0_1120, \
                         osg0_1122, osf_626, osf_749, osg1_1120, osg1_1122, qsf_746, \
                         qsf_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = f_7 * osf_749[k]
                    + f_3 * pc_x[k] * qsf_749[k];

        t_1120[k] = pa_x[k] * osg0_1120[k]
                    - f_6 * pc_x[k] * osg1_1120[k];

        t_1121[k] = f_15 * osf_626[k]
                    + f_3 * pc_z[k] * qsf_746[k];

        t_1122[k] = pa_x[k] * osg0_1122[k]
                    - f_6 * pc_x[k] * osg1_1122[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, pa_x, pc_x, pc_y, osg0_1124, \
                         osg0_1125, osf_639, osf_640, osf_750, osg1_1124, osg1_1125, qsf_749, \
                         qsf_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_14 * osf_639[k]
                    + f_3 * pc_y[k] * qsf_749[k];

        t_1124[k] = pa_x[k] * osg0_1124[k]
                    - f_6 * pc_x[k] * osg1_1124[k];

        t_1125[k] = pa_x[k] * osg0_1125[k]
                    + f_16 * osf_750[k]
                    - f_6 * pc_x[k] * osg1_1125[k];

        t_1126[k] = f_8 * osf_640[k]
                    + f_3 * pc_y[k] * qsf_750[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pa_x, pc_x, pc_y, pc_z, osg0_1128, osf_630, \
                         osf_642, osf_753, osg1_1128, qsf_750, \
                         qsf_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_13 * osf_630[k]
                    + f_3 * pc_z[k] * qsf_750[k];

        t_1128[k] = pa_x[k] * osg0_1128[k]
                    + f_8 * osf_753[k]
                    - f_6 * pc_x[k] * osg1_1128[k];

        t_1129[k] = f_8 * osf_642[k]
                    + f_3 * pc_y[k] * qsf_752[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, pa_x, pc_x, osg0_1130, osf_755, \
                         osf_756, osf_757, osf_758, osg1_1130, qsf_756, qsf_757, \
                         qsf_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pa_x[k] * osg0_1130[k]
                    + f_8 * osf_755[k]
                    - f_6 * pc_x[k] * osg1_1130[k];

        t_1131[k] = f_7 * osf_756[k]
                    + f_3 * pc_x[k] * qsf_756[k];

        t_1132[k] = f_7 * osf_757[k]
                    + f_3 * pc_x[k] * qsf_757[k];

        t_1133[k] = f_7 * osf_758[k]
                    + f_3 * pc_x[k] * qsf_758[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pa_x, pc_x, pc_z, osg0_1135, \
                         osg0_1137, osf_636, osf_759, osg1_1135, osg1_1137, qsf_756, \
                         qsf_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_7 * osf_759[k]
                    + f_3 * pc_x[k] * qsf_759[k];

        t_1135[k] = pa_x[k] * osg0_1135[k]
                    - f_6 * pc_x[k] * osg1_1135[k];

        t_1136[k] = f_13 * osf_636[k]
                    + f_3 * pc_z[k] * qsf_756[k];

        t_1137[k] = pa_x[k] * osg0_1137[k]
                    - f_6 * pc_x[k] * osg1_1137[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, pa_x, pa_y, pc_x, pc_y, osg0_975, \
                         osg0_1139, osf_649, osf_650, osg1_975, osg1_1139, qsf_759, \
                         qsf_760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_8 * osf_649[k]
                    + f_3 * pc_y[k] * qsf_759[k];

        t_1139[k] = pa_x[k] * osg0_1139[k]
                    - f_6 * pc_x[k] * osg1_1139[k];

        t_1140[k] = pa_y[k] * osg0_975[k]
                    - f_6 * pc_y[k] * osg1_975[k];

        t_1141[k] = f_7 * osf_650[k]
                    + f_3 * pc_y[k] * qsf_760[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pa_x, pc_x, pc_y, pc_z, osg0_1143, osf_640, \
                         osf_652, osf_763, osg1_1143, qsf_760, \
                         qsf_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_12 * osf_640[k]
                    + f_3 * pc_z[k] * qsf_760[k];

        t_1143[k] = pa_x[k] * osg0_1143[k]
                    + f_8 * osf_763[k]
                    - f_6 * pc_x[k] * osg1_1143[k];

        t_1144[k] = f_7 * osf_652[k]
                    + f_3 * pc_y[k] * qsf_762[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, pa_y, pc_x, pc_y, osg0_980, osf_766, \
                         osf_767, osf_768, osg1_980, qsf_766, qsf_767, \
                         qsf_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = pa_y[k] * osg0_980[k]
                    - f_6 * pc_y[k] * osg1_980[k];

        t_1146[k] = f_7 * osf_766[k]
                    + f_3 * pc_x[k] * qsf_766[k];

        t_1147[k] = f_7 * osf_767[k]
                    + f_3 * pc_x[k] * qsf_767[k];

        t_1148[k] = f_7 * osf_768[k]
                    + f_3 * pc_x[k] * qsf_768[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pa_x, pc_x, pc_z, osg0_1150, \
                         osg0_1152, osf_646, osf_769, osg1_1150, osg1_1152, qsf_766, \
                         qsf_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_7 * osf_769[k]
                    + f_3 * pc_x[k] * qsf_769[k];

        t_1150[k] = pa_x[k] * osg0_1150[k]
                    - f_6 * pc_x[k] * osg1_1150[k];

        t_1151[k] = f_12 * osf_646[k]
                    + f_3 * pc_z[k] * qsf_766[k];

        t_1152[k] = pa_x[k] * osg0_1152[k]
                    - f_6 * pc_x[k] * osg1_1152[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pa_x, pc_x, pc_y, osg0_1154, \
                         osg0_1155, osf_659, osf_770, osg1_1154, osg1_1155, qsf_769, \
                         qsf_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_7 * osf_659[k]
                    + f_3 * pc_y[k] * qsf_769[k];

        t_1154[k] = pa_x[k] * osg0_1154[k]
                    - f_6 * pc_x[k] * osg1_1154[k];

        t_1155[k] = pa_x[k] * osg0_1155[k]
                    + f_16 * osf_770[k]
                    - f_6 * pc_x[k] * osg1_1155[k];

        t_1156[k] = f_3 * pc_y[k] * qsf_770[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_y, pc_z, osf_650, qsd0_462, qsd1_462, \
                         qsf_770, qsf_771, qsf_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_9 * osf_650[k]
                    + f_3 * pc_z[k] * qsf_770[k];

        t_1158[k] = f_4 * qsd0_462[k]
                    - f_5 * qsd1_462[k]
                    + f_3 * pc_y[k] * qsf_771[k];

        t_1159[k] = f_3 * pc_y[k] * qsf_772[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, pa_x, pc_x, pc_y, osg0_1160, osf_775, \
                         osf_776, osf_777, osg1_1160, qsf_775, qsf_776, \
                         qsf_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = pa_x[k] * osg0_1160[k]
                    + f_8 * osf_775[k]
                    - f_6 * pc_x[k] * osg1_1160[k];

        t_1161[k] = f_7 * osf_776[k]
                    + f_3 * pc_x[k] * qsf_776[k];

        t_1162[k] = f_7 * osf_777[k]
                    + f_3 * pc_x[k] * qsf_777[k];

        t_1163[k] = f_3 * pc_y[k] * qsf_775[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, t_1167, t_1168, pa_x, pc_x, pc_y, osg0_1165, \
                         osg0_1166, osg0_1167, osf_779, osg1_1165, osg1_1166, osg1_1167, \
                         qsf_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_7 * osf_779[k]
                    + f_3 * pc_x[k] * qsf_779[k];

        t_1165[k] = pa_x[k] * osg0_1165[k]
                    - f_6 * pc_x[k] * osg1_1165[k];

        t_1166[k] = pa_x[k] * osg0_1166[k]
                    - f_6 * pc_x[k] * osg1_1166[k];

        t_1167[k] = pa_x[k] * osg0_1167[k]
                    - f_6 * pc_x[k] * osg1_1167[k];

        t_1168[k] = f_3 * pc_y[k] * qsf_779[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pa_x, pc_x, pc_z, osg0_1169, \
                         osg1_1169, qsd0_468, qsd0_469, qsd1_468, qsd1_469, qsf_780, \
                         qsf_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = pa_x[k] * osg0_1169[k]
                    - f_6 * pc_x[k] * osg1_1169[k];

        t_1170[k] = f_1 * qsd0_468[k]
                    - f_2 * qsd1_468[k]
                    + f_3 * pc_x[k] * qsf_780[k];

        t_1171[k] = f_10 * qsd0_469[k]
                    - f_11 * qsd1_469[k]
                    + f_3 * pc_x[k] * qsf_781[k];

        t_1172[k] = f_3 * pc_z[k] * qsf_780[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, t_1177, pc_x, pc_z, qsd0_471, \
                         qsd0_473, qsd1_471, qsd1_473, qsf_781, qsf_783, qsf_785, qsf_786, \
                         qsf_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_4 * qsd0_471[k]
                    - f_5 * qsd1_471[k]
                    + f_3 * pc_x[k] * qsf_783[k];

        t_1174[k] = f_3 * pc_z[k] * qsf_781[k];

        t_1175[k] = f_4 * qsd0_473[k]
                    - f_5 * qsd1_473[k]
                    + f_3 * pc_x[k] * qsf_785[k];

        t_1176[k] = f_3 * pc_x[k] * qsf_786[k];

        t_1177[k] = f_3 * pc_x[k] * qsf_787[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, t_1181, t_1182, pc_x, pc_y, pc_z, osf_666, \
                         qsd0_471, qsd1_471, qsf_786, qsf_787, qsf_788, \
                         qsf_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_3 * pc_x[k] * qsf_788[k];

        t_1179[k] = f_3 * pc_x[k] * qsf_789[k];

        t_1180[k] = f_0 * osf_666[k]
                    + f_1 * qsd0_471[k]
                    - f_2 * qsd1_471[k]
                    + f_3 * pc_y[k] * qsf_786[k];

        t_1181[k] = f_3 * pc_z[k] * qsf_786[k];

        t_1182[k] = f_4 * qsd0_471[k]
                    - f_5 * qsd1_471[k]
                    + f_3 * pc_z[k] * qsf_787[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, t_1186, pa_z, pc_y, pc_z, osg0_990, osg0_991, \
                         osf_669, osg1_990, osg1_991, qsd0_473, qsd1_473, \
                         qsf_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_0 * osf_669[k]
                    + f_3 * pc_y[k] * qsf_789[k];

        t_1184[k] = f_1 * qsd0_473[k]
                    - f_2 * qsd1_473[k]
                    + f_3 * pc_z[k] * qsf_789[k];

        t_1185[k] = pa_z[k] * osg0_990[k]
                    - f_6 * pc_z[k] * osg1_990[k];

        t_1186[k] = pa_z[k] * osg0_991[k]
                    - f_6 * pc_z[k] * osg1_991[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, pa_z, pc_x, pc_z, osg0_993, osg1_993, \
                         qsd0_476, qsd0_478, qsd1_476, qsd1_478, qsf_792, \
                         qsf_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = f_10 * qsd0_476[k]
                    - f_11 * qsd1_476[k]
                    + f_3 * pc_x[k] * qsf_792[k];

        t_1188[k] = pa_z[k] * osg0_993[k]
                    - f_6 * pc_z[k] * osg1_993[k];

        t_1189[k] = f_4 * qsd0_478[k]
                    - f_5 * qsd1_478[k]
                    + f_3 * pc_x[k] * qsf_794[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, t_1193, t_1194, pc_x, qsd0_479, qsd1_479, \
                         qsf_795, qsf_796, qsf_797, qsf_798, qsf_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = f_4 * qsd0_479[k]
                    - f_5 * qsd1_479[k]
                    + f_3 * pc_x[k] * qsf_795[k];

        t_1191[k] = f_3 * pc_x[k] * qsf_796[k];

        t_1192[k] = f_3 * pc_x[k] * qsf_797[k];

        t_1193[k] = f_3 * pc_x[k] * qsf_798[k];

        t_1194[k] = f_3 * pc_x[k] * qsf_799[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pa_z, pc_y, pc_z, osg0_1000, \
                         osg0_1002, osf_666, osf_667, osf_679, osg1_1000, osg1_1002, qsf_796, \
                         qsf_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = pa_z[k] * osg0_1000[k]
                    - f_6 * pc_z[k] * osg1_1000[k];

        t_1196[k] = f_7 * osf_666[k]
                    + f_3 * pc_z[k] * qsf_796[k];

        t_1197[k] = pa_z[k] * osg0_1002[k]
                    + f_8 * osf_667[k]
                    - f_6 * pc_z[k] * osg1_1002[k];

        t_1198[k] = f_9 * osf_679[k]
                    + f_3 * pc_y[k] * qsf_799[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t osf, const size_t qsd0,
                                                           const size_t qsd1, const size_t qsf,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf_669 = buffer.data(osf + 669);
    const auto *osf_676 = buffer.data(osf + 676);
    const auto *osf_679 = buffer.data(osf + 679);
    const auto *osf_686 = buffer.data(osf + 686);
    const auto *osf_688 = buffer.data(osf + 688);
    const auto *osf_689 = buffer.data(osf + 689);
    const auto *osf_696 = buffer.data(osf + 696);
    const auto *osf_698 = buffer.data(osf + 698);
    const auto *osf_699 = buffer.data(osf + 699);
    const auto *osf_706 = buffer.data(osf + 706);
    const auto *osf_708 = buffer.data(osf + 708);
    const auto *osf_709 = buffer.data(osf + 709);
    const auto *osf_716 = buffer.data(osf + 716);
    const auto *osf_718 = buffer.data(osf + 718);
    const auto *osf_719 = buffer.data(osf + 719);
    const auto *osf_726 = buffer.data(osf + 726);
    const auto *osf_728 = buffer.data(osf + 728);
    const auto *osf_729 = buffer.data(osf + 729);
    const auto *osf_736 = buffer.data(osf + 736);
    const auto *osf_738 = buffer.data(osf + 738);
    const auto *osf_739 = buffer.data(osf + 739);
    const auto *osf_746 = buffer.data(osf + 746);
    const auto *osf_748 = buffer.data(osf + 748);
    const auto *osf_749 = buffer.data(osf + 749);
    const auto *osf_756 = buffer.data(osf + 756);
    const auto *osf_758 = buffer.data(osf + 758);
    const auto *osf_759 = buffer.data(osf + 759);

    const auto *qsd0_479 = buffer.data(qsd0 + 479);
    const auto *qsd0_480 = buffer.data(qsd0 + 480);
    const auto *qsd0_481 = buffer.data(qsd0 + 481);
    const auto *qsd0_482 = buffer.data(qsd0 + 482);
    const auto *qsd0_483 = buffer.data(qsd0 + 483);
    const auto *qsd0_484 = buffer.data(qsd0 + 484);
    const auto *qsd0_485 = buffer.data(qsd0 + 485);
    const auto *qsd0_486 = buffer.data(qsd0 + 486);
    const auto *qsd0_487 = buffer.data(qsd0 + 487);
    const auto *qsd0_488 = buffer.data(qsd0 + 488);
    const auto *qsd0_489 = buffer.data(qsd0 + 489);
    const auto *qsd0_490 = buffer.data(qsd0 + 490);
    const auto *qsd0_491 = buffer.data(qsd0 + 491);
    const auto *qsd0_492 = buffer.data(qsd0 + 492);
    const auto *qsd0_493 = buffer.data(qsd0 + 493);
    const auto *qsd0_494 = buffer.data(qsd0 + 494);
    const auto *qsd0_495 = buffer.data(qsd0 + 495);
    const auto *qsd0_496 = buffer.data(qsd0 + 496);
    const auto *qsd0_497 = buffer.data(qsd0 + 497);
    const auto *qsd0_498 = buffer.data(qsd0 + 498);
    const auto *qsd0_499 = buffer.data(qsd0 + 499);
    const auto *qsd0_500 = buffer.data(qsd0 + 500);
    const auto *qsd0_501 = buffer.data(qsd0 + 501);
    const auto *qsd0_502 = buffer.data(qsd0 + 502);
    const auto *qsd0_503 = buffer.data(qsd0 + 503);
    const auto *qsd0_504 = buffer.data(qsd0 + 504);
    const auto *qsd0_505 = buffer.data(qsd0 + 505);
    const auto *qsd0_506 = buffer.data(qsd0 + 506);
    const auto *qsd0_507 = buffer.data(qsd0 + 507);
    const auto *qsd0_508 = buffer.data(qsd0 + 508);
    const auto *qsd0_509 = buffer.data(qsd0 + 509);
    const auto *qsd0_510 = buffer.data(qsd0 + 510);
    const auto *qsd0_511 = buffer.data(qsd0 + 511);
    const auto *qsd0_512 = buffer.data(qsd0 + 512);
    const auto *qsd0_513 = buffer.data(qsd0 + 513);
    const auto *qsd0_514 = buffer.data(qsd0 + 514);
    const auto *qsd0_515 = buffer.data(qsd0 + 515);
    const auto *qsd0_516 = buffer.data(qsd0 + 516);
    const auto *qsd0_517 = buffer.data(qsd0 + 517);
    const auto *qsd0_518 = buffer.data(qsd0 + 518);
    const auto *qsd0_519 = buffer.data(qsd0 + 519);
    const auto *qsd0_520 = buffer.data(qsd0 + 520);
    const auto *qsd0_521 = buffer.data(qsd0 + 521);
    const auto *qsd0_522 = buffer.data(qsd0 + 522);
    const auto *qsd0_523 = buffer.data(qsd0 + 523);
    const auto *qsd0_524 = buffer.data(qsd0 + 524);
    const auto *qsd0_525 = buffer.data(qsd0 + 525);
    const auto *qsd0_526 = buffer.data(qsd0 + 526);
    const auto *qsd0_527 = buffer.data(qsd0 + 527);

    const auto *qsd1_479 = buffer.data(qsd1 + 479);
    const auto *qsd1_480 = buffer.data(qsd1 + 480);
    const auto *qsd1_481 = buffer.data(qsd1 + 481);
    const auto *qsd1_482 = buffer.data(qsd1 + 482);
    const auto *qsd1_483 = buffer.data(qsd1 + 483);
    const auto *qsd1_484 = buffer.data(qsd1 + 484);
    const auto *qsd1_485 = buffer.data(qsd1 + 485);
    const auto *qsd1_486 = buffer.data(qsd1 + 486);
    const auto *qsd1_487 = buffer.data(qsd1 + 487);
    const auto *qsd1_488 = buffer.data(qsd1 + 488);
    const auto *qsd1_489 = buffer.data(qsd1 + 489);
    const auto *qsd1_490 = buffer.data(qsd1 + 490);
    const auto *qsd1_491 = buffer.data(qsd1 + 491);
    const auto *qsd1_492 = buffer.data(qsd1 + 492);
    const auto *qsd1_493 = buffer.data(qsd1 + 493);
    const auto *qsd1_494 = buffer.data(qsd1 + 494);
    const auto *qsd1_495 = buffer.data(qsd1 + 495);
    const auto *qsd1_496 = buffer.data(qsd1 + 496);
    const auto *qsd1_497 = buffer.data(qsd1 + 497);
    const auto *qsd1_498 = buffer.data(qsd1 + 498);
    const auto *qsd1_499 = buffer.data(qsd1 + 499);
    const auto *qsd1_500 = buffer.data(qsd1 + 500);
    const auto *qsd1_501 = buffer.data(qsd1 + 501);
    const auto *qsd1_502 = buffer.data(qsd1 + 502);
    const auto *qsd1_503 = buffer.data(qsd1 + 503);
    const auto *qsd1_504 = buffer.data(qsd1 + 504);
    const auto *qsd1_505 = buffer.data(qsd1 + 505);
    const auto *qsd1_506 = buffer.data(qsd1 + 506);
    const auto *qsd1_507 = buffer.data(qsd1 + 507);
    const auto *qsd1_508 = buffer.data(qsd1 + 508);
    const auto *qsd1_509 = buffer.data(qsd1 + 509);
    const auto *qsd1_510 = buffer.data(qsd1 + 510);
    const auto *qsd1_511 = buffer.data(qsd1 + 511);
    const auto *qsd1_512 = buffer.data(qsd1 + 512);
    const auto *qsd1_513 = buffer.data(qsd1 + 513);
    const auto *qsd1_514 = buffer.data(qsd1 + 514);
    const auto *qsd1_515 = buffer.data(qsd1 + 515);
    const auto *qsd1_516 = buffer.data(qsd1 + 516);
    const auto *qsd1_517 = buffer.data(qsd1 + 517);
    const auto *qsd1_518 = buffer.data(qsd1 + 518);
    const auto *qsd1_519 = buffer.data(qsd1 + 519);
    const auto *qsd1_520 = buffer.data(qsd1 + 520);
    const auto *qsd1_521 = buffer.data(qsd1 + 521);
    const auto *qsd1_522 = buffer.data(qsd1 + 522);
    const auto *qsd1_523 = buffer.data(qsd1 + 523);
    const auto *qsd1_524 = buffer.data(qsd1 + 524);
    const auto *qsd1_525 = buffer.data(qsd1 + 525);
    const auto *qsd1_526 = buffer.data(qsd1 + 526);
    const auto *qsd1_527 = buffer.data(qsd1 + 527);

    const auto *qsf_799 = buffer.data(qsf + 799);
    const auto *qsf_800 = buffer.data(qsf + 800);
    const auto *qsf_801 = buffer.data(qsf + 801);
    const auto *qsf_802 = buffer.data(qsf + 802);
    const auto *qsf_803 = buffer.data(qsf + 803);
    const auto *qsf_804 = buffer.data(qsf + 804);
    const auto *qsf_805 = buffer.data(qsf + 805);
    const auto *qsf_806 = buffer.data(qsf + 806);
    const auto *qsf_807 = buffer.data(qsf + 807);
    const auto *qsf_808 = buffer.data(qsf + 808);
    const auto *qsf_809 = buffer.data(qsf + 809);
    const auto *qsf_810 = buffer.data(qsf + 810);
    const auto *qsf_811 = buffer.data(qsf + 811);
    const auto *qsf_812 = buffer.data(qsf + 812);
    const auto *qsf_813 = buffer.data(qsf + 813);
    const auto *qsf_814 = buffer.data(qsf + 814);
    const auto *qsf_815 = buffer.data(qsf + 815);
    const auto *qsf_816 = buffer.data(qsf + 816);
    const auto *qsf_817 = buffer.data(qsf + 817);
    const auto *qsf_818 = buffer.data(qsf + 818);
    const auto *qsf_819 = buffer.data(qsf + 819);
    const auto *qsf_820 = buffer.data(qsf + 820);
    const auto *qsf_821 = buffer.data(qsf + 821);
    const auto *qsf_822 = buffer.data(qsf + 822);
    const auto *qsf_823 = buffer.data(qsf + 823);
    const auto *qsf_824 = buffer.data(qsf + 824);
    const auto *qsf_825 = buffer.data(qsf + 825);
    const auto *qsf_826 = buffer.data(qsf + 826);
    const auto *qsf_827 = buffer.data(qsf + 827);
    const auto *qsf_828 = buffer.data(qsf + 828);
    const auto *qsf_829 = buffer.data(qsf + 829);
    const auto *qsf_830 = buffer.data(qsf + 830);
    const auto *qsf_831 = buffer.data(qsf + 831);
    const auto *qsf_832 = buffer.data(qsf + 832);
    const auto *qsf_833 = buffer.data(qsf + 833);
    const auto *qsf_834 = buffer.data(qsf + 834);
    const auto *qsf_835 = buffer.data(qsf + 835);
    const auto *qsf_836 = buffer.data(qsf + 836);
    const auto *qsf_837 = buffer.data(qsf + 837);
    const auto *qsf_838 = buffer.data(qsf + 838);
    const auto *qsf_839 = buffer.data(qsf + 839);
    const auto *qsf_840 = buffer.data(qsf + 840);
    const auto *qsf_841 = buffer.data(qsf + 841);
    const auto *qsf_842 = buffer.data(qsf + 842);
    const auto *qsf_843 = buffer.data(qsf + 843);
    const auto *qsf_844 = buffer.data(qsf + 844);
    const auto *qsf_845 = buffer.data(qsf + 845);
    const auto *qsf_846 = buffer.data(qsf + 846);
    const auto *qsf_847 = buffer.data(qsf + 847);
    const auto *qsf_848 = buffer.data(qsf + 848);
    const auto *qsf_849 = buffer.data(qsf + 849);
    const auto *qsf_850 = buffer.data(qsf + 850);
    const auto *qsf_851 = buffer.data(qsf + 851);
    const auto *qsf_852 = buffer.data(qsf + 852);
    const auto *qsf_853 = buffer.data(qsf + 853);
    const auto *qsf_854 = buffer.data(qsf + 854);
    const auto *qsf_855 = buffer.data(qsf + 855);
    const auto *qsf_856 = buffer.data(qsf + 856);
    const auto *qsf_857 = buffer.data(qsf + 857);
    const auto *qsf_858 = buffer.data(qsf + 858);
    const auto *qsf_859 = buffer.data(qsf + 859);
    const auto *qsf_860 = buffer.data(qsf + 860);
    const auto *qsf_861 = buffer.data(qsf + 861);
    const auto *qsf_862 = buffer.data(qsf + 862);
    const auto *qsf_863 = buffer.data(qsf + 863);
    const auto *qsf_864 = buffer.data(qsf + 864);
    const auto *qsf_865 = buffer.data(qsf + 865);
    const auto *qsf_866 = buffer.data(qsf + 866);
    const auto *qsf_867 = buffer.data(qsf + 867);
    const auto *qsf_868 = buffer.data(qsf + 868);
    const auto *qsf_869 = buffer.data(qsf + 869);
    const auto *qsf_870 = buffer.data(qsf + 870);
    const auto *qsf_871 = buffer.data(qsf + 871);
    const auto *qsf_872 = buffer.data(qsf + 872);
    const auto *qsf_873 = buffer.data(qsf + 873);
    const auto *qsf_874 = buffer.data(qsf + 874);
    const auto *qsf_875 = buffer.data(qsf + 875);
    const auto *qsf_876 = buffer.data(qsf + 876);
    const auto *qsf_877 = buffer.data(qsf + 877);
    const auto *qsf_878 = buffer.data(qsf + 878);
    const auto *qsf_879 = buffer.data(qsf + 879);

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_z, osf_669, qsd0_479, qsd0_480, \
                         qsd0_481, qsd1_479, qsd1_480, qsd1_481, qsf_799, qsf_800, \
                         qsf_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_7 * osf_669[k]
                    + f_1 * qsd0_479[k]
                    - f_2 * qsd1_479[k]
                    + f_3 * pc_z[k] * qsf_799[k];

        t_1200[k] = f_1 * qsd0_480[k]
                    - f_2 * qsd1_480[k]
                    + f_3 * pc_x[k] * qsf_800[k];

        t_1201[k] = f_10 * qsd0_481[k]
                    - f_11 * qsd1_481[k]
                    + f_3 * pc_x[k] * qsf_801[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, qsd0_482, qsd0_483, qsd0_484, qsd1_482, \
                         qsd1_483, qsd1_484, qsf_802, qsf_803, \
                         qsf_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_10 * qsd0_482[k]
                    - f_11 * qsd1_482[k]
                    + f_3 * pc_x[k] * qsf_802[k];

        t_1203[k] = f_4 * qsd0_483[k]
                    - f_5 * qsd1_483[k]
                    + f_3 * pc_x[k] * qsf_803[k];

        t_1204[k] = f_4 * qsd0_484[k]
                    - f_5 * qsd1_484[k]
                    + f_3 * pc_x[k] * qsf_804[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, t_1209, pc_x, qsd0_485, qsd1_485, \
                         qsf_805, qsf_806, qsf_807, qsf_808, qsf_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_4 * qsd0_485[k]
                    - f_5 * qsd1_485[k]
                    + f_3 * pc_x[k] * qsf_805[k];

        t_1206[k] = f_3 * pc_x[k] * qsf_806[k];

        t_1207[k] = f_3 * pc_x[k] * qsf_807[k];

        t_1208[k] = f_3 * pc_x[k] * qsf_808[k];

        t_1209[k] = f_3 * pc_x[k] * qsf_809[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, pc_y, pc_z, osf_676, osf_686, osf_688, \
                         qsd0_483, qsd0_485, qsd1_483, qsd1_485, qsf_806, \
                         qsf_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_12 * osf_686[k]
                    + f_1 * qsd0_483[k]
                    - f_2 * qsd1_483[k]
                    + f_3 * pc_y[k] * qsf_806[k];

        t_1211[k] = f_8 * osf_676[k]
                    + f_3 * pc_z[k] * qsf_806[k];

        t_1212[k] = f_12 * osf_688[k]
                    + f_4 * qsd0_485[k]
                    - f_5 * qsd1_485[k]
                    + f_3 * pc_y[k] * qsf_808[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, pc_x, pc_y, pc_z, osf_679, osf_689, qsd0_485, \
                         qsd0_486, qsd1_485, qsd1_486, qsf_809, \
                         qsf_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = f_12 * osf_689[k]
                    + f_3 * pc_y[k] * qsf_809[k];

        t_1214[k] = f_8 * osf_679[k]
                    + f_1 * qsd0_485[k]
                    - f_2 * qsd1_485[k]
                    + f_3 * pc_z[k] * qsf_809[k];

        t_1215[k] = f_1 * qsd0_486[k]
                    - f_2 * qsd1_486[k]
                    + f_3 * pc_x[k] * qsf_810[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_x, qsd0_487, qsd0_488, qsd0_489, qsd1_487, \
                         qsd1_488, qsd1_489, qsf_811, qsf_812, \
                         qsf_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_10 * qsd0_487[k]
                    - f_11 * qsd1_487[k]
                    + f_3 * pc_x[k] * qsf_811[k];

        t_1217[k] = f_10 * qsd0_488[k]
                    - f_11 * qsd1_488[k]
                    + f_3 * pc_x[k] * qsf_812[k];

        t_1218[k] = f_4 * qsd0_489[k]
                    - f_5 * qsd1_489[k]
                    + f_3 * pc_x[k] * qsf_813[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, t_1222, t_1223, pc_x, qsd0_490, qsd0_491, \
                         qsd1_490, qsd1_491, qsf_814, qsf_815, qsf_816, qsf_817, \
                         qsf_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_4 * qsd0_490[k]
                    - f_5 * qsd1_490[k]
                    + f_3 * pc_x[k] * qsf_814[k];

        t_1220[k] = f_4 * qsd0_491[k]
                    - f_5 * qsd1_491[k]
                    + f_3 * pc_x[k] * qsf_815[k];

        t_1221[k] = f_3 * pc_x[k] * qsf_816[k];

        t_1222[k] = f_3 * pc_x[k] * qsf_817[k];

        t_1223[k] = f_3 * pc_x[k] * qsf_818[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pc_x, pc_y, pc_z, osf_686, osf_696, qsd0_489, \
                         qsd1_489, qsf_816, qsf_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_3 * pc_x[k] * qsf_819[k];

        t_1225[k] = f_13 * osf_696[k]
                    + f_1 * qsd0_489[k]
                    - f_2 * qsd1_489[k]
                    + f_3 * pc_y[k] * qsf_816[k];

        t_1226[k] = f_14 * osf_686[k]
                    + f_3 * pc_z[k] * qsf_816[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pc_y, pc_z, osf_689, osf_698, osf_699, \
                         qsd0_491, qsd1_491, qsf_818, qsf_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_13 * osf_698[k]
                    + f_4 * qsd0_491[k]
                    - f_5 * qsd1_491[k]
                    + f_3 * pc_y[k] * qsf_818[k];

        t_1228[k] = f_13 * osf_699[k]
                    + f_3 * pc_y[k] * qsf_819[k];

        t_1229[k] = f_14 * osf_689[k]
                    + f_1 * qsd0_491[k]
                    - f_2 * qsd1_491[k]
                    + f_3 * pc_z[k] * qsf_819[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pc_x, qsd0_492, qsd0_493, qsd0_494, qsd1_492, \
                         qsd1_493, qsd1_494, qsf_820, qsf_821, \
                         qsf_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_1 * qsd0_492[k]
                    - f_2 * qsd1_492[k]
                    + f_3 * pc_x[k] * qsf_820[k];

        t_1231[k] = f_10 * qsd0_493[k]
                    - f_11 * qsd1_493[k]
                    + f_3 * pc_x[k] * qsf_821[k];

        t_1232[k] = f_10 * qsd0_494[k]
                    - f_11 * qsd1_494[k]
                    + f_3 * pc_x[k] * qsf_822[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, t_1236, pc_x, qsd0_495, qsd0_496, qsd0_497, \
                         qsd1_495, qsd1_496, qsd1_497, qsf_823, qsf_824, qsf_825, \
                         qsf_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_4 * qsd0_495[k]
                    - f_5 * qsd1_495[k]
                    + f_3 * pc_x[k] * qsf_823[k];

        t_1234[k] = f_4 * qsd0_496[k]
                    - f_5 * qsd1_496[k]
                    + f_3 * pc_x[k] * qsf_824[k];

        t_1235[k] = f_4 * qsd0_497[k]
                    - f_5 * qsd1_497[k]
                    + f_3 * pc_x[k] * qsf_825[k];

        t_1236[k] = f_3 * pc_x[k] * qsf_826[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, t_1241, pc_x, pc_y, pc_z, osf_696, \
                         osf_706, qsd0_495, qsd1_495, qsf_826, qsf_827, qsf_828, \
                         qsf_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_3 * pc_x[k] * qsf_827[k];

        t_1238[k] = f_3 * pc_x[k] * qsf_828[k];

        t_1239[k] = f_3 * pc_x[k] * qsf_829[k];

        t_1240[k] = f_15 * osf_706[k]
                    + f_1 * qsd0_495[k]
                    - f_2 * qsd1_495[k]
                    + f_3 * pc_y[k] * qsf_826[k];

        t_1241[k] = f_16 * osf_696[k]
                    + f_3 * pc_z[k] * qsf_826[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, pc_y, pc_z, osf_699, osf_708, osf_709, \
                         qsd0_497, qsd1_497, qsf_828, qsf_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_15 * osf_708[k]
                    + f_4 * qsd0_497[k]
                    - f_5 * qsd1_497[k]
                    + f_3 * pc_y[k] * qsf_828[k];

        t_1243[k] = f_15 * osf_709[k]
                    + f_3 * pc_y[k] * qsf_829[k];

        t_1244[k] = f_16 * osf_699[k]
                    + f_1 * qsd0_497[k]
                    - f_2 * qsd1_497[k]
                    + f_3 * pc_z[k] * qsf_829[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, pc_x, qsd0_498, qsd0_499, qsd0_500, qsd1_498, \
                         qsd1_499, qsd1_500, qsf_830, qsf_831, \
                         qsf_832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_1 * qsd0_498[k]
                    - f_2 * qsd1_498[k]
                    + f_3 * pc_x[k] * qsf_830[k];

        t_1246[k] = f_10 * qsd0_499[k]
                    - f_11 * qsd1_499[k]
                    + f_3 * pc_x[k] * qsf_831[k];

        t_1247[k] = f_10 * qsd0_500[k]
                    - f_11 * qsd1_500[k]
                    + f_3 * pc_x[k] * qsf_832[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, pc_x, qsd0_501, qsd0_502, qsd0_503, \
                         qsd1_501, qsd1_502, qsd1_503, qsf_833, qsf_834, qsf_835, \
                         qsf_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_4 * qsd0_501[k]
                    - f_5 * qsd1_501[k]
                    + f_3 * pc_x[k] * qsf_833[k];

        t_1249[k] = f_4 * qsd0_502[k]
                    - f_5 * qsd1_502[k]
                    + f_3 * pc_x[k] * qsf_834[k];

        t_1250[k] = f_4 * qsd0_503[k]
                    - f_5 * qsd1_503[k]
                    + f_3 * pc_x[k] * qsf_835[k];

        t_1251[k] = f_3 * pc_x[k] * qsf_836[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, t_1256, pc_x, pc_y, pc_z, osf_706, \
                         osf_716, qsd0_501, qsd1_501, qsf_836, qsf_837, qsf_838, \
                         qsf_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_3 * pc_x[k] * qsf_837[k];

        t_1253[k] = f_3 * pc_x[k] * qsf_838[k];

        t_1254[k] = f_3 * pc_x[k] * qsf_839[k];

        t_1255[k] = f_17 * osf_716[k]
                    + f_1 * qsd0_501[k]
                    - f_2 * qsd1_501[k]
                    + f_3 * pc_y[k] * qsf_836[k];

        t_1256[k] = f_18 * osf_706[k]
                    + f_3 * pc_z[k] * qsf_836[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pc_y, pc_z, osf_709, osf_718, osf_719, \
                         qsd0_503, qsd1_503, qsf_838, qsf_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_17 * osf_718[k]
                    + f_4 * qsd0_503[k]
                    - f_5 * qsd1_503[k]
                    + f_3 * pc_y[k] * qsf_838[k];

        t_1258[k] = f_17 * osf_719[k]
                    + f_3 * pc_y[k] * qsf_839[k];

        t_1259[k] = f_18 * osf_709[k]
                    + f_1 * qsd0_503[k]
                    - f_2 * qsd1_503[k]
                    + f_3 * pc_z[k] * qsf_839[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, pc_x, qsd0_504, qsd0_505, qsd0_506, qsd1_504, \
                         qsd1_505, qsd1_506, qsf_840, qsf_841, \
                         qsf_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_1 * qsd0_504[k]
                    - f_2 * qsd1_504[k]
                    + f_3 * pc_x[k] * qsf_840[k];

        t_1261[k] = f_10 * qsd0_505[k]
                    - f_11 * qsd1_505[k]
                    + f_3 * pc_x[k] * qsf_841[k];

        t_1262[k] = f_10 * qsd0_506[k]
                    - f_11 * qsd1_506[k]
                    + f_3 * pc_x[k] * qsf_842[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, t_1266, pc_x, qsd0_507, qsd0_508, qsd0_509, \
                         qsd1_507, qsd1_508, qsd1_509, qsf_843, qsf_844, qsf_845, \
                         qsf_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_4 * qsd0_507[k]
                    - f_5 * qsd1_507[k]
                    + f_3 * pc_x[k] * qsf_843[k];

        t_1264[k] = f_4 * qsd0_508[k]
                    - f_5 * qsd1_508[k]
                    + f_3 * pc_x[k] * qsf_844[k];

        t_1265[k] = f_4 * qsd0_509[k]
                    - f_5 * qsd1_509[k]
                    + f_3 * pc_x[k] * qsf_845[k];

        t_1266[k] = f_3 * pc_x[k] * qsf_846[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, t_1270, t_1271, pc_x, pc_y, pc_z, osf_716, \
                         osf_726, qsd0_507, qsd1_507, qsf_846, qsf_847, qsf_848, \
                         qsf_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_3 * pc_x[k] * qsf_847[k];

        t_1268[k] = f_3 * pc_x[k] * qsf_848[k];

        t_1269[k] = f_3 * pc_x[k] * qsf_849[k];

        t_1270[k] = f_19 * osf_726[k]
                    + f_1 * qsd0_507[k]
                    - f_2 * qsd1_507[k]
                    + f_3 * pc_y[k] * qsf_846[k];

        t_1271[k] = f_19 * osf_716[k]
                    + f_3 * pc_z[k] * qsf_846[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pc_y, pc_z, osf_719, osf_728, osf_729, \
                         qsd0_509, qsd1_509, qsf_848, qsf_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_19 * osf_728[k]
                    + f_4 * qsd0_509[k]
                    - f_5 * qsd1_509[k]
                    + f_3 * pc_y[k] * qsf_848[k];

        t_1273[k] = f_19 * osf_729[k]
                    + f_3 * pc_y[k] * qsf_849[k];

        t_1274[k] = f_19 * osf_719[k]
                    + f_1 * qsd0_509[k]
                    - f_2 * qsd1_509[k]
                    + f_3 * pc_z[k] * qsf_849[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pc_x, qsd0_510, qsd0_511, qsd0_512, qsd1_510, \
                         qsd1_511, qsd1_512, qsf_850, qsf_851, \
                         qsf_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_1 * qsd0_510[k]
                    - f_2 * qsd1_510[k]
                    + f_3 * pc_x[k] * qsf_850[k];

        t_1276[k] = f_10 * qsd0_511[k]
                    - f_11 * qsd1_511[k]
                    + f_3 * pc_x[k] * qsf_851[k];

        t_1277[k] = f_10 * qsd0_512[k]
                    - f_11 * qsd1_512[k]
                    + f_3 * pc_x[k] * qsf_852[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, t_1281, pc_x, qsd0_513, qsd0_514, qsd0_515, \
                         qsd1_513, qsd1_514, qsd1_515, qsf_853, qsf_854, qsf_855, \
                         qsf_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_4 * qsd0_513[k]
                    - f_5 * qsd1_513[k]
                    + f_3 * pc_x[k] * qsf_853[k];

        t_1279[k] = f_4 * qsd0_514[k]
                    - f_5 * qsd1_514[k]
                    + f_3 * pc_x[k] * qsf_854[k];

        t_1280[k] = f_4 * qsd0_515[k]
                    - f_5 * qsd1_515[k]
                    + f_3 * pc_x[k] * qsf_855[k];

        t_1281[k] = f_3 * pc_x[k] * qsf_856[k];
    }

#pragma omp simd aligned(t_1282, t_1283, t_1284, t_1285, t_1286, pc_x, pc_y, pc_z, osf_726, \
                         osf_736, qsd0_513, qsd1_513, qsf_856, qsf_857, qsf_858, \
                         qsf_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1282[k] = f_3 * pc_x[k] * qsf_857[k];

        t_1283[k] = f_3 * pc_x[k] * qsf_858[k];

        t_1284[k] = f_3 * pc_x[k] * qsf_859[k];

        t_1285[k] = f_18 * osf_736[k]
                    + f_1 * qsd0_513[k]
                    - f_2 * qsd1_513[k]
                    + f_3 * pc_y[k] * qsf_856[k];

        t_1286[k] = f_17 * osf_726[k]
                    + f_3 * pc_z[k] * qsf_856[k];
    }

#pragma omp simd aligned(t_1287, t_1288, t_1289, pc_y, pc_z, osf_729, osf_738, osf_739, \
                         qsd0_515, qsd1_515, qsf_858, qsf_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1287[k] = f_18 * osf_738[k]
                    + f_4 * qsd0_515[k]
                    - f_5 * qsd1_515[k]
                    + f_3 * pc_y[k] * qsf_858[k];

        t_1288[k] = f_18 * osf_739[k]
                    + f_3 * pc_y[k] * qsf_859[k];

        t_1289[k] = f_17 * osf_729[k]
                    + f_1 * qsd0_515[k]
                    - f_2 * qsd1_515[k]
                    + f_3 * pc_z[k] * qsf_859[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, pc_x, qsd0_516, qsd0_517, qsd0_518, qsd1_516, \
                         qsd1_517, qsd1_518, qsf_860, qsf_861, \
                         qsf_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_1 * qsd0_516[k]
                    - f_2 * qsd1_516[k]
                    + f_3 * pc_x[k] * qsf_860[k];

        t_1291[k] = f_10 * qsd0_517[k]
                    - f_11 * qsd1_517[k]
                    + f_3 * pc_x[k] * qsf_861[k];

        t_1292[k] = f_10 * qsd0_518[k]
                    - f_11 * qsd1_518[k]
                    + f_3 * pc_x[k] * qsf_862[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pc_x, qsd0_519, qsd0_520, qsd0_521, \
                         qsd1_519, qsd1_520, qsd1_521, qsf_863, qsf_864, qsf_865, \
                         qsf_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_4 * qsd0_519[k]
                    - f_5 * qsd1_519[k]
                    + f_3 * pc_x[k] * qsf_863[k];

        t_1294[k] = f_4 * qsd0_520[k]
                    - f_5 * qsd1_520[k]
                    + f_3 * pc_x[k] * qsf_864[k];

        t_1295[k] = f_4 * qsd0_521[k]
                    - f_5 * qsd1_521[k]
                    + f_3 * pc_x[k] * qsf_865[k];

        t_1296[k] = f_3 * pc_x[k] * qsf_866[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, t_1300, t_1301, pc_x, pc_y, pc_z, osf_736, \
                         osf_746, qsd0_519, qsd1_519, qsf_866, qsf_867, qsf_868, \
                         qsf_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = f_3 * pc_x[k] * qsf_867[k];

        t_1298[k] = f_3 * pc_x[k] * qsf_868[k];

        t_1299[k] = f_3 * pc_x[k] * qsf_869[k];

        t_1300[k] = f_16 * osf_746[k]
                    + f_1 * qsd0_519[k]
                    - f_2 * qsd1_519[k]
                    + f_3 * pc_y[k] * qsf_866[k];

        t_1301[k] = f_15 * osf_736[k]
                    + f_3 * pc_z[k] * qsf_866[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, pc_y, pc_z, osf_739, osf_748, osf_749, \
                         qsd0_521, qsd1_521, qsf_868, qsf_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_16 * osf_748[k]
                    + f_4 * qsd0_521[k]
                    - f_5 * qsd1_521[k]
                    + f_3 * pc_y[k] * qsf_868[k];

        t_1303[k] = f_16 * osf_749[k]
                    + f_3 * pc_y[k] * qsf_869[k];

        t_1304[k] = f_15 * osf_739[k]
                    + f_1 * qsd0_521[k]
                    - f_2 * qsd1_521[k]
                    + f_3 * pc_z[k] * qsf_869[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, pc_x, qsd0_522, qsd0_523, qsd0_524, qsd1_522, \
                         qsd1_523, qsd1_524, qsf_870, qsf_871, \
                         qsf_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = f_1 * qsd0_522[k]
                    - f_2 * qsd1_522[k]
                    + f_3 * pc_x[k] * qsf_870[k];

        t_1306[k] = f_10 * qsd0_523[k]
                    - f_11 * qsd1_523[k]
                    + f_3 * pc_x[k] * qsf_871[k];

        t_1307[k] = f_10 * qsd0_524[k]
                    - f_11 * qsd1_524[k]
                    + f_3 * pc_x[k] * qsf_872[k];
    }

#pragma omp simd aligned(t_1308, t_1309, t_1310, t_1311, pc_x, qsd0_525, qsd0_526, qsd0_527, \
                         qsd1_525, qsd1_526, qsd1_527, qsf_873, qsf_874, qsf_875, \
                         qsf_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1308[k] = f_4 * qsd0_525[k]
                    - f_5 * qsd1_525[k]
                    + f_3 * pc_x[k] * qsf_873[k];

        t_1309[k] = f_4 * qsd0_526[k]
                    - f_5 * qsd1_526[k]
                    + f_3 * pc_x[k] * qsf_874[k];

        t_1310[k] = f_4 * qsd0_527[k]
                    - f_5 * qsd1_527[k]
                    + f_3 * pc_x[k] * qsf_875[k];

        t_1311[k] = f_3 * pc_x[k] * qsf_876[k];
    }

#pragma omp simd aligned(t_1312, t_1313, t_1314, t_1315, t_1316, pc_x, pc_y, pc_z, osf_746, \
                         osf_756, qsd0_525, qsd1_525, qsf_876, qsf_877, qsf_878, \
                         qsf_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1312[k] = f_3 * pc_x[k] * qsf_877[k];

        t_1313[k] = f_3 * pc_x[k] * qsf_878[k];

        t_1314[k] = f_3 * pc_x[k] * qsf_879[k];

        t_1315[k] = f_14 * osf_756[k]
                    + f_1 * qsd0_525[k]
                    - f_2 * qsd1_525[k]
                    + f_3 * pc_y[k] * qsf_876[k];

        t_1316[k] = f_13 * osf_746[k]
                    + f_3 * pc_z[k] * qsf_876[k];
    }

#pragma omp simd aligned(t_1317, t_1318, t_1319, pc_y, pc_z, osf_749, osf_758, osf_759, \
                         qsd0_527, qsd1_527, qsf_878, qsf_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1317[k] = f_14 * osf_758[k]
                    + f_4 * qsd0_527[k]
                    - f_5 * qsd1_527[k]
                    + f_3 * pc_y[k] * qsf_878[k];

        t_1318[k] = f_14 * osf_759[k]
                    + f_3 * pc_y[k] * qsf_879[k];

        t_1319[k] = f_13 * osf_749[k]
                    + f_1 * qsd0_527[k]
                    - f_2 * qsd1_527[k]
                    + f_3 * pc_z[k] * qsf_879[k];
    }
}

static auto
compute_prim_qsg_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osg0,
                                                           const size_t osf, const size_t osg1,
                                                           const size_t qsd0, const size_t qsd1,
                                                           const size_t qsf, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 5.0 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osg0_1155 = buffer.data(osg0 + 1155);
    const auto *osg0_1157 = buffer.data(osg0 + 1157);
    const auto *osg0_1160 = buffer.data(osg0 + 1160);
    const auto *osg0_1165 = buffer.data(osg0 + 1165);
    const auto *osg0_1167 = buffer.data(osg0 + 1167);
    const auto *osg0_1169 = buffer.data(osg0 + 1169);

    const auto *osf_756 = buffer.data(osf + 756);
    const auto *osf_759 = buffer.data(osf + 759);
    const auto *osf_766 = buffer.data(osf + 766);
    const auto *osf_768 = buffer.data(osf + 768);
    const auto *osf_769 = buffer.data(osf + 769);
    const auto *osf_776 = buffer.data(osf + 776);
    const auto *osf_778 = buffer.data(osf + 778);
    const auto *osf_779 = buffer.data(osf + 779);

    const auto *osg1_1155 = buffer.data(osg1 + 1155);
    const auto *osg1_1157 = buffer.data(osg1 + 1157);
    const auto *osg1_1160 = buffer.data(osg1 + 1160);
    const auto *osg1_1165 = buffer.data(osg1 + 1165);
    const auto *osg1_1167 = buffer.data(osg1 + 1167);
    const auto *osg1_1169 = buffer.data(osg1 + 1169);

    const auto *qsd0_528 = buffer.data(qsd0 + 528);
    const auto *qsd0_529 = buffer.data(qsd0 + 529);
    const auto *qsd0_530 = buffer.data(qsd0 + 530);
    const auto *qsd0_531 = buffer.data(qsd0 + 531);
    const auto *qsd0_532 = buffer.data(qsd0 + 532);
    const auto *qsd0_533 = buffer.data(qsd0 + 533);
    const auto *qsd0_535 = buffer.data(qsd0 + 535);
    const auto *qsd0_537 = buffer.data(qsd0 + 537);
    const auto *qsd0_538 = buffer.data(qsd0 + 538);
    const auto *qsd0_540 = buffer.data(qsd0 + 540);
    const auto *qsd0_542 = buffer.data(qsd0 + 542);
    const auto *qsd0_543 = buffer.data(qsd0 + 543);
    const auto *qsd0_544 = buffer.data(qsd0 + 544);
    const auto *qsd0_545 = buffer.data(qsd0 + 545);

    const auto *qsd1_528 = buffer.data(qsd1 + 528);
    const auto *qsd1_529 = buffer.data(qsd1 + 529);
    const auto *qsd1_530 = buffer.data(qsd1 + 530);
    const auto *qsd1_531 = buffer.data(qsd1 + 531);
    const auto *qsd1_532 = buffer.data(qsd1 + 532);
    const auto *qsd1_533 = buffer.data(qsd1 + 533);
    const auto *qsd1_535 = buffer.data(qsd1 + 535);
    const auto *qsd1_537 = buffer.data(qsd1 + 537);
    const auto *qsd1_538 = buffer.data(qsd1 + 538);
    const auto *qsd1_540 = buffer.data(qsd1 + 540);
    const auto *qsd1_542 = buffer.data(qsd1 + 542);
    const auto *qsd1_543 = buffer.data(qsd1 + 543);
    const auto *qsd1_544 = buffer.data(qsd1 + 544);
    const auto *qsd1_545 = buffer.data(qsd1 + 545);

    const auto *qsf_880 = buffer.data(qsf + 880);
    const auto *qsf_881 = buffer.data(qsf + 881);
    const auto *qsf_882 = buffer.data(qsf + 882);
    const auto *qsf_883 = buffer.data(qsf + 883);
    const auto *qsf_884 = buffer.data(qsf + 884);
    const auto *qsf_885 = buffer.data(qsf + 885);
    const auto *qsf_886 = buffer.data(qsf + 886);
    const auto *qsf_887 = buffer.data(qsf + 887);
    const auto *qsf_888 = buffer.data(qsf + 888);
    const auto *qsf_889 = buffer.data(qsf + 889);
    const auto *qsf_891 = buffer.data(qsf + 891);
    const auto *qsf_893 = buffer.data(qsf + 893);
    const auto *qsf_894 = buffer.data(qsf + 894);
    const auto *qsf_896 = buffer.data(qsf + 896);
    const auto *qsf_897 = buffer.data(qsf + 897);
    const auto *qsf_898 = buffer.data(qsf + 898);
    const auto *qsf_899 = buffer.data(qsf + 899);
    const auto *qsf_900 = buffer.data(qsf + 900);
    const auto *qsf_902 = buffer.data(qsf + 902);
    const auto *qsf_903 = buffer.data(qsf + 903);
    const auto *qsf_905 = buffer.data(qsf + 905);
    const auto *qsf_906 = buffer.data(qsf + 906);
    const auto *qsf_907 = buffer.data(qsf + 907);
    const auto *qsf_908 = buffer.data(qsf + 908);
    const auto *qsf_909 = buffer.data(qsf + 909);

#pragma omp simd aligned(t_1320, t_1321, t_1322, pc_x, qsd0_528, qsd0_529, qsd0_530, qsd1_528, \
                         qsd1_529, qsd1_530, qsf_880, qsf_881, \
                         qsf_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1320[k] = f_1 * qsd0_528[k]
                    - f_2 * qsd1_528[k]
                    + f_3 * pc_x[k] * qsf_880[k];

        t_1321[k] = f_10 * qsd0_529[k]
                    - f_11 * qsd1_529[k]
                    + f_3 * pc_x[k] * qsf_881[k];

        t_1322[k] = f_10 * qsd0_530[k]
                    - f_11 * qsd1_530[k]
                    + f_3 * pc_x[k] * qsf_882[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, t_1326, pc_x, qsd0_531, qsd0_532, qsd0_533, \
                         qsd1_531, qsd1_532, qsd1_533, qsf_883, qsf_884, qsf_885, \
                         qsf_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = f_4 * qsd0_531[k]
                    - f_5 * qsd1_531[k]
                    + f_3 * pc_x[k] * qsf_883[k];

        t_1324[k] = f_4 * qsd0_532[k]
                    - f_5 * qsd1_532[k]
                    + f_3 * pc_x[k] * qsf_884[k];

        t_1325[k] = f_4 * qsd0_533[k]
                    - f_5 * qsd1_533[k]
                    + f_3 * pc_x[k] * qsf_885[k];

        t_1326[k] = f_3 * pc_x[k] * qsf_886[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, t_1330, t_1331, pc_x, pc_y, pc_z, osf_756, \
                         osf_766, qsd0_531, qsd1_531, qsf_886, qsf_887, qsf_888, \
                         qsf_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_3 * pc_x[k] * qsf_887[k];

        t_1328[k] = f_3 * pc_x[k] * qsf_888[k];

        t_1329[k] = f_3 * pc_x[k] * qsf_889[k];

        t_1330[k] = f_8 * osf_766[k]
                    + f_1 * qsd0_531[k]
                    - f_2 * qsd1_531[k]
                    + f_3 * pc_y[k] * qsf_886[k];

        t_1331[k] = f_12 * osf_756[k]
                    + f_3 * pc_z[k] * qsf_886[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, pa_y, pc_y, pc_z, osg0_1155, osf_759, \
                         osf_768, osf_769, osg1_1155, qsd0_533, qsd1_533, qsf_888, \
                         qsf_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = f_8 * osf_768[k]
                    + f_4 * qsd0_533[k]
                    - f_5 * qsd1_533[k]
                    + f_3 * pc_y[k] * qsf_888[k];

        t_1333[k] = f_8 * osf_769[k]
                    + f_3 * pc_y[k] * qsf_889[k];

        t_1334[k] = f_12 * osf_759[k]
                    + f_1 * qsd0_533[k]
                    - f_2 * qsd1_533[k]
                    + f_3 * pc_z[k] * qsf_889[k];

        t_1335[k] = pa_y[k] * osg0_1155[k]
                    - f_6 * pc_y[k] * osg1_1155[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, pa_y, pc_x, pc_y, osg0_1157, osg1_1157, \
                         qsd0_535, qsd0_537, qsd1_535, qsd1_537, qsf_891, \
                         qsf_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_10 * qsd0_535[k]
                    - f_11 * qsd1_535[k]
                    + f_3 * pc_x[k] * qsf_891[k];

        t_1337[k] = pa_y[k] * osg0_1157[k]
                    - f_6 * pc_y[k] * osg1_1157[k];

        t_1338[k] = f_4 * qsd0_537[k]
                    - f_5 * qsd1_537[k]
                    + f_3 * pc_x[k] * qsf_893[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, t_1342, t_1343, pa_y, pc_x, pc_y, osg0_1160, \
                         osg1_1160, qsd0_538, qsd1_538, qsf_894, qsf_896, qsf_897, \
                         qsf_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_4 * qsd0_538[k]
                    - f_5 * qsd1_538[k]
                    + f_3 * pc_x[k] * qsf_894[k];

        t_1340[k] = pa_y[k] * osg0_1160[k]
                    - f_6 * pc_y[k] * osg1_1160[k];

        t_1341[k] = f_3 * pc_x[k] * qsf_896[k];

        t_1342[k] = f_3 * pc_x[k] * qsf_897[k];

        t_1343[k] = f_3 * pc_x[k] * qsf_898[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, pa_y, pc_x, pc_y, pc_z, osg0_1165, osf_766, \
                         osf_776, osg1_1165, qsf_896, qsf_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = f_3 * pc_x[k] * qsf_899[k];

        t_1345[k] = pa_y[k] * osg0_1165[k]
                    + f_16 * osf_776[k]
                    - f_6 * pc_y[k] * osg1_1165[k];

        t_1346[k] = f_9 * osf_766[k]
                    + f_3 * pc_z[k] * qsf_896[k];
    }

#pragma omp simd aligned(t_1347, t_1348, t_1349, pa_y, pc_y, osg0_1167, osg0_1169, osf_778, \
                         osf_779, osg1_1167, osg1_1169, qsf_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1347[k] = pa_y[k] * osg0_1167[k]
                    + f_8 * osf_778[k]
                    - f_6 * pc_y[k] * osg1_1167[k];

        t_1348[k] = f_7 * osf_779[k]
                    + f_3 * pc_y[k] * qsf_899[k];

        t_1349[k] = pa_y[k] * osg0_1169[k]
                    - f_6 * pc_y[k] * osg1_1169[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, t_1353, t_1354, pc_x, pc_y, qsd0_540, \
                         qsd0_542, qsd0_543, qsd1_540, qsd1_542, qsd1_543, qsf_900, qsf_902, \
                         qsf_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = f_1 * qsd0_540[k]
                    - f_2 * qsd1_540[k]
                    + f_3 * pc_x[k] * qsf_900[k];

        t_1351[k] = f_3 * pc_y[k] * qsf_900[k];

        t_1352[k] = f_10 * qsd0_542[k]
                    - f_11 * qsd1_542[k]
                    + f_3 * pc_x[k] * qsf_902[k];

        t_1353[k] = f_4 * qsd0_543[k]
                    - f_5 * qsd1_543[k]
                    + f_3 * pc_x[k] * qsf_903[k];

        t_1354[k] = f_3 * pc_y[k] * qsf_902[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, t_1358, t_1359, pc_x, qsd0_545, qsd1_545, \
                         qsf_905, qsf_906, qsf_907, qsf_908, qsf_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_4 * qsd0_545[k]
                    - f_5 * qsd1_545[k]
                    + f_3 * pc_x[k] * qsf_905[k];

        t_1356[k] = f_3 * pc_x[k] * qsf_906[k];

        t_1357[k] = f_3 * pc_x[k] * qsf_907[k];

        t_1358[k] = f_3 * pc_x[k] * qsf_908[k];

        t_1359[k] = f_3 * pc_x[k] * qsf_909[k];
    }

#pragma omp simd aligned(t_1360, t_1361, t_1362, t_1363, pc_y, qsd0_543, qsd0_544, qsd0_545, \
                         qsd1_543, qsd1_544, qsd1_545, qsf_906, qsf_907, qsf_908, \
                         qsf_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1360[k] = f_1 * qsd0_543[k]
                    - f_2 * qsd1_543[k]
                    + f_3 * pc_y[k] * qsf_906[k];

        t_1361[k] = f_10 * qsd0_544[k]
                    - f_11 * qsd1_544[k]
                    + f_3 * pc_y[k] * qsf_907[k];

        t_1362[k] = f_4 * qsd0_545[k]
                    - f_5 * qsd1_545[k]
                    + f_3 * pc_y[k] * qsf_908[k];

        t_1363[k] = f_3 * pc_y[k] * qsf_909[k];
    }

#pragma omp simd aligned(t_1364, pc_z, osf_779, qsd0_545, qsd1_545, \
                         qsf_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = f_0 * osf_779[k]
                    + f_1 * qsd0_545[k]
                    - f_2 * qsd1_545[k]
                    + f_3 * pc_z[k] * qsf_909[k];
    }
}

auto
compute_prim_qsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t osg0, const size_t osf,
                                                   const size_t osg1, const size_t qsd0,
                                                   const size_t qsd1, const size_t qsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_qsg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, osg0, osf,
                                                              osg1, qsd0, qsd1, qsf, ncols,
                                                              gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece10(buffer, target, pc, osf, qsd0,
                                                               qsd1, qsf, ncols, gamma, p, q);

    compute_prim_qsg_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, osg0,
                                                               osf, osg1, qsd0, qsd1, qsf,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
