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


#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qsd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osd0,
                                                          const size_t osp, const size_t osd1,
                                                          const size_t qss0, const size_t qss1,
                                                          const size_t qsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 5.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osd0_0 = buffer.data(osd0 + 0);
    const auto *osd0_3 = buffer.data(osd0 + 3);
    const auto *osd0_5 = buffer.data(osd0 + 5);
    const auto *osd0_9 = buffer.data(osd0 + 9);
    const auto *osd0_12 = buffer.data(osd0 + 12);
    const auto *osd0_17 = buffer.data(osd0 + 17);
    const auto *osd0_18 = buffer.data(osd0 + 18);
    const auto *osd0_21 = buffer.data(osd0 + 21);
    const auto *osd0_30 = buffer.data(osd0 + 30);
    const auto *osd0_35 = buffer.data(osd0 + 35);
    const auto *osd0_36 = buffer.data(osd0 + 36);
    const auto *osd0_39 = buffer.data(osd0 + 39);
    const auto *osd0_54 = buffer.data(osd0 + 54);
    const auto *osd0_59 = buffer.data(osd0 + 59);
    const auto *osd0_60 = buffer.data(osd0 + 60);
    const auto *osd0_63 = buffer.data(osd0 + 63);
    const auto *osd0_84 = buffer.data(osd0 + 84);
    const auto *osd0_89 = buffer.data(osd0 + 89);

    const auto *osp_0 = buffer.data(osp + 0);
    const auto *osp_1 = buffer.data(osp + 1);
    const auto *osp_2 = buffer.data(osp + 2);
    const auto *osp_4 = buffer.data(osp + 4);
    const auto *osp_8 = buffer.data(osp + 8);
    const auto *osp_9 = buffer.data(osp + 9);
    const auto *osp_10 = buffer.data(osp + 10);
    const auto *osp_11 = buffer.data(osp + 11);
    const auto *osp_13 = buffer.data(osp + 13);
    const auto *osp_14 = buffer.data(osp + 14);
    const auto *osp_15 = buffer.data(osp + 15);
    const auto *osp_16 = buffer.data(osp + 16);
    const auto *osp_17 = buffer.data(osp + 17);
    const auto *osp_18 = buffer.data(osp + 18);
    const auto *osp_19 = buffer.data(osp + 19);
    const auto *osp_20 = buffer.data(osp + 20);
    const auto *osp_22 = buffer.data(osp + 22);
    const auto *osp_23 = buffer.data(osp + 23);
    const auto *osp_25 = buffer.data(osp + 25);
    const auto *osp_26 = buffer.data(osp + 26);
    const auto *osp_27 = buffer.data(osp + 27);
    const auto *osp_28 = buffer.data(osp + 28);
    const auto *osp_29 = buffer.data(osp + 29);
    const auto *osp_30 = buffer.data(osp + 30);
    const auto *osp_31 = buffer.data(osp + 31);
    const auto *osp_32 = buffer.data(osp + 32);
    const auto *osp_34 = buffer.data(osp + 34);
    const auto *osp_35 = buffer.data(osp + 35);
    const auto *osp_36 = buffer.data(osp + 36);
    const auto *osp_37 = buffer.data(osp + 37);
    const auto *osp_38 = buffer.data(osp + 38);
    const auto *osp_40 = buffer.data(osp + 40);
    const auto *osp_41 = buffer.data(osp + 41);
    const auto *osp_42 = buffer.data(osp + 42);
    const auto *osp_43 = buffer.data(osp + 43);
    const auto *osp_44 = buffer.data(osp + 44);
    const auto *osp_45 = buffer.data(osp + 45);
    const auto *osp_46 = buffer.data(osp + 46);
    const auto *osp_49 = buffer.data(osp + 49);
    const auto *osp_50 = buffer.data(osp + 50);
    const auto *osp_51 = buffer.data(osp + 51);
    const auto *osp_52 = buffer.data(osp + 52);
    const auto *osp_53 = buffer.data(osp + 53);
    const auto *osp_54 = buffer.data(osp + 54);
    const auto *osp_55 = buffer.data(osp + 55);
    const auto *osp_56 = buffer.data(osp + 56);
    const auto *osp_58 = buffer.data(osp + 58);
    const auto *osp_59 = buffer.data(osp + 59);
    const auto *osp_60 = buffer.data(osp + 60);
    const auto *osp_62 = buffer.data(osp + 62);

    const auto *osd1_0 = buffer.data(osd1 + 0);
    const auto *osd1_3 = buffer.data(osd1 + 3);
    const auto *osd1_5 = buffer.data(osd1 + 5);
    const auto *osd1_9 = buffer.data(osd1 + 9);
    const auto *osd1_12 = buffer.data(osd1 + 12);
    const auto *osd1_17 = buffer.data(osd1 + 17);
    const auto *osd1_18 = buffer.data(osd1 + 18);
    const auto *osd1_21 = buffer.data(osd1 + 21);
    const auto *osd1_30 = buffer.data(osd1 + 30);
    const auto *osd1_35 = buffer.data(osd1 + 35);
    const auto *osd1_36 = buffer.data(osd1 + 36);
    const auto *osd1_39 = buffer.data(osd1 + 39);
    const auto *osd1_54 = buffer.data(osd1 + 54);
    const auto *osd1_59 = buffer.data(osd1 + 59);
    const auto *osd1_60 = buffer.data(osd1 + 60);
    const auto *osd1_63 = buffer.data(osd1 + 63);
    const auto *osd1_84 = buffer.data(osd1 + 84);
    const auto *osd1_89 = buffer.data(osd1 + 89);

    const auto *qss0_0 = buffer.data(qss0 + 0);
    const auto *qss0_1 = buffer.data(qss0 + 1);
    const auto *qss0_2 = buffer.data(qss0 + 2);
    const auto *qss0_3 = buffer.data(qss0 + 3);
    const auto *qss0_5 = buffer.data(qss0 + 5);
    const auto *qss0_6 = buffer.data(qss0 + 6);
    const auto *qss0_7 = buffer.data(qss0 + 7);
    const auto *qss0_8 = buffer.data(qss0 + 8);
    const auto *qss0_9 = buffer.data(qss0 + 9);
    const auto *qss0_10 = buffer.data(qss0 + 10);
    const auto *qss0_11 = buffer.data(qss0 + 11);
    const auto *qss0_12 = buffer.data(qss0 + 12);
    const auto *qss0_13 = buffer.data(qss0 + 13);
    const auto *qss0_14 = buffer.data(qss0 + 14);
    const auto *qss0_15 = buffer.data(qss0 + 15);
    const auto *qss0_16 = buffer.data(qss0 + 16);
    const auto *qss0_17 = buffer.data(qss0 + 17);
    const auto *qss0_18 = buffer.data(qss0 + 18);
    const auto *qss0_19 = buffer.data(qss0 + 19);
    const auto *qss0_20 = buffer.data(qss0 + 20);

    const auto *qss1_0 = buffer.data(qss1 + 0);
    const auto *qss1_1 = buffer.data(qss1 + 1);
    const auto *qss1_2 = buffer.data(qss1 + 2);
    const auto *qss1_3 = buffer.data(qss1 + 3);
    const auto *qss1_5 = buffer.data(qss1 + 5);
    const auto *qss1_6 = buffer.data(qss1 + 6);
    const auto *qss1_7 = buffer.data(qss1 + 7);
    const auto *qss1_8 = buffer.data(qss1 + 8);
    const auto *qss1_9 = buffer.data(qss1 + 9);
    const auto *qss1_10 = buffer.data(qss1 + 10);
    const auto *qss1_11 = buffer.data(qss1 + 11);
    const auto *qss1_12 = buffer.data(qss1 + 12);
    const auto *qss1_13 = buffer.data(qss1 + 13);
    const auto *qss1_14 = buffer.data(qss1 + 14);
    const auto *qss1_15 = buffer.data(qss1 + 15);
    const auto *qss1_16 = buffer.data(qss1 + 16);
    const auto *qss1_17 = buffer.data(qss1 + 17);
    const auto *qss1_18 = buffer.data(qss1 + 18);
    const auto *qss1_19 = buffer.data(qss1 + 19);
    const auto *qss1_20 = buffer.data(qss1 + 20);

    const auto *qsp_0 = buffer.data(qsp + 0);
    const auto *qsp_1 = buffer.data(qsp + 1);
    const auto *qsp_2 = buffer.data(qsp + 2);
    const auto *qsp_3 = buffer.data(qsp + 3);
    const auto *qsp_4 = buffer.data(qsp + 4);
    const auto *qsp_6 = buffer.data(qsp + 6);
    const auto *qsp_8 = buffer.data(qsp + 8);
    const auto *qsp_9 = buffer.data(qsp + 9);
    const auto *qsp_10 = buffer.data(qsp + 10);
    const auto *qsp_11 = buffer.data(qsp + 11);
    const auto *qsp_13 = buffer.data(qsp + 13);
    const auto *qsp_14 = buffer.data(qsp + 14);
    const auto *qsp_15 = buffer.data(qsp + 15);
    const auto *qsp_16 = buffer.data(qsp + 16);
    const auto *qsp_17 = buffer.data(qsp + 17);
    const auto *qsp_18 = buffer.data(qsp + 18);
    const auto *qsp_19 = buffer.data(qsp + 19);
    const auto *qsp_20 = buffer.data(qsp + 20);
    const auto *qsp_22 = buffer.data(qsp + 22);
    const auto *qsp_23 = buffer.data(qsp + 23);
    const auto *qsp_25 = buffer.data(qsp + 25);
    const auto *qsp_26 = buffer.data(qsp + 26);
    const auto *qsp_27 = buffer.data(qsp + 27);
    const auto *qsp_28 = buffer.data(qsp + 28);
    const auto *qsp_29 = buffer.data(qsp + 29);
    const auto *qsp_30 = buffer.data(qsp + 30);
    const auto *qsp_31 = buffer.data(qsp + 31);
    const auto *qsp_32 = buffer.data(qsp + 32);
    const auto *qsp_34 = buffer.data(qsp + 34);
    const auto *qsp_35 = buffer.data(qsp + 35);
    const auto *qsp_36 = buffer.data(qsp + 36);
    const auto *qsp_37 = buffer.data(qsp + 37);
    const auto *qsp_38 = buffer.data(qsp + 38);
    const auto *qsp_40 = buffer.data(qsp + 40);
    const auto *qsp_41 = buffer.data(qsp + 41);
    const auto *qsp_42 = buffer.data(qsp + 42);
    const auto *qsp_43 = buffer.data(qsp + 43);
    const auto *qsp_44 = buffer.data(qsp + 44);
    const auto *qsp_45 = buffer.data(qsp + 45);
    const auto *qsp_46 = buffer.data(qsp + 46);
    const auto *qsp_47 = buffer.data(qsp + 47);
    const auto *qsp_49 = buffer.data(qsp + 49);
    const auto *qsp_50 = buffer.data(qsp + 50);
    const auto *qsp_51 = buffer.data(qsp + 51);
    const auto *qsp_52 = buffer.data(qsp + 52);
    const auto *qsp_53 = buffer.data(qsp + 53);
    const auto *qsp_54 = buffer.data(qsp + 54);
    const auto *qsp_55 = buffer.data(qsp + 55);
    const auto *qsp_56 = buffer.data(qsp + 56);
    const auto *qsp_58 = buffer.data(qsp + 58);
    const auto *qsp_59 = buffer.data(qsp + 59);
    const auto *qsp_60 = buffer.data(qsp + 60);
    const auto *qsp_61 = buffer.data(qsp + 61);
    const auto *qsp_62 = buffer.data(qsp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, osp_0, qss0_0, \
                         qss1_0, qsp_0, qsp_1, qsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * osp_0[k]
                 + f_1 * qss0_0[k]
                 - f_2 * qss1_0[k]
                 + f_3 * pc_x[k] * qsp_0[k];

        t_1[k] = f_3 * pc_y[k] * qsp_0[k];

        t_2[k] = f_3 * pc_z[k] * qsp_0[k];

        t_3[k] = f_1 * qss0_0[k]
                 - f_2 * qss1_0[k]
                 + f_3 * pc_y[k] * qsp_1[k];

        t_4[k] = f_3 * pc_y[k] * qsp_2[k];

        t_5[k] = f_1 * qss0_0[k]
                 - f_2 * qss1_0[k]
                 + f_3 * pc_z[k] * qsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, osd0_0, osp_1, osp_4, \
                         osd1_0, qss0_1, qss1_1, qsp_3, qsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * osd0_0[k]
                 - f_4 * pc_y[k] * osd1_0[k];

        t_7[k] = f_5 * osp_4[k]
                 + f_3 * pc_x[k] * qsp_4[k];

        t_8[k] = f_3 * pc_z[k] * qsp_3[k];

        t_9[k] = f_6 * osp_1[k]
                 + f_1 * qss0_1[k]
                 - f_2 * qss1_1[k]
                 + f_3 * pc_y[k] * qsp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, osd0_0, osd0_5, \
                         osd1_0, osd1_5, qsp_4, qsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * qsp_4[k];

        t_11[k] = pa_y[k] * osd0_5[k]
                  - f_4 * pc_y[k] * osd1_5[k];

        t_12[k] = pa_z[k] * osd0_0[k]
                  - f_4 * pc_z[k] * osd1_0[k];

        t_13[k] = f_3 * pc_y[k] * qsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, osd0_3, osp_2, osp_8, \
                         osd1_3, qss0_2, qss1_2, qsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * osp_8[k]
                  + f_3 * pc_x[k] * qsp_8[k];

        t_15[k] = pa_z[k] * osd0_3[k]
                  - f_4 * pc_z[k] * osd1_3[k];

        t_16[k] = f_3 * pc_y[k] * qsp_8[k];

        t_17[k] = f_6 * osp_2[k]
                  + f_1 * qss0_2[k]
                  - f_2 * qss1_2[k]
                  + f_3 * pc_z[k] * qsp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, osp_4, osp_9, osp_10, \
                         qss0_3, qss1_3, qsp_9, qsp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * osp_9[k]
                  + f_1 * qss0_3[k]
                  - f_2 * qss1_3[k]
                  + f_3 * pc_x[k] * qsp_9[k];

        t_19[k] = f_7 * osp_10[k]
                  + f_3 * pc_x[k] * qsp_10[k];

        t_20[k] = f_3 * pc_z[k] * qsp_9[k];

        t_21[k] = f_8 * osp_4[k]
                  + f_1 * qss0_3[k]
                  - f_2 * qss1_3[k]
                  + f_3 * pc_y[k] * qsp_10[k];

        t_22[k] = f_3 * pc_z[k] * qsp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, osd0_12, osp_13, osd1_12, \
                         qss0_3, qss1_3, qsp_11, qsp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * qss0_3[k]
                  - f_2 * qss1_3[k]
                  + f_3 * pc_z[k] * qsp_11[k];

        t_24[k] = pa_y[k] * osd0_12[k]
                  - f_4 * pc_y[k] * osd1_12[k];

        t_25[k] = f_7 * osp_13[k]
                  + f_3 * pc_x[k] * qsp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, osd0_9, \
                         osd0_17, osp_8, osp_14, osd1_9, osd1_17, \
                         qsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * osp_14[k]
                  + f_3 * pc_x[k] * qsp_14[k];

        t_27[k] = pa_z[k] * osd0_9[k]
                  - f_4 * pc_z[k] * osd1_9[k];

        t_28[k] = f_6 * osp_8[k]
                  + f_3 * pc_y[k] * qsp_14[k];

        t_29[k] = pa_y[k] * osd0_17[k]
                  - f_4 * pc_y[k] * osd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, osp_15, osp_17, qss0_5, \
                         qss1_5, qsp_15, qsp_16, qsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * osp_15[k]
                  + f_1 * qss0_5[k]
                  - f_2 * qss1_5[k]
                  + f_3 * pc_x[k] * qsp_15[k];

        t_31[k] = f_3 * pc_y[k] * qsp_15[k];

        t_32[k] = f_7 * osp_17[k]
                  + f_3 * pc_x[k] * qsp_17[k];

        t_33[k] = f_1 * qss0_5[k]
                  - f_2 * qss1_5[k]
                  + f_3 * pc_y[k] * qsp_16[k];

        t_34[k] = f_3 * pc_y[k] * qsp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, osp_8, osp_18, osp_19, qss0_5, \
                         qss0_6, qss1_5, qss1_6, qsp_17, qsp_18, \
                         qsp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * osp_8[k]
                  + f_1 * qss0_5[k]
                  - f_2 * qss1_5[k]
                  + f_3 * pc_z[k] * qsp_17[k];

        t_36[k] = f_9 * osp_18[k]
                  + f_1 * qss0_6[k]
                  - f_2 * qss1_6[k]
                  + f_3 * pc_x[k] * qsp_18[k];

        t_37[k] = f_9 * osp_19[k]
                  + f_3 * pc_x[k] * qsp_19[k];

        t_38[k] = f_3 * pc_z[k] * qsp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, osd0_18, osp_10, osd1_18, \
                         qss0_6, qss1_6, qsp_19, qsp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * osp_10[k]
                  + f_1 * qss0_6[k]
                  - f_2 * qss1_6[k]
                  + f_3 * pc_y[k] * qsp_19[k];

        t_40[k] = f_3 * pc_z[k] * qsp_19[k];

        t_41[k] = f_1 * qss0_6[k]
                  - f_2 * qss1_6[k]
                  + f_3 * pc_z[k] * qsp_20[k];

        t_42[k] = pa_z[k] * osd0_18[k]
                  - f_4 * pc_z[k] * osd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, osd0_21, osp_14, \
                         osp_22, osp_23, osd1_21, qsp_22, qsp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * osp_22[k]
                  + f_3 * pc_x[k] * qsp_22[k];

        t_44[k] = f_9 * osp_23[k]
                  + f_3 * pc_x[k] * qsp_23[k];

        t_45[k] = pa_z[k] * osd0_21[k]
                  - f_4 * pc_z[k] * osd1_21[k];

        t_46[k] = f_8 * osp_14[k]
                  + f_3 * pc_y[k] * qsp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, osd0_30, osp_11, osp_25, \
                         osd1_30, qss0_7, qss1_7, qsp_23, qsp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * osp_11[k]
                  + f_1 * qss0_7[k]
                  - f_2 * qss1_7[k]
                  + f_3 * pc_z[k] * qsp_23[k];

        t_48[k] = pa_y[k] * osd0_30[k]
                  - f_4 * pc_y[k] * osd1_30[k];

        t_49[k] = f_9 * osp_25[k]
                  + f_3 * pc_x[k] * qsp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, osd0_35, osp_16, osp_17, \
                         osp_26, osd1_35, qss0_8, qss1_8, qsp_25, \
                         qsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * osp_26[k]
                  + f_3 * pc_x[k] * qsp_26[k];

        t_51[k] = f_6 * osp_16[k]
                  + f_1 * qss0_8[k]
                  - f_2 * qss1_8[k]
                  + f_3 * pc_y[k] * qsp_25[k];

        t_52[k] = f_6 * osp_17[k]
                  + f_3 * pc_y[k] * qsp_26[k];

        t_53[k] = pa_y[k] * osd0_35[k]
                  - f_4 * pc_y[k] * osd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, osp_27, osp_29, qss0_9, \
                         qss1_9, qsp_27, qsp_28, qsp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * osp_27[k]
                  + f_1 * qss0_9[k]
                  - f_2 * qss1_9[k]
                  + f_3 * pc_x[k] * qsp_27[k];

        t_55[k] = f_3 * pc_y[k] * qsp_27[k];

        t_56[k] = f_9 * osp_29[k]
                  + f_3 * pc_x[k] * qsp_29[k];

        t_57[k] = f_1 * qss0_9[k]
                  - f_2 * qss1_9[k]
                  + f_3 * pc_y[k] * qsp_28[k];

        t_58[k] = f_3 * pc_y[k] * qsp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, osp_17, osp_30, osp_31, qss0_9, \
                         qss0_10, qss1_9, qss1_10, qsp_29, qsp_30, \
                         qsp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * osp_17[k]
                  + f_1 * qss0_9[k]
                  - f_2 * qss1_9[k]
                  + f_3 * pc_z[k] * qsp_29[k];

        t_60[k] = f_11 * osp_30[k]
                  + f_1 * qss0_10[k]
                  - f_2 * qss1_10[k]
                  + f_3 * pc_x[k] * qsp_30[k];

        t_61[k] = f_11 * osp_31[k]
                  + f_3 * pc_x[k] * qsp_31[k];

        t_62[k] = f_3 * pc_z[k] * qsp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, osd0_36, osp_19, osd1_36, \
                         qss0_10, qss1_10, qsp_31, qsp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_12 * osp_19[k]
                  + f_1 * qss0_10[k]
                  - f_2 * qss1_10[k]
                  + f_3 * pc_y[k] * qsp_31[k];

        t_64[k] = f_3 * pc_z[k] * qsp_31[k];

        t_65[k] = f_1 * qss0_10[k]
                  - f_2 * qss1_10[k]
                  + f_3 * pc_z[k] * qsp_32[k];

        t_66[k] = pa_z[k] * osd0_36[k]
                  - f_4 * pc_z[k] * osd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, osd0_39, osp_23, \
                         osp_34, osp_35, osd1_39, qsp_34, qsp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * osp_34[k]
                  + f_3 * pc_x[k] * qsp_34[k];

        t_68[k] = f_11 * osp_35[k]
                  + f_3 * pc_x[k] * qsp_35[k];

        t_69[k] = pa_z[k] * osd0_39[k]
                  - f_4 * pc_z[k] * osd1_39[k];

        t_70[k] = f_10 * osp_23[k]
                  + f_3 * pc_y[k] * qsp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, osp_20, osp_36, osp_37, qss0_11, \
                         qss0_12, qss1_11, qss1_12, qsp_35, qsp_36, \
                         qsp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * osp_20[k]
                  + f_1 * qss0_11[k]
                  - f_2 * qss1_11[k]
                  + f_3 * pc_z[k] * qsp_35[k];

        t_72[k] = f_11 * osp_36[k]
                  + f_1 * qss0_12[k]
                  - f_2 * qss1_12[k]
                  + f_3 * pc_x[k] * qsp_36[k];

        t_73[k] = f_11 * osp_37[k]
                  + f_3 * pc_x[k] * qsp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, osp_23, osp_25, osp_26, \
                         osp_38, qss0_12, qss1_12, qsp_37, qsp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * osp_38[k]
                  + f_3 * pc_x[k] * qsp_38[k];

        t_75[k] = f_8 * osp_25[k]
                  + f_1 * qss0_12[k]
                  - f_2 * qss1_12[k]
                  + f_3 * pc_y[k] * qsp_37[k];

        t_76[k] = f_8 * osp_26[k]
                  + f_3 * pc_y[k] * qsp_38[k];

        t_77[k] = f_8 * osp_23[k]
                  + f_1 * qss0_12[k]
                  - f_2 * qss1_12[k]
                  + f_3 * pc_z[k] * qsp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, osd0_54, osp_28, osp_40, \
                         osp_41, osd1_54, qss0_13, qss1_13, qsp_40, \
                         qsp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * osd0_54[k]
                  - f_4 * pc_y[k] * osd1_54[k];

        t_79[k] = f_11 * osp_40[k]
                  + f_3 * pc_x[k] * qsp_40[k];

        t_80[k] = f_11 * osp_41[k]
                  + f_3 * pc_x[k] * qsp_41[k];

        t_81[k] = f_6 * osp_28[k]
                  + f_1 * qss0_13[k]
                  - f_2 * qss1_13[k]
                  + f_3 * pc_y[k] * qsp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, osd0_59, osp_29, osp_42, \
                         osd1_59, qss0_14, qss1_14, qsp_41, qsp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * osp_29[k]
                  + f_3 * pc_y[k] * qsp_41[k];

        t_83[k] = pa_y[k] * osd0_59[k]
                  - f_4 * pc_y[k] * osd1_59[k];

        t_84[k] = f_11 * osp_42[k]
                  + f_1 * qss0_14[k]
                  - f_2 * qss1_14[k]
                  + f_3 * pc_x[k] * qsp_42[k];

        t_85[k] = f_3 * pc_y[k] * qsp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, osp_29, osp_44, qss0_14, \
                         qss1_14, qsp_43, qsp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * osp_44[k]
                  + f_3 * pc_x[k] * qsp_44[k];

        t_87[k] = f_1 * qss0_14[k]
                  - f_2 * qss1_14[k]
                  + f_3 * pc_y[k] * qsp_43[k];

        t_88[k] = f_3 * pc_y[k] * qsp_44[k];

        t_89[k] = f_12 * osp_29[k]
                  + f_1 * qss0_14[k]
                  - f_2 * qss1_14[k]
                  + f_3 * pc_z[k] * qsp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, osp_31, osp_45, \
                         osp_46, qss0_15, qss1_15, qsp_45, qsp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_13 * osp_45[k]
                  + f_1 * qss0_15[k]
                  - f_2 * qss1_15[k]
                  + f_3 * pc_x[k] * qsp_45[k];

        t_91[k] = f_13 * osp_46[k]
                  + f_3 * pc_x[k] * qsp_46[k];

        t_92[k] = f_3 * pc_z[k] * qsp_45[k];

        t_93[k] = f_14 * osp_31[k]
                  + f_1 * qss0_15[k]
                  - f_2 * qss1_15[k]
                  + f_3 * pc_y[k] * qsp_46[k];

        t_94[k] = f_3 * pc_z[k] * qsp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pc_x, pc_z, osd0_60, osp_49, osp_50, \
                         osd1_60, qss0_15, qss1_15, qsp_47, qsp_49, \
                         qsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_1 * qss0_15[k]
                  - f_2 * qss1_15[k]
                  + f_3 * pc_z[k] * qsp_47[k];

        t_96[k] = pa_z[k] * osd0_60[k]
                  - f_4 * pc_z[k] * osd1_60[k];

        t_97[k] = f_13 * osp_49[k]
                  + f_3 * pc_x[k] * qsp_49[k];

        t_98[k] = f_13 * osp_50[k]
                  + f_3 * pc_x[k] * qsp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, osd0_63, osp_32, osp_35, \
                         osd1_63, qss0_16, qss1_16, qsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * osd0_63[k]
                  - f_4 * pc_z[k] * osd1_63[k];

        t_100[k] = f_12 * osp_35[k]
                   + f_3 * pc_y[k] * qsp_50[k];

        t_101[k] = f_6 * osp_32[k]
                   + f_1 * qss0_16[k]
                   - f_2 * qss1_16[k]
                   + f_3 * pc_z[k] * qsp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, osp_37, osp_51, osp_52, \
                         osp_53, qss0_17, qss1_17, qsp_51, qsp_52, \
                         qsp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_13 * osp_51[k]
                   + f_1 * qss0_17[k]
                   - f_2 * qss1_17[k]
                   + f_3 * pc_x[k] * qsp_51[k];

        t_103[k] = f_13 * osp_52[k]
                   + f_3 * pc_x[k] * qsp_52[k];

        t_104[k] = f_13 * osp_53[k]
                   + f_3 * pc_x[k] * qsp_53[k];

        t_105[k] = f_10 * osp_37[k]
                   + f_1 * qss0_17[k]
                   - f_2 * qss1_17[k]
                   + f_3 * pc_y[k] * qsp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, osp_35, osp_38, osp_54, \
                         qss0_17, qss0_18, qss1_17, qss1_18, qsp_53, \
                         qsp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * osp_38[k]
                   + f_3 * pc_y[k] * qsp_53[k];

        t_107[k] = f_8 * osp_35[k]
                   + f_1 * qss0_17[k]
                   - f_2 * qss1_17[k]
                   + f_3 * pc_z[k] * qsp_53[k];

        t_108[k] = f_13 * osp_54[k]
                   + f_1 * qss0_18[k]
                   - f_2 * qss1_18[k]
                   + f_3 * pc_x[k] * qsp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, osp_40, osp_41, osp_55, \
                         osp_56, qss0_18, qss1_18, qsp_55, qsp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_13 * osp_55[k]
                   + f_3 * pc_x[k] * qsp_55[k];

        t_110[k] = f_13 * osp_56[k]
                   + f_3 * pc_x[k] * qsp_56[k];

        t_111[k] = f_8 * osp_40[k]
                   + f_1 * qss0_18[k]
                   - f_2 * qss1_18[k]
                   + f_3 * pc_y[k] * qsp_55[k];

        t_112[k] = f_8 * osp_41[k]
                   + f_3 * pc_y[k] * qsp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, osd0_84, osp_38, osp_58, \
                         osd1_84, qss0_18, qss1_18, qsp_56, qsp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * osp_38[k]
                   + f_1 * qss0_18[k]
                   - f_2 * qss1_18[k]
                   + f_3 * pc_z[k] * qsp_56[k];

        t_114[k] = pa_y[k] * osd0_84[k]
                   - f_4 * pc_y[k] * osd1_84[k];

        t_115[k] = f_13 * osp_58[k]
                   + f_3 * pc_x[k] * qsp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, osd0_89, osp_43, \
                         osp_44, osp_59, osd1_89, qss0_19, qss1_19, qsp_58, \
                         qsp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_13 * osp_59[k]
                   + f_3 * pc_x[k] * qsp_59[k];

        t_117[k] = f_6 * osp_43[k]
                   + f_1 * qss0_19[k]
                   - f_2 * qss1_19[k]
                   + f_3 * pc_y[k] * qsp_58[k];

        t_118[k] = f_6 * osp_44[k]
                   + f_3 * pc_y[k] * qsp_59[k];

        t_119[k] = pa_y[k] * osd0_89[k]
                   - f_4 * pc_y[k] * osd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, osp_60, osp_62, \
                         qss0_20, qss1_20, qsp_60, qsp_61, qsp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_13 * osp_60[k]
                   + f_1 * qss0_20[k]
                   - f_2 * qss1_20[k]
                   + f_3 * pc_x[k] * qsp_60[k];

        t_121[k] = f_3 * pc_y[k] * qsp_60[k];

        t_122[k] = f_13 * osp_62[k]
                   + f_3 * pc_x[k] * qsp_62[k];

        t_123[k] = f_1 * qss0_20[k]
                   - f_2 * qss1_20[k]
                   + f_3 * pc_y[k] * qsp_61[k];

        t_124[k] = f_3 * pc_y[k] * qsp_62[k];
    }
}

static auto
compute_prim_qsd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osd0,
                                                          const size_t osp, const size_t osd1,
                                                          const size_t qss0, const size_t qss1,
                                                          const size_t qsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osd0_90 = buffer.data(osd0 + 90);
    const auto *osd0_93 = buffer.data(osd0 + 93);
    const auto *osd0_120 = buffer.data(osd0 + 120);
    const auto *osd0_125 = buffer.data(osd0 + 125);
    const auto *osd0_126 = buffer.data(osd0 + 126);
    const auto *osd0_129 = buffer.data(osd0 + 129);
    const auto *osd0_162 = buffer.data(osd0 + 162);
    const auto *osd0_167 = buffer.data(osd0 + 167);
    const auto *osd0_168 = buffer.data(osd0 + 168);
    const auto *osd0_171 = buffer.data(osd0 + 171);

    const auto *osp_44 = buffer.data(osp + 44);
    const auto *osp_46 = buffer.data(osp + 46);
    const auto *osp_47 = buffer.data(osp + 47);
    const auto *osp_50 = buffer.data(osp + 50);
    const auto *osp_52 = buffer.data(osp + 52);
    const auto *osp_53 = buffer.data(osp + 53);
    const auto *osp_55 = buffer.data(osp + 55);
    const auto *osp_56 = buffer.data(osp + 56);
    const auto *osp_58 = buffer.data(osp + 58);
    const auto *osp_59 = buffer.data(osp + 59);
    const auto *osp_61 = buffer.data(osp + 61);
    const auto *osp_62 = buffer.data(osp + 62);
    const auto *osp_63 = buffer.data(osp + 63);
    const auto *osp_64 = buffer.data(osp + 64);
    const auto *osp_65 = buffer.data(osp + 65);
    const auto *osp_67 = buffer.data(osp + 67);
    const auto *osp_68 = buffer.data(osp + 68);
    const auto *osp_69 = buffer.data(osp + 69);
    const auto *osp_70 = buffer.data(osp + 70);
    const auto *osp_71 = buffer.data(osp + 71);
    const auto *osp_72 = buffer.data(osp + 72);
    const auto *osp_73 = buffer.data(osp + 73);
    const auto *osp_74 = buffer.data(osp + 74);
    const auto *osp_75 = buffer.data(osp + 75);
    const auto *osp_76 = buffer.data(osp + 76);
    const auto *osp_77 = buffer.data(osp + 77);
    const auto *osp_79 = buffer.data(osp + 79);
    const auto *osp_80 = buffer.data(osp + 80);
    const auto *osp_81 = buffer.data(osp + 81);
    const auto *osp_82 = buffer.data(osp + 82);
    const auto *osp_83 = buffer.data(osp + 83);
    const auto *osp_84 = buffer.data(osp + 84);
    const auto *osp_85 = buffer.data(osp + 85);
    const auto *osp_86 = buffer.data(osp + 86);
    const auto *osp_88 = buffer.data(osp + 88);
    const auto *osp_89 = buffer.data(osp + 89);
    const auto *osp_90 = buffer.data(osp + 90);
    const auto *osp_91 = buffer.data(osp + 91);
    const auto *osp_92 = buffer.data(osp + 92);
    const auto *osp_93 = buffer.data(osp + 93);
    const auto *osp_94 = buffer.data(osp + 94);
    const auto *osp_95 = buffer.data(osp + 95);
    const auto *osp_96 = buffer.data(osp + 96);
    const auto *osp_97 = buffer.data(osp + 97);
    const auto *osp_98 = buffer.data(osp + 98);
    const auto *osp_99 = buffer.data(osp + 99);
    const auto *osp_100 = buffer.data(osp + 100);
    const auto *osp_101 = buffer.data(osp + 101);
    const auto *osp_103 = buffer.data(osp + 103);
    const auto *osp_104 = buffer.data(osp + 104);
    const auto *osp_105 = buffer.data(osp + 105);
    const auto *osp_107 = buffer.data(osp + 107);
    const auto *osp_108 = buffer.data(osp + 108);
    const auto *osp_109 = buffer.data(osp + 109);
    const auto *osp_112 = buffer.data(osp + 112);
    const auto *osp_113 = buffer.data(osp + 113);
    const auto *osp_114 = buffer.data(osp + 114);
    const auto *osp_115 = buffer.data(osp + 115);
    const auto *osp_116 = buffer.data(osp + 116);
    const auto *osp_117 = buffer.data(osp + 117);
    const auto *osp_118 = buffer.data(osp + 118);
    const auto *osp_119 = buffer.data(osp + 119);
    const auto *osp_120 = buffer.data(osp + 120);
    const auto *osp_121 = buffer.data(osp + 121);

    const auto *osd1_90 = buffer.data(osd1 + 90);
    const auto *osd1_93 = buffer.data(osd1 + 93);
    const auto *osd1_120 = buffer.data(osd1 + 120);
    const auto *osd1_125 = buffer.data(osd1 + 125);
    const auto *osd1_126 = buffer.data(osd1 + 126);
    const auto *osd1_129 = buffer.data(osd1 + 129);
    const auto *osd1_162 = buffer.data(osd1 + 162);
    const auto *osd1_167 = buffer.data(osd1 + 167);
    const auto *osd1_168 = buffer.data(osd1 + 168);
    const auto *osd1_171 = buffer.data(osd1 + 171);

    const auto *qss0_20 = buffer.data(qss0 + 20);
    const auto *qss0_21 = buffer.data(qss0 + 21);
    const auto *qss0_22 = buffer.data(qss0 + 22);
    const auto *qss0_23 = buffer.data(qss0 + 23);
    const auto *qss0_24 = buffer.data(qss0 + 24);
    const auto *qss0_25 = buffer.data(qss0 + 25);
    const auto *qss0_26 = buffer.data(qss0 + 26);
    const auto *qss0_27 = buffer.data(qss0 + 27);
    const auto *qss0_28 = buffer.data(qss0 + 28);
    const auto *qss0_29 = buffer.data(qss0 + 29);
    const auto *qss0_30 = buffer.data(qss0 + 30);
    const auto *qss0_31 = buffer.data(qss0 + 31);
    const auto *qss0_32 = buffer.data(qss0 + 32);
    const auto *qss0_33 = buffer.data(qss0 + 33);
    const auto *qss0_34 = buffer.data(qss0 + 34);
    const auto *qss0_35 = buffer.data(qss0 + 35);
    const auto *qss0_36 = buffer.data(qss0 + 36);
    const auto *qss0_37 = buffer.data(qss0 + 37);
    const auto *qss0_38 = buffer.data(qss0 + 38);
    const auto *qss0_39 = buffer.data(qss0 + 39);
    const auto *qss0_40 = buffer.data(qss0 + 40);

    const auto *qss1_20 = buffer.data(qss1 + 20);
    const auto *qss1_21 = buffer.data(qss1 + 21);
    const auto *qss1_22 = buffer.data(qss1 + 22);
    const auto *qss1_23 = buffer.data(qss1 + 23);
    const auto *qss1_24 = buffer.data(qss1 + 24);
    const auto *qss1_25 = buffer.data(qss1 + 25);
    const auto *qss1_26 = buffer.data(qss1 + 26);
    const auto *qss1_27 = buffer.data(qss1 + 27);
    const auto *qss1_28 = buffer.data(qss1 + 28);
    const auto *qss1_29 = buffer.data(qss1 + 29);
    const auto *qss1_30 = buffer.data(qss1 + 30);
    const auto *qss1_31 = buffer.data(qss1 + 31);
    const auto *qss1_32 = buffer.data(qss1 + 32);
    const auto *qss1_33 = buffer.data(qss1 + 33);
    const auto *qss1_34 = buffer.data(qss1 + 34);
    const auto *qss1_35 = buffer.data(qss1 + 35);
    const auto *qss1_36 = buffer.data(qss1 + 36);
    const auto *qss1_37 = buffer.data(qss1 + 37);
    const auto *qss1_38 = buffer.data(qss1 + 38);
    const auto *qss1_39 = buffer.data(qss1 + 39);
    const auto *qss1_40 = buffer.data(qss1 + 40);

    const auto *qsp_62 = buffer.data(qsp + 62);
    const auto *qsp_63 = buffer.data(qsp + 63);
    const auto *qsp_64 = buffer.data(qsp + 64);
    const auto *qsp_65 = buffer.data(qsp + 65);
    const auto *qsp_67 = buffer.data(qsp + 67);
    const auto *qsp_68 = buffer.data(qsp + 68);
    const auto *qsp_69 = buffer.data(qsp + 69);
    const auto *qsp_70 = buffer.data(qsp + 70);
    const auto *qsp_71 = buffer.data(qsp + 71);
    const auto *qsp_72 = buffer.data(qsp + 72);
    const auto *qsp_73 = buffer.data(qsp + 73);
    const auto *qsp_74 = buffer.data(qsp + 74);
    const auto *qsp_75 = buffer.data(qsp + 75);
    const auto *qsp_76 = buffer.data(qsp + 76);
    const auto *qsp_77 = buffer.data(qsp + 77);
    const auto *qsp_79 = buffer.data(qsp + 79);
    const auto *qsp_80 = buffer.data(qsp + 80);
    const auto *qsp_81 = buffer.data(qsp + 81);
    const auto *qsp_82 = buffer.data(qsp + 82);
    const auto *qsp_83 = buffer.data(qsp + 83);
    const auto *qsp_84 = buffer.data(qsp + 84);
    const auto *qsp_85 = buffer.data(qsp + 85);
    const auto *qsp_86 = buffer.data(qsp + 86);
    const auto *qsp_88 = buffer.data(qsp + 88);
    const auto *qsp_89 = buffer.data(qsp + 89);
    const auto *qsp_90 = buffer.data(qsp + 90);
    const auto *qsp_91 = buffer.data(qsp + 91);
    const auto *qsp_92 = buffer.data(qsp + 92);
    const auto *qsp_93 = buffer.data(qsp + 93);
    const auto *qsp_94 = buffer.data(qsp + 94);
    const auto *qsp_95 = buffer.data(qsp + 95);
    const auto *qsp_96 = buffer.data(qsp + 96);
    const auto *qsp_97 = buffer.data(qsp + 97);
    const auto *qsp_98 = buffer.data(qsp + 98);
    const auto *qsp_99 = buffer.data(qsp + 99);
    const auto *qsp_100 = buffer.data(qsp + 100);
    const auto *qsp_101 = buffer.data(qsp + 101);
    const auto *qsp_103 = buffer.data(qsp + 103);
    const auto *qsp_104 = buffer.data(qsp + 104);
    const auto *qsp_105 = buffer.data(qsp + 105);
    const auto *qsp_106 = buffer.data(qsp + 106);
    const auto *qsp_107 = buffer.data(qsp + 107);
    const auto *qsp_108 = buffer.data(qsp + 108);
    const auto *qsp_109 = buffer.data(qsp + 109);
    const auto *qsp_110 = buffer.data(qsp + 110);
    const auto *qsp_112 = buffer.data(qsp + 112);
    const auto *qsp_113 = buffer.data(qsp + 113);
    const auto *qsp_114 = buffer.data(qsp + 114);
    const auto *qsp_115 = buffer.data(qsp + 115);
    const auto *qsp_116 = buffer.data(qsp + 116);
    const auto *qsp_117 = buffer.data(qsp + 117);
    const auto *qsp_118 = buffer.data(qsp + 118);
    const auto *qsp_119 = buffer.data(qsp + 119);
    const auto *qsp_120 = buffer.data(qsp + 120);
    const auto *qsp_121 = buffer.data(qsp + 121);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_z, osp_44, osp_63, osp_64, \
                         qss0_20, qss0_21, qss1_20, qss1_21, qsp_62, qsp_63, \
                         qsp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_14 * osp_44[k]
                   + f_1 * qss0_20[k]
                   - f_2 * qss1_20[k]
                   + f_3 * pc_z[k] * qsp_62[k];

        t_126[k] = f_15 * osp_63[k]
                   + f_1 * qss0_21[k]
                   - f_2 * qss1_21[k]
                   + f_3 * pc_x[k] * qsp_63[k];

        t_127[k] = f_15 * osp_64[k]
                   + f_3 * pc_x[k] * qsp_64[k];

        t_128[k] = f_3 * pc_z[k] * qsp_63[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_z, pc_y, pc_z, osd0_90, osp_46, \
                         osd1_90, qss0_21, qss1_21, qsp_64, qsp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_15 * osp_46[k]
                   + f_1 * qss0_21[k]
                   - f_2 * qss1_21[k]
                   + f_3 * pc_y[k] * qsp_64[k];

        t_130[k] = f_3 * pc_z[k] * qsp_64[k];

        t_131[k] = f_1 * qss0_21[k]
                   - f_2 * qss1_21[k]
                   + f_3 * pc_z[k] * qsp_65[k];

        t_132[k] = pa_z[k] * osd0_90[k]
                   - f_4 * pc_z[k] * osd1_90[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_z, pc_x, pc_y, pc_z, osd0_93, osp_50, \
                         osp_67, osp_68, osd1_93, qsp_67, qsp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_15 * osp_67[k]
                   + f_3 * pc_x[k] * qsp_67[k];

        t_134[k] = f_15 * osp_68[k]
                   + f_3 * pc_x[k] * qsp_68[k];

        t_135[k] = pa_z[k] * osd0_93[k]
                   - f_4 * pc_z[k] * osd1_93[k];

        t_136[k] = f_14 * osp_50[k]
                   + f_3 * pc_y[k] * qsp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_x, pc_z, osp_47, osp_69, osp_70, qss0_22, \
                         qss0_23, qss1_22, qss1_23, qsp_68, qsp_69, \
                         qsp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * osp_47[k]
                   + f_1 * qss0_22[k]
                   - f_2 * qss1_22[k]
                   + f_3 * pc_z[k] * qsp_68[k];

        t_138[k] = f_15 * osp_69[k]
                   + f_1 * qss0_23[k]
                   - f_2 * qss1_23[k]
                   + f_3 * pc_x[k] * qsp_69[k];

        t_139[k] = f_15 * osp_70[k]
                   + f_3 * pc_x[k] * qsp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, osp_50, osp_52, osp_53, \
                         osp_71, qss0_23, qss1_23, qsp_70, qsp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * osp_71[k]
                   + f_3 * pc_x[k] * qsp_71[k];

        t_141[k] = f_12 * osp_52[k]
                   + f_1 * qss0_23[k]
                   - f_2 * qss1_23[k]
                   + f_3 * pc_y[k] * qsp_70[k];

        t_142[k] = f_12 * osp_53[k]
                   + f_3 * pc_y[k] * qsp_71[k];

        t_143[k] = f_8 * osp_50[k]
                   + f_1 * qss0_23[k]
                   - f_2 * qss1_23[k]
                   + f_3 * pc_z[k] * qsp_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pc_x, pc_y, osp_55, osp_72, osp_73, \
                         osp_74, qss0_24, qss1_24, qsp_72, qsp_73, \
                         qsp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_15 * osp_72[k]
                   + f_1 * qss0_24[k]
                   - f_2 * qss1_24[k]
                   + f_3 * pc_x[k] * qsp_72[k];

        t_145[k] = f_15 * osp_73[k]
                   + f_3 * pc_x[k] * qsp_73[k];

        t_146[k] = f_15 * osp_74[k]
                   + f_3 * pc_x[k] * qsp_74[k];

        t_147[k] = f_10 * osp_55[k]
                   + f_1 * qss0_24[k]
                   - f_2 * qss1_24[k]
                   + f_3 * pc_y[k] * qsp_73[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pc_x, pc_y, pc_z, osp_53, osp_56, osp_75, \
                         qss0_24, qss0_25, qss1_24, qss1_25, qsp_74, \
                         qsp_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * osp_56[k]
                   + f_3 * pc_y[k] * qsp_74[k];

        t_149[k] = f_10 * osp_53[k]
                   + f_1 * qss0_24[k]
                   - f_2 * qss1_24[k]
                   + f_3 * pc_z[k] * qsp_74[k];

        t_150[k] = f_15 * osp_75[k]
                   + f_1 * qss0_25[k]
                   - f_2 * qss1_25[k]
                   + f_3 * pc_x[k] * qsp_75[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, osp_58, osp_59, osp_76, \
                         osp_77, qss0_25, qss1_25, qsp_76, qsp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_15 * osp_76[k]
                   + f_3 * pc_x[k] * qsp_76[k];

        t_152[k] = f_15 * osp_77[k]
                   + f_3 * pc_x[k] * qsp_77[k];

        t_153[k] = f_8 * osp_58[k]
                   + f_1 * qss0_25[k]
                   - f_2 * qss1_25[k]
                   + f_3 * pc_y[k] * qsp_76[k];

        t_154[k] = f_8 * osp_59[k]
                   + f_3 * pc_y[k] * qsp_77[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pa_y, pc_x, pc_y, pc_z, osd0_120, osp_56, \
                         osp_79, osd1_120, qss0_25, qss1_25, qsp_77, \
                         qsp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * osp_56[k]
                   + f_1 * qss0_25[k]
                   - f_2 * qss1_25[k]
                   + f_3 * pc_z[k] * qsp_77[k];

        t_156[k] = pa_y[k] * osd0_120[k]
                   - f_4 * pc_y[k] * osd1_120[k];

        t_157[k] = f_15 * osp_79[k]
                   + f_3 * pc_x[k] * qsp_79[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_y, pc_x, pc_y, osd0_125, osp_61, \
                         osp_62, osp_80, osd1_125, qss0_26, qss1_26, qsp_79, \
                         qsp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_15 * osp_80[k]
                   + f_3 * pc_x[k] * qsp_80[k];

        t_159[k] = f_6 * osp_61[k]
                   + f_1 * qss0_26[k]
                   - f_2 * qss1_26[k]
                   + f_3 * pc_y[k] * qsp_79[k];

        t_160[k] = f_6 * osp_62[k]
                   + f_3 * pc_y[k] * qsp_80[k];

        t_161[k] = pa_y[k] * osd0_125[k]
                   - f_4 * pc_y[k] * osd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pc_x, pc_y, osp_81, osp_83, \
                         qss0_27, qss1_27, qsp_81, qsp_82, qsp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_15 * osp_81[k]
                   + f_1 * qss0_27[k]
                   - f_2 * qss1_27[k]
                   + f_3 * pc_x[k] * qsp_81[k];

        t_163[k] = f_3 * pc_y[k] * qsp_81[k];

        t_164[k] = f_15 * osp_83[k]
                   + f_3 * pc_x[k] * qsp_83[k];

        t_165[k] = f_1 * qss0_27[k]
                   - f_2 * qss1_27[k]
                   + f_3 * pc_y[k] * qsp_82[k];

        t_166[k] = f_3 * pc_y[k] * qsp_83[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pc_x, pc_z, osp_62, osp_84, osp_85, \
                         qss0_27, qss0_28, qss1_27, qss1_28, qsp_83, qsp_84, \
                         qsp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_15 * osp_62[k]
                   + f_1 * qss0_27[k]
                   - f_2 * qss1_27[k]
                   + f_3 * pc_z[k] * qsp_83[k];

        t_168[k] = f_14 * osp_84[k]
                   + f_1 * qss0_28[k]
                   - f_2 * qss1_28[k]
                   + f_3 * pc_x[k] * qsp_84[k];

        t_169[k] = f_14 * osp_85[k]
                   + f_3 * pc_x[k] * qsp_85[k];

        t_170[k] = f_3 * pc_z[k] * qsp_84[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_z, pc_y, pc_z, osd0_126, osp_64, \
                         osd1_126, qss0_28, qss1_28, qsp_85, qsp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_13 * osp_64[k]
                   + f_1 * qss0_28[k]
                   - f_2 * qss1_28[k]
                   + f_3 * pc_y[k] * qsp_85[k];

        t_172[k] = f_3 * pc_z[k] * qsp_85[k];

        t_173[k] = f_1 * qss0_28[k]
                   - f_2 * qss1_28[k]
                   + f_3 * pc_z[k] * qsp_86[k];

        t_174[k] = pa_z[k] * osd0_126[k]
                   - f_4 * pc_z[k] * osd1_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_z, pc_x, pc_y, pc_z, osd0_129, osp_68, \
                         osp_88, osp_89, osd1_129, qsp_88, qsp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_14 * osp_88[k]
                   + f_3 * pc_x[k] * qsp_88[k];

        t_176[k] = f_14 * osp_89[k]
                   + f_3 * pc_x[k] * qsp_89[k];

        t_177[k] = pa_z[k] * osd0_129[k]
                   - f_4 * pc_z[k] * osd1_129[k];

        t_178[k] = f_15 * osp_68[k]
                   + f_3 * pc_y[k] * qsp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_x, pc_z, osp_65, osp_90, osp_91, qss0_29, \
                         qss0_30, qss1_29, qss1_30, qsp_89, qsp_90, \
                         qsp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * osp_65[k]
                   + f_1 * qss0_29[k]
                   - f_2 * qss1_29[k]
                   + f_3 * pc_z[k] * qsp_89[k];

        t_180[k] = f_14 * osp_90[k]
                   + f_1 * qss0_30[k]
                   - f_2 * qss1_30[k]
                   + f_3 * pc_x[k] * qsp_90[k];

        t_181[k] = f_14 * osp_91[k]
                   + f_3 * pc_x[k] * qsp_91[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, osp_68, osp_70, osp_71, \
                         osp_92, qss0_30, qss1_30, qsp_91, qsp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_14 * osp_92[k]
                   + f_3 * pc_x[k] * qsp_92[k];

        t_183[k] = f_14 * osp_70[k]
                   + f_1 * qss0_30[k]
                   - f_2 * qss1_30[k]
                   + f_3 * pc_y[k] * qsp_91[k];

        t_184[k] = f_14 * osp_71[k]
                   + f_3 * pc_y[k] * qsp_92[k];

        t_185[k] = f_8 * osp_68[k]
                   + f_1 * qss0_30[k]
                   - f_2 * qss1_30[k]
                   + f_3 * pc_z[k] * qsp_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, pc_y, osp_73, osp_93, osp_94, \
                         osp_95, qss0_31, qss1_31, qsp_93, qsp_94, \
                         qsp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_14 * osp_93[k]
                   + f_1 * qss0_31[k]
                   - f_2 * qss1_31[k]
                   + f_3 * pc_x[k] * qsp_93[k];

        t_187[k] = f_14 * osp_94[k]
                   + f_3 * pc_x[k] * qsp_94[k];

        t_188[k] = f_14 * osp_95[k]
                   + f_3 * pc_x[k] * qsp_95[k];

        t_189[k] = f_12 * osp_73[k]
                   + f_1 * qss0_31[k]
                   - f_2 * qss1_31[k]
                   + f_3 * pc_y[k] * qsp_94[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_y, pc_z, osp_71, osp_74, osp_96, \
                         qss0_31, qss0_32, qss1_31, qss1_32, qsp_95, \
                         qsp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_12 * osp_74[k]
                   + f_3 * pc_y[k] * qsp_95[k];

        t_191[k] = f_10 * osp_71[k]
                   + f_1 * qss0_31[k]
                   - f_2 * qss1_31[k]
                   + f_3 * pc_z[k] * qsp_95[k];

        t_192[k] = f_14 * osp_96[k]
                   + f_1 * qss0_32[k]
                   - f_2 * qss1_32[k]
                   + f_3 * pc_x[k] * qsp_96[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, osp_76, osp_77, osp_97, \
                         osp_98, qss0_32, qss1_32, qsp_97, qsp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_14 * osp_97[k]
                   + f_3 * pc_x[k] * qsp_97[k];

        t_194[k] = f_14 * osp_98[k]
                   + f_3 * pc_x[k] * qsp_98[k];

        t_195[k] = f_10 * osp_76[k]
                   + f_1 * qss0_32[k]
                   - f_2 * qss1_32[k]
                   + f_3 * pc_y[k] * qsp_97[k];

        t_196[k] = f_10 * osp_77[k]
                   + f_3 * pc_y[k] * qsp_98[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pc_x, pc_z, osp_74, osp_99, osp_100, qss0_32, \
                         qss0_33, qss1_32, qss1_33, qsp_98, qsp_99, \
                         qsp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * osp_74[k]
                   + f_1 * qss0_32[k]
                   - f_2 * qss1_32[k]
                   + f_3 * pc_z[k] * qsp_98[k];

        t_198[k] = f_14 * osp_99[k]
                   + f_1 * qss0_33[k]
                   - f_2 * qss1_33[k]
                   + f_3 * pc_x[k] * qsp_99[k];

        t_199[k] = f_14 * osp_100[k]
                   + f_3 * pc_x[k] * qsp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, osp_77, osp_79, osp_80, \
                         osp_101, qss0_33, qss1_33, qsp_100, qsp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * osp_101[k]
                   + f_3 * pc_x[k] * qsp_101[k];

        t_201[k] = f_8 * osp_79[k]
                   + f_1 * qss0_33[k]
                   - f_2 * qss1_33[k]
                   + f_3 * pc_y[k] * qsp_100[k];

        t_202[k] = f_8 * osp_80[k]
                   + f_3 * pc_y[k] * qsp_101[k];

        t_203[k] = f_14 * osp_77[k]
                   + f_1 * qss0_33[k]
                   - f_2 * qss1_33[k]
                   + f_3 * pc_z[k] * qsp_101[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_x, pc_y, osd0_162, osp_82, \
                         osp_103, osp_104, osd1_162, qss0_34, qss1_34, qsp_103, \
                         qsp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * osd0_162[k]
                   - f_4 * pc_y[k] * osd1_162[k];

        t_205[k] = f_14 * osp_103[k]
                   + f_3 * pc_x[k] * qsp_103[k];

        t_206[k] = f_14 * osp_104[k]
                   + f_3 * pc_x[k] * qsp_104[k];

        t_207[k] = f_6 * osp_82[k]
                   + f_1 * qss0_34[k]
                   - f_2 * qss1_34[k]
                   + f_3 * pc_y[k] * qsp_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_y, pc_x, pc_y, osd0_167, osp_83, \
                         osp_105, osd1_167, qss0_35, qss1_35, qsp_104, \
                         qsp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_6 * osp_83[k]
                   + f_3 * pc_y[k] * qsp_104[k];

        t_209[k] = pa_y[k] * osd0_167[k]
                   - f_4 * pc_y[k] * osd1_167[k];

        t_210[k] = f_14 * osp_105[k]
                   + f_1 * qss0_35[k]
                   - f_2 * qss1_35[k]
                   + f_3 * pc_x[k] * qsp_105[k];

        t_211[k] = f_3 * pc_y[k] * qsp_105[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, osp_83, osp_107, \
                         qss0_35, qss1_35, qsp_106, qsp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_14 * osp_107[k]
                   + f_3 * pc_x[k] * qsp_107[k];

        t_213[k] = f_1 * qss0_35[k]
                   - f_2 * qss1_35[k]
                   + f_3 * pc_y[k] * qsp_106[k];

        t_214[k] = f_3 * pc_y[k] * qsp_107[k];

        t_215[k] = f_13 * osp_83[k]
                   + f_1 * qss0_35[k]
                   - f_2 * qss1_35[k]
                   + f_3 * pc_z[k] * qsp_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, osp_85, osp_108, \
                         osp_109, qss0_36, qss1_36, qsp_108, qsp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_12 * osp_108[k]
                   + f_1 * qss0_36[k]
                   - f_2 * qss1_36[k]
                   + f_3 * pc_x[k] * qsp_108[k];

        t_217[k] = f_12 * osp_109[k]
                   + f_3 * pc_x[k] * qsp_109[k];

        t_218[k] = f_3 * pc_z[k] * qsp_108[k];

        t_219[k] = f_11 * osp_85[k]
                   + f_1 * qss0_36[k]
                   - f_2 * qss1_36[k]
                   + f_3 * pc_y[k] * qsp_109[k];

        t_220[k] = f_3 * pc_z[k] * qsp_109[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_z, pc_x, pc_z, osd0_168, osp_112, \
                         osp_113, osd1_168, qss0_36, qss1_36, qsp_110, qsp_112, \
                         qsp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_1 * qss0_36[k]
                   - f_2 * qss1_36[k]
                   + f_3 * pc_z[k] * qsp_110[k];

        t_222[k] = pa_z[k] * osd0_168[k]
                   - f_4 * pc_z[k] * osd1_168[k];

        t_223[k] = f_12 * osp_112[k]
                   + f_3 * pc_x[k] * qsp_112[k];

        t_224[k] = f_12 * osp_113[k]
                   + f_3 * pc_x[k] * qsp_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_z, pc_y, pc_z, osd0_171, osp_86, osp_89, \
                         osd1_171, qss0_37, qss1_37, qsp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_z[k] * osd0_171[k]
                   - f_4 * pc_z[k] * osd1_171[k];

        t_226[k] = f_13 * osp_89[k]
                   + f_3 * pc_y[k] * qsp_113[k];

        t_227[k] = f_6 * osp_86[k]
                   + f_1 * qss0_37[k]
                   - f_2 * qss1_37[k]
                   + f_3 * pc_z[k] * qsp_113[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, osp_91, osp_114, osp_115, \
                         osp_116, qss0_38, qss1_38, qsp_114, qsp_115, \
                         qsp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_12 * osp_114[k]
                   + f_1 * qss0_38[k]
                   - f_2 * qss1_38[k]
                   + f_3 * pc_x[k] * qsp_114[k];

        t_229[k] = f_12 * osp_115[k]
                   + f_3 * pc_x[k] * qsp_115[k];

        t_230[k] = f_12 * osp_116[k]
                   + f_3 * pc_x[k] * qsp_116[k];

        t_231[k] = f_15 * osp_91[k]
                   + f_1 * qss0_38[k]
                   - f_2 * qss1_38[k]
                   + f_3 * pc_y[k] * qsp_115[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_y, pc_z, osp_89, osp_92, osp_117, \
                         qss0_38, qss0_39, qss1_38, qss1_39, qsp_116, \
                         qsp_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_15 * osp_92[k]
                   + f_3 * pc_y[k] * qsp_116[k];

        t_233[k] = f_8 * osp_89[k]
                   + f_1 * qss0_38[k]
                   - f_2 * qss1_38[k]
                   + f_3 * pc_z[k] * qsp_116[k];

        t_234[k] = f_12 * osp_117[k]
                   + f_1 * qss0_39[k]
                   - f_2 * qss1_39[k]
                   + f_3 * pc_x[k] * qsp_117[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, osp_94, osp_95, osp_118, \
                         osp_119, qss0_39, qss1_39, qsp_118, qsp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_12 * osp_118[k]
                   + f_3 * pc_x[k] * qsp_118[k];

        t_236[k] = f_12 * osp_119[k]
                   + f_3 * pc_x[k] * qsp_119[k];

        t_237[k] = f_14 * osp_94[k]
                   + f_1 * qss0_39[k]
                   - f_2 * qss1_39[k]
                   + f_3 * pc_y[k] * qsp_118[k];

        t_238[k] = f_14 * osp_95[k]
                   + f_3 * pc_y[k] * qsp_119[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, osp_92, osp_120, osp_121, qss0_39, \
                         qss0_40, qss1_39, qss1_40, qsp_119, qsp_120, \
                         qsp_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * osp_92[k]
                   + f_1 * qss0_39[k]
                   - f_2 * qss1_39[k]
                   + f_3 * pc_z[k] * qsp_119[k];

        t_240[k] = f_12 * osp_120[k]
                   + f_1 * qss0_40[k]
                   - f_2 * qss1_40[k]
                   + f_3 * pc_x[k] * qsp_120[k];

        t_241[k] = f_12 * osp_121[k]
                   + f_3 * pc_x[k] * qsp_121[k];
    }
}

static auto
compute_prim_qsd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osd0,
                                                          const size_t osp, const size_t osd1,
                                                          const size_t qss0, const size_t qss1,
                                                          const size_t qsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 5.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osd0_210 = buffer.data(osd0 + 210);
    const auto *osd0_215 = buffer.data(osd0 + 215);
    const auto *osd0_216 = buffer.data(osd0 + 216);
    const auto *osd0_219 = buffer.data(osd0 + 219);
    const auto *osd0_264 = buffer.data(osd0 + 264);
    const auto *osd0_269 = buffer.data(osd0 + 269);
    const auto *osd0_270 = buffer.data(osd0 + 270);
    const auto *osd0_273 = buffer.data(osd0 + 273);

    const auto *osp_95 = buffer.data(osp + 95);
    const auto *osp_97 = buffer.data(osp + 97);
    const auto *osp_98 = buffer.data(osp + 98);
    const auto *osp_100 = buffer.data(osp + 100);
    const auto *osp_101 = buffer.data(osp + 101);
    const auto *osp_103 = buffer.data(osp + 103);
    const auto *osp_104 = buffer.data(osp + 104);
    const auto *osp_106 = buffer.data(osp + 106);
    const auto *osp_107 = buffer.data(osp + 107);
    const auto *osp_109 = buffer.data(osp + 109);
    const auto *osp_110 = buffer.data(osp + 110);
    const auto *osp_113 = buffer.data(osp + 113);
    const auto *osp_115 = buffer.data(osp + 115);
    const auto *osp_116 = buffer.data(osp + 116);
    const auto *osp_118 = buffer.data(osp + 118);
    const auto *osp_119 = buffer.data(osp + 119);
    const auto *osp_121 = buffer.data(osp + 121);
    const auto *osp_122 = buffer.data(osp + 122);
    const auto *osp_123 = buffer.data(osp + 123);
    const auto *osp_124 = buffer.data(osp + 124);
    const auto *osp_125 = buffer.data(osp + 125);
    const auto *osp_126 = buffer.data(osp + 126);
    const auto *osp_127 = buffer.data(osp + 127);
    const auto *osp_128 = buffer.data(osp + 128);
    const auto *osp_130 = buffer.data(osp + 130);
    const auto *osp_131 = buffer.data(osp + 131);
    const auto *osp_132 = buffer.data(osp + 132);
    const auto *osp_133 = buffer.data(osp + 133);
    const auto *osp_134 = buffer.data(osp + 134);
    const auto *osp_135 = buffer.data(osp + 135);
    const auto *osp_136 = buffer.data(osp + 136);
    const auto *osp_137 = buffer.data(osp + 137);
    const auto *osp_139 = buffer.data(osp + 139);
    const auto *osp_140 = buffer.data(osp + 140);
    const auto *osp_141 = buffer.data(osp + 141);
    const auto *osp_142 = buffer.data(osp + 142);
    const auto *osp_143 = buffer.data(osp + 143);
    const auto *osp_144 = buffer.data(osp + 144);
    const auto *osp_145 = buffer.data(osp + 145);
    const auto *osp_146 = buffer.data(osp + 146);
    const auto *osp_147 = buffer.data(osp + 147);
    const auto *osp_148 = buffer.data(osp + 148);
    const auto *osp_149 = buffer.data(osp + 149);
    const auto *osp_150 = buffer.data(osp + 150);
    const auto *osp_151 = buffer.data(osp + 151);
    const auto *osp_152 = buffer.data(osp + 152);
    const auto *osp_153 = buffer.data(osp + 153);
    const auto *osp_154 = buffer.data(osp + 154);
    const auto *osp_155 = buffer.data(osp + 155);
    const auto *osp_156 = buffer.data(osp + 156);
    const auto *osp_157 = buffer.data(osp + 157);
    const auto *osp_158 = buffer.data(osp + 158);
    const auto *osp_160 = buffer.data(osp + 160);
    const auto *osp_161 = buffer.data(osp + 161);
    const auto *osp_162 = buffer.data(osp + 162);
    const auto *osp_164 = buffer.data(osp + 164);
    const auto *osp_165 = buffer.data(osp + 165);
    const auto *osp_166 = buffer.data(osp + 166);
    const auto *osp_169 = buffer.data(osp + 169);
    const auto *osp_170 = buffer.data(osp + 170);
    const auto *osp_171 = buffer.data(osp + 171);
    const auto *osp_172 = buffer.data(osp + 172);
    const auto *osp_173 = buffer.data(osp + 173);
    const auto *osp_174 = buffer.data(osp + 174);
    const auto *osp_175 = buffer.data(osp + 175);
    const auto *osp_176 = buffer.data(osp + 176);
    const auto *osp_177 = buffer.data(osp + 177);
    const auto *osp_178 = buffer.data(osp + 178);
    const auto *osp_179 = buffer.data(osp + 179);

    const auto *osd1_210 = buffer.data(osd1 + 210);
    const auto *osd1_215 = buffer.data(osd1 + 215);
    const auto *osd1_216 = buffer.data(osd1 + 216);
    const auto *osd1_219 = buffer.data(osd1 + 219);
    const auto *osd1_264 = buffer.data(osd1 + 264);
    const auto *osd1_269 = buffer.data(osd1 + 269);
    const auto *osd1_270 = buffer.data(osd1 + 270);
    const auto *osd1_273 = buffer.data(osd1 + 273);

    const auto *qss0_40 = buffer.data(qss0 + 40);
    const auto *qss0_41 = buffer.data(qss0 + 41);
    const auto *qss0_42 = buffer.data(qss0 + 42);
    const auto *qss0_43 = buffer.data(qss0 + 43);
    const auto *qss0_44 = buffer.data(qss0 + 44);
    const auto *qss0_45 = buffer.data(qss0 + 45);
    const auto *qss0_46 = buffer.data(qss0 + 46);
    const auto *qss0_47 = buffer.data(qss0 + 47);
    const auto *qss0_48 = buffer.data(qss0 + 48);
    const auto *qss0_49 = buffer.data(qss0 + 49);
    const auto *qss0_50 = buffer.data(qss0 + 50);
    const auto *qss0_51 = buffer.data(qss0 + 51);
    const auto *qss0_52 = buffer.data(qss0 + 52);
    const auto *qss0_53 = buffer.data(qss0 + 53);
    const auto *qss0_54 = buffer.data(qss0 + 54);
    const auto *qss0_55 = buffer.data(qss0 + 55);
    const auto *qss0_56 = buffer.data(qss0 + 56);
    const auto *qss0_57 = buffer.data(qss0 + 57);
    const auto *qss0_58 = buffer.data(qss0 + 58);
    const auto *qss0_59 = buffer.data(qss0 + 59);

    const auto *qss1_40 = buffer.data(qss1 + 40);
    const auto *qss1_41 = buffer.data(qss1 + 41);
    const auto *qss1_42 = buffer.data(qss1 + 42);
    const auto *qss1_43 = buffer.data(qss1 + 43);
    const auto *qss1_44 = buffer.data(qss1 + 44);
    const auto *qss1_45 = buffer.data(qss1 + 45);
    const auto *qss1_46 = buffer.data(qss1 + 46);
    const auto *qss1_47 = buffer.data(qss1 + 47);
    const auto *qss1_48 = buffer.data(qss1 + 48);
    const auto *qss1_49 = buffer.data(qss1 + 49);
    const auto *qss1_50 = buffer.data(qss1 + 50);
    const auto *qss1_51 = buffer.data(qss1 + 51);
    const auto *qss1_52 = buffer.data(qss1 + 52);
    const auto *qss1_53 = buffer.data(qss1 + 53);
    const auto *qss1_54 = buffer.data(qss1 + 54);
    const auto *qss1_55 = buffer.data(qss1 + 55);
    const auto *qss1_56 = buffer.data(qss1 + 56);
    const auto *qss1_57 = buffer.data(qss1 + 57);
    const auto *qss1_58 = buffer.data(qss1 + 58);
    const auto *qss1_59 = buffer.data(qss1 + 59);

    const auto *qsp_121 = buffer.data(qsp + 121);
    const auto *qsp_122 = buffer.data(qsp + 122);
    const auto *qsp_123 = buffer.data(qsp + 123);
    const auto *qsp_124 = buffer.data(qsp + 124);
    const auto *qsp_125 = buffer.data(qsp + 125);
    const auto *qsp_126 = buffer.data(qsp + 126);
    const auto *qsp_127 = buffer.data(qsp + 127);
    const auto *qsp_128 = buffer.data(qsp + 128);
    const auto *qsp_130 = buffer.data(qsp + 130);
    const auto *qsp_131 = buffer.data(qsp + 131);
    const auto *qsp_132 = buffer.data(qsp + 132);
    const auto *qsp_133 = buffer.data(qsp + 133);
    const auto *qsp_134 = buffer.data(qsp + 134);
    const auto *qsp_135 = buffer.data(qsp + 135);
    const auto *qsp_136 = buffer.data(qsp + 136);
    const auto *qsp_137 = buffer.data(qsp + 137);
    const auto *qsp_139 = buffer.data(qsp + 139);
    const auto *qsp_140 = buffer.data(qsp + 140);
    const auto *qsp_141 = buffer.data(qsp + 141);
    const auto *qsp_142 = buffer.data(qsp + 142);
    const auto *qsp_143 = buffer.data(qsp + 143);
    const auto *qsp_144 = buffer.data(qsp + 144);
    const auto *qsp_145 = buffer.data(qsp + 145);
    const auto *qsp_146 = buffer.data(qsp + 146);
    const auto *qsp_147 = buffer.data(qsp + 147);
    const auto *qsp_148 = buffer.data(qsp + 148);
    const auto *qsp_149 = buffer.data(qsp + 149);
    const auto *qsp_150 = buffer.data(qsp + 150);
    const auto *qsp_151 = buffer.data(qsp + 151);
    const auto *qsp_152 = buffer.data(qsp + 152);
    const auto *qsp_153 = buffer.data(qsp + 153);
    const auto *qsp_154 = buffer.data(qsp + 154);
    const auto *qsp_155 = buffer.data(qsp + 155);
    const auto *qsp_156 = buffer.data(qsp + 156);
    const auto *qsp_157 = buffer.data(qsp + 157);
    const auto *qsp_158 = buffer.data(qsp + 158);
    const auto *qsp_160 = buffer.data(qsp + 160);
    const auto *qsp_161 = buffer.data(qsp + 161);
    const auto *qsp_162 = buffer.data(qsp + 162);
    const auto *qsp_163 = buffer.data(qsp + 163);
    const auto *qsp_164 = buffer.data(qsp + 164);
    const auto *qsp_165 = buffer.data(qsp + 165);
    const auto *qsp_166 = buffer.data(qsp + 166);
    const auto *qsp_167 = buffer.data(qsp + 167);
    const auto *qsp_169 = buffer.data(qsp + 169);
    const auto *qsp_170 = buffer.data(qsp + 170);
    const auto *qsp_171 = buffer.data(qsp + 171);
    const auto *qsp_172 = buffer.data(qsp + 172);
    const auto *qsp_173 = buffer.data(qsp + 173);
    const auto *qsp_174 = buffer.data(qsp + 174);
    const auto *qsp_175 = buffer.data(qsp + 175);
    const auto *qsp_176 = buffer.data(qsp + 176);
    const auto *qsp_177 = buffer.data(qsp + 177);
    const auto *qsp_178 = buffer.data(qsp + 178);
    const auto *qsp_179 = buffer.data(qsp + 179);

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pc_x, pc_y, pc_z, osp_95, osp_97, osp_98, \
                         osp_122, qss0_40, qss1_40, qsp_121, qsp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_12 * osp_122[k]
                   + f_3 * pc_x[k] * qsp_122[k];

        t_243[k] = f_12 * osp_97[k]
                   + f_1 * qss0_40[k]
                   - f_2 * qss1_40[k]
                   + f_3 * pc_y[k] * qsp_121[k];

        t_244[k] = f_12 * osp_98[k]
                   + f_3 * pc_y[k] * qsp_122[k];

        t_245[k] = f_12 * osp_95[k]
                   + f_1 * qss0_40[k]
                   - f_2 * qss1_40[k]
                   + f_3 * pc_z[k] * qsp_122[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pc_y, osp_100, osp_123, osp_124, \
                         osp_125, qss0_41, qss1_41, qsp_123, qsp_124, \
                         qsp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * osp_123[k]
                   + f_1 * qss0_41[k]
                   - f_2 * qss1_41[k]
                   + f_3 * pc_x[k] * qsp_123[k];

        t_247[k] = f_12 * osp_124[k]
                   + f_3 * pc_x[k] * qsp_124[k];

        t_248[k] = f_12 * osp_125[k]
                   + f_3 * pc_x[k] * qsp_125[k];

        t_249[k] = f_10 * osp_100[k]
                   + f_1 * qss0_41[k]
                   - f_2 * qss1_41[k]
                   + f_3 * pc_y[k] * qsp_124[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, osp_98, osp_101, osp_126, \
                         qss0_41, qss0_42, qss1_41, qss1_42, qsp_125, \
                         qsp_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_10 * osp_101[k]
                   + f_3 * pc_y[k] * qsp_125[k];

        t_251[k] = f_14 * osp_98[k]
                   + f_1 * qss0_41[k]
                   - f_2 * qss1_41[k]
                   + f_3 * pc_z[k] * qsp_125[k];

        t_252[k] = f_12 * osp_126[k]
                   + f_1 * qss0_42[k]
                   - f_2 * qss1_42[k]
                   + f_3 * pc_x[k] * qsp_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, osp_103, osp_104, osp_127, \
                         osp_128, qss0_42, qss1_42, qsp_127, qsp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_12 * osp_127[k]
                   + f_3 * pc_x[k] * qsp_127[k];

        t_254[k] = f_12 * osp_128[k]
                   + f_3 * pc_x[k] * qsp_128[k];

        t_255[k] = f_8 * osp_103[k]
                   + f_1 * qss0_42[k]
                   - f_2 * qss1_42[k]
                   + f_3 * pc_y[k] * qsp_127[k];

        t_256[k] = f_8 * osp_104[k]
                   + f_3 * pc_y[k] * qsp_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_y, pc_x, pc_y, pc_z, osd0_210, osp_101, \
                         osp_130, osd1_210, qss0_42, qss1_42, qsp_128, \
                         qsp_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * osp_101[k]
                   + f_1 * qss0_42[k]
                   - f_2 * qss1_42[k]
                   + f_3 * pc_z[k] * qsp_128[k];

        t_258[k] = pa_y[k] * osd0_210[k]
                   - f_4 * pc_y[k] * osd1_210[k];

        t_259[k] = f_12 * osp_130[k]
                   + f_3 * pc_x[k] * qsp_130[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, osd0_215, osp_106, \
                         osp_107, osp_131, osd1_215, qss0_43, qss1_43, qsp_130, \
                         qsp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_12 * osp_131[k]
                   + f_3 * pc_x[k] * qsp_131[k];

        t_261[k] = f_6 * osp_106[k]
                   + f_1 * qss0_43[k]
                   - f_2 * qss1_43[k]
                   + f_3 * pc_y[k] * qsp_130[k];

        t_262[k] = f_6 * osp_107[k]
                   + f_3 * pc_y[k] * qsp_131[k];

        t_263[k] = pa_y[k] * osd0_215[k]
                   - f_4 * pc_y[k] * osd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, osp_132, osp_134, \
                         qss0_44, qss1_44, qsp_132, qsp_133, qsp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_12 * osp_132[k]
                   + f_1 * qss0_44[k]
                   - f_2 * qss1_44[k]
                   + f_3 * pc_x[k] * qsp_132[k];

        t_265[k] = f_3 * pc_y[k] * qsp_132[k];

        t_266[k] = f_12 * osp_134[k]
                   + f_3 * pc_x[k] * qsp_134[k];

        t_267[k] = f_1 * qss0_44[k]
                   - f_2 * qss1_44[k]
                   + f_3 * pc_y[k] * qsp_133[k];

        t_268[k] = f_3 * pc_y[k] * qsp_134[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pc_x, pc_z, osp_107, osp_135, osp_136, \
                         qss0_44, qss0_45, qss1_44, qss1_45, qsp_134, qsp_135, \
                         qsp_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_11 * osp_107[k]
                   + f_1 * qss0_44[k]
                   - f_2 * qss1_44[k]
                   + f_3 * pc_z[k] * qsp_134[k];

        t_270[k] = f_10 * osp_135[k]
                   + f_1 * qss0_45[k]
                   - f_2 * qss1_45[k]
                   + f_3 * pc_x[k] * qsp_135[k];

        t_271[k] = f_10 * osp_136[k]
                   + f_3 * pc_x[k] * qsp_136[k];

        t_272[k] = f_3 * pc_z[k] * qsp_135[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pa_z, pc_y, pc_z, osd0_216, osp_109, \
                         osd1_216, qss0_45, qss1_45, qsp_136, qsp_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_9 * osp_109[k]
                   + f_1 * qss0_45[k]
                   - f_2 * qss1_45[k]
                   + f_3 * pc_y[k] * qsp_136[k];

        t_274[k] = f_3 * pc_z[k] * qsp_136[k];

        t_275[k] = f_1 * qss0_45[k]
                   - f_2 * qss1_45[k]
                   + f_3 * pc_z[k] * qsp_137[k];

        t_276[k] = pa_z[k] * osd0_216[k]
                   - f_4 * pc_z[k] * osd1_216[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pa_z, pc_x, pc_y, pc_z, osd0_219, \
                         osp_113, osp_139, osp_140, osd1_219, qsp_139, \
                         qsp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_10 * osp_139[k]
                   + f_3 * pc_x[k] * qsp_139[k];

        t_278[k] = f_10 * osp_140[k]
                   + f_3 * pc_x[k] * qsp_140[k];

        t_279[k] = pa_z[k] * osd0_219[k]
                   - f_4 * pc_z[k] * osd1_219[k];

        t_280[k] = f_11 * osp_113[k]
                   + f_3 * pc_y[k] * qsp_140[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_z, osp_110, osp_141, osp_142, qss0_46, \
                         qss0_47, qss1_46, qss1_47, qsp_140, qsp_141, \
                         qsp_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_6 * osp_110[k]
                   + f_1 * qss0_46[k]
                   - f_2 * qss1_46[k]
                   + f_3 * pc_z[k] * qsp_140[k];

        t_282[k] = f_10 * osp_141[k]
                   + f_1 * qss0_47[k]
                   - f_2 * qss1_47[k]
                   + f_3 * pc_x[k] * qsp_141[k];

        t_283[k] = f_10 * osp_142[k]
                   + f_3 * pc_x[k] * qsp_142[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_y, pc_z, osp_113, osp_115, \
                         osp_116, osp_143, qss0_47, qss1_47, qsp_142, \
                         qsp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_10 * osp_143[k]
                   + f_3 * pc_x[k] * qsp_143[k];

        t_285[k] = f_13 * osp_115[k]
                   + f_1 * qss0_47[k]
                   - f_2 * qss1_47[k]
                   + f_3 * pc_y[k] * qsp_142[k];

        t_286[k] = f_13 * osp_116[k]
                   + f_3 * pc_y[k] * qsp_143[k];

        t_287[k] = f_8 * osp_113[k]
                   + f_1 * qss0_47[k]
                   - f_2 * qss1_47[k]
                   + f_3 * pc_z[k] * qsp_143[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, osp_118, osp_144, osp_145, \
                         osp_146, qss0_48, qss1_48, qsp_144, qsp_145, \
                         qsp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_10 * osp_144[k]
                   + f_1 * qss0_48[k]
                   - f_2 * qss1_48[k]
                   + f_3 * pc_x[k] * qsp_144[k];

        t_289[k] = f_10 * osp_145[k]
                   + f_3 * pc_x[k] * qsp_145[k];

        t_290[k] = f_10 * osp_146[k]
                   + f_3 * pc_x[k] * qsp_146[k];

        t_291[k] = f_15 * osp_118[k]
                   + f_1 * qss0_48[k]
                   - f_2 * qss1_48[k]
                   + f_3 * pc_y[k] * qsp_145[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, pc_x, pc_y, pc_z, osp_116, osp_119, osp_147, \
                         qss0_48, qss0_49, qss1_48, qss1_49, qsp_146, \
                         qsp_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_15 * osp_119[k]
                   + f_3 * pc_y[k] * qsp_146[k];

        t_293[k] = f_10 * osp_116[k]
                   + f_1 * qss0_48[k]
                   - f_2 * qss1_48[k]
                   + f_3 * pc_z[k] * qsp_146[k];

        t_294[k] = f_10 * osp_147[k]
                   + f_1 * qss0_49[k]
                   - f_2 * qss1_49[k]
                   + f_3 * pc_x[k] * qsp_147[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_x, pc_y, osp_121, osp_122, osp_148, \
                         osp_149, qss0_49, qss1_49, qsp_148, qsp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_10 * osp_148[k]
                   + f_3 * pc_x[k] * qsp_148[k];

        t_296[k] = f_10 * osp_149[k]
                   + f_3 * pc_x[k] * qsp_149[k];

        t_297[k] = f_14 * osp_121[k]
                   + f_1 * qss0_49[k]
                   - f_2 * qss1_49[k]
                   + f_3 * pc_y[k] * qsp_148[k];

        t_298[k] = f_14 * osp_122[k]
                   + f_3 * pc_y[k] * qsp_149[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pc_x, pc_z, osp_119, osp_150, osp_151, qss0_49, \
                         qss0_50, qss1_49, qss1_50, qsp_149, qsp_150, \
                         qsp_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_12 * osp_119[k]
                   + f_1 * qss0_49[k]
                   - f_2 * qss1_49[k]
                   + f_3 * pc_z[k] * qsp_149[k];

        t_300[k] = f_10 * osp_150[k]
                   + f_1 * qss0_50[k]
                   - f_2 * qss1_50[k]
                   + f_3 * pc_x[k] * qsp_150[k];

        t_301[k] = f_10 * osp_151[k]
                   + f_3 * pc_x[k] * qsp_151[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_y, pc_z, osp_122, osp_124, \
                         osp_125, osp_152, qss0_50, qss1_50, qsp_151, \
                         qsp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_10 * osp_152[k]
                   + f_3 * pc_x[k] * qsp_152[k];

        t_303[k] = f_12 * osp_124[k]
                   + f_1 * qss0_50[k]
                   - f_2 * qss1_50[k]
                   + f_3 * pc_y[k] * qsp_151[k];

        t_304[k] = f_12 * osp_125[k]
                   + f_3 * pc_y[k] * qsp_152[k];

        t_305[k] = f_14 * osp_122[k]
                   + f_1 * qss0_50[k]
                   - f_2 * qss1_50[k]
                   + f_3 * pc_z[k] * qsp_152[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_x, pc_y, osp_127, osp_153, osp_154, \
                         osp_155, qss0_51, qss1_51, qsp_153, qsp_154, \
                         qsp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_10 * osp_153[k]
                   + f_1 * qss0_51[k]
                   - f_2 * qss1_51[k]
                   + f_3 * pc_x[k] * qsp_153[k];

        t_307[k] = f_10 * osp_154[k]
                   + f_3 * pc_x[k] * qsp_154[k];

        t_308[k] = f_10 * osp_155[k]
                   + f_3 * pc_x[k] * qsp_155[k];

        t_309[k] = f_10 * osp_127[k]
                   + f_1 * qss0_51[k]
                   - f_2 * qss1_51[k]
                   + f_3 * pc_y[k] * qsp_154[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pc_x, pc_y, pc_z, osp_125, osp_128, osp_156, \
                         qss0_51, qss0_52, qss1_51, qss1_52, qsp_155, \
                         qsp_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_10 * osp_128[k]
                   + f_3 * pc_y[k] * qsp_155[k];

        t_311[k] = f_15 * osp_125[k]
                   + f_1 * qss0_51[k]
                   - f_2 * qss1_51[k]
                   + f_3 * pc_z[k] * qsp_155[k];

        t_312[k] = f_10 * osp_156[k]
                   + f_1 * qss0_52[k]
                   - f_2 * qss1_52[k]
                   + f_3 * pc_x[k] * qsp_156[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, osp_130, osp_131, osp_157, \
                         osp_158, qss0_52, qss1_52, qsp_157, qsp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_10 * osp_157[k]
                   + f_3 * pc_x[k] * qsp_157[k];

        t_314[k] = f_10 * osp_158[k]
                   + f_3 * pc_x[k] * qsp_158[k];

        t_315[k] = f_8 * osp_130[k]
                   + f_1 * qss0_52[k]
                   - f_2 * qss1_52[k]
                   + f_3 * pc_y[k] * qsp_157[k];

        t_316[k] = f_8 * osp_131[k]
                   + f_3 * pc_y[k] * qsp_158[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pa_y, pc_x, pc_y, pc_z, osd0_264, osp_128, \
                         osp_160, osd1_264, qss0_52, qss1_52, qsp_158, \
                         qsp_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_13 * osp_128[k]
                   + f_1 * qss0_52[k]
                   - f_2 * qss1_52[k]
                   + f_3 * pc_z[k] * qsp_158[k];

        t_318[k] = pa_y[k] * osd0_264[k]
                   - f_4 * pc_y[k] * osd1_264[k];

        t_319[k] = f_10 * osp_160[k]
                   + f_3 * pc_x[k] * qsp_160[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_y, pc_x, pc_y, osd0_269, osp_133, \
                         osp_134, osp_161, osd1_269, qss0_53, qss1_53, qsp_160, \
                         qsp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_10 * osp_161[k]
                   + f_3 * pc_x[k] * qsp_161[k];

        t_321[k] = f_6 * osp_133[k]
                   + f_1 * qss0_53[k]
                   - f_2 * qss1_53[k]
                   + f_3 * pc_y[k] * qsp_160[k];

        t_322[k] = f_6 * osp_134[k]
                   + f_3 * pc_y[k] * qsp_161[k];

        t_323[k] = pa_y[k] * osd0_269[k]
                   - f_4 * pc_y[k] * osd1_269[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, osp_162, osp_164, \
                         qss0_54, qss1_54, qsp_162, qsp_163, qsp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_10 * osp_162[k]
                   + f_1 * qss0_54[k]
                   - f_2 * qss1_54[k]
                   + f_3 * pc_x[k] * qsp_162[k];

        t_325[k] = f_3 * pc_y[k] * qsp_162[k];

        t_326[k] = f_10 * osp_164[k]
                   + f_3 * pc_x[k] * qsp_164[k];

        t_327[k] = f_1 * qss0_54[k]
                   - f_2 * qss1_54[k]
                   + f_3 * pc_y[k] * qsp_163[k];

        t_328[k] = f_3 * pc_y[k] * qsp_164[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_z, osp_134, osp_165, osp_166, \
                         qss0_54, qss0_55, qss1_54, qss1_55, qsp_164, qsp_165, \
                         qsp_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_9 * osp_134[k]
                   + f_1 * qss0_54[k]
                   - f_2 * qss1_54[k]
                   + f_3 * pc_z[k] * qsp_164[k];

        t_330[k] = f_8 * osp_165[k]
                   + f_1 * qss0_55[k]
                   - f_2 * qss1_55[k]
                   + f_3 * pc_x[k] * qsp_165[k];

        t_331[k] = f_8 * osp_166[k]
                   + f_3 * pc_x[k] * qsp_166[k];

        t_332[k] = f_3 * pc_z[k] * qsp_165[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pa_z, pc_y, pc_z, osd0_270, osp_136, \
                         osd1_270, qss0_55, qss1_55, qsp_166, qsp_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_7 * osp_136[k]
                   + f_1 * qss0_55[k]
                   - f_2 * qss1_55[k]
                   + f_3 * pc_y[k] * qsp_166[k];

        t_334[k] = f_3 * pc_z[k] * qsp_166[k];

        t_335[k] = f_1 * qss0_55[k]
                   - f_2 * qss1_55[k]
                   + f_3 * pc_z[k] * qsp_167[k];

        t_336[k] = pa_z[k] * osd0_270[k]
                   - f_4 * pc_z[k] * osd1_270[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pc_x, pc_y, pc_z, osd0_273, \
                         osp_140, osp_169, osp_170, osd1_273, qsp_169, \
                         qsp_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_8 * osp_169[k]
                   + f_3 * pc_x[k] * qsp_169[k];

        t_338[k] = f_8 * osp_170[k]
                   + f_3 * pc_x[k] * qsp_170[k];

        t_339[k] = pa_z[k] * osd0_273[k]
                   - f_4 * pc_z[k] * osd1_273[k];

        t_340[k] = f_9 * osp_140[k]
                   + f_3 * pc_y[k] * qsp_170[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pc_x, pc_z, osp_137, osp_171, osp_172, qss0_56, \
                         qss0_57, qss1_56, qss1_57, qsp_170, qsp_171, \
                         qsp_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_6 * osp_137[k]
                   + f_1 * qss0_56[k]
                   - f_2 * qss1_56[k]
                   + f_3 * pc_z[k] * qsp_170[k];

        t_342[k] = f_8 * osp_171[k]
                   + f_1 * qss0_57[k]
                   - f_2 * qss1_57[k]
                   + f_3 * pc_x[k] * qsp_171[k];

        t_343[k] = f_8 * osp_172[k]
                   + f_3 * pc_x[k] * qsp_172[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, pc_z, osp_140, osp_142, \
                         osp_143, osp_173, qss0_57, qss1_57, qsp_172, \
                         qsp_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_8 * osp_173[k]
                   + f_3 * pc_x[k] * qsp_173[k];

        t_345[k] = f_11 * osp_142[k]
                   + f_1 * qss0_57[k]
                   - f_2 * qss1_57[k]
                   + f_3 * pc_y[k] * qsp_172[k];

        t_346[k] = f_11 * osp_143[k]
                   + f_3 * pc_y[k] * qsp_173[k];

        t_347[k] = f_8 * osp_140[k]
                   + f_1 * qss0_57[k]
                   - f_2 * qss1_57[k]
                   + f_3 * pc_z[k] * qsp_173[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pc_x, pc_y, osp_145, osp_174, osp_175, \
                         osp_176, qss0_58, qss1_58, qsp_174, qsp_175, \
                         qsp_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_8 * osp_174[k]
                   + f_1 * qss0_58[k]
                   - f_2 * qss1_58[k]
                   + f_3 * pc_x[k] * qsp_174[k];

        t_349[k] = f_8 * osp_175[k]
                   + f_3 * pc_x[k] * qsp_175[k];

        t_350[k] = f_8 * osp_176[k]
                   + f_3 * pc_x[k] * qsp_176[k];

        t_351[k] = f_13 * osp_145[k]
                   + f_1 * qss0_58[k]
                   - f_2 * qss1_58[k]
                   + f_3 * pc_y[k] * qsp_175[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_x, pc_y, pc_z, osp_143, osp_146, osp_177, \
                         qss0_58, qss0_59, qss1_58, qss1_59, qsp_176, \
                         qsp_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_13 * osp_146[k]
                   + f_3 * pc_y[k] * qsp_176[k];

        t_353[k] = f_10 * osp_143[k]
                   + f_1 * qss0_58[k]
                   - f_2 * qss1_58[k]
                   + f_3 * pc_z[k] * qsp_176[k];

        t_354[k] = f_8 * osp_177[k]
                   + f_1 * qss0_59[k]
                   - f_2 * qss1_59[k]
                   + f_3 * pc_x[k] * qsp_177[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_x, pc_y, osp_148, osp_149, osp_178, \
                         osp_179, qss0_59, qss1_59, qsp_178, qsp_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_8 * osp_178[k]
                   + f_3 * pc_x[k] * qsp_178[k];

        t_356[k] = f_8 * osp_179[k]
                   + f_3 * pc_x[k] * qsp_179[k];

        t_357[k] = f_15 * osp_148[k]
                   + f_1 * qss0_59[k]
                   - f_2 * qss1_59[k]
                   + f_3 * pc_y[k] * qsp_178[k];

        t_358[k] = f_15 * osp_149[k]
                   + f_3 * pc_y[k] * qsp_179[k];
    }
}

static auto
compute_prim_qsd_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osd0,
                                                          const size_t osp, const size_t osd1,
                                                          const size_t qss0, const size_t qss1,
                                                          const size_t qsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 5.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osd0_324 = buffer.data(osd0 + 324);
    const auto *osd0_329 = buffer.data(osd0 + 329);
    const auto *osd0_330 = buffer.data(osd0 + 330);
    const auto *osd0_390 = buffer.data(osd0 + 390);
    const auto *osd0_396 = buffer.data(osd0 + 396);
    const auto *osd0_399 = buffer.data(osd0 + 399);
    const auto *osd0_401 = buffer.data(osd0 + 401);
    const auto *osd0_405 = buffer.data(osd0 + 405);
    const auto *osd0_407 = buffer.data(osd0 + 407);
    const auto *osd0_408 = buffer.data(osd0 + 408);
    const auto *osd0_411 = buffer.data(osd0 + 411);
    const auto *osd0_413 = buffer.data(osd0 + 413);
    const auto *osd0_414 = buffer.data(osd0 + 414);
    const auto *osd0_417 = buffer.data(osd0 + 417);
    const auto *osd0_419 = buffer.data(osd0 + 419);
    const auto *osd0_420 = buffer.data(osd0 + 420);
    const auto *osd0_423 = buffer.data(osd0 + 423);
    const auto *osd0_425 = buffer.data(osd0 + 425);
    const auto *osd0_426 = buffer.data(osd0 + 426);
    const auto *osd0_429 = buffer.data(osd0 + 429);
    const auto *osd0_431 = buffer.data(osd0 + 431);
    const auto *osd0_432 = buffer.data(osd0 + 432);
    const auto *osd0_435 = buffer.data(osd0 + 435);
    const auto *osd0_437 = buffer.data(osd0 + 437);
    const auto *osd0_438 = buffer.data(osd0 + 438);
    const auto *osd0_441 = buffer.data(osd0 + 441);
    const auto *osd0_443 = buffer.data(osd0 + 443);
    const auto *osd0_444 = buffer.data(osd0 + 444);
    const auto *osd0_447 = buffer.data(osd0 + 447);
    const auto *osd0_449 = buffer.data(osd0 + 449);
    const auto *osd0_450 = buffer.data(osd0 + 450);
    const auto *osd0_453 = buffer.data(osd0 + 453);
    const auto *osd0_455 = buffer.data(osd0 + 455);
    const auto *osd0_459 = buffer.data(osd0 + 459);
    const auto *osd0_461 = buffer.data(osd0 + 461);
    const auto *osd0_462 = buffer.data(osd0 + 462);
    const auto *osd0_465 = buffer.data(osd0 + 465);
    const auto *osd0_467 = buffer.data(osd0 + 467);

    const auto *osp_146 = buffer.data(osp + 146);
    const auto *osp_149 = buffer.data(osp + 149);
    const auto *osp_151 = buffer.data(osp + 151);
    const auto *osp_152 = buffer.data(osp + 152);
    const auto *osp_154 = buffer.data(osp + 154);
    const auto *osp_155 = buffer.data(osp + 155);
    const auto *osp_157 = buffer.data(osp + 157);
    const auto *osp_158 = buffer.data(osp + 158);
    const auto *osp_160 = buffer.data(osp + 160);
    const auto *osp_161 = buffer.data(osp + 161);
    const auto *osp_163 = buffer.data(osp + 163);
    const auto *osp_164 = buffer.data(osp + 164);
    const auto *osp_170 = buffer.data(osp + 170);
    const auto *osp_173 = buffer.data(osp + 173);
    const auto *osp_176 = buffer.data(osp + 176);
    const auto *osp_179 = buffer.data(osp + 179);
    const auto *osp_180 = buffer.data(osp + 180);
    const auto *osp_181 = buffer.data(osp + 181);
    const auto *osp_182 = buffer.data(osp + 182);
    const auto *osp_183 = buffer.data(osp + 183);
    const auto *osp_184 = buffer.data(osp + 184);
    const auto *osp_185 = buffer.data(osp + 185);
    const auto *osp_186 = buffer.data(osp + 186);
    const auto *osp_187 = buffer.data(osp + 187);
    const auto *osp_188 = buffer.data(osp + 188);
    const auto *osp_189 = buffer.data(osp + 189);
    const auto *osp_190 = buffer.data(osp + 190);
    const auto *osp_191 = buffer.data(osp + 191);
    const auto *osp_193 = buffer.data(osp + 193);
    const auto *osp_194 = buffer.data(osp + 194);
    const auto *osp_195 = buffer.data(osp + 195);
    const auto *osp_197 = buffer.data(osp + 197);
    const auto *osp_198 = buffer.data(osp + 198);
    const auto *osp_199 = buffer.data(osp + 199);
    const auto *osp_200 = buffer.data(osp + 200);
    const auto *osp_202 = buffer.data(osp + 202);
    const auto *osp_203 = buffer.data(osp + 203);
    const auto *osp_204 = buffer.data(osp + 204);
    const auto *osp_205 = buffer.data(osp + 205);
    const auto *osp_206 = buffer.data(osp + 206);
    const auto *osp_207 = buffer.data(osp + 207);
    const auto *osp_208 = buffer.data(osp + 208);
    const auto *osp_209 = buffer.data(osp + 209);
    const auto *osp_210 = buffer.data(osp + 210);
    const auto *osp_211 = buffer.data(osp + 211);
    const auto *osp_212 = buffer.data(osp + 212);
    const auto *osp_213 = buffer.data(osp + 213);
    const auto *osp_214 = buffer.data(osp + 214);
    const auto *osp_215 = buffer.data(osp + 215);
    const auto *osp_216 = buffer.data(osp + 216);
    const auto *osp_217 = buffer.data(osp + 217);
    const auto *osp_218 = buffer.data(osp + 218);
    const auto *osp_219 = buffer.data(osp + 219);
    const auto *osp_220 = buffer.data(osp + 220);
    const auto *osp_221 = buffer.data(osp + 221);
    const auto *osp_222 = buffer.data(osp + 222);
    const auto *osp_223 = buffer.data(osp + 223);
    const auto *osp_224 = buffer.data(osp + 224);
    const auto *osp_225 = buffer.data(osp + 225);
    const auto *osp_226 = buffer.data(osp + 226);
    const auto *osp_227 = buffer.data(osp + 227);
    const auto *osp_229 = buffer.data(osp + 229);
    const auto *osp_230 = buffer.data(osp + 230);
    const auto *osp_231 = buffer.data(osp + 231);
    const auto *osp_233 = buffer.data(osp + 233);

    const auto *osd1_324 = buffer.data(osd1 + 324);
    const auto *osd1_329 = buffer.data(osd1 + 329);
    const auto *osd1_330 = buffer.data(osd1 + 330);
    const auto *osd1_390 = buffer.data(osd1 + 390);
    const auto *osd1_396 = buffer.data(osd1 + 396);
    const auto *osd1_399 = buffer.data(osd1 + 399);
    const auto *osd1_401 = buffer.data(osd1 + 401);
    const auto *osd1_405 = buffer.data(osd1 + 405);
    const auto *osd1_407 = buffer.data(osd1 + 407);
    const auto *osd1_408 = buffer.data(osd1 + 408);
    const auto *osd1_411 = buffer.data(osd1 + 411);
    const auto *osd1_413 = buffer.data(osd1 + 413);
    const auto *osd1_414 = buffer.data(osd1 + 414);
    const auto *osd1_417 = buffer.data(osd1 + 417);
    const auto *osd1_419 = buffer.data(osd1 + 419);
    const auto *osd1_420 = buffer.data(osd1 + 420);
    const auto *osd1_423 = buffer.data(osd1 + 423);
    const auto *osd1_425 = buffer.data(osd1 + 425);
    const auto *osd1_426 = buffer.data(osd1 + 426);
    const auto *osd1_429 = buffer.data(osd1 + 429);
    const auto *osd1_431 = buffer.data(osd1 + 431);
    const auto *osd1_432 = buffer.data(osd1 + 432);
    const auto *osd1_435 = buffer.data(osd1 + 435);
    const auto *osd1_437 = buffer.data(osd1 + 437);
    const auto *osd1_438 = buffer.data(osd1 + 438);
    const auto *osd1_441 = buffer.data(osd1 + 441);
    const auto *osd1_443 = buffer.data(osd1 + 443);
    const auto *osd1_444 = buffer.data(osd1 + 444);
    const auto *osd1_447 = buffer.data(osd1 + 447);
    const auto *osd1_449 = buffer.data(osd1 + 449);
    const auto *osd1_450 = buffer.data(osd1 + 450);
    const auto *osd1_453 = buffer.data(osd1 + 453);
    const auto *osd1_455 = buffer.data(osd1 + 455);
    const auto *osd1_459 = buffer.data(osd1 + 459);
    const auto *osd1_461 = buffer.data(osd1 + 461);
    const auto *osd1_462 = buffer.data(osd1 + 462);
    const auto *osd1_465 = buffer.data(osd1 + 465);
    const auto *osd1_467 = buffer.data(osd1 + 467);

    const auto *qss0_59 = buffer.data(qss0 + 59);
    const auto *qss0_60 = buffer.data(qss0 + 60);
    const auto *qss0_61 = buffer.data(qss0 + 61);
    const auto *qss0_62 = buffer.data(qss0 + 62);
    const auto *qss0_63 = buffer.data(qss0 + 63);
    const auto *qss0_64 = buffer.data(qss0 + 64);
    const auto *qss0_65 = buffer.data(qss0 + 65);
    const auto *qss0_78 = buffer.data(qss0 + 78);
    const auto *qss0_79 = buffer.data(qss0 + 79);
    const auto *qss0_80 = buffer.data(qss0 + 80);

    const auto *qss1_59 = buffer.data(qss1 + 59);
    const auto *qss1_60 = buffer.data(qss1 + 60);
    const auto *qss1_61 = buffer.data(qss1 + 61);
    const auto *qss1_62 = buffer.data(qss1 + 62);
    const auto *qss1_63 = buffer.data(qss1 + 63);
    const auto *qss1_64 = buffer.data(qss1 + 64);
    const auto *qss1_65 = buffer.data(qss1 + 65);
    const auto *qss1_78 = buffer.data(qss1 + 78);
    const auto *qss1_79 = buffer.data(qss1 + 79);
    const auto *qss1_80 = buffer.data(qss1 + 80);

    const auto *qsp_179 = buffer.data(qsp + 179);
    const auto *qsp_180 = buffer.data(qsp + 180);
    const auto *qsp_181 = buffer.data(qsp + 181);
    const auto *qsp_182 = buffer.data(qsp + 182);
    const auto *qsp_183 = buffer.data(qsp + 183);
    const auto *qsp_184 = buffer.data(qsp + 184);
    const auto *qsp_185 = buffer.data(qsp + 185);
    const auto *qsp_186 = buffer.data(qsp + 186);
    const auto *qsp_187 = buffer.data(qsp + 187);
    const auto *qsp_188 = buffer.data(qsp + 188);
    const auto *qsp_189 = buffer.data(qsp + 189);
    const auto *qsp_190 = buffer.data(qsp + 190);
    const auto *qsp_191 = buffer.data(qsp + 191);
    const auto *qsp_193 = buffer.data(qsp + 193);
    const auto *qsp_194 = buffer.data(qsp + 194);
    const auto *qsp_195 = buffer.data(qsp + 195);
    const auto *qsp_196 = buffer.data(qsp + 196);
    const auto *qsp_197 = buffer.data(qsp + 197);
    const auto *qsp_198 = buffer.data(qsp + 198);
    const auto *qsp_199 = buffer.data(qsp + 199);
    const auto *qsp_202 = buffer.data(qsp + 202);
    const auto *qsp_203 = buffer.data(qsp + 203);
    const auto *qsp_205 = buffer.data(qsp + 205);
    const auto *qsp_206 = buffer.data(qsp + 206);
    const auto *qsp_208 = buffer.data(qsp + 208);
    const auto *qsp_209 = buffer.data(qsp + 209);
    const auto *qsp_211 = buffer.data(qsp + 211);
    const auto *qsp_212 = buffer.data(qsp + 212);
    const auto *qsp_214 = buffer.data(qsp + 214);
    const auto *qsp_215 = buffer.data(qsp + 215);
    const auto *qsp_217 = buffer.data(qsp + 217);
    const auto *qsp_218 = buffer.data(qsp + 218);
    const auto *qsp_220 = buffer.data(qsp + 220);
    const auto *qsp_221 = buffer.data(qsp + 221);
    const auto *qsp_223 = buffer.data(qsp + 223);
    const auto *qsp_224 = buffer.data(qsp + 224);
    const auto *qsp_226 = buffer.data(qsp + 226);
    const auto *qsp_227 = buffer.data(qsp + 227);
    const auto *qsp_229 = buffer.data(qsp + 229);
    const auto *qsp_230 = buffer.data(qsp + 230);
    const auto *qsp_231 = buffer.data(qsp + 231);
    const auto *qsp_233 = buffer.data(qsp + 233);
    const auto *qsp_234 = buffer.data(qsp + 234);
    const auto *qsp_235 = buffer.data(qsp + 235);
    const auto *qsp_236 = buffer.data(qsp + 236);
    const auto *qsp_238 = buffer.data(qsp + 238);
    const auto *qsp_239 = buffer.data(qsp + 239);
    const auto *qsp_240 = buffer.data(qsp + 240);
    const auto *qsp_241 = buffer.data(qsp + 241);
    const auto *qsp_242 = buffer.data(qsp + 242);

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_z, osp_146, osp_180, osp_181, qss0_59, \
                         qss0_60, qss1_59, qss1_60, qsp_179, qsp_180, \
                         qsp_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * osp_146[k]
                   + f_1 * qss0_59[k]
                   - f_2 * qss1_59[k]
                   + f_3 * pc_z[k] * qsp_179[k];

        t_360[k] = f_8 * osp_180[k]
                   + f_1 * qss0_60[k]
                   - f_2 * qss1_60[k]
                   + f_3 * pc_x[k] * qsp_180[k];

        t_361[k] = f_8 * osp_181[k]
                   + f_3 * pc_x[k] * qsp_181[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pc_x, pc_y, pc_z, osp_149, osp_151, \
                         osp_152, osp_182, qss0_60, qss1_60, qsp_181, \
                         qsp_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_8 * osp_182[k]
                   + f_3 * pc_x[k] * qsp_182[k];

        t_363[k] = f_14 * osp_151[k]
                   + f_1 * qss0_60[k]
                   - f_2 * qss1_60[k]
                   + f_3 * pc_y[k] * qsp_181[k];

        t_364[k] = f_14 * osp_152[k]
                   + f_3 * pc_y[k] * qsp_182[k];

        t_365[k] = f_14 * osp_149[k]
                   + f_1 * qss0_60[k]
                   - f_2 * qss1_60[k]
                   + f_3 * pc_z[k] * qsp_182[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, pc_y, osp_154, osp_183, osp_184, \
                         osp_185, qss0_61, qss1_61, qsp_183, qsp_184, \
                         qsp_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_8 * osp_183[k]
                   + f_1 * qss0_61[k]
                   - f_2 * qss1_61[k]
                   + f_3 * pc_x[k] * qsp_183[k];

        t_367[k] = f_8 * osp_184[k]
                   + f_3 * pc_x[k] * qsp_184[k];

        t_368[k] = f_8 * osp_185[k]
                   + f_3 * pc_x[k] * qsp_185[k];

        t_369[k] = f_12 * osp_154[k]
                   + f_1 * qss0_61[k]
                   - f_2 * qss1_61[k]
                   + f_3 * pc_y[k] * qsp_184[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_y, pc_z, osp_152, osp_155, osp_186, \
                         qss0_61, qss0_62, qss1_61, qss1_62, qsp_185, \
                         qsp_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_12 * osp_155[k]
                   + f_3 * pc_y[k] * qsp_185[k];

        t_371[k] = f_15 * osp_152[k]
                   + f_1 * qss0_61[k]
                   - f_2 * qss1_61[k]
                   + f_3 * pc_z[k] * qsp_185[k];

        t_372[k] = f_8 * osp_186[k]
                   + f_1 * qss0_62[k]
                   - f_2 * qss1_62[k]
                   + f_3 * pc_x[k] * qsp_186[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pc_x, pc_y, osp_157, osp_158, osp_187, \
                         osp_188, qss0_62, qss1_62, qsp_187, qsp_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_8 * osp_187[k]
                   + f_3 * pc_x[k] * qsp_187[k];

        t_374[k] = f_8 * osp_188[k]
                   + f_3 * pc_x[k] * qsp_188[k];

        t_375[k] = f_10 * osp_157[k]
                   + f_1 * qss0_62[k]
                   - f_2 * qss1_62[k]
                   + f_3 * pc_y[k] * qsp_187[k];

        t_376[k] = f_10 * osp_158[k]
                   + f_3 * pc_y[k] * qsp_188[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pc_x, pc_z, osp_155, osp_189, osp_190, qss0_62, \
                         qss0_63, qss1_62, qss1_63, qsp_188, qsp_189, \
                         qsp_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_13 * osp_155[k]
                   + f_1 * qss0_62[k]
                   - f_2 * qss1_62[k]
                   + f_3 * pc_z[k] * qsp_188[k];

        t_378[k] = f_8 * osp_189[k]
                   + f_1 * qss0_63[k]
                   - f_2 * qss1_63[k]
                   + f_3 * pc_x[k] * qsp_189[k];

        t_379[k] = f_8 * osp_190[k]
                   + f_3 * pc_x[k] * qsp_190[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_y, pc_z, osp_158, osp_160, \
                         osp_161, osp_191, qss0_63, qss1_63, qsp_190, \
                         qsp_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_8 * osp_191[k]
                   + f_3 * pc_x[k] * qsp_191[k];

        t_381[k] = f_8 * osp_160[k]
                   + f_1 * qss0_63[k]
                   - f_2 * qss1_63[k]
                   + f_3 * pc_y[k] * qsp_190[k];

        t_382[k] = f_8 * osp_161[k]
                   + f_3 * pc_y[k] * qsp_191[k];

        t_383[k] = f_11 * osp_158[k]
                   + f_1 * qss0_63[k]
                   - f_2 * qss1_63[k]
                   + f_3 * pc_z[k] * qsp_191[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pa_y, pc_x, pc_y, osd0_324, osp_163, \
                         osp_193, osp_194, osd1_324, qss0_64, qss1_64, qsp_193, \
                         qsp_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_y[k] * osd0_324[k]
                   - f_4 * pc_y[k] * osd1_324[k];

        t_385[k] = f_8 * osp_193[k]
                   + f_3 * pc_x[k] * qsp_193[k];

        t_386[k] = f_8 * osp_194[k]
                   + f_3 * pc_x[k] * qsp_194[k];

        t_387[k] = f_6 * osp_163[k]
                   + f_1 * qss0_64[k]
                   - f_2 * qss1_64[k]
                   + f_3 * pc_y[k] * qsp_193[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pa_y, pc_x, pc_y, osd0_329, osp_164, \
                         osp_195, osd1_329, qss0_65, qss1_65, qsp_194, \
                         qsp_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_6 * osp_164[k]
                   + f_3 * pc_y[k] * qsp_194[k];

        t_389[k] = pa_y[k] * osd0_329[k]
                   - f_4 * pc_y[k] * osd1_329[k];

        t_390[k] = f_8 * osp_195[k]
                   + f_1 * qss0_65[k]
                   - f_2 * qss1_65[k]
                   + f_3 * pc_x[k] * qsp_195[k];

        t_391[k] = f_3 * pc_y[k] * qsp_195[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pc_x, pc_y, pc_z, osp_164, osp_197, \
                         qss0_65, qss1_65, qsp_196, qsp_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_8 * osp_197[k]
                   + f_3 * pc_x[k] * qsp_197[k];

        t_393[k] = f_1 * qss0_65[k]
                   - f_2 * qss1_65[k]
                   + f_3 * pc_y[k] * qsp_196[k];

        t_394[k] = f_3 * pc_y[k] * qsp_197[k];

        t_395[k] = f_7 * osp_164[k]
                   + f_1 * qss0_65[k]
                   - f_2 * qss1_65[k]
                   + f_3 * pc_z[k] * qsp_197[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, pa_x, pc_x, pc_z, osd0_396, \
                         osd0_399, osp_198, osp_199, osd1_396, osd1_399, qsp_198, \
                         qsp_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pa_x[k] * osd0_396[k]
                   + f_8 * osp_198[k]
                   - f_4 * pc_x[k] * osd1_396[k];

        t_397[k] = f_6 * osp_199[k]
                   + f_3 * pc_x[k] * qsp_199[k];

        t_398[k] = f_3 * pc_z[k] * qsp_198[k];

        t_399[k] = pa_x[k] * osd0_399[k]
                   - f_4 * pc_x[k] * osd1_399[k];

        t_400[k] = f_3 * pc_z[k] * qsp_199[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_x, pa_z, pc_x, pc_z, osd0_330, \
                         osd0_401, osp_202, osp_203, osd1_330, osd1_401, qsp_202, \
                         qsp_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = pa_x[k] * osd0_401[k]
                   - f_4 * pc_x[k] * osd1_401[k];

        t_402[k] = pa_z[k] * osd0_330[k]
                   - f_4 * pc_z[k] * osd1_330[k];

        t_403[k] = f_6 * osp_202[k]
                   + f_3 * pc_x[k] * qsp_202[k];

        t_404[k] = f_6 * osp_203[k]
                   + f_3 * pc_x[k] * qsp_203[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_x, pc_x, pc_y, osd0_405, osd0_407, \
                         osd0_408, osp_170, osp_204, osd1_405, osd1_407, osd1_408, \
                         qsp_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_x[k] * osd0_405[k]
                   - f_4 * pc_x[k] * osd1_405[k];

        t_406[k] = f_7 * osp_170[k]
                   + f_3 * pc_y[k] * qsp_203[k];

        t_407[k] = pa_x[k] * osd0_407[k]
                   - f_4 * pc_x[k] * osd1_407[k];

        t_408[k] = pa_x[k] * osd0_408[k]
                   + f_8 * osp_204[k]
                   - f_4 * pc_x[k] * osd1_408[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, pa_x, pc_x, pc_y, osd0_411, osp_173, \
                         osp_205, osp_206, osd1_411, qsp_205, qsp_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_6 * osp_205[k]
                   + f_3 * pc_x[k] * qsp_205[k];

        t_410[k] = f_6 * osp_206[k]
                   + f_3 * pc_x[k] * qsp_206[k];

        t_411[k] = pa_x[k] * osd0_411[k]
                   - f_4 * pc_x[k] * osd1_411[k];

        t_412[k] = f_9 * osp_173[k]
                   + f_3 * pc_y[k] * qsp_206[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, pa_x, pc_x, osd0_413, osd0_414, osp_207, \
                         osp_208, osp_209, osd1_413, osd1_414, qsp_208, \
                         qsp_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_x[k] * osd0_413[k]
                   - f_4 * pc_x[k] * osd1_413[k];

        t_414[k] = pa_x[k] * osd0_414[k]
                   + f_8 * osp_207[k]
                   - f_4 * pc_x[k] * osd1_414[k];

        t_415[k] = f_6 * osp_208[k]
                   + f_3 * pc_x[k] * qsp_208[k];

        t_416[k] = f_6 * osp_209[k]
                   + f_3 * pc_x[k] * qsp_209[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pa_x, pc_x, pc_y, osd0_417, osd0_419, \
                         osd0_420, osp_176, osp_210, osd1_417, osd1_419, osd1_420, \
                         qsp_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = pa_x[k] * osd0_417[k]
                   - f_4 * pc_x[k] * osd1_417[k];

        t_418[k] = f_11 * osp_176[k]
                   + f_3 * pc_y[k] * qsp_209[k];

        t_419[k] = pa_x[k] * osd0_419[k]
                   - f_4 * pc_x[k] * osd1_419[k];

        t_420[k] = pa_x[k] * osd0_420[k]
                   + f_8 * osp_210[k]
                   - f_4 * pc_x[k] * osd1_420[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pa_x, pc_x, pc_y, osd0_423, osp_179, \
                         osp_211, osp_212, osd1_423, qsp_211, qsp_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_6 * osp_211[k]
                   + f_3 * pc_x[k] * qsp_211[k];

        t_422[k] = f_6 * osp_212[k]
                   + f_3 * pc_x[k] * qsp_212[k];

        t_423[k] = pa_x[k] * osd0_423[k]
                   - f_4 * pc_x[k] * osd1_423[k];

        t_424[k] = f_13 * osp_179[k]
                   + f_3 * pc_y[k] * qsp_212[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pa_x, pc_x, osd0_425, osd0_426, osp_213, \
                         osp_214, osp_215, osd1_425, osd1_426, qsp_214, \
                         qsp_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = pa_x[k] * osd0_425[k]
                   - f_4 * pc_x[k] * osd1_425[k];

        t_426[k] = pa_x[k] * osd0_426[k]
                   + f_8 * osp_213[k]
                   - f_4 * pc_x[k] * osd1_426[k];

        t_427[k] = f_6 * osp_214[k]
                   + f_3 * pc_x[k] * qsp_214[k];

        t_428[k] = f_6 * osp_215[k]
                   + f_3 * pc_x[k] * qsp_215[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_x, pc_x, pc_y, osd0_429, osd0_431, \
                         osd0_432, osp_182, osp_216, osd1_429, osd1_431, osd1_432, \
                         qsp_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = pa_x[k] * osd0_429[k]
                   - f_4 * pc_x[k] * osd1_429[k];

        t_430[k] = f_15 * osp_182[k]
                   + f_3 * pc_y[k] * qsp_215[k];

        t_431[k] = pa_x[k] * osd0_431[k]
                   - f_4 * pc_x[k] * osd1_431[k];

        t_432[k] = pa_x[k] * osd0_432[k]
                   + f_8 * osp_216[k]
                   - f_4 * pc_x[k] * osd1_432[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pa_x, pc_x, pc_y, osd0_435, osp_185, \
                         osp_217, osp_218, osd1_435, qsp_217, qsp_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_6 * osp_217[k]
                   + f_3 * pc_x[k] * qsp_217[k];

        t_434[k] = f_6 * osp_218[k]
                   + f_3 * pc_x[k] * qsp_218[k];

        t_435[k] = pa_x[k] * osd0_435[k]
                   - f_4 * pc_x[k] * osd1_435[k];

        t_436[k] = f_14 * osp_185[k]
                   + f_3 * pc_y[k] * qsp_218[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pa_x, pc_x, osd0_437, osd0_438, osp_219, \
                         osp_220, osp_221, osd1_437, osd1_438, qsp_220, \
                         qsp_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = pa_x[k] * osd0_437[k]
                   - f_4 * pc_x[k] * osd1_437[k];

        t_438[k] = pa_x[k] * osd0_438[k]
                   + f_8 * osp_219[k]
                   - f_4 * pc_x[k] * osd1_438[k];

        t_439[k] = f_6 * osp_220[k]
                   + f_3 * pc_x[k] * qsp_220[k];

        t_440[k] = f_6 * osp_221[k]
                   + f_3 * pc_x[k] * qsp_221[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pa_x, pc_x, pc_y, osd0_441, osd0_443, \
                         osd0_444, osp_188, osp_222, osd1_441, osd1_443, osd1_444, \
                         qsp_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pa_x[k] * osd0_441[k]
                   - f_4 * pc_x[k] * osd1_441[k];

        t_442[k] = f_12 * osp_188[k]
                   + f_3 * pc_y[k] * qsp_221[k];

        t_443[k] = pa_x[k] * osd0_443[k]
                   - f_4 * pc_x[k] * osd1_443[k];

        t_444[k] = pa_x[k] * osd0_444[k]
                   + f_8 * osp_222[k]
                   - f_4 * pc_x[k] * osd1_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_x, pc_x, pc_y, osd0_447, osp_191, \
                         osp_223, osp_224, osd1_447, qsp_223, qsp_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_6 * osp_223[k]
                   + f_3 * pc_x[k] * qsp_223[k];

        t_446[k] = f_6 * osp_224[k]
                   + f_3 * pc_x[k] * qsp_224[k];

        t_447[k] = pa_x[k] * osd0_447[k]
                   - f_4 * pc_x[k] * osd1_447[k];

        t_448[k] = f_10 * osp_191[k]
                   + f_3 * pc_y[k] * qsp_224[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_x, pc_x, osd0_449, osd0_450, osp_225, \
                         osp_226, osp_227, osd1_449, osd1_450, qsp_226, \
                         qsp_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_x[k] * osd0_449[k]
                   - f_4 * pc_x[k] * osd1_449[k];

        t_450[k] = pa_x[k] * osd0_450[k]
                   + f_8 * osp_225[k]
                   - f_4 * pc_x[k] * osd1_450[k];

        t_451[k] = f_6 * osp_226[k]
                   + f_3 * pc_x[k] * qsp_226[k];

        t_452[k] = f_6 * osp_227[k]
                   + f_3 * pc_x[k] * qsp_227[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pa_x, pa_y, pc_x, pc_y, osd0_390, \
                         osd0_453, osd0_455, osp_194, osd1_390, osd1_453, osd1_455, \
                         qsp_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = pa_x[k] * osd0_453[k]
                   - f_4 * pc_x[k] * osd1_453[k];

        t_454[k] = f_8 * osp_194[k]
                   + f_3 * pc_y[k] * qsp_227[k];

        t_455[k] = pa_x[k] * osd0_455[k]
                   - f_4 * pc_x[k] * osd1_455[k];

        t_456[k] = pa_y[k] * osd0_390[k]
                   - f_4 * pc_y[k] * osd1_390[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pa_x, pc_x, pc_y, osd0_459, osp_197, \
                         osp_229, osp_230, osd1_459, qsp_229, qsp_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_6 * osp_229[k]
                   + f_3 * pc_x[k] * qsp_229[k];

        t_458[k] = f_6 * osp_230[k]
                   + f_3 * pc_x[k] * qsp_230[k];

        t_459[k] = pa_x[k] * osd0_459[k]
                   - f_4 * pc_x[k] * osd1_459[k];

        t_460[k] = f_6 * osp_197[k]
                   + f_3 * pc_y[k] * qsp_230[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_x, pc_x, pc_y, osd0_461, osd0_462, \
                         osp_231, osp_233, osd1_461, osd1_462, qsp_231, \
                         qsp_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = pa_x[k] * osd0_461[k]
                   - f_4 * pc_x[k] * osd1_461[k];

        t_462[k] = pa_x[k] * osd0_462[k]
                   + f_8 * osp_231[k]
                   - f_4 * pc_x[k] * osd1_462[k];

        t_463[k] = f_3 * pc_y[k] * qsp_231[k];

        t_464[k] = f_6 * osp_233[k]
                   + f_3 * pc_x[k] * qsp_233[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pa_x, pc_x, pc_y, osd0_465, osd0_467, \
                         osd1_465, osd1_467, qss0_78, qss1_78, qsp_233, \
                         qsp_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_x[k] * osd0_465[k]
                   - f_4 * pc_x[k] * osd1_465[k];

        t_466[k] = f_3 * pc_y[k] * qsp_233[k];

        t_467[k] = pa_x[k] * osd0_467[k]
                   - f_4 * pc_x[k] * osd1_467[k];

        t_468[k] = f_1 * qss0_78[k]
                   - f_2 * qss1_78[k]
                   + f_3 * pc_x[k] * qsp_234[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, pc_x, pc_y, pc_z, osp_199, \
                         qss0_78, qss1_78, qsp_235, qsp_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_3 * pc_x[k] * qsp_235[k];

        t_470[k] = f_3 * pc_x[k] * qsp_236[k];

        t_471[k] = f_0 * osp_199[k]
                   + f_1 * qss0_78[k]
                   - f_2 * qss1_78[k]
                   + f_3 * pc_y[k] * qsp_235[k];

        t_472[k] = f_3 * pc_z[k] * qsp_235[k];

        t_473[k] = f_1 * qss0_78[k]
                   - f_2 * qss1_78[k]
                   + f_3 * pc_z[k] * qsp_236[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, t_478, pa_z, pc_x, pc_y, pc_z, osd0_396, \
                         osd0_399, osp_203, osd1_396, osd1_399, qsp_238, \
                         qsp_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_z[k] * osd0_396[k]
                   - f_4 * pc_z[k] * osd1_396[k];

        t_475[k] = f_3 * pc_x[k] * qsp_238[k];

        t_476[k] = f_3 * pc_x[k] * qsp_239[k];

        t_477[k] = pa_z[k] * osd0_399[k]
                   - f_4 * pc_z[k] * osd1_399[k];

        t_478[k] = f_5 * osp_203[k]
                   + f_3 * pc_y[k] * qsp_239[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pc_x, pc_z, osp_200, qss0_79, qss0_80, \
                         qss1_79, qss1_80, qsp_239, qsp_240, qsp_241, \
                         qsp_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_6 * osp_200[k]
                   + f_1 * qss0_79[k]
                   - f_2 * qss1_79[k]
                   + f_3 * pc_z[k] * qsp_239[k];

        t_480[k] = f_1 * qss0_80[k]
                   - f_2 * qss1_80[k]
                   + f_3 * pc_x[k] * qsp_240[k];

        t_481[k] = f_3 * pc_x[k] * qsp_241[k];

        t_482[k] = f_3 * pc_x[k] * qsp_242[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pc_y, pc_z, osp_203, osp_205, osp_206, qss0_80, \
                         qss1_80, qsp_241, qsp_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_7 * osp_205[k]
                   + f_1 * qss0_80[k]
                   - f_2 * qss1_80[k]
                   + f_3 * pc_y[k] * qsp_241[k];

        t_484[k] = f_7 * osp_206[k]
                   + f_3 * pc_y[k] * qsp_242[k];

        t_485[k] = f_8 * osp_203[k]
                   + f_1 * qss0_80[k]
                   - f_2 * qss1_80[k]
                   + f_3 * pc_z[k] * qsp_242[k];
    }
}

static auto
compute_prim_qsd_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osd0,
                                                          const size_t osp, const size_t osd1,
                                                          const size_t qss0, const size_t qss1,
                                                          const size_t qsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 5.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osd0_462 = buffer.data(osd0 + 462);
    const auto *osd0_465 = buffer.data(osd0 + 465);
    const auto *osd0_467 = buffer.data(osd0 + 467);

    const auto *osp_206 = buffer.data(osp + 206);
    const auto *osp_208 = buffer.data(osp + 208);
    const auto *osp_209 = buffer.data(osp + 209);
    const auto *osp_211 = buffer.data(osp + 211);
    const auto *osp_212 = buffer.data(osp + 212);
    const auto *osp_214 = buffer.data(osp + 214);
    const auto *osp_215 = buffer.data(osp + 215);
    const auto *osp_217 = buffer.data(osp + 217);
    const auto *osp_218 = buffer.data(osp + 218);
    const auto *osp_220 = buffer.data(osp + 220);
    const auto *osp_221 = buffer.data(osp + 221);
    const auto *osp_223 = buffer.data(osp + 223);
    const auto *osp_224 = buffer.data(osp + 224);
    const auto *osp_226 = buffer.data(osp + 226);
    const auto *osp_227 = buffer.data(osp + 227);
    const auto *osp_229 = buffer.data(osp + 229);
    const auto *osp_230 = buffer.data(osp + 230);
    const auto *osp_232 = buffer.data(osp + 232);
    const auto *osp_233 = buffer.data(osp + 233);

    const auto *osd1_462 = buffer.data(osd1 + 462);
    const auto *osd1_465 = buffer.data(osd1 + 465);
    const auto *osd1_467 = buffer.data(osd1 + 467);

    const auto *qss0_81 = buffer.data(qss0 + 81);
    const auto *qss0_82 = buffer.data(qss0 + 82);
    const auto *qss0_83 = buffer.data(qss0 + 83);
    const auto *qss0_84 = buffer.data(qss0 + 84);
    const auto *qss0_85 = buffer.data(qss0 + 85);
    const auto *qss0_86 = buffer.data(qss0 + 86);
    const auto *qss0_87 = buffer.data(qss0 + 87);
    const auto *qss0_88 = buffer.data(qss0 + 88);
    const auto *qss0_90 = buffer.data(qss0 + 90);

    const auto *qss1_81 = buffer.data(qss1 + 81);
    const auto *qss1_82 = buffer.data(qss1 + 82);
    const auto *qss1_83 = buffer.data(qss1 + 83);
    const auto *qss1_84 = buffer.data(qss1 + 84);
    const auto *qss1_85 = buffer.data(qss1 + 85);
    const auto *qss1_86 = buffer.data(qss1 + 86);
    const auto *qss1_87 = buffer.data(qss1 + 87);
    const auto *qss1_88 = buffer.data(qss1 + 88);
    const auto *qss1_90 = buffer.data(qss1 + 90);

    const auto *qsp_243 = buffer.data(qsp + 243);
    const auto *qsp_244 = buffer.data(qsp + 244);
    const auto *qsp_245 = buffer.data(qsp + 245);
    const auto *qsp_246 = buffer.data(qsp + 246);
    const auto *qsp_247 = buffer.data(qsp + 247);
    const auto *qsp_248 = buffer.data(qsp + 248);
    const auto *qsp_249 = buffer.data(qsp + 249);
    const auto *qsp_250 = buffer.data(qsp + 250);
    const auto *qsp_251 = buffer.data(qsp + 251);
    const auto *qsp_252 = buffer.data(qsp + 252);
    const auto *qsp_253 = buffer.data(qsp + 253);
    const auto *qsp_254 = buffer.data(qsp + 254);
    const auto *qsp_255 = buffer.data(qsp + 255);
    const auto *qsp_256 = buffer.data(qsp + 256);
    const auto *qsp_257 = buffer.data(qsp + 257);
    const auto *qsp_258 = buffer.data(qsp + 258);
    const auto *qsp_259 = buffer.data(qsp + 259);
    const auto *qsp_260 = buffer.data(qsp + 260);
    const auto *qsp_261 = buffer.data(qsp + 261);
    const auto *qsp_262 = buffer.data(qsp + 262);
    const auto *qsp_263 = buffer.data(qsp + 263);
    const auto *qsp_264 = buffer.data(qsp + 264);
    const auto *qsp_265 = buffer.data(qsp + 265);
    const auto *qsp_266 = buffer.data(qsp + 266);
    const auto *qsp_268 = buffer.data(qsp + 268);
    const auto *qsp_269 = buffer.data(qsp + 269);
    const auto *qsp_270 = buffer.data(qsp + 270);
    const auto *qsp_271 = buffer.data(qsp + 271);
    const auto *qsp_272 = buffer.data(qsp + 272);

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, pc_x, pc_y, osp_208, osp_209, \
                         qss0_81, qss1_81, qsp_243, qsp_244, qsp_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_1 * qss0_81[k]
                   - f_2 * qss1_81[k]
                   + f_3 * pc_x[k] * qsp_243[k];

        t_487[k] = f_3 * pc_x[k] * qsp_244[k];

        t_488[k] = f_3 * pc_x[k] * qsp_245[k];

        t_489[k] = f_9 * osp_208[k]
                   + f_1 * qss0_81[k]
                   - f_2 * qss1_81[k]
                   + f_3 * pc_y[k] * qsp_244[k];

        t_490[k] = f_9 * osp_209[k]
                   + f_3 * pc_y[k] * qsp_245[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_z, osp_206, qss0_81, qss0_82, \
                         qss1_81, qss1_82, qsp_245, qsp_246, qsp_247, \
                         qsp_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * osp_206[k]
                   + f_1 * qss0_81[k]
                   - f_2 * qss1_81[k]
                   + f_3 * pc_z[k] * qsp_245[k];

        t_492[k] = f_1 * qss0_82[k]
                   - f_2 * qss1_82[k]
                   + f_3 * pc_x[k] * qsp_246[k];

        t_493[k] = f_3 * pc_x[k] * qsp_247[k];

        t_494[k] = f_3 * pc_x[k] * qsp_248[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pc_y, pc_z, osp_209, osp_211, osp_212, qss0_82, \
                         qss1_82, qsp_247, qsp_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_11 * osp_211[k]
                   + f_1 * qss0_82[k]
                   - f_2 * qss1_82[k]
                   + f_3 * pc_y[k] * qsp_247[k];

        t_496[k] = f_11 * osp_212[k]
                   + f_3 * pc_y[k] * qsp_248[k];

        t_497[k] = f_12 * osp_209[k]
                   + f_1 * qss0_82[k]
                   - f_2 * qss1_82[k]
                   + f_3 * pc_z[k] * qsp_248[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, t_502, pc_x, pc_y, osp_214, osp_215, \
                         qss0_83, qss1_83, qsp_249, qsp_250, qsp_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_1 * qss0_83[k]
                   - f_2 * qss1_83[k]
                   + f_3 * pc_x[k] * qsp_249[k];

        t_499[k] = f_3 * pc_x[k] * qsp_250[k];

        t_500[k] = f_3 * pc_x[k] * qsp_251[k];

        t_501[k] = f_13 * osp_214[k]
                   + f_1 * qss0_83[k]
                   - f_2 * qss1_83[k]
                   + f_3 * pc_y[k] * qsp_250[k];

        t_502[k] = f_13 * osp_215[k]
                   + f_3 * pc_y[k] * qsp_251[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, pc_x, pc_z, osp_212, qss0_83, qss0_84, \
                         qss1_83, qss1_84, qsp_251, qsp_252, qsp_253, \
                         qsp_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_14 * osp_212[k]
                   + f_1 * qss0_83[k]
                   - f_2 * qss1_83[k]
                   + f_3 * pc_z[k] * qsp_251[k];

        t_504[k] = f_1 * qss0_84[k]
                   - f_2 * qss1_84[k]
                   + f_3 * pc_x[k] * qsp_252[k];

        t_505[k] = f_3 * pc_x[k] * qsp_253[k];

        t_506[k] = f_3 * pc_x[k] * qsp_254[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pc_y, pc_z, osp_215, osp_217, osp_218, qss0_84, \
                         qss1_84, qsp_253, qsp_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_15 * osp_217[k]
                   + f_1 * qss0_84[k]
                   - f_2 * qss1_84[k]
                   + f_3 * pc_y[k] * qsp_253[k];

        t_508[k] = f_15 * osp_218[k]
                   + f_3 * pc_y[k] * qsp_254[k];

        t_509[k] = f_15 * osp_215[k]
                   + f_1 * qss0_84[k]
                   - f_2 * qss1_84[k]
                   + f_3 * pc_z[k] * qsp_254[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, pc_x, pc_y, osp_220, osp_221, \
                         qss0_85, qss1_85, qsp_255, qsp_256, qsp_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_1 * qss0_85[k]
                   - f_2 * qss1_85[k]
                   + f_3 * pc_x[k] * qsp_255[k];

        t_511[k] = f_3 * pc_x[k] * qsp_256[k];

        t_512[k] = f_3 * pc_x[k] * qsp_257[k];

        t_513[k] = f_14 * osp_220[k]
                   + f_1 * qss0_85[k]
                   - f_2 * qss1_85[k]
                   + f_3 * pc_y[k] * qsp_256[k];

        t_514[k] = f_14 * osp_221[k]
                   + f_3 * pc_y[k] * qsp_257[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pc_x, pc_z, osp_218, qss0_85, qss0_86, \
                         qss1_85, qss1_86, qsp_257, qsp_258, qsp_259, \
                         qsp_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * osp_218[k]
                   + f_1 * qss0_85[k]
                   - f_2 * qss1_85[k]
                   + f_3 * pc_z[k] * qsp_257[k];

        t_516[k] = f_1 * qss0_86[k]
                   - f_2 * qss1_86[k]
                   + f_3 * pc_x[k] * qsp_258[k];

        t_517[k] = f_3 * pc_x[k] * qsp_259[k];

        t_518[k] = f_3 * pc_x[k] * qsp_260[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, pc_z, osp_221, osp_223, osp_224, qss0_86, \
                         qss1_86, qsp_259, qsp_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_12 * osp_223[k]
                   + f_1 * qss0_86[k]
                   - f_2 * qss1_86[k]
                   + f_3 * pc_y[k] * qsp_259[k];

        t_520[k] = f_12 * osp_224[k]
                   + f_3 * pc_y[k] * qsp_260[k];

        t_521[k] = f_11 * osp_221[k]
                   + f_1 * qss0_86[k]
                   - f_2 * qss1_86[k]
                   + f_3 * pc_z[k] * qsp_260[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pc_x, pc_y, osp_226, osp_227, \
                         qss0_87, qss1_87, qsp_261, qsp_262, qsp_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_1 * qss0_87[k]
                   - f_2 * qss1_87[k]
                   + f_3 * pc_x[k] * qsp_261[k];

        t_523[k] = f_3 * pc_x[k] * qsp_262[k];

        t_524[k] = f_3 * pc_x[k] * qsp_263[k];

        t_525[k] = f_10 * osp_226[k]
                   + f_1 * qss0_87[k]
                   - f_2 * qss1_87[k]
                   + f_3 * pc_y[k] * qsp_262[k];

        t_526[k] = f_10 * osp_227[k]
                   + f_3 * pc_y[k] * qsp_263[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, pc_x, pc_z, osp_224, qss0_87, qss0_88, \
                         qss1_87, qss1_88, qsp_263, qsp_264, qsp_265, \
                         qsp_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_9 * osp_224[k]
                   + f_1 * qss0_87[k]
                   - f_2 * qss1_87[k]
                   + f_3 * pc_z[k] * qsp_263[k];

        t_528[k] = f_1 * qss0_88[k]
                   - f_2 * qss1_88[k]
                   + f_3 * pc_x[k] * qsp_264[k];

        t_529[k] = f_3 * pc_x[k] * qsp_265[k];

        t_530[k] = f_3 * pc_x[k] * qsp_266[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, pa_y, pc_y, pc_z, osd0_462, osp_227, \
                         osp_229, osp_230, osd1_462, qss0_88, qss1_88, qsp_265, \
                         qsp_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_8 * osp_229[k]
                   + f_1 * qss0_88[k]
                   - f_2 * qss1_88[k]
                   + f_3 * pc_y[k] * qsp_265[k];

        t_532[k] = f_8 * osp_230[k]
                   + f_3 * pc_y[k] * qsp_266[k];

        t_533[k] = f_7 * osp_227[k]
                   + f_1 * qss0_88[k]
                   - f_2 * qss1_88[k]
                   + f_3 * pc_z[k] * qsp_266[k];

        t_534[k] = pa_y[k] * osd0_462[k]
                   - f_4 * pc_y[k] * osd1_462[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, pa_y, pc_x, pc_y, osd0_465, \
                         osd0_467, osp_232, osp_233, osd1_465, osd1_467, qsp_268, \
                         qsp_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_3 * pc_x[k] * qsp_268[k];

        t_536[k] = f_3 * pc_x[k] * qsp_269[k];

        t_537[k] = pa_y[k] * osd0_465[k]
                   + f_8 * osp_232[k]
                   - f_4 * pc_y[k] * osd1_465[k];

        t_538[k] = f_6 * osp_233[k]
                   + f_3 * pc_y[k] * qsp_269[k];

        t_539[k] = pa_y[k] * osd0_467[k]
                   - f_4 * pc_y[k] * osd1_467[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, t_545, pc_x, pc_y, pc_z, osp_233, \
                         qss0_90, qss1_90, qsp_270, qsp_271, qsp_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * qss0_90[k]
                   - f_2 * qss1_90[k]
                   + f_3 * pc_x[k] * qsp_270[k];

        t_541[k] = f_3 * pc_x[k] * qsp_271[k];

        t_542[k] = f_3 * pc_x[k] * qsp_272[k];

        t_543[k] = f_1 * qss0_90[k]
                   - f_2 * qss1_90[k]
                   + f_3 * pc_y[k] * qsp_271[k];

        t_544[k] = f_3 * pc_y[k] * qsp_272[k];

        t_545[k] = f_0 * osp_233[k]
                   + f_1 * qss0_90[k]
                   - f_2 * qss1_90[k]
                   + f_3 * pc_z[k] * qsp_272[k];
    }
}

auto
compute_prim_qsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t osd0, const size_t osp,
                                                   const size_t osd1, const size_t qss0,
                                                   const size_t qss1, const size_t qsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_qsd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, osd0, osp,
                                                              osd1, qss0, qss1, qsp, ncols,
                                                              gamma, p, q);

    compute_prim_qsd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, osd0, osp,
                                                              osd1, qss0, qss1, qsp, ncols,
                                                              gamma, p, q);

    compute_prim_qsd_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, osd0, osp,
                                                              osd1, qss0, qss1, qsp, ncols,
                                                              gamma, p, q);

    compute_prim_qsd_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, osd0, osp,
                                                              osd1, qss0, qss1, qsp, ncols,
                                                              gamma, p, q);

    compute_prim_qsd_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, osd0, osp,
                                                              osd1, qss0, qss1, qsp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
