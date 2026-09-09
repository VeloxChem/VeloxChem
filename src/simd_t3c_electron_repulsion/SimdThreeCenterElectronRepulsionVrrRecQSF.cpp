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


#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_0 = buffer.data(osf0 + 0);
    const auto *osf0_6 = buffer.data(osf0 + 6);
    const auto *osf0_9 = buffer.data(osf0 + 9);
    const auto *osf0_16 = buffer.data(osf0 + 16);
    const auto *osf0_20 = buffer.data(osf0 + 20);
    const auto *osf0_29 = buffer.data(osf0 + 29);
    const auto *osf0_30 = buffer.data(osf0 + 30);
    const auto *osf0_36 = buffer.data(osf0 + 36);
    const auto *osf0_50 = buffer.data(osf0 + 50);
    const auto *osf0_59 = buffer.data(osf0 + 59);
    const auto *osf0_60 = buffer.data(osf0 + 60);
    const auto *osf0_66 = buffer.data(osf0 + 66);
    const auto *osf0_90 = buffer.data(osf0 + 90);

    const auto *osd_0 = buffer.data(osd + 0);
    const auto *osd_3 = buffer.data(osd + 3);
    const auto *osd_5 = buffer.data(osd + 5);
    const auto *osd_6 = buffer.data(osd + 6);
    const auto *osd_9 = buffer.data(osd + 9);
    const auto *osd_11 = buffer.data(osd + 11);
    const auto *osd_12 = buffer.data(osd + 12);
    const auto *osd_15 = buffer.data(osd + 15);
    const auto *osd_17 = buffer.data(osd + 17);
    const auto *osd_18 = buffer.data(osd + 18);
    const auto *osd_21 = buffer.data(osd + 21);
    const auto *osd_23 = buffer.data(osd + 23);
    const auto *osd_24 = buffer.data(osd + 24);
    const auto *osd_27 = buffer.data(osd + 27);
    const auto *osd_28 = buffer.data(osd + 28);
    const auto *osd_29 = buffer.data(osd + 29);
    const auto *osd_30 = buffer.data(osd + 30);
    const auto *osd_33 = buffer.data(osd + 33);
    const auto *osd_35 = buffer.data(osd + 35);
    const auto *osd_36 = buffer.data(osd + 36);
    const auto *osd_39 = buffer.data(osd + 39);
    const auto *osd_41 = buffer.data(osd + 41);
    const auto *osd_42 = buffer.data(osd + 42);
    const auto *osd_45 = buffer.data(osd + 45);
    const auto *osd_46 = buffer.data(osd + 46);
    const auto *osd_47 = buffer.data(osd + 47);
    const auto *osd_48 = buffer.data(osd + 48);
    const auto *osd_51 = buffer.data(osd + 51);
    const auto *osd_52 = buffer.data(osd + 52);
    const auto *osd_53 = buffer.data(osd + 53);
    const auto *osd_54 = buffer.data(osd + 54);
    const auto *osd_57 = buffer.data(osd + 57);
    const auto *osd_59 = buffer.data(osd + 59);
    const auto *osd_60 = buffer.data(osd + 60);
    const auto *osd_63 = buffer.data(osd + 63);
    const auto *osd_65 = buffer.data(osd + 65);
    const auto *osd_69 = buffer.data(osd + 69);
    const auto *osd_70 = buffer.data(osd + 70);
    const auto *osd_71 = buffer.data(osd + 71);
    const auto *osd_72 = buffer.data(osd + 72);
    const auto *osd_75 = buffer.data(osd + 75);
    const auto *osd_76 = buffer.data(osd + 76);
    const auto *osd_77 = buffer.data(osd + 77);

    const auto *osf1_0 = buffer.data(osf1 + 0);
    const auto *osf1_6 = buffer.data(osf1 + 6);
    const auto *osf1_9 = buffer.data(osf1 + 9);
    const auto *osf1_16 = buffer.data(osf1 + 16);
    const auto *osf1_20 = buffer.data(osf1 + 20);
    const auto *osf1_29 = buffer.data(osf1 + 29);
    const auto *osf1_30 = buffer.data(osf1 + 30);
    const auto *osf1_36 = buffer.data(osf1 + 36);
    const auto *osf1_50 = buffer.data(osf1 + 50);
    const auto *osf1_59 = buffer.data(osf1 + 59);
    const auto *osf1_60 = buffer.data(osf1 + 60);
    const auto *osf1_66 = buffer.data(osf1 + 66);
    const auto *osf1_90 = buffer.data(osf1 + 90);

    const auto *qsp0_0 = buffer.data(qsp0 + 0);
    const auto *qsp0_1 = buffer.data(qsp0 + 1);
    const auto *qsp0_2 = buffer.data(qsp0 + 2);
    const auto *qsp0_4 = buffer.data(qsp0 + 4);
    const auto *qsp0_8 = buffer.data(qsp0 + 8);
    const auto *qsp0_9 = buffer.data(qsp0 + 9);
    const auto *qsp0_10 = buffer.data(qsp0 + 10);
    const auto *qsp0_11 = buffer.data(qsp0 + 11);
    const auto *qsp0_15 = buffer.data(qsp0 + 15);
    const auto *qsp0_16 = buffer.data(qsp0 + 16);
    const auto *qsp0_17 = buffer.data(qsp0 + 17);
    const auto *qsp0_18 = buffer.data(qsp0 + 18);
    const auto *qsp0_19 = buffer.data(qsp0 + 19);
    const auto *qsp0_20 = buffer.data(qsp0 + 20);
    const auto *qsp0_23 = buffer.data(qsp0 + 23);
    const auto *qsp0_25 = buffer.data(qsp0 + 25);
    const auto *qsp0_27 = buffer.data(qsp0 + 27);
    const auto *qsp0_28 = buffer.data(qsp0 + 28);
    const auto *qsp0_29 = buffer.data(qsp0 + 29);
    const auto *qsp0_30 = buffer.data(qsp0 + 30);
    const auto *qsp0_31 = buffer.data(qsp0 + 31);
    const auto *qsp0_32 = buffer.data(qsp0 + 32);
    const auto *qsp0_35 = buffer.data(qsp0 + 35);
    const auto *qsp0_36 = buffer.data(qsp0 + 36);
    const auto *qsp0_37 = buffer.data(qsp0 + 37);
    const auto *qsp0_38 = buffer.data(qsp0 + 38);

    const auto *qsp1_0 = buffer.data(qsp1 + 0);
    const auto *qsp1_1 = buffer.data(qsp1 + 1);
    const auto *qsp1_2 = buffer.data(qsp1 + 2);
    const auto *qsp1_4 = buffer.data(qsp1 + 4);
    const auto *qsp1_8 = buffer.data(qsp1 + 8);
    const auto *qsp1_9 = buffer.data(qsp1 + 9);
    const auto *qsp1_10 = buffer.data(qsp1 + 10);
    const auto *qsp1_11 = buffer.data(qsp1 + 11);
    const auto *qsp1_15 = buffer.data(qsp1 + 15);
    const auto *qsp1_16 = buffer.data(qsp1 + 16);
    const auto *qsp1_17 = buffer.data(qsp1 + 17);
    const auto *qsp1_18 = buffer.data(qsp1 + 18);
    const auto *qsp1_19 = buffer.data(qsp1 + 19);
    const auto *qsp1_20 = buffer.data(qsp1 + 20);
    const auto *qsp1_23 = buffer.data(qsp1 + 23);
    const auto *qsp1_25 = buffer.data(qsp1 + 25);
    const auto *qsp1_27 = buffer.data(qsp1 + 27);
    const auto *qsp1_28 = buffer.data(qsp1 + 28);
    const auto *qsp1_29 = buffer.data(qsp1 + 29);
    const auto *qsp1_30 = buffer.data(qsp1 + 30);
    const auto *qsp1_31 = buffer.data(qsp1 + 31);
    const auto *qsp1_32 = buffer.data(qsp1 + 32);
    const auto *qsp1_35 = buffer.data(qsp1 + 35);
    const auto *qsp1_36 = buffer.data(qsp1 + 36);
    const auto *qsp1_37 = buffer.data(qsp1 + 37);
    const auto *qsp1_38 = buffer.data(qsp1 + 38);

    const auto *qsd_0 = buffer.data(qsd + 0);
    const auto *qsd_2 = buffer.data(qsd + 2);
    const auto *qsd_3 = buffer.data(qsd + 3);
    const auto *qsd_5 = buffer.data(qsd + 5);
    const auto *qsd_6 = buffer.data(qsd + 6);
    const auto *qsd_7 = buffer.data(qsd + 7);
    const auto *qsd_9 = buffer.data(qsd + 9);
    const auto *qsd_11 = buffer.data(qsd + 11);
    const auto *qsd_12 = buffer.data(qsd + 12);
    const auto *qsd_14 = buffer.data(qsd + 14);
    const auto *qsd_15 = buffer.data(qsd + 15);
    const auto *qsd_16 = buffer.data(qsd + 16);
    const auto *qsd_17 = buffer.data(qsd + 17);
    const auto *qsd_18 = buffer.data(qsd + 18);
    const auto *qsd_19 = buffer.data(qsd + 19);
    const auto *qsd_21 = buffer.data(qsd + 21);
    const auto *qsd_23 = buffer.data(qsd + 23);
    const auto *qsd_24 = buffer.data(qsd + 24);
    const auto *qsd_27 = buffer.data(qsd + 27);
    const auto *qsd_28 = buffer.data(qsd + 28);
    const auto *qsd_29 = buffer.data(qsd + 29);
    const auto *qsd_30 = buffer.data(qsd + 30);
    const auto *qsd_32 = buffer.data(qsd + 32);
    const auto *qsd_33 = buffer.data(qsd + 33);
    const auto *qsd_34 = buffer.data(qsd + 34);
    const auto *qsd_35 = buffer.data(qsd + 35);
    const auto *qsd_36 = buffer.data(qsd + 36);
    const auto *qsd_37 = buffer.data(qsd + 37);
    const auto *qsd_39 = buffer.data(qsd + 39);
    const auto *qsd_41 = buffer.data(qsd + 41);
    const auto *qsd_42 = buffer.data(qsd + 42);
    const auto *qsd_45 = buffer.data(qsd + 45);
    const auto *qsd_46 = buffer.data(qsd + 46);
    const auto *qsd_47 = buffer.data(qsd + 47);
    const auto *qsd_48 = buffer.data(qsd + 48);
    const auto *qsd_51 = buffer.data(qsd + 51);
    const auto *qsd_52 = buffer.data(qsd + 52);
    const auto *qsd_53 = buffer.data(qsd + 53);
    const auto *qsd_54 = buffer.data(qsd + 54);
    const auto *qsd_56 = buffer.data(qsd + 56);
    const auto *qsd_57 = buffer.data(qsd + 57);
    const auto *qsd_58 = buffer.data(qsd + 58);
    const auto *qsd_59 = buffer.data(qsd + 59);
    const auto *qsd_60 = buffer.data(qsd + 60);
    const auto *qsd_61 = buffer.data(qsd + 61);
    const auto *qsd_63 = buffer.data(qsd + 63);
    const auto *qsd_65 = buffer.data(qsd + 65);
    const auto *qsd_66 = buffer.data(qsd + 66);
    const auto *qsd_69 = buffer.data(qsd + 69);
    const auto *qsd_70 = buffer.data(qsd + 70);
    const auto *qsd_71 = buffer.data(qsd + 71);
    const auto *qsd_72 = buffer.data(qsd + 72);
    const auto *qsd_75 = buffer.data(qsd + 75);
    const auto *qsd_76 = buffer.data(qsd + 76);
    const auto *qsd_77 = buffer.data(qsd + 77);
    const auto *qsd_78 = buffer.data(qsd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, osd_0, osd_3, qsp0_0, \
                         qsp1_0, qsd_0, qsd_2, qsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * osd_0[k]
                 + f_1 * qsp0_0[k]
                 - f_2 * qsp1_0[k]
                 + f_3 * pc_x[k] * qsd_0[k];

        t_1[k] = f_3 * pc_y[k] * qsd_0[k];

        t_2[k] = f_3 * pc_z[k] * qsd_0[k];

        t_3[k] = f_0 * osd_3[k]
                 + f_3 * pc_x[k] * qsd_3[k];

        t_4[k] = f_3 * pc_y[k] * qsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, osd_5, qsp0_1, qsp0_2, \
                         qsp1_1, qsp1_2, qsd_3, qsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * osd_5[k]
                 + f_3 * pc_x[k] * qsd_5[k];

        t_6[k] = f_1 * qsp0_1[k]
                 - f_2 * qsp1_1[k]
                 + f_3 * pc_y[k] * qsd_3[k];

        t_7[k] = f_3 * pc_z[k] * qsd_3[k];

        t_8[k] = f_3 * pc_y[k] * qsd_5[k];

        t_9[k] = f_1 * qsp0_2[k]
                 - f_2 * qsp1_2[k]
                 + f_3 * pc_z[k] * qsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, osf0_0, osd_0, \
                         osd_9, osf1_0, qsd_6, qsd_7, qsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * osf0_0[k]
                  - f_4 * pc_y[k] * osf1_0[k];

        t_11[k] = f_5 * osd_0[k]
                  + f_3 * pc_y[k] * qsd_6[k];

        t_12[k] = f_3 * pc_z[k] * qsd_6[k];

        t_13[k] = f_6 * osd_9[k]
                  + f_3 * pc_x[k] * qsd_9[k];

        t_14[k] = f_3 * pc_z[k] * qsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, osd_3, osd_5, osd_11, \
                         qsp0_4, qsp1_4, qsd_9, qsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * osd_11[k]
                  + f_3 * pc_x[k] * qsd_11[k];

        t_16[k] = f_5 * osd_3[k]
                  + f_1 * qsp0_4[k]
                  - f_2 * qsp1_4[k]
                  + f_3 * pc_y[k] * qsd_9[k];

        t_17[k] = f_3 * pc_z[k] * qsd_9[k];

        t_18[k] = f_5 * osd_5[k]
                  + f_3 * pc_y[k] * qsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, osf0_0, osf0_9, \
                         osd_0, osf1_0, osf1_9, qsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * osf0_9[k]
                  - f_4 * pc_y[k] * osf1_9[k];

        t_20[k] = pa_z[k] * osf0_0[k]
                  - f_4 * pc_z[k] * osf1_0[k];

        t_21[k] = f_3 * pc_y[k] * qsd_12[k];

        t_22[k] = f_5 * osd_0[k]
                  + f_3 * pc_z[k] * qsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, osf0_6, osd_15, \
                         osd_17, osf1_6, qsd_14, qsd_15, qsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * osd_15[k]
                  + f_3 * pc_x[k] * qsd_15[k];

        t_24[k] = f_3 * pc_y[k] * qsd_14[k];

        t_25[k] = f_6 * osd_17[k]
                  + f_3 * pc_x[k] * qsd_17[k];

        t_26[k] = pa_z[k] * osf0_6[k]
                  - f_4 * pc_z[k] * osf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, osd_5, osd_18, qsp0_8, \
                         qsp0_9, qsp1_8, qsp1_9, qsd_16, qsd_17, \
                         qsd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * qsp0_8[k]
                  - f_8 * qsp1_8[k]
                  + f_3 * pc_y[k] * qsd_16[k];

        t_28[k] = f_3 * pc_y[k] * qsd_17[k];

        t_29[k] = f_5 * osd_5[k]
                  + f_1 * qsp0_8[k]
                  - f_2 * qsp1_8[k]
                  + f_3 * pc_z[k] * qsd_17[k];

        t_30[k] = f_9 * osd_18[k]
                  + f_1 * qsp0_9[k]
                  - f_2 * qsp1_9[k]
                  + f_3 * pc_x[k] * qsd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, osd_6, osd_21, \
                         osd_23, qsd_18, qsd_19, qsd_21, qsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * osd_6[k]
                  + f_3 * pc_y[k] * qsd_18[k];

        t_32[k] = f_3 * pc_z[k] * qsd_18[k];

        t_33[k] = f_9 * osd_21[k]
                  + f_3 * pc_x[k] * qsd_21[k];

        t_34[k] = f_3 * pc_z[k] * qsd_19[k];

        t_35[k] = f_9 * osd_23[k]
                  + f_3 * pc_x[k] * qsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, osd_9, osd_11, qsp0_10, qsp0_11, \
                         qsp1_10, qsp1_11, qsd_21, qsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * osd_9[k]
                  + f_1 * qsp0_10[k]
                  - f_2 * qsp1_10[k]
                  + f_3 * pc_y[k] * qsd_21[k];

        t_37[k] = f_3 * pc_z[k] * qsd_21[k];

        t_38[k] = f_10 * osd_11[k]
                  + f_3 * pc_y[k] * qsd_23[k];

        t_39[k] = f_1 * qsp0_11[k]
                  - f_2 * qsp1_11[k]
                  + f_3 * pc_z[k] * qsd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, osf0_20, osd_6, \
                         osd_12, osd_27, osf1_20, qsd_24, qsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * osf0_20[k]
                  - f_4 * pc_y[k] * osf1_20[k];

        t_41[k] = f_5 * osd_12[k]
                  + f_3 * pc_y[k] * qsd_24[k];

        t_42[k] = f_5 * osd_6[k]
                  + f_3 * pc_z[k] * qsd_24[k];

        t_43[k] = f_9 * osd_27[k]
                  + f_3 * pc_x[k] * qsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, osf0_16, osd_9, osd_28, \
                         osd_29, osf1_16, qsd_27, qsd_28, qsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * osd_28[k]
                  + f_3 * pc_x[k] * qsd_28[k];

        t_45[k] = f_9 * osd_29[k]
                  + f_3 * pc_x[k] * qsd_29[k];

        t_46[k] = pa_z[k] * osf0_16[k]
                  - f_4 * pc_z[k] * osf1_16[k];

        t_47[k] = f_5 * osd_9[k]
                  + f_3 * pc_z[k] * qsd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, osf0_29, osd_17, osd_30, \
                         osf1_29, qsp0_15, qsp1_15, qsd_29, qsd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * osd_17[k]
                  + f_3 * pc_y[k] * qsd_29[k];

        t_49[k] = pa_y[k] * osf0_29[k]
                  - f_4 * pc_y[k] * osf1_29[k];

        t_50[k] = f_9 * osd_30[k]
                  + f_1 * qsp0_15[k]
                  - f_2 * qsp1_15[k]
                  + f_3 * pc_x[k] * qsd_30[k];

        t_51[k] = f_3 * pc_y[k] * qsd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, osd_12, osd_33, osd_35, \
                         qsd_30, qsd_32, qsd_33, qsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * osd_12[k]
                  + f_3 * pc_z[k] * qsd_30[k];

        t_53[k] = f_9 * osd_33[k]
                  + f_3 * pc_x[k] * qsd_33[k];

        t_54[k] = f_3 * pc_y[k] * qsd_32[k];

        t_55[k] = f_9 * osd_35[k]
                  + f_3 * pc_x[k] * qsd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, osd_17, qsp0_16, qsp0_17, \
                         qsp1_16, qsp1_17, qsd_33, qsd_34, qsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * qsp0_16[k]
                  - f_2 * qsp1_16[k]
                  + f_3 * pc_y[k] * qsd_33[k];

        t_57[k] = f_7 * qsp0_17[k]
                  - f_8 * qsp1_17[k]
                  + f_3 * pc_y[k] * qsd_34[k];

        t_58[k] = f_3 * pc_y[k] * qsd_35[k];

        t_59[k] = f_10 * osd_17[k]
                  + f_1 * qsp0_17[k]
                  - f_2 * qsp1_17[k]
                  + f_3 * pc_z[k] * qsd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, osd_18, osd_36, \
                         osd_39, qsp0_18, qsp1_18, qsd_36, qsd_37, \
                         qsd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * osd_36[k]
                  + f_1 * qsp0_18[k]
                  - f_2 * qsp1_18[k]
                  + f_3 * pc_x[k] * qsd_36[k];

        t_61[k] = f_12 * osd_18[k]
                  + f_3 * pc_y[k] * qsd_36[k];

        t_62[k] = f_3 * pc_z[k] * qsd_36[k];

        t_63[k] = f_11 * osd_39[k]
                  + f_3 * pc_x[k] * qsd_39[k];

        t_64[k] = f_3 * pc_z[k] * qsd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, osd_21, osd_23, osd_41, \
                         qsp0_19, qsp1_19, qsd_39, qsd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * osd_41[k]
                  + f_3 * pc_x[k] * qsd_41[k];

        t_66[k] = f_12 * osd_21[k]
                  + f_1 * qsp0_19[k]
                  - f_2 * qsp1_19[k]
                  + f_3 * pc_y[k] * qsd_39[k];

        t_67[k] = f_3 * pc_z[k] * qsd_39[k];

        t_68[k] = f_12 * osd_23[k]
                  + f_3 * pc_y[k] * qsd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, osf0_30, osd_18, osd_24, \
                         osf1_30, qsp0_20, qsp1_20, qsd_41, qsd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * qsp0_20[k]
                  - f_2 * qsp1_20[k]
                  + f_3 * pc_z[k] * qsd_41[k];

        t_70[k] = pa_z[k] * osf0_30[k]
                  - f_4 * pc_z[k] * osf1_30[k];

        t_71[k] = f_10 * osd_24[k]
                  + f_3 * pc_y[k] * qsd_42[k];

        t_72[k] = f_5 * osd_18[k]
                  + f_3 * pc_z[k] * qsd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, osf0_36, osd_45, osd_46, \
                         osd_47, osf1_36, qsd_45, qsd_46, qsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * osd_45[k]
                  + f_3 * pc_x[k] * qsd_45[k];

        t_74[k] = f_11 * osd_46[k]
                  + f_3 * pc_x[k] * qsd_46[k];

        t_75[k] = f_11 * osd_47[k]
                  + f_3 * pc_x[k] * qsd_47[k];

        t_76[k] = pa_z[k] * osf0_36[k]
                  - f_4 * pc_z[k] * osf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, osf0_50, osd_21, osd_23, \
                         osd_29, osf1_50, qsp0_23, qsp1_23, qsd_45, \
                         qsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * osd_21[k]
                  + f_3 * pc_z[k] * qsd_45[k];

        t_78[k] = f_10 * osd_29[k]
                  + f_3 * pc_y[k] * qsd_47[k];

        t_79[k] = f_5 * osd_23[k]
                  + f_1 * qsp0_23[k]
                  - f_2 * qsp1_23[k]
                  + f_3 * pc_z[k] * qsd_47[k];

        t_80[k] = pa_y[k] * osf0_50[k]
                  - f_4 * pc_y[k] * osf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, osd_24, osd_30, osd_51, \
                         osd_52, qsd_48, qsd_51, qsd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * osd_30[k]
                  + f_3 * pc_y[k] * qsd_48[k];

        t_82[k] = f_10 * osd_24[k]
                  + f_3 * pc_z[k] * qsd_48[k];

        t_83[k] = f_11 * osd_51[k]
                  + f_3 * pc_x[k] * qsd_51[k];

        t_84[k] = f_11 * osd_52[k]
                  + f_3 * pc_x[k] * qsd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, osd_27, osd_33, osd_35, \
                         osd_53, qsp0_25, qsp1_25, qsd_51, qsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * osd_53[k]
                  + f_3 * pc_x[k] * qsd_53[k];

        t_86[k] = f_5 * osd_33[k]
                  + f_1 * qsp0_25[k]
                  - f_2 * qsp1_25[k]
                  + f_3 * pc_y[k] * qsd_51[k];

        t_87[k] = f_10 * osd_27[k]
                  + f_3 * pc_z[k] * qsd_51[k];

        t_88[k] = f_5 * osd_35[k]
                  + f_3 * pc_y[k] * qsd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, osf0_59, osd_30, \
                         osd_54, osf1_59, qsp0_27, qsp1_27, qsd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * osf0_59[k]
                  - f_4 * pc_y[k] * osf1_59[k];

        t_90[k] = f_11 * osd_54[k]
                  + f_1 * qsp0_27[k]
                  - f_2 * qsp1_27[k]
                  + f_3 * pc_x[k] * qsd_54[k];

        t_91[k] = f_3 * pc_y[k] * qsd_54[k];

        t_92[k] = f_12 * osd_30[k]
                  + f_3 * pc_z[k] * qsd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, osd_57, osd_59, qsp0_28, qsp1_28, \
                         qsd_56, qsd_57, qsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * osd_57[k]
                  + f_3 * pc_x[k] * qsd_57[k];

        t_94[k] = f_3 * pc_y[k] * qsd_56[k];

        t_95[k] = f_11 * osd_59[k]
                  + f_3 * pc_x[k] * qsd_59[k];

        t_96[k] = f_1 * qsp0_28[k]
                  - f_2 * qsp1_28[k]
                  + f_3 * pc_y[k] * qsd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, osd_35, osd_60, qsp0_29, \
                         qsp0_30, qsp1_29, qsp1_30, qsd_58, qsd_59, \
                         qsd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * qsp0_29[k]
                  - f_8 * qsp1_29[k]
                  + f_3 * pc_y[k] * qsd_58[k];

        t_98[k] = f_3 * pc_y[k] * qsd_59[k];

        t_99[k] = f_12 * osd_35[k]
                  + f_1 * qsp0_29[k]
                  - f_2 * qsp1_29[k]
                  + f_3 * pc_z[k] * qsd_59[k];

        t_100[k] = f_13 * osd_60[k]
                   + f_1 * qsp0_30[k]
                   - f_2 * qsp1_30[k]
                   + f_3 * pc_x[k] * qsd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, osd_36, osd_63, \
                         osd_65, qsd_60, qsd_61, qsd_63, qsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_14 * osd_36[k]
                   + f_3 * pc_y[k] * qsd_60[k];

        t_102[k] = f_3 * pc_z[k] * qsd_60[k];

        t_103[k] = f_13 * osd_63[k]
                   + f_3 * pc_x[k] * qsd_63[k];

        t_104[k] = f_3 * pc_z[k] * qsd_61[k];

        t_105[k] = f_13 * osd_65[k]
                   + f_3 * pc_x[k] * qsd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, osd_39, osd_41, qsp0_31, \
                         qsp0_32, qsp1_31, qsp1_32, qsd_63, qsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_14 * osd_39[k]
                   + f_1 * qsp0_31[k]
                   - f_2 * qsp1_31[k]
                   + f_3 * pc_y[k] * qsd_63[k];

        t_107[k] = f_3 * pc_z[k] * qsd_63[k];

        t_108[k] = f_14 * osd_41[k]
                   + f_3 * pc_y[k] * qsd_65[k];

        t_109[k] = f_1 * qsp0_32[k]
                   - f_2 * qsp1_32[k]
                   + f_3 * pc_z[k] * qsd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, osf0_60, osd_36, \
                         osd_42, osd_69, osf1_60, qsd_66, qsd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * osf0_60[k]
                   - f_4 * pc_z[k] * osf1_60[k];

        t_111[k] = f_12 * osd_42[k]
                   + f_3 * pc_y[k] * qsd_66[k];

        t_112[k] = f_5 * osd_36[k]
                   + f_3 * pc_z[k] * qsd_66[k];

        t_113[k] = f_13 * osd_69[k]
                   + f_3 * pc_x[k] * qsd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, osf0_66, osd_39, \
                         osd_70, osd_71, osf1_66, qsd_69, qsd_70, \
                         qsd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * osd_70[k]
                   + f_3 * pc_x[k] * qsd_70[k];

        t_115[k] = f_13 * osd_71[k]
                   + f_3 * pc_x[k] * qsd_71[k];

        t_116[k] = pa_z[k] * osf0_66[k]
                   - f_4 * pc_z[k] * osf1_66[k];

        t_117[k] = f_5 * osd_39[k]
                   + f_3 * pc_z[k] * qsd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, osd_41, osd_47, osd_72, \
                         qsp0_35, qsp0_36, qsp1_35, qsp1_36, qsd_71, \
                         qsd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * osd_47[k]
                   + f_3 * pc_y[k] * qsd_71[k];

        t_119[k] = f_5 * osd_41[k]
                   + f_1 * qsp0_35[k]
                   - f_2 * qsp1_35[k]
                   + f_3 * pc_z[k] * qsd_71[k];

        t_120[k] = f_13 * osd_72[k]
                   + f_1 * qsp0_36[k]
                   - f_2 * qsp1_36[k]
                   + f_3 * pc_x[k] * qsd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, osd_42, osd_48, osd_75, \
                         osd_76, qsd_72, qsd_75, qsd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * osd_48[k]
                   + f_3 * pc_y[k] * qsd_72[k];

        t_122[k] = f_10 * osd_42[k]
                   + f_3 * pc_z[k] * qsd_72[k];

        t_123[k] = f_13 * osd_75[k]
                   + f_3 * pc_x[k] * qsd_75[k];

        t_124[k] = f_13 * osd_76[k]
                   + f_3 * pc_x[k] * qsd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, osd_45, osd_51, osd_53, \
                         osd_77, qsp0_37, qsp1_37, qsd_75, qsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * osd_77[k]
                   + f_3 * pc_x[k] * qsd_77[k];

        t_126[k] = f_10 * osd_51[k]
                   + f_1 * qsp0_37[k]
                   - f_2 * qsp1_37[k]
                   + f_3 * pc_y[k] * qsd_75[k];

        t_127[k] = f_10 * osd_45[k]
                   + f_3 * pc_z[k] * qsd_75[k];

        t_128[k] = f_10 * osd_53[k]
                   + f_3 * pc_y[k] * qsd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, osf0_90, osd_47, \
                         osd_48, osd_54, osf1_90, qsp0_38, qsp1_38, qsd_77, \
                         qsd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * osd_47[k]
                   + f_1 * qsp0_38[k]
                   - f_2 * qsp1_38[k]
                   + f_3 * pc_z[k] * qsd_77[k];

        t_130[k] = pa_y[k] * osf0_90[k]
                   - f_4 * pc_y[k] * osf1_90[k];

        t_131[k] = f_5 * osd_54[k]
                   + f_3 * pc_y[k] * qsd_78[k];

        t_132[k] = f_12 * osd_48[k]
                   + f_3 * pc_z[k] * qsd_78[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_99 = buffer.data(osf0 + 99);
    const auto *osf0_100 = buffer.data(osf0 + 100);
    const auto *osf0_106 = buffer.data(osf0 + 106);
    const auto *osf0_140 = buffer.data(osf0 + 140);
    const auto *osf0_149 = buffer.data(osf0 + 149);
    const auto *osf0_150 = buffer.data(osf0 + 150);
    const auto *osf0_156 = buffer.data(osf0 + 156);

    const auto *osd_51 = buffer.data(osd + 51);
    const auto *osd_54 = buffer.data(osd + 54);
    const auto *osd_57 = buffer.data(osd + 57);
    const auto *osd_59 = buffer.data(osd + 59);
    const auto *osd_60 = buffer.data(osd + 60);
    const auto *osd_63 = buffer.data(osd + 63);
    const auto *osd_65 = buffer.data(osd + 65);
    const auto *osd_66 = buffer.data(osd + 66);
    const auto *osd_69 = buffer.data(osd + 69);
    const auto *osd_71 = buffer.data(osd + 71);
    const auto *osd_72 = buffer.data(osd + 72);
    const auto *osd_75 = buffer.data(osd + 75);
    const auto *osd_77 = buffer.data(osd + 77);
    const auto *osd_78 = buffer.data(osd + 78);
    const auto *osd_81 = buffer.data(osd + 81);
    const auto *osd_82 = buffer.data(osd + 82);
    const auto *osd_83 = buffer.data(osd + 83);
    const auto *osd_84 = buffer.data(osd + 84);
    const auto *osd_87 = buffer.data(osd + 87);
    const auto *osd_89 = buffer.data(osd + 89);
    const auto *osd_90 = buffer.data(osd + 90);
    const auto *osd_93 = buffer.data(osd + 93);
    const auto *osd_95 = buffer.data(osd + 95);
    const auto *osd_96 = buffer.data(osd + 96);
    const auto *osd_99 = buffer.data(osd + 99);
    const auto *osd_100 = buffer.data(osd + 100);
    const auto *osd_101 = buffer.data(osd + 101);
    const auto *osd_102 = buffer.data(osd + 102);
    const auto *osd_105 = buffer.data(osd + 105);
    const auto *osd_106 = buffer.data(osd + 106);
    const auto *osd_107 = buffer.data(osd + 107);
    const auto *osd_108 = buffer.data(osd + 108);
    const auto *osd_111 = buffer.data(osd + 111);
    const auto *osd_112 = buffer.data(osd + 112);
    const auto *osd_113 = buffer.data(osd + 113);
    const auto *osd_114 = buffer.data(osd + 114);
    const auto *osd_117 = buffer.data(osd + 117);
    const auto *osd_118 = buffer.data(osd + 118);
    const auto *osd_119 = buffer.data(osd + 119);
    const auto *osd_120 = buffer.data(osd + 120);
    const auto *osd_123 = buffer.data(osd + 123);
    const auto *osd_125 = buffer.data(osd + 125);
    const auto *osd_126 = buffer.data(osd + 126);
    const auto *osd_129 = buffer.data(osd + 129);
    const auto *osd_131 = buffer.data(osd + 131);
    const auto *osd_135 = buffer.data(osd + 135);
    const auto *osd_136 = buffer.data(osd + 136);
    const auto *osd_137 = buffer.data(osd + 137);
    const auto *osd_138 = buffer.data(osd + 138);
    const auto *osd_141 = buffer.data(osd + 141);
    const auto *osd_142 = buffer.data(osd + 142);
    const auto *osd_143 = buffer.data(osd + 143);
    const auto *osd_144 = buffer.data(osd + 144);
    const auto *osd_147 = buffer.data(osd + 147);
    const auto *osd_148 = buffer.data(osd + 148);
    const auto *osd_149 = buffer.data(osd + 149);
    const auto *osd_150 = buffer.data(osd + 150);
    const auto *osd_153 = buffer.data(osd + 153);
    const auto *osd_154 = buffer.data(osd + 154);
    const auto *osd_155 = buffer.data(osd + 155);

    const auto *osf1_99 = buffer.data(osf1 + 99);
    const auto *osf1_100 = buffer.data(osf1 + 100);
    const auto *osf1_106 = buffer.data(osf1 + 106);
    const auto *osf1_140 = buffer.data(osf1 + 140);
    const auto *osf1_149 = buffer.data(osf1 + 149);
    const auto *osf1_150 = buffer.data(osf1 + 150);
    const auto *osf1_156 = buffer.data(osf1 + 156);

    const auto *qsp0_40 = buffer.data(qsp0 + 40);
    const auto *qsp0_42 = buffer.data(qsp0 + 42);
    const auto *qsp0_43 = buffer.data(qsp0 + 43);
    const auto *qsp0_44 = buffer.data(qsp0 + 44);
    const auto *qsp0_45 = buffer.data(qsp0 + 45);
    const auto *qsp0_46 = buffer.data(qsp0 + 46);
    const auto *qsp0_47 = buffer.data(qsp0 + 47);
    const auto *qsp0_50 = buffer.data(qsp0 + 50);
    const auto *qsp0_51 = buffer.data(qsp0 + 51);
    const auto *qsp0_52 = buffer.data(qsp0 + 52);
    const auto *qsp0_53 = buffer.data(qsp0 + 53);
    const auto *qsp0_54 = buffer.data(qsp0 + 54);
    const auto *qsp0_55 = buffer.data(qsp0 + 55);
    const auto *qsp0_56 = buffer.data(qsp0 + 56);
    const auto *qsp0_58 = buffer.data(qsp0 + 58);
    const auto *qsp0_60 = buffer.data(qsp0 + 60);
    const auto *qsp0_61 = buffer.data(qsp0 + 61);
    const auto *qsp0_62 = buffer.data(qsp0 + 62);
    const auto *qsp0_63 = buffer.data(qsp0 + 63);
    const auto *qsp0_64 = buffer.data(qsp0 + 64);
    const auto *qsp0_65 = buffer.data(qsp0 + 65);
    const auto *qsp0_68 = buffer.data(qsp0 + 68);
    const auto *qsp0_69 = buffer.data(qsp0 + 69);
    const auto *qsp0_70 = buffer.data(qsp0 + 70);
    const auto *qsp0_71 = buffer.data(qsp0 + 71);
    const auto *qsp0_72 = buffer.data(qsp0 + 72);
    const auto *qsp0_73 = buffer.data(qsp0 + 73);
    const auto *qsp0_74 = buffer.data(qsp0 + 74);
    const auto *qsp0_75 = buffer.data(qsp0 + 75);

    const auto *qsp1_40 = buffer.data(qsp1 + 40);
    const auto *qsp1_42 = buffer.data(qsp1 + 42);
    const auto *qsp1_43 = buffer.data(qsp1 + 43);
    const auto *qsp1_44 = buffer.data(qsp1 + 44);
    const auto *qsp1_45 = buffer.data(qsp1 + 45);
    const auto *qsp1_46 = buffer.data(qsp1 + 46);
    const auto *qsp1_47 = buffer.data(qsp1 + 47);
    const auto *qsp1_50 = buffer.data(qsp1 + 50);
    const auto *qsp1_51 = buffer.data(qsp1 + 51);
    const auto *qsp1_52 = buffer.data(qsp1 + 52);
    const auto *qsp1_53 = buffer.data(qsp1 + 53);
    const auto *qsp1_54 = buffer.data(qsp1 + 54);
    const auto *qsp1_55 = buffer.data(qsp1 + 55);
    const auto *qsp1_56 = buffer.data(qsp1 + 56);
    const auto *qsp1_58 = buffer.data(qsp1 + 58);
    const auto *qsp1_60 = buffer.data(qsp1 + 60);
    const auto *qsp1_61 = buffer.data(qsp1 + 61);
    const auto *qsp1_62 = buffer.data(qsp1 + 62);
    const auto *qsp1_63 = buffer.data(qsp1 + 63);
    const auto *qsp1_64 = buffer.data(qsp1 + 64);
    const auto *qsp1_65 = buffer.data(qsp1 + 65);
    const auto *qsp1_68 = buffer.data(qsp1 + 68);
    const auto *qsp1_69 = buffer.data(qsp1 + 69);
    const auto *qsp1_70 = buffer.data(qsp1 + 70);
    const auto *qsp1_71 = buffer.data(qsp1 + 71);
    const auto *qsp1_72 = buffer.data(qsp1 + 72);
    const auto *qsp1_73 = buffer.data(qsp1 + 73);
    const auto *qsp1_74 = buffer.data(qsp1 + 74);
    const auto *qsp1_75 = buffer.data(qsp1 + 75);

    const auto *qsd_81 = buffer.data(qsd + 81);
    const auto *qsd_82 = buffer.data(qsd + 82);
    const auto *qsd_83 = buffer.data(qsd + 83);
    const auto *qsd_84 = buffer.data(qsd + 84);
    const auto *qsd_86 = buffer.data(qsd + 86);
    const auto *qsd_87 = buffer.data(qsd + 87);
    const auto *qsd_88 = buffer.data(qsd + 88);
    const auto *qsd_89 = buffer.data(qsd + 89);
    const auto *qsd_90 = buffer.data(qsd + 90);
    const auto *qsd_91 = buffer.data(qsd + 91);
    const auto *qsd_93 = buffer.data(qsd + 93);
    const auto *qsd_95 = buffer.data(qsd + 95);
    const auto *qsd_96 = buffer.data(qsd + 96);
    const auto *qsd_99 = buffer.data(qsd + 99);
    const auto *qsd_100 = buffer.data(qsd + 100);
    const auto *qsd_101 = buffer.data(qsd + 101);
    const auto *qsd_102 = buffer.data(qsd + 102);
    const auto *qsd_105 = buffer.data(qsd + 105);
    const auto *qsd_106 = buffer.data(qsd + 106);
    const auto *qsd_107 = buffer.data(qsd + 107);
    const auto *qsd_108 = buffer.data(qsd + 108);
    const auto *qsd_111 = buffer.data(qsd + 111);
    const auto *qsd_112 = buffer.data(qsd + 112);
    const auto *qsd_113 = buffer.data(qsd + 113);
    const auto *qsd_114 = buffer.data(qsd + 114);
    const auto *qsd_117 = buffer.data(qsd + 117);
    const auto *qsd_118 = buffer.data(qsd + 118);
    const auto *qsd_119 = buffer.data(qsd + 119);
    const auto *qsd_120 = buffer.data(qsd + 120);
    const auto *qsd_122 = buffer.data(qsd + 122);
    const auto *qsd_123 = buffer.data(qsd + 123);
    const auto *qsd_124 = buffer.data(qsd + 124);
    const auto *qsd_125 = buffer.data(qsd + 125);
    const auto *qsd_126 = buffer.data(qsd + 126);
    const auto *qsd_127 = buffer.data(qsd + 127);
    const auto *qsd_129 = buffer.data(qsd + 129);
    const auto *qsd_131 = buffer.data(qsd + 131);
    const auto *qsd_132 = buffer.data(qsd + 132);
    const auto *qsd_135 = buffer.data(qsd + 135);
    const auto *qsd_136 = buffer.data(qsd + 136);
    const auto *qsd_137 = buffer.data(qsd + 137);
    const auto *qsd_138 = buffer.data(qsd + 138);
    const auto *qsd_141 = buffer.data(qsd + 141);
    const auto *qsd_142 = buffer.data(qsd + 142);
    const auto *qsd_143 = buffer.data(qsd + 143);
    const auto *qsd_144 = buffer.data(qsd + 144);
    const auto *qsd_147 = buffer.data(qsd + 147);
    const auto *qsd_148 = buffer.data(qsd + 148);
    const auto *qsd_149 = buffer.data(qsd + 149);
    const auto *qsd_150 = buffer.data(qsd + 150);
    const auto *qsd_153 = buffer.data(qsd + 153);
    const auto *qsd_154 = buffer.data(qsd + 154);
    const auto *qsd_155 = buffer.data(qsd + 155);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, osd_57, osd_81, osd_82, \
                         osd_83, qsp0_40, qsp1_40, qsd_81, qsd_82, \
                         qsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_13 * osd_81[k]
                   + f_3 * pc_x[k] * qsd_81[k];

        t_134[k] = f_13 * osd_82[k]
                   + f_3 * pc_x[k] * qsd_82[k];

        t_135[k] = f_13 * osd_83[k]
                   + f_3 * pc_x[k] * qsd_83[k];

        t_136[k] = f_5 * osd_57[k]
                   + f_1 * qsp0_40[k]
                   - f_2 * qsp1_40[k]
                   + f_3 * pc_y[k] * qsd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, osf0_99, osd_51, osd_59, \
                         osf1_99, qsd_81, qsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * osd_51[k]
                   + f_3 * pc_z[k] * qsd_81[k];

        t_138[k] = f_5 * osd_59[k]
                   + f_3 * pc_y[k] * qsd_83[k];

        t_139[k] = pa_y[k] * osf0_99[k]
                   - f_4 * pc_y[k] * osf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, osd_54, osd_84, \
                         osd_87, qsp0_42, qsp1_42, qsd_84, qsd_86, \
                         qsd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * osd_84[k]
                   + f_1 * qsp0_42[k]
                   - f_2 * qsp1_42[k]
                   + f_3 * pc_x[k] * qsd_84[k];

        t_141[k] = f_3 * pc_y[k] * qsd_84[k];

        t_142[k] = f_14 * osd_54[k]
                   + f_3 * pc_z[k] * qsd_84[k];

        t_143[k] = f_13 * osd_87[k]
                   + f_3 * pc_x[k] * qsd_87[k];

        t_144[k] = f_3 * pc_y[k] * qsd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, osd_89, qsp0_43, qsp0_44, \
                         qsp1_43, qsp1_44, qsd_87, qsd_88, qsd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * osd_89[k]
                   + f_3 * pc_x[k] * qsd_89[k];

        t_146[k] = f_1 * qsp0_43[k]
                   - f_2 * qsp1_43[k]
                   + f_3 * pc_y[k] * qsd_87[k];

        t_147[k] = f_7 * qsp0_44[k]
                   - f_8 * qsp1_44[k]
                   + f_3 * pc_y[k] * qsd_88[k];

        t_148[k] = f_3 * pc_y[k] * qsd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_y, pc_z, osd_59, osd_60, osd_90, \
                         qsp0_44, qsp0_45, qsp1_44, qsp1_45, qsd_89, \
                         qsd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_14 * osd_59[k]
                   + f_1 * qsp0_44[k]
                   - f_2 * qsp1_44[k]
                   + f_3 * pc_z[k] * qsd_89[k];

        t_150[k] = f_15 * osd_90[k]
                   + f_1 * qsp0_45[k]
                   - f_2 * qsp1_45[k]
                   + f_3 * pc_x[k] * qsd_90[k];

        t_151[k] = f_16 * osd_60[k]
                   + f_3 * pc_y[k] * qsd_90[k];

        t_152[k] = f_3 * pc_z[k] * qsd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, osd_63, osd_93, \
                         osd_95, qsp0_46, qsp1_46, qsd_91, qsd_93, \
                         qsd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_15 * osd_93[k]
                   + f_3 * pc_x[k] * qsd_93[k];

        t_154[k] = f_3 * pc_z[k] * qsd_91[k];

        t_155[k] = f_15 * osd_95[k]
                   + f_3 * pc_x[k] * qsd_95[k];

        t_156[k] = f_16 * osd_63[k]
                   + f_1 * qsp0_46[k]
                   - f_2 * qsp1_46[k]
                   + f_3 * pc_y[k] * qsd_93[k];

        t_157[k] = f_3 * pc_z[k] * qsd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pc_y, pc_z, osf0_100, osd_65, \
                         osd_66, osf1_100, qsp0_47, qsp1_47, qsd_95, \
                         qsd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_16 * osd_65[k]
                   + f_3 * pc_y[k] * qsd_95[k];

        t_159[k] = f_1 * qsp0_47[k]
                   - f_2 * qsp1_47[k]
                   + f_3 * pc_z[k] * qsd_95[k];

        t_160[k] = pa_z[k] * osf0_100[k]
                   - f_4 * pc_z[k] * osf1_100[k];

        t_161[k] = f_14 * osd_66[k]
                   + f_3 * pc_y[k] * qsd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, osd_60, osd_99, osd_100, \
                         osd_101, qsd_96, qsd_99, qsd_100, qsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * osd_60[k]
                   + f_3 * pc_z[k] * qsd_96[k];

        t_163[k] = f_15 * osd_99[k]
                   + f_3 * pc_x[k] * qsd_99[k];

        t_164[k] = f_15 * osd_100[k]
                   + f_3 * pc_x[k] * qsd_100[k];

        t_165[k] = f_15 * osd_101[k]
                   + f_3 * pc_x[k] * qsd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, osf0_106, osd_63, \
                         osd_65, osd_71, osf1_106, qsp0_50, qsp1_50, qsd_99, \
                         qsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * osf0_106[k]
                   - f_4 * pc_z[k] * osf1_106[k];

        t_167[k] = f_5 * osd_63[k]
                   + f_3 * pc_z[k] * qsd_99[k];

        t_168[k] = f_14 * osd_71[k]
                   + f_3 * pc_y[k] * qsd_101[k];

        t_169[k] = f_5 * osd_65[k]
                   + f_1 * qsp0_50[k]
                   - f_2 * qsp1_50[k]
                   + f_3 * pc_z[k] * qsd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, osd_66, osd_72, \
                         osd_102, osd_105, qsp0_51, qsp1_51, qsd_102, \
                         qsd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_15 * osd_102[k]
                   + f_1 * qsp0_51[k]
                   - f_2 * qsp1_51[k]
                   + f_3 * pc_x[k] * qsd_102[k];

        t_171[k] = f_12 * osd_72[k]
                   + f_3 * pc_y[k] * qsd_102[k];

        t_172[k] = f_10 * osd_66[k]
                   + f_3 * pc_z[k] * qsd_102[k];

        t_173[k] = f_15 * osd_105[k]
                   + f_3 * pc_x[k] * qsd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, osd_69, osd_75, \
                         osd_106, osd_107, qsp0_52, qsp1_52, qsd_105, qsd_106, \
                         qsd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_15 * osd_106[k]
                   + f_3 * pc_x[k] * qsd_106[k];

        t_175[k] = f_15 * osd_107[k]
                   + f_3 * pc_x[k] * qsd_107[k];

        t_176[k] = f_12 * osd_75[k]
                   + f_1 * qsp0_52[k]
                   - f_2 * qsp1_52[k]
                   + f_3 * pc_y[k] * qsd_105[k];

        t_177[k] = f_10 * osd_69[k]
                   + f_3 * pc_z[k] * qsd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, osd_71, osd_77, osd_108, \
                         qsp0_53, qsp0_54, qsp1_53, qsp1_54, qsd_107, \
                         qsd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * osd_77[k]
                   + f_3 * pc_y[k] * qsd_107[k];

        t_179[k] = f_10 * osd_71[k]
                   + f_1 * qsp0_53[k]
                   - f_2 * qsp1_53[k]
                   + f_3 * pc_z[k] * qsd_107[k];

        t_180[k] = f_15 * osd_108[k]
                   + f_1 * qsp0_54[k]
                   - f_2 * qsp1_54[k]
                   + f_3 * pc_x[k] * qsd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, osd_72, osd_78, \
                         osd_111, osd_112, qsd_108, qsd_111, qsd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * osd_78[k]
                   + f_3 * pc_y[k] * qsd_108[k];

        t_182[k] = f_12 * osd_72[k]
                   + f_3 * pc_z[k] * qsd_108[k];

        t_183[k] = f_15 * osd_111[k]
                   + f_3 * pc_x[k] * qsd_111[k];

        t_184[k] = f_15 * osd_112[k]
                   + f_3 * pc_x[k] * qsd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, osd_75, osd_81, osd_83, \
                         osd_113, qsp0_55, qsp1_55, qsd_111, qsd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * osd_113[k]
                   + f_3 * pc_x[k] * qsd_113[k];

        t_186[k] = f_10 * osd_81[k]
                   + f_1 * qsp0_55[k]
                   - f_2 * qsp1_55[k]
                   + f_3 * pc_y[k] * qsd_111[k];

        t_187[k] = f_12 * osd_75[k]
                   + f_3 * pc_z[k] * qsd_111[k];

        t_188[k] = f_10 * osd_83[k]
                   + f_3 * pc_y[k] * qsd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_y, pc_y, pc_z, osf0_140, osd_77, \
                         osd_78, osd_84, osf1_140, qsp0_56, qsp1_56, qsd_113, \
                         qsd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * osd_77[k]
                   + f_1 * qsp0_56[k]
                   - f_2 * qsp1_56[k]
                   + f_3 * pc_z[k] * qsd_113[k];

        t_190[k] = pa_y[k] * osf0_140[k]
                   - f_4 * pc_y[k] * osf1_140[k];

        t_191[k] = f_5 * osd_84[k]
                   + f_3 * pc_y[k] * qsd_114[k];

        t_192[k] = f_14 * osd_78[k]
                   + f_3 * pc_z[k] * qsd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, osd_87, osd_117, osd_118, \
                         osd_119, qsp0_58, qsp1_58, qsd_117, qsd_118, \
                         qsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_15 * osd_117[k]
                   + f_3 * pc_x[k] * qsd_117[k];

        t_194[k] = f_15 * osd_118[k]
                   + f_3 * pc_x[k] * qsd_118[k];

        t_195[k] = f_15 * osd_119[k]
                   + f_3 * pc_x[k] * qsd_119[k];

        t_196[k] = f_5 * osd_87[k]
                   + f_1 * qsp0_58[k]
                   - f_2 * qsp1_58[k]
                   + f_3 * pc_y[k] * qsd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, osf0_149, osd_81, osd_89, \
                         osf1_149, qsd_117, qsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * osd_81[k]
                   + f_3 * pc_z[k] * qsd_117[k];

        t_198[k] = f_5 * osd_89[k]
                   + f_3 * pc_y[k] * qsd_119[k];

        t_199[k] = pa_y[k] * osf0_149[k]
                   - f_4 * pc_y[k] * osf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, osd_84, osd_120, \
                         osd_123, qsp0_60, qsp1_60, qsd_120, qsd_122, \
                         qsd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_15 * osd_120[k]
                   + f_1 * qsp0_60[k]
                   - f_2 * qsp1_60[k]
                   + f_3 * pc_x[k] * qsd_120[k];

        t_201[k] = f_3 * pc_y[k] * qsd_120[k];

        t_202[k] = f_16 * osd_84[k]
                   + f_3 * pc_z[k] * qsd_120[k];

        t_203[k] = f_15 * osd_123[k]
                   + f_3 * pc_x[k] * qsd_123[k];

        t_204[k] = f_3 * pc_y[k] * qsd_122[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, osd_125, qsp0_61, qsp0_62, \
                         qsp1_61, qsp1_62, qsd_123, qsd_124, qsd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_15 * osd_125[k]
                   + f_3 * pc_x[k] * qsd_125[k];

        t_206[k] = f_1 * qsp0_61[k]
                   - f_2 * qsp1_61[k]
                   + f_3 * pc_y[k] * qsd_123[k];

        t_207[k] = f_7 * qsp0_62[k]
                   - f_8 * qsp1_62[k]
                   + f_3 * pc_y[k] * qsd_124[k];

        t_208[k] = f_3 * pc_y[k] * qsd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, osd_89, osd_90, \
                         osd_126, qsp0_62, qsp0_63, qsp1_62, qsp1_63, qsd_125, \
                         qsd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_16 * osd_89[k]
                   + f_1 * qsp0_62[k]
                   - f_2 * qsp1_62[k]
                   + f_3 * pc_z[k] * qsd_125[k];

        t_210[k] = f_17 * osd_126[k]
                   + f_1 * qsp0_63[k]
                   - f_2 * qsp1_63[k]
                   + f_3 * pc_x[k] * qsd_126[k];

        t_211[k] = f_17 * osd_90[k]
                   + f_3 * pc_y[k] * qsd_126[k];

        t_212[k] = f_3 * pc_z[k] * qsd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pc_x, pc_y, pc_z, osd_93, osd_129, \
                         osd_131, qsp0_64, qsp1_64, qsd_127, qsd_129, \
                         qsd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_17 * osd_129[k]
                   + f_3 * pc_x[k] * qsd_129[k];

        t_214[k] = f_3 * pc_z[k] * qsd_127[k];

        t_215[k] = f_17 * osd_131[k]
                   + f_3 * pc_x[k] * qsd_131[k];

        t_216[k] = f_17 * osd_93[k]
                   + f_1 * qsp0_64[k]
                   - f_2 * qsp1_64[k]
                   + f_3 * pc_y[k] * qsd_129[k];

        t_217[k] = f_3 * pc_z[k] * qsd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pc_y, pc_z, osf0_150, osd_95, \
                         osd_96, osf1_150, qsp0_65, qsp1_65, qsd_131, \
                         qsd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_17 * osd_95[k]
                   + f_3 * pc_y[k] * qsd_131[k];

        t_219[k] = f_1 * qsp0_65[k]
                   - f_2 * qsp1_65[k]
                   + f_3 * pc_z[k] * qsd_131[k];

        t_220[k] = pa_z[k] * osf0_150[k]
                   - f_4 * pc_z[k] * osf1_150[k];

        t_221[k] = f_16 * osd_96[k]
                   + f_3 * pc_y[k] * qsd_132[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_z, osd_90, osd_135, osd_136, \
                         osd_137, qsd_132, qsd_135, qsd_136, qsd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_5 * osd_90[k]
                   + f_3 * pc_z[k] * qsd_132[k];

        t_223[k] = f_17 * osd_135[k]
                   + f_3 * pc_x[k] * qsd_135[k];

        t_224[k] = f_17 * osd_136[k]
                   + f_3 * pc_x[k] * qsd_136[k];

        t_225[k] = f_17 * osd_137[k]
                   + f_3 * pc_x[k] * qsd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_z, pc_y, pc_z, osf0_156, osd_93, \
                         osd_95, osd_101, osf1_156, qsp0_68, qsp1_68, qsd_135, \
                         qsd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * osf0_156[k]
                   - f_4 * pc_z[k] * osf1_156[k];

        t_227[k] = f_5 * osd_93[k]
                   + f_3 * pc_z[k] * qsd_135[k];

        t_228[k] = f_16 * osd_101[k]
                   + f_3 * pc_y[k] * qsd_137[k];

        t_229[k] = f_5 * osd_95[k]
                   + f_1 * qsp0_68[k]
                   - f_2 * qsp1_68[k]
                   + f_3 * pc_z[k] * qsd_137[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, osd_96, osd_102, \
                         osd_138, osd_141, qsp0_69, qsp1_69, qsd_138, \
                         qsd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_17 * osd_138[k]
                   + f_1 * qsp0_69[k]
                   - f_2 * qsp1_69[k]
                   + f_3 * pc_x[k] * qsd_138[k];

        t_231[k] = f_14 * osd_102[k]
                   + f_3 * pc_y[k] * qsd_138[k];

        t_232[k] = f_10 * osd_96[k]
                   + f_3 * pc_z[k] * qsd_138[k];

        t_233[k] = f_17 * osd_141[k]
                   + f_3 * pc_x[k] * qsd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, osd_99, osd_105, \
                         osd_142, osd_143, qsp0_70, qsp1_70, qsd_141, qsd_142, \
                         qsd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_17 * osd_142[k]
                   + f_3 * pc_x[k] * qsd_142[k];

        t_235[k] = f_17 * osd_143[k]
                   + f_3 * pc_x[k] * qsd_143[k];

        t_236[k] = f_14 * osd_105[k]
                   + f_1 * qsp0_70[k]
                   - f_2 * qsp1_70[k]
                   + f_3 * pc_y[k] * qsd_141[k];

        t_237[k] = f_10 * osd_99[k]
                   + f_3 * pc_z[k] * qsd_141[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pc_x, pc_y, pc_z, osd_101, osd_107, osd_144, \
                         qsp0_71, qsp0_72, qsp1_71, qsp1_72, qsd_143, \
                         qsd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_14 * osd_107[k]
                   + f_3 * pc_y[k] * qsd_143[k];

        t_239[k] = f_10 * osd_101[k]
                   + f_1 * qsp0_71[k]
                   - f_2 * qsp1_71[k]
                   + f_3 * pc_z[k] * qsd_143[k];

        t_240[k] = f_17 * osd_144[k]
                   + f_1 * qsp0_72[k]
                   - f_2 * qsp1_72[k]
                   + f_3 * pc_x[k] * qsd_144[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, osd_102, osd_108, \
                         osd_147, osd_148, qsd_144, qsd_147, qsd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_12 * osd_108[k]
                   + f_3 * pc_y[k] * qsd_144[k];

        t_242[k] = f_12 * osd_102[k]
                   + f_3 * pc_z[k] * qsd_144[k];

        t_243[k] = f_17 * osd_147[k]
                   + f_3 * pc_x[k] * qsd_147[k];

        t_244[k] = f_17 * osd_148[k]
                   + f_3 * pc_x[k] * qsd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_y, pc_z, osd_105, osd_111, \
                         osd_113, osd_149, qsp0_73, qsp1_73, qsd_147, \
                         qsd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_17 * osd_149[k]
                   + f_3 * pc_x[k] * qsd_149[k];

        t_246[k] = f_12 * osd_111[k]
                   + f_1 * qsp0_73[k]
                   - f_2 * qsp1_73[k]
                   + f_3 * pc_y[k] * qsd_147[k];

        t_247[k] = f_12 * osd_105[k]
                   + f_3 * pc_z[k] * qsd_147[k];

        t_248[k] = f_12 * osd_113[k]
                   + f_3 * pc_y[k] * qsd_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_x, pc_y, pc_z, osd_107, osd_114, osd_150, \
                         qsp0_74, qsp0_75, qsp1_74, qsp1_75, qsd_149, \
                         qsd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * osd_107[k]
                   + f_1 * qsp0_74[k]
                   - f_2 * qsp1_74[k]
                   + f_3 * pc_z[k] * qsd_149[k];

        t_250[k] = f_17 * osd_150[k]
                   + f_1 * qsp0_75[k]
                   - f_2 * qsp1_75[k]
                   + f_3 * pc_x[k] * qsd_150[k];

        t_251[k] = f_10 * osd_114[k]
                   + f_3 * pc_y[k] * qsd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, osd_108, osd_153, osd_154, \
                         osd_155, qsd_150, qsd_153, qsd_154, qsd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * osd_108[k]
                   + f_3 * pc_z[k] * qsd_150[k];

        t_253[k] = f_17 * osd_153[k]
                   + f_3 * pc_x[k] * qsd_153[k];

        t_254[k] = f_17 * osd_154[k]
                   + f_3 * pc_x[k] * qsd_154[k];

        t_255[k] = f_17 * osd_155[k]
                   + f_3 * pc_x[k] * qsd_155[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.0 / q;

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
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_200 = buffer.data(osf0 + 200);
    const auto *osf0_209 = buffer.data(osf0 + 209);
    const auto *osf0_210 = buffer.data(osf0 + 210);
    const auto *osf0_216 = buffer.data(osf0 + 216);
    const auto *osf0_270 = buffer.data(osf0 + 270);
    const auto *osf0_279 = buffer.data(osf0 + 279);
    const auto *osf0_280 = buffer.data(osf0 + 280);
    const auto *osf0_286 = buffer.data(osf0 + 286);

    const auto *osd_111 = buffer.data(osd + 111);
    const auto *osd_113 = buffer.data(osd + 113);
    const auto *osd_114 = buffer.data(osd + 114);
    const auto *osd_117 = buffer.data(osd + 117);
    const auto *osd_119 = buffer.data(osd + 119);
    const auto *osd_120 = buffer.data(osd + 120);
    const auto *osd_123 = buffer.data(osd + 123);
    const auto *osd_125 = buffer.data(osd + 125);
    const auto *osd_126 = buffer.data(osd + 126);
    const auto *osd_129 = buffer.data(osd + 129);
    const auto *osd_131 = buffer.data(osd + 131);
    const auto *osd_132 = buffer.data(osd + 132);
    const auto *osd_135 = buffer.data(osd + 135);
    const auto *osd_137 = buffer.data(osd + 137);
    const auto *osd_138 = buffer.data(osd + 138);
    const auto *osd_141 = buffer.data(osd + 141);
    const auto *osd_143 = buffer.data(osd + 143);
    const auto *osd_144 = buffer.data(osd + 144);
    const auto *osd_147 = buffer.data(osd + 147);
    const auto *osd_149 = buffer.data(osd + 149);
    const auto *osd_150 = buffer.data(osd + 150);
    const auto *osd_153 = buffer.data(osd + 153);
    const auto *osd_155 = buffer.data(osd + 155);
    const auto *osd_156 = buffer.data(osd + 156);
    const auto *osd_159 = buffer.data(osd + 159);
    const auto *osd_160 = buffer.data(osd + 160);
    const auto *osd_161 = buffer.data(osd + 161);
    const auto *osd_162 = buffer.data(osd + 162);
    const auto *osd_165 = buffer.data(osd + 165);
    const auto *osd_167 = buffer.data(osd + 167);
    const auto *osd_168 = buffer.data(osd + 168);
    const auto *osd_171 = buffer.data(osd + 171);
    const auto *osd_173 = buffer.data(osd + 173);
    const auto *osd_174 = buffer.data(osd + 174);
    const auto *osd_177 = buffer.data(osd + 177);
    const auto *osd_178 = buffer.data(osd + 178);
    const auto *osd_179 = buffer.data(osd + 179);
    const auto *osd_180 = buffer.data(osd + 180);
    const auto *osd_183 = buffer.data(osd + 183);
    const auto *osd_184 = buffer.data(osd + 184);
    const auto *osd_185 = buffer.data(osd + 185);
    const auto *osd_186 = buffer.data(osd + 186);
    const auto *osd_189 = buffer.data(osd + 189);
    const auto *osd_190 = buffer.data(osd + 190);
    const auto *osd_191 = buffer.data(osd + 191);
    const auto *osd_192 = buffer.data(osd + 192);
    const auto *osd_195 = buffer.data(osd + 195);
    const auto *osd_196 = buffer.data(osd + 196);
    const auto *osd_197 = buffer.data(osd + 197);
    const auto *osd_198 = buffer.data(osd + 198);
    const auto *osd_201 = buffer.data(osd + 201);
    const auto *osd_202 = buffer.data(osd + 202);
    const auto *osd_203 = buffer.data(osd + 203);
    const auto *osd_207 = buffer.data(osd + 207);
    const auto *osd_208 = buffer.data(osd + 208);
    const auto *osd_209 = buffer.data(osd + 209);
    const auto *osd_210 = buffer.data(osd + 210);
    const auto *osd_213 = buffer.data(osd + 213);
    const auto *osd_215 = buffer.data(osd + 215);
    const auto *osd_216 = buffer.data(osd + 216);
    const auto *osd_219 = buffer.data(osd + 219);
    const auto *osd_221 = buffer.data(osd + 221);
    const auto *osd_225 = buffer.data(osd + 225);
    const auto *osd_226 = buffer.data(osd + 226);
    const auto *osd_227 = buffer.data(osd + 227);

    const auto *osf1_200 = buffer.data(osf1 + 200);
    const auto *osf1_209 = buffer.data(osf1 + 209);
    const auto *osf1_210 = buffer.data(osf1 + 210);
    const auto *osf1_216 = buffer.data(osf1 + 216);
    const auto *osf1_270 = buffer.data(osf1 + 270);
    const auto *osf1_279 = buffer.data(osf1 + 279);
    const auto *osf1_280 = buffer.data(osf1 + 280);
    const auto *osf1_286 = buffer.data(osf1 + 286);

    const auto *qsp0_76 = buffer.data(qsp0 + 76);
    const auto *qsp0_77 = buffer.data(qsp0 + 77);
    const auto *qsp0_79 = buffer.data(qsp0 + 79);
    const auto *qsp0_81 = buffer.data(qsp0 + 81);
    const auto *qsp0_82 = buffer.data(qsp0 + 82);
    const auto *qsp0_83 = buffer.data(qsp0 + 83);
    const auto *qsp0_84 = buffer.data(qsp0 + 84);
    const auto *qsp0_85 = buffer.data(qsp0 + 85);
    const auto *qsp0_86 = buffer.data(qsp0 + 86);
    const auto *qsp0_89 = buffer.data(qsp0 + 89);
    const auto *qsp0_90 = buffer.data(qsp0 + 90);
    const auto *qsp0_91 = buffer.data(qsp0 + 91);
    const auto *qsp0_92 = buffer.data(qsp0 + 92);
    const auto *qsp0_93 = buffer.data(qsp0 + 93);
    const auto *qsp0_94 = buffer.data(qsp0 + 94);
    const auto *qsp0_95 = buffer.data(qsp0 + 95);
    const auto *qsp0_96 = buffer.data(qsp0 + 96);
    const auto *qsp0_97 = buffer.data(qsp0 + 97);
    const auto *qsp0_98 = buffer.data(qsp0 + 98);
    const auto *qsp0_99 = buffer.data(qsp0 + 99);
    const auto *qsp0_100 = buffer.data(qsp0 + 100);
    const auto *qsp0_101 = buffer.data(qsp0 + 101);
    const auto *qsp0_103 = buffer.data(qsp0 + 103);
    const auto *qsp0_105 = buffer.data(qsp0 + 105);
    const auto *qsp0_106 = buffer.data(qsp0 + 106);
    const auto *qsp0_107 = buffer.data(qsp0 + 107);
    const auto *qsp0_108 = buffer.data(qsp0 + 108);
    const auto *qsp0_109 = buffer.data(qsp0 + 109);
    const auto *qsp0_110 = buffer.data(qsp0 + 110);
    const auto *qsp0_113 = buffer.data(qsp0 + 113);

    const auto *qsp1_76 = buffer.data(qsp1 + 76);
    const auto *qsp1_77 = buffer.data(qsp1 + 77);
    const auto *qsp1_79 = buffer.data(qsp1 + 79);
    const auto *qsp1_81 = buffer.data(qsp1 + 81);
    const auto *qsp1_82 = buffer.data(qsp1 + 82);
    const auto *qsp1_83 = buffer.data(qsp1 + 83);
    const auto *qsp1_84 = buffer.data(qsp1 + 84);
    const auto *qsp1_85 = buffer.data(qsp1 + 85);
    const auto *qsp1_86 = buffer.data(qsp1 + 86);
    const auto *qsp1_89 = buffer.data(qsp1 + 89);
    const auto *qsp1_90 = buffer.data(qsp1 + 90);
    const auto *qsp1_91 = buffer.data(qsp1 + 91);
    const auto *qsp1_92 = buffer.data(qsp1 + 92);
    const auto *qsp1_93 = buffer.data(qsp1 + 93);
    const auto *qsp1_94 = buffer.data(qsp1 + 94);
    const auto *qsp1_95 = buffer.data(qsp1 + 95);
    const auto *qsp1_96 = buffer.data(qsp1 + 96);
    const auto *qsp1_97 = buffer.data(qsp1 + 97);
    const auto *qsp1_98 = buffer.data(qsp1 + 98);
    const auto *qsp1_99 = buffer.data(qsp1 + 99);
    const auto *qsp1_100 = buffer.data(qsp1 + 100);
    const auto *qsp1_101 = buffer.data(qsp1 + 101);
    const auto *qsp1_103 = buffer.data(qsp1 + 103);
    const auto *qsp1_105 = buffer.data(qsp1 + 105);
    const auto *qsp1_106 = buffer.data(qsp1 + 106);
    const auto *qsp1_107 = buffer.data(qsp1 + 107);
    const auto *qsp1_108 = buffer.data(qsp1 + 108);
    const auto *qsp1_109 = buffer.data(qsp1 + 109);
    const auto *qsp1_110 = buffer.data(qsp1 + 110);
    const auto *qsp1_113 = buffer.data(qsp1 + 113);

    const auto *qsd_153 = buffer.data(qsd + 153);
    const auto *qsd_155 = buffer.data(qsd + 155);
    const auto *qsd_156 = buffer.data(qsd + 156);
    const auto *qsd_159 = buffer.data(qsd + 159);
    const auto *qsd_160 = buffer.data(qsd + 160);
    const auto *qsd_161 = buffer.data(qsd + 161);
    const auto *qsd_162 = buffer.data(qsd + 162);
    const auto *qsd_164 = buffer.data(qsd + 164);
    const auto *qsd_165 = buffer.data(qsd + 165);
    const auto *qsd_166 = buffer.data(qsd + 166);
    const auto *qsd_167 = buffer.data(qsd + 167);
    const auto *qsd_168 = buffer.data(qsd + 168);
    const auto *qsd_169 = buffer.data(qsd + 169);
    const auto *qsd_171 = buffer.data(qsd + 171);
    const auto *qsd_173 = buffer.data(qsd + 173);
    const auto *qsd_174 = buffer.data(qsd + 174);
    const auto *qsd_177 = buffer.data(qsd + 177);
    const auto *qsd_178 = buffer.data(qsd + 178);
    const auto *qsd_179 = buffer.data(qsd + 179);
    const auto *qsd_180 = buffer.data(qsd + 180);
    const auto *qsd_183 = buffer.data(qsd + 183);
    const auto *qsd_184 = buffer.data(qsd + 184);
    const auto *qsd_185 = buffer.data(qsd + 185);
    const auto *qsd_186 = buffer.data(qsd + 186);
    const auto *qsd_189 = buffer.data(qsd + 189);
    const auto *qsd_190 = buffer.data(qsd + 190);
    const auto *qsd_191 = buffer.data(qsd + 191);
    const auto *qsd_192 = buffer.data(qsd + 192);
    const auto *qsd_195 = buffer.data(qsd + 195);
    const auto *qsd_196 = buffer.data(qsd + 196);
    const auto *qsd_197 = buffer.data(qsd + 197);
    const auto *qsd_198 = buffer.data(qsd + 198);
    const auto *qsd_201 = buffer.data(qsd + 201);
    const auto *qsd_202 = buffer.data(qsd + 202);
    const auto *qsd_203 = buffer.data(qsd + 203);
    const auto *qsd_204 = buffer.data(qsd + 204);
    const auto *qsd_207 = buffer.data(qsd + 207);
    const auto *qsd_208 = buffer.data(qsd + 208);
    const auto *qsd_209 = buffer.data(qsd + 209);
    const auto *qsd_210 = buffer.data(qsd + 210);
    const auto *qsd_212 = buffer.data(qsd + 212);
    const auto *qsd_213 = buffer.data(qsd + 213);
    const auto *qsd_214 = buffer.data(qsd + 214);
    const auto *qsd_215 = buffer.data(qsd + 215);
    const auto *qsd_216 = buffer.data(qsd + 216);
    const auto *qsd_217 = buffer.data(qsd + 217);
    const auto *qsd_219 = buffer.data(qsd + 219);
    const auto *qsd_221 = buffer.data(qsd + 221);
    const auto *qsd_222 = buffer.data(qsd + 222);
    const auto *qsd_225 = buffer.data(qsd + 225);
    const auto *qsd_226 = buffer.data(qsd + 226);
    const auto *qsd_227 = buffer.data(qsd + 227);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_y, pc_z, osd_111, osd_113, osd_117, \
                         osd_119, qsp0_76, qsp0_77, qsp1_76, qsp1_77, qsd_153, \
                         qsd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_10 * osd_117[k]
                   + f_1 * qsp0_76[k]
                   - f_2 * qsp1_76[k]
                   + f_3 * pc_y[k] * qsd_153[k];

        t_257[k] = f_14 * osd_111[k]
                   + f_3 * pc_z[k] * qsd_153[k];

        t_258[k] = f_10 * osd_119[k]
                   + f_3 * pc_y[k] * qsd_155[k];

        t_259[k] = f_14 * osd_113[k]
                   + f_1 * qsp0_77[k]
                   - f_2 * qsp1_77[k]
                   + f_3 * pc_z[k] * qsd_155[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, pc_z, osf0_200, \
                         osd_114, osd_120, osd_159, osf1_200, qsd_156, \
                         qsd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * osf0_200[k]
                   - f_4 * pc_y[k] * osf1_200[k];

        t_261[k] = f_5 * osd_120[k]
                   + f_3 * pc_y[k] * qsd_156[k];

        t_262[k] = f_16 * osd_114[k]
                   + f_3 * pc_z[k] * qsd_156[k];

        t_263[k] = f_17 * osd_159[k]
                   + f_3 * pc_x[k] * qsd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, pc_z, osd_117, osd_123, \
                         osd_160, osd_161, qsp0_79, qsp1_79, qsd_159, qsd_160, \
                         qsd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * osd_160[k]
                   + f_3 * pc_x[k] * qsd_160[k];

        t_265[k] = f_17 * osd_161[k]
                   + f_3 * pc_x[k] * qsd_161[k];

        t_266[k] = f_5 * osd_123[k]
                   + f_1 * qsp0_79[k]
                   - f_2 * qsp1_79[k]
                   + f_3 * pc_y[k] * qsd_159[k];

        t_267[k] = f_16 * osd_117[k]
                   + f_3 * pc_z[k] * qsd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pc_x, pc_y, osf0_209, osd_125, \
                         osd_162, osf1_209, qsp0_81, qsp1_81, qsd_161, \
                         qsd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * osd_125[k]
                   + f_3 * pc_y[k] * qsd_161[k];

        t_269[k] = pa_y[k] * osf0_209[k]
                   - f_4 * pc_y[k] * osf1_209[k];

        t_270[k] = f_17 * osd_162[k]
                   + f_1 * qsp0_81[k]
                   - f_2 * qsp1_81[k]
                   + f_3 * pc_x[k] * qsd_162[k];

        t_271[k] = f_3 * pc_y[k] * qsd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, osd_120, osd_165, \
                         osd_167, qsd_162, qsd_164, qsd_165, qsd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * osd_120[k]
                   + f_3 * pc_z[k] * qsd_162[k];

        t_273[k] = f_17 * osd_165[k]
                   + f_3 * pc_x[k] * qsd_165[k];

        t_274[k] = f_3 * pc_y[k] * qsd_164[k];

        t_275[k] = f_17 * osd_167[k]
                   + f_3 * pc_x[k] * qsd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, osd_125, qsp0_82, qsp0_83, \
                         qsp1_82, qsp1_83, qsd_165, qsd_166, qsd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * qsp0_82[k]
                   - f_2 * qsp1_82[k]
                   + f_3 * pc_y[k] * qsd_165[k];

        t_277[k] = f_7 * qsp0_83[k]
                   - f_8 * qsp1_83[k]
                   + f_3 * pc_y[k] * qsd_166[k];

        t_278[k] = f_3 * pc_y[k] * qsd_167[k];

        t_279[k] = f_17 * osd_125[k]
                   + f_1 * qsp0_83[k]
                   - f_2 * qsp1_83[k]
                   + f_3 * pc_z[k] * qsd_167[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pc_x, pc_y, pc_z, osd_126, \
                         osd_168, osd_171, qsp0_84, qsp1_84, qsd_168, qsd_169, \
                         qsd_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_16 * osd_168[k]
                   + f_1 * qsp0_84[k]
                   - f_2 * qsp1_84[k]
                   + f_3 * pc_x[k] * qsd_168[k];

        t_281[k] = f_15 * osd_126[k]
                   + f_3 * pc_y[k] * qsd_168[k];

        t_282[k] = f_3 * pc_z[k] * qsd_168[k];

        t_283[k] = f_16 * osd_171[k]
                   + f_3 * pc_x[k] * qsd_171[k];

        t_284[k] = f_3 * pc_z[k] * qsd_169[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pc_x, pc_y, pc_z, osd_129, osd_131, \
                         osd_173, qsp0_85, qsp1_85, qsd_171, qsd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_16 * osd_173[k]
                   + f_3 * pc_x[k] * qsd_173[k];

        t_286[k] = f_15 * osd_129[k]
                   + f_1 * qsp0_85[k]
                   - f_2 * qsp1_85[k]
                   + f_3 * pc_y[k] * qsd_171[k];

        t_287[k] = f_3 * pc_z[k] * qsd_171[k];

        t_288[k] = f_15 * osd_131[k]
                   + f_3 * pc_y[k] * qsd_173[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_z, pc_y, pc_z, osf0_210, osd_126, \
                         osd_132, osf1_210, qsp0_86, qsp1_86, qsd_173, \
                         qsd_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_1 * qsp0_86[k]
                   - f_2 * qsp1_86[k]
                   + f_3 * pc_z[k] * qsd_173[k];

        t_290[k] = pa_z[k] * osf0_210[k]
                   - f_4 * pc_z[k] * osf1_210[k];

        t_291[k] = f_17 * osd_132[k]
                   + f_3 * pc_y[k] * qsd_174[k];

        t_292[k] = f_5 * osd_126[k]
                   + f_3 * pc_z[k] * qsd_174[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_z, pc_x, pc_z, osf0_216, osd_177, \
                         osd_178, osd_179, osf1_216, qsd_177, qsd_178, \
                         qsd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_16 * osd_177[k]
                   + f_3 * pc_x[k] * qsd_177[k];

        t_294[k] = f_16 * osd_178[k]
                   + f_3 * pc_x[k] * qsd_178[k];

        t_295[k] = f_16 * osd_179[k]
                   + f_3 * pc_x[k] * qsd_179[k];

        t_296[k] = pa_z[k] * osf0_216[k]
                   - f_4 * pc_z[k] * osf1_216[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, osd_129, osd_131, osd_137, qsp0_89, \
                         qsp1_89, qsd_177, qsd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * osd_129[k]
                   + f_3 * pc_z[k] * qsd_177[k];

        t_298[k] = f_17 * osd_137[k]
                   + f_3 * pc_y[k] * qsd_179[k];

        t_299[k] = f_5 * osd_131[k]
                   + f_1 * qsp0_89[k]
                   - f_2 * qsp1_89[k]
                   + f_3 * pc_z[k] * qsd_179[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, osd_132, osd_138, \
                         osd_180, osd_183, qsp0_90, qsp1_90, qsd_180, \
                         qsd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_16 * osd_180[k]
                   + f_1 * qsp0_90[k]
                   - f_2 * qsp1_90[k]
                   + f_3 * pc_x[k] * qsd_180[k];

        t_301[k] = f_16 * osd_138[k]
                   + f_3 * pc_y[k] * qsd_180[k];

        t_302[k] = f_10 * osd_132[k]
                   + f_3 * pc_z[k] * qsd_180[k];

        t_303[k] = f_16 * osd_183[k]
                   + f_3 * pc_x[k] * qsd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, pc_z, osd_135, osd_141, \
                         osd_184, osd_185, qsp0_91, qsp1_91, qsd_183, qsd_184, \
                         qsd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_16 * osd_184[k]
                   + f_3 * pc_x[k] * qsd_184[k];

        t_305[k] = f_16 * osd_185[k]
                   + f_3 * pc_x[k] * qsd_185[k];

        t_306[k] = f_16 * osd_141[k]
                   + f_1 * qsp0_91[k]
                   - f_2 * qsp1_91[k]
                   + f_3 * pc_y[k] * qsd_183[k];

        t_307[k] = f_10 * osd_135[k]
                   + f_3 * pc_z[k] * qsd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_x, pc_y, pc_z, osd_137, osd_143, osd_186, \
                         qsp0_92, qsp0_93, qsp1_92, qsp1_93, qsd_185, \
                         qsd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_16 * osd_143[k]
                   + f_3 * pc_y[k] * qsd_185[k];

        t_309[k] = f_10 * osd_137[k]
                   + f_1 * qsp0_92[k]
                   - f_2 * qsp1_92[k]
                   + f_3 * pc_z[k] * qsd_185[k];

        t_310[k] = f_16 * osd_186[k]
                   + f_1 * qsp0_93[k]
                   - f_2 * qsp1_93[k]
                   + f_3 * pc_x[k] * qsd_186[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_x, pc_y, pc_z, osd_138, osd_144, \
                         osd_189, osd_190, qsd_186, qsd_189, qsd_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_14 * osd_144[k]
                   + f_3 * pc_y[k] * qsd_186[k];

        t_312[k] = f_12 * osd_138[k]
                   + f_3 * pc_z[k] * qsd_186[k];

        t_313[k] = f_16 * osd_189[k]
                   + f_3 * pc_x[k] * qsd_189[k];

        t_314[k] = f_16 * osd_190[k]
                   + f_3 * pc_x[k] * qsd_190[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, osd_141, osd_147, \
                         osd_149, osd_191, qsp0_94, qsp1_94, qsd_189, \
                         qsd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_16 * osd_191[k]
                   + f_3 * pc_x[k] * qsd_191[k];

        t_316[k] = f_14 * osd_147[k]
                   + f_1 * qsp0_94[k]
                   - f_2 * qsp1_94[k]
                   + f_3 * pc_y[k] * qsd_189[k];

        t_317[k] = f_12 * osd_141[k]
                   + f_3 * pc_z[k] * qsd_189[k];

        t_318[k] = f_14 * osd_149[k]
                   + f_3 * pc_y[k] * qsd_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pc_x, pc_y, pc_z, osd_143, osd_150, osd_192, \
                         qsp0_95, qsp0_96, qsp1_95, qsp1_96, qsd_191, \
                         qsd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_12 * osd_143[k]
                   + f_1 * qsp0_95[k]
                   - f_2 * qsp1_95[k]
                   + f_3 * pc_z[k] * qsd_191[k];

        t_320[k] = f_16 * osd_192[k]
                   + f_1 * qsp0_96[k]
                   - f_2 * qsp1_96[k]
                   + f_3 * pc_x[k] * qsd_192[k];

        t_321[k] = f_12 * osd_150[k]
                   + f_3 * pc_y[k] * qsd_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, pc_z, osd_144, osd_195, osd_196, \
                         osd_197, qsd_192, qsd_195, qsd_196, qsd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * osd_144[k]
                   + f_3 * pc_z[k] * qsd_192[k];

        t_323[k] = f_16 * osd_195[k]
                   + f_3 * pc_x[k] * qsd_195[k];

        t_324[k] = f_16 * osd_196[k]
                   + f_3 * pc_x[k] * qsd_196[k];

        t_325[k] = f_16 * osd_197[k]
                   + f_3 * pc_x[k] * qsd_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_y, pc_z, osd_147, osd_149, osd_153, \
                         osd_155, qsp0_97, qsp0_98, qsp1_97, qsp1_98, qsd_195, \
                         qsd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * osd_153[k]
                   + f_1 * qsp0_97[k]
                   - f_2 * qsp1_97[k]
                   + f_3 * pc_y[k] * qsd_195[k];

        t_327[k] = f_14 * osd_147[k]
                   + f_3 * pc_z[k] * qsd_195[k];

        t_328[k] = f_12 * osd_155[k]
                   + f_3 * pc_y[k] * qsd_197[k];

        t_329[k] = f_14 * osd_149[k]
                   + f_1 * qsp0_98[k]
                   - f_2 * qsp1_98[k]
                   + f_3 * pc_z[k] * qsd_197[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, pc_y, pc_z, osd_150, osd_156, \
                         osd_198, osd_201, qsp0_99, qsp1_99, qsd_198, \
                         qsd_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_16 * osd_198[k]
                   + f_1 * qsp0_99[k]
                   - f_2 * qsp1_99[k]
                   + f_3 * pc_x[k] * qsd_198[k];

        t_331[k] = f_10 * osd_156[k]
                   + f_3 * pc_y[k] * qsd_198[k];

        t_332[k] = f_16 * osd_150[k]
                   + f_3 * pc_z[k] * qsd_198[k];

        t_333[k] = f_16 * osd_201[k]
                   + f_3 * pc_x[k] * qsd_201[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, osd_153, osd_159, \
                         osd_202, osd_203, qsp0_100, qsp1_100, qsd_201, qsd_202, \
                         qsd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_16 * osd_202[k]
                   + f_3 * pc_x[k] * qsd_202[k];

        t_335[k] = f_16 * osd_203[k]
                   + f_3 * pc_x[k] * qsd_203[k];

        t_336[k] = f_10 * osd_159[k]
                   + f_1 * qsp0_100[k]
                   - f_2 * qsp1_100[k]
                   + f_3 * pc_y[k] * qsd_201[k];

        t_337[k] = f_16 * osd_153[k]
                   + f_3 * pc_z[k] * qsd_201[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_y, pc_y, pc_z, osf0_270, osd_155, \
                         osd_161, osd_162, osf1_270, qsp0_101, qsp1_101, qsd_203, \
                         qsd_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_10 * osd_161[k]
                   + f_3 * pc_y[k] * qsd_203[k];

        t_339[k] = f_16 * osd_155[k]
                   + f_1 * qsp0_101[k]
                   - f_2 * qsp1_101[k]
                   + f_3 * pc_z[k] * qsd_203[k];

        t_340[k] = pa_y[k] * osf0_270[k]
                   - f_4 * pc_y[k] * osf1_270[k];

        t_341[k] = f_5 * osd_162[k]
                   + f_3 * pc_y[k] * qsd_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_z, osd_156, osd_207, osd_208, \
                         osd_209, qsd_204, qsd_207, qsd_208, qsd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_17 * osd_156[k]
                   + f_3 * pc_z[k] * qsd_204[k];

        t_343[k] = f_16 * osd_207[k]
                   + f_3 * pc_x[k] * qsd_207[k];

        t_344[k] = f_16 * osd_208[k]
                   + f_3 * pc_x[k] * qsd_208[k];

        t_345[k] = f_16 * osd_209[k]
                   + f_3 * pc_x[k] * qsd_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_y, pc_y, pc_z, osf0_279, osd_159, \
                         osd_165, osd_167, osf1_279, qsp0_103, qsp1_103, qsd_207, \
                         qsd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_5 * osd_165[k]
                   + f_1 * qsp0_103[k]
                   - f_2 * qsp1_103[k]
                   + f_3 * pc_y[k] * qsd_207[k];

        t_347[k] = f_17 * osd_159[k]
                   + f_3 * pc_z[k] * qsd_207[k];

        t_348[k] = f_5 * osd_167[k]
                   + f_3 * pc_y[k] * qsd_209[k];

        t_349[k] = pa_y[k] * osf0_279[k]
                   - f_4 * pc_y[k] * osf1_279[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, pc_x, pc_y, pc_z, osd_162, \
                         osd_210, osd_213, qsp0_105, qsp1_105, qsd_210, qsd_212, \
                         qsd_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_16 * osd_210[k]
                   + f_1 * qsp0_105[k]
                   - f_2 * qsp1_105[k]
                   + f_3 * pc_x[k] * qsd_210[k];

        t_351[k] = f_3 * pc_y[k] * qsd_210[k];

        t_352[k] = f_15 * osd_162[k]
                   + f_3 * pc_z[k] * qsd_210[k];

        t_353[k] = f_16 * osd_213[k]
                   + f_3 * pc_x[k] * qsd_213[k];

        t_354[k] = f_3 * pc_y[k] * qsd_212[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_x, pc_y, osd_215, qsp0_106, qsp0_107, \
                         qsp1_106, qsp1_107, qsd_213, qsd_214, \
                         qsd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_16 * osd_215[k]
                   + f_3 * pc_x[k] * qsd_215[k];

        t_356[k] = f_1 * qsp0_106[k]
                   - f_2 * qsp1_106[k]
                   + f_3 * pc_y[k] * qsd_213[k];

        t_357[k] = f_7 * qsp0_107[k]
                   - f_8 * qsp1_107[k]
                   + f_3 * pc_y[k] * qsd_214[k];

        t_358[k] = f_3 * pc_y[k] * qsd_215[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pc_x, pc_y, pc_z, osd_167, osd_168, \
                         osd_216, qsp0_107, qsp0_108, qsp1_107, qsp1_108, qsd_215, \
                         qsd_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * osd_167[k]
                   + f_1 * qsp0_107[k]
                   - f_2 * qsp1_107[k]
                   + f_3 * pc_z[k] * qsd_215[k];

        t_360[k] = f_14 * osd_216[k]
                   + f_1 * qsp0_108[k]
                   - f_2 * qsp1_108[k]
                   + f_3 * pc_x[k] * qsd_216[k];

        t_361[k] = f_13 * osd_168[k]
                   + f_3 * pc_y[k] * qsd_216[k];

        t_362[k] = f_3 * pc_z[k] * qsd_216[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, pc_x, pc_y, pc_z, osd_171, \
                         osd_219, osd_221, qsp0_109, qsp1_109, qsd_217, qsd_219, \
                         qsd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * osd_219[k]
                   + f_3 * pc_x[k] * qsd_219[k];

        t_364[k] = f_3 * pc_z[k] * qsd_217[k];

        t_365[k] = f_14 * osd_221[k]
                   + f_3 * pc_x[k] * qsd_221[k];

        t_366[k] = f_13 * osd_171[k]
                   + f_1 * qsp0_109[k]
                   - f_2 * qsp1_109[k]
                   + f_3 * pc_y[k] * qsd_219[k];

        t_367[k] = f_3 * pc_z[k] * qsd_219[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pa_z, pc_y, pc_z, osf0_280, osd_173, \
                         osd_174, osf1_280, qsp0_110, qsp1_110, qsd_221, \
                         qsd_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_13 * osd_173[k]
                   + f_3 * pc_y[k] * qsd_221[k];

        t_369[k] = f_1 * qsp0_110[k]
                   - f_2 * qsp1_110[k]
                   + f_3 * pc_z[k] * qsd_221[k];

        t_370[k] = pa_z[k] * osf0_280[k]
                   - f_4 * pc_z[k] * osf1_280[k];

        t_371[k] = f_15 * osd_174[k]
                   + f_3 * pc_y[k] * qsd_222[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pc_x, pc_z, osd_168, osd_225, osd_226, \
                         osd_227, qsd_222, qsd_225, qsd_226, qsd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_5 * osd_168[k]
                   + f_3 * pc_z[k] * qsd_222[k];

        t_373[k] = f_14 * osd_225[k]
                   + f_3 * pc_x[k] * qsd_225[k];

        t_374[k] = f_14 * osd_226[k]
                   + f_3 * pc_x[k] * qsd_226[k];

        t_375[k] = f_14 * osd_227[k]
                   + f_3 * pc_x[k] * qsd_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pc_y, pc_z, osf0_286, osd_171, \
                         osd_173, osd_179, osf1_286, qsp0_113, qsp1_113, qsd_225, \
                         qsd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_z[k] * osf0_286[k]
                   - f_4 * pc_z[k] * osf1_286[k];

        t_377[k] = f_5 * osd_171[k]
                   + f_3 * pc_z[k] * qsd_225[k];

        t_378[k] = f_15 * osd_179[k]
                   + f_3 * pc_y[k] * qsd_227[k];

        t_379[k] = f_5 * osd_173[k]
                   + f_1 * qsp0_113[k]
                   - f_2 * qsp1_113[k]
                   + f_3 * pc_z[k] * qsd_227[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_350 = buffer.data(osf0 + 350);
    const auto *osf0_359 = buffer.data(osf0 + 359);
    const auto *osf0_360 = buffer.data(osf0 + 360);
    const auto *osf0_366 = buffer.data(osf0 + 366);

    const auto *osd_174 = buffer.data(osd + 174);
    const auto *osd_177 = buffer.data(osd + 177);
    const auto *osd_179 = buffer.data(osd + 179);
    const auto *osd_180 = buffer.data(osd + 180);
    const auto *osd_183 = buffer.data(osd + 183);
    const auto *osd_185 = buffer.data(osd + 185);
    const auto *osd_186 = buffer.data(osd + 186);
    const auto *osd_189 = buffer.data(osd + 189);
    const auto *osd_191 = buffer.data(osd + 191);
    const auto *osd_192 = buffer.data(osd + 192);
    const auto *osd_195 = buffer.data(osd + 195);
    const auto *osd_197 = buffer.data(osd + 197);
    const auto *osd_198 = buffer.data(osd + 198);
    const auto *osd_201 = buffer.data(osd + 201);
    const auto *osd_203 = buffer.data(osd + 203);
    const auto *osd_204 = buffer.data(osd + 204);
    const auto *osd_207 = buffer.data(osd + 207);
    const auto *osd_209 = buffer.data(osd + 209);
    const auto *osd_210 = buffer.data(osd + 210);
    const auto *osd_213 = buffer.data(osd + 213);
    const auto *osd_215 = buffer.data(osd + 215);
    const auto *osd_216 = buffer.data(osd + 216);
    const auto *osd_219 = buffer.data(osd + 219);
    const auto *osd_221 = buffer.data(osd + 221);
    const auto *osd_222 = buffer.data(osd + 222);
    const auto *osd_225 = buffer.data(osd + 225);
    const auto *osd_227 = buffer.data(osd + 227);
    const auto *osd_228 = buffer.data(osd + 228);
    const auto *osd_231 = buffer.data(osd + 231);
    const auto *osd_232 = buffer.data(osd + 232);
    const auto *osd_233 = buffer.data(osd + 233);
    const auto *osd_234 = buffer.data(osd + 234);
    const auto *osd_237 = buffer.data(osd + 237);
    const auto *osd_238 = buffer.data(osd + 238);
    const auto *osd_239 = buffer.data(osd + 239);
    const auto *osd_240 = buffer.data(osd + 240);
    const auto *osd_243 = buffer.data(osd + 243);
    const auto *osd_244 = buffer.data(osd + 244);
    const auto *osd_245 = buffer.data(osd + 245);
    const auto *osd_246 = buffer.data(osd + 246);
    const auto *osd_249 = buffer.data(osd + 249);
    const auto *osd_250 = buffer.data(osd + 250);
    const auto *osd_251 = buffer.data(osd + 251);
    const auto *osd_252 = buffer.data(osd + 252);
    const auto *osd_255 = buffer.data(osd + 255);
    const auto *osd_256 = buffer.data(osd + 256);
    const auto *osd_257 = buffer.data(osd + 257);
    const auto *osd_261 = buffer.data(osd + 261);
    const auto *osd_262 = buffer.data(osd + 262);
    const auto *osd_263 = buffer.data(osd + 263);
    const auto *osd_264 = buffer.data(osd + 264);
    const auto *osd_267 = buffer.data(osd + 267);
    const auto *osd_269 = buffer.data(osd + 269);
    const auto *osd_270 = buffer.data(osd + 270);
    const auto *osd_273 = buffer.data(osd + 273);
    const auto *osd_275 = buffer.data(osd + 275);
    const auto *osd_279 = buffer.data(osd + 279);
    const auto *osd_280 = buffer.data(osd + 280);
    const auto *osd_281 = buffer.data(osd + 281);
    const auto *osd_282 = buffer.data(osd + 282);
    const auto *osd_285 = buffer.data(osd + 285);
    const auto *osd_286 = buffer.data(osd + 286);
    const auto *osd_287 = buffer.data(osd + 287);
    const auto *osd_288 = buffer.data(osd + 288);
    const auto *osd_291 = buffer.data(osd + 291);
    const auto *osd_292 = buffer.data(osd + 292);
    const auto *osd_293 = buffer.data(osd + 293);
    const auto *osd_294 = buffer.data(osd + 294);
    const auto *osd_297 = buffer.data(osd + 297);
    const auto *osd_298 = buffer.data(osd + 298);
    const auto *osd_299 = buffer.data(osd + 299);

    const auto *osf1_350 = buffer.data(osf1 + 350);
    const auto *osf1_359 = buffer.data(osf1 + 359);
    const auto *osf1_360 = buffer.data(osf1 + 360);
    const auto *osf1_366 = buffer.data(osf1 + 366);

    const auto *qsp0_114 = buffer.data(qsp0 + 114);
    const auto *qsp0_115 = buffer.data(qsp0 + 115);
    const auto *qsp0_116 = buffer.data(qsp0 + 116);
    const auto *qsp0_117 = buffer.data(qsp0 + 117);
    const auto *qsp0_118 = buffer.data(qsp0 + 118);
    const auto *qsp0_119 = buffer.data(qsp0 + 119);
    const auto *qsp0_120 = buffer.data(qsp0 + 120);
    const auto *qsp0_121 = buffer.data(qsp0 + 121);
    const auto *qsp0_122 = buffer.data(qsp0 + 122);
    const auto *qsp0_123 = buffer.data(qsp0 + 123);
    const auto *qsp0_124 = buffer.data(qsp0 + 124);
    const auto *qsp0_125 = buffer.data(qsp0 + 125);
    const auto *qsp0_126 = buffer.data(qsp0 + 126);
    const auto *qsp0_127 = buffer.data(qsp0 + 127);
    const auto *qsp0_128 = buffer.data(qsp0 + 128);
    const auto *qsp0_130 = buffer.data(qsp0 + 130);
    const auto *qsp0_132 = buffer.data(qsp0 + 132);
    const auto *qsp0_133 = buffer.data(qsp0 + 133);
    const auto *qsp0_134 = buffer.data(qsp0 + 134);
    const auto *qsp0_135 = buffer.data(qsp0 + 135);
    const auto *qsp0_136 = buffer.data(qsp0 + 136);
    const auto *qsp0_137 = buffer.data(qsp0 + 137);
    const auto *qsp0_140 = buffer.data(qsp0 + 140);
    const auto *qsp0_141 = buffer.data(qsp0 + 141);
    const auto *qsp0_142 = buffer.data(qsp0 + 142);
    const auto *qsp0_143 = buffer.data(qsp0 + 143);
    const auto *qsp0_144 = buffer.data(qsp0 + 144);
    const auto *qsp0_145 = buffer.data(qsp0 + 145);
    const auto *qsp0_146 = buffer.data(qsp0 + 146);
    const auto *qsp0_147 = buffer.data(qsp0 + 147);
    const auto *qsp0_148 = buffer.data(qsp0 + 148);
    const auto *qsp0_149 = buffer.data(qsp0 + 149);

    const auto *qsp1_114 = buffer.data(qsp1 + 114);
    const auto *qsp1_115 = buffer.data(qsp1 + 115);
    const auto *qsp1_116 = buffer.data(qsp1 + 116);
    const auto *qsp1_117 = buffer.data(qsp1 + 117);
    const auto *qsp1_118 = buffer.data(qsp1 + 118);
    const auto *qsp1_119 = buffer.data(qsp1 + 119);
    const auto *qsp1_120 = buffer.data(qsp1 + 120);
    const auto *qsp1_121 = buffer.data(qsp1 + 121);
    const auto *qsp1_122 = buffer.data(qsp1 + 122);
    const auto *qsp1_123 = buffer.data(qsp1 + 123);
    const auto *qsp1_124 = buffer.data(qsp1 + 124);
    const auto *qsp1_125 = buffer.data(qsp1 + 125);
    const auto *qsp1_126 = buffer.data(qsp1 + 126);
    const auto *qsp1_127 = buffer.data(qsp1 + 127);
    const auto *qsp1_128 = buffer.data(qsp1 + 128);
    const auto *qsp1_130 = buffer.data(qsp1 + 130);
    const auto *qsp1_132 = buffer.data(qsp1 + 132);
    const auto *qsp1_133 = buffer.data(qsp1 + 133);
    const auto *qsp1_134 = buffer.data(qsp1 + 134);
    const auto *qsp1_135 = buffer.data(qsp1 + 135);
    const auto *qsp1_136 = buffer.data(qsp1 + 136);
    const auto *qsp1_137 = buffer.data(qsp1 + 137);
    const auto *qsp1_140 = buffer.data(qsp1 + 140);
    const auto *qsp1_141 = buffer.data(qsp1 + 141);
    const auto *qsp1_142 = buffer.data(qsp1 + 142);
    const auto *qsp1_143 = buffer.data(qsp1 + 143);
    const auto *qsp1_144 = buffer.data(qsp1 + 144);
    const auto *qsp1_145 = buffer.data(qsp1 + 145);
    const auto *qsp1_146 = buffer.data(qsp1 + 146);
    const auto *qsp1_147 = buffer.data(qsp1 + 147);
    const auto *qsp1_148 = buffer.data(qsp1 + 148);
    const auto *qsp1_149 = buffer.data(qsp1 + 149);

    const auto *qsd_228 = buffer.data(qsd + 228);
    const auto *qsd_231 = buffer.data(qsd + 231);
    const auto *qsd_232 = buffer.data(qsd + 232);
    const auto *qsd_233 = buffer.data(qsd + 233);
    const auto *qsd_234 = buffer.data(qsd + 234);
    const auto *qsd_237 = buffer.data(qsd + 237);
    const auto *qsd_238 = buffer.data(qsd + 238);
    const auto *qsd_239 = buffer.data(qsd + 239);
    const auto *qsd_240 = buffer.data(qsd + 240);
    const auto *qsd_243 = buffer.data(qsd + 243);
    const auto *qsd_244 = buffer.data(qsd + 244);
    const auto *qsd_245 = buffer.data(qsd + 245);
    const auto *qsd_246 = buffer.data(qsd + 246);
    const auto *qsd_249 = buffer.data(qsd + 249);
    const auto *qsd_250 = buffer.data(qsd + 250);
    const auto *qsd_251 = buffer.data(qsd + 251);
    const auto *qsd_252 = buffer.data(qsd + 252);
    const auto *qsd_255 = buffer.data(qsd + 255);
    const auto *qsd_256 = buffer.data(qsd + 256);
    const auto *qsd_257 = buffer.data(qsd + 257);
    const auto *qsd_258 = buffer.data(qsd + 258);
    const auto *qsd_261 = buffer.data(qsd + 261);
    const auto *qsd_262 = buffer.data(qsd + 262);
    const auto *qsd_263 = buffer.data(qsd + 263);
    const auto *qsd_264 = buffer.data(qsd + 264);
    const auto *qsd_266 = buffer.data(qsd + 266);
    const auto *qsd_267 = buffer.data(qsd + 267);
    const auto *qsd_268 = buffer.data(qsd + 268);
    const auto *qsd_269 = buffer.data(qsd + 269);
    const auto *qsd_270 = buffer.data(qsd + 270);
    const auto *qsd_271 = buffer.data(qsd + 271);
    const auto *qsd_273 = buffer.data(qsd + 273);
    const auto *qsd_275 = buffer.data(qsd + 275);
    const auto *qsd_276 = buffer.data(qsd + 276);
    const auto *qsd_279 = buffer.data(qsd + 279);
    const auto *qsd_280 = buffer.data(qsd + 280);
    const auto *qsd_281 = buffer.data(qsd + 281);
    const auto *qsd_282 = buffer.data(qsd + 282);
    const auto *qsd_285 = buffer.data(qsd + 285);
    const auto *qsd_286 = buffer.data(qsd + 286);
    const auto *qsd_287 = buffer.data(qsd + 287);
    const auto *qsd_288 = buffer.data(qsd + 288);
    const auto *qsd_291 = buffer.data(qsd + 291);
    const auto *qsd_292 = buffer.data(qsd + 292);
    const auto *qsd_293 = buffer.data(qsd + 293);
    const auto *qsd_294 = buffer.data(qsd + 294);
    const auto *qsd_297 = buffer.data(qsd + 297);
    const auto *qsd_298 = buffer.data(qsd + 298);
    const auto *qsd_299 = buffer.data(qsd + 299);

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_y, pc_z, osd_174, osd_180, \
                         osd_228, osd_231, qsp0_114, qsp1_114, qsd_228, \
                         qsd_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_14 * osd_228[k]
                   + f_1 * qsp0_114[k]
                   - f_2 * qsp1_114[k]
                   + f_3 * pc_x[k] * qsd_228[k];

        t_381[k] = f_17 * osd_180[k]
                   + f_3 * pc_y[k] * qsd_228[k];

        t_382[k] = f_10 * osd_174[k]
                   + f_3 * pc_z[k] * qsd_228[k];

        t_383[k] = f_14 * osd_231[k]
                   + f_3 * pc_x[k] * qsd_231[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pc_x, pc_y, pc_z, osd_177, osd_183, \
                         osd_232, osd_233, qsp0_115, qsp1_115, qsd_231, qsd_232, \
                         qsd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_14 * osd_232[k]
                   + f_3 * pc_x[k] * qsd_232[k];

        t_385[k] = f_14 * osd_233[k]
                   + f_3 * pc_x[k] * qsd_233[k];

        t_386[k] = f_17 * osd_183[k]
                   + f_1 * qsp0_115[k]
                   - f_2 * qsp1_115[k]
                   + f_3 * pc_y[k] * qsd_231[k];

        t_387[k] = f_10 * osd_177[k]
                   + f_3 * pc_z[k] * qsd_231[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_x, pc_y, pc_z, osd_179, osd_185, osd_234, \
                         qsp0_116, qsp0_117, qsp1_116, qsp1_117, qsd_233, \
                         qsd_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_17 * osd_185[k]
                   + f_3 * pc_y[k] * qsd_233[k];

        t_389[k] = f_10 * osd_179[k]
                   + f_1 * qsp0_116[k]
                   - f_2 * qsp1_116[k]
                   + f_3 * pc_z[k] * qsd_233[k];

        t_390[k] = f_14 * osd_234[k]
                   + f_1 * qsp0_117[k]
                   - f_2 * qsp1_117[k]
                   + f_3 * pc_x[k] * qsd_234[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, osd_180, osd_186, \
                         osd_237, osd_238, qsd_234, qsd_237, qsd_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * osd_186[k]
                   + f_3 * pc_y[k] * qsd_234[k];

        t_392[k] = f_12 * osd_180[k]
                   + f_3 * pc_z[k] * qsd_234[k];

        t_393[k] = f_14 * osd_237[k]
                   + f_3 * pc_x[k] * qsd_237[k];

        t_394[k] = f_14 * osd_238[k]
                   + f_3 * pc_x[k] * qsd_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, osd_183, osd_189, \
                         osd_191, osd_239, qsp0_118, qsp1_118, qsd_237, \
                         qsd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_14 * osd_239[k]
                   + f_3 * pc_x[k] * qsd_239[k];

        t_396[k] = f_16 * osd_189[k]
                   + f_1 * qsp0_118[k]
                   - f_2 * qsp1_118[k]
                   + f_3 * pc_y[k] * qsd_237[k];

        t_397[k] = f_12 * osd_183[k]
                   + f_3 * pc_z[k] * qsd_237[k];

        t_398[k] = f_16 * osd_191[k]
                   + f_3 * pc_y[k] * qsd_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, osd_185, osd_192, osd_240, \
                         qsp0_119, qsp0_120, qsp1_119, qsp1_120, qsd_239, \
                         qsd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_12 * osd_185[k]
                   + f_1 * qsp0_119[k]
                   - f_2 * qsp1_119[k]
                   + f_3 * pc_z[k] * qsd_239[k];

        t_400[k] = f_14 * osd_240[k]
                   + f_1 * qsp0_120[k]
                   - f_2 * qsp1_120[k]
                   + f_3 * pc_x[k] * qsd_240[k];

        t_401[k] = f_14 * osd_192[k]
                   + f_3 * pc_y[k] * qsd_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, osd_186, osd_243, osd_244, \
                         osd_245, qsd_240, qsd_243, qsd_244, qsd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_14 * osd_186[k]
                   + f_3 * pc_z[k] * qsd_240[k];

        t_403[k] = f_14 * osd_243[k]
                   + f_3 * pc_x[k] * qsd_243[k];

        t_404[k] = f_14 * osd_244[k]
                   + f_3 * pc_x[k] * qsd_244[k];

        t_405[k] = f_14 * osd_245[k]
                   + f_3 * pc_x[k] * qsd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pc_y, pc_z, osd_189, osd_191, osd_195, \
                         osd_197, qsp0_121, qsp0_122, qsp1_121, qsp1_122, qsd_243, \
                         qsd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_14 * osd_195[k]
                   + f_1 * qsp0_121[k]
                   - f_2 * qsp1_121[k]
                   + f_3 * pc_y[k] * qsd_243[k];

        t_407[k] = f_14 * osd_189[k]
                   + f_3 * pc_z[k] * qsd_243[k];

        t_408[k] = f_14 * osd_197[k]
                   + f_3 * pc_y[k] * qsd_245[k];

        t_409[k] = f_14 * osd_191[k]
                   + f_1 * qsp0_122[k]
                   - f_2 * qsp1_122[k]
                   + f_3 * pc_z[k] * qsd_245[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, pc_z, osd_192, osd_198, \
                         osd_246, osd_249, qsp0_123, qsp1_123, qsd_246, \
                         qsd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_14 * osd_246[k]
                   + f_1 * qsp0_123[k]
                   - f_2 * qsp1_123[k]
                   + f_3 * pc_x[k] * qsd_246[k];

        t_411[k] = f_12 * osd_198[k]
                   + f_3 * pc_y[k] * qsd_246[k];

        t_412[k] = f_16 * osd_192[k]
                   + f_3 * pc_z[k] * qsd_246[k];

        t_413[k] = f_14 * osd_249[k]
                   + f_3 * pc_x[k] * qsd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, osd_195, osd_201, \
                         osd_250, osd_251, qsp0_124, qsp1_124, qsd_249, qsd_250, \
                         qsd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_14 * osd_250[k]
                   + f_3 * pc_x[k] * qsd_250[k];

        t_415[k] = f_14 * osd_251[k]
                   + f_3 * pc_x[k] * qsd_251[k];

        t_416[k] = f_12 * osd_201[k]
                   + f_1 * qsp0_124[k]
                   - f_2 * qsp1_124[k]
                   + f_3 * pc_y[k] * qsd_249[k];

        t_417[k] = f_16 * osd_195[k]
                   + f_3 * pc_z[k] * qsd_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pc_x, pc_y, pc_z, osd_197, osd_203, osd_252, \
                         qsp0_125, qsp0_126, qsp1_125, qsp1_126, qsd_251, \
                         qsd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_12 * osd_203[k]
                   + f_3 * pc_y[k] * qsd_251[k];

        t_419[k] = f_16 * osd_197[k]
                   + f_1 * qsp0_125[k]
                   - f_2 * qsp1_125[k]
                   + f_3 * pc_z[k] * qsd_251[k];

        t_420[k] = f_14 * osd_252[k]
                   + f_1 * qsp0_126[k]
                   - f_2 * qsp1_126[k]
                   + f_3 * pc_x[k] * qsd_252[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, osd_198, osd_204, \
                         osd_255, osd_256, qsd_252, qsd_255, qsd_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_10 * osd_204[k]
                   + f_3 * pc_y[k] * qsd_252[k];

        t_422[k] = f_17 * osd_198[k]
                   + f_3 * pc_z[k] * qsd_252[k];

        t_423[k] = f_14 * osd_255[k]
                   + f_3 * pc_x[k] * qsd_255[k];

        t_424[k] = f_14 * osd_256[k]
                   + f_3 * pc_x[k] * qsd_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, osd_201, osd_207, \
                         osd_209, osd_257, qsp0_127, qsp1_127, qsd_255, \
                         qsd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_14 * osd_257[k]
                   + f_3 * pc_x[k] * qsd_257[k];

        t_426[k] = f_10 * osd_207[k]
                   + f_1 * qsp0_127[k]
                   - f_2 * qsp1_127[k]
                   + f_3 * pc_y[k] * qsd_255[k];

        t_427[k] = f_17 * osd_201[k]
                   + f_3 * pc_z[k] * qsd_255[k];

        t_428[k] = f_10 * osd_209[k]
                   + f_3 * pc_y[k] * qsd_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_y, pc_y, pc_z, osf0_350, osd_203, \
                         osd_204, osd_210, osf1_350, qsp0_128, qsp1_128, qsd_257, \
                         qsd_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_17 * osd_203[k]
                   + f_1 * qsp0_128[k]
                   - f_2 * qsp1_128[k]
                   + f_3 * pc_z[k] * qsd_257[k];

        t_430[k] = pa_y[k] * osf0_350[k]
                   - f_4 * pc_y[k] * osf1_350[k];

        t_431[k] = f_5 * osd_210[k]
                   + f_3 * pc_y[k] * qsd_258[k];

        t_432[k] = f_15 * osd_204[k]
                   + f_3 * pc_z[k] * qsd_258[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, osd_213, osd_261, osd_262, \
                         osd_263, qsp0_130, qsp1_130, qsd_261, qsd_262, \
                         qsd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * osd_261[k]
                   + f_3 * pc_x[k] * qsd_261[k];

        t_434[k] = f_14 * osd_262[k]
                   + f_3 * pc_x[k] * qsd_262[k];

        t_435[k] = f_14 * osd_263[k]
                   + f_3 * pc_x[k] * qsd_263[k];

        t_436[k] = f_5 * osd_213[k]
                   + f_1 * qsp0_130[k]
                   - f_2 * qsp1_130[k]
                   + f_3 * pc_y[k] * qsd_261[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_y, pc_y, pc_z, osf0_359, osd_207, osd_215, \
                         osf1_359, qsd_261, qsd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_15 * osd_207[k]
                   + f_3 * pc_z[k] * qsd_261[k];

        t_438[k] = f_5 * osd_215[k]
                   + f_3 * pc_y[k] * qsd_263[k];

        t_439[k] = pa_y[k] * osf0_359[k]
                   - f_4 * pc_y[k] * osf1_359[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, osd_210, \
                         osd_264, osd_267, qsp0_132, qsp1_132, qsd_264, qsd_266, \
                         qsd_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * osd_264[k]
                   + f_1 * qsp0_132[k]
                   - f_2 * qsp1_132[k]
                   + f_3 * pc_x[k] * qsd_264[k];

        t_441[k] = f_3 * pc_y[k] * qsd_264[k];

        t_442[k] = f_13 * osd_210[k]
                   + f_3 * pc_z[k] * qsd_264[k];

        t_443[k] = f_14 * osd_267[k]
                   + f_3 * pc_x[k] * qsd_267[k];

        t_444[k] = f_3 * pc_y[k] * qsd_266[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_y, osd_269, qsp0_133, qsp0_134, \
                         qsp1_133, qsp1_134, qsd_267, qsd_268, \
                         qsd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_14 * osd_269[k]
                   + f_3 * pc_x[k] * qsd_269[k];

        t_446[k] = f_1 * qsp0_133[k]
                   - f_2 * qsp1_133[k]
                   + f_3 * pc_y[k] * qsd_267[k];

        t_447[k] = f_7 * qsp0_134[k]
                   - f_8 * qsp1_134[k]
                   + f_3 * pc_y[k] * qsd_268[k];

        t_448[k] = f_3 * pc_y[k] * qsd_269[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, osd_215, osd_216, \
                         osd_270, qsp0_134, qsp0_135, qsp1_134, qsp1_135, qsd_269, \
                         qsd_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_13 * osd_215[k]
                   + f_1 * qsp0_134[k]
                   - f_2 * qsp1_134[k]
                   + f_3 * pc_z[k] * qsd_269[k];

        t_450[k] = f_12 * osd_270[k]
                   + f_1 * qsp0_135[k]
                   - f_2 * qsp1_135[k]
                   + f_3 * pc_x[k] * qsd_270[k];

        t_451[k] = f_11 * osd_216[k]
                   + f_3 * pc_y[k] * qsd_270[k];

        t_452[k] = f_3 * pc_z[k] * qsd_270[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, pc_x, pc_y, pc_z, osd_219, \
                         osd_273, osd_275, qsp0_136, qsp1_136, qsd_271, qsd_273, \
                         qsd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_12 * osd_273[k]
                   + f_3 * pc_x[k] * qsd_273[k];

        t_454[k] = f_3 * pc_z[k] * qsd_271[k];

        t_455[k] = f_12 * osd_275[k]
                   + f_3 * pc_x[k] * qsd_275[k];

        t_456[k] = f_11 * osd_219[k]
                   + f_1 * qsp0_136[k]
                   - f_2 * qsp1_136[k]
                   + f_3 * pc_y[k] * qsd_273[k];

        t_457[k] = f_3 * pc_z[k] * qsd_273[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pa_z, pc_y, pc_z, osf0_360, osd_221, \
                         osd_222, osf1_360, qsp0_137, qsp1_137, qsd_275, \
                         qsd_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_11 * osd_221[k]
                   + f_3 * pc_y[k] * qsd_275[k];

        t_459[k] = f_1 * qsp0_137[k]
                   - f_2 * qsp1_137[k]
                   + f_3 * pc_z[k] * qsd_275[k];

        t_460[k] = pa_z[k] * osf0_360[k]
                   - f_4 * pc_z[k] * osf1_360[k];

        t_461[k] = f_13 * osd_222[k]
                   + f_3 * pc_y[k] * qsd_276[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, pc_z, osd_216, osd_279, osd_280, \
                         osd_281, qsd_276, qsd_279, qsd_280, qsd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_5 * osd_216[k]
                   + f_3 * pc_z[k] * qsd_276[k];

        t_463[k] = f_12 * osd_279[k]
                   + f_3 * pc_x[k] * qsd_279[k];

        t_464[k] = f_12 * osd_280[k]
                   + f_3 * pc_x[k] * qsd_280[k];

        t_465[k] = f_12 * osd_281[k]
                   + f_3 * pc_x[k] * qsd_281[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_y, pc_z, osf0_366, osd_219, \
                         osd_221, osd_227, osf1_366, qsp0_140, qsp1_140, qsd_279, \
                         qsd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = pa_z[k] * osf0_366[k]
                   - f_4 * pc_z[k] * osf1_366[k];

        t_467[k] = f_5 * osd_219[k]
                   + f_3 * pc_z[k] * qsd_279[k];

        t_468[k] = f_13 * osd_227[k]
                   + f_3 * pc_y[k] * qsd_281[k];

        t_469[k] = f_5 * osd_221[k]
                   + f_1 * qsp0_140[k]
                   - f_2 * qsp1_140[k]
                   + f_3 * pc_z[k] * qsd_281[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, pc_z, osd_222, osd_228, \
                         osd_282, osd_285, qsp0_141, qsp1_141, qsd_282, \
                         qsd_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_12 * osd_282[k]
                   + f_1 * qsp0_141[k]
                   - f_2 * qsp1_141[k]
                   + f_3 * pc_x[k] * qsd_282[k];

        t_471[k] = f_15 * osd_228[k]
                   + f_3 * pc_y[k] * qsd_282[k];

        t_472[k] = f_10 * osd_222[k]
                   + f_3 * pc_z[k] * qsd_282[k];

        t_473[k] = f_12 * osd_285[k]
                   + f_3 * pc_x[k] * qsd_285[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_x, pc_y, pc_z, osd_225, osd_231, \
                         osd_286, osd_287, qsp0_142, qsp1_142, qsd_285, qsd_286, \
                         qsd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_12 * osd_286[k]
                   + f_3 * pc_x[k] * qsd_286[k];

        t_475[k] = f_12 * osd_287[k]
                   + f_3 * pc_x[k] * qsd_287[k];

        t_476[k] = f_15 * osd_231[k]
                   + f_1 * qsp0_142[k]
                   - f_2 * qsp1_142[k]
                   + f_3 * pc_y[k] * qsd_285[k];

        t_477[k] = f_10 * osd_225[k]
                   + f_3 * pc_z[k] * qsd_285[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, osd_227, osd_233, osd_288, \
                         qsp0_143, qsp0_144, qsp1_143, qsp1_144, qsd_287, \
                         qsd_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_15 * osd_233[k]
                   + f_3 * pc_y[k] * qsd_287[k];

        t_479[k] = f_10 * osd_227[k]
                   + f_1 * qsp0_143[k]
                   - f_2 * qsp1_143[k]
                   + f_3 * pc_z[k] * qsd_287[k];

        t_480[k] = f_12 * osd_288[k]
                   + f_1 * qsp0_144[k]
                   - f_2 * qsp1_144[k]
                   + f_3 * pc_x[k] * qsd_288[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, osd_228, osd_234, \
                         osd_291, osd_292, qsd_288, qsd_291, qsd_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_17 * osd_234[k]
                   + f_3 * pc_y[k] * qsd_288[k];

        t_482[k] = f_12 * osd_228[k]
                   + f_3 * pc_z[k] * qsd_288[k];

        t_483[k] = f_12 * osd_291[k]
                   + f_3 * pc_x[k] * qsd_291[k];

        t_484[k] = f_12 * osd_292[k]
                   + f_3 * pc_x[k] * qsd_292[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, pc_y, pc_z, osd_231, osd_237, \
                         osd_239, osd_293, qsp0_145, qsp1_145, qsd_291, \
                         qsd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_12 * osd_293[k]
                   + f_3 * pc_x[k] * qsd_293[k];

        t_486[k] = f_17 * osd_237[k]
                   + f_1 * qsp0_145[k]
                   - f_2 * qsp1_145[k]
                   + f_3 * pc_y[k] * qsd_291[k];

        t_487[k] = f_12 * osd_231[k]
                   + f_3 * pc_z[k] * qsd_291[k];

        t_488[k] = f_17 * osd_239[k]
                   + f_3 * pc_y[k] * qsd_293[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, osd_233, osd_240, osd_294, \
                         qsp0_146, qsp0_147, qsp1_146, qsp1_147, qsd_293, \
                         qsd_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_12 * osd_233[k]
                   + f_1 * qsp0_146[k]
                   - f_2 * qsp1_146[k]
                   + f_3 * pc_z[k] * qsd_293[k];

        t_490[k] = f_12 * osd_294[k]
                   + f_1 * qsp0_147[k]
                   - f_2 * qsp1_147[k]
                   + f_3 * pc_x[k] * qsd_294[k];

        t_491[k] = f_16 * osd_240[k]
                   + f_3 * pc_y[k] * qsd_294[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pc_x, pc_z, osd_234, osd_297, osd_298, \
                         osd_299, qsd_294, qsd_297, qsd_298, qsd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * osd_234[k]
                   + f_3 * pc_z[k] * qsd_294[k];

        t_493[k] = f_12 * osd_297[k]
                   + f_3 * pc_x[k] * qsd_297[k];

        t_494[k] = f_12 * osd_298[k]
                   + f_3 * pc_x[k] * qsd_298[k];

        t_495[k] = f_12 * osd_299[k]
                   + f_3 * pc_x[k] * qsd_299[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_y, pc_z, osd_237, osd_239, osd_243, \
                         osd_245, qsp0_148, qsp0_149, qsp1_148, qsp1_149, qsd_297, \
                         qsd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_16 * osd_243[k]
                   + f_1 * qsp0_148[k]
                   - f_2 * qsp1_148[k]
                   + f_3 * pc_y[k] * qsd_297[k];

        t_497[k] = f_14 * osd_237[k]
                   + f_3 * pc_z[k] * qsd_297[k];

        t_498[k] = f_16 * osd_245[k]
                   + f_3 * pc_y[k] * qsd_299[k];

        t_499[k] = f_14 * osd_239[k]
                   + f_1 * qsp0_149[k]
                   - f_2 * qsp1_149[k]
                   + f_3 * pc_z[k] * qsd_299[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_440 = buffer.data(osf0 + 440);
    const auto *osf0_449 = buffer.data(osf0 + 449);
    const auto *osf0_450 = buffer.data(osf0 + 450);
    const auto *osf0_456 = buffer.data(osf0 + 456);

    const auto *osd_240 = buffer.data(osd + 240);
    const auto *osd_243 = buffer.data(osd + 243);
    const auto *osd_245 = buffer.data(osd + 245);
    const auto *osd_246 = buffer.data(osd + 246);
    const auto *osd_249 = buffer.data(osd + 249);
    const auto *osd_251 = buffer.data(osd + 251);
    const auto *osd_252 = buffer.data(osd + 252);
    const auto *osd_255 = buffer.data(osd + 255);
    const auto *osd_257 = buffer.data(osd + 257);
    const auto *osd_258 = buffer.data(osd + 258);
    const auto *osd_261 = buffer.data(osd + 261);
    const auto *osd_263 = buffer.data(osd + 263);
    const auto *osd_264 = buffer.data(osd + 264);
    const auto *osd_267 = buffer.data(osd + 267);
    const auto *osd_269 = buffer.data(osd + 269);
    const auto *osd_270 = buffer.data(osd + 270);
    const auto *osd_273 = buffer.data(osd + 273);
    const auto *osd_275 = buffer.data(osd + 275);
    const auto *osd_276 = buffer.data(osd + 276);
    const auto *osd_279 = buffer.data(osd + 279);
    const auto *osd_281 = buffer.data(osd + 281);
    const auto *osd_282 = buffer.data(osd + 282);
    const auto *osd_285 = buffer.data(osd + 285);
    const auto *osd_287 = buffer.data(osd + 287);
    const auto *osd_288 = buffer.data(osd + 288);
    const auto *osd_291 = buffer.data(osd + 291);
    const auto *osd_293 = buffer.data(osd + 293);
    const auto *osd_294 = buffer.data(osd + 294);
    const auto *osd_297 = buffer.data(osd + 297);
    const auto *osd_299 = buffer.data(osd + 299);
    const auto *osd_300 = buffer.data(osd + 300);
    const auto *osd_303 = buffer.data(osd + 303);
    const auto *osd_304 = buffer.data(osd + 304);
    const auto *osd_305 = buffer.data(osd + 305);
    const auto *osd_306 = buffer.data(osd + 306);
    const auto *osd_309 = buffer.data(osd + 309);
    const auto *osd_310 = buffer.data(osd + 310);
    const auto *osd_311 = buffer.data(osd + 311);
    const auto *osd_312 = buffer.data(osd + 312);
    const auto *osd_315 = buffer.data(osd + 315);
    const auto *osd_316 = buffer.data(osd + 316);
    const auto *osd_317 = buffer.data(osd + 317);
    const auto *osd_321 = buffer.data(osd + 321);
    const auto *osd_322 = buffer.data(osd + 322);
    const auto *osd_323 = buffer.data(osd + 323);
    const auto *osd_324 = buffer.data(osd + 324);
    const auto *osd_327 = buffer.data(osd + 327);
    const auto *osd_329 = buffer.data(osd + 329);
    const auto *osd_330 = buffer.data(osd + 330);
    const auto *osd_333 = buffer.data(osd + 333);
    const auto *osd_335 = buffer.data(osd + 335);
    const auto *osd_339 = buffer.data(osd + 339);
    const auto *osd_340 = buffer.data(osd + 340);
    const auto *osd_341 = buffer.data(osd + 341);
    const auto *osd_342 = buffer.data(osd + 342);
    const auto *osd_345 = buffer.data(osd + 345);
    const auto *osd_346 = buffer.data(osd + 346);
    const auto *osd_347 = buffer.data(osd + 347);
    const auto *osd_348 = buffer.data(osd + 348);
    const auto *osd_351 = buffer.data(osd + 351);
    const auto *osd_352 = buffer.data(osd + 352);
    const auto *osd_353 = buffer.data(osd + 353);
    const auto *osd_354 = buffer.data(osd + 354);
    const auto *osd_357 = buffer.data(osd + 357);
    const auto *osd_358 = buffer.data(osd + 358);
    const auto *osd_359 = buffer.data(osd + 359);
    const auto *osd_360 = buffer.data(osd + 360);
    const auto *osd_363 = buffer.data(osd + 363);
    const auto *osd_364 = buffer.data(osd + 364);
    const auto *osd_365 = buffer.data(osd + 365);
    const auto *osd_366 = buffer.data(osd + 366);
    const auto *osd_369 = buffer.data(osd + 369);
    const auto *osd_370 = buffer.data(osd + 370);
    const auto *osd_371 = buffer.data(osd + 371);
    const auto *osd_372 = buffer.data(osd + 372);

    const auto *osf1_440 = buffer.data(osf1 + 440);
    const auto *osf1_449 = buffer.data(osf1 + 449);
    const auto *osf1_450 = buffer.data(osf1 + 450);
    const auto *osf1_456 = buffer.data(osf1 + 456);

    const auto *qsp0_150 = buffer.data(qsp0 + 150);
    const auto *qsp0_151 = buffer.data(qsp0 + 151);
    const auto *qsp0_152 = buffer.data(qsp0 + 152);
    const auto *qsp0_153 = buffer.data(qsp0 + 153);
    const auto *qsp0_154 = buffer.data(qsp0 + 154);
    const auto *qsp0_155 = buffer.data(qsp0 + 155);
    const auto *qsp0_156 = buffer.data(qsp0 + 156);
    const auto *qsp0_157 = buffer.data(qsp0 + 157);
    const auto *qsp0_158 = buffer.data(qsp0 + 158);
    const auto *qsp0_160 = buffer.data(qsp0 + 160);
    const auto *qsp0_162 = buffer.data(qsp0 + 162);
    const auto *qsp0_163 = buffer.data(qsp0 + 163);
    const auto *qsp0_164 = buffer.data(qsp0 + 164);
    const auto *qsp0_165 = buffer.data(qsp0 + 165);
    const auto *qsp0_166 = buffer.data(qsp0 + 166);
    const auto *qsp0_167 = buffer.data(qsp0 + 167);
    const auto *qsp0_170 = buffer.data(qsp0 + 170);
    const auto *qsp0_171 = buffer.data(qsp0 + 171);
    const auto *qsp0_172 = buffer.data(qsp0 + 172);
    const auto *qsp0_173 = buffer.data(qsp0 + 173);
    const auto *qsp0_174 = buffer.data(qsp0 + 174);
    const auto *qsp0_175 = buffer.data(qsp0 + 175);
    const auto *qsp0_176 = buffer.data(qsp0 + 176);
    const auto *qsp0_177 = buffer.data(qsp0 + 177);
    const auto *qsp0_178 = buffer.data(qsp0 + 178);
    const auto *qsp0_179 = buffer.data(qsp0 + 179);
    const auto *qsp0_180 = buffer.data(qsp0 + 180);
    const auto *qsp0_181 = buffer.data(qsp0 + 181);
    const auto *qsp0_182 = buffer.data(qsp0 + 182);
    const auto *qsp0_183 = buffer.data(qsp0 + 183);
    const auto *qsp0_184 = buffer.data(qsp0 + 184);
    const auto *qsp0_185 = buffer.data(qsp0 + 185);
    const auto *qsp0_186 = buffer.data(qsp0 + 186);

    const auto *qsp1_150 = buffer.data(qsp1 + 150);
    const auto *qsp1_151 = buffer.data(qsp1 + 151);
    const auto *qsp1_152 = buffer.data(qsp1 + 152);
    const auto *qsp1_153 = buffer.data(qsp1 + 153);
    const auto *qsp1_154 = buffer.data(qsp1 + 154);
    const auto *qsp1_155 = buffer.data(qsp1 + 155);
    const auto *qsp1_156 = buffer.data(qsp1 + 156);
    const auto *qsp1_157 = buffer.data(qsp1 + 157);
    const auto *qsp1_158 = buffer.data(qsp1 + 158);
    const auto *qsp1_160 = buffer.data(qsp1 + 160);
    const auto *qsp1_162 = buffer.data(qsp1 + 162);
    const auto *qsp1_163 = buffer.data(qsp1 + 163);
    const auto *qsp1_164 = buffer.data(qsp1 + 164);
    const auto *qsp1_165 = buffer.data(qsp1 + 165);
    const auto *qsp1_166 = buffer.data(qsp1 + 166);
    const auto *qsp1_167 = buffer.data(qsp1 + 167);
    const auto *qsp1_170 = buffer.data(qsp1 + 170);
    const auto *qsp1_171 = buffer.data(qsp1 + 171);
    const auto *qsp1_172 = buffer.data(qsp1 + 172);
    const auto *qsp1_173 = buffer.data(qsp1 + 173);
    const auto *qsp1_174 = buffer.data(qsp1 + 174);
    const auto *qsp1_175 = buffer.data(qsp1 + 175);
    const auto *qsp1_176 = buffer.data(qsp1 + 176);
    const auto *qsp1_177 = buffer.data(qsp1 + 177);
    const auto *qsp1_178 = buffer.data(qsp1 + 178);
    const auto *qsp1_179 = buffer.data(qsp1 + 179);
    const auto *qsp1_180 = buffer.data(qsp1 + 180);
    const auto *qsp1_181 = buffer.data(qsp1 + 181);
    const auto *qsp1_182 = buffer.data(qsp1 + 182);
    const auto *qsp1_183 = buffer.data(qsp1 + 183);
    const auto *qsp1_184 = buffer.data(qsp1 + 184);
    const auto *qsp1_185 = buffer.data(qsp1 + 185);
    const auto *qsp1_186 = buffer.data(qsp1 + 186);

    const auto *qsd_300 = buffer.data(qsd + 300);
    const auto *qsd_303 = buffer.data(qsd + 303);
    const auto *qsd_304 = buffer.data(qsd + 304);
    const auto *qsd_305 = buffer.data(qsd + 305);
    const auto *qsd_306 = buffer.data(qsd + 306);
    const auto *qsd_309 = buffer.data(qsd + 309);
    const auto *qsd_310 = buffer.data(qsd + 310);
    const auto *qsd_311 = buffer.data(qsd + 311);
    const auto *qsd_312 = buffer.data(qsd + 312);
    const auto *qsd_315 = buffer.data(qsd + 315);
    const auto *qsd_316 = buffer.data(qsd + 316);
    const auto *qsd_317 = buffer.data(qsd + 317);
    const auto *qsd_318 = buffer.data(qsd + 318);
    const auto *qsd_321 = buffer.data(qsd + 321);
    const auto *qsd_322 = buffer.data(qsd + 322);
    const auto *qsd_323 = buffer.data(qsd + 323);
    const auto *qsd_324 = buffer.data(qsd + 324);
    const auto *qsd_326 = buffer.data(qsd + 326);
    const auto *qsd_327 = buffer.data(qsd + 327);
    const auto *qsd_328 = buffer.data(qsd + 328);
    const auto *qsd_329 = buffer.data(qsd + 329);
    const auto *qsd_330 = buffer.data(qsd + 330);
    const auto *qsd_331 = buffer.data(qsd + 331);
    const auto *qsd_333 = buffer.data(qsd + 333);
    const auto *qsd_335 = buffer.data(qsd + 335);
    const auto *qsd_336 = buffer.data(qsd + 336);
    const auto *qsd_339 = buffer.data(qsd + 339);
    const auto *qsd_340 = buffer.data(qsd + 340);
    const auto *qsd_341 = buffer.data(qsd + 341);
    const auto *qsd_342 = buffer.data(qsd + 342);
    const auto *qsd_345 = buffer.data(qsd + 345);
    const auto *qsd_346 = buffer.data(qsd + 346);
    const auto *qsd_347 = buffer.data(qsd + 347);
    const auto *qsd_348 = buffer.data(qsd + 348);
    const auto *qsd_351 = buffer.data(qsd + 351);
    const auto *qsd_352 = buffer.data(qsd + 352);
    const auto *qsd_353 = buffer.data(qsd + 353);
    const auto *qsd_354 = buffer.data(qsd + 354);
    const auto *qsd_357 = buffer.data(qsd + 357);
    const auto *qsd_358 = buffer.data(qsd + 358);
    const auto *qsd_359 = buffer.data(qsd + 359);
    const auto *qsd_360 = buffer.data(qsd + 360);
    const auto *qsd_363 = buffer.data(qsd + 363);
    const auto *qsd_364 = buffer.data(qsd + 364);
    const auto *qsd_365 = buffer.data(qsd + 365);
    const auto *qsd_366 = buffer.data(qsd + 366);
    const auto *qsd_369 = buffer.data(qsd + 369);
    const auto *qsd_370 = buffer.data(qsd + 370);
    const auto *qsd_371 = buffer.data(qsd + 371);
    const auto *qsd_372 = buffer.data(qsd + 372);

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_x, pc_y, pc_z, osd_240, osd_246, \
                         osd_300, osd_303, qsp0_150, qsp1_150, qsd_300, \
                         qsd_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_12 * osd_300[k]
                   + f_1 * qsp0_150[k]
                   - f_2 * qsp1_150[k]
                   + f_3 * pc_x[k] * qsd_300[k];

        t_501[k] = f_14 * osd_246[k]
                   + f_3 * pc_y[k] * qsd_300[k];

        t_502[k] = f_16 * osd_240[k]
                   + f_3 * pc_z[k] * qsd_300[k];

        t_503[k] = f_12 * osd_303[k]
                   + f_3 * pc_x[k] * qsd_303[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, osd_243, osd_249, \
                         osd_304, osd_305, qsp0_151, qsp1_151, qsd_303, qsd_304, \
                         qsd_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_12 * osd_304[k]
                   + f_3 * pc_x[k] * qsd_304[k];

        t_505[k] = f_12 * osd_305[k]
                   + f_3 * pc_x[k] * qsd_305[k];

        t_506[k] = f_14 * osd_249[k]
                   + f_1 * qsp0_151[k]
                   - f_2 * qsp1_151[k]
                   + f_3 * pc_y[k] * qsd_303[k];

        t_507[k] = f_16 * osd_243[k]
                   + f_3 * pc_z[k] * qsd_303[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, pc_z, osd_245, osd_251, osd_306, \
                         qsp0_152, qsp0_153, qsp1_152, qsp1_153, qsd_305, \
                         qsd_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_14 * osd_251[k]
                   + f_3 * pc_y[k] * qsd_305[k];

        t_509[k] = f_16 * osd_245[k]
                   + f_1 * qsp0_152[k]
                   - f_2 * qsp1_152[k]
                   + f_3 * pc_z[k] * qsd_305[k];

        t_510[k] = f_12 * osd_306[k]
                   + f_1 * qsp0_153[k]
                   - f_2 * qsp1_153[k]
                   + f_3 * pc_x[k] * qsd_306[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pc_x, pc_y, pc_z, osd_246, osd_252, \
                         osd_309, osd_310, qsd_306, qsd_309, qsd_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_12 * osd_252[k]
                   + f_3 * pc_y[k] * qsd_306[k];

        t_512[k] = f_17 * osd_246[k]
                   + f_3 * pc_z[k] * qsd_306[k];

        t_513[k] = f_12 * osd_309[k]
                   + f_3 * pc_x[k] * qsd_309[k];

        t_514[k] = f_12 * osd_310[k]
                   + f_3 * pc_x[k] * qsd_310[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pc_x, pc_y, pc_z, osd_249, osd_255, \
                         osd_257, osd_311, qsp0_154, qsp1_154, qsd_309, \
                         qsd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_12 * osd_311[k]
                   + f_3 * pc_x[k] * qsd_311[k];

        t_516[k] = f_12 * osd_255[k]
                   + f_1 * qsp0_154[k]
                   - f_2 * qsp1_154[k]
                   + f_3 * pc_y[k] * qsd_309[k];

        t_517[k] = f_17 * osd_249[k]
                   + f_3 * pc_z[k] * qsd_309[k];

        t_518[k] = f_12 * osd_257[k]
                   + f_3 * pc_y[k] * qsd_311[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, pc_z, osd_251, osd_258, osd_312, \
                         qsp0_155, qsp0_156, qsp1_155, qsp1_156, qsd_311, \
                         qsd_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_17 * osd_251[k]
                   + f_1 * qsp0_155[k]
                   - f_2 * qsp1_155[k]
                   + f_3 * pc_z[k] * qsd_311[k];

        t_520[k] = f_12 * osd_312[k]
                   + f_1 * qsp0_156[k]
                   - f_2 * qsp1_156[k]
                   + f_3 * pc_x[k] * qsd_312[k];

        t_521[k] = f_10 * osd_258[k]
                   + f_3 * pc_y[k] * qsd_312[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_z, osd_252, osd_315, osd_316, \
                         osd_317, qsd_312, qsd_315, qsd_316, qsd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_15 * osd_252[k]
                   + f_3 * pc_z[k] * qsd_312[k];

        t_523[k] = f_12 * osd_315[k]
                   + f_3 * pc_x[k] * qsd_315[k];

        t_524[k] = f_12 * osd_316[k]
                   + f_3 * pc_x[k] * qsd_316[k];

        t_525[k] = f_12 * osd_317[k]
                   + f_3 * pc_x[k] * qsd_317[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_y, pc_z, osd_255, osd_257, osd_261, \
                         osd_263, qsp0_157, qsp0_158, qsp1_157, qsp1_158, qsd_315, \
                         qsd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * osd_261[k]
                   + f_1 * qsp0_157[k]
                   - f_2 * qsp1_157[k]
                   + f_3 * pc_y[k] * qsd_315[k];

        t_527[k] = f_15 * osd_255[k]
                   + f_3 * pc_z[k] * qsd_315[k];

        t_528[k] = f_10 * osd_263[k]
                   + f_3 * pc_y[k] * qsd_317[k];

        t_529[k] = f_15 * osd_257[k]
                   + f_1 * qsp0_158[k]
                   - f_2 * qsp1_158[k]
                   + f_3 * pc_z[k] * qsd_317[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pa_y, pc_x, pc_y, pc_z, osf0_440, \
                         osd_258, osd_264, osd_321, osf1_440, qsd_318, \
                         qsd_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = pa_y[k] * osf0_440[k]
                   - f_4 * pc_y[k] * osf1_440[k];

        t_531[k] = f_5 * osd_264[k]
                   + f_3 * pc_y[k] * qsd_318[k];

        t_532[k] = f_13 * osd_258[k]
                   + f_3 * pc_z[k] * qsd_318[k];

        t_533[k] = f_12 * osd_321[k]
                   + f_3 * pc_x[k] * qsd_321[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pc_x, pc_y, pc_z, osd_261, osd_267, \
                         osd_322, osd_323, qsp0_160, qsp1_160, qsd_321, qsd_322, \
                         qsd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_12 * osd_322[k]
                   + f_3 * pc_x[k] * qsd_322[k];

        t_535[k] = f_12 * osd_323[k]
                   + f_3 * pc_x[k] * qsd_323[k];

        t_536[k] = f_5 * osd_267[k]
                   + f_1 * qsp0_160[k]
                   - f_2 * qsp1_160[k]
                   + f_3 * pc_y[k] * qsd_321[k];

        t_537[k] = f_13 * osd_261[k]
                   + f_3 * pc_z[k] * qsd_321[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_y, pc_x, pc_y, osf0_449, osd_269, \
                         osd_324, osf1_449, qsp0_162, qsp1_162, qsd_323, \
                         qsd_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_5 * osd_269[k]
                   + f_3 * pc_y[k] * qsd_323[k];

        t_539[k] = pa_y[k] * osf0_449[k]
                   - f_4 * pc_y[k] * osf1_449[k];

        t_540[k] = f_12 * osd_324[k]
                   + f_1 * qsp0_162[k]
                   - f_2 * qsp1_162[k]
                   + f_3 * pc_x[k] * qsd_324[k];

        t_541[k] = f_3 * pc_y[k] * qsd_324[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pc_x, pc_y, pc_z, osd_264, osd_327, \
                         osd_329, qsd_324, qsd_326, qsd_327, qsd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_11 * osd_264[k]
                   + f_3 * pc_z[k] * qsd_324[k];

        t_543[k] = f_12 * osd_327[k]
                   + f_3 * pc_x[k] * qsd_327[k];

        t_544[k] = f_3 * pc_y[k] * qsd_326[k];

        t_545[k] = f_12 * osd_329[k]
                   + f_3 * pc_x[k] * qsd_329[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pc_y, pc_z, osd_269, qsp0_163, qsp0_164, \
                         qsp1_163, qsp1_164, qsd_327, qsd_328, \
                         qsd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_1 * qsp0_163[k]
                   - f_2 * qsp1_163[k]
                   + f_3 * pc_y[k] * qsd_327[k];

        t_547[k] = f_7 * qsp0_164[k]
                   - f_8 * qsp1_164[k]
                   + f_3 * pc_y[k] * qsd_328[k];

        t_548[k] = f_3 * pc_y[k] * qsd_329[k];

        t_549[k] = f_11 * osd_269[k]
                   + f_1 * qsp0_164[k]
                   - f_2 * qsp1_164[k]
                   + f_3 * pc_z[k] * qsd_329[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, pc_x, pc_y, pc_z, osd_270, \
                         osd_330, osd_333, qsp0_165, qsp1_165, qsd_330, qsd_331, \
                         qsd_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_10 * osd_330[k]
                   + f_1 * qsp0_165[k]
                   - f_2 * qsp1_165[k]
                   + f_3 * pc_x[k] * qsd_330[k];

        t_551[k] = f_9 * osd_270[k]
                   + f_3 * pc_y[k] * qsd_330[k];

        t_552[k] = f_3 * pc_z[k] * qsd_330[k];

        t_553[k] = f_10 * osd_333[k]
                   + f_3 * pc_x[k] * qsd_333[k];

        t_554[k] = f_3 * pc_z[k] * qsd_331[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, pc_x, pc_y, pc_z, osd_273, osd_275, \
                         osd_335, qsp0_166, qsp1_166, qsd_333, \
                         qsd_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_10 * osd_335[k]
                   + f_3 * pc_x[k] * qsd_335[k];

        t_556[k] = f_9 * osd_273[k]
                   + f_1 * qsp0_166[k]
                   - f_2 * qsp1_166[k]
                   + f_3 * pc_y[k] * qsd_333[k];

        t_557[k] = f_3 * pc_z[k] * qsd_333[k];

        t_558[k] = f_9 * osd_275[k]
                   + f_3 * pc_y[k] * qsd_335[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_z, pc_y, pc_z, osf0_450, osd_270, \
                         osd_276, osf1_450, qsp0_167, qsp1_167, qsd_335, \
                         qsd_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_1 * qsp0_167[k]
                   - f_2 * qsp1_167[k]
                   + f_3 * pc_z[k] * qsd_335[k];

        t_560[k] = pa_z[k] * osf0_450[k]
                   - f_4 * pc_z[k] * osf1_450[k];

        t_561[k] = f_11 * osd_276[k]
                   + f_3 * pc_y[k] * qsd_336[k];

        t_562[k] = f_5 * osd_270[k]
                   + f_3 * pc_z[k] * qsd_336[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pa_z, pc_x, pc_z, osf0_456, osd_339, \
                         osd_340, osd_341, osf1_456, qsd_339, qsd_340, \
                         qsd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_10 * osd_339[k]
                   + f_3 * pc_x[k] * qsd_339[k];

        t_564[k] = f_10 * osd_340[k]
                   + f_3 * pc_x[k] * qsd_340[k];

        t_565[k] = f_10 * osd_341[k]
                   + f_3 * pc_x[k] * qsd_341[k];

        t_566[k] = pa_z[k] * osf0_456[k]
                   - f_4 * pc_z[k] * osf1_456[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pc_y, pc_z, osd_273, osd_275, osd_281, qsp0_170, \
                         qsp1_170, qsd_339, qsd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_5 * osd_273[k]
                   + f_3 * pc_z[k] * qsd_339[k];

        t_568[k] = f_11 * osd_281[k]
                   + f_3 * pc_y[k] * qsd_341[k];

        t_569[k] = f_5 * osd_275[k]
                   + f_1 * qsp0_170[k]
                   - f_2 * qsp1_170[k]
                   + f_3 * pc_z[k] * qsd_341[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pc_x, pc_y, pc_z, osd_276, osd_282, \
                         osd_342, osd_345, qsp0_171, qsp1_171, qsd_342, \
                         qsd_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_10 * osd_342[k]
                   + f_1 * qsp0_171[k]
                   - f_2 * qsp1_171[k]
                   + f_3 * pc_x[k] * qsd_342[k];

        t_571[k] = f_13 * osd_282[k]
                   + f_3 * pc_y[k] * qsd_342[k];

        t_572[k] = f_10 * osd_276[k]
                   + f_3 * pc_z[k] * qsd_342[k];

        t_573[k] = f_10 * osd_345[k]
                   + f_3 * pc_x[k] * qsd_345[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pc_x, pc_y, pc_z, osd_279, osd_285, \
                         osd_346, osd_347, qsp0_172, qsp1_172, qsd_345, qsd_346, \
                         qsd_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_10 * osd_346[k]
                   + f_3 * pc_x[k] * qsd_346[k];

        t_575[k] = f_10 * osd_347[k]
                   + f_3 * pc_x[k] * qsd_347[k];

        t_576[k] = f_13 * osd_285[k]
                   + f_1 * qsp0_172[k]
                   - f_2 * qsp1_172[k]
                   + f_3 * pc_y[k] * qsd_345[k];

        t_577[k] = f_10 * osd_279[k]
                   + f_3 * pc_z[k] * qsd_345[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_x, pc_y, pc_z, osd_281, osd_287, osd_348, \
                         qsp0_173, qsp0_174, qsp1_173, qsp1_174, qsd_347, \
                         qsd_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * osd_287[k]
                   + f_3 * pc_y[k] * qsd_347[k];

        t_579[k] = f_10 * osd_281[k]
                   + f_1 * qsp0_173[k]
                   - f_2 * qsp1_173[k]
                   + f_3 * pc_z[k] * qsd_347[k];

        t_580[k] = f_10 * osd_348[k]
                   + f_1 * qsp0_174[k]
                   - f_2 * qsp1_174[k]
                   + f_3 * pc_x[k] * qsd_348[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, t_584, pc_x, pc_y, pc_z, osd_282, osd_288, \
                         osd_351, osd_352, qsd_348, qsd_351, qsd_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_15 * osd_288[k]
                   + f_3 * pc_y[k] * qsd_348[k];

        t_582[k] = f_12 * osd_282[k]
                   + f_3 * pc_z[k] * qsd_348[k];

        t_583[k] = f_10 * osd_351[k]
                   + f_3 * pc_x[k] * qsd_351[k];

        t_584[k] = f_10 * osd_352[k]
                   + f_3 * pc_x[k] * qsd_352[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pc_x, pc_y, pc_z, osd_285, osd_291, \
                         osd_293, osd_353, qsp0_175, qsp1_175, qsd_351, \
                         qsd_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_10 * osd_353[k]
                   + f_3 * pc_x[k] * qsd_353[k];

        t_586[k] = f_15 * osd_291[k]
                   + f_1 * qsp0_175[k]
                   - f_2 * qsp1_175[k]
                   + f_3 * pc_y[k] * qsd_351[k];

        t_587[k] = f_12 * osd_285[k]
                   + f_3 * pc_z[k] * qsd_351[k];

        t_588[k] = f_15 * osd_293[k]
                   + f_3 * pc_y[k] * qsd_353[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, pc_x, pc_y, pc_z, osd_287, osd_294, osd_354, \
                         qsp0_176, qsp0_177, qsp1_176, qsp1_177, qsd_353, \
                         qsd_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_12 * osd_287[k]
                   + f_1 * qsp0_176[k]
                   - f_2 * qsp1_176[k]
                   + f_3 * pc_z[k] * qsd_353[k];

        t_590[k] = f_10 * osd_354[k]
                   + f_1 * qsp0_177[k]
                   - f_2 * qsp1_177[k]
                   + f_3 * pc_x[k] * qsd_354[k];

        t_591[k] = f_17 * osd_294[k]
                   + f_3 * pc_y[k] * qsd_354[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, osd_288, osd_357, osd_358, \
                         osd_359, qsd_354, qsd_357, qsd_358, qsd_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_14 * osd_288[k]
                   + f_3 * pc_z[k] * qsd_354[k];

        t_593[k] = f_10 * osd_357[k]
                   + f_3 * pc_x[k] * qsd_357[k];

        t_594[k] = f_10 * osd_358[k]
                   + f_3 * pc_x[k] * qsd_358[k];

        t_595[k] = f_10 * osd_359[k]
                   + f_3 * pc_x[k] * qsd_359[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_y, pc_z, osd_291, osd_293, osd_297, \
                         osd_299, qsp0_178, qsp0_179, qsp1_178, qsp1_179, qsd_357, \
                         qsd_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * osd_297[k]
                   + f_1 * qsp0_178[k]
                   - f_2 * qsp1_178[k]
                   + f_3 * pc_y[k] * qsd_357[k];

        t_597[k] = f_14 * osd_291[k]
                   + f_3 * pc_z[k] * qsd_357[k];

        t_598[k] = f_17 * osd_299[k]
                   + f_3 * pc_y[k] * qsd_359[k];

        t_599[k] = f_14 * osd_293[k]
                   + f_1 * qsp0_179[k]
                   - f_2 * qsp1_179[k]
                   + f_3 * pc_z[k] * qsd_359[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, osd_294, osd_300, \
                         osd_360, osd_363, qsp0_180, qsp1_180, qsd_360, \
                         qsd_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_10 * osd_360[k]
                   + f_1 * qsp0_180[k]
                   - f_2 * qsp1_180[k]
                   + f_3 * pc_x[k] * qsd_360[k];

        t_601[k] = f_16 * osd_300[k]
                   + f_3 * pc_y[k] * qsd_360[k];

        t_602[k] = f_16 * osd_294[k]
                   + f_3 * pc_z[k] * qsd_360[k];

        t_603[k] = f_10 * osd_363[k]
                   + f_3 * pc_x[k] * qsd_363[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pc_x, pc_y, pc_z, osd_297, osd_303, \
                         osd_364, osd_365, qsp0_181, qsp1_181, qsd_363, qsd_364, \
                         qsd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_10 * osd_364[k]
                   + f_3 * pc_x[k] * qsd_364[k];

        t_605[k] = f_10 * osd_365[k]
                   + f_3 * pc_x[k] * qsd_365[k];

        t_606[k] = f_16 * osd_303[k]
                   + f_1 * qsp0_181[k]
                   - f_2 * qsp1_181[k]
                   + f_3 * pc_y[k] * qsd_363[k];

        t_607[k] = f_16 * osd_297[k]
                   + f_3 * pc_z[k] * qsd_363[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, pc_x, pc_y, pc_z, osd_299, osd_305, osd_366, \
                         qsp0_182, qsp0_183, qsp1_182, qsp1_183, qsd_365, \
                         qsd_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_16 * osd_305[k]
                   + f_3 * pc_y[k] * qsd_365[k];

        t_609[k] = f_16 * osd_299[k]
                   + f_1 * qsp0_182[k]
                   - f_2 * qsp1_182[k]
                   + f_3 * pc_z[k] * qsd_365[k];

        t_610[k] = f_10 * osd_366[k]
                   + f_1 * qsp0_183[k]
                   - f_2 * qsp1_183[k]
                   + f_3 * pc_x[k] * qsd_366[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pc_x, pc_y, pc_z, osd_300, osd_306, \
                         osd_369, osd_370, qsd_366, qsd_369, qsd_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_14 * osd_306[k]
                   + f_3 * pc_y[k] * qsd_366[k];

        t_612[k] = f_17 * osd_300[k]
                   + f_3 * pc_z[k] * qsd_366[k];

        t_613[k] = f_10 * osd_369[k]
                   + f_3 * pc_x[k] * qsd_369[k];

        t_614[k] = f_10 * osd_370[k]
                   + f_3 * pc_x[k] * qsd_370[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pc_x, pc_y, pc_z, osd_303, osd_309, \
                         osd_311, osd_371, qsp0_184, qsp1_184, qsd_369, \
                         qsd_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_10 * osd_371[k]
                   + f_3 * pc_x[k] * qsd_371[k];

        t_616[k] = f_14 * osd_309[k]
                   + f_1 * qsp0_184[k]
                   - f_2 * qsp1_184[k]
                   + f_3 * pc_y[k] * qsd_369[k];

        t_617[k] = f_17 * osd_303[k]
                   + f_3 * pc_z[k] * qsd_369[k];

        t_618[k] = f_14 * osd_311[k]
                   + f_3 * pc_y[k] * qsd_371[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, pc_x, pc_y, pc_z, osd_305, osd_312, osd_372, \
                         qsp0_185, qsp0_186, qsp1_185, qsp1_186, qsd_371, \
                         qsd_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_17 * osd_305[k]
                   + f_1 * qsp0_185[k]
                   - f_2 * qsp1_185[k]
                   + f_3 * pc_z[k] * qsd_371[k];

        t_620[k] = f_10 * osd_372[k]
                   + f_1 * qsp0_186[k]
                   - f_2 * qsp1_186[k]
                   + f_3 * pc_x[k] * qsd_372[k];

        t_621[k] = f_12 * osd_312[k]
                   + f_3 * pc_y[k] * qsd_372[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_540 = buffer.data(osf0 + 540);
    const auto *osf0_549 = buffer.data(osf0 + 549);
    const auto *osf0_550 = buffer.data(osf0 + 550);
    const auto *osf0_660 = buffer.data(osf0 + 660);
    const auto *osf0_666 = buffer.data(osf0 + 666);
    const auto *osf0_669 = buffer.data(osf0 + 669);
    const auto *osf0_676 = buffer.data(osf0 + 676);
    const auto *osf0_679 = buffer.data(osf0 + 679);
    const auto *osf0_680 = buffer.data(osf0 + 680);
    const auto *osf0_686 = buffer.data(osf0 + 686);
    const auto *osf0_689 = buffer.data(osf0 + 689);
    const auto *osf0_690 = buffer.data(osf0 + 690);
    const auto *osf0_696 = buffer.data(osf0 + 696);
    const auto *osf0_699 = buffer.data(osf0 + 699);
    const auto *osf0_700 = buffer.data(osf0 + 700);
    const auto *osf0_706 = buffer.data(osf0 + 706);
    const auto *osf0_709 = buffer.data(osf0 + 709);
    const auto *osf0_710 = buffer.data(osf0 + 710);
    const auto *osf0_716 = buffer.data(osf0 + 716);
    const auto *osf0_719 = buffer.data(osf0 + 719);
    const auto *osf0_720 = buffer.data(osf0 + 720);
    const auto *osf0_726 = buffer.data(osf0 + 726);
    const auto *osf0_729 = buffer.data(osf0 + 729);
    const auto *osf0_730 = buffer.data(osf0 + 730);
    const auto *osf0_736 = buffer.data(osf0 + 736);
    const auto *osf0_739 = buffer.data(osf0 + 739);
    const auto *osf0_740 = buffer.data(osf0 + 740);
    const auto *osf0_746 = buffer.data(osf0 + 746);
    const auto *osf0_749 = buffer.data(osf0 + 749);

    const auto *osd_306 = buffer.data(osd + 306);
    const auto *osd_309 = buffer.data(osd + 309);
    const auto *osd_311 = buffer.data(osd + 311);
    const auto *osd_312 = buffer.data(osd + 312);
    const auto *osd_315 = buffer.data(osd + 315);
    const auto *osd_317 = buffer.data(osd + 317);
    const auto *osd_318 = buffer.data(osd + 318);
    const auto *osd_321 = buffer.data(osd + 321);
    const auto *osd_323 = buffer.data(osd + 323);
    const auto *osd_324 = buffer.data(osd + 324);
    const auto *osd_327 = buffer.data(osd + 327);
    const auto *osd_329 = buffer.data(osd + 329);
    const auto *osd_330 = buffer.data(osd + 330);
    const auto *osd_333 = buffer.data(osd + 333);
    const auto *osd_335 = buffer.data(osd + 335);
    const auto *osd_336 = buffer.data(osd + 336);
    const auto *osd_339 = buffer.data(osd + 339);
    const auto *osd_341 = buffer.data(osd + 341);
    const auto *osd_342 = buffer.data(osd + 342);
    const auto *osd_345 = buffer.data(osd + 345);
    const auto *osd_347 = buffer.data(osd + 347);
    const auto *osd_348 = buffer.data(osd + 348);
    const auto *osd_351 = buffer.data(osd + 351);
    const auto *osd_353 = buffer.data(osd + 353);
    const auto *osd_354 = buffer.data(osd + 354);
    const auto *osd_357 = buffer.data(osd + 357);
    const auto *osd_359 = buffer.data(osd + 359);
    const auto *osd_360 = buffer.data(osd + 360);
    const auto *osd_363 = buffer.data(osd + 363);
    const auto *osd_365 = buffer.data(osd + 365);
    const auto *osd_366 = buffer.data(osd + 366);
    const auto *osd_369 = buffer.data(osd + 369);
    const auto *osd_371 = buffer.data(osd + 371);
    const auto *osd_372 = buffer.data(osd + 372);
    const auto *osd_375 = buffer.data(osd + 375);
    const auto *osd_376 = buffer.data(osd + 376);
    const auto *osd_377 = buffer.data(osd + 377);
    const auto *osd_378 = buffer.data(osd + 378);
    const auto *osd_381 = buffer.data(osd + 381);
    const auto *osd_382 = buffer.data(osd + 382);
    const auto *osd_383 = buffer.data(osd + 383);
    const auto *osd_387 = buffer.data(osd + 387);
    const auto *osd_388 = buffer.data(osd + 388);
    const auto *osd_389 = buffer.data(osd + 389);
    const auto *osd_390 = buffer.data(osd + 390);
    const auto *osd_393 = buffer.data(osd + 393);
    const auto *osd_395 = buffer.data(osd + 395);
    const auto *osd_396 = buffer.data(osd + 396);
    const auto *osd_399 = buffer.data(osd + 399);
    const auto *osd_401 = buffer.data(osd + 401);
    const auto *osd_405 = buffer.data(osd + 405);
    const auto *osd_406 = buffer.data(osd + 406);
    const auto *osd_407 = buffer.data(osd + 407);
    const auto *osd_408 = buffer.data(osd + 408);
    const auto *osd_411 = buffer.data(osd + 411);
    const auto *osd_412 = buffer.data(osd + 412);
    const auto *osd_413 = buffer.data(osd + 413);
    const auto *osd_414 = buffer.data(osd + 414);
    const auto *osd_417 = buffer.data(osd + 417);
    const auto *osd_418 = buffer.data(osd + 418);
    const auto *osd_419 = buffer.data(osd + 419);
    const auto *osd_420 = buffer.data(osd + 420);
    const auto *osd_423 = buffer.data(osd + 423);
    const auto *osd_424 = buffer.data(osd + 424);
    const auto *osd_425 = buffer.data(osd + 425);
    const auto *osd_426 = buffer.data(osd + 426);
    const auto *osd_429 = buffer.data(osd + 429);
    const auto *osd_430 = buffer.data(osd + 430);
    const auto *osd_431 = buffer.data(osd + 431);
    const auto *osd_432 = buffer.data(osd + 432);
    const auto *osd_435 = buffer.data(osd + 435);
    const auto *osd_436 = buffer.data(osd + 436);
    const auto *osd_437 = buffer.data(osd + 437);
    const auto *osd_438 = buffer.data(osd + 438);
    const auto *osd_441 = buffer.data(osd + 441);
    const auto *osd_442 = buffer.data(osd + 442);
    const auto *osd_443 = buffer.data(osd + 443);
    const auto *osd_444 = buffer.data(osd + 444);
    const auto *osd_447 = buffer.data(osd + 447);
    const auto *osd_448 = buffer.data(osd + 448);
    const auto *osd_449 = buffer.data(osd + 449);

    const auto *osf1_540 = buffer.data(osf1 + 540);
    const auto *osf1_549 = buffer.data(osf1 + 549);
    const auto *osf1_550 = buffer.data(osf1 + 550);
    const auto *osf1_660 = buffer.data(osf1 + 660);
    const auto *osf1_666 = buffer.data(osf1 + 666);
    const auto *osf1_669 = buffer.data(osf1 + 669);
    const auto *osf1_676 = buffer.data(osf1 + 676);
    const auto *osf1_679 = buffer.data(osf1 + 679);
    const auto *osf1_680 = buffer.data(osf1 + 680);
    const auto *osf1_686 = buffer.data(osf1 + 686);
    const auto *osf1_689 = buffer.data(osf1 + 689);
    const auto *osf1_690 = buffer.data(osf1 + 690);
    const auto *osf1_696 = buffer.data(osf1 + 696);
    const auto *osf1_699 = buffer.data(osf1 + 699);
    const auto *osf1_700 = buffer.data(osf1 + 700);
    const auto *osf1_706 = buffer.data(osf1 + 706);
    const auto *osf1_709 = buffer.data(osf1 + 709);
    const auto *osf1_710 = buffer.data(osf1 + 710);
    const auto *osf1_716 = buffer.data(osf1 + 716);
    const auto *osf1_719 = buffer.data(osf1 + 719);
    const auto *osf1_720 = buffer.data(osf1 + 720);
    const auto *osf1_726 = buffer.data(osf1 + 726);
    const auto *osf1_729 = buffer.data(osf1 + 729);
    const auto *osf1_730 = buffer.data(osf1 + 730);
    const auto *osf1_736 = buffer.data(osf1 + 736);
    const auto *osf1_739 = buffer.data(osf1 + 739);
    const auto *osf1_740 = buffer.data(osf1 + 740);
    const auto *osf1_746 = buffer.data(osf1 + 746);
    const auto *osf1_749 = buffer.data(osf1 + 749);

    const auto *qsp0_187 = buffer.data(qsp0 + 187);
    const auto *qsp0_188 = buffer.data(qsp0 + 188);
    const auto *qsp0_189 = buffer.data(qsp0 + 189);
    const auto *qsp0_190 = buffer.data(qsp0 + 190);
    const auto *qsp0_191 = buffer.data(qsp0 + 191);
    const auto *qsp0_193 = buffer.data(qsp0 + 193);
    const auto *qsp0_195 = buffer.data(qsp0 + 195);
    const auto *qsp0_196 = buffer.data(qsp0 + 196);
    const auto *qsp0_197 = buffer.data(qsp0 + 197);

    const auto *qsp1_187 = buffer.data(qsp1 + 187);
    const auto *qsp1_188 = buffer.data(qsp1 + 188);
    const auto *qsp1_189 = buffer.data(qsp1 + 189);
    const auto *qsp1_190 = buffer.data(qsp1 + 190);
    const auto *qsp1_191 = buffer.data(qsp1 + 191);
    const auto *qsp1_193 = buffer.data(qsp1 + 193);
    const auto *qsp1_195 = buffer.data(qsp1 + 195);
    const auto *qsp1_196 = buffer.data(qsp1 + 196);
    const auto *qsp1_197 = buffer.data(qsp1 + 197);

    const auto *qsd_372 = buffer.data(qsd + 372);
    const auto *qsd_375 = buffer.data(qsd + 375);
    const auto *qsd_376 = buffer.data(qsd + 376);
    const auto *qsd_377 = buffer.data(qsd + 377);
    const auto *qsd_378 = buffer.data(qsd + 378);
    const auto *qsd_381 = buffer.data(qsd + 381);
    const auto *qsd_382 = buffer.data(qsd + 382);
    const auto *qsd_383 = buffer.data(qsd + 383);
    const auto *qsd_384 = buffer.data(qsd + 384);
    const auto *qsd_387 = buffer.data(qsd + 387);
    const auto *qsd_388 = buffer.data(qsd + 388);
    const auto *qsd_389 = buffer.data(qsd + 389);
    const auto *qsd_390 = buffer.data(qsd + 390);
    const auto *qsd_392 = buffer.data(qsd + 392);
    const auto *qsd_393 = buffer.data(qsd + 393);
    const auto *qsd_394 = buffer.data(qsd + 394);
    const auto *qsd_395 = buffer.data(qsd + 395);
    const auto *qsd_396 = buffer.data(qsd + 396);
    const auto *qsd_397 = buffer.data(qsd + 397);
    const auto *qsd_399 = buffer.data(qsd + 399);
    const auto *qsd_401 = buffer.data(qsd + 401);
    const auto *qsd_402 = buffer.data(qsd + 402);
    const auto *qsd_405 = buffer.data(qsd + 405);
    const auto *qsd_406 = buffer.data(qsd + 406);
    const auto *qsd_407 = buffer.data(qsd + 407);
    const auto *qsd_408 = buffer.data(qsd + 408);
    const auto *qsd_411 = buffer.data(qsd + 411);
    const auto *qsd_412 = buffer.data(qsd + 412);
    const auto *qsd_413 = buffer.data(qsd + 413);
    const auto *qsd_414 = buffer.data(qsd + 414);
    const auto *qsd_417 = buffer.data(qsd + 417);
    const auto *qsd_418 = buffer.data(qsd + 418);
    const auto *qsd_419 = buffer.data(qsd + 419);
    const auto *qsd_420 = buffer.data(qsd + 420);
    const auto *qsd_423 = buffer.data(qsd + 423);
    const auto *qsd_424 = buffer.data(qsd + 424);
    const auto *qsd_425 = buffer.data(qsd + 425);
    const auto *qsd_426 = buffer.data(qsd + 426);
    const auto *qsd_429 = buffer.data(qsd + 429);
    const auto *qsd_430 = buffer.data(qsd + 430);
    const auto *qsd_431 = buffer.data(qsd + 431);
    const auto *qsd_432 = buffer.data(qsd + 432);
    const auto *qsd_435 = buffer.data(qsd + 435);
    const auto *qsd_436 = buffer.data(qsd + 436);
    const auto *qsd_437 = buffer.data(qsd + 437);
    const auto *qsd_438 = buffer.data(qsd + 438);
    const auto *qsd_441 = buffer.data(qsd + 441);
    const auto *qsd_442 = buffer.data(qsd + 442);
    const auto *qsd_443 = buffer.data(qsd + 443);
    const auto *qsd_444 = buffer.data(qsd + 444);
    const auto *qsd_447 = buffer.data(qsd + 447);
    const auto *qsd_448 = buffer.data(qsd + 448);
    const auto *qsd_449 = buffer.data(qsd + 449);

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pc_x, pc_z, osd_306, osd_375, osd_376, \
                         osd_377, qsd_372, qsd_375, qsd_376, qsd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_15 * osd_306[k]
                   + f_3 * pc_z[k] * qsd_372[k];

        t_623[k] = f_10 * osd_375[k]
                   + f_3 * pc_x[k] * qsd_375[k];

        t_624[k] = f_10 * osd_376[k]
                   + f_3 * pc_x[k] * qsd_376[k];

        t_625[k] = f_10 * osd_377[k]
                   + f_3 * pc_x[k] * qsd_377[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pc_y, pc_z, osd_309, osd_311, osd_315, \
                         osd_317, qsp0_187, qsp0_188, qsp1_187, qsp1_188, qsd_375, \
                         qsd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_12 * osd_315[k]
                   + f_1 * qsp0_187[k]
                   - f_2 * qsp1_187[k]
                   + f_3 * pc_y[k] * qsd_375[k];

        t_627[k] = f_15 * osd_309[k]
                   + f_3 * pc_z[k] * qsd_375[k];

        t_628[k] = f_12 * osd_317[k]
                   + f_3 * pc_y[k] * qsd_377[k];

        t_629[k] = f_15 * osd_311[k]
                   + f_1 * qsp0_188[k]
                   - f_2 * qsp1_188[k]
                   + f_3 * pc_z[k] * qsd_377[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, pc_y, pc_z, osd_312, osd_318, \
                         osd_378, osd_381, qsp0_189, qsp1_189, qsd_378, \
                         qsd_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_10 * osd_378[k]
                   + f_1 * qsp0_189[k]
                   - f_2 * qsp1_189[k]
                   + f_3 * pc_x[k] * qsd_378[k];

        t_631[k] = f_10 * osd_318[k]
                   + f_3 * pc_y[k] * qsd_378[k];

        t_632[k] = f_13 * osd_312[k]
                   + f_3 * pc_z[k] * qsd_378[k];

        t_633[k] = f_10 * osd_381[k]
                   + f_3 * pc_x[k] * qsd_381[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pc_x, pc_y, pc_z, osd_315, osd_321, \
                         osd_382, osd_383, qsp0_190, qsp1_190, qsd_381, qsd_382, \
                         qsd_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_10 * osd_382[k]
                   + f_3 * pc_x[k] * qsd_382[k];

        t_635[k] = f_10 * osd_383[k]
                   + f_3 * pc_x[k] * qsd_383[k];

        t_636[k] = f_10 * osd_321[k]
                   + f_1 * qsp0_190[k]
                   - f_2 * qsp1_190[k]
                   + f_3 * pc_y[k] * qsd_381[k];

        t_637[k] = f_13 * osd_315[k]
                   + f_3 * pc_z[k] * qsd_381[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pa_y, pc_y, pc_z, osf0_540, osd_317, \
                         osd_323, osd_324, osf1_540, qsp0_191, qsp1_191, qsd_383, \
                         qsd_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_10 * osd_323[k]
                   + f_3 * pc_y[k] * qsd_383[k];

        t_639[k] = f_13 * osd_317[k]
                   + f_1 * qsp0_191[k]
                   - f_2 * qsp1_191[k]
                   + f_3 * pc_z[k] * qsd_383[k];

        t_640[k] = pa_y[k] * osf0_540[k]
                   - f_4 * pc_y[k] * osf1_540[k];

        t_641[k] = f_5 * osd_324[k]
                   + f_3 * pc_y[k] * qsd_384[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_z, osd_318, osd_387, osd_388, \
                         osd_389, qsd_384, qsd_387, qsd_388, qsd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_11 * osd_318[k]
                   + f_3 * pc_z[k] * qsd_384[k];

        t_643[k] = f_10 * osd_387[k]
                   + f_3 * pc_x[k] * qsd_387[k];

        t_644[k] = f_10 * osd_388[k]
                   + f_3 * pc_x[k] * qsd_388[k];

        t_645[k] = f_10 * osd_389[k]
                   + f_3 * pc_x[k] * qsd_389[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pa_y, pc_y, pc_z, osf0_549, osd_321, \
                         osd_327, osd_329, osf1_549, qsp0_193, qsp1_193, qsd_387, \
                         qsd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_5 * osd_327[k]
                   + f_1 * qsp0_193[k]
                   - f_2 * qsp1_193[k]
                   + f_3 * pc_y[k] * qsd_387[k];

        t_647[k] = f_11 * osd_321[k]
                   + f_3 * pc_z[k] * qsd_387[k];

        t_648[k] = f_5 * osd_329[k]
                   + f_3 * pc_y[k] * qsd_389[k];

        t_649[k] = pa_y[k] * osf0_549[k]
                   - f_4 * pc_y[k] * osf1_549[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, pc_x, pc_y, pc_z, osd_324, \
                         osd_390, osd_393, qsp0_195, qsp1_195, qsd_390, qsd_392, \
                         qsd_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_10 * osd_390[k]
                   + f_1 * qsp0_195[k]
                   - f_2 * qsp1_195[k]
                   + f_3 * pc_x[k] * qsd_390[k];

        t_651[k] = f_3 * pc_y[k] * qsd_390[k];

        t_652[k] = f_9 * osd_324[k]
                   + f_3 * pc_z[k] * qsd_390[k];

        t_653[k] = f_10 * osd_393[k]
                   + f_3 * pc_x[k] * qsd_393[k];

        t_654[k] = f_3 * pc_y[k] * qsd_392[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, pc_x, pc_y, osd_395, qsp0_196, qsp0_197, \
                         qsp1_196, qsp1_197, qsd_393, qsd_394, \
                         qsd_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_10 * osd_395[k]
                   + f_3 * pc_x[k] * qsd_395[k];

        t_656[k] = f_1 * qsp0_196[k]
                   - f_2 * qsp1_196[k]
                   + f_3 * pc_y[k] * qsd_393[k];

        t_657[k] = f_7 * qsp0_197[k]
                   - f_8 * qsp1_197[k]
                   + f_3 * pc_y[k] * qsd_394[k];

        t_658[k] = f_3 * pc_y[k] * qsd_395[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pa_x, pc_x, pc_y, pc_z, osf0_660, osd_329, \
                         osd_330, osd_396, osf1_660, qsp0_197, qsp1_197, qsd_395, \
                         qsd_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_9 * osd_329[k]
                   + f_1 * qsp0_197[k]
                   - f_2 * qsp1_197[k]
                   + f_3 * pc_z[k] * qsd_395[k];

        t_660[k] = pa_x[k] * osf0_660[k]
                   + f_12 * osd_396[k]
                   - f_4 * pc_x[k] * osf1_660[k];

        t_661[k] = f_6 * osd_330[k]
                   + f_3 * pc_y[k] * qsd_396[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, pa_x, pc_x, pc_z, osf0_666, \
                         osd_399, osd_401, osf1_666, qsd_396, qsd_397, qsd_399, \
                         qsd_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_3 * pc_z[k] * qsd_396[k];

        t_663[k] = f_5 * osd_399[k]
                   + f_3 * pc_x[k] * qsd_399[k];

        t_664[k] = f_3 * pc_z[k] * qsd_397[k];

        t_665[k] = f_5 * osd_401[k]
                   + f_3 * pc_x[k] * qsd_401[k];

        t_666[k] = pa_x[k] * osf0_666[k]
                   - f_4 * pc_x[k] * osf1_666[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, pa_x, pa_z, pc_x, pc_y, pc_z, osf0_550, \
                         osf0_669, osd_335, osf1_550, osf1_669, qsd_399, \
                         qsd_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_3 * pc_z[k] * qsd_399[k];

        t_668[k] = f_6 * osd_335[k]
                   + f_3 * pc_y[k] * qsd_401[k];

        t_669[k] = pa_x[k] * osf0_669[k]
                   - f_4 * pc_x[k] * osf1_669[k];

        t_670[k] = pa_z[k] * osf0_550[k]
                   - f_4 * pc_z[k] * osf1_550[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, pc_x, pc_y, pc_z, osd_330, osd_336, \
                         osd_405, osd_406, qsd_402, qsd_405, qsd_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = f_9 * osd_336[k]
                   + f_3 * pc_y[k] * qsd_402[k];

        t_672[k] = f_5 * osd_330[k]
                   + f_3 * pc_z[k] * qsd_402[k];

        t_673[k] = f_5 * osd_405[k]
                   + f_3 * pc_x[k] * qsd_405[k];

        t_674[k] = f_5 * osd_406[k]
                   + f_3 * pc_x[k] * qsd_406[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, pa_x, pc_x, pc_y, pc_z, osf0_676, \
                         osd_333, osd_341, osd_407, osf1_676, qsd_405, \
                         qsd_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = f_5 * osd_407[k]
                   + f_3 * pc_x[k] * qsd_407[k];

        t_676[k] = pa_x[k] * osf0_676[k]
                   - f_4 * pc_x[k] * osf1_676[k];

        t_677[k] = f_5 * osd_333[k]
                   + f_3 * pc_z[k] * qsd_405[k];

        t_678[k] = f_9 * osd_341[k]
                   + f_3 * pc_y[k] * qsd_407[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, pa_x, pc_x, pc_y, pc_z, osf0_679, \
                         osf0_680, osd_336, osd_342, osd_408, osf1_679, osf1_680, \
                         qsd_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = pa_x[k] * osf0_679[k]
                   - f_4 * pc_x[k] * osf1_679[k];

        t_680[k] = pa_x[k] * osf0_680[k]
                   + f_12 * osd_408[k]
                   - f_4 * pc_x[k] * osf1_680[k];

        t_681[k] = f_11 * osd_342[k]
                   + f_3 * pc_y[k] * qsd_408[k];

        t_682[k] = f_10 * osd_336[k]
                   + f_3 * pc_z[k] * qsd_408[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, pa_x, pc_x, osf0_686, osd_411, osd_412, \
                         osd_413, osf1_686, qsd_411, qsd_412, qsd_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_5 * osd_411[k]
                   + f_3 * pc_x[k] * qsd_411[k];

        t_684[k] = f_5 * osd_412[k]
                   + f_3 * pc_x[k] * qsd_412[k];

        t_685[k] = f_5 * osd_413[k]
                   + f_3 * pc_x[k] * qsd_413[k];

        t_686[k] = pa_x[k] * osf0_686[k]
                   - f_4 * pc_x[k] * osf1_686[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pa_x, pc_x, pc_y, pc_z, osf0_689, osd_339, \
                         osd_347, osf1_689, qsd_411, qsd_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_10 * osd_339[k]
                   + f_3 * pc_z[k] * qsd_411[k];

        t_688[k] = f_11 * osd_347[k]
                   + f_3 * pc_y[k] * qsd_413[k];

        t_689[k] = pa_x[k] * osf0_689[k]
                   - f_4 * pc_x[k] * osf1_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_x, pc_x, pc_y, pc_z, osf0_690, \
                         osd_342, osd_348, osd_414, osd_417, osf1_690, qsd_414, \
                         qsd_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_x[k] * osf0_690[k]
                   + f_12 * osd_414[k]
                   - f_4 * pc_x[k] * osf1_690[k];

        t_691[k] = f_13 * osd_348[k]
                   + f_3 * pc_y[k] * qsd_414[k];

        t_692[k] = f_12 * osd_342[k]
                   + f_3 * pc_z[k] * qsd_414[k];

        t_693[k] = f_5 * osd_417[k]
                   + f_3 * pc_x[k] * qsd_417[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pa_x, pc_x, pc_z, osf0_696, osd_345, \
                         osd_418, osd_419, osf1_696, qsd_417, qsd_418, \
                         qsd_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_5 * osd_418[k]
                   + f_3 * pc_x[k] * qsd_418[k];

        t_695[k] = f_5 * osd_419[k]
                   + f_3 * pc_x[k] * qsd_419[k];

        t_696[k] = pa_x[k] * osf0_696[k]
                   - f_4 * pc_x[k] * osf1_696[k];

        t_697[k] = f_12 * osd_345[k]
                   + f_3 * pc_z[k] * qsd_417[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, t_701, pa_x, pc_x, pc_y, osf0_699, osf0_700, \
                         osd_353, osd_354, osd_420, osf1_699, osf1_700, qsd_419, \
                         qsd_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_13 * osd_353[k]
                   + f_3 * pc_y[k] * qsd_419[k];

        t_699[k] = pa_x[k] * osf0_699[k]
                   - f_4 * pc_x[k] * osf1_699[k];

        t_700[k] = pa_x[k] * osf0_700[k]
                   + f_12 * osd_420[k]
                   - f_4 * pc_x[k] * osf1_700[k];

        t_701[k] = f_15 * osd_354[k]
                   + f_3 * pc_y[k] * qsd_420[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, pc_x, pc_z, osd_348, osd_423, osd_424, \
                         osd_425, qsd_420, qsd_423, qsd_424, qsd_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_14 * osd_348[k]
                   + f_3 * pc_z[k] * qsd_420[k];

        t_703[k] = f_5 * osd_423[k]
                   + f_3 * pc_x[k] * qsd_423[k];

        t_704[k] = f_5 * osd_424[k]
                   + f_3 * pc_x[k] * qsd_424[k];

        t_705[k] = f_5 * osd_425[k]
                   + f_3 * pc_x[k] * qsd_425[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, pa_x, pc_x, pc_y, pc_z, osf0_706, \
                         osf0_709, osd_351, osd_359, osf1_706, osf1_709, qsd_423, \
                         qsd_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = pa_x[k] * osf0_706[k]
                   - f_4 * pc_x[k] * osf1_706[k];

        t_707[k] = f_14 * osd_351[k]
                   + f_3 * pc_z[k] * qsd_423[k];

        t_708[k] = f_15 * osd_359[k]
                   + f_3 * pc_y[k] * qsd_425[k];

        t_709[k] = pa_x[k] * osf0_709[k]
                   - f_4 * pc_x[k] * osf1_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pa_x, pc_x, pc_y, pc_z, osf0_710, \
                         osd_354, osd_360, osd_426, osd_429, osf1_710, qsd_426, \
                         qsd_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = pa_x[k] * osf0_710[k]
                   + f_12 * osd_426[k]
                   - f_4 * pc_x[k] * osf1_710[k];

        t_711[k] = f_17 * osd_360[k]
                   + f_3 * pc_y[k] * qsd_426[k];

        t_712[k] = f_16 * osd_354[k]
                   + f_3 * pc_z[k] * qsd_426[k];

        t_713[k] = f_5 * osd_429[k]
                   + f_3 * pc_x[k] * qsd_429[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pa_x, pc_x, pc_z, osf0_716, osd_357, \
                         osd_430, osd_431, osf1_716, qsd_429, qsd_430, \
                         qsd_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_5 * osd_430[k]
                   + f_3 * pc_x[k] * qsd_430[k];

        t_715[k] = f_5 * osd_431[k]
                   + f_3 * pc_x[k] * qsd_431[k];

        t_716[k] = pa_x[k] * osf0_716[k]
                   - f_4 * pc_x[k] * osf1_716[k];

        t_717[k] = f_16 * osd_357[k]
                   + f_3 * pc_z[k] * qsd_429[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pa_x, pc_x, pc_y, osf0_719, osf0_720, \
                         osd_365, osd_366, osd_432, osf1_719, osf1_720, qsd_431, \
                         qsd_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_17 * osd_365[k]
                   + f_3 * pc_y[k] * qsd_431[k];

        t_719[k] = pa_x[k] * osf0_719[k]
                   - f_4 * pc_x[k] * osf1_719[k];

        t_720[k] = pa_x[k] * osf0_720[k]
                   + f_12 * osd_432[k]
                   - f_4 * pc_x[k] * osf1_720[k];

        t_721[k] = f_16 * osd_366[k]
                   + f_3 * pc_y[k] * qsd_432[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, pc_x, pc_z, osd_360, osd_435, osd_436, \
                         osd_437, qsd_432, qsd_435, qsd_436, qsd_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_17 * osd_360[k]
                   + f_3 * pc_z[k] * qsd_432[k];

        t_723[k] = f_5 * osd_435[k]
                   + f_3 * pc_x[k] * qsd_435[k];

        t_724[k] = f_5 * osd_436[k]
                   + f_3 * pc_x[k] * qsd_436[k];

        t_725[k] = f_5 * osd_437[k]
                   + f_3 * pc_x[k] * qsd_437[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_x, pc_x, pc_y, pc_z, osf0_726, \
                         osf0_729, osd_363, osd_371, osf1_726, osf1_729, qsd_435, \
                         qsd_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = pa_x[k] * osf0_726[k]
                   - f_4 * pc_x[k] * osf1_726[k];

        t_727[k] = f_17 * osd_363[k]
                   + f_3 * pc_z[k] * qsd_435[k];

        t_728[k] = f_16 * osd_371[k]
                   + f_3 * pc_y[k] * qsd_437[k];

        t_729[k] = pa_x[k] * osf0_729[k]
                   - f_4 * pc_x[k] * osf1_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_x, pc_x, pc_y, pc_z, osf0_730, \
                         osd_366, osd_372, osd_438, osd_441, osf1_730, qsd_438, \
                         qsd_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = pa_x[k] * osf0_730[k]
                   + f_12 * osd_438[k]
                   - f_4 * pc_x[k] * osf1_730[k];

        t_731[k] = f_14 * osd_372[k]
                   + f_3 * pc_y[k] * qsd_438[k];

        t_732[k] = f_15 * osd_366[k]
                   + f_3 * pc_z[k] * qsd_438[k];

        t_733[k] = f_5 * osd_441[k]
                   + f_3 * pc_x[k] * qsd_441[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_x, pc_x, pc_z, osf0_736, osd_369, \
                         osd_442, osd_443, osf1_736, qsd_441, qsd_442, \
                         qsd_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_5 * osd_442[k]
                   + f_3 * pc_x[k] * qsd_442[k];

        t_735[k] = f_5 * osd_443[k]
                   + f_3 * pc_x[k] * qsd_443[k];

        t_736[k] = pa_x[k] * osf0_736[k]
                   - f_4 * pc_x[k] * osf1_736[k];

        t_737[k] = f_15 * osd_369[k]
                   + f_3 * pc_z[k] * qsd_441[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, pa_x, pc_x, pc_y, osf0_739, osf0_740, \
                         osd_377, osd_378, osd_444, osf1_739, osf1_740, qsd_443, \
                         qsd_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_14 * osd_377[k]
                   + f_3 * pc_y[k] * qsd_443[k];

        t_739[k] = pa_x[k] * osf0_739[k]
                   - f_4 * pc_x[k] * osf1_739[k];

        t_740[k] = pa_x[k] * osf0_740[k]
                   + f_12 * osd_444[k]
                   - f_4 * pc_x[k] * osf1_740[k];

        t_741[k] = f_12 * osd_378[k]
                   + f_3 * pc_y[k] * qsd_444[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pc_x, pc_z, osd_372, osd_447, osd_448, \
                         osd_449, qsd_444, qsd_447, qsd_448, qsd_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_13 * osd_372[k]
                   + f_3 * pc_z[k] * qsd_444[k];

        t_743[k] = f_5 * osd_447[k]
                   + f_3 * pc_x[k] * qsd_447[k];

        t_744[k] = f_5 * osd_448[k]
                   + f_3 * pc_x[k] * qsd_448[k];

        t_745[k] = f_5 * osd_449[k]
                   + f_3 * pc_x[k] * qsd_449[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pa_x, pc_x, pc_y, pc_z, osf0_746, \
                         osf0_749, osd_375, osd_383, osf1_746, osf1_749, qsd_447, \
                         qsd_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = pa_x[k] * osf0_746[k]
                   - f_4 * pc_x[k] * osf1_746[k];

        t_747[k] = f_13 * osd_375[k]
                   + f_3 * pc_z[k] * qsd_447[k];

        t_748[k] = f_12 * osd_383[k]
                   + f_3 * pc_y[k] * qsd_449[k];

        t_749[k] = pa_x[k] * osf0_749[k]
                   - f_4 * pc_x[k] * osf1_749[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_650 = buffer.data(osf0 + 650);
    const auto *osf0_660 = buffer.data(osf0 + 660);
    const auto *osf0_661 = buffer.data(osf0 + 661);
    const auto *osf0_666 = buffer.data(osf0 + 666);
    const auto *osf0_750 = buffer.data(osf0 + 750);
    const auto *osf0_756 = buffer.data(osf0 + 756);
    const auto *osf0_759 = buffer.data(osf0 + 759);
    const auto *osf0_766 = buffer.data(osf0 + 766);
    const auto *osf0_769 = buffer.data(osf0 + 769);
    const auto *osf0_770 = buffer.data(osf0 + 770);
    const auto *osf0_776 = buffer.data(osf0 + 776);
    const auto *osf0_777 = buffer.data(osf0 + 777);
    const auto *osf0_779 = buffer.data(osf0 + 779);

    const auto *osd_378 = buffer.data(osd + 378);
    const auto *osd_381 = buffer.data(osd + 381);
    const auto *osd_384 = buffer.data(osd + 384);
    const auto *osd_387 = buffer.data(osd + 387);
    const auto *osd_389 = buffer.data(osd + 389);
    const auto *osd_390 = buffer.data(osd + 390);
    const auto *osd_395 = buffer.data(osd + 395);
    const auto *osd_399 = buffer.data(osd + 399);
    const auto *osd_401 = buffer.data(osd + 401);
    const auto *osd_405 = buffer.data(osd + 405);
    const auto *osd_407 = buffer.data(osd + 407);
    const auto *osd_411 = buffer.data(osd + 411);
    const auto *osd_413 = buffer.data(osd + 413);
    const auto *osd_417 = buffer.data(osd + 417);
    const auto *osd_419 = buffer.data(osd + 419);
    const auto *osd_423 = buffer.data(osd + 423);
    const auto *osd_425 = buffer.data(osd + 425);
    const auto *osd_429 = buffer.data(osd + 429);
    const auto *osd_431 = buffer.data(osd + 431);
    const auto *osd_435 = buffer.data(osd + 435);
    const auto *osd_437 = buffer.data(osd + 437);
    const auto *osd_441 = buffer.data(osd + 441);
    const auto *osd_443 = buffer.data(osd + 443);
    const auto *osd_447 = buffer.data(osd + 447);
    const auto *osd_449 = buffer.data(osd + 449);
    const auto *osd_450 = buffer.data(osd + 450);
    const auto *osd_453 = buffer.data(osd + 453);
    const auto *osd_454 = buffer.data(osd + 454);
    const auto *osd_455 = buffer.data(osd + 455);
    const auto *osd_459 = buffer.data(osd + 459);
    const auto *osd_460 = buffer.data(osd + 460);
    const auto *osd_461 = buffer.data(osd + 461);
    const auto *osd_462 = buffer.data(osd + 462);
    const auto *osd_465 = buffer.data(osd + 465);
    const auto *osd_467 = buffer.data(osd + 467);

    const auto *osf1_650 = buffer.data(osf1 + 650);
    const auto *osf1_660 = buffer.data(osf1 + 660);
    const auto *osf1_661 = buffer.data(osf1 + 661);
    const auto *osf1_666 = buffer.data(osf1 + 666);
    const auto *osf1_750 = buffer.data(osf1 + 750);
    const auto *osf1_756 = buffer.data(osf1 + 756);
    const auto *osf1_759 = buffer.data(osf1 + 759);
    const auto *osf1_766 = buffer.data(osf1 + 766);
    const auto *osf1_769 = buffer.data(osf1 + 769);
    const auto *osf1_770 = buffer.data(osf1 + 770);
    const auto *osf1_776 = buffer.data(osf1 + 776);
    const auto *osf1_777 = buffer.data(osf1 + 777);
    const auto *osf1_779 = buffer.data(osf1 + 779);

    const auto *qsp0_234 = buffer.data(qsp0 + 234);
    const auto *qsp0_235 = buffer.data(qsp0 + 235);
    const auto *qsp0_236 = buffer.data(qsp0 + 236);
    const auto *qsp0_239 = buffer.data(qsp0 + 239);
    const auto *qsp0_240 = buffer.data(qsp0 + 240);
    const auto *qsp0_241 = buffer.data(qsp0 + 241);
    const auto *qsp0_242 = buffer.data(qsp0 + 242);
    const auto *qsp0_243 = buffer.data(qsp0 + 243);
    const auto *qsp0_244 = buffer.data(qsp0 + 244);
    const auto *qsp0_245 = buffer.data(qsp0 + 245);
    const auto *qsp0_246 = buffer.data(qsp0 + 246);
    const auto *qsp0_247 = buffer.data(qsp0 + 247);
    const auto *qsp0_248 = buffer.data(qsp0 + 248);
    const auto *qsp0_249 = buffer.data(qsp0 + 249);
    const auto *qsp0_250 = buffer.data(qsp0 + 250);
    const auto *qsp0_251 = buffer.data(qsp0 + 251);
    const auto *qsp0_252 = buffer.data(qsp0 + 252);
    const auto *qsp0_253 = buffer.data(qsp0 + 253);
    const auto *qsp0_254 = buffer.data(qsp0 + 254);
    const auto *qsp0_255 = buffer.data(qsp0 + 255);
    const auto *qsp0_256 = buffer.data(qsp0 + 256);
    const auto *qsp0_257 = buffer.data(qsp0 + 257);
    const auto *qsp0_258 = buffer.data(qsp0 + 258);
    const auto *qsp0_259 = buffer.data(qsp0 + 259);
    const auto *qsp0_260 = buffer.data(qsp0 + 260);
    const auto *qsp0_261 = buffer.data(qsp0 + 261);
    const auto *qsp0_262 = buffer.data(qsp0 + 262);
    const auto *qsp0_263 = buffer.data(qsp0 + 263);

    const auto *qsp1_234 = buffer.data(qsp1 + 234);
    const auto *qsp1_235 = buffer.data(qsp1 + 235);
    const auto *qsp1_236 = buffer.data(qsp1 + 236);
    const auto *qsp1_239 = buffer.data(qsp1 + 239);
    const auto *qsp1_240 = buffer.data(qsp1 + 240);
    const auto *qsp1_241 = buffer.data(qsp1 + 241);
    const auto *qsp1_242 = buffer.data(qsp1 + 242);
    const auto *qsp1_243 = buffer.data(qsp1 + 243);
    const auto *qsp1_244 = buffer.data(qsp1 + 244);
    const auto *qsp1_245 = buffer.data(qsp1 + 245);
    const auto *qsp1_246 = buffer.data(qsp1 + 246);
    const auto *qsp1_247 = buffer.data(qsp1 + 247);
    const auto *qsp1_248 = buffer.data(qsp1 + 248);
    const auto *qsp1_249 = buffer.data(qsp1 + 249);
    const auto *qsp1_250 = buffer.data(qsp1 + 250);
    const auto *qsp1_251 = buffer.data(qsp1 + 251);
    const auto *qsp1_252 = buffer.data(qsp1 + 252);
    const auto *qsp1_253 = buffer.data(qsp1 + 253);
    const auto *qsp1_254 = buffer.data(qsp1 + 254);
    const auto *qsp1_255 = buffer.data(qsp1 + 255);
    const auto *qsp1_256 = buffer.data(qsp1 + 256);
    const auto *qsp1_257 = buffer.data(qsp1 + 257);
    const auto *qsp1_258 = buffer.data(qsp1 + 258);
    const auto *qsp1_259 = buffer.data(qsp1 + 259);
    const auto *qsp1_260 = buffer.data(qsp1 + 260);
    const auto *qsp1_261 = buffer.data(qsp1 + 261);
    const auto *qsp1_262 = buffer.data(qsp1 + 262);
    const auto *qsp1_263 = buffer.data(qsp1 + 263);

    const auto *qsd_450 = buffer.data(qsd + 450);
    const auto *qsd_453 = buffer.data(qsd + 453);
    const auto *qsd_454 = buffer.data(qsd + 454);
    const auto *qsd_455 = buffer.data(qsd + 455);
    const auto *qsd_456 = buffer.data(qsd + 456);
    const auto *qsd_459 = buffer.data(qsd + 459);
    const auto *qsd_460 = buffer.data(qsd + 460);
    const auto *qsd_461 = buffer.data(qsd + 461);
    const auto *qsd_462 = buffer.data(qsd + 462);
    const auto *qsd_464 = buffer.data(qsd + 464);
    const auto *qsd_465 = buffer.data(qsd + 465);
    const auto *qsd_467 = buffer.data(qsd + 467);
    const auto *qsd_468 = buffer.data(qsd + 468);
    const auto *qsd_469 = buffer.data(qsd + 469);
    const auto *qsd_471 = buffer.data(qsd + 471);
    const auto *qsd_472 = buffer.data(qsd + 472);
    const auto *qsd_473 = buffer.data(qsd + 473);
    const auto *qsd_476 = buffer.data(qsd + 476);
    const auto *qsd_477 = buffer.data(qsd + 477);
    const auto *qsd_478 = buffer.data(qsd + 478);
    const auto *qsd_479 = buffer.data(qsd + 479);
    const auto *qsd_480 = buffer.data(qsd + 480);
    const auto *qsd_481 = buffer.data(qsd + 481);
    const auto *qsd_482 = buffer.data(qsd + 482);
    const auto *qsd_483 = buffer.data(qsd + 483);
    const auto *qsd_484 = buffer.data(qsd + 484);
    const auto *qsd_485 = buffer.data(qsd + 485);
    const auto *qsd_486 = buffer.data(qsd + 486);
    const auto *qsd_487 = buffer.data(qsd + 487);
    const auto *qsd_488 = buffer.data(qsd + 488);
    const auto *qsd_489 = buffer.data(qsd + 489);
    const auto *qsd_490 = buffer.data(qsd + 490);
    const auto *qsd_491 = buffer.data(qsd + 491);
    const auto *qsd_492 = buffer.data(qsd + 492);
    const auto *qsd_493 = buffer.data(qsd + 493);
    const auto *qsd_494 = buffer.data(qsd + 494);
    const auto *qsd_495 = buffer.data(qsd + 495);
    const auto *qsd_496 = buffer.data(qsd + 496);
    const auto *qsd_497 = buffer.data(qsd + 497);
    const auto *qsd_498 = buffer.data(qsd + 498);
    const auto *qsd_499 = buffer.data(qsd + 499);
    const auto *qsd_500 = buffer.data(qsd + 500);
    const auto *qsd_501 = buffer.data(qsd + 501);
    const auto *qsd_502 = buffer.data(qsd + 502);
    const auto *qsd_503 = buffer.data(qsd + 503);
    const auto *qsd_504 = buffer.data(qsd + 504);
    const auto *qsd_505 = buffer.data(qsd + 505);
    const auto *qsd_506 = buffer.data(qsd + 506);
    const auto *qsd_507 = buffer.data(qsd + 507);
    const auto *qsd_508 = buffer.data(qsd + 508);
    const auto *qsd_509 = buffer.data(qsd + 509);
    const auto *qsd_510 = buffer.data(qsd + 510);
    const auto *qsd_511 = buffer.data(qsd + 511);
    const auto *qsd_512 = buffer.data(qsd + 512);
    const auto *qsd_513 = buffer.data(qsd + 513);
    const auto *qsd_514 = buffer.data(qsd + 514);
    const auto *qsd_515 = buffer.data(qsd + 515);
    const auto *qsd_516 = buffer.data(qsd + 516);
    const auto *qsd_517 = buffer.data(qsd + 517);
    const auto *qsd_518 = buffer.data(qsd + 518);
    const auto *qsd_519 = buffer.data(qsd + 519);
    const auto *qsd_520 = buffer.data(qsd + 520);
    const auto *qsd_521 = buffer.data(qsd + 521);
    const auto *qsd_522 = buffer.data(qsd + 522);
    const auto *qsd_523 = buffer.data(qsd + 523);
    const auto *qsd_524 = buffer.data(qsd + 524);
    const auto *qsd_525 = buffer.data(qsd + 525);
    const auto *qsd_526 = buffer.data(qsd + 526);
    const auto *qsd_527 = buffer.data(qsd + 527);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pa_x, pc_x, pc_y, pc_z, osf0_750, \
                         osd_378, osd_384, osd_450, osd_453, osf1_750, qsd_450, \
                         qsd_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = pa_x[k] * osf0_750[k]
                   + f_12 * osd_450[k]
                   - f_4 * pc_x[k] * osf1_750[k];

        t_751[k] = f_10 * osd_384[k]
                   + f_3 * pc_y[k] * qsd_450[k];

        t_752[k] = f_11 * osd_378[k]
                   + f_3 * pc_z[k] * qsd_450[k];

        t_753[k] = f_5 * osd_453[k]
                   + f_3 * pc_x[k] * qsd_453[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, pa_x, pc_x, pc_z, osf0_756, osd_381, \
                         osd_454, osd_455, osf1_756, qsd_453, qsd_454, \
                         qsd_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_5 * osd_454[k]
                   + f_3 * pc_x[k] * qsd_454[k];

        t_755[k] = f_5 * osd_455[k]
                   + f_3 * pc_x[k] * qsd_455[k];

        t_756[k] = pa_x[k] * osf0_756[k]
                   - f_4 * pc_x[k] * osf1_756[k];

        t_757[k] = f_11 * osd_381[k]
                   + f_3 * pc_z[k] * qsd_453[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, t_761, pa_x, pa_y, pc_x, pc_y, osf0_650, \
                         osf0_759, osd_389, osd_390, osf1_650, osf1_759, qsd_455, \
                         qsd_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_10 * osd_389[k]
                   + f_3 * pc_y[k] * qsd_455[k];

        t_759[k] = pa_x[k] * osf0_759[k]
                   - f_4 * pc_x[k] * osf1_759[k];

        t_760[k] = pa_y[k] * osf0_650[k]
                   - f_4 * pc_y[k] * osf1_650[k];

        t_761[k] = f_5 * osd_390[k]
                   + f_3 * pc_y[k] * qsd_456[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, pc_x, pc_z, osd_384, osd_459, osd_460, \
                         osd_461, qsd_456, qsd_459, qsd_460, qsd_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_9 * osd_384[k]
                   + f_3 * pc_z[k] * qsd_456[k];

        t_763[k] = f_5 * osd_459[k]
                   + f_3 * pc_x[k] * qsd_459[k];

        t_764[k] = f_5 * osd_460[k]
                   + f_3 * pc_x[k] * qsd_460[k];

        t_765[k] = f_5 * osd_461[k]
                   + f_3 * pc_x[k] * qsd_461[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, pa_x, pc_x, pc_y, pc_z, osf0_766, \
                         osf0_769, osd_387, osd_395, osf1_766, osf1_769, qsd_459, \
                         qsd_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = pa_x[k] * osf0_766[k]
                   - f_4 * pc_x[k] * osf1_766[k];

        t_767[k] = f_9 * osd_387[k]
                   + f_3 * pc_z[k] * qsd_459[k];

        t_768[k] = f_5 * osd_395[k]
                   + f_3 * pc_y[k] * qsd_461[k];

        t_769[k] = pa_x[k] * osf0_769[k]
                   - f_4 * pc_x[k] * osf1_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pa_x, pc_x, pc_y, pc_z, osf0_770, \
                         osd_390, osd_462, osd_465, osf1_770, qsd_462, \
                         qsd_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = pa_x[k] * osf0_770[k]
                   + f_12 * osd_462[k]
                   - f_4 * pc_x[k] * osf1_770[k];

        t_771[k] = f_3 * pc_y[k] * qsd_462[k];

        t_772[k] = f_6 * osd_390[k]
                   + f_3 * pc_z[k] * qsd_462[k];

        t_773[k] = f_5 * osd_465[k]
                   + f_3 * pc_x[k] * qsd_465[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, t_778, pa_x, pc_x, pc_y, osf0_776, \
                         osf0_777, osd_467, osf1_776, osf1_777, qsd_464, \
                         qsd_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_3 * pc_y[k] * qsd_464[k];

        t_775[k] = f_5 * osd_467[k]
                   + f_3 * pc_x[k] * qsd_467[k];

        t_776[k] = pa_x[k] * osf0_776[k]
                   - f_4 * pc_x[k] * osf1_776[k];

        t_777[k] = pa_x[k] * osf0_777[k]
                   - f_4 * pc_x[k] * osf1_777[k];

        t_778[k] = f_3 * pc_y[k] * qsd_467[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, pa_x, pc_x, pc_z, osf0_779, osf1_779, \
                         qsp0_234, qsp0_235, qsp1_234, qsp1_235, qsd_468, \
                         qsd_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = pa_x[k] * osf0_779[k]
                   - f_4 * pc_x[k] * osf1_779[k];

        t_780[k] = f_1 * qsp0_234[k]
                   - f_2 * qsp1_234[k]
                   + f_3 * pc_x[k] * qsd_468[k];

        t_781[k] = f_7 * qsp0_235[k]
                   - f_8 * qsp1_235[k]
                   + f_3 * pc_x[k] * qsd_469[k];

        t_782[k] = f_3 * pc_z[k] * qsd_468[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, t_787, t_788, pc_x, pc_y, pc_z, osd_399, \
                         osd_401, qsp0_235, qsp1_235, qsd_471, qsd_472, \
                         qsd_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_3 * pc_x[k] * qsd_471[k];

        t_784[k] = f_3 * pc_x[k] * qsd_472[k];

        t_785[k] = f_3 * pc_x[k] * qsd_473[k];

        t_786[k] = f_0 * osd_399[k]
                   + f_1 * qsp0_235[k]
                   - f_2 * qsp1_235[k]
                   + f_3 * pc_y[k] * qsd_471[k];

        t_787[k] = f_3 * pc_z[k] * qsd_471[k];

        t_788[k] = f_0 * osd_401[k]
                   + f_3 * pc_y[k] * qsd_473[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pa_z, pc_z, osf0_660, osf0_661, osf1_660, \
                         osf1_661, qsp0_236, qsp1_236, qsd_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_1 * qsp0_236[k]
                   - f_2 * qsp1_236[k]
                   + f_3 * pc_z[k] * qsd_473[k];

        t_790[k] = pa_z[k] * osf0_660[k]
                   - f_4 * pc_z[k] * osf1_660[k];

        t_791[k] = pa_z[k] * osf0_661[k]
                   - f_4 * pc_z[k] * osf1_661[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, t_796, pa_z, pc_x, pc_z, osf0_666, \
                         osf1_666, qsp0_239, qsp1_239, qsd_476, qsd_477, qsd_478, \
                         qsd_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_7 * qsp0_239[k]
                   - f_8 * qsp1_239[k]
                   + f_3 * pc_x[k] * qsd_476[k];

        t_793[k] = f_3 * pc_x[k] * qsd_477[k];

        t_794[k] = f_3 * pc_x[k] * qsd_478[k];

        t_795[k] = f_3 * pc_x[k] * qsd_479[k];

        t_796[k] = pa_z[k] * osf0_666[k]
                   - f_4 * pc_z[k] * osf1_666[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pc_y, pc_z, osd_399, osd_401, osd_407, qsp0_239, \
                         qsp1_239, qsd_477, qsd_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_5 * osd_399[k]
                   + f_3 * pc_z[k] * qsd_477[k];

        t_798[k] = f_6 * osd_407[k]
                   + f_3 * pc_y[k] * qsd_479[k];

        t_799[k] = f_5 * osd_401[k]
                   + f_1 * qsp0_239[k]
                   - f_2 * qsp1_239[k]
                   + f_3 * pc_z[k] * qsd_479[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, pc_x, qsp0_240, qsp0_241, qsp0_242, \
                         qsp1_240, qsp1_241, qsp1_242, qsd_480, qsd_481, qsd_482, \
                         qsd_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_1 * qsp0_240[k]
                   - f_2 * qsp1_240[k]
                   + f_3 * pc_x[k] * qsd_480[k];

        t_801[k] = f_7 * qsp0_241[k]
                   - f_8 * qsp1_241[k]
                   + f_3 * pc_x[k] * qsd_481[k];

        t_802[k] = f_7 * qsp0_242[k]
                   - f_8 * qsp1_242[k]
                   + f_3 * pc_x[k] * qsd_482[k];

        t_803[k] = f_3 * pc_x[k] * qsd_483[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, t_808, pc_x, pc_y, pc_z, osd_405, \
                         osd_411, osd_413, qsp0_241, qsp1_241, qsd_483, qsd_484, \
                         qsd_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_3 * pc_x[k] * qsd_484[k];

        t_805[k] = f_3 * pc_x[k] * qsd_485[k];

        t_806[k] = f_9 * osd_411[k]
                   + f_1 * qsp0_241[k]
                   - f_2 * qsp1_241[k]
                   + f_3 * pc_y[k] * qsd_483[k];

        t_807[k] = f_10 * osd_405[k]
                   + f_3 * pc_z[k] * qsd_483[k];

        t_808[k] = f_9 * osd_413[k]
                   + f_3 * pc_y[k] * qsd_485[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_z, osd_407, qsp0_242, qsp0_243, \
                         qsp0_244, qsp1_242, qsp1_243, qsp1_244, qsd_485, qsd_486, \
                         qsd_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_10 * osd_407[k]
                   + f_1 * qsp0_242[k]
                   - f_2 * qsp1_242[k]
                   + f_3 * pc_z[k] * qsd_485[k];

        t_810[k] = f_1 * qsp0_243[k]
                   - f_2 * qsp1_243[k]
                   + f_3 * pc_x[k] * qsd_486[k];

        t_811[k] = f_7 * qsp0_244[k]
                   - f_8 * qsp1_244[k]
                   + f_3 * pc_x[k] * qsd_487[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, t_816, pc_x, pc_y, osd_417, qsp0_244, \
                         qsp0_245, qsp1_244, qsp1_245, qsd_488, qsd_489, qsd_490, \
                         qsd_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_7 * qsp0_245[k]
                   - f_8 * qsp1_245[k]
                   + f_3 * pc_x[k] * qsd_488[k];

        t_813[k] = f_3 * pc_x[k] * qsd_489[k];

        t_814[k] = f_3 * pc_x[k] * qsd_490[k];

        t_815[k] = f_3 * pc_x[k] * qsd_491[k];

        t_816[k] = f_11 * osd_417[k]
                   + f_1 * qsp0_244[k]
                   - f_2 * qsp1_244[k]
                   + f_3 * pc_y[k] * qsd_489[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pc_y, pc_z, osd_411, osd_413, osd_419, qsp0_245, \
                         qsp1_245, qsd_489, qsd_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_12 * osd_411[k]
                   + f_3 * pc_z[k] * qsd_489[k];

        t_818[k] = f_11 * osd_419[k]
                   + f_3 * pc_y[k] * qsd_491[k];

        t_819[k] = f_12 * osd_413[k]
                   + f_1 * qsp0_245[k]
                   - f_2 * qsp1_245[k]
                   + f_3 * pc_z[k] * qsd_491[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_x, qsp0_246, qsp0_247, qsp0_248, \
                         qsp1_246, qsp1_247, qsp1_248, qsd_492, qsd_493, qsd_494, \
                         qsd_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_1 * qsp0_246[k]
                   - f_2 * qsp1_246[k]
                   + f_3 * pc_x[k] * qsd_492[k];

        t_821[k] = f_7 * qsp0_247[k]
                   - f_8 * qsp1_247[k]
                   + f_3 * pc_x[k] * qsd_493[k];

        t_822[k] = f_7 * qsp0_248[k]
                   - f_8 * qsp1_248[k]
                   + f_3 * pc_x[k] * qsd_494[k];

        t_823[k] = f_3 * pc_x[k] * qsd_495[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, t_827, t_828, pc_x, pc_y, pc_z, osd_417, \
                         osd_423, osd_425, qsp0_247, qsp1_247, qsd_495, qsd_496, \
                         qsd_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_3 * pc_x[k] * qsd_496[k];

        t_825[k] = f_3 * pc_x[k] * qsd_497[k];

        t_826[k] = f_13 * osd_423[k]
                   + f_1 * qsp0_247[k]
                   - f_2 * qsp1_247[k]
                   + f_3 * pc_y[k] * qsd_495[k];

        t_827[k] = f_14 * osd_417[k]
                   + f_3 * pc_z[k] * qsd_495[k];

        t_828[k] = f_13 * osd_425[k]
                   + f_3 * pc_y[k] * qsd_497[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, pc_x, pc_z, osd_419, qsp0_248, qsp0_249, \
                         qsp0_250, qsp1_248, qsp1_249, qsp1_250, qsd_497, qsd_498, \
                         qsd_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_14 * osd_419[k]
                   + f_1 * qsp0_248[k]
                   - f_2 * qsp1_248[k]
                   + f_3 * pc_z[k] * qsd_497[k];

        t_830[k] = f_1 * qsp0_249[k]
                   - f_2 * qsp1_249[k]
                   + f_3 * pc_x[k] * qsd_498[k];

        t_831[k] = f_7 * qsp0_250[k]
                   - f_8 * qsp1_250[k]
                   + f_3 * pc_x[k] * qsd_499[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, t_836, pc_x, pc_y, osd_429, qsp0_250, \
                         qsp0_251, qsp1_250, qsp1_251, qsd_500, qsd_501, qsd_502, \
                         qsd_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_7 * qsp0_251[k]
                   - f_8 * qsp1_251[k]
                   + f_3 * pc_x[k] * qsd_500[k];

        t_833[k] = f_3 * pc_x[k] * qsd_501[k];

        t_834[k] = f_3 * pc_x[k] * qsd_502[k];

        t_835[k] = f_3 * pc_x[k] * qsd_503[k];

        t_836[k] = f_15 * osd_429[k]
                   + f_1 * qsp0_250[k]
                   - f_2 * qsp1_250[k]
                   + f_3 * pc_y[k] * qsd_501[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, osd_423, osd_425, osd_431, qsp0_251, \
                         qsp1_251, qsd_501, qsd_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_16 * osd_423[k]
                   + f_3 * pc_z[k] * qsd_501[k];

        t_838[k] = f_15 * osd_431[k]
                   + f_3 * pc_y[k] * qsd_503[k];

        t_839[k] = f_16 * osd_425[k]
                   + f_1 * qsp0_251[k]
                   - f_2 * qsp1_251[k]
                   + f_3 * pc_z[k] * qsd_503[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, pc_x, qsp0_252, qsp0_253, qsp0_254, \
                         qsp1_252, qsp1_253, qsp1_254, qsd_504, qsd_505, qsd_506, \
                         qsd_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_1 * qsp0_252[k]
                   - f_2 * qsp1_252[k]
                   + f_3 * pc_x[k] * qsd_504[k];

        t_841[k] = f_7 * qsp0_253[k]
                   - f_8 * qsp1_253[k]
                   + f_3 * pc_x[k] * qsd_505[k];

        t_842[k] = f_7 * qsp0_254[k]
                   - f_8 * qsp1_254[k]
                   + f_3 * pc_x[k] * qsd_506[k];

        t_843[k] = f_3 * pc_x[k] * qsd_507[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, t_847, t_848, pc_x, pc_y, pc_z, osd_429, \
                         osd_435, osd_437, qsp0_253, qsp1_253, qsd_507, qsd_508, \
                         qsd_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = f_3 * pc_x[k] * qsd_508[k];

        t_845[k] = f_3 * pc_x[k] * qsd_509[k];

        t_846[k] = f_17 * osd_435[k]
                   + f_1 * qsp0_253[k]
                   - f_2 * qsp1_253[k]
                   + f_3 * pc_y[k] * qsd_507[k];

        t_847[k] = f_17 * osd_429[k]
                   + f_3 * pc_z[k] * qsd_507[k];

        t_848[k] = f_17 * osd_437[k]
                   + f_3 * pc_y[k] * qsd_509[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, osd_431, qsp0_254, qsp0_255, \
                         qsp0_256, qsp1_254, qsp1_255, qsp1_256, qsd_509, qsd_510, \
                         qsd_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_17 * osd_431[k]
                   + f_1 * qsp0_254[k]
                   - f_2 * qsp1_254[k]
                   + f_3 * pc_z[k] * qsd_509[k];

        t_850[k] = f_1 * qsp0_255[k]
                   - f_2 * qsp1_255[k]
                   + f_3 * pc_x[k] * qsd_510[k];

        t_851[k] = f_7 * qsp0_256[k]
                   - f_8 * qsp1_256[k]
                   + f_3 * pc_x[k] * qsd_511[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, t_856, pc_x, pc_y, osd_441, qsp0_256, \
                         qsp0_257, qsp1_256, qsp1_257, qsd_512, qsd_513, qsd_514, \
                         qsd_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_7 * qsp0_257[k]
                   - f_8 * qsp1_257[k]
                   + f_3 * pc_x[k] * qsd_512[k];

        t_853[k] = f_3 * pc_x[k] * qsd_513[k];

        t_854[k] = f_3 * pc_x[k] * qsd_514[k];

        t_855[k] = f_3 * pc_x[k] * qsd_515[k];

        t_856[k] = f_16 * osd_441[k]
                   + f_1 * qsp0_256[k]
                   - f_2 * qsp1_256[k]
                   + f_3 * pc_y[k] * qsd_513[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pc_y, pc_z, osd_435, osd_437, osd_443, qsp0_257, \
                         qsp1_257, qsd_513, qsd_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = f_15 * osd_435[k]
                   + f_3 * pc_z[k] * qsd_513[k];

        t_858[k] = f_16 * osd_443[k]
                   + f_3 * pc_y[k] * qsd_515[k];

        t_859[k] = f_15 * osd_437[k]
                   + f_1 * qsp0_257[k]
                   - f_2 * qsp1_257[k]
                   + f_3 * pc_z[k] * qsd_515[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, pc_x, qsp0_258, qsp0_259, qsp0_260, \
                         qsp1_258, qsp1_259, qsp1_260, qsd_516, qsd_517, qsd_518, \
                         qsd_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_1 * qsp0_258[k]
                   - f_2 * qsp1_258[k]
                   + f_3 * pc_x[k] * qsd_516[k];

        t_861[k] = f_7 * qsp0_259[k]
                   - f_8 * qsp1_259[k]
                   + f_3 * pc_x[k] * qsd_517[k];

        t_862[k] = f_7 * qsp0_260[k]
                   - f_8 * qsp1_260[k]
                   + f_3 * pc_x[k] * qsd_518[k];

        t_863[k] = f_3 * pc_x[k] * qsd_519[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, osd_441, \
                         osd_447, osd_449, qsp0_259, qsp1_259, qsd_519, qsd_520, \
                         qsd_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_3 * pc_x[k] * qsd_520[k];

        t_865[k] = f_3 * pc_x[k] * qsd_521[k];

        t_866[k] = f_14 * osd_447[k]
                   + f_1 * qsp0_259[k]
                   - f_2 * qsp1_259[k]
                   + f_3 * pc_y[k] * qsd_519[k];

        t_867[k] = f_13 * osd_441[k]
                   + f_3 * pc_z[k] * qsd_519[k];

        t_868[k] = f_14 * osd_449[k]
                   + f_3 * pc_y[k] * qsd_521[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, osd_443, qsp0_260, qsp0_261, \
                         qsp0_262, qsp1_260, qsp1_261, qsp1_262, qsd_521, qsd_522, \
                         qsd_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_13 * osd_443[k]
                   + f_1 * qsp0_260[k]
                   - f_2 * qsp1_260[k]
                   + f_3 * pc_z[k] * qsd_521[k];

        t_870[k] = f_1 * qsp0_261[k]
                   - f_2 * qsp1_261[k]
                   + f_3 * pc_x[k] * qsd_522[k];

        t_871[k] = f_7 * qsp0_262[k]
                   - f_8 * qsp1_262[k]
                   + f_3 * pc_x[k] * qsd_523[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, pc_x, pc_y, osd_453, qsp0_262, \
                         qsp0_263, qsp1_262, qsp1_263, qsd_524, qsd_525, qsd_526, \
                         qsd_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_7 * qsp0_263[k]
                   - f_8 * qsp1_263[k]
                   + f_3 * pc_x[k] * qsd_524[k];

        t_873[k] = f_3 * pc_x[k] * qsd_525[k];

        t_874[k] = f_3 * pc_x[k] * qsd_526[k];

        t_875[k] = f_3 * pc_x[k] * qsd_527[k];

        t_876[k] = f_12 * osd_453[k]
                   + f_1 * qsp0_262[k]
                   - f_2 * qsp1_262[k]
                   + f_3 * pc_y[k] * qsd_525[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, osd_447, osd_449, osd_455, qsp0_263, \
                         qsp1_263, qsd_525, qsd_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_11 * osd_447[k]
                   + f_3 * pc_z[k] * qsd_525[k];

        t_878[k] = f_12 * osd_455[k]
                   + f_3 * pc_y[k] * qsd_527[k];

        t_879[k] = f_11 * osd_449[k]
                   + f_1 * qsp0_263[k]
                   - f_2 * qsp1_263[k]
                   + f_3 * pc_z[k] * qsd_527[k];
    }
}

static auto
compute_prim_qsf_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osf0,
                                                          const size_t osd, const size_t osf1,
                                                          const size_t qsp0, const size_t qsp1,
                                                          const size_t qsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osf0_770 = buffer.data(osf0 + 770);
    const auto *osf0_772 = buffer.data(osf0 + 772);
    const auto *osf0_776 = buffer.data(osf0 + 776);
    const auto *osf0_779 = buffer.data(osf0 + 779);

    const auto *osd_453 = buffer.data(osd + 453);
    const auto *osd_455 = buffer.data(osd + 455);
    const auto *osd_459 = buffer.data(osd + 459);
    const auto *osd_461 = buffer.data(osd + 461);
    const auto *osd_465 = buffer.data(osd + 465);
    const auto *osd_467 = buffer.data(osd + 467);

    const auto *osf1_770 = buffer.data(osf1 + 770);
    const auto *osf1_772 = buffer.data(osf1 + 772);
    const auto *osf1_776 = buffer.data(osf1 + 776);
    const auto *osf1_779 = buffer.data(osf1 + 779);

    const auto *qsp0_264 = buffer.data(qsp0 + 264);
    const auto *qsp0_265 = buffer.data(qsp0 + 265);
    const auto *qsp0_266 = buffer.data(qsp0 + 266);
    const auto *qsp0_268 = buffer.data(qsp0 + 268);
    const auto *qsp0_270 = buffer.data(qsp0 + 270);
    const auto *qsp0_271 = buffer.data(qsp0 + 271);
    const auto *qsp0_272 = buffer.data(qsp0 + 272);

    const auto *qsp1_264 = buffer.data(qsp1 + 264);
    const auto *qsp1_265 = buffer.data(qsp1 + 265);
    const auto *qsp1_266 = buffer.data(qsp1 + 266);
    const auto *qsp1_268 = buffer.data(qsp1 + 268);
    const auto *qsp1_270 = buffer.data(qsp1 + 270);
    const auto *qsp1_271 = buffer.data(qsp1 + 271);
    const auto *qsp1_272 = buffer.data(qsp1 + 272);

    const auto *qsd_528 = buffer.data(qsd + 528);
    const auto *qsd_529 = buffer.data(qsd + 529);
    const auto *qsd_530 = buffer.data(qsd + 530);
    const auto *qsd_531 = buffer.data(qsd + 531);
    const auto *qsd_532 = buffer.data(qsd + 532);
    const auto *qsd_533 = buffer.data(qsd + 533);
    const auto *qsd_535 = buffer.data(qsd + 535);
    const auto *qsd_537 = buffer.data(qsd + 537);
    const auto *qsd_538 = buffer.data(qsd + 538);
    const auto *qsd_539 = buffer.data(qsd + 539);
    const auto *qsd_540 = buffer.data(qsd + 540);
    const auto *qsd_542 = buffer.data(qsd + 542);
    const auto *qsd_543 = buffer.data(qsd + 543);
    const auto *qsd_544 = buffer.data(qsd + 544);
    const auto *qsd_545 = buffer.data(qsd + 545);

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pc_x, qsp0_264, qsp0_265, qsp0_266, \
                         qsp1_264, qsp1_265, qsp1_266, qsd_528, qsd_529, qsd_530, \
                         qsd_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_1 * qsp0_264[k]
                   - f_2 * qsp1_264[k]
                   + f_3 * pc_x[k] * qsd_528[k];

        t_881[k] = f_7 * qsp0_265[k]
                   - f_8 * qsp1_265[k]
                   + f_3 * pc_x[k] * qsd_529[k];

        t_882[k] = f_7 * qsp0_266[k]
                   - f_8 * qsp1_266[k]
                   + f_3 * pc_x[k] * qsd_530[k];

        t_883[k] = f_3 * pc_x[k] * qsd_531[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, pc_y, pc_z, osd_453, \
                         osd_459, osd_461, qsp0_265, qsp1_265, qsd_531, qsd_532, \
                         qsd_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_3 * pc_x[k] * qsd_532[k];

        t_885[k] = f_3 * pc_x[k] * qsd_533[k];

        t_886[k] = f_10 * osd_459[k]
                   + f_1 * qsp0_265[k]
                   - f_2 * qsp1_265[k]
                   + f_3 * pc_y[k] * qsd_531[k];

        t_887[k] = f_9 * osd_453[k]
                   + f_3 * pc_z[k] * qsd_531[k];

        t_888[k] = f_10 * osd_461[k]
                   + f_3 * pc_y[k] * qsd_533[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, pa_y, pc_x, pc_y, pc_z, osf0_770, osd_455, \
                         osf1_770, qsp0_266, qsp0_268, qsp1_266, qsp1_268, qsd_533, \
                         qsd_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_9 * osd_455[k]
                   + f_1 * qsp0_266[k]
                   - f_2 * qsp1_266[k]
                   + f_3 * pc_z[k] * qsd_533[k];

        t_890[k] = pa_y[k] * osf0_770[k]
                   - f_4 * pc_y[k] * osf1_770[k];

        t_891[k] = f_7 * qsp0_268[k]
                   - f_8 * qsp1_268[k]
                   + f_3 * pc_x[k] * qsd_535[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, t_896, pa_y, pc_x, pc_y, osf0_772, \
                         osf0_776, osd_465, osf1_772, osf1_776, qsd_537, qsd_538, \
                         qsd_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = pa_y[k] * osf0_772[k]
                   - f_4 * pc_y[k] * osf1_772[k];

        t_893[k] = f_3 * pc_x[k] * qsd_537[k];

        t_894[k] = f_3 * pc_x[k] * qsd_538[k];

        t_895[k] = f_3 * pc_x[k] * qsd_539[k];

        t_896[k] = pa_y[k] * osf0_776[k]
                   + f_12 * osd_465[k]
                   - f_4 * pc_y[k] * osf1_776[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pa_y, pc_y, pc_z, osf0_779, osd_459, osd_467, \
                         osf1_779, qsd_537, qsd_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_6 * osd_459[k]
                   + f_3 * pc_z[k] * qsd_537[k];

        t_898[k] = f_5 * osd_467[k]
                   + f_3 * pc_y[k] * qsd_539[k];

        t_899[k] = pa_y[k] * osf0_779[k]
                   - f_4 * pc_y[k] * osf1_779[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, pc_x, pc_y, qsp0_270, qsp0_272, \
                         qsp1_270, qsp1_272, qsd_540, qsd_542, qsd_543, \
                         qsd_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_1 * qsp0_270[k]
                   - f_2 * qsp1_270[k]
                   + f_3 * pc_x[k] * qsd_540[k];

        t_901[k] = f_3 * pc_y[k] * qsd_540[k];

        t_902[k] = f_7 * qsp0_272[k]
                   - f_8 * qsp1_272[k]
                   + f_3 * pc_x[k] * qsd_542[k];

        t_903[k] = f_3 * pc_x[k] * qsd_543[k];

        t_904[k] = f_3 * pc_x[k] * qsd_544[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, pc_x, pc_y, pc_z, osd_467, \
                         qsp0_271, qsp0_272, qsp1_271, qsp1_272, qsd_543, qsd_544, \
                         qsd_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_3 * pc_x[k] * qsd_545[k];

        t_906[k] = f_1 * qsp0_271[k]
                   - f_2 * qsp1_271[k]
                   + f_3 * pc_y[k] * qsd_543[k];

        t_907[k] = f_7 * qsp0_272[k]
                   - f_8 * qsp1_272[k]
                   + f_3 * pc_y[k] * qsd_544[k];

        t_908[k] = f_3 * pc_y[k] * qsd_545[k];

        t_909[k] = f_0 * osd_467[k]
                   + f_1 * qsp0_272[k]
                   - f_2 * qsp1_272[k]
                   + f_3 * pc_z[k] * qsd_545[k];
    }
}

auto
compute_prim_qsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t osf0, const size_t osd,
                                                   const size_t osf1, const size_t qsp0,
                                                   const size_t qsp1, const size_t qsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_qsf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);

    compute_prim_qsf_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, osf0, osd,
                                                              osf1, qsp0, qsp1, qsd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
