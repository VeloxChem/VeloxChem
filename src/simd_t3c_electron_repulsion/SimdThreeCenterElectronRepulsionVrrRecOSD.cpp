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


#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_osd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsd0,
                                                          const size_t nsp, const size_t nsd1,
                                                          const size_t oss0, const size_t oss1,
                                                          const size_t osp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
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

    const auto *nsd0_0 = buffer.data(nsd0 + 0);
    const auto *nsd0_3 = buffer.data(nsd0 + 3);
    const auto *nsd0_5 = buffer.data(nsd0 + 5);
    const auto *nsd0_9 = buffer.data(nsd0 + 9);
    const auto *nsd0_12 = buffer.data(nsd0 + 12);
    const auto *nsd0_17 = buffer.data(nsd0 + 17);
    const auto *nsd0_18 = buffer.data(nsd0 + 18);
    const auto *nsd0_21 = buffer.data(nsd0 + 21);
    const auto *nsd0_30 = buffer.data(nsd0 + 30);
    const auto *nsd0_35 = buffer.data(nsd0 + 35);
    const auto *nsd0_36 = buffer.data(nsd0 + 36);
    const auto *nsd0_39 = buffer.data(nsd0 + 39);
    const auto *nsd0_54 = buffer.data(nsd0 + 54);
    const auto *nsd0_59 = buffer.data(nsd0 + 59);
    const auto *nsd0_60 = buffer.data(nsd0 + 60);
    const auto *nsd0_63 = buffer.data(nsd0 + 63);
    const auto *nsd0_84 = buffer.data(nsd0 + 84);
    const auto *nsd0_89 = buffer.data(nsd0 + 89);

    const auto *nsp_0 = buffer.data(nsp + 0);
    const auto *nsp_1 = buffer.data(nsp + 1);
    const auto *nsp_2 = buffer.data(nsp + 2);
    const auto *nsp_4 = buffer.data(nsp + 4);
    const auto *nsp_8 = buffer.data(nsp + 8);
    const auto *nsp_9 = buffer.data(nsp + 9);
    const auto *nsp_10 = buffer.data(nsp + 10);
    const auto *nsp_11 = buffer.data(nsp + 11);
    const auto *nsp_13 = buffer.data(nsp + 13);
    const auto *nsp_14 = buffer.data(nsp + 14);
    const auto *nsp_15 = buffer.data(nsp + 15);
    const auto *nsp_16 = buffer.data(nsp + 16);
    const auto *nsp_17 = buffer.data(nsp + 17);
    const auto *nsp_18 = buffer.data(nsp + 18);
    const auto *nsp_19 = buffer.data(nsp + 19);
    const auto *nsp_20 = buffer.data(nsp + 20);
    const auto *nsp_22 = buffer.data(nsp + 22);
    const auto *nsp_23 = buffer.data(nsp + 23);
    const auto *nsp_25 = buffer.data(nsp + 25);
    const auto *nsp_26 = buffer.data(nsp + 26);
    const auto *nsp_27 = buffer.data(nsp + 27);
    const auto *nsp_28 = buffer.data(nsp + 28);
    const auto *nsp_29 = buffer.data(nsp + 29);
    const auto *nsp_30 = buffer.data(nsp + 30);
    const auto *nsp_31 = buffer.data(nsp + 31);
    const auto *nsp_32 = buffer.data(nsp + 32);
    const auto *nsp_34 = buffer.data(nsp + 34);
    const auto *nsp_35 = buffer.data(nsp + 35);
    const auto *nsp_36 = buffer.data(nsp + 36);
    const auto *nsp_37 = buffer.data(nsp + 37);
    const auto *nsp_38 = buffer.data(nsp + 38);
    const auto *nsp_40 = buffer.data(nsp + 40);
    const auto *nsp_41 = buffer.data(nsp + 41);
    const auto *nsp_42 = buffer.data(nsp + 42);
    const auto *nsp_43 = buffer.data(nsp + 43);
    const auto *nsp_44 = buffer.data(nsp + 44);
    const auto *nsp_45 = buffer.data(nsp + 45);
    const auto *nsp_46 = buffer.data(nsp + 46);
    const auto *nsp_49 = buffer.data(nsp + 49);
    const auto *nsp_50 = buffer.data(nsp + 50);
    const auto *nsp_51 = buffer.data(nsp + 51);
    const auto *nsp_52 = buffer.data(nsp + 52);
    const auto *nsp_53 = buffer.data(nsp + 53);
    const auto *nsp_54 = buffer.data(nsp + 54);
    const auto *nsp_55 = buffer.data(nsp + 55);
    const auto *nsp_56 = buffer.data(nsp + 56);
    const auto *nsp_58 = buffer.data(nsp + 58);
    const auto *nsp_59 = buffer.data(nsp + 59);
    const auto *nsp_60 = buffer.data(nsp + 60);
    const auto *nsp_62 = buffer.data(nsp + 62);

    const auto *nsd1_0 = buffer.data(nsd1 + 0);
    const auto *nsd1_3 = buffer.data(nsd1 + 3);
    const auto *nsd1_5 = buffer.data(nsd1 + 5);
    const auto *nsd1_9 = buffer.data(nsd1 + 9);
    const auto *nsd1_12 = buffer.data(nsd1 + 12);
    const auto *nsd1_17 = buffer.data(nsd1 + 17);
    const auto *nsd1_18 = buffer.data(nsd1 + 18);
    const auto *nsd1_21 = buffer.data(nsd1 + 21);
    const auto *nsd1_30 = buffer.data(nsd1 + 30);
    const auto *nsd1_35 = buffer.data(nsd1 + 35);
    const auto *nsd1_36 = buffer.data(nsd1 + 36);
    const auto *nsd1_39 = buffer.data(nsd1 + 39);
    const auto *nsd1_54 = buffer.data(nsd1 + 54);
    const auto *nsd1_59 = buffer.data(nsd1 + 59);
    const auto *nsd1_60 = buffer.data(nsd1 + 60);
    const auto *nsd1_63 = buffer.data(nsd1 + 63);
    const auto *nsd1_84 = buffer.data(nsd1 + 84);
    const auto *nsd1_89 = buffer.data(nsd1 + 89);

    const auto *oss0_0 = buffer.data(oss0 + 0);
    const auto *oss0_1 = buffer.data(oss0 + 1);
    const auto *oss0_2 = buffer.data(oss0 + 2);
    const auto *oss0_3 = buffer.data(oss0 + 3);
    const auto *oss0_5 = buffer.data(oss0 + 5);
    const auto *oss0_6 = buffer.data(oss0 + 6);
    const auto *oss0_7 = buffer.data(oss0 + 7);
    const auto *oss0_8 = buffer.data(oss0 + 8);
    const auto *oss0_9 = buffer.data(oss0 + 9);
    const auto *oss0_10 = buffer.data(oss0 + 10);
    const auto *oss0_11 = buffer.data(oss0 + 11);
    const auto *oss0_12 = buffer.data(oss0 + 12);
    const auto *oss0_13 = buffer.data(oss0 + 13);
    const auto *oss0_14 = buffer.data(oss0 + 14);
    const auto *oss0_15 = buffer.data(oss0 + 15);
    const auto *oss0_16 = buffer.data(oss0 + 16);
    const auto *oss0_17 = buffer.data(oss0 + 17);
    const auto *oss0_18 = buffer.data(oss0 + 18);
    const auto *oss0_19 = buffer.data(oss0 + 19);
    const auto *oss0_20 = buffer.data(oss0 + 20);

    const auto *oss1_0 = buffer.data(oss1 + 0);
    const auto *oss1_1 = buffer.data(oss1 + 1);
    const auto *oss1_2 = buffer.data(oss1 + 2);
    const auto *oss1_3 = buffer.data(oss1 + 3);
    const auto *oss1_5 = buffer.data(oss1 + 5);
    const auto *oss1_6 = buffer.data(oss1 + 6);
    const auto *oss1_7 = buffer.data(oss1 + 7);
    const auto *oss1_8 = buffer.data(oss1 + 8);
    const auto *oss1_9 = buffer.data(oss1 + 9);
    const auto *oss1_10 = buffer.data(oss1 + 10);
    const auto *oss1_11 = buffer.data(oss1 + 11);
    const auto *oss1_12 = buffer.data(oss1 + 12);
    const auto *oss1_13 = buffer.data(oss1 + 13);
    const auto *oss1_14 = buffer.data(oss1 + 14);
    const auto *oss1_15 = buffer.data(oss1 + 15);
    const auto *oss1_16 = buffer.data(oss1 + 16);
    const auto *oss1_17 = buffer.data(oss1 + 17);
    const auto *oss1_18 = buffer.data(oss1 + 18);
    const auto *oss1_19 = buffer.data(oss1 + 19);
    const auto *oss1_20 = buffer.data(oss1 + 20);

    const auto *osp_0 = buffer.data(osp + 0);
    const auto *osp_1 = buffer.data(osp + 1);
    const auto *osp_2 = buffer.data(osp + 2);
    const auto *osp_3 = buffer.data(osp + 3);
    const auto *osp_4 = buffer.data(osp + 4);
    const auto *osp_6 = buffer.data(osp + 6);
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
    const auto *osp_47 = buffer.data(osp + 47);
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
    const auto *osp_61 = buffer.data(osp + 61);
    const auto *osp_62 = buffer.data(osp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, nsp_0, oss0_0, \
                         oss1_0, osp_0, osp_1, osp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * nsp_0[k]
                 + f_1 * oss0_0[k]
                 - f_2 * oss1_0[k]
                 + f_3 * pc_x[k] * osp_0[k];

        t_1[k] = f_3 * pc_y[k] * osp_0[k];

        t_2[k] = f_3 * pc_z[k] * osp_0[k];

        t_3[k] = f_1 * oss0_0[k]
                 - f_2 * oss1_0[k]
                 + f_3 * pc_y[k] * osp_1[k];

        t_4[k] = f_3 * pc_y[k] * osp_2[k];

        t_5[k] = f_1 * oss0_0[k]
                 - f_2 * oss1_0[k]
                 + f_3 * pc_z[k] * osp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, nsd0_0, nsp_1, nsp_4, \
                         nsd1_0, oss0_1, oss1_1, osp_3, osp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * nsd0_0[k]
                 - f_4 * pc_y[k] * nsd1_0[k];

        t_7[k] = f_5 * nsp_4[k]
                 + f_3 * pc_x[k] * osp_4[k];

        t_8[k] = f_3 * pc_z[k] * osp_3[k];

        t_9[k] = f_6 * nsp_1[k]
                 + f_1 * oss0_1[k]
                 - f_2 * oss1_1[k]
                 + f_3 * pc_y[k] * osp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, nsd0_0, nsd0_5, \
                         nsd1_0, nsd1_5, osp_4, osp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * osp_4[k];

        t_11[k] = pa_y[k] * nsd0_5[k]
                  - f_4 * pc_y[k] * nsd1_5[k];

        t_12[k] = pa_z[k] * nsd0_0[k]
                  - f_4 * pc_z[k] * nsd1_0[k];

        t_13[k] = f_3 * pc_y[k] * osp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, nsd0_3, nsp_2, nsp_8, \
                         nsd1_3, oss0_2, oss1_2, osp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * nsp_8[k]
                  + f_3 * pc_x[k] * osp_8[k];

        t_15[k] = pa_z[k] * nsd0_3[k]
                  - f_4 * pc_z[k] * nsd1_3[k];

        t_16[k] = f_3 * pc_y[k] * osp_8[k];

        t_17[k] = f_6 * nsp_2[k]
                  + f_1 * oss0_2[k]
                  - f_2 * oss1_2[k]
                  + f_3 * pc_z[k] * osp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, nsp_4, nsp_9, nsp_10, \
                         oss0_3, oss1_3, osp_9, osp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * nsp_9[k]
                  + f_1 * oss0_3[k]
                  - f_2 * oss1_3[k]
                  + f_3 * pc_x[k] * osp_9[k];

        t_19[k] = f_7 * nsp_10[k]
                  + f_3 * pc_x[k] * osp_10[k];

        t_20[k] = f_3 * pc_z[k] * osp_9[k];

        t_21[k] = f_8 * nsp_4[k]
                  + f_1 * oss0_3[k]
                  - f_2 * oss1_3[k]
                  + f_3 * pc_y[k] * osp_10[k];

        t_22[k] = f_3 * pc_z[k] * osp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, nsd0_12, nsp_13, nsd1_12, \
                         oss0_3, oss1_3, osp_11, osp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * oss0_3[k]
                  - f_2 * oss1_3[k]
                  + f_3 * pc_z[k] * osp_11[k];

        t_24[k] = pa_y[k] * nsd0_12[k]
                  - f_4 * pc_y[k] * nsd1_12[k];

        t_25[k] = f_7 * nsp_13[k]
                  + f_3 * pc_x[k] * osp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, nsd0_9, \
                         nsd0_17, nsp_8, nsp_14, nsd1_9, nsd1_17, \
                         osp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * nsp_14[k]
                  + f_3 * pc_x[k] * osp_14[k];

        t_27[k] = pa_z[k] * nsd0_9[k]
                  - f_4 * pc_z[k] * nsd1_9[k];

        t_28[k] = f_6 * nsp_8[k]
                  + f_3 * pc_y[k] * osp_14[k];

        t_29[k] = pa_y[k] * nsd0_17[k]
                  - f_4 * pc_y[k] * nsd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, nsp_15, nsp_17, oss0_5, \
                         oss1_5, osp_15, osp_16, osp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * nsp_15[k]
                  + f_1 * oss0_5[k]
                  - f_2 * oss1_5[k]
                  + f_3 * pc_x[k] * osp_15[k];

        t_31[k] = f_3 * pc_y[k] * osp_15[k];

        t_32[k] = f_7 * nsp_17[k]
                  + f_3 * pc_x[k] * osp_17[k];

        t_33[k] = f_1 * oss0_5[k]
                  - f_2 * oss1_5[k]
                  + f_3 * pc_y[k] * osp_16[k];

        t_34[k] = f_3 * pc_y[k] * osp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, nsp_8, nsp_18, nsp_19, oss0_5, \
                         oss0_6, oss1_5, oss1_6, osp_17, osp_18, \
                         osp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * nsp_8[k]
                  + f_1 * oss0_5[k]
                  - f_2 * oss1_5[k]
                  + f_3 * pc_z[k] * osp_17[k];

        t_36[k] = f_9 * nsp_18[k]
                  + f_1 * oss0_6[k]
                  - f_2 * oss1_6[k]
                  + f_3 * pc_x[k] * osp_18[k];

        t_37[k] = f_9 * nsp_19[k]
                  + f_3 * pc_x[k] * osp_19[k];

        t_38[k] = f_3 * pc_z[k] * osp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, nsd0_18, nsp_10, nsd1_18, \
                         oss0_6, oss1_6, osp_19, osp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * nsp_10[k]
                  + f_1 * oss0_6[k]
                  - f_2 * oss1_6[k]
                  + f_3 * pc_y[k] * osp_19[k];

        t_40[k] = f_3 * pc_z[k] * osp_19[k];

        t_41[k] = f_1 * oss0_6[k]
                  - f_2 * oss1_6[k]
                  + f_3 * pc_z[k] * osp_20[k];

        t_42[k] = pa_z[k] * nsd0_18[k]
                  - f_4 * pc_z[k] * nsd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, nsd0_21, nsp_14, \
                         nsp_22, nsp_23, nsd1_21, osp_22, osp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * nsp_22[k]
                  + f_3 * pc_x[k] * osp_22[k];

        t_44[k] = f_9 * nsp_23[k]
                  + f_3 * pc_x[k] * osp_23[k];

        t_45[k] = pa_z[k] * nsd0_21[k]
                  - f_4 * pc_z[k] * nsd1_21[k];

        t_46[k] = f_8 * nsp_14[k]
                  + f_3 * pc_y[k] * osp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, nsd0_30, nsp_11, nsp_25, \
                         nsd1_30, oss0_7, oss1_7, osp_23, osp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * nsp_11[k]
                  + f_1 * oss0_7[k]
                  - f_2 * oss1_7[k]
                  + f_3 * pc_z[k] * osp_23[k];

        t_48[k] = pa_y[k] * nsd0_30[k]
                  - f_4 * pc_y[k] * nsd1_30[k];

        t_49[k] = f_9 * nsp_25[k]
                  + f_3 * pc_x[k] * osp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, nsd0_35, nsp_16, nsp_17, \
                         nsp_26, nsd1_35, oss0_8, oss1_8, osp_25, \
                         osp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * nsp_26[k]
                  + f_3 * pc_x[k] * osp_26[k];

        t_51[k] = f_6 * nsp_16[k]
                  + f_1 * oss0_8[k]
                  - f_2 * oss1_8[k]
                  + f_3 * pc_y[k] * osp_25[k];

        t_52[k] = f_6 * nsp_17[k]
                  + f_3 * pc_y[k] * osp_26[k];

        t_53[k] = pa_y[k] * nsd0_35[k]
                  - f_4 * pc_y[k] * nsd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, nsp_27, nsp_29, oss0_9, \
                         oss1_9, osp_27, osp_28, osp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * nsp_27[k]
                  + f_1 * oss0_9[k]
                  - f_2 * oss1_9[k]
                  + f_3 * pc_x[k] * osp_27[k];

        t_55[k] = f_3 * pc_y[k] * osp_27[k];

        t_56[k] = f_9 * nsp_29[k]
                  + f_3 * pc_x[k] * osp_29[k];

        t_57[k] = f_1 * oss0_9[k]
                  - f_2 * oss1_9[k]
                  + f_3 * pc_y[k] * osp_28[k];

        t_58[k] = f_3 * pc_y[k] * osp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, nsp_17, nsp_30, nsp_31, oss0_9, \
                         oss0_10, oss1_9, oss1_10, osp_29, osp_30, \
                         osp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * nsp_17[k]
                  + f_1 * oss0_9[k]
                  - f_2 * oss1_9[k]
                  + f_3 * pc_z[k] * osp_29[k];

        t_60[k] = f_11 * nsp_30[k]
                  + f_1 * oss0_10[k]
                  - f_2 * oss1_10[k]
                  + f_3 * pc_x[k] * osp_30[k];

        t_61[k] = f_11 * nsp_31[k]
                  + f_3 * pc_x[k] * osp_31[k];

        t_62[k] = f_3 * pc_z[k] * osp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, nsd0_36, nsp_19, nsd1_36, \
                         oss0_10, oss1_10, osp_31, osp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_12 * nsp_19[k]
                  + f_1 * oss0_10[k]
                  - f_2 * oss1_10[k]
                  + f_3 * pc_y[k] * osp_31[k];

        t_64[k] = f_3 * pc_z[k] * osp_31[k];

        t_65[k] = f_1 * oss0_10[k]
                  - f_2 * oss1_10[k]
                  + f_3 * pc_z[k] * osp_32[k];

        t_66[k] = pa_z[k] * nsd0_36[k]
                  - f_4 * pc_z[k] * nsd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, nsd0_39, nsp_23, \
                         nsp_34, nsp_35, nsd1_39, osp_34, osp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * nsp_34[k]
                  + f_3 * pc_x[k] * osp_34[k];

        t_68[k] = f_11 * nsp_35[k]
                  + f_3 * pc_x[k] * osp_35[k];

        t_69[k] = pa_z[k] * nsd0_39[k]
                  - f_4 * pc_z[k] * nsd1_39[k];

        t_70[k] = f_10 * nsp_23[k]
                  + f_3 * pc_y[k] * osp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, nsp_20, nsp_36, nsp_37, oss0_11, \
                         oss0_12, oss1_11, oss1_12, osp_35, osp_36, \
                         osp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * nsp_20[k]
                  + f_1 * oss0_11[k]
                  - f_2 * oss1_11[k]
                  + f_3 * pc_z[k] * osp_35[k];

        t_72[k] = f_11 * nsp_36[k]
                  + f_1 * oss0_12[k]
                  - f_2 * oss1_12[k]
                  + f_3 * pc_x[k] * osp_36[k];

        t_73[k] = f_11 * nsp_37[k]
                  + f_3 * pc_x[k] * osp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, nsp_23, nsp_25, nsp_26, \
                         nsp_38, oss0_12, oss1_12, osp_37, osp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * nsp_38[k]
                  + f_3 * pc_x[k] * osp_38[k];

        t_75[k] = f_8 * nsp_25[k]
                  + f_1 * oss0_12[k]
                  - f_2 * oss1_12[k]
                  + f_3 * pc_y[k] * osp_37[k];

        t_76[k] = f_8 * nsp_26[k]
                  + f_3 * pc_y[k] * osp_38[k];

        t_77[k] = f_8 * nsp_23[k]
                  + f_1 * oss0_12[k]
                  - f_2 * oss1_12[k]
                  + f_3 * pc_z[k] * osp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, nsd0_54, nsp_28, nsp_40, \
                         nsp_41, nsd1_54, oss0_13, oss1_13, osp_40, \
                         osp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * nsd0_54[k]
                  - f_4 * pc_y[k] * nsd1_54[k];

        t_79[k] = f_11 * nsp_40[k]
                  + f_3 * pc_x[k] * osp_40[k];

        t_80[k] = f_11 * nsp_41[k]
                  + f_3 * pc_x[k] * osp_41[k];

        t_81[k] = f_6 * nsp_28[k]
                  + f_1 * oss0_13[k]
                  - f_2 * oss1_13[k]
                  + f_3 * pc_y[k] * osp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, nsd0_59, nsp_29, nsp_42, \
                         nsd1_59, oss0_14, oss1_14, osp_41, osp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * nsp_29[k]
                  + f_3 * pc_y[k] * osp_41[k];

        t_83[k] = pa_y[k] * nsd0_59[k]
                  - f_4 * pc_y[k] * nsd1_59[k];

        t_84[k] = f_11 * nsp_42[k]
                  + f_1 * oss0_14[k]
                  - f_2 * oss1_14[k]
                  + f_3 * pc_x[k] * osp_42[k];

        t_85[k] = f_3 * pc_y[k] * osp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, nsp_29, nsp_44, oss0_14, \
                         oss1_14, osp_43, osp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * nsp_44[k]
                  + f_3 * pc_x[k] * osp_44[k];

        t_87[k] = f_1 * oss0_14[k]
                  - f_2 * oss1_14[k]
                  + f_3 * pc_y[k] * osp_43[k];

        t_88[k] = f_3 * pc_y[k] * osp_44[k];

        t_89[k] = f_12 * nsp_29[k]
                  + f_1 * oss0_14[k]
                  - f_2 * oss1_14[k]
                  + f_3 * pc_z[k] * osp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, nsp_31, nsp_45, \
                         nsp_46, oss0_15, oss1_15, osp_45, osp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_13 * nsp_45[k]
                  + f_1 * oss0_15[k]
                  - f_2 * oss1_15[k]
                  + f_3 * pc_x[k] * osp_45[k];

        t_91[k] = f_13 * nsp_46[k]
                  + f_3 * pc_x[k] * osp_46[k];

        t_92[k] = f_3 * pc_z[k] * osp_45[k];

        t_93[k] = f_14 * nsp_31[k]
                  + f_1 * oss0_15[k]
                  - f_2 * oss1_15[k]
                  + f_3 * pc_y[k] * osp_46[k];

        t_94[k] = f_3 * pc_z[k] * osp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pc_x, pc_z, nsd0_60, nsp_49, nsp_50, \
                         nsd1_60, oss0_15, oss1_15, osp_47, osp_49, \
                         osp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_1 * oss0_15[k]
                  - f_2 * oss1_15[k]
                  + f_3 * pc_z[k] * osp_47[k];

        t_96[k] = pa_z[k] * nsd0_60[k]
                  - f_4 * pc_z[k] * nsd1_60[k];

        t_97[k] = f_13 * nsp_49[k]
                  + f_3 * pc_x[k] * osp_49[k];

        t_98[k] = f_13 * nsp_50[k]
                  + f_3 * pc_x[k] * osp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, nsd0_63, nsp_32, nsp_35, \
                         nsd1_63, oss0_16, oss1_16, osp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * nsd0_63[k]
                  - f_4 * pc_z[k] * nsd1_63[k];

        t_100[k] = f_12 * nsp_35[k]
                   + f_3 * pc_y[k] * osp_50[k];

        t_101[k] = f_6 * nsp_32[k]
                   + f_1 * oss0_16[k]
                   - f_2 * oss1_16[k]
                   + f_3 * pc_z[k] * osp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, nsp_37, nsp_51, nsp_52, \
                         nsp_53, oss0_17, oss1_17, osp_51, osp_52, \
                         osp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_13 * nsp_51[k]
                   + f_1 * oss0_17[k]
                   - f_2 * oss1_17[k]
                   + f_3 * pc_x[k] * osp_51[k];

        t_103[k] = f_13 * nsp_52[k]
                   + f_3 * pc_x[k] * osp_52[k];

        t_104[k] = f_13 * nsp_53[k]
                   + f_3 * pc_x[k] * osp_53[k];

        t_105[k] = f_10 * nsp_37[k]
                   + f_1 * oss0_17[k]
                   - f_2 * oss1_17[k]
                   + f_3 * pc_y[k] * osp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, nsp_35, nsp_38, nsp_54, \
                         oss0_17, oss0_18, oss1_17, oss1_18, osp_53, \
                         osp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * nsp_38[k]
                   + f_3 * pc_y[k] * osp_53[k];

        t_107[k] = f_8 * nsp_35[k]
                   + f_1 * oss0_17[k]
                   - f_2 * oss1_17[k]
                   + f_3 * pc_z[k] * osp_53[k];

        t_108[k] = f_13 * nsp_54[k]
                   + f_1 * oss0_18[k]
                   - f_2 * oss1_18[k]
                   + f_3 * pc_x[k] * osp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, nsp_40, nsp_41, nsp_55, \
                         nsp_56, oss0_18, oss1_18, osp_55, osp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_13 * nsp_55[k]
                   + f_3 * pc_x[k] * osp_55[k];

        t_110[k] = f_13 * nsp_56[k]
                   + f_3 * pc_x[k] * osp_56[k];

        t_111[k] = f_8 * nsp_40[k]
                   + f_1 * oss0_18[k]
                   - f_2 * oss1_18[k]
                   + f_3 * pc_y[k] * osp_55[k];

        t_112[k] = f_8 * nsp_41[k]
                   + f_3 * pc_y[k] * osp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, nsd0_84, nsp_38, nsp_58, \
                         nsd1_84, oss0_18, oss1_18, osp_56, osp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * nsp_38[k]
                   + f_1 * oss0_18[k]
                   - f_2 * oss1_18[k]
                   + f_3 * pc_z[k] * osp_56[k];

        t_114[k] = pa_y[k] * nsd0_84[k]
                   - f_4 * pc_y[k] * nsd1_84[k];

        t_115[k] = f_13 * nsp_58[k]
                   + f_3 * pc_x[k] * osp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, nsd0_89, nsp_43, \
                         nsp_44, nsp_59, nsd1_89, oss0_19, oss1_19, osp_58, \
                         osp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_13 * nsp_59[k]
                   + f_3 * pc_x[k] * osp_59[k];

        t_117[k] = f_6 * nsp_43[k]
                   + f_1 * oss0_19[k]
                   - f_2 * oss1_19[k]
                   + f_3 * pc_y[k] * osp_58[k];

        t_118[k] = f_6 * nsp_44[k]
                   + f_3 * pc_y[k] * osp_59[k];

        t_119[k] = pa_y[k] * nsd0_89[k]
                   - f_4 * pc_y[k] * nsd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, nsp_60, nsp_62, \
                         oss0_20, oss1_20, osp_60, osp_61, osp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_13 * nsp_60[k]
                   + f_1 * oss0_20[k]
                   - f_2 * oss1_20[k]
                   + f_3 * pc_x[k] * osp_60[k];

        t_121[k] = f_3 * pc_y[k] * osp_60[k];

        t_122[k] = f_13 * nsp_62[k]
                   + f_3 * pc_x[k] * osp_62[k];

        t_123[k] = f_1 * oss0_20[k]
                   - f_2 * oss1_20[k]
                   + f_3 * pc_y[k] * osp_61[k];

        t_124[k] = f_3 * pc_y[k] * osp_62[k];
    }
}

static auto
compute_prim_osd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsd0,
                                                          const size_t nsp, const size_t nsd1,
                                                          const size_t oss0, const size_t oss1,
                                                          const size_t osp, const size_t ncols,
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
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *nsd0_90 = buffer.data(nsd0 + 90);
    const auto *nsd0_93 = buffer.data(nsd0 + 93);
    const auto *nsd0_120 = buffer.data(nsd0 + 120);
    const auto *nsd0_125 = buffer.data(nsd0 + 125);
    const auto *nsd0_126 = buffer.data(nsd0 + 126);
    const auto *nsd0_129 = buffer.data(nsd0 + 129);
    const auto *nsd0_162 = buffer.data(nsd0 + 162);
    const auto *nsd0_167 = buffer.data(nsd0 + 167);
    const auto *nsd0_168 = buffer.data(nsd0 + 168);
    const auto *nsd0_171 = buffer.data(nsd0 + 171);

    const auto *nsp_44 = buffer.data(nsp + 44);
    const auto *nsp_46 = buffer.data(nsp + 46);
    const auto *nsp_47 = buffer.data(nsp + 47);
    const auto *nsp_50 = buffer.data(nsp + 50);
    const auto *nsp_52 = buffer.data(nsp + 52);
    const auto *nsp_53 = buffer.data(nsp + 53);
    const auto *nsp_55 = buffer.data(nsp + 55);
    const auto *nsp_56 = buffer.data(nsp + 56);
    const auto *nsp_58 = buffer.data(nsp + 58);
    const auto *nsp_59 = buffer.data(nsp + 59);
    const auto *nsp_61 = buffer.data(nsp + 61);
    const auto *nsp_62 = buffer.data(nsp + 62);
    const auto *nsp_63 = buffer.data(nsp + 63);
    const auto *nsp_64 = buffer.data(nsp + 64);
    const auto *nsp_65 = buffer.data(nsp + 65);
    const auto *nsp_67 = buffer.data(nsp + 67);
    const auto *nsp_68 = buffer.data(nsp + 68);
    const auto *nsp_69 = buffer.data(nsp + 69);
    const auto *nsp_70 = buffer.data(nsp + 70);
    const auto *nsp_71 = buffer.data(nsp + 71);
    const auto *nsp_72 = buffer.data(nsp + 72);
    const auto *nsp_73 = buffer.data(nsp + 73);
    const auto *nsp_74 = buffer.data(nsp + 74);
    const auto *nsp_75 = buffer.data(nsp + 75);
    const auto *nsp_76 = buffer.data(nsp + 76);
    const auto *nsp_77 = buffer.data(nsp + 77);
    const auto *nsp_79 = buffer.data(nsp + 79);
    const auto *nsp_80 = buffer.data(nsp + 80);
    const auto *nsp_81 = buffer.data(nsp + 81);
    const auto *nsp_82 = buffer.data(nsp + 82);
    const auto *nsp_83 = buffer.data(nsp + 83);
    const auto *nsp_84 = buffer.data(nsp + 84);
    const auto *nsp_85 = buffer.data(nsp + 85);
    const auto *nsp_86 = buffer.data(nsp + 86);
    const auto *nsp_88 = buffer.data(nsp + 88);
    const auto *nsp_89 = buffer.data(nsp + 89);
    const auto *nsp_90 = buffer.data(nsp + 90);
    const auto *nsp_91 = buffer.data(nsp + 91);
    const auto *nsp_92 = buffer.data(nsp + 92);
    const auto *nsp_93 = buffer.data(nsp + 93);
    const auto *nsp_94 = buffer.data(nsp + 94);
    const auto *nsp_95 = buffer.data(nsp + 95);
    const auto *nsp_96 = buffer.data(nsp + 96);
    const auto *nsp_97 = buffer.data(nsp + 97);
    const auto *nsp_98 = buffer.data(nsp + 98);
    const auto *nsp_99 = buffer.data(nsp + 99);
    const auto *nsp_100 = buffer.data(nsp + 100);
    const auto *nsp_101 = buffer.data(nsp + 101);
    const auto *nsp_103 = buffer.data(nsp + 103);
    const auto *nsp_104 = buffer.data(nsp + 104);
    const auto *nsp_105 = buffer.data(nsp + 105);
    const auto *nsp_107 = buffer.data(nsp + 107);
    const auto *nsp_108 = buffer.data(nsp + 108);
    const auto *nsp_109 = buffer.data(nsp + 109);
    const auto *nsp_112 = buffer.data(nsp + 112);
    const auto *nsp_113 = buffer.data(nsp + 113);
    const auto *nsp_114 = buffer.data(nsp + 114);
    const auto *nsp_115 = buffer.data(nsp + 115);
    const auto *nsp_116 = buffer.data(nsp + 116);
    const auto *nsp_117 = buffer.data(nsp + 117);
    const auto *nsp_118 = buffer.data(nsp + 118);
    const auto *nsp_119 = buffer.data(nsp + 119);
    const auto *nsp_120 = buffer.data(nsp + 120);
    const auto *nsp_121 = buffer.data(nsp + 121);

    const auto *nsd1_90 = buffer.data(nsd1 + 90);
    const auto *nsd1_93 = buffer.data(nsd1 + 93);
    const auto *nsd1_120 = buffer.data(nsd1 + 120);
    const auto *nsd1_125 = buffer.data(nsd1 + 125);
    const auto *nsd1_126 = buffer.data(nsd1 + 126);
    const auto *nsd1_129 = buffer.data(nsd1 + 129);
    const auto *nsd1_162 = buffer.data(nsd1 + 162);
    const auto *nsd1_167 = buffer.data(nsd1 + 167);
    const auto *nsd1_168 = buffer.data(nsd1 + 168);
    const auto *nsd1_171 = buffer.data(nsd1 + 171);

    const auto *oss0_20 = buffer.data(oss0 + 20);
    const auto *oss0_21 = buffer.data(oss0 + 21);
    const auto *oss0_22 = buffer.data(oss0 + 22);
    const auto *oss0_23 = buffer.data(oss0 + 23);
    const auto *oss0_24 = buffer.data(oss0 + 24);
    const auto *oss0_25 = buffer.data(oss0 + 25);
    const auto *oss0_26 = buffer.data(oss0 + 26);
    const auto *oss0_27 = buffer.data(oss0 + 27);
    const auto *oss0_28 = buffer.data(oss0 + 28);
    const auto *oss0_29 = buffer.data(oss0 + 29);
    const auto *oss0_30 = buffer.data(oss0 + 30);
    const auto *oss0_31 = buffer.data(oss0 + 31);
    const auto *oss0_32 = buffer.data(oss0 + 32);
    const auto *oss0_33 = buffer.data(oss0 + 33);
    const auto *oss0_34 = buffer.data(oss0 + 34);
    const auto *oss0_35 = buffer.data(oss0 + 35);
    const auto *oss0_36 = buffer.data(oss0 + 36);
    const auto *oss0_37 = buffer.data(oss0 + 37);
    const auto *oss0_38 = buffer.data(oss0 + 38);
    const auto *oss0_39 = buffer.data(oss0 + 39);
    const auto *oss0_40 = buffer.data(oss0 + 40);

    const auto *oss1_20 = buffer.data(oss1 + 20);
    const auto *oss1_21 = buffer.data(oss1 + 21);
    const auto *oss1_22 = buffer.data(oss1 + 22);
    const auto *oss1_23 = buffer.data(oss1 + 23);
    const auto *oss1_24 = buffer.data(oss1 + 24);
    const auto *oss1_25 = buffer.data(oss1 + 25);
    const auto *oss1_26 = buffer.data(oss1 + 26);
    const auto *oss1_27 = buffer.data(oss1 + 27);
    const auto *oss1_28 = buffer.data(oss1 + 28);
    const auto *oss1_29 = buffer.data(oss1 + 29);
    const auto *oss1_30 = buffer.data(oss1 + 30);
    const auto *oss1_31 = buffer.data(oss1 + 31);
    const auto *oss1_32 = buffer.data(oss1 + 32);
    const auto *oss1_33 = buffer.data(oss1 + 33);
    const auto *oss1_34 = buffer.data(oss1 + 34);
    const auto *oss1_35 = buffer.data(oss1 + 35);
    const auto *oss1_36 = buffer.data(oss1 + 36);
    const auto *oss1_37 = buffer.data(oss1 + 37);
    const auto *oss1_38 = buffer.data(oss1 + 38);
    const auto *oss1_39 = buffer.data(oss1 + 39);
    const auto *oss1_40 = buffer.data(oss1 + 40);

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
    const auto *osp_106 = buffer.data(osp + 106);
    const auto *osp_107 = buffer.data(osp + 107);
    const auto *osp_108 = buffer.data(osp + 108);
    const auto *osp_109 = buffer.data(osp + 109);
    const auto *osp_110 = buffer.data(osp + 110);
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

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_z, nsp_44, nsp_63, nsp_64, \
                         oss0_20, oss0_21, oss1_20, oss1_21, osp_62, osp_63, \
                         osp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_14 * nsp_44[k]
                   + f_1 * oss0_20[k]
                   - f_2 * oss1_20[k]
                   + f_3 * pc_z[k] * osp_62[k];

        t_126[k] = f_14 * nsp_63[k]
                   + f_1 * oss0_21[k]
                   - f_2 * oss1_21[k]
                   + f_3 * pc_x[k] * osp_63[k];

        t_127[k] = f_14 * nsp_64[k]
                   + f_3 * pc_x[k] * osp_64[k];

        t_128[k] = f_3 * pc_z[k] * osp_63[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_z, pc_y, pc_z, nsd0_90, nsp_46, \
                         nsd1_90, oss0_21, oss1_21, osp_64, osp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_13 * nsp_46[k]
                   + f_1 * oss0_21[k]
                   - f_2 * oss1_21[k]
                   + f_3 * pc_y[k] * osp_64[k];

        t_130[k] = f_3 * pc_z[k] * osp_64[k];

        t_131[k] = f_1 * oss0_21[k]
                   - f_2 * oss1_21[k]
                   + f_3 * pc_z[k] * osp_65[k];

        t_132[k] = pa_z[k] * nsd0_90[k]
                   - f_4 * pc_z[k] * nsd1_90[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_z, pc_x, pc_y, pc_z, nsd0_93, nsp_50, \
                         nsp_67, nsp_68, nsd1_93, osp_67, osp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_14 * nsp_67[k]
                   + f_3 * pc_x[k] * osp_67[k];

        t_134[k] = f_14 * nsp_68[k]
                   + f_3 * pc_x[k] * osp_68[k];

        t_135[k] = pa_z[k] * nsd0_93[k]
                   - f_4 * pc_z[k] * nsd1_93[k];

        t_136[k] = f_14 * nsp_50[k]
                   + f_3 * pc_y[k] * osp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_x, pc_z, nsp_47, nsp_69, nsp_70, oss0_22, \
                         oss0_23, oss1_22, oss1_23, osp_68, osp_69, \
                         osp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * nsp_47[k]
                   + f_1 * oss0_22[k]
                   - f_2 * oss1_22[k]
                   + f_3 * pc_z[k] * osp_68[k];

        t_138[k] = f_14 * nsp_69[k]
                   + f_1 * oss0_23[k]
                   - f_2 * oss1_23[k]
                   + f_3 * pc_x[k] * osp_69[k];

        t_139[k] = f_14 * nsp_70[k]
                   + f_3 * pc_x[k] * osp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, nsp_50, nsp_52, nsp_53, \
                         nsp_71, oss0_23, oss1_23, osp_70, osp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_14 * nsp_71[k]
                   + f_3 * pc_x[k] * osp_71[k];

        t_141[k] = f_12 * nsp_52[k]
                   + f_1 * oss0_23[k]
                   - f_2 * oss1_23[k]
                   + f_3 * pc_y[k] * osp_70[k];

        t_142[k] = f_12 * nsp_53[k]
                   + f_3 * pc_y[k] * osp_71[k];

        t_143[k] = f_8 * nsp_50[k]
                   + f_1 * oss0_23[k]
                   - f_2 * oss1_23[k]
                   + f_3 * pc_z[k] * osp_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pc_x, pc_y, nsp_55, nsp_72, nsp_73, \
                         nsp_74, oss0_24, oss1_24, osp_72, osp_73, \
                         osp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_14 * nsp_72[k]
                   + f_1 * oss0_24[k]
                   - f_2 * oss1_24[k]
                   + f_3 * pc_x[k] * osp_72[k];

        t_145[k] = f_14 * nsp_73[k]
                   + f_3 * pc_x[k] * osp_73[k];

        t_146[k] = f_14 * nsp_74[k]
                   + f_3 * pc_x[k] * osp_74[k];

        t_147[k] = f_10 * nsp_55[k]
                   + f_1 * oss0_24[k]
                   - f_2 * oss1_24[k]
                   + f_3 * pc_y[k] * osp_73[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pc_x, pc_y, pc_z, nsp_53, nsp_56, nsp_75, \
                         oss0_24, oss0_25, oss1_24, oss1_25, osp_74, \
                         osp_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * nsp_56[k]
                   + f_3 * pc_y[k] * osp_74[k];

        t_149[k] = f_10 * nsp_53[k]
                   + f_1 * oss0_24[k]
                   - f_2 * oss1_24[k]
                   + f_3 * pc_z[k] * osp_74[k];

        t_150[k] = f_14 * nsp_75[k]
                   + f_1 * oss0_25[k]
                   - f_2 * oss1_25[k]
                   + f_3 * pc_x[k] * osp_75[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, nsp_58, nsp_59, nsp_76, \
                         nsp_77, oss0_25, oss1_25, osp_76, osp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_14 * nsp_76[k]
                   + f_3 * pc_x[k] * osp_76[k];

        t_152[k] = f_14 * nsp_77[k]
                   + f_3 * pc_x[k] * osp_77[k];

        t_153[k] = f_8 * nsp_58[k]
                   + f_1 * oss0_25[k]
                   - f_2 * oss1_25[k]
                   + f_3 * pc_y[k] * osp_76[k];

        t_154[k] = f_8 * nsp_59[k]
                   + f_3 * pc_y[k] * osp_77[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pa_y, pc_x, pc_y, pc_z, nsd0_120, nsp_56, \
                         nsp_79, nsd1_120, oss0_25, oss1_25, osp_77, \
                         osp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * nsp_56[k]
                   + f_1 * oss0_25[k]
                   - f_2 * oss1_25[k]
                   + f_3 * pc_z[k] * osp_77[k];

        t_156[k] = pa_y[k] * nsd0_120[k]
                   - f_4 * pc_y[k] * nsd1_120[k];

        t_157[k] = f_14 * nsp_79[k]
                   + f_3 * pc_x[k] * osp_79[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_y, pc_x, pc_y, nsd0_125, nsp_61, \
                         nsp_62, nsp_80, nsd1_125, oss0_26, oss1_26, osp_79, \
                         osp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_14 * nsp_80[k]
                   + f_3 * pc_x[k] * osp_80[k];

        t_159[k] = f_6 * nsp_61[k]
                   + f_1 * oss0_26[k]
                   - f_2 * oss1_26[k]
                   + f_3 * pc_y[k] * osp_79[k];

        t_160[k] = f_6 * nsp_62[k]
                   + f_3 * pc_y[k] * osp_80[k];

        t_161[k] = pa_y[k] * nsd0_125[k]
                   - f_4 * pc_y[k] * nsd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pc_x, pc_y, nsp_81, nsp_83, \
                         oss0_27, oss1_27, osp_81, osp_82, osp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_14 * nsp_81[k]
                   + f_1 * oss0_27[k]
                   - f_2 * oss1_27[k]
                   + f_3 * pc_x[k] * osp_81[k];

        t_163[k] = f_3 * pc_y[k] * osp_81[k];

        t_164[k] = f_14 * nsp_83[k]
                   + f_3 * pc_x[k] * osp_83[k];

        t_165[k] = f_1 * oss0_27[k]
                   - f_2 * oss1_27[k]
                   + f_3 * pc_y[k] * osp_82[k];

        t_166[k] = f_3 * pc_y[k] * osp_83[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pc_x, pc_z, nsp_62, nsp_84, nsp_85, \
                         oss0_27, oss0_28, oss1_27, oss1_28, osp_83, osp_84, \
                         osp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_13 * nsp_62[k]
                   + f_1 * oss0_27[k]
                   - f_2 * oss1_27[k]
                   + f_3 * pc_z[k] * osp_83[k];

        t_168[k] = f_12 * nsp_84[k]
                   + f_1 * oss0_28[k]
                   - f_2 * oss1_28[k]
                   + f_3 * pc_x[k] * osp_84[k];

        t_169[k] = f_12 * nsp_85[k]
                   + f_3 * pc_x[k] * osp_85[k];

        t_170[k] = f_3 * pc_z[k] * osp_84[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_z, pc_y, pc_z, nsd0_126, nsp_64, \
                         nsd1_126, oss0_28, oss1_28, osp_85, osp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_11 * nsp_64[k]
                   + f_1 * oss0_28[k]
                   - f_2 * oss1_28[k]
                   + f_3 * pc_y[k] * osp_85[k];

        t_172[k] = f_3 * pc_z[k] * osp_85[k];

        t_173[k] = f_1 * oss0_28[k]
                   - f_2 * oss1_28[k]
                   + f_3 * pc_z[k] * osp_86[k];

        t_174[k] = pa_z[k] * nsd0_126[k]
                   - f_4 * pc_z[k] * nsd1_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_z, pc_x, pc_y, pc_z, nsd0_129, nsp_68, \
                         nsp_88, nsp_89, nsd1_129, osp_88, osp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_12 * nsp_88[k]
                   + f_3 * pc_x[k] * osp_88[k];

        t_176[k] = f_12 * nsp_89[k]
                   + f_3 * pc_x[k] * osp_89[k];

        t_177[k] = pa_z[k] * nsd0_129[k]
                   - f_4 * pc_z[k] * nsd1_129[k];

        t_178[k] = f_13 * nsp_68[k]
                   + f_3 * pc_y[k] * osp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_x, pc_z, nsp_65, nsp_90, nsp_91, oss0_29, \
                         oss0_30, oss1_29, oss1_30, osp_89, osp_90, \
                         osp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * nsp_65[k]
                   + f_1 * oss0_29[k]
                   - f_2 * oss1_29[k]
                   + f_3 * pc_z[k] * osp_89[k];

        t_180[k] = f_12 * nsp_90[k]
                   + f_1 * oss0_30[k]
                   - f_2 * oss1_30[k]
                   + f_3 * pc_x[k] * osp_90[k];

        t_181[k] = f_12 * nsp_91[k]
                   + f_3 * pc_x[k] * osp_91[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, nsp_68, nsp_70, nsp_71, \
                         nsp_92, oss0_30, oss1_30, osp_91, osp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_12 * nsp_92[k]
                   + f_3 * pc_x[k] * osp_92[k];

        t_183[k] = f_14 * nsp_70[k]
                   + f_1 * oss0_30[k]
                   - f_2 * oss1_30[k]
                   + f_3 * pc_y[k] * osp_91[k];

        t_184[k] = f_14 * nsp_71[k]
                   + f_3 * pc_y[k] * osp_92[k];

        t_185[k] = f_8 * nsp_68[k]
                   + f_1 * oss0_30[k]
                   - f_2 * oss1_30[k]
                   + f_3 * pc_z[k] * osp_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, pc_y, nsp_73, nsp_93, nsp_94, \
                         nsp_95, oss0_31, oss1_31, osp_93, osp_94, \
                         osp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_12 * nsp_93[k]
                   + f_1 * oss0_31[k]
                   - f_2 * oss1_31[k]
                   + f_3 * pc_x[k] * osp_93[k];

        t_187[k] = f_12 * nsp_94[k]
                   + f_3 * pc_x[k] * osp_94[k];

        t_188[k] = f_12 * nsp_95[k]
                   + f_3 * pc_x[k] * osp_95[k];

        t_189[k] = f_12 * nsp_73[k]
                   + f_1 * oss0_31[k]
                   - f_2 * oss1_31[k]
                   + f_3 * pc_y[k] * osp_94[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_y, pc_z, nsp_71, nsp_74, nsp_96, \
                         oss0_31, oss0_32, oss1_31, oss1_32, osp_95, \
                         osp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_12 * nsp_74[k]
                   + f_3 * pc_y[k] * osp_95[k];

        t_191[k] = f_10 * nsp_71[k]
                   + f_1 * oss0_31[k]
                   - f_2 * oss1_31[k]
                   + f_3 * pc_z[k] * osp_95[k];

        t_192[k] = f_12 * nsp_96[k]
                   + f_1 * oss0_32[k]
                   - f_2 * oss1_32[k]
                   + f_3 * pc_x[k] * osp_96[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, nsp_76, nsp_77, nsp_97, \
                         nsp_98, oss0_32, oss1_32, osp_97, osp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_12 * nsp_97[k]
                   + f_3 * pc_x[k] * osp_97[k];

        t_194[k] = f_12 * nsp_98[k]
                   + f_3 * pc_x[k] * osp_98[k];

        t_195[k] = f_10 * nsp_76[k]
                   + f_1 * oss0_32[k]
                   - f_2 * oss1_32[k]
                   + f_3 * pc_y[k] * osp_97[k];

        t_196[k] = f_10 * nsp_77[k]
                   + f_3 * pc_y[k] * osp_98[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pc_x, pc_z, nsp_74, nsp_99, nsp_100, oss0_32, \
                         oss0_33, oss1_32, oss1_33, osp_98, osp_99, \
                         osp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * nsp_74[k]
                   + f_1 * oss0_32[k]
                   - f_2 * oss1_32[k]
                   + f_3 * pc_z[k] * osp_98[k];

        t_198[k] = f_12 * nsp_99[k]
                   + f_1 * oss0_33[k]
                   - f_2 * oss1_33[k]
                   + f_3 * pc_x[k] * osp_99[k];

        t_199[k] = f_12 * nsp_100[k]
                   + f_3 * pc_x[k] * osp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, nsp_77, nsp_79, nsp_80, \
                         nsp_101, oss0_33, oss1_33, osp_100, osp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * nsp_101[k]
                   + f_3 * pc_x[k] * osp_101[k];

        t_201[k] = f_8 * nsp_79[k]
                   + f_1 * oss0_33[k]
                   - f_2 * oss1_33[k]
                   + f_3 * pc_y[k] * osp_100[k];

        t_202[k] = f_8 * nsp_80[k]
                   + f_3 * pc_y[k] * osp_101[k];

        t_203[k] = f_14 * nsp_77[k]
                   + f_1 * oss0_33[k]
                   - f_2 * oss1_33[k]
                   + f_3 * pc_z[k] * osp_101[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_x, pc_y, nsd0_162, nsp_82, \
                         nsp_103, nsp_104, nsd1_162, oss0_34, oss1_34, osp_103, \
                         osp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * nsd0_162[k]
                   - f_4 * pc_y[k] * nsd1_162[k];

        t_205[k] = f_12 * nsp_103[k]
                   + f_3 * pc_x[k] * osp_103[k];

        t_206[k] = f_12 * nsp_104[k]
                   + f_3 * pc_x[k] * osp_104[k];

        t_207[k] = f_6 * nsp_82[k]
                   + f_1 * oss0_34[k]
                   - f_2 * oss1_34[k]
                   + f_3 * pc_y[k] * osp_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_y, pc_x, pc_y, nsd0_167, nsp_83, \
                         nsp_105, nsd1_167, oss0_35, oss1_35, osp_104, \
                         osp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_6 * nsp_83[k]
                   + f_3 * pc_y[k] * osp_104[k];

        t_209[k] = pa_y[k] * nsd0_167[k]
                   - f_4 * pc_y[k] * nsd1_167[k];

        t_210[k] = f_12 * nsp_105[k]
                   + f_1 * oss0_35[k]
                   - f_2 * oss1_35[k]
                   + f_3 * pc_x[k] * osp_105[k];

        t_211[k] = f_3 * pc_y[k] * osp_105[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, nsp_83, nsp_107, \
                         oss0_35, oss1_35, osp_106, osp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_12 * nsp_107[k]
                   + f_3 * pc_x[k] * osp_107[k];

        t_213[k] = f_1 * oss0_35[k]
                   - f_2 * oss1_35[k]
                   + f_3 * pc_y[k] * osp_106[k];

        t_214[k] = f_3 * pc_y[k] * osp_107[k];

        t_215[k] = f_11 * nsp_83[k]
                   + f_1 * oss0_35[k]
                   - f_2 * oss1_35[k]
                   + f_3 * pc_z[k] * osp_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, nsp_85, nsp_108, \
                         nsp_109, oss0_36, oss1_36, osp_108, osp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_10 * nsp_108[k]
                   + f_1 * oss0_36[k]
                   - f_2 * oss1_36[k]
                   + f_3 * pc_x[k] * osp_108[k];

        t_217[k] = f_10 * nsp_109[k]
                   + f_3 * pc_x[k] * osp_109[k];

        t_218[k] = f_3 * pc_z[k] * osp_108[k];

        t_219[k] = f_9 * nsp_85[k]
                   + f_1 * oss0_36[k]
                   - f_2 * oss1_36[k]
                   + f_3 * pc_y[k] * osp_109[k];

        t_220[k] = f_3 * pc_z[k] * osp_109[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_z, pc_x, pc_z, nsd0_168, nsp_112, \
                         nsp_113, nsd1_168, oss0_36, oss1_36, osp_110, osp_112, \
                         osp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_1 * oss0_36[k]
                   - f_2 * oss1_36[k]
                   + f_3 * pc_z[k] * osp_110[k];

        t_222[k] = pa_z[k] * nsd0_168[k]
                   - f_4 * pc_z[k] * nsd1_168[k];

        t_223[k] = f_10 * nsp_112[k]
                   + f_3 * pc_x[k] * osp_112[k];

        t_224[k] = f_10 * nsp_113[k]
                   + f_3 * pc_x[k] * osp_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_z, pc_y, pc_z, nsd0_171, nsp_86, nsp_89, \
                         nsd1_171, oss0_37, oss1_37, osp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_z[k] * nsd0_171[k]
                   - f_4 * pc_z[k] * nsd1_171[k];

        t_226[k] = f_11 * nsp_89[k]
                   + f_3 * pc_y[k] * osp_113[k];

        t_227[k] = f_6 * nsp_86[k]
                   + f_1 * oss0_37[k]
                   - f_2 * oss1_37[k]
                   + f_3 * pc_z[k] * osp_113[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, nsp_91, nsp_114, nsp_115, \
                         nsp_116, oss0_38, oss1_38, osp_114, osp_115, \
                         osp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_10 * nsp_114[k]
                   + f_1 * oss0_38[k]
                   - f_2 * oss1_38[k]
                   + f_3 * pc_x[k] * osp_114[k];

        t_229[k] = f_10 * nsp_115[k]
                   + f_3 * pc_x[k] * osp_115[k];

        t_230[k] = f_10 * nsp_116[k]
                   + f_3 * pc_x[k] * osp_116[k];

        t_231[k] = f_13 * nsp_91[k]
                   + f_1 * oss0_38[k]
                   - f_2 * oss1_38[k]
                   + f_3 * pc_y[k] * osp_115[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_y, pc_z, nsp_89, nsp_92, nsp_117, \
                         oss0_38, oss0_39, oss1_38, oss1_39, osp_116, \
                         osp_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_13 * nsp_92[k]
                   + f_3 * pc_y[k] * osp_116[k];

        t_233[k] = f_8 * nsp_89[k]
                   + f_1 * oss0_38[k]
                   - f_2 * oss1_38[k]
                   + f_3 * pc_z[k] * osp_116[k];

        t_234[k] = f_10 * nsp_117[k]
                   + f_1 * oss0_39[k]
                   - f_2 * oss1_39[k]
                   + f_3 * pc_x[k] * osp_117[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, nsp_94, nsp_95, nsp_118, \
                         nsp_119, oss0_39, oss1_39, osp_118, osp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_10 * nsp_118[k]
                   + f_3 * pc_x[k] * osp_118[k];

        t_236[k] = f_10 * nsp_119[k]
                   + f_3 * pc_x[k] * osp_119[k];

        t_237[k] = f_14 * nsp_94[k]
                   + f_1 * oss0_39[k]
                   - f_2 * oss1_39[k]
                   + f_3 * pc_y[k] * osp_118[k];

        t_238[k] = f_14 * nsp_95[k]
                   + f_3 * pc_y[k] * osp_119[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, nsp_92, nsp_120, nsp_121, oss0_39, \
                         oss0_40, oss1_39, oss1_40, osp_119, osp_120, \
                         osp_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * nsp_92[k]
                   + f_1 * oss0_39[k]
                   - f_2 * oss1_39[k]
                   + f_3 * pc_z[k] * osp_119[k];

        t_240[k] = f_10 * nsp_120[k]
                   + f_1 * oss0_40[k]
                   - f_2 * oss1_40[k]
                   + f_3 * pc_x[k] * osp_120[k];

        t_241[k] = f_10 * nsp_121[k]
                   + f_3 * pc_x[k] * osp_121[k];
    }
}

static auto
compute_prim_osd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsd0,
                                                          const size_t nsp, const size_t nsd1,
                                                          const size_t oss0, const size_t oss1,
                                                          const size_t osp, const size_t ncols,
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
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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
    auto *t_359 = buffer.data(target + 359);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsd0_210 = buffer.data(nsd0 + 210);
    const auto *nsd0_215 = buffer.data(nsd0 + 215);
    const auto *nsd0_216 = buffer.data(nsd0 + 216);
    const auto *nsd0_219 = buffer.data(nsd0 + 219);
    const auto *nsd0_264 = buffer.data(nsd0 + 264);
    const auto *nsd0_269 = buffer.data(nsd0 + 269);
    const auto *nsd0_270 = buffer.data(nsd0 + 270);
    const auto *nsd0_330 = buffer.data(nsd0 + 330);
    const auto *nsd0_333 = buffer.data(nsd0 + 333);
    const auto *nsd0_335 = buffer.data(nsd0 + 335);
    const auto *nsd0_339 = buffer.data(nsd0 + 339);
    const auto *nsd0_341 = buffer.data(nsd0 + 341);
    const auto *nsd0_342 = buffer.data(nsd0 + 342);
    const auto *nsd0_345 = buffer.data(nsd0 + 345);
    const auto *nsd0_347 = buffer.data(nsd0 + 347);
    const auto *nsd0_348 = buffer.data(nsd0 + 348);
    const auto *nsd0_351 = buffer.data(nsd0 + 351);
    const auto *nsd0_353 = buffer.data(nsd0 + 353);
    const auto *nsd0_354 = buffer.data(nsd0 + 354);
    const auto *nsd0_357 = buffer.data(nsd0 + 357);
    const auto *nsd0_359 = buffer.data(nsd0 + 359);

    const auto *nsp_95 = buffer.data(nsp + 95);
    const auto *nsp_97 = buffer.data(nsp + 97);
    const auto *nsp_98 = buffer.data(nsp + 98);
    const auto *nsp_100 = buffer.data(nsp + 100);
    const auto *nsp_101 = buffer.data(nsp + 101);
    const auto *nsp_103 = buffer.data(nsp + 103);
    const auto *nsp_104 = buffer.data(nsp + 104);
    const auto *nsp_106 = buffer.data(nsp + 106);
    const auto *nsp_107 = buffer.data(nsp + 107);
    const auto *nsp_109 = buffer.data(nsp + 109);
    const auto *nsp_110 = buffer.data(nsp + 110);
    const auto *nsp_113 = buffer.data(nsp + 113);
    const auto *nsp_115 = buffer.data(nsp + 115);
    const auto *nsp_116 = buffer.data(nsp + 116);
    const auto *nsp_118 = buffer.data(nsp + 118);
    const auto *nsp_119 = buffer.data(nsp + 119);
    const auto *nsp_121 = buffer.data(nsp + 121);
    const auto *nsp_122 = buffer.data(nsp + 122);
    const auto *nsp_123 = buffer.data(nsp + 123);
    const auto *nsp_124 = buffer.data(nsp + 124);
    const auto *nsp_125 = buffer.data(nsp + 125);
    const auto *nsp_126 = buffer.data(nsp + 126);
    const auto *nsp_127 = buffer.data(nsp + 127);
    const auto *nsp_128 = buffer.data(nsp + 128);
    const auto *nsp_130 = buffer.data(nsp + 130);
    const auto *nsp_131 = buffer.data(nsp + 131);
    const auto *nsp_132 = buffer.data(nsp + 132);
    const auto *nsp_133 = buffer.data(nsp + 133);
    const auto *nsp_134 = buffer.data(nsp + 134);
    const auto *nsp_135 = buffer.data(nsp + 135);
    const auto *nsp_136 = buffer.data(nsp + 136);
    const auto *nsp_139 = buffer.data(nsp + 139);
    const auto *nsp_140 = buffer.data(nsp + 140);
    const auto *nsp_141 = buffer.data(nsp + 141);
    const auto *nsp_142 = buffer.data(nsp + 142);
    const auto *nsp_143 = buffer.data(nsp + 143);
    const auto *nsp_144 = buffer.data(nsp + 144);
    const auto *nsp_145 = buffer.data(nsp + 145);
    const auto *nsp_146 = buffer.data(nsp + 146);
    const auto *nsp_147 = buffer.data(nsp + 147);
    const auto *nsp_148 = buffer.data(nsp + 148);
    const auto *nsp_149 = buffer.data(nsp + 149);
    const auto *nsp_150 = buffer.data(nsp + 150);
    const auto *nsp_151 = buffer.data(nsp + 151);
    const auto *nsp_152 = buffer.data(nsp + 152);
    const auto *nsp_153 = buffer.data(nsp + 153);
    const auto *nsp_154 = buffer.data(nsp + 154);
    const auto *nsp_155 = buffer.data(nsp + 155);
    const auto *nsp_156 = buffer.data(nsp + 156);
    const auto *nsp_157 = buffer.data(nsp + 157);
    const auto *nsp_158 = buffer.data(nsp + 158);
    const auto *nsp_160 = buffer.data(nsp + 160);
    const auto *nsp_161 = buffer.data(nsp + 161);
    const auto *nsp_162 = buffer.data(nsp + 162);
    const auto *nsp_164 = buffer.data(nsp + 164);
    const auto *nsp_165 = buffer.data(nsp + 165);
    const auto *nsp_166 = buffer.data(nsp + 166);
    const auto *nsp_169 = buffer.data(nsp + 169);
    const auto *nsp_170 = buffer.data(nsp + 170);
    const auto *nsp_171 = buffer.data(nsp + 171);
    const auto *nsp_172 = buffer.data(nsp + 172);
    const auto *nsp_173 = buffer.data(nsp + 173);
    const auto *nsp_174 = buffer.data(nsp + 174);
    const auto *nsp_175 = buffer.data(nsp + 175);
    const auto *nsp_176 = buffer.data(nsp + 176);
    const auto *nsp_177 = buffer.data(nsp + 177);
    const auto *nsp_178 = buffer.data(nsp + 178);
    const auto *nsp_179 = buffer.data(nsp + 179);

    const auto *nsd1_210 = buffer.data(nsd1 + 210);
    const auto *nsd1_215 = buffer.data(nsd1 + 215);
    const auto *nsd1_216 = buffer.data(nsd1 + 216);
    const auto *nsd1_219 = buffer.data(nsd1 + 219);
    const auto *nsd1_264 = buffer.data(nsd1 + 264);
    const auto *nsd1_269 = buffer.data(nsd1 + 269);
    const auto *nsd1_270 = buffer.data(nsd1 + 270);
    const auto *nsd1_330 = buffer.data(nsd1 + 330);
    const auto *nsd1_333 = buffer.data(nsd1 + 333);
    const auto *nsd1_335 = buffer.data(nsd1 + 335);
    const auto *nsd1_339 = buffer.data(nsd1 + 339);
    const auto *nsd1_341 = buffer.data(nsd1 + 341);
    const auto *nsd1_342 = buffer.data(nsd1 + 342);
    const auto *nsd1_345 = buffer.data(nsd1 + 345);
    const auto *nsd1_347 = buffer.data(nsd1 + 347);
    const auto *nsd1_348 = buffer.data(nsd1 + 348);
    const auto *nsd1_351 = buffer.data(nsd1 + 351);
    const auto *nsd1_353 = buffer.data(nsd1 + 353);
    const auto *nsd1_354 = buffer.data(nsd1 + 354);
    const auto *nsd1_357 = buffer.data(nsd1 + 357);
    const auto *nsd1_359 = buffer.data(nsd1 + 359);

    const auto *oss0_40 = buffer.data(oss0 + 40);
    const auto *oss0_41 = buffer.data(oss0 + 41);
    const auto *oss0_42 = buffer.data(oss0 + 42);
    const auto *oss0_43 = buffer.data(oss0 + 43);
    const auto *oss0_44 = buffer.data(oss0 + 44);
    const auto *oss0_45 = buffer.data(oss0 + 45);
    const auto *oss0_46 = buffer.data(oss0 + 46);
    const auto *oss0_47 = buffer.data(oss0 + 47);
    const auto *oss0_48 = buffer.data(oss0 + 48);
    const auto *oss0_49 = buffer.data(oss0 + 49);
    const auto *oss0_50 = buffer.data(oss0 + 50);
    const auto *oss0_51 = buffer.data(oss0 + 51);
    const auto *oss0_52 = buffer.data(oss0 + 52);
    const auto *oss0_53 = buffer.data(oss0 + 53);
    const auto *oss0_54 = buffer.data(oss0 + 54);

    const auto *oss1_40 = buffer.data(oss1 + 40);
    const auto *oss1_41 = buffer.data(oss1 + 41);
    const auto *oss1_42 = buffer.data(oss1 + 42);
    const auto *oss1_43 = buffer.data(oss1 + 43);
    const auto *oss1_44 = buffer.data(oss1 + 44);
    const auto *oss1_45 = buffer.data(oss1 + 45);
    const auto *oss1_46 = buffer.data(oss1 + 46);
    const auto *oss1_47 = buffer.data(oss1 + 47);
    const auto *oss1_48 = buffer.data(oss1 + 48);
    const auto *oss1_49 = buffer.data(oss1 + 49);
    const auto *oss1_50 = buffer.data(oss1 + 50);
    const auto *oss1_51 = buffer.data(oss1 + 51);
    const auto *oss1_52 = buffer.data(oss1 + 52);
    const auto *oss1_53 = buffer.data(oss1 + 53);
    const auto *oss1_54 = buffer.data(oss1 + 54);

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
    const auto *osp_163 = buffer.data(osp + 163);
    const auto *osp_164 = buffer.data(osp + 164);
    const auto *osp_165 = buffer.data(osp + 165);
    const auto *osp_166 = buffer.data(osp + 166);
    const auto *osp_169 = buffer.data(osp + 169);
    const auto *osp_170 = buffer.data(osp + 170);
    const auto *osp_172 = buffer.data(osp + 172);
    const auto *osp_173 = buffer.data(osp + 173);
    const auto *osp_175 = buffer.data(osp + 175);
    const auto *osp_176 = buffer.data(osp + 176);
    const auto *osp_178 = buffer.data(osp + 178);
    const auto *osp_179 = buffer.data(osp + 179);

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pc_x, pc_y, pc_z, nsp_95, nsp_97, nsp_98, \
                         nsp_122, oss0_40, oss1_40, osp_121, osp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_10 * nsp_122[k]
                   + f_3 * pc_x[k] * osp_122[k];

        t_243[k] = f_12 * nsp_97[k]
                   + f_1 * oss0_40[k]
                   - f_2 * oss1_40[k]
                   + f_3 * pc_y[k] * osp_121[k];

        t_244[k] = f_12 * nsp_98[k]
                   + f_3 * pc_y[k] * osp_122[k];

        t_245[k] = f_12 * nsp_95[k]
                   + f_1 * oss0_40[k]
                   - f_2 * oss1_40[k]
                   + f_3 * pc_z[k] * osp_122[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pc_y, nsp_100, nsp_123, nsp_124, \
                         nsp_125, oss0_41, oss1_41, osp_123, osp_124, \
                         osp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_10 * nsp_123[k]
                   + f_1 * oss0_41[k]
                   - f_2 * oss1_41[k]
                   + f_3 * pc_x[k] * osp_123[k];

        t_247[k] = f_10 * nsp_124[k]
                   + f_3 * pc_x[k] * osp_124[k];

        t_248[k] = f_10 * nsp_125[k]
                   + f_3 * pc_x[k] * osp_125[k];

        t_249[k] = f_10 * nsp_100[k]
                   + f_1 * oss0_41[k]
                   - f_2 * oss1_41[k]
                   + f_3 * pc_y[k] * osp_124[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, nsp_98, nsp_101, nsp_126, \
                         oss0_41, oss0_42, oss1_41, oss1_42, osp_125, \
                         osp_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_10 * nsp_101[k]
                   + f_3 * pc_y[k] * osp_125[k];

        t_251[k] = f_14 * nsp_98[k]
                   + f_1 * oss0_41[k]
                   - f_2 * oss1_41[k]
                   + f_3 * pc_z[k] * osp_125[k];

        t_252[k] = f_10 * nsp_126[k]
                   + f_1 * oss0_42[k]
                   - f_2 * oss1_42[k]
                   + f_3 * pc_x[k] * osp_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, nsp_103, nsp_104, nsp_127, \
                         nsp_128, oss0_42, oss1_42, osp_127, osp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * nsp_127[k]
                   + f_3 * pc_x[k] * osp_127[k];

        t_254[k] = f_10 * nsp_128[k]
                   + f_3 * pc_x[k] * osp_128[k];

        t_255[k] = f_8 * nsp_103[k]
                   + f_1 * oss0_42[k]
                   - f_2 * oss1_42[k]
                   + f_3 * pc_y[k] * osp_127[k];

        t_256[k] = f_8 * nsp_104[k]
                   + f_3 * pc_y[k] * osp_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_y, pc_x, pc_y, pc_z, nsd0_210, nsp_101, \
                         nsp_130, nsd1_210, oss0_42, oss1_42, osp_128, \
                         osp_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_13 * nsp_101[k]
                   + f_1 * oss0_42[k]
                   - f_2 * oss1_42[k]
                   + f_3 * pc_z[k] * osp_128[k];

        t_258[k] = pa_y[k] * nsd0_210[k]
                   - f_4 * pc_y[k] * nsd1_210[k];

        t_259[k] = f_10 * nsp_130[k]
                   + f_3 * pc_x[k] * osp_130[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, nsd0_215, nsp_106, \
                         nsp_107, nsp_131, nsd1_215, oss0_43, oss1_43, osp_130, \
                         osp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * nsp_131[k]
                   + f_3 * pc_x[k] * osp_131[k];

        t_261[k] = f_6 * nsp_106[k]
                   + f_1 * oss0_43[k]
                   - f_2 * oss1_43[k]
                   + f_3 * pc_y[k] * osp_130[k];

        t_262[k] = f_6 * nsp_107[k]
                   + f_3 * pc_y[k] * osp_131[k];

        t_263[k] = pa_y[k] * nsd0_215[k]
                   - f_4 * pc_y[k] * nsd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, nsp_132, nsp_134, \
                         oss0_44, oss1_44, osp_132, osp_133, osp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * nsp_132[k]
                   + f_1 * oss0_44[k]
                   - f_2 * oss1_44[k]
                   + f_3 * pc_x[k] * osp_132[k];

        t_265[k] = f_3 * pc_y[k] * osp_132[k];

        t_266[k] = f_10 * nsp_134[k]
                   + f_3 * pc_x[k] * osp_134[k];

        t_267[k] = f_1 * oss0_44[k]
                   - f_2 * oss1_44[k]
                   + f_3 * pc_y[k] * osp_133[k];

        t_268[k] = f_3 * pc_y[k] * osp_134[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pc_x, pc_z, nsp_107, nsp_135, nsp_136, \
                         oss0_44, oss0_45, oss1_44, oss1_45, osp_134, osp_135, \
                         osp_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_9 * nsp_107[k]
                   + f_1 * oss0_44[k]
                   - f_2 * oss1_44[k]
                   + f_3 * pc_z[k] * osp_134[k];

        t_270[k] = f_8 * nsp_135[k]
                   + f_1 * oss0_45[k]
                   - f_2 * oss1_45[k]
                   + f_3 * pc_x[k] * osp_135[k];

        t_271[k] = f_8 * nsp_136[k]
                   + f_3 * pc_x[k] * osp_136[k];

        t_272[k] = f_3 * pc_z[k] * osp_135[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pa_z, pc_y, pc_z, nsd0_216, nsp_109, \
                         nsd1_216, oss0_45, oss1_45, osp_136, osp_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_7 * nsp_109[k]
                   + f_1 * oss0_45[k]
                   - f_2 * oss1_45[k]
                   + f_3 * pc_y[k] * osp_136[k];

        t_274[k] = f_3 * pc_z[k] * osp_136[k];

        t_275[k] = f_1 * oss0_45[k]
                   - f_2 * oss1_45[k]
                   + f_3 * pc_z[k] * osp_137[k];

        t_276[k] = pa_z[k] * nsd0_216[k]
                   - f_4 * pc_z[k] * nsd1_216[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pa_z, pc_x, pc_y, pc_z, nsd0_219, \
                         nsp_113, nsp_139, nsp_140, nsd1_219, osp_139, \
                         osp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_8 * nsp_139[k]
                   + f_3 * pc_x[k] * osp_139[k];

        t_278[k] = f_8 * nsp_140[k]
                   + f_3 * pc_x[k] * osp_140[k];

        t_279[k] = pa_z[k] * nsd0_219[k]
                   - f_4 * pc_z[k] * nsd1_219[k];

        t_280[k] = f_9 * nsp_113[k]
                   + f_3 * pc_y[k] * osp_140[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_z, nsp_110, nsp_141, nsp_142, oss0_46, \
                         oss0_47, oss1_46, oss1_47, osp_140, osp_141, \
                         osp_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_6 * nsp_110[k]
                   + f_1 * oss0_46[k]
                   - f_2 * oss1_46[k]
                   + f_3 * pc_z[k] * osp_140[k];

        t_282[k] = f_8 * nsp_141[k]
                   + f_1 * oss0_47[k]
                   - f_2 * oss1_47[k]
                   + f_3 * pc_x[k] * osp_141[k];

        t_283[k] = f_8 * nsp_142[k]
                   + f_3 * pc_x[k] * osp_142[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_y, pc_z, nsp_113, nsp_115, \
                         nsp_116, nsp_143, oss0_47, oss1_47, osp_142, \
                         osp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_8 * nsp_143[k]
                   + f_3 * pc_x[k] * osp_143[k];

        t_285[k] = f_11 * nsp_115[k]
                   + f_1 * oss0_47[k]
                   - f_2 * oss1_47[k]
                   + f_3 * pc_y[k] * osp_142[k];

        t_286[k] = f_11 * nsp_116[k]
                   + f_3 * pc_y[k] * osp_143[k];

        t_287[k] = f_8 * nsp_113[k]
                   + f_1 * oss0_47[k]
                   - f_2 * oss1_47[k]
                   + f_3 * pc_z[k] * osp_143[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, nsp_118, nsp_144, nsp_145, \
                         nsp_146, oss0_48, oss1_48, osp_144, osp_145, \
                         osp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_8 * nsp_144[k]
                   + f_1 * oss0_48[k]
                   - f_2 * oss1_48[k]
                   + f_3 * pc_x[k] * osp_144[k];

        t_289[k] = f_8 * nsp_145[k]
                   + f_3 * pc_x[k] * osp_145[k];

        t_290[k] = f_8 * nsp_146[k]
                   + f_3 * pc_x[k] * osp_146[k];

        t_291[k] = f_13 * nsp_118[k]
                   + f_1 * oss0_48[k]
                   - f_2 * oss1_48[k]
                   + f_3 * pc_y[k] * osp_145[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, pc_x, pc_y, pc_z, nsp_116, nsp_119, nsp_147, \
                         oss0_48, oss0_49, oss1_48, oss1_49, osp_146, \
                         osp_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_13 * nsp_119[k]
                   + f_3 * pc_y[k] * osp_146[k];

        t_293[k] = f_10 * nsp_116[k]
                   + f_1 * oss0_48[k]
                   - f_2 * oss1_48[k]
                   + f_3 * pc_z[k] * osp_146[k];

        t_294[k] = f_8 * nsp_147[k]
                   + f_1 * oss0_49[k]
                   - f_2 * oss1_49[k]
                   + f_3 * pc_x[k] * osp_147[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_x, pc_y, nsp_121, nsp_122, nsp_148, \
                         nsp_149, oss0_49, oss1_49, osp_148, osp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_8 * nsp_148[k]
                   + f_3 * pc_x[k] * osp_148[k];

        t_296[k] = f_8 * nsp_149[k]
                   + f_3 * pc_x[k] * osp_149[k];

        t_297[k] = f_14 * nsp_121[k]
                   + f_1 * oss0_49[k]
                   - f_2 * oss1_49[k]
                   + f_3 * pc_y[k] * osp_148[k];

        t_298[k] = f_14 * nsp_122[k]
                   + f_3 * pc_y[k] * osp_149[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pc_x, pc_z, nsp_119, nsp_150, nsp_151, oss0_49, \
                         oss0_50, oss1_49, oss1_50, osp_149, osp_150, \
                         osp_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_12 * nsp_119[k]
                   + f_1 * oss0_49[k]
                   - f_2 * oss1_49[k]
                   + f_3 * pc_z[k] * osp_149[k];

        t_300[k] = f_8 * nsp_150[k]
                   + f_1 * oss0_50[k]
                   - f_2 * oss1_50[k]
                   + f_3 * pc_x[k] * osp_150[k];

        t_301[k] = f_8 * nsp_151[k]
                   + f_3 * pc_x[k] * osp_151[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_y, pc_z, nsp_122, nsp_124, \
                         nsp_125, nsp_152, oss0_50, oss1_50, osp_151, \
                         osp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_8 * nsp_152[k]
                   + f_3 * pc_x[k] * osp_152[k];

        t_303[k] = f_12 * nsp_124[k]
                   + f_1 * oss0_50[k]
                   - f_2 * oss1_50[k]
                   + f_3 * pc_y[k] * osp_151[k];

        t_304[k] = f_12 * nsp_125[k]
                   + f_3 * pc_y[k] * osp_152[k];

        t_305[k] = f_14 * nsp_122[k]
                   + f_1 * oss0_50[k]
                   - f_2 * oss1_50[k]
                   + f_3 * pc_z[k] * osp_152[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_x, pc_y, nsp_127, nsp_153, nsp_154, \
                         nsp_155, oss0_51, oss1_51, osp_153, osp_154, \
                         osp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_8 * nsp_153[k]
                   + f_1 * oss0_51[k]
                   - f_2 * oss1_51[k]
                   + f_3 * pc_x[k] * osp_153[k];

        t_307[k] = f_8 * nsp_154[k]
                   + f_3 * pc_x[k] * osp_154[k];

        t_308[k] = f_8 * nsp_155[k]
                   + f_3 * pc_x[k] * osp_155[k];

        t_309[k] = f_10 * nsp_127[k]
                   + f_1 * oss0_51[k]
                   - f_2 * oss1_51[k]
                   + f_3 * pc_y[k] * osp_154[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pc_x, pc_y, pc_z, nsp_125, nsp_128, nsp_156, \
                         oss0_51, oss0_52, oss1_51, oss1_52, osp_155, \
                         osp_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_10 * nsp_128[k]
                   + f_3 * pc_y[k] * osp_155[k];

        t_311[k] = f_13 * nsp_125[k]
                   + f_1 * oss0_51[k]
                   - f_2 * oss1_51[k]
                   + f_3 * pc_z[k] * osp_155[k];

        t_312[k] = f_8 * nsp_156[k]
                   + f_1 * oss0_52[k]
                   - f_2 * oss1_52[k]
                   + f_3 * pc_x[k] * osp_156[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, nsp_130, nsp_131, nsp_157, \
                         nsp_158, oss0_52, oss1_52, osp_157, osp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_8 * nsp_157[k]
                   + f_3 * pc_x[k] * osp_157[k];

        t_314[k] = f_8 * nsp_158[k]
                   + f_3 * pc_x[k] * osp_158[k];

        t_315[k] = f_8 * nsp_130[k]
                   + f_1 * oss0_52[k]
                   - f_2 * oss1_52[k]
                   + f_3 * pc_y[k] * osp_157[k];

        t_316[k] = f_8 * nsp_131[k]
                   + f_3 * pc_y[k] * osp_158[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pa_y, pc_x, pc_y, pc_z, nsd0_264, nsp_128, \
                         nsp_160, nsd1_264, oss0_52, oss1_52, osp_158, \
                         osp_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_11 * nsp_128[k]
                   + f_1 * oss0_52[k]
                   - f_2 * oss1_52[k]
                   + f_3 * pc_z[k] * osp_158[k];

        t_318[k] = pa_y[k] * nsd0_264[k]
                   - f_4 * pc_y[k] * nsd1_264[k];

        t_319[k] = f_8 * nsp_160[k]
                   + f_3 * pc_x[k] * osp_160[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_y, pc_x, pc_y, nsd0_269, nsp_133, \
                         nsp_134, nsp_161, nsd1_269, oss0_53, oss1_53, osp_160, \
                         osp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_8 * nsp_161[k]
                   + f_3 * pc_x[k] * osp_161[k];

        t_321[k] = f_6 * nsp_133[k]
                   + f_1 * oss0_53[k]
                   - f_2 * oss1_53[k]
                   + f_3 * pc_y[k] * osp_160[k];

        t_322[k] = f_6 * nsp_134[k]
                   + f_3 * pc_y[k] * osp_161[k];

        t_323[k] = pa_y[k] * nsd0_269[k]
                   - f_4 * pc_y[k] * nsd1_269[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, nsp_162, nsp_164, \
                         oss0_54, oss1_54, osp_162, osp_163, osp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_8 * nsp_162[k]
                   + f_1 * oss0_54[k]
                   - f_2 * oss1_54[k]
                   + f_3 * pc_x[k] * osp_162[k];

        t_325[k] = f_3 * pc_y[k] * osp_162[k];

        t_326[k] = f_8 * nsp_164[k]
                   + f_3 * pc_x[k] * osp_164[k];

        t_327[k] = f_1 * oss0_54[k]
                   - f_2 * oss1_54[k]
                   + f_3 * pc_y[k] * osp_163[k];

        t_328[k] = f_3 * pc_y[k] * osp_164[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_x, pc_x, pc_z, nsd0_330, nsp_134, nsp_165, \
                         nsp_166, nsd1_330, oss0_54, oss1_54, osp_164, \
                         osp_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_7 * nsp_134[k]
                   + f_1 * oss0_54[k]
                   - f_2 * oss1_54[k]
                   + f_3 * pc_z[k] * osp_164[k];

        t_330[k] = pa_x[k] * nsd0_330[k]
                   + f_8 * nsp_165[k]
                   - f_4 * pc_x[k] * nsd1_330[k];

        t_331[k] = f_6 * nsp_166[k]
                   + f_3 * pc_x[k] * osp_166[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pa_x, pc_x, pc_z, nsd0_333, nsd0_335, \
                         nsd1_333, nsd1_335, osp_165, osp_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_3 * pc_z[k] * osp_165[k];

        t_333[k] = pa_x[k] * nsd0_333[k]
                   - f_4 * pc_x[k] * nsd1_333[k];

        t_334[k] = f_3 * pc_z[k] * osp_166[k];

        t_335[k] = pa_x[k] * nsd0_335[k]
                   - f_4 * pc_x[k] * nsd1_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_x, pa_z, pc_x, pc_z, nsd0_270, \
                         nsd0_339, nsp_169, nsp_170, nsd1_270, nsd1_339, osp_169, \
                         osp_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_z[k] * nsd0_270[k]
                   - f_4 * pc_z[k] * nsd1_270[k];

        t_337[k] = f_6 * nsp_169[k]
                   + f_3 * pc_x[k] * osp_169[k];

        t_338[k] = f_6 * nsp_170[k]
                   + f_3 * pc_x[k] * osp_170[k];

        t_339[k] = pa_x[k] * nsd0_339[k]
                   - f_4 * pc_x[k] * nsd1_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_x, pc_x, pc_y, nsd0_341, nsd0_342, \
                         nsp_140, nsp_171, nsp_172, nsd1_341, nsd1_342, osp_170, \
                         osp_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_7 * nsp_140[k]
                   + f_3 * pc_y[k] * osp_170[k];

        t_341[k] = pa_x[k] * nsd0_341[k]
                   - f_4 * pc_x[k] * nsd1_341[k];

        t_342[k] = pa_x[k] * nsd0_342[k]
                   + f_8 * nsp_171[k]
                   - f_4 * pc_x[k] * nsd1_342[k];

        t_343[k] = f_6 * nsp_172[k]
                   + f_3 * pc_x[k] * osp_172[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_x, pc_x, pc_y, nsd0_345, nsd0_347, \
                         nsp_143, nsp_173, nsd1_345, nsd1_347, \
                         osp_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_6 * nsp_173[k]
                   + f_3 * pc_x[k] * osp_173[k];

        t_345[k] = pa_x[k] * nsd0_345[k]
                   - f_4 * pc_x[k] * nsd1_345[k];

        t_346[k] = f_9 * nsp_143[k]
                   + f_3 * pc_y[k] * osp_173[k];

        t_347[k] = pa_x[k] * nsd0_347[k]
                   - f_4 * pc_x[k] * nsd1_347[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_x, pc_x, nsd0_348, nsd0_351, nsp_174, \
                         nsp_175, nsp_176, nsd1_348, nsd1_351, osp_175, \
                         osp_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = pa_x[k] * nsd0_348[k]
                   + f_8 * nsp_174[k]
                   - f_4 * pc_x[k] * nsd1_348[k];

        t_349[k] = f_6 * nsp_175[k]
                   + f_3 * pc_x[k] * osp_175[k];

        t_350[k] = f_6 * nsp_176[k]
                   + f_3 * pc_x[k] * osp_176[k];

        t_351[k] = pa_x[k] * nsd0_351[k]
                   - f_4 * pc_x[k] * nsd1_351[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_x, pc_x, pc_y, nsd0_353, nsd0_354, \
                         nsp_146, nsp_177, nsp_178, nsd1_353, nsd1_354, osp_176, \
                         osp_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_11 * nsp_146[k]
                   + f_3 * pc_y[k] * osp_176[k];

        t_353[k] = pa_x[k] * nsd0_353[k]
                   - f_4 * pc_x[k] * nsd1_353[k];

        t_354[k] = pa_x[k] * nsd0_354[k]
                   + f_8 * nsp_177[k]
                   - f_4 * pc_x[k] * nsd1_354[k];

        t_355[k] = f_6 * nsp_178[k]
                   + f_3 * pc_x[k] * osp_178[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pc_x, pc_y, nsd0_357, nsd0_359, \
                         nsp_149, nsp_179, nsd1_357, nsd1_359, \
                         osp_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_6 * nsp_179[k]
                   + f_3 * pc_x[k] * osp_179[k];

        t_357[k] = pa_x[k] * nsd0_357[k]
                   - f_4 * pc_x[k] * nsd1_357[k];

        t_358[k] = f_13 * nsp_149[k]
                   + f_3 * pc_y[k] * osp_179[k];

        t_359[k] = pa_x[k] * nsd0_359[k]
                   - f_4 * pc_x[k] * nsd1_359[k];
    }
}

static auto
compute_prim_osd_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsd0,
                                                          const size_t nsp, const size_t nsd1,
                                                          const size_t oss0, const size_t oss1,
                                                          const size_t osp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 5.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsd0_324 = buffer.data(nsd0 + 324);
    const auto *nsd0_330 = buffer.data(nsd0 + 330);
    const auto *nsd0_333 = buffer.data(nsd0 + 333);
    const auto *nsd0_360 = buffer.data(nsd0 + 360);
    const auto *nsd0_363 = buffer.data(nsd0 + 363);
    const auto *nsd0_365 = buffer.data(nsd0 + 365);
    const auto *nsd0_366 = buffer.data(nsd0 + 366);
    const auto *nsd0_369 = buffer.data(nsd0 + 369);
    const auto *nsd0_371 = buffer.data(nsd0 + 371);
    const auto *nsd0_372 = buffer.data(nsd0 + 372);
    const auto *nsd0_375 = buffer.data(nsd0 + 375);
    const auto *nsd0_377 = buffer.data(nsd0 + 377);
    const auto *nsd0_378 = buffer.data(nsd0 + 378);
    const auto *nsd0_381 = buffer.data(nsd0 + 381);
    const auto *nsd0_383 = buffer.data(nsd0 + 383);
    const auto *nsd0_387 = buffer.data(nsd0 + 387);
    const auto *nsd0_389 = buffer.data(nsd0 + 389);
    const auto *nsd0_390 = buffer.data(nsd0 + 390);
    const auto *nsd0_393 = buffer.data(nsd0 + 393);
    const auto *nsd0_395 = buffer.data(nsd0 + 395);

    const auto *nsp_152 = buffer.data(nsp + 152);
    const auto *nsp_155 = buffer.data(nsp + 155);
    const auto *nsp_158 = buffer.data(nsp + 158);
    const auto *nsp_161 = buffer.data(nsp + 161);
    const auto *nsp_164 = buffer.data(nsp + 164);
    const auto *nsp_166 = buffer.data(nsp + 166);
    const auto *nsp_167 = buffer.data(nsp + 167);
    const auto *nsp_170 = buffer.data(nsp + 170);
    const auto *nsp_172 = buffer.data(nsp + 172);
    const auto *nsp_173 = buffer.data(nsp + 173);
    const auto *nsp_175 = buffer.data(nsp + 175);
    const auto *nsp_176 = buffer.data(nsp + 176);
    const auto *nsp_178 = buffer.data(nsp + 178);
    const auto *nsp_179 = buffer.data(nsp + 179);
    const auto *nsp_180 = buffer.data(nsp + 180);
    const auto *nsp_181 = buffer.data(nsp + 181);
    const auto *nsp_182 = buffer.data(nsp + 182);
    const auto *nsp_183 = buffer.data(nsp + 183);
    const auto *nsp_184 = buffer.data(nsp + 184);
    const auto *nsp_185 = buffer.data(nsp + 185);
    const auto *nsp_186 = buffer.data(nsp + 186);
    const auto *nsp_187 = buffer.data(nsp + 187);
    const auto *nsp_188 = buffer.data(nsp + 188);
    const auto *nsp_189 = buffer.data(nsp + 189);
    const auto *nsp_190 = buffer.data(nsp + 190);
    const auto *nsp_191 = buffer.data(nsp + 191);
    const auto *nsp_193 = buffer.data(nsp + 193);
    const auto *nsp_194 = buffer.data(nsp + 194);
    const auto *nsp_195 = buffer.data(nsp + 195);
    const auto *nsp_196 = buffer.data(nsp + 196);
    const auto *nsp_197 = buffer.data(nsp + 197);

    const auto *nsd1_324 = buffer.data(nsd1 + 324);
    const auto *nsd1_330 = buffer.data(nsd1 + 330);
    const auto *nsd1_333 = buffer.data(nsd1 + 333);
    const auto *nsd1_360 = buffer.data(nsd1 + 360);
    const auto *nsd1_363 = buffer.data(nsd1 + 363);
    const auto *nsd1_365 = buffer.data(nsd1 + 365);
    const auto *nsd1_366 = buffer.data(nsd1 + 366);
    const auto *nsd1_369 = buffer.data(nsd1 + 369);
    const auto *nsd1_371 = buffer.data(nsd1 + 371);
    const auto *nsd1_372 = buffer.data(nsd1 + 372);
    const auto *nsd1_375 = buffer.data(nsd1 + 375);
    const auto *nsd1_377 = buffer.data(nsd1 + 377);
    const auto *nsd1_378 = buffer.data(nsd1 + 378);
    const auto *nsd1_381 = buffer.data(nsd1 + 381);
    const auto *nsd1_383 = buffer.data(nsd1 + 383);
    const auto *nsd1_387 = buffer.data(nsd1 + 387);
    const auto *nsd1_389 = buffer.data(nsd1 + 389);
    const auto *nsd1_390 = buffer.data(nsd1 + 390);
    const auto *nsd1_393 = buffer.data(nsd1 + 393);
    const auto *nsd1_395 = buffer.data(nsd1 + 395);

    const auto *oss0_66 = buffer.data(oss0 + 66);
    const auto *oss0_67 = buffer.data(oss0 + 67);
    const auto *oss0_68 = buffer.data(oss0 + 68);
    const auto *oss0_69 = buffer.data(oss0 + 69);
    const auto *oss0_70 = buffer.data(oss0 + 70);
    const auto *oss0_71 = buffer.data(oss0 + 71);
    const auto *oss0_72 = buffer.data(oss0 + 72);
    const auto *oss0_73 = buffer.data(oss0 + 73);
    const auto *oss0_74 = buffer.data(oss0 + 74);
    const auto *oss0_75 = buffer.data(oss0 + 75);
    const auto *oss0_77 = buffer.data(oss0 + 77);

    const auto *oss1_66 = buffer.data(oss1 + 66);
    const auto *oss1_67 = buffer.data(oss1 + 67);
    const auto *oss1_68 = buffer.data(oss1 + 68);
    const auto *oss1_69 = buffer.data(oss1 + 69);
    const auto *oss1_70 = buffer.data(oss1 + 70);
    const auto *oss1_71 = buffer.data(oss1 + 71);
    const auto *oss1_72 = buffer.data(oss1 + 72);
    const auto *oss1_73 = buffer.data(oss1 + 73);
    const auto *oss1_74 = buffer.data(oss1 + 74);
    const auto *oss1_75 = buffer.data(oss1 + 75);
    const auto *oss1_77 = buffer.data(oss1 + 77);

    const auto *osp_181 = buffer.data(osp + 181);
    const auto *osp_182 = buffer.data(osp + 182);
    const auto *osp_184 = buffer.data(osp + 184);
    const auto *osp_185 = buffer.data(osp + 185);
    const auto *osp_187 = buffer.data(osp + 187);
    const auto *osp_188 = buffer.data(osp + 188);
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
    const auto *osp_232 = buffer.data(osp + 232);
    const auto *osp_233 = buffer.data(osp + 233);

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pa_x, pc_x, nsd0_360, nsd0_363, nsp_180, \
                         nsp_181, nsp_182, nsd1_360, nsd1_363, osp_181, \
                         osp_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pa_x[k] * nsd0_360[k]
                   + f_8 * nsp_180[k]
                   - f_4 * pc_x[k] * nsd1_360[k];

        t_361[k] = f_6 * nsp_181[k]
                   + f_3 * pc_x[k] * osp_181[k];

        t_362[k] = f_6 * nsp_182[k]
                   + f_3 * pc_x[k] * osp_182[k];

        t_363[k] = pa_x[k] * nsd0_363[k]
                   - f_4 * pc_x[k] * nsd1_363[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pa_x, pc_x, pc_y, nsd0_365, nsd0_366, \
                         nsp_152, nsp_183, nsp_184, nsd1_365, nsd1_366, osp_182, \
                         osp_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_14 * nsp_152[k]
                   + f_3 * pc_y[k] * osp_182[k];

        t_365[k] = pa_x[k] * nsd0_365[k]
                   - f_4 * pc_x[k] * nsd1_365[k];

        t_366[k] = pa_x[k] * nsd0_366[k]
                   + f_8 * nsp_183[k]
                   - f_4 * pc_x[k] * nsd1_366[k];

        t_367[k] = f_6 * nsp_184[k]
                   + f_3 * pc_x[k] * osp_184[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pa_x, pc_x, pc_y, nsd0_369, nsd0_371, \
                         nsp_155, nsp_185, nsd1_369, nsd1_371, \
                         osp_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_6 * nsp_185[k]
                   + f_3 * pc_x[k] * osp_185[k];

        t_369[k] = pa_x[k] * nsd0_369[k]
                   - f_4 * pc_x[k] * nsd1_369[k];

        t_370[k] = f_12 * nsp_155[k]
                   + f_3 * pc_y[k] * osp_185[k];

        t_371[k] = pa_x[k] * nsd0_371[k]
                   - f_4 * pc_x[k] * nsd1_371[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pa_x, pc_x, nsd0_372, nsd0_375, nsp_186, \
                         nsp_187, nsp_188, nsd1_372, nsd1_375, osp_187, \
                         osp_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = pa_x[k] * nsd0_372[k]
                   + f_8 * nsp_186[k]
                   - f_4 * pc_x[k] * nsd1_372[k];

        t_373[k] = f_6 * nsp_187[k]
                   + f_3 * pc_x[k] * osp_187[k];

        t_374[k] = f_6 * nsp_188[k]
                   + f_3 * pc_x[k] * osp_188[k];

        t_375[k] = pa_x[k] * nsd0_375[k]
                   - f_4 * pc_x[k] * nsd1_375[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_x, pc_x, pc_y, nsd0_377, nsd0_378, \
                         nsp_158, nsp_189, nsp_190, nsd1_377, nsd1_378, osp_188, \
                         osp_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_10 * nsp_158[k]
                   + f_3 * pc_y[k] * osp_188[k];

        t_377[k] = pa_x[k] * nsd0_377[k]
                   - f_4 * pc_x[k] * nsd1_377[k];

        t_378[k] = pa_x[k] * nsd0_378[k]
                   + f_8 * nsp_189[k]
                   - f_4 * pc_x[k] * nsd1_378[k];

        t_379[k] = f_6 * nsp_190[k]
                   + f_3 * pc_x[k] * osp_190[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pa_x, pc_x, pc_y, nsd0_381, nsd0_383, \
                         nsp_161, nsp_191, nsd1_381, nsd1_383, \
                         osp_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_6 * nsp_191[k]
                   + f_3 * pc_x[k] * osp_191[k];

        t_381[k] = pa_x[k] * nsd0_381[k]
                   - f_4 * pc_x[k] * nsd1_381[k];

        t_382[k] = f_8 * nsp_161[k]
                   + f_3 * pc_y[k] * osp_191[k];

        t_383[k] = pa_x[k] * nsd0_383[k]
                   - f_4 * pc_x[k] * nsd1_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pa_x, pa_y, pc_x, pc_y, nsd0_324, \
                         nsd0_387, nsp_193, nsp_194, nsd1_324, nsd1_387, osp_193, \
                         osp_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_y[k] * nsd0_324[k]
                   - f_4 * pc_y[k] * nsd1_324[k];

        t_385[k] = f_6 * nsp_193[k]
                   + f_3 * pc_x[k] * osp_193[k];

        t_386[k] = f_6 * nsp_194[k]
                   + f_3 * pc_x[k] * osp_194[k];

        t_387[k] = pa_x[k] * nsd0_387[k]
                   - f_4 * pc_x[k] * nsd1_387[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pa_x, pc_x, pc_y, nsd0_389, nsd0_390, \
                         nsp_164, nsp_195, nsd1_389, nsd1_390, osp_194, \
                         osp_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_6 * nsp_164[k]
                   + f_3 * pc_y[k] * osp_194[k];

        t_389[k] = pa_x[k] * nsd0_389[k]
                   - f_4 * pc_x[k] * nsd1_389[k];

        t_390[k] = pa_x[k] * nsd0_390[k]
                   + f_8 * nsp_195[k]
                   - f_4 * pc_x[k] * nsd1_390[k];

        t_391[k] = f_3 * pc_y[k] * osp_195[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pa_x, pc_x, pc_y, nsd0_393, nsd0_395, \
                         nsp_197, nsd1_393, nsd1_395, osp_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_6 * nsp_197[k]
                   + f_3 * pc_x[k] * osp_197[k];

        t_393[k] = pa_x[k] * nsd0_393[k]
                   - f_4 * pc_x[k] * nsd1_393[k];

        t_394[k] = f_3 * pc_y[k] * osp_197[k];

        t_395[k] = pa_x[k] * nsd0_395[k]
                   - f_4 * pc_x[k] * nsd1_395[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, t_401, pc_x, pc_y, pc_z, nsp_166, \
                         oss0_66, oss1_66, osp_198, osp_199, osp_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_1 * oss0_66[k]
                   - f_2 * oss1_66[k]
                   + f_3 * pc_x[k] * osp_198[k];

        t_397[k] = f_3 * pc_x[k] * osp_199[k];

        t_398[k] = f_3 * pc_x[k] * osp_200[k];

        t_399[k] = f_0 * nsp_166[k]
                   + f_1 * oss0_66[k]
                   - f_2 * oss1_66[k]
                   + f_3 * pc_y[k] * osp_199[k];

        t_400[k] = f_3 * pc_z[k] * osp_199[k];

        t_401[k] = f_1 * oss0_66[k]
                   - f_2 * oss1_66[k]
                   + f_3 * pc_z[k] * osp_200[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_z, pc_x, pc_y, pc_z, nsd0_330, \
                         nsd0_333, nsp_170, nsd1_330, nsd1_333, osp_202, \
                         osp_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_z[k] * nsd0_330[k]
                   - f_4 * pc_z[k] * nsd1_330[k];

        t_403[k] = f_3 * pc_x[k] * osp_202[k];

        t_404[k] = f_3 * pc_x[k] * osp_203[k];

        t_405[k] = pa_z[k] * nsd0_333[k]
                   - f_4 * pc_z[k] * nsd1_333[k];

        t_406[k] = f_5 * nsp_170[k]
                   + f_3 * pc_y[k] * osp_203[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pc_x, pc_z, nsp_167, oss0_67, oss0_68, \
                         oss1_67, oss1_68, osp_203, osp_204, osp_205, \
                         osp_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_6 * nsp_167[k]
                   + f_1 * oss0_67[k]
                   - f_2 * oss1_67[k]
                   + f_3 * pc_z[k] * osp_203[k];

        t_408[k] = f_1 * oss0_68[k]
                   - f_2 * oss1_68[k]
                   + f_3 * pc_x[k] * osp_204[k];

        t_409[k] = f_3 * pc_x[k] * osp_205[k];

        t_410[k] = f_3 * pc_x[k] * osp_206[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, pc_z, nsp_170, nsp_172, nsp_173, oss0_68, \
                         oss1_68, osp_205, osp_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_7 * nsp_172[k]
                   + f_1 * oss0_68[k]
                   - f_2 * oss1_68[k]
                   + f_3 * pc_y[k] * osp_205[k];

        t_412[k] = f_7 * nsp_173[k]
                   + f_3 * pc_y[k] * osp_206[k];

        t_413[k] = f_8 * nsp_170[k]
                   + f_1 * oss0_68[k]
                   - f_2 * oss1_68[k]
                   + f_3 * pc_z[k] * osp_206[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, pc_x, pc_y, nsp_175, nsp_176, \
                         oss0_69, oss1_69, osp_207, osp_208, osp_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_1 * oss0_69[k]
                   - f_2 * oss1_69[k]
                   + f_3 * pc_x[k] * osp_207[k];

        t_415[k] = f_3 * pc_x[k] * osp_208[k];

        t_416[k] = f_3 * pc_x[k] * osp_209[k];

        t_417[k] = f_9 * nsp_175[k]
                   + f_1 * oss0_69[k]
                   - f_2 * oss1_69[k]
                   + f_3 * pc_y[k] * osp_208[k];

        t_418[k] = f_9 * nsp_176[k]
                   + f_3 * pc_y[k] * osp_209[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pc_x, pc_z, nsp_173, oss0_69, oss0_70, \
                         oss1_69, oss1_70, osp_209, osp_210, osp_211, \
                         osp_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_10 * nsp_173[k]
                   + f_1 * oss0_69[k]
                   - f_2 * oss1_69[k]
                   + f_3 * pc_z[k] * osp_209[k];

        t_420[k] = f_1 * oss0_70[k]
                   - f_2 * oss1_70[k]
                   + f_3 * pc_x[k] * osp_210[k];

        t_421[k] = f_3 * pc_x[k] * osp_211[k];

        t_422[k] = f_3 * pc_x[k] * osp_212[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_y, pc_z, nsp_176, nsp_178, nsp_179, oss0_70, \
                         oss1_70, osp_211, osp_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_11 * nsp_178[k]
                   + f_1 * oss0_70[k]
                   - f_2 * oss1_70[k]
                   + f_3 * pc_y[k] * osp_211[k];

        t_424[k] = f_11 * nsp_179[k]
                   + f_3 * pc_y[k] * osp_212[k];

        t_425[k] = f_12 * nsp_176[k]
                   + f_1 * oss0_70[k]
                   - f_2 * oss1_70[k]
                   + f_3 * pc_z[k] * osp_212[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, t_430, pc_x, pc_y, nsp_181, nsp_182, \
                         oss0_71, oss1_71, osp_213, osp_214, osp_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_1 * oss0_71[k]
                   - f_2 * oss1_71[k]
                   + f_3 * pc_x[k] * osp_213[k];

        t_427[k] = f_3 * pc_x[k] * osp_214[k];

        t_428[k] = f_3 * pc_x[k] * osp_215[k];

        t_429[k] = f_13 * nsp_181[k]
                   + f_1 * oss0_71[k]
                   - f_2 * oss1_71[k]
                   + f_3 * pc_y[k] * osp_214[k];

        t_430[k] = f_13 * nsp_182[k]
                   + f_3 * pc_y[k] * osp_215[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, pc_z, nsp_179, oss0_71, oss0_72, \
                         oss1_71, oss1_72, osp_215, osp_216, osp_217, \
                         osp_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_14 * nsp_179[k]
                   + f_1 * oss0_71[k]
                   - f_2 * oss1_71[k]
                   + f_3 * pc_z[k] * osp_215[k];

        t_432[k] = f_1 * oss0_72[k]
                   - f_2 * oss1_72[k]
                   + f_3 * pc_x[k] * osp_216[k];

        t_433[k] = f_3 * pc_x[k] * osp_217[k];

        t_434[k] = f_3 * pc_x[k] * osp_218[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pc_y, pc_z, nsp_182, nsp_184, nsp_185, oss0_72, \
                         oss1_72, osp_217, osp_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_14 * nsp_184[k]
                   + f_1 * oss0_72[k]
                   - f_2 * oss1_72[k]
                   + f_3 * pc_y[k] * osp_217[k];

        t_436[k] = f_14 * nsp_185[k]
                   + f_3 * pc_y[k] * osp_218[k];

        t_437[k] = f_13 * nsp_182[k]
                   + f_1 * oss0_72[k]
                   - f_2 * oss1_72[k]
                   + f_3 * pc_z[k] * osp_218[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pc_x, pc_y, nsp_187, nsp_188, \
                         oss0_73, oss1_73, osp_219, osp_220, osp_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_1 * oss0_73[k]
                   - f_2 * oss1_73[k]
                   + f_3 * pc_x[k] * osp_219[k];

        t_439[k] = f_3 * pc_x[k] * osp_220[k];

        t_440[k] = f_3 * pc_x[k] * osp_221[k];

        t_441[k] = f_12 * nsp_187[k]
                   + f_1 * oss0_73[k]
                   - f_2 * oss1_73[k]
                   + f_3 * pc_y[k] * osp_220[k];

        t_442[k] = f_12 * nsp_188[k]
                   + f_3 * pc_y[k] * osp_221[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pc_x, pc_z, nsp_185, oss0_73, oss0_74, \
                         oss1_73, oss1_74, osp_221, osp_222, osp_223, \
                         osp_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_11 * nsp_185[k]
                   + f_1 * oss0_73[k]
                   - f_2 * oss1_73[k]
                   + f_3 * pc_z[k] * osp_221[k];

        t_444[k] = f_1 * oss0_74[k]
                   - f_2 * oss1_74[k]
                   + f_3 * pc_x[k] * osp_222[k];

        t_445[k] = f_3 * pc_x[k] * osp_223[k];

        t_446[k] = f_3 * pc_x[k] * osp_224[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, nsp_188, nsp_190, nsp_191, oss0_74, \
                         oss1_74, osp_223, osp_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_10 * nsp_190[k]
                   + f_1 * oss0_74[k]
                   - f_2 * oss1_74[k]
                   + f_3 * pc_y[k] * osp_223[k];

        t_448[k] = f_10 * nsp_191[k]
                   + f_3 * pc_y[k] * osp_224[k];

        t_449[k] = f_9 * nsp_188[k]
                   + f_1 * oss0_74[k]
                   - f_2 * oss1_74[k]
                   + f_3 * pc_z[k] * osp_224[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, pc_x, pc_y, nsp_193, nsp_194, \
                         oss0_75, oss1_75, osp_225, osp_226, osp_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * oss0_75[k]
                   - f_2 * oss1_75[k]
                   + f_3 * pc_x[k] * osp_225[k];

        t_451[k] = f_3 * pc_x[k] * osp_226[k];

        t_452[k] = f_3 * pc_x[k] * osp_227[k];

        t_453[k] = f_8 * nsp_193[k]
                   + f_1 * oss0_75[k]
                   - f_2 * oss1_75[k]
                   + f_3 * pc_y[k] * osp_226[k];

        t_454[k] = f_8 * nsp_194[k]
                   + f_3 * pc_y[k] * osp_227[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pa_y, pc_x, pc_y, pc_z, nsd0_390, \
                         nsp_191, nsd1_390, oss0_75, oss1_75, osp_227, osp_229, \
                         osp_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_7 * nsp_191[k]
                   + f_1 * oss0_75[k]
                   - f_2 * oss1_75[k]
                   + f_3 * pc_z[k] * osp_227[k];

        t_456[k] = pa_y[k] * nsd0_390[k]
                   - f_4 * pc_y[k] * nsd1_390[k];

        t_457[k] = f_3 * pc_x[k] * osp_229[k];

        t_458[k] = f_3 * pc_x[k] * osp_230[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_y, pc_y, nsd0_393, nsd0_395, nsp_196, \
                         nsp_197, nsd1_393, nsd1_395, osp_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = pa_y[k] * nsd0_393[k]
                   + f_8 * nsp_196[k]
                   - f_4 * pc_y[k] * nsd1_393[k];

        t_460[k] = f_6 * nsp_197[k]
                   + f_3 * pc_y[k] * osp_230[k];

        t_461[k] = pa_y[k] * nsd0_395[k]
                   - f_4 * pc_y[k] * nsd1_395[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, t_467, pc_x, pc_y, pc_z, nsp_197, \
                         oss0_77, oss1_77, osp_231, osp_232, osp_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_1 * oss0_77[k]
                   - f_2 * oss1_77[k]
                   + f_3 * pc_x[k] * osp_231[k];

        t_463[k] = f_3 * pc_x[k] * osp_232[k];

        t_464[k] = f_3 * pc_x[k] * osp_233[k];

        t_465[k] = f_1 * oss0_77[k]
                   - f_2 * oss1_77[k]
                   + f_3 * pc_y[k] * osp_232[k];

        t_466[k] = f_3 * pc_y[k] * osp_233[k];

        t_467[k] = f_0 * nsp_197[k]
                   + f_1 * oss0_77[k]
                   - f_2 * oss1_77[k]
                   + f_3 * pc_z[k] * osp_233[k];
    }
}

auto
compute_prim_osd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t nsd0, const size_t nsp,
                                                   const size_t nsd1, const size_t oss0,
                                                   const size_t oss1, const size_t osp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_osd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, nsd0, nsp,
                                                              nsd1, oss0, oss1, osp, ncols,
                                                              gamma, p, q);

    compute_prim_osd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, nsd0, nsp,
                                                              nsd1, oss0, oss1, osp, ncols,
                                                              gamma, p, q);

    compute_prim_osd_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, nsd0, nsp,
                                                              nsd1, oss0, oss1, osp, ncols,
                                                              gamma, p, q);

    compute_prim_osd_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, nsd0, nsp,
                                                              nsd1, oss0, oss1, osp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
