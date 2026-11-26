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


#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fpd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppd0, const size_t ppd1,
                                                          const size_t dpd0, const size_t dpp,
                                                          const size_t dpd1, const size_t fsd0,
                                                          const size_t fsp, const size_t fsd1,
                                                          const size_t fps0, const size_t fps1,
                                                          const size_t fpp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppd0_27 = buffer.data(ppd0 + 27);
    const auto *ppd0_53 = buffer.data(ppd0 + 53);

    const auto *ppd1_27 = buffer.data(ppd1 + 27);
    const auto *ppd1_53 = buffer.data(ppd1 + 53);

    const auto *dpd0_0 = buffer.data(dpd0 + 0);
    const auto *dpd0_3 = buffer.data(dpd0 + 3);
    const auto *dpd0_5 = buffer.data(dpd0 + 5);
    const auto *dpd0_6 = buffer.data(dpd0 + 6);
    const auto *dpd0_9 = buffer.data(dpd0 + 9);
    const auto *dpd0_12 = buffer.data(dpd0 + 12);
    const auto *dpd0_17 = buffer.data(dpd0 + 17);
    const auto *dpd0_19 = buffer.data(dpd0 + 19);
    const auto *dpd0_27 = buffer.data(dpd0 + 27);
    const auto *dpd0_36 = buffer.data(dpd0 + 36);
    const auto *dpd0_38 = buffer.data(dpd0 + 38);
    const auto *dpd0_48 = buffer.data(dpd0 + 48);
    const auto *dpd0_53 = buffer.data(dpd0 + 53);
    const auto *dpd0_54 = buffer.data(dpd0 + 54);
    const auto *dpd0_60 = buffer.data(dpd0 + 60);
    const auto *dpd0_63 = buffer.data(dpd0 + 63);
    const auto *dpd0_65 = buffer.data(dpd0 + 65);
    const auto *dpd0_69 = buffer.data(dpd0 + 69);
    const auto *dpd0_71 = buffer.data(dpd0 + 71);
    const auto *dpd0_81 = buffer.data(dpd0 + 81);
    const auto *dpd0_82 = buffer.data(dpd0 + 82);
    const auto *dpd0_83 = buffer.data(dpd0 + 83);
    const auto *dpd0_87 = buffer.data(dpd0 + 87);
    const auto *dpd0_89 = buffer.data(dpd0 + 89);
    const auto *dpd0_99 = buffer.data(dpd0 + 99);
    const auto *dpd0_101 = buffer.data(dpd0 + 101);
    const auto *dpd0_102 = buffer.data(dpd0 + 102);
    const auto *dpd0_105 = buffer.data(dpd0 + 105);
    const auto *dpd0_107 = buffer.data(dpd0 + 107);

    const auto *dpp_0 = buffer.data(dpp + 0);
    const auto *dpp_1 = buffer.data(dpp + 1);
    const auto *dpp_2 = buffer.data(dpp + 2);
    const auto *dpp_10 = buffer.data(dpp + 10);
    const auto *dpp_11 = buffer.data(dpp + 11);
    const auto *dpp_12 = buffer.data(dpp + 12);
    const auto *dpp_13 = buffer.data(dpp + 13);
    const auto *dpp_16 = buffer.data(dpp + 16);
    const auto *dpp_19 = buffer.data(dpp + 19);
    const auto *dpp_20 = buffer.data(dpp + 20);
    const auto *dpp_23 = buffer.data(dpp + 23);
    const auto *dpp_24 = buffer.data(dpp + 24);
    const auto *dpp_26 = buffer.data(dpp + 26);
    const auto *dpp_27 = buffer.data(dpp + 27);
    const auto *dpp_28 = buffer.data(dpp + 28);
    const auto *dpp_30 = buffer.data(dpp + 30);
    const auto *dpp_31 = buffer.data(dpp + 31);
    const auto *dpp_34 = buffer.data(dpp + 34);
    const auto *dpp_39 = buffer.data(dpp + 39);
    const auto *dpp_40 = buffer.data(dpp + 40);
    const auto *dpp_41 = buffer.data(dpp + 41);
    const auto *dpp_43 = buffer.data(dpp + 43);
    const auto *dpp_44 = buffer.data(dpp + 44);
    const auto *dpp_45 = buffer.data(dpp + 45);
    const auto *dpp_47 = buffer.data(dpp + 47);
    const auto *dpp_50 = buffer.data(dpp + 50);
    const auto *dpp_51 = buffer.data(dpp + 51);
    const auto *dpp_53 = buffer.data(dpp + 53);

    const auto *dpd1_0 = buffer.data(dpd1 + 0);
    const auto *dpd1_3 = buffer.data(dpd1 + 3);
    const auto *dpd1_5 = buffer.data(dpd1 + 5);
    const auto *dpd1_6 = buffer.data(dpd1 + 6);
    const auto *dpd1_9 = buffer.data(dpd1 + 9);
    const auto *dpd1_12 = buffer.data(dpd1 + 12);
    const auto *dpd1_17 = buffer.data(dpd1 + 17);
    const auto *dpd1_19 = buffer.data(dpd1 + 19);
    const auto *dpd1_27 = buffer.data(dpd1 + 27);
    const auto *dpd1_36 = buffer.data(dpd1 + 36);
    const auto *dpd1_38 = buffer.data(dpd1 + 38);
    const auto *dpd1_48 = buffer.data(dpd1 + 48);
    const auto *dpd1_53 = buffer.data(dpd1 + 53);
    const auto *dpd1_54 = buffer.data(dpd1 + 54);
    const auto *dpd1_60 = buffer.data(dpd1 + 60);
    const auto *dpd1_63 = buffer.data(dpd1 + 63);
    const auto *dpd1_65 = buffer.data(dpd1 + 65);
    const auto *dpd1_69 = buffer.data(dpd1 + 69);
    const auto *dpd1_71 = buffer.data(dpd1 + 71);
    const auto *dpd1_81 = buffer.data(dpd1 + 81);
    const auto *dpd1_82 = buffer.data(dpd1 + 82);
    const auto *dpd1_83 = buffer.data(dpd1 + 83);
    const auto *dpd1_87 = buffer.data(dpd1 + 87);
    const auto *dpd1_89 = buffer.data(dpd1 + 89);
    const auto *dpd1_99 = buffer.data(dpd1 + 99);
    const auto *dpd1_101 = buffer.data(dpd1 + 101);
    const auto *dpd1_102 = buffer.data(dpd1 + 102);
    const auto *dpd1_105 = buffer.data(dpd1 + 105);
    const auto *dpd1_107 = buffer.data(dpd1 + 107);

    const auto *fsd0_0 = buffer.data(fsd0 + 0);
    const auto *fsd0_3 = buffer.data(fsd0 + 3);
    const auto *fsd0_5 = buffer.data(fsd0 + 5);
    const auto *fsd0_9 = buffer.data(fsd0 + 9);
    const auto *fsd0_17 = buffer.data(fsd0 + 17);
    const auto *fsd0_18 = buffer.data(fsd0 + 18);
    const auto *fsd0_30 = buffer.data(fsd0 + 30);
    const auto *fsd0_36 = buffer.data(fsd0 + 36);
    const auto *fsd0_39 = buffer.data(fsd0 + 39);
    const auto *fsd0_41 = buffer.data(fsd0 + 41);

    const auto *fsp_0 = buffer.data(fsp + 0);
    const auto *fsp_1 = buffer.data(fsp + 1);
    const auto *fsp_2 = buffer.data(fsp + 2);
    const auto *fsp_3 = buffer.data(fsp + 3);
    const auto *fsp_4 = buffer.data(fsp + 4);
    const auto *fsp_6 = buffer.data(fsp + 6);
    const auto *fsp_8 = buffer.data(fsp + 8);
    const auto *fsp_9 = buffer.data(fsp + 9);
    const auto *fsp_10 = buffer.data(fsp + 10);
    const auto *fsp_15 = buffer.data(fsp + 15);
    const auto *fsp_17 = buffer.data(fsp + 17);
    const auto *fsp_18 = buffer.data(fsp + 18);
    const auto *fsp_19 = buffer.data(fsp + 19);
    const auto *fsp_20 = buffer.data(fsp + 20);
    const auto *fsp_22 = buffer.data(fsp + 22);

    const auto *fsd1_0 = buffer.data(fsd1 + 0);
    const auto *fsd1_3 = buffer.data(fsd1 + 3);
    const auto *fsd1_5 = buffer.data(fsd1 + 5);
    const auto *fsd1_9 = buffer.data(fsd1 + 9);
    const auto *fsd1_17 = buffer.data(fsd1 + 17);
    const auto *fsd1_18 = buffer.data(fsd1 + 18);
    const auto *fsd1_30 = buffer.data(fsd1 + 30);
    const auto *fsd1_36 = buffer.data(fsd1 + 36);
    const auto *fsd1_39 = buffer.data(fsd1 + 39);
    const auto *fsd1_41 = buffer.data(fsd1 + 41);

    const auto *fps0_0 = buffer.data(fps0 + 0);
    const auto *fps0_3 = buffer.data(fps0 + 3);
    const auto *fps0_4 = buffer.data(fps0 + 4);
    const auto *fps0_6 = buffer.data(fps0 + 6);
    const auto *fps0_8 = buffer.data(fps0 + 8);
    const auto *fps0_9 = buffer.data(fps0 + 9);
    const auto *fps0_12 = buffer.data(fps0 + 12);
    const auto *fps0_13 = buffer.data(fps0 + 13);
    const auto *fps0_15 = buffer.data(fps0 + 15);
    const auto *fps0_19 = buffer.data(fps0 + 19);

    const auto *fps1_0 = buffer.data(fps1 + 0);
    const auto *fps1_3 = buffer.data(fps1 + 3);
    const auto *fps1_4 = buffer.data(fps1 + 4);
    const auto *fps1_6 = buffer.data(fps1 + 6);
    const auto *fps1_8 = buffer.data(fps1 + 8);
    const auto *fps1_9 = buffer.data(fps1 + 9);
    const auto *fps1_12 = buffer.data(fps1 + 12);
    const auto *fps1_13 = buffer.data(fps1 + 13);
    const auto *fps1_15 = buffer.data(fps1 + 15);
    const auto *fps1_19 = buffer.data(fps1 + 19);

    const auto *fpp_0 = buffer.data(fpp + 0);
    const auto *fpp_1 = buffer.data(fpp + 1);
    const auto *fpp_2 = buffer.data(fpp + 2);
    const auto *fpp_3 = buffer.data(fpp + 3);
    const auto *fpp_5 = buffer.data(fpp + 5);
    const auto *fpp_6 = buffer.data(fpp + 6);
    const auto *fpp_8 = buffer.data(fpp + 8);
    const auto *fpp_9 = buffer.data(fpp + 9);
    const auto *fpp_10 = buffer.data(fpp + 10);
    const auto *fpp_12 = buffer.data(fpp + 12);
    const auto *fpp_13 = buffer.data(fpp + 13);
    const auto *fpp_14 = buffer.data(fpp + 14);
    const auto *fpp_15 = buffer.data(fpp + 15);
    const auto *fpp_16 = buffer.data(fpp + 16);
    const auto *fpp_18 = buffer.data(fpp + 18);
    const auto *fpp_20 = buffer.data(fpp + 20);
    const auto *fpp_21 = buffer.data(fpp + 21);
    const auto *fpp_23 = buffer.data(fpp + 23);
    const auto *fpp_24 = buffer.data(fpp + 24);
    const auto *fpp_25 = buffer.data(fpp + 25);
    const auto *fpp_26 = buffer.data(fpp + 26);
    const auto *fpp_27 = buffer.data(fpp + 27);
    const auto *fpp_28 = buffer.data(fpp + 28);
    const auto *fpp_29 = buffer.data(fpp + 29);
    const auto *fpp_30 = buffer.data(fpp + 30);
    const auto *fpp_31 = buffer.data(fpp + 31);
    const auto *fpp_33 = buffer.data(fpp + 33);
    const auto *fpp_34 = buffer.data(fpp + 34);
    const auto *fpp_37 = buffer.data(fpp + 37);
    const auto *fpp_38 = buffer.data(fpp + 38);
    const auto *fpp_39 = buffer.data(fpp + 39);
    const auto *fpp_40 = buffer.data(fpp + 40);
    const auto *fpp_41 = buffer.data(fpp + 41);
    const auto *fpp_43 = buffer.data(fpp + 43);
    const auto *fpp_44 = buffer.data(fpp + 44);
    const auto *fpp_45 = buffer.data(fpp + 45);
    const auto *fpp_46 = buffer.data(fpp + 46);
    const auto *fpp_47 = buffer.data(fpp + 47);
    const auto *fpp_48 = buffer.data(fpp + 48);
    const auto *fpp_50 = buffer.data(fpp + 50);
    const auto *fpp_51 = buffer.data(fpp + 51);
    const auto *fpp_53 = buffer.data(fpp + 53);
    const auto *fpp_55 = buffer.data(fpp + 55);
    const auto *fpp_56 = buffer.data(fpp + 56);
    const auto *fpp_57 = buffer.data(fpp + 57);
    const auto *fpp_58 = buffer.data(fpp + 58);
    const auto *fpp_59 = buffer.data(fpp + 59);
    const auto *fpp_61 = buffer.data(fpp + 61);
    const auto *fpp_62 = buffer.data(fpp + 62);
    const auto *fpp_64 = buffer.data(fpp + 64);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dpp_0, fsp_0, fps0_0, \
                         fps1_0, fpp_0, fpp_1, fpp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dpp_0[k]
                 + f_1 * fsp_0[k]
                 + f_2 * fps0_0[k]
                 - f_3 * fps1_0[k]
                 + f_4 * pc_x[k] * fpp_0[k];

        t_1[k] = f_4 * pc_y[k] * fpp_0[k];

        t_2[k] = f_4 * pc_z[k] * fpp_0[k];

        t_3[k] = f_2 * fps0_0[k]
                 - f_3 * fps1_0[k]
                 + f_4 * pc_y[k] * fpp_1[k];

        t_4[k] = f_4 * pc_y[k] * fpp_2[k];

        t_5[k] = f_2 * fps0_0[k]
                 - f_3 * fps1_0[k]
                 + f_4 * pc_z[k] * fpp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pc_y, pc_z, fsd0_0, fsd0_3, fsp_0, fsp_1, \
                         fsd1_0, fsd1_3, fpp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * fsd0_0[k]
                 - f_5 * pc_y[k] * fsd1_0[k];

        t_7[k] = f_1 * fsp_0[k]
                 + f_4 * pc_y[k] * fpp_3[k];

        t_8[k] = f_4 * pc_z[k] * fpp_3[k];

        t_9[k] = pb_y[k] * fsd0_3[k]
                 + f_6 * fsp_1[k]
                 - f_5 * pc_y[k] * fsd1_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, pc_y, pc_z, fsd0_0, fsd0_5, \
                         fsp_2, fsd1_0, fsd1_5, fpp_5, fpp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * fsp_2[k]
                  + f_4 * pc_y[k] * fpp_5[k];

        t_11[k] = pb_y[k] * fsd0_5[k]
                  - f_5 * pc_y[k] * fsd1_5[k];

        t_12[k] = pb_z[k] * fsd0_0[k]
                  - f_5 * pc_z[k] * fsd1_0[k];

        t_13[k] = f_4 * pc_y[k] * fpp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_z, pc_y, pc_z, fsd0_3, fsd0_5, fsp_0, \
                         fsp_2, fsd1_3, fsd1_5, fpp_6, fpp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * fsp_0[k]
                  + f_4 * pc_z[k] * fpp_6[k];

        t_15[k] = pb_z[k] * fsd0_3[k]
                  - f_5 * pc_z[k] * fsd1_3[k];

        t_16[k] = f_4 * pc_y[k] * fpp_8[k];

        t_17[k] = pb_z[k] * fsd0_5[k]
                  + f_6 * fsp_2[k]
                  - f_5 * pc_z[k] * fsd1_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pc_x, pc_y, pc_z, dpd0_0, dpp_10, dpd1_0, \
                         fsp_4, fpp_9, fpp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_y[k] * dpd0_0[k]
                  - f_5 * pc_y[k] * dpd1_0[k];

        t_19[k] = f_6 * dpp_10[k]
                  + f_1 * fsp_4[k]
                  + f_4 * pc_x[k] * fpp_10[k];

        t_20[k] = f_4 * pc_z[k] * fpp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pc_y, pc_z, dpd0_5, dpp_1, dpd1_5, fps0_3, \
                         fps1_3, fpp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * dpp_1[k]
                  + f_2 * fps0_3[k]
                  - f_3 * fps1_3[k]
                  + f_4 * pc_y[k] * fpp_10[k];

        t_22[k] = f_4 * pc_z[k] * fpp_10[k];

        t_23[k] = pa_y[k] * dpd0_5[k]
                  - f_5 * pc_y[k] * dpd1_5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pc_x, pc_z, dpp_12, dpp_13, fps0_4, fps1_4, fpp_12, \
                         fpp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * dpp_12[k]
                  + f_2 * fps0_4[k]
                  - f_3 * fps1_4[k]
                  + f_4 * pc_x[k] * fpp_12[k];

        t_25[k] = f_6 * dpp_13[k]
                  + f_4 * pc_x[k] * fpp_13[k];

        t_26[k] = f_4 * pc_z[k] * fpp_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pc_x, pc_z, ppd0_27, ppd1_27, dpd0_27, \
                         dpd1_27, fps0_4, fps1_4, fpp_13, fpp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * ppd0_27[k]
                  - f_8 * ppd1_27[k]
                  + pa_x[k] * dpd0_27[k]
                  - f_5 * pc_x[k] * dpd1_27[k];

        t_28[k] = f_4 * pc_z[k] * fpp_13[k];

        t_29[k] = f_2 * fps0_4[k]
                  - f_3 * fps1_4[k]
                  + f_4 * pc_z[k] * fpp_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pc_x, pc_y, pc_z, dpd0_12, dpp_16, dpd1_12, \
                         fsp_3, fpp_15, fpp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * dpd0_12[k]
                  - f_5 * pc_y[k] * dpd1_12[k];

        t_31[k] = f_6 * dpp_16[k]
                  + f_4 * pc_x[k] * fpp_16[k];

        t_32[k] = f_1 * fsp_3[k]
                  + f_4 * pc_z[k] * fpp_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_z, pc_y, pc_z, dpd0_17, dpd1_17, fsd0_9, \
                         fsp_4, fsd1_9, fpp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * fsd0_9[k]
                  - f_5 * pc_z[k] * fsd1_9[k];

        t_34[k] = f_1 * fsp_4[k]
                  + f_4 * pc_z[k] * fpp_16[k];

        t_35[k] = pa_y[k] * dpd0_17[k]
                  - f_5 * pc_y[k] * dpd1_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pc_x, pc_y, pc_z, dpd0_0, dpd0_3, \
                         dpp_20, dpd1_0, dpd1_3, fsp_8, fpp_18, \
                         fpp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_z[k] * dpd0_0[k]
                  - f_5 * pc_z[k] * dpd1_0[k];

        t_37[k] = f_4 * pc_y[k] * fpp_18[k];

        t_38[k] = f_6 * dpp_20[k]
                  + f_1 * fsp_8[k]
                  + f_4 * pc_x[k] * fpp_20[k];

        t_39[k] = pa_z[k] * dpd0_3[k]
                  - f_5 * pc_z[k] * dpd1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pc_y, pc_z, dpd0_6, dpp_2, dpd1_6, \
                         fsp_6, fps0_6, fps1_6, fpp_20, fpp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_4 * pc_y[k] * fpp_20[k];

        t_41[k] = f_1 * dpp_2[k]
                  + f_2 * fps0_6[k]
                  - f_3 * fps1_6[k]
                  + f_4 * pc_z[k] * fpp_20[k];

        t_42[k] = pa_z[k] * dpd0_6[k]
                  - f_5 * pc_z[k] * dpd1_6[k];

        t_43[k] = f_1 * fsp_6[k]
                  + f_4 * pc_y[k] * fpp_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_y, pc_x, pc_y, pc_z, dpd0_9, dpp_23, \
                         dpd1_9, fsd0_17, fsp_8, fsd1_17, fpp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_6 * dpp_23[k]
                  + f_4 * pc_x[k] * fpp_23[k];

        t_45[k] = pa_z[k] * dpd0_9[k]
                  - f_5 * pc_z[k] * dpd1_9[k];

        t_46[k] = f_1 * fsp_8[k]
                  + f_4 * pc_y[k] * fpp_23[k];

        t_47[k] = pb_y[k] * fsd0_17[k]
                  - f_5 * pc_y[k] * fsd1_17[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pc_x, pc_y, dpp_24, dpp_26, fps0_8, \
                         fps1_8, fpp_24, fpp_25, fpp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * dpp_24[k]
                  + f_2 * fps0_8[k]
                  - f_3 * fps1_8[k]
                  + f_4 * pc_x[k] * fpp_24[k];

        t_49[k] = f_4 * pc_y[k] * fpp_24[k];

        t_50[k] = f_6 * dpp_26[k]
                  + f_4 * pc_x[k] * fpp_26[k];

        t_51[k] = f_2 * fps0_8[k]
                  - f_3 * fps1_8[k]
                  + f_4 * pc_y[k] * fpp_25[k];

        t_52[k] = f_4 * pc_y[k] * fpp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, pa_x, pc_x, ppd0_53, ppd1_53, dpd0_53, dpp_27, dpd1_53, \
                         fsp_9, fps0_9, fps1_9, fpp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * ppd0_53[k]
                  - f_8 * ppd1_53[k]
                  + pa_x[k] * dpd0_53[k]
                  - f_5 * pc_x[k] * dpd1_53[k];

        t_54[k] = f_1 * dpp_27[k]
                  + f_1 * fsp_9[k]
                  + f_2 * fps0_9[k]
                  - f_3 * fps1_9[k]
                  + f_4 * pc_x[k] * fpp_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, dpp_10, dpp_28, \
                         fsp_10, fps0_9, fps1_9, fpp_27, fpp_28, \
                         fpp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * dpp_28[k]
                  + f_1 * fsp_10[k]
                  + f_4 * pc_x[k] * fpp_28[k];

        t_56[k] = f_4 * pc_z[k] * fpp_27[k];

        t_57[k] = f_6 * dpp_10[k]
                  + f_2 * fps0_9[k]
                  - f_3 * fps1_9[k]
                  + f_4 * pc_y[k] * fpp_28[k];

        t_58[k] = f_4 * pc_z[k] * fpp_28[k];

        t_59[k] = f_2 * fps0_9[k]
                  - f_3 * fps1_9[k]
                  + f_4 * pc_z[k] * fpp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_x, pc_x, pc_z, dpd0_60, dpd0_63, \
                         dpp_30, dpp_31, dpd1_60, dpd1_63, fpp_30, \
                         fpp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * dpd0_60[k]
                  + f_6 * dpp_30[k]
                  - f_5 * pc_x[k] * dpd1_60[k];

        t_61[k] = f_1 * dpp_31[k]
                  + f_4 * pc_x[k] * fpp_31[k];

        t_62[k] = f_4 * pc_z[k] * fpp_30[k];

        t_63[k] = pa_x[k] * dpd0_63[k]
                  - f_5 * pc_x[k] * dpd1_63[k];

        t_64[k] = f_4 * pc_z[k] * fpp_31[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_z, pc_x, pc_z, dpd0_65, dpp_34, \
                         dpd1_65, fsd0_18, fsp_9, fsd1_18, fpp_33, \
                         fpp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_x[k] * dpd0_65[k]
                  - f_5 * pc_x[k] * dpd1_65[k];

        t_66[k] = pb_z[k] * fsd0_18[k]
                  - f_5 * pc_z[k] * fsd1_18[k];

        t_67[k] = f_1 * dpp_34[k]
                  + f_4 * pc_x[k] * fpp_34[k];

        t_68[k] = f_1 * fsp_9[k]
                  + f_4 * pc_z[k] * fpp_33[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_x, pc_x, pc_z, dpd0_69, dpd0_71, dpd1_69, \
                         dpd1_71, fsp_10, fpp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_x[k] * dpd0_69[k]
                  - f_5 * pc_x[k] * dpd1_69[k];

        t_70[k] = f_1 * fsp_10[k]
                  + f_4 * pc_z[k] * fpp_34[k];

        t_71[k] = pa_x[k] * dpd0_71[k]
                  - f_5 * pc_x[k] * dpd1_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pa_z, pc_y, pc_z, dpd0_19, dpd0_36, dpd0_38, \
                         dpd1_19, dpd1_36, dpd1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * dpd0_36[k]
                  - f_5 * pc_y[k] * dpd1_36[k];

        t_73[k] = pa_z[k] * dpd0_19[k]
                  - f_5 * pc_z[k] * dpd1_19[k];

        t_74[k] = pa_y[k] * dpd0_38[k]
                  - f_5 * pc_y[k] * dpd1_38[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pc_y, pc_z, dpp_11, dpp_19, dpp_20, fps0_12, \
                         fps1_12, fpp_37, fpp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * dpp_19[k]
                  + f_2 * fps0_12[k]
                  - f_3 * fps1_12[k]
                  + f_4 * pc_y[k] * fpp_37[k];

        t_76[k] = f_1 * dpp_20[k]
                  + f_4 * pc_y[k] * fpp_38[k];

        t_77[k] = f_1 * dpp_11[k]
                  + f_2 * fps0_12[k]
                  - f_3 * fps1_12[k]
                  + f_4 * pc_z[k] * fpp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pc_x, dpd0_81, dpp_39, dpp_40, dpp_41, \
                         dpd1_81, fps0_13, fps1_13, fpp_39, fpp_40, \
                         fpp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * dpp_39[k]
                  + f_2 * fps0_13[k]
                  - f_3 * fps1_13[k]
                  + f_4 * pc_x[k] * fpp_39[k];

        t_79[k] = f_1 * dpp_40[k]
                  + f_4 * pc_x[k] * fpp_40[k];

        t_80[k] = f_1 * dpp_41[k]
                  + f_4 * pc_x[k] * fpp_41[k];

        t_81[k] = pa_x[k] * dpd0_81[k]
                  - f_5 * pc_x[k] * dpd1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, pc_x, pc_y, dpd0_48, dpd0_82, \
                         dpd0_83, dpp_43, dpd1_48, dpd1_82, dpd1_83, \
                         fpp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_x[k] * dpd0_82[k]
                  - f_5 * pc_x[k] * dpd1_82[k];

        t_83[k] = pa_x[k] * dpd0_83[k]
                  - f_5 * pc_x[k] * dpd1_83[k];

        t_84[k] = pa_y[k] * dpd0_48[k]
                  - f_5 * pc_y[k] * dpd1_48[k];

        t_85[k] = f_1 * dpp_43[k]
                  + f_4 * pc_x[k] * fpp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pc_x, pc_y, dpd0_87, dpd0_89, dpp_26, \
                         dpp_44, dpd1_87, dpd1_89, fpp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_1 * dpp_44[k]
                  + f_4 * pc_x[k] * fpp_44[k];

        t_87[k] = pa_x[k] * dpd0_87[k]
                  - f_5 * pc_x[k] * dpd1_87[k];

        t_88[k] = f_1 * dpp_26[k]
                  + f_4 * pc_y[k] * fpp_44[k];

        t_89[k] = pa_x[k] * dpd0_89[k]
                  - f_5 * pc_x[k] * dpd1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, dpp_45, dpp_47, fsp_15, \
                         fsp_17, fps0_15, fps1_15, fpp_45, fpp_46, \
                         fpp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * dpp_45[k]
                  + f_1 * fsp_15[k]
                  + f_2 * fps0_15[k]
                  - f_3 * fps1_15[k]
                  + f_4 * pc_x[k] * fpp_45[k];

        t_91[k] = f_4 * pc_y[k] * fpp_45[k];

        t_92[k] = f_1 * dpp_47[k]
                  + f_1 * fsp_17[k]
                  + f_4 * pc_x[k] * fpp_47[k];

        t_93[k] = f_2 * fps0_15[k]
                  - f_3 * fps1_15[k]
                  + f_4 * pc_y[k] * fpp_46[k];

        t_94[k] = f_4 * pc_y[k] * fpp_47[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_y, pc_y, pc_z, dpp_20, fsd0_30, fsp_15, fsd1_30, \
                         fps0_15, fps1_15, fpp_47, fpp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_6 * dpp_20[k]
                  + f_2 * fps0_15[k]
                  - f_3 * fps1_15[k]
                  + f_4 * pc_z[k] * fpp_47[k];

        t_96[k] = pb_y[k] * fsd0_30[k]
                  - f_5 * pc_y[k] * fsd1_30[k];

        t_97[k] = f_1 * fsp_15[k]
                  + f_4 * pc_y[k] * fpp_48[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_x, pc_x, pc_y, dpd0_99, dpd0_101, \
                         dpp_50, dpd1_99, dpd1_101, fsp_17, fpp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_1 * dpp_50[k]
                  + f_4 * pc_x[k] * fpp_50[k];

        t_99[k] = pa_x[k] * dpd0_99[k]
                  - f_5 * pc_x[k] * dpd1_99[k];

        t_100[k] = f_1 * fsp_17[k]
                   + f_4 * pc_y[k] * fpp_50[k];

        t_101[k] = pa_x[k] * dpd0_101[k]
                   - f_5 * pc_x[k] * dpd1_101[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_x, pc_x, pc_y, dpd0_102, \
                         dpd0_105, dpp_51, dpp_53, dpd1_102, dpd1_105, fpp_51, \
                         fpp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pa_x[k] * dpd0_102[k]
                   + f_6 * dpp_51[k]
                   - f_5 * pc_x[k] * dpd1_102[k];

        t_103[k] = f_4 * pc_y[k] * fpp_51[k];

        t_104[k] = f_1 * dpp_53[k]
                   + f_4 * pc_x[k] * fpp_53[k];

        t_105[k] = pa_x[k] * dpd0_105[k]
                   - f_5 * pc_x[k] * dpd1_105[k];

        t_106[k] = f_4 * pc_y[k] * fpp_53[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, pc_x, dpd0_107, dpd1_107, \
                         fsd0_36, fsp_18, fsp_19, fsp_20, fsd1_36, fpp_55, \
                         fpp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * dpd0_107[k]
                   - f_5 * pc_x[k] * dpd1_107[k];

        t_108[k] = pb_x[k] * fsd0_36[k]
                   + f_6 * fsp_18[k]
                   - f_5 * pc_x[k] * fsd1_36[k];

        t_109[k] = f_1 * fsp_19[k]
                   + f_4 * pc_x[k] * fpp_55[k];

        t_110[k] = f_1 * fsp_20[k]
                   + f_4 * pc_x[k] * fpp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pc_x, pc_z, fsd0_39, fsd0_41, \
                         fsd1_39, fsd1_41, fps0_19, fps1_19, fpp_55, \
                         fpp_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pb_x[k] * fsd0_39[k]
                   - f_5 * pc_x[k] * fsd1_39[k];

        t_112[k] = f_4 * pc_z[k] * fpp_55[k];

        t_113[k] = pb_x[k] * fsd0_41[k]
                   - f_5 * pc_x[k] * fsd1_41[k];

        t_114[k] = f_2 * fps0_19[k]
                   - f_3 * fps1_19[k]
                   + f_4 * pc_x[k] * fpp_57[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pc_x, pc_y, pc_z, dpp_31, fsp_19, \
                         fps0_19, fps1_19, fpp_58, fpp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_4 * pc_x[k] * fpp_58[k];

        t_116[k] = f_4 * pc_x[k] * fpp_59[k];

        t_117[k] = f_0 * dpp_31[k]
                   + f_1 * fsp_19[k]
                   + f_2 * fps0_19[k]
                   - f_3 * fps1_19[k]
                   + f_4 * pc_y[k] * fpp_58[k];

        t_118[k] = f_4 * pc_z[k] * fpp_58[k];

        t_119[k] = f_2 * fps0_19[k]
                   - f_3 * fps1_19[k]
                   + f_4 * pc_z[k] * fpp_59[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pb_z, pc_x, pc_z, fsd0_36, \
                         fsd0_39, fsp_19, fsd1_36, fsd1_39, fpp_61, \
                         fpp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_z[k] * fsd0_36[k]
                   - f_5 * pc_z[k] * fsd1_36[k];

        t_121[k] = f_4 * pc_x[k] * fpp_61[k];

        t_122[k] = f_4 * pc_x[k] * fpp_62[k];

        t_123[k] = pb_z[k] * fsd0_39[k]
                   - f_5 * pc_z[k] * fsd1_39[k];

        t_124[k] = f_1 * fsp_19[k]
                   + f_4 * pc_z[k] * fpp_61[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_z, pb_z, pc_x, pc_z, dpd0_54, dpd1_54, \
                         fsd0_41, fsp_20, fsp_22, fsd1_41, fpp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_z[k] * fsd0_41[k]
                   + f_6 * fsp_20[k]
                   - f_5 * pc_z[k] * fsd1_41[k];

        t_126[k] = pa_z[k] * dpd0_54[k]
                   - f_5 * pc_z[k] * dpd1_54[k];

        t_127[k] = f_1 * fsp_22[k]
                   + f_4 * pc_x[k] * fpp_64[k];
    }
}

static auto
compute_prim_fpd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppd0, const size_t ppd1,
                                                          const size_t dpd0, const size_t dpp,
                                                          const size_t dpd1, const size_t fsd0,
                                                          const size_t fsp, const size_t fsd1,
                                                          const size_t fps0, const size_t fps1,
                                                          const size_t fpp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppd0_53 = buffer.data(ppd0 + 53);

    const auto *ppd1_53 = buffer.data(ppd1 + 53);

    const auto *dpd0_57 = buffer.data(dpd0 + 57);
    const auto *dpd0_60 = buffer.data(dpd0 + 60);
    const auto *dpd0_63 = buffer.data(dpd0 + 63);
    const auto *dpd0_89 = buffer.data(dpd0 + 89);
    const auto *dpd0_90 = buffer.data(dpd0 + 90);
    const auto *dpd0_95 = buffer.data(dpd0 + 95);
    const auto *dpd0_102 = buffer.data(dpd0 + 102);
    const auto *dpd0_105 = buffer.data(dpd0 + 105);
    const auto *dpd0_107 = buffer.data(dpd0 + 107);

    const auto *dpp_32 = buffer.data(dpp + 32);
    const auto *dpp_38 = buffer.data(dpp + 38);
    const auto *dpp_41 = buffer.data(dpp + 41);
    const auto *dpp_43 = buffer.data(dpp + 43);
    const auto *dpp_44 = buffer.data(dpp + 44);
    const auto *dpp_47 = buffer.data(dpp + 47);
    const auto *dpp_49 = buffer.data(dpp + 49);
    const auto *dpp_50 = buffer.data(dpp + 50);
    const auto *dpp_52 = buffer.data(dpp + 52);
    const auto *dpp_53 = buffer.data(dpp + 53);

    const auto *dpd1_57 = buffer.data(dpd1 + 57);
    const auto *dpd1_60 = buffer.data(dpd1 + 60);
    const auto *dpd1_63 = buffer.data(dpd1 + 63);
    const auto *dpd1_89 = buffer.data(dpd1 + 89);
    const auto *dpd1_90 = buffer.data(dpd1 + 90);
    const auto *dpd1_95 = buffer.data(dpd1 + 95);
    const auto *dpd1_102 = buffer.data(dpd1 + 102);
    const auto *dpd1_105 = buffer.data(dpd1 + 105);
    const auto *dpd1_107 = buffer.data(dpd1 + 107);

    const auto *fsd0_47 = buffer.data(fsd0 + 47);
    const auto *fsd0_51 = buffer.data(fsd0 + 51);
    const auto *fsd0_54 = buffer.data(fsd0 + 54);
    const auto *fsd0_57 = buffer.data(fsd0 + 57);
    const auto *fsd0_59 = buffer.data(fsd0 + 59);

    const auto *fsp_23 = buffer.data(fsp + 23);
    const auto *fsp_25 = buffer.data(fsp + 25);
    const auto *fsp_26 = buffer.data(fsp + 26);
    const auto *fsp_27 = buffer.data(fsp + 27);
    const auto *fsp_28 = buffer.data(fsp + 28);
    const auto *fsp_29 = buffer.data(fsp + 29);

    const auto *fsd1_47 = buffer.data(fsd1 + 47);
    const auto *fsd1_51 = buffer.data(fsd1 + 51);
    const auto *fsd1_54 = buffer.data(fsd1 + 54);
    const auto *fsd1_57 = buffer.data(fsd1 + 57);
    const auto *fsd1_59 = buffer.data(fsd1 + 59);

    const auto *fps0_22 = buffer.data(fps0 + 22);
    const auto *fps0_23 = buffer.data(fps0 + 23);
    const auto *fps0_25 = buffer.data(fps0 + 25);
    const auto *fps0_29 = buffer.data(fps0 + 29);

    const auto *fps1_22 = buffer.data(fps1 + 22);
    const auto *fps1_23 = buffer.data(fps1 + 23);
    const auto *fps1_25 = buffer.data(fps1 + 25);
    const auto *fps1_29 = buffer.data(fps1 + 29);

    const auto *fpp_65 = buffer.data(fpp + 65);
    const auto *fpp_67 = buffer.data(fpp + 67);
    const auto *fpp_68 = buffer.data(fpp + 68);
    const auto *fpp_69 = buffer.data(fpp + 69);
    const auto *fpp_70 = buffer.data(fpp + 70);
    const auto *fpp_71 = buffer.data(fpp + 71);
    const auto *fpp_73 = buffer.data(fpp + 73);
    const auto *fpp_74 = buffer.data(fpp + 74);
    const auto *fpp_75 = buffer.data(fpp + 75);
    const auto *fpp_76 = buffer.data(fpp + 76);
    const auto *fpp_77 = buffer.data(fpp + 77);
    const auto *fpp_79 = buffer.data(fpp + 79);
    const auto *fpp_80 = buffer.data(fpp + 80);
    const auto *fpp_82 = buffer.data(fpp + 82);
    const auto *fpp_83 = buffer.data(fpp + 83);
    const auto *fpp_85 = buffer.data(fpp + 85);
    const auto *fpp_86 = buffer.data(fpp + 86);
    const auto *fpp_87 = buffer.data(fpp + 87);
    const auto *fpp_88 = buffer.data(fpp + 88);
    const auto *fpp_89 = buffer.data(fpp + 89);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_z, pb_x, pc_x, pc_y, pc_z, dpd0_57, \
                         dpp_38, dpd1_57, fsd0_47, fsp_23, fsd1_47, \
                         fpp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * fsp_23[k]
                   + f_4 * pc_x[k] * fpp_65[k];

        t_129[k] = pa_z[k] * dpd0_57[k]
                   - f_5 * pc_z[k] * dpd1_57[k];

        t_130[k] = f_6 * dpp_38[k]
                   + f_4 * pc_y[k] * fpp_65[k];

        t_131[k] = pb_x[k] * fsd0_47[k]
                   - f_5 * pc_x[k] * fsd1_47[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_z, pc_x, pc_z, dpd0_60, dpd0_63, \
                         dpd1_60, dpd1_63, fpp_67, fpp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_z[k] * dpd0_60[k]
                   - f_5 * pc_z[k] * dpd1_60[k];

        t_133[k] = f_4 * pc_x[k] * fpp_67[k];

        t_134[k] = f_4 * pc_x[k] * fpp_68[k];

        t_135[k] = pa_z[k] * dpd0_63[k]
                   - f_5 * pc_z[k] * dpd1_63[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_x, pc_y, pc_z, dpp_32, dpp_41, fsp_23, \
                         fps0_22, fps0_23, fps1_22, fps1_23, fpp_68, \
                         fpp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_6 * dpp_41[k]
                   + f_1 * fsp_23[k]
                   + f_4 * pc_y[k] * fpp_68[k];

        t_137[k] = f_1 * dpp_32[k]
                   + f_2 * fps0_22[k]
                   - f_3 * fps1_22[k]
                   + f_4 * pc_z[k] * fpp_68[k];

        t_138[k] = f_2 * fps0_23[k]
                   - f_3 * fps1_23[k]
                   + f_4 * pc_x[k] * fpp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, dpp_43, dpp_44, fps0_23, \
                         fps1_23, fpp_70, fpp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_4 * pc_x[k] * fpp_70[k];

        t_140[k] = f_4 * pc_x[k] * fpp_71[k];

        t_141[k] = f_6 * dpp_43[k]
                   + f_2 * fps0_23[k]
                   - f_3 * fps1_23[k]
                   + f_4 * pc_y[k] * fpp_70[k];

        t_142[k] = f_6 * dpp_44[k]
                   + f_4 * pc_y[k] * fpp_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_y, pc_x, pc_y, ppd0_53, ppd1_53, dpd0_89, \
                         dpd0_90, dpd1_89, dpd1_90, fsp_25, fpp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * ppd0_53[k]
                   - f_8 * ppd1_53[k]
                   + pa_y[k] * dpd0_89[k]
                   - f_5 * pc_y[k] * dpd1_89[k];

        t_144[k] = pa_y[k] * dpd0_90[k]
                   - f_5 * pc_y[k] * dpd1_90[k];

        t_145[k] = f_1 * fsp_25[k]
                   + f_4 * pc_x[k] * fpp_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_y, pb_x, pc_x, pc_y, dpd0_95, dpp_47, \
                         dpd1_95, fsd0_51, fsp_26, fsd1_51, fpp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * fsp_26[k]
                   + f_4 * pc_x[k] * fpp_74[k];

        t_147[k] = pb_x[k] * fsd0_51[k]
                   - f_5 * pc_x[k] * fsd1_51[k];

        t_148[k] = f_1 * dpp_47[k]
                   + f_4 * pc_y[k] * fpp_74[k];

        t_149[k] = pa_y[k] * dpd0_95[k]
                   - f_5 * pc_y[k] * dpd1_95[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pc_x, pc_y, dpp_49, dpp_50, \
                         fsp_25, fsp_26, fps0_25, fps1_25, fpp_75, fpp_76, \
                         fpp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_2 * fps0_25[k]
                   - f_3 * fps1_25[k]
                   + f_4 * pc_x[k] * fpp_75[k];

        t_151[k] = f_4 * pc_x[k] * fpp_76[k];

        t_152[k] = f_4 * pc_x[k] * fpp_77[k];

        t_153[k] = f_1 * dpp_49[k]
                   + f_1 * fsp_25[k]
                   + f_2 * fps0_25[k]
                   - f_3 * fps1_25[k]
                   + f_4 * pc_y[k] * fpp_76[k];

        t_154[k] = f_1 * dpp_50[k]
                   + f_1 * fsp_26[k]
                   + f_4 * pc_y[k] * fpp_77[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pc_x, pc_y, pc_z, dpd0_102, dpp_41, \
                         dpd1_102, fps0_25, fps1_25, fpp_77, fpp_79, \
                         fpp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_6 * dpp_41[k]
                   + f_2 * fps0_25[k]
                   - f_3 * fps1_25[k]
                   + f_4 * pc_z[k] * fpp_77[k];

        t_156[k] = pa_y[k] * dpd0_102[k]
                   - f_5 * pc_y[k] * dpd1_102[k];

        t_157[k] = f_4 * pc_x[k] * fpp_79[k];

        t_158[k] = f_4 * pc_x[k] * fpp_80[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pc_y, dpd0_105, dpd0_107, dpp_52, dpp_53, \
                         dpd1_105, dpd1_107, fpp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_y[k] * dpd0_105[k]
                   + f_6 * dpp_52[k]
                   - f_5 * pc_y[k] * dpd1_105[k];

        t_160[k] = f_1 * dpp_53[k]
                   + f_4 * pc_y[k] * fpp_80[k];

        t_161[k] = pa_y[k] * dpd0_107[k]
                   - f_5 * pc_y[k] * dpd1_107[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_x, pc_x, fsd0_54, fsd0_57, fsp_27, \
                         fsp_28, fsp_29, fsd1_54, fsd1_57, fpp_82, \
                         fpp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pb_x[k] * fsd0_54[k]
                   + f_6 * fsp_27[k]
                   - f_5 * pc_x[k] * fsd1_54[k];

        t_163[k] = f_1 * fsp_28[k]
                   + f_4 * pc_x[k] * fpp_82[k];

        t_164[k] = f_1 * fsp_29[k]
                   + f_4 * pc_x[k] * fpp_83[k];

        t_165[k] = pb_x[k] * fsd0_57[k]
                   - f_5 * pc_x[k] * fsd1_57[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pb_x, pb_y, pc_x, pc_y, fsd0_54, \
                         fsd0_59, fsd1_54, fsd1_59, fpp_83, fpp_85, \
                         fpp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_4 * pc_y[k] * fpp_83[k];

        t_167[k] = pb_x[k] * fsd0_59[k]
                   - f_5 * pc_x[k] * fsd1_59[k];

        t_168[k] = pb_y[k] * fsd0_54[k]
                   - f_5 * pc_y[k] * fsd1_54[k];

        t_169[k] = f_4 * pc_x[k] * fpp_85[k];

        t_170[k] = f_4 * pc_x[k] * fpp_86[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_y, pc_y, fsd0_57, fsd0_59, fsp_28, fsp_29, \
                         fsd1_57, fsd1_59, fpp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pb_y[k] * fsd0_57[k]
                   + f_6 * fsp_28[k]
                   - f_5 * pc_y[k] * fsd1_57[k];

        t_172[k] = f_1 * fsp_29[k]
                   + f_4 * pc_y[k] * fpp_86[k];

        t_173[k] = pb_y[k] * fsd0_59[k]
                   - f_5 * pc_y[k] * fsd1_59[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, dpp_53, \
                         fsp_29, fps0_29, fps1_29, fpp_87, fpp_88, \
                         fpp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_2 * fps0_29[k]
                   - f_3 * fps1_29[k]
                   + f_4 * pc_x[k] * fpp_87[k];

        t_175[k] = f_4 * pc_x[k] * fpp_88[k];

        t_176[k] = f_4 * pc_x[k] * fpp_89[k];

        t_177[k] = f_2 * fps0_29[k]
                   - f_3 * fps1_29[k]
                   + f_4 * pc_y[k] * fpp_88[k];

        t_178[k] = f_4 * pc_y[k] * fpp_89[k];

        t_179[k] = f_0 * dpp_53[k]
                   + f_1 * fsp_29[k]
                   + f_2 * fps0_29[k]
                   - f_3 * fps1_29[k]
                   + f_4 * pc_z[k] * fpp_89[k];
    }
}

auto
compute_prim_fpd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppd0,
                                                   const size_t ppd1, const size_t dpd0,
                                                   const size_t dpp, const size_t dpd1,
                                                   const size_t fsd0, const size_t fsp,
                                                   const size_t fsd1, const size_t fps0,
                                                   const size_t fps1, const size_t fpp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fpd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, ppd0,
                                                              ppd1, dpd0, dpp, dpd1, fsd0, fsp,
                                                              fsd1, fps0, fps1, fpp, ncols,
                                                              gamma, p, q);

    compute_prim_fpd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppd0,
                                                              ppd1, dpd0, dpp, dpd1, fsd0, fsp,
                                                              fsd1, fps0, fps1, fpp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
