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


#include "SimdThreeCenterElectronRepulsionVrrRecDPF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dpf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppf0, const size_t ppd,
                                                          const size_t ppf1, const size_t dsf0,
                                                          const size_t dsd, const size_t dsf1,
                                                          const size_t dpp0, const size_t dpp1,
                                                          const size_t dpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);

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
    auto *t_133 = buffer.data(target + 133);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppf0_0 = buffer.data(ppf0 + 0);
    const auto *ppf0_3 = buffer.data(ppf0 + 3);
    const auto *ppf0_5 = buffer.data(ppf0 + 5);
    const auto *ppf0_10 = buffer.data(ppf0 + 10);
    const auto *ppf0_20 = buffer.data(ppf0 + 20);
    const auto *ppf0_31 = buffer.data(ppf0 + 31);
    const auto *ppf0_36 = buffer.data(ppf0 + 36);
    const auto *ppf0_41 = buffer.data(ppf0 + 41);
    const auto *ppf0_46 = buffer.data(ppf0 + 46);
    const auto *ppf0_48 = buffer.data(ppf0 + 48);
    const auto *ppf0_49 = buffer.data(ppf0 + 49);
    const auto *ppf0_56 = buffer.data(ppf0 + 56);
    const auto *ppf0_59 = buffer.data(ppf0 + 59);
    const auto *ppf0_60 = buffer.data(ppf0 + 60);
    const auto *ppf0_62 = buffer.data(ppf0 + 62);
    const auto *ppf0_69 = buffer.data(ppf0 + 69);
    const auto *ppf0_76 = buffer.data(ppf0 + 76);
    const auto *ppf0_77 = buffer.data(ppf0 + 77);
    const auto *ppf0_79 = buffer.data(ppf0 + 79);
    const auto *ppf0_86 = buffer.data(ppf0 + 86);
    const auto *ppf0_87 = buffer.data(ppf0 + 87);
    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppd_0 = buffer.data(ppd + 0);
    const auto *ppd_3 = buffer.data(ppd + 3);
    const auto *ppd_5 = buffer.data(ppd + 5);
    const auto *ppd_6 = buffer.data(ppd + 6);
    const auto *ppd_9 = buffer.data(ppd + 9);
    const auto *ppd_11 = buffer.data(ppd + 11);
    const auto *ppd_12 = buffer.data(ppd + 12);
    const auto *ppd_15 = buffer.data(ppd + 15);
    const auto *ppd_17 = buffer.data(ppd + 17);
    const auto *ppd_21 = buffer.data(ppd + 21);
    const auto *ppd_23 = buffer.data(ppd + 23);
    const auto *ppd_24 = buffer.data(ppd + 24);
    const auto *ppd_27 = buffer.data(ppd + 27);
    const auto *ppd_29 = buffer.data(ppd + 29);
    const auto *ppd_33 = buffer.data(ppd + 33);
    const auto *ppd_35 = buffer.data(ppd + 35);
    const auto *ppd_41 = buffer.data(ppd + 41);
    const auto *ppd_45 = buffer.data(ppd + 45);
    const auto *ppd_47 = buffer.data(ppd + 47);
    const auto *ppd_48 = buffer.data(ppd + 48);
    const auto *ppd_51 = buffer.data(ppd + 51);
    const auto *ppd_53 = buffer.data(ppd + 53);

    const auto *ppf1_0 = buffer.data(ppf1 + 0);
    const auto *ppf1_3 = buffer.data(ppf1 + 3);
    const auto *ppf1_5 = buffer.data(ppf1 + 5);
    const auto *ppf1_10 = buffer.data(ppf1 + 10);
    const auto *ppf1_20 = buffer.data(ppf1 + 20);
    const auto *ppf1_31 = buffer.data(ppf1 + 31);
    const auto *ppf1_36 = buffer.data(ppf1 + 36);
    const auto *ppf1_41 = buffer.data(ppf1 + 41);
    const auto *ppf1_46 = buffer.data(ppf1 + 46);
    const auto *ppf1_48 = buffer.data(ppf1 + 48);
    const auto *ppf1_49 = buffer.data(ppf1 + 49);
    const auto *ppf1_56 = buffer.data(ppf1 + 56);
    const auto *ppf1_59 = buffer.data(ppf1 + 59);
    const auto *ppf1_60 = buffer.data(ppf1 + 60);
    const auto *ppf1_62 = buffer.data(ppf1 + 62);
    const auto *ppf1_69 = buffer.data(ppf1 + 69);
    const auto *ppf1_76 = buffer.data(ppf1 + 76);
    const auto *ppf1_77 = buffer.data(ppf1 + 77);
    const auto *ppf1_79 = buffer.data(ppf1 + 79);
    const auto *ppf1_86 = buffer.data(ppf1 + 86);
    const auto *ppf1_87 = buffer.data(ppf1 + 87);
    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *dsf0_0 = buffer.data(dsf0 + 0);
    const auto *dsf0_6 = buffer.data(dsf0 + 6);
    const auto *dsf0_9 = buffer.data(dsf0 + 9);
    const auto *dsf0_30 = buffer.data(dsf0 + 30);
    const auto *dsf0_31 = buffer.data(dsf0 + 31);
    const auto *dsf0_36 = buffer.data(dsf0 + 36);
    const auto *dsf0_39 = buffer.data(dsf0 + 39);

    const auto *dsd_0 = buffer.data(dsd + 0);
    const auto *dsd_2 = buffer.data(dsd + 2);
    const auto *dsd_3 = buffer.data(dsd + 3);
    const auto *dsd_5 = buffer.data(dsd + 5);
    const auto *dsd_6 = buffer.data(dsd + 6);
    const auto *dsd_7 = buffer.data(dsd + 7);
    const auto *dsd_9 = buffer.data(dsd + 9);
    const auto *dsd_12 = buffer.data(dsd + 12);
    const auto *dsd_14 = buffer.data(dsd + 14);
    const auto *dsd_17 = buffer.data(dsd + 17);
    const auto *dsd_18 = buffer.data(dsd + 18);
    const auto *dsd_19 = buffer.data(dsd + 19);
    const auto *dsd_21 = buffer.data(dsd + 21);
    const auto *dsd_22 = buffer.data(dsd + 22);
    const auto *dsd_23 = buffer.data(dsd + 23);
    const auto *dsd_27 = buffer.data(dsd + 27);
    const auto *dsd_28 = buffer.data(dsd + 28);
    const auto *dsd_29 = buffer.data(dsd + 29);

    const auto *dsf1_0 = buffer.data(dsf1 + 0);
    const auto *dsf1_6 = buffer.data(dsf1 + 6);
    const auto *dsf1_9 = buffer.data(dsf1 + 9);
    const auto *dsf1_30 = buffer.data(dsf1 + 30);
    const auto *dsf1_31 = buffer.data(dsf1 + 31);
    const auto *dsf1_36 = buffer.data(dsf1 + 36);
    const auto *dsf1_39 = buffer.data(dsf1 + 39);

    const auto *dpp0_0 = buffer.data(dpp0 + 0);
    const auto *dpp0_1 = buffer.data(dpp0 + 1);
    const auto *dpp0_2 = buffer.data(dpp0 + 2);
    const auto *dpp0_10 = buffer.data(dpp0 + 10);
    const auto *dpp0_11 = buffer.data(dpp0 + 11);
    const auto *dpp0_12 = buffer.data(dpp0 + 12);
    const auto *dpp0_19 = buffer.data(dpp0 + 19);
    const auto *dpp0_20 = buffer.data(dpp0 + 20);
    const auto *dpp0_24 = buffer.data(dpp0 + 24);
    const auto *dpp0_30 = buffer.data(dpp0 + 30);
    const auto *dpp0_31 = buffer.data(dpp0 + 31);
    const auto *dpp0_32 = buffer.data(dpp0 + 32);
    const auto *dpp0_39 = buffer.data(dpp0 + 39);
    const auto *dpp0_41 = buffer.data(dpp0 + 41);

    const auto *dpp1_0 = buffer.data(dpp1 + 0);
    const auto *dpp1_1 = buffer.data(dpp1 + 1);
    const auto *dpp1_2 = buffer.data(dpp1 + 2);
    const auto *dpp1_10 = buffer.data(dpp1 + 10);
    const auto *dpp1_11 = buffer.data(dpp1 + 11);
    const auto *dpp1_12 = buffer.data(dpp1 + 12);
    const auto *dpp1_19 = buffer.data(dpp1 + 19);
    const auto *dpp1_20 = buffer.data(dpp1 + 20);
    const auto *dpp1_24 = buffer.data(dpp1 + 24);
    const auto *dpp1_30 = buffer.data(dpp1 + 30);
    const auto *dpp1_31 = buffer.data(dpp1 + 31);
    const auto *dpp1_32 = buffer.data(dpp1 + 32);
    const auto *dpp1_39 = buffer.data(dpp1 + 39);
    const auto *dpp1_41 = buffer.data(dpp1 + 41);

    const auto *dpd_0 = buffer.data(dpd + 0);
    const auto *dpd_2 = buffer.data(dpd + 2);
    const auto *dpd_3 = buffer.data(dpd + 3);
    const auto *dpd_5 = buffer.data(dpd + 5);
    const auto *dpd_6 = buffer.data(dpd + 6);
    const auto *dpd_8 = buffer.data(dpd + 8);
    const auto *dpd_9 = buffer.data(dpd + 9);
    const auto *dpd_11 = buffer.data(dpd + 11);
    const auto *dpd_12 = buffer.data(dpd + 12);
    const auto *dpd_14 = buffer.data(dpd + 14);
    const auto *dpd_15 = buffer.data(dpd + 15);
    const auto *dpd_17 = buffer.data(dpd + 17);
    const auto *dpd_18 = buffer.data(dpd + 18);
    const auto *dpd_19 = buffer.data(dpd + 19);
    const auto *dpd_21 = buffer.data(dpd + 21);
    const auto *dpd_23 = buffer.data(dpd + 23);
    const auto *dpd_24 = buffer.data(dpd + 24);
    const auto *dpd_25 = buffer.data(dpd + 25);
    const auto *dpd_27 = buffer.data(dpd + 27);
    const auto *dpd_29 = buffer.data(dpd + 29);
    const auto *dpd_30 = buffer.data(dpd + 30);
    const auto *dpd_31 = buffer.data(dpd + 31);
    const auto *dpd_33 = buffer.data(dpd + 33);
    const auto *dpd_35 = buffer.data(dpd + 35);
    const auto *dpd_36 = buffer.data(dpd + 36);
    const auto *dpd_38 = buffer.data(dpd + 38);
    const auto *dpd_39 = buffer.data(dpd + 39);
    const auto *dpd_40 = buffer.data(dpd + 40);
    const auto *dpd_41 = buffer.data(dpd + 41);
    const auto *dpd_42 = buffer.data(dpd + 42);
    const auto *dpd_44 = buffer.data(dpd + 44);
    const auto *dpd_45 = buffer.data(dpd + 45);
    const auto *dpd_47 = buffer.data(dpd + 47);
    const auto *dpd_48 = buffer.data(dpd + 48);
    const auto *dpd_50 = buffer.data(dpd + 50);
    const auto *dpd_51 = buffer.data(dpd + 51);
    const auto *dpd_53 = buffer.data(dpd + 53);
    const auto *dpd_54 = buffer.data(dpd + 54);
    const auto *dpd_57 = buffer.data(dpd + 57);
    const auto *dpd_58 = buffer.data(dpd + 58);
    const auto *dpd_59 = buffer.data(dpd + 59);
    const auto *dpd_60 = buffer.data(dpd + 60);
    const auto *dpd_61 = buffer.data(dpd + 61);
    const auto *dpd_63 = buffer.data(dpd + 63);
    const auto *dpd_64 = buffer.data(dpd + 64);
    const auto *dpd_65 = buffer.data(dpd + 65);
    const auto *dpd_66 = buffer.data(dpd + 66);
    const auto *dpd_69 = buffer.data(dpd + 69);
    const auto *dpd_70 = buffer.data(dpd + 70);
    const auto *dpd_71 = buffer.data(dpd + 71);
    const auto *dpd_75 = buffer.data(dpd + 75);
    const auto *dpd_76 = buffer.data(dpd + 76);
    const auto *dpd_77 = buffer.data(dpd + 77);
    const auto *dpd_78 = buffer.data(dpd + 78);
    const auto *dpd_80 = buffer.data(dpd + 80);
    const auto *dpd_81 = buffer.data(dpd + 81);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, ppd_0, ppd_3, dsd_0, dsd_3, \
                         dpp0_0, dpp1_0, dpd_0, dpd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ppd_0[k]
                 + f_1 * dsd_0[k]
                 + f_2 * dpp0_0[k]
                 - f_3 * dpp1_0[k]
                 + f_4 * pc_x[k] * dpd_0[k];

        t_1[k] = f_4 * pc_y[k] * dpd_0[k];

        t_2[k] = f_4 * pc_z[k] * dpd_0[k];

        t_3[k] = f_0 * ppd_3[k]
                 + f_1 * dsd_3[k]
                 + f_4 * pc_x[k] * dpd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, ppd_5, dsd_5, dpp0_1, \
                         dpp1_1, dpd_2, dpd_3, dpd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * pc_y[k] * dpd_2[k];

        t_5[k] = f_0 * ppd_5[k]
                 + f_1 * dsd_5[k]
                 + f_4 * pc_x[k] * dpd_5[k];

        t_6[k] = f_2 * dpp0_1[k]
                 - f_3 * dpp1_1[k]
                 + f_4 * pc_y[k] * dpd_3[k];

        t_7[k] = f_4 * pc_z[k] * dpd_3[k];

        t_8[k] = f_4 * pc_y[k] * dpd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, dsf0_0, dsd_0, dsf1_0, \
                         dpp0_2, dpp1_2, dpd_5, dpd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * dpp0_2[k]
                 - f_3 * dpp1_2[k]
                 + f_4 * pc_z[k] * dpd_5[k];

        t_10[k] = pb_y[k] * dsf0_0[k]
                  - f_5 * pc_y[k] * dsf1_0[k];

        t_11[k] = f_1 * dsd_0[k]
                  + f_4 * pc_y[k] * dpd_6[k];

        t_12[k] = f_4 * pc_z[k] * dpd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pc_x, pc_y, ppd_9, ppd_11, dsf0_6, \
                         dsd_2, dsd_3, dsf1_6, dpd_8, dpd_9, dpd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * ppd_9[k]
                  + f_4 * pc_x[k] * dpd_9[k];

        t_14[k] = f_1 * dsd_2[k]
                  + f_4 * pc_y[k] * dpd_8[k];

        t_15[k] = f_0 * ppd_11[k]
                  + f_4 * pc_x[k] * dpd_11[k];

        t_16[k] = pb_y[k] * dsf0_6[k]
                  + f_6 * dsd_3[k]
                  - f_5 * pc_y[k] * dsf1_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, pc_y, pc_z, dsf0_0, dsf0_9, \
                         dsd_5, dsf1_0, dsf1_9, dpd_9, dpd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * pc_z[k] * dpd_9[k];

        t_18[k] = f_1 * dsd_5[k]
                  + f_4 * pc_y[k] * dpd_11[k];

        t_19[k] = pb_y[k] * dsf0_9[k]
                  - f_5 * pc_y[k] * dsf1_9[k];

        t_20[k] = pb_z[k] * dsf0_0[k]
                  - f_5 * pc_z[k] * dsf1_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, ppd_15, ppd_17, \
                         dsd_0, dpd_12, dpd_14, dpd_15, dpd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * pc_y[k] * dpd_12[k];

        t_22[k] = f_1 * dsd_0[k]
                  + f_4 * pc_z[k] * dpd_12[k];

        t_23[k] = f_0 * ppd_15[k]
                  + f_4 * pc_x[k] * dpd_15[k];

        t_24[k] = f_4 * pc_y[k] * dpd_14[k];

        t_25[k] = f_0 * ppd_17[k]
                  + f_4 * pc_x[k] * dpd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, dsf0_6, dsf0_9, dsd_3, \
                         dsd_5, dsf1_6, dsf1_9, dpd_15, dpd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * dsf0_6[k]
                  - f_5 * pc_z[k] * dsf1_6[k];

        t_27[k] = f_1 * dsd_3[k]
                  + f_4 * pc_z[k] * dpd_15[k];

        t_28[k] = f_4 * pc_y[k] * dpd_17[k];

        t_29[k] = pb_z[k] * dsf0_9[k]
                  + f_6 * dsd_5[k]
                  - f_5 * pc_z[k] * dsf1_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_x, pc_y, pc_z, ppf0_0, ppd_0, \
                         ppd_21, ppf1_0, dsd_9, dpd_18, dpd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * ppf0_0[k]
                  - f_5 * pc_y[k] * ppf1_0[k];

        t_31[k] = f_1 * ppd_0[k]
                  + f_4 * pc_y[k] * dpd_18[k];

        t_32[k] = f_4 * pc_z[k] * dpd_18[k];

        t_33[k] = f_1 * ppd_21[k]
                  + f_1 * dsd_9[k]
                  + f_4 * pc_x[k] * dpd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, ppf0_5, ppd_3, ppf1_5, \
                         dpp0_10, dpp1_10, dpd_19, dpd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_4 * pc_z[k] * dpd_19[k];

        t_35[k] = pa_y[k] * ppf0_5[k]
                  - f_5 * pc_y[k] * ppf1_5[k];

        t_36[k] = f_1 * ppd_3[k]
                  + f_2 * dpp0_10[k]
                  - f_3 * dpp1_10[k]
                  + f_4 * pc_y[k] * dpd_21[k];

        t_37[k] = f_4 * pc_z[k] * dpd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pc_x, pc_y, pc_z, ppd_5, ppd_24, dpp0_11, dpp0_12, \
                         dpp1_11, dpp1_12, dpd_23, dpd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ppd_5[k]
                  + f_4 * pc_y[k] * dpd_23[k];

        t_39[k] = f_2 * dpp0_11[k]
                  - f_3 * dpp1_11[k]
                  + f_4 * pc_z[k] * dpd_23[k];

        t_40[k] = f_1 * ppd_24[k]
                  + f_2 * dpp0_12[k]
                  - f_3 * dpp1_12[k]
                  + f_4 * pc_x[k] * dpd_24[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pc_x, pc_y, pc_z, ppd_6, ppd_27, \
                         ppd_29, dsd_6, dpd_24, dpd_25, dpd_27, \
                         dpd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ppd_6[k]
                  + f_1 * dsd_6[k]
                  + f_4 * pc_y[k] * dpd_24[k];

        t_42[k] = f_4 * pc_z[k] * dpd_24[k];

        t_43[k] = f_1 * ppd_27[k]
                  + f_4 * pc_x[k] * dpd_27[k];

        t_44[k] = f_4 * pc_z[k] * dpd_25[k];

        t_45[k] = f_1 * ppd_29[k]
                  + f_4 * pc_x[k] * dpd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pc_x, pc_z, ppf0_46, ppf0_48, ppf0_49, \
                         ppf1_46, ppf1_48, ppf1_49, dpd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * ppf0_46[k]
                  - f_5 * pc_x[k] * ppf1_46[k];

        t_47[k] = f_4 * pc_z[k] * dpd_27[k];

        t_48[k] = pa_x[k] * ppf0_48[k]
                  - f_5 * pc_x[k] * ppf1_48[k];

        t_49[k] = pa_x[k] * ppf0_49[k]
                  - f_5 * pc_x[k] * ppf1_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, pc_z, ppf0_20, ppd_12, \
                         ppd_33, ppf1_20, dsd_6, dpd_30, dpd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * ppf0_20[k]
                  - f_5 * pc_y[k] * ppf1_20[k];

        t_51[k] = f_1 * ppd_12[k]
                  + f_4 * pc_y[k] * dpd_30[k];

        t_52[k] = f_1 * dsd_6[k]
                  + f_4 * pc_z[k] * dpd_30[k];

        t_53[k] = f_1 * ppd_33[k]
                  + f_4 * pc_x[k] * dpd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pc_x, pc_z, ppf0_56, ppd_35, ppf1_56, \
                         dsd_7, dsd_9, dpd_31, dpd_33, dpd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * dsd_7[k]
                  + f_4 * pc_z[k] * dpd_31[k];

        t_55[k] = f_1 * ppd_35[k]
                  + f_4 * pc_x[k] * dpd_35[k];

        t_56[k] = pa_x[k] * ppf0_56[k]
                  - f_5 * pc_x[k] * ppf1_56[k];

        t_57[k] = f_1 * dsd_9[k]
                  + f_4 * pc_z[k] * dpd_33[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pa_z, pc_x, pc_y, pc_z, ppf0_0, \
                         ppf0_59, ppd_17, ppf1_0, ppf1_59, dpd_35, \
                         dpd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * ppd_17[k]
                  + f_4 * pc_y[k] * dpd_35[k];

        t_59[k] = pa_x[k] * ppf0_59[k]
                  - f_5 * pc_x[k] * ppf1_59[k];

        t_60[k] = pa_z[k] * ppf0_0[k]
                  - f_5 * pc_z[k] * ppf1_0[k];

        t_61[k] = f_4 * pc_y[k] * dpd_36[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_x, pc_y, pc_z, ppf0_3, ppd_0, \
                         ppd_41, ppf1_3, dsd_17, dpd_36, dpd_38, \
                         dpd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_1 * ppd_0[k]
                  + f_4 * pc_z[k] * dpd_36[k];

        t_63[k] = pa_z[k] * ppf0_3[k]
                  - f_5 * pc_z[k] * ppf1_3[k];

        t_64[k] = f_4 * pc_y[k] * dpd_38[k];

        t_65[k] = f_1 * ppd_41[k]
                  + f_1 * dsd_17[k]
                  + f_4 * pc_x[k] * dpd_41[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_y, pc_z, ppd_5, dpp0_19, dpp0_20, dpp1_19, \
                         dpp1_20, dpd_39, dpd_40, dpd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_2 * dpp0_19[k]
                  - f_3 * dpp1_19[k]
                  + f_4 * pc_y[k] * dpd_39[k];

        t_67[k] = f_7 * dpp0_20[k]
                  - f_8 * dpp1_20[k]
                  + f_4 * pc_y[k] * dpd_40[k];

        t_68[k] = f_4 * pc_y[k] * dpd_41[k];

        t_69[k] = f_1 * ppd_5[k]
                  + f_2 * dpp0_20[k]
                  - f_3 * dpp1_20[k]
                  + f_4 * pc_z[k] * dpd_41[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_y, pc_z, ppf0_10, ppd_6, \
                         ppd_45, ppf1_10, dsd_12, dpd_42, dpd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * ppf0_10[k]
                  - f_5 * pc_z[k] * ppf1_10[k];

        t_71[k] = f_1 * dsd_12[k]
                  + f_4 * pc_y[k] * dpd_42[k];

        t_72[k] = f_1 * ppd_6[k]
                  + f_4 * pc_z[k] * dpd_42[k];

        t_73[k] = f_1 * ppd_45[k]
                  + f_4 * pc_x[k] * dpd_45[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_y, ppf0_76, ppf0_77, ppd_47, \
                         ppf1_76, ppf1_77, dsd_14, dpd_44, dpd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * dsd_14[k]
                  + f_4 * pc_y[k] * dpd_44[k];

        t_75[k] = f_1 * ppd_47[k]
                  + f_4 * pc_x[k] * dpd_47[k];

        t_76[k] = pa_x[k] * ppf0_76[k]
                  - f_5 * pc_x[k] * ppf1_76[k];

        t_77[k] = pa_x[k] * ppf0_77[k]
                  - f_5 * pc_x[k] * ppf1_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pc_x, pc_y, ppf0_79, ppd_48, ppf1_79, \
                         dsd_17, dpp0_24, dpp1_24, dpd_47, dpd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * dsd_17[k]
                  + f_4 * pc_y[k] * dpd_47[k];

        t_79[k] = pa_x[k] * ppf0_79[k]
                  - f_5 * pc_x[k] * ppf1_79[k];

        t_80[k] = f_1 * ppd_48[k]
                  + f_2 * dpp0_24[k]
                  - f_3 * dpp1_24[k]
                  + f_4 * pc_x[k] * dpd_48[k];

        t_81[k] = f_4 * pc_y[k] * dpd_48[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pc_x, pc_y, pc_z, ppd_12, ppd_51, ppd_53, \
                         dsd_12, dpd_48, dpd_50, dpd_51, dpd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_1 * ppd_12[k]
                  + f_1 * dsd_12[k]
                  + f_4 * pc_z[k] * dpd_48[k];

        t_83[k] = f_1 * ppd_51[k]
                  + f_4 * pc_x[k] * dpd_51[k];

        t_84[k] = f_4 * pc_y[k] * dpd_50[k];

        t_85[k] = f_1 * ppd_53[k]
                  + f_4 * pc_x[k] * dpd_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pc_x, pc_y, ppf0_86, ppf0_87, ppf0_89, \
                         ppf1_86, ppf1_87, ppf1_89, dpd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * ppf0_86[k]
                  - f_5 * pc_x[k] * ppf1_86[k];

        t_87[k] = pa_x[k] * ppf0_87[k]
                  - f_5 * pc_x[k] * ppf1_87[k];

        t_88[k] = f_4 * pc_y[k] * dpd_53[k];

        t_89[k] = pa_x[k] * ppf0_89[k]
                  - f_5 * pc_x[k] * ppf1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pc_x, pc_z, dsf0_30, dsf0_31, dsd_18, \
                         dsd_19, dsd_21, dsf1_30, dsf1_31, dpd_54, \
                         dpd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_x[k] * dsf0_30[k]
                  + f_6 * dsd_18[k]
                  - f_5 * pc_x[k] * dsf1_30[k];

        t_91[k] = pb_x[k] * dsf0_31[k]
                  + f_0 * dsd_19[k]
                  - f_5 * pc_x[k] * dsf1_31[k];

        t_92[k] = f_4 * pc_z[k] * dpd_54[k];

        t_93[k] = f_1 * dsd_21[k]
                  + f_4 * pc_x[k] * dpd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_x, pc_x, pc_z, dsf0_36, dsd_22, dsd_23, \
                         dsf1_36, dpd_57, dpd_58, dpd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_1 * dsd_22[k]
                  + f_4 * pc_x[k] * dpd_58[k];

        t_95[k] = f_1 * dsd_23[k]
                  + f_4 * pc_x[k] * dpd_59[k];

        t_96[k] = pb_x[k] * dsf0_36[k]
                  - f_5 * pc_x[k] * dsf1_36[k];

        t_97[k] = f_4 * pc_z[k] * dpd_57[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pb_x, pc_x, pc_y, ppd_23, dsf0_39, dsf1_39, \
                         dpp0_30, dpp1_30, dpd_59, dpd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * ppd_23[k]
                  + f_4 * pc_y[k] * dpd_59[k];

        t_99[k] = pb_x[k] * dsf0_39[k]
                  - f_5 * pc_x[k] * dsf1_39[k];

        t_100[k] = f_2 * dpp0_30[k]
                   - f_3 * dpp1_30[k]
                   + f_4 * pc_x[k] * dpd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_z, dpp0_31, dpp1_31, \
                         dpd_60, dpd_61, dpd_63, dpd_64, dpd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_7 * dpp0_31[k]
                   - f_8 * dpp1_31[k]
                   + f_4 * pc_x[k] * dpd_61[k];

        t_102[k] = f_4 * pc_z[k] * dpd_60[k];

        t_103[k] = f_4 * pc_x[k] * dpd_63[k];

        t_104[k] = f_4 * pc_x[k] * dpd_64[k];

        t_105[k] = f_4 * pc_x[k] * dpd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, ppd_27, ppd_29, dsd_21, \
                         dsd_23, dpp0_31, dpp0_32, dpp1_31, dpp1_32, dpd_63, \
                         dpd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * ppd_27[k]
                   + f_1 * dsd_21[k]
                   + f_2 * dpp0_31[k]
                   - f_3 * dpp1_31[k]
                   + f_4 * pc_y[k] * dpd_63[k];

        t_107[k] = f_4 * pc_z[k] * dpd_63[k];

        t_108[k] = f_0 * ppd_29[k]
                   + f_1 * dsd_23[k]
                   + f_4 * pc_y[k] * dpd_65[k];

        t_109[k] = f_2 * dpp0_32[k]
                   - f_3 * dpp1_32[k]
                   + f_4 * pc_z[k] * dpd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pb_z, pc_x, pc_z, dsf0_30, \
                         dsf0_31, dsd_18, dsf1_30, dsf1_31, dpd_66, dpd_69, \
                         dpd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_z[k] * dsf0_30[k]
                   - f_5 * pc_z[k] * dsf1_30[k];

        t_111[k] = pb_z[k] * dsf0_31[k]
                   - f_5 * pc_z[k] * dsf1_31[k];

        t_112[k] = f_1 * dsd_18[k]
                   + f_4 * pc_z[k] * dpd_66[k];

        t_113[k] = f_4 * pc_x[k] * dpd_69[k];

        t_114[k] = f_4 * pc_x[k] * dpd_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, ppd_35, dsf0_36, \
                         dsd_21, dsf1_36, dpd_69, dpd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_4 * pc_x[k] * dpd_71[k];

        t_116[k] = pb_z[k] * dsf0_36[k]
                   - f_5 * pc_z[k] * dsf1_36[k];

        t_117[k] = f_1 * dsd_21[k]
                   + f_4 * pc_z[k] * dpd_69[k];

        t_118[k] = f_0 * ppd_35[k]
                   + f_4 * pc_y[k] * dpd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pa_y, pa_z, pb_z, pc_y, pc_z, ppf0_31, ppf0_60, \
                         ppf1_31, ppf1_60, dsf0_39, dsd_23, dsf1_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_z[k] * dsf0_39[k]
                   + f_6 * dsd_23[k]
                   - f_5 * pc_z[k] * dsf1_39[k];

        t_120[k] = pa_y[k] * ppf0_60[k]
                   - f_5 * pc_y[k] * ppf1_60[k];

        t_121[k] = pa_z[k] * ppf0_31[k]
                   - f_5 * pc_z[k] * ppf1_31[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_x, pc_y, ppf0_62, ppf1_62, \
                         dsd_27, dsd_28, dsd_29, dpd_75, dpd_76, \
                         dpd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = pa_y[k] * ppf0_62[k]
                   - f_5 * pc_y[k] * ppf1_62[k];

        t_123[k] = f_1 * dsd_27[k]
                   + f_4 * pc_x[k] * dpd_75[k];

        t_124[k] = f_1 * dsd_28[k]
                   + f_4 * pc_x[k] * dpd_76[k];

        t_125[k] = f_1 * dsd_29[k]
                   + f_4 * pc_x[k] * dpd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_y, pa_z, pc_y, pc_z, ppf0_36, ppf0_69, \
                         ppd_21, ppd_41, ppf1_36, ppf1_69, dpd_75, \
                         dpd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * ppf0_36[k]
                   - f_5 * pc_z[k] * ppf1_36[k];

        t_127[k] = f_1 * ppd_21[k]
                   + f_4 * pc_z[k] * dpd_75[k];

        t_128[k] = f_1 * ppd_41[k]
                   + f_4 * pc_y[k] * dpd_77[k];

        t_129[k] = pa_y[k] * ppf0_69[k]
                   - f_5 * pc_y[k] * ppf1_69[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_z, pc_x, pc_z, ppf0_41, ppf1_41, \
                         dpp0_39, dpp0_41, dpp1_39, dpp1_41, dpd_78, dpd_80, \
                         dpd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * dpp0_39[k]
                   - f_3 * dpp1_39[k]
                   + f_4 * pc_x[k] * dpd_78[k];

        t_131[k] = pa_z[k] * ppf0_41[k]
                   - f_5 * pc_z[k] * ppf1_41[k];

        t_132[k] = f_7 * dpp0_41[k]
                   - f_8 * dpp1_41[k]
                   + f_4 * pc_x[k] * dpd_80[k];

        t_133[k] = f_4 * pc_x[k] * dpd_81[k];
    }
}

static auto
compute_prim_dpf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppf0, const size_t ppd,
                                                          const size_t ppf1, const size_t dsf0,
                                                          const size_t dsd, const size_t dsf1,
                                                          const size_t dpp0, const size_t dpp1,
                                                          const size_t dpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);

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

    const auto *ppf0_46 = buffer.data(ppf0 + 46);
    const auto *ppf0_80 = buffer.data(ppf0 + 80);
    const auto *ppf0_82 = buffer.data(ppf0 + 82);
    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppd_27 = buffer.data(ppd + 27);
    const auto *ppd_29 = buffer.data(ppd + 29);
    const auto *ppd_33 = buffer.data(ppd + 33);
    const auto *ppd_47 = buffer.data(ppd + 47);
    const auto *ppd_51 = buffer.data(ppd + 51);
    const auto *ppd_53 = buffer.data(ppd + 53);

    const auto *ppf1_46 = buffer.data(ppf1 + 46);
    const auto *ppf1_80 = buffer.data(ppf1 + 80);
    const auto *ppf1_82 = buffer.data(ppf1 + 82);
    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *dsf0_50 = buffer.data(dsf0 + 50);
    const auto *dsf0_52 = buffer.data(dsf0 + 52);
    const auto *dsf0_56 = buffer.data(dsf0 + 56);
    const auto *dsf0_57 = buffer.data(dsf0 + 57);
    const auto *dsf0_59 = buffer.data(dsf0 + 59);

    const auto *dsd_27 = buffer.data(dsd + 27);
    const auto *dsd_29 = buffer.data(dsd + 29);
    const auto *dsd_30 = buffer.data(dsd + 30);
    const auto *dsd_32 = buffer.data(dsd + 32);
    const auto *dsd_33 = buffer.data(dsd + 33);
    const auto *dsd_34 = buffer.data(dsd + 34);
    const auto *dsd_35 = buffer.data(dsd + 35);

    const auto *dsf1_50 = buffer.data(dsf1 + 50);
    const auto *dsf1_52 = buffer.data(dsf1 + 52);
    const auto *dsf1_56 = buffer.data(dsf1 + 56);
    const auto *dsf1_57 = buffer.data(dsf1 + 57);
    const auto *dsf1_59 = buffer.data(dsf1 + 59);

    const auto *dpp0_41 = buffer.data(dpp0 + 41);
    const auto *dpp0_43 = buffer.data(dpp0 + 43);
    const auto *dpp0_51 = buffer.data(dpp0 + 51);
    const auto *dpp0_52 = buffer.data(dpp0 + 52);
    const auto *dpp0_53 = buffer.data(dpp0 + 53);

    const auto *dpp1_41 = buffer.data(dpp1 + 41);
    const auto *dpp1_43 = buffer.data(dpp1 + 43);
    const auto *dpp1_51 = buffer.data(dpp1 + 51);
    const auto *dpp1_52 = buffer.data(dpp1 + 52);
    const auto *dpp1_53 = buffer.data(dpp1 + 53);

    const auto *dpd_81 = buffer.data(dpd + 81);
    const auto *dpd_82 = buffer.data(dpd + 82);
    const auto *dpd_83 = buffer.data(dpd + 83);
    const auto *dpd_85 = buffer.data(dpd + 85);
    const auto *dpd_87 = buffer.data(dpd + 87);
    const auto *dpd_88 = buffer.data(dpd + 88);
    const auto *dpd_89 = buffer.data(dpd + 89);
    const auto *dpd_90 = buffer.data(dpd + 90);
    const auto *dpd_93 = buffer.data(dpd + 93);
    const auto *dpd_94 = buffer.data(dpd + 94);
    const auto *dpd_95 = buffer.data(dpd + 95);
    const auto *dpd_96 = buffer.data(dpd + 96);
    const auto *dpd_99 = buffer.data(dpd + 99);
    const auto *dpd_100 = buffer.data(dpd + 100);
    const auto *dpd_101 = buffer.data(dpd + 101);
    const auto *dpd_102 = buffer.data(dpd + 102);
    const auto *dpd_104 = buffer.data(dpd + 104);
    const auto *dpd_105 = buffer.data(dpd + 105);
    const auto *dpd_106 = buffer.data(dpd + 106);
    const auto *dpd_107 = buffer.data(dpd + 107);

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_z, pc_x, pc_z, ppf0_46, ppd_27, \
                         ppf1_46, dpd_81, dpd_82, dpd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_4 * pc_x[k] * dpd_82[k];

        t_135[k] = f_4 * pc_x[k] * dpd_83[k];

        t_136[k] = pa_z[k] * ppf0_46[k]
                   - f_5 * pc_z[k] * ppf1_46[k];

        t_137[k] = f_1 * ppd_27[k]
                   + f_4 * pc_z[k] * dpd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pa_y, pc_y, pc_z, ppf0_80, ppd_29, ppd_47, \
                         ppf1_80, dsd_29, dpp0_41, dpp1_41, dpd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_1 * ppd_47[k]
                   + f_1 * dsd_29[k]
                   + f_4 * pc_y[k] * dpd_83[k];

        t_139[k] = f_1 * ppd_29[k]
                   + f_2 * dpp0_41[k]
                   - f_3 * dpp1_41[k]
                   + f_4 * pc_z[k] * dpd_83[k];

        t_140[k] = pa_y[k] * ppf0_80[k]
                   - f_5 * pc_y[k] * ppf1_80[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pa_y, pc_x, pc_y, ppf0_82, \
                         ppf1_82, dpp0_43, dpp1_43, dpd_85, dpd_87, dpd_88, \
                         dpd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_7 * dpp0_43[k]
                   - f_8 * dpp1_43[k]
                   + f_4 * pc_x[k] * dpd_85[k];

        t_142[k] = pa_y[k] * ppf0_82[k]
                   - f_5 * pc_y[k] * ppf1_82[k];

        t_143[k] = f_4 * pc_x[k] * dpd_87[k];

        t_144[k] = f_4 * pc_x[k] * dpd_88[k];

        t_145[k] = f_4 * pc_x[k] * dpd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, pc_z, ppd_33, ppd_51, ppd_53, dsd_27, \
                         dpp0_43, dpp1_43, dpd_87, dpd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * ppd_51[k]
                   + f_2 * dpp0_43[k]
                   - f_3 * dpp1_43[k]
                   + f_4 * pc_y[k] * dpd_87[k];

        t_147[k] = f_1 * ppd_33[k]
                   + f_1 * dsd_27[k]
                   + f_4 * pc_z[k] * dpd_87[k];

        t_148[k] = f_1 * ppd_53[k]
                   + f_4 * pc_y[k] * dpd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_y, pb_x, pc_x, pc_y, ppf0_89, ppf1_89, \
                         dsf0_50, dsd_30, dsf1_50, dpd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * ppf0_89[k]
                   - f_5 * pc_y[k] * ppf1_89[k];

        t_150[k] = pb_x[k] * dsf0_50[k]
                   + f_6 * dsd_30[k]
                   - f_5 * pc_x[k] * dsf1_50[k];

        t_151[k] = f_4 * pc_y[k] * dpd_90[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pc_x, dsf0_52, dsd_32, dsd_33, \
                         dsd_34, dsd_35, dsf1_52, dpd_93, dpd_94, \
                         dpd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_x[k] * dsf0_52[k]
                   + f_0 * dsd_32[k]
                   - f_5 * pc_x[k] * dsf1_52[k];

        t_153[k] = f_1 * dsd_33[k]
                   + f_4 * pc_x[k] * dpd_93[k];

        t_154[k] = f_1 * dsd_34[k]
                   + f_4 * pc_x[k] * dpd_94[k];

        t_155[k] = f_1 * dsd_35[k]
                   + f_4 * pc_x[k] * dpd_95[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pc_x, pc_y, dsf0_56, dsf0_57, \
                         dsf0_59, dsf1_56, dsf1_57, dsf1_59, dpd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_x[k] * dsf0_56[k]
                   - f_5 * pc_x[k] * dsf1_56[k];

        t_157[k] = pb_x[k] * dsf0_57[k]
                   - f_5 * pc_x[k] * dsf1_57[k];

        t_158[k] = f_4 * pc_y[k] * dpd_95[k];

        t_159[k] = pb_x[k] * dsf0_59[k]
                   - f_5 * pc_x[k] * dsf1_59[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pb_y, pc_x, pc_y, dsf0_50, \
                         dsf0_52, dsd_30, dsf1_50, dsf1_52, dpd_96, dpd_99, \
                         dpd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_y[k] * dsf0_50[k]
                   - f_5 * pc_y[k] * dsf1_50[k];

        t_161[k] = f_1 * dsd_30[k]
                   + f_4 * pc_y[k] * dpd_96[k];

        t_162[k] = pb_y[k] * dsf0_52[k]
                   - f_5 * pc_y[k] * dsf1_52[k];

        t_163[k] = f_4 * pc_x[k] * dpd_99[k];

        t_164[k] = f_4 * pc_x[k] * dpd_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_y, pc_x, pc_y, dsf0_56, dsf0_57, \
                         dsd_33, dsd_34, dsd_35, dsf1_56, dsf1_57, \
                         dpd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_4 * pc_x[k] * dpd_101[k];

        t_166[k] = pb_y[k] * dsf0_56[k]
                   + f_6 * dsd_33[k]
                   - f_5 * pc_y[k] * dsf1_56[k];

        t_167[k] = pb_y[k] * dsf0_57[k]
                   + f_0 * dsd_34[k]
                   - f_5 * pc_y[k] * dsf1_57[k];

        t_168[k] = f_1 * dsd_35[k]
                   + f_4 * pc_y[k] * dpd_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pb_y, pc_x, pc_y, dsf0_59, dsf1_59, \
                         dpp0_51, dpp0_53, dpp1_51, dpp1_53, dpd_102, \
                         dpd_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pb_y[k] * dsf0_59[k]
                   - f_5 * pc_y[k] * dsf1_59[k];

        t_170[k] = f_2 * dpp0_51[k]
                   - f_3 * dpp1_51[k]
                   + f_4 * pc_x[k] * dpd_102[k];

        t_171[k] = f_4 * pc_y[k] * dpd_102[k];

        t_172[k] = f_7 * dpp0_53[k]
                   - f_8 * dpp1_53[k]
                   + f_4 * pc_x[k] * dpd_104[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, pc_x, pc_y, dpp0_52, \
                         dpp0_53, dpp1_52, dpp1_53, dpd_105, dpd_106, \
                         dpd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_4 * pc_x[k] * dpd_105[k];

        t_174[k] = f_4 * pc_x[k] * dpd_106[k];

        t_175[k] = f_4 * pc_x[k] * dpd_107[k];

        t_176[k] = f_2 * dpp0_52[k]
                   - f_3 * dpp1_52[k]
                   + f_4 * pc_y[k] * dpd_105[k];

        t_177[k] = f_7 * dpp0_53[k]
                   - f_8 * dpp1_53[k]
                   + f_4 * pc_y[k] * dpd_106[k];

        t_178[k] = f_4 * pc_y[k] * dpd_107[k];
    }

#pragma omp simd aligned(t_179, pc_z, ppd_53, dsd_35, dpp0_53, dpp1_53, \
                         dpd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * ppd_53[k]
                   + f_1 * dsd_35[k]
                   + f_2 * dpp0_53[k]
                   - f_3 * dpp1_53[k]
                   + f_4 * pc_z[k] * dpd_107[k];
    }
}

auto
compute_prim_dpf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppf0,
                                                   const size_t ppd, const size_t ppf1,
                                                   const size_t dsf0, const size_t dsd,
                                                   const size_t dsf1, const size_t dpp0,
                                                   const size_t dpp1, const size_t dpd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dpf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, ppf0,
                                                              ppd, ppf1, dsf0, dsd, dsf1, dpp0,
                                                              dpp1, dpd, ncols, gamma, p, q);

    compute_prim_dpf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppf0,
                                                              ppd, ppf1, dsf0, dsd, dsf1, dpp0,
                                                              dpp1, dpd, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
