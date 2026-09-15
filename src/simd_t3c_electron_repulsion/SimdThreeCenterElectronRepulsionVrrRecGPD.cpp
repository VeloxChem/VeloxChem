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


#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gpd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpd0, const size_t dpd1,
                                                          const size_t fpd0, const size_t fpp,
                                                          const size_t fpd1, const size_t gsd0,
                                                          const size_t gsp, const size_t gsd1,
                                                          const size_t gps0, const size_t gps1,
                                                          const size_t gpp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 1.0 / p;
    const auto f_9 = gamma / (p * q);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpd0_0 = buffer.data(dpd0 + 0);
    const auto *dpd0_27 = buffer.data(dpd0 + 27);
    const auto *dpd0_53 = buffer.data(dpd0 + 53);
    const auto *dpd0_63 = buffer.data(dpd0 + 63);
    const auto *dpd0_107 = buffer.data(dpd0 + 107);

    const auto *dpd1_0 = buffer.data(dpd1 + 0);
    const auto *dpd1_27 = buffer.data(dpd1 + 27);
    const auto *dpd1_53 = buffer.data(dpd1 + 53);
    const auto *dpd1_63 = buffer.data(dpd1 + 63);
    const auto *dpd1_107 = buffer.data(dpd1 + 107);

    const auto *fpd0_0 = buffer.data(fpd0 + 0);
    const auto *fpd0_3 = buffer.data(fpd0 + 3);
    const auto *fpd0_5 = buffer.data(fpd0 + 5);
    const auto *fpd0_6 = buffer.data(fpd0 + 6);
    const auto *fpd0_9 = buffer.data(fpd0 + 9);
    const auto *fpd0_12 = buffer.data(fpd0 + 12);
    const auto *fpd0_17 = buffer.data(fpd0 + 17);
    const auto *fpd0_18 = buffer.data(fpd0 + 18);
    const auto *fpd0_19 = buffer.data(fpd0 + 19);
    const auto *fpd0_21 = buffer.data(fpd0 + 21);
    const auto *fpd0_27 = buffer.data(fpd0 + 27);
    const auto *fpd0_36 = buffer.data(fpd0 + 36);
    const auto *fpd0_38 = buffer.data(fpd0 + 38);
    const auto *fpd0_41 = buffer.data(fpd0 + 41);
    const auto *fpd0_48 = buffer.data(fpd0 + 48);
    const auto *fpd0_53 = buffer.data(fpd0 + 53);
    const auto *fpd0_63 = buffer.data(fpd0 + 63);
    const auto *fpd0_107 = buffer.data(fpd0 + 107);
    const auto *fpd0_114 = buffer.data(fpd0 + 114);
    const auto *fpd0_117 = buffer.data(fpd0 + 117);
    const auto *fpd0_119 = buffer.data(fpd0 + 119);
    const auto *fpd0_123 = buffer.data(fpd0 + 123);

    const auto *fpp_0 = buffer.data(fpp + 0);
    const auto *fpp_1 = buffer.data(fpp + 1);
    const auto *fpp_2 = buffer.data(fpp + 2);
    const auto *fpp_10 = buffer.data(fpp + 10);
    const auto *fpp_12 = buffer.data(fpp + 12);
    const auto *fpp_13 = buffer.data(fpp + 13);
    const auto *fpp_14 = buffer.data(fpp + 14);
    const auto *fpp_16 = buffer.data(fpp + 16);
    const auto *fpp_20 = buffer.data(fpp + 20);
    const auto *fpp_23 = buffer.data(fpp + 23);
    const auto *fpp_24 = buffer.data(fpp + 24);
    const auto *fpp_25 = buffer.data(fpp + 25);
    const auto *fpp_26 = buffer.data(fpp + 26);
    const auto *fpp_28 = buffer.data(fpp + 28);
    const auto *fpp_30 = buffer.data(fpp + 30);
    const auto *fpp_31 = buffer.data(fpp + 31);
    const auto *fpp_34 = buffer.data(fpp + 34);
    const auto *fpp_39 = buffer.data(fpp + 39);
    const auto *fpp_40 = buffer.data(fpp + 40);
    const auto *fpp_41 = buffer.data(fpp + 41);
    const auto *fpp_43 = buffer.data(fpp + 43);
    const auto *fpp_44 = buffer.data(fpp + 44);
    const auto *fpp_47 = buffer.data(fpp + 47);
    const auto *fpp_50 = buffer.data(fpp + 50);
    const auto *fpp_51 = buffer.data(fpp + 51);
    const auto *fpp_53 = buffer.data(fpp + 53);
    const auto *fpp_54 = buffer.data(fpp + 54);
    const auto *fpp_55 = buffer.data(fpp + 55);
    const auto *fpp_57 = buffer.data(fpp + 57);
    const auto *fpp_58 = buffer.data(fpp + 58);
    const auto *fpp_61 = buffer.data(fpp + 61);

    const auto *fpd1_0 = buffer.data(fpd1 + 0);
    const auto *fpd1_3 = buffer.data(fpd1 + 3);
    const auto *fpd1_5 = buffer.data(fpd1 + 5);
    const auto *fpd1_6 = buffer.data(fpd1 + 6);
    const auto *fpd1_9 = buffer.data(fpd1 + 9);
    const auto *fpd1_12 = buffer.data(fpd1 + 12);
    const auto *fpd1_17 = buffer.data(fpd1 + 17);
    const auto *fpd1_18 = buffer.data(fpd1 + 18);
    const auto *fpd1_19 = buffer.data(fpd1 + 19);
    const auto *fpd1_21 = buffer.data(fpd1 + 21);
    const auto *fpd1_27 = buffer.data(fpd1 + 27);
    const auto *fpd1_36 = buffer.data(fpd1 + 36);
    const auto *fpd1_38 = buffer.data(fpd1 + 38);
    const auto *fpd1_41 = buffer.data(fpd1 + 41);
    const auto *fpd1_48 = buffer.data(fpd1 + 48);
    const auto *fpd1_53 = buffer.data(fpd1 + 53);
    const auto *fpd1_63 = buffer.data(fpd1 + 63);
    const auto *fpd1_107 = buffer.data(fpd1 + 107);
    const auto *fpd1_114 = buffer.data(fpd1 + 114);
    const auto *fpd1_117 = buffer.data(fpd1 + 117);
    const auto *fpd1_119 = buffer.data(fpd1 + 119);
    const auto *fpd1_123 = buffer.data(fpd1 + 123);

    const auto *gsd0_0 = buffer.data(gsd0 + 0);
    const auto *gsd0_3 = buffer.data(gsd0 + 3);
    const auto *gsd0_5 = buffer.data(gsd0 + 5);
    const auto *gsd0_9 = buffer.data(gsd0 + 9);
    const auto *gsd0_17 = buffer.data(gsd0 + 17);
    const auto *gsd0_18 = buffer.data(gsd0 + 18);
    const auto *gsd0_21 = buffer.data(gsd0 + 21);
    const auto *gsd0_23 = buffer.data(gsd0 + 23);
    const auto *gsd0_30 = buffer.data(gsd0 + 30);
    const auto *gsd0_33 = buffer.data(gsd0 + 33);
    const auto *gsd0_35 = buffer.data(gsd0 + 35);
    const auto *gsd0_36 = buffer.data(gsd0 + 36);

    const auto *gsp_0 = buffer.data(gsp + 0);
    const auto *gsp_1 = buffer.data(gsp + 1);
    const auto *gsp_2 = buffer.data(gsp + 2);
    const auto *gsp_3 = buffer.data(gsp + 3);
    const auto *gsp_4 = buffer.data(gsp + 4);
    const auto *gsp_6 = buffer.data(gsp + 6);
    const auto *gsp_8 = buffer.data(gsp + 8);
    const auto *gsp_9 = buffer.data(gsp + 9);
    const auto *gsp_10 = buffer.data(gsp + 10);
    const auto *gsp_11 = buffer.data(gsp + 11);
    const auto *gsp_14 = buffer.data(gsp + 14);
    const auto *gsp_15 = buffer.data(gsp + 15);
    const auto *gsp_16 = buffer.data(gsp + 16);
    const auto *gsp_17 = buffer.data(gsp + 17);
    const auto *gsp_18 = buffer.data(gsp + 18);
    const auto *gsp_19 = buffer.data(gsp + 19);

    const auto *gsd1_0 = buffer.data(gsd1 + 0);
    const auto *gsd1_3 = buffer.data(gsd1 + 3);
    const auto *gsd1_5 = buffer.data(gsd1 + 5);
    const auto *gsd1_9 = buffer.data(gsd1 + 9);
    const auto *gsd1_17 = buffer.data(gsd1 + 17);
    const auto *gsd1_18 = buffer.data(gsd1 + 18);
    const auto *gsd1_21 = buffer.data(gsd1 + 21);
    const auto *gsd1_23 = buffer.data(gsd1 + 23);
    const auto *gsd1_30 = buffer.data(gsd1 + 30);
    const auto *gsd1_33 = buffer.data(gsd1 + 33);
    const auto *gsd1_35 = buffer.data(gsd1 + 35);
    const auto *gsd1_36 = buffer.data(gsd1 + 36);

    const auto *gps0_0 = buffer.data(gps0 + 0);
    const auto *gps0_3 = buffer.data(gps0 + 3);
    const auto *gps0_4 = buffer.data(gps0 + 4);
    const auto *gps0_6 = buffer.data(gps0 + 6);
    const auto *gps0_8 = buffer.data(gps0 + 8);
    const auto *gps0_9 = buffer.data(gps0 + 9);
    const auto *gps0_10 = buffer.data(gps0 + 10);
    const auto *gps0_13 = buffer.data(gps0 + 13);
    const auto *gps0_14 = buffer.data(gps0 + 14);
    const auto *gps0_15 = buffer.data(gps0 + 15);
    const auto *gps0_17 = buffer.data(gps0 + 17);
    const auto *gps0_18 = buffer.data(gps0 + 18);

    const auto *gps1_0 = buffer.data(gps1 + 0);
    const auto *gps1_3 = buffer.data(gps1 + 3);
    const auto *gps1_4 = buffer.data(gps1 + 4);
    const auto *gps1_6 = buffer.data(gps1 + 6);
    const auto *gps1_8 = buffer.data(gps1 + 8);
    const auto *gps1_9 = buffer.data(gps1 + 9);
    const auto *gps1_10 = buffer.data(gps1 + 10);
    const auto *gps1_13 = buffer.data(gps1 + 13);
    const auto *gps1_14 = buffer.data(gps1 + 14);
    const auto *gps1_15 = buffer.data(gps1 + 15);
    const auto *gps1_17 = buffer.data(gps1 + 17);
    const auto *gps1_18 = buffer.data(gps1 + 18);

    const auto *gpp_0 = buffer.data(gpp + 0);
    const auto *gpp_1 = buffer.data(gpp + 1);
    const auto *gpp_2 = buffer.data(gpp + 2);
    const auto *gpp_3 = buffer.data(gpp + 3);
    const auto *gpp_5 = buffer.data(gpp + 5);
    const auto *gpp_6 = buffer.data(gpp + 6);
    const auto *gpp_8 = buffer.data(gpp + 8);
    const auto *gpp_9 = buffer.data(gpp + 9);
    const auto *gpp_10 = buffer.data(gpp + 10);
    const auto *gpp_12 = buffer.data(gpp + 12);
    const auto *gpp_13 = buffer.data(gpp + 13);
    const auto *gpp_14 = buffer.data(gpp + 14);
    const auto *gpp_15 = buffer.data(gpp + 15);
    const auto *gpp_16 = buffer.data(gpp + 16);
    const auto *gpp_18 = buffer.data(gpp + 18);
    const auto *gpp_20 = buffer.data(gpp + 20);
    const auto *gpp_21 = buffer.data(gpp + 21);
    const auto *gpp_23 = buffer.data(gpp + 23);
    const auto *gpp_24 = buffer.data(gpp + 24);
    const auto *gpp_25 = buffer.data(gpp + 25);
    const auto *gpp_26 = buffer.data(gpp + 26);
    const auto *gpp_27 = buffer.data(gpp + 27);
    const auto *gpp_28 = buffer.data(gpp + 28);
    const auto *gpp_29 = buffer.data(gpp + 29);
    const auto *gpp_30 = buffer.data(gpp + 30);
    const auto *gpp_31 = buffer.data(gpp + 31);
    const auto *gpp_32 = buffer.data(gpp + 32);
    const auto *gpp_33 = buffer.data(gpp + 33);
    const auto *gpp_34 = buffer.data(gpp + 34);
    const auto *gpp_38 = buffer.data(gpp + 38);
    const auto *gpp_39 = buffer.data(gpp + 39);
    const auto *gpp_40 = buffer.data(gpp + 40);
    const auto *gpp_41 = buffer.data(gpp + 41);
    const auto *gpp_43 = buffer.data(gpp + 43);
    const auto *gpp_44 = buffer.data(gpp + 44);
    const auto *gpp_45 = buffer.data(gpp + 45);
    const auto *gpp_46 = buffer.data(gpp + 46);
    const auto *gpp_47 = buffer.data(gpp + 47);
    const auto *gpp_48 = buffer.data(gpp + 48);
    const auto *gpp_50 = buffer.data(gpp + 50);
    const auto *gpp_51 = buffer.data(gpp + 51);
    const auto *gpp_52 = buffer.data(gpp + 52);
    const auto *gpp_53 = buffer.data(gpp + 53);
    const auto *gpp_54 = buffer.data(gpp + 54);
    const auto *gpp_55 = buffer.data(gpp + 55);
    const auto *gpp_56 = buffer.data(gpp + 56);
    const auto *gpp_57 = buffer.data(gpp + 57);
    const auto *gpp_58 = buffer.data(gpp + 58);
    const auto *gpp_60 = buffer.data(gpp + 60);
    const auto *gpp_61 = buffer.data(gpp + 61);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fpp_0, gsp_0, gps0_0, \
                         gps1_0, gpp_0, gpp_1, gpp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fpp_0[k]
                 + f_1 * gsp_0[k]
                 + f_2 * gps0_0[k]
                 - f_3 * gps1_0[k]
                 + f_4 * pc_x[k] * gpp_0[k];

        t_1[k] = f_4 * pc_y[k] * gpp_0[k];

        t_2[k] = f_4 * pc_z[k] * gpp_0[k];

        t_3[k] = f_2 * gps0_0[k]
                 - f_3 * gps1_0[k]
                 + f_4 * pc_y[k] * gpp_1[k];

        t_4[k] = f_4 * pc_y[k] * gpp_2[k];

        t_5[k] = f_2 * gps0_0[k]
                 - f_3 * gps1_0[k]
                 + f_4 * pc_z[k] * gpp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pc_y, pc_z, gsd0_0, gsd0_3, gsp_0, gsp_1, \
                         gsd1_0, gsd1_3, gpp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * gsd0_0[k]
                 - f_5 * pc_y[k] * gsd1_0[k];

        t_7[k] = f_1 * gsp_0[k]
                 + f_4 * pc_y[k] * gpp_3[k];

        t_8[k] = f_4 * pc_z[k] * gpp_3[k];

        t_9[k] = pb_y[k] * gsd0_3[k]
                 + f_6 * gsp_1[k]
                 - f_5 * pc_y[k] * gsd1_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, pc_y, pc_z, gsd0_0, gsd0_5, \
                         gsp_2, gsd1_0, gsd1_5, gpp_5, gpp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * gsp_2[k]
                  + f_4 * pc_y[k] * gpp_5[k];

        t_11[k] = pb_y[k] * gsd0_5[k]
                  - f_5 * pc_y[k] * gsd1_5[k];

        t_12[k] = pb_z[k] * gsd0_0[k]
                  - f_5 * pc_z[k] * gsd1_0[k];

        t_13[k] = f_4 * pc_y[k] * gpp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_z, pc_y, pc_z, gsd0_3, gsd0_5, gsp_0, \
                         gsp_2, gsd1_3, gsd1_5, gpp_6, gpp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * gsp_0[k]
                  + f_4 * pc_z[k] * gpp_6[k];

        t_15[k] = pb_z[k] * gsd0_3[k]
                  - f_5 * pc_z[k] * gsd1_3[k];

        t_16[k] = f_4 * pc_y[k] * gpp_8[k];

        t_17[k] = pb_z[k] * gsd0_5[k]
                  + f_6 * gsp_2[k]
                  - f_5 * pc_z[k] * gsd1_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pc_x, pc_y, pc_z, fpd0_0, fpp_10, fpd1_0, \
                         gsp_4, gpp_9, gpp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_y[k] * fpd0_0[k]
                  - f_5 * pc_y[k] * fpd1_0[k];

        t_19[k] = f_7 * fpp_10[k]
                  + f_1 * gsp_4[k]
                  + f_4 * pc_x[k] * gpp_10[k];

        t_20[k] = f_4 * pc_z[k] * gpp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pc_y, pc_z, fpd0_5, fpp_1, fpd1_5, gps0_3, \
                         gps1_3, gpp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * fpp_1[k]
                  + f_2 * gps0_3[k]
                  - f_3 * gps1_3[k]
                  + f_4 * pc_y[k] * gpp_10[k];

        t_22[k] = f_4 * pc_z[k] * gpp_10[k];

        t_23[k] = pa_y[k] * fpd0_5[k]
                  - f_5 * pc_y[k] * fpd1_5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pc_x, pc_z, fpp_12, fpp_13, gps0_4, gps1_4, gpp_12, \
                         gpp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * fpp_12[k]
                  + f_2 * gps0_4[k]
                  - f_3 * gps1_4[k]
                  + f_4 * pc_x[k] * gpp_12[k];

        t_25[k] = f_7 * fpp_13[k]
                  + f_4 * pc_x[k] * gpp_13[k];

        t_26[k] = f_4 * pc_z[k] * gpp_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pc_x, pc_z, dpd0_27, dpd1_27, fpd0_27, \
                         fpd1_27, gps0_4, gps1_4, gpp_13, gpp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_8 * dpd0_27[k]
                  - f_9 * dpd1_27[k]
                  + pa_x[k] * fpd0_27[k]
                  - f_5 * pc_x[k] * fpd1_27[k];

        t_28[k] = f_4 * pc_z[k] * gpp_13[k];

        t_29[k] = f_2 * gps0_4[k]
                  - f_3 * gps1_4[k]
                  + f_4 * pc_z[k] * gpp_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pc_x, pc_y, pc_z, fpd0_12, fpp_16, fpd1_12, \
                         gsp_3, gpp_15, gpp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * fpd0_12[k]
                  - f_5 * pc_y[k] * fpd1_12[k];

        t_31[k] = f_7 * fpp_16[k]
                  + f_4 * pc_x[k] * gpp_16[k];

        t_32[k] = f_1 * gsp_3[k]
                  + f_4 * pc_z[k] * gpp_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_z, pc_y, pc_z, fpd0_17, fpd1_17, gsd0_9, \
                         gsp_4, gsd1_9, gpp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * gsd0_9[k]
                  - f_5 * pc_z[k] * gsd1_9[k];

        t_34[k] = f_1 * gsp_4[k]
                  + f_4 * pc_z[k] * gpp_16[k];

        t_35[k] = pa_y[k] * fpd0_17[k]
                  - f_5 * pc_y[k] * fpd1_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pc_x, pc_y, pc_z, fpd0_0, fpd0_3, \
                         fpp_20, fpd1_0, fpd1_3, gsp_8, gpp_18, \
                         gpp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_z[k] * fpd0_0[k]
                  - f_5 * pc_z[k] * fpd1_0[k];

        t_37[k] = f_4 * pc_y[k] * gpp_18[k];

        t_38[k] = f_7 * fpp_20[k]
                  + f_1 * gsp_8[k]
                  + f_4 * pc_x[k] * gpp_20[k];

        t_39[k] = pa_z[k] * fpd0_3[k]
                  - f_5 * pc_z[k] * fpd1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pc_y, pc_z, fpd0_6, fpp_2, fpd1_6, \
                         gsp_6, gps0_6, gps1_6, gpp_20, gpp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_4 * pc_y[k] * gpp_20[k];

        t_41[k] = f_1 * fpp_2[k]
                  + f_2 * gps0_6[k]
                  - f_3 * gps1_6[k]
                  + f_4 * pc_z[k] * gpp_20[k];

        t_42[k] = pa_z[k] * fpd0_6[k]
                  - f_5 * pc_z[k] * fpd1_6[k];

        t_43[k] = f_1 * gsp_6[k]
                  + f_4 * pc_y[k] * gpp_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_y, pc_x, pc_y, pc_z, fpd0_9, fpp_23, \
                         fpd1_9, gsd0_17, gsp_8, gsd1_17, gpp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_7 * fpp_23[k]
                  + f_4 * pc_x[k] * gpp_23[k];

        t_45[k] = pa_z[k] * fpd0_9[k]
                  - f_5 * pc_z[k] * fpd1_9[k];

        t_46[k] = f_1 * gsp_8[k]
                  + f_4 * pc_y[k] * gpp_23[k];

        t_47[k] = pb_y[k] * gsd0_17[k]
                  - f_5 * pc_y[k] * gsd1_17[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pc_x, pc_y, fpp_24, fpp_26, gps0_8, \
                         gps1_8, gpp_24, gpp_25, gpp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_7 * fpp_24[k]
                  + f_2 * gps0_8[k]
                  - f_3 * gps1_8[k]
                  + f_4 * pc_x[k] * gpp_24[k];

        t_49[k] = f_4 * pc_y[k] * gpp_24[k];

        t_50[k] = f_7 * fpp_26[k]
                  + f_4 * pc_x[k] * gpp_26[k];

        t_51[k] = f_2 * gps0_8[k]
                  - f_3 * gps1_8[k]
                  + f_4 * pc_y[k] * gpp_25[k];

        t_52[k] = f_4 * pc_y[k] * gpp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, pa_x, pa_y, pc_x, pc_y, dpd0_0, dpd0_53, dpd1_0, dpd1_53, \
                         fpd0_18, fpd0_53, fpd1_18, fpd1_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * dpd0_53[k]
                  - f_9 * dpd1_53[k]
                  + pa_x[k] * fpd0_53[k]
                  - f_5 * pc_x[k] * fpd1_53[k];

        t_54[k] = f_10 * dpd0_0[k]
                  - f_11 * dpd1_0[k]
                  + pa_y[k] * fpd0_18[k]
                  - f_5 * pc_y[k] * fpd1_18[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, fpp_10, fpp_28, \
                         gsp_10, gps0_9, gps1_9, gpp_27, gpp_28, \
                         gpp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_6 * fpp_28[k]
                  + f_1 * gsp_10[k]
                  + f_4 * pc_x[k] * gpp_28[k];

        t_56[k] = f_4 * pc_z[k] * gpp_27[k];

        t_57[k] = f_6 * fpp_10[k]
                  + f_2 * gps0_9[k]
                  - f_3 * gps1_9[k]
                  + f_4 * pc_y[k] * gpp_28[k];

        t_58[k] = f_4 * pc_z[k] * gpp_28[k];

        t_59[k] = f_2 * gps0_9[k]
                  - f_3 * gps1_9[k]
                  + f_4 * pc_z[k] * gpp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pc_x, pc_z, fpp_30, fpp_31, gps0_10, gps1_10, \
                         gpp_30, gpp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * fpp_30[k]
                  + f_2 * gps0_10[k]
                  - f_3 * gps1_10[k]
                  + f_4 * pc_x[k] * gpp_30[k];

        t_61[k] = f_6 * fpp_31[k]
                  + f_4 * pc_x[k] * gpp_31[k];

        t_62[k] = f_4 * pc_z[k] * gpp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_x, pc_x, pc_z, dpd0_63, dpd1_63, fpd0_63, \
                         fpd1_63, gps0_10, gps1_10, gpp_31, gpp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_10 * dpd0_63[k]
                  - f_11 * dpd1_63[k]
                  + pa_x[k] * fpd0_63[k]
                  - f_5 * pc_x[k] * fpd1_63[k];

        t_64[k] = f_4 * pc_z[k] * gpp_31[k];

        t_65[k] = f_2 * gps0_10[k]
                  - f_3 * gps1_10[k]
                  + f_4 * pc_z[k] * gpp_32[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_z, pc_x, pc_z, fpp_34, gsd0_18, gsd0_21, \
                         gsp_9, gsd1_18, gsd1_21, gpp_33, gpp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * gsd0_18[k]
                  - f_5 * pc_z[k] * gsd1_18[k];

        t_67[k] = f_6 * fpp_34[k]
                  + f_4 * pc_x[k] * gpp_34[k];

        t_68[k] = f_1 * gsp_9[k]
                  + f_4 * pc_z[k] * gpp_33[k];

        t_69[k] = pb_z[k] * gsd0_21[k]
                  - f_5 * pc_z[k] * gsd1_21[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pb_z, pc_y, pc_z, fpd0_36, fpd1_36, gsd0_23, \
                         gsp_10, gsp_11, gsd1_23, gpp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * gsp_10[k]
                  + f_4 * pc_z[k] * gpp_34[k];

        t_71[k] = pb_z[k] * gsd0_23[k]
                  + f_6 * gsp_11[k]
                  - f_5 * pc_z[k] * gsd1_23[k];

        t_72[k] = pa_y[k] * fpd0_36[k]
                  - f_5 * pc_y[k] * fpd1_36[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pa_z, pc_y, pc_z, fpd0_19, fpd0_21, \
                         fpd0_38, fpp_20, fpd1_19, fpd1_21, fpd1_38, \
                         gpp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * fpd0_19[k]
                  - f_5 * pc_z[k] * fpd1_19[k];

        t_74[k] = pa_y[k] * fpd0_38[k]
                  - f_5 * pc_y[k] * fpd1_38[k];

        t_75[k] = pa_z[k] * fpd0_21[k]
                  - f_5 * pc_z[k] * fpd1_21[k];

        t_76[k] = f_1 * fpp_20[k]
                  + f_4 * pc_y[k] * gpp_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_y, pc_x, pc_y, fpd0_41, fpp_39, fpp_40, fpd1_41, \
                         gps0_13, gps1_13, gpp_39, gpp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_y[k] * fpd0_41[k]
                  - f_5 * pc_y[k] * fpd1_41[k];

        t_78[k] = f_6 * fpp_39[k]
                  + f_2 * gps0_13[k]
                  - f_3 * gps1_13[k]
                  + f_4 * pc_x[k] * gpp_39[k];

        t_79[k] = f_6 * fpp_40[k]
                  + f_4 * pc_x[k] * gpp_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_z, pc_x, pc_y, pc_z, fpd0_27, fpp_23, fpp_41, \
                         fpd1_27, gsp_14, gpp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * fpp_41[k]
                  + f_4 * pc_x[k] * gpp_41[k];

        t_81[k] = pa_z[k] * fpd0_27[k]
                  - f_5 * pc_z[k] * fpd1_27[k];

        t_82[k] = f_1 * fpp_23[k]
                  + f_1 * gsp_14[k]
                  + f_4 * pc_y[k] * gpp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_y, pc_x, pc_y, pc_z, fpd0_48, fpp_14, fpp_43, \
                         fpd1_48, gps0_13, gps1_13, gpp_41, gpp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * fpp_14[k]
                  + f_2 * gps0_13[k]
                  - f_3 * gps1_13[k]
                  + f_4 * pc_z[k] * gpp_41[k];

        t_84[k] = pa_y[k] * fpd0_48[k]
                  - f_5 * pc_y[k] * fpd1_48[k];

        t_85[k] = f_6 * fpp_43[k]
                  + f_4 * pc_x[k] * gpp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pc_x, pc_y, fpd0_53, fpp_25, fpp_26, \
                         fpp_44, fpd1_53, gps0_14, gps1_14, gpp_43, \
                         gpp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_6 * fpp_44[k]
                  + f_4 * pc_x[k] * gpp_44[k];

        t_87[k] = f_1 * fpp_25[k]
                  + f_2 * gps0_14[k]
                  - f_3 * gps1_14[k]
                  + f_4 * pc_y[k] * gpp_43[k];

        t_88[k] = f_1 * fpp_26[k]
                  + f_4 * pc_y[k] * gpp_44[k];

        t_89[k] = pa_y[k] * fpd0_53[k]
                  - f_5 * pc_y[k] * fpd1_53[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_z, pc_x, pc_y, pc_z, dpd0_0, dpd1_0, fpd0_36, \
                         fpp_47, fpd1_36, gsp_17, gpp_45, gpp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_10 * dpd0_0[k]
                  - f_11 * dpd1_0[k]
                  + pa_z[k] * fpd0_36[k]
                  - f_5 * pc_z[k] * fpd1_36[k];

        t_91[k] = f_4 * pc_y[k] * gpp_45[k];

        t_92[k] = f_6 * fpp_47[k]
                  + f_1 * gsp_17[k]
                  + f_4 * pc_x[k] * gpp_47[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_y, pc_z, fpp_20, gsd0_30, gsd1_30, \
                         gps0_15, gps1_15, gpp_46, gpp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * gps0_15[k]
                  - f_3 * gps1_15[k]
                  + f_4 * pc_y[k] * gpp_46[k];

        t_94[k] = f_4 * pc_y[k] * gpp_47[k];

        t_95[k] = f_6 * fpp_20[k]
                  + f_2 * gps0_15[k]
                  - f_3 * gps1_15[k]
                  + f_4 * pc_z[k] * gpp_47[k];

        t_96[k] = pb_y[k] * gsd0_30[k]
                  - f_5 * pc_y[k] * gsd1_30[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_y, pc_x, pc_y, fpp_50, gsd0_33, gsp_15, \
                         gsp_16, gsp_17, gsd1_33, gpp_48, gpp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * gsp_15[k]
                  + f_4 * pc_y[k] * gpp_48[k];

        t_98[k] = f_6 * fpp_50[k]
                  + f_4 * pc_x[k] * gpp_50[k];

        t_99[k] = pb_y[k] * gsd0_33[k]
                  + f_6 * gsp_16[k]
                  - f_5 * pc_y[k] * gsd1_33[k];

        t_100[k] = f_1 * gsp_17[k]
                   + f_4 * pc_y[k] * gpp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_y, pc_x, pc_y, fpp_51, fpp_53, \
                         gsd0_35, gsd1_35, gps0_17, gps1_17, gpp_51, \
                         gpp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_y[k] * gsd0_35[k]
                   - f_5 * pc_y[k] * gsd1_35[k];

        t_102[k] = f_6 * fpp_51[k]
                   + f_2 * gps0_17[k]
                   - f_3 * gps1_17[k]
                   + f_4 * pc_x[k] * gpp_51[k];

        t_103[k] = f_4 * pc_y[k] * gpp_51[k];

        t_104[k] = f_6 * fpp_53[k]
                   + f_4 * pc_x[k] * gpp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pc_x, pc_y, dpd0_107, dpd1_107, fpd0_107, \
                         fpd1_107, gps0_17, gps1_17, gpp_52, gpp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_2 * gps0_17[k]
                   - f_3 * gps1_17[k]
                   + f_4 * pc_y[k] * gpp_52[k];

        t_106[k] = f_4 * pc_y[k] * gpp_53[k];

        t_107[k] = f_10 * dpd0_107[k]
                   - f_11 * dpd1_107[k]
                   + pa_x[k] * fpd0_107[k]
                   - f_5 * pc_x[k] * fpd1_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_x, pc_y, pc_z, fpp_28, fpp_54, fpp_55, \
                         gsp_18, gsp_19, gps0_18, gps1_18, gpp_54, \
                         gpp_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * fpp_54[k]
                   + f_1 * gsp_18[k]
                   + f_2 * gps0_18[k]
                   - f_3 * gps1_18[k]
                   + f_4 * pc_x[k] * gpp_54[k];

        t_109[k] = f_1 * fpp_55[k]
                   + f_1 * gsp_19[k]
                   + f_4 * pc_x[k] * gpp_55[k];

        t_110[k] = f_4 * pc_z[k] * gpp_54[k];

        t_111[k] = f_7 * fpp_28[k]
                   + f_2 * gps0_18[k]
                   - f_3 * gps1_18[k]
                   + f_4 * pc_y[k] * gpp_55[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_x, pc_x, pc_z, fpd0_114, fpp_57, \
                         fpp_58, fpd1_114, gps0_18, gps1_18, gpp_55, gpp_56, \
                         gpp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_4 * pc_z[k] * gpp_55[k];

        t_113[k] = f_2 * gps0_18[k]
                   - f_3 * gps1_18[k]
                   + f_4 * pc_z[k] * gpp_56[k];

        t_114[k] = pa_x[k] * fpd0_114[k]
                   + f_6 * fpp_57[k]
                   - f_5 * pc_x[k] * fpd1_114[k];

        t_115[k] = f_1 * fpp_58[k]
                   + f_4 * pc_x[k] * gpp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_x, pc_x, pc_z, fpd0_117, fpd0_119, \
                         fpd1_117, fpd1_119, gpp_57, gpp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_4 * pc_z[k] * gpp_57[k];

        t_117[k] = pa_x[k] * fpd0_117[k]
                   - f_5 * pc_x[k] * fpd1_117[k];

        t_118[k] = f_4 * pc_z[k] * gpp_58[k];

        t_119[k] = pa_x[k] * fpd0_119[k]
                   - f_5 * pc_x[k] * fpd1_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, pb_z, pc_x, pc_z, fpd0_123, fpp_61, \
                         fpd1_123, gsd0_36, gsp_18, gsd1_36, gpp_60, \
                         gpp_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_z[k] * gsd0_36[k]
                   - f_5 * pc_z[k] * gsd1_36[k];

        t_121[k] = f_1 * fpp_61[k]
                   + f_4 * pc_x[k] * gpp_61[k];

        t_122[k] = f_1 * gsp_18[k]
                   + f_4 * pc_z[k] * gpp_60[k];

        t_123[k] = pa_x[k] * fpd0_123[k]
                   - f_5 * pc_x[k] * fpd1_123[k];
    }
}

static auto
compute_prim_gpd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpd0, const size_t dpd1,
                                                          const size_t fpd0, const size_t fpp,
                                                          const size_t fpd1, const size_t gsd0,
                                                          const size_t gsp, const size_t gsd1,
                                                          const size_t gps0, const size_t gps1,
                                                          const size_t gpp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 1.0 / p;
    const auto f_9 = gamma / (p * q);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpd0_63 = buffer.data(dpd0 + 63);
    const auto *dpd0_89 = buffer.data(dpd0 + 89);
    const auto *dpd0_107 = buffer.data(dpd0 + 107);

    const auto *dpd1_63 = buffer.data(dpd1 + 63);
    const auto *dpd1_89 = buffer.data(dpd1 + 89);
    const auto *dpd1_107 = buffer.data(dpd1 + 107);

    const auto *fpd0_54 = buffer.data(fpd0 + 54);
    const auto *fpd0_55 = buffer.data(fpd0 + 55);
    const auto *fpd0_60 = buffer.data(fpd0 + 60);
    const auto *fpd0_90 = buffer.data(fpd0 + 90);
    const auto *fpd0_92 = buffer.data(fpd0 + 92);
    const auto *fpd0_102 = buffer.data(fpd0 + 102);
    const auto *fpd0_108 = buffer.data(fpd0 + 108);
    const auto *fpd0_111 = buffer.data(fpd0 + 111);
    const auto *fpd0_114 = buffer.data(fpd0 + 114);
    const auto *fpd0_117 = buffer.data(fpd0 + 117);
    const auto *fpd0_125 = buffer.data(fpd0 + 125);
    const auto *fpd0_135 = buffer.data(fpd0 + 135);
    const auto *fpd0_136 = buffer.data(fpd0 + 136);
    const auto *fpd0_137 = buffer.data(fpd0 + 137);
    const auto *fpd0_141 = buffer.data(fpd0 + 141);
    const auto *fpd0_143 = buffer.data(fpd0 + 143);
    const auto *fpd0_153 = buffer.data(fpd0 + 153);
    const auto *fpd0_154 = buffer.data(fpd0 + 154);
    const auto *fpd0_155 = buffer.data(fpd0 + 155);
    const auto *fpd0_159 = buffer.data(fpd0 + 159);
    const auto *fpd0_161 = buffer.data(fpd0 + 161);
    const auto *fpd0_162 = buffer.data(fpd0 + 162);
    const auto *fpd0_167 = buffer.data(fpd0 + 167);
    const auto *fpd0_171 = buffer.data(fpd0 + 171);
    const auto *fpd0_173 = buffer.data(fpd0 + 173);
    const auto *fpd0_174 = buffer.data(fpd0 + 174);
    const auto *fpd0_177 = buffer.data(fpd0 + 177);
    const auto *fpd0_179 = buffer.data(fpd0 + 179);

    const auto *fpp_29 = buffer.data(fpp + 29);
    const auto *fpp_37 = buffer.data(fpp + 37);
    const auto *fpp_38 = buffer.data(fpp + 38);
    const auto *fpp_44 = buffer.data(fpp + 44);
    const auto *fpp_46 = buffer.data(fpp + 46);
    const auto *fpp_47 = buffer.data(fpp + 47);
    const auto *fpp_53 = buffer.data(fpp + 53);
    const auto *fpp_58 = buffer.data(fpp + 58);
    const auto *fpp_59 = buffer.data(fpp + 59);
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
    const auto *fpp_81 = buffer.data(fpp + 81);
    const auto *fpp_83 = buffer.data(fpp + 83);
    const auto *fpp_85 = buffer.data(fpp + 85);
    const auto *fpp_86 = buffer.data(fpp + 86);
    const auto *fpp_87 = buffer.data(fpp + 87);
    const auto *fpp_89 = buffer.data(fpp + 89);

    const auto *fpd1_54 = buffer.data(fpd1 + 54);
    const auto *fpd1_55 = buffer.data(fpd1 + 55);
    const auto *fpd1_60 = buffer.data(fpd1 + 60);
    const auto *fpd1_90 = buffer.data(fpd1 + 90);
    const auto *fpd1_92 = buffer.data(fpd1 + 92);
    const auto *fpd1_102 = buffer.data(fpd1 + 102);
    const auto *fpd1_108 = buffer.data(fpd1 + 108);
    const auto *fpd1_111 = buffer.data(fpd1 + 111);
    const auto *fpd1_114 = buffer.data(fpd1 + 114);
    const auto *fpd1_117 = buffer.data(fpd1 + 117);
    const auto *fpd1_125 = buffer.data(fpd1 + 125);
    const auto *fpd1_135 = buffer.data(fpd1 + 135);
    const auto *fpd1_136 = buffer.data(fpd1 + 136);
    const auto *fpd1_137 = buffer.data(fpd1 + 137);
    const auto *fpd1_141 = buffer.data(fpd1 + 141);
    const auto *fpd1_143 = buffer.data(fpd1 + 143);
    const auto *fpd1_153 = buffer.data(fpd1 + 153);
    const auto *fpd1_154 = buffer.data(fpd1 + 154);
    const auto *fpd1_155 = buffer.data(fpd1 + 155);
    const auto *fpd1_159 = buffer.data(fpd1 + 159);
    const auto *fpd1_161 = buffer.data(fpd1 + 161);
    const auto *fpd1_162 = buffer.data(fpd1 + 162);
    const auto *fpd1_167 = buffer.data(fpd1 + 167);
    const auto *fpd1_171 = buffer.data(fpd1 + 171);
    const auto *fpd1_173 = buffer.data(fpd1 + 173);
    const auto *fpd1_174 = buffer.data(fpd1 + 174);
    const auto *fpd1_177 = buffer.data(fpd1 + 177);
    const auto *fpd1_179 = buffer.data(fpd1 + 179);

    const auto *gsd0_54 = buffer.data(gsd0 + 54);
    const auto *gsd0_60 = buffer.data(gsd0 + 60);
    const auto *gsd0_63 = buffer.data(gsd0 + 63);
    const auto *gsd0_65 = buffer.data(gsd0 + 65);
    const auto *gsd0_71 = buffer.data(gsd0 + 71);
    const auto *gsd0_72 = buffer.data(gsd0 + 72);
    const auto *gsd0_75 = buffer.data(gsd0 + 75);
    const auto *gsd0_77 = buffer.data(gsd0 + 77);
    const auto *gsd0_81 = buffer.data(gsd0 + 81);

    const auto *gsp_19 = buffer.data(gsp + 19);
    const auto *gsp_23 = buffer.data(gsp + 23);
    const auto *gsp_25 = buffer.data(gsp + 25);
    const auto *gsp_27 = buffer.data(gsp + 27);
    const auto *gsp_29 = buffer.data(gsp + 29);
    const auto *gsp_30 = buffer.data(gsp + 30);
    const auto *gsp_31 = buffer.data(gsp + 31);
    const auto *gsp_32 = buffer.data(gsp + 32);
    const auto *gsp_34 = buffer.data(gsp + 34);
    const auto *gsp_35 = buffer.data(gsp + 35);
    const auto *gsp_36 = buffer.data(gsp + 36);
    const auto *gsp_37 = buffer.data(gsp + 37);
    const auto *gsp_38 = buffer.data(gsp + 38);
    const auto *gsp_40 = buffer.data(gsp + 40);
    const auto *gsp_41 = buffer.data(gsp + 41);

    const auto *gsd1_54 = buffer.data(gsd1 + 54);
    const auto *gsd1_60 = buffer.data(gsd1 + 60);
    const auto *gsd1_63 = buffer.data(gsd1 + 63);
    const auto *gsd1_65 = buffer.data(gsd1 + 65);
    const auto *gsd1_71 = buffer.data(gsd1 + 71);
    const auto *gsd1_72 = buffer.data(gsd1 + 72);
    const auto *gsd1_75 = buffer.data(gsd1 + 75);
    const auto *gsd1_77 = buffer.data(gsd1 + 77);
    const auto *gsd1_81 = buffer.data(gsd1 + 81);

    const auto *gps0_21 = buffer.data(gps0 + 21);
    const auto *gps0_23 = buffer.data(gps0 + 23);
    const auto *gps0_24 = buffer.data(gps0 + 24);
    const auto *gps0_25 = buffer.data(gps0 + 25);
    const auto *gps0_27 = buffer.data(gps0 + 27);
    const auto *gps0_31 = buffer.data(gps0 + 31);
    const auto *gps0_34 = buffer.data(gps0 + 34);
    const auto *gps0_35 = buffer.data(gps0 + 35);
    const auto *gps0_37 = buffer.data(gps0 + 37);
    const auto *gps0_38 = buffer.data(gps0 + 38);
    const auto *gps0_40 = buffer.data(gps0 + 40);

    const auto *gps1_21 = buffer.data(gps1 + 21);
    const auto *gps1_23 = buffer.data(gps1 + 23);
    const auto *gps1_24 = buffer.data(gps1 + 24);
    const auto *gps1_25 = buffer.data(gps1 + 25);
    const auto *gps1_27 = buffer.data(gps1 + 27);
    const auto *gps1_31 = buffer.data(gps1 + 31);
    const auto *gps1_34 = buffer.data(gps1 + 34);
    const auto *gps1_35 = buffer.data(gps1 + 35);
    const auto *gps1_37 = buffer.data(gps1 + 37);
    const auto *gps1_38 = buffer.data(gps1 + 38);
    const auto *gps1_40 = buffer.data(gps1 + 40);

    const auto *gpp_61 = buffer.data(gpp + 61);
    const auto *gpp_64 = buffer.data(gpp + 64);
    const auto *gpp_65 = buffer.data(gpp + 65);
    const auto *gpp_67 = buffer.data(gpp + 67);
    const auto *gpp_68 = buffer.data(gpp + 68);
    const auto *gpp_69 = buffer.data(gpp + 69);
    const auto *gpp_70 = buffer.data(gpp + 70);
    const auto *gpp_71 = buffer.data(gpp + 71);
    const auto *gpp_73 = buffer.data(gpp + 73);
    const auto *gpp_74 = buffer.data(gpp + 74);
    const auto *gpp_75 = buffer.data(gpp + 75);
    const auto *gpp_76 = buffer.data(gpp + 76);
    const auto *gpp_77 = buffer.data(gpp + 77);
    const auto *gpp_79 = buffer.data(gpp + 79);
    const auto *gpp_80 = buffer.data(gpp + 80);
    const auto *gpp_81 = buffer.data(gpp + 81);
    const auto *gpp_82 = buffer.data(gpp + 82);
    const auto *gpp_83 = buffer.data(gpp + 83);
    const auto *gpp_84 = buffer.data(gpp + 84);
    const auto *gpp_86 = buffer.data(gpp + 86);
    const auto *gpp_87 = buffer.data(gpp + 87);
    const auto *gpp_89 = buffer.data(gpp + 89);
    const auto *gpp_91 = buffer.data(gpp + 91);
    const auto *gpp_92 = buffer.data(gpp + 92);
    const auto *gpp_93 = buffer.data(gpp + 93);
    const auto *gpp_94 = buffer.data(gpp + 94);
    const auto *gpp_95 = buffer.data(gpp + 95);
    const auto *gpp_97 = buffer.data(gpp + 97);
    const auto *gpp_98 = buffer.data(gpp + 98);
    const auto *gpp_100 = buffer.data(gpp + 100);
    const auto *gpp_101 = buffer.data(gpp + 101);
    const auto *gpp_103 = buffer.data(gpp + 103);
    const auto *gpp_104 = buffer.data(gpp + 104);
    const auto *gpp_105 = buffer.data(gpp + 105);
    const auto *gpp_106 = buffer.data(gpp + 106);
    const auto *gpp_107 = buffer.data(gpp + 107);
    const auto *gpp_109 = buffer.data(gpp + 109);
    const auto *gpp_110 = buffer.data(gpp + 110);
    const auto *gpp_111 = buffer.data(gpp + 111);
    const auto *gpp_112 = buffer.data(gpp + 112);
    const auto *gpp_113 = buffer.data(gpp + 113);
    const auto *gpp_114 = buffer.data(gpp + 114);
    const auto *gpp_115 = buffer.data(gpp + 115);
    const auto *gpp_116 = buffer.data(gpp + 116);
    const auto *gpp_118 = buffer.data(gpp + 118);
    const auto *gpp_119 = buffer.data(gpp + 119);
    const auto *gpp_120 = buffer.data(gpp + 120);
    const auto *gpp_121 = buffer.data(gpp + 121);
    const auto *gpp_122 = buffer.data(gpp + 122);
    const auto *gpp_124 = buffer.data(gpp + 124);
    const auto *gpp_125 = buffer.data(gpp + 125);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pa_z, pc_x, pc_z, fpd0_54, fpd0_55, \
                         fpd0_125, fpd1_54, fpd1_55, fpd1_125, gsp_19, \
                         gpp_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * gsp_19[k]
                   + f_4 * pc_z[k] * gpp_61[k];

        t_125[k] = pa_x[k] * fpd0_125[k]
                   - f_5 * pc_x[k] * fpd1_125[k];

        t_126[k] = pa_z[k] * fpd0_54[k]
                   - f_5 * pc_z[k] * fpd1_54[k];

        t_127[k] = pa_z[k] * fpd0_55[k]
                   - f_5 * pc_z[k] * fpd1_55[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, fpp_29, fpp_37, fpp_38, \
                         fpp_65, gsp_23, gps0_21, gps1_21, gpp_64, \
                         gpp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * fpp_65[k]
                   + f_1 * gsp_23[k]
                   + f_4 * pc_x[k] * gpp_65[k];

        t_129[k] = f_6 * fpp_37[k]
                   + f_2 * gps0_21[k]
                   - f_3 * gps1_21[k]
                   + f_4 * pc_y[k] * gpp_64[k];

        t_130[k] = f_6 * fpp_38[k]
                   + f_4 * pc_y[k] * gpp_65[k];

        t_131[k] = f_1 * fpp_29[k]
                   + f_2 * gps0_21[k]
                   - f_3 * gps1_21[k]
                   + f_4 * pc_z[k] * gpp_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pa_z, pc_x, pc_z, fpd0_60, \
                         fpd0_135, fpp_67, fpp_68, fpd1_60, fpd1_135, gpp_67, \
                         gpp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_z[k] * fpd0_60[k]
                   - f_5 * pc_z[k] * fpd1_60[k];

        t_133[k] = f_1 * fpp_67[k]
                   + f_4 * pc_x[k] * gpp_67[k];

        t_134[k] = f_1 * fpp_68[k]
                   + f_4 * pc_x[k] * gpp_68[k];

        t_135[k] = pa_x[k] * fpd0_135[k]
                   - f_5 * pc_x[k] * fpd1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, pc_x, fpd0_136, fpd0_137, fpp_69, \
                         fpp_70, fpd1_136, fpd1_137, gps0_23, gps1_23, gpp_69, \
                         gpp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * fpd0_136[k]
                   - f_5 * pc_x[k] * fpd1_136[k];

        t_137[k] = pa_x[k] * fpd0_137[k]
                   - f_5 * pc_x[k] * fpd1_137[k];

        t_138[k] = f_1 * fpp_69[k]
                   + f_2 * gps0_23[k]
                   - f_3 * gps1_23[k]
                   + f_4 * pc_x[k] * gpp_69[k];

        t_139[k] = f_1 * fpp_70[k]
                   + f_4 * pc_x[k] * gpp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pc_x, pc_y, fpd0_141, fpd0_143, \
                         fpp_44, fpp_71, fpd1_141, fpd1_143, gpp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_1 * fpp_71[k]
                   + f_4 * pc_x[k] * gpp_71[k];

        t_141[k] = pa_x[k] * fpd0_141[k]
                   - f_5 * pc_x[k] * fpd1_141[k];

        t_142[k] = f_6 * fpp_44[k]
                   + f_4 * pc_y[k] * gpp_71[k];

        t_143[k] = pa_x[k] * fpd0_143[k]
                   - f_5 * pc_x[k] * fpd1_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pa_y, pc_x, pc_y, fpd0_90, fpd0_92, fpp_73, \
                         fpd1_90, fpd1_92, gsp_25, gpp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pa_y[k] * fpd0_90[k]
                   - f_5 * pc_y[k] * fpd1_90[k];

        t_145[k] = f_1 * fpp_73[k]
                   + f_1 * gsp_25[k]
                   + f_4 * pc_x[k] * gpp_73[k];

        t_146[k] = pa_y[k] * fpd0_92[k]
                   - f_5 * pc_y[k] * fpd1_92[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, fpp_38, fpp_46, fpp_47, gps0_24, \
                         gps1_24, gpp_73, gpp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_1 * fpp_46[k]
                   + f_2 * gps0_24[k]
                   - f_3 * gps1_24[k]
                   + f_4 * pc_y[k] * gpp_73[k];

        t_148[k] = f_1 * fpp_47[k]
                   + f_4 * pc_y[k] * gpp_74[k];

        t_149[k] = f_6 * fpp_38[k]
                   + f_2 * gps0_24[k]
                   - f_3 * gps1_24[k]
                   + f_4 * pc_z[k] * gpp_74[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_x, pc_x, fpd0_153, fpp_75, fpp_76, \
                         fpp_77, fpd1_153, gps0_25, gps1_25, gpp_75, gpp_76, \
                         gpp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * fpp_75[k]
                   + f_2 * gps0_25[k]
                   - f_3 * gps1_25[k]
                   + f_4 * pc_x[k] * gpp_75[k];

        t_151[k] = f_1 * fpp_76[k]
                   + f_4 * pc_x[k] * gpp_76[k];

        t_152[k] = f_1 * fpp_77[k]
                   + f_4 * pc_x[k] * gpp_77[k];

        t_153[k] = pa_x[k] * fpd0_153[k]
                   - f_5 * pc_x[k] * fpd1_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pa_y, pc_x, pc_y, fpd0_102, \
                         fpd0_154, fpd0_155, fpp_79, fpd1_102, fpd1_154, fpd1_155, \
                         gpp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_x[k] * fpd0_154[k]
                   - f_5 * pc_x[k] * fpd1_154[k];

        t_155[k] = pa_x[k] * fpd0_155[k]
                   - f_5 * pc_x[k] * fpd1_155[k];

        t_156[k] = pa_y[k] * fpd0_102[k]
                   - f_5 * pc_y[k] * fpd1_102[k];

        t_157[k] = f_1 * fpp_79[k]
                   + f_4 * pc_x[k] * gpp_79[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_x, pc_x, pc_y, fpd0_159, fpd0_161, \
                         fpp_53, fpp_80, fpd1_159, fpd1_161, gpp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_1 * fpp_80[k]
                   + f_4 * pc_x[k] * gpp_80[k];

        t_159[k] = pa_x[k] * fpd0_159[k]
                   - f_5 * pc_x[k] * fpd1_159[k];

        t_160[k] = f_1 * fpp_53[k]
                   + f_4 * pc_y[k] * gpp_80[k];

        t_161[k] = pa_x[k] * fpd0_161[k]
                   - f_5 * pc_x[k] * fpd1_161[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pc_x, pc_y, fpp_81, fpp_83, \
                         gsp_27, gsp_29, gps0_27, gps1_27, gpp_81, gpp_82, \
                         gpp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_1 * fpp_81[k]
                   + f_1 * gsp_27[k]
                   + f_2 * gps0_27[k]
                   - f_3 * gps1_27[k]
                   + f_4 * pc_x[k] * gpp_81[k];

        t_163[k] = f_4 * pc_y[k] * gpp_81[k];

        t_164[k] = f_1 * fpp_83[k]
                   + f_1 * gsp_29[k]
                   + f_4 * pc_x[k] * gpp_83[k];

        t_165[k] = f_2 * gps0_27[k]
                   - f_3 * gps1_27[k]
                   + f_4 * pc_y[k] * gpp_82[k];

        t_166[k] = f_4 * pc_y[k] * gpp_83[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pb_y, pc_y, pc_z, fpp_47, gsd0_54, gsp_27, \
                         gsd1_54, gps0_27, gps1_27, gpp_83, gpp_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_7 * fpp_47[k]
                   + f_2 * gps0_27[k]
                   - f_3 * gps1_27[k]
                   + f_4 * pc_z[k] * gpp_83[k];

        t_168[k] = pb_y[k] * gsd0_54[k]
                   - f_5 * pc_y[k] * gsd1_54[k];

        t_169[k] = f_1 * gsp_27[k]
                   + f_4 * pc_y[k] * gpp_84[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, pc_x, pc_y, fpd0_171, fpd0_173, \
                         fpp_86, fpd1_171, fpd1_173, gsp_29, gpp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_1 * fpp_86[k]
                   + f_4 * pc_x[k] * gpp_86[k];

        t_171[k] = pa_x[k] * fpd0_171[k]
                   - f_5 * pc_x[k] * fpd1_171[k];

        t_172[k] = f_1 * gsp_29[k]
                   + f_4 * pc_y[k] * gpp_86[k];

        t_173[k] = pa_x[k] * fpd0_173[k]
                   - f_5 * pc_x[k] * fpd1_173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pa_x, pc_x, pc_y, fpd0_174, \
                         fpd0_177, fpp_87, fpp_89, fpd1_174, fpd1_177, gpp_87, \
                         gpp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_x[k] * fpd0_174[k]
                   + f_6 * fpp_87[k]
                   - f_5 * pc_x[k] * fpd1_174[k];

        t_175[k] = f_4 * pc_y[k] * gpp_87[k];

        t_176[k] = f_1 * fpp_89[k]
                   + f_4 * pc_x[k] * gpp_89[k];

        t_177[k] = pa_x[k] * fpd0_177[k]
                   - f_5 * pc_x[k] * fpd1_177[k];

        t_178[k] = f_4 * pc_y[k] * gpp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pa_x, pb_x, pc_x, fpd0_179, fpd1_179, \
                         gsd0_60, gsp_30, gsp_31, gsp_32, gsd1_60, gpp_91, \
                         gpp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_x[k] * fpd0_179[k]
                   - f_5 * pc_x[k] * fpd1_179[k];

        t_180[k] = pb_x[k] * gsd0_60[k]
                   + f_6 * gsp_30[k]
                   - f_5 * pc_x[k] * gsd1_60[k];

        t_181[k] = f_1 * gsp_31[k]
                   + f_4 * pc_x[k] * gpp_91[k];

        t_182[k] = f_1 * gsp_32[k]
                   + f_4 * pc_x[k] * gpp_92[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pb_x, pc_x, pc_z, gsd0_63, gsd0_65, \
                         gsd1_63, gsd1_65, gps0_31, gps1_31, gpp_91, \
                         gpp_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pb_x[k] * gsd0_63[k]
                   - f_5 * pc_x[k] * gsd1_63[k];

        t_184[k] = f_4 * pc_z[k] * gpp_91[k];

        t_185[k] = pb_x[k] * gsd0_65[k]
                   - f_5 * pc_x[k] * gsd1_65[k];

        t_186[k] = f_2 * gps0_31[k]
                   - f_3 * gps1_31[k]
                   + f_4 * pc_x[k] * gpp_93[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, pc_x, pc_y, pc_z, fpp_58, gsp_31, \
                         gps0_31, gps1_31, gpp_94, gpp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_4 * pc_x[k] * gpp_94[k];

        t_188[k] = f_4 * pc_x[k] * gpp_95[k];

        t_189[k] = f_0 * fpp_58[k]
                   + f_1 * gsp_31[k]
                   + f_2 * gps0_31[k]
                   - f_3 * gps1_31[k]
                   + f_4 * pc_y[k] * gpp_94[k];

        t_190[k] = f_4 * pc_z[k] * gpp_94[k];

        t_191[k] = f_2 * gps0_31[k]
                   - f_3 * gps1_31[k]
                   + f_4 * pc_z[k] * gpp_95[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, pb_z, pc_x, pc_z, gsd0_60, \
                         gsd0_63, gsp_31, gsd1_60, gsd1_63, gpp_97, \
                         gpp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_z[k] * gsd0_60[k]
                   - f_5 * pc_z[k] * gsd1_60[k];

        t_193[k] = f_4 * pc_x[k] * gpp_97[k];

        t_194[k] = f_4 * pc_x[k] * gpp_98[k];

        t_195[k] = pb_z[k] * gsd0_63[k]
                   - f_5 * pc_z[k] * gsd1_63[k];

        t_196[k] = f_1 * gsp_31[k]
                   + f_4 * pc_z[k] * gpp_97[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_z, pb_z, pc_x, pc_z, fpd0_108, fpd1_108, \
                         gsd0_65, gsp_32, gsp_34, gsd1_65, gpp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pb_z[k] * gsd0_65[k]
                   + f_6 * gsp_32[k]
                   - f_5 * pc_z[k] * gsd1_65[k];

        t_198[k] = pa_z[k] * fpd0_108[k]
                   - f_5 * pc_z[k] * fpd1_108[k];

        t_199[k] = f_1 * gsp_34[k]
                   + f_4 * pc_x[k] * gpp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_z, pb_x, pc_x, pc_y, pc_z, fpd0_111, \
                         fpp_65, fpd1_111, gsd0_71, gsp_35, gsd1_71, \
                         gpp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * gsp_35[k]
                   + f_4 * pc_x[k] * gpp_101[k];

        t_201[k] = pa_z[k] * fpd0_111[k]
                   - f_5 * pc_z[k] * fpd1_111[k];

        t_202[k] = f_7 * fpp_65[k]
                   + f_4 * pc_y[k] * gpp_101[k];

        t_203[k] = pb_x[k] * gsd0_71[k]
                   - f_5 * pc_x[k] * gsd1_71[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_z, pc_x, pc_z, fpd0_114, fpd0_117, \
                         fpd1_114, fpd1_117, gpp_103, gpp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_z[k] * fpd0_114[k]
                   - f_5 * pc_z[k] * fpd1_114[k];

        t_205[k] = f_4 * pc_x[k] * gpp_103[k];

        t_206[k] = f_4 * pc_x[k] * gpp_104[k];

        t_207[k] = pa_z[k] * fpd0_117[k]
                   - f_5 * pc_z[k] * fpd1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pc_x, pc_y, pc_z, fpp_59, fpp_68, gsp_35, \
                         gps0_34, gps0_35, gps1_34, gps1_35, gpp_104, \
                         gpp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * fpp_68[k]
                   + f_1 * gsp_35[k]
                   + f_4 * pc_y[k] * gpp_104[k];

        t_209[k] = f_1 * fpp_59[k]
                   + f_2 * gps0_34[k]
                   - f_3 * gps1_34[k]
                   + f_4 * pc_z[k] * gpp_104[k];

        t_210[k] = f_2 * gps0_35[k]
                   - f_3 * gps1_35[k]
                   + f_4 * pc_x[k] * gpp_105[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pc_x, pc_y, fpp_70, fpp_71, gps0_35, \
                         gps1_35, gpp_106, gpp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_4 * pc_x[k] * gpp_106[k];

        t_212[k] = f_4 * pc_x[k] * gpp_107[k];

        t_213[k] = f_7 * fpp_70[k]
                   + f_2 * gps0_35[k]
                   - f_3 * gps1_35[k]
                   + f_4 * pc_y[k] * gpp_106[k];

        t_214[k] = f_7 * fpp_71[k]
                   + f_4 * pc_y[k] * gpp_107[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pa_y, pb_x, pc_x, pc_y, dpd0_89, dpd1_89, \
                         fpd0_143, fpd1_143, gsd0_72, gsp_36, gsp_37, gsd1_72, \
                         gpp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_8 * dpd0_89[k]
                   - f_9 * dpd1_89[k]
                   + pa_y[k] * fpd0_143[k]
                   - f_5 * pc_y[k] * fpd1_143[k];

        t_216[k] = pb_x[k] * gsd0_72[k]
                   + f_6 * gsp_36[k]
                   - f_5 * pc_x[k] * gsd1_72[k];

        t_217[k] = f_1 * gsp_37[k]
                   + f_4 * pc_x[k] * gpp_109[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pb_x, pc_x, pc_y, fpp_74, gsd0_75, \
                         gsd0_77, gsp_38, gsd1_75, gsd1_77, gpp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_1 * gsp_38[k]
                   + f_4 * pc_x[k] * gpp_110[k];

        t_219[k] = pb_x[k] * gsd0_75[k]
                   - f_5 * pc_x[k] * gsd1_75[k];

        t_220[k] = f_6 * fpp_74[k]
                   + f_4 * pc_y[k] * gpp_110[k];

        t_221[k] = pb_x[k] * gsd0_77[k]
                   - f_5 * pc_x[k] * gsd1_77[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pa_z, pc_x, pc_z, dpd0_63, dpd1_63, \
                         fpd0_135, fpd1_135, gps0_37, gps1_37, gpp_111, gpp_112, \
                         gpp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_2 * gps0_37[k]
                   - f_3 * gps1_37[k]
                   + f_4 * pc_x[k] * gpp_111[k];

        t_223[k] = f_4 * pc_x[k] * gpp_112[k];

        t_224[k] = f_4 * pc_x[k] * gpp_113[k];

        t_225[k] = f_10 * dpd0_63[k]
                   - f_11 * dpd1_63[k]
                   + pa_z[k] * fpd0_135[k]
                   - f_5 * pc_z[k] * fpd1_135[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_y, pc_z, fpp_68, fpp_77, gsp_38, \
                         gps0_37, gps0_38, gps1_37, gps1_38, gpp_113, \
                         gpp_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_6 * fpp_77[k]
                   + f_1 * gsp_38[k]
                   + f_4 * pc_y[k] * gpp_113[k];

        t_227[k] = f_6 * fpp_68[k]
                   + f_2 * gps0_37[k]
                   - f_3 * gps1_37[k]
                   + f_4 * pc_z[k] * gpp_113[k];

        t_228[k] = f_2 * gps0_38[k]
                   - f_3 * gps1_38[k]
                   + f_4 * pc_x[k] * gpp_114[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pc_x, pc_y, fpp_79, fpp_80, gps0_38, \
                         gps1_38, gpp_115, gpp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_4 * pc_x[k] * gpp_115[k];

        t_230[k] = f_4 * pc_x[k] * gpp_116[k];

        t_231[k] = f_6 * fpp_79[k]
                   + f_2 * gps0_38[k]
                   - f_3 * gps1_38[k]
                   + f_4 * pc_y[k] * gpp_115[k];

        t_232[k] = f_6 * fpp_80[k]
                   + f_4 * pc_y[k] * gpp_116[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_x, pc_y, dpd0_107, dpd1_107, fpd0_161, \
                         fpd0_162, fpd1_161, fpd1_162, gsp_40, \
                         gpp_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_10 * dpd0_107[k]
                   - f_11 * dpd1_107[k]
                   + pa_y[k] * fpd0_161[k]
                   - f_5 * pc_y[k] * fpd1_161[k];

        t_234[k] = pa_y[k] * fpd0_162[k]
                   - f_5 * pc_y[k] * fpd1_162[k];

        t_235[k] = f_1 * gsp_40[k]
                   + f_4 * pc_x[k] * gpp_118[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pb_x, pc_x, pc_y, fpd0_167, fpp_83, \
                         fpd1_167, gsd0_81, gsp_41, gsd1_81, gpp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_1 * gsp_41[k]
                   + f_4 * pc_x[k] * gpp_119[k];

        t_237[k] = pb_x[k] * gsd0_81[k]
                   - f_5 * pc_x[k] * gsd1_81[k];

        t_238[k] = f_1 * fpp_83[k]
                   + f_4 * pc_y[k] * gpp_119[k];

        t_239[k] = pa_y[k] * fpd0_167[k]
                   - f_5 * pc_y[k] * fpd1_167[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_y, fpp_85, fpp_86, \
                         gsp_40, gsp_41, gps0_40, gps1_40, gpp_120, gpp_121, \
                         gpp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_2 * gps0_40[k]
                   - f_3 * gps1_40[k]
                   + f_4 * pc_x[k] * gpp_120[k];

        t_241[k] = f_4 * pc_x[k] * gpp_121[k];

        t_242[k] = f_4 * pc_x[k] * gpp_122[k];

        t_243[k] = f_1 * fpp_85[k]
                   + f_1 * gsp_40[k]
                   + f_2 * gps0_40[k]
                   - f_3 * gps1_40[k]
                   + f_4 * pc_y[k] * gpp_121[k];

        t_244[k] = f_1 * fpp_86[k]
                   + f_1 * gsp_41[k]
                   + f_4 * pc_y[k] * gpp_122[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_y, pc_x, pc_y, pc_z, fpd0_174, fpp_77, \
                         fpd1_174, gps0_40, gps1_40, gpp_122, gpp_124, \
                         gpp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_7 * fpp_77[k]
                   + f_2 * gps0_40[k]
                   - f_3 * gps1_40[k]
                   + f_4 * pc_z[k] * gpp_122[k];

        t_246[k] = pa_y[k] * fpd0_174[k]
                   - f_5 * pc_y[k] * fpd1_174[k];

        t_247[k] = f_4 * pc_x[k] * gpp_124[k];

        t_248[k] = f_4 * pc_x[k] * gpp_125[k];
    }
}

static auto
compute_prim_gpd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fpd0, const size_t fpp,
                                                          const size_t fpd1, const size_t gsd0,
                                                          const size_t gsp, const size_t gsd1,
                                                          const size_t gps0, const size_t gps1,
                                                          const size_t gpp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpd0_177 = buffer.data(fpd0 + 177);
    const auto *fpd0_179 = buffer.data(fpd0 + 179);

    const auto *fpp_88 = buffer.data(fpp + 88);
    const auto *fpp_89 = buffer.data(fpp + 89);

    const auto *fpd1_177 = buffer.data(fpd1 + 177);
    const auto *fpd1_179 = buffer.data(fpd1 + 179);

    const auto *gsd0_84 = buffer.data(gsd0 + 84);
    const auto *gsd0_87 = buffer.data(gsd0 + 87);
    const auto *gsd0_89 = buffer.data(gsd0 + 89);

    const auto *gsp_42 = buffer.data(gsp + 42);
    const auto *gsp_43 = buffer.data(gsp + 43);
    const auto *gsp_44 = buffer.data(gsp + 44);

    const auto *gsd1_84 = buffer.data(gsd1 + 84);
    const auto *gsd1_87 = buffer.data(gsd1 + 87);
    const auto *gsd1_89 = buffer.data(gsd1 + 89);

    const auto *gps0_44 = buffer.data(gps0 + 44);

    const auto *gps1_44 = buffer.data(gps1 + 44);

    const auto *gpp_125 = buffer.data(gpp + 125);
    const auto *gpp_127 = buffer.data(gpp + 127);
    const auto *gpp_128 = buffer.data(gpp + 128);
    const auto *gpp_130 = buffer.data(gpp + 130);
    const auto *gpp_131 = buffer.data(gpp + 131);
    const auto *gpp_132 = buffer.data(gpp + 132);
    const auto *gpp_133 = buffer.data(gpp + 133);
    const auto *gpp_134 = buffer.data(gpp + 134);

#pragma omp simd aligned(t_249, t_250, t_251, pa_y, pc_y, fpd0_177, fpd0_179, fpp_88, fpp_89, \
                         fpd1_177, fpd1_179, gpp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pa_y[k] * fpd0_177[k]
                   + f_6 * fpp_88[k]
                   - f_5 * pc_y[k] * fpd1_177[k];

        t_250[k] = f_1 * fpp_89[k]
                   + f_4 * pc_y[k] * gpp_125[k];

        t_251[k] = pa_y[k] * fpd0_179[k]
                   - f_5 * pc_y[k] * fpd1_179[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, pc_x, gsd0_84, gsd0_87, gsp_42, \
                         gsp_43, gsp_44, gsd1_84, gsd1_87, gpp_127, \
                         gpp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pb_x[k] * gsd0_84[k]
                   + f_6 * gsp_42[k]
                   - f_5 * pc_x[k] * gsd1_84[k];

        t_253[k] = f_1 * gsp_43[k]
                   + f_4 * pc_x[k] * gpp_127[k];

        t_254[k] = f_1 * gsp_44[k]
                   + f_4 * pc_x[k] * gpp_128[k];

        t_255[k] = pb_x[k] * gsd0_87[k]
                   - f_5 * pc_x[k] * gsd1_87[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, pc_x, pc_y, gsd0_84, \
                         gsd0_89, gsd1_84, gsd1_89, gpp_128, gpp_130, \
                         gpp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_4 * pc_y[k] * gpp_128[k];

        t_257[k] = pb_x[k] * gsd0_89[k]
                   - f_5 * pc_x[k] * gsd1_89[k];

        t_258[k] = pb_y[k] * gsd0_84[k]
                   - f_5 * pc_y[k] * gsd1_84[k];

        t_259[k] = f_4 * pc_x[k] * gpp_130[k];

        t_260[k] = f_4 * pc_x[k] * gpp_131[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pb_y, pc_y, gsd0_87, gsd0_89, gsp_43, gsp_44, \
                         gsd1_87, gsd1_89, gpp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pb_y[k] * gsd0_87[k]
                   + f_6 * gsp_43[k]
                   - f_5 * pc_y[k] * gsd1_87[k];

        t_262[k] = f_1 * gsp_44[k]
                   + f_4 * pc_y[k] * gpp_131[k];

        t_263[k] = pb_y[k] * gsd0_89[k]
                   - f_5 * pc_y[k] * gsd1_89[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, t_269, pc_x, pc_y, pc_z, fpp_89, \
                         gsp_44, gps0_44, gps1_44, gpp_132, gpp_133, \
                         gpp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_2 * gps0_44[k]
                   - f_3 * gps1_44[k]
                   + f_4 * pc_x[k] * gpp_132[k];

        t_265[k] = f_4 * pc_x[k] * gpp_133[k];

        t_266[k] = f_4 * pc_x[k] * gpp_134[k];

        t_267[k] = f_2 * gps0_44[k]
                   - f_3 * gps1_44[k]
                   + f_4 * pc_y[k] * gpp_133[k];

        t_268[k] = f_4 * pc_y[k] * gpp_134[k];

        t_269[k] = f_0 * fpp_89[k]
                   + f_1 * gsp_44[k]
                   + f_2 * gps0_44[k]
                   - f_3 * gps1_44[k]
                   + f_4 * pc_z[k] * gpp_134[k];
    }
}

auto
compute_prim_gpd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dpd0,
                                                   const size_t dpd1, const size_t fpd0,
                                                   const size_t fpp, const size_t fpd1,
                                                   const size_t gsd0, const size_t gsp,
                                                   const size_t gsd1, const size_t gps0,
                                                   const size_t gps1, const size_t gpp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gpd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, dpd0,
                                                              dpd1, fpd0, fpp, fpd1, gsd0, gsp,
                                                              gsd1, gps0, gps1, gpp, ncols,
                                                              gamma, p, q);

    compute_prim_gpd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, dpd0,
                                                              dpd1, fpd0, fpp, fpd1, gsd0, gsp,
                                                              gsd1, gps0, gps1, gpp, ncols,
                                                              gamma, p, q);

    compute_prim_gpd_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, fpd0,
                                                              fpp, fpd1, gsd0, gsp, gsd1, gps0,
                                                              gps1, gpp, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
