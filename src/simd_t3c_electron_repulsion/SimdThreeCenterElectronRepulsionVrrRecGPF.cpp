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


#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gpf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t fpf0, const size_t fpd,
                                                          const size_t fpf1, const size_t gsf0,
                                                          const size_t gsd, const size_t gsf1,
                                                          const size_t gpp0, const size_t gpp1,
                                                          const size_t gpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 1.0 / p;
    const auto f_8 = gamma / (p * q);
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);
    const auto f_11 = 1.0 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpf0_0 = buffer.data(dpf0 + 0);
    const auto *dpf0_46 = buffer.data(dpf0 + 46);
    const auto *dpf0_89 = buffer.data(dpf0 + 89);
    const auto *dpf0_106 = buffer.data(dpf0 + 106);

    const auto *dpf1_0 = buffer.data(dpf1 + 0);
    const auto *dpf1_46 = buffer.data(dpf1 + 46);
    const auto *dpf1_89 = buffer.data(dpf1 + 89);
    const auto *dpf1_106 = buffer.data(dpf1 + 106);

    const auto *fpf0_0 = buffer.data(fpf0 + 0);
    const auto *fpf0_3 = buffer.data(fpf0 + 3);
    const auto *fpf0_5 = buffer.data(fpf0 + 5);
    const auto *fpf0_6 = buffer.data(fpf0 + 6);
    const auto *fpf0_9 = buffer.data(fpf0 + 9);
    const auto *fpf0_10 = buffer.data(fpf0 + 10);
    const auto *fpf0_16 = buffer.data(fpf0 + 16);
    const auto *fpf0_20 = buffer.data(fpf0 + 20);
    const auto *fpf0_29 = buffer.data(fpf0 + 29);
    const auto *fpf0_30 = buffer.data(fpf0 + 30);
    const auto *fpf0_33 = buffer.data(fpf0 + 33);
    const auto *fpf0_46 = buffer.data(fpf0 + 46);
    const auto *fpf0_60 = buffer.data(fpf0 + 60);
    const auto *fpf0_65 = buffer.data(fpf0 + 65);
    const auto *fpf0_89 = buffer.data(fpf0 + 89);
    const auto *fpf0_106 = buffer.data(fpf0 + 106);

    const auto *fpd_0 = buffer.data(fpd + 0);
    const auto *fpd_3 = buffer.data(fpd + 3);
    const auto *fpd_5 = buffer.data(fpd + 5);
    const auto *fpd_6 = buffer.data(fpd + 6);
    const auto *fpd_9 = buffer.data(fpd + 9);
    const auto *fpd_11 = buffer.data(fpd + 11);
    const auto *fpd_12 = buffer.data(fpd + 12);
    const auto *fpd_15 = buffer.data(fpd + 15);
    const auto *fpd_17 = buffer.data(fpd + 17);
    const auto *fpd_18 = buffer.data(fpd + 18);
    const auto *fpd_21 = buffer.data(fpd + 21);
    const auto *fpd_23 = buffer.data(fpd + 23);
    const auto *fpd_24 = buffer.data(fpd + 24);
    const auto *fpd_27 = buffer.data(fpd + 27);
    const auto *fpd_29 = buffer.data(fpd + 29);
    const auto *fpd_30 = buffer.data(fpd + 30);
    const auto *fpd_33 = buffer.data(fpd + 33);
    const auto *fpd_35 = buffer.data(fpd + 35);
    const auto *fpd_36 = buffer.data(fpd + 36);
    const auto *fpd_41 = buffer.data(fpd + 41);
    const auto *fpd_45 = buffer.data(fpd + 45);
    const auto *fpd_47 = buffer.data(fpd + 47);
    const auto *fpd_48 = buffer.data(fpd + 48);
    const auto *fpd_51 = buffer.data(fpd + 51);
    const auto *fpd_53 = buffer.data(fpd + 53);
    const auto *fpd_57 = buffer.data(fpd + 57);
    const auto *fpd_59 = buffer.data(fpd + 59);
    const auto *fpd_60 = buffer.data(fpd + 60);
    const auto *fpd_63 = buffer.data(fpd + 63);
    const auto *fpd_65 = buffer.data(fpd + 65);
    const auto *fpd_69 = buffer.data(fpd + 69);
    const auto *fpd_71 = buffer.data(fpd + 71);
    const auto *fpd_76 = buffer.data(fpd + 76);

    const auto *fpf1_0 = buffer.data(fpf1 + 0);
    const auto *fpf1_3 = buffer.data(fpf1 + 3);
    const auto *fpf1_5 = buffer.data(fpf1 + 5);
    const auto *fpf1_6 = buffer.data(fpf1 + 6);
    const auto *fpf1_9 = buffer.data(fpf1 + 9);
    const auto *fpf1_10 = buffer.data(fpf1 + 10);
    const auto *fpf1_16 = buffer.data(fpf1 + 16);
    const auto *fpf1_20 = buffer.data(fpf1 + 20);
    const auto *fpf1_29 = buffer.data(fpf1 + 29);
    const auto *fpf1_30 = buffer.data(fpf1 + 30);
    const auto *fpf1_33 = buffer.data(fpf1 + 33);
    const auto *fpf1_46 = buffer.data(fpf1 + 46);
    const auto *fpf1_60 = buffer.data(fpf1 + 60);
    const auto *fpf1_65 = buffer.data(fpf1 + 65);
    const auto *fpf1_89 = buffer.data(fpf1 + 89);
    const auto *fpf1_106 = buffer.data(fpf1 + 106);

    const auto *gsf0_0 = buffer.data(gsf0 + 0);
    const auto *gsf0_6 = buffer.data(gsf0 + 6);
    const auto *gsf0_9 = buffer.data(gsf0 + 9);
    const auto *gsf0_16 = buffer.data(gsf0 + 16);
    const auto *gsf0_27 = buffer.data(gsf0 + 27);
    const auto *gsf0_29 = buffer.data(gsf0 + 29);
    const auto *gsf0_30 = buffer.data(gsf0 + 30);
    const auto *gsf0_36 = buffer.data(gsf0 + 36);
    const auto *gsf0_39 = buffer.data(gsf0 + 39);

    const auto *gsd_0 = buffer.data(gsd + 0);
    const auto *gsd_2 = buffer.data(gsd + 2);
    const auto *gsd_3 = buffer.data(gsd + 3);
    const auto *gsd_5 = buffer.data(gsd + 5);
    const auto *gsd_6 = buffer.data(gsd + 6);
    const auto *gsd_7 = buffer.data(gsd + 7);
    const auto *gsd_9 = buffer.data(gsd + 9);
    const auto *gsd_11 = buffer.data(gsd + 11);
    const auto *gsd_12 = buffer.data(gsd + 12);
    const auto *gsd_14 = buffer.data(gsd + 14);
    const auto *gsd_16 = buffer.data(gsd + 16);
    const auto *gsd_17 = buffer.data(gsd + 17);
    const auto *gsd_18 = buffer.data(gsd + 18);
    const auto *gsd_19 = buffer.data(gsd + 19);
    const auto *gsd_21 = buffer.data(gsd + 21);
    const auto *gsd_23 = buffer.data(gsd + 23);
    const auto *gsd_28 = buffer.data(gsd + 28);

    const auto *gsf1_0 = buffer.data(gsf1 + 0);
    const auto *gsf1_6 = buffer.data(gsf1 + 6);
    const auto *gsf1_9 = buffer.data(gsf1 + 9);
    const auto *gsf1_16 = buffer.data(gsf1 + 16);
    const auto *gsf1_27 = buffer.data(gsf1 + 27);
    const auto *gsf1_29 = buffer.data(gsf1 + 29);
    const auto *gsf1_30 = buffer.data(gsf1 + 30);
    const auto *gsf1_36 = buffer.data(gsf1 + 36);
    const auto *gsf1_39 = buffer.data(gsf1 + 39);

    const auto *gpp0_0 = buffer.data(gpp0 + 0);
    const auto *gpp0_1 = buffer.data(gpp0 + 1);
    const auto *gpp0_2 = buffer.data(gpp0 + 2);
    const auto *gpp0_10 = buffer.data(gpp0 + 10);
    const auto *gpp0_12 = buffer.data(gpp0 + 12);
    const auto *gpp0_14 = buffer.data(gpp0 + 14);
    const auto *gpp0_20 = buffer.data(gpp0 + 20);
    const auto *gpp0_24 = buffer.data(gpp0 + 24);
    const auto *gpp0_25 = buffer.data(gpp0 + 25);
    const auto *gpp0_26 = buffer.data(gpp0 + 26);
    const auto *gpp0_28 = buffer.data(gpp0 + 28);
    const auto *gpp0_29 = buffer.data(gpp0 + 29);
    const auto *gpp0_30 = buffer.data(gpp0 + 30);
    const auto *gpp0_32 = buffer.data(gpp0 + 32);

    const auto *gpp1_0 = buffer.data(gpp1 + 0);
    const auto *gpp1_1 = buffer.data(gpp1 + 1);
    const auto *gpp1_2 = buffer.data(gpp1 + 2);
    const auto *gpp1_10 = buffer.data(gpp1 + 10);
    const auto *gpp1_12 = buffer.data(gpp1 + 12);
    const auto *gpp1_14 = buffer.data(gpp1 + 14);
    const auto *gpp1_20 = buffer.data(gpp1 + 20);
    const auto *gpp1_24 = buffer.data(gpp1 + 24);
    const auto *gpp1_25 = buffer.data(gpp1 + 25);
    const auto *gpp1_26 = buffer.data(gpp1 + 26);
    const auto *gpp1_28 = buffer.data(gpp1 + 28);
    const auto *gpp1_29 = buffer.data(gpp1 + 29);
    const auto *gpp1_30 = buffer.data(gpp1 + 30);
    const auto *gpp1_32 = buffer.data(gpp1 + 32);

    const auto *gpd_0 = buffer.data(gpd + 0);
    const auto *gpd_2 = buffer.data(gpd + 2);
    const auto *gpd_3 = buffer.data(gpd + 3);
    const auto *gpd_5 = buffer.data(gpd + 5);
    const auto *gpd_6 = buffer.data(gpd + 6);
    const auto *gpd_8 = buffer.data(gpd + 8);
    const auto *gpd_9 = buffer.data(gpd + 9);
    const auto *gpd_11 = buffer.data(gpd + 11);
    const auto *gpd_12 = buffer.data(gpd + 12);
    const auto *gpd_14 = buffer.data(gpd + 14);
    const auto *gpd_15 = buffer.data(gpd + 15);
    const auto *gpd_17 = buffer.data(gpd + 17);
    const auto *gpd_18 = buffer.data(gpd + 18);
    const auto *gpd_19 = buffer.data(gpd + 19);
    const auto *gpd_21 = buffer.data(gpd + 21);
    const auto *gpd_23 = buffer.data(gpd + 23);
    const auto *gpd_24 = buffer.data(gpd + 24);
    const auto *gpd_25 = buffer.data(gpd + 25);
    const auto *gpd_27 = buffer.data(gpd + 27);
    const auto *gpd_29 = buffer.data(gpd + 29);
    const auto *gpd_30 = buffer.data(gpd + 30);
    const auto *gpd_31 = buffer.data(gpd + 31);
    const auto *gpd_33 = buffer.data(gpd + 33);
    const auto *gpd_35 = buffer.data(gpd + 35);
    const auto *gpd_36 = buffer.data(gpd + 36);
    const auto *gpd_38 = buffer.data(gpd + 38);
    const auto *gpd_40 = buffer.data(gpd + 40);
    const auto *gpd_41 = buffer.data(gpd + 41);
    const auto *gpd_42 = buffer.data(gpd + 42);
    const auto *gpd_44 = buffer.data(gpd + 44);
    const auto *gpd_45 = buffer.data(gpd + 45);
    const auto *gpd_47 = buffer.data(gpd + 47);
    const auto *gpd_48 = buffer.data(gpd + 48);
    const auto *gpd_50 = buffer.data(gpd + 50);
    const auto *gpd_51 = buffer.data(gpd + 51);
    const auto *gpd_52 = buffer.data(gpd + 52);
    const auto *gpd_53 = buffer.data(gpd + 53);
    const auto *gpd_54 = buffer.data(gpd + 54);
    const auto *gpd_55 = buffer.data(gpd + 55);
    const auto *gpd_57 = buffer.data(gpd + 57);
    const auto *gpd_59 = buffer.data(gpd + 59);
    const auto *gpd_60 = buffer.data(gpd + 60);
    const auto *gpd_61 = buffer.data(gpd + 61);
    const auto *gpd_63 = buffer.data(gpd + 63);
    const auto *gpd_65 = buffer.data(gpd + 65);
    const auto *gpd_66 = buffer.data(gpd + 66);
    const auto *gpd_67 = buffer.data(gpd + 67);
    const auto *gpd_69 = buffer.data(gpd + 69);
    const auto *gpd_71 = buffer.data(gpd + 71);
    const auto *gpd_72 = buffer.data(gpd + 72);
    const auto *gpd_76 = buffer.data(gpd + 76);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, fpd_0, fpd_3, gsd_0, gsd_3, \
                         gpp0_0, gpp1_0, gpd_0, gpd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fpd_0[k]
                 + f_1 * gsd_0[k]
                 + f_2 * gpp0_0[k]
                 - f_3 * gpp1_0[k]
                 + f_4 * pc_x[k] * gpd_0[k];

        t_1[k] = f_4 * pc_y[k] * gpd_0[k];

        t_2[k] = f_4 * pc_z[k] * gpd_0[k];

        t_3[k] = f_0 * fpd_3[k]
                 + f_1 * gsd_3[k]
                 + f_4 * pc_x[k] * gpd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, fpd_5, gsd_5, gpp0_1, \
                         gpp1_1, gpd_2, gpd_3, gpd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * pc_y[k] * gpd_2[k];

        t_5[k] = f_0 * fpd_5[k]
                 + f_1 * gsd_5[k]
                 + f_4 * pc_x[k] * gpd_5[k];

        t_6[k] = f_2 * gpp0_1[k]
                 - f_3 * gpp1_1[k]
                 + f_4 * pc_y[k] * gpd_3[k];

        t_7[k] = f_4 * pc_z[k] * gpd_3[k];

        t_8[k] = f_4 * pc_y[k] * gpd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, gsf0_0, gsd_0, gsf1_0, \
                         gpp0_2, gpp1_2, gpd_5, gpd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * gpp0_2[k]
                 - f_3 * gpp1_2[k]
                 + f_4 * pc_z[k] * gpd_5[k];

        t_10[k] = pb_y[k] * gsf0_0[k]
                  - f_5 * pc_y[k] * gsf1_0[k];

        t_11[k] = f_1 * gsd_0[k]
                  + f_4 * pc_y[k] * gpd_6[k];

        t_12[k] = f_4 * pc_z[k] * gpd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pc_x, pc_y, fpd_9, fpd_11, gsf0_6, \
                         gsd_2, gsd_3, gsf1_6, gpd_8, gpd_9, gpd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * fpd_9[k]
                  + f_4 * pc_x[k] * gpd_9[k];

        t_14[k] = f_1 * gsd_2[k]
                  + f_4 * pc_y[k] * gpd_8[k];

        t_15[k] = f_0 * fpd_11[k]
                  + f_4 * pc_x[k] * gpd_11[k];

        t_16[k] = pb_y[k] * gsf0_6[k]
                  + f_6 * gsd_3[k]
                  - f_5 * pc_y[k] * gsf1_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, pc_y, pc_z, gsf0_0, gsf0_9, \
                         gsd_5, gsf1_0, gsf1_9, gpd_9, gpd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * pc_z[k] * gpd_9[k];

        t_18[k] = f_1 * gsd_5[k]
                  + f_4 * pc_y[k] * gpd_11[k];

        t_19[k] = pb_y[k] * gsf0_9[k]
                  - f_5 * pc_y[k] * gsf1_9[k];

        t_20[k] = pb_z[k] * gsf0_0[k]
                  - f_5 * pc_z[k] * gsf1_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, fpd_15, fpd_17, \
                         gsd_0, gpd_12, gpd_14, gpd_15, gpd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * pc_y[k] * gpd_12[k];

        t_22[k] = f_1 * gsd_0[k]
                  + f_4 * pc_z[k] * gpd_12[k];

        t_23[k] = f_0 * fpd_15[k]
                  + f_4 * pc_x[k] * gpd_15[k];

        t_24[k] = f_4 * pc_y[k] * gpd_14[k];

        t_25[k] = f_0 * fpd_17[k]
                  + f_4 * pc_x[k] * gpd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, gsf0_6, gsf0_9, gsd_3, \
                         gsd_5, gsf1_6, gsf1_9, gpd_15, gpd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * gsf0_6[k]
                  - f_5 * pc_z[k] * gsf1_6[k];

        t_27[k] = f_1 * gsd_3[k]
                  + f_4 * pc_z[k] * gpd_15[k];

        t_28[k] = f_4 * pc_y[k] * gpd_17[k];

        t_29[k] = pb_z[k] * gsf0_9[k]
                  + f_6 * gsd_5[k]
                  - f_5 * pc_z[k] * gsf1_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_x, pc_y, pc_z, fpf0_0, fpd_0, \
                         fpd_21, fpf1_0, gsd_9, gpd_18, gpd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * fpf0_0[k]
                  - f_5 * pc_y[k] * fpf1_0[k];

        t_31[k] = f_1 * fpd_0[k]
                  + f_4 * pc_y[k] * gpd_18[k];

        t_32[k] = f_4 * pc_z[k] * gpd_18[k];

        t_33[k] = f_6 * fpd_21[k]
                  + f_1 * gsd_9[k]
                  + f_4 * pc_x[k] * gpd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, fpf0_5, fpd_3, fpf1_5, \
                         gpp0_10, gpp1_10, gpd_19, gpd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_4 * pc_z[k] * gpd_19[k];

        t_35[k] = pa_y[k] * fpf0_5[k]
                  - f_5 * pc_y[k] * fpf1_5[k];

        t_36[k] = f_1 * fpd_3[k]
                  + f_2 * gpp0_10[k]
                  - f_3 * gpp1_10[k]
                  + f_4 * pc_y[k] * gpd_21[k];

        t_37[k] = f_4 * pc_z[k] * gpd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pc_x, pc_y, fpf0_9, fpd_5, fpd_24, fpf1_9, \
                         gpp0_12, gpp1_12, gpd_23, gpd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * fpd_5[k]
                  + f_4 * pc_y[k] * gpd_23[k];

        t_39[k] = pa_y[k] * fpf0_9[k]
                  - f_5 * pc_y[k] * fpf1_9[k];

        t_40[k] = f_6 * fpd_24[k]
                  + f_2 * gpp0_12[k]
                  - f_3 * gpp1_12[k]
                  + f_4 * pc_x[k] * gpd_24[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pc_x, pc_y, pc_z, fpd_6, fpd_27, \
                         fpd_29, gsd_6, gpd_24, gpd_25, gpd_27, \
                         gpd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * fpd_6[k]
                  + f_1 * gsd_6[k]
                  + f_4 * pc_y[k] * gpd_24[k];

        t_42[k] = f_4 * pc_z[k] * gpd_24[k];

        t_43[k] = f_6 * fpd_27[k]
                  + f_4 * pc_x[k] * gpd_27[k];

        t_44[k] = f_4 * pc_z[k] * gpd_25[k];

        t_45[k] = f_6 * fpd_29[k]
                  + f_4 * pc_x[k] * gpd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pc_x, pc_y, pc_z, dpf0_46, dpf1_46, fpf0_46, \
                         fpd_11, fpf1_46, gsd_11, gpd_27, gpd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_7 * dpf0_46[k]
                  - f_8 * dpf1_46[k]
                  + pa_x[k] * fpf0_46[k]
                  - f_5 * pc_x[k] * fpf1_46[k];

        t_47[k] = f_4 * pc_z[k] * gpd_27[k];

        t_48[k] = f_1 * fpd_11[k]
                  + f_1 * gsd_11[k]
                  + f_4 * pc_y[k] * gpd_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pc_y, pc_z, fpf0_20, fpd_12, fpf1_20, \
                         gsd_6, gpp0_14, gpp1_14, gpd_29, gpd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * gpp0_14[k]
                  - f_3 * gpp1_14[k]
                  + f_4 * pc_z[k] * gpd_29[k];

        t_50[k] = pa_y[k] * fpf0_20[k]
                  - f_5 * pc_y[k] * fpf1_20[k];

        t_51[k] = f_1 * fpd_12[k]
                  + f_4 * pc_y[k] * gpd_30[k];

        t_52[k] = f_1 * gsd_6[k]
                  + f_4 * pc_z[k] * gpd_30[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_z, pc_x, pc_z, fpd_33, fpd_35, gsf0_16, \
                         gsd_7, gsf1_16, gpd_31, gpd_33, gpd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_6 * fpd_33[k]
                  + f_4 * pc_x[k] * gpd_33[k];

        t_54[k] = f_1 * gsd_7[k]
                  + f_4 * pc_z[k] * gpd_31[k];

        t_55[k] = f_6 * fpd_35[k]
                  + f_4 * pc_x[k] * gpd_35[k];

        t_56[k] = pb_z[k] * gsf0_16[k]
                  - f_5 * pc_z[k] * gsf1_16[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_y, pa_z, pc_y, pc_z, fpf0_0, fpf0_29, \
                         fpd_17, fpf1_0, fpf1_29, gsd_9, gpd_33, \
                         gpd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * gsd_9[k]
                  + f_4 * pc_z[k] * gpd_33[k];

        t_58[k] = f_1 * fpd_17[k]
                  + f_4 * pc_y[k] * gpd_35[k];

        t_59[k] = pa_y[k] * fpf0_29[k]
                  - f_5 * pc_y[k] * fpf1_29[k];

        t_60[k] = pa_z[k] * fpf0_0[k]
                  - f_5 * pc_z[k] * fpf1_0[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_z, pc_y, pc_z, fpf0_3, fpd_0, fpf1_3, \
                         gpd_36, gpd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_4 * pc_y[k] * gpd_36[k];

        t_62[k] = f_1 * fpd_0[k]
                  + f_4 * pc_z[k] * gpd_36[k];

        t_63[k] = pa_z[k] * fpf0_3[k]
                  - f_5 * pc_z[k] * fpf1_3[k];

        t_64[k] = f_4 * pc_y[k] * gpd_38[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pc_x, pc_y, pc_z, fpf0_6, fpd_41, \
                         fpf1_6, gsd_17, gpp0_20, gpp1_20, gpd_40, \
                         gpd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_6 * fpd_41[k]
                  + f_1 * gsd_17[k]
                  + f_4 * pc_x[k] * gpd_41[k];

        t_66[k] = pa_z[k] * fpf0_6[k]
                  - f_5 * pc_z[k] * fpf1_6[k];

        t_67[k] = f_9 * gpp0_20[k]
                  - f_10 * gpp1_20[k]
                  + f_4 * pc_y[k] * gpd_40[k];

        t_68[k] = f_4 * pc_y[k] * gpd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, fpf0_10, fpd_5, fpd_6, \
                         fpf1_10, gsd_12, gpp0_20, gpp1_20, gpd_41, \
                         gpd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * fpd_5[k]
                  + f_2 * gpp0_20[k]
                  - f_3 * gpp1_20[k]
                  + f_4 * pc_z[k] * gpd_41[k];

        t_70[k] = pa_z[k] * fpf0_10[k]
                  - f_5 * pc_z[k] * fpf1_10[k];

        t_71[k] = f_1 * gsd_12[k]
                  + f_4 * pc_y[k] * gpd_42[k];

        t_72[k] = f_1 * fpd_6[k]
                  + f_4 * pc_z[k] * gpd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_y, pc_z, fpf0_16, fpd_45, \
                         fpd_47, fpf1_16, gsd_14, gpd_44, gpd_45, \
                         gpd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * fpd_45[k]
                  + f_4 * pc_x[k] * gpd_45[k];

        t_74[k] = f_1 * gsd_14[k]
                  + f_4 * pc_y[k] * gpd_44[k];

        t_75[k] = f_6 * fpd_47[k]
                  + f_4 * pc_x[k] * gpd_47[k];

        t_76[k] = pa_z[k] * fpf0_16[k]
                  - f_5 * pc_z[k] * fpf1_16[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, pc_y, gsf0_27, gsf0_29, gsd_16, gsd_17, \
                         gsf1_27, gsf1_29, gpd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_y[k] * gsf0_27[k]
                  + f_11 * gsd_16[k]
                  - f_5 * pc_y[k] * gsf1_27[k];

        t_78[k] = f_1 * gsd_17[k]
                  + f_4 * pc_y[k] * gpd_47[k];

        t_79[k] = pb_y[k] * gsf0_29[k]
                  - f_5 * pc_y[k] * gsf1_29[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, fpd_12, fpd_48, fpd_51, \
                         gsd_12, gpp0_24, gpp1_24, gpd_48, gpd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * fpd_48[k]
                  + f_2 * gpp0_24[k]
                  - f_3 * gpp1_24[k]
                  + f_4 * pc_x[k] * gpd_48[k];

        t_81[k] = f_4 * pc_y[k] * gpd_48[k];

        t_82[k] = f_1 * fpd_12[k]
                  + f_1 * gsd_12[k]
                  + f_4 * pc_z[k] * gpd_48[k];

        t_83[k] = f_6 * fpd_51[k]
                  + f_4 * pc_x[k] * gpd_51[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pc_x, pc_y, fpd_53, gpp0_25, gpp0_26, \
                         gpp1_25, gpp1_26, gpd_50, gpd_51, gpd_52, \
                         gpd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * pc_y[k] * gpd_50[k];

        t_85[k] = f_6 * fpd_53[k]
                  + f_4 * pc_x[k] * gpd_53[k];

        t_86[k] = f_2 * gpp0_25[k]
                  - f_3 * gpp1_25[k]
                  + f_4 * pc_y[k] * gpd_51[k];

        t_87[k] = f_9 * gpp0_26[k]
                  - f_10 * gpp1_26[k]
                  + f_4 * pc_y[k] * gpd_52[k];

        t_88[k] = f_4 * pc_y[k] * gpd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, pa_x, pa_y, pc_x, pc_y, dpf0_0, dpf0_89, dpf1_0, dpf1_89, \
                         fpf0_30, fpf0_89, fpf1_30, fpf1_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_7 * dpf0_89[k]
                  - f_8 * dpf1_89[k]
                  + pa_x[k] * fpf0_89[k]
                  - f_5 * pc_x[k] * fpf1_89[k];

        t_90[k] = f_12 * dpf0_0[k]
                  - f_13 * dpf1_0[k]
                  + pa_y[k] * fpf0_30[k]
                  - f_5 * pc_y[k] * fpf1_30[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, fpd_18, fpd_57, gsd_21, \
                         gpd_54, gpd_55, gpd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_11 * fpd_18[k]
                  + f_4 * pc_y[k] * gpd_54[k];

        t_92[k] = f_4 * pc_z[k] * gpd_54[k];

        t_93[k] = f_11 * fpd_57[k]
                  + f_1 * gsd_21[k]
                  + f_4 * pc_x[k] * gpd_57[k];

        t_94[k] = f_4 * pc_z[k] * gpd_55[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, fpd_21, fpd_23, fpd_59, \
                         gsd_23, gpp0_28, gpp1_28, gpd_57, gpd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * fpd_59[k]
                  + f_1 * gsd_23[k]
                  + f_4 * pc_x[k] * gpd_59[k];

        t_96[k] = f_11 * fpd_21[k]
                  + f_2 * gpp0_28[k]
                  - f_3 * gpp1_28[k]
                  + f_4 * pc_y[k] * gpd_57[k];

        t_97[k] = f_4 * pc_z[k] * gpd_57[k];

        t_98[k] = f_11 * fpd_23[k]
                  + f_4 * pc_y[k] * gpd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, fpd_24, fpd_60, gsd_18, \
                         gpp0_29, gpp0_30, gpp1_29, gpp1_30, gpd_59, \
                         gpd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_2 * gpp0_29[k]
                  - f_3 * gpp1_29[k]
                  + f_4 * pc_z[k] * gpd_59[k];

        t_100[k] = f_11 * fpd_60[k]
                   + f_2 * gpp0_30[k]
                   - f_3 * gpp1_30[k]
                   + f_4 * pc_x[k] * gpd_60[k];

        t_101[k] = f_11 * fpd_24[k]
                   + f_1 * gsd_18[k]
                   + f_4 * pc_y[k] * gpd_60[k];

        t_102[k] = f_4 * pc_z[k] * gpd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_x, pc_x, pc_z, dpf0_106, dpf1_106, \
                         fpf0_106, fpd_63, fpd_65, fpf1_106, gpd_61, gpd_63, \
                         gpd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_11 * fpd_63[k]
                   + f_4 * pc_x[k] * gpd_63[k];

        t_104[k] = f_4 * pc_z[k] * gpd_61[k];

        t_105[k] = f_11 * fpd_65[k]
                   + f_4 * pc_x[k] * gpd_65[k];

        t_106[k] = f_12 * dpf0_106[k]
                   - f_13 * dpf1_106[k]
                   + pa_x[k] * fpf0_106[k]
                   - f_5 * pc_x[k] * fpf1_106[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, fpd_29, gsf0_30, \
                         gsd_23, gsf1_30, gpp0_32, gpp1_32, gpd_63, \
                         gpd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_4 * pc_z[k] * gpd_63[k];

        t_108[k] = f_11 * fpd_29[k]
                   + f_1 * gsd_23[k]
                   + f_4 * pc_y[k] * gpd_65[k];

        t_109[k] = f_2 * gpp0_32[k]
                   - f_3 * gpp1_32[k]
                   + f_4 * pc_z[k] * gpd_65[k];

        t_110[k] = pb_z[k] * gsf0_30[k]
                   - f_5 * pc_z[k] * gsf1_30[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, fpd_30, fpd_69, gsd_18, \
                         gsd_19, gpd_66, gpd_67, gpd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_11 * fpd_30[k]
                   + f_4 * pc_y[k] * gpd_66[k];

        t_112[k] = f_1 * gsd_18[k]
                   + f_4 * pc_z[k] * gpd_66[k];

        t_113[k] = f_11 * fpd_69[k]
                   + f_4 * pc_x[k] * gpd_69[k];

        t_114[k] = f_1 * gsd_19[k]
                   + f_4 * pc_z[k] * gpd_67[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, fpd_35, fpd_71, \
                         gsf0_36, gsd_21, gsf1_36, gpd_69, gpd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_11 * fpd_71[k]
                   + f_4 * pc_x[k] * gpd_71[k];

        t_116[k] = pb_z[k] * gsf0_36[k]
                   - f_5 * pc_z[k] * gsf1_36[k];

        t_117[k] = f_1 * gsd_21[k]
                   + f_4 * pc_z[k] * gpd_69[k];

        t_118[k] = f_11 * fpd_35[k]
                   + f_4 * pc_y[k] * gpd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_y, pb_z, pc_y, pc_z, fpf0_60, fpd_18, \
                         fpd_36, fpf1_60, gsf0_39, gsd_23, gsf1_39, \
                         gpd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_z[k] * gsf0_39[k]
                   + f_6 * gsd_23[k]
                   - f_5 * pc_z[k] * gsf1_39[k];

        t_120[k] = pa_y[k] * fpf0_60[k]
                   - f_5 * pc_y[k] * fpf1_60[k];

        t_121[k] = f_1 * fpd_36[k]
                   + f_4 * pc_y[k] * gpd_72[k];

        t_122[k] = f_1 * fpd_18[k]
                   + f_4 * pc_z[k] * gpd_72[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pa_y, pa_z, pc_x, pc_y, pc_z, fpf0_33, fpf0_65, \
                         fpd_76, fpf1_33, fpf1_65, gsd_28, gpd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = pa_z[k] * fpf0_33[k]
                   - f_5 * pc_z[k] * fpf1_33[k];

        t_124[k] = f_11 * fpd_76[k]
                   + f_1 * gsd_28[k]
                   + f_4 * pc_x[k] * gpd_76[k];

        t_125[k] = pa_y[k] * fpf0_65[k]
                   - f_5 * pc_y[k] * fpf1_65[k];
    }
}

static auto
compute_prim_gpf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t fpf0, const size_t fpd,
                                                          const size_t fpf1, const size_t gsf0,
                                                          const size_t gsd, const size_t gsf1,
                                                          const size_t gpp0, const size_t gpp1,
                                                          const size_t gpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);
    const auto f_11 = 1.0 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpf0_0 = buffer.data(dpf0 + 0);
    const auto *dpf0_179 = buffer.data(dpf0 + 179);

    const auto *dpf1_0 = buffer.data(dpf1 + 0);
    const auto *dpf1_179 = buffer.data(dpf1 + 179);

    const auto *fpf0_36 = buffer.data(fpf0 + 36);
    const auto *fpf0_41 = buffer.data(fpf0 + 41);
    const auto *fpf0_46 = buffer.data(fpf0 + 46);
    const auto *fpf0_60 = buffer.data(fpf0 + 60);
    const auto *fpf0_69 = buffer.data(fpf0 + 69);
    const auto *fpf0_80 = buffer.data(fpf0 + 80);
    const auto *fpf0_82 = buffer.data(fpf0 + 82);
    const auto *fpf0_89 = buffer.data(fpf0 + 89);
    const auto *fpf0_90 = buffer.data(fpf0 + 90);
    const auto *fpf0_93 = buffer.data(fpf0 + 93);
    const auto *fpf0_100 = buffer.data(fpf0 + 100);
    const auto *fpf0_101 = buffer.data(fpf0 + 101);
    const auto *fpf0_150 = buffer.data(fpf0 + 150);
    const auto *fpf0_155 = buffer.data(fpf0 + 155);
    const auto *fpf0_179 = buffer.data(fpf0 + 179);
    const auto *fpf0_190 = buffer.data(fpf0 + 190);
    const auto *fpf0_196 = buffer.data(fpf0 + 196);
    const auto *fpf0_198 = buffer.data(fpf0 + 198);
    const auto *fpf0_199 = buffer.data(fpf0 + 199);
    const auto *fpf0_206 = buffer.data(fpf0 + 206);
    const auto *fpf0_209 = buffer.data(fpf0 + 209);
    const auto *fpf0_226 = buffer.data(fpf0 + 226);
    const auto *fpf0_228 = buffer.data(fpf0 + 228);
    const auto *fpf0_229 = buffer.data(fpf0 + 229);
    const auto *fpf0_236 = buffer.data(fpf0 + 236);
    const auto *fpf0_237 = buffer.data(fpf0 + 237);
    const auto *fpf0_239 = buffer.data(fpf0 + 239);

    const auto *fpd_21 = buffer.data(fpd + 21);
    const auto *fpd_24 = buffer.data(fpd + 24);
    const auto *fpd_27 = buffer.data(fpd + 27);
    const auto *fpd_29 = buffer.data(fpd + 29);
    const auto *fpd_33 = buffer.data(fpd + 33);
    const auto *fpd_36 = buffer.data(fpd + 36);
    const auto *fpd_41 = buffer.data(fpd + 41);
    const auto *fpd_42 = buffer.data(fpd + 42);
    const auto *fpd_47 = buffer.data(fpd + 47);
    const auto *fpd_48 = buffer.data(fpd + 48);
    const auto *fpd_51 = buffer.data(fpd + 51);
    const auto *fpd_53 = buffer.data(fpd + 53);
    const auto *fpd_54 = buffer.data(fpd + 54);
    const auto *fpd_57 = buffer.data(fpd + 57);
    const auto *fpd_59 = buffer.data(fpd + 59);
    const auto *fpd_60 = buffer.data(fpd + 60);
    const auto *fpd_63 = buffer.data(fpd + 63);
    const auto *fpd_66 = buffer.data(fpd + 66);
    const auto *fpd_71 = buffer.data(fpd + 71);
    const auto *fpd_72 = buffer.data(fpd + 72);
    const auto *fpd_75 = buffer.data(fpd + 75);
    const auto *fpd_77 = buffer.data(fpd + 77);
    const auto *fpd_78 = buffer.data(fpd + 78);
    const auto *fpd_81 = buffer.data(fpd + 81);
    const auto *fpd_82 = buffer.data(fpd + 82);
    const auto *fpd_83 = buffer.data(fpd + 83);
    const auto *fpd_84 = buffer.data(fpd + 84);
    const auto *fpd_87 = buffer.data(fpd + 87);
    const auto *fpd_88 = buffer.data(fpd + 88);
    const auto *fpd_89 = buffer.data(fpd + 89);
    const auto *fpd_90 = buffer.data(fpd + 90);
    const auto *fpd_93 = buffer.data(fpd + 93);
    const auto *fpd_95 = buffer.data(fpd + 95);
    const auto *fpd_99 = buffer.data(fpd + 99);
    const auto *fpd_101 = buffer.data(fpd + 101);
    const auto *fpd_102 = buffer.data(fpd + 102);
    const auto *fpd_105 = buffer.data(fpd + 105);
    const auto *fpd_107 = buffer.data(fpd + 107);
    const auto *fpd_108 = buffer.data(fpd + 108);
    const auto *fpd_111 = buffer.data(fpd + 111);
    const auto *fpd_113 = buffer.data(fpd + 113);
    const auto *fpd_114 = buffer.data(fpd + 114);
    const auto *fpd_117 = buffer.data(fpd + 117);
    const auto *fpd_119 = buffer.data(fpd + 119);
    const auto *fpd_123 = buffer.data(fpd + 123);
    const auto *fpd_125 = buffer.data(fpd + 125);
    const auto *fpd_130 = buffer.data(fpd + 130);
    const auto *fpd_131 = buffer.data(fpd + 131);
    const auto *fpd_135 = buffer.data(fpd + 135);
    const auto *fpd_136 = buffer.data(fpd + 136);
    const auto *fpd_137 = buffer.data(fpd + 137);
    const auto *fpd_138 = buffer.data(fpd + 138);
    const auto *fpd_141 = buffer.data(fpd + 141);
    const auto *fpd_142 = buffer.data(fpd + 142);
    const auto *fpd_143 = buffer.data(fpd + 143);
    const auto *fpd_147 = buffer.data(fpd + 147);
    const auto *fpd_148 = buffer.data(fpd + 148);

    const auto *fpf1_36 = buffer.data(fpf1 + 36);
    const auto *fpf1_41 = buffer.data(fpf1 + 41);
    const auto *fpf1_46 = buffer.data(fpf1 + 46);
    const auto *fpf1_60 = buffer.data(fpf1 + 60);
    const auto *fpf1_69 = buffer.data(fpf1 + 69);
    const auto *fpf1_80 = buffer.data(fpf1 + 80);
    const auto *fpf1_82 = buffer.data(fpf1 + 82);
    const auto *fpf1_89 = buffer.data(fpf1 + 89);
    const auto *fpf1_90 = buffer.data(fpf1 + 90);
    const auto *fpf1_93 = buffer.data(fpf1 + 93);
    const auto *fpf1_100 = buffer.data(fpf1 + 100);
    const auto *fpf1_101 = buffer.data(fpf1 + 101);
    const auto *fpf1_150 = buffer.data(fpf1 + 150);
    const auto *fpf1_155 = buffer.data(fpf1 + 155);
    const auto *fpf1_179 = buffer.data(fpf1 + 179);
    const auto *fpf1_190 = buffer.data(fpf1 + 190);
    const auto *fpf1_196 = buffer.data(fpf1 + 196);
    const auto *fpf1_198 = buffer.data(fpf1 + 198);
    const auto *fpf1_199 = buffer.data(fpf1 + 199);
    const auto *fpf1_206 = buffer.data(fpf1 + 206);
    const auto *fpf1_209 = buffer.data(fpf1 + 209);
    const auto *fpf1_226 = buffer.data(fpf1 + 226);
    const auto *fpf1_228 = buffer.data(fpf1 + 228);
    const auto *fpf1_229 = buffer.data(fpf1 + 229);
    const auto *fpf1_236 = buffer.data(fpf1 + 236);
    const auto *fpf1_237 = buffer.data(fpf1 + 237);
    const auto *fpf1_239 = buffer.data(fpf1 + 239);

    const auto *gsf0_50 = buffer.data(gsf0 + 50);
    const auto *gsf0_56 = buffer.data(gsf0 + 56);
    const auto *gsf0_57 = buffer.data(gsf0 + 57);
    const auto *gsf0_59 = buffer.data(gsf0 + 59);
    const auto *gsf0_60 = buffer.data(gsf0 + 60);

    const auto *gsd_27 = buffer.data(gsd + 27);
    const auto *gsd_29 = buffer.data(gsd + 29);
    const auto *gsd_30 = buffer.data(gsd + 30);
    const auto *gsd_32 = buffer.data(gsd + 32);
    const auto *gsd_33 = buffer.data(gsd + 33);
    const auto *gsd_34 = buffer.data(gsd + 34);
    const auto *gsd_35 = buffer.data(gsd + 35);
    const auto *gsd_36 = buffer.data(gsd + 36);
    const auto *gsd_37 = buffer.data(gsd + 37);
    const auto *gsd_39 = buffer.data(gsd + 39);
    const auto *gsd_41 = buffer.data(gsd + 41);
    const auto *gsd_42 = buffer.data(gsd + 42);
    const auto *gsd_46 = buffer.data(gsd + 46);
    const auto *gsd_47 = buffer.data(gsd + 47);
    const auto *gsd_51 = buffer.data(gsd + 51);
    const auto *gsd_52 = buffer.data(gsd + 52);

    const auto *gsf1_50 = buffer.data(gsf1 + 50);
    const auto *gsf1_56 = buffer.data(gsf1 + 56);
    const auto *gsf1_57 = buffer.data(gsf1 + 57);
    const auto *gsf1_59 = buffer.data(gsf1 + 59);
    const auto *gsf1_60 = buffer.data(gsf1 + 60);

    const auto *gpp0_39 = buffer.data(gpp0 + 39);
    const auto *gpp0_41 = buffer.data(gpp0 + 41);
    const auto *gpp0_43 = buffer.data(gpp0 + 43);
    const auto *gpp0_46 = buffer.data(gpp0 + 46);
    const auto *gpp0_47 = buffer.data(gpp0 + 47);
    const auto *gpp0_51 = buffer.data(gpp0 + 51);
    const auto *gpp0_52 = buffer.data(gpp0 + 52);
    const auto *gpp0_53 = buffer.data(gpp0 + 53);
    const auto *gpp0_54 = buffer.data(gpp0 + 54);
    const auto *gpp0_55 = buffer.data(gpp0 + 55);
    const auto *gpp0_56 = buffer.data(gpp0 + 56);
    const auto *gpp0_64 = buffer.data(gpp0 + 64);
    const auto *gpp0_65 = buffer.data(gpp0 + 65);
    const auto *gpp0_69 = buffer.data(gpp0 + 69);
    const auto *gpp0_73 = buffer.data(gpp0 + 73);

    const auto *gpp1_39 = buffer.data(gpp1 + 39);
    const auto *gpp1_41 = buffer.data(gpp1 + 41);
    const auto *gpp1_43 = buffer.data(gpp1 + 43);
    const auto *gpp1_46 = buffer.data(gpp1 + 46);
    const auto *gpp1_47 = buffer.data(gpp1 + 47);
    const auto *gpp1_51 = buffer.data(gpp1 + 51);
    const auto *gpp1_52 = buffer.data(gpp1 + 52);
    const auto *gpp1_53 = buffer.data(gpp1 + 53);
    const auto *gpp1_54 = buffer.data(gpp1 + 54);
    const auto *gpp1_55 = buffer.data(gpp1 + 55);
    const auto *gpp1_56 = buffer.data(gpp1 + 56);
    const auto *gpp1_64 = buffer.data(gpp1 + 64);
    const auto *gpp1_65 = buffer.data(gpp1 + 65);
    const auto *gpp1_69 = buffer.data(gpp1 + 69);
    const auto *gpp1_73 = buffer.data(gpp1 + 73);

    const auto *gpd_75 = buffer.data(gpd + 75);
    const auto *gpd_77 = buffer.data(gpd + 77);
    const auto *gpd_78 = buffer.data(gpd + 78);
    const auto *gpd_81 = buffer.data(gpd + 81);
    const auto *gpd_82 = buffer.data(gpd + 82);
    const auto *gpd_83 = buffer.data(gpd + 83);
    const auto *gpd_84 = buffer.data(gpd + 84);
    const auto *gpd_87 = buffer.data(gpd + 87);
    const auto *gpd_88 = buffer.data(gpd + 88);
    const auto *gpd_89 = buffer.data(gpd + 89);
    const auto *gpd_90 = buffer.data(gpd + 90);
    const auto *gpd_92 = buffer.data(gpd + 92);
    const auto *gpd_93 = buffer.data(gpd + 93);
    const auto *gpd_94 = buffer.data(gpd + 94);
    const auto *gpd_95 = buffer.data(gpd + 95);
    const auto *gpd_96 = buffer.data(gpd + 96);
    const auto *gpd_98 = buffer.data(gpd + 98);
    const auto *gpd_99 = buffer.data(gpd + 99);
    const auto *gpd_101 = buffer.data(gpd + 101);
    const auto *gpd_102 = buffer.data(gpd + 102);
    const auto *gpd_104 = buffer.data(gpd + 104);
    const auto *gpd_105 = buffer.data(gpd + 105);
    const auto *gpd_106 = buffer.data(gpd + 106);
    const auto *gpd_107 = buffer.data(gpd + 107);
    const auto *gpd_108 = buffer.data(gpd + 108);
    const auto *gpd_109 = buffer.data(gpd + 109);
    const auto *gpd_111 = buffer.data(gpd + 111);
    const auto *gpd_113 = buffer.data(gpd + 113);
    const auto *gpd_114 = buffer.data(gpd + 114);
    const auto *gpd_115 = buffer.data(gpd + 115);
    const auto *gpd_117 = buffer.data(gpd + 117);
    const auto *gpd_119 = buffer.data(gpd + 119);
    const auto *gpd_120 = buffer.data(gpd + 120);
    const auto *gpd_121 = buffer.data(gpd + 121);
    const auto *gpd_123 = buffer.data(gpd + 123);
    const auto *gpd_125 = buffer.data(gpd + 125);
    const auto *gpd_126 = buffer.data(gpd + 126);
    const auto *gpd_129 = buffer.data(gpd + 129);
    const auto *gpd_130 = buffer.data(gpd + 130);
    const auto *gpd_131 = buffer.data(gpd + 131);
    const auto *gpd_132 = buffer.data(gpd + 132);
    const auto *gpd_135 = buffer.data(gpd + 135);
    const auto *gpd_136 = buffer.data(gpd + 136);
    const auto *gpd_137 = buffer.data(gpd + 137);
    const auto *gpd_138 = buffer.data(gpd + 138);
    const auto *gpd_141 = buffer.data(gpd + 141);
    const auto *gpd_142 = buffer.data(gpd + 142);
    const auto *gpd_143 = buffer.data(gpd + 143);
    const auto *gpd_144 = buffer.data(gpd + 144);
    const auto *gpd_147 = buffer.data(gpd + 147);
    const auto *gpd_148 = buffer.data(gpd + 148);
    const auto *gpd_149 = buffer.data(gpd + 149);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_y, pa_z, pc_y, pc_z, fpf0_36, fpf0_69, \
                         fpd_21, fpd_41, fpf1_36, fpf1_69, gpd_75, \
                         gpd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * fpf0_36[k]
                   - f_5 * pc_z[k] * fpf1_36[k];

        t_127[k] = f_1 * fpd_21[k]
                   + f_4 * pc_z[k] * gpd_75[k];

        t_128[k] = f_1 * fpd_41[k]
                   + f_4 * pc_y[k] * gpd_77[k];

        t_129[k] = pa_y[k] * fpf0_69[k]
                   - f_5 * pc_y[k] * fpf1_69[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_z, pc_x, pc_z, fpf0_41, fpd_24, \
                         fpd_78, fpd_81, fpf1_41, gpp0_39, gpp1_39, gpd_78, \
                         gpd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * fpd_78[k]
                   + f_2 * gpp0_39[k]
                   - f_3 * gpp1_39[k]
                   + f_4 * pc_x[k] * gpd_78[k];

        t_131[k] = pa_z[k] * fpf0_41[k]
                   - f_5 * pc_z[k] * fpf1_41[k];

        t_132[k] = f_1 * fpd_24[k]
                   + f_4 * pc_z[k] * gpd_78[k];

        t_133[k] = f_11 * fpd_81[k]
                   + f_4 * pc_x[k] * gpd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_z, pc_x, pc_z, fpf0_46, fpd_27, \
                         fpd_82, fpd_83, fpf1_46, gpd_81, gpd_82, \
                         gpd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * fpd_82[k]
                   + f_4 * pc_x[k] * gpd_82[k];

        t_135[k] = f_11 * fpd_83[k]
                   + f_4 * pc_x[k] * gpd_83[k];

        t_136[k] = pa_z[k] * fpf0_46[k]
                   - f_5 * pc_z[k] * fpf1_46[k];

        t_137[k] = f_1 * fpd_27[k]
                   + f_4 * pc_z[k] * gpd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pa_y, pc_y, pc_z, fpf0_80, fpd_29, fpd_47, \
                         fpf1_80, gsd_29, gpp0_41, gpp1_41, gpd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_1 * fpd_47[k]
                   + f_1 * gsd_29[k]
                   + f_4 * pc_y[k] * gpd_83[k];

        t_139[k] = f_1 * fpd_29[k]
                   + f_2 * gpp0_41[k]
                   - f_3 * gpp1_41[k]
                   + f_4 * pc_z[k] * gpd_83[k];

        t_140[k] = pa_y[k] * fpf0_80[k]
                   - f_5 * pc_y[k] * fpf1_80[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_x, pc_y, fpf0_82, fpd_48, \
                         fpd_87, fpd_88, fpf1_82, gpd_84, gpd_87, \
                         gpd_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_1 * fpd_48[k]
                   + f_4 * pc_y[k] * gpd_84[k];

        t_142[k] = pa_y[k] * fpf0_82[k]
                   - f_5 * pc_y[k] * fpf1_82[k];

        t_143[k] = f_11 * fpd_87[k]
                   + f_4 * pc_x[k] * gpd_87[k];

        t_144[k] = f_11 * fpd_88[k]
                   + f_4 * pc_x[k] * gpd_88[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, pc_z, fpd_33, fpd_51, fpd_53, \
                         fpd_89, gsd_27, gpp0_43, gpp1_43, gpd_87, \
                         gpd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_11 * fpd_89[k]
                   + f_4 * pc_x[k] * gpd_89[k];

        t_146[k] = f_1 * fpd_51[k]
                   + f_2 * gpp0_43[k]
                   - f_3 * gpp1_43[k]
                   + f_4 * pc_y[k] * gpd_87[k];

        t_147[k] = f_1 * fpd_33[k]
                   + f_1 * gsd_27[k]
                   + f_4 * pc_z[k] * gpd_87[k];

        t_148[k] = f_1 * fpd_53[k]
                   + f_4 * pc_y[k] * gpd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, dpf0_0, dpf1_0, \
                         fpf0_60, fpf0_89, fpd_36, fpf1_60, fpf1_89, \
                         gpd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * fpf0_89[k]
                   - f_5 * pc_y[k] * fpf1_89[k];

        t_150[k] = f_12 * dpf0_0[k]
                   - f_13 * dpf1_0[k]
                   + pa_z[k] * fpf0_60[k]
                   - f_5 * pc_z[k] * fpf1_60[k];

        t_151[k] = f_4 * pc_y[k] * gpd_90[k];

        t_152[k] = f_11 * fpd_36[k]
                   + f_4 * pc_z[k] * gpd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pc_x, pc_y, fpd_93, fpd_95, gsd_33, \
                         gsd_35, gpp0_46, gpp1_46, gpd_92, gpd_93, \
                         gpd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_11 * fpd_93[k]
                   + f_1 * gsd_33[k]
                   + f_4 * pc_x[k] * gpd_93[k];

        t_154[k] = f_4 * pc_y[k] * gpd_92[k];

        t_155[k] = f_11 * fpd_95[k]
                   + f_1 * gsd_35[k]
                   + f_4 * pc_x[k] * gpd_95[k];

        t_156[k] = f_2 * gpp0_46[k]
                   - f_3 * gpp1_46[k]
                   + f_4 * pc_y[k] * gpd_93[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pc_y, pc_z, fpd_41, gsf0_50, \
                         gsf1_50, gpp0_47, gpp1_47, gpd_94, gpd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_9 * gpp0_47[k]
                   - f_10 * gpp1_47[k]
                   + f_4 * pc_y[k] * gpd_94[k];

        t_158[k] = f_4 * pc_y[k] * gpd_95[k];

        t_159[k] = f_11 * fpd_41[k]
                   + f_2 * gpp0_47[k]
                   - f_3 * gpp1_47[k]
                   + f_4 * pc_z[k] * gpd_95[k];

        t_160[k] = pb_y[k] * gsf0_50[k]
                   - f_5 * pc_y[k] * gsf1_50[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pc_x, pc_y, pc_z, fpd_42, fpd_99, gsd_30, \
                         gsd_32, gpd_96, gpd_98, gpd_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * gsd_30[k]
                   + f_4 * pc_y[k] * gpd_96[k];

        t_162[k] = f_11 * fpd_42[k]
                   + f_4 * pc_z[k] * gpd_96[k];

        t_163[k] = f_11 * fpd_99[k]
                   + f_4 * pc_x[k] * gpd_99[k];

        t_164[k] = f_1 * gsd_32[k]
                   + f_4 * pc_y[k] * gpd_98[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_y, pc_x, pc_y, fpd_101, gsf0_56, \
                         gsf0_57, gsd_33, gsd_34, gsd_35, gsf1_56, gsf1_57, \
                         gpd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_11 * fpd_101[k]
                   + f_4 * pc_x[k] * gpd_101[k];

        t_166[k] = pb_y[k] * gsf0_56[k]
                   + f_6 * gsd_33[k]
                   - f_5 * pc_y[k] * gsf1_56[k];

        t_167[k] = pb_y[k] * gsf0_57[k]
                   + f_11 * gsd_34[k]
                   - f_5 * pc_y[k] * gsf1_57[k];

        t_168[k] = f_1 * gsd_35[k]
                   + f_4 * pc_y[k] * gpd_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pb_y, pc_x, pc_y, pc_z, fpd_48, fpd_102, \
                         gsf0_59, gsd_30, gsf1_59, gpp0_51, gpp1_51, \
                         gpd_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pb_y[k] * gsf0_59[k]
                   - f_5 * pc_y[k] * gsf1_59[k];

        t_170[k] = f_11 * fpd_102[k]
                   + f_2 * gpp0_51[k]
                   - f_3 * gpp1_51[k]
                   + f_4 * pc_x[k] * gpd_102[k];

        t_171[k] = f_4 * pc_y[k] * gpd_102[k];

        t_172[k] = f_11 * fpd_48[k]
                   + f_1 * gsd_30[k]
                   + f_4 * pc_z[k] * gpd_102[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pc_x, pc_y, fpd_105, fpd_107, gpp0_52, \
                         gpp1_52, gpd_104, gpd_105, gpd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_11 * fpd_105[k]
                   + f_4 * pc_x[k] * gpd_105[k];

        t_174[k] = f_4 * pc_y[k] * gpd_104[k];

        t_175[k] = f_11 * fpd_107[k]
                   + f_4 * pc_x[k] * gpd_107[k];

        t_176[k] = f_2 * gpp0_52[k]
                   - f_3 * gpp1_52[k]
                   + f_4 * pc_y[k] * gpd_105[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_x, pc_x, pc_y, dpf0_179, dpf1_179, fpf0_179, \
                         fpf1_179, gpp0_53, gpp1_53, gpd_106, gpd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_9 * gpp0_53[k]
                   - f_10 * gpp1_53[k]
                   + f_4 * pc_y[k] * gpd_106[k];

        t_178[k] = f_4 * pc_y[k] * gpd_107[k];

        t_179[k] = f_12 * dpf0_179[k]
                   - f_13 * dpf1_179[k]
                   + pa_x[k] * fpf0_179[k]
                   - f_5 * pc_x[k] * fpf1_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, fpd_54, fpd_108, \
                         fpd_111, gsd_36, gsd_39, gpp0_54, gpp1_54, gpd_108, \
                         gpd_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * fpd_108[k]
                   + f_1 * gsd_36[k]
                   + f_2 * gpp0_54[k]
                   - f_3 * gpp1_54[k]
                   + f_4 * pc_x[k] * gpd_108[k];

        t_181[k] = f_6 * fpd_54[k]
                   + f_4 * pc_y[k] * gpd_108[k];

        t_182[k] = f_4 * pc_z[k] * gpd_108[k];

        t_183[k] = f_1 * fpd_111[k]
                   + f_1 * gsd_39[k]
                   + f_4 * pc_x[k] * gpd_111[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pc_x, pc_y, pc_z, fpd_57, fpd_113, \
                         gsd_41, gpp0_55, gpp1_55, gpd_109, gpd_111, \
                         gpd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_4 * pc_z[k] * gpd_109[k];

        t_185[k] = f_1 * fpd_113[k]
                   + f_1 * gsd_41[k]
                   + f_4 * pc_x[k] * gpd_113[k];

        t_186[k] = f_6 * fpd_57[k]
                   + f_2 * gpp0_55[k]
                   - f_3 * gpp1_55[k]
                   + f_4 * pc_y[k] * gpd_111[k];

        t_187[k] = f_4 * pc_z[k] * gpd_111[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_x, pc_x, pc_y, pc_z, fpf0_190, fpd_59, \
                         fpd_114, fpf1_190, gpp0_56, gpp1_56, gpd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_6 * fpd_59[k]
                   + f_4 * pc_y[k] * gpd_113[k];

        t_189[k] = f_2 * gpp0_56[k]
                   - f_3 * gpp1_56[k]
                   + f_4 * pc_z[k] * gpd_113[k];

        t_190[k] = pa_x[k] * fpf0_190[k]
                   + f_6 * fpd_114[k]
                   - f_5 * pc_x[k] * fpf1_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pc_x, pc_y, pc_z, fpd_60, fpd_117, \
                         fpd_119, gsd_36, gpd_114, gpd_115, gpd_117, \
                         gpd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_6 * fpd_60[k]
                   + f_1 * gsd_36[k]
                   + f_4 * pc_y[k] * gpd_114[k];

        t_192[k] = f_4 * pc_z[k] * gpd_114[k];

        t_193[k] = f_1 * fpd_117[k]
                   + f_4 * pc_x[k] * gpd_117[k];

        t_194[k] = f_4 * pc_z[k] * gpd_115[k];

        t_195[k] = f_1 * fpd_119[k]
                   + f_4 * pc_x[k] * gpd_119[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_x, pc_x, pc_z, fpf0_196, fpf0_198, \
                         fpf0_199, fpf1_196, fpf1_198, fpf1_199, \
                         gpd_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * fpf0_196[k]
                   - f_5 * pc_x[k] * fpf1_196[k];

        t_197[k] = f_4 * pc_z[k] * gpd_117[k];

        t_198[k] = pa_x[k] * fpf0_198[k]
                   - f_5 * pc_x[k] * fpf1_198[k];

        t_199[k] = pa_x[k] * fpf0_199[k]
                   - f_5 * pc_x[k] * fpf1_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pb_z, pc_x, pc_y, pc_z, fpd_66, fpd_123, \
                         gsf0_60, gsd_36, gsf1_60, gpd_120, gpd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_z[k] * gsf0_60[k]
                   - f_5 * pc_z[k] * gsf1_60[k];

        t_201[k] = f_6 * fpd_66[k]
                   + f_4 * pc_y[k] * gpd_120[k];

        t_202[k] = f_1 * gsd_36[k]
                   + f_4 * pc_z[k] * gpd_120[k];

        t_203[k] = f_1 * fpd_123[k]
                   + f_4 * pc_x[k] * gpd_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, pc_x, pc_z, fpf0_206, fpd_125, \
                         fpf1_206, gsd_37, gsd_39, gpd_121, gpd_123, \
                         gpd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * gsd_37[k]
                   + f_4 * pc_z[k] * gpd_121[k];

        t_205[k] = f_1 * fpd_125[k]
                   + f_4 * pc_x[k] * gpd_125[k];

        t_206[k] = pa_x[k] * fpf0_206[k]
                   - f_5 * pc_x[k] * fpf1_206[k];

        t_207[k] = f_1 * gsd_39[k]
                   + f_4 * pc_z[k] * gpd_123[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_x, pa_z, pc_x, pc_y, pc_z, fpf0_90, fpf0_209, \
                         fpd_71, fpf1_90, fpf1_209, gpd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_6 * fpd_71[k]
                   + f_4 * pc_y[k] * gpd_125[k];

        t_209[k] = pa_x[k] * fpf0_209[k]
                   - f_5 * pc_x[k] * fpf1_209[k];

        t_210[k] = pa_z[k] * fpf0_90[k]
                   - f_5 * pc_z[k] * fpf1_90[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pa_z, pc_x, pc_y, pc_z, fpf0_93, fpd_54, \
                         fpd_72, fpd_130, fpf1_93, gsd_46, gpd_126, \
                         gpd_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_11 * fpd_72[k]
                   + f_4 * pc_y[k] * gpd_126[k];

        t_212[k] = f_1 * fpd_54[k]
                   + f_4 * pc_z[k] * gpd_126[k];

        t_213[k] = pa_z[k] * fpf0_93[k]
                   - f_5 * pc_z[k] * fpf1_93[k];

        t_214[k] = f_1 * fpd_130[k]
                   + f_1 * gsd_46[k]
                   + f_4 * pc_x[k] * gpd_130[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, pc_z, fpd_57, fpd_75, fpd_77, \
                         fpd_131, gsd_47, gpp0_64, gpp1_64, gpd_129, \
                         gpd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_1 * fpd_131[k]
                   + f_1 * gsd_47[k]
                   + f_4 * pc_x[k] * gpd_131[k];

        t_216[k] = f_11 * fpd_75[k]
                   + f_2 * gpp0_64[k]
                   - f_3 * gpp1_64[k]
                   + f_4 * pc_y[k] * gpd_129[k];

        t_217[k] = f_1 * fpd_57[k]
                   + f_4 * pc_z[k] * gpd_129[k];

        t_218[k] = f_11 * fpd_77[k]
                   + f_4 * pc_y[k] * gpd_131[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pa_z, pc_z, fpf0_100, fpf0_101, fpd_59, \
                         fpd_60, fpf1_100, fpf1_101, gpp0_65, gpp1_65, gpd_131, \
                         gpd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_1 * fpd_59[k]
                   + f_2 * gpp0_65[k]
                   - f_3 * gpp1_65[k]
                   + f_4 * pc_z[k] * gpd_131[k];

        t_220[k] = pa_z[k] * fpf0_100[k]
                   - f_5 * pc_z[k] * fpf1_100[k];

        t_221[k] = pa_z[k] * fpf0_101[k]
                   - f_5 * pc_z[k] * fpf1_101[k];

        t_222[k] = f_1 * fpd_60[k]
                   + f_4 * pc_z[k] * gpd_132[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pc_x, fpf0_226, fpd_135, fpd_136, \
                         fpd_137, fpf1_226, gpd_135, gpd_136, gpd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_1 * fpd_135[k]
                   + f_4 * pc_x[k] * gpd_135[k];

        t_224[k] = f_1 * fpd_136[k]
                   + f_4 * pc_x[k] * gpd_136[k];

        t_225[k] = f_1 * fpd_137[k]
                   + f_4 * pc_x[k] * gpd_137[k];

        t_226[k] = pa_x[k] * fpf0_226[k]
                   - f_5 * pc_x[k] * fpf1_226[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pa_x, pc_x, pc_z, fpf0_228, fpf0_229, fpd_63, \
                         fpf1_228, fpf1_229, gpd_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_1 * fpd_63[k]
                   + f_4 * pc_z[k] * gpd_135[k];

        t_228[k] = pa_x[k] * fpf0_228[k]
                   - f_5 * pc_x[k] * fpf1_228[k];

        t_229[k] = pa_x[k] * fpf0_229[k]
                   - f_5 * pc_x[k] * fpf1_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, fpd_66, fpd_84, \
                         fpd_138, fpd_141, gsd_42, gpp0_69, gpp1_69, gpd_138, \
                         gpd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * fpd_138[k]
                   + f_2 * gpp0_69[k]
                   - f_3 * gpp1_69[k]
                   + f_4 * pc_x[k] * gpd_138[k];

        t_231[k] = f_11 * fpd_84[k]
                   + f_4 * pc_y[k] * gpd_138[k];

        t_232[k] = f_1 * fpd_66[k]
                   + f_1 * gsd_42[k]
                   + f_4 * pc_z[k] * gpd_138[k];

        t_233[k] = f_1 * fpd_141[k]
                   + f_4 * pc_x[k] * gpd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_x, pc_x, fpf0_236, fpf0_237, fpd_142, \
                         fpd_143, fpf1_236, fpf1_237, gpd_142, \
                         gpd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_1 * fpd_142[k]
                   + f_4 * pc_x[k] * gpd_142[k];

        t_235[k] = f_1 * fpd_143[k]
                   + f_4 * pc_x[k] * gpd_143[k];

        t_236[k] = pa_x[k] * fpf0_236[k]
                   - f_5 * pc_x[k] * fpf1_236[k];

        t_237[k] = pa_x[k] * fpf0_237[k]
                   - f_5 * pc_x[k] * fpf1_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_x, pa_y, pc_x, pc_y, fpf0_150, \
                         fpf0_239, fpd_89, fpd_90, fpf1_150, fpf1_239, gpd_143, \
                         gpd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_11 * fpd_89[k]
                   + f_4 * pc_y[k] * gpd_143[k];

        t_239[k] = pa_x[k] * fpf0_239[k]
                   - f_5 * pc_x[k] * fpf1_239[k];

        t_240[k] = pa_y[k] * fpf0_150[k]
                   - f_5 * pc_y[k] * fpf1_150[k];

        t_241[k] = f_1 * fpd_90[k]
                   + f_4 * pc_y[k] * gpd_144[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_z, fpd_72, fpd_147, fpd_148, gsd_51, \
                         gsd_52, gpd_144, gpd_147, gpd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_11 * fpd_72[k]
                   + f_4 * pc_z[k] * gpd_144[k];

        t_243[k] = f_1 * fpd_147[k]
                   + f_1 * gsd_51[k]
                   + f_4 * pc_x[k] * gpd_147[k];

        t_244[k] = f_1 * fpd_148[k]
                   + f_1 * gsd_52[k]
                   + f_4 * pc_x[k] * gpd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_y, pc_y, pc_z, fpf0_155, fpd_75, \
                         fpd_93, fpd_95, fpf1_155, gpp0_73, gpp1_73, gpd_147, \
                         gpd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pa_y[k] * fpf0_155[k]
                   - f_5 * pc_y[k] * fpf1_155[k];

        t_246[k] = f_1 * fpd_93[k]
                   + f_2 * gpp0_73[k]
                   - f_3 * gpp1_73[k]
                   + f_4 * pc_y[k] * gpd_147[k];

        t_247[k] = f_11 * fpd_75[k]
                   + f_4 * pc_z[k] * gpd_147[k];

        t_248[k] = f_1 * fpd_95[k]
                   + f_4 * pc_y[k] * gpd_149[k];
    }
}

static auto
compute_prim_gpf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t fpf0, const size_t fpd,
                                                          const size_t fpf1, const size_t gsf0,
                                                          const size_t gsd, const size_t gsf1,
                                                          const size_t gpp0, const size_t gpp1,
                                                          const size_t gpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 1.0 / p;
    const auto f_8 = gamma / (p * q);
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);
    const auto f_11 = 1.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpf0_149 = buffer.data(dpf0 + 149);

    const auto *dpf1_149 = buffer.data(dpf1 + 149);

    const auto *fpf0_170 = buffer.data(fpf0 + 170);
    const auto *fpf0_172 = buffer.data(fpf0 + 172);
    const auto *fpf0_180 = buffer.data(fpf0 + 180);
    const auto *fpf0_181 = buffer.data(fpf0 + 181);
    const auto *fpf0_186 = buffer.data(fpf0 + 186);
    const auto *fpf0_190 = buffer.data(fpf0 + 190);
    const auto *fpf0_191 = buffer.data(fpf0 + 191);
    const auto *fpf0_196 = buffer.data(fpf0 + 196);
    const auto *fpf0_239 = buffer.data(fpf0 + 239);
    const auto *fpf0_256 = buffer.data(fpf0 + 256);
    const auto *fpf0_258 = buffer.data(fpf0 + 258);
    const auto *fpf0_259 = buffer.data(fpf0 + 259);
    const auto *fpf0_266 = buffer.data(fpf0 + 266);
    const auto *fpf0_267 = buffer.data(fpf0 + 267);
    const auto *fpf0_269 = buffer.data(fpf0 + 269);
    const auto *fpf0_286 = buffer.data(fpf0 + 286);
    const auto *fpf0_287 = buffer.data(fpf0 + 287);
    const auto *fpf0_289 = buffer.data(fpf0 + 289);
    const auto *fpf0_290 = buffer.data(fpf0 + 290);
    const auto *fpf0_296 = buffer.data(fpf0 + 296);
    const auto *fpf0_297 = buffer.data(fpf0 + 297);
    const auto *fpf0_299 = buffer.data(fpf0 + 299);

    const auto *fpd_77 = buffer.data(fpd + 77);
    const auto *fpd_78 = buffer.data(fpd + 78);
    const auto *fpd_81 = buffer.data(fpd + 81);
    const auto *fpd_90 = buffer.data(fpd + 90);
    const auto *fpd_95 = buffer.data(fpd + 95);
    const auto *fpd_96 = buffer.data(fpd + 96);
    const auto *fpd_102 = buffer.data(fpd + 102);
    const auto *fpd_107 = buffer.data(fpd + 107);
    const auto *fpd_111 = buffer.data(fpd + 111);
    const auto *fpd_113 = buffer.data(fpd + 113);
    const auto *fpd_117 = buffer.data(fpd + 117);
    const auto *fpd_119 = buffer.data(fpd + 119);
    const auto *fpd_123 = buffer.data(fpd + 123);
    const auto *fpd_125 = buffer.data(fpd + 125);
    const auto *fpd_129 = buffer.data(fpd + 129);
    const auto *fpd_131 = buffer.data(fpd + 131);
    const auto *fpd_137 = buffer.data(fpd + 137);
    const auto *fpd_141 = buffer.data(fpd + 141);
    const auto *fpd_143 = buffer.data(fpd + 143);
    const auto *fpd_149 = buffer.data(fpd + 149);
    const auto *fpd_150 = buffer.data(fpd + 150);
    const auto *fpd_153 = buffer.data(fpd + 153);
    const auto *fpd_154 = buffer.data(fpd + 154);
    const auto *fpd_155 = buffer.data(fpd + 155);
    const auto *fpd_159 = buffer.data(fpd + 159);
    const auto *fpd_160 = buffer.data(fpd + 160);
    const auto *fpd_161 = buffer.data(fpd + 161);
    const auto *fpd_162 = buffer.data(fpd + 162);
    const auto *fpd_165 = buffer.data(fpd + 165);
    const auto *fpd_167 = buffer.data(fpd + 167);
    const auto *fpd_171 = buffer.data(fpd + 171);
    const auto *fpd_173 = buffer.data(fpd + 173);
    const auto *fpd_174 = buffer.data(fpd + 174);
    const auto *fpd_177 = buffer.data(fpd + 177);
    const auto *fpd_179 = buffer.data(fpd + 179);

    const auto *fpf1_170 = buffer.data(fpf1 + 170);
    const auto *fpf1_172 = buffer.data(fpf1 + 172);
    const auto *fpf1_180 = buffer.data(fpf1 + 180);
    const auto *fpf1_181 = buffer.data(fpf1 + 181);
    const auto *fpf1_186 = buffer.data(fpf1 + 186);
    const auto *fpf1_190 = buffer.data(fpf1 + 190);
    const auto *fpf1_191 = buffer.data(fpf1 + 191);
    const auto *fpf1_196 = buffer.data(fpf1 + 196);
    const auto *fpf1_239 = buffer.data(fpf1 + 239);
    const auto *fpf1_256 = buffer.data(fpf1 + 256);
    const auto *fpf1_258 = buffer.data(fpf1 + 258);
    const auto *fpf1_259 = buffer.data(fpf1 + 259);
    const auto *fpf1_266 = buffer.data(fpf1 + 266);
    const auto *fpf1_267 = buffer.data(fpf1 + 267);
    const auto *fpf1_269 = buffer.data(fpf1 + 269);
    const auto *fpf1_286 = buffer.data(fpf1 + 286);
    const auto *fpf1_287 = buffer.data(fpf1 + 287);
    const auto *fpf1_289 = buffer.data(fpf1 + 289);
    const auto *fpf1_290 = buffer.data(fpf1 + 290);
    const auto *fpf1_296 = buffer.data(fpf1 + 296);
    const auto *fpf1_297 = buffer.data(fpf1 + 297);
    const auto *fpf1_299 = buffer.data(fpf1 + 299);

    const auto *gsf0_90 = buffer.data(gsf0 + 90);
    const auto *gsf0_100 = buffer.data(gsf0 + 100);
    const auto *gsf0_101 = buffer.data(gsf0 + 101);
    const auto *gsf0_106 = buffer.data(gsf0 + 106);
    const auto *gsf0_109 = buffer.data(gsf0 + 109);
    const auto *gsf0_112 = buffer.data(gsf0 + 112);
    const auto *gsf0_119 = buffer.data(gsf0 + 119);
    const auto *gsf0_120 = buffer.data(gsf0 + 120);
    const auto *gsf0_121 = buffer.data(gsf0 + 121);
    const auto *gsf0_122 = buffer.data(gsf0 + 122);
    const auto *gsf0_126 = buffer.data(gsf0 + 126);
    const auto *gsf0_129 = buffer.data(gsf0 + 129);

    const auto *gsd_48 = buffer.data(gsd + 48);
    const auto *gsd_54 = buffer.data(gsd + 54);
    const auto *gsd_56 = buffer.data(gsd + 56);
    const auto *gsd_57 = buffer.data(gsd + 57);
    const auto *gsd_59 = buffer.data(gsd + 59);
    const auto *gsd_60 = buffer.data(gsd + 60);
    const auto *gsd_61 = buffer.data(gsd + 61);
    const auto *gsd_63 = buffer.data(gsd + 63);
    const auto *gsd_64 = buffer.data(gsd + 64);
    const auto *gsd_65 = buffer.data(gsd + 65);
    const auto *gsd_68 = buffer.data(gsd + 68);
    const auto *gsd_69 = buffer.data(gsd + 69);
    const auto *gsd_70 = buffer.data(gsd + 70);
    const auto *gsd_71 = buffer.data(gsd + 71);
    const auto *gsd_72 = buffer.data(gsd + 72);
    const auto *gsd_73 = buffer.data(gsd + 73);
    const auto *gsd_74 = buffer.data(gsd + 74);
    const auto *gsd_75 = buffer.data(gsd + 75);
    const auto *gsd_76 = buffer.data(gsd + 76);
    const auto *gsd_77 = buffer.data(gsd + 77);

    const auto *gsf1_90 = buffer.data(gsf1 + 90);
    const auto *gsf1_100 = buffer.data(gsf1 + 100);
    const auto *gsf1_101 = buffer.data(gsf1 + 101);
    const auto *gsf1_106 = buffer.data(gsf1 + 106);
    const auto *gsf1_109 = buffer.data(gsf1 + 109);
    const auto *gsf1_112 = buffer.data(gsf1 + 112);
    const auto *gsf1_119 = buffer.data(gsf1 + 119);
    const auto *gsf1_120 = buffer.data(gsf1 + 120);
    const auto *gsf1_121 = buffer.data(gsf1 + 121);
    const auto *gsf1_122 = buffer.data(gsf1 + 122);
    const auto *gsf1_126 = buffer.data(gsf1 + 126);
    const auto *gsf1_129 = buffer.data(gsf1 + 129);

    const auto *gpp0_74 = buffer.data(gpp0 + 74);
    const auto *gpp0_75 = buffer.data(gpp0 + 75);
    const auto *gpp0_81 = buffer.data(gpp0 + 81);
    const auto *gpp0_82 = buffer.data(gpp0 + 82);
    const auto *gpp0_83 = buffer.data(gpp0 + 83);
    const auto *gpp0_93 = buffer.data(gpp0 + 93);
    const auto *gpp0_94 = buffer.data(gpp0 + 94);
    const auto *gpp0_95 = buffer.data(gpp0 + 95);
    const auto *gpp0_104 = buffer.data(gpp0 + 104);
    const auto *gpp0_105 = buffer.data(gpp0 + 105);
    const auto *gpp0_106 = buffer.data(gpp0 + 106);
    const auto *gpp0_107 = buffer.data(gpp0 + 107);
    const auto *gpp0_111 = buffer.data(gpp0 + 111);
    const auto *gpp0_112 = buffer.data(gpp0 + 112);
    const auto *gpp0_113 = buffer.data(gpp0 + 113);

    const auto *gpp1_74 = buffer.data(gpp1 + 74);
    const auto *gpp1_75 = buffer.data(gpp1 + 75);
    const auto *gpp1_81 = buffer.data(gpp1 + 81);
    const auto *gpp1_82 = buffer.data(gpp1 + 82);
    const auto *gpp1_83 = buffer.data(gpp1 + 83);
    const auto *gpp1_93 = buffer.data(gpp1 + 93);
    const auto *gpp1_94 = buffer.data(gpp1 + 94);
    const auto *gpp1_95 = buffer.data(gpp1 + 95);
    const auto *gpp1_104 = buffer.data(gpp1 + 104);
    const auto *gpp1_105 = buffer.data(gpp1 + 105);
    const auto *gpp1_106 = buffer.data(gpp1 + 106);
    const auto *gpp1_107 = buffer.data(gpp1 + 107);
    const auto *gpp1_111 = buffer.data(gpp1 + 111);
    const auto *gpp1_112 = buffer.data(gpp1 + 112);
    const auto *gpp1_113 = buffer.data(gpp1 + 113);

    const auto *gpd_149 = buffer.data(gpd + 149);
    const auto *gpd_150 = buffer.data(gpd + 150);
    const auto *gpd_153 = buffer.data(gpd + 153);
    const auto *gpd_154 = buffer.data(gpd + 154);
    const auto *gpd_155 = buffer.data(gpd + 155);
    const auto *gpd_156 = buffer.data(gpd + 156);
    const auto *gpd_159 = buffer.data(gpd + 159);
    const auto *gpd_160 = buffer.data(gpd + 160);
    const auto *gpd_161 = buffer.data(gpd + 161);
    const auto *gpd_162 = buffer.data(gpd + 162);
    const auto *gpd_164 = buffer.data(gpd + 164);
    const auto *gpd_165 = buffer.data(gpd + 165);
    const auto *gpd_166 = buffer.data(gpd + 166);
    const auto *gpd_167 = buffer.data(gpd + 167);
    const auto *gpd_168 = buffer.data(gpd + 168);
    const auto *gpd_170 = buffer.data(gpd + 170);
    const auto *gpd_171 = buffer.data(gpd + 171);
    const auto *gpd_173 = buffer.data(gpd + 173);
    const auto *gpd_174 = buffer.data(gpd + 174);
    const auto *gpd_176 = buffer.data(gpd + 176);
    const auto *gpd_177 = buffer.data(gpd + 177);
    const auto *gpd_179 = buffer.data(gpd + 179);
    const auto *gpd_180 = buffer.data(gpd + 180);
    const auto *gpd_183 = buffer.data(gpd + 183);
    const auto *gpd_184 = buffer.data(gpd + 184);
    const auto *gpd_185 = buffer.data(gpd + 185);
    const auto *gpd_186 = buffer.data(gpd + 186);
    const auto *gpd_187 = buffer.data(gpd + 187);
    const auto *gpd_189 = buffer.data(gpd + 189);
    const auto *gpd_190 = buffer.data(gpd + 190);
    const auto *gpd_191 = buffer.data(gpd + 191);
    const auto *gpd_192 = buffer.data(gpd + 192);
    const auto *gpd_195 = buffer.data(gpd + 195);
    const auto *gpd_196 = buffer.data(gpd + 196);
    const auto *gpd_197 = buffer.data(gpd + 197);
    const auto *gpd_201 = buffer.data(gpd + 201);
    const auto *gpd_202 = buffer.data(gpd + 202);
    const auto *gpd_203 = buffer.data(gpd + 203);
    const auto *gpd_206 = buffer.data(gpd + 206);
    const auto *gpd_207 = buffer.data(gpd + 207);
    const auto *gpd_208 = buffer.data(gpd + 208);
    const auto *gpd_209 = buffer.data(gpd + 209);
    const auto *gpd_210 = buffer.data(gpd + 210);
    const auto *gpd_211 = buffer.data(gpd + 211);
    const auto *gpd_212 = buffer.data(gpd + 212);
    const auto *gpd_213 = buffer.data(gpd + 213);
    const auto *gpd_214 = buffer.data(gpd + 214);
    const auto *gpd_215 = buffer.data(gpd + 215);
    const auto *gpd_219 = buffer.data(gpd + 219);
    const auto *gpd_220 = buffer.data(gpd + 220);
    const auto *gpd_221 = buffer.data(gpd + 221);
    const auto *gpd_222 = buffer.data(gpd + 222);
    const auto *gpd_223 = buffer.data(gpd + 223);
    const auto *gpd_224 = buffer.data(gpd + 224);
    const auto *gpd_225 = buffer.data(gpd + 225);
    const auto *gpd_226 = buffer.data(gpd + 226);
    const auto *gpd_227 = buffer.data(gpd + 227);

#pragma omp simd aligned(t_249, t_250, t_251, pc_x, pc_y, pc_z, fpd_77, fpd_96, fpd_150, \
                         gsd_48, gpp0_74, gpp0_75, gpp1_74, gpp1_75, gpd_149, \
                         gpd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_11 * fpd_77[k]
                   + f_2 * gpp0_74[k]
                   - f_3 * gpp1_74[k]
                   + f_4 * pc_z[k] * gpd_149[k];

        t_250[k] = f_1 * fpd_150[k]
                   + f_2 * gpp0_75[k]
                   - f_3 * gpp1_75[k]
                   + f_4 * pc_x[k] * gpd_150[k];

        t_251[k] = f_1 * fpd_96[k]
                   + f_1 * gsd_48[k]
                   + f_4 * pc_y[k] * gpd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, fpd_78, fpd_153, fpd_154, \
                         fpd_155, gpd_150, gpd_153, gpd_154, gpd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_11 * fpd_78[k]
                   + f_4 * pc_z[k] * gpd_150[k];

        t_253[k] = f_1 * fpd_153[k]
                   + f_4 * pc_x[k] * gpd_153[k];

        t_254[k] = f_1 * fpd_154[k]
                   + f_4 * pc_x[k] * gpd_154[k];

        t_255[k] = f_1 * fpd_155[k]
                   + f_4 * pc_x[k] * gpd_155[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_x, pc_x, pc_z, fpf0_256, fpf0_258, \
                         fpf0_259, fpd_81, fpf1_256, fpf1_258, fpf1_259, \
                         gpd_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pa_x[k] * fpf0_256[k]
                   - f_5 * pc_x[k] * fpf1_256[k];

        t_257[k] = f_11 * fpd_81[k]
                   + f_4 * pc_z[k] * gpd_153[k];

        t_258[k] = pa_x[k] * fpf0_258[k]
                   - f_5 * pc_x[k] * fpf1_258[k];

        t_259[k] = pa_x[k] * fpf0_259[k]
                   - f_5 * pc_x[k] * fpf1_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, fpf0_170, fpf0_172, \
                         fpd_102, fpd_159, fpf1_170, fpf1_172, gpd_156, \
                         gpd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * fpf0_170[k]
                   - f_5 * pc_y[k] * fpf1_170[k];

        t_261[k] = f_1 * fpd_102[k]
                   + f_4 * pc_y[k] * gpd_156[k];

        t_262[k] = pa_y[k] * fpf0_172[k]
                   - f_5 * pc_y[k] * fpf1_172[k];

        t_263[k] = f_1 * fpd_159[k]
                   + f_4 * pc_x[k] * gpd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_x, pc_x, fpf0_266, fpf0_267, fpd_160, \
                         fpd_161, fpf1_266, fpf1_267, gpd_160, \
                         gpd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_1 * fpd_160[k]
                   + f_4 * pc_x[k] * gpd_160[k];

        t_265[k] = f_1 * fpd_161[k]
                   + f_4 * pc_x[k] * gpd_161[k];

        t_266[k] = pa_x[k] * fpf0_266[k]
                   - f_5 * pc_x[k] * fpf1_266[k];

        t_267[k] = pa_x[k] * fpf0_267[k]
                   - f_5 * pc_x[k] * fpf1_267[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_x, pc_x, pc_y, fpf0_269, fpd_107, \
                         fpd_162, fpf1_269, gsd_54, gpp0_81, gpp1_81, gpd_161, \
                         gpd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_1 * fpd_107[k]
                   + f_4 * pc_y[k] * gpd_161[k];

        t_269[k] = pa_x[k] * fpf0_269[k]
                   - f_5 * pc_x[k] * fpf1_269[k];

        t_270[k] = f_1 * fpd_162[k]
                   + f_1 * gsd_54[k]
                   + f_2 * gpp0_81[k]
                   - f_3 * gpp1_81[k]
                   + f_4 * pc_x[k] * gpd_162[k];

        t_271[k] = f_4 * pc_y[k] * gpd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, fpd_90, fpd_165, \
                         fpd_167, gsd_57, gsd_59, gpd_162, gpd_164, gpd_165, \
                         gpd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_6 * fpd_90[k]
                   + f_4 * pc_z[k] * gpd_162[k];

        t_273[k] = f_1 * fpd_165[k]
                   + f_1 * gsd_57[k]
                   + f_4 * pc_x[k] * gpd_165[k];

        t_274[k] = f_4 * pc_y[k] * gpd_164[k];

        t_275[k] = f_1 * fpd_167[k]
                   + f_1 * gsd_59[k]
                   + f_4 * pc_x[k] * gpd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, fpd_95, gpp0_82, gpp0_83, \
                         gpp1_82, gpp1_83, gpd_165, gpd_166, gpd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_2 * gpp0_82[k]
                   - f_3 * gpp1_82[k]
                   + f_4 * pc_y[k] * gpd_165[k];

        t_277[k] = f_9 * gpp0_83[k]
                   - f_10 * gpp1_83[k]
                   + f_4 * pc_y[k] * gpd_166[k];

        t_278[k] = f_4 * pc_y[k] * gpd_167[k];

        t_279[k] = f_6 * fpd_95[k]
                   + f_2 * gpp0_83[k]
                   - f_3 * gpp1_83[k]
                   + f_4 * pc_z[k] * gpd_167[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pb_y, pc_x, pc_y, pc_z, fpd_96, fpd_171, \
                         gsf0_90, gsd_54, gsf1_90, gpd_168, gpd_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_y[k] * gsf0_90[k]
                   - f_5 * pc_y[k] * gsf1_90[k];

        t_281[k] = f_1 * gsd_54[k]
                   + f_4 * pc_y[k] * gpd_168[k];

        t_282[k] = f_6 * fpd_96[k]
                   + f_4 * pc_z[k] * gpd_168[k];

        t_283[k] = f_1 * fpd_171[k]
                   + f_4 * pc_x[k] * gpd_171[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pa_x, pc_x, pc_y, fpf0_286, fpf0_287, \
                         fpd_173, fpf1_286, fpf1_287, gsd_56, gpd_170, \
                         gpd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_1 * gsd_56[k]
                   + f_4 * pc_y[k] * gpd_170[k];

        t_285[k] = f_1 * fpd_173[k]
                   + f_4 * pc_x[k] * gpd_173[k];

        t_286[k] = pa_x[k] * fpf0_286[k]
                   - f_5 * pc_x[k] * fpf1_286[k];

        t_287[k] = pa_x[k] * fpf0_287[k]
                   - f_5 * pc_x[k] * fpf1_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_x, pc_x, pc_y, fpf0_289, fpf0_290, \
                         fpd_174, fpf1_289, fpf1_290, gsd_59, gpd_173, \
                         gpd_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_1 * gsd_59[k]
                   + f_4 * pc_y[k] * gpd_173[k];

        t_289[k] = pa_x[k] * fpf0_289[k]
                   - f_5 * pc_x[k] * fpf1_289[k];

        t_290[k] = pa_x[k] * fpf0_290[k]
                   + f_6 * fpd_174[k]
                   - f_5 * pc_x[k] * fpf1_290[k];

        t_291[k] = f_4 * pc_y[k] * gpd_174[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, fpd_102, fpd_177, \
                         fpd_179, gsd_54, gpd_174, gpd_176, gpd_177, \
                         gpd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_6 * fpd_102[k]
                   + f_1 * gsd_54[k]
                   + f_4 * pc_z[k] * gpd_174[k];

        t_293[k] = f_1 * fpd_177[k]
                   + f_4 * pc_x[k] * gpd_177[k];

        t_294[k] = f_4 * pc_y[k] * gpd_176[k];

        t_295[k] = f_1 * fpd_179[k]
                   + f_4 * pc_x[k] * gpd_179[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_x, pc_x, pc_y, fpf0_296, fpf0_297, \
                         fpf0_299, fpf1_296, fpf1_297, fpf1_299, \
                         gpd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = pa_x[k] * fpf0_296[k]
                   - f_5 * pc_x[k] * fpf1_296[k];

        t_297[k] = pa_x[k] * fpf0_297[k]
                   - f_5 * pc_x[k] * fpf1_297[k];

        t_298[k] = f_4 * pc_y[k] * gpd_179[k];

        t_299[k] = pa_x[k] * fpf0_299[k]
                   - f_5 * pc_x[k] * fpf1_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pb_x, pc_x, pc_z, gsf0_100, gsf0_101, \
                         gsd_60, gsd_61, gsd_63, gsf1_100, gsf1_101, gpd_180, \
                         gpd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pb_x[k] * gsf0_100[k]
                   + f_6 * gsd_60[k]
                   - f_5 * pc_x[k] * gsf1_100[k];

        t_301[k] = pb_x[k] * gsf0_101[k]
                   + f_11 * gsd_61[k]
                   - f_5 * pc_x[k] * gsf1_101[k];

        t_302[k] = f_4 * pc_z[k] * gpd_180[k];

        t_303[k] = f_1 * gsd_63[k]
                   + f_4 * pc_x[k] * gpd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pb_x, pc_x, pc_z, gsf0_106, gsd_64, \
                         gsd_65, gsf1_106, gpd_183, gpd_184, gpd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_1 * gsd_64[k]
                   + f_4 * pc_x[k] * gpd_184[k];

        t_305[k] = f_1 * gsd_65[k]
                   + f_4 * pc_x[k] * gpd_185[k];

        t_306[k] = pb_x[k] * gsf0_106[k]
                   - f_5 * pc_x[k] * gsf1_106[k];

        t_307[k] = f_4 * pc_z[k] * gpd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pb_x, pc_x, pc_y, fpd_113, gsf0_109, gsf1_109, \
                         gpp0_93, gpp1_93, gpd_185, gpd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * fpd_113[k]
                   + f_4 * pc_y[k] * gpd_185[k];

        t_309[k] = pb_x[k] * gsf0_109[k]
                   - f_5 * pc_x[k] * gsf1_109[k];

        t_310[k] = f_2 * gpp0_93[k]
                   - f_3 * gpp1_93[k]
                   + f_4 * pc_x[k] * gpd_186[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, t_315, pc_x, pc_z, gpp0_94, gpp1_94, \
                         gpd_186, gpd_187, gpd_189, gpd_190, gpd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_9 * gpp0_94[k]
                   - f_10 * gpp1_94[k]
                   + f_4 * pc_x[k] * gpd_187[k];

        t_312[k] = f_4 * pc_z[k] * gpd_186[k];

        t_313[k] = f_4 * pc_x[k] * gpd_189[k];

        t_314[k] = f_4 * pc_x[k] * gpd_190[k];

        t_315[k] = f_4 * pc_x[k] * gpd_191[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_y, pc_z, fpd_117, fpd_119, gsd_63, \
                         gsd_65, gpp0_94, gpp0_95, gpp1_94, gpp1_95, gpd_189, \
                         gpd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_0 * fpd_117[k]
                   + f_1 * gsd_63[k]
                   + f_2 * gpp0_94[k]
                   - f_3 * gpp1_94[k]
                   + f_4 * pc_y[k] * gpd_189[k];

        t_317[k] = f_4 * pc_z[k] * gpd_189[k];

        t_318[k] = f_0 * fpd_119[k]
                   + f_1 * gsd_65[k]
                   + f_4 * pc_y[k] * gpd_191[k];

        t_319[k] = f_2 * gpp0_95[k]
                   - f_3 * gpp1_95[k]
                   + f_4 * pc_z[k] * gpd_191[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, pb_z, pc_x, pc_z, gsf0_100, \
                         gsf0_101, gsd_60, gsf1_100, gsf1_101, gpd_192, gpd_195, \
                         gpd_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pb_z[k] * gsf0_100[k]
                   - f_5 * pc_z[k] * gsf1_100[k];

        t_321[k] = pb_z[k] * gsf0_101[k]
                   - f_5 * pc_z[k] * gsf1_101[k];

        t_322[k] = f_1 * gsd_60[k]
                   + f_4 * pc_z[k] * gpd_192[k];

        t_323[k] = f_4 * pc_x[k] * gpd_195[k];

        t_324[k] = f_4 * pc_x[k] * gpd_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pb_z, pc_x, pc_y, pc_z, fpd_125, \
                         gsf0_106, gsd_63, gsf1_106, gpd_195, gpd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_4 * pc_x[k] * gpd_197[k];

        t_326[k] = pb_z[k] * gsf0_106[k]
                   - f_5 * pc_z[k] * gsf1_106[k];

        t_327[k] = f_1 * gsd_63[k]
                   + f_4 * pc_z[k] * gpd_195[k];

        t_328[k] = f_0 * fpd_125[k]
                   + f_4 * pc_y[k] * gpd_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pb_z, pc_z, fpf0_180, fpf0_181, fpf1_180, \
                         fpf1_181, gsf0_109, gsd_65, gsf1_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pb_z[k] * gsf0_109[k]
                   + f_6 * gsd_65[k]
                   - f_5 * pc_z[k] * gsf1_109[k];

        t_330[k] = pa_z[k] * fpf0_180[k]
                   - f_5 * pc_z[k] * fpf1_180[k];

        t_331[k] = pa_z[k] * fpf0_181[k]
                   - f_5 * pc_z[k] * fpf1_181[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pb_x, pc_x, gsf0_112, gsd_68, gsd_69, \
                         gsd_70, gsd_71, gsf1_112, gpd_201, gpd_202, \
                         gpd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_x[k] * gsf0_112[k]
                   + f_11 * gsd_68[k]
                   - f_5 * pc_x[k] * gsf1_112[k];

        t_333[k] = f_1 * gsd_69[k]
                   + f_4 * pc_x[k] * gpd_201[k];

        t_334[k] = f_1 * gsd_70[k]
                   + f_4 * pc_x[k] * gpd_202[k];

        t_335[k] = f_1 * gsd_71[k]
                   + f_4 * pc_x[k] * gpd_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_z, pc_y, pc_z, fpf0_186, fpd_111, fpd_131, \
                         fpf1_186, gpd_201, gpd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_z[k] * fpf0_186[k]
                   - f_5 * pc_z[k] * fpf1_186[k];

        t_337[k] = f_1 * fpd_111[k]
                   + f_4 * pc_z[k] * gpd_201[k];

        t_338[k] = f_6 * fpd_131[k]
                   + f_4 * pc_y[k] * gpd_203[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pb_x, pc_x, pc_z, fpf0_190, fpf0_191, \
                         fpf1_190, fpf1_191, gsf0_119, gsf1_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pb_x[k] * gsf0_119[k]
                   - f_5 * pc_x[k] * gsf1_119[k];

        t_340[k] = pa_z[k] * fpf0_190[k]
                   - f_5 * pc_z[k] * fpf1_190[k];

        t_341[k] = pa_z[k] * fpf0_191[k]
                   - f_5 * pc_z[k] * fpf1_191[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, pa_z, pc_x, pc_z, fpf0_196, \
                         fpf1_196, gpp0_104, gpp1_104, gpd_206, gpd_207, gpd_208, \
                         gpd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_9 * gpp0_104[k]
                   - f_10 * gpp1_104[k]
                   + f_4 * pc_x[k] * gpd_206[k];

        t_343[k] = f_4 * pc_x[k] * gpd_207[k];

        t_344[k] = f_4 * pc_x[k] * gpd_208[k];

        t_345[k] = f_4 * pc_x[k] * gpd_209[k];

        t_346[k] = pa_z[k] * fpf0_196[k]
                   - f_5 * pc_z[k] * fpf1_196[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pc_y, pc_z, fpd_117, fpd_119, fpd_137, gsd_71, \
                         gpp0_104, gpp1_104, gpd_207, gpd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_1 * fpd_117[k]
                   + f_4 * pc_z[k] * gpd_207[k];

        t_348[k] = f_6 * fpd_137[k]
                   + f_1 * gsd_71[k]
                   + f_4 * pc_y[k] * gpd_209[k];

        t_349[k] = f_1 * fpd_119[k]
                   + f_2 * gpp0_104[k]
                   - f_3 * gpp1_104[k]
                   + f_4 * pc_z[k] * gpd_209[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pc_x, gpp0_105, gpp0_106, gpp0_107, \
                         gpp1_105, gpp1_106, gpp1_107, gpd_210, gpd_211, gpd_212, \
                         gpd_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_2 * gpp0_105[k]
                   - f_3 * gpp1_105[k]
                   + f_4 * pc_x[k] * gpd_210[k];

        t_351[k] = f_9 * gpp0_106[k]
                   - f_10 * gpp1_106[k]
                   + f_4 * pc_x[k] * gpd_211[k];

        t_352[k] = f_9 * gpp0_107[k]
                   - f_10 * gpp1_107[k]
                   + f_4 * pc_x[k] * gpd_212[k];

        t_353[k] = f_4 * pc_x[k] * gpd_213[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pc_x, pc_y, pc_z, fpd_123, fpd_141, \
                         gsd_69, gpp0_106, gpp1_106, gpd_213, gpd_214, \
                         gpd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_4 * pc_x[k] * gpd_214[k];

        t_355[k] = f_4 * pc_x[k] * gpd_215[k];

        t_356[k] = f_6 * fpd_141[k]
                   + f_2 * gpp0_106[k]
                   - f_3 * gpp1_106[k]
                   + f_4 * pc_y[k] * gpd_213[k];

        t_357[k] = f_1 * fpd_123[k]
                   + f_1 * gsd_69[k]
                   + f_4 * pc_z[k] * gpd_213[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pa_y, pb_x, pc_x, pc_y, dpf0_149, dpf1_149, \
                         fpf0_239, fpd_143, fpf1_239, gsf0_120, gsd_72, gsf1_120, \
                         gpd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_6 * fpd_143[k]
                   + f_4 * pc_y[k] * gpd_215[k];

        t_359[k] = f_7 * dpf0_149[k]
                   - f_8 * dpf1_149[k]
                   + pa_y[k] * fpf0_239[k]
                   - f_5 * pc_y[k] * fpf1_239[k];

        t_360[k] = pb_x[k] * gsf0_120[k]
                   + f_6 * gsd_72[k]
                   - f_5 * pc_x[k] * gsf1_120[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pb_x, pc_x, gsf0_121, gsf0_122, gsd_73, \
                         gsd_74, gsd_75, gsd_76, gsf1_121, gsf1_122, gpd_219, \
                         gpd_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = pb_x[k] * gsf0_121[k]
                   + f_11 * gsd_73[k]
                   - f_5 * pc_x[k] * gsf1_121[k];

        t_362[k] = pb_x[k] * gsf0_122[k]
                   + f_11 * gsd_74[k]
                   - f_5 * pc_x[k] * gsf1_122[k];

        t_363[k] = f_1 * gsd_75[k]
                   + f_4 * pc_x[k] * gpd_219[k];

        t_364[k] = f_1 * gsd_76[k]
                   + f_4 * pc_x[k] * gpd_220[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pb_x, pc_x, pc_y, pc_z, fpd_129, fpd_149, \
                         gsf0_126, gsd_77, gsf1_126, gpd_219, gpd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_1 * gsd_77[k]
                   + f_4 * pc_x[k] * gpd_221[k];

        t_366[k] = pb_x[k] * gsf0_126[k]
                   - f_5 * pc_x[k] * gsf1_126[k];

        t_367[k] = f_11 * fpd_129[k]
                   + f_4 * pc_z[k] * gpd_219[k];

        t_368[k] = f_11 * fpd_149[k]
                   + f_4 * pc_y[k] * gpd_221[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pb_x, pc_x, gsf0_129, gsf1_129, gpp0_111, \
                         gpp0_112, gpp1_111, gpp1_112, gpd_222, \
                         gpd_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_x[k] * gsf0_129[k]
                   - f_5 * pc_x[k] * gsf1_129[k];

        t_370[k] = f_2 * gpp0_111[k]
                   - f_3 * gpp1_111[k]
                   + f_4 * pc_x[k] * gpd_222[k];

        t_371[k] = f_9 * gpp0_112[k]
                   - f_10 * gpp1_112[k]
                   + f_4 * pc_x[k] * gpd_223[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pc_x, gpp0_113, gpp1_113, gpd_224, \
                         gpd_225, gpd_226, gpd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_9 * gpp0_113[k]
                   - f_10 * gpp1_113[k]
                   + f_4 * pc_x[k] * gpd_224[k];

        t_373[k] = f_4 * pc_x[k] * gpd_225[k];

        t_374[k] = f_4 * pc_x[k] * gpd_226[k];

        t_375[k] = f_4 * pc_x[k] * gpd_227[k];
    }
}

static auto
compute_prim_gpf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t fpf0, const size_t fpd,
                                                          const size_t fpf1, const size_t gsf0,
                                                          const size_t gsd, const size_t gsf1,
                                                          const size_t gpp0, const size_t gpp1,
                                                          const size_t gpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);
    const auto f_11 = 1.0 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpf0_106 = buffer.data(dpf0 + 106);
    const auto *dpf0_179 = buffer.data(dpf0 + 179);

    const auto *dpf1_106 = buffer.data(dpf1 + 106);
    const auto *dpf1_179 = buffer.data(dpf1 + 179);

    const auto *fpf0_226 = buffer.data(fpf0 + 226);
    const auto *fpf0_269 = buffer.data(fpf0 + 269);
    const auto *fpf0_270 = buffer.data(fpf0 + 270);
    const auto *fpf0_271 = buffer.data(fpf0 + 271);
    const auto *fpf0_272 = buffer.data(fpf0 + 272);
    const auto *fpf0_279 = buffer.data(fpf0 + 279);
    const auto *fpf0_290 = buffer.data(fpf0 + 290);
    const auto *fpf0_292 = buffer.data(fpf0 + 292);
    const auto *fpf0_296 = buffer.data(fpf0 + 296);
    const auto *fpf0_299 = buffer.data(fpf0 + 299);

    const auto *fpd_135 = buffer.data(fpd + 135);
    const auto *fpd_137 = buffer.data(fpd + 137);
    const auto *fpd_141 = buffer.data(fpd + 141);
    const auto *fpd_147 = buffer.data(fpd + 147);
    const auto *fpd_153 = buffer.data(fpd + 153);
    const auto *fpd_155 = buffer.data(fpd + 155);
    const auto *fpd_159 = buffer.data(fpd + 159);
    const auto *fpd_161 = buffer.data(fpd + 161);
    const auto *fpd_162 = buffer.data(fpd + 162);
    const auto *fpd_167 = buffer.data(fpd + 167);
    const auto *fpd_171 = buffer.data(fpd + 171);
    const auto *fpd_173 = buffer.data(fpd + 173);
    const auto *fpd_177 = buffer.data(fpd + 177);
    const auto *fpd_179 = buffer.data(fpd + 179);

    const auto *fpf1_226 = buffer.data(fpf1 + 226);
    const auto *fpf1_269 = buffer.data(fpf1 + 269);
    const auto *fpf1_270 = buffer.data(fpf1 + 270);
    const auto *fpf1_271 = buffer.data(fpf1 + 271);
    const auto *fpf1_272 = buffer.data(fpf1 + 272);
    const auto *fpf1_279 = buffer.data(fpf1 + 279);
    const auto *fpf1_290 = buffer.data(fpf1 + 290);
    const auto *fpf1_292 = buffer.data(fpf1 + 292);
    const auto *fpf1_296 = buffer.data(fpf1 + 296);
    const auto *fpf1_299 = buffer.data(fpf1 + 299);

    const auto *gsf0_136 = buffer.data(gsf0 + 136);
    const auto *gsf0_140 = buffer.data(gsf0 + 140);
    const auto *gsf0_142 = buffer.data(gsf0 + 142);
    const auto *gsf0_146 = buffer.data(gsf0 + 146);
    const auto *gsf0_147 = buffer.data(gsf0 + 147);
    const auto *gsf0_149 = buffer.data(gsf0 + 149);

    const auto *gsd_75 = buffer.data(gsd + 75);
    const auto *gsd_77 = buffer.data(gsd + 77);
    const auto *gsd_81 = buffer.data(gsd + 81);
    const auto *gsd_82 = buffer.data(gsd + 82);
    const auto *gsd_83 = buffer.data(gsd + 83);
    const auto *gsd_84 = buffer.data(gsd + 84);
    const auto *gsd_86 = buffer.data(gsd + 86);
    const auto *gsd_87 = buffer.data(gsd + 87);
    const auto *gsd_88 = buffer.data(gsd + 88);
    const auto *gsd_89 = buffer.data(gsd + 89);

    const auto *gsf1_136 = buffer.data(gsf1 + 136);
    const auto *gsf1_140 = buffer.data(gsf1 + 140);
    const auto *gsf1_142 = buffer.data(gsf1 + 142);
    const auto *gsf1_146 = buffer.data(gsf1 + 146);
    const auto *gsf1_147 = buffer.data(gsf1 + 147);
    const auto *gsf1_149 = buffer.data(gsf1 + 149);

    const auto *gpp0_113 = buffer.data(gpp0 + 113);
    const auto *gpp0_114 = buffer.data(gpp0 + 114);
    const auto *gpp0_115 = buffer.data(gpp0 + 115);
    const auto *gpp0_116 = buffer.data(gpp0 + 116);
    const auto *gpp0_120 = buffer.data(gpp0 + 120);
    const auto *gpp0_121 = buffer.data(gpp0 + 121);
    const auto *gpp0_122 = buffer.data(gpp0 + 122);
    const auto *gpp0_124 = buffer.data(gpp0 + 124);
    const auto *gpp0_132 = buffer.data(gpp0 + 132);
    const auto *gpp0_133 = buffer.data(gpp0 + 133);
    const auto *gpp0_134 = buffer.data(gpp0 + 134);

    const auto *gpp1_113 = buffer.data(gpp1 + 113);
    const auto *gpp1_114 = buffer.data(gpp1 + 114);
    const auto *gpp1_115 = buffer.data(gpp1 + 115);
    const auto *gpp1_116 = buffer.data(gpp1 + 116);
    const auto *gpp1_120 = buffer.data(gpp1 + 120);
    const auto *gpp1_121 = buffer.data(gpp1 + 121);
    const auto *gpp1_122 = buffer.data(gpp1 + 122);
    const auto *gpp1_124 = buffer.data(gpp1 + 124);
    const auto *gpp1_132 = buffer.data(gpp1 + 132);
    const auto *gpp1_133 = buffer.data(gpp1 + 133);
    const auto *gpp1_134 = buffer.data(gpp1 + 134);

    const auto *gpd_225 = buffer.data(gpd + 225);
    const auto *gpd_227 = buffer.data(gpd + 227);
    const auto *gpd_228 = buffer.data(gpd + 228);
    const auto *gpd_229 = buffer.data(gpd + 229);
    const auto *gpd_230 = buffer.data(gpd + 230);
    const auto *gpd_231 = buffer.data(gpd + 231);
    const auto *gpd_232 = buffer.data(gpd + 232);
    const auto *gpd_233 = buffer.data(gpd + 233);
    const auto *gpd_237 = buffer.data(gpd + 237);
    const auto *gpd_238 = buffer.data(gpd + 238);
    const auto *gpd_239 = buffer.data(gpd + 239);
    const auto *gpd_240 = buffer.data(gpd + 240);
    const auto *gpd_241 = buffer.data(gpd + 241);
    const auto *gpd_242 = buffer.data(gpd + 242);
    const auto *gpd_243 = buffer.data(gpd + 243);
    const auto *gpd_244 = buffer.data(gpd + 244);
    const auto *gpd_245 = buffer.data(gpd + 245);
    const auto *gpd_247 = buffer.data(gpd + 247);
    const auto *gpd_249 = buffer.data(gpd + 249);
    const auto *gpd_250 = buffer.data(gpd + 250);
    const auto *gpd_251 = buffer.data(gpd + 251);
    const auto *gpd_252 = buffer.data(gpd + 252);
    const auto *gpd_255 = buffer.data(gpd + 255);
    const auto *gpd_256 = buffer.data(gpd + 256);
    const auto *gpd_257 = buffer.data(gpd + 257);
    const auto *gpd_258 = buffer.data(gpd + 258);
    const auto *gpd_261 = buffer.data(gpd + 261);
    const auto *gpd_262 = buffer.data(gpd + 262);
    const auto *gpd_263 = buffer.data(gpd + 263);
    const auto *gpd_264 = buffer.data(gpd + 264);
    const auto *gpd_266 = buffer.data(gpd + 266);
    const auto *gpd_267 = buffer.data(gpd + 267);
    const auto *gpd_268 = buffer.data(gpd + 268);
    const auto *gpd_269 = buffer.data(gpd + 269);

#pragma omp simd aligned(t_376, t_377, t_378, pa_z, pc_y, pc_z, dpf0_106, dpf1_106, fpf0_226, \
                         fpd_135, fpd_155, fpf1_226, gsd_77, gpd_225, \
                         gpd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_12 * dpf0_106[k]
                   - f_13 * dpf1_106[k]
                   + pa_z[k] * fpf0_226[k]
                   - f_5 * pc_z[k] * fpf1_226[k];

        t_377[k] = f_11 * fpd_135[k]
                   + f_4 * pc_z[k] * gpd_225[k];

        t_378[k] = f_11 * fpd_155[k]
                   + f_1 * gsd_77[k]
                   + f_4 * pc_y[k] * gpd_227[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, pc_x, pc_z, fpd_137, gpp0_113, gpp0_114, \
                         gpp0_115, gpp1_113, gpp1_114, gpp1_115, gpd_227, gpd_228, \
                         gpd_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_11 * fpd_137[k]
                   + f_2 * gpp0_113[k]
                   - f_3 * gpp1_113[k]
                   + f_4 * pc_z[k] * gpd_227[k];

        t_380[k] = f_2 * gpp0_114[k]
                   - f_3 * gpp1_114[k]
                   + f_4 * pc_x[k] * gpd_228[k];

        t_381[k] = f_9 * gpp0_115[k]
                   - f_10 * gpp1_115[k]
                   + f_4 * pc_x[k] * gpd_229[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_y, fpd_159, gpp0_115, \
                         gpp0_116, gpp1_115, gpp1_116, gpd_230, gpd_231, gpd_232, \
                         gpd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_9 * gpp0_116[k]
                   - f_10 * gpp1_116[k]
                   + f_4 * pc_x[k] * gpd_230[k];

        t_383[k] = f_4 * pc_x[k] * gpd_231[k];

        t_384[k] = f_4 * pc_x[k] * gpd_232[k];

        t_385[k] = f_4 * pc_x[k] * gpd_233[k];

        t_386[k] = f_11 * fpd_159[k]
                   + f_2 * gpp0_115[k]
                   - f_3 * gpp1_115[k]
                   + f_4 * pc_y[k] * gpd_231[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, pa_y, pc_y, pc_z, dpf0_179, dpf1_179, fpf0_269, \
                         fpd_141, fpd_161, fpf1_269, gsd_75, gpd_231, \
                         gpd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_11 * fpd_141[k]
                   + f_1 * gsd_75[k]
                   + f_4 * pc_z[k] * gpd_231[k];

        t_388[k] = f_11 * fpd_161[k]
                   + f_4 * pc_y[k] * gpd_233[k];

        t_389[k] = f_12 * dpf0_179[k]
                   - f_13 * dpf1_179[k]
                   + pa_y[k] * fpf0_269[k]
                   - f_5 * pc_y[k] * fpf1_269[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pa_y, pc_x, pc_y, fpf0_270, fpf0_271, \
                         fpf0_272, fpd_162, fpf1_270, fpf1_271, fpf1_272, gsd_81, \
                         gpd_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = pa_y[k] * fpf0_270[k]
                   - f_5 * pc_y[k] * fpf1_270[k];

        t_391[k] = pa_y[k] * fpf0_271[k]
                   + f_1 * fpd_162[k]
                   - f_5 * pc_y[k] * fpf1_271[k];

        t_392[k] = pa_y[k] * fpf0_272[k]
                   - f_5 * pc_y[k] * fpf1_272[k];

        t_393[k] = f_1 * gsd_81[k]
                   + f_4 * pc_x[k] * gpd_237[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_x, pc_x, pc_z, fpd_147, gsf0_136, \
                         gsd_82, gsd_83, gsf1_136, gpd_237, gpd_238, \
                         gpd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_1 * gsd_82[k]
                   + f_4 * pc_x[k] * gpd_238[k];

        t_395[k] = f_1 * gsd_83[k]
                   + f_4 * pc_x[k] * gpd_239[k];

        t_396[k] = pb_x[k] * gsf0_136[k]
                   - f_5 * pc_x[k] * gsf1_136[k];

        t_397[k] = f_6 * fpd_147[k]
                   + f_4 * pc_z[k] * gpd_237[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_y, pc_x, pc_y, fpf0_279, fpd_167, fpf1_279, \
                         gpp0_120, gpp1_120, gpd_239, gpd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_1 * fpd_167[k]
                   + f_4 * pc_y[k] * gpd_239[k];

        t_399[k] = pa_y[k] * fpf0_279[k]
                   - f_5 * pc_y[k] * fpf1_279[k];

        t_400[k] = f_2 * gpp0_120[k]
                   - f_3 * gpp1_120[k]
                   + f_4 * pc_x[k] * gpd_240[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, pc_x, gpp0_121, gpp0_122, \
                         gpp1_121, gpp1_122, gpd_241, gpd_242, gpd_243, gpd_244, \
                         gpd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_9 * gpp0_121[k]
                   - f_10 * gpp1_121[k]
                   + f_4 * pc_x[k] * gpd_241[k];

        t_402[k] = f_9 * gpp0_122[k]
                   - f_10 * gpp1_122[k]
                   + f_4 * pc_x[k] * gpd_242[k];

        t_403[k] = f_4 * pc_x[k] * gpd_243[k];

        t_404[k] = f_4 * pc_x[k] * gpd_244[k];

        t_405[k] = f_4 * pc_x[k] * gpd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, pc_y, pc_z, fpd_153, fpd_171, fpd_173, gsd_81, \
                         gsd_83, gpp0_121, gpp1_121, gpd_243, gpd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_1 * fpd_171[k]
                   + f_1 * gsd_81[k]
                   + f_2 * gpp0_121[k]
                   - f_3 * gpp1_121[k]
                   + f_4 * pc_y[k] * gpd_243[k];

        t_407[k] = f_6 * fpd_153[k]
                   + f_4 * pc_z[k] * gpd_243[k];

        t_408[k] = f_1 * fpd_173[k]
                   + f_1 * gsd_83[k]
                   + f_4 * pc_y[k] * gpd_245[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, pa_y, pc_x, pc_y, pc_z, fpf0_290, fpd_155, \
                         fpf1_290, gpp0_122, gpp0_124, gpp1_122, gpp1_124, gpd_245, \
                         gpd_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_6 * fpd_155[k]
                   + f_2 * gpp0_122[k]
                   - f_3 * gpp1_122[k]
                   + f_4 * pc_z[k] * gpd_245[k];

        t_410[k] = pa_y[k] * fpf0_290[k]
                   - f_5 * pc_y[k] * fpf1_290[k];

        t_411[k] = f_9 * gpp0_124[k]
                   - f_10 * gpp1_124[k]
                   + f_4 * pc_x[k] * gpd_247[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, pa_y, pc_x, pc_y, fpf0_292, \
                         fpf0_296, fpd_177, fpf1_292, fpf1_296, gpd_249, gpd_250, \
                         gpd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = pa_y[k] * fpf0_292[k]
                   - f_5 * pc_y[k] * fpf1_292[k];

        t_413[k] = f_4 * pc_x[k] * gpd_249[k];

        t_414[k] = f_4 * pc_x[k] * gpd_250[k];

        t_415[k] = f_4 * pc_x[k] * gpd_251[k];

        t_416[k] = pa_y[k] * fpf0_296[k]
                   + f_6 * fpd_177[k]
                   - f_5 * pc_y[k] * fpf1_296[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, pc_z, fpf0_299, fpd_159, fpd_179, \
                         fpf1_299, gsd_81, gpd_249, gpd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_6 * fpd_159[k]
                   + f_1 * gsd_81[k]
                   + f_4 * pc_z[k] * gpd_249[k];

        t_418[k] = f_1 * fpd_179[k]
                   + f_4 * pc_y[k] * gpd_251[k];

        t_419[k] = pa_y[k] * fpf0_299[k]
                   - f_5 * pc_y[k] * fpf1_299[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_x, pc_x, pc_y, gsf0_140, gsf0_142, \
                         gsd_84, gsd_86, gsd_87, gsf1_140, gsf1_142, gpd_252, \
                         gpd_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pb_x[k] * gsf0_140[k]
                   + f_6 * gsd_84[k]
                   - f_5 * pc_x[k] * gsf1_140[k];

        t_421[k] = f_4 * pc_y[k] * gpd_252[k];

        t_422[k] = pb_x[k] * gsf0_142[k]
                   + f_11 * gsd_86[k]
                   - f_5 * pc_x[k] * gsf1_142[k];

        t_423[k] = f_1 * gsd_87[k]
                   + f_4 * pc_x[k] * gpd_255[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, pb_x, pc_x, pc_y, gsf0_146, \
                         gsf0_147, gsd_88, gsd_89, gsf1_146, gsf1_147, gpd_256, \
                         gpd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_1 * gsd_88[k]
                   + f_4 * pc_x[k] * gpd_256[k];

        t_425[k] = f_1 * gsd_89[k]
                   + f_4 * pc_x[k] * gpd_257[k];

        t_426[k] = pb_x[k] * gsf0_146[k]
                   - f_5 * pc_x[k] * gsf1_146[k];

        t_427[k] = pb_x[k] * gsf0_147[k]
                   - f_5 * pc_x[k] * gsf1_147[k];

        t_428[k] = f_4 * pc_y[k] * gpd_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pb_x, pb_y, pc_x, pc_y, gsf0_140, \
                         gsf0_142, gsf0_149, gsd_84, gsf1_140, gsf1_142, gsf1_149, \
                         gpd_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = pb_x[k] * gsf0_149[k]
                   - f_5 * pc_x[k] * gsf1_149[k];

        t_430[k] = pb_y[k] * gsf0_140[k]
                   - f_5 * pc_y[k] * gsf1_140[k];

        t_431[k] = f_1 * gsd_84[k]
                   + f_4 * pc_y[k] * gpd_258[k];

        t_432[k] = pb_y[k] * gsf0_142[k]
                   - f_5 * pc_y[k] * gsf1_142[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pb_y, pc_x, pc_y, gsf0_146, gsd_87, \
                         gsf1_146, gpd_261, gpd_262, gpd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_4 * pc_x[k] * gpd_261[k];

        t_434[k] = f_4 * pc_x[k] * gpd_262[k];

        t_435[k] = f_4 * pc_x[k] * gpd_263[k];

        t_436[k] = pb_y[k] * gsf0_146[k]
                   + f_6 * gsd_87[k]
                   - f_5 * pc_y[k] * gsf1_146[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pb_y, pc_y, gsf0_147, gsf0_149, gsd_88, gsd_89, \
                         gsf1_147, gsf1_149, gpd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = pb_y[k] * gsf0_147[k]
                   + f_11 * gsd_88[k]
                   - f_5 * pc_y[k] * gsf1_147[k];

        t_438[k] = f_1 * gsd_89[k]
                   + f_4 * pc_y[k] * gpd_263[k];

        t_439[k] = pb_y[k] * gsf0_149[k]
                   - f_5 * pc_y[k] * gsf1_149[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, pc_y, gpp0_132, gpp0_134, \
                         gpp1_132, gpp1_134, gpd_264, gpd_266, gpd_267, \
                         gpd_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_2 * gpp0_132[k]
                   - f_3 * gpp1_132[k]
                   + f_4 * pc_x[k] * gpd_264[k];

        t_441[k] = f_4 * pc_y[k] * gpd_264[k];

        t_442[k] = f_9 * gpp0_134[k]
                   - f_10 * gpp1_134[k]
                   + f_4 * pc_x[k] * gpd_266[k];

        t_443[k] = f_4 * pc_x[k] * gpd_267[k];

        t_444[k] = f_4 * pc_x[k] * gpd_268[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_y, gpp0_133, gpp0_134, gpp1_133, \
                         gpp1_134, gpd_267, gpd_268, gpd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_4 * pc_x[k] * gpd_269[k];

        t_446[k] = f_2 * gpp0_133[k]
                   - f_3 * gpp1_133[k]
                   + f_4 * pc_y[k] * gpd_267[k];

        t_447[k] = f_9 * gpp0_134[k]
                   - f_10 * gpp1_134[k]
                   + f_4 * pc_y[k] * gpd_268[k];

        t_448[k] = f_4 * pc_y[k] * gpd_269[k];
    }

#pragma omp simd aligned(t_449, pc_z, fpd_179, gsd_89, gpp0_134, gpp1_134, \
                         gpd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * fpd_179[k]
                   + f_1 * gsd_89[k]
                   + f_2 * gpp0_134[k]
                   - f_3 * gpp1_134[k]
                   + f_4 * pc_z[k] * gpd_269[k];
    }
}

auto
compute_prim_gpf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dpf0,
                                                   const size_t dpf1, const size_t fpf0,
                                                   const size_t fpd, const size_t fpf1,
                                                   const size_t gsf0, const size_t gsd,
                                                   const size_t gsf1, const size_t gpp0,
                                                   const size_t gpp1, const size_t gpd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gpf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, dpf0,
                                                              dpf1, fpf0, fpd, fpf1, gsf0, gsd,
                                                              gsf1, gpp0, gpp1, gpd, ncols,
                                                              gamma, p, q);

    compute_prim_gpf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, dpf0,
                                                              dpf1, fpf0, fpd, fpf1, gsf0, gsd,
                                                              gsf1, gpp0, gpp1, gpd, ncols,
                                                              gamma, p, q);

    compute_prim_gpf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, dpf0,
                                                              dpf1, fpf0, fpd, fpf1, gsf0, gsd,
                                                              gsf1, gpp0, gpp1, gpd, ncols,
                                                              gamma, p, q);

    compute_prim_gpf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, dpf0,
                                                              dpf1, fpf0, fpd, fpf1, gsf0, gsd,
                                                              gsf1, gpp0, gpp1, gpd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
