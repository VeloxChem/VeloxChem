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


#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fpf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppf0, const size_t ppf1,
                                                          const size_t dpf0, const size_t dpd,
                                                          const size_t dpf1, const size_t fsf0,
                                                          const size_t fsd, const size_t fsf1,
                                                          const size_t fpp0, const size_t fpp1,
                                                          const size_t fpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppf0_46 = buffer.data(ppf0 + 46);
    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppf1_46 = buffer.data(ppf1 + 46);
    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *dpf0_0 = buffer.data(dpf0 + 0);
    const auto *dpf0_3 = buffer.data(dpf0 + 3);
    const auto *dpf0_5 = buffer.data(dpf0 + 5);
    const auto *dpf0_6 = buffer.data(dpf0 + 6);
    const auto *dpf0_9 = buffer.data(dpf0 + 9);
    const auto *dpf0_10 = buffer.data(dpf0 + 10);
    const auto *dpf0_16 = buffer.data(dpf0 + 16);
    const auto *dpf0_20 = buffer.data(dpf0 + 20);
    const auto *dpf0_29 = buffer.data(dpf0 + 29);
    const auto *dpf0_33 = buffer.data(dpf0 + 33);
    const auto *dpf0_46 = buffer.data(dpf0 + 46);
    const auto *dpf0_60 = buffer.data(dpf0 + 60);
    const auto *dpf0_65 = buffer.data(dpf0 + 65);
    const auto *dpf0_89 = buffer.data(dpf0 + 89);
    const auto *dpf0_100 = buffer.data(dpf0 + 100);
    const auto *dpf0_106 = buffer.data(dpf0 + 106);
    const auto *dpf0_108 = buffer.data(dpf0 + 108);
    const auto *dpf0_109 = buffer.data(dpf0 + 109);
    const auto *dpf0_116 = buffer.data(dpf0 + 116);
    const auto *dpf0_119 = buffer.data(dpf0 + 119);

    const auto *dpd_0 = buffer.data(dpd + 0);
    const auto *dpd_3 = buffer.data(dpd + 3);
    const auto *dpd_5 = buffer.data(dpd + 5);
    const auto *dpd_6 = buffer.data(dpd + 6);
    const auto *dpd_9 = buffer.data(dpd + 9);
    const auto *dpd_11 = buffer.data(dpd + 11);
    const auto *dpd_12 = buffer.data(dpd + 12);
    const auto *dpd_15 = buffer.data(dpd + 15);
    const auto *dpd_17 = buffer.data(dpd + 17);
    const auto *dpd_18 = buffer.data(dpd + 18);
    const auto *dpd_21 = buffer.data(dpd + 21);
    const auto *dpd_23 = buffer.data(dpd + 23);
    const auto *dpd_24 = buffer.data(dpd + 24);
    const auto *dpd_27 = buffer.data(dpd + 27);
    const auto *dpd_29 = buffer.data(dpd + 29);
    const auto *dpd_30 = buffer.data(dpd + 30);
    const auto *dpd_33 = buffer.data(dpd + 33);
    const auto *dpd_35 = buffer.data(dpd + 35);
    const auto *dpd_36 = buffer.data(dpd + 36);
    const auto *dpd_39 = buffer.data(dpd + 39);
    const auto *dpd_41 = buffer.data(dpd + 41);
    const auto *dpd_45 = buffer.data(dpd + 45);
    const auto *dpd_47 = buffer.data(dpd + 47);
    const auto *dpd_48 = buffer.data(dpd + 48);
    const auto *dpd_51 = buffer.data(dpd + 51);
    const auto *dpd_53 = buffer.data(dpd + 53);
    const auto *dpd_54 = buffer.data(dpd + 54);
    const auto *dpd_57 = buffer.data(dpd + 57);
    const auto *dpd_59 = buffer.data(dpd + 59);
    const auto *dpd_60 = buffer.data(dpd + 60);
    const auto *dpd_63 = buffer.data(dpd + 63);
    const auto *dpd_65 = buffer.data(dpd + 65);
    const auto *dpd_69 = buffer.data(dpd + 69);
    const auto *dpd_71 = buffer.data(dpd + 71);
    const auto *dpd_76 = buffer.data(dpd + 76);

    const auto *dpf1_0 = buffer.data(dpf1 + 0);
    const auto *dpf1_3 = buffer.data(dpf1 + 3);
    const auto *dpf1_5 = buffer.data(dpf1 + 5);
    const auto *dpf1_6 = buffer.data(dpf1 + 6);
    const auto *dpf1_9 = buffer.data(dpf1 + 9);
    const auto *dpf1_10 = buffer.data(dpf1 + 10);
    const auto *dpf1_16 = buffer.data(dpf1 + 16);
    const auto *dpf1_20 = buffer.data(dpf1 + 20);
    const auto *dpf1_29 = buffer.data(dpf1 + 29);
    const auto *dpf1_33 = buffer.data(dpf1 + 33);
    const auto *dpf1_46 = buffer.data(dpf1 + 46);
    const auto *dpf1_60 = buffer.data(dpf1 + 60);
    const auto *dpf1_65 = buffer.data(dpf1 + 65);
    const auto *dpf1_89 = buffer.data(dpf1 + 89);
    const auto *dpf1_100 = buffer.data(dpf1 + 100);
    const auto *dpf1_106 = buffer.data(dpf1 + 106);
    const auto *dpf1_108 = buffer.data(dpf1 + 108);
    const auto *dpf1_109 = buffer.data(dpf1 + 109);
    const auto *dpf1_116 = buffer.data(dpf1 + 116);
    const auto *dpf1_119 = buffer.data(dpf1 + 119);

    const auto *fsf0_0 = buffer.data(fsf0 + 0);
    const auto *fsf0_6 = buffer.data(fsf0 + 6);
    const auto *fsf0_9 = buffer.data(fsf0 + 9);
    const auto *fsf0_16 = buffer.data(fsf0 + 16);
    const auto *fsf0_27 = buffer.data(fsf0 + 27);
    const auto *fsf0_29 = buffer.data(fsf0 + 29);
    const auto *fsf0_30 = buffer.data(fsf0 + 30);

    const auto *fsd_0 = buffer.data(fsd + 0);
    const auto *fsd_2 = buffer.data(fsd + 2);
    const auto *fsd_3 = buffer.data(fsd + 3);
    const auto *fsd_5 = buffer.data(fsd + 5);
    const auto *fsd_6 = buffer.data(fsd + 6);
    const auto *fsd_7 = buffer.data(fsd + 7);
    const auto *fsd_9 = buffer.data(fsd + 9);
    const auto *fsd_11 = buffer.data(fsd + 11);
    const auto *fsd_12 = buffer.data(fsd + 12);
    const auto *fsd_14 = buffer.data(fsd + 14);
    const auto *fsd_16 = buffer.data(fsd + 16);
    const auto *fsd_17 = buffer.data(fsd + 17);
    const auto *fsd_18 = buffer.data(fsd + 18);
    const auto *fsd_19 = buffer.data(fsd + 19);
    const auto *fsd_21 = buffer.data(fsd + 21);
    const auto *fsd_23 = buffer.data(fsd + 23);
    const auto *fsd_28 = buffer.data(fsd + 28);

    const auto *fsf1_0 = buffer.data(fsf1 + 0);
    const auto *fsf1_6 = buffer.data(fsf1 + 6);
    const auto *fsf1_9 = buffer.data(fsf1 + 9);
    const auto *fsf1_16 = buffer.data(fsf1 + 16);
    const auto *fsf1_27 = buffer.data(fsf1 + 27);
    const auto *fsf1_29 = buffer.data(fsf1 + 29);
    const auto *fsf1_30 = buffer.data(fsf1 + 30);

    const auto *fpp0_0 = buffer.data(fpp0 + 0);
    const auto *fpp0_1 = buffer.data(fpp0 + 1);
    const auto *fpp0_2 = buffer.data(fpp0 + 2);
    const auto *fpp0_10 = buffer.data(fpp0 + 10);
    const auto *fpp0_12 = buffer.data(fpp0 + 12);
    const auto *fpp0_14 = buffer.data(fpp0 + 14);
    const auto *fpp0_20 = buffer.data(fpp0 + 20);
    const auto *fpp0_24 = buffer.data(fpp0 + 24);
    const auto *fpp0_25 = buffer.data(fpp0 + 25);
    const auto *fpp0_26 = buffer.data(fpp0 + 26);
    const auto *fpp0_27 = buffer.data(fpp0 + 27);
    const auto *fpp0_28 = buffer.data(fpp0 + 28);
    const auto *fpp0_29 = buffer.data(fpp0 + 29);
    const auto *fpp0_37 = buffer.data(fpp0 + 37);
    const auto *fpp0_38 = buffer.data(fpp0 + 38);

    const auto *fpp1_0 = buffer.data(fpp1 + 0);
    const auto *fpp1_1 = buffer.data(fpp1 + 1);
    const auto *fpp1_2 = buffer.data(fpp1 + 2);
    const auto *fpp1_10 = buffer.data(fpp1 + 10);
    const auto *fpp1_12 = buffer.data(fpp1 + 12);
    const auto *fpp1_14 = buffer.data(fpp1 + 14);
    const auto *fpp1_20 = buffer.data(fpp1 + 20);
    const auto *fpp1_24 = buffer.data(fpp1 + 24);
    const auto *fpp1_25 = buffer.data(fpp1 + 25);
    const auto *fpp1_26 = buffer.data(fpp1 + 26);
    const auto *fpp1_27 = buffer.data(fpp1 + 27);
    const auto *fpp1_28 = buffer.data(fpp1 + 28);
    const auto *fpp1_29 = buffer.data(fpp1 + 29);
    const auto *fpp1_37 = buffer.data(fpp1 + 37);
    const auto *fpp1_38 = buffer.data(fpp1 + 38);

    const auto *fpd_0 = buffer.data(fpd + 0);
    const auto *fpd_2 = buffer.data(fpd + 2);
    const auto *fpd_3 = buffer.data(fpd + 3);
    const auto *fpd_5 = buffer.data(fpd + 5);
    const auto *fpd_6 = buffer.data(fpd + 6);
    const auto *fpd_8 = buffer.data(fpd + 8);
    const auto *fpd_9 = buffer.data(fpd + 9);
    const auto *fpd_11 = buffer.data(fpd + 11);
    const auto *fpd_12 = buffer.data(fpd + 12);
    const auto *fpd_14 = buffer.data(fpd + 14);
    const auto *fpd_15 = buffer.data(fpd + 15);
    const auto *fpd_17 = buffer.data(fpd + 17);
    const auto *fpd_18 = buffer.data(fpd + 18);
    const auto *fpd_19 = buffer.data(fpd + 19);
    const auto *fpd_21 = buffer.data(fpd + 21);
    const auto *fpd_23 = buffer.data(fpd + 23);
    const auto *fpd_24 = buffer.data(fpd + 24);
    const auto *fpd_25 = buffer.data(fpd + 25);
    const auto *fpd_27 = buffer.data(fpd + 27);
    const auto *fpd_29 = buffer.data(fpd + 29);
    const auto *fpd_30 = buffer.data(fpd + 30);
    const auto *fpd_31 = buffer.data(fpd + 31);
    const auto *fpd_33 = buffer.data(fpd + 33);
    const auto *fpd_35 = buffer.data(fpd + 35);
    const auto *fpd_36 = buffer.data(fpd + 36);
    const auto *fpd_38 = buffer.data(fpd + 38);
    const auto *fpd_40 = buffer.data(fpd + 40);
    const auto *fpd_41 = buffer.data(fpd + 41);
    const auto *fpd_42 = buffer.data(fpd + 42);
    const auto *fpd_44 = buffer.data(fpd + 44);
    const auto *fpd_45 = buffer.data(fpd + 45);
    const auto *fpd_47 = buffer.data(fpd + 47);
    const auto *fpd_48 = buffer.data(fpd + 48);
    const auto *fpd_50 = buffer.data(fpd + 50);
    const auto *fpd_51 = buffer.data(fpd + 51);
    const auto *fpd_52 = buffer.data(fpd + 52);
    const auto *fpd_53 = buffer.data(fpd + 53);
    const auto *fpd_54 = buffer.data(fpd + 54);
    const auto *fpd_55 = buffer.data(fpd + 55);
    const auto *fpd_57 = buffer.data(fpd + 57);
    const auto *fpd_59 = buffer.data(fpd + 59);
    const auto *fpd_60 = buffer.data(fpd + 60);
    const auto *fpd_61 = buffer.data(fpd + 61);
    const auto *fpd_63 = buffer.data(fpd + 63);
    const auto *fpd_65 = buffer.data(fpd + 65);
    const auto *fpd_66 = buffer.data(fpd + 66);
    const auto *fpd_67 = buffer.data(fpd + 67);
    const auto *fpd_69 = buffer.data(fpd + 69);
    const auto *fpd_71 = buffer.data(fpd + 71);
    const auto *fpd_72 = buffer.data(fpd + 72);
    const auto *fpd_75 = buffer.data(fpd + 75);
    const auto *fpd_76 = buffer.data(fpd + 76);
    const auto *fpd_77 = buffer.data(fpd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, dpd_0, dpd_3, fsd_0, fsd_3, \
                         fpp0_0, fpp1_0, fpd_0, fpd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dpd_0[k]
                 + f_1 * fsd_0[k]
                 + f_2 * fpp0_0[k]
                 - f_3 * fpp1_0[k]
                 + f_4 * pc_x[k] * fpd_0[k];

        t_1[k] = f_4 * pc_y[k] * fpd_0[k];

        t_2[k] = f_4 * pc_z[k] * fpd_0[k];

        t_3[k] = f_0 * dpd_3[k]
                 + f_1 * fsd_3[k]
                 + f_4 * pc_x[k] * fpd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, dpd_5, fsd_5, fpp0_1, \
                         fpp1_1, fpd_2, fpd_3, fpd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * pc_y[k] * fpd_2[k];

        t_5[k] = f_0 * dpd_5[k]
                 + f_1 * fsd_5[k]
                 + f_4 * pc_x[k] * fpd_5[k];

        t_6[k] = f_2 * fpp0_1[k]
                 - f_3 * fpp1_1[k]
                 + f_4 * pc_y[k] * fpd_3[k];

        t_7[k] = f_4 * pc_z[k] * fpd_3[k];

        t_8[k] = f_4 * pc_y[k] * fpd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, fsf0_0, fsd_0, fsf1_0, \
                         fpp0_2, fpp1_2, fpd_5, fpd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * fpp0_2[k]
                 - f_3 * fpp1_2[k]
                 + f_4 * pc_z[k] * fpd_5[k];

        t_10[k] = pb_y[k] * fsf0_0[k]
                  - f_5 * pc_y[k] * fsf1_0[k];

        t_11[k] = f_1 * fsd_0[k]
                  + f_4 * pc_y[k] * fpd_6[k];

        t_12[k] = f_4 * pc_z[k] * fpd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pc_x, pc_y, dpd_9, dpd_11, fsf0_6, \
                         fsd_2, fsd_3, fsf1_6, fpd_8, fpd_9, fpd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dpd_9[k]
                  + f_4 * pc_x[k] * fpd_9[k];

        t_14[k] = f_1 * fsd_2[k]
                  + f_4 * pc_y[k] * fpd_8[k];

        t_15[k] = f_0 * dpd_11[k]
                  + f_4 * pc_x[k] * fpd_11[k];

        t_16[k] = pb_y[k] * fsf0_6[k]
                  + f_0 * fsd_3[k]
                  - f_5 * pc_y[k] * fsf1_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, pc_y, pc_z, fsf0_0, fsf0_9, \
                         fsd_5, fsf1_0, fsf1_9, fpd_9, fpd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * pc_z[k] * fpd_9[k];

        t_18[k] = f_1 * fsd_5[k]
                  + f_4 * pc_y[k] * fpd_11[k];

        t_19[k] = pb_y[k] * fsf0_9[k]
                  - f_5 * pc_y[k] * fsf1_9[k];

        t_20[k] = pb_z[k] * fsf0_0[k]
                  - f_5 * pc_z[k] * fsf1_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, dpd_15, dpd_17, \
                         fsd_0, fpd_12, fpd_14, fpd_15, fpd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * pc_y[k] * fpd_12[k];

        t_22[k] = f_1 * fsd_0[k]
                  + f_4 * pc_z[k] * fpd_12[k];

        t_23[k] = f_0 * dpd_15[k]
                  + f_4 * pc_x[k] * fpd_15[k];

        t_24[k] = f_4 * pc_y[k] * fpd_14[k];

        t_25[k] = f_0 * dpd_17[k]
                  + f_4 * pc_x[k] * fpd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, fsf0_6, fsf0_9, fsd_3, \
                         fsd_5, fsf1_6, fsf1_9, fpd_15, fpd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * fsf0_6[k]
                  - f_5 * pc_z[k] * fsf1_6[k];

        t_27[k] = f_1 * fsd_3[k]
                  + f_4 * pc_z[k] * fpd_15[k];

        t_28[k] = f_4 * pc_y[k] * fpd_17[k];

        t_29[k] = pb_z[k] * fsf0_9[k]
                  + f_0 * fsd_5[k]
                  - f_5 * pc_z[k] * fsf1_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_x, pc_y, pc_z, dpf0_0, dpd_0, \
                         dpd_21, dpf1_0, fsd_9, fpd_18, fpd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * dpf0_0[k]
                  - f_5 * pc_y[k] * dpf1_0[k];

        t_31[k] = f_1 * dpd_0[k]
                  + f_4 * pc_y[k] * fpd_18[k];

        t_32[k] = f_4 * pc_z[k] * fpd_18[k];

        t_33[k] = f_6 * dpd_21[k]
                  + f_1 * fsd_9[k]
                  + f_4 * pc_x[k] * fpd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, dpf0_5, dpd_3, dpf1_5, \
                         fpp0_10, fpp1_10, fpd_19, fpd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_4 * pc_z[k] * fpd_19[k];

        t_35[k] = pa_y[k] * dpf0_5[k]
                  - f_5 * pc_y[k] * dpf1_5[k];

        t_36[k] = f_1 * dpd_3[k]
                  + f_2 * fpp0_10[k]
                  - f_3 * fpp1_10[k]
                  + f_4 * pc_y[k] * fpd_21[k];

        t_37[k] = f_4 * pc_z[k] * fpd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pc_x, pc_y, dpf0_9, dpd_5, dpd_24, dpf1_9, \
                         fpp0_12, fpp1_12, fpd_23, fpd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * dpd_5[k]
                  + f_4 * pc_y[k] * fpd_23[k];

        t_39[k] = pa_y[k] * dpf0_9[k]
                  - f_5 * pc_y[k] * dpf1_9[k];

        t_40[k] = f_6 * dpd_24[k]
                  + f_2 * fpp0_12[k]
                  - f_3 * fpp1_12[k]
                  + f_4 * pc_x[k] * fpd_24[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pc_x, pc_y, pc_z, dpd_6, dpd_27, \
                         dpd_29, fsd_6, fpd_24, fpd_25, fpd_27, \
                         fpd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * dpd_6[k]
                  + f_1 * fsd_6[k]
                  + f_4 * pc_y[k] * fpd_24[k];

        t_42[k] = f_4 * pc_z[k] * fpd_24[k];

        t_43[k] = f_6 * dpd_27[k]
                  + f_4 * pc_x[k] * fpd_27[k];

        t_44[k] = f_4 * pc_z[k] * fpd_25[k];

        t_45[k] = f_6 * dpd_29[k]
                  + f_4 * pc_x[k] * fpd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pc_x, pc_y, pc_z, ppf0_46, ppf1_46, dpf0_46, \
                         dpd_11, dpf1_46, fsd_11, fpd_27, fpd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_7 * ppf0_46[k]
                  - f_8 * ppf1_46[k]
                  + pa_x[k] * dpf0_46[k]
                  - f_5 * pc_x[k] * dpf1_46[k];

        t_47[k] = f_4 * pc_z[k] * fpd_27[k];

        t_48[k] = f_1 * dpd_11[k]
                  + f_1 * fsd_11[k]
                  + f_4 * pc_y[k] * fpd_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pc_y, pc_z, dpf0_20, dpd_12, dpf1_20, \
                         fsd_6, fpp0_14, fpp1_14, fpd_29, fpd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * fpp0_14[k]
                  - f_3 * fpp1_14[k]
                  + f_4 * pc_z[k] * fpd_29[k];

        t_50[k] = pa_y[k] * dpf0_20[k]
                  - f_5 * pc_y[k] * dpf1_20[k];

        t_51[k] = f_1 * dpd_12[k]
                  + f_4 * pc_y[k] * fpd_30[k];

        t_52[k] = f_1 * fsd_6[k]
                  + f_4 * pc_z[k] * fpd_30[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_z, pc_x, pc_z, dpd_33, dpd_35, fsf0_16, \
                         fsd_7, fsf1_16, fpd_31, fpd_33, fpd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_6 * dpd_33[k]
                  + f_4 * pc_x[k] * fpd_33[k];

        t_54[k] = f_1 * fsd_7[k]
                  + f_4 * pc_z[k] * fpd_31[k];

        t_55[k] = f_6 * dpd_35[k]
                  + f_4 * pc_x[k] * fpd_35[k];

        t_56[k] = pb_z[k] * fsf0_16[k]
                  - f_5 * pc_z[k] * fsf1_16[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_y, pa_z, pc_y, pc_z, dpf0_0, dpf0_29, \
                         dpd_17, dpf1_0, dpf1_29, fsd_9, fpd_33, \
                         fpd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * fsd_9[k]
                  + f_4 * pc_z[k] * fpd_33[k];

        t_58[k] = f_1 * dpd_17[k]
                  + f_4 * pc_y[k] * fpd_35[k];

        t_59[k] = pa_y[k] * dpf0_29[k]
                  - f_5 * pc_y[k] * dpf1_29[k];

        t_60[k] = pa_z[k] * dpf0_0[k]
                  - f_5 * pc_z[k] * dpf1_0[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_z, pc_y, pc_z, dpf0_3, dpd_0, dpf1_3, \
                         fpd_36, fpd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_4 * pc_y[k] * fpd_36[k];

        t_62[k] = f_1 * dpd_0[k]
                  + f_4 * pc_z[k] * fpd_36[k];

        t_63[k] = pa_z[k] * dpf0_3[k]
                  - f_5 * pc_z[k] * dpf1_3[k];

        t_64[k] = f_4 * pc_y[k] * fpd_38[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pc_x, pc_y, pc_z, dpf0_6, dpd_41, \
                         dpf1_6, fsd_17, fpp0_20, fpp1_20, fpd_40, \
                         fpd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_6 * dpd_41[k]
                  + f_1 * fsd_17[k]
                  + f_4 * pc_x[k] * fpd_41[k];

        t_66[k] = pa_z[k] * dpf0_6[k]
                  - f_5 * pc_z[k] * dpf1_6[k];

        t_67[k] = f_9 * fpp0_20[k]
                  - f_10 * fpp1_20[k]
                  + f_4 * pc_y[k] * fpd_40[k];

        t_68[k] = f_4 * pc_y[k] * fpd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, dpf0_10, dpd_5, dpd_6, \
                         dpf1_10, fsd_12, fpp0_20, fpp1_20, fpd_41, \
                         fpd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * dpd_5[k]
                  + f_2 * fpp0_20[k]
                  - f_3 * fpp1_20[k]
                  + f_4 * pc_z[k] * fpd_41[k];

        t_70[k] = pa_z[k] * dpf0_10[k]
                  - f_5 * pc_z[k] * dpf1_10[k];

        t_71[k] = f_1 * fsd_12[k]
                  + f_4 * pc_y[k] * fpd_42[k];

        t_72[k] = f_1 * dpd_6[k]
                  + f_4 * pc_z[k] * fpd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_y, pc_z, dpf0_16, dpd_45, \
                         dpd_47, dpf1_16, fsd_14, fpd_44, fpd_45, \
                         fpd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * dpd_45[k]
                  + f_4 * pc_x[k] * fpd_45[k];

        t_74[k] = f_1 * fsd_14[k]
                  + f_4 * pc_y[k] * fpd_44[k];

        t_75[k] = f_6 * dpd_47[k]
                  + f_4 * pc_x[k] * fpd_47[k];

        t_76[k] = pa_z[k] * dpf0_16[k]
                  - f_5 * pc_z[k] * dpf1_16[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, pc_y, fsf0_27, fsf0_29, fsd_16, fsd_17, \
                         fsf1_27, fsf1_29, fpd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_y[k] * fsf0_27[k]
                  + f_6 * fsd_16[k]
                  - f_5 * pc_y[k] * fsf1_27[k];

        t_78[k] = f_1 * fsd_17[k]
                  + f_4 * pc_y[k] * fpd_47[k];

        t_79[k] = pb_y[k] * fsf0_29[k]
                  - f_5 * pc_y[k] * fsf1_29[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, dpd_12, dpd_48, dpd_51, \
                         fsd_12, fpp0_24, fpp1_24, fpd_48, fpd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * dpd_48[k]
                  + f_2 * fpp0_24[k]
                  - f_3 * fpp1_24[k]
                  + f_4 * pc_x[k] * fpd_48[k];

        t_81[k] = f_4 * pc_y[k] * fpd_48[k];

        t_82[k] = f_1 * dpd_12[k]
                  + f_1 * fsd_12[k]
                  + f_4 * pc_z[k] * fpd_48[k];

        t_83[k] = f_6 * dpd_51[k]
                  + f_4 * pc_x[k] * fpd_51[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pc_x, pc_y, dpd_53, fpp0_25, fpp0_26, \
                         fpp1_25, fpp1_26, fpd_50, fpd_51, fpd_52, \
                         fpd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * pc_y[k] * fpd_50[k];

        t_85[k] = f_6 * dpd_53[k]
                  + f_4 * pc_x[k] * fpd_53[k];

        t_86[k] = f_2 * fpp0_25[k]
                  - f_3 * fpp1_25[k]
                  + f_4 * pc_y[k] * fpd_51[k];

        t_87[k] = f_9 * fpp0_26[k]
                  - f_10 * fpp1_26[k]
                  + f_4 * pc_y[k] * fpd_52[k];

        t_88[k] = f_4 * pc_y[k] * fpd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pc_x, pc_y, ppf0_89, ppf1_89, dpf0_89, \
                         dpd_18, dpd_54, dpf1_89, fsd_18, fpp0_27, fpp1_27, \
                         fpd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_7 * ppf0_89[k]
                  - f_8 * ppf1_89[k]
                  + pa_x[k] * dpf0_89[k]
                  - f_5 * pc_x[k] * dpf1_89[k];

        t_90[k] = f_1 * dpd_54[k]
                  + f_1 * fsd_18[k]
                  + f_2 * fpp0_27[k]
                  - f_3 * fpp1_27[k]
                  + f_4 * pc_x[k] * fpd_54[k];

        t_91[k] = f_6 * dpd_18[k]
                  + f_4 * pc_y[k] * fpd_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_z, dpd_57, dpd_59, fsd_21, fsd_23, \
                         fpd_54, fpd_55, fpd_57, fpd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * pc_z[k] * fpd_54[k];

        t_93[k] = f_1 * dpd_57[k]
                  + f_1 * fsd_21[k]
                  + f_4 * pc_x[k] * fpd_57[k];

        t_94[k] = f_4 * pc_z[k] * fpd_55[k];

        t_95[k] = f_1 * dpd_59[k]
                  + f_1 * fsd_23[k]
                  + f_4 * pc_x[k] * fpd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_y, pc_z, dpd_21, dpd_23, fpp0_28, fpp0_29, \
                         fpp1_28, fpp1_29, fpd_57, fpd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_6 * dpd_21[k]
                  + f_2 * fpp0_28[k]
                  - f_3 * fpp1_28[k]
                  + f_4 * pc_y[k] * fpd_57[k];

        t_97[k] = f_4 * pc_z[k] * fpd_57[k];

        t_98[k] = f_6 * dpd_23[k]
                  + f_4 * pc_y[k] * fpd_59[k];

        t_99[k] = f_2 * fpp0_29[k]
                  - f_3 * fpp1_29[k]
                  + f_4 * pc_z[k] * fpd_59[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pc_x, pc_y, pc_z, dpf0_100, dpd_24, \
                         dpd_60, dpd_63, dpf1_100, fsd_18, fpd_60, \
                         fpd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * dpf0_100[k]
                   + f_0 * dpd_60[k]
                   - f_5 * pc_x[k] * dpf1_100[k];

        t_101[k] = f_6 * dpd_24[k]
                   + f_1 * fsd_18[k]
                   + f_4 * pc_y[k] * fpd_60[k];

        t_102[k] = f_4 * pc_z[k] * fpd_60[k];

        t_103[k] = f_1 * dpd_63[k]
                   + f_4 * pc_x[k] * fpd_63[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_x, pc_x, pc_z, dpf0_106, \
                         dpf0_108, dpd_65, dpf1_106, dpf1_108, fpd_61, fpd_63, \
                         fpd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_4 * pc_z[k] * fpd_61[k];

        t_105[k] = f_1 * dpd_65[k]
                   + f_4 * pc_x[k] * fpd_65[k];

        t_106[k] = pa_x[k] * dpf0_106[k]
                   - f_5 * pc_x[k] * dpf1_106[k];

        t_107[k] = f_4 * pc_z[k] * fpd_63[k];

        t_108[k] = pa_x[k] * dpf0_108[k]
                   - f_5 * pc_x[k] * dpf1_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pb_z, pc_x, pc_y, pc_z, dpf0_109, \
                         dpd_30, dpf1_109, fsf0_30, fsd_18, fsf1_30, \
                         fpd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * dpf0_109[k]
                   - f_5 * pc_x[k] * dpf1_109[k];

        t_110[k] = pb_z[k] * fsf0_30[k]
                   - f_5 * pc_z[k] * fsf1_30[k];

        t_111[k] = f_6 * dpd_30[k]
                   + f_4 * pc_y[k] * fpd_66[k];

        t_112[k] = f_1 * fsd_18[k]
                   + f_4 * pc_z[k] * fpd_66[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pc_x, pc_z, dpf0_116, dpd_69, \
                         dpd_71, dpf1_116, fsd_19, fpd_67, fpd_69, \
                         fpd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_1 * dpd_69[k]
                   + f_4 * pc_x[k] * fpd_69[k];

        t_114[k] = f_1 * fsd_19[k]
                   + f_4 * pc_z[k] * fpd_67[k];

        t_115[k] = f_1 * dpd_71[k]
                   + f_4 * pc_x[k] * fpd_71[k];

        t_116[k] = pa_x[k] * dpf0_116[k]
                   - f_5 * pc_x[k] * dpf1_116[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_x, pc_x, pc_y, pc_z, dpf0_119, dpd_35, \
                         dpf1_119, fsd_21, fpd_69, fpd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_1 * fsd_21[k]
                   + f_4 * pc_z[k] * fpd_69[k];

        t_118[k] = f_6 * dpd_35[k]
                   + f_4 * pc_y[k] * fpd_71[k];

        t_119[k] = pa_x[k] * dpf0_119[k]
                   - f_5 * pc_x[k] * dpf1_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, dpf0_33, dpf0_60, \
                         dpd_18, dpd_36, dpf1_33, dpf1_60, fpd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * dpf0_60[k]
                   - f_5 * pc_y[k] * dpf1_60[k];

        t_121[k] = f_1 * dpd_36[k]
                   + f_4 * pc_y[k] * fpd_72[k];

        t_122[k] = f_1 * dpd_18[k]
                   + f_4 * pc_z[k] * fpd_72[k];

        t_123[k] = pa_z[k] * dpf0_33[k]
                   - f_5 * pc_z[k] * dpf1_33[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pc_x, pc_y, dpf0_65, dpd_39, dpd_76, \
                         dpf1_65, fsd_28, fpp0_37, fpp1_37, fpd_75, \
                         fpd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * dpd_76[k]
                   + f_1 * fsd_28[k]
                   + f_4 * pc_x[k] * fpd_76[k];

        t_125[k] = pa_y[k] * dpf0_65[k]
                   - f_5 * pc_y[k] * dpf1_65[k];

        t_126[k] = f_1 * dpd_39[k]
                   + f_2 * fpp0_37[k]
                   - f_3 * fpp1_37[k]
                   + f_4 * pc_y[k] * fpd_75[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, dpd_21, dpd_23, dpd_41, fpp0_38, \
                         fpp1_38, fpd_75, fpd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_1 * dpd_21[k]
                   + f_4 * pc_z[k] * fpd_75[k];

        t_128[k] = f_1 * dpd_41[k]
                   + f_4 * pc_y[k] * fpd_77[k];

        t_129[k] = f_1 * dpd_23[k]
                   + f_2 * fpp0_38[k]
                   - f_3 * fpp1_38[k]
                   + f_4 * pc_z[k] * fpd_77[k];
    }
}

static auto
compute_prim_fpf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppf0, const size_t ppf1,
                                                          const size_t dpf0, const size_t dpd,
                                                          const size_t dpf1, const size_t fsf0,
                                                          const size_t fsd, const size_t fsf1,
                                                          const size_t fpp0, const size_t fpp1,
                                                          const size_t fpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);

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
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *dpf0_41 = buffer.data(dpf0 + 41);
    const auto *dpf0_80 = buffer.data(dpf0 + 80);
    const auto *dpf0_82 = buffer.data(dpf0 + 82);
    const auto *dpf0_90 = buffer.data(dpf0 + 90);
    const auto *dpf0_91 = buffer.data(dpf0 + 91);
    const auto *dpf0_96 = buffer.data(dpf0 + 96);
    const auto *dpf0_100 = buffer.data(dpf0 + 100);
    const auto *dpf0_101 = buffer.data(dpf0 + 101);
    const auto *dpf0_106 = buffer.data(dpf0 + 106);
    const auto *dpf0_136 = buffer.data(dpf0 + 136);
    const auto *dpf0_138 = buffer.data(dpf0 + 138);
    const auto *dpf0_139 = buffer.data(dpf0 + 139);
    const auto *dpf0_146 = buffer.data(dpf0 + 146);
    const auto *dpf0_147 = buffer.data(dpf0 + 147);
    const auto *dpf0_149 = buffer.data(dpf0 + 149);
    const auto *dpf0_150 = buffer.data(dpf0 + 150);
    const auto *dpf0_151 = buffer.data(dpf0 + 151);
    const auto *dpf0_152 = buffer.data(dpf0 + 152);
    const auto *dpf0_159 = buffer.data(dpf0 + 159);
    const auto *dpf0_166 = buffer.data(dpf0 + 166);
    const auto *dpf0_167 = buffer.data(dpf0 + 167);
    const auto *dpf0_169 = buffer.data(dpf0 + 169);
    const auto *dpf0_170 = buffer.data(dpf0 + 170);
    const auto *dpf0_176 = buffer.data(dpf0 + 176);
    const auto *dpf0_177 = buffer.data(dpf0 + 177);
    const auto *dpf0_179 = buffer.data(dpf0 + 179);

    const auto *dpd_24 = buffer.data(dpd + 24);
    const auto *dpd_27 = buffer.data(dpd + 27);
    const auto *dpd_36 = buffer.data(dpd + 36);
    const auto *dpd_41 = buffer.data(dpd + 41);
    const auto *dpd_42 = buffer.data(dpd + 42);
    const auto *dpd_48 = buffer.data(dpd + 48);
    const auto *dpd_53 = buffer.data(dpd + 53);
    const auto *dpd_57 = buffer.data(dpd + 57);
    const auto *dpd_59 = buffer.data(dpd + 59);
    const auto *dpd_63 = buffer.data(dpd + 63);
    const auto *dpd_65 = buffer.data(dpd + 65);
    const auto *dpd_69 = buffer.data(dpd + 69);
    const auto *dpd_71 = buffer.data(dpd + 71);
    const auto *dpd_75 = buffer.data(dpd + 75);
    const auto *dpd_77 = buffer.data(dpd + 77);
    const auto *dpd_78 = buffer.data(dpd + 78);
    const auto *dpd_81 = buffer.data(dpd + 81);
    const auto *dpd_82 = buffer.data(dpd + 82);
    const auto *dpd_83 = buffer.data(dpd + 83);
    const auto *dpd_87 = buffer.data(dpd + 87);
    const auto *dpd_88 = buffer.data(dpd + 88);
    const auto *dpd_89 = buffer.data(dpd + 89);
    const auto *dpd_90 = buffer.data(dpd + 90);
    const auto *dpd_93 = buffer.data(dpd + 93);
    const auto *dpd_95 = buffer.data(dpd + 95);
    const auto *dpd_99 = buffer.data(dpd + 99);
    const auto *dpd_101 = buffer.data(dpd + 101);
    const auto *dpd_102 = buffer.data(dpd + 102);
    const auto *dpd_105 = buffer.data(dpd + 105);
    const auto *dpd_107 = buffer.data(dpd + 107);

    const auto *dpf1_41 = buffer.data(dpf1 + 41);
    const auto *dpf1_80 = buffer.data(dpf1 + 80);
    const auto *dpf1_82 = buffer.data(dpf1 + 82);
    const auto *dpf1_90 = buffer.data(dpf1 + 90);
    const auto *dpf1_91 = buffer.data(dpf1 + 91);
    const auto *dpf1_96 = buffer.data(dpf1 + 96);
    const auto *dpf1_100 = buffer.data(dpf1 + 100);
    const auto *dpf1_101 = buffer.data(dpf1 + 101);
    const auto *dpf1_106 = buffer.data(dpf1 + 106);
    const auto *dpf1_136 = buffer.data(dpf1 + 136);
    const auto *dpf1_138 = buffer.data(dpf1 + 138);
    const auto *dpf1_139 = buffer.data(dpf1 + 139);
    const auto *dpf1_146 = buffer.data(dpf1 + 146);
    const auto *dpf1_147 = buffer.data(dpf1 + 147);
    const auto *dpf1_149 = buffer.data(dpf1 + 149);
    const auto *dpf1_150 = buffer.data(dpf1 + 150);
    const auto *dpf1_151 = buffer.data(dpf1 + 151);
    const auto *dpf1_152 = buffer.data(dpf1 + 152);
    const auto *dpf1_159 = buffer.data(dpf1 + 159);
    const auto *dpf1_166 = buffer.data(dpf1 + 166);
    const auto *dpf1_167 = buffer.data(dpf1 + 167);
    const auto *dpf1_169 = buffer.data(dpf1 + 169);
    const auto *dpf1_170 = buffer.data(dpf1 + 170);
    const auto *dpf1_176 = buffer.data(dpf1 + 176);
    const auto *dpf1_177 = buffer.data(dpf1 + 177);
    const auto *dpf1_179 = buffer.data(dpf1 + 179);

    const auto *fsf0_50 = buffer.data(fsf0 + 50);
    const auto *fsf0_60 = buffer.data(fsf0 + 60);
    const auto *fsf0_61 = buffer.data(fsf0 + 61);
    const auto *fsf0_66 = buffer.data(fsf0 + 66);
    const auto *fsf0_69 = buffer.data(fsf0 + 69);
    const auto *fsf0_72 = buffer.data(fsf0 + 72);
    const auto *fsf0_79 = buffer.data(fsf0 + 79);
    const auto *fsf0_86 = buffer.data(fsf0 + 86);

    const auto *fsd_30 = buffer.data(fsd + 30);
    const auto *fsd_32 = buffer.data(fsd + 32);
    const auto *fsd_33 = buffer.data(fsd + 33);
    const auto *fsd_35 = buffer.data(fsd + 35);
    const auto *fsd_36 = buffer.data(fsd + 36);
    const auto *fsd_37 = buffer.data(fsd + 37);
    const auto *fsd_39 = buffer.data(fsd + 39);
    const auto *fsd_40 = buffer.data(fsd + 40);
    const auto *fsd_41 = buffer.data(fsd + 41);
    const auto *fsd_44 = buffer.data(fsd + 44);
    const auto *fsd_45 = buffer.data(fsd + 45);
    const auto *fsd_46 = buffer.data(fsd + 46);
    const auto *fsd_47 = buffer.data(fsd + 47);
    const auto *fsd_51 = buffer.data(fsd + 51);
    const auto *fsd_52 = buffer.data(fsd + 52);
    const auto *fsd_53 = buffer.data(fsd + 53);

    const auto *fsf1_50 = buffer.data(fsf1 + 50);
    const auto *fsf1_60 = buffer.data(fsf1 + 60);
    const auto *fsf1_61 = buffer.data(fsf1 + 61);
    const auto *fsf1_66 = buffer.data(fsf1 + 66);
    const auto *fsf1_69 = buffer.data(fsf1 + 69);
    const auto *fsf1_72 = buffer.data(fsf1 + 72);
    const auto *fsf1_79 = buffer.data(fsf1 + 79);
    const auto *fsf1_86 = buffer.data(fsf1 + 86);

    const auto *fpp0_39 = buffer.data(fpp0 + 39);
    const auto *fpp0_45 = buffer.data(fpp0 + 45);
    const auto *fpp0_46 = buffer.data(fpp0 + 46);
    const auto *fpp0_47 = buffer.data(fpp0 + 47);
    const auto *fpp0_57 = buffer.data(fpp0 + 57);
    const auto *fpp0_58 = buffer.data(fpp0 + 58);
    const auto *fpp0_59 = buffer.data(fpp0 + 59);
    const auto *fpp0_68 = buffer.data(fpp0 + 68);
    const auto *fpp0_69 = buffer.data(fpp0 + 69);
    const auto *fpp0_70 = buffer.data(fpp0 + 70);
    const auto *fpp0_71 = buffer.data(fpp0 + 71);
    const auto *fpp0_75 = buffer.data(fpp0 + 75);
    const auto *fpp0_76 = buffer.data(fpp0 + 76);
    const auto *fpp0_77 = buffer.data(fpp0 + 77);

    const auto *fpp1_39 = buffer.data(fpp1 + 39);
    const auto *fpp1_45 = buffer.data(fpp1 + 45);
    const auto *fpp1_46 = buffer.data(fpp1 + 46);
    const auto *fpp1_47 = buffer.data(fpp1 + 47);
    const auto *fpp1_57 = buffer.data(fpp1 + 57);
    const auto *fpp1_58 = buffer.data(fpp1 + 58);
    const auto *fpp1_59 = buffer.data(fpp1 + 59);
    const auto *fpp1_68 = buffer.data(fpp1 + 68);
    const auto *fpp1_69 = buffer.data(fpp1 + 69);
    const auto *fpp1_70 = buffer.data(fpp1 + 70);
    const auto *fpp1_71 = buffer.data(fpp1 + 71);
    const auto *fpp1_75 = buffer.data(fpp1 + 75);
    const auto *fpp1_76 = buffer.data(fpp1 + 76);
    const auto *fpp1_77 = buffer.data(fpp1 + 77);

    const auto *fpd_78 = buffer.data(fpd + 78);
    const auto *fpd_81 = buffer.data(fpd + 81);
    const auto *fpd_82 = buffer.data(fpd + 82);
    const auto *fpd_83 = buffer.data(fpd + 83);
    const auto *fpd_84 = buffer.data(fpd + 84);
    const auto *fpd_87 = buffer.data(fpd + 87);
    const auto *fpd_88 = buffer.data(fpd + 88);
    const auto *fpd_89 = buffer.data(fpd + 89);
    const auto *fpd_90 = buffer.data(fpd + 90);
    const auto *fpd_92 = buffer.data(fpd + 92);
    const auto *fpd_93 = buffer.data(fpd + 93);
    const auto *fpd_94 = buffer.data(fpd + 94);
    const auto *fpd_95 = buffer.data(fpd + 95);
    const auto *fpd_96 = buffer.data(fpd + 96);
    const auto *fpd_98 = buffer.data(fpd + 98);
    const auto *fpd_99 = buffer.data(fpd + 99);
    const auto *fpd_101 = buffer.data(fpd + 101);
    const auto *fpd_102 = buffer.data(fpd + 102);
    const auto *fpd_104 = buffer.data(fpd + 104);
    const auto *fpd_105 = buffer.data(fpd + 105);
    const auto *fpd_107 = buffer.data(fpd + 107);
    const auto *fpd_108 = buffer.data(fpd + 108);
    const auto *fpd_111 = buffer.data(fpd + 111);
    const auto *fpd_112 = buffer.data(fpd + 112);
    const auto *fpd_113 = buffer.data(fpd + 113);
    const auto *fpd_114 = buffer.data(fpd + 114);
    const auto *fpd_115 = buffer.data(fpd + 115);
    const auto *fpd_117 = buffer.data(fpd + 117);
    const auto *fpd_118 = buffer.data(fpd + 118);
    const auto *fpd_119 = buffer.data(fpd + 119);
    const auto *fpd_120 = buffer.data(fpd + 120);
    const auto *fpd_123 = buffer.data(fpd + 123);
    const auto *fpd_124 = buffer.data(fpd + 124);
    const auto *fpd_125 = buffer.data(fpd + 125);
    const auto *fpd_129 = buffer.data(fpd + 129);
    const auto *fpd_130 = buffer.data(fpd + 130);
    const auto *fpd_131 = buffer.data(fpd + 131);
    const auto *fpd_134 = buffer.data(fpd + 134);
    const auto *fpd_135 = buffer.data(fpd + 135);
    const auto *fpd_136 = buffer.data(fpd + 136);
    const auto *fpd_137 = buffer.data(fpd + 137);
    const auto *fpd_138 = buffer.data(fpd + 138);
    const auto *fpd_139 = buffer.data(fpd + 139);
    const auto *fpd_140 = buffer.data(fpd + 140);
    const auto *fpd_141 = buffer.data(fpd + 141);
    const auto *fpd_142 = buffer.data(fpd + 142);
    const auto *fpd_143 = buffer.data(fpd + 143);
    const auto *fpd_147 = buffer.data(fpd + 147);
    const auto *fpd_148 = buffer.data(fpd + 148);
    const auto *fpd_149 = buffer.data(fpd + 149);
    const auto *fpd_150 = buffer.data(fpd + 150);
    const auto *fpd_151 = buffer.data(fpd + 151);
    const auto *fpd_152 = buffer.data(fpd + 152);
    const auto *fpd_153 = buffer.data(fpd + 153);
    const auto *fpd_154 = buffer.data(fpd + 154);
    const auto *fpd_155 = buffer.data(fpd + 155);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_z, pc_x, pc_z, dpf0_41, dpd_24, \
                         dpd_78, dpd_81, dpf1_41, fpp0_39, fpp1_39, fpd_78, \
                         fpd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_1 * dpd_78[k]
                   + f_2 * fpp0_39[k]
                   - f_3 * fpp1_39[k]
                   + f_4 * pc_x[k] * fpd_78[k];

        t_131[k] = pa_z[k] * dpf0_41[k]
                   - f_5 * pc_z[k] * dpf1_41[k];

        t_132[k] = f_1 * dpd_24[k]
                   + f_4 * pc_z[k] * fpd_78[k];

        t_133[k] = f_1 * dpd_81[k]
                   + f_4 * pc_x[k] * fpd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pc_x, pc_z, dpf0_136, dpd_27, \
                         dpd_82, dpd_83, dpf1_136, fpd_81, fpd_82, \
                         fpd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_1 * dpd_82[k]
                   + f_4 * pc_x[k] * fpd_82[k];

        t_135[k] = f_1 * dpd_83[k]
                   + f_4 * pc_x[k] * fpd_83[k];

        t_136[k] = pa_x[k] * dpf0_136[k]
                   - f_5 * pc_x[k] * dpf1_136[k];

        t_137[k] = f_1 * dpd_27[k]
                   + f_4 * pc_z[k] * fpd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_x, pa_y, pc_x, pc_y, dpf0_80, \
                         dpf0_138, dpf0_139, dpd_48, dpf1_80, dpf1_138, dpf1_139, \
                         fpd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pa_x[k] * dpf0_138[k]
                   - f_5 * pc_x[k] * dpf1_138[k];

        t_139[k] = pa_x[k] * dpf0_139[k]
                   - f_5 * pc_x[k] * dpf1_139[k];

        t_140[k] = pa_y[k] * dpf0_80[k]
                   - f_5 * pc_y[k] * dpf1_80[k];

        t_141[k] = f_1 * dpd_48[k]
                   + f_4 * pc_y[k] * fpd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_y, pc_x, pc_y, dpf0_82, dpd_87, \
                         dpd_88, dpd_89, dpf1_82, fpd_87, fpd_88, \
                         fpd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = pa_y[k] * dpf0_82[k]
                   - f_5 * pc_y[k] * dpf1_82[k];

        t_143[k] = f_1 * dpd_87[k]
                   + f_4 * pc_x[k] * fpd_87[k];

        t_144[k] = f_1 * dpd_88[k]
                   + f_4 * pc_x[k] * fpd_88[k];

        t_145[k] = f_1 * dpd_89[k]
                   + f_4 * pc_x[k] * fpd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pc_x, pc_y, dpf0_146, dpf0_147, \
                         dpf0_149, dpd_53, dpf1_146, dpf1_147, dpf1_149, \
                         fpd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * dpf0_146[k]
                   - f_5 * pc_x[k] * dpf1_146[k];

        t_147[k] = pa_x[k] * dpf0_147[k]
                   - f_5 * pc_x[k] * dpf1_147[k];

        t_148[k] = f_1 * dpd_53[k]
                   + f_4 * pc_y[k] * fpd_89[k];

        t_149[k] = pa_x[k] * dpf0_149[k]
                   - f_5 * pc_x[k] * dpf1_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, pc_z, dpd_36, dpd_90, dpd_93, \
                         fsd_30, fsd_33, fpp0_45, fpp1_45, fpd_90, \
                         fpd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * dpd_90[k]
                   + f_1 * fsd_30[k]
                   + f_2 * fpp0_45[k]
                   - f_3 * fpp1_45[k]
                   + f_4 * pc_x[k] * fpd_90[k];

        t_151[k] = f_4 * pc_y[k] * fpd_90[k];

        t_152[k] = f_6 * dpd_36[k]
                   + f_4 * pc_z[k] * fpd_90[k];

        t_153[k] = f_1 * dpd_93[k]
                   + f_1 * fsd_33[k]
                   + f_4 * pc_x[k] * fpd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, dpd_95, fsd_35, fpp0_46, \
                         fpp0_47, fpp1_46, fpp1_47, fpd_92, fpd_93, fpd_94, \
                         fpd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_4 * pc_y[k] * fpd_92[k];

        t_155[k] = f_1 * dpd_95[k]
                   + f_1 * fsd_35[k]
                   + f_4 * pc_x[k] * fpd_95[k];

        t_156[k] = f_2 * fpp0_46[k]
                   - f_3 * fpp1_46[k]
                   + f_4 * pc_y[k] * fpd_93[k];

        t_157[k] = f_9 * fpp0_47[k]
                   - f_10 * fpp1_47[k]
                   + f_4 * pc_y[k] * fpd_94[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_y, pc_y, pc_z, dpd_41, fsf0_50, \
                         fsd_30, fsf1_50, fpp0_47, fpp1_47, fpd_95, \
                         fpd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_4 * pc_y[k] * fpd_95[k];

        t_159[k] = f_6 * dpd_41[k]
                   + f_2 * fpp0_47[k]
                   - f_3 * fpp1_47[k]
                   + f_4 * pc_z[k] * fpd_95[k];

        t_160[k] = pb_y[k] * fsf0_50[k]
                   - f_5 * pc_y[k] * fsf1_50[k];

        t_161[k] = f_1 * fsd_30[k]
                   + f_4 * pc_y[k] * fpd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_y, pc_z, dpd_42, dpd_99, \
                         dpd_101, fsd_32, fpd_96, fpd_98, fpd_99, \
                         fpd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_6 * dpd_42[k]
                   + f_4 * pc_z[k] * fpd_96[k];

        t_163[k] = f_1 * dpd_99[k]
                   + f_4 * pc_x[k] * fpd_99[k];

        t_164[k] = f_1 * fsd_32[k]
                   + f_4 * pc_y[k] * fpd_98[k];

        t_165[k] = f_1 * dpd_101[k]
                   + f_4 * pc_x[k] * fpd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_x, pc_x, pc_y, dpf0_166, dpf0_167, \
                         dpf0_169, dpf1_166, dpf1_167, dpf1_169, fsd_35, \
                         fpd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_x[k] * dpf0_166[k]
                   - f_5 * pc_x[k] * dpf1_166[k];

        t_167[k] = pa_x[k] * dpf0_167[k]
                   - f_5 * pc_x[k] * dpf1_167[k];

        t_168[k] = f_1 * fsd_35[k]
                   + f_4 * pc_y[k] * fpd_101[k];

        t_169[k] = pa_x[k] * dpf0_169[k]
                   - f_5 * pc_x[k] * dpf1_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, pc_x, pc_y, pc_z, dpf0_170, dpd_48, \
                         dpd_102, dpd_105, dpf1_170, fsd_30, fpd_102, \
                         fpd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pa_x[k] * dpf0_170[k]
                   + f_0 * dpd_102[k]
                   - f_5 * pc_x[k] * dpf1_170[k];

        t_171[k] = f_4 * pc_y[k] * fpd_102[k];

        t_172[k] = f_6 * dpd_48[k]
                   + f_1 * fsd_30[k]
                   + f_4 * pc_z[k] * fpd_102[k];

        t_173[k] = f_1 * dpd_105[k]
                   + f_4 * pc_x[k] * fpd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pa_x, pc_x, pc_y, dpf0_176, \
                         dpf0_177, dpd_107, dpf1_176, dpf1_177, fpd_104, \
                         fpd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_4 * pc_y[k] * fpd_104[k];

        t_175[k] = f_1 * dpd_107[k]
                   + f_4 * pc_x[k] * fpd_107[k];

        t_176[k] = pa_x[k] * dpf0_176[k]
                   - f_5 * pc_x[k] * dpf1_176[k];

        t_177[k] = pa_x[k] * dpf0_177[k]
                   - f_5 * pc_x[k] * dpf1_177[k];

        t_178[k] = f_4 * pc_y[k] * fpd_107[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pb_x, pc_x, dpf0_179, dpf1_179, fsf0_60, \
                         fsf0_61, fsd_36, fsd_37, fsf1_60, fsf1_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_x[k] * dpf0_179[k]
                   - f_5 * pc_x[k] * dpf1_179[k];

        t_180[k] = pb_x[k] * fsf0_60[k]
                   + f_0 * fsd_36[k]
                   - f_5 * pc_x[k] * fsf1_60[k];

        t_181[k] = pb_x[k] * fsf0_61[k]
                   + f_6 * fsd_37[k]
                   - f_5 * pc_x[k] * fsf1_61[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_z, fsd_39, fsd_40, fsd_41, \
                         fpd_108, fpd_111, fpd_112, fpd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_4 * pc_z[k] * fpd_108[k];

        t_183[k] = f_1 * fsd_39[k]
                   + f_4 * pc_x[k] * fpd_111[k];

        t_184[k] = f_1 * fsd_40[k]
                   + f_4 * pc_x[k] * fpd_112[k];

        t_185[k] = f_1 * fsd_41[k]
                   + f_4 * pc_x[k] * fpd_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pb_x, pc_x, pc_y, pc_z, dpd_59, fsf0_66, \
                         fsf0_69, fsf1_66, fsf1_69, fpd_111, fpd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pb_x[k] * fsf0_66[k]
                   - f_5 * pc_x[k] * fsf1_66[k];

        t_187[k] = f_4 * pc_z[k] * fpd_111[k];

        t_188[k] = f_0 * dpd_59[k]
                   + f_4 * pc_y[k] * fpd_113[k];

        t_189[k] = pb_x[k] * fsf0_69[k]
                   - f_5 * pc_x[k] * fsf1_69[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pc_x, pc_z, fpp0_57, fpp0_58, \
                         fpp1_57, fpp1_58, fpd_114, fpd_115, fpd_117, \
                         fpd_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_2 * fpp0_57[k]
                   - f_3 * fpp1_57[k]
                   + f_4 * pc_x[k] * fpd_114[k];

        t_191[k] = f_9 * fpp0_58[k]
                   - f_10 * fpp1_58[k]
                   + f_4 * pc_x[k] * fpd_115[k];

        t_192[k] = f_4 * pc_z[k] * fpd_114[k];

        t_193[k] = f_4 * pc_x[k] * fpd_117[k];

        t_194[k] = f_4 * pc_x[k] * fpd_118[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pc_x, pc_y, pc_z, dpd_63, dpd_65, fsd_39, \
                         fsd_41, fpp0_58, fpp1_58, fpd_117, fpd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_4 * pc_x[k] * fpd_119[k];

        t_196[k] = f_0 * dpd_63[k]
                   + f_1 * fsd_39[k]
                   + f_2 * fpp0_58[k]
                   - f_3 * fpp1_58[k]
                   + f_4 * pc_y[k] * fpd_117[k];

        t_197[k] = f_4 * pc_z[k] * fpd_117[k];

        t_198[k] = f_0 * dpd_65[k]
                   + f_1 * fsd_41[k]
                   + f_4 * pc_y[k] * fpd_119[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_z, pc_z, fsf0_60, fsf0_61, fsd_36, \
                         fsf1_60, fsf1_61, fpp0_59, fpp1_59, fpd_119, \
                         fpd_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_2 * fpp0_59[k]
                   - f_3 * fpp1_59[k]
                   + f_4 * pc_z[k] * fpd_119[k];

        t_200[k] = pb_z[k] * fsf0_60[k]
                   - f_5 * pc_z[k] * fsf1_60[k];

        t_201[k] = pb_z[k] * fsf0_61[k]
                   - f_5 * pc_z[k] * fsf1_61[k];

        t_202[k] = f_1 * fsd_36[k]
                   + f_4 * pc_z[k] * fpd_120[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_z, pc_x, pc_z, fsf0_66, fsd_39, \
                         fsf1_66, fpd_123, fpd_124, fpd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_4 * pc_x[k] * fpd_123[k];

        t_204[k] = f_4 * pc_x[k] * fpd_124[k];

        t_205[k] = f_4 * pc_x[k] * fpd_125[k];

        t_206[k] = pb_z[k] * fsf0_66[k]
                   - f_5 * pc_z[k] * fsf1_66[k];

        t_207[k] = f_1 * fsd_39[k]
                   + f_4 * pc_z[k] * fpd_123[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_z, pb_z, pc_y, pc_z, dpf0_90, dpd_71, \
                         dpf1_90, fsf0_69, fsd_41, fsf1_69, fpd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_0 * dpd_71[k]
                   + f_4 * pc_y[k] * fpd_125[k];

        t_209[k] = pb_z[k] * fsf0_69[k]
                   + f_0 * fsd_41[k]
                   - f_5 * pc_z[k] * fsf1_69[k];

        t_210[k] = pa_z[k] * dpf0_90[k]
                   - f_5 * pc_z[k] * dpf1_90[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pa_z, pb_x, pc_x, pc_z, dpf0_91, dpf1_91, \
                         fsf0_72, fsd_44, fsd_45, fsf1_72, fpd_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = pa_z[k] * dpf0_91[k]
                   - f_5 * pc_z[k] * dpf1_91[k];

        t_212[k] = pb_x[k] * fsf0_72[k]
                   + f_6 * fsd_44[k]
                   - f_5 * pc_x[k] * fsf1_72[k];

        t_213[k] = f_1 * fsd_45[k]
                   + f_4 * pc_x[k] * fpd_129[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, dpf0_96, dpd_57, \
                         dpf1_96, fsd_46, fsd_47, fpd_129, fpd_130, \
                         fpd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_1 * fsd_46[k]
                   + f_4 * pc_x[k] * fpd_130[k];

        t_215[k] = f_1 * fsd_47[k]
                   + f_4 * pc_x[k] * fpd_131[k];

        t_216[k] = pa_z[k] * dpf0_96[k]
                   - f_5 * pc_z[k] * dpf1_96[k];

        t_217[k] = f_1 * dpd_57[k]
                   + f_4 * pc_z[k] * fpd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pa_z, pb_x, pc_x, pc_y, pc_z, dpf0_100, dpd_77, \
                         dpf1_100, fsf0_79, fsf1_79, fpd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_6 * dpd_77[k]
                   + f_4 * pc_y[k] * fpd_131[k];

        t_219[k] = pb_x[k] * fsf0_79[k]
                   - f_5 * pc_x[k] * fsf1_79[k];

        t_220[k] = pa_z[k] * dpf0_100[k]
                   - f_5 * pc_z[k] * dpf1_100[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, pa_z, pc_x, pc_z, dpf0_101, \
                         dpf1_101, fpp0_68, fpp1_68, fpd_134, fpd_135, fpd_136, \
                         fpd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pa_z[k] * dpf0_101[k]
                   - f_5 * pc_z[k] * dpf1_101[k];

        t_222[k] = f_9 * fpp0_68[k]
                   - f_10 * fpp1_68[k]
                   + f_4 * pc_x[k] * fpd_134[k];

        t_223[k] = f_4 * pc_x[k] * fpd_135[k];

        t_224[k] = f_4 * pc_x[k] * fpd_136[k];

        t_225[k] = f_4 * pc_x[k] * fpd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pa_z, pc_y, pc_z, dpf0_106, dpd_63, dpd_83, \
                         dpf1_106, fsd_47, fpd_135, fpd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * dpf0_106[k]
                   - f_5 * pc_z[k] * dpf1_106[k];

        t_227[k] = f_1 * dpd_63[k]
                   + f_4 * pc_z[k] * fpd_135[k];

        t_228[k] = f_6 * dpd_83[k]
                   + f_1 * fsd_47[k]
                   + f_4 * pc_y[k] * fpd_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_z, dpd_65, fpp0_68, fpp0_69, fpp0_70, \
                         fpp1_68, fpp1_69, fpp1_70, fpd_137, fpd_138, \
                         fpd_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_1 * dpd_65[k]
                   + f_2 * fpp0_68[k]
                   - f_3 * fpp1_68[k]
                   + f_4 * pc_z[k] * fpd_137[k];

        t_230[k] = f_2 * fpp0_69[k]
                   - f_3 * fpp1_69[k]
                   + f_4 * pc_x[k] * fpd_138[k];

        t_231[k] = f_9 * fpp0_70[k]
                   - f_10 * fpp1_70[k]
                   + f_4 * pc_x[k] * fpd_139[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pc_x, pc_y, dpd_87, fpp0_70, \
                         fpp0_71, fpp1_70, fpp1_71, fpd_140, fpd_141, fpd_142, \
                         fpd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_9 * fpp0_71[k]
                   - f_10 * fpp1_71[k]
                   + f_4 * pc_x[k] * fpd_140[k];

        t_233[k] = f_4 * pc_x[k] * fpd_141[k];

        t_234[k] = f_4 * pc_x[k] * fpd_142[k];

        t_235[k] = f_4 * pc_x[k] * fpd_143[k];

        t_236[k] = f_6 * dpd_87[k]
                   + f_2 * fpp0_70[k]
                   - f_3 * fpp1_70[k]
                   + f_4 * pc_y[k] * fpd_141[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_y, pc_y, pc_z, ppf0_89, ppf1_89, dpf0_149, \
                         dpd_69, dpd_89, dpf1_149, fsd_45, fpd_141, \
                         fpd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_1 * dpd_69[k]
                   + f_1 * fsd_45[k]
                   + f_4 * pc_z[k] * fpd_141[k];

        t_238[k] = f_6 * dpd_89[k]
                   + f_4 * pc_y[k] * fpd_143[k];

        t_239[k] = f_7 * ppf0_89[k]
                   - f_8 * ppf1_89[k]
                   + pa_y[k] * dpf0_149[k]
                   - f_5 * pc_y[k] * dpf1_149[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pa_y, pc_x, pc_y, dpf0_150, dpf0_151, \
                         dpf0_152, dpd_90, dpf1_150, dpf1_151, dpf1_152, fsd_51, \
                         fpd_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_y[k] * dpf0_150[k]
                   - f_5 * pc_y[k] * dpf1_150[k];

        t_241[k] = pa_y[k] * dpf0_151[k]
                   + f_1 * dpd_90[k]
                   - f_5 * pc_y[k] * dpf1_151[k];

        t_242[k] = pa_y[k] * dpf0_152[k]
                   - f_5 * pc_y[k] * dpf1_152[k];

        t_243[k] = f_1 * fsd_51[k]
                   + f_4 * pc_x[k] * fpd_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pc_x, pc_z, dpd_75, fsf0_86, \
                         fsd_52, fsd_53, fsf1_86, fpd_147, fpd_148, \
                         fpd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_1 * fsd_52[k]
                   + f_4 * pc_x[k] * fpd_148[k];

        t_245[k] = f_1 * fsd_53[k]
                   + f_4 * pc_x[k] * fpd_149[k];

        t_246[k] = pb_x[k] * fsf0_86[k]
                   - f_5 * pc_x[k] * fsf1_86[k];

        t_247[k] = f_6 * dpd_75[k]
                   + f_4 * pc_z[k] * fpd_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pa_y, pc_x, pc_y, dpf0_159, dpd_95, dpf1_159, \
                         fpp0_75, fpp1_75, fpd_149, fpd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_1 * dpd_95[k]
                   + f_4 * pc_y[k] * fpd_149[k];

        t_249[k] = pa_y[k] * dpf0_159[k]
                   - f_5 * pc_y[k] * dpf1_159[k];

        t_250[k] = f_2 * fpp0_75[k]
                   - f_3 * fpp1_75[k]
                   + f_4 * pc_x[k] * fpd_150[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, pc_x, fpp0_76, fpp0_77, fpp1_76, \
                         fpp1_77, fpd_151, fpd_152, fpd_153, fpd_154, \
                         fpd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * fpp0_76[k]
                   - f_10 * fpp1_76[k]
                   + f_4 * pc_x[k] * fpd_151[k];

        t_252[k] = f_9 * fpp0_77[k]
                   - f_10 * fpp1_77[k]
                   + f_4 * pc_x[k] * fpd_152[k];

        t_253[k] = f_4 * pc_x[k] * fpd_153[k];

        t_254[k] = f_4 * pc_x[k] * fpd_154[k];

        t_255[k] = f_4 * pc_x[k] * fpd_155[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pc_y, pc_z, dpd_81, dpd_99, dpd_101, fsd_51, \
                         fsd_53, fpp0_76, fpp1_76, fpd_153, fpd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_1 * dpd_99[k]
                   + f_1 * fsd_51[k]
                   + f_2 * fpp0_76[k]
                   - f_3 * fpp1_76[k]
                   + f_4 * pc_y[k] * fpd_153[k];

        t_257[k] = f_6 * dpd_81[k]
                   + f_4 * pc_z[k] * fpd_153[k];

        t_258[k] = f_1 * dpd_101[k]
                   + f_1 * fsd_53[k]
                   + f_4 * pc_y[k] * fpd_155[k];
    }
}

static auto
compute_prim_fpf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpf0, const size_t dpd,
                                                          const size_t dpf1, const size_t fsf0,
                                                          const size_t fsd, const size_t fsf1,
                                                          const size_t fpp0, const size_t fpp1,
                                                          const size_t fpd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpf0_170 = buffer.data(dpf0 + 170);
    const auto *dpf0_172 = buffer.data(dpf0 + 172);
    const auto *dpf0_176 = buffer.data(dpf0 + 176);
    const auto *dpf0_179 = buffer.data(dpf0 + 179);

    const auto *dpd_83 = buffer.data(dpd + 83);
    const auto *dpd_87 = buffer.data(dpd + 87);
    const auto *dpd_105 = buffer.data(dpd + 105);
    const auto *dpd_107 = buffer.data(dpd + 107);

    const auto *dpf1_170 = buffer.data(dpf1 + 170);
    const auto *dpf1_172 = buffer.data(dpf1 + 172);
    const auto *dpf1_176 = buffer.data(dpf1 + 176);
    const auto *dpf1_179 = buffer.data(dpf1 + 179);

    const auto *fsf0_90 = buffer.data(fsf0 + 90);
    const auto *fsf0_92 = buffer.data(fsf0 + 92);
    const auto *fsf0_96 = buffer.data(fsf0 + 96);
    const auto *fsf0_97 = buffer.data(fsf0 + 97);
    const auto *fsf0_99 = buffer.data(fsf0 + 99);

    const auto *fsd_51 = buffer.data(fsd + 51);
    const auto *fsd_54 = buffer.data(fsd + 54);
    const auto *fsd_56 = buffer.data(fsd + 56);
    const auto *fsd_57 = buffer.data(fsd + 57);
    const auto *fsd_58 = buffer.data(fsd + 58);
    const auto *fsd_59 = buffer.data(fsd + 59);

    const auto *fsf1_90 = buffer.data(fsf1 + 90);
    const auto *fsf1_92 = buffer.data(fsf1 + 92);
    const auto *fsf1_96 = buffer.data(fsf1 + 96);
    const auto *fsf1_97 = buffer.data(fsf1 + 97);
    const auto *fsf1_99 = buffer.data(fsf1 + 99);

    const auto *fpp0_77 = buffer.data(fpp0 + 77);
    const auto *fpp0_79 = buffer.data(fpp0 + 79);
    const auto *fpp0_87 = buffer.data(fpp0 + 87);
    const auto *fpp0_88 = buffer.data(fpp0 + 88);
    const auto *fpp0_89 = buffer.data(fpp0 + 89);

    const auto *fpp1_77 = buffer.data(fpp1 + 77);
    const auto *fpp1_79 = buffer.data(fpp1 + 79);
    const auto *fpp1_87 = buffer.data(fpp1 + 87);
    const auto *fpp1_88 = buffer.data(fpp1 + 88);
    const auto *fpp1_89 = buffer.data(fpp1 + 89);

    const auto *fpd_155 = buffer.data(fpd + 155);
    const auto *fpd_157 = buffer.data(fpd + 157);
    const auto *fpd_159 = buffer.data(fpd + 159);
    const auto *fpd_160 = buffer.data(fpd + 160);
    const auto *fpd_161 = buffer.data(fpd + 161);
    const auto *fpd_162 = buffer.data(fpd + 162);
    const auto *fpd_165 = buffer.data(fpd + 165);
    const auto *fpd_166 = buffer.data(fpd + 166);
    const auto *fpd_167 = buffer.data(fpd + 167);
    const auto *fpd_168 = buffer.data(fpd + 168);
    const auto *fpd_171 = buffer.data(fpd + 171);
    const auto *fpd_172 = buffer.data(fpd + 172);
    const auto *fpd_173 = buffer.data(fpd + 173);
    const auto *fpd_174 = buffer.data(fpd + 174);
    const auto *fpd_176 = buffer.data(fpd + 176);
    const auto *fpd_177 = buffer.data(fpd + 177);
    const auto *fpd_178 = buffer.data(fpd + 178);
    const auto *fpd_179 = buffer.data(fpd + 179);

#pragma omp simd aligned(t_259, t_260, t_261, pa_y, pc_x, pc_y, pc_z, dpf0_170, dpd_83, \
                         dpf1_170, fpp0_77, fpp0_79, fpp1_77, fpp1_79, fpd_155, \
                         fpd_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_6 * dpd_83[k]
                   + f_2 * fpp0_77[k]
                   - f_3 * fpp1_77[k]
                   + f_4 * pc_z[k] * fpd_155[k];

        t_260[k] = pa_y[k] * dpf0_170[k]
                   - f_5 * pc_y[k] * dpf1_170[k];

        t_261[k] = f_9 * fpp0_79[k]
                   - f_10 * fpp1_79[k]
                   + f_4 * pc_x[k] * fpd_157[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pa_y, pc_x, pc_y, dpf0_172, \
                         dpf0_176, dpd_105, dpf1_172, dpf1_176, fpd_159, fpd_160, \
                         fpd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = pa_y[k] * dpf0_172[k]
                   - f_5 * pc_y[k] * dpf1_172[k];

        t_263[k] = f_4 * pc_x[k] * fpd_159[k];

        t_264[k] = f_4 * pc_x[k] * fpd_160[k];

        t_265[k] = f_4 * pc_x[k] * fpd_161[k];

        t_266[k] = pa_y[k] * dpf0_176[k]
                   + f_0 * dpd_105[k]
                   - f_5 * pc_y[k] * dpf1_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_y, pc_y, pc_z, dpf0_179, dpd_87, dpd_107, \
                         dpf1_179, fsd_51, fpd_159, fpd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_6 * dpd_87[k]
                   + f_1 * fsd_51[k]
                   + f_4 * pc_z[k] * fpd_159[k];

        t_268[k] = f_1 * dpd_107[k]
                   + f_4 * pc_y[k] * fpd_161[k];

        t_269[k] = pa_y[k] * dpf0_179[k]
                   - f_5 * pc_y[k] * dpf1_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_x, pc_x, pc_y, fsf0_90, fsf0_92, \
                         fsd_54, fsd_56, fsd_57, fsf1_90, fsf1_92, fpd_162, \
                         fpd_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = pb_x[k] * fsf0_90[k]
                   + f_0 * fsd_54[k]
                   - f_5 * pc_x[k] * fsf1_90[k];

        t_271[k] = f_4 * pc_y[k] * fpd_162[k];

        t_272[k] = pb_x[k] * fsf0_92[k]
                   + f_6 * fsd_56[k]
                   - f_5 * pc_x[k] * fsf1_92[k];

        t_273[k] = f_1 * fsd_57[k]
                   + f_4 * pc_x[k] * fpd_165[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_x, pc_x, pc_y, fsf0_96, \
                         fsf0_97, fsd_58, fsd_59, fsf1_96, fsf1_97, fpd_166, \
                         fpd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_1 * fsd_58[k]
                   + f_4 * pc_x[k] * fpd_166[k];

        t_275[k] = f_1 * fsd_59[k]
                   + f_4 * pc_x[k] * fpd_167[k];

        t_276[k] = pb_x[k] * fsf0_96[k]
                   - f_5 * pc_x[k] * fsf1_96[k];

        t_277[k] = pb_x[k] * fsf0_97[k]
                   - f_5 * pc_x[k] * fsf1_97[k];

        t_278[k] = f_4 * pc_y[k] * fpd_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_x, pb_y, pc_x, pc_y, fsf0_90, fsf0_92, \
                         fsf0_99, fsd_54, fsf1_90, fsf1_92, fsf1_99, \
                         fpd_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_x[k] * fsf0_99[k]
                   - f_5 * pc_x[k] * fsf1_99[k];

        t_280[k] = pb_y[k] * fsf0_90[k]
                   - f_5 * pc_y[k] * fsf1_90[k];

        t_281[k] = f_1 * fsd_54[k]
                   + f_4 * pc_y[k] * fpd_168[k];

        t_282[k] = pb_y[k] * fsf0_92[k]
                   - f_5 * pc_y[k] * fsf1_92[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_y, pc_x, pc_y, fsf0_96, fsd_57, \
                         fsf1_96, fpd_171, fpd_172, fpd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_4 * pc_x[k] * fpd_171[k];

        t_284[k] = f_4 * pc_x[k] * fpd_172[k];

        t_285[k] = f_4 * pc_x[k] * fpd_173[k];

        t_286[k] = pb_y[k] * fsf0_96[k]
                   + f_0 * fsd_57[k]
                   - f_5 * pc_y[k] * fsf1_96[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pb_y, pc_y, fsf0_97, fsf0_99, fsd_58, fsd_59, \
                         fsf1_97, fsf1_99, fpd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pb_y[k] * fsf0_97[k]
                   + f_6 * fsd_58[k]
                   - f_5 * pc_y[k] * fsf1_97[k];

        t_288[k] = f_1 * fsd_59[k]
                   + f_4 * pc_y[k] * fpd_173[k];

        t_289[k] = pb_y[k] * fsf0_99[k]
                   - f_5 * pc_y[k] * fsf1_99[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pc_x, pc_y, fpp0_87, fpp0_89, \
                         fpp1_87, fpp1_89, fpd_174, fpd_176, fpd_177, \
                         fpd_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_2 * fpp0_87[k]
                   - f_3 * fpp1_87[k]
                   + f_4 * pc_x[k] * fpd_174[k];

        t_291[k] = f_4 * pc_y[k] * fpd_174[k];

        t_292[k] = f_9 * fpp0_89[k]
                   - f_10 * fpp1_89[k]
                   + f_4 * pc_x[k] * fpd_176[k];

        t_293[k] = f_4 * pc_x[k] * fpd_177[k];

        t_294[k] = f_4 * pc_x[k] * fpd_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_x, pc_y, fpp0_88, fpp0_89, fpp1_88, \
                         fpp1_89, fpd_177, fpd_178, fpd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_4 * pc_x[k] * fpd_179[k];

        t_296[k] = f_2 * fpp0_88[k]
                   - f_3 * fpp1_88[k]
                   + f_4 * pc_y[k] * fpd_177[k];

        t_297[k] = f_9 * fpp0_89[k]
                   - f_10 * fpp1_89[k]
                   + f_4 * pc_y[k] * fpd_178[k];

        t_298[k] = f_4 * pc_y[k] * fpd_179[k];
    }

#pragma omp simd aligned(t_299, pc_z, dpd_107, fsd_59, fpp0_89, fpp1_89, \
                         fpd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * dpd_107[k]
                   + f_1 * fsd_59[k]
                   + f_2 * fpp0_89[k]
                   - f_3 * fpp1_89[k]
                   + f_4 * pc_z[k] * fpd_179[k];
    }
}

auto
compute_prim_fpf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppf0,
                                                   const size_t ppf1, const size_t dpf0,
                                                   const size_t dpd, const size_t dpf1,
                                                   const size_t fsf0, const size_t fsd,
                                                   const size_t fsf1, const size_t fpp0,
                                                   const size_t fpp1, const size_t fpd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fpf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, ppf0,
                                                              ppf1, dpf0, dpd, dpf1, fsf0, fsd,
                                                              fsf1, fpp0, fpp1, fpd, ncols,
                                                              gamma, p, q);

    compute_prim_fpf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppf0,
                                                              ppf1, dpf0, dpd, dpf1, fsf0, fsd,
                                                              fsf1, fpp0, fpp1, fpd, ncols,
                                                              gamma, p, q);

    compute_prim_fpf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, dpf0,
                                                              dpd, dpf1, fsf0, fsd, fsf1, fpp0,
                                                              fpp1, fpd, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
