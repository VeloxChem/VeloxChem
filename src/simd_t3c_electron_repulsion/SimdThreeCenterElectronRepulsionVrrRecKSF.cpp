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


#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isf0,
                                                          const size_t isd, const size_t isf1,
                                                          const size_t ksp0, const size_t ksp1,
                                                          const size_t ksd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *isf0_0 = buffer.data(isf0 + 0);
    const auto *isf0_6 = buffer.data(isf0 + 6);
    const auto *isf0_9 = buffer.data(isf0 + 9);
    const auto *isf0_16 = buffer.data(isf0 + 16);
    const auto *isf0_20 = buffer.data(isf0 + 20);
    const auto *isf0_29 = buffer.data(isf0 + 29);
    const auto *isf0_30 = buffer.data(isf0 + 30);
    const auto *isf0_36 = buffer.data(isf0 + 36);
    const auto *isf0_50 = buffer.data(isf0 + 50);
    const auto *isf0_59 = buffer.data(isf0 + 59);
    const auto *isf0_60 = buffer.data(isf0 + 60);
    const auto *isf0_66 = buffer.data(isf0 + 66);
    const auto *isf0_90 = buffer.data(isf0 + 90);

    const auto *isd_0 = buffer.data(isd + 0);
    const auto *isd_3 = buffer.data(isd + 3);
    const auto *isd_5 = buffer.data(isd + 5);
    const auto *isd_6 = buffer.data(isd + 6);
    const auto *isd_9 = buffer.data(isd + 9);
    const auto *isd_11 = buffer.data(isd + 11);
    const auto *isd_12 = buffer.data(isd + 12);
    const auto *isd_15 = buffer.data(isd + 15);
    const auto *isd_17 = buffer.data(isd + 17);
    const auto *isd_18 = buffer.data(isd + 18);
    const auto *isd_21 = buffer.data(isd + 21);
    const auto *isd_23 = buffer.data(isd + 23);
    const auto *isd_24 = buffer.data(isd + 24);
    const auto *isd_27 = buffer.data(isd + 27);
    const auto *isd_28 = buffer.data(isd + 28);
    const auto *isd_29 = buffer.data(isd + 29);
    const auto *isd_30 = buffer.data(isd + 30);
    const auto *isd_33 = buffer.data(isd + 33);
    const auto *isd_35 = buffer.data(isd + 35);
    const auto *isd_36 = buffer.data(isd + 36);
    const auto *isd_39 = buffer.data(isd + 39);
    const auto *isd_41 = buffer.data(isd + 41);
    const auto *isd_42 = buffer.data(isd + 42);
    const auto *isd_45 = buffer.data(isd + 45);
    const auto *isd_46 = buffer.data(isd + 46);
    const auto *isd_47 = buffer.data(isd + 47);
    const auto *isd_48 = buffer.data(isd + 48);
    const auto *isd_51 = buffer.data(isd + 51);
    const auto *isd_52 = buffer.data(isd + 52);
    const auto *isd_53 = buffer.data(isd + 53);
    const auto *isd_54 = buffer.data(isd + 54);
    const auto *isd_57 = buffer.data(isd + 57);
    const auto *isd_59 = buffer.data(isd + 59);
    const auto *isd_60 = buffer.data(isd + 60);
    const auto *isd_63 = buffer.data(isd + 63);
    const auto *isd_65 = buffer.data(isd + 65);
    const auto *isd_69 = buffer.data(isd + 69);
    const auto *isd_70 = buffer.data(isd + 70);
    const auto *isd_71 = buffer.data(isd + 71);
    const auto *isd_72 = buffer.data(isd + 72);
    const auto *isd_75 = buffer.data(isd + 75);
    const auto *isd_76 = buffer.data(isd + 76);
    const auto *isd_77 = buffer.data(isd + 77);

    const auto *isf1_0 = buffer.data(isf1 + 0);
    const auto *isf1_6 = buffer.data(isf1 + 6);
    const auto *isf1_9 = buffer.data(isf1 + 9);
    const auto *isf1_16 = buffer.data(isf1 + 16);
    const auto *isf1_20 = buffer.data(isf1 + 20);
    const auto *isf1_29 = buffer.data(isf1 + 29);
    const auto *isf1_30 = buffer.data(isf1 + 30);
    const auto *isf1_36 = buffer.data(isf1 + 36);
    const auto *isf1_50 = buffer.data(isf1 + 50);
    const auto *isf1_59 = buffer.data(isf1 + 59);
    const auto *isf1_60 = buffer.data(isf1 + 60);
    const auto *isf1_66 = buffer.data(isf1 + 66);
    const auto *isf1_90 = buffer.data(isf1 + 90);

    const auto *ksp0_0 = buffer.data(ksp0 + 0);
    const auto *ksp0_1 = buffer.data(ksp0 + 1);
    const auto *ksp0_2 = buffer.data(ksp0 + 2);
    const auto *ksp0_4 = buffer.data(ksp0 + 4);
    const auto *ksp0_8 = buffer.data(ksp0 + 8);
    const auto *ksp0_9 = buffer.data(ksp0 + 9);
    const auto *ksp0_10 = buffer.data(ksp0 + 10);
    const auto *ksp0_11 = buffer.data(ksp0 + 11);
    const auto *ksp0_15 = buffer.data(ksp0 + 15);
    const auto *ksp0_16 = buffer.data(ksp0 + 16);
    const auto *ksp0_17 = buffer.data(ksp0 + 17);
    const auto *ksp0_18 = buffer.data(ksp0 + 18);
    const auto *ksp0_19 = buffer.data(ksp0 + 19);
    const auto *ksp0_20 = buffer.data(ksp0 + 20);
    const auto *ksp0_23 = buffer.data(ksp0 + 23);
    const auto *ksp0_25 = buffer.data(ksp0 + 25);
    const auto *ksp0_27 = buffer.data(ksp0 + 27);
    const auto *ksp0_28 = buffer.data(ksp0 + 28);
    const auto *ksp0_29 = buffer.data(ksp0 + 29);
    const auto *ksp0_30 = buffer.data(ksp0 + 30);
    const auto *ksp0_31 = buffer.data(ksp0 + 31);
    const auto *ksp0_32 = buffer.data(ksp0 + 32);
    const auto *ksp0_35 = buffer.data(ksp0 + 35);
    const auto *ksp0_36 = buffer.data(ksp0 + 36);
    const auto *ksp0_37 = buffer.data(ksp0 + 37);
    const auto *ksp0_38 = buffer.data(ksp0 + 38);

    const auto *ksp1_0 = buffer.data(ksp1 + 0);
    const auto *ksp1_1 = buffer.data(ksp1 + 1);
    const auto *ksp1_2 = buffer.data(ksp1 + 2);
    const auto *ksp1_4 = buffer.data(ksp1 + 4);
    const auto *ksp1_8 = buffer.data(ksp1 + 8);
    const auto *ksp1_9 = buffer.data(ksp1 + 9);
    const auto *ksp1_10 = buffer.data(ksp1 + 10);
    const auto *ksp1_11 = buffer.data(ksp1 + 11);
    const auto *ksp1_15 = buffer.data(ksp1 + 15);
    const auto *ksp1_16 = buffer.data(ksp1 + 16);
    const auto *ksp1_17 = buffer.data(ksp1 + 17);
    const auto *ksp1_18 = buffer.data(ksp1 + 18);
    const auto *ksp1_19 = buffer.data(ksp1 + 19);
    const auto *ksp1_20 = buffer.data(ksp1 + 20);
    const auto *ksp1_23 = buffer.data(ksp1 + 23);
    const auto *ksp1_25 = buffer.data(ksp1 + 25);
    const auto *ksp1_27 = buffer.data(ksp1 + 27);
    const auto *ksp1_28 = buffer.data(ksp1 + 28);
    const auto *ksp1_29 = buffer.data(ksp1 + 29);
    const auto *ksp1_30 = buffer.data(ksp1 + 30);
    const auto *ksp1_31 = buffer.data(ksp1 + 31);
    const auto *ksp1_32 = buffer.data(ksp1 + 32);
    const auto *ksp1_35 = buffer.data(ksp1 + 35);
    const auto *ksp1_36 = buffer.data(ksp1 + 36);
    const auto *ksp1_37 = buffer.data(ksp1 + 37);
    const auto *ksp1_38 = buffer.data(ksp1 + 38);

    const auto *ksd_0 = buffer.data(ksd + 0);
    const auto *ksd_2 = buffer.data(ksd + 2);
    const auto *ksd_3 = buffer.data(ksd + 3);
    const auto *ksd_5 = buffer.data(ksd + 5);
    const auto *ksd_6 = buffer.data(ksd + 6);
    const auto *ksd_7 = buffer.data(ksd + 7);
    const auto *ksd_9 = buffer.data(ksd + 9);
    const auto *ksd_11 = buffer.data(ksd + 11);
    const auto *ksd_12 = buffer.data(ksd + 12);
    const auto *ksd_14 = buffer.data(ksd + 14);
    const auto *ksd_15 = buffer.data(ksd + 15);
    const auto *ksd_16 = buffer.data(ksd + 16);
    const auto *ksd_17 = buffer.data(ksd + 17);
    const auto *ksd_18 = buffer.data(ksd + 18);
    const auto *ksd_19 = buffer.data(ksd + 19);
    const auto *ksd_21 = buffer.data(ksd + 21);
    const auto *ksd_23 = buffer.data(ksd + 23);
    const auto *ksd_24 = buffer.data(ksd + 24);
    const auto *ksd_27 = buffer.data(ksd + 27);
    const auto *ksd_28 = buffer.data(ksd + 28);
    const auto *ksd_29 = buffer.data(ksd + 29);
    const auto *ksd_30 = buffer.data(ksd + 30);
    const auto *ksd_32 = buffer.data(ksd + 32);
    const auto *ksd_33 = buffer.data(ksd + 33);
    const auto *ksd_34 = buffer.data(ksd + 34);
    const auto *ksd_35 = buffer.data(ksd + 35);
    const auto *ksd_36 = buffer.data(ksd + 36);
    const auto *ksd_37 = buffer.data(ksd + 37);
    const auto *ksd_39 = buffer.data(ksd + 39);
    const auto *ksd_41 = buffer.data(ksd + 41);
    const auto *ksd_42 = buffer.data(ksd + 42);
    const auto *ksd_45 = buffer.data(ksd + 45);
    const auto *ksd_46 = buffer.data(ksd + 46);
    const auto *ksd_47 = buffer.data(ksd + 47);
    const auto *ksd_48 = buffer.data(ksd + 48);
    const auto *ksd_51 = buffer.data(ksd + 51);
    const auto *ksd_52 = buffer.data(ksd + 52);
    const auto *ksd_53 = buffer.data(ksd + 53);
    const auto *ksd_54 = buffer.data(ksd + 54);
    const auto *ksd_56 = buffer.data(ksd + 56);
    const auto *ksd_57 = buffer.data(ksd + 57);
    const auto *ksd_58 = buffer.data(ksd + 58);
    const auto *ksd_59 = buffer.data(ksd + 59);
    const auto *ksd_60 = buffer.data(ksd + 60);
    const auto *ksd_61 = buffer.data(ksd + 61);
    const auto *ksd_63 = buffer.data(ksd + 63);
    const auto *ksd_65 = buffer.data(ksd + 65);
    const auto *ksd_66 = buffer.data(ksd + 66);
    const auto *ksd_69 = buffer.data(ksd + 69);
    const auto *ksd_70 = buffer.data(ksd + 70);
    const auto *ksd_71 = buffer.data(ksd + 71);
    const auto *ksd_72 = buffer.data(ksd + 72);
    const auto *ksd_75 = buffer.data(ksd + 75);
    const auto *ksd_76 = buffer.data(ksd + 76);
    const auto *ksd_77 = buffer.data(ksd + 77);
    const auto *ksd_78 = buffer.data(ksd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, isd_0, isd_3, ksp0_0, \
                         ksp1_0, ksd_0, ksd_2, ksd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * isd_0[k]
                 + f_1 * ksp0_0[k]
                 - f_2 * ksp1_0[k]
                 + f_3 * pc_x[k] * ksd_0[k];

        t_1[k] = f_3 * pc_y[k] * ksd_0[k];

        t_2[k] = f_3 * pc_z[k] * ksd_0[k];

        t_3[k] = f_0 * isd_3[k]
                 + f_3 * pc_x[k] * ksd_3[k];

        t_4[k] = f_3 * pc_y[k] * ksd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, isd_5, ksp0_1, ksp0_2, \
                         ksp1_1, ksp1_2, ksd_3, ksd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * isd_5[k]
                 + f_3 * pc_x[k] * ksd_5[k];

        t_6[k] = f_1 * ksp0_1[k]
                 - f_2 * ksp1_1[k]
                 + f_3 * pc_y[k] * ksd_3[k];

        t_7[k] = f_3 * pc_z[k] * ksd_3[k];

        t_8[k] = f_3 * pc_y[k] * ksd_5[k];

        t_9[k] = f_1 * ksp0_2[k]
                 - f_2 * ksp1_2[k]
                 + f_3 * pc_z[k] * ksd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, isf0_0, isd_0, \
                         isd_9, isf1_0, ksd_6, ksd_7, ksd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * isf0_0[k]
                  - f_4 * pc_y[k] * isf1_0[k];

        t_11[k] = f_5 * isd_0[k]
                  + f_3 * pc_y[k] * ksd_6[k];

        t_12[k] = f_3 * pc_z[k] * ksd_6[k];

        t_13[k] = f_6 * isd_9[k]
                  + f_3 * pc_x[k] * ksd_9[k];

        t_14[k] = f_3 * pc_z[k] * ksd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, isd_3, isd_5, isd_11, \
                         ksp0_4, ksp1_4, ksd_9, ksd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * isd_11[k]
                  + f_3 * pc_x[k] * ksd_11[k];

        t_16[k] = f_5 * isd_3[k]
                  + f_1 * ksp0_4[k]
                  - f_2 * ksp1_4[k]
                  + f_3 * pc_y[k] * ksd_9[k];

        t_17[k] = f_3 * pc_z[k] * ksd_9[k];

        t_18[k] = f_5 * isd_5[k]
                  + f_3 * pc_y[k] * ksd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, isf0_0, isf0_9, \
                         isd_0, isf1_0, isf1_9, ksd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * isf0_9[k]
                  - f_4 * pc_y[k] * isf1_9[k];

        t_20[k] = pa_z[k] * isf0_0[k]
                  - f_4 * pc_z[k] * isf1_0[k];

        t_21[k] = f_3 * pc_y[k] * ksd_12[k];

        t_22[k] = f_5 * isd_0[k]
                  + f_3 * pc_z[k] * ksd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, isf0_6, isd_15, \
                         isd_17, isf1_6, ksd_14, ksd_15, ksd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * isd_15[k]
                  + f_3 * pc_x[k] * ksd_15[k];

        t_24[k] = f_3 * pc_y[k] * ksd_14[k];

        t_25[k] = f_6 * isd_17[k]
                  + f_3 * pc_x[k] * ksd_17[k];

        t_26[k] = pa_z[k] * isf0_6[k]
                  - f_4 * pc_z[k] * isf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, isd_5, isd_18, ksp0_8, \
                         ksp0_9, ksp1_8, ksp1_9, ksd_16, ksd_17, \
                         ksd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * ksp0_8[k]
                  - f_8 * ksp1_8[k]
                  + f_3 * pc_y[k] * ksd_16[k];

        t_28[k] = f_3 * pc_y[k] * ksd_17[k];

        t_29[k] = f_5 * isd_5[k]
                  + f_1 * ksp0_8[k]
                  - f_2 * ksp1_8[k]
                  + f_3 * pc_z[k] * ksd_17[k];

        t_30[k] = f_9 * isd_18[k]
                  + f_1 * ksp0_9[k]
                  - f_2 * ksp1_9[k]
                  + f_3 * pc_x[k] * ksd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, isd_6, isd_21, \
                         isd_23, ksd_18, ksd_19, ksd_21, ksd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * isd_6[k]
                  + f_3 * pc_y[k] * ksd_18[k];

        t_32[k] = f_3 * pc_z[k] * ksd_18[k];

        t_33[k] = f_9 * isd_21[k]
                  + f_3 * pc_x[k] * ksd_21[k];

        t_34[k] = f_3 * pc_z[k] * ksd_19[k];

        t_35[k] = f_9 * isd_23[k]
                  + f_3 * pc_x[k] * ksd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, isd_9, isd_11, ksp0_10, ksp0_11, \
                         ksp1_10, ksp1_11, ksd_21, ksd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * isd_9[k]
                  + f_1 * ksp0_10[k]
                  - f_2 * ksp1_10[k]
                  + f_3 * pc_y[k] * ksd_21[k];

        t_37[k] = f_3 * pc_z[k] * ksd_21[k];

        t_38[k] = f_10 * isd_11[k]
                  + f_3 * pc_y[k] * ksd_23[k];

        t_39[k] = f_1 * ksp0_11[k]
                  - f_2 * ksp1_11[k]
                  + f_3 * pc_z[k] * ksd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, isf0_20, isd_6, \
                         isd_12, isd_27, isf1_20, ksd_24, ksd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * isf0_20[k]
                  - f_4 * pc_y[k] * isf1_20[k];

        t_41[k] = f_5 * isd_12[k]
                  + f_3 * pc_y[k] * ksd_24[k];

        t_42[k] = f_5 * isd_6[k]
                  + f_3 * pc_z[k] * ksd_24[k];

        t_43[k] = f_9 * isd_27[k]
                  + f_3 * pc_x[k] * ksd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, isf0_16, isd_9, isd_28, \
                         isd_29, isf1_16, ksd_27, ksd_28, ksd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * isd_28[k]
                  + f_3 * pc_x[k] * ksd_28[k];

        t_45[k] = f_9 * isd_29[k]
                  + f_3 * pc_x[k] * ksd_29[k];

        t_46[k] = pa_z[k] * isf0_16[k]
                  - f_4 * pc_z[k] * isf1_16[k];

        t_47[k] = f_5 * isd_9[k]
                  + f_3 * pc_z[k] * ksd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, isf0_29, isd_17, isd_30, \
                         isf1_29, ksp0_15, ksp1_15, ksd_29, ksd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * isd_17[k]
                  + f_3 * pc_y[k] * ksd_29[k];

        t_49[k] = pa_y[k] * isf0_29[k]
                  - f_4 * pc_y[k] * isf1_29[k];

        t_50[k] = f_9 * isd_30[k]
                  + f_1 * ksp0_15[k]
                  - f_2 * ksp1_15[k]
                  + f_3 * pc_x[k] * ksd_30[k];

        t_51[k] = f_3 * pc_y[k] * ksd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, isd_12, isd_33, isd_35, \
                         ksd_30, ksd_32, ksd_33, ksd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * isd_12[k]
                  + f_3 * pc_z[k] * ksd_30[k];

        t_53[k] = f_9 * isd_33[k]
                  + f_3 * pc_x[k] * ksd_33[k];

        t_54[k] = f_3 * pc_y[k] * ksd_32[k];

        t_55[k] = f_9 * isd_35[k]
                  + f_3 * pc_x[k] * ksd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, isd_17, ksp0_16, ksp0_17, \
                         ksp1_16, ksp1_17, ksd_33, ksd_34, ksd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * ksp0_16[k]
                  - f_2 * ksp1_16[k]
                  + f_3 * pc_y[k] * ksd_33[k];

        t_57[k] = f_7 * ksp0_17[k]
                  - f_8 * ksp1_17[k]
                  + f_3 * pc_y[k] * ksd_34[k];

        t_58[k] = f_3 * pc_y[k] * ksd_35[k];

        t_59[k] = f_10 * isd_17[k]
                  + f_1 * ksp0_17[k]
                  - f_2 * ksp1_17[k]
                  + f_3 * pc_z[k] * ksd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, isd_18, isd_36, \
                         isd_39, ksp0_18, ksp1_18, ksd_36, ksd_37, \
                         ksd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * isd_36[k]
                  + f_1 * ksp0_18[k]
                  - f_2 * ksp1_18[k]
                  + f_3 * pc_x[k] * ksd_36[k];

        t_61[k] = f_12 * isd_18[k]
                  + f_3 * pc_y[k] * ksd_36[k];

        t_62[k] = f_3 * pc_z[k] * ksd_36[k];

        t_63[k] = f_11 * isd_39[k]
                  + f_3 * pc_x[k] * ksd_39[k];

        t_64[k] = f_3 * pc_z[k] * ksd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, isd_21, isd_23, isd_41, \
                         ksp0_19, ksp1_19, ksd_39, ksd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * isd_41[k]
                  + f_3 * pc_x[k] * ksd_41[k];

        t_66[k] = f_12 * isd_21[k]
                  + f_1 * ksp0_19[k]
                  - f_2 * ksp1_19[k]
                  + f_3 * pc_y[k] * ksd_39[k];

        t_67[k] = f_3 * pc_z[k] * ksd_39[k];

        t_68[k] = f_12 * isd_23[k]
                  + f_3 * pc_y[k] * ksd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, isf0_30, isd_18, isd_24, \
                         isf1_30, ksp0_20, ksp1_20, ksd_41, ksd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * ksp0_20[k]
                  - f_2 * ksp1_20[k]
                  + f_3 * pc_z[k] * ksd_41[k];

        t_70[k] = pa_z[k] * isf0_30[k]
                  - f_4 * pc_z[k] * isf1_30[k];

        t_71[k] = f_10 * isd_24[k]
                  + f_3 * pc_y[k] * ksd_42[k];

        t_72[k] = f_5 * isd_18[k]
                  + f_3 * pc_z[k] * ksd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, isf0_36, isd_45, isd_46, \
                         isd_47, isf1_36, ksd_45, ksd_46, ksd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * isd_45[k]
                  + f_3 * pc_x[k] * ksd_45[k];

        t_74[k] = f_11 * isd_46[k]
                  + f_3 * pc_x[k] * ksd_46[k];

        t_75[k] = f_11 * isd_47[k]
                  + f_3 * pc_x[k] * ksd_47[k];

        t_76[k] = pa_z[k] * isf0_36[k]
                  - f_4 * pc_z[k] * isf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, isf0_50, isd_21, isd_23, \
                         isd_29, isf1_50, ksp0_23, ksp1_23, ksd_45, \
                         ksd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * isd_21[k]
                  + f_3 * pc_z[k] * ksd_45[k];

        t_78[k] = f_10 * isd_29[k]
                  + f_3 * pc_y[k] * ksd_47[k];

        t_79[k] = f_5 * isd_23[k]
                  + f_1 * ksp0_23[k]
                  - f_2 * ksp1_23[k]
                  + f_3 * pc_z[k] * ksd_47[k];

        t_80[k] = pa_y[k] * isf0_50[k]
                  - f_4 * pc_y[k] * isf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, isd_24, isd_30, isd_51, \
                         isd_52, ksd_48, ksd_51, ksd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * isd_30[k]
                  + f_3 * pc_y[k] * ksd_48[k];

        t_82[k] = f_10 * isd_24[k]
                  + f_3 * pc_z[k] * ksd_48[k];

        t_83[k] = f_11 * isd_51[k]
                  + f_3 * pc_x[k] * ksd_51[k];

        t_84[k] = f_11 * isd_52[k]
                  + f_3 * pc_x[k] * ksd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, isd_27, isd_33, isd_35, \
                         isd_53, ksp0_25, ksp1_25, ksd_51, ksd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * isd_53[k]
                  + f_3 * pc_x[k] * ksd_53[k];

        t_86[k] = f_5 * isd_33[k]
                  + f_1 * ksp0_25[k]
                  - f_2 * ksp1_25[k]
                  + f_3 * pc_y[k] * ksd_51[k];

        t_87[k] = f_10 * isd_27[k]
                  + f_3 * pc_z[k] * ksd_51[k];

        t_88[k] = f_5 * isd_35[k]
                  + f_3 * pc_y[k] * ksd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, isf0_59, isd_30, \
                         isd_54, isf1_59, ksp0_27, ksp1_27, ksd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * isf0_59[k]
                  - f_4 * pc_y[k] * isf1_59[k];

        t_90[k] = f_11 * isd_54[k]
                  + f_1 * ksp0_27[k]
                  - f_2 * ksp1_27[k]
                  + f_3 * pc_x[k] * ksd_54[k];

        t_91[k] = f_3 * pc_y[k] * ksd_54[k];

        t_92[k] = f_12 * isd_30[k]
                  + f_3 * pc_z[k] * ksd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, isd_57, isd_59, ksp0_28, ksp1_28, \
                         ksd_56, ksd_57, ksd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * isd_57[k]
                  + f_3 * pc_x[k] * ksd_57[k];

        t_94[k] = f_3 * pc_y[k] * ksd_56[k];

        t_95[k] = f_11 * isd_59[k]
                  + f_3 * pc_x[k] * ksd_59[k];

        t_96[k] = f_1 * ksp0_28[k]
                  - f_2 * ksp1_28[k]
                  + f_3 * pc_y[k] * ksd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, isd_35, isd_60, ksp0_29, \
                         ksp0_30, ksp1_29, ksp1_30, ksd_58, ksd_59, \
                         ksd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * ksp0_29[k]
                  - f_8 * ksp1_29[k]
                  + f_3 * pc_y[k] * ksd_58[k];

        t_98[k] = f_3 * pc_y[k] * ksd_59[k];

        t_99[k] = f_12 * isd_35[k]
                  + f_1 * ksp0_29[k]
                  - f_2 * ksp1_29[k]
                  + f_3 * pc_z[k] * ksd_59[k];

        t_100[k] = f_12 * isd_60[k]
                   + f_1 * ksp0_30[k]
                   - f_2 * ksp1_30[k]
                   + f_3 * pc_x[k] * ksd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, isd_36, isd_63, \
                         isd_65, ksd_60, ksd_61, ksd_63, ksd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_11 * isd_36[k]
                   + f_3 * pc_y[k] * ksd_60[k];

        t_102[k] = f_3 * pc_z[k] * ksd_60[k];

        t_103[k] = f_12 * isd_63[k]
                   + f_3 * pc_x[k] * ksd_63[k];

        t_104[k] = f_3 * pc_z[k] * ksd_61[k];

        t_105[k] = f_12 * isd_65[k]
                   + f_3 * pc_x[k] * ksd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, isd_39, isd_41, ksp0_31, \
                         ksp0_32, ksp1_31, ksp1_32, ksd_63, ksd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_11 * isd_39[k]
                   + f_1 * ksp0_31[k]
                   - f_2 * ksp1_31[k]
                   + f_3 * pc_y[k] * ksd_63[k];

        t_107[k] = f_3 * pc_z[k] * ksd_63[k];

        t_108[k] = f_11 * isd_41[k]
                   + f_3 * pc_y[k] * ksd_65[k];

        t_109[k] = f_1 * ksp0_32[k]
                   - f_2 * ksp1_32[k]
                   + f_3 * pc_z[k] * ksd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, isf0_60, isd_36, \
                         isd_42, isd_69, isf1_60, ksd_66, ksd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * isf0_60[k]
                   - f_4 * pc_z[k] * isf1_60[k];

        t_111[k] = f_12 * isd_42[k]
                   + f_3 * pc_y[k] * ksd_66[k];

        t_112[k] = f_5 * isd_36[k]
                   + f_3 * pc_z[k] * ksd_66[k];

        t_113[k] = f_12 * isd_69[k]
                   + f_3 * pc_x[k] * ksd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, isf0_66, isd_39, \
                         isd_70, isd_71, isf1_66, ksd_69, ksd_70, \
                         ksd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_12 * isd_70[k]
                   + f_3 * pc_x[k] * ksd_70[k];

        t_115[k] = f_12 * isd_71[k]
                   + f_3 * pc_x[k] * ksd_71[k];

        t_116[k] = pa_z[k] * isf0_66[k]
                   - f_4 * pc_z[k] * isf1_66[k];

        t_117[k] = f_5 * isd_39[k]
                   + f_3 * pc_z[k] * ksd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, isd_41, isd_47, isd_72, \
                         ksp0_35, ksp0_36, ksp1_35, ksp1_36, ksd_71, \
                         ksd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * isd_47[k]
                   + f_3 * pc_y[k] * ksd_71[k];

        t_119[k] = f_5 * isd_41[k]
                   + f_1 * ksp0_35[k]
                   - f_2 * ksp1_35[k]
                   + f_3 * pc_z[k] * ksd_71[k];

        t_120[k] = f_12 * isd_72[k]
                   + f_1 * ksp0_36[k]
                   - f_2 * ksp1_36[k]
                   + f_3 * pc_x[k] * ksd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, isd_42, isd_48, isd_75, \
                         isd_76, ksd_72, ksd_75, ksd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * isd_48[k]
                   + f_3 * pc_y[k] * ksd_72[k];

        t_122[k] = f_10 * isd_42[k]
                   + f_3 * pc_z[k] * ksd_72[k];

        t_123[k] = f_12 * isd_75[k]
                   + f_3 * pc_x[k] * ksd_75[k];

        t_124[k] = f_12 * isd_76[k]
                   + f_3 * pc_x[k] * ksd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, isd_45, isd_51, isd_53, \
                         isd_77, ksp0_37, ksp1_37, ksd_75, ksd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * isd_77[k]
                   + f_3 * pc_x[k] * ksd_77[k];

        t_126[k] = f_10 * isd_51[k]
                   + f_1 * ksp0_37[k]
                   - f_2 * ksp1_37[k]
                   + f_3 * pc_y[k] * ksd_75[k];

        t_127[k] = f_10 * isd_45[k]
                   + f_3 * pc_z[k] * ksd_75[k];

        t_128[k] = f_10 * isd_53[k]
                   + f_3 * pc_y[k] * ksd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, isf0_90, isd_47, \
                         isd_48, isd_54, isf1_90, ksp0_38, ksp1_38, ksd_77, \
                         ksd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * isd_47[k]
                   + f_1 * ksp0_38[k]
                   - f_2 * ksp1_38[k]
                   + f_3 * pc_z[k] * ksd_77[k];

        t_130[k] = pa_y[k] * isf0_90[k]
                   - f_4 * pc_y[k] * isf1_90[k];

        t_131[k] = f_5 * isd_54[k]
                   + f_3 * pc_y[k] * ksd_78[k];

        t_132[k] = f_12 * isd_48[k]
                   + f_3 * pc_z[k] * ksd_78[k];
    }
}

static auto
compute_prim_ksf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isf0,
                                                          const size_t isd, const size_t isf1,
                                                          const size_t ksp0, const size_t ksp1,
                                                          const size_t ksd, const size_t ncols,
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
    const auto f_6 = 3.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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
    auto *t_259 = buffer.data(target + 259);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isf0_99 = buffer.data(isf0 + 99);
    const auto *isf0_100 = buffer.data(isf0 + 100);
    const auto *isf0_106 = buffer.data(isf0 + 106);
    const auto *isf0_140 = buffer.data(isf0 + 140);
    const auto *isf0_149 = buffer.data(isf0 + 149);
    const auto *isf0_150 = buffer.data(isf0 + 150);
    const auto *isf0_210 = buffer.data(isf0 + 210);
    const auto *isf0_216 = buffer.data(isf0 + 216);
    const auto *isf0_219 = buffer.data(isf0 + 219);
    const auto *isf0_226 = buffer.data(isf0 + 226);
    const auto *isf0_229 = buffer.data(isf0 + 229);
    const auto *isf0_230 = buffer.data(isf0 + 230);
    const auto *isf0_236 = buffer.data(isf0 + 236);
    const auto *isf0_239 = buffer.data(isf0 + 239);
    const auto *isf0_240 = buffer.data(isf0 + 240);
    const auto *isf0_246 = buffer.data(isf0 + 246);
    const auto *isf0_249 = buffer.data(isf0 + 249);
    const auto *isf0_250 = buffer.data(isf0 + 250);
    const auto *isf0_256 = buffer.data(isf0 + 256);
    const auto *isf0_259 = buffer.data(isf0 + 259);

    const auto *isd_51 = buffer.data(isd + 51);
    const auto *isd_54 = buffer.data(isd + 54);
    const auto *isd_57 = buffer.data(isd + 57);
    const auto *isd_59 = buffer.data(isd + 59);
    const auto *isd_60 = buffer.data(isd + 60);
    const auto *isd_63 = buffer.data(isd + 63);
    const auto *isd_65 = buffer.data(isd + 65);
    const auto *isd_66 = buffer.data(isd + 66);
    const auto *isd_69 = buffer.data(isd + 69);
    const auto *isd_71 = buffer.data(isd + 71);
    const auto *isd_72 = buffer.data(isd + 72);
    const auto *isd_75 = buffer.data(isd + 75);
    const auto *isd_77 = buffer.data(isd + 77);
    const auto *isd_78 = buffer.data(isd + 78);
    const auto *isd_81 = buffer.data(isd + 81);
    const auto *isd_82 = buffer.data(isd + 82);
    const auto *isd_83 = buffer.data(isd + 83);
    const auto *isd_84 = buffer.data(isd + 84);
    const auto *isd_87 = buffer.data(isd + 87);
    const auto *isd_89 = buffer.data(isd + 89);
    const auto *isd_90 = buffer.data(isd + 90);
    const auto *isd_93 = buffer.data(isd + 93);
    const auto *isd_95 = buffer.data(isd + 95);
    const auto *isd_96 = buffer.data(isd + 96);
    const auto *isd_99 = buffer.data(isd + 99);
    const auto *isd_100 = buffer.data(isd + 100);
    const auto *isd_101 = buffer.data(isd + 101);
    const auto *isd_102 = buffer.data(isd + 102);
    const auto *isd_105 = buffer.data(isd + 105);
    const auto *isd_106 = buffer.data(isd + 106);
    const auto *isd_107 = buffer.data(isd + 107);
    const auto *isd_108 = buffer.data(isd + 108);
    const auto *isd_111 = buffer.data(isd + 111);
    const auto *isd_112 = buffer.data(isd + 112);
    const auto *isd_113 = buffer.data(isd + 113);
    const auto *isd_114 = buffer.data(isd + 114);
    const auto *isd_117 = buffer.data(isd + 117);
    const auto *isd_118 = buffer.data(isd + 118);
    const auto *isd_119 = buffer.data(isd + 119);
    const auto *isd_120 = buffer.data(isd + 120);
    const auto *isd_123 = buffer.data(isd + 123);
    const auto *isd_125 = buffer.data(isd + 125);
    const auto *isd_126 = buffer.data(isd + 126);
    const auto *isd_129 = buffer.data(isd + 129);
    const auto *isd_131 = buffer.data(isd + 131);
    const auto *isd_135 = buffer.data(isd + 135);
    const auto *isd_136 = buffer.data(isd + 136);
    const auto *isd_137 = buffer.data(isd + 137);
    const auto *isd_138 = buffer.data(isd + 138);
    const auto *isd_141 = buffer.data(isd + 141);
    const auto *isd_142 = buffer.data(isd + 142);
    const auto *isd_143 = buffer.data(isd + 143);
    const auto *isd_144 = buffer.data(isd + 144);
    const auto *isd_147 = buffer.data(isd + 147);
    const auto *isd_148 = buffer.data(isd + 148);
    const auto *isd_149 = buffer.data(isd + 149);
    const auto *isd_150 = buffer.data(isd + 150);
    const auto *isd_153 = buffer.data(isd + 153);
    const auto *isd_154 = buffer.data(isd + 154);
    const auto *isd_155 = buffer.data(isd + 155);

    const auto *isf1_99 = buffer.data(isf1 + 99);
    const auto *isf1_100 = buffer.data(isf1 + 100);
    const auto *isf1_106 = buffer.data(isf1 + 106);
    const auto *isf1_140 = buffer.data(isf1 + 140);
    const auto *isf1_149 = buffer.data(isf1 + 149);
    const auto *isf1_150 = buffer.data(isf1 + 150);
    const auto *isf1_210 = buffer.data(isf1 + 210);
    const auto *isf1_216 = buffer.data(isf1 + 216);
    const auto *isf1_219 = buffer.data(isf1 + 219);
    const auto *isf1_226 = buffer.data(isf1 + 226);
    const auto *isf1_229 = buffer.data(isf1 + 229);
    const auto *isf1_230 = buffer.data(isf1 + 230);
    const auto *isf1_236 = buffer.data(isf1 + 236);
    const auto *isf1_239 = buffer.data(isf1 + 239);
    const auto *isf1_240 = buffer.data(isf1 + 240);
    const auto *isf1_246 = buffer.data(isf1 + 246);
    const auto *isf1_249 = buffer.data(isf1 + 249);
    const auto *isf1_250 = buffer.data(isf1 + 250);
    const auto *isf1_256 = buffer.data(isf1 + 256);
    const auto *isf1_259 = buffer.data(isf1 + 259);

    const auto *ksp0_40 = buffer.data(ksp0 + 40);
    const auto *ksp0_42 = buffer.data(ksp0 + 42);
    const auto *ksp0_43 = buffer.data(ksp0 + 43);
    const auto *ksp0_44 = buffer.data(ksp0 + 44);
    const auto *ksp0_45 = buffer.data(ksp0 + 45);
    const auto *ksp0_46 = buffer.data(ksp0 + 46);
    const auto *ksp0_47 = buffer.data(ksp0 + 47);
    const auto *ksp0_50 = buffer.data(ksp0 + 50);
    const auto *ksp0_51 = buffer.data(ksp0 + 51);
    const auto *ksp0_52 = buffer.data(ksp0 + 52);
    const auto *ksp0_53 = buffer.data(ksp0 + 53);
    const auto *ksp0_54 = buffer.data(ksp0 + 54);
    const auto *ksp0_55 = buffer.data(ksp0 + 55);
    const auto *ksp0_56 = buffer.data(ksp0 + 56);
    const auto *ksp0_58 = buffer.data(ksp0 + 58);
    const auto *ksp0_60 = buffer.data(ksp0 + 60);
    const auto *ksp0_61 = buffer.data(ksp0 + 61);
    const auto *ksp0_62 = buffer.data(ksp0 + 62);

    const auto *ksp1_40 = buffer.data(ksp1 + 40);
    const auto *ksp1_42 = buffer.data(ksp1 + 42);
    const auto *ksp1_43 = buffer.data(ksp1 + 43);
    const auto *ksp1_44 = buffer.data(ksp1 + 44);
    const auto *ksp1_45 = buffer.data(ksp1 + 45);
    const auto *ksp1_46 = buffer.data(ksp1 + 46);
    const auto *ksp1_47 = buffer.data(ksp1 + 47);
    const auto *ksp1_50 = buffer.data(ksp1 + 50);
    const auto *ksp1_51 = buffer.data(ksp1 + 51);
    const auto *ksp1_52 = buffer.data(ksp1 + 52);
    const auto *ksp1_53 = buffer.data(ksp1 + 53);
    const auto *ksp1_54 = buffer.data(ksp1 + 54);
    const auto *ksp1_55 = buffer.data(ksp1 + 55);
    const auto *ksp1_56 = buffer.data(ksp1 + 56);
    const auto *ksp1_58 = buffer.data(ksp1 + 58);
    const auto *ksp1_60 = buffer.data(ksp1 + 60);
    const auto *ksp1_61 = buffer.data(ksp1 + 61);
    const auto *ksp1_62 = buffer.data(ksp1 + 62);

    const auto *ksd_81 = buffer.data(ksd + 81);
    const auto *ksd_82 = buffer.data(ksd + 82);
    const auto *ksd_83 = buffer.data(ksd + 83);
    const auto *ksd_84 = buffer.data(ksd + 84);
    const auto *ksd_86 = buffer.data(ksd + 86);
    const auto *ksd_87 = buffer.data(ksd + 87);
    const auto *ksd_88 = buffer.data(ksd + 88);
    const auto *ksd_89 = buffer.data(ksd + 89);
    const auto *ksd_90 = buffer.data(ksd + 90);
    const auto *ksd_91 = buffer.data(ksd + 91);
    const auto *ksd_93 = buffer.data(ksd + 93);
    const auto *ksd_95 = buffer.data(ksd + 95);
    const auto *ksd_96 = buffer.data(ksd + 96);
    const auto *ksd_99 = buffer.data(ksd + 99);
    const auto *ksd_100 = buffer.data(ksd + 100);
    const auto *ksd_101 = buffer.data(ksd + 101);
    const auto *ksd_102 = buffer.data(ksd + 102);
    const auto *ksd_105 = buffer.data(ksd + 105);
    const auto *ksd_106 = buffer.data(ksd + 106);
    const auto *ksd_107 = buffer.data(ksd + 107);
    const auto *ksd_108 = buffer.data(ksd + 108);
    const auto *ksd_111 = buffer.data(ksd + 111);
    const auto *ksd_112 = buffer.data(ksd + 112);
    const auto *ksd_113 = buffer.data(ksd + 113);
    const auto *ksd_114 = buffer.data(ksd + 114);
    const auto *ksd_117 = buffer.data(ksd + 117);
    const auto *ksd_118 = buffer.data(ksd + 118);
    const auto *ksd_119 = buffer.data(ksd + 119);
    const auto *ksd_120 = buffer.data(ksd + 120);
    const auto *ksd_122 = buffer.data(ksd + 122);
    const auto *ksd_123 = buffer.data(ksd + 123);
    const auto *ksd_124 = buffer.data(ksd + 124);
    const auto *ksd_125 = buffer.data(ksd + 125);
    const auto *ksd_126 = buffer.data(ksd + 126);
    const auto *ksd_127 = buffer.data(ksd + 127);
    const auto *ksd_129 = buffer.data(ksd + 129);
    const auto *ksd_131 = buffer.data(ksd + 131);
    const auto *ksd_132 = buffer.data(ksd + 132);
    const auto *ksd_135 = buffer.data(ksd + 135);
    const auto *ksd_136 = buffer.data(ksd + 136);
    const auto *ksd_137 = buffer.data(ksd + 137);
    const auto *ksd_138 = buffer.data(ksd + 138);
    const auto *ksd_141 = buffer.data(ksd + 141);
    const auto *ksd_142 = buffer.data(ksd + 142);
    const auto *ksd_143 = buffer.data(ksd + 143);
    const auto *ksd_144 = buffer.data(ksd + 144);
    const auto *ksd_147 = buffer.data(ksd + 147);
    const auto *ksd_148 = buffer.data(ksd + 148);
    const auto *ksd_149 = buffer.data(ksd + 149);
    const auto *ksd_150 = buffer.data(ksd + 150);
    const auto *ksd_153 = buffer.data(ksd + 153);
    const auto *ksd_154 = buffer.data(ksd + 154);
    const auto *ksd_155 = buffer.data(ksd + 155);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, isd_57, isd_81, isd_82, \
                         isd_83, ksp0_40, ksp1_40, ksd_81, ksd_82, \
                         ksd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_12 * isd_81[k]
                   + f_3 * pc_x[k] * ksd_81[k];

        t_134[k] = f_12 * isd_82[k]
                   + f_3 * pc_x[k] * ksd_82[k];

        t_135[k] = f_12 * isd_83[k]
                   + f_3 * pc_x[k] * ksd_83[k];

        t_136[k] = f_5 * isd_57[k]
                   + f_1 * ksp0_40[k]
                   - f_2 * ksp1_40[k]
                   + f_3 * pc_y[k] * ksd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, isf0_99, isd_51, isd_59, \
                         isf1_99, ksd_81, ksd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * isd_51[k]
                   + f_3 * pc_z[k] * ksd_81[k];

        t_138[k] = f_5 * isd_59[k]
                   + f_3 * pc_y[k] * ksd_83[k];

        t_139[k] = pa_y[k] * isf0_99[k]
                   - f_4 * pc_y[k] * isf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, isd_54, isd_84, \
                         isd_87, ksp0_42, ksp1_42, ksd_84, ksd_86, \
                         ksd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_12 * isd_84[k]
                   + f_1 * ksp0_42[k]
                   - f_2 * ksp1_42[k]
                   + f_3 * pc_x[k] * ksd_84[k];

        t_141[k] = f_3 * pc_y[k] * ksd_84[k];

        t_142[k] = f_11 * isd_54[k]
                   + f_3 * pc_z[k] * ksd_84[k];

        t_143[k] = f_12 * isd_87[k]
                   + f_3 * pc_x[k] * ksd_87[k];

        t_144[k] = f_3 * pc_y[k] * ksd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, isd_89, ksp0_43, ksp0_44, \
                         ksp1_43, ksp1_44, ksd_87, ksd_88, ksd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_12 * isd_89[k]
                   + f_3 * pc_x[k] * ksd_89[k];

        t_146[k] = f_1 * ksp0_43[k]
                   - f_2 * ksp1_43[k]
                   + f_3 * pc_y[k] * ksd_87[k];

        t_147[k] = f_7 * ksp0_44[k]
                   - f_8 * ksp1_44[k]
                   + f_3 * pc_y[k] * ksd_88[k];

        t_148[k] = f_3 * pc_y[k] * ksd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_y, pc_z, isd_59, isd_60, isd_90, \
                         ksp0_44, ksp0_45, ksp1_44, ksp1_45, ksd_89, \
                         ksd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_11 * isd_59[k]
                   + f_1 * ksp0_44[k]
                   - f_2 * ksp1_44[k]
                   + f_3 * pc_z[k] * ksd_89[k];

        t_150[k] = f_10 * isd_90[k]
                   + f_1 * ksp0_45[k]
                   - f_2 * ksp1_45[k]
                   + f_3 * pc_x[k] * ksd_90[k];

        t_151[k] = f_9 * isd_60[k]
                   + f_3 * pc_y[k] * ksd_90[k];

        t_152[k] = f_3 * pc_z[k] * ksd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, isd_63, isd_93, \
                         isd_95, ksp0_46, ksp1_46, ksd_91, ksd_93, \
                         ksd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_10 * isd_93[k]
                   + f_3 * pc_x[k] * ksd_93[k];

        t_154[k] = f_3 * pc_z[k] * ksd_91[k];

        t_155[k] = f_10 * isd_95[k]
                   + f_3 * pc_x[k] * ksd_95[k];

        t_156[k] = f_9 * isd_63[k]
                   + f_1 * ksp0_46[k]
                   - f_2 * ksp1_46[k]
                   + f_3 * pc_y[k] * ksd_93[k];

        t_157[k] = f_3 * pc_z[k] * ksd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pc_y, pc_z, isf0_100, isd_65, \
                         isd_66, isf1_100, ksp0_47, ksp1_47, ksd_95, \
                         ksd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_9 * isd_65[k]
                   + f_3 * pc_y[k] * ksd_95[k];

        t_159[k] = f_1 * ksp0_47[k]
                   - f_2 * ksp1_47[k]
                   + f_3 * pc_z[k] * ksd_95[k];

        t_160[k] = pa_z[k] * isf0_100[k]
                   - f_4 * pc_z[k] * isf1_100[k];

        t_161[k] = f_11 * isd_66[k]
                   + f_3 * pc_y[k] * ksd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, isd_60, isd_99, isd_100, \
                         isd_101, ksd_96, ksd_99, ksd_100, ksd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * isd_60[k]
                   + f_3 * pc_z[k] * ksd_96[k];

        t_163[k] = f_10 * isd_99[k]
                   + f_3 * pc_x[k] * ksd_99[k];

        t_164[k] = f_10 * isd_100[k]
                   + f_3 * pc_x[k] * ksd_100[k];

        t_165[k] = f_10 * isd_101[k]
                   + f_3 * pc_x[k] * ksd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, isf0_106, isd_63, \
                         isd_65, isd_71, isf1_106, ksp0_50, ksp1_50, ksd_99, \
                         ksd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * isf0_106[k]
                   - f_4 * pc_z[k] * isf1_106[k];

        t_167[k] = f_5 * isd_63[k]
                   + f_3 * pc_z[k] * ksd_99[k];

        t_168[k] = f_11 * isd_71[k]
                   + f_3 * pc_y[k] * ksd_101[k];

        t_169[k] = f_5 * isd_65[k]
                   + f_1 * ksp0_50[k]
                   - f_2 * ksp1_50[k]
                   + f_3 * pc_z[k] * ksd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, isd_66, isd_72, \
                         isd_102, isd_105, ksp0_51, ksp1_51, ksd_102, \
                         ksd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * isd_102[k]
                   + f_1 * ksp0_51[k]
                   - f_2 * ksp1_51[k]
                   + f_3 * pc_x[k] * ksd_102[k];

        t_171[k] = f_12 * isd_72[k]
                   + f_3 * pc_y[k] * ksd_102[k];

        t_172[k] = f_10 * isd_66[k]
                   + f_3 * pc_z[k] * ksd_102[k];

        t_173[k] = f_10 * isd_105[k]
                   + f_3 * pc_x[k] * ksd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, isd_69, isd_75, \
                         isd_106, isd_107, ksp0_52, ksp1_52, ksd_105, ksd_106, \
                         ksd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_10 * isd_106[k]
                   + f_3 * pc_x[k] * ksd_106[k];

        t_175[k] = f_10 * isd_107[k]
                   + f_3 * pc_x[k] * ksd_107[k];

        t_176[k] = f_12 * isd_75[k]
                   + f_1 * ksp0_52[k]
                   - f_2 * ksp1_52[k]
                   + f_3 * pc_y[k] * ksd_105[k];

        t_177[k] = f_10 * isd_69[k]
                   + f_3 * pc_z[k] * ksd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, isd_71, isd_77, isd_108, \
                         ksp0_53, ksp0_54, ksp1_53, ksp1_54, ksd_107, \
                         ksd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * isd_77[k]
                   + f_3 * pc_y[k] * ksd_107[k];

        t_179[k] = f_10 * isd_71[k]
                   + f_1 * ksp0_53[k]
                   - f_2 * ksp1_53[k]
                   + f_3 * pc_z[k] * ksd_107[k];

        t_180[k] = f_10 * isd_108[k]
                   + f_1 * ksp0_54[k]
                   - f_2 * ksp1_54[k]
                   + f_3 * pc_x[k] * ksd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, isd_72, isd_78, \
                         isd_111, isd_112, ksd_108, ksd_111, ksd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * isd_78[k]
                   + f_3 * pc_y[k] * ksd_108[k];

        t_182[k] = f_12 * isd_72[k]
                   + f_3 * pc_z[k] * ksd_108[k];

        t_183[k] = f_10 * isd_111[k]
                   + f_3 * pc_x[k] * ksd_111[k];

        t_184[k] = f_10 * isd_112[k]
                   + f_3 * pc_x[k] * ksd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, isd_75, isd_81, isd_83, \
                         isd_113, ksp0_55, ksp1_55, ksd_111, ksd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_10 * isd_113[k]
                   + f_3 * pc_x[k] * ksd_113[k];

        t_186[k] = f_10 * isd_81[k]
                   + f_1 * ksp0_55[k]
                   - f_2 * ksp1_55[k]
                   + f_3 * pc_y[k] * ksd_111[k];

        t_187[k] = f_12 * isd_75[k]
                   + f_3 * pc_z[k] * ksd_111[k];

        t_188[k] = f_10 * isd_83[k]
                   + f_3 * pc_y[k] * ksd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_y, pc_y, pc_z, isf0_140, isd_77, \
                         isd_78, isd_84, isf1_140, ksp0_56, ksp1_56, ksd_113, \
                         ksd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * isd_77[k]
                   + f_1 * ksp0_56[k]
                   - f_2 * ksp1_56[k]
                   + f_3 * pc_z[k] * ksd_113[k];

        t_190[k] = pa_y[k] * isf0_140[k]
                   - f_4 * pc_y[k] * isf1_140[k];

        t_191[k] = f_5 * isd_84[k]
                   + f_3 * pc_y[k] * ksd_114[k];

        t_192[k] = f_11 * isd_78[k]
                   + f_3 * pc_z[k] * ksd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, isd_87, isd_117, isd_118, \
                         isd_119, ksp0_58, ksp1_58, ksd_117, ksd_118, \
                         ksd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_10 * isd_117[k]
                   + f_3 * pc_x[k] * ksd_117[k];

        t_194[k] = f_10 * isd_118[k]
                   + f_3 * pc_x[k] * ksd_118[k];

        t_195[k] = f_10 * isd_119[k]
                   + f_3 * pc_x[k] * ksd_119[k];

        t_196[k] = f_5 * isd_87[k]
                   + f_1 * ksp0_58[k]
                   - f_2 * ksp1_58[k]
                   + f_3 * pc_y[k] * ksd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, isf0_149, isd_81, isd_89, \
                         isf1_149, ksd_117, ksd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_11 * isd_81[k]
                   + f_3 * pc_z[k] * ksd_117[k];

        t_198[k] = f_5 * isd_89[k]
                   + f_3 * pc_y[k] * ksd_119[k];

        t_199[k] = pa_y[k] * isf0_149[k]
                   - f_4 * pc_y[k] * isf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, isd_84, isd_120, \
                         isd_123, ksp0_60, ksp1_60, ksd_120, ksd_122, \
                         ksd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_10 * isd_120[k]
                   + f_1 * ksp0_60[k]
                   - f_2 * ksp1_60[k]
                   + f_3 * pc_x[k] * ksd_120[k];

        t_201[k] = f_3 * pc_y[k] * ksd_120[k];

        t_202[k] = f_9 * isd_84[k]
                   + f_3 * pc_z[k] * ksd_120[k];

        t_203[k] = f_10 * isd_123[k]
                   + f_3 * pc_x[k] * ksd_123[k];

        t_204[k] = f_3 * pc_y[k] * ksd_122[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, isd_125, ksp0_61, ksp0_62, \
                         ksp1_61, ksp1_62, ksd_123, ksd_124, ksd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_10 * isd_125[k]
                   + f_3 * pc_x[k] * ksd_125[k];

        t_206[k] = f_1 * ksp0_61[k]
                   - f_2 * ksp1_61[k]
                   + f_3 * pc_y[k] * ksd_123[k];

        t_207[k] = f_7 * ksp0_62[k]
                   - f_8 * ksp1_62[k]
                   + f_3 * pc_y[k] * ksd_124[k];

        t_208[k] = f_3 * pc_y[k] * ksd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pa_x, pc_x, pc_y, pc_z, isf0_210, isd_89, \
                         isd_90, isd_126, isf1_210, ksp0_62, ksp1_62, ksd_125, \
                         ksd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_9 * isd_89[k]
                   + f_1 * ksp0_62[k]
                   - f_2 * ksp1_62[k]
                   + f_3 * pc_z[k] * ksd_125[k];

        t_210[k] = pa_x[k] * isf0_210[k]
                   + f_12 * isd_126[k]
                   - f_4 * pc_x[k] * isf1_210[k];

        t_211[k] = f_6 * isd_90[k]
                   + f_3 * pc_y[k] * ksd_126[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pc_x, pc_z, isf0_216, \
                         isd_129, isd_131, isf1_216, ksd_126, ksd_127, ksd_129, \
                         ksd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * ksd_126[k];

        t_213[k] = f_5 * isd_129[k]
                   + f_3 * pc_x[k] * ksd_129[k];

        t_214[k] = f_3 * pc_z[k] * ksd_127[k];

        t_215[k] = f_5 * isd_131[k]
                   + f_3 * pc_x[k] * ksd_131[k];

        t_216[k] = pa_x[k] * isf0_216[k]
                   - f_4 * pc_x[k] * isf1_216[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_x, pa_z, pc_x, pc_y, pc_z, isf0_150, \
                         isf0_219, isd_95, isf1_150, isf1_219, ksd_129, \
                         ksd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * pc_z[k] * ksd_129[k];

        t_218[k] = f_6 * isd_95[k]
                   + f_3 * pc_y[k] * ksd_131[k];

        t_219[k] = pa_x[k] * isf0_219[k]
                   - f_4 * pc_x[k] * isf1_219[k];

        t_220[k] = pa_z[k] * isf0_150[k]
                   - f_4 * pc_z[k] * isf1_150[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pc_x, pc_y, pc_z, isd_90, isd_96, \
                         isd_135, isd_136, ksd_132, ksd_135, ksd_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * isd_96[k]
                   + f_3 * pc_y[k] * ksd_132[k];

        t_222[k] = f_5 * isd_90[k]
                   + f_3 * pc_z[k] * ksd_132[k];

        t_223[k] = f_5 * isd_135[k]
                   + f_3 * pc_x[k] * ksd_135[k];

        t_224[k] = f_5 * isd_136[k]
                   + f_3 * pc_x[k] * ksd_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pc_x, pc_y, pc_z, isf0_226, isd_93, \
                         isd_101, isd_137, isf1_226, ksd_135, ksd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_5 * isd_137[k]
                   + f_3 * pc_x[k] * ksd_137[k];

        t_226[k] = pa_x[k] * isf0_226[k]
                   - f_4 * pc_x[k] * isf1_226[k];

        t_227[k] = f_5 * isd_93[k]
                   + f_3 * pc_z[k] * ksd_135[k];

        t_228[k] = f_9 * isd_101[k]
                   + f_3 * pc_y[k] * ksd_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_x, pc_x, pc_y, pc_z, isf0_229, \
                         isf0_230, isd_96, isd_102, isd_138, isf1_229, isf1_230, \
                         ksd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_x[k] * isf0_229[k]
                   - f_4 * pc_x[k] * isf1_229[k];

        t_230[k] = pa_x[k] * isf0_230[k]
                   + f_12 * isd_138[k]
                   - f_4 * pc_x[k] * isf1_230[k];

        t_231[k] = f_11 * isd_102[k]
                   + f_3 * pc_y[k] * ksd_138[k];

        t_232[k] = f_10 * isd_96[k]
                   + f_3 * pc_z[k] * ksd_138[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pc_x, isf0_236, isd_141, isd_142, \
                         isd_143, isf1_236, ksd_141, ksd_142, ksd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_5 * isd_141[k]
                   + f_3 * pc_x[k] * ksd_141[k];

        t_234[k] = f_5 * isd_142[k]
                   + f_3 * pc_x[k] * ksd_142[k];

        t_235[k] = f_5 * isd_143[k]
                   + f_3 * pc_x[k] * ksd_143[k];

        t_236[k] = pa_x[k] * isf0_236[k]
                   - f_4 * pc_x[k] * isf1_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_x, pc_x, pc_y, pc_z, isf0_239, isd_99, \
                         isd_107, isf1_239, ksd_141, ksd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_10 * isd_99[k]
                   + f_3 * pc_z[k] * ksd_141[k];

        t_238[k] = f_11 * isd_107[k]
                   + f_3 * pc_y[k] * ksd_143[k];

        t_239[k] = pa_x[k] * isf0_239[k]
                   - f_4 * pc_x[k] * isf1_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pa_x, pc_x, pc_y, pc_z, isf0_240, \
                         isd_102, isd_108, isd_144, isd_147, isf1_240, ksd_144, \
                         ksd_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_x[k] * isf0_240[k]
                   + f_12 * isd_144[k]
                   - f_4 * pc_x[k] * isf1_240[k];

        t_241[k] = f_12 * isd_108[k]
                   + f_3 * pc_y[k] * ksd_144[k];

        t_242[k] = f_12 * isd_102[k]
                   + f_3 * pc_z[k] * ksd_144[k];

        t_243[k] = f_5 * isd_147[k]
                   + f_3 * pc_x[k] * ksd_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_x, pc_x, pc_z, isf0_246, isd_105, \
                         isd_148, isd_149, isf1_246, ksd_147, ksd_148, \
                         ksd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_5 * isd_148[k]
                   + f_3 * pc_x[k] * ksd_148[k];

        t_245[k] = f_5 * isd_149[k]
                   + f_3 * pc_x[k] * ksd_149[k];

        t_246[k] = pa_x[k] * isf0_246[k]
                   - f_4 * pc_x[k] * isf1_246[k];

        t_247[k] = f_12 * isd_105[k]
                   + f_3 * pc_z[k] * ksd_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_x, pc_x, pc_y, isf0_249, isf0_250, \
                         isd_113, isd_114, isd_150, isf1_249, isf1_250, ksd_149, \
                         ksd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_12 * isd_113[k]
                   + f_3 * pc_y[k] * ksd_149[k];

        t_249[k] = pa_x[k] * isf0_249[k]
                   - f_4 * pc_x[k] * isf1_249[k];

        t_250[k] = pa_x[k] * isf0_250[k]
                   + f_12 * isd_150[k]
                   - f_4 * pc_x[k] * isf1_250[k];

        t_251[k] = f_10 * isd_114[k]
                   + f_3 * pc_y[k] * ksd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, isd_108, isd_153, isd_154, \
                         isd_155, ksd_150, ksd_153, ksd_154, ksd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_11 * isd_108[k]
                   + f_3 * pc_z[k] * ksd_150[k];

        t_253[k] = f_5 * isd_153[k]
                   + f_3 * pc_x[k] * ksd_153[k];

        t_254[k] = f_5 * isd_154[k]
                   + f_3 * pc_x[k] * ksd_154[k];

        t_255[k] = f_5 * isd_155[k]
                   + f_3 * pc_x[k] * ksd_155[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_x, pc_x, pc_y, pc_z, isf0_256, \
                         isf0_259, isd_111, isd_119, isf1_256, isf1_259, ksd_153, \
                         ksd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pa_x[k] * isf0_256[k]
                   - f_4 * pc_x[k] * isf1_256[k];

        t_257[k] = f_11 * isd_111[k]
                   + f_3 * pc_z[k] * ksd_153[k];

        t_258[k] = f_10 * isd_119[k]
                   + f_3 * pc_y[k] * ksd_155[k];

        t_259[k] = pa_x[k] * isf0_259[k]
                   - f_4 * pc_x[k] * isf1_259[k];
    }
}

static auto
compute_prim_ksf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isf0,
                                                          const size_t isd, const size_t isf1,
                                                          const size_t ksp0, const size_t ksp1,
                                                          const size_t ksd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *isf0_200 = buffer.data(isf0 + 200);
    const auto *isf0_210 = buffer.data(isf0 + 210);
    const auto *isf0_211 = buffer.data(isf0 + 211);
    const auto *isf0_216 = buffer.data(isf0 + 216);
    const auto *isf0_266 = buffer.data(isf0 + 266);
    const auto *isf0_269 = buffer.data(isf0 + 269);
    const auto *isf0_270 = buffer.data(isf0 + 270);
    const auto *isf0_272 = buffer.data(isf0 + 272);
    const auto *isf0_276 = buffer.data(isf0 + 276);
    const auto *isf0_277 = buffer.data(isf0 + 277);
    const auto *isf0_279 = buffer.data(isf0 + 279);

    const auto *isd_114 = buffer.data(isd + 114);
    const auto *isd_117 = buffer.data(isd + 117);
    const auto *isd_120 = buffer.data(isd + 120);
    const auto *isd_125 = buffer.data(isd + 125);
    const auto *isd_129 = buffer.data(isd + 129);
    const auto *isd_131 = buffer.data(isd + 131);
    const auto *isd_135 = buffer.data(isd + 135);
    const auto *isd_137 = buffer.data(isd + 137);
    const auto *isd_141 = buffer.data(isd + 141);
    const auto *isd_143 = buffer.data(isd + 143);
    const auto *isd_147 = buffer.data(isd + 147);
    const auto *isd_149 = buffer.data(isd + 149);
    const auto *isd_153 = buffer.data(isd + 153);
    const auto *isd_155 = buffer.data(isd + 155);
    const auto *isd_159 = buffer.data(isd + 159);
    const auto *isd_160 = buffer.data(isd + 160);
    const auto *isd_161 = buffer.data(isd + 161);
    const auto *isd_162 = buffer.data(isd + 162);
    const auto *isd_165 = buffer.data(isd + 165);
    const auto *isd_167 = buffer.data(isd + 167);

    const auto *isf1_200 = buffer.data(isf1 + 200);
    const auto *isf1_210 = buffer.data(isf1 + 210);
    const auto *isf1_211 = buffer.data(isf1 + 211);
    const auto *isf1_216 = buffer.data(isf1 + 216);
    const auto *isf1_266 = buffer.data(isf1 + 266);
    const auto *isf1_269 = buffer.data(isf1 + 269);
    const auto *isf1_270 = buffer.data(isf1 + 270);
    const auto *isf1_272 = buffer.data(isf1 + 272);
    const auto *isf1_276 = buffer.data(isf1 + 276);
    const auto *isf1_277 = buffer.data(isf1 + 277);
    const auto *isf1_279 = buffer.data(isf1 + 279);

    const auto *ksp0_84 = buffer.data(ksp0 + 84);
    const auto *ksp0_85 = buffer.data(ksp0 + 85);
    const auto *ksp0_86 = buffer.data(ksp0 + 86);
    const auto *ksp0_89 = buffer.data(ksp0 + 89);
    const auto *ksp0_90 = buffer.data(ksp0 + 90);
    const auto *ksp0_91 = buffer.data(ksp0 + 91);
    const auto *ksp0_92 = buffer.data(ksp0 + 92);
    const auto *ksp0_93 = buffer.data(ksp0 + 93);
    const auto *ksp0_94 = buffer.data(ksp0 + 94);
    const auto *ksp0_95 = buffer.data(ksp0 + 95);
    const auto *ksp0_96 = buffer.data(ksp0 + 96);
    const auto *ksp0_97 = buffer.data(ksp0 + 97);
    const auto *ksp0_98 = buffer.data(ksp0 + 98);
    const auto *ksp0_99 = buffer.data(ksp0 + 99);
    const auto *ksp0_100 = buffer.data(ksp0 + 100);
    const auto *ksp0_101 = buffer.data(ksp0 + 101);
    const auto *ksp0_103 = buffer.data(ksp0 + 103);
    const auto *ksp0_105 = buffer.data(ksp0 + 105);
    const auto *ksp0_106 = buffer.data(ksp0 + 106);
    const auto *ksp0_107 = buffer.data(ksp0 + 107);

    const auto *ksp1_84 = buffer.data(ksp1 + 84);
    const auto *ksp1_85 = buffer.data(ksp1 + 85);
    const auto *ksp1_86 = buffer.data(ksp1 + 86);
    const auto *ksp1_89 = buffer.data(ksp1 + 89);
    const auto *ksp1_90 = buffer.data(ksp1 + 90);
    const auto *ksp1_91 = buffer.data(ksp1 + 91);
    const auto *ksp1_92 = buffer.data(ksp1 + 92);
    const auto *ksp1_93 = buffer.data(ksp1 + 93);
    const auto *ksp1_94 = buffer.data(ksp1 + 94);
    const auto *ksp1_95 = buffer.data(ksp1 + 95);
    const auto *ksp1_96 = buffer.data(ksp1 + 96);
    const auto *ksp1_97 = buffer.data(ksp1 + 97);
    const auto *ksp1_98 = buffer.data(ksp1 + 98);
    const auto *ksp1_99 = buffer.data(ksp1 + 99);
    const auto *ksp1_100 = buffer.data(ksp1 + 100);
    const auto *ksp1_101 = buffer.data(ksp1 + 101);
    const auto *ksp1_103 = buffer.data(ksp1 + 103);
    const auto *ksp1_105 = buffer.data(ksp1 + 105);
    const auto *ksp1_106 = buffer.data(ksp1 + 106);
    const auto *ksp1_107 = buffer.data(ksp1 + 107);

    const auto *ksd_156 = buffer.data(ksd + 156);
    const auto *ksd_159 = buffer.data(ksd + 159);
    const auto *ksd_160 = buffer.data(ksd + 160);
    const auto *ksd_161 = buffer.data(ksd + 161);
    const auto *ksd_162 = buffer.data(ksd + 162);
    const auto *ksd_164 = buffer.data(ksd + 164);
    const auto *ksd_165 = buffer.data(ksd + 165);
    const auto *ksd_167 = buffer.data(ksd + 167);
    const auto *ksd_168 = buffer.data(ksd + 168);
    const auto *ksd_169 = buffer.data(ksd + 169);
    const auto *ksd_171 = buffer.data(ksd + 171);
    const auto *ksd_172 = buffer.data(ksd + 172);
    const auto *ksd_173 = buffer.data(ksd + 173);
    const auto *ksd_176 = buffer.data(ksd + 176);
    const auto *ksd_177 = buffer.data(ksd + 177);
    const auto *ksd_178 = buffer.data(ksd + 178);
    const auto *ksd_179 = buffer.data(ksd + 179);
    const auto *ksd_180 = buffer.data(ksd + 180);
    const auto *ksd_181 = buffer.data(ksd + 181);
    const auto *ksd_182 = buffer.data(ksd + 182);
    const auto *ksd_183 = buffer.data(ksd + 183);
    const auto *ksd_184 = buffer.data(ksd + 184);
    const auto *ksd_185 = buffer.data(ksd + 185);
    const auto *ksd_186 = buffer.data(ksd + 186);
    const auto *ksd_187 = buffer.data(ksd + 187);
    const auto *ksd_188 = buffer.data(ksd + 188);
    const auto *ksd_189 = buffer.data(ksd + 189);
    const auto *ksd_190 = buffer.data(ksd + 190);
    const auto *ksd_191 = buffer.data(ksd + 191);
    const auto *ksd_192 = buffer.data(ksd + 192);
    const auto *ksd_193 = buffer.data(ksd + 193);
    const auto *ksd_194 = buffer.data(ksd + 194);
    const auto *ksd_195 = buffer.data(ksd + 195);
    const auto *ksd_196 = buffer.data(ksd + 196);
    const auto *ksd_197 = buffer.data(ksd + 197);
    const auto *ksd_198 = buffer.data(ksd + 198);
    const auto *ksd_199 = buffer.data(ksd + 199);
    const auto *ksd_200 = buffer.data(ksd + 200);
    const auto *ksd_201 = buffer.data(ksd + 201);
    const auto *ksd_202 = buffer.data(ksd + 202);
    const auto *ksd_203 = buffer.data(ksd + 203);
    const auto *ksd_205 = buffer.data(ksd + 205);
    const auto *ksd_207 = buffer.data(ksd + 207);
    const auto *ksd_208 = buffer.data(ksd + 208);
    const auto *ksd_209 = buffer.data(ksd + 209);
    const auto *ksd_210 = buffer.data(ksd + 210);
    const auto *ksd_212 = buffer.data(ksd + 212);
    const auto *ksd_213 = buffer.data(ksd + 213);
    const auto *ksd_214 = buffer.data(ksd + 214);
    const auto *ksd_215 = buffer.data(ksd + 215);

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, pc_z, isf0_200, \
                         isd_114, isd_120, isd_159, isf1_200, ksd_156, \
                         ksd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * isf0_200[k]
                   - f_4 * pc_y[k] * isf1_200[k];

        t_261[k] = f_5 * isd_120[k]
                   + f_3 * pc_y[k] * ksd_156[k];

        t_262[k] = f_9 * isd_114[k]
                   + f_3 * pc_z[k] * ksd_156[k];

        t_263[k] = f_5 * isd_159[k]
                   + f_3 * pc_x[k] * ksd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_x, pc_x, pc_z, isf0_266, isd_117, \
                         isd_160, isd_161, isf1_266, ksd_159, ksd_160, \
                         ksd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_5 * isd_160[k]
                   + f_3 * pc_x[k] * ksd_160[k];

        t_265[k] = f_5 * isd_161[k]
                   + f_3 * pc_x[k] * ksd_161[k];

        t_266[k] = pa_x[k] * isf0_266[k]
                   - f_4 * pc_x[k] * isf1_266[k];

        t_267[k] = f_9 * isd_117[k]
                   + f_3 * pc_z[k] * ksd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_x, pc_x, pc_y, isf0_269, isf0_270, \
                         isd_125, isd_162, isf1_269, isf1_270, ksd_161, \
                         ksd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * isd_125[k]
                   + f_3 * pc_y[k] * ksd_161[k];

        t_269[k] = pa_x[k] * isf0_269[k]
                   - f_4 * pc_x[k] * isf1_269[k];

        t_270[k] = pa_x[k] * isf0_270[k]
                   + f_12 * isd_162[k]
                   - f_4 * pc_x[k] * isf1_270[k];

        t_271[k] = f_3 * pc_y[k] * ksd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, isd_120, isd_165, \
                         isd_167, ksd_162, ksd_164, ksd_165, ksd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_6 * isd_120[k]
                   + f_3 * pc_z[k] * ksd_162[k];

        t_273[k] = f_5 * isd_165[k]
                   + f_3 * pc_x[k] * ksd_165[k];

        t_274[k] = f_3 * pc_y[k] * ksd_164[k];

        t_275[k] = f_5 * isd_167[k]
                   + f_3 * pc_x[k] * ksd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pc_x, pc_y, isf0_276, isf0_277, \
                         isf0_279, isf1_276, isf1_277, isf1_279, \
                         ksd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_x[k] * isf0_276[k]
                   - f_4 * pc_x[k] * isf1_276[k];

        t_277[k] = pa_x[k] * isf0_277[k]
                   - f_4 * pc_x[k] * isf1_277[k];

        t_278[k] = f_3 * pc_y[k] * ksd_167[k];

        t_279[k] = pa_x[k] * isf0_279[k]
                   - f_4 * pc_x[k] * isf1_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pc_x, pc_z, ksp0_84, ksp0_85, \
                         ksp1_84, ksp1_85, ksd_168, ksd_169, ksd_171, \
                         ksd_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_1 * ksp0_84[k]
                   - f_2 * ksp1_84[k]
                   + f_3 * pc_x[k] * ksd_168[k];

        t_281[k] = f_7 * ksp0_85[k]
                   - f_8 * ksp1_85[k]
                   + f_3 * pc_x[k] * ksd_169[k];

        t_282[k] = f_3 * pc_z[k] * ksd_168[k];

        t_283[k] = f_3 * pc_x[k] * ksd_171[k];

        t_284[k] = f_3 * pc_x[k] * ksd_172[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, pc_x, pc_y, pc_z, isd_129, \
                         isd_131, ksp0_85, ksp0_86, ksp1_85, ksp1_86, ksd_171, \
                         ksd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_3 * pc_x[k] * ksd_173[k];

        t_286[k] = f_0 * isd_129[k]
                   + f_1 * ksp0_85[k]
                   - f_2 * ksp1_85[k]
                   + f_3 * pc_y[k] * ksd_171[k];

        t_287[k] = f_3 * pc_z[k] * ksd_171[k];

        t_288[k] = f_0 * isd_131[k]
                   + f_3 * pc_y[k] * ksd_173[k];

        t_289[k] = f_1 * ksp0_86[k]
                   - f_2 * ksp1_86[k]
                   + f_3 * pc_z[k] * ksd_173[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_z, pc_x, pc_z, isf0_210, isf0_211, \
                         isf1_210, isf1_211, ksp0_89, ksp1_89, ksd_176, \
                         ksd_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_z[k] * isf0_210[k]
                   - f_4 * pc_z[k] * isf1_210[k];

        t_291[k] = pa_z[k] * isf0_211[k]
                   - f_4 * pc_z[k] * isf1_211[k];

        t_292[k] = f_7 * ksp0_89[k]
                   - f_8 * ksp1_89[k]
                   + f_3 * pc_x[k] * ksd_176[k];

        t_293[k] = f_3 * pc_x[k] * ksd_177[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pa_z, pc_x, pc_y, pc_z, isf0_216, \
                         isd_129, isd_137, isf1_216, ksd_177, ksd_178, \
                         ksd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_3 * pc_x[k] * ksd_178[k];

        t_295[k] = f_3 * pc_x[k] * ksd_179[k];

        t_296[k] = pa_z[k] * isf0_216[k]
                   - f_4 * pc_z[k] * isf1_216[k];

        t_297[k] = f_5 * isd_129[k]
                   + f_3 * pc_z[k] * ksd_177[k];

        t_298[k] = f_6 * isd_137[k]
                   + f_3 * pc_y[k] * ksd_179[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pc_x, pc_z, isd_131, ksp0_89, ksp0_90, ksp0_91, \
                         ksp1_89, ksp1_90, ksp1_91, ksd_179, ksd_180, \
                         ksd_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_5 * isd_131[k]
                   + f_1 * ksp0_89[k]
                   - f_2 * ksp1_89[k]
                   + f_3 * pc_z[k] * ksd_179[k];

        t_300[k] = f_1 * ksp0_90[k]
                   - f_2 * ksp1_90[k]
                   + f_3 * pc_x[k] * ksd_180[k];

        t_301[k] = f_7 * ksp0_91[k]
                   - f_8 * ksp1_91[k]
                   + f_3 * pc_x[k] * ksd_181[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, pc_x, pc_y, isd_141, ksp0_91, \
                         ksp0_92, ksp1_91, ksp1_92, ksd_182, ksd_183, ksd_184, \
                         ksd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_7 * ksp0_92[k]
                   - f_8 * ksp1_92[k]
                   + f_3 * pc_x[k] * ksd_182[k];

        t_303[k] = f_3 * pc_x[k] * ksd_183[k];

        t_304[k] = f_3 * pc_x[k] * ksd_184[k];

        t_305[k] = f_3 * pc_x[k] * ksd_185[k];

        t_306[k] = f_9 * isd_141[k]
                   + f_1 * ksp0_91[k]
                   - f_2 * ksp1_91[k]
                   + f_3 * pc_y[k] * ksd_183[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, pc_y, pc_z, isd_135, isd_137, isd_143, ksp0_92, \
                         ksp1_92, ksd_183, ksd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_10 * isd_135[k]
                   + f_3 * pc_z[k] * ksd_183[k];

        t_308[k] = f_9 * isd_143[k]
                   + f_3 * pc_y[k] * ksd_185[k];

        t_309[k] = f_10 * isd_137[k]
                   + f_1 * ksp0_92[k]
                   - f_2 * ksp1_92[k]
                   + f_3 * pc_z[k] * ksd_185[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pc_x, ksp0_93, ksp0_94, ksp0_95, ksp1_93, \
                         ksp1_94, ksp1_95, ksd_186, ksd_187, ksd_188, \
                         ksd_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_1 * ksp0_93[k]
                   - f_2 * ksp1_93[k]
                   + f_3 * pc_x[k] * ksd_186[k];

        t_311[k] = f_7 * ksp0_94[k]
                   - f_8 * ksp1_94[k]
                   + f_3 * pc_x[k] * ksd_187[k];

        t_312[k] = f_7 * ksp0_95[k]
                   - f_8 * ksp1_95[k]
                   + f_3 * pc_x[k] * ksd_188[k];

        t_313[k] = f_3 * pc_x[k] * ksd_189[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, isd_141, \
                         isd_147, isd_149, ksp0_94, ksp1_94, ksd_189, ksd_190, \
                         ksd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_3 * pc_x[k] * ksd_190[k];

        t_315[k] = f_3 * pc_x[k] * ksd_191[k];

        t_316[k] = f_11 * isd_147[k]
                   + f_1 * ksp0_94[k]
                   - f_2 * ksp1_94[k]
                   + f_3 * pc_y[k] * ksd_189[k];

        t_317[k] = f_12 * isd_141[k]
                   + f_3 * pc_z[k] * ksd_189[k];

        t_318[k] = f_11 * isd_149[k]
                   + f_3 * pc_y[k] * ksd_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pc_x, pc_z, isd_143, ksp0_95, ksp0_96, ksp0_97, \
                         ksp1_95, ksp1_96, ksp1_97, ksd_191, ksd_192, \
                         ksd_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_12 * isd_143[k]
                   + f_1 * ksp0_95[k]
                   - f_2 * ksp1_95[k]
                   + f_3 * pc_z[k] * ksd_191[k];

        t_320[k] = f_1 * ksp0_96[k]
                   - f_2 * ksp1_96[k]
                   + f_3 * pc_x[k] * ksd_192[k];

        t_321[k] = f_7 * ksp0_97[k]
                   - f_8 * ksp1_97[k]
                   + f_3 * pc_x[k] * ksd_193[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, pc_x, pc_y, isd_153, ksp0_97, \
                         ksp0_98, ksp1_97, ksp1_98, ksd_194, ksd_195, ksd_196, \
                         ksd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_7 * ksp0_98[k]
                   - f_8 * ksp1_98[k]
                   + f_3 * pc_x[k] * ksd_194[k];

        t_323[k] = f_3 * pc_x[k] * ksd_195[k];

        t_324[k] = f_3 * pc_x[k] * ksd_196[k];

        t_325[k] = f_3 * pc_x[k] * ksd_197[k];

        t_326[k] = f_12 * isd_153[k]
                   + f_1 * ksp0_97[k]
                   - f_2 * ksp1_97[k]
                   + f_3 * pc_y[k] * ksd_195[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, pc_y, pc_z, isd_147, isd_149, isd_155, ksp0_98, \
                         ksp1_98, ksd_195, ksd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_11 * isd_147[k]
                   + f_3 * pc_z[k] * ksd_195[k];

        t_328[k] = f_12 * isd_155[k]
                   + f_3 * pc_y[k] * ksd_197[k];

        t_329[k] = f_11 * isd_149[k]
                   + f_1 * ksp0_98[k]
                   - f_2 * ksp1_98[k]
                   + f_3 * pc_z[k] * ksd_197[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, ksp0_99, ksp0_100, ksp0_101, \
                         ksp1_99, ksp1_100, ksp1_101, ksd_198, ksd_199, ksd_200, \
                         ksd_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_1 * ksp0_99[k]
                   - f_2 * ksp1_99[k]
                   + f_3 * pc_x[k] * ksd_198[k];

        t_331[k] = f_7 * ksp0_100[k]
                   - f_8 * ksp1_100[k]
                   + f_3 * pc_x[k] * ksd_199[k];

        t_332[k] = f_7 * ksp0_101[k]
                   - f_8 * ksp1_101[k]
                   + f_3 * pc_x[k] * ksd_200[k];

        t_333[k] = f_3 * pc_x[k] * ksd_201[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, pc_x, pc_y, pc_z, isd_153, \
                         isd_159, isd_161, ksp0_100, ksp1_100, ksd_201, ksd_202, \
                         ksd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_3 * pc_x[k] * ksd_202[k];

        t_335[k] = f_3 * pc_x[k] * ksd_203[k];

        t_336[k] = f_10 * isd_159[k]
                   + f_1 * ksp0_100[k]
                   - f_2 * ksp1_100[k]
                   + f_3 * pc_y[k] * ksd_201[k];

        t_337[k] = f_9 * isd_153[k]
                   + f_3 * pc_z[k] * ksd_201[k];

        t_338[k] = f_10 * isd_161[k]
                   + f_3 * pc_y[k] * ksd_203[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_y, pc_x, pc_y, pc_z, isf0_270, isd_155, \
                         isf1_270, ksp0_101, ksp0_103, ksp1_101, ksp1_103, ksd_203, \
                         ksd_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_9 * isd_155[k]
                   + f_1 * ksp0_101[k]
                   - f_2 * ksp1_101[k]
                   + f_3 * pc_z[k] * ksd_203[k];

        t_340[k] = pa_y[k] * isf0_270[k]
                   - f_4 * pc_y[k] * isf1_270[k];

        t_341[k] = f_7 * ksp0_103[k]
                   - f_8 * ksp1_103[k]
                   + f_3 * pc_x[k] * ksd_205[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, pa_y, pc_x, pc_y, isf0_272, \
                         isf0_276, isd_165, isf1_272, isf1_276, ksd_207, ksd_208, \
                         ksd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_y[k] * isf0_272[k]
                   - f_4 * pc_y[k] * isf1_272[k];

        t_343[k] = f_3 * pc_x[k] * ksd_207[k];

        t_344[k] = f_3 * pc_x[k] * ksd_208[k];

        t_345[k] = f_3 * pc_x[k] * ksd_209[k];

        t_346[k] = pa_y[k] * isf0_276[k]
                   + f_12 * isd_165[k]
                   - f_4 * pc_y[k] * isf1_276[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pa_y, pc_y, pc_z, isf0_279, isd_159, isd_167, \
                         isf1_279, ksd_207, ksd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_6 * isd_159[k]
                   + f_3 * pc_z[k] * ksd_207[k];

        t_348[k] = f_5 * isd_167[k]
                   + f_3 * pc_y[k] * ksd_209[k];

        t_349[k] = pa_y[k] * isf0_279[k]
                   - f_4 * pc_y[k] * isf1_279[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, pc_x, pc_y, ksp0_105, ksp0_107, \
                         ksp1_105, ksp1_107, ksd_210, ksd_212, ksd_213, \
                         ksd_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_1 * ksp0_105[k]
                   - f_2 * ksp1_105[k]
                   + f_3 * pc_x[k] * ksd_210[k];

        t_351[k] = f_3 * pc_y[k] * ksd_210[k];

        t_352[k] = f_7 * ksp0_107[k]
                   - f_8 * ksp1_107[k]
                   + f_3 * pc_x[k] * ksd_212[k];

        t_353[k] = f_3 * pc_x[k] * ksd_213[k];

        t_354[k] = f_3 * pc_x[k] * ksd_214[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, pc_x, pc_y, pc_z, isd_167, \
                         ksp0_106, ksp0_107, ksp1_106, ksp1_107, ksd_213, ksd_214, \
                         ksd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_3 * pc_x[k] * ksd_215[k];

        t_356[k] = f_1 * ksp0_106[k]
                   - f_2 * ksp1_106[k]
                   + f_3 * pc_y[k] * ksd_213[k];

        t_357[k] = f_7 * ksp0_107[k]
                   - f_8 * ksp1_107[k]
                   + f_3 * pc_y[k] * ksd_214[k];

        t_358[k] = f_3 * pc_y[k] * ksd_215[k];

        t_359[k] = f_0 * isd_167[k]
                   + f_1 * ksp0_107[k]
                   - f_2 * ksp1_107[k]
                   + f_3 * pc_z[k] * ksd_215[k];
    }
}

auto
compute_prim_ksf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t isf0, const size_t isd,
                                                   const size_t isf1, const size_t ksp0,
                                                   const size_t ksp1, const size_t ksd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, isf0, isd,
                                                              isf1, ksp0, ksp1, ksd, ncols,
                                                              gamma, p, q);

    compute_prim_ksf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, isf0, isd,
                                                              isf1, ksp0, ksp1, ksd, ncols,
                                                              gamma, p, q);

    compute_prim_ksf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, isf0, isd,
                                                              isf1, ksp0, ksp1, ksd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
