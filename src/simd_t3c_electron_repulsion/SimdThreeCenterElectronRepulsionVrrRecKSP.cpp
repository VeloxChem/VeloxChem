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


#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ksp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t iss,
                                                   const size_t kss, const size_t ncols,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = p / q;
    const auto f_2 = 3.0 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 2.5 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 2.0 / q;
    const auto f_7 = 1.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *iss_0 = buffer.data(iss + 0);
    const auto *iss_1 = buffer.data(iss + 1);
    const auto *iss_2 = buffer.data(iss + 2);
    const auto *iss_3 = buffer.data(iss + 3);
    const auto *iss_4 = buffer.data(iss + 4);
    const auto *iss_5 = buffer.data(iss + 5);
    const auto *iss_6 = buffer.data(iss + 6);
    const auto *iss_7 = buffer.data(iss + 7);
    const auto *iss_8 = buffer.data(iss + 8);
    const auto *iss_9 = buffer.data(iss + 9);
    const auto *iss_10 = buffer.data(iss + 10);
    const auto *iss_11 = buffer.data(iss + 11);
    const auto *iss_12 = buffer.data(iss + 12);
    const auto *iss_13 = buffer.data(iss + 13);
    const auto *iss_14 = buffer.data(iss + 14);
    const auto *iss_15 = buffer.data(iss + 15);
    const auto *iss_16 = buffer.data(iss + 16);
    const auto *iss_17 = buffer.data(iss + 17);
    const auto *iss_18 = buffer.data(iss + 18);
    const auto *iss_19 = buffer.data(iss + 19);
    const auto *iss_20 = buffer.data(iss + 20);
    const auto *iss_21 = buffer.data(iss + 21);
    const auto *iss_22 = buffer.data(iss + 22);
    const auto *iss_23 = buffer.data(iss + 23);
    const auto *iss_24 = buffer.data(iss + 24);
    const auto *iss_25 = buffer.data(iss + 25);
    const auto *iss_26 = buffer.data(iss + 26);
    const auto *iss_27 = buffer.data(iss + 27);

    const auto *kss_0 = buffer.data(kss + 0);
    const auto *kss_1 = buffer.data(kss + 1);
    const auto *kss_2 = buffer.data(kss + 2);
    const auto *kss_3 = buffer.data(kss + 3);
    const auto *kss_4 = buffer.data(kss + 4);
    const auto *kss_5 = buffer.data(kss + 5);
    const auto *kss_6 = buffer.data(kss + 6);
    const auto *kss_7 = buffer.data(kss + 7);
    const auto *kss_8 = buffer.data(kss + 8);
    const auto *kss_9 = buffer.data(kss + 9);
    const auto *kss_10 = buffer.data(kss + 10);
    const auto *kss_11 = buffer.data(kss + 11);
    const auto *kss_12 = buffer.data(kss + 12);
    const auto *kss_13 = buffer.data(kss + 13);
    const auto *kss_14 = buffer.data(kss + 14);
    const auto *kss_15 = buffer.data(kss + 15);
    const auto *kss_16 = buffer.data(kss + 16);
    const auto *kss_17 = buffer.data(kss + 17);
    const auto *kss_18 = buffer.data(kss + 18);
    const auto *kss_19 = buffer.data(kss + 19);
    const auto *kss_20 = buffer.data(kss + 20);
    const auto *kss_21 = buffer.data(kss + 21);
    const auto *kss_22 = buffer.data(kss + 22);
    const auto *kss_23 = buffer.data(kss + 23);
    const auto *kss_24 = buffer.data(kss + 24);
    const auto *kss_25 = buffer.data(kss + 25);
    const auto *kss_26 = buffer.data(kss + 26);
    const auto *kss_27 = buffer.data(kss + 27);
    const auto *kss_28 = buffer.data(kss + 28);
    const auto *kss_29 = buffer.data(kss + 29);
    const auto *kss_30 = buffer.data(kss + 30);
    const auto *kss_31 = buffer.data(kss + 31);
    const auto *kss_32 = buffer.data(kss + 32);
    const auto *kss_33 = buffer.data(kss + 33);
    const auto *kss_34 = buffer.data(kss + 34);
    const auto *kss_35 = buffer.data(kss + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, iss_0, iss_1, \
                         iss_2, kss_0, kss_1, kss_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * iss_0[k]
                 + f_1 * pc_x[k] * kss_0[k];

        t_1[k] = f_1 * pc_y[k] * kss_0[k];

        t_2[k] = f_1 * pc_z[k] * kss_0[k];

        t_3[k] = f_2 * iss_1[k]
                 + f_1 * pc_x[k] * kss_1[k];

        t_4[k] = f_3 * iss_0[k]
                 + f_1 * pc_y[k] * kss_1[k];

        t_5[k] = f_1 * pc_z[k] * kss_1[k];

        t_6[k] = f_2 * iss_2[k]
                 + f_1 * pc_x[k] * kss_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, iss_0, iss_1, \
                         iss_3, iss_4, kss_2, kss_3, kss_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * kss_2[k];

        t_8[k] = f_3 * iss_0[k]
                 + f_1 * pc_z[k] * kss_2[k];

        t_9[k] = f_4 * iss_3[k]
                 + f_1 * pc_x[k] * kss_3[k];

        t_10[k] = f_5 * iss_1[k]
                  + f_1 * pc_y[k] * kss_3[k];

        t_11[k] = f_1 * pc_z[k] * kss_3[k];

        t_12[k] = f_4 * iss_4[k]
                  + f_1 * pc_x[k] * kss_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, iss_1, iss_2, \
                         iss_5, iss_6, kss_4, kss_5, kss_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * iss_2[k]
                  + f_1 * pc_y[k] * kss_4[k];

        t_14[k] = f_3 * iss_1[k]
                  + f_1 * pc_z[k] * kss_4[k];

        t_15[k] = f_4 * iss_5[k]
                  + f_1 * pc_x[k] * kss_5[k];

        t_16[k] = f_1 * pc_y[k] * kss_5[k];

        t_17[k] = f_5 * iss_2[k]
                  + f_1 * pc_z[k] * kss_5[k];

        t_18[k] = f_6 * iss_6[k]
                  + f_1 * pc_x[k] * kss_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, iss_3, iss_4, \
                         iss_7, iss_8, kss_6, kss_7, kss_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * iss_3[k]
                  + f_1 * pc_y[k] * kss_6[k];

        t_20[k] = f_1 * pc_z[k] * kss_6[k];

        t_21[k] = f_6 * iss_7[k]
                  + f_1 * pc_x[k] * kss_7[k];

        t_22[k] = f_5 * iss_4[k]
                  + f_1 * pc_y[k] * kss_7[k];

        t_23[k] = f_3 * iss_3[k]
                  + f_1 * pc_z[k] * kss_7[k];

        t_24[k] = f_6 * iss_8[k]
                  + f_1 * pc_x[k] * kss_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, iss_4, iss_5, \
                         iss_9, iss_10, kss_8, kss_9, kss_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * iss_5[k]
                  + f_1 * pc_y[k] * kss_8[k];

        t_26[k] = f_5 * iss_4[k]
                  + f_1 * pc_z[k] * kss_8[k];

        t_27[k] = f_6 * iss_9[k]
                  + f_1 * pc_x[k] * kss_9[k];

        t_28[k] = f_1 * pc_y[k] * kss_9[k];

        t_29[k] = f_7 * iss_5[k]
                  + f_1 * pc_z[k] * kss_9[k];

        t_30[k] = f_7 * iss_10[k]
                  + f_1 * pc_x[k] * kss_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, iss_6, iss_7, \
                         iss_11, iss_12, kss_10, kss_11, kss_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * iss_6[k]
                  + f_1 * pc_y[k] * kss_10[k];

        t_32[k] = f_1 * pc_z[k] * kss_10[k];

        t_33[k] = f_7 * iss_11[k]
                  + f_1 * pc_x[k] * kss_11[k];

        t_34[k] = f_7 * iss_7[k]
                  + f_1 * pc_y[k] * kss_11[k];

        t_35[k] = f_3 * iss_6[k]
                  + f_1 * pc_z[k] * kss_11[k];

        t_36[k] = f_7 * iss_12[k]
                  + f_1 * pc_x[k] * kss_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, iss_7, iss_8, iss_9, \
                         iss_13, kss_12, kss_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * iss_8[k]
                  + f_1 * pc_y[k] * kss_12[k];

        t_38[k] = f_5 * iss_7[k]
                  + f_1 * pc_z[k] * kss_12[k];

        t_39[k] = f_7 * iss_13[k]
                  + f_1 * pc_x[k] * kss_13[k];

        t_40[k] = f_3 * iss_9[k]
                  + f_1 * pc_y[k] * kss_13[k];

        t_41[k] = f_7 * iss_8[k]
                  + f_1 * pc_z[k] * kss_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, iss_9, iss_10, \
                         iss_14, iss_15, kss_14, kss_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_7 * iss_14[k]
                  + f_1 * pc_x[k] * kss_14[k];

        t_43[k] = f_1 * pc_y[k] * kss_14[k];

        t_44[k] = f_6 * iss_9[k]
                  + f_1 * pc_z[k] * kss_14[k];

        t_45[k] = f_5 * iss_15[k]
                  + f_1 * pc_x[k] * kss_15[k];

        t_46[k] = f_4 * iss_10[k]
                  + f_1 * pc_y[k] * kss_15[k];

        t_47[k] = f_1 * pc_z[k] * kss_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, iss_10, iss_11, \
                         iss_12, iss_16, iss_17, kss_16, kss_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * iss_16[k]
                  + f_1 * pc_x[k] * kss_16[k];

        t_49[k] = f_6 * iss_11[k]
                  + f_1 * pc_y[k] * kss_16[k];

        t_50[k] = f_3 * iss_10[k]
                  + f_1 * pc_z[k] * kss_16[k];

        t_51[k] = f_5 * iss_17[k]
                  + f_1 * pc_x[k] * kss_17[k];

        t_52[k] = f_7 * iss_12[k]
                  + f_1 * pc_y[k] * kss_17[k];

        t_53[k] = f_5 * iss_11[k]
                  + f_1 * pc_z[k] * kss_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, iss_12, iss_13, \
                         iss_14, iss_18, iss_19, kss_18, kss_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * iss_18[k]
                  + f_1 * pc_x[k] * kss_18[k];

        t_55[k] = f_5 * iss_13[k]
                  + f_1 * pc_y[k] * kss_18[k];

        t_56[k] = f_7 * iss_12[k]
                  + f_1 * pc_z[k] * kss_18[k];

        t_57[k] = f_5 * iss_19[k]
                  + f_1 * pc_x[k] * kss_19[k];

        t_58[k] = f_3 * iss_14[k]
                  + f_1 * pc_y[k] * kss_19[k];

        t_59[k] = f_6 * iss_13[k]
                  + f_1 * pc_z[k] * kss_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, iss_14, iss_15, \
                         iss_20, iss_21, kss_20, kss_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * iss_20[k]
                  + f_1 * pc_x[k] * kss_20[k];

        t_61[k] = f_1 * pc_y[k] * kss_20[k];

        t_62[k] = f_4 * iss_14[k]
                  + f_1 * pc_z[k] * kss_20[k];

        t_63[k] = f_3 * iss_21[k]
                  + f_1 * pc_x[k] * kss_21[k];

        t_64[k] = f_2 * iss_15[k]
                  + f_1 * pc_y[k] * kss_21[k];

        t_65[k] = f_1 * pc_z[k] * kss_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pc_x, pc_y, pc_z, iss_15, iss_16, \
                         iss_17, iss_22, iss_23, kss_22, kss_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * iss_22[k]
                  + f_1 * pc_x[k] * kss_22[k];

        t_67[k] = f_4 * iss_16[k]
                  + f_1 * pc_y[k] * kss_22[k];

        t_68[k] = f_3 * iss_15[k]
                  + f_1 * pc_z[k] * kss_22[k];

        t_69[k] = f_3 * iss_23[k]
                  + f_1 * pc_x[k] * kss_23[k];

        t_70[k] = f_6 * iss_17[k]
                  + f_1 * pc_y[k] * kss_23[k];

        t_71[k] = f_5 * iss_16[k]
                  + f_1 * pc_z[k] * kss_23[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, iss_17, iss_18, \
                         iss_19, iss_24, iss_25, kss_24, kss_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * iss_24[k]
                  + f_1 * pc_x[k] * kss_24[k];

        t_73[k] = f_7 * iss_18[k]
                  + f_1 * pc_y[k] * kss_24[k];

        t_74[k] = f_7 * iss_17[k]
                  + f_1 * pc_z[k] * kss_24[k];

        t_75[k] = f_3 * iss_25[k]
                  + f_1 * pc_x[k] * kss_25[k];

        t_76[k] = f_5 * iss_19[k]
                  + f_1 * pc_y[k] * kss_25[k];

        t_77[k] = f_6 * iss_18[k]
                  + f_1 * pc_z[k] * kss_25[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, iss_19, iss_20, \
                         iss_26, iss_27, kss_26, kss_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_3 * iss_26[k]
                  + f_1 * pc_x[k] * kss_26[k];

        t_79[k] = f_3 * iss_20[k]
                  + f_1 * pc_y[k] * kss_26[k];

        t_80[k] = f_4 * iss_19[k]
                  + f_1 * pc_z[k] * kss_26[k];

        t_81[k] = f_3 * iss_27[k]
                  + f_1 * pc_x[k] * kss_27[k];

        t_82[k] = f_1 * pc_y[k] * kss_27[k];

        t_83[k] = f_2 * iss_20[k]
                  + f_1 * pc_z[k] * kss_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, iss_21, \
                         iss_22, kss_28, kss_29, kss_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * pc_x[k] * kss_28[k];

        t_85[k] = f_0 * iss_21[k]
                  + f_1 * pc_y[k] * kss_28[k];

        t_86[k] = f_1 * pc_z[k] * kss_28[k];

        t_87[k] = f_1 * pc_x[k] * kss_29[k];

        t_88[k] = f_2 * iss_22[k]
                  + f_1 * pc_y[k] * kss_29[k];

        t_89[k] = f_3 * iss_21[k]
                  + f_1 * pc_z[k] * kss_29[k];

        t_90[k] = f_1 * pc_x[k] * kss_30[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, t_96, pc_x, pc_y, pc_z, iss_22, iss_23, \
                         iss_24, kss_30, kss_31, kss_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_4 * iss_23[k]
                  + f_1 * pc_y[k] * kss_30[k];

        t_92[k] = f_5 * iss_22[k]
                  + f_1 * pc_z[k] * kss_30[k];

        t_93[k] = f_1 * pc_x[k] * kss_31[k];

        t_94[k] = f_6 * iss_24[k]
                  + f_1 * pc_y[k] * kss_31[k];

        t_95[k] = f_7 * iss_23[k]
                  + f_1 * pc_z[k] * kss_31[k];

        t_96[k] = f_1 * pc_x[k] * kss_32[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, iss_24, \
                         iss_25, iss_26, kss_32, kss_33, kss_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * iss_25[k]
                  + f_1 * pc_y[k] * kss_32[k];

        t_98[k] = f_6 * iss_24[k]
                  + f_1 * pc_z[k] * kss_32[k];

        t_99[k] = f_1 * pc_x[k] * kss_33[k];

        t_100[k] = f_5 * iss_26[k]
                   + f_1 * pc_y[k] * kss_33[k];

        t_101[k] = f_4 * iss_25[k]
                   + f_1 * pc_z[k] * kss_33[k];

        t_102[k] = f_1 * pc_x[k] * kss_34[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, iss_26, iss_27, \
                         kss_34, kss_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_3 * iss_27[k]
                   + f_1 * pc_y[k] * kss_34[k];

        t_104[k] = f_2 * iss_26[k]
                   + f_1 * pc_z[k] * kss_34[k];

        t_105[k] = f_1 * pc_x[k] * kss_35[k];

        t_106[k] = f_1 * pc_y[k] * kss_35[k];

        t_107[k] = f_0 * iss_27[k]
                   + f_1 * pc_z[k] * kss_35[k];
    }
}

}  // namespace simdt3ceri
