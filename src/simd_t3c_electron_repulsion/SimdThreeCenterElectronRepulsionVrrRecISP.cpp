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


#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_isp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t hss,
                                                   const size_t iss, const size_t ncols,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 2.5 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 2.0 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 1.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hss_0 = buffer.data(hss + 0);
    const auto *hss_1 = buffer.data(hss + 1);
    const auto *hss_2 = buffer.data(hss + 2);
    const auto *hss_3 = buffer.data(hss + 3);
    const auto *hss_4 = buffer.data(hss + 4);
    const auto *hss_5 = buffer.data(hss + 5);
    const auto *hss_6 = buffer.data(hss + 6);
    const auto *hss_7 = buffer.data(hss + 7);
    const auto *hss_8 = buffer.data(hss + 8);
    const auto *hss_9 = buffer.data(hss + 9);
    const auto *hss_10 = buffer.data(hss + 10);
    const auto *hss_11 = buffer.data(hss + 11);
    const auto *hss_12 = buffer.data(hss + 12);
    const auto *hss_13 = buffer.data(hss + 13);
    const auto *hss_14 = buffer.data(hss + 14);
    const auto *hss_15 = buffer.data(hss + 15);
    const auto *hss_16 = buffer.data(hss + 16);
    const auto *hss_17 = buffer.data(hss + 17);
    const auto *hss_18 = buffer.data(hss + 18);
    const auto *hss_19 = buffer.data(hss + 19);
    const auto *hss_20 = buffer.data(hss + 20);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, hss_0, hss_1, \
                         hss_2, iss_0, iss_1, iss_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hss_0[k]
                 + f_1 * pc_x[k] * iss_0[k];

        t_1[k] = f_1 * pc_y[k] * iss_0[k];

        t_2[k] = f_1 * pc_z[k] * iss_0[k];

        t_3[k] = f_2 * hss_1[k]
                 + f_1 * pc_x[k] * iss_1[k];

        t_4[k] = f_3 * hss_0[k]
                 + f_1 * pc_y[k] * iss_1[k];

        t_5[k] = f_1 * pc_z[k] * iss_1[k];

        t_6[k] = f_2 * hss_2[k]
                 + f_1 * pc_x[k] * iss_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, hss_0, hss_1, \
                         hss_3, hss_4, iss_2, iss_3, iss_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * iss_2[k];

        t_8[k] = f_3 * hss_0[k]
                 + f_1 * pc_z[k] * iss_2[k];

        t_9[k] = f_4 * hss_3[k]
                 + f_1 * pc_x[k] * iss_3[k];

        t_10[k] = f_5 * hss_1[k]
                  + f_1 * pc_y[k] * iss_3[k];

        t_11[k] = f_1 * pc_z[k] * iss_3[k];

        t_12[k] = f_4 * hss_4[k]
                  + f_1 * pc_x[k] * iss_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, hss_1, hss_2, \
                         hss_5, hss_6, iss_4, iss_5, iss_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * hss_2[k]
                  + f_1 * pc_y[k] * iss_4[k];

        t_14[k] = f_3 * hss_1[k]
                  + f_1 * pc_z[k] * iss_4[k];

        t_15[k] = f_4 * hss_5[k]
                  + f_1 * pc_x[k] * iss_5[k];

        t_16[k] = f_1 * pc_y[k] * iss_5[k];

        t_17[k] = f_5 * hss_2[k]
                  + f_1 * pc_z[k] * iss_5[k];

        t_18[k] = f_6 * hss_6[k]
                  + f_1 * pc_x[k] * iss_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, hss_3, hss_4, \
                         hss_7, hss_8, iss_6, iss_7, iss_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_6 * hss_3[k]
                  + f_1 * pc_y[k] * iss_6[k];

        t_20[k] = f_1 * pc_z[k] * iss_6[k];

        t_21[k] = f_6 * hss_7[k]
                  + f_1 * pc_x[k] * iss_7[k];

        t_22[k] = f_5 * hss_4[k]
                  + f_1 * pc_y[k] * iss_7[k];

        t_23[k] = f_3 * hss_3[k]
                  + f_1 * pc_z[k] * iss_7[k];

        t_24[k] = f_6 * hss_8[k]
                  + f_1 * pc_x[k] * iss_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, hss_4, hss_5, \
                         hss_9, hss_10, iss_8, iss_9, iss_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * hss_5[k]
                  + f_1 * pc_y[k] * iss_8[k];

        t_26[k] = f_5 * hss_4[k]
                  + f_1 * pc_z[k] * iss_8[k];

        t_27[k] = f_6 * hss_9[k]
                  + f_1 * pc_x[k] * iss_9[k];

        t_28[k] = f_1 * pc_y[k] * iss_9[k];

        t_29[k] = f_6 * hss_5[k]
                  + f_1 * pc_z[k] * iss_9[k];

        t_30[k] = f_5 * hss_10[k]
                  + f_1 * pc_x[k] * iss_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, hss_6, hss_7, \
                         hss_11, hss_12, iss_10, iss_11, iss_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_4 * hss_6[k]
                  + f_1 * pc_y[k] * iss_10[k];

        t_32[k] = f_1 * pc_z[k] * iss_10[k];

        t_33[k] = f_5 * hss_11[k]
                  + f_1 * pc_x[k] * iss_11[k];

        t_34[k] = f_6 * hss_7[k]
                  + f_1 * pc_y[k] * iss_11[k];

        t_35[k] = f_3 * hss_6[k]
                  + f_1 * pc_z[k] * iss_11[k];

        t_36[k] = f_5 * hss_12[k]
                  + f_1 * pc_x[k] * iss_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, hss_7, hss_8, hss_9, \
                         hss_13, iss_12, iss_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * hss_8[k]
                  + f_1 * pc_y[k] * iss_12[k];

        t_38[k] = f_5 * hss_7[k]
                  + f_1 * pc_z[k] * iss_12[k];

        t_39[k] = f_5 * hss_13[k]
                  + f_1 * pc_x[k] * iss_13[k];

        t_40[k] = f_3 * hss_9[k]
                  + f_1 * pc_y[k] * iss_13[k];

        t_41[k] = f_6 * hss_8[k]
                  + f_1 * pc_z[k] * iss_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, hss_9, hss_10, \
                         hss_14, hss_15, iss_14, iss_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * hss_14[k]
                  + f_1 * pc_x[k] * iss_14[k];

        t_43[k] = f_1 * pc_y[k] * iss_14[k];

        t_44[k] = f_4 * hss_9[k]
                  + f_1 * pc_z[k] * iss_14[k];

        t_45[k] = f_3 * hss_15[k]
                  + f_1 * pc_x[k] * iss_15[k];

        t_46[k] = f_2 * hss_10[k]
                  + f_1 * pc_y[k] * iss_15[k];

        t_47[k] = f_1 * pc_z[k] * iss_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, hss_10, hss_11, \
                         hss_12, hss_16, hss_17, iss_16, iss_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * hss_16[k]
                  + f_1 * pc_x[k] * iss_16[k];

        t_49[k] = f_4 * hss_11[k]
                  + f_1 * pc_y[k] * iss_16[k];

        t_50[k] = f_3 * hss_10[k]
                  + f_1 * pc_z[k] * iss_16[k];

        t_51[k] = f_3 * hss_17[k]
                  + f_1 * pc_x[k] * iss_17[k];

        t_52[k] = f_6 * hss_12[k]
                  + f_1 * pc_y[k] * iss_17[k];

        t_53[k] = f_5 * hss_11[k]
                  + f_1 * pc_z[k] * iss_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, hss_12, hss_13, \
                         hss_14, hss_18, hss_19, iss_18, iss_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * hss_18[k]
                  + f_1 * pc_x[k] * iss_18[k];

        t_55[k] = f_5 * hss_13[k]
                  + f_1 * pc_y[k] * iss_18[k];

        t_56[k] = f_6 * hss_12[k]
                  + f_1 * pc_z[k] * iss_18[k];

        t_57[k] = f_3 * hss_19[k]
                  + f_1 * pc_x[k] * iss_19[k];

        t_58[k] = f_3 * hss_14[k]
                  + f_1 * pc_y[k] * iss_19[k];

        t_59[k] = f_4 * hss_13[k]
                  + f_1 * pc_z[k] * iss_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, hss_14, \
                         hss_15, hss_20, iss_20, iss_21, iss_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * hss_20[k]
                  + f_1 * pc_x[k] * iss_20[k];

        t_61[k] = f_1 * pc_y[k] * iss_20[k];

        t_62[k] = f_2 * hss_14[k]
                  + f_1 * pc_z[k] * iss_20[k];

        t_63[k] = f_1 * pc_x[k] * iss_21[k];

        t_64[k] = f_0 * hss_15[k]
                  + f_1 * pc_y[k] * iss_21[k];

        t_65[k] = f_1 * pc_z[k] * iss_21[k];

        t_66[k] = f_1 * pc_x[k] * iss_22[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, t_72, pc_x, pc_y, pc_z, hss_15, hss_16, \
                         hss_17, iss_22, iss_23, iss_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * hss_16[k]
                  + f_1 * pc_y[k] * iss_22[k];

        t_68[k] = f_3 * hss_15[k]
                  + f_1 * pc_z[k] * iss_22[k];

        t_69[k] = f_1 * pc_x[k] * iss_23[k];

        t_70[k] = f_4 * hss_17[k]
                  + f_1 * pc_y[k] * iss_23[k];

        t_71[k] = f_5 * hss_16[k]
                  + f_1 * pc_z[k] * iss_23[k];

        t_72[k] = f_1 * pc_x[k] * iss_24[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, hss_17, hss_18, \
                         hss_19, iss_24, iss_25, iss_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * hss_18[k]
                  + f_1 * pc_y[k] * iss_24[k];

        t_74[k] = f_6 * hss_17[k]
                  + f_1 * pc_z[k] * iss_24[k];

        t_75[k] = f_1 * pc_x[k] * iss_25[k];

        t_76[k] = f_5 * hss_19[k]
                  + f_1 * pc_y[k] * iss_25[k];

        t_77[k] = f_4 * hss_18[k]
                  + f_1 * pc_z[k] * iss_25[k];

        t_78[k] = f_1 * pc_x[k] * iss_26[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, hss_19, hss_20, \
                         iss_26, iss_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * hss_20[k]
                  + f_1 * pc_y[k] * iss_26[k];

        t_80[k] = f_2 * hss_19[k]
                  + f_1 * pc_z[k] * iss_26[k];

        t_81[k] = f_1 * pc_x[k] * iss_27[k];

        t_82[k] = f_1 * pc_y[k] * iss_27[k];

        t_83[k] = f_0 * hss_20[k]
                  + f_1 * pc_z[k] * iss_27[k];
    }
}

}  // namespace simdt3ceri
