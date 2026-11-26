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


#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_fsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsf0, const size_t dsd,
                                                   const size_t dsf1, const size_t fsp0,
                                                   const size_t fsp1, const size_t fsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.0 / q;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsf0_0 = buffer.data(dsf0 + 0);
    const auto *dsf0_6 = buffer.data(dsf0 + 6);
    const auto *dsf0_9 = buffer.data(dsf0 + 9);
    const auto *dsf0_20 = buffer.data(dsf0 + 20);
    const auto *dsf0_30 = buffer.data(dsf0 + 30);
    const auto *dsf0_31 = buffer.data(dsf0 + 31);
    const auto *dsf0_36 = buffer.data(dsf0 + 36);
    const auto *dsf0_39 = buffer.data(dsf0 + 39);
    const auto *dsf0_46 = buffer.data(dsf0 + 46);
    const auto *dsf0_49 = buffer.data(dsf0 + 49);
    const auto *dsf0_50 = buffer.data(dsf0 + 50);
    const auto *dsf0_52 = buffer.data(dsf0 + 52);
    const auto *dsf0_56 = buffer.data(dsf0 + 56);
    const auto *dsf0_57 = buffer.data(dsf0 + 57);
    const auto *dsf0_59 = buffer.data(dsf0 + 59);

    const auto *dsd_0 = buffer.data(dsd + 0);
    const auto *dsd_3 = buffer.data(dsd + 3);
    const auto *dsd_5 = buffer.data(dsd + 5);
    const auto *dsd_6 = buffer.data(dsd + 6);
    const auto *dsd_9 = buffer.data(dsd + 9);
    const auto *dsd_11 = buffer.data(dsd + 11);
    const auto *dsd_12 = buffer.data(dsd + 12);
    const auto *dsd_15 = buffer.data(dsd + 15);
    const auto *dsd_17 = buffer.data(dsd + 17);
    const auto *dsd_18 = buffer.data(dsd + 18);
    const auto *dsd_21 = buffer.data(dsd + 21);
    const auto *dsd_23 = buffer.data(dsd + 23);
    const auto *dsd_27 = buffer.data(dsd + 27);
    const auto *dsd_28 = buffer.data(dsd + 28);
    const auto *dsd_29 = buffer.data(dsd + 29);
    const auto *dsd_30 = buffer.data(dsd + 30);
    const auto *dsd_33 = buffer.data(dsd + 33);
    const auto *dsd_35 = buffer.data(dsd + 35);

    const auto *dsf1_0 = buffer.data(dsf1 + 0);
    const auto *dsf1_6 = buffer.data(dsf1 + 6);
    const auto *dsf1_9 = buffer.data(dsf1 + 9);
    const auto *dsf1_20 = buffer.data(dsf1 + 20);
    const auto *dsf1_30 = buffer.data(dsf1 + 30);
    const auto *dsf1_31 = buffer.data(dsf1 + 31);
    const auto *dsf1_36 = buffer.data(dsf1 + 36);
    const auto *dsf1_39 = buffer.data(dsf1 + 39);
    const auto *dsf1_46 = buffer.data(dsf1 + 46);
    const auto *dsf1_49 = buffer.data(dsf1 + 49);
    const auto *dsf1_50 = buffer.data(dsf1 + 50);
    const auto *dsf1_52 = buffer.data(dsf1 + 52);
    const auto *dsf1_56 = buffer.data(dsf1 + 56);
    const auto *dsf1_57 = buffer.data(dsf1 + 57);
    const auto *dsf1_59 = buffer.data(dsf1 + 59);

    const auto *fsp0_0 = buffer.data(fsp0 + 0);
    const auto *fsp0_1 = buffer.data(fsp0 + 1);
    const auto *fsp0_2 = buffer.data(fsp0 + 2);
    const auto *fsp0_4 = buffer.data(fsp0 + 4);
    const auto *fsp0_8 = buffer.data(fsp0 + 8);
    const auto *fsp0_18 = buffer.data(fsp0 + 18);
    const auto *fsp0_19 = buffer.data(fsp0 + 19);
    const auto *fsp0_20 = buffer.data(fsp0 + 20);
    const auto *fsp0_23 = buffer.data(fsp0 + 23);
    const auto *fsp0_25 = buffer.data(fsp0 + 25);
    const auto *fsp0_27 = buffer.data(fsp0 + 27);
    const auto *fsp0_28 = buffer.data(fsp0 + 28);
    const auto *fsp0_29 = buffer.data(fsp0 + 29);

    const auto *fsp1_0 = buffer.data(fsp1 + 0);
    const auto *fsp1_1 = buffer.data(fsp1 + 1);
    const auto *fsp1_2 = buffer.data(fsp1 + 2);
    const auto *fsp1_4 = buffer.data(fsp1 + 4);
    const auto *fsp1_8 = buffer.data(fsp1 + 8);
    const auto *fsp1_18 = buffer.data(fsp1 + 18);
    const auto *fsp1_19 = buffer.data(fsp1 + 19);
    const auto *fsp1_20 = buffer.data(fsp1 + 20);
    const auto *fsp1_23 = buffer.data(fsp1 + 23);
    const auto *fsp1_25 = buffer.data(fsp1 + 25);
    const auto *fsp1_27 = buffer.data(fsp1 + 27);
    const auto *fsp1_28 = buffer.data(fsp1 + 28);
    const auto *fsp1_29 = buffer.data(fsp1 + 29);

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
    const auto *fsd_15 = buffer.data(fsd + 15);
    const auto *fsd_16 = buffer.data(fsd + 16);
    const auto *fsd_17 = buffer.data(fsd + 17);
    const auto *fsd_18 = buffer.data(fsd + 18);
    const auto *fsd_19 = buffer.data(fsd + 19);
    const auto *fsd_21 = buffer.data(fsd + 21);
    const auto *fsd_23 = buffer.data(fsd + 23);
    const auto *fsd_24 = buffer.data(fsd + 24);
    const auto *fsd_27 = buffer.data(fsd + 27);
    const auto *fsd_28 = buffer.data(fsd + 28);
    const auto *fsd_29 = buffer.data(fsd + 29);
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
    const auto *fsd_49 = buffer.data(fsd + 49);
    const auto *fsd_51 = buffer.data(fsd + 51);
    const auto *fsd_52 = buffer.data(fsd + 52);
    const auto *fsd_53 = buffer.data(fsd + 53);
    const auto *fsd_54 = buffer.data(fsd + 54);
    const auto *fsd_56 = buffer.data(fsd + 56);
    const auto *fsd_57 = buffer.data(fsd + 57);
    const auto *fsd_58 = buffer.data(fsd + 58);
    const auto *fsd_59 = buffer.data(fsd + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, dsd_0, dsd_3, fsp0_0, \
                         fsp1_0, fsd_0, fsd_2, fsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsd_0[k]
                 + f_1 * fsp0_0[k]
                 - f_2 * fsp1_0[k]
                 + f_3 * pc_x[k] * fsd_0[k];

        t_1[k] = f_3 * pc_y[k] * fsd_0[k];

        t_2[k] = f_3 * pc_z[k] * fsd_0[k];

        t_3[k] = f_0 * dsd_3[k]
                 + f_3 * pc_x[k] * fsd_3[k];

        t_4[k] = f_3 * pc_y[k] * fsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, dsd_5, fsp0_1, fsp0_2, \
                         fsp1_1, fsp1_2, fsd_3, fsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dsd_5[k]
                 + f_3 * pc_x[k] * fsd_5[k];

        t_6[k] = f_1 * fsp0_1[k]
                 - f_2 * fsp1_1[k]
                 + f_3 * pc_y[k] * fsd_3[k];

        t_7[k] = f_3 * pc_z[k] * fsd_3[k];

        t_8[k] = f_3 * pc_y[k] * fsd_5[k];

        t_9[k] = f_1 * fsp0_2[k]
                 - f_2 * fsp1_2[k]
                 + f_3 * pc_z[k] * fsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, dsf0_0, dsd_0, \
                         dsd_9, dsf1_0, fsd_6, fsd_7, fsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * dsf0_0[k]
                  - f_4 * pc_y[k] * dsf1_0[k];

        t_11[k] = f_5 * dsd_0[k]
                  + f_3 * pc_y[k] * fsd_6[k];

        t_12[k] = f_3 * pc_z[k] * fsd_6[k];

        t_13[k] = f_6 * dsd_9[k]
                  + f_3 * pc_x[k] * fsd_9[k];

        t_14[k] = f_3 * pc_z[k] * fsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, dsd_3, dsd_5, dsd_11, \
                         fsp0_4, fsp1_4, fsd_9, fsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * dsd_11[k]
                  + f_3 * pc_x[k] * fsd_11[k];

        t_16[k] = f_5 * dsd_3[k]
                  + f_1 * fsp0_4[k]
                  - f_2 * fsp1_4[k]
                  + f_3 * pc_y[k] * fsd_9[k];

        t_17[k] = f_3 * pc_z[k] * fsd_9[k];

        t_18[k] = f_5 * dsd_5[k]
                  + f_3 * pc_y[k] * fsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, dsf0_0, dsf0_9, \
                         dsd_0, dsf1_0, dsf1_9, fsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * dsf0_9[k]
                  - f_4 * pc_y[k] * dsf1_9[k];

        t_20[k] = pa_z[k] * dsf0_0[k]
                  - f_4 * pc_z[k] * dsf1_0[k];

        t_21[k] = f_3 * pc_y[k] * fsd_12[k];

        t_22[k] = f_5 * dsd_0[k]
                  + f_3 * pc_z[k] * fsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, dsf0_6, dsd_15, \
                         dsd_17, dsf1_6, fsd_14, fsd_15, fsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * dsd_15[k]
                  + f_3 * pc_x[k] * fsd_15[k];

        t_24[k] = f_3 * pc_y[k] * fsd_14[k];

        t_25[k] = f_6 * dsd_17[k]
                  + f_3 * pc_x[k] * fsd_17[k];

        t_26[k] = pa_z[k] * dsf0_6[k]
                  - f_4 * pc_z[k] * dsf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pc_x, pc_y, pc_z, dsf0_30, dsd_5, \
                         dsd_18, dsf1_30, fsp0_8, fsp1_8, fsd_16, \
                         fsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * fsp0_8[k]
                  - f_8 * fsp1_8[k]
                  + f_3 * pc_y[k] * fsd_16[k];

        t_28[k] = f_3 * pc_y[k] * fsd_17[k];

        t_29[k] = f_5 * dsd_5[k]
                  + f_1 * fsp0_8[k]
                  - f_2 * fsp1_8[k]
                  + f_3 * pc_z[k] * fsd_17[k];

        t_30[k] = pa_x[k] * dsf0_30[k]
                  + f_0 * dsd_18[k]
                  - f_4 * pc_x[k] * dsf1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, dsd_6, dsd_21, \
                         dsd_23, fsd_18, fsd_19, fsd_21, fsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * dsd_6[k]
                  + f_3 * pc_y[k] * fsd_18[k];

        t_32[k] = f_3 * pc_z[k] * fsd_18[k];

        t_33[k] = f_5 * dsd_21[k]
                  + f_3 * pc_x[k] * fsd_21[k];

        t_34[k] = f_3 * pc_z[k] * fsd_19[k];

        t_35[k] = f_5 * dsd_23[k]
                  + f_3 * pc_x[k] * fsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pc_x, pc_y, pc_z, dsf0_36, dsf0_39, \
                         dsd_11, dsf1_36, dsf1_39, fsd_21, fsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * dsf0_36[k]
                  - f_4 * pc_x[k] * dsf1_36[k];

        t_37[k] = f_3 * pc_z[k] * fsd_21[k];

        t_38[k] = f_6 * dsd_11[k]
                  + f_3 * pc_y[k] * fsd_23[k];

        t_39[k] = pa_x[k] * dsf0_39[k]
                  - f_4 * pc_x[k] * dsf1_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, dsf0_20, dsd_6, \
                         dsd_12, dsd_27, dsf1_20, fsd_24, fsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * dsf0_20[k]
                  - f_4 * pc_y[k] * dsf1_20[k];

        t_41[k] = f_5 * dsd_12[k]
                  + f_3 * pc_y[k] * fsd_24[k];

        t_42[k] = f_5 * dsd_6[k]
                  + f_3 * pc_z[k] * fsd_24[k];

        t_43[k] = f_5 * dsd_27[k]
                  + f_3 * pc_x[k] * fsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pc_x, pc_z, dsf0_46, dsd_9, dsd_28, \
                         dsd_29, dsf1_46, fsd_27, fsd_28, fsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * dsd_28[k]
                  + f_3 * pc_x[k] * fsd_28[k];

        t_45[k] = f_5 * dsd_29[k]
                  + f_3 * pc_x[k] * fsd_29[k];

        t_46[k] = pa_x[k] * dsf0_46[k]
                  - f_4 * pc_x[k] * dsf1_46[k];

        t_47[k] = f_5 * dsd_9[k]
                  + f_3 * pc_z[k] * fsd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pc_x, pc_y, dsf0_49, dsf0_50, dsd_17, \
                         dsd_30, dsf1_49, dsf1_50, fsd_29, fsd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * dsd_17[k]
                  + f_3 * pc_y[k] * fsd_29[k];

        t_49[k] = pa_x[k] * dsf0_49[k]
                  - f_4 * pc_x[k] * dsf1_49[k];

        t_50[k] = pa_x[k] * dsf0_50[k]
                  + f_0 * dsd_30[k]
                  - f_4 * pc_x[k] * dsf1_50[k];

        t_51[k] = f_3 * pc_y[k] * fsd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, dsd_12, dsd_33, dsd_35, \
                         fsd_30, fsd_32, fsd_33, fsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_6 * dsd_12[k]
                  + f_3 * pc_z[k] * fsd_30[k];

        t_53[k] = f_5 * dsd_33[k]
                  + f_3 * pc_x[k] * fsd_33[k];

        t_54[k] = f_3 * pc_y[k] * fsd_32[k];

        t_55[k] = f_5 * dsd_35[k]
                  + f_3 * pc_x[k] * fsd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pc_x, pc_y, dsf0_56, dsf0_57, dsf0_59, \
                         dsf1_56, dsf1_57, dsf1_59, fsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * dsf0_56[k]
                  - f_4 * pc_x[k] * dsf1_56[k];

        t_57[k] = pa_x[k] * dsf0_57[k]
                  - f_4 * pc_x[k] * dsf1_57[k];

        t_58[k] = f_3 * pc_y[k] * fsd_35[k];

        t_59[k] = pa_x[k] * dsf0_59[k]
                  - f_4 * pc_x[k] * dsf1_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_z, fsp0_18, fsp0_19, fsp1_18, \
                         fsp1_19, fsd_36, fsd_37, fsd_39, fsd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * fsp0_18[k]
                  - f_2 * fsp1_18[k]
                  + f_3 * pc_x[k] * fsd_36[k];

        t_61[k] = f_7 * fsp0_19[k]
                  - f_8 * fsp1_19[k]
                  + f_3 * pc_x[k] * fsd_37[k];

        t_62[k] = f_3 * pc_z[k] * fsd_36[k];

        t_63[k] = f_3 * pc_x[k] * fsd_39[k];

        t_64[k] = f_3 * pc_x[k] * fsd_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pc_x, pc_y, pc_z, dsd_21, dsd_23, \
                         fsp0_19, fsp0_20, fsp1_19, fsp1_20, fsd_39, \
                         fsd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * pc_x[k] * fsd_41[k];

        t_66[k] = f_0 * dsd_21[k]
                  + f_1 * fsp0_19[k]
                  - f_2 * fsp1_19[k]
                  + f_3 * pc_y[k] * fsd_39[k];

        t_67[k] = f_3 * pc_z[k] * fsd_39[k];

        t_68[k] = f_0 * dsd_23[k]
                  + f_3 * pc_y[k] * fsd_41[k];

        t_69[k] = f_1 * fsp0_20[k]
                  - f_2 * fsp1_20[k]
                  + f_3 * pc_z[k] * fsd_41[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, dsf0_30, dsf0_31, dsf1_30, \
                         dsf1_31, fsp0_23, fsp1_23, fsd_44, fsd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * dsf0_30[k]
                  - f_4 * pc_z[k] * dsf1_30[k];

        t_71[k] = pa_z[k] * dsf0_31[k]
                  - f_4 * pc_z[k] * dsf1_31[k];

        t_72[k] = f_7 * fsp0_23[k]
                  - f_8 * fsp1_23[k]
                  + f_3 * pc_x[k] * fsd_44[k];

        t_73[k] = f_3 * pc_x[k] * fsd_45[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pa_z, pc_x, pc_y, pc_z, dsf0_36, \
                         dsd_21, dsd_29, dsf1_36, fsd_45, fsd_46, \
                         fsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * pc_x[k] * fsd_46[k];

        t_75[k] = f_3 * pc_x[k] * fsd_47[k];

        t_76[k] = pa_z[k] * dsf0_36[k]
                  - f_4 * pc_z[k] * dsf1_36[k];

        t_77[k] = f_5 * dsd_21[k]
                  + f_3 * pc_z[k] * fsd_45[k];

        t_78[k] = f_6 * dsd_29[k]
                  + f_3 * pc_y[k] * fsd_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pa_y, pc_x, pc_y, pc_z, dsf0_50, dsd_23, dsf1_50, \
                         fsp0_23, fsp0_25, fsp1_23, fsp1_25, fsd_47, \
                         fsd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * dsd_23[k]
                  + f_1 * fsp0_23[k]
                  - f_2 * fsp1_23[k]
                  + f_3 * pc_z[k] * fsd_47[k];

        t_80[k] = pa_y[k] * dsf0_50[k]
                  - f_4 * pc_y[k] * dsf1_50[k];

        t_81[k] = f_7 * fsp0_25[k]
                  - f_8 * fsp1_25[k]
                  + f_3 * pc_x[k] * fsd_49[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pc_x, pc_y, dsf0_52, dsf0_56, \
                         dsd_33, dsf1_52, dsf1_56, fsd_51, fsd_52, \
                         fsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_y[k] * dsf0_52[k]
                  - f_4 * pc_y[k] * dsf1_52[k];

        t_83[k] = f_3 * pc_x[k] * fsd_51[k];

        t_84[k] = f_3 * pc_x[k] * fsd_52[k];

        t_85[k] = f_3 * pc_x[k] * fsd_53[k];

        t_86[k] = pa_y[k] * dsf0_56[k]
                  + f_0 * dsd_33[k]
                  - f_4 * pc_y[k] * dsf1_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pc_y, pc_z, dsf0_59, dsd_27, dsd_35, dsf1_59, \
                         fsd_51, fsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_6 * dsd_27[k]
                  + f_3 * pc_z[k] * fsd_51[k];

        t_88[k] = f_5 * dsd_35[k]
                  + f_3 * pc_y[k] * fsd_53[k];

        t_89[k] = pa_y[k] * dsf0_59[k]
                  - f_4 * pc_y[k] * dsf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, fsp0_27, fsp0_29, fsp1_27, \
                         fsp1_29, fsd_54, fsd_56, fsd_57, fsd_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * fsp0_27[k]
                  - f_2 * fsp1_27[k]
                  + f_3 * pc_x[k] * fsd_54[k];

        t_91[k] = f_3 * pc_y[k] * fsd_54[k];

        t_92[k] = f_7 * fsp0_29[k]
                  - f_8 * fsp1_29[k]
                  + f_3 * pc_x[k] * fsd_56[k];

        t_93[k] = f_3 * pc_x[k] * fsd_57[k];

        t_94[k] = f_3 * pc_x[k] * fsd_58[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, dsd_35, fsp0_28, \
                         fsp0_29, fsp1_28, fsp1_29, fsd_57, fsd_58, \
                         fsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * pc_x[k] * fsd_59[k];

        t_96[k] = f_1 * fsp0_28[k]
                  - f_2 * fsp1_28[k]
                  + f_3 * pc_y[k] * fsd_57[k];

        t_97[k] = f_7 * fsp0_29[k]
                  - f_8 * fsp1_29[k]
                  + f_3 * pc_y[k] * fsd_58[k];

        t_98[k] = f_3 * pc_y[k] * fsd_59[k];

        t_99[k] = f_0 * dsd_35[k]
                  + f_1 * fsp0_29[k]
                  - f_2 * fsp1_29[k]
                  + f_3 * pc_z[k] * fsd_59[k];
    }
}

}  // namespace simdt3ceri
