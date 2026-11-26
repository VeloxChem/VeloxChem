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


#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pfp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfs, const size_t pdp0,
                                                   const size_t pds, const size_t pdp1,
                                                   const size_t pfs, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = p / q;
    const auto f_3 = gamma / q;
    const auto f_4 = 1.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfs_0 = buffer.data(sfs + 0);
    const auto *sfs_1 = buffer.data(sfs + 1);
    const auto *sfs_2 = buffer.data(sfs + 2);
    const auto *sfs_3 = buffer.data(sfs + 3);
    const auto *sfs_5 = buffer.data(sfs + 5);
    const auto *sfs_6 = buffer.data(sfs + 6);
    const auto *sfs_7 = buffer.data(sfs + 7);
    const auto *sfs_8 = buffer.data(sfs + 8);
    const auto *sfs_9 = buffer.data(sfs + 9);

    const auto *pdp0_0 = buffer.data(pdp0 + 0);
    const auto *pdp0_6 = buffer.data(pdp0 + 6);
    const auto *pdp0_28 = buffer.data(pdp0 + 28);
    const auto *pdp0_31 = buffer.data(pdp0 + 31);
    const auto *pdp0_50 = buffer.data(pdp0 + 50);
    const auto *pdp0_53 = buffer.data(pdp0 + 53);

    const auto *pds_0 = buffer.data(pds + 0);
    const auto *pds_1 = buffer.data(pds + 1);
    const auto *pds_2 = buffer.data(pds + 2);
    const auto *pds_3 = buffer.data(pds + 3);
    const auto *pds_4 = buffer.data(pds + 4);
    const auto *pds_5 = buffer.data(pds + 5);
    const auto *pds_6 = buffer.data(pds + 6);
    const auto *pds_7 = buffer.data(pds + 7);
    const auto *pds_8 = buffer.data(pds + 8);
    const auto *pds_9 = buffer.data(pds + 9);
    const auto *pds_10 = buffer.data(pds + 10);
    const auto *pds_11 = buffer.data(pds + 11);
    const auto *pds_12 = buffer.data(pds + 12);
    const auto *pds_13 = buffer.data(pds + 13);
    const auto *pds_14 = buffer.data(pds + 14);
    const auto *pds_15 = buffer.data(pds + 15);
    const auto *pds_16 = buffer.data(pds + 16);
    const auto *pds_17 = buffer.data(pds + 17);

    const auto *pdp1_0 = buffer.data(pdp1 + 0);
    const auto *pdp1_6 = buffer.data(pdp1 + 6);
    const auto *pdp1_28 = buffer.data(pdp1 + 28);
    const auto *pdp1_31 = buffer.data(pdp1 + 31);
    const auto *pdp1_50 = buffer.data(pdp1 + 50);
    const auto *pdp1_53 = buffer.data(pdp1 + 53);

    const auto *pfs_0 = buffer.data(pfs + 0);
    const auto *pfs_1 = buffer.data(pfs + 1);
    const auto *pfs_2 = buffer.data(pfs + 2);
    const auto *pfs_3 = buffer.data(pfs + 3);
    const auto *pfs_4 = buffer.data(pfs + 4);
    const auto *pfs_5 = buffer.data(pfs + 5);
    const auto *pfs_6 = buffer.data(pfs + 6);
    const auto *pfs_7 = buffer.data(pfs + 7);
    const auto *pfs_8 = buffer.data(pfs + 8);
    const auto *pfs_9 = buffer.data(pfs + 9);
    const auto *pfs_10 = buffer.data(pfs + 10);
    const auto *pfs_11 = buffer.data(pfs + 11);
    const auto *pfs_12 = buffer.data(pfs + 12);
    const auto *pfs_13 = buffer.data(pfs + 13);
    const auto *pfs_14 = buffer.data(pfs + 14);
    const auto *pfs_15 = buffer.data(pfs + 15);
    const auto *pfs_16 = buffer.data(pfs + 16);
    const auto *pfs_17 = buffer.data(pfs + 17);
    const auto *pfs_18 = buffer.data(pfs + 18);
    const auto *pfs_19 = buffer.data(pfs + 19);
    const auto *pfs_20 = buffer.data(pfs + 20);
    const auto *pfs_21 = buffer.data(pfs + 21);
    const auto *pfs_22 = buffer.data(pfs + 22);
    const auto *pfs_23 = buffer.data(pfs + 23);
    const auto *pfs_24 = buffer.data(pfs + 24);
    const auto *pfs_25 = buffer.data(pfs + 25);
    const auto *pfs_26 = buffer.data(pfs + 26);
    const auto *pfs_27 = buffer.data(pfs + 27);
    const auto *pfs_28 = buffer.data(pfs + 28);
    const auto *pfs_29 = buffer.data(pfs + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_y, pc_x, pc_y, pc_z, sfs_0, pdp0_0, \
                         pds_0, pdp1_0, pfs_0, pfs_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfs_0[k]
                 + f_1 * pds_0[k]
                 + f_2 * pc_x[k] * pfs_0[k];

        t_1[k] = f_2 * pc_y[k] * pfs_0[k];

        t_2[k] = f_2 * pc_z[k] * pfs_0[k];

        t_3[k] = pb_y[k] * pdp0_0[k]
                 - f_3 * pc_y[k] * pdp1_0[k];

        t_4[k] = f_0 * pds_0[k]
                 + f_2 * pc_y[k] * pfs_1[k];

        t_5[k] = f_2 * pc_z[k] * pfs_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_z, pc_x, pc_y, pc_z, sfs_3, pdp0_0, pds_0, \
                         pds_3, pdp1_0, pfs_2, pfs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_z[k] * pdp0_0[k]
                 - f_3 * pc_z[k] * pdp1_0[k];

        t_7[k] = f_2 * pc_y[k] * pfs_2[k];

        t_8[k] = f_0 * pds_0[k]
                 + f_2 * pc_z[k] * pfs_2[k];

        t_9[k] = f_0 * sfs_3[k]
                 + f_0 * pds_3[k]
                 + f_2 * pc_x[k] * pfs_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pc_y, pc_z, pdp0_6, pds_1, pds_2, \
                         pdp1_6, pfs_3, pfs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * pds_1[k]
                  + f_2 * pc_y[k] * pfs_3[k];

        t_11[k] = f_2 * pc_z[k] * pfs_3[k];

        t_12[k] = pb_y[k] * pdp0_6[k]
                  - f_3 * pc_y[k] * pdp1_6[k];

        t_13[k] = f_0 * pds_2[k]
                  + f_2 * pc_y[k] * pfs_4[k];

        t_14[k] = f_0 * pds_1[k]
                  + f_2 * pc_z[k] * pfs_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, sfs_5, sfs_6, \
                         pds_2, pds_3, pds_5, pfs_5, pfs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * sfs_5[k]
                  + f_0 * pds_5[k]
                  + f_2 * pc_x[k] * pfs_5[k];

        t_16[k] = f_2 * pc_y[k] * pfs_5[k];

        t_17[k] = f_4 * pds_2[k]
                  + f_2 * pc_z[k] * pfs_5[k];

        t_18[k] = f_0 * sfs_6[k]
                  + f_2 * pc_x[k] * pfs_6[k];

        t_19[k] = f_1 * pds_3[k]
                  + f_2 * pc_y[k] * pfs_6[k];

        t_20[k] = f_2 * pc_z[k] * pfs_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pc_x, pc_y, pc_z, sfs_7, sfs_8, \
                         pds_3, pds_4, pds_5, pfs_7, pfs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * sfs_7[k]
                  + f_2 * pc_x[k] * pfs_7[k];

        t_22[k] = f_4 * pds_4[k]
                  + f_2 * pc_y[k] * pfs_7[k];

        t_23[k] = f_0 * pds_3[k]
                  + f_2 * pc_z[k] * pfs_7[k];

        t_24[k] = f_0 * sfs_8[k]
                  + f_2 * pc_x[k] * pfs_8[k];

        t_25[k] = f_0 * pds_5[k]
                  + f_2 * pc_y[k] * pfs_8[k];

        t_26[k] = f_4 * pds_4[k]
                  + f_2 * pc_z[k] * pfs_8[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pc_x, pc_y, pc_z, sfs_0, sfs_9, \
                         pds_5, pds_6, pfs_9, pfs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * sfs_9[k]
                  + f_2 * pc_x[k] * pfs_9[k];

        t_28[k] = f_2 * pc_y[k] * pfs_9[k];

        t_29[k] = f_1 * pds_5[k]
                  + f_2 * pc_z[k] * pfs_9[k];

        t_30[k] = f_1 * pds_6[k]
                  + f_2 * pc_x[k] * pfs_10[k];

        t_31[k] = f_0 * sfs_0[k]
                  + f_2 * pc_y[k] * pfs_10[k];

        t_32[k] = f_2 * pc_z[k] * pfs_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, sfs_1, sfs_2, \
                         pds_6, pds_7, pds_8, pfs_11, pfs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * pds_7[k]
                  + f_2 * pc_x[k] * pfs_11[k];

        t_34[k] = f_0 * sfs_1[k]
                  + f_0 * pds_6[k]
                  + f_2 * pc_y[k] * pfs_11[k];

        t_35[k] = f_2 * pc_z[k] * pfs_11[k];

        t_36[k] = f_4 * pds_8[k]
                  + f_2 * pc_x[k] * pfs_12[k];

        t_37[k] = f_0 * sfs_2[k]
                  + f_2 * pc_y[k] * pfs_12[k];

        t_38[k] = f_0 * pds_6[k]
                  + f_2 * pc_z[k] * pfs_12[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pb_x, pc_x, pc_z, pdp0_28, pdp0_31, \
                         pds_9, pds_10, pdp1_28, pdp1_31, pfs_13, \
                         pfs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * pds_9[k]
                  + f_2 * pc_x[k] * pfs_13[k];

        t_40[k] = pb_x[k] * pdp0_28[k]
                  - f_3 * pc_x[k] * pdp1_28[k];

        t_41[k] = f_2 * pc_z[k] * pfs_13[k];

        t_42[k] = f_0 * pds_10[k]
                  + f_2 * pc_x[k] * pfs_14[k];

        t_43[k] = pb_x[k] * pdp0_31[k]
                  - f_3 * pc_x[k] * pdp1_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, sfs_5, pds_7, pds_8, \
                         pds_11, pfs_14, pfs_15, pfs_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * pds_7[k]
                  + f_2 * pc_z[k] * pfs_14[k];

        t_45[k] = f_0 * pds_11[k]
                  + f_2 * pc_x[k] * pfs_15[k];

        t_46[k] = f_0 * sfs_5[k]
                  + f_2 * pc_y[k] * pfs_15[k];

        t_47[k] = f_4 * pds_8[k]
                  + f_2 * pc_z[k] * pfs_15[k];

        t_48[k] = f_2 * pc_x[k] * pfs_16[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pb_z, pc_x, pc_y, pc_z, sfs_6, pdp0_28, \
                         pds_9, pdp1_28, pfs_16, pfs_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * sfs_6[k]
                  + f_1 * pds_9[k]
                  + f_2 * pc_y[k] * pfs_16[k];

        t_50[k] = f_2 * pc_z[k] * pfs_16[k];

        t_51[k] = f_2 * pc_x[k] * pfs_17[k];

        t_52[k] = pb_z[k] * pdp0_28[k]
                  - f_3 * pc_z[k] * pdp1_28[k];

        t_53[k] = f_0 * pds_9[k]
                  + f_2 * pc_z[k] * pfs_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sfs_8, sfs_9, \
                         pds_10, pds_11, pfs_18, pfs_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_2 * pc_x[k] * pfs_18[k];

        t_55[k] = f_0 * sfs_8[k]
                  + f_0 * pds_11[k]
                  + f_2 * pc_y[k] * pfs_18[k];

        t_56[k] = f_4 * pds_10[k]
                  + f_2 * pc_z[k] * pfs_18[k];

        t_57[k] = f_2 * pc_x[k] * pfs_19[k];

        t_58[k] = f_0 * sfs_9[k]
                  + f_2 * pc_y[k] * pfs_19[k];

        t_59[k] = f_1 * pds_11[k]
                  + f_2 * pc_z[k] * pfs_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, sfs_0, sfs_1, \
                         pds_12, pds_13, pfs_20, pfs_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * pds_12[k]
                  + f_2 * pc_x[k] * pfs_20[k];

        t_61[k] = f_2 * pc_y[k] * pfs_20[k];

        t_62[k] = f_0 * sfs_0[k]
                  + f_2 * pc_z[k] * pfs_20[k];

        t_63[k] = f_4 * pds_13[k]
                  + f_2 * pc_x[k] * pfs_21[k];

        t_64[k] = f_0 * pds_12[k]
                  + f_2 * pc_y[k] * pfs_21[k];

        t_65[k] = f_0 * sfs_1[k]
                  + f_2 * pc_z[k] * pfs_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pc_x, pc_y, pc_z, sfs_2, pds_12, \
                         pds_13, pds_14, pds_15, pfs_22, pfs_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * pds_14[k]
                  + f_2 * pc_x[k] * pfs_22[k];

        t_67[k] = f_2 * pc_y[k] * pfs_22[k];

        t_68[k] = f_0 * sfs_2[k]
                  + f_0 * pds_12[k]
                  + f_2 * pc_z[k] * pfs_22[k];

        t_69[k] = f_0 * pds_15[k]
                  + f_2 * pc_x[k] * pfs_23[k];

        t_70[k] = f_4 * pds_13[k]
                  + f_2 * pc_y[k] * pfs_23[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pb_x, pc_x, pc_y, pc_z, sfs_3, pdp0_50, \
                         pds_14, pds_16, pdp1_50, pfs_23, pfs_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * sfs_3[k]
                  + f_2 * pc_z[k] * pfs_23[k];

        t_72[k] = f_0 * pds_16[k]
                  + f_2 * pc_x[k] * pfs_24[k];

        t_73[k] = f_0 * pds_14[k]
                  + f_2 * pc_y[k] * pfs_24[k];

        t_74[k] = pb_x[k] * pdp0_50[k]
                  - f_3 * pc_x[k] * pdp1_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pc_x, pc_y, pdp0_53, pds_15, \
                         pds_17, pdp1_53, pfs_25, pfs_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * pds_17[k]
                  + f_2 * pc_x[k] * pfs_25[k];

        t_76[k] = f_2 * pc_y[k] * pfs_25[k];

        t_77[k] = pb_x[k] * pdp0_53[k]
                  - f_3 * pc_x[k] * pdp1_53[k];

        t_78[k] = f_2 * pc_x[k] * pfs_26[k];

        t_79[k] = f_1 * pds_15[k]
                  + f_2 * pc_y[k] * pfs_26[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, sfs_6, sfs_7, pds_15, \
                         pds_16, pfs_26, pfs_27, pfs_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * sfs_6[k]
                  + f_2 * pc_z[k] * pfs_26[k];

        t_81[k] = f_2 * pc_x[k] * pfs_27[k];

        t_82[k] = f_4 * pds_16[k]
                  + f_2 * pc_y[k] * pfs_27[k];

        t_83[k] = f_0 * sfs_7[k]
                  + f_0 * pds_15[k]
                  + f_2 * pc_z[k] * pfs_27[k];

        t_84[k] = f_2 * pc_x[k] * pfs_28[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pb_y, pc_x, pc_y, pc_z, sfs_9, pdp0_53, \
                         pds_17, pdp1_53, pfs_28, pfs_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_0 * pds_17[k]
                  + f_2 * pc_y[k] * pfs_28[k];

        t_86[k] = pb_y[k] * pdp0_53[k]
                  - f_3 * pc_y[k] * pdp1_53[k];

        t_87[k] = f_2 * pc_x[k] * pfs_29[k];

        t_88[k] = f_2 * pc_y[k] * pfs_29[k];

        t_89[k] = f_0 * sfs_9[k]
                  + f_1 * pds_17[k]
                  + f_2 * pc_z[k] * pfs_29[k];
    }
}

}  // namespace simdt3ceri
