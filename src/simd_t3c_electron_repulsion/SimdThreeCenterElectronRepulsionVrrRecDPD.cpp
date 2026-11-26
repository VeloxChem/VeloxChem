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


#include "SimdThreeCenterElectronRepulsionVrrRecDPD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_dpd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppd0,
                                                   const size_t ppp, const size_t ppd1,
                                                   const size_t dsd0, const size_t dsp,
                                                   const size_t dsd1, const size_t dps0,
                                                   const size_t dps1, const size_t dpp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppd0_0 = buffer.data(ppd0 + 0);
    const auto *ppd0_6 = buffer.data(ppd0 + 6);
    const auto *ppd0_12 = buffer.data(ppd0 + 12);
    const auto *ppd0_21 = buffer.data(ppd0 + 21);
    const auto *ppd0_27 = buffer.data(ppd0 + 27);
    const auto *ppd0_29 = buffer.data(ppd0 + 29);
    const auto *ppd0_33 = buffer.data(ppd0 + 33);
    const auto *ppd0_35 = buffer.data(ppd0 + 35);
    const auto *ppd0_36 = buffer.data(ppd0 + 36);
    const auto *ppd0_41 = buffer.data(ppd0 + 41);
    const auto *ppd0_45 = buffer.data(ppd0 + 45);
    const auto *ppd0_47 = buffer.data(ppd0 + 47);
    const auto *ppd0_48 = buffer.data(ppd0 + 48);
    const auto *ppd0_51 = buffer.data(ppd0 + 51);
    const auto *ppd0_53 = buffer.data(ppd0 + 53);

    const auto *ppp_0 = buffer.data(ppp + 0);
    const auto *ppp_1 = buffer.data(ppp + 1);
    const auto *ppp_2 = buffer.data(ppp + 2);
    const auto *ppp_10 = buffer.data(ppp + 10);
    const auto *ppp_12 = buffer.data(ppp + 12);
    const auto *ppp_13 = buffer.data(ppp + 13);
    const auto *ppp_14 = buffer.data(ppp + 14);
    const auto *ppp_16 = buffer.data(ppp + 16);
    const auto *ppp_20 = buffer.data(ppp + 20);
    const auto *ppp_23 = buffer.data(ppp + 23);
    const auto *ppp_24 = buffer.data(ppp + 24);
    const auto *ppp_25 = buffer.data(ppp + 25);
    const auto *ppp_26 = buffer.data(ppp + 26);

    const auto *ppd1_0 = buffer.data(ppd1 + 0);
    const auto *ppd1_6 = buffer.data(ppd1 + 6);
    const auto *ppd1_12 = buffer.data(ppd1 + 12);
    const auto *ppd1_21 = buffer.data(ppd1 + 21);
    const auto *ppd1_27 = buffer.data(ppd1 + 27);
    const auto *ppd1_29 = buffer.data(ppd1 + 29);
    const auto *ppd1_33 = buffer.data(ppd1 + 33);
    const auto *ppd1_35 = buffer.data(ppd1 + 35);
    const auto *ppd1_36 = buffer.data(ppd1 + 36);
    const auto *ppd1_41 = buffer.data(ppd1 + 41);
    const auto *ppd1_45 = buffer.data(ppd1 + 45);
    const auto *ppd1_47 = buffer.data(ppd1 + 47);
    const auto *ppd1_48 = buffer.data(ppd1 + 48);
    const auto *ppd1_51 = buffer.data(ppd1 + 51);
    const auto *ppd1_53 = buffer.data(ppd1 + 53);

    const auto *dsd0_0 = buffer.data(dsd0 + 0);
    const auto *dsd0_3 = buffer.data(dsd0 + 3);
    const auto *dsd0_5 = buffer.data(dsd0 + 5);
    const auto *dsd0_18 = buffer.data(dsd0 + 18);
    const auto *dsd0_21 = buffer.data(dsd0 + 21);
    const auto *dsd0_23 = buffer.data(dsd0 + 23);
    const auto *dsd0_30 = buffer.data(dsd0 + 30);
    const auto *dsd0_33 = buffer.data(dsd0 + 33);
    const auto *dsd0_35 = buffer.data(dsd0 + 35);

    const auto *dsp_0 = buffer.data(dsp + 0);
    const auto *dsp_1 = buffer.data(dsp + 1);
    const auto *dsp_2 = buffer.data(dsp + 2);
    const auto *dsp_3 = buffer.data(dsp + 3);
    const auto *dsp_4 = buffer.data(dsp + 4);
    const auto *dsp_6 = buffer.data(dsp + 6);
    const auto *dsp_8 = buffer.data(dsp + 8);
    const auto *dsp_9 = buffer.data(dsp + 9);
    const auto *dsp_10 = buffer.data(dsp + 10);
    const auto *dsp_11 = buffer.data(dsp + 11);
    const auto *dsp_13 = buffer.data(dsp + 13);
    const auto *dsp_14 = buffer.data(dsp + 14);
    const auto *dsp_15 = buffer.data(dsp + 15);
    const auto *dsp_16 = buffer.data(dsp + 16);
    const auto *dsp_17 = buffer.data(dsp + 17);

    const auto *dsd1_0 = buffer.data(dsd1 + 0);
    const auto *dsd1_3 = buffer.data(dsd1 + 3);
    const auto *dsd1_5 = buffer.data(dsd1 + 5);
    const auto *dsd1_18 = buffer.data(dsd1 + 18);
    const auto *dsd1_21 = buffer.data(dsd1 + 21);
    const auto *dsd1_23 = buffer.data(dsd1 + 23);
    const auto *dsd1_30 = buffer.data(dsd1 + 30);
    const auto *dsd1_33 = buffer.data(dsd1 + 33);
    const auto *dsd1_35 = buffer.data(dsd1 + 35);

    const auto *dps0_0 = buffer.data(dps0 + 0);
    const auto *dps0_3 = buffer.data(dps0 + 3);
    const auto *dps0_4 = buffer.data(dps0 + 4);
    const auto *dps0_6 = buffer.data(dps0 + 6);
    const auto *dps0_8 = buffer.data(dps0 + 8);
    const auto *dps0_10 = buffer.data(dps0 + 10);
    const auto *dps0_13 = buffer.data(dps0 + 13);
    const auto *dps0_14 = buffer.data(dps0 + 14);
    const auto *dps0_17 = buffer.data(dps0 + 17);

    const auto *dps1_0 = buffer.data(dps1 + 0);
    const auto *dps1_3 = buffer.data(dps1 + 3);
    const auto *dps1_4 = buffer.data(dps1 + 4);
    const auto *dps1_6 = buffer.data(dps1 + 6);
    const auto *dps1_8 = buffer.data(dps1 + 8);
    const auto *dps1_10 = buffer.data(dps1 + 10);
    const auto *dps1_13 = buffer.data(dps1 + 13);
    const auto *dps1_14 = buffer.data(dps1 + 14);
    const auto *dps1_17 = buffer.data(dps1 + 17);

    const auto *dpp_0 = buffer.data(dpp + 0);
    const auto *dpp_1 = buffer.data(dpp + 1);
    const auto *dpp_2 = buffer.data(dpp + 2);
    const auto *dpp_3 = buffer.data(dpp + 3);
    const auto *dpp_5 = buffer.data(dpp + 5);
    const auto *dpp_6 = buffer.data(dpp + 6);
    const auto *dpp_8 = buffer.data(dpp + 8);
    const auto *dpp_9 = buffer.data(dpp + 9);
    const auto *dpp_10 = buffer.data(dpp + 10);
    const auto *dpp_11 = buffer.data(dpp + 11);
    const auto *dpp_12 = buffer.data(dpp + 12);
    const auto *dpp_13 = buffer.data(dpp + 13);
    const auto *dpp_15 = buffer.data(dpp + 15);
    const auto *dpp_16 = buffer.data(dpp + 16);
    const auto *dpp_18 = buffer.data(dpp + 18);
    const auto *dpp_19 = buffer.data(dpp + 19);
    const auto *dpp_20 = buffer.data(dpp + 20);
    const auto *dpp_21 = buffer.data(dpp + 21);
    const auto *dpp_23 = buffer.data(dpp + 23);
    const auto *dpp_24 = buffer.data(dpp + 24);
    const auto *dpp_26 = buffer.data(dpp + 26);
    const auto *dpp_28 = buffer.data(dpp + 28);
    const auto *dpp_29 = buffer.data(dpp + 29);
    const auto *dpp_30 = buffer.data(dpp + 30);
    const auto *dpp_31 = buffer.data(dpp + 31);
    const auto *dpp_32 = buffer.data(dpp + 32);
    const auto *dpp_34 = buffer.data(dpp + 34);
    const auto *dpp_35 = buffer.data(dpp + 35);
    const auto *dpp_37 = buffer.data(dpp + 37);
    const auto *dpp_38 = buffer.data(dpp + 38);
    const auto *dpp_39 = buffer.data(dpp + 39);
    const auto *dpp_40 = buffer.data(dpp + 40);
    const auto *dpp_41 = buffer.data(dpp + 41);
    const auto *dpp_43 = buffer.data(dpp + 43);
    const auto *dpp_44 = buffer.data(dpp + 44);
    const auto *dpp_46 = buffer.data(dpp + 46);
    const auto *dpp_47 = buffer.data(dpp + 47);
    const auto *dpp_49 = buffer.data(dpp + 49);
    const auto *dpp_50 = buffer.data(dpp + 50);
    const auto *dpp_51 = buffer.data(dpp + 51);
    const auto *dpp_52 = buffer.data(dpp + 52);
    const auto *dpp_53 = buffer.data(dpp + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ppp_0, dsp_0, dps0_0, \
                         dps1_0, dpp_0, dpp_1, dpp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ppp_0[k]
                 + f_1 * dsp_0[k]
                 + f_2 * dps0_0[k]
                 - f_3 * dps1_0[k]
                 + f_4 * pc_x[k] * dpp_0[k];

        t_1[k] = f_4 * pc_y[k] * dpp_0[k];

        t_2[k] = f_4 * pc_z[k] * dpp_0[k];

        t_3[k] = f_2 * dps0_0[k]
                 - f_3 * dps1_0[k]
                 + f_4 * pc_y[k] * dpp_1[k];

        t_4[k] = f_4 * pc_y[k] * dpp_2[k];

        t_5[k] = f_2 * dps0_0[k]
                 - f_3 * dps1_0[k]
                 + f_4 * pc_z[k] * dpp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pc_y, pc_z, dsd0_0, dsd0_3, dsp_0, dsp_1, \
                         dsd1_0, dsd1_3, dpp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * dsd0_0[k]
                 - f_5 * pc_y[k] * dsd1_0[k];

        t_7[k] = f_1 * dsp_0[k]
                 + f_4 * pc_y[k] * dpp_3[k];

        t_8[k] = f_4 * pc_z[k] * dpp_3[k];

        t_9[k] = pb_y[k] * dsd0_3[k]
                 + f_0 * dsp_1[k]
                 - f_5 * pc_y[k] * dsd1_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, pc_y, pc_z, dsd0_0, dsd0_5, \
                         dsp_2, dsd1_0, dsd1_5, dpp_5, dpp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dsp_2[k]
                  + f_4 * pc_y[k] * dpp_5[k];

        t_11[k] = pb_y[k] * dsd0_5[k]
                  - f_5 * pc_y[k] * dsd1_5[k];

        t_12[k] = pb_z[k] * dsd0_0[k]
                  - f_5 * pc_z[k] * dsd1_0[k];

        t_13[k] = f_4 * pc_y[k] * dpp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_z, pc_y, pc_z, dsd0_3, dsd0_5, dsp_0, \
                         dsp_2, dsd1_3, dsd1_5, dpp_6, dpp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * dsp_0[k]
                  + f_4 * pc_z[k] * dpp_6[k];

        t_15[k] = pb_z[k] * dsd0_3[k]
                  - f_5 * pc_z[k] * dsd1_3[k];

        t_16[k] = f_4 * pc_y[k] * dpp_8[k];

        t_17[k] = pb_z[k] * dsd0_5[k]
                  + f_0 * dsp_2[k]
                  - f_5 * pc_z[k] * dsd1_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pc_x, pc_y, pc_z, ppd0_0, ppp_10, ppd1_0, \
                         dsp_4, dpp_9, dpp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_y[k] * ppd0_0[k]
                  - f_5 * pc_y[k] * ppd1_0[k];

        t_19[k] = f_1 * ppp_10[k]
                  + f_1 * dsp_4[k]
                  + f_4 * pc_x[k] * dpp_10[k];

        t_20[k] = f_4 * pc_z[k] * dpp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, ppp_1, ppp_12, dps0_3, \
                         dps0_4, dps1_3, dps1_4, dpp_10, dpp_11, \
                         dpp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * ppp_1[k]
                  + f_2 * dps0_3[k]
                  - f_3 * dps1_3[k]
                  + f_4 * pc_y[k] * dpp_10[k];

        t_22[k] = f_4 * pc_z[k] * dpp_10[k];

        t_23[k] = f_2 * dps0_3[k]
                  - f_3 * dps1_3[k]
                  + f_4 * pc_z[k] * dpp_11[k];

        t_24[k] = f_1 * ppp_12[k]
                  + f_2 * dps0_4[k]
                  - f_3 * dps1_4[k]
                  + f_4 * pc_x[k] * dpp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pc_x, pc_z, ppd0_27, ppd0_29, \
                         ppp_13, ppd1_27, ppd1_29, dpp_12, dpp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ppp_13[k]
                  + f_4 * pc_x[k] * dpp_13[k];

        t_26[k] = f_4 * pc_z[k] * dpp_12[k];

        t_27[k] = pa_x[k] * ppd0_27[k]
                  - f_5 * pc_x[k] * ppd1_27[k];

        t_28[k] = f_4 * pc_z[k] * dpp_13[k];

        t_29[k] = pa_x[k] * ppd0_29[k]
                  - f_5 * pc_x[k] * ppd1_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pc_x, pc_y, pc_z, ppd0_12, ppp_16, ppd1_12, \
                         dsp_3, dpp_15, dpp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * ppd0_12[k]
                  - f_5 * pc_y[k] * ppd1_12[k];

        t_31[k] = f_1 * ppp_16[k]
                  + f_4 * pc_x[k] * dpp_16[k];

        t_32[k] = f_1 * dsp_3[k]
                  + f_4 * pc_z[k] * dpp_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_x, pa_z, pc_x, pc_z, ppd0_0, ppd0_33, \
                         ppd0_35, ppd1_0, ppd1_33, ppd1_35, dsp_4, \
                         dpp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * ppd0_33[k]
                  - f_5 * pc_x[k] * ppd1_33[k];

        t_34[k] = f_1 * dsp_4[k]
                  + f_4 * pc_z[k] * dpp_16[k];

        t_35[k] = pa_x[k] * ppd0_35[k]
                  - f_5 * pc_x[k] * ppd1_35[k];

        t_36[k] = pa_z[k] * ppd0_0[k]
                  - f_5 * pc_z[k] * ppd1_0[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, ppp_2, ppp_20, dsp_8, \
                         dps0_6, dps1_6, dpp_18, dpp_19, dpp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_4 * pc_y[k] * dpp_18[k];

        t_38[k] = f_1 * ppp_20[k]
                  + f_1 * dsp_8[k]
                  + f_4 * pc_x[k] * dpp_20[k];

        t_39[k] = f_2 * dps0_6[k]
                  - f_3 * dps1_6[k]
                  + f_4 * pc_y[k] * dpp_19[k];

        t_40[k] = f_4 * pc_y[k] * dpp_20[k];

        t_41[k] = f_1 * ppp_2[k]
                  + f_2 * dps0_6[k]
                  - f_3 * dps1_6[k]
                  + f_4 * pc_z[k] * dpp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_z, pc_x, pc_y, pc_z, ppd0_6, ppp_23, ppd1_6, \
                         dsp_6, dpp_21, dpp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * ppd0_6[k]
                  - f_5 * pc_z[k] * ppd1_6[k];

        t_43[k] = f_1 * dsp_6[k]
                  + f_4 * pc_y[k] * dpp_21[k];

        t_44[k] = f_1 * ppp_23[k]
                  + f_4 * pc_x[k] * dpp_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_x, pc_x, pc_y, ppd0_45, ppd0_47, ppd1_45, \
                         ppd1_47, dsp_8, dpp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_x[k] * ppd0_45[k]
                  - f_5 * pc_x[k] * ppd1_45[k];

        t_46[k] = f_1 * dsp_8[k]
                  + f_4 * pc_y[k] * dpp_23[k];

        t_47[k] = pa_x[k] * ppd0_47[k]
                  - f_5 * pc_x[k] * ppd1_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pc_x, pc_y, ppd0_51, ppp_24, \
                         ppp_26, ppd1_51, dps0_8, dps1_8, dpp_24, \
                         dpp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ppp_24[k]
                  + f_2 * dps0_8[k]
                  - f_3 * dps1_8[k]
                  + f_4 * pc_x[k] * dpp_24[k];

        t_49[k] = f_4 * pc_y[k] * dpp_24[k];

        t_50[k] = f_1 * ppp_26[k]
                  + f_4 * pc_x[k] * dpp_26[k];

        t_51[k] = pa_x[k] * ppd0_51[k]
                  - f_5 * pc_x[k] * ppd1_51[k];

        t_52[k] = f_4 * pc_y[k] * dpp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pb_x, pc_x, ppd0_53, ppd1_53, dsd0_18, \
                         dsp_9, dsp_10, dsp_11, dsd1_18, dpp_28, \
                         dpp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * ppd0_53[k]
                  - f_5 * pc_x[k] * ppd1_53[k];

        t_54[k] = pb_x[k] * dsd0_18[k]
                  + f_0 * dsp_9[k]
                  - f_5 * pc_x[k] * dsd1_18[k];

        t_55[k] = f_1 * dsp_10[k]
                  + f_4 * pc_x[k] * dpp_28[k];

        t_56[k] = f_1 * dsp_11[k]
                  + f_4 * pc_x[k] * dpp_29[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_x, pc_x, pc_z, dsd0_21, dsd0_23, dsd1_21, \
                         dsd1_23, dps0_10, dps1_10, dpp_28, dpp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_x[k] * dsd0_21[k]
                  - f_5 * pc_x[k] * dsd1_21[k];

        t_58[k] = f_4 * pc_z[k] * dpp_28[k];

        t_59[k] = pb_x[k] * dsd0_23[k]
                  - f_5 * pc_x[k] * dsd1_23[k];

        t_60[k] = f_2 * dps0_10[k]
                  - f_3 * dps1_10[k]
                  + f_4 * pc_x[k] * dpp_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, ppp_13, dsp_10, \
                         dps0_10, dps1_10, dpp_31, dpp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_4 * pc_x[k] * dpp_31[k];

        t_62[k] = f_4 * pc_x[k] * dpp_32[k];

        t_63[k] = f_0 * ppp_13[k]
                  + f_1 * dsp_10[k]
                  + f_2 * dps0_10[k]
                  - f_3 * dps1_10[k]
                  + f_4 * pc_y[k] * dpp_31[k];

        t_64[k] = f_4 * pc_z[k] * dpp_31[k];

        t_65[k] = f_2 * dps0_10[k]
                  - f_3 * dps1_10[k]
                  + f_4 * pc_z[k] * dpp_32[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pb_z, pc_x, pc_z, dsd0_18, dsd0_21, \
                         dsp_10, dsd1_18, dsd1_21, dpp_34, dpp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * dsd0_18[k]
                  - f_5 * pc_z[k] * dsd1_18[k];

        t_67[k] = f_4 * pc_x[k] * dpp_34[k];

        t_68[k] = f_4 * pc_x[k] * dpp_35[k];

        t_69[k] = pb_z[k] * dsd0_21[k]
                  - f_5 * pc_z[k] * dsd1_21[k];

        t_70[k] = f_1 * dsp_10[k]
                  + f_4 * pc_z[k] * dpp_34[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_y, pb_z, pc_x, pc_y, pc_z, ppd0_36, ppd1_36, \
                         dsd0_23, dsp_11, dsp_13, dsd1_23, dpp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_z[k] * dsd0_23[k]
                  + f_0 * dsp_11[k]
                  - f_5 * pc_z[k] * dsd1_23[k];

        t_72[k] = pa_y[k] * ppd0_36[k]
                  - f_5 * pc_y[k] * ppd1_36[k];

        t_73[k] = f_1 * dsp_13[k]
                  + f_4 * pc_x[k] * dpp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_y, pa_z, pc_x, pc_y, pc_z, ppd0_21, \
                         ppd0_41, ppp_20, ppd1_21, ppd1_41, dsp_14, \
                         dpp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * dsp_14[k]
                  + f_4 * pc_x[k] * dpp_38[k];

        t_75[k] = pa_z[k] * ppd0_21[k]
                  - f_5 * pc_z[k] * ppd1_21[k];

        t_76[k] = f_1 * ppp_20[k]
                  + f_4 * pc_y[k] * dpp_38[k];

        t_77[k] = pa_y[k] * ppd0_41[k]
                  - f_5 * pc_y[k] * ppd1_41[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_z, pc_x, pc_z, ppd0_27, ppd1_27, dps0_13, \
                         dps1_13, dpp_39, dpp_40, dpp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_2 * dps0_13[k]
                  - f_3 * dps1_13[k]
                  + f_4 * pc_x[k] * dpp_39[k];

        t_79[k] = f_4 * pc_x[k] * dpp_40[k];

        t_80[k] = f_4 * pc_x[k] * dpp_41[k];

        t_81[k] = pa_z[k] * ppd0_27[k]
                  - f_5 * pc_z[k] * ppd1_27[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pc_y, pc_z, ppd0_48, ppp_14, ppp_23, ppd1_48, \
                         dsp_14, dps0_13, dps1_13, dpp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_1 * ppp_23[k]
                  + f_1 * dsp_14[k]
                  + f_4 * pc_y[k] * dpp_41[k];

        t_83[k] = f_1 * ppp_14[k]
                  + f_2 * dps0_13[k]
                  - f_3 * dps1_13[k]
                  + f_4 * pc_z[k] * dpp_41[k];

        t_84[k] = pa_y[k] * ppd0_48[k]
                  - f_5 * pc_y[k] * ppd1_48[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_y, pc_x, pc_y, ppd0_53, ppp_25, \
                         ppp_26, ppd1_53, dps0_14, dps1_14, dpp_43, \
                         dpp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_4 * pc_x[k] * dpp_43[k];

        t_86[k] = f_4 * pc_x[k] * dpp_44[k];

        t_87[k] = f_1 * ppp_25[k]
                  + f_2 * dps0_14[k]
                  - f_3 * dps1_14[k]
                  + f_4 * pc_y[k] * dpp_43[k];

        t_88[k] = f_1 * ppp_26[k]
                  + f_4 * pc_y[k] * dpp_44[k];

        t_89[k] = pa_y[k] * ppd0_53[k]
                  - f_5 * pc_y[k] * ppd1_53[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pc_x, dsd0_30, dsd0_33, dsp_15, dsp_16, \
                         dsp_17, dsd1_30, dsd1_33, dpp_46, dpp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_x[k] * dsd0_30[k]
                  + f_0 * dsp_15[k]
                  - f_5 * pc_x[k] * dsd1_30[k];

        t_91[k] = f_1 * dsp_16[k]
                  + f_4 * pc_x[k] * dpp_46[k];

        t_92[k] = f_1 * dsp_17[k]
                  + f_4 * pc_x[k] * dpp_47[k];

        t_93[k] = pb_x[k] * dsd0_33[k]
                  - f_5 * pc_x[k] * dsd1_33[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, pc_x, pc_y, dsd0_30, \
                         dsd0_35, dsd1_30, dsd1_35, dpp_47, dpp_49, \
                         dpp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_4 * pc_y[k] * dpp_47[k];

        t_95[k] = pb_x[k] * dsd0_35[k]
                  - f_5 * pc_x[k] * dsd1_35[k];

        t_96[k] = pb_y[k] * dsd0_30[k]
                  - f_5 * pc_y[k] * dsd1_30[k];

        t_97[k] = f_4 * pc_x[k] * dpp_49[k];

        t_98[k] = f_4 * pc_x[k] * dpp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_y, pc_y, dsd0_33, dsd0_35, dsp_16, dsp_17, \
                         dsd1_33, dsd1_35, dpp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pb_y[k] * dsd0_33[k]
                  + f_0 * dsp_16[k]
                  - f_5 * pc_y[k] * dsd1_33[k];

        t_100[k] = f_1 * dsp_17[k]
                   + f_4 * pc_y[k] * dpp_50[k];

        t_101[k] = pb_y[k] * dsd0_35[k]
                   - f_5 * pc_y[k] * dsd1_35[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, ppp_26, \
                         dsp_17, dps0_17, dps1_17, dpp_51, dpp_52, \
                         dpp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_2 * dps0_17[k]
                   - f_3 * dps1_17[k]
                   + f_4 * pc_x[k] * dpp_51[k];

        t_103[k] = f_4 * pc_x[k] * dpp_52[k];

        t_104[k] = f_4 * pc_x[k] * dpp_53[k];

        t_105[k] = f_2 * dps0_17[k]
                   - f_3 * dps1_17[k]
                   + f_4 * pc_y[k] * dpp_52[k];

        t_106[k] = f_4 * pc_y[k] * dpp_53[k];

        t_107[k] = f_0 * ppp_26[k]
                   + f_1 * dsp_17[k]
                   + f_2 * dps0_17[k]
                   - f_3 * dps1_17[k]
                   + f_4 * pc_z[k] * dpp_53[k];
    }
}

}  // namespace simdt3ceri
