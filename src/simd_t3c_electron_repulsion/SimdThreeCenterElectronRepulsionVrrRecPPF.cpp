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


#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ppf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t spf0,
                                                   const size_t spd, const size_t spf1,
                                                   const size_t psf0, const size_t psd,
                                                   const size_t psf1, const size_t ppp0,
                                                   const size_t ppp1, const size_t ppd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spf0_0 = buffer.data(spf0 + 0);
    const auto *spf0_1 = buffer.data(spf0 + 1);
    const auto *spf0_2 = buffer.data(spf0 + 2);
    const auto *spf0_6 = buffer.data(spf0 + 6);
    const auto *spf0_9 = buffer.data(spf0 + 9);
    const auto *spf0_10 = buffer.data(spf0 + 10);
    const auto *spf0_16 = buffer.data(spf0 + 16);
    const auto *spf0_17 = buffer.data(spf0 + 17);
    const auto *spf0_19 = buffer.data(spf0 + 19);
    const auto *spf0_20 = buffer.data(spf0 + 20);
    const auto *spf0_26 = buffer.data(spf0 + 26);
    const auto *spf0_29 = buffer.data(spf0 + 29);

    const auto *spd_0 = buffer.data(spd + 0);
    const auto *spd_3 = buffer.data(spd + 3);
    const auto *spd_5 = buffer.data(spd + 5);
    const auto *spd_9 = buffer.data(spd + 9);
    const auto *spd_11 = buffer.data(spd + 11);
    const auto *spd_15 = buffer.data(spd + 15);
    const auto *spd_17 = buffer.data(spd + 17);

    const auto *spf1_0 = buffer.data(spf1 + 0);
    const auto *spf1_1 = buffer.data(spf1 + 1);
    const auto *spf1_2 = buffer.data(spf1 + 2);
    const auto *spf1_6 = buffer.data(spf1 + 6);
    const auto *spf1_9 = buffer.data(spf1 + 9);
    const auto *spf1_10 = buffer.data(spf1 + 10);
    const auto *spf1_16 = buffer.data(spf1 + 16);
    const auto *spf1_17 = buffer.data(spf1 + 17);
    const auto *spf1_19 = buffer.data(spf1 + 19);
    const auto *spf1_20 = buffer.data(spf1 + 20);
    const auto *spf1_26 = buffer.data(spf1 + 26);
    const auto *spf1_29 = buffer.data(spf1 + 29);

    const auto *psf0_0 = buffer.data(psf0 + 0);
    const auto *psf0_11 = buffer.data(psf0 + 11);
    const auto *psf0_16 = buffer.data(psf0 + 16);
    const auto *psf0_22 = buffer.data(psf0 + 22);
    const auto *psf0_27 = buffer.data(psf0 + 27);
    const auto *psf0_29 = buffer.data(psf0 + 29);

    const auto *psd_0 = buffer.data(psd + 0);
    const auto *psd_2 = buffer.data(psd + 2);
    const auto *psd_3 = buffer.data(psd + 3);
    const auto *psd_5 = buffer.data(psd + 5);
    const auto *psd_6 = buffer.data(psd + 6);
    const auto *psd_9 = buffer.data(psd + 9);
    const auto *psd_10 = buffer.data(psd + 10);
    const auto *psd_11 = buffer.data(psd + 11);
    const auto *psd_12 = buffer.data(psd + 12);
    const auto *psd_15 = buffer.data(psd + 15);
    const auto *psd_16 = buffer.data(psd + 16);
    const auto *psd_17 = buffer.data(psd + 17);

    const auto *psf1_0 = buffer.data(psf1 + 0);
    const auto *psf1_11 = buffer.data(psf1 + 11);
    const auto *psf1_16 = buffer.data(psf1 + 16);
    const auto *psf1_22 = buffer.data(psf1 + 22);
    const auto *psf1_27 = buffer.data(psf1 + 27);
    const auto *psf1_29 = buffer.data(psf1 + 29);

    const auto *ppp0_0 = buffer.data(ppp0 + 0);
    const auto *ppp0_1 = buffer.data(ppp0 + 1);
    const auto *ppp0_2 = buffer.data(ppp0 + 2);
    const auto *ppp0_12 = buffer.data(ppp0 + 12);
    const auto *ppp0_13 = buffer.data(ppp0 + 13);
    const auto *ppp0_14 = buffer.data(ppp0 + 14);
    const auto *ppp0_24 = buffer.data(ppp0 + 24);
    const auto *ppp0_25 = buffer.data(ppp0 + 25);
    const auto *ppp0_26 = buffer.data(ppp0 + 26);

    const auto *ppp1_0 = buffer.data(ppp1 + 0);
    const auto *ppp1_1 = buffer.data(ppp1 + 1);
    const auto *ppp1_2 = buffer.data(ppp1 + 2);
    const auto *ppp1_12 = buffer.data(ppp1 + 12);
    const auto *ppp1_13 = buffer.data(ppp1 + 13);
    const auto *ppp1_14 = buffer.data(ppp1 + 14);
    const auto *ppp1_24 = buffer.data(ppp1 + 24);
    const auto *ppp1_25 = buffer.data(ppp1 + 25);
    const auto *ppp1_26 = buffer.data(ppp1 + 26);

    const auto *ppd_0 = buffer.data(ppd + 0);
    const auto *ppd_2 = buffer.data(ppd + 2);
    const auto *ppd_3 = buffer.data(ppd + 3);
    const auto *ppd_5 = buffer.data(ppd + 5);
    const auto *ppd_6 = buffer.data(ppd + 6);
    const auto *ppd_8 = buffer.data(ppd + 8);
    const auto *ppd_9 = buffer.data(ppd + 9);
    const auto *ppd_11 = buffer.data(ppd + 11);
    const auto *ppd_12 = buffer.data(ppd + 12);
    const auto *ppd_14 = buffer.data(ppd + 14);
    const auto *ppd_15 = buffer.data(ppd + 15);
    const auto *ppd_17 = buffer.data(ppd + 17);
    const auto *ppd_18 = buffer.data(ppd + 18);
    const auto *ppd_21 = buffer.data(ppd + 21);
    const auto *ppd_22 = buffer.data(ppd + 22);
    const auto *ppd_23 = buffer.data(ppd + 23);
    const auto *ppd_24 = buffer.data(ppd + 24);
    const auto *ppd_25 = buffer.data(ppd + 25);
    const auto *ppd_27 = buffer.data(ppd + 27);
    const auto *ppd_28 = buffer.data(ppd + 28);
    const auto *ppd_29 = buffer.data(ppd + 29);
    const auto *ppd_30 = buffer.data(ppd + 30);
    const auto *ppd_33 = buffer.data(ppd + 33);
    const auto *ppd_34 = buffer.data(ppd + 34);
    const auto *ppd_35 = buffer.data(ppd + 35);
    const auto *ppd_36 = buffer.data(ppd + 36);
    const auto *ppd_39 = buffer.data(ppd + 39);
    const auto *ppd_40 = buffer.data(ppd + 40);
    const auto *ppd_41 = buffer.data(ppd + 41);
    const auto *ppd_42 = buffer.data(ppd + 42);
    const auto *ppd_45 = buffer.data(ppd + 45);
    const auto *ppd_46 = buffer.data(ppd + 46);
    const auto *ppd_47 = buffer.data(ppd + 47);
    const auto *ppd_48 = buffer.data(ppd + 48);
    const auto *ppd_50 = buffer.data(ppd + 50);
    const auto *ppd_51 = buffer.data(ppd + 51);
    const auto *ppd_52 = buffer.data(ppd + 52);
    const auto *ppd_53 = buffer.data(ppd + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, spd_0, spd_3, psd_0, psd_3, \
                         ppp0_0, ppp1_0, ppd_0, ppd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spd_0[k]
                 + f_0 * psd_0[k]
                 + f_1 * ppp0_0[k]
                 - f_2 * ppp1_0[k]
                 + f_3 * pc_x[k] * ppd_0[k];

        t_1[k] = f_3 * pc_y[k] * ppd_0[k];

        t_2[k] = f_3 * pc_z[k] * ppd_0[k];

        t_3[k] = f_0 * spd_3[k]
                 + f_0 * psd_3[k]
                 + f_3 * pc_x[k] * ppd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, spd_5, psd_5, ppp0_1, \
                         ppp1_1, ppd_2, ppd_3, ppd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * ppd_2[k];

        t_5[k] = f_0 * spd_5[k]
                 + f_0 * psd_5[k]
                 + f_3 * pc_x[k] * ppd_5[k];

        t_6[k] = f_1 * ppp0_1[k]
                 - f_2 * ppp1_1[k]
                 + f_3 * pc_y[k] * ppd_3[k];

        t_7[k] = f_3 * pc_z[k] * ppd_3[k];

        t_8[k] = f_3 * pc_y[k] * ppd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, psf0_0, psd_0, psf1_0, \
                         ppp0_2, ppp1_2, ppd_5, ppd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * ppp0_2[k]
                 - f_2 * ppp1_2[k]
                 + f_3 * pc_z[k] * ppd_5[k];

        t_10[k] = pb_y[k] * psf0_0[k]
                  - f_4 * pc_y[k] * psf1_0[k];

        t_11[k] = f_0 * psd_0[k]
                  + f_3 * pc_y[k] * ppd_6[k];

        t_12[k] = f_3 * pc_z[k] * ppd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pc_x, pc_y, spf0_16, spd_9, spd_11, \
                         spf1_16, psd_2, ppd_8, ppd_9, ppd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * spd_9[k]
                  + f_3 * pc_x[k] * ppd_9[k];

        t_14[k] = f_0 * psd_2[k]
                  + f_3 * pc_y[k] * ppd_8[k];

        t_15[k] = f_0 * spd_11[k]
                  + f_3 * pc_x[k] * ppd_11[k];

        t_16[k] = pa_x[k] * spf0_16[k]
                  - f_4 * pc_x[k] * spf1_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pb_z, pc_x, pc_y, pc_z, spf0_19, \
                         spf1_19, psf0_0, psd_5, psf1_0, ppd_9, \
                         ppd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * ppd_9[k];

        t_18[k] = f_0 * psd_5[k]
                  + f_3 * pc_y[k] * ppd_11[k];

        t_19[k] = pa_x[k] * spf0_19[k]
                  - f_4 * pc_x[k] * spf1_19[k];

        t_20[k] = pb_z[k] * psf0_0[k]
                  - f_4 * pc_z[k] * psf1_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, spd_15, spd_17, \
                         psd_0, ppd_12, ppd_14, ppd_15, ppd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * pc_y[k] * ppd_12[k];

        t_22[k] = f_0 * psd_0[k]
                  + f_3 * pc_z[k] * ppd_12[k];

        t_23[k] = f_0 * spd_15[k]
                  + f_3 * pc_x[k] * ppd_15[k];

        t_24[k] = f_3 * pc_y[k] * ppd_14[k];

        t_25[k] = f_0 * spd_17[k]
                  + f_3 * pc_x[k] * ppd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pc_x, pc_y, pc_z, spf0_26, spf0_29, \
                         spf1_26, spf1_29, psd_3, ppd_15, ppd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * spf0_26[k]
                  - f_4 * pc_x[k] * spf1_26[k];

        t_27[k] = f_0 * psd_3[k]
                  + f_3 * pc_z[k] * ppd_15[k];

        t_28[k] = f_3 * pc_y[k] * ppd_17[k];

        t_29[k] = pa_x[k] * spf0_29[k]
                  - f_4 * pc_x[k] * spf1_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_x, pc_y, pc_z, spf0_0, spf0_1, \
                         spd_0, spf1_0, spf1_1, psd_9, ppd_18, ppd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * spf0_0[k]
                  - f_4 * pc_y[k] * spf1_0[k];

        t_31[k] = pa_y[k] * spf0_1[k]
                  + f_0 * spd_0[k]
                  - f_4 * pc_y[k] * spf1_1[k];

        t_32[k] = f_3 * pc_z[k] * ppd_18[k];

        t_33[k] = f_0 * psd_9[k]
                  + f_3 * pc_x[k] * ppd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, pc_x, pc_z, psf0_16, psd_10, psd_11, \
                         psf1_16, ppd_21, ppd_22, ppd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * psd_10[k]
                  + f_3 * pc_x[k] * ppd_22[k];

        t_35[k] = f_0 * psd_11[k]
                  + f_3 * pc_x[k] * ppd_23[k];

        t_36[k] = pb_x[k] * psf0_16[k]
                  - f_4 * pc_x[k] * psf1_16[k];

        t_37[k] = f_3 * pc_z[k] * ppd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pc_x, pc_y, spf0_9, spd_5, spf1_9, ppp0_12, \
                         ppp1_12, ppd_23, ppd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * spd_5[k]
                  + f_3 * pc_y[k] * ppd_23[k];

        t_39[k] = pa_y[k] * spf0_9[k]
                  - f_4 * pc_y[k] * spf1_9[k];

        t_40[k] = f_1 * ppp0_12[k]
                  - f_2 * ppp1_12[k]
                  + f_3 * pc_x[k] * ppd_24[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pc_x, pc_z, ppp0_13, ppp1_13, ppd_24, \
                         ppd_25, ppd_27, ppd_28, ppd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * ppp0_13[k]
                  - f_6 * ppp1_13[k]
                  + f_3 * pc_x[k] * ppd_25[k];

        t_42[k] = f_3 * pc_z[k] * ppd_24[k];

        t_43[k] = f_3 * pc_x[k] * ppd_27[k];

        t_44[k] = f_3 * pc_x[k] * ppd_28[k];

        t_45[k] = f_3 * pc_x[k] * ppd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_y, pc_z, spd_9, spd_11, psd_9, psd_11, \
                         ppp0_13, ppp0_14, ppp1_13, ppp1_14, ppd_27, \
                         ppd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_0 * spd_9[k]
                  + f_0 * psd_9[k]
                  + f_1 * ppp0_13[k]
                  - f_2 * ppp1_13[k]
                  + f_3 * pc_y[k] * ppd_27[k];

        t_47[k] = f_3 * pc_z[k] * ppd_27[k];

        t_48[k] = f_0 * spd_11[k]
                  + f_0 * psd_11[k]
                  + f_3 * pc_y[k] * ppd_29[k];

        t_49[k] = f_1 * ppp0_14[k]
                  - f_2 * ppp1_14[k]
                  + f_3 * pc_z[k] * ppd_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_z, pc_x, pc_y, pc_z, spf0_20, \
                         spf1_20, psf0_11, psd_6, psf1_11, ppd_30, \
                         ppd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * spf0_20[k]
                  - f_4 * pc_y[k] * spf1_20[k];

        t_51[k] = pb_z[k] * psf0_11[k]
                  - f_4 * pc_z[k] * psf1_11[k];

        t_52[k] = f_0 * psd_6[k]
                  + f_3 * pc_z[k] * ppd_30[k];

        t_53[k] = f_3 * pc_x[k] * ppd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pb_z, pc_x, pc_y, pc_z, spd_17, \
                         psf0_16, psd_9, psf1_16, ppd_33, ppd_34, \
                         ppd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * pc_x[k] * ppd_34[k];

        t_55[k] = f_3 * pc_x[k] * ppd_35[k];

        t_56[k] = pb_z[k] * psf0_16[k]
                  - f_4 * pc_z[k] * psf1_16[k];

        t_57[k] = f_0 * psd_9[k]
                  + f_3 * pc_z[k] * ppd_33[k];

        t_58[k] = f_0 * spd_17[k]
                  + f_3 * pc_y[k] * ppd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pa_z, pc_y, pc_z, spf0_0, spf0_2, \
                         spf0_29, spd_0, spf1_0, spf1_2, spf1_29, \
                         ppd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * spf0_29[k]
                  - f_4 * pc_y[k] * spf1_29[k];

        t_60[k] = pa_z[k] * spf0_0[k]
                  - f_4 * pc_z[k] * spf1_0[k];

        t_61[k] = f_3 * pc_y[k] * ppd_36[k];

        t_62[k] = pa_z[k] * spf0_2[k]
                  + f_0 * spd_0[k]
                  - f_4 * pc_z[k] * spf1_2[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_x, pc_z, spf0_6, spf1_6, psd_15, \
                         psd_16, psd_17, ppd_39, ppd_40, ppd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * psd_15[k]
                  + f_3 * pc_x[k] * ppd_39[k];

        t_64[k] = f_0 * psd_16[k]
                  + f_3 * pc_x[k] * ppd_40[k];

        t_65[k] = f_0 * psd_17[k]
                  + f_3 * pc_x[k] * ppd_41[k];

        t_66[k] = pa_z[k] * spf0_6[k]
                  - f_4 * pc_z[k] * spf1_6[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, pc_x, pc_y, pc_z, spf0_10, \
                         spf1_10, psf0_27, psf0_29, psf1_27, psf1_29, \
                         ppd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * psf0_27[k]
                  - f_4 * pc_x[k] * psf1_27[k];

        t_68[k] = f_3 * pc_y[k] * ppd_41[k];

        t_69[k] = pb_x[k] * psf0_29[k]
                  - f_4 * pc_x[k] * psf1_29[k];

        t_70[k] = pa_z[k] * spf0_10[k]
                  - f_4 * pc_z[k] * spf1_10[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pb_y, pc_x, pc_y, psf0_22, psd_12, \
                         psf1_22, ppd_42, ppd_45, ppd_46, ppd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * psd_12[k]
                  + f_3 * pc_y[k] * ppd_42[k];

        t_72[k] = pb_y[k] * psf0_22[k]
                  - f_4 * pc_y[k] * psf1_22[k];

        t_73[k] = f_3 * pc_x[k] * ppd_45[k];

        t_74[k] = f_3 * pc_x[k] * ppd_46[k];

        t_75[k] = f_3 * pc_x[k] * ppd_47[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pc_y, pc_z, spf0_16, spf0_17, spd_9, spf1_16, \
                         spf1_17, psd_17, ppd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_z[k] * spf0_16[k]
                  - f_4 * pc_z[k] * spf1_16[k];

        t_77[k] = pa_z[k] * spf0_17[k]
                  + f_0 * spd_9[k]
                  - f_4 * pc_z[k] * spf1_17[k];

        t_78[k] = f_0 * psd_17[k]
                  + f_3 * pc_y[k] * ppd_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_x, pc_y, psf0_29, psf1_29, ppp0_24, \
                         ppp0_26, ppp1_24, ppp1_26, ppd_48, ppd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * psf0_29[k]
                  - f_4 * pc_y[k] * psf1_29[k];

        t_80[k] = f_1 * ppp0_24[k]
                  - f_2 * ppp1_24[k]
                  + f_3 * pc_x[k] * ppd_48[k];

        t_81[k] = f_3 * pc_y[k] * ppd_48[k];

        t_82[k] = f_5 * ppp0_26[k]
                  - f_6 * ppp1_26[k]
                  + f_3 * pc_x[k] * ppd_50[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, pc_x, pc_y, ppp0_25, ppp0_26, \
                         ppp1_25, ppp1_26, ppd_51, ppd_52, ppd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * pc_x[k] * ppd_51[k];

        t_84[k] = f_3 * pc_x[k] * ppd_52[k];

        t_85[k] = f_3 * pc_x[k] * ppd_53[k];

        t_86[k] = f_1 * ppp0_25[k]
                  - f_2 * ppp1_25[k]
                  + f_3 * pc_y[k] * ppd_51[k];

        t_87[k] = f_5 * ppp0_26[k]
                  - f_6 * ppp1_26[k]
                  + f_3 * pc_y[k] * ppd_52[k];

        t_88[k] = f_3 * pc_y[k] * ppd_53[k];
    }

#pragma omp simd aligned(t_89, pc_z, spd_17, psd_17, ppp0_26, ppp1_26, \
                         ppd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * spd_17[k]
                  + f_0 * psd_17[k]
                  + f_1 * ppp0_26[k]
                  - f_2 * ppp1_26[k]
                  + f_3 * pc_z[k] * ppd_53[k];
    }
}

}  // namespace simdt3ceri
