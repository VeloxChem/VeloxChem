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


#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_psh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ssh0, const size_t ssg,
                                                   const size_t ssh1, const size_t psf0,
                                                   const size_t psf1, const size_t psg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / gamma;
    const auto f_4 = 0.5 * p / (gamma * q);
    const auto f_5 = 1.0 / gamma;
    const auto f_6 = p / (gamma * q);
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssh0_0 = buffer.data(ssh0 + 0);
    const auto *ssh0_3 = buffer.data(ssh0 + 3);
    const auto *ssh0_5 = buffer.data(ssh0 + 5);
    const auto *ssh0_6 = buffer.data(ssh0 + 6);
    const auto *ssh0_9 = buffer.data(ssh0 + 9);
    const auto *ssh0_15 = buffer.data(ssh0 + 15);
    const auto *ssh0_17 = buffer.data(ssh0 + 17);
    const auto *ssh0_18 = buffer.data(ssh0 + 18);
    const auto *ssh0_20 = buffer.data(ssh0 + 20);

    const auto *ssg_0 = buffer.data(ssg + 0);
    const auto *ssg_10 = buffer.data(ssg + 10);
    const auto *ssg_12 = buffer.data(ssg + 12);
    const auto *ssg_14 = buffer.data(ssg + 14);

    const auto *ssh1_0 = buffer.data(ssh1 + 0);
    const auto *ssh1_3 = buffer.data(ssh1 + 3);
    const auto *ssh1_5 = buffer.data(ssh1 + 5);
    const auto *ssh1_6 = buffer.data(ssh1 + 6);
    const auto *ssh1_9 = buffer.data(ssh1 + 9);
    const auto *ssh1_15 = buffer.data(ssh1 + 15);
    const auto *ssh1_17 = buffer.data(ssh1 + 17);
    const auto *ssh1_18 = buffer.data(ssh1 + 18);
    const auto *ssh1_20 = buffer.data(ssh1 + 20);

    const auto *psf0_0 = buffer.data(psf0 + 0);
    const auto *psf0_1 = buffer.data(psf0 + 1);
    const auto *psf0_2 = buffer.data(psf0 + 2);
    const auto *psf0_11 = buffer.data(psf0 + 11);
    const auto *psf0_13 = buffer.data(psf0 + 13);
    const auto *psf0_16 = buffer.data(psf0 + 16);
    const auto *psf0_17 = buffer.data(psf0 + 17);
    const auto *psf0_18 = buffer.data(psf0 + 18);
    const auto *psf0_22 = buffer.data(psf0 + 22);
    const auto *psf0_25 = buffer.data(psf0 + 25);
    const auto *psf0_27 = buffer.data(psf0 + 27);
    const auto *psf0_28 = buffer.data(psf0 + 28);
    const auto *psf0_29 = buffer.data(psf0 + 29);

    const auto *psf1_0 = buffer.data(psf1 + 0);
    const auto *psf1_1 = buffer.data(psf1 + 1);
    const auto *psf1_2 = buffer.data(psf1 + 2);
    const auto *psf1_11 = buffer.data(psf1 + 11);
    const auto *psf1_13 = buffer.data(psf1 + 13);
    const auto *psf1_16 = buffer.data(psf1 + 16);
    const auto *psf1_17 = buffer.data(psf1 + 17);
    const auto *psf1_18 = buffer.data(psf1 + 18);
    const auto *psf1_22 = buffer.data(psf1 + 22);
    const auto *psf1_25 = buffer.data(psf1 + 25);
    const auto *psf1_27 = buffer.data(psf1 + 27);
    const auto *psf1_28 = buffer.data(psf1 + 28);
    const auto *psf1_29 = buffer.data(psf1 + 29);

    const auto *psg_0 = buffer.data(psg + 0);
    const auto *psg_1 = buffer.data(psg + 1);
    const auto *psg_2 = buffer.data(psg + 2);
    const auto *psg_3 = buffer.data(psg + 3);
    const auto *psg_5 = buffer.data(psg + 5);
    const auto *psg_6 = buffer.data(psg + 6);
    const auto *psg_9 = buffer.data(psg + 9);
    const auto *psg_10 = buffer.data(psg + 10);
    const auto *psg_12 = buffer.data(psg + 12);
    const auto *psg_14 = buffer.data(psg + 14);
    const auto *psg_15 = buffer.data(psg + 15);
    const auto *psg_16 = buffer.data(psg + 16);
    const auto *psg_18 = buffer.data(psg + 18);
    const auto *psg_21 = buffer.data(psg + 21);
    const auto *psg_23 = buffer.data(psg + 23);
    const auto *psg_25 = buffer.data(psg + 25);
    const auto *psg_26 = buffer.data(psg + 26);
    const auto *psg_27 = buffer.data(psg + 27);
    const auto *psg_28 = buffer.data(psg + 28);
    const auto *psg_29 = buffer.data(psg + 29);
    const auto *psg_30 = buffer.data(psg + 30);
    const auto *psg_32 = buffer.data(psg + 32);
    const auto *psg_35 = buffer.data(psg + 35);
    const auto *psg_37 = buffer.data(psg + 37);
    const auto *psg_39 = buffer.data(psg + 39);
    const auto *psg_40 = buffer.data(psg + 40);
    const auto *psg_41 = buffer.data(psg + 41);
    const auto *psg_42 = buffer.data(psg + 42);
    const auto *psg_43 = buffer.data(psg + 43);
    const auto *psg_44 = buffer.data(psg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pc_x, pc_y, pc_z, ssh0_0, ssg_0, ssh1_0, \
                         psf0_0, psf1_0, psg_0, psg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pa_x[k] * ssh0_0[k]
                 + f_0 * ssg_0[k]
                 - f_1 * pc_x[k] * ssh1_0[k];

        t_1[k] = f_2 * pc_y[k] * psg_0[k];

        t_2[k] = f_2 * pc_z[k] * psg_0[k];

        t_3[k] = f_3 * psf0_0[k]
                 - f_4 * psf1_0[k]
                 + f_2 * pc_y[k] * psg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_y, pc_z, psf0_0, psf0_1, psf1_0, psf1_1, \
                         psg_2, psg_3, psg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * psg_2[k];

        t_5[k] = f_3 * psf0_0[k]
                 - f_4 * psf1_0[k]
                 + f_2 * pc_z[k] * psg_2[k];

        t_6[k] = f_5 * psf0_1[k]
                 - f_6 * psf1_1[k]
                 + f_2 * pc_y[k] * psg_3[k];

        t_7[k] = f_2 * pc_z[k] * psg_3[k];

        t_8[k] = f_2 * pc_y[k] * psg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_x, pc_z, ssg_10, ssg_12, psf0_2, psf1_2, \
                         psg_5, psg_6, psg_10, psg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * psf0_2[k]
                 - f_6 * psf1_2[k]
                 + f_2 * pc_z[k] * psg_5[k];

        t_10[k] = f_7 * ssg_10[k]
                  + f_2 * pc_x[k] * psg_10[k];

        t_11[k] = f_2 * pc_z[k] * psg_6[k];

        t_12[k] = f_7 * ssg_12[k]
                  + f_2 * pc_x[k] * psg_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pc_x, pc_y, pc_z, ssh0_15, ssg_14, \
                         ssh1_15, psg_9, psg_10, psg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * psg_9[k];

        t_14[k] = f_7 * ssg_14[k]
                  + f_2 * pc_x[k] * psg_14[k];

        t_15[k] = pa_x[k] * ssh0_15[k]
                  - f_1 * pc_x[k] * ssh1_15[k];

        t_16[k] = f_2 * pc_z[k] * psg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pc_x, pc_y, ssh0_17, ssh0_18, ssh0_20, \
                         ssh1_17, ssh1_18, ssh1_20, psg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_x[k] * ssh0_17[k]
                  - f_1 * pc_x[k] * ssh1_17[k];

        t_18[k] = pa_x[k] * ssh0_18[k]
                  - f_1 * pc_x[k] * ssh1_18[k];

        t_19[k] = f_2 * pc_y[k] * psg_14[k];

        t_20[k] = pa_x[k] * ssh0_20[k]
                  - f_1 * pc_x[k] * ssh1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pc_x, pc_y, pc_z, ssh0_0, ssh1_0, psf0_11, \
                         psf1_11, psg_15, psg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_y[k] * ssh0_0[k]
                  - f_1 * pc_y[k] * ssh1_0[k];

        t_22[k] = f_8 * psf0_11[k]
                  - f_9 * psf1_11[k]
                  + f_2 * pc_x[k] * psg_16[k];

        t_23[k] = f_2 * pc_z[k] * psg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pc_x, pc_y, pc_z, ssh0_5, ssh1_5, psf0_13, \
                         psf1_13, psg_16, psg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * psf0_13[k]
                  - f_6 * psf1_13[k]
                  + f_2 * pc_x[k] * psg_18[k];

        t_25[k] = f_2 * pc_z[k] * psg_16[k];

        t_26[k] = pa_y[k] * ssh0_5[k]
                  - f_1 * pc_y[k] * ssh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pc_x, pc_z, psf0_16, psf0_18, psf1_16, psf1_18, \
                         psg_18, psg_21, psg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * psf0_16[k]
                  - f_4 * psf1_16[k]
                  + f_2 * pc_x[k] * psg_21[k];

        t_28[k] = f_2 * pc_z[k] * psg_18[k];

        t_29[k] = f_3 * psf0_18[k]
                  - f_4 * psf1_18[k]
                  + f_2 * pc_x[k] * psg_23[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_y, pc_x, pc_y, ssh0_9, ssh1_9, \
                         psg_25, psg_26, psg_27, psg_28, psg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * ssh0_9[k]
                  - f_1 * pc_y[k] * ssh1_9[k];

        t_31[k] = f_2 * pc_x[k] * psg_25[k];

        t_32[k] = f_2 * pc_x[k] * psg_26[k];

        t_33[k] = f_2 * pc_x[k] * psg_27[k];

        t_34[k] = f_2 * pc_x[k] * psg_28[k];

        t_35[k] = f_2 * pc_x[k] * psg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pc_y, pc_z, ssh0_15, ssg_10, ssh1_15, \
                         psf0_16, psf1_16, psg_25, psg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_y[k] * ssh0_15[k]
                  + f_0 * ssg_10[k]
                  - f_1 * pc_y[k] * ssh1_15[k];

        t_37[k] = f_2 * pc_z[k] * psg_25[k];

        t_38[k] = f_3 * psf0_16[k]
                  - f_4 * psf1_16[k]
                  + f_2 * pc_z[k] * psg_26[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pc_y, pc_z, ssh0_20, ssg_14, ssh1_20, \
                         psf0_17, psf1_17, psg_27, psg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * psf0_17[k]
                  - f_6 * psf1_17[k]
                  + f_2 * pc_z[k] * psg_27[k];

        t_40[k] = f_7 * ssg_14[k]
                  + f_2 * pc_y[k] * psg_29[k];

        t_41[k] = pa_y[k] * ssh0_20[k]
                  - f_1 * pc_y[k] * ssh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pc_x, pc_y, pc_z, ssh0_0, ssh0_3, \
                         ssh1_0, ssh1_3, psf0_22, psf1_22, psg_30, \
                         psg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * ssh0_0[k]
                  - f_1 * pc_z[k] * ssh1_0[k];

        t_43[k] = f_2 * pc_y[k] * psg_30[k];

        t_44[k] = f_8 * psf0_22[k]
                  - f_9 * psf1_22[k]
                  + f_2 * pc_x[k] * psg_32[k];

        t_45[k] = pa_z[k] * ssh0_3[k]
                  - f_1 * pc_z[k] * ssh1_3[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_z, pc_x, pc_y, pc_z, ssh0_6, ssh1_6, psf0_25, \
                         psf1_25, psg_32, psg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * pc_y[k] * psg_32[k];

        t_47[k] = f_5 * psf0_25[k]
                  - f_6 * psf1_25[k]
                  + f_2 * pc_x[k] * psg_35[k];

        t_48[k] = pa_z[k] * ssh0_6[k]
                  - f_1 * pc_z[k] * ssh1_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, psf0_27, psf0_29, psf1_27, \
                         psf1_29, psg_35, psg_37, psg_39, psg_40, \
                         psg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * psf0_27[k]
                  - f_4 * psf1_27[k]
                  + f_2 * pc_x[k] * psg_37[k];

        t_50[k] = f_2 * pc_y[k] * psg_35[k];

        t_51[k] = f_3 * psf0_29[k]
                  - f_4 * psf1_29[k]
                  + f_2 * pc_x[k] * psg_39[k];

        t_52[k] = f_2 * pc_x[k] * psg_40[k];

        t_53[k] = f_2 * pc_x[k] * psg_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pc_x, pc_z, ssh0_15, ssh1_15, psg_42, \
                         psg_43, psg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_2 * pc_x[k] * psg_42[k];

        t_55[k] = f_2 * pc_x[k] * psg_43[k];

        t_56[k] = f_2 * pc_x[k] * psg_44[k];

        t_57[k] = pa_z[k] * ssh0_15[k]
                  - f_1 * pc_z[k] * ssh1_15[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pc_y, psf0_27, psf0_28, psf0_29, psf1_27, \
                         psf1_28, psf1_29, psg_41, psg_42, psg_43, \
                         psg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_8 * psf0_27[k]
                  - f_9 * psf1_27[k]
                  + f_2 * pc_y[k] * psg_41[k];

        t_59[k] = f_5 * psf0_28[k]
                  - f_6 * psf1_28[k]
                  + f_2 * pc_y[k] * psg_42[k];

        t_60[k] = f_3 * psf0_29[k]
                  - f_4 * psf1_29[k]
                  + f_2 * pc_y[k] * psg_43[k];

        t_61[k] = f_2 * pc_y[k] * psg_44[k];
    }

#pragma omp simd aligned(t_62, pa_z, pc_z, ssh0_20, ssg_14, ssh1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * ssh0_20[k]
                  + f_0 * ssg_14[k]
                  - f_1 * pc_z[k] * ssh1_20[k];
    }
}

}  // namespace simdt3ceri
