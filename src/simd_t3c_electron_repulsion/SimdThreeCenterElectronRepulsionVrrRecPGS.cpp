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


#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pgs_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sfs0,
                                                   const size_t sfs1, const size_t sgs0,
                                                   const size_t sgs1, const size_t pfs0,
                                                   const size_t pfs1, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 1.0 / p;
    const auto f_4 = gamma / (p * q);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfs0_0 = buffer.data(sfs0 + 0);
    const auto *sfs0_1 = buffer.data(sfs0 + 1);
    const auto *sfs0_2 = buffer.data(sfs0 + 2);
    const auto *sfs0_3 = buffer.data(sfs0 + 3);
    const auto *sfs0_5 = buffer.data(sfs0 + 5);
    const auto *sfs0_6 = buffer.data(sfs0 + 6);
    const auto *sfs0_7 = buffer.data(sfs0 + 7);
    const auto *sfs0_8 = buffer.data(sfs0 + 8);
    const auto *sfs0_9 = buffer.data(sfs0 + 9);

    const auto *sfs1_0 = buffer.data(sfs1 + 0);
    const auto *sfs1_1 = buffer.data(sfs1 + 1);
    const auto *sfs1_2 = buffer.data(sfs1 + 2);
    const auto *sfs1_3 = buffer.data(sfs1 + 3);
    const auto *sfs1_5 = buffer.data(sfs1 + 5);
    const auto *sfs1_6 = buffer.data(sfs1 + 6);
    const auto *sfs1_7 = buffer.data(sfs1 + 7);
    const auto *sfs1_8 = buffer.data(sfs1 + 8);
    const auto *sfs1_9 = buffer.data(sfs1 + 9);

    const auto *sgs0_0 = buffer.data(sgs0 + 0);
    const auto *sgs0_1 = buffer.data(sgs0 + 1);
    const auto *sgs0_2 = buffer.data(sgs0 + 2);
    const auto *sgs0_3 = buffer.data(sgs0 + 3);
    const auto *sgs0_4 = buffer.data(sgs0 + 4);
    const auto *sgs0_5 = buffer.data(sgs0 + 5);
    const auto *sgs0_6 = buffer.data(sgs0 + 6);
    const auto *sgs0_7 = buffer.data(sgs0 + 7);
    const auto *sgs0_8 = buffer.data(sgs0 + 8);
    const auto *sgs0_9 = buffer.data(sgs0 + 9);
    const auto *sgs0_10 = buffer.data(sgs0 + 10);
    const auto *sgs0_11 = buffer.data(sgs0 + 11);
    const auto *sgs0_12 = buffer.data(sgs0 + 12);
    const auto *sgs0_13 = buffer.data(sgs0 + 13);
    const auto *sgs0_14 = buffer.data(sgs0 + 14);

    const auto *sgs1_0 = buffer.data(sgs1 + 0);
    const auto *sgs1_1 = buffer.data(sgs1 + 1);
    const auto *sgs1_2 = buffer.data(sgs1 + 2);
    const auto *sgs1_3 = buffer.data(sgs1 + 3);
    const auto *sgs1_4 = buffer.data(sgs1 + 4);
    const auto *sgs1_5 = buffer.data(sgs1 + 5);
    const auto *sgs1_6 = buffer.data(sgs1 + 6);
    const auto *sgs1_7 = buffer.data(sgs1 + 7);
    const auto *sgs1_8 = buffer.data(sgs1 + 8);
    const auto *sgs1_9 = buffer.data(sgs1 + 9);
    const auto *sgs1_10 = buffer.data(sgs1 + 10);
    const auto *sgs1_11 = buffer.data(sgs1 + 11);
    const auto *sgs1_12 = buffer.data(sgs1 + 12);
    const auto *sgs1_13 = buffer.data(sgs1 + 13);
    const auto *sgs1_14 = buffer.data(sgs1 + 14);

    const auto *pfs0_0 = buffer.data(pfs0 + 0);
    const auto *pfs0_2 = buffer.data(pfs0 + 2);
    const auto *pfs0_16 = buffer.data(pfs0 + 16);
    const auto *pfs0_17 = buffer.data(pfs0 + 17);
    const auto *pfs0_18 = buffer.data(pfs0 + 18);
    const auto *pfs0_22 = buffer.data(pfs0 + 22);
    const auto *pfs0_27 = buffer.data(pfs0 + 27);
    const auto *pfs0_28 = buffer.data(pfs0 + 28);
    const auto *pfs0_29 = buffer.data(pfs0 + 29);

    const auto *pfs1_0 = buffer.data(pfs1 + 0);
    const auto *pfs1_2 = buffer.data(pfs1 + 2);
    const auto *pfs1_16 = buffer.data(pfs1 + 16);
    const auto *pfs1_17 = buffer.data(pfs1 + 17);
    const auto *pfs1_18 = buffer.data(pfs1 + 18);
    const auto *pfs1_22 = buffer.data(pfs1 + 22);
    const auto *pfs1_27 = buffer.data(pfs1 + 27);
    const auto *pfs1_28 = buffer.data(pfs1 + 28);
    const auto *pfs1_29 = buffer.data(pfs1 + 29);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pb_y, pb_z, pc_x, pc_y, pc_z, sfs0_0, sfs1_0, \
                         sgs0_0, sgs1_0, pfs0_0, pfs1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfs0_0[k]
                 - f_1 * sfs1_0[k]
                 + pa_x[k] * sgs0_0[k]
                 - f_2 * pc_x[k] * sgs1_0[k];

        t_1[k] = pb_y[k] * pfs0_0[k]
                 - f_2 * pc_y[k] * pfs1_0[k];

        t_2[k] = pb_z[k] * pfs0_0[k]
                 - f_2 * pc_z[k] * pfs1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pb_y, pc_x, pc_y, sfs0_3, sfs1_3, sgs0_3, sgs1_3, \
                         pfs0_2, pfs1_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sfs0_3[k]
                 - f_4 * sfs1_3[k]
                 + pa_x[k] * sgs0_3[k]
                 - f_2 * pc_x[k] * sgs1_3[k];

        t_4[k] = pb_y[k] * pfs0_2[k]
                 - f_2 * pc_y[k] * pfs1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, sfs0_5, sfs0_6, sfs1_5, sfs1_6, sgs0_5, sgs0_6, \
                         sgs1_5, sgs1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sfs0_5[k]
                 - f_4 * sfs1_5[k]
                 + pa_x[k] * sgs0_5[k]
                 - f_2 * pc_x[k] * sgs1_5[k];

        t_6[k] = f_5 * sfs0_6[k]
                 - f_6 * sfs1_6[k]
                 + pa_x[k] * sgs0_6[k]
                 - f_2 * pc_x[k] * sgs1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_x, pc_x, sfs0_7, sfs0_8, sfs1_7, sfs1_8, sgs0_7, sgs0_8, \
                         sgs1_7, sgs1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * sfs0_7[k]
                 - f_6 * sfs1_7[k]
                 + pa_x[k] * sgs0_7[k]
                 - f_2 * pc_x[k] * sgs1_7[k];

        t_8[k] = f_5 * sfs0_8[k]
                 - f_6 * sfs1_8[k]
                 + pa_x[k] * sgs0_8[k]
                 - f_2 * pc_x[k] * sgs1_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pc_x, sfs0_9, sfs1_9, sgs0_9, sgs0_10, \
                         sgs0_11, sgs0_12, sgs1_9, sgs1_10, sgs1_11, \
                         sgs1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sfs0_9[k]
                 - f_6 * sfs1_9[k]
                 + pa_x[k] * sgs0_9[k]
                 - f_2 * pc_x[k] * sgs1_9[k];

        t_10[k] = pa_x[k] * sgs0_10[k]
                  - f_2 * pc_x[k] * sgs1_10[k];

        t_11[k] = pa_x[k] * sgs0_11[k]
                  - f_2 * pc_x[k] * sgs1_11[k];

        t_12[k] = pa_x[k] * sgs0_12[k]
                  - f_2 * pc_x[k] * sgs1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pc_x, pc_y, sgs0_0, sgs0_13, sgs0_14, \
                         sgs1_0, sgs1_13, sgs1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * sgs0_13[k]
                  - f_2 * pc_x[k] * sgs1_13[k];

        t_14[k] = pa_x[k] * sgs0_14[k]
                  - f_2 * pc_x[k] * sgs1_14[k];

        t_15[k] = pa_y[k] * sgs0_0[k]
                  - f_2 * pc_y[k] * sgs1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pc_y, sfs0_0, sfs0_1, sfs1_0, sfs1_1, sgs0_1, \
                         sgs0_2, sgs0_3, sgs1_1, sgs1_2, sgs1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * sfs0_0[k]
                  - f_6 * sfs1_0[k]
                  + pa_y[k] * sgs0_1[k]
                  - f_2 * pc_y[k] * sgs1_1[k];

        t_17[k] = pa_y[k] * sgs0_2[k]
                  - f_2 * pc_y[k] * sgs1_2[k];

        t_18[k] = f_3 * sfs0_1[k]
                  - f_4 * sfs1_1[k]
                  + pa_y[k] * sgs0_3[k]
                  - f_2 * pc_y[k] * sgs1_3[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pb_x, pc_x, pc_y, sfs0_2, sfs1_2, sgs0_4, \
                         sgs0_5, sgs1_4, sgs1_5, pfs0_16, pfs1_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * sfs0_2[k]
                  - f_6 * sfs1_2[k]
                  + pa_y[k] * sgs0_4[k]
                  - f_2 * pc_y[k] * sgs1_4[k];

        t_20[k] = pa_y[k] * sgs0_5[k]
                  - f_2 * pc_y[k] * sgs1_5[k];

        t_21[k] = pb_x[k] * pfs0_16[k]
                  - f_2 * pc_x[k] * pfs1_16[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_x, pc_x, pc_y, sgs0_9, sgs1_9, pfs0_17, \
                         pfs0_18, pfs1_17, pfs1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_x[k] * pfs0_17[k]
                  - f_2 * pc_x[k] * pfs1_17[k];

        t_23[k] = pb_x[k] * pfs0_18[k]
                  - f_2 * pc_x[k] * pfs1_18[k];

        t_24[k] = pa_y[k] * sgs0_9[k]
                  - f_2 * pc_y[k] * sgs1_9[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_y, pb_z, pc_y, pc_z, sfs0_6, sfs1_6, sgs0_10, sgs1_10, \
                         pfs0_16, pfs1_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sfs0_6[k]
                  - f_1 * sfs1_6[k]
                  + pa_y[k] * sgs0_10[k]
                  - f_2 * pc_y[k] * sgs1_10[k];

        t_26[k] = pb_z[k] * pfs0_16[k]
                  - f_2 * pc_z[k] * pfs1_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pc_y, sfs0_8, sfs0_9, sfs1_8, sfs1_9, \
                         sgs0_12, sgs0_13, sgs0_14, sgs1_12, sgs1_13, \
                         sgs1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * sfs0_8[k]
                  - f_4 * sfs1_8[k]
                  + pa_y[k] * sgs0_12[k]
                  - f_2 * pc_y[k] * sgs1_12[k];

        t_28[k] = f_5 * sfs0_9[k]
                  - f_6 * sfs1_9[k]
                  + pa_y[k] * sgs0_13[k]
                  - f_2 * pc_y[k] * sgs1_13[k];

        t_29[k] = pa_y[k] * sgs0_14[k]
                  - f_2 * pc_y[k] * sgs1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pc_z, sfs0_0, sfs1_0, sgs0_0, sgs0_1, \
                         sgs0_2, sgs0_3, sgs1_0, sgs1_1, sgs1_2, \
                         sgs1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * sgs0_0[k]
                  - f_2 * pc_z[k] * sgs1_0[k];

        t_31[k] = pa_z[k] * sgs0_1[k]
                  - f_2 * pc_z[k] * sgs1_1[k];

        t_32[k] = f_5 * sfs0_0[k]
                  - f_6 * sfs1_0[k]
                  + pa_z[k] * sgs0_2[k]
                  - f_2 * pc_z[k] * sgs1_2[k];

        t_33[k] = pa_z[k] * sgs0_3[k]
                  - f_2 * pc_z[k] * sgs1_3[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_y, pc_y, pc_z, sfs0_2, sfs1_2, sgs0_5, \
                         sgs0_6, sgs1_5, sgs1_6, pfs0_22, pfs1_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * pfs0_22[k]
                  - f_2 * pc_y[k] * pfs1_22[k];

        t_35[k] = f_3 * sfs0_2[k]
                  - f_4 * sfs1_2[k]
                  + pa_z[k] * sgs0_5[k]
                  - f_2 * pc_z[k] * sgs1_5[k];

        t_36[k] = pa_z[k] * sgs0_6[k]
                  - f_2 * pc_z[k] * sgs1_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pc_x, pc_z, sgs0_10, sgs1_10, \
                         pfs0_27, pfs0_28, pfs0_29, pfs1_27, pfs1_28, \
                         pfs1_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_x[k] * pfs0_27[k]
                  - f_2 * pc_x[k] * pfs1_27[k];

        t_38[k] = pb_x[k] * pfs0_28[k]
                  - f_2 * pc_x[k] * pfs1_28[k];

        t_39[k] = pb_x[k] * pfs0_29[k]
                  - f_2 * pc_x[k] * pfs1_29[k];

        t_40[k] = pa_z[k] * sgs0_10[k]
                  - f_2 * pc_z[k] * sgs1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_z, pc_z, sfs0_6, sfs0_7, sfs1_6, sfs1_7, sgs0_11, \
                         sgs0_12, sgs1_11, sgs1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * sfs0_6[k]
                  - f_6 * sfs1_6[k]
                  + pa_z[k] * sgs0_11[k]
                  - f_2 * pc_z[k] * sgs1_11[k];

        t_42[k] = f_3 * sfs0_7[k]
                  - f_4 * sfs1_7[k]
                  + pa_z[k] * sgs0_12[k]
                  - f_2 * pc_z[k] * sgs1_12[k];
    }

#pragma omp simd aligned(t_43, t_44, pa_z, pb_y, pc_y, pc_z, sfs0_9, sfs1_9, sgs0_14, sgs1_14, \
                         pfs0_29, pfs1_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * pfs0_29[k]
                  - f_2 * pc_y[k] * pfs1_29[k];

        t_44[k] = f_0 * sfs0_9[k]
                  - f_1 * sfs1_9[k]
                  + pa_z[k] * sgs0_14[k]
                  - f_2 * pc_z[k] * sgs1_14[k];
    }
}

}  // namespace simdt3ceri
