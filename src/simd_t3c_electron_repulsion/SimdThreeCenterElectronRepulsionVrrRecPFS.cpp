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


#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pfs_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sds0,
                                                   const size_t sds1, const size_t sfs0,
                                                   const size_t sfs1, const size_t pds0,
                                                   const size_t pds1, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sds0_0 = buffer.data(sds0 + 0);
    const auto *sds0_3 = buffer.data(sds0 + 3);
    const auto *sds0_4 = buffer.data(sds0 + 4);
    const auto *sds0_5 = buffer.data(sds0 + 5);

    const auto *sds1_0 = buffer.data(sds1 + 0);
    const auto *sds1_3 = buffer.data(sds1 + 3);
    const auto *sds1_4 = buffer.data(sds1 + 4);
    const auto *sds1_5 = buffer.data(sds1 + 5);

    const auto *sfs0_0 = buffer.data(sfs0 + 0);
    const auto *sfs0_1 = buffer.data(sfs0 + 1);
    const auto *sfs0_2 = buffer.data(sfs0 + 2);
    const auto *sfs0_3 = buffer.data(sfs0 + 3);
    const auto *sfs0_4 = buffer.data(sfs0 + 4);
    const auto *sfs0_5 = buffer.data(sfs0 + 5);
    const auto *sfs0_6 = buffer.data(sfs0 + 6);
    const auto *sfs0_7 = buffer.data(sfs0 + 7);
    const auto *sfs0_8 = buffer.data(sfs0 + 8);
    const auto *sfs0_9 = buffer.data(sfs0 + 9);

    const auto *sfs1_0 = buffer.data(sfs1 + 0);
    const auto *sfs1_1 = buffer.data(sfs1 + 1);
    const auto *sfs1_2 = buffer.data(sfs1 + 2);
    const auto *sfs1_3 = buffer.data(sfs1 + 3);
    const auto *sfs1_4 = buffer.data(sfs1 + 4);
    const auto *sfs1_5 = buffer.data(sfs1 + 5);
    const auto *sfs1_6 = buffer.data(sfs1 + 6);
    const auto *sfs1_7 = buffer.data(sfs1 + 7);
    const auto *sfs1_8 = buffer.data(sfs1 + 8);
    const auto *sfs1_9 = buffer.data(sfs1 + 9);

    const auto *pds0_0 = buffer.data(pds0 + 0);
    const auto *pds0_9 = buffer.data(pds0 + 9);
    const auto *pds0_10 = buffer.data(pds0 + 10);
    const auto *pds0_16 = buffer.data(pds0 + 16);
    const auto *pds0_17 = buffer.data(pds0 + 17);

    const auto *pds1_0 = buffer.data(pds1 + 0);
    const auto *pds1_9 = buffer.data(pds1 + 9);
    const auto *pds1_10 = buffer.data(pds1 + 10);
    const auto *pds1_16 = buffer.data(pds1 + 16);
    const auto *pds1_17 = buffer.data(pds1 + 17);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pb_y, pb_z, pc_x, pc_y, pc_z, sds0_0, sds1_0, \
                         sfs0_0, sfs1_0, pds0_0, pds1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sds0_0[k]
                 - f_1 * sds1_0[k]
                 + pa_x[k] * sfs0_0[k]
                 - f_2 * pc_x[k] * sfs1_0[k];

        t_1[k] = pb_y[k] * pds0_0[k]
                 - f_2 * pc_y[k] * pds1_0[k];

        t_2[k] = pb_z[k] * pds0_0[k]
                 - f_2 * pc_z[k] * pds1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pc_x, sds0_3, sds0_4, sds1_3, sds1_4, sfs0_3, sfs0_4, \
                         sfs1_3, sfs1_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sds0_3[k]
                 - f_4 * sds1_3[k]
                 + pa_x[k] * sfs0_3[k]
                 - f_2 * pc_x[k] * sfs1_3[k];

        t_4[k] = f_3 * sds0_4[k]
                 - f_4 * sds1_4[k]
                 + pa_x[k] * sfs0_4[k]
                 - f_2 * pc_x[k] * sfs1_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pc_x, sds0_5, sds1_5, sfs0_5, sfs0_6, \
                         sfs0_7, sfs0_8, sfs1_5, sfs1_6, sfs1_7, \
                         sfs1_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sds0_5[k]
                 - f_4 * sds1_5[k]
                 + pa_x[k] * sfs0_5[k]
                 - f_2 * pc_x[k] * sfs1_5[k];

        t_6[k] = pa_x[k] * sfs0_6[k]
                 - f_2 * pc_x[k] * sfs1_6[k];

        t_7[k] = pa_x[k] * sfs0_7[k]
                 - f_2 * pc_x[k] * sfs1_7[k];

        t_8[k] = pa_x[k] * sfs0_8[k]
                 - f_2 * pc_x[k] * sfs1_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pc_x, pc_y, sds0_0, sds1_0, sfs0_0, \
                         sfs0_1, sfs0_9, sfs1_0, sfs1_1, sfs1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * sfs0_9[k]
                 - f_2 * pc_x[k] * sfs1_9[k];

        t_10[k] = pa_y[k] * sfs0_0[k]
                  - f_2 * pc_y[k] * sfs1_0[k];

        t_11[k] = f_3 * sds0_0[k]
                  - f_4 * sds1_0[k]
                  + pa_y[k] * sfs0_1[k]
                  - f_2 * pc_y[k] * sfs1_1[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pb_x, pc_x, pc_y, sfs0_2, sfs0_5, \
                         sfs1_2, sfs1_5, pds0_9, pds0_10, pds1_9, \
                         pds1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * sfs0_2[k]
                  - f_2 * pc_y[k] * sfs1_2[k];

        t_13[k] = pb_x[k] * pds0_9[k]
                  - f_2 * pc_x[k] * pds1_9[k];

        t_14[k] = pb_x[k] * pds0_10[k]
                  - f_2 * pc_x[k] * pds1_10[k];

        t_15[k] = pa_y[k] * sfs0_5[k]
                  - f_2 * pc_y[k] * sfs1_5[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_y, pb_z, pc_y, pc_z, sds0_3, sds1_3, sfs0_6, sfs1_6, \
                         pds0_9, pds1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * sds0_3[k]
                  - f_1 * sds1_3[k]
                  + pa_y[k] * sfs0_6[k]
                  - f_2 * pc_y[k] * sfs1_6[k];

        t_17[k] = pb_z[k] * pds0_9[k]
                  - f_2 * pc_z[k] * pds1_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pc_y, pc_z, sds0_5, sds1_5, sfs0_0, \
                         sfs0_8, sfs0_9, sfs1_0, sfs1_8, sfs1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * sds0_5[k]
                  - f_4 * sds1_5[k]
                  + pa_y[k] * sfs0_8[k]
                  - f_2 * pc_y[k] * sfs1_8[k];

        t_19[k] = pa_y[k] * sfs0_9[k]
                  - f_2 * pc_y[k] * sfs1_9[k];

        t_20[k] = pa_z[k] * sfs0_0[k]
                  - f_2 * pc_z[k] * sfs1_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, pc_z, sds0_0, sds1_0, sfs0_1, sfs0_2, sfs0_3, \
                         sfs1_1, sfs1_2, sfs1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_z[k] * sfs0_1[k]
                  - f_2 * pc_z[k] * sfs1_1[k];

        t_22[k] = f_3 * sds0_0[k]
                  - f_4 * sds1_0[k]
                  + pa_z[k] * sfs0_2[k]
                  - f_2 * pc_z[k] * sfs1_2[k];

        t_23[k] = pa_z[k] * sfs0_3[k]
                  - f_2 * pc_z[k] * sfs1_3[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_z, pb_x, pc_x, pc_z, sfs0_6, sfs1_6, pds0_16, \
                         pds0_17, pds1_16, pds1_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_x[k] * pds0_16[k]
                  - f_2 * pc_x[k] * pds1_16[k];

        t_25[k] = pb_x[k] * pds0_17[k]
                  - f_2 * pc_x[k] * pds1_17[k];

        t_26[k] = pa_z[k] * sfs0_6[k]
                  - f_2 * pc_z[k] * sfs1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_z, pb_y, pc_y, pc_z, sds0_3, sds1_3, sfs0_7, sfs1_7, \
                         pds0_17, pds1_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * sds0_3[k]
                  - f_4 * sds1_3[k]
                  + pa_z[k] * sfs0_7[k]
                  - f_2 * pc_z[k] * sfs1_7[k];

        t_28[k] = pb_y[k] * pds0_17[k]
                  - f_2 * pc_y[k] * pds1_17[k];
    }

#pragma omp simd aligned(t_29, pa_z, pc_z, sds0_5, sds1_5, sfs0_9, \
                         sfs1_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * sds0_5[k]
                  - f_1 * sds1_5[k]
                  + pa_z[k] * sfs0_9[k]
                  - f_2 * pc_z[k] * sfs1_9[k];
    }
}

}  // namespace simdt3ceri
