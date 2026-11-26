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


#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_shs_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfs0, const size_t sfs1,
                                                   const size_t sgs0, const size_t sgs1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfs0_0 = buffer.data(sfs0 + 0);
    const auto *sfs0_3 = buffer.data(sfs0 + 3);
    const auto *sfs0_5 = buffer.data(sfs0 + 5);
    const auto *sfs0_6 = buffer.data(sfs0 + 6);
    const auto *sfs0_8 = buffer.data(sfs0 + 8);
    const auto *sfs0_9 = buffer.data(sfs0 + 9);

    const auto *sfs1_0 = buffer.data(sfs1 + 0);
    const auto *sfs1_3 = buffer.data(sfs1 + 3);
    const auto *sfs1_5 = buffer.data(sfs1 + 5);
    const auto *sfs1_6 = buffer.data(sfs1 + 6);
    const auto *sfs1_8 = buffer.data(sfs1 + 8);
    const auto *sfs1_9 = buffer.data(sfs1 + 9);

    const auto *sgs0_0 = buffer.data(sgs0 + 0);
    const auto *sgs0_2 = buffer.data(sgs0 + 2);
    const auto *sgs0_3 = buffer.data(sgs0 + 3);
    const auto *sgs0_5 = buffer.data(sgs0 + 5);
    const auto *sgs0_6 = buffer.data(sgs0 + 6);
    const auto *sgs0_9 = buffer.data(sgs0 + 9);
    const auto *sgs0_10 = buffer.data(sgs0 + 10);
    const auto *sgs0_11 = buffer.data(sgs0 + 11);
    const auto *sgs0_12 = buffer.data(sgs0 + 12);
    const auto *sgs0_13 = buffer.data(sgs0 + 13);
    const auto *sgs0_14 = buffer.data(sgs0 + 14);

    const auto *sgs1_0 = buffer.data(sgs1 + 0);
    const auto *sgs1_2 = buffer.data(sgs1 + 2);
    const auto *sgs1_3 = buffer.data(sgs1 + 3);
    const auto *sgs1_5 = buffer.data(sgs1 + 5);
    const auto *sgs1_6 = buffer.data(sgs1 + 6);
    const auto *sgs1_9 = buffer.data(sgs1 + 9);
    const auto *sgs1_10 = buffer.data(sgs1 + 10);
    const auto *sgs1_11 = buffer.data(sgs1 + 11);
    const auto *sgs1_12 = buffer.data(sgs1 + 12);
    const auto *sgs1_13 = buffer.data(sgs1 + 13);
    const auto *sgs1_14 = buffer.data(sgs1 + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sfs0_0, sfs1_0, \
                         sgs0_0, sgs1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfs0_0[k]
                 - f_1 * sfs1_0[k]
                 + pb_x[k] * sgs0_0[k]
                 - f_2 * pc_x[k] * sgs1_0[k];

        t_1[k] = pb_y[k] * sgs0_0[k]
                 - f_2 * pc_y[k] * sgs1_0[k];

        t_2[k] = pb_z[k] * sgs0_0[k]
                 - f_2 * pc_z[k] * sgs1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sfs0_3, sfs1_3, sgs0_2, sgs0_3, \
                         sgs1_2, sgs1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sfs0_3[k]
                 - f_4 * sfs1_3[k]
                 + pb_x[k] * sgs0_3[k]
                 - f_2 * pc_x[k] * sgs1_3[k];

        t_4[k] = pb_y[k] * sgs0_2[k]
                 - f_2 * pc_y[k] * sgs1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sfs0_5, sfs0_6, sfs1_5, sfs1_6, sgs0_5, sgs0_6, \
                         sgs1_5, sgs1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sfs0_5[k]
                 - f_4 * sfs1_5[k]
                 + pb_x[k] * sgs0_5[k]
                 - f_2 * pc_x[k] * sgs1_5[k];

        t_6[k] = f_5 * sfs0_6[k]
                 - f_6 * sfs1_6[k]
                 + pb_x[k] * sgs0_6[k]
                 - f_2 * pc_x[k] * sgs1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sgs0_3, sgs0_5, sgs1_3, \
                         sgs1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sgs0_3[k]
                 - f_2 * pc_z[k] * sgs1_3[k];

        t_8[k] = pb_y[k] * sgs0_5[k]
                 - f_2 * pc_y[k] * sgs1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, pc_x, sfs0_9, sfs1_9, sgs0_9, sgs0_10, \
                         sgs0_11, sgs0_12, sgs1_9, sgs1_10, sgs1_11, \
                         sgs1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sfs0_9[k]
                 - f_6 * sfs1_9[k]
                 + pb_x[k] * sgs0_9[k]
                 - f_2 * pc_x[k] * sgs1_9[k];

        t_10[k] = pb_x[k] * sgs0_10[k]
                  - f_2 * pc_x[k] * sgs1_10[k];

        t_11[k] = pb_x[k] * sgs0_11[k]
                  - f_2 * pc_x[k] * sgs1_11[k];

        t_12[k] = pb_x[k] * sgs0_12[k]
                  - f_2 * pc_x[k] * sgs1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, pc_x, pc_y, sfs0_6, sfs1_6, sgs0_10, \
                         sgs0_13, sgs0_14, sgs1_10, sgs1_13, sgs1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * sgs0_13[k]
                  - f_2 * pc_x[k] * sgs1_13[k];

        t_14[k] = pb_x[k] * sgs0_14[k]
                  - f_2 * pc_x[k] * sgs1_14[k];

        t_15[k] = f_0 * sfs0_6[k]
                  - f_1 * sfs1_6[k]
                  + pb_y[k] * sgs0_10[k]
                  - f_2 * pc_y[k] * sgs1_10[k];
    }

#pragma omp simd aligned(t_16, t_17, pb_y, pb_z, pc_y, pc_z, sfs0_8, sfs1_8, sgs0_10, sgs0_12, \
                         sgs1_10, sgs1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * sgs0_10[k]
                  - f_2 * pc_z[k] * sgs1_10[k];

        t_17[k] = f_3 * sfs0_8[k]
                  - f_4 * sfs1_8[k]
                  + pb_y[k] * sgs0_12[k]
                  - f_2 * pc_y[k] * sgs1_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, pb_z, pc_y, pc_z, sfs0_9, sfs1_9, sgs0_13, \
                         sgs0_14, sgs1_13, sgs1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sfs0_9[k]
                  - f_6 * sfs1_9[k]
                  + pb_y[k] * sgs0_13[k]
                  - f_2 * pc_y[k] * sgs1_13[k];

        t_19[k] = pb_y[k] * sgs0_14[k]
                  - f_2 * pc_y[k] * sgs1_14[k];

        t_20[k] = f_0 * sfs0_9[k]
                  - f_1 * sfs1_9[k]
                  + pb_z[k] * sgs0_14[k]
                  - f_2 * pc_z[k] * sgs1_14[k];
    }
}

}  // namespace simdt3ceri
