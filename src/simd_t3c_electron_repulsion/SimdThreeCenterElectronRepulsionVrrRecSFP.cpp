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


#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sfp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t sds,
                                                   const size_t sfs, const size_t ncols,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = p / q;
    const auto f_2 = 1.0 / q;
    const auto f_3 = 0.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sds_0 = buffer.data(sds + 0);
    const auto *sds_1 = buffer.data(sds + 1);
    const auto *sds_2 = buffer.data(sds + 2);
    const auto *sds_3 = buffer.data(sds + 3);
    const auto *sds_4 = buffer.data(sds + 4);
    const auto *sds_5 = buffer.data(sds + 5);

    const auto *sfs_0 = buffer.data(sfs + 0);
    const auto *sfs_1 = buffer.data(sfs + 1);
    const auto *sfs_2 = buffer.data(sfs + 2);
    const auto *sfs_3 = buffer.data(sfs + 3);
    const auto *sfs_4 = buffer.data(sfs + 4);
    const auto *sfs_5 = buffer.data(sfs + 5);
    const auto *sfs_6 = buffer.data(sfs + 6);
    const auto *sfs_7 = buffer.data(sfs + 7);
    const auto *sfs_8 = buffer.data(sfs + 8);
    const auto *sfs_9 = buffer.data(sfs + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, sds_0, sds_1, \
                         sds_2, sfs_0, sfs_1, sfs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sds_0[k]
                 + f_1 * pc_x[k] * sfs_0[k];

        t_1[k] = f_1 * pc_y[k] * sfs_0[k];

        t_2[k] = f_1 * pc_z[k] * sfs_0[k];

        t_3[k] = f_2 * sds_1[k]
                 + f_1 * pc_x[k] * sfs_1[k];

        t_4[k] = f_3 * sds_0[k]
                 + f_1 * pc_y[k] * sfs_1[k];

        t_5[k] = f_1 * pc_z[k] * sfs_1[k];

        t_6[k] = f_2 * sds_2[k]
                 + f_1 * pc_x[k] * sfs_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, sds_0, sds_1, \
                         sds_3, sds_4, sfs_2, sfs_3, sfs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * sfs_2[k];

        t_8[k] = f_3 * sds_0[k]
                 + f_1 * pc_z[k] * sfs_2[k];

        t_9[k] = f_3 * sds_3[k]
                 + f_1 * pc_x[k] * sfs_3[k];

        t_10[k] = f_2 * sds_1[k]
                  + f_1 * pc_y[k] * sfs_3[k];

        t_11[k] = f_1 * pc_z[k] * sfs_3[k];

        t_12[k] = f_3 * sds_4[k]
                  + f_1 * pc_x[k] * sfs_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, sds_1, sds_2, \
                         sds_5, sfs_4, sfs_5, sfs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * sds_2[k]
                  + f_1 * pc_y[k] * sfs_4[k];

        t_14[k] = f_3 * sds_1[k]
                  + f_1 * pc_z[k] * sfs_4[k];

        t_15[k] = f_3 * sds_5[k]
                  + f_1 * pc_x[k] * sfs_5[k];

        t_16[k] = f_1 * pc_y[k] * sfs_5[k];

        t_17[k] = f_2 * sds_2[k]
                  + f_1 * pc_z[k] * sfs_5[k];

        t_18[k] = f_1 * pc_x[k] * sfs_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, sds_3, \
                         sds_4, sds_5, sfs_6, sfs_7, sfs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * sds_3[k]
                  + f_1 * pc_y[k] * sfs_6[k];

        t_20[k] = f_1 * pc_z[k] * sfs_6[k];

        t_21[k] = f_1 * pc_x[k] * sfs_7[k];

        t_22[k] = f_2 * sds_4[k]
                  + f_1 * pc_y[k] * sfs_7[k];

        t_23[k] = f_3 * sds_3[k]
                  + f_1 * pc_z[k] * sfs_7[k];

        t_24[k] = f_1 * pc_x[k] * sfs_8[k];

        t_25[k] = f_3 * sds_5[k]
                  + f_1 * pc_y[k] * sfs_8[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, sds_4, sds_5, sfs_8, \
                         sfs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * sds_4[k]
                  + f_1 * pc_z[k] * sfs_8[k];

        t_27[k] = f_1 * pc_x[k] * sfs_9[k];

        t_28[k] = f_1 * pc_y[k] * sfs_9[k];

        t_29[k] = f_0 * sds_5[k]
                  + f_1 * pc_z[k] * sfs_9[k];
    }
}

}  // namespace simdt3ceri
