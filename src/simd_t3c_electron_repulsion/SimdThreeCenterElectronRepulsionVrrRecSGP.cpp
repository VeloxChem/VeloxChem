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


#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sgp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t sfs,
                                                   const size_t sgs, const size_t ncols,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 1.5 / q;
    const auto f_3 = 0.5 / q;
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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

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

    const auto *sgs_0 = buffer.data(sgs + 0);
    const auto *sgs_1 = buffer.data(sgs + 1);
    const auto *sgs_2 = buffer.data(sgs + 2);
    const auto *sgs_3 = buffer.data(sgs + 3);
    const auto *sgs_4 = buffer.data(sgs + 4);
    const auto *sgs_5 = buffer.data(sgs + 5);
    const auto *sgs_6 = buffer.data(sgs + 6);
    const auto *sgs_7 = buffer.data(sgs + 7);
    const auto *sgs_8 = buffer.data(sgs + 8);
    const auto *sgs_9 = buffer.data(sgs + 9);
    const auto *sgs_10 = buffer.data(sgs + 10);
    const auto *sgs_11 = buffer.data(sgs + 11);
    const auto *sgs_12 = buffer.data(sgs + 12);
    const auto *sgs_13 = buffer.data(sgs + 13);
    const auto *sgs_14 = buffer.data(sgs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, sfs_0, sfs_1, \
                         sfs_2, sgs_0, sgs_1, sgs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfs_0[k]
                 + f_1 * pc_x[k] * sgs_0[k];

        t_1[k] = f_1 * pc_y[k] * sgs_0[k];

        t_2[k] = f_1 * pc_z[k] * sgs_0[k];

        t_3[k] = f_2 * sfs_1[k]
                 + f_1 * pc_x[k] * sgs_1[k];

        t_4[k] = f_3 * sfs_0[k]
                 + f_1 * pc_y[k] * sgs_1[k];

        t_5[k] = f_1 * pc_z[k] * sgs_1[k];

        t_6[k] = f_2 * sfs_2[k]
                 + f_1 * pc_x[k] * sgs_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, sfs_0, sfs_1, \
                         sfs_3, sfs_4, sgs_2, sgs_3, sgs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * sgs_2[k];

        t_8[k] = f_3 * sfs_0[k]
                 + f_1 * pc_z[k] * sgs_2[k];

        t_9[k] = f_4 * sfs_3[k]
                 + f_1 * pc_x[k] * sgs_3[k];

        t_10[k] = f_4 * sfs_1[k]
                  + f_1 * pc_y[k] * sgs_3[k];

        t_11[k] = f_1 * pc_z[k] * sgs_3[k];

        t_12[k] = f_4 * sfs_4[k]
                  + f_1 * pc_x[k] * sgs_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, sfs_1, sfs_2, \
                         sfs_5, sfs_6, sgs_4, sgs_5, sgs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * sfs_2[k]
                  + f_1 * pc_y[k] * sgs_4[k];

        t_14[k] = f_3 * sfs_1[k]
                  + f_1 * pc_z[k] * sgs_4[k];

        t_15[k] = f_4 * sfs_5[k]
                  + f_1 * pc_x[k] * sgs_5[k];

        t_16[k] = f_1 * pc_y[k] * sgs_5[k];

        t_17[k] = f_4 * sfs_2[k]
                  + f_1 * pc_z[k] * sgs_5[k];

        t_18[k] = f_3 * sfs_6[k]
                  + f_1 * pc_x[k] * sgs_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, sfs_3, sfs_4, \
                         sfs_7, sfs_8, sgs_6, sgs_7, sgs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * sfs_3[k]
                  + f_1 * pc_y[k] * sgs_6[k];

        t_20[k] = f_1 * pc_z[k] * sgs_6[k];

        t_21[k] = f_3 * sfs_7[k]
                  + f_1 * pc_x[k] * sgs_7[k];

        t_22[k] = f_4 * sfs_4[k]
                  + f_1 * pc_y[k] * sgs_7[k];

        t_23[k] = f_3 * sfs_3[k]
                  + f_1 * pc_z[k] * sgs_7[k];

        t_24[k] = f_3 * sfs_8[k]
                  + f_1 * pc_x[k] * sgs_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, sfs_4, sfs_5, \
                         sfs_9, sgs_8, sgs_9, sgs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * sfs_5[k]
                  + f_1 * pc_y[k] * sgs_8[k];

        t_26[k] = f_4 * sfs_4[k]
                  + f_1 * pc_z[k] * sgs_8[k];

        t_27[k] = f_3 * sfs_9[k]
                  + f_1 * pc_x[k] * sgs_9[k];

        t_28[k] = f_1 * pc_y[k] * sgs_9[k];

        t_29[k] = f_2 * sfs_5[k]
                  + f_1 * pc_z[k] * sgs_9[k];

        t_30[k] = f_1 * pc_x[k] * sgs_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, sfs_6, \
                         sfs_7, sfs_8, sgs_10, sgs_11, sgs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * sfs_6[k]
                  + f_1 * pc_y[k] * sgs_10[k];

        t_32[k] = f_1 * pc_z[k] * sgs_10[k];

        t_33[k] = f_1 * pc_x[k] * sgs_11[k];

        t_34[k] = f_2 * sfs_7[k]
                  + f_1 * pc_y[k] * sgs_11[k];

        t_35[k] = f_3 * sfs_6[k]
                  + f_1 * pc_z[k] * sgs_11[k];

        t_36[k] = f_1 * pc_x[k] * sgs_12[k];

        t_37[k] = f_4 * sfs_8[k]
                  + f_1 * pc_y[k] * sgs_12[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, pc_x, pc_y, pc_z, sfs_7, \
                         sfs_8, sfs_9, sgs_12, sgs_13, sgs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_4 * sfs_7[k]
                  + f_1 * pc_z[k] * sgs_12[k];

        t_39[k] = f_1 * pc_x[k] * sgs_13[k];

        t_40[k] = f_3 * sfs_9[k]
                  + f_1 * pc_y[k] * sgs_13[k];

        t_41[k] = f_2 * sfs_8[k]
                  + f_1 * pc_z[k] * sgs_13[k];

        t_42[k] = f_1 * pc_x[k] * sgs_14[k];

        t_43[k] = f_1 * pc_y[k] * sgs_14[k];

        t_44[k] = f_0 * sfs_9[k]
                  + f_1 * pc_z[k] * sgs_14[k];
    }
}

}  // namespace simdt3ceri
