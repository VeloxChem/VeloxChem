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


#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_fsp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t dss,
                                                   const size_t fss, const size_t ncols,
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

    const auto *dss_0 = buffer.data(dss + 0);
    const auto *dss_1 = buffer.data(dss + 1);
    const auto *dss_2 = buffer.data(dss + 2);
    const auto *dss_3 = buffer.data(dss + 3);
    const auto *dss_4 = buffer.data(dss + 4);
    const auto *dss_5 = buffer.data(dss + 5);

    const auto *fss_0 = buffer.data(fss + 0);
    const auto *fss_1 = buffer.data(fss + 1);
    const auto *fss_2 = buffer.data(fss + 2);
    const auto *fss_3 = buffer.data(fss + 3);
    const auto *fss_4 = buffer.data(fss + 4);
    const auto *fss_5 = buffer.data(fss + 5);
    const auto *fss_6 = buffer.data(fss + 6);
    const auto *fss_7 = buffer.data(fss + 7);
    const auto *fss_8 = buffer.data(fss + 8);
    const auto *fss_9 = buffer.data(fss + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, dss_0, dss_1, \
                         dss_2, fss_0, fss_1, fss_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dss_0[k]
                 + f_1 * pc_x[k] * fss_0[k];

        t_1[k] = f_1 * pc_y[k] * fss_0[k];

        t_2[k] = f_1 * pc_z[k] * fss_0[k];

        t_3[k] = f_2 * dss_1[k]
                 + f_1 * pc_x[k] * fss_1[k];

        t_4[k] = f_3 * dss_0[k]
                 + f_1 * pc_y[k] * fss_1[k];

        t_5[k] = f_1 * pc_z[k] * fss_1[k];

        t_6[k] = f_2 * dss_2[k]
                 + f_1 * pc_x[k] * fss_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, dss_0, dss_1, \
                         dss_3, dss_4, fss_2, fss_3, fss_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * fss_2[k];

        t_8[k] = f_3 * dss_0[k]
                 + f_1 * pc_z[k] * fss_2[k];

        t_9[k] = f_3 * dss_3[k]
                 + f_1 * pc_x[k] * fss_3[k];

        t_10[k] = f_2 * dss_1[k]
                  + f_1 * pc_y[k] * fss_3[k];

        t_11[k] = f_1 * pc_z[k] * fss_3[k];

        t_12[k] = f_3 * dss_4[k]
                  + f_1 * pc_x[k] * fss_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, dss_1, dss_2, \
                         dss_5, fss_4, fss_5, fss_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * dss_2[k]
                  + f_1 * pc_y[k] * fss_4[k];

        t_14[k] = f_3 * dss_1[k]
                  + f_1 * pc_z[k] * fss_4[k];

        t_15[k] = f_3 * dss_5[k]
                  + f_1 * pc_x[k] * fss_5[k];

        t_16[k] = f_1 * pc_y[k] * fss_5[k];

        t_17[k] = f_2 * dss_2[k]
                  + f_1 * pc_z[k] * fss_5[k];

        t_18[k] = f_1 * pc_x[k] * fss_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, pc_x, pc_y, pc_z, dss_3, \
                         dss_4, dss_5, fss_6, fss_7, fss_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dss_3[k]
                  + f_1 * pc_y[k] * fss_6[k];

        t_20[k] = f_1 * pc_z[k] * fss_6[k];

        t_21[k] = f_1 * pc_x[k] * fss_7[k];

        t_22[k] = f_2 * dss_4[k]
                  + f_1 * pc_y[k] * fss_7[k];

        t_23[k] = f_3 * dss_3[k]
                  + f_1 * pc_z[k] * fss_7[k];

        t_24[k] = f_1 * pc_x[k] * fss_8[k];

        t_25[k] = f_3 * dss_5[k]
                  + f_1 * pc_y[k] * fss_8[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, dss_4, dss_5, fss_8, \
                         fss_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * dss_4[k]
                  + f_1 * pc_z[k] * fss_8[k];

        t_27[k] = f_1 * pc_x[k] * fss_9[k];

        t_28[k] = f_1 * pc_y[k] * fss_9[k];

        t_29[k] = f_0 * dss_5[k]
                  + f_1 * pc_z[k] * fss_9[k];
    }
}

}  // namespace simdt3ceri
