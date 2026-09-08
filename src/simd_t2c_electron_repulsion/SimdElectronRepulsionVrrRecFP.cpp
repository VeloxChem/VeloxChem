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


#include "SimdElectronRepulsionVrrRecFP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_fp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ds, const size_t dp,
                                     const size_t fs, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 0.5 / p;

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

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_17 = buffer.data(dp + 17);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_y, pa_z, pb_x, pb_y, pb_z, \
                         ds_0, dp_0, fs_0, fs_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds_0[k]
                 + pb_x[k] * fs_0[k];

        t_1[k] = pb_y[k] * fs_0[k];

        t_2[k] = pb_z[k] * fs_0[k];

        t_3[k] = pa_y[k] * dp_0[k];

        t_4[k] = f_1 * ds_0[k]
                 + pb_y[k] * fs_1[k];

        t_5[k] = pb_z[k] * fs_1[k];

        t_6[k] = pa_z[k] * dp_0[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pa_x, pb_x, pb_y, pb_z, ds_0, ds_3, dp_10, \
                         fs_2, fs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * fs_2[k];

        t_8[k] = f_1 * ds_0[k]
                 + pb_z[k] * fs_2[k];

        t_9[k] = f_1 * ds_3[k]
                 + pb_x[k] * fs_3[k];

        t_10[k] = pa_x[k] * dp_10[k];

        t_11[k] = pb_z[k] * fs_3[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_x, pa_y, pb_x, pb_y, ds_5, \
                         dp_6, dp_13, dp_14, dp_17, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * dp_6[k];

        t_13[k] = pa_x[k] * dp_13[k];

        t_14[k] = pa_x[k] * dp_14[k];

        t_15[k] = f_1 * ds_5[k]
                  + pb_x[k] * fs_5[k];

        t_16[k] = pb_y[k] * fs_5[k];

        t_17[k] = pa_x[k] * dp_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, t_24, pa_z, pb_x, pb_y, pb_z, \
                         ds_3, dp_10, fs_6, fs_7, fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_x[k] * fs_6[k];

        t_19[k] = f_0 * ds_3[k]
                  + pb_y[k] * fs_6[k];

        t_20[k] = pb_z[k] * fs_6[k];

        t_21[k] = pb_x[k] * fs_7[k];

        t_22[k] = pa_z[k] * dp_10[k];

        t_23[k] = f_1 * ds_3[k]
                  + pb_z[k] * fs_7[k];

        t_24[k] = pb_x[k] * fs_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_y, pb_z, ds_5, dp_17, \
                         fs_8, fs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ds_5[k]
                  + pb_y[k] * fs_8[k];

        t_26[k] = pa_y[k] * dp_17[k];

        t_27[k] = pb_x[k] * fs_9[k];

        t_28[k] = pb_y[k] * fs_9[k];

        t_29[k] = f_0 * ds_5[k]
                  + pb_z[k] * fs_9[k];
    }
}

}  // namespace simdt2ceri
