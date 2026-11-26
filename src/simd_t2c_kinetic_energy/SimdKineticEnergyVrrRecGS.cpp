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


#include "SimdKineticEnergyVrrRecGS.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_gs_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t ds_s, const size_t ds, const size_t fs,
                                 const size_t gs_s, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * beta / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = beta / p;
    const auto f_4 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds_s_0 = buffer.data(ds_s + 0);
    const auto *ds_s_3 = buffer.data(ds_s + 3);
    const auto *ds_s_5 = buffer.data(ds_s + 5);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_1 = buffer.data(gs_s + 1);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_5 = buffer.data(gs_s + 5);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_9 = buffer.data(gs_s + 9);
    const auto *gs_s_10 = buffer.data(gs_s + 10);
    const auto *gs_s_11 = buffer.data(gs_s + 11);
    const auto *gs_s_12 = buffer.data(gs_s + 12);
    const auto *gs_s_13 = buffer.data(gs_s + 13);
    const auto *gs_s_14 = buffer.data(gs_s + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, ds_s_0, ds_0, fs_0, gs_s_0, gs_s_1, \
                         gs_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ds_s_0[k]
                 + f_1 * ds_0[k]
                 + pa_x[k] * fs_0[k]
                 + f_2 * gs_s_0[k];

        t_1[k] = pa_y[k] * fs_0[k]
                 + f_2 * gs_s_1[k];

        t_2[k] = pa_z[k] * fs_0[k]
                 + f_2 * gs_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, ds_s_3, ds_s_5, ds_3, ds_5, fs_2, fs_3, \
                         fs_5, gs_s_3, gs_s_4, gs_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * ds_s_3[k]
                 + f_4 * ds_3[k]
                 + pa_x[k] * fs_3[k]
                 + f_2 * gs_s_3[k];

        t_4[k] = pa_y[k] * fs_2[k]
                 + f_2 * gs_s_4[k];

        t_5[k] = -f_3 * ds_s_5[k]
                 + f_4 * ds_5[k]
                 + pa_x[k] * fs_5[k]
                 + f_2 * gs_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, fs_6, fs_7, fs_8, fs_9, gs_s_6, gs_s_7, \
                         gs_s_8, gs_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * fs_6[k]
                 + f_2 * gs_s_6[k];

        t_7[k] = pa_x[k] * fs_7[k]
                 + f_2 * gs_s_7[k];

        t_8[k] = pa_x[k] * fs_8[k]
                 + f_2 * gs_s_8[k];

        t_9[k] = pa_x[k] * fs_9[k]
                 + f_2 * gs_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pa_z, ds_s_3, ds_s_5, ds_3, ds_5, fs_6, fs_8, \
                         gs_s_10, gs_s_11, gs_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_0 * ds_s_3[k]
                  + f_1 * ds_3[k]
                  + pa_y[k] * fs_6[k]
                  + f_2 * gs_s_10[k];

        t_11[k] = pa_z[k] * fs_6[k]
                  + f_2 * gs_s_11[k];

        t_12[k] = -f_3 * ds_s_5[k]
                  + f_4 * ds_5[k]
                  + pa_y[k] * fs_8[k]
                  + f_2 * gs_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pa_z, ds_s_5, ds_5, fs_9, gs_s_13, \
                         gs_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fs_9[k]
                  + f_2 * gs_s_13[k];

        t_14[k] = -f_0 * ds_s_5[k]
                  + f_1 * ds_5[k]
                  + pa_z[k] * fs_9[k]
                  + f_2 * gs_s_14[k];
    }
}

}  // namespace simdkin
