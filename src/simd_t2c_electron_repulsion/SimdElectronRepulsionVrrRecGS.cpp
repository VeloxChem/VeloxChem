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


#include "SimdElectronRepulsionVrrRecGS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_gs_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, ds0_0, ds0_1, ds1_0, ds1_1, \
                         fs_0, fs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = pa_y[k] * fs_0[k];

        t_2[k] = pa_z[k] * fs_0[k];

        t_3[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_4, fs_5, fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_4[k];

        t_5[k] = pa_x[k] * fs_5[k];

        t_6[k] = pa_x[k] * fs_8[k];

        t_7[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_5[k];

        t_8[k] = pa_z[k] * fs_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, ds0_2, ds1_2, fs_7, \
                         fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_7[k];

        t_10[k] = pa_y[k] * fs_8[k];

        t_11[k] = f_0 * ds0_2[k]
                  - f_1 * ds1_2[k]
                  + pa_z[k] * fs_8[k];
    }
}

auto
compute_prim_gs_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, \
                         ds1_2, fs_0, fs_1, fs_2, fs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_1[k];

        t_2[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_2[k];

        t_3[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_3[k];
    }

#pragma omp simd aligned(t_4, t_5, pa_y, pa_z, ds0_2, ds1_2, fs_4, \
                         fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_4[k];

        t_5[k] = f_0 * ds0_2[k]
                 - f_1 * ds1_2[k]
                 + pa_z[k] * fs_5[k];
    }
}

auto
compute_prim_gs_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, \
                         ds1_2, fs_0, fs_3, fs_4, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_3[k];

        t_2[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_4[k];

        t_3[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_5[k];
    }

#pragma omp simd aligned(t_4, t_5, pa_y, pa_z, ds0_2, ds1_2, fs_7, \
                         fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_7[k];

        t_5[k] = f_0 * ds0_2[k]
                 - f_1 * ds1_2[k]
                 + pa_z[k] * fs_8[k];
    }
}

auto
compute_prim_gs_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, ds0_0, ds0_1, ds1_0, ds1_1, \
                         fs_0, fs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = pa_y[k] * fs_0[k];

        t_2[k] = pa_z[k] * fs_0[k];

        t_3[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_3, fs_4, fs_6, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_3[k];

        t_5[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_4[k];

        t_6[k] = pa_z[k] * fs_4[k];

        t_7[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_6[k];

        t_8[k] = f_0 * ds0_2[k]
                 - f_1 * ds1_2[k]
                 + pa_z[k] * fs_7[k];
    }
}

auto
compute_prim_gs_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_z, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, \
                         ds1_2, fs_0, fs_1, fs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = pa_z[k] * fs_0[k];

        t_2[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_1[k];

        t_3[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_3, fs_4, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * fs_3[k];

        t_5[k] = pa_x[k] * fs_5[k];

        t_6[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_3[k];

        t_7[k] = pa_z[k] * fs_3[k];

        t_8[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_4[k];

        t_9[k] = pa_y[k] * fs_5[k];
    }

#pragma omp simd aligned(t_10, pa_z, ds0_2, ds1_2, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ds0_2[k]
                  - f_1 * ds1_2[k]
                  + pa_z[k] * fs_5[k];
    }
}

auto
compute_prim_gs_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, ds0_0, ds0_1, ds1_0, ds1_1, \
                         fs_0, fs_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = pa_y[k] * fs_0[k];

        t_2[k] = pa_z[k] * fs_0[k];

        t_3[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_2, fs_3, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_2[k];

        t_5[k] = pa_x[k] * fs_3[k];

        t_6[k] = pa_x[k] * fs_5[k];

        t_7[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_3[k];

        t_8[k] = pa_z[k] * fs_3[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, ds0_2, ds1_2, fs_4, \
                         fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_4[k];

        t_10[k] = pa_y[k] * fs_5[k];

        t_11[k] = f_0 * ds0_2[k]
                  - f_1 * ds1_2[k]
                  + pa_z[k] * fs_5[k];
    }
}

auto
compute_prim_gs_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, ds1_2, \
                         fs_0, fs_3, fs_4, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_3[k];

        t_2[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_4[k];

        t_3[k] = pa_x[k] * fs_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_5, fs_7, fs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * fs_8[k];

        t_5[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_5[k];

        t_6[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_7[k];

        t_7[k] = pa_y[k] * fs_8[k];

        t_8[k] = f_0 * ds0_2[k]
                 - f_1 * ds1_2[k]
                 + pa_z[k] * fs_8[k];
    }
}

auto
compute_prim_gs_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, ds0_0, ds0_1, ds1_0, ds1_1, \
                         fs_0, fs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = pa_y[k] * fs_0[k];

        t_2[k] = pa_z[k] * fs_0[k];

        t_3[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_3, fs_4, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_3[k];

        t_5[k] = pa_x[k] * fs_4[k];

        t_6[k] = pa_x[k] * fs_7[k];

        t_7[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_4[k];

        t_8[k] = pa_z[k] * fs_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, ds0_2, ds1_2, fs_6, \
                         fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_6[k];

        t_10[k] = pa_y[k] * fs_7[k];

        t_11[k] = f_0 * ds0_2[k]
                  - f_1 * ds1_2[k]
                  + pa_z[k] * fs_7[k];
    }
}

auto
compute_prim_gs_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, ds1_2, \
                         fs_0, fs_1, fs_2, fs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_1[k];

        t_2[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_2[k];

        t_3[k] = pa_x[k] * fs_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_3, fs_4, fs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * fs_5[k];

        t_5[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_3[k];

        t_6[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_4[k];

        t_7[k] = pa_y[k] * fs_5[k];

        t_8[k] = f_0 * ds0_2[k]
                 - f_1 * ds1_2[k]
                 + pa_z[k] * fs_5[k];
    }
}

auto
compute_prim_gs_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ds0, const size_t ds1, const size_t fs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, ds1_2, \
                         fs_0, fs_2, fs_3, fs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_2[k];

        t_2[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_3[k];

        t_3[k] = pa_x[k] * fs_4[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_4, fs_6, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * fs_7[k];

        t_5[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_4[k];

        t_6[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_6[k];

        t_7[k] = pa_y[k] * fs_7[k];

        t_8[k] = f_0 * ds0_2[k]
                 - f_1 * ds1_2[k]
                 + pa_z[k] * fs_7[k];
    }
}

auto
compute_prim_gs_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t ds0, const size_t ds1, const size_t fs,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / alpha;
    const auto f_1 = 1.5 * beta / (alpha * p);
    const auto f_2 = 0.5 / alpha;
    const auto f_3 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_z, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, \
                         ds1_2, fs_0, fs_2, fs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ds0_0[k]
                 - f_1 * ds1_0[k]
                 + pa_x[k] * fs_0[k];

        t_1[k] = pa_z[k] * fs_0[k];

        t_2[k] = f_2 * ds0_1[k]
                 - f_3 * ds1_1[k]
                 + pa_x[k] * fs_2[k];

        t_3[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_x[k] * fs_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, ds0_1, ds0_2, ds1_1, \
                         ds1_2, fs_4, fs_6, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * fs_4[k];

        t_5[k] = pa_x[k] * fs_7[k];

        t_6[k] = f_0 * ds0_1[k]
                 - f_1 * ds1_1[k]
                 + pa_y[k] * fs_4[k];

        t_7[k] = pa_z[k] * fs_4[k];

        t_8[k] = f_2 * ds0_2[k]
                 - f_3 * ds1_2[k]
                 + pa_y[k] * fs_6[k];

        t_9[k] = pa_y[k] * fs_7[k];
    }

#pragma omp simd aligned(t_10, pa_z, ds0_2, ds1_2, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ds0_2[k]
                  - f_1 * ds1_2[k]
                  + pa_z[k] * fs_7[k];
    }
}

}  // namespace simdt2ceri
