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


#include "SimdOverlapVrrRecHS.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_hs_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t fs, const size_t gs, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pa_z, fs_0, fs_3, fs_5, \
                         gs_0, gs_2, gs_3, gs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pa_x[k] * gs_0[k];

        t_1[k] = pa_y[k] * gs_0[k];

        t_2[k] = pa_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_3[k]
                 + pa_x[k] * gs_3[k];

        t_4[k] = pa_y[k] * gs_2[k];

        t_5[k] = f_1 * fs_5[k]
                 + pa_x[k] * gs_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, fs_6, fs_9, gs_3, gs_5, \
                         gs_6, gs_9, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * fs_6[k]
                 + pa_x[k] * gs_6[k];

        t_7[k] = pa_z[k] * gs_3[k];

        t_8[k] = pa_y[k] * gs_5[k];

        t_9[k] = f_2 * fs_9[k]
                 + pa_x[k] * gs_9[k];

        t_10[k] = pa_x[k] * gs_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, fs_6, gs_10, \
                         gs_11, gs_12, gs_13, gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_x[k] * gs_11[k];

        t_12[k] = pa_x[k] * gs_12[k];

        t_13[k] = pa_x[k] * gs_13[k];

        t_14[k] = pa_x[k] * gs_14[k];

        t_15[k] = f_0 * fs_6[k]
                  + pa_y[k] * gs_10[k];

        t_16[k] = pa_z[k] * gs_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, fs_8, fs_9, gs_12, gs_13, \
                         gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * fs_8[k]
                  + pa_y[k] * gs_12[k];

        t_18[k] = f_2 * fs_9[k]
                  + pa_y[k] * gs_13[k];

        t_19[k] = pa_y[k] * gs_14[k];

        t_20[k] = f_0 * fs_9[k]
                  + pa_z[k] * gs_14[k];
    }
}

auto
compute_prim_hs_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t fs, const size_t gs, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pa_z, fs_0, fs_2, fs_3, \
                         gs_0, gs_2, gs_3, gs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pa_x[k] * gs_0[k];

        t_1[k] = pa_y[k] * gs_0[k];

        t_2[k] = pa_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_2[k]
                 + pa_x[k] * gs_3[k];

        t_4[k] = pa_y[k] * gs_2[k];

        t_5[k] = f_1 * fs_3[k]
                 + pa_x[k] * gs_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, fs_4, fs_7, gs_3, gs_5, \
                         gs_6, gs_9, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * fs_4[k]
                 + pa_x[k] * gs_6[k];

        t_7[k] = pa_z[k] * gs_3[k];

        t_8[k] = pa_y[k] * gs_5[k];

        t_9[k] = f_2 * fs_7[k]
                 + pa_x[k] * gs_9[k];

        t_10[k] = pa_x[k] * gs_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, fs_4, gs_10, \
                         gs_11, gs_12, gs_13, gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_x[k] * gs_11[k];

        t_12[k] = pa_x[k] * gs_12[k];

        t_13[k] = pa_x[k] * gs_13[k];

        t_14[k] = pa_x[k] * gs_14[k];

        t_15[k] = f_0 * fs_4[k]
                  + pa_y[k] * gs_10[k];

        t_16[k] = pa_z[k] * gs_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, fs_6, fs_7, gs_12, gs_13, \
                         gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * fs_6[k]
                  + pa_y[k] * gs_12[k];

        t_18[k] = f_2 * fs_7[k]
                  + pa_y[k] * gs_13[k];

        t_19[k] = pa_y[k] * gs_14[k];

        t_20[k] = f_0 * fs_7[k]
                  + pa_z[k] * gs_14[k];
    }
}

auto
compute_prim_hs_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t fs, const size_t gs, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pa_z, fs_0, fs_1, fs_2, \
                         gs_0, gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pa_x[k] * gs_0[k];

        t_1[k] = pa_y[k] * gs_0[k];

        t_2[k] = pa_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_1[k]
                 + pa_x[k] * gs_2[k];

        t_4[k] = pa_y[k] * gs_1[k];

        t_5[k] = f_1 * fs_2[k]
                 + pa_x[k] * gs_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, fs_3, fs_5, gs_2, gs_3, \
                         gs_4, gs_5, gs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * fs_3[k]
                 + pa_x[k] * gs_4[k];

        t_7[k] = pa_z[k] * gs_2[k];

        t_8[k] = pa_y[k] * gs_3[k];

        t_9[k] = f_2 * fs_5[k]
                 + pa_x[k] * gs_5[k];

        t_10[k] = pa_x[k] * gs_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, fs_3, gs_6, \
                         gs_7, gs_8, gs_9, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_x[k] * gs_7[k];

        t_12[k] = pa_x[k] * gs_8[k];

        t_13[k] = pa_x[k] * gs_9[k];

        t_14[k] = pa_x[k] * gs_10[k];

        t_15[k] = f_0 * fs_3[k]
                  + pa_y[k] * gs_6[k];

        t_16[k] = pa_z[k] * gs_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, fs_4, fs_5, gs_8, gs_9, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * fs_4[k]
                  + pa_y[k] * gs_8[k];

        t_18[k] = f_2 * fs_5[k]
                  + pa_y[k] * gs_9[k];

        t_19[k] = pa_y[k] * gs_10[k];

        t_20[k] = f_0 * fs_5[k]
                  + pa_z[k] * gs_10[k];
    }
}

auto
compute_prim_hs_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t fs, const size_t gs, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_z, fs_0, fs_1, fs_2, fs_3, gs_0, \
                         gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pa_x[k] * gs_0[k];

        t_1[k] = pa_z[k] * gs_0[k];

        t_2[k] = f_1 * fs_1[k]
                 + pa_x[k] * gs_1[k];

        t_3[k] = f_1 * fs_2[k]
                 + pa_x[k] * gs_2[k];

        t_4[k] = f_2 * fs_3[k]
                 + pa_x[k] * gs_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, fs_3, fs_5, gs_4, \
                         gs_5, gs_6, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * fs_5[k]
                 + pa_x[k] * gs_4[k];

        t_6[k] = pa_x[k] * gs_5[k];

        t_7[k] = pa_x[k] * gs_6[k];

        t_8[k] = pa_x[k] * gs_8[k];

        t_9[k] = f_0 * fs_3[k]
                 + pa_y[k] * gs_5[k];

        t_10[k] = pa_z[k] * gs_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_y, pa_z, fs_4, fs_5, gs_6, gs_7, \
                         gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * fs_4[k]
                  + pa_y[k] * gs_6[k];

        t_12[k] = f_2 * fs_5[k]
                  + pa_y[k] * gs_7[k];

        t_13[k] = pa_y[k] * gs_8[k];

        t_14[k] = f_0 * fs_5[k]
                  + pa_z[k] * gs_8[k];
    }
}

}  // namespace simdovl
