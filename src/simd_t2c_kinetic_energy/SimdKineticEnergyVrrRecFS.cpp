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


#include "SimdKineticEnergyVrrRecFS.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_fs_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t ps_s, const size_t ps, const size_t ds,
                                 const size_t fs_s, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * beta / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 2.0 * alpha * beta / p;

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

    const auto *ps_s_0 = buffer.data(ps_s + 0);
    const auto *ps_s_1 = buffer.data(ps_s + 1);
    const auto *ps_s_2 = buffer.data(ps_s + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_6 = buffer.data(fs_s + 6);
    const auto *fs_s_7 = buffer.data(fs_s + 7);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, ps_s_0, ps_0, ds_0, ds_1, \
                         fs_s_0, fs_s_1, fs_s_2, fs_s_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ps_s_0[k]
                 + f_1 * ps_0[k]
                 + pa_x[k] * ds_0[k]
                 + f_2 * fs_s_0[k];

        t_1[k] = pa_y[k] * ds_0[k]
                 + f_2 * fs_s_1[k];

        t_2[k] = pa_z[k] * ds_0[k]
                 + f_2 * fs_s_2[k];

        t_3[k] = pa_x[k] * ds_1[k]
                 + f_2 * fs_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, ps_s_1, ps_1, ds_1, ds_2, \
                         fs_s_4, fs_s_5, fs_s_6, fs_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * ds_2[k]
                 + f_2 * fs_s_4[k];

        t_5[k] = -f_0 * ps_s_1[k]
                 + f_1 * ps_1[k]
                 + pa_y[k] * ds_1[k]
                 + f_2 * fs_s_5[k];

        t_6[k] = pa_z[k] * ds_1[k]
                 + f_2 * fs_s_6[k];

        t_7[k] = pa_y[k] * ds_2[k]
                 + f_2 * fs_s_7[k];
    }

#pragma omp simd aligned(t_8, pa_z, ps_s_2, ps_2, ds_2, fs_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_0 * ps_s_2[k]
                 + f_1 * ps_2[k]
                 + pa_z[k] * ds_2[k]
                 + f_2 * fs_s_8[k];
    }
}

auto
compute_prim_fs_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t ps_s, const size_t ps, const size_t ds,
                                 const size_t fs_s, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * beta / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_s_0 = buffer.data(ps_s + 0);
    const auto *ps_s_1 = buffer.data(ps_s + 1);
    const auto *ps_s_2 = buffer.data(ps_s + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_6 = buffer.data(fs_s + 6);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, ps_s_0, ps_0, ds_0, fs_s_0, fs_s_1, \
                         fs_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ps_s_0[k]
                 + f_1 * ps_0[k]
                 + pa_x[k] * ds_0[k]
                 + f_2 * fs_s_0[k];

        t_1[k] = pa_y[k] * ds_0[k]
                 + f_2 * fs_s_1[k];

        t_2[k] = pa_z[k] * ds_0[k]
                 + f_2 * fs_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pa_z, ps_s_1, ps_s_2, ps_1, ps_2, ds_1, ds_2, \
                         fs_s_3, fs_s_4, fs_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_0 * ps_s_1[k]
                 + f_1 * ps_1[k]
                 + pa_y[k] * ds_1[k]
                 + f_2 * fs_s_3[k];

        t_4[k] = pa_z[k] * ds_1[k]
                 + f_2 * fs_s_4[k];

        t_5[k] = -f_0 * ps_s_2[k]
                 + f_1 * ps_2[k]
                 + pa_z[k] * ds_2[k]
                 + f_2 * fs_s_6[k];
    }
}

auto
compute_prim_fs_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t ps_s, const size_t ps, const size_t ds,
                                 const size_t fs_s, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * beta / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_s_0 = buffer.data(ps_s + 0);
    const auto *ps_s_1 = buffer.data(ps_s + 1);
    const auto *ps_s_2 = buffer.data(ps_s + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, ps_s_0, ps_s_1, ps_0, ps_1, ds_0, \
                         ds_1, fs_s_0, fs_s_1, fs_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ps_s_0[k]
                 + f_1 * ps_0[k]
                 + pa_x[k] * ds_0[k]
                 + f_2 * fs_s_0[k];

        t_1[k] = pa_z[k] * ds_0[k]
                 + f_2 * fs_s_1[k];

        t_2[k] = -f_0 * ps_s_1[k]
                 + f_1 * ps_1[k]
                 + pa_y[k] * ds_1[k]
                 + f_2 * fs_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_z, ps_s_2, ps_2, ds_1, ds_2, fs_s_3, \
                         fs_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_z[k] * ds_1[k]
                 + f_2 * fs_s_3[k];

        t_4[k] = -f_0 * ps_s_2[k]
                 + f_1 * ps_2[k]
                 + pa_z[k] * ds_2[k]
                 + f_2 * fs_s_5[k];
    }
}

auto
compute_prim_fs_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t ps_s, const size_t ps, const size_t ds,
                                 const size_t fs_s, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * beta / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_s_0 = buffer.data(ps_s + 0);
    const auto *ps_s_1 = buffer.data(ps_s + 1);
    const auto *ps_s_2 = buffer.data(ps_s + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_6 = buffer.data(fs_s + 6);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_z, ps_s_0, ps_0, ds_0, ds_1, ds_2, \
                         fs_s_0, fs_s_1, fs_s_2, fs_s_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ps_s_0[k]
                 + f_1 * ps_0[k]
                 + pa_x[k] * ds_0[k]
                 + f_2 * fs_s_0[k];

        t_1[k] = pa_z[k] * ds_0[k]
                 + f_2 * fs_s_1[k];

        t_2[k] = pa_x[k] * ds_1[k]
                 + f_2 * fs_s_2[k];

        t_3[k] = pa_x[k] * ds_2[k]
                 + f_2 * fs_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, ps_s_1, ps_s_2, ps_1, ps_2, ds_1, \
                         ds_2, fs_s_4, fs_s_5, fs_s_6, fs_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_0 * ps_s_1[k]
                 + f_1 * ps_1[k]
                 + pa_y[k] * ds_1[k]
                 + f_2 * fs_s_4[k];

        t_5[k] = pa_z[k] * ds_1[k]
                 + f_2 * fs_s_5[k];

        t_6[k] = pa_y[k] * ds_2[k]
                 + f_2 * fs_s_6[k];

        t_7[k] = -f_0 * ps_s_2[k]
                 + f_1 * ps_2[k]
                 + pa_z[k] * ds_2[k]
                 + f_2 * fs_s_7[k];
    }
}

auto
compute_prim_fs_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t ps_s, const size_t ps, const size_t ds,
                                 const size_t fs_s, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * beta / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_s_0 = buffer.data(ps_s + 0);
    const auto *ps_s_1 = buffer.data(ps_s + 1);
    const auto *ps_s_2 = buffer.data(ps_s + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, ps_s_0, ps_0, ds_0, ds_1, ds_2, fs_s_0, fs_s_1, \
                         fs_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * ps_s_0[k]
                 + f_1 * ps_0[k]
                 + pa_x[k] * ds_0[k]
                 + f_2 * fs_s_0[k];

        t_1[k] = pa_x[k] * ds_1[k]
                 + f_2 * fs_s_1[k];

        t_2[k] = pa_x[k] * ds_2[k]
                 + f_2 * fs_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pa_z, ps_s_1, ps_s_2, ps_1, ps_2, ds_1, ds_2, \
                         fs_s_3, fs_s_4, fs_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_0 * ps_s_1[k]
                 + f_1 * ps_1[k]
                 + pa_y[k] * ds_1[k]
                 + f_2 * fs_s_3[k];

        t_4[k] = pa_y[k] * ds_2[k]
                 + f_2 * fs_s_4[k];

        t_5[k] = -f_0 * ps_s_2[k]
                 + f_1 * ps_2[k]
                 + pa_z[k] * ds_2[k]
                 + f_2 * fs_s_5[k];
    }
}

}  // namespace simdkin
