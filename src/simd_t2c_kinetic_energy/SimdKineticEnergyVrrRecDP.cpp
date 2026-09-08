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


#include "SimdKineticEnergyVrrRecDP.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_dp_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ps, const size_t pp,
                                 const size_t dp_s, const size_t ds, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_8 = buffer.data(pp + 8);

    const auto *dp_s_0 = buffer.data(dp_s + 0);
    const auto *dp_s_1 = buffer.data(dp_s + 1);
    const auto *dp_s_2 = buffer.data(dp_s + 2);
    const auto *dp_s_3 = buffer.data(dp_s + 3);
    const auto *dp_s_4 = buffer.data(dp_s + 4);
    const auto *dp_s_5 = buffer.data(dp_s + 5);
    const auto *dp_s_6 = buffer.data(dp_s + 6);
    const auto *dp_s_7 = buffer.data(dp_s + 7);
    const auto *dp_s_8 = buffer.data(dp_s + 8);
    const auto *dp_s_9 = buffer.data(dp_s + 9);
    const auto *dp_s_10 = buffer.data(dp_s + 10);
    const auto *dp_s_11 = buffer.data(dp_s + 11);
    const auto *dp_s_12 = buffer.data(dp_s + 12);
    const auto *dp_s_13 = buffer.data(dp_s + 13);
    const auto *dp_s_14 = buffer.data(dp_s + 14);
    const auto *dp_s_15 = buffer.data(dp_s + 15);
    const auto *dp_s_16 = buffer.data(dp_s + 16);
    const auto *dp_s_17 = buffer.data(dp_s + 17);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_4 = buffer.data(ds + 4);
    const auto *ds_5 = buffer.data(ds + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ps_0, pp_0, dp_s_0, \
                         dp_s_1, dp_s_2, dp_s_3, ds_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ps_0[k]
                 + f_1 * dp_s_0[k]
                 + pb_x[k] * ds_0[k];

        t_1[k] = f_1 * dp_s_1[k]
                 + pb_y[k] * ds_0[k];

        t_2[k] = f_1 * dp_s_2[k]
                 + pb_z[k] * ds_0[k];

        t_3[k] = pa_y[k] * pp_0[k]
                 + f_1 * dp_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, pb_y, pb_z, pp_0, pp_4, dp_s_4, \
                         dp_s_5, dp_s_6, dp_s_7, ds_1, ds_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * pp_4[k]
                 + f_1 * dp_s_4[k];

        t_5[k] = f_1 * dp_s_5[k]
                 + pb_z[k] * ds_1[k];

        t_6[k] = pa_z[k] * pp_0[k]
                 + f_1 * dp_s_6[k];

        t_7[k] = f_1 * dp_s_7[k]
                 + pb_y[k] * ds_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pb_x, pb_y, pb_z, ps_1, pp_8, dp_s_8, \
                         dp_s_9, dp_s_10, dp_s_11, ds_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_x[k] * pp_8[k]
                 + f_1 * dp_s_8[k];

        t_9[k] = f_1 * dp_s_9[k]
                 + pb_x[k] * ds_3[k];

        t_10[k] = f_0 * ps_1[k]
                  + f_1 * dp_s_10[k]
                  + pb_y[k] * ds_3[k];

        t_11[k] = f_1 * dp_s_11[k]
                  + pb_z[k] * ds_3[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, pb_x, pp_4, pp_8, dp_s_12, \
                         dp_s_13, dp_s_14, dp_s_15, ds_4, ds_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * dp_s_12[k]
                  + pb_x[k] * ds_4[k];

        t_13[k] = pa_z[k] * pp_4[k]
                  + f_1 * dp_s_13[k];

        t_14[k] = pa_y[k] * pp_8[k]
                  + f_1 * dp_s_14[k];

        t_15[k] = f_1 * dp_s_15[k]
                  + pb_x[k] * ds_5[k];
    }

#pragma omp simd aligned(t_16, t_17, pb_y, pb_z, ps_2, dp_s_16, dp_s_17, \
                         ds_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * dp_s_16[k]
                  + pb_y[k] * ds_5[k];

        t_17[k] = f_0 * ps_2[k]
                  + f_1 * dp_s_17[k]
                  + pb_z[k] * ds_5[k];
    }
}

}  // namespace simdkin
