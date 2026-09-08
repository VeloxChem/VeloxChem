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


#include "SimdKineticEnergyVrrRecPD.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_pd_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sp, const size_t sd,
                                 const size_t pd_s, const size_t pp, const size_t ncols,
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

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_3 = buffer.data(pd_s + 3);
    const auto *pd_s_4 = buffer.data(pd_s + 4);
    const auto *pd_s_5 = buffer.data(pd_s + 5);
    const auto *pd_s_6 = buffer.data(pd_s + 6);
    const auto *pd_s_7 = buffer.data(pd_s + 7);
    const auto *pd_s_8 = buffer.data(pd_s + 8);
    const auto *pd_s_9 = buffer.data(pd_s + 9);
    const auto *pd_s_10 = buffer.data(pd_s + 10);
    const auto *pd_s_11 = buffer.data(pd_s + 11);
    const auto *pd_s_12 = buffer.data(pd_s + 12);
    const auto *pd_s_13 = buffer.data(pd_s + 13);
    const auto *pd_s_14 = buffer.data(pd_s + 14);
    const auto *pd_s_15 = buffer.data(pd_s + 15);
    const auto *pd_s_16 = buffer.data(pd_s + 16);
    const auto *pd_s_17 = buffer.data(pd_s + 17);

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_1 = buffer.data(pp + 1);
    const auto *pp_2 = buffer.data(pp + 2);
    const auto *pp_3 = buffer.data(pp + 3);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_5 = buffer.data(pp + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_y, pb_z, sp_0, sd_0, sd_1, pd_s_0, \
                         pd_s_1, pd_s_2, pd_s_3, pp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_1 * pd_s_1[k]
                 + pb_y[k] * pp_0[k];

        t_2[k] = f_1 * pd_s_2[k]
                 + pb_z[k] * pp_0[k];

        t_3[k] = pa_x[k] * sd_1[k]
                 + f_1 * pd_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pb_x, pb_y, sd_0, sd_2, pd_s_4, \
                         pd_s_5, pd_s_6, pd_s_7, pp_1, pp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * pd_s_4[k]
                 + pb_y[k] * pp_1[k];

        t_5[k] = pa_x[k] * sd_2[k]
                 + f_1 * pd_s_5[k];

        t_6[k] = pa_y[k] * sd_0[k]
                 + f_1 * pd_s_6[k];

        t_7[k] = f_1 * pd_s_7[k]
                 + pb_x[k] * pp_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, sp_1, sd_1, sd_2, pd_s_8, \
                         pd_s_9, pd_s_10, pd_s_11, pp_2, pp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * pd_s_8[k]
                 + pb_x[k] * pp_3[k];

        t_9[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_9[k];

        t_10[k] = f_1 * pd_s_10[k]
                  + pb_z[k] * pp_2[k];

        t_11[k] = pa_y[k] * sd_2[k]
                  + f_1 * pd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, sd_0, sd_1, pd_s_12, pd_s_13, \
                         pd_s_14, pd_s_15, pp_4, pp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * sd_0[k]
                  + f_1 * pd_s_12[k];

        t_13[k] = f_1 * pd_s_13[k]
                  + pb_x[k] * pp_4[k];

        t_14[k] = f_1 * pd_s_14[k]
                  + pb_x[k] * pp_5[k];

        t_15[k] = pa_z[k] * sd_1[k]
                  + f_1 * pd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_z, pb_y, sp_2, sd_2, pd_s_16, pd_s_17, \
                         pp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * pd_s_16[k]
                  + pb_y[k] * pp_5[k];

        t_17[k] = f_0 * sp_2[k]
                  + pa_z[k] * sd_2[k]
                  + f_1 * pd_s_17[k];
    }
}

auto
compute_prim_pd_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sp, const size_t sd,
                                 const size_t pd_s, const size_t pp, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_3 = buffer.data(pd_s + 3);
    const auto *pd_s_4 = buffer.data(pd_s + 4);
    const auto *pd_s_6 = buffer.data(pd_s + 6);
    const auto *pd_s_7 = buffer.data(pd_s + 7);
    const auto *pd_s_8 = buffer.data(pd_s + 8);
    const auto *pd_s_9 = buffer.data(pd_s + 9);
    const auto *pd_s_11 = buffer.data(pd_s + 11);
    const auto *pd_s_12 = buffer.data(pd_s + 12);
    const auto *pd_s_13 = buffer.data(pd_s + 13);

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_1 = buffer.data(pp + 1);
    const auto *pp_2 = buffer.data(pp + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_z, sp_0, sd_0, sd_1, sd_2, pd_s_0, \
                         pd_s_1, pd_s_2, pd_s_3, pp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_1 * pd_s_1[k]
                 + pb_z[k] * pp_0[k];

        t_2[k] = pa_x[k] * sd_1[k]
                 + f_1 * pd_s_2[k];

        t_3[k] = pa_x[k] * sd_2[k]
                 + f_1 * pd_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_z, sp_1, sd_0, sd_1, sd_2, pd_s_4, \
                         pd_s_6, pd_s_7, pd_s_8, pp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_y[k] * sd_0[k]
                 + f_1 * pd_s_4[k];

        t_5[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_6[k];

        t_6[k] = f_1 * pd_s_7[k]
                 + pb_z[k] * pp_1[k];

        t_7[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_z, pb_y, sp_2, sd_0, sd_1, sd_2, pd_s_9, \
                         pd_s_11, pd_s_12, pd_s_13, pp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_z[k] * sd_0[k]
                 + f_1 * pd_s_9[k];

        t_9[k] = pa_z[k] * sd_1[k]
                 + f_1 * pd_s_11[k];

        t_10[k] = f_1 * pd_s_12[k]
                  + pb_y[k] * pp_2[k];

        t_11[k] = f_0 * sp_2[k]
                  + pa_z[k] * sd_2[k]
                  + f_1 * pd_s_13[k];
    }
}

auto
compute_prim_pd_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_4 = buffer.data(pd_s + 4);
    const auto *pd_s_5 = buffer.data(pd_s + 5);
    const auto *pd_s_8 = buffer.data(pd_s + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, sp_0, sp_1, sd_0, sd_1, sd_2, pd_s_0, \
                         pd_s_1, pd_s_2, pd_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = pa_x[k] * sd_1[k]
                 + f_1 * pd_s_1[k];

        t_2[k] = pa_x[k] * sd_2[k]
                 + f_1 * pd_s_2[k];

        t_3[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_4[k];
    }

#pragma omp simd aligned(t_4, t_5, pa_y, pa_z, sp_2, sd_2, pd_s_5, \
                         pd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_5[k];

        t_5[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_8[k];
    }
}

auto
compute_prim_pd_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_3 = buffer.data(pd_s + 3);
    const auto *pd_s_4 = buffer.data(pd_s + 4);
    const auto *pd_s_5 = buffer.data(pd_s + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, sp_0, sp_1, sd_0, sd_1, sd_2, \
                         pd_s_0, pd_s_1, pd_s_2, pd_s_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_1[k];

        t_2[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_2[k];

        t_3[k] = pa_z[k] * sd_0[k]
                 + f_1 * pd_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, pa_z, sp_2, sd_1, sd_2, pd_s_4, \
                         pd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * sd_1[k]
                 + f_1 * pd_s_4[k];

        t_5[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_5[k];
    }
}

auto
compute_prim_pd_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_3 = buffer.data(pd_s + 3);
    const auto *pd_s_6 = buffer.data(pd_s + 6);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, sp_0, sp_1, sd_0, sd_1, sd_2, pd_s_0, \
                         pd_s_2, pd_s_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_2[k];

        t_2[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_3[k];
    }

#pragma omp simd aligned(t_3, pa_z, sp_2, sd_2, pd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_6[k];
    }
}

auto
compute_prim_pd_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, sp_0, sp_1, sp_2, sd_0, sd_1, sd_2, \
                         pd_s_0, pd_s_1, pd_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_1[k];

        t_2[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_2[k];
    }
}

auto
compute_prim_pd_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_3 = buffer.data(pd_s + 3);
    const auto *pd_s_4 = buffer.data(pd_s + 4);
    const auto *pd_s_5 = buffer.data(pd_s + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, sp_0, sp_1, sd_0, sd_1, sd_2, pd_s_0, \
                         pd_s_1, pd_s_2, pd_s_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = pa_x[k] * sd_1[k]
                 + f_1 * pd_s_1[k];

        t_2[k] = pa_x[k] * sd_2[k]
                 + f_1 * pd_s_2[k];

        t_3[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, pa_y, pa_z, sp_2, sd_2, pd_s_4, \
                         pd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_4[k];

        t_5[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_5[k];
    }
}

auto
compute_prim_pd_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_4 = buffer.data(pd_s + 4);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, sp_0, sp_1, sd_0, sd_1, sd_2, pd_s_0, \
                         pd_s_1, pd_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_1[k];

        t_2[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_2[k];
    }

#pragma omp simd aligned(t_3, pa_z, sp_2, sd_2, pd_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_4[k];
    }
}

auto
compute_prim_pd_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);

#pragma omp simd aligned(t_0, t_1, pa_y, pa_z, sp_0, sp_1, sd_0, sd_1, pd_s_0, \
                         pd_s_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_y[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_0 * sp_1[k]
                 + pa_z[k] * sd_1[k]
                 + f_1 * pd_s_1[k];
    }
}

auto
compute_prim_pd_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sp, const size_t sd, const size_t pd_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);
    const auto *pd_s_3 = buffer.data(pd_s + 3);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, sp_0, sp_1, sd_0, sd_1, sd_2, pd_s_0, \
                         pd_s_1, pd_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp_0[k]
                 + pa_x[k] * sd_0[k]
                 + f_1 * pd_s_0[k];

        t_1[k] = f_0 * sp_1[k]
                 + pa_y[k] * sd_1[k]
                 + f_1 * pd_s_1[k];

        t_2[k] = pa_y[k] * sd_2[k]
                 + f_1 * pd_s_2[k];
    }

#pragma omp simd aligned(t_3, pa_z, sp_2, sd_2, pd_s_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * sp_2[k]
                 + pa_z[k] * sd_2[k]
                 + f_1 * pd_s_3[k];
    }
}

}  // namespace simdkin
