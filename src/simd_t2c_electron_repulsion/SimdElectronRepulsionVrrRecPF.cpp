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


#include "SimdElectronRepulsionVrrRecPF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_pf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sd, const size_t sf,
                                     const size_t pd, const size_t ncols,
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

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_x, pb_y, pb_z, sd_0, sd_1, sf_0, \
                         pd_0, pd_1, pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * pd_0[k];

        t_2[k] = pb_z[k] * pd_0[k];

        t_3[k] = f_1 * sd_1[k]
                 + pb_x[k] * pd_2[k];

        t_4[k] = pb_y[k] * pd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_x, pb_y, pb_z, sd_2, sf_1, sf_2, \
                         pd_2, pd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sd_2[k]
                 + pb_x[k] * pd_3[k];

        t_6[k] = pa_x[k] * sf_1[k];

        t_7[k] = pb_z[k] * pd_2[k];

        t_8[k] = pb_y[k] * pd_3[k];

        t_9[k] = pa_x[k] * sf_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, pb_z, sd_0, \
                         sf_0, pd_4, pd_5, pd_6, pd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * sf_0[k];

        t_11[k] = f_1 * sd_0[k]
                  + pb_y[k] * pd_4[k];

        t_12[k] = pb_z[k] * pd_4[k];

        t_13[k] = pb_x[k] * pd_5[k];

        t_14[k] = pb_x[k] * pd_6[k];

        t_15[k] = pb_x[k] * pd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pa_z, pb_y, pb_z, sd_1, sd_2, \
                         sf_0, sf_1, sf_2, pd_5, pd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * sd_1[k]
                  + pa_y[k] * sf_1[k];

        t_17[k] = pb_z[k] * pd_5[k];

        t_18[k] = f_1 * sd_2[k]
                  + pb_y[k] * pd_7[k];

        t_19[k] = pa_y[k] * sf_2[k];

        t_20[k] = pa_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, pb_z, sd_0, \
                         sf_1, pd_8, pd_9, pd_10, pd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_y[k] * pd_8[k];

        t_22[k] = f_1 * sd_0[k]
                  + pb_z[k] * pd_8[k];

        t_23[k] = pb_x[k] * pd_9[k];

        t_24[k] = pb_x[k] * pd_10[k];

        t_25[k] = pb_x[k] * pd_11[k];

        t_26[k] = pa_z[k] * sf_1[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_z, pb_y, pb_z, sd_1, sd_2, sf_2, pd_9, \
                         pd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * sd_1[k]
                  + pb_z[k] * pd_9[k];

        t_28[k] = pb_y[k] * pd_11[k];

        t_29[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_2[k];
    }
}

auto
compute_prim_pf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sd, const size_t sf,
                                     const size_t pd, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_x, pb_z, sd_0, sd_1, sd_2, sf_0, \
                         sf_3, pd_0, pd_1, pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pb_z[k] * pd_0[k];

        t_2[k] = f_1 * sd_1[k]
                 + pb_x[k] * pd_1[k];

        t_3[k] = f_1 * sd_2[k]
                 + pb_x[k] * pd_2[k];

        t_4[k] = pa_x[k] * sf_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pb_y, pb_z, sd_0, sd_1, sf_0, \
                         sf_3, sf_5, pd_3, pd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * sf_5[k];

        t_6[k] = pa_y[k] * sf_0[k];

        t_7[k] = f_1 * sd_0[k]
                 + pb_y[k] * pd_3[k];

        t_8[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];

        t_9[k] = pb_z[k] * pd_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_y, pb_z, sd_0, sd_2, \
                         sf_0, sf_3, sf_5, pd_5, pd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * sd_2[k]
                  + pb_y[k] * pd_5[k];

        t_11[k] = pa_y[k] * sf_5[k];

        t_12[k] = pa_z[k] * sf_0[k];

        t_13[k] = f_1 * sd_0[k]
                  + pb_z[k] * pd_6[k];

        t_14[k] = pa_z[k] * sf_3[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_z, pb_y, pb_z, sd_1, sd_2, sf_5, pd_7, \
                         pd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * sd_1[k]
                  + pb_z[k] * pd_7[k];

        t_16[k] = pb_y[k] * pd_8[k];

        t_17[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_5[k];
    }
}

auto
compute_prim_pf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sd, const size_t sf, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

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

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, sd_0, sd_1, \
                         sf_0, sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pa_x[k] * sf_3[k];

        t_2[k] = pa_x[k] * sf_5[k];

        t_3[k] = pa_y[k] * sf_0[k];

        t_4[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];

        t_5[k] = pa_y[k] * sf_5[k];

        t_6[k] = pa_z[k] * sf_0[k];

        t_7[k] = pa_z[k] * sf_3[k];
    }

#pragma omp simd aligned(t_8, pa_z, sd_2, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_5[k];
    }
}

auto
compute_prim_pf_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sd, const size_t sf,
                                     const size_t pd, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_x, pb_y, sd_0, sd_1, sd_2, sf_0, pd_1, \
                         pd_2, pd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_1 * sd_1[k]
                 + pb_x[k] * pd_1[k];

        t_2[k] = f_1 * sd_2[k]
                 + pb_x[k] * pd_2[k];

        t_3[k] = f_1 * sd_0[k]
                 + pb_y[k] * pd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pb_y, sd_1, sd_2, sf_0, sf_1, sf_2, \
                         pd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_1[k];

        t_5[k] = f_1 * sd_2[k]
                 + pb_y[k] * pd_5[k];

        t_6[k] = pa_y[k] * sf_2[k];

        t_7[k] = pa_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_z, pb_z, sd_0, sd_1, sd_2, sf_1, sf_2, pd_6, \
                         pd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * sd_0[k]
                 + pb_z[k] * pd_6[k];

        t_9[k] = pa_z[k] * sf_1[k];

        t_10[k] = f_1 * sd_1[k]
                  + pb_z[k] * pd_7[k];

        t_11[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_2[k];
    }
}

auto
compute_prim_pf_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sd, const size_t sf,
                                     const size_t pd, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pb_x, sd_0, sd_1, sf_0, \
                         sf_5, sf_8, pd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pa_x[k] * sf_5[k];

        t_2[k] = pa_x[k] * sf_8[k];

        t_3[k] = pa_y[k] * sf_0[k];

        t_4[k] = pb_x[k] * pd_4[k];

        t_5[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pa_z, pb_y, pb_z, sd_0, sd_2, sf_0, sf_8, \
                         pd_5, pd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sd_2[k]
                 + pb_y[k] * pd_5[k];

        t_7[k] = pa_y[k] * sf_8[k];

        t_8[k] = pa_z[k] * sf_0[k];

        t_9[k] = f_1 * sd_0[k]
                 + pb_z[k] * pd_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_z, pb_x, pb_y, pb_z, sd_1, sd_2, \
                         sf_5, sf_8, pd_7, pd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * pd_8[k];

        t_11[k] = pa_z[k] * sf_5[k];

        t_12[k] = f_1 * sd_1[k]
                  + pb_z[k] * pd_7[k];

        t_13[k] = pb_y[k] * pd_8[k];

        t_14[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_8[k];
    }
}

auto
compute_prim_pf_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sd, const size_t sf,
                                     const size_t pd, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pb_y, sd_0, sd_1, sd_2, \
                         sf_0, sf_3, sf_5, pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pa_x[k] * sf_3[k];

        t_2[k] = pa_x[k] * sf_5[k];

        t_3[k] = pa_y[k] * sf_0[k];

        t_4[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];

        t_5[k] = f_1 * sd_2[k]
                 + pb_y[k] * pd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pa_z, pb_z, sd_0, sd_1, sf_0, sf_3, \
                         sf_5, pd_3, pd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * sf_5[k];

        t_7[k] = pa_z[k] * sf_0[k];

        t_8[k] = f_1 * sd_0[k]
                 + pb_z[k] * pd_3[k];

        t_9[k] = pa_z[k] * sf_3[k];

        t_10[k] = f_1 * sd_1[k]
                  + pb_z[k] * pd_4[k];
    }

#pragma omp simd aligned(t_11, pa_z, sd_2, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_5[k];
    }
}

auto
compute_prim_pf_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sd, const size_t sf,
                                     const size_t pd, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_x, pa_y, pb_x, sd_0, sd_1, \
                         sf_0, sf_5, sf_8, pd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pa_x[k] * sf_5[k];

        t_2[k] = pa_x[k] * sf_8[k];

        t_3[k] = pa_y[k] * sf_0[k];

        t_4[k] = pb_x[k] * pd_1[k];

        t_5[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_5[k];

        t_6[k] = pa_y[k] * sf_8[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pa_z, pb_x, pb_y, sd_2, sf_0, sf_5, sf_8, \
                         pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * sf_0[k];

        t_8[k] = pb_x[k] * pd_2[k];

        t_9[k] = pa_z[k] * sf_5[k];

        t_10[k] = pb_y[k] * pd_2[k];

        t_11[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_8[k];
    }
}

auto
compute_prim_pf_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sd, const size_t sf, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pa_z, sd_0, sd_1, sd_2, \
                         sf_0, sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];

        t_2[k] = pa_y[k] * sf_5[k];

        t_3[k] = pa_z[k] * sf_0[k];

        t_4[k] = pa_z[k] * sf_3[k];

        t_5[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_5[k];
    }
}

auto
compute_prim_pf_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sd, const size_t sf, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, sd_0, sd_1, sd_2, sf_0, sf_3, \
                         sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];

        t_2[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_5[k];
    }
}

auto
compute_prim_pf_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sd, const size_t sf, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, sd_0, sd_1, sd_2, sf_0, sf_5, \
                         sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_5[k];

        t_2[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_8[k];
    }
}

auto
compute_prim_pf_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t sd, const size_t sf, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, sd_0, sd_1, sd_2, sf_0, sf_1, \
                         sf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_1[k];

        t_2[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_2[k];
    }
}

auto
compute_prim_pf_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t sd, const size_t sf,
                                      const size_t pd, const size_t ncols,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pb_x, sd_0, sd_1, sf_0, \
                         sf_3, sf_5, pd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pa_x[k] * sf_3[k];

        t_2[k] = pa_x[k] * sf_5[k];

        t_3[k] = pa_y[k] * sf_0[k];

        t_4[k] = pb_x[k] * pd_1[k];

        t_5[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pa_z, pb_y, pb_z, sd_0, sd_2, sf_0, sf_5, \
                         pd_2, pd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sd_2[k]
                 + pb_y[k] * pd_2[k];

        t_7[k] = pa_y[k] * sf_5[k];

        t_8[k] = pa_z[k] * sf_0[k];

        t_9[k] = f_1 * sd_0[k]
                 + pb_z[k] * pd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_z, pb_x, pb_y, pb_z, sd_1, sd_2, \
                         sf_3, sf_5, pd_4, pd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * pd_5[k];

        t_11[k] = pa_z[k] * sf_3[k];

        t_12[k] = f_1 * sd_1[k]
                  + pb_z[k] * pd_4[k];

        t_13[k] = pb_y[k] * pd_5[k];

        t_14[k] = f_0 * sd_2[k]
                  + pa_z[k] * sf_5[k];
    }
}

auto
compute_prim_pf_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t sd, const size_t sf,
                                      const size_t pd, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, sd_0, sd_1, sf_0, \
                         sf_5, sf_8, pd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pb_x[k] * pd_1[k];

        t_2[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_5[k];

        t_3[k] = pa_y[k] * sf_8[k];

        t_4[k] = pa_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_z, pb_x, pb_y, sd_2, sf_5, sf_8, \
                         pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_x[k] * pd_2[k];

        t_6[k] = pa_z[k] * sf_5[k];

        t_7[k] = pb_y[k] * pd_2[k];

        t_8[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_8[k];
    }
}

auto
compute_prim_pf_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t sd, const size_t sf,
                                      const size_t pd, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *pd_2 = buffer.data(pd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_y, sd_0, sd_1, sd_2, sf_0, \
                         sf_5, sf_8, pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_5[k];

        t_2[k] = pb_y[k] * pd_2[k];

        t_3[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_8[k];
    }
}

auto
compute_prim_pf_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t sd, const size_t sf, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_7 = buffer.data(sf + 7);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, sd_0, sd_1, sd_2, sf_0, sf_4, \
                         sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_4[k];

        t_2[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_7[k];
    }
}

auto
compute_prim_pf_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t sd, const size_t sf,
                                      const size_t pd, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, sd_0, sd_1, sf_0, \
                         sf_3, sf_5, pd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k];

        t_1[k] = pb_x[k] * pd_1[k];

        t_2[k] = f_0 * sd_1[k]
                 + pa_y[k] * sf_3[k];

        t_3[k] = pa_y[k] * sf_5[k];

        t_4[k] = pa_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_z, pb_x, pb_y, sd_2, sf_3, sf_5, \
                         pd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_x[k] * pd_2[k];

        t_6[k] = pa_z[k] * sf_3[k];

        t_7[k] = pb_y[k] * pd_2[k];

        t_8[k] = f_0 * sd_2[k]
                 + pa_z[k] * sf_5[k];
    }
}

}  // namespace simdt2ceri
