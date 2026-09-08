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


#include "SimdElectronRepulsionVrrRecSI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_si_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, sg0_0, sg0_1, sg0_2, sg1_0, sg1_1, \
                         sg1_2, sh_0, sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_z[k] * sh_0[k];

        t_2[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];

        t_3[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, sg0_3, sg0_4, sg0_5, sg1_3, sg1_4, sg1_5, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_5[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_6[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pb_x, sg0_6, sg0_8, sg1_6, sg1_8, sh_6, \
                         sh_7, sh_8, sh_9, sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * sg0_6[k]
                 - f_7 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_8[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_7[k];

        t_9[k] = pb_x[k] * sh_8[k];

        t_10[k] = pb_x[k] * sh_9[k];

        t_11[k] = pb_x[k] * sh_10[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sg0_5, sg0_6, sg1_5, sg1_6, \
                         sh_8, sh_9, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * sh_12[k];

        t_13[k] = f_0 * sg0_5[k]
                  - f_1 * sg1_5[k]
                  + pb_y[k] * sh_8[k];

        t_14[k] = pb_z[k] * sh_8[k];

        t_15[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, sg0_7, sg0_8, sg1_7, sg1_8, \
                         sh_10, sh_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_10[k];

        t_17[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_11[k];

        t_18[k] = pb_y[k] * sh_12[k];

        t_19[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sg0_0, sg0_1, sg0_2, sg1_0, sg1_1, sg1_2, sh_0, \
                         sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];

        t_2[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sg0_3, sg0_4, sg0_5, sg1_3, sg1_4, sg1_5, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, sg0_6, sg0_8, sg1_6, sg1_8, sh_6, \
                         sh_7, sh_8, sh_9, sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sg0_6[k]
                 - f_7 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_7[k];

        t_8[k] = pb_x[k] * sh_8[k];

        t_9[k] = pb_x[k] * sh_9[k];

        t_10[k] = pb_x[k] * sh_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, sg0_5, sg0_6, sg0_7, sg1_5, \
                         sg1_6, sg1_7, sh_8, sh_9, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_x[k] * sh_12[k];

        t_12[k] = f_0 * sg0_5[k]
                  - f_1 * sg1_5[k]
                  + pb_y[k] * sh_8[k];

        t_13[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_9[k];

        t_14[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, sg0_8, sg1_8, sh_11, \
                         sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_11[k];

        t_16[k] = pb_y[k] * sh_12[k];

        t_17[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_9 = buffer.data(sg1 + 9);
    const auto *sg1_10 = buffer.data(sg1 + 10);
    const auto *sg1_11 = buffer.data(sg1 + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sg0_0, sg0_1, sg0_2, sg1_0, sg1_3, sg1_4, sh_0, \
                         sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_3[k]
                 + pb_x[k] * sh_1[k];

        t_2[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_4[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sg0_3, sg0_4, sg0_5, sg1_5, sg1_6, sg1_7, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_5[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_6[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_7[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_y, sg0_5, sg0_6, sg0_8, sg1_7, sg1_9, \
                         sg1_11, sh_6, sh_7, sh_8, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sg0_6[k]
                 - f_7 * sg1_9[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_11[k]
                 + pb_x[k] * sh_7[k];

        t_8[k] = f_0 * sg0_5[k]
                 - f_1 * sg1_7[k]
                 + pb_y[k] * sh_8[k];

        t_9[k] = f_2 * sg0_6[k]
                 - f_3 * sg1_9[k]
                 + pb_y[k] * sh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, sg0_7, sg0_8, sg1_10, sg1_11, sh_10, \
                         sh_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_10[k]
                  + pb_y[k] * sh_10[k];

        t_11[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_11[k]
                  + pb_y[k] * sh_11[k];

        t_12[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_11[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg1_0, sg1_1, \
                         sh_0, sh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sg0_2, sg0_3, sg1_2, sg1_3, \
                         sh_3, sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_5[k];

        t_6[k] = pb_z[k] * sh_3[k];

        t_7[k] = pb_y[k] * sh_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sg0_4, sg0_5, sg0_6, sg1_4, sg1_5, \
                         sg1_6, sh_5, sh_6, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_6[k];

        t_9[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_7[k];

        t_10[k] = pb_z[k] * sh_5[k];

        t_11[k] = f_6 * sg0_6[k]
                  - f_7 * sg1_6[k]
                  + pb_x[k] * sh_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, sg0_8, sg1_8, sh_6, \
                         sh_9, sh_10, sh_12, sh_13, sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_y[k] * sh_6[k];

        t_13[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_x[k] * sh_9[k];

        t_14[k] = pb_x[k] * sh_10[k];

        t_15[k] = pb_x[k] * sh_12[k];

        t_16[k] = pb_x[k] * sh_13[k];

        t_17[k] = pb_x[k] * sh_15[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, sg0_5, sg0_6, sg0_7, sg1_5, \
                         sg1_6, sg1_7, sh_10, sh_12, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * sg0_5[k]
                  - f_1 * sg1_5[k]
                  + pb_y[k] * sh_10[k];

        t_19[k] = pb_z[k] * sh_10[k];

        t_20[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_12[k];

        t_21[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_13[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_y, pb_z, sg0_8, sg1_8, sh_14, \
                         sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_14[k];

        t_23[k] = pb_y[k] * sh_15[k];

        t_24[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_15[k];
    }
}

auto
compute_prim_si_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg1_0, sg1_1, \
                         sh_0, sh_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, sg0_2, sg0_3, sg0_4, sg1_2, sg1_3, sg1_4, sh_2, \
                         sh_3, sh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];

        t_5[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_6[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, sg0_5, sg0_6, sg0_8, sg1_5, sg1_6, sg1_8, \
                         sh_5, sh_6, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];

        t_8[k] = f_6 * sg0_6[k]
                 - f_7 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_9[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_7[k];

        t_10[k] = pb_x[k] * sh_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sg0_5, sg1_5, sh_8, \
                         sh_9, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_x[k] * sh_9[k];

        t_12[k] = pb_x[k] * sh_10[k];

        t_13[k] = pb_x[k] * sh_12[k];

        t_14[k] = f_0 * sg0_5[k]
                  - f_1 * sg1_5[k]
                  + pb_y[k] * sh_8[k];

        t_15[k] = pb_z[k] * sh_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, sg0_6, sg0_7, sg0_8, sg1_6, sg1_7, \
                         sg1_8, sh_9, sh_10, sh_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_9[k];

        t_17[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_10[k];

        t_18[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_11[k];

        t_19[k] = pb_y[k] * sh_12[k];
    }

#pragma omp simd aligned(t_20, pb_z, sg0_8, sg1_8, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_9 = buffer.data(sg0 + 9);
    const auto *sg0_10 = buffer.data(sg0 + 10);
    const auto *sg0_11 = buffer.data(sg0 + 11);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_9 = buffer.data(sg1 + 9);
    const auto *sg1_10 = buffer.data(sg1 + 10);
    const auto *sg1_11 = buffer.data(sg1 + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sg0_0, sg0_3, sg0_4, sg1_0, sg1_3, sg1_4, sh_0, \
                         sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_2 * sg0_3[k]
                 - f_3 * sg1_3[k]
                 + pb_x[k] * sh_1[k];

        t_2[k] = f_2 * sg0_4[k]
                 - f_3 * sg1_4[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sg0_5, sg0_6, sg0_7, sg1_5, sg1_6, sg1_7, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sg0_5[k]
                 - f_5 * sg1_5[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_4 * sg0_6[k]
                 - f_5 * sg1_6[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_6 * sg0_7[k]
                 - f_7 * sg1_7[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_y, sg0_7, sg0_9, sg0_11, sg1_7, sg1_9, \
                         sg1_11, sh_6, sh_7, sh_8, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sg0_9[k]
                 - f_7 * sg1_9[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = f_6 * sg0_11[k]
                 - f_7 * sg1_11[k]
                 + pb_x[k] * sh_7[k];

        t_8[k] = f_0 * sg0_7[k]
                 - f_1 * sg1_7[k]
                 + pb_y[k] * sh_8[k];

        t_9[k] = f_2 * sg0_9[k]
                 - f_3 * sg1_9[k]
                 + pb_y[k] * sh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, sg0_10, sg0_11, sg1_10, sg1_11, sh_10, \
                         sh_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * sg0_10[k]
                  - f_5 * sg1_10[k]
                  + pb_y[k] * sh_10[k];

        t_11[k] = f_6 * sg0_11[k]
                  - f_7 * sg1_11[k]
                  + pb_y[k] * sh_11[k];

        t_12[k] = f_0 * sg0_11[k]
                  - f_1 * sg1_11[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_9 = buffer.data(sg1 + 9);
    const auto *sg1_10 = buffer.data(sg1 + 10);
    const auto *sg1_11 = buffer.data(sg1 + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg1_0, sg1_3, \
                         sh_0, sh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_3[k]
                 + pb_x[k] * sh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sg0_2, sg0_3, sg1_4, sg1_5, \
                         sh_3, sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_5[k]
                 + pb_x[k] * sh_5[k];

        t_6[k] = pb_z[k] * sh_3[k];

        t_7[k] = pb_y[k] * sh_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sg0_4, sg0_5, sg0_6, sg1_6, sg1_7, \
                         sg1_9, sh_5, sh_6, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_9[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_7[k]
                 + pb_x[k] * sh_7[k];

        t_10[k] = pb_z[k] * sh_5[k];

        t_11[k] = f_6 * sg0_6[k]
                  - f_7 * sg1_9[k]
                  + pb_x[k] * sh_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, sg0_5, sg0_8, sg1_7, \
                         sg1_11, sh_6, sh_9, sh_10, sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_y[k] * sh_6[k];

        t_13[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_11[k]
                  + pb_x[k] * sh_9[k];

        t_14[k] = pb_x[k] * sh_10[k];

        t_15[k] = pb_x[k] * sh_15[k];

        t_16[k] = f_0 * sg0_5[k]
                  - f_1 * sg1_7[k]
                  + pb_y[k] * sh_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, sg0_6, sg0_7, sg0_8, sg1_9, \
                         sg1_10, sg1_11, sh_10, sh_12, sh_13, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pb_z[k] * sh_10[k];

        t_18[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_9[k]
                  + pb_y[k] * sh_12[k];

        t_19[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_10[k]
                  + pb_y[k] * sh_13[k];

        t_20[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_11[k]
                  + pb_y[k] * sh_14[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_y, pb_z, sg0_8, sg1_11, sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_y[k] * sh_15[k];

        t_22[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_11[k]
                  + pb_z[k] * sh_15[k];
    }
}

auto
compute_prim_si_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_9 = buffer.data(sg0 + 9);
    const auto *sg0_10 = buffer.data(sg0 + 10);
    const auto *sg0_11 = buffer.data(sg0 + 11);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sg0_0, sg0_3, sg0_4, sg1_0, sg1_1, sg1_2, sh_0, \
                         sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_2 * sg0_3[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];

        t_2[k] = f_2 * sg0_4[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sg0_5, sg0_6, sg0_7, sg1_3, sg1_4, sg1_5, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sg0_5[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_4 * sg0_6[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_6 * sg0_7[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, sg0_9, sg0_11, sg1_6, sg1_8, sh_6, \
                         sh_7, sh_8, sh_9, sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sg0_9[k]
                 - f_7 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = f_6 * sg0_11[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_7[k];

        t_8[k] = pb_x[k] * sh_8[k];

        t_9[k] = pb_x[k] * sh_9[k];

        t_10[k] = pb_x[k] * sh_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, sg0_7, sg0_9, sg0_10, sg1_5, \
                         sg1_6, sg1_7, sh_8, sh_9, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_x[k] * sh_12[k];

        t_12[k] = f_0 * sg0_7[k]
                  - f_1 * sg1_5[k]
                  + pb_y[k] * sh_8[k];

        t_13[k] = f_2 * sg0_9[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_9[k];

        t_14[k] = f_4 * sg0_10[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, sg0_11, sg1_8, sh_11, \
                         sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * sg0_11[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_11[k];

        t_16[k] = pb_y[k] * sh_12[k];

        t_17[k] = f_0 * sg0_11[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_9 = buffer.data(sg0 + 9);
    const auto *sg0_10 = buffer.data(sg0 + 10);
    const auto *sg0_11 = buffer.data(sg0 + 11);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_9 = buffer.data(sg1 + 9);
    const auto *sg1_10 = buffer.data(sg1 + 10);
    const auto *sg1_11 = buffer.data(sg1 + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_3, sg1_0, sg1_3, \
                         sh_0, sh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_3[k]
                 - f_3 * sg1_3[k]
                 + pb_x[k] * sh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_z, sg0_4, sg0_5, sg0_6, sg1_4, sg1_5, \
                         sg1_6, sh_3, sh_4, sh_5, sh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_4[k]
                 - f_3 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_4 * sg0_5[k]
                 - f_5 * sg1_5[k]
                 + pb_x[k] * sh_5[k];

        t_6[k] = pb_z[k] * sh_3[k];

        t_7[k] = f_4 * sg0_6[k]
                 - f_5 * sg1_6[k]
                 + pb_x[k] * sh_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sg0_7, sg0_9, sg0_11, sg1_7, sg1_9, \
                         sg1_11, sh_5, sh_7, sh_8, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * sg0_7[k]
                 - f_7 * sg1_7[k]
                 + pb_x[k] * sh_7[k];

        t_9[k] = pb_z[k] * sh_5[k];

        t_10[k] = f_6 * sg0_9[k]
                  - f_7 * sg1_9[k]
                  + pb_x[k] * sh_8[k];

        t_11[k] = f_6 * sg0_11[k]
                  - f_7 * sg1_11[k]
                  + pb_x[k] * sh_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_y, pb_z, sg0_7, sg0_9, sg0_10, sg1_7, \
                         sg1_9, sg1_10, sh_10, sh_12, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sg0_7[k]
                  - f_1 * sg1_7[k]
                  + pb_y[k] * sh_10[k];

        t_13[k] = pb_z[k] * sh_10[k];

        t_14[k] = f_2 * sg0_9[k]
                  - f_3 * sg1_9[k]
                  + pb_y[k] * sh_12[k];

        t_15[k] = f_4 * sg0_10[k]
                  - f_5 * sg1_10[k]
                  + pb_y[k] * sh_13[k];
    }

#pragma omp simd aligned(t_16, t_17, pb_y, pb_z, sg0_11, sg1_11, sh_14, \
                         sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * sg0_11[k]
                  - f_7 * sg1_11[k]
                  + pb_y[k] * sh_14[k];

        t_17[k] = f_0 * sg0_11[k]
                  - f_1 * sg1_11[k]
                  + pb_z[k] * sh_15[k];
    }
}

auto
compute_prim_si_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sg0, const size_t sg1, const size_t sh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_8 = buffer.data(sg1 + 8);
    const auto *sg1_9 = buffer.data(sg1 + 9);
    const auto *sg1_10 = buffer.data(sg1 + 10);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg1_0, sg1_2, \
                         sh_0, sh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_z, sg0_2, sg0_3, sg0_4, sg1_3, sg1_4, \
                         sg1_5, sh_3, sh_4, sh_5, sh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_3[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_5[k];

        t_6[k] = pb_z[k] * sh_3[k];

        t_7[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_5[k]
                 + pb_x[k] * sh_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sg0_5, sg0_6, sg0_8, sg1_6, sg1_8, \
                         sg1_10, sh_5, sh_7, sh_8, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_6[k]
                 + pb_x[k] * sh_7[k];

        t_9[k] = pb_z[k] * sh_5[k];

        t_10[k] = f_6 * sg0_6[k]
                  - f_7 * sg1_8[k]
                  + pb_x[k] * sh_8[k];

        t_11[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_10[k]
                  + pb_x[k] * sh_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, sg0_5, sg1_6, \
                         sh_10, sh_12, sh_13, sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * sh_10[k];

        t_13[k] = pb_x[k] * sh_12[k];

        t_14[k] = pb_x[k] * sh_13[k];

        t_15[k] = pb_x[k] * sh_15[k];

        t_16[k] = f_0 * sg0_5[k]
                  - f_1 * sg1_6[k]
                  + pb_y[k] * sh_10[k];

        t_17[k] = pb_z[k] * sh_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, sg0_6, sg0_7, sg0_8, sg1_8, sg1_9, \
                         sg1_10, sh_12, sh_13, sh_14, sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_8[k]
                  + pb_y[k] * sh_12[k];

        t_19[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_9[k]
                  + pb_y[k] * sh_13[k];

        t_20[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_10[k]
                  + pb_y[k] * sh_14[k];

        t_21[k] = pb_y[k] * sh_15[k];
    }

#pragma omp simd aligned(t_22, pb_z, sg0_8, sg1_10, sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_10[k]
                  + pb_z[k] * sh_15[k];
    }
}

auto
compute_prim_si_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_8 = buffer.data(sg0 + 8);
    const auto *sg0_9 = buffer.data(sg0 + 9);
    const auto *sg0_10 = buffer.data(sg0 + 10);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_2, sg1_0, sg1_1, \
                         sh_0, sh_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, sg0_3, sg0_4, sg0_5, sg1_2, sg1_3, sg1_4, sh_2, \
                         sh_3, sh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_3[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];

        t_5[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_6[k] = f_4 * sg0_5[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, sg0_6, sg0_8, sg0_10, sg1_5, sg1_6, sg1_8, \
                         sh_5, sh_6, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * sg0_6[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];

        t_8[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_9[k] = f_6 * sg0_10[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_7[k];

        t_10[k] = pb_x[k] * sh_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sg0_6, sg1_5, sh_8, \
                         sh_9, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_x[k] * sh_9[k];

        t_12[k] = pb_x[k] * sh_10[k];

        t_13[k] = pb_x[k] * sh_12[k];

        t_14[k] = f_0 * sg0_6[k]
                  - f_1 * sg1_5[k]
                  + pb_y[k] * sh_8[k];

        t_15[k] = pb_z[k] * sh_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, sg0_8, sg0_9, sg0_10, sg1_6, sg1_7, \
                         sg1_8, sh_9, sh_10, sh_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * sg0_8[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_9[k];

        t_17[k] = f_4 * sg0_9[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_10[k];

        t_18[k] = f_6 * sg0_10[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_11[k];

        t_19[k] = pb_y[k] * sh_12[k];
    }

#pragma omp simd aligned(t_20, pb_z, sg0_10, sg1_8, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * sg0_10[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_12[k];
    }
}

auto
compute_prim_si_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg0_2, sg1_0, sg1_1, \
                         sg1_2, sh_0, sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_0 * sg0_1[k]
                 - f_1 * sg1_1[k]
                 + pb_y[k] * sh_1[k];

        t_2[k] = f_0 * sg0_2[k]
                 - f_1 * sg1_2[k]
                 + pb_z[k] * sh_2[k];
    }
}

auto
compute_prim_si_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_11 = buffer.data(sg1 + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_11 = buffer.data(sh + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg0_2, sg1_0, sg1_7, \
                         sg1_11, sh_0, sh_7, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_0 * sg0_1[k]
                 - f_1 * sg1_7[k]
                 + pb_y[k] * sh_7[k];

        t_2[k] = f_0 * sg0_2[k]
                 - f_1 * sg1_11[k]
                 + pb_z[k] * sh_11[k];
    }
}

auto
compute_prim_si_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_9 = buffer.data(sg0 + 9);
    const auto *sg0_10 = buffer.data(sg0 + 10);
    const auto *sg0_11 = buffer.data(sg0 + 11);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sg0_0, sg0_3, sg0_4, sg1_0, sg1_1, sg1_2, sh_0, \
                         sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_2 * sg0_3[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];

        t_2[k] = f_2 * sg0_4[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sg0_5, sg0_6, sg0_7, sg1_3, sg1_4, sg1_5, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sg0_5[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_4 * sg0_6[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_6 * sg0_7[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, sg0_7, sg0_9, sg0_11, sg1_5, sg1_6, sg1_8, \
                         sh_6, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sg0_11[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = f_0 * sg0_7[k]
                 - f_1 * sg1_5[k]
                 + pb_y[k] * sh_7[k];

        t_8[k] = f_2 * sg0_9[k]
                 - f_3 * sg1_6[k]
                 + pb_y[k] * sh_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sg0_10, sg0_11, sg1_7, sg1_8, sh_9, \
                         sh_10, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * sg0_10[k]
                 - f_5 * sg1_7[k]
                 + pb_y[k] * sh_9[k];

        t_10[k] = f_6 * sg0_11[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_10[k];

        t_11[k] = f_0 * sg0_11[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_11[k];
    }
}

auto
compute_prim_si_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_14 = buffer.data(sh + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg0_2, sg1_0, sg1_5, \
                         sg1_8, sh_0, sh_9, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_0 * sg0_1[k]
                 - f_1 * sg1_5[k]
                 + pb_y[k] * sh_9[k];

        t_2[k] = f_0 * sg0_2[k]
                 - f_1 * sg1_8[k]
                 + pb_z[k] * sh_14[k];
    }
}

auto
compute_prim_si_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sg0_0, sg0_1, sg1_0, sg1_1, \
                         sh_0, sh_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * sh_0[k];

        t_2[k] = pb_z[k] * sh_0[k];

        t_3[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, sg0_2, sg0_3, sg0_4, sg1_2, sg1_3, sg1_4, sh_2, \
                         sh_3, sh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];

        t_5[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_6[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, sg0_5, sg0_8, sg1_5, sg1_8, \
                         sh_5, sh_6, sh_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];

        t_8[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_6[k];

        t_9[k] = f_0 * sg0_5[k]
                 - f_1 * sg1_5[k]
                 + pb_y[k] * sh_7[k];

        t_10[k] = pb_z[k] * sh_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_y, pb_z, sg0_6, sg0_7, sg0_8, sg1_6, \
                         sg1_7, sg1_8, sh_8, sh_9, sh_10, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_8[k];

        t_12[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_9[k];

        t_13[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_10[k];

        t_14[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_11[k];
    }
}

auto
compute_prim_si_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sg0, const size_t sg1, const size_t sh,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / beta;
    const auto f_1 = 2.5 * alpha / (beta * p);
    const auto f_2 = 1.5 / beta;
    const auto f_3 = 1.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / beta;
    const auto f_5 = alpha / (beta * p);
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_1 = buffer.data(sg0 + 1);
    const auto *sg0_2 = buffer.data(sg0 + 2);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_4 = buffer.data(sg0 + 4);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_7 = buffer.data(sg0 + 7);
    const auto *sg0_8 = buffer.data(sg0 + 8);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_1 = buffer.data(sg1 + 1);
    const auto *sg1_2 = buffer.data(sg1 + 2);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_4 = buffer.data(sg1 + 4);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_7 = buffer.data(sg1 + 7);
    const auto *sg1_8 = buffer.data(sg1 + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sg0_0, sg0_1, sg0_2, sg1_0, sg1_1, sg1_2, sh_0, \
                         sh_1, sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg0_0[k]
                 - f_1 * sg1_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_2 * sg0_1[k]
                 - f_3 * sg1_1[k]
                 + pb_x[k] * sh_1[k];

        t_2[k] = f_2 * sg0_2[k]
                 - f_3 * sg1_2[k]
                 + pb_x[k] * sh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sg0_3, sg0_4, sg0_5, sg1_3, sg1_4, sg1_5, sh_3, \
                         sh_4, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sg0_3[k]
                 - f_5 * sg1_3[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_4 * sg0_4[k]
                 - f_5 * sg1_4[k]
                 + pb_x[k] * sh_4[k];

        t_5[k] = f_6 * sg0_5[k]
                 - f_7 * sg1_5[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_y, sg0_5, sg0_8, sg1_5, sg1_8, sh_6, \
                         sh_7, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sg0_8[k]
                 - f_7 * sg1_8[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = pb_x[k] * sh_7[k];

        t_8[k] = pb_x[k] * sh_11[k];

        t_9[k] = f_0 * sg0_5[k]
                 - f_1 * sg1_5[k]
                 + pb_y[k] * sh_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, sg0_6, sg0_7, sg0_8, sg1_6, sg1_7, \
                         sg1_8, sh_8, sh_9, sh_10, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * sg0_6[k]
                  - f_3 * sg1_6[k]
                  + pb_y[k] * sh_8[k];

        t_11[k] = f_4 * sg0_7[k]
                  - f_5 * sg1_7[k]
                  + pb_y[k] * sh_9[k];

        t_12[k] = f_6 * sg0_8[k]
                  - f_7 * sg1_8[k]
                  + pb_y[k] * sh_10[k];

        t_13[k] = pb_y[k] * sh_11[k];
    }

#pragma omp simd aligned(t_14, pb_z, sg0_8, sg1_8, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * sg0_8[k]
                  - f_1 * sg1_8[k]
                  + pb_z[k] * sh_11[k];
    }
}

}  // namespace simdt2ceri
