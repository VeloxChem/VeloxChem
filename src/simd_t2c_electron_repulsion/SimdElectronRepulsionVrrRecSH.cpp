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


#include "SimdElectronRepulsionVrrRecSH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_sh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

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

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_1 = buffer.data(sf1 + 1);
    const auto *sf1_2 = buffer.data(sf1 + 2);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, sf0_0, sf0_1, sf0_2, sf1_0, sf1_1, \
                         sf1_2, sg_0, sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = pb_z[k] * sg_0[k];

        t_2[k] = f_2 * sf0_1[k]
                 - f_3 * sf1_1[k]
                 + pb_x[k] * sg_1[k];

        t_3[k] = f_2 * sf0_2[k]
                 - f_3 * sf1_2[k]
                 + pb_x[k] * sg_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, sf0_3, sf0_5, sf1_3, sf1_5, sg_3, \
                         sg_4, sg_5, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * sf0_3[k]
                 - f_5 * sf1_3[k]
                 + pb_x[k] * sg_3[k];

        t_5[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_4[k];

        t_6[k] = pb_x[k] * sg_5[k];

        t_7[k] = pb_x[k] * sg_6[k];

        t_8[k] = pb_x[k] * sg_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, sf0_3, sf0_4, sf0_5, sf1_3, sf1_4, \
                         sf1_5, sg_5, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * sf0_3[k]
                 - f_1 * sf1_3[k]
                 + pb_y[k] * sg_5[k];

        t_10[k] = pb_z[k] * sg_5[k];

        t_11[k] = f_2 * sf0_4[k]
                  - f_3 * sf1_4[k]
                  + pb_y[k] * sg_6[k];

        t_12[k] = f_4 * sf0_5[k]
                  - f_5 * sf1_5[k]
                  + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_y, pb_z, sf0_5, sf1_5, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sg_8[k];

        t_14[k] = f_0 * sf0_5[k]
                  - f_1 * sf1_5[k]
                  + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);
    const auto *sf1_7 = buffer.data(sf1 + 7);
    const auto *sf1_8 = buffer.data(sf1 + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf0_0, sf0_1, sf0_2, sf1_0, sf1_3, sf1_4, sg_0, \
                         sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_2 * sf0_1[k]
                 - f_3 * sf1_3[k]
                 + pb_x[k] * sg_1[k];

        t_2[k] = f_2 * sf0_2[k]
                 - f_3 * sf1_4[k]
                 + pb_x[k] * sg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pb_x, pb_y, sf0_3, sf0_4, sf0_5, sf1_5, sf1_7, \
                         sf1_8, sg_3, sg_4, sg_5, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sf0_3[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_8[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = f_0 * sf0_3[k]
                 - f_1 * sf1_5[k]
                 + pb_y[k] * sg_5[k];

        t_6[k] = f_2 * sf0_4[k]
                 - f_3 * sf1_7[k]
                 + pb_y[k] * sg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, sf0_5, sf1_8, sg_7, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_8[k]
                 + pb_y[k] * sg_7[k];

        t_8[k] = f_0 * sf0_5[k]
                 - f_1 * sf1_8[k]
                 + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

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

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_1 = buffer.data(sf1 + 1);
    const auto *sf1_2 = buffer.data(sf1 + 2);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sf0_0, sf0_1, sf1_0, sf1_1, \
                         sg_0, sg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * sg_0[k];

        t_2[k] = pb_z[k] * sg_0[k];

        t_3[k] = f_2 * sf0_1[k]
                 - f_3 * sf1_1[k]
                 + pb_x[k] * sg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sf0_2, sf0_3, sf1_2, sf1_3, \
                         sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sf0_2[k]
                 - f_3 * sf1_2[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = f_4 * sf0_3[k]
                 - f_5 * sf1_3[k]
                 + pb_x[k] * sg_5[k];

        t_6[k] = pb_z[k] * sg_3[k];

        t_7[k] = pb_y[k] * sg_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pb_x, pb_y, sf0_3, sf0_5, sf1_3, sf1_5, \
                         sg_6, sg_7, sg_9, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_6[k];

        t_9[k] = pb_x[k] * sg_7[k];

        t_10[k] = pb_x[k] * sg_9[k];

        t_11[k] = pb_x[k] * sg_11[k];

        t_12[k] = f_0 * sf0_3[k]
                  - f_1 * sf1_3[k]
                  + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pb_y, pb_z, sf0_4, sf0_5, sf1_4, sf1_5, \
                         sg_7, sg_9, sg_10, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_z[k] * sg_7[k];

        t_14[k] = f_2 * sf0_4[k]
                  - f_3 * sf1_4[k]
                  + pb_y[k] * sg_9[k];

        t_15[k] = f_4 * sf0_5[k]
                  - f_5 * sf1_5[k]
                  + pb_y[k] * sg_10[k];

        t_16[k] = pb_y[k] * sg_11[k];

        t_17[k] = f_0 * sf0_5[k]
                  - f_1 * sf1_5[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

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

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_1 = buffer.data(sf1 + 1);
    const auto *sf1_2 = buffer.data(sf1 + 2);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf0_0, sf0_1, sf0_2, sf1_0, sf1_1, sf1_2, sg_0, \
                         sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_2 * sf0_1[k]
                 - f_3 * sf1_1[k]
                 + pb_x[k] * sg_1[k];

        t_2[k] = f_2 * sf0_2[k]
                 - f_3 * sf1_2[k]
                 + pb_x[k] * sg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, t_7, pb_x, sf0_3, sf0_5, sf1_3, sf1_5, sg_3, \
                         sg_4, sg_5, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sf0_3[k]
                 - f_5 * sf1_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = pb_x[k] * sg_5[k];

        t_6[k] = pb_x[k] * sg_6[k];

        t_7[k] = pb_x[k] * sg_8[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, sf0_3, sf0_4, sf0_5, sf1_3, sf1_4, sf1_5, \
                         sg_5, sg_6, sg_7, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sf0_3[k]
                 - f_1 * sf1_3[k]
                 + pb_y[k] * sg_5[k];

        t_9[k] = f_2 * sf0_4[k]
                 - f_3 * sf1_4[k]
                 + pb_y[k] * sg_6[k];

        t_10[k] = f_4 * sf0_5[k]
                  - f_5 * sf1_5[k]
                  + pb_y[k] * sg_7[k];

        t_11[k] = pb_y[k] * sg_8[k];
    }

#pragma omp simd aligned(t_12, pb_z, sf0_5, sf1_5, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sf0_5[k]
                  - f_1 * sf1_5[k]
                  + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);
    const auto *sf0_7 = buffer.data(sf0 + 7);
    const auto *sf0_8 = buffer.data(sf0 + 8);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);
    const auto *sf1_7 = buffer.data(sf1 + 7);
    const auto *sf1_8 = buffer.data(sf1 + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf0_0, sf0_3, sf0_4, sf1_0, sf1_3, sf1_4, sg_0, \
                         sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_2 * sf0_3[k]
                 - f_3 * sf1_3[k]
                 + pb_x[k] * sg_1[k];

        t_2[k] = f_2 * sf0_4[k]
                 - f_3 * sf1_4[k]
                 + pb_x[k] * sg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pb_x, pb_y, sf0_5, sf0_7, sf0_8, sf1_5, sf1_7, \
                         sf1_8, sg_3, sg_4, sg_5, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = f_4 * sf0_8[k]
                 - f_5 * sf1_8[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = f_0 * sf0_5[k]
                 - f_1 * sf1_5[k]
                 + pb_y[k] * sg_5[k];

        t_6[k] = f_2 * sf0_7[k]
                 - f_3 * sf1_7[k]
                 + pb_y[k] * sg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, sf0_8, sf1_8, sg_7, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * sf0_8[k]
                 - f_5 * sf1_8[k]
                 + pb_y[k] * sg_7[k];

        t_8[k] = f_0 * sf0_8[k]
                 - f_1 * sf1_8[k]
                 + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);
    const auto *sf1_7 = buffer.data(sf1 + 7);
    const auto *sf1_8 = buffer.data(sf1 + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sf0_0, sf0_1, sf1_0, sf1_3, \
                         sg_0, sg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * sg_0[k];

        t_2[k] = pb_z[k] * sg_0[k];

        t_3[k] = f_2 * sf0_1[k]
                 - f_3 * sf1_3[k]
                 + pb_x[k] * sg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sf0_2, sf0_3, sf1_4, sf1_5, \
                         sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sf0_2[k]
                 - f_3 * sf1_4[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = f_4 * sf0_3[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_5[k];

        t_6[k] = pb_z[k] * sg_3[k];

        t_7[k] = pb_y[k] * sg_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pb_x, pb_y, pb_z, sf0_3, sf0_5, sf1_5, \
                         sf1_8, sg_6, sg_7, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_8[k]
                 + pb_x[k] * sg_6[k];

        t_9[k] = pb_x[k] * sg_7[k];

        t_10[k] = pb_x[k] * sg_11[k];

        t_11[k] = f_0 * sf0_3[k]
                  - f_1 * sf1_5[k]
                  + pb_y[k] * sg_7[k];

        t_12[k] = pb_z[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, sf0_4, sf0_5, sf1_7, sf1_8, sg_9, \
                         sg_10, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * sf0_4[k]
                  - f_3 * sf1_7[k]
                  + pb_y[k] * sg_9[k];

        t_14[k] = f_4 * sf0_5[k]
                  - f_5 * sf1_8[k]
                  + pb_y[k] * sg_10[k];

        t_15[k] = pb_y[k] * sg_11[k];

        t_16[k] = f_0 * sf0_5[k]
                  - f_1 * sf1_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

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

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);
    const auto *sf0_7 = buffer.data(sf0 + 7);
    const auto *sf0_8 = buffer.data(sf0 + 8);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_1 = buffer.data(sf1 + 1);
    const auto *sf1_2 = buffer.data(sf1 + 2);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf0_0, sf0_3, sf0_4, sf1_0, sf1_1, sf1_2, sg_0, \
                         sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_2 * sf0_3[k]
                 - f_3 * sf1_1[k]
                 + pb_x[k] * sg_1[k];

        t_2[k] = f_2 * sf0_4[k]
                 - f_3 * sf1_2[k]
                 + pb_x[k] * sg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, t_7, pb_x, sf0_5, sf0_8, sf1_3, sf1_5, sg_3, \
                         sg_4, sg_5, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = f_4 * sf0_8[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = pb_x[k] * sg_5[k];

        t_6[k] = pb_x[k] * sg_6[k];

        t_7[k] = pb_x[k] * sg_8[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, sf0_5, sf0_7, sf0_8, sf1_3, sf1_4, sf1_5, \
                         sg_5, sg_6, sg_7, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sf0_5[k]
                 - f_1 * sf1_3[k]
                 + pb_y[k] * sg_5[k];

        t_9[k] = f_2 * sf0_7[k]
                 - f_3 * sf1_4[k]
                 + pb_y[k] * sg_6[k];

        t_10[k] = f_4 * sf0_8[k]
                  - f_5 * sf1_5[k]
                  + pb_y[k] * sg_7[k];

        t_11[k] = pb_y[k] * sg_8[k];
    }

#pragma omp simd aligned(t_12, pb_z, sf0_8, sf1_5, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sf0_8[k]
                  - f_1 * sf1_5[k]
                  + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

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

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);
    const auto *sf0_7 = buffer.data(sf0 + 7);
    const auto *sf0_8 = buffer.data(sf0 + 8);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);
    const auto *sf1_7 = buffer.data(sf1 + 7);
    const auto *sf1_8 = buffer.data(sf1 + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sf0_0, sf0_3, sf1_0, sf1_3, \
                         sg_0, sg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * sg_0[k];

        t_2[k] = pb_z[k] * sg_0[k];

        t_3[k] = f_2 * sf0_3[k]
                 - f_3 * sf1_3[k]
                 + pb_x[k] * sg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_z, sf0_4, sf0_5, sf0_8, sf1_4, sf1_5, \
                         sf1_8, sg_3, sg_4, sg_5, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sf0_4[k]
                 - f_3 * sf1_4[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_5[k];

        t_6[k] = pb_z[k] * sg_3[k];

        t_7[k] = f_4 * sf0_8[k]
                 - f_5 * sf1_8[k]
                 + pb_x[k] * sg_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, sf0_5, sf0_7, sf0_8, sf1_5, sf1_7, \
                         sf1_8, sg_7, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sf0_5[k]
                 - f_1 * sf1_5[k]
                 + pb_y[k] * sg_7[k];

        t_9[k] = pb_z[k] * sg_7[k];

        t_10[k] = f_2 * sf0_7[k]
                  - f_3 * sf1_7[k]
                  + pb_y[k] * sg_9[k];

        t_11[k] = f_4 * sf0_8[k]
                  - f_5 * sf1_8[k]
                  + pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(t_12, pb_z, sf0_8, sf1_8, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sf0_8[k]
                  - f_1 * sf1_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_1 = buffer.data(sf1 + 1);
    const auto *sf1_2 = buffer.data(sf1 + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf0_0, sf0_1, sf0_2, sf1_0, sf1_1, \
                         sf1_2, sg_0, sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_0 * sf0_1[k]
                 - f_1 * sf1_1[k]
                 + pb_y[k] * sg_1[k];

        t_2[k] = f_0 * sf0_2[k]
                 - f_1 * sf1_2[k]
                 + pb_z[k] * sg_2[k];
    }
}

auto
compute_prim_sh_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sf0, const size_t sf1, const size_t sg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_1 = buffer.data(sf0 + 1);
    const auto *sf0_2 = buffer.data(sf0 + 2);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_5 = buffer.data(sf1 + 5);
    const auto *sf1_8 = buffer.data(sf1 + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf0_0, sf0_1, sf0_2, sf1_0, sf1_5, \
                         sf1_8, sg_0, sg_5, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_0 * sf0_1[k]
                 - f_1 * sf1_5[k]
                 + pb_y[k] * sg_5[k];

        t_2[k] = f_0 * sf0_2[k]
                 - f_1 * sf1_8[k]
                 + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sf0, const size_t sf1, const size_t sg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / beta;
    const auto f_1 = 2.0 * alpha / (beta * p);
    const auto f_2 = 1.0 / beta;
    const auto f_3 = alpha / (beta * p);
    const auto f_4 = 0.5 / beta;
    const auto f_5 = 0.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf0_0 = buffer.data(sf0 + 0);
    const auto *sf0_3 = buffer.data(sf0 + 3);
    const auto *sf0_4 = buffer.data(sf0 + 4);
    const auto *sf0_5 = buffer.data(sf0 + 5);
    const auto *sf0_7 = buffer.data(sf0 + 7);
    const auto *sf0_8 = buffer.data(sf0 + 8);

    const auto *sf1_0 = buffer.data(sf1 + 0);
    const auto *sf1_1 = buffer.data(sf1 + 1);
    const auto *sf1_2 = buffer.data(sf1 + 2);
    const auto *sf1_3 = buffer.data(sf1 + 3);
    const auto *sf1_4 = buffer.data(sf1 + 4);
    const auto *sf1_5 = buffer.data(sf1 + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf0_0, sf0_3, sf0_4, sf1_0, sf1_1, sf1_2, sg_0, \
                         sg_1, sg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf0_0[k]
                 - f_1 * sf1_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_2 * sf0_3[k]
                 - f_3 * sf1_1[k]
                 + pb_x[k] * sg_1[k];

        t_2[k] = f_2 * sf0_4[k]
                 - f_3 * sf1_2[k]
                 + pb_x[k] * sg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pb_x, pb_y, sf0_5, sf0_7, sf0_8, sf1_3, sf1_4, \
                         sf1_5, sg_3, sg_4, sg_5, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sf0_5[k]
                 - f_5 * sf1_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = f_4 * sf0_8[k]
                 - f_5 * sf1_5[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = f_0 * sf0_5[k]
                 - f_1 * sf1_3[k]
                 + pb_y[k] * sg_5[k];

        t_6[k] = f_2 * sf0_7[k]
                 - f_3 * sf1_4[k]
                 + pb_y[k] * sg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, sf0_8, sf1_5, sg_7, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * sf0_8[k]
                 - f_5 * sf1_5[k]
                 + pb_y[k] * sg_7[k];

        t_8[k] = f_0 * sf0_8[k]
                 - f_1 * sf1_5[k]
                 + pb_z[k] * sg_8[k];
    }
}

}  // namespace simdt2ceri
