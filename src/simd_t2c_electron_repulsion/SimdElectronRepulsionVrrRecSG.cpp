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


#include "SimdElectronRepulsionVrrRecSG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_sg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, sd0_0, sd0_1, sd0_2, sd1_0, sd1_1, \
                         sd1_2, sf_0, sf_1, sf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_z[k] * sf_0[k];

        t_2[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_1[k];

        t_3[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_3, sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pb_x[k] * sf_3[k];

        t_5[k] = pb_x[k] * sf_5[k];

        t_6[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_3[k];

        t_7[k] = pb_z[k] * sf_3[k];

        t_8[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_4[k];

        t_9[k] = pb_y[k] * sf_5[k];
    }

#pragma omp simd aligned(t_10, pb_z, sd0_2, sd1_2, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sd0_2[k]
                  - f_1 * sd1_2[k]
                  + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, sd0_0, sd0_1, sd0_2, sd1_0, sd1_1, sd1_2, \
                         sf_0, sf_1, sf_2, sf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_1[k];

        t_2[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_2[k];

        t_3[k] = pb_x[k] * sf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_3, sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pb_x[k] * sf_5[k];

        t_5[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_3[k];

        t_6[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_4[k];

        t_7[k] = pb_y[k] * sf_5[k];

        t_8[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, sd0_0, sd0_1, sd0_2, sd1_0, sd1_1, \
                         sd1_2, sf_0, sf_1, sf_2, sf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_1[k];

        t_2[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_2[k];

        t_3[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, sd0_2, sd1_2, sf_4, \
                         sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_4[k];

        t_5[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd1_0, sd1_1, \
                         sf_0, sf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * sf_0[k];

        t_2[k] = pb_z[k] * sf_0[k];

        t_3[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_4, sf_5, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_4[k];

        t_5[k] = pb_x[k] * sf_5[k];

        t_6[k] = pb_x[k] * sf_8[k];

        t_7[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_5[k];

        t_8[k] = pb_z[k] * sf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sd0_2, sd1_2, sf_7, \
                         sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_7[k];

        t_10[k] = pb_y[k] * sf_8[k];

        t_11[k] = f_0 * sd0_2[k]
                  - f_1 * sd1_2[k]
                  + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd1_0, sd1_1, \
                         sf_0, sf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * sf_0[k];

        t_2[k] = pb_z[k] * sf_0[k];

        t_3[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_2, sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_2[k];

        t_5[k] = pb_x[k] * sf_3[k];

        t_6[k] = pb_x[k] * sf_5[k];

        t_7[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_3[k];

        t_8[k] = pb_z[k] * sf_3[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sd0_2, sd1_2, sf_4, \
                         sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_4[k];

        t_10[k] = pb_y[k] * sf_5[k];

        t_11[k] = f_0 * sd0_2[k]
                  - f_1 * sd1_2[k]
                  + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd1_0, sd1_1, \
                         sf_0, sf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * sf_0[k];

        t_2[k] = pb_z[k] * sf_0[k];

        t_3[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_4, sf_5, sf_7, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_4[k];

        t_5[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_5[k];

        t_6[k] = pb_z[k] * sf_5[k];

        t_7[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_7[k];

        t_8[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd1_0, sd1_1, \
                         sf_0, sf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * sf_0[k];

        t_2[k] = pb_z[k] * sf_0[k];

        t_3[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_3, sf_4, sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_3[k];

        t_5[k] = pb_x[k] * sf_4[k];

        t_6[k] = pb_x[k] * sf_7[k];

        t_7[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_4[k];

        t_8[k] = pb_z[k] * sf_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sd0_2, sd1_2, sf_6, \
                         sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_6[k];

        t_10[k] = pb_y[k] * sf_7[k];

        t_11[k] = f_0 * sd0_2[k]
                  - f_1 * sd1_2[k]
                  + pb_z[k] * sf_7[k];
    }
}

auto
compute_prim_sg_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd0_2, sd1_0, sd1_1, \
                         sd1_2, sf_0, sf_1, sf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_1[k];

        t_2[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_2[k];
    }
}

auto
compute_prim_sg_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd0_2, sd1_0, sd1_1, \
                         sd1_2, sf_0, sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_3[k];

        t_2[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sd0, const size_t sd1, const size_t sf,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd0_2, sd1_0, sd1_1, \
                         sd1_2, sf_0, sf_5, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_5[k];

        t_2[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sd0, const size_t sd1, const size_t sf,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / beta;
    const auto f_1 = 1.5 * alpha / (beta * p);
    const auto f_2 = 0.5 / beta;
    const auto f_3 = 0.5 * alpha / (beta * p);

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

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_1 = buffer.data(sd0 + 1);
    const auto *sd0_2 = buffer.data(sd0 + 2);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_1 = buffer.data(sd1 + 1);
    const auto *sd1_2 = buffer.data(sd1 + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sd0_0, sd0_1, sd1_0, sd1_1, \
                         sf_0, sf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * sf_0[k];

        t_2[k] = pb_z[k] * sf_0[k];

        t_3[k] = f_2 * sd0_1[k]
                 - f_3 * sd1_1[k]
                 + pb_x[k] * sf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sd0_1, sd0_2, sd1_1, \
                         sd1_2, sf_2, sf_3, sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_x[k] * sf_2[k];

        t_5[k] = f_0 * sd0_1[k]
                 - f_1 * sd1_1[k]
                 + pb_y[k] * sf_3[k];

        t_6[k] = pb_z[k] * sf_3[k];

        t_7[k] = f_2 * sd0_2[k]
                 - f_3 * sd1_2[k]
                 + pb_y[k] * sf_4[k];

        t_8[k] = f_0 * sd0_2[k]
                 - f_1 * sd1_2[k]
                 + pb_z[k] * sf_5[k];
    }
}

}  // namespace simdt2ceri
