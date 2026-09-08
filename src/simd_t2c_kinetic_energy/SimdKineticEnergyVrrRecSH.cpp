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


#include "SimdKineticEnergyVrrRecSH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_sh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_2 = buffer.data(sf_s + 2);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);
    const auto *sh_s_13 = buffer.data(sh_s + 13);
    const auto *sh_s_14 = buffer.data(sh_s + 14);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, sf_s_0, sf_s_1, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sf_1, sg_0, sg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_1 * sh_s_1[k]
                 + pb_z[k] * sg_0[k];

        t_2[k] = -f_3 * sf_s_1[k]
                 + f_1 * sh_s_2[k]
                 + f_4 * sf_1[k]
                 + pb_x[k] * sg_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sf_s_2, sf_s_3, sf_s_5, sh_s_3, sh_s_4, sh_s_5, \
                         sf_2, sf_3, sf_5, sg_2, sg_3, sg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_2[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_2[k]
                 + pb_x[k] * sg_2[k];

        t_4[k] = -f_5 * sf_s_3[k]
                 + f_1 * sh_s_4[k]
                 + f_6 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_5[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_y, sf_s_3, sh_s_6, sh_s_7, sh_s_8, \
                         sh_s_9, sf_3, sg_5, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_x[k] * sg_5[k];

        t_7[k] = f_1 * sh_s_7[k]
                 + pb_x[k] * sg_6[k];

        t_8[k] = f_1 * sh_s_8[k]
                 + pb_x[k] * sg_8[k];

        t_9[k] = -f_0 * sf_s_3[k]
                 + f_1 * sh_s_9[k]
                 + f_2 * sf_3[k]
                 + pb_y[k] * sg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, sf_s_4, sf_s_5, sh_s_10, sh_s_11, \
                         sh_s_12, sf_4, sf_5, sg_5, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * sh_s_10[k]
                  + pb_z[k] * sg_5[k];

        t_11[k] = -f_3 * sf_s_4[k]
                  + f_1 * sh_s_11[k]
                  + f_4 * sf_4[k]
                  + pb_y[k] * sg_6[k];

        t_12[k] = -f_5 * sf_s_5[k]
                  + f_1 * sh_s_12[k]
                  + f_6 * sf_5[k]
                  + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_y, pb_z, sf_s_5, sh_s_13, sh_s_14, sf_5, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * sh_s_13[k]
                  + pb_y[k] * sg_8[k];

        t_14[k] = -f_0 * sf_s_5[k]
                  + f_1 * sh_s_14[k]
                  + f_2 * sf_5[k]
                  + pb_z[k] * sg_8[k];
    }
}

auto
compute_prim_sh_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);
    const auto *sh_s_13 = buffer.data(sh_s + 13);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf_s_0, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_1 * sh_s_1[k]
                 + pb_y[k] * sg_0[k];

        t_2[k] = f_1 * sh_s_2[k]
                 + pb_z[k] * sg_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sf_s_3, sf_s_4, sf_s_5, sh_s_3, sh_s_4, sh_s_5, \
                         sf_3, sf_4, sf_5, sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = -f_3 * sf_s_4[k]
                 + f_1 * sh_s_4[k]
                 + f_4 * sf_4[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sf_s_8, sh_s_6, sh_s_7, sh_s_8, \
                         sf_8, sg_3, sg_4, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_z[k] * sg_3[k];

        t_7[k] = f_1 * sh_s_7[k]
                 + pb_y[k] * sg_4[k];

        t_8[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_8[k]
                 + f_6 * sf_8[k]
                 + pb_x[k] * sg_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sf_s_5, sf_s_7, sh_s_9, sh_s_10, \
                         sh_s_11, sf_5, sf_7, sg_7, sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_0 * sf_s_5[k]
                 + f_1 * sh_s_9[k]
                 + f_2 * sf_5[k]
                 + pb_y[k] * sg_7[k];

        t_10[k] = f_1 * sh_s_10[k]
                  + pb_z[k] * sg_7[k];

        t_11[k] = -f_3 * sf_s_7[k]
                  + f_1 * sh_s_11[k]
                  + f_4 * sf_7[k]
                  + pb_y[k] * sg_9[k];
    }

#pragma omp simd aligned(t_12, t_13, pb_y, pb_z, sf_s_8, sh_s_12, sh_s_13, sf_8, sg_10, \
                         sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_5 * sf_s_8[k]
                  + f_1 * sh_s_12[k]
                  + f_6 * sf_8[k]
                  + pb_y[k] * sg_10[k];

        t_13[k] = -f_0 * sf_s_8[k]
                  + f_1 * sh_s_13[k]
                  + f_2 * sf_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_2 = buffer.data(sf_s + 2);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);
    const auto *sh_s_13 = buffer.data(sh_s + 13);
    const auto *sh_s_14 = buffer.data(sh_s + 14);
    const auto *sh_s_15 = buffer.data(sh_s + 15);
    const auto *sh_s_16 = buffer.data(sh_s + 16);
    const auto *sh_s_17 = buffer.data(sh_s + 17);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf_s_0, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_1 * sh_s_1[k]
                 + pb_y[k] * sg_0[k];

        t_2[k] = f_1 * sh_s_2[k]
                 + pb_z[k] * sg_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sf_s_1, sf_s_2, sf_s_3, sh_s_3, sh_s_4, sh_s_5, \
                         sf_1, sf_2, sf_3, sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_1[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_1[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = -f_3 * sf_s_2[k]
                 + f_1 * sh_s_4[k]
                 + f_4 * sf_2[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = -f_5 * sf_s_3[k]
                 + f_1 * sh_s_5[k]
                 + f_6 * sf_3[k]
                 + pb_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sf_s_5, sh_s_6, sh_s_7, sh_s_8, \
                         sf_5, sg_3, sg_4, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_z[k] * sg_3[k];

        t_7[k] = f_1 * sh_s_7[k]
                 + pb_y[k] * sg_4[k];

        t_8[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_8[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, pb_y, sf_s_3, sh_s_9, sh_s_10, sh_s_11, \
                         sh_s_12, sf_3, sg_7, sg_9, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * sh_s_9[k]
                 + pb_x[k] * sg_7[k];

        t_10[k] = f_1 * sh_s_10[k]
                  + pb_x[k] * sg_9[k];

        t_11[k] = f_1 * sh_s_11[k]
                  + pb_x[k] * sg_11[k];

        t_12[k] = -f_0 * sf_s_3[k]
                  + f_1 * sh_s_12[k]
                  + f_2 * sf_3[k]
                  + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, pb_z, sf_s_4, sf_s_5, sh_s_13, sh_s_14, \
                         sh_s_15, sf_4, sf_5, sg_7, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * sh_s_13[k]
                  + pb_z[k] * sg_7[k];

        t_14[k] = -f_3 * sf_s_4[k]
                  + f_1 * sh_s_14[k]
                  + f_4 * sf_4[k]
                  + pb_y[k] * sg_9[k];

        t_15[k] = -f_5 * sf_s_5[k]
                  + f_1 * sh_s_15[k]
                  + f_6 * sf_5[k]
                  + pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(t_16, t_17, pb_y, pb_z, sf_s_5, sh_s_16, sh_s_17, sf_5, \
                         sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * sh_s_16[k]
                  + pb_y[k] * sg_11[k];

        t_17[k] = -f_0 * sf_s_5[k]
                  + f_1 * sh_s_17[k]
                  + f_2 * sf_5[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf_s_0, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_1 * sh_s_1[k]
                 + pb_y[k] * sg_0[k];

        t_2[k] = f_1 * sh_s_2[k]
                 + pb_z[k] * sg_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sf_s_3, sf_s_4, sf_s_5, sh_s_3, sh_s_4, sh_s_5, \
                         sf_3, sf_4, sf_5, sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = -f_3 * sf_s_4[k]
                 + f_1 * sh_s_4[k]
                 + f_4 * sf_4[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sf_s_5, sf_s_8, sh_s_6, sh_s_7, \
                         sh_s_8, sf_5, sf_8, sg_3, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_z[k] * sg_3[k];

        t_7[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_7[k]
                 + f_6 * sf_8[k]
                 + pb_x[k] * sg_6[k];

        t_8[k] = -f_0 * sf_s_5[k]
                 + f_1 * sh_s_8[k]
                 + f_2 * sf_5[k]
                 + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sf_s_7, sf_s_8, sh_s_9, sh_s_10, \
                         sh_s_11, sf_7, sf_8, sg_7, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * sh_s_9[k]
                 + pb_z[k] * sg_7[k];

        t_10[k] = -f_3 * sf_s_7[k]
                  + f_1 * sh_s_10[k]
                  + f_4 * sf_7[k]
                  + pb_y[k] * sg_9[k];

        t_11[k] = -f_5 * sf_s_8[k]
                  + f_1 * sh_s_11[k]
                  + f_6 * sf_8[k]
                  + pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(t_12, pb_z, sf_s_8, sh_s_12, sf_8, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_0 * sf_s_8[k]
                  + f_1 * sh_s_12[k]
                  + f_2 * sf_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);
    const auto *sh_s_13 = buffer.data(sh_s + 13);
    const auto *sh_s_14 = buffer.data(sh_s + 14);
    const auto *sh_s_15 = buffer.data(sh_s + 15);
    const auto *sh_s_16 = buffer.data(sh_s + 16);
    const auto *sh_s_17 = buffer.data(sh_s + 17);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf_s_0, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_1 * sh_s_1[k]
                 + pb_y[k] * sg_0[k];

        t_2[k] = f_1 * sh_s_2[k]
                 + pb_z[k] * sg_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sf_s_3, sf_s_4, sf_s_5, sh_s_3, sh_s_4, sh_s_5, \
                         sf_3, sf_4, sf_5, sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = -f_3 * sf_s_4[k]
                 + f_1 * sh_s_4[k]
                 + f_4 * sf_4[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sf_s_8, sh_s_6, sh_s_7, sh_s_8, \
                         sf_8, sg_3, sg_4, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_z[k] * sg_3[k];

        t_7[k] = f_1 * sh_s_7[k]
                 + pb_y[k] * sg_4[k];

        t_8[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_8[k]
                 + f_6 * sf_8[k]
                 + pb_x[k] * sg_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, pb_y, sf_s_5, sh_s_9, sh_s_10, sh_s_11, \
                         sh_s_12, sf_5, sg_7, sg_9, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * sh_s_9[k]
                 + pb_x[k] * sg_7[k];

        t_10[k] = f_1 * sh_s_10[k]
                  + pb_x[k] * sg_9[k];

        t_11[k] = f_1 * sh_s_11[k]
                  + pb_x[k] * sg_11[k];

        t_12[k] = -f_0 * sf_s_5[k]
                  + f_1 * sh_s_12[k]
                  + f_2 * sf_5[k]
                  + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, pb_z, sf_s_7, sf_s_8, sh_s_13, sh_s_14, \
                         sh_s_15, sf_7, sf_8, sg_7, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * sh_s_13[k]
                  + pb_z[k] * sg_7[k];

        t_14[k] = -f_3 * sf_s_7[k]
                  + f_1 * sh_s_14[k]
                  + f_4 * sf_7[k]
                  + pb_y[k] * sg_9[k];

        t_15[k] = -f_5 * sf_s_8[k]
                  + f_1 * sh_s_15[k]
                  + f_6 * sf_8[k]
                  + pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(t_16, t_17, pb_y, pb_z, sf_s_8, sh_s_16, sh_s_17, sf_8, \
                         sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * sh_s_16[k]
                  + pb_y[k] * sg_11[k];

        t_17[k] = -f_0 * sf_s_8[k]
                  + f_1 * sh_s_17[k]
                  + f_2 * sf_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf_s_0, sf_s_3, sf_s_4, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sf_3, sf_4, sg_0, sg_3, sg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_1[k]
                 + f_4 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_2[k] = -f_3 * sf_s_4[k]
                 + f_1 * sh_s_2[k]
                 + f_4 * sf_4[k]
                 + pb_x[k] * sg_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, sf_s_5, sf_s_8, sh_s_3, sh_s_4, sh_s_5, \
                         sf_5, sf_8, sg_5, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_3[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_5[k];

        t_4[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_4[k]
                 + f_6 * sf_8[k]
                 + pb_x[k] * sg_6[k];

        t_5[k] = -f_0 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_2 * sf_5[k]
                 + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, sf_s_7, sf_s_8, sh_s_6, sh_s_7, sh_s_8, \
                         sf_7, sf_8, sg_9, sg_10, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_3 * sf_s_7[k]
                 + f_1 * sh_s_6[k]
                 + f_4 * sf_7[k]
                 + pb_y[k] * sg_9[k];

        t_7[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_7[k]
                 + f_6 * sf_8[k]
                 + pb_y[k] * sg_10[k];

        t_8[k] = -f_0 * sf_s_8[k]
                 + f_1 * sh_s_8[k]
                 + f_2 * sf_8[k]
                 + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);
    const auto *sh_s_13 = buffer.data(sh_s + 13);
    const auto *sh_s_14 = buffer.data(sh_s + 14);
    const auto *sh_s_15 = buffer.data(sh_s + 15);
    const auto *sh_s_16 = buffer.data(sh_s + 16);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf_s_0, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = f_1 * sh_s_1[k]
                 + pb_y[k] * sg_0[k];

        t_2[k] = f_1 * sh_s_2[k]
                 + pb_z[k] * sg_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sf_s_3, sf_s_4, sf_s_5, sh_s_3, sh_s_4, sh_s_5, \
                         sf_3, sf_4, sf_5, sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = -f_3 * sf_s_4[k]
                 + f_1 * sh_s_4[k]
                 + f_4 * sf_4[k]
                 + pb_x[k] * sg_4[k];

        t_5[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_z, sf_s_8, sh_s_6, sh_s_7, sh_s_8, \
                         sh_s_9, sf_8, sg_3, sg_6, sg_7, sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_z[k] * sg_3[k];

        t_7[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_7[k]
                 + f_6 * sf_8[k]
                 + pb_x[k] * sg_6[k];

        t_8[k] = f_1 * sh_s_8[k]
                 + pb_x[k] * sg_7[k];

        t_9[k] = f_1 * sh_s_9[k]
                 + pb_x[k] * sg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, pb_y, pb_z, sf_s_5, sh_s_10, sh_s_11, \
                         sh_s_12, sf_5, sg_7, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * sh_s_10[k]
                  + pb_x[k] * sg_11[k];

        t_11[k] = -f_0 * sf_s_5[k]
                  + f_1 * sh_s_11[k]
                  + f_2 * sf_5[k]
                  + pb_y[k] * sg_7[k];

        t_12[k] = f_1 * sh_s_12[k]
                  + pb_z[k] * sg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, sf_s_7, sf_s_8, sh_s_13, sh_s_14, sh_s_15, \
                         sf_7, sf_8, sg_9, sg_10, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_3 * sf_s_7[k]
                  + f_1 * sh_s_13[k]
                  + f_4 * sf_7[k]
                  + pb_y[k] * sg_9[k];

        t_14[k] = -f_5 * sf_s_8[k]
                  + f_1 * sh_s_14[k]
                  + f_6 * sf_8[k]
                  + pb_y[k] * sg_10[k];

        t_15[k] = f_1 * sh_s_15[k]
                  + pb_y[k] * sg_11[k];
    }

#pragma omp simd aligned(t_16, pb_z, sf_s_8, sh_s_16, sf_8, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_0 * sf_s_8[k]
                  + f_1 * sh_s_16[k]
                  + f_2 * sf_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);

    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);

#pragma omp simd aligned(t_0, t_1, t_2, pb_y, sf_s_5, sf_s_7, sf_s_8, sh_s_0, sh_s_1, sh_s_2, \
                         sf_5, sf_7, sf_8, sg_3, sg_5, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_5[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_5[k]
                 + pb_y[k] * sg_3[k];

        t_1[k] = -f_3 * sf_s_7[k]
                 + f_1 * sh_s_1[k]
                 + f_4 * sf_7[k]
                 + pb_y[k] * sg_5[k];

        t_2[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_2[k]
                 + f_6 * sf_8[k]
                 + pb_y[k] * sg_6[k];
    }

#pragma omp simd aligned(t_3, pb_z, sf_s_8, sh_s_3, sf_8, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_0 * sf_s_8[k]
                 + f_1 * sh_s_3[k]
                 + f_2 * sf_8[k]
                 + pb_z[k] * sg_7[k];
    }
}

auto
compute_prim_sh_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sf_s_0, sf_s_3, sf_s_4, sh_s_0, sh_s_1, sh_s_2, \
                         sf_0, sf_3, sf_4, sg_0, sg_3, sg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_1[k]
                 + f_4 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_2[k] = -f_3 * sf_s_4[k]
                 + f_1 * sh_s_2[k]
                 + f_4 * sf_4[k]
                 + pb_x[k] * sg_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, sf_s_5, sf_s_8, sh_s_3, sh_s_4, sh_s_5, \
                         sf_5, sf_8, sg_5, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_5 * sf_s_5[k]
                 + f_1 * sh_s_3[k]
                 + f_6 * sf_5[k]
                 + pb_x[k] * sg_5[k];

        t_4[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_4[k]
                 + f_6 * sf_8[k]
                 + pb_x[k] * sg_6[k];

        t_5[k] = -f_0 * sf_s_5[k]
                 + f_1 * sh_s_5[k]
                 + f_2 * sf_5[k]
                 + pb_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, sf_s_7, sf_s_8, sh_s_6, sh_s_7, sh_s_8, \
                         sf_7, sf_8, sg_7, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sh_s_6[k]
                 + pb_z[k] * sg_7[k];

        t_7[k] = -f_3 * sf_s_7[k]
                 + f_1 * sh_s_7[k]
                 + f_4 * sf_7[k]
                 + pb_y[k] * sg_9[k];

        t_8[k] = -f_5 * sf_s_8[k]
                 + f_1 * sh_s_8[k]
                 + f_6 * sf_8[k]
                 + pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_y, pb_z, sf_s_8, sh_s_9, sh_s_10, sf_8, \
                         sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * sh_s_9[k]
                 + pb_y[k] * sg_11[k];

        t_10[k] = -f_0 * sf_s_8[k]
                  + f_1 * sh_s_10[k]
                  + f_2 * sf_8[k]
                  + pb_z[k] * sg_11[k];
    }
}

auto
compute_prim_sh_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sf_s, const size_t sh_s, const size_t sf,
                                 const size_t sg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);

    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);

    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_y, sf_s_1, sf_s_3, sf_s_4, sh_s_0, sh_s_1, sh_s_2, \
                         sf_1, sf_3, sf_4, sg_1, sg_3, sg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_1[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_1[k]
                 + pb_y[k] * sg_1[k];

        t_1[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_1[k]
                 + f_4 * sf_3[k]
                 + pb_y[k] * sg_3[k];

        t_2[k] = -f_5 * sf_s_4[k]
                 + f_1 * sh_s_2[k]
                 + f_6 * sf_4[k]
                 + pb_y[k] * sg_4[k];
    }

#pragma omp simd aligned(t_3, pb_z, sf_s_4, sh_s_3, sf_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_0 * sf_s_4[k]
                 + f_1 * sh_s_3[k]
                 + f_2 * sf_4[k]
                 + pb_z[k] * sg_5[k];
    }
}

auto
compute_prim_sh_kinetic_energy_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                  const size_t sf_s, const size_t sh_s, const size_t sf,
                                  const size_t sg, const size_t ncols, const double alpha,
                                  const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 2.0 * alpha / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = alpha / p;
    const auto f_6 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sf_s_0, sf_s_1, sh_s_0, sh_s_1, \
                         sh_s_2, sf_0, sf_1, sg_0, sg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sf_s_0[k]
                 + f_1 * sh_s_0[k]
                 + f_2 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = -f_0 * sf_s_1[k]
                 + f_1 * sh_s_1[k]
                 + f_2 * sf_1[k]
                 + pb_y[k] * sg_1[k];

        t_2[k] = f_1 * sh_s_2[k]
                 + pb_z[k] * sg_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, sf_s_3, sf_s_4, sh_s_3, sh_s_4, sh_s_5, sf_3, \
                         sf_4, sg_3, sg_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sf_s_3[k]
                 + f_1 * sh_s_3[k]
                 + f_4 * sf_3[k]
                 + pb_y[k] * sg_3[k];

        t_4[k] = -f_5 * sf_s_4[k]
                 + f_1 * sh_s_4[k]
                 + f_6 * sf_4[k]
                 + pb_y[k] * sg_4[k];

        t_5[k] = f_1 * sh_s_5[k]
                 + pb_y[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, pb_z, sf_s_4, sh_s_6, sf_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_0 * sf_s_4[k]
                 + f_1 * sh_s_6[k]
                 + f_2 * sf_4[k]
                 + pb_z[k] * sg_5[k];
    }
}

}  // namespace simdkin
