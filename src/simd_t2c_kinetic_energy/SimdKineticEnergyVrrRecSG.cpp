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


#include "SimdKineticEnergyVrrRecSG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_sg_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

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

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);
    const auto *sg_s_9 = buffer.data(sg_s + 9);
    const auto *sg_s_10 = buffer.data(sg_s + 10);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, sd_s_0, sd_s_1, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sd_1, sf_0, sf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_1 * sg_s_1[k]
                 + pb_z[k] * sf_0[k];

        t_2[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_2[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sd_s_2, sg_s_3, sg_s_4, sg_s_5, sd_2, sf_2, \
                         sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_3[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_2[k];

        t_4[k] = f_1 * sg_s_4[k]
                 + pb_x[k] * sf_3[k];

        t_5[k] = f_1 * sg_s_5[k]
                 + pb_x[k] * sf_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, sd_s_1, sd_s_2, sg_s_6, sg_s_7, sg_s_8, \
                         sd_1, sd_2, sf_3, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_6[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_3[k];

        t_7[k] = f_1 * sg_s_7[k]
                 + pb_z[k] * sf_3[k];

        t_8[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_8[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_4[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_y, pb_z, sd_s_2, sg_s_9, sg_s_10, sd_2, \
                         sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * sg_s_9[k]
                 + pb_y[k] * sf_5[k];

        t_10[k] = -f_0 * sd_s_2[k]
                  + f_1 * sg_s_10[k]
                  + f_2 * sd_2[k]
                  + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

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

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sd_s_0, sd_s_1, sd_s_2, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sd_1, sd_2, sf_0, sf_1, sf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_1[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_1[k];

        t_2[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_2[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, sd_s_1, sg_s_3, sg_s_4, sg_s_5, sd_1, \
                         sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * sg_s_3[k]
                 + pb_x[k] * sf_3[k];

        t_4[k] = f_1 * sg_s_4[k]
                 + pb_x[k] * sf_5[k];

        t_5[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_5[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, sd_s_2, sg_s_6, sg_s_7, sg_s_8, sd_2, \
                         sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_6[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_4[k];

        t_7[k] = f_1 * sg_s_7[k]
                 + pb_y[k] * sf_5[k];

        t_8[k] = -f_0 * sd_s_2[k]
                 + f_1 * sg_s_8[k]
                 + f_2 * sd_2[k]
                 + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

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

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd_s_0, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_1 * sg_s_1[k]
                 + pb_y[k] * sf_0[k];

        t_2[k] = f_1 * sg_s_2[k]
                 + pb_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, sd_s_1, sd_s_2, sg_s_3, sg_s_4, sg_s_5, \
                         sd_1, sd_2, sf_3, sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_3[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_3[k];

        t_4[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_4[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_4[k];

        t_5[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_5[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, sd_s_2, sg_s_6, sg_s_7, sg_s_8, sd_2, \
                         sf_5, sf_7, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sg_s_6[k]
                 + pb_z[k] * sf_5[k];

        t_7[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_7[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_7[k];

        t_8[k] = -f_0 * sd_s_2[k]
                 + f_1 * sg_s_8[k]
                 + f_2 * sd_2[k]
                 + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

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

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);
    const auto *sg_s_9 = buffer.data(sg_s + 9);
    const auto *sg_s_10 = buffer.data(sg_s + 10);
    const auto *sg_s_11 = buffer.data(sg_s + 11);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd_s_0, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_1 * sg_s_1[k]
                 + pb_y[k] * sf_0[k];

        t_2[k] = f_1 * sg_s_2[k]
                 + pb_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sd_s_1, sd_s_2, sg_s_3, sg_s_4, sg_s_5, sd_1, \
                         sd_2, sf_3, sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_3[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_3[k];

        t_4[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_4[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_4[k];

        t_5[k] = f_1 * sg_s_5[k]
                 + pb_x[k] * sf_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sd_s_1, sg_s_6, sg_s_7, sg_s_8, \
                         sd_1, sf_5, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sg_s_6[k]
                 + pb_x[k] * sf_8[k];

        t_7[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_7[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_5[k];

        t_8[k] = f_1 * sg_s_8[k]
                 + pb_z[k] * sf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sd_s_2, sg_s_9, sg_s_10, sg_s_11, sd_2, \
                         sf_7, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_9[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_7[k];

        t_10[k] = f_1 * sg_s_10[k]
                  + pb_y[k] * sf_8[k];

        t_11[k] = -f_0 * sd_s_2[k]
                  + f_1 * sg_s_11[k]
                  + f_2 * sd_2[k]
                  + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

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

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);
    const auto *sg_s_9 = buffer.data(sg_s + 9);
    const auto *sg_s_10 = buffer.data(sg_s + 10);
    const auto *sg_s_11 = buffer.data(sg_s + 11);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd_s_0, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = f_1 * sg_s_1[k]
                 + pb_y[k] * sf_0[k];

        t_2[k] = f_1 * sg_s_2[k]
                 + pb_z[k] * sf_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sd_s_1, sd_s_2, sg_s_3, sg_s_4, sg_s_5, sd_1, \
                         sd_2, sf_1, sf_2, sf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_3[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_1[k];

        t_4[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_4[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_2[k];

        t_5[k] = f_1 * sg_s_5[k]
                 + pb_x[k] * sf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sd_s_1, sg_s_6, sg_s_7, sg_s_8, \
                         sd_1, sf_3, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sg_s_6[k]
                 + pb_x[k] * sf_5[k];

        t_7[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_7[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_3[k];

        t_8[k] = f_1 * sg_s_8[k]
                 + pb_z[k] * sf_3[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, sd_s_2, sg_s_9, sg_s_10, sg_s_11, sd_2, \
                         sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_9[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_4[k];

        t_10[k] = f_1 * sg_s_10[k]
                  + pb_y[k] * sf_5[k];

        t_11[k] = -f_0 * sd_s_2[k]
                  + f_1 * sg_s_11[k]
                  + f_2 * sd_2[k]
                  + pb_z[k] * sf_5[k];
    }
}

auto
compute_prim_sg_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sd_s_0, sd_s_1, sd_s_2, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sd_1, sd_2, sf_0, sf_3, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_1[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_3[k];

        t_2[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_2[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, sd_s_1, sd_s_2, sg_s_3, sg_s_4, sg_s_5, \
                         sd_1, sd_2, sf_5, sf_7, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_3[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_5[k];

        t_4[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_4[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_7[k];

        t_5[k] = -f_0 * sd_s_2[k]
                 + f_1 * sg_s_5[k]
                 + f_2 * sd_2[k]
                 + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);

    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);

#pragma omp simd aligned(t_0, t_1, t_2, pb_y, pb_z, sd_s_1, sd_s_2, sg_s_0, sg_s_1, sg_s_2, \
                         sd_1, sd_2, sf_1, sf_3, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_1[k];

        t_1[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_1[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_3[k];

        t_2[k] = -f_0 * sd_s_2[k]
                 + f_1 * sg_s_2[k]
                 + f_2 * sd_2[k]
                 + pb_z[k] * sf_4[k];
    }
}

auto
compute_prim_sg_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sd_s_0, sd_s_1, sd_s_2, sg_s_0, sg_s_1, sg_s_2, \
                         sd_0, sd_1, sd_2, sf_0, sf_3, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = -f_3 * sd_s_1[k]
                 + f_1 * sg_s_1[k]
                 + f_4 * sd_1[k]
                 + pb_x[k] * sf_3[k];

        t_2[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_2[k]
                 + f_4 * sd_2[k]
                 + pb_x[k] * sf_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, sd_s_1, sd_s_2, sg_s_3, sg_s_4, sg_s_5, \
                         sd_1, sd_2, sf_5, sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_3[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_5[k];

        t_4[k] = f_1 * sg_s_4[k]
                 + pb_z[k] * sf_5[k];

        t_5[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_5[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_7[k];
    }

#pragma omp simd aligned(t_6, t_7, pb_y, pb_z, sd_s_2, sg_s_6, sg_s_7, sd_2, \
                         sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sg_s_6[k]
                 + pb_y[k] * sf_8[k];

        t_7[k] = -f_0 * sd_s_2[k]
                 + f_1 * sg_s_7[k]
                 + f_2 * sd_2[k]
                 + pb_z[k] * sf_8[k];
    }
}

auto
compute_prim_sg_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sd_s, const size_t sg_s, const size_t sd,
                                 const size_t sf, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sd_s_0, sd_s_1, sg_s_0, sg_s_1, \
                         sg_s_2, sd_0, sd_1, sf_0, sf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sd_s_0[k]
                 + f_1 * sg_s_0[k]
                 + f_2 * sd_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = -f_0 * sd_s_1[k]
                 + f_1 * sg_s_1[k]
                 + f_2 * sd_1[k]
                 + pb_y[k] * sf_1[k];

        t_2[k] = f_1 * sg_s_2[k]
                 + pb_z[k] * sf_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, sd_s_2, sg_s_3, sg_s_4, sg_s_5, sd_2, \
                         sf_3, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sd_s_2[k]
                 + f_1 * sg_s_3[k]
                 + f_4 * sd_2[k]
                 + pb_y[k] * sf_3[k];

        t_4[k] = f_1 * sg_s_4[k]
                 + pb_y[k] * sf_4[k];

        t_5[k] = -f_0 * sd_s_2[k]
                 + f_1 * sg_s_5[k]
                 + f_2 * sd_2[k]
                 + pb_z[k] * sf_4[k];
    }
}

}  // namespace simdkin
