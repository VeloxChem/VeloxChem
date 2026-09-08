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

}  // namespace simdt2ceri
