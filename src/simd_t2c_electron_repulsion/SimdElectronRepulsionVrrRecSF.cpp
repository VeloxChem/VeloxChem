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


#include "SimdElectronRepulsionVrrRecSF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_sf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sp0, const size_t sp1, const size_t sd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / beta;
    const auto f_1 = alpha / (beta * p);

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

    const auto *sp0_0 = buffer.data(sp0 + 0);
    const auto *sp0_1 = buffer.data(sp0 + 1);
    const auto *sp0_2 = buffer.data(sp0 + 2);

    const auto *sp1_0 = buffer.data(sp1 + 0);
    const auto *sp1_1 = buffer.data(sp1 + 1);
    const auto *sp1_2 = buffer.data(sp1 + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, sp0_0, sp0_1, sp1_0, \
                         sp1_1, sd_0, sd_1, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp0_0[k]
                 - f_1 * sp1_0[k]
                 + pb_x[k] * sd_0[k];

        t_1[k] = pb_z[k] * sd_0[k];

        t_2[k] = pb_x[k] * sd_1[k];

        t_3[k] = pb_x[k] * sd_2[k];

        t_4[k] = f_0 * sp0_1[k]
                 - f_1 * sp1_1[k]
                 + pb_y[k] * sd_1[k];

        t_5[k] = pb_z[k] * sd_1[k];
    }

#pragma omp simd aligned(t_6, t_7, pb_y, pb_z, sp0_2, sp1_2, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * sd_2[k];

        t_7[k] = f_0 * sp0_2[k]
                 - f_1 * sp1_2[k]
                 + pb_z[k] * sd_2[k];
    }
}

auto
compute_prim_sf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sp0, const size_t sp1, const size_t sd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / beta;
    const auto f_1 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sp0_0 = buffer.data(sp0 + 0);
    const auto *sp0_1 = buffer.data(sp0 + 1);
    const auto *sp0_2 = buffer.data(sp0 + 2);

    const auto *sp1_0 = buffer.data(sp1 + 0);
    const auto *sp1_1 = buffer.data(sp1 + 1);
    const auto *sp1_2 = buffer.data(sp1 + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, sp0_0, sp0_1, sp1_0, sp1_1, \
                         sd_0, sd_1, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp0_0[k]
                 - f_1 * sp1_0[k]
                 + pb_x[k] * sd_0[k];

        t_1[k] = pb_x[k] * sd_1[k];

        t_2[k] = pb_x[k] * sd_2[k];

        t_3[k] = f_0 * sp0_1[k]
                 - f_1 * sp1_1[k]
                 + pb_y[k] * sd_1[k];

        t_4[k] = pb_y[k] * sd_2[k];
    }

#pragma omp simd aligned(t_5, pb_z, sp0_2, sp1_2, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sp0_2[k]
                 - f_1 * sp1_2[k]
                 + pb_z[k] * sd_2[k];
    }
}

auto
compute_prim_sf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sp0, const size_t sp1, const size_t sd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / beta;
    const auto f_1 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sp0_0 = buffer.data(sp0 + 0);
    const auto *sp0_1 = buffer.data(sp0 + 1);
    const auto *sp0_2 = buffer.data(sp0 + 2);

    const auto *sp1_0 = buffer.data(sp1 + 0);
    const auto *sp1_1 = buffer.data(sp1 + 1);
    const auto *sp1_2 = buffer.data(sp1 + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sp0_0, sp0_1, sp0_2, sp1_0, sp1_1, \
                         sp1_2, sd_0, sd_1, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp0_0[k]
                 - f_1 * sp1_0[k]
                 + pb_x[k] * sd_0[k];

        t_1[k] = f_0 * sp0_1[k]
                 - f_1 * sp1_1[k]
                 + pb_y[k] * sd_1[k];

        t_2[k] = f_0 * sp0_2[k]
                 - f_1 * sp1_2[k]
                 + pb_z[k] * sd_2[k];
    }
}

auto
compute_prim_sf_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sp0, const size_t sp1, const size_t sd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / beta;
    const auto f_1 = alpha / (beta * p);

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

    const auto *sp0_0 = buffer.data(sp0 + 0);
    const auto *sp0_1 = buffer.data(sp0 + 1);
    const auto *sp0_2 = buffer.data(sp0 + 2);

    const auto *sp1_0 = buffer.data(sp1 + 0);
    const auto *sp1_1 = buffer.data(sp1 + 1);
    const auto *sp1_2 = buffer.data(sp1 + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, sp0_0, sp0_1, sp1_0, \
                         sp1_1, sd_0, sd_1, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp0_0[k]
                 - f_1 * sp1_0[k]
                 + pb_x[k] * sd_0[k];

        t_1[k] = pb_y[k] * sd_0[k];

        t_2[k] = pb_z[k] * sd_0[k];

        t_3[k] = pb_x[k] * sd_1[k];

        t_4[k] = pb_x[k] * sd_2[k];

        t_5[k] = f_0 * sp0_1[k]
                 - f_1 * sp1_1[k]
                 + pb_y[k] * sd_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, sp0_2, sp1_2, sd_1, \
                         sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_z[k] * sd_1[k];

        t_7[k] = pb_y[k] * sd_2[k];

        t_8[k] = f_0 * sp0_2[k]
                 - f_1 * sp1_2[k]
                 + pb_z[k] * sd_2[k];
    }
}

auto
compute_prim_sf_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sp0, const size_t sp1, const size_t sd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / beta;
    const auto f_1 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sp0_0 = buffer.data(sp0 + 0);
    const auto *sp0_1 = buffer.data(sp0 + 1);
    const auto *sp0_2 = buffer.data(sp0 + 2);

    const auto *sp1_0 = buffer.data(sp1 + 0);
    const auto *sp1_1 = buffer.data(sp1 + 1);
    const auto *sp1_2 = buffer.data(sp1 + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, sp0_0, sp0_1, sp1_0, \
                         sp1_1, sd_0, sd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sp0_0[k]
                 - f_1 * sp1_0[k]
                 + pb_x[k] * sd_0[k];

        t_1[k] = pb_y[k] * sd_0[k];

        t_2[k] = pb_z[k] * sd_0[k];

        t_3[k] = f_0 * sp0_1[k]
                 - f_1 * sp1_1[k]
                 + pb_y[k] * sd_1[k];

        t_4[k] = pb_z[k] * sd_1[k];
    }

#pragma omp simd aligned(t_5, pb_z, sp0_2, sp1_2, sd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sp0_2[k]
                 - f_1 * sp1_2[k]
                 + pb_z[k] * sd_2[k];
    }
}

}  // namespace simdt2ceri
