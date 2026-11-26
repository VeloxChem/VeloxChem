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
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd0_0 = buffer.data(sd0 + 0);
    const auto *sd0_3 = buffer.data(sd0 + 3);
    const auto *sd0_5 = buffer.data(sd0 + 5);

    const auto *sd1_0 = buffer.data(sd1 + 0);
    const auto *sd1_3 = buffer.data(sd1 + 3);
    const auto *sd1_5 = buffer.data(sd1 + 5);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, sd0_0, sd0_3, sd1_0, \
                         sd1_3, sf_0, sf_2, sf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd0_0[k]
                 - f_1 * sd1_0[k]
                 + pb_x[k] * sf_0[k];

        t_1[k] = pb_y[k] * sf_0[k];

        t_2[k] = pb_z[k] * sf_0[k];

        t_3[k] = f_2 * sd0_3[k]
                 - f_3 * sd1_3[k]
                 + pb_x[k] * sf_3[k];

        t_4[k] = pb_y[k] * sf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, sd0_5, sd1_5, sf_5, sf_6, sf_7, sf_8, \
                         sf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * sd0_5[k]
                 - f_3 * sd1_5[k]
                 + pb_x[k] * sf_5[k];

        t_6[k] = pb_x[k] * sf_6[k];

        t_7[k] = pb_x[k] * sf_7[k];

        t_8[k] = pb_x[k] * sf_8[k];

        t_9[k] = pb_x[k] * sf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pb_z, sd0_3, sd0_5, sd1_3, sd1_5, \
                         sf_6, sf_8, sf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sd0_3[k]
                  - f_1 * sd1_3[k]
                  + pb_y[k] * sf_6[k];

        t_11[k] = pb_z[k] * sf_6[k];

        t_12[k] = f_2 * sd0_5[k]
                  - f_3 * sd1_5[k]
                  + pb_y[k] * sf_8[k];

        t_13[k] = pb_y[k] * sf_9[k];

        t_14[k] = f_0 * sd0_5[k]
                  - f_1 * sd1_5[k]
                  + pb_z[k] * sf_9[k];
    }
}

}  // namespace simdt2ceri
