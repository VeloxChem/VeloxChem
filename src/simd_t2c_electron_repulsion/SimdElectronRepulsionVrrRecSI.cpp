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
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg0_0 = buffer.data(sg0 + 0);
    const auto *sg0_3 = buffer.data(sg0 + 3);
    const auto *sg0_5 = buffer.data(sg0 + 5);
    const auto *sg0_6 = buffer.data(sg0 + 6);
    const auto *sg0_9 = buffer.data(sg0 + 9);
    const auto *sg0_10 = buffer.data(sg0 + 10);
    const auto *sg0_12 = buffer.data(sg0 + 12);
    const auto *sg0_13 = buffer.data(sg0 + 13);
    const auto *sg0_14 = buffer.data(sg0 + 14);

    const auto *sg1_0 = buffer.data(sg1 + 0);
    const auto *sg1_3 = buffer.data(sg1 + 3);
    const auto *sg1_5 = buffer.data(sg1 + 5);
    const auto *sg1_6 = buffer.data(sg1 + 6);
    const auto *sg1_9 = buffer.data(sg1 + 9);
    const auto *sg1_10 = buffer.data(sg1 + 10);
    const auto *sg1_12 = buffer.data(sg1 + 12);
    const auto *sg1_13 = buffer.data(sg1 + 13);
    const auto *sg1_14 = buffer.data(sg1 + 14);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, sg0_0, sg0_3, sg1_0, \
                         sg1_3, sh_0, sh_2, sh_3 : simd::cache_line_size())
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

        t_4[k] = pb_y[k] * sh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, sg0_5, sg0_6, sg1_5, sg1_6, \
                         sh_3, sh_5, sh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * sg0_5[k]
                 - f_3 * sg1_5[k]
                 + pb_x[k] * sh_5[k];

        t_6[k] = f_4 * sg0_6[k]
                 - f_5 * sg1_6[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = pb_z[k] * sh_3[k];

        t_8[k] = pb_y[k] * sh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, pb_z, sg0_9, sg0_10, sg0_12, sg1_9, \
                         sg1_10, sg1_12, sh_6, sh_9, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * sg0_9[k]
                 - f_5 * sg1_9[k]
                 + pb_x[k] * sh_9[k];

        t_10[k] = f_6 * sg0_10[k]
                  - f_7 * sg1_10[k]
                  + pb_x[k] * sh_10[k];

        t_11[k] = pb_z[k] * sh_6[k];

        t_12[k] = f_6 * sg0_12[k]
                  - f_7 * sg1_12[k]
                  + pb_x[k] * sh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, sg0_14, sg1_14, sh_9, \
                         sh_14, sh_15, sh_16, sh_17, sh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sh_9[k];

        t_14[k] = f_6 * sg0_14[k]
                  - f_7 * sg1_14[k]
                  + pb_x[k] * sh_14[k];

        t_15[k] = pb_x[k] * sh_15[k];

        t_16[k] = pb_x[k] * sh_16[k];

        t_17[k] = pb_x[k] * sh_17[k];

        t_18[k] = pb_x[k] * sh_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, sg0_10, sg0_12, \
                         sg1_10, sg1_12, sh_15, sh_17, sh_19, sh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_x[k] * sh_19[k];

        t_20[k] = pb_x[k] * sh_20[k];

        t_21[k] = f_0 * sg0_10[k]
                  - f_1 * sg1_10[k]
                  + pb_y[k] * sh_15[k];

        t_22[k] = pb_z[k] * sh_15[k];

        t_23[k] = f_2 * sg0_12[k]
                  - f_3 * sg1_12[k]
                  + pb_y[k] * sh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pb_z, sg0_13, sg0_14, sg1_13, sg1_14, \
                         sh_18, sh_19, sh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_4 * sg0_13[k]
                  - f_5 * sg1_13[k]
                  + pb_y[k] * sh_18[k];

        t_25[k] = f_6 * sg0_14[k]
                  - f_7 * sg1_14[k]
                  + pb_y[k] * sh_19[k];

        t_26[k] = pb_y[k] * sh_20[k];

        t_27[k] = f_0 * sg0_14[k]
                  - f_1 * sg1_14[k]
                  + pb_z[k] * sh_20[k];
    }
}

}  // namespace simdt2ceri
