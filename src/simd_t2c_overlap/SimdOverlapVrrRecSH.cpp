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


#include "SimdOverlapVrrRecSH.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_sh_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                          const size_t sf, const size_t sg, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, sf_0, sf_3, sf_5, \
                         sg_0, sg_2, sg_3, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pb_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * sg_0[k];

        t_2[k] = pb_z[k] * sg_0[k];

        t_3[k] = f_1 * sf_3[k]
                 + pb_x[k] * sg_3[k];

        t_4[k] = pb_y[k] * sg_2[k];

        t_5[k] = f_1 * sf_5[k]
                 + pb_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, sf_6, sf_9, sg_3, sg_5, \
                         sg_6, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * sf_6[k]
                 + pb_x[k] * sg_6[k];

        t_7[k] = pb_z[k] * sg_3[k];

        t_8[k] = pb_y[k] * sg_5[k];

        t_9[k] = f_2 * sf_9[k]
                 + pb_x[k] * sg_9[k];

        t_10[k] = pb_x[k] * sg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, sf_6, sg_10, \
                         sg_11, sg_12, sg_13, sg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_x[k] * sg_11[k];

        t_12[k] = pb_x[k] * sg_12[k];

        t_13[k] = pb_x[k] * sg_13[k];

        t_14[k] = pb_x[k] * sg_14[k];

        t_15[k] = f_0 * sf_6[k]
                  + pb_y[k] * sg_10[k];

        t_16[k] = pb_z[k] * sg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, sf_8, sf_9, sg_12, sg_13, \
                         sg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * sf_8[k]
                  + pb_y[k] * sg_12[k];

        t_18[k] = f_2 * sf_9[k]
                  + pb_y[k] * sg_13[k];

        t_19[k] = pb_y[k] * sg_14[k];

        t_20[k] = f_0 * sf_9[k]
                  + pb_z[k] * sg_14[k];
    }
}

}  // namespace simdovl
