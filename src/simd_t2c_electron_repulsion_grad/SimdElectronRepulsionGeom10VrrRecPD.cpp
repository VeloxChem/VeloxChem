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


#include "SimdElectronRepulsionGeom10VrrRecPD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sd, const size_t dd,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_4 = buffer.data(sd + 4);
    const auto *sd_5 = buffer.data(sd + 5);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, sd_0, sd_1, sd_2, sd_3, sd_4, dd_0, dd_1, \
                         dd_2, dd_3, dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sd_0[k]
                 + f_0 * dd_0[k];

        t_1[k] = -sd_1[k]
                 + f_0 * dd_1[k];

        t_2[k] = -sd_2[k]
                 + f_0 * dd_2[k];

        t_3[k] = -sd_3[k]
                 + f_0 * dd_3[k];

        t_4[k] = -sd_4[k]
                 + f_0 * dd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, t_11, sd_5, dd_5, dd_6, dd_7, dd_8, \
                         dd_9, dd_10, dd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -sd_5[k]
                 + f_0 * dd_5[k];

        t_6[k] = f_0 * dd_6[k];

        t_7[k] = f_0 * dd_7[k];

        t_8[k] = f_0 * dd_8[k];

        t_9[k] = f_0 * dd_9[k];

        t_10[k] = f_0 * dd_10[k];

        t_11[k] = f_0 * dd_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, dd_12, dd_13, dd_14, dd_15, \
                         dd_16, dd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * dd_12[k];

        t_13[k] = f_0 * dd_13[k];

        t_14[k] = f_0 * dd_14[k];

        t_15[k] = f_0 * dd_15[k];

        t_16[k] = f_0 * dd_16[k];

        t_17[k] = f_0 * dd_17[k];
    }
}

auto
compute_prim_geom_10_pd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sd, const size_t dd,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_4 = buffer.data(sd + 4);
    const auto *sd_5 = buffer.data(sd + 5);

    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_20 = buffer.data(dd + 20);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_22 = buffer.data(dd + 22);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_25 = buffer.data(dd + 25);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, sd_0, dd_6, dd_7, dd_8, dd_9, \
                         dd_10, dd_11, dd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_6[k];

        t_1[k] = f_0 * dd_7[k];

        t_2[k] = f_0 * dd_8[k];

        t_3[k] = f_0 * dd_9[k];

        t_4[k] = f_0 * dd_10[k];

        t_5[k] = f_0 * dd_11[k];

        t_6[k] = -sd_0[k]
                 + f_0 * dd_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, sd_1, sd_2, sd_3, sd_4, sd_5, dd_19, \
                         dd_20, dd_21, dd_22, dd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -sd_1[k]
                 + f_0 * dd_19[k];

        t_8[k] = -sd_2[k]
                 + f_0 * dd_20[k];

        t_9[k] = -sd_3[k]
                 + f_0 * dd_21[k];

        t_10[k] = -sd_4[k]
                  + f_0 * dd_22[k];

        t_11[k] = -sd_5[k]
                  + f_0 * dd_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, dd_24, dd_25, dd_26, dd_27, \
                         dd_28, dd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * dd_24[k];

        t_13[k] = f_0 * dd_25[k];

        t_14[k] = f_0 * dd_26[k];

        t_15[k] = f_0 * dd_27[k];

        t_16[k] = f_0 * dd_28[k];

        t_17[k] = f_0 * dd_29[k];
    }
}

auto
compute_prim_geom_10_pd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sd, const size_t dd,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_4 = buffer.data(sd + 4);
    const auto *sd_5 = buffer.data(sd + 5);

    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_25 = buffer.data(dd + 25);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_31 = buffer.data(dd + 31);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dd_12, dd_13, dd_14, dd_15, \
                         dd_16, dd_17, dd_24, dd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_12[k];

        t_1[k] = f_0 * dd_13[k];

        t_2[k] = f_0 * dd_14[k];

        t_3[k] = f_0 * dd_15[k];

        t_4[k] = f_0 * dd_16[k];

        t_5[k] = f_0 * dd_17[k];

        t_6[k] = f_0 * dd_24[k];

        t_7[k] = f_0 * dd_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, sd_0, sd_1, dd_26, dd_27, dd_28, \
                         dd_29, dd_30, dd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dd_26[k];

        t_9[k] = f_0 * dd_27[k];

        t_10[k] = f_0 * dd_28[k];

        t_11[k] = f_0 * dd_29[k];

        t_12[k] = -sd_0[k]
                  + f_0 * dd_30[k];

        t_13[k] = -sd_1[k]
                  + f_0 * dd_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, sd_2, sd_3, sd_4, sd_5, dd_32, dd_33, dd_34, \
                         dd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -sd_2[k]
                  + f_0 * dd_32[k];

        t_15[k] = -sd_3[k]
                  + f_0 * dd_33[k];

        t_16[k] = -sd_4[k]
                  + f_0 * dd_34[k];

        t_17[k] = -sd_5[k]
                  + f_0 * dd_35[k];
    }
}

}  // namespace simdt2ceri
