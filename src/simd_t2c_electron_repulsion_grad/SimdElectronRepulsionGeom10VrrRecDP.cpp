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


#include "SimdElectronRepulsionGeom10VrrRecDP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_dp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pp, const size_t fp,
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

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_1 = buffer.data(pp + 1);
    const auto *pp_2 = buffer.data(pp + 2);
    const auto *pp_3 = buffer.data(pp + 3);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_5 = buffer.data(pp + 5);
    const auto *pp_6 = buffer.data(pp + 6);
    const auto *pp_7 = buffer.data(pp + 7);
    const auto *pp_8 = buffer.data(pp + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pp_0, pp_1, pp_2, pp_3, pp_4, fp_0, fp_1, \
                         fp_2, fp_3, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pp_0[k]
                 + f_0 * fp_0[k];

        t_1[k] = -2.0 * pp_1[k]
                 + f_0 * fp_1[k];

        t_2[k] = -2.0 * pp_2[k]
                 + f_0 * fp_2[k];

        t_3[k] = -pp_3[k]
                 + f_0 * fp_3[k];

        t_4[k] = -pp_4[k]
                 + f_0 * fp_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pp_5, pp_6, pp_7, pp_8, fp_5, fp_6, \
                         fp_7, fp_8, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -pp_5[k]
                 + f_0 * fp_5[k];

        t_6[k] = -pp_6[k]
                 + f_0 * fp_6[k];

        t_7[k] = -pp_7[k]
                 + f_0 * fp_7[k];

        t_8[k] = -pp_8[k]
                 + f_0 * fp_8[k];

        t_9[k] = f_0 * fp_9[k];

        t_10[k] = f_0 * fp_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, t_17, fp_11, fp_12, fp_13, fp_14, \
                         fp_15, fp_16, fp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * fp_11[k];

        t_12[k] = f_0 * fp_12[k];

        t_13[k] = f_0 * fp_13[k];

        t_14[k] = f_0 * fp_14[k];

        t_15[k] = f_0 * fp_15[k];

        t_16[k] = f_0 * fp_16[k];

        t_17[k] = f_0 * fp_17[k];
    }
}

auto
compute_prim_geom_10_dp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pp, const size_t fp,
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

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_1 = buffer.data(pp + 1);
    const auto *pp_2 = buffer.data(pp + 2);
    const auto *pp_3 = buffer.data(pp + 3);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_5 = buffer.data(pp + 5);
    const auto *pp_6 = buffer.data(pp + 6);
    const auto *pp_7 = buffer.data(pp + 7);
    const auto *pp_8 = buffer.data(pp + 8);

    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pp_0, pp_1, pp_2, fp_3, fp_4, fp_5, \
                         fp_9, fp_10, fp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_3[k];

        t_1[k] = f_0 * fp_4[k];

        t_2[k] = f_0 * fp_5[k];

        t_3[k] = -pp_0[k]
                 + f_0 * fp_9[k];

        t_4[k] = -pp_1[k]
                 + f_0 * fp_10[k];

        t_5[k] = -pp_2[k]
                 + f_0 * fp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pp_3, pp_4, pp_5, fp_12, fp_13, \
                         fp_14, fp_18, fp_19, fp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * fp_12[k];

        t_7[k] = f_0 * fp_13[k];

        t_8[k] = f_0 * fp_14[k];

        t_9[k] = -2.0 * pp_3[k]
                 + f_0 * fp_18[k];

        t_10[k] = -2.0 * pp_4[k]
                  + f_0 * fp_19[k];

        t_11[k] = -2.0 * pp_5[k]
                  + f_0 * fp_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pp_6, pp_7, pp_8, fp_21, fp_22, \
                         fp_23, fp_24, fp_25, fp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -pp_6[k]
                  + f_0 * fp_21[k];

        t_13[k] = -pp_7[k]
                  + f_0 * fp_22[k];

        t_14[k] = -pp_8[k]
                  + f_0 * fp_23[k];

        t_15[k] = f_0 * fp_24[k];

        t_16[k] = f_0 * fp_25[k];

        t_17[k] = f_0 * fp_26[k];
    }
}

auto
compute_prim_geom_10_dp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pp, const size_t fp,
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

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_1 = buffer.data(pp + 1);
    const auto *pp_2 = buffer.data(pp + 2);
    const auto *pp_3 = buffer.data(pp + 3);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_5 = buffer.data(pp + 5);
    const auto *pp_6 = buffer.data(pp + 6);
    const auto *pp_7 = buffer.data(pp + 7);
    const auto *pp_8 = buffer.data(pp + 8);

    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pp_0, fp_6, fp_7, fp_8, fp_12, \
                         fp_13, fp_14, fp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_6[k];

        t_1[k] = f_0 * fp_7[k];

        t_2[k] = f_0 * fp_8[k];

        t_3[k] = f_0 * fp_12[k];

        t_4[k] = f_0 * fp_13[k];

        t_5[k] = f_0 * fp_14[k];

        t_6[k] = -pp_0[k]
                 + f_0 * fp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pp_1, pp_2, pp_3, fp_16, fp_17, \
                         fp_21, fp_22, fp_23, fp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -pp_1[k]
                 + f_0 * fp_16[k];

        t_8[k] = -pp_2[k]
                 + f_0 * fp_17[k];

        t_9[k] = f_0 * fp_21[k];

        t_10[k] = f_0 * fp_22[k];

        t_11[k] = f_0 * fp_23[k];

        t_12[k] = -pp_3[k]
                  + f_0 * fp_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pp_4, pp_5, pp_6, pp_7, pp_8, fp_25, \
                         fp_26, fp_27, fp_28, fp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -pp_4[k]
                  + f_0 * fp_25[k];

        t_14[k] = -pp_5[k]
                  + f_0 * fp_26[k];

        t_15[k] = -2.0 * pp_6[k]
                  + f_0 * fp_27[k];

        t_16[k] = -2.0 * pp_7[k]
                  + f_0 * fp_28[k];

        t_17[k] = -2.0 * pp_8[k]
                  + f_0 * fp_29[k];
    }
}

}  // namespace simdt2ceri
