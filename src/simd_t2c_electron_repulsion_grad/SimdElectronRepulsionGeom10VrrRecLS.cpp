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


#include "SimdElectronRepulsionGeom10VrrRecLS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_ls_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ks, const size_t ms,
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
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

    const auto *ms_0 = buffer.data(ms + 0);
    const auto *ms_1 = buffer.data(ms + 1);
    const auto *ms_2 = buffer.data(ms + 2);
    const auto *ms_3 = buffer.data(ms + 3);
    const auto *ms_4 = buffer.data(ms + 4);
    const auto *ms_5 = buffer.data(ms + 5);
    const auto *ms_6 = buffer.data(ms + 6);
    const auto *ms_7 = buffer.data(ms + 7);
    const auto *ms_8 = buffer.data(ms + 8);
    const auto *ms_9 = buffer.data(ms + 9);
    const auto *ms_10 = buffer.data(ms + 10);
    const auto *ms_11 = buffer.data(ms + 11);
    const auto *ms_12 = buffer.data(ms + 12);
    const auto *ms_13 = buffer.data(ms + 13);
    const auto *ms_14 = buffer.data(ms + 14);
    const auto *ms_15 = buffer.data(ms + 15);
    const auto *ms_16 = buffer.data(ms + 16);
    const auto *ms_17 = buffer.data(ms + 17);
    const auto *ms_18 = buffer.data(ms + 18);
    const auto *ms_19 = buffer.data(ms + 19);
    const auto *ms_20 = buffer.data(ms + 20);
    const auto *ms_21 = buffer.data(ms + 21);
    const auto *ms_22 = buffer.data(ms + 22);
    const auto *ms_23 = buffer.data(ms + 23);
    const auto *ms_24 = buffer.data(ms + 24);
    const auto *ms_25 = buffer.data(ms + 25);
    const auto *ms_26 = buffer.data(ms + 26);
    const auto *ms_27 = buffer.data(ms + 27);
    const auto *ms_28 = buffer.data(ms + 28);
    const auto *ms_29 = buffer.data(ms + 29);
    const auto *ms_30 = buffer.data(ms + 30);
    const auto *ms_31 = buffer.data(ms + 31);
    const auto *ms_32 = buffer.data(ms + 32);
    const auto *ms_33 = buffer.data(ms + 33);
    const auto *ms_34 = buffer.data(ms + 34);
    const auto *ms_35 = buffer.data(ms + 35);
    const auto *ms_36 = buffer.data(ms + 36);
    const auto *ms_37 = buffer.data(ms + 37);
    const auto *ms_38 = buffer.data(ms + 38);
    const auto *ms_39 = buffer.data(ms + 39);
    const auto *ms_40 = buffer.data(ms + 40);
    const auto *ms_41 = buffer.data(ms + 41);
    const auto *ms_42 = buffer.data(ms + 42);
    const auto *ms_43 = buffer.data(ms + 43);
    const auto *ms_44 = buffer.data(ms + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ks_0, ks_1, ks_2, ks_3, ks_4, ms_0, ms_1, \
                         ms_2, ms_3, ms_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * ks_0[k]
                 + f_0 * ms_0[k];

        t_1[k] = -7.0 * ks_1[k]
                 + f_0 * ms_1[k];

        t_2[k] = -7.0 * ks_2[k]
                 + f_0 * ms_2[k];

        t_3[k] = -6.0 * ks_3[k]
                 + f_0 * ms_3[k];

        t_4[k] = -6.0 * ks_4[k]
                 + f_0 * ms_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ks_5, ks_6, ks_7, ks_8, ks_9, ms_5, ms_6, \
                         ms_7, ms_8, ms_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * ks_5[k]
                 + f_0 * ms_5[k];

        t_6[k] = -5.0 * ks_6[k]
                 + f_0 * ms_6[k];

        t_7[k] = -5.0 * ks_7[k]
                 + f_0 * ms_7[k];

        t_8[k] = -5.0 * ks_8[k]
                 + f_0 * ms_8[k];

        t_9[k] = -5.0 * ks_9[k]
                 + f_0 * ms_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ks_10, ks_11, ks_12, ks_13, ks_14, \
                         ms_10, ms_11, ms_12, ms_13, ms_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * ks_10[k]
                  + f_0 * ms_10[k];

        t_11[k] = -4.0 * ks_11[k]
                  + f_0 * ms_11[k];

        t_12[k] = -4.0 * ks_12[k]
                  + f_0 * ms_12[k];

        t_13[k] = -4.0 * ks_13[k]
                  + f_0 * ms_13[k];

        t_14[k] = -4.0 * ks_14[k]
                  + f_0 * ms_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ks_15, ks_16, ks_17, ks_18, ks_19, \
                         ms_15, ms_16, ms_17, ms_18, ms_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * ks_15[k]
                  + f_0 * ms_15[k];

        t_16[k] = -3.0 * ks_16[k]
                  + f_0 * ms_16[k];

        t_17[k] = -3.0 * ks_17[k]
                  + f_0 * ms_17[k];

        t_18[k] = -3.0 * ks_18[k]
                  + f_0 * ms_18[k];

        t_19[k] = -3.0 * ks_19[k]
                  + f_0 * ms_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ks_20, ks_21, ks_22, ks_23, ks_24, \
                         ms_20, ms_21, ms_22, ms_23, ms_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * ks_20[k]
                  + f_0 * ms_20[k];

        t_21[k] = -2.0 * ks_21[k]
                  + f_0 * ms_21[k];

        t_22[k] = -2.0 * ks_22[k]
                  + f_0 * ms_22[k];

        t_23[k] = -2.0 * ks_23[k]
                  + f_0 * ms_23[k];

        t_24[k] = -2.0 * ks_24[k]
                  + f_0 * ms_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ks_25, ks_26, ks_27, ks_28, ks_29, \
                         ms_25, ms_26, ms_27, ms_28, ms_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * ks_25[k]
                  + f_0 * ms_25[k];

        t_26[k] = -2.0 * ks_26[k]
                  + f_0 * ms_26[k];

        t_27[k] = -2.0 * ks_27[k]
                  + f_0 * ms_27[k];

        t_28[k] = -ks_28[k]
                  + f_0 * ms_28[k];

        t_29[k] = -ks_29[k]
                  + f_0 * ms_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ks_30, ks_31, ks_32, ks_33, ks_34, \
                         ms_30, ms_31, ms_32, ms_33, ms_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -ks_30[k]
                  + f_0 * ms_30[k];

        t_31[k] = -ks_31[k]
                  + f_0 * ms_31[k];

        t_32[k] = -ks_32[k]
                  + f_0 * ms_32[k];

        t_33[k] = -ks_33[k]
                  + f_0 * ms_33[k];

        t_34[k] = -ks_34[k]
                  + f_0 * ms_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, t_41, ks_35, ms_35, ms_36, ms_37, \
                         ms_38, ms_39, ms_40, ms_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -ks_35[k]
                  + f_0 * ms_35[k];

        t_36[k] = f_0 * ms_36[k];

        t_37[k] = f_0 * ms_37[k];

        t_38[k] = f_0 * ms_38[k];

        t_39[k] = f_0 * ms_39[k];

        t_40[k] = f_0 * ms_40[k];

        t_41[k] = f_0 * ms_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, ms_42, ms_43, ms_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_0 * ms_42[k];

        t_43[k] = f_0 * ms_43[k];

        t_44[k] = f_0 * ms_44[k];
    }
}

auto
compute_prim_geom_10_ls_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ks, const size_t ms,
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
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

    const auto *ms_1 = buffer.data(ms + 1);
    const auto *ms_3 = buffer.data(ms + 3);
    const auto *ms_4 = buffer.data(ms + 4);
    const auto *ms_6 = buffer.data(ms + 6);
    const auto *ms_7 = buffer.data(ms + 7);
    const auto *ms_8 = buffer.data(ms + 8);
    const auto *ms_10 = buffer.data(ms + 10);
    const auto *ms_11 = buffer.data(ms + 11);
    const auto *ms_12 = buffer.data(ms + 12);
    const auto *ms_13 = buffer.data(ms + 13);
    const auto *ms_15 = buffer.data(ms + 15);
    const auto *ms_16 = buffer.data(ms + 16);
    const auto *ms_17 = buffer.data(ms + 17);
    const auto *ms_18 = buffer.data(ms + 18);
    const auto *ms_19 = buffer.data(ms + 19);
    const auto *ms_21 = buffer.data(ms + 21);
    const auto *ms_22 = buffer.data(ms + 22);
    const auto *ms_23 = buffer.data(ms + 23);
    const auto *ms_24 = buffer.data(ms + 24);
    const auto *ms_25 = buffer.data(ms + 25);
    const auto *ms_26 = buffer.data(ms + 26);
    const auto *ms_28 = buffer.data(ms + 28);
    const auto *ms_29 = buffer.data(ms + 29);
    const auto *ms_30 = buffer.data(ms + 30);
    const auto *ms_31 = buffer.data(ms + 31);
    const auto *ms_32 = buffer.data(ms + 32);
    const auto *ms_33 = buffer.data(ms + 33);
    const auto *ms_34 = buffer.data(ms + 34);
    const auto *ms_36 = buffer.data(ms + 36);
    const auto *ms_37 = buffer.data(ms + 37);
    const auto *ms_38 = buffer.data(ms + 38);
    const auto *ms_39 = buffer.data(ms + 39);
    const auto *ms_40 = buffer.data(ms + 40);
    const auto *ms_41 = buffer.data(ms + 41);
    const auto *ms_42 = buffer.data(ms + 42);
    const auto *ms_43 = buffer.data(ms + 43);
    const auto *ms_45 = buffer.data(ms + 45);
    const auto *ms_46 = buffer.data(ms + 46);
    const auto *ms_47 = buffer.data(ms + 47);
    const auto *ms_48 = buffer.data(ms + 48);
    const auto *ms_49 = buffer.data(ms + 49);
    const auto *ms_50 = buffer.data(ms + 50);
    const auto *ms_51 = buffer.data(ms + 51);
    const auto *ms_52 = buffer.data(ms + 52);
    const auto *ms_53 = buffer.data(ms + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ks_0, ks_1, ks_2, ms_1, ms_3, ms_4, \
                         ms_6, ms_7, ms_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ms_1[k];

        t_1[k] = -ks_0[k]
                 + f_0 * ms_3[k];

        t_2[k] = f_0 * ms_4[k];

        t_3[k] = -2.0 * ks_1[k]
                 + f_0 * ms_6[k];

        t_4[k] = -ks_2[k]
                 + f_0 * ms_7[k];

        t_5[k] = f_0 * ms_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ks_3, ks_4, ks_5, ks_6, ms_10, ms_11, \
                         ms_12, ms_13, ms_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -3.0 * ks_3[k]
                 + f_0 * ms_10[k];

        t_7[k] = -2.0 * ks_4[k]
                 + f_0 * ms_11[k];

        t_8[k] = -ks_5[k]
                 + f_0 * ms_12[k];

        t_9[k] = f_0 * ms_13[k];

        t_10[k] = -4.0 * ks_6[k]
                  + f_0 * ms_15[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, ks_7, ks_8, ks_9, ks_10, ms_16, ms_17, \
                         ms_18, ms_19, ms_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -3.0 * ks_7[k]
                  + f_0 * ms_16[k];

        t_12[k] = -2.0 * ks_8[k]
                  + f_0 * ms_17[k];

        t_13[k] = -ks_9[k]
                  + f_0 * ms_18[k];

        t_14[k] = f_0 * ms_19[k];

        t_15[k] = -5.0 * ks_10[k]
                  + f_0 * ms_21[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, ks_11, ks_12, ks_13, ks_14, ms_22, \
                         ms_23, ms_24, ms_25, ms_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -4.0 * ks_11[k]
                  + f_0 * ms_22[k];

        t_17[k] = -3.0 * ks_12[k]
                  + f_0 * ms_23[k];

        t_18[k] = -2.0 * ks_13[k]
                  + f_0 * ms_24[k];

        t_19[k] = -ks_14[k]
                  + f_0 * ms_25[k];

        t_20[k] = f_0 * ms_26[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, ks_15, ks_16, ks_17, ks_18, ks_19, \
                         ms_28, ms_29, ms_30, ms_31, ms_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -6.0 * ks_15[k]
                  + f_0 * ms_28[k];

        t_22[k] = -5.0 * ks_16[k]
                  + f_0 * ms_29[k];

        t_23[k] = -4.0 * ks_17[k]
                  + f_0 * ms_30[k];

        t_24[k] = -3.0 * ks_18[k]
                  + f_0 * ms_31[k];

        t_25[k] = -2.0 * ks_19[k]
                  + f_0 * ms_32[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, ks_20, ks_21, ks_22, ks_23, ms_33, \
                         ms_34, ms_36, ms_37, ms_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -ks_20[k]
                  + f_0 * ms_33[k];

        t_27[k] = f_0 * ms_34[k];

        t_28[k] = -7.0 * ks_21[k]
                  + f_0 * ms_36[k];

        t_29[k] = -6.0 * ks_22[k]
                  + f_0 * ms_37[k];

        t_30[k] = -5.0 * ks_23[k]
                  + f_0 * ms_38[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ks_24, ks_25, ks_26, ks_27, ms_39, \
                         ms_40, ms_41, ms_42, ms_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -4.0 * ks_24[k]
                  + f_0 * ms_39[k];

        t_32[k] = -3.0 * ks_25[k]
                  + f_0 * ms_40[k];

        t_33[k] = -2.0 * ks_26[k]
                  + f_0 * ms_41[k];

        t_34[k] = -ks_27[k]
                  + f_0 * ms_42[k];

        t_35[k] = f_0 * ms_43[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ks_28, ks_29, ks_30, ks_31, ks_32, \
                         ms_45, ms_46, ms_47, ms_48, ms_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -8.0 * ks_28[k]
                  + f_0 * ms_45[k];

        t_37[k] = -7.0 * ks_29[k]
                  + f_0 * ms_46[k];

        t_38[k] = -6.0 * ks_30[k]
                  + f_0 * ms_47[k];

        t_39[k] = -5.0 * ks_31[k]
                  + f_0 * ms_48[k];

        t_40[k] = -4.0 * ks_32[k]
                  + f_0 * ms_49[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ks_33, ks_34, ks_35, ms_50, ms_51, ms_52, \
                         ms_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * ks_33[k]
                  + f_0 * ms_50[k];

        t_42[k] = -2.0 * ks_34[k]
                  + f_0 * ms_51[k];

        t_43[k] = -ks_35[k]
                  + f_0 * ms_52[k];

        t_44[k] = f_0 * ms_53[k];
    }
}

auto
compute_prim_geom_10_ls_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ks, const size_t ms,
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
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

    const auto *ms_2 = buffer.data(ms + 2);
    const auto *ms_4 = buffer.data(ms + 4);
    const auto *ms_5 = buffer.data(ms + 5);
    const auto *ms_7 = buffer.data(ms + 7);
    const auto *ms_8 = buffer.data(ms + 8);
    const auto *ms_9 = buffer.data(ms + 9);
    const auto *ms_11 = buffer.data(ms + 11);
    const auto *ms_12 = buffer.data(ms + 12);
    const auto *ms_13 = buffer.data(ms + 13);
    const auto *ms_14 = buffer.data(ms + 14);
    const auto *ms_16 = buffer.data(ms + 16);
    const auto *ms_17 = buffer.data(ms + 17);
    const auto *ms_18 = buffer.data(ms + 18);
    const auto *ms_19 = buffer.data(ms + 19);
    const auto *ms_20 = buffer.data(ms + 20);
    const auto *ms_22 = buffer.data(ms + 22);
    const auto *ms_23 = buffer.data(ms + 23);
    const auto *ms_24 = buffer.data(ms + 24);
    const auto *ms_25 = buffer.data(ms + 25);
    const auto *ms_26 = buffer.data(ms + 26);
    const auto *ms_27 = buffer.data(ms + 27);
    const auto *ms_29 = buffer.data(ms + 29);
    const auto *ms_30 = buffer.data(ms + 30);
    const auto *ms_31 = buffer.data(ms + 31);
    const auto *ms_32 = buffer.data(ms + 32);
    const auto *ms_33 = buffer.data(ms + 33);
    const auto *ms_34 = buffer.data(ms + 34);
    const auto *ms_35 = buffer.data(ms + 35);
    const auto *ms_37 = buffer.data(ms + 37);
    const auto *ms_38 = buffer.data(ms + 38);
    const auto *ms_39 = buffer.data(ms + 39);
    const auto *ms_40 = buffer.data(ms + 40);
    const auto *ms_41 = buffer.data(ms + 41);
    const auto *ms_42 = buffer.data(ms + 42);
    const auto *ms_43 = buffer.data(ms + 43);
    const auto *ms_44 = buffer.data(ms + 44);
    const auto *ms_46 = buffer.data(ms + 46);
    const auto *ms_47 = buffer.data(ms + 47);
    const auto *ms_48 = buffer.data(ms + 48);
    const auto *ms_49 = buffer.data(ms + 49);
    const auto *ms_50 = buffer.data(ms + 50);
    const auto *ms_51 = buffer.data(ms + 51);
    const auto *ms_52 = buffer.data(ms + 52);
    const auto *ms_53 = buffer.data(ms + 53);
    const auto *ms_54 = buffer.data(ms + 54);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ks_0, ks_1, ks_2, ms_2, ms_4, ms_5, \
                         ms_7, ms_8, ms_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ms_2[k];

        t_1[k] = f_0 * ms_4[k];

        t_2[k] = -ks_0[k]
                 + f_0 * ms_5[k];

        t_3[k] = f_0 * ms_7[k];

        t_4[k] = -ks_1[k]
                 + f_0 * ms_8[k];

        t_5[k] = -2.0 * ks_2[k]
                 + f_0 * ms_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, ks_3, ks_4, ks_5, ks_6, ms_11, ms_12, \
                         ms_13, ms_14, ms_16, ms_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ms_11[k];

        t_7[k] = -ks_3[k]
                 + f_0 * ms_12[k];

        t_8[k] = -2.0 * ks_4[k]
                 + f_0 * ms_13[k];

        t_9[k] = -3.0 * ks_5[k]
                 + f_0 * ms_14[k];

        t_10[k] = f_0 * ms_16[k];

        t_11[k] = -ks_6[k]
                  + f_0 * ms_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, ks_7, ks_8, ks_9, ks_10, ms_18, ms_19, \
                         ms_20, ms_22, ms_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -2.0 * ks_7[k]
                  + f_0 * ms_18[k];

        t_13[k] = -3.0 * ks_8[k]
                  + f_0 * ms_19[k];

        t_14[k] = -4.0 * ks_9[k]
                  + f_0 * ms_20[k];

        t_15[k] = f_0 * ms_22[k];

        t_16[k] = -ks_10[k]
                  + f_0 * ms_23[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ks_11, ks_12, ks_13, ks_14, ms_24, \
                         ms_25, ms_26, ms_27, ms_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -2.0 * ks_11[k]
                  + f_0 * ms_24[k];

        t_18[k] = -3.0 * ks_12[k]
                  + f_0 * ms_25[k];

        t_19[k] = -4.0 * ks_13[k]
                  + f_0 * ms_26[k];

        t_20[k] = -5.0 * ks_14[k]
                  + f_0 * ms_27[k];

        t_21[k] = f_0 * ms_29[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ks_15, ks_16, ks_17, ks_18, ks_19, \
                         ms_30, ms_31, ms_32, ms_33, ms_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -ks_15[k]
                  + f_0 * ms_30[k];

        t_23[k] = -2.0 * ks_16[k]
                  + f_0 * ms_31[k];

        t_24[k] = -3.0 * ks_17[k]
                  + f_0 * ms_32[k];

        t_25[k] = -4.0 * ks_18[k]
                  + f_0 * ms_33[k];

        t_26[k] = -5.0 * ks_19[k]
                  + f_0 * ms_34[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ks_20, ks_21, ks_22, ks_23, ms_35, \
                         ms_37, ms_38, ms_39, ms_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -6.0 * ks_20[k]
                  + f_0 * ms_35[k];

        t_28[k] = f_0 * ms_37[k];

        t_29[k] = -ks_21[k]
                  + f_0 * ms_38[k];

        t_30[k] = -2.0 * ks_22[k]
                  + f_0 * ms_39[k];

        t_31[k] = -3.0 * ks_23[k]
                  + f_0 * ms_40[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, ks_24, ks_25, ks_26, ks_27, ms_41, \
                         ms_42, ms_43, ms_44, ms_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -4.0 * ks_24[k]
                  + f_0 * ms_41[k];

        t_33[k] = -5.0 * ks_25[k]
                  + f_0 * ms_42[k];

        t_34[k] = -6.0 * ks_26[k]
                  + f_0 * ms_43[k];

        t_35[k] = -7.0 * ks_27[k]
                  + f_0 * ms_44[k];

        t_36[k] = f_0 * ms_46[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, ks_28, ks_29, ks_30, ks_31, ks_32, \
                         ms_47, ms_48, ms_49, ms_50, ms_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -ks_28[k]
                  + f_0 * ms_47[k];

        t_38[k] = -2.0 * ks_29[k]
                  + f_0 * ms_48[k];

        t_39[k] = -3.0 * ks_30[k]
                  + f_0 * ms_49[k];

        t_40[k] = -4.0 * ks_31[k]
                  + f_0 * ms_50[k];

        t_41[k] = -5.0 * ks_32[k]
                  + f_0 * ms_51[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, ks_33, ks_34, ks_35, ms_52, ms_53, \
                         ms_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -6.0 * ks_33[k]
                  + f_0 * ms_52[k];

        t_43[k] = -7.0 * ks_34[k]
                  + f_0 * ms_53[k];

        t_44[k] = -8.0 * ks_35[k]
                  + f_0 * ms_54[k];
    }
}

}  // namespace simdt2ceri
