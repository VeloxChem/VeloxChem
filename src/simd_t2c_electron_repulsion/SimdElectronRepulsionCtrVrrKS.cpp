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


#include "SimdElectronRepulsionCtrVrrKS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_ks_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pa, const size_t hs0, const size_t hs1,
                                    const size_t is, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.640625 * std::sqrt(429.0) / alpha;
    const auto f_1 = 0.234375 * std::sqrt(429.0) / alpha;
    const auto f_2 = 1.640625 * std::sqrt(429.0) * beta / (alpha * p);
    const auto f_3 = 0.234375 * std::sqrt(429.0) * beta / (alpha * p);
    const auto f_4 = 0.21875 * std::sqrt(429.0);
    const auto f_5 = 1.09375 * std::sqrt(429.0);
    const auto f_6 = 0.65625 * std::sqrt(429.0);
    const auto f_7 = 0.03125 * std::sqrt(429.0);
    const auto f_8 = 0.1875 * std::sqrt(6006.0);
    const auto f_9 = 0.625 * std::sqrt(6006.0);
    const auto f_10 = 0.234375 * std::sqrt(231.0) / alpha;
    const auto f_11 = 0.046875 * std::sqrt(231.0) / alpha;
    const auto f_12 = 1.125 * std::sqrt(231.0) / alpha;
    const auto f_13 = 0.234375 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_14 = 0.046875 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_15 = 1.125 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_16 = 0.15625 * std::sqrt(231.0);
    const auto f_17 = 1.875 * std::sqrt(231.0);
    const auto f_18 = 0.28125 * std::sqrt(231.0);
    const auto f_19 = 3.75 * std::sqrt(231.0);
    const auto f_20 = 0.03125 * std::sqrt(231.0);
    const auto f_21 = 0.375 * std::sqrt(231.0);
    const auto f_22 = 0.75 * std::sqrt(231.0);
    const auto f_23 = 2.5 * std::sqrt(231.0);
    const auto f_24 = 0.703125 * std::sqrt(21.0) / alpha;
    const auto f_25 = 0.234375 * std::sqrt(21.0) / alpha;
    const auto f_26 = 1.875 * std::sqrt(21.0) / alpha;
    const auto f_27 = 2.5 * std::sqrt(21.0) / alpha;
    const auto f_28 = 0.703125 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_29 = 0.234375 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_30 = 1.875 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_31 = 2.5 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_32 = 0.28125 * std::sqrt(21.0);
    const auto f_33 = 5.625 * std::sqrt(21.0);
    const auto f_34 = 0.46875 * std::sqrt(21.0);
    const auto f_35 = 7.5 * std::sqrt(21.0);
    const auto f_36 = 0.09375 * std::sqrt(21.0);
    const auto f_37 = 3.75 * std::sqrt(21.0);
    const auto f_38 = 1.875 * std::sqrt(21.0);
    const auto f_39 = 2.5 * std::sqrt(21.0);
    const auto f_40 = 0.9375 * std::sqrt(42.0);
    const auto f_41 = 1.875 * std::sqrt(42.0);
    const auto f_42 = 5.0 * std::sqrt(42.0);
    const auto f_43 = 3.0 * std::sqrt(42.0);
    const auto f_44 = 0.703125 * std::sqrt(7.0) / alpha;
    const auto f_45 = 11.25 * std::sqrt(7.0) / alpha;
    const auto f_46 = 7.5 * std::sqrt(7.0) / alpha;
    const auto f_47 = 0.703125 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_48 = 11.25 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_49 = 7.5 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_50 = 0.15625 * std::sqrt(7.0);
    const auto f_51 = 3.75 * std::sqrt(7.0);
    const auto f_52 = 0.46875 * std::sqrt(7.0);
    const auto f_53 = 7.5 * std::sqrt(7.0);
    const auto f_54 = 2.0 * std::sqrt(7.0);
    const auto f_55 = 19.6875 / alpha;
    const auto f_56 = 32.8125 / alpha;
    const auto f_57 = 7.5 / alpha;
    const auto f_58 = 19.6875 * beta / (alpha * p);
    const auto f_59 = 32.8125 * beta / (alpha * p);
    const auto f_60 = 7.5 * beta / (alpha * p);
    const auto f_61 = 0.46875 * std::sqrt(7.0) / alpha;
    const auto f_62 = 0.9375 * std::sqrt(7.0) / alpha;
    const auto f_63 = 0.46875 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_64 = 0.9375 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_65 = 3.75 * std::sqrt(42.0) / alpha;
    const auto f_66 = 3.75 * std::sqrt(42.0) * beta / (alpha * p);
    const auto f_67 = 0.46875 * std::sqrt(42.0);
    const auto f_68 = 2.5 * std::sqrt(42.0);
    const auto f_69 = 1.5 * std::sqrt(42.0);
    const auto f_70 = 0.28125 * std::sqrt(21.0) / alpha;
    const auto f_71 = 0.1875 * std::sqrt(21.0) / alpha;
    const auto f_72 = 3.75 * std::sqrt(21.0) / alpha;
    const auto f_73 = 0.46875 * std::sqrt(21.0) / alpha;
    const auto f_74 = 0.28125 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_75 = 0.1875 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_76 = 3.75 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_77 = 0.46875 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_78 = 0.9375 * std::sqrt(231.0) / alpha;
    const auto f_79 = 0.9375 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_80 = 0.1875 * std::sqrt(231.0);
    const auto f_81 = 0.9375 * std::sqrt(231.0);
    const auto f_82 = 0.625 * std::sqrt(231.0);
    const auto f_83 = 0.09375 * std::sqrt(231.0) / alpha;
    const auto f_84 = 0.5625 * std::sqrt(231.0) / alpha;
    const auto f_85 = 0.75 * std::sqrt(231.0) / alpha;
    const auto f_86 = 0.15625 * std::sqrt(231.0) / alpha;
    const auto f_87 = 3.75 * std::sqrt(231.0) / alpha;
    const auto f_88 = 0.09375 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_89 = 0.5625 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_90 = 0.75 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_91 = 0.15625 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_92 = 3.75 * std::sqrt(231.0) * beta / (alpha * p);
    const auto f_93 = 0.03125 * std::sqrt(6006.0);
    const auto f_94 = 0.46875 * std::sqrt(6006.0);
    const auto f_95 = 0.09375 * std::sqrt(429.0) / alpha;
    const auto f_96 = 1.3125 * std::sqrt(429.0) / alpha;
    const auto f_97 = 1.09375 * std::sqrt(429.0) / alpha;
    const auto f_98 = 0.09375 * std::sqrt(429.0) * beta / (alpha * p);
    const auto f_99 = 1.3125 * std::sqrt(429.0) * beta / (alpha * p);
    const auto f_100 = 1.09375 * std::sqrt(429.0) * beta / (alpha * p);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);

#pragma omp simd aligned(pa_x, pa_y, pa_z, hs0_3, hs0_8, hs1_3, hs1_8, is_0, is_1, is_4, is_9, \
                         is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * hs0_3[k]
                  + f_1 * hs0_8[k]
                  + f_2 * hs1_3[k]
                  - f_3 * hs1_8[k]
                  + f_4 * pa_y[k] * is_0[k]
                  - f_5 * pa_x[k] * is_4[k]
                  + f_6 * pa_x[k] * is_9[k]
                  - f_7 * pa_y[k] * is_13[k];

        g_1[k] += f_8 * pa_y[k] * is_1[k]
                  - f_9 * pa_z[k] * is_4[k]
                  + f_8 * pa_x[k] * is_14[k];
    }

#pragma omp simd aligned(pa_x, pa_y, hs0_3, hs0_8, hs0_9, hs1_3, hs1_8, hs1_9, is_0, is_3, \
                         is_4, is_9, is_10, is_13, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += f_10 * hs0_3[k]
                  + f_11 * hs0_8[k]
                  - f_12 * hs0_9[k]
                  - f_13 * hs1_3[k]
                  - f_14 * hs1_8[k]
                  + f_15 * hs1_9[k]
                  - f_16 * pa_y[k] * is_0[k]
                  + f_17 * pa_y[k] * is_3[k]
                  + f_16 * pa_x[k] * is_4[k]
                  + f_18 * pa_x[k] * is_9[k]
                  - f_19 * pa_x[k] * is_10[k]
                  - f_20 * pa_y[k] * is_13[k]
                  + f_21 * pa_y[k] * is_15[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is_1, is_5, is_14, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_22 * pa_y[k] * is_1[k]
                  + f_23 * pa_y[k] * is_5[k]
                  + f_22 * pa_x[k] * is_14[k]
                  - f_23 * pa_x[k] * is_16[k];
    }

#pragma omp simd aligned(pa_x, pa_y, hs0_3, hs0_8, hs0_9, hs0_11, hs1_3, hs1_8, hs1_9, hs1_11, \
                         is_0, is_3, is_4, is_8, is_9, is_10, is_13, is_15, \
                         is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_24 * hs0_3[k]
                  - f_25 * hs0_8[k]
                  + f_26 * hs0_9[k]
                  - f_27 * hs0_11[k]
                  - f_28 * hs1_3[k]
                  + f_29 * hs1_8[k]
                  - f_30 * hs1_9[k]
                  + f_31 * hs1_11[k]
                  + f_32 * pa_y[k] * is_0[k]
                  - f_33 * pa_y[k] * is_3[k]
                  + f_34 * pa_x[k] * is_4[k]
                  + f_35 * pa_y[k] * is_8[k]
                  + f_36 * pa_x[k] * is_9[k]
                  - f_37 * pa_x[k] * is_10[k]
                  - f_36 * pa_y[k] * is_13[k]
                  + f_38 * pa_y[k] * is_15[k]
                  - f_39 * pa_y[k] * is_17[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, is_1, is_4, is_5, is_14, is_16, \
                         is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_40 * pa_y[k] * is_1[k]
                  + f_41 * pa_z[k] * is_4[k]
                  - f_42 * pa_y[k] * is_5[k]
                  + f_40 * pa_x[k] * is_14[k]
                  - f_42 * pa_x[k] * is_16[k]
                  + f_43 * pa_x[k] * is_18[k];
    }

#pragma omp simd aligned(pa_x, pa_y, hs0_3, hs0_8, hs0_9, hs0_11, hs1_3, hs1_8, hs1_9, hs1_11, \
                         is_0, is_3, is_4, is_8, is_9, is_10, is_13, is_15, is_17, \
                         is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_44 * hs0_3[k]
                  - f_44 * hs0_8[k]
                  + f_45 * hs0_9[k]
                  - f_46 * hs0_11[k]
                  + f_47 * hs1_3[k]
                  + f_47 * hs1_8[k]
                  - f_48 * hs1_9[k]
                  + f_49 * hs1_11[k]
                  - f_50 * pa_y[k] * is_0[k]
                  + f_51 * pa_y[k] * is_3[k]
                  - f_52 * pa_x[k] * is_4[k]
                  - f_53 * pa_y[k] * is_8[k]
                  - f_52 * pa_x[k] * is_9[k]
                  + f_53 * pa_x[k] * is_10[k]
                  - f_50 * pa_y[k] * is_13[k]
                  + f_51 * pa_y[k] * is_15[k]
                  - f_53 * pa_y[k] * is_17[k]
                  + f_54 * pa_y[k] * is_19[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, hs0_4, hs0_10, hs0_12, hs1_4, hs1_10, hs1_12, is_0, \
                         is_2, is_5, is_6, is_11, is_12, is_13, is_16, is_18, \
                         is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_55 * hs0_4[k]
                  + f_56 * hs0_10[k]
                  - f_57 * hs0_12[k]
                  - f_58 * hs1_4[k]
                  - f_59 * hs1_10[k]
                  + f_60 * hs1_12[k]
                  - 2.1875 * pa_z[k] * is_0[k]
                  - 6.5625 * pa_z[k] * is_2[k]
                  + 13.125 * pa_x[k] * is_5[k]
                  - 6.5625 * pa_z[k] * is_6[k]
                  + 26.25 * pa_x[k] * is_11[k]
                  - 10.5 * pa_x[k] * is_12[k]
                  - 2.1875 * pa_z[k] * is_13[k]
                  + 13.125 * pa_y[k] * is_16[k]
                  - 10.5 * pa_y[k] * is_18[k]
                  + pa_z[k] * is_19[k];
    }

#pragma omp simd aligned(pa_x, hs0_0, hs0_1, hs0_2, hs0_5, hs0_6, hs0_7, hs1_0, hs1_1, hs1_2, \
                         hs1_5, hs1_6, hs1_7, is_0, is_2, is_3, is_6, is_7, is_8, is_13, \
                         is_15, is_17, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_61 * hs0_0[k]
                  - f_62 * hs0_1[k]
                  + f_46 * hs0_2[k]
                  - f_61 * hs0_5[k]
                  + f_46 * hs0_6[k]
                  - f_46 * hs0_7[k]
                  + f_63 * hs1_0[k]
                  + f_64 * hs1_1[k]
                  - f_49 * hs1_2[k]
                  + f_63 * hs1_5[k]
                  - f_49 * hs1_6[k]
                  + f_49 * hs1_7[k]
                  - f_50 * pa_x[k] * is_0[k]
                  - f_52 * pa_x[k] * is_2[k]
                  + f_51 * pa_x[k] * is_3[k]
                  - f_52 * pa_x[k] * is_6[k]
                  + f_53 * pa_x[k] * is_7[k]
                  - f_53 * pa_x[k] * is_8[k]
                  - f_50 * pa_x[k] * is_13[k]
                  + f_51 * pa_x[k] * is_15[k]
                  - f_53 * pa_x[k] * is_17[k]
                  + f_54 * pa_x[k] * is_19[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, hs0_4, hs0_10, hs1_4, hs1_10, is_0, is_2, is_5, \
                         is_6, is_12, is_13, is_16, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += -f_65 * hs0_4[k]
                  + f_65 * hs0_10[k]
                  + f_66 * hs1_4[k]
                  - f_66 * hs1_10[k]
                  + f_67 * pa_z[k] * is_0[k]
                  + f_67 * pa_z[k] * is_2[k]
                  - f_68 * pa_x[k] * is_5[k]
                  - f_67 * pa_z[k] * is_6[k]
                  + f_69 * pa_x[k] * is_12[k]
                  - f_67 * pa_z[k] * is_13[k]
                  + f_68 * pa_y[k] * is_16[k]
                  - f_69 * pa_y[k] * is_18[k];
    }

#pragma omp simd aligned(pa_x, hs0_0, hs0_1, hs0_2, hs0_5, hs0_6, hs0_7, hs1_0, hs1_1, hs1_2, \
                         hs1_5, hs1_6, hs1_7, is_0, is_2, is_3, is_6, is_7, is_8, is_13, \
                         is_15, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += f_70 * hs0_0[k]
                   - f_71 * hs0_1[k]
                   - f_72 * hs0_2[k]
                   - f_73 * hs0_5[k]
                   + f_72 * hs0_6[k]
                   + f_27 * hs0_7[k]
                   - f_74 * hs1_0[k]
                   + f_75 * hs1_1[k]
                   + f_76 * hs1_2[k]
                   + f_77 * hs1_5[k]
                   - f_76 * hs1_6[k]
                   - f_31 * hs1_7[k]
                   + f_36 * pa_x[k] * is_0[k]
                   - f_36 * pa_x[k] * is_2[k]
                   - f_38 * pa_x[k] * is_3[k]
                   - f_34 * pa_x[k] * is_6[k]
                   + f_37 * pa_x[k] * is_7[k]
                   + f_39 * pa_x[k] * is_8[k]
                   - f_32 * pa_x[k] * is_13[k]
                   + f_33 * pa_x[k] * is_15[k]
                   - f_35 * pa_x[k] * is_17[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, hs0_4, hs0_10, hs1_4, hs1_10, is_0, is_2, is_5, \
                         is_6, is_11, is_13, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_78 * hs0_4[k]
                   - f_78 * hs0_10[k]
                   - f_79 * hs1_4[k]
                   + f_79 * hs1_10[k]
                   - f_80 * pa_z[k] * is_0[k]
                   + f_81 * pa_z[k] * is_2[k]
                   + f_82 * pa_x[k] * is_5[k]
                   + f_81 * pa_z[k] * is_6[k]
                   - f_19 * pa_x[k] * is_11[k]
                   - f_80 * pa_z[k] * is_13[k]
                   + f_82 * pa_y[k] * is_16[k];
    }

#pragma omp simd aligned(pa_x, hs0_0, hs0_1, hs0_2, hs0_5, hs0_6, hs1_0, hs1_1, hs1_2, hs1_5, \
                         hs1_6, is_0, is_2, is_3, is_6, is_7, is_13, \
                         is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_12[k] += -f_83 * hs0_0[k]
                   + f_84 * hs0_1[k]
                   + f_85 * hs0_2[k]
                   + f_86 * hs0_5[k]
                   - f_87 * hs0_6[k]
                   + f_88 * hs1_0[k]
                   - f_89 * hs1_1[k]
                   - f_90 * hs1_2[k]
                   - f_91 * hs1_5[k]
                   + f_92 * hs1_6[k]
                   - f_20 * pa_x[k] * is_0[k]
                   + f_18 * pa_x[k] * is_2[k]
                   + f_21 * pa_x[k] * is_3[k]
                   + f_16 * pa_x[k] * is_6[k]
                   - f_19 * pa_x[k] * is_7[k]
                   - f_16 * pa_x[k] * is_13[k]
                   + f_17 * pa_x[k] * is_15[k];
    }

#pragma omp simd aligned(pa_x, pa_z, hs0_0, hs0_1, hs0_5, hs1_0, hs1_1, hs1_5, is_0, is_2, \
                         is_6, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_13[k] += f_93 * pa_z[k] * is_0[k]
                   - f_94 * pa_z[k] * is_2[k]
                   + f_94 * pa_z[k] * is_6[k]
                   - f_93 * pa_z[k] * is_13[k];

        g_14[k] += f_95 * hs0_0[k]
                   - f_96 * hs0_1[k]
                   + f_97 * hs0_5[k]
                   - f_98 * hs1_0[k]
                   + f_99 * hs1_1[k]
                   - f_100 * hs1_5[k]
                   + f_7 * pa_x[k] * is_0[k]
                   - f_6 * pa_x[k] * is_2[k]
                   + f_5 * pa_x[k] * is_6[k]
                   - f_4 * pa_x[k] * is_13[k];
    }
}

}  // namespace simdt2ceri
