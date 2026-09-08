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


#include "SimdElectronRepulsionCtrVrrLS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_ls_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pa, const size_t is0, const size_t is1,
                                    const size_t ks, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 2.625 * std::sqrt(715.0) / alpha;
    const auto f_1 = 1.3125 * std::sqrt(715.0) / alpha;
    const auto f_2 = 2.625 * std::sqrt(715.0) * beta / (alpha * p);
    const auto f_3 = 1.3125 * std::sqrt(715.0) * beta / (alpha * p);
    const auto f_4 = 0.1875 * std::sqrt(715.0);
    const auto f_5 = 1.3125 * std::sqrt(715.0);
    const auto f_6 = 0.65625 * std::sqrt(715.0);
    const auto f_7 = 3.28125 * std::sqrt(715.0);
    const auto f_8 = 1.96875 * std::sqrt(715.0);
    const auto f_9 = 0.09375 * std::sqrt(715.0);
    const auto f_10 = 0.4375 * std::sqrt(858.0) / alpha;
    const auto f_11 = 0.21875 * std::sqrt(858.0) / alpha;
    const auto f_12 = 4.375 * std::sqrt(858.0) / alpha;
    const auto f_13 = 0.4375 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_14 = 0.21875 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_15 = 4.375 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_16 = 0.09375 * std::sqrt(858.0);
    const auto f_17 = 1.3125 * std::sqrt(858.0);
    const auto f_18 = 0.21875 * std::sqrt(858.0);
    const auto f_19 = 4.375 * std::sqrt(858.0);
    const auto f_20 = 1.125 * std::sqrt(1001.0) / alpha;
    const auto f_21 = 1.125 * std::sqrt(1001.0) * beta / (alpha * p);
    const auto f_22 = 0.46875 * std::sqrt(1001.0);
    const auto f_23 = 1.875 * std::sqrt(1001.0);
    const auto f_24 = 0.84375 * std::sqrt(1001.0);
    const auto f_25 = 3.75 * std::sqrt(1001.0);
    const auto f_26 = 0.09375 * std::sqrt(1001.0);
    const auto f_27 = 0.375 * std::sqrt(1001.0);
    const auto f_28 = 0.375 * std::sqrt(77.0) / alpha;
    const auto f_29 = 0.1875 * std::sqrt(77.0) / alpha;
    const auto f_30 = 0.375 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_31 = 0.1875 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_32 = 0.1875 * std::sqrt(77.0);
    const auto f_33 = 4.5 * std::sqrt(77.0);
    const auto f_34 = 7.5 * std::sqrt(77.0);
    const auto f_35 = 0.625 * std::sqrt(1155.0) / alpha;
    const auto f_36 = 0.5 * std::sqrt(1155.0) / alpha;
    const auto f_37 = 0.625 * std::sqrt(1155.0) * beta / (alpha * p);
    const auto f_38 = 0.5 * std::sqrt(1155.0) * beta / (alpha * p);
    const auto f_39 = 0.28125 * std::sqrt(1155.0);
    const auto f_40 = 0.46875 * std::sqrt(1155.0);
    const auto f_41 = 1.875 * std::sqrt(1155.0);
    const auto f_42 = 0.09375 * std::sqrt(1155.0);
    const auto f_43 = 1.5 * std::sqrt(1155.0);
    const auto f_44 = 1.25 * std::sqrt(1155.0);
    const auto f_45 = 0.625 * std::sqrt(1155.0);
    const auto f_46 = 0.5 * std::sqrt(1155.0);
    const auto f_47 = 0.5625 * std::sqrt(70.0) / alpha;
    const auto f_48 = 0.28125 * std::sqrt(70.0) / alpha;
    const auto f_49 = 5.625 * std::sqrt(70.0) / alpha;
    const auto f_50 = 0.5625 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_51 = 0.28125 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_52 = 5.625 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_53 = 0.09375 * std::sqrt(70.0);
    const auto f_54 = 2.8125 * std::sqrt(70.0);
    const auto f_55 = 0.28125 * std::sqrt(70.0);
    const auto f_56 = 7.5 * std::sqrt(70.0);
    const auto f_57 = 5.625 * std::sqrt(70.0);
    const auto f_58 = 3.0 * std::sqrt(70.0);
    const auto f_59 = 78.75 / alpha;
    const auto f_60 = 31.5 / alpha;
    const auto f_61 = 78.75 * beta / (alpha * p);
    const auto f_62 = 31.5 * beta / (alpha * p);
    const auto f_63 = 0.95703125 / alpha;
    const auto f_64 = 2.734375 / alpha;
    const auto f_65 = 21.875 / alpha;
    const auto f_66 = 2.4609375 / alpha;
    const auto f_67 = 39.375 / alpha;
    const auto f_68 = 1.50390625 / alpha;
    const auto f_69 = 35.0 / alpha;
    const auto f_70 = 65.625 / alpha;
    const auto f_71 = 10.5 / alpha;
    const auto f_72 = 0.95703125 * beta / (alpha * p);
    const auto f_73 = 2.734375 * beta / (alpha * p);
    const auto f_74 = 21.875 * beta / (alpha * p);
    const auto f_75 = 2.4609375 * beta / (alpha * p);
    const auto f_76 = 39.375 * beta / (alpha * p);
    const auto f_77 = 1.50390625 * beta / (alpha * p);
    const auto f_78 = 35.0 * beta / (alpha * p);
    const auto f_79 = 65.625 * beta / (alpha * p);
    const auto f_80 = 10.5 * beta / (alpha * p);
    const auto f_81 = 52.5 / alpha;
    const auto f_82 = 52.5 * beta / (alpha * p);
    const auto f_83 = 0.1640625 * std::sqrt(70.0) / alpha;
    const auto f_84 = 0.234375 * std::sqrt(70.0) / alpha;
    const auto f_85 = 3.515625 * std::sqrt(70.0) / alpha;
    const auto f_86 = 2.109375 * std::sqrt(70.0) / alpha;
    const auto f_87 = 0.2109375 * std::sqrt(70.0) / alpha;
    const auto f_88 = 4.21875 * std::sqrt(70.0) / alpha;
    const auto f_89 = 0.1640625 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_90 = 0.234375 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_91 = 3.515625 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_92 = 2.109375 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_93 = 0.2109375 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_94 = 4.21875 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_95 = 0.046875 * std::sqrt(70.0);
    const auto f_96 = 1.40625 * std::sqrt(70.0);
    const auto f_97 = 3.75 * std::sqrt(70.0);
    const auto f_98 = 1.5 * std::sqrt(70.0);
    const auto f_99 = 1.25 * std::sqrt(1155.0) / alpha;
    const auto f_100 = 1.25 * std::sqrt(1155.0) * beta / (alpha * p);
    const auto f_101 = 0.1640625 * std::sqrt(77.0) / alpha;
    const auto f_102 = 0.46875 * std::sqrt(77.0) / alpha;
    const auto f_103 = 2.8125 * std::sqrt(77.0) / alpha;
    const auto f_104 = 0.703125 * std::sqrt(77.0) / alpha;
    const auto f_105 = 8.4375 * std::sqrt(77.0) / alpha;
    const auto f_106 = 0.0703125 * std::sqrt(77.0) / alpha;
    const auto f_107 = 0.1640625 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_108 = 0.46875 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_109 = 2.8125 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_110 = 0.703125 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_111 = 8.4375 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_112 = 0.0703125 * std::sqrt(77.0) * beta / (alpha * p);
    const auto f_113 = 0.046875 * std::sqrt(77.0);
    const auto f_114 = 1.125 * std::sqrt(77.0);
    const auto f_115 = 0.46875 * std::sqrt(77.0);
    const auto f_116 = 5.625 * std::sqrt(77.0);
    const auto f_117 = 1.875 * std::sqrt(77.0);
    const auto f_118 = 11.25 * std::sqrt(77.0);
    const auto f_119 = 0.75 * std::sqrt(1001.0) / alpha;
    const auto f_120 = 3.75 * std::sqrt(1001.0) / alpha;
    const auto f_121 = 0.75 * std::sqrt(1001.0) * beta / (alpha * p);
    const auto f_122 = 3.75 * std::sqrt(1001.0) * beta / (alpha * p);
    const auto f_123 = 0.0546875 * std::sqrt(858.0) / alpha;
    const auto f_124 = 0.546875 * std::sqrt(858.0) / alpha;
    const auto f_125 = 4.921875 * std::sqrt(858.0) / alpha;
    const auto f_126 = 1.09375 * std::sqrt(858.0) / alpha;
    const auto f_127 = 0.0546875 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_128 = 0.546875 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_129 = 4.921875 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_130 = 1.09375 * std::sqrt(858.0) * beta / (alpha * p);
    const auto f_131 = 0.015625 * std::sqrt(858.0);
    const auto f_132 = 3.28125 * std::sqrt(858.0);
    const auto f_133 = 0.08203125 * std::sqrt(715.0) / alpha;
    const auto f_134 = 1.640625 * std::sqrt(715.0) / alpha;
    const auto f_135 = 2.4609375 * std::sqrt(715.0) / alpha;
    const auto f_136 = 0.24609375 * std::sqrt(715.0) / alpha;
    const auto f_137 = 0.08203125 * std::sqrt(715.0) * beta / (alpha * p);
    const auto f_138 = 1.640625 * std::sqrt(715.0) * beta / (alpha * p);
    const auto f_139 = 2.4609375 * std::sqrt(715.0) * beta / (alpha * p);
    const auto f_140 = 0.24609375 * std::sqrt(715.0) * beta / (alpha * p);
    const auto f_141 = 0.0234375 * std::sqrt(715.0);
    const auto f_142 = 1.640625 * std::sqrt(715.0);

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
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_13 = buffer.data(is0 + 13);
    const auto *is0_14 = buffer.data(is0 + 14);
    const auto *is0_15 = buffer.data(is0 + 15);
    const auto *is0_16 = buffer.data(is0 + 16);
    const auto *is0_17 = buffer.data(is0 + 17);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_13 = buffer.data(is1 + 13);
    const auto *is1_14 = buffer.data(is1 + 14);
    const auto *is1_15 = buffer.data(is1 + 15);
    const auto *is1_16 = buffer.data(is1 + 16);
    const auto *is1_17 = buffer.data(is1 + 17);

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

#pragma omp simd aligned(pa_x, pa_y, pa_z, is0_3, is0_8, is1_3, is1_8, ks_0, ks_1, ks_4, ks_9, \
                         ks_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * is0_3[k]
                  + f_1 * is0_8[k]
                  + f_2 * is1_3[k]
                  - f_3 * is1_8[k]
                  + f_4 * pa_y[k] * ks_0[k]
                  - f_5 * pa_x[k] * ks_4[k]
                  + f_5 * pa_x[k] * ks_9[k]
                  - f_4 * pa_x[k] * ks_18[k];

        g_1[k] += f_6 * pa_y[k] * ks_1[k]
                  - f_7 * pa_z[k] * ks_4[k]
                  + f_8 * pa_z[k] * ks_9[k]
                  - f_9 * pa_z[k] * ks_18[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_3, is0_8, is0_9, is1_3, is1_8, is1_9, ks_0, ks_3, \
                         ks_4, ks_9, ks_10, ks_18, ks_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += f_10 * is0_3[k]
                  + f_11 * is0_8[k]
                  - f_12 * is0_9[k]
                  - f_13 * is1_3[k]
                  - f_14 * is1_8[k]
                  + f_15 * is1_9[k]
                  - f_16 * pa_y[k] * ks_0[k]
                  + f_17 * pa_y[k] * ks_3[k]
                  + f_18 * pa_x[k] * ks_4[k]
                  + f_18 * pa_x[k] * ks_9[k]
                  - f_19 * pa_x[k] * ks_10[k]
                  - f_16 * pa_x[k] * ks_18[k]
                  + f_17 * pa_x[k] * ks_20[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, is0_14, is1_14, ks_1, ks_4, ks_5, ks_9, ks_15, \
                         ks_18, ks_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_20 * is0_14[k]
                  + f_21 * is1_14[k]
                  - f_22 * pa_y[k] * ks_1[k]
                  + f_22 * pa_z[k] * ks_4[k]
                  + f_23 * pa_y[k] * ks_5[k]
                  + f_24 * pa_z[k] * ks_9[k]
                  - f_25 * pa_x[k] * ks_15[k]
                  - f_26 * pa_z[k] * ks_18[k]
                  + f_27 * pa_y[k] * ks_21[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_3, is0_8, is1_3, is1_8, ks_0, ks_3, ks_4, ks_8, ks_9, \
                         ks_18, ks_20, ks_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_28 * is0_3[k]
                  - f_29 * is0_8[k]
                  - f_30 * is1_3[k]
                  + f_31 * is1_8[k]
                  + f_32 * pa_y[k] * ks_0[k]
                  - f_33 * pa_y[k] * ks_3[k]
                  + f_32 * pa_x[k] * ks_4[k]
                  + f_34 * pa_y[k] * ks_8[k]
                  - f_32 * pa_x[k] * ks_9[k]
                  - f_32 * pa_x[k] * ks_18[k]
                  + f_33 * pa_x[k] * ks_20[k]
                  - f_34 * pa_x[k] * ks_22[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, is0_14, is0_16, is1_14, is1_16, ks_1, ks_4, ks_5, \
                         ks_9, ks_12, ks_15, ks_18, ks_21, ks_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_35 * is0_14[k]
                  - f_36 * is0_16[k]
                  - f_37 * is1_14[k]
                  + f_38 * is1_16[k]
                  + f_39 * pa_y[k] * ks_1[k]
                  + f_40 * pa_z[k] * ks_4[k]
                  - f_41 * pa_y[k] * ks_5[k]
                  + f_42 * pa_z[k] * ks_9[k]
                  + f_43 * pa_y[k] * ks_12[k]
                  - f_44 * pa_x[k] * ks_15[k]
                  - f_42 * pa_z[k] * ks_18[k]
                  + f_45 * pa_y[k] * ks_21[k]
                  - f_46 * pa_y[k] * ks_23[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_3, is0_8, is0_9, is1_3, is1_8, is1_9, ks_0, ks_3, \
                         ks_4, ks_8, ks_9, ks_10, ks_18, ks_20, ks_22, \
                         ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_47 * is0_3[k]
                  - f_48 * is0_8[k]
                  + f_49 * is0_9[k]
                  + f_50 * is1_3[k]
                  + f_51 * is1_8[k]
                  - f_52 * is1_9[k]
                  - f_53 * pa_y[k] * ks_0[k]
                  + f_54 * pa_y[k] * ks_3[k]
                  - f_55 * pa_x[k] * ks_4[k]
                  - f_56 * pa_y[k] * ks_8[k]
                  - f_55 * pa_x[k] * ks_9[k]
                  + f_57 * pa_x[k] * ks_10[k]
                  - f_53 * pa_x[k] * ks_18[k]
                  + f_54 * pa_x[k] * ks_20[k]
                  - f_56 * pa_x[k] * ks_22[k]
                  + f_58 * pa_x[k] * ks_24[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, is0_14, is0_16, is1_14, is1_16, ks_1, ks_4, ks_5, \
                         ks_9, ks_12, ks_15, ks_18, ks_21, ks_23, \
                         ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_59 * is0_14[k]
                  - f_60 * is0_16[k]
                  - f_61 * is1_14[k]
                  + f_62 * is1_16[k]
                  - 3.28125 * pa_y[k] * ks_1[k]
                  - 9.84375 * pa_z[k] * ks_4[k]
                  + 26.25 * pa_y[k] * ks_5[k]
                  - 9.84375 * pa_z[k] * ks_9[k]
                  - 31.5 * pa_y[k] * ks_12[k]
                  + 52.5 * pa_x[k] * ks_15[k]
                  - 3.28125 * pa_z[k] * ks_18[k]
                  + 26.25 * pa_y[k] * ks_21[k]
                  - 31.5 * pa_y[k] * ks_23[k]
                  + 6.0 * pa_y[k] * ks_25[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, is0_0, is0_1, is0_2, is0_5, is0_6, is0_7, is0_12, \
                         is0_13, is0_15, is0_17, is1_0, is1_1, is1_2, is1_5, is1_6, is1_7, \
                         is1_12, is1_13, is1_15, is1_17, ks_0, ks_2, ks_3, ks_6, ks_7, ks_8, \
                         ks_13, ks_14, ks_16, ks_17, ks_18, ks_20, ks_22, ks_24, \
                         ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_63 * is0_0[k]
                  + f_64 * is0_1[k]
                  - f_65 * is0_2[k]
                  + f_66 * is0_5[k]
                  - f_67 * is0_6[k]
                  + f_67 * is0_7[k]
                  + f_68 * is0_12[k]
                  - f_69 * is0_13[k]
                  + f_70 * is0_15[k]
                  - f_71 * is0_17[k]
                  - f_72 * is1_0[k]
                  - f_73 * is1_1[k]
                  + f_74 * is1_2[k]
                  - f_75 * is1_5[k]
                  + f_76 * is1_6[k]
                  - f_76 * is1_7[k]
                  - f_77 * is1_12[k]
                  + f_78 * is1_13[k]
                  - f_79 * is1_15[k]
                  + f_80 * is1_17[k]
                  + 0.2734375 * pa_x[k] * ks_0[k]
                  + 1.09375 * pa_x[k] * ks_2[k]
                  - 8.75 * pa_x[k] * ks_3[k]
                  + 1.640625 * pa_x[k] * ks_6[k]
                  - 26.25 * pa_x[k] * ks_7[k]
                  + 26.25 * pa_x[k] * ks_8[k]
                  + 1.09375 * pa_x[k] * ks_13[k]
                  - 26.25 * pa_x[k] * ks_14[k]
                  + 52.5 * pa_x[k] * ks_16[k]
                  - 14.0 * pa_x[k] * ks_17[k]
                  + 0.2734375 * pa_y[k] * ks_18[k]
                  - 8.75 * pa_y[k] * ks_20[k]
                  + 26.25 * pa_y[k] * ks_22[k]
                  - 14.0 * pa_y[k] * ks_24[k]
                  + pa_z[k] * ks_25[k];
    }

#pragma omp simd aligned(pa_x, pa_z, is0_4, is0_10, is0_11, is1_4, is1_10, is1_11, ks_0, ks_2, \
                         ks_5, ks_6, ks_11, ks_12, ks_19, ks_21, ks_23, \
                         ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_81 * is0_4[k]
                  + f_81 * is0_10[k]
                  - f_60 * is0_11[k]
                  - f_82 * is1_4[k]
                  - f_82 * is1_10[k]
                  + f_62 * is1_11[k]
                  - 3.28125 * pa_z[k] * ks_0[k]
                  - 9.84375 * pa_z[k] * ks_2[k]
                  + 26.25 * pa_x[k] * ks_5[k]
                  - 9.84375 * pa_z[k] * ks_6[k]
                  + 52.5 * pa_x[k] * ks_11[k]
                  - 31.5 * pa_x[k] * ks_12[k]
                  - 3.28125 * pa_x[k] * ks_19[k]
                  + 26.25 * pa_x[k] * ks_21[k]
                  - 31.5 * pa_x[k] * ks_23[k]
                  + 6.0 * pa_x[k] * ks_25[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_0, is0_1, is0_2, is0_6, is0_7, is0_12, is0_13, \
                         is0_15, is1_0, is1_1, is1_2, is1_6, is1_7, is1_12, is1_13, is1_15, \
                         ks_0, ks_2, ks_3, ks_7, ks_8, ks_13, ks_14, ks_17, ks_18, ks_20, \
                         ks_22, ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_83 * is0_0[k]
                   - f_84 * is0_1[k]
                   + f_85 * is0_2[k]
                   + f_86 * is0_6[k]
                   - f_49 * is0_7[k]
                   + f_87 * is0_12[k]
                   - f_88 * is0_13[k]
                   + f_49 * is0_15[k]
                   + f_89 * is1_0[k]
                   + f_90 * is1_1[k]
                   - f_91 * is1_2[k]
                   - f_92 * is1_6[k]
                   + f_52 * is1_7[k]
                   - f_93 * is1_12[k]
                   + f_94 * is1_13[k]
                   - f_52 * is1_15[k]
                   - f_95 * pa_x[k] * ks_0[k]
                   - f_53 * pa_x[k] * ks_2[k]
                   + f_96 * pa_x[k] * ks_3[k]
                   + f_96 * pa_x[k] * ks_7[k]
                   - f_97 * pa_x[k] * ks_8[k]
                   + f_53 * pa_x[k] * ks_13[k]
                   - f_96 * pa_x[k] * ks_14[k]
                   + f_98 * pa_x[k] * ks_17[k]
                   + f_95 * pa_y[k] * ks_18[k]
                   - f_96 * pa_y[k] * ks_20[k]
                   + f_97 * pa_y[k] * ks_22[k]
                   - f_98 * pa_y[k] * ks_24[k];
    }

#pragma omp simd aligned(pa_x, pa_z, is0_4, is0_10, is0_11, is1_4, is1_10, is1_11, ks_0, ks_2, \
                         ks_5, ks_6, ks_11, ks_12, ks_19, ks_21, \
                         ks_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += -f_99 * is0_4[k]
                   + f_99 * is0_10[k]
                   + f_36 * is0_11[k]
                   + f_100 * is1_4[k]
                   - f_100 * is1_10[k]
                   - f_38 * is1_11[k]
                   + f_42 * pa_z[k] * ks_0[k]
                   - f_42 * pa_z[k] * ks_2[k]
                   - f_45 * pa_x[k] * ks_5[k]
                   - f_40 * pa_z[k] * ks_6[k]
                   + f_44 * pa_x[k] * ks_11[k]
                   + f_46 * pa_x[k] * ks_12[k]
                   - f_39 * pa_x[k] * ks_19[k]
                   + f_41 * pa_x[k] * ks_21[k]
                   - f_43 * pa_x[k] * ks_23[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_0, is0_1, is0_2, is0_5, is0_6, is0_7, is0_12, is0_15, \
                         is1_0, is1_1, is1_2, is1_5, is1_6, is1_7, is1_12, is1_15, ks_0, ks_2, \
                         ks_3, ks_6, ks_7, ks_8, ks_13, ks_14, ks_16, ks_18, ks_20, \
                         ks_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_12[k] += f_101 * is0_0[k]
                   - f_102 * is0_1[k]
                   - f_103 * is0_2[k]
                   - f_104 * is0_5[k]
                   + f_105 * is0_6[k]
                   + f_103 * is0_7[k]
                   + f_106 * is0_12[k]
                   - f_103 * is0_15[k]
                   - f_107 * is1_0[k]
                   + f_108 * is1_1[k]
                   + f_109 * is1_2[k]
                   + f_110 * is1_5[k]
                   - f_111 * is1_6[k]
                   - f_109 * is1_7[k]
                   - f_112 * is1_12[k]
                   + f_109 * is1_15[k]
                   + f_113 * pa_x[k] * ks_0[k]
                   - f_32 * pa_x[k] * ks_2[k]
                   - f_114 * pa_x[k] * ks_3[k]
                   - f_115 * pa_x[k] * ks_6[k]
                   + f_116 * pa_x[k] * ks_7[k]
                   + f_117 * pa_x[k] * ks_8[k]
                   - f_32 * pa_x[k] * ks_13[k]
                   + f_116 * pa_x[k] * ks_14[k]
                   - f_118 * pa_x[k] * ks_16[k]
                   + f_113 * pa_y[k] * ks_18[k]
                   - f_114 * pa_y[k] * ks_20[k]
                   + f_117 * pa_y[k] * ks_22[k];
    }

#pragma omp simd aligned(pa_x, pa_z, is0_4, is0_10, is1_4, is1_10, ks_0, ks_2, ks_5, ks_6, \
                         ks_11, ks_19, ks_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_13[k] += f_119 * is0_4[k]
                   - f_120 * is0_10[k]
                   - f_121 * is1_4[k]
                   + f_122 * is1_10[k]
                   - f_26 * pa_z[k] * ks_0[k]
                   + f_24 * pa_z[k] * ks_2[k]
                   + f_27 * pa_x[k] * ks_5[k]
                   + f_22 * pa_z[k] * ks_6[k]
                   - f_25 * pa_x[k] * ks_11[k]
                   - f_22 * pa_x[k] * ks_19[k]
                   + f_23 * pa_x[k] * ks_21[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_0, is0_1, is0_2, is0_6, is0_12, is0_13, is1_0, is1_1, \
                         is1_2, is1_6, is1_12, is1_13, ks_0, ks_2, ks_3, ks_7, ks_13, ks_14, \
                         ks_18, ks_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_14[k] += -f_123 * is0_0[k]
                   + f_124 * is0_1[k]
                   + f_124 * is0_2[k]
                   - f_125 * is0_6[k]
                   - f_123 * is0_12[k]
                   + f_126 * is0_13[k]
                   + f_127 * is1_0[k]
                   - f_128 * is1_1[k]
                   - f_128 * is1_2[k]
                   + f_129 * is1_6[k]
                   + f_127 * is1_12[k]
                   - f_130 * is1_13[k]
                   - f_131 * pa_x[k] * ks_0[k]
                   + f_18 * pa_x[k] * ks_2[k]
                   + f_18 * pa_x[k] * ks_3[k]
                   - f_132 * pa_x[k] * ks_7[k]
                   - f_18 * pa_x[k] * ks_13[k]
                   + f_132 * pa_x[k] * ks_14[k]
                   + f_131 * pa_y[k] * ks_18[k]
                   - f_18 * pa_y[k] * ks_20[k];
    }

#pragma omp simd aligned(pa_x, pa_z, ks_0, ks_2, ks_6, ks_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_15[k] += f_9 * pa_z[k] * ks_0[k]
                   - f_8 * pa_z[k] * ks_2[k]
                   + f_7 * pa_z[k] * ks_6[k]
                   - f_6 * pa_x[k] * ks_19[k];
    }

#pragma omp simd aligned(pa_x, pa_y, is0_0, is0_1, is0_5, is0_12, is1_0, is1_1, is1_5, is1_12, \
                         ks_0, ks_2, ks_6, ks_13, ks_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_16[k] += f_133 * is0_0[k]
                   - f_134 * is0_1[k]
                   + f_135 * is0_5[k]
                   - f_136 * is0_12[k]
                   - f_137 * is1_0[k]
                   + f_138 * is1_1[k]
                   - f_139 * is1_5[k]
                   + f_140 * is1_12[k]
                   + f_141 * pa_x[k] * ks_0[k]
                   - f_6 * pa_x[k] * ks_2[k]
                   + f_142 * pa_x[k] * ks_6[k]
                   - f_6 * pa_x[k] * ks_13[k]
                   + f_141 * pa_y[k] * ks_18[k];
    }
}

}  // namespace simdt2ceri
