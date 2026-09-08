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


#include "SimdElectronRepulsionCtrVrrSL.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_sl_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pb, const size_t si0, const size_t si1,
                                    const size_t sk, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 2.625 * std::sqrt(715.0) / beta;
    const auto f_1 = 1.3125 * std::sqrt(715.0) / beta;
    const auto f_2 = 2.625 * std::sqrt(715.0) * alpha / (beta * p);
    const auto f_3 = 1.3125 * std::sqrt(715.0) * alpha / (beta * p);
    const auto f_4 = 0.1875 * std::sqrt(715.0);
    const auto f_5 = 1.3125 * std::sqrt(715.0);
    const auto f_6 = 0.65625 * std::sqrt(715.0);
    const auto f_7 = 3.28125 * std::sqrt(715.0);
    const auto f_8 = 1.96875 * std::sqrt(715.0);
    const auto f_9 = 0.09375 * std::sqrt(715.0);
    const auto f_10 = 0.4375 * std::sqrt(858.0) / beta;
    const auto f_11 = 0.21875 * std::sqrt(858.0) / beta;
    const auto f_12 = 4.375 * std::sqrt(858.0) / beta;
    const auto f_13 = 0.4375 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_14 = 0.21875 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_15 = 4.375 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_16 = 0.09375 * std::sqrt(858.0);
    const auto f_17 = 1.3125 * std::sqrt(858.0);
    const auto f_18 = 0.21875 * std::sqrt(858.0);
    const auto f_19 = 4.375 * std::sqrt(858.0);
    const auto f_20 = 1.125 * std::sqrt(1001.0) / beta;
    const auto f_21 = 1.125 * std::sqrt(1001.0) * alpha / (beta * p);
    const auto f_22 = 0.46875 * std::sqrt(1001.0);
    const auto f_23 = 1.875 * std::sqrt(1001.0);
    const auto f_24 = 0.84375 * std::sqrt(1001.0);
    const auto f_25 = 3.75 * std::sqrt(1001.0);
    const auto f_26 = 0.09375 * std::sqrt(1001.0);
    const auto f_27 = 0.375 * std::sqrt(1001.0);
    const auto f_28 = 0.375 * std::sqrt(77.0) / beta;
    const auto f_29 = 0.1875 * std::sqrt(77.0) / beta;
    const auto f_30 = 0.375 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_31 = 0.1875 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_32 = 0.1875 * std::sqrt(77.0);
    const auto f_33 = 4.5 * std::sqrt(77.0);
    const auto f_34 = 7.5 * std::sqrt(77.0);
    const auto f_35 = 0.625 * std::sqrt(1155.0) / beta;
    const auto f_36 = 0.5 * std::sqrt(1155.0) / beta;
    const auto f_37 = 0.625 * std::sqrt(1155.0) * alpha / (beta * p);
    const auto f_38 = 0.5 * std::sqrt(1155.0) * alpha / (beta * p);
    const auto f_39 = 0.28125 * std::sqrt(1155.0);
    const auto f_40 = 0.46875 * std::sqrt(1155.0);
    const auto f_41 = 1.875 * std::sqrt(1155.0);
    const auto f_42 = 0.09375 * std::sqrt(1155.0);
    const auto f_43 = 1.5 * std::sqrt(1155.0);
    const auto f_44 = 1.25 * std::sqrt(1155.0);
    const auto f_45 = 0.625 * std::sqrt(1155.0);
    const auto f_46 = 0.5 * std::sqrt(1155.0);
    const auto f_47 = 0.5625 * std::sqrt(70.0) / beta;
    const auto f_48 = 0.28125 * std::sqrt(70.0) / beta;
    const auto f_49 = 5.625 * std::sqrt(70.0) / beta;
    const auto f_50 = 0.5625 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_51 = 0.28125 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_52 = 5.625 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_53 = 0.09375 * std::sqrt(70.0);
    const auto f_54 = 2.8125 * std::sqrt(70.0);
    const auto f_55 = 0.28125 * std::sqrt(70.0);
    const auto f_56 = 7.5 * std::sqrt(70.0);
    const auto f_57 = 5.625 * std::sqrt(70.0);
    const auto f_58 = 3.0 * std::sqrt(70.0);
    const auto f_59 = 78.75 / beta;
    const auto f_60 = 31.5 / beta;
    const auto f_61 = 78.75 * alpha / (beta * p);
    const auto f_62 = 31.5 * alpha / (beta * p);
    const auto f_63 = 0.95703125 / beta;
    const auto f_64 = 2.734375 / beta;
    const auto f_65 = 21.875 / beta;
    const auto f_66 = 2.4609375 / beta;
    const auto f_67 = 39.375 / beta;
    const auto f_68 = 1.50390625 / beta;
    const auto f_69 = 35.0 / beta;
    const auto f_70 = 65.625 / beta;
    const auto f_71 = 10.5 / beta;
    const auto f_72 = 0.95703125 * alpha / (beta * p);
    const auto f_73 = 2.734375 * alpha / (beta * p);
    const auto f_74 = 21.875 * alpha / (beta * p);
    const auto f_75 = 2.4609375 * alpha / (beta * p);
    const auto f_76 = 39.375 * alpha / (beta * p);
    const auto f_77 = 1.50390625 * alpha / (beta * p);
    const auto f_78 = 35.0 * alpha / (beta * p);
    const auto f_79 = 65.625 * alpha / (beta * p);
    const auto f_80 = 10.5 * alpha / (beta * p);
    const auto f_81 = 52.5 / beta;
    const auto f_82 = 52.5 * alpha / (beta * p);
    const auto f_83 = 0.1640625 * std::sqrt(70.0) / beta;
    const auto f_84 = 0.234375 * std::sqrt(70.0) / beta;
    const auto f_85 = 3.515625 * std::sqrt(70.0) / beta;
    const auto f_86 = 2.109375 * std::sqrt(70.0) / beta;
    const auto f_87 = 0.2109375 * std::sqrt(70.0) / beta;
    const auto f_88 = 4.21875 * std::sqrt(70.0) / beta;
    const auto f_89 = 0.1640625 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_90 = 0.234375 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_91 = 3.515625 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_92 = 2.109375 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_93 = 0.2109375 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_94 = 4.21875 * std::sqrt(70.0) * alpha / (beta * p);
    const auto f_95 = 0.046875 * std::sqrt(70.0);
    const auto f_96 = 1.40625 * std::sqrt(70.0);
    const auto f_97 = 3.75 * std::sqrt(70.0);
    const auto f_98 = 1.5 * std::sqrt(70.0);
    const auto f_99 = 1.25 * std::sqrt(1155.0) / beta;
    const auto f_100 = 1.25 * std::sqrt(1155.0) * alpha / (beta * p);
    const auto f_101 = 0.1640625 * std::sqrt(77.0) / beta;
    const auto f_102 = 0.46875 * std::sqrt(77.0) / beta;
    const auto f_103 = 2.8125 * std::sqrt(77.0) / beta;
    const auto f_104 = 0.703125 * std::sqrt(77.0) / beta;
    const auto f_105 = 8.4375 * std::sqrt(77.0) / beta;
    const auto f_106 = 0.0703125 * std::sqrt(77.0) / beta;
    const auto f_107 = 0.1640625 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_108 = 0.46875 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_109 = 2.8125 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_110 = 0.703125 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_111 = 8.4375 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_112 = 0.0703125 * std::sqrt(77.0) * alpha / (beta * p);
    const auto f_113 = 0.046875 * std::sqrt(77.0);
    const auto f_114 = 1.125 * std::sqrt(77.0);
    const auto f_115 = 0.46875 * std::sqrt(77.0);
    const auto f_116 = 5.625 * std::sqrt(77.0);
    const auto f_117 = 1.875 * std::sqrt(77.0);
    const auto f_118 = 11.25 * std::sqrt(77.0);
    const auto f_119 = 0.75 * std::sqrt(1001.0) / beta;
    const auto f_120 = 3.75 * std::sqrt(1001.0) / beta;
    const auto f_121 = 0.75 * std::sqrt(1001.0) * alpha / (beta * p);
    const auto f_122 = 3.75 * std::sqrt(1001.0) * alpha / (beta * p);
    const auto f_123 = 0.0546875 * std::sqrt(858.0) / beta;
    const auto f_124 = 0.546875 * std::sqrt(858.0) / beta;
    const auto f_125 = 4.921875 * std::sqrt(858.0) / beta;
    const auto f_126 = 1.09375 * std::sqrt(858.0) / beta;
    const auto f_127 = 0.0546875 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_128 = 0.546875 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_129 = 4.921875 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_130 = 1.09375 * std::sqrt(858.0) * alpha / (beta * p);
    const auto f_131 = 0.015625 * std::sqrt(858.0);
    const auto f_132 = 3.28125 * std::sqrt(858.0);
    const auto f_133 = 0.08203125 * std::sqrt(715.0) / beta;
    const auto f_134 = 1.640625 * std::sqrt(715.0) / beta;
    const auto f_135 = 2.4609375 * std::sqrt(715.0) / beta;
    const auto f_136 = 0.24609375 * std::sqrt(715.0) / beta;
    const auto f_137 = 0.08203125 * std::sqrt(715.0) * alpha / (beta * p);
    const auto f_138 = 1.640625 * std::sqrt(715.0) * alpha / (beta * p);
    const auto f_139 = 2.4609375 * std::sqrt(715.0) * alpha / (beta * p);
    const auto f_140 = 0.24609375 * std::sqrt(715.0) * alpha / (beta * p);
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si0_0 = buffer.data(si0 + 0);
    const auto *si0_1 = buffer.data(si0 + 1);
    const auto *si0_2 = buffer.data(si0 + 2);
    const auto *si0_3 = buffer.data(si0 + 3);
    const auto *si0_4 = buffer.data(si0 + 4);
    const auto *si0_5 = buffer.data(si0 + 5);
    const auto *si0_6 = buffer.data(si0 + 6);
    const auto *si0_7 = buffer.data(si0 + 7);
    const auto *si0_8 = buffer.data(si0 + 8);
    const auto *si0_9 = buffer.data(si0 + 9);
    const auto *si0_10 = buffer.data(si0 + 10);
    const auto *si0_11 = buffer.data(si0 + 11);
    const auto *si0_12 = buffer.data(si0 + 12);
    const auto *si0_13 = buffer.data(si0 + 13);
    const auto *si0_14 = buffer.data(si0 + 14);
    const auto *si0_15 = buffer.data(si0 + 15);
    const auto *si0_16 = buffer.data(si0 + 16);
    const auto *si0_17 = buffer.data(si0 + 17);

    const auto *si1_0 = buffer.data(si1 + 0);
    const auto *si1_1 = buffer.data(si1 + 1);
    const auto *si1_2 = buffer.data(si1 + 2);
    const auto *si1_3 = buffer.data(si1 + 3);
    const auto *si1_4 = buffer.data(si1 + 4);
    const auto *si1_5 = buffer.data(si1 + 5);
    const auto *si1_6 = buffer.data(si1 + 6);
    const auto *si1_7 = buffer.data(si1 + 7);
    const auto *si1_8 = buffer.data(si1 + 8);
    const auto *si1_9 = buffer.data(si1 + 9);
    const auto *si1_10 = buffer.data(si1 + 10);
    const auto *si1_11 = buffer.data(si1 + 11);
    const auto *si1_12 = buffer.data(si1 + 12);
    const auto *si1_13 = buffer.data(si1 + 13);
    const auto *si1_14 = buffer.data(si1 + 14);
    const auto *si1_15 = buffer.data(si1 + 15);
    const auto *si1_16 = buffer.data(si1 + 16);
    const auto *si1_17 = buffer.data(si1 + 17);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_6 = buffer.data(sk + 6);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_13 = buffer.data(sk + 13);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_16 = buffer.data(sk + 16);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);

#pragma omp simd aligned(pb_x, pb_y, pb_z, si0_3, si0_8, si1_3, si1_8, sk_0, sk_1, sk_4, sk_9, \
                         sk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * si0_3[k]
                  + f_1 * si0_8[k]
                  + f_2 * si1_3[k]
                  - f_3 * si1_8[k]
                  + f_4 * pb_y[k] * sk_0[k]
                  - f_5 * pb_x[k] * sk_4[k]
                  + f_5 * pb_x[k] * sk_9[k]
                  - f_4 * pb_x[k] * sk_18[k];

        g_1[k] += f_6 * pb_y[k] * sk_1[k]
                  - f_7 * pb_z[k] * sk_4[k]
                  + f_8 * pb_z[k] * sk_9[k]
                  - f_9 * pb_z[k] * sk_18[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_3, si0_8, si0_9, si1_3, si1_8, si1_9, sk_0, sk_3, \
                         sk_4, sk_9, sk_10, sk_18, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += f_10 * si0_3[k]
                  + f_11 * si0_8[k]
                  - f_12 * si0_9[k]
                  - f_13 * si1_3[k]
                  - f_14 * si1_8[k]
                  + f_15 * si1_9[k]
                  - f_16 * pb_y[k] * sk_0[k]
                  + f_17 * pb_y[k] * sk_3[k]
                  + f_18 * pb_x[k] * sk_4[k]
                  + f_18 * pb_x[k] * sk_9[k]
                  - f_19 * pb_x[k] * sk_10[k]
                  - f_16 * pb_x[k] * sk_18[k]
                  + f_17 * pb_x[k] * sk_20[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, si0_14, si1_14, sk_1, sk_4, sk_5, sk_9, sk_15, \
                         sk_18, sk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_20 * si0_14[k]
                  + f_21 * si1_14[k]
                  - f_22 * pb_y[k] * sk_1[k]
                  + f_22 * pb_z[k] * sk_4[k]
                  + f_23 * pb_y[k] * sk_5[k]
                  + f_24 * pb_z[k] * sk_9[k]
                  - f_25 * pb_x[k] * sk_15[k]
                  - f_26 * pb_z[k] * sk_18[k]
                  + f_27 * pb_y[k] * sk_21[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_3, si0_8, si1_3, si1_8, sk_0, sk_3, sk_4, sk_8, sk_9, \
                         sk_18, sk_20, sk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_28 * si0_3[k]
                  - f_29 * si0_8[k]
                  - f_30 * si1_3[k]
                  + f_31 * si1_8[k]
                  + f_32 * pb_y[k] * sk_0[k]
                  - f_33 * pb_y[k] * sk_3[k]
                  + f_32 * pb_x[k] * sk_4[k]
                  + f_34 * pb_y[k] * sk_8[k]
                  - f_32 * pb_x[k] * sk_9[k]
                  - f_32 * pb_x[k] * sk_18[k]
                  + f_33 * pb_x[k] * sk_20[k]
                  - f_34 * pb_x[k] * sk_22[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, si0_14, si0_16, si1_14, si1_16, sk_1, sk_4, sk_5, \
                         sk_9, sk_12, sk_15, sk_18, sk_21, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_35 * si0_14[k]
                  - f_36 * si0_16[k]
                  - f_37 * si1_14[k]
                  + f_38 * si1_16[k]
                  + f_39 * pb_y[k] * sk_1[k]
                  + f_40 * pb_z[k] * sk_4[k]
                  - f_41 * pb_y[k] * sk_5[k]
                  + f_42 * pb_z[k] * sk_9[k]
                  + f_43 * pb_y[k] * sk_12[k]
                  - f_44 * pb_x[k] * sk_15[k]
                  - f_42 * pb_z[k] * sk_18[k]
                  + f_45 * pb_y[k] * sk_21[k]
                  - f_46 * pb_y[k] * sk_23[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_3, si0_8, si0_9, si1_3, si1_8, si1_9, sk_0, sk_3, \
                         sk_4, sk_8, sk_9, sk_10, sk_18, sk_20, sk_22, \
                         sk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_47 * si0_3[k]
                  - f_48 * si0_8[k]
                  + f_49 * si0_9[k]
                  + f_50 * si1_3[k]
                  + f_51 * si1_8[k]
                  - f_52 * si1_9[k]
                  - f_53 * pb_y[k] * sk_0[k]
                  + f_54 * pb_y[k] * sk_3[k]
                  - f_55 * pb_x[k] * sk_4[k]
                  - f_56 * pb_y[k] * sk_8[k]
                  - f_55 * pb_x[k] * sk_9[k]
                  + f_57 * pb_x[k] * sk_10[k]
                  - f_53 * pb_x[k] * sk_18[k]
                  + f_54 * pb_x[k] * sk_20[k]
                  - f_56 * pb_x[k] * sk_22[k]
                  + f_58 * pb_x[k] * sk_24[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, si0_14, si0_16, si1_14, si1_16, sk_1, sk_4, sk_5, \
                         sk_9, sk_12, sk_15, sk_18, sk_21, sk_23, \
                         sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_59 * si0_14[k]
                  - f_60 * si0_16[k]
                  - f_61 * si1_14[k]
                  + f_62 * si1_16[k]
                  - 3.28125 * pb_y[k] * sk_1[k]
                  - 9.84375 * pb_z[k] * sk_4[k]
                  + 26.25 * pb_y[k] * sk_5[k]
                  - 9.84375 * pb_z[k] * sk_9[k]
                  - 31.5 * pb_y[k] * sk_12[k]
                  + 52.5 * pb_x[k] * sk_15[k]
                  - 3.28125 * pb_z[k] * sk_18[k]
                  + 26.25 * pb_y[k] * sk_21[k]
                  - 31.5 * pb_y[k] * sk_23[k]
                  + 6.0 * pb_y[k] * sk_25[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, si0_0, si0_1, si0_2, si0_5, si0_6, si0_7, si0_12, \
                         si0_13, si0_15, si0_17, si1_0, si1_1, si1_2, si1_5, si1_6, si1_7, \
                         si1_12, si1_13, si1_15, si1_17, sk_0, sk_2, sk_3, sk_6, sk_7, sk_8, \
                         sk_13, sk_14, sk_16, sk_17, sk_18, sk_20, sk_22, sk_24, \
                         sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_63 * si0_0[k]
                  + f_64 * si0_1[k]
                  - f_65 * si0_2[k]
                  + f_66 * si0_5[k]
                  - f_67 * si0_6[k]
                  + f_67 * si0_7[k]
                  + f_68 * si0_12[k]
                  - f_69 * si0_13[k]
                  + f_70 * si0_15[k]
                  - f_71 * si0_17[k]
                  - f_72 * si1_0[k]
                  - f_73 * si1_1[k]
                  + f_74 * si1_2[k]
                  - f_75 * si1_5[k]
                  + f_76 * si1_6[k]
                  - f_76 * si1_7[k]
                  - f_77 * si1_12[k]
                  + f_78 * si1_13[k]
                  - f_79 * si1_15[k]
                  + f_80 * si1_17[k]
                  + 0.2734375 * pb_x[k] * sk_0[k]
                  + 1.09375 * pb_x[k] * sk_2[k]
                  - 8.75 * pb_x[k] * sk_3[k]
                  + 1.640625 * pb_x[k] * sk_6[k]
                  - 26.25 * pb_x[k] * sk_7[k]
                  + 26.25 * pb_x[k] * sk_8[k]
                  + 1.09375 * pb_x[k] * sk_13[k]
                  - 26.25 * pb_x[k] * sk_14[k]
                  + 52.5 * pb_x[k] * sk_16[k]
                  - 14.0 * pb_x[k] * sk_17[k]
                  + 0.2734375 * pb_y[k] * sk_18[k]
                  - 8.75 * pb_y[k] * sk_20[k]
                  + 26.25 * pb_y[k] * sk_22[k]
                  - 14.0 * pb_y[k] * sk_24[k]
                  + pb_z[k] * sk_25[k];
    }

#pragma omp simd aligned(pb_x, pb_z, si0_4, si0_10, si0_11, si1_4, si1_10, si1_11, sk_0, sk_2, \
                         sk_5, sk_6, sk_11, sk_12, sk_19, sk_21, sk_23, \
                         sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_81 * si0_4[k]
                  + f_81 * si0_10[k]
                  - f_60 * si0_11[k]
                  - f_82 * si1_4[k]
                  - f_82 * si1_10[k]
                  + f_62 * si1_11[k]
                  - 3.28125 * pb_z[k] * sk_0[k]
                  - 9.84375 * pb_z[k] * sk_2[k]
                  + 26.25 * pb_x[k] * sk_5[k]
                  - 9.84375 * pb_z[k] * sk_6[k]
                  + 52.5 * pb_x[k] * sk_11[k]
                  - 31.5 * pb_x[k] * sk_12[k]
                  - 3.28125 * pb_x[k] * sk_19[k]
                  + 26.25 * pb_x[k] * sk_21[k]
                  - 31.5 * pb_x[k] * sk_23[k]
                  + 6.0 * pb_x[k] * sk_25[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_0, si0_1, si0_2, si0_6, si0_7, si0_12, si0_13, \
                         si0_15, si1_0, si1_1, si1_2, si1_6, si1_7, si1_12, si1_13, si1_15, \
                         sk_0, sk_2, sk_3, sk_7, sk_8, sk_13, sk_14, sk_17, sk_18, sk_20, \
                         sk_22, sk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_83 * si0_0[k]
                   - f_84 * si0_1[k]
                   + f_85 * si0_2[k]
                   + f_86 * si0_6[k]
                   - f_49 * si0_7[k]
                   + f_87 * si0_12[k]
                   - f_88 * si0_13[k]
                   + f_49 * si0_15[k]
                   + f_89 * si1_0[k]
                   + f_90 * si1_1[k]
                   - f_91 * si1_2[k]
                   - f_92 * si1_6[k]
                   + f_52 * si1_7[k]
                   - f_93 * si1_12[k]
                   + f_94 * si1_13[k]
                   - f_52 * si1_15[k]
                   - f_95 * pb_x[k] * sk_0[k]
                   - f_53 * pb_x[k] * sk_2[k]
                   + f_96 * pb_x[k] * sk_3[k]
                   + f_96 * pb_x[k] * sk_7[k]
                   - f_97 * pb_x[k] * sk_8[k]
                   + f_53 * pb_x[k] * sk_13[k]
                   - f_96 * pb_x[k] * sk_14[k]
                   + f_98 * pb_x[k] * sk_17[k]
                   + f_95 * pb_y[k] * sk_18[k]
                   - f_96 * pb_y[k] * sk_20[k]
                   + f_97 * pb_y[k] * sk_22[k]
                   - f_98 * pb_y[k] * sk_24[k];
    }

#pragma omp simd aligned(pb_x, pb_z, si0_4, si0_10, si0_11, si1_4, si1_10, si1_11, sk_0, sk_2, \
                         sk_5, sk_6, sk_11, sk_12, sk_19, sk_21, \
                         sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += -f_99 * si0_4[k]
                   + f_99 * si0_10[k]
                   + f_36 * si0_11[k]
                   + f_100 * si1_4[k]
                   - f_100 * si1_10[k]
                   - f_38 * si1_11[k]
                   + f_42 * pb_z[k] * sk_0[k]
                   - f_42 * pb_z[k] * sk_2[k]
                   - f_45 * pb_x[k] * sk_5[k]
                   - f_40 * pb_z[k] * sk_6[k]
                   + f_44 * pb_x[k] * sk_11[k]
                   + f_46 * pb_x[k] * sk_12[k]
                   - f_39 * pb_x[k] * sk_19[k]
                   + f_41 * pb_x[k] * sk_21[k]
                   - f_43 * pb_x[k] * sk_23[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_0, si0_1, si0_2, si0_5, si0_6, si0_7, si0_12, si0_15, \
                         si1_0, si1_1, si1_2, si1_5, si1_6, si1_7, si1_12, si1_15, sk_0, sk_2, \
                         sk_3, sk_6, sk_7, sk_8, sk_13, sk_14, sk_16, sk_18, sk_20, \
                         sk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_12[k] += f_101 * si0_0[k]
                   - f_102 * si0_1[k]
                   - f_103 * si0_2[k]
                   - f_104 * si0_5[k]
                   + f_105 * si0_6[k]
                   + f_103 * si0_7[k]
                   + f_106 * si0_12[k]
                   - f_103 * si0_15[k]
                   - f_107 * si1_0[k]
                   + f_108 * si1_1[k]
                   + f_109 * si1_2[k]
                   + f_110 * si1_5[k]
                   - f_111 * si1_6[k]
                   - f_109 * si1_7[k]
                   - f_112 * si1_12[k]
                   + f_109 * si1_15[k]
                   + f_113 * pb_x[k] * sk_0[k]
                   - f_32 * pb_x[k] * sk_2[k]
                   - f_114 * pb_x[k] * sk_3[k]
                   - f_115 * pb_x[k] * sk_6[k]
                   + f_116 * pb_x[k] * sk_7[k]
                   + f_117 * pb_x[k] * sk_8[k]
                   - f_32 * pb_x[k] * sk_13[k]
                   + f_116 * pb_x[k] * sk_14[k]
                   - f_118 * pb_x[k] * sk_16[k]
                   + f_113 * pb_y[k] * sk_18[k]
                   - f_114 * pb_y[k] * sk_20[k]
                   + f_117 * pb_y[k] * sk_22[k];
    }

#pragma omp simd aligned(pb_x, pb_z, si0_4, si0_10, si1_4, si1_10, sk_0, sk_2, sk_5, sk_6, \
                         sk_11, sk_19, sk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_13[k] += f_119 * si0_4[k]
                   - f_120 * si0_10[k]
                   - f_121 * si1_4[k]
                   + f_122 * si1_10[k]
                   - f_26 * pb_z[k] * sk_0[k]
                   + f_24 * pb_z[k] * sk_2[k]
                   + f_27 * pb_x[k] * sk_5[k]
                   + f_22 * pb_z[k] * sk_6[k]
                   - f_25 * pb_x[k] * sk_11[k]
                   - f_22 * pb_x[k] * sk_19[k]
                   + f_23 * pb_x[k] * sk_21[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_0, si0_1, si0_2, si0_6, si0_12, si0_13, si1_0, si1_1, \
                         si1_2, si1_6, si1_12, si1_13, sk_0, sk_2, sk_3, sk_7, sk_13, sk_14, \
                         sk_18, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_14[k] += -f_123 * si0_0[k]
                   + f_124 * si0_1[k]
                   + f_124 * si0_2[k]
                   - f_125 * si0_6[k]
                   - f_123 * si0_12[k]
                   + f_126 * si0_13[k]
                   + f_127 * si1_0[k]
                   - f_128 * si1_1[k]
                   - f_128 * si1_2[k]
                   + f_129 * si1_6[k]
                   + f_127 * si1_12[k]
                   - f_130 * si1_13[k]
                   - f_131 * pb_x[k] * sk_0[k]
                   + f_18 * pb_x[k] * sk_2[k]
                   + f_18 * pb_x[k] * sk_3[k]
                   - f_132 * pb_x[k] * sk_7[k]
                   - f_18 * pb_x[k] * sk_13[k]
                   + f_132 * pb_x[k] * sk_14[k]
                   + f_131 * pb_y[k] * sk_18[k]
                   - f_18 * pb_y[k] * sk_20[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sk_0, sk_2, sk_6, sk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_15[k] += f_9 * pb_z[k] * sk_0[k]
                   - f_8 * pb_z[k] * sk_2[k]
                   + f_7 * pb_z[k] * sk_6[k]
                   - f_6 * pb_x[k] * sk_19[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si0_0, si0_1, si0_5, si0_12, si1_0, si1_1, si1_5, si1_12, \
                         sk_0, sk_2, sk_6, sk_13, sk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_16[k] += f_133 * si0_0[k]
                   - f_134 * si0_1[k]
                   + f_135 * si0_5[k]
                   - f_136 * si0_12[k]
                   - f_137 * si1_0[k]
                   + f_138 * si1_1[k]
                   - f_139 * si1_5[k]
                   + f_140 * si1_12[k]
                   + f_141 * pb_x[k] * sk_0[k]
                   - f_6 * pb_x[k] * sk_2[k]
                   + f_142 * pb_x[k] * sk_6[k]
                   - f_6 * pb_x[k] * sk_13[k]
                   + f_141 * pb_y[k] * sk_18[k];
    }
}

}  // namespace simdt2ceri
