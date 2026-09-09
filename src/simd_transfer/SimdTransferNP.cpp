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


#include "SimdTransferNP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_np_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ns, const size_t os,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_0 = buffer.data(target + 0 * ncomps + c);
        auto *t_1 = buffer.data(target + 1 * ncomps + c);
        auto *t_2 = buffer.data(target + 2 * ncomps + c);
        auto *t_3 = buffer.data(target + 3 * ncomps + c);
        auto *t_4 = buffer.data(target + 4 * ncomps + c);
        auto *t_5 = buffer.data(target + 5 * ncomps + c);
        auto *t_6 = buffer.data(target + 6 * ncomps + c);
        auto *t_7 = buffer.data(target + 7 * ncomps + c);
        auto *t_8 = buffer.data(target + 8 * ncomps + c);
        auto *t_9 = buffer.data(target + 9 * ncomps + c);
        auto *t_10 = buffer.data(target + 10 * ncomps + c);
        auto *t_11 = buffer.data(target + 11 * ncomps + c);
        auto *t_12 = buffer.data(target + 12 * ncomps + c);
        auto *t_13 = buffer.data(target + 13 * ncomps + c);
        auto *t_14 = buffer.data(target + 14 * ncomps + c);
        auto *t_15 = buffer.data(target + 15 * ncomps + c);
        auto *t_16 = buffer.data(target + 16 * ncomps + c);
        auto *t_17 = buffer.data(target + 17 * ncomps + c);
        auto *t_18 = buffer.data(target + 18 * ncomps + c);
        auto *t_19 = buffer.data(target + 19 * ncomps + c);
        auto *t_20 = buffer.data(target + 20 * ncomps + c);
        auto *t_21 = buffer.data(target + 21 * ncomps + c);
        auto *t_22 = buffer.data(target + 22 * ncomps + c);
        auto *t_23 = buffer.data(target + 23 * ncomps + c);
        auto *t_24 = buffer.data(target + 24 * ncomps + c);
        auto *t_25 = buffer.data(target + 25 * ncomps + c);
        auto *t_26 = buffer.data(target + 26 * ncomps + c);
        auto *t_27 = buffer.data(target + 27 * ncomps + c);
        auto *t_28 = buffer.data(target + 28 * ncomps + c);
        auto *t_29 = buffer.data(target + 29 * ncomps + c);
        auto *t_30 = buffer.data(target + 30 * ncomps + c);
        auto *t_31 = buffer.data(target + 31 * ncomps + c);
        auto *t_32 = buffer.data(target + 32 * ncomps + c);
        auto *t_33 = buffer.data(target + 33 * ncomps + c);
        auto *t_34 = buffer.data(target + 34 * ncomps + c);
        auto *t_35 = buffer.data(target + 35 * ncomps + c);
        auto *t_36 = buffer.data(target + 36 * ncomps + c);
        auto *t_37 = buffer.data(target + 37 * ncomps + c);
        auto *t_38 = buffer.data(target + 38 * ncomps + c);
        auto *t_39 = buffer.data(target + 39 * ncomps + c);
        auto *t_40 = buffer.data(target + 40 * ncomps + c);
        auto *t_41 = buffer.data(target + 41 * ncomps + c);
        auto *t_42 = buffer.data(target + 42 * ncomps + c);
        auto *t_43 = buffer.data(target + 43 * ncomps + c);
        auto *t_44 = buffer.data(target + 44 * ncomps + c);
        auto *t_45 = buffer.data(target + 45 * ncomps + c);
        auto *t_46 = buffer.data(target + 46 * ncomps + c);
        auto *t_47 = buffer.data(target + 47 * ncomps + c);
        auto *t_48 = buffer.data(target + 48 * ncomps + c);
        auto *t_49 = buffer.data(target + 49 * ncomps + c);
        auto *t_50 = buffer.data(target + 50 * ncomps + c);
        auto *t_51 = buffer.data(target + 51 * ncomps + c);
        auto *t_52 = buffer.data(target + 52 * ncomps + c);
        auto *t_53 = buffer.data(target + 53 * ncomps + c);
        auto *t_54 = buffer.data(target + 54 * ncomps + c);
        auto *t_55 = buffer.data(target + 55 * ncomps + c);
        auto *t_56 = buffer.data(target + 56 * ncomps + c);
        auto *t_57 = buffer.data(target + 57 * ncomps + c);
        auto *t_58 = buffer.data(target + 58 * ncomps + c);
        auto *t_59 = buffer.data(target + 59 * ncomps + c);
        auto *t_60 = buffer.data(target + 60 * ncomps + c);
        auto *t_61 = buffer.data(target + 61 * ncomps + c);
        auto *t_62 = buffer.data(target + 62 * ncomps + c);
        auto *t_63 = buffer.data(target + 63 * ncomps + c);
        auto *t_64 = buffer.data(target + 64 * ncomps + c);
        auto *t_65 = buffer.data(target + 65 * ncomps + c);
        auto *t_66 = buffer.data(target + 66 * ncomps + c);
        auto *t_67 = buffer.data(target + 67 * ncomps + c);
        auto *t_68 = buffer.data(target + 68 * ncomps + c);
        auto *t_69 = buffer.data(target + 69 * ncomps + c);
        auto *t_70 = buffer.data(target + 70 * ncomps + c);
        auto *t_71 = buffer.data(target + 71 * ncomps + c);
        auto *t_72 = buffer.data(target + 72 * ncomps + c);
        auto *t_73 = buffer.data(target + 73 * ncomps + c);
        auto *t_74 = buffer.data(target + 74 * ncomps + c);
        auto *t_75 = buffer.data(target + 75 * ncomps + c);
        auto *t_76 = buffer.data(target + 76 * ncomps + c);
        auto *t_77 = buffer.data(target + 77 * ncomps + c);
        auto *t_78 = buffer.data(target + 78 * ncomps + c);
        auto *t_79 = buffer.data(target + 79 * ncomps + c);
        auto *t_80 = buffer.data(target + 80 * ncomps + c);
        auto *t_81 = buffer.data(target + 81 * ncomps + c);
        auto *t_82 = buffer.data(target + 82 * ncomps + c);
        auto *t_83 = buffer.data(target + 83 * ncomps + c);
        auto *t_84 = buffer.data(target + 84 * ncomps + c);
        auto *t_85 = buffer.data(target + 85 * ncomps + c);
        auto *t_86 = buffer.data(target + 86 * ncomps + c);
        auto *t_87 = buffer.data(target + 87 * ncomps + c);
        auto *t_88 = buffer.data(target + 88 * ncomps + c);
        auto *t_89 = buffer.data(target + 89 * ncomps + c);
        auto *t_90 = buffer.data(target + 90 * ncomps + c);
        auto *t_91 = buffer.data(target + 91 * ncomps + c);
        auto *t_92 = buffer.data(target + 92 * ncomps + c);
        auto *t_93 = buffer.data(target + 93 * ncomps + c);
        auto *t_94 = buffer.data(target + 94 * ncomps + c);
        auto *t_95 = buffer.data(target + 95 * ncomps + c);
        auto *t_96 = buffer.data(target + 96 * ncomps + c);
        auto *t_97 = buffer.data(target + 97 * ncomps + c);
        auto *t_98 = buffer.data(target + 98 * ncomps + c);
        auto *t_99 = buffer.data(target + 99 * ncomps + c);
        auto *t_100 = buffer.data(target + 100 * ncomps + c);
        auto *t_101 = buffer.data(target + 101 * ncomps + c);
        auto *t_102 = buffer.data(target + 102 * ncomps + c);
        auto *t_103 = buffer.data(target + 103 * ncomps + c);
        auto *t_104 = buffer.data(target + 104 * ncomps + c);
        auto *t_105 = buffer.data(target + 105 * ncomps + c);
        auto *t_106 = buffer.data(target + 106 * ncomps + c);
        auto *t_107 = buffer.data(target + 107 * ncomps + c);
        auto *t_108 = buffer.data(target + 108 * ncomps + c);
        auto *t_109 = buffer.data(target + 109 * ncomps + c);
        auto *t_110 = buffer.data(target + 110 * ncomps + c);
        auto *t_111 = buffer.data(target + 111 * ncomps + c);
        auto *t_112 = buffer.data(target + 112 * ncomps + c);
        auto *t_113 = buffer.data(target + 113 * ncomps + c);
        auto *t_114 = buffer.data(target + 114 * ncomps + c);
        auto *t_115 = buffer.data(target + 115 * ncomps + c);
        auto *t_116 = buffer.data(target + 116 * ncomps + c);
        auto *t_117 = buffer.data(target + 117 * ncomps + c);
        auto *t_118 = buffer.data(target + 118 * ncomps + c);
        auto *t_119 = buffer.data(target + 119 * ncomps + c);
        auto *t_120 = buffer.data(target + 120 * ncomps + c);
        auto *t_121 = buffer.data(target + 121 * ncomps + c);
        auto *t_122 = buffer.data(target + 122 * ncomps + c);
        auto *t_123 = buffer.data(target + 123 * ncomps + c);
        auto *t_124 = buffer.data(target + 124 * ncomps + c);
        auto *t_125 = buffer.data(target + 125 * ncomps + c);
        auto *t_126 = buffer.data(target + 126 * ncomps + c);
        auto *t_127 = buffer.data(target + 127 * ncomps + c);
        auto *t_128 = buffer.data(target + 128 * ncomps + c);
        auto *t_129 = buffer.data(target + 129 * ncomps + c);
        auto *t_130 = buffer.data(target + 130 * ncomps + c);
        auto *t_131 = buffer.data(target + 131 * ncomps + c);
        auto *t_132 = buffer.data(target + 132 * ncomps + c);
        auto *t_133 = buffer.data(target + 133 * ncomps + c);
        auto *t_134 = buffer.data(target + 134 * ncomps + c);
        auto *t_135 = buffer.data(target + 135 * ncomps + c);
        auto *t_136 = buffer.data(target + 136 * ncomps + c);
        auto *t_137 = buffer.data(target + 137 * ncomps + c);
        auto *t_138 = buffer.data(target + 138 * ncomps + c);
        auto *t_139 = buffer.data(target + 139 * ncomps + c);
        auto *t_140 = buffer.data(target + 140 * ncomps + c);
        auto *t_141 = buffer.data(target + 141 * ncomps + c);
        auto *t_142 = buffer.data(target + 142 * ncomps + c);
        auto *t_143 = buffer.data(target + 143 * ncomps + c);
        auto *t_144 = buffer.data(target + 144 * ncomps + c);
        auto *t_145 = buffer.data(target + 145 * ncomps + c);
        auto *t_146 = buffer.data(target + 146 * ncomps + c);
        auto *t_147 = buffer.data(target + 147 * ncomps + c);
        auto *t_148 = buffer.data(target + 148 * ncomps + c);
        auto *t_149 = buffer.data(target + 149 * ncomps + c);
        auto *t_150 = buffer.data(target + 150 * ncomps + c);
        auto *t_151 = buffer.data(target + 151 * ncomps + c);
        auto *t_152 = buffer.data(target + 152 * ncomps + c);
        auto *t_153 = buffer.data(target + 153 * ncomps + c);
        auto *t_154 = buffer.data(target + 154 * ncomps + c);
        auto *t_155 = buffer.data(target + 155 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ns_0 = buffer.data(ns + 0 * ncomps + c);
        const auto *ns_1 = buffer.data(ns + 1 * ncomps + c);
        const auto *ns_2 = buffer.data(ns + 2 * ncomps + c);
        const auto *ns_3 = buffer.data(ns + 3 * ncomps + c);
        const auto *ns_4 = buffer.data(ns + 4 * ncomps + c);
        const auto *ns_5 = buffer.data(ns + 5 * ncomps + c);
        const auto *ns_6 = buffer.data(ns + 6 * ncomps + c);
        const auto *ns_7 = buffer.data(ns + 7 * ncomps + c);
        const auto *ns_8 = buffer.data(ns + 8 * ncomps + c);
        const auto *ns_9 = buffer.data(ns + 9 * ncomps + c);
        const auto *ns_10 = buffer.data(ns + 10 * ncomps + c);
        const auto *ns_11 = buffer.data(ns + 11 * ncomps + c);
        const auto *ns_12 = buffer.data(ns + 12 * ncomps + c);
        const auto *ns_13 = buffer.data(ns + 13 * ncomps + c);
        const auto *ns_14 = buffer.data(ns + 14 * ncomps + c);
        const auto *ns_15 = buffer.data(ns + 15 * ncomps + c);
        const auto *ns_16 = buffer.data(ns + 16 * ncomps + c);
        const auto *ns_17 = buffer.data(ns + 17 * ncomps + c);
        const auto *ns_18 = buffer.data(ns + 18 * ncomps + c);
        const auto *ns_19 = buffer.data(ns + 19 * ncomps + c);
        const auto *ns_20 = buffer.data(ns + 20 * ncomps + c);
        const auto *ns_21 = buffer.data(ns + 21 * ncomps + c);
        const auto *ns_22 = buffer.data(ns + 22 * ncomps + c);
        const auto *ns_23 = buffer.data(ns + 23 * ncomps + c);
        const auto *ns_24 = buffer.data(ns + 24 * ncomps + c);
        const auto *ns_25 = buffer.data(ns + 25 * ncomps + c);
        const auto *ns_26 = buffer.data(ns + 26 * ncomps + c);
        const auto *ns_27 = buffer.data(ns + 27 * ncomps + c);
        const auto *ns_28 = buffer.data(ns + 28 * ncomps + c);
        const auto *ns_29 = buffer.data(ns + 29 * ncomps + c);
        const auto *ns_30 = buffer.data(ns + 30 * ncomps + c);
        const auto *ns_31 = buffer.data(ns + 31 * ncomps + c);
        const auto *ns_32 = buffer.data(ns + 32 * ncomps + c);
        const auto *ns_33 = buffer.data(ns + 33 * ncomps + c);
        const auto *ns_34 = buffer.data(ns + 34 * ncomps + c);
        const auto *ns_35 = buffer.data(ns + 35 * ncomps + c);
        const auto *ns_36 = buffer.data(ns + 36 * ncomps + c);
        const auto *ns_37 = buffer.data(ns + 37 * ncomps + c);
        const auto *ns_38 = buffer.data(ns + 38 * ncomps + c);
        const auto *ns_39 = buffer.data(ns + 39 * ncomps + c);
        const auto *ns_40 = buffer.data(ns + 40 * ncomps + c);
        const auto *ns_41 = buffer.data(ns + 41 * ncomps + c);
        const auto *ns_42 = buffer.data(ns + 42 * ncomps + c);
        const auto *ns_43 = buffer.data(ns + 43 * ncomps + c);
        const auto *ns_44 = buffer.data(ns + 44 * ncomps + c);
        const auto *ns_45 = buffer.data(ns + 45 * ncomps + c);
        const auto *ns_46 = buffer.data(ns + 46 * ncomps + c);
        const auto *ns_47 = buffer.data(ns + 47 * ncomps + c);
        const auto *ns_48 = buffer.data(ns + 48 * ncomps + c);
        const auto *ns_49 = buffer.data(ns + 49 * ncomps + c);
        const auto *ns_50 = buffer.data(ns + 50 * ncomps + c);
        const auto *ns_51 = buffer.data(ns + 51 * ncomps + c);

        const auto *os_0 = buffer.data(os + 0 * ncomps + c);
        const auto *os_1 = buffer.data(os + 1 * ncomps + c);
        const auto *os_2 = buffer.data(os + 2 * ncomps + c);
        const auto *os_3 = buffer.data(os + 3 * ncomps + c);
        const auto *os_4 = buffer.data(os + 4 * ncomps + c);
        const auto *os_5 = buffer.data(os + 5 * ncomps + c);
        const auto *os_6 = buffer.data(os + 6 * ncomps + c);
        const auto *os_7 = buffer.data(os + 7 * ncomps + c);
        const auto *os_8 = buffer.data(os + 8 * ncomps + c);
        const auto *os_9 = buffer.data(os + 9 * ncomps + c);
        const auto *os_10 = buffer.data(os + 10 * ncomps + c);
        const auto *os_11 = buffer.data(os + 11 * ncomps + c);
        const auto *os_12 = buffer.data(os + 12 * ncomps + c);
        const auto *os_13 = buffer.data(os + 13 * ncomps + c);
        const auto *os_14 = buffer.data(os + 14 * ncomps + c);
        const auto *os_15 = buffer.data(os + 15 * ncomps + c);
        const auto *os_16 = buffer.data(os + 16 * ncomps + c);
        const auto *os_17 = buffer.data(os + 17 * ncomps + c);
        const auto *os_18 = buffer.data(os + 18 * ncomps + c);
        const auto *os_19 = buffer.data(os + 19 * ncomps + c);
        const auto *os_20 = buffer.data(os + 20 * ncomps + c);
        const auto *os_21 = buffer.data(os + 21 * ncomps + c);
        const auto *os_22 = buffer.data(os + 22 * ncomps + c);
        const auto *os_23 = buffer.data(os + 23 * ncomps + c);
        const auto *os_24 = buffer.data(os + 24 * ncomps + c);
        const auto *os_25 = buffer.data(os + 25 * ncomps + c);
        const auto *os_26 = buffer.data(os + 26 * ncomps + c);
        const auto *os_27 = buffer.data(os + 27 * ncomps + c);
        const auto *os_28 = buffer.data(os + 28 * ncomps + c);
        const auto *os_29 = buffer.data(os + 29 * ncomps + c);
        const auto *os_30 = buffer.data(os + 30 * ncomps + c);
        const auto *os_31 = buffer.data(os + 31 * ncomps + c);
        const auto *os_32 = buffer.data(os + 32 * ncomps + c);
        const auto *os_33 = buffer.data(os + 33 * ncomps + c);
        const auto *os_34 = buffer.data(os + 34 * ncomps + c);
        const auto *os_35 = buffer.data(os + 35 * ncomps + c);
        const auto *os_36 = buffer.data(os + 36 * ncomps + c);
        const auto *os_37 = buffer.data(os + 37 * ncomps + c);
        const auto *os_38 = buffer.data(os + 38 * ncomps + c);
        const auto *os_39 = buffer.data(os + 39 * ncomps + c);
        const auto *os_40 = buffer.data(os + 40 * ncomps + c);
        const auto *os_41 = buffer.data(os + 41 * ncomps + c);
        const auto *os_42 = buffer.data(os + 42 * ncomps + c);
        const auto *os_43 = buffer.data(os + 43 * ncomps + c);
        const auto *os_44 = buffer.data(os + 44 * ncomps + c);
        const auto *os_45 = buffer.data(os + 45 * ncomps + c);
        const auto *os_46 = buffer.data(os + 46 * ncomps + c);
        const auto *os_47 = buffer.data(os + 47 * ncomps + c);
        const auto *os_48 = buffer.data(os + 48 * ncomps + c);
        const auto *os_49 = buffer.data(os + 49 * ncomps + c);
        const auto *os_50 = buffer.data(os + 50 * ncomps + c);
        const auto *os_51 = buffer.data(os + 51 * ncomps + c);
        const auto *os_52 = buffer.data(os + 52 * ncomps + c);
        const auto *os_53 = buffer.data(os + 53 * ncomps + c);
        const auto *os_54 = buffer.data(os + 54 * ncomps + c);
        const auto *os_55 = buffer.data(os + 55 * ncomps + c);
        const auto *os_56 = buffer.data(os + 56 * ncomps + c);
        const auto *os_57 = buffer.data(os + 57 * ncomps + c);
        const auto *os_58 = buffer.data(os + 58 * ncomps + c);
        const auto *os_59 = buffer.data(os + 59 * ncomps + c);
        const auto *os_60 = buffer.data(os + 60 * ncomps + c);
        const auto *os_61 = buffer.data(os + 61 * ncomps + c);
        const auto *os_62 = buffer.data(os + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, ns_0, ns_1, os_0, \
                         os_1, os_2, os_3, os_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ns_0[k]
                     + os_0[k];

            t_1[k] = ab_y[k] * ns_0[k]
                     + os_1[k];

            t_2[k] = ab_z[k] * ns_0[k]
                     + os_2[k];

            t_3[k] = ab_x[k] * ns_1[k]
                     + os_1[k];

            t_4[k] = ab_y[k] * ns_1[k]
                     + os_3[k];

            t_5[k] = ab_z[k] * ns_1[k]
                     + os_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, ns_2, ns_3, os_2, os_3, \
                         os_4, os_5, os_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * ns_2[k]
                     + os_2[k];

            t_7[k] = ab_y[k] * ns_2[k]
                     + os_4[k];

            t_8[k] = ab_z[k] * ns_2[k]
                     + os_5[k];

            t_9[k] = ab_x[k] * ns_3[k]
                     + os_3[k];

            t_10[k] = ab_y[k] * ns_3[k]
                      + os_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, ns_3, ns_4, \
                         ns_5, os_4, os_5, os_7, os_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * ns_3[k]
                      + os_7[k];

            t_12[k] = ab_x[k] * ns_4[k]
                      + os_4[k];

            t_13[k] = ab_y[k] * ns_4[k]
                      + os_7[k];

            t_14[k] = ab_z[k] * ns_4[k]
                      + os_8[k];

            t_15[k] = ab_x[k] * ns_5[k]
                      + os_5[k];

            t_16[k] = ab_y[k] * ns_5[k]
                      + os_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, ns_5, ns_6, ns_7, \
                         os_6, os_7, os_9, os_10, os_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * ns_5[k]
                      + os_9[k];

            t_18[k] = ab_x[k] * ns_6[k]
                      + os_6[k];

            t_19[k] = ab_y[k] * ns_6[k]
                      + os_10[k];

            t_20[k] = ab_z[k] * ns_6[k]
                      + os_11[k];

            t_21[k] = ab_x[k] * ns_7[k]
                      + os_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, ns_7, ns_8, os_8, \
                         os_11, os_12, os_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * ns_7[k]
                      + os_11[k];

            t_23[k] = ab_z[k] * ns_7[k]
                      + os_12[k];

            t_24[k] = ab_x[k] * ns_8[k]
                      + os_8[k];

            t_25[k] = ab_y[k] * ns_8[k]
                      + os_12[k];

            t_26[k] = ab_z[k] * ns_8[k]
                      + os_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ns_9, ns_10, os_9, \
                         os_10, os_13, os_14, os_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * ns_9[k]
                      + os_9[k];

            t_28[k] = ab_y[k] * ns_9[k]
                      + os_13[k];

            t_29[k] = ab_z[k] * ns_9[k]
                      + os_14[k];

            t_30[k] = ab_x[k] * ns_10[k]
                      + os_10[k];

            t_31[k] = ab_y[k] * ns_10[k]
                      + os_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, ns_10, ns_11, \
                         ns_12, os_11, os_12, os_16, os_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * ns_10[k]
                      + os_16[k];

            t_33[k] = ab_x[k] * ns_11[k]
                      + os_11[k];

            t_34[k] = ab_y[k] * ns_11[k]
                      + os_16[k];

            t_35[k] = ab_z[k] * ns_11[k]
                      + os_17[k];

            t_36[k] = ab_x[k] * ns_12[k]
                      + os_12[k];

            t_37[k] = ab_y[k] * ns_12[k]
                      + os_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, ns_12, ns_13, \
                         ns_14, os_13, os_14, os_18, os_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * ns_12[k]
                      + os_18[k];

            t_39[k] = ab_x[k] * ns_13[k]
                      + os_13[k];

            t_40[k] = ab_y[k] * ns_13[k]
                      + os_18[k];

            t_41[k] = ab_z[k] * ns_13[k]
                      + os_19[k];

            t_42[k] = ab_x[k] * ns_14[k]
                      + os_14[k];

            t_43[k] = ab_y[k] * ns_14[k]
                      + os_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, ns_14, ns_15, ns_16, \
                         os_15, os_16, os_20, os_21, os_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * ns_14[k]
                      + os_20[k];

            t_45[k] = ab_x[k] * ns_15[k]
                      + os_15[k];

            t_46[k] = ab_y[k] * ns_15[k]
                      + os_21[k];

            t_47[k] = ab_z[k] * ns_15[k]
                      + os_22[k];

            t_48[k] = ab_x[k] * ns_16[k]
                      + os_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ns_16, ns_17, os_17, \
                         os_22, os_23, os_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * ns_16[k]
                      + os_22[k];

            t_50[k] = ab_z[k] * ns_16[k]
                      + os_23[k];

            t_51[k] = ab_x[k] * ns_17[k]
                      + os_17[k];

            t_52[k] = ab_y[k] * ns_17[k]
                      + os_23[k];

            t_53[k] = ab_z[k] * ns_17[k]
                      + os_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ns_18, ns_19, \
                         os_18, os_19, os_24, os_25, os_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * ns_18[k]
                      + os_18[k];

            t_55[k] = ab_y[k] * ns_18[k]
                      + os_24[k];

            t_56[k] = ab_z[k] * ns_18[k]
                      + os_25[k];

            t_57[k] = ab_x[k] * ns_19[k]
                      + os_19[k];

            t_58[k] = ab_y[k] * ns_19[k]
                      + os_25[k];

            t_59[k] = ab_z[k] * ns_19[k]
                      + os_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, ns_20, ns_21, os_20, \
                         os_21, os_26, os_27, os_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ns_20[k]
                      + os_20[k];

            t_61[k] = ab_y[k] * ns_20[k]
                      + os_26[k];

            t_62[k] = ab_z[k] * ns_20[k]
                      + os_27[k];

            t_63[k] = ab_x[k] * ns_21[k]
                      + os_21[k];

            t_64[k] = ab_y[k] * ns_21[k]
                      + os_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, ns_21, ns_22, \
                         ns_23, os_22, os_23, os_29, os_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * ns_21[k]
                      + os_29[k];

            t_66[k] = ab_x[k] * ns_22[k]
                      + os_22[k];

            t_67[k] = ab_y[k] * ns_22[k]
                      + os_29[k];

            t_68[k] = ab_z[k] * ns_22[k]
                      + os_30[k];

            t_69[k] = ab_x[k] * ns_23[k]
                      + os_23[k];

            t_70[k] = ab_y[k] * ns_23[k]
                      + os_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, ns_23, ns_24, \
                         ns_25, os_24, os_25, os_31, os_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * ns_23[k]
                      + os_31[k];

            t_72[k] = ab_x[k] * ns_24[k]
                      + os_24[k];

            t_73[k] = ab_y[k] * ns_24[k]
                      + os_31[k];

            t_74[k] = ab_z[k] * ns_24[k]
                      + os_32[k];

            t_75[k] = ab_x[k] * ns_25[k]
                      + os_25[k];

            t_76[k] = ab_y[k] * ns_25[k]
                      + os_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, ns_25, ns_26, \
                         ns_27, os_26, os_27, os_33, os_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * ns_25[k]
                      + os_33[k];

            t_78[k] = ab_x[k] * ns_26[k]
                      + os_26[k];

            t_79[k] = ab_y[k] * ns_26[k]
                      + os_33[k];

            t_80[k] = ab_z[k] * ns_26[k]
                      + os_34[k];

            t_81[k] = ab_x[k] * ns_27[k]
                      + os_27[k];

            t_82[k] = ab_y[k] * ns_27[k]
                      + os_34[k];
        }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ab_z, ns_27, ns_28, ns_29, \
                         os_28, os_29, os_35, os_36, os_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * ns_27[k]
                      + os_35[k];

            t_84[k] = ab_x[k] * ns_28[k]
                      + os_28[k];

            t_85[k] = ab_y[k] * ns_28[k]
                      + os_36[k];

            t_86[k] = ab_z[k] * ns_28[k]
                      + os_37[k];

            t_87[k] = ab_x[k] * ns_29[k]
                      + os_29[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ab_x, ab_y, ab_z, ns_29, ns_30, os_30, \
                         os_37, os_38, os_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * ns_29[k]
                      + os_37[k];

            t_89[k] = ab_z[k] * ns_29[k]
                      + os_38[k];

            t_90[k] = ab_x[k] * ns_30[k]
                      + os_30[k];

            t_91[k] = ab_y[k] * ns_30[k]
                      + os_38[k];

            t_92[k] = ab_z[k] * ns_30[k]
                      + os_39[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, ns_31, ns_32, \
                         os_31, os_32, os_39, os_40, os_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = ab_x[k] * ns_31[k]
                      + os_31[k];

            t_94[k] = ab_y[k] * ns_31[k]
                      + os_39[k];

            t_95[k] = ab_z[k] * ns_31[k]
                      + os_40[k];

            t_96[k] = ab_x[k] * ns_32[k]
                      + os_32[k];

            t_97[k] = ab_y[k] * ns_32[k]
                      + os_40[k];

            t_98[k] = ab_z[k] * ns_32[k]
                      + os_41[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, ab_x, ab_y, ab_z, ns_33, \
                         ns_34, os_33, os_34, os_41, os_42, os_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * ns_33[k]
                      + os_33[k];

            t_100[k] = ab_y[k] * ns_33[k]
                       + os_41[k];

            t_101[k] = ab_z[k] * ns_33[k]
                       + os_42[k];

            t_102[k] = ab_x[k] * ns_34[k]
                       + os_34[k];

            t_103[k] = ab_y[k] * ns_34[k]
                       + os_42[k];

            t_104[k] = ab_z[k] * ns_34[k]
                       + os_43[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, ns_35, ns_36, \
                         os_35, os_36, os_43, os_44, os_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * ns_35[k]
                       + os_35[k];

            t_106[k] = ab_y[k] * ns_35[k]
                       + os_43[k];

            t_107[k] = ab_z[k] * ns_35[k]
                       + os_44[k];

            t_108[k] = ab_x[k] * ns_36[k]
                       + os_36[k];

            t_109[k] = ab_y[k] * ns_36[k]
                       + os_45[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, ab_x, ab_y, ab_z, ns_36, \
                         ns_37, ns_38, os_37, os_38, os_46, os_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_z[k] * ns_36[k]
                       + os_46[k];

            t_111[k] = ab_x[k] * ns_37[k]
                       + os_37[k];

            t_112[k] = ab_y[k] * ns_37[k]
                       + os_46[k];

            t_113[k] = ab_z[k] * ns_37[k]
                       + os_47[k];

            t_114[k] = ab_x[k] * ns_38[k]
                       + os_38[k];

            t_115[k] = ab_y[k] * ns_38[k]
                       + os_47[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, ab_x, ab_y, ab_z, ns_38, \
                         ns_39, ns_40, os_39, os_40, os_48, os_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_z[k] * ns_38[k]
                       + os_48[k];

            t_117[k] = ab_x[k] * ns_39[k]
                       + os_39[k];

            t_118[k] = ab_y[k] * ns_39[k]
                       + os_48[k];

            t_119[k] = ab_z[k] * ns_39[k]
                       + os_49[k];

            t_120[k] = ab_x[k] * ns_40[k]
                       + os_40[k];

            t_121[k] = ab_y[k] * ns_40[k]
                       + os_49[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, ab_x, ab_y, ab_z, ns_40, \
                         ns_41, ns_42, os_41, os_42, os_50, os_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_z[k] * ns_40[k]
                       + os_50[k];

            t_123[k] = ab_x[k] * ns_41[k]
                       + os_41[k];

            t_124[k] = ab_y[k] * ns_41[k]
                       + os_50[k];

            t_125[k] = ab_z[k] * ns_41[k]
                       + os_51[k];

            t_126[k] = ab_x[k] * ns_42[k]
                       + os_42[k];

            t_127[k] = ab_y[k] * ns_42[k]
                       + os_51[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, ab_x, ab_y, ab_z, ns_42, \
                         ns_43, ns_44, os_43, os_44, os_52, os_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = ab_z[k] * ns_42[k]
                       + os_52[k];

            t_129[k] = ab_x[k] * ns_43[k]
                       + os_43[k];

            t_130[k] = ab_y[k] * ns_43[k]
                       + os_52[k];

            t_131[k] = ab_z[k] * ns_43[k]
                       + os_53[k];

            t_132[k] = ab_x[k] * ns_44[k]
                       + os_44[k];

            t_133[k] = ab_y[k] * ns_44[k]
                       + os_53[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, ns_44, ns_45, \
                         ns_46, os_45, os_46, os_54, os_55, os_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = ab_z[k] * ns_44[k]
                       + os_54[k];

            t_135[k] = ab_x[k] * ns_45[k]
                       + os_45[k];

            t_136[k] = ab_y[k] * ns_45[k]
                       + os_55[k];

            t_137[k] = ab_z[k] * ns_45[k]
                       + os_56[k];

            t_138[k] = ab_x[k] * ns_46[k]
                       + os_46[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, ns_46, ns_47, \
                         os_47, os_56, os_57, os_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_y[k] * ns_46[k]
                       + os_56[k];

            t_140[k] = ab_z[k] * ns_46[k]
                       + os_57[k];

            t_141[k] = ab_x[k] * ns_47[k]
                       + os_47[k];

            t_142[k] = ab_y[k] * ns_47[k]
                       + os_57[k];

            t_143[k] = ab_z[k] * ns_47[k]
                       + os_58[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, ns_48, \
                         ns_49, os_48, os_49, os_58, os_59, os_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * ns_48[k]
                       + os_48[k];

            t_145[k] = ab_y[k] * ns_48[k]
                       + os_58[k];

            t_146[k] = ab_z[k] * ns_48[k]
                       + os_59[k];

            t_147[k] = ab_x[k] * ns_49[k]
                       + os_49[k];

            t_148[k] = ab_y[k] * ns_49[k]
                       + os_59[k];

            t_149[k] = ab_z[k] * ns_49[k]
                       + os_60[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, ab_x, ab_y, ab_z, ns_50, \
                         ns_51, os_50, os_51, os_60, os_61, os_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ns_50[k]
                       + os_50[k];

            t_151[k] = ab_y[k] * ns_50[k]
                       + os_60[k];

            t_152[k] = ab_z[k] * ns_50[k]
                       + os_61[k];

            t_153[k] = ab_x[k] * ns_51[k]
                       + os_51[k];

            t_154[k] = ab_y[k] * ns_51[k]
                       + os_61[k];

            t_155[k] = ab_z[k] * ns_51[k]
                       + os_62[k];
        }
    }
}

static auto
compute_hrr_np_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ns, const size_t os,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_156 = buffer.data(target + 156 * ncomps + c);
        auto *t_157 = buffer.data(target + 157 * ncomps + c);
        auto *t_158 = buffer.data(target + 158 * ncomps + c);
        auto *t_159 = buffer.data(target + 159 * ncomps + c);
        auto *t_160 = buffer.data(target + 160 * ncomps + c);
        auto *t_161 = buffer.data(target + 161 * ncomps + c);
        auto *t_162 = buffer.data(target + 162 * ncomps + c);
        auto *t_163 = buffer.data(target + 163 * ncomps + c);
        auto *t_164 = buffer.data(target + 164 * ncomps + c);
        auto *t_165 = buffer.data(target + 165 * ncomps + c);
        auto *t_166 = buffer.data(target + 166 * ncomps + c);
        auto *t_167 = buffer.data(target + 167 * ncomps + c);
        auto *t_168 = buffer.data(target + 168 * ncomps + c);
        auto *t_169 = buffer.data(target + 169 * ncomps + c);
        auto *t_170 = buffer.data(target + 170 * ncomps + c);
        auto *t_171 = buffer.data(target + 171 * ncomps + c);
        auto *t_172 = buffer.data(target + 172 * ncomps + c);
        auto *t_173 = buffer.data(target + 173 * ncomps + c);
        auto *t_174 = buffer.data(target + 174 * ncomps + c);
        auto *t_175 = buffer.data(target + 175 * ncomps + c);
        auto *t_176 = buffer.data(target + 176 * ncomps + c);
        auto *t_177 = buffer.data(target + 177 * ncomps + c);
        auto *t_178 = buffer.data(target + 178 * ncomps + c);
        auto *t_179 = buffer.data(target + 179 * ncomps + c);
        auto *t_180 = buffer.data(target + 180 * ncomps + c);
        auto *t_181 = buffer.data(target + 181 * ncomps + c);
        auto *t_182 = buffer.data(target + 182 * ncomps + c);
        auto *t_183 = buffer.data(target + 183 * ncomps + c);
        auto *t_184 = buffer.data(target + 184 * ncomps + c);
        auto *t_185 = buffer.data(target + 185 * ncomps + c);
        auto *t_186 = buffer.data(target + 186 * ncomps + c);
        auto *t_187 = buffer.data(target + 187 * ncomps + c);
        auto *t_188 = buffer.data(target + 188 * ncomps + c);
        auto *t_189 = buffer.data(target + 189 * ncomps + c);
        auto *t_190 = buffer.data(target + 190 * ncomps + c);
        auto *t_191 = buffer.data(target + 191 * ncomps + c);
        auto *t_192 = buffer.data(target + 192 * ncomps + c);
        auto *t_193 = buffer.data(target + 193 * ncomps + c);
        auto *t_194 = buffer.data(target + 194 * ncomps + c);
        auto *t_195 = buffer.data(target + 195 * ncomps + c);
        auto *t_196 = buffer.data(target + 196 * ncomps + c);
        auto *t_197 = buffer.data(target + 197 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ns_52 = buffer.data(ns + 52 * ncomps + c);
        const auto *ns_53 = buffer.data(ns + 53 * ncomps + c);
        const auto *ns_54 = buffer.data(ns + 54 * ncomps + c);
        const auto *ns_55 = buffer.data(ns + 55 * ncomps + c);
        const auto *ns_56 = buffer.data(ns + 56 * ncomps + c);
        const auto *ns_57 = buffer.data(ns + 57 * ncomps + c);
        const auto *ns_58 = buffer.data(ns + 58 * ncomps + c);
        const auto *ns_59 = buffer.data(ns + 59 * ncomps + c);
        const auto *ns_60 = buffer.data(ns + 60 * ncomps + c);
        const auto *ns_61 = buffer.data(ns + 61 * ncomps + c);
        const auto *ns_62 = buffer.data(ns + 62 * ncomps + c);
        const auto *ns_63 = buffer.data(ns + 63 * ncomps + c);
        const auto *ns_64 = buffer.data(ns + 64 * ncomps + c);
        const auto *ns_65 = buffer.data(ns + 65 * ncomps + c);

        const auto *os_52 = buffer.data(os + 52 * ncomps + c);
        const auto *os_53 = buffer.data(os + 53 * ncomps + c);
        const auto *os_54 = buffer.data(os + 54 * ncomps + c);
        const auto *os_55 = buffer.data(os + 55 * ncomps + c);
        const auto *os_56 = buffer.data(os + 56 * ncomps + c);
        const auto *os_57 = buffer.data(os + 57 * ncomps + c);
        const auto *os_58 = buffer.data(os + 58 * ncomps + c);
        const auto *os_59 = buffer.data(os + 59 * ncomps + c);
        const auto *os_60 = buffer.data(os + 60 * ncomps + c);
        const auto *os_61 = buffer.data(os + 61 * ncomps + c);
        const auto *os_62 = buffer.data(os + 62 * ncomps + c);
        const auto *os_63 = buffer.data(os + 63 * ncomps + c);
        const auto *os_64 = buffer.data(os + 64 * ncomps + c);
        const auto *os_65 = buffer.data(os + 65 * ncomps + c);
        const auto *os_66 = buffer.data(os + 66 * ncomps + c);
        const auto *os_67 = buffer.data(os + 67 * ncomps + c);
        const auto *os_68 = buffer.data(os + 68 * ncomps + c);
        const auto *os_69 = buffer.data(os + 69 * ncomps + c);
        const auto *os_70 = buffer.data(os + 70 * ncomps + c);
        const auto *os_71 = buffer.data(os + 71 * ncomps + c);
        const auto *os_72 = buffer.data(os + 72 * ncomps + c);
        const auto *os_73 = buffer.data(os + 73 * ncomps + c);
        const auto *os_74 = buffer.data(os + 74 * ncomps + c);
        const auto *os_75 = buffer.data(os + 75 * ncomps + c);
        const auto *os_76 = buffer.data(os + 76 * ncomps + c);
        const auto *os_77 = buffer.data(os + 77 * ncomps + c);

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, ns_52, \
                         ns_53, os_52, os_53, os_62, os_63, os_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = ab_x[k] * ns_52[k]
                       + os_52[k];

            t_157[k] = ab_y[k] * ns_52[k]
                       + os_62[k];

            t_158[k] = ab_z[k] * ns_52[k]
                       + os_63[k];

            t_159[k] = ab_x[k] * ns_53[k]
                       + os_53[k];

            t_160[k] = ab_y[k] * ns_53[k]
                       + os_63[k];

            t_161[k] = ab_z[k] * ns_53[k]
                       + os_64[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, ab_z, ns_54, ns_55, \
                         os_54, os_55, os_64, os_65, os_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * ns_54[k]
                       + os_54[k];

            t_163[k] = ab_y[k] * ns_54[k]
                       + os_64[k];

            t_164[k] = ab_z[k] * ns_54[k]
                       + os_65[k];

            t_165[k] = ab_x[k] * ns_55[k]
                       + os_55[k];

            t_166[k] = ab_y[k] * ns_55[k]
                       + os_66[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, t_172, ab_x, ab_y, ab_z, ns_55, \
                         ns_56, ns_57, os_56, os_57, os_67, os_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = ab_z[k] * ns_55[k]
                       + os_67[k];

            t_168[k] = ab_x[k] * ns_56[k]
                       + os_56[k];

            t_169[k] = ab_y[k] * ns_56[k]
                       + os_67[k];

            t_170[k] = ab_z[k] * ns_56[k]
                       + os_68[k];

            t_171[k] = ab_x[k] * ns_57[k]
                       + os_57[k];

            t_172[k] = ab_y[k] * ns_57[k]
                       + os_68[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, ab_x, ab_y, ab_z, ns_57, \
                         ns_58, ns_59, os_58, os_59, os_69, os_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_z[k] * ns_57[k]
                       + os_69[k];

            t_174[k] = ab_x[k] * ns_58[k]
                       + os_58[k];

            t_175[k] = ab_y[k] * ns_58[k]
                       + os_69[k];

            t_176[k] = ab_z[k] * ns_58[k]
                       + os_70[k];

            t_177[k] = ab_x[k] * ns_59[k]
                       + os_59[k];

            t_178[k] = ab_y[k] * ns_59[k]
                       + os_70[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, ab_z, ns_59, \
                         ns_60, ns_61, os_60, os_61, os_71, os_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = ab_z[k] * ns_59[k]
                       + os_71[k];

            t_180[k] = ab_x[k] * ns_60[k]
                       + os_60[k];

            t_181[k] = ab_y[k] * ns_60[k]
                       + os_71[k];

            t_182[k] = ab_z[k] * ns_60[k]
                       + os_72[k];

            t_183[k] = ab_x[k] * ns_61[k]
                       + os_61[k];

            t_184[k] = ab_y[k] * ns_61[k]
                       + os_72[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, ab_x, ab_y, ab_z, ns_61, \
                         ns_62, ns_63, os_62, os_63, os_73, os_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_z[k] * ns_61[k]
                       + os_73[k];

            t_186[k] = ab_x[k] * ns_62[k]
                       + os_62[k];

            t_187[k] = ab_y[k] * ns_62[k]
                       + os_73[k];

            t_188[k] = ab_z[k] * ns_62[k]
                       + os_74[k];

            t_189[k] = ab_x[k] * ns_63[k]
                       + os_63[k];

            t_190[k] = ab_y[k] * ns_63[k]
                       + os_74[k];
        }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, ab_x, ab_y, ab_z, ns_63, \
                         ns_64, ns_65, os_64, os_65, os_75, os_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_191[k] = ab_z[k] * ns_63[k]
                       + os_75[k];

            t_192[k] = ab_x[k] * ns_64[k]
                       + os_64[k];

            t_193[k] = ab_y[k] * ns_64[k]
                       + os_75[k];

            t_194[k] = ab_z[k] * ns_64[k]
                       + os_76[k];

            t_195[k] = ab_x[k] * ns_65[k]
                       + os_65[k];

            t_196[k] = ab_y[k] * ns_65[k]
                       + os_76[k];
        }

#pragma omp simd aligned(t_197, ab_z, ns_65, os_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_z[k] * ns_65[k]
                       + os_77[k];
        }
    }
}

auto
compute_hrr_np_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t ns, const size_t os,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_np_out_of_first_piece0(buffer, coordinates, target, ns, os, ncomps, nmax);

    compute_hrr_np_out_of_first_piece1(buffer, coordinates, target, ns, os, ncomps, nmax);
}

static auto
compute_hrr_np_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ns, const size_t os, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_0 = buffer.data(target + 0 * ncomps + c);
        auto *t_1 = buffer.data(target + 1 * ncomps + c);
        auto *t_2 = buffer.data(target + 2 * ncomps + c);
        auto *t_3 = buffer.data(target + 3 * ncomps + c);
        auto *t_4 = buffer.data(target + 4 * ncomps + c);
        auto *t_5 = buffer.data(target + 5 * ncomps + c);
        auto *t_6 = buffer.data(target + 6 * ncomps + c);
        auto *t_7 = buffer.data(target + 7 * ncomps + c);
        auto *t_8 = buffer.data(target + 8 * ncomps + c);
        auto *t_9 = buffer.data(target + 9 * ncomps + c);
        auto *t_10 = buffer.data(target + 10 * ncomps + c);
        auto *t_11 = buffer.data(target + 11 * ncomps + c);
        auto *t_12 = buffer.data(target + 12 * ncomps + c);
        auto *t_13 = buffer.data(target + 13 * ncomps + c);
        auto *t_14 = buffer.data(target + 14 * ncomps + c);
        auto *t_15 = buffer.data(target + 15 * ncomps + c);
        auto *t_16 = buffer.data(target + 16 * ncomps + c);
        auto *t_17 = buffer.data(target + 17 * ncomps + c);
        auto *t_18 = buffer.data(target + 18 * ncomps + c);
        auto *t_19 = buffer.data(target + 19 * ncomps + c);
        auto *t_20 = buffer.data(target + 20 * ncomps + c);
        auto *t_21 = buffer.data(target + 21 * ncomps + c);
        auto *t_22 = buffer.data(target + 22 * ncomps + c);
        auto *t_23 = buffer.data(target + 23 * ncomps + c);
        auto *t_24 = buffer.data(target + 24 * ncomps + c);
        auto *t_25 = buffer.data(target + 25 * ncomps + c);
        auto *t_26 = buffer.data(target + 26 * ncomps + c);
        auto *t_27 = buffer.data(target + 27 * ncomps + c);
        auto *t_28 = buffer.data(target + 28 * ncomps + c);
        auto *t_29 = buffer.data(target + 29 * ncomps + c);
        auto *t_30 = buffer.data(target + 30 * ncomps + c);
        auto *t_31 = buffer.data(target + 31 * ncomps + c);
        auto *t_32 = buffer.data(target + 32 * ncomps + c);
        auto *t_33 = buffer.data(target + 33 * ncomps + c);
        auto *t_34 = buffer.data(target + 34 * ncomps + c);
        auto *t_35 = buffer.data(target + 35 * ncomps + c);
        auto *t_36 = buffer.data(target + 36 * ncomps + c);
        auto *t_37 = buffer.data(target + 37 * ncomps + c);
        auto *t_38 = buffer.data(target + 38 * ncomps + c);
        auto *t_39 = buffer.data(target + 39 * ncomps + c);
        auto *t_40 = buffer.data(target + 40 * ncomps + c);
        auto *t_41 = buffer.data(target + 41 * ncomps + c);
        auto *t_42 = buffer.data(target + 42 * ncomps + c);
        auto *t_43 = buffer.data(target + 43 * ncomps + c);
        auto *t_44 = buffer.data(target + 44 * ncomps + c);
        auto *t_45 = buffer.data(target + 45 * ncomps + c);
        auto *t_46 = buffer.data(target + 46 * ncomps + c);
        auto *t_47 = buffer.data(target + 47 * ncomps + c);
        auto *t_48 = buffer.data(target + 48 * ncomps + c);
        auto *t_49 = buffer.data(target + 49 * ncomps + c);
        auto *t_50 = buffer.data(target + 50 * ncomps + c);
        auto *t_51 = buffer.data(target + 51 * ncomps + c);
        auto *t_52 = buffer.data(target + 52 * ncomps + c);
        auto *t_53 = buffer.data(target + 53 * ncomps + c);
        auto *t_54 = buffer.data(target + 54 * ncomps + c);
        auto *t_55 = buffer.data(target + 55 * ncomps + c);
        auto *t_56 = buffer.data(target + 56 * ncomps + c);
        auto *t_57 = buffer.data(target + 57 * ncomps + c);
        auto *t_58 = buffer.data(target + 58 * ncomps + c);
        auto *t_59 = buffer.data(target + 59 * ncomps + c);
        auto *t_60 = buffer.data(target + 60 * ncomps + c);
        auto *t_61 = buffer.data(target + 61 * ncomps + c);
        auto *t_62 = buffer.data(target + 62 * ncomps + c);
        auto *t_63 = buffer.data(target + 63 * ncomps + c);
        auto *t_64 = buffer.data(target + 64 * ncomps + c);
        auto *t_65 = buffer.data(target + 65 * ncomps + c);
        auto *t_66 = buffer.data(target + 66 * ncomps + c);
        auto *t_67 = buffer.data(target + 67 * ncomps + c);
        auto *t_68 = buffer.data(target + 68 * ncomps + c);
        auto *t_69 = buffer.data(target + 69 * ncomps + c);
        auto *t_70 = buffer.data(target + 70 * ncomps + c);
        auto *t_71 = buffer.data(target + 71 * ncomps + c);
        auto *t_72 = buffer.data(target + 72 * ncomps + c);
        auto *t_73 = buffer.data(target + 73 * ncomps + c);
        auto *t_74 = buffer.data(target + 74 * ncomps + c);
        auto *t_75 = buffer.data(target + 75 * ncomps + c);
        auto *t_76 = buffer.data(target + 76 * ncomps + c);
        auto *t_77 = buffer.data(target + 77 * ncomps + c);
        auto *t_78 = buffer.data(target + 78 * ncomps + c);
        auto *t_79 = buffer.data(target + 79 * ncomps + c);
        auto *t_80 = buffer.data(target + 80 * ncomps + c);
        auto *t_81 = buffer.data(target + 81 * ncomps + c);
        auto *t_82 = buffer.data(target + 82 * ncomps + c);
        auto *t_83 = buffer.data(target + 83 * ncomps + c);
        auto *t_84 = buffer.data(target + 84 * ncomps + c);
        auto *t_85 = buffer.data(target + 85 * ncomps + c);
        auto *t_86 = buffer.data(target + 86 * ncomps + c);
        auto *t_87 = buffer.data(target + 87 * ncomps + c);
        auto *t_88 = buffer.data(target + 88 * ncomps + c);
        auto *t_89 = buffer.data(target + 89 * ncomps + c);
        auto *t_90 = buffer.data(target + 90 * ncomps + c);
        auto *t_91 = buffer.data(target + 91 * ncomps + c);
        auto *t_92 = buffer.data(target + 92 * ncomps + c);
        auto *t_93 = buffer.data(target + 93 * ncomps + c);
        auto *t_94 = buffer.data(target + 94 * ncomps + c);
        auto *t_95 = buffer.data(target + 95 * ncomps + c);
        auto *t_96 = buffer.data(target + 96 * ncomps + c);
        auto *t_97 = buffer.data(target + 97 * ncomps + c);
        auto *t_98 = buffer.data(target + 98 * ncomps + c);
        auto *t_99 = buffer.data(target + 99 * ncomps + c);
        auto *t_100 = buffer.data(target + 100 * ncomps + c);
        auto *t_101 = buffer.data(target + 101 * ncomps + c);
        auto *t_102 = buffer.data(target + 102 * ncomps + c);
        auto *t_103 = buffer.data(target + 103 * ncomps + c);
        auto *t_104 = buffer.data(target + 104 * ncomps + c);
        auto *t_105 = buffer.data(target + 105 * ncomps + c);
        auto *t_106 = buffer.data(target + 106 * ncomps + c);
        auto *t_107 = buffer.data(target + 107 * ncomps + c);
        auto *t_108 = buffer.data(target + 108 * ncomps + c);
        auto *t_109 = buffer.data(target + 109 * ncomps + c);
        auto *t_110 = buffer.data(target + 110 * ncomps + c);
        auto *t_111 = buffer.data(target + 111 * ncomps + c);
        auto *t_112 = buffer.data(target + 112 * ncomps + c);
        auto *t_113 = buffer.data(target + 113 * ncomps + c);
        auto *t_114 = buffer.data(target + 114 * ncomps + c);
        auto *t_115 = buffer.data(target + 115 * ncomps + c);
        auto *t_116 = buffer.data(target + 116 * ncomps + c);
        auto *t_117 = buffer.data(target + 117 * ncomps + c);
        auto *t_118 = buffer.data(target + 118 * ncomps + c);
        auto *t_119 = buffer.data(target + 119 * ncomps + c);
        auto *t_120 = buffer.data(target + 120 * ncomps + c);
        auto *t_121 = buffer.data(target + 121 * ncomps + c);
        auto *t_122 = buffer.data(target + 122 * ncomps + c);
        auto *t_123 = buffer.data(target + 123 * ncomps + c);
        auto *t_124 = buffer.data(target + 124 * ncomps + c);
        auto *t_125 = buffer.data(target + 125 * ncomps + c);
        auto *t_126 = buffer.data(target + 126 * ncomps + c);
        auto *t_127 = buffer.data(target + 127 * ncomps + c);
        auto *t_128 = buffer.data(target + 128 * ncomps + c);
        auto *t_129 = buffer.data(target + 129 * ncomps + c);
        auto *t_130 = buffer.data(target + 130 * ncomps + c);
        auto *t_131 = buffer.data(target + 131 * ncomps + c);
        auto *t_132 = buffer.data(target + 132 * ncomps + c);
        auto *t_133 = buffer.data(target + 133 * ncomps + c);
        auto *t_134 = buffer.data(target + 134 * ncomps + c);
        auto *t_135 = buffer.data(target + 135 * ncomps + c);
        auto *t_136 = buffer.data(target + 136 * ncomps + c);
        auto *t_137 = buffer.data(target + 137 * ncomps + c);
        auto *t_138 = buffer.data(target + 138 * ncomps + c);
        auto *t_139 = buffer.data(target + 139 * ncomps + c);
        auto *t_140 = buffer.data(target + 140 * ncomps + c);
        auto *t_141 = buffer.data(target + 141 * ncomps + c);
        auto *t_142 = buffer.data(target + 142 * ncomps + c);
        auto *t_143 = buffer.data(target + 143 * ncomps + c);
        auto *t_144 = buffer.data(target + 144 * ncomps + c);
        auto *t_145 = buffer.data(target + 145 * ncomps + c);
        auto *t_146 = buffer.data(target + 146 * ncomps + c);
        auto *t_147 = buffer.data(target + 147 * ncomps + c);
        auto *t_148 = buffer.data(target + 148 * ncomps + c);
        auto *t_149 = buffer.data(target + 149 * ncomps + c);
        auto *t_150 = buffer.data(target + 150 * ncomps + c);
        auto *t_151 = buffer.data(target + 151 * ncomps + c);
        auto *t_152 = buffer.data(target + 152 * ncomps + c);
        auto *t_153 = buffer.data(target + 153 * ncomps + c);
        auto *t_154 = buffer.data(target + 154 * ncomps + c);
        auto *t_155 = buffer.data(target + 155 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ns_0 = buffer.data(ns + 0 * ncomps + c);
        const auto *ns_1 = buffer.data(ns + 1 * ncomps + c);
        const auto *ns_2 = buffer.data(ns + 2 * ncomps + c);
        const auto *ns_3 = buffer.data(ns + 3 * ncomps + c);
        const auto *ns_4 = buffer.data(ns + 4 * ncomps + c);
        const auto *ns_5 = buffer.data(ns + 5 * ncomps + c);
        const auto *ns_6 = buffer.data(ns + 6 * ncomps + c);
        const auto *ns_7 = buffer.data(ns + 7 * ncomps + c);
        const auto *ns_8 = buffer.data(ns + 8 * ncomps + c);
        const auto *ns_9 = buffer.data(ns + 9 * ncomps + c);
        const auto *ns_10 = buffer.data(ns + 10 * ncomps + c);
        const auto *ns_11 = buffer.data(ns + 11 * ncomps + c);
        const auto *ns_12 = buffer.data(ns + 12 * ncomps + c);
        const auto *ns_13 = buffer.data(ns + 13 * ncomps + c);
        const auto *ns_14 = buffer.data(ns + 14 * ncomps + c);
        const auto *ns_15 = buffer.data(ns + 15 * ncomps + c);
        const auto *ns_16 = buffer.data(ns + 16 * ncomps + c);
        const auto *ns_17 = buffer.data(ns + 17 * ncomps + c);
        const auto *ns_18 = buffer.data(ns + 18 * ncomps + c);
        const auto *ns_19 = buffer.data(ns + 19 * ncomps + c);
        const auto *ns_20 = buffer.data(ns + 20 * ncomps + c);
        const auto *ns_21 = buffer.data(ns + 21 * ncomps + c);
        const auto *ns_22 = buffer.data(ns + 22 * ncomps + c);
        const auto *ns_23 = buffer.data(ns + 23 * ncomps + c);
        const auto *ns_24 = buffer.data(ns + 24 * ncomps + c);
        const auto *ns_25 = buffer.data(ns + 25 * ncomps + c);
        const auto *ns_26 = buffer.data(ns + 26 * ncomps + c);
        const auto *ns_27 = buffer.data(ns + 27 * ncomps + c);
        const auto *ns_28 = buffer.data(ns + 28 * ncomps + c);
        const auto *ns_29 = buffer.data(ns + 29 * ncomps + c);
        const auto *ns_30 = buffer.data(ns + 30 * ncomps + c);
        const auto *ns_31 = buffer.data(ns + 31 * ncomps + c);
        const auto *ns_32 = buffer.data(ns + 32 * ncomps + c);
        const auto *ns_33 = buffer.data(ns + 33 * ncomps + c);
        const auto *ns_34 = buffer.data(ns + 34 * ncomps + c);
        const auto *ns_35 = buffer.data(ns + 35 * ncomps + c);
        const auto *ns_36 = buffer.data(ns + 36 * ncomps + c);
        const auto *ns_37 = buffer.data(ns + 37 * ncomps + c);
        const auto *ns_38 = buffer.data(ns + 38 * ncomps + c);
        const auto *ns_39 = buffer.data(ns + 39 * ncomps + c);
        const auto *ns_40 = buffer.data(ns + 40 * ncomps + c);
        const auto *ns_41 = buffer.data(ns + 41 * ncomps + c);
        const auto *ns_42 = buffer.data(ns + 42 * ncomps + c);
        const auto *ns_43 = buffer.data(ns + 43 * ncomps + c);
        const auto *ns_44 = buffer.data(ns + 44 * ncomps + c);
        const auto *ns_45 = buffer.data(ns + 45 * ncomps + c);
        const auto *ns_46 = buffer.data(ns + 46 * ncomps + c);
        const auto *ns_47 = buffer.data(ns + 47 * ncomps + c);
        const auto *ns_48 = buffer.data(ns + 48 * ncomps + c);
        const auto *ns_49 = buffer.data(ns + 49 * ncomps + c);
        const auto *ns_50 = buffer.data(ns + 50 * ncomps + c);
        const auto *ns_51 = buffer.data(ns + 51 * ncomps + c);

        const auto *os_0 = buffer.data(os + 0 * ncomps + c);
        const auto *os_1 = buffer.data(os + 1 * ncomps + c);
        const auto *os_2 = buffer.data(os + 2 * ncomps + c);
        const auto *os_3 = buffer.data(os + 3 * ncomps + c);
        const auto *os_4 = buffer.data(os + 4 * ncomps + c);
        const auto *os_5 = buffer.data(os + 5 * ncomps + c);
        const auto *os_6 = buffer.data(os + 6 * ncomps + c);
        const auto *os_7 = buffer.data(os + 7 * ncomps + c);
        const auto *os_8 = buffer.data(os + 8 * ncomps + c);
        const auto *os_9 = buffer.data(os + 9 * ncomps + c);
        const auto *os_10 = buffer.data(os + 10 * ncomps + c);
        const auto *os_11 = buffer.data(os + 11 * ncomps + c);
        const auto *os_12 = buffer.data(os + 12 * ncomps + c);
        const auto *os_13 = buffer.data(os + 13 * ncomps + c);
        const auto *os_14 = buffer.data(os + 14 * ncomps + c);
        const auto *os_15 = buffer.data(os + 15 * ncomps + c);
        const auto *os_16 = buffer.data(os + 16 * ncomps + c);
        const auto *os_17 = buffer.data(os + 17 * ncomps + c);
        const auto *os_18 = buffer.data(os + 18 * ncomps + c);
        const auto *os_19 = buffer.data(os + 19 * ncomps + c);
        const auto *os_20 = buffer.data(os + 20 * ncomps + c);
        const auto *os_21 = buffer.data(os + 21 * ncomps + c);
        const auto *os_22 = buffer.data(os + 22 * ncomps + c);
        const auto *os_23 = buffer.data(os + 23 * ncomps + c);
        const auto *os_24 = buffer.data(os + 24 * ncomps + c);
        const auto *os_25 = buffer.data(os + 25 * ncomps + c);
        const auto *os_26 = buffer.data(os + 26 * ncomps + c);
        const auto *os_27 = buffer.data(os + 27 * ncomps + c);
        const auto *os_28 = buffer.data(os + 28 * ncomps + c);
        const auto *os_29 = buffer.data(os + 29 * ncomps + c);
        const auto *os_30 = buffer.data(os + 30 * ncomps + c);
        const auto *os_31 = buffer.data(os + 31 * ncomps + c);
        const auto *os_32 = buffer.data(os + 32 * ncomps + c);
        const auto *os_33 = buffer.data(os + 33 * ncomps + c);
        const auto *os_34 = buffer.data(os + 34 * ncomps + c);
        const auto *os_35 = buffer.data(os + 35 * ncomps + c);
        const auto *os_36 = buffer.data(os + 36 * ncomps + c);
        const auto *os_37 = buffer.data(os + 37 * ncomps + c);
        const auto *os_38 = buffer.data(os + 38 * ncomps + c);
        const auto *os_39 = buffer.data(os + 39 * ncomps + c);
        const auto *os_40 = buffer.data(os + 40 * ncomps + c);
        const auto *os_41 = buffer.data(os + 41 * ncomps + c);
        const auto *os_42 = buffer.data(os + 42 * ncomps + c);
        const auto *os_43 = buffer.data(os + 43 * ncomps + c);
        const auto *os_44 = buffer.data(os + 44 * ncomps + c);
        const auto *os_45 = buffer.data(os + 45 * ncomps + c);
        const auto *os_46 = buffer.data(os + 46 * ncomps + c);
        const auto *os_47 = buffer.data(os + 47 * ncomps + c);
        const auto *os_48 = buffer.data(os + 48 * ncomps + c);
        const auto *os_49 = buffer.data(os + 49 * ncomps + c);
        const auto *os_50 = buffer.data(os + 50 * ncomps + c);
        const auto *os_51 = buffer.data(os + 51 * ncomps + c);
        const auto *os_52 = buffer.data(os + 52 * ncomps + c);
        const auto *os_53 = buffer.data(os + 53 * ncomps + c);
        const auto *os_54 = buffer.data(os + 54 * ncomps + c);
        const auto *os_55 = buffer.data(os + 55 * ncomps + c);
        const auto *os_56 = buffer.data(os + 56 * ncomps + c);
        const auto *os_57 = buffer.data(os + 57 * ncomps + c);
        const auto *os_58 = buffer.data(os + 58 * ncomps + c);
        const auto *os_59 = buffer.data(os + 59 * ncomps + c);
        const auto *os_60 = buffer.data(os + 60 * ncomps + c);
        const auto *os_61 = buffer.data(os + 61 * ncomps + c);
        const auto *os_62 = buffer.data(os + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, ns_0, ns_1, os_0, \
                         os_1, os_2, os_3, os_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ns_0[k]
                     + os_0[k];

            t_1[k] = ab_y[k] * ns_0[k]
                     + os_1[k];

            t_2[k] = ab_z[k] * ns_0[k]
                     + os_2[k];

            t_3[k] = ab_x[k] * ns_1[k]
                     + os_1[k];

            t_4[k] = ab_y[k] * ns_1[k]
                     + os_3[k];

            t_5[k] = ab_z[k] * ns_1[k]
                     + os_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, ns_2, ns_3, os_2, os_3, \
                         os_4, os_5, os_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * ns_2[k]
                     + os_2[k];

            t_7[k] = ab_y[k] * ns_2[k]
                     + os_4[k];

            t_8[k] = ab_z[k] * ns_2[k]
                     + os_5[k];

            t_9[k] = ab_x[k] * ns_3[k]
                     + os_3[k];

            t_10[k] = ab_y[k] * ns_3[k]
                      + os_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, ns_3, ns_4, \
                         ns_5, os_4, os_5, os_7, os_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * ns_3[k]
                      + os_7[k];

            t_12[k] = ab_x[k] * ns_4[k]
                      + os_4[k];

            t_13[k] = ab_y[k] * ns_4[k]
                      + os_7[k];

            t_14[k] = ab_z[k] * ns_4[k]
                      + os_8[k];

            t_15[k] = ab_x[k] * ns_5[k]
                      + os_5[k];

            t_16[k] = ab_y[k] * ns_5[k]
                      + os_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, ns_5, ns_6, ns_7, \
                         os_6, os_7, os_9, os_10, os_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * ns_5[k]
                      + os_9[k];

            t_18[k] = ab_x[k] * ns_6[k]
                      + os_6[k];

            t_19[k] = ab_y[k] * ns_6[k]
                      + os_10[k];

            t_20[k] = ab_z[k] * ns_6[k]
                      + os_11[k];

            t_21[k] = ab_x[k] * ns_7[k]
                      + os_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, ns_7, ns_8, os_8, \
                         os_11, os_12, os_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * ns_7[k]
                      + os_11[k];

            t_23[k] = ab_z[k] * ns_7[k]
                      + os_12[k];

            t_24[k] = ab_x[k] * ns_8[k]
                      + os_8[k];

            t_25[k] = ab_y[k] * ns_8[k]
                      + os_12[k];

            t_26[k] = ab_z[k] * ns_8[k]
                      + os_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ns_9, ns_10, os_9, \
                         os_10, os_13, os_14, os_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * ns_9[k]
                      + os_9[k];

            t_28[k] = ab_y[k] * ns_9[k]
                      + os_13[k];

            t_29[k] = ab_z[k] * ns_9[k]
                      + os_14[k];

            t_30[k] = ab_x[k] * ns_10[k]
                      + os_10[k];

            t_31[k] = ab_y[k] * ns_10[k]
                      + os_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, ns_10, ns_11, \
                         ns_12, os_11, os_12, os_16, os_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * ns_10[k]
                      + os_16[k];

            t_33[k] = ab_x[k] * ns_11[k]
                      + os_11[k];

            t_34[k] = ab_y[k] * ns_11[k]
                      + os_16[k];

            t_35[k] = ab_z[k] * ns_11[k]
                      + os_17[k];

            t_36[k] = ab_x[k] * ns_12[k]
                      + os_12[k];

            t_37[k] = ab_y[k] * ns_12[k]
                      + os_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, ns_12, ns_13, \
                         ns_14, os_13, os_14, os_18, os_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * ns_12[k]
                      + os_18[k];

            t_39[k] = ab_x[k] * ns_13[k]
                      + os_13[k];

            t_40[k] = ab_y[k] * ns_13[k]
                      + os_18[k];

            t_41[k] = ab_z[k] * ns_13[k]
                      + os_19[k];

            t_42[k] = ab_x[k] * ns_14[k]
                      + os_14[k];

            t_43[k] = ab_y[k] * ns_14[k]
                      + os_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, ns_14, ns_15, ns_16, \
                         os_15, os_16, os_20, os_21, os_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * ns_14[k]
                      + os_20[k];

            t_45[k] = ab_x[k] * ns_15[k]
                      + os_15[k];

            t_46[k] = ab_y[k] * ns_15[k]
                      + os_21[k];

            t_47[k] = ab_z[k] * ns_15[k]
                      + os_22[k];

            t_48[k] = ab_x[k] * ns_16[k]
                      + os_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ns_16, ns_17, os_17, \
                         os_22, os_23, os_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * ns_16[k]
                      + os_22[k];

            t_50[k] = ab_z[k] * ns_16[k]
                      + os_23[k];

            t_51[k] = ab_x[k] * ns_17[k]
                      + os_17[k];

            t_52[k] = ab_y[k] * ns_17[k]
                      + os_23[k];

            t_53[k] = ab_z[k] * ns_17[k]
                      + os_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ns_18, ns_19, \
                         os_18, os_19, os_24, os_25, os_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * ns_18[k]
                      + os_18[k];

            t_55[k] = ab_y[k] * ns_18[k]
                      + os_24[k];

            t_56[k] = ab_z[k] * ns_18[k]
                      + os_25[k];

            t_57[k] = ab_x[k] * ns_19[k]
                      + os_19[k];

            t_58[k] = ab_y[k] * ns_19[k]
                      + os_25[k];

            t_59[k] = ab_z[k] * ns_19[k]
                      + os_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, ns_20, ns_21, os_20, \
                         os_21, os_26, os_27, os_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ns_20[k]
                      + os_20[k];

            t_61[k] = ab_y[k] * ns_20[k]
                      + os_26[k];

            t_62[k] = ab_z[k] * ns_20[k]
                      + os_27[k];

            t_63[k] = ab_x[k] * ns_21[k]
                      + os_21[k];

            t_64[k] = ab_y[k] * ns_21[k]
                      + os_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, ns_21, ns_22, \
                         ns_23, os_22, os_23, os_29, os_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * ns_21[k]
                      + os_29[k];

            t_66[k] = ab_x[k] * ns_22[k]
                      + os_22[k];

            t_67[k] = ab_y[k] * ns_22[k]
                      + os_29[k];

            t_68[k] = ab_z[k] * ns_22[k]
                      + os_30[k];

            t_69[k] = ab_x[k] * ns_23[k]
                      + os_23[k];

            t_70[k] = ab_y[k] * ns_23[k]
                      + os_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, ns_23, ns_24, \
                         ns_25, os_24, os_25, os_31, os_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * ns_23[k]
                      + os_31[k];

            t_72[k] = ab_x[k] * ns_24[k]
                      + os_24[k];

            t_73[k] = ab_y[k] * ns_24[k]
                      + os_31[k];

            t_74[k] = ab_z[k] * ns_24[k]
                      + os_32[k];

            t_75[k] = ab_x[k] * ns_25[k]
                      + os_25[k];

            t_76[k] = ab_y[k] * ns_25[k]
                      + os_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, ns_25, ns_26, \
                         ns_27, os_26, os_27, os_33, os_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * ns_25[k]
                      + os_33[k];

            t_78[k] = ab_x[k] * ns_26[k]
                      + os_26[k];

            t_79[k] = ab_y[k] * ns_26[k]
                      + os_33[k];

            t_80[k] = ab_z[k] * ns_26[k]
                      + os_34[k];

            t_81[k] = ab_x[k] * ns_27[k]
                      + os_27[k];

            t_82[k] = ab_y[k] * ns_27[k]
                      + os_34[k];
        }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ab_z, ns_27, ns_28, ns_29, \
                         os_28, os_29, os_35, os_36, os_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * ns_27[k]
                      + os_35[k];

            t_84[k] = ab_x[k] * ns_28[k]
                      + os_28[k];

            t_85[k] = ab_y[k] * ns_28[k]
                      + os_36[k];

            t_86[k] = ab_z[k] * ns_28[k]
                      + os_37[k];

            t_87[k] = ab_x[k] * ns_29[k]
                      + os_29[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ab_x, ab_y, ab_z, ns_29, ns_30, os_30, \
                         os_37, os_38, os_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * ns_29[k]
                      + os_37[k];

            t_89[k] = ab_z[k] * ns_29[k]
                      + os_38[k];

            t_90[k] = ab_x[k] * ns_30[k]
                      + os_30[k];

            t_91[k] = ab_y[k] * ns_30[k]
                      + os_38[k];

            t_92[k] = ab_z[k] * ns_30[k]
                      + os_39[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, ns_31, ns_32, \
                         os_31, os_32, os_39, os_40, os_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = ab_x[k] * ns_31[k]
                      + os_31[k];

            t_94[k] = ab_y[k] * ns_31[k]
                      + os_39[k];

            t_95[k] = ab_z[k] * ns_31[k]
                      + os_40[k];

            t_96[k] = ab_x[k] * ns_32[k]
                      + os_32[k];

            t_97[k] = ab_y[k] * ns_32[k]
                      + os_40[k];

            t_98[k] = ab_z[k] * ns_32[k]
                      + os_41[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, ab_x, ab_y, ab_z, ns_33, \
                         ns_34, os_33, os_34, os_41, os_42, os_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * ns_33[k]
                      + os_33[k];

            t_100[k] = ab_y[k] * ns_33[k]
                       + os_41[k];

            t_101[k] = ab_z[k] * ns_33[k]
                       + os_42[k];

            t_102[k] = ab_x[k] * ns_34[k]
                       + os_34[k];

            t_103[k] = ab_y[k] * ns_34[k]
                       + os_42[k];

            t_104[k] = ab_z[k] * ns_34[k]
                       + os_43[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, ns_35, ns_36, \
                         os_35, os_36, os_43, os_44, os_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * ns_35[k]
                       + os_35[k];

            t_106[k] = ab_y[k] * ns_35[k]
                       + os_43[k];

            t_107[k] = ab_z[k] * ns_35[k]
                       + os_44[k];

            t_108[k] = ab_x[k] * ns_36[k]
                       + os_36[k];

            t_109[k] = ab_y[k] * ns_36[k]
                       + os_45[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, ab_x, ab_y, ab_z, ns_36, \
                         ns_37, ns_38, os_37, os_38, os_46, os_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_z[k] * ns_36[k]
                       + os_46[k];

            t_111[k] = ab_x[k] * ns_37[k]
                       + os_37[k];

            t_112[k] = ab_y[k] * ns_37[k]
                       + os_46[k];

            t_113[k] = ab_z[k] * ns_37[k]
                       + os_47[k];

            t_114[k] = ab_x[k] * ns_38[k]
                       + os_38[k];

            t_115[k] = ab_y[k] * ns_38[k]
                       + os_47[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, ab_x, ab_y, ab_z, ns_38, \
                         ns_39, ns_40, os_39, os_40, os_48, os_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_z[k] * ns_38[k]
                       + os_48[k];

            t_117[k] = ab_x[k] * ns_39[k]
                       + os_39[k];

            t_118[k] = ab_y[k] * ns_39[k]
                       + os_48[k];

            t_119[k] = ab_z[k] * ns_39[k]
                       + os_49[k];

            t_120[k] = ab_x[k] * ns_40[k]
                       + os_40[k];

            t_121[k] = ab_y[k] * ns_40[k]
                       + os_49[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, ab_x, ab_y, ab_z, ns_40, \
                         ns_41, ns_42, os_41, os_42, os_50, os_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_z[k] * ns_40[k]
                       + os_50[k];

            t_123[k] = ab_x[k] * ns_41[k]
                       + os_41[k];

            t_124[k] = ab_y[k] * ns_41[k]
                       + os_50[k];

            t_125[k] = ab_z[k] * ns_41[k]
                       + os_51[k];

            t_126[k] = ab_x[k] * ns_42[k]
                       + os_42[k];

            t_127[k] = ab_y[k] * ns_42[k]
                       + os_51[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, ab_x, ab_y, ab_z, ns_42, \
                         ns_43, ns_44, os_43, os_44, os_52, os_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = ab_z[k] * ns_42[k]
                       + os_52[k];

            t_129[k] = ab_x[k] * ns_43[k]
                       + os_43[k];

            t_130[k] = ab_y[k] * ns_43[k]
                       + os_52[k];

            t_131[k] = ab_z[k] * ns_43[k]
                       + os_53[k];

            t_132[k] = ab_x[k] * ns_44[k]
                       + os_44[k];

            t_133[k] = ab_y[k] * ns_44[k]
                       + os_53[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, ns_44, ns_45, \
                         ns_46, os_45, os_46, os_54, os_55, os_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = ab_z[k] * ns_44[k]
                       + os_54[k];

            t_135[k] = ab_x[k] * ns_45[k]
                       + os_45[k];

            t_136[k] = ab_y[k] * ns_45[k]
                       + os_55[k];

            t_137[k] = ab_z[k] * ns_45[k]
                       + os_56[k];

            t_138[k] = ab_x[k] * ns_46[k]
                       + os_46[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, ns_46, ns_47, \
                         os_47, os_56, os_57, os_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_y[k] * ns_46[k]
                       + os_56[k];

            t_140[k] = ab_z[k] * ns_46[k]
                       + os_57[k];

            t_141[k] = ab_x[k] * ns_47[k]
                       + os_47[k];

            t_142[k] = ab_y[k] * ns_47[k]
                       + os_57[k];

            t_143[k] = ab_z[k] * ns_47[k]
                       + os_58[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, ns_48, \
                         ns_49, os_48, os_49, os_58, os_59, os_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * ns_48[k]
                       + os_48[k];

            t_145[k] = ab_y[k] * ns_48[k]
                       + os_58[k];

            t_146[k] = ab_z[k] * ns_48[k]
                       + os_59[k];

            t_147[k] = ab_x[k] * ns_49[k]
                       + os_49[k];

            t_148[k] = ab_y[k] * ns_49[k]
                       + os_59[k];

            t_149[k] = ab_z[k] * ns_49[k]
                       + os_60[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, ab_x, ab_y, ab_z, ns_50, \
                         ns_51, os_50, os_51, os_60, os_61, os_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ns_50[k]
                       + os_50[k];

            t_151[k] = ab_y[k] * ns_50[k]
                       + os_60[k];

            t_152[k] = ab_z[k] * ns_50[k]
                       + os_61[k];

            t_153[k] = ab_x[k] * ns_51[k]
                       + os_51[k];

            t_154[k] = ab_y[k] * ns_51[k]
                       + os_61[k];

            t_155[k] = ab_z[k] * ns_51[k]
                       + os_62[k];
        }
    }
}

static auto
compute_hrr_np_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ns, const size_t os, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_156 = buffer.data(target + 156 * ncomps + c);
        auto *t_157 = buffer.data(target + 157 * ncomps + c);
        auto *t_158 = buffer.data(target + 158 * ncomps + c);
        auto *t_159 = buffer.data(target + 159 * ncomps + c);
        auto *t_160 = buffer.data(target + 160 * ncomps + c);
        auto *t_161 = buffer.data(target + 161 * ncomps + c);
        auto *t_162 = buffer.data(target + 162 * ncomps + c);
        auto *t_163 = buffer.data(target + 163 * ncomps + c);
        auto *t_164 = buffer.data(target + 164 * ncomps + c);
        auto *t_165 = buffer.data(target + 165 * ncomps + c);
        auto *t_166 = buffer.data(target + 166 * ncomps + c);
        auto *t_167 = buffer.data(target + 167 * ncomps + c);
        auto *t_168 = buffer.data(target + 168 * ncomps + c);
        auto *t_169 = buffer.data(target + 169 * ncomps + c);
        auto *t_170 = buffer.data(target + 170 * ncomps + c);
        auto *t_171 = buffer.data(target + 171 * ncomps + c);
        auto *t_172 = buffer.data(target + 172 * ncomps + c);
        auto *t_173 = buffer.data(target + 173 * ncomps + c);
        auto *t_174 = buffer.data(target + 174 * ncomps + c);
        auto *t_175 = buffer.data(target + 175 * ncomps + c);
        auto *t_176 = buffer.data(target + 176 * ncomps + c);
        auto *t_177 = buffer.data(target + 177 * ncomps + c);
        auto *t_178 = buffer.data(target + 178 * ncomps + c);
        auto *t_179 = buffer.data(target + 179 * ncomps + c);
        auto *t_180 = buffer.data(target + 180 * ncomps + c);
        auto *t_181 = buffer.data(target + 181 * ncomps + c);
        auto *t_182 = buffer.data(target + 182 * ncomps + c);
        auto *t_183 = buffer.data(target + 183 * ncomps + c);
        auto *t_184 = buffer.data(target + 184 * ncomps + c);
        auto *t_185 = buffer.data(target + 185 * ncomps + c);
        auto *t_186 = buffer.data(target + 186 * ncomps + c);
        auto *t_187 = buffer.data(target + 187 * ncomps + c);
        auto *t_188 = buffer.data(target + 188 * ncomps + c);
        auto *t_189 = buffer.data(target + 189 * ncomps + c);
        auto *t_190 = buffer.data(target + 190 * ncomps + c);
        auto *t_191 = buffer.data(target + 191 * ncomps + c);
        auto *t_192 = buffer.data(target + 192 * ncomps + c);
        auto *t_193 = buffer.data(target + 193 * ncomps + c);
        auto *t_194 = buffer.data(target + 194 * ncomps + c);
        auto *t_195 = buffer.data(target + 195 * ncomps + c);
        auto *t_196 = buffer.data(target + 196 * ncomps + c);
        auto *t_197 = buffer.data(target + 197 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ns_52 = buffer.data(ns + 52 * ncomps + c);
        const auto *ns_53 = buffer.data(ns + 53 * ncomps + c);
        const auto *ns_54 = buffer.data(ns + 54 * ncomps + c);
        const auto *ns_55 = buffer.data(ns + 55 * ncomps + c);
        const auto *ns_56 = buffer.data(ns + 56 * ncomps + c);
        const auto *ns_57 = buffer.data(ns + 57 * ncomps + c);
        const auto *ns_58 = buffer.data(ns + 58 * ncomps + c);
        const auto *ns_59 = buffer.data(ns + 59 * ncomps + c);
        const auto *ns_60 = buffer.data(ns + 60 * ncomps + c);
        const auto *ns_61 = buffer.data(ns + 61 * ncomps + c);
        const auto *ns_62 = buffer.data(ns + 62 * ncomps + c);
        const auto *ns_63 = buffer.data(ns + 63 * ncomps + c);
        const auto *ns_64 = buffer.data(ns + 64 * ncomps + c);
        const auto *ns_65 = buffer.data(ns + 65 * ncomps + c);

        const auto *os_52 = buffer.data(os + 52 * ncomps + c);
        const auto *os_53 = buffer.data(os + 53 * ncomps + c);
        const auto *os_54 = buffer.data(os + 54 * ncomps + c);
        const auto *os_55 = buffer.data(os + 55 * ncomps + c);
        const auto *os_56 = buffer.data(os + 56 * ncomps + c);
        const auto *os_57 = buffer.data(os + 57 * ncomps + c);
        const auto *os_58 = buffer.data(os + 58 * ncomps + c);
        const auto *os_59 = buffer.data(os + 59 * ncomps + c);
        const auto *os_60 = buffer.data(os + 60 * ncomps + c);
        const auto *os_61 = buffer.data(os + 61 * ncomps + c);
        const auto *os_62 = buffer.data(os + 62 * ncomps + c);
        const auto *os_63 = buffer.data(os + 63 * ncomps + c);
        const auto *os_64 = buffer.data(os + 64 * ncomps + c);
        const auto *os_65 = buffer.data(os + 65 * ncomps + c);
        const auto *os_66 = buffer.data(os + 66 * ncomps + c);
        const auto *os_67 = buffer.data(os + 67 * ncomps + c);
        const auto *os_68 = buffer.data(os + 68 * ncomps + c);
        const auto *os_69 = buffer.data(os + 69 * ncomps + c);
        const auto *os_70 = buffer.data(os + 70 * ncomps + c);
        const auto *os_71 = buffer.data(os + 71 * ncomps + c);
        const auto *os_72 = buffer.data(os + 72 * ncomps + c);
        const auto *os_73 = buffer.data(os + 73 * ncomps + c);
        const auto *os_74 = buffer.data(os + 74 * ncomps + c);
        const auto *os_75 = buffer.data(os + 75 * ncomps + c);
        const auto *os_76 = buffer.data(os + 76 * ncomps + c);
        const auto *os_77 = buffer.data(os + 77 * ncomps + c);

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, ns_52, \
                         ns_53, os_52, os_53, os_62, os_63, os_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = ab_x[k] * ns_52[k]
                       + os_52[k];

            t_157[k] = ab_y[k] * ns_52[k]
                       + os_62[k];

            t_158[k] = ab_z[k] * ns_52[k]
                       + os_63[k];

            t_159[k] = ab_x[k] * ns_53[k]
                       + os_53[k];

            t_160[k] = ab_y[k] * ns_53[k]
                       + os_63[k];

            t_161[k] = ab_z[k] * ns_53[k]
                       + os_64[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, ab_z, ns_54, ns_55, \
                         os_54, os_55, os_64, os_65, os_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * ns_54[k]
                       + os_54[k];

            t_163[k] = ab_y[k] * ns_54[k]
                       + os_64[k];

            t_164[k] = ab_z[k] * ns_54[k]
                       + os_65[k];

            t_165[k] = ab_x[k] * ns_55[k]
                       + os_55[k];

            t_166[k] = ab_y[k] * ns_55[k]
                       + os_66[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, t_172, ab_x, ab_y, ab_z, ns_55, \
                         ns_56, ns_57, os_56, os_57, os_67, os_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = ab_z[k] * ns_55[k]
                       + os_67[k];

            t_168[k] = ab_x[k] * ns_56[k]
                       + os_56[k];

            t_169[k] = ab_y[k] * ns_56[k]
                       + os_67[k];

            t_170[k] = ab_z[k] * ns_56[k]
                       + os_68[k];

            t_171[k] = ab_x[k] * ns_57[k]
                       + os_57[k];

            t_172[k] = ab_y[k] * ns_57[k]
                       + os_68[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, ab_x, ab_y, ab_z, ns_57, \
                         ns_58, ns_59, os_58, os_59, os_69, os_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_z[k] * ns_57[k]
                       + os_69[k];

            t_174[k] = ab_x[k] * ns_58[k]
                       + os_58[k];

            t_175[k] = ab_y[k] * ns_58[k]
                       + os_69[k];

            t_176[k] = ab_z[k] * ns_58[k]
                       + os_70[k];

            t_177[k] = ab_x[k] * ns_59[k]
                       + os_59[k];

            t_178[k] = ab_y[k] * ns_59[k]
                       + os_70[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, ab_z, ns_59, \
                         ns_60, ns_61, os_60, os_61, os_71, os_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = ab_z[k] * ns_59[k]
                       + os_71[k];

            t_180[k] = ab_x[k] * ns_60[k]
                       + os_60[k];

            t_181[k] = ab_y[k] * ns_60[k]
                       + os_71[k];

            t_182[k] = ab_z[k] * ns_60[k]
                       + os_72[k];

            t_183[k] = ab_x[k] * ns_61[k]
                       + os_61[k];

            t_184[k] = ab_y[k] * ns_61[k]
                       + os_72[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, ab_x, ab_y, ab_z, ns_61, \
                         ns_62, ns_63, os_62, os_63, os_73, os_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_z[k] * ns_61[k]
                       + os_73[k];

            t_186[k] = ab_x[k] * ns_62[k]
                       + os_62[k];

            t_187[k] = ab_y[k] * ns_62[k]
                       + os_73[k];

            t_188[k] = ab_z[k] * ns_62[k]
                       + os_74[k];

            t_189[k] = ab_x[k] * ns_63[k]
                       + os_63[k];

            t_190[k] = ab_y[k] * ns_63[k]
                       + os_74[k];
        }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, ab_x, ab_y, ab_z, ns_63, \
                         ns_64, ns_65, os_64, os_65, os_75, os_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_191[k] = ab_z[k] * ns_63[k]
                       + os_75[k];

            t_192[k] = ab_x[k] * ns_64[k]
                       + os_64[k];

            t_193[k] = ab_y[k] * ns_64[k]
                       + os_75[k];

            t_194[k] = ab_z[k] * ns_64[k]
                       + os_76[k];

            t_195[k] = ab_x[k] * ns_65[k]
                       + os_65[k];

            t_196[k] = ab_y[k] * ns_65[k]
                       + os_76[k];
        }

#pragma omp simd aligned(t_197, ab_z, ns_65, os_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_z[k] * ns_65[k]
                       + os_77[k];
        }
    }
}

auto
compute_hrr_np(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ns, const size_t os, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_np_piece0(buffer, coordinates, target, ns, os, ncomps, nmax);

    compute_hrr_np_piece1(buffer, coordinates, target, ns, os, ncomps, nmax);
}

}  // namespace simdtrf
