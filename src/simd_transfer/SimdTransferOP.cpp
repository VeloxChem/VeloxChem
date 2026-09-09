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


#include "SimdTransferOP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_op_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t os, const size_t qs,
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

        const auto *qs_0 = buffer.data(qs + 0 * ncomps + c);
        const auto *qs_1 = buffer.data(qs + 1 * ncomps + c);
        const auto *qs_2 = buffer.data(qs + 2 * ncomps + c);
        const auto *qs_3 = buffer.data(qs + 3 * ncomps + c);
        const auto *qs_4 = buffer.data(qs + 4 * ncomps + c);
        const auto *qs_5 = buffer.data(qs + 5 * ncomps + c);
        const auto *qs_6 = buffer.data(qs + 6 * ncomps + c);
        const auto *qs_7 = buffer.data(qs + 7 * ncomps + c);
        const auto *qs_8 = buffer.data(qs + 8 * ncomps + c);
        const auto *qs_9 = buffer.data(qs + 9 * ncomps + c);
        const auto *qs_10 = buffer.data(qs + 10 * ncomps + c);
        const auto *qs_11 = buffer.data(qs + 11 * ncomps + c);
        const auto *qs_12 = buffer.data(qs + 12 * ncomps + c);
        const auto *qs_13 = buffer.data(qs + 13 * ncomps + c);
        const auto *qs_14 = buffer.data(qs + 14 * ncomps + c);
        const auto *qs_15 = buffer.data(qs + 15 * ncomps + c);
        const auto *qs_16 = buffer.data(qs + 16 * ncomps + c);
        const auto *qs_17 = buffer.data(qs + 17 * ncomps + c);
        const auto *qs_18 = buffer.data(qs + 18 * ncomps + c);
        const auto *qs_19 = buffer.data(qs + 19 * ncomps + c);
        const auto *qs_20 = buffer.data(qs + 20 * ncomps + c);
        const auto *qs_21 = buffer.data(qs + 21 * ncomps + c);
        const auto *qs_22 = buffer.data(qs + 22 * ncomps + c);
        const auto *qs_23 = buffer.data(qs + 23 * ncomps + c);
        const auto *qs_24 = buffer.data(qs + 24 * ncomps + c);
        const auto *qs_25 = buffer.data(qs + 25 * ncomps + c);
        const auto *qs_26 = buffer.data(qs + 26 * ncomps + c);
        const auto *qs_27 = buffer.data(qs + 27 * ncomps + c);
        const auto *qs_28 = buffer.data(qs + 28 * ncomps + c);
        const auto *qs_29 = buffer.data(qs + 29 * ncomps + c);
        const auto *qs_30 = buffer.data(qs + 30 * ncomps + c);
        const auto *qs_31 = buffer.data(qs + 31 * ncomps + c);
        const auto *qs_32 = buffer.data(qs + 32 * ncomps + c);
        const auto *qs_33 = buffer.data(qs + 33 * ncomps + c);
        const auto *qs_34 = buffer.data(qs + 34 * ncomps + c);
        const auto *qs_35 = buffer.data(qs + 35 * ncomps + c);
        const auto *qs_36 = buffer.data(qs + 36 * ncomps + c);
        const auto *qs_37 = buffer.data(qs + 37 * ncomps + c);
        const auto *qs_38 = buffer.data(qs + 38 * ncomps + c);
        const auto *qs_39 = buffer.data(qs + 39 * ncomps + c);
        const auto *qs_40 = buffer.data(qs + 40 * ncomps + c);
        const auto *qs_41 = buffer.data(qs + 41 * ncomps + c);
        const auto *qs_42 = buffer.data(qs + 42 * ncomps + c);
        const auto *qs_43 = buffer.data(qs + 43 * ncomps + c);
        const auto *qs_44 = buffer.data(qs + 44 * ncomps + c);
        const auto *qs_45 = buffer.data(qs + 45 * ncomps + c);
        const auto *qs_46 = buffer.data(qs + 46 * ncomps + c);
        const auto *qs_47 = buffer.data(qs + 47 * ncomps + c);
        const auto *qs_48 = buffer.data(qs + 48 * ncomps + c);
        const auto *qs_49 = buffer.data(qs + 49 * ncomps + c);
        const auto *qs_50 = buffer.data(qs + 50 * ncomps + c);
        const auto *qs_51 = buffer.data(qs + 51 * ncomps + c);
        const auto *qs_52 = buffer.data(qs + 52 * ncomps + c);
        const auto *qs_53 = buffer.data(qs + 53 * ncomps + c);
        const auto *qs_54 = buffer.data(qs + 54 * ncomps + c);
        const auto *qs_55 = buffer.data(qs + 55 * ncomps + c);
        const auto *qs_56 = buffer.data(qs + 56 * ncomps + c);
        const auto *qs_57 = buffer.data(qs + 57 * ncomps + c);
        const auto *qs_58 = buffer.data(qs + 58 * ncomps + c);
        const auto *qs_59 = buffer.data(qs + 59 * ncomps + c);
        const auto *qs_60 = buffer.data(qs + 60 * ncomps + c);
        const auto *qs_61 = buffer.data(qs + 61 * ncomps + c);
        const auto *qs_62 = buffer.data(qs + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, os_0, os_1, qs_0, \
                         qs_1, qs_2, qs_3, qs_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * os_0[k]
                     + qs_0[k];

            t_1[k] = ab_y[k] * os_0[k]
                     + qs_1[k];

            t_2[k] = ab_z[k] * os_0[k]
                     + qs_2[k];

            t_3[k] = ab_x[k] * os_1[k]
                     + qs_1[k];

            t_4[k] = ab_y[k] * os_1[k]
                     + qs_3[k];

            t_5[k] = ab_z[k] * os_1[k]
                     + qs_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, os_2, os_3, qs_2, qs_3, \
                         qs_4, qs_5, qs_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * os_2[k]
                     + qs_2[k];

            t_7[k] = ab_y[k] * os_2[k]
                     + qs_4[k];

            t_8[k] = ab_z[k] * os_2[k]
                     + qs_5[k];

            t_9[k] = ab_x[k] * os_3[k]
                     + qs_3[k];

            t_10[k] = ab_y[k] * os_3[k]
                      + qs_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, os_3, os_4, \
                         os_5, qs_4, qs_5, qs_7, qs_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * os_3[k]
                      + qs_7[k];

            t_12[k] = ab_x[k] * os_4[k]
                      + qs_4[k];

            t_13[k] = ab_y[k] * os_4[k]
                      + qs_7[k];

            t_14[k] = ab_z[k] * os_4[k]
                      + qs_8[k];

            t_15[k] = ab_x[k] * os_5[k]
                      + qs_5[k];

            t_16[k] = ab_y[k] * os_5[k]
                      + qs_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, os_5, os_6, os_7, \
                         qs_6, qs_7, qs_9, qs_10, qs_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * os_5[k]
                      + qs_9[k];

            t_18[k] = ab_x[k] * os_6[k]
                      + qs_6[k];

            t_19[k] = ab_y[k] * os_6[k]
                      + qs_10[k];

            t_20[k] = ab_z[k] * os_6[k]
                      + qs_11[k];

            t_21[k] = ab_x[k] * os_7[k]
                      + qs_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, os_7, os_8, qs_8, \
                         qs_11, qs_12, qs_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * os_7[k]
                      + qs_11[k];

            t_23[k] = ab_z[k] * os_7[k]
                      + qs_12[k];

            t_24[k] = ab_x[k] * os_8[k]
                      + qs_8[k];

            t_25[k] = ab_y[k] * os_8[k]
                      + qs_12[k];

            t_26[k] = ab_z[k] * os_8[k]
                      + qs_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, os_9, os_10, qs_9, \
                         qs_10, qs_13, qs_14, qs_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * os_9[k]
                      + qs_9[k];

            t_28[k] = ab_y[k] * os_9[k]
                      + qs_13[k];

            t_29[k] = ab_z[k] * os_9[k]
                      + qs_14[k];

            t_30[k] = ab_x[k] * os_10[k]
                      + qs_10[k];

            t_31[k] = ab_y[k] * os_10[k]
                      + qs_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, os_10, os_11, \
                         os_12, qs_11, qs_12, qs_16, qs_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * os_10[k]
                      + qs_16[k];

            t_33[k] = ab_x[k] * os_11[k]
                      + qs_11[k];

            t_34[k] = ab_y[k] * os_11[k]
                      + qs_16[k];

            t_35[k] = ab_z[k] * os_11[k]
                      + qs_17[k];

            t_36[k] = ab_x[k] * os_12[k]
                      + qs_12[k];

            t_37[k] = ab_y[k] * os_12[k]
                      + qs_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, os_12, os_13, \
                         os_14, qs_13, qs_14, qs_18, qs_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * os_12[k]
                      + qs_18[k];

            t_39[k] = ab_x[k] * os_13[k]
                      + qs_13[k];

            t_40[k] = ab_y[k] * os_13[k]
                      + qs_18[k];

            t_41[k] = ab_z[k] * os_13[k]
                      + qs_19[k];

            t_42[k] = ab_x[k] * os_14[k]
                      + qs_14[k];

            t_43[k] = ab_y[k] * os_14[k]
                      + qs_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, os_14, os_15, os_16, \
                         qs_15, qs_16, qs_20, qs_21, qs_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * os_14[k]
                      + qs_20[k];

            t_45[k] = ab_x[k] * os_15[k]
                      + qs_15[k];

            t_46[k] = ab_y[k] * os_15[k]
                      + qs_21[k];

            t_47[k] = ab_z[k] * os_15[k]
                      + qs_22[k];

            t_48[k] = ab_x[k] * os_16[k]
                      + qs_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, os_16, os_17, qs_17, \
                         qs_22, qs_23, qs_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * os_16[k]
                      + qs_22[k];

            t_50[k] = ab_z[k] * os_16[k]
                      + qs_23[k];

            t_51[k] = ab_x[k] * os_17[k]
                      + qs_17[k];

            t_52[k] = ab_y[k] * os_17[k]
                      + qs_23[k];

            t_53[k] = ab_z[k] * os_17[k]
                      + qs_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, os_18, os_19, \
                         qs_18, qs_19, qs_24, qs_25, qs_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * os_18[k]
                      + qs_18[k];

            t_55[k] = ab_y[k] * os_18[k]
                      + qs_24[k];

            t_56[k] = ab_z[k] * os_18[k]
                      + qs_25[k];

            t_57[k] = ab_x[k] * os_19[k]
                      + qs_19[k];

            t_58[k] = ab_y[k] * os_19[k]
                      + qs_25[k];

            t_59[k] = ab_z[k] * os_19[k]
                      + qs_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, os_20, os_21, qs_20, \
                         qs_21, qs_26, qs_27, qs_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * os_20[k]
                      + qs_20[k];

            t_61[k] = ab_y[k] * os_20[k]
                      + qs_26[k];

            t_62[k] = ab_z[k] * os_20[k]
                      + qs_27[k];

            t_63[k] = ab_x[k] * os_21[k]
                      + qs_21[k];

            t_64[k] = ab_y[k] * os_21[k]
                      + qs_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, os_21, os_22, \
                         os_23, qs_22, qs_23, qs_29, qs_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * os_21[k]
                      + qs_29[k];

            t_66[k] = ab_x[k] * os_22[k]
                      + qs_22[k];

            t_67[k] = ab_y[k] * os_22[k]
                      + qs_29[k];

            t_68[k] = ab_z[k] * os_22[k]
                      + qs_30[k];

            t_69[k] = ab_x[k] * os_23[k]
                      + qs_23[k];

            t_70[k] = ab_y[k] * os_23[k]
                      + qs_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, os_23, os_24, \
                         os_25, qs_24, qs_25, qs_31, qs_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * os_23[k]
                      + qs_31[k];

            t_72[k] = ab_x[k] * os_24[k]
                      + qs_24[k];

            t_73[k] = ab_y[k] * os_24[k]
                      + qs_31[k];

            t_74[k] = ab_z[k] * os_24[k]
                      + qs_32[k];

            t_75[k] = ab_x[k] * os_25[k]
                      + qs_25[k];

            t_76[k] = ab_y[k] * os_25[k]
                      + qs_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, os_25, os_26, \
                         os_27, qs_26, qs_27, qs_33, qs_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * os_25[k]
                      + qs_33[k];

            t_78[k] = ab_x[k] * os_26[k]
                      + qs_26[k];

            t_79[k] = ab_y[k] * os_26[k]
                      + qs_33[k];

            t_80[k] = ab_z[k] * os_26[k]
                      + qs_34[k];

            t_81[k] = ab_x[k] * os_27[k]
                      + qs_27[k];

            t_82[k] = ab_y[k] * os_27[k]
                      + qs_34[k];
        }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ab_z, os_27, os_28, os_29, \
                         qs_28, qs_29, qs_35, qs_36, qs_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * os_27[k]
                      + qs_35[k];

            t_84[k] = ab_x[k] * os_28[k]
                      + qs_28[k];

            t_85[k] = ab_y[k] * os_28[k]
                      + qs_36[k];

            t_86[k] = ab_z[k] * os_28[k]
                      + qs_37[k];

            t_87[k] = ab_x[k] * os_29[k]
                      + qs_29[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ab_x, ab_y, ab_z, os_29, os_30, qs_30, \
                         qs_37, qs_38, qs_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * os_29[k]
                      + qs_37[k];

            t_89[k] = ab_z[k] * os_29[k]
                      + qs_38[k];

            t_90[k] = ab_x[k] * os_30[k]
                      + qs_30[k];

            t_91[k] = ab_y[k] * os_30[k]
                      + qs_38[k];

            t_92[k] = ab_z[k] * os_30[k]
                      + qs_39[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, os_31, os_32, \
                         qs_31, qs_32, qs_39, qs_40, qs_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = ab_x[k] * os_31[k]
                      + qs_31[k];

            t_94[k] = ab_y[k] * os_31[k]
                      + qs_39[k];

            t_95[k] = ab_z[k] * os_31[k]
                      + qs_40[k];

            t_96[k] = ab_x[k] * os_32[k]
                      + qs_32[k];

            t_97[k] = ab_y[k] * os_32[k]
                      + qs_40[k];

            t_98[k] = ab_z[k] * os_32[k]
                      + qs_41[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, ab_x, ab_y, ab_z, os_33, \
                         os_34, qs_33, qs_34, qs_41, qs_42, qs_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * os_33[k]
                      + qs_33[k];

            t_100[k] = ab_y[k] * os_33[k]
                       + qs_41[k];

            t_101[k] = ab_z[k] * os_33[k]
                       + qs_42[k];

            t_102[k] = ab_x[k] * os_34[k]
                       + qs_34[k];

            t_103[k] = ab_y[k] * os_34[k]
                       + qs_42[k];

            t_104[k] = ab_z[k] * os_34[k]
                       + qs_43[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, os_35, os_36, \
                         qs_35, qs_36, qs_43, qs_44, qs_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * os_35[k]
                       + qs_35[k];

            t_106[k] = ab_y[k] * os_35[k]
                       + qs_43[k];

            t_107[k] = ab_z[k] * os_35[k]
                       + qs_44[k];

            t_108[k] = ab_x[k] * os_36[k]
                       + qs_36[k];

            t_109[k] = ab_y[k] * os_36[k]
                       + qs_45[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, ab_x, ab_y, ab_z, os_36, \
                         os_37, os_38, qs_37, qs_38, qs_46, qs_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_z[k] * os_36[k]
                       + qs_46[k];

            t_111[k] = ab_x[k] * os_37[k]
                       + qs_37[k];

            t_112[k] = ab_y[k] * os_37[k]
                       + qs_46[k];

            t_113[k] = ab_z[k] * os_37[k]
                       + qs_47[k];

            t_114[k] = ab_x[k] * os_38[k]
                       + qs_38[k];

            t_115[k] = ab_y[k] * os_38[k]
                       + qs_47[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, ab_x, ab_y, ab_z, os_38, \
                         os_39, os_40, qs_39, qs_40, qs_48, qs_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_z[k] * os_38[k]
                       + qs_48[k];

            t_117[k] = ab_x[k] * os_39[k]
                       + qs_39[k];

            t_118[k] = ab_y[k] * os_39[k]
                       + qs_48[k];

            t_119[k] = ab_z[k] * os_39[k]
                       + qs_49[k];

            t_120[k] = ab_x[k] * os_40[k]
                       + qs_40[k];

            t_121[k] = ab_y[k] * os_40[k]
                       + qs_49[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, ab_x, ab_y, ab_z, os_40, \
                         os_41, os_42, qs_41, qs_42, qs_50, qs_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_z[k] * os_40[k]
                       + qs_50[k];

            t_123[k] = ab_x[k] * os_41[k]
                       + qs_41[k];

            t_124[k] = ab_y[k] * os_41[k]
                       + qs_50[k];

            t_125[k] = ab_z[k] * os_41[k]
                       + qs_51[k];

            t_126[k] = ab_x[k] * os_42[k]
                       + qs_42[k];

            t_127[k] = ab_y[k] * os_42[k]
                       + qs_51[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, ab_x, ab_y, ab_z, os_42, \
                         os_43, os_44, qs_43, qs_44, qs_52, qs_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = ab_z[k] * os_42[k]
                       + qs_52[k];

            t_129[k] = ab_x[k] * os_43[k]
                       + qs_43[k];

            t_130[k] = ab_y[k] * os_43[k]
                       + qs_52[k];

            t_131[k] = ab_z[k] * os_43[k]
                       + qs_53[k];

            t_132[k] = ab_x[k] * os_44[k]
                       + qs_44[k];

            t_133[k] = ab_y[k] * os_44[k]
                       + qs_53[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, os_44, os_45, \
                         os_46, qs_45, qs_46, qs_54, qs_55, qs_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = ab_z[k] * os_44[k]
                       + qs_54[k];

            t_135[k] = ab_x[k] * os_45[k]
                       + qs_45[k];

            t_136[k] = ab_y[k] * os_45[k]
                       + qs_55[k];

            t_137[k] = ab_z[k] * os_45[k]
                       + qs_56[k];

            t_138[k] = ab_x[k] * os_46[k]
                       + qs_46[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, os_46, os_47, \
                         qs_47, qs_56, qs_57, qs_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_y[k] * os_46[k]
                       + qs_56[k];

            t_140[k] = ab_z[k] * os_46[k]
                       + qs_57[k];

            t_141[k] = ab_x[k] * os_47[k]
                       + qs_47[k];

            t_142[k] = ab_y[k] * os_47[k]
                       + qs_57[k];

            t_143[k] = ab_z[k] * os_47[k]
                       + qs_58[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, os_48, \
                         os_49, qs_48, qs_49, qs_58, qs_59, qs_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * os_48[k]
                       + qs_48[k];

            t_145[k] = ab_y[k] * os_48[k]
                       + qs_58[k];

            t_146[k] = ab_z[k] * os_48[k]
                       + qs_59[k];

            t_147[k] = ab_x[k] * os_49[k]
                       + qs_49[k];

            t_148[k] = ab_y[k] * os_49[k]
                       + qs_59[k];

            t_149[k] = ab_z[k] * os_49[k]
                       + qs_60[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, ab_x, ab_y, ab_z, os_50, \
                         os_51, qs_50, qs_51, qs_60, qs_61, qs_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * os_50[k]
                       + qs_50[k];

            t_151[k] = ab_y[k] * os_50[k]
                       + qs_60[k];

            t_152[k] = ab_z[k] * os_50[k]
                       + qs_61[k];

            t_153[k] = ab_x[k] * os_51[k]
                       + qs_51[k];

            t_154[k] = ab_y[k] * os_51[k]
                       + qs_61[k];

            t_155[k] = ab_z[k] * os_51[k]
                       + qs_62[k];
        }
    }
}

static auto
compute_hrr_op_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t os, const size_t qs,
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
        auto *t_198 = buffer.data(target + 198 * ncomps + c);
        auto *t_199 = buffer.data(target + 199 * ncomps + c);
        auto *t_200 = buffer.data(target + 200 * ncomps + c);
        auto *t_201 = buffer.data(target + 201 * ncomps + c);
        auto *t_202 = buffer.data(target + 202 * ncomps + c);
        auto *t_203 = buffer.data(target + 203 * ncomps + c);
        auto *t_204 = buffer.data(target + 204 * ncomps + c);
        auto *t_205 = buffer.data(target + 205 * ncomps + c);
        auto *t_206 = buffer.data(target + 206 * ncomps + c);
        auto *t_207 = buffer.data(target + 207 * ncomps + c);
        auto *t_208 = buffer.data(target + 208 * ncomps + c);
        auto *t_209 = buffer.data(target + 209 * ncomps + c);
        auto *t_210 = buffer.data(target + 210 * ncomps + c);
        auto *t_211 = buffer.data(target + 211 * ncomps + c);
        auto *t_212 = buffer.data(target + 212 * ncomps + c);
        auto *t_213 = buffer.data(target + 213 * ncomps + c);
        auto *t_214 = buffer.data(target + 214 * ncomps + c);
        auto *t_215 = buffer.data(target + 215 * ncomps + c);
        auto *t_216 = buffer.data(target + 216 * ncomps + c);
        auto *t_217 = buffer.data(target + 217 * ncomps + c);
        auto *t_218 = buffer.data(target + 218 * ncomps + c);
        auto *t_219 = buffer.data(target + 219 * ncomps + c);
        auto *t_220 = buffer.data(target + 220 * ncomps + c);
        auto *t_221 = buffer.data(target + 221 * ncomps + c);
        auto *t_222 = buffer.data(target + 222 * ncomps + c);
        auto *t_223 = buffer.data(target + 223 * ncomps + c);
        auto *t_224 = buffer.data(target + 224 * ncomps + c);
        auto *t_225 = buffer.data(target + 225 * ncomps + c);
        auto *t_226 = buffer.data(target + 226 * ncomps + c);
        auto *t_227 = buffer.data(target + 227 * ncomps + c);
        auto *t_228 = buffer.data(target + 228 * ncomps + c);
        auto *t_229 = buffer.data(target + 229 * ncomps + c);
        auto *t_230 = buffer.data(target + 230 * ncomps + c);
        auto *t_231 = buffer.data(target + 231 * ncomps + c);
        auto *t_232 = buffer.data(target + 232 * ncomps + c);
        auto *t_233 = buffer.data(target + 233 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

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

        const auto *qs_52 = buffer.data(qs + 52 * ncomps + c);
        const auto *qs_53 = buffer.data(qs + 53 * ncomps + c);
        const auto *qs_54 = buffer.data(qs + 54 * ncomps + c);
        const auto *qs_55 = buffer.data(qs + 55 * ncomps + c);
        const auto *qs_56 = buffer.data(qs + 56 * ncomps + c);
        const auto *qs_57 = buffer.data(qs + 57 * ncomps + c);
        const auto *qs_58 = buffer.data(qs + 58 * ncomps + c);
        const auto *qs_59 = buffer.data(qs + 59 * ncomps + c);
        const auto *qs_60 = buffer.data(qs + 60 * ncomps + c);
        const auto *qs_61 = buffer.data(qs + 61 * ncomps + c);
        const auto *qs_62 = buffer.data(qs + 62 * ncomps + c);
        const auto *qs_63 = buffer.data(qs + 63 * ncomps + c);
        const auto *qs_64 = buffer.data(qs + 64 * ncomps + c);
        const auto *qs_65 = buffer.data(qs + 65 * ncomps + c);
        const auto *qs_66 = buffer.data(qs + 66 * ncomps + c);
        const auto *qs_67 = buffer.data(qs + 67 * ncomps + c);
        const auto *qs_68 = buffer.data(qs + 68 * ncomps + c);
        const auto *qs_69 = buffer.data(qs + 69 * ncomps + c);
        const auto *qs_70 = buffer.data(qs + 70 * ncomps + c);
        const auto *qs_71 = buffer.data(qs + 71 * ncomps + c);
        const auto *qs_72 = buffer.data(qs + 72 * ncomps + c);
        const auto *qs_73 = buffer.data(qs + 73 * ncomps + c);
        const auto *qs_74 = buffer.data(qs + 74 * ncomps + c);
        const auto *qs_75 = buffer.data(qs + 75 * ncomps + c);
        const auto *qs_76 = buffer.data(qs + 76 * ncomps + c);
        const auto *qs_77 = buffer.data(qs + 77 * ncomps + c);
        const auto *qs_78 = buffer.data(qs + 78 * ncomps + c);
        const auto *qs_79 = buffer.data(qs + 79 * ncomps + c);
        const auto *qs_80 = buffer.data(qs + 80 * ncomps + c);
        const auto *qs_81 = buffer.data(qs + 81 * ncomps + c);
        const auto *qs_82 = buffer.data(qs + 82 * ncomps + c);
        const auto *qs_83 = buffer.data(qs + 83 * ncomps + c);
        const auto *qs_84 = buffer.data(qs + 84 * ncomps + c);
        const auto *qs_85 = buffer.data(qs + 85 * ncomps + c);
        const auto *qs_86 = buffer.data(qs + 86 * ncomps + c);
        const auto *qs_87 = buffer.data(qs + 87 * ncomps + c);
        const auto *qs_88 = buffer.data(qs + 88 * ncomps + c);
        const auto *qs_89 = buffer.data(qs + 89 * ncomps + c);
        const auto *qs_90 = buffer.data(qs + 90 * ncomps + c);

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, os_52, \
                         os_53, qs_52, qs_53, qs_62, qs_63, qs_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = ab_x[k] * os_52[k]
                       + qs_52[k];

            t_157[k] = ab_y[k] * os_52[k]
                       + qs_62[k];

            t_158[k] = ab_z[k] * os_52[k]
                       + qs_63[k];

            t_159[k] = ab_x[k] * os_53[k]
                       + qs_53[k];

            t_160[k] = ab_y[k] * os_53[k]
                       + qs_63[k];

            t_161[k] = ab_z[k] * os_53[k]
                       + qs_64[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, ab_z, os_54, os_55, \
                         qs_54, qs_55, qs_64, qs_65, qs_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * os_54[k]
                       + qs_54[k];

            t_163[k] = ab_y[k] * os_54[k]
                       + qs_64[k];

            t_164[k] = ab_z[k] * os_54[k]
                       + qs_65[k];

            t_165[k] = ab_x[k] * os_55[k]
                       + qs_55[k];

            t_166[k] = ab_y[k] * os_55[k]
                       + qs_66[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, t_172, ab_x, ab_y, ab_z, os_55, \
                         os_56, os_57, qs_56, qs_57, qs_67, qs_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = ab_z[k] * os_55[k]
                       + qs_67[k];

            t_168[k] = ab_x[k] * os_56[k]
                       + qs_56[k];

            t_169[k] = ab_y[k] * os_56[k]
                       + qs_67[k];

            t_170[k] = ab_z[k] * os_56[k]
                       + qs_68[k];

            t_171[k] = ab_x[k] * os_57[k]
                       + qs_57[k];

            t_172[k] = ab_y[k] * os_57[k]
                       + qs_68[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, ab_x, ab_y, ab_z, os_57, \
                         os_58, os_59, qs_58, qs_59, qs_69, qs_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_z[k] * os_57[k]
                       + qs_69[k];

            t_174[k] = ab_x[k] * os_58[k]
                       + qs_58[k];

            t_175[k] = ab_y[k] * os_58[k]
                       + qs_69[k];

            t_176[k] = ab_z[k] * os_58[k]
                       + qs_70[k];

            t_177[k] = ab_x[k] * os_59[k]
                       + qs_59[k];

            t_178[k] = ab_y[k] * os_59[k]
                       + qs_70[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, ab_z, os_59, \
                         os_60, os_61, qs_60, qs_61, qs_71, qs_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = ab_z[k] * os_59[k]
                       + qs_71[k];

            t_180[k] = ab_x[k] * os_60[k]
                       + qs_60[k];

            t_181[k] = ab_y[k] * os_60[k]
                       + qs_71[k];

            t_182[k] = ab_z[k] * os_60[k]
                       + qs_72[k];

            t_183[k] = ab_x[k] * os_61[k]
                       + qs_61[k];

            t_184[k] = ab_y[k] * os_61[k]
                       + qs_72[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, ab_x, ab_y, ab_z, os_61, \
                         os_62, os_63, qs_62, qs_63, qs_73, qs_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_z[k] * os_61[k]
                       + qs_73[k];

            t_186[k] = ab_x[k] * os_62[k]
                       + qs_62[k];

            t_187[k] = ab_y[k] * os_62[k]
                       + qs_73[k];

            t_188[k] = ab_z[k] * os_62[k]
                       + qs_74[k];

            t_189[k] = ab_x[k] * os_63[k]
                       + qs_63[k];

            t_190[k] = ab_y[k] * os_63[k]
                       + qs_74[k];
        }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, ab_x, ab_y, ab_z, os_63, \
                         os_64, os_65, qs_64, qs_65, qs_75, qs_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_191[k] = ab_z[k] * os_63[k]
                       + qs_75[k];

            t_192[k] = ab_x[k] * os_64[k]
                       + qs_64[k];

            t_193[k] = ab_y[k] * os_64[k]
                       + qs_75[k];

            t_194[k] = ab_z[k] * os_64[k]
                       + qs_76[k];

            t_195[k] = ab_x[k] * os_65[k]
                       + qs_65[k];

            t_196[k] = ab_y[k] * os_65[k]
                       + qs_76[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_x, ab_y, ab_z, os_65, os_66, \
                         os_67, qs_66, qs_67, qs_77, qs_78, qs_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_z[k] * os_65[k]
                       + qs_77[k];

            t_198[k] = ab_x[k] * os_66[k]
                       + qs_66[k];

            t_199[k] = ab_y[k] * os_66[k]
                       + qs_78[k];

            t_200[k] = ab_z[k] * os_66[k]
                       + qs_79[k];

            t_201[k] = ab_x[k] * os_67[k]
                       + qs_67[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, ab_x, ab_y, ab_z, os_67, os_68, \
                         qs_68, qs_79, qs_80, qs_81 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_y[k] * os_67[k]
                       + qs_79[k];

            t_203[k] = ab_z[k] * os_67[k]
                       + qs_80[k];

            t_204[k] = ab_x[k] * os_68[k]
                       + qs_68[k];

            t_205[k] = ab_y[k] * os_68[k]
                       + qs_80[k];

            t_206[k] = ab_z[k] * os_68[k]
                       + qs_81[k];
        }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, ab_x, ab_y, ab_z, os_69, \
                         os_70, qs_69, qs_70, qs_81, qs_82, qs_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_207[k] = ab_x[k] * os_69[k]
                       + qs_69[k];

            t_208[k] = ab_y[k] * os_69[k]
                       + qs_81[k];

            t_209[k] = ab_z[k] * os_69[k]
                       + qs_82[k];

            t_210[k] = ab_x[k] * os_70[k]
                       + qs_70[k];

            t_211[k] = ab_y[k] * os_70[k]
                       + qs_82[k];

            t_212[k] = ab_z[k] * os_70[k]
                       + qs_83[k];
        }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, ab_x, ab_y, ab_z, os_71, \
                         os_72, qs_71, qs_72, qs_83, qs_84, qs_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_213[k] = ab_x[k] * os_71[k]
                       + qs_71[k];

            t_214[k] = ab_y[k] * os_71[k]
                       + qs_83[k];

            t_215[k] = ab_z[k] * os_71[k]
                       + qs_84[k];

            t_216[k] = ab_x[k] * os_72[k]
                       + qs_72[k];

            t_217[k] = ab_y[k] * os_72[k]
                       + qs_84[k];

            t_218[k] = ab_z[k] * os_72[k]
                       + qs_85[k];
        }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, ab_x, ab_y, ab_z, os_73, \
                         os_74, qs_73, qs_74, qs_85, qs_86, qs_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_219[k] = ab_x[k] * os_73[k]
                       + qs_73[k];

            t_220[k] = ab_y[k] * os_73[k]
                       + qs_85[k];

            t_221[k] = ab_z[k] * os_73[k]
                       + qs_86[k];

            t_222[k] = ab_x[k] * os_74[k]
                       + qs_74[k];

            t_223[k] = ab_y[k] * os_74[k]
                       + qs_86[k];

            t_224[k] = ab_z[k] * os_74[k]
                       + qs_87[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, t_230, ab_x, ab_y, ab_z, os_75, \
                         os_76, qs_75, qs_76, qs_87, qs_88, qs_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * os_75[k]
                       + qs_75[k];

            t_226[k] = ab_y[k] * os_75[k]
                       + qs_87[k];

            t_227[k] = ab_z[k] * os_75[k]
                       + qs_88[k];

            t_228[k] = ab_x[k] * os_76[k]
                       + qs_76[k];

            t_229[k] = ab_y[k] * os_76[k]
                       + qs_88[k];

            t_230[k] = ab_z[k] * os_76[k]
                       + qs_89[k];
        }

#pragma omp simd aligned(t_231, t_232, t_233, ab_x, ab_y, ab_z, os_77, qs_77, qs_89, \
                         qs_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_231[k] = ab_x[k] * os_77[k]
                       + qs_77[k];

            t_232[k] = ab_y[k] * os_77[k]
                       + qs_89[k];

            t_233[k] = ab_z[k] * os_77[k]
                       + qs_90[k];
        }
    }
}

auto
compute_hrr_op_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t os, const size_t qs,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_op_out_of_first_piece0(buffer, coordinates, target, os, qs, ncomps, nmax);

    compute_hrr_op_out_of_first_piece1(buffer, coordinates, target, os, qs, ncomps, nmax);
}

}  // namespace simdtrf
