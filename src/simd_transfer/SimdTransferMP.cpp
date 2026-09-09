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


#include "SimdTransferMP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_mp_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ms, const size_t ns,
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

        const auto *ms_0 = buffer.data(ms + 0 * ncomps + c);
        const auto *ms_1 = buffer.data(ms + 1 * ncomps + c);
        const auto *ms_2 = buffer.data(ms + 2 * ncomps + c);
        const auto *ms_3 = buffer.data(ms + 3 * ncomps + c);
        const auto *ms_4 = buffer.data(ms + 4 * ncomps + c);
        const auto *ms_5 = buffer.data(ms + 5 * ncomps + c);
        const auto *ms_6 = buffer.data(ms + 6 * ncomps + c);
        const auto *ms_7 = buffer.data(ms + 7 * ncomps + c);
        const auto *ms_8 = buffer.data(ms + 8 * ncomps + c);
        const auto *ms_9 = buffer.data(ms + 9 * ncomps + c);
        const auto *ms_10 = buffer.data(ms + 10 * ncomps + c);
        const auto *ms_11 = buffer.data(ms + 11 * ncomps + c);
        const auto *ms_12 = buffer.data(ms + 12 * ncomps + c);
        const auto *ms_13 = buffer.data(ms + 13 * ncomps + c);
        const auto *ms_14 = buffer.data(ms + 14 * ncomps + c);
        const auto *ms_15 = buffer.data(ms + 15 * ncomps + c);
        const auto *ms_16 = buffer.data(ms + 16 * ncomps + c);
        const auto *ms_17 = buffer.data(ms + 17 * ncomps + c);
        const auto *ms_18 = buffer.data(ms + 18 * ncomps + c);
        const auto *ms_19 = buffer.data(ms + 19 * ncomps + c);
        const auto *ms_20 = buffer.data(ms + 20 * ncomps + c);
        const auto *ms_21 = buffer.data(ms + 21 * ncomps + c);
        const auto *ms_22 = buffer.data(ms + 22 * ncomps + c);
        const auto *ms_23 = buffer.data(ms + 23 * ncomps + c);
        const auto *ms_24 = buffer.data(ms + 24 * ncomps + c);
        const auto *ms_25 = buffer.data(ms + 25 * ncomps + c);
        const auto *ms_26 = buffer.data(ms + 26 * ncomps + c);
        const auto *ms_27 = buffer.data(ms + 27 * ncomps + c);
        const auto *ms_28 = buffer.data(ms + 28 * ncomps + c);
        const auto *ms_29 = buffer.data(ms + 29 * ncomps + c);
        const auto *ms_30 = buffer.data(ms + 30 * ncomps + c);
        const auto *ms_31 = buffer.data(ms + 31 * ncomps + c);
        const auto *ms_32 = buffer.data(ms + 32 * ncomps + c);
        const auto *ms_33 = buffer.data(ms + 33 * ncomps + c);
        const auto *ms_34 = buffer.data(ms + 34 * ncomps + c);
        const auto *ms_35 = buffer.data(ms + 35 * ncomps + c);
        const auto *ms_36 = buffer.data(ms + 36 * ncomps + c);
        const auto *ms_37 = buffer.data(ms + 37 * ncomps + c);
        const auto *ms_38 = buffer.data(ms + 38 * ncomps + c);
        const auto *ms_39 = buffer.data(ms + 39 * ncomps + c);
        const auto *ms_40 = buffer.data(ms + 40 * ncomps + c);
        const auto *ms_41 = buffer.data(ms + 41 * ncomps + c);
        const auto *ms_42 = buffer.data(ms + 42 * ncomps + c);
        const auto *ms_43 = buffer.data(ms + 43 * ncomps + c);
        const auto *ms_44 = buffer.data(ms + 44 * ncomps + c);
        const auto *ms_45 = buffer.data(ms + 45 * ncomps + c);
        const auto *ms_46 = buffer.data(ms + 46 * ncomps + c);
        const auto *ms_47 = buffer.data(ms + 47 * ncomps + c);
        const auto *ms_48 = buffer.data(ms + 48 * ncomps + c);
        const auto *ms_49 = buffer.data(ms + 49 * ncomps + c);
        const auto *ms_50 = buffer.data(ms + 50 * ncomps + c);
        const auto *ms_51 = buffer.data(ms + 51 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, ms_0, ms_1, ns_0, \
                         ns_1, ns_2, ns_3, ns_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ms_0[k]
                     + ns_0[k];

            t_1[k] = ab_y[k] * ms_0[k]
                     + ns_1[k];

            t_2[k] = ab_z[k] * ms_0[k]
                     + ns_2[k];

            t_3[k] = ab_x[k] * ms_1[k]
                     + ns_1[k];

            t_4[k] = ab_y[k] * ms_1[k]
                     + ns_3[k];

            t_5[k] = ab_z[k] * ms_1[k]
                     + ns_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, ms_2, ms_3, ns_2, ns_3, \
                         ns_4, ns_5, ns_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * ms_2[k]
                     + ns_2[k];

            t_7[k] = ab_y[k] * ms_2[k]
                     + ns_4[k];

            t_8[k] = ab_z[k] * ms_2[k]
                     + ns_5[k];

            t_9[k] = ab_x[k] * ms_3[k]
                     + ns_3[k];

            t_10[k] = ab_y[k] * ms_3[k]
                      + ns_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, ms_3, ms_4, \
                         ms_5, ns_4, ns_5, ns_7, ns_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * ms_3[k]
                      + ns_7[k];

            t_12[k] = ab_x[k] * ms_4[k]
                      + ns_4[k];

            t_13[k] = ab_y[k] * ms_4[k]
                      + ns_7[k];

            t_14[k] = ab_z[k] * ms_4[k]
                      + ns_8[k];

            t_15[k] = ab_x[k] * ms_5[k]
                      + ns_5[k];

            t_16[k] = ab_y[k] * ms_5[k]
                      + ns_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, ms_5, ms_6, ms_7, \
                         ns_6, ns_7, ns_9, ns_10, ns_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * ms_5[k]
                      + ns_9[k];

            t_18[k] = ab_x[k] * ms_6[k]
                      + ns_6[k];

            t_19[k] = ab_y[k] * ms_6[k]
                      + ns_10[k];

            t_20[k] = ab_z[k] * ms_6[k]
                      + ns_11[k];

            t_21[k] = ab_x[k] * ms_7[k]
                      + ns_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, ms_7, ms_8, ns_8, \
                         ns_11, ns_12, ns_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * ms_7[k]
                      + ns_11[k];

            t_23[k] = ab_z[k] * ms_7[k]
                      + ns_12[k];

            t_24[k] = ab_x[k] * ms_8[k]
                      + ns_8[k];

            t_25[k] = ab_y[k] * ms_8[k]
                      + ns_12[k];

            t_26[k] = ab_z[k] * ms_8[k]
                      + ns_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ms_9, ms_10, ns_9, \
                         ns_10, ns_13, ns_14, ns_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * ms_9[k]
                      + ns_9[k];

            t_28[k] = ab_y[k] * ms_9[k]
                      + ns_13[k];

            t_29[k] = ab_z[k] * ms_9[k]
                      + ns_14[k];

            t_30[k] = ab_x[k] * ms_10[k]
                      + ns_10[k];

            t_31[k] = ab_y[k] * ms_10[k]
                      + ns_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, ms_10, ms_11, \
                         ms_12, ns_11, ns_12, ns_16, ns_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * ms_10[k]
                      + ns_16[k];

            t_33[k] = ab_x[k] * ms_11[k]
                      + ns_11[k];

            t_34[k] = ab_y[k] * ms_11[k]
                      + ns_16[k];

            t_35[k] = ab_z[k] * ms_11[k]
                      + ns_17[k];

            t_36[k] = ab_x[k] * ms_12[k]
                      + ns_12[k];

            t_37[k] = ab_y[k] * ms_12[k]
                      + ns_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, ms_12, ms_13, \
                         ms_14, ns_13, ns_14, ns_18, ns_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * ms_12[k]
                      + ns_18[k];

            t_39[k] = ab_x[k] * ms_13[k]
                      + ns_13[k];

            t_40[k] = ab_y[k] * ms_13[k]
                      + ns_18[k];

            t_41[k] = ab_z[k] * ms_13[k]
                      + ns_19[k];

            t_42[k] = ab_x[k] * ms_14[k]
                      + ns_14[k];

            t_43[k] = ab_y[k] * ms_14[k]
                      + ns_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, ms_14, ms_15, ms_16, \
                         ns_15, ns_16, ns_20, ns_21, ns_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * ms_14[k]
                      + ns_20[k];

            t_45[k] = ab_x[k] * ms_15[k]
                      + ns_15[k];

            t_46[k] = ab_y[k] * ms_15[k]
                      + ns_21[k];

            t_47[k] = ab_z[k] * ms_15[k]
                      + ns_22[k];

            t_48[k] = ab_x[k] * ms_16[k]
                      + ns_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ms_16, ms_17, ns_17, \
                         ns_22, ns_23, ns_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * ms_16[k]
                      + ns_22[k];

            t_50[k] = ab_z[k] * ms_16[k]
                      + ns_23[k];

            t_51[k] = ab_x[k] * ms_17[k]
                      + ns_17[k];

            t_52[k] = ab_y[k] * ms_17[k]
                      + ns_23[k];

            t_53[k] = ab_z[k] * ms_17[k]
                      + ns_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ms_18, ms_19, \
                         ns_18, ns_19, ns_24, ns_25, ns_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * ms_18[k]
                      + ns_18[k];

            t_55[k] = ab_y[k] * ms_18[k]
                      + ns_24[k];

            t_56[k] = ab_z[k] * ms_18[k]
                      + ns_25[k];

            t_57[k] = ab_x[k] * ms_19[k]
                      + ns_19[k];

            t_58[k] = ab_y[k] * ms_19[k]
                      + ns_25[k];

            t_59[k] = ab_z[k] * ms_19[k]
                      + ns_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, ms_20, ms_21, ns_20, \
                         ns_21, ns_26, ns_27, ns_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ms_20[k]
                      + ns_20[k];

            t_61[k] = ab_y[k] * ms_20[k]
                      + ns_26[k];

            t_62[k] = ab_z[k] * ms_20[k]
                      + ns_27[k];

            t_63[k] = ab_x[k] * ms_21[k]
                      + ns_21[k];

            t_64[k] = ab_y[k] * ms_21[k]
                      + ns_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, ms_21, ms_22, \
                         ms_23, ns_22, ns_23, ns_29, ns_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * ms_21[k]
                      + ns_29[k];

            t_66[k] = ab_x[k] * ms_22[k]
                      + ns_22[k];

            t_67[k] = ab_y[k] * ms_22[k]
                      + ns_29[k];

            t_68[k] = ab_z[k] * ms_22[k]
                      + ns_30[k];

            t_69[k] = ab_x[k] * ms_23[k]
                      + ns_23[k];

            t_70[k] = ab_y[k] * ms_23[k]
                      + ns_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, ms_23, ms_24, \
                         ms_25, ns_24, ns_25, ns_31, ns_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * ms_23[k]
                      + ns_31[k];

            t_72[k] = ab_x[k] * ms_24[k]
                      + ns_24[k];

            t_73[k] = ab_y[k] * ms_24[k]
                      + ns_31[k];

            t_74[k] = ab_z[k] * ms_24[k]
                      + ns_32[k];

            t_75[k] = ab_x[k] * ms_25[k]
                      + ns_25[k];

            t_76[k] = ab_y[k] * ms_25[k]
                      + ns_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, ms_25, ms_26, \
                         ms_27, ns_26, ns_27, ns_33, ns_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * ms_25[k]
                      + ns_33[k];

            t_78[k] = ab_x[k] * ms_26[k]
                      + ns_26[k];

            t_79[k] = ab_y[k] * ms_26[k]
                      + ns_33[k];

            t_80[k] = ab_z[k] * ms_26[k]
                      + ns_34[k];

            t_81[k] = ab_x[k] * ms_27[k]
                      + ns_27[k];

            t_82[k] = ab_y[k] * ms_27[k]
                      + ns_34[k];
        }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ab_z, ms_27, ms_28, ms_29, \
                         ns_28, ns_29, ns_35, ns_36, ns_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * ms_27[k]
                      + ns_35[k];

            t_84[k] = ab_x[k] * ms_28[k]
                      + ns_28[k];

            t_85[k] = ab_y[k] * ms_28[k]
                      + ns_36[k];

            t_86[k] = ab_z[k] * ms_28[k]
                      + ns_37[k];

            t_87[k] = ab_x[k] * ms_29[k]
                      + ns_29[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ab_x, ab_y, ab_z, ms_29, ms_30, ns_30, \
                         ns_37, ns_38, ns_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * ms_29[k]
                      + ns_37[k];

            t_89[k] = ab_z[k] * ms_29[k]
                      + ns_38[k];

            t_90[k] = ab_x[k] * ms_30[k]
                      + ns_30[k];

            t_91[k] = ab_y[k] * ms_30[k]
                      + ns_38[k];

            t_92[k] = ab_z[k] * ms_30[k]
                      + ns_39[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, ms_31, ms_32, \
                         ns_31, ns_32, ns_39, ns_40, ns_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = ab_x[k] * ms_31[k]
                      + ns_31[k];

            t_94[k] = ab_y[k] * ms_31[k]
                      + ns_39[k];

            t_95[k] = ab_z[k] * ms_31[k]
                      + ns_40[k];

            t_96[k] = ab_x[k] * ms_32[k]
                      + ns_32[k];

            t_97[k] = ab_y[k] * ms_32[k]
                      + ns_40[k];

            t_98[k] = ab_z[k] * ms_32[k]
                      + ns_41[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, ab_x, ab_y, ab_z, ms_33, \
                         ms_34, ns_33, ns_34, ns_41, ns_42, ns_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * ms_33[k]
                      + ns_33[k];

            t_100[k] = ab_y[k] * ms_33[k]
                       + ns_41[k];

            t_101[k] = ab_z[k] * ms_33[k]
                       + ns_42[k];

            t_102[k] = ab_x[k] * ms_34[k]
                       + ns_34[k];

            t_103[k] = ab_y[k] * ms_34[k]
                       + ns_42[k];

            t_104[k] = ab_z[k] * ms_34[k]
                       + ns_43[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, ms_35, ms_36, \
                         ns_35, ns_36, ns_43, ns_44, ns_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * ms_35[k]
                       + ns_35[k];

            t_106[k] = ab_y[k] * ms_35[k]
                       + ns_43[k];

            t_107[k] = ab_z[k] * ms_35[k]
                       + ns_44[k];

            t_108[k] = ab_x[k] * ms_36[k]
                       + ns_36[k];

            t_109[k] = ab_y[k] * ms_36[k]
                       + ns_45[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, ab_x, ab_y, ab_z, ms_36, \
                         ms_37, ms_38, ns_37, ns_38, ns_46, ns_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_z[k] * ms_36[k]
                       + ns_46[k];

            t_111[k] = ab_x[k] * ms_37[k]
                       + ns_37[k];

            t_112[k] = ab_y[k] * ms_37[k]
                       + ns_46[k];

            t_113[k] = ab_z[k] * ms_37[k]
                       + ns_47[k];

            t_114[k] = ab_x[k] * ms_38[k]
                       + ns_38[k];

            t_115[k] = ab_y[k] * ms_38[k]
                       + ns_47[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, ab_x, ab_y, ab_z, ms_38, \
                         ms_39, ms_40, ns_39, ns_40, ns_48, ns_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_z[k] * ms_38[k]
                       + ns_48[k];

            t_117[k] = ab_x[k] * ms_39[k]
                       + ns_39[k];

            t_118[k] = ab_y[k] * ms_39[k]
                       + ns_48[k];

            t_119[k] = ab_z[k] * ms_39[k]
                       + ns_49[k];

            t_120[k] = ab_x[k] * ms_40[k]
                       + ns_40[k];

            t_121[k] = ab_y[k] * ms_40[k]
                       + ns_49[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, ab_x, ab_y, ab_z, ms_40, \
                         ms_41, ms_42, ns_41, ns_42, ns_50, ns_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_z[k] * ms_40[k]
                       + ns_50[k];

            t_123[k] = ab_x[k] * ms_41[k]
                       + ns_41[k];

            t_124[k] = ab_y[k] * ms_41[k]
                       + ns_50[k];

            t_125[k] = ab_z[k] * ms_41[k]
                       + ns_51[k];

            t_126[k] = ab_x[k] * ms_42[k]
                       + ns_42[k];

            t_127[k] = ab_y[k] * ms_42[k]
                       + ns_51[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, ab_x, ab_y, ab_z, ms_42, \
                         ms_43, ms_44, ns_43, ns_44, ns_52, ns_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = ab_z[k] * ms_42[k]
                       + ns_52[k];

            t_129[k] = ab_x[k] * ms_43[k]
                       + ns_43[k];

            t_130[k] = ab_y[k] * ms_43[k]
                       + ns_52[k];

            t_131[k] = ab_z[k] * ms_43[k]
                       + ns_53[k];

            t_132[k] = ab_x[k] * ms_44[k]
                       + ns_44[k];

            t_133[k] = ab_y[k] * ms_44[k]
                       + ns_53[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, ms_44, ms_45, \
                         ms_46, ns_45, ns_46, ns_54, ns_55, ns_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = ab_z[k] * ms_44[k]
                       + ns_54[k];

            t_135[k] = ab_x[k] * ms_45[k]
                       + ns_45[k];

            t_136[k] = ab_y[k] * ms_45[k]
                       + ns_55[k];

            t_137[k] = ab_z[k] * ms_45[k]
                       + ns_56[k];

            t_138[k] = ab_x[k] * ms_46[k]
                       + ns_46[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, ms_46, ms_47, \
                         ns_47, ns_56, ns_57, ns_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_y[k] * ms_46[k]
                       + ns_56[k];

            t_140[k] = ab_z[k] * ms_46[k]
                       + ns_57[k];

            t_141[k] = ab_x[k] * ms_47[k]
                       + ns_47[k];

            t_142[k] = ab_y[k] * ms_47[k]
                       + ns_57[k];

            t_143[k] = ab_z[k] * ms_47[k]
                       + ns_58[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, ms_48, \
                         ms_49, ns_48, ns_49, ns_58, ns_59, ns_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * ms_48[k]
                       + ns_48[k];

            t_145[k] = ab_y[k] * ms_48[k]
                       + ns_58[k];

            t_146[k] = ab_z[k] * ms_48[k]
                       + ns_59[k];

            t_147[k] = ab_x[k] * ms_49[k]
                       + ns_49[k];

            t_148[k] = ab_y[k] * ms_49[k]
                       + ns_59[k];

            t_149[k] = ab_z[k] * ms_49[k]
                       + ns_60[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, ab_x, ab_y, ab_z, ms_50, \
                         ms_51, ns_50, ns_51, ns_60, ns_61, ns_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ms_50[k]
                       + ns_50[k];

            t_151[k] = ab_y[k] * ms_50[k]
                       + ns_60[k];

            t_152[k] = ab_z[k] * ms_50[k]
                       + ns_61[k];

            t_153[k] = ab_x[k] * ms_51[k]
                       + ns_51[k];

            t_154[k] = ab_y[k] * ms_51[k]
                       + ns_61[k];

            t_155[k] = ab_z[k] * ms_51[k]
                       + ns_62[k];
        }
    }
}

static auto
compute_hrr_mp_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ms, const size_t ns,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ms_52 = buffer.data(ms + 52 * ncomps + c);
        const auto *ms_53 = buffer.data(ms + 53 * ncomps + c);
        const auto *ms_54 = buffer.data(ms + 54 * ncomps + c);

        const auto *ns_52 = buffer.data(ns + 52 * ncomps + c);
        const auto *ns_53 = buffer.data(ns + 53 * ncomps + c);
        const auto *ns_54 = buffer.data(ns + 54 * ncomps + c);
        const auto *ns_62 = buffer.data(ns + 62 * ncomps + c);
        const auto *ns_63 = buffer.data(ns + 63 * ncomps + c);
        const auto *ns_64 = buffer.data(ns + 64 * ncomps + c);
        const auto *ns_65 = buffer.data(ns + 65 * ncomps + c);

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, ms_52, \
                         ms_53, ns_52, ns_53, ns_62, ns_63, ns_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = ab_x[k] * ms_52[k]
                       + ns_52[k];

            t_157[k] = ab_y[k] * ms_52[k]
                       + ns_62[k];

            t_158[k] = ab_z[k] * ms_52[k]
                       + ns_63[k];

            t_159[k] = ab_x[k] * ms_53[k]
                       + ns_53[k];

            t_160[k] = ab_y[k] * ms_53[k]
                       + ns_63[k];

            t_161[k] = ab_z[k] * ms_53[k]
                       + ns_64[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, ab_x, ab_y, ab_z, ms_54, ns_54, ns_64, \
                         ns_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * ms_54[k]
                       + ns_54[k];

            t_163[k] = ab_y[k] * ms_54[k]
                       + ns_64[k];

            t_164[k] = ab_z[k] * ms_54[k]
                       + ns_65[k];
        }
    }
}

auto
compute_hrr_mp_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t ms, const size_t ns,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_mp_out_of_first_piece0(buffer, coordinates, target, ms, ns, ncomps, nmax);

    compute_hrr_mp_out_of_first_piece1(buffer, coordinates, target, ms, ns, ncomps, nmax);
}

static auto
compute_hrr_mp_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ms, const size_t ns, const size_t ncomps,
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

        const auto *ms_0 = buffer.data(ms + 0 * ncomps + c);
        const auto *ms_1 = buffer.data(ms + 1 * ncomps + c);
        const auto *ms_2 = buffer.data(ms + 2 * ncomps + c);
        const auto *ms_3 = buffer.data(ms + 3 * ncomps + c);
        const auto *ms_4 = buffer.data(ms + 4 * ncomps + c);
        const auto *ms_5 = buffer.data(ms + 5 * ncomps + c);
        const auto *ms_6 = buffer.data(ms + 6 * ncomps + c);
        const auto *ms_7 = buffer.data(ms + 7 * ncomps + c);
        const auto *ms_8 = buffer.data(ms + 8 * ncomps + c);
        const auto *ms_9 = buffer.data(ms + 9 * ncomps + c);
        const auto *ms_10 = buffer.data(ms + 10 * ncomps + c);
        const auto *ms_11 = buffer.data(ms + 11 * ncomps + c);
        const auto *ms_12 = buffer.data(ms + 12 * ncomps + c);
        const auto *ms_13 = buffer.data(ms + 13 * ncomps + c);
        const auto *ms_14 = buffer.data(ms + 14 * ncomps + c);
        const auto *ms_15 = buffer.data(ms + 15 * ncomps + c);
        const auto *ms_16 = buffer.data(ms + 16 * ncomps + c);
        const auto *ms_17 = buffer.data(ms + 17 * ncomps + c);
        const auto *ms_18 = buffer.data(ms + 18 * ncomps + c);
        const auto *ms_19 = buffer.data(ms + 19 * ncomps + c);
        const auto *ms_20 = buffer.data(ms + 20 * ncomps + c);
        const auto *ms_21 = buffer.data(ms + 21 * ncomps + c);
        const auto *ms_22 = buffer.data(ms + 22 * ncomps + c);
        const auto *ms_23 = buffer.data(ms + 23 * ncomps + c);
        const auto *ms_24 = buffer.data(ms + 24 * ncomps + c);
        const auto *ms_25 = buffer.data(ms + 25 * ncomps + c);
        const auto *ms_26 = buffer.data(ms + 26 * ncomps + c);
        const auto *ms_27 = buffer.data(ms + 27 * ncomps + c);
        const auto *ms_28 = buffer.data(ms + 28 * ncomps + c);
        const auto *ms_29 = buffer.data(ms + 29 * ncomps + c);
        const auto *ms_30 = buffer.data(ms + 30 * ncomps + c);
        const auto *ms_31 = buffer.data(ms + 31 * ncomps + c);
        const auto *ms_32 = buffer.data(ms + 32 * ncomps + c);
        const auto *ms_33 = buffer.data(ms + 33 * ncomps + c);
        const auto *ms_34 = buffer.data(ms + 34 * ncomps + c);
        const auto *ms_35 = buffer.data(ms + 35 * ncomps + c);
        const auto *ms_36 = buffer.data(ms + 36 * ncomps + c);
        const auto *ms_37 = buffer.data(ms + 37 * ncomps + c);
        const auto *ms_38 = buffer.data(ms + 38 * ncomps + c);
        const auto *ms_39 = buffer.data(ms + 39 * ncomps + c);
        const auto *ms_40 = buffer.data(ms + 40 * ncomps + c);
        const auto *ms_41 = buffer.data(ms + 41 * ncomps + c);
        const auto *ms_42 = buffer.data(ms + 42 * ncomps + c);
        const auto *ms_43 = buffer.data(ms + 43 * ncomps + c);
        const auto *ms_44 = buffer.data(ms + 44 * ncomps + c);
        const auto *ms_45 = buffer.data(ms + 45 * ncomps + c);
        const auto *ms_46 = buffer.data(ms + 46 * ncomps + c);
        const auto *ms_47 = buffer.data(ms + 47 * ncomps + c);
        const auto *ms_48 = buffer.data(ms + 48 * ncomps + c);
        const auto *ms_49 = buffer.data(ms + 49 * ncomps + c);
        const auto *ms_50 = buffer.data(ms + 50 * ncomps + c);
        const auto *ms_51 = buffer.data(ms + 51 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, ms_0, ms_1, ns_0, \
                         ns_1, ns_2, ns_3, ns_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ms_0[k]
                     + ns_0[k];

            t_1[k] = ab_y[k] * ms_0[k]
                     + ns_1[k];

            t_2[k] = ab_z[k] * ms_0[k]
                     + ns_2[k];

            t_3[k] = ab_x[k] * ms_1[k]
                     + ns_1[k];

            t_4[k] = ab_y[k] * ms_1[k]
                     + ns_3[k];

            t_5[k] = ab_z[k] * ms_1[k]
                     + ns_4[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, ms_2, ms_3, ns_2, ns_3, \
                         ns_4, ns_5, ns_6 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = ab_x[k] * ms_2[k]
                     + ns_2[k];

            t_7[k] = ab_y[k] * ms_2[k]
                     + ns_4[k];

            t_8[k] = ab_z[k] * ms_2[k]
                     + ns_5[k];

            t_9[k] = ab_x[k] * ms_3[k]
                     + ns_3[k];

            t_10[k] = ab_y[k] * ms_3[k]
                      + ns_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, ms_3, ms_4, \
                         ms_5, ns_4, ns_5, ns_7, ns_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_z[k] * ms_3[k]
                      + ns_7[k];

            t_12[k] = ab_x[k] * ms_4[k]
                      + ns_4[k];

            t_13[k] = ab_y[k] * ms_4[k]
                      + ns_7[k];

            t_14[k] = ab_z[k] * ms_4[k]
                      + ns_8[k];

            t_15[k] = ab_x[k] * ms_5[k]
                      + ns_5[k];

            t_16[k] = ab_y[k] * ms_5[k]
                      + ns_8[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, ms_5, ms_6, ms_7, \
                         ns_6, ns_7, ns_9, ns_10, ns_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_z[k] * ms_5[k]
                      + ns_9[k];

            t_18[k] = ab_x[k] * ms_6[k]
                      + ns_6[k];

            t_19[k] = ab_y[k] * ms_6[k]
                      + ns_10[k];

            t_20[k] = ab_z[k] * ms_6[k]
                      + ns_11[k];

            t_21[k] = ab_x[k] * ms_7[k]
                      + ns_7[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, ms_7, ms_8, ns_8, \
                         ns_11, ns_12, ns_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_y[k] * ms_7[k]
                      + ns_11[k];

            t_23[k] = ab_z[k] * ms_7[k]
                      + ns_12[k];

            t_24[k] = ab_x[k] * ms_8[k]
                      + ns_8[k];

            t_25[k] = ab_y[k] * ms_8[k]
                      + ns_12[k];

            t_26[k] = ab_z[k] * ms_8[k]
                      + ns_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ms_9, ms_10, ns_9, \
                         ns_10, ns_13, ns_14, ns_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * ms_9[k]
                      + ns_9[k];

            t_28[k] = ab_y[k] * ms_9[k]
                      + ns_13[k];

            t_29[k] = ab_z[k] * ms_9[k]
                      + ns_14[k];

            t_30[k] = ab_x[k] * ms_10[k]
                      + ns_10[k];

            t_31[k] = ab_y[k] * ms_10[k]
                      + ns_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, ms_10, ms_11, \
                         ms_12, ns_11, ns_12, ns_16, ns_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * ms_10[k]
                      + ns_16[k];

            t_33[k] = ab_x[k] * ms_11[k]
                      + ns_11[k];

            t_34[k] = ab_y[k] * ms_11[k]
                      + ns_16[k];

            t_35[k] = ab_z[k] * ms_11[k]
                      + ns_17[k];

            t_36[k] = ab_x[k] * ms_12[k]
                      + ns_12[k];

            t_37[k] = ab_y[k] * ms_12[k]
                      + ns_17[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, ms_12, ms_13, \
                         ms_14, ns_13, ns_14, ns_18, ns_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_z[k] * ms_12[k]
                      + ns_18[k];

            t_39[k] = ab_x[k] * ms_13[k]
                      + ns_13[k];

            t_40[k] = ab_y[k] * ms_13[k]
                      + ns_18[k];

            t_41[k] = ab_z[k] * ms_13[k]
                      + ns_19[k];

            t_42[k] = ab_x[k] * ms_14[k]
                      + ns_14[k];

            t_43[k] = ab_y[k] * ms_14[k]
                      + ns_19[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, ms_14, ms_15, ms_16, \
                         ns_15, ns_16, ns_20, ns_21, ns_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_z[k] * ms_14[k]
                      + ns_20[k];

            t_45[k] = ab_x[k] * ms_15[k]
                      + ns_15[k];

            t_46[k] = ab_y[k] * ms_15[k]
                      + ns_21[k];

            t_47[k] = ab_z[k] * ms_15[k]
                      + ns_22[k];

            t_48[k] = ab_x[k] * ms_16[k]
                      + ns_16[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ms_16, ms_17, ns_17, \
                         ns_22, ns_23, ns_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_y[k] * ms_16[k]
                      + ns_22[k];

            t_50[k] = ab_z[k] * ms_16[k]
                      + ns_23[k];

            t_51[k] = ab_x[k] * ms_17[k]
                      + ns_17[k];

            t_52[k] = ab_y[k] * ms_17[k]
                      + ns_23[k];

            t_53[k] = ab_z[k] * ms_17[k]
                      + ns_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ms_18, ms_19, \
                         ns_18, ns_19, ns_24, ns_25, ns_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * ms_18[k]
                      + ns_18[k];

            t_55[k] = ab_y[k] * ms_18[k]
                      + ns_24[k];

            t_56[k] = ab_z[k] * ms_18[k]
                      + ns_25[k];

            t_57[k] = ab_x[k] * ms_19[k]
                      + ns_19[k];

            t_58[k] = ab_y[k] * ms_19[k]
                      + ns_25[k];

            t_59[k] = ab_z[k] * ms_19[k]
                      + ns_26[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, ms_20, ms_21, ns_20, \
                         ns_21, ns_26, ns_27, ns_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ms_20[k]
                      + ns_20[k];

            t_61[k] = ab_y[k] * ms_20[k]
                      + ns_26[k];

            t_62[k] = ab_z[k] * ms_20[k]
                      + ns_27[k];

            t_63[k] = ab_x[k] * ms_21[k]
                      + ns_21[k];

            t_64[k] = ab_y[k] * ms_21[k]
                      + ns_28[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, ms_21, ms_22, \
                         ms_23, ns_22, ns_23, ns_29, ns_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_z[k] * ms_21[k]
                      + ns_29[k];

            t_66[k] = ab_x[k] * ms_22[k]
                      + ns_22[k];

            t_67[k] = ab_y[k] * ms_22[k]
                      + ns_29[k];

            t_68[k] = ab_z[k] * ms_22[k]
                      + ns_30[k];

            t_69[k] = ab_x[k] * ms_23[k]
                      + ns_23[k];

            t_70[k] = ab_y[k] * ms_23[k]
                      + ns_30[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, ms_23, ms_24, \
                         ms_25, ns_24, ns_25, ns_31, ns_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_z[k] * ms_23[k]
                      + ns_31[k];

            t_72[k] = ab_x[k] * ms_24[k]
                      + ns_24[k];

            t_73[k] = ab_y[k] * ms_24[k]
                      + ns_31[k];

            t_74[k] = ab_z[k] * ms_24[k]
                      + ns_32[k];

            t_75[k] = ab_x[k] * ms_25[k]
                      + ns_25[k];

            t_76[k] = ab_y[k] * ms_25[k]
                      + ns_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, ms_25, ms_26, \
                         ms_27, ns_26, ns_27, ns_33, ns_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * ms_25[k]
                      + ns_33[k];

            t_78[k] = ab_x[k] * ms_26[k]
                      + ns_26[k];

            t_79[k] = ab_y[k] * ms_26[k]
                      + ns_33[k];

            t_80[k] = ab_z[k] * ms_26[k]
                      + ns_34[k];

            t_81[k] = ab_x[k] * ms_27[k]
                      + ns_27[k];

            t_82[k] = ab_y[k] * ms_27[k]
                      + ns_34[k];
        }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ab_z, ms_27, ms_28, ms_29, \
                         ns_28, ns_29, ns_35, ns_36, ns_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_z[k] * ms_27[k]
                      + ns_35[k];

            t_84[k] = ab_x[k] * ms_28[k]
                      + ns_28[k];

            t_85[k] = ab_y[k] * ms_28[k]
                      + ns_36[k];

            t_86[k] = ab_z[k] * ms_28[k]
                      + ns_37[k];

            t_87[k] = ab_x[k] * ms_29[k]
                      + ns_29[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ab_x, ab_y, ab_z, ms_29, ms_30, ns_30, \
                         ns_37, ns_38, ns_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * ms_29[k]
                      + ns_37[k];

            t_89[k] = ab_z[k] * ms_29[k]
                      + ns_38[k];

            t_90[k] = ab_x[k] * ms_30[k]
                      + ns_30[k];

            t_91[k] = ab_y[k] * ms_30[k]
                      + ns_38[k];

            t_92[k] = ab_z[k] * ms_30[k]
                      + ns_39[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, ms_31, ms_32, \
                         ns_31, ns_32, ns_39, ns_40, ns_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = ab_x[k] * ms_31[k]
                      + ns_31[k];

            t_94[k] = ab_y[k] * ms_31[k]
                      + ns_39[k];

            t_95[k] = ab_z[k] * ms_31[k]
                      + ns_40[k];

            t_96[k] = ab_x[k] * ms_32[k]
                      + ns_32[k];

            t_97[k] = ab_y[k] * ms_32[k]
                      + ns_40[k];

            t_98[k] = ab_z[k] * ms_32[k]
                      + ns_41[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, ab_x, ab_y, ab_z, ms_33, \
                         ms_34, ns_33, ns_34, ns_41, ns_42, ns_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * ms_33[k]
                      + ns_33[k];

            t_100[k] = ab_y[k] * ms_33[k]
                       + ns_41[k];

            t_101[k] = ab_z[k] * ms_33[k]
                       + ns_42[k];

            t_102[k] = ab_x[k] * ms_34[k]
                       + ns_34[k];

            t_103[k] = ab_y[k] * ms_34[k]
                       + ns_42[k];

            t_104[k] = ab_z[k] * ms_34[k]
                       + ns_43[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, ms_35, ms_36, \
                         ns_35, ns_36, ns_43, ns_44, ns_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * ms_35[k]
                       + ns_35[k];

            t_106[k] = ab_y[k] * ms_35[k]
                       + ns_43[k];

            t_107[k] = ab_z[k] * ms_35[k]
                       + ns_44[k];

            t_108[k] = ab_x[k] * ms_36[k]
                       + ns_36[k];

            t_109[k] = ab_y[k] * ms_36[k]
                       + ns_45[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, ab_x, ab_y, ab_z, ms_36, \
                         ms_37, ms_38, ns_37, ns_38, ns_46, ns_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_z[k] * ms_36[k]
                       + ns_46[k];

            t_111[k] = ab_x[k] * ms_37[k]
                       + ns_37[k];

            t_112[k] = ab_y[k] * ms_37[k]
                       + ns_46[k];

            t_113[k] = ab_z[k] * ms_37[k]
                       + ns_47[k];

            t_114[k] = ab_x[k] * ms_38[k]
                       + ns_38[k];

            t_115[k] = ab_y[k] * ms_38[k]
                       + ns_47[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, ab_x, ab_y, ab_z, ms_38, \
                         ms_39, ms_40, ns_39, ns_40, ns_48, ns_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_z[k] * ms_38[k]
                       + ns_48[k];

            t_117[k] = ab_x[k] * ms_39[k]
                       + ns_39[k];

            t_118[k] = ab_y[k] * ms_39[k]
                       + ns_48[k];

            t_119[k] = ab_z[k] * ms_39[k]
                       + ns_49[k];

            t_120[k] = ab_x[k] * ms_40[k]
                       + ns_40[k];

            t_121[k] = ab_y[k] * ms_40[k]
                       + ns_49[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, ab_x, ab_y, ab_z, ms_40, \
                         ms_41, ms_42, ns_41, ns_42, ns_50, ns_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_z[k] * ms_40[k]
                       + ns_50[k];

            t_123[k] = ab_x[k] * ms_41[k]
                       + ns_41[k];

            t_124[k] = ab_y[k] * ms_41[k]
                       + ns_50[k];

            t_125[k] = ab_z[k] * ms_41[k]
                       + ns_51[k];

            t_126[k] = ab_x[k] * ms_42[k]
                       + ns_42[k];

            t_127[k] = ab_y[k] * ms_42[k]
                       + ns_51[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, ab_x, ab_y, ab_z, ms_42, \
                         ms_43, ms_44, ns_43, ns_44, ns_52, ns_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = ab_z[k] * ms_42[k]
                       + ns_52[k];

            t_129[k] = ab_x[k] * ms_43[k]
                       + ns_43[k];

            t_130[k] = ab_y[k] * ms_43[k]
                       + ns_52[k];

            t_131[k] = ab_z[k] * ms_43[k]
                       + ns_53[k];

            t_132[k] = ab_x[k] * ms_44[k]
                       + ns_44[k];

            t_133[k] = ab_y[k] * ms_44[k]
                       + ns_53[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, ms_44, ms_45, \
                         ms_46, ns_45, ns_46, ns_54, ns_55, ns_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = ab_z[k] * ms_44[k]
                       + ns_54[k];

            t_135[k] = ab_x[k] * ms_45[k]
                       + ns_45[k];

            t_136[k] = ab_y[k] * ms_45[k]
                       + ns_55[k];

            t_137[k] = ab_z[k] * ms_45[k]
                       + ns_56[k];

            t_138[k] = ab_x[k] * ms_46[k]
                       + ns_46[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, ms_46, ms_47, \
                         ns_47, ns_56, ns_57, ns_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_y[k] * ms_46[k]
                       + ns_56[k];

            t_140[k] = ab_z[k] * ms_46[k]
                       + ns_57[k];

            t_141[k] = ab_x[k] * ms_47[k]
                       + ns_47[k];

            t_142[k] = ab_y[k] * ms_47[k]
                       + ns_57[k];

            t_143[k] = ab_z[k] * ms_47[k]
                       + ns_58[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, ms_48, \
                         ms_49, ns_48, ns_49, ns_58, ns_59, ns_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * ms_48[k]
                       + ns_48[k];

            t_145[k] = ab_y[k] * ms_48[k]
                       + ns_58[k];

            t_146[k] = ab_z[k] * ms_48[k]
                       + ns_59[k];

            t_147[k] = ab_x[k] * ms_49[k]
                       + ns_49[k];

            t_148[k] = ab_y[k] * ms_49[k]
                       + ns_59[k];

            t_149[k] = ab_z[k] * ms_49[k]
                       + ns_60[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, ab_x, ab_y, ab_z, ms_50, \
                         ms_51, ns_50, ns_51, ns_60, ns_61, ns_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ms_50[k]
                       + ns_50[k];

            t_151[k] = ab_y[k] * ms_50[k]
                       + ns_60[k];

            t_152[k] = ab_z[k] * ms_50[k]
                       + ns_61[k];

            t_153[k] = ab_x[k] * ms_51[k]
                       + ns_51[k];

            t_154[k] = ab_y[k] * ms_51[k]
                       + ns_61[k];

            t_155[k] = ab_z[k] * ms_51[k]
                       + ns_62[k];
        }
    }
}

static auto
compute_hrr_mp_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ms, const size_t ns, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ms_52 = buffer.data(ms + 52 * ncomps + c);
        const auto *ms_53 = buffer.data(ms + 53 * ncomps + c);
        const auto *ms_54 = buffer.data(ms + 54 * ncomps + c);

        const auto *ns_52 = buffer.data(ns + 52 * ncomps + c);
        const auto *ns_53 = buffer.data(ns + 53 * ncomps + c);
        const auto *ns_54 = buffer.data(ns + 54 * ncomps + c);
        const auto *ns_62 = buffer.data(ns + 62 * ncomps + c);
        const auto *ns_63 = buffer.data(ns + 63 * ncomps + c);
        const auto *ns_64 = buffer.data(ns + 64 * ncomps + c);
        const auto *ns_65 = buffer.data(ns + 65 * ncomps + c);

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, ms_52, \
                         ms_53, ns_52, ns_53, ns_62, ns_63, ns_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = ab_x[k] * ms_52[k]
                       + ns_52[k];

            t_157[k] = ab_y[k] * ms_52[k]
                       + ns_62[k];

            t_158[k] = ab_z[k] * ms_52[k]
                       + ns_63[k];

            t_159[k] = ab_x[k] * ms_53[k]
                       + ns_53[k];

            t_160[k] = ab_y[k] * ms_53[k]
                       + ns_63[k];

            t_161[k] = ab_z[k] * ms_53[k]
                       + ns_64[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, ab_x, ab_y, ab_z, ms_54, ns_54, ns_64, \
                         ns_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * ms_54[k]
                       + ns_54[k];

            t_163[k] = ab_y[k] * ms_54[k]
                       + ns_64[k];

            t_164[k] = ab_z[k] * ms_54[k]
                       + ns_65[k];
        }
    }
}

auto
compute_hrr_mp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ms, const size_t ns, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_mp_piece0(buffer, coordinates, target, ms, ns, ncomps, nmax);

    compute_hrr_mp_piece1(buffer, coordinates, target, ms, ns, ncomps, nmax);
}

}  // namespace simdtrf
