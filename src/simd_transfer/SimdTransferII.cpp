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


#include "SimdTransferII.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_ii_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ih, const size_t kh,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ih_0 = buffer.data(ih + 0 * ncomps + c);
        const auto *ih_1 = buffer.data(ih + 1 * ncomps + c);
        const auto *ih_2 = buffer.data(ih + 2 * ncomps + c);
        const auto *ih_3 = buffer.data(ih + 3 * ncomps + c);
        const auto *ih_4 = buffer.data(ih + 4 * ncomps + c);
        const auto *ih_5 = buffer.data(ih + 5 * ncomps + c);
        const auto *ih_6 = buffer.data(ih + 6 * ncomps + c);
        const auto *ih_7 = buffer.data(ih + 7 * ncomps + c);
        const auto *ih_8 = buffer.data(ih + 8 * ncomps + c);
        const auto *ih_9 = buffer.data(ih + 9 * ncomps + c);
        const auto *ih_10 = buffer.data(ih + 10 * ncomps + c);
        const auto *ih_11 = buffer.data(ih + 11 * ncomps + c);
        const auto *ih_12 = buffer.data(ih + 12 * ncomps + c);
        const auto *ih_13 = buffer.data(ih + 13 * ncomps + c);
        const auto *ih_14 = buffer.data(ih + 14 * ncomps + c);
        const auto *ih_15 = buffer.data(ih + 15 * ncomps + c);
        const auto *ih_16 = buffer.data(ih + 16 * ncomps + c);
        const auto *ih_17 = buffer.data(ih + 17 * ncomps + c);
        const auto *ih_18 = buffer.data(ih + 18 * ncomps + c);
        const auto *ih_19 = buffer.data(ih + 19 * ncomps + c);
        const auto *ih_20 = buffer.data(ih + 20 * ncomps + c);
        const auto *ih_21 = buffer.data(ih + 21 * ncomps + c);
        const auto *ih_22 = buffer.data(ih + 22 * ncomps + c);
        const auto *ih_23 = buffer.data(ih + 23 * ncomps + c);
        const auto *ih_24 = buffer.data(ih + 24 * ncomps + c);
        const auto *ih_25 = buffer.data(ih + 25 * ncomps + c);
        const auto *ih_26 = buffer.data(ih + 26 * ncomps + c);
        const auto *ih_27 = buffer.data(ih + 27 * ncomps + c);
        const auto *ih_28 = buffer.data(ih + 28 * ncomps + c);
        const auto *ih_29 = buffer.data(ih + 29 * ncomps + c);
        const auto *ih_30 = buffer.data(ih + 30 * ncomps + c);
        const auto *ih_31 = buffer.data(ih + 31 * ncomps + c);
        const auto *ih_32 = buffer.data(ih + 32 * ncomps + c);
        const auto *ih_33 = buffer.data(ih + 33 * ncomps + c);
        const auto *ih_34 = buffer.data(ih + 34 * ncomps + c);
        const auto *ih_35 = buffer.data(ih + 35 * ncomps + c);
        const auto *ih_36 = buffer.data(ih + 36 * ncomps + c);
        const auto *ih_37 = buffer.data(ih + 37 * ncomps + c);
        const auto *ih_38 = buffer.data(ih + 38 * ncomps + c);
        const auto *ih_39 = buffer.data(ih + 39 * ncomps + c);
        const auto *ih_40 = buffer.data(ih + 40 * ncomps + c);
        const auto *ih_41 = buffer.data(ih + 41 * ncomps + c);
        const auto *ih_42 = buffer.data(ih + 42 * ncomps + c);
        const auto *ih_43 = buffer.data(ih + 43 * ncomps + c);
        const auto *ih_44 = buffer.data(ih + 44 * ncomps + c);
        const auto *ih_45 = buffer.data(ih + 45 * ncomps + c);
        const auto *ih_46 = buffer.data(ih + 46 * ncomps + c);
        const auto *ih_47 = buffer.data(ih + 47 * ncomps + c);
        const auto *ih_48 = buffer.data(ih + 48 * ncomps + c);
        const auto *ih_49 = buffer.data(ih + 49 * ncomps + c);
        const auto *ih_50 = buffer.data(ih + 50 * ncomps + c);
        const auto *ih_51 = buffer.data(ih + 51 * ncomps + c);
        const auto *ih_52 = buffer.data(ih + 52 * ncomps + c);
        const auto *ih_53 = buffer.data(ih + 53 * ncomps + c);
        const auto *ih_54 = buffer.data(ih + 54 * ncomps + c);
        const auto *ih_55 = buffer.data(ih + 55 * ncomps + c);
        const auto *ih_56 = buffer.data(ih + 56 * ncomps + c);
        const auto *ih_57 = buffer.data(ih + 57 * ncomps + c);
        const auto *ih_58 = buffer.data(ih + 58 * ncomps + c);
        const auto *ih_59 = buffer.data(ih + 59 * ncomps + c);
        const auto *ih_60 = buffer.data(ih + 60 * ncomps + c);
        const auto *ih_61 = buffer.data(ih + 61 * ncomps + c);
        const auto *ih_62 = buffer.data(ih + 62 * ncomps + c);
        const auto *ih_63 = buffer.data(ih + 63 * ncomps + c);
        const auto *ih_64 = buffer.data(ih + 64 * ncomps + c);
        const auto *ih_65 = buffer.data(ih + 65 * ncomps + c);
        const auto *ih_66 = buffer.data(ih + 66 * ncomps + c);
        const auto *ih_67 = buffer.data(ih + 67 * ncomps + c);
        const auto *ih_68 = buffer.data(ih + 68 * ncomps + c);
        const auto *ih_69 = buffer.data(ih + 69 * ncomps + c);
        const auto *ih_70 = buffer.data(ih + 70 * ncomps + c);
        const auto *ih_71 = buffer.data(ih + 71 * ncomps + c);
        const auto *ih_72 = buffer.data(ih + 72 * ncomps + c);
        const auto *ih_73 = buffer.data(ih + 73 * ncomps + c);
        const auto *ih_74 = buffer.data(ih + 74 * ncomps + c);
        const auto *ih_75 = buffer.data(ih + 75 * ncomps + c);
        const auto *ih_76 = buffer.data(ih + 76 * ncomps + c);
        const auto *ih_77 = buffer.data(ih + 77 * ncomps + c);
        const auto *ih_78 = buffer.data(ih + 78 * ncomps + c);
        const auto *ih_79 = buffer.data(ih + 79 * ncomps + c);
        const auto *ih_80 = buffer.data(ih + 80 * ncomps + c);
        const auto *ih_81 = buffer.data(ih + 81 * ncomps + c);
        const auto *ih_82 = buffer.data(ih + 82 * ncomps + c);
        const auto *ih_83 = buffer.data(ih + 83 * ncomps + c);
        const auto *ih_84 = buffer.data(ih + 84 * ncomps + c);
        const auto *ih_85 = buffer.data(ih + 85 * ncomps + c);
        const auto *ih_86 = buffer.data(ih + 86 * ncomps + c);
        const auto *ih_87 = buffer.data(ih + 87 * ncomps + c);
        const auto *ih_88 = buffer.data(ih + 88 * ncomps + c);
        const auto *ih_89 = buffer.data(ih + 89 * ncomps + c);
        const auto *ih_90 = buffer.data(ih + 90 * ncomps + c);
        const auto *ih_91 = buffer.data(ih + 91 * ncomps + c);
        const auto *ih_92 = buffer.data(ih + 92 * ncomps + c);
        const auto *ih_93 = buffer.data(ih + 93 * ncomps + c);
        const auto *ih_94 = buffer.data(ih + 94 * ncomps + c);
        const auto *ih_95 = buffer.data(ih + 95 * ncomps + c);
        const auto *ih_96 = buffer.data(ih + 96 * ncomps + c);
        const auto *ih_97 = buffer.data(ih + 97 * ncomps + c);
        const auto *ih_98 = buffer.data(ih + 98 * ncomps + c);
        const auto *ih_99 = buffer.data(ih + 99 * ncomps + c);
        const auto *ih_100 = buffer.data(ih + 100 * ncomps + c);
        const auto *ih_101 = buffer.data(ih + 101 * ncomps + c);
        const auto *ih_102 = buffer.data(ih + 102 * ncomps + c);
        const auto *ih_103 = buffer.data(ih + 103 * ncomps + c);
        const auto *ih_104 = buffer.data(ih + 104 * ncomps + c);
        const auto *ih_105 = buffer.data(ih + 105 * ncomps + c);
        const auto *ih_106 = buffer.data(ih + 106 * ncomps + c);
        const auto *ih_107 = buffer.data(ih + 107 * ncomps + c);
        const auto *ih_108 = buffer.data(ih + 108 * ncomps + c);
        const auto *ih_109 = buffer.data(ih + 109 * ncomps + c);

        const auto *kh_0 = buffer.data(kh + 0 * ncomps + c);
        const auto *kh_1 = buffer.data(kh + 1 * ncomps + c);
        const auto *kh_2 = buffer.data(kh + 2 * ncomps + c);
        const auto *kh_3 = buffer.data(kh + 3 * ncomps + c);
        const auto *kh_4 = buffer.data(kh + 4 * ncomps + c);
        const auto *kh_5 = buffer.data(kh + 5 * ncomps + c);
        const auto *kh_6 = buffer.data(kh + 6 * ncomps + c);
        const auto *kh_7 = buffer.data(kh + 7 * ncomps + c);
        const auto *kh_8 = buffer.data(kh + 8 * ncomps + c);
        const auto *kh_9 = buffer.data(kh + 9 * ncomps + c);
        const auto *kh_10 = buffer.data(kh + 10 * ncomps + c);
        const auto *kh_11 = buffer.data(kh + 11 * ncomps + c);
        const auto *kh_12 = buffer.data(kh + 12 * ncomps + c);
        const auto *kh_13 = buffer.data(kh + 13 * ncomps + c);
        const auto *kh_14 = buffer.data(kh + 14 * ncomps + c);
        const auto *kh_15 = buffer.data(kh + 15 * ncomps + c);
        const auto *kh_16 = buffer.data(kh + 16 * ncomps + c);
        const auto *kh_17 = buffer.data(kh + 17 * ncomps + c);
        const auto *kh_18 = buffer.data(kh + 18 * ncomps + c);
        const auto *kh_19 = buffer.data(kh + 19 * ncomps + c);
        const auto *kh_20 = buffer.data(kh + 20 * ncomps + c);
        const auto *kh_21 = buffer.data(kh + 21 * ncomps + c);
        const auto *kh_22 = buffer.data(kh + 22 * ncomps + c);
        const auto *kh_23 = buffer.data(kh + 23 * ncomps + c);
        const auto *kh_24 = buffer.data(kh + 24 * ncomps + c);
        const auto *kh_25 = buffer.data(kh + 25 * ncomps + c);
        const auto *kh_26 = buffer.data(kh + 26 * ncomps + c);
        const auto *kh_27 = buffer.data(kh + 27 * ncomps + c);
        const auto *kh_28 = buffer.data(kh + 28 * ncomps + c);
        const auto *kh_29 = buffer.data(kh + 29 * ncomps + c);
        const auto *kh_30 = buffer.data(kh + 30 * ncomps + c);
        const auto *kh_31 = buffer.data(kh + 31 * ncomps + c);
        const auto *kh_32 = buffer.data(kh + 32 * ncomps + c);
        const auto *kh_33 = buffer.data(kh + 33 * ncomps + c);
        const auto *kh_34 = buffer.data(kh + 34 * ncomps + c);
        const auto *kh_35 = buffer.data(kh + 35 * ncomps + c);
        const auto *kh_36 = buffer.data(kh + 36 * ncomps + c);
        const auto *kh_37 = buffer.data(kh + 37 * ncomps + c);
        const auto *kh_38 = buffer.data(kh + 38 * ncomps + c);
        const auto *kh_39 = buffer.data(kh + 39 * ncomps + c);
        const auto *kh_40 = buffer.data(kh + 40 * ncomps + c);
        const auto *kh_41 = buffer.data(kh + 41 * ncomps + c);
        const auto *kh_42 = buffer.data(kh + 42 * ncomps + c);
        const auto *kh_43 = buffer.data(kh + 43 * ncomps + c);
        const auto *kh_44 = buffer.data(kh + 44 * ncomps + c);
        const auto *kh_45 = buffer.data(kh + 45 * ncomps + c);
        const auto *kh_46 = buffer.data(kh + 46 * ncomps + c);
        const auto *kh_47 = buffer.data(kh + 47 * ncomps + c);
        const auto *kh_48 = buffer.data(kh + 48 * ncomps + c);
        const auto *kh_49 = buffer.data(kh + 49 * ncomps + c);
        const auto *kh_50 = buffer.data(kh + 50 * ncomps + c);
        const auto *kh_51 = buffer.data(kh + 51 * ncomps + c);
        const auto *kh_52 = buffer.data(kh + 52 * ncomps + c);
        const auto *kh_53 = buffer.data(kh + 53 * ncomps + c);
        const auto *kh_54 = buffer.data(kh + 54 * ncomps + c);
        const auto *kh_55 = buffer.data(kh + 55 * ncomps + c);
        const auto *kh_56 = buffer.data(kh + 56 * ncomps + c);
        const auto *kh_57 = buffer.data(kh + 57 * ncomps + c);
        const auto *kh_58 = buffer.data(kh + 58 * ncomps + c);
        const auto *kh_59 = buffer.data(kh + 59 * ncomps + c);
        const auto *kh_60 = buffer.data(kh + 60 * ncomps + c);
        const auto *kh_61 = buffer.data(kh + 61 * ncomps + c);
        const auto *kh_62 = buffer.data(kh + 62 * ncomps + c);
        const auto *kh_63 = buffer.data(kh + 63 * ncomps + c);
        const auto *kh_64 = buffer.data(kh + 64 * ncomps + c);
        const auto *kh_65 = buffer.data(kh + 65 * ncomps + c);
        const auto *kh_66 = buffer.data(kh + 66 * ncomps + c);
        const auto *kh_67 = buffer.data(kh + 67 * ncomps + c);
        const auto *kh_68 = buffer.data(kh + 68 * ncomps + c);
        const auto *kh_69 = buffer.data(kh + 69 * ncomps + c);
        const auto *kh_70 = buffer.data(kh + 70 * ncomps + c);
        const auto *kh_71 = buffer.data(kh + 71 * ncomps + c);
        const auto *kh_72 = buffer.data(kh + 72 * ncomps + c);
        const auto *kh_73 = buffer.data(kh + 73 * ncomps + c);
        const auto *kh_74 = buffer.data(kh + 74 * ncomps + c);
        const auto *kh_75 = buffer.data(kh + 75 * ncomps + c);
        const auto *kh_76 = buffer.data(kh + 76 * ncomps + c);
        const auto *kh_77 = buffer.data(kh + 77 * ncomps + c);
        const auto *kh_78 = buffer.data(kh + 78 * ncomps + c);
        const auto *kh_79 = buffer.data(kh + 79 * ncomps + c);
        const auto *kh_80 = buffer.data(kh + 80 * ncomps + c);
        const auto *kh_81 = buffer.data(kh + 81 * ncomps + c);
        const auto *kh_82 = buffer.data(kh + 82 * ncomps + c);
        const auto *kh_83 = buffer.data(kh + 83 * ncomps + c);
        const auto *kh_84 = buffer.data(kh + 84 * ncomps + c);
        const auto *kh_85 = buffer.data(kh + 85 * ncomps + c);
        const auto *kh_86 = buffer.data(kh + 86 * ncomps + c);
        const auto *kh_87 = buffer.data(kh + 87 * ncomps + c);
        const auto *kh_88 = buffer.data(kh + 88 * ncomps + c);
        const auto *kh_89 = buffer.data(kh + 89 * ncomps + c);
        const auto *kh_90 = buffer.data(kh + 90 * ncomps + c);
        const auto *kh_91 = buffer.data(kh + 91 * ncomps + c);
        const auto *kh_92 = buffer.data(kh + 92 * ncomps + c);
        const auto *kh_93 = buffer.data(kh + 93 * ncomps + c);
        const auto *kh_94 = buffer.data(kh + 94 * ncomps + c);
        const auto *kh_95 = buffer.data(kh + 95 * ncomps + c);
        const auto *kh_96 = buffer.data(kh + 96 * ncomps + c);
        const auto *kh_97 = buffer.data(kh + 97 * ncomps + c);
        const auto *kh_98 = buffer.data(kh + 98 * ncomps + c);
        const auto *kh_99 = buffer.data(kh + 99 * ncomps + c);
        const auto *kh_100 = buffer.data(kh + 100 * ncomps + c);
        const auto *kh_101 = buffer.data(kh + 101 * ncomps + c);
        const auto *kh_102 = buffer.data(kh + 102 * ncomps + c);
        const auto *kh_103 = buffer.data(kh + 103 * ncomps + c);
        const auto *kh_104 = buffer.data(kh + 104 * ncomps + c);
        const auto *kh_105 = buffer.data(kh + 105 * ncomps + c);
        const auto *kh_106 = buffer.data(kh + 106 * ncomps + c);
        const auto *kh_107 = buffer.data(kh + 107 * ncomps + c);
        const auto *kh_108 = buffer.data(kh + 108 * ncomps + c);
        const auto *kh_109 = buffer.data(kh + 109 * ncomps + c);
        const auto *kh_125 = buffer.data(kh + 125 * ncomps + c);
        const auto *kh_141 = buffer.data(kh + 141 * ncomps + c);
        const auto *kh_142 = buffer.data(kh + 142 * ncomps + c);
        const auto *kh_143 = buffer.data(kh + 143 * ncomps + c);
        const auto *kh_144 = buffer.data(kh + 144 * ncomps + c);
        const auto *kh_145 = buffer.data(kh + 145 * ncomps + c);
        const auto *kh_146 = buffer.data(kh + 146 * ncomps + c);
        const auto *kh_162 = buffer.data(kh + 162 * ncomps + c);
        const auto *kh_163 = buffer.data(kh + 163 * ncomps + c);
        const auto *kh_164 = buffer.data(kh + 164 * ncomps + c);
        const auto *kh_165 = buffer.data(kh + 165 * ncomps + c);
        const auto *kh_166 = buffer.data(kh + 166 * ncomps + c);
        const auto *kh_167 = buffer.data(kh + 167 * ncomps + c);
        const auto *kh_188 = buffer.data(kh + 188 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ih_0, ih_1, ih_2, ih_3, ih_4, kh_0, \
                         kh_1, kh_2, kh_3, kh_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ih_0[k]
                     + kh_0[k];

            t_1[k] = ab_x[k] * ih_1[k]
                     + kh_1[k];

            t_2[k] = ab_x[k] * ih_2[k]
                     + kh_2[k];

            t_3[k] = ab_x[k] * ih_3[k]
                     + kh_3[k];

            t_4[k] = ab_x[k] * ih_4[k]
                     + kh_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ih_5, ih_6, ih_7, ih_8, ih_9, kh_5, \
                         kh_6, kh_7, kh_8, kh_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * ih_5[k]
                     + kh_5[k];

            t_6[k] = ab_x[k] * ih_6[k]
                     + kh_6[k];

            t_7[k] = ab_x[k] * ih_7[k]
                     + kh_7[k];

            t_8[k] = ab_x[k] * ih_8[k]
                     + kh_8[k];

            t_9[k] = ab_x[k] * ih_9[k]
                     + kh_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ih_10, ih_11, ih_12, ih_13, \
                         ih_14, kh_10, kh_11, kh_12, kh_13, kh_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * ih_10[k]
                      + kh_10[k];

            t_11[k] = ab_x[k] * ih_11[k]
                      + kh_11[k];

            t_12[k] = ab_x[k] * ih_12[k]
                      + kh_12[k];

            t_13[k] = ab_x[k] * ih_13[k]
                      + kh_13[k];

            t_14[k] = ab_x[k] * ih_14[k]
                      + kh_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ih_15, ih_16, ih_17, ih_18, \
                         ih_19, kh_15, kh_16, kh_17, kh_18, kh_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * ih_15[k]
                      + kh_15[k];

            t_16[k] = ab_x[k] * ih_16[k]
                      + kh_16[k];

            t_17[k] = ab_x[k] * ih_17[k]
                      + kh_17[k];

            t_18[k] = ab_x[k] * ih_18[k]
                      + kh_18[k];

            t_19[k] = ab_x[k] * ih_19[k]
                      + kh_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_y, ih_15, ih_16, ih_17, ih_20, \
                         kh_20, kh_36, kh_37, kh_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * ih_20[k]
                      + kh_20[k];

            t_21[k] = ab_y[k] * ih_15[k]
                      + kh_36[k];

            t_22[k] = ab_y[k] * ih_16[k]
                      + kh_37[k];

            t_23[k] = ab_y[k] * ih_17[k]
                      + kh_38[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, ab_y, ab_z, ih_18, ih_19, ih_20, kh_39, \
                         kh_40, kh_41, kh_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_y[k] * ih_18[k]
                      + kh_39[k];

            t_25[k] = ab_y[k] * ih_19[k]
                      + kh_40[k];

            t_26[k] = ab_y[k] * ih_20[k]
                      + kh_41[k];

            t_27[k] = ab_z[k] * ih_20[k]
                      + kh_62[k];
        }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, ab_x, ih_21, ih_22, ih_23, ih_24, \
                         ih_25, kh_21, kh_22, kh_23, kh_24, kh_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_28[k] = ab_x[k] * ih_21[k]
                      + kh_21[k];

            t_29[k] = ab_x[k] * ih_22[k]
                      + kh_22[k];

            t_30[k] = ab_x[k] * ih_23[k]
                      + kh_23[k];

            t_31[k] = ab_x[k] * ih_24[k]
                      + kh_24[k];

            t_32[k] = ab_x[k] * ih_25[k]
                      + kh_25[k];
        }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, ab_x, ih_26, ih_27, ih_28, ih_29, \
                         ih_30, kh_26, kh_27, kh_28, kh_29, kh_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_33[k] = ab_x[k] * ih_26[k]
                      + kh_26[k];

            t_34[k] = ab_x[k] * ih_27[k]
                      + kh_27[k];

            t_35[k] = ab_x[k] * ih_28[k]
                      + kh_28[k];

            t_36[k] = ab_x[k] * ih_29[k]
                      + kh_29[k];

            t_37[k] = ab_x[k] * ih_30[k]
                      + kh_30[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, ab_x, ih_31, ih_32, ih_33, ih_34, \
                         ih_35, kh_31, kh_32, kh_33, kh_34, kh_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_x[k] * ih_31[k]
                      + kh_31[k];

            t_39[k] = ab_x[k] * ih_32[k]
                      + kh_32[k];

            t_40[k] = ab_x[k] * ih_33[k]
                      + kh_33[k];

            t_41[k] = ab_x[k] * ih_34[k]
                      + kh_34[k];

            t_42[k] = ab_x[k] * ih_35[k]
                      + kh_35[k];
        }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, ab_x, ih_36, ih_37, ih_38, ih_39, \
                         ih_40, kh_36, kh_37, kh_38, kh_39, kh_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_43[k] = ab_x[k] * ih_36[k]
                      + kh_36[k];

            t_44[k] = ab_x[k] * ih_37[k]
                      + kh_37[k];

            t_45[k] = ab_x[k] * ih_38[k]
                      + kh_38[k];

            t_46[k] = ab_x[k] * ih_39[k]
                      + kh_39[k];

            t_47[k] = ab_x[k] * ih_40[k]
                      + kh_40[k];
        }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, ab_x, ab_y, ih_36, ih_37, ih_38, ih_41, \
                         kh_41, kh_78, kh_79, kh_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_48[k] = ab_x[k] * ih_41[k]
                      + kh_41[k];

            t_49[k] = ab_y[k] * ih_36[k]
                      + kh_78[k];

            t_50[k] = ab_y[k] * ih_37[k]
                      + kh_79[k];

            t_51[k] = ab_y[k] * ih_38[k]
                      + kh_80[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, ab_y, ab_z, ih_39, ih_40, ih_41, kh_81, \
                         kh_82, kh_83, kh_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_y[k] * ih_39[k]
                      + kh_81[k];

            t_53[k] = ab_y[k] * ih_40[k]
                      + kh_82[k];

            t_54[k] = ab_y[k] * ih_41[k]
                      + kh_83[k];

            t_55[k] = ab_z[k] * ih_41[k]
                      + kh_104[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, ab_x, ih_42, ih_43, ih_44, ih_45, \
                         ih_46, kh_42, kh_43, kh_44, kh_45, kh_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = ab_x[k] * ih_42[k]
                      + kh_42[k];

            t_57[k] = ab_x[k] * ih_43[k]
                      + kh_43[k];

            t_58[k] = ab_x[k] * ih_44[k]
                      + kh_44[k];

            t_59[k] = ab_x[k] * ih_45[k]
                      + kh_45[k];

            t_60[k] = ab_x[k] * ih_46[k]
                      + kh_46[k];
        }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, ab_x, ih_47, ih_48, ih_49, ih_50, \
                         ih_51, kh_47, kh_48, kh_49, kh_50, kh_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_61[k] = ab_x[k] * ih_47[k]
                      + kh_47[k];

            t_62[k] = ab_x[k] * ih_48[k]
                      + kh_48[k];

            t_63[k] = ab_x[k] * ih_49[k]
                      + kh_49[k];

            t_64[k] = ab_x[k] * ih_50[k]
                      + kh_50[k];

            t_65[k] = ab_x[k] * ih_51[k]
                      + kh_51[k];
        }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ab_x, ih_52, ih_53, ih_54, ih_55, \
                         ih_56, kh_52, kh_53, kh_54, kh_55, kh_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_66[k] = ab_x[k] * ih_52[k]
                      + kh_52[k];

            t_67[k] = ab_x[k] * ih_53[k]
                      + kh_53[k];

            t_68[k] = ab_x[k] * ih_54[k]
                      + kh_54[k];

            t_69[k] = ab_x[k] * ih_55[k]
                      + kh_55[k];

            t_70[k] = ab_x[k] * ih_56[k]
                      + kh_56[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ab_x, ih_57, ih_58, ih_59, ih_60, \
                         ih_61, kh_57, kh_58, kh_59, kh_60, kh_61 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_x[k] * ih_57[k]
                      + kh_57[k];

            t_72[k] = ab_x[k] * ih_58[k]
                      + kh_58[k];

            t_73[k] = ab_x[k] * ih_59[k]
                      + kh_59[k];

            t_74[k] = ab_x[k] * ih_60[k]
                      + kh_60[k];

            t_75[k] = ab_x[k] * ih_61[k]
                      + kh_61[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, ih_57, ih_58, ih_59, ih_62, \
                         kh_62, kh_99, kh_100, kh_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_x[k] * ih_62[k]
                      + kh_62[k];

            t_77[k] = ab_y[k] * ih_57[k]
                      + kh_99[k];

            t_78[k] = ab_y[k] * ih_58[k]
                      + kh_100[k];

            t_79[k] = ab_y[k] * ih_59[k]
                      + kh_101[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_y, ab_z, ih_60, ih_61, ih_62, kh_102, \
                         kh_103, kh_104, kh_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_y[k] * ih_60[k]
                      + kh_102[k];

            t_81[k] = ab_y[k] * ih_61[k]
                      + kh_103[k];

            t_82[k] = ab_y[k] * ih_62[k]
                      + kh_104[k];

            t_83[k] = ab_z[k] * ih_62[k]
                      + kh_125[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_x, ih_63, ih_64, ih_65, ih_66, \
                         ih_67, kh_63, kh_64, kh_65, kh_66, kh_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * ih_63[k]
                      + kh_63[k];

            t_85[k] = ab_x[k] * ih_64[k]
                      + kh_64[k];

            t_86[k] = ab_x[k] * ih_65[k]
                      + kh_65[k];

            t_87[k] = ab_x[k] * ih_66[k]
                      + kh_66[k];

            t_88[k] = ab_x[k] * ih_67[k]
                      + kh_67[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_x, ih_68, ih_69, ih_70, ih_71, \
                         ih_72, kh_68, kh_69, kh_70, kh_71, kh_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_x[k] * ih_68[k]
                      + kh_68[k];

            t_90[k] = ab_x[k] * ih_69[k]
                      + kh_69[k];

            t_91[k] = ab_x[k] * ih_70[k]
                      + kh_70[k];

            t_92[k] = ab_x[k] * ih_71[k]
                      + kh_71[k];

            t_93[k] = ab_x[k] * ih_72[k]
                      + kh_72[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_x, ih_73, ih_74, ih_75, ih_76, \
                         ih_77, kh_73, kh_74, kh_75, kh_76, kh_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = ab_x[k] * ih_73[k]
                      + kh_73[k];

            t_95[k] = ab_x[k] * ih_74[k]
                      + kh_74[k];

            t_96[k] = ab_x[k] * ih_75[k]
                      + kh_75[k];

            t_97[k] = ab_x[k] * ih_76[k]
                      + kh_76[k];

            t_98[k] = ab_x[k] * ih_77[k]
                      + kh_77[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_x, ih_78, ih_79, ih_80, ih_81, \
                         ih_82, kh_78, kh_79, kh_80, kh_81, kh_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * ih_78[k]
                      + kh_78[k];

            t_100[k] = ab_x[k] * ih_79[k]
                       + kh_79[k];

            t_101[k] = ab_x[k] * ih_80[k]
                       + kh_80[k];

            t_102[k] = ab_x[k] * ih_81[k]
                       + kh_81[k];

            t_103[k] = ab_x[k] * ih_82[k]
                       + kh_82[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_y, ih_78, ih_79, ih_80, ih_83, \
                         kh_83, kh_141, kh_142, kh_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_x[k] * ih_83[k]
                       + kh_83[k];

            t_105[k] = ab_y[k] * ih_78[k]
                       + kh_141[k];

            t_106[k] = ab_y[k] * ih_79[k]
                       + kh_142[k];

            t_107[k] = ab_y[k] * ih_80[k]
                       + kh_143[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, ab_y, ab_z, ih_81, ih_82, ih_83, kh_144, \
                         kh_145, kh_146, kh_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_y[k] * ih_81[k]
                       + kh_144[k];

            t_109[k] = ab_y[k] * ih_82[k]
                       + kh_145[k];

            t_110[k] = ab_y[k] * ih_83[k]
                       + kh_146[k];

            t_111[k] = ab_z[k] * ih_83[k]
                       + kh_167[k];
        }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, ab_x, ih_84, ih_85, ih_86, ih_87, \
                         ih_88, kh_84, kh_85, kh_86, kh_87, kh_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_112[k] = ab_x[k] * ih_84[k]
                       + kh_84[k];

            t_113[k] = ab_x[k] * ih_85[k]
                       + kh_85[k];

            t_114[k] = ab_x[k] * ih_86[k]
                       + kh_86[k];

            t_115[k] = ab_x[k] * ih_87[k]
                       + kh_87[k];

            t_116[k] = ab_x[k] * ih_88[k]
                       + kh_88[k];
        }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, ab_x, ih_89, ih_90, ih_91, ih_92, \
                         ih_93, kh_89, kh_90, kh_91, kh_92, kh_93 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_117[k] = ab_x[k] * ih_89[k]
                       + kh_89[k];

            t_118[k] = ab_x[k] * ih_90[k]
                       + kh_90[k];

            t_119[k] = ab_x[k] * ih_91[k]
                       + kh_91[k];

            t_120[k] = ab_x[k] * ih_92[k]
                       + kh_92[k];

            t_121[k] = ab_x[k] * ih_93[k]
                       + kh_93[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, ab_x, ih_94, ih_95, ih_96, ih_97, \
                         ih_98, kh_94, kh_95, kh_96, kh_97, kh_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_x[k] * ih_94[k]
                       + kh_94[k];

            t_123[k] = ab_x[k] * ih_95[k]
                       + kh_95[k];

            t_124[k] = ab_x[k] * ih_96[k]
                       + kh_96[k];

            t_125[k] = ab_x[k] * ih_97[k]
                       + kh_97[k];

            t_126[k] = ab_x[k] * ih_98[k]
                       + kh_98[k];
        }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, ab_x, ih_99, ih_100, ih_101, \
                         ih_102, ih_103, kh_99, kh_100, kh_101, kh_102, \
                         kh_103 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_127[k] = ab_x[k] * ih_99[k]
                       + kh_99[k];

            t_128[k] = ab_x[k] * ih_100[k]
                       + kh_100[k];

            t_129[k] = ab_x[k] * ih_101[k]
                       + kh_101[k];

            t_130[k] = ab_x[k] * ih_102[k]
                       + kh_102[k];

            t_131[k] = ab_x[k] * ih_103[k]
                       + kh_103[k];
        }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, ab_x, ab_y, ih_99, ih_100, ih_101, \
                         ih_104, kh_104, kh_162, kh_163, kh_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_132[k] = ab_x[k] * ih_104[k]
                       + kh_104[k];

            t_133[k] = ab_y[k] * ih_99[k]
                       + kh_162[k];

            t_134[k] = ab_y[k] * ih_100[k]
                       + kh_163[k];

            t_135[k] = ab_y[k] * ih_101[k]
                       + kh_164[k];
        }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, ab_y, ab_z, ih_102, ih_103, ih_104, \
                         kh_165, kh_166, kh_167, kh_188 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_136[k] = ab_y[k] * ih_102[k]
                       + kh_165[k];

            t_137[k] = ab_y[k] * ih_103[k]
                       + kh_166[k];

            t_138[k] = ab_y[k] * ih_104[k]
                       + kh_167[k];

            t_139[k] = ab_z[k] * ih_104[k]
                       + kh_188[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, ih_105, ih_106, ih_107, \
                         ih_108, ih_109, kh_105, kh_106, kh_107, kh_108, \
                         kh_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * ih_105[k]
                       + kh_105[k];

            t_141[k] = ab_x[k] * ih_106[k]
                       + kh_106[k];

            t_142[k] = ab_x[k] * ih_107[k]
                       + kh_107[k];

            t_143[k] = ab_x[k] * ih_108[k]
                       + kh_108[k];

            t_144[k] = ab_x[k] * ih_109[k]
                       + kh_109[k];
        }
    }
}

static auto
compute_hrr_ii_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ih, const size_t kh,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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
        auto *t_234 = buffer.data(target + 234 * ncomps + c);
        auto *t_235 = buffer.data(target + 235 * ncomps + c);
        auto *t_236 = buffer.data(target + 236 * ncomps + c);
        auto *t_237 = buffer.data(target + 237 * ncomps + c);
        auto *t_238 = buffer.data(target + 238 * ncomps + c);
        auto *t_239 = buffer.data(target + 239 * ncomps + c);
        auto *t_240 = buffer.data(target + 240 * ncomps + c);
        auto *t_241 = buffer.data(target + 241 * ncomps + c);
        auto *t_242 = buffer.data(target + 242 * ncomps + c);
        auto *t_243 = buffer.data(target + 243 * ncomps + c);
        auto *t_244 = buffer.data(target + 244 * ncomps + c);
        auto *t_245 = buffer.data(target + 245 * ncomps + c);
        auto *t_246 = buffer.data(target + 246 * ncomps + c);
        auto *t_247 = buffer.data(target + 247 * ncomps + c);
        auto *t_248 = buffer.data(target + 248 * ncomps + c);
        auto *t_249 = buffer.data(target + 249 * ncomps + c);
        auto *t_250 = buffer.data(target + 250 * ncomps + c);
        auto *t_251 = buffer.data(target + 251 * ncomps + c);
        auto *t_252 = buffer.data(target + 252 * ncomps + c);
        auto *t_253 = buffer.data(target + 253 * ncomps + c);
        auto *t_254 = buffer.data(target + 254 * ncomps + c);
        auto *t_255 = buffer.data(target + 255 * ncomps + c);
        auto *t_256 = buffer.data(target + 256 * ncomps + c);
        auto *t_257 = buffer.data(target + 257 * ncomps + c);
        auto *t_258 = buffer.data(target + 258 * ncomps + c);
        auto *t_259 = buffer.data(target + 259 * ncomps + c);
        auto *t_260 = buffer.data(target + 260 * ncomps + c);
        auto *t_261 = buffer.data(target + 261 * ncomps + c);
        auto *t_262 = buffer.data(target + 262 * ncomps + c);
        auto *t_263 = buffer.data(target + 263 * ncomps + c);
        auto *t_264 = buffer.data(target + 264 * ncomps + c);
        auto *t_265 = buffer.data(target + 265 * ncomps + c);
        auto *t_266 = buffer.data(target + 266 * ncomps + c);
        auto *t_267 = buffer.data(target + 267 * ncomps + c);
        auto *t_268 = buffer.data(target + 268 * ncomps + c);
        auto *t_269 = buffer.data(target + 269 * ncomps + c);
        auto *t_270 = buffer.data(target + 270 * ncomps + c);
        auto *t_271 = buffer.data(target + 271 * ncomps + c);
        auto *t_272 = buffer.data(target + 272 * ncomps + c);
        auto *t_273 = buffer.data(target + 273 * ncomps + c);
        auto *t_274 = buffer.data(target + 274 * ncomps + c);
        auto *t_275 = buffer.data(target + 275 * ncomps + c);
        auto *t_276 = buffer.data(target + 276 * ncomps + c);
        auto *t_277 = buffer.data(target + 277 * ncomps + c);
        auto *t_278 = buffer.data(target + 278 * ncomps + c);
        auto *t_279 = buffer.data(target + 279 * ncomps + c);
        auto *t_280 = buffer.data(target + 280 * ncomps + c);
        auto *t_281 = buffer.data(target + 281 * ncomps + c);
        auto *t_282 = buffer.data(target + 282 * ncomps + c);
        auto *t_283 = buffer.data(target + 283 * ncomps + c);
        auto *t_284 = buffer.data(target + 284 * ncomps + c);
        auto *t_285 = buffer.data(target + 285 * ncomps + c);
        auto *t_286 = buffer.data(target + 286 * ncomps + c);
        auto *t_287 = buffer.data(target + 287 * ncomps + c);
        auto *t_288 = buffer.data(target + 288 * ncomps + c);
        auto *t_289 = buffer.data(target + 289 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ih_110 = buffer.data(ih + 110 * ncomps + c);
        const auto *ih_111 = buffer.data(ih + 111 * ncomps + c);
        const auto *ih_112 = buffer.data(ih + 112 * ncomps + c);
        const auto *ih_113 = buffer.data(ih + 113 * ncomps + c);
        const auto *ih_114 = buffer.data(ih + 114 * ncomps + c);
        const auto *ih_115 = buffer.data(ih + 115 * ncomps + c);
        const auto *ih_116 = buffer.data(ih + 116 * ncomps + c);
        const auto *ih_117 = buffer.data(ih + 117 * ncomps + c);
        const auto *ih_118 = buffer.data(ih + 118 * ncomps + c);
        const auto *ih_119 = buffer.data(ih + 119 * ncomps + c);
        const auto *ih_120 = buffer.data(ih + 120 * ncomps + c);
        const auto *ih_121 = buffer.data(ih + 121 * ncomps + c);
        const auto *ih_122 = buffer.data(ih + 122 * ncomps + c);
        const auto *ih_123 = buffer.data(ih + 123 * ncomps + c);
        const auto *ih_124 = buffer.data(ih + 124 * ncomps + c);
        const auto *ih_125 = buffer.data(ih + 125 * ncomps + c);
        const auto *ih_126 = buffer.data(ih + 126 * ncomps + c);
        const auto *ih_127 = buffer.data(ih + 127 * ncomps + c);
        const auto *ih_128 = buffer.data(ih + 128 * ncomps + c);
        const auto *ih_129 = buffer.data(ih + 129 * ncomps + c);
        const auto *ih_130 = buffer.data(ih + 130 * ncomps + c);
        const auto *ih_131 = buffer.data(ih + 131 * ncomps + c);
        const auto *ih_132 = buffer.data(ih + 132 * ncomps + c);
        const auto *ih_133 = buffer.data(ih + 133 * ncomps + c);
        const auto *ih_134 = buffer.data(ih + 134 * ncomps + c);
        const auto *ih_135 = buffer.data(ih + 135 * ncomps + c);
        const auto *ih_136 = buffer.data(ih + 136 * ncomps + c);
        const auto *ih_137 = buffer.data(ih + 137 * ncomps + c);
        const auto *ih_138 = buffer.data(ih + 138 * ncomps + c);
        const auto *ih_139 = buffer.data(ih + 139 * ncomps + c);
        const auto *ih_140 = buffer.data(ih + 140 * ncomps + c);
        const auto *ih_141 = buffer.data(ih + 141 * ncomps + c);
        const auto *ih_142 = buffer.data(ih + 142 * ncomps + c);
        const auto *ih_143 = buffer.data(ih + 143 * ncomps + c);
        const auto *ih_144 = buffer.data(ih + 144 * ncomps + c);
        const auto *ih_145 = buffer.data(ih + 145 * ncomps + c);
        const auto *ih_146 = buffer.data(ih + 146 * ncomps + c);
        const auto *ih_147 = buffer.data(ih + 147 * ncomps + c);
        const auto *ih_148 = buffer.data(ih + 148 * ncomps + c);
        const auto *ih_149 = buffer.data(ih + 149 * ncomps + c);
        const auto *ih_150 = buffer.data(ih + 150 * ncomps + c);
        const auto *ih_151 = buffer.data(ih + 151 * ncomps + c);
        const auto *ih_152 = buffer.data(ih + 152 * ncomps + c);
        const auto *ih_153 = buffer.data(ih + 153 * ncomps + c);
        const auto *ih_154 = buffer.data(ih + 154 * ncomps + c);
        const auto *ih_155 = buffer.data(ih + 155 * ncomps + c);
        const auto *ih_156 = buffer.data(ih + 156 * ncomps + c);
        const auto *ih_157 = buffer.data(ih + 157 * ncomps + c);
        const auto *ih_158 = buffer.data(ih + 158 * ncomps + c);
        const auto *ih_159 = buffer.data(ih + 159 * ncomps + c);
        const auto *ih_160 = buffer.data(ih + 160 * ncomps + c);
        const auto *ih_161 = buffer.data(ih + 161 * ncomps + c);
        const auto *ih_162 = buffer.data(ih + 162 * ncomps + c);
        const auto *ih_163 = buffer.data(ih + 163 * ncomps + c);
        const auto *ih_164 = buffer.data(ih + 164 * ncomps + c);
        const auto *ih_165 = buffer.data(ih + 165 * ncomps + c);
        const auto *ih_166 = buffer.data(ih + 166 * ncomps + c);
        const auto *ih_167 = buffer.data(ih + 167 * ncomps + c);
        const auto *ih_168 = buffer.data(ih + 168 * ncomps + c);
        const auto *ih_169 = buffer.data(ih + 169 * ncomps + c);
        const auto *ih_170 = buffer.data(ih + 170 * ncomps + c);
        const auto *ih_171 = buffer.data(ih + 171 * ncomps + c);
        const auto *ih_172 = buffer.data(ih + 172 * ncomps + c);
        const auto *ih_173 = buffer.data(ih + 173 * ncomps + c);
        const auto *ih_174 = buffer.data(ih + 174 * ncomps + c);
        const auto *ih_175 = buffer.data(ih + 175 * ncomps + c);
        const auto *ih_176 = buffer.data(ih + 176 * ncomps + c);
        const auto *ih_177 = buffer.data(ih + 177 * ncomps + c);
        const auto *ih_178 = buffer.data(ih + 178 * ncomps + c);
        const auto *ih_179 = buffer.data(ih + 179 * ncomps + c);
        const auto *ih_180 = buffer.data(ih + 180 * ncomps + c);
        const auto *ih_181 = buffer.data(ih + 181 * ncomps + c);
        const auto *ih_182 = buffer.data(ih + 182 * ncomps + c);
        const auto *ih_183 = buffer.data(ih + 183 * ncomps + c);
        const auto *ih_184 = buffer.data(ih + 184 * ncomps + c);
        const auto *ih_185 = buffer.data(ih + 185 * ncomps + c);
        const auto *ih_186 = buffer.data(ih + 186 * ncomps + c);
        const auto *ih_187 = buffer.data(ih + 187 * ncomps + c);
        const auto *ih_188 = buffer.data(ih + 188 * ncomps + c);
        const auto *ih_189 = buffer.data(ih + 189 * ncomps + c);
        const auto *ih_190 = buffer.data(ih + 190 * ncomps + c);
        const auto *ih_191 = buffer.data(ih + 191 * ncomps + c);
        const auto *ih_192 = buffer.data(ih + 192 * ncomps + c);
        const auto *ih_193 = buffer.data(ih + 193 * ncomps + c);
        const auto *ih_194 = buffer.data(ih + 194 * ncomps + c);
        const auto *ih_195 = buffer.data(ih + 195 * ncomps + c);
        const auto *ih_196 = buffer.data(ih + 196 * ncomps + c);
        const auto *ih_197 = buffer.data(ih + 197 * ncomps + c);
        const auto *ih_198 = buffer.data(ih + 198 * ncomps + c);
        const auto *ih_199 = buffer.data(ih + 199 * ncomps + c);
        const auto *ih_200 = buffer.data(ih + 200 * ncomps + c);
        const auto *ih_201 = buffer.data(ih + 201 * ncomps + c);
        const auto *ih_202 = buffer.data(ih + 202 * ncomps + c);
        const auto *ih_203 = buffer.data(ih + 203 * ncomps + c);
        const auto *ih_204 = buffer.data(ih + 204 * ncomps + c);
        const auto *ih_205 = buffer.data(ih + 205 * ncomps + c);
        const auto *ih_206 = buffer.data(ih + 206 * ncomps + c);
        const auto *ih_207 = buffer.data(ih + 207 * ncomps + c);
        const auto *ih_208 = buffer.data(ih + 208 * ncomps + c);
        const auto *ih_209 = buffer.data(ih + 209 * ncomps + c);
        const auto *ih_210 = buffer.data(ih + 210 * ncomps + c);
        const auto *ih_211 = buffer.data(ih + 211 * ncomps + c);
        const auto *ih_212 = buffer.data(ih + 212 * ncomps + c);
        const auto *ih_213 = buffer.data(ih + 213 * ncomps + c);
        const auto *ih_214 = buffer.data(ih + 214 * ncomps + c);
        const auto *ih_215 = buffer.data(ih + 215 * ncomps + c);
        const auto *ih_216 = buffer.data(ih + 216 * ncomps + c);
        const auto *ih_217 = buffer.data(ih + 217 * ncomps + c);
        const auto *ih_218 = buffer.data(ih + 218 * ncomps + c);
        const auto *ih_219 = buffer.data(ih + 219 * ncomps + c);

        const auto *kh_110 = buffer.data(kh + 110 * ncomps + c);
        const auto *kh_111 = buffer.data(kh + 111 * ncomps + c);
        const auto *kh_112 = buffer.data(kh + 112 * ncomps + c);
        const auto *kh_113 = buffer.data(kh + 113 * ncomps + c);
        const auto *kh_114 = buffer.data(kh + 114 * ncomps + c);
        const auto *kh_115 = buffer.data(kh + 115 * ncomps + c);
        const auto *kh_116 = buffer.data(kh + 116 * ncomps + c);
        const auto *kh_117 = buffer.data(kh + 117 * ncomps + c);
        const auto *kh_118 = buffer.data(kh + 118 * ncomps + c);
        const auto *kh_119 = buffer.data(kh + 119 * ncomps + c);
        const auto *kh_120 = buffer.data(kh + 120 * ncomps + c);
        const auto *kh_121 = buffer.data(kh + 121 * ncomps + c);
        const auto *kh_122 = buffer.data(kh + 122 * ncomps + c);
        const auto *kh_123 = buffer.data(kh + 123 * ncomps + c);
        const auto *kh_124 = buffer.data(kh + 124 * ncomps + c);
        const auto *kh_125 = buffer.data(kh + 125 * ncomps + c);
        const auto *kh_126 = buffer.data(kh + 126 * ncomps + c);
        const auto *kh_127 = buffer.data(kh + 127 * ncomps + c);
        const auto *kh_128 = buffer.data(kh + 128 * ncomps + c);
        const auto *kh_129 = buffer.data(kh + 129 * ncomps + c);
        const auto *kh_130 = buffer.data(kh + 130 * ncomps + c);
        const auto *kh_131 = buffer.data(kh + 131 * ncomps + c);
        const auto *kh_132 = buffer.data(kh + 132 * ncomps + c);
        const auto *kh_133 = buffer.data(kh + 133 * ncomps + c);
        const auto *kh_134 = buffer.data(kh + 134 * ncomps + c);
        const auto *kh_135 = buffer.data(kh + 135 * ncomps + c);
        const auto *kh_136 = buffer.data(kh + 136 * ncomps + c);
        const auto *kh_137 = buffer.data(kh + 137 * ncomps + c);
        const auto *kh_138 = buffer.data(kh + 138 * ncomps + c);
        const auto *kh_139 = buffer.data(kh + 139 * ncomps + c);
        const auto *kh_140 = buffer.data(kh + 140 * ncomps + c);
        const auto *kh_141 = buffer.data(kh + 141 * ncomps + c);
        const auto *kh_142 = buffer.data(kh + 142 * ncomps + c);
        const auto *kh_143 = buffer.data(kh + 143 * ncomps + c);
        const auto *kh_144 = buffer.data(kh + 144 * ncomps + c);
        const auto *kh_145 = buffer.data(kh + 145 * ncomps + c);
        const auto *kh_146 = buffer.data(kh + 146 * ncomps + c);
        const auto *kh_147 = buffer.data(kh + 147 * ncomps + c);
        const auto *kh_148 = buffer.data(kh + 148 * ncomps + c);
        const auto *kh_149 = buffer.data(kh + 149 * ncomps + c);
        const auto *kh_150 = buffer.data(kh + 150 * ncomps + c);
        const auto *kh_151 = buffer.data(kh + 151 * ncomps + c);
        const auto *kh_152 = buffer.data(kh + 152 * ncomps + c);
        const auto *kh_153 = buffer.data(kh + 153 * ncomps + c);
        const auto *kh_154 = buffer.data(kh + 154 * ncomps + c);
        const auto *kh_155 = buffer.data(kh + 155 * ncomps + c);
        const auto *kh_156 = buffer.data(kh + 156 * ncomps + c);
        const auto *kh_157 = buffer.data(kh + 157 * ncomps + c);
        const auto *kh_158 = buffer.data(kh + 158 * ncomps + c);
        const auto *kh_159 = buffer.data(kh + 159 * ncomps + c);
        const auto *kh_160 = buffer.data(kh + 160 * ncomps + c);
        const auto *kh_161 = buffer.data(kh + 161 * ncomps + c);
        const auto *kh_162 = buffer.data(kh + 162 * ncomps + c);
        const auto *kh_163 = buffer.data(kh + 163 * ncomps + c);
        const auto *kh_164 = buffer.data(kh + 164 * ncomps + c);
        const auto *kh_165 = buffer.data(kh + 165 * ncomps + c);
        const auto *kh_166 = buffer.data(kh + 166 * ncomps + c);
        const auto *kh_167 = buffer.data(kh + 167 * ncomps + c);
        const auto *kh_168 = buffer.data(kh + 168 * ncomps + c);
        const auto *kh_169 = buffer.data(kh + 169 * ncomps + c);
        const auto *kh_170 = buffer.data(kh + 170 * ncomps + c);
        const auto *kh_171 = buffer.data(kh + 171 * ncomps + c);
        const auto *kh_172 = buffer.data(kh + 172 * ncomps + c);
        const auto *kh_173 = buffer.data(kh + 173 * ncomps + c);
        const auto *kh_174 = buffer.data(kh + 174 * ncomps + c);
        const auto *kh_175 = buffer.data(kh + 175 * ncomps + c);
        const auto *kh_176 = buffer.data(kh + 176 * ncomps + c);
        const auto *kh_177 = buffer.data(kh + 177 * ncomps + c);
        const auto *kh_178 = buffer.data(kh + 178 * ncomps + c);
        const auto *kh_179 = buffer.data(kh + 179 * ncomps + c);
        const auto *kh_180 = buffer.data(kh + 180 * ncomps + c);
        const auto *kh_181 = buffer.data(kh + 181 * ncomps + c);
        const auto *kh_182 = buffer.data(kh + 182 * ncomps + c);
        const auto *kh_183 = buffer.data(kh + 183 * ncomps + c);
        const auto *kh_184 = buffer.data(kh + 184 * ncomps + c);
        const auto *kh_185 = buffer.data(kh + 185 * ncomps + c);
        const auto *kh_186 = buffer.data(kh + 186 * ncomps + c);
        const auto *kh_187 = buffer.data(kh + 187 * ncomps + c);
        const auto *kh_188 = buffer.data(kh + 188 * ncomps + c);
        const auto *kh_189 = buffer.data(kh + 189 * ncomps + c);
        const auto *kh_190 = buffer.data(kh + 190 * ncomps + c);
        const auto *kh_191 = buffer.data(kh + 191 * ncomps + c);
        const auto *kh_192 = buffer.data(kh + 192 * ncomps + c);
        const auto *kh_193 = buffer.data(kh + 193 * ncomps + c);
        const auto *kh_194 = buffer.data(kh + 194 * ncomps + c);
        const auto *kh_195 = buffer.data(kh + 195 * ncomps + c);
        const auto *kh_196 = buffer.data(kh + 196 * ncomps + c);
        const auto *kh_197 = buffer.data(kh + 197 * ncomps + c);
        const auto *kh_198 = buffer.data(kh + 198 * ncomps + c);
        const auto *kh_199 = buffer.data(kh + 199 * ncomps + c);
        const auto *kh_200 = buffer.data(kh + 200 * ncomps + c);
        const auto *kh_201 = buffer.data(kh + 201 * ncomps + c);
        const auto *kh_202 = buffer.data(kh + 202 * ncomps + c);
        const auto *kh_203 = buffer.data(kh + 203 * ncomps + c);
        const auto *kh_204 = buffer.data(kh + 204 * ncomps + c);
        const auto *kh_205 = buffer.data(kh + 205 * ncomps + c);
        const auto *kh_206 = buffer.data(kh + 206 * ncomps + c);
        const auto *kh_207 = buffer.data(kh + 207 * ncomps + c);
        const auto *kh_208 = buffer.data(kh + 208 * ncomps + c);
        const auto *kh_209 = buffer.data(kh + 209 * ncomps + c);
        const auto *kh_210 = buffer.data(kh + 210 * ncomps + c);
        const auto *kh_211 = buffer.data(kh + 211 * ncomps + c);
        const auto *kh_212 = buffer.data(kh + 212 * ncomps + c);
        const auto *kh_213 = buffer.data(kh + 213 * ncomps + c);
        const auto *kh_214 = buffer.data(kh + 214 * ncomps + c);
        const auto *kh_215 = buffer.data(kh + 215 * ncomps + c);
        const auto *kh_216 = buffer.data(kh + 216 * ncomps + c);
        const auto *kh_217 = buffer.data(kh + 217 * ncomps + c);
        const auto *kh_218 = buffer.data(kh + 218 * ncomps + c);
        const auto *kh_219 = buffer.data(kh + 219 * ncomps + c);
        const auto *kh_225 = buffer.data(kh + 225 * ncomps + c);
        const auto *kh_226 = buffer.data(kh + 226 * ncomps + c);
        const auto *kh_227 = buffer.data(kh + 227 * ncomps + c);
        const auto *kh_228 = buffer.data(kh + 228 * ncomps + c);
        const auto *kh_229 = buffer.data(kh + 229 * ncomps + c);
        const auto *kh_230 = buffer.data(kh + 230 * ncomps + c);
        const auto *kh_246 = buffer.data(kh + 246 * ncomps + c);
        const auto *kh_247 = buffer.data(kh + 247 * ncomps + c);
        const auto *kh_248 = buffer.data(kh + 248 * ncomps + c);
        const auto *kh_249 = buffer.data(kh + 249 * ncomps + c);
        const auto *kh_250 = buffer.data(kh + 250 * ncomps + c);
        const auto *kh_251 = buffer.data(kh + 251 * ncomps + c);
        const auto *kh_267 = buffer.data(kh + 267 * ncomps + c);
        const auto *kh_268 = buffer.data(kh + 268 * ncomps + c);
        const auto *kh_269 = buffer.data(kh + 269 * ncomps + c);
        const auto *kh_270 = buffer.data(kh + 270 * ncomps + c);
        const auto *kh_271 = buffer.data(kh + 271 * ncomps + c);
        const auto *kh_272 = buffer.data(kh + 272 * ncomps + c);
        const auto *kh_288 = buffer.data(kh + 288 * ncomps + c);
        const auto *kh_289 = buffer.data(kh + 289 * ncomps + c);
        const auto *kh_290 = buffer.data(kh + 290 * ncomps + c);
        const auto *kh_291 = buffer.data(kh + 291 * ncomps + c);
        const auto *kh_292 = buffer.data(kh + 292 * ncomps + c);
        const auto *kh_293 = buffer.data(kh + 293 * ncomps + c);
        const auto *kh_314 = buffer.data(kh + 314 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ih_110, ih_111, ih_112, \
                         ih_113, ih_114, kh_110, kh_111, kh_112, kh_113, \
                         kh_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * ih_110[k]
                       + kh_110[k];

            t_146[k] = ab_x[k] * ih_111[k]
                       + kh_111[k];

            t_147[k] = ab_x[k] * ih_112[k]
                       + kh_112[k];

            t_148[k] = ab_x[k] * ih_113[k]
                       + kh_113[k];

            t_149[k] = ab_x[k] * ih_114[k]
                       + kh_114[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, ih_115, ih_116, ih_117, \
                         ih_118, ih_119, kh_115, kh_116, kh_117, kh_118, \
                         kh_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ih_115[k]
                       + kh_115[k];

            t_151[k] = ab_x[k] * ih_116[k]
                       + kh_116[k];

            t_152[k] = ab_x[k] * ih_117[k]
                       + kh_117[k];

            t_153[k] = ab_x[k] * ih_118[k]
                       + kh_118[k];

            t_154[k] = ab_x[k] * ih_119[k]
                       + kh_119[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ih_120, ih_121, ih_122, \
                         ih_123, ih_124, kh_120, kh_121, kh_122, kh_123, \
                         kh_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * ih_120[k]
                       + kh_120[k];

            t_156[k] = ab_x[k] * ih_121[k]
                       + kh_121[k];

            t_157[k] = ab_x[k] * ih_122[k]
                       + kh_122[k];

            t_158[k] = ab_x[k] * ih_123[k]
                       + kh_123[k];

            t_159[k] = ab_x[k] * ih_124[k]
                       + kh_124[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, ab_x, ab_y, ih_120, ih_121, ih_122, \
                         ih_125, kh_125, kh_183, kh_184, kh_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * ih_125[k]
                       + kh_125[k];

            t_161[k] = ab_y[k] * ih_120[k]
                       + kh_183[k];

            t_162[k] = ab_y[k] * ih_121[k]
                       + kh_184[k];

            t_163[k] = ab_y[k] * ih_122[k]
                       + kh_185[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, ab_y, ab_z, ih_123, ih_124, ih_125, \
                         kh_186, kh_187, kh_188, kh_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = ab_y[k] * ih_123[k]
                       + kh_186[k];

            t_165[k] = ab_y[k] * ih_124[k]
                       + kh_187[k];

            t_166[k] = ab_y[k] * ih_125[k]
                       + kh_188[k];

            t_167[k] = ab_z[k] * ih_125[k]
                       + kh_209[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_x, ih_126, ih_127, ih_128, \
                         ih_129, ih_130, kh_126, kh_127, kh_128, kh_129, \
                         kh_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = ab_x[k] * ih_126[k]
                       + kh_126[k];

            t_169[k] = ab_x[k] * ih_127[k]
                       + kh_127[k];

            t_170[k] = ab_x[k] * ih_128[k]
                       + kh_128[k];

            t_171[k] = ab_x[k] * ih_129[k]
                       + kh_129[k];

            t_172[k] = ab_x[k] * ih_130[k]
                       + kh_130[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_x, ih_131, ih_132, ih_133, \
                         ih_134, ih_135, kh_131, kh_132, kh_133, kh_134, \
                         kh_135 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_x[k] * ih_131[k]
                       + kh_131[k];

            t_174[k] = ab_x[k] * ih_132[k]
                       + kh_132[k];

            t_175[k] = ab_x[k] * ih_133[k]
                       + kh_133[k];

            t_176[k] = ab_x[k] * ih_134[k]
                       + kh_134[k];

            t_177[k] = ab_x[k] * ih_135[k]
                       + kh_135[k];
        }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_x, ih_136, ih_137, ih_138, \
                         ih_139, ih_140, kh_136, kh_137, kh_138, kh_139, \
                         kh_140 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_178[k] = ab_x[k] * ih_136[k]
                       + kh_136[k];

            t_179[k] = ab_x[k] * ih_137[k]
                       + kh_137[k];

            t_180[k] = ab_x[k] * ih_138[k]
                       + kh_138[k];

            t_181[k] = ab_x[k] * ih_139[k]
                       + kh_139[k];

            t_182[k] = ab_x[k] * ih_140[k]
                       + kh_140[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_x, ih_141, ih_142, ih_143, \
                         ih_144, ih_145, kh_141, kh_142, kh_143, kh_144, \
                         kh_145 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = ab_x[k] * ih_141[k]
                       + kh_141[k];

            t_184[k] = ab_x[k] * ih_142[k]
                       + kh_142[k];

            t_185[k] = ab_x[k] * ih_143[k]
                       + kh_143[k];

            t_186[k] = ab_x[k] * ih_144[k]
                       + kh_144[k];

            t_187[k] = ab_x[k] * ih_145[k]
                       + kh_145[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, ab_x, ab_y, ih_141, ih_142, ih_143, \
                         ih_146, kh_146, kh_225, kh_226, kh_227 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = ab_x[k] * ih_146[k]
                       + kh_146[k];

            t_189[k] = ab_y[k] * ih_141[k]
                       + kh_225[k];

            t_190[k] = ab_y[k] * ih_142[k]
                       + kh_226[k];

            t_191[k] = ab_y[k] * ih_143[k]
                       + kh_227[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, ab_y, ab_z, ih_144, ih_145, ih_146, \
                         kh_228, kh_229, kh_230, kh_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = ab_y[k] * ih_144[k]
                       + kh_228[k];

            t_193[k] = ab_y[k] * ih_145[k]
                       + kh_229[k];

            t_194[k] = ab_y[k] * ih_146[k]
                       + kh_230[k];

            t_195[k] = ab_z[k] * ih_146[k]
                       + kh_251[k];
        }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, ab_x, ih_147, ih_148, ih_149, \
                         ih_150, ih_151, kh_147, kh_148, kh_149, kh_150, \
                         kh_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_196[k] = ab_x[k] * ih_147[k]
                       + kh_147[k];

            t_197[k] = ab_x[k] * ih_148[k]
                       + kh_148[k];

            t_198[k] = ab_x[k] * ih_149[k]
                       + kh_149[k];

            t_199[k] = ab_x[k] * ih_150[k]
                       + kh_150[k];

            t_200[k] = ab_x[k] * ih_151[k]
                       + kh_151[k];
        }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, ab_x, ih_152, ih_153, ih_154, \
                         ih_155, ih_156, kh_152, kh_153, kh_154, kh_155, \
                         kh_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_201[k] = ab_x[k] * ih_152[k]
                       + kh_152[k];

            t_202[k] = ab_x[k] * ih_153[k]
                       + kh_153[k];

            t_203[k] = ab_x[k] * ih_154[k]
                       + kh_154[k];

            t_204[k] = ab_x[k] * ih_155[k]
                       + kh_155[k];

            t_205[k] = ab_x[k] * ih_156[k]
                       + kh_156[k];
        }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, ab_x, ih_157, ih_158, ih_159, \
                         ih_160, ih_161, kh_157, kh_158, kh_159, kh_160, \
                         kh_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_206[k] = ab_x[k] * ih_157[k]
                       + kh_157[k];

            t_207[k] = ab_x[k] * ih_158[k]
                       + kh_158[k];

            t_208[k] = ab_x[k] * ih_159[k]
                       + kh_159[k];

            t_209[k] = ab_x[k] * ih_160[k]
                       + kh_160[k];

            t_210[k] = ab_x[k] * ih_161[k]
                       + kh_161[k];
        }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ab_x, ih_162, ih_163, ih_164, \
                         ih_165, ih_166, kh_162, kh_163, kh_164, kh_165, \
                         kh_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_211[k] = ab_x[k] * ih_162[k]
                       + kh_162[k];

            t_212[k] = ab_x[k] * ih_163[k]
                       + kh_163[k];

            t_213[k] = ab_x[k] * ih_164[k]
                       + kh_164[k];

            t_214[k] = ab_x[k] * ih_165[k]
                       + kh_165[k];

            t_215[k] = ab_x[k] * ih_166[k]
                       + kh_166[k];
        }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, ab_x, ab_y, ih_162, ih_163, ih_164, \
                         ih_167, kh_167, kh_246, kh_247, kh_248 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_216[k] = ab_x[k] * ih_167[k]
                       + kh_167[k];

            t_217[k] = ab_y[k] * ih_162[k]
                       + kh_246[k];

            t_218[k] = ab_y[k] * ih_163[k]
                       + kh_247[k];

            t_219[k] = ab_y[k] * ih_164[k]
                       + kh_248[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, ab_y, ab_z, ih_165, ih_166, ih_167, \
                         kh_249, kh_250, kh_251, kh_272 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * ih_165[k]
                       + kh_249[k];

            t_221[k] = ab_y[k] * ih_166[k]
                       + kh_250[k];

            t_222[k] = ab_y[k] * ih_167[k]
                       + kh_251[k];

            t_223[k] = ab_z[k] * ih_167[k]
                       + kh_272[k];
        }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, ab_x, ih_168, ih_169, ih_170, \
                         ih_171, ih_172, kh_168, kh_169, kh_170, kh_171, \
                         kh_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_224[k] = ab_x[k] * ih_168[k]
                       + kh_168[k];

            t_225[k] = ab_x[k] * ih_169[k]
                       + kh_169[k];

            t_226[k] = ab_x[k] * ih_170[k]
                       + kh_170[k];

            t_227[k] = ab_x[k] * ih_171[k]
                       + kh_171[k];

            t_228[k] = ab_x[k] * ih_172[k]
                       + kh_172[k];
        }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, ab_x, ih_173, ih_174, ih_175, \
                         ih_176, ih_177, kh_173, kh_174, kh_175, kh_176, \
                         kh_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_229[k] = ab_x[k] * ih_173[k]
                       + kh_173[k];

            t_230[k] = ab_x[k] * ih_174[k]
                       + kh_174[k];

            t_231[k] = ab_x[k] * ih_175[k]
                       + kh_175[k];

            t_232[k] = ab_x[k] * ih_176[k]
                       + kh_176[k];

            t_233[k] = ab_x[k] * ih_177[k]
                       + kh_177[k];
        }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, ih_178, ih_179, ih_180, \
                         ih_181, ih_182, kh_178, kh_179, kh_180, kh_181, \
                         kh_182 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_234[k] = ab_x[k] * ih_178[k]
                       + kh_178[k];

            t_235[k] = ab_x[k] * ih_179[k]
                       + kh_179[k];

            t_236[k] = ab_x[k] * ih_180[k]
                       + kh_180[k];

            t_237[k] = ab_x[k] * ih_181[k]
                       + kh_181[k];

            t_238[k] = ab_x[k] * ih_182[k]
                       + kh_182[k];
        }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_x, ih_183, ih_184, ih_185, \
                         ih_186, ih_187, kh_183, kh_184, kh_185, kh_186, \
                         kh_187 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_239[k] = ab_x[k] * ih_183[k]
                       + kh_183[k];

            t_240[k] = ab_x[k] * ih_184[k]
                       + kh_184[k];

            t_241[k] = ab_x[k] * ih_185[k]
                       + kh_185[k];

            t_242[k] = ab_x[k] * ih_186[k]
                       + kh_186[k];

            t_243[k] = ab_x[k] * ih_187[k]
                       + kh_187[k];
        }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, ab_x, ab_y, ih_183, ih_184, ih_185, \
                         ih_188, kh_188, kh_267, kh_268, kh_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_244[k] = ab_x[k] * ih_188[k]
                       + kh_188[k];

            t_245[k] = ab_y[k] * ih_183[k]
                       + kh_267[k];

            t_246[k] = ab_y[k] * ih_184[k]
                       + kh_268[k];

            t_247[k] = ab_y[k] * ih_185[k]
                       + kh_269[k];
        }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, ab_y, ab_z, ih_186, ih_187, ih_188, \
                         kh_270, kh_271, kh_272, kh_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_248[k] = ab_y[k] * ih_186[k]
                       + kh_270[k];

            t_249[k] = ab_y[k] * ih_187[k]
                       + kh_271[k];

            t_250[k] = ab_y[k] * ih_188[k]
                       + kh_272[k];

            t_251[k] = ab_z[k] * ih_188[k]
                       + kh_293[k];
        }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, ih_189, ih_190, ih_191, \
                         ih_192, ih_193, kh_189, kh_190, kh_191, kh_192, \
                         kh_193 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_252[k] = ab_x[k] * ih_189[k]
                       + kh_189[k];

            t_253[k] = ab_x[k] * ih_190[k]
                       + kh_190[k];

            t_254[k] = ab_x[k] * ih_191[k]
                       + kh_191[k];

            t_255[k] = ab_x[k] * ih_192[k]
                       + kh_192[k];

            t_256[k] = ab_x[k] * ih_193[k]
                       + kh_193[k];
        }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ab_x, ih_194, ih_195, ih_196, \
                         ih_197, ih_198, kh_194, kh_195, kh_196, kh_197, \
                         kh_198 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_257[k] = ab_x[k] * ih_194[k]
                       + kh_194[k];

            t_258[k] = ab_x[k] * ih_195[k]
                       + kh_195[k];

            t_259[k] = ab_x[k] * ih_196[k]
                       + kh_196[k];

            t_260[k] = ab_x[k] * ih_197[k]
                       + kh_197[k];

            t_261[k] = ab_x[k] * ih_198[k]
                       + kh_198[k];
        }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ab_x, ih_199, ih_200, ih_201, \
                         ih_202, ih_203, kh_199, kh_200, kh_201, kh_202, \
                         kh_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_262[k] = ab_x[k] * ih_199[k]
                       + kh_199[k];

            t_263[k] = ab_x[k] * ih_200[k]
                       + kh_200[k];

            t_264[k] = ab_x[k] * ih_201[k]
                       + kh_201[k];

            t_265[k] = ab_x[k] * ih_202[k]
                       + kh_202[k];

            t_266[k] = ab_x[k] * ih_203[k]
                       + kh_203[k];
        }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ab_x, ih_204, ih_205, ih_206, \
                         ih_207, ih_208, kh_204, kh_205, kh_206, kh_207, \
                         kh_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_267[k] = ab_x[k] * ih_204[k]
                       + kh_204[k];

            t_268[k] = ab_x[k] * ih_205[k]
                       + kh_205[k];

            t_269[k] = ab_x[k] * ih_206[k]
                       + kh_206[k];

            t_270[k] = ab_x[k] * ih_207[k]
                       + kh_207[k];

            t_271[k] = ab_x[k] * ih_208[k]
                       + kh_208[k];
        }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, ab_x, ab_y, ih_204, ih_205, ih_206, \
                         ih_209, kh_209, kh_288, kh_289, kh_290 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_272[k] = ab_x[k] * ih_209[k]
                       + kh_209[k];

            t_273[k] = ab_y[k] * ih_204[k]
                       + kh_288[k];

            t_274[k] = ab_y[k] * ih_205[k]
                       + kh_289[k];

            t_275[k] = ab_y[k] * ih_206[k]
                       + kh_290[k];
        }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, ab_y, ab_z, ih_207, ih_208, ih_209, \
                         kh_291, kh_292, kh_293, kh_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_276[k] = ab_y[k] * ih_207[k]
                       + kh_291[k];

            t_277[k] = ab_y[k] * ih_208[k]
                       + kh_292[k];

            t_278[k] = ab_y[k] * ih_209[k]
                       + kh_293[k];

            t_279[k] = ab_z[k] * ih_209[k]
                       + kh_314[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, ih_210, ih_211, ih_212, \
                         ih_213, ih_214, kh_210, kh_211, kh_212, kh_213, \
                         kh_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_x[k] * ih_210[k]
                       + kh_210[k];

            t_281[k] = ab_x[k] * ih_211[k]
                       + kh_211[k];

            t_282[k] = ab_x[k] * ih_212[k]
                       + kh_212[k];

            t_283[k] = ab_x[k] * ih_213[k]
                       + kh_213[k];

            t_284[k] = ab_x[k] * ih_214[k]
                       + kh_214[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, ih_215, ih_216, ih_217, \
                         ih_218, ih_219, kh_215, kh_216, kh_217, kh_218, \
                         kh_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * ih_215[k]
                       + kh_215[k];

            t_286[k] = ab_x[k] * ih_216[k]
                       + kh_216[k];

            t_287[k] = ab_x[k] * ih_217[k]
                       + kh_217[k];

            t_288[k] = ab_x[k] * ih_218[k]
                       + kh_218[k];

            t_289[k] = ab_x[k] * ih_219[k]
                       + kh_219[k];
        }
    }
}

static auto
compute_hrr_ii_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ih, const size_t kh,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_290 = buffer.data(target + 290 * ncomps + c);
        auto *t_291 = buffer.data(target + 291 * ncomps + c);
        auto *t_292 = buffer.data(target + 292 * ncomps + c);
        auto *t_293 = buffer.data(target + 293 * ncomps + c);
        auto *t_294 = buffer.data(target + 294 * ncomps + c);
        auto *t_295 = buffer.data(target + 295 * ncomps + c);
        auto *t_296 = buffer.data(target + 296 * ncomps + c);
        auto *t_297 = buffer.data(target + 297 * ncomps + c);
        auto *t_298 = buffer.data(target + 298 * ncomps + c);
        auto *t_299 = buffer.data(target + 299 * ncomps + c);
        auto *t_300 = buffer.data(target + 300 * ncomps + c);
        auto *t_301 = buffer.data(target + 301 * ncomps + c);
        auto *t_302 = buffer.data(target + 302 * ncomps + c);
        auto *t_303 = buffer.data(target + 303 * ncomps + c);
        auto *t_304 = buffer.data(target + 304 * ncomps + c);
        auto *t_305 = buffer.data(target + 305 * ncomps + c);
        auto *t_306 = buffer.data(target + 306 * ncomps + c);
        auto *t_307 = buffer.data(target + 307 * ncomps + c);
        auto *t_308 = buffer.data(target + 308 * ncomps + c);
        auto *t_309 = buffer.data(target + 309 * ncomps + c);
        auto *t_310 = buffer.data(target + 310 * ncomps + c);
        auto *t_311 = buffer.data(target + 311 * ncomps + c);
        auto *t_312 = buffer.data(target + 312 * ncomps + c);
        auto *t_313 = buffer.data(target + 313 * ncomps + c);
        auto *t_314 = buffer.data(target + 314 * ncomps + c);
        auto *t_315 = buffer.data(target + 315 * ncomps + c);
        auto *t_316 = buffer.data(target + 316 * ncomps + c);
        auto *t_317 = buffer.data(target + 317 * ncomps + c);
        auto *t_318 = buffer.data(target + 318 * ncomps + c);
        auto *t_319 = buffer.data(target + 319 * ncomps + c);
        auto *t_320 = buffer.data(target + 320 * ncomps + c);
        auto *t_321 = buffer.data(target + 321 * ncomps + c);
        auto *t_322 = buffer.data(target + 322 * ncomps + c);
        auto *t_323 = buffer.data(target + 323 * ncomps + c);
        auto *t_324 = buffer.data(target + 324 * ncomps + c);
        auto *t_325 = buffer.data(target + 325 * ncomps + c);
        auto *t_326 = buffer.data(target + 326 * ncomps + c);
        auto *t_327 = buffer.data(target + 327 * ncomps + c);
        auto *t_328 = buffer.data(target + 328 * ncomps + c);
        auto *t_329 = buffer.data(target + 329 * ncomps + c);
        auto *t_330 = buffer.data(target + 330 * ncomps + c);
        auto *t_331 = buffer.data(target + 331 * ncomps + c);
        auto *t_332 = buffer.data(target + 332 * ncomps + c);
        auto *t_333 = buffer.data(target + 333 * ncomps + c);
        auto *t_334 = buffer.data(target + 334 * ncomps + c);
        auto *t_335 = buffer.data(target + 335 * ncomps + c);
        auto *t_336 = buffer.data(target + 336 * ncomps + c);
        auto *t_337 = buffer.data(target + 337 * ncomps + c);
        auto *t_338 = buffer.data(target + 338 * ncomps + c);
        auto *t_339 = buffer.data(target + 339 * ncomps + c);
        auto *t_340 = buffer.data(target + 340 * ncomps + c);
        auto *t_341 = buffer.data(target + 341 * ncomps + c);
        auto *t_342 = buffer.data(target + 342 * ncomps + c);
        auto *t_343 = buffer.data(target + 343 * ncomps + c);
        auto *t_344 = buffer.data(target + 344 * ncomps + c);
        auto *t_345 = buffer.data(target + 345 * ncomps + c);
        auto *t_346 = buffer.data(target + 346 * ncomps + c);
        auto *t_347 = buffer.data(target + 347 * ncomps + c);
        auto *t_348 = buffer.data(target + 348 * ncomps + c);
        auto *t_349 = buffer.data(target + 349 * ncomps + c);
        auto *t_350 = buffer.data(target + 350 * ncomps + c);
        auto *t_351 = buffer.data(target + 351 * ncomps + c);
        auto *t_352 = buffer.data(target + 352 * ncomps + c);
        auto *t_353 = buffer.data(target + 353 * ncomps + c);
        auto *t_354 = buffer.data(target + 354 * ncomps + c);
        auto *t_355 = buffer.data(target + 355 * ncomps + c);
        auto *t_356 = buffer.data(target + 356 * ncomps + c);
        auto *t_357 = buffer.data(target + 357 * ncomps + c);
        auto *t_358 = buffer.data(target + 358 * ncomps + c);
        auto *t_359 = buffer.data(target + 359 * ncomps + c);
        auto *t_360 = buffer.data(target + 360 * ncomps + c);
        auto *t_361 = buffer.data(target + 361 * ncomps + c);
        auto *t_362 = buffer.data(target + 362 * ncomps + c);
        auto *t_363 = buffer.data(target + 363 * ncomps + c);
        auto *t_364 = buffer.data(target + 364 * ncomps + c);
        auto *t_365 = buffer.data(target + 365 * ncomps + c);
        auto *t_366 = buffer.data(target + 366 * ncomps + c);
        auto *t_367 = buffer.data(target + 367 * ncomps + c);
        auto *t_368 = buffer.data(target + 368 * ncomps + c);
        auto *t_369 = buffer.data(target + 369 * ncomps + c);
        auto *t_370 = buffer.data(target + 370 * ncomps + c);
        auto *t_371 = buffer.data(target + 371 * ncomps + c);
        auto *t_372 = buffer.data(target + 372 * ncomps + c);
        auto *t_373 = buffer.data(target + 373 * ncomps + c);
        auto *t_374 = buffer.data(target + 374 * ncomps + c);
        auto *t_375 = buffer.data(target + 375 * ncomps + c);
        auto *t_376 = buffer.data(target + 376 * ncomps + c);
        auto *t_377 = buffer.data(target + 377 * ncomps + c);
        auto *t_378 = buffer.data(target + 378 * ncomps + c);
        auto *t_379 = buffer.data(target + 379 * ncomps + c);
        auto *t_380 = buffer.data(target + 380 * ncomps + c);
        auto *t_381 = buffer.data(target + 381 * ncomps + c);
        auto *t_382 = buffer.data(target + 382 * ncomps + c);
        auto *t_383 = buffer.data(target + 383 * ncomps + c);
        auto *t_384 = buffer.data(target + 384 * ncomps + c);
        auto *t_385 = buffer.data(target + 385 * ncomps + c);
        auto *t_386 = buffer.data(target + 386 * ncomps + c);
        auto *t_387 = buffer.data(target + 387 * ncomps + c);
        auto *t_388 = buffer.data(target + 388 * ncomps + c);
        auto *t_389 = buffer.data(target + 389 * ncomps + c);
        auto *t_390 = buffer.data(target + 390 * ncomps + c);
        auto *t_391 = buffer.data(target + 391 * ncomps + c);
        auto *t_392 = buffer.data(target + 392 * ncomps + c);
        auto *t_393 = buffer.data(target + 393 * ncomps + c);
        auto *t_394 = buffer.data(target + 394 * ncomps + c);
        auto *t_395 = buffer.data(target + 395 * ncomps + c);
        auto *t_396 = buffer.data(target + 396 * ncomps + c);
        auto *t_397 = buffer.data(target + 397 * ncomps + c);
        auto *t_398 = buffer.data(target + 398 * ncomps + c);
        auto *t_399 = buffer.data(target + 399 * ncomps + c);
        auto *t_400 = buffer.data(target + 400 * ncomps + c);
        auto *t_401 = buffer.data(target + 401 * ncomps + c);
        auto *t_402 = buffer.data(target + 402 * ncomps + c);
        auto *t_403 = buffer.data(target + 403 * ncomps + c);
        auto *t_404 = buffer.data(target + 404 * ncomps + c);
        auto *t_405 = buffer.data(target + 405 * ncomps + c);
        auto *t_406 = buffer.data(target + 406 * ncomps + c);
        auto *t_407 = buffer.data(target + 407 * ncomps + c);
        auto *t_408 = buffer.data(target + 408 * ncomps + c);
        auto *t_409 = buffer.data(target + 409 * ncomps + c);
        auto *t_410 = buffer.data(target + 410 * ncomps + c);
        auto *t_411 = buffer.data(target + 411 * ncomps + c);
        auto *t_412 = buffer.data(target + 412 * ncomps + c);
        auto *t_413 = buffer.data(target + 413 * ncomps + c);
        auto *t_414 = buffer.data(target + 414 * ncomps + c);
        auto *t_415 = buffer.data(target + 415 * ncomps + c);
        auto *t_416 = buffer.data(target + 416 * ncomps + c);
        auto *t_417 = buffer.data(target + 417 * ncomps + c);
        auto *t_418 = buffer.data(target + 418 * ncomps + c);
        auto *t_419 = buffer.data(target + 419 * ncomps + c);
        auto *t_420 = buffer.data(target + 420 * ncomps + c);
        auto *t_421 = buffer.data(target + 421 * ncomps + c);
        auto *t_422 = buffer.data(target + 422 * ncomps + c);
        auto *t_423 = buffer.data(target + 423 * ncomps + c);
        auto *t_424 = buffer.data(target + 424 * ncomps + c);
        auto *t_425 = buffer.data(target + 425 * ncomps + c);
        auto *t_426 = buffer.data(target + 426 * ncomps + c);
        auto *t_427 = buffer.data(target + 427 * ncomps + c);
        auto *t_428 = buffer.data(target + 428 * ncomps + c);
        auto *t_429 = buffer.data(target + 429 * ncomps + c);
        auto *t_430 = buffer.data(target + 430 * ncomps + c);
        auto *t_431 = buffer.data(target + 431 * ncomps + c);
        auto *t_432 = buffer.data(target + 432 * ncomps + c);
        auto *t_433 = buffer.data(target + 433 * ncomps + c);
        auto *t_434 = buffer.data(target + 434 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ih_220 = buffer.data(ih + 220 * ncomps + c);
        const auto *ih_221 = buffer.data(ih + 221 * ncomps + c);
        const auto *ih_222 = buffer.data(ih + 222 * ncomps + c);
        const auto *ih_223 = buffer.data(ih + 223 * ncomps + c);
        const auto *ih_224 = buffer.data(ih + 224 * ncomps + c);
        const auto *ih_225 = buffer.data(ih + 225 * ncomps + c);
        const auto *ih_226 = buffer.data(ih + 226 * ncomps + c);
        const auto *ih_227 = buffer.data(ih + 227 * ncomps + c);
        const auto *ih_228 = buffer.data(ih + 228 * ncomps + c);
        const auto *ih_229 = buffer.data(ih + 229 * ncomps + c);
        const auto *ih_230 = buffer.data(ih + 230 * ncomps + c);
        const auto *ih_231 = buffer.data(ih + 231 * ncomps + c);
        const auto *ih_232 = buffer.data(ih + 232 * ncomps + c);
        const auto *ih_233 = buffer.data(ih + 233 * ncomps + c);
        const auto *ih_234 = buffer.data(ih + 234 * ncomps + c);
        const auto *ih_235 = buffer.data(ih + 235 * ncomps + c);
        const auto *ih_236 = buffer.data(ih + 236 * ncomps + c);
        const auto *ih_237 = buffer.data(ih + 237 * ncomps + c);
        const auto *ih_238 = buffer.data(ih + 238 * ncomps + c);
        const auto *ih_239 = buffer.data(ih + 239 * ncomps + c);
        const auto *ih_240 = buffer.data(ih + 240 * ncomps + c);
        const auto *ih_241 = buffer.data(ih + 241 * ncomps + c);
        const auto *ih_242 = buffer.data(ih + 242 * ncomps + c);
        const auto *ih_243 = buffer.data(ih + 243 * ncomps + c);
        const auto *ih_244 = buffer.data(ih + 244 * ncomps + c);
        const auto *ih_245 = buffer.data(ih + 245 * ncomps + c);
        const auto *ih_246 = buffer.data(ih + 246 * ncomps + c);
        const auto *ih_247 = buffer.data(ih + 247 * ncomps + c);
        const auto *ih_248 = buffer.data(ih + 248 * ncomps + c);
        const auto *ih_249 = buffer.data(ih + 249 * ncomps + c);
        const auto *ih_250 = buffer.data(ih + 250 * ncomps + c);
        const auto *ih_251 = buffer.data(ih + 251 * ncomps + c);
        const auto *ih_252 = buffer.data(ih + 252 * ncomps + c);
        const auto *ih_253 = buffer.data(ih + 253 * ncomps + c);
        const auto *ih_254 = buffer.data(ih + 254 * ncomps + c);
        const auto *ih_255 = buffer.data(ih + 255 * ncomps + c);
        const auto *ih_256 = buffer.data(ih + 256 * ncomps + c);
        const auto *ih_257 = buffer.data(ih + 257 * ncomps + c);
        const auto *ih_258 = buffer.data(ih + 258 * ncomps + c);
        const auto *ih_259 = buffer.data(ih + 259 * ncomps + c);
        const auto *ih_260 = buffer.data(ih + 260 * ncomps + c);
        const auto *ih_261 = buffer.data(ih + 261 * ncomps + c);
        const auto *ih_262 = buffer.data(ih + 262 * ncomps + c);
        const auto *ih_263 = buffer.data(ih + 263 * ncomps + c);
        const auto *ih_264 = buffer.data(ih + 264 * ncomps + c);
        const auto *ih_265 = buffer.data(ih + 265 * ncomps + c);
        const auto *ih_266 = buffer.data(ih + 266 * ncomps + c);
        const auto *ih_267 = buffer.data(ih + 267 * ncomps + c);
        const auto *ih_268 = buffer.data(ih + 268 * ncomps + c);
        const auto *ih_269 = buffer.data(ih + 269 * ncomps + c);
        const auto *ih_270 = buffer.data(ih + 270 * ncomps + c);
        const auto *ih_271 = buffer.data(ih + 271 * ncomps + c);
        const auto *ih_272 = buffer.data(ih + 272 * ncomps + c);
        const auto *ih_273 = buffer.data(ih + 273 * ncomps + c);
        const auto *ih_274 = buffer.data(ih + 274 * ncomps + c);
        const auto *ih_275 = buffer.data(ih + 275 * ncomps + c);
        const auto *ih_276 = buffer.data(ih + 276 * ncomps + c);
        const auto *ih_277 = buffer.data(ih + 277 * ncomps + c);
        const auto *ih_278 = buffer.data(ih + 278 * ncomps + c);
        const auto *ih_279 = buffer.data(ih + 279 * ncomps + c);
        const auto *ih_280 = buffer.data(ih + 280 * ncomps + c);
        const auto *ih_281 = buffer.data(ih + 281 * ncomps + c);
        const auto *ih_282 = buffer.data(ih + 282 * ncomps + c);
        const auto *ih_283 = buffer.data(ih + 283 * ncomps + c);
        const auto *ih_284 = buffer.data(ih + 284 * ncomps + c);
        const auto *ih_285 = buffer.data(ih + 285 * ncomps + c);
        const auto *ih_286 = buffer.data(ih + 286 * ncomps + c);
        const auto *ih_287 = buffer.data(ih + 287 * ncomps + c);
        const auto *ih_288 = buffer.data(ih + 288 * ncomps + c);
        const auto *ih_289 = buffer.data(ih + 289 * ncomps + c);
        const auto *ih_290 = buffer.data(ih + 290 * ncomps + c);
        const auto *ih_291 = buffer.data(ih + 291 * ncomps + c);
        const auto *ih_292 = buffer.data(ih + 292 * ncomps + c);
        const auto *ih_293 = buffer.data(ih + 293 * ncomps + c);
        const auto *ih_294 = buffer.data(ih + 294 * ncomps + c);
        const auto *ih_295 = buffer.data(ih + 295 * ncomps + c);
        const auto *ih_296 = buffer.data(ih + 296 * ncomps + c);
        const auto *ih_297 = buffer.data(ih + 297 * ncomps + c);
        const auto *ih_298 = buffer.data(ih + 298 * ncomps + c);
        const auto *ih_299 = buffer.data(ih + 299 * ncomps + c);
        const auto *ih_300 = buffer.data(ih + 300 * ncomps + c);
        const auto *ih_301 = buffer.data(ih + 301 * ncomps + c);
        const auto *ih_302 = buffer.data(ih + 302 * ncomps + c);
        const auto *ih_303 = buffer.data(ih + 303 * ncomps + c);
        const auto *ih_304 = buffer.data(ih + 304 * ncomps + c);
        const auto *ih_305 = buffer.data(ih + 305 * ncomps + c);
        const auto *ih_306 = buffer.data(ih + 306 * ncomps + c);
        const auto *ih_307 = buffer.data(ih + 307 * ncomps + c);
        const auto *ih_308 = buffer.data(ih + 308 * ncomps + c);
        const auto *ih_309 = buffer.data(ih + 309 * ncomps + c);
        const auto *ih_310 = buffer.data(ih + 310 * ncomps + c);
        const auto *ih_311 = buffer.data(ih + 311 * ncomps + c);
        const auto *ih_312 = buffer.data(ih + 312 * ncomps + c);
        const auto *ih_313 = buffer.data(ih + 313 * ncomps + c);
        const auto *ih_314 = buffer.data(ih + 314 * ncomps + c);
        const auto *ih_315 = buffer.data(ih + 315 * ncomps + c);
        const auto *ih_316 = buffer.data(ih + 316 * ncomps + c);
        const auto *ih_317 = buffer.data(ih + 317 * ncomps + c);
        const auto *ih_318 = buffer.data(ih + 318 * ncomps + c);
        const auto *ih_319 = buffer.data(ih + 319 * ncomps + c);
        const auto *ih_320 = buffer.data(ih + 320 * ncomps + c);
        const auto *ih_321 = buffer.data(ih + 321 * ncomps + c);
        const auto *ih_322 = buffer.data(ih + 322 * ncomps + c);
        const auto *ih_323 = buffer.data(ih + 323 * ncomps + c);
        const auto *ih_324 = buffer.data(ih + 324 * ncomps + c);
        const auto *ih_325 = buffer.data(ih + 325 * ncomps + c);
        const auto *ih_326 = buffer.data(ih + 326 * ncomps + c);
        const auto *ih_327 = buffer.data(ih + 327 * ncomps + c);
        const auto *ih_328 = buffer.data(ih + 328 * ncomps + c);
        const auto *ih_329 = buffer.data(ih + 329 * ncomps + c);

        const auto *kh_220 = buffer.data(kh + 220 * ncomps + c);
        const auto *kh_221 = buffer.data(kh + 221 * ncomps + c);
        const auto *kh_222 = buffer.data(kh + 222 * ncomps + c);
        const auto *kh_223 = buffer.data(kh + 223 * ncomps + c);
        const auto *kh_224 = buffer.data(kh + 224 * ncomps + c);
        const auto *kh_225 = buffer.data(kh + 225 * ncomps + c);
        const auto *kh_226 = buffer.data(kh + 226 * ncomps + c);
        const auto *kh_227 = buffer.data(kh + 227 * ncomps + c);
        const auto *kh_228 = buffer.data(kh + 228 * ncomps + c);
        const auto *kh_229 = buffer.data(kh + 229 * ncomps + c);
        const auto *kh_230 = buffer.data(kh + 230 * ncomps + c);
        const auto *kh_231 = buffer.data(kh + 231 * ncomps + c);
        const auto *kh_232 = buffer.data(kh + 232 * ncomps + c);
        const auto *kh_233 = buffer.data(kh + 233 * ncomps + c);
        const auto *kh_234 = buffer.data(kh + 234 * ncomps + c);
        const auto *kh_235 = buffer.data(kh + 235 * ncomps + c);
        const auto *kh_236 = buffer.data(kh + 236 * ncomps + c);
        const auto *kh_237 = buffer.data(kh + 237 * ncomps + c);
        const auto *kh_238 = buffer.data(kh + 238 * ncomps + c);
        const auto *kh_239 = buffer.data(kh + 239 * ncomps + c);
        const auto *kh_240 = buffer.data(kh + 240 * ncomps + c);
        const auto *kh_241 = buffer.data(kh + 241 * ncomps + c);
        const auto *kh_242 = buffer.data(kh + 242 * ncomps + c);
        const auto *kh_243 = buffer.data(kh + 243 * ncomps + c);
        const auto *kh_244 = buffer.data(kh + 244 * ncomps + c);
        const auto *kh_245 = buffer.data(kh + 245 * ncomps + c);
        const auto *kh_246 = buffer.data(kh + 246 * ncomps + c);
        const auto *kh_247 = buffer.data(kh + 247 * ncomps + c);
        const auto *kh_248 = buffer.data(kh + 248 * ncomps + c);
        const auto *kh_249 = buffer.data(kh + 249 * ncomps + c);
        const auto *kh_250 = buffer.data(kh + 250 * ncomps + c);
        const auto *kh_251 = buffer.data(kh + 251 * ncomps + c);
        const auto *kh_252 = buffer.data(kh + 252 * ncomps + c);
        const auto *kh_253 = buffer.data(kh + 253 * ncomps + c);
        const auto *kh_254 = buffer.data(kh + 254 * ncomps + c);
        const auto *kh_255 = buffer.data(kh + 255 * ncomps + c);
        const auto *kh_256 = buffer.data(kh + 256 * ncomps + c);
        const auto *kh_257 = buffer.data(kh + 257 * ncomps + c);
        const auto *kh_258 = buffer.data(kh + 258 * ncomps + c);
        const auto *kh_259 = buffer.data(kh + 259 * ncomps + c);
        const auto *kh_260 = buffer.data(kh + 260 * ncomps + c);
        const auto *kh_261 = buffer.data(kh + 261 * ncomps + c);
        const auto *kh_262 = buffer.data(kh + 262 * ncomps + c);
        const auto *kh_263 = buffer.data(kh + 263 * ncomps + c);
        const auto *kh_264 = buffer.data(kh + 264 * ncomps + c);
        const auto *kh_265 = buffer.data(kh + 265 * ncomps + c);
        const auto *kh_266 = buffer.data(kh + 266 * ncomps + c);
        const auto *kh_267 = buffer.data(kh + 267 * ncomps + c);
        const auto *kh_268 = buffer.data(kh + 268 * ncomps + c);
        const auto *kh_269 = buffer.data(kh + 269 * ncomps + c);
        const auto *kh_270 = buffer.data(kh + 270 * ncomps + c);
        const auto *kh_271 = buffer.data(kh + 271 * ncomps + c);
        const auto *kh_272 = buffer.data(kh + 272 * ncomps + c);
        const auto *kh_273 = buffer.data(kh + 273 * ncomps + c);
        const auto *kh_274 = buffer.data(kh + 274 * ncomps + c);
        const auto *kh_275 = buffer.data(kh + 275 * ncomps + c);
        const auto *kh_276 = buffer.data(kh + 276 * ncomps + c);
        const auto *kh_277 = buffer.data(kh + 277 * ncomps + c);
        const auto *kh_278 = buffer.data(kh + 278 * ncomps + c);
        const auto *kh_279 = buffer.data(kh + 279 * ncomps + c);
        const auto *kh_280 = buffer.data(kh + 280 * ncomps + c);
        const auto *kh_281 = buffer.data(kh + 281 * ncomps + c);
        const auto *kh_282 = buffer.data(kh + 282 * ncomps + c);
        const auto *kh_283 = buffer.data(kh + 283 * ncomps + c);
        const auto *kh_284 = buffer.data(kh + 284 * ncomps + c);
        const auto *kh_285 = buffer.data(kh + 285 * ncomps + c);
        const auto *kh_286 = buffer.data(kh + 286 * ncomps + c);
        const auto *kh_287 = buffer.data(kh + 287 * ncomps + c);
        const auto *kh_288 = buffer.data(kh + 288 * ncomps + c);
        const auto *kh_289 = buffer.data(kh + 289 * ncomps + c);
        const auto *kh_290 = buffer.data(kh + 290 * ncomps + c);
        const auto *kh_291 = buffer.data(kh + 291 * ncomps + c);
        const auto *kh_292 = buffer.data(kh + 292 * ncomps + c);
        const auto *kh_293 = buffer.data(kh + 293 * ncomps + c);
        const auto *kh_294 = buffer.data(kh + 294 * ncomps + c);
        const auto *kh_295 = buffer.data(kh + 295 * ncomps + c);
        const auto *kh_296 = buffer.data(kh + 296 * ncomps + c);
        const auto *kh_297 = buffer.data(kh + 297 * ncomps + c);
        const auto *kh_298 = buffer.data(kh + 298 * ncomps + c);
        const auto *kh_299 = buffer.data(kh + 299 * ncomps + c);
        const auto *kh_300 = buffer.data(kh + 300 * ncomps + c);
        const auto *kh_301 = buffer.data(kh + 301 * ncomps + c);
        const auto *kh_302 = buffer.data(kh + 302 * ncomps + c);
        const auto *kh_303 = buffer.data(kh + 303 * ncomps + c);
        const auto *kh_304 = buffer.data(kh + 304 * ncomps + c);
        const auto *kh_305 = buffer.data(kh + 305 * ncomps + c);
        const auto *kh_306 = buffer.data(kh + 306 * ncomps + c);
        const auto *kh_307 = buffer.data(kh + 307 * ncomps + c);
        const auto *kh_308 = buffer.data(kh + 308 * ncomps + c);
        const auto *kh_309 = buffer.data(kh + 309 * ncomps + c);
        const auto *kh_310 = buffer.data(kh + 310 * ncomps + c);
        const auto *kh_311 = buffer.data(kh + 311 * ncomps + c);
        const auto *kh_312 = buffer.data(kh + 312 * ncomps + c);
        const auto *kh_313 = buffer.data(kh + 313 * ncomps + c);
        const auto *kh_314 = buffer.data(kh + 314 * ncomps + c);
        const auto *kh_315 = buffer.data(kh + 315 * ncomps + c);
        const auto *kh_316 = buffer.data(kh + 316 * ncomps + c);
        const auto *kh_317 = buffer.data(kh + 317 * ncomps + c);
        const auto *kh_318 = buffer.data(kh + 318 * ncomps + c);
        const auto *kh_319 = buffer.data(kh + 319 * ncomps + c);
        const auto *kh_320 = buffer.data(kh + 320 * ncomps + c);
        const auto *kh_321 = buffer.data(kh + 321 * ncomps + c);
        const auto *kh_322 = buffer.data(kh + 322 * ncomps + c);
        const auto *kh_323 = buffer.data(kh + 323 * ncomps + c);
        const auto *kh_324 = buffer.data(kh + 324 * ncomps + c);
        const auto *kh_325 = buffer.data(kh + 325 * ncomps + c);
        const auto *kh_326 = buffer.data(kh + 326 * ncomps + c);
        const auto *kh_327 = buffer.data(kh + 327 * ncomps + c);
        const auto *kh_328 = buffer.data(kh + 328 * ncomps + c);
        const auto *kh_329 = buffer.data(kh + 329 * ncomps + c);
        const auto *kh_330 = buffer.data(kh + 330 * ncomps + c);
        const auto *kh_331 = buffer.data(kh + 331 * ncomps + c);
        const auto *kh_332 = buffer.data(kh + 332 * ncomps + c);
        const auto *kh_333 = buffer.data(kh + 333 * ncomps + c);
        const auto *kh_334 = buffer.data(kh + 334 * ncomps + c);
        const auto *kh_335 = buffer.data(kh + 335 * ncomps + c);
        const auto *kh_351 = buffer.data(kh + 351 * ncomps + c);
        const auto *kh_352 = buffer.data(kh + 352 * ncomps + c);
        const auto *kh_353 = buffer.data(kh + 353 * ncomps + c);
        const auto *kh_354 = buffer.data(kh + 354 * ncomps + c);
        const auto *kh_355 = buffer.data(kh + 355 * ncomps + c);
        const auto *kh_356 = buffer.data(kh + 356 * ncomps + c);
        const auto *kh_372 = buffer.data(kh + 372 * ncomps + c);
        const auto *kh_373 = buffer.data(kh + 373 * ncomps + c);
        const auto *kh_374 = buffer.data(kh + 374 * ncomps + c);
        const auto *kh_375 = buffer.data(kh + 375 * ncomps + c);
        const auto *kh_376 = buffer.data(kh + 376 * ncomps + c);
        const auto *kh_377 = buffer.data(kh + 377 * ncomps + c);
        const auto *kh_393 = buffer.data(kh + 393 * ncomps + c);
        const auto *kh_394 = buffer.data(kh + 394 * ncomps + c);
        const auto *kh_395 = buffer.data(kh + 395 * ncomps + c);
        const auto *kh_396 = buffer.data(kh + 396 * ncomps + c);
        const auto *kh_397 = buffer.data(kh + 397 * ncomps + c);
        const auto *kh_398 = buffer.data(kh + 398 * ncomps + c);
        const auto *kh_414 = buffer.data(kh + 414 * ncomps + c);
        const auto *kh_415 = buffer.data(kh + 415 * ncomps + c);
        const auto *kh_416 = buffer.data(kh + 416 * ncomps + c);
        const auto *kh_417 = buffer.data(kh + 417 * ncomps + c);
        const auto *kh_418 = buffer.data(kh + 418 * ncomps + c);
        const auto *kh_419 = buffer.data(kh + 419 * ncomps + c);
        const auto *kh_440 = buffer.data(kh + 440 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, ih_220, ih_221, ih_222, \
                         ih_223, ih_224, kh_220, kh_221, kh_222, kh_223, \
                         kh_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * ih_220[k]
                       + kh_220[k];

            t_291[k] = ab_x[k] * ih_221[k]
                       + kh_221[k];

            t_292[k] = ab_x[k] * ih_222[k]
                       + kh_222[k];

            t_293[k] = ab_x[k] * ih_223[k]
                       + kh_223[k];

            t_294[k] = ab_x[k] * ih_224[k]
                       + kh_224[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, ih_225, ih_226, ih_227, \
                         ih_228, ih_229, kh_225, kh_226, kh_227, kh_228, \
                         kh_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_x[k] * ih_225[k]
                       + kh_225[k];

            t_296[k] = ab_x[k] * ih_226[k]
                       + kh_226[k];

            t_297[k] = ab_x[k] * ih_227[k]
                       + kh_227[k];

            t_298[k] = ab_x[k] * ih_228[k]
                       + kh_228[k];

            t_299[k] = ab_x[k] * ih_229[k]
                       + kh_229[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, ab_x, ab_y, ih_225, ih_226, ih_227, \
                         ih_230, kh_230, kh_330, kh_331, kh_332 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * ih_230[k]
                       + kh_230[k];

            t_301[k] = ab_y[k] * ih_225[k]
                       + kh_330[k];

            t_302[k] = ab_y[k] * ih_226[k]
                       + kh_331[k];

            t_303[k] = ab_y[k] * ih_227[k]
                       + kh_332[k];
        }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, ab_y, ab_z, ih_228, ih_229, ih_230, \
                         kh_333, kh_334, kh_335, kh_356 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_304[k] = ab_y[k] * ih_228[k]
                       + kh_333[k];

            t_305[k] = ab_y[k] * ih_229[k]
                       + kh_334[k];

            t_306[k] = ab_y[k] * ih_230[k]
                       + kh_335[k];

            t_307[k] = ab_z[k] * ih_230[k]
                       + kh_356[k];
        }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, ab_x, ih_231, ih_232, ih_233, \
                         ih_234, ih_235, kh_231, kh_232, kh_233, kh_234, \
                         kh_235 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_308[k] = ab_x[k] * ih_231[k]
                       + kh_231[k];

            t_309[k] = ab_x[k] * ih_232[k]
                       + kh_232[k];

            t_310[k] = ab_x[k] * ih_233[k]
                       + kh_233[k];

            t_311[k] = ab_x[k] * ih_234[k]
                       + kh_234[k];

            t_312[k] = ab_x[k] * ih_235[k]
                       + kh_235[k];
        }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, ab_x, ih_236, ih_237, ih_238, \
                         ih_239, ih_240, kh_236, kh_237, kh_238, kh_239, \
                         kh_240 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_313[k] = ab_x[k] * ih_236[k]
                       + kh_236[k];

            t_314[k] = ab_x[k] * ih_237[k]
                       + kh_237[k];

            t_315[k] = ab_x[k] * ih_238[k]
                       + kh_238[k];

            t_316[k] = ab_x[k] * ih_239[k]
                       + kh_239[k];

            t_317[k] = ab_x[k] * ih_240[k]
                       + kh_240[k];
        }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, ab_x, ih_241, ih_242, ih_243, \
                         ih_244, ih_245, kh_241, kh_242, kh_243, kh_244, \
                         kh_245 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_318[k] = ab_x[k] * ih_241[k]
                       + kh_241[k];

            t_319[k] = ab_x[k] * ih_242[k]
                       + kh_242[k];

            t_320[k] = ab_x[k] * ih_243[k]
                       + kh_243[k];

            t_321[k] = ab_x[k] * ih_244[k]
                       + kh_244[k];

            t_322[k] = ab_x[k] * ih_245[k]
                       + kh_245[k];
        }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, ab_x, ih_246, ih_247, ih_248, \
                         ih_249, ih_250, kh_246, kh_247, kh_248, kh_249, \
                         kh_250 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_323[k] = ab_x[k] * ih_246[k]
                       + kh_246[k];

            t_324[k] = ab_x[k] * ih_247[k]
                       + kh_247[k];

            t_325[k] = ab_x[k] * ih_248[k]
                       + kh_248[k];

            t_326[k] = ab_x[k] * ih_249[k]
                       + kh_249[k];

            t_327[k] = ab_x[k] * ih_250[k]
                       + kh_250[k];
        }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, ab_x, ab_y, ih_246, ih_247, ih_248, \
                         ih_251, kh_251, kh_351, kh_352, kh_353 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_328[k] = ab_x[k] * ih_251[k]
                       + kh_251[k];

            t_329[k] = ab_y[k] * ih_246[k]
                       + kh_351[k];

            t_330[k] = ab_y[k] * ih_247[k]
                       + kh_352[k];

            t_331[k] = ab_y[k] * ih_248[k]
                       + kh_353[k];
        }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, ab_y, ab_z, ih_249, ih_250, ih_251, \
                         kh_354, kh_355, kh_356, kh_377 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_332[k] = ab_y[k] * ih_249[k]
                       + kh_354[k];

            t_333[k] = ab_y[k] * ih_250[k]
                       + kh_355[k];

            t_334[k] = ab_y[k] * ih_251[k]
                       + kh_356[k];

            t_335[k] = ab_z[k] * ih_251[k]
                       + kh_377[k];
        }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, ab_x, ih_252, ih_253, ih_254, \
                         ih_255, ih_256, kh_252, kh_253, kh_254, kh_255, \
                         kh_256 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_336[k] = ab_x[k] * ih_252[k]
                       + kh_252[k];

            t_337[k] = ab_x[k] * ih_253[k]
                       + kh_253[k];

            t_338[k] = ab_x[k] * ih_254[k]
                       + kh_254[k];

            t_339[k] = ab_x[k] * ih_255[k]
                       + kh_255[k];

            t_340[k] = ab_x[k] * ih_256[k]
                       + kh_256[k];
        }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, ab_x, ih_257, ih_258, ih_259, \
                         ih_260, ih_261, kh_257, kh_258, kh_259, kh_260, \
                         kh_261 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_341[k] = ab_x[k] * ih_257[k]
                       + kh_257[k];

            t_342[k] = ab_x[k] * ih_258[k]
                       + kh_258[k];

            t_343[k] = ab_x[k] * ih_259[k]
                       + kh_259[k];

            t_344[k] = ab_x[k] * ih_260[k]
                       + kh_260[k];

            t_345[k] = ab_x[k] * ih_261[k]
                       + kh_261[k];
        }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, ab_x, ih_262, ih_263, ih_264, \
                         ih_265, ih_266, kh_262, kh_263, kh_264, kh_265, \
                         kh_266 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_346[k] = ab_x[k] * ih_262[k]
                       + kh_262[k];

            t_347[k] = ab_x[k] * ih_263[k]
                       + kh_263[k];

            t_348[k] = ab_x[k] * ih_264[k]
                       + kh_264[k];

            t_349[k] = ab_x[k] * ih_265[k]
                       + kh_265[k];

            t_350[k] = ab_x[k] * ih_266[k]
                       + kh_266[k];
        }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, ab_x, ih_267, ih_268, ih_269, \
                         ih_270, ih_271, kh_267, kh_268, kh_269, kh_270, \
                         kh_271 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_351[k] = ab_x[k] * ih_267[k]
                       + kh_267[k];

            t_352[k] = ab_x[k] * ih_268[k]
                       + kh_268[k];

            t_353[k] = ab_x[k] * ih_269[k]
                       + kh_269[k];

            t_354[k] = ab_x[k] * ih_270[k]
                       + kh_270[k];

            t_355[k] = ab_x[k] * ih_271[k]
                       + kh_271[k];
        }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, ab_x, ab_y, ih_267, ih_268, ih_269, \
                         ih_272, kh_272, kh_372, kh_373, kh_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_356[k] = ab_x[k] * ih_272[k]
                       + kh_272[k];

            t_357[k] = ab_y[k] * ih_267[k]
                       + kh_372[k];

            t_358[k] = ab_y[k] * ih_268[k]
                       + kh_373[k];

            t_359[k] = ab_y[k] * ih_269[k]
                       + kh_374[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, ab_y, ab_z, ih_270, ih_271, ih_272, \
                         kh_375, kh_376, kh_377, kh_398 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_y[k] * ih_270[k]
                       + kh_375[k];

            t_361[k] = ab_y[k] * ih_271[k]
                       + kh_376[k];

            t_362[k] = ab_y[k] * ih_272[k]
                       + kh_377[k];

            t_363[k] = ab_z[k] * ih_272[k]
                       + kh_398[k];
        }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, ab_x, ih_273, ih_274, ih_275, \
                         ih_276, ih_277, kh_273, kh_274, kh_275, kh_276, \
                         kh_277 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_364[k] = ab_x[k] * ih_273[k]
                       + kh_273[k];

            t_365[k] = ab_x[k] * ih_274[k]
                       + kh_274[k];

            t_366[k] = ab_x[k] * ih_275[k]
                       + kh_275[k];

            t_367[k] = ab_x[k] * ih_276[k]
                       + kh_276[k];

            t_368[k] = ab_x[k] * ih_277[k]
                       + kh_277[k];
        }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, t_373, ab_x, ih_278, ih_279, ih_280, \
                         ih_281, ih_282, kh_278, kh_279, kh_280, kh_281, \
                         kh_282 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_369[k] = ab_x[k] * ih_278[k]
                       + kh_278[k];

            t_370[k] = ab_x[k] * ih_279[k]
                       + kh_279[k];

            t_371[k] = ab_x[k] * ih_280[k]
                       + kh_280[k];

            t_372[k] = ab_x[k] * ih_281[k]
                       + kh_281[k];

            t_373[k] = ab_x[k] * ih_282[k]
                       + kh_282[k];
        }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, ab_x, ih_283, ih_284, ih_285, \
                         ih_286, ih_287, kh_283, kh_284, kh_285, kh_286, \
                         kh_287 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_374[k] = ab_x[k] * ih_283[k]
                       + kh_283[k];

            t_375[k] = ab_x[k] * ih_284[k]
                       + kh_284[k];

            t_376[k] = ab_x[k] * ih_285[k]
                       + kh_285[k];

            t_377[k] = ab_x[k] * ih_286[k]
                       + kh_286[k];

            t_378[k] = ab_x[k] * ih_287[k]
                       + kh_287[k];
        }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, ab_x, ih_288, ih_289, ih_290, \
                         ih_291, ih_292, kh_288, kh_289, kh_290, kh_291, \
                         kh_292 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_379[k] = ab_x[k] * ih_288[k]
                       + kh_288[k];

            t_380[k] = ab_x[k] * ih_289[k]
                       + kh_289[k];

            t_381[k] = ab_x[k] * ih_290[k]
                       + kh_290[k];

            t_382[k] = ab_x[k] * ih_291[k]
                       + kh_291[k];

            t_383[k] = ab_x[k] * ih_292[k]
                       + kh_292[k];
        }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, ab_x, ab_y, ih_288, ih_289, ih_290, \
                         ih_293, kh_293, kh_393, kh_394, kh_395 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_384[k] = ab_x[k] * ih_293[k]
                       + kh_293[k];

            t_385[k] = ab_y[k] * ih_288[k]
                       + kh_393[k];

            t_386[k] = ab_y[k] * ih_289[k]
                       + kh_394[k];

            t_387[k] = ab_y[k] * ih_290[k]
                       + kh_395[k];
        }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, ab_y, ab_z, ih_291, ih_292, ih_293, \
                         kh_396, kh_397, kh_398, kh_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_388[k] = ab_y[k] * ih_291[k]
                       + kh_396[k];

            t_389[k] = ab_y[k] * ih_292[k]
                       + kh_397[k];

            t_390[k] = ab_y[k] * ih_293[k]
                       + kh_398[k];

            t_391[k] = ab_z[k] * ih_293[k]
                       + kh_419[k];
        }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, ab_x, ih_294, ih_295, ih_296, \
                         ih_297, ih_298, kh_294, kh_295, kh_296, kh_297, \
                         kh_298 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_392[k] = ab_x[k] * ih_294[k]
                       + kh_294[k];

            t_393[k] = ab_x[k] * ih_295[k]
                       + kh_295[k];

            t_394[k] = ab_x[k] * ih_296[k]
                       + kh_296[k];

            t_395[k] = ab_x[k] * ih_297[k]
                       + kh_297[k];

            t_396[k] = ab_x[k] * ih_298[k]
                       + kh_298[k];
        }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, ab_x, ih_299, ih_300, ih_301, \
                         ih_302, ih_303, kh_299, kh_300, kh_301, kh_302, \
                         kh_303 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_397[k] = ab_x[k] * ih_299[k]
                       + kh_299[k];

            t_398[k] = ab_x[k] * ih_300[k]
                       + kh_300[k];

            t_399[k] = ab_x[k] * ih_301[k]
                       + kh_301[k];

            t_400[k] = ab_x[k] * ih_302[k]
                       + kh_302[k];

            t_401[k] = ab_x[k] * ih_303[k]
                       + kh_303[k];
        }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ab_x, ih_304, ih_305, ih_306, \
                         ih_307, ih_308, kh_304, kh_305, kh_306, kh_307, \
                         kh_308 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_402[k] = ab_x[k] * ih_304[k]
                       + kh_304[k];

            t_403[k] = ab_x[k] * ih_305[k]
                       + kh_305[k];

            t_404[k] = ab_x[k] * ih_306[k]
                       + kh_306[k];

            t_405[k] = ab_x[k] * ih_307[k]
                       + kh_307[k];

            t_406[k] = ab_x[k] * ih_308[k]
                       + kh_308[k];
        }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ab_x, ih_309, ih_310, ih_311, \
                         ih_312, ih_313, kh_309, kh_310, kh_311, kh_312, \
                         kh_313 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_407[k] = ab_x[k] * ih_309[k]
                       + kh_309[k];

            t_408[k] = ab_x[k] * ih_310[k]
                       + kh_310[k];

            t_409[k] = ab_x[k] * ih_311[k]
                       + kh_311[k];

            t_410[k] = ab_x[k] * ih_312[k]
                       + kh_312[k];

            t_411[k] = ab_x[k] * ih_313[k]
                       + kh_313[k];
        }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, ab_x, ab_y, ih_309, ih_310, ih_311, \
                         ih_314, kh_314, kh_414, kh_415, kh_416 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_412[k] = ab_x[k] * ih_314[k]
                       + kh_314[k];

            t_413[k] = ab_y[k] * ih_309[k]
                       + kh_414[k];

            t_414[k] = ab_y[k] * ih_310[k]
                       + kh_415[k];

            t_415[k] = ab_y[k] * ih_311[k]
                       + kh_416[k];
        }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, ab_y, ab_z, ih_312, ih_313, ih_314, \
                         kh_417, kh_418, kh_419, kh_440 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_416[k] = ab_y[k] * ih_312[k]
                       + kh_417[k];

            t_417[k] = ab_y[k] * ih_313[k]
                       + kh_418[k];

            t_418[k] = ab_y[k] * ih_314[k]
                       + kh_419[k];

            t_419[k] = ab_z[k] * ih_314[k]
                       + kh_440[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, ih_315, ih_316, ih_317, \
                         ih_318, ih_319, kh_315, kh_316, kh_317, kh_318, \
                         kh_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * ih_315[k]
                       + kh_315[k];

            t_421[k] = ab_x[k] * ih_316[k]
                       + kh_316[k];

            t_422[k] = ab_x[k] * ih_317[k]
                       + kh_317[k];

            t_423[k] = ab_x[k] * ih_318[k]
                       + kh_318[k];

            t_424[k] = ab_x[k] * ih_319[k]
                       + kh_319[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, ih_320, ih_321, ih_322, \
                         ih_323, ih_324, kh_320, kh_321, kh_322, kh_323, \
                         kh_324 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * ih_320[k]
                       + kh_320[k];

            t_426[k] = ab_x[k] * ih_321[k]
                       + kh_321[k];

            t_427[k] = ab_x[k] * ih_322[k]
                       + kh_322[k];

            t_428[k] = ab_x[k] * ih_323[k]
                       + kh_323[k];

            t_429[k] = ab_x[k] * ih_324[k]
                       + kh_324[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, ih_325, ih_326, ih_327, \
                         ih_328, ih_329, kh_325, kh_326, kh_327, kh_328, \
                         kh_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * ih_325[k]
                       + kh_325[k];

            t_431[k] = ab_x[k] * ih_326[k]
                       + kh_326[k];

            t_432[k] = ab_x[k] * ih_327[k]
                       + kh_327[k];

            t_433[k] = ab_x[k] * ih_328[k]
                       + kh_328[k];

            t_434[k] = ab_x[k] * ih_329[k]
                       + kh_329[k];
        }
    }
}

static auto
compute_hrr_ii_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ih, const size_t kh,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_435 = buffer.data(target + 435 * ncomps + c);
        auto *t_436 = buffer.data(target + 436 * ncomps + c);
        auto *t_437 = buffer.data(target + 437 * ncomps + c);
        auto *t_438 = buffer.data(target + 438 * ncomps + c);
        auto *t_439 = buffer.data(target + 439 * ncomps + c);
        auto *t_440 = buffer.data(target + 440 * ncomps + c);
        auto *t_441 = buffer.data(target + 441 * ncomps + c);
        auto *t_442 = buffer.data(target + 442 * ncomps + c);
        auto *t_443 = buffer.data(target + 443 * ncomps + c);
        auto *t_444 = buffer.data(target + 444 * ncomps + c);
        auto *t_445 = buffer.data(target + 445 * ncomps + c);
        auto *t_446 = buffer.data(target + 446 * ncomps + c);
        auto *t_447 = buffer.data(target + 447 * ncomps + c);
        auto *t_448 = buffer.data(target + 448 * ncomps + c);
        auto *t_449 = buffer.data(target + 449 * ncomps + c);
        auto *t_450 = buffer.data(target + 450 * ncomps + c);
        auto *t_451 = buffer.data(target + 451 * ncomps + c);
        auto *t_452 = buffer.data(target + 452 * ncomps + c);
        auto *t_453 = buffer.data(target + 453 * ncomps + c);
        auto *t_454 = buffer.data(target + 454 * ncomps + c);
        auto *t_455 = buffer.data(target + 455 * ncomps + c);
        auto *t_456 = buffer.data(target + 456 * ncomps + c);
        auto *t_457 = buffer.data(target + 457 * ncomps + c);
        auto *t_458 = buffer.data(target + 458 * ncomps + c);
        auto *t_459 = buffer.data(target + 459 * ncomps + c);
        auto *t_460 = buffer.data(target + 460 * ncomps + c);
        auto *t_461 = buffer.data(target + 461 * ncomps + c);
        auto *t_462 = buffer.data(target + 462 * ncomps + c);
        auto *t_463 = buffer.data(target + 463 * ncomps + c);
        auto *t_464 = buffer.data(target + 464 * ncomps + c);
        auto *t_465 = buffer.data(target + 465 * ncomps + c);
        auto *t_466 = buffer.data(target + 466 * ncomps + c);
        auto *t_467 = buffer.data(target + 467 * ncomps + c);
        auto *t_468 = buffer.data(target + 468 * ncomps + c);
        auto *t_469 = buffer.data(target + 469 * ncomps + c);
        auto *t_470 = buffer.data(target + 470 * ncomps + c);
        auto *t_471 = buffer.data(target + 471 * ncomps + c);
        auto *t_472 = buffer.data(target + 472 * ncomps + c);
        auto *t_473 = buffer.data(target + 473 * ncomps + c);
        auto *t_474 = buffer.data(target + 474 * ncomps + c);
        auto *t_475 = buffer.data(target + 475 * ncomps + c);
        auto *t_476 = buffer.data(target + 476 * ncomps + c);
        auto *t_477 = buffer.data(target + 477 * ncomps + c);
        auto *t_478 = buffer.data(target + 478 * ncomps + c);
        auto *t_479 = buffer.data(target + 479 * ncomps + c);
        auto *t_480 = buffer.data(target + 480 * ncomps + c);
        auto *t_481 = buffer.data(target + 481 * ncomps + c);
        auto *t_482 = buffer.data(target + 482 * ncomps + c);
        auto *t_483 = buffer.data(target + 483 * ncomps + c);
        auto *t_484 = buffer.data(target + 484 * ncomps + c);
        auto *t_485 = buffer.data(target + 485 * ncomps + c);
        auto *t_486 = buffer.data(target + 486 * ncomps + c);
        auto *t_487 = buffer.data(target + 487 * ncomps + c);
        auto *t_488 = buffer.data(target + 488 * ncomps + c);
        auto *t_489 = buffer.data(target + 489 * ncomps + c);
        auto *t_490 = buffer.data(target + 490 * ncomps + c);
        auto *t_491 = buffer.data(target + 491 * ncomps + c);
        auto *t_492 = buffer.data(target + 492 * ncomps + c);
        auto *t_493 = buffer.data(target + 493 * ncomps + c);
        auto *t_494 = buffer.data(target + 494 * ncomps + c);
        auto *t_495 = buffer.data(target + 495 * ncomps + c);
        auto *t_496 = buffer.data(target + 496 * ncomps + c);
        auto *t_497 = buffer.data(target + 497 * ncomps + c);
        auto *t_498 = buffer.data(target + 498 * ncomps + c);
        auto *t_499 = buffer.data(target + 499 * ncomps + c);
        auto *t_500 = buffer.data(target + 500 * ncomps + c);
        auto *t_501 = buffer.data(target + 501 * ncomps + c);
        auto *t_502 = buffer.data(target + 502 * ncomps + c);
        auto *t_503 = buffer.data(target + 503 * ncomps + c);
        auto *t_504 = buffer.data(target + 504 * ncomps + c);
        auto *t_505 = buffer.data(target + 505 * ncomps + c);
        auto *t_506 = buffer.data(target + 506 * ncomps + c);
        auto *t_507 = buffer.data(target + 507 * ncomps + c);
        auto *t_508 = buffer.data(target + 508 * ncomps + c);
        auto *t_509 = buffer.data(target + 509 * ncomps + c);
        auto *t_510 = buffer.data(target + 510 * ncomps + c);
        auto *t_511 = buffer.data(target + 511 * ncomps + c);
        auto *t_512 = buffer.data(target + 512 * ncomps + c);
        auto *t_513 = buffer.data(target + 513 * ncomps + c);
        auto *t_514 = buffer.data(target + 514 * ncomps + c);
        auto *t_515 = buffer.data(target + 515 * ncomps + c);
        auto *t_516 = buffer.data(target + 516 * ncomps + c);
        auto *t_517 = buffer.data(target + 517 * ncomps + c);
        auto *t_518 = buffer.data(target + 518 * ncomps + c);
        auto *t_519 = buffer.data(target + 519 * ncomps + c);
        auto *t_520 = buffer.data(target + 520 * ncomps + c);
        auto *t_521 = buffer.data(target + 521 * ncomps + c);
        auto *t_522 = buffer.data(target + 522 * ncomps + c);
        auto *t_523 = buffer.data(target + 523 * ncomps + c);
        auto *t_524 = buffer.data(target + 524 * ncomps + c);
        auto *t_525 = buffer.data(target + 525 * ncomps + c);
        auto *t_526 = buffer.data(target + 526 * ncomps + c);
        auto *t_527 = buffer.data(target + 527 * ncomps + c);
        auto *t_528 = buffer.data(target + 528 * ncomps + c);
        auto *t_529 = buffer.data(target + 529 * ncomps + c);
        auto *t_530 = buffer.data(target + 530 * ncomps + c);
        auto *t_531 = buffer.data(target + 531 * ncomps + c);
        auto *t_532 = buffer.data(target + 532 * ncomps + c);
        auto *t_533 = buffer.data(target + 533 * ncomps + c);
        auto *t_534 = buffer.data(target + 534 * ncomps + c);
        auto *t_535 = buffer.data(target + 535 * ncomps + c);
        auto *t_536 = buffer.data(target + 536 * ncomps + c);
        auto *t_537 = buffer.data(target + 537 * ncomps + c);
        auto *t_538 = buffer.data(target + 538 * ncomps + c);
        auto *t_539 = buffer.data(target + 539 * ncomps + c);
        auto *t_540 = buffer.data(target + 540 * ncomps + c);
        auto *t_541 = buffer.data(target + 541 * ncomps + c);
        auto *t_542 = buffer.data(target + 542 * ncomps + c);
        auto *t_543 = buffer.data(target + 543 * ncomps + c);
        auto *t_544 = buffer.data(target + 544 * ncomps + c);
        auto *t_545 = buffer.data(target + 545 * ncomps + c);
        auto *t_546 = buffer.data(target + 546 * ncomps + c);
        auto *t_547 = buffer.data(target + 547 * ncomps + c);
        auto *t_548 = buffer.data(target + 548 * ncomps + c);
        auto *t_549 = buffer.data(target + 549 * ncomps + c);
        auto *t_550 = buffer.data(target + 550 * ncomps + c);
        auto *t_551 = buffer.data(target + 551 * ncomps + c);
        auto *t_552 = buffer.data(target + 552 * ncomps + c);
        auto *t_553 = buffer.data(target + 553 * ncomps + c);
        auto *t_554 = buffer.data(target + 554 * ncomps + c);
        auto *t_555 = buffer.data(target + 555 * ncomps + c);
        auto *t_556 = buffer.data(target + 556 * ncomps + c);
        auto *t_557 = buffer.data(target + 557 * ncomps + c);
        auto *t_558 = buffer.data(target + 558 * ncomps + c);
        auto *t_559 = buffer.data(target + 559 * ncomps + c);
        auto *t_560 = buffer.data(target + 560 * ncomps + c);
        auto *t_561 = buffer.data(target + 561 * ncomps + c);
        auto *t_562 = buffer.data(target + 562 * ncomps + c);
        auto *t_563 = buffer.data(target + 563 * ncomps + c);
        auto *t_564 = buffer.data(target + 564 * ncomps + c);
        auto *t_565 = buffer.data(target + 565 * ncomps + c);
        auto *t_566 = buffer.data(target + 566 * ncomps + c);
        auto *t_567 = buffer.data(target + 567 * ncomps + c);
        auto *t_568 = buffer.data(target + 568 * ncomps + c);
        auto *t_569 = buffer.data(target + 569 * ncomps + c);
        auto *t_570 = buffer.data(target + 570 * ncomps + c);
        auto *t_571 = buffer.data(target + 571 * ncomps + c);
        auto *t_572 = buffer.data(target + 572 * ncomps + c);
        auto *t_573 = buffer.data(target + 573 * ncomps + c);
        auto *t_574 = buffer.data(target + 574 * ncomps + c);
        auto *t_575 = buffer.data(target + 575 * ncomps + c);
        auto *t_576 = buffer.data(target + 576 * ncomps + c);
        auto *t_577 = buffer.data(target + 577 * ncomps + c);
        auto *t_578 = buffer.data(target + 578 * ncomps + c);
        auto *t_579 = buffer.data(target + 579 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ih_330 = buffer.data(ih + 330 * ncomps + c);
        const auto *ih_331 = buffer.data(ih + 331 * ncomps + c);
        const auto *ih_332 = buffer.data(ih + 332 * ncomps + c);
        const auto *ih_333 = buffer.data(ih + 333 * ncomps + c);
        const auto *ih_334 = buffer.data(ih + 334 * ncomps + c);
        const auto *ih_335 = buffer.data(ih + 335 * ncomps + c);
        const auto *ih_336 = buffer.data(ih + 336 * ncomps + c);
        const auto *ih_337 = buffer.data(ih + 337 * ncomps + c);
        const auto *ih_338 = buffer.data(ih + 338 * ncomps + c);
        const auto *ih_339 = buffer.data(ih + 339 * ncomps + c);
        const auto *ih_340 = buffer.data(ih + 340 * ncomps + c);
        const auto *ih_341 = buffer.data(ih + 341 * ncomps + c);
        const auto *ih_342 = buffer.data(ih + 342 * ncomps + c);
        const auto *ih_343 = buffer.data(ih + 343 * ncomps + c);
        const auto *ih_344 = buffer.data(ih + 344 * ncomps + c);
        const auto *ih_345 = buffer.data(ih + 345 * ncomps + c);
        const auto *ih_346 = buffer.data(ih + 346 * ncomps + c);
        const auto *ih_347 = buffer.data(ih + 347 * ncomps + c);
        const auto *ih_348 = buffer.data(ih + 348 * ncomps + c);
        const auto *ih_349 = buffer.data(ih + 349 * ncomps + c);
        const auto *ih_350 = buffer.data(ih + 350 * ncomps + c);
        const auto *ih_351 = buffer.data(ih + 351 * ncomps + c);
        const auto *ih_352 = buffer.data(ih + 352 * ncomps + c);
        const auto *ih_353 = buffer.data(ih + 353 * ncomps + c);
        const auto *ih_354 = buffer.data(ih + 354 * ncomps + c);
        const auto *ih_355 = buffer.data(ih + 355 * ncomps + c);
        const auto *ih_356 = buffer.data(ih + 356 * ncomps + c);
        const auto *ih_357 = buffer.data(ih + 357 * ncomps + c);
        const auto *ih_358 = buffer.data(ih + 358 * ncomps + c);
        const auto *ih_359 = buffer.data(ih + 359 * ncomps + c);
        const auto *ih_360 = buffer.data(ih + 360 * ncomps + c);
        const auto *ih_361 = buffer.data(ih + 361 * ncomps + c);
        const auto *ih_362 = buffer.data(ih + 362 * ncomps + c);
        const auto *ih_363 = buffer.data(ih + 363 * ncomps + c);
        const auto *ih_364 = buffer.data(ih + 364 * ncomps + c);
        const auto *ih_365 = buffer.data(ih + 365 * ncomps + c);
        const auto *ih_366 = buffer.data(ih + 366 * ncomps + c);
        const auto *ih_367 = buffer.data(ih + 367 * ncomps + c);
        const auto *ih_368 = buffer.data(ih + 368 * ncomps + c);
        const auto *ih_369 = buffer.data(ih + 369 * ncomps + c);
        const auto *ih_370 = buffer.data(ih + 370 * ncomps + c);
        const auto *ih_371 = buffer.data(ih + 371 * ncomps + c);
        const auto *ih_372 = buffer.data(ih + 372 * ncomps + c);
        const auto *ih_373 = buffer.data(ih + 373 * ncomps + c);
        const auto *ih_374 = buffer.data(ih + 374 * ncomps + c);
        const auto *ih_375 = buffer.data(ih + 375 * ncomps + c);
        const auto *ih_376 = buffer.data(ih + 376 * ncomps + c);
        const auto *ih_377 = buffer.data(ih + 377 * ncomps + c);
        const auto *ih_378 = buffer.data(ih + 378 * ncomps + c);
        const auto *ih_379 = buffer.data(ih + 379 * ncomps + c);
        const auto *ih_380 = buffer.data(ih + 380 * ncomps + c);
        const auto *ih_381 = buffer.data(ih + 381 * ncomps + c);
        const auto *ih_382 = buffer.data(ih + 382 * ncomps + c);
        const auto *ih_383 = buffer.data(ih + 383 * ncomps + c);
        const auto *ih_384 = buffer.data(ih + 384 * ncomps + c);
        const auto *ih_385 = buffer.data(ih + 385 * ncomps + c);
        const auto *ih_386 = buffer.data(ih + 386 * ncomps + c);
        const auto *ih_387 = buffer.data(ih + 387 * ncomps + c);
        const auto *ih_388 = buffer.data(ih + 388 * ncomps + c);
        const auto *ih_389 = buffer.data(ih + 389 * ncomps + c);
        const auto *ih_390 = buffer.data(ih + 390 * ncomps + c);
        const auto *ih_391 = buffer.data(ih + 391 * ncomps + c);
        const auto *ih_392 = buffer.data(ih + 392 * ncomps + c);
        const auto *ih_393 = buffer.data(ih + 393 * ncomps + c);
        const auto *ih_394 = buffer.data(ih + 394 * ncomps + c);
        const auto *ih_395 = buffer.data(ih + 395 * ncomps + c);
        const auto *ih_396 = buffer.data(ih + 396 * ncomps + c);
        const auto *ih_397 = buffer.data(ih + 397 * ncomps + c);
        const auto *ih_398 = buffer.data(ih + 398 * ncomps + c);
        const auto *ih_399 = buffer.data(ih + 399 * ncomps + c);
        const auto *ih_400 = buffer.data(ih + 400 * ncomps + c);
        const auto *ih_401 = buffer.data(ih + 401 * ncomps + c);
        const auto *ih_402 = buffer.data(ih + 402 * ncomps + c);
        const auto *ih_403 = buffer.data(ih + 403 * ncomps + c);
        const auto *ih_404 = buffer.data(ih + 404 * ncomps + c);
        const auto *ih_405 = buffer.data(ih + 405 * ncomps + c);
        const auto *ih_406 = buffer.data(ih + 406 * ncomps + c);
        const auto *ih_407 = buffer.data(ih + 407 * ncomps + c);
        const auto *ih_408 = buffer.data(ih + 408 * ncomps + c);
        const auto *ih_409 = buffer.data(ih + 409 * ncomps + c);
        const auto *ih_410 = buffer.data(ih + 410 * ncomps + c);
        const auto *ih_411 = buffer.data(ih + 411 * ncomps + c);
        const auto *ih_412 = buffer.data(ih + 412 * ncomps + c);
        const auto *ih_413 = buffer.data(ih + 413 * ncomps + c);
        const auto *ih_414 = buffer.data(ih + 414 * ncomps + c);
        const auto *ih_415 = buffer.data(ih + 415 * ncomps + c);
        const auto *ih_416 = buffer.data(ih + 416 * ncomps + c);
        const auto *ih_417 = buffer.data(ih + 417 * ncomps + c);
        const auto *ih_418 = buffer.data(ih + 418 * ncomps + c);
        const auto *ih_419 = buffer.data(ih + 419 * ncomps + c);
        const auto *ih_420 = buffer.data(ih + 420 * ncomps + c);
        const auto *ih_421 = buffer.data(ih + 421 * ncomps + c);
        const auto *ih_422 = buffer.data(ih + 422 * ncomps + c);
        const auto *ih_423 = buffer.data(ih + 423 * ncomps + c);
        const auto *ih_424 = buffer.data(ih + 424 * ncomps + c);
        const auto *ih_425 = buffer.data(ih + 425 * ncomps + c);
        const auto *ih_426 = buffer.data(ih + 426 * ncomps + c);
        const auto *ih_427 = buffer.data(ih + 427 * ncomps + c);
        const auto *ih_428 = buffer.data(ih + 428 * ncomps + c);
        const auto *ih_429 = buffer.data(ih + 429 * ncomps + c);
        const auto *ih_430 = buffer.data(ih + 430 * ncomps + c);
        const auto *ih_431 = buffer.data(ih + 431 * ncomps + c);
        const auto *ih_432 = buffer.data(ih + 432 * ncomps + c);
        const auto *ih_433 = buffer.data(ih + 433 * ncomps + c);
        const auto *ih_434 = buffer.data(ih + 434 * ncomps + c);
        const auto *ih_435 = buffer.data(ih + 435 * ncomps + c);
        const auto *ih_436 = buffer.data(ih + 436 * ncomps + c);
        const auto *ih_437 = buffer.data(ih + 437 * ncomps + c);
        const auto *ih_438 = buffer.data(ih + 438 * ncomps + c);
        const auto *ih_439 = buffer.data(ih + 439 * ncomps + c);

        const auto *kh_330 = buffer.data(kh + 330 * ncomps + c);
        const auto *kh_331 = buffer.data(kh + 331 * ncomps + c);
        const auto *kh_332 = buffer.data(kh + 332 * ncomps + c);
        const auto *kh_333 = buffer.data(kh + 333 * ncomps + c);
        const auto *kh_334 = buffer.data(kh + 334 * ncomps + c);
        const auto *kh_335 = buffer.data(kh + 335 * ncomps + c);
        const auto *kh_336 = buffer.data(kh + 336 * ncomps + c);
        const auto *kh_337 = buffer.data(kh + 337 * ncomps + c);
        const auto *kh_338 = buffer.data(kh + 338 * ncomps + c);
        const auto *kh_339 = buffer.data(kh + 339 * ncomps + c);
        const auto *kh_340 = buffer.data(kh + 340 * ncomps + c);
        const auto *kh_341 = buffer.data(kh + 341 * ncomps + c);
        const auto *kh_342 = buffer.data(kh + 342 * ncomps + c);
        const auto *kh_343 = buffer.data(kh + 343 * ncomps + c);
        const auto *kh_344 = buffer.data(kh + 344 * ncomps + c);
        const auto *kh_345 = buffer.data(kh + 345 * ncomps + c);
        const auto *kh_346 = buffer.data(kh + 346 * ncomps + c);
        const auto *kh_347 = buffer.data(kh + 347 * ncomps + c);
        const auto *kh_348 = buffer.data(kh + 348 * ncomps + c);
        const auto *kh_349 = buffer.data(kh + 349 * ncomps + c);
        const auto *kh_350 = buffer.data(kh + 350 * ncomps + c);
        const auto *kh_351 = buffer.data(kh + 351 * ncomps + c);
        const auto *kh_352 = buffer.data(kh + 352 * ncomps + c);
        const auto *kh_353 = buffer.data(kh + 353 * ncomps + c);
        const auto *kh_354 = buffer.data(kh + 354 * ncomps + c);
        const auto *kh_355 = buffer.data(kh + 355 * ncomps + c);
        const auto *kh_356 = buffer.data(kh + 356 * ncomps + c);
        const auto *kh_357 = buffer.data(kh + 357 * ncomps + c);
        const auto *kh_358 = buffer.data(kh + 358 * ncomps + c);
        const auto *kh_359 = buffer.data(kh + 359 * ncomps + c);
        const auto *kh_360 = buffer.data(kh + 360 * ncomps + c);
        const auto *kh_361 = buffer.data(kh + 361 * ncomps + c);
        const auto *kh_362 = buffer.data(kh + 362 * ncomps + c);
        const auto *kh_363 = buffer.data(kh + 363 * ncomps + c);
        const auto *kh_364 = buffer.data(kh + 364 * ncomps + c);
        const auto *kh_365 = buffer.data(kh + 365 * ncomps + c);
        const auto *kh_366 = buffer.data(kh + 366 * ncomps + c);
        const auto *kh_367 = buffer.data(kh + 367 * ncomps + c);
        const auto *kh_368 = buffer.data(kh + 368 * ncomps + c);
        const auto *kh_369 = buffer.data(kh + 369 * ncomps + c);
        const auto *kh_370 = buffer.data(kh + 370 * ncomps + c);
        const auto *kh_371 = buffer.data(kh + 371 * ncomps + c);
        const auto *kh_372 = buffer.data(kh + 372 * ncomps + c);
        const auto *kh_373 = buffer.data(kh + 373 * ncomps + c);
        const auto *kh_374 = buffer.data(kh + 374 * ncomps + c);
        const auto *kh_375 = buffer.data(kh + 375 * ncomps + c);
        const auto *kh_376 = buffer.data(kh + 376 * ncomps + c);
        const auto *kh_377 = buffer.data(kh + 377 * ncomps + c);
        const auto *kh_378 = buffer.data(kh + 378 * ncomps + c);
        const auto *kh_379 = buffer.data(kh + 379 * ncomps + c);
        const auto *kh_380 = buffer.data(kh + 380 * ncomps + c);
        const auto *kh_381 = buffer.data(kh + 381 * ncomps + c);
        const auto *kh_382 = buffer.data(kh + 382 * ncomps + c);
        const auto *kh_383 = buffer.data(kh + 383 * ncomps + c);
        const auto *kh_384 = buffer.data(kh + 384 * ncomps + c);
        const auto *kh_385 = buffer.data(kh + 385 * ncomps + c);
        const auto *kh_386 = buffer.data(kh + 386 * ncomps + c);
        const auto *kh_387 = buffer.data(kh + 387 * ncomps + c);
        const auto *kh_388 = buffer.data(kh + 388 * ncomps + c);
        const auto *kh_389 = buffer.data(kh + 389 * ncomps + c);
        const auto *kh_390 = buffer.data(kh + 390 * ncomps + c);
        const auto *kh_391 = buffer.data(kh + 391 * ncomps + c);
        const auto *kh_392 = buffer.data(kh + 392 * ncomps + c);
        const auto *kh_393 = buffer.data(kh + 393 * ncomps + c);
        const auto *kh_394 = buffer.data(kh + 394 * ncomps + c);
        const auto *kh_395 = buffer.data(kh + 395 * ncomps + c);
        const auto *kh_396 = buffer.data(kh + 396 * ncomps + c);
        const auto *kh_397 = buffer.data(kh + 397 * ncomps + c);
        const auto *kh_398 = buffer.data(kh + 398 * ncomps + c);
        const auto *kh_399 = buffer.data(kh + 399 * ncomps + c);
        const auto *kh_400 = buffer.data(kh + 400 * ncomps + c);
        const auto *kh_401 = buffer.data(kh + 401 * ncomps + c);
        const auto *kh_402 = buffer.data(kh + 402 * ncomps + c);
        const auto *kh_403 = buffer.data(kh + 403 * ncomps + c);
        const auto *kh_404 = buffer.data(kh + 404 * ncomps + c);
        const auto *kh_405 = buffer.data(kh + 405 * ncomps + c);
        const auto *kh_406 = buffer.data(kh + 406 * ncomps + c);
        const auto *kh_407 = buffer.data(kh + 407 * ncomps + c);
        const auto *kh_408 = buffer.data(kh + 408 * ncomps + c);
        const auto *kh_409 = buffer.data(kh + 409 * ncomps + c);
        const auto *kh_410 = buffer.data(kh + 410 * ncomps + c);
        const auto *kh_411 = buffer.data(kh + 411 * ncomps + c);
        const auto *kh_412 = buffer.data(kh + 412 * ncomps + c);
        const auto *kh_413 = buffer.data(kh + 413 * ncomps + c);
        const auto *kh_414 = buffer.data(kh + 414 * ncomps + c);
        const auto *kh_415 = buffer.data(kh + 415 * ncomps + c);
        const auto *kh_416 = buffer.data(kh + 416 * ncomps + c);
        const auto *kh_417 = buffer.data(kh + 417 * ncomps + c);
        const auto *kh_418 = buffer.data(kh + 418 * ncomps + c);
        const auto *kh_419 = buffer.data(kh + 419 * ncomps + c);
        const auto *kh_420 = buffer.data(kh + 420 * ncomps + c);
        const auto *kh_421 = buffer.data(kh + 421 * ncomps + c);
        const auto *kh_422 = buffer.data(kh + 422 * ncomps + c);
        const auto *kh_423 = buffer.data(kh + 423 * ncomps + c);
        const auto *kh_424 = buffer.data(kh + 424 * ncomps + c);
        const auto *kh_425 = buffer.data(kh + 425 * ncomps + c);
        const auto *kh_426 = buffer.data(kh + 426 * ncomps + c);
        const auto *kh_427 = buffer.data(kh + 427 * ncomps + c);
        const auto *kh_428 = buffer.data(kh + 428 * ncomps + c);
        const auto *kh_429 = buffer.data(kh + 429 * ncomps + c);
        const auto *kh_430 = buffer.data(kh + 430 * ncomps + c);
        const auto *kh_431 = buffer.data(kh + 431 * ncomps + c);
        const auto *kh_432 = buffer.data(kh + 432 * ncomps + c);
        const auto *kh_433 = buffer.data(kh + 433 * ncomps + c);
        const auto *kh_434 = buffer.data(kh + 434 * ncomps + c);
        const auto *kh_435 = buffer.data(kh + 435 * ncomps + c);
        const auto *kh_436 = buffer.data(kh + 436 * ncomps + c);
        const auto *kh_437 = buffer.data(kh + 437 * ncomps + c);
        const auto *kh_438 = buffer.data(kh + 438 * ncomps + c);
        const auto *kh_439 = buffer.data(kh + 439 * ncomps + c);
        const auto *kh_456 = buffer.data(kh + 456 * ncomps + c);
        const auto *kh_457 = buffer.data(kh + 457 * ncomps + c);
        const auto *kh_458 = buffer.data(kh + 458 * ncomps + c);
        const auto *kh_459 = buffer.data(kh + 459 * ncomps + c);
        const auto *kh_460 = buffer.data(kh + 460 * ncomps + c);
        const auto *kh_461 = buffer.data(kh + 461 * ncomps + c);
        const auto *kh_477 = buffer.data(kh + 477 * ncomps + c);
        const auto *kh_478 = buffer.data(kh + 478 * ncomps + c);
        const auto *kh_479 = buffer.data(kh + 479 * ncomps + c);
        const auto *kh_480 = buffer.data(kh + 480 * ncomps + c);
        const auto *kh_481 = buffer.data(kh + 481 * ncomps + c);
        const auto *kh_482 = buffer.data(kh + 482 * ncomps + c);
        const auto *kh_498 = buffer.data(kh + 498 * ncomps + c);
        const auto *kh_499 = buffer.data(kh + 499 * ncomps + c);
        const auto *kh_500 = buffer.data(kh + 500 * ncomps + c);
        const auto *kh_501 = buffer.data(kh + 501 * ncomps + c);
        const auto *kh_502 = buffer.data(kh + 502 * ncomps + c);
        const auto *kh_503 = buffer.data(kh + 503 * ncomps + c);
        const auto *kh_519 = buffer.data(kh + 519 * ncomps + c);
        const auto *kh_520 = buffer.data(kh + 520 * ncomps + c);
        const auto *kh_521 = buffer.data(kh + 521 * ncomps + c);
        const auto *kh_522 = buffer.data(kh + 522 * ncomps + c);
        const auto *kh_523 = buffer.data(kh + 523 * ncomps + c);
        const auto *kh_524 = buffer.data(kh + 524 * ncomps + c);
        const auto *kh_540 = buffer.data(kh + 540 * ncomps + c);
        const auto *kh_541 = buffer.data(kh + 541 * ncomps + c);
        const auto *kh_542 = buffer.data(kh + 542 * ncomps + c);
        const auto *kh_543 = buffer.data(kh + 543 * ncomps + c);
        const auto *kh_544 = buffer.data(kh + 544 * ncomps + c);
        const auto *kh_545 = buffer.data(kh + 545 * ncomps + c);
        const auto *kh_566 = buffer.data(kh + 566 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, ih_330, ih_331, ih_332, \
                         ih_333, ih_334, kh_330, kh_331, kh_332, kh_333, \
                         kh_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * ih_330[k]
                       + kh_330[k];

            t_436[k] = ab_x[k] * ih_331[k]
                       + kh_331[k];

            t_437[k] = ab_x[k] * ih_332[k]
                       + kh_332[k];

            t_438[k] = ab_x[k] * ih_333[k]
                       + kh_333[k];

            t_439[k] = ab_x[k] * ih_334[k]
                       + kh_334[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, ab_x, ab_y, ih_330, ih_331, ih_332, \
                         ih_335, kh_335, kh_456, kh_457, kh_458 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * ih_335[k]
                       + kh_335[k];

            t_441[k] = ab_y[k] * ih_330[k]
                       + kh_456[k];

            t_442[k] = ab_y[k] * ih_331[k]
                       + kh_457[k];

            t_443[k] = ab_y[k] * ih_332[k]
                       + kh_458[k];
        }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, ab_y, ab_z, ih_333, ih_334, ih_335, \
                         kh_459, kh_460, kh_461, kh_482 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_444[k] = ab_y[k] * ih_333[k]
                       + kh_459[k];

            t_445[k] = ab_y[k] * ih_334[k]
                       + kh_460[k];

            t_446[k] = ab_y[k] * ih_335[k]
                       + kh_461[k];

            t_447[k] = ab_z[k] * ih_335[k]
                       + kh_482[k];
        }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, ab_x, ih_336, ih_337, ih_338, \
                         ih_339, ih_340, kh_336, kh_337, kh_338, kh_339, \
                         kh_340 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_448[k] = ab_x[k] * ih_336[k]
                       + kh_336[k];

            t_449[k] = ab_x[k] * ih_337[k]
                       + kh_337[k];

            t_450[k] = ab_x[k] * ih_338[k]
                       + kh_338[k];

            t_451[k] = ab_x[k] * ih_339[k]
                       + kh_339[k];

            t_452[k] = ab_x[k] * ih_340[k]
                       + kh_340[k];
        }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, ab_x, ih_341, ih_342, ih_343, \
                         ih_344, ih_345, kh_341, kh_342, kh_343, kh_344, \
                         kh_345 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_453[k] = ab_x[k] * ih_341[k]
                       + kh_341[k];

            t_454[k] = ab_x[k] * ih_342[k]
                       + kh_342[k];

            t_455[k] = ab_x[k] * ih_343[k]
                       + kh_343[k];

            t_456[k] = ab_x[k] * ih_344[k]
                       + kh_344[k];

            t_457[k] = ab_x[k] * ih_345[k]
                       + kh_345[k];
        }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, ab_x, ih_346, ih_347, ih_348, \
                         ih_349, ih_350, kh_346, kh_347, kh_348, kh_349, \
                         kh_350 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_458[k] = ab_x[k] * ih_346[k]
                       + kh_346[k];

            t_459[k] = ab_x[k] * ih_347[k]
                       + kh_347[k];

            t_460[k] = ab_x[k] * ih_348[k]
                       + kh_348[k];

            t_461[k] = ab_x[k] * ih_349[k]
                       + kh_349[k];

            t_462[k] = ab_x[k] * ih_350[k]
                       + kh_350[k];
        }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, ab_x, ih_351, ih_352, ih_353, \
                         ih_354, ih_355, kh_351, kh_352, kh_353, kh_354, \
                         kh_355 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_463[k] = ab_x[k] * ih_351[k]
                       + kh_351[k];

            t_464[k] = ab_x[k] * ih_352[k]
                       + kh_352[k];

            t_465[k] = ab_x[k] * ih_353[k]
                       + kh_353[k];

            t_466[k] = ab_x[k] * ih_354[k]
                       + kh_354[k];

            t_467[k] = ab_x[k] * ih_355[k]
                       + kh_355[k];
        }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, ab_x, ab_y, ih_351, ih_352, ih_353, \
                         ih_356, kh_356, kh_477, kh_478, kh_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_468[k] = ab_x[k] * ih_356[k]
                       + kh_356[k];

            t_469[k] = ab_y[k] * ih_351[k]
                       + kh_477[k];

            t_470[k] = ab_y[k] * ih_352[k]
                       + kh_478[k];

            t_471[k] = ab_y[k] * ih_353[k]
                       + kh_479[k];
        }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, ab_y, ab_z, ih_354, ih_355, ih_356, \
                         kh_480, kh_481, kh_482, kh_503 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_472[k] = ab_y[k] * ih_354[k]
                       + kh_480[k];

            t_473[k] = ab_y[k] * ih_355[k]
                       + kh_481[k];

            t_474[k] = ab_y[k] * ih_356[k]
                       + kh_482[k];

            t_475[k] = ab_z[k] * ih_356[k]
                       + kh_503[k];
        }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, ab_x, ih_357, ih_358, ih_359, \
                         ih_360, ih_361, kh_357, kh_358, kh_359, kh_360, \
                         kh_361 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_476[k] = ab_x[k] * ih_357[k]
                       + kh_357[k];

            t_477[k] = ab_x[k] * ih_358[k]
                       + kh_358[k];

            t_478[k] = ab_x[k] * ih_359[k]
                       + kh_359[k];

            t_479[k] = ab_x[k] * ih_360[k]
                       + kh_360[k];

            t_480[k] = ab_x[k] * ih_361[k]
                       + kh_361[k];
        }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, ab_x, ih_362, ih_363, ih_364, \
                         ih_365, ih_366, kh_362, kh_363, kh_364, kh_365, \
                         kh_366 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_481[k] = ab_x[k] * ih_362[k]
                       + kh_362[k];

            t_482[k] = ab_x[k] * ih_363[k]
                       + kh_363[k];

            t_483[k] = ab_x[k] * ih_364[k]
                       + kh_364[k];

            t_484[k] = ab_x[k] * ih_365[k]
                       + kh_365[k];

            t_485[k] = ab_x[k] * ih_366[k]
                       + kh_366[k];
        }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, ab_x, ih_367, ih_368, ih_369, \
                         ih_370, ih_371, kh_367, kh_368, kh_369, kh_370, \
                         kh_371 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_486[k] = ab_x[k] * ih_367[k]
                       + kh_367[k];

            t_487[k] = ab_x[k] * ih_368[k]
                       + kh_368[k];

            t_488[k] = ab_x[k] * ih_369[k]
                       + kh_369[k];

            t_489[k] = ab_x[k] * ih_370[k]
                       + kh_370[k];

            t_490[k] = ab_x[k] * ih_371[k]
                       + kh_371[k];
        }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, ab_x, ih_372, ih_373, ih_374, \
                         ih_375, ih_376, kh_372, kh_373, kh_374, kh_375, \
                         kh_376 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_491[k] = ab_x[k] * ih_372[k]
                       + kh_372[k];

            t_492[k] = ab_x[k] * ih_373[k]
                       + kh_373[k];

            t_493[k] = ab_x[k] * ih_374[k]
                       + kh_374[k];

            t_494[k] = ab_x[k] * ih_375[k]
                       + kh_375[k];

            t_495[k] = ab_x[k] * ih_376[k]
                       + kh_376[k];
        }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, ab_x, ab_y, ih_372, ih_373, ih_374, \
                         ih_377, kh_377, kh_498, kh_499, kh_500 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_496[k] = ab_x[k] * ih_377[k]
                       + kh_377[k];

            t_497[k] = ab_y[k] * ih_372[k]
                       + kh_498[k];

            t_498[k] = ab_y[k] * ih_373[k]
                       + kh_499[k];

            t_499[k] = ab_y[k] * ih_374[k]
                       + kh_500[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, ab_y, ab_z, ih_375, ih_376, ih_377, \
                         kh_501, kh_502, kh_503, kh_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_y[k] * ih_375[k]
                       + kh_501[k];

            t_501[k] = ab_y[k] * ih_376[k]
                       + kh_502[k];

            t_502[k] = ab_y[k] * ih_377[k]
                       + kh_503[k];

            t_503[k] = ab_z[k] * ih_377[k]
                       + kh_524[k];
        }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, ab_x, ih_378, ih_379, ih_380, \
                         ih_381, ih_382, kh_378, kh_379, kh_380, kh_381, \
                         kh_382 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_504[k] = ab_x[k] * ih_378[k]
                       + kh_378[k];

            t_505[k] = ab_x[k] * ih_379[k]
                       + kh_379[k];

            t_506[k] = ab_x[k] * ih_380[k]
                       + kh_380[k];

            t_507[k] = ab_x[k] * ih_381[k]
                       + kh_381[k];

            t_508[k] = ab_x[k] * ih_382[k]
                       + kh_382[k];
        }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, ab_x, ih_383, ih_384, ih_385, \
                         ih_386, ih_387, kh_383, kh_384, kh_385, kh_386, \
                         kh_387 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_509[k] = ab_x[k] * ih_383[k]
                       + kh_383[k];

            t_510[k] = ab_x[k] * ih_384[k]
                       + kh_384[k];

            t_511[k] = ab_x[k] * ih_385[k]
                       + kh_385[k];

            t_512[k] = ab_x[k] * ih_386[k]
                       + kh_386[k];

            t_513[k] = ab_x[k] * ih_387[k]
                       + kh_387[k];
        }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, ab_x, ih_388, ih_389, ih_390, \
                         ih_391, ih_392, kh_388, kh_389, kh_390, kh_391, \
                         kh_392 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_514[k] = ab_x[k] * ih_388[k]
                       + kh_388[k];

            t_515[k] = ab_x[k] * ih_389[k]
                       + kh_389[k];

            t_516[k] = ab_x[k] * ih_390[k]
                       + kh_390[k];

            t_517[k] = ab_x[k] * ih_391[k]
                       + kh_391[k];

            t_518[k] = ab_x[k] * ih_392[k]
                       + kh_392[k];
        }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, ab_x, ih_393, ih_394, ih_395, \
                         ih_396, ih_397, kh_393, kh_394, kh_395, kh_396, \
                         kh_397 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_519[k] = ab_x[k] * ih_393[k]
                       + kh_393[k];

            t_520[k] = ab_x[k] * ih_394[k]
                       + kh_394[k];

            t_521[k] = ab_x[k] * ih_395[k]
                       + kh_395[k];

            t_522[k] = ab_x[k] * ih_396[k]
                       + kh_396[k];

            t_523[k] = ab_x[k] * ih_397[k]
                       + kh_397[k];
        }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, ab_x, ab_y, ih_393, ih_394, ih_395, \
                         ih_398, kh_398, kh_519, kh_520, kh_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_524[k] = ab_x[k] * ih_398[k]
                       + kh_398[k];

            t_525[k] = ab_y[k] * ih_393[k]
                       + kh_519[k];

            t_526[k] = ab_y[k] * ih_394[k]
                       + kh_520[k];

            t_527[k] = ab_y[k] * ih_395[k]
                       + kh_521[k];
        }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, ab_y, ab_z, ih_396, ih_397, ih_398, \
                         kh_522, kh_523, kh_524, kh_545 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_528[k] = ab_y[k] * ih_396[k]
                       + kh_522[k];

            t_529[k] = ab_y[k] * ih_397[k]
                       + kh_523[k];

            t_530[k] = ab_y[k] * ih_398[k]
                       + kh_524[k];

            t_531[k] = ab_z[k] * ih_398[k]
                       + kh_545[k];
        }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, ab_x, ih_399, ih_400, ih_401, \
                         ih_402, ih_403, kh_399, kh_400, kh_401, kh_402, \
                         kh_403 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_532[k] = ab_x[k] * ih_399[k]
                       + kh_399[k];

            t_533[k] = ab_x[k] * ih_400[k]
                       + kh_400[k];

            t_534[k] = ab_x[k] * ih_401[k]
                       + kh_401[k];

            t_535[k] = ab_x[k] * ih_402[k]
                       + kh_402[k];

            t_536[k] = ab_x[k] * ih_403[k]
                       + kh_403[k];
        }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, ab_x, ih_404, ih_405, ih_406, \
                         ih_407, ih_408, kh_404, kh_405, kh_406, kh_407, \
                         kh_408 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_537[k] = ab_x[k] * ih_404[k]
                       + kh_404[k];

            t_538[k] = ab_x[k] * ih_405[k]
                       + kh_405[k];

            t_539[k] = ab_x[k] * ih_406[k]
                       + kh_406[k];

            t_540[k] = ab_x[k] * ih_407[k]
                       + kh_407[k];

            t_541[k] = ab_x[k] * ih_408[k]
                       + kh_408[k];
        }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, ab_x, ih_409, ih_410, ih_411, \
                         ih_412, ih_413, kh_409, kh_410, kh_411, kh_412, \
                         kh_413 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_542[k] = ab_x[k] * ih_409[k]
                       + kh_409[k];

            t_543[k] = ab_x[k] * ih_410[k]
                       + kh_410[k];

            t_544[k] = ab_x[k] * ih_411[k]
                       + kh_411[k];

            t_545[k] = ab_x[k] * ih_412[k]
                       + kh_412[k];

            t_546[k] = ab_x[k] * ih_413[k]
                       + kh_413[k];
        }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, ab_x, ih_414, ih_415, ih_416, \
                         ih_417, ih_418, kh_414, kh_415, kh_416, kh_417, \
                         kh_418 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_547[k] = ab_x[k] * ih_414[k]
                       + kh_414[k];

            t_548[k] = ab_x[k] * ih_415[k]
                       + kh_415[k];

            t_549[k] = ab_x[k] * ih_416[k]
                       + kh_416[k];

            t_550[k] = ab_x[k] * ih_417[k]
                       + kh_417[k];

            t_551[k] = ab_x[k] * ih_418[k]
                       + kh_418[k];
        }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, ab_x, ab_y, ih_414, ih_415, ih_416, \
                         ih_419, kh_419, kh_540, kh_541, kh_542 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_552[k] = ab_x[k] * ih_419[k]
                       + kh_419[k];

            t_553[k] = ab_y[k] * ih_414[k]
                       + kh_540[k];

            t_554[k] = ab_y[k] * ih_415[k]
                       + kh_541[k];

            t_555[k] = ab_y[k] * ih_416[k]
                       + kh_542[k];
        }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, ab_y, ab_z, ih_417, ih_418, ih_419, \
                         kh_543, kh_544, kh_545, kh_566 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_556[k] = ab_y[k] * ih_417[k]
                       + kh_543[k];

            t_557[k] = ab_y[k] * ih_418[k]
                       + kh_544[k];

            t_558[k] = ab_y[k] * ih_419[k]
                       + kh_545[k];

            t_559[k] = ab_z[k] * ih_419[k]
                       + kh_566[k];
        }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_x, ih_420, ih_421, ih_422, \
                         ih_423, ih_424, kh_420, kh_421, kh_422, kh_423, \
                         kh_424 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_560[k] = ab_x[k] * ih_420[k]
                       + kh_420[k];

            t_561[k] = ab_x[k] * ih_421[k]
                       + kh_421[k];

            t_562[k] = ab_x[k] * ih_422[k]
                       + kh_422[k];

            t_563[k] = ab_x[k] * ih_423[k]
                       + kh_423[k];

            t_564[k] = ab_x[k] * ih_424[k]
                       + kh_424[k];
        }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_x, ih_425, ih_426, ih_427, \
                         ih_428, ih_429, kh_425, kh_426, kh_427, kh_428, \
                         kh_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_565[k] = ab_x[k] * ih_425[k]
                       + kh_425[k];

            t_566[k] = ab_x[k] * ih_426[k]
                       + kh_426[k];

            t_567[k] = ab_x[k] * ih_427[k]
                       + kh_427[k];

            t_568[k] = ab_x[k] * ih_428[k]
                       + kh_428[k];

            t_569[k] = ab_x[k] * ih_429[k]
                       + kh_429[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_x, ih_430, ih_431, ih_432, \
                         ih_433, ih_434, kh_430, kh_431, kh_432, kh_433, \
                         kh_434 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = ab_x[k] * ih_430[k]
                       + kh_430[k];

            t_571[k] = ab_x[k] * ih_431[k]
                       + kh_431[k];

            t_572[k] = ab_x[k] * ih_432[k]
                       + kh_432[k];

            t_573[k] = ab_x[k] * ih_433[k]
                       + kh_433[k];

            t_574[k] = ab_x[k] * ih_434[k]
                       + kh_434[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_x, ih_435, ih_436, ih_437, \
                         ih_438, ih_439, kh_435, kh_436, kh_437, kh_438, \
                         kh_439 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = ab_x[k] * ih_435[k]
                       + kh_435[k];

            t_576[k] = ab_x[k] * ih_436[k]
                       + kh_436[k];

            t_577[k] = ab_x[k] * ih_437[k]
                       + kh_437[k];

            t_578[k] = ab_x[k] * ih_438[k]
                       + kh_438[k];

            t_579[k] = ab_x[k] * ih_439[k]
                       + kh_439[k];
        }
    }
}

static auto
compute_hrr_ii_out_of_first_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ih, const size_t kh,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_580 = buffer.data(target + 580 * ncomps + c);
        auto *t_581 = buffer.data(target + 581 * ncomps + c);
        auto *t_582 = buffer.data(target + 582 * ncomps + c);
        auto *t_583 = buffer.data(target + 583 * ncomps + c);
        auto *t_584 = buffer.data(target + 584 * ncomps + c);
        auto *t_585 = buffer.data(target + 585 * ncomps + c);
        auto *t_586 = buffer.data(target + 586 * ncomps + c);
        auto *t_587 = buffer.data(target + 587 * ncomps + c);
        auto *t_588 = buffer.data(target + 588 * ncomps + c);
        auto *t_589 = buffer.data(target + 589 * ncomps + c);
        auto *t_590 = buffer.data(target + 590 * ncomps + c);
        auto *t_591 = buffer.data(target + 591 * ncomps + c);
        auto *t_592 = buffer.data(target + 592 * ncomps + c);
        auto *t_593 = buffer.data(target + 593 * ncomps + c);
        auto *t_594 = buffer.data(target + 594 * ncomps + c);
        auto *t_595 = buffer.data(target + 595 * ncomps + c);
        auto *t_596 = buffer.data(target + 596 * ncomps + c);
        auto *t_597 = buffer.data(target + 597 * ncomps + c);
        auto *t_598 = buffer.data(target + 598 * ncomps + c);
        auto *t_599 = buffer.data(target + 599 * ncomps + c);
        auto *t_600 = buffer.data(target + 600 * ncomps + c);
        auto *t_601 = buffer.data(target + 601 * ncomps + c);
        auto *t_602 = buffer.data(target + 602 * ncomps + c);
        auto *t_603 = buffer.data(target + 603 * ncomps + c);
        auto *t_604 = buffer.data(target + 604 * ncomps + c);
        auto *t_605 = buffer.data(target + 605 * ncomps + c);
        auto *t_606 = buffer.data(target + 606 * ncomps + c);
        auto *t_607 = buffer.data(target + 607 * ncomps + c);
        auto *t_608 = buffer.data(target + 608 * ncomps + c);
        auto *t_609 = buffer.data(target + 609 * ncomps + c);
        auto *t_610 = buffer.data(target + 610 * ncomps + c);
        auto *t_611 = buffer.data(target + 611 * ncomps + c);
        auto *t_612 = buffer.data(target + 612 * ncomps + c);
        auto *t_613 = buffer.data(target + 613 * ncomps + c);
        auto *t_614 = buffer.data(target + 614 * ncomps + c);
        auto *t_615 = buffer.data(target + 615 * ncomps + c);
        auto *t_616 = buffer.data(target + 616 * ncomps + c);
        auto *t_617 = buffer.data(target + 617 * ncomps + c);
        auto *t_618 = buffer.data(target + 618 * ncomps + c);
        auto *t_619 = buffer.data(target + 619 * ncomps + c);
        auto *t_620 = buffer.data(target + 620 * ncomps + c);
        auto *t_621 = buffer.data(target + 621 * ncomps + c);
        auto *t_622 = buffer.data(target + 622 * ncomps + c);
        auto *t_623 = buffer.data(target + 623 * ncomps + c);
        auto *t_624 = buffer.data(target + 624 * ncomps + c);
        auto *t_625 = buffer.data(target + 625 * ncomps + c);
        auto *t_626 = buffer.data(target + 626 * ncomps + c);
        auto *t_627 = buffer.data(target + 627 * ncomps + c);
        auto *t_628 = buffer.data(target + 628 * ncomps + c);
        auto *t_629 = buffer.data(target + 629 * ncomps + c);
        auto *t_630 = buffer.data(target + 630 * ncomps + c);
        auto *t_631 = buffer.data(target + 631 * ncomps + c);
        auto *t_632 = buffer.data(target + 632 * ncomps + c);
        auto *t_633 = buffer.data(target + 633 * ncomps + c);
        auto *t_634 = buffer.data(target + 634 * ncomps + c);
        auto *t_635 = buffer.data(target + 635 * ncomps + c);
        auto *t_636 = buffer.data(target + 636 * ncomps + c);
        auto *t_637 = buffer.data(target + 637 * ncomps + c);
        auto *t_638 = buffer.data(target + 638 * ncomps + c);
        auto *t_639 = buffer.data(target + 639 * ncomps + c);
        auto *t_640 = buffer.data(target + 640 * ncomps + c);
        auto *t_641 = buffer.data(target + 641 * ncomps + c);
        auto *t_642 = buffer.data(target + 642 * ncomps + c);
        auto *t_643 = buffer.data(target + 643 * ncomps + c);
        auto *t_644 = buffer.data(target + 644 * ncomps + c);
        auto *t_645 = buffer.data(target + 645 * ncomps + c);
        auto *t_646 = buffer.data(target + 646 * ncomps + c);
        auto *t_647 = buffer.data(target + 647 * ncomps + c);
        auto *t_648 = buffer.data(target + 648 * ncomps + c);
        auto *t_649 = buffer.data(target + 649 * ncomps + c);
        auto *t_650 = buffer.data(target + 650 * ncomps + c);
        auto *t_651 = buffer.data(target + 651 * ncomps + c);
        auto *t_652 = buffer.data(target + 652 * ncomps + c);
        auto *t_653 = buffer.data(target + 653 * ncomps + c);
        auto *t_654 = buffer.data(target + 654 * ncomps + c);
        auto *t_655 = buffer.data(target + 655 * ncomps + c);
        auto *t_656 = buffer.data(target + 656 * ncomps + c);
        auto *t_657 = buffer.data(target + 657 * ncomps + c);
        auto *t_658 = buffer.data(target + 658 * ncomps + c);
        auto *t_659 = buffer.data(target + 659 * ncomps + c);
        auto *t_660 = buffer.data(target + 660 * ncomps + c);
        auto *t_661 = buffer.data(target + 661 * ncomps + c);
        auto *t_662 = buffer.data(target + 662 * ncomps + c);
        auto *t_663 = buffer.data(target + 663 * ncomps + c);
        auto *t_664 = buffer.data(target + 664 * ncomps + c);
        auto *t_665 = buffer.data(target + 665 * ncomps + c);
        auto *t_666 = buffer.data(target + 666 * ncomps + c);
        auto *t_667 = buffer.data(target + 667 * ncomps + c);
        auto *t_668 = buffer.data(target + 668 * ncomps + c);
        auto *t_669 = buffer.data(target + 669 * ncomps + c);
        auto *t_670 = buffer.data(target + 670 * ncomps + c);
        auto *t_671 = buffer.data(target + 671 * ncomps + c);
        auto *t_672 = buffer.data(target + 672 * ncomps + c);
        auto *t_673 = buffer.data(target + 673 * ncomps + c);
        auto *t_674 = buffer.data(target + 674 * ncomps + c);
        auto *t_675 = buffer.data(target + 675 * ncomps + c);
        auto *t_676 = buffer.data(target + 676 * ncomps + c);
        auto *t_677 = buffer.data(target + 677 * ncomps + c);
        auto *t_678 = buffer.data(target + 678 * ncomps + c);
        auto *t_679 = buffer.data(target + 679 * ncomps + c);
        auto *t_680 = buffer.data(target + 680 * ncomps + c);
        auto *t_681 = buffer.data(target + 681 * ncomps + c);
        auto *t_682 = buffer.data(target + 682 * ncomps + c);
        auto *t_683 = buffer.data(target + 683 * ncomps + c);
        auto *t_684 = buffer.data(target + 684 * ncomps + c);
        auto *t_685 = buffer.data(target + 685 * ncomps + c);
        auto *t_686 = buffer.data(target + 686 * ncomps + c);
        auto *t_687 = buffer.data(target + 687 * ncomps + c);
        auto *t_688 = buffer.data(target + 688 * ncomps + c);
        auto *t_689 = buffer.data(target + 689 * ncomps + c);
        auto *t_690 = buffer.data(target + 690 * ncomps + c);
        auto *t_691 = buffer.data(target + 691 * ncomps + c);
        auto *t_692 = buffer.data(target + 692 * ncomps + c);
        auto *t_693 = buffer.data(target + 693 * ncomps + c);
        auto *t_694 = buffer.data(target + 694 * ncomps + c);
        auto *t_695 = buffer.data(target + 695 * ncomps + c);
        auto *t_696 = buffer.data(target + 696 * ncomps + c);
        auto *t_697 = buffer.data(target + 697 * ncomps + c);
        auto *t_698 = buffer.data(target + 698 * ncomps + c);
        auto *t_699 = buffer.data(target + 699 * ncomps + c);
        auto *t_700 = buffer.data(target + 700 * ncomps + c);
        auto *t_701 = buffer.data(target + 701 * ncomps + c);
        auto *t_702 = buffer.data(target + 702 * ncomps + c);
        auto *t_703 = buffer.data(target + 703 * ncomps + c);
        auto *t_704 = buffer.data(target + 704 * ncomps + c);
        auto *t_705 = buffer.data(target + 705 * ncomps + c);
        auto *t_706 = buffer.data(target + 706 * ncomps + c);
        auto *t_707 = buffer.data(target + 707 * ncomps + c);
        auto *t_708 = buffer.data(target + 708 * ncomps + c);
        auto *t_709 = buffer.data(target + 709 * ncomps + c);
        auto *t_710 = buffer.data(target + 710 * ncomps + c);
        auto *t_711 = buffer.data(target + 711 * ncomps + c);
        auto *t_712 = buffer.data(target + 712 * ncomps + c);
        auto *t_713 = buffer.data(target + 713 * ncomps + c);
        auto *t_714 = buffer.data(target + 714 * ncomps + c);
        auto *t_715 = buffer.data(target + 715 * ncomps + c);
        auto *t_716 = buffer.data(target + 716 * ncomps + c);
        auto *t_717 = buffer.data(target + 717 * ncomps + c);
        auto *t_718 = buffer.data(target + 718 * ncomps + c);
        auto *t_719 = buffer.data(target + 719 * ncomps + c);
        auto *t_720 = buffer.data(target + 720 * ncomps + c);
        auto *t_721 = buffer.data(target + 721 * ncomps + c);
        auto *t_722 = buffer.data(target + 722 * ncomps + c);
        auto *t_723 = buffer.data(target + 723 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ih_435 = buffer.data(ih + 435 * ncomps + c);
        const auto *ih_436 = buffer.data(ih + 436 * ncomps + c);
        const auto *ih_437 = buffer.data(ih + 437 * ncomps + c);
        const auto *ih_438 = buffer.data(ih + 438 * ncomps + c);
        const auto *ih_439 = buffer.data(ih + 439 * ncomps + c);
        const auto *ih_440 = buffer.data(ih + 440 * ncomps + c);
        const auto *ih_441 = buffer.data(ih + 441 * ncomps + c);
        const auto *ih_442 = buffer.data(ih + 442 * ncomps + c);
        const auto *ih_443 = buffer.data(ih + 443 * ncomps + c);
        const auto *ih_444 = buffer.data(ih + 444 * ncomps + c);
        const auto *ih_445 = buffer.data(ih + 445 * ncomps + c);
        const auto *ih_446 = buffer.data(ih + 446 * ncomps + c);
        const auto *ih_447 = buffer.data(ih + 447 * ncomps + c);
        const auto *ih_448 = buffer.data(ih + 448 * ncomps + c);
        const auto *ih_449 = buffer.data(ih + 449 * ncomps + c);
        const auto *ih_450 = buffer.data(ih + 450 * ncomps + c);
        const auto *ih_451 = buffer.data(ih + 451 * ncomps + c);
        const auto *ih_452 = buffer.data(ih + 452 * ncomps + c);
        const auto *ih_453 = buffer.data(ih + 453 * ncomps + c);
        const auto *ih_454 = buffer.data(ih + 454 * ncomps + c);
        const auto *ih_455 = buffer.data(ih + 455 * ncomps + c);
        const auto *ih_456 = buffer.data(ih + 456 * ncomps + c);
        const auto *ih_457 = buffer.data(ih + 457 * ncomps + c);
        const auto *ih_458 = buffer.data(ih + 458 * ncomps + c);
        const auto *ih_459 = buffer.data(ih + 459 * ncomps + c);
        const auto *ih_460 = buffer.data(ih + 460 * ncomps + c);
        const auto *ih_461 = buffer.data(ih + 461 * ncomps + c);
        const auto *ih_462 = buffer.data(ih + 462 * ncomps + c);
        const auto *ih_463 = buffer.data(ih + 463 * ncomps + c);
        const auto *ih_464 = buffer.data(ih + 464 * ncomps + c);
        const auto *ih_465 = buffer.data(ih + 465 * ncomps + c);
        const auto *ih_466 = buffer.data(ih + 466 * ncomps + c);
        const auto *ih_467 = buffer.data(ih + 467 * ncomps + c);
        const auto *ih_468 = buffer.data(ih + 468 * ncomps + c);
        const auto *ih_469 = buffer.data(ih + 469 * ncomps + c);
        const auto *ih_470 = buffer.data(ih + 470 * ncomps + c);
        const auto *ih_471 = buffer.data(ih + 471 * ncomps + c);
        const auto *ih_472 = buffer.data(ih + 472 * ncomps + c);
        const auto *ih_473 = buffer.data(ih + 473 * ncomps + c);
        const auto *ih_474 = buffer.data(ih + 474 * ncomps + c);
        const auto *ih_475 = buffer.data(ih + 475 * ncomps + c);
        const auto *ih_476 = buffer.data(ih + 476 * ncomps + c);
        const auto *ih_477 = buffer.data(ih + 477 * ncomps + c);
        const auto *ih_478 = buffer.data(ih + 478 * ncomps + c);
        const auto *ih_479 = buffer.data(ih + 479 * ncomps + c);
        const auto *ih_480 = buffer.data(ih + 480 * ncomps + c);
        const auto *ih_481 = buffer.data(ih + 481 * ncomps + c);
        const auto *ih_482 = buffer.data(ih + 482 * ncomps + c);
        const auto *ih_483 = buffer.data(ih + 483 * ncomps + c);
        const auto *ih_484 = buffer.data(ih + 484 * ncomps + c);
        const auto *ih_485 = buffer.data(ih + 485 * ncomps + c);
        const auto *ih_486 = buffer.data(ih + 486 * ncomps + c);
        const auto *ih_487 = buffer.data(ih + 487 * ncomps + c);
        const auto *ih_488 = buffer.data(ih + 488 * ncomps + c);
        const auto *ih_489 = buffer.data(ih + 489 * ncomps + c);
        const auto *ih_490 = buffer.data(ih + 490 * ncomps + c);
        const auto *ih_491 = buffer.data(ih + 491 * ncomps + c);
        const auto *ih_492 = buffer.data(ih + 492 * ncomps + c);
        const auto *ih_493 = buffer.data(ih + 493 * ncomps + c);
        const auto *ih_494 = buffer.data(ih + 494 * ncomps + c);
        const auto *ih_495 = buffer.data(ih + 495 * ncomps + c);
        const auto *ih_496 = buffer.data(ih + 496 * ncomps + c);
        const auto *ih_497 = buffer.data(ih + 497 * ncomps + c);
        const auto *ih_498 = buffer.data(ih + 498 * ncomps + c);
        const auto *ih_499 = buffer.data(ih + 499 * ncomps + c);
        const auto *ih_500 = buffer.data(ih + 500 * ncomps + c);
        const auto *ih_501 = buffer.data(ih + 501 * ncomps + c);
        const auto *ih_502 = buffer.data(ih + 502 * ncomps + c);
        const auto *ih_503 = buffer.data(ih + 503 * ncomps + c);
        const auto *ih_504 = buffer.data(ih + 504 * ncomps + c);
        const auto *ih_505 = buffer.data(ih + 505 * ncomps + c);
        const auto *ih_506 = buffer.data(ih + 506 * ncomps + c);
        const auto *ih_507 = buffer.data(ih + 507 * ncomps + c);
        const auto *ih_508 = buffer.data(ih + 508 * ncomps + c);
        const auto *ih_509 = buffer.data(ih + 509 * ncomps + c);
        const auto *ih_510 = buffer.data(ih + 510 * ncomps + c);
        const auto *ih_511 = buffer.data(ih + 511 * ncomps + c);
        const auto *ih_512 = buffer.data(ih + 512 * ncomps + c);
        const auto *ih_513 = buffer.data(ih + 513 * ncomps + c);
        const auto *ih_514 = buffer.data(ih + 514 * ncomps + c);
        const auto *ih_515 = buffer.data(ih + 515 * ncomps + c);
        const auto *ih_516 = buffer.data(ih + 516 * ncomps + c);
        const auto *ih_517 = buffer.data(ih + 517 * ncomps + c);
        const auto *ih_518 = buffer.data(ih + 518 * ncomps + c);
        const auto *ih_519 = buffer.data(ih + 519 * ncomps + c);
        const auto *ih_520 = buffer.data(ih + 520 * ncomps + c);
        const auto *ih_521 = buffer.data(ih + 521 * ncomps + c);
        const auto *ih_522 = buffer.data(ih + 522 * ncomps + c);
        const auto *ih_523 = buffer.data(ih + 523 * ncomps + c);
        const auto *ih_524 = buffer.data(ih + 524 * ncomps + c);
        const auto *ih_525 = buffer.data(ih + 525 * ncomps + c);
        const auto *ih_526 = buffer.data(ih + 526 * ncomps + c);
        const auto *ih_527 = buffer.data(ih + 527 * ncomps + c);
        const auto *ih_528 = buffer.data(ih + 528 * ncomps + c);
        const auto *ih_529 = buffer.data(ih + 529 * ncomps + c);
        const auto *ih_530 = buffer.data(ih + 530 * ncomps + c);
        const auto *ih_531 = buffer.data(ih + 531 * ncomps + c);
        const auto *ih_532 = buffer.data(ih + 532 * ncomps + c);
        const auto *ih_533 = buffer.data(ih + 533 * ncomps + c);
        const auto *ih_534 = buffer.data(ih + 534 * ncomps + c);
        const auto *ih_535 = buffer.data(ih + 535 * ncomps + c);
        const auto *ih_536 = buffer.data(ih + 536 * ncomps + c);
        const auto *ih_537 = buffer.data(ih + 537 * ncomps + c);
        const auto *ih_538 = buffer.data(ih + 538 * ncomps + c);
        const auto *ih_539 = buffer.data(ih + 539 * ncomps + c);
        const auto *ih_540 = buffer.data(ih + 540 * ncomps + c);
        const auto *ih_541 = buffer.data(ih + 541 * ncomps + c);
        const auto *ih_542 = buffer.data(ih + 542 * ncomps + c);
        const auto *ih_543 = buffer.data(ih + 543 * ncomps + c);
        const auto *ih_544 = buffer.data(ih + 544 * ncomps + c);
        const auto *ih_545 = buffer.data(ih + 545 * ncomps + c);

        const auto *kh_440 = buffer.data(kh + 440 * ncomps + c);
        const auto *kh_441 = buffer.data(kh + 441 * ncomps + c);
        const auto *kh_442 = buffer.data(kh + 442 * ncomps + c);
        const auto *kh_443 = buffer.data(kh + 443 * ncomps + c);
        const auto *kh_444 = buffer.data(kh + 444 * ncomps + c);
        const auto *kh_445 = buffer.data(kh + 445 * ncomps + c);
        const auto *kh_446 = buffer.data(kh + 446 * ncomps + c);
        const auto *kh_447 = buffer.data(kh + 447 * ncomps + c);
        const auto *kh_448 = buffer.data(kh + 448 * ncomps + c);
        const auto *kh_449 = buffer.data(kh + 449 * ncomps + c);
        const auto *kh_450 = buffer.data(kh + 450 * ncomps + c);
        const auto *kh_451 = buffer.data(kh + 451 * ncomps + c);
        const auto *kh_452 = buffer.data(kh + 452 * ncomps + c);
        const auto *kh_453 = buffer.data(kh + 453 * ncomps + c);
        const auto *kh_454 = buffer.data(kh + 454 * ncomps + c);
        const auto *kh_455 = buffer.data(kh + 455 * ncomps + c);
        const auto *kh_456 = buffer.data(kh + 456 * ncomps + c);
        const auto *kh_457 = buffer.data(kh + 457 * ncomps + c);
        const auto *kh_458 = buffer.data(kh + 458 * ncomps + c);
        const auto *kh_459 = buffer.data(kh + 459 * ncomps + c);
        const auto *kh_460 = buffer.data(kh + 460 * ncomps + c);
        const auto *kh_461 = buffer.data(kh + 461 * ncomps + c);
        const auto *kh_462 = buffer.data(kh + 462 * ncomps + c);
        const auto *kh_463 = buffer.data(kh + 463 * ncomps + c);
        const auto *kh_464 = buffer.data(kh + 464 * ncomps + c);
        const auto *kh_465 = buffer.data(kh + 465 * ncomps + c);
        const auto *kh_466 = buffer.data(kh + 466 * ncomps + c);
        const auto *kh_467 = buffer.data(kh + 467 * ncomps + c);
        const auto *kh_468 = buffer.data(kh + 468 * ncomps + c);
        const auto *kh_469 = buffer.data(kh + 469 * ncomps + c);
        const auto *kh_470 = buffer.data(kh + 470 * ncomps + c);
        const auto *kh_471 = buffer.data(kh + 471 * ncomps + c);
        const auto *kh_472 = buffer.data(kh + 472 * ncomps + c);
        const auto *kh_473 = buffer.data(kh + 473 * ncomps + c);
        const auto *kh_474 = buffer.data(kh + 474 * ncomps + c);
        const auto *kh_475 = buffer.data(kh + 475 * ncomps + c);
        const auto *kh_476 = buffer.data(kh + 476 * ncomps + c);
        const auto *kh_477 = buffer.data(kh + 477 * ncomps + c);
        const auto *kh_478 = buffer.data(kh + 478 * ncomps + c);
        const auto *kh_479 = buffer.data(kh + 479 * ncomps + c);
        const auto *kh_480 = buffer.data(kh + 480 * ncomps + c);
        const auto *kh_481 = buffer.data(kh + 481 * ncomps + c);
        const auto *kh_482 = buffer.data(kh + 482 * ncomps + c);
        const auto *kh_483 = buffer.data(kh + 483 * ncomps + c);
        const auto *kh_484 = buffer.data(kh + 484 * ncomps + c);
        const auto *kh_485 = buffer.data(kh + 485 * ncomps + c);
        const auto *kh_486 = buffer.data(kh + 486 * ncomps + c);
        const auto *kh_487 = buffer.data(kh + 487 * ncomps + c);
        const auto *kh_488 = buffer.data(kh + 488 * ncomps + c);
        const auto *kh_489 = buffer.data(kh + 489 * ncomps + c);
        const auto *kh_490 = buffer.data(kh + 490 * ncomps + c);
        const auto *kh_491 = buffer.data(kh + 491 * ncomps + c);
        const auto *kh_492 = buffer.data(kh + 492 * ncomps + c);
        const auto *kh_493 = buffer.data(kh + 493 * ncomps + c);
        const auto *kh_494 = buffer.data(kh + 494 * ncomps + c);
        const auto *kh_495 = buffer.data(kh + 495 * ncomps + c);
        const auto *kh_496 = buffer.data(kh + 496 * ncomps + c);
        const auto *kh_497 = buffer.data(kh + 497 * ncomps + c);
        const auto *kh_498 = buffer.data(kh + 498 * ncomps + c);
        const auto *kh_499 = buffer.data(kh + 499 * ncomps + c);
        const auto *kh_500 = buffer.data(kh + 500 * ncomps + c);
        const auto *kh_501 = buffer.data(kh + 501 * ncomps + c);
        const auto *kh_502 = buffer.data(kh + 502 * ncomps + c);
        const auto *kh_503 = buffer.data(kh + 503 * ncomps + c);
        const auto *kh_504 = buffer.data(kh + 504 * ncomps + c);
        const auto *kh_505 = buffer.data(kh + 505 * ncomps + c);
        const auto *kh_506 = buffer.data(kh + 506 * ncomps + c);
        const auto *kh_507 = buffer.data(kh + 507 * ncomps + c);
        const auto *kh_508 = buffer.data(kh + 508 * ncomps + c);
        const auto *kh_509 = buffer.data(kh + 509 * ncomps + c);
        const auto *kh_510 = buffer.data(kh + 510 * ncomps + c);
        const auto *kh_511 = buffer.data(kh + 511 * ncomps + c);
        const auto *kh_512 = buffer.data(kh + 512 * ncomps + c);
        const auto *kh_513 = buffer.data(kh + 513 * ncomps + c);
        const auto *kh_514 = buffer.data(kh + 514 * ncomps + c);
        const auto *kh_515 = buffer.data(kh + 515 * ncomps + c);
        const auto *kh_516 = buffer.data(kh + 516 * ncomps + c);
        const auto *kh_517 = buffer.data(kh + 517 * ncomps + c);
        const auto *kh_518 = buffer.data(kh + 518 * ncomps + c);
        const auto *kh_519 = buffer.data(kh + 519 * ncomps + c);
        const auto *kh_520 = buffer.data(kh + 520 * ncomps + c);
        const auto *kh_521 = buffer.data(kh + 521 * ncomps + c);
        const auto *kh_522 = buffer.data(kh + 522 * ncomps + c);
        const auto *kh_523 = buffer.data(kh + 523 * ncomps + c);
        const auto *kh_524 = buffer.data(kh + 524 * ncomps + c);
        const auto *kh_525 = buffer.data(kh + 525 * ncomps + c);
        const auto *kh_526 = buffer.data(kh + 526 * ncomps + c);
        const auto *kh_527 = buffer.data(kh + 527 * ncomps + c);
        const auto *kh_528 = buffer.data(kh + 528 * ncomps + c);
        const auto *kh_529 = buffer.data(kh + 529 * ncomps + c);
        const auto *kh_530 = buffer.data(kh + 530 * ncomps + c);
        const auto *kh_531 = buffer.data(kh + 531 * ncomps + c);
        const auto *kh_532 = buffer.data(kh + 532 * ncomps + c);
        const auto *kh_533 = buffer.data(kh + 533 * ncomps + c);
        const auto *kh_534 = buffer.data(kh + 534 * ncomps + c);
        const auto *kh_535 = buffer.data(kh + 535 * ncomps + c);
        const auto *kh_536 = buffer.data(kh + 536 * ncomps + c);
        const auto *kh_537 = buffer.data(kh + 537 * ncomps + c);
        const auto *kh_538 = buffer.data(kh + 538 * ncomps + c);
        const auto *kh_539 = buffer.data(kh + 539 * ncomps + c);
        const auto *kh_540 = buffer.data(kh + 540 * ncomps + c);
        const auto *kh_541 = buffer.data(kh + 541 * ncomps + c);
        const auto *kh_542 = buffer.data(kh + 542 * ncomps + c);
        const auto *kh_543 = buffer.data(kh + 543 * ncomps + c);
        const auto *kh_544 = buffer.data(kh + 544 * ncomps + c);
        const auto *kh_545 = buffer.data(kh + 545 * ncomps + c);
        const auto *kh_561 = buffer.data(kh + 561 * ncomps + c);
        const auto *kh_562 = buffer.data(kh + 562 * ncomps + c);
        const auto *kh_563 = buffer.data(kh + 563 * ncomps + c);
        const auto *kh_564 = buffer.data(kh + 564 * ncomps + c);
        const auto *kh_565 = buffer.data(kh + 565 * ncomps + c);
        const auto *kh_566 = buffer.data(kh + 566 * ncomps + c);
        const auto *kh_587 = buffer.data(kh + 587 * ncomps + c);
        const auto *kh_603 = buffer.data(kh + 603 * ncomps + c);
        const auto *kh_604 = buffer.data(kh + 604 * ncomps + c);
        const auto *kh_605 = buffer.data(kh + 605 * ncomps + c);
        const auto *kh_606 = buffer.data(kh + 606 * ncomps + c);
        const auto *kh_607 = buffer.data(kh + 607 * ncomps + c);
        const auto *kh_608 = buffer.data(kh + 608 * ncomps + c);
        const auto *kh_624 = buffer.data(kh + 624 * ncomps + c);
        const auto *kh_625 = buffer.data(kh + 625 * ncomps + c);
        const auto *kh_626 = buffer.data(kh + 626 * ncomps + c);
        const auto *kh_627 = buffer.data(kh + 627 * ncomps + c);
        const auto *kh_628 = buffer.data(kh + 628 * ncomps + c);
        const auto *kh_629 = buffer.data(kh + 629 * ncomps + c);
        const auto *kh_645 = buffer.data(kh + 645 * ncomps + c);
        const auto *kh_646 = buffer.data(kh + 646 * ncomps + c);
        const auto *kh_647 = buffer.data(kh + 647 * ncomps + c);
        const auto *kh_648 = buffer.data(kh + 648 * ncomps + c);
        const auto *kh_649 = buffer.data(kh + 649 * ncomps + c);
        const auto *kh_650 = buffer.data(kh + 650 * ncomps + c);
        const auto *kh_666 = buffer.data(kh + 666 * ncomps + c);
        const auto *kh_667 = buffer.data(kh + 667 * ncomps + c);
        const auto *kh_668 = buffer.data(kh + 668 * ncomps + c);
        const auto *kh_669 = buffer.data(kh + 669 * ncomps + c);
        const auto *kh_670 = buffer.data(kh + 670 * ncomps + c);
        const auto *kh_671 = buffer.data(kh + 671 * ncomps + c);
        const auto *kh_687 = buffer.data(kh + 687 * ncomps + c);
        const auto *kh_688 = buffer.data(kh + 688 * ncomps + c);
        const auto *kh_689 = buffer.data(kh + 689 * ncomps + c);
        const auto *kh_692 = buffer.data(kh + 692 * ncomps + c);

#pragma omp simd aligned(t_580, t_581, t_582, t_583, ab_x, ab_y, ih_435, ih_436, ih_437, \
                         ih_440, kh_440, kh_561, kh_562, kh_563 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = ab_x[k] * ih_440[k]
                       + kh_440[k];

            t_581[k] = ab_y[k] * ih_435[k]
                       + kh_561[k];

            t_582[k] = ab_y[k] * ih_436[k]
                       + kh_562[k];

            t_583[k] = ab_y[k] * ih_437[k]
                       + kh_563[k];
        }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, ab_y, ab_z, ih_438, ih_439, ih_440, \
                         kh_564, kh_565, kh_566, kh_587 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_584[k] = ab_y[k] * ih_438[k]
                       + kh_564[k];

            t_585[k] = ab_y[k] * ih_439[k]
                       + kh_565[k];

            t_586[k] = ab_y[k] * ih_440[k]
                       + kh_566[k];

            t_587[k] = ab_z[k] * ih_440[k]
                       + kh_587[k];
        }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, t_592, ab_x, ih_441, ih_442, ih_443, \
                         ih_444, ih_445, kh_441, kh_442, kh_443, kh_444, \
                         kh_445 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_588[k] = ab_x[k] * ih_441[k]
                       + kh_441[k];

            t_589[k] = ab_x[k] * ih_442[k]
                       + kh_442[k];

            t_590[k] = ab_x[k] * ih_443[k]
                       + kh_443[k];

            t_591[k] = ab_x[k] * ih_444[k]
                       + kh_444[k];

            t_592[k] = ab_x[k] * ih_445[k]
                       + kh_445[k];
        }

#pragma omp simd aligned(t_593, t_594, t_595, t_596, t_597, ab_x, ih_446, ih_447, ih_448, \
                         ih_449, ih_450, kh_446, kh_447, kh_448, kh_449, \
                         kh_450 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_593[k] = ab_x[k] * ih_446[k]
                       + kh_446[k];

            t_594[k] = ab_x[k] * ih_447[k]
                       + kh_447[k];

            t_595[k] = ab_x[k] * ih_448[k]
                       + kh_448[k];

            t_596[k] = ab_x[k] * ih_449[k]
                       + kh_449[k];

            t_597[k] = ab_x[k] * ih_450[k]
                       + kh_450[k];
        }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, t_602, ab_x, ih_451, ih_452, ih_453, \
                         ih_454, ih_455, kh_451, kh_452, kh_453, kh_454, \
                         kh_455 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_598[k] = ab_x[k] * ih_451[k]
                       + kh_451[k];

            t_599[k] = ab_x[k] * ih_452[k]
                       + kh_452[k];

            t_600[k] = ab_x[k] * ih_453[k]
                       + kh_453[k];

            t_601[k] = ab_x[k] * ih_454[k]
                       + kh_454[k];

            t_602[k] = ab_x[k] * ih_455[k]
                       + kh_455[k];
        }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, t_607, ab_x, ih_456, ih_457, ih_458, \
                         ih_459, ih_460, kh_456, kh_457, kh_458, kh_459, \
                         kh_460 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_603[k] = ab_x[k] * ih_456[k]
                       + kh_456[k];

            t_604[k] = ab_x[k] * ih_457[k]
                       + kh_457[k];

            t_605[k] = ab_x[k] * ih_458[k]
                       + kh_458[k];

            t_606[k] = ab_x[k] * ih_459[k]
                       + kh_459[k];

            t_607[k] = ab_x[k] * ih_460[k]
                       + kh_460[k];
        }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, ab_x, ab_y, ih_456, ih_457, ih_458, \
                         ih_461, kh_461, kh_603, kh_604, kh_605 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_608[k] = ab_x[k] * ih_461[k]
                       + kh_461[k];

            t_609[k] = ab_y[k] * ih_456[k]
                       + kh_603[k];

            t_610[k] = ab_y[k] * ih_457[k]
                       + kh_604[k];

            t_611[k] = ab_y[k] * ih_458[k]
                       + kh_605[k];
        }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, ab_y, ab_z, ih_459, ih_460, ih_461, \
                         kh_606, kh_607, kh_608, kh_629 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_612[k] = ab_y[k] * ih_459[k]
                       + kh_606[k];

            t_613[k] = ab_y[k] * ih_460[k]
                       + kh_607[k];

            t_614[k] = ab_y[k] * ih_461[k]
                       + kh_608[k];

            t_615[k] = ab_z[k] * ih_461[k]
                       + kh_629[k];
        }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, t_620, ab_x, ih_462, ih_463, ih_464, \
                         ih_465, ih_466, kh_462, kh_463, kh_464, kh_465, \
                         kh_466 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_616[k] = ab_x[k] * ih_462[k]
                       + kh_462[k];

            t_617[k] = ab_x[k] * ih_463[k]
                       + kh_463[k];

            t_618[k] = ab_x[k] * ih_464[k]
                       + kh_464[k];

            t_619[k] = ab_x[k] * ih_465[k]
                       + kh_465[k];

            t_620[k] = ab_x[k] * ih_466[k]
                       + kh_466[k];
        }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, t_625, ab_x, ih_467, ih_468, ih_469, \
                         ih_470, ih_471, kh_467, kh_468, kh_469, kh_470, \
                         kh_471 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_621[k] = ab_x[k] * ih_467[k]
                       + kh_467[k];

            t_622[k] = ab_x[k] * ih_468[k]
                       + kh_468[k];

            t_623[k] = ab_x[k] * ih_469[k]
                       + kh_469[k];

            t_624[k] = ab_x[k] * ih_470[k]
                       + kh_470[k];

            t_625[k] = ab_x[k] * ih_471[k]
                       + kh_471[k];
        }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, t_630, ab_x, ih_472, ih_473, ih_474, \
                         ih_475, ih_476, kh_472, kh_473, kh_474, kh_475, \
                         kh_476 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_626[k] = ab_x[k] * ih_472[k]
                       + kh_472[k];

            t_627[k] = ab_x[k] * ih_473[k]
                       + kh_473[k];

            t_628[k] = ab_x[k] * ih_474[k]
                       + kh_474[k];

            t_629[k] = ab_x[k] * ih_475[k]
                       + kh_475[k];

            t_630[k] = ab_x[k] * ih_476[k]
                       + kh_476[k];
        }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, ab_x, ih_477, ih_478, ih_479, \
                         ih_480, ih_481, kh_477, kh_478, kh_479, kh_480, \
                         kh_481 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_631[k] = ab_x[k] * ih_477[k]
                       + kh_477[k];

            t_632[k] = ab_x[k] * ih_478[k]
                       + kh_478[k];

            t_633[k] = ab_x[k] * ih_479[k]
                       + kh_479[k];

            t_634[k] = ab_x[k] * ih_480[k]
                       + kh_480[k];

            t_635[k] = ab_x[k] * ih_481[k]
                       + kh_481[k];
        }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, ab_x, ab_y, ih_477, ih_478, ih_479, \
                         ih_482, kh_482, kh_624, kh_625, kh_626 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_636[k] = ab_x[k] * ih_482[k]
                       + kh_482[k];

            t_637[k] = ab_y[k] * ih_477[k]
                       + kh_624[k];

            t_638[k] = ab_y[k] * ih_478[k]
                       + kh_625[k];

            t_639[k] = ab_y[k] * ih_479[k]
                       + kh_626[k];
        }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, ab_y, ab_z, ih_480, ih_481, ih_482, \
                         kh_627, kh_628, kh_629, kh_650 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_640[k] = ab_y[k] * ih_480[k]
                       + kh_627[k];

            t_641[k] = ab_y[k] * ih_481[k]
                       + kh_628[k];

            t_642[k] = ab_y[k] * ih_482[k]
                       + kh_629[k];

            t_643[k] = ab_z[k] * ih_482[k]
                       + kh_650[k];
        }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, ab_x, ih_483, ih_484, ih_485, \
                         ih_486, ih_487, kh_483, kh_484, kh_485, kh_486, \
                         kh_487 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_644[k] = ab_x[k] * ih_483[k]
                       + kh_483[k];

            t_645[k] = ab_x[k] * ih_484[k]
                       + kh_484[k];

            t_646[k] = ab_x[k] * ih_485[k]
                       + kh_485[k];

            t_647[k] = ab_x[k] * ih_486[k]
                       + kh_486[k];

            t_648[k] = ab_x[k] * ih_487[k]
                       + kh_487[k];
        }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, ab_x, ih_488, ih_489, ih_490, \
                         ih_491, ih_492, kh_488, kh_489, kh_490, kh_491, \
                         kh_492 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_649[k] = ab_x[k] * ih_488[k]
                       + kh_488[k];

            t_650[k] = ab_x[k] * ih_489[k]
                       + kh_489[k];

            t_651[k] = ab_x[k] * ih_490[k]
                       + kh_490[k];

            t_652[k] = ab_x[k] * ih_491[k]
                       + kh_491[k];

            t_653[k] = ab_x[k] * ih_492[k]
                       + kh_492[k];
        }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, ab_x, ih_493, ih_494, ih_495, \
                         ih_496, ih_497, kh_493, kh_494, kh_495, kh_496, \
                         kh_497 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_654[k] = ab_x[k] * ih_493[k]
                       + kh_493[k];

            t_655[k] = ab_x[k] * ih_494[k]
                       + kh_494[k];

            t_656[k] = ab_x[k] * ih_495[k]
                       + kh_495[k];

            t_657[k] = ab_x[k] * ih_496[k]
                       + kh_496[k];

            t_658[k] = ab_x[k] * ih_497[k]
                       + kh_497[k];
        }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, ab_x, ih_498, ih_499, ih_500, \
                         ih_501, ih_502, kh_498, kh_499, kh_500, kh_501, \
                         kh_502 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_659[k] = ab_x[k] * ih_498[k]
                       + kh_498[k];

            t_660[k] = ab_x[k] * ih_499[k]
                       + kh_499[k];

            t_661[k] = ab_x[k] * ih_500[k]
                       + kh_500[k];

            t_662[k] = ab_x[k] * ih_501[k]
                       + kh_501[k];

            t_663[k] = ab_x[k] * ih_502[k]
                       + kh_502[k];
        }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, ab_x, ab_y, ih_498, ih_499, ih_500, \
                         ih_503, kh_503, kh_645, kh_646, kh_647 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_664[k] = ab_x[k] * ih_503[k]
                       + kh_503[k];

            t_665[k] = ab_y[k] * ih_498[k]
                       + kh_645[k];

            t_666[k] = ab_y[k] * ih_499[k]
                       + kh_646[k];

            t_667[k] = ab_y[k] * ih_500[k]
                       + kh_647[k];
        }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, ab_y, ab_z, ih_501, ih_502, ih_503, \
                         kh_648, kh_649, kh_650, kh_671 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_668[k] = ab_y[k] * ih_501[k]
                       + kh_648[k];

            t_669[k] = ab_y[k] * ih_502[k]
                       + kh_649[k];

            t_670[k] = ab_y[k] * ih_503[k]
                       + kh_650[k];

            t_671[k] = ab_z[k] * ih_503[k]
                       + kh_671[k];
        }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, ab_x, ih_504, ih_505, ih_506, \
                         ih_507, ih_508, kh_504, kh_505, kh_506, kh_507, \
                         kh_508 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_672[k] = ab_x[k] * ih_504[k]
                       + kh_504[k];

            t_673[k] = ab_x[k] * ih_505[k]
                       + kh_505[k];

            t_674[k] = ab_x[k] * ih_506[k]
                       + kh_506[k];

            t_675[k] = ab_x[k] * ih_507[k]
                       + kh_507[k];

            t_676[k] = ab_x[k] * ih_508[k]
                       + kh_508[k];
        }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, ab_x, ih_509, ih_510, ih_511, \
                         ih_512, ih_513, kh_509, kh_510, kh_511, kh_512, \
                         kh_513 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_677[k] = ab_x[k] * ih_509[k]
                       + kh_509[k];

            t_678[k] = ab_x[k] * ih_510[k]
                       + kh_510[k];

            t_679[k] = ab_x[k] * ih_511[k]
                       + kh_511[k];

            t_680[k] = ab_x[k] * ih_512[k]
                       + kh_512[k];

            t_681[k] = ab_x[k] * ih_513[k]
                       + kh_513[k];
        }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, ab_x, ih_514, ih_515, ih_516, \
                         ih_517, ih_518, kh_514, kh_515, kh_516, kh_517, \
                         kh_518 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_682[k] = ab_x[k] * ih_514[k]
                       + kh_514[k];

            t_683[k] = ab_x[k] * ih_515[k]
                       + kh_515[k];

            t_684[k] = ab_x[k] * ih_516[k]
                       + kh_516[k];

            t_685[k] = ab_x[k] * ih_517[k]
                       + kh_517[k];

            t_686[k] = ab_x[k] * ih_518[k]
                       + kh_518[k];
        }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, ab_x, ih_519, ih_520, ih_521, \
                         ih_522, ih_523, kh_519, kh_520, kh_521, kh_522, \
                         kh_523 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_687[k] = ab_x[k] * ih_519[k]
                       + kh_519[k];

            t_688[k] = ab_x[k] * ih_520[k]
                       + kh_520[k];

            t_689[k] = ab_x[k] * ih_521[k]
                       + kh_521[k];

            t_690[k] = ab_x[k] * ih_522[k]
                       + kh_522[k];

            t_691[k] = ab_x[k] * ih_523[k]
                       + kh_523[k];
        }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, ab_x, ab_y, ih_519, ih_520, ih_521, \
                         ih_524, kh_524, kh_666, kh_667, kh_668 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_692[k] = ab_x[k] * ih_524[k]
                       + kh_524[k];

            t_693[k] = ab_y[k] * ih_519[k]
                       + kh_666[k];

            t_694[k] = ab_y[k] * ih_520[k]
                       + kh_667[k];

            t_695[k] = ab_y[k] * ih_521[k]
                       + kh_668[k];
        }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, ab_y, ab_z, ih_522, ih_523, ih_524, \
                         kh_669, kh_670, kh_671, kh_692 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_696[k] = ab_y[k] * ih_522[k]
                       + kh_669[k];

            t_697[k] = ab_y[k] * ih_523[k]
                       + kh_670[k];

            t_698[k] = ab_y[k] * ih_524[k]
                       + kh_671[k];

            t_699[k] = ab_z[k] * ih_524[k]
                       + kh_692[k];
        }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ab_x, ih_525, ih_526, ih_527, \
                         ih_528, ih_529, kh_525, kh_526, kh_527, kh_528, \
                         kh_529 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_700[k] = ab_x[k] * ih_525[k]
                       + kh_525[k];

            t_701[k] = ab_x[k] * ih_526[k]
                       + kh_526[k];

            t_702[k] = ab_x[k] * ih_527[k]
                       + kh_527[k];

            t_703[k] = ab_x[k] * ih_528[k]
                       + kh_528[k];

            t_704[k] = ab_x[k] * ih_529[k]
                       + kh_529[k];
        }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ab_x, ih_530, ih_531, ih_532, \
                         ih_533, ih_534, kh_530, kh_531, kh_532, kh_533, \
                         kh_534 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_705[k] = ab_x[k] * ih_530[k]
                       + kh_530[k];

            t_706[k] = ab_x[k] * ih_531[k]
                       + kh_531[k];

            t_707[k] = ab_x[k] * ih_532[k]
                       + kh_532[k];

            t_708[k] = ab_x[k] * ih_533[k]
                       + kh_533[k];

            t_709[k] = ab_x[k] * ih_534[k]
                       + kh_534[k];
        }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ab_x, ih_535, ih_536, ih_537, \
                         ih_538, ih_539, kh_535, kh_536, kh_537, kh_538, \
                         kh_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_710[k] = ab_x[k] * ih_535[k]
                       + kh_535[k];

            t_711[k] = ab_x[k] * ih_536[k]
                       + kh_536[k];

            t_712[k] = ab_x[k] * ih_537[k]
                       + kh_537[k];

            t_713[k] = ab_x[k] * ih_538[k]
                       + kh_538[k];

            t_714[k] = ab_x[k] * ih_539[k]
                       + kh_539[k];
        }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ab_x, ih_540, ih_541, ih_542, \
                         ih_543, ih_544, kh_540, kh_541, kh_542, kh_543, \
                         kh_544 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_715[k] = ab_x[k] * ih_540[k]
                       + kh_540[k];

            t_716[k] = ab_x[k] * ih_541[k]
                       + kh_541[k];

            t_717[k] = ab_x[k] * ih_542[k]
                       + kh_542[k];

            t_718[k] = ab_x[k] * ih_543[k]
                       + kh_543[k];

            t_719[k] = ab_x[k] * ih_544[k]
                       + kh_544[k];
        }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, ab_x, ab_y, ih_540, ih_541, ih_542, \
                         ih_545, kh_545, kh_687, kh_688, kh_689 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_720[k] = ab_x[k] * ih_545[k]
                       + kh_545[k];

            t_721[k] = ab_y[k] * ih_540[k]
                       + kh_687[k];

            t_722[k] = ab_y[k] * ih_541[k]
                       + kh_688[k];

            t_723[k] = ab_y[k] * ih_542[k]
                       + kh_689[k];
        }
    }
}

static auto
compute_hrr_ii_out_of_first_piece5(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ih, const size_t kh,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_724 = buffer.data(target + 724 * ncomps + c);
        auto *t_725 = buffer.data(target + 725 * ncomps + c);
        auto *t_726 = buffer.data(target + 726 * ncomps + c);
        auto *t_727 = buffer.data(target + 727 * ncomps + c);
        auto *t_728 = buffer.data(target + 728 * ncomps + c);
        auto *t_729 = buffer.data(target + 729 * ncomps + c);
        auto *t_730 = buffer.data(target + 730 * ncomps + c);
        auto *t_731 = buffer.data(target + 731 * ncomps + c);
        auto *t_732 = buffer.data(target + 732 * ncomps + c);
        auto *t_733 = buffer.data(target + 733 * ncomps + c);
        auto *t_734 = buffer.data(target + 734 * ncomps + c);
        auto *t_735 = buffer.data(target + 735 * ncomps + c);
        auto *t_736 = buffer.data(target + 736 * ncomps + c);
        auto *t_737 = buffer.data(target + 737 * ncomps + c);
        auto *t_738 = buffer.data(target + 738 * ncomps + c);
        auto *t_739 = buffer.data(target + 739 * ncomps + c);
        auto *t_740 = buffer.data(target + 740 * ncomps + c);
        auto *t_741 = buffer.data(target + 741 * ncomps + c);
        auto *t_742 = buffer.data(target + 742 * ncomps + c);
        auto *t_743 = buffer.data(target + 743 * ncomps + c);
        auto *t_744 = buffer.data(target + 744 * ncomps + c);
        auto *t_745 = buffer.data(target + 745 * ncomps + c);
        auto *t_746 = buffer.data(target + 746 * ncomps + c);
        auto *t_747 = buffer.data(target + 747 * ncomps + c);
        auto *t_748 = buffer.data(target + 748 * ncomps + c);
        auto *t_749 = buffer.data(target + 749 * ncomps + c);
        auto *t_750 = buffer.data(target + 750 * ncomps + c);
        auto *t_751 = buffer.data(target + 751 * ncomps + c);
        auto *t_752 = buffer.data(target + 752 * ncomps + c);
        auto *t_753 = buffer.data(target + 753 * ncomps + c);
        auto *t_754 = buffer.data(target + 754 * ncomps + c);
        auto *t_755 = buffer.data(target + 755 * ncomps + c);
        auto *t_756 = buffer.data(target + 756 * ncomps + c);
        auto *t_757 = buffer.data(target + 757 * ncomps + c);
        auto *t_758 = buffer.data(target + 758 * ncomps + c);
        auto *t_759 = buffer.data(target + 759 * ncomps + c);
        auto *t_760 = buffer.data(target + 760 * ncomps + c);
        auto *t_761 = buffer.data(target + 761 * ncomps + c);
        auto *t_762 = buffer.data(target + 762 * ncomps + c);
        auto *t_763 = buffer.data(target + 763 * ncomps + c);
        auto *t_764 = buffer.data(target + 764 * ncomps + c);
        auto *t_765 = buffer.data(target + 765 * ncomps + c);
        auto *t_766 = buffer.data(target + 766 * ncomps + c);
        auto *t_767 = buffer.data(target + 767 * ncomps + c);
        auto *t_768 = buffer.data(target + 768 * ncomps + c);
        auto *t_769 = buffer.data(target + 769 * ncomps + c);
        auto *t_770 = buffer.data(target + 770 * ncomps + c);
        auto *t_771 = buffer.data(target + 771 * ncomps + c);
        auto *t_772 = buffer.data(target + 772 * ncomps + c);
        auto *t_773 = buffer.data(target + 773 * ncomps + c);
        auto *t_774 = buffer.data(target + 774 * ncomps + c);
        auto *t_775 = buffer.data(target + 775 * ncomps + c);
        auto *t_776 = buffer.data(target + 776 * ncomps + c);
        auto *t_777 = buffer.data(target + 777 * ncomps + c);
        auto *t_778 = buffer.data(target + 778 * ncomps + c);
        auto *t_779 = buffer.data(target + 779 * ncomps + c);
        auto *t_780 = buffer.data(target + 780 * ncomps + c);
        auto *t_781 = buffer.data(target + 781 * ncomps + c);
        auto *t_782 = buffer.data(target + 782 * ncomps + c);
        auto *t_783 = buffer.data(target + 783 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ih_543 = buffer.data(ih + 543 * ncomps + c);
        const auto *ih_544 = buffer.data(ih + 544 * ncomps + c);
        const auto *ih_545 = buffer.data(ih + 545 * ncomps + c);
        const auto *ih_546 = buffer.data(ih + 546 * ncomps + c);
        const auto *ih_547 = buffer.data(ih + 547 * ncomps + c);
        const auto *ih_548 = buffer.data(ih + 548 * ncomps + c);
        const auto *ih_549 = buffer.data(ih + 549 * ncomps + c);
        const auto *ih_550 = buffer.data(ih + 550 * ncomps + c);
        const auto *ih_551 = buffer.data(ih + 551 * ncomps + c);
        const auto *ih_552 = buffer.data(ih + 552 * ncomps + c);
        const auto *ih_553 = buffer.data(ih + 553 * ncomps + c);
        const auto *ih_554 = buffer.data(ih + 554 * ncomps + c);
        const auto *ih_555 = buffer.data(ih + 555 * ncomps + c);
        const auto *ih_556 = buffer.data(ih + 556 * ncomps + c);
        const auto *ih_557 = buffer.data(ih + 557 * ncomps + c);
        const auto *ih_558 = buffer.data(ih + 558 * ncomps + c);
        const auto *ih_559 = buffer.data(ih + 559 * ncomps + c);
        const auto *ih_560 = buffer.data(ih + 560 * ncomps + c);
        const auto *ih_561 = buffer.data(ih + 561 * ncomps + c);
        const auto *ih_562 = buffer.data(ih + 562 * ncomps + c);
        const auto *ih_563 = buffer.data(ih + 563 * ncomps + c);
        const auto *ih_564 = buffer.data(ih + 564 * ncomps + c);
        const auto *ih_565 = buffer.data(ih + 565 * ncomps + c);
        const auto *ih_566 = buffer.data(ih + 566 * ncomps + c);
        const auto *ih_567 = buffer.data(ih + 567 * ncomps + c);
        const auto *ih_568 = buffer.data(ih + 568 * ncomps + c);
        const auto *ih_569 = buffer.data(ih + 569 * ncomps + c);
        const auto *ih_570 = buffer.data(ih + 570 * ncomps + c);
        const auto *ih_571 = buffer.data(ih + 571 * ncomps + c);
        const auto *ih_572 = buffer.data(ih + 572 * ncomps + c);
        const auto *ih_573 = buffer.data(ih + 573 * ncomps + c);
        const auto *ih_574 = buffer.data(ih + 574 * ncomps + c);
        const auto *ih_575 = buffer.data(ih + 575 * ncomps + c);
        const auto *ih_576 = buffer.data(ih + 576 * ncomps + c);
        const auto *ih_577 = buffer.data(ih + 577 * ncomps + c);
        const auto *ih_578 = buffer.data(ih + 578 * ncomps + c);
        const auto *ih_579 = buffer.data(ih + 579 * ncomps + c);
        const auto *ih_580 = buffer.data(ih + 580 * ncomps + c);
        const auto *ih_581 = buffer.data(ih + 581 * ncomps + c);
        const auto *ih_582 = buffer.data(ih + 582 * ncomps + c);
        const auto *ih_583 = buffer.data(ih + 583 * ncomps + c);
        const auto *ih_584 = buffer.data(ih + 584 * ncomps + c);
        const auto *ih_585 = buffer.data(ih + 585 * ncomps + c);
        const auto *ih_586 = buffer.data(ih + 586 * ncomps + c);
        const auto *ih_587 = buffer.data(ih + 587 * ncomps + c);

        const auto *kh_546 = buffer.data(kh + 546 * ncomps + c);
        const auto *kh_547 = buffer.data(kh + 547 * ncomps + c);
        const auto *kh_548 = buffer.data(kh + 548 * ncomps + c);
        const auto *kh_549 = buffer.data(kh + 549 * ncomps + c);
        const auto *kh_550 = buffer.data(kh + 550 * ncomps + c);
        const auto *kh_551 = buffer.data(kh + 551 * ncomps + c);
        const auto *kh_552 = buffer.data(kh + 552 * ncomps + c);
        const auto *kh_553 = buffer.data(kh + 553 * ncomps + c);
        const auto *kh_554 = buffer.data(kh + 554 * ncomps + c);
        const auto *kh_555 = buffer.data(kh + 555 * ncomps + c);
        const auto *kh_556 = buffer.data(kh + 556 * ncomps + c);
        const auto *kh_557 = buffer.data(kh + 557 * ncomps + c);
        const auto *kh_558 = buffer.data(kh + 558 * ncomps + c);
        const auto *kh_559 = buffer.data(kh + 559 * ncomps + c);
        const auto *kh_560 = buffer.data(kh + 560 * ncomps + c);
        const auto *kh_561 = buffer.data(kh + 561 * ncomps + c);
        const auto *kh_562 = buffer.data(kh + 562 * ncomps + c);
        const auto *kh_563 = buffer.data(kh + 563 * ncomps + c);
        const auto *kh_564 = buffer.data(kh + 564 * ncomps + c);
        const auto *kh_565 = buffer.data(kh + 565 * ncomps + c);
        const auto *kh_566 = buffer.data(kh + 566 * ncomps + c);
        const auto *kh_567 = buffer.data(kh + 567 * ncomps + c);
        const auto *kh_568 = buffer.data(kh + 568 * ncomps + c);
        const auto *kh_569 = buffer.data(kh + 569 * ncomps + c);
        const auto *kh_570 = buffer.data(kh + 570 * ncomps + c);
        const auto *kh_571 = buffer.data(kh + 571 * ncomps + c);
        const auto *kh_572 = buffer.data(kh + 572 * ncomps + c);
        const auto *kh_573 = buffer.data(kh + 573 * ncomps + c);
        const auto *kh_574 = buffer.data(kh + 574 * ncomps + c);
        const auto *kh_575 = buffer.data(kh + 575 * ncomps + c);
        const auto *kh_576 = buffer.data(kh + 576 * ncomps + c);
        const auto *kh_577 = buffer.data(kh + 577 * ncomps + c);
        const auto *kh_578 = buffer.data(kh + 578 * ncomps + c);
        const auto *kh_579 = buffer.data(kh + 579 * ncomps + c);
        const auto *kh_580 = buffer.data(kh + 580 * ncomps + c);
        const auto *kh_581 = buffer.data(kh + 581 * ncomps + c);
        const auto *kh_582 = buffer.data(kh + 582 * ncomps + c);
        const auto *kh_583 = buffer.data(kh + 583 * ncomps + c);
        const auto *kh_584 = buffer.data(kh + 584 * ncomps + c);
        const auto *kh_585 = buffer.data(kh + 585 * ncomps + c);
        const auto *kh_586 = buffer.data(kh + 586 * ncomps + c);
        const auto *kh_587 = buffer.data(kh + 587 * ncomps + c);
        const auto *kh_690 = buffer.data(kh + 690 * ncomps + c);
        const auto *kh_691 = buffer.data(kh + 691 * ncomps + c);
        const auto *kh_692 = buffer.data(kh + 692 * ncomps + c);
        const auto *kh_708 = buffer.data(kh + 708 * ncomps + c);
        const auto *kh_709 = buffer.data(kh + 709 * ncomps + c);
        const auto *kh_710 = buffer.data(kh + 710 * ncomps + c);
        const auto *kh_711 = buffer.data(kh + 711 * ncomps + c);
        const auto *kh_712 = buffer.data(kh + 712 * ncomps + c);
        const auto *kh_713 = buffer.data(kh + 713 * ncomps + c);
        const auto *kh_729 = buffer.data(kh + 729 * ncomps + c);
        const auto *kh_730 = buffer.data(kh + 730 * ncomps + c);
        const auto *kh_731 = buffer.data(kh + 731 * ncomps + c);
        const auto *kh_732 = buffer.data(kh + 732 * ncomps + c);
        const auto *kh_733 = buffer.data(kh + 733 * ncomps + c);
        const auto *kh_734 = buffer.data(kh + 734 * ncomps + c);
        const auto *kh_755 = buffer.data(kh + 755 * ncomps + c);

#pragma omp simd aligned(t_724, t_725, t_726, t_727, ab_y, ab_z, ih_543, ih_544, ih_545, \
                         kh_690, kh_691, kh_692, kh_713 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_724[k] = ab_y[k] * ih_543[k]
                       + kh_690[k];

            t_725[k] = ab_y[k] * ih_544[k]
                       + kh_691[k];

            t_726[k] = ab_y[k] * ih_545[k]
                       + kh_692[k];

            t_727[k] = ab_z[k] * ih_545[k]
                       + kh_713[k];
        }

#pragma omp simd aligned(t_728, t_729, t_730, t_731, t_732, ab_x, ih_546, ih_547, ih_548, \
                         ih_549, ih_550, kh_546, kh_547, kh_548, kh_549, \
                         kh_550 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_728[k] = ab_x[k] * ih_546[k]
                       + kh_546[k];

            t_729[k] = ab_x[k] * ih_547[k]
                       + kh_547[k];

            t_730[k] = ab_x[k] * ih_548[k]
                       + kh_548[k];

            t_731[k] = ab_x[k] * ih_549[k]
                       + kh_549[k];

            t_732[k] = ab_x[k] * ih_550[k]
                       + kh_550[k];
        }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, t_737, ab_x, ih_551, ih_552, ih_553, \
                         ih_554, ih_555, kh_551, kh_552, kh_553, kh_554, \
                         kh_555 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_733[k] = ab_x[k] * ih_551[k]
                       + kh_551[k];

            t_734[k] = ab_x[k] * ih_552[k]
                       + kh_552[k];

            t_735[k] = ab_x[k] * ih_553[k]
                       + kh_553[k];

            t_736[k] = ab_x[k] * ih_554[k]
                       + kh_554[k];

            t_737[k] = ab_x[k] * ih_555[k]
                       + kh_555[k];
        }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, t_742, ab_x, ih_556, ih_557, ih_558, \
                         ih_559, ih_560, kh_556, kh_557, kh_558, kh_559, \
                         kh_560 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_738[k] = ab_x[k] * ih_556[k]
                       + kh_556[k];

            t_739[k] = ab_x[k] * ih_557[k]
                       + kh_557[k];

            t_740[k] = ab_x[k] * ih_558[k]
                       + kh_558[k];

            t_741[k] = ab_x[k] * ih_559[k]
                       + kh_559[k];

            t_742[k] = ab_x[k] * ih_560[k]
                       + kh_560[k];
        }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, ab_x, ih_561, ih_562, ih_563, \
                         ih_564, ih_565, kh_561, kh_562, kh_563, kh_564, \
                         kh_565 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_743[k] = ab_x[k] * ih_561[k]
                       + kh_561[k];

            t_744[k] = ab_x[k] * ih_562[k]
                       + kh_562[k];

            t_745[k] = ab_x[k] * ih_563[k]
                       + kh_563[k];

            t_746[k] = ab_x[k] * ih_564[k]
                       + kh_564[k];

            t_747[k] = ab_x[k] * ih_565[k]
                       + kh_565[k];
        }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, ab_x, ab_y, ih_561, ih_562, ih_563, \
                         ih_566, kh_566, kh_708, kh_709, kh_710 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_748[k] = ab_x[k] * ih_566[k]
                       + kh_566[k];

            t_749[k] = ab_y[k] * ih_561[k]
                       + kh_708[k];

            t_750[k] = ab_y[k] * ih_562[k]
                       + kh_709[k];

            t_751[k] = ab_y[k] * ih_563[k]
                       + kh_710[k];
        }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, ab_y, ab_z, ih_564, ih_565, ih_566, \
                         kh_711, kh_712, kh_713, kh_734 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_752[k] = ab_y[k] * ih_564[k]
                       + kh_711[k];

            t_753[k] = ab_y[k] * ih_565[k]
                       + kh_712[k];

            t_754[k] = ab_y[k] * ih_566[k]
                       + kh_713[k];

            t_755[k] = ab_z[k] * ih_566[k]
                       + kh_734[k];
        }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, ab_x, ih_567, ih_568, ih_569, \
                         ih_570, ih_571, kh_567, kh_568, kh_569, kh_570, \
                         kh_571 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_756[k] = ab_x[k] * ih_567[k]
                       + kh_567[k];

            t_757[k] = ab_x[k] * ih_568[k]
                       + kh_568[k];

            t_758[k] = ab_x[k] * ih_569[k]
                       + kh_569[k];

            t_759[k] = ab_x[k] * ih_570[k]
                       + kh_570[k];

            t_760[k] = ab_x[k] * ih_571[k]
                       + kh_571[k];
        }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, t_765, ab_x, ih_572, ih_573, ih_574, \
                         ih_575, ih_576, kh_572, kh_573, kh_574, kh_575, \
                         kh_576 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_761[k] = ab_x[k] * ih_572[k]
                       + kh_572[k];

            t_762[k] = ab_x[k] * ih_573[k]
                       + kh_573[k];

            t_763[k] = ab_x[k] * ih_574[k]
                       + kh_574[k];

            t_764[k] = ab_x[k] * ih_575[k]
                       + kh_575[k];

            t_765[k] = ab_x[k] * ih_576[k]
                       + kh_576[k];
        }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, t_770, ab_x, ih_577, ih_578, ih_579, \
                         ih_580, ih_581, kh_577, kh_578, kh_579, kh_580, \
                         kh_581 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_766[k] = ab_x[k] * ih_577[k]
                       + kh_577[k];

            t_767[k] = ab_x[k] * ih_578[k]
                       + kh_578[k];

            t_768[k] = ab_x[k] * ih_579[k]
                       + kh_579[k];

            t_769[k] = ab_x[k] * ih_580[k]
                       + kh_580[k];

            t_770[k] = ab_x[k] * ih_581[k]
                       + kh_581[k];
        }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, t_775, ab_x, ih_582, ih_583, ih_584, \
                         ih_585, ih_586, kh_582, kh_583, kh_584, kh_585, \
                         kh_586 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_771[k] = ab_x[k] * ih_582[k]
                       + kh_582[k];

            t_772[k] = ab_x[k] * ih_583[k]
                       + kh_583[k];

            t_773[k] = ab_x[k] * ih_584[k]
                       + kh_584[k];

            t_774[k] = ab_x[k] * ih_585[k]
                       + kh_585[k];

            t_775[k] = ab_x[k] * ih_586[k]
                       + kh_586[k];
        }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, ab_x, ab_y, ih_582, ih_583, ih_584, \
                         ih_587, kh_587, kh_729, kh_730, kh_731 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_776[k] = ab_x[k] * ih_587[k]
                       + kh_587[k];

            t_777[k] = ab_y[k] * ih_582[k]
                       + kh_729[k];

            t_778[k] = ab_y[k] * ih_583[k]
                       + kh_730[k];

            t_779[k] = ab_y[k] * ih_584[k]
                       + kh_731[k];
        }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, ab_y, ab_z, ih_585, ih_586, ih_587, \
                         kh_732, kh_733, kh_734, kh_755 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_780[k] = ab_y[k] * ih_585[k]
                       + kh_732[k];

            t_781[k] = ab_y[k] * ih_586[k]
                       + kh_733[k];

            t_782[k] = ab_y[k] * ih_587[k]
                       + kh_734[k];

            t_783[k] = ab_z[k] * ih_587[k]
                       + kh_755[k];
        }
    }
}

auto
compute_hrr_ii_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t ih, const size_t kh,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_ii_out_of_first_piece0(buffer, coordinates, target, ih, kh, ncomps, nmax);

    compute_hrr_ii_out_of_first_piece1(buffer, coordinates, target, ih, kh, ncomps, nmax);

    compute_hrr_ii_out_of_first_piece2(buffer, coordinates, target, ih, kh, ncomps, nmax);

    compute_hrr_ii_out_of_first_piece3(buffer, coordinates, target, ih, kh, ncomps, nmax);

    compute_hrr_ii_out_of_first_piece4(buffer, coordinates, target, ih, kh, ncomps, nmax);

    compute_hrr_ii_out_of_first_piece5(buffer, coordinates, target, ih, kh, ncomps, nmax);
}

static auto
compute_hrr_ii_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hi, const size_t hk, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);

        const auto *hi_0 = buffer.data(hi + 0 * ncomps + c);
        const auto *hi_1 = buffer.data(hi + 1 * ncomps + c);
        const auto *hi_2 = buffer.data(hi + 2 * ncomps + c);
        const auto *hi_3 = buffer.data(hi + 3 * ncomps + c);
        const auto *hi_4 = buffer.data(hi + 4 * ncomps + c);
        const auto *hi_5 = buffer.data(hi + 5 * ncomps + c);
        const auto *hi_6 = buffer.data(hi + 6 * ncomps + c);
        const auto *hi_7 = buffer.data(hi + 7 * ncomps + c);
        const auto *hi_8 = buffer.data(hi + 8 * ncomps + c);
        const auto *hi_9 = buffer.data(hi + 9 * ncomps + c);
        const auto *hi_10 = buffer.data(hi + 10 * ncomps + c);
        const auto *hi_11 = buffer.data(hi + 11 * ncomps + c);
        const auto *hi_12 = buffer.data(hi + 12 * ncomps + c);
        const auto *hi_13 = buffer.data(hi + 13 * ncomps + c);
        const auto *hi_14 = buffer.data(hi + 14 * ncomps + c);
        const auto *hi_15 = buffer.data(hi + 15 * ncomps + c);
        const auto *hi_16 = buffer.data(hi + 16 * ncomps + c);
        const auto *hi_17 = buffer.data(hi + 17 * ncomps + c);
        const auto *hi_18 = buffer.data(hi + 18 * ncomps + c);
        const auto *hi_19 = buffer.data(hi + 19 * ncomps + c);
        const auto *hi_20 = buffer.data(hi + 20 * ncomps + c);
        const auto *hi_21 = buffer.data(hi + 21 * ncomps + c);
        const auto *hi_22 = buffer.data(hi + 22 * ncomps + c);
        const auto *hi_23 = buffer.data(hi + 23 * ncomps + c);
        const auto *hi_24 = buffer.data(hi + 24 * ncomps + c);
        const auto *hi_25 = buffer.data(hi + 25 * ncomps + c);
        const auto *hi_26 = buffer.data(hi + 26 * ncomps + c);
        const auto *hi_27 = buffer.data(hi + 27 * ncomps + c);
        const auto *hi_28 = buffer.data(hi + 28 * ncomps + c);
        const auto *hi_29 = buffer.data(hi + 29 * ncomps + c);
        const auto *hi_30 = buffer.data(hi + 30 * ncomps + c);
        const auto *hi_31 = buffer.data(hi + 31 * ncomps + c);
        const auto *hi_32 = buffer.data(hi + 32 * ncomps + c);
        const auto *hi_33 = buffer.data(hi + 33 * ncomps + c);
        const auto *hi_34 = buffer.data(hi + 34 * ncomps + c);
        const auto *hi_35 = buffer.data(hi + 35 * ncomps + c);
        const auto *hi_36 = buffer.data(hi + 36 * ncomps + c);
        const auto *hi_37 = buffer.data(hi + 37 * ncomps + c);
        const auto *hi_38 = buffer.data(hi + 38 * ncomps + c);
        const auto *hi_39 = buffer.data(hi + 39 * ncomps + c);
        const auto *hi_40 = buffer.data(hi + 40 * ncomps + c);
        const auto *hi_41 = buffer.data(hi + 41 * ncomps + c);
        const auto *hi_42 = buffer.data(hi + 42 * ncomps + c);
        const auto *hi_43 = buffer.data(hi + 43 * ncomps + c);
        const auto *hi_44 = buffer.data(hi + 44 * ncomps + c);
        const auto *hi_45 = buffer.data(hi + 45 * ncomps + c);
        const auto *hi_46 = buffer.data(hi + 46 * ncomps + c);
        const auto *hi_47 = buffer.data(hi + 47 * ncomps + c);
        const auto *hi_48 = buffer.data(hi + 48 * ncomps + c);
        const auto *hi_49 = buffer.data(hi + 49 * ncomps + c);
        const auto *hi_50 = buffer.data(hi + 50 * ncomps + c);
        const auto *hi_51 = buffer.data(hi + 51 * ncomps + c);
        const auto *hi_52 = buffer.data(hi + 52 * ncomps + c);
        const auto *hi_53 = buffer.data(hi + 53 * ncomps + c);
        const auto *hi_54 = buffer.data(hi + 54 * ncomps + c);
        const auto *hi_55 = buffer.data(hi + 55 * ncomps + c);
        const auto *hi_56 = buffer.data(hi + 56 * ncomps + c);
        const auto *hi_57 = buffer.data(hi + 57 * ncomps + c);
        const auto *hi_58 = buffer.data(hi + 58 * ncomps + c);
        const auto *hi_59 = buffer.data(hi + 59 * ncomps + c);
        const auto *hi_60 = buffer.data(hi + 60 * ncomps + c);
        const auto *hi_61 = buffer.data(hi + 61 * ncomps + c);
        const auto *hi_62 = buffer.data(hi + 62 * ncomps + c);
        const auto *hi_63 = buffer.data(hi + 63 * ncomps + c);
        const auto *hi_64 = buffer.data(hi + 64 * ncomps + c);
        const auto *hi_65 = buffer.data(hi + 65 * ncomps + c);
        const auto *hi_66 = buffer.data(hi + 66 * ncomps + c);
        const auto *hi_67 = buffer.data(hi + 67 * ncomps + c);
        const auto *hi_68 = buffer.data(hi + 68 * ncomps + c);
        const auto *hi_69 = buffer.data(hi + 69 * ncomps + c);
        const auto *hi_70 = buffer.data(hi + 70 * ncomps + c);
        const auto *hi_71 = buffer.data(hi + 71 * ncomps + c);
        const auto *hi_72 = buffer.data(hi + 72 * ncomps + c);
        const auto *hi_73 = buffer.data(hi + 73 * ncomps + c);
        const auto *hi_74 = buffer.data(hi + 74 * ncomps + c);
        const auto *hi_75 = buffer.data(hi + 75 * ncomps + c);
        const auto *hi_76 = buffer.data(hi + 76 * ncomps + c);
        const auto *hi_77 = buffer.data(hi + 77 * ncomps + c);
        const auto *hi_78 = buffer.data(hi + 78 * ncomps + c);
        const auto *hi_79 = buffer.data(hi + 79 * ncomps + c);
        const auto *hi_80 = buffer.data(hi + 80 * ncomps + c);
        const auto *hi_81 = buffer.data(hi + 81 * ncomps + c);
        const auto *hi_82 = buffer.data(hi + 82 * ncomps + c);
        const auto *hi_83 = buffer.data(hi + 83 * ncomps + c);
        const auto *hi_84 = buffer.data(hi + 84 * ncomps + c);
        const auto *hi_85 = buffer.data(hi + 85 * ncomps + c);
        const auto *hi_86 = buffer.data(hi + 86 * ncomps + c);
        const auto *hi_87 = buffer.data(hi + 87 * ncomps + c);
        const auto *hi_88 = buffer.data(hi + 88 * ncomps + c);
        const auto *hi_89 = buffer.data(hi + 89 * ncomps + c);
        const auto *hi_90 = buffer.data(hi + 90 * ncomps + c);
        const auto *hi_91 = buffer.data(hi + 91 * ncomps + c);
        const auto *hi_92 = buffer.data(hi + 92 * ncomps + c);
        const auto *hi_93 = buffer.data(hi + 93 * ncomps + c);
        const auto *hi_94 = buffer.data(hi + 94 * ncomps + c);
        const auto *hi_95 = buffer.data(hi + 95 * ncomps + c);
        const auto *hi_96 = buffer.data(hi + 96 * ncomps + c);
        const auto *hi_97 = buffer.data(hi + 97 * ncomps + c);
        const auto *hi_98 = buffer.data(hi + 98 * ncomps + c);
        const auto *hi_99 = buffer.data(hi + 99 * ncomps + c);
        const auto *hi_100 = buffer.data(hi + 100 * ncomps + c);
        const auto *hi_101 = buffer.data(hi + 101 * ncomps + c);
        const auto *hi_102 = buffer.data(hi + 102 * ncomps + c);
        const auto *hi_103 = buffer.data(hi + 103 * ncomps + c);
        const auto *hi_104 = buffer.data(hi + 104 * ncomps + c);
        const auto *hi_105 = buffer.data(hi + 105 * ncomps + c);
        const auto *hi_106 = buffer.data(hi + 106 * ncomps + c);
        const auto *hi_107 = buffer.data(hi + 107 * ncomps + c);
        const auto *hi_108 = buffer.data(hi + 108 * ncomps + c);
        const auto *hi_109 = buffer.data(hi + 109 * ncomps + c);
        const auto *hi_110 = buffer.data(hi + 110 * ncomps + c);
        const auto *hi_111 = buffer.data(hi + 111 * ncomps + c);
        const auto *hi_112 = buffer.data(hi + 112 * ncomps + c);
        const auto *hi_113 = buffer.data(hi + 113 * ncomps + c);
        const auto *hi_114 = buffer.data(hi + 114 * ncomps + c);
        const auto *hi_115 = buffer.data(hi + 115 * ncomps + c);
        const auto *hi_116 = buffer.data(hi + 116 * ncomps + c);
        const auto *hi_117 = buffer.data(hi + 117 * ncomps + c);
        const auto *hi_118 = buffer.data(hi + 118 * ncomps + c);
        const auto *hi_119 = buffer.data(hi + 119 * ncomps + c);
        const auto *hi_120 = buffer.data(hi + 120 * ncomps + c);
        const auto *hi_121 = buffer.data(hi + 121 * ncomps + c);
        const auto *hi_122 = buffer.data(hi + 122 * ncomps + c);
        const auto *hi_123 = buffer.data(hi + 123 * ncomps + c);
        const auto *hi_124 = buffer.data(hi + 124 * ncomps + c);
        const auto *hi_125 = buffer.data(hi + 125 * ncomps + c);
        const auto *hi_126 = buffer.data(hi + 126 * ncomps + c);
        const auto *hi_127 = buffer.data(hi + 127 * ncomps + c);
        const auto *hi_128 = buffer.data(hi + 128 * ncomps + c);
        const auto *hi_129 = buffer.data(hi + 129 * ncomps + c);
        const auto *hi_130 = buffer.data(hi + 130 * ncomps + c);
        const auto *hi_131 = buffer.data(hi + 131 * ncomps + c);
        const auto *hi_132 = buffer.data(hi + 132 * ncomps + c);
        const auto *hi_133 = buffer.data(hi + 133 * ncomps + c);
        const auto *hi_134 = buffer.data(hi + 134 * ncomps + c);
        const auto *hi_135 = buffer.data(hi + 135 * ncomps + c);
        const auto *hi_136 = buffer.data(hi + 136 * ncomps + c);
        const auto *hi_137 = buffer.data(hi + 137 * ncomps + c);
        const auto *hi_138 = buffer.data(hi + 138 * ncomps + c);
        const auto *hi_139 = buffer.data(hi + 139 * ncomps + c);
        const auto *hi_140 = buffer.data(hi + 140 * ncomps + c);
        const auto *hi_141 = buffer.data(hi + 141 * ncomps + c);
        const auto *hi_142 = buffer.data(hi + 142 * ncomps + c);
        const auto *hi_143 = buffer.data(hi + 143 * ncomps + c);
        const auto *hi_144 = buffer.data(hi + 144 * ncomps + c);

        const auto *hk_0 = buffer.data(hk + 0 * ncomps + c);
        const auto *hk_1 = buffer.data(hk + 1 * ncomps + c);
        const auto *hk_2 = buffer.data(hk + 2 * ncomps + c);
        const auto *hk_3 = buffer.data(hk + 3 * ncomps + c);
        const auto *hk_4 = buffer.data(hk + 4 * ncomps + c);
        const auto *hk_5 = buffer.data(hk + 5 * ncomps + c);
        const auto *hk_6 = buffer.data(hk + 6 * ncomps + c);
        const auto *hk_7 = buffer.data(hk + 7 * ncomps + c);
        const auto *hk_8 = buffer.data(hk + 8 * ncomps + c);
        const auto *hk_9 = buffer.data(hk + 9 * ncomps + c);
        const auto *hk_10 = buffer.data(hk + 10 * ncomps + c);
        const auto *hk_11 = buffer.data(hk + 11 * ncomps + c);
        const auto *hk_12 = buffer.data(hk + 12 * ncomps + c);
        const auto *hk_13 = buffer.data(hk + 13 * ncomps + c);
        const auto *hk_14 = buffer.data(hk + 14 * ncomps + c);
        const auto *hk_15 = buffer.data(hk + 15 * ncomps + c);
        const auto *hk_16 = buffer.data(hk + 16 * ncomps + c);
        const auto *hk_17 = buffer.data(hk + 17 * ncomps + c);
        const auto *hk_18 = buffer.data(hk + 18 * ncomps + c);
        const auto *hk_19 = buffer.data(hk + 19 * ncomps + c);
        const auto *hk_20 = buffer.data(hk + 20 * ncomps + c);
        const auto *hk_21 = buffer.data(hk + 21 * ncomps + c);
        const auto *hk_22 = buffer.data(hk + 22 * ncomps + c);
        const auto *hk_23 = buffer.data(hk + 23 * ncomps + c);
        const auto *hk_24 = buffer.data(hk + 24 * ncomps + c);
        const auto *hk_25 = buffer.data(hk + 25 * ncomps + c);
        const auto *hk_26 = buffer.data(hk + 26 * ncomps + c);
        const auto *hk_27 = buffer.data(hk + 27 * ncomps + c);
        const auto *hk_36 = buffer.data(hk + 36 * ncomps + c);
        const auto *hk_37 = buffer.data(hk + 37 * ncomps + c);
        const auto *hk_38 = buffer.data(hk + 38 * ncomps + c);
        const auto *hk_39 = buffer.data(hk + 39 * ncomps + c);
        const auto *hk_40 = buffer.data(hk + 40 * ncomps + c);
        const auto *hk_41 = buffer.data(hk + 41 * ncomps + c);
        const auto *hk_42 = buffer.data(hk + 42 * ncomps + c);
        const auto *hk_43 = buffer.data(hk + 43 * ncomps + c);
        const auto *hk_44 = buffer.data(hk + 44 * ncomps + c);
        const auto *hk_45 = buffer.data(hk + 45 * ncomps + c);
        const auto *hk_46 = buffer.data(hk + 46 * ncomps + c);
        const auto *hk_47 = buffer.data(hk + 47 * ncomps + c);
        const auto *hk_48 = buffer.data(hk + 48 * ncomps + c);
        const auto *hk_49 = buffer.data(hk + 49 * ncomps + c);
        const auto *hk_50 = buffer.data(hk + 50 * ncomps + c);
        const auto *hk_51 = buffer.data(hk + 51 * ncomps + c);
        const auto *hk_52 = buffer.data(hk + 52 * ncomps + c);
        const auto *hk_53 = buffer.data(hk + 53 * ncomps + c);
        const auto *hk_54 = buffer.data(hk + 54 * ncomps + c);
        const auto *hk_55 = buffer.data(hk + 55 * ncomps + c);
        const auto *hk_56 = buffer.data(hk + 56 * ncomps + c);
        const auto *hk_57 = buffer.data(hk + 57 * ncomps + c);
        const auto *hk_58 = buffer.data(hk + 58 * ncomps + c);
        const auto *hk_59 = buffer.data(hk + 59 * ncomps + c);
        const auto *hk_60 = buffer.data(hk + 60 * ncomps + c);
        const auto *hk_61 = buffer.data(hk + 61 * ncomps + c);
        const auto *hk_62 = buffer.data(hk + 62 * ncomps + c);
        const auto *hk_63 = buffer.data(hk + 63 * ncomps + c);
        const auto *hk_72 = buffer.data(hk + 72 * ncomps + c);
        const auto *hk_73 = buffer.data(hk + 73 * ncomps + c);
        const auto *hk_74 = buffer.data(hk + 74 * ncomps + c);
        const auto *hk_75 = buffer.data(hk + 75 * ncomps + c);
        const auto *hk_76 = buffer.data(hk + 76 * ncomps + c);
        const auto *hk_77 = buffer.data(hk + 77 * ncomps + c);
        const auto *hk_78 = buffer.data(hk + 78 * ncomps + c);
        const auto *hk_79 = buffer.data(hk + 79 * ncomps + c);
        const auto *hk_80 = buffer.data(hk + 80 * ncomps + c);
        const auto *hk_81 = buffer.data(hk + 81 * ncomps + c);
        const auto *hk_82 = buffer.data(hk + 82 * ncomps + c);
        const auto *hk_83 = buffer.data(hk + 83 * ncomps + c);
        const auto *hk_84 = buffer.data(hk + 84 * ncomps + c);
        const auto *hk_85 = buffer.data(hk + 85 * ncomps + c);
        const auto *hk_86 = buffer.data(hk + 86 * ncomps + c);
        const auto *hk_87 = buffer.data(hk + 87 * ncomps + c);
        const auto *hk_88 = buffer.data(hk + 88 * ncomps + c);
        const auto *hk_89 = buffer.data(hk + 89 * ncomps + c);
        const auto *hk_90 = buffer.data(hk + 90 * ncomps + c);
        const auto *hk_91 = buffer.data(hk + 91 * ncomps + c);
        const auto *hk_92 = buffer.data(hk + 92 * ncomps + c);
        const auto *hk_93 = buffer.data(hk + 93 * ncomps + c);
        const auto *hk_94 = buffer.data(hk + 94 * ncomps + c);
        const auto *hk_95 = buffer.data(hk + 95 * ncomps + c);
        const auto *hk_96 = buffer.data(hk + 96 * ncomps + c);
        const auto *hk_97 = buffer.data(hk + 97 * ncomps + c);
        const auto *hk_98 = buffer.data(hk + 98 * ncomps + c);
        const auto *hk_99 = buffer.data(hk + 99 * ncomps + c);
        const auto *hk_108 = buffer.data(hk + 108 * ncomps + c);
        const auto *hk_109 = buffer.data(hk + 109 * ncomps + c);
        const auto *hk_110 = buffer.data(hk + 110 * ncomps + c);
        const auto *hk_111 = buffer.data(hk + 111 * ncomps + c);
        const auto *hk_112 = buffer.data(hk + 112 * ncomps + c);
        const auto *hk_113 = buffer.data(hk + 113 * ncomps + c);
        const auto *hk_114 = buffer.data(hk + 114 * ncomps + c);
        const auto *hk_115 = buffer.data(hk + 115 * ncomps + c);
        const auto *hk_116 = buffer.data(hk + 116 * ncomps + c);
        const auto *hk_117 = buffer.data(hk + 117 * ncomps + c);
        const auto *hk_118 = buffer.data(hk + 118 * ncomps + c);
        const auto *hk_119 = buffer.data(hk + 119 * ncomps + c);
        const auto *hk_120 = buffer.data(hk + 120 * ncomps + c);
        const auto *hk_121 = buffer.data(hk + 121 * ncomps + c);
        const auto *hk_122 = buffer.data(hk + 122 * ncomps + c);
        const auto *hk_123 = buffer.data(hk + 123 * ncomps + c);
        const auto *hk_124 = buffer.data(hk + 124 * ncomps + c);
        const auto *hk_125 = buffer.data(hk + 125 * ncomps + c);
        const auto *hk_126 = buffer.data(hk + 126 * ncomps + c);
        const auto *hk_127 = buffer.data(hk + 127 * ncomps + c);
        const auto *hk_128 = buffer.data(hk + 128 * ncomps + c);
        const auto *hk_129 = buffer.data(hk + 129 * ncomps + c);
        const auto *hk_130 = buffer.data(hk + 130 * ncomps + c);
        const auto *hk_131 = buffer.data(hk + 131 * ncomps + c);
        const auto *hk_132 = buffer.data(hk + 132 * ncomps + c);
        const auto *hk_133 = buffer.data(hk + 133 * ncomps + c);
        const auto *hk_134 = buffer.data(hk + 134 * ncomps + c);
        const auto *hk_135 = buffer.data(hk + 135 * ncomps + c);
        const auto *hk_144 = buffer.data(hk + 144 * ncomps + c);
        const auto *hk_145 = buffer.data(hk + 145 * ncomps + c);
        const auto *hk_146 = buffer.data(hk + 146 * ncomps + c);
        const auto *hk_147 = buffer.data(hk + 147 * ncomps + c);
        const auto *hk_148 = buffer.data(hk + 148 * ncomps + c);
        const auto *hk_149 = buffer.data(hk + 149 * ncomps + c);
        const auto *hk_150 = buffer.data(hk + 150 * ncomps + c);
        const auto *hk_151 = buffer.data(hk + 151 * ncomps + c);
        const auto *hk_152 = buffer.data(hk + 152 * ncomps + c);
        const auto *hk_153 = buffer.data(hk + 153 * ncomps + c);
        const auto *hk_154 = buffer.data(hk + 154 * ncomps + c);
        const auto *hk_155 = buffer.data(hk + 155 * ncomps + c);
        const auto *hk_156 = buffer.data(hk + 156 * ncomps + c);
        const auto *hk_157 = buffer.data(hk + 157 * ncomps + c);
        const auto *hk_158 = buffer.data(hk + 158 * ncomps + c);
        const auto *hk_159 = buffer.data(hk + 159 * ncomps + c);
        const auto *hk_160 = buffer.data(hk + 160 * ncomps + c);
        const auto *hk_161 = buffer.data(hk + 161 * ncomps + c);
        const auto *hk_162 = buffer.data(hk + 162 * ncomps + c);
        const auto *hk_163 = buffer.data(hk + 163 * ncomps + c);
        const auto *hk_164 = buffer.data(hk + 164 * ncomps + c);
        const auto *hk_165 = buffer.data(hk + 165 * ncomps + c);
        const auto *hk_166 = buffer.data(hk + 166 * ncomps + c);
        const auto *hk_167 = buffer.data(hk + 167 * ncomps + c);
        const auto *hk_168 = buffer.data(hk + 168 * ncomps + c);
        const auto *hk_169 = buffer.data(hk + 169 * ncomps + c);
        const auto *hk_170 = buffer.data(hk + 170 * ncomps + c);
        const auto *hk_171 = buffer.data(hk + 171 * ncomps + c);
        const auto *hk_180 = buffer.data(hk + 180 * ncomps + c);
        const auto *hk_181 = buffer.data(hk + 181 * ncomps + c);
        const auto *hk_182 = buffer.data(hk + 182 * ncomps + c);
        const auto *hk_183 = buffer.data(hk + 183 * ncomps + c);
        const auto *hk_184 = buffer.data(hk + 184 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, hi_0, hi_1, hi_2, hi_3, hi_4, hk_0, \
                         hk_1, hk_2, hk_3, hk_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * hi_0[k]
                     + hk_0[k];

            t_1[k] = -ab_x[k] * hi_1[k]
                     + hk_1[k];

            t_2[k] = -ab_x[k] * hi_2[k]
                     + hk_2[k];

            t_3[k] = -ab_x[k] * hi_3[k]
                     + hk_3[k];

            t_4[k] = -ab_x[k] * hi_4[k]
                     + hk_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, hi_5, hi_6, hi_7, hi_8, hi_9, hk_5, \
                         hk_6, hk_7, hk_8, hk_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * hi_5[k]
                     + hk_5[k];

            t_6[k] = -ab_x[k] * hi_6[k]
                     + hk_6[k];

            t_7[k] = -ab_x[k] * hi_7[k]
                     + hk_7[k];

            t_8[k] = -ab_x[k] * hi_8[k]
                     + hk_8[k];

            t_9[k] = -ab_x[k] * hi_9[k]
                     + hk_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, hi_10, hi_11, hi_12, hi_13, \
                         hi_14, hk_10, hk_11, hk_12, hk_13, hk_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * hi_10[k]
                      + hk_10[k];

            t_11[k] = -ab_x[k] * hi_11[k]
                      + hk_11[k];

            t_12[k] = -ab_x[k] * hi_12[k]
                      + hk_12[k];

            t_13[k] = -ab_x[k] * hi_13[k]
                      + hk_13[k];

            t_14[k] = -ab_x[k] * hi_14[k]
                      + hk_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, hi_15, hi_16, hi_17, hi_18, \
                         hi_19, hk_15, hk_16, hk_17, hk_18, hk_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * hi_15[k]
                      + hk_15[k];

            t_16[k] = -ab_x[k] * hi_16[k]
                      + hk_16[k];

            t_17[k] = -ab_x[k] * hi_17[k]
                      + hk_17[k];

            t_18[k] = -ab_x[k] * hi_18[k]
                      + hk_18[k];

            t_19[k] = -ab_x[k] * hi_19[k]
                      + hk_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, hi_20, hi_21, hi_22, hi_23, \
                         hi_24, hk_20, hk_21, hk_22, hk_23, hk_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * hi_20[k]
                      + hk_20[k];

            t_21[k] = -ab_x[k] * hi_21[k]
                      + hk_21[k];

            t_22[k] = -ab_x[k] * hi_22[k]
                      + hk_22[k];

            t_23[k] = -ab_x[k] * hi_23[k]
                      + hk_23[k];

            t_24[k] = -ab_x[k] * hi_24[k]
                      + hk_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, hi_25, hi_26, hi_27, hi_28, \
                         hi_29, hk_25, hk_26, hk_27, hk_36, hk_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * hi_25[k]
                      + hk_25[k];

            t_26[k] = -ab_x[k] * hi_26[k]
                      + hk_26[k];

            t_27[k] = -ab_x[k] * hi_27[k]
                      + hk_27[k];

            t_28[k] = -ab_x[k] * hi_28[k]
                      + hk_36[k];

            t_29[k] = -ab_x[k] * hi_29[k]
                      + hk_37[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, hi_30, hi_31, hi_32, hi_33, \
                         hi_34, hk_38, hk_39, hk_40, hk_41, hk_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * hi_30[k]
                      + hk_38[k];

            t_31[k] = -ab_x[k] * hi_31[k]
                      + hk_39[k];

            t_32[k] = -ab_x[k] * hi_32[k]
                      + hk_40[k];

            t_33[k] = -ab_x[k] * hi_33[k]
                      + hk_41[k];

            t_34[k] = -ab_x[k] * hi_34[k]
                      + hk_42[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, hi_35, hi_36, hi_37, hi_38, \
                         hi_39, hk_43, hk_44, hk_45, hk_46, hk_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * hi_35[k]
                      + hk_43[k];

            t_36[k] = -ab_x[k] * hi_36[k]
                      + hk_44[k];

            t_37[k] = -ab_x[k] * hi_37[k]
                      + hk_45[k];

            t_38[k] = -ab_x[k] * hi_38[k]
                      + hk_46[k];

            t_39[k] = -ab_x[k] * hi_39[k]
                      + hk_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, hi_40, hi_41, hi_42, hi_43, \
                         hi_44, hk_48, hk_49, hk_50, hk_51, hk_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * hi_40[k]
                      + hk_48[k];

            t_41[k] = -ab_x[k] * hi_41[k]
                      + hk_49[k];

            t_42[k] = -ab_x[k] * hi_42[k]
                      + hk_50[k];

            t_43[k] = -ab_x[k] * hi_43[k]
                      + hk_51[k];

            t_44[k] = -ab_x[k] * hi_44[k]
                      + hk_52[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, hi_45, hi_46, hi_47, hi_48, \
                         hi_49, hk_53, hk_54, hk_55, hk_56, hk_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * hi_45[k]
                      + hk_53[k];

            t_46[k] = -ab_x[k] * hi_46[k]
                      + hk_54[k];

            t_47[k] = -ab_x[k] * hi_47[k]
                      + hk_55[k];

            t_48[k] = -ab_x[k] * hi_48[k]
                      + hk_56[k];

            t_49[k] = -ab_x[k] * hi_49[k]
                      + hk_57[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, hi_50, hi_51, hi_52, hi_53, \
                         hi_54, hk_58, hk_59, hk_60, hk_61, hk_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * hi_50[k]
                      + hk_58[k];

            t_51[k] = -ab_x[k] * hi_51[k]
                      + hk_59[k];

            t_52[k] = -ab_x[k] * hi_52[k]
                      + hk_60[k];

            t_53[k] = -ab_x[k] * hi_53[k]
                      + hk_61[k];

            t_54[k] = -ab_x[k] * hi_54[k]
                      + hk_62[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, hi_55, hi_56, hi_57, hi_58, \
                         hi_59, hk_63, hk_72, hk_73, hk_74, hk_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * hi_55[k]
                      + hk_63[k];

            t_56[k] = -ab_x[k] * hi_56[k]
                      + hk_72[k];

            t_57[k] = -ab_x[k] * hi_57[k]
                      + hk_73[k];

            t_58[k] = -ab_x[k] * hi_58[k]
                      + hk_74[k];

            t_59[k] = -ab_x[k] * hi_59[k]
                      + hk_75[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, hi_60, hi_61, hi_62, hi_63, \
                         hi_64, hk_76, hk_77, hk_78, hk_79, hk_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * hi_60[k]
                      + hk_76[k];

            t_61[k] = -ab_x[k] * hi_61[k]
                      + hk_77[k];

            t_62[k] = -ab_x[k] * hi_62[k]
                      + hk_78[k];

            t_63[k] = -ab_x[k] * hi_63[k]
                      + hk_79[k];

            t_64[k] = -ab_x[k] * hi_64[k]
                      + hk_80[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, hi_65, hi_66, hi_67, hi_68, \
                         hi_69, hk_81, hk_82, hk_83, hk_84, hk_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * hi_65[k]
                      + hk_81[k];

            t_66[k] = -ab_x[k] * hi_66[k]
                      + hk_82[k];

            t_67[k] = -ab_x[k] * hi_67[k]
                      + hk_83[k];

            t_68[k] = -ab_x[k] * hi_68[k]
                      + hk_84[k];

            t_69[k] = -ab_x[k] * hi_69[k]
                      + hk_85[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, hi_70, hi_71, hi_72, hi_73, \
                         hi_74, hk_86, hk_87, hk_88, hk_89, hk_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * hi_70[k]
                      + hk_86[k];

            t_71[k] = -ab_x[k] * hi_71[k]
                      + hk_87[k];

            t_72[k] = -ab_x[k] * hi_72[k]
                      + hk_88[k];

            t_73[k] = -ab_x[k] * hi_73[k]
                      + hk_89[k];

            t_74[k] = -ab_x[k] * hi_74[k]
                      + hk_90[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, hi_75, hi_76, hi_77, hi_78, \
                         hi_79, hk_91, hk_92, hk_93, hk_94, hk_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * hi_75[k]
                      + hk_91[k];

            t_76[k] = -ab_x[k] * hi_76[k]
                      + hk_92[k];

            t_77[k] = -ab_x[k] * hi_77[k]
                      + hk_93[k];

            t_78[k] = -ab_x[k] * hi_78[k]
                      + hk_94[k];

            t_79[k] = -ab_x[k] * hi_79[k]
                      + hk_95[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, hi_80, hi_81, hi_82, hi_83, \
                         hi_84, hk_96, hk_97, hk_98, hk_99, hk_108 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * hi_80[k]
                      + hk_96[k];

            t_81[k] = -ab_x[k] * hi_81[k]
                      + hk_97[k];

            t_82[k] = -ab_x[k] * hi_82[k]
                      + hk_98[k];

            t_83[k] = -ab_x[k] * hi_83[k]
                      + hk_99[k];

            t_84[k] = -ab_x[k] * hi_84[k]
                      + hk_108[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, hi_85, hi_86, hi_87, hi_88, \
                         hi_89, hk_109, hk_110, hk_111, hk_112, \
                         hk_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * hi_85[k]
                      + hk_109[k];

            t_86[k] = -ab_x[k] * hi_86[k]
                      + hk_110[k];

            t_87[k] = -ab_x[k] * hi_87[k]
                      + hk_111[k];

            t_88[k] = -ab_x[k] * hi_88[k]
                      + hk_112[k];

            t_89[k] = -ab_x[k] * hi_89[k]
                      + hk_113[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, hi_90, hi_91, hi_92, hi_93, \
                         hi_94, hk_114, hk_115, hk_116, hk_117, \
                         hk_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * hi_90[k]
                      + hk_114[k];

            t_91[k] = -ab_x[k] * hi_91[k]
                      + hk_115[k];

            t_92[k] = -ab_x[k] * hi_92[k]
                      + hk_116[k];

            t_93[k] = -ab_x[k] * hi_93[k]
                      + hk_117[k];

            t_94[k] = -ab_x[k] * hi_94[k]
                      + hk_118[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, hi_95, hi_96, hi_97, hi_98, \
                         hi_99, hk_119, hk_120, hk_121, hk_122, \
                         hk_123 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * hi_95[k]
                      + hk_119[k];

            t_96[k] = -ab_x[k] * hi_96[k]
                      + hk_120[k];

            t_97[k] = -ab_x[k] * hi_97[k]
                      + hk_121[k];

            t_98[k] = -ab_x[k] * hi_98[k]
                      + hk_122[k];

            t_99[k] = -ab_x[k] * hi_99[k]
                      + hk_123[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, hi_100, hi_101, hi_102, \
                         hi_103, hi_104, hk_124, hk_125, hk_126, hk_127, \
                         hk_128 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * hi_100[k]
                       + hk_124[k];

            t_101[k] = -ab_x[k] * hi_101[k]
                       + hk_125[k];

            t_102[k] = -ab_x[k] * hi_102[k]
                       + hk_126[k];

            t_103[k] = -ab_x[k] * hi_103[k]
                       + hk_127[k];

            t_104[k] = -ab_x[k] * hi_104[k]
                       + hk_128[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, hi_105, hi_106, hi_107, \
                         hi_108, hi_109, hk_129, hk_130, hk_131, hk_132, \
                         hk_133 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * hi_105[k]
                       + hk_129[k];

            t_106[k] = -ab_x[k] * hi_106[k]
                       + hk_130[k];

            t_107[k] = -ab_x[k] * hi_107[k]
                       + hk_131[k];

            t_108[k] = -ab_x[k] * hi_108[k]
                       + hk_132[k];

            t_109[k] = -ab_x[k] * hi_109[k]
                       + hk_133[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, hi_110, hi_111, hi_112, \
                         hi_113, hi_114, hk_134, hk_135, hk_144, hk_145, \
                         hk_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * hi_110[k]
                       + hk_134[k];

            t_111[k] = -ab_x[k] * hi_111[k]
                       + hk_135[k];

            t_112[k] = -ab_x[k] * hi_112[k]
                       + hk_144[k];

            t_113[k] = -ab_x[k] * hi_113[k]
                       + hk_145[k];

            t_114[k] = -ab_x[k] * hi_114[k]
                       + hk_146[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, hi_115, hi_116, hi_117, \
                         hi_118, hi_119, hk_147, hk_148, hk_149, hk_150, \
                         hk_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * hi_115[k]
                       + hk_147[k];

            t_116[k] = -ab_x[k] * hi_116[k]
                       + hk_148[k];

            t_117[k] = -ab_x[k] * hi_117[k]
                       + hk_149[k];

            t_118[k] = -ab_x[k] * hi_118[k]
                       + hk_150[k];

            t_119[k] = -ab_x[k] * hi_119[k]
                       + hk_151[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, hi_120, hi_121, hi_122, \
                         hi_123, hi_124, hk_152, hk_153, hk_154, hk_155, \
                         hk_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * hi_120[k]
                       + hk_152[k];

            t_121[k] = -ab_x[k] * hi_121[k]
                       + hk_153[k];

            t_122[k] = -ab_x[k] * hi_122[k]
                       + hk_154[k];

            t_123[k] = -ab_x[k] * hi_123[k]
                       + hk_155[k];

            t_124[k] = -ab_x[k] * hi_124[k]
                       + hk_156[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, hi_125, hi_126, hi_127, \
                         hi_128, hi_129, hk_157, hk_158, hk_159, hk_160, \
                         hk_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * hi_125[k]
                       + hk_157[k];

            t_126[k] = -ab_x[k] * hi_126[k]
                       + hk_158[k];

            t_127[k] = -ab_x[k] * hi_127[k]
                       + hk_159[k];

            t_128[k] = -ab_x[k] * hi_128[k]
                       + hk_160[k];

            t_129[k] = -ab_x[k] * hi_129[k]
                       + hk_161[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, hi_130, hi_131, hi_132, \
                         hi_133, hi_134, hk_162, hk_163, hk_164, hk_165, \
                         hk_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * hi_130[k]
                       + hk_162[k];

            t_131[k] = -ab_x[k] * hi_131[k]
                       + hk_163[k];

            t_132[k] = -ab_x[k] * hi_132[k]
                       + hk_164[k];

            t_133[k] = -ab_x[k] * hi_133[k]
                       + hk_165[k];

            t_134[k] = -ab_x[k] * hi_134[k]
                       + hk_166[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, hi_135, hi_136, hi_137, \
                         hi_138, hi_139, hk_167, hk_168, hk_169, hk_170, \
                         hk_171 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * hi_135[k]
                       + hk_167[k];

            t_136[k] = -ab_x[k] * hi_136[k]
                       + hk_168[k];

            t_137[k] = -ab_x[k] * hi_137[k]
                       + hk_169[k];

            t_138[k] = -ab_x[k] * hi_138[k]
                       + hk_170[k];

            t_139[k] = -ab_x[k] * hi_139[k]
                       + hk_171[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, hi_140, hi_141, hi_142, \
                         hi_143, hi_144, hk_180, hk_181, hk_182, hk_183, \
                         hk_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * hi_140[k]
                       + hk_180[k];

            t_141[k] = -ab_x[k] * hi_141[k]
                       + hk_181[k];

            t_142[k] = -ab_x[k] * hi_142[k]
                       + hk_182[k];

            t_143[k] = -ab_x[k] * hi_143[k]
                       + hk_183[k];

            t_144[k] = -ab_x[k] * hi_144[k]
                       + hk_184[k];
        }
    }
}

static auto
compute_hrr_ii_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hi, const size_t hk, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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
        auto *t_234 = buffer.data(target + 234 * ncomps + c);
        auto *t_235 = buffer.data(target + 235 * ncomps + c);
        auto *t_236 = buffer.data(target + 236 * ncomps + c);
        auto *t_237 = buffer.data(target + 237 * ncomps + c);
        auto *t_238 = buffer.data(target + 238 * ncomps + c);
        auto *t_239 = buffer.data(target + 239 * ncomps + c);
        auto *t_240 = buffer.data(target + 240 * ncomps + c);
        auto *t_241 = buffer.data(target + 241 * ncomps + c);
        auto *t_242 = buffer.data(target + 242 * ncomps + c);
        auto *t_243 = buffer.data(target + 243 * ncomps + c);
        auto *t_244 = buffer.data(target + 244 * ncomps + c);
        auto *t_245 = buffer.data(target + 245 * ncomps + c);
        auto *t_246 = buffer.data(target + 246 * ncomps + c);
        auto *t_247 = buffer.data(target + 247 * ncomps + c);
        auto *t_248 = buffer.data(target + 248 * ncomps + c);
        auto *t_249 = buffer.data(target + 249 * ncomps + c);
        auto *t_250 = buffer.data(target + 250 * ncomps + c);
        auto *t_251 = buffer.data(target + 251 * ncomps + c);
        auto *t_252 = buffer.data(target + 252 * ncomps + c);
        auto *t_253 = buffer.data(target + 253 * ncomps + c);
        auto *t_254 = buffer.data(target + 254 * ncomps + c);
        auto *t_255 = buffer.data(target + 255 * ncomps + c);
        auto *t_256 = buffer.data(target + 256 * ncomps + c);
        auto *t_257 = buffer.data(target + 257 * ncomps + c);
        auto *t_258 = buffer.data(target + 258 * ncomps + c);
        auto *t_259 = buffer.data(target + 259 * ncomps + c);
        auto *t_260 = buffer.data(target + 260 * ncomps + c);
        auto *t_261 = buffer.data(target + 261 * ncomps + c);
        auto *t_262 = buffer.data(target + 262 * ncomps + c);
        auto *t_263 = buffer.data(target + 263 * ncomps + c);
        auto *t_264 = buffer.data(target + 264 * ncomps + c);
        auto *t_265 = buffer.data(target + 265 * ncomps + c);
        auto *t_266 = buffer.data(target + 266 * ncomps + c);
        auto *t_267 = buffer.data(target + 267 * ncomps + c);
        auto *t_268 = buffer.data(target + 268 * ncomps + c);
        auto *t_269 = buffer.data(target + 269 * ncomps + c);
        auto *t_270 = buffer.data(target + 270 * ncomps + c);
        auto *t_271 = buffer.data(target + 271 * ncomps + c);
        auto *t_272 = buffer.data(target + 272 * ncomps + c);
        auto *t_273 = buffer.data(target + 273 * ncomps + c);
        auto *t_274 = buffer.data(target + 274 * ncomps + c);
        auto *t_275 = buffer.data(target + 275 * ncomps + c);
        auto *t_276 = buffer.data(target + 276 * ncomps + c);
        auto *t_277 = buffer.data(target + 277 * ncomps + c);
        auto *t_278 = buffer.data(target + 278 * ncomps + c);
        auto *t_279 = buffer.data(target + 279 * ncomps + c);
        auto *t_280 = buffer.data(target + 280 * ncomps + c);
        auto *t_281 = buffer.data(target + 281 * ncomps + c);
        auto *t_282 = buffer.data(target + 282 * ncomps + c);
        auto *t_283 = buffer.data(target + 283 * ncomps + c);
        auto *t_284 = buffer.data(target + 284 * ncomps + c);
        auto *t_285 = buffer.data(target + 285 * ncomps + c);
        auto *t_286 = buffer.data(target + 286 * ncomps + c);
        auto *t_287 = buffer.data(target + 287 * ncomps + c);
        auto *t_288 = buffer.data(target + 288 * ncomps + c);
        auto *t_289 = buffer.data(target + 289 * ncomps + c);

        const auto *ab_x = coordinates.data(6);

        const auto *hi_145 = buffer.data(hi + 145 * ncomps + c);
        const auto *hi_146 = buffer.data(hi + 146 * ncomps + c);
        const auto *hi_147 = buffer.data(hi + 147 * ncomps + c);
        const auto *hi_148 = buffer.data(hi + 148 * ncomps + c);
        const auto *hi_149 = buffer.data(hi + 149 * ncomps + c);
        const auto *hi_150 = buffer.data(hi + 150 * ncomps + c);
        const auto *hi_151 = buffer.data(hi + 151 * ncomps + c);
        const auto *hi_152 = buffer.data(hi + 152 * ncomps + c);
        const auto *hi_153 = buffer.data(hi + 153 * ncomps + c);
        const auto *hi_154 = buffer.data(hi + 154 * ncomps + c);
        const auto *hi_155 = buffer.data(hi + 155 * ncomps + c);
        const auto *hi_156 = buffer.data(hi + 156 * ncomps + c);
        const auto *hi_157 = buffer.data(hi + 157 * ncomps + c);
        const auto *hi_158 = buffer.data(hi + 158 * ncomps + c);
        const auto *hi_159 = buffer.data(hi + 159 * ncomps + c);
        const auto *hi_160 = buffer.data(hi + 160 * ncomps + c);
        const auto *hi_161 = buffer.data(hi + 161 * ncomps + c);
        const auto *hi_162 = buffer.data(hi + 162 * ncomps + c);
        const auto *hi_163 = buffer.data(hi + 163 * ncomps + c);
        const auto *hi_164 = buffer.data(hi + 164 * ncomps + c);
        const auto *hi_165 = buffer.data(hi + 165 * ncomps + c);
        const auto *hi_166 = buffer.data(hi + 166 * ncomps + c);
        const auto *hi_167 = buffer.data(hi + 167 * ncomps + c);
        const auto *hi_168 = buffer.data(hi + 168 * ncomps + c);
        const auto *hi_169 = buffer.data(hi + 169 * ncomps + c);
        const auto *hi_170 = buffer.data(hi + 170 * ncomps + c);
        const auto *hi_171 = buffer.data(hi + 171 * ncomps + c);
        const auto *hi_172 = buffer.data(hi + 172 * ncomps + c);
        const auto *hi_173 = buffer.data(hi + 173 * ncomps + c);
        const auto *hi_174 = buffer.data(hi + 174 * ncomps + c);
        const auto *hi_175 = buffer.data(hi + 175 * ncomps + c);
        const auto *hi_176 = buffer.data(hi + 176 * ncomps + c);
        const auto *hi_177 = buffer.data(hi + 177 * ncomps + c);
        const auto *hi_178 = buffer.data(hi + 178 * ncomps + c);
        const auto *hi_179 = buffer.data(hi + 179 * ncomps + c);
        const auto *hi_180 = buffer.data(hi + 180 * ncomps + c);
        const auto *hi_181 = buffer.data(hi + 181 * ncomps + c);
        const auto *hi_182 = buffer.data(hi + 182 * ncomps + c);
        const auto *hi_183 = buffer.data(hi + 183 * ncomps + c);
        const auto *hi_184 = buffer.data(hi + 184 * ncomps + c);
        const auto *hi_185 = buffer.data(hi + 185 * ncomps + c);
        const auto *hi_186 = buffer.data(hi + 186 * ncomps + c);
        const auto *hi_187 = buffer.data(hi + 187 * ncomps + c);
        const auto *hi_188 = buffer.data(hi + 188 * ncomps + c);
        const auto *hi_189 = buffer.data(hi + 189 * ncomps + c);
        const auto *hi_190 = buffer.data(hi + 190 * ncomps + c);
        const auto *hi_191 = buffer.data(hi + 191 * ncomps + c);
        const auto *hi_192 = buffer.data(hi + 192 * ncomps + c);
        const auto *hi_193 = buffer.data(hi + 193 * ncomps + c);
        const auto *hi_194 = buffer.data(hi + 194 * ncomps + c);
        const auto *hi_195 = buffer.data(hi + 195 * ncomps + c);
        const auto *hi_196 = buffer.data(hi + 196 * ncomps + c);
        const auto *hi_197 = buffer.data(hi + 197 * ncomps + c);
        const auto *hi_198 = buffer.data(hi + 198 * ncomps + c);
        const auto *hi_199 = buffer.data(hi + 199 * ncomps + c);
        const auto *hi_200 = buffer.data(hi + 200 * ncomps + c);
        const auto *hi_201 = buffer.data(hi + 201 * ncomps + c);
        const auto *hi_202 = buffer.data(hi + 202 * ncomps + c);
        const auto *hi_203 = buffer.data(hi + 203 * ncomps + c);
        const auto *hi_204 = buffer.data(hi + 204 * ncomps + c);
        const auto *hi_205 = buffer.data(hi + 205 * ncomps + c);
        const auto *hi_206 = buffer.data(hi + 206 * ncomps + c);
        const auto *hi_207 = buffer.data(hi + 207 * ncomps + c);
        const auto *hi_208 = buffer.data(hi + 208 * ncomps + c);
        const auto *hi_209 = buffer.data(hi + 209 * ncomps + c);
        const auto *hi_210 = buffer.data(hi + 210 * ncomps + c);
        const auto *hi_211 = buffer.data(hi + 211 * ncomps + c);
        const auto *hi_212 = buffer.data(hi + 212 * ncomps + c);
        const auto *hi_213 = buffer.data(hi + 213 * ncomps + c);
        const auto *hi_214 = buffer.data(hi + 214 * ncomps + c);
        const auto *hi_215 = buffer.data(hi + 215 * ncomps + c);
        const auto *hi_216 = buffer.data(hi + 216 * ncomps + c);
        const auto *hi_217 = buffer.data(hi + 217 * ncomps + c);
        const auto *hi_218 = buffer.data(hi + 218 * ncomps + c);
        const auto *hi_219 = buffer.data(hi + 219 * ncomps + c);
        const auto *hi_220 = buffer.data(hi + 220 * ncomps + c);
        const auto *hi_221 = buffer.data(hi + 221 * ncomps + c);
        const auto *hi_222 = buffer.data(hi + 222 * ncomps + c);
        const auto *hi_223 = buffer.data(hi + 223 * ncomps + c);
        const auto *hi_224 = buffer.data(hi + 224 * ncomps + c);
        const auto *hi_225 = buffer.data(hi + 225 * ncomps + c);
        const auto *hi_226 = buffer.data(hi + 226 * ncomps + c);
        const auto *hi_227 = buffer.data(hi + 227 * ncomps + c);
        const auto *hi_228 = buffer.data(hi + 228 * ncomps + c);
        const auto *hi_229 = buffer.data(hi + 229 * ncomps + c);
        const auto *hi_230 = buffer.data(hi + 230 * ncomps + c);
        const auto *hi_231 = buffer.data(hi + 231 * ncomps + c);
        const auto *hi_232 = buffer.data(hi + 232 * ncomps + c);
        const auto *hi_233 = buffer.data(hi + 233 * ncomps + c);
        const auto *hi_234 = buffer.data(hi + 234 * ncomps + c);
        const auto *hi_235 = buffer.data(hi + 235 * ncomps + c);
        const auto *hi_236 = buffer.data(hi + 236 * ncomps + c);
        const auto *hi_237 = buffer.data(hi + 237 * ncomps + c);
        const auto *hi_238 = buffer.data(hi + 238 * ncomps + c);
        const auto *hi_239 = buffer.data(hi + 239 * ncomps + c);
        const auto *hi_240 = buffer.data(hi + 240 * ncomps + c);
        const auto *hi_241 = buffer.data(hi + 241 * ncomps + c);
        const auto *hi_242 = buffer.data(hi + 242 * ncomps + c);
        const auto *hi_243 = buffer.data(hi + 243 * ncomps + c);
        const auto *hi_244 = buffer.data(hi + 244 * ncomps + c);
        const auto *hi_245 = buffer.data(hi + 245 * ncomps + c);
        const auto *hi_246 = buffer.data(hi + 246 * ncomps + c);
        const auto *hi_247 = buffer.data(hi + 247 * ncomps + c);
        const auto *hi_248 = buffer.data(hi + 248 * ncomps + c);
        const auto *hi_249 = buffer.data(hi + 249 * ncomps + c);
        const auto *hi_250 = buffer.data(hi + 250 * ncomps + c);
        const auto *hi_251 = buffer.data(hi + 251 * ncomps + c);
        const auto *hi_252 = buffer.data(hi + 252 * ncomps + c);
        const auto *hi_253 = buffer.data(hi + 253 * ncomps + c);
        const auto *hi_254 = buffer.data(hi + 254 * ncomps + c);
        const auto *hi_255 = buffer.data(hi + 255 * ncomps + c);
        const auto *hi_256 = buffer.data(hi + 256 * ncomps + c);
        const auto *hi_257 = buffer.data(hi + 257 * ncomps + c);
        const auto *hi_258 = buffer.data(hi + 258 * ncomps + c);
        const auto *hi_259 = buffer.data(hi + 259 * ncomps + c);
        const auto *hi_260 = buffer.data(hi + 260 * ncomps + c);
        const auto *hi_261 = buffer.data(hi + 261 * ncomps + c);
        const auto *hi_262 = buffer.data(hi + 262 * ncomps + c);
        const auto *hi_263 = buffer.data(hi + 263 * ncomps + c);
        const auto *hi_264 = buffer.data(hi + 264 * ncomps + c);
        const auto *hi_265 = buffer.data(hi + 265 * ncomps + c);
        const auto *hi_266 = buffer.data(hi + 266 * ncomps + c);
        const auto *hi_267 = buffer.data(hi + 267 * ncomps + c);
        const auto *hi_268 = buffer.data(hi + 268 * ncomps + c);
        const auto *hi_269 = buffer.data(hi + 269 * ncomps + c);
        const auto *hi_270 = buffer.data(hi + 270 * ncomps + c);
        const auto *hi_271 = buffer.data(hi + 271 * ncomps + c);
        const auto *hi_272 = buffer.data(hi + 272 * ncomps + c);
        const auto *hi_273 = buffer.data(hi + 273 * ncomps + c);
        const auto *hi_274 = buffer.data(hi + 274 * ncomps + c);
        const auto *hi_275 = buffer.data(hi + 275 * ncomps + c);
        const auto *hi_276 = buffer.data(hi + 276 * ncomps + c);
        const auto *hi_277 = buffer.data(hi + 277 * ncomps + c);
        const auto *hi_278 = buffer.data(hi + 278 * ncomps + c);
        const auto *hi_279 = buffer.data(hi + 279 * ncomps + c);
        const auto *hi_280 = buffer.data(hi + 280 * ncomps + c);
        const auto *hi_281 = buffer.data(hi + 281 * ncomps + c);
        const auto *hi_282 = buffer.data(hi + 282 * ncomps + c);
        const auto *hi_283 = buffer.data(hi + 283 * ncomps + c);
        const auto *hi_284 = buffer.data(hi + 284 * ncomps + c);
        const auto *hi_285 = buffer.data(hi + 285 * ncomps + c);
        const auto *hi_286 = buffer.data(hi + 286 * ncomps + c);
        const auto *hi_287 = buffer.data(hi + 287 * ncomps + c);
        const auto *hi_288 = buffer.data(hi + 288 * ncomps + c);
        const auto *hi_289 = buffer.data(hi + 289 * ncomps + c);

        const auto *hk_185 = buffer.data(hk + 185 * ncomps + c);
        const auto *hk_186 = buffer.data(hk + 186 * ncomps + c);
        const auto *hk_187 = buffer.data(hk + 187 * ncomps + c);
        const auto *hk_188 = buffer.data(hk + 188 * ncomps + c);
        const auto *hk_189 = buffer.data(hk + 189 * ncomps + c);
        const auto *hk_190 = buffer.data(hk + 190 * ncomps + c);
        const auto *hk_191 = buffer.data(hk + 191 * ncomps + c);
        const auto *hk_192 = buffer.data(hk + 192 * ncomps + c);
        const auto *hk_193 = buffer.data(hk + 193 * ncomps + c);
        const auto *hk_194 = buffer.data(hk + 194 * ncomps + c);
        const auto *hk_195 = buffer.data(hk + 195 * ncomps + c);
        const auto *hk_196 = buffer.data(hk + 196 * ncomps + c);
        const auto *hk_197 = buffer.data(hk + 197 * ncomps + c);
        const auto *hk_198 = buffer.data(hk + 198 * ncomps + c);
        const auto *hk_199 = buffer.data(hk + 199 * ncomps + c);
        const auto *hk_200 = buffer.data(hk + 200 * ncomps + c);
        const auto *hk_201 = buffer.data(hk + 201 * ncomps + c);
        const auto *hk_202 = buffer.data(hk + 202 * ncomps + c);
        const auto *hk_203 = buffer.data(hk + 203 * ncomps + c);
        const auto *hk_204 = buffer.data(hk + 204 * ncomps + c);
        const auto *hk_205 = buffer.data(hk + 205 * ncomps + c);
        const auto *hk_206 = buffer.data(hk + 206 * ncomps + c);
        const auto *hk_207 = buffer.data(hk + 207 * ncomps + c);
        const auto *hk_216 = buffer.data(hk + 216 * ncomps + c);
        const auto *hk_217 = buffer.data(hk + 217 * ncomps + c);
        const auto *hk_218 = buffer.data(hk + 218 * ncomps + c);
        const auto *hk_219 = buffer.data(hk + 219 * ncomps + c);
        const auto *hk_220 = buffer.data(hk + 220 * ncomps + c);
        const auto *hk_221 = buffer.data(hk + 221 * ncomps + c);
        const auto *hk_222 = buffer.data(hk + 222 * ncomps + c);
        const auto *hk_223 = buffer.data(hk + 223 * ncomps + c);
        const auto *hk_224 = buffer.data(hk + 224 * ncomps + c);
        const auto *hk_225 = buffer.data(hk + 225 * ncomps + c);
        const auto *hk_226 = buffer.data(hk + 226 * ncomps + c);
        const auto *hk_227 = buffer.data(hk + 227 * ncomps + c);
        const auto *hk_228 = buffer.data(hk + 228 * ncomps + c);
        const auto *hk_229 = buffer.data(hk + 229 * ncomps + c);
        const auto *hk_230 = buffer.data(hk + 230 * ncomps + c);
        const auto *hk_231 = buffer.data(hk + 231 * ncomps + c);
        const auto *hk_232 = buffer.data(hk + 232 * ncomps + c);
        const auto *hk_233 = buffer.data(hk + 233 * ncomps + c);
        const auto *hk_234 = buffer.data(hk + 234 * ncomps + c);
        const auto *hk_235 = buffer.data(hk + 235 * ncomps + c);
        const auto *hk_236 = buffer.data(hk + 236 * ncomps + c);
        const auto *hk_237 = buffer.data(hk + 237 * ncomps + c);
        const auto *hk_238 = buffer.data(hk + 238 * ncomps + c);
        const auto *hk_239 = buffer.data(hk + 239 * ncomps + c);
        const auto *hk_240 = buffer.data(hk + 240 * ncomps + c);
        const auto *hk_241 = buffer.data(hk + 241 * ncomps + c);
        const auto *hk_242 = buffer.data(hk + 242 * ncomps + c);
        const auto *hk_243 = buffer.data(hk + 243 * ncomps + c);
        const auto *hk_252 = buffer.data(hk + 252 * ncomps + c);
        const auto *hk_253 = buffer.data(hk + 253 * ncomps + c);
        const auto *hk_254 = buffer.data(hk + 254 * ncomps + c);
        const auto *hk_255 = buffer.data(hk + 255 * ncomps + c);
        const auto *hk_256 = buffer.data(hk + 256 * ncomps + c);
        const auto *hk_257 = buffer.data(hk + 257 * ncomps + c);
        const auto *hk_258 = buffer.data(hk + 258 * ncomps + c);
        const auto *hk_259 = buffer.data(hk + 259 * ncomps + c);
        const auto *hk_260 = buffer.data(hk + 260 * ncomps + c);
        const auto *hk_261 = buffer.data(hk + 261 * ncomps + c);
        const auto *hk_262 = buffer.data(hk + 262 * ncomps + c);
        const auto *hk_263 = buffer.data(hk + 263 * ncomps + c);
        const auto *hk_264 = buffer.data(hk + 264 * ncomps + c);
        const auto *hk_265 = buffer.data(hk + 265 * ncomps + c);
        const auto *hk_266 = buffer.data(hk + 266 * ncomps + c);
        const auto *hk_267 = buffer.data(hk + 267 * ncomps + c);
        const auto *hk_268 = buffer.data(hk + 268 * ncomps + c);
        const auto *hk_269 = buffer.data(hk + 269 * ncomps + c);
        const auto *hk_270 = buffer.data(hk + 270 * ncomps + c);
        const auto *hk_271 = buffer.data(hk + 271 * ncomps + c);
        const auto *hk_272 = buffer.data(hk + 272 * ncomps + c);
        const auto *hk_273 = buffer.data(hk + 273 * ncomps + c);
        const auto *hk_274 = buffer.data(hk + 274 * ncomps + c);
        const auto *hk_275 = buffer.data(hk + 275 * ncomps + c);
        const auto *hk_276 = buffer.data(hk + 276 * ncomps + c);
        const auto *hk_277 = buffer.data(hk + 277 * ncomps + c);
        const auto *hk_278 = buffer.data(hk + 278 * ncomps + c);
        const auto *hk_279 = buffer.data(hk + 279 * ncomps + c);
        const auto *hk_288 = buffer.data(hk + 288 * ncomps + c);
        const auto *hk_289 = buffer.data(hk + 289 * ncomps + c);
        const auto *hk_290 = buffer.data(hk + 290 * ncomps + c);
        const auto *hk_291 = buffer.data(hk + 291 * ncomps + c);
        const auto *hk_292 = buffer.data(hk + 292 * ncomps + c);
        const auto *hk_293 = buffer.data(hk + 293 * ncomps + c);
        const auto *hk_294 = buffer.data(hk + 294 * ncomps + c);
        const auto *hk_295 = buffer.data(hk + 295 * ncomps + c);
        const auto *hk_296 = buffer.data(hk + 296 * ncomps + c);
        const auto *hk_297 = buffer.data(hk + 297 * ncomps + c);
        const auto *hk_298 = buffer.data(hk + 298 * ncomps + c);
        const auto *hk_299 = buffer.data(hk + 299 * ncomps + c);
        const auto *hk_300 = buffer.data(hk + 300 * ncomps + c);
        const auto *hk_301 = buffer.data(hk + 301 * ncomps + c);
        const auto *hk_302 = buffer.data(hk + 302 * ncomps + c);
        const auto *hk_303 = buffer.data(hk + 303 * ncomps + c);
        const auto *hk_304 = buffer.data(hk + 304 * ncomps + c);
        const auto *hk_305 = buffer.data(hk + 305 * ncomps + c);
        const auto *hk_306 = buffer.data(hk + 306 * ncomps + c);
        const auto *hk_307 = buffer.data(hk + 307 * ncomps + c);
        const auto *hk_308 = buffer.data(hk + 308 * ncomps + c);
        const auto *hk_309 = buffer.data(hk + 309 * ncomps + c);
        const auto *hk_310 = buffer.data(hk + 310 * ncomps + c);
        const auto *hk_311 = buffer.data(hk + 311 * ncomps + c);
        const auto *hk_312 = buffer.data(hk + 312 * ncomps + c);
        const auto *hk_313 = buffer.data(hk + 313 * ncomps + c);
        const auto *hk_314 = buffer.data(hk + 314 * ncomps + c);
        const auto *hk_315 = buffer.data(hk + 315 * ncomps + c);
        const auto *hk_324 = buffer.data(hk + 324 * ncomps + c);
        const auto *hk_325 = buffer.data(hk + 325 * ncomps + c);
        const auto *hk_326 = buffer.data(hk + 326 * ncomps + c);
        const auto *hk_327 = buffer.data(hk + 327 * ncomps + c);
        const auto *hk_328 = buffer.data(hk + 328 * ncomps + c);
        const auto *hk_329 = buffer.data(hk + 329 * ncomps + c);
        const auto *hk_330 = buffer.data(hk + 330 * ncomps + c);
        const auto *hk_331 = buffer.data(hk + 331 * ncomps + c);
        const auto *hk_332 = buffer.data(hk + 332 * ncomps + c);
        const auto *hk_333 = buffer.data(hk + 333 * ncomps + c);
        const auto *hk_334 = buffer.data(hk + 334 * ncomps + c);
        const auto *hk_335 = buffer.data(hk + 335 * ncomps + c);
        const auto *hk_336 = buffer.data(hk + 336 * ncomps + c);
        const auto *hk_337 = buffer.data(hk + 337 * ncomps + c);
        const auto *hk_338 = buffer.data(hk + 338 * ncomps + c);
        const auto *hk_339 = buffer.data(hk + 339 * ncomps + c);
        const auto *hk_340 = buffer.data(hk + 340 * ncomps + c);
        const auto *hk_341 = buffer.data(hk + 341 * ncomps + c);
        const auto *hk_342 = buffer.data(hk + 342 * ncomps + c);
        const auto *hk_343 = buffer.data(hk + 343 * ncomps + c);
        const auto *hk_344 = buffer.data(hk + 344 * ncomps + c);
        const auto *hk_345 = buffer.data(hk + 345 * ncomps + c);
        const auto *hk_346 = buffer.data(hk + 346 * ncomps + c);
        const auto *hk_347 = buffer.data(hk + 347 * ncomps + c);
        const auto *hk_348 = buffer.data(hk + 348 * ncomps + c);
        const auto *hk_349 = buffer.data(hk + 349 * ncomps + c);
        const auto *hk_350 = buffer.data(hk + 350 * ncomps + c);
        const auto *hk_351 = buffer.data(hk + 351 * ncomps + c);
        const auto *hk_360 = buffer.data(hk + 360 * ncomps + c);
        const auto *hk_361 = buffer.data(hk + 361 * ncomps + c);
        const auto *hk_362 = buffer.data(hk + 362 * ncomps + c);
        const auto *hk_363 = buffer.data(hk + 363 * ncomps + c);
        const auto *hk_364 = buffer.data(hk + 364 * ncomps + c);
        const auto *hk_365 = buffer.data(hk + 365 * ncomps + c);
        const auto *hk_366 = buffer.data(hk + 366 * ncomps + c);
        const auto *hk_367 = buffer.data(hk + 367 * ncomps + c);
        const auto *hk_368 = buffer.data(hk + 368 * ncomps + c);
        const auto *hk_369 = buffer.data(hk + 369 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, hi_145, hi_146, hi_147, \
                         hi_148, hi_149, hk_185, hk_186, hk_187, hk_188, \
                         hk_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * hi_145[k]
                       + hk_185[k];

            t_146[k] = -ab_x[k] * hi_146[k]
                       + hk_186[k];

            t_147[k] = -ab_x[k] * hi_147[k]
                       + hk_187[k];

            t_148[k] = -ab_x[k] * hi_148[k]
                       + hk_188[k];

            t_149[k] = -ab_x[k] * hi_149[k]
                       + hk_189[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, hi_150, hi_151, hi_152, \
                         hi_153, hi_154, hk_190, hk_191, hk_192, hk_193, \
                         hk_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * hi_150[k]
                       + hk_190[k];

            t_151[k] = -ab_x[k] * hi_151[k]
                       + hk_191[k];

            t_152[k] = -ab_x[k] * hi_152[k]
                       + hk_192[k];

            t_153[k] = -ab_x[k] * hi_153[k]
                       + hk_193[k];

            t_154[k] = -ab_x[k] * hi_154[k]
                       + hk_194[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, hi_155, hi_156, hi_157, \
                         hi_158, hi_159, hk_195, hk_196, hk_197, hk_198, \
                         hk_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * hi_155[k]
                       + hk_195[k];

            t_156[k] = -ab_x[k] * hi_156[k]
                       + hk_196[k];

            t_157[k] = -ab_x[k] * hi_157[k]
                       + hk_197[k];

            t_158[k] = -ab_x[k] * hi_158[k]
                       + hk_198[k];

            t_159[k] = -ab_x[k] * hi_159[k]
                       + hk_199[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, hi_160, hi_161, hi_162, \
                         hi_163, hi_164, hk_200, hk_201, hk_202, hk_203, \
                         hk_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * hi_160[k]
                       + hk_200[k];

            t_161[k] = -ab_x[k] * hi_161[k]
                       + hk_201[k];

            t_162[k] = -ab_x[k] * hi_162[k]
                       + hk_202[k];

            t_163[k] = -ab_x[k] * hi_163[k]
                       + hk_203[k];

            t_164[k] = -ab_x[k] * hi_164[k]
                       + hk_204[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, hi_165, hi_166, hi_167, \
                         hi_168, hi_169, hk_205, hk_206, hk_207, hk_216, \
                         hk_217 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * hi_165[k]
                       + hk_205[k];

            t_166[k] = -ab_x[k] * hi_166[k]
                       + hk_206[k];

            t_167[k] = -ab_x[k] * hi_167[k]
                       + hk_207[k];

            t_168[k] = -ab_x[k] * hi_168[k]
                       + hk_216[k];

            t_169[k] = -ab_x[k] * hi_169[k]
                       + hk_217[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, hi_170, hi_171, hi_172, \
                         hi_173, hi_174, hk_218, hk_219, hk_220, hk_221, \
                         hk_222 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * hi_170[k]
                       + hk_218[k];

            t_171[k] = -ab_x[k] * hi_171[k]
                       + hk_219[k];

            t_172[k] = -ab_x[k] * hi_172[k]
                       + hk_220[k];

            t_173[k] = -ab_x[k] * hi_173[k]
                       + hk_221[k];

            t_174[k] = -ab_x[k] * hi_174[k]
                       + hk_222[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, hi_175, hi_176, hi_177, \
                         hi_178, hi_179, hk_223, hk_224, hk_225, hk_226, \
                         hk_227 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * hi_175[k]
                       + hk_223[k];

            t_176[k] = -ab_x[k] * hi_176[k]
                       + hk_224[k];

            t_177[k] = -ab_x[k] * hi_177[k]
                       + hk_225[k];

            t_178[k] = -ab_x[k] * hi_178[k]
                       + hk_226[k];

            t_179[k] = -ab_x[k] * hi_179[k]
                       + hk_227[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, hi_180, hi_181, hi_182, \
                         hi_183, hi_184, hk_228, hk_229, hk_230, hk_231, \
                         hk_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * hi_180[k]
                       + hk_228[k];

            t_181[k] = -ab_x[k] * hi_181[k]
                       + hk_229[k];

            t_182[k] = -ab_x[k] * hi_182[k]
                       + hk_230[k];

            t_183[k] = -ab_x[k] * hi_183[k]
                       + hk_231[k];

            t_184[k] = -ab_x[k] * hi_184[k]
                       + hk_232[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, hi_185, hi_186, hi_187, \
                         hi_188, hi_189, hk_233, hk_234, hk_235, hk_236, \
                         hk_237 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * hi_185[k]
                       + hk_233[k];

            t_186[k] = -ab_x[k] * hi_186[k]
                       + hk_234[k];

            t_187[k] = -ab_x[k] * hi_187[k]
                       + hk_235[k];

            t_188[k] = -ab_x[k] * hi_188[k]
                       + hk_236[k];

            t_189[k] = -ab_x[k] * hi_189[k]
                       + hk_237[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, hi_190, hi_191, hi_192, \
                         hi_193, hi_194, hk_238, hk_239, hk_240, hk_241, \
                         hk_242 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * hi_190[k]
                       + hk_238[k];

            t_191[k] = -ab_x[k] * hi_191[k]
                       + hk_239[k];

            t_192[k] = -ab_x[k] * hi_192[k]
                       + hk_240[k];

            t_193[k] = -ab_x[k] * hi_193[k]
                       + hk_241[k];

            t_194[k] = -ab_x[k] * hi_194[k]
                       + hk_242[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, hi_195, hi_196, hi_197, \
                         hi_198, hi_199, hk_243, hk_252, hk_253, hk_254, \
                         hk_255 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * hi_195[k]
                       + hk_243[k];

            t_196[k] = -ab_x[k] * hi_196[k]
                       + hk_252[k];

            t_197[k] = -ab_x[k] * hi_197[k]
                       + hk_253[k];

            t_198[k] = -ab_x[k] * hi_198[k]
                       + hk_254[k];

            t_199[k] = -ab_x[k] * hi_199[k]
                       + hk_255[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, hi_200, hi_201, hi_202, \
                         hi_203, hi_204, hk_256, hk_257, hk_258, hk_259, \
                         hk_260 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * hi_200[k]
                       + hk_256[k];

            t_201[k] = -ab_x[k] * hi_201[k]
                       + hk_257[k];

            t_202[k] = -ab_x[k] * hi_202[k]
                       + hk_258[k];

            t_203[k] = -ab_x[k] * hi_203[k]
                       + hk_259[k];

            t_204[k] = -ab_x[k] * hi_204[k]
                       + hk_260[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, hi_205, hi_206, hi_207, \
                         hi_208, hi_209, hk_261, hk_262, hk_263, hk_264, \
                         hk_265 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * hi_205[k]
                       + hk_261[k];

            t_206[k] = -ab_x[k] * hi_206[k]
                       + hk_262[k];

            t_207[k] = -ab_x[k] * hi_207[k]
                       + hk_263[k];

            t_208[k] = -ab_x[k] * hi_208[k]
                       + hk_264[k];

            t_209[k] = -ab_x[k] * hi_209[k]
                       + hk_265[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, hi_210, hi_211, hi_212, \
                         hi_213, hi_214, hk_266, hk_267, hk_268, hk_269, \
                         hk_270 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * hi_210[k]
                       + hk_266[k];

            t_211[k] = -ab_x[k] * hi_211[k]
                       + hk_267[k];

            t_212[k] = -ab_x[k] * hi_212[k]
                       + hk_268[k];

            t_213[k] = -ab_x[k] * hi_213[k]
                       + hk_269[k];

            t_214[k] = -ab_x[k] * hi_214[k]
                       + hk_270[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, hi_215, hi_216, hi_217, \
                         hi_218, hi_219, hk_271, hk_272, hk_273, hk_274, \
                         hk_275 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * hi_215[k]
                       + hk_271[k];

            t_216[k] = -ab_x[k] * hi_216[k]
                       + hk_272[k];

            t_217[k] = -ab_x[k] * hi_217[k]
                       + hk_273[k];

            t_218[k] = -ab_x[k] * hi_218[k]
                       + hk_274[k];

            t_219[k] = -ab_x[k] * hi_219[k]
                       + hk_275[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, hi_220, hi_221, hi_222, \
                         hi_223, hi_224, hk_276, hk_277, hk_278, hk_279, \
                         hk_288 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_x[k] * hi_220[k]
                       + hk_276[k];

            t_221[k] = -ab_x[k] * hi_221[k]
                       + hk_277[k];

            t_222[k] = -ab_x[k] * hi_222[k]
                       + hk_278[k];

            t_223[k] = -ab_x[k] * hi_223[k]
                       + hk_279[k];

            t_224[k] = -ab_x[k] * hi_224[k]
                       + hk_288[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, hi_225, hi_226, hi_227, \
                         hi_228, hi_229, hk_289, hk_290, hk_291, hk_292, \
                         hk_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_x[k] * hi_225[k]
                       + hk_289[k];

            t_226[k] = -ab_x[k] * hi_226[k]
                       + hk_290[k];

            t_227[k] = -ab_x[k] * hi_227[k]
                       + hk_291[k];

            t_228[k] = -ab_x[k] * hi_228[k]
                       + hk_292[k];

            t_229[k] = -ab_x[k] * hi_229[k]
                       + hk_293[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, hi_230, hi_231, hi_232, \
                         hi_233, hi_234, hk_294, hk_295, hk_296, hk_297, \
                         hk_298 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_x[k] * hi_230[k]
                       + hk_294[k];

            t_231[k] = -ab_x[k] * hi_231[k]
                       + hk_295[k];

            t_232[k] = -ab_x[k] * hi_232[k]
                       + hk_296[k];

            t_233[k] = -ab_x[k] * hi_233[k]
                       + hk_297[k];

            t_234[k] = -ab_x[k] * hi_234[k]
                       + hk_298[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, hi_235, hi_236, hi_237, \
                         hi_238, hi_239, hk_299, hk_300, hk_301, hk_302, \
                         hk_303 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_x[k] * hi_235[k]
                       + hk_299[k];

            t_236[k] = -ab_x[k] * hi_236[k]
                       + hk_300[k];

            t_237[k] = -ab_x[k] * hi_237[k]
                       + hk_301[k];

            t_238[k] = -ab_x[k] * hi_238[k]
                       + hk_302[k];

            t_239[k] = -ab_x[k] * hi_239[k]
                       + hk_303[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, hi_240, hi_241, hi_242, \
                         hi_243, hi_244, hk_304, hk_305, hk_306, hk_307, \
                         hk_308 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_x[k] * hi_240[k]
                       + hk_304[k];

            t_241[k] = -ab_x[k] * hi_241[k]
                       + hk_305[k];

            t_242[k] = -ab_x[k] * hi_242[k]
                       + hk_306[k];

            t_243[k] = -ab_x[k] * hi_243[k]
                       + hk_307[k];

            t_244[k] = -ab_x[k] * hi_244[k]
                       + hk_308[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, hi_245, hi_246, hi_247, \
                         hi_248, hi_249, hk_309, hk_310, hk_311, hk_312, \
                         hk_313 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_x[k] * hi_245[k]
                       + hk_309[k];

            t_246[k] = -ab_x[k] * hi_246[k]
                       + hk_310[k];

            t_247[k] = -ab_x[k] * hi_247[k]
                       + hk_311[k];

            t_248[k] = -ab_x[k] * hi_248[k]
                       + hk_312[k];

            t_249[k] = -ab_x[k] * hi_249[k]
                       + hk_313[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, hi_250, hi_251, hi_252, \
                         hi_253, hi_254, hk_314, hk_315, hk_324, hk_325, \
                         hk_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_x[k] * hi_250[k]
                       + hk_314[k];

            t_251[k] = -ab_x[k] * hi_251[k]
                       + hk_315[k];

            t_252[k] = -ab_x[k] * hi_252[k]
                       + hk_324[k];

            t_253[k] = -ab_x[k] * hi_253[k]
                       + hk_325[k];

            t_254[k] = -ab_x[k] * hi_254[k]
                       + hk_326[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, hi_255, hi_256, hi_257, \
                         hi_258, hi_259, hk_327, hk_328, hk_329, hk_330, \
                         hk_331 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_x[k] * hi_255[k]
                       + hk_327[k];

            t_256[k] = -ab_x[k] * hi_256[k]
                       + hk_328[k];

            t_257[k] = -ab_x[k] * hi_257[k]
                       + hk_329[k];

            t_258[k] = -ab_x[k] * hi_258[k]
                       + hk_330[k];

            t_259[k] = -ab_x[k] * hi_259[k]
                       + hk_331[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, hi_260, hi_261, hi_262, \
                         hi_263, hi_264, hk_332, hk_333, hk_334, hk_335, \
                         hk_336 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_x[k] * hi_260[k]
                       + hk_332[k];

            t_261[k] = -ab_x[k] * hi_261[k]
                       + hk_333[k];

            t_262[k] = -ab_x[k] * hi_262[k]
                       + hk_334[k];

            t_263[k] = -ab_x[k] * hi_263[k]
                       + hk_335[k];

            t_264[k] = -ab_x[k] * hi_264[k]
                       + hk_336[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, hi_265, hi_266, hi_267, \
                         hi_268, hi_269, hk_337, hk_338, hk_339, hk_340, \
                         hk_341 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_x[k] * hi_265[k]
                       + hk_337[k];

            t_266[k] = -ab_x[k] * hi_266[k]
                       + hk_338[k];

            t_267[k] = -ab_x[k] * hi_267[k]
                       + hk_339[k];

            t_268[k] = -ab_x[k] * hi_268[k]
                       + hk_340[k];

            t_269[k] = -ab_x[k] * hi_269[k]
                       + hk_341[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, hi_270, hi_271, hi_272, \
                         hi_273, hi_274, hk_342, hk_343, hk_344, hk_345, \
                         hk_346 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_x[k] * hi_270[k]
                       + hk_342[k];

            t_271[k] = -ab_x[k] * hi_271[k]
                       + hk_343[k];

            t_272[k] = -ab_x[k] * hi_272[k]
                       + hk_344[k];

            t_273[k] = -ab_x[k] * hi_273[k]
                       + hk_345[k];

            t_274[k] = -ab_x[k] * hi_274[k]
                       + hk_346[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, hi_275, hi_276, hi_277, \
                         hi_278, hi_279, hk_347, hk_348, hk_349, hk_350, \
                         hk_351 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_x[k] * hi_275[k]
                       + hk_347[k];

            t_276[k] = -ab_x[k] * hi_276[k]
                       + hk_348[k];

            t_277[k] = -ab_x[k] * hi_277[k]
                       + hk_349[k];

            t_278[k] = -ab_x[k] * hi_278[k]
                       + hk_350[k];

            t_279[k] = -ab_x[k] * hi_279[k]
                       + hk_351[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, hi_280, hi_281, hi_282, \
                         hi_283, hi_284, hk_360, hk_361, hk_362, hk_363, \
                         hk_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_x[k] * hi_280[k]
                       + hk_360[k];

            t_281[k] = -ab_x[k] * hi_281[k]
                       + hk_361[k];

            t_282[k] = -ab_x[k] * hi_282[k]
                       + hk_362[k];

            t_283[k] = -ab_x[k] * hi_283[k]
                       + hk_363[k];

            t_284[k] = -ab_x[k] * hi_284[k]
                       + hk_364[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, hi_285, hi_286, hi_287, \
                         hi_288, hi_289, hk_365, hk_366, hk_367, hk_368, \
                         hk_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_x[k] * hi_285[k]
                       + hk_365[k];

            t_286[k] = -ab_x[k] * hi_286[k]
                       + hk_366[k];

            t_287[k] = -ab_x[k] * hi_287[k]
                       + hk_367[k];

            t_288[k] = -ab_x[k] * hi_288[k]
                       + hk_368[k];

            t_289[k] = -ab_x[k] * hi_289[k]
                       + hk_369[k];
        }
    }
}

static auto
compute_hrr_ii_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hi, const size_t hk, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_290 = buffer.data(target + 290 * ncomps + c);
        auto *t_291 = buffer.data(target + 291 * ncomps + c);
        auto *t_292 = buffer.data(target + 292 * ncomps + c);
        auto *t_293 = buffer.data(target + 293 * ncomps + c);
        auto *t_294 = buffer.data(target + 294 * ncomps + c);
        auto *t_295 = buffer.data(target + 295 * ncomps + c);
        auto *t_296 = buffer.data(target + 296 * ncomps + c);
        auto *t_297 = buffer.data(target + 297 * ncomps + c);
        auto *t_298 = buffer.data(target + 298 * ncomps + c);
        auto *t_299 = buffer.data(target + 299 * ncomps + c);
        auto *t_300 = buffer.data(target + 300 * ncomps + c);
        auto *t_301 = buffer.data(target + 301 * ncomps + c);
        auto *t_302 = buffer.data(target + 302 * ncomps + c);
        auto *t_303 = buffer.data(target + 303 * ncomps + c);
        auto *t_304 = buffer.data(target + 304 * ncomps + c);
        auto *t_305 = buffer.data(target + 305 * ncomps + c);
        auto *t_306 = buffer.data(target + 306 * ncomps + c);
        auto *t_307 = buffer.data(target + 307 * ncomps + c);
        auto *t_308 = buffer.data(target + 308 * ncomps + c);
        auto *t_309 = buffer.data(target + 309 * ncomps + c);
        auto *t_310 = buffer.data(target + 310 * ncomps + c);
        auto *t_311 = buffer.data(target + 311 * ncomps + c);
        auto *t_312 = buffer.data(target + 312 * ncomps + c);
        auto *t_313 = buffer.data(target + 313 * ncomps + c);
        auto *t_314 = buffer.data(target + 314 * ncomps + c);
        auto *t_315 = buffer.data(target + 315 * ncomps + c);
        auto *t_316 = buffer.data(target + 316 * ncomps + c);
        auto *t_317 = buffer.data(target + 317 * ncomps + c);
        auto *t_318 = buffer.data(target + 318 * ncomps + c);
        auto *t_319 = buffer.data(target + 319 * ncomps + c);
        auto *t_320 = buffer.data(target + 320 * ncomps + c);
        auto *t_321 = buffer.data(target + 321 * ncomps + c);
        auto *t_322 = buffer.data(target + 322 * ncomps + c);
        auto *t_323 = buffer.data(target + 323 * ncomps + c);
        auto *t_324 = buffer.data(target + 324 * ncomps + c);
        auto *t_325 = buffer.data(target + 325 * ncomps + c);
        auto *t_326 = buffer.data(target + 326 * ncomps + c);
        auto *t_327 = buffer.data(target + 327 * ncomps + c);
        auto *t_328 = buffer.data(target + 328 * ncomps + c);
        auto *t_329 = buffer.data(target + 329 * ncomps + c);
        auto *t_330 = buffer.data(target + 330 * ncomps + c);
        auto *t_331 = buffer.data(target + 331 * ncomps + c);
        auto *t_332 = buffer.data(target + 332 * ncomps + c);
        auto *t_333 = buffer.data(target + 333 * ncomps + c);
        auto *t_334 = buffer.data(target + 334 * ncomps + c);
        auto *t_335 = buffer.data(target + 335 * ncomps + c);
        auto *t_336 = buffer.data(target + 336 * ncomps + c);
        auto *t_337 = buffer.data(target + 337 * ncomps + c);
        auto *t_338 = buffer.data(target + 338 * ncomps + c);
        auto *t_339 = buffer.data(target + 339 * ncomps + c);
        auto *t_340 = buffer.data(target + 340 * ncomps + c);
        auto *t_341 = buffer.data(target + 341 * ncomps + c);
        auto *t_342 = buffer.data(target + 342 * ncomps + c);
        auto *t_343 = buffer.data(target + 343 * ncomps + c);
        auto *t_344 = buffer.data(target + 344 * ncomps + c);
        auto *t_345 = buffer.data(target + 345 * ncomps + c);
        auto *t_346 = buffer.data(target + 346 * ncomps + c);
        auto *t_347 = buffer.data(target + 347 * ncomps + c);
        auto *t_348 = buffer.data(target + 348 * ncomps + c);
        auto *t_349 = buffer.data(target + 349 * ncomps + c);
        auto *t_350 = buffer.data(target + 350 * ncomps + c);
        auto *t_351 = buffer.data(target + 351 * ncomps + c);
        auto *t_352 = buffer.data(target + 352 * ncomps + c);
        auto *t_353 = buffer.data(target + 353 * ncomps + c);
        auto *t_354 = buffer.data(target + 354 * ncomps + c);
        auto *t_355 = buffer.data(target + 355 * ncomps + c);
        auto *t_356 = buffer.data(target + 356 * ncomps + c);
        auto *t_357 = buffer.data(target + 357 * ncomps + c);
        auto *t_358 = buffer.data(target + 358 * ncomps + c);
        auto *t_359 = buffer.data(target + 359 * ncomps + c);
        auto *t_360 = buffer.data(target + 360 * ncomps + c);
        auto *t_361 = buffer.data(target + 361 * ncomps + c);
        auto *t_362 = buffer.data(target + 362 * ncomps + c);
        auto *t_363 = buffer.data(target + 363 * ncomps + c);
        auto *t_364 = buffer.data(target + 364 * ncomps + c);
        auto *t_365 = buffer.data(target + 365 * ncomps + c);
        auto *t_366 = buffer.data(target + 366 * ncomps + c);
        auto *t_367 = buffer.data(target + 367 * ncomps + c);
        auto *t_368 = buffer.data(target + 368 * ncomps + c);
        auto *t_369 = buffer.data(target + 369 * ncomps + c);
        auto *t_370 = buffer.data(target + 370 * ncomps + c);
        auto *t_371 = buffer.data(target + 371 * ncomps + c);
        auto *t_372 = buffer.data(target + 372 * ncomps + c);
        auto *t_373 = buffer.data(target + 373 * ncomps + c);
        auto *t_374 = buffer.data(target + 374 * ncomps + c);
        auto *t_375 = buffer.data(target + 375 * ncomps + c);
        auto *t_376 = buffer.data(target + 376 * ncomps + c);
        auto *t_377 = buffer.data(target + 377 * ncomps + c);
        auto *t_378 = buffer.data(target + 378 * ncomps + c);
        auto *t_379 = buffer.data(target + 379 * ncomps + c);
        auto *t_380 = buffer.data(target + 380 * ncomps + c);
        auto *t_381 = buffer.data(target + 381 * ncomps + c);
        auto *t_382 = buffer.data(target + 382 * ncomps + c);
        auto *t_383 = buffer.data(target + 383 * ncomps + c);
        auto *t_384 = buffer.data(target + 384 * ncomps + c);
        auto *t_385 = buffer.data(target + 385 * ncomps + c);
        auto *t_386 = buffer.data(target + 386 * ncomps + c);
        auto *t_387 = buffer.data(target + 387 * ncomps + c);
        auto *t_388 = buffer.data(target + 388 * ncomps + c);
        auto *t_389 = buffer.data(target + 389 * ncomps + c);
        auto *t_390 = buffer.data(target + 390 * ncomps + c);
        auto *t_391 = buffer.data(target + 391 * ncomps + c);
        auto *t_392 = buffer.data(target + 392 * ncomps + c);
        auto *t_393 = buffer.data(target + 393 * ncomps + c);
        auto *t_394 = buffer.data(target + 394 * ncomps + c);
        auto *t_395 = buffer.data(target + 395 * ncomps + c);
        auto *t_396 = buffer.data(target + 396 * ncomps + c);
        auto *t_397 = buffer.data(target + 397 * ncomps + c);
        auto *t_398 = buffer.data(target + 398 * ncomps + c);
        auto *t_399 = buffer.data(target + 399 * ncomps + c);
        auto *t_400 = buffer.data(target + 400 * ncomps + c);
        auto *t_401 = buffer.data(target + 401 * ncomps + c);
        auto *t_402 = buffer.data(target + 402 * ncomps + c);
        auto *t_403 = buffer.data(target + 403 * ncomps + c);
        auto *t_404 = buffer.data(target + 404 * ncomps + c);
        auto *t_405 = buffer.data(target + 405 * ncomps + c);
        auto *t_406 = buffer.data(target + 406 * ncomps + c);
        auto *t_407 = buffer.data(target + 407 * ncomps + c);
        auto *t_408 = buffer.data(target + 408 * ncomps + c);
        auto *t_409 = buffer.data(target + 409 * ncomps + c);
        auto *t_410 = buffer.data(target + 410 * ncomps + c);
        auto *t_411 = buffer.data(target + 411 * ncomps + c);
        auto *t_412 = buffer.data(target + 412 * ncomps + c);
        auto *t_413 = buffer.data(target + 413 * ncomps + c);
        auto *t_414 = buffer.data(target + 414 * ncomps + c);
        auto *t_415 = buffer.data(target + 415 * ncomps + c);
        auto *t_416 = buffer.data(target + 416 * ncomps + c);
        auto *t_417 = buffer.data(target + 417 * ncomps + c);
        auto *t_418 = buffer.data(target + 418 * ncomps + c);
        auto *t_419 = buffer.data(target + 419 * ncomps + c);
        auto *t_420 = buffer.data(target + 420 * ncomps + c);
        auto *t_421 = buffer.data(target + 421 * ncomps + c);
        auto *t_422 = buffer.data(target + 422 * ncomps + c);
        auto *t_423 = buffer.data(target + 423 * ncomps + c);
        auto *t_424 = buffer.data(target + 424 * ncomps + c);
        auto *t_425 = buffer.data(target + 425 * ncomps + c);
        auto *t_426 = buffer.data(target + 426 * ncomps + c);
        auto *t_427 = buffer.data(target + 427 * ncomps + c);
        auto *t_428 = buffer.data(target + 428 * ncomps + c);
        auto *t_429 = buffer.data(target + 429 * ncomps + c);
        auto *t_430 = buffer.data(target + 430 * ncomps + c);
        auto *t_431 = buffer.data(target + 431 * ncomps + c);
        auto *t_432 = buffer.data(target + 432 * ncomps + c);
        auto *t_433 = buffer.data(target + 433 * ncomps + c);
        auto *t_434 = buffer.data(target + 434 * ncomps + c);

        const auto *ab_x = coordinates.data(6);

        const auto *hi_290 = buffer.data(hi + 290 * ncomps + c);
        const auto *hi_291 = buffer.data(hi + 291 * ncomps + c);
        const auto *hi_292 = buffer.data(hi + 292 * ncomps + c);
        const auto *hi_293 = buffer.data(hi + 293 * ncomps + c);
        const auto *hi_294 = buffer.data(hi + 294 * ncomps + c);
        const auto *hi_295 = buffer.data(hi + 295 * ncomps + c);
        const auto *hi_296 = buffer.data(hi + 296 * ncomps + c);
        const auto *hi_297 = buffer.data(hi + 297 * ncomps + c);
        const auto *hi_298 = buffer.data(hi + 298 * ncomps + c);
        const auto *hi_299 = buffer.data(hi + 299 * ncomps + c);
        const auto *hi_300 = buffer.data(hi + 300 * ncomps + c);
        const auto *hi_301 = buffer.data(hi + 301 * ncomps + c);
        const auto *hi_302 = buffer.data(hi + 302 * ncomps + c);
        const auto *hi_303 = buffer.data(hi + 303 * ncomps + c);
        const auto *hi_304 = buffer.data(hi + 304 * ncomps + c);
        const auto *hi_305 = buffer.data(hi + 305 * ncomps + c);
        const auto *hi_306 = buffer.data(hi + 306 * ncomps + c);
        const auto *hi_307 = buffer.data(hi + 307 * ncomps + c);
        const auto *hi_308 = buffer.data(hi + 308 * ncomps + c);
        const auto *hi_309 = buffer.data(hi + 309 * ncomps + c);
        const auto *hi_310 = buffer.data(hi + 310 * ncomps + c);
        const auto *hi_311 = buffer.data(hi + 311 * ncomps + c);
        const auto *hi_312 = buffer.data(hi + 312 * ncomps + c);
        const auto *hi_313 = buffer.data(hi + 313 * ncomps + c);
        const auto *hi_314 = buffer.data(hi + 314 * ncomps + c);
        const auto *hi_315 = buffer.data(hi + 315 * ncomps + c);
        const auto *hi_316 = buffer.data(hi + 316 * ncomps + c);
        const auto *hi_317 = buffer.data(hi + 317 * ncomps + c);
        const auto *hi_318 = buffer.data(hi + 318 * ncomps + c);
        const auto *hi_319 = buffer.data(hi + 319 * ncomps + c);
        const auto *hi_320 = buffer.data(hi + 320 * ncomps + c);
        const auto *hi_321 = buffer.data(hi + 321 * ncomps + c);
        const auto *hi_322 = buffer.data(hi + 322 * ncomps + c);
        const auto *hi_323 = buffer.data(hi + 323 * ncomps + c);
        const auto *hi_324 = buffer.data(hi + 324 * ncomps + c);
        const auto *hi_325 = buffer.data(hi + 325 * ncomps + c);
        const auto *hi_326 = buffer.data(hi + 326 * ncomps + c);
        const auto *hi_327 = buffer.data(hi + 327 * ncomps + c);
        const auto *hi_328 = buffer.data(hi + 328 * ncomps + c);
        const auto *hi_329 = buffer.data(hi + 329 * ncomps + c);
        const auto *hi_330 = buffer.data(hi + 330 * ncomps + c);
        const auto *hi_331 = buffer.data(hi + 331 * ncomps + c);
        const auto *hi_332 = buffer.data(hi + 332 * ncomps + c);
        const auto *hi_333 = buffer.data(hi + 333 * ncomps + c);
        const auto *hi_334 = buffer.data(hi + 334 * ncomps + c);
        const auto *hi_335 = buffer.data(hi + 335 * ncomps + c);
        const auto *hi_336 = buffer.data(hi + 336 * ncomps + c);
        const auto *hi_337 = buffer.data(hi + 337 * ncomps + c);
        const auto *hi_338 = buffer.data(hi + 338 * ncomps + c);
        const auto *hi_339 = buffer.data(hi + 339 * ncomps + c);
        const auto *hi_340 = buffer.data(hi + 340 * ncomps + c);
        const auto *hi_341 = buffer.data(hi + 341 * ncomps + c);
        const auto *hi_342 = buffer.data(hi + 342 * ncomps + c);
        const auto *hi_343 = buffer.data(hi + 343 * ncomps + c);
        const auto *hi_344 = buffer.data(hi + 344 * ncomps + c);
        const auto *hi_345 = buffer.data(hi + 345 * ncomps + c);
        const auto *hi_346 = buffer.data(hi + 346 * ncomps + c);
        const auto *hi_347 = buffer.data(hi + 347 * ncomps + c);
        const auto *hi_348 = buffer.data(hi + 348 * ncomps + c);
        const auto *hi_349 = buffer.data(hi + 349 * ncomps + c);
        const auto *hi_350 = buffer.data(hi + 350 * ncomps + c);
        const auto *hi_351 = buffer.data(hi + 351 * ncomps + c);
        const auto *hi_352 = buffer.data(hi + 352 * ncomps + c);
        const auto *hi_353 = buffer.data(hi + 353 * ncomps + c);
        const auto *hi_354 = buffer.data(hi + 354 * ncomps + c);
        const auto *hi_355 = buffer.data(hi + 355 * ncomps + c);
        const auto *hi_356 = buffer.data(hi + 356 * ncomps + c);
        const auto *hi_357 = buffer.data(hi + 357 * ncomps + c);
        const auto *hi_358 = buffer.data(hi + 358 * ncomps + c);
        const auto *hi_359 = buffer.data(hi + 359 * ncomps + c);
        const auto *hi_360 = buffer.data(hi + 360 * ncomps + c);
        const auto *hi_361 = buffer.data(hi + 361 * ncomps + c);
        const auto *hi_362 = buffer.data(hi + 362 * ncomps + c);
        const auto *hi_363 = buffer.data(hi + 363 * ncomps + c);
        const auto *hi_364 = buffer.data(hi + 364 * ncomps + c);
        const auto *hi_365 = buffer.data(hi + 365 * ncomps + c);
        const auto *hi_366 = buffer.data(hi + 366 * ncomps + c);
        const auto *hi_367 = buffer.data(hi + 367 * ncomps + c);
        const auto *hi_368 = buffer.data(hi + 368 * ncomps + c);
        const auto *hi_369 = buffer.data(hi + 369 * ncomps + c);
        const auto *hi_370 = buffer.data(hi + 370 * ncomps + c);
        const auto *hi_371 = buffer.data(hi + 371 * ncomps + c);
        const auto *hi_372 = buffer.data(hi + 372 * ncomps + c);
        const auto *hi_373 = buffer.data(hi + 373 * ncomps + c);
        const auto *hi_374 = buffer.data(hi + 374 * ncomps + c);
        const auto *hi_375 = buffer.data(hi + 375 * ncomps + c);
        const auto *hi_376 = buffer.data(hi + 376 * ncomps + c);
        const auto *hi_377 = buffer.data(hi + 377 * ncomps + c);
        const auto *hi_378 = buffer.data(hi + 378 * ncomps + c);
        const auto *hi_379 = buffer.data(hi + 379 * ncomps + c);
        const auto *hi_380 = buffer.data(hi + 380 * ncomps + c);
        const auto *hi_381 = buffer.data(hi + 381 * ncomps + c);
        const auto *hi_382 = buffer.data(hi + 382 * ncomps + c);
        const auto *hi_383 = buffer.data(hi + 383 * ncomps + c);
        const auto *hi_384 = buffer.data(hi + 384 * ncomps + c);
        const auto *hi_385 = buffer.data(hi + 385 * ncomps + c);
        const auto *hi_386 = buffer.data(hi + 386 * ncomps + c);
        const auto *hi_387 = buffer.data(hi + 387 * ncomps + c);
        const auto *hi_388 = buffer.data(hi + 388 * ncomps + c);
        const auto *hi_389 = buffer.data(hi + 389 * ncomps + c);
        const auto *hi_390 = buffer.data(hi + 390 * ncomps + c);
        const auto *hi_391 = buffer.data(hi + 391 * ncomps + c);
        const auto *hi_392 = buffer.data(hi + 392 * ncomps + c);
        const auto *hi_393 = buffer.data(hi + 393 * ncomps + c);
        const auto *hi_394 = buffer.data(hi + 394 * ncomps + c);
        const auto *hi_395 = buffer.data(hi + 395 * ncomps + c);
        const auto *hi_396 = buffer.data(hi + 396 * ncomps + c);
        const auto *hi_397 = buffer.data(hi + 397 * ncomps + c);
        const auto *hi_398 = buffer.data(hi + 398 * ncomps + c);
        const auto *hi_399 = buffer.data(hi + 399 * ncomps + c);
        const auto *hi_400 = buffer.data(hi + 400 * ncomps + c);
        const auto *hi_401 = buffer.data(hi + 401 * ncomps + c);
        const auto *hi_402 = buffer.data(hi + 402 * ncomps + c);
        const auto *hi_403 = buffer.data(hi + 403 * ncomps + c);
        const auto *hi_404 = buffer.data(hi + 404 * ncomps + c);
        const auto *hi_405 = buffer.data(hi + 405 * ncomps + c);
        const auto *hi_406 = buffer.data(hi + 406 * ncomps + c);
        const auto *hi_407 = buffer.data(hi + 407 * ncomps + c);
        const auto *hi_408 = buffer.data(hi + 408 * ncomps + c);
        const auto *hi_409 = buffer.data(hi + 409 * ncomps + c);
        const auto *hi_410 = buffer.data(hi + 410 * ncomps + c);
        const auto *hi_411 = buffer.data(hi + 411 * ncomps + c);
        const auto *hi_412 = buffer.data(hi + 412 * ncomps + c);
        const auto *hi_413 = buffer.data(hi + 413 * ncomps + c);
        const auto *hi_414 = buffer.data(hi + 414 * ncomps + c);
        const auto *hi_415 = buffer.data(hi + 415 * ncomps + c);
        const auto *hi_416 = buffer.data(hi + 416 * ncomps + c);
        const auto *hi_417 = buffer.data(hi + 417 * ncomps + c);
        const auto *hi_418 = buffer.data(hi + 418 * ncomps + c);
        const auto *hi_419 = buffer.data(hi + 419 * ncomps + c);
        const auto *hi_420 = buffer.data(hi + 420 * ncomps + c);
        const auto *hi_421 = buffer.data(hi + 421 * ncomps + c);
        const auto *hi_422 = buffer.data(hi + 422 * ncomps + c);
        const auto *hi_423 = buffer.data(hi + 423 * ncomps + c);
        const auto *hi_424 = buffer.data(hi + 424 * ncomps + c);
        const auto *hi_425 = buffer.data(hi + 425 * ncomps + c);
        const auto *hi_426 = buffer.data(hi + 426 * ncomps + c);
        const auto *hi_427 = buffer.data(hi + 427 * ncomps + c);
        const auto *hi_428 = buffer.data(hi + 428 * ncomps + c);
        const auto *hi_429 = buffer.data(hi + 429 * ncomps + c);
        const auto *hi_430 = buffer.data(hi + 430 * ncomps + c);
        const auto *hi_431 = buffer.data(hi + 431 * ncomps + c);
        const auto *hi_432 = buffer.data(hi + 432 * ncomps + c);
        const auto *hi_433 = buffer.data(hi + 433 * ncomps + c);
        const auto *hi_434 = buffer.data(hi + 434 * ncomps + c);

        const auto *hk_370 = buffer.data(hk + 370 * ncomps + c);
        const auto *hk_371 = buffer.data(hk + 371 * ncomps + c);
        const auto *hk_372 = buffer.data(hk + 372 * ncomps + c);
        const auto *hk_373 = buffer.data(hk + 373 * ncomps + c);
        const auto *hk_374 = buffer.data(hk + 374 * ncomps + c);
        const auto *hk_375 = buffer.data(hk + 375 * ncomps + c);
        const auto *hk_376 = buffer.data(hk + 376 * ncomps + c);
        const auto *hk_377 = buffer.data(hk + 377 * ncomps + c);
        const auto *hk_378 = buffer.data(hk + 378 * ncomps + c);
        const auto *hk_379 = buffer.data(hk + 379 * ncomps + c);
        const auto *hk_380 = buffer.data(hk + 380 * ncomps + c);
        const auto *hk_381 = buffer.data(hk + 381 * ncomps + c);
        const auto *hk_382 = buffer.data(hk + 382 * ncomps + c);
        const auto *hk_383 = buffer.data(hk + 383 * ncomps + c);
        const auto *hk_384 = buffer.data(hk + 384 * ncomps + c);
        const auto *hk_385 = buffer.data(hk + 385 * ncomps + c);
        const auto *hk_386 = buffer.data(hk + 386 * ncomps + c);
        const auto *hk_387 = buffer.data(hk + 387 * ncomps + c);
        const auto *hk_396 = buffer.data(hk + 396 * ncomps + c);
        const auto *hk_397 = buffer.data(hk + 397 * ncomps + c);
        const auto *hk_398 = buffer.data(hk + 398 * ncomps + c);
        const auto *hk_399 = buffer.data(hk + 399 * ncomps + c);
        const auto *hk_400 = buffer.data(hk + 400 * ncomps + c);
        const auto *hk_401 = buffer.data(hk + 401 * ncomps + c);
        const auto *hk_402 = buffer.data(hk + 402 * ncomps + c);
        const auto *hk_403 = buffer.data(hk + 403 * ncomps + c);
        const auto *hk_404 = buffer.data(hk + 404 * ncomps + c);
        const auto *hk_405 = buffer.data(hk + 405 * ncomps + c);
        const auto *hk_406 = buffer.data(hk + 406 * ncomps + c);
        const auto *hk_407 = buffer.data(hk + 407 * ncomps + c);
        const auto *hk_408 = buffer.data(hk + 408 * ncomps + c);
        const auto *hk_409 = buffer.data(hk + 409 * ncomps + c);
        const auto *hk_410 = buffer.data(hk + 410 * ncomps + c);
        const auto *hk_411 = buffer.data(hk + 411 * ncomps + c);
        const auto *hk_412 = buffer.data(hk + 412 * ncomps + c);
        const auto *hk_413 = buffer.data(hk + 413 * ncomps + c);
        const auto *hk_414 = buffer.data(hk + 414 * ncomps + c);
        const auto *hk_415 = buffer.data(hk + 415 * ncomps + c);
        const auto *hk_416 = buffer.data(hk + 416 * ncomps + c);
        const auto *hk_417 = buffer.data(hk + 417 * ncomps + c);
        const auto *hk_418 = buffer.data(hk + 418 * ncomps + c);
        const auto *hk_419 = buffer.data(hk + 419 * ncomps + c);
        const auto *hk_420 = buffer.data(hk + 420 * ncomps + c);
        const auto *hk_421 = buffer.data(hk + 421 * ncomps + c);
        const auto *hk_422 = buffer.data(hk + 422 * ncomps + c);
        const auto *hk_423 = buffer.data(hk + 423 * ncomps + c);
        const auto *hk_432 = buffer.data(hk + 432 * ncomps + c);
        const auto *hk_433 = buffer.data(hk + 433 * ncomps + c);
        const auto *hk_434 = buffer.data(hk + 434 * ncomps + c);
        const auto *hk_435 = buffer.data(hk + 435 * ncomps + c);
        const auto *hk_436 = buffer.data(hk + 436 * ncomps + c);
        const auto *hk_437 = buffer.data(hk + 437 * ncomps + c);
        const auto *hk_438 = buffer.data(hk + 438 * ncomps + c);
        const auto *hk_439 = buffer.data(hk + 439 * ncomps + c);
        const auto *hk_440 = buffer.data(hk + 440 * ncomps + c);
        const auto *hk_441 = buffer.data(hk + 441 * ncomps + c);
        const auto *hk_442 = buffer.data(hk + 442 * ncomps + c);
        const auto *hk_443 = buffer.data(hk + 443 * ncomps + c);
        const auto *hk_444 = buffer.data(hk + 444 * ncomps + c);
        const auto *hk_445 = buffer.data(hk + 445 * ncomps + c);
        const auto *hk_446 = buffer.data(hk + 446 * ncomps + c);
        const auto *hk_447 = buffer.data(hk + 447 * ncomps + c);
        const auto *hk_448 = buffer.data(hk + 448 * ncomps + c);
        const auto *hk_449 = buffer.data(hk + 449 * ncomps + c);
        const auto *hk_450 = buffer.data(hk + 450 * ncomps + c);
        const auto *hk_451 = buffer.data(hk + 451 * ncomps + c);
        const auto *hk_452 = buffer.data(hk + 452 * ncomps + c);
        const auto *hk_453 = buffer.data(hk + 453 * ncomps + c);
        const auto *hk_454 = buffer.data(hk + 454 * ncomps + c);
        const auto *hk_455 = buffer.data(hk + 455 * ncomps + c);
        const auto *hk_456 = buffer.data(hk + 456 * ncomps + c);
        const auto *hk_457 = buffer.data(hk + 457 * ncomps + c);
        const auto *hk_458 = buffer.data(hk + 458 * ncomps + c);
        const auto *hk_459 = buffer.data(hk + 459 * ncomps + c);
        const auto *hk_468 = buffer.data(hk + 468 * ncomps + c);
        const auto *hk_469 = buffer.data(hk + 469 * ncomps + c);
        const auto *hk_470 = buffer.data(hk + 470 * ncomps + c);
        const auto *hk_471 = buffer.data(hk + 471 * ncomps + c);
        const auto *hk_472 = buffer.data(hk + 472 * ncomps + c);
        const auto *hk_473 = buffer.data(hk + 473 * ncomps + c);
        const auto *hk_474 = buffer.data(hk + 474 * ncomps + c);
        const auto *hk_475 = buffer.data(hk + 475 * ncomps + c);
        const auto *hk_476 = buffer.data(hk + 476 * ncomps + c);
        const auto *hk_477 = buffer.data(hk + 477 * ncomps + c);
        const auto *hk_478 = buffer.data(hk + 478 * ncomps + c);
        const auto *hk_479 = buffer.data(hk + 479 * ncomps + c);
        const auto *hk_480 = buffer.data(hk + 480 * ncomps + c);
        const auto *hk_481 = buffer.data(hk + 481 * ncomps + c);
        const auto *hk_482 = buffer.data(hk + 482 * ncomps + c);
        const auto *hk_483 = buffer.data(hk + 483 * ncomps + c);
        const auto *hk_484 = buffer.data(hk + 484 * ncomps + c);
        const auto *hk_485 = buffer.data(hk + 485 * ncomps + c);
        const auto *hk_486 = buffer.data(hk + 486 * ncomps + c);
        const auto *hk_487 = buffer.data(hk + 487 * ncomps + c);
        const auto *hk_488 = buffer.data(hk + 488 * ncomps + c);
        const auto *hk_489 = buffer.data(hk + 489 * ncomps + c);
        const auto *hk_490 = buffer.data(hk + 490 * ncomps + c);
        const auto *hk_491 = buffer.data(hk + 491 * ncomps + c);
        const auto *hk_492 = buffer.data(hk + 492 * ncomps + c);
        const auto *hk_493 = buffer.data(hk + 493 * ncomps + c);
        const auto *hk_494 = buffer.data(hk + 494 * ncomps + c);
        const auto *hk_495 = buffer.data(hk + 495 * ncomps + c);
        const auto *hk_504 = buffer.data(hk + 504 * ncomps + c);
        const auto *hk_505 = buffer.data(hk + 505 * ncomps + c);
        const auto *hk_506 = buffer.data(hk + 506 * ncomps + c);
        const auto *hk_507 = buffer.data(hk + 507 * ncomps + c);
        const auto *hk_508 = buffer.data(hk + 508 * ncomps + c);
        const auto *hk_509 = buffer.data(hk + 509 * ncomps + c);
        const auto *hk_510 = buffer.data(hk + 510 * ncomps + c);
        const auto *hk_511 = buffer.data(hk + 511 * ncomps + c);
        const auto *hk_512 = buffer.data(hk + 512 * ncomps + c);
        const auto *hk_513 = buffer.data(hk + 513 * ncomps + c);
        const auto *hk_514 = buffer.data(hk + 514 * ncomps + c);
        const auto *hk_515 = buffer.data(hk + 515 * ncomps + c);
        const auto *hk_516 = buffer.data(hk + 516 * ncomps + c);
        const auto *hk_517 = buffer.data(hk + 517 * ncomps + c);
        const auto *hk_518 = buffer.data(hk + 518 * ncomps + c);
        const auto *hk_519 = buffer.data(hk + 519 * ncomps + c);
        const auto *hk_520 = buffer.data(hk + 520 * ncomps + c);
        const auto *hk_521 = buffer.data(hk + 521 * ncomps + c);
        const auto *hk_522 = buffer.data(hk + 522 * ncomps + c);
        const auto *hk_523 = buffer.data(hk + 523 * ncomps + c);
        const auto *hk_524 = buffer.data(hk + 524 * ncomps + c);
        const auto *hk_525 = buffer.data(hk + 525 * ncomps + c);
        const auto *hk_526 = buffer.data(hk + 526 * ncomps + c);
        const auto *hk_527 = buffer.data(hk + 527 * ncomps + c);
        const auto *hk_528 = buffer.data(hk + 528 * ncomps + c);
        const auto *hk_529 = buffer.data(hk + 529 * ncomps + c);
        const auto *hk_530 = buffer.data(hk + 530 * ncomps + c);
        const auto *hk_531 = buffer.data(hk + 531 * ncomps + c);
        const auto *hk_540 = buffer.data(hk + 540 * ncomps + c);
        const auto *hk_541 = buffer.data(hk + 541 * ncomps + c);
        const auto *hk_542 = buffer.data(hk + 542 * ncomps + c);
        const auto *hk_543 = buffer.data(hk + 543 * ncomps + c);
        const auto *hk_544 = buffer.data(hk + 544 * ncomps + c);
        const auto *hk_545 = buffer.data(hk + 545 * ncomps + c);
        const auto *hk_546 = buffer.data(hk + 546 * ncomps + c);
        const auto *hk_547 = buffer.data(hk + 547 * ncomps + c);
        const auto *hk_548 = buffer.data(hk + 548 * ncomps + c);
        const auto *hk_549 = buffer.data(hk + 549 * ncomps + c);
        const auto *hk_550 = buffer.data(hk + 550 * ncomps + c);
        const auto *hk_551 = buffer.data(hk + 551 * ncomps + c);
        const auto *hk_552 = buffer.data(hk + 552 * ncomps + c);
        const auto *hk_553 = buffer.data(hk + 553 * ncomps + c);
        const auto *hk_554 = buffer.data(hk + 554 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, hi_290, hi_291, hi_292, \
                         hi_293, hi_294, hk_370, hk_371, hk_372, hk_373, \
                         hk_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_x[k] * hi_290[k]
                       + hk_370[k];

            t_291[k] = -ab_x[k] * hi_291[k]
                       + hk_371[k];

            t_292[k] = -ab_x[k] * hi_292[k]
                       + hk_372[k];

            t_293[k] = -ab_x[k] * hi_293[k]
                       + hk_373[k];

            t_294[k] = -ab_x[k] * hi_294[k]
                       + hk_374[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, hi_295, hi_296, hi_297, \
                         hi_298, hi_299, hk_375, hk_376, hk_377, hk_378, \
                         hk_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_x[k] * hi_295[k]
                       + hk_375[k];

            t_296[k] = -ab_x[k] * hi_296[k]
                       + hk_376[k];

            t_297[k] = -ab_x[k] * hi_297[k]
                       + hk_377[k];

            t_298[k] = -ab_x[k] * hi_298[k]
                       + hk_378[k];

            t_299[k] = -ab_x[k] * hi_299[k]
                       + hk_379[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, hi_300, hi_301, hi_302, \
                         hi_303, hi_304, hk_380, hk_381, hk_382, hk_383, \
                         hk_384 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_x[k] * hi_300[k]
                       + hk_380[k];

            t_301[k] = -ab_x[k] * hi_301[k]
                       + hk_381[k];

            t_302[k] = -ab_x[k] * hi_302[k]
                       + hk_382[k];

            t_303[k] = -ab_x[k] * hi_303[k]
                       + hk_383[k];

            t_304[k] = -ab_x[k] * hi_304[k]
                       + hk_384[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, hi_305, hi_306, hi_307, \
                         hi_308, hi_309, hk_385, hk_386, hk_387, hk_396, \
                         hk_397 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_x[k] * hi_305[k]
                       + hk_385[k];

            t_306[k] = -ab_x[k] * hi_306[k]
                       + hk_386[k];

            t_307[k] = -ab_x[k] * hi_307[k]
                       + hk_387[k];

            t_308[k] = -ab_x[k] * hi_308[k]
                       + hk_396[k];

            t_309[k] = -ab_x[k] * hi_309[k]
                       + hk_397[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, hi_310, hi_311, hi_312, \
                         hi_313, hi_314, hk_398, hk_399, hk_400, hk_401, \
                         hk_402 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_x[k] * hi_310[k]
                       + hk_398[k];

            t_311[k] = -ab_x[k] * hi_311[k]
                       + hk_399[k];

            t_312[k] = -ab_x[k] * hi_312[k]
                       + hk_400[k];

            t_313[k] = -ab_x[k] * hi_313[k]
                       + hk_401[k];

            t_314[k] = -ab_x[k] * hi_314[k]
                       + hk_402[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, hi_315, hi_316, hi_317, \
                         hi_318, hi_319, hk_403, hk_404, hk_405, hk_406, \
                         hk_407 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_x[k] * hi_315[k]
                       + hk_403[k];

            t_316[k] = -ab_x[k] * hi_316[k]
                       + hk_404[k];

            t_317[k] = -ab_x[k] * hi_317[k]
                       + hk_405[k];

            t_318[k] = -ab_x[k] * hi_318[k]
                       + hk_406[k];

            t_319[k] = -ab_x[k] * hi_319[k]
                       + hk_407[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, hi_320, hi_321, hi_322, \
                         hi_323, hi_324, hk_408, hk_409, hk_410, hk_411, \
                         hk_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_x[k] * hi_320[k]
                       + hk_408[k];

            t_321[k] = -ab_x[k] * hi_321[k]
                       + hk_409[k];

            t_322[k] = -ab_x[k] * hi_322[k]
                       + hk_410[k];

            t_323[k] = -ab_x[k] * hi_323[k]
                       + hk_411[k];

            t_324[k] = -ab_x[k] * hi_324[k]
                       + hk_412[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, hi_325, hi_326, hi_327, \
                         hi_328, hi_329, hk_413, hk_414, hk_415, hk_416, \
                         hk_417 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_x[k] * hi_325[k]
                       + hk_413[k];

            t_326[k] = -ab_x[k] * hi_326[k]
                       + hk_414[k];

            t_327[k] = -ab_x[k] * hi_327[k]
                       + hk_415[k];

            t_328[k] = -ab_x[k] * hi_328[k]
                       + hk_416[k];

            t_329[k] = -ab_x[k] * hi_329[k]
                       + hk_417[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, hi_330, hi_331, hi_332, \
                         hi_333, hi_334, hk_418, hk_419, hk_420, hk_421, \
                         hk_422 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = -ab_x[k] * hi_330[k]
                       + hk_418[k];

            t_331[k] = -ab_x[k] * hi_331[k]
                       + hk_419[k];

            t_332[k] = -ab_x[k] * hi_332[k]
                       + hk_420[k];

            t_333[k] = -ab_x[k] * hi_333[k]
                       + hk_421[k];

            t_334[k] = -ab_x[k] * hi_334[k]
                       + hk_422[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, hi_335, hi_336, hi_337, \
                         hi_338, hi_339, hk_423, hk_432, hk_433, hk_434, \
                         hk_435 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = -ab_x[k] * hi_335[k]
                       + hk_423[k];

            t_336[k] = -ab_x[k] * hi_336[k]
                       + hk_432[k];

            t_337[k] = -ab_x[k] * hi_337[k]
                       + hk_433[k];

            t_338[k] = -ab_x[k] * hi_338[k]
                       + hk_434[k];

            t_339[k] = -ab_x[k] * hi_339[k]
                       + hk_435[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, hi_340, hi_341, hi_342, \
                         hi_343, hi_344, hk_436, hk_437, hk_438, hk_439, \
                         hk_440 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = -ab_x[k] * hi_340[k]
                       + hk_436[k];

            t_341[k] = -ab_x[k] * hi_341[k]
                       + hk_437[k];

            t_342[k] = -ab_x[k] * hi_342[k]
                       + hk_438[k];

            t_343[k] = -ab_x[k] * hi_343[k]
                       + hk_439[k];

            t_344[k] = -ab_x[k] * hi_344[k]
                       + hk_440[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, hi_345, hi_346, hi_347, \
                         hi_348, hi_349, hk_441, hk_442, hk_443, hk_444, \
                         hk_445 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = -ab_x[k] * hi_345[k]
                       + hk_441[k];

            t_346[k] = -ab_x[k] * hi_346[k]
                       + hk_442[k];

            t_347[k] = -ab_x[k] * hi_347[k]
                       + hk_443[k];

            t_348[k] = -ab_x[k] * hi_348[k]
                       + hk_444[k];

            t_349[k] = -ab_x[k] * hi_349[k]
                       + hk_445[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, hi_350, hi_351, hi_352, \
                         hi_353, hi_354, hk_446, hk_447, hk_448, hk_449, \
                         hk_450 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = -ab_x[k] * hi_350[k]
                       + hk_446[k];

            t_351[k] = -ab_x[k] * hi_351[k]
                       + hk_447[k];

            t_352[k] = -ab_x[k] * hi_352[k]
                       + hk_448[k];

            t_353[k] = -ab_x[k] * hi_353[k]
                       + hk_449[k];

            t_354[k] = -ab_x[k] * hi_354[k]
                       + hk_450[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, hi_355, hi_356, hi_357, \
                         hi_358, hi_359, hk_451, hk_452, hk_453, hk_454, \
                         hk_455 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = -ab_x[k] * hi_355[k]
                       + hk_451[k];

            t_356[k] = -ab_x[k] * hi_356[k]
                       + hk_452[k];

            t_357[k] = -ab_x[k] * hi_357[k]
                       + hk_453[k];

            t_358[k] = -ab_x[k] * hi_358[k]
                       + hk_454[k];

            t_359[k] = -ab_x[k] * hi_359[k]
                       + hk_455[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, hi_360, hi_361, hi_362, \
                         hi_363, hi_364, hk_456, hk_457, hk_458, hk_459, \
                         hk_468 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = -ab_x[k] * hi_360[k]
                       + hk_456[k];

            t_361[k] = -ab_x[k] * hi_361[k]
                       + hk_457[k];

            t_362[k] = -ab_x[k] * hi_362[k]
                       + hk_458[k];

            t_363[k] = -ab_x[k] * hi_363[k]
                       + hk_459[k];

            t_364[k] = -ab_x[k] * hi_364[k]
                       + hk_468[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, hi_365, hi_366, hi_367, \
                         hi_368, hi_369, hk_469, hk_470, hk_471, hk_472, \
                         hk_473 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = -ab_x[k] * hi_365[k]
                       + hk_469[k];

            t_366[k] = -ab_x[k] * hi_366[k]
                       + hk_470[k];

            t_367[k] = -ab_x[k] * hi_367[k]
                       + hk_471[k];

            t_368[k] = -ab_x[k] * hi_368[k]
                       + hk_472[k];

            t_369[k] = -ab_x[k] * hi_369[k]
                       + hk_473[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, hi_370, hi_371, hi_372, \
                         hi_373, hi_374, hk_474, hk_475, hk_476, hk_477, \
                         hk_478 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = -ab_x[k] * hi_370[k]
                       + hk_474[k];

            t_371[k] = -ab_x[k] * hi_371[k]
                       + hk_475[k];

            t_372[k] = -ab_x[k] * hi_372[k]
                       + hk_476[k];

            t_373[k] = -ab_x[k] * hi_373[k]
                       + hk_477[k];

            t_374[k] = -ab_x[k] * hi_374[k]
                       + hk_478[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, hi_375, hi_376, hi_377, \
                         hi_378, hi_379, hk_479, hk_480, hk_481, hk_482, \
                         hk_483 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = -ab_x[k] * hi_375[k]
                       + hk_479[k];

            t_376[k] = -ab_x[k] * hi_376[k]
                       + hk_480[k];

            t_377[k] = -ab_x[k] * hi_377[k]
                       + hk_481[k];

            t_378[k] = -ab_x[k] * hi_378[k]
                       + hk_482[k];

            t_379[k] = -ab_x[k] * hi_379[k]
                       + hk_483[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, hi_380, hi_381, hi_382, \
                         hi_383, hi_384, hk_484, hk_485, hk_486, hk_487, \
                         hk_488 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = -ab_x[k] * hi_380[k]
                       + hk_484[k];

            t_381[k] = -ab_x[k] * hi_381[k]
                       + hk_485[k];

            t_382[k] = -ab_x[k] * hi_382[k]
                       + hk_486[k];

            t_383[k] = -ab_x[k] * hi_383[k]
                       + hk_487[k];

            t_384[k] = -ab_x[k] * hi_384[k]
                       + hk_488[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, hi_385, hi_386, hi_387, \
                         hi_388, hi_389, hk_489, hk_490, hk_491, hk_492, \
                         hk_493 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = -ab_x[k] * hi_385[k]
                       + hk_489[k];

            t_386[k] = -ab_x[k] * hi_386[k]
                       + hk_490[k];

            t_387[k] = -ab_x[k] * hi_387[k]
                       + hk_491[k];

            t_388[k] = -ab_x[k] * hi_388[k]
                       + hk_492[k];

            t_389[k] = -ab_x[k] * hi_389[k]
                       + hk_493[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, hi_390, hi_391, hi_392, \
                         hi_393, hi_394, hk_494, hk_495, hk_504, hk_505, \
                         hk_506 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = -ab_x[k] * hi_390[k]
                       + hk_494[k];

            t_391[k] = -ab_x[k] * hi_391[k]
                       + hk_495[k];

            t_392[k] = -ab_x[k] * hi_392[k]
                       + hk_504[k];

            t_393[k] = -ab_x[k] * hi_393[k]
                       + hk_505[k];

            t_394[k] = -ab_x[k] * hi_394[k]
                       + hk_506[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, hi_395, hi_396, hi_397, \
                         hi_398, hi_399, hk_507, hk_508, hk_509, hk_510, \
                         hk_511 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = -ab_x[k] * hi_395[k]
                       + hk_507[k];

            t_396[k] = -ab_x[k] * hi_396[k]
                       + hk_508[k];

            t_397[k] = -ab_x[k] * hi_397[k]
                       + hk_509[k];

            t_398[k] = -ab_x[k] * hi_398[k]
                       + hk_510[k];

            t_399[k] = -ab_x[k] * hi_399[k]
                       + hk_511[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, hi_400, hi_401, hi_402, \
                         hi_403, hi_404, hk_512, hk_513, hk_514, hk_515, \
                         hk_516 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = -ab_x[k] * hi_400[k]
                       + hk_512[k];

            t_401[k] = -ab_x[k] * hi_401[k]
                       + hk_513[k];

            t_402[k] = -ab_x[k] * hi_402[k]
                       + hk_514[k];

            t_403[k] = -ab_x[k] * hi_403[k]
                       + hk_515[k];

            t_404[k] = -ab_x[k] * hi_404[k]
                       + hk_516[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, hi_405, hi_406, hi_407, \
                         hi_408, hi_409, hk_517, hk_518, hk_519, hk_520, \
                         hk_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = -ab_x[k] * hi_405[k]
                       + hk_517[k];

            t_406[k] = -ab_x[k] * hi_406[k]
                       + hk_518[k];

            t_407[k] = -ab_x[k] * hi_407[k]
                       + hk_519[k];

            t_408[k] = -ab_x[k] * hi_408[k]
                       + hk_520[k];

            t_409[k] = -ab_x[k] * hi_409[k]
                       + hk_521[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, hi_410, hi_411, hi_412, \
                         hi_413, hi_414, hk_522, hk_523, hk_524, hk_525, \
                         hk_526 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = -ab_x[k] * hi_410[k]
                       + hk_522[k];

            t_411[k] = -ab_x[k] * hi_411[k]
                       + hk_523[k];

            t_412[k] = -ab_x[k] * hi_412[k]
                       + hk_524[k];

            t_413[k] = -ab_x[k] * hi_413[k]
                       + hk_525[k];

            t_414[k] = -ab_x[k] * hi_414[k]
                       + hk_526[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, hi_415, hi_416, hi_417, \
                         hi_418, hi_419, hk_527, hk_528, hk_529, hk_530, \
                         hk_531 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = -ab_x[k] * hi_415[k]
                       + hk_527[k];

            t_416[k] = -ab_x[k] * hi_416[k]
                       + hk_528[k];

            t_417[k] = -ab_x[k] * hi_417[k]
                       + hk_529[k];

            t_418[k] = -ab_x[k] * hi_418[k]
                       + hk_530[k];

            t_419[k] = -ab_x[k] * hi_419[k]
                       + hk_531[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, hi_420, hi_421, hi_422, \
                         hi_423, hi_424, hk_540, hk_541, hk_542, hk_543, \
                         hk_544 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = -ab_x[k] * hi_420[k]
                       + hk_540[k];

            t_421[k] = -ab_x[k] * hi_421[k]
                       + hk_541[k];

            t_422[k] = -ab_x[k] * hi_422[k]
                       + hk_542[k];

            t_423[k] = -ab_x[k] * hi_423[k]
                       + hk_543[k];

            t_424[k] = -ab_x[k] * hi_424[k]
                       + hk_544[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, hi_425, hi_426, hi_427, \
                         hi_428, hi_429, hk_545, hk_546, hk_547, hk_548, \
                         hk_549 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = -ab_x[k] * hi_425[k]
                       + hk_545[k];

            t_426[k] = -ab_x[k] * hi_426[k]
                       + hk_546[k];

            t_427[k] = -ab_x[k] * hi_427[k]
                       + hk_547[k];

            t_428[k] = -ab_x[k] * hi_428[k]
                       + hk_548[k];

            t_429[k] = -ab_x[k] * hi_429[k]
                       + hk_549[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, hi_430, hi_431, hi_432, \
                         hi_433, hi_434, hk_550, hk_551, hk_552, hk_553, \
                         hk_554 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = -ab_x[k] * hi_430[k]
                       + hk_550[k];

            t_431[k] = -ab_x[k] * hi_431[k]
                       + hk_551[k];

            t_432[k] = -ab_x[k] * hi_432[k]
                       + hk_552[k];

            t_433[k] = -ab_x[k] * hi_433[k]
                       + hk_553[k];

            t_434[k] = -ab_x[k] * hi_434[k]
                       + hk_554[k];
        }
    }
}

static auto
compute_hrr_ii_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hi, const size_t hk, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_435 = buffer.data(target + 435 * ncomps + c);
        auto *t_436 = buffer.data(target + 436 * ncomps + c);
        auto *t_437 = buffer.data(target + 437 * ncomps + c);
        auto *t_438 = buffer.data(target + 438 * ncomps + c);
        auto *t_439 = buffer.data(target + 439 * ncomps + c);
        auto *t_440 = buffer.data(target + 440 * ncomps + c);
        auto *t_441 = buffer.data(target + 441 * ncomps + c);
        auto *t_442 = buffer.data(target + 442 * ncomps + c);
        auto *t_443 = buffer.data(target + 443 * ncomps + c);
        auto *t_444 = buffer.data(target + 444 * ncomps + c);
        auto *t_445 = buffer.data(target + 445 * ncomps + c);
        auto *t_446 = buffer.data(target + 446 * ncomps + c);
        auto *t_447 = buffer.data(target + 447 * ncomps + c);
        auto *t_448 = buffer.data(target + 448 * ncomps + c);
        auto *t_449 = buffer.data(target + 449 * ncomps + c);
        auto *t_450 = buffer.data(target + 450 * ncomps + c);
        auto *t_451 = buffer.data(target + 451 * ncomps + c);
        auto *t_452 = buffer.data(target + 452 * ncomps + c);
        auto *t_453 = buffer.data(target + 453 * ncomps + c);
        auto *t_454 = buffer.data(target + 454 * ncomps + c);
        auto *t_455 = buffer.data(target + 455 * ncomps + c);
        auto *t_456 = buffer.data(target + 456 * ncomps + c);
        auto *t_457 = buffer.data(target + 457 * ncomps + c);
        auto *t_458 = buffer.data(target + 458 * ncomps + c);
        auto *t_459 = buffer.data(target + 459 * ncomps + c);
        auto *t_460 = buffer.data(target + 460 * ncomps + c);
        auto *t_461 = buffer.data(target + 461 * ncomps + c);
        auto *t_462 = buffer.data(target + 462 * ncomps + c);
        auto *t_463 = buffer.data(target + 463 * ncomps + c);
        auto *t_464 = buffer.data(target + 464 * ncomps + c);
        auto *t_465 = buffer.data(target + 465 * ncomps + c);
        auto *t_466 = buffer.data(target + 466 * ncomps + c);
        auto *t_467 = buffer.data(target + 467 * ncomps + c);
        auto *t_468 = buffer.data(target + 468 * ncomps + c);
        auto *t_469 = buffer.data(target + 469 * ncomps + c);
        auto *t_470 = buffer.data(target + 470 * ncomps + c);
        auto *t_471 = buffer.data(target + 471 * ncomps + c);
        auto *t_472 = buffer.data(target + 472 * ncomps + c);
        auto *t_473 = buffer.data(target + 473 * ncomps + c);
        auto *t_474 = buffer.data(target + 474 * ncomps + c);
        auto *t_475 = buffer.data(target + 475 * ncomps + c);
        auto *t_476 = buffer.data(target + 476 * ncomps + c);
        auto *t_477 = buffer.data(target + 477 * ncomps + c);
        auto *t_478 = buffer.data(target + 478 * ncomps + c);
        auto *t_479 = buffer.data(target + 479 * ncomps + c);
        auto *t_480 = buffer.data(target + 480 * ncomps + c);
        auto *t_481 = buffer.data(target + 481 * ncomps + c);
        auto *t_482 = buffer.data(target + 482 * ncomps + c);
        auto *t_483 = buffer.data(target + 483 * ncomps + c);
        auto *t_484 = buffer.data(target + 484 * ncomps + c);
        auto *t_485 = buffer.data(target + 485 * ncomps + c);
        auto *t_486 = buffer.data(target + 486 * ncomps + c);
        auto *t_487 = buffer.data(target + 487 * ncomps + c);
        auto *t_488 = buffer.data(target + 488 * ncomps + c);
        auto *t_489 = buffer.data(target + 489 * ncomps + c);
        auto *t_490 = buffer.data(target + 490 * ncomps + c);
        auto *t_491 = buffer.data(target + 491 * ncomps + c);
        auto *t_492 = buffer.data(target + 492 * ncomps + c);
        auto *t_493 = buffer.data(target + 493 * ncomps + c);
        auto *t_494 = buffer.data(target + 494 * ncomps + c);
        auto *t_495 = buffer.data(target + 495 * ncomps + c);
        auto *t_496 = buffer.data(target + 496 * ncomps + c);
        auto *t_497 = buffer.data(target + 497 * ncomps + c);
        auto *t_498 = buffer.data(target + 498 * ncomps + c);
        auto *t_499 = buffer.data(target + 499 * ncomps + c);
        auto *t_500 = buffer.data(target + 500 * ncomps + c);
        auto *t_501 = buffer.data(target + 501 * ncomps + c);
        auto *t_502 = buffer.data(target + 502 * ncomps + c);
        auto *t_503 = buffer.data(target + 503 * ncomps + c);
        auto *t_504 = buffer.data(target + 504 * ncomps + c);
        auto *t_505 = buffer.data(target + 505 * ncomps + c);
        auto *t_506 = buffer.data(target + 506 * ncomps + c);
        auto *t_507 = buffer.data(target + 507 * ncomps + c);
        auto *t_508 = buffer.data(target + 508 * ncomps + c);
        auto *t_509 = buffer.data(target + 509 * ncomps + c);
        auto *t_510 = buffer.data(target + 510 * ncomps + c);
        auto *t_511 = buffer.data(target + 511 * ncomps + c);
        auto *t_512 = buffer.data(target + 512 * ncomps + c);
        auto *t_513 = buffer.data(target + 513 * ncomps + c);
        auto *t_514 = buffer.data(target + 514 * ncomps + c);
        auto *t_515 = buffer.data(target + 515 * ncomps + c);
        auto *t_516 = buffer.data(target + 516 * ncomps + c);
        auto *t_517 = buffer.data(target + 517 * ncomps + c);
        auto *t_518 = buffer.data(target + 518 * ncomps + c);
        auto *t_519 = buffer.data(target + 519 * ncomps + c);
        auto *t_520 = buffer.data(target + 520 * ncomps + c);
        auto *t_521 = buffer.data(target + 521 * ncomps + c);
        auto *t_522 = buffer.data(target + 522 * ncomps + c);
        auto *t_523 = buffer.data(target + 523 * ncomps + c);
        auto *t_524 = buffer.data(target + 524 * ncomps + c);
        auto *t_525 = buffer.data(target + 525 * ncomps + c);
        auto *t_526 = buffer.data(target + 526 * ncomps + c);
        auto *t_527 = buffer.data(target + 527 * ncomps + c);
        auto *t_528 = buffer.data(target + 528 * ncomps + c);
        auto *t_529 = buffer.data(target + 529 * ncomps + c);
        auto *t_530 = buffer.data(target + 530 * ncomps + c);
        auto *t_531 = buffer.data(target + 531 * ncomps + c);
        auto *t_532 = buffer.data(target + 532 * ncomps + c);
        auto *t_533 = buffer.data(target + 533 * ncomps + c);
        auto *t_534 = buffer.data(target + 534 * ncomps + c);
        auto *t_535 = buffer.data(target + 535 * ncomps + c);
        auto *t_536 = buffer.data(target + 536 * ncomps + c);
        auto *t_537 = buffer.data(target + 537 * ncomps + c);
        auto *t_538 = buffer.data(target + 538 * ncomps + c);
        auto *t_539 = buffer.data(target + 539 * ncomps + c);
        auto *t_540 = buffer.data(target + 540 * ncomps + c);
        auto *t_541 = buffer.data(target + 541 * ncomps + c);
        auto *t_542 = buffer.data(target + 542 * ncomps + c);
        auto *t_543 = buffer.data(target + 543 * ncomps + c);
        auto *t_544 = buffer.data(target + 544 * ncomps + c);
        auto *t_545 = buffer.data(target + 545 * ncomps + c);
        auto *t_546 = buffer.data(target + 546 * ncomps + c);
        auto *t_547 = buffer.data(target + 547 * ncomps + c);
        auto *t_548 = buffer.data(target + 548 * ncomps + c);
        auto *t_549 = buffer.data(target + 549 * ncomps + c);
        auto *t_550 = buffer.data(target + 550 * ncomps + c);
        auto *t_551 = buffer.data(target + 551 * ncomps + c);
        auto *t_552 = buffer.data(target + 552 * ncomps + c);
        auto *t_553 = buffer.data(target + 553 * ncomps + c);
        auto *t_554 = buffer.data(target + 554 * ncomps + c);
        auto *t_555 = buffer.data(target + 555 * ncomps + c);
        auto *t_556 = buffer.data(target + 556 * ncomps + c);
        auto *t_557 = buffer.data(target + 557 * ncomps + c);
        auto *t_558 = buffer.data(target + 558 * ncomps + c);
        auto *t_559 = buffer.data(target + 559 * ncomps + c);
        auto *t_560 = buffer.data(target + 560 * ncomps + c);
        auto *t_561 = buffer.data(target + 561 * ncomps + c);
        auto *t_562 = buffer.data(target + 562 * ncomps + c);
        auto *t_563 = buffer.data(target + 563 * ncomps + c);
        auto *t_564 = buffer.data(target + 564 * ncomps + c);
        auto *t_565 = buffer.data(target + 565 * ncomps + c);
        auto *t_566 = buffer.data(target + 566 * ncomps + c);
        auto *t_567 = buffer.data(target + 567 * ncomps + c);
        auto *t_568 = buffer.data(target + 568 * ncomps + c);
        auto *t_569 = buffer.data(target + 569 * ncomps + c);
        auto *t_570 = buffer.data(target + 570 * ncomps + c);
        auto *t_571 = buffer.data(target + 571 * ncomps + c);
        auto *t_572 = buffer.data(target + 572 * ncomps + c);
        auto *t_573 = buffer.data(target + 573 * ncomps + c);
        auto *t_574 = buffer.data(target + 574 * ncomps + c);
        auto *t_575 = buffer.data(target + 575 * ncomps + c);
        auto *t_576 = buffer.data(target + 576 * ncomps + c);
        auto *t_577 = buffer.data(target + 577 * ncomps + c);
        auto *t_578 = buffer.data(target + 578 * ncomps + c);
        auto *t_579 = buffer.data(target + 579 * ncomps + c);

        const auto *ab_x = coordinates.data(6);

        const auto *hi_435 = buffer.data(hi + 435 * ncomps + c);
        const auto *hi_436 = buffer.data(hi + 436 * ncomps + c);
        const auto *hi_437 = buffer.data(hi + 437 * ncomps + c);
        const auto *hi_438 = buffer.data(hi + 438 * ncomps + c);
        const auto *hi_439 = buffer.data(hi + 439 * ncomps + c);
        const auto *hi_440 = buffer.data(hi + 440 * ncomps + c);
        const auto *hi_441 = buffer.data(hi + 441 * ncomps + c);
        const auto *hi_442 = buffer.data(hi + 442 * ncomps + c);
        const auto *hi_443 = buffer.data(hi + 443 * ncomps + c);
        const auto *hi_444 = buffer.data(hi + 444 * ncomps + c);
        const auto *hi_445 = buffer.data(hi + 445 * ncomps + c);
        const auto *hi_446 = buffer.data(hi + 446 * ncomps + c);
        const auto *hi_447 = buffer.data(hi + 447 * ncomps + c);
        const auto *hi_448 = buffer.data(hi + 448 * ncomps + c);
        const auto *hi_449 = buffer.data(hi + 449 * ncomps + c);
        const auto *hi_450 = buffer.data(hi + 450 * ncomps + c);
        const auto *hi_451 = buffer.data(hi + 451 * ncomps + c);
        const auto *hi_452 = buffer.data(hi + 452 * ncomps + c);
        const auto *hi_453 = buffer.data(hi + 453 * ncomps + c);
        const auto *hi_454 = buffer.data(hi + 454 * ncomps + c);
        const auto *hi_455 = buffer.data(hi + 455 * ncomps + c);
        const auto *hi_456 = buffer.data(hi + 456 * ncomps + c);
        const auto *hi_457 = buffer.data(hi + 457 * ncomps + c);
        const auto *hi_458 = buffer.data(hi + 458 * ncomps + c);
        const auto *hi_459 = buffer.data(hi + 459 * ncomps + c);
        const auto *hi_460 = buffer.data(hi + 460 * ncomps + c);
        const auto *hi_461 = buffer.data(hi + 461 * ncomps + c);
        const auto *hi_462 = buffer.data(hi + 462 * ncomps + c);
        const auto *hi_463 = buffer.data(hi + 463 * ncomps + c);
        const auto *hi_464 = buffer.data(hi + 464 * ncomps + c);
        const auto *hi_465 = buffer.data(hi + 465 * ncomps + c);
        const auto *hi_466 = buffer.data(hi + 466 * ncomps + c);
        const auto *hi_467 = buffer.data(hi + 467 * ncomps + c);
        const auto *hi_468 = buffer.data(hi + 468 * ncomps + c);
        const auto *hi_469 = buffer.data(hi + 469 * ncomps + c);
        const auto *hi_470 = buffer.data(hi + 470 * ncomps + c);
        const auto *hi_471 = buffer.data(hi + 471 * ncomps + c);
        const auto *hi_472 = buffer.data(hi + 472 * ncomps + c);
        const auto *hi_473 = buffer.data(hi + 473 * ncomps + c);
        const auto *hi_474 = buffer.data(hi + 474 * ncomps + c);
        const auto *hi_475 = buffer.data(hi + 475 * ncomps + c);
        const auto *hi_476 = buffer.data(hi + 476 * ncomps + c);
        const auto *hi_477 = buffer.data(hi + 477 * ncomps + c);
        const auto *hi_478 = buffer.data(hi + 478 * ncomps + c);
        const auto *hi_479 = buffer.data(hi + 479 * ncomps + c);
        const auto *hi_480 = buffer.data(hi + 480 * ncomps + c);
        const auto *hi_481 = buffer.data(hi + 481 * ncomps + c);
        const auto *hi_482 = buffer.data(hi + 482 * ncomps + c);
        const auto *hi_483 = buffer.data(hi + 483 * ncomps + c);
        const auto *hi_484 = buffer.data(hi + 484 * ncomps + c);
        const auto *hi_485 = buffer.data(hi + 485 * ncomps + c);
        const auto *hi_486 = buffer.data(hi + 486 * ncomps + c);
        const auto *hi_487 = buffer.data(hi + 487 * ncomps + c);
        const auto *hi_488 = buffer.data(hi + 488 * ncomps + c);
        const auto *hi_489 = buffer.data(hi + 489 * ncomps + c);
        const auto *hi_490 = buffer.data(hi + 490 * ncomps + c);
        const auto *hi_491 = buffer.data(hi + 491 * ncomps + c);
        const auto *hi_492 = buffer.data(hi + 492 * ncomps + c);
        const auto *hi_493 = buffer.data(hi + 493 * ncomps + c);
        const auto *hi_494 = buffer.data(hi + 494 * ncomps + c);
        const auto *hi_495 = buffer.data(hi + 495 * ncomps + c);
        const auto *hi_496 = buffer.data(hi + 496 * ncomps + c);
        const auto *hi_497 = buffer.data(hi + 497 * ncomps + c);
        const auto *hi_498 = buffer.data(hi + 498 * ncomps + c);
        const auto *hi_499 = buffer.data(hi + 499 * ncomps + c);
        const auto *hi_500 = buffer.data(hi + 500 * ncomps + c);
        const auto *hi_501 = buffer.data(hi + 501 * ncomps + c);
        const auto *hi_502 = buffer.data(hi + 502 * ncomps + c);
        const auto *hi_503 = buffer.data(hi + 503 * ncomps + c);
        const auto *hi_504 = buffer.data(hi + 504 * ncomps + c);
        const auto *hi_505 = buffer.data(hi + 505 * ncomps + c);
        const auto *hi_506 = buffer.data(hi + 506 * ncomps + c);
        const auto *hi_507 = buffer.data(hi + 507 * ncomps + c);
        const auto *hi_508 = buffer.data(hi + 508 * ncomps + c);
        const auto *hi_509 = buffer.data(hi + 509 * ncomps + c);
        const auto *hi_510 = buffer.data(hi + 510 * ncomps + c);
        const auto *hi_511 = buffer.data(hi + 511 * ncomps + c);
        const auto *hi_512 = buffer.data(hi + 512 * ncomps + c);
        const auto *hi_513 = buffer.data(hi + 513 * ncomps + c);
        const auto *hi_514 = buffer.data(hi + 514 * ncomps + c);
        const auto *hi_515 = buffer.data(hi + 515 * ncomps + c);
        const auto *hi_516 = buffer.data(hi + 516 * ncomps + c);
        const auto *hi_517 = buffer.data(hi + 517 * ncomps + c);
        const auto *hi_518 = buffer.data(hi + 518 * ncomps + c);
        const auto *hi_519 = buffer.data(hi + 519 * ncomps + c);
        const auto *hi_520 = buffer.data(hi + 520 * ncomps + c);
        const auto *hi_521 = buffer.data(hi + 521 * ncomps + c);
        const auto *hi_522 = buffer.data(hi + 522 * ncomps + c);
        const auto *hi_523 = buffer.data(hi + 523 * ncomps + c);
        const auto *hi_524 = buffer.data(hi + 524 * ncomps + c);
        const auto *hi_525 = buffer.data(hi + 525 * ncomps + c);
        const auto *hi_526 = buffer.data(hi + 526 * ncomps + c);
        const auto *hi_527 = buffer.data(hi + 527 * ncomps + c);
        const auto *hi_528 = buffer.data(hi + 528 * ncomps + c);
        const auto *hi_529 = buffer.data(hi + 529 * ncomps + c);
        const auto *hi_530 = buffer.data(hi + 530 * ncomps + c);
        const auto *hi_531 = buffer.data(hi + 531 * ncomps + c);
        const auto *hi_532 = buffer.data(hi + 532 * ncomps + c);
        const auto *hi_533 = buffer.data(hi + 533 * ncomps + c);
        const auto *hi_534 = buffer.data(hi + 534 * ncomps + c);
        const auto *hi_535 = buffer.data(hi + 535 * ncomps + c);
        const auto *hi_536 = buffer.data(hi + 536 * ncomps + c);
        const auto *hi_537 = buffer.data(hi + 537 * ncomps + c);
        const auto *hi_538 = buffer.data(hi + 538 * ncomps + c);
        const auto *hi_539 = buffer.data(hi + 539 * ncomps + c);
        const auto *hi_540 = buffer.data(hi + 540 * ncomps + c);
        const auto *hi_541 = buffer.data(hi + 541 * ncomps + c);
        const auto *hi_542 = buffer.data(hi + 542 * ncomps + c);
        const auto *hi_543 = buffer.data(hi + 543 * ncomps + c);
        const auto *hi_544 = buffer.data(hi + 544 * ncomps + c);
        const auto *hi_545 = buffer.data(hi + 545 * ncomps + c);
        const auto *hi_546 = buffer.data(hi + 546 * ncomps + c);
        const auto *hi_547 = buffer.data(hi + 547 * ncomps + c);
        const auto *hi_548 = buffer.data(hi + 548 * ncomps + c);
        const auto *hi_549 = buffer.data(hi + 549 * ncomps + c);
        const auto *hi_550 = buffer.data(hi + 550 * ncomps + c);
        const auto *hi_551 = buffer.data(hi + 551 * ncomps + c);
        const auto *hi_552 = buffer.data(hi + 552 * ncomps + c);
        const auto *hi_553 = buffer.data(hi + 553 * ncomps + c);
        const auto *hi_554 = buffer.data(hi + 554 * ncomps + c);
        const auto *hi_555 = buffer.data(hi + 555 * ncomps + c);
        const auto *hi_556 = buffer.data(hi + 556 * ncomps + c);
        const auto *hi_557 = buffer.data(hi + 557 * ncomps + c);
        const auto *hi_558 = buffer.data(hi + 558 * ncomps + c);
        const auto *hi_559 = buffer.data(hi + 559 * ncomps + c);
        const auto *hi_560 = buffer.data(hi + 560 * ncomps + c);
        const auto *hi_561 = buffer.data(hi + 561 * ncomps + c);
        const auto *hi_562 = buffer.data(hi + 562 * ncomps + c);
        const auto *hi_563 = buffer.data(hi + 563 * ncomps + c);
        const auto *hi_564 = buffer.data(hi + 564 * ncomps + c);
        const auto *hi_565 = buffer.data(hi + 565 * ncomps + c);
        const auto *hi_566 = buffer.data(hi + 566 * ncomps + c);
        const auto *hi_567 = buffer.data(hi + 567 * ncomps + c);
        const auto *hi_568 = buffer.data(hi + 568 * ncomps + c);
        const auto *hi_569 = buffer.data(hi + 569 * ncomps + c);
        const auto *hi_570 = buffer.data(hi + 570 * ncomps + c);
        const auto *hi_571 = buffer.data(hi + 571 * ncomps + c);
        const auto *hi_572 = buffer.data(hi + 572 * ncomps + c);
        const auto *hi_573 = buffer.data(hi + 573 * ncomps + c);
        const auto *hi_574 = buffer.data(hi + 574 * ncomps + c);
        const auto *hi_575 = buffer.data(hi + 575 * ncomps + c);
        const auto *hi_576 = buffer.data(hi + 576 * ncomps + c);
        const auto *hi_577 = buffer.data(hi + 577 * ncomps + c);
        const auto *hi_578 = buffer.data(hi + 578 * ncomps + c);
        const auto *hi_579 = buffer.data(hi + 579 * ncomps + c);

        const auto *hk_555 = buffer.data(hk + 555 * ncomps + c);
        const auto *hk_556 = buffer.data(hk + 556 * ncomps + c);
        const auto *hk_557 = buffer.data(hk + 557 * ncomps + c);
        const auto *hk_558 = buffer.data(hk + 558 * ncomps + c);
        const auto *hk_559 = buffer.data(hk + 559 * ncomps + c);
        const auto *hk_560 = buffer.data(hk + 560 * ncomps + c);
        const auto *hk_561 = buffer.data(hk + 561 * ncomps + c);
        const auto *hk_562 = buffer.data(hk + 562 * ncomps + c);
        const auto *hk_563 = buffer.data(hk + 563 * ncomps + c);
        const auto *hk_564 = buffer.data(hk + 564 * ncomps + c);
        const auto *hk_565 = buffer.data(hk + 565 * ncomps + c);
        const auto *hk_566 = buffer.data(hk + 566 * ncomps + c);
        const auto *hk_567 = buffer.data(hk + 567 * ncomps + c);
        const auto *hk_576 = buffer.data(hk + 576 * ncomps + c);
        const auto *hk_577 = buffer.data(hk + 577 * ncomps + c);
        const auto *hk_578 = buffer.data(hk + 578 * ncomps + c);
        const auto *hk_579 = buffer.data(hk + 579 * ncomps + c);
        const auto *hk_580 = buffer.data(hk + 580 * ncomps + c);
        const auto *hk_581 = buffer.data(hk + 581 * ncomps + c);
        const auto *hk_582 = buffer.data(hk + 582 * ncomps + c);
        const auto *hk_583 = buffer.data(hk + 583 * ncomps + c);
        const auto *hk_584 = buffer.data(hk + 584 * ncomps + c);
        const auto *hk_585 = buffer.data(hk + 585 * ncomps + c);
        const auto *hk_586 = buffer.data(hk + 586 * ncomps + c);
        const auto *hk_587 = buffer.data(hk + 587 * ncomps + c);
        const auto *hk_588 = buffer.data(hk + 588 * ncomps + c);
        const auto *hk_589 = buffer.data(hk + 589 * ncomps + c);
        const auto *hk_590 = buffer.data(hk + 590 * ncomps + c);
        const auto *hk_591 = buffer.data(hk + 591 * ncomps + c);
        const auto *hk_592 = buffer.data(hk + 592 * ncomps + c);
        const auto *hk_593 = buffer.data(hk + 593 * ncomps + c);
        const auto *hk_594 = buffer.data(hk + 594 * ncomps + c);
        const auto *hk_595 = buffer.data(hk + 595 * ncomps + c);
        const auto *hk_596 = buffer.data(hk + 596 * ncomps + c);
        const auto *hk_597 = buffer.data(hk + 597 * ncomps + c);
        const auto *hk_598 = buffer.data(hk + 598 * ncomps + c);
        const auto *hk_599 = buffer.data(hk + 599 * ncomps + c);
        const auto *hk_600 = buffer.data(hk + 600 * ncomps + c);
        const auto *hk_601 = buffer.data(hk + 601 * ncomps + c);
        const auto *hk_602 = buffer.data(hk + 602 * ncomps + c);
        const auto *hk_603 = buffer.data(hk + 603 * ncomps + c);
        const auto *hk_612 = buffer.data(hk + 612 * ncomps + c);
        const auto *hk_613 = buffer.data(hk + 613 * ncomps + c);
        const auto *hk_614 = buffer.data(hk + 614 * ncomps + c);
        const auto *hk_615 = buffer.data(hk + 615 * ncomps + c);
        const auto *hk_616 = buffer.data(hk + 616 * ncomps + c);
        const auto *hk_617 = buffer.data(hk + 617 * ncomps + c);
        const auto *hk_618 = buffer.data(hk + 618 * ncomps + c);
        const auto *hk_619 = buffer.data(hk + 619 * ncomps + c);
        const auto *hk_620 = buffer.data(hk + 620 * ncomps + c);
        const auto *hk_621 = buffer.data(hk + 621 * ncomps + c);
        const auto *hk_622 = buffer.data(hk + 622 * ncomps + c);
        const auto *hk_623 = buffer.data(hk + 623 * ncomps + c);
        const auto *hk_624 = buffer.data(hk + 624 * ncomps + c);
        const auto *hk_625 = buffer.data(hk + 625 * ncomps + c);
        const auto *hk_626 = buffer.data(hk + 626 * ncomps + c);
        const auto *hk_627 = buffer.data(hk + 627 * ncomps + c);
        const auto *hk_628 = buffer.data(hk + 628 * ncomps + c);
        const auto *hk_629 = buffer.data(hk + 629 * ncomps + c);
        const auto *hk_630 = buffer.data(hk + 630 * ncomps + c);
        const auto *hk_631 = buffer.data(hk + 631 * ncomps + c);
        const auto *hk_632 = buffer.data(hk + 632 * ncomps + c);
        const auto *hk_633 = buffer.data(hk + 633 * ncomps + c);
        const auto *hk_634 = buffer.data(hk + 634 * ncomps + c);
        const auto *hk_635 = buffer.data(hk + 635 * ncomps + c);
        const auto *hk_636 = buffer.data(hk + 636 * ncomps + c);
        const auto *hk_637 = buffer.data(hk + 637 * ncomps + c);
        const auto *hk_638 = buffer.data(hk + 638 * ncomps + c);
        const auto *hk_639 = buffer.data(hk + 639 * ncomps + c);
        const auto *hk_648 = buffer.data(hk + 648 * ncomps + c);
        const auto *hk_649 = buffer.data(hk + 649 * ncomps + c);
        const auto *hk_650 = buffer.data(hk + 650 * ncomps + c);
        const auto *hk_651 = buffer.data(hk + 651 * ncomps + c);
        const auto *hk_652 = buffer.data(hk + 652 * ncomps + c);
        const auto *hk_653 = buffer.data(hk + 653 * ncomps + c);
        const auto *hk_654 = buffer.data(hk + 654 * ncomps + c);
        const auto *hk_655 = buffer.data(hk + 655 * ncomps + c);
        const auto *hk_656 = buffer.data(hk + 656 * ncomps + c);
        const auto *hk_657 = buffer.data(hk + 657 * ncomps + c);
        const auto *hk_658 = buffer.data(hk + 658 * ncomps + c);
        const auto *hk_659 = buffer.data(hk + 659 * ncomps + c);
        const auto *hk_660 = buffer.data(hk + 660 * ncomps + c);
        const auto *hk_661 = buffer.data(hk + 661 * ncomps + c);
        const auto *hk_662 = buffer.data(hk + 662 * ncomps + c);
        const auto *hk_663 = buffer.data(hk + 663 * ncomps + c);
        const auto *hk_664 = buffer.data(hk + 664 * ncomps + c);
        const auto *hk_665 = buffer.data(hk + 665 * ncomps + c);
        const auto *hk_666 = buffer.data(hk + 666 * ncomps + c);
        const auto *hk_667 = buffer.data(hk + 667 * ncomps + c);
        const auto *hk_668 = buffer.data(hk + 668 * ncomps + c);
        const auto *hk_669 = buffer.data(hk + 669 * ncomps + c);
        const auto *hk_670 = buffer.data(hk + 670 * ncomps + c);
        const auto *hk_671 = buffer.data(hk + 671 * ncomps + c);
        const auto *hk_672 = buffer.data(hk + 672 * ncomps + c);
        const auto *hk_673 = buffer.data(hk + 673 * ncomps + c);
        const auto *hk_674 = buffer.data(hk + 674 * ncomps + c);
        const auto *hk_675 = buffer.data(hk + 675 * ncomps + c);
        const auto *hk_684 = buffer.data(hk + 684 * ncomps + c);
        const auto *hk_685 = buffer.data(hk + 685 * ncomps + c);
        const auto *hk_686 = buffer.data(hk + 686 * ncomps + c);
        const auto *hk_687 = buffer.data(hk + 687 * ncomps + c);
        const auto *hk_688 = buffer.data(hk + 688 * ncomps + c);
        const auto *hk_689 = buffer.data(hk + 689 * ncomps + c);
        const auto *hk_690 = buffer.data(hk + 690 * ncomps + c);
        const auto *hk_691 = buffer.data(hk + 691 * ncomps + c);
        const auto *hk_692 = buffer.data(hk + 692 * ncomps + c);
        const auto *hk_693 = buffer.data(hk + 693 * ncomps + c);
        const auto *hk_694 = buffer.data(hk + 694 * ncomps + c);
        const auto *hk_695 = buffer.data(hk + 695 * ncomps + c);
        const auto *hk_696 = buffer.data(hk + 696 * ncomps + c);
        const auto *hk_697 = buffer.data(hk + 697 * ncomps + c);
        const auto *hk_698 = buffer.data(hk + 698 * ncomps + c);
        const auto *hk_699 = buffer.data(hk + 699 * ncomps + c);
        const auto *hk_700 = buffer.data(hk + 700 * ncomps + c);
        const auto *hk_701 = buffer.data(hk + 701 * ncomps + c);
        const auto *hk_702 = buffer.data(hk + 702 * ncomps + c);
        const auto *hk_703 = buffer.data(hk + 703 * ncomps + c);
        const auto *hk_704 = buffer.data(hk + 704 * ncomps + c);
        const auto *hk_705 = buffer.data(hk + 705 * ncomps + c);
        const auto *hk_706 = buffer.data(hk + 706 * ncomps + c);
        const auto *hk_707 = buffer.data(hk + 707 * ncomps + c);
        const auto *hk_708 = buffer.data(hk + 708 * ncomps + c);
        const auto *hk_709 = buffer.data(hk + 709 * ncomps + c);
        const auto *hk_710 = buffer.data(hk + 710 * ncomps + c);
        const auto *hk_711 = buffer.data(hk + 711 * ncomps + c);
        const auto *hk_720 = buffer.data(hk + 720 * ncomps + c);
        const auto *hk_721 = buffer.data(hk + 721 * ncomps + c);
        const auto *hk_722 = buffer.data(hk + 722 * ncomps + c);
        const auto *hk_723 = buffer.data(hk + 723 * ncomps + c);
        const auto *hk_724 = buffer.data(hk + 724 * ncomps + c);
        const auto *hk_725 = buffer.data(hk + 725 * ncomps + c);
        const auto *hk_726 = buffer.data(hk + 726 * ncomps + c);
        const auto *hk_727 = buffer.data(hk + 727 * ncomps + c);
        const auto *hk_728 = buffer.data(hk + 728 * ncomps + c);
        const auto *hk_729 = buffer.data(hk + 729 * ncomps + c);
        const auto *hk_730 = buffer.data(hk + 730 * ncomps + c);
        const auto *hk_731 = buffer.data(hk + 731 * ncomps + c);
        const auto *hk_732 = buffer.data(hk + 732 * ncomps + c);
        const auto *hk_733 = buffer.data(hk + 733 * ncomps + c);
        const auto *hk_734 = buffer.data(hk + 734 * ncomps + c);
        const auto *hk_735 = buffer.data(hk + 735 * ncomps + c);
        const auto *hk_736 = buffer.data(hk + 736 * ncomps + c);
        const auto *hk_737 = buffer.data(hk + 737 * ncomps + c);
        const auto *hk_738 = buffer.data(hk + 738 * ncomps + c);
        const auto *hk_739 = buffer.data(hk + 739 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, hi_435, hi_436, hi_437, \
                         hi_438, hi_439, hk_555, hk_556, hk_557, hk_558, \
                         hk_559 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = -ab_x[k] * hi_435[k]
                       + hk_555[k];

            t_436[k] = -ab_x[k] * hi_436[k]
                       + hk_556[k];

            t_437[k] = -ab_x[k] * hi_437[k]
                       + hk_557[k];

            t_438[k] = -ab_x[k] * hi_438[k]
                       + hk_558[k];

            t_439[k] = -ab_x[k] * hi_439[k]
                       + hk_559[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, hi_440, hi_441, hi_442, \
                         hi_443, hi_444, hk_560, hk_561, hk_562, hk_563, \
                         hk_564 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = -ab_x[k] * hi_440[k]
                       + hk_560[k];

            t_441[k] = -ab_x[k] * hi_441[k]
                       + hk_561[k];

            t_442[k] = -ab_x[k] * hi_442[k]
                       + hk_562[k];

            t_443[k] = -ab_x[k] * hi_443[k]
                       + hk_563[k];

            t_444[k] = -ab_x[k] * hi_444[k]
                       + hk_564[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_x, hi_445, hi_446, hi_447, \
                         hi_448, hi_449, hk_565, hk_566, hk_567, hk_576, \
                         hk_577 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = -ab_x[k] * hi_445[k]
                       + hk_565[k];

            t_446[k] = -ab_x[k] * hi_446[k]
                       + hk_566[k];

            t_447[k] = -ab_x[k] * hi_447[k]
                       + hk_567[k];

            t_448[k] = -ab_x[k] * hi_448[k]
                       + hk_576[k];

            t_449[k] = -ab_x[k] * hi_449[k]
                       + hk_577[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_x, hi_450, hi_451, hi_452, \
                         hi_453, hi_454, hk_578, hk_579, hk_580, hk_581, \
                         hk_582 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = -ab_x[k] * hi_450[k]
                       + hk_578[k];

            t_451[k] = -ab_x[k] * hi_451[k]
                       + hk_579[k];

            t_452[k] = -ab_x[k] * hi_452[k]
                       + hk_580[k];

            t_453[k] = -ab_x[k] * hi_453[k]
                       + hk_581[k];

            t_454[k] = -ab_x[k] * hi_454[k]
                       + hk_582[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_x, hi_455, hi_456, hi_457, \
                         hi_458, hi_459, hk_583, hk_584, hk_585, hk_586, \
                         hk_587 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = -ab_x[k] * hi_455[k]
                       + hk_583[k];

            t_456[k] = -ab_x[k] * hi_456[k]
                       + hk_584[k];

            t_457[k] = -ab_x[k] * hi_457[k]
                       + hk_585[k];

            t_458[k] = -ab_x[k] * hi_458[k]
                       + hk_586[k];

            t_459[k] = -ab_x[k] * hi_459[k]
                       + hk_587[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_x, hi_460, hi_461, hi_462, \
                         hi_463, hi_464, hk_588, hk_589, hk_590, hk_591, \
                         hk_592 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = -ab_x[k] * hi_460[k]
                       + hk_588[k];

            t_461[k] = -ab_x[k] * hi_461[k]
                       + hk_589[k];

            t_462[k] = -ab_x[k] * hi_462[k]
                       + hk_590[k];

            t_463[k] = -ab_x[k] * hi_463[k]
                       + hk_591[k];

            t_464[k] = -ab_x[k] * hi_464[k]
                       + hk_592[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_x, hi_465, hi_466, hi_467, \
                         hi_468, hi_469, hk_593, hk_594, hk_595, hk_596, \
                         hk_597 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = -ab_x[k] * hi_465[k]
                       + hk_593[k];

            t_466[k] = -ab_x[k] * hi_466[k]
                       + hk_594[k];

            t_467[k] = -ab_x[k] * hi_467[k]
                       + hk_595[k];

            t_468[k] = -ab_x[k] * hi_468[k]
                       + hk_596[k];

            t_469[k] = -ab_x[k] * hi_469[k]
                       + hk_597[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_x, hi_470, hi_471, hi_472, \
                         hi_473, hi_474, hk_598, hk_599, hk_600, hk_601, \
                         hk_602 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = -ab_x[k] * hi_470[k]
                       + hk_598[k];

            t_471[k] = -ab_x[k] * hi_471[k]
                       + hk_599[k];

            t_472[k] = -ab_x[k] * hi_472[k]
                       + hk_600[k];

            t_473[k] = -ab_x[k] * hi_473[k]
                       + hk_601[k];

            t_474[k] = -ab_x[k] * hi_474[k]
                       + hk_602[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_x, hi_475, hi_476, hi_477, \
                         hi_478, hi_479, hk_603, hk_612, hk_613, hk_614, \
                         hk_615 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = -ab_x[k] * hi_475[k]
                       + hk_603[k];

            t_476[k] = -ab_x[k] * hi_476[k]
                       + hk_612[k];

            t_477[k] = -ab_x[k] * hi_477[k]
                       + hk_613[k];

            t_478[k] = -ab_x[k] * hi_478[k]
                       + hk_614[k];

            t_479[k] = -ab_x[k] * hi_479[k]
                       + hk_615[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_x, hi_480, hi_481, hi_482, \
                         hi_483, hi_484, hk_616, hk_617, hk_618, hk_619, \
                         hk_620 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = -ab_x[k] * hi_480[k]
                       + hk_616[k];

            t_481[k] = -ab_x[k] * hi_481[k]
                       + hk_617[k];

            t_482[k] = -ab_x[k] * hi_482[k]
                       + hk_618[k];

            t_483[k] = -ab_x[k] * hi_483[k]
                       + hk_619[k];

            t_484[k] = -ab_x[k] * hi_484[k]
                       + hk_620[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_x, hi_485, hi_486, hi_487, \
                         hi_488, hi_489, hk_621, hk_622, hk_623, hk_624, \
                         hk_625 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = -ab_x[k] * hi_485[k]
                       + hk_621[k];

            t_486[k] = -ab_x[k] * hi_486[k]
                       + hk_622[k];

            t_487[k] = -ab_x[k] * hi_487[k]
                       + hk_623[k];

            t_488[k] = -ab_x[k] * hi_488[k]
                       + hk_624[k];

            t_489[k] = -ab_x[k] * hi_489[k]
                       + hk_625[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_x, hi_490, hi_491, hi_492, \
                         hi_493, hi_494, hk_626, hk_627, hk_628, hk_629, \
                         hk_630 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = -ab_x[k] * hi_490[k]
                       + hk_626[k];

            t_491[k] = -ab_x[k] * hi_491[k]
                       + hk_627[k];

            t_492[k] = -ab_x[k] * hi_492[k]
                       + hk_628[k];

            t_493[k] = -ab_x[k] * hi_493[k]
                       + hk_629[k];

            t_494[k] = -ab_x[k] * hi_494[k]
                       + hk_630[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_x, hi_495, hi_496, hi_497, \
                         hi_498, hi_499, hk_631, hk_632, hk_633, hk_634, \
                         hk_635 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = -ab_x[k] * hi_495[k]
                       + hk_631[k];

            t_496[k] = -ab_x[k] * hi_496[k]
                       + hk_632[k];

            t_497[k] = -ab_x[k] * hi_497[k]
                       + hk_633[k];

            t_498[k] = -ab_x[k] * hi_498[k]
                       + hk_634[k];

            t_499[k] = -ab_x[k] * hi_499[k]
                       + hk_635[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_x, hi_500, hi_501, hi_502, \
                         hi_503, hi_504, hk_636, hk_637, hk_638, hk_639, \
                         hk_648 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = -ab_x[k] * hi_500[k]
                       + hk_636[k];

            t_501[k] = -ab_x[k] * hi_501[k]
                       + hk_637[k];

            t_502[k] = -ab_x[k] * hi_502[k]
                       + hk_638[k];

            t_503[k] = -ab_x[k] * hi_503[k]
                       + hk_639[k];

            t_504[k] = -ab_x[k] * hi_504[k]
                       + hk_648[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_x, hi_505, hi_506, hi_507, \
                         hi_508, hi_509, hk_649, hk_650, hk_651, hk_652, \
                         hk_653 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = -ab_x[k] * hi_505[k]
                       + hk_649[k];

            t_506[k] = -ab_x[k] * hi_506[k]
                       + hk_650[k];

            t_507[k] = -ab_x[k] * hi_507[k]
                       + hk_651[k];

            t_508[k] = -ab_x[k] * hi_508[k]
                       + hk_652[k];

            t_509[k] = -ab_x[k] * hi_509[k]
                       + hk_653[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_x, hi_510, hi_511, hi_512, \
                         hi_513, hi_514, hk_654, hk_655, hk_656, hk_657, \
                         hk_658 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = -ab_x[k] * hi_510[k]
                       + hk_654[k];

            t_511[k] = -ab_x[k] * hi_511[k]
                       + hk_655[k];

            t_512[k] = -ab_x[k] * hi_512[k]
                       + hk_656[k];

            t_513[k] = -ab_x[k] * hi_513[k]
                       + hk_657[k];

            t_514[k] = -ab_x[k] * hi_514[k]
                       + hk_658[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_x, hi_515, hi_516, hi_517, \
                         hi_518, hi_519, hk_659, hk_660, hk_661, hk_662, \
                         hk_663 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = -ab_x[k] * hi_515[k]
                       + hk_659[k];

            t_516[k] = -ab_x[k] * hi_516[k]
                       + hk_660[k];

            t_517[k] = -ab_x[k] * hi_517[k]
                       + hk_661[k];

            t_518[k] = -ab_x[k] * hi_518[k]
                       + hk_662[k];

            t_519[k] = -ab_x[k] * hi_519[k]
                       + hk_663[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_x, hi_520, hi_521, hi_522, \
                         hi_523, hi_524, hk_664, hk_665, hk_666, hk_667, \
                         hk_668 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = -ab_x[k] * hi_520[k]
                       + hk_664[k];

            t_521[k] = -ab_x[k] * hi_521[k]
                       + hk_665[k];

            t_522[k] = -ab_x[k] * hi_522[k]
                       + hk_666[k];

            t_523[k] = -ab_x[k] * hi_523[k]
                       + hk_667[k];

            t_524[k] = -ab_x[k] * hi_524[k]
                       + hk_668[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_x, hi_525, hi_526, hi_527, \
                         hi_528, hi_529, hk_669, hk_670, hk_671, hk_672, \
                         hk_673 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = -ab_x[k] * hi_525[k]
                       + hk_669[k];

            t_526[k] = -ab_x[k] * hi_526[k]
                       + hk_670[k];

            t_527[k] = -ab_x[k] * hi_527[k]
                       + hk_671[k];

            t_528[k] = -ab_x[k] * hi_528[k]
                       + hk_672[k];

            t_529[k] = -ab_x[k] * hi_529[k]
                       + hk_673[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_x, hi_530, hi_531, hi_532, \
                         hi_533, hi_534, hk_674, hk_675, hk_684, hk_685, \
                         hk_686 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = -ab_x[k] * hi_530[k]
                       + hk_674[k];

            t_531[k] = -ab_x[k] * hi_531[k]
                       + hk_675[k];

            t_532[k] = -ab_x[k] * hi_532[k]
                       + hk_684[k];

            t_533[k] = -ab_x[k] * hi_533[k]
                       + hk_685[k];

            t_534[k] = -ab_x[k] * hi_534[k]
                       + hk_686[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_x, hi_535, hi_536, hi_537, \
                         hi_538, hi_539, hk_687, hk_688, hk_689, hk_690, \
                         hk_691 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = -ab_x[k] * hi_535[k]
                       + hk_687[k];

            t_536[k] = -ab_x[k] * hi_536[k]
                       + hk_688[k];

            t_537[k] = -ab_x[k] * hi_537[k]
                       + hk_689[k];

            t_538[k] = -ab_x[k] * hi_538[k]
                       + hk_690[k];

            t_539[k] = -ab_x[k] * hi_539[k]
                       + hk_691[k];
        }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_x, hi_540, hi_541, hi_542, \
                         hi_543, hi_544, hk_692, hk_693, hk_694, hk_695, \
                         hk_696 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_540[k] = -ab_x[k] * hi_540[k]
                       + hk_692[k];

            t_541[k] = -ab_x[k] * hi_541[k]
                       + hk_693[k];

            t_542[k] = -ab_x[k] * hi_542[k]
                       + hk_694[k];

            t_543[k] = -ab_x[k] * hi_543[k]
                       + hk_695[k];

            t_544[k] = -ab_x[k] * hi_544[k]
                       + hk_696[k];
        }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_x, hi_545, hi_546, hi_547, \
                         hi_548, hi_549, hk_697, hk_698, hk_699, hk_700, \
                         hk_701 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_545[k] = -ab_x[k] * hi_545[k]
                       + hk_697[k];

            t_546[k] = -ab_x[k] * hi_546[k]
                       + hk_698[k];

            t_547[k] = -ab_x[k] * hi_547[k]
                       + hk_699[k];

            t_548[k] = -ab_x[k] * hi_548[k]
                       + hk_700[k];

            t_549[k] = -ab_x[k] * hi_549[k]
                       + hk_701[k];
        }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ab_x, hi_550, hi_551, hi_552, \
                         hi_553, hi_554, hk_702, hk_703, hk_704, hk_705, \
                         hk_706 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_550[k] = -ab_x[k] * hi_550[k]
                       + hk_702[k];

            t_551[k] = -ab_x[k] * hi_551[k]
                       + hk_703[k];

            t_552[k] = -ab_x[k] * hi_552[k]
                       + hk_704[k];

            t_553[k] = -ab_x[k] * hi_553[k]
                       + hk_705[k];

            t_554[k] = -ab_x[k] * hi_554[k]
                       + hk_706[k];
        }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ab_x, hi_555, hi_556, hi_557, \
                         hi_558, hi_559, hk_707, hk_708, hk_709, hk_710, \
                         hk_711 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_555[k] = -ab_x[k] * hi_555[k]
                       + hk_707[k];

            t_556[k] = -ab_x[k] * hi_556[k]
                       + hk_708[k];

            t_557[k] = -ab_x[k] * hi_557[k]
                       + hk_709[k];

            t_558[k] = -ab_x[k] * hi_558[k]
                       + hk_710[k];

            t_559[k] = -ab_x[k] * hi_559[k]
                       + hk_711[k];
        }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_x, hi_560, hi_561, hi_562, \
                         hi_563, hi_564, hk_720, hk_721, hk_722, hk_723, \
                         hk_724 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_560[k] = -ab_x[k] * hi_560[k]
                       + hk_720[k];

            t_561[k] = -ab_x[k] * hi_561[k]
                       + hk_721[k];

            t_562[k] = -ab_x[k] * hi_562[k]
                       + hk_722[k];

            t_563[k] = -ab_x[k] * hi_563[k]
                       + hk_723[k];

            t_564[k] = -ab_x[k] * hi_564[k]
                       + hk_724[k];
        }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_x, hi_565, hi_566, hi_567, \
                         hi_568, hi_569, hk_725, hk_726, hk_727, hk_728, \
                         hk_729 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_565[k] = -ab_x[k] * hi_565[k]
                       + hk_725[k];

            t_566[k] = -ab_x[k] * hi_566[k]
                       + hk_726[k];

            t_567[k] = -ab_x[k] * hi_567[k]
                       + hk_727[k];

            t_568[k] = -ab_x[k] * hi_568[k]
                       + hk_728[k];

            t_569[k] = -ab_x[k] * hi_569[k]
                       + hk_729[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_x, hi_570, hi_571, hi_572, \
                         hi_573, hi_574, hk_730, hk_731, hk_732, hk_733, \
                         hk_734 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = -ab_x[k] * hi_570[k]
                       + hk_730[k];

            t_571[k] = -ab_x[k] * hi_571[k]
                       + hk_731[k];

            t_572[k] = -ab_x[k] * hi_572[k]
                       + hk_732[k];

            t_573[k] = -ab_x[k] * hi_573[k]
                       + hk_733[k];

            t_574[k] = -ab_x[k] * hi_574[k]
                       + hk_734[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_x, hi_575, hi_576, hi_577, \
                         hi_578, hi_579, hk_735, hk_736, hk_737, hk_738, \
                         hk_739 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = -ab_x[k] * hi_575[k]
                       + hk_735[k];

            t_576[k] = -ab_x[k] * hi_576[k]
                       + hk_736[k];

            t_577[k] = -ab_x[k] * hi_577[k]
                       + hk_737[k];

            t_578[k] = -ab_x[k] * hi_578[k]
                       + hk_738[k];

            t_579[k] = -ab_x[k] * hi_579[k]
                       + hk_739[k];
        }
    }
}

static auto
compute_hrr_ii_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hi, const size_t hk, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_580 = buffer.data(target + 580 * ncomps + c);
        auto *t_581 = buffer.data(target + 581 * ncomps + c);
        auto *t_582 = buffer.data(target + 582 * ncomps + c);
        auto *t_583 = buffer.data(target + 583 * ncomps + c);
        auto *t_584 = buffer.data(target + 584 * ncomps + c);
        auto *t_585 = buffer.data(target + 585 * ncomps + c);
        auto *t_586 = buffer.data(target + 586 * ncomps + c);
        auto *t_587 = buffer.data(target + 587 * ncomps + c);
        auto *t_588 = buffer.data(target + 588 * ncomps + c);
        auto *t_589 = buffer.data(target + 589 * ncomps + c);
        auto *t_590 = buffer.data(target + 590 * ncomps + c);
        auto *t_591 = buffer.data(target + 591 * ncomps + c);
        auto *t_592 = buffer.data(target + 592 * ncomps + c);
        auto *t_593 = buffer.data(target + 593 * ncomps + c);
        auto *t_594 = buffer.data(target + 594 * ncomps + c);
        auto *t_595 = buffer.data(target + 595 * ncomps + c);
        auto *t_596 = buffer.data(target + 596 * ncomps + c);
        auto *t_597 = buffer.data(target + 597 * ncomps + c);
        auto *t_598 = buffer.data(target + 598 * ncomps + c);
        auto *t_599 = buffer.data(target + 599 * ncomps + c);
        auto *t_600 = buffer.data(target + 600 * ncomps + c);
        auto *t_601 = buffer.data(target + 601 * ncomps + c);
        auto *t_602 = buffer.data(target + 602 * ncomps + c);
        auto *t_603 = buffer.data(target + 603 * ncomps + c);
        auto *t_604 = buffer.data(target + 604 * ncomps + c);
        auto *t_605 = buffer.data(target + 605 * ncomps + c);
        auto *t_606 = buffer.data(target + 606 * ncomps + c);
        auto *t_607 = buffer.data(target + 607 * ncomps + c);
        auto *t_608 = buffer.data(target + 608 * ncomps + c);
        auto *t_609 = buffer.data(target + 609 * ncomps + c);
        auto *t_610 = buffer.data(target + 610 * ncomps + c);
        auto *t_611 = buffer.data(target + 611 * ncomps + c);
        auto *t_612 = buffer.data(target + 612 * ncomps + c);
        auto *t_613 = buffer.data(target + 613 * ncomps + c);
        auto *t_614 = buffer.data(target + 614 * ncomps + c);
        auto *t_615 = buffer.data(target + 615 * ncomps + c);
        auto *t_616 = buffer.data(target + 616 * ncomps + c);
        auto *t_617 = buffer.data(target + 617 * ncomps + c);
        auto *t_618 = buffer.data(target + 618 * ncomps + c);
        auto *t_619 = buffer.data(target + 619 * ncomps + c);
        auto *t_620 = buffer.data(target + 620 * ncomps + c);
        auto *t_621 = buffer.data(target + 621 * ncomps + c);
        auto *t_622 = buffer.data(target + 622 * ncomps + c);
        auto *t_623 = buffer.data(target + 623 * ncomps + c);
        auto *t_624 = buffer.data(target + 624 * ncomps + c);
        auto *t_625 = buffer.data(target + 625 * ncomps + c);
        auto *t_626 = buffer.data(target + 626 * ncomps + c);
        auto *t_627 = buffer.data(target + 627 * ncomps + c);
        auto *t_628 = buffer.data(target + 628 * ncomps + c);
        auto *t_629 = buffer.data(target + 629 * ncomps + c);
        auto *t_630 = buffer.data(target + 630 * ncomps + c);
        auto *t_631 = buffer.data(target + 631 * ncomps + c);
        auto *t_632 = buffer.data(target + 632 * ncomps + c);
        auto *t_633 = buffer.data(target + 633 * ncomps + c);
        auto *t_634 = buffer.data(target + 634 * ncomps + c);
        auto *t_635 = buffer.data(target + 635 * ncomps + c);
        auto *t_636 = buffer.data(target + 636 * ncomps + c);
        auto *t_637 = buffer.data(target + 637 * ncomps + c);
        auto *t_638 = buffer.data(target + 638 * ncomps + c);
        auto *t_639 = buffer.data(target + 639 * ncomps + c);
        auto *t_640 = buffer.data(target + 640 * ncomps + c);
        auto *t_641 = buffer.data(target + 641 * ncomps + c);
        auto *t_642 = buffer.data(target + 642 * ncomps + c);
        auto *t_643 = buffer.data(target + 643 * ncomps + c);
        auto *t_644 = buffer.data(target + 644 * ncomps + c);
        auto *t_645 = buffer.data(target + 645 * ncomps + c);
        auto *t_646 = buffer.data(target + 646 * ncomps + c);
        auto *t_647 = buffer.data(target + 647 * ncomps + c);
        auto *t_648 = buffer.data(target + 648 * ncomps + c);
        auto *t_649 = buffer.data(target + 649 * ncomps + c);
        auto *t_650 = buffer.data(target + 650 * ncomps + c);
        auto *t_651 = buffer.data(target + 651 * ncomps + c);
        auto *t_652 = buffer.data(target + 652 * ncomps + c);
        auto *t_653 = buffer.data(target + 653 * ncomps + c);
        auto *t_654 = buffer.data(target + 654 * ncomps + c);
        auto *t_655 = buffer.data(target + 655 * ncomps + c);
        auto *t_656 = buffer.data(target + 656 * ncomps + c);
        auto *t_657 = buffer.data(target + 657 * ncomps + c);
        auto *t_658 = buffer.data(target + 658 * ncomps + c);
        auto *t_659 = buffer.data(target + 659 * ncomps + c);
        auto *t_660 = buffer.data(target + 660 * ncomps + c);
        auto *t_661 = buffer.data(target + 661 * ncomps + c);
        auto *t_662 = buffer.data(target + 662 * ncomps + c);
        auto *t_663 = buffer.data(target + 663 * ncomps + c);
        auto *t_664 = buffer.data(target + 664 * ncomps + c);
        auto *t_665 = buffer.data(target + 665 * ncomps + c);
        auto *t_666 = buffer.data(target + 666 * ncomps + c);
        auto *t_667 = buffer.data(target + 667 * ncomps + c);
        auto *t_668 = buffer.data(target + 668 * ncomps + c);
        auto *t_669 = buffer.data(target + 669 * ncomps + c);
        auto *t_670 = buffer.data(target + 670 * ncomps + c);
        auto *t_671 = buffer.data(target + 671 * ncomps + c);
        auto *t_672 = buffer.data(target + 672 * ncomps + c);
        auto *t_673 = buffer.data(target + 673 * ncomps + c);
        auto *t_674 = buffer.data(target + 674 * ncomps + c);
        auto *t_675 = buffer.data(target + 675 * ncomps + c);
        auto *t_676 = buffer.data(target + 676 * ncomps + c);
        auto *t_677 = buffer.data(target + 677 * ncomps + c);
        auto *t_678 = buffer.data(target + 678 * ncomps + c);
        auto *t_679 = buffer.data(target + 679 * ncomps + c);
        auto *t_680 = buffer.data(target + 680 * ncomps + c);
        auto *t_681 = buffer.data(target + 681 * ncomps + c);
        auto *t_682 = buffer.data(target + 682 * ncomps + c);
        auto *t_683 = buffer.data(target + 683 * ncomps + c);
        auto *t_684 = buffer.data(target + 684 * ncomps + c);
        auto *t_685 = buffer.data(target + 685 * ncomps + c);
        auto *t_686 = buffer.data(target + 686 * ncomps + c);
        auto *t_687 = buffer.data(target + 687 * ncomps + c);
        auto *t_688 = buffer.data(target + 688 * ncomps + c);
        auto *t_689 = buffer.data(target + 689 * ncomps + c);
        auto *t_690 = buffer.data(target + 690 * ncomps + c);
        auto *t_691 = buffer.data(target + 691 * ncomps + c);
        auto *t_692 = buffer.data(target + 692 * ncomps + c);
        auto *t_693 = buffer.data(target + 693 * ncomps + c);
        auto *t_694 = buffer.data(target + 694 * ncomps + c);
        auto *t_695 = buffer.data(target + 695 * ncomps + c);
        auto *t_696 = buffer.data(target + 696 * ncomps + c);
        auto *t_697 = buffer.data(target + 697 * ncomps + c);
        auto *t_698 = buffer.data(target + 698 * ncomps + c);
        auto *t_699 = buffer.data(target + 699 * ncomps + c);
        auto *t_700 = buffer.data(target + 700 * ncomps + c);
        auto *t_701 = buffer.data(target + 701 * ncomps + c);
        auto *t_702 = buffer.data(target + 702 * ncomps + c);
        auto *t_703 = buffer.data(target + 703 * ncomps + c);
        auto *t_704 = buffer.data(target + 704 * ncomps + c);
        auto *t_705 = buffer.data(target + 705 * ncomps + c);
        auto *t_706 = buffer.data(target + 706 * ncomps + c);
        auto *t_707 = buffer.data(target + 707 * ncomps + c);
        auto *t_708 = buffer.data(target + 708 * ncomps + c);
        auto *t_709 = buffer.data(target + 709 * ncomps + c);
        auto *t_710 = buffer.data(target + 710 * ncomps + c);
        auto *t_711 = buffer.data(target + 711 * ncomps + c);
        auto *t_712 = buffer.data(target + 712 * ncomps + c);
        auto *t_713 = buffer.data(target + 713 * ncomps + c);
        auto *t_714 = buffer.data(target + 714 * ncomps + c);
        auto *t_715 = buffer.data(target + 715 * ncomps + c);
        auto *t_716 = buffer.data(target + 716 * ncomps + c);
        auto *t_717 = buffer.data(target + 717 * ncomps + c);
        auto *t_718 = buffer.data(target + 718 * ncomps + c);
        auto *t_719 = buffer.data(target + 719 * ncomps + c);
        auto *t_720 = buffer.data(target + 720 * ncomps + c);
        auto *t_721 = buffer.data(target + 721 * ncomps + c);
        auto *t_722 = buffer.data(target + 722 * ncomps + c);
        auto *t_723 = buffer.data(target + 723 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);

        const auto *hi_420 = buffer.data(hi + 420 * ncomps + c);
        const auto *hi_421 = buffer.data(hi + 421 * ncomps + c);
        const auto *hi_422 = buffer.data(hi + 422 * ncomps + c);
        const auto *hi_423 = buffer.data(hi + 423 * ncomps + c);
        const auto *hi_424 = buffer.data(hi + 424 * ncomps + c);
        const auto *hi_425 = buffer.data(hi + 425 * ncomps + c);
        const auto *hi_426 = buffer.data(hi + 426 * ncomps + c);
        const auto *hi_427 = buffer.data(hi + 427 * ncomps + c);
        const auto *hi_428 = buffer.data(hi + 428 * ncomps + c);
        const auto *hi_429 = buffer.data(hi + 429 * ncomps + c);
        const auto *hi_430 = buffer.data(hi + 430 * ncomps + c);
        const auto *hi_431 = buffer.data(hi + 431 * ncomps + c);
        const auto *hi_432 = buffer.data(hi + 432 * ncomps + c);
        const auto *hi_433 = buffer.data(hi + 433 * ncomps + c);
        const auto *hi_434 = buffer.data(hi + 434 * ncomps + c);
        const auto *hi_435 = buffer.data(hi + 435 * ncomps + c);
        const auto *hi_436 = buffer.data(hi + 436 * ncomps + c);
        const auto *hi_437 = buffer.data(hi + 437 * ncomps + c);
        const auto *hi_438 = buffer.data(hi + 438 * ncomps + c);
        const auto *hi_439 = buffer.data(hi + 439 * ncomps + c);
        const auto *hi_440 = buffer.data(hi + 440 * ncomps + c);
        const auto *hi_441 = buffer.data(hi + 441 * ncomps + c);
        const auto *hi_442 = buffer.data(hi + 442 * ncomps + c);
        const auto *hi_443 = buffer.data(hi + 443 * ncomps + c);
        const auto *hi_444 = buffer.data(hi + 444 * ncomps + c);
        const auto *hi_445 = buffer.data(hi + 445 * ncomps + c);
        const auto *hi_446 = buffer.data(hi + 446 * ncomps + c);
        const auto *hi_447 = buffer.data(hi + 447 * ncomps + c);
        const auto *hi_448 = buffer.data(hi + 448 * ncomps + c);
        const auto *hi_449 = buffer.data(hi + 449 * ncomps + c);
        const auto *hi_450 = buffer.data(hi + 450 * ncomps + c);
        const auto *hi_451 = buffer.data(hi + 451 * ncomps + c);
        const auto *hi_452 = buffer.data(hi + 452 * ncomps + c);
        const auto *hi_453 = buffer.data(hi + 453 * ncomps + c);
        const auto *hi_454 = buffer.data(hi + 454 * ncomps + c);
        const auto *hi_455 = buffer.data(hi + 455 * ncomps + c);
        const auto *hi_456 = buffer.data(hi + 456 * ncomps + c);
        const auto *hi_457 = buffer.data(hi + 457 * ncomps + c);
        const auto *hi_458 = buffer.data(hi + 458 * ncomps + c);
        const auto *hi_459 = buffer.data(hi + 459 * ncomps + c);
        const auto *hi_460 = buffer.data(hi + 460 * ncomps + c);
        const auto *hi_461 = buffer.data(hi + 461 * ncomps + c);
        const auto *hi_462 = buffer.data(hi + 462 * ncomps + c);
        const auto *hi_463 = buffer.data(hi + 463 * ncomps + c);
        const auto *hi_464 = buffer.data(hi + 464 * ncomps + c);
        const auto *hi_465 = buffer.data(hi + 465 * ncomps + c);
        const auto *hi_466 = buffer.data(hi + 466 * ncomps + c);
        const auto *hi_467 = buffer.data(hi + 467 * ncomps + c);
        const auto *hi_468 = buffer.data(hi + 468 * ncomps + c);
        const auto *hi_469 = buffer.data(hi + 469 * ncomps + c);
        const auto *hi_470 = buffer.data(hi + 470 * ncomps + c);
        const auto *hi_471 = buffer.data(hi + 471 * ncomps + c);
        const auto *hi_472 = buffer.data(hi + 472 * ncomps + c);
        const auto *hi_473 = buffer.data(hi + 473 * ncomps + c);
        const auto *hi_474 = buffer.data(hi + 474 * ncomps + c);
        const auto *hi_475 = buffer.data(hi + 475 * ncomps + c);
        const auto *hi_476 = buffer.data(hi + 476 * ncomps + c);
        const auto *hi_477 = buffer.data(hi + 477 * ncomps + c);
        const auto *hi_478 = buffer.data(hi + 478 * ncomps + c);
        const auto *hi_479 = buffer.data(hi + 479 * ncomps + c);
        const auto *hi_480 = buffer.data(hi + 480 * ncomps + c);
        const auto *hi_481 = buffer.data(hi + 481 * ncomps + c);
        const auto *hi_482 = buffer.data(hi + 482 * ncomps + c);
        const auto *hi_483 = buffer.data(hi + 483 * ncomps + c);
        const auto *hi_484 = buffer.data(hi + 484 * ncomps + c);
        const auto *hi_485 = buffer.data(hi + 485 * ncomps + c);
        const auto *hi_486 = buffer.data(hi + 486 * ncomps + c);
        const auto *hi_487 = buffer.data(hi + 487 * ncomps + c);
        const auto *hi_488 = buffer.data(hi + 488 * ncomps + c);
        const auto *hi_489 = buffer.data(hi + 489 * ncomps + c);
        const auto *hi_490 = buffer.data(hi + 490 * ncomps + c);
        const auto *hi_491 = buffer.data(hi + 491 * ncomps + c);
        const auto *hi_492 = buffer.data(hi + 492 * ncomps + c);
        const auto *hi_493 = buffer.data(hi + 493 * ncomps + c);
        const auto *hi_494 = buffer.data(hi + 494 * ncomps + c);
        const auto *hi_495 = buffer.data(hi + 495 * ncomps + c);
        const auto *hi_496 = buffer.data(hi + 496 * ncomps + c);
        const auto *hi_497 = buffer.data(hi + 497 * ncomps + c);
        const auto *hi_498 = buffer.data(hi + 498 * ncomps + c);
        const auto *hi_499 = buffer.data(hi + 499 * ncomps + c);
        const auto *hi_500 = buffer.data(hi + 500 * ncomps + c);
        const auto *hi_501 = buffer.data(hi + 501 * ncomps + c);
        const auto *hi_502 = buffer.data(hi + 502 * ncomps + c);
        const auto *hi_503 = buffer.data(hi + 503 * ncomps + c);
        const auto *hi_504 = buffer.data(hi + 504 * ncomps + c);
        const auto *hi_505 = buffer.data(hi + 505 * ncomps + c);
        const auto *hi_506 = buffer.data(hi + 506 * ncomps + c);
        const auto *hi_507 = buffer.data(hi + 507 * ncomps + c);
        const auto *hi_508 = buffer.data(hi + 508 * ncomps + c);
        const auto *hi_509 = buffer.data(hi + 509 * ncomps + c);
        const auto *hi_510 = buffer.data(hi + 510 * ncomps + c);
        const auto *hi_511 = buffer.data(hi + 511 * ncomps + c);
        const auto *hi_512 = buffer.data(hi + 512 * ncomps + c);
        const auto *hi_513 = buffer.data(hi + 513 * ncomps + c);
        const auto *hi_514 = buffer.data(hi + 514 * ncomps + c);
        const auto *hi_515 = buffer.data(hi + 515 * ncomps + c);
        const auto *hi_516 = buffer.data(hi + 516 * ncomps + c);
        const auto *hi_517 = buffer.data(hi + 517 * ncomps + c);
        const auto *hi_518 = buffer.data(hi + 518 * ncomps + c);
        const auto *hi_519 = buffer.data(hi + 519 * ncomps + c);
        const auto *hi_520 = buffer.data(hi + 520 * ncomps + c);
        const auto *hi_521 = buffer.data(hi + 521 * ncomps + c);
        const auto *hi_522 = buffer.data(hi + 522 * ncomps + c);
        const auto *hi_523 = buffer.data(hi + 523 * ncomps + c);
        const auto *hi_524 = buffer.data(hi + 524 * ncomps + c);
        const auto *hi_525 = buffer.data(hi + 525 * ncomps + c);
        const auto *hi_526 = buffer.data(hi + 526 * ncomps + c);
        const auto *hi_527 = buffer.data(hi + 527 * ncomps + c);
        const auto *hi_528 = buffer.data(hi + 528 * ncomps + c);
        const auto *hi_529 = buffer.data(hi + 529 * ncomps + c);
        const auto *hi_530 = buffer.data(hi + 530 * ncomps + c);
        const auto *hi_531 = buffer.data(hi + 531 * ncomps + c);
        const auto *hi_532 = buffer.data(hi + 532 * ncomps + c);
        const auto *hi_533 = buffer.data(hi + 533 * ncomps + c);
        const auto *hi_534 = buffer.data(hi + 534 * ncomps + c);
        const auto *hi_535 = buffer.data(hi + 535 * ncomps + c);
        const auto *hi_536 = buffer.data(hi + 536 * ncomps + c);
        const auto *hi_537 = buffer.data(hi + 537 * ncomps + c);
        const auto *hi_538 = buffer.data(hi + 538 * ncomps + c);
        const auto *hi_539 = buffer.data(hi + 539 * ncomps + c);
        const auto *hi_540 = buffer.data(hi + 540 * ncomps + c);
        const auto *hi_541 = buffer.data(hi + 541 * ncomps + c);
        const auto *hi_542 = buffer.data(hi + 542 * ncomps + c);
        const auto *hi_543 = buffer.data(hi + 543 * ncomps + c);
        const auto *hi_544 = buffer.data(hi + 544 * ncomps + c);
        const auto *hi_545 = buffer.data(hi + 545 * ncomps + c);
        const auto *hi_546 = buffer.data(hi + 546 * ncomps + c);
        const auto *hi_547 = buffer.data(hi + 547 * ncomps + c);
        const auto *hi_548 = buffer.data(hi + 548 * ncomps + c);
        const auto *hi_549 = buffer.data(hi + 549 * ncomps + c);
        const auto *hi_550 = buffer.data(hi + 550 * ncomps + c);
        const auto *hi_551 = buffer.data(hi + 551 * ncomps + c);
        const auto *hi_552 = buffer.data(hi + 552 * ncomps + c);
        const auto *hi_553 = buffer.data(hi + 553 * ncomps + c);
        const auto *hi_554 = buffer.data(hi + 554 * ncomps + c);
        const auto *hi_555 = buffer.data(hi + 555 * ncomps + c);
        const auto *hi_580 = buffer.data(hi + 580 * ncomps + c);
        const auto *hi_581 = buffer.data(hi + 581 * ncomps + c);
        const auto *hi_582 = buffer.data(hi + 582 * ncomps + c);
        const auto *hi_583 = buffer.data(hi + 583 * ncomps + c);
        const auto *hi_584 = buffer.data(hi + 584 * ncomps + c);
        const auto *hi_585 = buffer.data(hi + 585 * ncomps + c);
        const auto *hi_586 = buffer.data(hi + 586 * ncomps + c);
        const auto *hi_587 = buffer.data(hi + 587 * ncomps + c);

        const auto *hk_541 = buffer.data(hk + 541 * ncomps + c);
        const auto *hk_543 = buffer.data(hk + 543 * ncomps + c);
        const auto *hk_544 = buffer.data(hk + 544 * ncomps + c);
        const auto *hk_546 = buffer.data(hk + 546 * ncomps + c);
        const auto *hk_547 = buffer.data(hk + 547 * ncomps + c);
        const auto *hk_548 = buffer.data(hk + 548 * ncomps + c);
        const auto *hk_550 = buffer.data(hk + 550 * ncomps + c);
        const auto *hk_551 = buffer.data(hk + 551 * ncomps + c);
        const auto *hk_552 = buffer.data(hk + 552 * ncomps + c);
        const auto *hk_553 = buffer.data(hk + 553 * ncomps + c);
        const auto *hk_555 = buffer.data(hk + 555 * ncomps + c);
        const auto *hk_556 = buffer.data(hk + 556 * ncomps + c);
        const auto *hk_557 = buffer.data(hk + 557 * ncomps + c);
        const auto *hk_558 = buffer.data(hk + 558 * ncomps + c);
        const auto *hk_559 = buffer.data(hk + 559 * ncomps + c);
        const auto *hk_561 = buffer.data(hk + 561 * ncomps + c);
        const auto *hk_562 = buffer.data(hk + 562 * ncomps + c);
        const auto *hk_563 = buffer.data(hk + 563 * ncomps + c);
        const auto *hk_564 = buffer.data(hk + 564 * ncomps + c);
        const auto *hk_565 = buffer.data(hk + 565 * ncomps + c);
        const auto *hk_566 = buffer.data(hk + 566 * ncomps + c);
        const auto *hk_568 = buffer.data(hk + 568 * ncomps + c);
        const auto *hk_569 = buffer.data(hk + 569 * ncomps + c);
        const auto *hk_570 = buffer.data(hk + 570 * ncomps + c);
        const auto *hk_571 = buffer.data(hk + 571 * ncomps + c);
        const auto *hk_572 = buffer.data(hk + 572 * ncomps + c);
        const auto *hk_573 = buffer.data(hk + 573 * ncomps + c);
        const auto *hk_574 = buffer.data(hk + 574 * ncomps + c);
        const auto *hk_577 = buffer.data(hk + 577 * ncomps + c);
        const auto *hk_579 = buffer.data(hk + 579 * ncomps + c);
        const auto *hk_580 = buffer.data(hk + 580 * ncomps + c);
        const auto *hk_582 = buffer.data(hk + 582 * ncomps + c);
        const auto *hk_583 = buffer.data(hk + 583 * ncomps + c);
        const auto *hk_584 = buffer.data(hk + 584 * ncomps + c);
        const auto *hk_586 = buffer.data(hk + 586 * ncomps + c);
        const auto *hk_587 = buffer.data(hk + 587 * ncomps + c);
        const auto *hk_588 = buffer.data(hk + 588 * ncomps + c);
        const auto *hk_589 = buffer.data(hk + 589 * ncomps + c);
        const auto *hk_591 = buffer.data(hk + 591 * ncomps + c);
        const auto *hk_592 = buffer.data(hk + 592 * ncomps + c);
        const auto *hk_593 = buffer.data(hk + 593 * ncomps + c);
        const auto *hk_594 = buffer.data(hk + 594 * ncomps + c);
        const auto *hk_595 = buffer.data(hk + 595 * ncomps + c);
        const auto *hk_597 = buffer.data(hk + 597 * ncomps + c);
        const auto *hk_598 = buffer.data(hk + 598 * ncomps + c);
        const auto *hk_599 = buffer.data(hk + 599 * ncomps + c);
        const auto *hk_600 = buffer.data(hk + 600 * ncomps + c);
        const auto *hk_601 = buffer.data(hk + 601 * ncomps + c);
        const auto *hk_602 = buffer.data(hk + 602 * ncomps + c);
        const auto *hk_604 = buffer.data(hk + 604 * ncomps + c);
        const auto *hk_605 = buffer.data(hk + 605 * ncomps + c);
        const auto *hk_606 = buffer.data(hk + 606 * ncomps + c);
        const auto *hk_607 = buffer.data(hk + 607 * ncomps + c);
        const auto *hk_608 = buffer.data(hk + 608 * ncomps + c);
        const auto *hk_609 = buffer.data(hk + 609 * ncomps + c);
        const auto *hk_610 = buffer.data(hk + 610 * ncomps + c);
        const auto *hk_613 = buffer.data(hk + 613 * ncomps + c);
        const auto *hk_615 = buffer.data(hk + 615 * ncomps + c);
        const auto *hk_616 = buffer.data(hk + 616 * ncomps + c);
        const auto *hk_618 = buffer.data(hk + 618 * ncomps + c);
        const auto *hk_619 = buffer.data(hk + 619 * ncomps + c);
        const auto *hk_620 = buffer.data(hk + 620 * ncomps + c);
        const auto *hk_622 = buffer.data(hk + 622 * ncomps + c);
        const auto *hk_623 = buffer.data(hk + 623 * ncomps + c);
        const auto *hk_624 = buffer.data(hk + 624 * ncomps + c);
        const auto *hk_625 = buffer.data(hk + 625 * ncomps + c);
        const auto *hk_627 = buffer.data(hk + 627 * ncomps + c);
        const auto *hk_628 = buffer.data(hk + 628 * ncomps + c);
        const auto *hk_629 = buffer.data(hk + 629 * ncomps + c);
        const auto *hk_630 = buffer.data(hk + 630 * ncomps + c);
        const auto *hk_631 = buffer.data(hk + 631 * ncomps + c);
        const auto *hk_633 = buffer.data(hk + 633 * ncomps + c);
        const auto *hk_634 = buffer.data(hk + 634 * ncomps + c);
        const auto *hk_635 = buffer.data(hk + 635 * ncomps + c);
        const auto *hk_636 = buffer.data(hk + 636 * ncomps + c);
        const auto *hk_637 = buffer.data(hk + 637 * ncomps + c);
        const auto *hk_638 = buffer.data(hk + 638 * ncomps + c);
        const auto *hk_640 = buffer.data(hk + 640 * ncomps + c);
        const auto *hk_641 = buffer.data(hk + 641 * ncomps + c);
        const auto *hk_642 = buffer.data(hk + 642 * ncomps + c);
        const auto *hk_643 = buffer.data(hk + 643 * ncomps + c);
        const auto *hk_644 = buffer.data(hk + 644 * ncomps + c);
        const auto *hk_645 = buffer.data(hk + 645 * ncomps + c);
        const auto *hk_646 = buffer.data(hk + 646 * ncomps + c);
        const auto *hk_649 = buffer.data(hk + 649 * ncomps + c);
        const auto *hk_651 = buffer.data(hk + 651 * ncomps + c);
        const auto *hk_652 = buffer.data(hk + 652 * ncomps + c);
        const auto *hk_654 = buffer.data(hk + 654 * ncomps + c);
        const auto *hk_655 = buffer.data(hk + 655 * ncomps + c);
        const auto *hk_656 = buffer.data(hk + 656 * ncomps + c);
        const auto *hk_658 = buffer.data(hk + 658 * ncomps + c);
        const auto *hk_659 = buffer.data(hk + 659 * ncomps + c);
        const auto *hk_660 = buffer.data(hk + 660 * ncomps + c);
        const auto *hk_661 = buffer.data(hk + 661 * ncomps + c);
        const auto *hk_663 = buffer.data(hk + 663 * ncomps + c);
        const auto *hk_664 = buffer.data(hk + 664 * ncomps + c);
        const auto *hk_665 = buffer.data(hk + 665 * ncomps + c);
        const auto *hk_666 = buffer.data(hk + 666 * ncomps + c);
        const auto *hk_667 = buffer.data(hk + 667 * ncomps + c);
        const auto *hk_669 = buffer.data(hk + 669 * ncomps + c);
        const auto *hk_670 = buffer.data(hk + 670 * ncomps + c);
        const auto *hk_671 = buffer.data(hk + 671 * ncomps + c);
        const auto *hk_672 = buffer.data(hk + 672 * ncomps + c);
        const auto *hk_673 = buffer.data(hk + 673 * ncomps + c);
        const auto *hk_674 = buffer.data(hk + 674 * ncomps + c);
        const auto *hk_676 = buffer.data(hk + 676 * ncomps + c);
        const auto *hk_677 = buffer.data(hk + 677 * ncomps + c);
        const auto *hk_678 = buffer.data(hk + 678 * ncomps + c);
        const auto *hk_679 = buffer.data(hk + 679 * ncomps + c);
        const auto *hk_680 = buffer.data(hk + 680 * ncomps + c);
        const auto *hk_681 = buffer.data(hk + 681 * ncomps + c);
        const auto *hk_682 = buffer.data(hk + 682 * ncomps + c);
        const auto *hk_685 = buffer.data(hk + 685 * ncomps + c);
        const auto *hk_687 = buffer.data(hk + 687 * ncomps + c);
        const auto *hk_688 = buffer.data(hk + 688 * ncomps + c);
        const auto *hk_690 = buffer.data(hk + 690 * ncomps + c);
        const auto *hk_691 = buffer.data(hk + 691 * ncomps + c);
        const auto *hk_692 = buffer.data(hk + 692 * ncomps + c);
        const auto *hk_694 = buffer.data(hk + 694 * ncomps + c);
        const auto *hk_695 = buffer.data(hk + 695 * ncomps + c);
        const auto *hk_696 = buffer.data(hk + 696 * ncomps + c);
        const auto *hk_697 = buffer.data(hk + 697 * ncomps + c);
        const auto *hk_699 = buffer.data(hk + 699 * ncomps + c);
        const auto *hk_700 = buffer.data(hk + 700 * ncomps + c);
        const auto *hk_701 = buffer.data(hk + 701 * ncomps + c);
        const auto *hk_702 = buffer.data(hk + 702 * ncomps + c);
        const auto *hk_703 = buffer.data(hk + 703 * ncomps + c);
        const auto *hk_705 = buffer.data(hk + 705 * ncomps + c);
        const auto *hk_706 = buffer.data(hk + 706 * ncomps + c);
        const auto *hk_707 = buffer.data(hk + 707 * ncomps + c);
        const auto *hk_708 = buffer.data(hk + 708 * ncomps + c);
        const auto *hk_709 = buffer.data(hk + 709 * ncomps + c);
        const auto *hk_710 = buffer.data(hk + 710 * ncomps + c);
        const auto *hk_712 = buffer.data(hk + 712 * ncomps + c);
        const auto *hk_713 = buffer.data(hk + 713 * ncomps + c);
        const auto *hk_714 = buffer.data(hk + 714 * ncomps + c);
        const auto *hk_740 = buffer.data(hk + 740 * ncomps + c);
        const auto *hk_741 = buffer.data(hk + 741 * ncomps + c);
        const auto *hk_742 = buffer.data(hk + 742 * ncomps + c);
        const auto *hk_743 = buffer.data(hk + 743 * ncomps + c);
        const auto *hk_744 = buffer.data(hk + 744 * ncomps + c);
        const auto *hk_745 = buffer.data(hk + 745 * ncomps + c);
        const auto *hk_746 = buffer.data(hk + 746 * ncomps + c);
        const auto *hk_747 = buffer.data(hk + 747 * ncomps + c);

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ab_x, hi_580, hi_581, hi_582, \
                         hi_583, hi_584, hk_740, hk_741, hk_742, hk_743, \
                         hk_744 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = -ab_x[k] * hi_580[k]
                       + hk_740[k];

            t_581[k] = -ab_x[k] * hi_581[k]
                       + hk_741[k];

            t_582[k] = -ab_x[k] * hi_582[k]
                       + hk_742[k];

            t_583[k] = -ab_x[k] * hi_583[k]
                       + hk_743[k];

            t_584[k] = -ab_x[k] * hi_584[k]
                       + hk_744[k];
        }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, ab_x, ab_y, hi_420, hi_585, hi_586, \
                         hi_587, hk_541, hk_745, hk_746, hk_747 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_585[k] = -ab_x[k] * hi_585[k]
                       + hk_745[k];

            t_586[k] = -ab_x[k] * hi_586[k]
                       + hk_746[k];

            t_587[k] = -ab_x[k] * hi_587[k]
                       + hk_747[k];

            t_588[k] = -ab_y[k] * hi_420[k]
                       + hk_541[k];
        }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, ab_y, hi_421, hi_422, hi_423, \
                         hi_424, hi_425, hk_543, hk_544, hk_546, hk_547, \
                         hk_548 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_589[k] = -ab_y[k] * hi_421[k]
                       + hk_543[k];

            t_590[k] = -ab_y[k] * hi_422[k]
                       + hk_544[k];

            t_591[k] = -ab_y[k] * hi_423[k]
                       + hk_546[k];

            t_592[k] = -ab_y[k] * hi_424[k]
                       + hk_547[k];

            t_593[k] = -ab_y[k] * hi_425[k]
                       + hk_548[k];
        }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, ab_y, hi_426, hi_427, hi_428, \
                         hi_429, hi_430, hk_550, hk_551, hk_552, hk_553, \
                         hk_555 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_594[k] = -ab_y[k] * hi_426[k]
                       + hk_550[k];

            t_595[k] = -ab_y[k] * hi_427[k]
                       + hk_551[k];

            t_596[k] = -ab_y[k] * hi_428[k]
                       + hk_552[k];

            t_597[k] = -ab_y[k] * hi_429[k]
                       + hk_553[k];

            t_598[k] = -ab_y[k] * hi_430[k]
                       + hk_555[k];
        }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, ab_y, hi_431, hi_432, hi_433, \
                         hi_434, hi_435, hk_556, hk_557, hk_558, hk_559, \
                         hk_561 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_599[k] = -ab_y[k] * hi_431[k]
                       + hk_556[k];

            t_600[k] = -ab_y[k] * hi_432[k]
                       + hk_557[k];

            t_601[k] = -ab_y[k] * hi_433[k]
                       + hk_558[k];

            t_602[k] = -ab_y[k] * hi_434[k]
                       + hk_559[k];

            t_603[k] = -ab_y[k] * hi_435[k]
                       + hk_561[k];
        }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, ab_y, hi_436, hi_437, hi_438, \
                         hi_439, hi_440, hk_562, hk_563, hk_564, hk_565, \
                         hk_566 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_604[k] = -ab_y[k] * hi_436[k]
                       + hk_562[k];

            t_605[k] = -ab_y[k] * hi_437[k]
                       + hk_563[k];

            t_606[k] = -ab_y[k] * hi_438[k]
                       + hk_564[k];

            t_607[k] = -ab_y[k] * hi_439[k]
                       + hk_565[k];

            t_608[k] = -ab_y[k] * hi_440[k]
                       + hk_566[k];
        }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, ab_y, hi_441, hi_442, hi_443, \
                         hi_444, hi_445, hk_568, hk_569, hk_570, hk_571, \
                         hk_572 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_609[k] = -ab_y[k] * hi_441[k]
                       + hk_568[k];

            t_610[k] = -ab_y[k] * hi_442[k]
                       + hk_569[k];

            t_611[k] = -ab_y[k] * hi_443[k]
                       + hk_570[k];

            t_612[k] = -ab_y[k] * hi_444[k]
                       + hk_571[k];

            t_613[k] = -ab_y[k] * hi_445[k]
                       + hk_572[k];
        }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, ab_y, hi_446, hi_447, hi_448, \
                         hi_449, hi_450, hk_573, hk_574, hk_577, hk_579, \
                         hk_580 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_614[k] = -ab_y[k] * hi_446[k]
                       + hk_573[k];

            t_615[k] = -ab_y[k] * hi_447[k]
                       + hk_574[k];

            t_616[k] = -ab_y[k] * hi_448[k]
                       + hk_577[k];

            t_617[k] = -ab_y[k] * hi_449[k]
                       + hk_579[k];

            t_618[k] = -ab_y[k] * hi_450[k]
                       + hk_580[k];
        }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, ab_y, hi_451, hi_452, hi_453, \
                         hi_454, hi_455, hk_582, hk_583, hk_584, hk_586, \
                         hk_587 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_619[k] = -ab_y[k] * hi_451[k]
                       + hk_582[k];

            t_620[k] = -ab_y[k] * hi_452[k]
                       + hk_583[k];

            t_621[k] = -ab_y[k] * hi_453[k]
                       + hk_584[k];

            t_622[k] = -ab_y[k] * hi_454[k]
                       + hk_586[k];

            t_623[k] = -ab_y[k] * hi_455[k]
                       + hk_587[k];
        }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, ab_y, hi_456, hi_457, hi_458, \
                         hi_459, hi_460, hk_588, hk_589, hk_591, hk_592, \
                         hk_593 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_624[k] = -ab_y[k] * hi_456[k]
                       + hk_588[k];

            t_625[k] = -ab_y[k] * hi_457[k]
                       + hk_589[k];

            t_626[k] = -ab_y[k] * hi_458[k]
                       + hk_591[k];

            t_627[k] = -ab_y[k] * hi_459[k]
                       + hk_592[k];

            t_628[k] = -ab_y[k] * hi_460[k]
                       + hk_593[k];
        }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, ab_y, hi_461, hi_462, hi_463, \
                         hi_464, hi_465, hk_594, hk_595, hk_597, hk_598, \
                         hk_599 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_629[k] = -ab_y[k] * hi_461[k]
                       + hk_594[k];

            t_630[k] = -ab_y[k] * hi_462[k]
                       + hk_595[k];

            t_631[k] = -ab_y[k] * hi_463[k]
                       + hk_597[k];

            t_632[k] = -ab_y[k] * hi_464[k]
                       + hk_598[k];

            t_633[k] = -ab_y[k] * hi_465[k]
                       + hk_599[k];
        }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, ab_y, hi_466, hi_467, hi_468, \
                         hi_469, hi_470, hk_600, hk_601, hk_602, hk_604, \
                         hk_605 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_634[k] = -ab_y[k] * hi_466[k]
                       + hk_600[k];

            t_635[k] = -ab_y[k] * hi_467[k]
                       + hk_601[k];

            t_636[k] = -ab_y[k] * hi_468[k]
                       + hk_602[k];

            t_637[k] = -ab_y[k] * hi_469[k]
                       + hk_604[k];

            t_638[k] = -ab_y[k] * hi_470[k]
                       + hk_605[k];
        }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, ab_y, hi_471, hi_472, hi_473, \
                         hi_474, hi_475, hk_606, hk_607, hk_608, hk_609, \
                         hk_610 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_639[k] = -ab_y[k] * hi_471[k]
                       + hk_606[k];

            t_640[k] = -ab_y[k] * hi_472[k]
                       + hk_607[k];

            t_641[k] = -ab_y[k] * hi_473[k]
                       + hk_608[k];

            t_642[k] = -ab_y[k] * hi_474[k]
                       + hk_609[k];

            t_643[k] = -ab_y[k] * hi_475[k]
                       + hk_610[k];
        }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, ab_y, hi_476, hi_477, hi_478, \
                         hi_479, hi_480, hk_613, hk_615, hk_616, hk_618, \
                         hk_619 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_644[k] = -ab_y[k] * hi_476[k]
                       + hk_613[k];

            t_645[k] = -ab_y[k] * hi_477[k]
                       + hk_615[k];

            t_646[k] = -ab_y[k] * hi_478[k]
                       + hk_616[k];

            t_647[k] = -ab_y[k] * hi_479[k]
                       + hk_618[k];

            t_648[k] = -ab_y[k] * hi_480[k]
                       + hk_619[k];
        }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, ab_y, hi_481, hi_482, hi_483, \
                         hi_484, hi_485, hk_620, hk_622, hk_623, hk_624, \
                         hk_625 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_649[k] = -ab_y[k] * hi_481[k]
                       + hk_620[k];

            t_650[k] = -ab_y[k] * hi_482[k]
                       + hk_622[k];

            t_651[k] = -ab_y[k] * hi_483[k]
                       + hk_623[k];

            t_652[k] = -ab_y[k] * hi_484[k]
                       + hk_624[k];

            t_653[k] = -ab_y[k] * hi_485[k]
                       + hk_625[k];
        }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, ab_y, hi_486, hi_487, hi_488, \
                         hi_489, hi_490, hk_627, hk_628, hk_629, hk_630, \
                         hk_631 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_654[k] = -ab_y[k] * hi_486[k]
                       + hk_627[k];

            t_655[k] = -ab_y[k] * hi_487[k]
                       + hk_628[k];

            t_656[k] = -ab_y[k] * hi_488[k]
                       + hk_629[k];

            t_657[k] = -ab_y[k] * hi_489[k]
                       + hk_630[k];

            t_658[k] = -ab_y[k] * hi_490[k]
                       + hk_631[k];
        }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, ab_y, hi_491, hi_492, hi_493, \
                         hi_494, hi_495, hk_633, hk_634, hk_635, hk_636, \
                         hk_637 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_659[k] = -ab_y[k] * hi_491[k]
                       + hk_633[k];

            t_660[k] = -ab_y[k] * hi_492[k]
                       + hk_634[k];

            t_661[k] = -ab_y[k] * hi_493[k]
                       + hk_635[k];

            t_662[k] = -ab_y[k] * hi_494[k]
                       + hk_636[k];

            t_663[k] = -ab_y[k] * hi_495[k]
                       + hk_637[k];
        }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, ab_y, hi_496, hi_497, hi_498, \
                         hi_499, hi_500, hk_638, hk_640, hk_641, hk_642, \
                         hk_643 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_664[k] = -ab_y[k] * hi_496[k]
                       + hk_638[k];

            t_665[k] = -ab_y[k] * hi_497[k]
                       + hk_640[k];

            t_666[k] = -ab_y[k] * hi_498[k]
                       + hk_641[k];

            t_667[k] = -ab_y[k] * hi_499[k]
                       + hk_642[k];

            t_668[k] = -ab_y[k] * hi_500[k]
                       + hk_643[k];
        }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, ab_y, hi_501, hi_502, hi_503, \
                         hi_504, hi_505, hk_644, hk_645, hk_646, hk_649, \
                         hk_651 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_669[k] = -ab_y[k] * hi_501[k]
                       + hk_644[k];

            t_670[k] = -ab_y[k] * hi_502[k]
                       + hk_645[k];

            t_671[k] = -ab_y[k] * hi_503[k]
                       + hk_646[k];

            t_672[k] = -ab_y[k] * hi_504[k]
                       + hk_649[k];

            t_673[k] = -ab_y[k] * hi_505[k]
                       + hk_651[k];
        }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, t_678, ab_y, hi_506, hi_507, hi_508, \
                         hi_509, hi_510, hk_652, hk_654, hk_655, hk_656, \
                         hk_658 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_674[k] = -ab_y[k] * hi_506[k]
                       + hk_652[k];

            t_675[k] = -ab_y[k] * hi_507[k]
                       + hk_654[k];

            t_676[k] = -ab_y[k] * hi_508[k]
                       + hk_655[k];

            t_677[k] = -ab_y[k] * hi_509[k]
                       + hk_656[k];

            t_678[k] = -ab_y[k] * hi_510[k]
                       + hk_658[k];
        }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, ab_y, hi_511, hi_512, hi_513, \
                         hi_514, hi_515, hk_659, hk_660, hk_661, hk_663, \
                         hk_664 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_679[k] = -ab_y[k] * hi_511[k]
                       + hk_659[k];

            t_680[k] = -ab_y[k] * hi_512[k]
                       + hk_660[k];

            t_681[k] = -ab_y[k] * hi_513[k]
                       + hk_661[k];

            t_682[k] = -ab_y[k] * hi_514[k]
                       + hk_663[k];

            t_683[k] = -ab_y[k] * hi_515[k]
                       + hk_664[k];
        }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, ab_y, hi_516, hi_517, hi_518, \
                         hi_519, hi_520, hk_665, hk_666, hk_667, hk_669, \
                         hk_670 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_684[k] = -ab_y[k] * hi_516[k]
                       + hk_665[k];

            t_685[k] = -ab_y[k] * hi_517[k]
                       + hk_666[k];

            t_686[k] = -ab_y[k] * hi_518[k]
                       + hk_667[k];

            t_687[k] = -ab_y[k] * hi_519[k]
                       + hk_669[k];

            t_688[k] = -ab_y[k] * hi_520[k]
                       + hk_670[k];
        }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, ab_y, hi_521, hi_522, hi_523, \
                         hi_524, hi_525, hk_671, hk_672, hk_673, hk_674, \
                         hk_676 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_689[k] = -ab_y[k] * hi_521[k]
                       + hk_671[k];

            t_690[k] = -ab_y[k] * hi_522[k]
                       + hk_672[k];

            t_691[k] = -ab_y[k] * hi_523[k]
                       + hk_673[k];

            t_692[k] = -ab_y[k] * hi_524[k]
                       + hk_674[k];

            t_693[k] = -ab_y[k] * hi_525[k]
                       + hk_676[k];
        }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, t_698, ab_y, hi_526, hi_527, hi_528, \
                         hi_529, hi_530, hk_677, hk_678, hk_679, hk_680, \
                         hk_681 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_694[k] = -ab_y[k] * hi_526[k]
                       + hk_677[k];

            t_695[k] = -ab_y[k] * hi_527[k]
                       + hk_678[k];

            t_696[k] = -ab_y[k] * hi_528[k]
                       + hk_679[k];

            t_697[k] = -ab_y[k] * hi_529[k]
                       + hk_680[k];

            t_698[k] = -ab_y[k] * hi_530[k]
                       + hk_681[k];
        }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, t_703, ab_y, hi_531, hi_532, hi_533, \
                         hi_534, hi_535, hk_682, hk_685, hk_687, hk_688, \
                         hk_690 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_699[k] = -ab_y[k] * hi_531[k]
                       + hk_682[k];

            t_700[k] = -ab_y[k] * hi_532[k]
                       + hk_685[k];

            t_701[k] = -ab_y[k] * hi_533[k]
                       + hk_687[k];

            t_702[k] = -ab_y[k] * hi_534[k]
                       + hk_688[k];

            t_703[k] = -ab_y[k] * hi_535[k]
                       + hk_690[k];
        }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, ab_y, hi_536, hi_537, hi_538, \
                         hi_539, hi_540, hk_691, hk_692, hk_694, hk_695, \
                         hk_696 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_704[k] = -ab_y[k] * hi_536[k]
                       + hk_691[k];

            t_705[k] = -ab_y[k] * hi_537[k]
                       + hk_692[k];

            t_706[k] = -ab_y[k] * hi_538[k]
                       + hk_694[k];

            t_707[k] = -ab_y[k] * hi_539[k]
                       + hk_695[k];

            t_708[k] = -ab_y[k] * hi_540[k]
                       + hk_696[k];
        }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, ab_y, hi_541, hi_542, hi_543, \
                         hi_544, hi_545, hk_697, hk_699, hk_700, hk_701, \
                         hk_702 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_709[k] = -ab_y[k] * hi_541[k]
                       + hk_697[k];

            t_710[k] = -ab_y[k] * hi_542[k]
                       + hk_699[k];

            t_711[k] = -ab_y[k] * hi_543[k]
                       + hk_700[k];

            t_712[k] = -ab_y[k] * hi_544[k]
                       + hk_701[k];

            t_713[k] = -ab_y[k] * hi_545[k]
                       + hk_702[k];
        }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, ab_y, hi_546, hi_547, hi_548, \
                         hi_549, hi_550, hk_703, hk_705, hk_706, hk_707, \
                         hk_708 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_714[k] = -ab_y[k] * hi_546[k]
                       + hk_703[k];

            t_715[k] = -ab_y[k] * hi_547[k]
                       + hk_705[k];

            t_716[k] = -ab_y[k] * hi_548[k]
                       + hk_706[k];

            t_717[k] = -ab_y[k] * hi_549[k]
                       + hk_707[k];

            t_718[k] = -ab_y[k] * hi_550[k]
                       + hk_708[k];
        }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, ab_y, hi_551, hi_552, hi_553, \
                         hi_554, hi_555, hk_709, hk_710, hk_712, hk_713, \
                         hk_714 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_719[k] = -ab_y[k] * hi_551[k]
                       + hk_709[k];

            t_720[k] = -ab_y[k] * hi_552[k]
                       + hk_710[k];

            t_721[k] = -ab_y[k] * hi_553[k]
                       + hk_712[k];

            t_722[k] = -ab_y[k] * hi_554[k]
                       + hk_713[k];

            t_723[k] = -ab_y[k] * hi_555[k]
                       + hk_714[k];
        }
    }
}

static auto
compute_hrr_ii_piece5(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hi, const size_t hk, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_724 = buffer.data(target + 724 * ncomps + c);
        auto *t_725 = buffer.data(target + 725 * ncomps + c);
        auto *t_726 = buffer.data(target + 726 * ncomps + c);
        auto *t_727 = buffer.data(target + 727 * ncomps + c);
        auto *t_728 = buffer.data(target + 728 * ncomps + c);
        auto *t_729 = buffer.data(target + 729 * ncomps + c);
        auto *t_730 = buffer.data(target + 730 * ncomps + c);
        auto *t_731 = buffer.data(target + 731 * ncomps + c);
        auto *t_732 = buffer.data(target + 732 * ncomps + c);
        auto *t_733 = buffer.data(target + 733 * ncomps + c);
        auto *t_734 = buffer.data(target + 734 * ncomps + c);
        auto *t_735 = buffer.data(target + 735 * ncomps + c);
        auto *t_736 = buffer.data(target + 736 * ncomps + c);
        auto *t_737 = buffer.data(target + 737 * ncomps + c);
        auto *t_738 = buffer.data(target + 738 * ncomps + c);
        auto *t_739 = buffer.data(target + 739 * ncomps + c);
        auto *t_740 = buffer.data(target + 740 * ncomps + c);
        auto *t_741 = buffer.data(target + 741 * ncomps + c);
        auto *t_742 = buffer.data(target + 742 * ncomps + c);
        auto *t_743 = buffer.data(target + 743 * ncomps + c);
        auto *t_744 = buffer.data(target + 744 * ncomps + c);
        auto *t_745 = buffer.data(target + 745 * ncomps + c);
        auto *t_746 = buffer.data(target + 746 * ncomps + c);
        auto *t_747 = buffer.data(target + 747 * ncomps + c);
        auto *t_748 = buffer.data(target + 748 * ncomps + c);
        auto *t_749 = buffer.data(target + 749 * ncomps + c);
        auto *t_750 = buffer.data(target + 750 * ncomps + c);
        auto *t_751 = buffer.data(target + 751 * ncomps + c);
        auto *t_752 = buffer.data(target + 752 * ncomps + c);
        auto *t_753 = buffer.data(target + 753 * ncomps + c);
        auto *t_754 = buffer.data(target + 754 * ncomps + c);
        auto *t_755 = buffer.data(target + 755 * ncomps + c);
        auto *t_756 = buffer.data(target + 756 * ncomps + c);
        auto *t_757 = buffer.data(target + 757 * ncomps + c);
        auto *t_758 = buffer.data(target + 758 * ncomps + c);
        auto *t_759 = buffer.data(target + 759 * ncomps + c);
        auto *t_760 = buffer.data(target + 760 * ncomps + c);
        auto *t_761 = buffer.data(target + 761 * ncomps + c);
        auto *t_762 = buffer.data(target + 762 * ncomps + c);
        auto *t_763 = buffer.data(target + 763 * ncomps + c);
        auto *t_764 = buffer.data(target + 764 * ncomps + c);
        auto *t_765 = buffer.data(target + 765 * ncomps + c);
        auto *t_766 = buffer.data(target + 766 * ncomps + c);
        auto *t_767 = buffer.data(target + 767 * ncomps + c);
        auto *t_768 = buffer.data(target + 768 * ncomps + c);
        auto *t_769 = buffer.data(target + 769 * ncomps + c);
        auto *t_770 = buffer.data(target + 770 * ncomps + c);
        auto *t_771 = buffer.data(target + 771 * ncomps + c);
        auto *t_772 = buffer.data(target + 772 * ncomps + c);
        auto *t_773 = buffer.data(target + 773 * ncomps + c);
        auto *t_774 = buffer.data(target + 774 * ncomps + c);
        auto *t_775 = buffer.data(target + 775 * ncomps + c);
        auto *t_776 = buffer.data(target + 776 * ncomps + c);
        auto *t_777 = buffer.data(target + 777 * ncomps + c);
        auto *t_778 = buffer.data(target + 778 * ncomps + c);
        auto *t_779 = buffer.data(target + 779 * ncomps + c);
        auto *t_780 = buffer.data(target + 780 * ncomps + c);
        auto *t_781 = buffer.data(target + 781 * ncomps + c);
        auto *t_782 = buffer.data(target + 782 * ncomps + c);
        auto *t_783 = buffer.data(target + 783 * ncomps + c);

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hi_556 = buffer.data(hi + 556 * ncomps + c);
        const auto *hi_557 = buffer.data(hi + 557 * ncomps + c);
        const auto *hi_558 = buffer.data(hi + 558 * ncomps + c);
        const auto *hi_559 = buffer.data(hi + 559 * ncomps + c);
        const auto *hi_560 = buffer.data(hi + 560 * ncomps + c);
        const auto *hi_561 = buffer.data(hi + 561 * ncomps + c);
        const auto *hi_562 = buffer.data(hi + 562 * ncomps + c);
        const auto *hi_563 = buffer.data(hi + 563 * ncomps + c);
        const auto *hi_564 = buffer.data(hi + 564 * ncomps + c);
        const auto *hi_565 = buffer.data(hi + 565 * ncomps + c);
        const auto *hi_566 = buffer.data(hi + 566 * ncomps + c);
        const auto *hi_567 = buffer.data(hi + 567 * ncomps + c);
        const auto *hi_568 = buffer.data(hi + 568 * ncomps + c);
        const auto *hi_569 = buffer.data(hi + 569 * ncomps + c);
        const auto *hi_570 = buffer.data(hi + 570 * ncomps + c);
        const auto *hi_571 = buffer.data(hi + 571 * ncomps + c);
        const auto *hi_572 = buffer.data(hi + 572 * ncomps + c);
        const auto *hi_573 = buffer.data(hi + 573 * ncomps + c);
        const auto *hi_574 = buffer.data(hi + 574 * ncomps + c);
        const auto *hi_575 = buffer.data(hi + 575 * ncomps + c);
        const auto *hi_576 = buffer.data(hi + 576 * ncomps + c);
        const auto *hi_577 = buffer.data(hi + 577 * ncomps + c);
        const auto *hi_578 = buffer.data(hi + 578 * ncomps + c);
        const auto *hi_579 = buffer.data(hi + 579 * ncomps + c);
        const auto *hi_580 = buffer.data(hi + 580 * ncomps + c);
        const auto *hi_581 = buffer.data(hi + 581 * ncomps + c);
        const auto *hi_582 = buffer.data(hi + 582 * ncomps + c);
        const auto *hi_583 = buffer.data(hi + 583 * ncomps + c);
        const auto *hi_584 = buffer.data(hi + 584 * ncomps + c);
        const auto *hi_585 = buffer.data(hi + 585 * ncomps + c);
        const auto *hi_586 = buffer.data(hi + 586 * ncomps + c);
        const auto *hi_587 = buffer.data(hi + 587 * ncomps + c);

        const auto *hk_715 = buffer.data(hk + 715 * ncomps + c);
        const auto *hk_716 = buffer.data(hk + 716 * ncomps + c);
        const auto *hk_717 = buffer.data(hk + 717 * ncomps + c);
        const auto *hk_718 = buffer.data(hk + 718 * ncomps + c);
        const auto *hk_721 = buffer.data(hk + 721 * ncomps + c);
        const auto *hk_722 = buffer.data(hk + 722 * ncomps + c);
        const auto *hk_723 = buffer.data(hk + 723 * ncomps + c);
        const auto *hk_724 = buffer.data(hk + 724 * ncomps + c);
        const auto *hk_725 = buffer.data(hk + 725 * ncomps + c);
        const auto *hk_726 = buffer.data(hk + 726 * ncomps + c);
        const auto *hk_727 = buffer.data(hk + 727 * ncomps + c);
        const auto *hk_728 = buffer.data(hk + 728 * ncomps + c);
        const auto *hk_729 = buffer.data(hk + 729 * ncomps + c);
        const auto *hk_730 = buffer.data(hk + 730 * ncomps + c);
        const auto *hk_731 = buffer.data(hk + 731 * ncomps + c);
        const auto *hk_732 = buffer.data(hk + 732 * ncomps + c);
        const auto *hk_733 = buffer.data(hk + 733 * ncomps + c);
        const auto *hk_734 = buffer.data(hk + 734 * ncomps + c);
        const auto *hk_735 = buffer.data(hk + 735 * ncomps + c);
        const auto *hk_736 = buffer.data(hk + 736 * ncomps + c);
        const auto *hk_737 = buffer.data(hk + 737 * ncomps + c);
        const auto *hk_738 = buffer.data(hk + 738 * ncomps + c);
        const auto *hk_739 = buffer.data(hk + 739 * ncomps + c);
        const auto *hk_740 = buffer.data(hk + 740 * ncomps + c);
        const auto *hk_741 = buffer.data(hk + 741 * ncomps + c);
        const auto *hk_742 = buffer.data(hk + 742 * ncomps + c);
        const auto *hk_743 = buffer.data(hk + 743 * ncomps + c);
        const auto *hk_744 = buffer.data(hk + 744 * ncomps + c);
        const auto *hk_745 = buffer.data(hk + 745 * ncomps + c);
        const auto *hk_746 = buffer.data(hk + 746 * ncomps + c);
        const auto *hk_747 = buffer.data(hk + 747 * ncomps + c);
        const auto *hk_748 = buffer.data(hk + 748 * ncomps + c);
        const auto *hk_749 = buffer.data(hk + 749 * ncomps + c);
        const auto *hk_750 = buffer.data(hk + 750 * ncomps + c);
        const auto *hk_751 = buffer.data(hk + 751 * ncomps + c);
        const auto *hk_752 = buffer.data(hk + 752 * ncomps + c);
        const auto *hk_753 = buffer.data(hk + 753 * ncomps + c);
        const auto *hk_754 = buffer.data(hk + 754 * ncomps + c);
        const auto *hk_755 = buffer.data(hk + 755 * ncomps + c);

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, ab_y, hi_556, hi_557, hi_558, \
                         hi_559, hi_560, hk_715, hk_716, hk_717, hk_718, \
                         hk_721 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_724[k] = -ab_y[k] * hi_556[k]
                       + hk_715[k];

            t_725[k] = -ab_y[k] * hi_557[k]
                       + hk_716[k];

            t_726[k] = -ab_y[k] * hi_558[k]
                       + hk_717[k];

            t_727[k] = -ab_y[k] * hi_559[k]
                       + hk_718[k];

            t_728[k] = -ab_y[k] * hi_560[k]
                       + hk_721[k];
        }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, ab_y, hi_561, hi_562, hi_563, \
                         hi_564, hi_565, hk_723, hk_724, hk_726, hk_727, \
                         hk_728 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_729[k] = -ab_y[k] * hi_561[k]
                       + hk_723[k];

            t_730[k] = -ab_y[k] * hi_562[k]
                       + hk_724[k];

            t_731[k] = -ab_y[k] * hi_563[k]
                       + hk_726[k];

            t_732[k] = -ab_y[k] * hi_564[k]
                       + hk_727[k];

            t_733[k] = -ab_y[k] * hi_565[k]
                       + hk_728[k];
        }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, t_738, ab_y, hi_566, hi_567, hi_568, \
                         hi_569, hi_570, hk_730, hk_731, hk_732, hk_733, \
                         hk_735 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_734[k] = -ab_y[k] * hi_566[k]
                       + hk_730[k];

            t_735[k] = -ab_y[k] * hi_567[k]
                       + hk_731[k];

            t_736[k] = -ab_y[k] * hi_568[k]
                       + hk_732[k];

            t_737[k] = -ab_y[k] * hi_569[k]
                       + hk_733[k];

            t_738[k] = -ab_y[k] * hi_570[k]
                       + hk_735[k];
        }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, ab_y, hi_571, hi_572, hi_573, \
                         hi_574, hi_575, hk_736, hk_737, hk_738, hk_739, \
                         hk_741 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_739[k] = -ab_y[k] * hi_571[k]
                       + hk_736[k];

            t_740[k] = -ab_y[k] * hi_572[k]
                       + hk_737[k];

            t_741[k] = -ab_y[k] * hi_573[k]
                       + hk_738[k];

            t_742[k] = -ab_y[k] * hi_574[k]
                       + hk_739[k];

            t_743[k] = -ab_y[k] * hi_575[k]
                       + hk_741[k];
        }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, ab_y, hi_576, hi_577, hi_578, \
                         hi_579, hi_580, hk_742, hk_743, hk_744, hk_745, \
                         hk_746 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_744[k] = -ab_y[k] * hi_576[k]
                       + hk_742[k];

            t_745[k] = -ab_y[k] * hi_577[k]
                       + hk_743[k];

            t_746[k] = -ab_y[k] * hi_578[k]
                       + hk_744[k];

            t_747[k] = -ab_y[k] * hi_579[k]
                       + hk_745[k];

            t_748[k] = -ab_y[k] * hi_580[k]
                       + hk_746[k];
        }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, t_753, ab_y, hi_581, hi_582, hi_583, \
                         hi_584, hi_585, hk_748, hk_749, hk_750, hk_751, \
                         hk_752 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_749[k] = -ab_y[k] * hi_581[k]
                       + hk_748[k];

            t_750[k] = -ab_y[k] * hi_582[k]
                       + hk_749[k];

            t_751[k] = -ab_y[k] * hi_583[k]
                       + hk_750[k];

            t_752[k] = -ab_y[k] * hi_584[k]
                       + hk_751[k];

            t_753[k] = -ab_y[k] * hi_585[k]
                       + hk_752[k];
        }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, ab_y, ab_z, hi_560, hi_561, hi_586, \
                         hi_587, hk_722, hk_724, hk_753, hk_754 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_754[k] = -ab_y[k] * hi_586[k]
                       + hk_753[k];

            t_755[k] = -ab_y[k] * hi_587[k]
                       + hk_754[k];

            t_756[k] = -ab_z[k] * hi_560[k]
                       + hk_722[k];

            t_757[k] = -ab_z[k] * hi_561[k]
                       + hk_724[k];
        }

#pragma omp simd aligned(t_758, t_759, t_760, t_761, t_762, ab_z, hi_562, hi_563, hi_564, \
                         hi_565, hi_566, hk_725, hk_727, hk_728, hk_729, \
                         hk_731 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_758[k] = -ab_z[k] * hi_562[k]
                       + hk_725[k];

            t_759[k] = -ab_z[k] * hi_563[k]
                       + hk_727[k];

            t_760[k] = -ab_z[k] * hi_564[k]
                       + hk_728[k];

            t_761[k] = -ab_z[k] * hi_565[k]
                       + hk_729[k];

            t_762[k] = -ab_z[k] * hi_566[k]
                       + hk_731[k];
        }

#pragma omp simd aligned(t_763, t_764, t_765, t_766, t_767, ab_z, hi_567, hi_568, hi_569, \
                         hi_570, hi_571, hk_732, hk_733, hk_734, hk_736, \
                         hk_737 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_763[k] = -ab_z[k] * hi_567[k]
                       + hk_732[k];

            t_764[k] = -ab_z[k] * hi_568[k]
                       + hk_733[k];

            t_765[k] = -ab_z[k] * hi_569[k]
                       + hk_734[k];

            t_766[k] = -ab_z[k] * hi_570[k]
                       + hk_736[k];

            t_767[k] = -ab_z[k] * hi_571[k]
                       + hk_737[k];
        }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, t_772, ab_z, hi_572, hi_573, hi_574, \
                         hi_575, hi_576, hk_738, hk_739, hk_740, hk_742, \
                         hk_743 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_768[k] = -ab_z[k] * hi_572[k]
                       + hk_738[k];

            t_769[k] = -ab_z[k] * hi_573[k]
                       + hk_739[k];

            t_770[k] = -ab_z[k] * hi_574[k]
                       + hk_740[k];

            t_771[k] = -ab_z[k] * hi_575[k]
                       + hk_742[k];

            t_772[k] = -ab_z[k] * hi_576[k]
                       + hk_743[k];
        }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, t_777, ab_z, hi_577, hi_578, hi_579, \
                         hi_580, hi_581, hk_744, hk_745, hk_746, hk_747, \
                         hk_749 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_773[k] = -ab_z[k] * hi_577[k]
                       + hk_744[k];

            t_774[k] = -ab_z[k] * hi_578[k]
                       + hk_745[k];

            t_775[k] = -ab_z[k] * hi_579[k]
                       + hk_746[k];

            t_776[k] = -ab_z[k] * hi_580[k]
                       + hk_747[k];

            t_777[k] = -ab_z[k] * hi_581[k]
                       + hk_749[k];
        }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, ab_z, hi_582, hi_583, hi_584, \
                         hi_585, hi_586, hk_750, hk_751, hk_752, hk_753, \
                         hk_754 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_778[k] = -ab_z[k] * hi_582[k]
                       + hk_750[k];

            t_779[k] = -ab_z[k] * hi_583[k]
                       + hk_751[k];

            t_780[k] = -ab_z[k] * hi_584[k]
                       + hk_752[k];

            t_781[k] = -ab_z[k] * hi_585[k]
                       + hk_753[k];

            t_782[k] = -ab_z[k] * hi_586[k]
                       + hk_754[k];
        }

#pragma omp simd aligned(t_783, ab_z, hi_587, hk_755 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_783[k] = -ab_z[k] * hi_587[k]
                       + hk_755[k];
        }
    }
}

auto
compute_hrr_ii(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t hi, const size_t hk, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_ii_piece0(buffer, coordinates, target, hi, hk, ncomps, nmax);

    compute_hrr_ii_piece1(buffer, coordinates, target, hi, hk, ncomps, nmax);

    compute_hrr_ii_piece2(buffer, coordinates, target, hi, hk, ncomps, nmax);

    compute_hrr_ii_piece3(buffer, coordinates, target, hi, hk, ncomps, nmax);

    compute_hrr_ii_piece4(buffer, coordinates, target, hi, hk, ncomps, nmax);

    compute_hrr_ii_piece5(buffer, coordinates, target, hi, hk, ncomps, nmax);
}

}  // namespace simdtrf
