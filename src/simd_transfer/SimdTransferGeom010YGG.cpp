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


#include "SimdTransferGeom010YGG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_010y_gg_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t fg_1, const size_t fh_1,
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

        const auto *fg_1_0 = buffer.data(fg_1 + 0 * ncomps + c);
        const auto *fg_1_1 = buffer.data(fg_1 + 1 * ncomps + c);
        const auto *fg_1_2 = buffer.data(fg_1 + 2 * ncomps + c);
        const auto *fg_1_3 = buffer.data(fg_1 + 3 * ncomps + c);
        const auto *fg_1_4 = buffer.data(fg_1 + 4 * ncomps + c);
        const auto *fg_1_5 = buffer.data(fg_1 + 5 * ncomps + c);
        const auto *fg_1_6 = buffer.data(fg_1 + 6 * ncomps + c);
        const auto *fg_1_7 = buffer.data(fg_1 + 7 * ncomps + c);
        const auto *fg_1_8 = buffer.data(fg_1 + 8 * ncomps + c);
        const auto *fg_1_9 = buffer.data(fg_1 + 9 * ncomps + c);
        const auto *fg_1_10 = buffer.data(fg_1 + 10 * ncomps + c);
        const auto *fg_1_11 = buffer.data(fg_1 + 11 * ncomps + c);
        const auto *fg_1_12 = buffer.data(fg_1 + 12 * ncomps + c);
        const auto *fg_1_13 = buffer.data(fg_1 + 13 * ncomps + c);
        const auto *fg_1_14 = buffer.data(fg_1 + 14 * ncomps + c);
        const auto *fg_1_15 = buffer.data(fg_1 + 15 * ncomps + c);
        const auto *fg_1_16 = buffer.data(fg_1 + 16 * ncomps + c);
        const auto *fg_1_17 = buffer.data(fg_1 + 17 * ncomps + c);
        const auto *fg_1_18 = buffer.data(fg_1 + 18 * ncomps + c);
        const auto *fg_1_19 = buffer.data(fg_1 + 19 * ncomps + c);
        const auto *fg_1_20 = buffer.data(fg_1 + 20 * ncomps + c);
        const auto *fg_1_21 = buffer.data(fg_1 + 21 * ncomps + c);
        const auto *fg_1_22 = buffer.data(fg_1 + 22 * ncomps + c);
        const auto *fg_1_23 = buffer.data(fg_1 + 23 * ncomps + c);
        const auto *fg_1_24 = buffer.data(fg_1 + 24 * ncomps + c);
        const auto *fg_1_25 = buffer.data(fg_1 + 25 * ncomps + c);
        const auto *fg_1_26 = buffer.data(fg_1 + 26 * ncomps + c);
        const auto *fg_1_27 = buffer.data(fg_1 + 27 * ncomps + c);
        const auto *fg_1_28 = buffer.data(fg_1 + 28 * ncomps + c);
        const auto *fg_1_29 = buffer.data(fg_1 + 29 * ncomps + c);
        const auto *fg_1_30 = buffer.data(fg_1 + 30 * ncomps + c);
        const auto *fg_1_31 = buffer.data(fg_1 + 31 * ncomps + c);
        const auto *fg_1_32 = buffer.data(fg_1 + 32 * ncomps + c);
        const auto *fg_1_33 = buffer.data(fg_1 + 33 * ncomps + c);
        const auto *fg_1_34 = buffer.data(fg_1 + 34 * ncomps + c);
        const auto *fg_1_35 = buffer.data(fg_1 + 35 * ncomps + c);
        const auto *fg_1_36 = buffer.data(fg_1 + 36 * ncomps + c);
        const auto *fg_1_37 = buffer.data(fg_1 + 37 * ncomps + c);
        const auto *fg_1_38 = buffer.data(fg_1 + 38 * ncomps + c);
        const auto *fg_1_39 = buffer.data(fg_1 + 39 * ncomps + c);
        const auto *fg_1_40 = buffer.data(fg_1 + 40 * ncomps + c);
        const auto *fg_1_41 = buffer.data(fg_1 + 41 * ncomps + c);
        const auto *fg_1_42 = buffer.data(fg_1 + 42 * ncomps + c);
        const auto *fg_1_43 = buffer.data(fg_1 + 43 * ncomps + c);
        const auto *fg_1_44 = buffer.data(fg_1 + 44 * ncomps + c);
        const auto *fg_1_45 = buffer.data(fg_1 + 45 * ncomps + c);
        const auto *fg_1_46 = buffer.data(fg_1 + 46 * ncomps + c);
        const auto *fg_1_47 = buffer.data(fg_1 + 47 * ncomps + c);
        const auto *fg_1_48 = buffer.data(fg_1 + 48 * ncomps + c);
        const auto *fg_1_49 = buffer.data(fg_1 + 49 * ncomps + c);
        const auto *fg_1_50 = buffer.data(fg_1 + 50 * ncomps + c);
        const auto *fg_1_51 = buffer.data(fg_1 + 51 * ncomps + c);
        const auto *fg_1_52 = buffer.data(fg_1 + 52 * ncomps + c);
        const auto *fg_1_53 = buffer.data(fg_1 + 53 * ncomps + c);
        const auto *fg_1_54 = buffer.data(fg_1 + 54 * ncomps + c);
        const auto *fg_1_55 = buffer.data(fg_1 + 55 * ncomps + c);
        const auto *fg_1_56 = buffer.data(fg_1 + 56 * ncomps + c);
        const auto *fg_1_57 = buffer.data(fg_1 + 57 * ncomps + c);
        const auto *fg_1_58 = buffer.data(fg_1 + 58 * ncomps + c);
        const auto *fg_1_59 = buffer.data(fg_1 + 59 * ncomps + c);
        const auto *fg_1_60 = buffer.data(fg_1 + 60 * ncomps + c);
        const auto *fg_1_61 = buffer.data(fg_1 + 61 * ncomps + c);
        const auto *fg_1_62 = buffer.data(fg_1 + 62 * ncomps + c);
        const auto *fg_1_63 = buffer.data(fg_1 + 63 * ncomps + c);
        const auto *fg_1_64 = buffer.data(fg_1 + 64 * ncomps + c);
        const auto *fg_1_65 = buffer.data(fg_1 + 65 * ncomps + c);
        const auto *fg_1_66 = buffer.data(fg_1 + 66 * ncomps + c);
        const auto *fg_1_67 = buffer.data(fg_1 + 67 * ncomps + c);
        const auto *fg_1_68 = buffer.data(fg_1 + 68 * ncomps + c);
        const auto *fg_1_69 = buffer.data(fg_1 + 69 * ncomps + c);
        const auto *fg_1_70 = buffer.data(fg_1 + 70 * ncomps + c);
        const auto *fg_1_71 = buffer.data(fg_1 + 71 * ncomps + c);
        const auto *fg_1_72 = buffer.data(fg_1 + 72 * ncomps + c);
        const auto *fg_1_73 = buffer.data(fg_1 + 73 * ncomps + c);
        const auto *fg_1_74 = buffer.data(fg_1 + 74 * ncomps + c);
        const auto *fg_1_75 = buffer.data(fg_1 + 75 * ncomps + c);
        const auto *fg_1_76 = buffer.data(fg_1 + 76 * ncomps + c);
        const auto *fg_1_77 = buffer.data(fg_1 + 77 * ncomps + c);
        const auto *fg_1_78 = buffer.data(fg_1 + 78 * ncomps + c);
        const auto *fg_1_79 = buffer.data(fg_1 + 79 * ncomps + c);
        const auto *fg_1_80 = buffer.data(fg_1 + 80 * ncomps + c);
        const auto *fg_1_81 = buffer.data(fg_1 + 81 * ncomps + c);
        const auto *fg_1_82 = buffer.data(fg_1 + 82 * ncomps + c);
        const auto *fg_1_83 = buffer.data(fg_1 + 83 * ncomps + c);
        const auto *fg_1_84 = buffer.data(fg_1 + 84 * ncomps + c);
        const auto *fg_1_85 = buffer.data(fg_1 + 85 * ncomps + c);
        const auto *fg_1_86 = buffer.data(fg_1 + 86 * ncomps + c);
        const auto *fg_1_87 = buffer.data(fg_1 + 87 * ncomps + c);
        const auto *fg_1_88 = buffer.data(fg_1 + 88 * ncomps + c);
        const auto *fg_1_89 = buffer.data(fg_1 + 89 * ncomps + c);
        const auto *fg_1_90 = buffer.data(fg_1 + 90 * ncomps + c);
        const auto *fg_1_91 = buffer.data(fg_1 + 91 * ncomps + c);
        const auto *fg_1_92 = buffer.data(fg_1 + 92 * ncomps + c);
        const auto *fg_1_93 = buffer.data(fg_1 + 93 * ncomps + c);
        const auto *fg_1_94 = buffer.data(fg_1 + 94 * ncomps + c);
        const auto *fg_1_95 = buffer.data(fg_1 + 95 * ncomps + c);
        const auto *fg_1_96 = buffer.data(fg_1 + 96 * ncomps + c);
        const auto *fg_1_97 = buffer.data(fg_1 + 97 * ncomps + c);
        const auto *fg_1_98 = buffer.data(fg_1 + 98 * ncomps + c);
        const auto *fg_1_99 = buffer.data(fg_1 + 99 * ncomps + c);
        const auto *fg_1_100 = buffer.data(fg_1 + 100 * ncomps + c);
        const auto *fg_1_101 = buffer.data(fg_1 + 101 * ncomps + c);
        const auto *fg_1_102 = buffer.data(fg_1 + 102 * ncomps + c);
        const auto *fg_1_103 = buffer.data(fg_1 + 103 * ncomps + c);
        const auto *fg_1_104 = buffer.data(fg_1 + 104 * ncomps + c);
        const auto *fg_1_105 = buffer.data(fg_1 + 105 * ncomps + c);
        const auto *fg_1_106 = buffer.data(fg_1 + 106 * ncomps + c);
        const auto *fg_1_107 = buffer.data(fg_1 + 107 * ncomps + c);
        const auto *fg_1_108 = buffer.data(fg_1 + 108 * ncomps + c);
        const auto *fg_1_109 = buffer.data(fg_1 + 109 * ncomps + c);
        const auto *fg_1_110 = buffer.data(fg_1 + 110 * ncomps + c);
        const auto *fg_1_111 = buffer.data(fg_1 + 111 * ncomps + c);
        const auto *fg_1_112 = buffer.data(fg_1 + 112 * ncomps + c);
        const auto *fg_1_113 = buffer.data(fg_1 + 113 * ncomps + c);
        const auto *fg_1_114 = buffer.data(fg_1 + 114 * ncomps + c);
        const auto *fg_1_115 = buffer.data(fg_1 + 115 * ncomps + c);
        const auto *fg_1_116 = buffer.data(fg_1 + 116 * ncomps + c);
        const auto *fg_1_117 = buffer.data(fg_1 + 117 * ncomps + c);
        const auto *fg_1_118 = buffer.data(fg_1 + 118 * ncomps + c);
        const auto *fg_1_119 = buffer.data(fg_1 + 119 * ncomps + c);
        const auto *fg_1_120 = buffer.data(fg_1 + 120 * ncomps + c);
        const auto *fg_1_121 = buffer.data(fg_1 + 121 * ncomps + c);
        const auto *fg_1_122 = buffer.data(fg_1 + 122 * ncomps + c);
        const auto *fg_1_123 = buffer.data(fg_1 + 123 * ncomps + c);
        const auto *fg_1_124 = buffer.data(fg_1 + 124 * ncomps + c);
        const auto *fg_1_125 = buffer.data(fg_1 + 125 * ncomps + c);
        const auto *fg_1_126 = buffer.data(fg_1 + 126 * ncomps + c);
        const auto *fg_1_127 = buffer.data(fg_1 + 127 * ncomps + c);
        const auto *fg_1_128 = buffer.data(fg_1 + 128 * ncomps + c);
        const auto *fg_1_129 = buffer.data(fg_1 + 129 * ncomps + c);
        const auto *fg_1_130 = buffer.data(fg_1 + 130 * ncomps + c);
        const auto *fg_1_131 = buffer.data(fg_1 + 131 * ncomps + c);
        const auto *fg_1_132 = buffer.data(fg_1 + 132 * ncomps + c);
        const auto *fg_1_133 = buffer.data(fg_1 + 133 * ncomps + c);
        const auto *fg_1_134 = buffer.data(fg_1 + 134 * ncomps + c);
        const auto *fg_1_135 = buffer.data(fg_1 + 135 * ncomps + c);
        const auto *fg_1_136 = buffer.data(fg_1 + 136 * ncomps + c);
        const auto *fg_1_137 = buffer.data(fg_1 + 137 * ncomps + c);
        const auto *fg_1_138 = buffer.data(fg_1 + 138 * ncomps + c);
        const auto *fg_1_139 = buffer.data(fg_1 + 139 * ncomps + c);
        const auto *fg_1_140 = buffer.data(fg_1 + 140 * ncomps + c);
        const auto *fg_1_141 = buffer.data(fg_1 + 141 * ncomps + c);
        const auto *fg_1_142 = buffer.data(fg_1 + 142 * ncomps + c);
        const auto *fg_1_143 = buffer.data(fg_1 + 143 * ncomps + c);
        const auto *fg_1_144 = buffer.data(fg_1 + 144 * ncomps + c);

        const auto *fh_1_0 = buffer.data(fh_1 + 0 * ncomps + c);
        const auto *fh_1_1 = buffer.data(fh_1 + 1 * ncomps + c);
        const auto *fh_1_2 = buffer.data(fh_1 + 2 * ncomps + c);
        const auto *fh_1_3 = buffer.data(fh_1 + 3 * ncomps + c);
        const auto *fh_1_4 = buffer.data(fh_1 + 4 * ncomps + c);
        const auto *fh_1_5 = buffer.data(fh_1 + 5 * ncomps + c);
        const auto *fh_1_6 = buffer.data(fh_1 + 6 * ncomps + c);
        const auto *fh_1_7 = buffer.data(fh_1 + 7 * ncomps + c);
        const auto *fh_1_8 = buffer.data(fh_1 + 8 * ncomps + c);
        const auto *fh_1_9 = buffer.data(fh_1 + 9 * ncomps + c);
        const auto *fh_1_10 = buffer.data(fh_1 + 10 * ncomps + c);
        const auto *fh_1_11 = buffer.data(fh_1 + 11 * ncomps + c);
        const auto *fh_1_12 = buffer.data(fh_1 + 12 * ncomps + c);
        const auto *fh_1_13 = buffer.data(fh_1 + 13 * ncomps + c);
        const auto *fh_1_14 = buffer.data(fh_1 + 14 * ncomps + c);
        const auto *fh_1_21 = buffer.data(fh_1 + 21 * ncomps + c);
        const auto *fh_1_22 = buffer.data(fh_1 + 22 * ncomps + c);
        const auto *fh_1_23 = buffer.data(fh_1 + 23 * ncomps + c);
        const auto *fh_1_24 = buffer.data(fh_1 + 24 * ncomps + c);
        const auto *fh_1_25 = buffer.data(fh_1 + 25 * ncomps + c);
        const auto *fh_1_26 = buffer.data(fh_1 + 26 * ncomps + c);
        const auto *fh_1_27 = buffer.data(fh_1 + 27 * ncomps + c);
        const auto *fh_1_28 = buffer.data(fh_1 + 28 * ncomps + c);
        const auto *fh_1_29 = buffer.data(fh_1 + 29 * ncomps + c);
        const auto *fh_1_30 = buffer.data(fh_1 + 30 * ncomps + c);
        const auto *fh_1_31 = buffer.data(fh_1 + 31 * ncomps + c);
        const auto *fh_1_32 = buffer.data(fh_1 + 32 * ncomps + c);
        const auto *fh_1_33 = buffer.data(fh_1 + 33 * ncomps + c);
        const auto *fh_1_34 = buffer.data(fh_1 + 34 * ncomps + c);
        const auto *fh_1_35 = buffer.data(fh_1 + 35 * ncomps + c);
        const auto *fh_1_42 = buffer.data(fh_1 + 42 * ncomps + c);
        const auto *fh_1_43 = buffer.data(fh_1 + 43 * ncomps + c);
        const auto *fh_1_44 = buffer.data(fh_1 + 44 * ncomps + c);
        const auto *fh_1_45 = buffer.data(fh_1 + 45 * ncomps + c);
        const auto *fh_1_46 = buffer.data(fh_1 + 46 * ncomps + c);
        const auto *fh_1_47 = buffer.data(fh_1 + 47 * ncomps + c);
        const auto *fh_1_48 = buffer.data(fh_1 + 48 * ncomps + c);
        const auto *fh_1_49 = buffer.data(fh_1 + 49 * ncomps + c);
        const auto *fh_1_50 = buffer.data(fh_1 + 50 * ncomps + c);
        const auto *fh_1_51 = buffer.data(fh_1 + 51 * ncomps + c);
        const auto *fh_1_52 = buffer.data(fh_1 + 52 * ncomps + c);
        const auto *fh_1_53 = buffer.data(fh_1 + 53 * ncomps + c);
        const auto *fh_1_54 = buffer.data(fh_1 + 54 * ncomps + c);
        const auto *fh_1_55 = buffer.data(fh_1 + 55 * ncomps + c);
        const auto *fh_1_56 = buffer.data(fh_1 + 56 * ncomps + c);
        const auto *fh_1_63 = buffer.data(fh_1 + 63 * ncomps + c);
        const auto *fh_1_64 = buffer.data(fh_1 + 64 * ncomps + c);
        const auto *fh_1_65 = buffer.data(fh_1 + 65 * ncomps + c);
        const auto *fh_1_66 = buffer.data(fh_1 + 66 * ncomps + c);
        const auto *fh_1_67 = buffer.data(fh_1 + 67 * ncomps + c);
        const auto *fh_1_68 = buffer.data(fh_1 + 68 * ncomps + c);
        const auto *fh_1_69 = buffer.data(fh_1 + 69 * ncomps + c);
        const auto *fh_1_70 = buffer.data(fh_1 + 70 * ncomps + c);
        const auto *fh_1_71 = buffer.data(fh_1 + 71 * ncomps + c);
        const auto *fh_1_72 = buffer.data(fh_1 + 72 * ncomps + c);
        const auto *fh_1_73 = buffer.data(fh_1 + 73 * ncomps + c);
        const auto *fh_1_74 = buffer.data(fh_1 + 74 * ncomps + c);
        const auto *fh_1_75 = buffer.data(fh_1 + 75 * ncomps + c);
        const auto *fh_1_76 = buffer.data(fh_1 + 76 * ncomps + c);
        const auto *fh_1_77 = buffer.data(fh_1 + 77 * ncomps + c);
        const auto *fh_1_84 = buffer.data(fh_1 + 84 * ncomps + c);
        const auto *fh_1_85 = buffer.data(fh_1 + 85 * ncomps + c);
        const auto *fh_1_86 = buffer.data(fh_1 + 86 * ncomps + c);
        const auto *fh_1_87 = buffer.data(fh_1 + 87 * ncomps + c);
        const auto *fh_1_88 = buffer.data(fh_1 + 88 * ncomps + c);
        const auto *fh_1_89 = buffer.data(fh_1 + 89 * ncomps + c);
        const auto *fh_1_90 = buffer.data(fh_1 + 90 * ncomps + c);
        const auto *fh_1_91 = buffer.data(fh_1 + 91 * ncomps + c);
        const auto *fh_1_92 = buffer.data(fh_1 + 92 * ncomps + c);
        const auto *fh_1_93 = buffer.data(fh_1 + 93 * ncomps + c);
        const auto *fh_1_94 = buffer.data(fh_1 + 94 * ncomps + c);
        const auto *fh_1_95 = buffer.data(fh_1 + 95 * ncomps + c);
        const auto *fh_1_96 = buffer.data(fh_1 + 96 * ncomps + c);
        const auto *fh_1_97 = buffer.data(fh_1 + 97 * ncomps + c);
        const auto *fh_1_98 = buffer.data(fh_1 + 98 * ncomps + c);
        const auto *fh_1_105 = buffer.data(fh_1 + 105 * ncomps + c);
        const auto *fh_1_106 = buffer.data(fh_1 + 106 * ncomps + c);
        const auto *fh_1_107 = buffer.data(fh_1 + 107 * ncomps + c);
        const auto *fh_1_108 = buffer.data(fh_1 + 108 * ncomps + c);
        const auto *fh_1_109 = buffer.data(fh_1 + 109 * ncomps + c);
        const auto *fh_1_110 = buffer.data(fh_1 + 110 * ncomps + c);
        const auto *fh_1_111 = buffer.data(fh_1 + 111 * ncomps + c);
        const auto *fh_1_112 = buffer.data(fh_1 + 112 * ncomps + c);
        const auto *fh_1_113 = buffer.data(fh_1 + 113 * ncomps + c);
        const auto *fh_1_114 = buffer.data(fh_1 + 114 * ncomps + c);
        const auto *fh_1_115 = buffer.data(fh_1 + 115 * ncomps + c);
        const auto *fh_1_116 = buffer.data(fh_1 + 116 * ncomps + c);
        const auto *fh_1_117 = buffer.data(fh_1 + 117 * ncomps + c);
        const auto *fh_1_118 = buffer.data(fh_1 + 118 * ncomps + c);
        const auto *fh_1_119 = buffer.data(fh_1 + 119 * ncomps + c);
        const auto *fh_1_126 = buffer.data(fh_1 + 126 * ncomps + c);
        const auto *fh_1_127 = buffer.data(fh_1 + 127 * ncomps + c);
        const auto *fh_1_128 = buffer.data(fh_1 + 128 * ncomps + c);
        const auto *fh_1_129 = buffer.data(fh_1 + 129 * ncomps + c);
        const auto *fh_1_130 = buffer.data(fh_1 + 130 * ncomps + c);
        const auto *fh_1_131 = buffer.data(fh_1 + 131 * ncomps + c);
        const auto *fh_1_132 = buffer.data(fh_1 + 132 * ncomps + c);
        const auto *fh_1_133 = buffer.data(fh_1 + 133 * ncomps + c);
        const auto *fh_1_134 = buffer.data(fh_1 + 134 * ncomps + c);
        const auto *fh_1_135 = buffer.data(fh_1 + 135 * ncomps + c);
        const auto *fh_1_136 = buffer.data(fh_1 + 136 * ncomps + c);
        const auto *fh_1_137 = buffer.data(fh_1 + 137 * ncomps + c);
        const auto *fh_1_138 = buffer.data(fh_1 + 138 * ncomps + c);
        const auto *fh_1_139 = buffer.data(fh_1 + 139 * ncomps + c);
        const auto *fh_1_140 = buffer.data(fh_1 + 140 * ncomps + c);
        const auto *fh_1_147 = buffer.data(fh_1 + 147 * ncomps + c);
        const auto *fh_1_148 = buffer.data(fh_1 + 148 * ncomps + c);
        const auto *fh_1_149 = buffer.data(fh_1 + 149 * ncomps + c);
        const auto *fh_1_150 = buffer.data(fh_1 + 150 * ncomps + c);
        const auto *fh_1_151 = buffer.data(fh_1 + 151 * ncomps + c);
        const auto *fh_1_152 = buffer.data(fh_1 + 152 * ncomps + c);
        const auto *fh_1_153 = buffer.data(fh_1 + 153 * ncomps + c);
        const auto *fh_1_154 = buffer.data(fh_1 + 154 * ncomps + c);
        const auto *fh_1_155 = buffer.data(fh_1 + 155 * ncomps + c);
        const auto *fh_1_156 = buffer.data(fh_1 + 156 * ncomps + c);
        const auto *fh_1_157 = buffer.data(fh_1 + 157 * ncomps + c);
        const auto *fh_1_158 = buffer.data(fh_1 + 158 * ncomps + c);
        const auto *fh_1_159 = buffer.data(fh_1 + 159 * ncomps + c);
        const auto *fh_1_160 = buffer.data(fh_1 + 160 * ncomps + c);
        const auto *fh_1_161 = buffer.data(fh_1 + 161 * ncomps + c);
        const auto *fh_1_168 = buffer.data(fh_1 + 168 * ncomps + c);
        const auto *fh_1_169 = buffer.data(fh_1 + 169 * ncomps + c);
        const auto *fh_1_170 = buffer.data(fh_1 + 170 * ncomps + c);
        const auto *fh_1_171 = buffer.data(fh_1 + 171 * ncomps + c);
        const auto *fh_1_172 = buffer.data(fh_1 + 172 * ncomps + c);
        const auto *fh_1_173 = buffer.data(fh_1 + 173 * ncomps + c);
        const auto *fh_1_174 = buffer.data(fh_1 + 174 * ncomps + c);
        const auto *fh_1_175 = buffer.data(fh_1 + 175 * ncomps + c);
        const auto *fh_1_176 = buffer.data(fh_1 + 176 * ncomps + c);
        const auto *fh_1_177 = buffer.data(fh_1 + 177 * ncomps + c);
        const auto *fh_1_178 = buffer.data(fh_1 + 178 * ncomps + c);
        const auto *fh_1_179 = buffer.data(fh_1 + 179 * ncomps + c);
        const auto *fh_1_180 = buffer.data(fh_1 + 180 * ncomps + c);
        const auto *fh_1_181 = buffer.data(fh_1 + 181 * ncomps + c);
        const auto *fh_1_182 = buffer.data(fh_1 + 182 * ncomps + c);
        const auto *fh_1_189 = buffer.data(fh_1 + 189 * ncomps + c);
        const auto *fh_1_190 = buffer.data(fh_1 + 190 * ncomps + c);
        const auto *fh_1_191 = buffer.data(fh_1 + 191 * ncomps + c);
        const auto *fh_1_192 = buffer.data(fh_1 + 192 * ncomps + c);
        const auto *fh_1_193 = buffer.data(fh_1 + 193 * ncomps + c);
        const auto *fh_1_194 = buffer.data(fh_1 + 194 * ncomps + c);
        const auto *fh_1_195 = buffer.data(fh_1 + 195 * ncomps + c);
        const auto *fh_1_196 = buffer.data(fh_1 + 196 * ncomps + c);
        const auto *fh_1_197 = buffer.data(fh_1 + 197 * ncomps + c);
        const auto *fh_1_198 = buffer.data(fh_1 + 198 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fg_1_0, fg_1_1, fg_1_2, fg_1_3, \
                         fg_1_4, fh_1_0, fh_1_1, fh_1_2, fh_1_3, \
                         fh_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * fg_1_0[k]
                     + fh_1_0[k];

            t_1[k] = -ab_x[k] * fg_1_1[k]
                     + fh_1_1[k];

            t_2[k] = -ab_x[k] * fg_1_2[k]
                     + fh_1_2[k];

            t_3[k] = -ab_x[k] * fg_1_3[k]
                     + fh_1_3[k];

            t_4[k] = -ab_x[k] * fg_1_4[k]
                     + fh_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fg_1_5, fg_1_6, fg_1_7, fg_1_8, \
                         fg_1_9, fh_1_5, fh_1_6, fh_1_7, fh_1_8, \
                         fh_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * fg_1_5[k]
                     + fh_1_5[k];

            t_6[k] = -ab_x[k] * fg_1_6[k]
                     + fh_1_6[k];

            t_7[k] = -ab_x[k] * fg_1_7[k]
                     + fh_1_7[k];

            t_8[k] = -ab_x[k] * fg_1_8[k]
                     + fh_1_8[k];

            t_9[k] = -ab_x[k] * fg_1_9[k]
                     + fh_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fg_1_10, fg_1_11, fg_1_12, \
                         fg_1_13, fg_1_14, fh_1_10, fh_1_11, fh_1_12, fh_1_13, \
                         fh_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * fg_1_10[k]
                      + fh_1_10[k];

            t_11[k] = -ab_x[k] * fg_1_11[k]
                      + fh_1_11[k];

            t_12[k] = -ab_x[k] * fg_1_12[k]
                      + fh_1_12[k];

            t_13[k] = -ab_x[k] * fg_1_13[k]
                      + fh_1_13[k];

            t_14[k] = -ab_x[k] * fg_1_14[k]
                      + fh_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fg_1_15, fg_1_16, fg_1_17, \
                         fg_1_18, fg_1_19, fh_1_21, fh_1_22, fh_1_23, fh_1_24, \
                         fh_1_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * fg_1_15[k]
                      + fh_1_21[k];

            t_16[k] = -ab_x[k] * fg_1_16[k]
                      + fh_1_22[k];

            t_17[k] = -ab_x[k] * fg_1_17[k]
                      + fh_1_23[k];

            t_18[k] = -ab_x[k] * fg_1_18[k]
                      + fh_1_24[k];

            t_19[k] = -ab_x[k] * fg_1_19[k]
                      + fh_1_25[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fg_1_20, fg_1_21, fg_1_22, \
                         fg_1_23, fg_1_24, fh_1_26, fh_1_27, fh_1_28, fh_1_29, \
                         fh_1_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * fg_1_20[k]
                      + fh_1_26[k];

            t_21[k] = -ab_x[k] * fg_1_21[k]
                      + fh_1_27[k];

            t_22[k] = -ab_x[k] * fg_1_22[k]
                      + fh_1_28[k];

            t_23[k] = -ab_x[k] * fg_1_23[k]
                      + fh_1_29[k];

            t_24[k] = -ab_x[k] * fg_1_24[k]
                      + fh_1_30[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fg_1_25, fg_1_26, fg_1_27, \
                         fg_1_28, fg_1_29, fh_1_31, fh_1_32, fh_1_33, fh_1_34, \
                         fh_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * fg_1_25[k]
                      + fh_1_31[k];

            t_26[k] = -ab_x[k] * fg_1_26[k]
                      + fh_1_32[k];

            t_27[k] = -ab_x[k] * fg_1_27[k]
                      + fh_1_33[k];

            t_28[k] = -ab_x[k] * fg_1_28[k]
                      + fh_1_34[k];

            t_29[k] = -ab_x[k] * fg_1_29[k]
                      + fh_1_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fg_1_30, fg_1_31, fg_1_32, \
                         fg_1_33, fg_1_34, fh_1_42, fh_1_43, fh_1_44, fh_1_45, \
                         fh_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * fg_1_30[k]
                      + fh_1_42[k];

            t_31[k] = -ab_x[k] * fg_1_31[k]
                      + fh_1_43[k];

            t_32[k] = -ab_x[k] * fg_1_32[k]
                      + fh_1_44[k];

            t_33[k] = -ab_x[k] * fg_1_33[k]
                      + fh_1_45[k];

            t_34[k] = -ab_x[k] * fg_1_34[k]
                      + fh_1_46[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, fg_1_35, fg_1_36, fg_1_37, \
                         fg_1_38, fg_1_39, fh_1_47, fh_1_48, fh_1_49, fh_1_50, \
                         fh_1_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * fg_1_35[k]
                      + fh_1_47[k];

            t_36[k] = -ab_x[k] * fg_1_36[k]
                      + fh_1_48[k];

            t_37[k] = -ab_x[k] * fg_1_37[k]
                      + fh_1_49[k];

            t_38[k] = -ab_x[k] * fg_1_38[k]
                      + fh_1_50[k];

            t_39[k] = -ab_x[k] * fg_1_39[k]
                      + fh_1_51[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fg_1_40, fg_1_41, fg_1_42, \
                         fg_1_43, fg_1_44, fh_1_52, fh_1_53, fh_1_54, fh_1_55, \
                         fh_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * fg_1_40[k]
                      + fh_1_52[k];

            t_41[k] = -ab_x[k] * fg_1_41[k]
                      + fh_1_53[k];

            t_42[k] = -ab_x[k] * fg_1_42[k]
                      + fh_1_54[k];

            t_43[k] = -ab_x[k] * fg_1_43[k]
                      + fh_1_55[k];

            t_44[k] = -ab_x[k] * fg_1_44[k]
                      + fh_1_56[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, fg_1_45, fg_1_46, fg_1_47, \
                         fg_1_48, fg_1_49, fh_1_63, fh_1_64, fh_1_65, fh_1_66, \
                         fh_1_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * fg_1_45[k]
                      + fh_1_63[k];

            t_46[k] = -ab_x[k] * fg_1_46[k]
                      + fh_1_64[k];

            t_47[k] = -ab_x[k] * fg_1_47[k]
                      + fh_1_65[k];

            t_48[k] = -ab_x[k] * fg_1_48[k]
                      + fh_1_66[k];

            t_49[k] = -ab_x[k] * fg_1_49[k]
                      + fh_1_67[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fg_1_50, fg_1_51, fg_1_52, \
                         fg_1_53, fg_1_54, fh_1_68, fh_1_69, fh_1_70, fh_1_71, \
                         fh_1_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * fg_1_50[k]
                      + fh_1_68[k];

            t_51[k] = -ab_x[k] * fg_1_51[k]
                      + fh_1_69[k];

            t_52[k] = -ab_x[k] * fg_1_52[k]
                      + fh_1_70[k];

            t_53[k] = -ab_x[k] * fg_1_53[k]
                      + fh_1_71[k];

            t_54[k] = -ab_x[k] * fg_1_54[k]
                      + fh_1_72[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, fg_1_55, fg_1_56, fg_1_57, \
                         fg_1_58, fg_1_59, fh_1_73, fh_1_74, fh_1_75, fh_1_76, \
                         fh_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * fg_1_55[k]
                      + fh_1_73[k];

            t_56[k] = -ab_x[k] * fg_1_56[k]
                      + fh_1_74[k];

            t_57[k] = -ab_x[k] * fg_1_57[k]
                      + fh_1_75[k];

            t_58[k] = -ab_x[k] * fg_1_58[k]
                      + fh_1_76[k];

            t_59[k] = -ab_x[k] * fg_1_59[k]
                      + fh_1_77[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fg_1_60, fg_1_61, fg_1_62, \
                         fg_1_63, fg_1_64, fh_1_84, fh_1_85, fh_1_86, fh_1_87, \
                         fh_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * fg_1_60[k]
                      + fh_1_84[k];

            t_61[k] = -ab_x[k] * fg_1_61[k]
                      + fh_1_85[k];

            t_62[k] = -ab_x[k] * fg_1_62[k]
                      + fh_1_86[k];

            t_63[k] = -ab_x[k] * fg_1_63[k]
                      + fh_1_87[k];

            t_64[k] = -ab_x[k] * fg_1_64[k]
                      + fh_1_88[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, fg_1_65, fg_1_66, fg_1_67, \
                         fg_1_68, fg_1_69, fh_1_89, fh_1_90, fh_1_91, fh_1_92, \
                         fh_1_93 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * fg_1_65[k]
                      + fh_1_89[k];

            t_66[k] = -ab_x[k] * fg_1_66[k]
                      + fh_1_90[k];

            t_67[k] = -ab_x[k] * fg_1_67[k]
                      + fh_1_91[k];

            t_68[k] = -ab_x[k] * fg_1_68[k]
                      + fh_1_92[k];

            t_69[k] = -ab_x[k] * fg_1_69[k]
                      + fh_1_93[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fg_1_70, fg_1_71, fg_1_72, \
                         fg_1_73, fg_1_74, fh_1_94, fh_1_95, fh_1_96, fh_1_97, \
                         fh_1_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * fg_1_70[k]
                      + fh_1_94[k];

            t_71[k] = -ab_x[k] * fg_1_71[k]
                      + fh_1_95[k];

            t_72[k] = -ab_x[k] * fg_1_72[k]
                      + fh_1_96[k];

            t_73[k] = -ab_x[k] * fg_1_73[k]
                      + fh_1_97[k];

            t_74[k] = -ab_x[k] * fg_1_74[k]
                      + fh_1_98[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, fg_1_75, fg_1_76, fg_1_77, \
                         fg_1_78, fg_1_79, fh_1_105, fh_1_106, fh_1_107, fh_1_108, \
                         fh_1_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * fg_1_75[k]
                      + fh_1_105[k];

            t_76[k] = -ab_x[k] * fg_1_76[k]
                      + fh_1_106[k];

            t_77[k] = -ab_x[k] * fg_1_77[k]
                      + fh_1_107[k];

            t_78[k] = -ab_x[k] * fg_1_78[k]
                      + fh_1_108[k];

            t_79[k] = -ab_x[k] * fg_1_79[k]
                      + fh_1_109[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fg_1_80, fg_1_81, fg_1_82, \
                         fg_1_83, fg_1_84, fh_1_110, fh_1_111, fh_1_112, fh_1_113, \
                         fh_1_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * fg_1_80[k]
                      + fh_1_110[k];

            t_81[k] = -ab_x[k] * fg_1_81[k]
                      + fh_1_111[k];

            t_82[k] = -ab_x[k] * fg_1_82[k]
                      + fh_1_112[k];

            t_83[k] = -ab_x[k] * fg_1_83[k]
                      + fh_1_113[k];

            t_84[k] = -ab_x[k] * fg_1_84[k]
                      + fh_1_114[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, fg_1_85, fg_1_86, fg_1_87, \
                         fg_1_88, fg_1_89, fh_1_115, fh_1_116, fh_1_117, fh_1_118, \
                         fh_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * fg_1_85[k]
                      + fh_1_115[k];

            t_86[k] = -ab_x[k] * fg_1_86[k]
                      + fh_1_116[k];

            t_87[k] = -ab_x[k] * fg_1_87[k]
                      + fh_1_117[k];

            t_88[k] = -ab_x[k] * fg_1_88[k]
                      + fh_1_118[k];

            t_89[k] = -ab_x[k] * fg_1_89[k]
                      + fh_1_119[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fg_1_90, fg_1_91, fg_1_92, \
                         fg_1_93, fg_1_94, fh_1_126, fh_1_127, fh_1_128, fh_1_129, \
                         fh_1_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * fg_1_90[k]
                      + fh_1_126[k];

            t_91[k] = -ab_x[k] * fg_1_91[k]
                      + fh_1_127[k];

            t_92[k] = -ab_x[k] * fg_1_92[k]
                      + fh_1_128[k];

            t_93[k] = -ab_x[k] * fg_1_93[k]
                      + fh_1_129[k];

            t_94[k] = -ab_x[k] * fg_1_94[k]
                      + fh_1_130[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, fg_1_95, fg_1_96, fg_1_97, \
                         fg_1_98, fg_1_99, fh_1_131, fh_1_132, fh_1_133, fh_1_134, \
                         fh_1_135 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * fg_1_95[k]
                      + fh_1_131[k];

            t_96[k] = -ab_x[k] * fg_1_96[k]
                      + fh_1_132[k];

            t_97[k] = -ab_x[k] * fg_1_97[k]
                      + fh_1_133[k];

            t_98[k] = -ab_x[k] * fg_1_98[k]
                      + fh_1_134[k];

            t_99[k] = -ab_x[k] * fg_1_99[k]
                      + fh_1_135[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, fg_1_100, fg_1_101, \
                         fg_1_102, fg_1_103, fg_1_104, fh_1_136, fh_1_137, fh_1_138, fh_1_139, \
                         fh_1_140 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * fg_1_100[k]
                       + fh_1_136[k];

            t_101[k] = -ab_x[k] * fg_1_101[k]
                       + fh_1_137[k];

            t_102[k] = -ab_x[k] * fg_1_102[k]
                       + fh_1_138[k];

            t_103[k] = -ab_x[k] * fg_1_103[k]
                       + fh_1_139[k];

            t_104[k] = -ab_x[k] * fg_1_104[k]
                       + fh_1_140[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, fg_1_105, fg_1_106, \
                         fg_1_107, fg_1_108, fg_1_109, fh_1_147, fh_1_148, fh_1_149, fh_1_150, \
                         fh_1_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * fg_1_105[k]
                       + fh_1_147[k];

            t_106[k] = -ab_x[k] * fg_1_106[k]
                       + fh_1_148[k];

            t_107[k] = -ab_x[k] * fg_1_107[k]
                       + fh_1_149[k];

            t_108[k] = -ab_x[k] * fg_1_108[k]
                       + fh_1_150[k];

            t_109[k] = -ab_x[k] * fg_1_109[k]
                       + fh_1_151[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, fg_1_110, fg_1_111, \
                         fg_1_112, fg_1_113, fg_1_114, fh_1_152, fh_1_153, fh_1_154, fh_1_155, \
                         fh_1_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * fg_1_110[k]
                       + fh_1_152[k];

            t_111[k] = -ab_x[k] * fg_1_111[k]
                       + fh_1_153[k];

            t_112[k] = -ab_x[k] * fg_1_112[k]
                       + fh_1_154[k];

            t_113[k] = -ab_x[k] * fg_1_113[k]
                       + fh_1_155[k];

            t_114[k] = -ab_x[k] * fg_1_114[k]
                       + fh_1_156[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, fg_1_115, fg_1_116, \
                         fg_1_117, fg_1_118, fg_1_119, fh_1_157, fh_1_158, fh_1_159, fh_1_160, \
                         fh_1_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * fg_1_115[k]
                       + fh_1_157[k];

            t_116[k] = -ab_x[k] * fg_1_116[k]
                       + fh_1_158[k];

            t_117[k] = -ab_x[k] * fg_1_117[k]
                       + fh_1_159[k];

            t_118[k] = -ab_x[k] * fg_1_118[k]
                       + fh_1_160[k];

            t_119[k] = -ab_x[k] * fg_1_119[k]
                       + fh_1_161[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, fg_1_120, fg_1_121, \
                         fg_1_122, fg_1_123, fg_1_124, fh_1_168, fh_1_169, fh_1_170, fh_1_171, \
                         fh_1_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * fg_1_120[k]
                       + fh_1_168[k];

            t_121[k] = -ab_x[k] * fg_1_121[k]
                       + fh_1_169[k];

            t_122[k] = -ab_x[k] * fg_1_122[k]
                       + fh_1_170[k];

            t_123[k] = -ab_x[k] * fg_1_123[k]
                       + fh_1_171[k];

            t_124[k] = -ab_x[k] * fg_1_124[k]
                       + fh_1_172[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, fg_1_125, fg_1_126, \
                         fg_1_127, fg_1_128, fg_1_129, fh_1_173, fh_1_174, fh_1_175, fh_1_176, \
                         fh_1_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * fg_1_125[k]
                       + fh_1_173[k];

            t_126[k] = -ab_x[k] * fg_1_126[k]
                       + fh_1_174[k];

            t_127[k] = -ab_x[k] * fg_1_127[k]
                       + fh_1_175[k];

            t_128[k] = -ab_x[k] * fg_1_128[k]
                       + fh_1_176[k];

            t_129[k] = -ab_x[k] * fg_1_129[k]
                       + fh_1_177[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, fg_1_130, fg_1_131, \
                         fg_1_132, fg_1_133, fg_1_134, fh_1_178, fh_1_179, fh_1_180, fh_1_181, \
                         fh_1_182 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * fg_1_130[k]
                       + fh_1_178[k];

            t_131[k] = -ab_x[k] * fg_1_131[k]
                       + fh_1_179[k];

            t_132[k] = -ab_x[k] * fg_1_132[k]
                       + fh_1_180[k];

            t_133[k] = -ab_x[k] * fg_1_133[k]
                       + fh_1_181[k];

            t_134[k] = -ab_x[k] * fg_1_134[k]
                       + fh_1_182[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, fg_1_135, fg_1_136, \
                         fg_1_137, fg_1_138, fg_1_139, fh_1_189, fh_1_190, fh_1_191, fh_1_192, \
                         fh_1_193 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * fg_1_135[k]
                       + fh_1_189[k];

            t_136[k] = -ab_x[k] * fg_1_136[k]
                       + fh_1_190[k];

            t_137[k] = -ab_x[k] * fg_1_137[k]
                       + fh_1_191[k];

            t_138[k] = -ab_x[k] * fg_1_138[k]
                       + fh_1_192[k];

            t_139[k] = -ab_x[k] * fg_1_139[k]
                       + fh_1_193[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, fg_1_140, fg_1_141, \
                         fg_1_142, fg_1_143, fg_1_144, fh_1_194, fh_1_195, fh_1_196, fh_1_197, \
                         fh_1_198 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * fg_1_140[k]
                       + fh_1_194[k];

            t_141[k] = -ab_x[k] * fg_1_141[k]
                       + fh_1_195[k];

            t_142[k] = -ab_x[k] * fg_1_142[k]
                       + fh_1_196[k];

            t_143[k] = -ab_x[k] * fg_1_143[k]
                       + fh_1_197[k];

            t_144[k] = -ab_x[k] * fg_1_144[k]
                       + fh_1_198[k];
        }
    }
}

static auto
compute_hrr_geom_010y_gg_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t fg_1, const size_t fg_0,
                                const size_t fh_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fg_1_90 = buffer.data(fg_1 + 90 * ncomps + c);
        const auto *fg_1_91 = buffer.data(fg_1 + 91 * ncomps + c);
        const auto *fg_1_92 = buffer.data(fg_1 + 92 * ncomps + c);
        const auto *fg_1_93 = buffer.data(fg_1 + 93 * ncomps + c);
        const auto *fg_1_94 = buffer.data(fg_1 + 94 * ncomps + c);
        const auto *fg_1_95 = buffer.data(fg_1 + 95 * ncomps + c);
        const auto *fg_1_96 = buffer.data(fg_1 + 96 * ncomps + c);
        const auto *fg_1_97 = buffer.data(fg_1 + 97 * ncomps + c);
        const auto *fg_1_98 = buffer.data(fg_1 + 98 * ncomps + c);
        const auto *fg_1_99 = buffer.data(fg_1 + 99 * ncomps + c);
        const auto *fg_1_100 = buffer.data(fg_1 + 100 * ncomps + c);
        const auto *fg_1_101 = buffer.data(fg_1 + 101 * ncomps + c);
        const auto *fg_1_102 = buffer.data(fg_1 + 102 * ncomps + c);
        const auto *fg_1_103 = buffer.data(fg_1 + 103 * ncomps + c);
        const auto *fg_1_104 = buffer.data(fg_1 + 104 * ncomps + c);
        const auto *fg_1_105 = buffer.data(fg_1 + 105 * ncomps + c);
        const auto *fg_1_106 = buffer.data(fg_1 + 106 * ncomps + c);
        const auto *fg_1_107 = buffer.data(fg_1 + 107 * ncomps + c);
        const auto *fg_1_108 = buffer.data(fg_1 + 108 * ncomps + c);
        const auto *fg_1_109 = buffer.data(fg_1 + 109 * ncomps + c);
        const auto *fg_1_110 = buffer.data(fg_1 + 110 * ncomps + c);
        const auto *fg_1_111 = buffer.data(fg_1 + 111 * ncomps + c);
        const auto *fg_1_112 = buffer.data(fg_1 + 112 * ncomps + c);
        const auto *fg_1_113 = buffer.data(fg_1 + 113 * ncomps + c);
        const auto *fg_1_114 = buffer.data(fg_1 + 114 * ncomps + c);
        const auto *fg_1_115 = buffer.data(fg_1 + 115 * ncomps + c);
        const auto *fg_1_116 = buffer.data(fg_1 + 116 * ncomps + c);
        const auto *fg_1_117 = buffer.data(fg_1 + 117 * ncomps + c);
        const auto *fg_1_118 = buffer.data(fg_1 + 118 * ncomps + c);
        const auto *fg_1_119 = buffer.data(fg_1 + 119 * ncomps + c);
        const auto *fg_1_120 = buffer.data(fg_1 + 120 * ncomps + c);
        const auto *fg_1_121 = buffer.data(fg_1 + 121 * ncomps + c);
        const auto *fg_1_122 = buffer.data(fg_1 + 122 * ncomps + c);
        const auto *fg_1_123 = buffer.data(fg_1 + 123 * ncomps + c);
        const auto *fg_1_124 = buffer.data(fg_1 + 124 * ncomps + c);
        const auto *fg_1_125 = buffer.data(fg_1 + 125 * ncomps + c);
        const auto *fg_1_126 = buffer.data(fg_1 + 126 * ncomps + c);
        const auto *fg_1_127 = buffer.data(fg_1 + 127 * ncomps + c);
        const auto *fg_1_128 = buffer.data(fg_1 + 128 * ncomps + c);
        const auto *fg_1_129 = buffer.data(fg_1 + 129 * ncomps + c);
        const auto *fg_1_130 = buffer.data(fg_1 + 130 * ncomps + c);
        const auto *fg_1_131 = buffer.data(fg_1 + 131 * ncomps + c);
        const auto *fg_1_132 = buffer.data(fg_1 + 132 * ncomps + c);
        const auto *fg_1_133 = buffer.data(fg_1 + 133 * ncomps + c);
        const auto *fg_1_134 = buffer.data(fg_1 + 134 * ncomps + c);
        const auto *fg_1_135 = buffer.data(fg_1 + 135 * ncomps + c);
        const auto *fg_1_136 = buffer.data(fg_1 + 136 * ncomps + c);
        const auto *fg_1_137 = buffer.data(fg_1 + 137 * ncomps + c);
        const auto *fg_1_138 = buffer.data(fg_1 + 138 * ncomps + c);
        const auto *fg_1_139 = buffer.data(fg_1 + 139 * ncomps + c);
        const auto *fg_1_140 = buffer.data(fg_1 + 140 * ncomps + c);
        const auto *fg_1_141 = buffer.data(fg_1 + 141 * ncomps + c);
        const auto *fg_1_142 = buffer.data(fg_1 + 142 * ncomps + c);
        const auto *fg_1_143 = buffer.data(fg_1 + 143 * ncomps + c);
        const auto *fg_1_144 = buffer.data(fg_1 + 144 * ncomps + c);
        const auto *fg_1_145 = buffer.data(fg_1 + 145 * ncomps + c);
        const auto *fg_1_146 = buffer.data(fg_1 + 146 * ncomps + c);
        const auto *fg_1_147 = buffer.data(fg_1 + 147 * ncomps + c);
        const auto *fg_1_148 = buffer.data(fg_1 + 148 * ncomps + c);
        const auto *fg_1_149 = buffer.data(fg_1 + 149 * ncomps + c);

        const auto *fg_0_90 = buffer.data(fg_0 + 90 * ncomps + c);
        const auto *fg_0_91 = buffer.data(fg_0 + 91 * ncomps + c);
        const auto *fg_0_92 = buffer.data(fg_0 + 92 * ncomps + c);
        const auto *fg_0_93 = buffer.data(fg_0 + 93 * ncomps + c);
        const auto *fg_0_94 = buffer.data(fg_0 + 94 * ncomps + c);
        const auto *fg_0_95 = buffer.data(fg_0 + 95 * ncomps + c);
        const auto *fg_0_96 = buffer.data(fg_0 + 96 * ncomps + c);
        const auto *fg_0_97 = buffer.data(fg_0 + 97 * ncomps + c);
        const auto *fg_0_98 = buffer.data(fg_0 + 98 * ncomps + c);
        const auto *fg_0_99 = buffer.data(fg_0 + 99 * ncomps + c);
        const auto *fg_0_100 = buffer.data(fg_0 + 100 * ncomps + c);
        const auto *fg_0_101 = buffer.data(fg_0 + 101 * ncomps + c);
        const auto *fg_0_102 = buffer.data(fg_0 + 102 * ncomps + c);
        const auto *fg_0_103 = buffer.data(fg_0 + 103 * ncomps + c);
        const auto *fg_0_104 = buffer.data(fg_0 + 104 * ncomps + c);
        const auto *fg_0_105 = buffer.data(fg_0 + 105 * ncomps + c);
        const auto *fg_0_106 = buffer.data(fg_0 + 106 * ncomps + c);
        const auto *fg_0_107 = buffer.data(fg_0 + 107 * ncomps + c);
        const auto *fg_0_108 = buffer.data(fg_0 + 108 * ncomps + c);
        const auto *fg_0_109 = buffer.data(fg_0 + 109 * ncomps + c);
        const auto *fg_0_110 = buffer.data(fg_0 + 110 * ncomps + c);
        const auto *fg_0_111 = buffer.data(fg_0 + 111 * ncomps + c);
        const auto *fg_0_112 = buffer.data(fg_0 + 112 * ncomps + c);
        const auto *fg_0_113 = buffer.data(fg_0 + 113 * ncomps + c);
        const auto *fg_0_114 = buffer.data(fg_0 + 114 * ncomps + c);
        const auto *fg_0_115 = buffer.data(fg_0 + 115 * ncomps + c);
        const auto *fg_0_116 = buffer.data(fg_0 + 116 * ncomps + c);
        const auto *fg_0_117 = buffer.data(fg_0 + 117 * ncomps + c);
        const auto *fg_0_118 = buffer.data(fg_0 + 118 * ncomps + c);
        const auto *fg_0_119 = buffer.data(fg_0 + 119 * ncomps + c);
        const auto *fg_0_120 = buffer.data(fg_0 + 120 * ncomps + c);
        const auto *fg_0_121 = buffer.data(fg_0 + 121 * ncomps + c);
        const auto *fg_0_122 = buffer.data(fg_0 + 122 * ncomps + c);
        const auto *fg_0_123 = buffer.data(fg_0 + 123 * ncomps + c);
        const auto *fg_0_124 = buffer.data(fg_0 + 124 * ncomps + c);
        const auto *fg_0_125 = buffer.data(fg_0 + 125 * ncomps + c);
        const auto *fg_0_126 = buffer.data(fg_0 + 126 * ncomps + c);
        const auto *fg_0_127 = buffer.data(fg_0 + 127 * ncomps + c);
        const auto *fg_0_128 = buffer.data(fg_0 + 128 * ncomps + c);
        const auto *fg_0_129 = buffer.data(fg_0 + 129 * ncomps + c);
        const auto *fg_0_130 = buffer.data(fg_0 + 130 * ncomps + c);
        const auto *fg_0_131 = buffer.data(fg_0 + 131 * ncomps + c);
        const auto *fg_0_132 = buffer.data(fg_0 + 132 * ncomps + c);
        const auto *fg_0_133 = buffer.data(fg_0 + 133 * ncomps + c);
        const auto *fg_0_134 = buffer.data(fg_0 + 134 * ncomps + c);
        const auto *fg_0_135 = buffer.data(fg_0 + 135 * ncomps + c);
        const auto *fg_0_136 = buffer.data(fg_0 + 136 * ncomps + c);
        const auto *fg_0_137 = buffer.data(fg_0 + 137 * ncomps + c);
        const auto *fg_0_138 = buffer.data(fg_0 + 138 * ncomps + c);
        const auto *fg_0_139 = buffer.data(fg_0 + 139 * ncomps + c);
        const auto *fg_0_140 = buffer.data(fg_0 + 140 * ncomps + c);
        const auto *fg_0_141 = buffer.data(fg_0 + 141 * ncomps + c);
        const auto *fg_0_142 = buffer.data(fg_0 + 142 * ncomps + c);
        const auto *fg_0_143 = buffer.data(fg_0 + 143 * ncomps + c);
        const auto *fg_0_144 = buffer.data(fg_0 + 144 * ncomps + c);
        const auto *fg_0_145 = buffer.data(fg_0 + 145 * ncomps + c);
        const auto *fg_0_146 = buffer.data(fg_0 + 146 * ncomps + c);
        const auto *fg_0_147 = buffer.data(fg_0 + 147 * ncomps + c);
        const auto *fg_0_148 = buffer.data(fg_0 + 148 * ncomps + c);
        const auto *fg_0_149 = buffer.data(fg_0 + 149 * ncomps + c);

        const auto *fh_1_127 = buffer.data(fh_1 + 127 * ncomps + c);
        const auto *fh_1_129 = buffer.data(fh_1 + 129 * ncomps + c);
        const auto *fh_1_130 = buffer.data(fh_1 + 130 * ncomps + c);
        const auto *fh_1_132 = buffer.data(fh_1 + 132 * ncomps + c);
        const auto *fh_1_133 = buffer.data(fh_1 + 133 * ncomps + c);
        const auto *fh_1_134 = buffer.data(fh_1 + 134 * ncomps + c);
        const auto *fh_1_136 = buffer.data(fh_1 + 136 * ncomps + c);
        const auto *fh_1_137 = buffer.data(fh_1 + 137 * ncomps + c);
        const auto *fh_1_138 = buffer.data(fh_1 + 138 * ncomps + c);
        const auto *fh_1_139 = buffer.data(fh_1 + 139 * ncomps + c);
        const auto *fh_1_141 = buffer.data(fh_1 + 141 * ncomps + c);
        const auto *fh_1_142 = buffer.data(fh_1 + 142 * ncomps + c);
        const auto *fh_1_143 = buffer.data(fh_1 + 143 * ncomps + c);
        const auto *fh_1_144 = buffer.data(fh_1 + 144 * ncomps + c);
        const auto *fh_1_145 = buffer.data(fh_1 + 145 * ncomps + c);
        const auto *fh_1_148 = buffer.data(fh_1 + 148 * ncomps + c);
        const auto *fh_1_150 = buffer.data(fh_1 + 150 * ncomps + c);
        const auto *fh_1_151 = buffer.data(fh_1 + 151 * ncomps + c);
        const auto *fh_1_153 = buffer.data(fh_1 + 153 * ncomps + c);
        const auto *fh_1_154 = buffer.data(fh_1 + 154 * ncomps + c);
        const auto *fh_1_155 = buffer.data(fh_1 + 155 * ncomps + c);
        const auto *fh_1_157 = buffer.data(fh_1 + 157 * ncomps + c);
        const auto *fh_1_158 = buffer.data(fh_1 + 158 * ncomps + c);
        const auto *fh_1_159 = buffer.data(fh_1 + 159 * ncomps + c);
        const auto *fh_1_160 = buffer.data(fh_1 + 160 * ncomps + c);
        const auto *fh_1_162 = buffer.data(fh_1 + 162 * ncomps + c);
        const auto *fh_1_163 = buffer.data(fh_1 + 163 * ncomps + c);
        const auto *fh_1_164 = buffer.data(fh_1 + 164 * ncomps + c);
        const auto *fh_1_165 = buffer.data(fh_1 + 165 * ncomps + c);
        const auto *fh_1_166 = buffer.data(fh_1 + 166 * ncomps + c);
        const auto *fh_1_169 = buffer.data(fh_1 + 169 * ncomps + c);
        const auto *fh_1_171 = buffer.data(fh_1 + 171 * ncomps + c);
        const auto *fh_1_172 = buffer.data(fh_1 + 172 * ncomps + c);
        const auto *fh_1_174 = buffer.data(fh_1 + 174 * ncomps + c);
        const auto *fh_1_175 = buffer.data(fh_1 + 175 * ncomps + c);
        const auto *fh_1_176 = buffer.data(fh_1 + 176 * ncomps + c);
        const auto *fh_1_178 = buffer.data(fh_1 + 178 * ncomps + c);
        const auto *fh_1_179 = buffer.data(fh_1 + 179 * ncomps + c);
        const auto *fh_1_180 = buffer.data(fh_1 + 180 * ncomps + c);
        const auto *fh_1_181 = buffer.data(fh_1 + 181 * ncomps + c);
        const auto *fh_1_183 = buffer.data(fh_1 + 183 * ncomps + c);
        const auto *fh_1_184 = buffer.data(fh_1 + 184 * ncomps + c);
        const auto *fh_1_185 = buffer.data(fh_1 + 185 * ncomps + c);
        const auto *fh_1_186 = buffer.data(fh_1 + 186 * ncomps + c);
        const auto *fh_1_187 = buffer.data(fh_1 + 187 * ncomps + c);
        const auto *fh_1_190 = buffer.data(fh_1 + 190 * ncomps + c);
        const auto *fh_1_191 = buffer.data(fh_1 + 191 * ncomps + c);
        const auto *fh_1_192 = buffer.data(fh_1 + 192 * ncomps + c);
        const auto *fh_1_193 = buffer.data(fh_1 + 193 * ncomps + c);
        const auto *fh_1_194 = buffer.data(fh_1 + 194 * ncomps + c);
        const auto *fh_1_195 = buffer.data(fh_1 + 195 * ncomps + c);
        const auto *fh_1_196 = buffer.data(fh_1 + 196 * ncomps + c);
        const auto *fh_1_197 = buffer.data(fh_1 + 197 * ncomps + c);
        const auto *fh_1_198 = buffer.data(fh_1 + 198 * ncomps + c);
        const auto *fh_1_199 = buffer.data(fh_1 + 199 * ncomps + c);
        const auto *fh_1_200 = buffer.data(fh_1 + 200 * ncomps + c);
        const auto *fh_1_201 = buffer.data(fh_1 + 201 * ncomps + c);
        const auto *fh_1_202 = buffer.data(fh_1 + 202 * ncomps + c);
        const auto *fh_1_203 = buffer.data(fh_1 + 203 * ncomps + c);
        const auto *fh_1_204 = buffer.data(fh_1 + 204 * ncomps + c);
        const auto *fh_1_205 = buffer.data(fh_1 + 205 * ncomps + c);
        const auto *fh_1_206 = buffer.data(fh_1 + 206 * ncomps + c);
        const auto *fh_1_207 = buffer.data(fh_1 + 207 * ncomps + c);
        const auto *fh_1_208 = buffer.data(fh_1 + 208 * ncomps + c);
        const auto *fh_1_209 = buffer.data(fh_1 + 209 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, fg_1_145, fg_1_146, \
                         fg_1_147, fg_1_148, fg_1_149, fh_1_199, fh_1_200, fh_1_201, fh_1_202, \
                         fh_1_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * fg_1_145[k]
                       + fh_1_199[k];

            t_146[k] = -ab_x[k] * fg_1_146[k]
                       + fh_1_200[k];

            t_147[k] = -ab_x[k] * fg_1_147[k]
                       + fh_1_201[k];

            t_148[k] = -ab_x[k] * fg_1_148[k]
                       + fh_1_202[k];

            t_149[k] = -ab_x[k] * fg_1_149[k]
                       + fh_1_203[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, ab_y, fg_1_90, fg_1_91, fg_1_92, fg_0_90, \
                         fg_0_91, fg_0_92, fh_1_127, fh_1_129, \
                         fh_1_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_y[k] * fg_1_90[k]
                       + fg_0_90[k]
                       + fh_1_127[k];

            t_151[k] = -ab_y[k] * fg_1_91[k]
                       + fg_0_91[k]
                       + fh_1_129[k];

            t_152[k] = -ab_y[k] * fg_1_92[k]
                       + fg_0_92[k]
                       + fh_1_130[k];
        }

#pragma omp simd aligned(t_153, t_154, t_155, ab_y, fg_1_93, fg_1_94, fg_1_95, fg_0_93, \
                         fg_0_94, fg_0_95, fh_1_132, fh_1_133, \
                         fh_1_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_153[k] = -ab_y[k] * fg_1_93[k]
                       + fg_0_93[k]
                       + fh_1_132[k];

            t_154[k] = -ab_y[k] * fg_1_94[k]
                       + fg_0_94[k]
                       + fh_1_133[k];

            t_155[k] = -ab_y[k] * fg_1_95[k]
                       + fg_0_95[k]
                       + fh_1_134[k];
        }

#pragma omp simd aligned(t_156, t_157, t_158, ab_y, fg_1_96, fg_1_97, fg_1_98, fg_0_96, \
                         fg_0_97, fg_0_98, fh_1_136, fh_1_137, \
                         fh_1_138 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = -ab_y[k] * fg_1_96[k]
                       + fg_0_96[k]
                       + fh_1_136[k];

            t_157[k] = -ab_y[k] * fg_1_97[k]
                       + fg_0_97[k]
                       + fh_1_137[k];

            t_158[k] = -ab_y[k] * fg_1_98[k]
                       + fg_0_98[k]
                       + fh_1_138[k];
        }

#pragma omp simd aligned(t_159, t_160, t_161, ab_y, fg_1_99, fg_1_100, fg_1_101, fg_0_99, \
                         fg_0_100, fg_0_101, fh_1_139, fh_1_141, \
                         fh_1_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_159[k] = -ab_y[k] * fg_1_99[k]
                       + fg_0_99[k]
                       + fh_1_139[k];

            t_160[k] = -ab_y[k] * fg_1_100[k]
                       + fg_0_100[k]
                       + fh_1_141[k];

            t_161[k] = -ab_y[k] * fg_1_101[k]
                       + fg_0_101[k]
                       + fh_1_142[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, ab_y, fg_1_102, fg_1_103, fg_1_104, fg_0_102, \
                         fg_0_103, fg_0_104, fh_1_143, fh_1_144, \
                         fh_1_145 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = -ab_y[k] * fg_1_102[k]
                       + fg_0_102[k]
                       + fh_1_143[k];

            t_163[k] = -ab_y[k] * fg_1_103[k]
                       + fg_0_103[k]
                       + fh_1_144[k];

            t_164[k] = -ab_y[k] * fg_1_104[k]
                       + fg_0_104[k]
                       + fh_1_145[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, ab_y, fg_1_105, fg_1_106, fg_1_107, fg_0_105, \
                         fg_0_106, fg_0_107, fh_1_148, fh_1_150, \
                         fh_1_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_y[k] * fg_1_105[k]
                       + fg_0_105[k]
                       + fh_1_148[k];

            t_166[k] = -ab_y[k] * fg_1_106[k]
                       + fg_0_106[k]
                       + fh_1_150[k];

            t_167[k] = -ab_y[k] * fg_1_107[k]
                       + fg_0_107[k]
                       + fh_1_151[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, ab_y, fg_1_108, fg_1_109, fg_1_110, fg_0_108, \
                         fg_0_109, fg_0_110, fh_1_153, fh_1_154, \
                         fh_1_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = -ab_y[k] * fg_1_108[k]
                       + fg_0_108[k]
                       + fh_1_153[k];

            t_169[k] = -ab_y[k] * fg_1_109[k]
                       + fg_0_109[k]
                       + fh_1_154[k];

            t_170[k] = -ab_y[k] * fg_1_110[k]
                       + fg_0_110[k]
                       + fh_1_155[k];
        }

#pragma omp simd aligned(t_171, t_172, t_173, ab_y, fg_1_111, fg_1_112, fg_1_113, fg_0_111, \
                         fg_0_112, fg_0_113, fh_1_157, fh_1_158, \
                         fh_1_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_171[k] = -ab_y[k] * fg_1_111[k]
                       + fg_0_111[k]
                       + fh_1_157[k];

            t_172[k] = -ab_y[k] * fg_1_112[k]
                       + fg_0_112[k]
                       + fh_1_158[k];

            t_173[k] = -ab_y[k] * fg_1_113[k]
                       + fg_0_113[k]
                       + fh_1_159[k];
        }

#pragma omp simd aligned(t_174, t_175, t_176, ab_y, fg_1_114, fg_1_115, fg_1_116, fg_0_114, \
                         fg_0_115, fg_0_116, fh_1_160, fh_1_162, \
                         fh_1_163 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_174[k] = -ab_y[k] * fg_1_114[k]
                       + fg_0_114[k]
                       + fh_1_160[k];

            t_175[k] = -ab_y[k] * fg_1_115[k]
                       + fg_0_115[k]
                       + fh_1_162[k];

            t_176[k] = -ab_y[k] * fg_1_116[k]
                       + fg_0_116[k]
                       + fh_1_163[k];
        }

#pragma omp simd aligned(t_177, t_178, t_179, ab_y, fg_1_117, fg_1_118, fg_1_119, fg_0_117, \
                         fg_0_118, fg_0_119, fh_1_164, fh_1_165, \
                         fh_1_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_177[k] = -ab_y[k] * fg_1_117[k]
                       + fg_0_117[k]
                       + fh_1_164[k];

            t_178[k] = -ab_y[k] * fg_1_118[k]
                       + fg_0_118[k]
                       + fh_1_165[k];

            t_179[k] = -ab_y[k] * fg_1_119[k]
                       + fg_0_119[k]
                       + fh_1_166[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, ab_y, fg_1_120, fg_1_121, fg_1_122, fg_0_120, \
                         fg_0_121, fg_0_122, fh_1_169, fh_1_171, \
                         fh_1_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_y[k] * fg_1_120[k]
                       + fg_0_120[k]
                       + fh_1_169[k];

            t_181[k] = -ab_y[k] * fg_1_121[k]
                       + fg_0_121[k]
                       + fh_1_171[k];

            t_182[k] = -ab_y[k] * fg_1_122[k]
                       + fg_0_122[k]
                       + fh_1_172[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, ab_y, fg_1_123, fg_1_124, fg_1_125, fg_0_123, \
                         fg_0_124, fg_0_125, fh_1_174, fh_1_175, \
                         fh_1_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = -ab_y[k] * fg_1_123[k]
                       + fg_0_123[k]
                       + fh_1_174[k];

            t_184[k] = -ab_y[k] * fg_1_124[k]
                       + fg_0_124[k]
                       + fh_1_175[k];

            t_185[k] = -ab_y[k] * fg_1_125[k]
                       + fg_0_125[k]
                       + fh_1_176[k];
        }

#pragma omp simd aligned(t_186, t_187, t_188, ab_y, fg_1_126, fg_1_127, fg_1_128, fg_0_126, \
                         fg_0_127, fg_0_128, fh_1_178, fh_1_179, \
                         fh_1_180 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_186[k] = -ab_y[k] * fg_1_126[k]
                       + fg_0_126[k]
                       + fh_1_178[k];

            t_187[k] = -ab_y[k] * fg_1_127[k]
                       + fg_0_127[k]
                       + fh_1_179[k];

            t_188[k] = -ab_y[k] * fg_1_128[k]
                       + fg_0_128[k]
                       + fh_1_180[k];
        }

#pragma omp simd aligned(t_189, t_190, t_191, ab_y, fg_1_129, fg_1_130, fg_1_131, fg_0_129, \
                         fg_0_130, fg_0_131, fh_1_181, fh_1_183, \
                         fh_1_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_189[k] = -ab_y[k] * fg_1_129[k]
                       + fg_0_129[k]
                       + fh_1_181[k];

            t_190[k] = -ab_y[k] * fg_1_130[k]
                       + fg_0_130[k]
                       + fh_1_183[k];

            t_191[k] = -ab_y[k] * fg_1_131[k]
                       + fg_0_131[k]
                       + fh_1_184[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, ab_y, fg_1_132, fg_1_133, fg_1_134, fg_0_132, \
                         fg_0_133, fg_0_134, fh_1_185, fh_1_186, \
                         fh_1_187 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = -ab_y[k] * fg_1_132[k]
                       + fg_0_132[k]
                       + fh_1_185[k];

            t_193[k] = -ab_y[k] * fg_1_133[k]
                       + fg_0_133[k]
                       + fh_1_186[k];

            t_194[k] = -ab_y[k] * fg_1_134[k]
                       + fg_0_134[k]
                       + fh_1_187[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, ab_y, fg_1_135, fg_1_136, fg_1_137, fg_0_135, \
                         fg_0_136, fg_0_137, fh_1_190, fh_1_192, \
                         fh_1_193 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_y[k] * fg_1_135[k]
                       + fg_0_135[k]
                       + fh_1_190[k];

            t_196[k] = -ab_y[k] * fg_1_136[k]
                       + fg_0_136[k]
                       + fh_1_192[k];

            t_197[k] = -ab_y[k] * fg_1_137[k]
                       + fg_0_137[k]
                       + fh_1_193[k];
        }

#pragma omp simd aligned(t_198, t_199, t_200, ab_y, fg_1_138, fg_1_139, fg_1_140, fg_0_138, \
                         fg_0_139, fg_0_140, fh_1_195, fh_1_196, \
                         fh_1_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_198[k] = -ab_y[k] * fg_1_138[k]
                       + fg_0_138[k]
                       + fh_1_195[k];

            t_199[k] = -ab_y[k] * fg_1_139[k]
                       + fg_0_139[k]
                       + fh_1_196[k];

            t_200[k] = -ab_y[k] * fg_1_140[k]
                       + fg_0_140[k]
                       + fh_1_197[k];
        }

#pragma omp simd aligned(t_201, t_202, t_203, ab_y, fg_1_141, fg_1_142, fg_1_143, fg_0_141, \
                         fg_0_142, fg_0_143, fh_1_199, fh_1_200, \
                         fh_1_201 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_201[k] = -ab_y[k] * fg_1_141[k]
                       + fg_0_141[k]
                       + fh_1_199[k];

            t_202[k] = -ab_y[k] * fg_1_142[k]
                       + fg_0_142[k]
                       + fh_1_200[k];

            t_203[k] = -ab_y[k] * fg_1_143[k]
                       + fg_0_143[k]
                       + fh_1_201[k];
        }

#pragma omp simd aligned(t_204, t_205, t_206, ab_y, fg_1_144, fg_1_145, fg_1_146, fg_0_144, \
                         fg_0_145, fg_0_146, fh_1_202, fh_1_204, \
                         fh_1_205 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_204[k] = -ab_y[k] * fg_1_144[k]
                       + fg_0_144[k]
                       + fh_1_202[k];

            t_205[k] = -ab_y[k] * fg_1_145[k]
                       + fg_0_145[k]
                       + fh_1_204[k];

            t_206[k] = -ab_y[k] * fg_1_146[k]
                       + fg_0_146[k]
                       + fh_1_205[k];
        }

#pragma omp simd aligned(t_207, t_208, t_209, ab_y, fg_1_147, fg_1_148, fg_1_149, fg_0_147, \
                         fg_0_148, fg_0_149, fh_1_206, fh_1_207, \
                         fh_1_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_207[k] = -ab_y[k] * fg_1_147[k]
                       + fg_0_147[k]
                       + fh_1_206[k];

            t_208[k] = -ab_y[k] * fg_1_148[k]
                       + fg_0_148[k]
                       + fh_1_207[k];

            t_209[k] = -ab_y[k] * fg_1_149[k]
                       + fg_0_149[k]
                       + fh_1_208[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_z, fg_1_135, fg_1_136, \
                         fg_1_137, fg_1_138, fg_1_139, fh_1_191, fh_1_193, fh_1_194, fh_1_196, \
                         fh_1_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_z[k] * fg_1_135[k]
                       + fh_1_191[k];

            t_211[k] = -ab_z[k] * fg_1_136[k]
                       + fh_1_193[k];

            t_212[k] = -ab_z[k] * fg_1_137[k]
                       + fh_1_194[k];

            t_213[k] = -ab_z[k] * fg_1_138[k]
                       + fh_1_196[k];

            t_214[k] = -ab_z[k] * fg_1_139[k]
                       + fh_1_197[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_z, fg_1_140, fg_1_141, \
                         fg_1_142, fg_1_143, fg_1_144, fh_1_198, fh_1_200, fh_1_201, fh_1_202, \
                         fh_1_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_z[k] * fg_1_140[k]
                       + fh_1_198[k];

            t_216[k] = -ab_z[k] * fg_1_141[k]
                       + fh_1_200[k];

            t_217[k] = -ab_z[k] * fg_1_142[k]
                       + fh_1_201[k];

            t_218[k] = -ab_z[k] * fg_1_143[k]
                       + fh_1_202[k];

            t_219[k] = -ab_z[k] * fg_1_144[k]
                       + fh_1_203[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_z, fg_1_145, fg_1_146, \
                         fg_1_147, fg_1_148, fg_1_149, fh_1_205, fh_1_206, fh_1_207, fh_1_208, \
                         fh_1_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_z[k] * fg_1_145[k]
                       + fh_1_205[k];

            t_221[k] = -ab_z[k] * fg_1_146[k]
                       + fh_1_206[k];

            t_222[k] = -ab_z[k] * fg_1_147[k]
                       + fh_1_207[k];

            t_223[k] = -ab_z[k] * fg_1_148[k]
                       + fh_1_208[k];

            t_224[k] = -ab_z[k] * fg_1_149[k]
                       + fh_1_209[k];
        }
    }
}

auto
compute_hrr_geom_010y_gg(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t fg_1, const size_t fg_0,
                         const size_t fh_1, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_geom_010y_gg_piece0(buffer, coordinates, target, fg_1, fh_1, ncomps, nmax);

    compute_hrr_geom_010y_gg_piece1(buffer, coordinates, target, fg_1, fg_0, fh_1, ncomps,
                                    nmax);
}

}  // namespace simdtrf
