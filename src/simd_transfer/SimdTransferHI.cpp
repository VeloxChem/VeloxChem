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


#include "SimdTransferHI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_hi_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gi, const size_t gk, const size_t ncomps,
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

        const auto *gi_0 = buffer.data(gi + 0 * ncomps + c);
        const auto *gi_1 = buffer.data(gi + 1 * ncomps + c);
        const auto *gi_2 = buffer.data(gi + 2 * ncomps + c);
        const auto *gi_3 = buffer.data(gi + 3 * ncomps + c);
        const auto *gi_4 = buffer.data(gi + 4 * ncomps + c);
        const auto *gi_5 = buffer.data(gi + 5 * ncomps + c);
        const auto *gi_6 = buffer.data(gi + 6 * ncomps + c);
        const auto *gi_7 = buffer.data(gi + 7 * ncomps + c);
        const auto *gi_8 = buffer.data(gi + 8 * ncomps + c);
        const auto *gi_9 = buffer.data(gi + 9 * ncomps + c);
        const auto *gi_10 = buffer.data(gi + 10 * ncomps + c);
        const auto *gi_11 = buffer.data(gi + 11 * ncomps + c);
        const auto *gi_12 = buffer.data(gi + 12 * ncomps + c);
        const auto *gi_13 = buffer.data(gi + 13 * ncomps + c);
        const auto *gi_14 = buffer.data(gi + 14 * ncomps + c);
        const auto *gi_15 = buffer.data(gi + 15 * ncomps + c);
        const auto *gi_16 = buffer.data(gi + 16 * ncomps + c);
        const auto *gi_17 = buffer.data(gi + 17 * ncomps + c);
        const auto *gi_18 = buffer.data(gi + 18 * ncomps + c);
        const auto *gi_19 = buffer.data(gi + 19 * ncomps + c);
        const auto *gi_20 = buffer.data(gi + 20 * ncomps + c);
        const auto *gi_21 = buffer.data(gi + 21 * ncomps + c);
        const auto *gi_22 = buffer.data(gi + 22 * ncomps + c);
        const auto *gi_23 = buffer.data(gi + 23 * ncomps + c);
        const auto *gi_24 = buffer.data(gi + 24 * ncomps + c);
        const auto *gi_25 = buffer.data(gi + 25 * ncomps + c);
        const auto *gi_26 = buffer.data(gi + 26 * ncomps + c);
        const auto *gi_27 = buffer.data(gi + 27 * ncomps + c);
        const auto *gi_28 = buffer.data(gi + 28 * ncomps + c);
        const auto *gi_29 = buffer.data(gi + 29 * ncomps + c);
        const auto *gi_30 = buffer.data(gi + 30 * ncomps + c);
        const auto *gi_31 = buffer.data(gi + 31 * ncomps + c);
        const auto *gi_32 = buffer.data(gi + 32 * ncomps + c);
        const auto *gi_33 = buffer.data(gi + 33 * ncomps + c);
        const auto *gi_34 = buffer.data(gi + 34 * ncomps + c);
        const auto *gi_35 = buffer.data(gi + 35 * ncomps + c);
        const auto *gi_36 = buffer.data(gi + 36 * ncomps + c);
        const auto *gi_37 = buffer.data(gi + 37 * ncomps + c);
        const auto *gi_38 = buffer.data(gi + 38 * ncomps + c);
        const auto *gi_39 = buffer.data(gi + 39 * ncomps + c);
        const auto *gi_40 = buffer.data(gi + 40 * ncomps + c);
        const auto *gi_41 = buffer.data(gi + 41 * ncomps + c);
        const auto *gi_42 = buffer.data(gi + 42 * ncomps + c);
        const auto *gi_43 = buffer.data(gi + 43 * ncomps + c);
        const auto *gi_44 = buffer.data(gi + 44 * ncomps + c);
        const auto *gi_45 = buffer.data(gi + 45 * ncomps + c);
        const auto *gi_46 = buffer.data(gi + 46 * ncomps + c);
        const auto *gi_47 = buffer.data(gi + 47 * ncomps + c);
        const auto *gi_48 = buffer.data(gi + 48 * ncomps + c);
        const auto *gi_49 = buffer.data(gi + 49 * ncomps + c);
        const auto *gi_50 = buffer.data(gi + 50 * ncomps + c);
        const auto *gi_51 = buffer.data(gi + 51 * ncomps + c);
        const auto *gi_52 = buffer.data(gi + 52 * ncomps + c);
        const auto *gi_53 = buffer.data(gi + 53 * ncomps + c);
        const auto *gi_54 = buffer.data(gi + 54 * ncomps + c);
        const auto *gi_55 = buffer.data(gi + 55 * ncomps + c);
        const auto *gi_56 = buffer.data(gi + 56 * ncomps + c);
        const auto *gi_57 = buffer.data(gi + 57 * ncomps + c);
        const auto *gi_58 = buffer.data(gi + 58 * ncomps + c);
        const auto *gi_59 = buffer.data(gi + 59 * ncomps + c);
        const auto *gi_60 = buffer.data(gi + 60 * ncomps + c);
        const auto *gi_61 = buffer.data(gi + 61 * ncomps + c);
        const auto *gi_62 = buffer.data(gi + 62 * ncomps + c);
        const auto *gi_63 = buffer.data(gi + 63 * ncomps + c);
        const auto *gi_64 = buffer.data(gi + 64 * ncomps + c);
        const auto *gi_65 = buffer.data(gi + 65 * ncomps + c);
        const auto *gi_66 = buffer.data(gi + 66 * ncomps + c);
        const auto *gi_67 = buffer.data(gi + 67 * ncomps + c);
        const auto *gi_68 = buffer.data(gi + 68 * ncomps + c);
        const auto *gi_69 = buffer.data(gi + 69 * ncomps + c);
        const auto *gi_70 = buffer.data(gi + 70 * ncomps + c);
        const auto *gi_71 = buffer.data(gi + 71 * ncomps + c);
        const auto *gi_72 = buffer.data(gi + 72 * ncomps + c);
        const auto *gi_73 = buffer.data(gi + 73 * ncomps + c);
        const auto *gi_74 = buffer.data(gi + 74 * ncomps + c);
        const auto *gi_75 = buffer.data(gi + 75 * ncomps + c);
        const auto *gi_76 = buffer.data(gi + 76 * ncomps + c);
        const auto *gi_77 = buffer.data(gi + 77 * ncomps + c);
        const auto *gi_78 = buffer.data(gi + 78 * ncomps + c);
        const auto *gi_79 = buffer.data(gi + 79 * ncomps + c);
        const auto *gi_80 = buffer.data(gi + 80 * ncomps + c);
        const auto *gi_81 = buffer.data(gi + 81 * ncomps + c);
        const auto *gi_82 = buffer.data(gi + 82 * ncomps + c);
        const auto *gi_83 = buffer.data(gi + 83 * ncomps + c);
        const auto *gi_84 = buffer.data(gi + 84 * ncomps + c);
        const auto *gi_85 = buffer.data(gi + 85 * ncomps + c);
        const auto *gi_86 = buffer.data(gi + 86 * ncomps + c);
        const auto *gi_87 = buffer.data(gi + 87 * ncomps + c);
        const auto *gi_88 = buffer.data(gi + 88 * ncomps + c);
        const auto *gi_89 = buffer.data(gi + 89 * ncomps + c);
        const auto *gi_90 = buffer.data(gi + 90 * ncomps + c);
        const auto *gi_91 = buffer.data(gi + 91 * ncomps + c);
        const auto *gi_92 = buffer.data(gi + 92 * ncomps + c);
        const auto *gi_93 = buffer.data(gi + 93 * ncomps + c);
        const auto *gi_94 = buffer.data(gi + 94 * ncomps + c);
        const auto *gi_95 = buffer.data(gi + 95 * ncomps + c);
        const auto *gi_96 = buffer.data(gi + 96 * ncomps + c);
        const auto *gi_97 = buffer.data(gi + 97 * ncomps + c);
        const auto *gi_98 = buffer.data(gi + 98 * ncomps + c);
        const auto *gi_99 = buffer.data(gi + 99 * ncomps + c);
        const auto *gi_100 = buffer.data(gi + 100 * ncomps + c);
        const auto *gi_101 = buffer.data(gi + 101 * ncomps + c);
        const auto *gi_102 = buffer.data(gi + 102 * ncomps + c);
        const auto *gi_103 = buffer.data(gi + 103 * ncomps + c);
        const auto *gi_104 = buffer.data(gi + 104 * ncomps + c);
        const auto *gi_105 = buffer.data(gi + 105 * ncomps + c);
        const auto *gi_106 = buffer.data(gi + 106 * ncomps + c);
        const auto *gi_107 = buffer.data(gi + 107 * ncomps + c);
        const auto *gi_108 = buffer.data(gi + 108 * ncomps + c);
        const auto *gi_109 = buffer.data(gi + 109 * ncomps + c);
        const auto *gi_110 = buffer.data(gi + 110 * ncomps + c);
        const auto *gi_111 = buffer.data(gi + 111 * ncomps + c);
        const auto *gi_112 = buffer.data(gi + 112 * ncomps + c);
        const auto *gi_113 = buffer.data(gi + 113 * ncomps + c);
        const auto *gi_114 = buffer.data(gi + 114 * ncomps + c);
        const auto *gi_115 = buffer.data(gi + 115 * ncomps + c);
        const auto *gi_116 = buffer.data(gi + 116 * ncomps + c);
        const auto *gi_117 = buffer.data(gi + 117 * ncomps + c);
        const auto *gi_118 = buffer.data(gi + 118 * ncomps + c);
        const auto *gi_119 = buffer.data(gi + 119 * ncomps + c);
        const auto *gi_120 = buffer.data(gi + 120 * ncomps + c);
        const auto *gi_121 = buffer.data(gi + 121 * ncomps + c);
        const auto *gi_122 = buffer.data(gi + 122 * ncomps + c);
        const auto *gi_123 = buffer.data(gi + 123 * ncomps + c);
        const auto *gi_124 = buffer.data(gi + 124 * ncomps + c);
        const auto *gi_125 = buffer.data(gi + 125 * ncomps + c);
        const auto *gi_126 = buffer.data(gi + 126 * ncomps + c);
        const auto *gi_127 = buffer.data(gi + 127 * ncomps + c);
        const auto *gi_128 = buffer.data(gi + 128 * ncomps + c);
        const auto *gi_129 = buffer.data(gi + 129 * ncomps + c);
        const auto *gi_130 = buffer.data(gi + 130 * ncomps + c);
        const auto *gi_131 = buffer.data(gi + 131 * ncomps + c);
        const auto *gi_132 = buffer.data(gi + 132 * ncomps + c);
        const auto *gi_133 = buffer.data(gi + 133 * ncomps + c);
        const auto *gi_134 = buffer.data(gi + 134 * ncomps + c);
        const auto *gi_135 = buffer.data(gi + 135 * ncomps + c);
        const auto *gi_136 = buffer.data(gi + 136 * ncomps + c);
        const auto *gi_137 = buffer.data(gi + 137 * ncomps + c);
        const auto *gi_138 = buffer.data(gi + 138 * ncomps + c);
        const auto *gi_139 = buffer.data(gi + 139 * ncomps + c);
        const auto *gi_140 = buffer.data(gi + 140 * ncomps + c);
        const auto *gi_141 = buffer.data(gi + 141 * ncomps + c);
        const auto *gi_142 = buffer.data(gi + 142 * ncomps + c);
        const auto *gi_143 = buffer.data(gi + 143 * ncomps + c);
        const auto *gi_144 = buffer.data(gi + 144 * ncomps + c);

        const auto *gk_0 = buffer.data(gk + 0 * ncomps + c);
        const auto *gk_1 = buffer.data(gk + 1 * ncomps + c);
        const auto *gk_2 = buffer.data(gk + 2 * ncomps + c);
        const auto *gk_3 = buffer.data(gk + 3 * ncomps + c);
        const auto *gk_4 = buffer.data(gk + 4 * ncomps + c);
        const auto *gk_5 = buffer.data(gk + 5 * ncomps + c);
        const auto *gk_6 = buffer.data(gk + 6 * ncomps + c);
        const auto *gk_7 = buffer.data(gk + 7 * ncomps + c);
        const auto *gk_8 = buffer.data(gk + 8 * ncomps + c);
        const auto *gk_9 = buffer.data(gk + 9 * ncomps + c);
        const auto *gk_10 = buffer.data(gk + 10 * ncomps + c);
        const auto *gk_11 = buffer.data(gk + 11 * ncomps + c);
        const auto *gk_12 = buffer.data(gk + 12 * ncomps + c);
        const auto *gk_13 = buffer.data(gk + 13 * ncomps + c);
        const auto *gk_14 = buffer.data(gk + 14 * ncomps + c);
        const auto *gk_15 = buffer.data(gk + 15 * ncomps + c);
        const auto *gk_16 = buffer.data(gk + 16 * ncomps + c);
        const auto *gk_17 = buffer.data(gk + 17 * ncomps + c);
        const auto *gk_18 = buffer.data(gk + 18 * ncomps + c);
        const auto *gk_19 = buffer.data(gk + 19 * ncomps + c);
        const auto *gk_20 = buffer.data(gk + 20 * ncomps + c);
        const auto *gk_21 = buffer.data(gk + 21 * ncomps + c);
        const auto *gk_22 = buffer.data(gk + 22 * ncomps + c);
        const auto *gk_23 = buffer.data(gk + 23 * ncomps + c);
        const auto *gk_24 = buffer.data(gk + 24 * ncomps + c);
        const auto *gk_25 = buffer.data(gk + 25 * ncomps + c);
        const auto *gk_26 = buffer.data(gk + 26 * ncomps + c);
        const auto *gk_27 = buffer.data(gk + 27 * ncomps + c);
        const auto *gk_36 = buffer.data(gk + 36 * ncomps + c);
        const auto *gk_37 = buffer.data(gk + 37 * ncomps + c);
        const auto *gk_38 = buffer.data(gk + 38 * ncomps + c);
        const auto *gk_39 = buffer.data(gk + 39 * ncomps + c);
        const auto *gk_40 = buffer.data(gk + 40 * ncomps + c);
        const auto *gk_41 = buffer.data(gk + 41 * ncomps + c);
        const auto *gk_42 = buffer.data(gk + 42 * ncomps + c);
        const auto *gk_43 = buffer.data(gk + 43 * ncomps + c);
        const auto *gk_44 = buffer.data(gk + 44 * ncomps + c);
        const auto *gk_45 = buffer.data(gk + 45 * ncomps + c);
        const auto *gk_46 = buffer.data(gk + 46 * ncomps + c);
        const auto *gk_47 = buffer.data(gk + 47 * ncomps + c);
        const auto *gk_48 = buffer.data(gk + 48 * ncomps + c);
        const auto *gk_49 = buffer.data(gk + 49 * ncomps + c);
        const auto *gk_50 = buffer.data(gk + 50 * ncomps + c);
        const auto *gk_51 = buffer.data(gk + 51 * ncomps + c);
        const auto *gk_52 = buffer.data(gk + 52 * ncomps + c);
        const auto *gk_53 = buffer.data(gk + 53 * ncomps + c);
        const auto *gk_54 = buffer.data(gk + 54 * ncomps + c);
        const auto *gk_55 = buffer.data(gk + 55 * ncomps + c);
        const auto *gk_56 = buffer.data(gk + 56 * ncomps + c);
        const auto *gk_57 = buffer.data(gk + 57 * ncomps + c);
        const auto *gk_58 = buffer.data(gk + 58 * ncomps + c);
        const auto *gk_59 = buffer.data(gk + 59 * ncomps + c);
        const auto *gk_60 = buffer.data(gk + 60 * ncomps + c);
        const auto *gk_61 = buffer.data(gk + 61 * ncomps + c);
        const auto *gk_62 = buffer.data(gk + 62 * ncomps + c);
        const auto *gk_63 = buffer.data(gk + 63 * ncomps + c);
        const auto *gk_72 = buffer.data(gk + 72 * ncomps + c);
        const auto *gk_73 = buffer.data(gk + 73 * ncomps + c);
        const auto *gk_74 = buffer.data(gk + 74 * ncomps + c);
        const auto *gk_75 = buffer.data(gk + 75 * ncomps + c);
        const auto *gk_76 = buffer.data(gk + 76 * ncomps + c);
        const auto *gk_77 = buffer.data(gk + 77 * ncomps + c);
        const auto *gk_78 = buffer.data(gk + 78 * ncomps + c);
        const auto *gk_79 = buffer.data(gk + 79 * ncomps + c);
        const auto *gk_80 = buffer.data(gk + 80 * ncomps + c);
        const auto *gk_81 = buffer.data(gk + 81 * ncomps + c);
        const auto *gk_82 = buffer.data(gk + 82 * ncomps + c);
        const auto *gk_83 = buffer.data(gk + 83 * ncomps + c);
        const auto *gk_84 = buffer.data(gk + 84 * ncomps + c);
        const auto *gk_85 = buffer.data(gk + 85 * ncomps + c);
        const auto *gk_86 = buffer.data(gk + 86 * ncomps + c);
        const auto *gk_87 = buffer.data(gk + 87 * ncomps + c);
        const auto *gk_88 = buffer.data(gk + 88 * ncomps + c);
        const auto *gk_89 = buffer.data(gk + 89 * ncomps + c);
        const auto *gk_90 = buffer.data(gk + 90 * ncomps + c);
        const auto *gk_91 = buffer.data(gk + 91 * ncomps + c);
        const auto *gk_92 = buffer.data(gk + 92 * ncomps + c);
        const auto *gk_93 = buffer.data(gk + 93 * ncomps + c);
        const auto *gk_94 = buffer.data(gk + 94 * ncomps + c);
        const auto *gk_95 = buffer.data(gk + 95 * ncomps + c);
        const auto *gk_96 = buffer.data(gk + 96 * ncomps + c);
        const auto *gk_97 = buffer.data(gk + 97 * ncomps + c);
        const auto *gk_98 = buffer.data(gk + 98 * ncomps + c);
        const auto *gk_99 = buffer.data(gk + 99 * ncomps + c);
        const auto *gk_108 = buffer.data(gk + 108 * ncomps + c);
        const auto *gk_109 = buffer.data(gk + 109 * ncomps + c);
        const auto *gk_110 = buffer.data(gk + 110 * ncomps + c);
        const auto *gk_111 = buffer.data(gk + 111 * ncomps + c);
        const auto *gk_112 = buffer.data(gk + 112 * ncomps + c);
        const auto *gk_113 = buffer.data(gk + 113 * ncomps + c);
        const auto *gk_114 = buffer.data(gk + 114 * ncomps + c);
        const auto *gk_115 = buffer.data(gk + 115 * ncomps + c);
        const auto *gk_116 = buffer.data(gk + 116 * ncomps + c);
        const auto *gk_117 = buffer.data(gk + 117 * ncomps + c);
        const auto *gk_118 = buffer.data(gk + 118 * ncomps + c);
        const auto *gk_119 = buffer.data(gk + 119 * ncomps + c);
        const auto *gk_120 = buffer.data(gk + 120 * ncomps + c);
        const auto *gk_121 = buffer.data(gk + 121 * ncomps + c);
        const auto *gk_122 = buffer.data(gk + 122 * ncomps + c);
        const auto *gk_123 = buffer.data(gk + 123 * ncomps + c);
        const auto *gk_124 = buffer.data(gk + 124 * ncomps + c);
        const auto *gk_125 = buffer.data(gk + 125 * ncomps + c);
        const auto *gk_126 = buffer.data(gk + 126 * ncomps + c);
        const auto *gk_127 = buffer.data(gk + 127 * ncomps + c);
        const auto *gk_128 = buffer.data(gk + 128 * ncomps + c);
        const auto *gk_129 = buffer.data(gk + 129 * ncomps + c);
        const auto *gk_130 = buffer.data(gk + 130 * ncomps + c);
        const auto *gk_131 = buffer.data(gk + 131 * ncomps + c);
        const auto *gk_132 = buffer.data(gk + 132 * ncomps + c);
        const auto *gk_133 = buffer.data(gk + 133 * ncomps + c);
        const auto *gk_134 = buffer.data(gk + 134 * ncomps + c);
        const auto *gk_135 = buffer.data(gk + 135 * ncomps + c);
        const auto *gk_144 = buffer.data(gk + 144 * ncomps + c);
        const auto *gk_145 = buffer.data(gk + 145 * ncomps + c);
        const auto *gk_146 = buffer.data(gk + 146 * ncomps + c);
        const auto *gk_147 = buffer.data(gk + 147 * ncomps + c);
        const auto *gk_148 = buffer.data(gk + 148 * ncomps + c);
        const auto *gk_149 = buffer.data(gk + 149 * ncomps + c);
        const auto *gk_150 = buffer.data(gk + 150 * ncomps + c);
        const auto *gk_151 = buffer.data(gk + 151 * ncomps + c);
        const auto *gk_152 = buffer.data(gk + 152 * ncomps + c);
        const auto *gk_153 = buffer.data(gk + 153 * ncomps + c);
        const auto *gk_154 = buffer.data(gk + 154 * ncomps + c);
        const auto *gk_155 = buffer.data(gk + 155 * ncomps + c);
        const auto *gk_156 = buffer.data(gk + 156 * ncomps + c);
        const auto *gk_157 = buffer.data(gk + 157 * ncomps + c);
        const auto *gk_158 = buffer.data(gk + 158 * ncomps + c);
        const auto *gk_159 = buffer.data(gk + 159 * ncomps + c);
        const auto *gk_160 = buffer.data(gk + 160 * ncomps + c);
        const auto *gk_161 = buffer.data(gk + 161 * ncomps + c);
        const auto *gk_162 = buffer.data(gk + 162 * ncomps + c);
        const auto *gk_163 = buffer.data(gk + 163 * ncomps + c);
        const auto *gk_164 = buffer.data(gk + 164 * ncomps + c);
        const auto *gk_165 = buffer.data(gk + 165 * ncomps + c);
        const auto *gk_166 = buffer.data(gk + 166 * ncomps + c);
        const auto *gk_167 = buffer.data(gk + 167 * ncomps + c);
        const auto *gk_168 = buffer.data(gk + 168 * ncomps + c);
        const auto *gk_169 = buffer.data(gk + 169 * ncomps + c);
        const auto *gk_170 = buffer.data(gk + 170 * ncomps + c);
        const auto *gk_171 = buffer.data(gk + 171 * ncomps + c);
        const auto *gk_180 = buffer.data(gk + 180 * ncomps + c);
        const auto *gk_181 = buffer.data(gk + 181 * ncomps + c);
        const auto *gk_182 = buffer.data(gk + 182 * ncomps + c);
        const auto *gk_183 = buffer.data(gk + 183 * ncomps + c);
        const auto *gk_184 = buffer.data(gk + 184 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gi_0, gi_1, gi_2, gi_3, gi_4, gk_0, \
                         gk_1, gk_2, gk_3, gk_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * gi_0[k]
                     + gk_0[k];

            t_1[k] = -ab_x[k] * gi_1[k]
                     + gk_1[k];

            t_2[k] = -ab_x[k] * gi_2[k]
                     + gk_2[k];

            t_3[k] = -ab_x[k] * gi_3[k]
                     + gk_3[k];

            t_4[k] = -ab_x[k] * gi_4[k]
                     + gk_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gi_5, gi_6, gi_7, gi_8, gi_9, gk_5, \
                         gk_6, gk_7, gk_8, gk_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * gi_5[k]
                     + gk_5[k];

            t_6[k] = -ab_x[k] * gi_6[k]
                     + gk_6[k];

            t_7[k] = -ab_x[k] * gi_7[k]
                     + gk_7[k];

            t_8[k] = -ab_x[k] * gi_8[k]
                     + gk_8[k];

            t_9[k] = -ab_x[k] * gi_9[k]
                     + gk_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, gi_10, gi_11, gi_12, gi_13, \
                         gi_14, gk_10, gk_11, gk_12, gk_13, gk_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * gi_10[k]
                      + gk_10[k];

            t_11[k] = -ab_x[k] * gi_11[k]
                      + gk_11[k];

            t_12[k] = -ab_x[k] * gi_12[k]
                      + gk_12[k];

            t_13[k] = -ab_x[k] * gi_13[k]
                      + gk_13[k];

            t_14[k] = -ab_x[k] * gi_14[k]
                      + gk_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, gi_15, gi_16, gi_17, gi_18, \
                         gi_19, gk_15, gk_16, gk_17, gk_18, gk_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * gi_15[k]
                      + gk_15[k];

            t_16[k] = -ab_x[k] * gi_16[k]
                      + gk_16[k];

            t_17[k] = -ab_x[k] * gi_17[k]
                      + gk_17[k];

            t_18[k] = -ab_x[k] * gi_18[k]
                      + gk_18[k];

            t_19[k] = -ab_x[k] * gi_19[k]
                      + gk_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gi_20, gi_21, gi_22, gi_23, \
                         gi_24, gk_20, gk_21, gk_22, gk_23, gk_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * gi_20[k]
                      + gk_20[k];

            t_21[k] = -ab_x[k] * gi_21[k]
                      + gk_21[k];

            t_22[k] = -ab_x[k] * gi_22[k]
                      + gk_22[k];

            t_23[k] = -ab_x[k] * gi_23[k]
                      + gk_23[k];

            t_24[k] = -ab_x[k] * gi_24[k]
                      + gk_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, gi_25, gi_26, gi_27, gi_28, \
                         gi_29, gk_25, gk_26, gk_27, gk_36, gk_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * gi_25[k]
                      + gk_25[k];

            t_26[k] = -ab_x[k] * gi_26[k]
                      + gk_26[k];

            t_27[k] = -ab_x[k] * gi_27[k]
                      + gk_27[k];

            t_28[k] = -ab_x[k] * gi_28[k]
                      + gk_36[k];

            t_29[k] = -ab_x[k] * gi_29[k]
                      + gk_37[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gi_30, gi_31, gi_32, gi_33, \
                         gi_34, gk_38, gk_39, gk_40, gk_41, gk_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * gi_30[k]
                      + gk_38[k];

            t_31[k] = -ab_x[k] * gi_31[k]
                      + gk_39[k];

            t_32[k] = -ab_x[k] * gi_32[k]
                      + gk_40[k];

            t_33[k] = -ab_x[k] * gi_33[k]
                      + gk_41[k];

            t_34[k] = -ab_x[k] * gi_34[k]
                      + gk_42[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gi_35, gi_36, gi_37, gi_38, \
                         gi_39, gk_43, gk_44, gk_45, gk_46, gk_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * gi_35[k]
                      + gk_43[k];

            t_36[k] = -ab_x[k] * gi_36[k]
                      + gk_44[k];

            t_37[k] = -ab_x[k] * gi_37[k]
                      + gk_45[k];

            t_38[k] = -ab_x[k] * gi_38[k]
                      + gk_46[k];

            t_39[k] = -ab_x[k] * gi_39[k]
                      + gk_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gi_40, gi_41, gi_42, gi_43, \
                         gi_44, gk_48, gk_49, gk_50, gk_51, gk_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * gi_40[k]
                      + gk_48[k];

            t_41[k] = -ab_x[k] * gi_41[k]
                      + gk_49[k];

            t_42[k] = -ab_x[k] * gi_42[k]
                      + gk_50[k];

            t_43[k] = -ab_x[k] * gi_43[k]
                      + gk_51[k];

            t_44[k] = -ab_x[k] * gi_44[k]
                      + gk_52[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, gi_45, gi_46, gi_47, gi_48, \
                         gi_49, gk_53, gk_54, gk_55, gk_56, gk_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * gi_45[k]
                      + gk_53[k];

            t_46[k] = -ab_x[k] * gi_46[k]
                      + gk_54[k];

            t_47[k] = -ab_x[k] * gi_47[k]
                      + gk_55[k];

            t_48[k] = -ab_x[k] * gi_48[k]
                      + gk_56[k];

            t_49[k] = -ab_x[k] * gi_49[k]
                      + gk_57[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gi_50, gi_51, gi_52, gi_53, \
                         gi_54, gk_58, gk_59, gk_60, gk_61, gk_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * gi_50[k]
                      + gk_58[k];

            t_51[k] = -ab_x[k] * gi_51[k]
                      + gk_59[k];

            t_52[k] = -ab_x[k] * gi_52[k]
                      + gk_60[k];

            t_53[k] = -ab_x[k] * gi_53[k]
                      + gk_61[k];

            t_54[k] = -ab_x[k] * gi_54[k]
                      + gk_62[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, gi_55, gi_56, gi_57, gi_58, \
                         gi_59, gk_63, gk_72, gk_73, gk_74, gk_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * gi_55[k]
                      + gk_63[k];

            t_56[k] = -ab_x[k] * gi_56[k]
                      + gk_72[k];

            t_57[k] = -ab_x[k] * gi_57[k]
                      + gk_73[k];

            t_58[k] = -ab_x[k] * gi_58[k]
                      + gk_74[k];

            t_59[k] = -ab_x[k] * gi_59[k]
                      + gk_75[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gi_60, gi_61, gi_62, gi_63, \
                         gi_64, gk_76, gk_77, gk_78, gk_79, gk_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * gi_60[k]
                      + gk_76[k];

            t_61[k] = -ab_x[k] * gi_61[k]
                      + gk_77[k];

            t_62[k] = -ab_x[k] * gi_62[k]
                      + gk_78[k];

            t_63[k] = -ab_x[k] * gi_63[k]
                      + gk_79[k];

            t_64[k] = -ab_x[k] * gi_64[k]
                      + gk_80[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gi_65, gi_66, gi_67, gi_68, \
                         gi_69, gk_81, gk_82, gk_83, gk_84, gk_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * gi_65[k]
                      + gk_81[k];

            t_66[k] = -ab_x[k] * gi_66[k]
                      + gk_82[k];

            t_67[k] = -ab_x[k] * gi_67[k]
                      + gk_83[k];

            t_68[k] = -ab_x[k] * gi_68[k]
                      + gk_84[k];

            t_69[k] = -ab_x[k] * gi_69[k]
                      + gk_85[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, gi_70, gi_71, gi_72, gi_73, \
                         gi_74, gk_86, gk_87, gk_88, gk_89, gk_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * gi_70[k]
                      + gk_86[k];

            t_71[k] = -ab_x[k] * gi_71[k]
                      + gk_87[k];

            t_72[k] = -ab_x[k] * gi_72[k]
                      + gk_88[k];

            t_73[k] = -ab_x[k] * gi_73[k]
                      + gk_89[k];

            t_74[k] = -ab_x[k] * gi_74[k]
                      + gk_90[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, gi_75, gi_76, gi_77, gi_78, \
                         gi_79, gk_91, gk_92, gk_93, gk_94, gk_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * gi_75[k]
                      + gk_91[k];

            t_76[k] = -ab_x[k] * gi_76[k]
                      + gk_92[k];

            t_77[k] = -ab_x[k] * gi_77[k]
                      + gk_93[k];

            t_78[k] = -ab_x[k] * gi_78[k]
                      + gk_94[k];

            t_79[k] = -ab_x[k] * gi_79[k]
                      + gk_95[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gi_80, gi_81, gi_82, gi_83, \
                         gi_84, gk_96, gk_97, gk_98, gk_99, gk_108 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * gi_80[k]
                      + gk_96[k];

            t_81[k] = -ab_x[k] * gi_81[k]
                      + gk_97[k];

            t_82[k] = -ab_x[k] * gi_82[k]
                      + gk_98[k];

            t_83[k] = -ab_x[k] * gi_83[k]
                      + gk_99[k];

            t_84[k] = -ab_x[k] * gi_84[k]
                      + gk_108[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, gi_85, gi_86, gi_87, gi_88, \
                         gi_89, gk_109, gk_110, gk_111, gk_112, \
                         gk_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * gi_85[k]
                      + gk_109[k];

            t_86[k] = -ab_x[k] * gi_86[k]
                      + gk_110[k];

            t_87[k] = -ab_x[k] * gi_87[k]
                      + gk_111[k];

            t_88[k] = -ab_x[k] * gi_88[k]
                      + gk_112[k];

            t_89[k] = -ab_x[k] * gi_89[k]
                      + gk_113[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gi_90, gi_91, gi_92, gi_93, \
                         gi_94, gk_114, gk_115, gk_116, gk_117, \
                         gk_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * gi_90[k]
                      + gk_114[k];

            t_91[k] = -ab_x[k] * gi_91[k]
                      + gk_115[k];

            t_92[k] = -ab_x[k] * gi_92[k]
                      + gk_116[k];

            t_93[k] = -ab_x[k] * gi_93[k]
                      + gk_117[k];

            t_94[k] = -ab_x[k] * gi_94[k]
                      + gk_118[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gi_95, gi_96, gi_97, gi_98, \
                         gi_99, gk_119, gk_120, gk_121, gk_122, \
                         gk_123 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * gi_95[k]
                      + gk_119[k];

            t_96[k] = -ab_x[k] * gi_96[k]
                      + gk_120[k];

            t_97[k] = -ab_x[k] * gi_97[k]
                      + gk_121[k];

            t_98[k] = -ab_x[k] * gi_98[k]
                      + gk_122[k];

            t_99[k] = -ab_x[k] * gi_99[k]
                      + gk_123[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, gi_100, gi_101, gi_102, \
                         gi_103, gi_104, gk_124, gk_125, gk_126, gk_127, \
                         gk_128 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * gi_100[k]
                       + gk_124[k];

            t_101[k] = -ab_x[k] * gi_101[k]
                       + gk_125[k];

            t_102[k] = -ab_x[k] * gi_102[k]
                       + gk_126[k];

            t_103[k] = -ab_x[k] * gi_103[k]
                       + gk_127[k];

            t_104[k] = -ab_x[k] * gi_104[k]
                       + gk_128[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, gi_105, gi_106, gi_107, \
                         gi_108, gi_109, gk_129, gk_130, gk_131, gk_132, \
                         gk_133 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * gi_105[k]
                       + gk_129[k];

            t_106[k] = -ab_x[k] * gi_106[k]
                       + gk_130[k];

            t_107[k] = -ab_x[k] * gi_107[k]
                       + gk_131[k];

            t_108[k] = -ab_x[k] * gi_108[k]
                       + gk_132[k];

            t_109[k] = -ab_x[k] * gi_109[k]
                       + gk_133[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gi_110, gi_111, gi_112, \
                         gi_113, gi_114, gk_134, gk_135, gk_144, gk_145, \
                         gk_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * gi_110[k]
                       + gk_134[k];

            t_111[k] = -ab_x[k] * gi_111[k]
                       + gk_135[k];

            t_112[k] = -ab_x[k] * gi_112[k]
                       + gk_144[k];

            t_113[k] = -ab_x[k] * gi_113[k]
                       + gk_145[k];

            t_114[k] = -ab_x[k] * gi_114[k]
                       + gk_146[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, gi_115, gi_116, gi_117, \
                         gi_118, gi_119, gk_147, gk_148, gk_149, gk_150, \
                         gk_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * gi_115[k]
                       + gk_147[k];

            t_116[k] = -ab_x[k] * gi_116[k]
                       + gk_148[k];

            t_117[k] = -ab_x[k] * gi_117[k]
                       + gk_149[k];

            t_118[k] = -ab_x[k] * gi_118[k]
                       + gk_150[k];

            t_119[k] = -ab_x[k] * gi_119[k]
                       + gk_151[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gi_120, gi_121, gi_122, \
                         gi_123, gi_124, gk_152, gk_153, gk_154, gk_155, \
                         gk_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * gi_120[k]
                       + gk_152[k];

            t_121[k] = -ab_x[k] * gi_121[k]
                       + gk_153[k];

            t_122[k] = -ab_x[k] * gi_122[k]
                       + gk_154[k];

            t_123[k] = -ab_x[k] * gi_123[k]
                       + gk_155[k];

            t_124[k] = -ab_x[k] * gi_124[k]
                       + gk_156[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gi_125, gi_126, gi_127, \
                         gi_128, gi_129, gk_157, gk_158, gk_159, gk_160, \
                         gk_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * gi_125[k]
                       + gk_157[k];

            t_126[k] = -ab_x[k] * gi_126[k]
                       + gk_158[k];

            t_127[k] = -ab_x[k] * gi_127[k]
                       + gk_159[k];

            t_128[k] = -ab_x[k] * gi_128[k]
                       + gk_160[k];

            t_129[k] = -ab_x[k] * gi_129[k]
                       + gk_161[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, gi_130, gi_131, gi_132, \
                         gi_133, gi_134, gk_162, gk_163, gk_164, gk_165, \
                         gk_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * gi_130[k]
                       + gk_162[k];

            t_131[k] = -ab_x[k] * gi_131[k]
                       + gk_163[k];

            t_132[k] = -ab_x[k] * gi_132[k]
                       + gk_164[k];

            t_133[k] = -ab_x[k] * gi_133[k]
                       + gk_165[k];

            t_134[k] = -ab_x[k] * gi_134[k]
                       + gk_166[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, gi_135, gi_136, gi_137, \
                         gi_138, gi_139, gk_167, gk_168, gk_169, gk_170, \
                         gk_171 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * gi_135[k]
                       + gk_167[k];

            t_136[k] = -ab_x[k] * gi_136[k]
                       + gk_168[k];

            t_137[k] = -ab_x[k] * gi_137[k]
                       + gk_169[k];

            t_138[k] = -ab_x[k] * gi_138[k]
                       + gk_170[k];

            t_139[k] = -ab_x[k] * gi_139[k]
                       + gk_171[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gi_140, gi_141, gi_142, \
                         gi_143, gi_144, gk_180, gk_181, gk_182, gk_183, \
                         gk_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * gi_140[k]
                       + gk_180[k];

            t_141[k] = -ab_x[k] * gi_141[k]
                       + gk_181[k];

            t_142[k] = -ab_x[k] * gi_142[k]
                       + gk_182[k];

            t_143[k] = -ab_x[k] * gi_143[k]
                       + gk_183[k];

            t_144[k] = -ab_x[k] * gi_144[k]
                       + gk_184[k];
        }
    }
}

static auto
compute_hrr_hi_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gi, const size_t gk, const size_t ncomps,
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

        const auto *gi_145 = buffer.data(gi + 145 * ncomps + c);
        const auto *gi_146 = buffer.data(gi + 146 * ncomps + c);
        const auto *gi_147 = buffer.data(gi + 147 * ncomps + c);
        const auto *gi_148 = buffer.data(gi + 148 * ncomps + c);
        const auto *gi_149 = buffer.data(gi + 149 * ncomps + c);
        const auto *gi_150 = buffer.data(gi + 150 * ncomps + c);
        const auto *gi_151 = buffer.data(gi + 151 * ncomps + c);
        const auto *gi_152 = buffer.data(gi + 152 * ncomps + c);
        const auto *gi_153 = buffer.data(gi + 153 * ncomps + c);
        const auto *gi_154 = buffer.data(gi + 154 * ncomps + c);
        const auto *gi_155 = buffer.data(gi + 155 * ncomps + c);
        const auto *gi_156 = buffer.data(gi + 156 * ncomps + c);
        const auto *gi_157 = buffer.data(gi + 157 * ncomps + c);
        const auto *gi_158 = buffer.data(gi + 158 * ncomps + c);
        const auto *gi_159 = buffer.data(gi + 159 * ncomps + c);
        const auto *gi_160 = buffer.data(gi + 160 * ncomps + c);
        const auto *gi_161 = buffer.data(gi + 161 * ncomps + c);
        const auto *gi_162 = buffer.data(gi + 162 * ncomps + c);
        const auto *gi_163 = buffer.data(gi + 163 * ncomps + c);
        const auto *gi_164 = buffer.data(gi + 164 * ncomps + c);
        const auto *gi_165 = buffer.data(gi + 165 * ncomps + c);
        const auto *gi_166 = buffer.data(gi + 166 * ncomps + c);
        const auto *gi_167 = buffer.data(gi + 167 * ncomps + c);
        const auto *gi_168 = buffer.data(gi + 168 * ncomps + c);
        const auto *gi_169 = buffer.data(gi + 169 * ncomps + c);
        const auto *gi_170 = buffer.data(gi + 170 * ncomps + c);
        const auto *gi_171 = buffer.data(gi + 171 * ncomps + c);
        const auto *gi_172 = buffer.data(gi + 172 * ncomps + c);
        const auto *gi_173 = buffer.data(gi + 173 * ncomps + c);
        const auto *gi_174 = buffer.data(gi + 174 * ncomps + c);
        const auto *gi_175 = buffer.data(gi + 175 * ncomps + c);
        const auto *gi_176 = buffer.data(gi + 176 * ncomps + c);
        const auto *gi_177 = buffer.data(gi + 177 * ncomps + c);
        const auto *gi_178 = buffer.data(gi + 178 * ncomps + c);
        const auto *gi_179 = buffer.data(gi + 179 * ncomps + c);
        const auto *gi_180 = buffer.data(gi + 180 * ncomps + c);
        const auto *gi_181 = buffer.data(gi + 181 * ncomps + c);
        const auto *gi_182 = buffer.data(gi + 182 * ncomps + c);
        const auto *gi_183 = buffer.data(gi + 183 * ncomps + c);
        const auto *gi_184 = buffer.data(gi + 184 * ncomps + c);
        const auto *gi_185 = buffer.data(gi + 185 * ncomps + c);
        const auto *gi_186 = buffer.data(gi + 186 * ncomps + c);
        const auto *gi_187 = buffer.data(gi + 187 * ncomps + c);
        const auto *gi_188 = buffer.data(gi + 188 * ncomps + c);
        const auto *gi_189 = buffer.data(gi + 189 * ncomps + c);
        const auto *gi_190 = buffer.data(gi + 190 * ncomps + c);
        const auto *gi_191 = buffer.data(gi + 191 * ncomps + c);
        const auto *gi_192 = buffer.data(gi + 192 * ncomps + c);
        const auto *gi_193 = buffer.data(gi + 193 * ncomps + c);
        const auto *gi_194 = buffer.data(gi + 194 * ncomps + c);
        const auto *gi_195 = buffer.data(gi + 195 * ncomps + c);
        const auto *gi_196 = buffer.data(gi + 196 * ncomps + c);
        const auto *gi_197 = buffer.data(gi + 197 * ncomps + c);
        const auto *gi_198 = buffer.data(gi + 198 * ncomps + c);
        const auto *gi_199 = buffer.data(gi + 199 * ncomps + c);
        const auto *gi_200 = buffer.data(gi + 200 * ncomps + c);
        const auto *gi_201 = buffer.data(gi + 201 * ncomps + c);
        const auto *gi_202 = buffer.data(gi + 202 * ncomps + c);
        const auto *gi_203 = buffer.data(gi + 203 * ncomps + c);
        const auto *gi_204 = buffer.data(gi + 204 * ncomps + c);
        const auto *gi_205 = buffer.data(gi + 205 * ncomps + c);
        const auto *gi_206 = buffer.data(gi + 206 * ncomps + c);
        const auto *gi_207 = buffer.data(gi + 207 * ncomps + c);
        const auto *gi_208 = buffer.data(gi + 208 * ncomps + c);
        const auto *gi_209 = buffer.data(gi + 209 * ncomps + c);
        const auto *gi_210 = buffer.data(gi + 210 * ncomps + c);
        const auto *gi_211 = buffer.data(gi + 211 * ncomps + c);
        const auto *gi_212 = buffer.data(gi + 212 * ncomps + c);
        const auto *gi_213 = buffer.data(gi + 213 * ncomps + c);
        const auto *gi_214 = buffer.data(gi + 214 * ncomps + c);
        const auto *gi_215 = buffer.data(gi + 215 * ncomps + c);
        const auto *gi_216 = buffer.data(gi + 216 * ncomps + c);
        const auto *gi_217 = buffer.data(gi + 217 * ncomps + c);
        const auto *gi_218 = buffer.data(gi + 218 * ncomps + c);
        const auto *gi_219 = buffer.data(gi + 219 * ncomps + c);
        const auto *gi_220 = buffer.data(gi + 220 * ncomps + c);
        const auto *gi_221 = buffer.data(gi + 221 * ncomps + c);
        const auto *gi_222 = buffer.data(gi + 222 * ncomps + c);
        const auto *gi_223 = buffer.data(gi + 223 * ncomps + c);
        const auto *gi_224 = buffer.data(gi + 224 * ncomps + c);
        const auto *gi_225 = buffer.data(gi + 225 * ncomps + c);
        const auto *gi_226 = buffer.data(gi + 226 * ncomps + c);
        const auto *gi_227 = buffer.data(gi + 227 * ncomps + c);
        const auto *gi_228 = buffer.data(gi + 228 * ncomps + c);
        const auto *gi_229 = buffer.data(gi + 229 * ncomps + c);
        const auto *gi_230 = buffer.data(gi + 230 * ncomps + c);
        const auto *gi_231 = buffer.data(gi + 231 * ncomps + c);
        const auto *gi_232 = buffer.data(gi + 232 * ncomps + c);
        const auto *gi_233 = buffer.data(gi + 233 * ncomps + c);
        const auto *gi_234 = buffer.data(gi + 234 * ncomps + c);
        const auto *gi_235 = buffer.data(gi + 235 * ncomps + c);
        const auto *gi_236 = buffer.data(gi + 236 * ncomps + c);
        const auto *gi_237 = buffer.data(gi + 237 * ncomps + c);
        const auto *gi_238 = buffer.data(gi + 238 * ncomps + c);
        const auto *gi_239 = buffer.data(gi + 239 * ncomps + c);
        const auto *gi_240 = buffer.data(gi + 240 * ncomps + c);
        const auto *gi_241 = buffer.data(gi + 241 * ncomps + c);
        const auto *gi_242 = buffer.data(gi + 242 * ncomps + c);
        const auto *gi_243 = buffer.data(gi + 243 * ncomps + c);
        const auto *gi_244 = buffer.data(gi + 244 * ncomps + c);
        const auto *gi_245 = buffer.data(gi + 245 * ncomps + c);
        const auto *gi_246 = buffer.data(gi + 246 * ncomps + c);
        const auto *gi_247 = buffer.data(gi + 247 * ncomps + c);
        const auto *gi_248 = buffer.data(gi + 248 * ncomps + c);
        const auto *gi_249 = buffer.data(gi + 249 * ncomps + c);
        const auto *gi_250 = buffer.data(gi + 250 * ncomps + c);
        const auto *gi_251 = buffer.data(gi + 251 * ncomps + c);
        const auto *gi_252 = buffer.data(gi + 252 * ncomps + c);
        const auto *gi_253 = buffer.data(gi + 253 * ncomps + c);
        const auto *gi_254 = buffer.data(gi + 254 * ncomps + c);
        const auto *gi_255 = buffer.data(gi + 255 * ncomps + c);
        const auto *gi_256 = buffer.data(gi + 256 * ncomps + c);
        const auto *gi_257 = buffer.data(gi + 257 * ncomps + c);
        const auto *gi_258 = buffer.data(gi + 258 * ncomps + c);
        const auto *gi_259 = buffer.data(gi + 259 * ncomps + c);
        const auto *gi_260 = buffer.data(gi + 260 * ncomps + c);
        const auto *gi_261 = buffer.data(gi + 261 * ncomps + c);
        const auto *gi_262 = buffer.data(gi + 262 * ncomps + c);
        const auto *gi_263 = buffer.data(gi + 263 * ncomps + c);
        const auto *gi_264 = buffer.data(gi + 264 * ncomps + c);
        const auto *gi_265 = buffer.data(gi + 265 * ncomps + c);
        const auto *gi_266 = buffer.data(gi + 266 * ncomps + c);
        const auto *gi_267 = buffer.data(gi + 267 * ncomps + c);
        const auto *gi_268 = buffer.data(gi + 268 * ncomps + c);
        const auto *gi_269 = buffer.data(gi + 269 * ncomps + c);
        const auto *gi_270 = buffer.data(gi + 270 * ncomps + c);
        const auto *gi_271 = buffer.data(gi + 271 * ncomps + c);
        const auto *gi_272 = buffer.data(gi + 272 * ncomps + c);
        const auto *gi_273 = buffer.data(gi + 273 * ncomps + c);
        const auto *gi_274 = buffer.data(gi + 274 * ncomps + c);
        const auto *gi_275 = buffer.data(gi + 275 * ncomps + c);
        const auto *gi_276 = buffer.data(gi + 276 * ncomps + c);
        const auto *gi_277 = buffer.data(gi + 277 * ncomps + c);
        const auto *gi_278 = buffer.data(gi + 278 * ncomps + c);
        const auto *gi_279 = buffer.data(gi + 279 * ncomps + c);
        const auto *gi_280 = buffer.data(gi + 280 * ncomps + c);
        const auto *gi_281 = buffer.data(gi + 281 * ncomps + c);
        const auto *gi_282 = buffer.data(gi + 282 * ncomps + c);
        const auto *gi_283 = buffer.data(gi + 283 * ncomps + c);
        const auto *gi_284 = buffer.data(gi + 284 * ncomps + c);
        const auto *gi_285 = buffer.data(gi + 285 * ncomps + c);
        const auto *gi_286 = buffer.data(gi + 286 * ncomps + c);
        const auto *gi_287 = buffer.data(gi + 287 * ncomps + c);
        const auto *gi_288 = buffer.data(gi + 288 * ncomps + c);
        const auto *gi_289 = buffer.data(gi + 289 * ncomps + c);

        const auto *gk_185 = buffer.data(gk + 185 * ncomps + c);
        const auto *gk_186 = buffer.data(gk + 186 * ncomps + c);
        const auto *gk_187 = buffer.data(gk + 187 * ncomps + c);
        const auto *gk_188 = buffer.data(gk + 188 * ncomps + c);
        const auto *gk_189 = buffer.data(gk + 189 * ncomps + c);
        const auto *gk_190 = buffer.data(gk + 190 * ncomps + c);
        const auto *gk_191 = buffer.data(gk + 191 * ncomps + c);
        const auto *gk_192 = buffer.data(gk + 192 * ncomps + c);
        const auto *gk_193 = buffer.data(gk + 193 * ncomps + c);
        const auto *gk_194 = buffer.data(gk + 194 * ncomps + c);
        const auto *gk_195 = buffer.data(gk + 195 * ncomps + c);
        const auto *gk_196 = buffer.data(gk + 196 * ncomps + c);
        const auto *gk_197 = buffer.data(gk + 197 * ncomps + c);
        const auto *gk_198 = buffer.data(gk + 198 * ncomps + c);
        const auto *gk_199 = buffer.data(gk + 199 * ncomps + c);
        const auto *gk_200 = buffer.data(gk + 200 * ncomps + c);
        const auto *gk_201 = buffer.data(gk + 201 * ncomps + c);
        const auto *gk_202 = buffer.data(gk + 202 * ncomps + c);
        const auto *gk_203 = buffer.data(gk + 203 * ncomps + c);
        const auto *gk_204 = buffer.data(gk + 204 * ncomps + c);
        const auto *gk_205 = buffer.data(gk + 205 * ncomps + c);
        const auto *gk_206 = buffer.data(gk + 206 * ncomps + c);
        const auto *gk_207 = buffer.data(gk + 207 * ncomps + c);
        const auto *gk_216 = buffer.data(gk + 216 * ncomps + c);
        const auto *gk_217 = buffer.data(gk + 217 * ncomps + c);
        const auto *gk_218 = buffer.data(gk + 218 * ncomps + c);
        const auto *gk_219 = buffer.data(gk + 219 * ncomps + c);
        const auto *gk_220 = buffer.data(gk + 220 * ncomps + c);
        const auto *gk_221 = buffer.data(gk + 221 * ncomps + c);
        const auto *gk_222 = buffer.data(gk + 222 * ncomps + c);
        const auto *gk_223 = buffer.data(gk + 223 * ncomps + c);
        const auto *gk_224 = buffer.data(gk + 224 * ncomps + c);
        const auto *gk_225 = buffer.data(gk + 225 * ncomps + c);
        const auto *gk_226 = buffer.data(gk + 226 * ncomps + c);
        const auto *gk_227 = buffer.data(gk + 227 * ncomps + c);
        const auto *gk_228 = buffer.data(gk + 228 * ncomps + c);
        const auto *gk_229 = buffer.data(gk + 229 * ncomps + c);
        const auto *gk_230 = buffer.data(gk + 230 * ncomps + c);
        const auto *gk_231 = buffer.data(gk + 231 * ncomps + c);
        const auto *gk_232 = buffer.data(gk + 232 * ncomps + c);
        const auto *gk_233 = buffer.data(gk + 233 * ncomps + c);
        const auto *gk_234 = buffer.data(gk + 234 * ncomps + c);
        const auto *gk_235 = buffer.data(gk + 235 * ncomps + c);
        const auto *gk_236 = buffer.data(gk + 236 * ncomps + c);
        const auto *gk_237 = buffer.data(gk + 237 * ncomps + c);
        const auto *gk_238 = buffer.data(gk + 238 * ncomps + c);
        const auto *gk_239 = buffer.data(gk + 239 * ncomps + c);
        const auto *gk_240 = buffer.data(gk + 240 * ncomps + c);
        const auto *gk_241 = buffer.data(gk + 241 * ncomps + c);
        const auto *gk_242 = buffer.data(gk + 242 * ncomps + c);
        const auto *gk_243 = buffer.data(gk + 243 * ncomps + c);
        const auto *gk_252 = buffer.data(gk + 252 * ncomps + c);
        const auto *gk_253 = buffer.data(gk + 253 * ncomps + c);
        const auto *gk_254 = buffer.data(gk + 254 * ncomps + c);
        const auto *gk_255 = buffer.data(gk + 255 * ncomps + c);
        const auto *gk_256 = buffer.data(gk + 256 * ncomps + c);
        const auto *gk_257 = buffer.data(gk + 257 * ncomps + c);
        const auto *gk_258 = buffer.data(gk + 258 * ncomps + c);
        const auto *gk_259 = buffer.data(gk + 259 * ncomps + c);
        const auto *gk_260 = buffer.data(gk + 260 * ncomps + c);
        const auto *gk_261 = buffer.data(gk + 261 * ncomps + c);
        const auto *gk_262 = buffer.data(gk + 262 * ncomps + c);
        const auto *gk_263 = buffer.data(gk + 263 * ncomps + c);
        const auto *gk_264 = buffer.data(gk + 264 * ncomps + c);
        const auto *gk_265 = buffer.data(gk + 265 * ncomps + c);
        const auto *gk_266 = buffer.data(gk + 266 * ncomps + c);
        const auto *gk_267 = buffer.data(gk + 267 * ncomps + c);
        const auto *gk_268 = buffer.data(gk + 268 * ncomps + c);
        const auto *gk_269 = buffer.data(gk + 269 * ncomps + c);
        const auto *gk_270 = buffer.data(gk + 270 * ncomps + c);
        const auto *gk_271 = buffer.data(gk + 271 * ncomps + c);
        const auto *gk_272 = buffer.data(gk + 272 * ncomps + c);
        const auto *gk_273 = buffer.data(gk + 273 * ncomps + c);
        const auto *gk_274 = buffer.data(gk + 274 * ncomps + c);
        const auto *gk_275 = buffer.data(gk + 275 * ncomps + c);
        const auto *gk_276 = buffer.data(gk + 276 * ncomps + c);
        const auto *gk_277 = buffer.data(gk + 277 * ncomps + c);
        const auto *gk_278 = buffer.data(gk + 278 * ncomps + c);
        const auto *gk_279 = buffer.data(gk + 279 * ncomps + c);
        const auto *gk_288 = buffer.data(gk + 288 * ncomps + c);
        const auto *gk_289 = buffer.data(gk + 289 * ncomps + c);
        const auto *gk_290 = buffer.data(gk + 290 * ncomps + c);
        const auto *gk_291 = buffer.data(gk + 291 * ncomps + c);
        const auto *gk_292 = buffer.data(gk + 292 * ncomps + c);
        const auto *gk_293 = buffer.data(gk + 293 * ncomps + c);
        const auto *gk_294 = buffer.data(gk + 294 * ncomps + c);
        const auto *gk_295 = buffer.data(gk + 295 * ncomps + c);
        const auto *gk_296 = buffer.data(gk + 296 * ncomps + c);
        const auto *gk_297 = buffer.data(gk + 297 * ncomps + c);
        const auto *gk_298 = buffer.data(gk + 298 * ncomps + c);
        const auto *gk_299 = buffer.data(gk + 299 * ncomps + c);
        const auto *gk_300 = buffer.data(gk + 300 * ncomps + c);
        const auto *gk_301 = buffer.data(gk + 301 * ncomps + c);
        const auto *gk_302 = buffer.data(gk + 302 * ncomps + c);
        const auto *gk_303 = buffer.data(gk + 303 * ncomps + c);
        const auto *gk_304 = buffer.data(gk + 304 * ncomps + c);
        const auto *gk_305 = buffer.data(gk + 305 * ncomps + c);
        const auto *gk_306 = buffer.data(gk + 306 * ncomps + c);
        const auto *gk_307 = buffer.data(gk + 307 * ncomps + c);
        const auto *gk_308 = buffer.data(gk + 308 * ncomps + c);
        const auto *gk_309 = buffer.data(gk + 309 * ncomps + c);
        const auto *gk_310 = buffer.data(gk + 310 * ncomps + c);
        const auto *gk_311 = buffer.data(gk + 311 * ncomps + c);
        const auto *gk_312 = buffer.data(gk + 312 * ncomps + c);
        const auto *gk_313 = buffer.data(gk + 313 * ncomps + c);
        const auto *gk_314 = buffer.data(gk + 314 * ncomps + c);
        const auto *gk_315 = buffer.data(gk + 315 * ncomps + c);
        const auto *gk_324 = buffer.data(gk + 324 * ncomps + c);
        const auto *gk_325 = buffer.data(gk + 325 * ncomps + c);
        const auto *gk_326 = buffer.data(gk + 326 * ncomps + c);
        const auto *gk_327 = buffer.data(gk + 327 * ncomps + c);
        const auto *gk_328 = buffer.data(gk + 328 * ncomps + c);
        const auto *gk_329 = buffer.data(gk + 329 * ncomps + c);
        const auto *gk_330 = buffer.data(gk + 330 * ncomps + c);
        const auto *gk_331 = buffer.data(gk + 331 * ncomps + c);
        const auto *gk_332 = buffer.data(gk + 332 * ncomps + c);
        const auto *gk_333 = buffer.data(gk + 333 * ncomps + c);
        const auto *gk_334 = buffer.data(gk + 334 * ncomps + c);
        const auto *gk_335 = buffer.data(gk + 335 * ncomps + c);
        const auto *gk_336 = buffer.data(gk + 336 * ncomps + c);
        const auto *gk_337 = buffer.data(gk + 337 * ncomps + c);
        const auto *gk_338 = buffer.data(gk + 338 * ncomps + c);
        const auto *gk_339 = buffer.data(gk + 339 * ncomps + c);
        const auto *gk_340 = buffer.data(gk + 340 * ncomps + c);
        const auto *gk_341 = buffer.data(gk + 341 * ncomps + c);
        const auto *gk_342 = buffer.data(gk + 342 * ncomps + c);
        const auto *gk_343 = buffer.data(gk + 343 * ncomps + c);
        const auto *gk_344 = buffer.data(gk + 344 * ncomps + c);
        const auto *gk_345 = buffer.data(gk + 345 * ncomps + c);
        const auto *gk_346 = buffer.data(gk + 346 * ncomps + c);
        const auto *gk_347 = buffer.data(gk + 347 * ncomps + c);
        const auto *gk_348 = buffer.data(gk + 348 * ncomps + c);
        const auto *gk_349 = buffer.data(gk + 349 * ncomps + c);
        const auto *gk_350 = buffer.data(gk + 350 * ncomps + c);
        const auto *gk_351 = buffer.data(gk + 351 * ncomps + c);
        const auto *gk_360 = buffer.data(gk + 360 * ncomps + c);
        const auto *gk_361 = buffer.data(gk + 361 * ncomps + c);
        const auto *gk_362 = buffer.data(gk + 362 * ncomps + c);
        const auto *gk_363 = buffer.data(gk + 363 * ncomps + c);
        const auto *gk_364 = buffer.data(gk + 364 * ncomps + c);
        const auto *gk_365 = buffer.data(gk + 365 * ncomps + c);
        const auto *gk_366 = buffer.data(gk + 366 * ncomps + c);
        const auto *gk_367 = buffer.data(gk + 367 * ncomps + c);
        const auto *gk_368 = buffer.data(gk + 368 * ncomps + c);
        const auto *gk_369 = buffer.data(gk + 369 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, gi_145, gi_146, gi_147, \
                         gi_148, gi_149, gk_185, gk_186, gk_187, gk_188, \
                         gk_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * gi_145[k]
                       + gk_185[k];

            t_146[k] = -ab_x[k] * gi_146[k]
                       + gk_186[k];

            t_147[k] = -ab_x[k] * gi_147[k]
                       + gk_187[k];

            t_148[k] = -ab_x[k] * gi_148[k]
                       + gk_188[k];

            t_149[k] = -ab_x[k] * gi_149[k]
                       + gk_189[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gi_150, gi_151, gi_152, \
                         gi_153, gi_154, gk_190, gk_191, gk_192, gk_193, \
                         gk_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * gi_150[k]
                       + gk_190[k];

            t_151[k] = -ab_x[k] * gi_151[k]
                       + gk_191[k];

            t_152[k] = -ab_x[k] * gi_152[k]
                       + gk_192[k];

            t_153[k] = -ab_x[k] * gi_153[k]
                       + gk_193[k];

            t_154[k] = -ab_x[k] * gi_154[k]
                       + gk_194[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gi_155, gi_156, gi_157, \
                         gi_158, gi_159, gk_195, gk_196, gk_197, gk_198, \
                         gk_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * gi_155[k]
                       + gk_195[k];

            t_156[k] = -ab_x[k] * gi_156[k]
                       + gk_196[k];

            t_157[k] = -ab_x[k] * gi_157[k]
                       + gk_197[k];

            t_158[k] = -ab_x[k] * gi_158[k]
                       + gk_198[k];

            t_159[k] = -ab_x[k] * gi_159[k]
                       + gk_199[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, gi_160, gi_161, gi_162, \
                         gi_163, gi_164, gk_200, gk_201, gk_202, gk_203, \
                         gk_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * gi_160[k]
                       + gk_200[k];

            t_161[k] = -ab_x[k] * gi_161[k]
                       + gk_201[k];

            t_162[k] = -ab_x[k] * gi_162[k]
                       + gk_202[k];

            t_163[k] = -ab_x[k] * gi_163[k]
                       + gk_203[k];

            t_164[k] = -ab_x[k] * gi_164[k]
                       + gk_204[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, gi_165, gi_166, gi_167, \
                         gi_168, gi_169, gk_205, gk_206, gk_207, gk_216, \
                         gk_217 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * gi_165[k]
                       + gk_205[k];

            t_166[k] = -ab_x[k] * gi_166[k]
                       + gk_206[k];

            t_167[k] = -ab_x[k] * gi_167[k]
                       + gk_207[k];

            t_168[k] = -ab_x[k] * gi_168[k]
                       + gk_216[k];

            t_169[k] = -ab_x[k] * gi_169[k]
                       + gk_217[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, gi_170, gi_171, gi_172, \
                         gi_173, gi_174, gk_218, gk_219, gk_220, gk_221, \
                         gk_222 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * gi_170[k]
                       + gk_218[k];

            t_171[k] = -ab_x[k] * gi_171[k]
                       + gk_219[k];

            t_172[k] = -ab_x[k] * gi_172[k]
                       + gk_220[k];

            t_173[k] = -ab_x[k] * gi_173[k]
                       + gk_221[k];

            t_174[k] = -ab_x[k] * gi_174[k]
                       + gk_222[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, gi_175, gi_176, gi_177, \
                         gi_178, gi_179, gk_223, gk_224, gk_225, gk_226, \
                         gk_227 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * gi_175[k]
                       + gk_223[k];

            t_176[k] = -ab_x[k] * gi_176[k]
                       + gk_224[k];

            t_177[k] = -ab_x[k] * gi_177[k]
                       + gk_225[k];

            t_178[k] = -ab_x[k] * gi_178[k]
                       + gk_226[k];

            t_179[k] = -ab_x[k] * gi_179[k]
                       + gk_227[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gi_180, gi_181, gi_182, \
                         gi_183, gi_184, gk_228, gk_229, gk_230, gk_231, \
                         gk_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * gi_180[k]
                       + gk_228[k];

            t_181[k] = -ab_x[k] * gi_181[k]
                       + gk_229[k];

            t_182[k] = -ab_x[k] * gi_182[k]
                       + gk_230[k];

            t_183[k] = -ab_x[k] * gi_183[k]
                       + gk_231[k];

            t_184[k] = -ab_x[k] * gi_184[k]
                       + gk_232[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gi_185, gi_186, gi_187, \
                         gi_188, gi_189, gk_233, gk_234, gk_235, gk_236, \
                         gk_237 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * gi_185[k]
                       + gk_233[k];

            t_186[k] = -ab_x[k] * gi_186[k]
                       + gk_234[k];

            t_187[k] = -ab_x[k] * gi_187[k]
                       + gk_235[k];

            t_188[k] = -ab_x[k] * gi_188[k]
                       + gk_236[k];

            t_189[k] = -ab_x[k] * gi_189[k]
                       + gk_237[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, gi_190, gi_191, gi_192, \
                         gi_193, gi_194, gk_238, gk_239, gk_240, gk_241, \
                         gk_242 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * gi_190[k]
                       + gk_238[k];

            t_191[k] = -ab_x[k] * gi_191[k]
                       + gk_239[k];

            t_192[k] = -ab_x[k] * gi_192[k]
                       + gk_240[k];

            t_193[k] = -ab_x[k] * gi_193[k]
                       + gk_241[k];

            t_194[k] = -ab_x[k] * gi_194[k]
                       + gk_242[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, gi_195, gi_196, gi_197, \
                         gi_198, gi_199, gk_243, gk_252, gk_253, gk_254, \
                         gk_255 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * gi_195[k]
                       + gk_243[k];

            t_196[k] = -ab_x[k] * gi_196[k]
                       + gk_252[k];

            t_197[k] = -ab_x[k] * gi_197[k]
                       + gk_253[k];

            t_198[k] = -ab_x[k] * gi_198[k]
                       + gk_254[k];

            t_199[k] = -ab_x[k] * gi_199[k]
                       + gk_255[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, gi_200, gi_201, gi_202, \
                         gi_203, gi_204, gk_256, gk_257, gk_258, gk_259, \
                         gk_260 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * gi_200[k]
                       + gk_256[k];

            t_201[k] = -ab_x[k] * gi_201[k]
                       + gk_257[k];

            t_202[k] = -ab_x[k] * gi_202[k]
                       + gk_258[k];

            t_203[k] = -ab_x[k] * gi_203[k]
                       + gk_259[k];

            t_204[k] = -ab_x[k] * gi_204[k]
                       + gk_260[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, gi_205, gi_206, gi_207, \
                         gi_208, gi_209, gk_261, gk_262, gk_263, gk_264, \
                         gk_265 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * gi_205[k]
                       + gk_261[k];

            t_206[k] = -ab_x[k] * gi_206[k]
                       + gk_262[k];

            t_207[k] = -ab_x[k] * gi_207[k]
                       + gk_263[k];

            t_208[k] = -ab_x[k] * gi_208[k]
                       + gk_264[k];

            t_209[k] = -ab_x[k] * gi_209[k]
                       + gk_265[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gi_210, gi_211, gi_212, \
                         gi_213, gi_214, gk_266, gk_267, gk_268, gk_269, \
                         gk_270 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * gi_210[k]
                       + gk_266[k];

            t_211[k] = -ab_x[k] * gi_211[k]
                       + gk_267[k];

            t_212[k] = -ab_x[k] * gi_212[k]
                       + gk_268[k];

            t_213[k] = -ab_x[k] * gi_213[k]
                       + gk_269[k];

            t_214[k] = -ab_x[k] * gi_214[k]
                       + gk_270[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gi_215, gi_216, gi_217, \
                         gi_218, gi_219, gk_271, gk_272, gk_273, gk_274, \
                         gk_275 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * gi_215[k]
                       + gk_271[k];

            t_216[k] = -ab_x[k] * gi_216[k]
                       + gk_272[k];

            t_217[k] = -ab_x[k] * gi_217[k]
                       + gk_273[k];

            t_218[k] = -ab_x[k] * gi_218[k]
                       + gk_274[k];

            t_219[k] = -ab_x[k] * gi_219[k]
                       + gk_275[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, gi_220, gi_221, gi_222, \
                         gi_223, gi_224, gk_276, gk_277, gk_278, gk_279, \
                         gk_288 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_x[k] * gi_220[k]
                       + gk_276[k];

            t_221[k] = -ab_x[k] * gi_221[k]
                       + gk_277[k];

            t_222[k] = -ab_x[k] * gi_222[k]
                       + gk_278[k];

            t_223[k] = -ab_x[k] * gi_223[k]
                       + gk_279[k];

            t_224[k] = -ab_x[k] * gi_224[k]
                       + gk_288[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, gi_225, gi_226, gi_227, \
                         gi_228, gi_229, gk_289, gk_290, gk_291, gk_292, \
                         gk_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_x[k] * gi_225[k]
                       + gk_289[k];

            t_226[k] = -ab_x[k] * gi_226[k]
                       + gk_290[k];

            t_227[k] = -ab_x[k] * gi_227[k]
                       + gk_291[k];

            t_228[k] = -ab_x[k] * gi_228[k]
                       + gk_292[k];

            t_229[k] = -ab_x[k] * gi_229[k]
                       + gk_293[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, gi_230, gi_231, gi_232, \
                         gi_233, gi_234, gk_294, gk_295, gk_296, gk_297, \
                         gk_298 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_x[k] * gi_230[k]
                       + gk_294[k];

            t_231[k] = -ab_x[k] * gi_231[k]
                       + gk_295[k];

            t_232[k] = -ab_x[k] * gi_232[k]
                       + gk_296[k];

            t_233[k] = -ab_x[k] * gi_233[k]
                       + gk_297[k];

            t_234[k] = -ab_x[k] * gi_234[k]
                       + gk_298[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, gi_235, gi_236, gi_237, \
                         gi_238, gi_239, gk_299, gk_300, gk_301, gk_302, \
                         gk_303 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_x[k] * gi_235[k]
                       + gk_299[k];

            t_236[k] = -ab_x[k] * gi_236[k]
                       + gk_300[k];

            t_237[k] = -ab_x[k] * gi_237[k]
                       + gk_301[k];

            t_238[k] = -ab_x[k] * gi_238[k]
                       + gk_302[k];

            t_239[k] = -ab_x[k] * gi_239[k]
                       + gk_303[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, gi_240, gi_241, gi_242, \
                         gi_243, gi_244, gk_304, gk_305, gk_306, gk_307, \
                         gk_308 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_x[k] * gi_240[k]
                       + gk_304[k];

            t_241[k] = -ab_x[k] * gi_241[k]
                       + gk_305[k];

            t_242[k] = -ab_x[k] * gi_242[k]
                       + gk_306[k];

            t_243[k] = -ab_x[k] * gi_243[k]
                       + gk_307[k];

            t_244[k] = -ab_x[k] * gi_244[k]
                       + gk_308[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, gi_245, gi_246, gi_247, \
                         gi_248, gi_249, gk_309, gk_310, gk_311, gk_312, \
                         gk_313 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_x[k] * gi_245[k]
                       + gk_309[k];

            t_246[k] = -ab_x[k] * gi_246[k]
                       + gk_310[k];

            t_247[k] = -ab_x[k] * gi_247[k]
                       + gk_311[k];

            t_248[k] = -ab_x[k] * gi_248[k]
                       + gk_312[k];

            t_249[k] = -ab_x[k] * gi_249[k]
                       + gk_313[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, gi_250, gi_251, gi_252, \
                         gi_253, gi_254, gk_314, gk_315, gk_324, gk_325, \
                         gk_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_x[k] * gi_250[k]
                       + gk_314[k];

            t_251[k] = -ab_x[k] * gi_251[k]
                       + gk_315[k];

            t_252[k] = -ab_x[k] * gi_252[k]
                       + gk_324[k];

            t_253[k] = -ab_x[k] * gi_253[k]
                       + gk_325[k];

            t_254[k] = -ab_x[k] * gi_254[k]
                       + gk_326[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, gi_255, gi_256, gi_257, \
                         gi_258, gi_259, gk_327, gk_328, gk_329, gk_330, \
                         gk_331 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_x[k] * gi_255[k]
                       + gk_327[k];

            t_256[k] = -ab_x[k] * gi_256[k]
                       + gk_328[k];

            t_257[k] = -ab_x[k] * gi_257[k]
                       + gk_329[k];

            t_258[k] = -ab_x[k] * gi_258[k]
                       + gk_330[k];

            t_259[k] = -ab_x[k] * gi_259[k]
                       + gk_331[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, gi_260, gi_261, gi_262, \
                         gi_263, gi_264, gk_332, gk_333, gk_334, gk_335, \
                         gk_336 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_x[k] * gi_260[k]
                       + gk_332[k];

            t_261[k] = -ab_x[k] * gi_261[k]
                       + gk_333[k];

            t_262[k] = -ab_x[k] * gi_262[k]
                       + gk_334[k];

            t_263[k] = -ab_x[k] * gi_263[k]
                       + gk_335[k];

            t_264[k] = -ab_x[k] * gi_264[k]
                       + gk_336[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, gi_265, gi_266, gi_267, \
                         gi_268, gi_269, gk_337, gk_338, gk_339, gk_340, \
                         gk_341 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_x[k] * gi_265[k]
                       + gk_337[k];

            t_266[k] = -ab_x[k] * gi_266[k]
                       + gk_338[k];

            t_267[k] = -ab_x[k] * gi_267[k]
                       + gk_339[k];

            t_268[k] = -ab_x[k] * gi_268[k]
                       + gk_340[k];

            t_269[k] = -ab_x[k] * gi_269[k]
                       + gk_341[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, gi_270, gi_271, gi_272, \
                         gi_273, gi_274, gk_342, gk_343, gk_344, gk_345, \
                         gk_346 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_x[k] * gi_270[k]
                       + gk_342[k];

            t_271[k] = -ab_x[k] * gi_271[k]
                       + gk_343[k];

            t_272[k] = -ab_x[k] * gi_272[k]
                       + gk_344[k];

            t_273[k] = -ab_x[k] * gi_273[k]
                       + gk_345[k];

            t_274[k] = -ab_x[k] * gi_274[k]
                       + gk_346[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, gi_275, gi_276, gi_277, \
                         gi_278, gi_279, gk_347, gk_348, gk_349, gk_350, \
                         gk_351 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_x[k] * gi_275[k]
                       + gk_347[k];

            t_276[k] = -ab_x[k] * gi_276[k]
                       + gk_348[k];

            t_277[k] = -ab_x[k] * gi_277[k]
                       + gk_349[k];

            t_278[k] = -ab_x[k] * gi_278[k]
                       + gk_350[k];

            t_279[k] = -ab_x[k] * gi_279[k]
                       + gk_351[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, gi_280, gi_281, gi_282, \
                         gi_283, gi_284, gk_360, gk_361, gk_362, gk_363, \
                         gk_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_x[k] * gi_280[k]
                       + gk_360[k];

            t_281[k] = -ab_x[k] * gi_281[k]
                       + gk_361[k];

            t_282[k] = -ab_x[k] * gi_282[k]
                       + gk_362[k];

            t_283[k] = -ab_x[k] * gi_283[k]
                       + gk_363[k];

            t_284[k] = -ab_x[k] * gi_284[k]
                       + gk_364[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, gi_285, gi_286, gi_287, \
                         gi_288, gi_289, gk_365, gk_366, gk_367, gk_368, \
                         gk_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_x[k] * gi_285[k]
                       + gk_365[k];

            t_286[k] = -ab_x[k] * gi_286[k]
                       + gk_366[k];

            t_287[k] = -ab_x[k] * gi_287[k]
                       + gk_367[k];

            t_288[k] = -ab_x[k] * gi_288[k]
                       + gk_368[k];

            t_289[k] = -ab_x[k] * gi_289[k]
                       + gk_369[k];
        }
    }
}

static auto
compute_hrr_hi_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gi, const size_t gk, const size_t ncomps,
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
        const auto *ab_y = coordinates.data(7);

        const auto *gi_280 = buffer.data(gi + 280 * ncomps + c);
        const auto *gi_281 = buffer.data(gi + 281 * ncomps + c);
        const auto *gi_282 = buffer.data(gi + 282 * ncomps + c);
        const auto *gi_283 = buffer.data(gi + 283 * ncomps + c);
        const auto *gi_284 = buffer.data(gi + 284 * ncomps + c);
        const auto *gi_285 = buffer.data(gi + 285 * ncomps + c);
        const auto *gi_286 = buffer.data(gi + 286 * ncomps + c);
        const auto *gi_287 = buffer.data(gi + 287 * ncomps + c);
        const auto *gi_288 = buffer.data(gi + 288 * ncomps + c);
        const auto *gi_289 = buffer.data(gi + 289 * ncomps + c);
        const auto *gi_290 = buffer.data(gi + 290 * ncomps + c);
        const auto *gi_291 = buffer.data(gi + 291 * ncomps + c);
        const auto *gi_292 = buffer.data(gi + 292 * ncomps + c);
        const auto *gi_293 = buffer.data(gi + 293 * ncomps + c);
        const auto *gi_294 = buffer.data(gi + 294 * ncomps + c);
        const auto *gi_295 = buffer.data(gi + 295 * ncomps + c);
        const auto *gi_296 = buffer.data(gi + 296 * ncomps + c);
        const auto *gi_297 = buffer.data(gi + 297 * ncomps + c);
        const auto *gi_298 = buffer.data(gi + 298 * ncomps + c);
        const auto *gi_299 = buffer.data(gi + 299 * ncomps + c);
        const auto *gi_300 = buffer.data(gi + 300 * ncomps + c);
        const auto *gi_301 = buffer.data(gi + 301 * ncomps + c);
        const auto *gi_302 = buffer.data(gi + 302 * ncomps + c);
        const auto *gi_303 = buffer.data(gi + 303 * ncomps + c);
        const auto *gi_304 = buffer.data(gi + 304 * ncomps + c);
        const auto *gi_305 = buffer.data(gi + 305 * ncomps + c);
        const auto *gi_306 = buffer.data(gi + 306 * ncomps + c);
        const auto *gi_307 = buffer.data(gi + 307 * ncomps + c);
        const auto *gi_308 = buffer.data(gi + 308 * ncomps + c);
        const auto *gi_309 = buffer.data(gi + 309 * ncomps + c);
        const auto *gi_310 = buffer.data(gi + 310 * ncomps + c);
        const auto *gi_311 = buffer.data(gi + 311 * ncomps + c);
        const auto *gi_312 = buffer.data(gi + 312 * ncomps + c);
        const auto *gi_313 = buffer.data(gi + 313 * ncomps + c);
        const auto *gi_314 = buffer.data(gi + 314 * ncomps + c);
        const auto *gi_315 = buffer.data(gi + 315 * ncomps + c);
        const auto *gi_316 = buffer.data(gi + 316 * ncomps + c);
        const auto *gi_317 = buffer.data(gi + 317 * ncomps + c);
        const auto *gi_318 = buffer.data(gi + 318 * ncomps + c);
        const auto *gi_319 = buffer.data(gi + 319 * ncomps + c);
        const auto *gi_320 = buffer.data(gi + 320 * ncomps + c);
        const auto *gi_321 = buffer.data(gi + 321 * ncomps + c);
        const auto *gi_322 = buffer.data(gi + 322 * ncomps + c);
        const auto *gi_323 = buffer.data(gi + 323 * ncomps + c);
        const auto *gi_324 = buffer.data(gi + 324 * ncomps + c);
        const auto *gi_325 = buffer.data(gi + 325 * ncomps + c);
        const auto *gi_326 = buffer.data(gi + 326 * ncomps + c);
        const auto *gi_327 = buffer.data(gi + 327 * ncomps + c);
        const auto *gi_328 = buffer.data(gi + 328 * ncomps + c);
        const auto *gi_329 = buffer.data(gi + 329 * ncomps + c);
        const auto *gi_330 = buffer.data(gi + 330 * ncomps + c);
        const auto *gi_331 = buffer.data(gi + 331 * ncomps + c);
        const auto *gi_332 = buffer.data(gi + 332 * ncomps + c);
        const auto *gi_333 = buffer.data(gi + 333 * ncomps + c);
        const auto *gi_334 = buffer.data(gi + 334 * ncomps + c);
        const auto *gi_335 = buffer.data(gi + 335 * ncomps + c);
        const auto *gi_336 = buffer.data(gi + 336 * ncomps + c);
        const auto *gi_337 = buffer.data(gi + 337 * ncomps + c);
        const auto *gi_338 = buffer.data(gi + 338 * ncomps + c);
        const auto *gi_339 = buffer.data(gi + 339 * ncomps + c);
        const auto *gi_340 = buffer.data(gi + 340 * ncomps + c);
        const auto *gi_341 = buffer.data(gi + 341 * ncomps + c);
        const auto *gi_342 = buffer.data(gi + 342 * ncomps + c);
        const auto *gi_343 = buffer.data(gi + 343 * ncomps + c);
        const auto *gi_344 = buffer.data(gi + 344 * ncomps + c);
        const auto *gi_345 = buffer.data(gi + 345 * ncomps + c);
        const auto *gi_346 = buffer.data(gi + 346 * ncomps + c);
        const auto *gi_347 = buffer.data(gi + 347 * ncomps + c);
        const auto *gi_348 = buffer.data(gi + 348 * ncomps + c);
        const auto *gi_349 = buffer.data(gi + 349 * ncomps + c);
        const auto *gi_350 = buffer.data(gi + 350 * ncomps + c);
        const auto *gi_351 = buffer.data(gi + 351 * ncomps + c);
        const auto *gi_352 = buffer.data(gi + 352 * ncomps + c);
        const auto *gi_353 = buffer.data(gi + 353 * ncomps + c);
        const auto *gi_354 = buffer.data(gi + 354 * ncomps + c);
        const auto *gi_355 = buffer.data(gi + 355 * ncomps + c);
        const auto *gi_356 = buffer.data(gi + 356 * ncomps + c);
        const auto *gi_357 = buffer.data(gi + 357 * ncomps + c);
        const auto *gi_358 = buffer.data(gi + 358 * ncomps + c);
        const auto *gi_359 = buffer.data(gi + 359 * ncomps + c);
        const auto *gi_360 = buffer.data(gi + 360 * ncomps + c);
        const auto *gi_361 = buffer.data(gi + 361 * ncomps + c);
        const auto *gi_362 = buffer.data(gi + 362 * ncomps + c);
        const auto *gi_363 = buffer.data(gi + 363 * ncomps + c);
        const auto *gi_364 = buffer.data(gi + 364 * ncomps + c);
        const auto *gi_365 = buffer.data(gi + 365 * ncomps + c);
        const auto *gi_366 = buffer.data(gi + 366 * ncomps + c);
        const auto *gi_367 = buffer.data(gi + 367 * ncomps + c);
        const auto *gi_368 = buffer.data(gi + 368 * ncomps + c);
        const auto *gi_369 = buffer.data(gi + 369 * ncomps + c);
        const auto *gi_370 = buffer.data(gi + 370 * ncomps + c);
        const auto *gi_371 = buffer.data(gi + 371 * ncomps + c);
        const auto *gi_372 = buffer.data(gi + 372 * ncomps + c);
        const auto *gi_373 = buffer.data(gi + 373 * ncomps + c);
        const auto *gi_374 = buffer.data(gi + 374 * ncomps + c);
        const auto *gi_375 = buffer.data(gi + 375 * ncomps + c);
        const auto *gi_376 = buffer.data(gi + 376 * ncomps + c);
        const auto *gi_377 = buffer.data(gi + 377 * ncomps + c);
        const auto *gi_378 = buffer.data(gi + 378 * ncomps + c);
        const auto *gi_379 = buffer.data(gi + 379 * ncomps + c);
        const auto *gi_380 = buffer.data(gi + 380 * ncomps + c);
        const auto *gi_381 = buffer.data(gi + 381 * ncomps + c);
        const auto *gi_382 = buffer.data(gi + 382 * ncomps + c);
        const auto *gi_383 = buffer.data(gi + 383 * ncomps + c);
        const auto *gi_384 = buffer.data(gi + 384 * ncomps + c);
        const auto *gi_385 = buffer.data(gi + 385 * ncomps + c);
        const auto *gi_386 = buffer.data(gi + 386 * ncomps + c);
        const auto *gi_387 = buffer.data(gi + 387 * ncomps + c);
        const auto *gi_388 = buffer.data(gi + 388 * ncomps + c);
        const auto *gi_389 = buffer.data(gi + 389 * ncomps + c);
        const auto *gi_390 = buffer.data(gi + 390 * ncomps + c);
        const auto *gi_391 = buffer.data(gi + 391 * ncomps + c);
        const auto *gi_392 = buffer.data(gi + 392 * ncomps + c);
        const auto *gi_393 = buffer.data(gi + 393 * ncomps + c);
        const auto *gi_394 = buffer.data(gi + 394 * ncomps + c);
        const auto *gi_395 = buffer.data(gi + 395 * ncomps + c);
        const auto *gi_396 = buffer.data(gi + 396 * ncomps + c);
        const auto *gi_397 = buffer.data(gi + 397 * ncomps + c);
        const auto *gi_398 = buffer.data(gi + 398 * ncomps + c);
        const auto *gi_399 = buffer.data(gi + 399 * ncomps + c);
        const auto *gi_400 = buffer.data(gi + 400 * ncomps + c);
        const auto *gi_401 = buffer.data(gi + 401 * ncomps + c);
        const auto *gi_402 = buffer.data(gi + 402 * ncomps + c);
        const auto *gi_403 = buffer.data(gi + 403 * ncomps + c);
        const auto *gi_404 = buffer.data(gi + 404 * ncomps + c);
        const auto *gi_405 = buffer.data(gi + 405 * ncomps + c);
        const auto *gi_406 = buffer.data(gi + 406 * ncomps + c);
        const auto *gi_407 = buffer.data(gi + 407 * ncomps + c);
        const auto *gi_408 = buffer.data(gi + 408 * ncomps + c);
        const auto *gi_409 = buffer.data(gi + 409 * ncomps + c);
        const auto *gi_410 = buffer.data(gi + 410 * ncomps + c);
        const auto *gi_411 = buffer.data(gi + 411 * ncomps + c);
        const auto *gi_412 = buffer.data(gi + 412 * ncomps + c);
        const auto *gi_413 = buffer.data(gi + 413 * ncomps + c);
        const auto *gi_414 = buffer.data(gi + 414 * ncomps + c);
        const auto *gi_415 = buffer.data(gi + 415 * ncomps + c);
        const auto *gi_416 = buffer.data(gi + 416 * ncomps + c);
        const auto *gi_417 = buffer.data(gi + 417 * ncomps + c);
        const auto *gi_418 = buffer.data(gi + 418 * ncomps + c);
        const auto *gi_419 = buffer.data(gi + 419 * ncomps + c);

        const auto *gk_361 = buffer.data(gk + 361 * ncomps + c);
        const auto *gk_363 = buffer.data(gk + 363 * ncomps + c);
        const auto *gk_364 = buffer.data(gk + 364 * ncomps + c);
        const auto *gk_366 = buffer.data(gk + 366 * ncomps + c);
        const auto *gk_367 = buffer.data(gk + 367 * ncomps + c);
        const auto *gk_368 = buffer.data(gk + 368 * ncomps + c);
        const auto *gk_370 = buffer.data(gk + 370 * ncomps + c);
        const auto *gk_371 = buffer.data(gk + 371 * ncomps + c);
        const auto *gk_372 = buffer.data(gk + 372 * ncomps + c);
        const auto *gk_373 = buffer.data(gk + 373 * ncomps + c);
        const auto *gk_374 = buffer.data(gk + 374 * ncomps + c);
        const auto *gk_375 = buffer.data(gk + 375 * ncomps + c);
        const auto *gk_376 = buffer.data(gk + 376 * ncomps + c);
        const auto *gk_377 = buffer.data(gk + 377 * ncomps + c);
        const auto *gk_378 = buffer.data(gk + 378 * ncomps + c);
        const auto *gk_379 = buffer.data(gk + 379 * ncomps + c);
        const auto *gk_380 = buffer.data(gk + 380 * ncomps + c);
        const auto *gk_381 = buffer.data(gk + 381 * ncomps + c);
        const auto *gk_382 = buffer.data(gk + 382 * ncomps + c);
        const auto *gk_383 = buffer.data(gk + 383 * ncomps + c);
        const auto *gk_384 = buffer.data(gk + 384 * ncomps + c);
        const auto *gk_385 = buffer.data(gk + 385 * ncomps + c);
        const auto *gk_386 = buffer.data(gk + 386 * ncomps + c);
        const auto *gk_387 = buffer.data(gk + 387 * ncomps + c);
        const auto *gk_396 = buffer.data(gk + 396 * ncomps + c);
        const auto *gk_397 = buffer.data(gk + 397 * ncomps + c);
        const auto *gk_398 = buffer.data(gk + 398 * ncomps + c);
        const auto *gk_399 = buffer.data(gk + 399 * ncomps + c);
        const auto *gk_400 = buffer.data(gk + 400 * ncomps + c);
        const auto *gk_401 = buffer.data(gk + 401 * ncomps + c);
        const auto *gk_402 = buffer.data(gk + 402 * ncomps + c);
        const auto *gk_403 = buffer.data(gk + 403 * ncomps + c);
        const auto *gk_404 = buffer.data(gk + 404 * ncomps + c);
        const auto *gk_405 = buffer.data(gk + 405 * ncomps + c);
        const auto *gk_406 = buffer.data(gk + 406 * ncomps + c);
        const auto *gk_407 = buffer.data(gk + 407 * ncomps + c);
        const auto *gk_408 = buffer.data(gk + 408 * ncomps + c);
        const auto *gk_409 = buffer.data(gk + 409 * ncomps + c);
        const auto *gk_410 = buffer.data(gk + 410 * ncomps + c);
        const auto *gk_411 = buffer.data(gk + 411 * ncomps + c);
        const auto *gk_412 = buffer.data(gk + 412 * ncomps + c);
        const auto *gk_413 = buffer.data(gk + 413 * ncomps + c);
        const auto *gk_414 = buffer.data(gk + 414 * ncomps + c);
        const auto *gk_415 = buffer.data(gk + 415 * ncomps + c);
        const auto *gk_416 = buffer.data(gk + 416 * ncomps + c);
        const auto *gk_417 = buffer.data(gk + 417 * ncomps + c);
        const auto *gk_418 = buffer.data(gk + 418 * ncomps + c);
        const auto *gk_419 = buffer.data(gk + 419 * ncomps + c);
        const auto *gk_420 = buffer.data(gk + 420 * ncomps + c);
        const auto *gk_421 = buffer.data(gk + 421 * ncomps + c);
        const auto *gk_422 = buffer.data(gk + 422 * ncomps + c);
        const auto *gk_423 = buffer.data(gk + 423 * ncomps + c);
        const auto *gk_432 = buffer.data(gk + 432 * ncomps + c);
        const auto *gk_433 = buffer.data(gk + 433 * ncomps + c);
        const auto *gk_434 = buffer.data(gk + 434 * ncomps + c);
        const auto *gk_435 = buffer.data(gk + 435 * ncomps + c);
        const auto *gk_436 = buffer.data(gk + 436 * ncomps + c);
        const auto *gk_437 = buffer.data(gk + 437 * ncomps + c);
        const auto *gk_438 = buffer.data(gk + 438 * ncomps + c);
        const auto *gk_439 = buffer.data(gk + 439 * ncomps + c);
        const auto *gk_440 = buffer.data(gk + 440 * ncomps + c);
        const auto *gk_441 = buffer.data(gk + 441 * ncomps + c);
        const auto *gk_442 = buffer.data(gk + 442 * ncomps + c);
        const auto *gk_443 = buffer.data(gk + 443 * ncomps + c);
        const auto *gk_444 = buffer.data(gk + 444 * ncomps + c);
        const auto *gk_445 = buffer.data(gk + 445 * ncomps + c);
        const auto *gk_446 = buffer.data(gk + 446 * ncomps + c);
        const auto *gk_447 = buffer.data(gk + 447 * ncomps + c);
        const auto *gk_448 = buffer.data(gk + 448 * ncomps + c);
        const auto *gk_449 = buffer.data(gk + 449 * ncomps + c);
        const auto *gk_450 = buffer.data(gk + 450 * ncomps + c);
        const auto *gk_451 = buffer.data(gk + 451 * ncomps + c);
        const auto *gk_452 = buffer.data(gk + 452 * ncomps + c);
        const auto *gk_453 = buffer.data(gk + 453 * ncomps + c);
        const auto *gk_454 = buffer.data(gk + 454 * ncomps + c);
        const auto *gk_455 = buffer.data(gk + 455 * ncomps + c);
        const auto *gk_456 = buffer.data(gk + 456 * ncomps + c);
        const auto *gk_457 = buffer.data(gk + 457 * ncomps + c);
        const auto *gk_458 = buffer.data(gk + 458 * ncomps + c);
        const auto *gk_459 = buffer.data(gk + 459 * ncomps + c);
        const auto *gk_468 = buffer.data(gk + 468 * ncomps + c);
        const auto *gk_469 = buffer.data(gk + 469 * ncomps + c);
        const auto *gk_470 = buffer.data(gk + 470 * ncomps + c);
        const auto *gk_471 = buffer.data(gk + 471 * ncomps + c);
        const auto *gk_472 = buffer.data(gk + 472 * ncomps + c);
        const auto *gk_473 = buffer.data(gk + 473 * ncomps + c);
        const auto *gk_474 = buffer.data(gk + 474 * ncomps + c);
        const auto *gk_475 = buffer.data(gk + 475 * ncomps + c);
        const auto *gk_476 = buffer.data(gk + 476 * ncomps + c);
        const auto *gk_477 = buffer.data(gk + 477 * ncomps + c);
        const auto *gk_478 = buffer.data(gk + 478 * ncomps + c);
        const auto *gk_479 = buffer.data(gk + 479 * ncomps + c);
        const auto *gk_480 = buffer.data(gk + 480 * ncomps + c);
        const auto *gk_481 = buffer.data(gk + 481 * ncomps + c);
        const auto *gk_482 = buffer.data(gk + 482 * ncomps + c);
        const auto *gk_483 = buffer.data(gk + 483 * ncomps + c);
        const auto *gk_484 = buffer.data(gk + 484 * ncomps + c);
        const auto *gk_485 = buffer.data(gk + 485 * ncomps + c);
        const auto *gk_486 = buffer.data(gk + 486 * ncomps + c);
        const auto *gk_487 = buffer.data(gk + 487 * ncomps + c);
        const auto *gk_488 = buffer.data(gk + 488 * ncomps + c);
        const auto *gk_489 = buffer.data(gk + 489 * ncomps + c);
        const auto *gk_490 = buffer.data(gk + 490 * ncomps + c);
        const auto *gk_491 = buffer.data(gk + 491 * ncomps + c);
        const auto *gk_492 = buffer.data(gk + 492 * ncomps + c);
        const auto *gk_493 = buffer.data(gk + 493 * ncomps + c);
        const auto *gk_494 = buffer.data(gk + 494 * ncomps + c);
        const auto *gk_495 = buffer.data(gk + 495 * ncomps + c);
        const auto *gk_504 = buffer.data(gk + 504 * ncomps + c);
        const auto *gk_505 = buffer.data(gk + 505 * ncomps + c);
        const auto *gk_506 = buffer.data(gk + 506 * ncomps + c);
        const auto *gk_507 = buffer.data(gk + 507 * ncomps + c);
        const auto *gk_508 = buffer.data(gk + 508 * ncomps + c);
        const auto *gk_509 = buffer.data(gk + 509 * ncomps + c);
        const auto *gk_510 = buffer.data(gk + 510 * ncomps + c);
        const auto *gk_511 = buffer.data(gk + 511 * ncomps + c);
        const auto *gk_512 = buffer.data(gk + 512 * ncomps + c);
        const auto *gk_513 = buffer.data(gk + 513 * ncomps + c);
        const auto *gk_514 = buffer.data(gk + 514 * ncomps + c);
        const auto *gk_515 = buffer.data(gk + 515 * ncomps + c);
        const auto *gk_516 = buffer.data(gk + 516 * ncomps + c);
        const auto *gk_517 = buffer.data(gk + 517 * ncomps + c);
        const auto *gk_518 = buffer.data(gk + 518 * ncomps + c);
        const auto *gk_519 = buffer.data(gk + 519 * ncomps + c);
        const auto *gk_520 = buffer.data(gk + 520 * ncomps + c);
        const auto *gk_521 = buffer.data(gk + 521 * ncomps + c);
        const auto *gk_522 = buffer.data(gk + 522 * ncomps + c);
        const auto *gk_523 = buffer.data(gk + 523 * ncomps + c);
        const auto *gk_524 = buffer.data(gk + 524 * ncomps + c);
        const auto *gk_525 = buffer.data(gk + 525 * ncomps + c);
        const auto *gk_526 = buffer.data(gk + 526 * ncomps + c);
        const auto *gk_527 = buffer.data(gk + 527 * ncomps + c);
        const auto *gk_528 = buffer.data(gk + 528 * ncomps + c);
        const auto *gk_529 = buffer.data(gk + 529 * ncomps + c);
        const auto *gk_530 = buffer.data(gk + 530 * ncomps + c);
        const auto *gk_531 = buffer.data(gk + 531 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, gi_290, gi_291, gi_292, \
                         gi_293, gi_294, gk_370, gk_371, gk_372, gk_373, \
                         gk_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_x[k] * gi_290[k]
                       + gk_370[k];

            t_291[k] = -ab_x[k] * gi_291[k]
                       + gk_371[k];

            t_292[k] = -ab_x[k] * gi_292[k]
                       + gk_372[k];

            t_293[k] = -ab_x[k] * gi_293[k]
                       + gk_373[k];

            t_294[k] = -ab_x[k] * gi_294[k]
                       + gk_374[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, gi_295, gi_296, gi_297, \
                         gi_298, gi_299, gk_375, gk_376, gk_377, gk_378, \
                         gk_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_x[k] * gi_295[k]
                       + gk_375[k];

            t_296[k] = -ab_x[k] * gi_296[k]
                       + gk_376[k];

            t_297[k] = -ab_x[k] * gi_297[k]
                       + gk_377[k];

            t_298[k] = -ab_x[k] * gi_298[k]
                       + gk_378[k];

            t_299[k] = -ab_x[k] * gi_299[k]
                       + gk_379[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, gi_300, gi_301, gi_302, \
                         gi_303, gi_304, gk_380, gk_381, gk_382, gk_383, \
                         gk_384 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_x[k] * gi_300[k]
                       + gk_380[k];

            t_301[k] = -ab_x[k] * gi_301[k]
                       + gk_381[k];

            t_302[k] = -ab_x[k] * gi_302[k]
                       + gk_382[k];

            t_303[k] = -ab_x[k] * gi_303[k]
                       + gk_383[k];

            t_304[k] = -ab_x[k] * gi_304[k]
                       + gk_384[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, gi_305, gi_306, gi_307, \
                         gi_308, gi_309, gk_385, gk_386, gk_387, gk_396, \
                         gk_397 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_x[k] * gi_305[k]
                       + gk_385[k];

            t_306[k] = -ab_x[k] * gi_306[k]
                       + gk_386[k];

            t_307[k] = -ab_x[k] * gi_307[k]
                       + gk_387[k];

            t_308[k] = -ab_x[k] * gi_308[k]
                       + gk_396[k];

            t_309[k] = -ab_x[k] * gi_309[k]
                       + gk_397[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, gi_310, gi_311, gi_312, \
                         gi_313, gi_314, gk_398, gk_399, gk_400, gk_401, \
                         gk_402 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_x[k] * gi_310[k]
                       + gk_398[k];

            t_311[k] = -ab_x[k] * gi_311[k]
                       + gk_399[k];

            t_312[k] = -ab_x[k] * gi_312[k]
                       + gk_400[k];

            t_313[k] = -ab_x[k] * gi_313[k]
                       + gk_401[k];

            t_314[k] = -ab_x[k] * gi_314[k]
                       + gk_402[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, gi_315, gi_316, gi_317, \
                         gi_318, gi_319, gk_403, gk_404, gk_405, gk_406, \
                         gk_407 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_x[k] * gi_315[k]
                       + gk_403[k];

            t_316[k] = -ab_x[k] * gi_316[k]
                       + gk_404[k];

            t_317[k] = -ab_x[k] * gi_317[k]
                       + gk_405[k];

            t_318[k] = -ab_x[k] * gi_318[k]
                       + gk_406[k];

            t_319[k] = -ab_x[k] * gi_319[k]
                       + gk_407[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, gi_320, gi_321, gi_322, \
                         gi_323, gi_324, gk_408, gk_409, gk_410, gk_411, \
                         gk_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_x[k] * gi_320[k]
                       + gk_408[k];

            t_321[k] = -ab_x[k] * gi_321[k]
                       + gk_409[k];

            t_322[k] = -ab_x[k] * gi_322[k]
                       + gk_410[k];

            t_323[k] = -ab_x[k] * gi_323[k]
                       + gk_411[k];

            t_324[k] = -ab_x[k] * gi_324[k]
                       + gk_412[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, gi_325, gi_326, gi_327, \
                         gi_328, gi_329, gk_413, gk_414, gk_415, gk_416, \
                         gk_417 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_x[k] * gi_325[k]
                       + gk_413[k];

            t_326[k] = -ab_x[k] * gi_326[k]
                       + gk_414[k];

            t_327[k] = -ab_x[k] * gi_327[k]
                       + gk_415[k];

            t_328[k] = -ab_x[k] * gi_328[k]
                       + gk_416[k];

            t_329[k] = -ab_x[k] * gi_329[k]
                       + gk_417[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, gi_330, gi_331, gi_332, \
                         gi_333, gi_334, gk_418, gk_419, gk_420, gk_421, \
                         gk_422 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = -ab_x[k] * gi_330[k]
                       + gk_418[k];

            t_331[k] = -ab_x[k] * gi_331[k]
                       + gk_419[k];

            t_332[k] = -ab_x[k] * gi_332[k]
                       + gk_420[k];

            t_333[k] = -ab_x[k] * gi_333[k]
                       + gk_421[k];

            t_334[k] = -ab_x[k] * gi_334[k]
                       + gk_422[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, gi_335, gi_336, gi_337, \
                         gi_338, gi_339, gk_423, gk_432, gk_433, gk_434, \
                         gk_435 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = -ab_x[k] * gi_335[k]
                       + gk_423[k];

            t_336[k] = -ab_x[k] * gi_336[k]
                       + gk_432[k];

            t_337[k] = -ab_x[k] * gi_337[k]
                       + gk_433[k];

            t_338[k] = -ab_x[k] * gi_338[k]
                       + gk_434[k];

            t_339[k] = -ab_x[k] * gi_339[k]
                       + gk_435[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, gi_340, gi_341, gi_342, \
                         gi_343, gi_344, gk_436, gk_437, gk_438, gk_439, \
                         gk_440 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = -ab_x[k] * gi_340[k]
                       + gk_436[k];

            t_341[k] = -ab_x[k] * gi_341[k]
                       + gk_437[k];

            t_342[k] = -ab_x[k] * gi_342[k]
                       + gk_438[k];

            t_343[k] = -ab_x[k] * gi_343[k]
                       + gk_439[k];

            t_344[k] = -ab_x[k] * gi_344[k]
                       + gk_440[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, gi_345, gi_346, gi_347, \
                         gi_348, gi_349, gk_441, gk_442, gk_443, gk_444, \
                         gk_445 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = -ab_x[k] * gi_345[k]
                       + gk_441[k];

            t_346[k] = -ab_x[k] * gi_346[k]
                       + gk_442[k];

            t_347[k] = -ab_x[k] * gi_347[k]
                       + gk_443[k];

            t_348[k] = -ab_x[k] * gi_348[k]
                       + gk_444[k];

            t_349[k] = -ab_x[k] * gi_349[k]
                       + gk_445[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, gi_350, gi_351, gi_352, \
                         gi_353, gi_354, gk_446, gk_447, gk_448, gk_449, \
                         gk_450 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = -ab_x[k] * gi_350[k]
                       + gk_446[k];

            t_351[k] = -ab_x[k] * gi_351[k]
                       + gk_447[k];

            t_352[k] = -ab_x[k] * gi_352[k]
                       + gk_448[k];

            t_353[k] = -ab_x[k] * gi_353[k]
                       + gk_449[k];

            t_354[k] = -ab_x[k] * gi_354[k]
                       + gk_450[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, gi_355, gi_356, gi_357, \
                         gi_358, gi_359, gk_451, gk_452, gk_453, gk_454, \
                         gk_455 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = -ab_x[k] * gi_355[k]
                       + gk_451[k];

            t_356[k] = -ab_x[k] * gi_356[k]
                       + gk_452[k];

            t_357[k] = -ab_x[k] * gi_357[k]
                       + gk_453[k];

            t_358[k] = -ab_x[k] * gi_358[k]
                       + gk_454[k];

            t_359[k] = -ab_x[k] * gi_359[k]
                       + gk_455[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, gi_360, gi_361, gi_362, \
                         gi_363, gi_364, gk_456, gk_457, gk_458, gk_459, \
                         gk_468 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = -ab_x[k] * gi_360[k]
                       + gk_456[k];

            t_361[k] = -ab_x[k] * gi_361[k]
                       + gk_457[k];

            t_362[k] = -ab_x[k] * gi_362[k]
                       + gk_458[k];

            t_363[k] = -ab_x[k] * gi_363[k]
                       + gk_459[k];

            t_364[k] = -ab_x[k] * gi_364[k]
                       + gk_468[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, gi_365, gi_366, gi_367, \
                         gi_368, gi_369, gk_469, gk_470, gk_471, gk_472, \
                         gk_473 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = -ab_x[k] * gi_365[k]
                       + gk_469[k];

            t_366[k] = -ab_x[k] * gi_366[k]
                       + gk_470[k];

            t_367[k] = -ab_x[k] * gi_367[k]
                       + gk_471[k];

            t_368[k] = -ab_x[k] * gi_368[k]
                       + gk_472[k];

            t_369[k] = -ab_x[k] * gi_369[k]
                       + gk_473[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, gi_370, gi_371, gi_372, \
                         gi_373, gi_374, gk_474, gk_475, gk_476, gk_477, \
                         gk_478 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = -ab_x[k] * gi_370[k]
                       + gk_474[k];

            t_371[k] = -ab_x[k] * gi_371[k]
                       + gk_475[k];

            t_372[k] = -ab_x[k] * gi_372[k]
                       + gk_476[k];

            t_373[k] = -ab_x[k] * gi_373[k]
                       + gk_477[k];

            t_374[k] = -ab_x[k] * gi_374[k]
                       + gk_478[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, gi_375, gi_376, gi_377, \
                         gi_378, gi_379, gk_479, gk_480, gk_481, gk_482, \
                         gk_483 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = -ab_x[k] * gi_375[k]
                       + gk_479[k];

            t_376[k] = -ab_x[k] * gi_376[k]
                       + gk_480[k];

            t_377[k] = -ab_x[k] * gi_377[k]
                       + gk_481[k];

            t_378[k] = -ab_x[k] * gi_378[k]
                       + gk_482[k];

            t_379[k] = -ab_x[k] * gi_379[k]
                       + gk_483[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, gi_380, gi_381, gi_382, \
                         gi_383, gi_384, gk_484, gk_485, gk_486, gk_487, \
                         gk_488 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = -ab_x[k] * gi_380[k]
                       + gk_484[k];

            t_381[k] = -ab_x[k] * gi_381[k]
                       + gk_485[k];

            t_382[k] = -ab_x[k] * gi_382[k]
                       + gk_486[k];

            t_383[k] = -ab_x[k] * gi_383[k]
                       + gk_487[k];

            t_384[k] = -ab_x[k] * gi_384[k]
                       + gk_488[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, gi_385, gi_386, gi_387, \
                         gi_388, gi_389, gk_489, gk_490, gk_491, gk_492, \
                         gk_493 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = -ab_x[k] * gi_385[k]
                       + gk_489[k];

            t_386[k] = -ab_x[k] * gi_386[k]
                       + gk_490[k];

            t_387[k] = -ab_x[k] * gi_387[k]
                       + gk_491[k];

            t_388[k] = -ab_x[k] * gi_388[k]
                       + gk_492[k];

            t_389[k] = -ab_x[k] * gi_389[k]
                       + gk_493[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, gi_390, gi_391, gi_392, \
                         gi_393, gi_394, gk_494, gk_495, gk_504, gk_505, \
                         gk_506 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = -ab_x[k] * gi_390[k]
                       + gk_494[k];

            t_391[k] = -ab_x[k] * gi_391[k]
                       + gk_495[k];

            t_392[k] = -ab_x[k] * gi_392[k]
                       + gk_504[k];

            t_393[k] = -ab_x[k] * gi_393[k]
                       + gk_505[k];

            t_394[k] = -ab_x[k] * gi_394[k]
                       + gk_506[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, gi_395, gi_396, gi_397, \
                         gi_398, gi_399, gk_507, gk_508, gk_509, gk_510, \
                         gk_511 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = -ab_x[k] * gi_395[k]
                       + gk_507[k];

            t_396[k] = -ab_x[k] * gi_396[k]
                       + gk_508[k];

            t_397[k] = -ab_x[k] * gi_397[k]
                       + gk_509[k];

            t_398[k] = -ab_x[k] * gi_398[k]
                       + gk_510[k];

            t_399[k] = -ab_x[k] * gi_399[k]
                       + gk_511[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, gi_400, gi_401, gi_402, \
                         gi_403, gi_404, gk_512, gk_513, gk_514, gk_515, \
                         gk_516 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = -ab_x[k] * gi_400[k]
                       + gk_512[k];

            t_401[k] = -ab_x[k] * gi_401[k]
                       + gk_513[k];

            t_402[k] = -ab_x[k] * gi_402[k]
                       + gk_514[k];

            t_403[k] = -ab_x[k] * gi_403[k]
                       + gk_515[k];

            t_404[k] = -ab_x[k] * gi_404[k]
                       + gk_516[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, gi_405, gi_406, gi_407, \
                         gi_408, gi_409, gk_517, gk_518, gk_519, gk_520, \
                         gk_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = -ab_x[k] * gi_405[k]
                       + gk_517[k];

            t_406[k] = -ab_x[k] * gi_406[k]
                       + gk_518[k];

            t_407[k] = -ab_x[k] * gi_407[k]
                       + gk_519[k];

            t_408[k] = -ab_x[k] * gi_408[k]
                       + gk_520[k];

            t_409[k] = -ab_x[k] * gi_409[k]
                       + gk_521[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, gi_410, gi_411, gi_412, \
                         gi_413, gi_414, gk_522, gk_523, gk_524, gk_525, \
                         gk_526 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = -ab_x[k] * gi_410[k]
                       + gk_522[k];

            t_411[k] = -ab_x[k] * gi_411[k]
                       + gk_523[k];

            t_412[k] = -ab_x[k] * gi_412[k]
                       + gk_524[k];

            t_413[k] = -ab_x[k] * gi_413[k]
                       + gk_525[k];

            t_414[k] = -ab_x[k] * gi_414[k]
                       + gk_526[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, gi_415, gi_416, gi_417, \
                         gi_418, gi_419, gk_527, gk_528, gk_529, gk_530, \
                         gk_531 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = -ab_x[k] * gi_415[k]
                       + gk_527[k];

            t_416[k] = -ab_x[k] * gi_416[k]
                       + gk_528[k];

            t_417[k] = -ab_x[k] * gi_417[k]
                       + gk_529[k];

            t_418[k] = -ab_x[k] * gi_418[k]
                       + gk_530[k];

            t_419[k] = -ab_x[k] * gi_419[k]
                       + gk_531[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_y, gi_280, gi_281, gi_282, \
                         gi_283, gi_284, gk_361, gk_363, gk_364, gk_366, \
                         gk_367 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = -ab_y[k] * gi_280[k]
                       + gk_361[k];

            t_421[k] = -ab_y[k] * gi_281[k]
                       + gk_363[k];

            t_422[k] = -ab_y[k] * gi_282[k]
                       + gk_364[k];

            t_423[k] = -ab_y[k] * gi_283[k]
                       + gk_366[k];

            t_424[k] = -ab_y[k] * gi_284[k]
                       + gk_367[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_y, gi_285, gi_286, gi_287, \
                         gi_288, gi_289, gk_368, gk_370, gk_371, gk_372, \
                         gk_373 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = -ab_y[k] * gi_285[k]
                       + gk_368[k];

            t_426[k] = -ab_y[k] * gi_286[k]
                       + gk_370[k];

            t_427[k] = -ab_y[k] * gi_287[k]
                       + gk_371[k];

            t_428[k] = -ab_y[k] * gi_288[k]
                       + gk_372[k];

            t_429[k] = -ab_y[k] * gi_289[k]
                       + gk_373[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_y, gi_290, gi_291, gi_292, \
                         gi_293, gi_294, gk_375, gk_376, gk_377, gk_378, \
                         gk_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = -ab_y[k] * gi_290[k]
                       + gk_375[k];

            t_431[k] = -ab_y[k] * gi_291[k]
                       + gk_376[k];

            t_432[k] = -ab_y[k] * gi_292[k]
                       + gk_377[k];

            t_433[k] = -ab_y[k] * gi_293[k]
                       + gk_378[k];

            t_434[k] = -ab_y[k] * gi_294[k]
                       + gk_379[k];
        }
    }
}

static auto
compute_hrr_hi_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gi, const size_t gk, const size_t ncomps,
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gi_295 = buffer.data(gi + 295 * ncomps + c);
        const auto *gi_296 = buffer.data(gi + 296 * ncomps + c);
        const auto *gi_297 = buffer.data(gi + 297 * ncomps + c);
        const auto *gi_298 = buffer.data(gi + 298 * ncomps + c);
        const auto *gi_299 = buffer.data(gi + 299 * ncomps + c);
        const auto *gi_300 = buffer.data(gi + 300 * ncomps + c);
        const auto *gi_301 = buffer.data(gi + 301 * ncomps + c);
        const auto *gi_302 = buffer.data(gi + 302 * ncomps + c);
        const auto *gi_303 = buffer.data(gi + 303 * ncomps + c);
        const auto *gi_304 = buffer.data(gi + 304 * ncomps + c);
        const auto *gi_305 = buffer.data(gi + 305 * ncomps + c);
        const auto *gi_306 = buffer.data(gi + 306 * ncomps + c);
        const auto *gi_307 = buffer.data(gi + 307 * ncomps + c);
        const auto *gi_308 = buffer.data(gi + 308 * ncomps + c);
        const auto *gi_309 = buffer.data(gi + 309 * ncomps + c);
        const auto *gi_310 = buffer.data(gi + 310 * ncomps + c);
        const auto *gi_311 = buffer.data(gi + 311 * ncomps + c);
        const auto *gi_312 = buffer.data(gi + 312 * ncomps + c);
        const auto *gi_313 = buffer.data(gi + 313 * ncomps + c);
        const auto *gi_314 = buffer.data(gi + 314 * ncomps + c);
        const auto *gi_315 = buffer.data(gi + 315 * ncomps + c);
        const auto *gi_316 = buffer.data(gi + 316 * ncomps + c);
        const auto *gi_317 = buffer.data(gi + 317 * ncomps + c);
        const auto *gi_318 = buffer.data(gi + 318 * ncomps + c);
        const auto *gi_319 = buffer.data(gi + 319 * ncomps + c);
        const auto *gi_320 = buffer.data(gi + 320 * ncomps + c);
        const auto *gi_321 = buffer.data(gi + 321 * ncomps + c);
        const auto *gi_322 = buffer.data(gi + 322 * ncomps + c);
        const auto *gi_323 = buffer.data(gi + 323 * ncomps + c);
        const auto *gi_324 = buffer.data(gi + 324 * ncomps + c);
        const auto *gi_325 = buffer.data(gi + 325 * ncomps + c);
        const auto *gi_326 = buffer.data(gi + 326 * ncomps + c);
        const auto *gi_327 = buffer.data(gi + 327 * ncomps + c);
        const auto *gi_328 = buffer.data(gi + 328 * ncomps + c);
        const auto *gi_329 = buffer.data(gi + 329 * ncomps + c);
        const auto *gi_330 = buffer.data(gi + 330 * ncomps + c);
        const auto *gi_331 = buffer.data(gi + 331 * ncomps + c);
        const auto *gi_332 = buffer.data(gi + 332 * ncomps + c);
        const auto *gi_333 = buffer.data(gi + 333 * ncomps + c);
        const auto *gi_334 = buffer.data(gi + 334 * ncomps + c);
        const auto *gi_335 = buffer.data(gi + 335 * ncomps + c);
        const auto *gi_336 = buffer.data(gi + 336 * ncomps + c);
        const auto *gi_337 = buffer.data(gi + 337 * ncomps + c);
        const auto *gi_338 = buffer.data(gi + 338 * ncomps + c);
        const auto *gi_339 = buffer.data(gi + 339 * ncomps + c);
        const auto *gi_340 = buffer.data(gi + 340 * ncomps + c);
        const auto *gi_341 = buffer.data(gi + 341 * ncomps + c);
        const auto *gi_342 = buffer.data(gi + 342 * ncomps + c);
        const auto *gi_343 = buffer.data(gi + 343 * ncomps + c);
        const auto *gi_344 = buffer.data(gi + 344 * ncomps + c);
        const auto *gi_345 = buffer.data(gi + 345 * ncomps + c);
        const auto *gi_346 = buffer.data(gi + 346 * ncomps + c);
        const auto *gi_347 = buffer.data(gi + 347 * ncomps + c);
        const auto *gi_348 = buffer.data(gi + 348 * ncomps + c);
        const auto *gi_349 = buffer.data(gi + 349 * ncomps + c);
        const auto *gi_350 = buffer.data(gi + 350 * ncomps + c);
        const auto *gi_351 = buffer.data(gi + 351 * ncomps + c);
        const auto *gi_352 = buffer.data(gi + 352 * ncomps + c);
        const auto *gi_353 = buffer.data(gi + 353 * ncomps + c);
        const auto *gi_354 = buffer.data(gi + 354 * ncomps + c);
        const auto *gi_355 = buffer.data(gi + 355 * ncomps + c);
        const auto *gi_356 = buffer.data(gi + 356 * ncomps + c);
        const auto *gi_357 = buffer.data(gi + 357 * ncomps + c);
        const auto *gi_358 = buffer.data(gi + 358 * ncomps + c);
        const auto *gi_359 = buffer.data(gi + 359 * ncomps + c);
        const auto *gi_360 = buffer.data(gi + 360 * ncomps + c);
        const auto *gi_361 = buffer.data(gi + 361 * ncomps + c);
        const auto *gi_362 = buffer.data(gi + 362 * ncomps + c);
        const auto *gi_363 = buffer.data(gi + 363 * ncomps + c);
        const auto *gi_364 = buffer.data(gi + 364 * ncomps + c);
        const auto *gi_365 = buffer.data(gi + 365 * ncomps + c);
        const auto *gi_366 = buffer.data(gi + 366 * ncomps + c);
        const auto *gi_367 = buffer.data(gi + 367 * ncomps + c);
        const auto *gi_368 = buffer.data(gi + 368 * ncomps + c);
        const auto *gi_369 = buffer.data(gi + 369 * ncomps + c);
        const auto *gi_370 = buffer.data(gi + 370 * ncomps + c);
        const auto *gi_371 = buffer.data(gi + 371 * ncomps + c);
        const auto *gi_372 = buffer.data(gi + 372 * ncomps + c);
        const auto *gi_373 = buffer.data(gi + 373 * ncomps + c);
        const auto *gi_374 = buffer.data(gi + 374 * ncomps + c);
        const auto *gi_375 = buffer.data(gi + 375 * ncomps + c);
        const auto *gi_376 = buffer.data(gi + 376 * ncomps + c);
        const auto *gi_377 = buffer.data(gi + 377 * ncomps + c);
        const auto *gi_378 = buffer.data(gi + 378 * ncomps + c);
        const auto *gi_379 = buffer.data(gi + 379 * ncomps + c);
        const auto *gi_380 = buffer.data(gi + 380 * ncomps + c);
        const auto *gi_381 = buffer.data(gi + 381 * ncomps + c);
        const auto *gi_382 = buffer.data(gi + 382 * ncomps + c);
        const auto *gi_383 = buffer.data(gi + 383 * ncomps + c);
        const auto *gi_384 = buffer.data(gi + 384 * ncomps + c);
        const auto *gi_385 = buffer.data(gi + 385 * ncomps + c);
        const auto *gi_386 = buffer.data(gi + 386 * ncomps + c);
        const auto *gi_387 = buffer.data(gi + 387 * ncomps + c);
        const auto *gi_388 = buffer.data(gi + 388 * ncomps + c);
        const auto *gi_389 = buffer.data(gi + 389 * ncomps + c);
        const auto *gi_390 = buffer.data(gi + 390 * ncomps + c);
        const auto *gi_391 = buffer.data(gi + 391 * ncomps + c);
        const auto *gi_392 = buffer.data(gi + 392 * ncomps + c);
        const auto *gi_393 = buffer.data(gi + 393 * ncomps + c);
        const auto *gi_394 = buffer.data(gi + 394 * ncomps + c);
        const auto *gi_395 = buffer.data(gi + 395 * ncomps + c);
        const auto *gi_396 = buffer.data(gi + 396 * ncomps + c);
        const auto *gi_397 = buffer.data(gi + 397 * ncomps + c);
        const auto *gi_398 = buffer.data(gi + 398 * ncomps + c);
        const auto *gi_399 = buffer.data(gi + 399 * ncomps + c);
        const auto *gi_400 = buffer.data(gi + 400 * ncomps + c);
        const auto *gi_401 = buffer.data(gi + 401 * ncomps + c);
        const auto *gi_402 = buffer.data(gi + 402 * ncomps + c);
        const auto *gi_403 = buffer.data(gi + 403 * ncomps + c);
        const auto *gi_404 = buffer.data(gi + 404 * ncomps + c);
        const auto *gi_405 = buffer.data(gi + 405 * ncomps + c);
        const auto *gi_406 = buffer.data(gi + 406 * ncomps + c);
        const auto *gi_407 = buffer.data(gi + 407 * ncomps + c);
        const auto *gi_408 = buffer.data(gi + 408 * ncomps + c);
        const auto *gi_409 = buffer.data(gi + 409 * ncomps + c);
        const auto *gi_410 = buffer.data(gi + 410 * ncomps + c);
        const auto *gi_411 = buffer.data(gi + 411 * ncomps + c);
        const auto *gi_412 = buffer.data(gi + 412 * ncomps + c);
        const auto *gi_413 = buffer.data(gi + 413 * ncomps + c);
        const auto *gi_414 = buffer.data(gi + 414 * ncomps + c);
        const auto *gi_415 = buffer.data(gi + 415 * ncomps + c);
        const auto *gi_416 = buffer.data(gi + 416 * ncomps + c);
        const auto *gi_417 = buffer.data(gi + 417 * ncomps + c);
        const auto *gi_418 = buffer.data(gi + 418 * ncomps + c);
        const auto *gi_419 = buffer.data(gi + 419 * ncomps + c);

        const auto *gk_381 = buffer.data(gk + 381 * ncomps + c);
        const auto *gk_382 = buffer.data(gk + 382 * ncomps + c);
        const auto *gk_383 = buffer.data(gk + 383 * ncomps + c);
        const auto *gk_384 = buffer.data(gk + 384 * ncomps + c);
        const auto *gk_385 = buffer.data(gk + 385 * ncomps + c);
        const auto *gk_386 = buffer.data(gk + 386 * ncomps + c);
        const auto *gk_388 = buffer.data(gk + 388 * ncomps + c);
        const auto *gk_389 = buffer.data(gk + 389 * ncomps + c);
        const auto *gk_390 = buffer.data(gk + 390 * ncomps + c);
        const auto *gk_391 = buffer.data(gk + 391 * ncomps + c);
        const auto *gk_392 = buffer.data(gk + 392 * ncomps + c);
        const auto *gk_393 = buffer.data(gk + 393 * ncomps + c);
        const auto *gk_394 = buffer.data(gk + 394 * ncomps + c);
        const auto *gk_397 = buffer.data(gk + 397 * ncomps + c);
        const auto *gk_399 = buffer.data(gk + 399 * ncomps + c);
        const auto *gk_400 = buffer.data(gk + 400 * ncomps + c);
        const auto *gk_402 = buffer.data(gk + 402 * ncomps + c);
        const auto *gk_403 = buffer.data(gk + 403 * ncomps + c);
        const auto *gk_404 = buffer.data(gk + 404 * ncomps + c);
        const auto *gk_406 = buffer.data(gk + 406 * ncomps + c);
        const auto *gk_407 = buffer.data(gk + 407 * ncomps + c);
        const auto *gk_408 = buffer.data(gk + 408 * ncomps + c);
        const auto *gk_409 = buffer.data(gk + 409 * ncomps + c);
        const auto *gk_411 = buffer.data(gk + 411 * ncomps + c);
        const auto *gk_412 = buffer.data(gk + 412 * ncomps + c);
        const auto *gk_413 = buffer.data(gk + 413 * ncomps + c);
        const auto *gk_414 = buffer.data(gk + 414 * ncomps + c);
        const auto *gk_415 = buffer.data(gk + 415 * ncomps + c);
        const auto *gk_417 = buffer.data(gk + 417 * ncomps + c);
        const auto *gk_418 = buffer.data(gk + 418 * ncomps + c);
        const auto *gk_419 = buffer.data(gk + 419 * ncomps + c);
        const auto *gk_420 = buffer.data(gk + 420 * ncomps + c);
        const auto *gk_421 = buffer.data(gk + 421 * ncomps + c);
        const auto *gk_422 = buffer.data(gk + 422 * ncomps + c);
        const auto *gk_424 = buffer.data(gk + 424 * ncomps + c);
        const auto *gk_425 = buffer.data(gk + 425 * ncomps + c);
        const auto *gk_426 = buffer.data(gk + 426 * ncomps + c);
        const auto *gk_427 = buffer.data(gk + 427 * ncomps + c);
        const auto *gk_428 = buffer.data(gk + 428 * ncomps + c);
        const auto *gk_429 = buffer.data(gk + 429 * ncomps + c);
        const auto *gk_430 = buffer.data(gk + 430 * ncomps + c);
        const auto *gk_433 = buffer.data(gk + 433 * ncomps + c);
        const auto *gk_435 = buffer.data(gk + 435 * ncomps + c);
        const auto *gk_436 = buffer.data(gk + 436 * ncomps + c);
        const auto *gk_438 = buffer.data(gk + 438 * ncomps + c);
        const auto *gk_439 = buffer.data(gk + 439 * ncomps + c);
        const auto *gk_440 = buffer.data(gk + 440 * ncomps + c);
        const auto *gk_442 = buffer.data(gk + 442 * ncomps + c);
        const auto *gk_443 = buffer.data(gk + 443 * ncomps + c);
        const auto *gk_444 = buffer.data(gk + 444 * ncomps + c);
        const auto *gk_445 = buffer.data(gk + 445 * ncomps + c);
        const auto *gk_447 = buffer.data(gk + 447 * ncomps + c);
        const auto *gk_448 = buffer.data(gk + 448 * ncomps + c);
        const auto *gk_449 = buffer.data(gk + 449 * ncomps + c);
        const auto *gk_450 = buffer.data(gk + 450 * ncomps + c);
        const auto *gk_451 = buffer.data(gk + 451 * ncomps + c);
        const auto *gk_453 = buffer.data(gk + 453 * ncomps + c);
        const auto *gk_454 = buffer.data(gk + 454 * ncomps + c);
        const auto *gk_455 = buffer.data(gk + 455 * ncomps + c);
        const auto *gk_456 = buffer.data(gk + 456 * ncomps + c);
        const auto *gk_457 = buffer.data(gk + 457 * ncomps + c);
        const auto *gk_458 = buffer.data(gk + 458 * ncomps + c);
        const auto *gk_460 = buffer.data(gk + 460 * ncomps + c);
        const auto *gk_461 = buffer.data(gk + 461 * ncomps + c);
        const auto *gk_462 = buffer.data(gk + 462 * ncomps + c);
        const auto *gk_463 = buffer.data(gk + 463 * ncomps + c);
        const auto *gk_464 = buffer.data(gk + 464 * ncomps + c);
        const auto *gk_465 = buffer.data(gk + 465 * ncomps + c);
        const auto *gk_466 = buffer.data(gk + 466 * ncomps + c);
        const auto *gk_469 = buffer.data(gk + 469 * ncomps + c);
        const auto *gk_471 = buffer.data(gk + 471 * ncomps + c);
        const auto *gk_472 = buffer.data(gk + 472 * ncomps + c);
        const auto *gk_474 = buffer.data(gk + 474 * ncomps + c);
        const auto *gk_475 = buffer.data(gk + 475 * ncomps + c);
        const auto *gk_476 = buffer.data(gk + 476 * ncomps + c);
        const auto *gk_478 = buffer.data(gk + 478 * ncomps + c);
        const auto *gk_479 = buffer.data(gk + 479 * ncomps + c);
        const auto *gk_480 = buffer.data(gk + 480 * ncomps + c);
        const auto *gk_481 = buffer.data(gk + 481 * ncomps + c);
        const auto *gk_483 = buffer.data(gk + 483 * ncomps + c);
        const auto *gk_484 = buffer.data(gk + 484 * ncomps + c);
        const auto *gk_485 = buffer.data(gk + 485 * ncomps + c);
        const auto *gk_486 = buffer.data(gk + 486 * ncomps + c);
        const auto *gk_487 = buffer.data(gk + 487 * ncomps + c);
        const auto *gk_489 = buffer.data(gk + 489 * ncomps + c);
        const auto *gk_490 = buffer.data(gk + 490 * ncomps + c);
        const auto *gk_491 = buffer.data(gk + 491 * ncomps + c);
        const auto *gk_492 = buffer.data(gk + 492 * ncomps + c);
        const auto *gk_493 = buffer.data(gk + 493 * ncomps + c);
        const auto *gk_494 = buffer.data(gk + 494 * ncomps + c);
        const auto *gk_496 = buffer.data(gk + 496 * ncomps + c);
        const auto *gk_497 = buffer.data(gk + 497 * ncomps + c);
        const auto *gk_498 = buffer.data(gk + 498 * ncomps + c);
        const auto *gk_499 = buffer.data(gk + 499 * ncomps + c);
        const auto *gk_500 = buffer.data(gk + 500 * ncomps + c);
        const auto *gk_501 = buffer.data(gk + 501 * ncomps + c);
        const auto *gk_502 = buffer.data(gk + 502 * ncomps + c);
        const auto *gk_505 = buffer.data(gk + 505 * ncomps + c);
        const auto *gk_506 = buffer.data(gk + 506 * ncomps + c);
        const auto *gk_507 = buffer.data(gk + 507 * ncomps + c);
        const auto *gk_508 = buffer.data(gk + 508 * ncomps + c);
        const auto *gk_509 = buffer.data(gk + 509 * ncomps + c);
        const auto *gk_510 = buffer.data(gk + 510 * ncomps + c);
        const auto *gk_511 = buffer.data(gk + 511 * ncomps + c);
        const auto *gk_512 = buffer.data(gk + 512 * ncomps + c);
        const auto *gk_513 = buffer.data(gk + 513 * ncomps + c);
        const auto *gk_514 = buffer.data(gk + 514 * ncomps + c);
        const auto *gk_515 = buffer.data(gk + 515 * ncomps + c);
        const auto *gk_516 = buffer.data(gk + 516 * ncomps + c);
        const auto *gk_517 = buffer.data(gk + 517 * ncomps + c);
        const auto *gk_518 = buffer.data(gk + 518 * ncomps + c);
        const auto *gk_519 = buffer.data(gk + 519 * ncomps + c);
        const auto *gk_520 = buffer.data(gk + 520 * ncomps + c);
        const auto *gk_521 = buffer.data(gk + 521 * ncomps + c);
        const auto *gk_522 = buffer.data(gk + 522 * ncomps + c);
        const auto *gk_523 = buffer.data(gk + 523 * ncomps + c);
        const auto *gk_524 = buffer.data(gk + 524 * ncomps + c);
        const auto *gk_525 = buffer.data(gk + 525 * ncomps + c);
        const auto *gk_526 = buffer.data(gk + 526 * ncomps + c);
        const auto *gk_527 = buffer.data(gk + 527 * ncomps + c);
        const auto *gk_528 = buffer.data(gk + 528 * ncomps + c);
        const auto *gk_529 = buffer.data(gk + 529 * ncomps + c);
        const auto *gk_530 = buffer.data(gk + 530 * ncomps + c);
        const auto *gk_532 = buffer.data(gk + 532 * ncomps + c);
        const auto *gk_533 = buffer.data(gk + 533 * ncomps + c);
        const auto *gk_534 = buffer.data(gk + 534 * ncomps + c);
        const auto *gk_535 = buffer.data(gk + 535 * ncomps + c);
        const auto *gk_536 = buffer.data(gk + 536 * ncomps + c);
        const auto *gk_537 = buffer.data(gk + 537 * ncomps + c);
        const auto *gk_538 = buffer.data(gk + 538 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, gi_295, gi_296, gi_297, \
                         gi_298, gi_299, gk_381, gk_382, gk_383, gk_384, \
                         gk_385 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = -ab_y[k] * gi_295[k]
                       + gk_381[k];

            t_436[k] = -ab_y[k] * gi_296[k]
                       + gk_382[k];

            t_437[k] = -ab_y[k] * gi_297[k]
                       + gk_383[k];

            t_438[k] = -ab_y[k] * gi_298[k]
                       + gk_384[k];

            t_439[k] = -ab_y[k] * gi_299[k]
                       + gk_385[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_y, gi_300, gi_301, gi_302, \
                         gi_303, gi_304, gk_386, gk_388, gk_389, gk_390, \
                         gk_391 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = -ab_y[k] * gi_300[k]
                       + gk_386[k];

            t_441[k] = -ab_y[k] * gi_301[k]
                       + gk_388[k];

            t_442[k] = -ab_y[k] * gi_302[k]
                       + gk_389[k];

            t_443[k] = -ab_y[k] * gi_303[k]
                       + gk_390[k];

            t_444[k] = -ab_y[k] * gi_304[k]
                       + gk_391[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_y, gi_305, gi_306, gi_307, \
                         gi_308, gi_309, gk_392, gk_393, gk_394, gk_397, \
                         gk_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = -ab_y[k] * gi_305[k]
                       + gk_392[k];

            t_446[k] = -ab_y[k] * gi_306[k]
                       + gk_393[k];

            t_447[k] = -ab_y[k] * gi_307[k]
                       + gk_394[k];

            t_448[k] = -ab_y[k] * gi_308[k]
                       + gk_397[k];

            t_449[k] = -ab_y[k] * gi_309[k]
                       + gk_399[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_y, gi_310, gi_311, gi_312, \
                         gi_313, gi_314, gk_400, gk_402, gk_403, gk_404, \
                         gk_406 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = -ab_y[k] * gi_310[k]
                       + gk_400[k];

            t_451[k] = -ab_y[k] * gi_311[k]
                       + gk_402[k];

            t_452[k] = -ab_y[k] * gi_312[k]
                       + gk_403[k];

            t_453[k] = -ab_y[k] * gi_313[k]
                       + gk_404[k];

            t_454[k] = -ab_y[k] * gi_314[k]
                       + gk_406[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_y, gi_315, gi_316, gi_317, \
                         gi_318, gi_319, gk_407, gk_408, gk_409, gk_411, \
                         gk_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = -ab_y[k] * gi_315[k]
                       + gk_407[k];

            t_456[k] = -ab_y[k] * gi_316[k]
                       + gk_408[k];

            t_457[k] = -ab_y[k] * gi_317[k]
                       + gk_409[k];

            t_458[k] = -ab_y[k] * gi_318[k]
                       + gk_411[k];

            t_459[k] = -ab_y[k] * gi_319[k]
                       + gk_412[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, gi_320, gi_321, gi_322, \
                         gi_323, gi_324, gk_413, gk_414, gk_415, gk_417, \
                         gk_418 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = -ab_y[k] * gi_320[k]
                       + gk_413[k];

            t_461[k] = -ab_y[k] * gi_321[k]
                       + gk_414[k];

            t_462[k] = -ab_y[k] * gi_322[k]
                       + gk_415[k];

            t_463[k] = -ab_y[k] * gi_323[k]
                       + gk_417[k];

            t_464[k] = -ab_y[k] * gi_324[k]
                       + gk_418[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_y, gi_325, gi_326, gi_327, \
                         gi_328, gi_329, gk_419, gk_420, gk_421, gk_422, \
                         gk_424 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = -ab_y[k] * gi_325[k]
                       + gk_419[k];

            t_466[k] = -ab_y[k] * gi_326[k]
                       + gk_420[k];

            t_467[k] = -ab_y[k] * gi_327[k]
                       + gk_421[k];

            t_468[k] = -ab_y[k] * gi_328[k]
                       + gk_422[k];

            t_469[k] = -ab_y[k] * gi_329[k]
                       + gk_424[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_y, gi_330, gi_331, gi_332, \
                         gi_333, gi_334, gk_425, gk_426, gk_427, gk_428, \
                         gk_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = -ab_y[k] * gi_330[k]
                       + gk_425[k];

            t_471[k] = -ab_y[k] * gi_331[k]
                       + gk_426[k];

            t_472[k] = -ab_y[k] * gi_332[k]
                       + gk_427[k];

            t_473[k] = -ab_y[k] * gi_333[k]
                       + gk_428[k];

            t_474[k] = -ab_y[k] * gi_334[k]
                       + gk_429[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, gi_335, gi_336, gi_337, \
                         gi_338, gi_339, gk_430, gk_433, gk_435, gk_436, \
                         gk_438 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = -ab_y[k] * gi_335[k]
                       + gk_430[k];

            t_476[k] = -ab_y[k] * gi_336[k]
                       + gk_433[k];

            t_477[k] = -ab_y[k] * gi_337[k]
                       + gk_435[k];

            t_478[k] = -ab_y[k] * gi_338[k]
                       + gk_436[k];

            t_479[k] = -ab_y[k] * gi_339[k]
                       + gk_438[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_y, gi_340, gi_341, gi_342, \
                         gi_343, gi_344, gk_439, gk_440, gk_442, gk_443, \
                         gk_444 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = -ab_y[k] * gi_340[k]
                       + gk_439[k];

            t_481[k] = -ab_y[k] * gi_341[k]
                       + gk_440[k];

            t_482[k] = -ab_y[k] * gi_342[k]
                       + gk_442[k];

            t_483[k] = -ab_y[k] * gi_343[k]
                       + gk_443[k];

            t_484[k] = -ab_y[k] * gi_344[k]
                       + gk_444[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_y, gi_345, gi_346, gi_347, \
                         gi_348, gi_349, gk_445, gk_447, gk_448, gk_449, \
                         gk_450 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = -ab_y[k] * gi_345[k]
                       + gk_445[k];

            t_486[k] = -ab_y[k] * gi_346[k]
                       + gk_447[k];

            t_487[k] = -ab_y[k] * gi_347[k]
                       + gk_448[k];

            t_488[k] = -ab_y[k] * gi_348[k]
                       + gk_449[k];

            t_489[k] = -ab_y[k] * gi_349[k]
                       + gk_450[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, gi_350, gi_351, gi_352, \
                         gi_353, gi_354, gk_451, gk_453, gk_454, gk_455, \
                         gk_456 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = -ab_y[k] * gi_350[k]
                       + gk_451[k];

            t_491[k] = -ab_y[k] * gi_351[k]
                       + gk_453[k];

            t_492[k] = -ab_y[k] * gi_352[k]
                       + gk_454[k];

            t_493[k] = -ab_y[k] * gi_353[k]
                       + gk_455[k];

            t_494[k] = -ab_y[k] * gi_354[k]
                       + gk_456[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_y, gi_355, gi_356, gi_357, \
                         gi_358, gi_359, gk_457, gk_458, gk_460, gk_461, \
                         gk_462 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = -ab_y[k] * gi_355[k]
                       + gk_457[k];

            t_496[k] = -ab_y[k] * gi_356[k]
                       + gk_458[k];

            t_497[k] = -ab_y[k] * gi_357[k]
                       + gk_460[k];

            t_498[k] = -ab_y[k] * gi_358[k]
                       + gk_461[k];

            t_499[k] = -ab_y[k] * gi_359[k]
                       + gk_462[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_y, gi_360, gi_361, gi_362, \
                         gi_363, gi_364, gk_463, gk_464, gk_465, gk_466, \
                         gk_469 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = -ab_y[k] * gi_360[k]
                       + gk_463[k];

            t_501[k] = -ab_y[k] * gi_361[k]
                       + gk_464[k];

            t_502[k] = -ab_y[k] * gi_362[k]
                       + gk_465[k];

            t_503[k] = -ab_y[k] * gi_363[k]
                       + gk_466[k];

            t_504[k] = -ab_y[k] * gi_364[k]
                       + gk_469[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_y, gi_365, gi_366, gi_367, \
                         gi_368, gi_369, gk_471, gk_472, gk_474, gk_475, \
                         gk_476 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = -ab_y[k] * gi_365[k]
                       + gk_471[k];

            t_506[k] = -ab_y[k] * gi_366[k]
                       + gk_472[k];

            t_507[k] = -ab_y[k] * gi_367[k]
                       + gk_474[k];

            t_508[k] = -ab_y[k] * gi_368[k]
                       + gk_475[k];

            t_509[k] = -ab_y[k] * gi_369[k]
                       + gk_476[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_y, gi_370, gi_371, gi_372, \
                         gi_373, gi_374, gk_478, gk_479, gk_480, gk_481, \
                         gk_483 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = -ab_y[k] * gi_370[k]
                       + gk_478[k];

            t_511[k] = -ab_y[k] * gi_371[k]
                       + gk_479[k];

            t_512[k] = -ab_y[k] * gi_372[k]
                       + gk_480[k];

            t_513[k] = -ab_y[k] * gi_373[k]
                       + gk_481[k];

            t_514[k] = -ab_y[k] * gi_374[k]
                       + gk_483[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_y, gi_375, gi_376, gi_377, \
                         gi_378, gi_379, gk_484, gk_485, gk_486, gk_487, \
                         gk_489 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = -ab_y[k] * gi_375[k]
                       + gk_484[k];

            t_516[k] = -ab_y[k] * gi_376[k]
                       + gk_485[k];

            t_517[k] = -ab_y[k] * gi_377[k]
                       + gk_486[k];

            t_518[k] = -ab_y[k] * gi_378[k]
                       + gk_487[k];

            t_519[k] = -ab_y[k] * gi_379[k]
                       + gk_489[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_y, gi_380, gi_381, gi_382, \
                         gi_383, gi_384, gk_490, gk_491, gk_492, gk_493, \
                         gk_494 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = -ab_y[k] * gi_380[k]
                       + gk_490[k];

            t_521[k] = -ab_y[k] * gi_381[k]
                       + gk_491[k];

            t_522[k] = -ab_y[k] * gi_382[k]
                       + gk_492[k];

            t_523[k] = -ab_y[k] * gi_383[k]
                       + gk_493[k];

            t_524[k] = -ab_y[k] * gi_384[k]
                       + gk_494[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_y, gi_385, gi_386, gi_387, \
                         gi_388, gi_389, gk_496, gk_497, gk_498, gk_499, \
                         gk_500 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = -ab_y[k] * gi_385[k]
                       + gk_496[k];

            t_526[k] = -ab_y[k] * gi_386[k]
                       + gk_497[k];

            t_527[k] = -ab_y[k] * gi_387[k]
                       + gk_498[k];

            t_528[k] = -ab_y[k] * gi_388[k]
                       + gk_499[k];

            t_529[k] = -ab_y[k] * gi_389[k]
                       + gk_500[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_y, gi_390, gi_391, gi_392, \
                         gi_393, gi_394, gk_501, gk_502, gk_505, gk_507, \
                         gk_508 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = -ab_y[k] * gi_390[k]
                       + gk_501[k];

            t_531[k] = -ab_y[k] * gi_391[k]
                       + gk_502[k];

            t_532[k] = -ab_y[k] * gi_392[k]
                       + gk_505[k];

            t_533[k] = -ab_y[k] * gi_393[k]
                       + gk_507[k];

            t_534[k] = -ab_y[k] * gi_394[k]
                       + gk_508[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_y, gi_395, gi_396, gi_397, \
                         gi_398, gi_399, gk_510, gk_511, gk_512, gk_514, \
                         gk_515 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = -ab_y[k] * gi_395[k]
                       + gk_510[k];

            t_536[k] = -ab_y[k] * gi_396[k]
                       + gk_511[k];

            t_537[k] = -ab_y[k] * gi_397[k]
                       + gk_512[k];

            t_538[k] = -ab_y[k] * gi_398[k]
                       + gk_514[k];

            t_539[k] = -ab_y[k] * gi_399[k]
                       + gk_515[k];
        }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_y, gi_400, gi_401, gi_402, \
                         gi_403, gi_404, gk_516, gk_517, gk_519, gk_520, \
                         gk_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_540[k] = -ab_y[k] * gi_400[k]
                       + gk_516[k];

            t_541[k] = -ab_y[k] * gi_401[k]
                       + gk_517[k];

            t_542[k] = -ab_y[k] * gi_402[k]
                       + gk_519[k];

            t_543[k] = -ab_y[k] * gi_403[k]
                       + gk_520[k];

            t_544[k] = -ab_y[k] * gi_404[k]
                       + gk_521[k];
        }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_y, gi_405, gi_406, gi_407, \
                         gi_408, gi_409, gk_522, gk_523, gk_525, gk_526, \
                         gk_527 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_545[k] = -ab_y[k] * gi_405[k]
                       + gk_522[k];

            t_546[k] = -ab_y[k] * gi_406[k]
                       + gk_523[k];

            t_547[k] = -ab_y[k] * gi_407[k]
                       + gk_525[k];

            t_548[k] = -ab_y[k] * gi_408[k]
                       + gk_526[k];

            t_549[k] = -ab_y[k] * gi_409[k]
                       + gk_527[k];
        }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ab_y, gi_410, gi_411, gi_412, \
                         gi_413, gi_414, gk_528, gk_529, gk_530, gk_532, \
                         gk_533 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_550[k] = -ab_y[k] * gi_410[k]
                       + gk_528[k];

            t_551[k] = -ab_y[k] * gi_411[k]
                       + gk_529[k];

            t_552[k] = -ab_y[k] * gi_412[k]
                       + gk_530[k];

            t_553[k] = -ab_y[k] * gi_413[k]
                       + gk_532[k];

            t_554[k] = -ab_y[k] * gi_414[k]
                       + gk_533[k];
        }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ab_y, gi_415, gi_416, gi_417, \
                         gi_418, gi_419, gk_534, gk_535, gk_536, gk_537, \
                         gk_538 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_555[k] = -ab_y[k] * gi_415[k]
                       + gk_534[k];

            t_556[k] = -ab_y[k] * gi_416[k]
                       + gk_535[k];

            t_557[k] = -ab_y[k] * gi_417[k]
                       + gk_536[k];

            t_558[k] = -ab_y[k] * gi_418[k]
                       + gk_537[k];

            t_559[k] = -ab_y[k] * gi_419[k]
                       + gk_538[k];
        }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_z, gi_392, gi_393, gi_394, \
                         gi_395, gi_396, gk_506, gk_508, gk_509, gk_511, \
                         gk_512 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_560[k] = -ab_z[k] * gi_392[k]
                       + gk_506[k];

            t_561[k] = -ab_z[k] * gi_393[k]
                       + gk_508[k];

            t_562[k] = -ab_z[k] * gi_394[k]
                       + gk_509[k];

            t_563[k] = -ab_z[k] * gi_395[k]
                       + gk_511[k];

            t_564[k] = -ab_z[k] * gi_396[k]
                       + gk_512[k];
        }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_z, gi_397, gi_398, gi_399, \
                         gi_400, gi_401, gk_513, gk_515, gk_516, gk_517, \
                         gk_518 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_565[k] = -ab_z[k] * gi_397[k]
                       + gk_513[k];

            t_566[k] = -ab_z[k] * gi_398[k]
                       + gk_515[k];

            t_567[k] = -ab_z[k] * gi_399[k]
                       + gk_516[k];

            t_568[k] = -ab_z[k] * gi_400[k]
                       + gk_517[k];

            t_569[k] = -ab_z[k] * gi_401[k]
                       + gk_518[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_z, gi_402, gi_403, gi_404, \
                         gi_405, gi_406, gk_520, gk_521, gk_522, gk_523, \
                         gk_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = -ab_z[k] * gi_402[k]
                       + gk_520[k];

            t_571[k] = -ab_z[k] * gi_403[k]
                       + gk_521[k];

            t_572[k] = -ab_z[k] * gi_404[k]
                       + gk_522[k];

            t_573[k] = -ab_z[k] * gi_405[k]
                       + gk_523[k];

            t_574[k] = -ab_z[k] * gi_406[k]
                       + gk_524[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_z, gi_407, gi_408, gi_409, \
                         gi_410, gi_411, gk_526, gk_527, gk_528, gk_529, \
                         gk_530 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = -ab_z[k] * gi_407[k]
                       + gk_526[k];

            t_576[k] = -ab_z[k] * gi_408[k]
                       + gk_527[k];

            t_577[k] = -ab_z[k] * gi_409[k]
                       + gk_528[k];

            t_578[k] = -ab_z[k] * gi_410[k]
                       + gk_529[k];

            t_579[k] = -ab_z[k] * gi_411[k]
                       + gk_530[k];
        }
    }
}

static auto
compute_hrr_hi_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gi, const size_t gk, const size_t ncomps,
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

        const auto *ab_z = coordinates.data(8);

        const auto *gi_412 = buffer.data(gi + 412 * ncomps + c);
        const auto *gi_413 = buffer.data(gi + 413 * ncomps + c);
        const auto *gi_414 = buffer.data(gi + 414 * ncomps + c);
        const auto *gi_415 = buffer.data(gi + 415 * ncomps + c);
        const auto *gi_416 = buffer.data(gi + 416 * ncomps + c);
        const auto *gi_417 = buffer.data(gi + 417 * ncomps + c);
        const auto *gi_418 = buffer.data(gi + 418 * ncomps + c);
        const auto *gi_419 = buffer.data(gi + 419 * ncomps + c);

        const auto *gk_531 = buffer.data(gk + 531 * ncomps + c);
        const auto *gk_533 = buffer.data(gk + 533 * ncomps + c);
        const auto *gk_534 = buffer.data(gk + 534 * ncomps + c);
        const auto *gk_535 = buffer.data(gk + 535 * ncomps + c);
        const auto *gk_536 = buffer.data(gk + 536 * ncomps + c);
        const auto *gk_537 = buffer.data(gk + 537 * ncomps + c);
        const auto *gk_538 = buffer.data(gk + 538 * ncomps + c);
        const auto *gk_539 = buffer.data(gk + 539 * ncomps + c);

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ab_z, gi_412, gi_413, gi_414, \
                         gi_415, gi_416, gk_531, gk_533, gk_534, gk_535, \
                         gk_536 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = -ab_z[k] * gi_412[k]
                       + gk_531[k];

            t_581[k] = -ab_z[k] * gi_413[k]
                       + gk_533[k];

            t_582[k] = -ab_z[k] * gi_414[k]
                       + gk_534[k];

            t_583[k] = -ab_z[k] * gi_415[k]
                       + gk_535[k];

            t_584[k] = -ab_z[k] * gi_416[k]
                       + gk_536[k];
        }

#pragma omp simd aligned(t_585, t_586, t_587, ab_z, gi_417, gi_418, gi_419, gk_537, gk_538, \
                         gk_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_585[k] = -ab_z[k] * gi_417[k]
                       + gk_537[k];

            t_586[k] = -ab_z[k] * gi_418[k]
                       + gk_538[k];

            t_587[k] = -ab_z[k] * gi_419[k]
                       + gk_539[k];
        }
    }
}

auto
compute_hrr_hi(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gi, const size_t gk, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_hi_piece0(buffer, coordinates, target, gi, gk, ncomps, nmax);

    compute_hrr_hi_piece1(buffer, coordinates, target, gi, gk, ncomps, nmax);

    compute_hrr_hi_piece2(buffer, coordinates, target, gi, gk, ncomps, nmax);

    compute_hrr_hi_piece3(buffer, coordinates, target, gi, gk, ncomps, nmax);

    compute_hrr_hi_piece4(buffer, coordinates, target, gi, gk, ncomps, nmax);
}

}  // namespace simdtrf
