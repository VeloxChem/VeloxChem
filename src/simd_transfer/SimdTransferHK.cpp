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


#include "SimdTransferHK.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_hk_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gk, const size_t gl, const size_t ncomps,
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
        const auto *gk_28 = buffer.data(gk + 28 * ncomps + c);
        const auto *gk_29 = buffer.data(gk + 29 * ncomps + c);
        const auto *gk_30 = buffer.data(gk + 30 * ncomps + c);
        const auto *gk_31 = buffer.data(gk + 31 * ncomps + c);
        const auto *gk_32 = buffer.data(gk + 32 * ncomps + c);
        const auto *gk_33 = buffer.data(gk + 33 * ncomps + c);
        const auto *gk_34 = buffer.data(gk + 34 * ncomps + c);
        const auto *gk_35 = buffer.data(gk + 35 * ncomps + c);
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
        const auto *gk_64 = buffer.data(gk + 64 * ncomps + c);
        const auto *gk_65 = buffer.data(gk + 65 * ncomps + c);
        const auto *gk_66 = buffer.data(gk + 66 * ncomps + c);
        const auto *gk_67 = buffer.data(gk + 67 * ncomps + c);
        const auto *gk_68 = buffer.data(gk + 68 * ncomps + c);
        const auto *gk_69 = buffer.data(gk + 69 * ncomps + c);
        const auto *gk_70 = buffer.data(gk + 70 * ncomps + c);
        const auto *gk_71 = buffer.data(gk + 71 * ncomps + c);
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
        const auto *gk_100 = buffer.data(gk + 100 * ncomps + c);
        const auto *gk_101 = buffer.data(gk + 101 * ncomps + c);
        const auto *gk_102 = buffer.data(gk + 102 * ncomps + c);
        const auto *gk_103 = buffer.data(gk + 103 * ncomps + c);
        const auto *gk_104 = buffer.data(gk + 104 * ncomps + c);
        const auto *gk_105 = buffer.data(gk + 105 * ncomps + c);
        const auto *gk_106 = buffer.data(gk + 106 * ncomps + c);
        const auto *gk_107 = buffer.data(gk + 107 * ncomps + c);
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
        const auto *gk_136 = buffer.data(gk + 136 * ncomps + c);
        const auto *gk_137 = buffer.data(gk + 137 * ncomps + c);
        const auto *gk_138 = buffer.data(gk + 138 * ncomps + c);
        const auto *gk_139 = buffer.data(gk + 139 * ncomps + c);
        const auto *gk_140 = buffer.data(gk + 140 * ncomps + c);
        const auto *gk_141 = buffer.data(gk + 141 * ncomps + c);
        const auto *gk_142 = buffer.data(gk + 142 * ncomps + c);
        const auto *gk_143 = buffer.data(gk + 143 * ncomps + c);
        const auto *gk_144 = buffer.data(gk + 144 * ncomps + c);

        const auto *gl_0 = buffer.data(gl + 0 * ncomps + c);
        const auto *gl_1 = buffer.data(gl + 1 * ncomps + c);
        const auto *gl_2 = buffer.data(gl + 2 * ncomps + c);
        const auto *gl_3 = buffer.data(gl + 3 * ncomps + c);
        const auto *gl_4 = buffer.data(gl + 4 * ncomps + c);
        const auto *gl_5 = buffer.data(gl + 5 * ncomps + c);
        const auto *gl_6 = buffer.data(gl + 6 * ncomps + c);
        const auto *gl_7 = buffer.data(gl + 7 * ncomps + c);
        const auto *gl_8 = buffer.data(gl + 8 * ncomps + c);
        const auto *gl_9 = buffer.data(gl + 9 * ncomps + c);
        const auto *gl_10 = buffer.data(gl + 10 * ncomps + c);
        const auto *gl_11 = buffer.data(gl + 11 * ncomps + c);
        const auto *gl_12 = buffer.data(gl + 12 * ncomps + c);
        const auto *gl_13 = buffer.data(gl + 13 * ncomps + c);
        const auto *gl_14 = buffer.data(gl + 14 * ncomps + c);
        const auto *gl_15 = buffer.data(gl + 15 * ncomps + c);
        const auto *gl_16 = buffer.data(gl + 16 * ncomps + c);
        const auto *gl_17 = buffer.data(gl + 17 * ncomps + c);
        const auto *gl_18 = buffer.data(gl + 18 * ncomps + c);
        const auto *gl_19 = buffer.data(gl + 19 * ncomps + c);
        const auto *gl_20 = buffer.data(gl + 20 * ncomps + c);
        const auto *gl_21 = buffer.data(gl + 21 * ncomps + c);
        const auto *gl_22 = buffer.data(gl + 22 * ncomps + c);
        const auto *gl_23 = buffer.data(gl + 23 * ncomps + c);
        const auto *gl_24 = buffer.data(gl + 24 * ncomps + c);
        const auto *gl_25 = buffer.data(gl + 25 * ncomps + c);
        const auto *gl_26 = buffer.data(gl + 26 * ncomps + c);
        const auto *gl_27 = buffer.data(gl + 27 * ncomps + c);
        const auto *gl_28 = buffer.data(gl + 28 * ncomps + c);
        const auto *gl_29 = buffer.data(gl + 29 * ncomps + c);
        const auto *gl_30 = buffer.data(gl + 30 * ncomps + c);
        const auto *gl_31 = buffer.data(gl + 31 * ncomps + c);
        const auto *gl_32 = buffer.data(gl + 32 * ncomps + c);
        const auto *gl_33 = buffer.data(gl + 33 * ncomps + c);
        const auto *gl_34 = buffer.data(gl + 34 * ncomps + c);
        const auto *gl_35 = buffer.data(gl + 35 * ncomps + c);
        const auto *gl_45 = buffer.data(gl + 45 * ncomps + c);
        const auto *gl_46 = buffer.data(gl + 46 * ncomps + c);
        const auto *gl_47 = buffer.data(gl + 47 * ncomps + c);
        const auto *gl_48 = buffer.data(gl + 48 * ncomps + c);
        const auto *gl_49 = buffer.data(gl + 49 * ncomps + c);
        const auto *gl_50 = buffer.data(gl + 50 * ncomps + c);
        const auto *gl_51 = buffer.data(gl + 51 * ncomps + c);
        const auto *gl_52 = buffer.data(gl + 52 * ncomps + c);
        const auto *gl_53 = buffer.data(gl + 53 * ncomps + c);
        const auto *gl_54 = buffer.data(gl + 54 * ncomps + c);
        const auto *gl_55 = buffer.data(gl + 55 * ncomps + c);
        const auto *gl_56 = buffer.data(gl + 56 * ncomps + c);
        const auto *gl_57 = buffer.data(gl + 57 * ncomps + c);
        const auto *gl_58 = buffer.data(gl + 58 * ncomps + c);
        const auto *gl_59 = buffer.data(gl + 59 * ncomps + c);
        const auto *gl_60 = buffer.data(gl + 60 * ncomps + c);
        const auto *gl_61 = buffer.data(gl + 61 * ncomps + c);
        const auto *gl_62 = buffer.data(gl + 62 * ncomps + c);
        const auto *gl_63 = buffer.data(gl + 63 * ncomps + c);
        const auto *gl_64 = buffer.data(gl + 64 * ncomps + c);
        const auto *gl_65 = buffer.data(gl + 65 * ncomps + c);
        const auto *gl_66 = buffer.data(gl + 66 * ncomps + c);
        const auto *gl_67 = buffer.data(gl + 67 * ncomps + c);
        const auto *gl_68 = buffer.data(gl + 68 * ncomps + c);
        const auto *gl_69 = buffer.data(gl + 69 * ncomps + c);
        const auto *gl_70 = buffer.data(gl + 70 * ncomps + c);
        const auto *gl_71 = buffer.data(gl + 71 * ncomps + c);
        const auto *gl_72 = buffer.data(gl + 72 * ncomps + c);
        const auto *gl_73 = buffer.data(gl + 73 * ncomps + c);
        const auto *gl_74 = buffer.data(gl + 74 * ncomps + c);
        const auto *gl_75 = buffer.data(gl + 75 * ncomps + c);
        const auto *gl_76 = buffer.data(gl + 76 * ncomps + c);
        const auto *gl_77 = buffer.data(gl + 77 * ncomps + c);
        const auto *gl_78 = buffer.data(gl + 78 * ncomps + c);
        const auto *gl_79 = buffer.data(gl + 79 * ncomps + c);
        const auto *gl_80 = buffer.data(gl + 80 * ncomps + c);
        const auto *gl_90 = buffer.data(gl + 90 * ncomps + c);
        const auto *gl_91 = buffer.data(gl + 91 * ncomps + c);
        const auto *gl_92 = buffer.data(gl + 92 * ncomps + c);
        const auto *gl_93 = buffer.data(gl + 93 * ncomps + c);
        const auto *gl_94 = buffer.data(gl + 94 * ncomps + c);
        const auto *gl_95 = buffer.data(gl + 95 * ncomps + c);
        const auto *gl_96 = buffer.data(gl + 96 * ncomps + c);
        const auto *gl_97 = buffer.data(gl + 97 * ncomps + c);
        const auto *gl_98 = buffer.data(gl + 98 * ncomps + c);
        const auto *gl_99 = buffer.data(gl + 99 * ncomps + c);
        const auto *gl_100 = buffer.data(gl + 100 * ncomps + c);
        const auto *gl_101 = buffer.data(gl + 101 * ncomps + c);
        const auto *gl_102 = buffer.data(gl + 102 * ncomps + c);
        const auto *gl_103 = buffer.data(gl + 103 * ncomps + c);
        const auto *gl_104 = buffer.data(gl + 104 * ncomps + c);
        const auto *gl_105 = buffer.data(gl + 105 * ncomps + c);
        const auto *gl_106 = buffer.data(gl + 106 * ncomps + c);
        const auto *gl_107 = buffer.data(gl + 107 * ncomps + c);
        const auto *gl_108 = buffer.data(gl + 108 * ncomps + c);
        const auto *gl_109 = buffer.data(gl + 109 * ncomps + c);
        const auto *gl_110 = buffer.data(gl + 110 * ncomps + c);
        const auto *gl_111 = buffer.data(gl + 111 * ncomps + c);
        const auto *gl_112 = buffer.data(gl + 112 * ncomps + c);
        const auto *gl_113 = buffer.data(gl + 113 * ncomps + c);
        const auto *gl_114 = buffer.data(gl + 114 * ncomps + c);
        const auto *gl_115 = buffer.data(gl + 115 * ncomps + c);
        const auto *gl_116 = buffer.data(gl + 116 * ncomps + c);
        const auto *gl_117 = buffer.data(gl + 117 * ncomps + c);
        const auto *gl_118 = buffer.data(gl + 118 * ncomps + c);
        const auto *gl_119 = buffer.data(gl + 119 * ncomps + c);
        const auto *gl_120 = buffer.data(gl + 120 * ncomps + c);
        const auto *gl_121 = buffer.data(gl + 121 * ncomps + c);
        const auto *gl_122 = buffer.data(gl + 122 * ncomps + c);
        const auto *gl_123 = buffer.data(gl + 123 * ncomps + c);
        const auto *gl_124 = buffer.data(gl + 124 * ncomps + c);
        const auto *gl_125 = buffer.data(gl + 125 * ncomps + c);
        const auto *gl_135 = buffer.data(gl + 135 * ncomps + c);
        const auto *gl_136 = buffer.data(gl + 136 * ncomps + c);
        const auto *gl_137 = buffer.data(gl + 137 * ncomps + c);
        const auto *gl_138 = buffer.data(gl + 138 * ncomps + c);
        const auto *gl_139 = buffer.data(gl + 139 * ncomps + c);
        const auto *gl_140 = buffer.data(gl + 140 * ncomps + c);
        const auto *gl_141 = buffer.data(gl + 141 * ncomps + c);
        const auto *gl_142 = buffer.data(gl + 142 * ncomps + c);
        const auto *gl_143 = buffer.data(gl + 143 * ncomps + c);
        const auto *gl_144 = buffer.data(gl + 144 * ncomps + c);
        const auto *gl_145 = buffer.data(gl + 145 * ncomps + c);
        const auto *gl_146 = buffer.data(gl + 146 * ncomps + c);
        const auto *gl_147 = buffer.data(gl + 147 * ncomps + c);
        const auto *gl_148 = buffer.data(gl + 148 * ncomps + c);
        const auto *gl_149 = buffer.data(gl + 149 * ncomps + c);
        const auto *gl_150 = buffer.data(gl + 150 * ncomps + c);
        const auto *gl_151 = buffer.data(gl + 151 * ncomps + c);
        const auto *gl_152 = buffer.data(gl + 152 * ncomps + c);
        const auto *gl_153 = buffer.data(gl + 153 * ncomps + c);
        const auto *gl_154 = buffer.data(gl + 154 * ncomps + c);
        const auto *gl_155 = buffer.data(gl + 155 * ncomps + c);
        const auto *gl_156 = buffer.data(gl + 156 * ncomps + c);
        const auto *gl_157 = buffer.data(gl + 157 * ncomps + c);
        const auto *gl_158 = buffer.data(gl + 158 * ncomps + c);
        const auto *gl_159 = buffer.data(gl + 159 * ncomps + c);
        const auto *gl_160 = buffer.data(gl + 160 * ncomps + c);
        const auto *gl_161 = buffer.data(gl + 161 * ncomps + c);
        const auto *gl_162 = buffer.data(gl + 162 * ncomps + c);
        const auto *gl_163 = buffer.data(gl + 163 * ncomps + c);
        const auto *gl_164 = buffer.data(gl + 164 * ncomps + c);
        const auto *gl_165 = buffer.data(gl + 165 * ncomps + c);
        const auto *gl_166 = buffer.data(gl + 166 * ncomps + c);
        const auto *gl_167 = buffer.data(gl + 167 * ncomps + c);
        const auto *gl_168 = buffer.data(gl + 168 * ncomps + c);
        const auto *gl_169 = buffer.data(gl + 169 * ncomps + c);
        const auto *gl_170 = buffer.data(gl + 170 * ncomps + c);
        const auto *gl_180 = buffer.data(gl + 180 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gk_0, gk_1, gk_2, gk_3, gk_4, gl_0, \
                         gl_1, gl_2, gl_3, gl_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * gk_0[k]
                     + gl_0[k];

            t_1[k] = -ab_x[k] * gk_1[k]
                     + gl_1[k];

            t_2[k] = -ab_x[k] * gk_2[k]
                     + gl_2[k];

            t_3[k] = -ab_x[k] * gk_3[k]
                     + gl_3[k];

            t_4[k] = -ab_x[k] * gk_4[k]
                     + gl_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gk_5, gk_6, gk_7, gk_8, gk_9, gl_5, \
                         gl_6, gl_7, gl_8, gl_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * gk_5[k]
                     + gl_5[k];

            t_6[k] = -ab_x[k] * gk_6[k]
                     + gl_6[k];

            t_7[k] = -ab_x[k] * gk_7[k]
                     + gl_7[k];

            t_8[k] = -ab_x[k] * gk_8[k]
                     + gl_8[k];

            t_9[k] = -ab_x[k] * gk_9[k]
                     + gl_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, gk_10, gk_11, gk_12, gk_13, \
                         gk_14, gl_10, gl_11, gl_12, gl_13, gl_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * gk_10[k]
                      + gl_10[k];

            t_11[k] = -ab_x[k] * gk_11[k]
                      + gl_11[k];

            t_12[k] = -ab_x[k] * gk_12[k]
                      + gl_12[k];

            t_13[k] = -ab_x[k] * gk_13[k]
                      + gl_13[k];

            t_14[k] = -ab_x[k] * gk_14[k]
                      + gl_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, gk_15, gk_16, gk_17, gk_18, \
                         gk_19, gl_15, gl_16, gl_17, gl_18, gl_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * gk_15[k]
                      + gl_15[k];

            t_16[k] = -ab_x[k] * gk_16[k]
                      + gl_16[k];

            t_17[k] = -ab_x[k] * gk_17[k]
                      + gl_17[k];

            t_18[k] = -ab_x[k] * gk_18[k]
                      + gl_18[k];

            t_19[k] = -ab_x[k] * gk_19[k]
                      + gl_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gk_20, gk_21, gk_22, gk_23, \
                         gk_24, gl_20, gl_21, gl_22, gl_23, gl_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * gk_20[k]
                      + gl_20[k];

            t_21[k] = -ab_x[k] * gk_21[k]
                      + gl_21[k];

            t_22[k] = -ab_x[k] * gk_22[k]
                      + gl_22[k];

            t_23[k] = -ab_x[k] * gk_23[k]
                      + gl_23[k];

            t_24[k] = -ab_x[k] * gk_24[k]
                      + gl_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, gk_25, gk_26, gk_27, gk_28, \
                         gk_29, gl_25, gl_26, gl_27, gl_28, gl_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * gk_25[k]
                      + gl_25[k];

            t_26[k] = -ab_x[k] * gk_26[k]
                      + gl_26[k];

            t_27[k] = -ab_x[k] * gk_27[k]
                      + gl_27[k];

            t_28[k] = -ab_x[k] * gk_28[k]
                      + gl_28[k];

            t_29[k] = -ab_x[k] * gk_29[k]
                      + gl_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gk_30, gk_31, gk_32, gk_33, \
                         gk_34, gl_30, gl_31, gl_32, gl_33, gl_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * gk_30[k]
                      + gl_30[k];

            t_31[k] = -ab_x[k] * gk_31[k]
                      + gl_31[k];

            t_32[k] = -ab_x[k] * gk_32[k]
                      + gl_32[k];

            t_33[k] = -ab_x[k] * gk_33[k]
                      + gl_33[k];

            t_34[k] = -ab_x[k] * gk_34[k]
                      + gl_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gk_35, gk_36, gk_37, gk_38, \
                         gk_39, gl_35, gl_45, gl_46, gl_47, gl_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * gk_35[k]
                      + gl_35[k];

            t_36[k] = -ab_x[k] * gk_36[k]
                      + gl_45[k];

            t_37[k] = -ab_x[k] * gk_37[k]
                      + gl_46[k];

            t_38[k] = -ab_x[k] * gk_38[k]
                      + gl_47[k];

            t_39[k] = -ab_x[k] * gk_39[k]
                      + gl_48[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gk_40, gk_41, gk_42, gk_43, \
                         gk_44, gl_49, gl_50, gl_51, gl_52, gl_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * gk_40[k]
                      + gl_49[k];

            t_41[k] = -ab_x[k] * gk_41[k]
                      + gl_50[k];

            t_42[k] = -ab_x[k] * gk_42[k]
                      + gl_51[k];

            t_43[k] = -ab_x[k] * gk_43[k]
                      + gl_52[k];

            t_44[k] = -ab_x[k] * gk_44[k]
                      + gl_53[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, gk_45, gk_46, gk_47, gk_48, \
                         gk_49, gl_54, gl_55, gl_56, gl_57, gl_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * gk_45[k]
                      + gl_54[k];

            t_46[k] = -ab_x[k] * gk_46[k]
                      + gl_55[k];

            t_47[k] = -ab_x[k] * gk_47[k]
                      + gl_56[k];

            t_48[k] = -ab_x[k] * gk_48[k]
                      + gl_57[k];

            t_49[k] = -ab_x[k] * gk_49[k]
                      + gl_58[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gk_50, gk_51, gk_52, gk_53, \
                         gk_54, gl_59, gl_60, gl_61, gl_62, gl_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * gk_50[k]
                      + gl_59[k];

            t_51[k] = -ab_x[k] * gk_51[k]
                      + gl_60[k];

            t_52[k] = -ab_x[k] * gk_52[k]
                      + gl_61[k];

            t_53[k] = -ab_x[k] * gk_53[k]
                      + gl_62[k];

            t_54[k] = -ab_x[k] * gk_54[k]
                      + gl_63[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, gk_55, gk_56, gk_57, gk_58, \
                         gk_59, gl_64, gl_65, gl_66, gl_67, gl_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * gk_55[k]
                      + gl_64[k];

            t_56[k] = -ab_x[k] * gk_56[k]
                      + gl_65[k];

            t_57[k] = -ab_x[k] * gk_57[k]
                      + gl_66[k];

            t_58[k] = -ab_x[k] * gk_58[k]
                      + gl_67[k];

            t_59[k] = -ab_x[k] * gk_59[k]
                      + gl_68[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gk_60, gk_61, gk_62, gk_63, \
                         gk_64, gl_69, gl_70, gl_71, gl_72, gl_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * gk_60[k]
                      + gl_69[k];

            t_61[k] = -ab_x[k] * gk_61[k]
                      + gl_70[k];

            t_62[k] = -ab_x[k] * gk_62[k]
                      + gl_71[k];

            t_63[k] = -ab_x[k] * gk_63[k]
                      + gl_72[k];

            t_64[k] = -ab_x[k] * gk_64[k]
                      + gl_73[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gk_65, gk_66, gk_67, gk_68, \
                         gk_69, gl_74, gl_75, gl_76, gl_77, gl_78 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * gk_65[k]
                      + gl_74[k];

            t_66[k] = -ab_x[k] * gk_66[k]
                      + gl_75[k];

            t_67[k] = -ab_x[k] * gk_67[k]
                      + gl_76[k];

            t_68[k] = -ab_x[k] * gk_68[k]
                      + gl_77[k];

            t_69[k] = -ab_x[k] * gk_69[k]
                      + gl_78[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, gk_70, gk_71, gk_72, gk_73, \
                         gk_74, gl_79, gl_80, gl_90, gl_91, gl_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * gk_70[k]
                      + gl_79[k];

            t_71[k] = -ab_x[k] * gk_71[k]
                      + gl_80[k];

            t_72[k] = -ab_x[k] * gk_72[k]
                      + gl_90[k];

            t_73[k] = -ab_x[k] * gk_73[k]
                      + gl_91[k];

            t_74[k] = -ab_x[k] * gk_74[k]
                      + gl_92[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, gk_75, gk_76, gk_77, gk_78, \
                         gk_79, gl_93, gl_94, gl_95, gl_96, gl_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * gk_75[k]
                      + gl_93[k];

            t_76[k] = -ab_x[k] * gk_76[k]
                      + gl_94[k];

            t_77[k] = -ab_x[k] * gk_77[k]
                      + gl_95[k];

            t_78[k] = -ab_x[k] * gk_78[k]
                      + gl_96[k];

            t_79[k] = -ab_x[k] * gk_79[k]
                      + gl_97[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gk_80, gk_81, gk_82, gk_83, \
                         gk_84, gl_98, gl_99, gl_100, gl_101, gl_102 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * gk_80[k]
                      + gl_98[k];

            t_81[k] = -ab_x[k] * gk_81[k]
                      + gl_99[k];

            t_82[k] = -ab_x[k] * gk_82[k]
                      + gl_100[k];

            t_83[k] = -ab_x[k] * gk_83[k]
                      + gl_101[k];

            t_84[k] = -ab_x[k] * gk_84[k]
                      + gl_102[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, gk_85, gk_86, gk_87, gk_88, \
                         gk_89, gl_103, gl_104, gl_105, gl_106, \
                         gl_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * gk_85[k]
                      + gl_103[k];

            t_86[k] = -ab_x[k] * gk_86[k]
                      + gl_104[k];

            t_87[k] = -ab_x[k] * gk_87[k]
                      + gl_105[k];

            t_88[k] = -ab_x[k] * gk_88[k]
                      + gl_106[k];

            t_89[k] = -ab_x[k] * gk_89[k]
                      + gl_107[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gk_90, gk_91, gk_92, gk_93, \
                         gk_94, gl_108, gl_109, gl_110, gl_111, \
                         gl_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * gk_90[k]
                      + gl_108[k];

            t_91[k] = -ab_x[k] * gk_91[k]
                      + gl_109[k];

            t_92[k] = -ab_x[k] * gk_92[k]
                      + gl_110[k];

            t_93[k] = -ab_x[k] * gk_93[k]
                      + gl_111[k];

            t_94[k] = -ab_x[k] * gk_94[k]
                      + gl_112[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gk_95, gk_96, gk_97, gk_98, \
                         gk_99, gl_113, gl_114, gl_115, gl_116, \
                         gl_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * gk_95[k]
                      + gl_113[k];

            t_96[k] = -ab_x[k] * gk_96[k]
                      + gl_114[k];

            t_97[k] = -ab_x[k] * gk_97[k]
                      + gl_115[k];

            t_98[k] = -ab_x[k] * gk_98[k]
                      + gl_116[k];

            t_99[k] = -ab_x[k] * gk_99[k]
                      + gl_117[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, gk_100, gk_101, gk_102, \
                         gk_103, gk_104, gl_118, gl_119, gl_120, gl_121, \
                         gl_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * gk_100[k]
                       + gl_118[k];

            t_101[k] = -ab_x[k] * gk_101[k]
                       + gl_119[k];

            t_102[k] = -ab_x[k] * gk_102[k]
                       + gl_120[k];

            t_103[k] = -ab_x[k] * gk_103[k]
                       + gl_121[k];

            t_104[k] = -ab_x[k] * gk_104[k]
                       + gl_122[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, gk_105, gk_106, gk_107, \
                         gk_108, gk_109, gl_123, gl_124, gl_125, gl_135, \
                         gl_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * gk_105[k]
                       + gl_123[k];

            t_106[k] = -ab_x[k] * gk_106[k]
                       + gl_124[k];

            t_107[k] = -ab_x[k] * gk_107[k]
                       + gl_125[k];

            t_108[k] = -ab_x[k] * gk_108[k]
                       + gl_135[k];

            t_109[k] = -ab_x[k] * gk_109[k]
                       + gl_136[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gk_110, gk_111, gk_112, \
                         gk_113, gk_114, gl_137, gl_138, gl_139, gl_140, \
                         gl_141 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * gk_110[k]
                       + gl_137[k];

            t_111[k] = -ab_x[k] * gk_111[k]
                       + gl_138[k];

            t_112[k] = -ab_x[k] * gk_112[k]
                       + gl_139[k];

            t_113[k] = -ab_x[k] * gk_113[k]
                       + gl_140[k];

            t_114[k] = -ab_x[k] * gk_114[k]
                       + gl_141[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, gk_115, gk_116, gk_117, \
                         gk_118, gk_119, gl_142, gl_143, gl_144, gl_145, \
                         gl_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * gk_115[k]
                       + gl_142[k];

            t_116[k] = -ab_x[k] * gk_116[k]
                       + gl_143[k];

            t_117[k] = -ab_x[k] * gk_117[k]
                       + gl_144[k];

            t_118[k] = -ab_x[k] * gk_118[k]
                       + gl_145[k];

            t_119[k] = -ab_x[k] * gk_119[k]
                       + gl_146[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gk_120, gk_121, gk_122, \
                         gk_123, gk_124, gl_147, gl_148, gl_149, gl_150, \
                         gl_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * gk_120[k]
                       + gl_147[k];

            t_121[k] = -ab_x[k] * gk_121[k]
                       + gl_148[k];

            t_122[k] = -ab_x[k] * gk_122[k]
                       + gl_149[k];

            t_123[k] = -ab_x[k] * gk_123[k]
                       + gl_150[k];

            t_124[k] = -ab_x[k] * gk_124[k]
                       + gl_151[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gk_125, gk_126, gk_127, \
                         gk_128, gk_129, gl_152, gl_153, gl_154, gl_155, \
                         gl_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * gk_125[k]
                       + gl_152[k];

            t_126[k] = -ab_x[k] * gk_126[k]
                       + gl_153[k];

            t_127[k] = -ab_x[k] * gk_127[k]
                       + gl_154[k];

            t_128[k] = -ab_x[k] * gk_128[k]
                       + gl_155[k];

            t_129[k] = -ab_x[k] * gk_129[k]
                       + gl_156[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, gk_130, gk_131, gk_132, \
                         gk_133, gk_134, gl_157, gl_158, gl_159, gl_160, \
                         gl_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * gk_130[k]
                       + gl_157[k];

            t_131[k] = -ab_x[k] * gk_131[k]
                       + gl_158[k];

            t_132[k] = -ab_x[k] * gk_132[k]
                       + gl_159[k];

            t_133[k] = -ab_x[k] * gk_133[k]
                       + gl_160[k];

            t_134[k] = -ab_x[k] * gk_134[k]
                       + gl_161[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, gk_135, gk_136, gk_137, \
                         gk_138, gk_139, gl_162, gl_163, gl_164, gl_165, \
                         gl_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * gk_135[k]
                       + gl_162[k];

            t_136[k] = -ab_x[k] * gk_136[k]
                       + gl_163[k];

            t_137[k] = -ab_x[k] * gk_137[k]
                       + gl_164[k];

            t_138[k] = -ab_x[k] * gk_138[k]
                       + gl_165[k];

            t_139[k] = -ab_x[k] * gk_139[k]
                       + gl_166[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gk_140, gk_141, gk_142, \
                         gk_143, gk_144, gl_167, gl_168, gl_169, gl_170, \
                         gl_180 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * gk_140[k]
                       + gl_167[k];

            t_141[k] = -ab_x[k] * gk_141[k]
                       + gl_168[k];

            t_142[k] = -ab_x[k] * gk_142[k]
                       + gl_169[k];

            t_143[k] = -ab_x[k] * gk_143[k]
                       + gl_170[k];

            t_144[k] = -ab_x[k] * gk_144[k]
                       + gl_180[k];
        }
    }
}

static auto
compute_hrr_hk_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gk, const size_t gl, const size_t ncomps,
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
        const auto *gk_172 = buffer.data(gk + 172 * ncomps + c);
        const auto *gk_173 = buffer.data(gk + 173 * ncomps + c);
        const auto *gk_174 = buffer.data(gk + 174 * ncomps + c);
        const auto *gk_175 = buffer.data(gk + 175 * ncomps + c);
        const auto *gk_176 = buffer.data(gk + 176 * ncomps + c);
        const auto *gk_177 = buffer.data(gk + 177 * ncomps + c);
        const auto *gk_178 = buffer.data(gk + 178 * ncomps + c);
        const auto *gk_179 = buffer.data(gk + 179 * ncomps + c);
        const auto *gk_180 = buffer.data(gk + 180 * ncomps + c);
        const auto *gk_181 = buffer.data(gk + 181 * ncomps + c);
        const auto *gk_182 = buffer.data(gk + 182 * ncomps + c);
        const auto *gk_183 = buffer.data(gk + 183 * ncomps + c);
        const auto *gk_184 = buffer.data(gk + 184 * ncomps + c);
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
        const auto *gk_208 = buffer.data(gk + 208 * ncomps + c);
        const auto *gk_209 = buffer.data(gk + 209 * ncomps + c);
        const auto *gk_210 = buffer.data(gk + 210 * ncomps + c);
        const auto *gk_211 = buffer.data(gk + 211 * ncomps + c);
        const auto *gk_212 = buffer.data(gk + 212 * ncomps + c);
        const auto *gk_213 = buffer.data(gk + 213 * ncomps + c);
        const auto *gk_214 = buffer.data(gk + 214 * ncomps + c);
        const auto *gk_215 = buffer.data(gk + 215 * ncomps + c);
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
        const auto *gk_244 = buffer.data(gk + 244 * ncomps + c);
        const auto *gk_245 = buffer.data(gk + 245 * ncomps + c);
        const auto *gk_246 = buffer.data(gk + 246 * ncomps + c);
        const auto *gk_247 = buffer.data(gk + 247 * ncomps + c);
        const auto *gk_248 = buffer.data(gk + 248 * ncomps + c);
        const auto *gk_249 = buffer.data(gk + 249 * ncomps + c);
        const auto *gk_250 = buffer.data(gk + 250 * ncomps + c);
        const auto *gk_251 = buffer.data(gk + 251 * ncomps + c);
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
        const auto *gk_280 = buffer.data(gk + 280 * ncomps + c);
        const auto *gk_281 = buffer.data(gk + 281 * ncomps + c);
        const auto *gk_282 = buffer.data(gk + 282 * ncomps + c);
        const auto *gk_283 = buffer.data(gk + 283 * ncomps + c);
        const auto *gk_284 = buffer.data(gk + 284 * ncomps + c);
        const auto *gk_285 = buffer.data(gk + 285 * ncomps + c);
        const auto *gk_286 = buffer.data(gk + 286 * ncomps + c);
        const auto *gk_287 = buffer.data(gk + 287 * ncomps + c);
        const auto *gk_288 = buffer.data(gk + 288 * ncomps + c);
        const auto *gk_289 = buffer.data(gk + 289 * ncomps + c);

        const auto *gl_181 = buffer.data(gl + 181 * ncomps + c);
        const auto *gl_182 = buffer.data(gl + 182 * ncomps + c);
        const auto *gl_183 = buffer.data(gl + 183 * ncomps + c);
        const auto *gl_184 = buffer.data(gl + 184 * ncomps + c);
        const auto *gl_185 = buffer.data(gl + 185 * ncomps + c);
        const auto *gl_186 = buffer.data(gl + 186 * ncomps + c);
        const auto *gl_187 = buffer.data(gl + 187 * ncomps + c);
        const auto *gl_188 = buffer.data(gl + 188 * ncomps + c);
        const auto *gl_189 = buffer.data(gl + 189 * ncomps + c);
        const auto *gl_190 = buffer.data(gl + 190 * ncomps + c);
        const auto *gl_191 = buffer.data(gl + 191 * ncomps + c);
        const auto *gl_192 = buffer.data(gl + 192 * ncomps + c);
        const auto *gl_193 = buffer.data(gl + 193 * ncomps + c);
        const auto *gl_194 = buffer.data(gl + 194 * ncomps + c);
        const auto *gl_195 = buffer.data(gl + 195 * ncomps + c);
        const auto *gl_196 = buffer.data(gl + 196 * ncomps + c);
        const auto *gl_197 = buffer.data(gl + 197 * ncomps + c);
        const auto *gl_198 = buffer.data(gl + 198 * ncomps + c);
        const auto *gl_199 = buffer.data(gl + 199 * ncomps + c);
        const auto *gl_200 = buffer.data(gl + 200 * ncomps + c);
        const auto *gl_201 = buffer.data(gl + 201 * ncomps + c);
        const auto *gl_202 = buffer.data(gl + 202 * ncomps + c);
        const auto *gl_203 = buffer.data(gl + 203 * ncomps + c);
        const auto *gl_204 = buffer.data(gl + 204 * ncomps + c);
        const auto *gl_205 = buffer.data(gl + 205 * ncomps + c);
        const auto *gl_206 = buffer.data(gl + 206 * ncomps + c);
        const auto *gl_207 = buffer.data(gl + 207 * ncomps + c);
        const auto *gl_208 = buffer.data(gl + 208 * ncomps + c);
        const auto *gl_209 = buffer.data(gl + 209 * ncomps + c);
        const auto *gl_210 = buffer.data(gl + 210 * ncomps + c);
        const auto *gl_211 = buffer.data(gl + 211 * ncomps + c);
        const auto *gl_212 = buffer.data(gl + 212 * ncomps + c);
        const auto *gl_213 = buffer.data(gl + 213 * ncomps + c);
        const auto *gl_214 = buffer.data(gl + 214 * ncomps + c);
        const auto *gl_215 = buffer.data(gl + 215 * ncomps + c);
        const auto *gl_225 = buffer.data(gl + 225 * ncomps + c);
        const auto *gl_226 = buffer.data(gl + 226 * ncomps + c);
        const auto *gl_227 = buffer.data(gl + 227 * ncomps + c);
        const auto *gl_228 = buffer.data(gl + 228 * ncomps + c);
        const auto *gl_229 = buffer.data(gl + 229 * ncomps + c);
        const auto *gl_230 = buffer.data(gl + 230 * ncomps + c);
        const auto *gl_231 = buffer.data(gl + 231 * ncomps + c);
        const auto *gl_232 = buffer.data(gl + 232 * ncomps + c);
        const auto *gl_233 = buffer.data(gl + 233 * ncomps + c);
        const auto *gl_234 = buffer.data(gl + 234 * ncomps + c);
        const auto *gl_235 = buffer.data(gl + 235 * ncomps + c);
        const auto *gl_236 = buffer.data(gl + 236 * ncomps + c);
        const auto *gl_237 = buffer.data(gl + 237 * ncomps + c);
        const auto *gl_238 = buffer.data(gl + 238 * ncomps + c);
        const auto *gl_239 = buffer.data(gl + 239 * ncomps + c);
        const auto *gl_240 = buffer.data(gl + 240 * ncomps + c);
        const auto *gl_241 = buffer.data(gl + 241 * ncomps + c);
        const auto *gl_242 = buffer.data(gl + 242 * ncomps + c);
        const auto *gl_243 = buffer.data(gl + 243 * ncomps + c);
        const auto *gl_244 = buffer.data(gl + 244 * ncomps + c);
        const auto *gl_245 = buffer.data(gl + 245 * ncomps + c);
        const auto *gl_246 = buffer.data(gl + 246 * ncomps + c);
        const auto *gl_247 = buffer.data(gl + 247 * ncomps + c);
        const auto *gl_248 = buffer.data(gl + 248 * ncomps + c);
        const auto *gl_249 = buffer.data(gl + 249 * ncomps + c);
        const auto *gl_250 = buffer.data(gl + 250 * ncomps + c);
        const auto *gl_251 = buffer.data(gl + 251 * ncomps + c);
        const auto *gl_252 = buffer.data(gl + 252 * ncomps + c);
        const auto *gl_253 = buffer.data(gl + 253 * ncomps + c);
        const auto *gl_254 = buffer.data(gl + 254 * ncomps + c);
        const auto *gl_255 = buffer.data(gl + 255 * ncomps + c);
        const auto *gl_256 = buffer.data(gl + 256 * ncomps + c);
        const auto *gl_257 = buffer.data(gl + 257 * ncomps + c);
        const auto *gl_258 = buffer.data(gl + 258 * ncomps + c);
        const auto *gl_259 = buffer.data(gl + 259 * ncomps + c);
        const auto *gl_260 = buffer.data(gl + 260 * ncomps + c);
        const auto *gl_270 = buffer.data(gl + 270 * ncomps + c);
        const auto *gl_271 = buffer.data(gl + 271 * ncomps + c);
        const auto *gl_272 = buffer.data(gl + 272 * ncomps + c);
        const auto *gl_273 = buffer.data(gl + 273 * ncomps + c);
        const auto *gl_274 = buffer.data(gl + 274 * ncomps + c);
        const auto *gl_275 = buffer.data(gl + 275 * ncomps + c);
        const auto *gl_276 = buffer.data(gl + 276 * ncomps + c);
        const auto *gl_277 = buffer.data(gl + 277 * ncomps + c);
        const auto *gl_278 = buffer.data(gl + 278 * ncomps + c);
        const auto *gl_279 = buffer.data(gl + 279 * ncomps + c);
        const auto *gl_280 = buffer.data(gl + 280 * ncomps + c);
        const auto *gl_281 = buffer.data(gl + 281 * ncomps + c);
        const auto *gl_282 = buffer.data(gl + 282 * ncomps + c);
        const auto *gl_283 = buffer.data(gl + 283 * ncomps + c);
        const auto *gl_284 = buffer.data(gl + 284 * ncomps + c);
        const auto *gl_285 = buffer.data(gl + 285 * ncomps + c);
        const auto *gl_286 = buffer.data(gl + 286 * ncomps + c);
        const auto *gl_287 = buffer.data(gl + 287 * ncomps + c);
        const auto *gl_288 = buffer.data(gl + 288 * ncomps + c);
        const auto *gl_289 = buffer.data(gl + 289 * ncomps + c);
        const auto *gl_290 = buffer.data(gl + 290 * ncomps + c);
        const auto *gl_291 = buffer.data(gl + 291 * ncomps + c);
        const auto *gl_292 = buffer.data(gl + 292 * ncomps + c);
        const auto *gl_293 = buffer.data(gl + 293 * ncomps + c);
        const auto *gl_294 = buffer.data(gl + 294 * ncomps + c);
        const auto *gl_295 = buffer.data(gl + 295 * ncomps + c);
        const auto *gl_296 = buffer.data(gl + 296 * ncomps + c);
        const auto *gl_297 = buffer.data(gl + 297 * ncomps + c);
        const auto *gl_298 = buffer.data(gl + 298 * ncomps + c);
        const auto *gl_299 = buffer.data(gl + 299 * ncomps + c);
        const auto *gl_300 = buffer.data(gl + 300 * ncomps + c);
        const auto *gl_301 = buffer.data(gl + 301 * ncomps + c);
        const auto *gl_302 = buffer.data(gl + 302 * ncomps + c);
        const auto *gl_303 = buffer.data(gl + 303 * ncomps + c);
        const auto *gl_304 = buffer.data(gl + 304 * ncomps + c);
        const auto *gl_305 = buffer.data(gl + 305 * ncomps + c);
        const auto *gl_315 = buffer.data(gl + 315 * ncomps + c);
        const auto *gl_316 = buffer.data(gl + 316 * ncomps + c);
        const auto *gl_317 = buffer.data(gl + 317 * ncomps + c);
        const auto *gl_318 = buffer.data(gl + 318 * ncomps + c);
        const auto *gl_319 = buffer.data(gl + 319 * ncomps + c);
        const auto *gl_320 = buffer.data(gl + 320 * ncomps + c);
        const auto *gl_321 = buffer.data(gl + 321 * ncomps + c);
        const auto *gl_322 = buffer.data(gl + 322 * ncomps + c);
        const auto *gl_323 = buffer.data(gl + 323 * ncomps + c);
        const auto *gl_324 = buffer.data(gl + 324 * ncomps + c);
        const auto *gl_325 = buffer.data(gl + 325 * ncomps + c);
        const auto *gl_326 = buffer.data(gl + 326 * ncomps + c);
        const auto *gl_327 = buffer.data(gl + 327 * ncomps + c);
        const auto *gl_328 = buffer.data(gl + 328 * ncomps + c);
        const auto *gl_329 = buffer.data(gl + 329 * ncomps + c);
        const auto *gl_330 = buffer.data(gl + 330 * ncomps + c);
        const auto *gl_331 = buffer.data(gl + 331 * ncomps + c);
        const auto *gl_332 = buffer.data(gl + 332 * ncomps + c);
        const auto *gl_333 = buffer.data(gl + 333 * ncomps + c);
        const auto *gl_334 = buffer.data(gl + 334 * ncomps + c);
        const auto *gl_335 = buffer.data(gl + 335 * ncomps + c);
        const auto *gl_336 = buffer.data(gl + 336 * ncomps + c);
        const auto *gl_337 = buffer.data(gl + 337 * ncomps + c);
        const auto *gl_338 = buffer.data(gl + 338 * ncomps + c);
        const auto *gl_339 = buffer.data(gl + 339 * ncomps + c);
        const auto *gl_340 = buffer.data(gl + 340 * ncomps + c);
        const auto *gl_341 = buffer.data(gl + 341 * ncomps + c);
        const auto *gl_342 = buffer.data(gl + 342 * ncomps + c);
        const auto *gl_343 = buffer.data(gl + 343 * ncomps + c);
        const auto *gl_344 = buffer.data(gl + 344 * ncomps + c);
        const auto *gl_345 = buffer.data(gl + 345 * ncomps + c);
        const auto *gl_346 = buffer.data(gl + 346 * ncomps + c);
        const auto *gl_347 = buffer.data(gl + 347 * ncomps + c);
        const auto *gl_348 = buffer.data(gl + 348 * ncomps + c);
        const auto *gl_349 = buffer.data(gl + 349 * ncomps + c);
        const auto *gl_350 = buffer.data(gl + 350 * ncomps + c);
        const auto *gl_360 = buffer.data(gl + 360 * ncomps + c);
        const auto *gl_361 = buffer.data(gl + 361 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, gk_145, gk_146, gk_147, \
                         gk_148, gk_149, gl_181, gl_182, gl_183, gl_184, \
                         gl_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * gk_145[k]
                       + gl_181[k];

            t_146[k] = -ab_x[k] * gk_146[k]
                       + gl_182[k];

            t_147[k] = -ab_x[k] * gk_147[k]
                       + gl_183[k];

            t_148[k] = -ab_x[k] * gk_148[k]
                       + gl_184[k];

            t_149[k] = -ab_x[k] * gk_149[k]
                       + gl_185[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gk_150, gk_151, gk_152, \
                         gk_153, gk_154, gl_186, gl_187, gl_188, gl_189, \
                         gl_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * gk_150[k]
                       + gl_186[k];

            t_151[k] = -ab_x[k] * gk_151[k]
                       + gl_187[k];

            t_152[k] = -ab_x[k] * gk_152[k]
                       + gl_188[k];

            t_153[k] = -ab_x[k] * gk_153[k]
                       + gl_189[k];

            t_154[k] = -ab_x[k] * gk_154[k]
                       + gl_190[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gk_155, gk_156, gk_157, \
                         gk_158, gk_159, gl_191, gl_192, gl_193, gl_194, \
                         gl_195 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * gk_155[k]
                       + gl_191[k];

            t_156[k] = -ab_x[k] * gk_156[k]
                       + gl_192[k];

            t_157[k] = -ab_x[k] * gk_157[k]
                       + gl_193[k];

            t_158[k] = -ab_x[k] * gk_158[k]
                       + gl_194[k];

            t_159[k] = -ab_x[k] * gk_159[k]
                       + gl_195[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, gk_160, gk_161, gk_162, \
                         gk_163, gk_164, gl_196, gl_197, gl_198, gl_199, \
                         gl_200 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * gk_160[k]
                       + gl_196[k];

            t_161[k] = -ab_x[k] * gk_161[k]
                       + gl_197[k];

            t_162[k] = -ab_x[k] * gk_162[k]
                       + gl_198[k];

            t_163[k] = -ab_x[k] * gk_163[k]
                       + gl_199[k];

            t_164[k] = -ab_x[k] * gk_164[k]
                       + gl_200[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, gk_165, gk_166, gk_167, \
                         gk_168, gk_169, gl_201, gl_202, gl_203, gl_204, \
                         gl_205 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * gk_165[k]
                       + gl_201[k];

            t_166[k] = -ab_x[k] * gk_166[k]
                       + gl_202[k];

            t_167[k] = -ab_x[k] * gk_167[k]
                       + gl_203[k];

            t_168[k] = -ab_x[k] * gk_168[k]
                       + gl_204[k];

            t_169[k] = -ab_x[k] * gk_169[k]
                       + gl_205[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, gk_170, gk_171, gk_172, \
                         gk_173, gk_174, gl_206, gl_207, gl_208, gl_209, \
                         gl_210 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * gk_170[k]
                       + gl_206[k];

            t_171[k] = -ab_x[k] * gk_171[k]
                       + gl_207[k];

            t_172[k] = -ab_x[k] * gk_172[k]
                       + gl_208[k];

            t_173[k] = -ab_x[k] * gk_173[k]
                       + gl_209[k];

            t_174[k] = -ab_x[k] * gk_174[k]
                       + gl_210[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, gk_175, gk_176, gk_177, \
                         gk_178, gk_179, gl_211, gl_212, gl_213, gl_214, \
                         gl_215 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * gk_175[k]
                       + gl_211[k];

            t_176[k] = -ab_x[k] * gk_176[k]
                       + gl_212[k];

            t_177[k] = -ab_x[k] * gk_177[k]
                       + gl_213[k];

            t_178[k] = -ab_x[k] * gk_178[k]
                       + gl_214[k];

            t_179[k] = -ab_x[k] * gk_179[k]
                       + gl_215[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gk_180, gk_181, gk_182, \
                         gk_183, gk_184, gl_225, gl_226, gl_227, gl_228, \
                         gl_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * gk_180[k]
                       + gl_225[k];

            t_181[k] = -ab_x[k] * gk_181[k]
                       + gl_226[k];

            t_182[k] = -ab_x[k] * gk_182[k]
                       + gl_227[k];

            t_183[k] = -ab_x[k] * gk_183[k]
                       + gl_228[k];

            t_184[k] = -ab_x[k] * gk_184[k]
                       + gl_229[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gk_185, gk_186, gk_187, \
                         gk_188, gk_189, gl_230, gl_231, gl_232, gl_233, \
                         gl_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * gk_185[k]
                       + gl_230[k];

            t_186[k] = -ab_x[k] * gk_186[k]
                       + gl_231[k];

            t_187[k] = -ab_x[k] * gk_187[k]
                       + gl_232[k];

            t_188[k] = -ab_x[k] * gk_188[k]
                       + gl_233[k];

            t_189[k] = -ab_x[k] * gk_189[k]
                       + gl_234[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, gk_190, gk_191, gk_192, \
                         gk_193, gk_194, gl_235, gl_236, gl_237, gl_238, \
                         gl_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * gk_190[k]
                       + gl_235[k];

            t_191[k] = -ab_x[k] * gk_191[k]
                       + gl_236[k];

            t_192[k] = -ab_x[k] * gk_192[k]
                       + gl_237[k];

            t_193[k] = -ab_x[k] * gk_193[k]
                       + gl_238[k];

            t_194[k] = -ab_x[k] * gk_194[k]
                       + gl_239[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, gk_195, gk_196, gk_197, \
                         gk_198, gk_199, gl_240, gl_241, gl_242, gl_243, \
                         gl_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * gk_195[k]
                       + gl_240[k];

            t_196[k] = -ab_x[k] * gk_196[k]
                       + gl_241[k];

            t_197[k] = -ab_x[k] * gk_197[k]
                       + gl_242[k];

            t_198[k] = -ab_x[k] * gk_198[k]
                       + gl_243[k];

            t_199[k] = -ab_x[k] * gk_199[k]
                       + gl_244[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, gk_200, gk_201, gk_202, \
                         gk_203, gk_204, gl_245, gl_246, gl_247, gl_248, \
                         gl_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * gk_200[k]
                       + gl_245[k];

            t_201[k] = -ab_x[k] * gk_201[k]
                       + gl_246[k];

            t_202[k] = -ab_x[k] * gk_202[k]
                       + gl_247[k];

            t_203[k] = -ab_x[k] * gk_203[k]
                       + gl_248[k];

            t_204[k] = -ab_x[k] * gk_204[k]
                       + gl_249[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, gk_205, gk_206, gk_207, \
                         gk_208, gk_209, gl_250, gl_251, gl_252, gl_253, \
                         gl_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * gk_205[k]
                       + gl_250[k];

            t_206[k] = -ab_x[k] * gk_206[k]
                       + gl_251[k];

            t_207[k] = -ab_x[k] * gk_207[k]
                       + gl_252[k];

            t_208[k] = -ab_x[k] * gk_208[k]
                       + gl_253[k];

            t_209[k] = -ab_x[k] * gk_209[k]
                       + gl_254[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gk_210, gk_211, gk_212, \
                         gk_213, gk_214, gl_255, gl_256, gl_257, gl_258, \
                         gl_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * gk_210[k]
                       + gl_255[k];

            t_211[k] = -ab_x[k] * gk_211[k]
                       + gl_256[k];

            t_212[k] = -ab_x[k] * gk_212[k]
                       + gl_257[k];

            t_213[k] = -ab_x[k] * gk_213[k]
                       + gl_258[k];

            t_214[k] = -ab_x[k] * gk_214[k]
                       + gl_259[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gk_215, gk_216, gk_217, \
                         gk_218, gk_219, gl_260, gl_270, gl_271, gl_272, \
                         gl_273 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * gk_215[k]
                       + gl_260[k];

            t_216[k] = -ab_x[k] * gk_216[k]
                       + gl_270[k];

            t_217[k] = -ab_x[k] * gk_217[k]
                       + gl_271[k];

            t_218[k] = -ab_x[k] * gk_218[k]
                       + gl_272[k];

            t_219[k] = -ab_x[k] * gk_219[k]
                       + gl_273[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, gk_220, gk_221, gk_222, \
                         gk_223, gk_224, gl_274, gl_275, gl_276, gl_277, \
                         gl_278 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_x[k] * gk_220[k]
                       + gl_274[k];

            t_221[k] = -ab_x[k] * gk_221[k]
                       + gl_275[k];

            t_222[k] = -ab_x[k] * gk_222[k]
                       + gl_276[k];

            t_223[k] = -ab_x[k] * gk_223[k]
                       + gl_277[k];

            t_224[k] = -ab_x[k] * gk_224[k]
                       + gl_278[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, gk_225, gk_226, gk_227, \
                         gk_228, gk_229, gl_279, gl_280, gl_281, gl_282, \
                         gl_283 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_x[k] * gk_225[k]
                       + gl_279[k];

            t_226[k] = -ab_x[k] * gk_226[k]
                       + gl_280[k];

            t_227[k] = -ab_x[k] * gk_227[k]
                       + gl_281[k];

            t_228[k] = -ab_x[k] * gk_228[k]
                       + gl_282[k];

            t_229[k] = -ab_x[k] * gk_229[k]
                       + gl_283[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, gk_230, gk_231, gk_232, \
                         gk_233, gk_234, gl_284, gl_285, gl_286, gl_287, \
                         gl_288 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_x[k] * gk_230[k]
                       + gl_284[k];

            t_231[k] = -ab_x[k] * gk_231[k]
                       + gl_285[k];

            t_232[k] = -ab_x[k] * gk_232[k]
                       + gl_286[k];

            t_233[k] = -ab_x[k] * gk_233[k]
                       + gl_287[k];

            t_234[k] = -ab_x[k] * gk_234[k]
                       + gl_288[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, gk_235, gk_236, gk_237, \
                         gk_238, gk_239, gl_289, gl_290, gl_291, gl_292, \
                         gl_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_x[k] * gk_235[k]
                       + gl_289[k];

            t_236[k] = -ab_x[k] * gk_236[k]
                       + gl_290[k];

            t_237[k] = -ab_x[k] * gk_237[k]
                       + gl_291[k];

            t_238[k] = -ab_x[k] * gk_238[k]
                       + gl_292[k];

            t_239[k] = -ab_x[k] * gk_239[k]
                       + gl_293[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, gk_240, gk_241, gk_242, \
                         gk_243, gk_244, gl_294, gl_295, gl_296, gl_297, \
                         gl_298 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_x[k] * gk_240[k]
                       + gl_294[k];

            t_241[k] = -ab_x[k] * gk_241[k]
                       + gl_295[k];

            t_242[k] = -ab_x[k] * gk_242[k]
                       + gl_296[k];

            t_243[k] = -ab_x[k] * gk_243[k]
                       + gl_297[k];

            t_244[k] = -ab_x[k] * gk_244[k]
                       + gl_298[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, gk_245, gk_246, gk_247, \
                         gk_248, gk_249, gl_299, gl_300, gl_301, gl_302, \
                         gl_303 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_x[k] * gk_245[k]
                       + gl_299[k];

            t_246[k] = -ab_x[k] * gk_246[k]
                       + gl_300[k];

            t_247[k] = -ab_x[k] * gk_247[k]
                       + gl_301[k];

            t_248[k] = -ab_x[k] * gk_248[k]
                       + gl_302[k];

            t_249[k] = -ab_x[k] * gk_249[k]
                       + gl_303[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, gk_250, gk_251, gk_252, \
                         gk_253, gk_254, gl_304, gl_305, gl_315, gl_316, \
                         gl_317 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_x[k] * gk_250[k]
                       + gl_304[k];

            t_251[k] = -ab_x[k] * gk_251[k]
                       + gl_305[k];

            t_252[k] = -ab_x[k] * gk_252[k]
                       + gl_315[k];

            t_253[k] = -ab_x[k] * gk_253[k]
                       + gl_316[k];

            t_254[k] = -ab_x[k] * gk_254[k]
                       + gl_317[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, gk_255, gk_256, gk_257, \
                         gk_258, gk_259, gl_318, gl_319, gl_320, gl_321, \
                         gl_322 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_x[k] * gk_255[k]
                       + gl_318[k];

            t_256[k] = -ab_x[k] * gk_256[k]
                       + gl_319[k];

            t_257[k] = -ab_x[k] * gk_257[k]
                       + gl_320[k];

            t_258[k] = -ab_x[k] * gk_258[k]
                       + gl_321[k];

            t_259[k] = -ab_x[k] * gk_259[k]
                       + gl_322[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, gk_260, gk_261, gk_262, \
                         gk_263, gk_264, gl_323, gl_324, gl_325, gl_326, \
                         gl_327 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_x[k] * gk_260[k]
                       + gl_323[k];

            t_261[k] = -ab_x[k] * gk_261[k]
                       + gl_324[k];

            t_262[k] = -ab_x[k] * gk_262[k]
                       + gl_325[k];

            t_263[k] = -ab_x[k] * gk_263[k]
                       + gl_326[k];

            t_264[k] = -ab_x[k] * gk_264[k]
                       + gl_327[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, gk_265, gk_266, gk_267, \
                         gk_268, gk_269, gl_328, gl_329, gl_330, gl_331, \
                         gl_332 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_x[k] * gk_265[k]
                       + gl_328[k];

            t_266[k] = -ab_x[k] * gk_266[k]
                       + gl_329[k];

            t_267[k] = -ab_x[k] * gk_267[k]
                       + gl_330[k];

            t_268[k] = -ab_x[k] * gk_268[k]
                       + gl_331[k];

            t_269[k] = -ab_x[k] * gk_269[k]
                       + gl_332[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, gk_270, gk_271, gk_272, \
                         gk_273, gk_274, gl_333, gl_334, gl_335, gl_336, \
                         gl_337 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_x[k] * gk_270[k]
                       + gl_333[k];

            t_271[k] = -ab_x[k] * gk_271[k]
                       + gl_334[k];

            t_272[k] = -ab_x[k] * gk_272[k]
                       + gl_335[k];

            t_273[k] = -ab_x[k] * gk_273[k]
                       + gl_336[k];

            t_274[k] = -ab_x[k] * gk_274[k]
                       + gl_337[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, gk_275, gk_276, gk_277, \
                         gk_278, gk_279, gl_338, gl_339, gl_340, gl_341, \
                         gl_342 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_x[k] * gk_275[k]
                       + gl_338[k];

            t_276[k] = -ab_x[k] * gk_276[k]
                       + gl_339[k];

            t_277[k] = -ab_x[k] * gk_277[k]
                       + gl_340[k];

            t_278[k] = -ab_x[k] * gk_278[k]
                       + gl_341[k];

            t_279[k] = -ab_x[k] * gk_279[k]
                       + gl_342[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, gk_280, gk_281, gk_282, \
                         gk_283, gk_284, gl_343, gl_344, gl_345, gl_346, \
                         gl_347 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_x[k] * gk_280[k]
                       + gl_343[k];

            t_281[k] = -ab_x[k] * gk_281[k]
                       + gl_344[k];

            t_282[k] = -ab_x[k] * gk_282[k]
                       + gl_345[k];

            t_283[k] = -ab_x[k] * gk_283[k]
                       + gl_346[k];

            t_284[k] = -ab_x[k] * gk_284[k]
                       + gl_347[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, gk_285, gk_286, gk_287, \
                         gk_288, gk_289, gl_348, gl_349, gl_350, gl_360, \
                         gl_361 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_x[k] * gk_285[k]
                       + gl_348[k];

            t_286[k] = -ab_x[k] * gk_286[k]
                       + gl_349[k];

            t_287[k] = -ab_x[k] * gk_287[k]
                       + gl_350[k];

            t_288[k] = -ab_x[k] * gk_288[k]
                       + gl_360[k];

            t_289[k] = -ab_x[k] * gk_289[k]
                       + gl_361[k];
        }
    }
}

static auto
compute_hrr_hk_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gk, const size_t gl, const size_t ncomps,
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
        const auto *gk_316 = buffer.data(gk + 316 * ncomps + c);
        const auto *gk_317 = buffer.data(gk + 317 * ncomps + c);
        const auto *gk_318 = buffer.data(gk + 318 * ncomps + c);
        const auto *gk_319 = buffer.data(gk + 319 * ncomps + c);
        const auto *gk_320 = buffer.data(gk + 320 * ncomps + c);
        const auto *gk_321 = buffer.data(gk + 321 * ncomps + c);
        const auto *gk_322 = buffer.data(gk + 322 * ncomps + c);
        const auto *gk_323 = buffer.data(gk + 323 * ncomps + c);
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
        const auto *gk_352 = buffer.data(gk + 352 * ncomps + c);
        const auto *gk_353 = buffer.data(gk + 353 * ncomps + c);
        const auto *gk_354 = buffer.data(gk + 354 * ncomps + c);
        const auto *gk_355 = buffer.data(gk + 355 * ncomps + c);
        const auto *gk_356 = buffer.data(gk + 356 * ncomps + c);
        const auto *gk_357 = buffer.data(gk + 357 * ncomps + c);
        const auto *gk_358 = buffer.data(gk + 358 * ncomps + c);
        const auto *gk_359 = buffer.data(gk + 359 * ncomps + c);
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
        const auto *gk_388 = buffer.data(gk + 388 * ncomps + c);
        const auto *gk_389 = buffer.data(gk + 389 * ncomps + c);
        const auto *gk_390 = buffer.data(gk + 390 * ncomps + c);
        const auto *gk_391 = buffer.data(gk + 391 * ncomps + c);
        const auto *gk_392 = buffer.data(gk + 392 * ncomps + c);
        const auto *gk_393 = buffer.data(gk + 393 * ncomps + c);
        const auto *gk_394 = buffer.data(gk + 394 * ncomps + c);
        const auto *gk_395 = buffer.data(gk + 395 * ncomps + c);
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
        const auto *gk_424 = buffer.data(gk + 424 * ncomps + c);
        const auto *gk_425 = buffer.data(gk + 425 * ncomps + c);
        const auto *gk_426 = buffer.data(gk + 426 * ncomps + c);
        const auto *gk_427 = buffer.data(gk + 427 * ncomps + c);
        const auto *gk_428 = buffer.data(gk + 428 * ncomps + c);
        const auto *gk_429 = buffer.data(gk + 429 * ncomps + c);
        const auto *gk_430 = buffer.data(gk + 430 * ncomps + c);
        const auto *gk_431 = buffer.data(gk + 431 * ncomps + c);
        const auto *gk_432 = buffer.data(gk + 432 * ncomps + c);
        const auto *gk_433 = buffer.data(gk + 433 * ncomps + c);
        const auto *gk_434 = buffer.data(gk + 434 * ncomps + c);

        const auto *gl_362 = buffer.data(gl + 362 * ncomps + c);
        const auto *gl_363 = buffer.data(gl + 363 * ncomps + c);
        const auto *gl_364 = buffer.data(gl + 364 * ncomps + c);
        const auto *gl_365 = buffer.data(gl + 365 * ncomps + c);
        const auto *gl_366 = buffer.data(gl + 366 * ncomps + c);
        const auto *gl_367 = buffer.data(gl + 367 * ncomps + c);
        const auto *gl_368 = buffer.data(gl + 368 * ncomps + c);
        const auto *gl_369 = buffer.data(gl + 369 * ncomps + c);
        const auto *gl_370 = buffer.data(gl + 370 * ncomps + c);
        const auto *gl_371 = buffer.data(gl + 371 * ncomps + c);
        const auto *gl_372 = buffer.data(gl + 372 * ncomps + c);
        const auto *gl_373 = buffer.data(gl + 373 * ncomps + c);
        const auto *gl_374 = buffer.data(gl + 374 * ncomps + c);
        const auto *gl_375 = buffer.data(gl + 375 * ncomps + c);
        const auto *gl_376 = buffer.data(gl + 376 * ncomps + c);
        const auto *gl_377 = buffer.data(gl + 377 * ncomps + c);
        const auto *gl_378 = buffer.data(gl + 378 * ncomps + c);
        const auto *gl_379 = buffer.data(gl + 379 * ncomps + c);
        const auto *gl_380 = buffer.data(gl + 380 * ncomps + c);
        const auto *gl_381 = buffer.data(gl + 381 * ncomps + c);
        const auto *gl_382 = buffer.data(gl + 382 * ncomps + c);
        const auto *gl_383 = buffer.data(gl + 383 * ncomps + c);
        const auto *gl_384 = buffer.data(gl + 384 * ncomps + c);
        const auto *gl_385 = buffer.data(gl + 385 * ncomps + c);
        const auto *gl_386 = buffer.data(gl + 386 * ncomps + c);
        const auto *gl_387 = buffer.data(gl + 387 * ncomps + c);
        const auto *gl_388 = buffer.data(gl + 388 * ncomps + c);
        const auto *gl_389 = buffer.data(gl + 389 * ncomps + c);
        const auto *gl_390 = buffer.data(gl + 390 * ncomps + c);
        const auto *gl_391 = buffer.data(gl + 391 * ncomps + c);
        const auto *gl_392 = buffer.data(gl + 392 * ncomps + c);
        const auto *gl_393 = buffer.data(gl + 393 * ncomps + c);
        const auto *gl_394 = buffer.data(gl + 394 * ncomps + c);
        const auto *gl_395 = buffer.data(gl + 395 * ncomps + c);
        const auto *gl_405 = buffer.data(gl + 405 * ncomps + c);
        const auto *gl_406 = buffer.data(gl + 406 * ncomps + c);
        const auto *gl_407 = buffer.data(gl + 407 * ncomps + c);
        const auto *gl_408 = buffer.data(gl + 408 * ncomps + c);
        const auto *gl_409 = buffer.data(gl + 409 * ncomps + c);
        const auto *gl_410 = buffer.data(gl + 410 * ncomps + c);
        const auto *gl_411 = buffer.data(gl + 411 * ncomps + c);
        const auto *gl_412 = buffer.data(gl + 412 * ncomps + c);
        const auto *gl_413 = buffer.data(gl + 413 * ncomps + c);
        const auto *gl_414 = buffer.data(gl + 414 * ncomps + c);
        const auto *gl_415 = buffer.data(gl + 415 * ncomps + c);
        const auto *gl_416 = buffer.data(gl + 416 * ncomps + c);
        const auto *gl_417 = buffer.data(gl + 417 * ncomps + c);
        const auto *gl_418 = buffer.data(gl + 418 * ncomps + c);
        const auto *gl_419 = buffer.data(gl + 419 * ncomps + c);
        const auto *gl_420 = buffer.data(gl + 420 * ncomps + c);
        const auto *gl_421 = buffer.data(gl + 421 * ncomps + c);
        const auto *gl_422 = buffer.data(gl + 422 * ncomps + c);
        const auto *gl_423 = buffer.data(gl + 423 * ncomps + c);
        const auto *gl_424 = buffer.data(gl + 424 * ncomps + c);
        const auto *gl_425 = buffer.data(gl + 425 * ncomps + c);
        const auto *gl_426 = buffer.data(gl + 426 * ncomps + c);
        const auto *gl_427 = buffer.data(gl + 427 * ncomps + c);
        const auto *gl_428 = buffer.data(gl + 428 * ncomps + c);
        const auto *gl_429 = buffer.data(gl + 429 * ncomps + c);
        const auto *gl_430 = buffer.data(gl + 430 * ncomps + c);
        const auto *gl_431 = buffer.data(gl + 431 * ncomps + c);
        const auto *gl_432 = buffer.data(gl + 432 * ncomps + c);
        const auto *gl_433 = buffer.data(gl + 433 * ncomps + c);
        const auto *gl_434 = buffer.data(gl + 434 * ncomps + c);
        const auto *gl_435 = buffer.data(gl + 435 * ncomps + c);
        const auto *gl_436 = buffer.data(gl + 436 * ncomps + c);
        const auto *gl_437 = buffer.data(gl + 437 * ncomps + c);
        const auto *gl_438 = buffer.data(gl + 438 * ncomps + c);
        const auto *gl_439 = buffer.data(gl + 439 * ncomps + c);
        const auto *gl_440 = buffer.data(gl + 440 * ncomps + c);
        const auto *gl_450 = buffer.data(gl + 450 * ncomps + c);
        const auto *gl_451 = buffer.data(gl + 451 * ncomps + c);
        const auto *gl_452 = buffer.data(gl + 452 * ncomps + c);
        const auto *gl_453 = buffer.data(gl + 453 * ncomps + c);
        const auto *gl_454 = buffer.data(gl + 454 * ncomps + c);
        const auto *gl_455 = buffer.data(gl + 455 * ncomps + c);
        const auto *gl_456 = buffer.data(gl + 456 * ncomps + c);
        const auto *gl_457 = buffer.data(gl + 457 * ncomps + c);
        const auto *gl_458 = buffer.data(gl + 458 * ncomps + c);
        const auto *gl_459 = buffer.data(gl + 459 * ncomps + c);
        const auto *gl_460 = buffer.data(gl + 460 * ncomps + c);
        const auto *gl_461 = buffer.data(gl + 461 * ncomps + c);
        const auto *gl_462 = buffer.data(gl + 462 * ncomps + c);
        const auto *gl_463 = buffer.data(gl + 463 * ncomps + c);
        const auto *gl_464 = buffer.data(gl + 464 * ncomps + c);
        const auto *gl_465 = buffer.data(gl + 465 * ncomps + c);
        const auto *gl_466 = buffer.data(gl + 466 * ncomps + c);
        const auto *gl_467 = buffer.data(gl + 467 * ncomps + c);
        const auto *gl_468 = buffer.data(gl + 468 * ncomps + c);
        const auto *gl_469 = buffer.data(gl + 469 * ncomps + c);
        const auto *gl_470 = buffer.data(gl + 470 * ncomps + c);
        const auto *gl_471 = buffer.data(gl + 471 * ncomps + c);
        const auto *gl_472 = buffer.data(gl + 472 * ncomps + c);
        const auto *gl_473 = buffer.data(gl + 473 * ncomps + c);
        const auto *gl_474 = buffer.data(gl + 474 * ncomps + c);
        const auto *gl_475 = buffer.data(gl + 475 * ncomps + c);
        const auto *gl_476 = buffer.data(gl + 476 * ncomps + c);
        const auto *gl_477 = buffer.data(gl + 477 * ncomps + c);
        const auto *gl_478 = buffer.data(gl + 478 * ncomps + c);
        const auto *gl_479 = buffer.data(gl + 479 * ncomps + c);
        const auto *gl_480 = buffer.data(gl + 480 * ncomps + c);
        const auto *gl_481 = buffer.data(gl + 481 * ncomps + c);
        const auto *gl_482 = buffer.data(gl + 482 * ncomps + c);
        const auto *gl_483 = buffer.data(gl + 483 * ncomps + c);
        const auto *gl_484 = buffer.data(gl + 484 * ncomps + c);
        const auto *gl_485 = buffer.data(gl + 485 * ncomps + c);
        const auto *gl_495 = buffer.data(gl + 495 * ncomps + c);
        const auto *gl_496 = buffer.data(gl + 496 * ncomps + c);
        const auto *gl_497 = buffer.data(gl + 497 * ncomps + c);
        const auto *gl_498 = buffer.data(gl + 498 * ncomps + c);
        const auto *gl_499 = buffer.data(gl + 499 * ncomps + c);
        const auto *gl_500 = buffer.data(gl + 500 * ncomps + c);
        const auto *gl_501 = buffer.data(gl + 501 * ncomps + c);
        const auto *gl_502 = buffer.data(gl + 502 * ncomps + c);
        const auto *gl_503 = buffer.data(gl + 503 * ncomps + c);
        const auto *gl_504 = buffer.data(gl + 504 * ncomps + c);
        const auto *gl_505 = buffer.data(gl + 505 * ncomps + c);
        const auto *gl_506 = buffer.data(gl + 506 * ncomps + c);
        const auto *gl_507 = buffer.data(gl + 507 * ncomps + c);
        const auto *gl_508 = buffer.data(gl + 508 * ncomps + c);
        const auto *gl_509 = buffer.data(gl + 509 * ncomps + c);
        const auto *gl_510 = buffer.data(gl + 510 * ncomps + c);
        const auto *gl_511 = buffer.data(gl + 511 * ncomps + c);
        const auto *gl_512 = buffer.data(gl + 512 * ncomps + c);
        const auto *gl_513 = buffer.data(gl + 513 * ncomps + c);
        const auto *gl_514 = buffer.data(gl + 514 * ncomps + c);
        const auto *gl_515 = buffer.data(gl + 515 * ncomps + c);
        const auto *gl_516 = buffer.data(gl + 516 * ncomps + c);
        const auto *gl_517 = buffer.data(gl + 517 * ncomps + c);
        const auto *gl_518 = buffer.data(gl + 518 * ncomps + c);
        const auto *gl_519 = buffer.data(gl + 519 * ncomps + c);
        const auto *gl_520 = buffer.data(gl + 520 * ncomps + c);
        const auto *gl_521 = buffer.data(gl + 521 * ncomps + c);
        const auto *gl_522 = buffer.data(gl + 522 * ncomps + c);
        const auto *gl_523 = buffer.data(gl + 523 * ncomps + c);
        const auto *gl_524 = buffer.data(gl + 524 * ncomps + c);
        const auto *gl_525 = buffer.data(gl + 525 * ncomps + c);
        const auto *gl_526 = buffer.data(gl + 526 * ncomps + c);
        const auto *gl_527 = buffer.data(gl + 527 * ncomps + c);
        const auto *gl_528 = buffer.data(gl + 528 * ncomps + c);
        const auto *gl_529 = buffer.data(gl + 529 * ncomps + c);
        const auto *gl_530 = buffer.data(gl + 530 * ncomps + c);
        const auto *gl_540 = buffer.data(gl + 540 * ncomps + c);
        const auto *gl_541 = buffer.data(gl + 541 * ncomps + c);
        const auto *gl_542 = buffer.data(gl + 542 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, gk_290, gk_291, gk_292, \
                         gk_293, gk_294, gl_362, gl_363, gl_364, gl_365, \
                         gl_366 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_x[k] * gk_290[k]
                       + gl_362[k];

            t_291[k] = -ab_x[k] * gk_291[k]
                       + gl_363[k];

            t_292[k] = -ab_x[k] * gk_292[k]
                       + gl_364[k];

            t_293[k] = -ab_x[k] * gk_293[k]
                       + gl_365[k];

            t_294[k] = -ab_x[k] * gk_294[k]
                       + gl_366[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, gk_295, gk_296, gk_297, \
                         gk_298, gk_299, gl_367, gl_368, gl_369, gl_370, \
                         gl_371 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_x[k] * gk_295[k]
                       + gl_367[k];

            t_296[k] = -ab_x[k] * gk_296[k]
                       + gl_368[k];

            t_297[k] = -ab_x[k] * gk_297[k]
                       + gl_369[k];

            t_298[k] = -ab_x[k] * gk_298[k]
                       + gl_370[k];

            t_299[k] = -ab_x[k] * gk_299[k]
                       + gl_371[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, gk_300, gk_301, gk_302, \
                         gk_303, gk_304, gl_372, gl_373, gl_374, gl_375, \
                         gl_376 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_x[k] * gk_300[k]
                       + gl_372[k];

            t_301[k] = -ab_x[k] * gk_301[k]
                       + gl_373[k];

            t_302[k] = -ab_x[k] * gk_302[k]
                       + gl_374[k];

            t_303[k] = -ab_x[k] * gk_303[k]
                       + gl_375[k];

            t_304[k] = -ab_x[k] * gk_304[k]
                       + gl_376[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, gk_305, gk_306, gk_307, \
                         gk_308, gk_309, gl_377, gl_378, gl_379, gl_380, \
                         gl_381 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_x[k] * gk_305[k]
                       + gl_377[k];

            t_306[k] = -ab_x[k] * gk_306[k]
                       + gl_378[k];

            t_307[k] = -ab_x[k] * gk_307[k]
                       + gl_379[k];

            t_308[k] = -ab_x[k] * gk_308[k]
                       + gl_380[k];

            t_309[k] = -ab_x[k] * gk_309[k]
                       + gl_381[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, gk_310, gk_311, gk_312, \
                         gk_313, gk_314, gl_382, gl_383, gl_384, gl_385, \
                         gl_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_x[k] * gk_310[k]
                       + gl_382[k];

            t_311[k] = -ab_x[k] * gk_311[k]
                       + gl_383[k];

            t_312[k] = -ab_x[k] * gk_312[k]
                       + gl_384[k];

            t_313[k] = -ab_x[k] * gk_313[k]
                       + gl_385[k];

            t_314[k] = -ab_x[k] * gk_314[k]
                       + gl_386[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, gk_315, gk_316, gk_317, \
                         gk_318, gk_319, gl_387, gl_388, gl_389, gl_390, \
                         gl_391 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_x[k] * gk_315[k]
                       + gl_387[k];

            t_316[k] = -ab_x[k] * gk_316[k]
                       + gl_388[k];

            t_317[k] = -ab_x[k] * gk_317[k]
                       + gl_389[k];

            t_318[k] = -ab_x[k] * gk_318[k]
                       + gl_390[k];

            t_319[k] = -ab_x[k] * gk_319[k]
                       + gl_391[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, gk_320, gk_321, gk_322, \
                         gk_323, gk_324, gl_392, gl_393, gl_394, gl_395, \
                         gl_405 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_x[k] * gk_320[k]
                       + gl_392[k];

            t_321[k] = -ab_x[k] * gk_321[k]
                       + gl_393[k];

            t_322[k] = -ab_x[k] * gk_322[k]
                       + gl_394[k];

            t_323[k] = -ab_x[k] * gk_323[k]
                       + gl_395[k];

            t_324[k] = -ab_x[k] * gk_324[k]
                       + gl_405[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, gk_325, gk_326, gk_327, \
                         gk_328, gk_329, gl_406, gl_407, gl_408, gl_409, \
                         gl_410 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_x[k] * gk_325[k]
                       + gl_406[k];

            t_326[k] = -ab_x[k] * gk_326[k]
                       + gl_407[k];

            t_327[k] = -ab_x[k] * gk_327[k]
                       + gl_408[k];

            t_328[k] = -ab_x[k] * gk_328[k]
                       + gl_409[k];

            t_329[k] = -ab_x[k] * gk_329[k]
                       + gl_410[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, gk_330, gk_331, gk_332, \
                         gk_333, gk_334, gl_411, gl_412, gl_413, gl_414, \
                         gl_415 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = -ab_x[k] * gk_330[k]
                       + gl_411[k];

            t_331[k] = -ab_x[k] * gk_331[k]
                       + gl_412[k];

            t_332[k] = -ab_x[k] * gk_332[k]
                       + gl_413[k];

            t_333[k] = -ab_x[k] * gk_333[k]
                       + gl_414[k];

            t_334[k] = -ab_x[k] * gk_334[k]
                       + gl_415[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, gk_335, gk_336, gk_337, \
                         gk_338, gk_339, gl_416, gl_417, gl_418, gl_419, \
                         gl_420 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = -ab_x[k] * gk_335[k]
                       + gl_416[k];

            t_336[k] = -ab_x[k] * gk_336[k]
                       + gl_417[k];

            t_337[k] = -ab_x[k] * gk_337[k]
                       + gl_418[k];

            t_338[k] = -ab_x[k] * gk_338[k]
                       + gl_419[k];

            t_339[k] = -ab_x[k] * gk_339[k]
                       + gl_420[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, gk_340, gk_341, gk_342, \
                         gk_343, gk_344, gl_421, gl_422, gl_423, gl_424, \
                         gl_425 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = -ab_x[k] * gk_340[k]
                       + gl_421[k];

            t_341[k] = -ab_x[k] * gk_341[k]
                       + gl_422[k];

            t_342[k] = -ab_x[k] * gk_342[k]
                       + gl_423[k];

            t_343[k] = -ab_x[k] * gk_343[k]
                       + gl_424[k];

            t_344[k] = -ab_x[k] * gk_344[k]
                       + gl_425[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, gk_345, gk_346, gk_347, \
                         gk_348, gk_349, gl_426, gl_427, gl_428, gl_429, \
                         gl_430 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = -ab_x[k] * gk_345[k]
                       + gl_426[k];

            t_346[k] = -ab_x[k] * gk_346[k]
                       + gl_427[k];

            t_347[k] = -ab_x[k] * gk_347[k]
                       + gl_428[k];

            t_348[k] = -ab_x[k] * gk_348[k]
                       + gl_429[k];

            t_349[k] = -ab_x[k] * gk_349[k]
                       + gl_430[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, gk_350, gk_351, gk_352, \
                         gk_353, gk_354, gl_431, gl_432, gl_433, gl_434, \
                         gl_435 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = -ab_x[k] * gk_350[k]
                       + gl_431[k];

            t_351[k] = -ab_x[k] * gk_351[k]
                       + gl_432[k];

            t_352[k] = -ab_x[k] * gk_352[k]
                       + gl_433[k];

            t_353[k] = -ab_x[k] * gk_353[k]
                       + gl_434[k];

            t_354[k] = -ab_x[k] * gk_354[k]
                       + gl_435[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, gk_355, gk_356, gk_357, \
                         gk_358, gk_359, gl_436, gl_437, gl_438, gl_439, \
                         gl_440 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = -ab_x[k] * gk_355[k]
                       + gl_436[k];

            t_356[k] = -ab_x[k] * gk_356[k]
                       + gl_437[k];

            t_357[k] = -ab_x[k] * gk_357[k]
                       + gl_438[k];

            t_358[k] = -ab_x[k] * gk_358[k]
                       + gl_439[k];

            t_359[k] = -ab_x[k] * gk_359[k]
                       + gl_440[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, gk_360, gk_361, gk_362, \
                         gk_363, gk_364, gl_450, gl_451, gl_452, gl_453, \
                         gl_454 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = -ab_x[k] * gk_360[k]
                       + gl_450[k];

            t_361[k] = -ab_x[k] * gk_361[k]
                       + gl_451[k];

            t_362[k] = -ab_x[k] * gk_362[k]
                       + gl_452[k];

            t_363[k] = -ab_x[k] * gk_363[k]
                       + gl_453[k];

            t_364[k] = -ab_x[k] * gk_364[k]
                       + gl_454[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, gk_365, gk_366, gk_367, \
                         gk_368, gk_369, gl_455, gl_456, gl_457, gl_458, \
                         gl_459 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = -ab_x[k] * gk_365[k]
                       + gl_455[k];

            t_366[k] = -ab_x[k] * gk_366[k]
                       + gl_456[k];

            t_367[k] = -ab_x[k] * gk_367[k]
                       + gl_457[k];

            t_368[k] = -ab_x[k] * gk_368[k]
                       + gl_458[k];

            t_369[k] = -ab_x[k] * gk_369[k]
                       + gl_459[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, gk_370, gk_371, gk_372, \
                         gk_373, gk_374, gl_460, gl_461, gl_462, gl_463, \
                         gl_464 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = -ab_x[k] * gk_370[k]
                       + gl_460[k];

            t_371[k] = -ab_x[k] * gk_371[k]
                       + gl_461[k];

            t_372[k] = -ab_x[k] * gk_372[k]
                       + gl_462[k];

            t_373[k] = -ab_x[k] * gk_373[k]
                       + gl_463[k];

            t_374[k] = -ab_x[k] * gk_374[k]
                       + gl_464[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, gk_375, gk_376, gk_377, \
                         gk_378, gk_379, gl_465, gl_466, gl_467, gl_468, \
                         gl_469 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = -ab_x[k] * gk_375[k]
                       + gl_465[k];

            t_376[k] = -ab_x[k] * gk_376[k]
                       + gl_466[k];

            t_377[k] = -ab_x[k] * gk_377[k]
                       + gl_467[k];

            t_378[k] = -ab_x[k] * gk_378[k]
                       + gl_468[k];

            t_379[k] = -ab_x[k] * gk_379[k]
                       + gl_469[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, gk_380, gk_381, gk_382, \
                         gk_383, gk_384, gl_470, gl_471, gl_472, gl_473, \
                         gl_474 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = -ab_x[k] * gk_380[k]
                       + gl_470[k];

            t_381[k] = -ab_x[k] * gk_381[k]
                       + gl_471[k];

            t_382[k] = -ab_x[k] * gk_382[k]
                       + gl_472[k];

            t_383[k] = -ab_x[k] * gk_383[k]
                       + gl_473[k];

            t_384[k] = -ab_x[k] * gk_384[k]
                       + gl_474[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, gk_385, gk_386, gk_387, \
                         gk_388, gk_389, gl_475, gl_476, gl_477, gl_478, \
                         gl_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = -ab_x[k] * gk_385[k]
                       + gl_475[k];

            t_386[k] = -ab_x[k] * gk_386[k]
                       + gl_476[k];

            t_387[k] = -ab_x[k] * gk_387[k]
                       + gl_477[k];

            t_388[k] = -ab_x[k] * gk_388[k]
                       + gl_478[k];

            t_389[k] = -ab_x[k] * gk_389[k]
                       + gl_479[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, gk_390, gk_391, gk_392, \
                         gk_393, gk_394, gl_480, gl_481, gl_482, gl_483, \
                         gl_484 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = -ab_x[k] * gk_390[k]
                       + gl_480[k];

            t_391[k] = -ab_x[k] * gk_391[k]
                       + gl_481[k];

            t_392[k] = -ab_x[k] * gk_392[k]
                       + gl_482[k];

            t_393[k] = -ab_x[k] * gk_393[k]
                       + gl_483[k];

            t_394[k] = -ab_x[k] * gk_394[k]
                       + gl_484[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, gk_395, gk_396, gk_397, \
                         gk_398, gk_399, gl_485, gl_495, gl_496, gl_497, \
                         gl_498 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = -ab_x[k] * gk_395[k]
                       + gl_485[k];

            t_396[k] = -ab_x[k] * gk_396[k]
                       + gl_495[k];

            t_397[k] = -ab_x[k] * gk_397[k]
                       + gl_496[k];

            t_398[k] = -ab_x[k] * gk_398[k]
                       + gl_497[k];

            t_399[k] = -ab_x[k] * gk_399[k]
                       + gl_498[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, gk_400, gk_401, gk_402, \
                         gk_403, gk_404, gl_499, gl_500, gl_501, gl_502, \
                         gl_503 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = -ab_x[k] * gk_400[k]
                       + gl_499[k];

            t_401[k] = -ab_x[k] * gk_401[k]
                       + gl_500[k];

            t_402[k] = -ab_x[k] * gk_402[k]
                       + gl_501[k];

            t_403[k] = -ab_x[k] * gk_403[k]
                       + gl_502[k];

            t_404[k] = -ab_x[k] * gk_404[k]
                       + gl_503[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, gk_405, gk_406, gk_407, \
                         gk_408, gk_409, gl_504, gl_505, gl_506, gl_507, \
                         gl_508 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = -ab_x[k] * gk_405[k]
                       + gl_504[k];

            t_406[k] = -ab_x[k] * gk_406[k]
                       + gl_505[k];

            t_407[k] = -ab_x[k] * gk_407[k]
                       + gl_506[k];

            t_408[k] = -ab_x[k] * gk_408[k]
                       + gl_507[k];

            t_409[k] = -ab_x[k] * gk_409[k]
                       + gl_508[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, gk_410, gk_411, gk_412, \
                         gk_413, gk_414, gl_509, gl_510, gl_511, gl_512, \
                         gl_513 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = -ab_x[k] * gk_410[k]
                       + gl_509[k];

            t_411[k] = -ab_x[k] * gk_411[k]
                       + gl_510[k];

            t_412[k] = -ab_x[k] * gk_412[k]
                       + gl_511[k];

            t_413[k] = -ab_x[k] * gk_413[k]
                       + gl_512[k];

            t_414[k] = -ab_x[k] * gk_414[k]
                       + gl_513[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, gk_415, gk_416, gk_417, \
                         gk_418, gk_419, gl_514, gl_515, gl_516, gl_517, \
                         gl_518 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = -ab_x[k] * gk_415[k]
                       + gl_514[k];

            t_416[k] = -ab_x[k] * gk_416[k]
                       + gl_515[k];

            t_417[k] = -ab_x[k] * gk_417[k]
                       + gl_516[k];

            t_418[k] = -ab_x[k] * gk_418[k]
                       + gl_517[k];

            t_419[k] = -ab_x[k] * gk_419[k]
                       + gl_518[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, gk_420, gk_421, gk_422, \
                         gk_423, gk_424, gl_519, gl_520, gl_521, gl_522, \
                         gl_523 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = -ab_x[k] * gk_420[k]
                       + gl_519[k];

            t_421[k] = -ab_x[k] * gk_421[k]
                       + gl_520[k];

            t_422[k] = -ab_x[k] * gk_422[k]
                       + gl_521[k];

            t_423[k] = -ab_x[k] * gk_423[k]
                       + gl_522[k];

            t_424[k] = -ab_x[k] * gk_424[k]
                       + gl_523[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, gk_425, gk_426, gk_427, \
                         gk_428, gk_429, gl_524, gl_525, gl_526, gl_527, \
                         gl_528 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = -ab_x[k] * gk_425[k]
                       + gl_524[k];

            t_426[k] = -ab_x[k] * gk_426[k]
                       + gl_525[k];

            t_427[k] = -ab_x[k] * gk_427[k]
                       + gl_526[k];

            t_428[k] = -ab_x[k] * gk_428[k]
                       + gl_527[k];

            t_429[k] = -ab_x[k] * gk_429[k]
                       + gl_528[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, gk_430, gk_431, gk_432, \
                         gk_433, gk_434, gl_529, gl_530, gl_540, gl_541, \
                         gl_542 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = -ab_x[k] * gk_430[k]
                       + gl_529[k];

            t_431[k] = -ab_x[k] * gk_431[k]
                       + gl_530[k];

            t_432[k] = -ab_x[k] * gk_432[k]
                       + gl_540[k];

            t_433[k] = -ab_x[k] * gk_433[k]
                       + gl_541[k];

            t_434[k] = -ab_x[k] * gk_434[k]
                       + gl_542[k];
        }
    }
}

static auto
compute_hrr_hk_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gk, const size_t gl, const size_t ncomps,
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
        const auto *ab_y = coordinates.data(7);

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
        const auto *gk_388 = buffer.data(gk + 388 * ncomps + c);
        const auto *gk_389 = buffer.data(gk + 389 * ncomps + c);
        const auto *gk_390 = buffer.data(gk + 390 * ncomps + c);
        const auto *gk_391 = buffer.data(gk + 391 * ncomps + c);
        const auto *gk_392 = buffer.data(gk + 392 * ncomps + c);
        const auto *gk_393 = buffer.data(gk + 393 * ncomps + c);
        const auto *gk_394 = buffer.data(gk + 394 * ncomps + c);
        const auto *gk_395 = buffer.data(gk + 395 * ncomps + c);
        const auto *gk_396 = buffer.data(gk + 396 * ncomps + c);
        const auto *gk_397 = buffer.data(gk + 397 * ncomps + c);
        const auto *gk_398 = buffer.data(gk + 398 * ncomps + c);
        const auto *gk_399 = buffer.data(gk + 399 * ncomps + c);
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
        const auto *gk_460 = buffer.data(gk + 460 * ncomps + c);
        const auto *gk_461 = buffer.data(gk + 461 * ncomps + c);
        const auto *gk_462 = buffer.data(gk + 462 * ncomps + c);
        const auto *gk_463 = buffer.data(gk + 463 * ncomps + c);
        const auto *gk_464 = buffer.data(gk + 464 * ncomps + c);
        const auto *gk_465 = buffer.data(gk + 465 * ncomps + c);
        const auto *gk_466 = buffer.data(gk + 466 * ncomps + c);
        const auto *gk_467 = buffer.data(gk + 467 * ncomps + c);
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
        const auto *gk_496 = buffer.data(gk + 496 * ncomps + c);
        const auto *gk_497 = buffer.data(gk + 497 * ncomps + c);
        const auto *gk_498 = buffer.data(gk + 498 * ncomps + c);
        const auto *gk_499 = buffer.data(gk + 499 * ncomps + c);
        const auto *gk_500 = buffer.data(gk + 500 * ncomps + c);
        const auto *gk_501 = buffer.data(gk + 501 * ncomps + c);
        const auto *gk_502 = buffer.data(gk + 502 * ncomps + c);
        const auto *gk_503 = buffer.data(gk + 503 * ncomps + c);
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
        const auto *gk_532 = buffer.data(gk + 532 * ncomps + c);
        const auto *gk_533 = buffer.data(gk + 533 * ncomps + c);
        const auto *gk_534 = buffer.data(gk + 534 * ncomps + c);
        const auto *gk_535 = buffer.data(gk + 535 * ncomps + c);
        const auto *gk_536 = buffer.data(gk + 536 * ncomps + c);
        const auto *gk_537 = buffer.data(gk + 537 * ncomps + c);
        const auto *gk_538 = buffer.data(gk + 538 * ncomps + c);
        const auto *gk_539 = buffer.data(gk + 539 * ncomps + c);

        const auto *gl_451 = buffer.data(gl + 451 * ncomps + c);
        const auto *gl_453 = buffer.data(gl + 453 * ncomps + c);
        const auto *gl_454 = buffer.data(gl + 454 * ncomps + c);
        const auto *gl_456 = buffer.data(gl + 456 * ncomps + c);
        const auto *gl_457 = buffer.data(gl + 457 * ncomps + c);
        const auto *gl_458 = buffer.data(gl + 458 * ncomps + c);
        const auto *gl_460 = buffer.data(gl + 460 * ncomps + c);
        const auto *gl_461 = buffer.data(gl + 461 * ncomps + c);
        const auto *gl_462 = buffer.data(gl + 462 * ncomps + c);
        const auto *gl_463 = buffer.data(gl + 463 * ncomps + c);
        const auto *gl_465 = buffer.data(gl + 465 * ncomps + c);
        const auto *gl_466 = buffer.data(gl + 466 * ncomps + c);
        const auto *gl_467 = buffer.data(gl + 467 * ncomps + c);
        const auto *gl_468 = buffer.data(gl + 468 * ncomps + c);
        const auto *gl_469 = buffer.data(gl + 469 * ncomps + c);
        const auto *gl_471 = buffer.data(gl + 471 * ncomps + c);
        const auto *gl_472 = buffer.data(gl + 472 * ncomps + c);
        const auto *gl_473 = buffer.data(gl + 473 * ncomps + c);
        const auto *gl_474 = buffer.data(gl + 474 * ncomps + c);
        const auto *gl_475 = buffer.data(gl + 475 * ncomps + c);
        const auto *gl_476 = buffer.data(gl + 476 * ncomps + c);
        const auto *gl_478 = buffer.data(gl + 478 * ncomps + c);
        const auto *gl_479 = buffer.data(gl + 479 * ncomps + c);
        const auto *gl_480 = buffer.data(gl + 480 * ncomps + c);
        const auto *gl_481 = buffer.data(gl + 481 * ncomps + c);
        const auto *gl_482 = buffer.data(gl + 482 * ncomps + c);
        const auto *gl_483 = buffer.data(gl + 483 * ncomps + c);
        const auto *gl_484 = buffer.data(gl + 484 * ncomps + c);
        const auto *gl_486 = buffer.data(gl + 486 * ncomps + c);
        const auto *gl_487 = buffer.data(gl + 487 * ncomps + c);
        const auto *gl_488 = buffer.data(gl + 488 * ncomps + c);
        const auto *gl_489 = buffer.data(gl + 489 * ncomps + c);
        const auto *gl_490 = buffer.data(gl + 490 * ncomps + c);
        const auto *gl_491 = buffer.data(gl + 491 * ncomps + c);
        const auto *gl_492 = buffer.data(gl + 492 * ncomps + c);
        const auto *gl_493 = buffer.data(gl + 493 * ncomps + c);
        const auto *gl_496 = buffer.data(gl + 496 * ncomps + c);
        const auto *gl_498 = buffer.data(gl + 498 * ncomps + c);
        const auto *gl_499 = buffer.data(gl + 499 * ncomps + c);
        const auto *gl_501 = buffer.data(gl + 501 * ncomps + c);
        const auto *gl_543 = buffer.data(gl + 543 * ncomps + c);
        const auto *gl_544 = buffer.data(gl + 544 * ncomps + c);
        const auto *gl_545 = buffer.data(gl + 545 * ncomps + c);
        const auto *gl_546 = buffer.data(gl + 546 * ncomps + c);
        const auto *gl_547 = buffer.data(gl + 547 * ncomps + c);
        const auto *gl_548 = buffer.data(gl + 548 * ncomps + c);
        const auto *gl_549 = buffer.data(gl + 549 * ncomps + c);
        const auto *gl_550 = buffer.data(gl + 550 * ncomps + c);
        const auto *gl_551 = buffer.data(gl + 551 * ncomps + c);
        const auto *gl_552 = buffer.data(gl + 552 * ncomps + c);
        const auto *gl_553 = buffer.data(gl + 553 * ncomps + c);
        const auto *gl_554 = buffer.data(gl + 554 * ncomps + c);
        const auto *gl_555 = buffer.data(gl + 555 * ncomps + c);
        const auto *gl_556 = buffer.data(gl + 556 * ncomps + c);
        const auto *gl_557 = buffer.data(gl + 557 * ncomps + c);
        const auto *gl_558 = buffer.data(gl + 558 * ncomps + c);
        const auto *gl_559 = buffer.data(gl + 559 * ncomps + c);
        const auto *gl_560 = buffer.data(gl + 560 * ncomps + c);
        const auto *gl_561 = buffer.data(gl + 561 * ncomps + c);
        const auto *gl_562 = buffer.data(gl + 562 * ncomps + c);
        const auto *gl_563 = buffer.data(gl + 563 * ncomps + c);
        const auto *gl_564 = buffer.data(gl + 564 * ncomps + c);
        const auto *gl_565 = buffer.data(gl + 565 * ncomps + c);
        const auto *gl_566 = buffer.data(gl + 566 * ncomps + c);
        const auto *gl_567 = buffer.data(gl + 567 * ncomps + c);
        const auto *gl_568 = buffer.data(gl + 568 * ncomps + c);
        const auto *gl_569 = buffer.data(gl + 569 * ncomps + c);
        const auto *gl_570 = buffer.data(gl + 570 * ncomps + c);
        const auto *gl_571 = buffer.data(gl + 571 * ncomps + c);
        const auto *gl_572 = buffer.data(gl + 572 * ncomps + c);
        const auto *gl_573 = buffer.data(gl + 573 * ncomps + c);
        const auto *gl_574 = buffer.data(gl + 574 * ncomps + c);
        const auto *gl_575 = buffer.data(gl + 575 * ncomps + c);
        const auto *gl_585 = buffer.data(gl + 585 * ncomps + c);
        const auto *gl_586 = buffer.data(gl + 586 * ncomps + c);
        const auto *gl_587 = buffer.data(gl + 587 * ncomps + c);
        const auto *gl_588 = buffer.data(gl + 588 * ncomps + c);
        const auto *gl_589 = buffer.data(gl + 589 * ncomps + c);
        const auto *gl_590 = buffer.data(gl + 590 * ncomps + c);
        const auto *gl_591 = buffer.data(gl + 591 * ncomps + c);
        const auto *gl_592 = buffer.data(gl + 592 * ncomps + c);
        const auto *gl_593 = buffer.data(gl + 593 * ncomps + c);
        const auto *gl_594 = buffer.data(gl + 594 * ncomps + c);
        const auto *gl_595 = buffer.data(gl + 595 * ncomps + c);
        const auto *gl_596 = buffer.data(gl + 596 * ncomps + c);
        const auto *gl_597 = buffer.data(gl + 597 * ncomps + c);
        const auto *gl_598 = buffer.data(gl + 598 * ncomps + c);
        const auto *gl_599 = buffer.data(gl + 599 * ncomps + c);
        const auto *gl_600 = buffer.data(gl + 600 * ncomps + c);
        const auto *gl_601 = buffer.data(gl + 601 * ncomps + c);
        const auto *gl_602 = buffer.data(gl + 602 * ncomps + c);
        const auto *gl_603 = buffer.data(gl + 603 * ncomps + c);
        const auto *gl_604 = buffer.data(gl + 604 * ncomps + c);
        const auto *gl_605 = buffer.data(gl + 605 * ncomps + c);
        const auto *gl_606 = buffer.data(gl + 606 * ncomps + c);
        const auto *gl_607 = buffer.data(gl + 607 * ncomps + c);
        const auto *gl_608 = buffer.data(gl + 608 * ncomps + c);
        const auto *gl_609 = buffer.data(gl + 609 * ncomps + c);
        const auto *gl_610 = buffer.data(gl + 610 * ncomps + c);
        const auto *gl_611 = buffer.data(gl + 611 * ncomps + c);
        const auto *gl_612 = buffer.data(gl + 612 * ncomps + c);
        const auto *gl_613 = buffer.data(gl + 613 * ncomps + c);
        const auto *gl_614 = buffer.data(gl + 614 * ncomps + c);
        const auto *gl_615 = buffer.data(gl + 615 * ncomps + c);
        const auto *gl_616 = buffer.data(gl + 616 * ncomps + c);
        const auto *gl_617 = buffer.data(gl + 617 * ncomps + c);
        const auto *gl_618 = buffer.data(gl + 618 * ncomps + c);
        const auto *gl_619 = buffer.data(gl + 619 * ncomps + c);
        const auto *gl_620 = buffer.data(gl + 620 * ncomps + c);
        const auto *gl_630 = buffer.data(gl + 630 * ncomps + c);
        const auto *gl_631 = buffer.data(gl + 631 * ncomps + c);
        const auto *gl_632 = buffer.data(gl + 632 * ncomps + c);
        const auto *gl_633 = buffer.data(gl + 633 * ncomps + c);
        const auto *gl_634 = buffer.data(gl + 634 * ncomps + c);
        const auto *gl_635 = buffer.data(gl + 635 * ncomps + c);
        const auto *gl_636 = buffer.data(gl + 636 * ncomps + c);
        const auto *gl_637 = buffer.data(gl + 637 * ncomps + c);
        const auto *gl_638 = buffer.data(gl + 638 * ncomps + c);
        const auto *gl_639 = buffer.data(gl + 639 * ncomps + c);
        const auto *gl_640 = buffer.data(gl + 640 * ncomps + c);
        const auto *gl_641 = buffer.data(gl + 641 * ncomps + c);
        const auto *gl_642 = buffer.data(gl + 642 * ncomps + c);
        const auto *gl_643 = buffer.data(gl + 643 * ncomps + c);
        const auto *gl_644 = buffer.data(gl + 644 * ncomps + c);
        const auto *gl_645 = buffer.data(gl + 645 * ncomps + c);
        const auto *gl_646 = buffer.data(gl + 646 * ncomps + c);
        const auto *gl_647 = buffer.data(gl + 647 * ncomps + c);
        const auto *gl_648 = buffer.data(gl + 648 * ncomps + c);
        const auto *gl_649 = buffer.data(gl + 649 * ncomps + c);
        const auto *gl_650 = buffer.data(gl + 650 * ncomps + c);
        const auto *gl_651 = buffer.data(gl + 651 * ncomps + c);
        const auto *gl_652 = buffer.data(gl + 652 * ncomps + c);
        const auto *gl_653 = buffer.data(gl + 653 * ncomps + c);
        const auto *gl_654 = buffer.data(gl + 654 * ncomps + c);
        const auto *gl_655 = buffer.data(gl + 655 * ncomps + c);
        const auto *gl_656 = buffer.data(gl + 656 * ncomps + c);
        const auto *gl_657 = buffer.data(gl + 657 * ncomps + c);
        const auto *gl_658 = buffer.data(gl + 658 * ncomps + c);
        const auto *gl_659 = buffer.data(gl + 659 * ncomps + c);
        const auto *gl_660 = buffer.data(gl + 660 * ncomps + c);
        const auto *gl_661 = buffer.data(gl + 661 * ncomps + c);
        const auto *gl_662 = buffer.data(gl + 662 * ncomps + c);
        const auto *gl_663 = buffer.data(gl + 663 * ncomps + c);
        const auto *gl_664 = buffer.data(gl + 664 * ncomps + c);
        const auto *gl_665 = buffer.data(gl + 665 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, gk_435, gk_436, gk_437, \
                         gk_438, gk_439, gl_543, gl_544, gl_545, gl_546, \
                         gl_547 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = -ab_x[k] * gk_435[k]
                       + gl_543[k];

            t_436[k] = -ab_x[k] * gk_436[k]
                       + gl_544[k];

            t_437[k] = -ab_x[k] * gk_437[k]
                       + gl_545[k];

            t_438[k] = -ab_x[k] * gk_438[k]
                       + gl_546[k];

            t_439[k] = -ab_x[k] * gk_439[k]
                       + gl_547[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, gk_440, gk_441, gk_442, \
                         gk_443, gk_444, gl_548, gl_549, gl_550, gl_551, \
                         gl_552 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = -ab_x[k] * gk_440[k]
                       + gl_548[k];

            t_441[k] = -ab_x[k] * gk_441[k]
                       + gl_549[k];

            t_442[k] = -ab_x[k] * gk_442[k]
                       + gl_550[k];

            t_443[k] = -ab_x[k] * gk_443[k]
                       + gl_551[k];

            t_444[k] = -ab_x[k] * gk_444[k]
                       + gl_552[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_x, gk_445, gk_446, gk_447, \
                         gk_448, gk_449, gl_553, gl_554, gl_555, gl_556, \
                         gl_557 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = -ab_x[k] * gk_445[k]
                       + gl_553[k];

            t_446[k] = -ab_x[k] * gk_446[k]
                       + gl_554[k];

            t_447[k] = -ab_x[k] * gk_447[k]
                       + gl_555[k];

            t_448[k] = -ab_x[k] * gk_448[k]
                       + gl_556[k];

            t_449[k] = -ab_x[k] * gk_449[k]
                       + gl_557[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_x, gk_450, gk_451, gk_452, \
                         gk_453, gk_454, gl_558, gl_559, gl_560, gl_561, \
                         gl_562 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = -ab_x[k] * gk_450[k]
                       + gl_558[k];

            t_451[k] = -ab_x[k] * gk_451[k]
                       + gl_559[k];

            t_452[k] = -ab_x[k] * gk_452[k]
                       + gl_560[k];

            t_453[k] = -ab_x[k] * gk_453[k]
                       + gl_561[k];

            t_454[k] = -ab_x[k] * gk_454[k]
                       + gl_562[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_x, gk_455, gk_456, gk_457, \
                         gk_458, gk_459, gl_563, gl_564, gl_565, gl_566, \
                         gl_567 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = -ab_x[k] * gk_455[k]
                       + gl_563[k];

            t_456[k] = -ab_x[k] * gk_456[k]
                       + gl_564[k];

            t_457[k] = -ab_x[k] * gk_457[k]
                       + gl_565[k];

            t_458[k] = -ab_x[k] * gk_458[k]
                       + gl_566[k];

            t_459[k] = -ab_x[k] * gk_459[k]
                       + gl_567[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_x, gk_460, gk_461, gk_462, \
                         gk_463, gk_464, gl_568, gl_569, gl_570, gl_571, \
                         gl_572 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = -ab_x[k] * gk_460[k]
                       + gl_568[k];

            t_461[k] = -ab_x[k] * gk_461[k]
                       + gl_569[k];

            t_462[k] = -ab_x[k] * gk_462[k]
                       + gl_570[k];

            t_463[k] = -ab_x[k] * gk_463[k]
                       + gl_571[k];

            t_464[k] = -ab_x[k] * gk_464[k]
                       + gl_572[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_x, gk_465, gk_466, gk_467, \
                         gk_468, gk_469, gl_573, gl_574, gl_575, gl_585, \
                         gl_586 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = -ab_x[k] * gk_465[k]
                       + gl_573[k];

            t_466[k] = -ab_x[k] * gk_466[k]
                       + gl_574[k];

            t_467[k] = -ab_x[k] * gk_467[k]
                       + gl_575[k];

            t_468[k] = -ab_x[k] * gk_468[k]
                       + gl_585[k];

            t_469[k] = -ab_x[k] * gk_469[k]
                       + gl_586[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_x, gk_470, gk_471, gk_472, \
                         gk_473, gk_474, gl_587, gl_588, gl_589, gl_590, \
                         gl_591 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = -ab_x[k] * gk_470[k]
                       + gl_587[k];

            t_471[k] = -ab_x[k] * gk_471[k]
                       + gl_588[k];

            t_472[k] = -ab_x[k] * gk_472[k]
                       + gl_589[k];

            t_473[k] = -ab_x[k] * gk_473[k]
                       + gl_590[k];

            t_474[k] = -ab_x[k] * gk_474[k]
                       + gl_591[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_x, gk_475, gk_476, gk_477, \
                         gk_478, gk_479, gl_592, gl_593, gl_594, gl_595, \
                         gl_596 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = -ab_x[k] * gk_475[k]
                       + gl_592[k];

            t_476[k] = -ab_x[k] * gk_476[k]
                       + gl_593[k];

            t_477[k] = -ab_x[k] * gk_477[k]
                       + gl_594[k];

            t_478[k] = -ab_x[k] * gk_478[k]
                       + gl_595[k];

            t_479[k] = -ab_x[k] * gk_479[k]
                       + gl_596[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_x, gk_480, gk_481, gk_482, \
                         gk_483, gk_484, gl_597, gl_598, gl_599, gl_600, \
                         gl_601 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = -ab_x[k] * gk_480[k]
                       + gl_597[k];

            t_481[k] = -ab_x[k] * gk_481[k]
                       + gl_598[k];

            t_482[k] = -ab_x[k] * gk_482[k]
                       + gl_599[k];

            t_483[k] = -ab_x[k] * gk_483[k]
                       + gl_600[k];

            t_484[k] = -ab_x[k] * gk_484[k]
                       + gl_601[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_x, gk_485, gk_486, gk_487, \
                         gk_488, gk_489, gl_602, gl_603, gl_604, gl_605, \
                         gl_606 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = -ab_x[k] * gk_485[k]
                       + gl_602[k];

            t_486[k] = -ab_x[k] * gk_486[k]
                       + gl_603[k];

            t_487[k] = -ab_x[k] * gk_487[k]
                       + gl_604[k];

            t_488[k] = -ab_x[k] * gk_488[k]
                       + gl_605[k];

            t_489[k] = -ab_x[k] * gk_489[k]
                       + gl_606[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_x, gk_490, gk_491, gk_492, \
                         gk_493, gk_494, gl_607, gl_608, gl_609, gl_610, \
                         gl_611 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = -ab_x[k] * gk_490[k]
                       + gl_607[k];

            t_491[k] = -ab_x[k] * gk_491[k]
                       + gl_608[k];

            t_492[k] = -ab_x[k] * gk_492[k]
                       + gl_609[k];

            t_493[k] = -ab_x[k] * gk_493[k]
                       + gl_610[k];

            t_494[k] = -ab_x[k] * gk_494[k]
                       + gl_611[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_x, gk_495, gk_496, gk_497, \
                         gk_498, gk_499, gl_612, gl_613, gl_614, gl_615, \
                         gl_616 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = -ab_x[k] * gk_495[k]
                       + gl_612[k];

            t_496[k] = -ab_x[k] * gk_496[k]
                       + gl_613[k];

            t_497[k] = -ab_x[k] * gk_497[k]
                       + gl_614[k];

            t_498[k] = -ab_x[k] * gk_498[k]
                       + gl_615[k];

            t_499[k] = -ab_x[k] * gk_499[k]
                       + gl_616[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_x, gk_500, gk_501, gk_502, \
                         gk_503, gk_504, gl_617, gl_618, gl_619, gl_620, \
                         gl_630 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = -ab_x[k] * gk_500[k]
                       + gl_617[k];

            t_501[k] = -ab_x[k] * gk_501[k]
                       + gl_618[k];

            t_502[k] = -ab_x[k] * gk_502[k]
                       + gl_619[k];

            t_503[k] = -ab_x[k] * gk_503[k]
                       + gl_620[k];

            t_504[k] = -ab_x[k] * gk_504[k]
                       + gl_630[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_x, gk_505, gk_506, gk_507, \
                         gk_508, gk_509, gl_631, gl_632, gl_633, gl_634, \
                         gl_635 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = -ab_x[k] * gk_505[k]
                       + gl_631[k];

            t_506[k] = -ab_x[k] * gk_506[k]
                       + gl_632[k];

            t_507[k] = -ab_x[k] * gk_507[k]
                       + gl_633[k];

            t_508[k] = -ab_x[k] * gk_508[k]
                       + gl_634[k];

            t_509[k] = -ab_x[k] * gk_509[k]
                       + gl_635[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_x, gk_510, gk_511, gk_512, \
                         gk_513, gk_514, gl_636, gl_637, gl_638, gl_639, \
                         gl_640 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = -ab_x[k] * gk_510[k]
                       + gl_636[k];

            t_511[k] = -ab_x[k] * gk_511[k]
                       + gl_637[k];

            t_512[k] = -ab_x[k] * gk_512[k]
                       + gl_638[k];

            t_513[k] = -ab_x[k] * gk_513[k]
                       + gl_639[k];

            t_514[k] = -ab_x[k] * gk_514[k]
                       + gl_640[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_x, gk_515, gk_516, gk_517, \
                         gk_518, gk_519, gl_641, gl_642, gl_643, gl_644, \
                         gl_645 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = -ab_x[k] * gk_515[k]
                       + gl_641[k];

            t_516[k] = -ab_x[k] * gk_516[k]
                       + gl_642[k];

            t_517[k] = -ab_x[k] * gk_517[k]
                       + gl_643[k];

            t_518[k] = -ab_x[k] * gk_518[k]
                       + gl_644[k];

            t_519[k] = -ab_x[k] * gk_519[k]
                       + gl_645[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_x, gk_520, gk_521, gk_522, \
                         gk_523, gk_524, gl_646, gl_647, gl_648, gl_649, \
                         gl_650 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = -ab_x[k] * gk_520[k]
                       + gl_646[k];

            t_521[k] = -ab_x[k] * gk_521[k]
                       + gl_647[k];

            t_522[k] = -ab_x[k] * gk_522[k]
                       + gl_648[k];

            t_523[k] = -ab_x[k] * gk_523[k]
                       + gl_649[k];

            t_524[k] = -ab_x[k] * gk_524[k]
                       + gl_650[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_x, gk_525, gk_526, gk_527, \
                         gk_528, gk_529, gl_651, gl_652, gl_653, gl_654, \
                         gl_655 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = -ab_x[k] * gk_525[k]
                       + gl_651[k];

            t_526[k] = -ab_x[k] * gk_526[k]
                       + gl_652[k];

            t_527[k] = -ab_x[k] * gk_527[k]
                       + gl_653[k];

            t_528[k] = -ab_x[k] * gk_528[k]
                       + gl_654[k];

            t_529[k] = -ab_x[k] * gk_529[k]
                       + gl_655[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_x, gk_530, gk_531, gk_532, \
                         gk_533, gk_534, gl_656, gl_657, gl_658, gl_659, \
                         gl_660 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = -ab_x[k] * gk_530[k]
                       + gl_656[k];

            t_531[k] = -ab_x[k] * gk_531[k]
                       + gl_657[k];

            t_532[k] = -ab_x[k] * gk_532[k]
                       + gl_658[k];

            t_533[k] = -ab_x[k] * gk_533[k]
                       + gl_659[k];

            t_534[k] = -ab_x[k] * gk_534[k]
                       + gl_660[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_x, gk_535, gk_536, gk_537, \
                         gk_538, gk_539, gl_661, gl_662, gl_663, gl_664, \
                         gl_665 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = -ab_x[k] * gk_535[k]
                       + gl_661[k];

            t_536[k] = -ab_x[k] * gk_536[k]
                       + gl_662[k];

            t_537[k] = -ab_x[k] * gk_537[k]
                       + gl_663[k];

            t_538[k] = -ab_x[k] * gk_538[k]
                       + gl_664[k];

            t_539[k] = -ab_x[k] * gk_539[k]
                       + gl_665[k];
        }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_y, gk_360, gk_361, gk_362, \
                         gk_363, gk_364, gl_451, gl_453, gl_454, gl_456, \
                         gl_457 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_540[k] = -ab_y[k] * gk_360[k]
                       + gl_451[k];

            t_541[k] = -ab_y[k] * gk_361[k]
                       + gl_453[k];

            t_542[k] = -ab_y[k] * gk_362[k]
                       + gl_454[k];

            t_543[k] = -ab_y[k] * gk_363[k]
                       + gl_456[k];

            t_544[k] = -ab_y[k] * gk_364[k]
                       + gl_457[k];
        }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_y, gk_365, gk_366, gk_367, \
                         gk_368, gk_369, gl_458, gl_460, gl_461, gl_462, \
                         gl_463 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_545[k] = -ab_y[k] * gk_365[k]
                       + gl_458[k];

            t_546[k] = -ab_y[k] * gk_366[k]
                       + gl_460[k];

            t_547[k] = -ab_y[k] * gk_367[k]
                       + gl_461[k];

            t_548[k] = -ab_y[k] * gk_368[k]
                       + gl_462[k];

            t_549[k] = -ab_y[k] * gk_369[k]
                       + gl_463[k];
        }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ab_y, gk_370, gk_371, gk_372, \
                         gk_373, gk_374, gl_465, gl_466, gl_467, gl_468, \
                         gl_469 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_550[k] = -ab_y[k] * gk_370[k]
                       + gl_465[k];

            t_551[k] = -ab_y[k] * gk_371[k]
                       + gl_466[k];

            t_552[k] = -ab_y[k] * gk_372[k]
                       + gl_467[k];

            t_553[k] = -ab_y[k] * gk_373[k]
                       + gl_468[k];

            t_554[k] = -ab_y[k] * gk_374[k]
                       + gl_469[k];
        }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ab_y, gk_375, gk_376, gk_377, \
                         gk_378, gk_379, gl_471, gl_472, gl_473, gl_474, \
                         gl_475 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_555[k] = -ab_y[k] * gk_375[k]
                       + gl_471[k];

            t_556[k] = -ab_y[k] * gk_376[k]
                       + gl_472[k];

            t_557[k] = -ab_y[k] * gk_377[k]
                       + gl_473[k];

            t_558[k] = -ab_y[k] * gk_378[k]
                       + gl_474[k];

            t_559[k] = -ab_y[k] * gk_379[k]
                       + gl_475[k];
        }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_y, gk_380, gk_381, gk_382, \
                         gk_383, gk_384, gl_476, gl_478, gl_479, gl_480, \
                         gl_481 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_560[k] = -ab_y[k] * gk_380[k]
                       + gl_476[k];

            t_561[k] = -ab_y[k] * gk_381[k]
                       + gl_478[k];

            t_562[k] = -ab_y[k] * gk_382[k]
                       + gl_479[k];

            t_563[k] = -ab_y[k] * gk_383[k]
                       + gl_480[k];

            t_564[k] = -ab_y[k] * gk_384[k]
                       + gl_481[k];
        }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_y, gk_385, gk_386, gk_387, \
                         gk_388, gk_389, gl_482, gl_483, gl_484, gl_486, \
                         gl_487 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_565[k] = -ab_y[k] * gk_385[k]
                       + gl_482[k];

            t_566[k] = -ab_y[k] * gk_386[k]
                       + gl_483[k];

            t_567[k] = -ab_y[k] * gk_387[k]
                       + gl_484[k];

            t_568[k] = -ab_y[k] * gk_388[k]
                       + gl_486[k];

            t_569[k] = -ab_y[k] * gk_389[k]
                       + gl_487[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_y, gk_390, gk_391, gk_392, \
                         gk_393, gk_394, gl_488, gl_489, gl_490, gl_491, \
                         gl_492 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = -ab_y[k] * gk_390[k]
                       + gl_488[k];

            t_571[k] = -ab_y[k] * gk_391[k]
                       + gl_489[k];

            t_572[k] = -ab_y[k] * gk_392[k]
                       + gl_490[k];

            t_573[k] = -ab_y[k] * gk_393[k]
                       + gl_491[k];

            t_574[k] = -ab_y[k] * gk_394[k]
                       + gl_492[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_y, gk_395, gk_396, gk_397, \
                         gk_398, gk_399, gl_493, gl_496, gl_498, gl_499, \
                         gl_501 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = -ab_y[k] * gk_395[k]
                       + gl_493[k];

            t_576[k] = -ab_y[k] * gk_396[k]
                       + gl_496[k];

            t_577[k] = -ab_y[k] * gk_397[k]
                       + gl_498[k];

            t_578[k] = -ab_y[k] * gk_398[k]
                       + gl_499[k];

            t_579[k] = -ab_y[k] * gk_399[k]
                       + gl_501[k];
        }
    }
}

static auto
compute_hrr_hk_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gk, const size_t gl, const size_t ncomps,
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
        auto *t_724 = buffer.data(target + 724 * ncomps + c);

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

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
        const auto *gk_424 = buffer.data(gk + 424 * ncomps + c);
        const auto *gk_425 = buffer.data(gk + 425 * ncomps + c);
        const auto *gk_426 = buffer.data(gk + 426 * ncomps + c);
        const auto *gk_427 = buffer.data(gk + 427 * ncomps + c);
        const auto *gk_428 = buffer.data(gk + 428 * ncomps + c);
        const auto *gk_429 = buffer.data(gk + 429 * ncomps + c);
        const auto *gk_430 = buffer.data(gk + 430 * ncomps + c);
        const auto *gk_431 = buffer.data(gk + 431 * ncomps + c);
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
        const auto *gk_460 = buffer.data(gk + 460 * ncomps + c);
        const auto *gk_461 = buffer.data(gk + 461 * ncomps + c);
        const auto *gk_462 = buffer.data(gk + 462 * ncomps + c);
        const auto *gk_463 = buffer.data(gk + 463 * ncomps + c);
        const auto *gk_464 = buffer.data(gk + 464 * ncomps + c);
        const auto *gk_465 = buffer.data(gk + 465 * ncomps + c);
        const auto *gk_466 = buffer.data(gk + 466 * ncomps + c);
        const auto *gk_467 = buffer.data(gk + 467 * ncomps + c);
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
        const auto *gk_496 = buffer.data(gk + 496 * ncomps + c);
        const auto *gk_497 = buffer.data(gk + 497 * ncomps + c);
        const auto *gk_498 = buffer.data(gk + 498 * ncomps + c);
        const auto *gk_499 = buffer.data(gk + 499 * ncomps + c);
        const auto *gk_500 = buffer.data(gk + 500 * ncomps + c);
        const auto *gk_501 = buffer.data(gk + 501 * ncomps + c);
        const auto *gk_502 = buffer.data(gk + 502 * ncomps + c);
        const auto *gk_503 = buffer.data(gk + 503 * ncomps + c);
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
        const auto *gk_532 = buffer.data(gk + 532 * ncomps + c);
        const auto *gk_533 = buffer.data(gk + 533 * ncomps + c);
        const auto *gk_534 = buffer.data(gk + 534 * ncomps + c);
        const auto *gk_535 = buffer.data(gk + 535 * ncomps + c);
        const auto *gk_536 = buffer.data(gk + 536 * ncomps + c);
        const auto *gk_537 = buffer.data(gk + 537 * ncomps + c);
        const auto *gk_538 = buffer.data(gk + 538 * ncomps + c);
        const auto *gk_539 = buffer.data(gk + 539 * ncomps + c);

        const auto *gl_502 = buffer.data(gl + 502 * ncomps + c);
        const auto *gl_503 = buffer.data(gl + 503 * ncomps + c);
        const auto *gl_505 = buffer.data(gl + 505 * ncomps + c);
        const auto *gl_506 = buffer.data(gl + 506 * ncomps + c);
        const auto *gl_507 = buffer.data(gl + 507 * ncomps + c);
        const auto *gl_508 = buffer.data(gl + 508 * ncomps + c);
        const auto *gl_510 = buffer.data(gl + 510 * ncomps + c);
        const auto *gl_511 = buffer.data(gl + 511 * ncomps + c);
        const auto *gl_512 = buffer.data(gl + 512 * ncomps + c);
        const auto *gl_513 = buffer.data(gl + 513 * ncomps + c);
        const auto *gl_514 = buffer.data(gl + 514 * ncomps + c);
        const auto *gl_516 = buffer.data(gl + 516 * ncomps + c);
        const auto *gl_517 = buffer.data(gl + 517 * ncomps + c);
        const auto *gl_518 = buffer.data(gl + 518 * ncomps + c);
        const auto *gl_519 = buffer.data(gl + 519 * ncomps + c);
        const auto *gl_520 = buffer.data(gl + 520 * ncomps + c);
        const auto *gl_521 = buffer.data(gl + 521 * ncomps + c);
        const auto *gl_523 = buffer.data(gl + 523 * ncomps + c);
        const auto *gl_524 = buffer.data(gl + 524 * ncomps + c);
        const auto *gl_525 = buffer.data(gl + 525 * ncomps + c);
        const auto *gl_526 = buffer.data(gl + 526 * ncomps + c);
        const auto *gl_527 = buffer.data(gl + 527 * ncomps + c);
        const auto *gl_528 = buffer.data(gl + 528 * ncomps + c);
        const auto *gl_529 = buffer.data(gl + 529 * ncomps + c);
        const auto *gl_531 = buffer.data(gl + 531 * ncomps + c);
        const auto *gl_532 = buffer.data(gl + 532 * ncomps + c);
        const auto *gl_533 = buffer.data(gl + 533 * ncomps + c);
        const auto *gl_534 = buffer.data(gl + 534 * ncomps + c);
        const auto *gl_535 = buffer.data(gl + 535 * ncomps + c);
        const auto *gl_536 = buffer.data(gl + 536 * ncomps + c);
        const auto *gl_537 = buffer.data(gl + 537 * ncomps + c);
        const auto *gl_538 = buffer.data(gl + 538 * ncomps + c);
        const auto *gl_541 = buffer.data(gl + 541 * ncomps + c);
        const auto *gl_543 = buffer.data(gl + 543 * ncomps + c);
        const auto *gl_544 = buffer.data(gl + 544 * ncomps + c);
        const auto *gl_546 = buffer.data(gl + 546 * ncomps + c);
        const auto *gl_547 = buffer.data(gl + 547 * ncomps + c);
        const auto *gl_548 = buffer.data(gl + 548 * ncomps + c);
        const auto *gl_550 = buffer.data(gl + 550 * ncomps + c);
        const auto *gl_551 = buffer.data(gl + 551 * ncomps + c);
        const auto *gl_552 = buffer.data(gl + 552 * ncomps + c);
        const auto *gl_553 = buffer.data(gl + 553 * ncomps + c);
        const auto *gl_555 = buffer.data(gl + 555 * ncomps + c);
        const auto *gl_556 = buffer.data(gl + 556 * ncomps + c);
        const auto *gl_557 = buffer.data(gl + 557 * ncomps + c);
        const auto *gl_558 = buffer.data(gl + 558 * ncomps + c);
        const auto *gl_559 = buffer.data(gl + 559 * ncomps + c);
        const auto *gl_561 = buffer.data(gl + 561 * ncomps + c);
        const auto *gl_562 = buffer.data(gl + 562 * ncomps + c);
        const auto *gl_563 = buffer.data(gl + 563 * ncomps + c);
        const auto *gl_564 = buffer.data(gl + 564 * ncomps + c);
        const auto *gl_565 = buffer.data(gl + 565 * ncomps + c);
        const auto *gl_566 = buffer.data(gl + 566 * ncomps + c);
        const auto *gl_568 = buffer.data(gl + 568 * ncomps + c);
        const auto *gl_569 = buffer.data(gl + 569 * ncomps + c);
        const auto *gl_570 = buffer.data(gl + 570 * ncomps + c);
        const auto *gl_571 = buffer.data(gl + 571 * ncomps + c);
        const auto *gl_572 = buffer.data(gl + 572 * ncomps + c);
        const auto *gl_573 = buffer.data(gl + 573 * ncomps + c);
        const auto *gl_574 = buffer.data(gl + 574 * ncomps + c);
        const auto *gl_576 = buffer.data(gl + 576 * ncomps + c);
        const auto *gl_577 = buffer.data(gl + 577 * ncomps + c);
        const auto *gl_578 = buffer.data(gl + 578 * ncomps + c);
        const auto *gl_579 = buffer.data(gl + 579 * ncomps + c);
        const auto *gl_580 = buffer.data(gl + 580 * ncomps + c);
        const auto *gl_581 = buffer.data(gl + 581 * ncomps + c);
        const auto *gl_582 = buffer.data(gl + 582 * ncomps + c);
        const auto *gl_583 = buffer.data(gl + 583 * ncomps + c);
        const auto *gl_586 = buffer.data(gl + 586 * ncomps + c);
        const auto *gl_588 = buffer.data(gl + 588 * ncomps + c);
        const auto *gl_589 = buffer.data(gl + 589 * ncomps + c);
        const auto *gl_591 = buffer.data(gl + 591 * ncomps + c);
        const auto *gl_592 = buffer.data(gl + 592 * ncomps + c);
        const auto *gl_593 = buffer.data(gl + 593 * ncomps + c);
        const auto *gl_595 = buffer.data(gl + 595 * ncomps + c);
        const auto *gl_596 = buffer.data(gl + 596 * ncomps + c);
        const auto *gl_597 = buffer.data(gl + 597 * ncomps + c);
        const auto *gl_598 = buffer.data(gl + 598 * ncomps + c);
        const auto *gl_600 = buffer.data(gl + 600 * ncomps + c);
        const auto *gl_601 = buffer.data(gl + 601 * ncomps + c);
        const auto *gl_602 = buffer.data(gl + 602 * ncomps + c);
        const auto *gl_603 = buffer.data(gl + 603 * ncomps + c);
        const auto *gl_604 = buffer.data(gl + 604 * ncomps + c);
        const auto *gl_606 = buffer.data(gl + 606 * ncomps + c);
        const auto *gl_607 = buffer.data(gl + 607 * ncomps + c);
        const auto *gl_608 = buffer.data(gl + 608 * ncomps + c);
        const auto *gl_609 = buffer.data(gl + 609 * ncomps + c);
        const auto *gl_610 = buffer.data(gl + 610 * ncomps + c);
        const auto *gl_611 = buffer.data(gl + 611 * ncomps + c);
        const auto *gl_613 = buffer.data(gl + 613 * ncomps + c);
        const auto *gl_614 = buffer.data(gl + 614 * ncomps + c);
        const auto *gl_615 = buffer.data(gl + 615 * ncomps + c);
        const auto *gl_616 = buffer.data(gl + 616 * ncomps + c);
        const auto *gl_617 = buffer.data(gl + 617 * ncomps + c);
        const auto *gl_618 = buffer.data(gl + 618 * ncomps + c);
        const auto *gl_619 = buffer.data(gl + 619 * ncomps + c);
        const auto *gl_621 = buffer.data(gl + 621 * ncomps + c);
        const auto *gl_622 = buffer.data(gl + 622 * ncomps + c);
        const auto *gl_623 = buffer.data(gl + 623 * ncomps + c);
        const auto *gl_624 = buffer.data(gl + 624 * ncomps + c);
        const auto *gl_625 = buffer.data(gl + 625 * ncomps + c);
        const auto *gl_626 = buffer.data(gl + 626 * ncomps + c);
        const auto *gl_627 = buffer.data(gl + 627 * ncomps + c);
        const auto *gl_628 = buffer.data(gl + 628 * ncomps + c);
        const auto *gl_631 = buffer.data(gl + 631 * ncomps + c);
        const auto *gl_632 = buffer.data(gl + 632 * ncomps + c);
        const auto *gl_633 = buffer.data(gl + 633 * ncomps + c);
        const auto *gl_634 = buffer.data(gl + 634 * ncomps + c);
        const auto *gl_635 = buffer.data(gl + 635 * ncomps + c);
        const auto *gl_636 = buffer.data(gl + 636 * ncomps + c);
        const auto *gl_637 = buffer.data(gl + 637 * ncomps + c);
        const auto *gl_638 = buffer.data(gl + 638 * ncomps + c);
        const auto *gl_640 = buffer.data(gl + 640 * ncomps + c);
        const auto *gl_641 = buffer.data(gl + 641 * ncomps + c);
        const auto *gl_642 = buffer.data(gl + 642 * ncomps + c);
        const auto *gl_643 = buffer.data(gl + 643 * ncomps + c);
        const auto *gl_645 = buffer.data(gl + 645 * ncomps + c);
        const auto *gl_646 = buffer.data(gl + 646 * ncomps + c);
        const auto *gl_647 = buffer.data(gl + 647 * ncomps + c);
        const auto *gl_648 = buffer.data(gl + 648 * ncomps + c);
        const auto *gl_649 = buffer.data(gl + 649 * ncomps + c);
        const auto *gl_651 = buffer.data(gl + 651 * ncomps + c);
        const auto *gl_652 = buffer.data(gl + 652 * ncomps + c);
        const auto *gl_653 = buffer.data(gl + 653 * ncomps + c);
        const auto *gl_654 = buffer.data(gl + 654 * ncomps + c);
        const auto *gl_655 = buffer.data(gl + 655 * ncomps + c);
        const auto *gl_656 = buffer.data(gl + 656 * ncomps + c);
        const auto *gl_658 = buffer.data(gl + 658 * ncomps + c);
        const auto *gl_659 = buffer.data(gl + 659 * ncomps + c);
        const auto *gl_660 = buffer.data(gl + 660 * ncomps + c);
        const auto *gl_661 = buffer.data(gl + 661 * ncomps + c);
        const auto *gl_662 = buffer.data(gl + 662 * ncomps + c);
        const auto *gl_663 = buffer.data(gl + 663 * ncomps + c);
        const auto *gl_664 = buffer.data(gl + 664 * ncomps + c);
        const auto *gl_666 = buffer.data(gl + 666 * ncomps + c);
        const auto *gl_667 = buffer.data(gl + 667 * ncomps + c);
        const auto *gl_668 = buffer.data(gl + 668 * ncomps + c);
        const auto *gl_669 = buffer.data(gl + 669 * ncomps + c);
        const auto *gl_670 = buffer.data(gl + 670 * ncomps + c);
        const auto *gl_671 = buffer.data(gl + 671 * ncomps + c);
        const auto *gl_672 = buffer.data(gl + 672 * ncomps + c);
        const auto *gl_673 = buffer.data(gl + 673 * ncomps + c);

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ab_y, gk_400, gk_401, gk_402, \
                         gk_403, gk_404, gl_502, gl_503, gl_505, gl_506, \
                         gl_507 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = -ab_y[k] * gk_400[k]
                       + gl_502[k];

            t_581[k] = -ab_y[k] * gk_401[k]
                       + gl_503[k];

            t_582[k] = -ab_y[k] * gk_402[k]
                       + gl_505[k];

            t_583[k] = -ab_y[k] * gk_403[k]
                       + gl_506[k];

            t_584[k] = -ab_y[k] * gk_404[k]
                       + gl_507[k];
        }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ab_y, gk_405, gk_406, gk_407, \
                         gk_408, gk_409, gl_508, gl_510, gl_511, gl_512, \
                         gl_513 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_585[k] = -ab_y[k] * gk_405[k]
                       + gl_508[k];

            t_586[k] = -ab_y[k] * gk_406[k]
                       + gl_510[k];

            t_587[k] = -ab_y[k] * gk_407[k]
                       + gl_511[k];

            t_588[k] = -ab_y[k] * gk_408[k]
                       + gl_512[k];

            t_589[k] = -ab_y[k] * gk_409[k]
                       + gl_513[k];
        }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ab_y, gk_410, gk_411, gk_412, \
                         gk_413, gk_414, gl_514, gl_516, gl_517, gl_518, \
                         gl_519 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_590[k] = -ab_y[k] * gk_410[k]
                       + gl_514[k];

            t_591[k] = -ab_y[k] * gk_411[k]
                       + gl_516[k];

            t_592[k] = -ab_y[k] * gk_412[k]
                       + gl_517[k];

            t_593[k] = -ab_y[k] * gk_413[k]
                       + gl_518[k];

            t_594[k] = -ab_y[k] * gk_414[k]
                       + gl_519[k];
        }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ab_y, gk_415, gk_416, gk_417, \
                         gk_418, gk_419, gl_520, gl_521, gl_523, gl_524, \
                         gl_525 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_595[k] = -ab_y[k] * gk_415[k]
                       + gl_520[k];

            t_596[k] = -ab_y[k] * gk_416[k]
                       + gl_521[k];

            t_597[k] = -ab_y[k] * gk_417[k]
                       + gl_523[k];

            t_598[k] = -ab_y[k] * gk_418[k]
                       + gl_524[k];

            t_599[k] = -ab_y[k] * gk_419[k]
                       + gl_525[k];
        }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ab_y, gk_420, gk_421, gk_422, \
                         gk_423, gk_424, gl_526, gl_527, gl_528, gl_529, \
                         gl_531 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_600[k] = -ab_y[k] * gk_420[k]
                       + gl_526[k];

            t_601[k] = -ab_y[k] * gk_421[k]
                       + gl_527[k];

            t_602[k] = -ab_y[k] * gk_422[k]
                       + gl_528[k];

            t_603[k] = -ab_y[k] * gk_423[k]
                       + gl_529[k];

            t_604[k] = -ab_y[k] * gk_424[k]
                       + gl_531[k];
        }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ab_y, gk_425, gk_426, gk_427, \
                         gk_428, gk_429, gl_532, gl_533, gl_534, gl_535, \
                         gl_536 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_605[k] = -ab_y[k] * gk_425[k]
                       + gl_532[k];

            t_606[k] = -ab_y[k] * gk_426[k]
                       + gl_533[k];

            t_607[k] = -ab_y[k] * gk_427[k]
                       + gl_534[k];

            t_608[k] = -ab_y[k] * gk_428[k]
                       + gl_535[k];

            t_609[k] = -ab_y[k] * gk_429[k]
                       + gl_536[k];
        }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ab_y, gk_430, gk_431, gk_432, \
                         gk_433, gk_434, gl_537, gl_538, gl_541, gl_543, \
                         gl_544 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_610[k] = -ab_y[k] * gk_430[k]
                       + gl_537[k];

            t_611[k] = -ab_y[k] * gk_431[k]
                       + gl_538[k];

            t_612[k] = -ab_y[k] * gk_432[k]
                       + gl_541[k];

            t_613[k] = -ab_y[k] * gk_433[k]
                       + gl_543[k];

            t_614[k] = -ab_y[k] * gk_434[k]
                       + gl_544[k];
        }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ab_y, gk_435, gk_436, gk_437, \
                         gk_438, gk_439, gl_546, gl_547, gl_548, gl_550, \
                         gl_551 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_615[k] = -ab_y[k] * gk_435[k]
                       + gl_546[k];

            t_616[k] = -ab_y[k] * gk_436[k]
                       + gl_547[k];

            t_617[k] = -ab_y[k] * gk_437[k]
                       + gl_548[k];

            t_618[k] = -ab_y[k] * gk_438[k]
                       + gl_550[k];

            t_619[k] = -ab_y[k] * gk_439[k]
                       + gl_551[k];
        }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ab_y, gk_440, gk_441, gk_442, \
                         gk_443, gk_444, gl_552, gl_553, gl_555, gl_556, \
                         gl_557 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_620[k] = -ab_y[k] * gk_440[k]
                       + gl_552[k];

            t_621[k] = -ab_y[k] * gk_441[k]
                       + gl_553[k];

            t_622[k] = -ab_y[k] * gk_442[k]
                       + gl_555[k];

            t_623[k] = -ab_y[k] * gk_443[k]
                       + gl_556[k];

            t_624[k] = -ab_y[k] * gk_444[k]
                       + gl_557[k];
        }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ab_y, gk_445, gk_446, gk_447, \
                         gk_448, gk_449, gl_558, gl_559, gl_561, gl_562, \
                         gl_563 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_625[k] = -ab_y[k] * gk_445[k]
                       + gl_558[k];

            t_626[k] = -ab_y[k] * gk_446[k]
                       + gl_559[k];

            t_627[k] = -ab_y[k] * gk_447[k]
                       + gl_561[k];

            t_628[k] = -ab_y[k] * gk_448[k]
                       + gl_562[k];

            t_629[k] = -ab_y[k] * gk_449[k]
                       + gl_563[k];
        }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ab_y, gk_450, gk_451, gk_452, \
                         gk_453, gk_454, gl_564, gl_565, gl_566, gl_568, \
                         gl_569 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_630[k] = -ab_y[k] * gk_450[k]
                       + gl_564[k];

            t_631[k] = -ab_y[k] * gk_451[k]
                       + gl_565[k];

            t_632[k] = -ab_y[k] * gk_452[k]
                       + gl_566[k];

            t_633[k] = -ab_y[k] * gk_453[k]
                       + gl_568[k];

            t_634[k] = -ab_y[k] * gk_454[k]
                       + gl_569[k];
        }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ab_y, gk_455, gk_456, gk_457, \
                         gk_458, gk_459, gl_570, gl_571, gl_572, gl_573, \
                         gl_574 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_635[k] = -ab_y[k] * gk_455[k]
                       + gl_570[k];

            t_636[k] = -ab_y[k] * gk_456[k]
                       + gl_571[k];

            t_637[k] = -ab_y[k] * gk_457[k]
                       + gl_572[k];

            t_638[k] = -ab_y[k] * gk_458[k]
                       + gl_573[k];

            t_639[k] = -ab_y[k] * gk_459[k]
                       + gl_574[k];
        }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ab_y, gk_460, gk_461, gk_462, \
                         gk_463, gk_464, gl_576, gl_577, gl_578, gl_579, \
                         gl_580 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_640[k] = -ab_y[k] * gk_460[k]
                       + gl_576[k];

            t_641[k] = -ab_y[k] * gk_461[k]
                       + gl_577[k];

            t_642[k] = -ab_y[k] * gk_462[k]
                       + gl_578[k];

            t_643[k] = -ab_y[k] * gk_463[k]
                       + gl_579[k];

            t_644[k] = -ab_y[k] * gk_464[k]
                       + gl_580[k];
        }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ab_y, gk_465, gk_466, gk_467, \
                         gk_468, gk_469, gl_581, gl_582, gl_583, gl_586, \
                         gl_588 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_645[k] = -ab_y[k] * gk_465[k]
                       + gl_581[k];

            t_646[k] = -ab_y[k] * gk_466[k]
                       + gl_582[k];

            t_647[k] = -ab_y[k] * gk_467[k]
                       + gl_583[k];

            t_648[k] = -ab_y[k] * gk_468[k]
                       + gl_586[k];

            t_649[k] = -ab_y[k] * gk_469[k]
                       + gl_588[k];
        }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ab_y, gk_470, gk_471, gk_472, \
                         gk_473, gk_474, gl_589, gl_591, gl_592, gl_593, \
                         gl_595 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_650[k] = -ab_y[k] * gk_470[k]
                       + gl_589[k];

            t_651[k] = -ab_y[k] * gk_471[k]
                       + gl_591[k];

            t_652[k] = -ab_y[k] * gk_472[k]
                       + gl_592[k];

            t_653[k] = -ab_y[k] * gk_473[k]
                       + gl_593[k];

            t_654[k] = -ab_y[k] * gk_474[k]
                       + gl_595[k];
        }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ab_y, gk_475, gk_476, gk_477, \
                         gk_478, gk_479, gl_596, gl_597, gl_598, gl_600, \
                         gl_601 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_655[k] = -ab_y[k] * gk_475[k]
                       + gl_596[k];

            t_656[k] = -ab_y[k] * gk_476[k]
                       + gl_597[k];

            t_657[k] = -ab_y[k] * gk_477[k]
                       + gl_598[k];

            t_658[k] = -ab_y[k] * gk_478[k]
                       + gl_600[k];

            t_659[k] = -ab_y[k] * gk_479[k]
                       + gl_601[k];
        }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ab_y, gk_480, gk_481, gk_482, \
                         gk_483, gk_484, gl_602, gl_603, gl_604, gl_606, \
                         gl_607 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_660[k] = -ab_y[k] * gk_480[k]
                       + gl_602[k];

            t_661[k] = -ab_y[k] * gk_481[k]
                       + gl_603[k];

            t_662[k] = -ab_y[k] * gk_482[k]
                       + gl_604[k];

            t_663[k] = -ab_y[k] * gk_483[k]
                       + gl_606[k];

            t_664[k] = -ab_y[k] * gk_484[k]
                       + gl_607[k];
        }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ab_y, gk_485, gk_486, gk_487, \
                         gk_488, gk_489, gl_608, gl_609, gl_610, gl_611, \
                         gl_613 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_665[k] = -ab_y[k] * gk_485[k]
                       + gl_608[k];

            t_666[k] = -ab_y[k] * gk_486[k]
                       + gl_609[k];

            t_667[k] = -ab_y[k] * gk_487[k]
                       + gl_610[k];

            t_668[k] = -ab_y[k] * gk_488[k]
                       + gl_611[k];

            t_669[k] = -ab_y[k] * gk_489[k]
                       + gl_613[k];
        }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ab_y, gk_490, gk_491, gk_492, \
                         gk_493, gk_494, gl_614, gl_615, gl_616, gl_617, \
                         gl_618 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_670[k] = -ab_y[k] * gk_490[k]
                       + gl_614[k];

            t_671[k] = -ab_y[k] * gk_491[k]
                       + gl_615[k];

            t_672[k] = -ab_y[k] * gk_492[k]
                       + gl_616[k];

            t_673[k] = -ab_y[k] * gk_493[k]
                       + gl_617[k];

            t_674[k] = -ab_y[k] * gk_494[k]
                       + gl_618[k];
        }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ab_y, gk_495, gk_496, gk_497, \
                         gk_498, gk_499, gl_619, gl_621, gl_622, gl_623, \
                         gl_624 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_675[k] = -ab_y[k] * gk_495[k]
                       + gl_619[k];

            t_676[k] = -ab_y[k] * gk_496[k]
                       + gl_621[k];

            t_677[k] = -ab_y[k] * gk_497[k]
                       + gl_622[k];

            t_678[k] = -ab_y[k] * gk_498[k]
                       + gl_623[k];

            t_679[k] = -ab_y[k] * gk_499[k]
                       + gl_624[k];
        }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ab_y, gk_500, gk_501, gk_502, \
                         gk_503, gk_504, gl_625, gl_626, gl_627, gl_628, \
                         gl_631 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_680[k] = -ab_y[k] * gk_500[k]
                       + gl_625[k];

            t_681[k] = -ab_y[k] * gk_501[k]
                       + gl_626[k];

            t_682[k] = -ab_y[k] * gk_502[k]
                       + gl_627[k];

            t_683[k] = -ab_y[k] * gk_503[k]
                       + gl_628[k];

            t_684[k] = -ab_y[k] * gk_504[k]
                       + gl_631[k];
        }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ab_y, gk_505, gk_506, gk_507, \
                         gk_508, gk_509, gl_633, gl_634, gl_636, gl_637, \
                         gl_638 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_685[k] = -ab_y[k] * gk_505[k]
                       + gl_633[k];

            t_686[k] = -ab_y[k] * gk_506[k]
                       + gl_634[k];

            t_687[k] = -ab_y[k] * gk_507[k]
                       + gl_636[k];

            t_688[k] = -ab_y[k] * gk_508[k]
                       + gl_637[k];

            t_689[k] = -ab_y[k] * gk_509[k]
                       + gl_638[k];
        }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ab_y, gk_510, gk_511, gk_512, \
                         gk_513, gk_514, gl_640, gl_641, gl_642, gl_643, \
                         gl_645 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_690[k] = -ab_y[k] * gk_510[k]
                       + gl_640[k];

            t_691[k] = -ab_y[k] * gk_511[k]
                       + gl_641[k];

            t_692[k] = -ab_y[k] * gk_512[k]
                       + gl_642[k];

            t_693[k] = -ab_y[k] * gk_513[k]
                       + gl_643[k];

            t_694[k] = -ab_y[k] * gk_514[k]
                       + gl_645[k];
        }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ab_y, gk_515, gk_516, gk_517, \
                         gk_518, gk_519, gl_646, gl_647, gl_648, gl_649, \
                         gl_651 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_695[k] = -ab_y[k] * gk_515[k]
                       + gl_646[k];

            t_696[k] = -ab_y[k] * gk_516[k]
                       + gl_647[k];

            t_697[k] = -ab_y[k] * gk_517[k]
                       + gl_648[k];

            t_698[k] = -ab_y[k] * gk_518[k]
                       + gl_649[k];

            t_699[k] = -ab_y[k] * gk_519[k]
                       + gl_651[k];
        }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ab_y, gk_520, gk_521, gk_522, \
                         gk_523, gk_524, gl_652, gl_653, gl_654, gl_655, \
                         gl_656 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_700[k] = -ab_y[k] * gk_520[k]
                       + gl_652[k];

            t_701[k] = -ab_y[k] * gk_521[k]
                       + gl_653[k];

            t_702[k] = -ab_y[k] * gk_522[k]
                       + gl_654[k];

            t_703[k] = -ab_y[k] * gk_523[k]
                       + gl_655[k];

            t_704[k] = -ab_y[k] * gk_524[k]
                       + gl_656[k];
        }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ab_y, gk_525, gk_526, gk_527, \
                         gk_528, gk_529, gl_658, gl_659, gl_660, gl_661, \
                         gl_662 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_705[k] = -ab_y[k] * gk_525[k]
                       + gl_658[k];

            t_706[k] = -ab_y[k] * gk_526[k]
                       + gl_659[k];

            t_707[k] = -ab_y[k] * gk_527[k]
                       + gl_660[k];

            t_708[k] = -ab_y[k] * gk_528[k]
                       + gl_661[k];

            t_709[k] = -ab_y[k] * gk_529[k]
                       + gl_662[k];
        }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ab_y, gk_530, gk_531, gk_532, \
                         gk_533, gk_534, gl_663, gl_664, gl_666, gl_667, \
                         gl_668 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_710[k] = -ab_y[k] * gk_530[k]
                       + gl_663[k];

            t_711[k] = -ab_y[k] * gk_531[k]
                       + gl_664[k];

            t_712[k] = -ab_y[k] * gk_532[k]
                       + gl_666[k];

            t_713[k] = -ab_y[k] * gk_533[k]
                       + gl_667[k];

            t_714[k] = -ab_y[k] * gk_534[k]
                       + gl_668[k];
        }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ab_y, gk_535, gk_536, gk_537, \
                         gk_538, gk_539, gl_669, gl_670, gl_671, gl_672, \
                         gl_673 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_715[k] = -ab_y[k] * gk_535[k]
                       + gl_669[k];

            t_716[k] = -ab_y[k] * gk_536[k]
                       + gl_670[k];

            t_717[k] = -ab_y[k] * gk_537[k]
                       + gl_671[k];

            t_718[k] = -ab_y[k] * gk_538[k]
                       + gl_672[k];

            t_719[k] = -ab_y[k] * gk_539[k]
                       + gl_673[k];
        }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ab_z, gk_504, gk_505, gk_506, \
                         gk_507, gk_508, gl_632, gl_634, gl_635, gl_637, \
                         gl_638 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_720[k] = -ab_z[k] * gk_504[k]
                       + gl_632[k];

            t_721[k] = -ab_z[k] * gk_505[k]
                       + gl_634[k];

            t_722[k] = -ab_z[k] * gk_506[k]
                       + gl_635[k];

            t_723[k] = -ab_z[k] * gk_507[k]
                       + gl_637[k];

            t_724[k] = -ab_z[k] * gk_508[k]
                       + gl_638[k];
        }
    }
}

static auto
compute_hrr_hk_piece5(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gk, const size_t gl, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_z = coordinates.data(8);

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
        const auto *gk_532 = buffer.data(gk + 532 * ncomps + c);
        const auto *gk_533 = buffer.data(gk + 533 * ncomps + c);
        const auto *gk_534 = buffer.data(gk + 534 * ncomps + c);
        const auto *gk_535 = buffer.data(gk + 535 * ncomps + c);
        const auto *gk_536 = buffer.data(gk + 536 * ncomps + c);
        const auto *gk_537 = buffer.data(gk + 537 * ncomps + c);
        const auto *gk_538 = buffer.data(gk + 538 * ncomps + c);
        const auto *gk_539 = buffer.data(gk + 539 * ncomps + c);

        const auto *gl_639 = buffer.data(gl + 639 * ncomps + c);
        const auto *gl_641 = buffer.data(gl + 641 * ncomps + c);
        const auto *gl_642 = buffer.data(gl + 642 * ncomps + c);
        const auto *gl_643 = buffer.data(gl + 643 * ncomps + c);
        const auto *gl_644 = buffer.data(gl + 644 * ncomps + c);
        const auto *gl_646 = buffer.data(gl + 646 * ncomps + c);
        const auto *gl_647 = buffer.data(gl + 647 * ncomps + c);
        const auto *gl_648 = buffer.data(gl + 648 * ncomps + c);
        const auto *gl_649 = buffer.data(gl + 649 * ncomps + c);
        const auto *gl_650 = buffer.data(gl + 650 * ncomps + c);
        const auto *gl_652 = buffer.data(gl + 652 * ncomps + c);
        const auto *gl_653 = buffer.data(gl + 653 * ncomps + c);
        const auto *gl_654 = buffer.data(gl + 654 * ncomps + c);
        const auto *gl_655 = buffer.data(gl + 655 * ncomps + c);
        const auto *gl_656 = buffer.data(gl + 656 * ncomps + c);
        const auto *gl_657 = buffer.data(gl + 657 * ncomps + c);
        const auto *gl_659 = buffer.data(gl + 659 * ncomps + c);
        const auto *gl_660 = buffer.data(gl + 660 * ncomps + c);
        const auto *gl_661 = buffer.data(gl + 661 * ncomps + c);
        const auto *gl_662 = buffer.data(gl + 662 * ncomps + c);
        const auto *gl_663 = buffer.data(gl + 663 * ncomps + c);
        const auto *gl_664 = buffer.data(gl + 664 * ncomps + c);
        const auto *gl_665 = buffer.data(gl + 665 * ncomps + c);
        const auto *gl_667 = buffer.data(gl + 667 * ncomps + c);
        const auto *gl_668 = buffer.data(gl + 668 * ncomps + c);
        const auto *gl_669 = buffer.data(gl + 669 * ncomps + c);
        const auto *gl_670 = buffer.data(gl + 670 * ncomps + c);
        const auto *gl_671 = buffer.data(gl + 671 * ncomps + c);
        const auto *gl_672 = buffer.data(gl + 672 * ncomps + c);
        const auto *gl_673 = buffer.data(gl + 673 * ncomps + c);
        const auto *gl_674 = buffer.data(gl + 674 * ncomps + c);

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ab_z, gk_509, gk_510, gk_511, \
                         gk_512, gk_513, gl_639, gl_641, gl_642, gl_643, \
                         gl_644 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_725[k] = -ab_z[k] * gk_509[k]
                       + gl_639[k];

            t_726[k] = -ab_z[k] * gk_510[k]
                       + gl_641[k];

            t_727[k] = -ab_z[k] * gk_511[k]
                       + gl_642[k];

            t_728[k] = -ab_z[k] * gk_512[k]
                       + gl_643[k];

            t_729[k] = -ab_z[k] * gk_513[k]
                       + gl_644[k];
        }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ab_z, gk_514, gk_515, gk_516, \
                         gk_517, gk_518, gl_646, gl_647, gl_648, gl_649, \
                         gl_650 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_730[k] = -ab_z[k] * gk_514[k]
                       + gl_646[k];

            t_731[k] = -ab_z[k] * gk_515[k]
                       + gl_647[k];

            t_732[k] = -ab_z[k] * gk_516[k]
                       + gl_648[k];

            t_733[k] = -ab_z[k] * gk_517[k]
                       + gl_649[k];

            t_734[k] = -ab_z[k] * gk_518[k]
                       + gl_650[k];
        }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ab_z, gk_519, gk_520, gk_521, \
                         gk_522, gk_523, gl_652, gl_653, gl_654, gl_655, \
                         gl_656 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_735[k] = -ab_z[k] * gk_519[k]
                       + gl_652[k];

            t_736[k] = -ab_z[k] * gk_520[k]
                       + gl_653[k];

            t_737[k] = -ab_z[k] * gk_521[k]
                       + gl_654[k];

            t_738[k] = -ab_z[k] * gk_522[k]
                       + gl_655[k];

            t_739[k] = -ab_z[k] * gk_523[k]
                       + gl_656[k];
        }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ab_z, gk_524, gk_525, gk_526, \
                         gk_527, gk_528, gl_657, gl_659, gl_660, gl_661, \
                         gl_662 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_740[k] = -ab_z[k] * gk_524[k]
                       + gl_657[k];

            t_741[k] = -ab_z[k] * gk_525[k]
                       + gl_659[k];

            t_742[k] = -ab_z[k] * gk_526[k]
                       + gl_660[k];

            t_743[k] = -ab_z[k] * gk_527[k]
                       + gl_661[k];

            t_744[k] = -ab_z[k] * gk_528[k]
                       + gl_662[k];
        }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ab_z, gk_529, gk_530, gk_531, \
                         gk_532, gk_533, gl_663, gl_664, gl_665, gl_667, \
                         gl_668 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_745[k] = -ab_z[k] * gk_529[k]
                       + gl_663[k];

            t_746[k] = -ab_z[k] * gk_530[k]
                       + gl_664[k];

            t_747[k] = -ab_z[k] * gk_531[k]
                       + gl_665[k];

            t_748[k] = -ab_z[k] * gk_532[k]
                       + gl_667[k];

            t_749[k] = -ab_z[k] * gk_533[k]
                       + gl_668[k];
        }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ab_z, gk_534, gk_535, gk_536, \
                         gk_537, gk_538, gl_669, gl_670, gl_671, gl_672, \
                         gl_673 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_750[k] = -ab_z[k] * gk_534[k]
                       + gl_669[k];

            t_751[k] = -ab_z[k] * gk_535[k]
                       + gl_670[k];

            t_752[k] = -ab_z[k] * gk_536[k]
                       + gl_671[k];

            t_753[k] = -ab_z[k] * gk_537[k]
                       + gl_672[k];

            t_754[k] = -ab_z[k] * gk_538[k]
                       + gl_673[k];
        }

#pragma omp simd aligned(t_755, ab_z, gk_539, gl_674 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_755[k] = -ab_z[k] * gk_539[k]
                       + gl_674[k];
        }
    }
}

auto
compute_hrr_hk(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gk, const size_t gl, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_hk_piece0(buffer, coordinates, target, gk, gl, ncomps, nmax);

    compute_hrr_hk_piece1(buffer, coordinates, target, gk, gl, ncomps, nmax);

    compute_hrr_hk_piece2(buffer, coordinates, target, gk, gl, ncomps, nmax);

    compute_hrr_hk_piece3(buffer, coordinates, target, gk, gl, ncomps, nmax);

    compute_hrr_hk_piece4(buffer, coordinates, target, gk, gl, ncomps, nmax);

    compute_hrr_hk_piece5(buffer, coordinates, target, gk, gl, ncomps, nmax);
}

}  // namespace simdtrf
