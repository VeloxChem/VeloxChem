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


#include "SimdTransferGK.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_gk_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fk, const size_t fl, const size_t ncomps,
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

        const auto *fk_0 = buffer.data(fk + 0 * ncomps + c);
        const auto *fk_1 = buffer.data(fk + 1 * ncomps + c);
        const auto *fk_2 = buffer.data(fk + 2 * ncomps + c);
        const auto *fk_3 = buffer.data(fk + 3 * ncomps + c);
        const auto *fk_4 = buffer.data(fk + 4 * ncomps + c);
        const auto *fk_5 = buffer.data(fk + 5 * ncomps + c);
        const auto *fk_6 = buffer.data(fk + 6 * ncomps + c);
        const auto *fk_7 = buffer.data(fk + 7 * ncomps + c);
        const auto *fk_8 = buffer.data(fk + 8 * ncomps + c);
        const auto *fk_9 = buffer.data(fk + 9 * ncomps + c);
        const auto *fk_10 = buffer.data(fk + 10 * ncomps + c);
        const auto *fk_11 = buffer.data(fk + 11 * ncomps + c);
        const auto *fk_12 = buffer.data(fk + 12 * ncomps + c);
        const auto *fk_13 = buffer.data(fk + 13 * ncomps + c);
        const auto *fk_14 = buffer.data(fk + 14 * ncomps + c);
        const auto *fk_15 = buffer.data(fk + 15 * ncomps + c);
        const auto *fk_16 = buffer.data(fk + 16 * ncomps + c);
        const auto *fk_17 = buffer.data(fk + 17 * ncomps + c);
        const auto *fk_18 = buffer.data(fk + 18 * ncomps + c);
        const auto *fk_19 = buffer.data(fk + 19 * ncomps + c);
        const auto *fk_20 = buffer.data(fk + 20 * ncomps + c);
        const auto *fk_21 = buffer.data(fk + 21 * ncomps + c);
        const auto *fk_22 = buffer.data(fk + 22 * ncomps + c);
        const auto *fk_23 = buffer.data(fk + 23 * ncomps + c);
        const auto *fk_24 = buffer.data(fk + 24 * ncomps + c);
        const auto *fk_25 = buffer.data(fk + 25 * ncomps + c);
        const auto *fk_26 = buffer.data(fk + 26 * ncomps + c);
        const auto *fk_27 = buffer.data(fk + 27 * ncomps + c);
        const auto *fk_28 = buffer.data(fk + 28 * ncomps + c);
        const auto *fk_29 = buffer.data(fk + 29 * ncomps + c);
        const auto *fk_30 = buffer.data(fk + 30 * ncomps + c);
        const auto *fk_31 = buffer.data(fk + 31 * ncomps + c);
        const auto *fk_32 = buffer.data(fk + 32 * ncomps + c);
        const auto *fk_33 = buffer.data(fk + 33 * ncomps + c);
        const auto *fk_34 = buffer.data(fk + 34 * ncomps + c);
        const auto *fk_35 = buffer.data(fk + 35 * ncomps + c);
        const auto *fk_36 = buffer.data(fk + 36 * ncomps + c);
        const auto *fk_37 = buffer.data(fk + 37 * ncomps + c);
        const auto *fk_38 = buffer.data(fk + 38 * ncomps + c);
        const auto *fk_39 = buffer.data(fk + 39 * ncomps + c);
        const auto *fk_40 = buffer.data(fk + 40 * ncomps + c);
        const auto *fk_41 = buffer.data(fk + 41 * ncomps + c);
        const auto *fk_42 = buffer.data(fk + 42 * ncomps + c);
        const auto *fk_43 = buffer.data(fk + 43 * ncomps + c);
        const auto *fk_44 = buffer.data(fk + 44 * ncomps + c);
        const auto *fk_45 = buffer.data(fk + 45 * ncomps + c);
        const auto *fk_46 = buffer.data(fk + 46 * ncomps + c);
        const auto *fk_47 = buffer.data(fk + 47 * ncomps + c);
        const auto *fk_48 = buffer.data(fk + 48 * ncomps + c);
        const auto *fk_49 = buffer.data(fk + 49 * ncomps + c);
        const auto *fk_50 = buffer.data(fk + 50 * ncomps + c);
        const auto *fk_51 = buffer.data(fk + 51 * ncomps + c);
        const auto *fk_52 = buffer.data(fk + 52 * ncomps + c);
        const auto *fk_53 = buffer.data(fk + 53 * ncomps + c);
        const auto *fk_54 = buffer.data(fk + 54 * ncomps + c);
        const auto *fk_55 = buffer.data(fk + 55 * ncomps + c);
        const auto *fk_56 = buffer.data(fk + 56 * ncomps + c);
        const auto *fk_57 = buffer.data(fk + 57 * ncomps + c);
        const auto *fk_58 = buffer.data(fk + 58 * ncomps + c);
        const auto *fk_59 = buffer.data(fk + 59 * ncomps + c);
        const auto *fk_60 = buffer.data(fk + 60 * ncomps + c);
        const auto *fk_61 = buffer.data(fk + 61 * ncomps + c);
        const auto *fk_62 = buffer.data(fk + 62 * ncomps + c);
        const auto *fk_63 = buffer.data(fk + 63 * ncomps + c);
        const auto *fk_64 = buffer.data(fk + 64 * ncomps + c);
        const auto *fk_65 = buffer.data(fk + 65 * ncomps + c);
        const auto *fk_66 = buffer.data(fk + 66 * ncomps + c);
        const auto *fk_67 = buffer.data(fk + 67 * ncomps + c);
        const auto *fk_68 = buffer.data(fk + 68 * ncomps + c);
        const auto *fk_69 = buffer.data(fk + 69 * ncomps + c);
        const auto *fk_70 = buffer.data(fk + 70 * ncomps + c);
        const auto *fk_71 = buffer.data(fk + 71 * ncomps + c);
        const auto *fk_72 = buffer.data(fk + 72 * ncomps + c);
        const auto *fk_73 = buffer.data(fk + 73 * ncomps + c);
        const auto *fk_74 = buffer.data(fk + 74 * ncomps + c);
        const auto *fk_75 = buffer.data(fk + 75 * ncomps + c);
        const auto *fk_76 = buffer.data(fk + 76 * ncomps + c);
        const auto *fk_77 = buffer.data(fk + 77 * ncomps + c);
        const auto *fk_78 = buffer.data(fk + 78 * ncomps + c);
        const auto *fk_79 = buffer.data(fk + 79 * ncomps + c);
        const auto *fk_80 = buffer.data(fk + 80 * ncomps + c);
        const auto *fk_81 = buffer.data(fk + 81 * ncomps + c);
        const auto *fk_82 = buffer.data(fk + 82 * ncomps + c);
        const auto *fk_83 = buffer.data(fk + 83 * ncomps + c);
        const auto *fk_84 = buffer.data(fk + 84 * ncomps + c);
        const auto *fk_85 = buffer.data(fk + 85 * ncomps + c);
        const auto *fk_86 = buffer.data(fk + 86 * ncomps + c);
        const auto *fk_87 = buffer.data(fk + 87 * ncomps + c);
        const auto *fk_88 = buffer.data(fk + 88 * ncomps + c);
        const auto *fk_89 = buffer.data(fk + 89 * ncomps + c);
        const auto *fk_90 = buffer.data(fk + 90 * ncomps + c);
        const auto *fk_91 = buffer.data(fk + 91 * ncomps + c);
        const auto *fk_92 = buffer.data(fk + 92 * ncomps + c);
        const auto *fk_93 = buffer.data(fk + 93 * ncomps + c);
        const auto *fk_94 = buffer.data(fk + 94 * ncomps + c);
        const auto *fk_95 = buffer.data(fk + 95 * ncomps + c);
        const auto *fk_96 = buffer.data(fk + 96 * ncomps + c);
        const auto *fk_97 = buffer.data(fk + 97 * ncomps + c);
        const auto *fk_98 = buffer.data(fk + 98 * ncomps + c);
        const auto *fk_99 = buffer.data(fk + 99 * ncomps + c);
        const auto *fk_100 = buffer.data(fk + 100 * ncomps + c);
        const auto *fk_101 = buffer.data(fk + 101 * ncomps + c);
        const auto *fk_102 = buffer.data(fk + 102 * ncomps + c);
        const auto *fk_103 = buffer.data(fk + 103 * ncomps + c);
        const auto *fk_104 = buffer.data(fk + 104 * ncomps + c);
        const auto *fk_105 = buffer.data(fk + 105 * ncomps + c);
        const auto *fk_106 = buffer.data(fk + 106 * ncomps + c);
        const auto *fk_107 = buffer.data(fk + 107 * ncomps + c);
        const auto *fk_108 = buffer.data(fk + 108 * ncomps + c);
        const auto *fk_109 = buffer.data(fk + 109 * ncomps + c);
        const auto *fk_110 = buffer.data(fk + 110 * ncomps + c);
        const auto *fk_111 = buffer.data(fk + 111 * ncomps + c);
        const auto *fk_112 = buffer.data(fk + 112 * ncomps + c);
        const auto *fk_113 = buffer.data(fk + 113 * ncomps + c);
        const auto *fk_114 = buffer.data(fk + 114 * ncomps + c);
        const auto *fk_115 = buffer.data(fk + 115 * ncomps + c);
        const auto *fk_116 = buffer.data(fk + 116 * ncomps + c);
        const auto *fk_117 = buffer.data(fk + 117 * ncomps + c);
        const auto *fk_118 = buffer.data(fk + 118 * ncomps + c);
        const auto *fk_119 = buffer.data(fk + 119 * ncomps + c);
        const auto *fk_120 = buffer.data(fk + 120 * ncomps + c);
        const auto *fk_121 = buffer.data(fk + 121 * ncomps + c);
        const auto *fk_122 = buffer.data(fk + 122 * ncomps + c);
        const auto *fk_123 = buffer.data(fk + 123 * ncomps + c);
        const auto *fk_124 = buffer.data(fk + 124 * ncomps + c);
        const auto *fk_125 = buffer.data(fk + 125 * ncomps + c);
        const auto *fk_126 = buffer.data(fk + 126 * ncomps + c);
        const auto *fk_127 = buffer.data(fk + 127 * ncomps + c);
        const auto *fk_128 = buffer.data(fk + 128 * ncomps + c);
        const auto *fk_129 = buffer.data(fk + 129 * ncomps + c);
        const auto *fk_130 = buffer.data(fk + 130 * ncomps + c);
        const auto *fk_131 = buffer.data(fk + 131 * ncomps + c);
        const auto *fk_132 = buffer.data(fk + 132 * ncomps + c);
        const auto *fk_133 = buffer.data(fk + 133 * ncomps + c);
        const auto *fk_134 = buffer.data(fk + 134 * ncomps + c);
        const auto *fk_135 = buffer.data(fk + 135 * ncomps + c);
        const auto *fk_136 = buffer.data(fk + 136 * ncomps + c);
        const auto *fk_137 = buffer.data(fk + 137 * ncomps + c);
        const auto *fk_138 = buffer.data(fk + 138 * ncomps + c);
        const auto *fk_139 = buffer.data(fk + 139 * ncomps + c);
        const auto *fk_140 = buffer.data(fk + 140 * ncomps + c);
        const auto *fk_141 = buffer.data(fk + 141 * ncomps + c);
        const auto *fk_142 = buffer.data(fk + 142 * ncomps + c);
        const auto *fk_143 = buffer.data(fk + 143 * ncomps + c);
        const auto *fk_144 = buffer.data(fk + 144 * ncomps + c);

        const auto *fl_0 = buffer.data(fl + 0 * ncomps + c);
        const auto *fl_1 = buffer.data(fl + 1 * ncomps + c);
        const auto *fl_2 = buffer.data(fl + 2 * ncomps + c);
        const auto *fl_3 = buffer.data(fl + 3 * ncomps + c);
        const auto *fl_4 = buffer.data(fl + 4 * ncomps + c);
        const auto *fl_5 = buffer.data(fl + 5 * ncomps + c);
        const auto *fl_6 = buffer.data(fl + 6 * ncomps + c);
        const auto *fl_7 = buffer.data(fl + 7 * ncomps + c);
        const auto *fl_8 = buffer.data(fl + 8 * ncomps + c);
        const auto *fl_9 = buffer.data(fl + 9 * ncomps + c);
        const auto *fl_10 = buffer.data(fl + 10 * ncomps + c);
        const auto *fl_11 = buffer.data(fl + 11 * ncomps + c);
        const auto *fl_12 = buffer.data(fl + 12 * ncomps + c);
        const auto *fl_13 = buffer.data(fl + 13 * ncomps + c);
        const auto *fl_14 = buffer.data(fl + 14 * ncomps + c);
        const auto *fl_15 = buffer.data(fl + 15 * ncomps + c);
        const auto *fl_16 = buffer.data(fl + 16 * ncomps + c);
        const auto *fl_17 = buffer.data(fl + 17 * ncomps + c);
        const auto *fl_18 = buffer.data(fl + 18 * ncomps + c);
        const auto *fl_19 = buffer.data(fl + 19 * ncomps + c);
        const auto *fl_20 = buffer.data(fl + 20 * ncomps + c);
        const auto *fl_21 = buffer.data(fl + 21 * ncomps + c);
        const auto *fl_22 = buffer.data(fl + 22 * ncomps + c);
        const auto *fl_23 = buffer.data(fl + 23 * ncomps + c);
        const auto *fl_24 = buffer.data(fl + 24 * ncomps + c);
        const auto *fl_25 = buffer.data(fl + 25 * ncomps + c);
        const auto *fl_26 = buffer.data(fl + 26 * ncomps + c);
        const auto *fl_27 = buffer.data(fl + 27 * ncomps + c);
        const auto *fl_28 = buffer.data(fl + 28 * ncomps + c);
        const auto *fl_29 = buffer.data(fl + 29 * ncomps + c);
        const auto *fl_30 = buffer.data(fl + 30 * ncomps + c);
        const auto *fl_31 = buffer.data(fl + 31 * ncomps + c);
        const auto *fl_32 = buffer.data(fl + 32 * ncomps + c);
        const auto *fl_33 = buffer.data(fl + 33 * ncomps + c);
        const auto *fl_34 = buffer.data(fl + 34 * ncomps + c);
        const auto *fl_35 = buffer.data(fl + 35 * ncomps + c);
        const auto *fl_45 = buffer.data(fl + 45 * ncomps + c);
        const auto *fl_46 = buffer.data(fl + 46 * ncomps + c);
        const auto *fl_47 = buffer.data(fl + 47 * ncomps + c);
        const auto *fl_48 = buffer.data(fl + 48 * ncomps + c);
        const auto *fl_49 = buffer.data(fl + 49 * ncomps + c);
        const auto *fl_50 = buffer.data(fl + 50 * ncomps + c);
        const auto *fl_51 = buffer.data(fl + 51 * ncomps + c);
        const auto *fl_52 = buffer.data(fl + 52 * ncomps + c);
        const auto *fl_53 = buffer.data(fl + 53 * ncomps + c);
        const auto *fl_54 = buffer.data(fl + 54 * ncomps + c);
        const auto *fl_55 = buffer.data(fl + 55 * ncomps + c);
        const auto *fl_56 = buffer.data(fl + 56 * ncomps + c);
        const auto *fl_57 = buffer.data(fl + 57 * ncomps + c);
        const auto *fl_58 = buffer.data(fl + 58 * ncomps + c);
        const auto *fl_59 = buffer.data(fl + 59 * ncomps + c);
        const auto *fl_60 = buffer.data(fl + 60 * ncomps + c);
        const auto *fl_61 = buffer.data(fl + 61 * ncomps + c);
        const auto *fl_62 = buffer.data(fl + 62 * ncomps + c);
        const auto *fl_63 = buffer.data(fl + 63 * ncomps + c);
        const auto *fl_64 = buffer.data(fl + 64 * ncomps + c);
        const auto *fl_65 = buffer.data(fl + 65 * ncomps + c);
        const auto *fl_66 = buffer.data(fl + 66 * ncomps + c);
        const auto *fl_67 = buffer.data(fl + 67 * ncomps + c);
        const auto *fl_68 = buffer.data(fl + 68 * ncomps + c);
        const auto *fl_69 = buffer.data(fl + 69 * ncomps + c);
        const auto *fl_70 = buffer.data(fl + 70 * ncomps + c);
        const auto *fl_71 = buffer.data(fl + 71 * ncomps + c);
        const auto *fl_72 = buffer.data(fl + 72 * ncomps + c);
        const auto *fl_73 = buffer.data(fl + 73 * ncomps + c);
        const auto *fl_74 = buffer.data(fl + 74 * ncomps + c);
        const auto *fl_75 = buffer.data(fl + 75 * ncomps + c);
        const auto *fl_76 = buffer.data(fl + 76 * ncomps + c);
        const auto *fl_77 = buffer.data(fl + 77 * ncomps + c);
        const auto *fl_78 = buffer.data(fl + 78 * ncomps + c);
        const auto *fl_79 = buffer.data(fl + 79 * ncomps + c);
        const auto *fl_80 = buffer.data(fl + 80 * ncomps + c);
        const auto *fl_90 = buffer.data(fl + 90 * ncomps + c);
        const auto *fl_91 = buffer.data(fl + 91 * ncomps + c);
        const auto *fl_92 = buffer.data(fl + 92 * ncomps + c);
        const auto *fl_93 = buffer.data(fl + 93 * ncomps + c);
        const auto *fl_94 = buffer.data(fl + 94 * ncomps + c);
        const auto *fl_95 = buffer.data(fl + 95 * ncomps + c);
        const auto *fl_96 = buffer.data(fl + 96 * ncomps + c);
        const auto *fl_97 = buffer.data(fl + 97 * ncomps + c);
        const auto *fl_98 = buffer.data(fl + 98 * ncomps + c);
        const auto *fl_99 = buffer.data(fl + 99 * ncomps + c);
        const auto *fl_100 = buffer.data(fl + 100 * ncomps + c);
        const auto *fl_101 = buffer.data(fl + 101 * ncomps + c);
        const auto *fl_102 = buffer.data(fl + 102 * ncomps + c);
        const auto *fl_103 = buffer.data(fl + 103 * ncomps + c);
        const auto *fl_104 = buffer.data(fl + 104 * ncomps + c);
        const auto *fl_105 = buffer.data(fl + 105 * ncomps + c);
        const auto *fl_106 = buffer.data(fl + 106 * ncomps + c);
        const auto *fl_107 = buffer.data(fl + 107 * ncomps + c);
        const auto *fl_108 = buffer.data(fl + 108 * ncomps + c);
        const auto *fl_109 = buffer.data(fl + 109 * ncomps + c);
        const auto *fl_110 = buffer.data(fl + 110 * ncomps + c);
        const auto *fl_111 = buffer.data(fl + 111 * ncomps + c);
        const auto *fl_112 = buffer.data(fl + 112 * ncomps + c);
        const auto *fl_113 = buffer.data(fl + 113 * ncomps + c);
        const auto *fl_114 = buffer.data(fl + 114 * ncomps + c);
        const auto *fl_115 = buffer.data(fl + 115 * ncomps + c);
        const auto *fl_116 = buffer.data(fl + 116 * ncomps + c);
        const auto *fl_117 = buffer.data(fl + 117 * ncomps + c);
        const auto *fl_118 = buffer.data(fl + 118 * ncomps + c);
        const auto *fl_119 = buffer.data(fl + 119 * ncomps + c);
        const auto *fl_120 = buffer.data(fl + 120 * ncomps + c);
        const auto *fl_121 = buffer.data(fl + 121 * ncomps + c);
        const auto *fl_122 = buffer.data(fl + 122 * ncomps + c);
        const auto *fl_123 = buffer.data(fl + 123 * ncomps + c);
        const auto *fl_124 = buffer.data(fl + 124 * ncomps + c);
        const auto *fl_125 = buffer.data(fl + 125 * ncomps + c);
        const auto *fl_135 = buffer.data(fl + 135 * ncomps + c);
        const auto *fl_136 = buffer.data(fl + 136 * ncomps + c);
        const auto *fl_137 = buffer.data(fl + 137 * ncomps + c);
        const auto *fl_138 = buffer.data(fl + 138 * ncomps + c);
        const auto *fl_139 = buffer.data(fl + 139 * ncomps + c);
        const auto *fl_140 = buffer.data(fl + 140 * ncomps + c);
        const auto *fl_141 = buffer.data(fl + 141 * ncomps + c);
        const auto *fl_142 = buffer.data(fl + 142 * ncomps + c);
        const auto *fl_143 = buffer.data(fl + 143 * ncomps + c);
        const auto *fl_144 = buffer.data(fl + 144 * ncomps + c);
        const auto *fl_145 = buffer.data(fl + 145 * ncomps + c);
        const auto *fl_146 = buffer.data(fl + 146 * ncomps + c);
        const auto *fl_147 = buffer.data(fl + 147 * ncomps + c);
        const auto *fl_148 = buffer.data(fl + 148 * ncomps + c);
        const auto *fl_149 = buffer.data(fl + 149 * ncomps + c);
        const auto *fl_150 = buffer.data(fl + 150 * ncomps + c);
        const auto *fl_151 = buffer.data(fl + 151 * ncomps + c);
        const auto *fl_152 = buffer.data(fl + 152 * ncomps + c);
        const auto *fl_153 = buffer.data(fl + 153 * ncomps + c);
        const auto *fl_154 = buffer.data(fl + 154 * ncomps + c);
        const auto *fl_155 = buffer.data(fl + 155 * ncomps + c);
        const auto *fl_156 = buffer.data(fl + 156 * ncomps + c);
        const auto *fl_157 = buffer.data(fl + 157 * ncomps + c);
        const auto *fl_158 = buffer.data(fl + 158 * ncomps + c);
        const auto *fl_159 = buffer.data(fl + 159 * ncomps + c);
        const auto *fl_160 = buffer.data(fl + 160 * ncomps + c);
        const auto *fl_161 = buffer.data(fl + 161 * ncomps + c);
        const auto *fl_162 = buffer.data(fl + 162 * ncomps + c);
        const auto *fl_163 = buffer.data(fl + 163 * ncomps + c);
        const auto *fl_164 = buffer.data(fl + 164 * ncomps + c);
        const auto *fl_165 = buffer.data(fl + 165 * ncomps + c);
        const auto *fl_166 = buffer.data(fl + 166 * ncomps + c);
        const auto *fl_167 = buffer.data(fl + 167 * ncomps + c);
        const auto *fl_168 = buffer.data(fl + 168 * ncomps + c);
        const auto *fl_169 = buffer.data(fl + 169 * ncomps + c);
        const auto *fl_170 = buffer.data(fl + 170 * ncomps + c);
        const auto *fl_180 = buffer.data(fl + 180 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fk_0, fk_1, fk_2, fk_3, fk_4, fl_0, \
                         fl_1, fl_2, fl_3, fl_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * fk_0[k]
                     + fl_0[k];

            t_1[k] = -ab_x[k] * fk_1[k]
                     + fl_1[k];

            t_2[k] = -ab_x[k] * fk_2[k]
                     + fl_2[k];

            t_3[k] = -ab_x[k] * fk_3[k]
                     + fl_3[k];

            t_4[k] = -ab_x[k] * fk_4[k]
                     + fl_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fk_5, fk_6, fk_7, fk_8, fk_9, fl_5, \
                         fl_6, fl_7, fl_8, fl_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * fk_5[k]
                     + fl_5[k];

            t_6[k] = -ab_x[k] * fk_6[k]
                     + fl_6[k];

            t_7[k] = -ab_x[k] * fk_7[k]
                     + fl_7[k];

            t_8[k] = -ab_x[k] * fk_8[k]
                     + fl_8[k];

            t_9[k] = -ab_x[k] * fk_9[k]
                     + fl_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fk_10, fk_11, fk_12, fk_13, \
                         fk_14, fl_10, fl_11, fl_12, fl_13, fl_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * fk_10[k]
                      + fl_10[k];

            t_11[k] = -ab_x[k] * fk_11[k]
                      + fl_11[k];

            t_12[k] = -ab_x[k] * fk_12[k]
                      + fl_12[k];

            t_13[k] = -ab_x[k] * fk_13[k]
                      + fl_13[k];

            t_14[k] = -ab_x[k] * fk_14[k]
                      + fl_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fk_15, fk_16, fk_17, fk_18, \
                         fk_19, fl_15, fl_16, fl_17, fl_18, fl_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * fk_15[k]
                      + fl_15[k];

            t_16[k] = -ab_x[k] * fk_16[k]
                      + fl_16[k];

            t_17[k] = -ab_x[k] * fk_17[k]
                      + fl_17[k];

            t_18[k] = -ab_x[k] * fk_18[k]
                      + fl_18[k];

            t_19[k] = -ab_x[k] * fk_19[k]
                      + fl_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fk_20, fk_21, fk_22, fk_23, \
                         fk_24, fl_20, fl_21, fl_22, fl_23, fl_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * fk_20[k]
                      + fl_20[k];

            t_21[k] = -ab_x[k] * fk_21[k]
                      + fl_21[k];

            t_22[k] = -ab_x[k] * fk_22[k]
                      + fl_22[k];

            t_23[k] = -ab_x[k] * fk_23[k]
                      + fl_23[k];

            t_24[k] = -ab_x[k] * fk_24[k]
                      + fl_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fk_25, fk_26, fk_27, fk_28, \
                         fk_29, fl_25, fl_26, fl_27, fl_28, fl_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * fk_25[k]
                      + fl_25[k];

            t_26[k] = -ab_x[k] * fk_26[k]
                      + fl_26[k];

            t_27[k] = -ab_x[k] * fk_27[k]
                      + fl_27[k];

            t_28[k] = -ab_x[k] * fk_28[k]
                      + fl_28[k];

            t_29[k] = -ab_x[k] * fk_29[k]
                      + fl_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fk_30, fk_31, fk_32, fk_33, \
                         fk_34, fl_30, fl_31, fl_32, fl_33, fl_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * fk_30[k]
                      + fl_30[k];

            t_31[k] = -ab_x[k] * fk_31[k]
                      + fl_31[k];

            t_32[k] = -ab_x[k] * fk_32[k]
                      + fl_32[k];

            t_33[k] = -ab_x[k] * fk_33[k]
                      + fl_33[k];

            t_34[k] = -ab_x[k] * fk_34[k]
                      + fl_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, fk_35, fk_36, fk_37, fk_38, \
                         fk_39, fl_35, fl_45, fl_46, fl_47, fl_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * fk_35[k]
                      + fl_35[k];

            t_36[k] = -ab_x[k] * fk_36[k]
                      + fl_45[k];

            t_37[k] = -ab_x[k] * fk_37[k]
                      + fl_46[k];

            t_38[k] = -ab_x[k] * fk_38[k]
                      + fl_47[k];

            t_39[k] = -ab_x[k] * fk_39[k]
                      + fl_48[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fk_40, fk_41, fk_42, fk_43, \
                         fk_44, fl_49, fl_50, fl_51, fl_52, fl_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * fk_40[k]
                      + fl_49[k];

            t_41[k] = -ab_x[k] * fk_41[k]
                      + fl_50[k];

            t_42[k] = -ab_x[k] * fk_42[k]
                      + fl_51[k];

            t_43[k] = -ab_x[k] * fk_43[k]
                      + fl_52[k];

            t_44[k] = -ab_x[k] * fk_44[k]
                      + fl_53[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, fk_45, fk_46, fk_47, fk_48, \
                         fk_49, fl_54, fl_55, fl_56, fl_57, fl_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * fk_45[k]
                      + fl_54[k];

            t_46[k] = -ab_x[k] * fk_46[k]
                      + fl_55[k];

            t_47[k] = -ab_x[k] * fk_47[k]
                      + fl_56[k];

            t_48[k] = -ab_x[k] * fk_48[k]
                      + fl_57[k];

            t_49[k] = -ab_x[k] * fk_49[k]
                      + fl_58[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fk_50, fk_51, fk_52, fk_53, \
                         fk_54, fl_59, fl_60, fl_61, fl_62, fl_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * fk_50[k]
                      + fl_59[k];

            t_51[k] = -ab_x[k] * fk_51[k]
                      + fl_60[k];

            t_52[k] = -ab_x[k] * fk_52[k]
                      + fl_61[k];

            t_53[k] = -ab_x[k] * fk_53[k]
                      + fl_62[k];

            t_54[k] = -ab_x[k] * fk_54[k]
                      + fl_63[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, fk_55, fk_56, fk_57, fk_58, \
                         fk_59, fl_64, fl_65, fl_66, fl_67, fl_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * fk_55[k]
                      + fl_64[k];

            t_56[k] = -ab_x[k] * fk_56[k]
                      + fl_65[k];

            t_57[k] = -ab_x[k] * fk_57[k]
                      + fl_66[k];

            t_58[k] = -ab_x[k] * fk_58[k]
                      + fl_67[k];

            t_59[k] = -ab_x[k] * fk_59[k]
                      + fl_68[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fk_60, fk_61, fk_62, fk_63, \
                         fk_64, fl_69, fl_70, fl_71, fl_72, fl_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * fk_60[k]
                      + fl_69[k];

            t_61[k] = -ab_x[k] * fk_61[k]
                      + fl_70[k];

            t_62[k] = -ab_x[k] * fk_62[k]
                      + fl_71[k];

            t_63[k] = -ab_x[k] * fk_63[k]
                      + fl_72[k];

            t_64[k] = -ab_x[k] * fk_64[k]
                      + fl_73[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, fk_65, fk_66, fk_67, fk_68, \
                         fk_69, fl_74, fl_75, fl_76, fl_77, fl_78 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * fk_65[k]
                      + fl_74[k];

            t_66[k] = -ab_x[k] * fk_66[k]
                      + fl_75[k];

            t_67[k] = -ab_x[k] * fk_67[k]
                      + fl_76[k];

            t_68[k] = -ab_x[k] * fk_68[k]
                      + fl_77[k];

            t_69[k] = -ab_x[k] * fk_69[k]
                      + fl_78[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fk_70, fk_71, fk_72, fk_73, \
                         fk_74, fl_79, fl_80, fl_90, fl_91, fl_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * fk_70[k]
                      + fl_79[k];

            t_71[k] = -ab_x[k] * fk_71[k]
                      + fl_80[k];

            t_72[k] = -ab_x[k] * fk_72[k]
                      + fl_90[k];

            t_73[k] = -ab_x[k] * fk_73[k]
                      + fl_91[k];

            t_74[k] = -ab_x[k] * fk_74[k]
                      + fl_92[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, fk_75, fk_76, fk_77, fk_78, \
                         fk_79, fl_93, fl_94, fl_95, fl_96, fl_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * fk_75[k]
                      + fl_93[k];

            t_76[k] = -ab_x[k] * fk_76[k]
                      + fl_94[k];

            t_77[k] = -ab_x[k] * fk_77[k]
                      + fl_95[k];

            t_78[k] = -ab_x[k] * fk_78[k]
                      + fl_96[k];

            t_79[k] = -ab_x[k] * fk_79[k]
                      + fl_97[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fk_80, fk_81, fk_82, fk_83, \
                         fk_84, fl_98, fl_99, fl_100, fl_101, fl_102 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * fk_80[k]
                      + fl_98[k];

            t_81[k] = -ab_x[k] * fk_81[k]
                      + fl_99[k];

            t_82[k] = -ab_x[k] * fk_82[k]
                      + fl_100[k];

            t_83[k] = -ab_x[k] * fk_83[k]
                      + fl_101[k];

            t_84[k] = -ab_x[k] * fk_84[k]
                      + fl_102[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, fk_85, fk_86, fk_87, fk_88, \
                         fk_89, fl_103, fl_104, fl_105, fl_106, \
                         fl_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * fk_85[k]
                      + fl_103[k];

            t_86[k] = -ab_x[k] * fk_86[k]
                      + fl_104[k];

            t_87[k] = -ab_x[k] * fk_87[k]
                      + fl_105[k];

            t_88[k] = -ab_x[k] * fk_88[k]
                      + fl_106[k];

            t_89[k] = -ab_x[k] * fk_89[k]
                      + fl_107[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fk_90, fk_91, fk_92, fk_93, \
                         fk_94, fl_108, fl_109, fl_110, fl_111, \
                         fl_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * fk_90[k]
                      + fl_108[k];

            t_91[k] = -ab_x[k] * fk_91[k]
                      + fl_109[k];

            t_92[k] = -ab_x[k] * fk_92[k]
                      + fl_110[k];

            t_93[k] = -ab_x[k] * fk_93[k]
                      + fl_111[k];

            t_94[k] = -ab_x[k] * fk_94[k]
                      + fl_112[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, fk_95, fk_96, fk_97, fk_98, \
                         fk_99, fl_113, fl_114, fl_115, fl_116, \
                         fl_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * fk_95[k]
                      + fl_113[k];

            t_96[k] = -ab_x[k] * fk_96[k]
                      + fl_114[k];

            t_97[k] = -ab_x[k] * fk_97[k]
                      + fl_115[k];

            t_98[k] = -ab_x[k] * fk_98[k]
                      + fl_116[k];

            t_99[k] = -ab_x[k] * fk_99[k]
                      + fl_117[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, fk_100, fk_101, fk_102, \
                         fk_103, fk_104, fl_118, fl_119, fl_120, fl_121, \
                         fl_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * fk_100[k]
                       + fl_118[k];

            t_101[k] = -ab_x[k] * fk_101[k]
                       + fl_119[k];

            t_102[k] = -ab_x[k] * fk_102[k]
                       + fl_120[k];

            t_103[k] = -ab_x[k] * fk_103[k]
                       + fl_121[k];

            t_104[k] = -ab_x[k] * fk_104[k]
                       + fl_122[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, fk_105, fk_106, fk_107, \
                         fk_108, fk_109, fl_123, fl_124, fl_125, fl_135, \
                         fl_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * fk_105[k]
                       + fl_123[k];

            t_106[k] = -ab_x[k] * fk_106[k]
                       + fl_124[k];

            t_107[k] = -ab_x[k] * fk_107[k]
                       + fl_125[k];

            t_108[k] = -ab_x[k] * fk_108[k]
                       + fl_135[k];

            t_109[k] = -ab_x[k] * fk_109[k]
                       + fl_136[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, fk_110, fk_111, fk_112, \
                         fk_113, fk_114, fl_137, fl_138, fl_139, fl_140, \
                         fl_141 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * fk_110[k]
                       + fl_137[k];

            t_111[k] = -ab_x[k] * fk_111[k]
                       + fl_138[k];

            t_112[k] = -ab_x[k] * fk_112[k]
                       + fl_139[k];

            t_113[k] = -ab_x[k] * fk_113[k]
                       + fl_140[k];

            t_114[k] = -ab_x[k] * fk_114[k]
                       + fl_141[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, fk_115, fk_116, fk_117, \
                         fk_118, fk_119, fl_142, fl_143, fl_144, fl_145, \
                         fl_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * fk_115[k]
                       + fl_142[k];

            t_116[k] = -ab_x[k] * fk_116[k]
                       + fl_143[k];

            t_117[k] = -ab_x[k] * fk_117[k]
                       + fl_144[k];

            t_118[k] = -ab_x[k] * fk_118[k]
                       + fl_145[k];

            t_119[k] = -ab_x[k] * fk_119[k]
                       + fl_146[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, fk_120, fk_121, fk_122, \
                         fk_123, fk_124, fl_147, fl_148, fl_149, fl_150, \
                         fl_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * fk_120[k]
                       + fl_147[k];

            t_121[k] = -ab_x[k] * fk_121[k]
                       + fl_148[k];

            t_122[k] = -ab_x[k] * fk_122[k]
                       + fl_149[k];

            t_123[k] = -ab_x[k] * fk_123[k]
                       + fl_150[k];

            t_124[k] = -ab_x[k] * fk_124[k]
                       + fl_151[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, fk_125, fk_126, fk_127, \
                         fk_128, fk_129, fl_152, fl_153, fl_154, fl_155, \
                         fl_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * fk_125[k]
                       + fl_152[k];

            t_126[k] = -ab_x[k] * fk_126[k]
                       + fl_153[k];

            t_127[k] = -ab_x[k] * fk_127[k]
                       + fl_154[k];

            t_128[k] = -ab_x[k] * fk_128[k]
                       + fl_155[k];

            t_129[k] = -ab_x[k] * fk_129[k]
                       + fl_156[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, fk_130, fk_131, fk_132, \
                         fk_133, fk_134, fl_157, fl_158, fl_159, fl_160, \
                         fl_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * fk_130[k]
                       + fl_157[k];

            t_131[k] = -ab_x[k] * fk_131[k]
                       + fl_158[k];

            t_132[k] = -ab_x[k] * fk_132[k]
                       + fl_159[k];

            t_133[k] = -ab_x[k] * fk_133[k]
                       + fl_160[k];

            t_134[k] = -ab_x[k] * fk_134[k]
                       + fl_161[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, fk_135, fk_136, fk_137, \
                         fk_138, fk_139, fl_162, fl_163, fl_164, fl_165, \
                         fl_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * fk_135[k]
                       + fl_162[k];

            t_136[k] = -ab_x[k] * fk_136[k]
                       + fl_163[k];

            t_137[k] = -ab_x[k] * fk_137[k]
                       + fl_164[k];

            t_138[k] = -ab_x[k] * fk_138[k]
                       + fl_165[k];

            t_139[k] = -ab_x[k] * fk_139[k]
                       + fl_166[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, fk_140, fk_141, fk_142, \
                         fk_143, fk_144, fl_167, fl_168, fl_169, fl_170, \
                         fl_180 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * fk_140[k]
                       + fl_167[k];

            t_141[k] = -ab_x[k] * fk_141[k]
                       + fl_168[k];

            t_142[k] = -ab_x[k] * fk_142[k]
                       + fl_169[k];

            t_143[k] = -ab_x[k] * fk_143[k]
                       + fl_170[k];

            t_144[k] = -ab_x[k] * fk_144[k]
                       + fl_180[k];
        }
    }
}

static auto
compute_hrr_gk_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fk, const size_t fl, const size_t ncomps,
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

        const auto *fk_145 = buffer.data(fk + 145 * ncomps + c);
        const auto *fk_146 = buffer.data(fk + 146 * ncomps + c);
        const auto *fk_147 = buffer.data(fk + 147 * ncomps + c);
        const auto *fk_148 = buffer.data(fk + 148 * ncomps + c);
        const auto *fk_149 = buffer.data(fk + 149 * ncomps + c);
        const auto *fk_150 = buffer.data(fk + 150 * ncomps + c);
        const auto *fk_151 = buffer.data(fk + 151 * ncomps + c);
        const auto *fk_152 = buffer.data(fk + 152 * ncomps + c);
        const auto *fk_153 = buffer.data(fk + 153 * ncomps + c);
        const auto *fk_154 = buffer.data(fk + 154 * ncomps + c);
        const auto *fk_155 = buffer.data(fk + 155 * ncomps + c);
        const auto *fk_156 = buffer.data(fk + 156 * ncomps + c);
        const auto *fk_157 = buffer.data(fk + 157 * ncomps + c);
        const auto *fk_158 = buffer.data(fk + 158 * ncomps + c);
        const auto *fk_159 = buffer.data(fk + 159 * ncomps + c);
        const auto *fk_160 = buffer.data(fk + 160 * ncomps + c);
        const auto *fk_161 = buffer.data(fk + 161 * ncomps + c);
        const auto *fk_162 = buffer.data(fk + 162 * ncomps + c);
        const auto *fk_163 = buffer.data(fk + 163 * ncomps + c);
        const auto *fk_164 = buffer.data(fk + 164 * ncomps + c);
        const auto *fk_165 = buffer.data(fk + 165 * ncomps + c);
        const auto *fk_166 = buffer.data(fk + 166 * ncomps + c);
        const auto *fk_167 = buffer.data(fk + 167 * ncomps + c);
        const auto *fk_168 = buffer.data(fk + 168 * ncomps + c);
        const auto *fk_169 = buffer.data(fk + 169 * ncomps + c);
        const auto *fk_170 = buffer.data(fk + 170 * ncomps + c);
        const auto *fk_171 = buffer.data(fk + 171 * ncomps + c);
        const auto *fk_172 = buffer.data(fk + 172 * ncomps + c);
        const auto *fk_173 = buffer.data(fk + 173 * ncomps + c);
        const auto *fk_174 = buffer.data(fk + 174 * ncomps + c);
        const auto *fk_175 = buffer.data(fk + 175 * ncomps + c);
        const auto *fk_176 = buffer.data(fk + 176 * ncomps + c);
        const auto *fk_177 = buffer.data(fk + 177 * ncomps + c);
        const auto *fk_178 = buffer.data(fk + 178 * ncomps + c);
        const auto *fk_179 = buffer.data(fk + 179 * ncomps + c);
        const auto *fk_180 = buffer.data(fk + 180 * ncomps + c);
        const auto *fk_181 = buffer.data(fk + 181 * ncomps + c);
        const auto *fk_182 = buffer.data(fk + 182 * ncomps + c);
        const auto *fk_183 = buffer.data(fk + 183 * ncomps + c);
        const auto *fk_184 = buffer.data(fk + 184 * ncomps + c);
        const auto *fk_185 = buffer.data(fk + 185 * ncomps + c);
        const auto *fk_186 = buffer.data(fk + 186 * ncomps + c);
        const auto *fk_187 = buffer.data(fk + 187 * ncomps + c);
        const auto *fk_188 = buffer.data(fk + 188 * ncomps + c);
        const auto *fk_189 = buffer.data(fk + 189 * ncomps + c);
        const auto *fk_190 = buffer.data(fk + 190 * ncomps + c);
        const auto *fk_191 = buffer.data(fk + 191 * ncomps + c);
        const auto *fk_192 = buffer.data(fk + 192 * ncomps + c);
        const auto *fk_193 = buffer.data(fk + 193 * ncomps + c);
        const auto *fk_194 = buffer.data(fk + 194 * ncomps + c);
        const auto *fk_195 = buffer.data(fk + 195 * ncomps + c);
        const auto *fk_196 = buffer.data(fk + 196 * ncomps + c);
        const auto *fk_197 = buffer.data(fk + 197 * ncomps + c);
        const auto *fk_198 = buffer.data(fk + 198 * ncomps + c);
        const auto *fk_199 = buffer.data(fk + 199 * ncomps + c);
        const auto *fk_200 = buffer.data(fk + 200 * ncomps + c);
        const auto *fk_201 = buffer.data(fk + 201 * ncomps + c);
        const auto *fk_202 = buffer.data(fk + 202 * ncomps + c);
        const auto *fk_203 = buffer.data(fk + 203 * ncomps + c);
        const auto *fk_204 = buffer.data(fk + 204 * ncomps + c);
        const auto *fk_205 = buffer.data(fk + 205 * ncomps + c);
        const auto *fk_206 = buffer.data(fk + 206 * ncomps + c);
        const auto *fk_207 = buffer.data(fk + 207 * ncomps + c);
        const auto *fk_208 = buffer.data(fk + 208 * ncomps + c);
        const auto *fk_209 = buffer.data(fk + 209 * ncomps + c);
        const auto *fk_210 = buffer.data(fk + 210 * ncomps + c);
        const auto *fk_211 = buffer.data(fk + 211 * ncomps + c);
        const auto *fk_212 = buffer.data(fk + 212 * ncomps + c);
        const auto *fk_213 = buffer.data(fk + 213 * ncomps + c);
        const auto *fk_214 = buffer.data(fk + 214 * ncomps + c);
        const auto *fk_215 = buffer.data(fk + 215 * ncomps + c);
        const auto *fk_216 = buffer.data(fk + 216 * ncomps + c);
        const auto *fk_217 = buffer.data(fk + 217 * ncomps + c);
        const auto *fk_218 = buffer.data(fk + 218 * ncomps + c);
        const auto *fk_219 = buffer.data(fk + 219 * ncomps + c);
        const auto *fk_220 = buffer.data(fk + 220 * ncomps + c);
        const auto *fk_221 = buffer.data(fk + 221 * ncomps + c);
        const auto *fk_222 = buffer.data(fk + 222 * ncomps + c);
        const auto *fk_223 = buffer.data(fk + 223 * ncomps + c);
        const auto *fk_224 = buffer.data(fk + 224 * ncomps + c);
        const auto *fk_225 = buffer.data(fk + 225 * ncomps + c);
        const auto *fk_226 = buffer.data(fk + 226 * ncomps + c);
        const auto *fk_227 = buffer.data(fk + 227 * ncomps + c);
        const auto *fk_228 = buffer.data(fk + 228 * ncomps + c);
        const auto *fk_229 = buffer.data(fk + 229 * ncomps + c);
        const auto *fk_230 = buffer.data(fk + 230 * ncomps + c);
        const auto *fk_231 = buffer.data(fk + 231 * ncomps + c);
        const auto *fk_232 = buffer.data(fk + 232 * ncomps + c);
        const auto *fk_233 = buffer.data(fk + 233 * ncomps + c);
        const auto *fk_234 = buffer.data(fk + 234 * ncomps + c);
        const auto *fk_235 = buffer.data(fk + 235 * ncomps + c);
        const auto *fk_236 = buffer.data(fk + 236 * ncomps + c);
        const auto *fk_237 = buffer.data(fk + 237 * ncomps + c);
        const auto *fk_238 = buffer.data(fk + 238 * ncomps + c);
        const auto *fk_239 = buffer.data(fk + 239 * ncomps + c);
        const auto *fk_240 = buffer.data(fk + 240 * ncomps + c);
        const auto *fk_241 = buffer.data(fk + 241 * ncomps + c);
        const auto *fk_242 = buffer.data(fk + 242 * ncomps + c);
        const auto *fk_243 = buffer.data(fk + 243 * ncomps + c);
        const auto *fk_244 = buffer.data(fk + 244 * ncomps + c);
        const auto *fk_245 = buffer.data(fk + 245 * ncomps + c);
        const auto *fk_246 = buffer.data(fk + 246 * ncomps + c);
        const auto *fk_247 = buffer.data(fk + 247 * ncomps + c);
        const auto *fk_248 = buffer.data(fk + 248 * ncomps + c);
        const auto *fk_249 = buffer.data(fk + 249 * ncomps + c);
        const auto *fk_250 = buffer.data(fk + 250 * ncomps + c);
        const auto *fk_251 = buffer.data(fk + 251 * ncomps + c);
        const auto *fk_252 = buffer.data(fk + 252 * ncomps + c);
        const auto *fk_253 = buffer.data(fk + 253 * ncomps + c);
        const auto *fk_254 = buffer.data(fk + 254 * ncomps + c);
        const auto *fk_255 = buffer.data(fk + 255 * ncomps + c);
        const auto *fk_256 = buffer.data(fk + 256 * ncomps + c);
        const auto *fk_257 = buffer.data(fk + 257 * ncomps + c);
        const auto *fk_258 = buffer.data(fk + 258 * ncomps + c);
        const auto *fk_259 = buffer.data(fk + 259 * ncomps + c);
        const auto *fk_260 = buffer.data(fk + 260 * ncomps + c);
        const auto *fk_261 = buffer.data(fk + 261 * ncomps + c);
        const auto *fk_262 = buffer.data(fk + 262 * ncomps + c);
        const auto *fk_263 = buffer.data(fk + 263 * ncomps + c);
        const auto *fk_264 = buffer.data(fk + 264 * ncomps + c);
        const auto *fk_265 = buffer.data(fk + 265 * ncomps + c);
        const auto *fk_266 = buffer.data(fk + 266 * ncomps + c);
        const auto *fk_267 = buffer.data(fk + 267 * ncomps + c);
        const auto *fk_268 = buffer.data(fk + 268 * ncomps + c);
        const auto *fk_269 = buffer.data(fk + 269 * ncomps + c);
        const auto *fk_270 = buffer.data(fk + 270 * ncomps + c);
        const auto *fk_271 = buffer.data(fk + 271 * ncomps + c);
        const auto *fk_272 = buffer.data(fk + 272 * ncomps + c);
        const auto *fk_273 = buffer.data(fk + 273 * ncomps + c);
        const auto *fk_274 = buffer.data(fk + 274 * ncomps + c);
        const auto *fk_275 = buffer.data(fk + 275 * ncomps + c);
        const auto *fk_276 = buffer.data(fk + 276 * ncomps + c);
        const auto *fk_277 = buffer.data(fk + 277 * ncomps + c);
        const auto *fk_278 = buffer.data(fk + 278 * ncomps + c);
        const auto *fk_279 = buffer.data(fk + 279 * ncomps + c);
        const auto *fk_280 = buffer.data(fk + 280 * ncomps + c);
        const auto *fk_281 = buffer.data(fk + 281 * ncomps + c);
        const auto *fk_282 = buffer.data(fk + 282 * ncomps + c);
        const auto *fk_283 = buffer.data(fk + 283 * ncomps + c);
        const auto *fk_284 = buffer.data(fk + 284 * ncomps + c);
        const auto *fk_285 = buffer.data(fk + 285 * ncomps + c);
        const auto *fk_286 = buffer.data(fk + 286 * ncomps + c);
        const auto *fk_287 = buffer.data(fk + 287 * ncomps + c);
        const auto *fk_288 = buffer.data(fk + 288 * ncomps + c);
        const auto *fk_289 = buffer.data(fk + 289 * ncomps + c);

        const auto *fl_181 = buffer.data(fl + 181 * ncomps + c);
        const auto *fl_182 = buffer.data(fl + 182 * ncomps + c);
        const auto *fl_183 = buffer.data(fl + 183 * ncomps + c);
        const auto *fl_184 = buffer.data(fl + 184 * ncomps + c);
        const auto *fl_185 = buffer.data(fl + 185 * ncomps + c);
        const auto *fl_186 = buffer.data(fl + 186 * ncomps + c);
        const auto *fl_187 = buffer.data(fl + 187 * ncomps + c);
        const auto *fl_188 = buffer.data(fl + 188 * ncomps + c);
        const auto *fl_189 = buffer.data(fl + 189 * ncomps + c);
        const auto *fl_190 = buffer.data(fl + 190 * ncomps + c);
        const auto *fl_191 = buffer.data(fl + 191 * ncomps + c);
        const auto *fl_192 = buffer.data(fl + 192 * ncomps + c);
        const auto *fl_193 = buffer.data(fl + 193 * ncomps + c);
        const auto *fl_194 = buffer.data(fl + 194 * ncomps + c);
        const auto *fl_195 = buffer.data(fl + 195 * ncomps + c);
        const auto *fl_196 = buffer.data(fl + 196 * ncomps + c);
        const auto *fl_197 = buffer.data(fl + 197 * ncomps + c);
        const auto *fl_198 = buffer.data(fl + 198 * ncomps + c);
        const auto *fl_199 = buffer.data(fl + 199 * ncomps + c);
        const auto *fl_200 = buffer.data(fl + 200 * ncomps + c);
        const auto *fl_201 = buffer.data(fl + 201 * ncomps + c);
        const auto *fl_202 = buffer.data(fl + 202 * ncomps + c);
        const auto *fl_203 = buffer.data(fl + 203 * ncomps + c);
        const auto *fl_204 = buffer.data(fl + 204 * ncomps + c);
        const auto *fl_205 = buffer.data(fl + 205 * ncomps + c);
        const auto *fl_206 = buffer.data(fl + 206 * ncomps + c);
        const auto *fl_207 = buffer.data(fl + 207 * ncomps + c);
        const auto *fl_208 = buffer.data(fl + 208 * ncomps + c);
        const auto *fl_209 = buffer.data(fl + 209 * ncomps + c);
        const auto *fl_210 = buffer.data(fl + 210 * ncomps + c);
        const auto *fl_211 = buffer.data(fl + 211 * ncomps + c);
        const auto *fl_212 = buffer.data(fl + 212 * ncomps + c);
        const auto *fl_213 = buffer.data(fl + 213 * ncomps + c);
        const auto *fl_214 = buffer.data(fl + 214 * ncomps + c);
        const auto *fl_215 = buffer.data(fl + 215 * ncomps + c);
        const auto *fl_225 = buffer.data(fl + 225 * ncomps + c);
        const auto *fl_226 = buffer.data(fl + 226 * ncomps + c);
        const auto *fl_227 = buffer.data(fl + 227 * ncomps + c);
        const auto *fl_228 = buffer.data(fl + 228 * ncomps + c);
        const auto *fl_229 = buffer.data(fl + 229 * ncomps + c);
        const auto *fl_230 = buffer.data(fl + 230 * ncomps + c);
        const auto *fl_231 = buffer.data(fl + 231 * ncomps + c);
        const auto *fl_232 = buffer.data(fl + 232 * ncomps + c);
        const auto *fl_233 = buffer.data(fl + 233 * ncomps + c);
        const auto *fl_234 = buffer.data(fl + 234 * ncomps + c);
        const auto *fl_235 = buffer.data(fl + 235 * ncomps + c);
        const auto *fl_236 = buffer.data(fl + 236 * ncomps + c);
        const auto *fl_237 = buffer.data(fl + 237 * ncomps + c);
        const auto *fl_238 = buffer.data(fl + 238 * ncomps + c);
        const auto *fl_239 = buffer.data(fl + 239 * ncomps + c);
        const auto *fl_240 = buffer.data(fl + 240 * ncomps + c);
        const auto *fl_241 = buffer.data(fl + 241 * ncomps + c);
        const auto *fl_242 = buffer.data(fl + 242 * ncomps + c);
        const auto *fl_243 = buffer.data(fl + 243 * ncomps + c);
        const auto *fl_244 = buffer.data(fl + 244 * ncomps + c);
        const auto *fl_245 = buffer.data(fl + 245 * ncomps + c);
        const auto *fl_246 = buffer.data(fl + 246 * ncomps + c);
        const auto *fl_247 = buffer.data(fl + 247 * ncomps + c);
        const auto *fl_248 = buffer.data(fl + 248 * ncomps + c);
        const auto *fl_249 = buffer.data(fl + 249 * ncomps + c);
        const auto *fl_250 = buffer.data(fl + 250 * ncomps + c);
        const auto *fl_251 = buffer.data(fl + 251 * ncomps + c);
        const auto *fl_252 = buffer.data(fl + 252 * ncomps + c);
        const auto *fl_253 = buffer.data(fl + 253 * ncomps + c);
        const auto *fl_254 = buffer.data(fl + 254 * ncomps + c);
        const auto *fl_255 = buffer.data(fl + 255 * ncomps + c);
        const auto *fl_256 = buffer.data(fl + 256 * ncomps + c);
        const auto *fl_257 = buffer.data(fl + 257 * ncomps + c);
        const auto *fl_258 = buffer.data(fl + 258 * ncomps + c);
        const auto *fl_259 = buffer.data(fl + 259 * ncomps + c);
        const auto *fl_260 = buffer.data(fl + 260 * ncomps + c);
        const auto *fl_270 = buffer.data(fl + 270 * ncomps + c);
        const auto *fl_271 = buffer.data(fl + 271 * ncomps + c);
        const auto *fl_272 = buffer.data(fl + 272 * ncomps + c);
        const auto *fl_273 = buffer.data(fl + 273 * ncomps + c);
        const auto *fl_274 = buffer.data(fl + 274 * ncomps + c);
        const auto *fl_275 = buffer.data(fl + 275 * ncomps + c);
        const auto *fl_276 = buffer.data(fl + 276 * ncomps + c);
        const auto *fl_277 = buffer.data(fl + 277 * ncomps + c);
        const auto *fl_278 = buffer.data(fl + 278 * ncomps + c);
        const auto *fl_279 = buffer.data(fl + 279 * ncomps + c);
        const auto *fl_280 = buffer.data(fl + 280 * ncomps + c);
        const auto *fl_281 = buffer.data(fl + 281 * ncomps + c);
        const auto *fl_282 = buffer.data(fl + 282 * ncomps + c);
        const auto *fl_283 = buffer.data(fl + 283 * ncomps + c);
        const auto *fl_284 = buffer.data(fl + 284 * ncomps + c);
        const auto *fl_285 = buffer.data(fl + 285 * ncomps + c);
        const auto *fl_286 = buffer.data(fl + 286 * ncomps + c);
        const auto *fl_287 = buffer.data(fl + 287 * ncomps + c);
        const auto *fl_288 = buffer.data(fl + 288 * ncomps + c);
        const auto *fl_289 = buffer.data(fl + 289 * ncomps + c);
        const auto *fl_290 = buffer.data(fl + 290 * ncomps + c);
        const auto *fl_291 = buffer.data(fl + 291 * ncomps + c);
        const auto *fl_292 = buffer.data(fl + 292 * ncomps + c);
        const auto *fl_293 = buffer.data(fl + 293 * ncomps + c);
        const auto *fl_294 = buffer.data(fl + 294 * ncomps + c);
        const auto *fl_295 = buffer.data(fl + 295 * ncomps + c);
        const auto *fl_296 = buffer.data(fl + 296 * ncomps + c);
        const auto *fl_297 = buffer.data(fl + 297 * ncomps + c);
        const auto *fl_298 = buffer.data(fl + 298 * ncomps + c);
        const auto *fl_299 = buffer.data(fl + 299 * ncomps + c);
        const auto *fl_300 = buffer.data(fl + 300 * ncomps + c);
        const auto *fl_301 = buffer.data(fl + 301 * ncomps + c);
        const auto *fl_302 = buffer.data(fl + 302 * ncomps + c);
        const auto *fl_303 = buffer.data(fl + 303 * ncomps + c);
        const auto *fl_304 = buffer.data(fl + 304 * ncomps + c);
        const auto *fl_305 = buffer.data(fl + 305 * ncomps + c);
        const auto *fl_315 = buffer.data(fl + 315 * ncomps + c);
        const auto *fl_316 = buffer.data(fl + 316 * ncomps + c);
        const auto *fl_317 = buffer.data(fl + 317 * ncomps + c);
        const auto *fl_318 = buffer.data(fl + 318 * ncomps + c);
        const auto *fl_319 = buffer.data(fl + 319 * ncomps + c);
        const auto *fl_320 = buffer.data(fl + 320 * ncomps + c);
        const auto *fl_321 = buffer.data(fl + 321 * ncomps + c);
        const auto *fl_322 = buffer.data(fl + 322 * ncomps + c);
        const auto *fl_323 = buffer.data(fl + 323 * ncomps + c);
        const auto *fl_324 = buffer.data(fl + 324 * ncomps + c);
        const auto *fl_325 = buffer.data(fl + 325 * ncomps + c);
        const auto *fl_326 = buffer.data(fl + 326 * ncomps + c);
        const auto *fl_327 = buffer.data(fl + 327 * ncomps + c);
        const auto *fl_328 = buffer.data(fl + 328 * ncomps + c);
        const auto *fl_329 = buffer.data(fl + 329 * ncomps + c);
        const auto *fl_330 = buffer.data(fl + 330 * ncomps + c);
        const auto *fl_331 = buffer.data(fl + 331 * ncomps + c);
        const auto *fl_332 = buffer.data(fl + 332 * ncomps + c);
        const auto *fl_333 = buffer.data(fl + 333 * ncomps + c);
        const auto *fl_334 = buffer.data(fl + 334 * ncomps + c);
        const auto *fl_335 = buffer.data(fl + 335 * ncomps + c);
        const auto *fl_336 = buffer.data(fl + 336 * ncomps + c);
        const auto *fl_337 = buffer.data(fl + 337 * ncomps + c);
        const auto *fl_338 = buffer.data(fl + 338 * ncomps + c);
        const auto *fl_339 = buffer.data(fl + 339 * ncomps + c);
        const auto *fl_340 = buffer.data(fl + 340 * ncomps + c);
        const auto *fl_341 = buffer.data(fl + 341 * ncomps + c);
        const auto *fl_342 = buffer.data(fl + 342 * ncomps + c);
        const auto *fl_343 = buffer.data(fl + 343 * ncomps + c);
        const auto *fl_344 = buffer.data(fl + 344 * ncomps + c);
        const auto *fl_345 = buffer.data(fl + 345 * ncomps + c);
        const auto *fl_346 = buffer.data(fl + 346 * ncomps + c);
        const auto *fl_347 = buffer.data(fl + 347 * ncomps + c);
        const auto *fl_348 = buffer.data(fl + 348 * ncomps + c);
        const auto *fl_349 = buffer.data(fl + 349 * ncomps + c);
        const auto *fl_350 = buffer.data(fl + 350 * ncomps + c);
        const auto *fl_360 = buffer.data(fl + 360 * ncomps + c);
        const auto *fl_361 = buffer.data(fl + 361 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, fk_145, fk_146, fk_147, \
                         fk_148, fk_149, fl_181, fl_182, fl_183, fl_184, \
                         fl_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * fk_145[k]
                       + fl_181[k];

            t_146[k] = -ab_x[k] * fk_146[k]
                       + fl_182[k];

            t_147[k] = -ab_x[k] * fk_147[k]
                       + fl_183[k];

            t_148[k] = -ab_x[k] * fk_148[k]
                       + fl_184[k];

            t_149[k] = -ab_x[k] * fk_149[k]
                       + fl_185[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, fk_150, fk_151, fk_152, \
                         fk_153, fk_154, fl_186, fl_187, fl_188, fl_189, \
                         fl_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * fk_150[k]
                       + fl_186[k];

            t_151[k] = -ab_x[k] * fk_151[k]
                       + fl_187[k];

            t_152[k] = -ab_x[k] * fk_152[k]
                       + fl_188[k];

            t_153[k] = -ab_x[k] * fk_153[k]
                       + fl_189[k];

            t_154[k] = -ab_x[k] * fk_154[k]
                       + fl_190[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, fk_155, fk_156, fk_157, \
                         fk_158, fk_159, fl_191, fl_192, fl_193, fl_194, \
                         fl_195 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * fk_155[k]
                       + fl_191[k];

            t_156[k] = -ab_x[k] * fk_156[k]
                       + fl_192[k];

            t_157[k] = -ab_x[k] * fk_157[k]
                       + fl_193[k];

            t_158[k] = -ab_x[k] * fk_158[k]
                       + fl_194[k];

            t_159[k] = -ab_x[k] * fk_159[k]
                       + fl_195[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, fk_160, fk_161, fk_162, \
                         fk_163, fk_164, fl_196, fl_197, fl_198, fl_199, \
                         fl_200 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * fk_160[k]
                       + fl_196[k];

            t_161[k] = -ab_x[k] * fk_161[k]
                       + fl_197[k];

            t_162[k] = -ab_x[k] * fk_162[k]
                       + fl_198[k];

            t_163[k] = -ab_x[k] * fk_163[k]
                       + fl_199[k];

            t_164[k] = -ab_x[k] * fk_164[k]
                       + fl_200[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, fk_165, fk_166, fk_167, \
                         fk_168, fk_169, fl_201, fl_202, fl_203, fl_204, \
                         fl_205 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * fk_165[k]
                       + fl_201[k];

            t_166[k] = -ab_x[k] * fk_166[k]
                       + fl_202[k];

            t_167[k] = -ab_x[k] * fk_167[k]
                       + fl_203[k];

            t_168[k] = -ab_x[k] * fk_168[k]
                       + fl_204[k];

            t_169[k] = -ab_x[k] * fk_169[k]
                       + fl_205[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, fk_170, fk_171, fk_172, \
                         fk_173, fk_174, fl_206, fl_207, fl_208, fl_209, \
                         fl_210 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * fk_170[k]
                       + fl_206[k];

            t_171[k] = -ab_x[k] * fk_171[k]
                       + fl_207[k];

            t_172[k] = -ab_x[k] * fk_172[k]
                       + fl_208[k];

            t_173[k] = -ab_x[k] * fk_173[k]
                       + fl_209[k];

            t_174[k] = -ab_x[k] * fk_174[k]
                       + fl_210[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, fk_175, fk_176, fk_177, \
                         fk_178, fk_179, fl_211, fl_212, fl_213, fl_214, \
                         fl_215 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * fk_175[k]
                       + fl_211[k];

            t_176[k] = -ab_x[k] * fk_176[k]
                       + fl_212[k];

            t_177[k] = -ab_x[k] * fk_177[k]
                       + fl_213[k];

            t_178[k] = -ab_x[k] * fk_178[k]
                       + fl_214[k];

            t_179[k] = -ab_x[k] * fk_179[k]
                       + fl_215[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, fk_180, fk_181, fk_182, \
                         fk_183, fk_184, fl_225, fl_226, fl_227, fl_228, \
                         fl_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * fk_180[k]
                       + fl_225[k];

            t_181[k] = -ab_x[k] * fk_181[k]
                       + fl_226[k];

            t_182[k] = -ab_x[k] * fk_182[k]
                       + fl_227[k];

            t_183[k] = -ab_x[k] * fk_183[k]
                       + fl_228[k];

            t_184[k] = -ab_x[k] * fk_184[k]
                       + fl_229[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, fk_185, fk_186, fk_187, \
                         fk_188, fk_189, fl_230, fl_231, fl_232, fl_233, \
                         fl_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * fk_185[k]
                       + fl_230[k];

            t_186[k] = -ab_x[k] * fk_186[k]
                       + fl_231[k];

            t_187[k] = -ab_x[k] * fk_187[k]
                       + fl_232[k];

            t_188[k] = -ab_x[k] * fk_188[k]
                       + fl_233[k];

            t_189[k] = -ab_x[k] * fk_189[k]
                       + fl_234[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, fk_190, fk_191, fk_192, \
                         fk_193, fk_194, fl_235, fl_236, fl_237, fl_238, \
                         fl_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * fk_190[k]
                       + fl_235[k];

            t_191[k] = -ab_x[k] * fk_191[k]
                       + fl_236[k];

            t_192[k] = -ab_x[k] * fk_192[k]
                       + fl_237[k];

            t_193[k] = -ab_x[k] * fk_193[k]
                       + fl_238[k];

            t_194[k] = -ab_x[k] * fk_194[k]
                       + fl_239[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, fk_195, fk_196, fk_197, \
                         fk_198, fk_199, fl_240, fl_241, fl_242, fl_243, \
                         fl_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * fk_195[k]
                       + fl_240[k];

            t_196[k] = -ab_x[k] * fk_196[k]
                       + fl_241[k];

            t_197[k] = -ab_x[k] * fk_197[k]
                       + fl_242[k];

            t_198[k] = -ab_x[k] * fk_198[k]
                       + fl_243[k];

            t_199[k] = -ab_x[k] * fk_199[k]
                       + fl_244[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, fk_200, fk_201, fk_202, \
                         fk_203, fk_204, fl_245, fl_246, fl_247, fl_248, \
                         fl_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * fk_200[k]
                       + fl_245[k];

            t_201[k] = -ab_x[k] * fk_201[k]
                       + fl_246[k];

            t_202[k] = -ab_x[k] * fk_202[k]
                       + fl_247[k];

            t_203[k] = -ab_x[k] * fk_203[k]
                       + fl_248[k];

            t_204[k] = -ab_x[k] * fk_204[k]
                       + fl_249[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, fk_205, fk_206, fk_207, \
                         fk_208, fk_209, fl_250, fl_251, fl_252, fl_253, \
                         fl_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * fk_205[k]
                       + fl_250[k];

            t_206[k] = -ab_x[k] * fk_206[k]
                       + fl_251[k];

            t_207[k] = -ab_x[k] * fk_207[k]
                       + fl_252[k];

            t_208[k] = -ab_x[k] * fk_208[k]
                       + fl_253[k];

            t_209[k] = -ab_x[k] * fk_209[k]
                       + fl_254[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, fk_210, fk_211, fk_212, \
                         fk_213, fk_214, fl_255, fl_256, fl_257, fl_258, \
                         fl_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * fk_210[k]
                       + fl_255[k];

            t_211[k] = -ab_x[k] * fk_211[k]
                       + fl_256[k];

            t_212[k] = -ab_x[k] * fk_212[k]
                       + fl_257[k];

            t_213[k] = -ab_x[k] * fk_213[k]
                       + fl_258[k];

            t_214[k] = -ab_x[k] * fk_214[k]
                       + fl_259[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, fk_215, fk_216, fk_217, \
                         fk_218, fk_219, fl_260, fl_270, fl_271, fl_272, \
                         fl_273 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * fk_215[k]
                       + fl_260[k];

            t_216[k] = -ab_x[k] * fk_216[k]
                       + fl_270[k];

            t_217[k] = -ab_x[k] * fk_217[k]
                       + fl_271[k];

            t_218[k] = -ab_x[k] * fk_218[k]
                       + fl_272[k];

            t_219[k] = -ab_x[k] * fk_219[k]
                       + fl_273[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, fk_220, fk_221, fk_222, \
                         fk_223, fk_224, fl_274, fl_275, fl_276, fl_277, \
                         fl_278 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_x[k] * fk_220[k]
                       + fl_274[k];

            t_221[k] = -ab_x[k] * fk_221[k]
                       + fl_275[k];

            t_222[k] = -ab_x[k] * fk_222[k]
                       + fl_276[k];

            t_223[k] = -ab_x[k] * fk_223[k]
                       + fl_277[k];

            t_224[k] = -ab_x[k] * fk_224[k]
                       + fl_278[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, fk_225, fk_226, fk_227, \
                         fk_228, fk_229, fl_279, fl_280, fl_281, fl_282, \
                         fl_283 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_x[k] * fk_225[k]
                       + fl_279[k];

            t_226[k] = -ab_x[k] * fk_226[k]
                       + fl_280[k];

            t_227[k] = -ab_x[k] * fk_227[k]
                       + fl_281[k];

            t_228[k] = -ab_x[k] * fk_228[k]
                       + fl_282[k];

            t_229[k] = -ab_x[k] * fk_229[k]
                       + fl_283[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, fk_230, fk_231, fk_232, \
                         fk_233, fk_234, fl_284, fl_285, fl_286, fl_287, \
                         fl_288 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_x[k] * fk_230[k]
                       + fl_284[k];

            t_231[k] = -ab_x[k] * fk_231[k]
                       + fl_285[k];

            t_232[k] = -ab_x[k] * fk_232[k]
                       + fl_286[k];

            t_233[k] = -ab_x[k] * fk_233[k]
                       + fl_287[k];

            t_234[k] = -ab_x[k] * fk_234[k]
                       + fl_288[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, fk_235, fk_236, fk_237, \
                         fk_238, fk_239, fl_289, fl_290, fl_291, fl_292, \
                         fl_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_x[k] * fk_235[k]
                       + fl_289[k];

            t_236[k] = -ab_x[k] * fk_236[k]
                       + fl_290[k];

            t_237[k] = -ab_x[k] * fk_237[k]
                       + fl_291[k];

            t_238[k] = -ab_x[k] * fk_238[k]
                       + fl_292[k];

            t_239[k] = -ab_x[k] * fk_239[k]
                       + fl_293[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, fk_240, fk_241, fk_242, \
                         fk_243, fk_244, fl_294, fl_295, fl_296, fl_297, \
                         fl_298 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_x[k] * fk_240[k]
                       + fl_294[k];

            t_241[k] = -ab_x[k] * fk_241[k]
                       + fl_295[k];

            t_242[k] = -ab_x[k] * fk_242[k]
                       + fl_296[k];

            t_243[k] = -ab_x[k] * fk_243[k]
                       + fl_297[k];

            t_244[k] = -ab_x[k] * fk_244[k]
                       + fl_298[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, fk_245, fk_246, fk_247, \
                         fk_248, fk_249, fl_299, fl_300, fl_301, fl_302, \
                         fl_303 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_x[k] * fk_245[k]
                       + fl_299[k];

            t_246[k] = -ab_x[k] * fk_246[k]
                       + fl_300[k];

            t_247[k] = -ab_x[k] * fk_247[k]
                       + fl_301[k];

            t_248[k] = -ab_x[k] * fk_248[k]
                       + fl_302[k];

            t_249[k] = -ab_x[k] * fk_249[k]
                       + fl_303[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, fk_250, fk_251, fk_252, \
                         fk_253, fk_254, fl_304, fl_305, fl_315, fl_316, \
                         fl_317 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_x[k] * fk_250[k]
                       + fl_304[k];

            t_251[k] = -ab_x[k] * fk_251[k]
                       + fl_305[k];

            t_252[k] = -ab_x[k] * fk_252[k]
                       + fl_315[k];

            t_253[k] = -ab_x[k] * fk_253[k]
                       + fl_316[k];

            t_254[k] = -ab_x[k] * fk_254[k]
                       + fl_317[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, fk_255, fk_256, fk_257, \
                         fk_258, fk_259, fl_318, fl_319, fl_320, fl_321, \
                         fl_322 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_x[k] * fk_255[k]
                       + fl_318[k];

            t_256[k] = -ab_x[k] * fk_256[k]
                       + fl_319[k];

            t_257[k] = -ab_x[k] * fk_257[k]
                       + fl_320[k];

            t_258[k] = -ab_x[k] * fk_258[k]
                       + fl_321[k];

            t_259[k] = -ab_x[k] * fk_259[k]
                       + fl_322[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, fk_260, fk_261, fk_262, \
                         fk_263, fk_264, fl_323, fl_324, fl_325, fl_326, \
                         fl_327 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_x[k] * fk_260[k]
                       + fl_323[k];

            t_261[k] = -ab_x[k] * fk_261[k]
                       + fl_324[k];

            t_262[k] = -ab_x[k] * fk_262[k]
                       + fl_325[k];

            t_263[k] = -ab_x[k] * fk_263[k]
                       + fl_326[k];

            t_264[k] = -ab_x[k] * fk_264[k]
                       + fl_327[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, fk_265, fk_266, fk_267, \
                         fk_268, fk_269, fl_328, fl_329, fl_330, fl_331, \
                         fl_332 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_x[k] * fk_265[k]
                       + fl_328[k];

            t_266[k] = -ab_x[k] * fk_266[k]
                       + fl_329[k];

            t_267[k] = -ab_x[k] * fk_267[k]
                       + fl_330[k];

            t_268[k] = -ab_x[k] * fk_268[k]
                       + fl_331[k];

            t_269[k] = -ab_x[k] * fk_269[k]
                       + fl_332[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, fk_270, fk_271, fk_272, \
                         fk_273, fk_274, fl_333, fl_334, fl_335, fl_336, \
                         fl_337 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_x[k] * fk_270[k]
                       + fl_333[k];

            t_271[k] = -ab_x[k] * fk_271[k]
                       + fl_334[k];

            t_272[k] = -ab_x[k] * fk_272[k]
                       + fl_335[k];

            t_273[k] = -ab_x[k] * fk_273[k]
                       + fl_336[k];

            t_274[k] = -ab_x[k] * fk_274[k]
                       + fl_337[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, fk_275, fk_276, fk_277, \
                         fk_278, fk_279, fl_338, fl_339, fl_340, fl_341, \
                         fl_342 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_x[k] * fk_275[k]
                       + fl_338[k];

            t_276[k] = -ab_x[k] * fk_276[k]
                       + fl_339[k];

            t_277[k] = -ab_x[k] * fk_277[k]
                       + fl_340[k];

            t_278[k] = -ab_x[k] * fk_278[k]
                       + fl_341[k];

            t_279[k] = -ab_x[k] * fk_279[k]
                       + fl_342[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, fk_280, fk_281, fk_282, \
                         fk_283, fk_284, fl_343, fl_344, fl_345, fl_346, \
                         fl_347 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_x[k] * fk_280[k]
                       + fl_343[k];

            t_281[k] = -ab_x[k] * fk_281[k]
                       + fl_344[k];

            t_282[k] = -ab_x[k] * fk_282[k]
                       + fl_345[k];

            t_283[k] = -ab_x[k] * fk_283[k]
                       + fl_346[k];

            t_284[k] = -ab_x[k] * fk_284[k]
                       + fl_347[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, fk_285, fk_286, fk_287, \
                         fk_288, fk_289, fl_348, fl_349, fl_350, fl_360, \
                         fl_361 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_x[k] * fk_285[k]
                       + fl_348[k];

            t_286[k] = -ab_x[k] * fk_286[k]
                       + fl_349[k];

            t_287[k] = -ab_x[k] * fk_287[k]
                       + fl_350[k];

            t_288[k] = -ab_x[k] * fk_288[k]
                       + fl_360[k];

            t_289[k] = -ab_x[k] * fk_289[k]
                       + fl_361[k];
        }
    }
}

static auto
compute_hrr_gk_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fk, const size_t fl, const size_t ncomps,
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

        const auto *fk_216 = buffer.data(fk + 216 * ncomps + c);
        const auto *fk_217 = buffer.data(fk + 217 * ncomps + c);
        const auto *fk_218 = buffer.data(fk + 218 * ncomps + c);
        const auto *fk_219 = buffer.data(fk + 219 * ncomps + c);
        const auto *fk_220 = buffer.data(fk + 220 * ncomps + c);
        const auto *fk_221 = buffer.data(fk + 221 * ncomps + c);
        const auto *fk_222 = buffer.data(fk + 222 * ncomps + c);
        const auto *fk_223 = buffer.data(fk + 223 * ncomps + c);
        const auto *fk_224 = buffer.data(fk + 224 * ncomps + c);
        const auto *fk_225 = buffer.data(fk + 225 * ncomps + c);
        const auto *fk_226 = buffer.data(fk + 226 * ncomps + c);
        const auto *fk_227 = buffer.data(fk + 227 * ncomps + c);
        const auto *fk_228 = buffer.data(fk + 228 * ncomps + c);
        const auto *fk_229 = buffer.data(fk + 229 * ncomps + c);
        const auto *fk_230 = buffer.data(fk + 230 * ncomps + c);
        const auto *fk_231 = buffer.data(fk + 231 * ncomps + c);
        const auto *fk_232 = buffer.data(fk + 232 * ncomps + c);
        const auto *fk_233 = buffer.data(fk + 233 * ncomps + c);
        const auto *fk_234 = buffer.data(fk + 234 * ncomps + c);
        const auto *fk_235 = buffer.data(fk + 235 * ncomps + c);
        const auto *fk_236 = buffer.data(fk + 236 * ncomps + c);
        const auto *fk_237 = buffer.data(fk + 237 * ncomps + c);
        const auto *fk_238 = buffer.data(fk + 238 * ncomps + c);
        const auto *fk_239 = buffer.data(fk + 239 * ncomps + c);
        const auto *fk_240 = buffer.data(fk + 240 * ncomps + c);
        const auto *fk_241 = buffer.data(fk + 241 * ncomps + c);
        const auto *fk_242 = buffer.data(fk + 242 * ncomps + c);
        const auto *fk_243 = buffer.data(fk + 243 * ncomps + c);
        const auto *fk_244 = buffer.data(fk + 244 * ncomps + c);
        const auto *fk_245 = buffer.data(fk + 245 * ncomps + c);
        const auto *fk_246 = buffer.data(fk + 246 * ncomps + c);
        const auto *fk_247 = buffer.data(fk + 247 * ncomps + c);
        const auto *fk_248 = buffer.data(fk + 248 * ncomps + c);
        const auto *fk_249 = buffer.data(fk + 249 * ncomps + c);
        const auto *fk_250 = buffer.data(fk + 250 * ncomps + c);
        const auto *fk_251 = buffer.data(fk + 251 * ncomps + c);
        const auto *fk_252 = buffer.data(fk + 252 * ncomps + c);
        const auto *fk_253 = buffer.data(fk + 253 * ncomps + c);
        const auto *fk_254 = buffer.data(fk + 254 * ncomps + c);
        const auto *fk_255 = buffer.data(fk + 255 * ncomps + c);
        const auto *fk_256 = buffer.data(fk + 256 * ncomps + c);
        const auto *fk_257 = buffer.data(fk + 257 * ncomps + c);
        const auto *fk_258 = buffer.data(fk + 258 * ncomps + c);
        const auto *fk_259 = buffer.data(fk + 259 * ncomps + c);
        const auto *fk_260 = buffer.data(fk + 260 * ncomps + c);
        const auto *fk_261 = buffer.data(fk + 261 * ncomps + c);
        const auto *fk_262 = buffer.data(fk + 262 * ncomps + c);
        const auto *fk_263 = buffer.data(fk + 263 * ncomps + c);
        const auto *fk_264 = buffer.data(fk + 264 * ncomps + c);
        const auto *fk_265 = buffer.data(fk + 265 * ncomps + c);
        const auto *fk_266 = buffer.data(fk + 266 * ncomps + c);
        const auto *fk_267 = buffer.data(fk + 267 * ncomps + c);
        const auto *fk_268 = buffer.data(fk + 268 * ncomps + c);
        const auto *fk_269 = buffer.data(fk + 269 * ncomps + c);
        const auto *fk_270 = buffer.data(fk + 270 * ncomps + c);
        const auto *fk_271 = buffer.data(fk + 271 * ncomps + c);
        const auto *fk_272 = buffer.data(fk + 272 * ncomps + c);
        const auto *fk_273 = buffer.data(fk + 273 * ncomps + c);
        const auto *fk_274 = buffer.data(fk + 274 * ncomps + c);
        const auto *fk_275 = buffer.data(fk + 275 * ncomps + c);
        const auto *fk_276 = buffer.data(fk + 276 * ncomps + c);
        const auto *fk_277 = buffer.data(fk + 277 * ncomps + c);
        const auto *fk_278 = buffer.data(fk + 278 * ncomps + c);
        const auto *fk_279 = buffer.data(fk + 279 * ncomps + c);
        const auto *fk_280 = buffer.data(fk + 280 * ncomps + c);
        const auto *fk_281 = buffer.data(fk + 281 * ncomps + c);
        const auto *fk_282 = buffer.data(fk + 282 * ncomps + c);
        const auto *fk_283 = buffer.data(fk + 283 * ncomps + c);
        const auto *fk_284 = buffer.data(fk + 284 * ncomps + c);
        const auto *fk_285 = buffer.data(fk + 285 * ncomps + c);
        const auto *fk_286 = buffer.data(fk + 286 * ncomps + c);
        const auto *fk_287 = buffer.data(fk + 287 * ncomps + c);
        const auto *fk_288 = buffer.data(fk + 288 * ncomps + c);
        const auto *fk_289 = buffer.data(fk + 289 * ncomps + c);
        const auto *fk_290 = buffer.data(fk + 290 * ncomps + c);
        const auto *fk_291 = buffer.data(fk + 291 * ncomps + c);
        const auto *fk_292 = buffer.data(fk + 292 * ncomps + c);
        const auto *fk_293 = buffer.data(fk + 293 * ncomps + c);
        const auto *fk_294 = buffer.data(fk + 294 * ncomps + c);
        const auto *fk_295 = buffer.data(fk + 295 * ncomps + c);
        const auto *fk_296 = buffer.data(fk + 296 * ncomps + c);
        const auto *fk_297 = buffer.data(fk + 297 * ncomps + c);
        const auto *fk_298 = buffer.data(fk + 298 * ncomps + c);
        const auto *fk_299 = buffer.data(fk + 299 * ncomps + c);
        const auto *fk_300 = buffer.data(fk + 300 * ncomps + c);
        const auto *fk_301 = buffer.data(fk + 301 * ncomps + c);
        const auto *fk_302 = buffer.data(fk + 302 * ncomps + c);
        const auto *fk_303 = buffer.data(fk + 303 * ncomps + c);
        const auto *fk_304 = buffer.data(fk + 304 * ncomps + c);
        const auto *fk_305 = buffer.data(fk + 305 * ncomps + c);
        const auto *fk_306 = buffer.data(fk + 306 * ncomps + c);
        const auto *fk_307 = buffer.data(fk + 307 * ncomps + c);
        const auto *fk_308 = buffer.data(fk + 308 * ncomps + c);
        const auto *fk_309 = buffer.data(fk + 309 * ncomps + c);
        const auto *fk_310 = buffer.data(fk + 310 * ncomps + c);
        const auto *fk_311 = buffer.data(fk + 311 * ncomps + c);
        const auto *fk_312 = buffer.data(fk + 312 * ncomps + c);
        const auto *fk_313 = buffer.data(fk + 313 * ncomps + c);
        const auto *fk_314 = buffer.data(fk + 314 * ncomps + c);
        const auto *fk_315 = buffer.data(fk + 315 * ncomps + c);
        const auto *fk_316 = buffer.data(fk + 316 * ncomps + c);
        const auto *fk_317 = buffer.data(fk + 317 * ncomps + c);
        const auto *fk_318 = buffer.data(fk + 318 * ncomps + c);
        const auto *fk_319 = buffer.data(fk + 319 * ncomps + c);
        const auto *fk_320 = buffer.data(fk + 320 * ncomps + c);
        const auto *fk_321 = buffer.data(fk + 321 * ncomps + c);
        const auto *fk_322 = buffer.data(fk + 322 * ncomps + c);
        const auto *fk_323 = buffer.data(fk + 323 * ncomps + c);
        const auto *fk_324 = buffer.data(fk + 324 * ncomps + c);
        const auto *fk_325 = buffer.data(fk + 325 * ncomps + c);
        const auto *fk_326 = buffer.data(fk + 326 * ncomps + c);
        const auto *fk_327 = buffer.data(fk + 327 * ncomps + c);
        const auto *fk_328 = buffer.data(fk + 328 * ncomps + c);
        const auto *fk_329 = buffer.data(fk + 329 * ncomps + c);
        const auto *fk_330 = buffer.data(fk + 330 * ncomps + c);
        const auto *fk_331 = buffer.data(fk + 331 * ncomps + c);
        const auto *fk_332 = buffer.data(fk + 332 * ncomps + c);
        const auto *fk_333 = buffer.data(fk + 333 * ncomps + c);
        const auto *fk_334 = buffer.data(fk + 334 * ncomps + c);
        const auto *fk_335 = buffer.data(fk + 335 * ncomps + c);
        const auto *fk_336 = buffer.data(fk + 336 * ncomps + c);
        const auto *fk_337 = buffer.data(fk + 337 * ncomps + c);
        const auto *fk_338 = buffer.data(fk + 338 * ncomps + c);
        const auto *fk_339 = buffer.data(fk + 339 * ncomps + c);
        const auto *fk_340 = buffer.data(fk + 340 * ncomps + c);
        const auto *fk_341 = buffer.data(fk + 341 * ncomps + c);
        const auto *fk_342 = buffer.data(fk + 342 * ncomps + c);
        const auto *fk_343 = buffer.data(fk + 343 * ncomps + c);
        const auto *fk_344 = buffer.data(fk + 344 * ncomps + c);
        const auto *fk_345 = buffer.data(fk + 345 * ncomps + c);
        const auto *fk_346 = buffer.data(fk + 346 * ncomps + c);
        const auto *fk_347 = buffer.data(fk + 347 * ncomps + c);
        const auto *fk_348 = buffer.data(fk + 348 * ncomps + c);
        const auto *fk_349 = buffer.data(fk + 349 * ncomps + c);
        const auto *fk_350 = buffer.data(fk + 350 * ncomps + c);
        const auto *fk_351 = buffer.data(fk + 351 * ncomps + c);
        const auto *fk_352 = buffer.data(fk + 352 * ncomps + c);
        const auto *fk_353 = buffer.data(fk + 353 * ncomps + c);
        const auto *fk_354 = buffer.data(fk + 354 * ncomps + c);
        const auto *fk_355 = buffer.data(fk + 355 * ncomps + c);
        const auto *fk_356 = buffer.data(fk + 356 * ncomps + c);
        const auto *fk_357 = buffer.data(fk + 357 * ncomps + c);
        const auto *fk_358 = buffer.data(fk + 358 * ncomps + c);
        const auto *fk_359 = buffer.data(fk + 359 * ncomps + c);

        const auto *fl_271 = buffer.data(fl + 271 * ncomps + c);
        const auto *fl_273 = buffer.data(fl + 273 * ncomps + c);
        const auto *fl_274 = buffer.data(fl + 274 * ncomps + c);
        const auto *fl_276 = buffer.data(fl + 276 * ncomps + c);
        const auto *fl_277 = buffer.data(fl + 277 * ncomps + c);
        const auto *fl_278 = buffer.data(fl + 278 * ncomps + c);
        const auto *fl_280 = buffer.data(fl + 280 * ncomps + c);
        const auto *fl_281 = buffer.data(fl + 281 * ncomps + c);
        const auto *fl_282 = buffer.data(fl + 282 * ncomps + c);
        const auto *fl_283 = buffer.data(fl + 283 * ncomps + c);
        const auto *fl_285 = buffer.data(fl + 285 * ncomps + c);
        const auto *fl_286 = buffer.data(fl + 286 * ncomps + c);
        const auto *fl_287 = buffer.data(fl + 287 * ncomps + c);
        const auto *fl_288 = buffer.data(fl + 288 * ncomps + c);
        const auto *fl_289 = buffer.data(fl + 289 * ncomps + c);
        const auto *fl_291 = buffer.data(fl + 291 * ncomps + c);
        const auto *fl_292 = buffer.data(fl + 292 * ncomps + c);
        const auto *fl_293 = buffer.data(fl + 293 * ncomps + c);
        const auto *fl_294 = buffer.data(fl + 294 * ncomps + c);
        const auto *fl_295 = buffer.data(fl + 295 * ncomps + c);
        const auto *fl_296 = buffer.data(fl + 296 * ncomps + c);
        const auto *fl_298 = buffer.data(fl + 298 * ncomps + c);
        const auto *fl_299 = buffer.data(fl + 299 * ncomps + c);
        const auto *fl_300 = buffer.data(fl + 300 * ncomps + c);
        const auto *fl_301 = buffer.data(fl + 301 * ncomps + c);
        const auto *fl_302 = buffer.data(fl + 302 * ncomps + c);
        const auto *fl_303 = buffer.data(fl + 303 * ncomps + c);
        const auto *fl_304 = buffer.data(fl + 304 * ncomps + c);
        const auto *fl_306 = buffer.data(fl + 306 * ncomps + c);
        const auto *fl_307 = buffer.data(fl + 307 * ncomps + c);
        const auto *fl_308 = buffer.data(fl + 308 * ncomps + c);
        const auto *fl_309 = buffer.data(fl + 309 * ncomps + c);
        const auto *fl_310 = buffer.data(fl + 310 * ncomps + c);
        const auto *fl_311 = buffer.data(fl + 311 * ncomps + c);
        const auto *fl_312 = buffer.data(fl + 312 * ncomps + c);
        const auto *fl_313 = buffer.data(fl + 313 * ncomps + c);
        const auto *fl_316 = buffer.data(fl + 316 * ncomps + c);
        const auto *fl_318 = buffer.data(fl + 318 * ncomps + c);
        const auto *fl_319 = buffer.data(fl + 319 * ncomps + c);
        const auto *fl_321 = buffer.data(fl + 321 * ncomps + c);
        const auto *fl_322 = buffer.data(fl + 322 * ncomps + c);
        const auto *fl_323 = buffer.data(fl + 323 * ncomps + c);
        const auto *fl_325 = buffer.data(fl + 325 * ncomps + c);
        const auto *fl_326 = buffer.data(fl + 326 * ncomps + c);
        const auto *fl_327 = buffer.data(fl + 327 * ncomps + c);
        const auto *fl_328 = buffer.data(fl + 328 * ncomps + c);
        const auto *fl_330 = buffer.data(fl + 330 * ncomps + c);
        const auto *fl_331 = buffer.data(fl + 331 * ncomps + c);
        const auto *fl_332 = buffer.data(fl + 332 * ncomps + c);
        const auto *fl_333 = buffer.data(fl + 333 * ncomps + c);
        const auto *fl_334 = buffer.data(fl + 334 * ncomps + c);
        const auto *fl_336 = buffer.data(fl + 336 * ncomps + c);
        const auto *fl_337 = buffer.data(fl + 337 * ncomps + c);
        const auto *fl_338 = buffer.data(fl + 338 * ncomps + c);
        const auto *fl_339 = buffer.data(fl + 339 * ncomps + c);
        const auto *fl_340 = buffer.data(fl + 340 * ncomps + c);
        const auto *fl_341 = buffer.data(fl + 341 * ncomps + c);
        const auto *fl_343 = buffer.data(fl + 343 * ncomps + c);
        const auto *fl_344 = buffer.data(fl + 344 * ncomps + c);
        const auto *fl_345 = buffer.data(fl + 345 * ncomps + c);
        const auto *fl_346 = buffer.data(fl + 346 * ncomps + c);
        const auto *fl_347 = buffer.data(fl + 347 * ncomps + c);
        const auto *fl_348 = buffer.data(fl + 348 * ncomps + c);
        const auto *fl_349 = buffer.data(fl + 349 * ncomps + c);
        const auto *fl_351 = buffer.data(fl + 351 * ncomps + c);
        const auto *fl_352 = buffer.data(fl + 352 * ncomps + c);
        const auto *fl_353 = buffer.data(fl + 353 * ncomps + c);
        const auto *fl_354 = buffer.data(fl + 354 * ncomps + c);
        const auto *fl_355 = buffer.data(fl + 355 * ncomps + c);
        const auto *fl_356 = buffer.data(fl + 356 * ncomps + c);
        const auto *fl_357 = buffer.data(fl + 357 * ncomps + c);
        const auto *fl_358 = buffer.data(fl + 358 * ncomps + c);
        const auto *fl_361 = buffer.data(fl + 361 * ncomps + c);
        const auto *fl_362 = buffer.data(fl + 362 * ncomps + c);
        const auto *fl_363 = buffer.data(fl + 363 * ncomps + c);
        const auto *fl_364 = buffer.data(fl + 364 * ncomps + c);
        const auto *fl_365 = buffer.data(fl + 365 * ncomps + c);
        const auto *fl_366 = buffer.data(fl + 366 * ncomps + c);
        const auto *fl_367 = buffer.data(fl + 367 * ncomps + c);
        const auto *fl_368 = buffer.data(fl + 368 * ncomps + c);
        const auto *fl_369 = buffer.data(fl + 369 * ncomps + c);
        const auto *fl_370 = buffer.data(fl + 370 * ncomps + c);
        const auto *fl_371 = buffer.data(fl + 371 * ncomps + c);
        const auto *fl_372 = buffer.data(fl + 372 * ncomps + c);
        const auto *fl_373 = buffer.data(fl + 373 * ncomps + c);
        const auto *fl_374 = buffer.data(fl + 374 * ncomps + c);
        const auto *fl_375 = buffer.data(fl + 375 * ncomps + c);
        const auto *fl_376 = buffer.data(fl + 376 * ncomps + c);
        const auto *fl_377 = buffer.data(fl + 377 * ncomps + c);
        const auto *fl_378 = buffer.data(fl + 378 * ncomps + c);
        const auto *fl_379 = buffer.data(fl + 379 * ncomps + c);
        const auto *fl_380 = buffer.data(fl + 380 * ncomps + c);
        const auto *fl_381 = buffer.data(fl + 381 * ncomps + c);
        const auto *fl_382 = buffer.data(fl + 382 * ncomps + c);
        const auto *fl_383 = buffer.data(fl + 383 * ncomps + c);
        const auto *fl_384 = buffer.data(fl + 384 * ncomps + c);
        const auto *fl_385 = buffer.data(fl + 385 * ncomps + c);
        const auto *fl_386 = buffer.data(fl + 386 * ncomps + c);
        const auto *fl_387 = buffer.data(fl + 387 * ncomps + c);
        const auto *fl_388 = buffer.data(fl + 388 * ncomps + c);
        const auto *fl_389 = buffer.data(fl + 389 * ncomps + c);
        const auto *fl_390 = buffer.data(fl + 390 * ncomps + c);
        const auto *fl_391 = buffer.data(fl + 391 * ncomps + c);
        const auto *fl_392 = buffer.data(fl + 392 * ncomps + c);
        const auto *fl_393 = buffer.data(fl + 393 * ncomps + c);
        const auto *fl_394 = buffer.data(fl + 394 * ncomps + c);
        const auto *fl_395 = buffer.data(fl + 395 * ncomps + c);
        const auto *fl_405 = buffer.data(fl + 405 * ncomps + c);
        const auto *fl_406 = buffer.data(fl + 406 * ncomps + c);
        const auto *fl_407 = buffer.data(fl + 407 * ncomps + c);
        const auto *fl_408 = buffer.data(fl + 408 * ncomps + c);
        const auto *fl_409 = buffer.data(fl + 409 * ncomps + c);
        const auto *fl_410 = buffer.data(fl + 410 * ncomps + c);
        const auto *fl_411 = buffer.data(fl + 411 * ncomps + c);
        const auto *fl_412 = buffer.data(fl + 412 * ncomps + c);
        const auto *fl_413 = buffer.data(fl + 413 * ncomps + c);
        const auto *fl_414 = buffer.data(fl + 414 * ncomps + c);
        const auto *fl_415 = buffer.data(fl + 415 * ncomps + c);
        const auto *fl_416 = buffer.data(fl + 416 * ncomps + c);
        const auto *fl_417 = buffer.data(fl + 417 * ncomps + c);
        const auto *fl_418 = buffer.data(fl + 418 * ncomps + c);
        const auto *fl_419 = buffer.data(fl + 419 * ncomps + c);
        const auto *fl_420 = buffer.data(fl + 420 * ncomps + c);
        const auto *fl_421 = buffer.data(fl + 421 * ncomps + c);
        const auto *fl_422 = buffer.data(fl + 422 * ncomps + c);
        const auto *fl_423 = buffer.data(fl + 423 * ncomps + c);
        const auto *fl_424 = buffer.data(fl + 424 * ncomps + c);
        const auto *fl_425 = buffer.data(fl + 425 * ncomps + c);
        const auto *fl_426 = buffer.data(fl + 426 * ncomps + c);
        const auto *fl_427 = buffer.data(fl + 427 * ncomps + c);
        const auto *fl_428 = buffer.data(fl + 428 * ncomps + c);
        const auto *fl_429 = buffer.data(fl + 429 * ncomps + c);
        const auto *fl_430 = buffer.data(fl + 430 * ncomps + c);
        const auto *fl_431 = buffer.data(fl + 431 * ncomps + c);
        const auto *fl_432 = buffer.data(fl + 432 * ncomps + c);
        const auto *fl_433 = buffer.data(fl + 433 * ncomps + c);
        const auto *fl_434 = buffer.data(fl + 434 * ncomps + c);
        const auto *fl_435 = buffer.data(fl + 435 * ncomps + c);
        const auto *fl_436 = buffer.data(fl + 436 * ncomps + c);
        const auto *fl_437 = buffer.data(fl + 437 * ncomps + c);
        const auto *fl_438 = buffer.data(fl + 438 * ncomps + c);
        const auto *fl_439 = buffer.data(fl + 439 * ncomps + c);
        const auto *fl_440 = buffer.data(fl + 440 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, fk_290, fk_291, fk_292, \
                         fk_293, fk_294, fl_362, fl_363, fl_364, fl_365, \
                         fl_366 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_x[k] * fk_290[k]
                       + fl_362[k];

            t_291[k] = -ab_x[k] * fk_291[k]
                       + fl_363[k];

            t_292[k] = -ab_x[k] * fk_292[k]
                       + fl_364[k];

            t_293[k] = -ab_x[k] * fk_293[k]
                       + fl_365[k];

            t_294[k] = -ab_x[k] * fk_294[k]
                       + fl_366[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, fk_295, fk_296, fk_297, \
                         fk_298, fk_299, fl_367, fl_368, fl_369, fl_370, \
                         fl_371 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_x[k] * fk_295[k]
                       + fl_367[k];

            t_296[k] = -ab_x[k] * fk_296[k]
                       + fl_368[k];

            t_297[k] = -ab_x[k] * fk_297[k]
                       + fl_369[k];

            t_298[k] = -ab_x[k] * fk_298[k]
                       + fl_370[k];

            t_299[k] = -ab_x[k] * fk_299[k]
                       + fl_371[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, fk_300, fk_301, fk_302, \
                         fk_303, fk_304, fl_372, fl_373, fl_374, fl_375, \
                         fl_376 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_x[k] * fk_300[k]
                       + fl_372[k];

            t_301[k] = -ab_x[k] * fk_301[k]
                       + fl_373[k];

            t_302[k] = -ab_x[k] * fk_302[k]
                       + fl_374[k];

            t_303[k] = -ab_x[k] * fk_303[k]
                       + fl_375[k];

            t_304[k] = -ab_x[k] * fk_304[k]
                       + fl_376[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, fk_305, fk_306, fk_307, \
                         fk_308, fk_309, fl_377, fl_378, fl_379, fl_380, \
                         fl_381 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_x[k] * fk_305[k]
                       + fl_377[k];

            t_306[k] = -ab_x[k] * fk_306[k]
                       + fl_378[k];

            t_307[k] = -ab_x[k] * fk_307[k]
                       + fl_379[k];

            t_308[k] = -ab_x[k] * fk_308[k]
                       + fl_380[k];

            t_309[k] = -ab_x[k] * fk_309[k]
                       + fl_381[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, fk_310, fk_311, fk_312, \
                         fk_313, fk_314, fl_382, fl_383, fl_384, fl_385, \
                         fl_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_x[k] * fk_310[k]
                       + fl_382[k];

            t_311[k] = -ab_x[k] * fk_311[k]
                       + fl_383[k];

            t_312[k] = -ab_x[k] * fk_312[k]
                       + fl_384[k];

            t_313[k] = -ab_x[k] * fk_313[k]
                       + fl_385[k];

            t_314[k] = -ab_x[k] * fk_314[k]
                       + fl_386[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, fk_315, fk_316, fk_317, \
                         fk_318, fk_319, fl_387, fl_388, fl_389, fl_390, \
                         fl_391 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_x[k] * fk_315[k]
                       + fl_387[k];

            t_316[k] = -ab_x[k] * fk_316[k]
                       + fl_388[k];

            t_317[k] = -ab_x[k] * fk_317[k]
                       + fl_389[k];

            t_318[k] = -ab_x[k] * fk_318[k]
                       + fl_390[k];

            t_319[k] = -ab_x[k] * fk_319[k]
                       + fl_391[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, fk_320, fk_321, fk_322, \
                         fk_323, fk_324, fl_392, fl_393, fl_394, fl_395, \
                         fl_405 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_x[k] * fk_320[k]
                       + fl_392[k];

            t_321[k] = -ab_x[k] * fk_321[k]
                       + fl_393[k];

            t_322[k] = -ab_x[k] * fk_322[k]
                       + fl_394[k];

            t_323[k] = -ab_x[k] * fk_323[k]
                       + fl_395[k];

            t_324[k] = -ab_x[k] * fk_324[k]
                       + fl_405[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, fk_325, fk_326, fk_327, \
                         fk_328, fk_329, fl_406, fl_407, fl_408, fl_409, \
                         fl_410 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_x[k] * fk_325[k]
                       + fl_406[k];

            t_326[k] = -ab_x[k] * fk_326[k]
                       + fl_407[k];

            t_327[k] = -ab_x[k] * fk_327[k]
                       + fl_408[k];

            t_328[k] = -ab_x[k] * fk_328[k]
                       + fl_409[k];

            t_329[k] = -ab_x[k] * fk_329[k]
                       + fl_410[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, fk_330, fk_331, fk_332, \
                         fk_333, fk_334, fl_411, fl_412, fl_413, fl_414, \
                         fl_415 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = -ab_x[k] * fk_330[k]
                       + fl_411[k];

            t_331[k] = -ab_x[k] * fk_331[k]
                       + fl_412[k];

            t_332[k] = -ab_x[k] * fk_332[k]
                       + fl_413[k];

            t_333[k] = -ab_x[k] * fk_333[k]
                       + fl_414[k];

            t_334[k] = -ab_x[k] * fk_334[k]
                       + fl_415[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, fk_335, fk_336, fk_337, \
                         fk_338, fk_339, fl_416, fl_417, fl_418, fl_419, \
                         fl_420 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = -ab_x[k] * fk_335[k]
                       + fl_416[k];

            t_336[k] = -ab_x[k] * fk_336[k]
                       + fl_417[k];

            t_337[k] = -ab_x[k] * fk_337[k]
                       + fl_418[k];

            t_338[k] = -ab_x[k] * fk_338[k]
                       + fl_419[k];

            t_339[k] = -ab_x[k] * fk_339[k]
                       + fl_420[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, fk_340, fk_341, fk_342, \
                         fk_343, fk_344, fl_421, fl_422, fl_423, fl_424, \
                         fl_425 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = -ab_x[k] * fk_340[k]
                       + fl_421[k];

            t_341[k] = -ab_x[k] * fk_341[k]
                       + fl_422[k];

            t_342[k] = -ab_x[k] * fk_342[k]
                       + fl_423[k];

            t_343[k] = -ab_x[k] * fk_343[k]
                       + fl_424[k];

            t_344[k] = -ab_x[k] * fk_344[k]
                       + fl_425[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, fk_345, fk_346, fk_347, \
                         fk_348, fk_349, fl_426, fl_427, fl_428, fl_429, \
                         fl_430 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = -ab_x[k] * fk_345[k]
                       + fl_426[k];

            t_346[k] = -ab_x[k] * fk_346[k]
                       + fl_427[k];

            t_347[k] = -ab_x[k] * fk_347[k]
                       + fl_428[k];

            t_348[k] = -ab_x[k] * fk_348[k]
                       + fl_429[k];

            t_349[k] = -ab_x[k] * fk_349[k]
                       + fl_430[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, fk_350, fk_351, fk_352, \
                         fk_353, fk_354, fl_431, fl_432, fl_433, fl_434, \
                         fl_435 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = -ab_x[k] * fk_350[k]
                       + fl_431[k];

            t_351[k] = -ab_x[k] * fk_351[k]
                       + fl_432[k];

            t_352[k] = -ab_x[k] * fk_352[k]
                       + fl_433[k];

            t_353[k] = -ab_x[k] * fk_353[k]
                       + fl_434[k];

            t_354[k] = -ab_x[k] * fk_354[k]
                       + fl_435[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, fk_355, fk_356, fk_357, \
                         fk_358, fk_359, fl_436, fl_437, fl_438, fl_439, \
                         fl_440 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = -ab_x[k] * fk_355[k]
                       + fl_436[k];

            t_356[k] = -ab_x[k] * fk_356[k]
                       + fl_437[k];

            t_357[k] = -ab_x[k] * fk_357[k]
                       + fl_438[k];

            t_358[k] = -ab_x[k] * fk_358[k]
                       + fl_439[k];

            t_359[k] = -ab_x[k] * fk_359[k]
                       + fl_440[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_y, fk_216, fk_217, fk_218, \
                         fk_219, fk_220, fl_271, fl_273, fl_274, fl_276, \
                         fl_277 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = -ab_y[k] * fk_216[k]
                       + fl_271[k];

            t_361[k] = -ab_y[k] * fk_217[k]
                       + fl_273[k];

            t_362[k] = -ab_y[k] * fk_218[k]
                       + fl_274[k];

            t_363[k] = -ab_y[k] * fk_219[k]
                       + fl_276[k];

            t_364[k] = -ab_y[k] * fk_220[k]
                       + fl_277[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_y, fk_221, fk_222, fk_223, \
                         fk_224, fk_225, fl_278, fl_280, fl_281, fl_282, \
                         fl_283 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = -ab_y[k] * fk_221[k]
                       + fl_278[k];

            t_366[k] = -ab_y[k] * fk_222[k]
                       + fl_280[k];

            t_367[k] = -ab_y[k] * fk_223[k]
                       + fl_281[k];

            t_368[k] = -ab_y[k] * fk_224[k]
                       + fl_282[k];

            t_369[k] = -ab_y[k] * fk_225[k]
                       + fl_283[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, fk_226, fk_227, fk_228, \
                         fk_229, fk_230, fl_285, fl_286, fl_287, fl_288, \
                         fl_289 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = -ab_y[k] * fk_226[k]
                       + fl_285[k];

            t_371[k] = -ab_y[k] * fk_227[k]
                       + fl_286[k];

            t_372[k] = -ab_y[k] * fk_228[k]
                       + fl_287[k];

            t_373[k] = -ab_y[k] * fk_229[k]
                       + fl_288[k];

            t_374[k] = -ab_y[k] * fk_230[k]
                       + fl_289[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_y, fk_231, fk_232, fk_233, \
                         fk_234, fk_235, fl_291, fl_292, fl_293, fl_294, \
                         fl_295 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = -ab_y[k] * fk_231[k]
                       + fl_291[k];

            t_376[k] = -ab_y[k] * fk_232[k]
                       + fl_292[k];

            t_377[k] = -ab_y[k] * fk_233[k]
                       + fl_293[k];

            t_378[k] = -ab_y[k] * fk_234[k]
                       + fl_294[k];

            t_379[k] = -ab_y[k] * fk_235[k]
                       + fl_295[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_y, fk_236, fk_237, fk_238, \
                         fk_239, fk_240, fl_296, fl_298, fl_299, fl_300, \
                         fl_301 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = -ab_y[k] * fk_236[k]
                       + fl_296[k];

            t_381[k] = -ab_y[k] * fk_237[k]
                       + fl_298[k];

            t_382[k] = -ab_y[k] * fk_238[k]
                       + fl_299[k];

            t_383[k] = -ab_y[k] * fk_239[k]
                       + fl_300[k];

            t_384[k] = -ab_y[k] * fk_240[k]
                       + fl_301[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, fk_241, fk_242, fk_243, \
                         fk_244, fk_245, fl_302, fl_303, fl_304, fl_306, \
                         fl_307 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = -ab_y[k] * fk_241[k]
                       + fl_302[k];

            t_386[k] = -ab_y[k] * fk_242[k]
                       + fl_303[k];

            t_387[k] = -ab_y[k] * fk_243[k]
                       + fl_304[k];

            t_388[k] = -ab_y[k] * fk_244[k]
                       + fl_306[k];

            t_389[k] = -ab_y[k] * fk_245[k]
                       + fl_307[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_y, fk_246, fk_247, fk_248, \
                         fk_249, fk_250, fl_308, fl_309, fl_310, fl_311, \
                         fl_312 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = -ab_y[k] * fk_246[k]
                       + fl_308[k];

            t_391[k] = -ab_y[k] * fk_247[k]
                       + fl_309[k];

            t_392[k] = -ab_y[k] * fk_248[k]
                       + fl_310[k];

            t_393[k] = -ab_y[k] * fk_249[k]
                       + fl_311[k];

            t_394[k] = -ab_y[k] * fk_250[k]
                       + fl_312[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_y, fk_251, fk_252, fk_253, \
                         fk_254, fk_255, fl_313, fl_316, fl_318, fl_319, \
                         fl_321 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = -ab_y[k] * fk_251[k]
                       + fl_313[k];

            t_396[k] = -ab_y[k] * fk_252[k]
                       + fl_316[k];

            t_397[k] = -ab_y[k] * fk_253[k]
                       + fl_318[k];

            t_398[k] = -ab_y[k] * fk_254[k]
                       + fl_319[k];

            t_399[k] = -ab_y[k] * fk_255[k]
                       + fl_321[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, fk_256, fk_257, fk_258, \
                         fk_259, fk_260, fl_322, fl_323, fl_325, fl_326, \
                         fl_327 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = -ab_y[k] * fk_256[k]
                       + fl_322[k];

            t_401[k] = -ab_y[k] * fk_257[k]
                       + fl_323[k];

            t_402[k] = -ab_y[k] * fk_258[k]
                       + fl_325[k];

            t_403[k] = -ab_y[k] * fk_259[k]
                       + fl_326[k];

            t_404[k] = -ab_y[k] * fk_260[k]
                       + fl_327[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_y, fk_261, fk_262, fk_263, \
                         fk_264, fk_265, fl_328, fl_330, fl_331, fl_332, \
                         fl_333 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = -ab_y[k] * fk_261[k]
                       + fl_328[k];

            t_406[k] = -ab_y[k] * fk_262[k]
                       + fl_330[k];

            t_407[k] = -ab_y[k] * fk_263[k]
                       + fl_331[k];

            t_408[k] = -ab_y[k] * fk_264[k]
                       + fl_332[k];

            t_409[k] = -ab_y[k] * fk_265[k]
                       + fl_333[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_y, fk_266, fk_267, fk_268, \
                         fk_269, fk_270, fl_334, fl_336, fl_337, fl_338, \
                         fl_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = -ab_y[k] * fk_266[k]
                       + fl_334[k];

            t_411[k] = -ab_y[k] * fk_267[k]
                       + fl_336[k];

            t_412[k] = -ab_y[k] * fk_268[k]
                       + fl_337[k];

            t_413[k] = -ab_y[k] * fk_269[k]
                       + fl_338[k];

            t_414[k] = -ab_y[k] * fk_270[k]
                       + fl_339[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, fk_271, fk_272, fk_273, \
                         fk_274, fk_275, fl_340, fl_341, fl_343, fl_344, \
                         fl_345 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = -ab_y[k] * fk_271[k]
                       + fl_340[k];

            t_416[k] = -ab_y[k] * fk_272[k]
                       + fl_341[k];

            t_417[k] = -ab_y[k] * fk_273[k]
                       + fl_343[k];

            t_418[k] = -ab_y[k] * fk_274[k]
                       + fl_344[k];

            t_419[k] = -ab_y[k] * fk_275[k]
                       + fl_345[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_y, fk_276, fk_277, fk_278, \
                         fk_279, fk_280, fl_346, fl_347, fl_348, fl_349, \
                         fl_351 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = -ab_y[k] * fk_276[k]
                       + fl_346[k];

            t_421[k] = -ab_y[k] * fk_277[k]
                       + fl_347[k];

            t_422[k] = -ab_y[k] * fk_278[k]
                       + fl_348[k];

            t_423[k] = -ab_y[k] * fk_279[k]
                       + fl_349[k];

            t_424[k] = -ab_y[k] * fk_280[k]
                       + fl_351[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_y, fk_281, fk_282, fk_283, \
                         fk_284, fk_285, fl_352, fl_353, fl_354, fl_355, \
                         fl_356 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = -ab_y[k] * fk_281[k]
                       + fl_352[k];

            t_426[k] = -ab_y[k] * fk_282[k]
                       + fl_353[k];

            t_427[k] = -ab_y[k] * fk_283[k]
                       + fl_354[k];

            t_428[k] = -ab_y[k] * fk_284[k]
                       + fl_355[k];

            t_429[k] = -ab_y[k] * fk_285[k]
                       + fl_356[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_y, fk_286, fk_287, fk_288, \
                         fk_289, fk_290, fl_357, fl_358, fl_361, fl_363, \
                         fl_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = -ab_y[k] * fk_286[k]
                       + fl_357[k];

            t_431[k] = -ab_y[k] * fk_287[k]
                       + fl_358[k];

            t_432[k] = -ab_y[k] * fk_288[k]
                       + fl_361[k];

            t_433[k] = -ab_y[k] * fk_289[k]
                       + fl_363[k];

            t_434[k] = -ab_y[k] * fk_290[k]
                       + fl_364[k];
        }
    }
}

static auto
compute_hrr_gk_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fk, const size_t fl, const size_t ncomps,
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fk_291 = buffer.data(fk + 291 * ncomps + c);
        const auto *fk_292 = buffer.data(fk + 292 * ncomps + c);
        const auto *fk_293 = buffer.data(fk + 293 * ncomps + c);
        const auto *fk_294 = buffer.data(fk + 294 * ncomps + c);
        const auto *fk_295 = buffer.data(fk + 295 * ncomps + c);
        const auto *fk_296 = buffer.data(fk + 296 * ncomps + c);
        const auto *fk_297 = buffer.data(fk + 297 * ncomps + c);
        const auto *fk_298 = buffer.data(fk + 298 * ncomps + c);
        const auto *fk_299 = buffer.data(fk + 299 * ncomps + c);
        const auto *fk_300 = buffer.data(fk + 300 * ncomps + c);
        const auto *fk_301 = buffer.data(fk + 301 * ncomps + c);
        const auto *fk_302 = buffer.data(fk + 302 * ncomps + c);
        const auto *fk_303 = buffer.data(fk + 303 * ncomps + c);
        const auto *fk_304 = buffer.data(fk + 304 * ncomps + c);
        const auto *fk_305 = buffer.data(fk + 305 * ncomps + c);
        const auto *fk_306 = buffer.data(fk + 306 * ncomps + c);
        const auto *fk_307 = buffer.data(fk + 307 * ncomps + c);
        const auto *fk_308 = buffer.data(fk + 308 * ncomps + c);
        const auto *fk_309 = buffer.data(fk + 309 * ncomps + c);
        const auto *fk_310 = buffer.data(fk + 310 * ncomps + c);
        const auto *fk_311 = buffer.data(fk + 311 * ncomps + c);
        const auto *fk_312 = buffer.data(fk + 312 * ncomps + c);
        const auto *fk_313 = buffer.data(fk + 313 * ncomps + c);
        const auto *fk_314 = buffer.data(fk + 314 * ncomps + c);
        const auto *fk_315 = buffer.data(fk + 315 * ncomps + c);
        const auto *fk_316 = buffer.data(fk + 316 * ncomps + c);
        const auto *fk_317 = buffer.data(fk + 317 * ncomps + c);
        const auto *fk_318 = buffer.data(fk + 318 * ncomps + c);
        const auto *fk_319 = buffer.data(fk + 319 * ncomps + c);
        const auto *fk_320 = buffer.data(fk + 320 * ncomps + c);
        const auto *fk_321 = buffer.data(fk + 321 * ncomps + c);
        const auto *fk_322 = buffer.data(fk + 322 * ncomps + c);
        const auto *fk_323 = buffer.data(fk + 323 * ncomps + c);
        const auto *fk_324 = buffer.data(fk + 324 * ncomps + c);
        const auto *fk_325 = buffer.data(fk + 325 * ncomps + c);
        const auto *fk_326 = buffer.data(fk + 326 * ncomps + c);
        const auto *fk_327 = buffer.data(fk + 327 * ncomps + c);
        const auto *fk_328 = buffer.data(fk + 328 * ncomps + c);
        const auto *fk_329 = buffer.data(fk + 329 * ncomps + c);
        const auto *fk_330 = buffer.data(fk + 330 * ncomps + c);
        const auto *fk_331 = buffer.data(fk + 331 * ncomps + c);
        const auto *fk_332 = buffer.data(fk + 332 * ncomps + c);
        const auto *fk_333 = buffer.data(fk + 333 * ncomps + c);
        const auto *fk_334 = buffer.data(fk + 334 * ncomps + c);
        const auto *fk_335 = buffer.data(fk + 335 * ncomps + c);
        const auto *fk_336 = buffer.data(fk + 336 * ncomps + c);
        const auto *fk_337 = buffer.data(fk + 337 * ncomps + c);
        const auto *fk_338 = buffer.data(fk + 338 * ncomps + c);
        const auto *fk_339 = buffer.data(fk + 339 * ncomps + c);
        const auto *fk_340 = buffer.data(fk + 340 * ncomps + c);
        const auto *fk_341 = buffer.data(fk + 341 * ncomps + c);
        const auto *fk_342 = buffer.data(fk + 342 * ncomps + c);
        const auto *fk_343 = buffer.data(fk + 343 * ncomps + c);
        const auto *fk_344 = buffer.data(fk + 344 * ncomps + c);
        const auto *fk_345 = buffer.data(fk + 345 * ncomps + c);
        const auto *fk_346 = buffer.data(fk + 346 * ncomps + c);
        const auto *fk_347 = buffer.data(fk + 347 * ncomps + c);
        const auto *fk_348 = buffer.data(fk + 348 * ncomps + c);
        const auto *fk_349 = buffer.data(fk + 349 * ncomps + c);
        const auto *fk_350 = buffer.data(fk + 350 * ncomps + c);
        const auto *fk_351 = buffer.data(fk + 351 * ncomps + c);
        const auto *fk_352 = buffer.data(fk + 352 * ncomps + c);
        const auto *fk_353 = buffer.data(fk + 353 * ncomps + c);
        const auto *fk_354 = buffer.data(fk + 354 * ncomps + c);
        const auto *fk_355 = buffer.data(fk + 355 * ncomps + c);
        const auto *fk_356 = buffer.data(fk + 356 * ncomps + c);
        const auto *fk_357 = buffer.data(fk + 357 * ncomps + c);
        const auto *fk_358 = buffer.data(fk + 358 * ncomps + c);
        const auto *fk_359 = buffer.data(fk + 359 * ncomps + c);

        const auto *fl_366 = buffer.data(fl + 366 * ncomps + c);
        const auto *fl_367 = buffer.data(fl + 367 * ncomps + c);
        const auto *fl_368 = buffer.data(fl + 368 * ncomps + c);
        const auto *fl_370 = buffer.data(fl + 370 * ncomps + c);
        const auto *fl_371 = buffer.data(fl + 371 * ncomps + c);
        const auto *fl_372 = buffer.data(fl + 372 * ncomps + c);
        const auto *fl_373 = buffer.data(fl + 373 * ncomps + c);
        const auto *fl_375 = buffer.data(fl + 375 * ncomps + c);
        const auto *fl_376 = buffer.data(fl + 376 * ncomps + c);
        const auto *fl_377 = buffer.data(fl + 377 * ncomps + c);
        const auto *fl_378 = buffer.data(fl + 378 * ncomps + c);
        const auto *fl_379 = buffer.data(fl + 379 * ncomps + c);
        const auto *fl_381 = buffer.data(fl + 381 * ncomps + c);
        const auto *fl_382 = buffer.data(fl + 382 * ncomps + c);
        const auto *fl_383 = buffer.data(fl + 383 * ncomps + c);
        const auto *fl_384 = buffer.data(fl + 384 * ncomps + c);
        const auto *fl_385 = buffer.data(fl + 385 * ncomps + c);
        const auto *fl_386 = buffer.data(fl + 386 * ncomps + c);
        const auto *fl_388 = buffer.data(fl + 388 * ncomps + c);
        const auto *fl_389 = buffer.data(fl + 389 * ncomps + c);
        const auto *fl_390 = buffer.data(fl + 390 * ncomps + c);
        const auto *fl_391 = buffer.data(fl + 391 * ncomps + c);
        const auto *fl_392 = buffer.data(fl + 392 * ncomps + c);
        const auto *fl_393 = buffer.data(fl + 393 * ncomps + c);
        const auto *fl_394 = buffer.data(fl + 394 * ncomps + c);
        const auto *fl_396 = buffer.data(fl + 396 * ncomps + c);
        const auto *fl_397 = buffer.data(fl + 397 * ncomps + c);
        const auto *fl_398 = buffer.data(fl + 398 * ncomps + c);
        const auto *fl_399 = buffer.data(fl + 399 * ncomps + c);
        const auto *fl_400 = buffer.data(fl + 400 * ncomps + c);
        const auto *fl_401 = buffer.data(fl + 401 * ncomps + c);
        const auto *fl_402 = buffer.data(fl + 402 * ncomps + c);
        const auto *fl_403 = buffer.data(fl + 403 * ncomps + c);
        const auto *fl_406 = buffer.data(fl + 406 * ncomps + c);
        const auto *fl_407 = buffer.data(fl + 407 * ncomps + c);
        const auto *fl_408 = buffer.data(fl + 408 * ncomps + c);
        const auto *fl_409 = buffer.data(fl + 409 * ncomps + c);
        const auto *fl_410 = buffer.data(fl + 410 * ncomps + c);
        const auto *fl_411 = buffer.data(fl + 411 * ncomps + c);
        const auto *fl_412 = buffer.data(fl + 412 * ncomps + c);
        const auto *fl_413 = buffer.data(fl + 413 * ncomps + c);
        const auto *fl_414 = buffer.data(fl + 414 * ncomps + c);
        const auto *fl_415 = buffer.data(fl + 415 * ncomps + c);
        const auto *fl_416 = buffer.data(fl + 416 * ncomps + c);
        const auto *fl_417 = buffer.data(fl + 417 * ncomps + c);
        const auto *fl_418 = buffer.data(fl + 418 * ncomps + c);
        const auto *fl_419 = buffer.data(fl + 419 * ncomps + c);
        const auto *fl_420 = buffer.data(fl + 420 * ncomps + c);
        const auto *fl_421 = buffer.data(fl + 421 * ncomps + c);
        const auto *fl_422 = buffer.data(fl + 422 * ncomps + c);
        const auto *fl_423 = buffer.data(fl + 423 * ncomps + c);
        const auto *fl_424 = buffer.data(fl + 424 * ncomps + c);
        const auto *fl_425 = buffer.data(fl + 425 * ncomps + c);
        const auto *fl_426 = buffer.data(fl + 426 * ncomps + c);
        const auto *fl_427 = buffer.data(fl + 427 * ncomps + c);
        const auto *fl_428 = buffer.data(fl + 428 * ncomps + c);
        const auto *fl_429 = buffer.data(fl + 429 * ncomps + c);
        const auto *fl_430 = buffer.data(fl + 430 * ncomps + c);
        const auto *fl_431 = buffer.data(fl + 431 * ncomps + c);
        const auto *fl_432 = buffer.data(fl + 432 * ncomps + c);
        const auto *fl_433 = buffer.data(fl + 433 * ncomps + c);
        const auto *fl_434 = buffer.data(fl + 434 * ncomps + c);
        const auto *fl_435 = buffer.data(fl + 435 * ncomps + c);
        const auto *fl_436 = buffer.data(fl + 436 * ncomps + c);
        const auto *fl_437 = buffer.data(fl + 437 * ncomps + c);
        const auto *fl_438 = buffer.data(fl + 438 * ncomps + c);
        const auto *fl_439 = buffer.data(fl + 439 * ncomps + c);
        const auto *fl_440 = buffer.data(fl + 440 * ncomps + c);
        const auto *fl_441 = buffer.data(fl + 441 * ncomps + c);
        const auto *fl_442 = buffer.data(fl + 442 * ncomps + c);
        const auto *fl_443 = buffer.data(fl + 443 * ncomps + c);
        const auto *fl_444 = buffer.data(fl + 444 * ncomps + c);
        const auto *fl_445 = buffer.data(fl + 445 * ncomps + c);
        const auto *fl_446 = buffer.data(fl + 446 * ncomps + c);
        const auto *fl_447 = buffer.data(fl + 447 * ncomps + c);
        const auto *fl_448 = buffer.data(fl + 448 * ncomps + c);
        const auto *fl_449 = buffer.data(fl + 449 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, fk_291, fk_292, fk_293, \
                         fk_294, fk_295, fl_366, fl_367, fl_368, fl_370, \
                         fl_371 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = -ab_y[k] * fk_291[k]
                       + fl_366[k];

            t_436[k] = -ab_y[k] * fk_292[k]
                       + fl_367[k];

            t_437[k] = -ab_y[k] * fk_293[k]
                       + fl_368[k];

            t_438[k] = -ab_y[k] * fk_294[k]
                       + fl_370[k];

            t_439[k] = -ab_y[k] * fk_295[k]
                       + fl_371[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_y, fk_296, fk_297, fk_298, \
                         fk_299, fk_300, fl_372, fl_373, fl_375, fl_376, \
                         fl_377 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = -ab_y[k] * fk_296[k]
                       + fl_372[k];

            t_441[k] = -ab_y[k] * fk_297[k]
                       + fl_373[k];

            t_442[k] = -ab_y[k] * fk_298[k]
                       + fl_375[k];

            t_443[k] = -ab_y[k] * fk_299[k]
                       + fl_376[k];

            t_444[k] = -ab_y[k] * fk_300[k]
                       + fl_377[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_y, fk_301, fk_302, fk_303, \
                         fk_304, fk_305, fl_378, fl_379, fl_381, fl_382, \
                         fl_383 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = -ab_y[k] * fk_301[k]
                       + fl_378[k];

            t_446[k] = -ab_y[k] * fk_302[k]
                       + fl_379[k];

            t_447[k] = -ab_y[k] * fk_303[k]
                       + fl_381[k];

            t_448[k] = -ab_y[k] * fk_304[k]
                       + fl_382[k];

            t_449[k] = -ab_y[k] * fk_305[k]
                       + fl_383[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_y, fk_306, fk_307, fk_308, \
                         fk_309, fk_310, fl_384, fl_385, fl_386, fl_388, \
                         fl_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = -ab_y[k] * fk_306[k]
                       + fl_384[k];

            t_451[k] = -ab_y[k] * fk_307[k]
                       + fl_385[k];

            t_452[k] = -ab_y[k] * fk_308[k]
                       + fl_386[k];

            t_453[k] = -ab_y[k] * fk_309[k]
                       + fl_388[k];

            t_454[k] = -ab_y[k] * fk_310[k]
                       + fl_389[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_y, fk_311, fk_312, fk_313, \
                         fk_314, fk_315, fl_390, fl_391, fl_392, fl_393, \
                         fl_394 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = -ab_y[k] * fk_311[k]
                       + fl_390[k];

            t_456[k] = -ab_y[k] * fk_312[k]
                       + fl_391[k];

            t_457[k] = -ab_y[k] * fk_313[k]
                       + fl_392[k];

            t_458[k] = -ab_y[k] * fk_314[k]
                       + fl_393[k];

            t_459[k] = -ab_y[k] * fk_315[k]
                       + fl_394[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, fk_316, fk_317, fk_318, \
                         fk_319, fk_320, fl_396, fl_397, fl_398, fl_399, \
                         fl_400 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = -ab_y[k] * fk_316[k]
                       + fl_396[k];

            t_461[k] = -ab_y[k] * fk_317[k]
                       + fl_397[k];

            t_462[k] = -ab_y[k] * fk_318[k]
                       + fl_398[k];

            t_463[k] = -ab_y[k] * fk_319[k]
                       + fl_399[k];

            t_464[k] = -ab_y[k] * fk_320[k]
                       + fl_400[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_y, fk_321, fk_322, fk_323, \
                         fk_324, fk_325, fl_401, fl_402, fl_403, fl_406, \
                         fl_408 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = -ab_y[k] * fk_321[k]
                       + fl_401[k];

            t_466[k] = -ab_y[k] * fk_322[k]
                       + fl_402[k];

            t_467[k] = -ab_y[k] * fk_323[k]
                       + fl_403[k];

            t_468[k] = -ab_y[k] * fk_324[k]
                       + fl_406[k];

            t_469[k] = -ab_y[k] * fk_325[k]
                       + fl_408[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_y, fk_326, fk_327, fk_328, \
                         fk_329, fk_330, fl_409, fl_411, fl_412, fl_413, \
                         fl_415 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = -ab_y[k] * fk_326[k]
                       + fl_409[k];

            t_471[k] = -ab_y[k] * fk_327[k]
                       + fl_411[k];

            t_472[k] = -ab_y[k] * fk_328[k]
                       + fl_412[k];

            t_473[k] = -ab_y[k] * fk_329[k]
                       + fl_413[k];

            t_474[k] = -ab_y[k] * fk_330[k]
                       + fl_415[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, fk_331, fk_332, fk_333, \
                         fk_334, fk_335, fl_416, fl_417, fl_418, fl_420, \
                         fl_421 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = -ab_y[k] * fk_331[k]
                       + fl_416[k];

            t_476[k] = -ab_y[k] * fk_332[k]
                       + fl_417[k];

            t_477[k] = -ab_y[k] * fk_333[k]
                       + fl_418[k];

            t_478[k] = -ab_y[k] * fk_334[k]
                       + fl_420[k];

            t_479[k] = -ab_y[k] * fk_335[k]
                       + fl_421[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_y, fk_336, fk_337, fk_338, \
                         fk_339, fk_340, fl_422, fl_423, fl_424, fl_426, \
                         fl_427 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = -ab_y[k] * fk_336[k]
                       + fl_422[k];

            t_481[k] = -ab_y[k] * fk_337[k]
                       + fl_423[k];

            t_482[k] = -ab_y[k] * fk_338[k]
                       + fl_424[k];

            t_483[k] = -ab_y[k] * fk_339[k]
                       + fl_426[k];

            t_484[k] = -ab_y[k] * fk_340[k]
                       + fl_427[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_y, fk_341, fk_342, fk_343, \
                         fk_344, fk_345, fl_428, fl_429, fl_430, fl_431, \
                         fl_433 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = -ab_y[k] * fk_341[k]
                       + fl_428[k];

            t_486[k] = -ab_y[k] * fk_342[k]
                       + fl_429[k];

            t_487[k] = -ab_y[k] * fk_343[k]
                       + fl_430[k];

            t_488[k] = -ab_y[k] * fk_344[k]
                       + fl_431[k];

            t_489[k] = -ab_y[k] * fk_345[k]
                       + fl_433[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, fk_346, fk_347, fk_348, \
                         fk_349, fk_350, fl_434, fl_435, fl_436, fl_437, \
                         fl_438 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = -ab_y[k] * fk_346[k]
                       + fl_434[k];

            t_491[k] = -ab_y[k] * fk_347[k]
                       + fl_435[k];

            t_492[k] = -ab_y[k] * fk_348[k]
                       + fl_436[k];

            t_493[k] = -ab_y[k] * fk_349[k]
                       + fl_437[k];

            t_494[k] = -ab_y[k] * fk_350[k]
                       + fl_438[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_y, fk_351, fk_352, fk_353, \
                         fk_354, fk_355, fl_439, fl_441, fl_442, fl_443, \
                         fl_444 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = -ab_y[k] * fk_351[k]
                       + fl_439[k];

            t_496[k] = -ab_y[k] * fk_352[k]
                       + fl_441[k];

            t_497[k] = -ab_y[k] * fk_353[k]
                       + fl_442[k];

            t_498[k] = -ab_y[k] * fk_354[k]
                       + fl_443[k];

            t_499[k] = -ab_y[k] * fk_355[k]
                       + fl_444[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, ab_y, fk_356, fk_357, fk_358, fk_359, \
                         fl_445, fl_446, fl_447, fl_448 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = -ab_y[k] * fk_356[k]
                       + fl_445[k];

            t_501[k] = -ab_y[k] * fk_357[k]
                       + fl_446[k];

            t_502[k] = -ab_y[k] * fk_358[k]
                       + fl_447[k];

            t_503[k] = -ab_y[k] * fk_359[k]
                       + fl_448[k];
        }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, ab_z, fk_324, fk_325, fk_326, \
                         fk_327, fk_328, fl_407, fl_409, fl_410, fl_412, \
                         fl_413 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_504[k] = -ab_z[k] * fk_324[k]
                       + fl_407[k];

            t_505[k] = -ab_z[k] * fk_325[k]
                       + fl_409[k];

            t_506[k] = -ab_z[k] * fk_326[k]
                       + fl_410[k];

            t_507[k] = -ab_z[k] * fk_327[k]
                       + fl_412[k];

            t_508[k] = -ab_z[k] * fk_328[k]
                       + fl_413[k];
        }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, ab_z, fk_329, fk_330, fk_331, \
                         fk_332, fk_333, fl_414, fl_416, fl_417, fl_418, \
                         fl_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_509[k] = -ab_z[k] * fk_329[k]
                       + fl_414[k];

            t_510[k] = -ab_z[k] * fk_330[k]
                       + fl_416[k];

            t_511[k] = -ab_z[k] * fk_331[k]
                       + fl_417[k];

            t_512[k] = -ab_z[k] * fk_332[k]
                       + fl_418[k];

            t_513[k] = -ab_z[k] * fk_333[k]
                       + fl_419[k];
        }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, ab_z, fk_334, fk_335, fk_336, \
                         fk_337, fk_338, fl_421, fl_422, fl_423, fl_424, \
                         fl_425 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_514[k] = -ab_z[k] * fk_334[k]
                       + fl_421[k];

            t_515[k] = -ab_z[k] * fk_335[k]
                       + fl_422[k];

            t_516[k] = -ab_z[k] * fk_336[k]
                       + fl_423[k];

            t_517[k] = -ab_z[k] * fk_337[k]
                       + fl_424[k];

            t_518[k] = -ab_z[k] * fk_338[k]
                       + fl_425[k];
        }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, ab_z, fk_339, fk_340, fk_341, \
                         fk_342, fk_343, fl_427, fl_428, fl_429, fl_430, \
                         fl_431 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_519[k] = -ab_z[k] * fk_339[k]
                       + fl_427[k];

            t_520[k] = -ab_z[k] * fk_340[k]
                       + fl_428[k];

            t_521[k] = -ab_z[k] * fk_341[k]
                       + fl_429[k];

            t_522[k] = -ab_z[k] * fk_342[k]
                       + fl_430[k];

            t_523[k] = -ab_z[k] * fk_343[k]
                       + fl_431[k];
        }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, ab_z, fk_344, fk_345, fk_346, \
                         fk_347, fk_348, fl_432, fl_434, fl_435, fl_436, \
                         fl_437 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_524[k] = -ab_z[k] * fk_344[k]
                       + fl_432[k];

            t_525[k] = -ab_z[k] * fk_345[k]
                       + fl_434[k];

            t_526[k] = -ab_z[k] * fk_346[k]
                       + fl_435[k];

            t_527[k] = -ab_z[k] * fk_347[k]
                       + fl_436[k];

            t_528[k] = -ab_z[k] * fk_348[k]
                       + fl_437[k];
        }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, ab_z, fk_349, fk_350, fk_351, \
                         fk_352, fk_353, fl_438, fl_439, fl_440, fl_442, \
                         fl_443 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_529[k] = -ab_z[k] * fk_349[k]
                       + fl_438[k];

            t_530[k] = -ab_z[k] * fk_350[k]
                       + fl_439[k];

            t_531[k] = -ab_z[k] * fk_351[k]
                       + fl_440[k];

            t_532[k] = -ab_z[k] * fk_352[k]
                       + fl_442[k];

            t_533[k] = -ab_z[k] * fk_353[k]
                       + fl_443[k];
        }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, ab_z, fk_354, fk_355, fk_356, \
                         fk_357, fk_358, fl_444, fl_445, fl_446, fl_447, \
                         fl_448 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_534[k] = -ab_z[k] * fk_354[k]
                       + fl_444[k];

            t_535[k] = -ab_z[k] * fk_355[k]
                       + fl_445[k];

            t_536[k] = -ab_z[k] * fk_356[k]
                       + fl_446[k];

            t_537[k] = -ab_z[k] * fk_357[k]
                       + fl_447[k];

            t_538[k] = -ab_z[k] * fk_358[k]
                       + fl_448[k];
        }

#pragma omp simd aligned(t_539, ab_z, fk_359, fl_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_539[k] = -ab_z[k] * fk_359[k]
                       + fl_449[k];
        }
    }
}

auto
compute_hrr_gk(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fk, const size_t fl, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_gk_piece0(buffer, coordinates, target, fk, fl, ncomps, nmax);

    compute_hrr_gk_piece1(buffer, coordinates, target, fk, fl, ncomps, nmax);

    compute_hrr_gk_piece2(buffer, coordinates, target, fk, fl, ncomps, nmax);

    compute_hrr_gk_piece3(buffer, coordinates, target, fk, fl, ncomps, nmax);
}

}  // namespace simdtrf
