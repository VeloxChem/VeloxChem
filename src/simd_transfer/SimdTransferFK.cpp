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


#include "SimdTransferFK.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_fk_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dk, const size_t dl, const size_t ncomps,
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

        const auto *dk_0 = buffer.data(dk + 0 * ncomps + c);
        const auto *dk_1 = buffer.data(dk + 1 * ncomps + c);
        const auto *dk_2 = buffer.data(dk + 2 * ncomps + c);
        const auto *dk_3 = buffer.data(dk + 3 * ncomps + c);
        const auto *dk_4 = buffer.data(dk + 4 * ncomps + c);
        const auto *dk_5 = buffer.data(dk + 5 * ncomps + c);
        const auto *dk_6 = buffer.data(dk + 6 * ncomps + c);
        const auto *dk_7 = buffer.data(dk + 7 * ncomps + c);
        const auto *dk_8 = buffer.data(dk + 8 * ncomps + c);
        const auto *dk_9 = buffer.data(dk + 9 * ncomps + c);
        const auto *dk_10 = buffer.data(dk + 10 * ncomps + c);
        const auto *dk_11 = buffer.data(dk + 11 * ncomps + c);
        const auto *dk_12 = buffer.data(dk + 12 * ncomps + c);
        const auto *dk_13 = buffer.data(dk + 13 * ncomps + c);
        const auto *dk_14 = buffer.data(dk + 14 * ncomps + c);
        const auto *dk_15 = buffer.data(dk + 15 * ncomps + c);
        const auto *dk_16 = buffer.data(dk + 16 * ncomps + c);
        const auto *dk_17 = buffer.data(dk + 17 * ncomps + c);
        const auto *dk_18 = buffer.data(dk + 18 * ncomps + c);
        const auto *dk_19 = buffer.data(dk + 19 * ncomps + c);
        const auto *dk_20 = buffer.data(dk + 20 * ncomps + c);
        const auto *dk_21 = buffer.data(dk + 21 * ncomps + c);
        const auto *dk_22 = buffer.data(dk + 22 * ncomps + c);
        const auto *dk_23 = buffer.data(dk + 23 * ncomps + c);
        const auto *dk_24 = buffer.data(dk + 24 * ncomps + c);
        const auto *dk_25 = buffer.data(dk + 25 * ncomps + c);
        const auto *dk_26 = buffer.data(dk + 26 * ncomps + c);
        const auto *dk_27 = buffer.data(dk + 27 * ncomps + c);
        const auto *dk_28 = buffer.data(dk + 28 * ncomps + c);
        const auto *dk_29 = buffer.data(dk + 29 * ncomps + c);
        const auto *dk_30 = buffer.data(dk + 30 * ncomps + c);
        const auto *dk_31 = buffer.data(dk + 31 * ncomps + c);
        const auto *dk_32 = buffer.data(dk + 32 * ncomps + c);
        const auto *dk_33 = buffer.data(dk + 33 * ncomps + c);
        const auto *dk_34 = buffer.data(dk + 34 * ncomps + c);
        const auto *dk_35 = buffer.data(dk + 35 * ncomps + c);
        const auto *dk_36 = buffer.data(dk + 36 * ncomps + c);
        const auto *dk_37 = buffer.data(dk + 37 * ncomps + c);
        const auto *dk_38 = buffer.data(dk + 38 * ncomps + c);
        const auto *dk_39 = buffer.data(dk + 39 * ncomps + c);
        const auto *dk_40 = buffer.data(dk + 40 * ncomps + c);
        const auto *dk_41 = buffer.data(dk + 41 * ncomps + c);
        const auto *dk_42 = buffer.data(dk + 42 * ncomps + c);
        const auto *dk_43 = buffer.data(dk + 43 * ncomps + c);
        const auto *dk_44 = buffer.data(dk + 44 * ncomps + c);
        const auto *dk_45 = buffer.data(dk + 45 * ncomps + c);
        const auto *dk_46 = buffer.data(dk + 46 * ncomps + c);
        const auto *dk_47 = buffer.data(dk + 47 * ncomps + c);
        const auto *dk_48 = buffer.data(dk + 48 * ncomps + c);
        const auto *dk_49 = buffer.data(dk + 49 * ncomps + c);
        const auto *dk_50 = buffer.data(dk + 50 * ncomps + c);
        const auto *dk_51 = buffer.data(dk + 51 * ncomps + c);
        const auto *dk_52 = buffer.data(dk + 52 * ncomps + c);
        const auto *dk_53 = buffer.data(dk + 53 * ncomps + c);
        const auto *dk_54 = buffer.data(dk + 54 * ncomps + c);
        const auto *dk_55 = buffer.data(dk + 55 * ncomps + c);
        const auto *dk_56 = buffer.data(dk + 56 * ncomps + c);
        const auto *dk_57 = buffer.data(dk + 57 * ncomps + c);
        const auto *dk_58 = buffer.data(dk + 58 * ncomps + c);
        const auto *dk_59 = buffer.data(dk + 59 * ncomps + c);
        const auto *dk_60 = buffer.data(dk + 60 * ncomps + c);
        const auto *dk_61 = buffer.data(dk + 61 * ncomps + c);
        const auto *dk_62 = buffer.data(dk + 62 * ncomps + c);
        const auto *dk_63 = buffer.data(dk + 63 * ncomps + c);
        const auto *dk_64 = buffer.data(dk + 64 * ncomps + c);
        const auto *dk_65 = buffer.data(dk + 65 * ncomps + c);
        const auto *dk_66 = buffer.data(dk + 66 * ncomps + c);
        const auto *dk_67 = buffer.data(dk + 67 * ncomps + c);
        const auto *dk_68 = buffer.data(dk + 68 * ncomps + c);
        const auto *dk_69 = buffer.data(dk + 69 * ncomps + c);
        const auto *dk_70 = buffer.data(dk + 70 * ncomps + c);
        const auto *dk_71 = buffer.data(dk + 71 * ncomps + c);
        const auto *dk_72 = buffer.data(dk + 72 * ncomps + c);
        const auto *dk_73 = buffer.data(dk + 73 * ncomps + c);
        const auto *dk_74 = buffer.data(dk + 74 * ncomps + c);
        const auto *dk_75 = buffer.data(dk + 75 * ncomps + c);
        const auto *dk_76 = buffer.data(dk + 76 * ncomps + c);
        const auto *dk_77 = buffer.data(dk + 77 * ncomps + c);
        const auto *dk_78 = buffer.data(dk + 78 * ncomps + c);
        const auto *dk_79 = buffer.data(dk + 79 * ncomps + c);
        const auto *dk_80 = buffer.data(dk + 80 * ncomps + c);
        const auto *dk_81 = buffer.data(dk + 81 * ncomps + c);
        const auto *dk_82 = buffer.data(dk + 82 * ncomps + c);
        const auto *dk_83 = buffer.data(dk + 83 * ncomps + c);
        const auto *dk_84 = buffer.data(dk + 84 * ncomps + c);
        const auto *dk_85 = buffer.data(dk + 85 * ncomps + c);
        const auto *dk_86 = buffer.data(dk + 86 * ncomps + c);
        const auto *dk_87 = buffer.data(dk + 87 * ncomps + c);
        const auto *dk_88 = buffer.data(dk + 88 * ncomps + c);
        const auto *dk_89 = buffer.data(dk + 89 * ncomps + c);
        const auto *dk_90 = buffer.data(dk + 90 * ncomps + c);
        const auto *dk_91 = buffer.data(dk + 91 * ncomps + c);
        const auto *dk_92 = buffer.data(dk + 92 * ncomps + c);
        const auto *dk_93 = buffer.data(dk + 93 * ncomps + c);
        const auto *dk_94 = buffer.data(dk + 94 * ncomps + c);
        const auto *dk_95 = buffer.data(dk + 95 * ncomps + c);
        const auto *dk_96 = buffer.data(dk + 96 * ncomps + c);
        const auto *dk_97 = buffer.data(dk + 97 * ncomps + c);
        const auto *dk_98 = buffer.data(dk + 98 * ncomps + c);
        const auto *dk_99 = buffer.data(dk + 99 * ncomps + c);
        const auto *dk_100 = buffer.data(dk + 100 * ncomps + c);
        const auto *dk_101 = buffer.data(dk + 101 * ncomps + c);
        const auto *dk_102 = buffer.data(dk + 102 * ncomps + c);
        const auto *dk_103 = buffer.data(dk + 103 * ncomps + c);
        const auto *dk_104 = buffer.data(dk + 104 * ncomps + c);
        const auto *dk_105 = buffer.data(dk + 105 * ncomps + c);
        const auto *dk_106 = buffer.data(dk + 106 * ncomps + c);
        const auto *dk_107 = buffer.data(dk + 107 * ncomps + c);
        const auto *dk_108 = buffer.data(dk + 108 * ncomps + c);
        const auto *dk_109 = buffer.data(dk + 109 * ncomps + c);
        const auto *dk_110 = buffer.data(dk + 110 * ncomps + c);
        const auto *dk_111 = buffer.data(dk + 111 * ncomps + c);
        const auto *dk_112 = buffer.data(dk + 112 * ncomps + c);
        const auto *dk_113 = buffer.data(dk + 113 * ncomps + c);
        const auto *dk_114 = buffer.data(dk + 114 * ncomps + c);
        const auto *dk_115 = buffer.data(dk + 115 * ncomps + c);
        const auto *dk_116 = buffer.data(dk + 116 * ncomps + c);
        const auto *dk_117 = buffer.data(dk + 117 * ncomps + c);
        const auto *dk_118 = buffer.data(dk + 118 * ncomps + c);
        const auto *dk_119 = buffer.data(dk + 119 * ncomps + c);
        const auto *dk_120 = buffer.data(dk + 120 * ncomps + c);
        const auto *dk_121 = buffer.data(dk + 121 * ncomps + c);
        const auto *dk_122 = buffer.data(dk + 122 * ncomps + c);
        const auto *dk_123 = buffer.data(dk + 123 * ncomps + c);
        const auto *dk_124 = buffer.data(dk + 124 * ncomps + c);
        const auto *dk_125 = buffer.data(dk + 125 * ncomps + c);
        const auto *dk_126 = buffer.data(dk + 126 * ncomps + c);
        const auto *dk_127 = buffer.data(dk + 127 * ncomps + c);
        const auto *dk_128 = buffer.data(dk + 128 * ncomps + c);
        const auto *dk_129 = buffer.data(dk + 129 * ncomps + c);
        const auto *dk_130 = buffer.data(dk + 130 * ncomps + c);
        const auto *dk_131 = buffer.data(dk + 131 * ncomps + c);
        const auto *dk_132 = buffer.data(dk + 132 * ncomps + c);
        const auto *dk_133 = buffer.data(dk + 133 * ncomps + c);
        const auto *dk_134 = buffer.data(dk + 134 * ncomps + c);
        const auto *dk_135 = buffer.data(dk + 135 * ncomps + c);
        const auto *dk_136 = buffer.data(dk + 136 * ncomps + c);
        const auto *dk_137 = buffer.data(dk + 137 * ncomps + c);
        const auto *dk_138 = buffer.data(dk + 138 * ncomps + c);
        const auto *dk_139 = buffer.data(dk + 139 * ncomps + c);
        const auto *dk_140 = buffer.data(dk + 140 * ncomps + c);
        const auto *dk_141 = buffer.data(dk + 141 * ncomps + c);
        const auto *dk_142 = buffer.data(dk + 142 * ncomps + c);
        const auto *dk_143 = buffer.data(dk + 143 * ncomps + c);
        const auto *dk_144 = buffer.data(dk + 144 * ncomps + c);

        const auto *dl_0 = buffer.data(dl + 0 * ncomps + c);
        const auto *dl_1 = buffer.data(dl + 1 * ncomps + c);
        const auto *dl_2 = buffer.data(dl + 2 * ncomps + c);
        const auto *dl_3 = buffer.data(dl + 3 * ncomps + c);
        const auto *dl_4 = buffer.data(dl + 4 * ncomps + c);
        const auto *dl_5 = buffer.data(dl + 5 * ncomps + c);
        const auto *dl_6 = buffer.data(dl + 6 * ncomps + c);
        const auto *dl_7 = buffer.data(dl + 7 * ncomps + c);
        const auto *dl_8 = buffer.data(dl + 8 * ncomps + c);
        const auto *dl_9 = buffer.data(dl + 9 * ncomps + c);
        const auto *dl_10 = buffer.data(dl + 10 * ncomps + c);
        const auto *dl_11 = buffer.data(dl + 11 * ncomps + c);
        const auto *dl_12 = buffer.data(dl + 12 * ncomps + c);
        const auto *dl_13 = buffer.data(dl + 13 * ncomps + c);
        const auto *dl_14 = buffer.data(dl + 14 * ncomps + c);
        const auto *dl_15 = buffer.data(dl + 15 * ncomps + c);
        const auto *dl_16 = buffer.data(dl + 16 * ncomps + c);
        const auto *dl_17 = buffer.data(dl + 17 * ncomps + c);
        const auto *dl_18 = buffer.data(dl + 18 * ncomps + c);
        const auto *dl_19 = buffer.data(dl + 19 * ncomps + c);
        const auto *dl_20 = buffer.data(dl + 20 * ncomps + c);
        const auto *dl_21 = buffer.data(dl + 21 * ncomps + c);
        const auto *dl_22 = buffer.data(dl + 22 * ncomps + c);
        const auto *dl_23 = buffer.data(dl + 23 * ncomps + c);
        const auto *dl_24 = buffer.data(dl + 24 * ncomps + c);
        const auto *dl_25 = buffer.data(dl + 25 * ncomps + c);
        const auto *dl_26 = buffer.data(dl + 26 * ncomps + c);
        const auto *dl_27 = buffer.data(dl + 27 * ncomps + c);
        const auto *dl_28 = buffer.data(dl + 28 * ncomps + c);
        const auto *dl_29 = buffer.data(dl + 29 * ncomps + c);
        const auto *dl_30 = buffer.data(dl + 30 * ncomps + c);
        const auto *dl_31 = buffer.data(dl + 31 * ncomps + c);
        const auto *dl_32 = buffer.data(dl + 32 * ncomps + c);
        const auto *dl_33 = buffer.data(dl + 33 * ncomps + c);
        const auto *dl_34 = buffer.data(dl + 34 * ncomps + c);
        const auto *dl_35 = buffer.data(dl + 35 * ncomps + c);
        const auto *dl_45 = buffer.data(dl + 45 * ncomps + c);
        const auto *dl_46 = buffer.data(dl + 46 * ncomps + c);
        const auto *dl_47 = buffer.data(dl + 47 * ncomps + c);
        const auto *dl_48 = buffer.data(dl + 48 * ncomps + c);
        const auto *dl_49 = buffer.data(dl + 49 * ncomps + c);
        const auto *dl_50 = buffer.data(dl + 50 * ncomps + c);
        const auto *dl_51 = buffer.data(dl + 51 * ncomps + c);
        const auto *dl_52 = buffer.data(dl + 52 * ncomps + c);
        const auto *dl_53 = buffer.data(dl + 53 * ncomps + c);
        const auto *dl_54 = buffer.data(dl + 54 * ncomps + c);
        const auto *dl_55 = buffer.data(dl + 55 * ncomps + c);
        const auto *dl_56 = buffer.data(dl + 56 * ncomps + c);
        const auto *dl_57 = buffer.data(dl + 57 * ncomps + c);
        const auto *dl_58 = buffer.data(dl + 58 * ncomps + c);
        const auto *dl_59 = buffer.data(dl + 59 * ncomps + c);
        const auto *dl_60 = buffer.data(dl + 60 * ncomps + c);
        const auto *dl_61 = buffer.data(dl + 61 * ncomps + c);
        const auto *dl_62 = buffer.data(dl + 62 * ncomps + c);
        const auto *dl_63 = buffer.data(dl + 63 * ncomps + c);
        const auto *dl_64 = buffer.data(dl + 64 * ncomps + c);
        const auto *dl_65 = buffer.data(dl + 65 * ncomps + c);
        const auto *dl_66 = buffer.data(dl + 66 * ncomps + c);
        const auto *dl_67 = buffer.data(dl + 67 * ncomps + c);
        const auto *dl_68 = buffer.data(dl + 68 * ncomps + c);
        const auto *dl_69 = buffer.data(dl + 69 * ncomps + c);
        const auto *dl_70 = buffer.data(dl + 70 * ncomps + c);
        const auto *dl_71 = buffer.data(dl + 71 * ncomps + c);
        const auto *dl_72 = buffer.data(dl + 72 * ncomps + c);
        const auto *dl_73 = buffer.data(dl + 73 * ncomps + c);
        const auto *dl_74 = buffer.data(dl + 74 * ncomps + c);
        const auto *dl_75 = buffer.data(dl + 75 * ncomps + c);
        const auto *dl_76 = buffer.data(dl + 76 * ncomps + c);
        const auto *dl_77 = buffer.data(dl + 77 * ncomps + c);
        const auto *dl_78 = buffer.data(dl + 78 * ncomps + c);
        const auto *dl_79 = buffer.data(dl + 79 * ncomps + c);
        const auto *dl_80 = buffer.data(dl + 80 * ncomps + c);
        const auto *dl_90 = buffer.data(dl + 90 * ncomps + c);
        const auto *dl_91 = buffer.data(dl + 91 * ncomps + c);
        const auto *dl_92 = buffer.data(dl + 92 * ncomps + c);
        const auto *dl_93 = buffer.data(dl + 93 * ncomps + c);
        const auto *dl_94 = buffer.data(dl + 94 * ncomps + c);
        const auto *dl_95 = buffer.data(dl + 95 * ncomps + c);
        const auto *dl_96 = buffer.data(dl + 96 * ncomps + c);
        const auto *dl_97 = buffer.data(dl + 97 * ncomps + c);
        const auto *dl_98 = buffer.data(dl + 98 * ncomps + c);
        const auto *dl_99 = buffer.data(dl + 99 * ncomps + c);
        const auto *dl_100 = buffer.data(dl + 100 * ncomps + c);
        const auto *dl_101 = buffer.data(dl + 101 * ncomps + c);
        const auto *dl_102 = buffer.data(dl + 102 * ncomps + c);
        const auto *dl_103 = buffer.data(dl + 103 * ncomps + c);
        const auto *dl_104 = buffer.data(dl + 104 * ncomps + c);
        const auto *dl_105 = buffer.data(dl + 105 * ncomps + c);
        const auto *dl_106 = buffer.data(dl + 106 * ncomps + c);
        const auto *dl_107 = buffer.data(dl + 107 * ncomps + c);
        const auto *dl_108 = buffer.data(dl + 108 * ncomps + c);
        const auto *dl_109 = buffer.data(dl + 109 * ncomps + c);
        const auto *dl_110 = buffer.data(dl + 110 * ncomps + c);
        const auto *dl_111 = buffer.data(dl + 111 * ncomps + c);
        const auto *dl_112 = buffer.data(dl + 112 * ncomps + c);
        const auto *dl_113 = buffer.data(dl + 113 * ncomps + c);
        const auto *dl_114 = buffer.data(dl + 114 * ncomps + c);
        const auto *dl_115 = buffer.data(dl + 115 * ncomps + c);
        const auto *dl_116 = buffer.data(dl + 116 * ncomps + c);
        const auto *dl_117 = buffer.data(dl + 117 * ncomps + c);
        const auto *dl_118 = buffer.data(dl + 118 * ncomps + c);
        const auto *dl_119 = buffer.data(dl + 119 * ncomps + c);
        const auto *dl_120 = buffer.data(dl + 120 * ncomps + c);
        const auto *dl_121 = buffer.data(dl + 121 * ncomps + c);
        const auto *dl_122 = buffer.data(dl + 122 * ncomps + c);
        const auto *dl_123 = buffer.data(dl + 123 * ncomps + c);
        const auto *dl_124 = buffer.data(dl + 124 * ncomps + c);
        const auto *dl_125 = buffer.data(dl + 125 * ncomps + c);
        const auto *dl_135 = buffer.data(dl + 135 * ncomps + c);
        const auto *dl_136 = buffer.data(dl + 136 * ncomps + c);
        const auto *dl_137 = buffer.data(dl + 137 * ncomps + c);
        const auto *dl_138 = buffer.data(dl + 138 * ncomps + c);
        const auto *dl_139 = buffer.data(dl + 139 * ncomps + c);
        const auto *dl_140 = buffer.data(dl + 140 * ncomps + c);
        const auto *dl_141 = buffer.data(dl + 141 * ncomps + c);
        const auto *dl_142 = buffer.data(dl + 142 * ncomps + c);
        const auto *dl_143 = buffer.data(dl + 143 * ncomps + c);
        const auto *dl_144 = buffer.data(dl + 144 * ncomps + c);
        const auto *dl_145 = buffer.data(dl + 145 * ncomps + c);
        const auto *dl_146 = buffer.data(dl + 146 * ncomps + c);
        const auto *dl_147 = buffer.data(dl + 147 * ncomps + c);
        const auto *dl_148 = buffer.data(dl + 148 * ncomps + c);
        const auto *dl_149 = buffer.data(dl + 149 * ncomps + c);
        const auto *dl_150 = buffer.data(dl + 150 * ncomps + c);
        const auto *dl_151 = buffer.data(dl + 151 * ncomps + c);
        const auto *dl_152 = buffer.data(dl + 152 * ncomps + c);
        const auto *dl_153 = buffer.data(dl + 153 * ncomps + c);
        const auto *dl_154 = buffer.data(dl + 154 * ncomps + c);
        const auto *dl_155 = buffer.data(dl + 155 * ncomps + c);
        const auto *dl_156 = buffer.data(dl + 156 * ncomps + c);
        const auto *dl_157 = buffer.data(dl + 157 * ncomps + c);
        const auto *dl_158 = buffer.data(dl + 158 * ncomps + c);
        const auto *dl_159 = buffer.data(dl + 159 * ncomps + c);
        const auto *dl_160 = buffer.data(dl + 160 * ncomps + c);
        const auto *dl_161 = buffer.data(dl + 161 * ncomps + c);
        const auto *dl_162 = buffer.data(dl + 162 * ncomps + c);
        const auto *dl_163 = buffer.data(dl + 163 * ncomps + c);
        const auto *dl_164 = buffer.data(dl + 164 * ncomps + c);
        const auto *dl_165 = buffer.data(dl + 165 * ncomps + c);
        const auto *dl_166 = buffer.data(dl + 166 * ncomps + c);
        const auto *dl_167 = buffer.data(dl + 167 * ncomps + c);
        const auto *dl_168 = buffer.data(dl + 168 * ncomps + c);
        const auto *dl_169 = buffer.data(dl + 169 * ncomps + c);
        const auto *dl_170 = buffer.data(dl + 170 * ncomps + c);
        const auto *dl_180 = buffer.data(dl + 180 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dk_0, dk_1, dk_2, dk_3, dk_4, dl_0, \
                         dl_1, dl_2, dl_3, dl_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dk_0[k]
                     + dl_0[k];

            t_1[k] = -ab_x[k] * dk_1[k]
                     + dl_1[k];

            t_2[k] = -ab_x[k] * dk_2[k]
                     + dl_2[k];

            t_3[k] = -ab_x[k] * dk_3[k]
                     + dl_3[k];

            t_4[k] = -ab_x[k] * dk_4[k]
                     + dl_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dk_5, dk_6, dk_7, dk_8, dk_9, dl_5, \
                         dl_6, dl_7, dl_8, dl_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dk_5[k]
                     + dl_5[k];

            t_6[k] = -ab_x[k] * dk_6[k]
                     + dl_6[k];

            t_7[k] = -ab_x[k] * dk_7[k]
                     + dl_7[k];

            t_8[k] = -ab_x[k] * dk_8[k]
                     + dl_8[k];

            t_9[k] = -ab_x[k] * dk_9[k]
                     + dl_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dk_10, dk_11, dk_12, dk_13, \
                         dk_14, dl_10, dl_11, dl_12, dl_13, dl_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dk_10[k]
                      + dl_10[k];

            t_11[k] = -ab_x[k] * dk_11[k]
                      + dl_11[k];

            t_12[k] = -ab_x[k] * dk_12[k]
                      + dl_12[k];

            t_13[k] = -ab_x[k] * dk_13[k]
                      + dl_13[k];

            t_14[k] = -ab_x[k] * dk_14[k]
                      + dl_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dk_15, dk_16, dk_17, dk_18, \
                         dk_19, dl_15, dl_16, dl_17, dl_18, dl_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dk_15[k]
                      + dl_15[k];

            t_16[k] = -ab_x[k] * dk_16[k]
                      + dl_16[k];

            t_17[k] = -ab_x[k] * dk_17[k]
                      + dl_17[k];

            t_18[k] = -ab_x[k] * dk_18[k]
                      + dl_18[k];

            t_19[k] = -ab_x[k] * dk_19[k]
                      + dl_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dk_20, dk_21, dk_22, dk_23, \
                         dk_24, dl_20, dl_21, dl_22, dl_23, dl_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * dk_20[k]
                      + dl_20[k];

            t_21[k] = -ab_x[k] * dk_21[k]
                      + dl_21[k];

            t_22[k] = -ab_x[k] * dk_22[k]
                      + dl_22[k];

            t_23[k] = -ab_x[k] * dk_23[k]
                      + dl_23[k];

            t_24[k] = -ab_x[k] * dk_24[k]
                      + dl_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dk_25, dk_26, dk_27, dk_28, \
                         dk_29, dl_25, dl_26, dl_27, dl_28, dl_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * dk_25[k]
                      + dl_25[k];

            t_26[k] = -ab_x[k] * dk_26[k]
                      + dl_26[k];

            t_27[k] = -ab_x[k] * dk_27[k]
                      + dl_27[k];

            t_28[k] = -ab_x[k] * dk_28[k]
                      + dl_28[k];

            t_29[k] = -ab_x[k] * dk_29[k]
                      + dl_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dk_30, dk_31, dk_32, dk_33, \
                         dk_34, dl_30, dl_31, dl_32, dl_33, dl_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dk_30[k]
                      + dl_30[k];

            t_31[k] = -ab_x[k] * dk_31[k]
                      + dl_31[k];

            t_32[k] = -ab_x[k] * dk_32[k]
                      + dl_32[k];

            t_33[k] = -ab_x[k] * dk_33[k]
                      + dl_33[k];

            t_34[k] = -ab_x[k] * dk_34[k]
                      + dl_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dk_35, dk_36, dk_37, dk_38, \
                         dk_39, dl_35, dl_45, dl_46, dl_47, dl_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * dk_35[k]
                      + dl_35[k];

            t_36[k] = -ab_x[k] * dk_36[k]
                      + dl_45[k];

            t_37[k] = -ab_x[k] * dk_37[k]
                      + dl_46[k];

            t_38[k] = -ab_x[k] * dk_38[k]
                      + dl_47[k];

            t_39[k] = -ab_x[k] * dk_39[k]
                      + dl_48[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dk_40, dk_41, dk_42, dk_43, \
                         dk_44, dl_49, dl_50, dl_51, dl_52, dl_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * dk_40[k]
                      + dl_49[k];

            t_41[k] = -ab_x[k] * dk_41[k]
                      + dl_50[k];

            t_42[k] = -ab_x[k] * dk_42[k]
                      + dl_51[k];

            t_43[k] = -ab_x[k] * dk_43[k]
                      + dl_52[k];

            t_44[k] = -ab_x[k] * dk_44[k]
                      + dl_53[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dk_45, dk_46, dk_47, dk_48, \
                         dk_49, dl_54, dl_55, dl_56, dl_57, dl_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * dk_45[k]
                      + dl_54[k];

            t_46[k] = -ab_x[k] * dk_46[k]
                      + dl_55[k];

            t_47[k] = -ab_x[k] * dk_47[k]
                      + dl_56[k];

            t_48[k] = -ab_x[k] * dk_48[k]
                      + dl_57[k];

            t_49[k] = -ab_x[k] * dk_49[k]
                      + dl_58[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dk_50, dk_51, dk_52, dk_53, \
                         dk_54, dl_59, dl_60, dl_61, dl_62, dl_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * dk_50[k]
                      + dl_59[k];

            t_51[k] = -ab_x[k] * dk_51[k]
                      + dl_60[k];

            t_52[k] = -ab_x[k] * dk_52[k]
                      + dl_61[k];

            t_53[k] = -ab_x[k] * dk_53[k]
                      + dl_62[k];

            t_54[k] = -ab_x[k] * dk_54[k]
                      + dl_63[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dk_55, dk_56, dk_57, dk_58, \
                         dk_59, dl_64, dl_65, dl_66, dl_67, dl_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * dk_55[k]
                      + dl_64[k];

            t_56[k] = -ab_x[k] * dk_56[k]
                      + dl_65[k];

            t_57[k] = -ab_x[k] * dk_57[k]
                      + dl_66[k];

            t_58[k] = -ab_x[k] * dk_58[k]
                      + dl_67[k];

            t_59[k] = -ab_x[k] * dk_59[k]
                      + dl_68[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dk_60, dk_61, dk_62, dk_63, \
                         dk_64, dl_69, dl_70, dl_71, dl_72, dl_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * dk_60[k]
                      + dl_69[k];

            t_61[k] = -ab_x[k] * dk_61[k]
                      + dl_70[k];

            t_62[k] = -ab_x[k] * dk_62[k]
                      + dl_71[k];

            t_63[k] = -ab_x[k] * dk_63[k]
                      + dl_72[k];

            t_64[k] = -ab_x[k] * dk_64[k]
                      + dl_73[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dk_65, dk_66, dk_67, dk_68, \
                         dk_69, dl_74, dl_75, dl_76, dl_77, dl_78 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * dk_65[k]
                      + dl_74[k];

            t_66[k] = -ab_x[k] * dk_66[k]
                      + dl_75[k];

            t_67[k] = -ab_x[k] * dk_67[k]
                      + dl_76[k];

            t_68[k] = -ab_x[k] * dk_68[k]
                      + dl_77[k];

            t_69[k] = -ab_x[k] * dk_69[k]
                      + dl_78[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dk_70, dk_71, dk_72, dk_73, \
                         dk_74, dl_79, dl_80, dl_90, dl_91, dl_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * dk_70[k]
                      + dl_79[k];

            t_71[k] = -ab_x[k] * dk_71[k]
                      + dl_80[k];

            t_72[k] = -ab_x[k] * dk_72[k]
                      + dl_90[k];

            t_73[k] = -ab_x[k] * dk_73[k]
                      + dl_91[k];

            t_74[k] = -ab_x[k] * dk_74[k]
                      + dl_92[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dk_75, dk_76, dk_77, dk_78, \
                         dk_79, dl_93, dl_94, dl_95, dl_96, dl_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * dk_75[k]
                      + dl_93[k];

            t_76[k] = -ab_x[k] * dk_76[k]
                      + dl_94[k];

            t_77[k] = -ab_x[k] * dk_77[k]
                      + dl_95[k];

            t_78[k] = -ab_x[k] * dk_78[k]
                      + dl_96[k];

            t_79[k] = -ab_x[k] * dk_79[k]
                      + dl_97[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dk_80, dk_81, dk_82, dk_83, \
                         dk_84, dl_98, dl_99, dl_100, dl_101, dl_102 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * dk_80[k]
                      + dl_98[k];

            t_81[k] = -ab_x[k] * dk_81[k]
                      + dl_99[k];

            t_82[k] = -ab_x[k] * dk_82[k]
                      + dl_100[k];

            t_83[k] = -ab_x[k] * dk_83[k]
                      + dl_101[k];

            t_84[k] = -ab_x[k] * dk_84[k]
                      + dl_102[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dk_85, dk_86, dk_87, dk_88, \
                         dk_89, dl_103, dl_104, dl_105, dl_106, \
                         dl_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * dk_85[k]
                      + dl_103[k];

            t_86[k] = -ab_x[k] * dk_86[k]
                      + dl_104[k];

            t_87[k] = -ab_x[k] * dk_87[k]
                      + dl_105[k];

            t_88[k] = -ab_x[k] * dk_88[k]
                      + dl_106[k];

            t_89[k] = -ab_x[k] * dk_89[k]
                      + dl_107[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, dk_90, dk_91, dk_92, dk_93, \
                         dk_94, dl_108, dl_109, dl_110, dl_111, \
                         dl_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * dk_90[k]
                      + dl_108[k];

            t_91[k] = -ab_x[k] * dk_91[k]
                      + dl_109[k];

            t_92[k] = -ab_x[k] * dk_92[k]
                      + dl_110[k];

            t_93[k] = -ab_x[k] * dk_93[k]
                      + dl_111[k];

            t_94[k] = -ab_x[k] * dk_94[k]
                      + dl_112[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, dk_95, dk_96, dk_97, dk_98, \
                         dk_99, dl_113, dl_114, dl_115, dl_116, \
                         dl_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * dk_95[k]
                      + dl_113[k];

            t_96[k] = -ab_x[k] * dk_96[k]
                      + dl_114[k];

            t_97[k] = -ab_x[k] * dk_97[k]
                      + dl_115[k];

            t_98[k] = -ab_x[k] * dk_98[k]
                      + dl_116[k];

            t_99[k] = -ab_x[k] * dk_99[k]
                      + dl_117[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, dk_100, dk_101, dk_102, \
                         dk_103, dk_104, dl_118, dl_119, dl_120, dl_121, \
                         dl_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * dk_100[k]
                       + dl_118[k];

            t_101[k] = -ab_x[k] * dk_101[k]
                       + dl_119[k];

            t_102[k] = -ab_x[k] * dk_102[k]
                       + dl_120[k];

            t_103[k] = -ab_x[k] * dk_103[k]
                       + dl_121[k];

            t_104[k] = -ab_x[k] * dk_104[k]
                       + dl_122[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, dk_105, dk_106, dk_107, \
                         dk_108, dk_109, dl_123, dl_124, dl_125, dl_135, \
                         dl_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * dk_105[k]
                       + dl_123[k];

            t_106[k] = -ab_x[k] * dk_106[k]
                       + dl_124[k];

            t_107[k] = -ab_x[k] * dk_107[k]
                       + dl_125[k];

            t_108[k] = -ab_x[k] * dk_108[k]
                       + dl_135[k];

            t_109[k] = -ab_x[k] * dk_109[k]
                       + dl_136[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, dk_110, dk_111, dk_112, \
                         dk_113, dk_114, dl_137, dl_138, dl_139, dl_140, \
                         dl_141 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * dk_110[k]
                       + dl_137[k];

            t_111[k] = -ab_x[k] * dk_111[k]
                       + dl_138[k];

            t_112[k] = -ab_x[k] * dk_112[k]
                       + dl_139[k];

            t_113[k] = -ab_x[k] * dk_113[k]
                       + dl_140[k];

            t_114[k] = -ab_x[k] * dk_114[k]
                       + dl_141[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, dk_115, dk_116, dk_117, \
                         dk_118, dk_119, dl_142, dl_143, dl_144, dl_145, \
                         dl_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * dk_115[k]
                       + dl_142[k];

            t_116[k] = -ab_x[k] * dk_116[k]
                       + dl_143[k];

            t_117[k] = -ab_x[k] * dk_117[k]
                       + dl_144[k];

            t_118[k] = -ab_x[k] * dk_118[k]
                       + dl_145[k];

            t_119[k] = -ab_x[k] * dk_119[k]
                       + dl_146[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, dk_120, dk_121, dk_122, \
                         dk_123, dk_124, dl_147, dl_148, dl_149, dl_150, \
                         dl_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * dk_120[k]
                       + dl_147[k];

            t_121[k] = -ab_x[k] * dk_121[k]
                       + dl_148[k];

            t_122[k] = -ab_x[k] * dk_122[k]
                       + dl_149[k];

            t_123[k] = -ab_x[k] * dk_123[k]
                       + dl_150[k];

            t_124[k] = -ab_x[k] * dk_124[k]
                       + dl_151[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, dk_125, dk_126, dk_127, \
                         dk_128, dk_129, dl_152, dl_153, dl_154, dl_155, \
                         dl_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * dk_125[k]
                       + dl_152[k];

            t_126[k] = -ab_x[k] * dk_126[k]
                       + dl_153[k];

            t_127[k] = -ab_x[k] * dk_127[k]
                       + dl_154[k];

            t_128[k] = -ab_x[k] * dk_128[k]
                       + dl_155[k];

            t_129[k] = -ab_x[k] * dk_129[k]
                       + dl_156[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, dk_130, dk_131, dk_132, \
                         dk_133, dk_134, dl_157, dl_158, dl_159, dl_160, \
                         dl_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * dk_130[k]
                       + dl_157[k];

            t_131[k] = -ab_x[k] * dk_131[k]
                       + dl_158[k];

            t_132[k] = -ab_x[k] * dk_132[k]
                       + dl_159[k];

            t_133[k] = -ab_x[k] * dk_133[k]
                       + dl_160[k];

            t_134[k] = -ab_x[k] * dk_134[k]
                       + dl_161[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, dk_135, dk_136, dk_137, \
                         dk_138, dk_139, dl_162, dl_163, dl_164, dl_165, \
                         dl_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * dk_135[k]
                       + dl_162[k];

            t_136[k] = -ab_x[k] * dk_136[k]
                       + dl_163[k];

            t_137[k] = -ab_x[k] * dk_137[k]
                       + dl_164[k];

            t_138[k] = -ab_x[k] * dk_138[k]
                       + dl_165[k];

            t_139[k] = -ab_x[k] * dk_139[k]
                       + dl_166[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, dk_140, dk_141, dk_142, \
                         dk_143, dk_144, dl_167, dl_168, dl_169, dl_170, \
                         dl_180 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * dk_140[k]
                       + dl_167[k];

            t_141[k] = -ab_x[k] * dk_141[k]
                       + dl_168[k];

            t_142[k] = -ab_x[k] * dk_142[k]
                       + dl_169[k];

            t_143[k] = -ab_x[k] * dk_143[k]
                       + dl_170[k];

            t_144[k] = -ab_x[k] * dk_144[k]
                       + dl_180[k];
        }
    }
}

static auto
compute_hrr_fk_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dk, const size_t dl, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);

        const auto *dk_108 = buffer.data(dk + 108 * ncomps + c);
        const auto *dk_109 = buffer.data(dk + 109 * ncomps + c);
        const auto *dk_110 = buffer.data(dk + 110 * ncomps + c);
        const auto *dk_111 = buffer.data(dk + 111 * ncomps + c);
        const auto *dk_112 = buffer.data(dk + 112 * ncomps + c);
        const auto *dk_113 = buffer.data(dk + 113 * ncomps + c);
        const auto *dk_114 = buffer.data(dk + 114 * ncomps + c);
        const auto *dk_115 = buffer.data(dk + 115 * ncomps + c);
        const auto *dk_116 = buffer.data(dk + 116 * ncomps + c);
        const auto *dk_117 = buffer.data(dk + 117 * ncomps + c);
        const auto *dk_118 = buffer.data(dk + 118 * ncomps + c);
        const auto *dk_119 = buffer.data(dk + 119 * ncomps + c);
        const auto *dk_120 = buffer.data(dk + 120 * ncomps + c);
        const auto *dk_121 = buffer.data(dk + 121 * ncomps + c);
        const auto *dk_122 = buffer.data(dk + 122 * ncomps + c);
        const auto *dk_123 = buffer.data(dk + 123 * ncomps + c);
        const auto *dk_124 = buffer.data(dk + 124 * ncomps + c);
        const auto *dk_125 = buffer.data(dk + 125 * ncomps + c);
        const auto *dk_126 = buffer.data(dk + 126 * ncomps + c);
        const auto *dk_127 = buffer.data(dk + 127 * ncomps + c);
        const auto *dk_128 = buffer.data(dk + 128 * ncomps + c);
        const auto *dk_129 = buffer.data(dk + 129 * ncomps + c);
        const auto *dk_130 = buffer.data(dk + 130 * ncomps + c);
        const auto *dk_131 = buffer.data(dk + 131 * ncomps + c);
        const auto *dk_132 = buffer.data(dk + 132 * ncomps + c);
        const auto *dk_133 = buffer.data(dk + 133 * ncomps + c);
        const auto *dk_134 = buffer.data(dk + 134 * ncomps + c);
        const auto *dk_135 = buffer.data(dk + 135 * ncomps + c);
        const auto *dk_136 = buffer.data(dk + 136 * ncomps + c);
        const auto *dk_137 = buffer.data(dk + 137 * ncomps + c);
        const auto *dk_138 = buffer.data(dk + 138 * ncomps + c);
        const auto *dk_139 = buffer.data(dk + 139 * ncomps + c);
        const auto *dk_140 = buffer.data(dk + 140 * ncomps + c);
        const auto *dk_141 = buffer.data(dk + 141 * ncomps + c);
        const auto *dk_142 = buffer.data(dk + 142 * ncomps + c);
        const auto *dk_143 = buffer.data(dk + 143 * ncomps + c);
        const auto *dk_144 = buffer.data(dk + 144 * ncomps + c);
        const auto *dk_145 = buffer.data(dk + 145 * ncomps + c);
        const auto *dk_146 = buffer.data(dk + 146 * ncomps + c);
        const auto *dk_147 = buffer.data(dk + 147 * ncomps + c);
        const auto *dk_148 = buffer.data(dk + 148 * ncomps + c);
        const auto *dk_149 = buffer.data(dk + 149 * ncomps + c);
        const auto *dk_150 = buffer.data(dk + 150 * ncomps + c);
        const auto *dk_151 = buffer.data(dk + 151 * ncomps + c);
        const auto *dk_152 = buffer.data(dk + 152 * ncomps + c);
        const auto *dk_153 = buffer.data(dk + 153 * ncomps + c);
        const auto *dk_154 = buffer.data(dk + 154 * ncomps + c);
        const auto *dk_155 = buffer.data(dk + 155 * ncomps + c);
        const auto *dk_156 = buffer.data(dk + 156 * ncomps + c);
        const auto *dk_157 = buffer.data(dk + 157 * ncomps + c);
        const auto *dk_158 = buffer.data(dk + 158 * ncomps + c);
        const auto *dk_159 = buffer.data(dk + 159 * ncomps + c);
        const auto *dk_160 = buffer.data(dk + 160 * ncomps + c);
        const auto *dk_161 = buffer.data(dk + 161 * ncomps + c);
        const auto *dk_162 = buffer.data(dk + 162 * ncomps + c);
        const auto *dk_163 = buffer.data(dk + 163 * ncomps + c);
        const auto *dk_164 = buffer.data(dk + 164 * ncomps + c);
        const auto *dk_165 = buffer.data(dk + 165 * ncomps + c);
        const auto *dk_166 = buffer.data(dk + 166 * ncomps + c);
        const auto *dk_167 = buffer.data(dk + 167 * ncomps + c);
        const auto *dk_168 = buffer.data(dk + 168 * ncomps + c);
        const auto *dk_169 = buffer.data(dk + 169 * ncomps + c);
        const auto *dk_170 = buffer.data(dk + 170 * ncomps + c);
        const auto *dk_171 = buffer.data(dk + 171 * ncomps + c);
        const auto *dk_172 = buffer.data(dk + 172 * ncomps + c);
        const auto *dk_173 = buffer.data(dk + 173 * ncomps + c);
        const auto *dk_174 = buffer.data(dk + 174 * ncomps + c);
        const auto *dk_175 = buffer.data(dk + 175 * ncomps + c);
        const auto *dk_176 = buffer.data(dk + 176 * ncomps + c);
        const auto *dk_177 = buffer.data(dk + 177 * ncomps + c);
        const auto *dk_178 = buffer.data(dk + 178 * ncomps + c);
        const auto *dk_179 = buffer.data(dk + 179 * ncomps + c);
        const auto *dk_180 = buffer.data(dk + 180 * ncomps + c);
        const auto *dk_181 = buffer.data(dk + 181 * ncomps + c);
        const auto *dk_182 = buffer.data(dk + 182 * ncomps + c);
        const auto *dk_183 = buffer.data(dk + 183 * ncomps + c);
        const auto *dk_184 = buffer.data(dk + 184 * ncomps + c);
        const auto *dk_185 = buffer.data(dk + 185 * ncomps + c);
        const auto *dk_186 = buffer.data(dk + 186 * ncomps + c);
        const auto *dk_187 = buffer.data(dk + 187 * ncomps + c);
        const auto *dk_188 = buffer.data(dk + 188 * ncomps + c);
        const auto *dk_189 = buffer.data(dk + 189 * ncomps + c);
        const auto *dk_190 = buffer.data(dk + 190 * ncomps + c);
        const auto *dk_191 = buffer.data(dk + 191 * ncomps + c);
        const auto *dk_192 = buffer.data(dk + 192 * ncomps + c);
        const auto *dk_193 = buffer.data(dk + 193 * ncomps + c);
        const auto *dk_194 = buffer.data(dk + 194 * ncomps + c);
        const auto *dk_195 = buffer.data(dk + 195 * ncomps + c);
        const auto *dk_196 = buffer.data(dk + 196 * ncomps + c);
        const auto *dk_197 = buffer.data(dk + 197 * ncomps + c);
        const auto *dk_198 = buffer.data(dk + 198 * ncomps + c);
        const auto *dk_199 = buffer.data(dk + 199 * ncomps + c);
        const auto *dk_200 = buffer.data(dk + 200 * ncomps + c);
        const auto *dk_201 = buffer.data(dk + 201 * ncomps + c);
        const auto *dk_202 = buffer.data(dk + 202 * ncomps + c);
        const auto *dk_203 = buffer.data(dk + 203 * ncomps + c);
        const auto *dk_204 = buffer.data(dk + 204 * ncomps + c);
        const auto *dk_205 = buffer.data(dk + 205 * ncomps + c);
        const auto *dk_206 = buffer.data(dk + 206 * ncomps + c);
        const auto *dk_207 = buffer.data(dk + 207 * ncomps + c);
        const auto *dk_208 = buffer.data(dk + 208 * ncomps + c);
        const auto *dk_209 = buffer.data(dk + 209 * ncomps + c);
        const auto *dk_210 = buffer.data(dk + 210 * ncomps + c);
        const auto *dk_211 = buffer.data(dk + 211 * ncomps + c);
        const auto *dk_212 = buffer.data(dk + 212 * ncomps + c);
        const auto *dk_213 = buffer.data(dk + 213 * ncomps + c);
        const auto *dk_214 = buffer.data(dk + 214 * ncomps + c);
        const auto *dk_215 = buffer.data(dk + 215 * ncomps + c);

        const auto *dl_136 = buffer.data(dl + 136 * ncomps + c);
        const auto *dl_138 = buffer.data(dl + 138 * ncomps + c);
        const auto *dl_139 = buffer.data(dl + 139 * ncomps + c);
        const auto *dl_141 = buffer.data(dl + 141 * ncomps + c);
        const auto *dl_142 = buffer.data(dl + 142 * ncomps + c);
        const auto *dl_143 = buffer.data(dl + 143 * ncomps + c);
        const auto *dl_145 = buffer.data(dl + 145 * ncomps + c);
        const auto *dl_146 = buffer.data(dl + 146 * ncomps + c);
        const auto *dl_147 = buffer.data(dl + 147 * ncomps + c);
        const auto *dl_148 = buffer.data(dl + 148 * ncomps + c);
        const auto *dl_150 = buffer.data(dl + 150 * ncomps + c);
        const auto *dl_151 = buffer.data(dl + 151 * ncomps + c);
        const auto *dl_152 = buffer.data(dl + 152 * ncomps + c);
        const auto *dl_153 = buffer.data(dl + 153 * ncomps + c);
        const auto *dl_154 = buffer.data(dl + 154 * ncomps + c);
        const auto *dl_156 = buffer.data(dl + 156 * ncomps + c);
        const auto *dl_157 = buffer.data(dl + 157 * ncomps + c);
        const auto *dl_158 = buffer.data(dl + 158 * ncomps + c);
        const auto *dl_159 = buffer.data(dl + 159 * ncomps + c);
        const auto *dl_160 = buffer.data(dl + 160 * ncomps + c);
        const auto *dl_161 = buffer.data(dl + 161 * ncomps + c);
        const auto *dl_163 = buffer.data(dl + 163 * ncomps + c);
        const auto *dl_164 = buffer.data(dl + 164 * ncomps + c);
        const auto *dl_165 = buffer.data(dl + 165 * ncomps + c);
        const auto *dl_166 = buffer.data(dl + 166 * ncomps + c);
        const auto *dl_167 = buffer.data(dl + 167 * ncomps + c);
        const auto *dl_168 = buffer.data(dl + 168 * ncomps + c);
        const auto *dl_169 = buffer.data(dl + 169 * ncomps + c);
        const auto *dl_171 = buffer.data(dl + 171 * ncomps + c);
        const auto *dl_172 = buffer.data(dl + 172 * ncomps + c);
        const auto *dl_173 = buffer.data(dl + 173 * ncomps + c);
        const auto *dl_174 = buffer.data(dl + 174 * ncomps + c);
        const auto *dl_175 = buffer.data(dl + 175 * ncomps + c);
        const auto *dl_176 = buffer.data(dl + 176 * ncomps + c);
        const auto *dl_177 = buffer.data(dl + 177 * ncomps + c);
        const auto *dl_178 = buffer.data(dl + 178 * ncomps + c);
        const auto *dl_181 = buffer.data(dl + 181 * ncomps + c);
        const auto *dl_182 = buffer.data(dl + 182 * ncomps + c);
        const auto *dl_183 = buffer.data(dl + 183 * ncomps + c);
        const auto *dl_184 = buffer.data(dl + 184 * ncomps + c);
        const auto *dl_185 = buffer.data(dl + 185 * ncomps + c);
        const auto *dl_186 = buffer.data(dl + 186 * ncomps + c);
        const auto *dl_187 = buffer.data(dl + 187 * ncomps + c);
        const auto *dl_188 = buffer.data(dl + 188 * ncomps + c);
        const auto *dl_189 = buffer.data(dl + 189 * ncomps + c);
        const auto *dl_190 = buffer.data(dl + 190 * ncomps + c);
        const auto *dl_191 = buffer.data(dl + 191 * ncomps + c);
        const auto *dl_192 = buffer.data(dl + 192 * ncomps + c);
        const auto *dl_193 = buffer.data(dl + 193 * ncomps + c);
        const auto *dl_194 = buffer.data(dl + 194 * ncomps + c);
        const auto *dl_195 = buffer.data(dl + 195 * ncomps + c);
        const auto *dl_196 = buffer.data(dl + 196 * ncomps + c);
        const auto *dl_197 = buffer.data(dl + 197 * ncomps + c);
        const auto *dl_198 = buffer.data(dl + 198 * ncomps + c);
        const auto *dl_199 = buffer.data(dl + 199 * ncomps + c);
        const auto *dl_200 = buffer.data(dl + 200 * ncomps + c);
        const auto *dl_201 = buffer.data(dl + 201 * ncomps + c);
        const auto *dl_202 = buffer.data(dl + 202 * ncomps + c);
        const auto *dl_203 = buffer.data(dl + 203 * ncomps + c);
        const auto *dl_204 = buffer.data(dl + 204 * ncomps + c);
        const auto *dl_205 = buffer.data(dl + 205 * ncomps + c);
        const auto *dl_206 = buffer.data(dl + 206 * ncomps + c);
        const auto *dl_207 = buffer.data(dl + 207 * ncomps + c);
        const auto *dl_208 = buffer.data(dl + 208 * ncomps + c);
        const auto *dl_209 = buffer.data(dl + 209 * ncomps + c);
        const auto *dl_210 = buffer.data(dl + 210 * ncomps + c);
        const auto *dl_211 = buffer.data(dl + 211 * ncomps + c);
        const auto *dl_212 = buffer.data(dl + 212 * ncomps + c);
        const auto *dl_213 = buffer.data(dl + 213 * ncomps + c);
        const auto *dl_214 = buffer.data(dl + 214 * ncomps + c);
        const auto *dl_215 = buffer.data(dl + 215 * ncomps + c);
        const auto *dl_216 = buffer.data(dl + 216 * ncomps + c);
        const auto *dl_217 = buffer.data(dl + 217 * ncomps + c);
        const auto *dl_218 = buffer.data(dl + 218 * ncomps + c);
        const auto *dl_219 = buffer.data(dl + 219 * ncomps + c);
        const auto *dl_220 = buffer.data(dl + 220 * ncomps + c);
        const auto *dl_221 = buffer.data(dl + 221 * ncomps + c);
        const auto *dl_222 = buffer.data(dl + 222 * ncomps + c);
        const auto *dl_223 = buffer.data(dl + 223 * ncomps + c);
        const auto *dl_225 = buffer.data(dl + 225 * ncomps + c);
        const auto *dl_226 = buffer.data(dl + 226 * ncomps + c);
        const auto *dl_227 = buffer.data(dl + 227 * ncomps + c);
        const auto *dl_228 = buffer.data(dl + 228 * ncomps + c);
        const auto *dl_229 = buffer.data(dl + 229 * ncomps + c);
        const auto *dl_230 = buffer.data(dl + 230 * ncomps + c);
        const auto *dl_231 = buffer.data(dl + 231 * ncomps + c);
        const auto *dl_232 = buffer.data(dl + 232 * ncomps + c);
        const auto *dl_233 = buffer.data(dl + 233 * ncomps + c);
        const auto *dl_234 = buffer.data(dl + 234 * ncomps + c);
        const auto *dl_235 = buffer.data(dl + 235 * ncomps + c);
        const auto *dl_236 = buffer.data(dl + 236 * ncomps + c);
        const auto *dl_237 = buffer.data(dl + 237 * ncomps + c);
        const auto *dl_238 = buffer.data(dl + 238 * ncomps + c);
        const auto *dl_239 = buffer.data(dl + 239 * ncomps + c);
        const auto *dl_240 = buffer.data(dl + 240 * ncomps + c);
        const auto *dl_241 = buffer.data(dl + 241 * ncomps + c);
        const auto *dl_242 = buffer.data(dl + 242 * ncomps + c);
        const auto *dl_243 = buffer.data(dl + 243 * ncomps + c);
        const auto *dl_244 = buffer.data(dl + 244 * ncomps + c);
        const auto *dl_245 = buffer.data(dl + 245 * ncomps + c);
        const auto *dl_246 = buffer.data(dl + 246 * ncomps + c);
        const auto *dl_247 = buffer.data(dl + 247 * ncomps + c);
        const auto *dl_248 = buffer.data(dl + 248 * ncomps + c);
        const auto *dl_249 = buffer.data(dl + 249 * ncomps + c);
        const auto *dl_250 = buffer.data(dl + 250 * ncomps + c);
        const auto *dl_251 = buffer.data(dl + 251 * ncomps + c);
        const auto *dl_252 = buffer.data(dl + 252 * ncomps + c);
        const auto *dl_253 = buffer.data(dl + 253 * ncomps + c);
        const auto *dl_254 = buffer.data(dl + 254 * ncomps + c);
        const auto *dl_255 = buffer.data(dl + 255 * ncomps + c);
        const auto *dl_256 = buffer.data(dl + 256 * ncomps + c);
        const auto *dl_257 = buffer.data(dl + 257 * ncomps + c);
        const auto *dl_258 = buffer.data(dl + 258 * ncomps + c);
        const auto *dl_259 = buffer.data(dl + 259 * ncomps + c);
        const auto *dl_260 = buffer.data(dl + 260 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, dk_145, dk_146, dk_147, \
                         dk_148, dk_149, dl_181, dl_182, dl_183, dl_184, \
                         dl_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * dk_145[k]
                       + dl_181[k];

            t_146[k] = -ab_x[k] * dk_146[k]
                       + dl_182[k];

            t_147[k] = -ab_x[k] * dk_147[k]
                       + dl_183[k];

            t_148[k] = -ab_x[k] * dk_148[k]
                       + dl_184[k];

            t_149[k] = -ab_x[k] * dk_149[k]
                       + dl_185[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, dk_150, dk_151, dk_152, \
                         dk_153, dk_154, dl_186, dl_187, dl_188, dl_189, \
                         dl_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * dk_150[k]
                       + dl_186[k];

            t_151[k] = -ab_x[k] * dk_151[k]
                       + dl_187[k];

            t_152[k] = -ab_x[k] * dk_152[k]
                       + dl_188[k];

            t_153[k] = -ab_x[k] * dk_153[k]
                       + dl_189[k];

            t_154[k] = -ab_x[k] * dk_154[k]
                       + dl_190[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, dk_155, dk_156, dk_157, \
                         dk_158, dk_159, dl_191, dl_192, dl_193, dl_194, \
                         dl_195 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * dk_155[k]
                       + dl_191[k];

            t_156[k] = -ab_x[k] * dk_156[k]
                       + dl_192[k];

            t_157[k] = -ab_x[k] * dk_157[k]
                       + dl_193[k];

            t_158[k] = -ab_x[k] * dk_158[k]
                       + dl_194[k];

            t_159[k] = -ab_x[k] * dk_159[k]
                       + dl_195[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, dk_160, dk_161, dk_162, \
                         dk_163, dk_164, dl_196, dl_197, dl_198, dl_199, \
                         dl_200 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * dk_160[k]
                       + dl_196[k];

            t_161[k] = -ab_x[k] * dk_161[k]
                       + dl_197[k];

            t_162[k] = -ab_x[k] * dk_162[k]
                       + dl_198[k];

            t_163[k] = -ab_x[k] * dk_163[k]
                       + dl_199[k];

            t_164[k] = -ab_x[k] * dk_164[k]
                       + dl_200[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, dk_165, dk_166, dk_167, \
                         dk_168, dk_169, dl_201, dl_202, dl_203, dl_204, \
                         dl_205 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * dk_165[k]
                       + dl_201[k];

            t_166[k] = -ab_x[k] * dk_166[k]
                       + dl_202[k];

            t_167[k] = -ab_x[k] * dk_167[k]
                       + dl_203[k];

            t_168[k] = -ab_x[k] * dk_168[k]
                       + dl_204[k];

            t_169[k] = -ab_x[k] * dk_169[k]
                       + dl_205[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, dk_170, dk_171, dk_172, \
                         dk_173, dk_174, dl_206, dl_207, dl_208, dl_209, \
                         dl_210 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * dk_170[k]
                       + dl_206[k];

            t_171[k] = -ab_x[k] * dk_171[k]
                       + dl_207[k];

            t_172[k] = -ab_x[k] * dk_172[k]
                       + dl_208[k];

            t_173[k] = -ab_x[k] * dk_173[k]
                       + dl_209[k];

            t_174[k] = -ab_x[k] * dk_174[k]
                       + dl_210[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, dk_175, dk_176, dk_177, \
                         dk_178, dk_179, dl_211, dl_212, dl_213, dl_214, \
                         dl_215 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * dk_175[k]
                       + dl_211[k];

            t_176[k] = -ab_x[k] * dk_176[k]
                       + dl_212[k];

            t_177[k] = -ab_x[k] * dk_177[k]
                       + dl_213[k];

            t_178[k] = -ab_x[k] * dk_178[k]
                       + dl_214[k];

            t_179[k] = -ab_x[k] * dk_179[k]
                       + dl_215[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, dk_180, dk_181, dk_182, \
                         dk_183, dk_184, dl_225, dl_226, dl_227, dl_228, \
                         dl_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * dk_180[k]
                       + dl_225[k];

            t_181[k] = -ab_x[k] * dk_181[k]
                       + dl_226[k];

            t_182[k] = -ab_x[k] * dk_182[k]
                       + dl_227[k];

            t_183[k] = -ab_x[k] * dk_183[k]
                       + dl_228[k];

            t_184[k] = -ab_x[k] * dk_184[k]
                       + dl_229[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, dk_185, dk_186, dk_187, \
                         dk_188, dk_189, dl_230, dl_231, dl_232, dl_233, \
                         dl_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * dk_185[k]
                       + dl_230[k];

            t_186[k] = -ab_x[k] * dk_186[k]
                       + dl_231[k];

            t_187[k] = -ab_x[k] * dk_187[k]
                       + dl_232[k];

            t_188[k] = -ab_x[k] * dk_188[k]
                       + dl_233[k];

            t_189[k] = -ab_x[k] * dk_189[k]
                       + dl_234[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, dk_190, dk_191, dk_192, \
                         dk_193, dk_194, dl_235, dl_236, dl_237, dl_238, \
                         dl_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * dk_190[k]
                       + dl_235[k];

            t_191[k] = -ab_x[k] * dk_191[k]
                       + dl_236[k];

            t_192[k] = -ab_x[k] * dk_192[k]
                       + dl_237[k];

            t_193[k] = -ab_x[k] * dk_193[k]
                       + dl_238[k];

            t_194[k] = -ab_x[k] * dk_194[k]
                       + dl_239[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, dk_195, dk_196, dk_197, \
                         dk_198, dk_199, dl_240, dl_241, dl_242, dl_243, \
                         dl_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * dk_195[k]
                       + dl_240[k];

            t_196[k] = -ab_x[k] * dk_196[k]
                       + dl_241[k];

            t_197[k] = -ab_x[k] * dk_197[k]
                       + dl_242[k];

            t_198[k] = -ab_x[k] * dk_198[k]
                       + dl_243[k];

            t_199[k] = -ab_x[k] * dk_199[k]
                       + dl_244[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, dk_200, dk_201, dk_202, \
                         dk_203, dk_204, dl_245, dl_246, dl_247, dl_248, \
                         dl_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * dk_200[k]
                       + dl_245[k];

            t_201[k] = -ab_x[k] * dk_201[k]
                       + dl_246[k];

            t_202[k] = -ab_x[k] * dk_202[k]
                       + dl_247[k];

            t_203[k] = -ab_x[k] * dk_203[k]
                       + dl_248[k];

            t_204[k] = -ab_x[k] * dk_204[k]
                       + dl_249[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, dk_205, dk_206, dk_207, \
                         dk_208, dk_209, dl_250, dl_251, dl_252, dl_253, \
                         dl_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * dk_205[k]
                       + dl_250[k];

            t_206[k] = -ab_x[k] * dk_206[k]
                       + dl_251[k];

            t_207[k] = -ab_x[k] * dk_207[k]
                       + dl_252[k];

            t_208[k] = -ab_x[k] * dk_208[k]
                       + dl_253[k];

            t_209[k] = -ab_x[k] * dk_209[k]
                       + dl_254[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, dk_210, dk_211, dk_212, \
                         dk_213, dk_214, dl_255, dl_256, dl_257, dl_258, \
                         dl_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * dk_210[k]
                       + dl_255[k];

            t_211[k] = -ab_x[k] * dk_211[k]
                       + dl_256[k];

            t_212[k] = -ab_x[k] * dk_212[k]
                       + dl_257[k];

            t_213[k] = -ab_x[k] * dk_213[k]
                       + dl_258[k];

            t_214[k] = -ab_x[k] * dk_214[k]
                       + dl_259[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, ab_x, ab_y, dk_108, dk_109, dk_110, \
                         dk_215, dl_136, dl_138, dl_139, dl_260 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * dk_215[k]
                       + dl_260[k];

            t_216[k] = -ab_y[k] * dk_108[k]
                       + dl_136[k];

            t_217[k] = -ab_y[k] * dk_109[k]
                       + dl_138[k];

            t_218[k] = -ab_y[k] * dk_110[k]
                       + dl_139[k];
        }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, ab_y, dk_111, dk_112, dk_113, \
                         dk_114, dk_115, dl_141, dl_142, dl_143, dl_145, \
                         dl_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_219[k] = -ab_y[k] * dk_111[k]
                       + dl_141[k];

            t_220[k] = -ab_y[k] * dk_112[k]
                       + dl_142[k];

            t_221[k] = -ab_y[k] * dk_113[k]
                       + dl_143[k];

            t_222[k] = -ab_y[k] * dk_114[k]
                       + dl_145[k];

            t_223[k] = -ab_y[k] * dk_115[k]
                       + dl_146[k];
        }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, ab_y, dk_116, dk_117, dk_118, \
                         dk_119, dk_120, dl_147, dl_148, dl_150, dl_151, \
                         dl_152 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_224[k] = -ab_y[k] * dk_116[k]
                       + dl_147[k];

            t_225[k] = -ab_y[k] * dk_117[k]
                       + dl_148[k];

            t_226[k] = -ab_y[k] * dk_118[k]
                       + dl_150[k];

            t_227[k] = -ab_y[k] * dk_119[k]
                       + dl_151[k];

            t_228[k] = -ab_y[k] * dk_120[k]
                       + dl_152[k];
        }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, ab_y, dk_121, dk_122, dk_123, \
                         dk_124, dk_125, dl_153, dl_154, dl_156, dl_157, \
                         dl_158 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_229[k] = -ab_y[k] * dk_121[k]
                       + dl_153[k];

            t_230[k] = -ab_y[k] * dk_122[k]
                       + dl_154[k];

            t_231[k] = -ab_y[k] * dk_123[k]
                       + dl_156[k];

            t_232[k] = -ab_y[k] * dk_124[k]
                       + dl_157[k];

            t_233[k] = -ab_y[k] * dk_125[k]
                       + dl_158[k];
        }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_y, dk_126, dk_127, dk_128, \
                         dk_129, dk_130, dl_159, dl_160, dl_161, dl_163, \
                         dl_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_234[k] = -ab_y[k] * dk_126[k]
                       + dl_159[k];

            t_235[k] = -ab_y[k] * dk_127[k]
                       + dl_160[k];

            t_236[k] = -ab_y[k] * dk_128[k]
                       + dl_161[k];

            t_237[k] = -ab_y[k] * dk_129[k]
                       + dl_163[k];

            t_238[k] = -ab_y[k] * dk_130[k]
                       + dl_164[k];
        }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_y, dk_131, dk_132, dk_133, \
                         dk_134, dk_135, dl_165, dl_166, dl_167, dl_168, \
                         dl_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_239[k] = -ab_y[k] * dk_131[k]
                       + dl_165[k];

            t_240[k] = -ab_y[k] * dk_132[k]
                       + dl_166[k];

            t_241[k] = -ab_y[k] * dk_133[k]
                       + dl_167[k];

            t_242[k] = -ab_y[k] * dk_134[k]
                       + dl_168[k];

            t_243[k] = -ab_y[k] * dk_135[k]
                       + dl_169[k];
        }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, ab_y, dk_136, dk_137, dk_138, \
                         dk_139, dk_140, dl_171, dl_172, dl_173, dl_174, \
                         dl_175 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_244[k] = -ab_y[k] * dk_136[k]
                       + dl_171[k];

            t_245[k] = -ab_y[k] * dk_137[k]
                       + dl_172[k];

            t_246[k] = -ab_y[k] * dk_138[k]
                       + dl_173[k];

            t_247[k] = -ab_y[k] * dk_139[k]
                       + dl_174[k];

            t_248[k] = -ab_y[k] * dk_140[k]
                       + dl_175[k];
        }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, ab_y, dk_141, dk_142, dk_143, \
                         dk_144, dk_145, dl_176, dl_177, dl_178, dl_181, \
                         dl_183 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_249[k] = -ab_y[k] * dk_141[k]
                       + dl_176[k];

            t_250[k] = -ab_y[k] * dk_142[k]
                       + dl_177[k];

            t_251[k] = -ab_y[k] * dk_143[k]
                       + dl_178[k];

            t_252[k] = -ab_y[k] * dk_144[k]
                       + dl_181[k];

            t_253[k] = -ab_y[k] * dk_145[k]
                       + dl_183[k];
        }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, ab_y, dk_146, dk_147, dk_148, \
                         dk_149, dk_150, dl_184, dl_186, dl_187, dl_188, \
                         dl_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_254[k] = -ab_y[k] * dk_146[k]
                       + dl_184[k];

            t_255[k] = -ab_y[k] * dk_147[k]
                       + dl_186[k];

            t_256[k] = -ab_y[k] * dk_148[k]
                       + dl_187[k];

            t_257[k] = -ab_y[k] * dk_149[k]
                       + dl_188[k];

            t_258[k] = -ab_y[k] * dk_150[k]
                       + dl_190[k];
        }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, ab_y, dk_151, dk_152, dk_153, \
                         dk_154, dk_155, dl_191, dl_192, dl_193, dl_195, \
                         dl_196 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_259[k] = -ab_y[k] * dk_151[k]
                       + dl_191[k];

            t_260[k] = -ab_y[k] * dk_152[k]
                       + dl_192[k];

            t_261[k] = -ab_y[k] * dk_153[k]
                       + dl_193[k];

            t_262[k] = -ab_y[k] * dk_154[k]
                       + dl_195[k];

            t_263[k] = -ab_y[k] * dk_155[k]
                       + dl_196[k];
        }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, ab_y, dk_156, dk_157, dk_158, \
                         dk_159, dk_160, dl_197, dl_198, dl_199, dl_201, \
                         dl_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_264[k] = -ab_y[k] * dk_156[k]
                       + dl_197[k];

            t_265[k] = -ab_y[k] * dk_157[k]
                       + dl_198[k];

            t_266[k] = -ab_y[k] * dk_158[k]
                       + dl_199[k];

            t_267[k] = -ab_y[k] * dk_159[k]
                       + dl_201[k];

            t_268[k] = -ab_y[k] * dk_160[k]
                       + dl_202[k];
        }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, ab_y, dk_161, dk_162, dk_163, \
                         dk_164, dk_165, dl_203, dl_204, dl_205, dl_206, \
                         dl_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_269[k] = -ab_y[k] * dk_161[k]
                       + dl_203[k];

            t_270[k] = -ab_y[k] * dk_162[k]
                       + dl_204[k];

            t_271[k] = -ab_y[k] * dk_163[k]
                       + dl_205[k];

            t_272[k] = -ab_y[k] * dk_164[k]
                       + dl_206[k];

            t_273[k] = -ab_y[k] * dk_165[k]
                       + dl_208[k];
        }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, ab_y, dk_166, dk_167, dk_168, \
                         dk_169, dk_170, dl_209, dl_210, dl_211, dl_212, \
                         dl_213 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_274[k] = -ab_y[k] * dk_166[k]
                       + dl_209[k];

            t_275[k] = -ab_y[k] * dk_167[k]
                       + dl_210[k];

            t_276[k] = -ab_y[k] * dk_168[k]
                       + dl_211[k];

            t_277[k] = -ab_y[k] * dk_169[k]
                       + dl_212[k];

            t_278[k] = -ab_y[k] * dk_170[k]
                       + dl_213[k];
        }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, t_283, ab_y, dk_171, dk_172, dk_173, \
                         dk_174, dk_175, dl_214, dl_216, dl_217, dl_218, \
                         dl_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_279[k] = -ab_y[k] * dk_171[k]
                       + dl_214[k];

            t_280[k] = -ab_y[k] * dk_172[k]
                       + dl_216[k];

            t_281[k] = -ab_y[k] * dk_173[k]
                       + dl_217[k];

            t_282[k] = -ab_y[k] * dk_174[k]
                       + dl_218[k];

            t_283[k] = -ab_y[k] * dk_175[k]
                       + dl_219[k];
        }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, ab_y, dk_176, dk_177, dk_178, \
                         dk_179, dk_180, dl_220, dl_221, dl_222, dl_223, \
                         dl_226 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_284[k] = -ab_y[k] * dk_176[k]
                       + dl_220[k];

            t_285[k] = -ab_y[k] * dk_177[k]
                       + dl_221[k];

            t_286[k] = -ab_y[k] * dk_178[k]
                       + dl_222[k];

            t_287[k] = -ab_y[k] * dk_179[k]
                       + dl_223[k];

            t_288[k] = -ab_y[k] * dk_180[k]
                       + dl_226[k];
        }
    }
}

static auto
compute_hrr_fk_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dk, const size_t dl, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_289 = buffer.data(target + 289 * ncomps + c);
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dk_180 = buffer.data(dk + 180 * ncomps + c);
        const auto *dk_181 = buffer.data(dk + 181 * ncomps + c);
        const auto *dk_182 = buffer.data(dk + 182 * ncomps + c);
        const auto *dk_183 = buffer.data(dk + 183 * ncomps + c);
        const auto *dk_184 = buffer.data(dk + 184 * ncomps + c);
        const auto *dk_185 = buffer.data(dk + 185 * ncomps + c);
        const auto *dk_186 = buffer.data(dk + 186 * ncomps + c);
        const auto *dk_187 = buffer.data(dk + 187 * ncomps + c);
        const auto *dk_188 = buffer.data(dk + 188 * ncomps + c);
        const auto *dk_189 = buffer.data(dk + 189 * ncomps + c);
        const auto *dk_190 = buffer.data(dk + 190 * ncomps + c);
        const auto *dk_191 = buffer.data(dk + 191 * ncomps + c);
        const auto *dk_192 = buffer.data(dk + 192 * ncomps + c);
        const auto *dk_193 = buffer.data(dk + 193 * ncomps + c);
        const auto *dk_194 = buffer.data(dk + 194 * ncomps + c);
        const auto *dk_195 = buffer.data(dk + 195 * ncomps + c);
        const auto *dk_196 = buffer.data(dk + 196 * ncomps + c);
        const auto *dk_197 = buffer.data(dk + 197 * ncomps + c);
        const auto *dk_198 = buffer.data(dk + 198 * ncomps + c);
        const auto *dk_199 = buffer.data(dk + 199 * ncomps + c);
        const auto *dk_200 = buffer.data(dk + 200 * ncomps + c);
        const auto *dk_201 = buffer.data(dk + 201 * ncomps + c);
        const auto *dk_202 = buffer.data(dk + 202 * ncomps + c);
        const auto *dk_203 = buffer.data(dk + 203 * ncomps + c);
        const auto *dk_204 = buffer.data(dk + 204 * ncomps + c);
        const auto *dk_205 = buffer.data(dk + 205 * ncomps + c);
        const auto *dk_206 = buffer.data(dk + 206 * ncomps + c);
        const auto *dk_207 = buffer.data(dk + 207 * ncomps + c);
        const auto *dk_208 = buffer.data(dk + 208 * ncomps + c);
        const auto *dk_209 = buffer.data(dk + 209 * ncomps + c);
        const auto *dk_210 = buffer.data(dk + 210 * ncomps + c);
        const auto *dk_211 = buffer.data(dk + 211 * ncomps + c);
        const auto *dk_212 = buffer.data(dk + 212 * ncomps + c);
        const auto *dk_213 = buffer.data(dk + 213 * ncomps + c);
        const auto *dk_214 = buffer.data(dk + 214 * ncomps + c);
        const auto *dk_215 = buffer.data(dk + 215 * ncomps + c);

        const auto *dl_227 = buffer.data(dl + 227 * ncomps + c);
        const auto *dl_228 = buffer.data(dl + 228 * ncomps + c);
        const auto *dl_229 = buffer.data(dl + 229 * ncomps + c);
        const auto *dl_230 = buffer.data(dl + 230 * ncomps + c);
        const auto *dl_231 = buffer.data(dl + 231 * ncomps + c);
        const auto *dl_232 = buffer.data(dl + 232 * ncomps + c);
        const auto *dl_233 = buffer.data(dl + 233 * ncomps + c);
        const auto *dl_234 = buffer.data(dl + 234 * ncomps + c);
        const auto *dl_235 = buffer.data(dl + 235 * ncomps + c);
        const auto *dl_236 = buffer.data(dl + 236 * ncomps + c);
        const auto *dl_237 = buffer.data(dl + 237 * ncomps + c);
        const auto *dl_238 = buffer.data(dl + 238 * ncomps + c);
        const auto *dl_239 = buffer.data(dl + 239 * ncomps + c);
        const auto *dl_240 = buffer.data(dl + 240 * ncomps + c);
        const auto *dl_241 = buffer.data(dl + 241 * ncomps + c);
        const auto *dl_242 = buffer.data(dl + 242 * ncomps + c);
        const auto *dl_243 = buffer.data(dl + 243 * ncomps + c);
        const auto *dl_244 = buffer.data(dl + 244 * ncomps + c);
        const auto *dl_245 = buffer.data(dl + 245 * ncomps + c);
        const auto *dl_246 = buffer.data(dl + 246 * ncomps + c);
        const auto *dl_247 = buffer.data(dl + 247 * ncomps + c);
        const auto *dl_248 = buffer.data(dl + 248 * ncomps + c);
        const auto *dl_249 = buffer.data(dl + 249 * ncomps + c);
        const auto *dl_250 = buffer.data(dl + 250 * ncomps + c);
        const auto *dl_251 = buffer.data(dl + 251 * ncomps + c);
        const auto *dl_252 = buffer.data(dl + 252 * ncomps + c);
        const auto *dl_253 = buffer.data(dl + 253 * ncomps + c);
        const auto *dl_254 = buffer.data(dl + 254 * ncomps + c);
        const auto *dl_255 = buffer.data(dl + 255 * ncomps + c);
        const auto *dl_256 = buffer.data(dl + 256 * ncomps + c);
        const auto *dl_257 = buffer.data(dl + 257 * ncomps + c);
        const auto *dl_258 = buffer.data(dl + 258 * ncomps + c);
        const auto *dl_259 = buffer.data(dl + 259 * ncomps + c);
        const auto *dl_260 = buffer.data(dl + 260 * ncomps + c);
        const auto *dl_261 = buffer.data(dl + 261 * ncomps + c);
        const auto *dl_262 = buffer.data(dl + 262 * ncomps + c);
        const auto *dl_263 = buffer.data(dl + 263 * ncomps + c);
        const auto *dl_264 = buffer.data(dl + 264 * ncomps + c);
        const auto *dl_265 = buffer.data(dl + 265 * ncomps + c);
        const auto *dl_266 = buffer.data(dl + 266 * ncomps + c);
        const auto *dl_267 = buffer.data(dl + 267 * ncomps + c);
        const auto *dl_268 = buffer.data(dl + 268 * ncomps + c);
        const auto *dl_269 = buffer.data(dl + 269 * ncomps + c);

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, ab_y, dk_181, dk_182, dk_183, \
                         dk_184, dk_185, dl_228, dl_229, dl_231, dl_232, \
                         dl_233 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_289[k] = -ab_y[k] * dk_181[k]
                       + dl_228[k];

            t_290[k] = -ab_y[k] * dk_182[k]
                       + dl_229[k];

            t_291[k] = -ab_y[k] * dk_183[k]
                       + dl_231[k];

            t_292[k] = -ab_y[k] * dk_184[k]
                       + dl_232[k];

            t_293[k] = -ab_y[k] * dk_185[k]
                       + dl_233[k];
        }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_y, dk_186, dk_187, dk_188, \
                         dk_189, dk_190, dl_235, dl_236, dl_237, dl_238, \
                         dl_240 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_294[k] = -ab_y[k] * dk_186[k]
                       + dl_235[k];

            t_295[k] = -ab_y[k] * dk_187[k]
                       + dl_236[k];

            t_296[k] = -ab_y[k] * dk_188[k]
                       + dl_237[k];

            t_297[k] = -ab_y[k] * dk_189[k]
                       + dl_238[k];

            t_298[k] = -ab_y[k] * dk_190[k]
                       + dl_240[k];
        }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_y, dk_191, dk_192, dk_193, \
                         dk_194, dk_195, dl_241, dl_242, dl_243, dl_244, \
                         dl_246 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_299[k] = -ab_y[k] * dk_191[k]
                       + dl_241[k];

            t_300[k] = -ab_y[k] * dk_192[k]
                       + dl_242[k];

            t_301[k] = -ab_y[k] * dk_193[k]
                       + dl_243[k];

            t_302[k] = -ab_y[k] * dk_194[k]
                       + dl_244[k];

            t_303[k] = -ab_y[k] * dk_195[k]
                       + dl_246[k];
        }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_y, dk_196, dk_197, dk_198, \
                         dk_199, dk_200, dl_247, dl_248, dl_249, dl_250, \
                         dl_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_304[k] = -ab_y[k] * dk_196[k]
                       + dl_247[k];

            t_305[k] = -ab_y[k] * dk_197[k]
                       + dl_248[k];

            t_306[k] = -ab_y[k] * dk_198[k]
                       + dl_249[k];

            t_307[k] = -ab_y[k] * dk_199[k]
                       + dl_250[k];

            t_308[k] = -ab_y[k] * dk_200[k]
                       + dl_251[k];
        }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_y, dk_201, dk_202, dk_203, \
                         dk_204, dk_205, dl_253, dl_254, dl_255, dl_256, \
                         dl_257 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_309[k] = -ab_y[k] * dk_201[k]
                       + dl_253[k];

            t_310[k] = -ab_y[k] * dk_202[k]
                       + dl_254[k];

            t_311[k] = -ab_y[k] * dk_203[k]
                       + dl_255[k];

            t_312[k] = -ab_y[k] * dk_204[k]
                       + dl_256[k];

            t_313[k] = -ab_y[k] * dk_205[k]
                       + dl_257[k];
        }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, ab_y, dk_206, dk_207, dk_208, \
                         dk_209, dk_210, dl_258, dl_259, dl_261, dl_262, \
                         dl_263 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_314[k] = -ab_y[k] * dk_206[k]
                       + dl_258[k];

            t_315[k] = -ab_y[k] * dk_207[k]
                       + dl_259[k];

            t_316[k] = -ab_y[k] * dk_208[k]
                       + dl_261[k];

            t_317[k] = -ab_y[k] * dk_209[k]
                       + dl_262[k];

            t_318[k] = -ab_y[k] * dk_210[k]
                       + dl_263[k];
        }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, ab_y, dk_211, dk_212, dk_213, \
                         dk_214, dk_215, dl_264, dl_265, dl_266, dl_267, \
                         dl_268 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_319[k] = -ab_y[k] * dk_211[k]
                       + dl_264[k];

            t_320[k] = -ab_y[k] * dk_212[k]
                       + dl_265[k];

            t_321[k] = -ab_y[k] * dk_213[k]
                       + dl_266[k];

            t_322[k] = -ab_y[k] * dk_214[k]
                       + dl_267[k];

            t_323[k] = -ab_y[k] * dk_215[k]
                       + dl_268[k];
        }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, ab_z, dk_180, dk_181, dk_182, \
                         dk_183, dk_184, dl_227, dl_229, dl_230, dl_232, \
                         dl_233 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_324[k] = -ab_z[k] * dk_180[k]
                       + dl_227[k];

            t_325[k] = -ab_z[k] * dk_181[k]
                       + dl_229[k];

            t_326[k] = -ab_z[k] * dk_182[k]
                       + dl_230[k];

            t_327[k] = -ab_z[k] * dk_183[k]
                       + dl_232[k];

            t_328[k] = -ab_z[k] * dk_184[k]
                       + dl_233[k];
        }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, ab_z, dk_185, dk_186, dk_187, \
                         dk_188, dk_189, dl_234, dl_236, dl_237, dl_238, \
                         dl_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_329[k] = -ab_z[k] * dk_185[k]
                       + dl_234[k];

            t_330[k] = -ab_z[k] * dk_186[k]
                       + dl_236[k];

            t_331[k] = -ab_z[k] * dk_187[k]
                       + dl_237[k];

            t_332[k] = -ab_z[k] * dk_188[k]
                       + dl_238[k];

            t_333[k] = -ab_z[k] * dk_189[k]
                       + dl_239[k];
        }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, ab_z, dk_190, dk_191, dk_192, \
                         dk_193, dk_194, dl_241, dl_242, dl_243, dl_244, \
                         dl_245 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_334[k] = -ab_z[k] * dk_190[k]
                       + dl_241[k];

            t_335[k] = -ab_z[k] * dk_191[k]
                       + dl_242[k];

            t_336[k] = -ab_z[k] * dk_192[k]
                       + dl_243[k];

            t_337[k] = -ab_z[k] * dk_193[k]
                       + dl_244[k];

            t_338[k] = -ab_z[k] * dk_194[k]
                       + dl_245[k];
        }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, ab_z, dk_195, dk_196, dk_197, \
                         dk_198, dk_199, dl_247, dl_248, dl_249, dl_250, \
                         dl_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_339[k] = -ab_z[k] * dk_195[k]
                       + dl_247[k];

            t_340[k] = -ab_z[k] * dk_196[k]
                       + dl_248[k];

            t_341[k] = -ab_z[k] * dk_197[k]
                       + dl_249[k];

            t_342[k] = -ab_z[k] * dk_198[k]
                       + dl_250[k];

            t_343[k] = -ab_z[k] * dk_199[k]
                       + dl_251[k];
        }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, ab_z, dk_200, dk_201, dk_202, \
                         dk_203, dk_204, dl_252, dl_254, dl_255, dl_256, \
                         dl_257 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_344[k] = -ab_z[k] * dk_200[k]
                       + dl_252[k];

            t_345[k] = -ab_z[k] * dk_201[k]
                       + dl_254[k];

            t_346[k] = -ab_z[k] * dk_202[k]
                       + dl_255[k];

            t_347[k] = -ab_z[k] * dk_203[k]
                       + dl_256[k];

            t_348[k] = -ab_z[k] * dk_204[k]
                       + dl_257[k];
        }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, ab_z, dk_205, dk_206, dk_207, \
                         dk_208, dk_209, dl_258, dl_259, dl_260, dl_262, \
                         dl_263 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_349[k] = -ab_z[k] * dk_205[k]
                       + dl_258[k];

            t_350[k] = -ab_z[k] * dk_206[k]
                       + dl_259[k];

            t_351[k] = -ab_z[k] * dk_207[k]
                       + dl_260[k];

            t_352[k] = -ab_z[k] * dk_208[k]
                       + dl_262[k];

            t_353[k] = -ab_z[k] * dk_209[k]
                       + dl_263[k];
        }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, ab_z, dk_210, dk_211, dk_212, \
                         dk_213, dk_214, dl_264, dl_265, dl_266, dl_267, \
                         dl_268 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_354[k] = -ab_z[k] * dk_210[k]
                       + dl_264[k];

            t_355[k] = -ab_z[k] * dk_211[k]
                       + dl_265[k];

            t_356[k] = -ab_z[k] * dk_212[k]
                       + dl_266[k];

            t_357[k] = -ab_z[k] * dk_213[k]
                       + dl_267[k];

            t_358[k] = -ab_z[k] * dk_214[k]
                       + dl_268[k];
        }

#pragma omp simd aligned(t_359, ab_z, dk_215, dl_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_359[k] = -ab_z[k] * dk_215[k]
                       + dl_269[k];
        }
    }
}

auto
compute_hrr_fk(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t dk, const size_t dl, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_fk_piece0(buffer, coordinates, target, dk, dl, ncomps, nmax);

    compute_hrr_fk_piece1(buffer, coordinates, target, dk, dl, ncomps, nmax);

    compute_hrr_fk_piece2(buffer, coordinates, target, dk, dl, ncomps, nmax);
}

}  // namespace simdtrf
