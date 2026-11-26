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


#include "SimdTransferGL.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_gl_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fl, const size_t fm, const size_t ncomps,
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
        const auto *fl_36 = buffer.data(fl + 36 * ncomps + c);
        const auto *fl_37 = buffer.data(fl + 37 * ncomps + c);
        const auto *fl_38 = buffer.data(fl + 38 * ncomps + c);
        const auto *fl_39 = buffer.data(fl + 39 * ncomps + c);
        const auto *fl_40 = buffer.data(fl + 40 * ncomps + c);
        const auto *fl_41 = buffer.data(fl + 41 * ncomps + c);
        const auto *fl_42 = buffer.data(fl + 42 * ncomps + c);
        const auto *fl_43 = buffer.data(fl + 43 * ncomps + c);
        const auto *fl_44 = buffer.data(fl + 44 * ncomps + c);
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
        const auto *fl_81 = buffer.data(fl + 81 * ncomps + c);
        const auto *fl_82 = buffer.data(fl + 82 * ncomps + c);
        const auto *fl_83 = buffer.data(fl + 83 * ncomps + c);
        const auto *fl_84 = buffer.data(fl + 84 * ncomps + c);
        const auto *fl_85 = buffer.data(fl + 85 * ncomps + c);
        const auto *fl_86 = buffer.data(fl + 86 * ncomps + c);
        const auto *fl_87 = buffer.data(fl + 87 * ncomps + c);
        const auto *fl_88 = buffer.data(fl + 88 * ncomps + c);
        const auto *fl_89 = buffer.data(fl + 89 * ncomps + c);
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
        const auto *fl_126 = buffer.data(fl + 126 * ncomps + c);
        const auto *fl_127 = buffer.data(fl + 127 * ncomps + c);
        const auto *fl_128 = buffer.data(fl + 128 * ncomps + c);
        const auto *fl_129 = buffer.data(fl + 129 * ncomps + c);
        const auto *fl_130 = buffer.data(fl + 130 * ncomps + c);
        const auto *fl_131 = buffer.data(fl + 131 * ncomps + c);
        const auto *fl_132 = buffer.data(fl + 132 * ncomps + c);
        const auto *fl_133 = buffer.data(fl + 133 * ncomps + c);
        const auto *fl_134 = buffer.data(fl + 134 * ncomps + c);
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

        const auto *fm_0 = buffer.data(fm + 0 * ncomps + c);
        const auto *fm_1 = buffer.data(fm + 1 * ncomps + c);
        const auto *fm_2 = buffer.data(fm + 2 * ncomps + c);
        const auto *fm_3 = buffer.data(fm + 3 * ncomps + c);
        const auto *fm_4 = buffer.data(fm + 4 * ncomps + c);
        const auto *fm_5 = buffer.data(fm + 5 * ncomps + c);
        const auto *fm_6 = buffer.data(fm + 6 * ncomps + c);
        const auto *fm_7 = buffer.data(fm + 7 * ncomps + c);
        const auto *fm_8 = buffer.data(fm + 8 * ncomps + c);
        const auto *fm_9 = buffer.data(fm + 9 * ncomps + c);
        const auto *fm_10 = buffer.data(fm + 10 * ncomps + c);
        const auto *fm_11 = buffer.data(fm + 11 * ncomps + c);
        const auto *fm_12 = buffer.data(fm + 12 * ncomps + c);
        const auto *fm_13 = buffer.data(fm + 13 * ncomps + c);
        const auto *fm_14 = buffer.data(fm + 14 * ncomps + c);
        const auto *fm_15 = buffer.data(fm + 15 * ncomps + c);
        const auto *fm_16 = buffer.data(fm + 16 * ncomps + c);
        const auto *fm_17 = buffer.data(fm + 17 * ncomps + c);
        const auto *fm_18 = buffer.data(fm + 18 * ncomps + c);
        const auto *fm_19 = buffer.data(fm + 19 * ncomps + c);
        const auto *fm_20 = buffer.data(fm + 20 * ncomps + c);
        const auto *fm_21 = buffer.data(fm + 21 * ncomps + c);
        const auto *fm_22 = buffer.data(fm + 22 * ncomps + c);
        const auto *fm_23 = buffer.data(fm + 23 * ncomps + c);
        const auto *fm_24 = buffer.data(fm + 24 * ncomps + c);
        const auto *fm_25 = buffer.data(fm + 25 * ncomps + c);
        const auto *fm_26 = buffer.data(fm + 26 * ncomps + c);
        const auto *fm_27 = buffer.data(fm + 27 * ncomps + c);
        const auto *fm_28 = buffer.data(fm + 28 * ncomps + c);
        const auto *fm_29 = buffer.data(fm + 29 * ncomps + c);
        const auto *fm_30 = buffer.data(fm + 30 * ncomps + c);
        const auto *fm_31 = buffer.data(fm + 31 * ncomps + c);
        const auto *fm_32 = buffer.data(fm + 32 * ncomps + c);
        const auto *fm_33 = buffer.data(fm + 33 * ncomps + c);
        const auto *fm_34 = buffer.data(fm + 34 * ncomps + c);
        const auto *fm_35 = buffer.data(fm + 35 * ncomps + c);
        const auto *fm_36 = buffer.data(fm + 36 * ncomps + c);
        const auto *fm_37 = buffer.data(fm + 37 * ncomps + c);
        const auto *fm_38 = buffer.data(fm + 38 * ncomps + c);
        const auto *fm_39 = buffer.data(fm + 39 * ncomps + c);
        const auto *fm_40 = buffer.data(fm + 40 * ncomps + c);
        const auto *fm_41 = buffer.data(fm + 41 * ncomps + c);
        const auto *fm_42 = buffer.data(fm + 42 * ncomps + c);
        const auto *fm_43 = buffer.data(fm + 43 * ncomps + c);
        const auto *fm_44 = buffer.data(fm + 44 * ncomps + c);
        const auto *fm_55 = buffer.data(fm + 55 * ncomps + c);
        const auto *fm_56 = buffer.data(fm + 56 * ncomps + c);
        const auto *fm_57 = buffer.data(fm + 57 * ncomps + c);
        const auto *fm_58 = buffer.data(fm + 58 * ncomps + c);
        const auto *fm_59 = buffer.data(fm + 59 * ncomps + c);
        const auto *fm_60 = buffer.data(fm + 60 * ncomps + c);
        const auto *fm_61 = buffer.data(fm + 61 * ncomps + c);
        const auto *fm_62 = buffer.data(fm + 62 * ncomps + c);
        const auto *fm_63 = buffer.data(fm + 63 * ncomps + c);
        const auto *fm_64 = buffer.data(fm + 64 * ncomps + c);
        const auto *fm_65 = buffer.data(fm + 65 * ncomps + c);
        const auto *fm_66 = buffer.data(fm + 66 * ncomps + c);
        const auto *fm_67 = buffer.data(fm + 67 * ncomps + c);
        const auto *fm_68 = buffer.data(fm + 68 * ncomps + c);
        const auto *fm_69 = buffer.data(fm + 69 * ncomps + c);
        const auto *fm_70 = buffer.data(fm + 70 * ncomps + c);
        const auto *fm_71 = buffer.data(fm + 71 * ncomps + c);
        const auto *fm_72 = buffer.data(fm + 72 * ncomps + c);
        const auto *fm_73 = buffer.data(fm + 73 * ncomps + c);
        const auto *fm_74 = buffer.data(fm + 74 * ncomps + c);
        const auto *fm_75 = buffer.data(fm + 75 * ncomps + c);
        const auto *fm_76 = buffer.data(fm + 76 * ncomps + c);
        const auto *fm_77 = buffer.data(fm + 77 * ncomps + c);
        const auto *fm_78 = buffer.data(fm + 78 * ncomps + c);
        const auto *fm_79 = buffer.data(fm + 79 * ncomps + c);
        const auto *fm_80 = buffer.data(fm + 80 * ncomps + c);
        const auto *fm_81 = buffer.data(fm + 81 * ncomps + c);
        const auto *fm_82 = buffer.data(fm + 82 * ncomps + c);
        const auto *fm_83 = buffer.data(fm + 83 * ncomps + c);
        const auto *fm_84 = buffer.data(fm + 84 * ncomps + c);
        const auto *fm_85 = buffer.data(fm + 85 * ncomps + c);
        const auto *fm_86 = buffer.data(fm + 86 * ncomps + c);
        const auto *fm_87 = buffer.data(fm + 87 * ncomps + c);
        const auto *fm_88 = buffer.data(fm + 88 * ncomps + c);
        const auto *fm_89 = buffer.data(fm + 89 * ncomps + c);
        const auto *fm_90 = buffer.data(fm + 90 * ncomps + c);
        const auto *fm_91 = buffer.data(fm + 91 * ncomps + c);
        const auto *fm_92 = buffer.data(fm + 92 * ncomps + c);
        const auto *fm_93 = buffer.data(fm + 93 * ncomps + c);
        const auto *fm_94 = buffer.data(fm + 94 * ncomps + c);
        const auto *fm_95 = buffer.data(fm + 95 * ncomps + c);
        const auto *fm_96 = buffer.data(fm + 96 * ncomps + c);
        const auto *fm_97 = buffer.data(fm + 97 * ncomps + c);
        const auto *fm_98 = buffer.data(fm + 98 * ncomps + c);
        const auto *fm_99 = buffer.data(fm + 99 * ncomps + c);
        const auto *fm_110 = buffer.data(fm + 110 * ncomps + c);
        const auto *fm_111 = buffer.data(fm + 111 * ncomps + c);
        const auto *fm_112 = buffer.data(fm + 112 * ncomps + c);
        const auto *fm_113 = buffer.data(fm + 113 * ncomps + c);
        const auto *fm_114 = buffer.data(fm + 114 * ncomps + c);
        const auto *fm_115 = buffer.data(fm + 115 * ncomps + c);
        const auto *fm_116 = buffer.data(fm + 116 * ncomps + c);
        const auto *fm_117 = buffer.data(fm + 117 * ncomps + c);
        const auto *fm_118 = buffer.data(fm + 118 * ncomps + c);
        const auto *fm_119 = buffer.data(fm + 119 * ncomps + c);
        const auto *fm_120 = buffer.data(fm + 120 * ncomps + c);
        const auto *fm_121 = buffer.data(fm + 121 * ncomps + c);
        const auto *fm_122 = buffer.data(fm + 122 * ncomps + c);
        const auto *fm_123 = buffer.data(fm + 123 * ncomps + c);
        const auto *fm_124 = buffer.data(fm + 124 * ncomps + c);
        const auto *fm_125 = buffer.data(fm + 125 * ncomps + c);
        const auto *fm_126 = buffer.data(fm + 126 * ncomps + c);
        const auto *fm_127 = buffer.data(fm + 127 * ncomps + c);
        const auto *fm_128 = buffer.data(fm + 128 * ncomps + c);
        const auto *fm_129 = buffer.data(fm + 129 * ncomps + c);
        const auto *fm_130 = buffer.data(fm + 130 * ncomps + c);
        const auto *fm_131 = buffer.data(fm + 131 * ncomps + c);
        const auto *fm_132 = buffer.data(fm + 132 * ncomps + c);
        const auto *fm_133 = buffer.data(fm + 133 * ncomps + c);
        const auto *fm_134 = buffer.data(fm + 134 * ncomps + c);
        const auto *fm_135 = buffer.data(fm + 135 * ncomps + c);
        const auto *fm_136 = buffer.data(fm + 136 * ncomps + c);
        const auto *fm_137 = buffer.data(fm + 137 * ncomps + c);
        const auto *fm_138 = buffer.data(fm + 138 * ncomps + c);
        const auto *fm_139 = buffer.data(fm + 139 * ncomps + c);
        const auto *fm_140 = buffer.data(fm + 140 * ncomps + c);
        const auto *fm_141 = buffer.data(fm + 141 * ncomps + c);
        const auto *fm_142 = buffer.data(fm + 142 * ncomps + c);
        const auto *fm_143 = buffer.data(fm + 143 * ncomps + c);
        const auto *fm_144 = buffer.data(fm + 144 * ncomps + c);
        const auto *fm_145 = buffer.data(fm + 145 * ncomps + c);
        const auto *fm_146 = buffer.data(fm + 146 * ncomps + c);
        const auto *fm_147 = buffer.data(fm + 147 * ncomps + c);
        const auto *fm_148 = buffer.data(fm + 148 * ncomps + c);
        const auto *fm_149 = buffer.data(fm + 149 * ncomps + c);
        const auto *fm_150 = buffer.data(fm + 150 * ncomps + c);
        const auto *fm_151 = buffer.data(fm + 151 * ncomps + c);
        const auto *fm_152 = buffer.data(fm + 152 * ncomps + c);
        const auto *fm_153 = buffer.data(fm + 153 * ncomps + c);
        const auto *fm_154 = buffer.data(fm + 154 * ncomps + c);
        const auto *fm_165 = buffer.data(fm + 165 * ncomps + c);
        const auto *fm_166 = buffer.data(fm + 166 * ncomps + c);
        const auto *fm_167 = buffer.data(fm + 167 * ncomps + c);
        const auto *fm_168 = buffer.data(fm + 168 * ncomps + c);
        const auto *fm_169 = buffer.data(fm + 169 * ncomps + c);
        const auto *fm_170 = buffer.data(fm + 170 * ncomps + c);
        const auto *fm_171 = buffer.data(fm + 171 * ncomps + c);
        const auto *fm_172 = buffer.data(fm + 172 * ncomps + c);
        const auto *fm_173 = buffer.data(fm + 173 * ncomps + c);
        const auto *fm_174 = buffer.data(fm + 174 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fl_0, fl_1, fl_2, fl_3, fl_4, fm_0, \
                         fm_1, fm_2, fm_3, fm_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * fl_0[k]
                     + fm_0[k];

            t_1[k] = -ab_x[k] * fl_1[k]
                     + fm_1[k];

            t_2[k] = -ab_x[k] * fl_2[k]
                     + fm_2[k];

            t_3[k] = -ab_x[k] * fl_3[k]
                     + fm_3[k];

            t_4[k] = -ab_x[k] * fl_4[k]
                     + fm_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fl_5, fl_6, fl_7, fl_8, fl_9, fm_5, \
                         fm_6, fm_7, fm_8, fm_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * fl_5[k]
                     + fm_5[k];

            t_6[k] = -ab_x[k] * fl_6[k]
                     + fm_6[k];

            t_7[k] = -ab_x[k] * fl_7[k]
                     + fm_7[k];

            t_8[k] = -ab_x[k] * fl_8[k]
                     + fm_8[k];

            t_9[k] = -ab_x[k] * fl_9[k]
                     + fm_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fl_10, fl_11, fl_12, fl_13, \
                         fl_14, fm_10, fm_11, fm_12, fm_13, fm_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * fl_10[k]
                      + fm_10[k];

            t_11[k] = -ab_x[k] * fl_11[k]
                      + fm_11[k];

            t_12[k] = -ab_x[k] * fl_12[k]
                      + fm_12[k];

            t_13[k] = -ab_x[k] * fl_13[k]
                      + fm_13[k];

            t_14[k] = -ab_x[k] * fl_14[k]
                      + fm_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fl_15, fl_16, fl_17, fl_18, \
                         fl_19, fm_15, fm_16, fm_17, fm_18, fm_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * fl_15[k]
                      + fm_15[k];

            t_16[k] = -ab_x[k] * fl_16[k]
                      + fm_16[k];

            t_17[k] = -ab_x[k] * fl_17[k]
                      + fm_17[k];

            t_18[k] = -ab_x[k] * fl_18[k]
                      + fm_18[k];

            t_19[k] = -ab_x[k] * fl_19[k]
                      + fm_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fl_20, fl_21, fl_22, fl_23, \
                         fl_24, fm_20, fm_21, fm_22, fm_23, fm_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * fl_20[k]
                      + fm_20[k];

            t_21[k] = -ab_x[k] * fl_21[k]
                      + fm_21[k];

            t_22[k] = -ab_x[k] * fl_22[k]
                      + fm_22[k];

            t_23[k] = -ab_x[k] * fl_23[k]
                      + fm_23[k];

            t_24[k] = -ab_x[k] * fl_24[k]
                      + fm_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fl_25, fl_26, fl_27, fl_28, \
                         fl_29, fm_25, fm_26, fm_27, fm_28, fm_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * fl_25[k]
                      + fm_25[k];

            t_26[k] = -ab_x[k] * fl_26[k]
                      + fm_26[k];

            t_27[k] = -ab_x[k] * fl_27[k]
                      + fm_27[k];

            t_28[k] = -ab_x[k] * fl_28[k]
                      + fm_28[k];

            t_29[k] = -ab_x[k] * fl_29[k]
                      + fm_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fl_30, fl_31, fl_32, fl_33, \
                         fl_34, fm_30, fm_31, fm_32, fm_33, fm_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * fl_30[k]
                      + fm_30[k];

            t_31[k] = -ab_x[k] * fl_31[k]
                      + fm_31[k];

            t_32[k] = -ab_x[k] * fl_32[k]
                      + fm_32[k];

            t_33[k] = -ab_x[k] * fl_33[k]
                      + fm_33[k];

            t_34[k] = -ab_x[k] * fl_34[k]
                      + fm_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, fl_35, fl_36, fl_37, fl_38, \
                         fl_39, fm_35, fm_36, fm_37, fm_38, fm_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * fl_35[k]
                      + fm_35[k];

            t_36[k] = -ab_x[k] * fl_36[k]
                      + fm_36[k];

            t_37[k] = -ab_x[k] * fl_37[k]
                      + fm_37[k];

            t_38[k] = -ab_x[k] * fl_38[k]
                      + fm_38[k];

            t_39[k] = -ab_x[k] * fl_39[k]
                      + fm_39[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fl_40, fl_41, fl_42, fl_43, \
                         fl_44, fm_40, fm_41, fm_42, fm_43, fm_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * fl_40[k]
                      + fm_40[k];

            t_41[k] = -ab_x[k] * fl_41[k]
                      + fm_41[k];

            t_42[k] = -ab_x[k] * fl_42[k]
                      + fm_42[k];

            t_43[k] = -ab_x[k] * fl_43[k]
                      + fm_43[k];

            t_44[k] = -ab_x[k] * fl_44[k]
                      + fm_44[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, fl_45, fl_46, fl_47, fl_48, \
                         fl_49, fm_55, fm_56, fm_57, fm_58, fm_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * fl_45[k]
                      + fm_55[k];

            t_46[k] = -ab_x[k] * fl_46[k]
                      + fm_56[k];

            t_47[k] = -ab_x[k] * fl_47[k]
                      + fm_57[k];

            t_48[k] = -ab_x[k] * fl_48[k]
                      + fm_58[k];

            t_49[k] = -ab_x[k] * fl_49[k]
                      + fm_59[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fl_50, fl_51, fl_52, fl_53, \
                         fl_54, fm_60, fm_61, fm_62, fm_63, fm_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * fl_50[k]
                      + fm_60[k];

            t_51[k] = -ab_x[k] * fl_51[k]
                      + fm_61[k];

            t_52[k] = -ab_x[k] * fl_52[k]
                      + fm_62[k];

            t_53[k] = -ab_x[k] * fl_53[k]
                      + fm_63[k];

            t_54[k] = -ab_x[k] * fl_54[k]
                      + fm_64[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, fl_55, fl_56, fl_57, fl_58, \
                         fl_59, fm_65, fm_66, fm_67, fm_68, fm_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * fl_55[k]
                      + fm_65[k];

            t_56[k] = -ab_x[k] * fl_56[k]
                      + fm_66[k];

            t_57[k] = -ab_x[k] * fl_57[k]
                      + fm_67[k];

            t_58[k] = -ab_x[k] * fl_58[k]
                      + fm_68[k];

            t_59[k] = -ab_x[k] * fl_59[k]
                      + fm_69[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fl_60, fl_61, fl_62, fl_63, \
                         fl_64, fm_70, fm_71, fm_72, fm_73, fm_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * fl_60[k]
                      + fm_70[k];

            t_61[k] = -ab_x[k] * fl_61[k]
                      + fm_71[k];

            t_62[k] = -ab_x[k] * fl_62[k]
                      + fm_72[k];

            t_63[k] = -ab_x[k] * fl_63[k]
                      + fm_73[k];

            t_64[k] = -ab_x[k] * fl_64[k]
                      + fm_74[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, fl_65, fl_66, fl_67, fl_68, \
                         fl_69, fm_75, fm_76, fm_77, fm_78, fm_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * fl_65[k]
                      + fm_75[k];

            t_66[k] = -ab_x[k] * fl_66[k]
                      + fm_76[k];

            t_67[k] = -ab_x[k] * fl_67[k]
                      + fm_77[k];

            t_68[k] = -ab_x[k] * fl_68[k]
                      + fm_78[k];

            t_69[k] = -ab_x[k] * fl_69[k]
                      + fm_79[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fl_70, fl_71, fl_72, fl_73, \
                         fl_74, fm_80, fm_81, fm_82, fm_83, fm_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * fl_70[k]
                      + fm_80[k];

            t_71[k] = -ab_x[k] * fl_71[k]
                      + fm_81[k];

            t_72[k] = -ab_x[k] * fl_72[k]
                      + fm_82[k];

            t_73[k] = -ab_x[k] * fl_73[k]
                      + fm_83[k];

            t_74[k] = -ab_x[k] * fl_74[k]
                      + fm_84[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, fl_75, fl_76, fl_77, fl_78, \
                         fl_79, fm_85, fm_86, fm_87, fm_88, fm_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * fl_75[k]
                      + fm_85[k];

            t_76[k] = -ab_x[k] * fl_76[k]
                      + fm_86[k];

            t_77[k] = -ab_x[k] * fl_77[k]
                      + fm_87[k];

            t_78[k] = -ab_x[k] * fl_78[k]
                      + fm_88[k];

            t_79[k] = -ab_x[k] * fl_79[k]
                      + fm_89[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fl_80, fl_81, fl_82, fl_83, \
                         fl_84, fm_90, fm_91, fm_92, fm_93, fm_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * fl_80[k]
                      + fm_90[k];

            t_81[k] = -ab_x[k] * fl_81[k]
                      + fm_91[k];

            t_82[k] = -ab_x[k] * fl_82[k]
                      + fm_92[k];

            t_83[k] = -ab_x[k] * fl_83[k]
                      + fm_93[k];

            t_84[k] = -ab_x[k] * fl_84[k]
                      + fm_94[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, fl_85, fl_86, fl_87, fl_88, \
                         fl_89, fm_95, fm_96, fm_97, fm_98, fm_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * fl_85[k]
                      + fm_95[k];

            t_86[k] = -ab_x[k] * fl_86[k]
                      + fm_96[k];

            t_87[k] = -ab_x[k] * fl_87[k]
                      + fm_97[k];

            t_88[k] = -ab_x[k] * fl_88[k]
                      + fm_98[k];

            t_89[k] = -ab_x[k] * fl_89[k]
                      + fm_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fl_90, fl_91, fl_92, fl_93, \
                         fl_94, fm_110, fm_111, fm_112, fm_113, \
                         fm_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * fl_90[k]
                      + fm_110[k];

            t_91[k] = -ab_x[k] * fl_91[k]
                      + fm_111[k];

            t_92[k] = -ab_x[k] * fl_92[k]
                      + fm_112[k];

            t_93[k] = -ab_x[k] * fl_93[k]
                      + fm_113[k];

            t_94[k] = -ab_x[k] * fl_94[k]
                      + fm_114[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, fl_95, fl_96, fl_97, fl_98, \
                         fl_99, fm_115, fm_116, fm_117, fm_118, \
                         fm_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * fl_95[k]
                      + fm_115[k];

            t_96[k] = -ab_x[k] * fl_96[k]
                      + fm_116[k];

            t_97[k] = -ab_x[k] * fl_97[k]
                      + fm_117[k];

            t_98[k] = -ab_x[k] * fl_98[k]
                      + fm_118[k];

            t_99[k] = -ab_x[k] * fl_99[k]
                      + fm_119[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, fl_100, fl_101, fl_102, \
                         fl_103, fl_104, fm_120, fm_121, fm_122, fm_123, \
                         fm_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * fl_100[k]
                       + fm_120[k];

            t_101[k] = -ab_x[k] * fl_101[k]
                       + fm_121[k];

            t_102[k] = -ab_x[k] * fl_102[k]
                       + fm_122[k];

            t_103[k] = -ab_x[k] * fl_103[k]
                       + fm_123[k];

            t_104[k] = -ab_x[k] * fl_104[k]
                       + fm_124[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, fl_105, fl_106, fl_107, \
                         fl_108, fl_109, fm_125, fm_126, fm_127, fm_128, \
                         fm_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * fl_105[k]
                       + fm_125[k];

            t_106[k] = -ab_x[k] * fl_106[k]
                       + fm_126[k];

            t_107[k] = -ab_x[k] * fl_107[k]
                       + fm_127[k];

            t_108[k] = -ab_x[k] * fl_108[k]
                       + fm_128[k];

            t_109[k] = -ab_x[k] * fl_109[k]
                       + fm_129[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, fl_110, fl_111, fl_112, \
                         fl_113, fl_114, fm_130, fm_131, fm_132, fm_133, \
                         fm_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * fl_110[k]
                       + fm_130[k];

            t_111[k] = -ab_x[k] * fl_111[k]
                       + fm_131[k];

            t_112[k] = -ab_x[k] * fl_112[k]
                       + fm_132[k];

            t_113[k] = -ab_x[k] * fl_113[k]
                       + fm_133[k];

            t_114[k] = -ab_x[k] * fl_114[k]
                       + fm_134[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, fl_115, fl_116, fl_117, \
                         fl_118, fl_119, fm_135, fm_136, fm_137, fm_138, \
                         fm_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * fl_115[k]
                       + fm_135[k];

            t_116[k] = -ab_x[k] * fl_116[k]
                       + fm_136[k];

            t_117[k] = -ab_x[k] * fl_117[k]
                       + fm_137[k];

            t_118[k] = -ab_x[k] * fl_118[k]
                       + fm_138[k];

            t_119[k] = -ab_x[k] * fl_119[k]
                       + fm_139[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, fl_120, fl_121, fl_122, \
                         fl_123, fl_124, fm_140, fm_141, fm_142, fm_143, \
                         fm_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * fl_120[k]
                       + fm_140[k];

            t_121[k] = -ab_x[k] * fl_121[k]
                       + fm_141[k];

            t_122[k] = -ab_x[k] * fl_122[k]
                       + fm_142[k];

            t_123[k] = -ab_x[k] * fl_123[k]
                       + fm_143[k];

            t_124[k] = -ab_x[k] * fl_124[k]
                       + fm_144[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, fl_125, fl_126, fl_127, \
                         fl_128, fl_129, fm_145, fm_146, fm_147, fm_148, \
                         fm_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * fl_125[k]
                       + fm_145[k];

            t_126[k] = -ab_x[k] * fl_126[k]
                       + fm_146[k];

            t_127[k] = -ab_x[k] * fl_127[k]
                       + fm_147[k];

            t_128[k] = -ab_x[k] * fl_128[k]
                       + fm_148[k];

            t_129[k] = -ab_x[k] * fl_129[k]
                       + fm_149[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, fl_130, fl_131, fl_132, \
                         fl_133, fl_134, fm_150, fm_151, fm_152, fm_153, \
                         fm_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * fl_130[k]
                       + fm_150[k];

            t_131[k] = -ab_x[k] * fl_131[k]
                       + fm_151[k];

            t_132[k] = -ab_x[k] * fl_132[k]
                       + fm_152[k];

            t_133[k] = -ab_x[k] * fl_133[k]
                       + fm_153[k];

            t_134[k] = -ab_x[k] * fl_134[k]
                       + fm_154[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, fl_135, fl_136, fl_137, \
                         fl_138, fl_139, fm_165, fm_166, fm_167, fm_168, \
                         fm_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * fl_135[k]
                       + fm_165[k];

            t_136[k] = -ab_x[k] * fl_136[k]
                       + fm_166[k];

            t_137[k] = -ab_x[k] * fl_137[k]
                       + fm_167[k];

            t_138[k] = -ab_x[k] * fl_138[k]
                       + fm_168[k];

            t_139[k] = -ab_x[k] * fl_139[k]
                       + fm_169[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, fl_140, fl_141, fl_142, \
                         fl_143, fl_144, fm_170, fm_171, fm_172, fm_173, \
                         fm_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * fl_140[k]
                       + fm_170[k];

            t_141[k] = -ab_x[k] * fl_141[k]
                       + fm_171[k];

            t_142[k] = -ab_x[k] * fl_142[k]
                       + fm_172[k];

            t_143[k] = -ab_x[k] * fl_143[k]
                       + fm_173[k];

            t_144[k] = -ab_x[k] * fl_144[k]
                       + fm_174[k];
        }
    }
}

static auto
compute_hrr_gl_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fl, const size_t fm, const size_t ncomps,
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
        const auto *fl_171 = buffer.data(fl + 171 * ncomps + c);
        const auto *fl_172 = buffer.data(fl + 172 * ncomps + c);
        const auto *fl_173 = buffer.data(fl + 173 * ncomps + c);
        const auto *fl_174 = buffer.data(fl + 174 * ncomps + c);
        const auto *fl_175 = buffer.data(fl + 175 * ncomps + c);
        const auto *fl_176 = buffer.data(fl + 176 * ncomps + c);
        const auto *fl_177 = buffer.data(fl + 177 * ncomps + c);
        const auto *fl_178 = buffer.data(fl + 178 * ncomps + c);
        const auto *fl_179 = buffer.data(fl + 179 * ncomps + c);
        const auto *fl_180 = buffer.data(fl + 180 * ncomps + c);
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
        const auto *fl_216 = buffer.data(fl + 216 * ncomps + c);
        const auto *fl_217 = buffer.data(fl + 217 * ncomps + c);
        const auto *fl_218 = buffer.data(fl + 218 * ncomps + c);
        const auto *fl_219 = buffer.data(fl + 219 * ncomps + c);
        const auto *fl_220 = buffer.data(fl + 220 * ncomps + c);
        const auto *fl_221 = buffer.data(fl + 221 * ncomps + c);
        const auto *fl_222 = buffer.data(fl + 222 * ncomps + c);
        const auto *fl_223 = buffer.data(fl + 223 * ncomps + c);
        const auto *fl_224 = buffer.data(fl + 224 * ncomps + c);
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
        const auto *fl_261 = buffer.data(fl + 261 * ncomps + c);
        const auto *fl_262 = buffer.data(fl + 262 * ncomps + c);
        const auto *fl_263 = buffer.data(fl + 263 * ncomps + c);
        const auto *fl_264 = buffer.data(fl + 264 * ncomps + c);
        const auto *fl_265 = buffer.data(fl + 265 * ncomps + c);
        const auto *fl_266 = buffer.data(fl + 266 * ncomps + c);
        const auto *fl_267 = buffer.data(fl + 267 * ncomps + c);
        const auto *fl_268 = buffer.data(fl + 268 * ncomps + c);
        const auto *fl_269 = buffer.data(fl + 269 * ncomps + c);
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

        const auto *fm_175 = buffer.data(fm + 175 * ncomps + c);
        const auto *fm_176 = buffer.data(fm + 176 * ncomps + c);
        const auto *fm_177 = buffer.data(fm + 177 * ncomps + c);
        const auto *fm_178 = buffer.data(fm + 178 * ncomps + c);
        const auto *fm_179 = buffer.data(fm + 179 * ncomps + c);
        const auto *fm_180 = buffer.data(fm + 180 * ncomps + c);
        const auto *fm_181 = buffer.data(fm + 181 * ncomps + c);
        const auto *fm_182 = buffer.data(fm + 182 * ncomps + c);
        const auto *fm_183 = buffer.data(fm + 183 * ncomps + c);
        const auto *fm_184 = buffer.data(fm + 184 * ncomps + c);
        const auto *fm_185 = buffer.data(fm + 185 * ncomps + c);
        const auto *fm_186 = buffer.data(fm + 186 * ncomps + c);
        const auto *fm_187 = buffer.data(fm + 187 * ncomps + c);
        const auto *fm_188 = buffer.data(fm + 188 * ncomps + c);
        const auto *fm_189 = buffer.data(fm + 189 * ncomps + c);
        const auto *fm_190 = buffer.data(fm + 190 * ncomps + c);
        const auto *fm_191 = buffer.data(fm + 191 * ncomps + c);
        const auto *fm_192 = buffer.data(fm + 192 * ncomps + c);
        const auto *fm_193 = buffer.data(fm + 193 * ncomps + c);
        const auto *fm_194 = buffer.data(fm + 194 * ncomps + c);
        const auto *fm_195 = buffer.data(fm + 195 * ncomps + c);
        const auto *fm_196 = buffer.data(fm + 196 * ncomps + c);
        const auto *fm_197 = buffer.data(fm + 197 * ncomps + c);
        const auto *fm_198 = buffer.data(fm + 198 * ncomps + c);
        const auto *fm_199 = buffer.data(fm + 199 * ncomps + c);
        const auto *fm_200 = buffer.data(fm + 200 * ncomps + c);
        const auto *fm_201 = buffer.data(fm + 201 * ncomps + c);
        const auto *fm_202 = buffer.data(fm + 202 * ncomps + c);
        const auto *fm_203 = buffer.data(fm + 203 * ncomps + c);
        const auto *fm_204 = buffer.data(fm + 204 * ncomps + c);
        const auto *fm_205 = buffer.data(fm + 205 * ncomps + c);
        const auto *fm_206 = buffer.data(fm + 206 * ncomps + c);
        const auto *fm_207 = buffer.data(fm + 207 * ncomps + c);
        const auto *fm_208 = buffer.data(fm + 208 * ncomps + c);
        const auto *fm_209 = buffer.data(fm + 209 * ncomps + c);
        const auto *fm_220 = buffer.data(fm + 220 * ncomps + c);
        const auto *fm_221 = buffer.data(fm + 221 * ncomps + c);
        const auto *fm_222 = buffer.data(fm + 222 * ncomps + c);
        const auto *fm_223 = buffer.data(fm + 223 * ncomps + c);
        const auto *fm_224 = buffer.data(fm + 224 * ncomps + c);
        const auto *fm_225 = buffer.data(fm + 225 * ncomps + c);
        const auto *fm_226 = buffer.data(fm + 226 * ncomps + c);
        const auto *fm_227 = buffer.data(fm + 227 * ncomps + c);
        const auto *fm_228 = buffer.data(fm + 228 * ncomps + c);
        const auto *fm_229 = buffer.data(fm + 229 * ncomps + c);
        const auto *fm_230 = buffer.data(fm + 230 * ncomps + c);
        const auto *fm_231 = buffer.data(fm + 231 * ncomps + c);
        const auto *fm_232 = buffer.data(fm + 232 * ncomps + c);
        const auto *fm_233 = buffer.data(fm + 233 * ncomps + c);
        const auto *fm_234 = buffer.data(fm + 234 * ncomps + c);
        const auto *fm_235 = buffer.data(fm + 235 * ncomps + c);
        const auto *fm_236 = buffer.data(fm + 236 * ncomps + c);
        const auto *fm_237 = buffer.data(fm + 237 * ncomps + c);
        const auto *fm_238 = buffer.data(fm + 238 * ncomps + c);
        const auto *fm_239 = buffer.data(fm + 239 * ncomps + c);
        const auto *fm_240 = buffer.data(fm + 240 * ncomps + c);
        const auto *fm_241 = buffer.data(fm + 241 * ncomps + c);
        const auto *fm_242 = buffer.data(fm + 242 * ncomps + c);
        const auto *fm_243 = buffer.data(fm + 243 * ncomps + c);
        const auto *fm_244 = buffer.data(fm + 244 * ncomps + c);
        const auto *fm_245 = buffer.data(fm + 245 * ncomps + c);
        const auto *fm_246 = buffer.data(fm + 246 * ncomps + c);
        const auto *fm_247 = buffer.data(fm + 247 * ncomps + c);
        const auto *fm_248 = buffer.data(fm + 248 * ncomps + c);
        const auto *fm_249 = buffer.data(fm + 249 * ncomps + c);
        const auto *fm_250 = buffer.data(fm + 250 * ncomps + c);
        const auto *fm_251 = buffer.data(fm + 251 * ncomps + c);
        const auto *fm_252 = buffer.data(fm + 252 * ncomps + c);
        const auto *fm_253 = buffer.data(fm + 253 * ncomps + c);
        const auto *fm_254 = buffer.data(fm + 254 * ncomps + c);
        const auto *fm_255 = buffer.data(fm + 255 * ncomps + c);
        const auto *fm_256 = buffer.data(fm + 256 * ncomps + c);
        const auto *fm_257 = buffer.data(fm + 257 * ncomps + c);
        const auto *fm_258 = buffer.data(fm + 258 * ncomps + c);
        const auto *fm_259 = buffer.data(fm + 259 * ncomps + c);
        const auto *fm_260 = buffer.data(fm + 260 * ncomps + c);
        const auto *fm_261 = buffer.data(fm + 261 * ncomps + c);
        const auto *fm_262 = buffer.data(fm + 262 * ncomps + c);
        const auto *fm_263 = buffer.data(fm + 263 * ncomps + c);
        const auto *fm_264 = buffer.data(fm + 264 * ncomps + c);
        const auto *fm_275 = buffer.data(fm + 275 * ncomps + c);
        const auto *fm_276 = buffer.data(fm + 276 * ncomps + c);
        const auto *fm_277 = buffer.data(fm + 277 * ncomps + c);
        const auto *fm_278 = buffer.data(fm + 278 * ncomps + c);
        const auto *fm_279 = buffer.data(fm + 279 * ncomps + c);
        const auto *fm_280 = buffer.data(fm + 280 * ncomps + c);
        const auto *fm_281 = buffer.data(fm + 281 * ncomps + c);
        const auto *fm_282 = buffer.data(fm + 282 * ncomps + c);
        const auto *fm_283 = buffer.data(fm + 283 * ncomps + c);
        const auto *fm_284 = buffer.data(fm + 284 * ncomps + c);
        const auto *fm_285 = buffer.data(fm + 285 * ncomps + c);
        const auto *fm_286 = buffer.data(fm + 286 * ncomps + c);
        const auto *fm_287 = buffer.data(fm + 287 * ncomps + c);
        const auto *fm_288 = buffer.data(fm + 288 * ncomps + c);
        const auto *fm_289 = buffer.data(fm + 289 * ncomps + c);
        const auto *fm_290 = buffer.data(fm + 290 * ncomps + c);
        const auto *fm_291 = buffer.data(fm + 291 * ncomps + c);
        const auto *fm_292 = buffer.data(fm + 292 * ncomps + c);
        const auto *fm_293 = buffer.data(fm + 293 * ncomps + c);
        const auto *fm_294 = buffer.data(fm + 294 * ncomps + c);
        const auto *fm_295 = buffer.data(fm + 295 * ncomps + c);
        const auto *fm_296 = buffer.data(fm + 296 * ncomps + c);
        const auto *fm_297 = buffer.data(fm + 297 * ncomps + c);
        const auto *fm_298 = buffer.data(fm + 298 * ncomps + c);
        const auto *fm_299 = buffer.data(fm + 299 * ncomps + c);
        const auto *fm_300 = buffer.data(fm + 300 * ncomps + c);
        const auto *fm_301 = buffer.data(fm + 301 * ncomps + c);
        const auto *fm_302 = buffer.data(fm + 302 * ncomps + c);
        const auto *fm_303 = buffer.data(fm + 303 * ncomps + c);
        const auto *fm_304 = buffer.data(fm + 304 * ncomps + c);
        const auto *fm_305 = buffer.data(fm + 305 * ncomps + c);
        const auto *fm_306 = buffer.data(fm + 306 * ncomps + c);
        const auto *fm_307 = buffer.data(fm + 307 * ncomps + c);
        const auto *fm_308 = buffer.data(fm + 308 * ncomps + c);
        const auto *fm_309 = buffer.data(fm + 309 * ncomps + c);
        const auto *fm_310 = buffer.data(fm + 310 * ncomps + c);
        const auto *fm_311 = buffer.data(fm + 311 * ncomps + c);
        const auto *fm_312 = buffer.data(fm + 312 * ncomps + c);
        const auto *fm_313 = buffer.data(fm + 313 * ncomps + c);
        const auto *fm_314 = buffer.data(fm + 314 * ncomps + c);
        const auto *fm_315 = buffer.data(fm + 315 * ncomps + c);
        const auto *fm_316 = buffer.data(fm + 316 * ncomps + c);
        const auto *fm_317 = buffer.data(fm + 317 * ncomps + c);
        const auto *fm_318 = buffer.data(fm + 318 * ncomps + c);
        const auto *fm_319 = buffer.data(fm + 319 * ncomps + c);
        const auto *fm_330 = buffer.data(fm + 330 * ncomps + c);
        const auto *fm_331 = buffer.data(fm + 331 * ncomps + c);
        const auto *fm_332 = buffer.data(fm + 332 * ncomps + c);
        const auto *fm_333 = buffer.data(fm + 333 * ncomps + c);
        const auto *fm_334 = buffer.data(fm + 334 * ncomps + c);
        const auto *fm_335 = buffer.data(fm + 335 * ncomps + c);
        const auto *fm_336 = buffer.data(fm + 336 * ncomps + c);
        const auto *fm_337 = buffer.data(fm + 337 * ncomps + c);
        const auto *fm_338 = buffer.data(fm + 338 * ncomps + c);
        const auto *fm_339 = buffer.data(fm + 339 * ncomps + c);
        const auto *fm_340 = buffer.data(fm + 340 * ncomps + c);
        const auto *fm_341 = buffer.data(fm + 341 * ncomps + c);
        const auto *fm_342 = buffer.data(fm + 342 * ncomps + c);
        const auto *fm_343 = buffer.data(fm + 343 * ncomps + c);
        const auto *fm_344 = buffer.data(fm + 344 * ncomps + c);
        const auto *fm_345 = buffer.data(fm + 345 * ncomps + c);
        const auto *fm_346 = buffer.data(fm + 346 * ncomps + c);
        const auto *fm_347 = buffer.data(fm + 347 * ncomps + c);
        const auto *fm_348 = buffer.data(fm + 348 * ncomps + c);
        const auto *fm_349 = buffer.data(fm + 349 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, fl_145, fl_146, fl_147, \
                         fl_148, fl_149, fm_175, fm_176, fm_177, fm_178, \
                         fm_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * fl_145[k]
                       + fm_175[k];

            t_146[k] = -ab_x[k] * fl_146[k]
                       + fm_176[k];

            t_147[k] = -ab_x[k] * fl_147[k]
                       + fm_177[k];

            t_148[k] = -ab_x[k] * fl_148[k]
                       + fm_178[k];

            t_149[k] = -ab_x[k] * fl_149[k]
                       + fm_179[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, fl_150, fl_151, fl_152, \
                         fl_153, fl_154, fm_180, fm_181, fm_182, fm_183, \
                         fm_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * fl_150[k]
                       + fm_180[k];

            t_151[k] = -ab_x[k] * fl_151[k]
                       + fm_181[k];

            t_152[k] = -ab_x[k] * fl_152[k]
                       + fm_182[k];

            t_153[k] = -ab_x[k] * fl_153[k]
                       + fm_183[k];

            t_154[k] = -ab_x[k] * fl_154[k]
                       + fm_184[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, fl_155, fl_156, fl_157, \
                         fl_158, fl_159, fm_185, fm_186, fm_187, fm_188, \
                         fm_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * fl_155[k]
                       + fm_185[k];

            t_156[k] = -ab_x[k] * fl_156[k]
                       + fm_186[k];

            t_157[k] = -ab_x[k] * fl_157[k]
                       + fm_187[k];

            t_158[k] = -ab_x[k] * fl_158[k]
                       + fm_188[k];

            t_159[k] = -ab_x[k] * fl_159[k]
                       + fm_189[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, fl_160, fl_161, fl_162, \
                         fl_163, fl_164, fm_190, fm_191, fm_192, fm_193, \
                         fm_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * fl_160[k]
                       + fm_190[k];

            t_161[k] = -ab_x[k] * fl_161[k]
                       + fm_191[k];

            t_162[k] = -ab_x[k] * fl_162[k]
                       + fm_192[k];

            t_163[k] = -ab_x[k] * fl_163[k]
                       + fm_193[k];

            t_164[k] = -ab_x[k] * fl_164[k]
                       + fm_194[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, fl_165, fl_166, fl_167, \
                         fl_168, fl_169, fm_195, fm_196, fm_197, fm_198, \
                         fm_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * fl_165[k]
                       + fm_195[k];

            t_166[k] = -ab_x[k] * fl_166[k]
                       + fm_196[k];

            t_167[k] = -ab_x[k] * fl_167[k]
                       + fm_197[k];

            t_168[k] = -ab_x[k] * fl_168[k]
                       + fm_198[k];

            t_169[k] = -ab_x[k] * fl_169[k]
                       + fm_199[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, fl_170, fl_171, fl_172, \
                         fl_173, fl_174, fm_200, fm_201, fm_202, fm_203, \
                         fm_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * fl_170[k]
                       + fm_200[k];

            t_171[k] = -ab_x[k] * fl_171[k]
                       + fm_201[k];

            t_172[k] = -ab_x[k] * fl_172[k]
                       + fm_202[k];

            t_173[k] = -ab_x[k] * fl_173[k]
                       + fm_203[k];

            t_174[k] = -ab_x[k] * fl_174[k]
                       + fm_204[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, fl_175, fl_176, fl_177, \
                         fl_178, fl_179, fm_205, fm_206, fm_207, fm_208, \
                         fm_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * fl_175[k]
                       + fm_205[k];

            t_176[k] = -ab_x[k] * fl_176[k]
                       + fm_206[k];

            t_177[k] = -ab_x[k] * fl_177[k]
                       + fm_207[k];

            t_178[k] = -ab_x[k] * fl_178[k]
                       + fm_208[k];

            t_179[k] = -ab_x[k] * fl_179[k]
                       + fm_209[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, fl_180, fl_181, fl_182, \
                         fl_183, fl_184, fm_220, fm_221, fm_222, fm_223, \
                         fm_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * fl_180[k]
                       + fm_220[k];

            t_181[k] = -ab_x[k] * fl_181[k]
                       + fm_221[k];

            t_182[k] = -ab_x[k] * fl_182[k]
                       + fm_222[k];

            t_183[k] = -ab_x[k] * fl_183[k]
                       + fm_223[k];

            t_184[k] = -ab_x[k] * fl_184[k]
                       + fm_224[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, fl_185, fl_186, fl_187, \
                         fl_188, fl_189, fm_225, fm_226, fm_227, fm_228, \
                         fm_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * fl_185[k]
                       + fm_225[k];

            t_186[k] = -ab_x[k] * fl_186[k]
                       + fm_226[k];

            t_187[k] = -ab_x[k] * fl_187[k]
                       + fm_227[k];

            t_188[k] = -ab_x[k] * fl_188[k]
                       + fm_228[k];

            t_189[k] = -ab_x[k] * fl_189[k]
                       + fm_229[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, fl_190, fl_191, fl_192, \
                         fl_193, fl_194, fm_230, fm_231, fm_232, fm_233, \
                         fm_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * fl_190[k]
                       + fm_230[k];

            t_191[k] = -ab_x[k] * fl_191[k]
                       + fm_231[k];

            t_192[k] = -ab_x[k] * fl_192[k]
                       + fm_232[k];

            t_193[k] = -ab_x[k] * fl_193[k]
                       + fm_233[k];

            t_194[k] = -ab_x[k] * fl_194[k]
                       + fm_234[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, fl_195, fl_196, fl_197, \
                         fl_198, fl_199, fm_235, fm_236, fm_237, fm_238, \
                         fm_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * fl_195[k]
                       + fm_235[k];

            t_196[k] = -ab_x[k] * fl_196[k]
                       + fm_236[k];

            t_197[k] = -ab_x[k] * fl_197[k]
                       + fm_237[k];

            t_198[k] = -ab_x[k] * fl_198[k]
                       + fm_238[k];

            t_199[k] = -ab_x[k] * fl_199[k]
                       + fm_239[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, fl_200, fl_201, fl_202, \
                         fl_203, fl_204, fm_240, fm_241, fm_242, fm_243, \
                         fm_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * fl_200[k]
                       + fm_240[k];

            t_201[k] = -ab_x[k] * fl_201[k]
                       + fm_241[k];

            t_202[k] = -ab_x[k] * fl_202[k]
                       + fm_242[k];

            t_203[k] = -ab_x[k] * fl_203[k]
                       + fm_243[k];

            t_204[k] = -ab_x[k] * fl_204[k]
                       + fm_244[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, fl_205, fl_206, fl_207, \
                         fl_208, fl_209, fm_245, fm_246, fm_247, fm_248, \
                         fm_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * fl_205[k]
                       + fm_245[k];

            t_206[k] = -ab_x[k] * fl_206[k]
                       + fm_246[k];

            t_207[k] = -ab_x[k] * fl_207[k]
                       + fm_247[k];

            t_208[k] = -ab_x[k] * fl_208[k]
                       + fm_248[k];

            t_209[k] = -ab_x[k] * fl_209[k]
                       + fm_249[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, fl_210, fl_211, fl_212, \
                         fl_213, fl_214, fm_250, fm_251, fm_252, fm_253, \
                         fm_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * fl_210[k]
                       + fm_250[k];

            t_211[k] = -ab_x[k] * fl_211[k]
                       + fm_251[k];

            t_212[k] = -ab_x[k] * fl_212[k]
                       + fm_252[k];

            t_213[k] = -ab_x[k] * fl_213[k]
                       + fm_253[k];

            t_214[k] = -ab_x[k] * fl_214[k]
                       + fm_254[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, fl_215, fl_216, fl_217, \
                         fl_218, fl_219, fm_255, fm_256, fm_257, fm_258, \
                         fm_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * fl_215[k]
                       + fm_255[k];

            t_216[k] = -ab_x[k] * fl_216[k]
                       + fm_256[k];

            t_217[k] = -ab_x[k] * fl_217[k]
                       + fm_257[k];

            t_218[k] = -ab_x[k] * fl_218[k]
                       + fm_258[k];

            t_219[k] = -ab_x[k] * fl_219[k]
                       + fm_259[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, fl_220, fl_221, fl_222, \
                         fl_223, fl_224, fm_260, fm_261, fm_262, fm_263, \
                         fm_264 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_x[k] * fl_220[k]
                       + fm_260[k];

            t_221[k] = -ab_x[k] * fl_221[k]
                       + fm_261[k];

            t_222[k] = -ab_x[k] * fl_222[k]
                       + fm_262[k];

            t_223[k] = -ab_x[k] * fl_223[k]
                       + fm_263[k];

            t_224[k] = -ab_x[k] * fl_224[k]
                       + fm_264[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, fl_225, fl_226, fl_227, \
                         fl_228, fl_229, fm_275, fm_276, fm_277, fm_278, \
                         fm_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_x[k] * fl_225[k]
                       + fm_275[k];

            t_226[k] = -ab_x[k] * fl_226[k]
                       + fm_276[k];

            t_227[k] = -ab_x[k] * fl_227[k]
                       + fm_277[k];

            t_228[k] = -ab_x[k] * fl_228[k]
                       + fm_278[k];

            t_229[k] = -ab_x[k] * fl_229[k]
                       + fm_279[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, fl_230, fl_231, fl_232, \
                         fl_233, fl_234, fm_280, fm_281, fm_282, fm_283, \
                         fm_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_x[k] * fl_230[k]
                       + fm_280[k];

            t_231[k] = -ab_x[k] * fl_231[k]
                       + fm_281[k];

            t_232[k] = -ab_x[k] * fl_232[k]
                       + fm_282[k];

            t_233[k] = -ab_x[k] * fl_233[k]
                       + fm_283[k];

            t_234[k] = -ab_x[k] * fl_234[k]
                       + fm_284[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, fl_235, fl_236, fl_237, \
                         fl_238, fl_239, fm_285, fm_286, fm_287, fm_288, \
                         fm_289 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_x[k] * fl_235[k]
                       + fm_285[k];

            t_236[k] = -ab_x[k] * fl_236[k]
                       + fm_286[k];

            t_237[k] = -ab_x[k] * fl_237[k]
                       + fm_287[k];

            t_238[k] = -ab_x[k] * fl_238[k]
                       + fm_288[k];

            t_239[k] = -ab_x[k] * fl_239[k]
                       + fm_289[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, fl_240, fl_241, fl_242, \
                         fl_243, fl_244, fm_290, fm_291, fm_292, fm_293, \
                         fm_294 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_x[k] * fl_240[k]
                       + fm_290[k];

            t_241[k] = -ab_x[k] * fl_241[k]
                       + fm_291[k];

            t_242[k] = -ab_x[k] * fl_242[k]
                       + fm_292[k];

            t_243[k] = -ab_x[k] * fl_243[k]
                       + fm_293[k];

            t_244[k] = -ab_x[k] * fl_244[k]
                       + fm_294[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, fl_245, fl_246, fl_247, \
                         fl_248, fl_249, fm_295, fm_296, fm_297, fm_298, \
                         fm_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_x[k] * fl_245[k]
                       + fm_295[k];

            t_246[k] = -ab_x[k] * fl_246[k]
                       + fm_296[k];

            t_247[k] = -ab_x[k] * fl_247[k]
                       + fm_297[k];

            t_248[k] = -ab_x[k] * fl_248[k]
                       + fm_298[k];

            t_249[k] = -ab_x[k] * fl_249[k]
                       + fm_299[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, fl_250, fl_251, fl_252, \
                         fl_253, fl_254, fm_300, fm_301, fm_302, fm_303, \
                         fm_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_x[k] * fl_250[k]
                       + fm_300[k];

            t_251[k] = -ab_x[k] * fl_251[k]
                       + fm_301[k];

            t_252[k] = -ab_x[k] * fl_252[k]
                       + fm_302[k];

            t_253[k] = -ab_x[k] * fl_253[k]
                       + fm_303[k];

            t_254[k] = -ab_x[k] * fl_254[k]
                       + fm_304[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, fl_255, fl_256, fl_257, \
                         fl_258, fl_259, fm_305, fm_306, fm_307, fm_308, \
                         fm_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_x[k] * fl_255[k]
                       + fm_305[k];

            t_256[k] = -ab_x[k] * fl_256[k]
                       + fm_306[k];

            t_257[k] = -ab_x[k] * fl_257[k]
                       + fm_307[k];

            t_258[k] = -ab_x[k] * fl_258[k]
                       + fm_308[k];

            t_259[k] = -ab_x[k] * fl_259[k]
                       + fm_309[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, fl_260, fl_261, fl_262, \
                         fl_263, fl_264, fm_310, fm_311, fm_312, fm_313, \
                         fm_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_x[k] * fl_260[k]
                       + fm_310[k];

            t_261[k] = -ab_x[k] * fl_261[k]
                       + fm_311[k];

            t_262[k] = -ab_x[k] * fl_262[k]
                       + fm_312[k];

            t_263[k] = -ab_x[k] * fl_263[k]
                       + fm_313[k];

            t_264[k] = -ab_x[k] * fl_264[k]
                       + fm_314[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, fl_265, fl_266, fl_267, \
                         fl_268, fl_269, fm_315, fm_316, fm_317, fm_318, \
                         fm_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_x[k] * fl_265[k]
                       + fm_315[k];

            t_266[k] = -ab_x[k] * fl_266[k]
                       + fm_316[k];

            t_267[k] = -ab_x[k] * fl_267[k]
                       + fm_317[k];

            t_268[k] = -ab_x[k] * fl_268[k]
                       + fm_318[k];

            t_269[k] = -ab_x[k] * fl_269[k]
                       + fm_319[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, fl_270, fl_271, fl_272, \
                         fl_273, fl_274, fm_330, fm_331, fm_332, fm_333, \
                         fm_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_x[k] * fl_270[k]
                       + fm_330[k];

            t_271[k] = -ab_x[k] * fl_271[k]
                       + fm_331[k];

            t_272[k] = -ab_x[k] * fl_272[k]
                       + fm_332[k];

            t_273[k] = -ab_x[k] * fl_273[k]
                       + fm_333[k];

            t_274[k] = -ab_x[k] * fl_274[k]
                       + fm_334[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, fl_275, fl_276, fl_277, \
                         fl_278, fl_279, fm_335, fm_336, fm_337, fm_338, \
                         fm_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_x[k] * fl_275[k]
                       + fm_335[k];

            t_276[k] = -ab_x[k] * fl_276[k]
                       + fm_336[k];

            t_277[k] = -ab_x[k] * fl_277[k]
                       + fm_337[k];

            t_278[k] = -ab_x[k] * fl_278[k]
                       + fm_338[k];

            t_279[k] = -ab_x[k] * fl_279[k]
                       + fm_339[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, fl_280, fl_281, fl_282, \
                         fl_283, fl_284, fm_340, fm_341, fm_342, fm_343, \
                         fm_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_x[k] * fl_280[k]
                       + fm_340[k];

            t_281[k] = -ab_x[k] * fl_281[k]
                       + fm_341[k];

            t_282[k] = -ab_x[k] * fl_282[k]
                       + fm_342[k];

            t_283[k] = -ab_x[k] * fl_283[k]
                       + fm_343[k];

            t_284[k] = -ab_x[k] * fl_284[k]
                       + fm_344[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, fl_285, fl_286, fl_287, \
                         fl_288, fl_289, fm_345, fm_346, fm_347, fm_348, \
                         fm_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_x[k] * fl_285[k]
                       + fm_345[k];

            t_286[k] = -ab_x[k] * fl_286[k]
                       + fm_346[k];

            t_287[k] = -ab_x[k] * fl_287[k]
                       + fm_347[k];

            t_288[k] = -ab_x[k] * fl_288[k]
                       + fm_348[k];

            t_289[k] = -ab_x[k] * fl_289[k]
                       + fm_349[k];
        }
    }
}

static auto
compute_hrr_gl_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fl, const size_t fm, const size_t ncomps,
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
        const auto *fl_306 = buffer.data(fl + 306 * ncomps + c);
        const auto *fl_307 = buffer.data(fl + 307 * ncomps + c);
        const auto *fl_308 = buffer.data(fl + 308 * ncomps + c);
        const auto *fl_309 = buffer.data(fl + 309 * ncomps + c);
        const auto *fl_310 = buffer.data(fl + 310 * ncomps + c);
        const auto *fl_311 = buffer.data(fl + 311 * ncomps + c);
        const auto *fl_312 = buffer.data(fl + 312 * ncomps + c);
        const auto *fl_313 = buffer.data(fl + 313 * ncomps + c);
        const auto *fl_314 = buffer.data(fl + 314 * ncomps + c);
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
        const auto *fl_351 = buffer.data(fl + 351 * ncomps + c);
        const auto *fl_352 = buffer.data(fl + 352 * ncomps + c);
        const auto *fl_353 = buffer.data(fl + 353 * ncomps + c);
        const auto *fl_354 = buffer.data(fl + 354 * ncomps + c);
        const auto *fl_355 = buffer.data(fl + 355 * ncomps + c);
        const auto *fl_356 = buffer.data(fl + 356 * ncomps + c);
        const auto *fl_357 = buffer.data(fl + 357 * ncomps + c);
        const auto *fl_358 = buffer.data(fl + 358 * ncomps + c);
        const auto *fl_359 = buffer.data(fl + 359 * ncomps + c);
        const auto *fl_360 = buffer.data(fl + 360 * ncomps + c);
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
        const auto *fl_396 = buffer.data(fl + 396 * ncomps + c);
        const auto *fl_397 = buffer.data(fl + 397 * ncomps + c);
        const auto *fl_398 = buffer.data(fl + 398 * ncomps + c);
        const auto *fl_399 = buffer.data(fl + 399 * ncomps + c);
        const auto *fl_400 = buffer.data(fl + 400 * ncomps + c);
        const auto *fl_401 = buffer.data(fl + 401 * ncomps + c);
        const auto *fl_402 = buffer.data(fl + 402 * ncomps + c);
        const auto *fl_403 = buffer.data(fl + 403 * ncomps + c);
        const auto *fl_404 = buffer.data(fl + 404 * ncomps + c);
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

        const auto *fm_350 = buffer.data(fm + 350 * ncomps + c);
        const auto *fm_351 = buffer.data(fm + 351 * ncomps + c);
        const auto *fm_352 = buffer.data(fm + 352 * ncomps + c);
        const auto *fm_353 = buffer.data(fm + 353 * ncomps + c);
        const auto *fm_354 = buffer.data(fm + 354 * ncomps + c);
        const auto *fm_355 = buffer.data(fm + 355 * ncomps + c);
        const auto *fm_356 = buffer.data(fm + 356 * ncomps + c);
        const auto *fm_357 = buffer.data(fm + 357 * ncomps + c);
        const auto *fm_358 = buffer.data(fm + 358 * ncomps + c);
        const auto *fm_359 = buffer.data(fm + 359 * ncomps + c);
        const auto *fm_360 = buffer.data(fm + 360 * ncomps + c);
        const auto *fm_361 = buffer.data(fm + 361 * ncomps + c);
        const auto *fm_362 = buffer.data(fm + 362 * ncomps + c);
        const auto *fm_363 = buffer.data(fm + 363 * ncomps + c);
        const auto *fm_364 = buffer.data(fm + 364 * ncomps + c);
        const auto *fm_365 = buffer.data(fm + 365 * ncomps + c);
        const auto *fm_366 = buffer.data(fm + 366 * ncomps + c);
        const auto *fm_367 = buffer.data(fm + 367 * ncomps + c);
        const auto *fm_368 = buffer.data(fm + 368 * ncomps + c);
        const auto *fm_369 = buffer.data(fm + 369 * ncomps + c);
        const auto *fm_370 = buffer.data(fm + 370 * ncomps + c);
        const auto *fm_371 = buffer.data(fm + 371 * ncomps + c);
        const auto *fm_372 = buffer.data(fm + 372 * ncomps + c);
        const auto *fm_373 = buffer.data(fm + 373 * ncomps + c);
        const auto *fm_374 = buffer.data(fm + 374 * ncomps + c);
        const auto *fm_385 = buffer.data(fm + 385 * ncomps + c);
        const auto *fm_386 = buffer.data(fm + 386 * ncomps + c);
        const auto *fm_387 = buffer.data(fm + 387 * ncomps + c);
        const auto *fm_388 = buffer.data(fm + 388 * ncomps + c);
        const auto *fm_389 = buffer.data(fm + 389 * ncomps + c);
        const auto *fm_390 = buffer.data(fm + 390 * ncomps + c);
        const auto *fm_391 = buffer.data(fm + 391 * ncomps + c);
        const auto *fm_392 = buffer.data(fm + 392 * ncomps + c);
        const auto *fm_393 = buffer.data(fm + 393 * ncomps + c);
        const auto *fm_394 = buffer.data(fm + 394 * ncomps + c);
        const auto *fm_395 = buffer.data(fm + 395 * ncomps + c);
        const auto *fm_396 = buffer.data(fm + 396 * ncomps + c);
        const auto *fm_397 = buffer.data(fm + 397 * ncomps + c);
        const auto *fm_398 = buffer.data(fm + 398 * ncomps + c);
        const auto *fm_399 = buffer.data(fm + 399 * ncomps + c);
        const auto *fm_400 = buffer.data(fm + 400 * ncomps + c);
        const auto *fm_401 = buffer.data(fm + 401 * ncomps + c);
        const auto *fm_402 = buffer.data(fm + 402 * ncomps + c);
        const auto *fm_403 = buffer.data(fm + 403 * ncomps + c);
        const auto *fm_404 = buffer.data(fm + 404 * ncomps + c);
        const auto *fm_405 = buffer.data(fm + 405 * ncomps + c);
        const auto *fm_406 = buffer.data(fm + 406 * ncomps + c);
        const auto *fm_407 = buffer.data(fm + 407 * ncomps + c);
        const auto *fm_408 = buffer.data(fm + 408 * ncomps + c);
        const auto *fm_409 = buffer.data(fm + 409 * ncomps + c);
        const auto *fm_410 = buffer.data(fm + 410 * ncomps + c);
        const auto *fm_411 = buffer.data(fm + 411 * ncomps + c);
        const auto *fm_412 = buffer.data(fm + 412 * ncomps + c);
        const auto *fm_413 = buffer.data(fm + 413 * ncomps + c);
        const auto *fm_414 = buffer.data(fm + 414 * ncomps + c);
        const auto *fm_415 = buffer.data(fm + 415 * ncomps + c);
        const auto *fm_416 = buffer.data(fm + 416 * ncomps + c);
        const auto *fm_417 = buffer.data(fm + 417 * ncomps + c);
        const auto *fm_418 = buffer.data(fm + 418 * ncomps + c);
        const auto *fm_419 = buffer.data(fm + 419 * ncomps + c);
        const auto *fm_420 = buffer.data(fm + 420 * ncomps + c);
        const auto *fm_421 = buffer.data(fm + 421 * ncomps + c);
        const auto *fm_422 = buffer.data(fm + 422 * ncomps + c);
        const auto *fm_423 = buffer.data(fm + 423 * ncomps + c);
        const auto *fm_424 = buffer.data(fm + 424 * ncomps + c);
        const auto *fm_425 = buffer.data(fm + 425 * ncomps + c);
        const auto *fm_426 = buffer.data(fm + 426 * ncomps + c);
        const auto *fm_427 = buffer.data(fm + 427 * ncomps + c);
        const auto *fm_428 = buffer.data(fm + 428 * ncomps + c);
        const auto *fm_429 = buffer.data(fm + 429 * ncomps + c);
        const auto *fm_440 = buffer.data(fm + 440 * ncomps + c);
        const auto *fm_441 = buffer.data(fm + 441 * ncomps + c);
        const auto *fm_442 = buffer.data(fm + 442 * ncomps + c);
        const auto *fm_443 = buffer.data(fm + 443 * ncomps + c);
        const auto *fm_444 = buffer.data(fm + 444 * ncomps + c);
        const auto *fm_445 = buffer.data(fm + 445 * ncomps + c);
        const auto *fm_446 = buffer.data(fm + 446 * ncomps + c);
        const auto *fm_447 = buffer.data(fm + 447 * ncomps + c);
        const auto *fm_448 = buffer.data(fm + 448 * ncomps + c);
        const auto *fm_449 = buffer.data(fm + 449 * ncomps + c);
        const auto *fm_450 = buffer.data(fm + 450 * ncomps + c);
        const auto *fm_451 = buffer.data(fm + 451 * ncomps + c);
        const auto *fm_452 = buffer.data(fm + 452 * ncomps + c);
        const auto *fm_453 = buffer.data(fm + 453 * ncomps + c);
        const auto *fm_454 = buffer.data(fm + 454 * ncomps + c);
        const auto *fm_455 = buffer.data(fm + 455 * ncomps + c);
        const auto *fm_456 = buffer.data(fm + 456 * ncomps + c);
        const auto *fm_457 = buffer.data(fm + 457 * ncomps + c);
        const auto *fm_458 = buffer.data(fm + 458 * ncomps + c);
        const auto *fm_459 = buffer.data(fm + 459 * ncomps + c);
        const auto *fm_460 = buffer.data(fm + 460 * ncomps + c);
        const auto *fm_461 = buffer.data(fm + 461 * ncomps + c);
        const auto *fm_462 = buffer.data(fm + 462 * ncomps + c);
        const auto *fm_463 = buffer.data(fm + 463 * ncomps + c);
        const auto *fm_464 = buffer.data(fm + 464 * ncomps + c);
        const auto *fm_465 = buffer.data(fm + 465 * ncomps + c);
        const auto *fm_466 = buffer.data(fm + 466 * ncomps + c);
        const auto *fm_467 = buffer.data(fm + 467 * ncomps + c);
        const auto *fm_468 = buffer.data(fm + 468 * ncomps + c);
        const auto *fm_469 = buffer.data(fm + 469 * ncomps + c);
        const auto *fm_470 = buffer.data(fm + 470 * ncomps + c);
        const auto *fm_471 = buffer.data(fm + 471 * ncomps + c);
        const auto *fm_472 = buffer.data(fm + 472 * ncomps + c);
        const auto *fm_473 = buffer.data(fm + 473 * ncomps + c);
        const auto *fm_474 = buffer.data(fm + 474 * ncomps + c);
        const auto *fm_475 = buffer.data(fm + 475 * ncomps + c);
        const auto *fm_476 = buffer.data(fm + 476 * ncomps + c);
        const auto *fm_477 = buffer.data(fm + 477 * ncomps + c);
        const auto *fm_478 = buffer.data(fm + 478 * ncomps + c);
        const auto *fm_479 = buffer.data(fm + 479 * ncomps + c);
        const auto *fm_480 = buffer.data(fm + 480 * ncomps + c);
        const auto *fm_481 = buffer.data(fm + 481 * ncomps + c);
        const auto *fm_482 = buffer.data(fm + 482 * ncomps + c);
        const auto *fm_483 = buffer.data(fm + 483 * ncomps + c);
        const auto *fm_484 = buffer.data(fm + 484 * ncomps + c);
        const auto *fm_495 = buffer.data(fm + 495 * ncomps + c);
        const auto *fm_496 = buffer.data(fm + 496 * ncomps + c);
        const auto *fm_497 = buffer.data(fm + 497 * ncomps + c);
        const auto *fm_498 = buffer.data(fm + 498 * ncomps + c);
        const auto *fm_499 = buffer.data(fm + 499 * ncomps + c);
        const auto *fm_500 = buffer.data(fm + 500 * ncomps + c);
        const auto *fm_501 = buffer.data(fm + 501 * ncomps + c);
        const auto *fm_502 = buffer.data(fm + 502 * ncomps + c);
        const auto *fm_503 = buffer.data(fm + 503 * ncomps + c);
        const auto *fm_504 = buffer.data(fm + 504 * ncomps + c);
        const auto *fm_505 = buffer.data(fm + 505 * ncomps + c);
        const auto *fm_506 = buffer.data(fm + 506 * ncomps + c);
        const auto *fm_507 = buffer.data(fm + 507 * ncomps + c);
        const auto *fm_508 = buffer.data(fm + 508 * ncomps + c);
        const auto *fm_509 = buffer.data(fm + 509 * ncomps + c);
        const auto *fm_510 = buffer.data(fm + 510 * ncomps + c);
        const auto *fm_511 = buffer.data(fm + 511 * ncomps + c);
        const auto *fm_512 = buffer.data(fm + 512 * ncomps + c);
        const auto *fm_513 = buffer.data(fm + 513 * ncomps + c);
        const auto *fm_514 = buffer.data(fm + 514 * ncomps + c);
        const auto *fm_515 = buffer.data(fm + 515 * ncomps + c);
        const auto *fm_516 = buffer.data(fm + 516 * ncomps + c);
        const auto *fm_517 = buffer.data(fm + 517 * ncomps + c);
        const auto *fm_518 = buffer.data(fm + 518 * ncomps + c);
        const auto *fm_519 = buffer.data(fm + 519 * ncomps + c);
        const auto *fm_520 = buffer.data(fm + 520 * ncomps + c);
        const auto *fm_521 = buffer.data(fm + 521 * ncomps + c);
        const auto *fm_522 = buffer.data(fm + 522 * ncomps + c);
        const auto *fm_523 = buffer.data(fm + 523 * ncomps + c);
        const auto *fm_524 = buffer.data(fm + 524 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, fl_290, fl_291, fl_292, \
                         fl_293, fl_294, fm_350, fm_351, fm_352, fm_353, \
                         fm_354 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_x[k] * fl_290[k]
                       + fm_350[k];

            t_291[k] = -ab_x[k] * fl_291[k]
                       + fm_351[k];

            t_292[k] = -ab_x[k] * fl_292[k]
                       + fm_352[k];

            t_293[k] = -ab_x[k] * fl_293[k]
                       + fm_353[k];

            t_294[k] = -ab_x[k] * fl_294[k]
                       + fm_354[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, fl_295, fl_296, fl_297, \
                         fl_298, fl_299, fm_355, fm_356, fm_357, fm_358, \
                         fm_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_x[k] * fl_295[k]
                       + fm_355[k];

            t_296[k] = -ab_x[k] * fl_296[k]
                       + fm_356[k];

            t_297[k] = -ab_x[k] * fl_297[k]
                       + fm_357[k];

            t_298[k] = -ab_x[k] * fl_298[k]
                       + fm_358[k];

            t_299[k] = -ab_x[k] * fl_299[k]
                       + fm_359[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, fl_300, fl_301, fl_302, \
                         fl_303, fl_304, fm_360, fm_361, fm_362, fm_363, \
                         fm_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_x[k] * fl_300[k]
                       + fm_360[k];

            t_301[k] = -ab_x[k] * fl_301[k]
                       + fm_361[k];

            t_302[k] = -ab_x[k] * fl_302[k]
                       + fm_362[k];

            t_303[k] = -ab_x[k] * fl_303[k]
                       + fm_363[k];

            t_304[k] = -ab_x[k] * fl_304[k]
                       + fm_364[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, fl_305, fl_306, fl_307, \
                         fl_308, fl_309, fm_365, fm_366, fm_367, fm_368, \
                         fm_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_x[k] * fl_305[k]
                       + fm_365[k];

            t_306[k] = -ab_x[k] * fl_306[k]
                       + fm_366[k];

            t_307[k] = -ab_x[k] * fl_307[k]
                       + fm_367[k];

            t_308[k] = -ab_x[k] * fl_308[k]
                       + fm_368[k];

            t_309[k] = -ab_x[k] * fl_309[k]
                       + fm_369[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, fl_310, fl_311, fl_312, \
                         fl_313, fl_314, fm_370, fm_371, fm_372, fm_373, \
                         fm_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_x[k] * fl_310[k]
                       + fm_370[k];

            t_311[k] = -ab_x[k] * fl_311[k]
                       + fm_371[k];

            t_312[k] = -ab_x[k] * fl_312[k]
                       + fm_372[k];

            t_313[k] = -ab_x[k] * fl_313[k]
                       + fm_373[k];

            t_314[k] = -ab_x[k] * fl_314[k]
                       + fm_374[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, fl_315, fl_316, fl_317, \
                         fl_318, fl_319, fm_385, fm_386, fm_387, fm_388, \
                         fm_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_x[k] * fl_315[k]
                       + fm_385[k];

            t_316[k] = -ab_x[k] * fl_316[k]
                       + fm_386[k];

            t_317[k] = -ab_x[k] * fl_317[k]
                       + fm_387[k];

            t_318[k] = -ab_x[k] * fl_318[k]
                       + fm_388[k];

            t_319[k] = -ab_x[k] * fl_319[k]
                       + fm_389[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, fl_320, fl_321, fl_322, \
                         fl_323, fl_324, fm_390, fm_391, fm_392, fm_393, \
                         fm_394 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_x[k] * fl_320[k]
                       + fm_390[k];

            t_321[k] = -ab_x[k] * fl_321[k]
                       + fm_391[k];

            t_322[k] = -ab_x[k] * fl_322[k]
                       + fm_392[k];

            t_323[k] = -ab_x[k] * fl_323[k]
                       + fm_393[k];

            t_324[k] = -ab_x[k] * fl_324[k]
                       + fm_394[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, fl_325, fl_326, fl_327, \
                         fl_328, fl_329, fm_395, fm_396, fm_397, fm_398, \
                         fm_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_x[k] * fl_325[k]
                       + fm_395[k];

            t_326[k] = -ab_x[k] * fl_326[k]
                       + fm_396[k];

            t_327[k] = -ab_x[k] * fl_327[k]
                       + fm_397[k];

            t_328[k] = -ab_x[k] * fl_328[k]
                       + fm_398[k];

            t_329[k] = -ab_x[k] * fl_329[k]
                       + fm_399[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, fl_330, fl_331, fl_332, \
                         fl_333, fl_334, fm_400, fm_401, fm_402, fm_403, \
                         fm_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = -ab_x[k] * fl_330[k]
                       + fm_400[k];

            t_331[k] = -ab_x[k] * fl_331[k]
                       + fm_401[k];

            t_332[k] = -ab_x[k] * fl_332[k]
                       + fm_402[k];

            t_333[k] = -ab_x[k] * fl_333[k]
                       + fm_403[k];

            t_334[k] = -ab_x[k] * fl_334[k]
                       + fm_404[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, fl_335, fl_336, fl_337, \
                         fl_338, fl_339, fm_405, fm_406, fm_407, fm_408, \
                         fm_409 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = -ab_x[k] * fl_335[k]
                       + fm_405[k];

            t_336[k] = -ab_x[k] * fl_336[k]
                       + fm_406[k];

            t_337[k] = -ab_x[k] * fl_337[k]
                       + fm_407[k];

            t_338[k] = -ab_x[k] * fl_338[k]
                       + fm_408[k];

            t_339[k] = -ab_x[k] * fl_339[k]
                       + fm_409[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, fl_340, fl_341, fl_342, \
                         fl_343, fl_344, fm_410, fm_411, fm_412, fm_413, \
                         fm_414 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = -ab_x[k] * fl_340[k]
                       + fm_410[k];

            t_341[k] = -ab_x[k] * fl_341[k]
                       + fm_411[k];

            t_342[k] = -ab_x[k] * fl_342[k]
                       + fm_412[k];

            t_343[k] = -ab_x[k] * fl_343[k]
                       + fm_413[k];

            t_344[k] = -ab_x[k] * fl_344[k]
                       + fm_414[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, fl_345, fl_346, fl_347, \
                         fl_348, fl_349, fm_415, fm_416, fm_417, fm_418, \
                         fm_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = -ab_x[k] * fl_345[k]
                       + fm_415[k];

            t_346[k] = -ab_x[k] * fl_346[k]
                       + fm_416[k];

            t_347[k] = -ab_x[k] * fl_347[k]
                       + fm_417[k];

            t_348[k] = -ab_x[k] * fl_348[k]
                       + fm_418[k];

            t_349[k] = -ab_x[k] * fl_349[k]
                       + fm_419[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, fl_350, fl_351, fl_352, \
                         fl_353, fl_354, fm_420, fm_421, fm_422, fm_423, \
                         fm_424 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = -ab_x[k] * fl_350[k]
                       + fm_420[k];

            t_351[k] = -ab_x[k] * fl_351[k]
                       + fm_421[k];

            t_352[k] = -ab_x[k] * fl_352[k]
                       + fm_422[k];

            t_353[k] = -ab_x[k] * fl_353[k]
                       + fm_423[k];

            t_354[k] = -ab_x[k] * fl_354[k]
                       + fm_424[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, fl_355, fl_356, fl_357, \
                         fl_358, fl_359, fm_425, fm_426, fm_427, fm_428, \
                         fm_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = -ab_x[k] * fl_355[k]
                       + fm_425[k];

            t_356[k] = -ab_x[k] * fl_356[k]
                       + fm_426[k];

            t_357[k] = -ab_x[k] * fl_357[k]
                       + fm_427[k];

            t_358[k] = -ab_x[k] * fl_358[k]
                       + fm_428[k];

            t_359[k] = -ab_x[k] * fl_359[k]
                       + fm_429[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, fl_360, fl_361, fl_362, \
                         fl_363, fl_364, fm_440, fm_441, fm_442, fm_443, \
                         fm_444 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = -ab_x[k] * fl_360[k]
                       + fm_440[k];

            t_361[k] = -ab_x[k] * fl_361[k]
                       + fm_441[k];

            t_362[k] = -ab_x[k] * fl_362[k]
                       + fm_442[k];

            t_363[k] = -ab_x[k] * fl_363[k]
                       + fm_443[k];

            t_364[k] = -ab_x[k] * fl_364[k]
                       + fm_444[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, fl_365, fl_366, fl_367, \
                         fl_368, fl_369, fm_445, fm_446, fm_447, fm_448, \
                         fm_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = -ab_x[k] * fl_365[k]
                       + fm_445[k];

            t_366[k] = -ab_x[k] * fl_366[k]
                       + fm_446[k];

            t_367[k] = -ab_x[k] * fl_367[k]
                       + fm_447[k];

            t_368[k] = -ab_x[k] * fl_368[k]
                       + fm_448[k];

            t_369[k] = -ab_x[k] * fl_369[k]
                       + fm_449[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, fl_370, fl_371, fl_372, \
                         fl_373, fl_374, fm_450, fm_451, fm_452, fm_453, \
                         fm_454 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = -ab_x[k] * fl_370[k]
                       + fm_450[k];

            t_371[k] = -ab_x[k] * fl_371[k]
                       + fm_451[k];

            t_372[k] = -ab_x[k] * fl_372[k]
                       + fm_452[k];

            t_373[k] = -ab_x[k] * fl_373[k]
                       + fm_453[k];

            t_374[k] = -ab_x[k] * fl_374[k]
                       + fm_454[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, fl_375, fl_376, fl_377, \
                         fl_378, fl_379, fm_455, fm_456, fm_457, fm_458, \
                         fm_459 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = -ab_x[k] * fl_375[k]
                       + fm_455[k];

            t_376[k] = -ab_x[k] * fl_376[k]
                       + fm_456[k];

            t_377[k] = -ab_x[k] * fl_377[k]
                       + fm_457[k];

            t_378[k] = -ab_x[k] * fl_378[k]
                       + fm_458[k];

            t_379[k] = -ab_x[k] * fl_379[k]
                       + fm_459[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, fl_380, fl_381, fl_382, \
                         fl_383, fl_384, fm_460, fm_461, fm_462, fm_463, \
                         fm_464 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = -ab_x[k] * fl_380[k]
                       + fm_460[k];

            t_381[k] = -ab_x[k] * fl_381[k]
                       + fm_461[k];

            t_382[k] = -ab_x[k] * fl_382[k]
                       + fm_462[k];

            t_383[k] = -ab_x[k] * fl_383[k]
                       + fm_463[k];

            t_384[k] = -ab_x[k] * fl_384[k]
                       + fm_464[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, fl_385, fl_386, fl_387, \
                         fl_388, fl_389, fm_465, fm_466, fm_467, fm_468, \
                         fm_469 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = -ab_x[k] * fl_385[k]
                       + fm_465[k];

            t_386[k] = -ab_x[k] * fl_386[k]
                       + fm_466[k];

            t_387[k] = -ab_x[k] * fl_387[k]
                       + fm_467[k];

            t_388[k] = -ab_x[k] * fl_388[k]
                       + fm_468[k];

            t_389[k] = -ab_x[k] * fl_389[k]
                       + fm_469[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, fl_390, fl_391, fl_392, \
                         fl_393, fl_394, fm_470, fm_471, fm_472, fm_473, \
                         fm_474 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = -ab_x[k] * fl_390[k]
                       + fm_470[k];

            t_391[k] = -ab_x[k] * fl_391[k]
                       + fm_471[k];

            t_392[k] = -ab_x[k] * fl_392[k]
                       + fm_472[k];

            t_393[k] = -ab_x[k] * fl_393[k]
                       + fm_473[k];

            t_394[k] = -ab_x[k] * fl_394[k]
                       + fm_474[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, fl_395, fl_396, fl_397, \
                         fl_398, fl_399, fm_475, fm_476, fm_477, fm_478, \
                         fm_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = -ab_x[k] * fl_395[k]
                       + fm_475[k];

            t_396[k] = -ab_x[k] * fl_396[k]
                       + fm_476[k];

            t_397[k] = -ab_x[k] * fl_397[k]
                       + fm_477[k];

            t_398[k] = -ab_x[k] * fl_398[k]
                       + fm_478[k];

            t_399[k] = -ab_x[k] * fl_399[k]
                       + fm_479[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, fl_400, fl_401, fl_402, \
                         fl_403, fl_404, fm_480, fm_481, fm_482, fm_483, \
                         fm_484 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = -ab_x[k] * fl_400[k]
                       + fm_480[k];

            t_401[k] = -ab_x[k] * fl_401[k]
                       + fm_481[k];

            t_402[k] = -ab_x[k] * fl_402[k]
                       + fm_482[k];

            t_403[k] = -ab_x[k] * fl_403[k]
                       + fm_483[k];

            t_404[k] = -ab_x[k] * fl_404[k]
                       + fm_484[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, fl_405, fl_406, fl_407, \
                         fl_408, fl_409, fm_495, fm_496, fm_497, fm_498, \
                         fm_499 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = -ab_x[k] * fl_405[k]
                       + fm_495[k];

            t_406[k] = -ab_x[k] * fl_406[k]
                       + fm_496[k];

            t_407[k] = -ab_x[k] * fl_407[k]
                       + fm_497[k];

            t_408[k] = -ab_x[k] * fl_408[k]
                       + fm_498[k];

            t_409[k] = -ab_x[k] * fl_409[k]
                       + fm_499[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, fl_410, fl_411, fl_412, \
                         fl_413, fl_414, fm_500, fm_501, fm_502, fm_503, \
                         fm_504 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = -ab_x[k] * fl_410[k]
                       + fm_500[k];

            t_411[k] = -ab_x[k] * fl_411[k]
                       + fm_501[k];

            t_412[k] = -ab_x[k] * fl_412[k]
                       + fm_502[k];

            t_413[k] = -ab_x[k] * fl_413[k]
                       + fm_503[k];

            t_414[k] = -ab_x[k] * fl_414[k]
                       + fm_504[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, fl_415, fl_416, fl_417, \
                         fl_418, fl_419, fm_505, fm_506, fm_507, fm_508, \
                         fm_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = -ab_x[k] * fl_415[k]
                       + fm_505[k];

            t_416[k] = -ab_x[k] * fl_416[k]
                       + fm_506[k];

            t_417[k] = -ab_x[k] * fl_417[k]
                       + fm_507[k];

            t_418[k] = -ab_x[k] * fl_418[k]
                       + fm_508[k];

            t_419[k] = -ab_x[k] * fl_419[k]
                       + fm_509[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, fl_420, fl_421, fl_422, \
                         fl_423, fl_424, fm_510, fm_511, fm_512, fm_513, \
                         fm_514 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = -ab_x[k] * fl_420[k]
                       + fm_510[k];

            t_421[k] = -ab_x[k] * fl_421[k]
                       + fm_511[k];

            t_422[k] = -ab_x[k] * fl_422[k]
                       + fm_512[k];

            t_423[k] = -ab_x[k] * fl_423[k]
                       + fm_513[k];

            t_424[k] = -ab_x[k] * fl_424[k]
                       + fm_514[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, fl_425, fl_426, fl_427, \
                         fl_428, fl_429, fm_515, fm_516, fm_517, fm_518, \
                         fm_519 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = -ab_x[k] * fl_425[k]
                       + fm_515[k];

            t_426[k] = -ab_x[k] * fl_426[k]
                       + fm_516[k];

            t_427[k] = -ab_x[k] * fl_427[k]
                       + fm_517[k];

            t_428[k] = -ab_x[k] * fl_428[k]
                       + fm_518[k];

            t_429[k] = -ab_x[k] * fl_429[k]
                       + fm_519[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, fl_430, fl_431, fl_432, \
                         fl_433, fl_434, fm_520, fm_521, fm_522, fm_523, \
                         fm_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = -ab_x[k] * fl_430[k]
                       + fm_520[k];

            t_431[k] = -ab_x[k] * fl_431[k]
                       + fm_521[k];

            t_432[k] = -ab_x[k] * fl_432[k]
                       + fm_522[k];

            t_433[k] = -ab_x[k] * fl_433[k]
                       + fm_523[k];

            t_434[k] = -ab_x[k] * fl_434[k]
                       + fm_524[k];
        }
    }
}

static auto
compute_hrr_gl_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fl, const size_t fm, const size_t ncomps,
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
        const auto *fl_306 = buffer.data(fl + 306 * ncomps + c);
        const auto *fl_307 = buffer.data(fl + 307 * ncomps + c);
        const auto *fl_308 = buffer.data(fl + 308 * ncomps + c);
        const auto *fl_309 = buffer.data(fl + 309 * ncomps + c);
        const auto *fl_310 = buffer.data(fl + 310 * ncomps + c);
        const auto *fl_311 = buffer.data(fl + 311 * ncomps + c);
        const auto *fl_312 = buffer.data(fl + 312 * ncomps + c);
        const auto *fl_313 = buffer.data(fl + 313 * ncomps + c);
        const auto *fl_314 = buffer.data(fl + 314 * ncomps + c);
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
        const auto *fl_351 = buffer.data(fl + 351 * ncomps + c);
        const auto *fl_352 = buffer.data(fl + 352 * ncomps + c);
        const auto *fl_353 = buffer.data(fl + 353 * ncomps + c);
        const auto *fl_354 = buffer.data(fl + 354 * ncomps + c);
        const auto *fl_355 = buffer.data(fl + 355 * ncomps + c);
        const auto *fl_356 = buffer.data(fl + 356 * ncomps + c);
        const auto *fl_357 = buffer.data(fl + 357 * ncomps + c);
        const auto *fl_358 = buffer.data(fl + 358 * ncomps + c);
        const auto *fl_359 = buffer.data(fl + 359 * ncomps + c);
        const auto *fl_360 = buffer.data(fl + 360 * ncomps + c);
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
        const auto *fl_396 = buffer.data(fl + 396 * ncomps + c);
        const auto *fl_397 = buffer.data(fl + 397 * ncomps + c);
        const auto *fl_398 = buffer.data(fl + 398 * ncomps + c);
        const auto *fl_399 = buffer.data(fl + 399 * ncomps + c);
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

        const auto *fm_331 = buffer.data(fm + 331 * ncomps + c);
        const auto *fm_333 = buffer.data(fm + 333 * ncomps + c);
        const auto *fm_334 = buffer.data(fm + 334 * ncomps + c);
        const auto *fm_336 = buffer.data(fm + 336 * ncomps + c);
        const auto *fm_337 = buffer.data(fm + 337 * ncomps + c);
        const auto *fm_338 = buffer.data(fm + 338 * ncomps + c);
        const auto *fm_340 = buffer.data(fm + 340 * ncomps + c);
        const auto *fm_341 = buffer.data(fm + 341 * ncomps + c);
        const auto *fm_342 = buffer.data(fm + 342 * ncomps + c);
        const auto *fm_343 = buffer.data(fm + 343 * ncomps + c);
        const auto *fm_345 = buffer.data(fm + 345 * ncomps + c);
        const auto *fm_346 = buffer.data(fm + 346 * ncomps + c);
        const auto *fm_347 = buffer.data(fm + 347 * ncomps + c);
        const auto *fm_348 = buffer.data(fm + 348 * ncomps + c);
        const auto *fm_349 = buffer.data(fm + 349 * ncomps + c);
        const auto *fm_351 = buffer.data(fm + 351 * ncomps + c);
        const auto *fm_352 = buffer.data(fm + 352 * ncomps + c);
        const auto *fm_353 = buffer.data(fm + 353 * ncomps + c);
        const auto *fm_354 = buffer.data(fm + 354 * ncomps + c);
        const auto *fm_355 = buffer.data(fm + 355 * ncomps + c);
        const auto *fm_356 = buffer.data(fm + 356 * ncomps + c);
        const auto *fm_358 = buffer.data(fm + 358 * ncomps + c);
        const auto *fm_359 = buffer.data(fm + 359 * ncomps + c);
        const auto *fm_360 = buffer.data(fm + 360 * ncomps + c);
        const auto *fm_361 = buffer.data(fm + 361 * ncomps + c);
        const auto *fm_362 = buffer.data(fm + 362 * ncomps + c);
        const auto *fm_363 = buffer.data(fm + 363 * ncomps + c);
        const auto *fm_364 = buffer.data(fm + 364 * ncomps + c);
        const auto *fm_366 = buffer.data(fm + 366 * ncomps + c);
        const auto *fm_367 = buffer.data(fm + 367 * ncomps + c);
        const auto *fm_368 = buffer.data(fm + 368 * ncomps + c);
        const auto *fm_369 = buffer.data(fm + 369 * ncomps + c);
        const auto *fm_370 = buffer.data(fm + 370 * ncomps + c);
        const auto *fm_371 = buffer.data(fm + 371 * ncomps + c);
        const auto *fm_372 = buffer.data(fm + 372 * ncomps + c);
        const auto *fm_373 = buffer.data(fm + 373 * ncomps + c);
        const auto *fm_375 = buffer.data(fm + 375 * ncomps + c);
        const auto *fm_376 = buffer.data(fm + 376 * ncomps + c);
        const auto *fm_377 = buffer.data(fm + 377 * ncomps + c);
        const auto *fm_378 = buffer.data(fm + 378 * ncomps + c);
        const auto *fm_379 = buffer.data(fm + 379 * ncomps + c);
        const auto *fm_380 = buffer.data(fm + 380 * ncomps + c);
        const auto *fm_381 = buffer.data(fm + 381 * ncomps + c);
        const auto *fm_382 = buffer.data(fm + 382 * ncomps + c);
        const auto *fm_383 = buffer.data(fm + 383 * ncomps + c);
        const auto *fm_386 = buffer.data(fm + 386 * ncomps + c);
        const auto *fm_388 = buffer.data(fm + 388 * ncomps + c);
        const auto *fm_389 = buffer.data(fm + 389 * ncomps + c);
        const auto *fm_391 = buffer.data(fm + 391 * ncomps + c);
        const auto *fm_392 = buffer.data(fm + 392 * ncomps + c);
        const auto *fm_393 = buffer.data(fm + 393 * ncomps + c);
        const auto *fm_395 = buffer.data(fm + 395 * ncomps + c);
        const auto *fm_396 = buffer.data(fm + 396 * ncomps + c);
        const auto *fm_397 = buffer.data(fm + 397 * ncomps + c);
        const auto *fm_398 = buffer.data(fm + 398 * ncomps + c);
        const auto *fm_400 = buffer.data(fm + 400 * ncomps + c);
        const auto *fm_401 = buffer.data(fm + 401 * ncomps + c);
        const auto *fm_402 = buffer.data(fm + 402 * ncomps + c);
        const auto *fm_403 = buffer.data(fm + 403 * ncomps + c);
        const auto *fm_404 = buffer.data(fm + 404 * ncomps + c);
        const auto *fm_406 = buffer.data(fm + 406 * ncomps + c);
        const auto *fm_407 = buffer.data(fm + 407 * ncomps + c);
        const auto *fm_408 = buffer.data(fm + 408 * ncomps + c);
        const auto *fm_409 = buffer.data(fm + 409 * ncomps + c);
        const auto *fm_410 = buffer.data(fm + 410 * ncomps + c);
        const auto *fm_411 = buffer.data(fm + 411 * ncomps + c);
        const auto *fm_413 = buffer.data(fm + 413 * ncomps + c);
        const auto *fm_414 = buffer.data(fm + 414 * ncomps + c);
        const auto *fm_415 = buffer.data(fm + 415 * ncomps + c);
        const auto *fm_416 = buffer.data(fm + 416 * ncomps + c);
        const auto *fm_417 = buffer.data(fm + 417 * ncomps + c);
        const auto *fm_418 = buffer.data(fm + 418 * ncomps + c);
        const auto *fm_419 = buffer.data(fm + 419 * ncomps + c);
        const auto *fm_421 = buffer.data(fm + 421 * ncomps + c);
        const auto *fm_422 = buffer.data(fm + 422 * ncomps + c);
        const auto *fm_423 = buffer.data(fm + 423 * ncomps + c);
        const auto *fm_424 = buffer.data(fm + 424 * ncomps + c);
        const auto *fm_425 = buffer.data(fm + 425 * ncomps + c);
        const auto *fm_426 = buffer.data(fm + 426 * ncomps + c);
        const auto *fm_427 = buffer.data(fm + 427 * ncomps + c);
        const auto *fm_428 = buffer.data(fm + 428 * ncomps + c);
        const auto *fm_430 = buffer.data(fm + 430 * ncomps + c);
        const auto *fm_431 = buffer.data(fm + 431 * ncomps + c);
        const auto *fm_432 = buffer.data(fm + 432 * ncomps + c);
        const auto *fm_433 = buffer.data(fm + 433 * ncomps + c);
        const auto *fm_434 = buffer.data(fm + 434 * ncomps + c);
        const auto *fm_435 = buffer.data(fm + 435 * ncomps + c);
        const auto *fm_436 = buffer.data(fm + 436 * ncomps + c);
        const auto *fm_437 = buffer.data(fm + 437 * ncomps + c);
        const auto *fm_438 = buffer.data(fm + 438 * ncomps + c);
        const auto *fm_441 = buffer.data(fm + 441 * ncomps + c);
        const auto *fm_443 = buffer.data(fm + 443 * ncomps + c);
        const auto *fm_444 = buffer.data(fm + 444 * ncomps + c);
        const auto *fm_446 = buffer.data(fm + 446 * ncomps + c);
        const auto *fm_447 = buffer.data(fm + 447 * ncomps + c);
        const auto *fm_448 = buffer.data(fm + 448 * ncomps + c);
        const auto *fm_450 = buffer.data(fm + 450 * ncomps + c);
        const auto *fm_451 = buffer.data(fm + 451 * ncomps + c);
        const auto *fm_452 = buffer.data(fm + 452 * ncomps + c);
        const auto *fm_453 = buffer.data(fm + 453 * ncomps + c);
        const auto *fm_455 = buffer.data(fm + 455 * ncomps + c);
        const auto *fm_456 = buffer.data(fm + 456 * ncomps + c);
        const auto *fm_457 = buffer.data(fm + 457 * ncomps + c);
        const auto *fm_458 = buffer.data(fm + 458 * ncomps + c);
        const auto *fm_459 = buffer.data(fm + 459 * ncomps + c);
        const auto *fm_461 = buffer.data(fm + 461 * ncomps + c);
        const auto *fm_462 = buffer.data(fm + 462 * ncomps + c);
        const auto *fm_463 = buffer.data(fm + 463 * ncomps + c);
        const auto *fm_464 = buffer.data(fm + 464 * ncomps + c);
        const auto *fm_465 = buffer.data(fm + 465 * ncomps + c);
        const auto *fm_466 = buffer.data(fm + 466 * ncomps + c);
        const auto *fm_468 = buffer.data(fm + 468 * ncomps + c);
        const auto *fm_469 = buffer.data(fm + 469 * ncomps + c);
        const auto *fm_470 = buffer.data(fm + 470 * ncomps + c);
        const auto *fm_471 = buffer.data(fm + 471 * ncomps + c);
        const auto *fm_472 = buffer.data(fm + 472 * ncomps + c);
        const auto *fm_473 = buffer.data(fm + 473 * ncomps + c);
        const auto *fm_474 = buffer.data(fm + 474 * ncomps + c);
        const auto *fm_476 = buffer.data(fm + 476 * ncomps + c);
        const auto *fm_477 = buffer.data(fm + 477 * ncomps + c);
        const auto *fm_478 = buffer.data(fm + 478 * ncomps + c);
        const auto *fm_479 = buffer.data(fm + 479 * ncomps + c);
        const auto *fm_480 = buffer.data(fm + 480 * ncomps + c);
        const auto *fm_481 = buffer.data(fm + 481 * ncomps + c);
        const auto *fm_482 = buffer.data(fm + 482 * ncomps + c);
        const auto *fm_483 = buffer.data(fm + 483 * ncomps + c);
        const auto *fm_485 = buffer.data(fm + 485 * ncomps + c);
        const auto *fm_486 = buffer.data(fm + 486 * ncomps + c);
        const auto *fm_487 = buffer.data(fm + 487 * ncomps + c);
        const auto *fm_488 = buffer.data(fm + 488 * ncomps + c);
        const auto *fm_525 = buffer.data(fm + 525 * ncomps + c);
        const auto *fm_526 = buffer.data(fm + 526 * ncomps + c);
        const auto *fm_527 = buffer.data(fm + 527 * ncomps + c);
        const auto *fm_528 = buffer.data(fm + 528 * ncomps + c);
        const auto *fm_529 = buffer.data(fm + 529 * ncomps + c);
        const auto *fm_530 = buffer.data(fm + 530 * ncomps + c);
        const auto *fm_531 = buffer.data(fm + 531 * ncomps + c);
        const auto *fm_532 = buffer.data(fm + 532 * ncomps + c);
        const auto *fm_533 = buffer.data(fm + 533 * ncomps + c);
        const auto *fm_534 = buffer.data(fm + 534 * ncomps + c);
        const auto *fm_535 = buffer.data(fm + 535 * ncomps + c);
        const auto *fm_536 = buffer.data(fm + 536 * ncomps + c);
        const auto *fm_537 = buffer.data(fm + 537 * ncomps + c);
        const auto *fm_538 = buffer.data(fm + 538 * ncomps + c);
        const auto *fm_539 = buffer.data(fm + 539 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, fl_435, fl_436, fl_437, \
                         fl_438, fl_439, fm_525, fm_526, fm_527, fm_528, \
                         fm_529 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = -ab_x[k] * fl_435[k]
                       + fm_525[k];

            t_436[k] = -ab_x[k] * fl_436[k]
                       + fm_526[k];

            t_437[k] = -ab_x[k] * fl_437[k]
                       + fm_527[k];

            t_438[k] = -ab_x[k] * fl_438[k]
                       + fm_528[k];

            t_439[k] = -ab_x[k] * fl_439[k]
                       + fm_529[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, fl_440, fl_441, fl_442, \
                         fl_443, fl_444, fm_530, fm_531, fm_532, fm_533, \
                         fm_534 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = -ab_x[k] * fl_440[k]
                       + fm_530[k];

            t_441[k] = -ab_x[k] * fl_441[k]
                       + fm_531[k];

            t_442[k] = -ab_x[k] * fl_442[k]
                       + fm_532[k];

            t_443[k] = -ab_x[k] * fl_443[k]
                       + fm_533[k];

            t_444[k] = -ab_x[k] * fl_444[k]
                       + fm_534[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_x, fl_445, fl_446, fl_447, \
                         fl_448, fl_449, fm_535, fm_536, fm_537, fm_538, \
                         fm_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = -ab_x[k] * fl_445[k]
                       + fm_535[k];

            t_446[k] = -ab_x[k] * fl_446[k]
                       + fm_536[k];

            t_447[k] = -ab_x[k] * fl_447[k]
                       + fm_537[k];

            t_448[k] = -ab_x[k] * fl_448[k]
                       + fm_538[k];

            t_449[k] = -ab_x[k] * fl_449[k]
                       + fm_539[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_y, fl_270, fl_271, fl_272, \
                         fl_273, fl_274, fm_331, fm_333, fm_334, fm_336, \
                         fm_337 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = -ab_y[k] * fl_270[k]
                       + fm_331[k];

            t_451[k] = -ab_y[k] * fl_271[k]
                       + fm_333[k];

            t_452[k] = -ab_y[k] * fl_272[k]
                       + fm_334[k];

            t_453[k] = -ab_y[k] * fl_273[k]
                       + fm_336[k];

            t_454[k] = -ab_y[k] * fl_274[k]
                       + fm_337[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_y, fl_275, fl_276, fl_277, \
                         fl_278, fl_279, fm_338, fm_340, fm_341, fm_342, \
                         fm_343 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = -ab_y[k] * fl_275[k]
                       + fm_338[k];

            t_456[k] = -ab_y[k] * fl_276[k]
                       + fm_340[k];

            t_457[k] = -ab_y[k] * fl_277[k]
                       + fm_341[k];

            t_458[k] = -ab_y[k] * fl_278[k]
                       + fm_342[k];

            t_459[k] = -ab_y[k] * fl_279[k]
                       + fm_343[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, fl_280, fl_281, fl_282, \
                         fl_283, fl_284, fm_345, fm_346, fm_347, fm_348, \
                         fm_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = -ab_y[k] * fl_280[k]
                       + fm_345[k];

            t_461[k] = -ab_y[k] * fl_281[k]
                       + fm_346[k];

            t_462[k] = -ab_y[k] * fl_282[k]
                       + fm_347[k];

            t_463[k] = -ab_y[k] * fl_283[k]
                       + fm_348[k];

            t_464[k] = -ab_y[k] * fl_284[k]
                       + fm_349[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_y, fl_285, fl_286, fl_287, \
                         fl_288, fl_289, fm_351, fm_352, fm_353, fm_354, \
                         fm_355 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = -ab_y[k] * fl_285[k]
                       + fm_351[k];

            t_466[k] = -ab_y[k] * fl_286[k]
                       + fm_352[k];

            t_467[k] = -ab_y[k] * fl_287[k]
                       + fm_353[k];

            t_468[k] = -ab_y[k] * fl_288[k]
                       + fm_354[k];

            t_469[k] = -ab_y[k] * fl_289[k]
                       + fm_355[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_y, fl_290, fl_291, fl_292, \
                         fl_293, fl_294, fm_356, fm_358, fm_359, fm_360, \
                         fm_361 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = -ab_y[k] * fl_290[k]
                       + fm_356[k];

            t_471[k] = -ab_y[k] * fl_291[k]
                       + fm_358[k];

            t_472[k] = -ab_y[k] * fl_292[k]
                       + fm_359[k];

            t_473[k] = -ab_y[k] * fl_293[k]
                       + fm_360[k];

            t_474[k] = -ab_y[k] * fl_294[k]
                       + fm_361[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, fl_295, fl_296, fl_297, \
                         fl_298, fl_299, fm_362, fm_363, fm_364, fm_366, \
                         fm_367 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = -ab_y[k] * fl_295[k]
                       + fm_362[k];

            t_476[k] = -ab_y[k] * fl_296[k]
                       + fm_363[k];

            t_477[k] = -ab_y[k] * fl_297[k]
                       + fm_364[k];

            t_478[k] = -ab_y[k] * fl_298[k]
                       + fm_366[k];

            t_479[k] = -ab_y[k] * fl_299[k]
                       + fm_367[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_y, fl_300, fl_301, fl_302, \
                         fl_303, fl_304, fm_368, fm_369, fm_370, fm_371, \
                         fm_372 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = -ab_y[k] * fl_300[k]
                       + fm_368[k];

            t_481[k] = -ab_y[k] * fl_301[k]
                       + fm_369[k];

            t_482[k] = -ab_y[k] * fl_302[k]
                       + fm_370[k];

            t_483[k] = -ab_y[k] * fl_303[k]
                       + fm_371[k];

            t_484[k] = -ab_y[k] * fl_304[k]
                       + fm_372[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_y, fl_305, fl_306, fl_307, \
                         fl_308, fl_309, fm_373, fm_375, fm_376, fm_377, \
                         fm_378 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = -ab_y[k] * fl_305[k]
                       + fm_373[k];

            t_486[k] = -ab_y[k] * fl_306[k]
                       + fm_375[k];

            t_487[k] = -ab_y[k] * fl_307[k]
                       + fm_376[k];

            t_488[k] = -ab_y[k] * fl_308[k]
                       + fm_377[k];

            t_489[k] = -ab_y[k] * fl_309[k]
                       + fm_378[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, fl_310, fl_311, fl_312, \
                         fl_313, fl_314, fm_379, fm_380, fm_381, fm_382, \
                         fm_383 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = -ab_y[k] * fl_310[k]
                       + fm_379[k];

            t_491[k] = -ab_y[k] * fl_311[k]
                       + fm_380[k];

            t_492[k] = -ab_y[k] * fl_312[k]
                       + fm_381[k];

            t_493[k] = -ab_y[k] * fl_313[k]
                       + fm_382[k];

            t_494[k] = -ab_y[k] * fl_314[k]
                       + fm_383[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_y, fl_315, fl_316, fl_317, \
                         fl_318, fl_319, fm_386, fm_388, fm_389, fm_391, \
                         fm_392 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = -ab_y[k] * fl_315[k]
                       + fm_386[k];

            t_496[k] = -ab_y[k] * fl_316[k]
                       + fm_388[k];

            t_497[k] = -ab_y[k] * fl_317[k]
                       + fm_389[k];

            t_498[k] = -ab_y[k] * fl_318[k]
                       + fm_391[k];

            t_499[k] = -ab_y[k] * fl_319[k]
                       + fm_392[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_y, fl_320, fl_321, fl_322, \
                         fl_323, fl_324, fm_393, fm_395, fm_396, fm_397, \
                         fm_398 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = -ab_y[k] * fl_320[k]
                       + fm_393[k];

            t_501[k] = -ab_y[k] * fl_321[k]
                       + fm_395[k];

            t_502[k] = -ab_y[k] * fl_322[k]
                       + fm_396[k];

            t_503[k] = -ab_y[k] * fl_323[k]
                       + fm_397[k];

            t_504[k] = -ab_y[k] * fl_324[k]
                       + fm_398[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_y, fl_325, fl_326, fl_327, \
                         fl_328, fl_329, fm_400, fm_401, fm_402, fm_403, \
                         fm_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = -ab_y[k] * fl_325[k]
                       + fm_400[k];

            t_506[k] = -ab_y[k] * fl_326[k]
                       + fm_401[k];

            t_507[k] = -ab_y[k] * fl_327[k]
                       + fm_402[k];

            t_508[k] = -ab_y[k] * fl_328[k]
                       + fm_403[k];

            t_509[k] = -ab_y[k] * fl_329[k]
                       + fm_404[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_y, fl_330, fl_331, fl_332, \
                         fl_333, fl_334, fm_406, fm_407, fm_408, fm_409, \
                         fm_410 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = -ab_y[k] * fl_330[k]
                       + fm_406[k];

            t_511[k] = -ab_y[k] * fl_331[k]
                       + fm_407[k];

            t_512[k] = -ab_y[k] * fl_332[k]
                       + fm_408[k];

            t_513[k] = -ab_y[k] * fl_333[k]
                       + fm_409[k];

            t_514[k] = -ab_y[k] * fl_334[k]
                       + fm_410[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_y, fl_335, fl_336, fl_337, \
                         fl_338, fl_339, fm_411, fm_413, fm_414, fm_415, \
                         fm_416 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = -ab_y[k] * fl_335[k]
                       + fm_411[k];

            t_516[k] = -ab_y[k] * fl_336[k]
                       + fm_413[k];

            t_517[k] = -ab_y[k] * fl_337[k]
                       + fm_414[k];

            t_518[k] = -ab_y[k] * fl_338[k]
                       + fm_415[k];

            t_519[k] = -ab_y[k] * fl_339[k]
                       + fm_416[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_y, fl_340, fl_341, fl_342, \
                         fl_343, fl_344, fm_417, fm_418, fm_419, fm_421, \
                         fm_422 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = -ab_y[k] * fl_340[k]
                       + fm_417[k];

            t_521[k] = -ab_y[k] * fl_341[k]
                       + fm_418[k];

            t_522[k] = -ab_y[k] * fl_342[k]
                       + fm_419[k];

            t_523[k] = -ab_y[k] * fl_343[k]
                       + fm_421[k];

            t_524[k] = -ab_y[k] * fl_344[k]
                       + fm_422[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_y, fl_345, fl_346, fl_347, \
                         fl_348, fl_349, fm_423, fm_424, fm_425, fm_426, \
                         fm_427 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = -ab_y[k] * fl_345[k]
                       + fm_423[k];

            t_526[k] = -ab_y[k] * fl_346[k]
                       + fm_424[k];

            t_527[k] = -ab_y[k] * fl_347[k]
                       + fm_425[k];

            t_528[k] = -ab_y[k] * fl_348[k]
                       + fm_426[k];

            t_529[k] = -ab_y[k] * fl_349[k]
                       + fm_427[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_y, fl_350, fl_351, fl_352, \
                         fl_353, fl_354, fm_428, fm_430, fm_431, fm_432, \
                         fm_433 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = -ab_y[k] * fl_350[k]
                       + fm_428[k];

            t_531[k] = -ab_y[k] * fl_351[k]
                       + fm_430[k];

            t_532[k] = -ab_y[k] * fl_352[k]
                       + fm_431[k];

            t_533[k] = -ab_y[k] * fl_353[k]
                       + fm_432[k];

            t_534[k] = -ab_y[k] * fl_354[k]
                       + fm_433[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_y, fl_355, fl_356, fl_357, \
                         fl_358, fl_359, fm_434, fm_435, fm_436, fm_437, \
                         fm_438 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = -ab_y[k] * fl_355[k]
                       + fm_434[k];

            t_536[k] = -ab_y[k] * fl_356[k]
                       + fm_435[k];

            t_537[k] = -ab_y[k] * fl_357[k]
                       + fm_436[k];

            t_538[k] = -ab_y[k] * fl_358[k]
                       + fm_437[k];

            t_539[k] = -ab_y[k] * fl_359[k]
                       + fm_438[k];
        }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_y, fl_360, fl_361, fl_362, \
                         fl_363, fl_364, fm_441, fm_443, fm_444, fm_446, \
                         fm_447 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_540[k] = -ab_y[k] * fl_360[k]
                       + fm_441[k];

            t_541[k] = -ab_y[k] * fl_361[k]
                       + fm_443[k];

            t_542[k] = -ab_y[k] * fl_362[k]
                       + fm_444[k];

            t_543[k] = -ab_y[k] * fl_363[k]
                       + fm_446[k];

            t_544[k] = -ab_y[k] * fl_364[k]
                       + fm_447[k];
        }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_y, fl_365, fl_366, fl_367, \
                         fl_368, fl_369, fm_448, fm_450, fm_451, fm_452, \
                         fm_453 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_545[k] = -ab_y[k] * fl_365[k]
                       + fm_448[k];

            t_546[k] = -ab_y[k] * fl_366[k]
                       + fm_450[k];

            t_547[k] = -ab_y[k] * fl_367[k]
                       + fm_451[k];

            t_548[k] = -ab_y[k] * fl_368[k]
                       + fm_452[k];

            t_549[k] = -ab_y[k] * fl_369[k]
                       + fm_453[k];
        }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ab_y, fl_370, fl_371, fl_372, \
                         fl_373, fl_374, fm_455, fm_456, fm_457, fm_458, \
                         fm_459 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_550[k] = -ab_y[k] * fl_370[k]
                       + fm_455[k];

            t_551[k] = -ab_y[k] * fl_371[k]
                       + fm_456[k];

            t_552[k] = -ab_y[k] * fl_372[k]
                       + fm_457[k];

            t_553[k] = -ab_y[k] * fl_373[k]
                       + fm_458[k];

            t_554[k] = -ab_y[k] * fl_374[k]
                       + fm_459[k];
        }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ab_y, fl_375, fl_376, fl_377, \
                         fl_378, fl_379, fm_461, fm_462, fm_463, fm_464, \
                         fm_465 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_555[k] = -ab_y[k] * fl_375[k]
                       + fm_461[k];

            t_556[k] = -ab_y[k] * fl_376[k]
                       + fm_462[k];

            t_557[k] = -ab_y[k] * fl_377[k]
                       + fm_463[k];

            t_558[k] = -ab_y[k] * fl_378[k]
                       + fm_464[k];

            t_559[k] = -ab_y[k] * fl_379[k]
                       + fm_465[k];
        }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_y, fl_380, fl_381, fl_382, \
                         fl_383, fl_384, fm_466, fm_468, fm_469, fm_470, \
                         fm_471 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_560[k] = -ab_y[k] * fl_380[k]
                       + fm_466[k];

            t_561[k] = -ab_y[k] * fl_381[k]
                       + fm_468[k];

            t_562[k] = -ab_y[k] * fl_382[k]
                       + fm_469[k];

            t_563[k] = -ab_y[k] * fl_383[k]
                       + fm_470[k];

            t_564[k] = -ab_y[k] * fl_384[k]
                       + fm_471[k];
        }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_y, fl_385, fl_386, fl_387, \
                         fl_388, fl_389, fm_472, fm_473, fm_474, fm_476, \
                         fm_477 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_565[k] = -ab_y[k] * fl_385[k]
                       + fm_472[k];

            t_566[k] = -ab_y[k] * fl_386[k]
                       + fm_473[k];

            t_567[k] = -ab_y[k] * fl_387[k]
                       + fm_474[k];

            t_568[k] = -ab_y[k] * fl_388[k]
                       + fm_476[k];

            t_569[k] = -ab_y[k] * fl_389[k]
                       + fm_477[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_y, fl_390, fl_391, fl_392, \
                         fl_393, fl_394, fm_478, fm_479, fm_480, fm_481, \
                         fm_482 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = -ab_y[k] * fl_390[k]
                       + fm_478[k];

            t_571[k] = -ab_y[k] * fl_391[k]
                       + fm_479[k];

            t_572[k] = -ab_y[k] * fl_392[k]
                       + fm_480[k];

            t_573[k] = -ab_y[k] * fl_393[k]
                       + fm_481[k];

            t_574[k] = -ab_y[k] * fl_394[k]
                       + fm_482[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_y, fl_395, fl_396, fl_397, \
                         fl_398, fl_399, fm_483, fm_485, fm_486, fm_487, \
                         fm_488 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = -ab_y[k] * fl_395[k]
                       + fm_483[k];

            t_576[k] = -ab_y[k] * fl_396[k]
                       + fm_485[k];

            t_577[k] = -ab_y[k] * fl_397[k]
                       + fm_486[k];

            t_578[k] = -ab_y[k] * fl_398[k]
                       + fm_487[k];

            t_579[k] = -ab_y[k] * fl_399[k]
                       + fm_488[k];
        }
    }
}

static auto
compute_hrr_gl_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fl, const size_t fm, const size_t ncomps,
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fl_400 = buffer.data(fl + 400 * ncomps + c);
        const auto *fl_401 = buffer.data(fl + 401 * ncomps + c);
        const auto *fl_402 = buffer.data(fl + 402 * ncomps + c);
        const auto *fl_403 = buffer.data(fl + 403 * ncomps + c);
        const auto *fl_404 = buffer.data(fl + 404 * ncomps + c);
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
        const auto *fl_441 = buffer.data(fl + 441 * ncomps + c);
        const auto *fl_442 = buffer.data(fl + 442 * ncomps + c);
        const auto *fl_443 = buffer.data(fl + 443 * ncomps + c);
        const auto *fl_444 = buffer.data(fl + 444 * ncomps + c);
        const auto *fl_445 = buffer.data(fl + 445 * ncomps + c);
        const auto *fl_446 = buffer.data(fl + 446 * ncomps + c);
        const auto *fl_447 = buffer.data(fl + 447 * ncomps + c);
        const auto *fl_448 = buffer.data(fl + 448 * ncomps + c);
        const auto *fl_449 = buffer.data(fl + 449 * ncomps + c);

        const auto *fm_489 = buffer.data(fm + 489 * ncomps + c);
        const auto *fm_490 = buffer.data(fm + 490 * ncomps + c);
        const auto *fm_491 = buffer.data(fm + 491 * ncomps + c);
        const auto *fm_492 = buffer.data(fm + 492 * ncomps + c);
        const auto *fm_493 = buffer.data(fm + 493 * ncomps + c);
        const auto *fm_496 = buffer.data(fm + 496 * ncomps + c);
        const auto *fm_497 = buffer.data(fm + 497 * ncomps + c);
        const auto *fm_498 = buffer.data(fm + 498 * ncomps + c);
        const auto *fm_499 = buffer.data(fm + 499 * ncomps + c);
        const auto *fm_500 = buffer.data(fm + 500 * ncomps + c);
        const auto *fm_501 = buffer.data(fm + 501 * ncomps + c);
        const auto *fm_502 = buffer.data(fm + 502 * ncomps + c);
        const auto *fm_503 = buffer.data(fm + 503 * ncomps + c);
        const auto *fm_504 = buffer.data(fm + 504 * ncomps + c);
        const auto *fm_505 = buffer.data(fm + 505 * ncomps + c);
        const auto *fm_506 = buffer.data(fm + 506 * ncomps + c);
        const auto *fm_507 = buffer.data(fm + 507 * ncomps + c);
        const auto *fm_508 = buffer.data(fm + 508 * ncomps + c);
        const auto *fm_509 = buffer.data(fm + 509 * ncomps + c);
        const auto *fm_510 = buffer.data(fm + 510 * ncomps + c);
        const auto *fm_511 = buffer.data(fm + 511 * ncomps + c);
        const auto *fm_512 = buffer.data(fm + 512 * ncomps + c);
        const auto *fm_513 = buffer.data(fm + 513 * ncomps + c);
        const auto *fm_514 = buffer.data(fm + 514 * ncomps + c);
        const auto *fm_515 = buffer.data(fm + 515 * ncomps + c);
        const auto *fm_516 = buffer.data(fm + 516 * ncomps + c);
        const auto *fm_517 = buffer.data(fm + 517 * ncomps + c);
        const auto *fm_518 = buffer.data(fm + 518 * ncomps + c);
        const auto *fm_519 = buffer.data(fm + 519 * ncomps + c);
        const auto *fm_520 = buffer.data(fm + 520 * ncomps + c);
        const auto *fm_521 = buffer.data(fm + 521 * ncomps + c);
        const auto *fm_522 = buffer.data(fm + 522 * ncomps + c);
        const auto *fm_523 = buffer.data(fm + 523 * ncomps + c);
        const auto *fm_524 = buffer.data(fm + 524 * ncomps + c);
        const auto *fm_525 = buffer.data(fm + 525 * ncomps + c);
        const auto *fm_526 = buffer.data(fm + 526 * ncomps + c);
        const auto *fm_527 = buffer.data(fm + 527 * ncomps + c);
        const auto *fm_528 = buffer.data(fm + 528 * ncomps + c);
        const auto *fm_529 = buffer.data(fm + 529 * ncomps + c);
        const auto *fm_530 = buffer.data(fm + 530 * ncomps + c);
        const auto *fm_531 = buffer.data(fm + 531 * ncomps + c);
        const auto *fm_532 = buffer.data(fm + 532 * ncomps + c);
        const auto *fm_533 = buffer.data(fm + 533 * ncomps + c);
        const auto *fm_534 = buffer.data(fm + 534 * ncomps + c);
        const auto *fm_535 = buffer.data(fm + 535 * ncomps + c);
        const auto *fm_536 = buffer.data(fm + 536 * ncomps + c);
        const auto *fm_537 = buffer.data(fm + 537 * ncomps + c);
        const auto *fm_538 = buffer.data(fm + 538 * ncomps + c);
        const auto *fm_539 = buffer.data(fm + 539 * ncomps + c);
        const auto *fm_540 = buffer.data(fm + 540 * ncomps + c);
        const auto *fm_541 = buffer.data(fm + 541 * ncomps + c);
        const auto *fm_542 = buffer.data(fm + 542 * ncomps + c);
        const auto *fm_543 = buffer.data(fm + 543 * ncomps + c);
        const auto *fm_544 = buffer.data(fm + 544 * ncomps + c);
        const auto *fm_545 = buffer.data(fm + 545 * ncomps + c);
        const auto *fm_546 = buffer.data(fm + 546 * ncomps + c);
        const auto *fm_547 = buffer.data(fm + 547 * ncomps + c);
        const auto *fm_548 = buffer.data(fm + 548 * ncomps + c);
        const auto *fm_549 = buffer.data(fm + 549 * ncomps + c);

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ab_y, fl_400, fl_401, fl_402, \
                         fl_403, fl_404, fm_489, fm_490, fm_491, fm_492, \
                         fm_493 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = -ab_y[k] * fl_400[k]
                       + fm_489[k];

            t_581[k] = -ab_y[k] * fl_401[k]
                       + fm_490[k];

            t_582[k] = -ab_y[k] * fl_402[k]
                       + fm_491[k];

            t_583[k] = -ab_y[k] * fl_403[k]
                       + fm_492[k];

            t_584[k] = -ab_y[k] * fl_404[k]
                       + fm_493[k];
        }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ab_y, fl_405, fl_406, fl_407, \
                         fl_408, fl_409, fm_496, fm_498, fm_499, fm_501, \
                         fm_502 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_585[k] = -ab_y[k] * fl_405[k]
                       + fm_496[k];

            t_586[k] = -ab_y[k] * fl_406[k]
                       + fm_498[k];

            t_587[k] = -ab_y[k] * fl_407[k]
                       + fm_499[k];

            t_588[k] = -ab_y[k] * fl_408[k]
                       + fm_501[k];

            t_589[k] = -ab_y[k] * fl_409[k]
                       + fm_502[k];
        }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ab_y, fl_410, fl_411, fl_412, \
                         fl_413, fl_414, fm_503, fm_505, fm_506, fm_507, \
                         fm_508 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_590[k] = -ab_y[k] * fl_410[k]
                       + fm_503[k];

            t_591[k] = -ab_y[k] * fl_411[k]
                       + fm_505[k];

            t_592[k] = -ab_y[k] * fl_412[k]
                       + fm_506[k];

            t_593[k] = -ab_y[k] * fl_413[k]
                       + fm_507[k];

            t_594[k] = -ab_y[k] * fl_414[k]
                       + fm_508[k];
        }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ab_y, fl_415, fl_416, fl_417, \
                         fl_418, fl_419, fm_510, fm_511, fm_512, fm_513, \
                         fm_514 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_595[k] = -ab_y[k] * fl_415[k]
                       + fm_510[k];

            t_596[k] = -ab_y[k] * fl_416[k]
                       + fm_511[k];

            t_597[k] = -ab_y[k] * fl_417[k]
                       + fm_512[k];

            t_598[k] = -ab_y[k] * fl_418[k]
                       + fm_513[k];

            t_599[k] = -ab_y[k] * fl_419[k]
                       + fm_514[k];
        }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ab_y, fl_420, fl_421, fl_422, \
                         fl_423, fl_424, fm_516, fm_517, fm_518, fm_519, \
                         fm_520 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_600[k] = -ab_y[k] * fl_420[k]
                       + fm_516[k];

            t_601[k] = -ab_y[k] * fl_421[k]
                       + fm_517[k];

            t_602[k] = -ab_y[k] * fl_422[k]
                       + fm_518[k];

            t_603[k] = -ab_y[k] * fl_423[k]
                       + fm_519[k];

            t_604[k] = -ab_y[k] * fl_424[k]
                       + fm_520[k];
        }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ab_y, fl_425, fl_426, fl_427, \
                         fl_428, fl_429, fm_521, fm_523, fm_524, fm_525, \
                         fm_526 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_605[k] = -ab_y[k] * fl_425[k]
                       + fm_521[k];

            t_606[k] = -ab_y[k] * fl_426[k]
                       + fm_523[k];

            t_607[k] = -ab_y[k] * fl_427[k]
                       + fm_524[k];

            t_608[k] = -ab_y[k] * fl_428[k]
                       + fm_525[k];

            t_609[k] = -ab_y[k] * fl_429[k]
                       + fm_526[k];
        }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ab_y, fl_430, fl_431, fl_432, \
                         fl_433, fl_434, fm_527, fm_528, fm_529, fm_531, \
                         fm_532 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_610[k] = -ab_y[k] * fl_430[k]
                       + fm_527[k];

            t_611[k] = -ab_y[k] * fl_431[k]
                       + fm_528[k];

            t_612[k] = -ab_y[k] * fl_432[k]
                       + fm_529[k];

            t_613[k] = -ab_y[k] * fl_433[k]
                       + fm_531[k];

            t_614[k] = -ab_y[k] * fl_434[k]
                       + fm_532[k];
        }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ab_y, fl_435, fl_436, fl_437, \
                         fl_438, fl_439, fm_533, fm_534, fm_535, fm_536, \
                         fm_537 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_615[k] = -ab_y[k] * fl_435[k]
                       + fm_533[k];

            t_616[k] = -ab_y[k] * fl_436[k]
                       + fm_534[k];

            t_617[k] = -ab_y[k] * fl_437[k]
                       + fm_535[k];

            t_618[k] = -ab_y[k] * fl_438[k]
                       + fm_536[k];

            t_619[k] = -ab_y[k] * fl_439[k]
                       + fm_537[k];
        }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ab_y, fl_440, fl_441, fl_442, \
                         fl_443, fl_444, fm_538, fm_540, fm_541, fm_542, \
                         fm_543 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_620[k] = -ab_y[k] * fl_440[k]
                       + fm_538[k];

            t_621[k] = -ab_y[k] * fl_441[k]
                       + fm_540[k];

            t_622[k] = -ab_y[k] * fl_442[k]
                       + fm_541[k];

            t_623[k] = -ab_y[k] * fl_443[k]
                       + fm_542[k];

            t_624[k] = -ab_y[k] * fl_444[k]
                       + fm_543[k];
        }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ab_y, fl_445, fl_446, fl_447, \
                         fl_448, fl_449, fm_544, fm_545, fm_546, fm_547, \
                         fm_548 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_625[k] = -ab_y[k] * fl_445[k]
                       + fm_544[k];

            t_626[k] = -ab_y[k] * fl_446[k]
                       + fm_545[k];

            t_627[k] = -ab_y[k] * fl_447[k]
                       + fm_546[k];

            t_628[k] = -ab_y[k] * fl_448[k]
                       + fm_547[k];

            t_629[k] = -ab_y[k] * fl_449[k]
                       + fm_548[k];
        }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ab_z, fl_405, fl_406, fl_407, \
                         fl_408, fl_409, fm_497, fm_499, fm_500, fm_502, \
                         fm_503 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_630[k] = -ab_z[k] * fl_405[k]
                       + fm_497[k];

            t_631[k] = -ab_z[k] * fl_406[k]
                       + fm_499[k];

            t_632[k] = -ab_z[k] * fl_407[k]
                       + fm_500[k];

            t_633[k] = -ab_z[k] * fl_408[k]
                       + fm_502[k];

            t_634[k] = -ab_z[k] * fl_409[k]
                       + fm_503[k];
        }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ab_z, fl_410, fl_411, fl_412, \
                         fl_413, fl_414, fm_504, fm_506, fm_507, fm_508, \
                         fm_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_635[k] = -ab_z[k] * fl_410[k]
                       + fm_504[k];

            t_636[k] = -ab_z[k] * fl_411[k]
                       + fm_506[k];

            t_637[k] = -ab_z[k] * fl_412[k]
                       + fm_507[k];

            t_638[k] = -ab_z[k] * fl_413[k]
                       + fm_508[k];

            t_639[k] = -ab_z[k] * fl_414[k]
                       + fm_509[k];
        }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ab_z, fl_415, fl_416, fl_417, \
                         fl_418, fl_419, fm_511, fm_512, fm_513, fm_514, \
                         fm_515 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_640[k] = -ab_z[k] * fl_415[k]
                       + fm_511[k];

            t_641[k] = -ab_z[k] * fl_416[k]
                       + fm_512[k];

            t_642[k] = -ab_z[k] * fl_417[k]
                       + fm_513[k];

            t_643[k] = -ab_z[k] * fl_418[k]
                       + fm_514[k];

            t_644[k] = -ab_z[k] * fl_419[k]
                       + fm_515[k];
        }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ab_z, fl_420, fl_421, fl_422, \
                         fl_423, fl_424, fm_517, fm_518, fm_519, fm_520, \
                         fm_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_645[k] = -ab_z[k] * fl_420[k]
                       + fm_517[k];

            t_646[k] = -ab_z[k] * fl_421[k]
                       + fm_518[k];

            t_647[k] = -ab_z[k] * fl_422[k]
                       + fm_519[k];

            t_648[k] = -ab_z[k] * fl_423[k]
                       + fm_520[k];

            t_649[k] = -ab_z[k] * fl_424[k]
                       + fm_521[k];
        }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ab_z, fl_425, fl_426, fl_427, \
                         fl_428, fl_429, fm_522, fm_524, fm_525, fm_526, \
                         fm_527 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_650[k] = -ab_z[k] * fl_425[k]
                       + fm_522[k];

            t_651[k] = -ab_z[k] * fl_426[k]
                       + fm_524[k];

            t_652[k] = -ab_z[k] * fl_427[k]
                       + fm_525[k];

            t_653[k] = -ab_z[k] * fl_428[k]
                       + fm_526[k];

            t_654[k] = -ab_z[k] * fl_429[k]
                       + fm_527[k];
        }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ab_z, fl_430, fl_431, fl_432, \
                         fl_433, fl_434, fm_528, fm_529, fm_530, fm_532, \
                         fm_533 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_655[k] = -ab_z[k] * fl_430[k]
                       + fm_528[k];

            t_656[k] = -ab_z[k] * fl_431[k]
                       + fm_529[k];

            t_657[k] = -ab_z[k] * fl_432[k]
                       + fm_530[k];

            t_658[k] = -ab_z[k] * fl_433[k]
                       + fm_532[k];

            t_659[k] = -ab_z[k] * fl_434[k]
                       + fm_533[k];
        }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ab_z, fl_435, fl_436, fl_437, \
                         fl_438, fl_439, fm_534, fm_535, fm_536, fm_537, \
                         fm_538 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_660[k] = -ab_z[k] * fl_435[k]
                       + fm_534[k];

            t_661[k] = -ab_z[k] * fl_436[k]
                       + fm_535[k];

            t_662[k] = -ab_z[k] * fl_437[k]
                       + fm_536[k];

            t_663[k] = -ab_z[k] * fl_438[k]
                       + fm_537[k];

            t_664[k] = -ab_z[k] * fl_439[k]
                       + fm_538[k];
        }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ab_z, fl_440, fl_441, fl_442, \
                         fl_443, fl_444, fm_539, fm_541, fm_542, fm_543, \
                         fm_544 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_665[k] = -ab_z[k] * fl_440[k]
                       + fm_539[k];

            t_666[k] = -ab_z[k] * fl_441[k]
                       + fm_541[k];

            t_667[k] = -ab_z[k] * fl_442[k]
                       + fm_542[k];

            t_668[k] = -ab_z[k] * fl_443[k]
                       + fm_543[k];

            t_669[k] = -ab_z[k] * fl_444[k]
                       + fm_544[k];
        }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ab_z, fl_445, fl_446, fl_447, \
                         fl_448, fl_449, fm_545, fm_546, fm_547, fm_548, \
                         fm_549 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_670[k] = -ab_z[k] * fl_445[k]
                       + fm_545[k];

            t_671[k] = -ab_z[k] * fl_446[k]
                       + fm_546[k];

            t_672[k] = -ab_z[k] * fl_447[k]
                       + fm_547[k];

            t_673[k] = -ab_z[k] * fl_448[k]
                       + fm_548[k];

            t_674[k] = -ab_z[k] * fl_449[k]
                       + fm_549[k];
        }
    }
}

auto
compute_hrr_gl(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fl, const size_t fm, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_gl_piece0(buffer, coordinates, target, fl, fm, ncomps, nmax);

    compute_hrr_gl_piece1(buffer, coordinates, target, fl, fm, ncomps, nmax);

    compute_hrr_gl_piece2(buffer, coordinates, target, fl, fm, ncomps, nmax);

    compute_hrr_gl_piece3(buffer, coordinates, target, fl, fm, ncomps, nmax);

    compute_hrr_gl_piece4(buffer, coordinates, target, fl, fm, ncomps, nmax);
}

}  // namespace simdtrf
