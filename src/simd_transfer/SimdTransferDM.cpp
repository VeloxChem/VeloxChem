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


#include "SimdTransferDM.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_dm_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pm, const size_t pn, const size_t ncomps,
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

        const auto *pm_0 = buffer.data(pm + 0 * ncomps + c);
        const auto *pm_1 = buffer.data(pm + 1 * ncomps + c);
        const auto *pm_2 = buffer.data(pm + 2 * ncomps + c);
        const auto *pm_3 = buffer.data(pm + 3 * ncomps + c);
        const auto *pm_4 = buffer.data(pm + 4 * ncomps + c);
        const auto *pm_5 = buffer.data(pm + 5 * ncomps + c);
        const auto *pm_6 = buffer.data(pm + 6 * ncomps + c);
        const auto *pm_7 = buffer.data(pm + 7 * ncomps + c);
        const auto *pm_8 = buffer.data(pm + 8 * ncomps + c);
        const auto *pm_9 = buffer.data(pm + 9 * ncomps + c);
        const auto *pm_10 = buffer.data(pm + 10 * ncomps + c);
        const auto *pm_11 = buffer.data(pm + 11 * ncomps + c);
        const auto *pm_12 = buffer.data(pm + 12 * ncomps + c);
        const auto *pm_13 = buffer.data(pm + 13 * ncomps + c);
        const auto *pm_14 = buffer.data(pm + 14 * ncomps + c);
        const auto *pm_15 = buffer.data(pm + 15 * ncomps + c);
        const auto *pm_16 = buffer.data(pm + 16 * ncomps + c);
        const auto *pm_17 = buffer.data(pm + 17 * ncomps + c);
        const auto *pm_18 = buffer.data(pm + 18 * ncomps + c);
        const auto *pm_19 = buffer.data(pm + 19 * ncomps + c);
        const auto *pm_20 = buffer.data(pm + 20 * ncomps + c);
        const auto *pm_21 = buffer.data(pm + 21 * ncomps + c);
        const auto *pm_22 = buffer.data(pm + 22 * ncomps + c);
        const auto *pm_23 = buffer.data(pm + 23 * ncomps + c);
        const auto *pm_24 = buffer.data(pm + 24 * ncomps + c);
        const auto *pm_25 = buffer.data(pm + 25 * ncomps + c);
        const auto *pm_26 = buffer.data(pm + 26 * ncomps + c);
        const auto *pm_27 = buffer.data(pm + 27 * ncomps + c);
        const auto *pm_28 = buffer.data(pm + 28 * ncomps + c);
        const auto *pm_29 = buffer.data(pm + 29 * ncomps + c);
        const auto *pm_30 = buffer.data(pm + 30 * ncomps + c);
        const auto *pm_31 = buffer.data(pm + 31 * ncomps + c);
        const auto *pm_32 = buffer.data(pm + 32 * ncomps + c);
        const auto *pm_33 = buffer.data(pm + 33 * ncomps + c);
        const auto *pm_34 = buffer.data(pm + 34 * ncomps + c);
        const auto *pm_35 = buffer.data(pm + 35 * ncomps + c);
        const auto *pm_36 = buffer.data(pm + 36 * ncomps + c);
        const auto *pm_37 = buffer.data(pm + 37 * ncomps + c);
        const auto *pm_38 = buffer.data(pm + 38 * ncomps + c);
        const auto *pm_39 = buffer.data(pm + 39 * ncomps + c);
        const auto *pm_40 = buffer.data(pm + 40 * ncomps + c);
        const auto *pm_41 = buffer.data(pm + 41 * ncomps + c);
        const auto *pm_42 = buffer.data(pm + 42 * ncomps + c);
        const auto *pm_43 = buffer.data(pm + 43 * ncomps + c);
        const auto *pm_44 = buffer.data(pm + 44 * ncomps + c);
        const auto *pm_45 = buffer.data(pm + 45 * ncomps + c);
        const auto *pm_46 = buffer.data(pm + 46 * ncomps + c);
        const auto *pm_47 = buffer.data(pm + 47 * ncomps + c);
        const auto *pm_48 = buffer.data(pm + 48 * ncomps + c);
        const auto *pm_49 = buffer.data(pm + 49 * ncomps + c);
        const auto *pm_50 = buffer.data(pm + 50 * ncomps + c);
        const auto *pm_51 = buffer.data(pm + 51 * ncomps + c);
        const auto *pm_52 = buffer.data(pm + 52 * ncomps + c);
        const auto *pm_53 = buffer.data(pm + 53 * ncomps + c);
        const auto *pm_54 = buffer.data(pm + 54 * ncomps + c);
        const auto *pm_55 = buffer.data(pm + 55 * ncomps + c);
        const auto *pm_56 = buffer.data(pm + 56 * ncomps + c);
        const auto *pm_57 = buffer.data(pm + 57 * ncomps + c);
        const auto *pm_58 = buffer.data(pm + 58 * ncomps + c);
        const auto *pm_59 = buffer.data(pm + 59 * ncomps + c);
        const auto *pm_60 = buffer.data(pm + 60 * ncomps + c);
        const auto *pm_61 = buffer.data(pm + 61 * ncomps + c);
        const auto *pm_62 = buffer.data(pm + 62 * ncomps + c);
        const auto *pm_63 = buffer.data(pm + 63 * ncomps + c);
        const auto *pm_64 = buffer.data(pm + 64 * ncomps + c);
        const auto *pm_65 = buffer.data(pm + 65 * ncomps + c);
        const auto *pm_66 = buffer.data(pm + 66 * ncomps + c);
        const auto *pm_67 = buffer.data(pm + 67 * ncomps + c);
        const auto *pm_68 = buffer.data(pm + 68 * ncomps + c);
        const auto *pm_69 = buffer.data(pm + 69 * ncomps + c);
        const auto *pm_70 = buffer.data(pm + 70 * ncomps + c);
        const auto *pm_71 = buffer.data(pm + 71 * ncomps + c);
        const auto *pm_72 = buffer.data(pm + 72 * ncomps + c);
        const auto *pm_73 = buffer.data(pm + 73 * ncomps + c);
        const auto *pm_74 = buffer.data(pm + 74 * ncomps + c);
        const auto *pm_75 = buffer.data(pm + 75 * ncomps + c);
        const auto *pm_76 = buffer.data(pm + 76 * ncomps + c);
        const auto *pm_77 = buffer.data(pm + 77 * ncomps + c);
        const auto *pm_78 = buffer.data(pm + 78 * ncomps + c);
        const auto *pm_79 = buffer.data(pm + 79 * ncomps + c);
        const auto *pm_80 = buffer.data(pm + 80 * ncomps + c);
        const auto *pm_81 = buffer.data(pm + 81 * ncomps + c);
        const auto *pm_82 = buffer.data(pm + 82 * ncomps + c);
        const auto *pm_83 = buffer.data(pm + 83 * ncomps + c);
        const auto *pm_84 = buffer.data(pm + 84 * ncomps + c);
        const auto *pm_85 = buffer.data(pm + 85 * ncomps + c);
        const auto *pm_86 = buffer.data(pm + 86 * ncomps + c);
        const auto *pm_87 = buffer.data(pm + 87 * ncomps + c);
        const auto *pm_88 = buffer.data(pm + 88 * ncomps + c);
        const auto *pm_89 = buffer.data(pm + 89 * ncomps + c);
        const auto *pm_90 = buffer.data(pm + 90 * ncomps + c);
        const auto *pm_91 = buffer.data(pm + 91 * ncomps + c);
        const auto *pm_92 = buffer.data(pm + 92 * ncomps + c);
        const auto *pm_93 = buffer.data(pm + 93 * ncomps + c);
        const auto *pm_94 = buffer.data(pm + 94 * ncomps + c);
        const auto *pm_95 = buffer.data(pm + 95 * ncomps + c);
        const auto *pm_96 = buffer.data(pm + 96 * ncomps + c);
        const auto *pm_97 = buffer.data(pm + 97 * ncomps + c);
        const auto *pm_98 = buffer.data(pm + 98 * ncomps + c);
        const auto *pm_99 = buffer.data(pm + 99 * ncomps + c);
        const auto *pm_100 = buffer.data(pm + 100 * ncomps + c);
        const auto *pm_101 = buffer.data(pm + 101 * ncomps + c);
        const auto *pm_102 = buffer.data(pm + 102 * ncomps + c);
        const auto *pm_103 = buffer.data(pm + 103 * ncomps + c);
        const auto *pm_104 = buffer.data(pm + 104 * ncomps + c);
        const auto *pm_105 = buffer.data(pm + 105 * ncomps + c);
        const auto *pm_106 = buffer.data(pm + 106 * ncomps + c);
        const auto *pm_107 = buffer.data(pm + 107 * ncomps + c);
        const auto *pm_108 = buffer.data(pm + 108 * ncomps + c);
        const auto *pm_109 = buffer.data(pm + 109 * ncomps + c);
        const auto *pm_110 = buffer.data(pm + 110 * ncomps + c);
        const auto *pm_111 = buffer.data(pm + 111 * ncomps + c);
        const auto *pm_112 = buffer.data(pm + 112 * ncomps + c);
        const auto *pm_113 = buffer.data(pm + 113 * ncomps + c);
        const auto *pm_114 = buffer.data(pm + 114 * ncomps + c);
        const auto *pm_115 = buffer.data(pm + 115 * ncomps + c);
        const auto *pm_116 = buffer.data(pm + 116 * ncomps + c);
        const auto *pm_117 = buffer.data(pm + 117 * ncomps + c);
        const auto *pm_118 = buffer.data(pm + 118 * ncomps + c);
        const auto *pm_119 = buffer.data(pm + 119 * ncomps + c);
        const auto *pm_120 = buffer.data(pm + 120 * ncomps + c);
        const auto *pm_121 = buffer.data(pm + 121 * ncomps + c);
        const auto *pm_122 = buffer.data(pm + 122 * ncomps + c);
        const auto *pm_123 = buffer.data(pm + 123 * ncomps + c);
        const auto *pm_124 = buffer.data(pm + 124 * ncomps + c);
        const auto *pm_125 = buffer.data(pm + 125 * ncomps + c);
        const auto *pm_126 = buffer.data(pm + 126 * ncomps + c);
        const auto *pm_127 = buffer.data(pm + 127 * ncomps + c);
        const auto *pm_128 = buffer.data(pm + 128 * ncomps + c);
        const auto *pm_129 = buffer.data(pm + 129 * ncomps + c);
        const auto *pm_130 = buffer.data(pm + 130 * ncomps + c);
        const auto *pm_131 = buffer.data(pm + 131 * ncomps + c);
        const auto *pm_132 = buffer.data(pm + 132 * ncomps + c);
        const auto *pm_133 = buffer.data(pm + 133 * ncomps + c);
        const auto *pm_134 = buffer.data(pm + 134 * ncomps + c);
        const auto *pm_135 = buffer.data(pm + 135 * ncomps + c);
        const auto *pm_136 = buffer.data(pm + 136 * ncomps + c);
        const auto *pm_137 = buffer.data(pm + 137 * ncomps + c);
        const auto *pm_138 = buffer.data(pm + 138 * ncomps + c);
        const auto *pm_139 = buffer.data(pm + 139 * ncomps + c);
        const auto *pm_140 = buffer.data(pm + 140 * ncomps + c);
        const auto *pm_141 = buffer.data(pm + 141 * ncomps + c);
        const auto *pm_142 = buffer.data(pm + 142 * ncomps + c);
        const auto *pm_143 = buffer.data(pm + 143 * ncomps + c);
        const auto *pm_144 = buffer.data(pm + 144 * ncomps + c);

        const auto *pn_0 = buffer.data(pn + 0 * ncomps + c);
        const auto *pn_1 = buffer.data(pn + 1 * ncomps + c);
        const auto *pn_2 = buffer.data(pn + 2 * ncomps + c);
        const auto *pn_3 = buffer.data(pn + 3 * ncomps + c);
        const auto *pn_4 = buffer.data(pn + 4 * ncomps + c);
        const auto *pn_5 = buffer.data(pn + 5 * ncomps + c);
        const auto *pn_6 = buffer.data(pn + 6 * ncomps + c);
        const auto *pn_7 = buffer.data(pn + 7 * ncomps + c);
        const auto *pn_8 = buffer.data(pn + 8 * ncomps + c);
        const auto *pn_9 = buffer.data(pn + 9 * ncomps + c);
        const auto *pn_10 = buffer.data(pn + 10 * ncomps + c);
        const auto *pn_11 = buffer.data(pn + 11 * ncomps + c);
        const auto *pn_12 = buffer.data(pn + 12 * ncomps + c);
        const auto *pn_13 = buffer.data(pn + 13 * ncomps + c);
        const auto *pn_14 = buffer.data(pn + 14 * ncomps + c);
        const auto *pn_15 = buffer.data(pn + 15 * ncomps + c);
        const auto *pn_16 = buffer.data(pn + 16 * ncomps + c);
        const auto *pn_17 = buffer.data(pn + 17 * ncomps + c);
        const auto *pn_18 = buffer.data(pn + 18 * ncomps + c);
        const auto *pn_19 = buffer.data(pn + 19 * ncomps + c);
        const auto *pn_20 = buffer.data(pn + 20 * ncomps + c);
        const auto *pn_21 = buffer.data(pn + 21 * ncomps + c);
        const auto *pn_22 = buffer.data(pn + 22 * ncomps + c);
        const auto *pn_23 = buffer.data(pn + 23 * ncomps + c);
        const auto *pn_24 = buffer.data(pn + 24 * ncomps + c);
        const auto *pn_25 = buffer.data(pn + 25 * ncomps + c);
        const auto *pn_26 = buffer.data(pn + 26 * ncomps + c);
        const auto *pn_27 = buffer.data(pn + 27 * ncomps + c);
        const auto *pn_28 = buffer.data(pn + 28 * ncomps + c);
        const auto *pn_29 = buffer.data(pn + 29 * ncomps + c);
        const auto *pn_30 = buffer.data(pn + 30 * ncomps + c);
        const auto *pn_31 = buffer.data(pn + 31 * ncomps + c);
        const auto *pn_32 = buffer.data(pn + 32 * ncomps + c);
        const auto *pn_33 = buffer.data(pn + 33 * ncomps + c);
        const auto *pn_34 = buffer.data(pn + 34 * ncomps + c);
        const auto *pn_35 = buffer.data(pn + 35 * ncomps + c);
        const auto *pn_36 = buffer.data(pn + 36 * ncomps + c);
        const auto *pn_37 = buffer.data(pn + 37 * ncomps + c);
        const auto *pn_38 = buffer.data(pn + 38 * ncomps + c);
        const auto *pn_39 = buffer.data(pn + 39 * ncomps + c);
        const auto *pn_40 = buffer.data(pn + 40 * ncomps + c);
        const auto *pn_41 = buffer.data(pn + 41 * ncomps + c);
        const auto *pn_42 = buffer.data(pn + 42 * ncomps + c);
        const auto *pn_43 = buffer.data(pn + 43 * ncomps + c);
        const auto *pn_44 = buffer.data(pn + 44 * ncomps + c);
        const auto *pn_45 = buffer.data(pn + 45 * ncomps + c);
        const auto *pn_46 = buffer.data(pn + 46 * ncomps + c);
        const auto *pn_47 = buffer.data(pn + 47 * ncomps + c);
        const auto *pn_48 = buffer.data(pn + 48 * ncomps + c);
        const auto *pn_49 = buffer.data(pn + 49 * ncomps + c);
        const auto *pn_50 = buffer.data(pn + 50 * ncomps + c);
        const auto *pn_51 = buffer.data(pn + 51 * ncomps + c);
        const auto *pn_52 = buffer.data(pn + 52 * ncomps + c);
        const auto *pn_53 = buffer.data(pn + 53 * ncomps + c);
        const auto *pn_54 = buffer.data(pn + 54 * ncomps + c);
        const auto *pn_66 = buffer.data(pn + 66 * ncomps + c);
        const auto *pn_67 = buffer.data(pn + 67 * ncomps + c);
        const auto *pn_68 = buffer.data(pn + 68 * ncomps + c);
        const auto *pn_69 = buffer.data(pn + 69 * ncomps + c);
        const auto *pn_70 = buffer.data(pn + 70 * ncomps + c);
        const auto *pn_71 = buffer.data(pn + 71 * ncomps + c);
        const auto *pn_72 = buffer.data(pn + 72 * ncomps + c);
        const auto *pn_73 = buffer.data(pn + 73 * ncomps + c);
        const auto *pn_74 = buffer.data(pn + 74 * ncomps + c);
        const auto *pn_75 = buffer.data(pn + 75 * ncomps + c);
        const auto *pn_76 = buffer.data(pn + 76 * ncomps + c);
        const auto *pn_77 = buffer.data(pn + 77 * ncomps + c);
        const auto *pn_78 = buffer.data(pn + 78 * ncomps + c);
        const auto *pn_79 = buffer.data(pn + 79 * ncomps + c);
        const auto *pn_80 = buffer.data(pn + 80 * ncomps + c);
        const auto *pn_81 = buffer.data(pn + 81 * ncomps + c);
        const auto *pn_82 = buffer.data(pn + 82 * ncomps + c);
        const auto *pn_83 = buffer.data(pn + 83 * ncomps + c);
        const auto *pn_84 = buffer.data(pn + 84 * ncomps + c);
        const auto *pn_85 = buffer.data(pn + 85 * ncomps + c);
        const auto *pn_86 = buffer.data(pn + 86 * ncomps + c);
        const auto *pn_87 = buffer.data(pn + 87 * ncomps + c);
        const auto *pn_88 = buffer.data(pn + 88 * ncomps + c);
        const auto *pn_89 = buffer.data(pn + 89 * ncomps + c);
        const auto *pn_90 = buffer.data(pn + 90 * ncomps + c);
        const auto *pn_91 = buffer.data(pn + 91 * ncomps + c);
        const auto *pn_92 = buffer.data(pn + 92 * ncomps + c);
        const auto *pn_93 = buffer.data(pn + 93 * ncomps + c);
        const auto *pn_94 = buffer.data(pn + 94 * ncomps + c);
        const auto *pn_95 = buffer.data(pn + 95 * ncomps + c);
        const auto *pn_96 = buffer.data(pn + 96 * ncomps + c);
        const auto *pn_97 = buffer.data(pn + 97 * ncomps + c);
        const auto *pn_98 = buffer.data(pn + 98 * ncomps + c);
        const auto *pn_99 = buffer.data(pn + 99 * ncomps + c);
        const auto *pn_100 = buffer.data(pn + 100 * ncomps + c);
        const auto *pn_101 = buffer.data(pn + 101 * ncomps + c);
        const auto *pn_102 = buffer.data(pn + 102 * ncomps + c);
        const auto *pn_103 = buffer.data(pn + 103 * ncomps + c);
        const auto *pn_104 = buffer.data(pn + 104 * ncomps + c);
        const auto *pn_105 = buffer.data(pn + 105 * ncomps + c);
        const auto *pn_106 = buffer.data(pn + 106 * ncomps + c);
        const auto *pn_107 = buffer.data(pn + 107 * ncomps + c);
        const auto *pn_108 = buffer.data(pn + 108 * ncomps + c);
        const auto *pn_109 = buffer.data(pn + 109 * ncomps + c);
        const auto *pn_110 = buffer.data(pn + 110 * ncomps + c);
        const auto *pn_111 = buffer.data(pn + 111 * ncomps + c);
        const auto *pn_112 = buffer.data(pn + 112 * ncomps + c);
        const auto *pn_113 = buffer.data(pn + 113 * ncomps + c);
        const auto *pn_114 = buffer.data(pn + 114 * ncomps + c);
        const auto *pn_115 = buffer.data(pn + 115 * ncomps + c);
        const auto *pn_116 = buffer.data(pn + 116 * ncomps + c);
        const auto *pn_117 = buffer.data(pn + 117 * ncomps + c);
        const auto *pn_118 = buffer.data(pn + 118 * ncomps + c);
        const auto *pn_119 = buffer.data(pn + 119 * ncomps + c);
        const auto *pn_120 = buffer.data(pn + 120 * ncomps + c);
        const auto *pn_132 = buffer.data(pn + 132 * ncomps + c);
        const auto *pn_133 = buffer.data(pn + 133 * ncomps + c);
        const auto *pn_134 = buffer.data(pn + 134 * ncomps + c);
        const auto *pn_135 = buffer.data(pn + 135 * ncomps + c);
        const auto *pn_136 = buffer.data(pn + 136 * ncomps + c);
        const auto *pn_137 = buffer.data(pn + 137 * ncomps + c);
        const auto *pn_138 = buffer.data(pn + 138 * ncomps + c);
        const auto *pn_139 = buffer.data(pn + 139 * ncomps + c);
        const auto *pn_140 = buffer.data(pn + 140 * ncomps + c);
        const auto *pn_141 = buffer.data(pn + 141 * ncomps + c);
        const auto *pn_142 = buffer.data(pn + 142 * ncomps + c);
        const auto *pn_143 = buffer.data(pn + 143 * ncomps + c);
        const auto *pn_144 = buffer.data(pn + 144 * ncomps + c);
        const auto *pn_145 = buffer.data(pn + 145 * ncomps + c);
        const auto *pn_146 = buffer.data(pn + 146 * ncomps + c);
        const auto *pn_147 = buffer.data(pn + 147 * ncomps + c);
        const auto *pn_148 = buffer.data(pn + 148 * ncomps + c);
        const auto *pn_149 = buffer.data(pn + 149 * ncomps + c);
        const auto *pn_150 = buffer.data(pn + 150 * ncomps + c);
        const auto *pn_151 = buffer.data(pn + 151 * ncomps + c);
        const auto *pn_152 = buffer.data(pn + 152 * ncomps + c);
        const auto *pn_153 = buffer.data(pn + 153 * ncomps + c);
        const auto *pn_154 = buffer.data(pn + 154 * ncomps + c);
        const auto *pn_155 = buffer.data(pn + 155 * ncomps + c);
        const auto *pn_156 = buffer.data(pn + 156 * ncomps + c);
        const auto *pn_157 = buffer.data(pn + 157 * ncomps + c);
        const auto *pn_158 = buffer.data(pn + 158 * ncomps + c);
        const auto *pn_159 = buffer.data(pn + 159 * ncomps + c);
        const auto *pn_160 = buffer.data(pn + 160 * ncomps + c);
        const auto *pn_161 = buffer.data(pn + 161 * ncomps + c);
        const auto *pn_162 = buffer.data(pn + 162 * ncomps + c);
        const auto *pn_163 = buffer.data(pn + 163 * ncomps + c);
        const auto *pn_164 = buffer.data(pn + 164 * ncomps + c);
        const auto *pn_165 = buffer.data(pn + 165 * ncomps + c);
        const auto *pn_166 = buffer.data(pn + 166 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pm_0, pm_1, pm_2, pm_3, pm_4, pn_0, \
                         pn_1, pn_2, pn_3, pn_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pm_0[k]
                     + pn_0[k];

            t_1[k] = -ab_x[k] * pm_1[k]
                     + pn_1[k];

            t_2[k] = -ab_x[k] * pm_2[k]
                     + pn_2[k];

            t_3[k] = -ab_x[k] * pm_3[k]
                     + pn_3[k];

            t_4[k] = -ab_x[k] * pm_4[k]
                     + pn_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pm_5, pm_6, pm_7, pm_8, pm_9, pn_5, \
                         pn_6, pn_7, pn_8, pn_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pm_5[k]
                     + pn_5[k];

            t_6[k] = -ab_x[k] * pm_6[k]
                     + pn_6[k];

            t_7[k] = -ab_x[k] * pm_7[k]
                     + pn_7[k];

            t_8[k] = -ab_x[k] * pm_8[k]
                     + pn_8[k];

            t_9[k] = -ab_x[k] * pm_9[k]
                     + pn_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pm_10, pm_11, pm_12, pm_13, \
                         pm_14, pn_10, pn_11, pn_12, pn_13, pn_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pm_10[k]
                      + pn_10[k];

            t_11[k] = -ab_x[k] * pm_11[k]
                      + pn_11[k];

            t_12[k] = -ab_x[k] * pm_12[k]
                      + pn_12[k];

            t_13[k] = -ab_x[k] * pm_13[k]
                      + pn_13[k];

            t_14[k] = -ab_x[k] * pm_14[k]
                      + pn_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pm_15, pm_16, pm_17, pm_18, \
                         pm_19, pn_15, pn_16, pn_17, pn_18, pn_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pm_15[k]
                      + pn_15[k];

            t_16[k] = -ab_x[k] * pm_16[k]
                      + pn_16[k];

            t_17[k] = -ab_x[k] * pm_17[k]
                      + pn_17[k];

            t_18[k] = -ab_x[k] * pm_18[k]
                      + pn_18[k];

            t_19[k] = -ab_x[k] * pm_19[k]
                      + pn_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pm_20, pm_21, pm_22, pm_23, \
                         pm_24, pn_20, pn_21, pn_22, pn_23, pn_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * pm_20[k]
                      + pn_20[k];

            t_21[k] = -ab_x[k] * pm_21[k]
                      + pn_21[k];

            t_22[k] = -ab_x[k] * pm_22[k]
                      + pn_22[k];

            t_23[k] = -ab_x[k] * pm_23[k]
                      + pn_23[k];

            t_24[k] = -ab_x[k] * pm_24[k]
                      + pn_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pm_25, pm_26, pm_27, pm_28, \
                         pm_29, pn_25, pn_26, pn_27, pn_28, pn_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * pm_25[k]
                      + pn_25[k];

            t_26[k] = -ab_x[k] * pm_26[k]
                      + pn_26[k];

            t_27[k] = -ab_x[k] * pm_27[k]
                      + pn_27[k];

            t_28[k] = -ab_x[k] * pm_28[k]
                      + pn_28[k];

            t_29[k] = -ab_x[k] * pm_29[k]
                      + pn_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pm_30, pm_31, pm_32, pm_33, \
                         pm_34, pn_30, pn_31, pn_32, pn_33, pn_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * pm_30[k]
                      + pn_30[k];

            t_31[k] = -ab_x[k] * pm_31[k]
                      + pn_31[k];

            t_32[k] = -ab_x[k] * pm_32[k]
                      + pn_32[k];

            t_33[k] = -ab_x[k] * pm_33[k]
                      + pn_33[k];

            t_34[k] = -ab_x[k] * pm_34[k]
                      + pn_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pm_35, pm_36, pm_37, pm_38, \
                         pm_39, pn_35, pn_36, pn_37, pn_38, pn_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * pm_35[k]
                      + pn_35[k];

            t_36[k] = -ab_x[k] * pm_36[k]
                      + pn_36[k];

            t_37[k] = -ab_x[k] * pm_37[k]
                      + pn_37[k];

            t_38[k] = -ab_x[k] * pm_38[k]
                      + pn_38[k];

            t_39[k] = -ab_x[k] * pm_39[k]
                      + pn_39[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pm_40, pm_41, pm_42, pm_43, \
                         pm_44, pn_40, pn_41, pn_42, pn_43, pn_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * pm_40[k]
                      + pn_40[k];

            t_41[k] = -ab_x[k] * pm_41[k]
                      + pn_41[k];

            t_42[k] = -ab_x[k] * pm_42[k]
                      + pn_42[k];

            t_43[k] = -ab_x[k] * pm_43[k]
                      + pn_43[k];

            t_44[k] = -ab_x[k] * pm_44[k]
                      + pn_44[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, pm_45, pm_46, pm_47, pm_48, \
                         pm_49, pn_45, pn_46, pn_47, pn_48, pn_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * pm_45[k]
                      + pn_45[k];

            t_46[k] = -ab_x[k] * pm_46[k]
                      + pn_46[k];

            t_47[k] = -ab_x[k] * pm_47[k]
                      + pn_47[k];

            t_48[k] = -ab_x[k] * pm_48[k]
                      + pn_48[k];

            t_49[k] = -ab_x[k] * pm_49[k]
                      + pn_49[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, pm_50, pm_51, pm_52, pm_53, \
                         pm_54, pn_50, pn_51, pn_52, pn_53, pn_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * pm_50[k]
                      + pn_50[k];

            t_51[k] = -ab_x[k] * pm_51[k]
                      + pn_51[k];

            t_52[k] = -ab_x[k] * pm_52[k]
                      + pn_52[k];

            t_53[k] = -ab_x[k] * pm_53[k]
                      + pn_53[k];

            t_54[k] = -ab_x[k] * pm_54[k]
                      + pn_54[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, pm_55, pm_56, pm_57, pm_58, \
                         pm_59, pn_66, pn_67, pn_68, pn_69, pn_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * pm_55[k]
                      + pn_66[k];

            t_56[k] = -ab_x[k] * pm_56[k]
                      + pn_67[k];

            t_57[k] = -ab_x[k] * pm_57[k]
                      + pn_68[k];

            t_58[k] = -ab_x[k] * pm_58[k]
                      + pn_69[k];

            t_59[k] = -ab_x[k] * pm_59[k]
                      + pn_70[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, pm_60, pm_61, pm_62, pm_63, \
                         pm_64, pn_71, pn_72, pn_73, pn_74, pn_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * pm_60[k]
                      + pn_71[k];

            t_61[k] = -ab_x[k] * pm_61[k]
                      + pn_72[k];

            t_62[k] = -ab_x[k] * pm_62[k]
                      + pn_73[k];

            t_63[k] = -ab_x[k] * pm_63[k]
                      + pn_74[k];

            t_64[k] = -ab_x[k] * pm_64[k]
                      + pn_75[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, pm_65, pm_66, pm_67, pm_68, \
                         pm_69, pn_76, pn_77, pn_78, pn_79, pn_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * pm_65[k]
                      + pn_76[k];

            t_66[k] = -ab_x[k] * pm_66[k]
                      + pn_77[k];

            t_67[k] = -ab_x[k] * pm_67[k]
                      + pn_78[k];

            t_68[k] = -ab_x[k] * pm_68[k]
                      + pn_79[k];

            t_69[k] = -ab_x[k] * pm_69[k]
                      + pn_80[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, pm_70, pm_71, pm_72, pm_73, \
                         pm_74, pn_81, pn_82, pn_83, pn_84, pn_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * pm_70[k]
                      + pn_81[k];

            t_71[k] = -ab_x[k] * pm_71[k]
                      + pn_82[k];

            t_72[k] = -ab_x[k] * pm_72[k]
                      + pn_83[k];

            t_73[k] = -ab_x[k] * pm_73[k]
                      + pn_84[k];

            t_74[k] = -ab_x[k] * pm_74[k]
                      + pn_85[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, pm_75, pm_76, pm_77, pm_78, \
                         pm_79, pn_86, pn_87, pn_88, pn_89, pn_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * pm_75[k]
                      + pn_86[k];

            t_76[k] = -ab_x[k] * pm_76[k]
                      + pn_87[k];

            t_77[k] = -ab_x[k] * pm_77[k]
                      + pn_88[k];

            t_78[k] = -ab_x[k] * pm_78[k]
                      + pn_89[k];

            t_79[k] = -ab_x[k] * pm_79[k]
                      + pn_90[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, pm_80, pm_81, pm_82, pm_83, \
                         pm_84, pn_91, pn_92, pn_93, pn_94, pn_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * pm_80[k]
                      + pn_91[k];

            t_81[k] = -ab_x[k] * pm_81[k]
                      + pn_92[k];

            t_82[k] = -ab_x[k] * pm_82[k]
                      + pn_93[k];

            t_83[k] = -ab_x[k] * pm_83[k]
                      + pn_94[k];

            t_84[k] = -ab_x[k] * pm_84[k]
                      + pn_95[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, pm_85, pm_86, pm_87, pm_88, \
                         pm_89, pn_96, pn_97, pn_98, pn_99, pn_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * pm_85[k]
                      + pn_96[k];

            t_86[k] = -ab_x[k] * pm_86[k]
                      + pn_97[k];

            t_87[k] = -ab_x[k] * pm_87[k]
                      + pn_98[k];

            t_88[k] = -ab_x[k] * pm_88[k]
                      + pn_99[k];

            t_89[k] = -ab_x[k] * pm_89[k]
                      + pn_100[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, pm_90, pm_91, pm_92, pm_93, \
                         pm_94, pn_101, pn_102, pn_103, pn_104, \
                         pn_105 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * pm_90[k]
                      + pn_101[k];

            t_91[k] = -ab_x[k] * pm_91[k]
                      + pn_102[k];

            t_92[k] = -ab_x[k] * pm_92[k]
                      + pn_103[k];

            t_93[k] = -ab_x[k] * pm_93[k]
                      + pn_104[k];

            t_94[k] = -ab_x[k] * pm_94[k]
                      + pn_105[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, pm_95, pm_96, pm_97, pm_98, \
                         pm_99, pn_106, pn_107, pn_108, pn_109, \
                         pn_110 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * pm_95[k]
                      + pn_106[k];

            t_96[k] = -ab_x[k] * pm_96[k]
                      + pn_107[k];

            t_97[k] = -ab_x[k] * pm_97[k]
                      + pn_108[k];

            t_98[k] = -ab_x[k] * pm_98[k]
                      + pn_109[k];

            t_99[k] = -ab_x[k] * pm_99[k]
                      + pn_110[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, pm_100, pm_101, pm_102, \
                         pm_103, pm_104, pn_111, pn_112, pn_113, pn_114, \
                         pn_115 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * pm_100[k]
                       + pn_111[k];

            t_101[k] = -ab_x[k] * pm_101[k]
                       + pn_112[k];

            t_102[k] = -ab_x[k] * pm_102[k]
                       + pn_113[k];

            t_103[k] = -ab_x[k] * pm_103[k]
                       + pn_114[k];

            t_104[k] = -ab_x[k] * pm_104[k]
                       + pn_115[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, pm_105, pm_106, pm_107, \
                         pm_108, pm_109, pn_116, pn_117, pn_118, pn_119, \
                         pn_120 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * pm_105[k]
                       + pn_116[k];

            t_106[k] = -ab_x[k] * pm_106[k]
                       + pn_117[k];

            t_107[k] = -ab_x[k] * pm_107[k]
                       + pn_118[k];

            t_108[k] = -ab_x[k] * pm_108[k]
                       + pn_119[k];

            t_109[k] = -ab_x[k] * pm_109[k]
                       + pn_120[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, pm_110, pm_111, pm_112, \
                         pm_113, pm_114, pn_132, pn_133, pn_134, pn_135, \
                         pn_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * pm_110[k]
                       + pn_132[k];

            t_111[k] = -ab_x[k] * pm_111[k]
                       + pn_133[k];

            t_112[k] = -ab_x[k] * pm_112[k]
                       + pn_134[k];

            t_113[k] = -ab_x[k] * pm_113[k]
                       + pn_135[k];

            t_114[k] = -ab_x[k] * pm_114[k]
                       + pn_136[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, pm_115, pm_116, pm_117, \
                         pm_118, pm_119, pn_137, pn_138, pn_139, pn_140, \
                         pn_141 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * pm_115[k]
                       + pn_137[k];

            t_116[k] = -ab_x[k] * pm_116[k]
                       + pn_138[k];

            t_117[k] = -ab_x[k] * pm_117[k]
                       + pn_139[k];

            t_118[k] = -ab_x[k] * pm_118[k]
                       + pn_140[k];

            t_119[k] = -ab_x[k] * pm_119[k]
                       + pn_141[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, pm_120, pm_121, pm_122, \
                         pm_123, pm_124, pn_142, pn_143, pn_144, pn_145, \
                         pn_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * pm_120[k]
                       + pn_142[k];

            t_121[k] = -ab_x[k] * pm_121[k]
                       + pn_143[k];

            t_122[k] = -ab_x[k] * pm_122[k]
                       + pn_144[k];

            t_123[k] = -ab_x[k] * pm_123[k]
                       + pn_145[k];

            t_124[k] = -ab_x[k] * pm_124[k]
                       + pn_146[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, pm_125, pm_126, pm_127, \
                         pm_128, pm_129, pn_147, pn_148, pn_149, pn_150, \
                         pn_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * pm_125[k]
                       + pn_147[k];

            t_126[k] = -ab_x[k] * pm_126[k]
                       + pn_148[k];

            t_127[k] = -ab_x[k] * pm_127[k]
                       + pn_149[k];

            t_128[k] = -ab_x[k] * pm_128[k]
                       + pn_150[k];

            t_129[k] = -ab_x[k] * pm_129[k]
                       + pn_151[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, pm_130, pm_131, pm_132, \
                         pm_133, pm_134, pn_152, pn_153, pn_154, pn_155, \
                         pn_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * pm_130[k]
                       + pn_152[k];

            t_131[k] = -ab_x[k] * pm_131[k]
                       + pn_153[k];

            t_132[k] = -ab_x[k] * pm_132[k]
                       + pn_154[k];

            t_133[k] = -ab_x[k] * pm_133[k]
                       + pn_155[k];

            t_134[k] = -ab_x[k] * pm_134[k]
                       + pn_156[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, pm_135, pm_136, pm_137, \
                         pm_138, pm_139, pn_157, pn_158, pn_159, pn_160, \
                         pn_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * pm_135[k]
                       + pn_157[k];

            t_136[k] = -ab_x[k] * pm_136[k]
                       + pn_158[k];

            t_137[k] = -ab_x[k] * pm_137[k]
                       + pn_159[k];

            t_138[k] = -ab_x[k] * pm_138[k]
                       + pn_160[k];

            t_139[k] = -ab_x[k] * pm_139[k]
                       + pn_161[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, pm_140, pm_141, pm_142, \
                         pm_143, pm_144, pn_162, pn_163, pn_164, pn_165, \
                         pn_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * pm_140[k]
                       + pn_162[k];

            t_141[k] = -ab_x[k] * pm_141[k]
                       + pn_163[k];

            t_142[k] = -ab_x[k] * pm_142[k]
                       + pn_164[k];

            t_143[k] = -ab_x[k] * pm_143[k]
                       + pn_165[k];

            t_144[k] = -ab_x[k] * pm_144[k]
                       + pn_166[k];
        }
    }
}

static auto
compute_hrr_dm_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pm, const size_t pn, const size_t ncomps,
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
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pm_55 = buffer.data(pm + 55 * ncomps + c);
        const auto *pm_56 = buffer.data(pm + 56 * ncomps + c);
        const auto *pm_57 = buffer.data(pm + 57 * ncomps + c);
        const auto *pm_58 = buffer.data(pm + 58 * ncomps + c);
        const auto *pm_59 = buffer.data(pm + 59 * ncomps + c);
        const auto *pm_60 = buffer.data(pm + 60 * ncomps + c);
        const auto *pm_61 = buffer.data(pm + 61 * ncomps + c);
        const auto *pm_62 = buffer.data(pm + 62 * ncomps + c);
        const auto *pm_63 = buffer.data(pm + 63 * ncomps + c);
        const auto *pm_64 = buffer.data(pm + 64 * ncomps + c);
        const auto *pm_65 = buffer.data(pm + 65 * ncomps + c);
        const auto *pm_66 = buffer.data(pm + 66 * ncomps + c);
        const auto *pm_67 = buffer.data(pm + 67 * ncomps + c);
        const auto *pm_68 = buffer.data(pm + 68 * ncomps + c);
        const auto *pm_69 = buffer.data(pm + 69 * ncomps + c);
        const auto *pm_70 = buffer.data(pm + 70 * ncomps + c);
        const auto *pm_71 = buffer.data(pm + 71 * ncomps + c);
        const auto *pm_72 = buffer.data(pm + 72 * ncomps + c);
        const auto *pm_73 = buffer.data(pm + 73 * ncomps + c);
        const auto *pm_74 = buffer.data(pm + 74 * ncomps + c);
        const auto *pm_75 = buffer.data(pm + 75 * ncomps + c);
        const auto *pm_76 = buffer.data(pm + 76 * ncomps + c);
        const auto *pm_77 = buffer.data(pm + 77 * ncomps + c);
        const auto *pm_78 = buffer.data(pm + 78 * ncomps + c);
        const auto *pm_79 = buffer.data(pm + 79 * ncomps + c);
        const auto *pm_80 = buffer.data(pm + 80 * ncomps + c);
        const auto *pm_81 = buffer.data(pm + 81 * ncomps + c);
        const auto *pm_82 = buffer.data(pm + 82 * ncomps + c);
        const auto *pm_83 = buffer.data(pm + 83 * ncomps + c);
        const auto *pm_84 = buffer.data(pm + 84 * ncomps + c);
        const auto *pm_85 = buffer.data(pm + 85 * ncomps + c);
        const auto *pm_86 = buffer.data(pm + 86 * ncomps + c);
        const auto *pm_87 = buffer.data(pm + 87 * ncomps + c);
        const auto *pm_88 = buffer.data(pm + 88 * ncomps + c);
        const auto *pm_89 = buffer.data(pm + 89 * ncomps + c);
        const auto *pm_90 = buffer.data(pm + 90 * ncomps + c);
        const auto *pm_91 = buffer.data(pm + 91 * ncomps + c);
        const auto *pm_92 = buffer.data(pm + 92 * ncomps + c);
        const auto *pm_93 = buffer.data(pm + 93 * ncomps + c);
        const auto *pm_94 = buffer.data(pm + 94 * ncomps + c);
        const auto *pm_95 = buffer.data(pm + 95 * ncomps + c);
        const auto *pm_96 = buffer.data(pm + 96 * ncomps + c);
        const auto *pm_97 = buffer.data(pm + 97 * ncomps + c);
        const auto *pm_98 = buffer.data(pm + 98 * ncomps + c);
        const auto *pm_99 = buffer.data(pm + 99 * ncomps + c);
        const auto *pm_100 = buffer.data(pm + 100 * ncomps + c);
        const auto *pm_101 = buffer.data(pm + 101 * ncomps + c);
        const auto *pm_102 = buffer.data(pm + 102 * ncomps + c);
        const auto *pm_103 = buffer.data(pm + 103 * ncomps + c);
        const auto *pm_104 = buffer.data(pm + 104 * ncomps + c);
        const auto *pm_105 = buffer.data(pm + 105 * ncomps + c);
        const auto *pm_106 = buffer.data(pm + 106 * ncomps + c);
        const auto *pm_107 = buffer.data(pm + 107 * ncomps + c);
        const auto *pm_108 = buffer.data(pm + 108 * ncomps + c);
        const auto *pm_109 = buffer.data(pm + 109 * ncomps + c);
        const auto *pm_110 = buffer.data(pm + 110 * ncomps + c);
        const auto *pm_111 = buffer.data(pm + 111 * ncomps + c);
        const auto *pm_112 = buffer.data(pm + 112 * ncomps + c);
        const auto *pm_113 = buffer.data(pm + 113 * ncomps + c);
        const auto *pm_114 = buffer.data(pm + 114 * ncomps + c);
        const auto *pm_115 = buffer.data(pm + 115 * ncomps + c);
        const auto *pm_116 = buffer.data(pm + 116 * ncomps + c);
        const auto *pm_117 = buffer.data(pm + 117 * ncomps + c);
        const auto *pm_118 = buffer.data(pm + 118 * ncomps + c);
        const auto *pm_119 = buffer.data(pm + 119 * ncomps + c);
        const auto *pm_120 = buffer.data(pm + 120 * ncomps + c);
        const auto *pm_121 = buffer.data(pm + 121 * ncomps + c);
        const auto *pm_122 = buffer.data(pm + 122 * ncomps + c);
        const auto *pm_123 = buffer.data(pm + 123 * ncomps + c);
        const auto *pm_124 = buffer.data(pm + 124 * ncomps + c);
        const auto *pm_125 = buffer.data(pm + 125 * ncomps + c);
        const auto *pm_126 = buffer.data(pm + 126 * ncomps + c);
        const auto *pm_127 = buffer.data(pm + 127 * ncomps + c);
        const auto *pm_128 = buffer.data(pm + 128 * ncomps + c);
        const auto *pm_129 = buffer.data(pm + 129 * ncomps + c);
        const auto *pm_130 = buffer.data(pm + 130 * ncomps + c);
        const auto *pm_131 = buffer.data(pm + 131 * ncomps + c);
        const auto *pm_132 = buffer.data(pm + 132 * ncomps + c);
        const auto *pm_133 = buffer.data(pm + 133 * ncomps + c);
        const auto *pm_134 = buffer.data(pm + 134 * ncomps + c);
        const auto *pm_135 = buffer.data(pm + 135 * ncomps + c);
        const auto *pm_136 = buffer.data(pm + 136 * ncomps + c);
        const auto *pm_137 = buffer.data(pm + 137 * ncomps + c);
        const auto *pm_138 = buffer.data(pm + 138 * ncomps + c);
        const auto *pm_139 = buffer.data(pm + 139 * ncomps + c);
        const auto *pm_140 = buffer.data(pm + 140 * ncomps + c);
        const auto *pm_141 = buffer.data(pm + 141 * ncomps + c);
        const auto *pm_142 = buffer.data(pm + 142 * ncomps + c);
        const auto *pm_143 = buffer.data(pm + 143 * ncomps + c);
        const auto *pm_144 = buffer.data(pm + 144 * ncomps + c);
        const auto *pm_145 = buffer.data(pm + 145 * ncomps + c);
        const auto *pm_146 = buffer.data(pm + 146 * ncomps + c);
        const auto *pm_147 = buffer.data(pm + 147 * ncomps + c);
        const auto *pm_148 = buffer.data(pm + 148 * ncomps + c);
        const auto *pm_149 = buffer.data(pm + 149 * ncomps + c);
        const auto *pm_150 = buffer.data(pm + 150 * ncomps + c);
        const auto *pm_151 = buffer.data(pm + 151 * ncomps + c);
        const auto *pm_152 = buffer.data(pm + 152 * ncomps + c);
        const auto *pm_153 = buffer.data(pm + 153 * ncomps + c);
        const auto *pm_154 = buffer.data(pm + 154 * ncomps + c);
        const auto *pm_155 = buffer.data(pm + 155 * ncomps + c);
        const auto *pm_156 = buffer.data(pm + 156 * ncomps + c);
        const auto *pm_157 = buffer.data(pm + 157 * ncomps + c);
        const auto *pm_158 = buffer.data(pm + 158 * ncomps + c);
        const auto *pm_159 = buffer.data(pm + 159 * ncomps + c);
        const auto *pm_160 = buffer.data(pm + 160 * ncomps + c);
        const auto *pm_161 = buffer.data(pm + 161 * ncomps + c);
        const auto *pm_162 = buffer.data(pm + 162 * ncomps + c);
        const auto *pm_163 = buffer.data(pm + 163 * ncomps + c);
        const auto *pm_164 = buffer.data(pm + 164 * ncomps + c);

        const auto *pn_67 = buffer.data(pn + 67 * ncomps + c);
        const auto *pn_69 = buffer.data(pn + 69 * ncomps + c);
        const auto *pn_70 = buffer.data(pn + 70 * ncomps + c);
        const auto *pn_72 = buffer.data(pn + 72 * ncomps + c);
        const auto *pn_73 = buffer.data(pn + 73 * ncomps + c);
        const auto *pn_74 = buffer.data(pn + 74 * ncomps + c);
        const auto *pn_76 = buffer.data(pn + 76 * ncomps + c);
        const auto *pn_77 = buffer.data(pn + 77 * ncomps + c);
        const auto *pn_78 = buffer.data(pn + 78 * ncomps + c);
        const auto *pn_79 = buffer.data(pn + 79 * ncomps + c);
        const auto *pn_81 = buffer.data(pn + 81 * ncomps + c);
        const auto *pn_82 = buffer.data(pn + 82 * ncomps + c);
        const auto *pn_83 = buffer.data(pn + 83 * ncomps + c);
        const auto *pn_84 = buffer.data(pn + 84 * ncomps + c);
        const auto *pn_85 = buffer.data(pn + 85 * ncomps + c);
        const auto *pn_87 = buffer.data(pn + 87 * ncomps + c);
        const auto *pn_88 = buffer.data(pn + 88 * ncomps + c);
        const auto *pn_89 = buffer.data(pn + 89 * ncomps + c);
        const auto *pn_90 = buffer.data(pn + 90 * ncomps + c);
        const auto *pn_91 = buffer.data(pn + 91 * ncomps + c);
        const auto *pn_92 = buffer.data(pn + 92 * ncomps + c);
        const auto *pn_94 = buffer.data(pn + 94 * ncomps + c);
        const auto *pn_95 = buffer.data(pn + 95 * ncomps + c);
        const auto *pn_96 = buffer.data(pn + 96 * ncomps + c);
        const auto *pn_97 = buffer.data(pn + 97 * ncomps + c);
        const auto *pn_98 = buffer.data(pn + 98 * ncomps + c);
        const auto *pn_99 = buffer.data(pn + 99 * ncomps + c);
        const auto *pn_100 = buffer.data(pn + 100 * ncomps + c);
        const auto *pn_102 = buffer.data(pn + 102 * ncomps + c);
        const auto *pn_103 = buffer.data(pn + 103 * ncomps + c);
        const auto *pn_104 = buffer.data(pn + 104 * ncomps + c);
        const auto *pn_105 = buffer.data(pn + 105 * ncomps + c);
        const auto *pn_106 = buffer.data(pn + 106 * ncomps + c);
        const auto *pn_107 = buffer.data(pn + 107 * ncomps + c);
        const auto *pn_108 = buffer.data(pn + 108 * ncomps + c);
        const auto *pn_109 = buffer.data(pn + 109 * ncomps + c);
        const auto *pn_111 = buffer.data(pn + 111 * ncomps + c);
        const auto *pn_112 = buffer.data(pn + 112 * ncomps + c);
        const auto *pn_113 = buffer.data(pn + 113 * ncomps + c);
        const auto *pn_114 = buffer.data(pn + 114 * ncomps + c);
        const auto *pn_115 = buffer.data(pn + 115 * ncomps + c);
        const auto *pn_116 = buffer.data(pn + 116 * ncomps + c);
        const auto *pn_117 = buffer.data(pn + 117 * ncomps + c);
        const auto *pn_118 = buffer.data(pn + 118 * ncomps + c);
        const auto *pn_119 = buffer.data(pn + 119 * ncomps + c);
        const auto *pn_121 = buffer.data(pn + 121 * ncomps + c);
        const auto *pn_122 = buffer.data(pn + 122 * ncomps + c);
        const auto *pn_123 = buffer.data(pn + 123 * ncomps + c);
        const auto *pn_124 = buffer.data(pn + 124 * ncomps + c);
        const auto *pn_125 = buffer.data(pn + 125 * ncomps + c);
        const auto *pn_126 = buffer.data(pn + 126 * ncomps + c);
        const auto *pn_127 = buffer.data(pn + 127 * ncomps + c);
        const auto *pn_128 = buffer.data(pn + 128 * ncomps + c);
        const auto *pn_129 = buffer.data(pn + 129 * ncomps + c);
        const auto *pn_130 = buffer.data(pn + 130 * ncomps + c);
        const auto *pn_133 = buffer.data(pn + 133 * ncomps + c);
        const auto *pn_134 = buffer.data(pn + 134 * ncomps + c);
        const auto *pn_135 = buffer.data(pn + 135 * ncomps + c);
        const auto *pn_136 = buffer.data(pn + 136 * ncomps + c);
        const auto *pn_137 = buffer.data(pn + 137 * ncomps + c);
        const auto *pn_138 = buffer.data(pn + 138 * ncomps + c);
        const auto *pn_139 = buffer.data(pn + 139 * ncomps + c);
        const auto *pn_140 = buffer.data(pn + 140 * ncomps + c);
        const auto *pn_141 = buffer.data(pn + 141 * ncomps + c);
        const auto *pn_142 = buffer.data(pn + 142 * ncomps + c);
        const auto *pn_143 = buffer.data(pn + 143 * ncomps + c);
        const auto *pn_144 = buffer.data(pn + 144 * ncomps + c);
        const auto *pn_145 = buffer.data(pn + 145 * ncomps + c);
        const auto *pn_146 = buffer.data(pn + 146 * ncomps + c);
        const auto *pn_147 = buffer.data(pn + 147 * ncomps + c);
        const auto *pn_148 = buffer.data(pn + 148 * ncomps + c);
        const auto *pn_149 = buffer.data(pn + 149 * ncomps + c);
        const auto *pn_150 = buffer.data(pn + 150 * ncomps + c);
        const auto *pn_151 = buffer.data(pn + 151 * ncomps + c);
        const auto *pn_152 = buffer.data(pn + 152 * ncomps + c);
        const auto *pn_153 = buffer.data(pn + 153 * ncomps + c);
        const auto *pn_154 = buffer.data(pn + 154 * ncomps + c);
        const auto *pn_155 = buffer.data(pn + 155 * ncomps + c);
        const auto *pn_156 = buffer.data(pn + 156 * ncomps + c);
        const auto *pn_157 = buffer.data(pn + 157 * ncomps + c);
        const auto *pn_158 = buffer.data(pn + 158 * ncomps + c);
        const auto *pn_160 = buffer.data(pn + 160 * ncomps + c);
        const auto *pn_161 = buffer.data(pn + 161 * ncomps + c);
        const auto *pn_162 = buffer.data(pn + 162 * ncomps + c);
        const auto *pn_163 = buffer.data(pn + 163 * ncomps + c);
        const auto *pn_164 = buffer.data(pn + 164 * ncomps + c);
        const auto *pn_165 = buffer.data(pn + 165 * ncomps + c);
        const auto *pn_166 = buffer.data(pn + 166 * ncomps + c);
        const auto *pn_167 = buffer.data(pn + 167 * ncomps + c);
        const auto *pn_168 = buffer.data(pn + 168 * ncomps + c);
        const auto *pn_169 = buffer.data(pn + 169 * ncomps + c);
        const auto *pn_170 = buffer.data(pn + 170 * ncomps + c);
        const auto *pn_171 = buffer.data(pn + 171 * ncomps + c);
        const auto *pn_172 = buffer.data(pn + 172 * ncomps + c);
        const auto *pn_173 = buffer.data(pn + 173 * ncomps + c);
        const auto *pn_174 = buffer.data(pn + 174 * ncomps + c);
        const auto *pn_175 = buffer.data(pn + 175 * ncomps + c);
        const auto *pn_176 = buffer.data(pn + 176 * ncomps + c);
        const auto *pn_177 = buffer.data(pn + 177 * ncomps + c);
        const auto *pn_178 = buffer.data(pn + 178 * ncomps + c);
        const auto *pn_179 = buffer.data(pn + 179 * ncomps + c);
        const auto *pn_180 = buffer.data(pn + 180 * ncomps + c);
        const auto *pn_181 = buffer.data(pn + 181 * ncomps + c);
        const auto *pn_182 = buffer.data(pn + 182 * ncomps + c);
        const auto *pn_183 = buffer.data(pn + 183 * ncomps + c);
        const auto *pn_184 = buffer.data(pn + 184 * ncomps + c);
        const auto *pn_185 = buffer.data(pn + 185 * ncomps + c);
        const auto *pn_186 = buffer.data(pn + 186 * ncomps + c);
        const auto *pn_187 = buffer.data(pn + 187 * ncomps + c);
        const auto *pn_188 = buffer.data(pn + 188 * ncomps + c);
        const auto *pn_189 = buffer.data(pn + 189 * ncomps + c);
        const auto *pn_190 = buffer.data(pn + 190 * ncomps + c);
        const auto *pn_191 = buffer.data(pn + 191 * ncomps + c);
        const auto *pn_192 = buffer.data(pn + 192 * ncomps + c);
        const auto *pn_193 = buffer.data(pn + 193 * ncomps + c);
        const auto *pn_194 = buffer.data(pn + 194 * ncomps + c);
        const auto *pn_195 = buffer.data(pn + 195 * ncomps + c);
        const auto *pn_196 = buffer.data(pn + 196 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, pm_145, pm_146, pm_147, \
                         pm_148, pm_149, pn_167, pn_168, pn_169, pn_170, \
                         pn_171 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * pm_145[k]
                       + pn_167[k];

            t_146[k] = -ab_x[k] * pm_146[k]
                       + pn_168[k];

            t_147[k] = -ab_x[k] * pm_147[k]
                       + pn_169[k];

            t_148[k] = -ab_x[k] * pm_148[k]
                       + pn_170[k];

            t_149[k] = -ab_x[k] * pm_149[k]
                       + pn_171[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, pm_150, pm_151, pm_152, \
                         pm_153, pm_154, pn_172, pn_173, pn_174, pn_175, \
                         pn_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * pm_150[k]
                       + pn_172[k];

            t_151[k] = -ab_x[k] * pm_151[k]
                       + pn_173[k];

            t_152[k] = -ab_x[k] * pm_152[k]
                       + pn_174[k];

            t_153[k] = -ab_x[k] * pm_153[k]
                       + pn_175[k];

            t_154[k] = -ab_x[k] * pm_154[k]
                       + pn_176[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, pm_155, pm_156, pm_157, \
                         pm_158, pm_159, pn_177, pn_178, pn_179, pn_180, \
                         pn_181 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * pm_155[k]
                       + pn_177[k];

            t_156[k] = -ab_x[k] * pm_156[k]
                       + pn_178[k];

            t_157[k] = -ab_x[k] * pm_157[k]
                       + pn_179[k];

            t_158[k] = -ab_x[k] * pm_158[k]
                       + pn_180[k];

            t_159[k] = -ab_x[k] * pm_159[k]
                       + pn_181[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, pm_160, pm_161, pm_162, \
                         pm_163, pm_164, pn_182, pn_183, pn_184, pn_185, \
                         pn_186 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * pm_160[k]
                       + pn_182[k];

            t_161[k] = -ab_x[k] * pm_161[k]
                       + pn_183[k];

            t_162[k] = -ab_x[k] * pm_162[k]
                       + pn_184[k];

            t_163[k] = -ab_x[k] * pm_163[k]
                       + pn_185[k];

            t_164[k] = -ab_x[k] * pm_164[k]
                       + pn_186[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_y, pm_55, pm_56, pm_57, pm_58, \
                         pm_59, pn_67, pn_69, pn_70, pn_72, pn_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_y[k] * pm_55[k]
                       + pn_67[k];

            t_166[k] = -ab_y[k] * pm_56[k]
                       + pn_69[k];

            t_167[k] = -ab_y[k] * pm_57[k]
                       + pn_70[k];

            t_168[k] = -ab_y[k] * pm_58[k]
                       + pn_72[k];

            t_169[k] = -ab_y[k] * pm_59[k]
                       + pn_73[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_y, pm_60, pm_61, pm_62, pm_63, \
                         pm_64, pn_74, pn_76, pn_77, pn_78, pn_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_y[k] * pm_60[k]
                       + pn_74[k];

            t_171[k] = -ab_y[k] * pm_61[k]
                       + pn_76[k];

            t_172[k] = -ab_y[k] * pm_62[k]
                       + pn_77[k];

            t_173[k] = -ab_y[k] * pm_63[k]
                       + pn_78[k];

            t_174[k] = -ab_y[k] * pm_64[k]
                       + pn_79[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, pm_65, pm_66, pm_67, pm_68, \
                         pm_69, pn_81, pn_82, pn_83, pn_84, pn_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_y[k] * pm_65[k]
                       + pn_81[k];

            t_176[k] = -ab_y[k] * pm_66[k]
                       + pn_82[k];

            t_177[k] = -ab_y[k] * pm_67[k]
                       + pn_83[k];

            t_178[k] = -ab_y[k] * pm_68[k]
                       + pn_84[k];

            t_179[k] = -ab_y[k] * pm_69[k]
                       + pn_85[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_y, pm_70, pm_71, pm_72, pm_73, \
                         pm_74, pn_87, pn_88, pn_89, pn_90, pn_91 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_y[k] * pm_70[k]
                       + pn_87[k];

            t_181[k] = -ab_y[k] * pm_71[k]
                       + pn_88[k];

            t_182[k] = -ab_y[k] * pm_72[k]
                       + pn_89[k];

            t_183[k] = -ab_y[k] * pm_73[k]
                       + pn_90[k];

            t_184[k] = -ab_y[k] * pm_74[k]
                       + pn_91[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_y, pm_75, pm_76, pm_77, pm_78, \
                         pm_79, pn_92, pn_94, pn_95, pn_96, pn_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_y[k] * pm_75[k]
                       + pn_92[k];

            t_186[k] = -ab_y[k] * pm_76[k]
                       + pn_94[k];

            t_187[k] = -ab_y[k] * pm_77[k]
                       + pn_95[k];

            t_188[k] = -ab_y[k] * pm_78[k]
                       + pn_96[k];

            t_189[k] = -ab_y[k] * pm_79[k]
                       + pn_97[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, pm_80, pm_81, pm_82, pm_83, \
                         pm_84, pn_98, pn_99, pn_100, pn_102, pn_103 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_y[k] * pm_80[k]
                       + pn_98[k];

            t_191[k] = -ab_y[k] * pm_81[k]
                       + pn_99[k];

            t_192[k] = -ab_y[k] * pm_82[k]
                       + pn_100[k];

            t_193[k] = -ab_y[k] * pm_83[k]
                       + pn_102[k];

            t_194[k] = -ab_y[k] * pm_84[k]
                       + pn_103[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_y, pm_85, pm_86, pm_87, pm_88, \
                         pm_89, pn_104, pn_105, pn_106, pn_107, \
                         pn_108 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_y[k] * pm_85[k]
                       + pn_104[k];

            t_196[k] = -ab_y[k] * pm_86[k]
                       + pn_105[k];

            t_197[k] = -ab_y[k] * pm_87[k]
                       + pn_106[k];

            t_198[k] = -ab_y[k] * pm_88[k]
                       + pn_107[k];

            t_199[k] = -ab_y[k] * pm_89[k]
                       + pn_108[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_y, pm_90, pm_91, pm_92, pm_93, \
                         pm_94, pn_109, pn_111, pn_112, pn_113, \
                         pn_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_y[k] * pm_90[k]
                       + pn_109[k];

            t_201[k] = -ab_y[k] * pm_91[k]
                       + pn_111[k];

            t_202[k] = -ab_y[k] * pm_92[k]
                       + pn_112[k];

            t_203[k] = -ab_y[k] * pm_93[k]
                       + pn_113[k];

            t_204[k] = -ab_y[k] * pm_94[k]
                       + pn_114[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, pm_95, pm_96, pm_97, pm_98, \
                         pm_99, pn_115, pn_116, pn_117, pn_118, \
                         pn_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_y[k] * pm_95[k]
                       + pn_115[k];

            t_206[k] = -ab_y[k] * pm_96[k]
                       + pn_116[k];

            t_207[k] = -ab_y[k] * pm_97[k]
                       + pn_117[k];

            t_208[k] = -ab_y[k] * pm_98[k]
                       + pn_118[k];

            t_209[k] = -ab_y[k] * pm_99[k]
                       + pn_119[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_y, pm_100, pm_101, pm_102, \
                         pm_103, pm_104, pn_121, pn_122, pn_123, pn_124, \
                         pn_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_y[k] * pm_100[k]
                       + pn_121[k];

            t_211[k] = -ab_y[k] * pm_101[k]
                       + pn_122[k];

            t_212[k] = -ab_y[k] * pm_102[k]
                       + pn_123[k];

            t_213[k] = -ab_y[k] * pm_103[k]
                       + pn_124[k];

            t_214[k] = -ab_y[k] * pm_104[k]
                       + pn_125[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_y, pm_105, pm_106, pm_107, \
                         pm_108, pm_109, pn_126, pn_127, pn_128, pn_129, \
                         pn_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_y[k] * pm_105[k]
                       + pn_126[k];

            t_216[k] = -ab_y[k] * pm_106[k]
                       + pn_127[k];

            t_217[k] = -ab_y[k] * pm_107[k]
                       + pn_128[k];

            t_218[k] = -ab_y[k] * pm_108[k]
                       + pn_129[k];

            t_219[k] = -ab_y[k] * pm_109[k]
                       + pn_130[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, pm_110, pm_111, pm_112, \
                         pm_113, pm_114, pn_133, pn_135, pn_136, pn_138, \
                         pn_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_y[k] * pm_110[k]
                       + pn_133[k];

            t_221[k] = -ab_y[k] * pm_111[k]
                       + pn_135[k];

            t_222[k] = -ab_y[k] * pm_112[k]
                       + pn_136[k];

            t_223[k] = -ab_y[k] * pm_113[k]
                       + pn_138[k];

            t_224[k] = -ab_y[k] * pm_114[k]
                       + pn_139[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_y, pm_115, pm_116, pm_117, \
                         pm_118, pm_119, pn_140, pn_142, pn_143, pn_144, \
                         pn_145 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_y[k] * pm_115[k]
                       + pn_140[k];

            t_226[k] = -ab_y[k] * pm_116[k]
                       + pn_142[k];

            t_227[k] = -ab_y[k] * pm_117[k]
                       + pn_143[k];

            t_228[k] = -ab_y[k] * pm_118[k]
                       + pn_144[k];

            t_229[k] = -ab_y[k] * pm_119[k]
                       + pn_145[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_y, pm_120, pm_121, pm_122, \
                         pm_123, pm_124, pn_147, pn_148, pn_149, pn_150, \
                         pn_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_y[k] * pm_120[k]
                       + pn_147[k];

            t_231[k] = -ab_y[k] * pm_121[k]
                       + pn_148[k];

            t_232[k] = -ab_y[k] * pm_122[k]
                       + pn_149[k];

            t_233[k] = -ab_y[k] * pm_123[k]
                       + pn_150[k];

            t_234[k] = -ab_y[k] * pm_124[k]
                       + pn_151[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, pm_125, pm_126, pm_127, \
                         pm_128, pm_129, pn_153, pn_154, pn_155, pn_156, \
                         pn_157 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_y[k] * pm_125[k]
                       + pn_153[k];

            t_236[k] = -ab_y[k] * pm_126[k]
                       + pn_154[k];

            t_237[k] = -ab_y[k] * pm_127[k]
                       + pn_155[k];

            t_238[k] = -ab_y[k] * pm_128[k]
                       + pn_156[k];

            t_239[k] = -ab_y[k] * pm_129[k]
                       + pn_157[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_y, pm_130, pm_131, pm_132, \
                         pm_133, pm_134, pn_158, pn_160, pn_161, pn_162, \
                         pn_163 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_y[k] * pm_130[k]
                       + pn_158[k];

            t_241[k] = -ab_y[k] * pm_131[k]
                       + pn_160[k];

            t_242[k] = -ab_y[k] * pm_132[k]
                       + pn_161[k];

            t_243[k] = -ab_y[k] * pm_133[k]
                       + pn_162[k];

            t_244[k] = -ab_y[k] * pm_134[k]
                       + pn_163[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_y, pm_135, pm_136, pm_137, \
                         pm_138, pm_139, pn_164, pn_165, pn_166, pn_168, \
                         pn_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_y[k] * pm_135[k]
                       + pn_164[k];

            t_246[k] = -ab_y[k] * pm_136[k]
                       + pn_165[k];

            t_247[k] = -ab_y[k] * pm_137[k]
                       + pn_166[k];

            t_248[k] = -ab_y[k] * pm_138[k]
                       + pn_168[k];

            t_249[k] = -ab_y[k] * pm_139[k]
                       + pn_169[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, pm_140, pm_141, pm_142, \
                         pm_143, pm_144, pn_170, pn_171, pn_172, pn_173, \
                         pn_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_y[k] * pm_140[k]
                       + pn_170[k];

            t_251[k] = -ab_y[k] * pm_141[k]
                       + pn_171[k];

            t_252[k] = -ab_y[k] * pm_142[k]
                       + pn_172[k];

            t_253[k] = -ab_y[k] * pm_143[k]
                       + pn_173[k];

            t_254[k] = -ab_y[k] * pm_144[k]
                       + pn_174[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_y, pm_145, pm_146, pm_147, \
                         pm_148, pm_149, pn_175, pn_177, pn_178, pn_179, \
                         pn_180 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_y[k] * pm_145[k]
                       + pn_175[k];

            t_256[k] = -ab_y[k] * pm_146[k]
                       + pn_177[k];

            t_257[k] = -ab_y[k] * pm_147[k]
                       + pn_178[k];

            t_258[k] = -ab_y[k] * pm_148[k]
                       + pn_179[k];

            t_259[k] = -ab_y[k] * pm_149[k]
                       + pn_180[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_y, pm_150, pm_151, pm_152, \
                         pm_153, pm_154, pn_181, pn_182, pn_183, pn_184, \
                         pn_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_y[k] * pm_150[k]
                       + pn_181[k];

            t_261[k] = -ab_y[k] * pm_151[k]
                       + pn_182[k];

            t_262[k] = -ab_y[k] * pm_152[k]
                       + pn_183[k];

            t_263[k] = -ab_y[k] * pm_153[k]
                       + pn_184[k];

            t_264[k] = -ab_y[k] * pm_154[k]
                       + pn_185[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, pm_155, pm_156, pm_157, \
                         pm_158, pm_159, pn_187, pn_188, pn_189, pn_190, \
                         pn_191 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_y[k] * pm_155[k]
                       + pn_187[k];

            t_266[k] = -ab_y[k] * pm_156[k]
                       + pn_188[k];

            t_267[k] = -ab_y[k] * pm_157[k]
                       + pn_189[k];

            t_268[k] = -ab_y[k] * pm_158[k]
                       + pn_190[k];

            t_269[k] = -ab_y[k] * pm_159[k]
                       + pn_191[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_y, pm_160, pm_161, pm_162, \
                         pm_163, pm_164, pn_192, pn_193, pn_194, pn_195, \
                         pn_196 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_y[k] * pm_160[k]
                       + pn_192[k];

            t_271[k] = -ab_y[k] * pm_161[k]
                       + pn_193[k];

            t_272[k] = -ab_y[k] * pm_162[k]
                       + pn_194[k];

            t_273[k] = -ab_y[k] * pm_163[k]
                       + pn_195[k];

            t_274[k] = -ab_y[k] * pm_164[k]
                       + pn_196[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_z, pm_110, pm_111, pm_112, \
                         pm_113, pm_114, pn_134, pn_136, pn_137, pn_139, \
                         pn_140 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_z[k] * pm_110[k]
                       + pn_134[k];

            t_276[k] = -ab_z[k] * pm_111[k]
                       + pn_136[k];

            t_277[k] = -ab_z[k] * pm_112[k]
                       + pn_137[k];

            t_278[k] = -ab_z[k] * pm_113[k]
                       + pn_139[k];

            t_279[k] = -ab_z[k] * pm_114[k]
                       + pn_140[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_z, pm_115, pm_116, pm_117, \
                         pm_118, pm_119, pn_141, pn_143, pn_144, pn_145, \
                         pn_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_z[k] * pm_115[k]
                       + pn_141[k];

            t_281[k] = -ab_z[k] * pm_116[k]
                       + pn_143[k];

            t_282[k] = -ab_z[k] * pm_117[k]
                       + pn_144[k];

            t_283[k] = -ab_z[k] * pm_118[k]
                       + pn_145[k];

            t_284[k] = -ab_z[k] * pm_119[k]
                       + pn_146[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_z, pm_120, pm_121, pm_122, \
                         pm_123, pm_124, pn_148, pn_149, pn_150, pn_151, \
                         pn_152 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_z[k] * pm_120[k]
                       + pn_148[k];

            t_286[k] = -ab_z[k] * pm_121[k]
                       + pn_149[k];

            t_287[k] = -ab_z[k] * pm_122[k]
                       + pn_150[k];

            t_288[k] = -ab_z[k] * pm_123[k]
                       + pn_151[k];

            t_289[k] = -ab_z[k] * pm_124[k]
                       + pn_152[k];
        }
    }
}

static auto
compute_hrr_dm_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pm, const size_t pn, const size_t ncomps,
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

        const auto *ab_z = coordinates.data(8);

        const auto *pm_125 = buffer.data(pm + 125 * ncomps + c);
        const auto *pm_126 = buffer.data(pm + 126 * ncomps + c);
        const auto *pm_127 = buffer.data(pm + 127 * ncomps + c);
        const auto *pm_128 = buffer.data(pm + 128 * ncomps + c);
        const auto *pm_129 = buffer.data(pm + 129 * ncomps + c);
        const auto *pm_130 = buffer.data(pm + 130 * ncomps + c);
        const auto *pm_131 = buffer.data(pm + 131 * ncomps + c);
        const auto *pm_132 = buffer.data(pm + 132 * ncomps + c);
        const auto *pm_133 = buffer.data(pm + 133 * ncomps + c);
        const auto *pm_134 = buffer.data(pm + 134 * ncomps + c);
        const auto *pm_135 = buffer.data(pm + 135 * ncomps + c);
        const auto *pm_136 = buffer.data(pm + 136 * ncomps + c);
        const auto *pm_137 = buffer.data(pm + 137 * ncomps + c);
        const auto *pm_138 = buffer.data(pm + 138 * ncomps + c);
        const auto *pm_139 = buffer.data(pm + 139 * ncomps + c);
        const auto *pm_140 = buffer.data(pm + 140 * ncomps + c);
        const auto *pm_141 = buffer.data(pm + 141 * ncomps + c);
        const auto *pm_142 = buffer.data(pm + 142 * ncomps + c);
        const auto *pm_143 = buffer.data(pm + 143 * ncomps + c);
        const auto *pm_144 = buffer.data(pm + 144 * ncomps + c);
        const auto *pm_145 = buffer.data(pm + 145 * ncomps + c);
        const auto *pm_146 = buffer.data(pm + 146 * ncomps + c);
        const auto *pm_147 = buffer.data(pm + 147 * ncomps + c);
        const auto *pm_148 = buffer.data(pm + 148 * ncomps + c);
        const auto *pm_149 = buffer.data(pm + 149 * ncomps + c);
        const auto *pm_150 = buffer.data(pm + 150 * ncomps + c);
        const auto *pm_151 = buffer.data(pm + 151 * ncomps + c);
        const auto *pm_152 = buffer.data(pm + 152 * ncomps + c);
        const auto *pm_153 = buffer.data(pm + 153 * ncomps + c);
        const auto *pm_154 = buffer.data(pm + 154 * ncomps + c);
        const auto *pm_155 = buffer.data(pm + 155 * ncomps + c);
        const auto *pm_156 = buffer.data(pm + 156 * ncomps + c);
        const auto *pm_157 = buffer.data(pm + 157 * ncomps + c);
        const auto *pm_158 = buffer.data(pm + 158 * ncomps + c);
        const auto *pm_159 = buffer.data(pm + 159 * ncomps + c);
        const auto *pm_160 = buffer.data(pm + 160 * ncomps + c);
        const auto *pm_161 = buffer.data(pm + 161 * ncomps + c);
        const auto *pm_162 = buffer.data(pm + 162 * ncomps + c);
        const auto *pm_163 = buffer.data(pm + 163 * ncomps + c);
        const auto *pm_164 = buffer.data(pm + 164 * ncomps + c);

        const auto *pn_154 = buffer.data(pn + 154 * ncomps + c);
        const auto *pn_155 = buffer.data(pn + 155 * ncomps + c);
        const auto *pn_156 = buffer.data(pn + 156 * ncomps + c);
        const auto *pn_157 = buffer.data(pn + 157 * ncomps + c);
        const auto *pn_158 = buffer.data(pn + 158 * ncomps + c);
        const auto *pn_159 = buffer.data(pn + 159 * ncomps + c);
        const auto *pn_161 = buffer.data(pn + 161 * ncomps + c);
        const auto *pn_162 = buffer.data(pn + 162 * ncomps + c);
        const auto *pn_163 = buffer.data(pn + 163 * ncomps + c);
        const auto *pn_164 = buffer.data(pn + 164 * ncomps + c);
        const auto *pn_165 = buffer.data(pn + 165 * ncomps + c);
        const auto *pn_166 = buffer.data(pn + 166 * ncomps + c);
        const auto *pn_167 = buffer.data(pn + 167 * ncomps + c);
        const auto *pn_169 = buffer.data(pn + 169 * ncomps + c);
        const auto *pn_170 = buffer.data(pn + 170 * ncomps + c);
        const auto *pn_171 = buffer.data(pn + 171 * ncomps + c);
        const auto *pn_172 = buffer.data(pn + 172 * ncomps + c);
        const auto *pn_173 = buffer.data(pn + 173 * ncomps + c);
        const auto *pn_174 = buffer.data(pn + 174 * ncomps + c);
        const auto *pn_175 = buffer.data(pn + 175 * ncomps + c);
        const auto *pn_176 = buffer.data(pn + 176 * ncomps + c);
        const auto *pn_178 = buffer.data(pn + 178 * ncomps + c);
        const auto *pn_179 = buffer.data(pn + 179 * ncomps + c);
        const auto *pn_180 = buffer.data(pn + 180 * ncomps + c);
        const auto *pn_181 = buffer.data(pn + 181 * ncomps + c);
        const auto *pn_182 = buffer.data(pn + 182 * ncomps + c);
        const auto *pn_183 = buffer.data(pn + 183 * ncomps + c);
        const auto *pn_184 = buffer.data(pn + 184 * ncomps + c);
        const auto *pn_185 = buffer.data(pn + 185 * ncomps + c);
        const auto *pn_186 = buffer.data(pn + 186 * ncomps + c);
        const auto *pn_188 = buffer.data(pn + 188 * ncomps + c);
        const auto *pn_189 = buffer.data(pn + 189 * ncomps + c);
        const auto *pn_190 = buffer.data(pn + 190 * ncomps + c);
        const auto *pn_191 = buffer.data(pn + 191 * ncomps + c);
        const auto *pn_192 = buffer.data(pn + 192 * ncomps + c);
        const auto *pn_193 = buffer.data(pn + 193 * ncomps + c);
        const auto *pn_194 = buffer.data(pn + 194 * ncomps + c);
        const auto *pn_195 = buffer.data(pn + 195 * ncomps + c);
        const auto *pn_196 = buffer.data(pn + 196 * ncomps + c);
        const auto *pn_197 = buffer.data(pn + 197 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_z, pm_125, pm_126, pm_127, \
                         pm_128, pm_129, pn_154, pn_155, pn_156, pn_157, \
                         pn_158 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_z[k] * pm_125[k]
                       + pn_154[k];

            t_291[k] = -ab_z[k] * pm_126[k]
                       + pn_155[k];

            t_292[k] = -ab_z[k] * pm_127[k]
                       + pn_156[k];

            t_293[k] = -ab_z[k] * pm_128[k]
                       + pn_157[k];

            t_294[k] = -ab_z[k] * pm_129[k]
                       + pn_158[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_z, pm_130, pm_131, pm_132, \
                         pm_133, pm_134, pn_159, pn_161, pn_162, pn_163, \
                         pn_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_z[k] * pm_130[k]
                       + pn_159[k];

            t_296[k] = -ab_z[k] * pm_131[k]
                       + pn_161[k];

            t_297[k] = -ab_z[k] * pm_132[k]
                       + pn_162[k];

            t_298[k] = -ab_z[k] * pm_133[k]
                       + pn_163[k];

            t_299[k] = -ab_z[k] * pm_134[k]
                       + pn_164[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_z, pm_135, pm_136, pm_137, \
                         pm_138, pm_139, pn_165, pn_166, pn_167, pn_169, \
                         pn_170 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_z[k] * pm_135[k]
                       + pn_165[k];

            t_301[k] = -ab_z[k] * pm_136[k]
                       + pn_166[k];

            t_302[k] = -ab_z[k] * pm_137[k]
                       + pn_167[k];

            t_303[k] = -ab_z[k] * pm_138[k]
                       + pn_169[k];

            t_304[k] = -ab_z[k] * pm_139[k]
                       + pn_170[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_z, pm_140, pm_141, pm_142, \
                         pm_143, pm_144, pn_171, pn_172, pn_173, pn_174, \
                         pn_175 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_z[k] * pm_140[k]
                       + pn_171[k];

            t_306[k] = -ab_z[k] * pm_141[k]
                       + pn_172[k];

            t_307[k] = -ab_z[k] * pm_142[k]
                       + pn_173[k];

            t_308[k] = -ab_z[k] * pm_143[k]
                       + pn_174[k];

            t_309[k] = -ab_z[k] * pm_144[k]
                       + pn_175[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_z, pm_145, pm_146, pm_147, \
                         pm_148, pm_149, pn_176, pn_178, pn_179, pn_180, \
                         pn_181 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_z[k] * pm_145[k]
                       + pn_176[k];

            t_311[k] = -ab_z[k] * pm_146[k]
                       + pn_178[k];

            t_312[k] = -ab_z[k] * pm_147[k]
                       + pn_179[k];

            t_313[k] = -ab_z[k] * pm_148[k]
                       + pn_180[k];

            t_314[k] = -ab_z[k] * pm_149[k]
                       + pn_181[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_z, pm_150, pm_151, pm_152, \
                         pm_153, pm_154, pn_182, pn_183, pn_184, pn_185, \
                         pn_186 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_z[k] * pm_150[k]
                       + pn_182[k];

            t_316[k] = -ab_z[k] * pm_151[k]
                       + pn_183[k];

            t_317[k] = -ab_z[k] * pm_152[k]
                       + pn_184[k];

            t_318[k] = -ab_z[k] * pm_153[k]
                       + pn_185[k];

            t_319[k] = -ab_z[k] * pm_154[k]
                       + pn_186[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_z, pm_155, pm_156, pm_157, \
                         pm_158, pm_159, pn_188, pn_189, pn_190, pn_191, \
                         pn_192 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_z[k] * pm_155[k]
                       + pn_188[k];

            t_321[k] = -ab_z[k] * pm_156[k]
                       + pn_189[k];

            t_322[k] = -ab_z[k] * pm_157[k]
                       + pn_190[k];

            t_323[k] = -ab_z[k] * pm_158[k]
                       + pn_191[k];

            t_324[k] = -ab_z[k] * pm_159[k]
                       + pn_192[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_z, pm_160, pm_161, pm_162, \
                         pm_163, pm_164, pn_193, pn_194, pn_195, pn_196, \
                         pn_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_z[k] * pm_160[k]
                       + pn_193[k];

            t_326[k] = -ab_z[k] * pm_161[k]
                       + pn_194[k];

            t_327[k] = -ab_z[k] * pm_162[k]
                       + pn_195[k];

            t_328[k] = -ab_z[k] * pm_163[k]
                       + pn_196[k];

            t_329[k] = -ab_z[k] * pm_164[k]
                       + pn_197[k];
        }
    }
}

auto
compute_hrr_dm(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pm, const size_t pn, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_dm_piece0(buffer, coordinates, target, pm, pn, ncomps, nmax);

    compute_hrr_dm_piece1(buffer, coordinates, target, pm, pn, ncomps, nmax);

    compute_hrr_dm_piece2(buffer, coordinates, target, pm, pn, ncomps, nmax);
}

}  // namespace simdtrf
