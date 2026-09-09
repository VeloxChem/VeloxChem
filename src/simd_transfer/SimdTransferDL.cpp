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


#include "SimdTransferDL.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_dl_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pl, const size_t pm, const size_t ncomps,
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
        const auto *ab_y = coordinates.data(7);

        const auto *pl_0 = buffer.data(pl + 0 * ncomps + c);
        const auto *pl_1 = buffer.data(pl + 1 * ncomps + c);
        const auto *pl_2 = buffer.data(pl + 2 * ncomps + c);
        const auto *pl_3 = buffer.data(pl + 3 * ncomps + c);
        const auto *pl_4 = buffer.data(pl + 4 * ncomps + c);
        const auto *pl_5 = buffer.data(pl + 5 * ncomps + c);
        const auto *pl_6 = buffer.data(pl + 6 * ncomps + c);
        const auto *pl_7 = buffer.data(pl + 7 * ncomps + c);
        const auto *pl_8 = buffer.data(pl + 8 * ncomps + c);
        const auto *pl_9 = buffer.data(pl + 9 * ncomps + c);
        const auto *pl_10 = buffer.data(pl + 10 * ncomps + c);
        const auto *pl_11 = buffer.data(pl + 11 * ncomps + c);
        const auto *pl_12 = buffer.data(pl + 12 * ncomps + c);
        const auto *pl_13 = buffer.data(pl + 13 * ncomps + c);
        const auto *pl_14 = buffer.data(pl + 14 * ncomps + c);
        const auto *pl_15 = buffer.data(pl + 15 * ncomps + c);
        const auto *pl_16 = buffer.data(pl + 16 * ncomps + c);
        const auto *pl_17 = buffer.data(pl + 17 * ncomps + c);
        const auto *pl_18 = buffer.data(pl + 18 * ncomps + c);
        const auto *pl_19 = buffer.data(pl + 19 * ncomps + c);
        const auto *pl_20 = buffer.data(pl + 20 * ncomps + c);
        const auto *pl_21 = buffer.data(pl + 21 * ncomps + c);
        const auto *pl_22 = buffer.data(pl + 22 * ncomps + c);
        const auto *pl_23 = buffer.data(pl + 23 * ncomps + c);
        const auto *pl_24 = buffer.data(pl + 24 * ncomps + c);
        const auto *pl_25 = buffer.data(pl + 25 * ncomps + c);
        const auto *pl_26 = buffer.data(pl + 26 * ncomps + c);
        const auto *pl_27 = buffer.data(pl + 27 * ncomps + c);
        const auto *pl_28 = buffer.data(pl + 28 * ncomps + c);
        const auto *pl_29 = buffer.data(pl + 29 * ncomps + c);
        const auto *pl_30 = buffer.data(pl + 30 * ncomps + c);
        const auto *pl_31 = buffer.data(pl + 31 * ncomps + c);
        const auto *pl_32 = buffer.data(pl + 32 * ncomps + c);
        const auto *pl_33 = buffer.data(pl + 33 * ncomps + c);
        const auto *pl_34 = buffer.data(pl + 34 * ncomps + c);
        const auto *pl_35 = buffer.data(pl + 35 * ncomps + c);
        const auto *pl_36 = buffer.data(pl + 36 * ncomps + c);
        const auto *pl_37 = buffer.data(pl + 37 * ncomps + c);
        const auto *pl_38 = buffer.data(pl + 38 * ncomps + c);
        const auto *pl_39 = buffer.data(pl + 39 * ncomps + c);
        const auto *pl_40 = buffer.data(pl + 40 * ncomps + c);
        const auto *pl_41 = buffer.data(pl + 41 * ncomps + c);
        const auto *pl_42 = buffer.data(pl + 42 * ncomps + c);
        const auto *pl_43 = buffer.data(pl + 43 * ncomps + c);
        const auto *pl_44 = buffer.data(pl + 44 * ncomps + c);
        const auto *pl_45 = buffer.data(pl + 45 * ncomps + c);
        const auto *pl_46 = buffer.data(pl + 46 * ncomps + c);
        const auto *pl_47 = buffer.data(pl + 47 * ncomps + c);
        const auto *pl_48 = buffer.data(pl + 48 * ncomps + c);
        const auto *pl_49 = buffer.data(pl + 49 * ncomps + c);
        const auto *pl_50 = buffer.data(pl + 50 * ncomps + c);
        const auto *pl_51 = buffer.data(pl + 51 * ncomps + c);
        const auto *pl_52 = buffer.data(pl + 52 * ncomps + c);
        const auto *pl_53 = buffer.data(pl + 53 * ncomps + c);
        const auto *pl_54 = buffer.data(pl + 54 * ncomps + c);
        const auto *pl_55 = buffer.data(pl + 55 * ncomps + c);
        const auto *pl_56 = buffer.data(pl + 56 * ncomps + c);
        const auto *pl_57 = buffer.data(pl + 57 * ncomps + c);
        const auto *pl_58 = buffer.data(pl + 58 * ncomps + c);
        const auto *pl_59 = buffer.data(pl + 59 * ncomps + c);
        const auto *pl_60 = buffer.data(pl + 60 * ncomps + c);
        const auto *pl_61 = buffer.data(pl + 61 * ncomps + c);
        const auto *pl_62 = buffer.data(pl + 62 * ncomps + c);
        const auto *pl_63 = buffer.data(pl + 63 * ncomps + c);
        const auto *pl_64 = buffer.data(pl + 64 * ncomps + c);
        const auto *pl_65 = buffer.data(pl + 65 * ncomps + c);
        const auto *pl_66 = buffer.data(pl + 66 * ncomps + c);
        const auto *pl_67 = buffer.data(pl + 67 * ncomps + c);
        const auto *pl_68 = buffer.data(pl + 68 * ncomps + c);
        const auto *pl_69 = buffer.data(pl + 69 * ncomps + c);
        const auto *pl_70 = buffer.data(pl + 70 * ncomps + c);
        const auto *pl_71 = buffer.data(pl + 71 * ncomps + c);
        const auto *pl_72 = buffer.data(pl + 72 * ncomps + c);
        const auto *pl_73 = buffer.data(pl + 73 * ncomps + c);
        const auto *pl_74 = buffer.data(pl + 74 * ncomps + c);
        const auto *pl_75 = buffer.data(pl + 75 * ncomps + c);
        const auto *pl_76 = buffer.data(pl + 76 * ncomps + c);
        const auto *pl_77 = buffer.data(pl + 77 * ncomps + c);
        const auto *pl_78 = buffer.data(pl + 78 * ncomps + c);
        const auto *pl_79 = buffer.data(pl + 79 * ncomps + c);
        const auto *pl_80 = buffer.data(pl + 80 * ncomps + c);
        const auto *pl_81 = buffer.data(pl + 81 * ncomps + c);
        const auto *pl_82 = buffer.data(pl + 82 * ncomps + c);
        const auto *pl_83 = buffer.data(pl + 83 * ncomps + c);
        const auto *pl_84 = buffer.data(pl + 84 * ncomps + c);
        const auto *pl_85 = buffer.data(pl + 85 * ncomps + c);
        const auto *pl_86 = buffer.data(pl + 86 * ncomps + c);
        const auto *pl_87 = buffer.data(pl + 87 * ncomps + c);
        const auto *pl_88 = buffer.data(pl + 88 * ncomps + c);
        const auto *pl_89 = buffer.data(pl + 89 * ncomps + c);
        const auto *pl_90 = buffer.data(pl + 90 * ncomps + c);
        const auto *pl_91 = buffer.data(pl + 91 * ncomps + c);
        const auto *pl_92 = buffer.data(pl + 92 * ncomps + c);
        const auto *pl_93 = buffer.data(pl + 93 * ncomps + c);
        const auto *pl_94 = buffer.data(pl + 94 * ncomps + c);
        const auto *pl_95 = buffer.data(pl + 95 * ncomps + c);
        const auto *pl_96 = buffer.data(pl + 96 * ncomps + c);
        const auto *pl_97 = buffer.data(pl + 97 * ncomps + c);
        const auto *pl_98 = buffer.data(pl + 98 * ncomps + c);
        const auto *pl_99 = buffer.data(pl + 99 * ncomps + c);
        const auto *pl_100 = buffer.data(pl + 100 * ncomps + c);
        const auto *pl_101 = buffer.data(pl + 101 * ncomps + c);
        const auto *pl_102 = buffer.data(pl + 102 * ncomps + c);
        const auto *pl_103 = buffer.data(pl + 103 * ncomps + c);
        const auto *pl_104 = buffer.data(pl + 104 * ncomps + c);
        const auto *pl_105 = buffer.data(pl + 105 * ncomps + c);
        const auto *pl_106 = buffer.data(pl + 106 * ncomps + c);
        const auto *pl_107 = buffer.data(pl + 107 * ncomps + c);
        const auto *pl_108 = buffer.data(pl + 108 * ncomps + c);
        const auto *pl_109 = buffer.data(pl + 109 * ncomps + c);
        const auto *pl_110 = buffer.data(pl + 110 * ncomps + c);
        const auto *pl_111 = buffer.data(pl + 111 * ncomps + c);
        const auto *pl_112 = buffer.data(pl + 112 * ncomps + c);
        const auto *pl_113 = buffer.data(pl + 113 * ncomps + c);
        const auto *pl_114 = buffer.data(pl + 114 * ncomps + c);
        const auto *pl_115 = buffer.data(pl + 115 * ncomps + c);
        const auto *pl_116 = buffer.data(pl + 116 * ncomps + c);
        const auto *pl_117 = buffer.data(pl + 117 * ncomps + c);
        const auto *pl_118 = buffer.data(pl + 118 * ncomps + c);
        const auto *pl_119 = buffer.data(pl + 119 * ncomps + c);
        const auto *pl_120 = buffer.data(pl + 120 * ncomps + c);
        const auto *pl_121 = buffer.data(pl + 121 * ncomps + c);
        const auto *pl_122 = buffer.data(pl + 122 * ncomps + c);
        const auto *pl_123 = buffer.data(pl + 123 * ncomps + c);
        const auto *pl_124 = buffer.data(pl + 124 * ncomps + c);
        const auto *pl_125 = buffer.data(pl + 125 * ncomps + c);
        const auto *pl_126 = buffer.data(pl + 126 * ncomps + c);
        const auto *pl_127 = buffer.data(pl + 127 * ncomps + c);
        const auto *pl_128 = buffer.data(pl + 128 * ncomps + c);
        const auto *pl_129 = buffer.data(pl + 129 * ncomps + c);
        const auto *pl_130 = buffer.data(pl + 130 * ncomps + c);
        const auto *pl_131 = buffer.data(pl + 131 * ncomps + c);
        const auto *pl_132 = buffer.data(pl + 132 * ncomps + c);
        const auto *pl_133 = buffer.data(pl + 133 * ncomps + c);
        const auto *pl_134 = buffer.data(pl + 134 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pl_0, pl_1, pl_2, pl_3, pl_4, pm_0, \
                         pm_1, pm_2, pm_3, pm_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pl_0[k]
                     + pm_0[k];

            t_1[k] = -ab_x[k] * pl_1[k]
                     + pm_1[k];

            t_2[k] = -ab_x[k] * pl_2[k]
                     + pm_2[k];

            t_3[k] = -ab_x[k] * pl_3[k]
                     + pm_3[k];

            t_4[k] = -ab_x[k] * pl_4[k]
                     + pm_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pl_5, pl_6, pl_7, pl_8, pl_9, pm_5, \
                         pm_6, pm_7, pm_8, pm_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pl_5[k]
                     + pm_5[k];

            t_6[k] = -ab_x[k] * pl_6[k]
                     + pm_6[k];

            t_7[k] = -ab_x[k] * pl_7[k]
                     + pm_7[k];

            t_8[k] = -ab_x[k] * pl_8[k]
                     + pm_8[k];

            t_9[k] = -ab_x[k] * pl_9[k]
                     + pm_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pl_10, pl_11, pl_12, pl_13, \
                         pl_14, pm_10, pm_11, pm_12, pm_13, pm_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pl_10[k]
                      + pm_10[k];

            t_11[k] = -ab_x[k] * pl_11[k]
                      + pm_11[k];

            t_12[k] = -ab_x[k] * pl_12[k]
                      + pm_12[k];

            t_13[k] = -ab_x[k] * pl_13[k]
                      + pm_13[k];

            t_14[k] = -ab_x[k] * pl_14[k]
                      + pm_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pl_15, pl_16, pl_17, pl_18, \
                         pl_19, pm_15, pm_16, pm_17, pm_18, pm_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pl_15[k]
                      + pm_15[k];

            t_16[k] = -ab_x[k] * pl_16[k]
                      + pm_16[k];

            t_17[k] = -ab_x[k] * pl_17[k]
                      + pm_17[k];

            t_18[k] = -ab_x[k] * pl_18[k]
                      + pm_18[k];

            t_19[k] = -ab_x[k] * pl_19[k]
                      + pm_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pl_20, pl_21, pl_22, pl_23, \
                         pl_24, pm_20, pm_21, pm_22, pm_23, pm_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * pl_20[k]
                      + pm_20[k];

            t_21[k] = -ab_x[k] * pl_21[k]
                      + pm_21[k];

            t_22[k] = -ab_x[k] * pl_22[k]
                      + pm_22[k];

            t_23[k] = -ab_x[k] * pl_23[k]
                      + pm_23[k];

            t_24[k] = -ab_x[k] * pl_24[k]
                      + pm_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pl_25, pl_26, pl_27, pl_28, \
                         pl_29, pm_25, pm_26, pm_27, pm_28, pm_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * pl_25[k]
                      + pm_25[k];

            t_26[k] = -ab_x[k] * pl_26[k]
                      + pm_26[k];

            t_27[k] = -ab_x[k] * pl_27[k]
                      + pm_27[k];

            t_28[k] = -ab_x[k] * pl_28[k]
                      + pm_28[k];

            t_29[k] = -ab_x[k] * pl_29[k]
                      + pm_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pl_30, pl_31, pl_32, pl_33, \
                         pl_34, pm_30, pm_31, pm_32, pm_33, pm_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * pl_30[k]
                      + pm_30[k];

            t_31[k] = -ab_x[k] * pl_31[k]
                      + pm_31[k];

            t_32[k] = -ab_x[k] * pl_32[k]
                      + pm_32[k];

            t_33[k] = -ab_x[k] * pl_33[k]
                      + pm_33[k];

            t_34[k] = -ab_x[k] * pl_34[k]
                      + pm_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pl_35, pl_36, pl_37, pl_38, \
                         pl_39, pm_35, pm_36, pm_37, pm_38, pm_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * pl_35[k]
                      + pm_35[k];

            t_36[k] = -ab_x[k] * pl_36[k]
                      + pm_36[k];

            t_37[k] = -ab_x[k] * pl_37[k]
                      + pm_37[k];

            t_38[k] = -ab_x[k] * pl_38[k]
                      + pm_38[k];

            t_39[k] = -ab_x[k] * pl_39[k]
                      + pm_39[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pl_40, pl_41, pl_42, pl_43, \
                         pl_44, pm_40, pm_41, pm_42, pm_43, pm_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * pl_40[k]
                      + pm_40[k];

            t_41[k] = -ab_x[k] * pl_41[k]
                      + pm_41[k];

            t_42[k] = -ab_x[k] * pl_42[k]
                      + pm_42[k];

            t_43[k] = -ab_x[k] * pl_43[k]
                      + pm_43[k];

            t_44[k] = -ab_x[k] * pl_44[k]
                      + pm_44[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, pl_45, pl_46, pl_47, pl_48, \
                         pl_49, pm_55, pm_56, pm_57, pm_58, pm_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * pl_45[k]
                      + pm_55[k];

            t_46[k] = -ab_x[k] * pl_46[k]
                      + pm_56[k];

            t_47[k] = -ab_x[k] * pl_47[k]
                      + pm_57[k];

            t_48[k] = -ab_x[k] * pl_48[k]
                      + pm_58[k];

            t_49[k] = -ab_x[k] * pl_49[k]
                      + pm_59[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, pl_50, pl_51, pl_52, pl_53, \
                         pl_54, pm_60, pm_61, pm_62, pm_63, pm_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * pl_50[k]
                      + pm_60[k];

            t_51[k] = -ab_x[k] * pl_51[k]
                      + pm_61[k];

            t_52[k] = -ab_x[k] * pl_52[k]
                      + pm_62[k];

            t_53[k] = -ab_x[k] * pl_53[k]
                      + pm_63[k];

            t_54[k] = -ab_x[k] * pl_54[k]
                      + pm_64[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, pl_55, pl_56, pl_57, pl_58, \
                         pl_59, pm_65, pm_66, pm_67, pm_68, pm_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * pl_55[k]
                      + pm_65[k];

            t_56[k] = -ab_x[k] * pl_56[k]
                      + pm_66[k];

            t_57[k] = -ab_x[k] * pl_57[k]
                      + pm_67[k];

            t_58[k] = -ab_x[k] * pl_58[k]
                      + pm_68[k];

            t_59[k] = -ab_x[k] * pl_59[k]
                      + pm_69[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, pl_60, pl_61, pl_62, pl_63, \
                         pl_64, pm_70, pm_71, pm_72, pm_73, pm_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * pl_60[k]
                      + pm_70[k];

            t_61[k] = -ab_x[k] * pl_61[k]
                      + pm_71[k];

            t_62[k] = -ab_x[k] * pl_62[k]
                      + pm_72[k];

            t_63[k] = -ab_x[k] * pl_63[k]
                      + pm_73[k];

            t_64[k] = -ab_x[k] * pl_64[k]
                      + pm_74[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, pl_65, pl_66, pl_67, pl_68, \
                         pl_69, pm_75, pm_76, pm_77, pm_78, pm_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * pl_65[k]
                      + pm_75[k];

            t_66[k] = -ab_x[k] * pl_66[k]
                      + pm_76[k];

            t_67[k] = -ab_x[k] * pl_67[k]
                      + pm_77[k];

            t_68[k] = -ab_x[k] * pl_68[k]
                      + pm_78[k];

            t_69[k] = -ab_x[k] * pl_69[k]
                      + pm_79[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, pl_70, pl_71, pl_72, pl_73, \
                         pl_74, pm_80, pm_81, pm_82, pm_83, pm_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * pl_70[k]
                      + pm_80[k];

            t_71[k] = -ab_x[k] * pl_71[k]
                      + pm_81[k];

            t_72[k] = -ab_x[k] * pl_72[k]
                      + pm_82[k];

            t_73[k] = -ab_x[k] * pl_73[k]
                      + pm_83[k];

            t_74[k] = -ab_x[k] * pl_74[k]
                      + pm_84[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, pl_75, pl_76, pl_77, pl_78, \
                         pl_79, pm_85, pm_86, pm_87, pm_88, pm_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * pl_75[k]
                      + pm_85[k];

            t_76[k] = -ab_x[k] * pl_76[k]
                      + pm_86[k];

            t_77[k] = -ab_x[k] * pl_77[k]
                      + pm_87[k];

            t_78[k] = -ab_x[k] * pl_78[k]
                      + pm_88[k];

            t_79[k] = -ab_x[k] * pl_79[k]
                      + pm_89[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, pl_80, pl_81, pl_82, pl_83, \
                         pl_84, pm_90, pm_91, pm_92, pm_93, pm_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * pl_80[k]
                      + pm_90[k];

            t_81[k] = -ab_x[k] * pl_81[k]
                      + pm_91[k];

            t_82[k] = -ab_x[k] * pl_82[k]
                      + pm_92[k];

            t_83[k] = -ab_x[k] * pl_83[k]
                      + pm_93[k];

            t_84[k] = -ab_x[k] * pl_84[k]
                      + pm_94[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, pl_85, pl_86, pl_87, pl_88, \
                         pl_89, pm_95, pm_96, pm_97, pm_98, pm_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * pl_85[k]
                      + pm_95[k];

            t_86[k] = -ab_x[k] * pl_86[k]
                      + pm_96[k];

            t_87[k] = -ab_x[k] * pl_87[k]
                      + pm_97[k];

            t_88[k] = -ab_x[k] * pl_88[k]
                      + pm_98[k];

            t_89[k] = -ab_x[k] * pl_89[k]
                      + pm_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, pl_90, pl_91, pl_92, pl_93, \
                         pl_94, pm_110, pm_111, pm_112, pm_113, \
                         pm_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * pl_90[k]
                      + pm_110[k];

            t_91[k] = -ab_x[k] * pl_91[k]
                      + pm_111[k];

            t_92[k] = -ab_x[k] * pl_92[k]
                      + pm_112[k];

            t_93[k] = -ab_x[k] * pl_93[k]
                      + pm_113[k];

            t_94[k] = -ab_x[k] * pl_94[k]
                      + pm_114[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, pl_95, pl_96, pl_97, pl_98, \
                         pl_99, pm_115, pm_116, pm_117, pm_118, \
                         pm_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * pl_95[k]
                      + pm_115[k];

            t_96[k] = -ab_x[k] * pl_96[k]
                      + pm_116[k];

            t_97[k] = -ab_x[k] * pl_97[k]
                      + pm_117[k];

            t_98[k] = -ab_x[k] * pl_98[k]
                      + pm_118[k];

            t_99[k] = -ab_x[k] * pl_99[k]
                      + pm_119[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, pl_100, pl_101, pl_102, \
                         pl_103, pl_104, pm_120, pm_121, pm_122, pm_123, \
                         pm_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * pl_100[k]
                       + pm_120[k];

            t_101[k] = -ab_x[k] * pl_101[k]
                       + pm_121[k];

            t_102[k] = -ab_x[k] * pl_102[k]
                       + pm_122[k];

            t_103[k] = -ab_x[k] * pl_103[k]
                       + pm_123[k];

            t_104[k] = -ab_x[k] * pl_104[k]
                       + pm_124[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, pl_105, pl_106, pl_107, \
                         pl_108, pl_109, pm_125, pm_126, pm_127, pm_128, \
                         pm_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * pl_105[k]
                       + pm_125[k];

            t_106[k] = -ab_x[k] * pl_106[k]
                       + pm_126[k];

            t_107[k] = -ab_x[k] * pl_107[k]
                       + pm_127[k];

            t_108[k] = -ab_x[k] * pl_108[k]
                       + pm_128[k];

            t_109[k] = -ab_x[k] * pl_109[k]
                       + pm_129[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, pl_110, pl_111, pl_112, \
                         pl_113, pl_114, pm_130, pm_131, pm_132, pm_133, \
                         pm_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * pl_110[k]
                       + pm_130[k];

            t_111[k] = -ab_x[k] * pl_111[k]
                       + pm_131[k];

            t_112[k] = -ab_x[k] * pl_112[k]
                       + pm_132[k];

            t_113[k] = -ab_x[k] * pl_113[k]
                       + pm_133[k];

            t_114[k] = -ab_x[k] * pl_114[k]
                       + pm_134[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, pl_115, pl_116, pl_117, \
                         pl_118, pl_119, pm_135, pm_136, pm_137, pm_138, \
                         pm_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * pl_115[k]
                       + pm_135[k];

            t_116[k] = -ab_x[k] * pl_116[k]
                       + pm_136[k];

            t_117[k] = -ab_x[k] * pl_117[k]
                       + pm_137[k];

            t_118[k] = -ab_x[k] * pl_118[k]
                       + pm_138[k];

            t_119[k] = -ab_x[k] * pl_119[k]
                       + pm_139[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, pl_120, pl_121, pl_122, \
                         pl_123, pl_124, pm_140, pm_141, pm_142, pm_143, \
                         pm_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * pl_120[k]
                       + pm_140[k];

            t_121[k] = -ab_x[k] * pl_121[k]
                       + pm_141[k];

            t_122[k] = -ab_x[k] * pl_122[k]
                       + pm_142[k];

            t_123[k] = -ab_x[k] * pl_123[k]
                       + pm_143[k];

            t_124[k] = -ab_x[k] * pl_124[k]
                       + pm_144[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, pl_125, pl_126, pl_127, \
                         pl_128, pl_129, pm_145, pm_146, pm_147, pm_148, \
                         pm_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * pl_125[k]
                       + pm_145[k];

            t_126[k] = -ab_x[k] * pl_126[k]
                       + pm_146[k];

            t_127[k] = -ab_x[k] * pl_127[k]
                       + pm_147[k];

            t_128[k] = -ab_x[k] * pl_128[k]
                       + pm_148[k];

            t_129[k] = -ab_x[k] * pl_129[k]
                       + pm_149[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, pl_130, pl_131, pl_132, \
                         pl_133, pl_134, pm_150, pm_151, pm_152, pm_153, \
                         pm_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * pl_130[k]
                       + pm_150[k];

            t_131[k] = -ab_x[k] * pl_131[k]
                       + pm_151[k];

            t_132[k] = -ab_x[k] * pl_132[k]
                       + pm_152[k];

            t_133[k] = -ab_x[k] * pl_133[k]
                       + pm_153[k];

            t_134[k] = -ab_x[k] * pl_134[k]
                       + pm_154[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_y, pl_45, pl_46, pl_47, pl_48, \
                         pl_49, pm_56, pm_58, pm_59, pm_61, pm_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_y[k] * pl_45[k]
                       + pm_56[k];

            t_136[k] = -ab_y[k] * pl_46[k]
                       + pm_58[k];

            t_137[k] = -ab_y[k] * pl_47[k]
                       + pm_59[k];

            t_138[k] = -ab_y[k] * pl_48[k]
                       + pm_61[k];

            t_139[k] = -ab_y[k] * pl_49[k]
                       + pm_62[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_y, pl_50, pl_51, pl_52, pl_53, \
                         pl_54, pm_63, pm_65, pm_66, pm_67, pm_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_y[k] * pl_50[k]
                       + pm_63[k];

            t_141[k] = -ab_y[k] * pl_51[k]
                       + pm_65[k];

            t_142[k] = -ab_y[k] * pl_52[k]
                       + pm_66[k];

            t_143[k] = -ab_y[k] * pl_53[k]
                       + pm_67[k];

            t_144[k] = -ab_y[k] * pl_54[k]
                       + pm_68[k];
        }
    }
}

static auto
compute_hrr_dl_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pl, const size_t pm, const size_t ncomps,
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pl_55 = buffer.data(pl + 55 * ncomps + c);
        const auto *pl_56 = buffer.data(pl + 56 * ncomps + c);
        const auto *pl_57 = buffer.data(pl + 57 * ncomps + c);
        const auto *pl_58 = buffer.data(pl + 58 * ncomps + c);
        const auto *pl_59 = buffer.data(pl + 59 * ncomps + c);
        const auto *pl_60 = buffer.data(pl + 60 * ncomps + c);
        const auto *pl_61 = buffer.data(pl + 61 * ncomps + c);
        const auto *pl_62 = buffer.data(pl + 62 * ncomps + c);
        const auto *pl_63 = buffer.data(pl + 63 * ncomps + c);
        const auto *pl_64 = buffer.data(pl + 64 * ncomps + c);
        const auto *pl_65 = buffer.data(pl + 65 * ncomps + c);
        const auto *pl_66 = buffer.data(pl + 66 * ncomps + c);
        const auto *pl_67 = buffer.data(pl + 67 * ncomps + c);
        const auto *pl_68 = buffer.data(pl + 68 * ncomps + c);
        const auto *pl_69 = buffer.data(pl + 69 * ncomps + c);
        const auto *pl_70 = buffer.data(pl + 70 * ncomps + c);
        const auto *pl_71 = buffer.data(pl + 71 * ncomps + c);
        const auto *pl_72 = buffer.data(pl + 72 * ncomps + c);
        const auto *pl_73 = buffer.data(pl + 73 * ncomps + c);
        const auto *pl_74 = buffer.data(pl + 74 * ncomps + c);
        const auto *pl_75 = buffer.data(pl + 75 * ncomps + c);
        const auto *pl_76 = buffer.data(pl + 76 * ncomps + c);
        const auto *pl_77 = buffer.data(pl + 77 * ncomps + c);
        const auto *pl_78 = buffer.data(pl + 78 * ncomps + c);
        const auto *pl_79 = buffer.data(pl + 79 * ncomps + c);
        const auto *pl_80 = buffer.data(pl + 80 * ncomps + c);
        const auto *pl_81 = buffer.data(pl + 81 * ncomps + c);
        const auto *pl_82 = buffer.data(pl + 82 * ncomps + c);
        const auto *pl_83 = buffer.data(pl + 83 * ncomps + c);
        const auto *pl_84 = buffer.data(pl + 84 * ncomps + c);
        const auto *pl_85 = buffer.data(pl + 85 * ncomps + c);
        const auto *pl_86 = buffer.data(pl + 86 * ncomps + c);
        const auto *pl_87 = buffer.data(pl + 87 * ncomps + c);
        const auto *pl_88 = buffer.data(pl + 88 * ncomps + c);
        const auto *pl_89 = buffer.data(pl + 89 * ncomps + c);
        const auto *pl_90 = buffer.data(pl + 90 * ncomps + c);
        const auto *pl_91 = buffer.data(pl + 91 * ncomps + c);
        const auto *pl_92 = buffer.data(pl + 92 * ncomps + c);
        const auto *pl_93 = buffer.data(pl + 93 * ncomps + c);
        const auto *pl_94 = buffer.data(pl + 94 * ncomps + c);
        const auto *pl_95 = buffer.data(pl + 95 * ncomps + c);
        const auto *pl_96 = buffer.data(pl + 96 * ncomps + c);
        const auto *pl_97 = buffer.data(pl + 97 * ncomps + c);
        const auto *pl_98 = buffer.data(pl + 98 * ncomps + c);
        const auto *pl_99 = buffer.data(pl + 99 * ncomps + c);
        const auto *pl_100 = buffer.data(pl + 100 * ncomps + c);
        const auto *pl_101 = buffer.data(pl + 101 * ncomps + c);
        const auto *pl_102 = buffer.data(pl + 102 * ncomps + c);
        const auto *pl_103 = buffer.data(pl + 103 * ncomps + c);
        const auto *pl_104 = buffer.data(pl + 104 * ncomps + c);
        const auto *pl_105 = buffer.data(pl + 105 * ncomps + c);
        const auto *pl_106 = buffer.data(pl + 106 * ncomps + c);
        const auto *pl_107 = buffer.data(pl + 107 * ncomps + c);
        const auto *pl_108 = buffer.data(pl + 108 * ncomps + c);
        const auto *pl_109 = buffer.data(pl + 109 * ncomps + c);
        const auto *pl_110 = buffer.data(pl + 110 * ncomps + c);
        const auto *pl_111 = buffer.data(pl + 111 * ncomps + c);
        const auto *pl_112 = buffer.data(pl + 112 * ncomps + c);
        const auto *pl_113 = buffer.data(pl + 113 * ncomps + c);
        const auto *pl_114 = buffer.data(pl + 114 * ncomps + c);
        const auto *pl_115 = buffer.data(pl + 115 * ncomps + c);
        const auto *pl_116 = buffer.data(pl + 116 * ncomps + c);
        const auto *pl_117 = buffer.data(pl + 117 * ncomps + c);
        const auto *pl_118 = buffer.data(pl + 118 * ncomps + c);
        const auto *pl_119 = buffer.data(pl + 119 * ncomps + c);
        const auto *pl_120 = buffer.data(pl + 120 * ncomps + c);
        const auto *pl_121 = buffer.data(pl + 121 * ncomps + c);
        const auto *pl_122 = buffer.data(pl + 122 * ncomps + c);
        const auto *pl_123 = buffer.data(pl + 123 * ncomps + c);
        const auto *pl_124 = buffer.data(pl + 124 * ncomps + c);
        const auto *pl_125 = buffer.data(pl + 125 * ncomps + c);
        const auto *pl_126 = buffer.data(pl + 126 * ncomps + c);
        const auto *pl_127 = buffer.data(pl + 127 * ncomps + c);
        const auto *pl_128 = buffer.data(pl + 128 * ncomps + c);
        const auto *pl_129 = buffer.data(pl + 129 * ncomps + c);
        const auto *pl_130 = buffer.data(pl + 130 * ncomps + c);
        const auto *pl_131 = buffer.data(pl + 131 * ncomps + c);
        const auto *pl_132 = buffer.data(pl + 132 * ncomps + c);
        const auto *pl_133 = buffer.data(pl + 133 * ncomps + c);
        const auto *pl_134 = buffer.data(pl + 134 * ncomps + c);

        const auto *pm_70 = buffer.data(pm + 70 * ncomps + c);
        const auto *pm_71 = buffer.data(pm + 71 * ncomps + c);
        const auto *pm_72 = buffer.data(pm + 72 * ncomps + c);
        const auto *pm_73 = buffer.data(pm + 73 * ncomps + c);
        const auto *pm_74 = buffer.data(pm + 74 * ncomps + c);
        const auto *pm_76 = buffer.data(pm + 76 * ncomps + c);
        const auto *pm_77 = buffer.data(pm + 77 * ncomps + c);
        const auto *pm_78 = buffer.data(pm + 78 * ncomps + c);
        const auto *pm_79 = buffer.data(pm + 79 * ncomps + c);
        const auto *pm_80 = buffer.data(pm + 80 * ncomps + c);
        const auto *pm_81 = buffer.data(pm + 81 * ncomps + c);
        const auto *pm_83 = buffer.data(pm + 83 * ncomps + c);
        const auto *pm_84 = buffer.data(pm + 84 * ncomps + c);
        const auto *pm_85 = buffer.data(pm + 85 * ncomps + c);
        const auto *pm_86 = buffer.data(pm + 86 * ncomps + c);
        const auto *pm_87 = buffer.data(pm + 87 * ncomps + c);
        const auto *pm_88 = buffer.data(pm + 88 * ncomps + c);
        const auto *pm_89 = buffer.data(pm + 89 * ncomps + c);
        const auto *pm_91 = buffer.data(pm + 91 * ncomps + c);
        const auto *pm_92 = buffer.data(pm + 92 * ncomps + c);
        const auto *pm_93 = buffer.data(pm + 93 * ncomps + c);
        const auto *pm_94 = buffer.data(pm + 94 * ncomps + c);
        const auto *pm_95 = buffer.data(pm + 95 * ncomps + c);
        const auto *pm_96 = buffer.data(pm + 96 * ncomps + c);
        const auto *pm_97 = buffer.data(pm + 97 * ncomps + c);
        const auto *pm_98 = buffer.data(pm + 98 * ncomps + c);
        const auto *pm_100 = buffer.data(pm + 100 * ncomps + c);
        const auto *pm_101 = buffer.data(pm + 101 * ncomps + c);
        const auto *pm_102 = buffer.data(pm + 102 * ncomps + c);
        const auto *pm_103 = buffer.data(pm + 103 * ncomps + c);
        const auto *pm_104 = buffer.data(pm + 104 * ncomps + c);
        const auto *pm_105 = buffer.data(pm + 105 * ncomps + c);
        const auto *pm_106 = buffer.data(pm + 106 * ncomps + c);
        const auto *pm_107 = buffer.data(pm + 107 * ncomps + c);
        const auto *pm_108 = buffer.data(pm + 108 * ncomps + c);
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

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, pl_55, pl_56, pl_57, pl_58, \
                         pl_59, pm_70, pm_71, pm_72, pm_73, pm_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_y[k] * pl_55[k]
                       + pm_70[k];

            t_146[k] = -ab_y[k] * pl_56[k]
                       + pm_71[k];

            t_147[k] = -ab_y[k] * pl_57[k]
                       + pm_72[k];

            t_148[k] = -ab_y[k] * pl_58[k]
                       + pm_73[k];

            t_149[k] = -ab_y[k] * pl_59[k]
                       + pm_74[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_y, pl_60, pl_61, pl_62, pl_63, \
                         pl_64, pm_76, pm_77, pm_78, pm_79, pm_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_y[k] * pl_60[k]
                       + pm_76[k];

            t_151[k] = -ab_y[k] * pl_61[k]
                       + pm_77[k];

            t_152[k] = -ab_y[k] * pl_62[k]
                       + pm_78[k];

            t_153[k] = -ab_y[k] * pl_63[k]
                       + pm_79[k];

            t_154[k] = -ab_y[k] * pl_64[k]
                       + pm_80[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_y, pl_65, pl_66, pl_67, pl_68, \
                         pl_69, pm_81, pm_83, pm_84, pm_85, pm_86 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_y[k] * pl_65[k]
                       + pm_81[k];

            t_156[k] = -ab_y[k] * pl_66[k]
                       + pm_83[k];

            t_157[k] = -ab_y[k] * pl_67[k]
                       + pm_84[k];

            t_158[k] = -ab_y[k] * pl_68[k]
                       + pm_85[k];

            t_159[k] = -ab_y[k] * pl_69[k]
                       + pm_86[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, pl_70, pl_71, pl_72, pl_73, \
                         pl_74, pm_87, pm_88, pm_89, pm_91, pm_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_y[k] * pl_70[k]
                       + pm_87[k];

            t_161[k] = -ab_y[k] * pl_71[k]
                       + pm_88[k];

            t_162[k] = -ab_y[k] * pl_72[k]
                       + pm_89[k];

            t_163[k] = -ab_y[k] * pl_73[k]
                       + pm_91[k];

            t_164[k] = -ab_y[k] * pl_74[k]
                       + pm_92[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_y, pl_75, pl_76, pl_77, pl_78, \
                         pl_79, pm_93, pm_94, pm_95, pm_96, pm_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_y[k] * pl_75[k]
                       + pm_93[k];

            t_166[k] = -ab_y[k] * pl_76[k]
                       + pm_94[k];

            t_167[k] = -ab_y[k] * pl_77[k]
                       + pm_95[k];

            t_168[k] = -ab_y[k] * pl_78[k]
                       + pm_96[k];

            t_169[k] = -ab_y[k] * pl_79[k]
                       + pm_97[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_y, pl_80, pl_81, pl_82, pl_83, \
                         pl_84, pm_98, pm_100, pm_101, pm_102, pm_103 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_y[k] * pl_80[k]
                       + pm_98[k];

            t_171[k] = -ab_y[k] * pl_81[k]
                       + pm_100[k];

            t_172[k] = -ab_y[k] * pl_82[k]
                       + pm_101[k];

            t_173[k] = -ab_y[k] * pl_83[k]
                       + pm_102[k];

            t_174[k] = -ab_y[k] * pl_84[k]
                       + pm_103[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, pl_85, pl_86, pl_87, pl_88, \
                         pl_89, pm_104, pm_105, pm_106, pm_107, \
                         pm_108 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_y[k] * pl_85[k]
                       + pm_104[k];

            t_176[k] = -ab_y[k] * pl_86[k]
                       + pm_105[k];

            t_177[k] = -ab_y[k] * pl_87[k]
                       + pm_106[k];

            t_178[k] = -ab_y[k] * pl_88[k]
                       + pm_107[k];

            t_179[k] = -ab_y[k] * pl_89[k]
                       + pm_108[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_y, pl_90, pl_91, pl_92, pl_93, \
                         pl_94, pm_111, pm_113, pm_114, pm_116, \
                         pm_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_y[k] * pl_90[k]
                       + pm_111[k];

            t_181[k] = -ab_y[k] * pl_91[k]
                       + pm_113[k];

            t_182[k] = -ab_y[k] * pl_92[k]
                       + pm_114[k];

            t_183[k] = -ab_y[k] * pl_93[k]
                       + pm_116[k];

            t_184[k] = -ab_y[k] * pl_94[k]
                       + pm_117[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_y, pl_95, pl_96, pl_97, pl_98, \
                         pl_99, pm_118, pm_120, pm_121, pm_122, \
                         pm_123 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_y[k] * pl_95[k]
                       + pm_118[k];

            t_186[k] = -ab_y[k] * pl_96[k]
                       + pm_120[k];

            t_187[k] = -ab_y[k] * pl_97[k]
                       + pm_121[k];

            t_188[k] = -ab_y[k] * pl_98[k]
                       + pm_122[k];

            t_189[k] = -ab_y[k] * pl_99[k]
                       + pm_123[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, pl_100, pl_101, pl_102, \
                         pl_103, pl_104, pm_125, pm_126, pm_127, pm_128, \
                         pm_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_y[k] * pl_100[k]
                       + pm_125[k];

            t_191[k] = -ab_y[k] * pl_101[k]
                       + pm_126[k];

            t_192[k] = -ab_y[k] * pl_102[k]
                       + pm_127[k];

            t_193[k] = -ab_y[k] * pl_103[k]
                       + pm_128[k];

            t_194[k] = -ab_y[k] * pl_104[k]
                       + pm_129[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_y, pl_105, pl_106, pl_107, \
                         pl_108, pl_109, pm_131, pm_132, pm_133, pm_134, \
                         pm_135 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_y[k] * pl_105[k]
                       + pm_131[k];

            t_196[k] = -ab_y[k] * pl_106[k]
                       + pm_132[k];

            t_197[k] = -ab_y[k] * pl_107[k]
                       + pm_133[k];

            t_198[k] = -ab_y[k] * pl_108[k]
                       + pm_134[k];

            t_199[k] = -ab_y[k] * pl_109[k]
                       + pm_135[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_y, pl_110, pl_111, pl_112, \
                         pl_113, pl_114, pm_136, pm_138, pm_139, pm_140, \
                         pm_141 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_y[k] * pl_110[k]
                       + pm_136[k];

            t_201[k] = -ab_y[k] * pl_111[k]
                       + pm_138[k];

            t_202[k] = -ab_y[k] * pl_112[k]
                       + pm_139[k];

            t_203[k] = -ab_y[k] * pl_113[k]
                       + pm_140[k];

            t_204[k] = -ab_y[k] * pl_114[k]
                       + pm_141[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, pl_115, pl_116, pl_117, \
                         pl_118, pl_119, pm_142, pm_143, pm_144, pm_146, \
                         pm_147 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_y[k] * pl_115[k]
                       + pm_142[k];

            t_206[k] = -ab_y[k] * pl_116[k]
                       + pm_143[k];

            t_207[k] = -ab_y[k] * pl_117[k]
                       + pm_144[k];

            t_208[k] = -ab_y[k] * pl_118[k]
                       + pm_146[k];

            t_209[k] = -ab_y[k] * pl_119[k]
                       + pm_147[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_y, pl_120, pl_121, pl_122, \
                         pl_123, pl_124, pm_148, pm_149, pm_150, pm_151, \
                         pm_152 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_y[k] * pl_120[k]
                       + pm_148[k];

            t_211[k] = -ab_y[k] * pl_121[k]
                       + pm_149[k];

            t_212[k] = -ab_y[k] * pl_122[k]
                       + pm_150[k];

            t_213[k] = -ab_y[k] * pl_123[k]
                       + pm_151[k];

            t_214[k] = -ab_y[k] * pl_124[k]
                       + pm_152[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_y, pl_125, pl_126, pl_127, \
                         pl_128, pl_129, pm_153, pm_155, pm_156, pm_157, \
                         pm_158 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_y[k] * pl_125[k]
                       + pm_153[k];

            t_216[k] = -ab_y[k] * pl_126[k]
                       + pm_155[k];

            t_217[k] = -ab_y[k] * pl_127[k]
                       + pm_156[k];

            t_218[k] = -ab_y[k] * pl_128[k]
                       + pm_157[k];

            t_219[k] = -ab_y[k] * pl_129[k]
                       + pm_158[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, pl_130, pl_131, pl_132, \
                         pl_133, pl_134, pm_159, pm_160, pm_161, pm_162, \
                         pm_163 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_y[k] * pl_130[k]
                       + pm_159[k];

            t_221[k] = -ab_y[k] * pl_131[k]
                       + pm_160[k];

            t_222[k] = -ab_y[k] * pl_132[k]
                       + pm_161[k];

            t_223[k] = -ab_y[k] * pl_133[k]
                       + pm_162[k];

            t_224[k] = -ab_y[k] * pl_134[k]
                       + pm_163[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_z, pl_90, pl_91, pl_92, pl_93, \
                         pl_94, pm_112, pm_114, pm_115, pm_117, \
                         pm_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_z[k] * pl_90[k]
                       + pm_112[k];

            t_226[k] = -ab_z[k] * pl_91[k]
                       + pm_114[k];

            t_227[k] = -ab_z[k] * pl_92[k]
                       + pm_115[k];

            t_228[k] = -ab_z[k] * pl_93[k]
                       + pm_117[k];

            t_229[k] = -ab_z[k] * pl_94[k]
                       + pm_118[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_z, pl_95, pl_96, pl_97, pl_98, \
                         pl_99, pm_119, pm_121, pm_122, pm_123, \
                         pm_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_z[k] * pl_95[k]
                       + pm_119[k];

            t_231[k] = -ab_z[k] * pl_96[k]
                       + pm_121[k];

            t_232[k] = -ab_z[k] * pl_97[k]
                       + pm_122[k];

            t_233[k] = -ab_z[k] * pl_98[k]
                       + pm_123[k];

            t_234[k] = -ab_z[k] * pl_99[k]
                       + pm_124[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_z, pl_100, pl_101, pl_102, \
                         pl_103, pl_104, pm_126, pm_127, pm_128, pm_129, \
                         pm_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_z[k] * pl_100[k]
                       + pm_126[k];

            t_236[k] = -ab_z[k] * pl_101[k]
                       + pm_127[k];

            t_237[k] = -ab_z[k] * pl_102[k]
                       + pm_128[k];

            t_238[k] = -ab_z[k] * pl_103[k]
                       + pm_129[k];

            t_239[k] = -ab_z[k] * pl_104[k]
                       + pm_130[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_z, pl_105, pl_106, pl_107, \
                         pl_108, pl_109, pm_132, pm_133, pm_134, pm_135, \
                         pm_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_z[k] * pl_105[k]
                       + pm_132[k];

            t_241[k] = -ab_z[k] * pl_106[k]
                       + pm_133[k];

            t_242[k] = -ab_z[k] * pl_107[k]
                       + pm_134[k];

            t_243[k] = -ab_z[k] * pl_108[k]
                       + pm_135[k];

            t_244[k] = -ab_z[k] * pl_109[k]
                       + pm_136[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_z, pl_110, pl_111, pl_112, \
                         pl_113, pl_114, pm_137, pm_139, pm_140, pm_141, \
                         pm_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_z[k] * pl_110[k]
                       + pm_137[k];

            t_246[k] = -ab_z[k] * pl_111[k]
                       + pm_139[k];

            t_247[k] = -ab_z[k] * pl_112[k]
                       + pm_140[k];

            t_248[k] = -ab_z[k] * pl_113[k]
                       + pm_141[k];

            t_249[k] = -ab_z[k] * pl_114[k]
                       + pm_142[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_z, pl_115, pl_116, pl_117, \
                         pl_118, pl_119, pm_143, pm_144, pm_145, pm_147, \
                         pm_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_z[k] * pl_115[k]
                       + pm_143[k];

            t_251[k] = -ab_z[k] * pl_116[k]
                       + pm_144[k];

            t_252[k] = -ab_z[k] * pl_117[k]
                       + pm_145[k];

            t_253[k] = -ab_z[k] * pl_118[k]
                       + pm_147[k];

            t_254[k] = -ab_z[k] * pl_119[k]
                       + pm_148[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_z, pl_120, pl_121, pl_122, \
                         pl_123, pl_124, pm_149, pm_150, pm_151, pm_152, \
                         pm_153 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_z[k] * pl_120[k]
                       + pm_149[k];

            t_256[k] = -ab_z[k] * pl_121[k]
                       + pm_150[k];

            t_257[k] = -ab_z[k] * pl_122[k]
                       + pm_151[k];

            t_258[k] = -ab_z[k] * pl_123[k]
                       + pm_152[k];

            t_259[k] = -ab_z[k] * pl_124[k]
                       + pm_153[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_z, pl_125, pl_126, pl_127, \
                         pl_128, pl_129, pm_154, pm_156, pm_157, pm_158, \
                         pm_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_z[k] * pl_125[k]
                       + pm_154[k];

            t_261[k] = -ab_z[k] * pl_126[k]
                       + pm_156[k];

            t_262[k] = -ab_z[k] * pl_127[k]
                       + pm_157[k];

            t_263[k] = -ab_z[k] * pl_128[k]
                       + pm_158[k];

            t_264[k] = -ab_z[k] * pl_129[k]
                       + pm_159[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_z, pl_130, pl_131, pl_132, \
                         pl_133, pl_134, pm_160, pm_161, pm_162, pm_163, \
                         pm_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_z[k] * pl_130[k]
                       + pm_160[k];

            t_266[k] = -ab_z[k] * pl_131[k]
                       + pm_161[k];

            t_267[k] = -ab_z[k] * pl_132[k]
                       + pm_162[k];

            t_268[k] = -ab_z[k] * pl_133[k]
                       + pm_163[k];

            t_269[k] = -ab_z[k] * pl_134[k]
                       + pm_164[k];
        }
    }
}

auto
compute_hrr_dl(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pl, const size_t pm, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_dl_piece0(buffer, coordinates, target, pl, pm, ncomps, nmax);

    compute_hrr_dl_piece1(buffer, coordinates, target, pl, pm, ncomps, nmax);
}

}  // namespace simdtrf
