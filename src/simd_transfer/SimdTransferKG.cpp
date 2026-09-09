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


#include "SimdTransferKG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_kg_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kf, const size_t lf,
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

        const auto *kf_0 = buffer.data(kf + 0 * ncomps + c);
        const auto *kf_1 = buffer.data(kf + 1 * ncomps + c);
        const auto *kf_2 = buffer.data(kf + 2 * ncomps + c);
        const auto *kf_3 = buffer.data(kf + 3 * ncomps + c);
        const auto *kf_4 = buffer.data(kf + 4 * ncomps + c);
        const auto *kf_5 = buffer.data(kf + 5 * ncomps + c);
        const auto *kf_6 = buffer.data(kf + 6 * ncomps + c);
        const auto *kf_7 = buffer.data(kf + 7 * ncomps + c);
        const auto *kf_8 = buffer.data(kf + 8 * ncomps + c);
        const auto *kf_9 = buffer.data(kf + 9 * ncomps + c);
        const auto *kf_10 = buffer.data(kf + 10 * ncomps + c);
        const auto *kf_11 = buffer.data(kf + 11 * ncomps + c);
        const auto *kf_12 = buffer.data(kf + 12 * ncomps + c);
        const auto *kf_13 = buffer.data(kf + 13 * ncomps + c);
        const auto *kf_14 = buffer.data(kf + 14 * ncomps + c);
        const auto *kf_15 = buffer.data(kf + 15 * ncomps + c);
        const auto *kf_16 = buffer.data(kf + 16 * ncomps + c);
        const auto *kf_17 = buffer.data(kf + 17 * ncomps + c);
        const auto *kf_18 = buffer.data(kf + 18 * ncomps + c);
        const auto *kf_19 = buffer.data(kf + 19 * ncomps + c);
        const auto *kf_20 = buffer.data(kf + 20 * ncomps + c);
        const auto *kf_21 = buffer.data(kf + 21 * ncomps + c);
        const auto *kf_22 = buffer.data(kf + 22 * ncomps + c);
        const auto *kf_23 = buffer.data(kf + 23 * ncomps + c);
        const auto *kf_24 = buffer.data(kf + 24 * ncomps + c);
        const auto *kf_25 = buffer.data(kf + 25 * ncomps + c);
        const auto *kf_26 = buffer.data(kf + 26 * ncomps + c);
        const auto *kf_27 = buffer.data(kf + 27 * ncomps + c);
        const auto *kf_28 = buffer.data(kf + 28 * ncomps + c);
        const auto *kf_29 = buffer.data(kf + 29 * ncomps + c);
        const auto *kf_30 = buffer.data(kf + 30 * ncomps + c);
        const auto *kf_31 = buffer.data(kf + 31 * ncomps + c);
        const auto *kf_32 = buffer.data(kf + 32 * ncomps + c);
        const auto *kf_33 = buffer.data(kf + 33 * ncomps + c);
        const auto *kf_34 = buffer.data(kf + 34 * ncomps + c);
        const auto *kf_35 = buffer.data(kf + 35 * ncomps + c);
        const auto *kf_36 = buffer.data(kf + 36 * ncomps + c);
        const auto *kf_37 = buffer.data(kf + 37 * ncomps + c);
        const auto *kf_38 = buffer.data(kf + 38 * ncomps + c);
        const auto *kf_39 = buffer.data(kf + 39 * ncomps + c);
        const auto *kf_40 = buffer.data(kf + 40 * ncomps + c);
        const auto *kf_41 = buffer.data(kf + 41 * ncomps + c);
        const auto *kf_42 = buffer.data(kf + 42 * ncomps + c);
        const auto *kf_43 = buffer.data(kf + 43 * ncomps + c);
        const auto *kf_44 = buffer.data(kf + 44 * ncomps + c);
        const auto *kf_45 = buffer.data(kf + 45 * ncomps + c);
        const auto *kf_46 = buffer.data(kf + 46 * ncomps + c);
        const auto *kf_47 = buffer.data(kf + 47 * ncomps + c);
        const auto *kf_48 = buffer.data(kf + 48 * ncomps + c);
        const auto *kf_49 = buffer.data(kf + 49 * ncomps + c);
        const auto *kf_50 = buffer.data(kf + 50 * ncomps + c);
        const auto *kf_51 = buffer.data(kf + 51 * ncomps + c);
        const auto *kf_52 = buffer.data(kf + 52 * ncomps + c);
        const auto *kf_53 = buffer.data(kf + 53 * ncomps + c);
        const auto *kf_54 = buffer.data(kf + 54 * ncomps + c);
        const auto *kf_55 = buffer.data(kf + 55 * ncomps + c);
        const auto *kf_56 = buffer.data(kf + 56 * ncomps + c);
        const auto *kf_57 = buffer.data(kf + 57 * ncomps + c);
        const auto *kf_58 = buffer.data(kf + 58 * ncomps + c);
        const auto *kf_59 = buffer.data(kf + 59 * ncomps + c);
        const auto *kf_60 = buffer.data(kf + 60 * ncomps + c);
        const auto *kf_61 = buffer.data(kf + 61 * ncomps + c);
        const auto *kf_62 = buffer.data(kf + 62 * ncomps + c);
        const auto *kf_63 = buffer.data(kf + 63 * ncomps + c);
        const auto *kf_64 = buffer.data(kf + 64 * ncomps + c);
        const auto *kf_65 = buffer.data(kf + 65 * ncomps + c);
        const auto *kf_66 = buffer.data(kf + 66 * ncomps + c);
        const auto *kf_67 = buffer.data(kf + 67 * ncomps + c);
        const auto *kf_68 = buffer.data(kf + 68 * ncomps + c);
        const auto *kf_69 = buffer.data(kf + 69 * ncomps + c);
        const auto *kf_70 = buffer.data(kf + 70 * ncomps + c);
        const auto *kf_71 = buffer.data(kf + 71 * ncomps + c);
        const auto *kf_72 = buffer.data(kf + 72 * ncomps + c);
        const auto *kf_73 = buffer.data(kf + 73 * ncomps + c);
        const auto *kf_74 = buffer.data(kf + 74 * ncomps + c);
        const auto *kf_75 = buffer.data(kf + 75 * ncomps + c);
        const auto *kf_76 = buffer.data(kf + 76 * ncomps + c);
        const auto *kf_77 = buffer.data(kf + 77 * ncomps + c);
        const auto *kf_78 = buffer.data(kf + 78 * ncomps + c);
        const auto *kf_79 = buffer.data(kf + 79 * ncomps + c);
        const auto *kf_80 = buffer.data(kf + 80 * ncomps + c);
        const auto *kf_81 = buffer.data(kf + 81 * ncomps + c);
        const auto *kf_82 = buffer.data(kf + 82 * ncomps + c);
        const auto *kf_83 = buffer.data(kf + 83 * ncomps + c);
        const auto *kf_84 = buffer.data(kf + 84 * ncomps + c);
        const auto *kf_85 = buffer.data(kf + 85 * ncomps + c);
        const auto *kf_86 = buffer.data(kf + 86 * ncomps + c);
        const auto *kf_87 = buffer.data(kf + 87 * ncomps + c);
        const auto *kf_88 = buffer.data(kf + 88 * ncomps + c);
        const auto *kf_89 = buffer.data(kf + 89 * ncomps + c);
        const auto *kf_90 = buffer.data(kf + 90 * ncomps + c);
        const auto *kf_91 = buffer.data(kf + 91 * ncomps + c);
        const auto *kf_92 = buffer.data(kf + 92 * ncomps + c);
        const auto *kf_93 = buffer.data(kf + 93 * ncomps + c);
        const auto *kf_94 = buffer.data(kf + 94 * ncomps + c);
        const auto *kf_95 = buffer.data(kf + 95 * ncomps + c);
        const auto *kf_96 = buffer.data(kf + 96 * ncomps + c);
        const auto *kf_97 = buffer.data(kf + 97 * ncomps + c);
        const auto *kf_98 = buffer.data(kf + 98 * ncomps + c);
        const auto *kf_99 = buffer.data(kf + 99 * ncomps + c);

        const auto *lf_0 = buffer.data(lf + 0 * ncomps + c);
        const auto *lf_1 = buffer.data(lf + 1 * ncomps + c);
        const auto *lf_2 = buffer.data(lf + 2 * ncomps + c);
        const auto *lf_3 = buffer.data(lf + 3 * ncomps + c);
        const auto *lf_4 = buffer.data(lf + 4 * ncomps + c);
        const auto *lf_5 = buffer.data(lf + 5 * ncomps + c);
        const auto *lf_6 = buffer.data(lf + 6 * ncomps + c);
        const auto *lf_7 = buffer.data(lf + 7 * ncomps + c);
        const auto *lf_8 = buffer.data(lf + 8 * ncomps + c);
        const auto *lf_9 = buffer.data(lf + 9 * ncomps + c);
        const auto *lf_10 = buffer.data(lf + 10 * ncomps + c);
        const auto *lf_11 = buffer.data(lf + 11 * ncomps + c);
        const auto *lf_12 = buffer.data(lf + 12 * ncomps + c);
        const auto *lf_13 = buffer.data(lf + 13 * ncomps + c);
        const auto *lf_14 = buffer.data(lf + 14 * ncomps + c);
        const auto *lf_15 = buffer.data(lf + 15 * ncomps + c);
        const auto *lf_16 = buffer.data(lf + 16 * ncomps + c);
        const auto *lf_17 = buffer.data(lf + 17 * ncomps + c);
        const auto *lf_18 = buffer.data(lf + 18 * ncomps + c);
        const auto *lf_19 = buffer.data(lf + 19 * ncomps + c);
        const auto *lf_20 = buffer.data(lf + 20 * ncomps + c);
        const auto *lf_21 = buffer.data(lf + 21 * ncomps + c);
        const auto *lf_22 = buffer.data(lf + 22 * ncomps + c);
        const auto *lf_23 = buffer.data(lf + 23 * ncomps + c);
        const auto *lf_24 = buffer.data(lf + 24 * ncomps + c);
        const auto *lf_25 = buffer.data(lf + 25 * ncomps + c);
        const auto *lf_26 = buffer.data(lf + 26 * ncomps + c);
        const auto *lf_27 = buffer.data(lf + 27 * ncomps + c);
        const auto *lf_28 = buffer.data(lf + 28 * ncomps + c);
        const auto *lf_29 = buffer.data(lf + 29 * ncomps + c);
        const auto *lf_30 = buffer.data(lf + 30 * ncomps + c);
        const auto *lf_31 = buffer.data(lf + 31 * ncomps + c);
        const auto *lf_32 = buffer.data(lf + 32 * ncomps + c);
        const auto *lf_33 = buffer.data(lf + 33 * ncomps + c);
        const auto *lf_34 = buffer.data(lf + 34 * ncomps + c);
        const auto *lf_35 = buffer.data(lf + 35 * ncomps + c);
        const auto *lf_36 = buffer.data(lf + 36 * ncomps + c);
        const auto *lf_37 = buffer.data(lf + 37 * ncomps + c);
        const auto *lf_38 = buffer.data(lf + 38 * ncomps + c);
        const auto *lf_39 = buffer.data(lf + 39 * ncomps + c);
        const auto *lf_40 = buffer.data(lf + 40 * ncomps + c);
        const auto *lf_41 = buffer.data(lf + 41 * ncomps + c);
        const auto *lf_42 = buffer.data(lf + 42 * ncomps + c);
        const auto *lf_43 = buffer.data(lf + 43 * ncomps + c);
        const auto *lf_44 = buffer.data(lf + 44 * ncomps + c);
        const auto *lf_45 = buffer.data(lf + 45 * ncomps + c);
        const auto *lf_46 = buffer.data(lf + 46 * ncomps + c);
        const auto *lf_47 = buffer.data(lf + 47 * ncomps + c);
        const auto *lf_48 = buffer.data(lf + 48 * ncomps + c);
        const auto *lf_49 = buffer.data(lf + 49 * ncomps + c);
        const auto *lf_50 = buffer.data(lf + 50 * ncomps + c);
        const auto *lf_51 = buffer.data(lf + 51 * ncomps + c);
        const auto *lf_52 = buffer.data(lf + 52 * ncomps + c);
        const auto *lf_53 = buffer.data(lf + 53 * ncomps + c);
        const auto *lf_54 = buffer.data(lf + 54 * ncomps + c);
        const auto *lf_55 = buffer.data(lf + 55 * ncomps + c);
        const auto *lf_56 = buffer.data(lf + 56 * ncomps + c);
        const auto *lf_57 = buffer.data(lf + 57 * ncomps + c);
        const auto *lf_58 = buffer.data(lf + 58 * ncomps + c);
        const auto *lf_59 = buffer.data(lf + 59 * ncomps + c);
        const auto *lf_60 = buffer.data(lf + 60 * ncomps + c);
        const auto *lf_61 = buffer.data(lf + 61 * ncomps + c);
        const auto *lf_62 = buffer.data(lf + 62 * ncomps + c);
        const auto *lf_63 = buffer.data(lf + 63 * ncomps + c);
        const auto *lf_64 = buffer.data(lf + 64 * ncomps + c);
        const auto *lf_65 = buffer.data(lf + 65 * ncomps + c);
        const auto *lf_66 = buffer.data(lf + 66 * ncomps + c);
        const auto *lf_67 = buffer.data(lf + 67 * ncomps + c);
        const auto *lf_68 = buffer.data(lf + 68 * ncomps + c);
        const auto *lf_69 = buffer.data(lf + 69 * ncomps + c);
        const auto *lf_70 = buffer.data(lf + 70 * ncomps + c);
        const auto *lf_71 = buffer.data(lf + 71 * ncomps + c);
        const auto *lf_72 = buffer.data(lf + 72 * ncomps + c);
        const auto *lf_73 = buffer.data(lf + 73 * ncomps + c);
        const auto *lf_74 = buffer.data(lf + 74 * ncomps + c);
        const auto *lf_75 = buffer.data(lf + 75 * ncomps + c);
        const auto *lf_76 = buffer.data(lf + 76 * ncomps + c);
        const auto *lf_77 = buffer.data(lf + 77 * ncomps + c);
        const auto *lf_78 = buffer.data(lf + 78 * ncomps + c);
        const auto *lf_79 = buffer.data(lf + 79 * ncomps + c);
        const auto *lf_80 = buffer.data(lf + 80 * ncomps + c);
        const auto *lf_81 = buffer.data(lf + 81 * ncomps + c);
        const auto *lf_82 = buffer.data(lf + 82 * ncomps + c);
        const auto *lf_83 = buffer.data(lf + 83 * ncomps + c);
        const auto *lf_84 = buffer.data(lf + 84 * ncomps + c);
        const auto *lf_85 = buffer.data(lf + 85 * ncomps + c);
        const auto *lf_86 = buffer.data(lf + 86 * ncomps + c);
        const auto *lf_87 = buffer.data(lf + 87 * ncomps + c);
        const auto *lf_88 = buffer.data(lf + 88 * ncomps + c);
        const auto *lf_89 = buffer.data(lf + 89 * ncomps + c);
        const auto *lf_90 = buffer.data(lf + 90 * ncomps + c);
        const auto *lf_91 = buffer.data(lf + 91 * ncomps + c);
        const auto *lf_92 = buffer.data(lf + 92 * ncomps + c);
        const auto *lf_93 = buffer.data(lf + 93 * ncomps + c);
        const auto *lf_94 = buffer.data(lf + 94 * ncomps + c);
        const auto *lf_95 = buffer.data(lf + 95 * ncomps + c);
        const auto *lf_96 = buffer.data(lf + 96 * ncomps + c);
        const auto *lf_97 = buffer.data(lf + 97 * ncomps + c);
        const auto *lf_98 = buffer.data(lf + 98 * ncomps + c);
        const auto *lf_99 = buffer.data(lf + 99 * ncomps + c);
        const auto *lf_106 = buffer.data(lf + 106 * ncomps + c);
        const auto *lf_107 = buffer.data(lf + 107 * ncomps + c);
        const auto *lf_108 = buffer.data(lf + 108 * ncomps + c);
        const auto *lf_109 = buffer.data(lf + 109 * ncomps + c);
        const auto *lf_116 = buffer.data(lf + 116 * ncomps + c);
        const auto *lf_117 = buffer.data(lf + 117 * ncomps + c);
        const auto *lf_118 = buffer.data(lf + 118 * ncomps + c);
        const auto *lf_119 = buffer.data(lf + 119 * ncomps + c);
        const auto *lf_126 = buffer.data(lf + 126 * ncomps + c);
        const auto *lf_127 = buffer.data(lf + 127 * ncomps + c);
        const auto *lf_128 = buffer.data(lf + 128 * ncomps + c);
        const auto *lf_129 = buffer.data(lf + 129 * ncomps + c);
        const auto *lf_139 = buffer.data(lf + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, kf_0, kf_1, kf_2, kf_3, kf_4, lf_0, \
                         lf_1, lf_2, lf_3, lf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * kf_0[k]
                     + lf_0[k];

            t_1[k] = ab_x[k] * kf_1[k]
                     + lf_1[k];

            t_2[k] = ab_x[k] * kf_2[k]
                     + lf_2[k];

            t_3[k] = ab_x[k] * kf_3[k]
                     + lf_3[k];

            t_4[k] = ab_x[k] * kf_4[k]
                     + lf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, kf_5, kf_6, kf_7, kf_8, kf_9, lf_5, \
                         lf_6, lf_7, lf_8, lf_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * kf_5[k]
                     + lf_5[k];

            t_6[k] = ab_x[k] * kf_6[k]
                     + lf_6[k];

            t_7[k] = ab_x[k] * kf_7[k]
                     + lf_7[k];

            t_8[k] = ab_x[k] * kf_8[k]
                     + lf_8[k];

            t_9[k] = ab_x[k] * kf_9[k]
                     + lf_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, kf_6, kf_7, kf_8, kf_9, \
                         lf_16, lf_17, lf_18, lf_19, lf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * kf_6[k]
                      + lf_16[k];

            t_11[k] = ab_y[k] * kf_7[k]
                      + lf_17[k];

            t_12[k] = ab_y[k] * kf_8[k]
                      + lf_18[k];

            t_13[k] = ab_y[k] * kf_9[k]
                      + lf_19[k];

            t_14[k] = ab_z[k] * kf_9[k]
                      + lf_29[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, kf_10, kf_11, kf_12, kf_13, \
                         kf_14, lf_10, lf_11, lf_12, lf_13, lf_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * kf_10[k]
                      + lf_10[k];

            t_16[k] = ab_x[k] * kf_11[k]
                      + lf_11[k];

            t_17[k] = ab_x[k] * kf_12[k]
                      + lf_12[k];

            t_18[k] = ab_x[k] * kf_13[k]
                      + lf_13[k];

            t_19[k] = ab_x[k] * kf_14[k]
                      + lf_14[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, kf_15, kf_16, kf_17, kf_18, \
                         kf_19, lf_15, lf_16, lf_17, lf_18, lf_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * kf_15[k]
                      + lf_15[k];

            t_21[k] = ab_x[k] * kf_16[k]
                      + lf_16[k];

            t_22[k] = ab_x[k] * kf_17[k]
                      + lf_17[k];

            t_23[k] = ab_x[k] * kf_18[k]
                      + lf_18[k];

            t_24[k] = ab_x[k] * kf_19[k]
                      + lf_19[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, kf_16, kf_17, kf_18, kf_19, \
                         lf_36, lf_37, lf_38, lf_39, lf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_y[k] * kf_16[k]
                      + lf_36[k];

            t_26[k] = ab_y[k] * kf_17[k]
                      + lf_37[k];

            t_27[k] = ab_y[k] * kf_18[k]
                      + lf_38[k];

            t_28[k] = ab_y[k] * kf_19[k]
                      + lf_39[k];

            t_29[k] = ab_z[k] * kf_19[k]
                      + lf_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, kf_20, kf_21, kf_22, kf_23, \
                         kf_24, lf_20, lf_21, lf_22, lf_23, lf_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * kf_20[k]
                      + lf_20[k];

            t_31[k] = ab_x[k] * kf_21[k]
                      + lf_21[k];

            t_32[k] = ab_x[k] * kf_22[k]
                      + lf_22[k];

            t_33[k] = ab_x[k] * kf_23[k]
                      + lf_23[k];

            t_34[k] = ab_x[k] * kf_24[k]
                      + lf_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, kf_25, kf_26, kf_27, kf_28, \
                         kf_29, lf_25, lf_26, lf_27, lf_28, lf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * kf_25[k]
                      + lf_25[k];

            t_36[k] = ab_x[k] * kf_26[k]
                      + lf_26[k];

            t_37[k] = ab_x[k] * kf_27[k]
                      + lf_27[k];

            t_38[k] = ab_x[k] * kf_28[k]
                      + lf_28[k];

            t_39[k] = ab_x[k] * kf_29[k]
                      + lf_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, kf_26, kf_27, kf_28, kf_29, \
                         lf_46, lf_47, lf_48, lf_49, lf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * kf_26[k]
                      + lf_46[k];

            t_41[k] = ab_y[k] * kf_27[k]
                      + lf_47[k];

            t_42[k] = ab_y[k] * kf_28[k]
                      + lf_48[k];

            t_43[k] = ab_y[k] * kf_29[k]
                      + lf_49[k];

            t_44[k] = ab_z[k] * kf_29[k]
                      + lf_59[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, kf_30, kf_31, kf_32, kf_33, \
                         kf_34, lf_30, lf_31, lf_32, lf_33, lf_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * kf_30[k]
                      + lf_30[k];

            t_46[k] = ab_x[k] * kf_31[k]
                      + lf_31[k];

            t_47[k] = ab_x[k] * kf_32[k]
                      + lf_32[k];

            t_48[k] = ab_x[k] * kf_33[k]
                      + lf_33[k];

            t_49[k] = ab_x[k] * kf_34[k]
                      + lf_34[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, kf_35, kf_36, kf_37, kf_38, \
                         kf_39, lf_35, lf_36, lf_37, lf_38, lf_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * kf_35[k]
                      + lf_35[k];

            t_51[k] = ab_x[k] * kf_36[k]
                      + lf_36[k];

            t_52[k] = ab_x[k] * kf_37[k]
                      + lf_37[k];

            t_53[k] = ab_x[k] * kf_38[k]
                      + lf_38[k];

            t_54[k] = ab_x[k] * kf_39[k]
                      + lf_39[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, kf_36, kf_37, kf_38, kf_39, \
                         lf_66, lf_67, lf_68, lf_69, lf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_y[k] * kf_36[k]
                      + lf_66[k];

            t_56[k] = ab_y[k] * kf_37[k]
                      + lf_67[k];

            t_57[k] = ab_y[k] * kf_38[k]
                      + lf_68[k];

            t_58[k] = ab_y[k] * kf_39[k]
                      + lf_69[k];

            t_59[k] = ab_z[k] * kf_39[k]
                      + lf_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, kf_40, kf_41, kf_42, kf_43, \
                         kf_44, lf_40, lf_41, lf_42, lf_43, lf_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * kf_40[k]
                      + lf_40[k];

            t_61[k] = ab_x[k] * kf_41[k]
                      + lf_41[k];

            t_62[k] = ab_x[k] * kf_42[k]
                      + lf_42[k];

            t_63[k] = ab_x[k] * kf_43[k]
                      + lf_43[k];

            t_64[k] = ab_x[k] * kf_44[k]
                      + lf_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, kf_45, kf_46, kf_47, kf_48, \
                         kf_49, lf_45, lf_46, lf_47, lf_48, lf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * kf_45[k]
                      + lf_45[k];

            t_66[k] = ab_x[k] * kf_46[k]
                      + lf_46[k];

            t_67[k] = ab_x[k] * kf_47[k]
                      + lf_47[k];

            t_68[k] = ab_x[k] * kf_48[k]
                      + lf_48[k];

            t_69[k] = ab_x[k] * kf_49[k]
                      + lf_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, kf_46, kf_47, kf_48, kf_49, \
                         lf_76, lf_77, lf_78, lf_79, lf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * kf_46[k]
                      + lf_76[k];

            t_71[k] = ab_y[k] * kf_47[k]
                      + lf_77[k];

            t_72[k] = ab_y[k] * kf_48[k]
                      + lf_78[k];

            t_73[k] = ab_y[k] * kf_49[k]
                      + lf_79[k];

            t_74[k] = ab_z[k] * kf_49[k]
                      + lf_89[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, kf_50, kf_51, kf_52, kf_53, \
                         kf_54, lf_50, lf_51, lf_52, lf_53, lf_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * kf_50[k]
                      + lf_50[k];

            t_76[k] = ab_x[k] * kf_51[k]
                      + lf_51[k];

            t_77[k] = ab_x[k] * kf_52[k]
                      + lf_52[k];

            t_78[k] = ab_x[k] * kf_53[k]
                      + lf_53[k];

            t_79[k] = ab_x[k] * kf_54[k]
                      + lf_54[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, kf_55, kf_56, kf_57, kf_58, \
                         kf_59, lf_55, lf_56, lf_57, lf_58, lf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * kf_55[k]
                      + lf_55[k];

            t_81[k] = ab_x[k] * kf_56[k]
                      + lf_56[k];

            t_82[k] = ab_x[k] * kf_57[k]
                      + lf_57[k];

            t_83[k] = ab_x[k] * kf_58[k]
                      + lf_58[k];

            t_84[k] = ab_x[k] * kf_59[k]
                      + lf_59[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, kf_56, kf_57, kf_58, kf_59, \
                         lf_86, lf_87, lf_88, lf_89, lf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_y[k] * kf_56[k]
                      + lf_86[k];

            t_86[k] = ab_y[k] * kf_57[k]
                      + lf_87[k];

            t_87[k] = ab_y[k] * kf_58[k]
                      + lf_88[k];

            t_88[k] = ab_y[k] * kf_59[k]
                      + lf_89[k];

            t_89[k] = ab_z[k] * kf_59[k]
                      + lf_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, kf_60, kf_61, kf_62, kf_63, \
                         kf_64, lf_60, lf_61, lf_62, lf_63, lf_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * kf_60[k]
                      + lf_60[k];

            t_91[k] = ab_x[k] * kf_61[k]
                      + lf_61[k];

            t_92[k] = ab_x[k] * kf_62[k]
                      + lf_62[k];

            t_93[k] = ab_x[k] * kf_63[k]
                      + lf_63[k];

            t_94[k] = ab_x[k] * kf_64[k]
                      + lf_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, kf_65, kf_66, kf_67, kf_68, \
                         kf_69, lf_65, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * kf_65[k]
                      + lf_65[k];

            t_96[k] = ab_x[k] * kf_66[k]
                      + lf_66[k];

            t_97[k] = ab_x[k] * kf_67[k]
                      + lf_67[k];

            t_98[k] = ab_x[k] * kf_68[k]
                      + lf_68[k];

            t_99[k] = ab_x[k] * kf_69[k]
                      + lf_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, kf_66, kf_67, kf_68, \
                         kf_69, lf_106, lf_107, lf_108, lf_109, \
                         lf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * kf_66[k]
                       + lf_106[k];

            t_101[k] = ab_y[k] * kf_67[k]
                       + lf_107[k];

            t_102[k] = ab_y[k] * kf_68[k]
                       + lf_108[k];

            t_103[k] = ab_y[k] * kf_69[k]
                       + lf_109[k];

            t_104[k] = ab_z[k] * kf_69[k]
                       + lf_119[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, kf_70, kf_71, kf_72, kf_73, \
                         kf_74, lf_70, lf_71, lf_72, lf_73, lf_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * kf_70[k]
                       + lf_70[k];

            t_106[k] = ab_x[k] * kf_71[k]
                       + lf_71[k];

            t_107[k] = ab_x[k] * kf_72[k]
                       + lf_72[k];

            t_108[k] = ab_x[k] * kf_73[k]
                       + lf_73[k];

            t_109[k] = ab_x[k] * kf_74[k]
                       + lf_74[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, kf_75, kf_76, kf_77, kf_78, \
                         kf_79, lf_75, lf_76, lf_77, lf_78, lf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * kf_75[k]
                       + lf_75[k];

            t_111[k] = ab_x[k] * kf_76[k]
                       + lf_76[k];

            t_112[k] = ab_x[k] * kf_77[k]
                       + lf_77[k];

            t_113[k] = ab_x[k] * kf_78[k]
                       + lf_78[k];

            t_114[k] = ab_x[k] * kf_79[k]
                       + lf_79[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, kf_76, kf_77, kf_78, \
                         kf_79, lf_116, lf_117, lf_118, lf_119, \
                         lf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_y[k] * kf_76[k]
                       + lf_116[k];

            t_116[k] = ab_y[k] * kf_77[k]
                       + lf_117[k];

            t_117[k] = ab_y[k] * kf_78[k]
                       + lf_118[k];

            t_118[k] = ab_y[k] * kf_79[k]
                       + lf_119[k];

            t_119[k] = ab_z[k] * kf_79[k]
                       + lf_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, kf_80, kf_81, kf_82, kf_83, \
                         kf_84, lf_80, lf_81, lf_82, lf_83, lf_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * kf_80[k]
                       + lf_80[k];

            t_121[k] = ab_x[k] * kf_81[k]
                       + lf_81[k];

            t_122[k] = ab_x[k] * kf_82[k]
                       + lf_82[k];

            t_123[k] = ab_x[k] * kf_83[k]
                       + lf_83[k];

            t_124[k] = ab_x[k] * kf_84[k]
                       + lf_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, kf_85, kf_86, kf_87, kf_88, \
                         kf_89, lf_85, lf_86, lf_87, lf_88, lf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * kf_85[k]
                       + lf_85[k];

            t_126[k] = ab_x[k] * kf_86[k]
                       + lf_86[k];

            t_127[k] = ab_x[k] * kf_87[k]
                       + lf_87[k];

            t_128[k] = ab_x[k] * kf_88[k]
                       + lf_88[k];

            t_129[k] = ab_x[k] * kf_89[k]
                       + lf_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, kf_86, kf_87, kf_88, \
                         kf_89, lf_126, lf_127, lf_128, lf_129, \
                         lf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * kf_86[k]
                       + lf_126[k];

            t_131[k] = ab_y[k] * kf_87[k]
                       + lf_127[k];

            t_132[k] = ab_y[k] * kf_88[k]
                       + lf_128[k];

            t_133[k] = ab_y[k] * kf_89[k]
                       + lf_129[k];

            t_134[k] = ab_z[k] * kf_89[k]
                       + lf_139[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, kf_90, kf_91, kf_92, kf_93, \
                         kf_94, lf_90, lf_91, lf_92, lf_93, lf_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * kf_90[k]
                       + lf_90[k];

            t_136[k] = ab_x[k] * kf_91[k]
                       + lf_91[k];

            t_137[k] = ab_x[k] * kf_92[k]
                       + lf_92[k];

            t_138[k] = ab_x[k] * kf_93[k]
                       + lf_93[k];

            t_139[k] = ab_x[k] * kf_94[k]
                       + lf_94[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, kf_95, kf_96, kf_97, kf_98, \
                         kf_99, lf_95, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * kf_95[k]
                       + lf_95[k];

            t_141[k] = ab_x[k] * kf_96[k]
                       + lf_96[k];

            t_142[k] = ab_x[k] * kf_97[k]
                       + lf_97[k];

            t_143[k] = ab_x[k] * kf_98[k]
                       + lf_98[k];

            t_144[k] = ab_x[k] * kf_99[k]
                       + lf_99[k];
        }
    }
}

static auto
compute_hrr_kg_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kf, const size_t lf,
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

        const auto *kf_96 = buffer.data(kf + 96 * ncomps + c);
        const auto *kf_97 = buffer.data(kf + 97 * ncomps + c);
        const auto *kf_98 = buffer.data(kf + 98 * ncomps + c);
        const auto *kf_99 = buffer.data(kf + 99 * ncomps + c);
        const auto *kf_100 = buffer.data(kf + 100 * ncomps + c);
        const auto *kf_101 = buffer.data(kf + 101 * ncomps + c);
        const auto *kf_102 = buffer.data(kf + 102 * ncomps + c);
        const auto *kf_103 = buffer.data(kf + 103 * ncomps + c);
        const auto *kf_104 = buffer.data(kf + 104 * ncomps + c);
        const auto *kf_105 = buffer.data(kf + 105 * ncomps + c);
        const auto *kf_106 = buffer.data(kf + 106 * ncomps + c);
        const auto *kf_107 = buffer.data(kf + 107 * ncomps + c);
        const auto *kf_108 = buffer.data(kf + 108 * ncomps + c);
        const auto *kf_109 = buffer.data(kf + 109 * ncomps + c);
        const auto *kf_110 = buffer.data(kf + 110 * ncomps + c);
        const auto *kf_111 = buffer.data(kf + 111 * ncomps + c);
        const auto *kf_112 = buffer.data(kf + 112 * ncomps + c);
        const auto *kf_113 = buffer.data(kf + 113 * ncomps + c);
        const auto *kf_114 = buffer.data(kf + 114 * ncomps + c);
        const auto *kf_115 = buffer.data(kf + 115 * ncomps + c);
        const auto *kf_116 = buffer.data(kf + 116 * ncomps + c);
        const auto *kf_117 = buffer.data(kf + 117 * ncomps + c);
        const auto *kf_118 = buffer.data(kf + 118 * ncomps + c);
        const auto *kf_119 = buffer.data(kf + 119 * ncomps + c);
        const auto *kf_120 = buffer.data(kf + 120 * ncomps + c);
        const auto *kf_121 = buffer.data(kf + 121 * ncomps + c);
        const auto *kf_122 = buffer.data(kf + 122 * ncomps + c);
        const auto *kf_123 = buffer.data(kf + 123 * ncomps + c);
        const auto *kf_124 = buffer.data(kf + 124 * ncomps + c);
        const auto *kf_125 = buffer.data(kf + 125 * ncomps + c);
        const auto *kf_126 = buffer.data(kf + 126 * ncomps + c);
        const auto *kf_127 = buffer.data(kf + 127 * ncomps + c);
        const auto *kf_128 = buffer.data(kf + 128 * ncomps + c);
        const auto *kf_129 = buffer.data(kf + 129 * ncomps + c);
        const auto *kf_130 = buffer.data(kf + 130 * ncomps + c);
        const auto *kf_131 = buffer.data(kf + 131 * ncomps + c);
        const auto *kf_132 = buffer.data(kf + 132 * ncomps + c);
        const auto *kf_133 = buffer.data(kf + 133 * ncomps + c);
        const auto *kf_134 = buffer.data(kf + 134 * ncomps + c);
        const auto *kf_135 = buffer.data(kf + 135 * ncomps + c);
        const auto *kf_136 = buffer.data(kf + 136 * ncomps + c);
        const auto *kf_137 = buffer.data(kf + 137 * ncomps + c);
        const auto *kf_138 = buffer.data(kf + 138 * ncomps + c);
        const auto *kf_139 = buffer.data(kf + 139 * ncomps + c);
        const auto *kf_140 = buffer.data(kf + 140 * ncomps + c);
        const auto *kf_141 = buffer.data(kf + 141 * ncomps + c);
        const auto *kf_142 = buffer.data(kf + 142 * ncomps + c);
        const auto *kf_143 = buffer.data(kf + 143 * ncomps + c);
        const auto *kf_144 = buffer.data(kf + 144 * ncomps + c);
        const auto *kf_145 = buffer.data(kf + 145 * ncomps + c);
        const auto *kf_146 = buffer.data(kf + 146 * ncomps + c);
        const auto *kf_147 = buffer.data(kf + 147 * ncomps + c);
        const auto *kf_148 = buffer.data(kf + 148 * ncomps + c);
        const auto *kf_149 = buffer.data(kf + 149 * ncomps + c);
        const auto *kf_150 = buffer.data(kf + 150 * ncomps + c);
        const auto *kf_151 = buffer.data(kf + 151 * ncomps + c);
        const auto *kf_152 = buffer.data(kf + 152 * ncomps + c);
        const auto *kf_153 = buffer.data(kf + 153 * ncomps + c);
        const auto *kf_154 = buffer.data(kf + 154 * ncomps + c);
        const auto *kf_155 = buffer.data(kf + 155 * ncomps + c);
        const auto *kf_156 = buffer.data(kf + 156 * ncomps + c);
        const auto *kf_157 = buffer.data(kf + 157 * ncomps + c);
        const auto *kf_158 = buffer.data(kf + 158 * ncomps + c);
        const auto *kf_159 = buffer.data(kf + 159 * ncomps + c);
        const auto *kf_160 = buffer.data(kf + 160 * ncomps + c);
        const auto *kf_161 = buffer.data(kf + 161 * ncomps + c);
        const auto *kf_162 = buffer.data(kf + 162 * ncomps + c);
        const auto *kf_163 = buffer.data(kf + 163 * ncomps + c);
        const auto *kf_164 = buffer.data(kf + 164 * ncomps + c);
        const auto *kf_165 = buffer.data(kf + 165 * ncomps + c);
        const auto *kf_166 = buffer.data(kf + 166 * ncomps + c);
        const auto *kf_167 = buffer.data(kf + 167 * ncomps + c);
        const auto *kf_168 = buffer.data(kf + 168 * ncomps + c);
        const auto *kf_169 = buffer.data(kf + 169 * ncomps + c);
        const auto *kf_170 = buffer.data(kf + 170 * ncomps + c);
        const auto *kf_171 = buffer.data(kf + 171 * ncomps + c);
        const auto *kf_172 = buffer.data(kf + 172 * ncomps + c);
        const auto *kf_173 = buffer.data(kf + 173 * ncomps + c);
        const auto *kf_174 = buffer.data(kf + 174 * ncomps + c);
        const auto *kf_175 = buffer.data(kf + 175 * ncomps + c);
        const auto *kf_176 = buffer.data(kf + 176 * ncomps + c);
        const auto *kf_177 = buffer.data(kf + 177 * ncomps + c);
        const auto *kf_178 = buffer.data(kf + 178 * ncomps + c);
        const auto *kf_179 = buffer.data(kf + 179 * ncomps + c);
        const auto *kf_180 = buffer.data(kf + 180 * ncomps + c);
        const auto *kf_181 = buffer.data(kf + 181 * ncomps + c);
        const auto *kf_182 = buffer.data(kf + 182 * ncomps + c);
        const auto *kf_183 = buffer.data(kf + 183 * ncomps + c);
        const auto *kf_184 = buffer.data(kf + 184 * ncomps + c);
        const auto *kf_185 = buffer.data(kf + 185 * ncomps + c);
        const auto *kf_186 = buffer.data(kf + 186 * ncomps + c);
        const auto *kf_187 = buffer.data(kf + 187 * ncomps + c);
        const auto *kf_188 = buffer.data(kf + 188 * ncomps + c);
        const auto *kf_189 = buffer.data(kf + 189 * ncomps + c);
        const auto *kf_190 = buffer.data(kf + 190 * ncomps + c);
        const auto *kf_191 = buffer.data(kf + 191 * ncomps + c);
        const auto *kf_192 = buffer.data(kf + 192 * ncomps + c);
        const auto *kf_193 = buffer.data(kf + 193 * ncomps + c);
        const auto *kf_194 = buffer.data(kf + 194 * ncomps + c);

        const auto *lf_100 = buffer.data(lf + 100 * ncomps + c);
        const auto *lf_101 = buffer.data(lf + 101 * ncomps + c);
        const auto *lf_102 = buffer.data(lf + 102 * ncomps + c);
        const auto *lf_103 = buffer.data(lf + 103 * ncomps + c);
        const auto *lf_104 = buffer.data(lf + 104 * ncomps + c);
        const auto *lf_105 = buffer.data(lf + 105 * ncomps + c);
        const auto *lf_106 = buffer.data(lf + 106 * ncomps + c);
        const auto *lf_107 = buffer.data(lf + 107 * ncomps + c);
        const auto *lf_108 = buffer.data(lf + 108 * ncomps + c);
        const auto *lf_109 = buffer.data(lf + 109 * ncomps + c);
        const auto *lf_110 = buffer.data(lf + 110 * ncomps + c);
        const auto *lf_111 = buffer.data(lf + 111 * ncomps + c);
        const auto *lf_112 = buffer.data(lf + 112 * ncomps + c);
        const auto *lf_113 = buffer.data(lf + 113 * ncomps + c);
        const auto *lf_114 = buffer.data(lf + 114 * ncomps + c);
        const auto *lf_115 = buffer.data(lf + 115 * ncomps + c);
        const auto *lf_116 = buffer.data(lf + 116 * ncomps + c);
        const auto *lf_117 = buffer.data(lf + 117 * ncomps + c);
        const auto *lf_118 = buffer.data(lf + 118 * ncomps + c);
        const auto *lf_119 = buffer.data(lf + 119 * ncomps + c);
        const auto *lf_120 = buffer.data(lf + 120 * ncomps + c);
        const auto *lf_121 = buffer.data(lf + 121 * ncomps + c);
        const auto *lf_122 = buffer.data(lf + 122 * ncomps + c);
        const auto *lf_123 = buffer.data(lf + 123 * ncomps + c);
        const auto *lf_124 = buffer.data(lf + 124 * ncomps + c);
        const auto *lf_125 = buffer.data(lf + 125 * ncomps + c);
        const auto *lf_126 = buffer.data(lf + 126 * ncomps + c);
        const auto *lf_127 = buffer.data(lf + 127 * ncomps + c);
        const auto *lf_128 = buffer.data(lf + 128 * ncomps + c);
        const auto *lf_129 = buffer.data(lf + 129 * ncomps + c);
        const auto *lf_130 = buffer.data(lf + 130 * ncomps + c);
        const auto *lf_131 = buffer.data(lf + 131 * ncomps + c);
        const auto *lf_132 = buffer.data(lf + 132 * ncomps + c);
        const auto *lf_133 = buffer.data(lf + 133 * ncomps + c);
        const auto *lf_134 = buffer.data(lf + 134 * ncomps + c);
        const auto *lf_135 = buffer.data(lf + 135 * ncomps + c);
        const auto *lf_136 = buffer.data(lf + 136 * ncomps + c);
        const auto *lf_137 = buffer.data(lf + 137 * ncomps + c);
        const auto *lf_138 = buffer.data(lf + 138 * ncomps + c);
        const auto *lf_139 = buffer.data(lf + 139 * ncomps + c);
        const auto *lf_140 = buffer.data(lf + 140 * ncomps + c);
        const auto *lf_141 = buffer.data(lf + 141 * ncomps + c);
        const auto *lf_142 = buffer.data(lf + 142 * ncomps + c);
        const auto *lf_143 = buffer.data(lf + 143 * ncomps + c);
        const auto *lf_144 = buffer.data(lf + 144 * ncomps + c);
        const auto *lf_145 = buffer.data(lf + 145 * ncomps + c);
        const auto *lf_146 = buffer.data(lf + 146 * ncomps + c);
        const auto *lf_147 = buffer.data(lf + 147 * ncomps + c);
        const auto *lf_148 = buffer.data(lf + 148 * ncomps + c);
        const auto *lf_149 = buffer.data(lf + 149 * ncomps + c);
        const auto *lf_150 = buffer.data(lf + 150 * ncomps + c);
        const auto *lf_151 = buffer.data(lf + 151 * ncomps + c);
        const auto *lf_152 = buffer.data(lf + 152 * ncomps + c);
        const auto *lf_153 = buffer.data(lf + 153 * ncomps + c);
        const auto *lf_154 = buffer.data(lf + 154 * ncomps + c);
        const auto *lf_155 = buffer.data(lf + 155 * ncomps + c);
        const auto *lf_156 = buffer.data(lf + 156 * ncomps + c);
        const auto *lf_157 = buffer.data(lf + 157 * ncomps + c);
        const auto *lf_158 = buffer.data(lf + 158 * ncomps + c);
        const auto *lf_159 = buffer.data(lf + 159 * ncomps + c);
        const auto *lf_160 = buffer.data(lf + 160 * ncomps + c);
        const auto *lf_161 = buffer.data(lf + 161 * ncomps + c);
        const auto *lf_162 = buffer.data(lf + 162 * ncomps + c);
        const auto *lf_163 = buffer.data(lf + 163 * ncomps + c);
        const auto *lf_164 = buffer.data(lf + 164 * ncomps + c);
        const auto *lf_165 = buffer.data(lf + 165 * ncomps + c);
        const auto *lf_166 = buffer.data(lf + 166 * ncomps + c);
        const auto *lf_167 = buffer.data(lf + 167 * ncomps + c);
        const auto *lf_168 = buffer.data(lf + 168 * ncomps + c);
        const auto *lf_169 = buffer.data(lf + 169 * ncomps + c);
        const auto *lf_170 = buffer.data(lf + 170 * ncomps + c);
        const auto *lf_171 = buffer.data(lf + 171 * ncomps + c);
        const auto *lf_172 = buffer.data(lf + 172 * ncomps + c);
        const auto *lf_173 = buffer.data(lf + 173 * ncomps + c);
        const auto *lf_174 = buffer.data(lf + 174 * ncomps + c);
        const auto *lf_175 = buffer.data(lf + 175 * ncomps + c);
        const auto *lf_176 = buffer.data(lf + 176 * ncomps + c);
        const auto *lf_177 = buffer.data(lf + 177 * ncomps + c);
        const auto *lf_178 = buffer.data(lf + 178 * ncomps + c);
        const auto *lf_179 = buffer.data(lf + 179 * ncomps + c);
        const auto *lf_180 = buffer.data(lf + 180 * ncomps + c);
        const auto *lf_181 = buffer.data(lf + 181 * ncomps + c);
        const auto *lf_182 = buffer.data(lf + 182 * ncomps + c);
        const auto *lf_183 = buffer.data(lf + 183 * ncomps + c);
        const auto *lf_184 = buffer.data(lf + 184 * ncomps + c);
        const auto *lf_185 = buffer.data(lf + 185 * ncomps + c);
        const auto *lf_186 = buffer.data(lf + 186 * ncomps + c);
        const auto *lf_187 = buffer.data(lf + 187 * ncomps + c);
        const auto *lf_188 = buffer.data(lf + 188 * ncomps + c);
        const auto *lf_189 = buffer.data(lf + 189 * ncomps + c);
        const auto *lf_190 = buffer.data(lf + 190 * ncomps + c);
        const auto *lf_191 = buffer.data(lf + 191 * ncomps + c);
        const auto *lf_192 = buffer.data(lf + 192 * ncomps + c);
        const auto *lf_193 = buffer.data(lf + 193 * ncomps + c);
        const auto *lf_194 = buffer.data(lf + 194 * ncomps + c);
        const auto *lf_196 = buffer.data(lf + 196 * ncomps + c);
        const auto *lf_197 = buffer.data(lf + 197 * ncomps + c);
        const auto *lf_198 = buffer.data(lf + 198 * ncomps + c);
        const auto *lf_199 = buffer.data(lf + 199 * ncomps + c);
        const auto *lf_209 = buffer.data(lf + 209 * ncomps + c);
        const auto *lf_216 = buffer.data(lf + 216 * ncomps + c);
        const auto *lf_217 = buffer.data(lf + 217 * ncomps + c);
        const auto *lf_218 = buffer.data(lf + 218 * ncomps + c);
        const auto *lf_219 = buffer.data(lf + 219 * ncomps + c);
        const auto *lf_226 = buffer.data(lf + 226 * ncomps + c);
        const auto *lf_227 = buffer.data(lf + 227 * ncomps + c);
        const auto *lf_228 = buffer.data(lf + 228 * ncomps + c);
        const auto *lf_229 = buffer.data(lf + 229 * ncomps + c);
        const auto *lf_236 = buffer.data(lf + 236 * ncomps + c);
        const auto *lf_237 = buffer.data(lf + 237 * ncomps + c);
        const auto *lf_238 = buffer.data(lf + 238 * ncomps + c);
        const auto *lf_239 = buffer.data(lf + 239 * ncomps + c);
        const auto *lf_246 = buffer.data(lf + 246 * ncomps + c);
        const auto *lf_247 = buffer.data(lf + 247 * ncomps + c);
        const auto *lf_248 = buffer.data(lf + 248 * ncomps + c);
        const auto *lf_249 = buffer.data(lf + 249 * ncomps + c);
        const auto *lf_259 = buffer.data(lf + 259 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, kf_96, kf_97, kf_98, \
                         kf_99, lf_136, lf_137, lf_138, lf_139, \
                         lf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_y[k] * kf_96[k]
                       + lf_136[k];

            t_146[k] = ab_y[k] * kf_97[k]
                       + lf_137[k];

            t_147[k] = ab_y[k] * kf_98[k]
                       + lf_138[k];

            t_148[k] = ab_y[k] * kf_99[k]
                       + lf_139[k];

            t_149[k] = ab_z[k] * kf_99[k]
                       + lf_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, kf_100, kf_101, kf_102, \
                         kf_103, kf_104, lf_100, lf_101, lf_102, lf_103, \
                         lf_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * kf_100[k]
                       + lf_100[k];

            t_151[k] = ab_x[k] * kf_101[k]
                       + lf_101[k];

            t_152[k] = ab_x[k] * kf_102[k]
                       + lf_102[k];

            t_153[k] = ab_x[k] * kf_103[k]
                       + lf_103[k];

            t_154[k] = ab_x[k] * kf_104[k]
                       + lf_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, kf_105, kf_106, kf_107, \
                         kf_108, kf_109, lf_105, lf_106, lf_107, lf_108, \
                         lf_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * kf_105[k]
                       + lf_105[k];

            t_156[k] = ab_x[k] * kf_106[k]
                       + lf_106[k];

            t_157[k] = ab_x[k] * kf_107[k]
                       + lf_107[k];

            t_158[k] = ab_x[k] * kf_108[k]
                       + lf_108[k];

            t_159[k] = ab_x[k] * kf_109[k]
                       + lf_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, kf_106, kf_107, \
                         kf_108, kf_109, lf_156, lf_157, lf_158, lf_159, \
                         lf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * kf_106[k]
                       + lf_156[k];

            t_161[k] = ab_y[k] * kf_107[k]
                       + lf_157[k];

            t_162[k] = ab_y[k] * kf_108[k]
                       + lf_158[k];

            t_163[k] = ab_y[k] * kf_109[k]
                       + lf_159[k];

            t_164[k] = ab_z[k] * kf_109[k]
                       + lf_169[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, kf_110, kf_111, kf_112, \
                         kf_113, kf_114, lf_110, lf_111, lf_112, lf_113, \
                         lf_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * kf_110[k]
                       + lf_110[k];

            t_166[k] = ab_x[k] * kf_111[k]
                       + lf_111[k];

            t_167[k] = ab_x[k] * kf_112[k]
                       + lf_112[k];

            t_168[k] = ab_x[k] * kf_113[k]
                       + lf_113[k];

            t_169[k] = ab_x[k] * kf_114[k]
                       + lf_114[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, kf_115, kf_116, kf_117, \
                         kf_118, kf_119, lf_115, lf_116, lf_117, lf_118, \
                         lf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * kf_115[k]
                       + lf_115[k];

            t_171[k] = ab_x[k] * kf_116[k]
                       + lf_116[k];

            t_172[k] = ab_x[k] * kf_117[k]
                       + lf_117[k];

            t_173[k] = ab_x[k] * kf_118[k]
                       + lf_118[k];

            t_174[k] = ab_x[k] * kf_119[k]
                       + lf_119[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, kf_116, kf_117, \
                         kf_118, kf_119, lf_166, lf_167, lf_168, lf_169, \
                         lf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_y[k] * kf_116[k]
                       + lf_166[k];

            t_176[k] = ab_y[k] * kf_117[k]
                       + lf_167[k];

            t_177[k] = ab_y[k] * kf_118[k]
                       + lf_168[k];

            t_178[k] = ab_y[k] * kf_119[k]
                       + lf_169[k];

            t_179[k] = ab_z[k] * kf_119[k]
                       + lf_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, kf_120, kf_121, kf_122, \
                         kf_123, kf_124, lf_120, lf_121, lf_122, lf_123, \
                         lf_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * kf_120[k]
                       + lf_120[k];

            t_181[k] = ab_x[k] * kf_121[k]
                       + lf_121[k];

            t_182[k] = ab_x[k] * kf_122[k]
                       + lf_122[k];

            t_183[k] = ab_x[k] * kf_123[k]
                       + lf_123[k];

            t_184[k] = ab_x[k] * kf_124[k]
                       + lf_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, kf_125, kf_126, kf_127, \
                         kf_128, kf_129, lf_125, lf_126, lf_127, lf_128, \
                         lf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * kf_125[k]
                       + lf_125[k];

            t_186[k] = ab_x[k] * kf_126[k]
                       + lf_126[k];

            t_187[k] = ab_x[k] * kf_127[k]
                       + lf_127[k];

            t_188[k] = ab_x[k] * kf_128[k]
                       + lf_128[k];

            t_189[k] = ab_x[k] * kf_129[k]
                       + lf_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, kf_126, kf_127, \
                         kf_128, kf_129, lf_176, lf_177, lf_178, lf_179, \
                         lf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * kf_126[k]
                       + lf_176[k];

            t_191[k] = ab_y[k] * kf_127[k]
                       + lf_177[k];

            t_192[k] = ab_y[k] * kf_128[k]
                       + lf_178[k];

            t_193[k] = ab_y[k] * kf_129[k]
                       + lf_179[k];

            t_194[k] = ab_z[k] * kf_129[k]
                       + lf_189[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, kf_130, kf_131, kf_132, \
                         kf_133, kf_134, lf_130, lf_131, lf_132, lf_133, \
                         lf_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * kf_130[k]
                       + lf_130[k];

            t_196[k] = ab_x[k] * kf_131[k]
                       + lf_131[k];

            t_197[k] = ab_x[k] * kf_132[k]
                       + lf_132[k];

            t_198[k] = ab_x[k] * kf_133[k]
                       + lf_133[k];

            t_199[k] = ab_x[k] * kf_134[k]
                       + lf_134[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, kf_135, kf_136, kf_137, \
                         kf_138, kf_139, lf_135, lf_136, lf_137, lf_138, \
                         lf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * kf_135[k]
                       + lf_135[k];

            t_201[k] = ab_x[k] * kf_136[k]
                       + lf_136[k];

            t_202[k] = ab_x[k] * kf_137[k]
                       + lf_137[k];

            t_203[k] = ab_x[k] * kf_138[k]
                       + lf_138[k];

            t_204[k] = ab_x[k] * kf_139[k]
                       + lf_139[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, kf_136, kf_137, \
                         kf_138, kf_139, lf_186, lf_187, lf_188, lf_189, \
                         lf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_y[k] * kf_136[k]
                       + lf_186[k];

            t_206[k] = ab_y[k] * kf_137[k]
                       + lf_187[k];

            t_207[k] = ab_y[k] * kf_138[k]
                       + lf_188[k];

            t_208[k] = ab_y[k] * kf_139[k]
                       + lf_189[k];

            t_209[k] = ab_z[k] * kf_139[k]
                       + lf_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, kf_140, kf_141, kf_142, \
                         kf_143, kf_144, lf_140, lf_141, lf_142, lf_143, \
                         lf_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * kf_140[k]
                       + lf_140[k];

            t_211[k] = ab_x[k] * kf_141[k]
                       + lf_141[k];

            t_212[k] = ab_x[k] * kf_142[k]
                       + lf_142[k];

            t_213[k] = ab_x[k] * kf_143[k]
                       + lf_143[k];

            t_214[k] = ab_x[k] * kf_144[k]
                       + lf_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, kf_145, kf_146, kf_147, \
                         kf_148, kf_149, lf_145, lf_146, lf_147, lf_148, \
                         lf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * kf_145[k]
                       + lf_145[k];

            t_216[k] = ab_x[k] * kf_146[k]
                       + lf_146[k];

            t_217[k] = ab_x[k] * kf_147[k]
                       + lf_147[k];

            t_218[k] = ab_x[k] * kf_148[k]
                       + lf_148[k];

            t_219[k] = ab_x[k] * kf_149[k]
                       + lf_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, kf_146, kf_147, \
                         kf_148, kf_149, lf_196, lf_197, lf_198, lf_199, \
                         lf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * kf_146[k]
                       + lf_196[k];

            t_221[k] = ab_y[k] * kf_147[k]
                       + lf_197[k];

            t_222[k] = ab_y[k] * kf_148[k]
                       + lf_198[k];

            t_223[k] = ab_y[k] * kf_149[k]
                       + lf_199[k];

            t_224[k] = ab_z[k] * kf_149[k]
                       + lf_209[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, kf_150, kf_151, kf_152, \
                         kf_153, kf_154, lf_150, lf_151, lf_152, lf_153, \
                         lf_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * kf_150[k]
                       + lf_150[k];

            t_226[k] = ab_x[k] * kf_151[k]
                       + lf_151[k];

            t_227[k] = ab_x[k] * kf_152[k]
                       + lf_152[k];

            t_228[k] = ab_x[k] * kf_153[k]
                       + lf_153[k];

            t_229[k] = ab_x[k] * kf_154[k]
                       + lf_154[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, kf_155, kf_156, kf_157, \
                         kf_158, kf_159, lf_155, lf_156, lf_157, lf_158, \
                         lf_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * kf_155[k]
                       + lf_155[k];

            t_231[k] = ab_x[k] * kf_156[k]
                       + lf_156[k];

            t_232[k] = ab_x[k] * kf_157[k]
                       + lf_157[k];

            t_233[k] = ab_x[k] * kf_158[k]
                       + lf_158[k];

            t_234[k] = ab_x[k] * kf_159[k]
                       + lf_159[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, ab_z, kf_156, kf_157, \
                         kf_158, kf_159, lf_216, lf_217, lf_218, lf_219, \
                         lf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_y[k] * kf_156[k]
                       + lf_216[k];

            t_236[k] = ab_y[k] * kf_157[k]
                       + lf_217[k];

            t_237[k] = ab_y[k] * kf_158[k]
                       + lf_218[k];

            t_238[k] = ab_y[k] * kf_159[k]
                       + lf_219[k];

            t_239[k] = ab_z[k] * kf_159[k]
                       + lf_229[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, kf_160, kf_161, kf_162, \
                         kf_163, kf_164, lf_160, lf_161, lf_162, lf_163, \
                         lf_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * kf_160[k]
                       + lf_160[k];

            t_241[k] = ab_x[k] * kf_161[k]
                       + lf_161[k];

            t_242[k] = ab_x[k] * kf_162[k]
                       + lf_162[k];

            t_243[k] = ab_x[k] * kf_163[k]
                       + lf_163[k];

            t_244[k] = ab_x[k] * kf_164[k]
                       + lf_164[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, kf_165, kf_166, kf_167, \
                         kf_168, kf_169, lf_165, lf_166, lf_167, lf_168, \
                         lf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * kf_165[k]
                       + lf_165[k];

            t_246[k] = ab_x[k] * kf_166[k]
                       + lf_166[k];

            t_247[k] = ab_x[k] * kf_167[k]
                       + lf_167[k];

            t_248[k] = ab_x[k] * kf_168[k]
                       + lf_168[k];

            t_249[k] = ab_x[k] * kf_169[k]
                       + lf_169[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, ab_z, kf_166, kf_167, \
                         kf_168, kf_169, lf_226, lf_227, lf_228, lf_229, \
                         lf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_y[k] * kf_166[k]
                       + lf_226[k];

            t_251[k] = ab_y[k] * kf_167[k]
                       + lf_227[k];

            t_252[k] = ab_y[k] * kf_168[k]
                       + lf_228[k];

            t_253[k] = ab_y[k] * kf_169[k]
                       + lf_229[k];

            t_254[k] = ab_z[k] * kf_169[k]
                       + lf_239[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, kf_170, kf_171, kf_172, \
                         kf_173, kf_174, lf_170, lf_171, lf_172, lf_173, \
                         lf_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * kf_170[k]
                       + lf_170[k];

            t_256[k] = ab_x[k] * kf_171[k]
                       + lf_171[k];

            t_257[k] = ab_x[k] * kf_172[k]
                       + lf_172[k];

            t_258[k] = ab_x[k] * kf_173[k]
                       + lf_173[k];

            t_259[k] = ab_x[k] * kf_174[k]
                       + lf_174[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, kf_175, kf_176, kf_177, \
                         kf_178, kf_179, lf_175, lf_176, lf_177, lf_178, \
                         lf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * kf_175[k]
                       + lf_175[k];

            t_261[k] = ab_x[k] * kf_176[k]
                       + lf_176[k];

            t_262[k] = ab_x[k] * kf_177[k]
                       + lf_177[k];

            t_263[k] = ab_x[k] * kf_178[k]
                       + lf_178[k];

            t_264[k] = ab_x[k] * kf_179[k]
                       + lf_179[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, ab_z, kf_176, kf_177, \
                         kf_178, kf_179, lf_236, lf_237, lf_238, lf_239, \
                         lf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_y[k] * kf_176[k]
                       + lf_236[k];

            t_266[k] = ab_y[k] * kf_177[k]
                       + lf_237[k];

            t_267[k] = ab_y[k] * kf_178[k]
                       + lf_238[k];

            t_268[k] = ab_y[k] * kf_179[k]
                       + lf_239[k];

            t_269[k] = ab_z[k] * kf_179[k]
                       + lf_249[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, kf_180, kf_181, kf_182, \
                         kf_183, kf_184, lf_180, lf_181, lf_182, lf_183, \
                         lf_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * kf_180[k]
                       + lf_180[k];

            t_271[k] = ab_x[k] * kf_181[k]
                       + lf_181[k];

            t_272[k] = ab_x[k] * kf_182[k]
                       + lf_182[k];

            t_273[k] = ab_x[k] * kf_183[k]
                       + lf_183[k];

            t_274[k] = ab_x[k] * kf_184[k]
                       + lf_184[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, kf_185, kf_186, kf_187, \
                         kf_188, kf_189, lf_185, lf_186, lf_187, lf_188, \
                         lf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * kf_185[k]
                       + lf_185[k];

            t_276[k] = ab_x[k] * kf_186[k]
                       + lf_186[k];

            t_277[k] = ab_x[k] * kf_187[k]
                       + lf_187[k];

            t_278[k] = ab_x[k] * kf_188[k]
                       + lf_188[k];

            t_279[k] = ab_x[k] * kf_189[k]
                       + lf_189[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, ab_z, kf_186, kf_187, \
                         kf_188, kf_189, lf_246, lf_247, lf_248, lf_249, \
                         lf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_y[k] * kf_186[k]
                       + lf_246[k];

            t_281[k] = ab_y[k] * kf_187[k]
                       + lf_247[k];

            t_282[k] = ab_y[k] * kf_188[k]
                       + lf_248[k];

            t_283[k] = ab_y[k] * kf_189[k]
                       + lf_249[k];

            t_284[k] = ab_z[k] * kf_189[k]
                       + lf_259[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, kf_190, kf_191, kf_192, \
                         kf_193, kf_194, lf_190, lf_191, lf_192, lf_193, \
                         lf_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * kf_190[k]
                       + lf_190[k];

            t_286[k] = ab_x[k] * kf_191[k]
                       + lf_191[k];

            t_287[k] = ab_x[k] * kf_192[k]
                       + lf_192[k];

            t_288[k] = ab_x[k] * kf_193[k]
                       + lf_193[k];

            t_289[k] = ab_x[k] * kf_194[k]
                       + lf_194[k];
        }
    }
}

static auto
compute_hrr_kg_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kf, const size_t lf,
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

        const auto *kf_195 = buffer.data(kf + 195 * ncomps + c);
        const auto *kf_196 = buffer.data(kf + 196 * ncomps + c);
        const auto *kf_197 = buffer.data(kf + 197 * ncomps + c);
        const auto *kf_198 = buffer.data(kf + 198 * ncomps + c);
        const auto *kf_199 = buffer.data(kf + 199 * ncomps + c);
        const auto *kf_200 = buffer.data(kf + 200 * ncomps + c);
        const auto *kf_201 = buffer.data(kf + 201 * ncomps + c);
        const auto *kf_202 = buffer.data(kf + 202 * ncomps + c);
        const auto *kf_203 = buffer.data(kf + 203 * ncomps + c);
        const auto *kf_204 = buffer.data(kf + 204 * ncomps + c);
        const auto *kf_205 = buffer.data(kf + 205 * ncomps + c);
        const auto *kf_206 = buffer.data(kf + 206 * ncomps + c);
        const auto *kf_207 = buffer.data(kf + 207 * ncomps + c);
        const auto *kf_208 = buffer.data(kf + 208 * ncomps + c);
        const auto *kf_209 = buffer.data(kf + 209 * ncomps + c);
        const auto *kf_210 = buffer.data(kf + 210 * ncomps + c);
        const auto *kf_211 = buffer.data(kf + 211 * ncomps + c);
        const auto *kf_212 = buffer.data(kf + 212 * ncomps + c);
        const auto *kf_213 = buffer.data(kf + 213 * ncomps + c);
        const auto *kf_214 = buffer.data(kf + 214 * ncomps + c);
        const auto *kf_215 = buffer.data(kf + 215 * ncomps + c);
        const auto *kf_216 = buffer.data(kf + 216 * ncomps + c);
        const auto *kf_217 = buffer.data(kf + 217 * ncomps + c);
        const auto *kf_218 = buffer.data(kf + 218 * ncomps + c);
        const auto *kf_219 = buffer.data(kf + 219 * ncomps + c);
        const auto *kf_220 = buffer.data(kf + 220 * ncomps + c);
        const auto *kf_221 = buffer.data(kf + 221 * ncomps + c);
        const auto *kf_222 = buffer.data(kf + 222 * ncomps + c);
        const auto *kf_223 = buffer.data(kf + 223 * ncomps + c);
        const auto *kf_224 = buffer.data(kf + 224 * ncomps + c);
        const auto *kf_225 = buffer.data(kf + 225 * ncomps + c);
        const auto *kf_226 = buffer.data(kf + 226 * ncomps + c);
        const auto *kf_227 = buffer.data(kf + 227 * ncomps + c);
        const auto *kf_228 = buffer.data(kf + 228 * ncomps + c);
        const auto *kf_229 = buffer.data(kf + 229 * ncomps + c);
        const auto *kf_230 = buffer.data(kf + 230 * ncomps + c);
        const auto *kf_231 = buffer.data(kf + 231 * ncomps + c);
        const auto *kf_232 = buffer.data(kf + 232 * ncomps + c);
        const auto *kf_233 = buffer.data(kf + 233 * ncomps + c);
        const auto *kf_234 = buffer.data(kf + 234 * ncomps + c);
        const auto *kf_235 = buffer.data(kf + 235 * ncomps + c);
        const auto *kf_236 = buffer.data(kf + 236 * ncomps + c);
        const auto *kf_237 = buffer.data(kf + 237 * ncomps + c);
        const auto *kf_238 = buffer.data(kf + 238 * ncomps + c);
        const auto *kf_239 = buffer.data(kf + 239 * ncomps + c);
        const auto *kf_240 = buffer.data(kf + 240 * ncomps + c);
        const auto *kf_241 = buffer.data(kf + 241 * ncomps + c);
        const auto *kf_242 = buffer.data(kf + 242 * ncomps + c);
        const auto *kf_243 = buffer.data(kf + 243 * ncomps + c);
        const auto *kf_244 = buffer.data(kf + 244 * ncomps + c);
        const auto *kf_245 = buffer.data(kf + 245 * ncomps + c);
        const auto *kf_246 = buffer.data(kf + 246 * ncomps + c);
        const auto *kf_247 = buffer.data(kf + 247 * ncomps + c);
        const auto *kf_248 = buffer.data(kf + 248 * ncomps + c);
        const auto *kf_249 = buffer.data(kf + 249 * ncomps + c);
        const auto *kf_250 = buffer.data(kf + 250 * ncomps + c);
        const auto *kf_251 = buffer.data(kf + 251 * ncomps + c);
        const auto *kf_252 = buffer.data(kf + 252 * ncomps + c);
        const auto *kf_253 = buffer.data(kf + 253 * ncomps + c);
        const auto *kf_254 = buffer.data(kf + 254 * ncomps + c);
        const auto *kf_255 = buffer.data(kf + 255 * ncomps + c);
        const auto *kf_256 = buffer.data(kf + 256 * ncomps + c);
        const auto *kf_257 = buffer.data(kf + 257 * ncomps + c);
        const auto *kf_258 = buffer.data(kf + 258 * ncomps + c);
        const auto *kf_259 = buffer.data(kf + 259 * ncomps + c);
        const auto *kf_260 = buffer.data(kf + 260 * ncomps + c);
        const auto *kf_261 = buffer.data(kf + 261 * ncomps + c);
        const auto *kf_262 = buffer.data(kf + 262 * ncomps + c);
        const auto *kf_263 = buffer.data(kf + 263 * ncomps + c);
        const auto *kf_264 = buffer.data(kf + 264 * ncomps + c);
        const auto *kf_265 = buffer.data(kf + 265 * ncomps + c);
        const auto *kf_266 = buffer.data(kf + 266 * ncomps + c);
        const auto *kf_267 = buffer.data(kf + 267 * ncomps + c);
        const auto *kf_268 = buffer.data(kf + 268 * ncomps + c);
        const auto *kf_269 = buffer.data(kf + 269 * ncomps + c);
        const auto *kf_270 = buffer.data(kf + 270 * ncomps + c);
        const auto *kf_271 = buffer.data(kf + 271 * ncomps + c);
        const auto *kf_272 = buffer.data(kf + 272 * ncomps + c);
        const auto *kf_273 = buffer.data(kf + 273 * ncomps + c);
        const auto *kf_274 = buffer.data(kf + 274 * ncomps + c);
        const auto *kf_275 = buffer.data(kf + 275 * ncomps + c);
        const auto *kf_276 = buffer.data(kf + 276 * ncomps + c);
        const auto *kf_277 = buffer.data(kf + 277 * ncomps + c);
        const auto *kf_278 = buffer.data(kf + 278 * ncomps + c);
        const auto *kf_279 = buffer.data(kf + 279 * ncomps + c);
        const auto *kf_280 = buffer.data(kf + 280 * ncomps + c);
        const auto *kf_281 = buffer.data(kf + 281 * ncomps + c);
        const auto *kf_282 = buffer.data(kf + 282 * ncomps + c);
        const auto *kf_283 = buffer.data(kf + 283 * ncomps + c);
        const auto *kf_284 = buffer.data(kf + 284 * ncomps + c);
        const auto *kf_285 = buffer.data(kf + 285 * ncomps + c);
        const auto *kf_286 = buffer.data(kf + 286 * ncomps + c);
        const auto *kf_287 = buffer.data(kf + 287 * ncomps + c);
        const auto *kf_288 = buffer.data(kf + 288 * ncomps + c);
        const auto *kf_289 = buffer.data(kf + 289 * ncomps + c);

        const auto *lf_195 = buffer.data(lf + 195 * ncomps + c);
        const auto *lf_196 = buffer.data(lf + 196 * ncomps + c);
        const auto *lf_197 = buffer.data(lf + 197 * ncomps + c);
        const auto *lf_198 = buffer.data(lf + 198 * ncomps + c);
        const auto *lf_199 = buffer.data(lf + 199 * ncomps + c);
        const auto *lf_200 = buffer.data(lf + 200 * ncomps + c);
        const auto *lf_201 = buffer.data(lf + 201 * ncomps + c);
        const auto *lf_202 = buffer.data(lf + 202 * ncomps + c);
        const auto *lf_203 = buffer.data(lf + 203 * ncomps + c);
        const auto *lf_204 = buffer.data(lf + 204 * ncomps + c);
        const auto *lf_205 = buffer.data(lf + 205 * ncomps + c);
        const auto *lf_206 = buffer.data(lf + 206 * ncomps + c);
        const auto *lf_207 = buffer.data(lf + 207 * ncomps + c);
        const auto *lf_208 = buffer.data(lf + 208 * ncomps + c);
        const auto *lf_209 = buffer.data(lf + 209 * ncomps + c);
        const auto *lf_210 = buffer.data(lf + 210 * ncomps + c);
        const auto *lf_211 = buffer.data(lf + 211 * ncomps + c);
        const auto *lf_212 = buffer.data(lf + 212 * ncomps + c);
        const auto *lf_213 = buffer.data(lf + 213 * ncomps + c);
        const auto *lf_214 = buffer.data(lf + 214 * ncomps + c);
        const auto *lf_215 = buffer.data(lf + 215 * ncomps + c);
        const auto *lf_216 = buffer.data(lf + 216 * ncomps + c);
        const auto *lf_217 = buffer.data(lf + 217 * ncomps + c);
        const auto *lf_218 = buffer.data(lf + 218 * ncomps + c);
        const auto *lf_219 = buffer.data(lf + 219 * ncomps + c);
        const auto *lf_220 = buffer.data(lf + 220 * ncomps + c);
        const auto *lf_221 = buffer.data(lf + 221 * ncomps + c);
        const auto *lf_222 = buffer.data(lf + 222 * ncomps + c);
        const auto *lf_223 = buffer.data(lf + 223 * ncomps + c);
        const auto *lf_224 = buffer.data(lf + 224 * ncomps + c);
        const auto *lf_225 = buffer.data(lf + 225 * ncomps + c);
        const auto *lf_226 = buffer.data(lf + 226 * ncomps + c);
        const auto *lf_227 = buffer.data(lf + 227 * ncomps + c);
        const auto *lf_228 = buffer.data(lf + 228 * ncomps + c);
        const auto *lf_229 = buffer.data(lf + 229 * ncomps + c);
        const auto *lf_230 = buffer.data(lf + 230 * ncomps + c);
        const auto *lf_231 = buffer.data(lf + 231 * ncomps + c);
        const auto *lf_232 = buffer.data(lf + 232 * ncomps + c);
        const auto *lf_233 = buffer.data(lf + 233 * ncomps + c);
        const auto *lf_234 = buffer.data(lf + 234 * ncomps + c);
        const auto *lf_235 = buffer.data(lf + 235 * ncomps + c);
        const auto *lf_236 = buffer.data(lf + 236 * ncomps + c);
        const auto *lf_237 = buffer.data(lf + 237 * ncomps + c);
        const auto *lf_238 = buffer.data(lf + 238 * ncomps + c);
        const auto *lf_239 = buffer.data(lf + 239 * ncomps + c);
        const auto *lf_240 = buffer.data(lf + 240 * ncomps + c);
        const auto *lf_241 = buffer.data(lf + 241 * ncomps + c);
        const auto *lf_242 = buffer.data(lf + 242 * ncomps + c);
        const auto *lf_243 = buffer.data(lf + 243 * ncomps + c);
        const auto *lf_244 = buffer.data(lf + 244 * ncomps + c);
        const auto *lf_245 = buffer.data(lf + 245 * ncomps + c);
        const auto *lf_246 = buffer.data(lf + 246 * ncomps + c);
        const auto *lf_247 = buffer.data(lf + 247 * ncomps + c);
        const auto *lf_248 = buffer.data(lf + 248 * ncomps + c);
        const auto *lf_249 = buffer.data(lf + 249 * ncomps + c);
        const auto *lf_250 = buffer.data(lf + 250 * ncomps + c);
        const auto *lf_251 = buffer.data(lf + 251 * ncomps + c);
        const auto *lf_252 = buffer.data(lf + 252 * ncomps + c);
        const auto *lf_253 = buffer.data(lf + 253 * ncomps + c);
        const auto *lf_254 = buffer.data(lf + 254 * ncomps + c);
        const auto *lf_255 = buffer.data(lf + 255 * ncomps + c);
        const auto *lf_256 = buffer.data(lf + 256 * ncomps + c);
        const auto *lf_257 = buffer.data(lf + 257 * ncomps + c);
        const auto *lf_258 = buffer.data(lf + 258 * ncomps + c);
        const auto *lf_259 = buffer.data(lf + 259 * ncomps + c);
        const auto *lf_260 = buffer.data(lf + 260 * ncomps + c);
        const auto *lf_261 = buffer.data(lf + 261 * ncomps + c);
        const auto *lf_262 = buffer.data(lf + 262 * ncomps + c);
        const auto *lf_263 = buffer.data(lf + 263 * ncomps + c);
        const auto *lf_264 = buffer.data(lf + 264 * ncomps + c);
        const auto *lf_265 = buffer.data(lf + 265 * ncomps + c);
        const auto *lf_266 = buffer.data(lf + 266 * ncomps + c);
        const auto *lf_267 = buffer.data(lf + 267 * ncomps + c);
        const auto *lf_268 = buffer.data(lf + 268 * ncomps + c);
        const auto *lf_269 = buffer.data(lf + 269 * ncomps + c);
        const auto *lf_270 = buffer.data(lf + 270 * ncomps + c);
        const auto *lf_271 = buffer.data(lf + 271 * ncomps + c);
        const auto *lf_272 = buffer.data(lf + 272 * ncomps + c);
        const auto *lf_273 = buffer.data(lf + 273 * ncomps + c);
        const auto *lf_274 = buffer.data(lf + 274 * ncomps + c);
        const auto *lf_275 = buffer.data(lf + 275 * ncomps + c);
        const auto *lf_276 = buffer.data(lf + 276 * ncomps + c);
        const auto *lf_277 = buffer.data(lf + 277 * ncomps + c);
        const auto *lf_278 = buffer.data(lf + 278 * ncomps + c);
        const auto *lf_279 = buffer.data(lf + 279 * ncomps + c);
        const auto *lf_280 = buffer.data(lf + 280 * ncomps + c);
        const auto *lf_281 = buffer.data(lf + 281 * ncomps + c);
        const auto *lf_282 = buffer.data(lf + 282 * ncomps + c);
        const auto *lf_283 = buffer.data(lf + 283 * ncomps + c);
        const auto *lf_284 = buffer.data(lf + 284 * ncomps + c);
        const auto *lf_285 = buffer.data(lf + 285 * ncomps + c);
        const auto *lf_286 = buffer.data(lf + 286 * ncomps + c);
        const auto *lf_287 = buffer.data(lf + 287 * ncomps + c);
        const auto *lf_288 = buffer.data(lf + 288 * ncomps + c);
        const auto *lf_289 = buffer.data(lf + 289 * ncomps + c);
        const auto *lf_296 = buffer.data(lf + 296 * ncomps + c);
        const auto *lf_297 = buffer.data(lf + 297 * ncomps + c);
        const auto *lf_298 = buffer.data(lf + 298 * ncomps + c);
        const auto *lf_299 = buffer.data(lf + 299 * ncomps + c);
        const auto *lf_306 = buffer.data(lf + 306 * ncomps + c);
        const auto *lf_307 = buffer.data(lf + 307 * ncomps + c);
        const auto *lf_308 = buffer.data(lf + 308 * ncomps + c);
        const auto *lf_309 = buffer.data(lf + 309 * ncomps + c);
        const auto *lf_316 = buffer.data(lf + 316 * ncomps + c);
        const auto *lf_317 = buffer.data(lf + 317 * ncomps + c);
        const auto *lf_318 = buffer.data(lf + 318 * ncomps + c);
        const auto *lf_319 = buffer.data(lf + 319 * ncomps + c);
        const auto *lf_326 = buffer.data(lf + 326 * ncomps + c);
        const auto *lf_327 = buffer.data(lf + 327 * ncomps + c);
        const auto *lf_328 = buffer.data(lf + 328 * ncomps + c);
        const auto *lf_329 = buffer.data(lf + 329 * ncomps + c);
        const auto *lf_336 = buffer.data(lf + 336 * ncomps + c);
        const auto *lf_337 = buffer.data(lf + 337 * ncomps + c);
        const auto *lf_338 = buffer.data(lf + 338 * ncomps + c);
        const auto *lf_339 = buffer.data(lf + 339 * ncomps + c);
        const auto *lf_346 = buffer.data(lf + 346 * ncomps + c);
        const auto *lf_347 = buffer.data(lf + 347 * ncomps + c);
        const auto *lf_348 = buffer.data(lf + 348 * ncomps + c);
        const auto *lf_349 = buffer.data(lf + 349 * ncomps + c);
        const auto *lf_359 = buffer.data(lf + 359 * ncomps + c);
        const auto *lf_366 = buffer.data(lf + 366 * ncomps + c);
        const auto *lf_367 = buffer.data(lf + 367 * ncomps + c);
        const auto *lf_368 = buffer.data(lf + 368 * ncomps + c);
        const auto *lf_369 = buffer.data(lf + 369 * ncomps + c);
        const auto *lf_379 = buffer.data(lf + 379 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, kf_195, kf_196, kf_197, \
                         kf_198, kf_199, lf_195, lf_196, lf_197, lf_198, \
                         lf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * kf_195[k]
                       + lf_195[k];

            t_291[k] = ab_x[k] * kf_196[k]
                       + lf_196[k];

            t_292[k] = ab_x[k] * kf_197[k]
                       + lf_197[k];

            t_293[k] = ab_x[k] * kf_198[k]
                       + lf_198[k];

            t_294[k] = ab_x[k] * kf_199[k]
                       + lf_199[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, ab_z, kf_196, kf_197, \
                         kf_198, kf_199, lf_256, lf_257, lf_258, lf_259, \
                         lf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_y[k] * kf_196[k]
                       + lf_256[k];

            t_296[k] = ab_y[k] * kf_197[k]
                       + lf_257[k];

            t_297[k] = ab_y[k] * kf_198[k]
                       + lf_258[k];

            t_298[k] = ab_y[k] * kf_199[k]
                       + lf_259[k];

            t_299[k] = ab_z[k] * kf_199[k]
                       + lf_269[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, kf_200, kf_201, kf_202, \
                         kf_203, kf_204, lf_200, lf_201, lf_202, lf_203, \
                         lf_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * kf_200[k]
                       + lf_200[k];

            t_301[k] = ab_x[k] * kf_201[k]
                       + lf_201[k];

            t_302[k] = ab_x[k] * kf_202[k]
                       + lf_202[k];

            t_303[k] = ab_x[k] * kf_203[k]
                       + lf_203[k];

            t_304[k] = ab_x[k] * kf_204[k]
                       + lf_204[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, kf_205, kf_206, kf_207, \
                         kf_208, kf_209, lf_205, lf_206, lf_207, lf_208, \
                         lf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * kf_205[k]
                       + lf_205[k];

            t_306[k] = ab_x[k] * kf_206[k]
                       + lf_206[k];

            t_307[k] = ab_x[k] * kf_207[k]
                       + lf_207[k];

            t_308[k] = ab_x[k] * kf_208[k]
                       + lf_208[k];

            t_309[k] = ab_x[k] * kf_209[k]
                       + lf_209[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, ab_z, kf_206, kf_207, \
                         kf_208, kf_209, lf_266, lf_267, lf_268, lf_269, \
                         lf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_y[k] * kf_206[k]
                       + lf_266[k];

            t_311[k] = ab_y[k] * kf_207[k]
                       + lf_267[k];

            t_312[k] = ab_y[k] * kf_208[k]
                       + lf_268[k];

            t_313[k] = ab_y[k] * kf_209[k]
                       + lf_269[k];

            t_314[k] = ab_z[k] * kf_209[k]
                       + lf_279[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, kf_210, kf_211, kf_212, \
                         kf_213, kf_214, lf_210, lf_211, lf_212, lf_213, \
                         lf_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * kf_210[k]
                       + lf_210[k];

            t_316[k] = ab_x[k] * kf_211[k]
                       + lf_211[k];

            t_317[k] = ab_x[k] * kf_212[k]
                       + lf_212[k];

            t_318[k] = ab_x[k] * kf_213[k]
                       + lf_213[k];

            t_319[k] = ab_x[k] * kf_214[k]
                       + lf_214[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, kf_215, kf_216, kf_217, \
                         kf_218, kf_219, lf_215, lf_216, lf_217, lf_218, \
                         lf_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * kf_215[k]
                       + lf_215[k];

            t_321[k] = ab_x[k] * kf_216[k]
                       + lf_216[k];

            t_322[k] = ab_x[k] * kf_217[k]
                       + lf_217[k];

            t_323[k] = ab_x[k] * kf_218[k]
                       + lf_218[k];

            t_324[k] = ab_x[k] * kf_219[k]
                       + lf_219[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, ab_z, kf_216, kf_217, \
                         kf_218, kf_219, lf_286, lf_287, lf_288, lf_289, \
                         lf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_y[k] * kf_216[k]
                       + lf_286[k];

            t_326[k] = ab_y[k] * kf_217[k]
                       + lf_287[k];

            t_327[k] = ab_y[k] * kf_218[k]
                       + lf_288[k];

            t_328[k] = ab_y[k] * kf_219[k]
                       + lf_289[k];

            t_329[k] = ab_z[k] * kf_219[k]
                       + lf_299[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, kf_220, kf_221, kf_222, \
                         kf_223, kf_224, lf_220, lf_221, lf_222, lf_223, \
                         lf_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * kf_220[k]
                       + lf_220[k];

            t_331[k] = ab_x[k] * kf_221[k]
                       + lf_221[k];

            t_332[k] = ab_x[k] * kf_222[k]
                       + lf_222[k];

            t_333[k] = ab_x[k] * kf_223[k]
                       + lf_223[k];

            t_334[k] = ab_x[k] * kf_224[k]
                       + lf_224[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, kf_225, kf_226, kf_227, \
                         kf_228, kf_229, lf_225, lf_226, lf_227, lf_228, \
                         lf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * kf_225[k]
                       + lf_225[k];

            t_336[k] = ab_x[k] * kf_226[k]
                       + lf_226[k];

            t_337[k] = ab_x[k] * kf_227[k]
                       + lf_227[k];

            t_338[k] = ab_x[k] * kf_228[k]
                       + lf_228[k];

            t_339[k] = ab_x[k] * kf_229[k]
                       + lf_229[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, ab_z, kf_226, kf_227, \
                         kf_228, kf_229, lf_296, lf_297, lf_298, lf_299, \
                         lf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_y[k] * kf_226[k]
                       + lf_296[k];

            t_341[k] = ab_y[k] * kf_227[k]
                       + lf_297[k];

            t_342[k] = ab_y[k] * kf_228[k]
                       + lf_298[k];

            t_343[k] = ab_y[k] * kf_229[k]
                       + lf_299[k];

            t_344[k] = ab_z[k] * kf_229[k]
                       + lf_309[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, kf_230, kf_231, kf_232, \
                         kf_233, kf_234, lf_230, lf_231, lf_232, lf_233, \
                         lf_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * kf_230[k]
                       + lf_230[k];

            t_346[k] = ab_x[k] * kf_231[k]
                       + lf_231[k];

            t_347[k] = ab_x[k] * kf_232[k]
                       + lf_232[k];

            t_348[k] = ab_x[k] * kf_233[k]
                       + lf_233[k];

            t_349[k] = ab_x[k] * kf_234[k]
                       + lf_234[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, kf_235, kf_236, kf_237, \
                         kf_238, kf_239, lf_235, lf_236, lf_237, lf_238, \
                         lf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * kf_235[k]
                       + lf_235[k];

            t_351[k] = ab_x[k] * kf_236[k]
                       + lf_236[k];

            t_352[k] = ab_x[k] * kf_237[k]
                       + lf_237[k];

            t_353[k] = ab_x[k] * kf_238[k]
                       + lf_238[k];

            t_354[k] = ab_x[k] * kf_239[k]
                       + lf_239[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, ab_z, kf_236, kf_237, \
                         kf_238, kf_239, lf_306, lf_307, lf_308, lf_309, \
                         lf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_y[k] * kf_236[k]
                       + lf_306[k];

            t_356[k] = ab_y[k] * kf_237[k]
                       + lf_307[k];

            t_357[k] = ab_y[k] * kf_238[k]
                       + lf_308[k];

            t_358[k] = ab_y[k] * kf_239[k]
                       + lf_309[k];

            t_359[k] = ab_z[k] * kf_239[k]
                       + lf_319[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, kf_240, kf_241, kf_242, \
                         kf_243, kf_244, lf_240, lf_241, lf_242, lf_243, \
                         lf_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * kf_240[k]
                       + lf_240[k];

            t_361[k] = ab_x[k] * kf_241[k]
                       + lf_241[k];

            t_362[k] = ab_x[k] * kf_242[k]
                       + lf_242[k];

            t_363[k] = ab_x[k] * kf_243[k]
                       + lf_243[k];

            t_364[k] = ab_x[k] * kf_244[k]
                       + lf_244[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, kf_245, kf_246, kf_247, \
                         kf_248, kf_249, lf_245, lf_246, lf_247, lf_248, \
                         lf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * kf_245[k]
                       + lf_245[k];

            t_366[k] = ab_x[k] * kf_246[k]
                       + lf_246[k];

            t_367[k] = ab_x[k] * kf_247[k]
                       + lf_247[k];

            t_368[k] = ab_x[k] * kf_248[k]
                       + lf_248[k];

            t_369[k] = ab_x[k] * kf_249[k]
                       + lf_249[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, ab_z, kf_246, kf_247, \
                         kf_248, kf_249, lf_316, lf_317, lf_318, lf_319, \
                         lf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_y[k] * kf_246[k]
                       + lf_316[k];

            t_371[k] = ab_y[k] * kf_247[k]
                       + lf_317[k];

            t_372[k] = ab_y[k] * kf_248[k]
                       + lf_318[k];

            t_373[k] = ab_y[k] * kf_249[k]
                       + lf_319[k];

            t_374[k] = ab_z[k] * kf_249[k]
                       + lf_329[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, kf_250, kf_251, kf_252, \
                         kf_253, kf_254, lf_250, lf_251, lf_252, lf_253, \
                         lf_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * kf_250[k]
                       + lf_250[k];

            t_376[k] = ab_x[k] * kf_251[k]
                       + lf_251[k];

            t_377[k] = ab_x[k] * kf_252[k]
                       + lf_252[k];

            t_378[k] = ab_x[k] * kf_253[k]
                       + lf_253[k];

            t_379[k] = ab_x[k] * kf_254[k]
                       + lf_254[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, kf_255, kf_256, kf_257, \
                         kf_258, kf_259, lf_255, lf_256, lf_257, lf_258, \
                         lf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * kf_255[k]
                       + lf_255[k];

            t_381[k] = ab_x[k] * kf_256[k]
                       + lf_256[k];

            t_382[k] = ab_x[k] * kf_257[k]
                       + lf_257[k];

            t_383[k] = ab_x[k] * kf_258[k]
                       + lf_258[k];

            t_384[k] = ab_x[k] * kf_259[k]
                       + lf_259[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, ab_z, kf_256, kf_257, \
                         kf_258, kf_259, lf_326, lf_327, lf_328, lf_329, \
                         lf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_y[k] * kf_256[k]
                       + lf_326[k];

            t_386[k] = ab_y[k] * kf_257[k]
                       + lf_327[k];

            t_387[k] = ab_y[k] * kf_258[k]
                       + lf_328[k];

            t_388[k] = ab_y[k] * kf_259[k]
                       + lf_329[k];

            t_389[k] = ab_z[k] * kf_259[k]
                       + lf_339[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, kf_260, kf_261, kf_262, \
                         kf_263, kf_264, lf_260, lf_261, lf_262, lf_263, \
                         lf_264 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * kf_260[k]
                       + lf_260[k];

            t_391[k] = ab_x[k] * kf_261[k]
                       + lf_261[k];

            t_392[k] = ab_x[k] * kf_262[k]
                       + lf_262[k];

            t_393[k] = ab_x[k] * kf_263[k]
                       + lf_263[k];

            t_394[k] = ab_x[k] * kf_264[k]
                       + lf_264[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, kf_265, kf_266, kf_267, \
                         kf_268, kf_269, lf_265, lf_266, lf_267, lf_268, \
                         lf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * kf_265[k]
                       + lf_265[k];

            t_396[k] = ab_x[k] * kf_266[k]
                       + lf_266[k];

            t_397[k] = ab_x[k] * kf_267[k]
                       + lf_267[k];

            t_398[k] = ab_x[k] * kf_268[k]
                       + lf_268[k];

            t_399[k] = ab_x[k] * kf_269[k]
                       + lf_269[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, ab_z, kf_266, kf_267, \
                         kf_268, kf_269, lf_336, lf_337, lf_338, lf_339, \
                         lf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_y[k] * kf_266[k]
                       + lf_336[k];

            t_401[k] = ab_y[k] * kf_267[k]
                       + lf_337[k];

            t_402[k] = ab_y[k] * kf_268[k]
                       + lf_338[k];

            t_403[k] = ab_y[k] * kf_269[k]
                       + lf_339[k];

            t_404[k] = ab_z[k] * kf_269[k]
                       + lf_349[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, kf_270, kf_271, kf_272, \
                         kf_273, kf_274, lf_270, lf_271, lf_272, lf_273, \
                         lf_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * kf_270[k]
                       + lf_270[k];

            t_406[k] = ab_x[k] * kf_271[k]
                       + lf_271[k];

            t_407[k] = ab_x[k] * kf_272[k]
                       + lf_272[k];

            t_408[k] = ab_x[k] * kf_273[k]
                       + lf_273[k];

            t_409[k] = ab_x[k] * kf_274[k]
                       + lf_274[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, kf_275, kf_276, kf_277, \
                         kf_278, kf_279, lf_275, lf_276, lf_277, lf_278, \
                         lf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * kf_275[k]
                       + lf_275[k];

            t_411[k] = ab_x[k] * kf_276[k]
                       + lf_276[k];

            t_412[k] = ab_x[k] * kf_277[k]
                       + lf_277[k];

            t_413[k] = ab_x[k] * kf_278[k]
                       + lf_278[k];

            t_414[k] = ab_x[k] * kf_279[k]
                       + lf_279[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, ab_z, kf_276, kf_277, \
                         kf_278, kf_279, lf_346, lf_347, lf_348, lf_349, \
                         lf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_y[k] * kf_276[k]
                       + lf_346[k];

            t_416[k] = ab_y[k] * kf_277[k]
                       + lf_347[k];

            t_417[k] = ab_y[k] * kf_278[k]
                       + lf_348[k];

            t_418[k] = ab_y[k] * kf_279[k]
                       + lf_349[k];

            t_419[k] = ab_z[k] * kf_279[k]
                       + lf_359[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, kf_280, kf_281, kf_282, \
                         kf_283, kf_284, lf_280, lf_281, lf_282, lf_283, \
                         lf_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * kf_280[k]
                       + lf_280[k];

            t_421[k] = ab_x[k] * kf_281[k]
                       + lf_281[k];

            t_422[k] = ab_x[k] * kf_282[k]
                       + lf_282[k];

            t_423[k] = ab_x[k] * kf_283[k]
                       + lf_283[k];

            t_424[k] = ab_x[k] * kf_284[k]
                       + lf_284[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, kf_285, kf_286, kf_287, \
                         kf_288, kf_289, lf_285, lf_286, lf_287, lf_288, \
                         lf_289 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * kf_285[k]
                       + lf_285[k];

            t_426[k] = ab_x[k] * kf_286[k]
                       + lf_286[k];

            t_427[k] = ab_x[k] * kf_287[k]
                       + lf_287[k];

            t_428[k] = ab_x[k] * kf_288[k]
                       + lf_288[k];

            t_429[k] = ab_x[k] * kf_289[k]
                       + lf_289[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_y, ab_z, kf_286, kf_287, \
                         kf_288, kf_289, lf_366, lf_367, lf_368, lf_369, \
                         lf_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_y[k] * kf_286[k]
                       + lf_366[k];

            t_431[k] = ab_y[k] * kf_287[k]
                       + lf_367[k];

            t_432[k] = ab_y[k] * kf_288[k]
                       + lf_368[k];

            t_433[k] = ab_y[k] * kf_289[k]
                       + lf_369[k];

            t_434[k] = ab_z[k] * kf_289[k]
                       + lf_379[k];
        }
    }
}

static auto
compute_hrr_kg_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kf, const size_t lf,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kf_290 = buffer.data(kf + 290 * ncomps + c);
        const auto *kf_291 = buffer.data(kf + 291 * ncomps + c);
        const auto *kf_292 = buffer.data(kf + 292 * ncomps + c);
        const auto *kf_293 = buffer.data(kf + 293 * ncomps + c);
        const auto *kf_294 = buffer.data(kf + 294 * ncomps + c);
        const auto *kf_295 = buffer.data(kf + 295 * ncomps + c);
        const auto *kf_296 = buffer.data(kf + 296 * ncomps + c);
        const auto *kf_297 = buffer.data(kf + 297 * ncomps + c);
        const auto *kf_298 = buffer.data(kf + 298 * ncomps + c);
        const auto *kf_299 = buffer.data(kf + 299 * ncomps + c);
        const auto *kf_300 = buffer.data(kf + 300 * ncomps + c);
        const auto *kf_301 = buffer.data(kf + 301 * ncomps + c);
        const auto *kf_302 = buffer.data(kf + 302 * ncomps + c);
        const auto *kf_303 = buffer.data(kf + 303 * ncomps + c);
        const auto *kf_304 = buffer.data(kf + 304 * ncomps + c);
        const auto *kf_305 = buffer.data(kf + 305 * ncomps + c);
        const auto *kf_306 = buffer.data(kf + 306 * ncomps + c);
        const auto *kf_307 = buffer.data(kf + 307 * ncomps + c);
        const auto *kf_308 = buffer.data(kf + 308 * ncomps + c);
        const auto *kf_309 = buffer.data(kf + 309 * ncomps + c);
        const auto *kf_310 = buffer.data(kf + 310 * ncomps + c);
        const auto *kf_311 = buffer.data(kf + 311 * ncomps + c);
        const auto *kf_312 = buffer.data(kf + 312 * ncomps + c);
        const auto *kf_313 = buffer.data(kf + 313 * ncomps + c);
        const auto *kf_314 = buffer.data(kf + 314 * ncomps + c);
        const auto *kf_315 = buffer.data(kf + 315 * ncomps + c);
        const auto *kf_316 = buffer.data(kf + 316 * ncomps + c);
        const auto *kf_317 = buffer.data(kf + 317 * ncomps + c);
        const auto *kf_318 = buffer.data(kf + 318 * ncomps + c);
        const auto *kf_319 = buffer.data(kf + 319 * ncomps + c);
        const auto *kf_320 = buffer.data(kf + 320 * ncomps + c);
        const auto *kf_321 = buffer.data(kf + 321 * ncomps + c);
        const auto *kf_322 = buffer.data(kf + 322 * ncomps + c);
        const auto *kf_323 = buffer.data(kf + 323 * ncomps + c);
        const auto *kf_324 = buffer.data(kf + 324 * ncomps + c);
        const auto *kf_325 = buffer.data(kf + 325 * ncomps + c);
        const auto *kf_326 = buffer.data(kf + 326 * ncomps + c);
        const auto *kf_327 = buffer.data(kf + 327 * ncomps + c);
        const auto *kf_328 = buffer.data(kf + 328 * ncomps + c);
        const auto *kf_329 = buffer.data(kf + 329 * ncomps + c);
        const auto *kf_330 = buffer.data(kf + 330 * ncomps + c);
        const auto *kf_331 = buffer.data(kf + 331 * ncomps + c);
        const auto *kf_332 = buffer.data(kf + 332 * ncomps + c);
        const auto *kf_333 = buffer.data(kf + 333 * ncomps + c);
        const auto *kf_334 = buffer.data(kf + 334 * ncomps + c);
        const auto *kf_335 = buffer.data(kf + 335 * ncomps + c);
        const auto *kf_336 = buffer.data(kf + 336 * ncomps + c);
        const auto *kf_337 = buffer.data(kf + 337 * ncomps + c);
        const auto *kf_338 = buffer.data(kf + 338 * ncomps + c);
        const auto *kf_339 = buffer.data(kf + 339 * ncomps + c);
        const auto *kf_340 = buffer.data(kf + 340 * ncomps + c);
        const auto *kf_341 = buffer.data(kf + 341 * ncomps + c);
        const auto *kf_342 = buffer.data(kf + 342 * ncomps + c);
        const auto *kf_343 = buffer.data(kf + 343 * ncomps + c);
        const auto *kf_344 = buffer.data(kf + 344 * ncomps + c);
        const auto *kf_345 = buffer.data(kf + 345 * ncomps + c);
        const auto *kf_346 = buffer.data(kf + 346 * ncomps + c);
        const auto *kf_347 = buffer.data(kf + 347 * ncomps + c);
        const auto *kf_348 = buffer.data(kf + 348 * ncomps + c);
        const auto *kf_349 = buffer.data(kf + 349 * ncomps + c);
        const auto *kf_350 = buffer.data(kf + 350 * ncomps + c);
        const auto *kf_351 = buffer.data(kf + 351 * ncomps + c);
        const auto *kf_352 = buffer.data(kf + 352 * ncomps + c);
        const auto *kf_353 = buffer.data(kf + 353 * ncomps + c);
        const auto *kf_354 = buffer.data(kf + 354 * ncomps + c);
        const auto *kf_355 = buffer.data(kf + 355 * ncomps + c);
        const auto *kf_356 = buffer.data(kf + 356 * ncomps + c);
        const auto *kf_357 = buffer.data(kf + 357 * ncomps + c);
        const auto *kf_358 = buffer.data(kf + 358 * ncomps + c);
        const auto *kf_359 = buffer.data(kf + 359 * ncomps + c);

        const auto *lf_290 = buffer.data(lf + 290 * ncomps + c);
        const auto *lf_291 = buffer.data(lf + 291 * ncomps + c);
        const auto *lf_292 = buffer.data(lf + 292 * ncomps + c);
        const auto *lf_293 = buffer.data(lf + 293 * ncomps + c);
        const auto *lf_294 = buffer.data(lf + 294 * ncomps + c);
        const auto *lf_295 = buffer.data(lf + 295 * ncomps + c);
        const auto *lf_296 = buffer.data(lf + 296 * ncomps + c);
        const auto *lf_297 = buffer.data(lf + 297 * ncomps + c);
        const auto *lf_298 = buffer.data(lf + 298 * ncomps + c);
        const auto *lf_299 = buffer.data(lf + 299 * ncomps + c);
        const auto *lf_300 = buffer.data(lf + 300 * ncomps + c);
        const auto *lf_301 = buffer.data(lf + 301 * ncomps + c);
        const auto *lf_302 = buffer.data(lf + 302 * ncomps + c);
        const auto *lf_303 = buffer.data(lf + 303 * ncomps + c);
        const auto *lf_304 = buffer.data(lf + 304 * ncomps + c);
        const auto *lf_305 = buffer.data(lf + 305 * ncomps + c);
        const auto *lf_306 = buffer.data(lf + 306 * ncomps + c);
        const auto *lf_307 = buffer.data(lf + 307 * ncomps + c);
        const auto *lf_308 = buffer.data(lf + 308 * ncomps + c);
        const auto *lf_309 = buffer.data(lf + 309 * ncomps + c);
        const auto *lf_310 = buffer.data(lf + 310 * ncomps + c);
        const auto *lf_311 = buffer.data(lf + 311 * ncomps + c);
        const auto *lf_312 = buffer.data(lf + 312 * ncomps + c);
        const auto *lf_313 = buffer.data(lf + 313 * ncomps + c);
        const auto *lf_314 = buffer.data(lf + 314 * ncomps + c);
        const auto *lf_315 = buffer.data(lf + 315 * ncomps + c);
        const auto *lf_316 = buffer.data(lf + 316 * ncomps + c);
        const auto *lf_317 = buffer.data(lf + 317 * ncomps + c);
        const auto *lf_318 = buffer.data(lf + 318 * ncomps + c);
        const auto *lf_319 = buffer.data(lf + 319 * ncomps + c);
        const auto *lf_320 = buffer.data(lf + 320 * ncomps + c);
        const auto *lf_321 = buffer.data(lf + 321 * ncomps + c);
        const auto *lf_322 = buffer.data(lf + 322 * ncomps + c);
        const auto *lf_323 = buffer.data(lf + 323 * ncomps + c);
        const auto *lf_324 = buffer.data(lf + 324 * ncomps + c);
        const auto *lf_325 = buffer.data(lf + 325 * ncomps + c);
        const auto *lf_326 = buffer.data(lf + 326 * ncomps + c);
        const auto *lf_327 = buffer.data(lf + 327 * ncomps + c);
        const auto *lf_328 = buffer.data(lf + 328 * ncomps + c);
        const auto *lf_329 = buffer.data(lf + 329 * ncomps + c);
        const auto *lf_330 = buffer.data(lf + 330 * ncomps + c);
        const auto *lf_331 = buffer.data(lf + 331 * ncomps + c);
        const auto *lf_332 = buffer.data(lf + 332 * ncomps + c);
        const auto *lf_333 = buffer.data(lf + 333 * ncomps + c);
        const auto *lf_334 = buffer.data(lf + 334 * ncomps + c);
        const auto *lf_335 = buffer.data(lf + 335 * ncomps + c);
        const auto *lf_336 = buffer.data(lf + 336 * ncomps + c);
        const auto *lf_337 = buffer.data(lf + 337 * ncomps + c);
        const auto *lf_338 = buffer.data(lf + 338 * ncomps + c);
        const auto *lf_339 = buffer.data(lf + 339 * ncomps + c);
        const auto *lf_340 = buffer.data(lf + 340 * ncomps + c);
        const auto *lf_341 = buffer.data(lf + 341 * ncomps + c);
        const auto *lf_342 = buffer.data(lf + 342 * ncomps + c);
        const auto *lf_343 = buffer.data(lf + 343 * ncomps + c);
        const auto *lf_344 = buffer.data(lf + 344 * ncomps + c);
        const auto *lf_345 = buffer.data(lf + 345 * ncomps + c);
        const auto *lf_346 = buffer.data(lf + 346 * ncomps + c);
        const auto *lf_347 = buffer.data(lf + 347 * ncomps + c);
        const auto *lf_348 = buffer.data(lf + 348 * ncomps + c);
        const auto *lf_349 = buffer.data(lf + 349 * ncomps + c);
        const auto *lf_350 = buffer.data(lf + 350 * ncomps + c);
        const auto *lf_351 = buffer.data(lf + 351 * ncomps + c);
        const auto *lf_352 = buffer.data(lf + 352 * ncomps + c);
        const auto *lf_353 = buffer.data(lf + 353 * ncomps + c);
        const auto *lf_354 = buffer.data(lf + 354 * ncomps + c);
        const auto *lf_355 = buffer.data(lf + 355 * ncomps + c);
        const auto *lf_356 = buffer.data(lf + 356 * ncomps + c);
        const auto *lf_357 = buffer.data(lf + 357 * ncomps + c);
        const auto *lf_358 = buffer.data(lf + 358 * ncomps + c);
        const auto *lf_359 = buffer.data(lf + 359 * ncomps + c);
        const auto *lf_376 = buffer.data(lf + 376 * ncomps + c);
        const auto *lf_377 = buffer.data(lf + 377 * ncomps + c);
        const auto *lf_378 = buffer.data(lf + 378 * ncomps + c);
        const auto *lf_379 = buffer.data(lf + 379 * ncomps + c);
        const auto *lf_386 = buffer.data(lf + 386 * ncomps + c);
        const auto *lf_387 = buffer.data(lf + 387 * ncomps + c);
        const auto *lf_388 = buffer.data(lf + 388 * ncomps + c);
        const auto *lf_389 = buffer.data(lf + 389 * ncomps + c);
        const auto *lf_396 = buffer.data(lf + 396 * ncomps + c);
        const auto *lf_397 = buffer.data(lf + 397 * ncomps + c);
        const auto *lf_398 = buffer.data(lf + 398 * ncomps + c);
        const auto *lf_399 = buffer.data(lf + 399 * ncomps + c);
        const auto *lf_406 = buffer.data(lf + 406 * ncomps + c);
        const auto *lf_407 = buffer.data(lf + 407 * ncomps + c);
        const auto *lf_408 = buffer.data(lf + 408 * ncomps + c);
        const auto *lf_409 = buffer.data(lf + 409 * ncomps + c);
        const auto *lf_416 = buffer.data(lf + 416 * ncomps + c);
        const auto *lf_417 = buffer.data(lf + 417 * ncomps + c);
        const auto *lf_418 = buffer.data(lf + 418 * ncomps + c);
        const auto *lf_419 = buffer.data(lf + 419 * ncomps + c);
        const auto *lf_426 = buffer.data(lf + 426 * ncomps + c);
        const auto *lf_427 = buffer.data(lf + 427 * ncomps + c);
        const auto *lf_428 = buffer.data(lf + 428 * ncomps + c);
        const auto *lf_429 = buffer.data(lf + 429 * ncomps + c);
        const auto *lf_436 = buffer.data(lf + 436 * ncomps + c);
        const auto *lf_437 = buffer.data(lf + 437 * ncomps + c);
        const auto *lf_438 = buffer.data(lf + 438 * ncomps + c);
        const auto *lf_439 = buffer.data(lf + 439 * ncomps + c);
        const auto *lf_449 = buffer.data(lf + 449 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, kf_290, kf_291, kf_292, \
                         kf_293, kf_294, lf_290, lf_291, lf_292, lf_293, \
                         lf_294 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * kf_290[k]
                       + lf_290[k];

            t_436[k] = ab_x[k] * kf_291[k]
                       + lf_291[k];

            t_437[k] = ab_x[k] * kf_292[k]
                       + lf_292[k];

            t_438[k] = ab_x[k] * kf_293[k]
                       + lf_293[k];

            t_439[k] = ab_x[k] * kf_294[k]
                       + lf_294[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, kf_295, kf_296, kf_297, \
                         kf_298, kf_299, lf_295, lf_296, lf_297, lf_298, \
                         lf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * kf_295[k]
                       + lf_295[k];

            t_441[k] = ab_x[k] * kf_296[k]
                       + lf_296[k];

            t_442[k] = ab_x[k] * kf_297[k]
                       + lf_297[k];

            t_443[k] = ab_x[k] * kf_298[k]
                       + lf_298[k];

            t_444[k] = ab_x[k] * kf_299[k]
                       + lf_299[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_y, ab_z, kf_296, kf_297, \
                         kf_298, kf_299, lf_376, lf_377, lf_378, lf_379, \
                         lf_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = ab_y[k] * kf_296[k]
                       + lf_376[k];

            t_446[k] = ab_y[k] * kf_297[k]
                       + lf_377[k];

            t_447[k] = ab_y[k] * kf_298[k]
                       + lf_378[k];

            t_448[k] = ab_y[k] * kf_299[k]
                       + lf_379[k];

            t_449[k] = ab_z[k] * kf_299[k]
                       + lf_389[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_x, kf_300, kf_301, kf_302, \
                         kf_303, kf_304, lf_300, lf_301, lf_302, lf_303, \
                         lf_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = ab_x[k] * kf_300[k]
                       + lf_300[k];

            t_451[k] = ab_x[k] * kf_301[k]
                       + lf_301[k];

            t_452[k] = ab_x[k] * kf_302[k]
                       + lf_302[k];

            t_453[k] = ab_x[k] * kf_303[k]
                       + lf_303[k];

            t_454[k] = ab_x[k] * kf_304[k]
                       + lf_304[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_x, kf_305, kf_306, kf_307, \
                         kf_308, kf_309, lf_305, lf_306, lf_307, lf_308, \
                         lf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = ab_x[k] * kf_305[k]
                       + lf_305[k];

            t_456[k] = ab_x[k] * kf_306[k]
                       + lf_306[k];

            t_457[k] = ab_x[k] * kf_307[k]
                       + lf_307[k];

            t_458[k] = ab_x[k] * kf_308[k]
                       + lf_308[k];

            t_459[k] = ab_x[k] * kf_309[k]
                       + lf_309[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, ab_z, kf_306, kf_307, \
                         kf_308, kf_309, lf_386, lf_387, lf_388, lf_389, \
                         lf_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = ab_y[k] * kf_306[k]
                       + lf_386[k];

            t_461[k] = ab_y[k] * kf_307[k]
                       + lf_387[k];

            t_462[k] = ab_y[k] * kf_308[k]
                       + lf_388[k];

            t_463[k] = ab_y[k] * kf_309[k]
                       + lf_389[k];

            t_464[k] = ab_z[k] * kf_309[k]
                       + lf_399[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_x, kf_310, kf_311, kf_312, \
                         kf_313, kf_314, lf_310, lf_311, lf_312, lf_313, \
                         lf_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = ab_x[k] * kf_310[k]
                       + lf_310[k];

            t_466[k] = ab_x[k] * kf_311[k]
                       + lf_311[k];

            t_467[k] = ab_x[k] * kf_312[k]
                       + lf_312[k];

            t_468[k] = ab_x[k] * kf_313[k]
                       + lf_313[k];

            t_469[k] = ab_x[k] * kf_314[k]
                       + lf_314[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_x, kf_315, kf_316, kf_317, \
                         kf_318, kf_319, lf_315, lf_316, lf_317, lf_318, \
                         lf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = ab_x[k] * kf_315[k]
                       + lf_315[k];

            t_471[k] = ab_x[k] * kf_316[k]
                       + lf_316[k];

            t_472[k] = ab_x[k] * kf_317[k]
                       + lf_317[k];

            t_473[k] = ab_x[k] * kf_318[k]
                       + lf_318[k];

            t_474[k] = ab_x[k] * kf_319[k]
                       + lf_319[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, ab_z, kf_316, kf_317, \
                         kf_318, kf_319, lf_396, lf_397, lf_398, lf_399, \
                         lf_409 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = ab_y[k] * kf_316[k]
                       + lf_396[k];

            t_476[k] = ab_y[k] * kf_317[k]
                       + lf_397[k];

            t_477[k] = ab_y[k] * kf_318[k]
                       + lf_398[k];

            t_478[k] = ab_y[k] * kf_319[k]
                       + lf_399[k];

            t_479[k] = ab_z[k] * kf_319[k]
                       + lf_409[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_x, kf_320, kf_321, kf_322, \
                         kf_323, kf_324, lf_320, lf_321, lf_322, lf_323, \
                         lf_324 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = ab_x[k] * kf_320[k]
                       + lf_320[k];

            t_481[k] = ab_x[k] * kf_321[k]
                       + lf_321[k];

            t_482[k] = ab_x[k] * kf_322[k]
                       + lf_322[k];

            t_483[k] = ab_x[k] * kf_323[k]
                       + lf_323[k];

            t_484[k] = ab_x[k] * kf_324[k]
                       + lf_324[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_x, kf_325, kf_326, kf_327, \
                         kf_328, kf_329, lf_325, lf_326, lf_327, lf_328, \
                         lf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = ab_x[k] * kf_325[k]
                       + lf_325[k];

            t_486[k] = ab_x[k] * kf_326[k]
                       + lf_326[k];

            t_487[k] = ab_x[k] * kf_327[k]
                       + lf_327[k];

            t_488[k] = ab_x[k] * kf_328[k]
                       + lf_328[k];

            t_489[k] = ab_x[k] * kf_329[k]
                       + lf_329[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, ab_z, kf_326, kf_327, \
                         kf_328, kf_329, lf_406, lf_407, lf_408, lf_409, \
                         lf_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = ab_y[k] * kf_326[k]
                       + lf_406[k];

            t_491[k] = ab_y[k] * kf_327[k]
                       + lf_407[k];

            t_492[k] = ab_y[k] * kf_328[k]
                       + lf_408[k];

            t_493[k] = ab_y[k] * kf_329[k]
                       + lf_409[k];

            t_494[k] = ab_z[k] * kf_329[k]
                       + lf_419[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_x, kf_330, kf_331, kf_332, \
                         kf_333, kf_334, lf_330, lf_331, lf_332, lf_333, \
                         lf_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = ab_x[k] * kf_330[k]
                       + lf_330[k];

            t_496[k] = ab_x[k] * kf_331[k]
                       + lf_331[k];

            t_497[k] = ab_x[k] * kf_332[k]
                       + lf_332[k];

            t_498[k] = ab_x[k] * kf_333[k]
                       + lf_333[k];

            t_499[k] = ab_x[k] * kf_334[k]
                       + lf_334[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_x, kf_335, kf_336, kf_337, \
                         kf_338, kf_339, lf_335, lf_336, lf_337, lf_338, \
                         lf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_x[k] * kf_335[k]
                       + lf_335[k];

            t_501[k] = ab_x[k] * kf_336[k]
                       + lf_336[k];

            t_502[k] = ab_x[k] * kf_337[k]
                       + lf_337[k];

            t_503[k] = ab_x[k] * kf_338[k]
                       + lf_338[k];

            t_504[k] = ab_x[k] * kf_339[k]
                       + lf_339[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_y, ab_z, kf_336, kf_337, \
                         kf_338, kf_339, lf_416, lf_417, lf_418, lf_419, \
                         lf_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = ab_y[k] * kf_336[k]
                       + lf_416[k];

            t_506[k] = ab_y[k] * kf_337[k]
                       + lf_417[k];

            t_507[k] = ab_y[k] * kf_338[k]
                       + lf_418[k];

            t_508[k] = ab_y[k] * kf_339[k]
                       + lf_419[k];

            t_509[k] = ab_z[k] * kf_339[k]
                       + lf_429[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_x, kf_340, kf_341, kf_342, \
                         kf_343, kf_344, lf_340, lf_341, lf_342, lf_343, \
                         lf_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = ab_x[k] * kf_340[k]
                       + lf_340[k];

            t_511[k] = ab_x[k] * kf_341[k]
                       + lf_341[k];

            t_512[k] = ab_x[k] * kf_342[k]
                       + lf_342[k];

            t_513[k] = ab_x[k] * kf_343[k]
                       + lf_343[k];

            t_514[k] = ab_x[k] * kf_344[k]
                       + lf_344[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_x, kf_345, kf_346, kf_347, \
                         kf_348, kf_349, lf_345, lf_346, lf_347, lf_348, \
                         lf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = ab_x[k] * kf_345[k]
                       + lf_345[k];

            t_516[k] = ab_x[k] * kf_346[k]
                       + lf_346[k];

            t_517[k] = ab_x[k] * kf_347[k]
                       + lf_347[k];

            t_518[k] = ab_x[k] * kf_348[k]
                       + lf_348[k];

            t_519[k] = ab_x[k] * kf_349[k]
                       + lf_349[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_y, ab_z, kf_346, kf_347, \
                         kf_348, kf_349, lf_426, lf_427, lf_428, lf_429, \
                         lf_439 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = ab_y[k] * kf_346[k]
                       + lf_426[k];

            t_521[k] = ab_y[k] * kf_347[k]
                       + lf_427[k];

            t_522[k] = ab_y[k] * kf_348[k]
                       + lf_428[k];

            t_523[k] = ab_y[k] * kf_349[k]
                       + lf_429[k];

            t_524[k] = ab_z[k] * kf_349[k]
                       + lf_439[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_x, kf_350, kf_351, kf_352, \
                         kf_353, kf_354, lf_350, lf_351, lf_352, lf_353, \
                         lf_354 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = ab_x[k] * kf_350[k]
                       + lf_350[k];

            t_526[k] = ab_x[k] * kf_351[k]
                       + lf_351[k];

            t_527[k] = ab_x[k] * kf_352[k]
                       + lf_352[k];

            t_528[k] = ab_x[k] * kf_353[k]
                       + lf_353[k];

            t_529[k] = ab_x[k] * kf_354[k]
                       + lf_354[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_x, kf_355, kf_356, kf_357, \
                         kf_358, kf_359, lf_355, lf_356, lf_357, lf_358, \
                         lf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = ab_x[k] * kf_355[k]
                       + lf_355[k];

            t_531[k] = ab_x[k] * kf_356[k]
                       + lf_356[k];

            t_532[k] = ab_x[k] * kf_357[k]
                       + lf_357[k];

            t_533[k] = ab_x[k] * kf_358[k]
                       + lf_358[k];

            t_534[k] = ab_x[k] * kf_359[k]
                       + lf_359[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_y, ab_z, kf_356, kf_357, \
                         kf_358, kf_359, lf_436, lf_437, lf_438, lf_439, \
                         lf_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = ab_y[k] * kf_356[k]
                       + lf_436[k];

            t_536[k] = ab_y[k] * kf_357[k]
                       + lf_437[k];

            t_537[k] = ab_y[k] * kf_358[k]
                       + lf_438[k];

            t_538[k] = ab_y[k] * kf_359[k]
                       + lf_439[k];

            t_539[k] = ab_z[k] * kf_359[k]
                       + lf_449[k];
        }
    }
}

auto
compute_hrr_kg_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t kf, const size_t lf,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_kg_out_of_first_piece0(buffer, coordinates, target, kf, lf, ncomps, nmax);

    compute_hrr_kg_out_of_first_piece1(buffer, coordinates, target, kf, lf, ncomps, nmax);

    compute_hrr_kg_out_of_first_piece2(buffer, coordinates, target, kf, lf, ncomps, nmax);

    compute_hrr_kg_out_of_first_piece3(buffer, coordinates, target, kf, lf, ncomps, nmax);
}

static auto
compute_hrr_kg_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t kf, const size_t lf, const size_t ncomps,
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
        const auto *ab_z = coordinates.data(8);

        const auto *kf_0 = buffer.data(kf + 0 * ncomps + c);
        const auto *kf_1 = buffer.data(kf + 1 * ncomps + c);
        const auto *kf_2 = buffer.data(kf + 2 * ncomps + c);
        const auto *kf_3 = buffer.data(kf + 3 * ncomps + c);
        const auto *kf_4 = buffer.data(kf + 4 * ncomps + c);
        const auto *kf_5 = buffer.data(kf + 5 * ncomps + c);
        const auto *kf_6 = buffer.data(kf + 6 * ncomps + c);
        const auto *kf_7 = buffer.data(kf + 7 * ncomps + c);
        const auto *kf_8 = buffer.data(kf + 8 * ncomps + c);
        const auto *kf_9 = buffer.data(kf + 9 * ncomps + c);
        const auto *kf_10 = buffer.data(kf + 10 * ncomps + c);
        const auto *kf_11 = buffer.data(kf + 11 * ncomps + c);
        const auto *kf_12 = buffer.data(kf + 12 * ncomps + c);
        const auto *kf_13 = buffer.data(kf + 13 * ncomps + c);
        const auto *kf_14 = buffer.data(kf + 14 * ncomps + c);
        const auto *kf_15 = buffer.data(kf + 15 * ncomps + c);
        const auto *kf_16 = buffer.data(kf + 16 * ncomps + c);
        const auto *kf_17 = buffer.data(kf + 17 * ncomps + c);
        const auto *kf_18 = buffer.data(kf + 18 * ncomps + c);
        const auto *kf_19 = buffer.data(kf + 19 * ncomps + c);
        const auto *kf_20 = buffer.data(kf + 20 * ncomps + c);
        const auto *kf_21 = buffer.data(kf + 21 * ncomps + c);
        const auto *kf_22 = buffer.data(kf + 22 * ncomps + c);
        const auto *kf_23 = buffer.data(kf + 23 * ncomps + c);
        const auto *kf_24 = buffer.data(kf + 24 * ncomps + c);
        const auto *kf_25 = buffer.data(kf + 25 * ncomps + c);
        const auto *kf_26 = buffer.data(kf + 26 * ncomps + c);
        const auto *kf_27 = buffer.data(kf + 27 * ncomps + c);
        const auto *kf_28 = buffer.data(kf + 28 * ncomps + c);
        const auto *kf_29 = buffer.data(kf + 29 * ncomps + c);
        const auto *kf_30 = buffer.data(kf + 30 * ncomps + c);
        const auto *kf_31 = buffer.data(kf + 31 * ncomps + c);
        const auto *kf_32 = buffer.data(kf + 32 * ncomps + c);
        const auto *kf_33 = buffer.data(kf + 33 * ncomps + c);
        const auto *kf_34 = buffer.data(kf + 34 * ncomps + c);
        const auto *kf_35 = buffer.data(kf + 35 * ncomps + c);
        const auto *kf_36 = buffer.data(kf + 36 * ncomps + c);
        const auto *kf_37 = buffer.data(kf + 37 * ncomps + c);
        const auto *kf_38 = buffer.data(kf + 38 * ncomps + c);
        const auto *kf_39 = buffer.data(kf + 39 * ncomps + c);
        const auto *kf_40 = buffer.data(kf + 40 * ncomps + c);
        const auto *kf_41 = buffer.data(kf + 41 * ncomps + c);
        const auto *kf_42 = buffer.data(kf + 42 * ncomps + c);
        const auto *kf_43 = buffer.data(kf + 43 * ncomps + c);
        const auto *kf_44 = buffer.data(kf + 44 * ncomps + c);
        const auto *kf_45 = buffer.data(kf + 45 * ncomps + c);
        const auto *kf_46 = buffer.data(kf + 46 * ncomps + c);
        const auto *kf_47 = buffer.data(kf + 47 * ncomps + c);
        const auto *kf_48 = buffer.data(kf + 48 * ncomps + c);
        const auto *kf_49 = buffer.data(kf + 49 * ncomps + c);
        const auto *kf_50 = buffer.data(kf + 50 * ncomps + c);
        const auto *kf_51 = buffer.data(kf + 51 * ncomps + c);
        const auto *kf_52 = buffer.data(kf + 52 * ncomps + c);
        const auto *kf_53 = buffer.data(kf + 53 * ncomps + c);
        const auto *kf_54 = buffer.data(kf + 54 * ncomps + c);
        const auto *kf_55 = buffer.data(kf + 55 * ncomps + c);
        const auto *kf_56 = buffer.data(kf + 56 * ncomps + c);
        const auto *kf_57 = buffer.data(kf + 57 * ncomps + c);
        const auto *kf_58 = buffer.data(kf + 58 * ncomps + c);
        const auto *kf_59 = buffer.data(kf + 59 * ncomps + c);
        const auto *kf_60 = buffer.data(kf + 60 * ncomps + c);
        const auto *kf_61 = buffer.data(kf + 61 * ncomps + c);
        const auto *kf_62 = buffer.data(kf + 62 * ncomps + c);
        const auto *kf_63 = buffer.data(kf + 63 * ncomps + c);
        const auto *kf_64 = buffer.data(kf + 64 * ncomps + c);
        const auto *kf_65 = buffer.data(kf + 65 * ncomps + c);
        const auto *kf_66 = buffer.data(kf + 66 * ncomps + c);
        const auto *kf_67 = buffer.data(kf + 67 * ncomps + c);
        const auto *kf_68 = buffer.data(kf + 68 * ncomps + c);
        const auto *kf_69 = buffer.data(kf + 69 * ncomps + c);
        const auto *kf_70 = buffer.data(kf + 70 * ncomps + c);
        const auto *kf_71 = buffer.data(kf + 71 * ncomps + c);
        const auto *kf_72 = buffer.data(kf + 72 * ncomps + c);
        const auto *kf_73 = buffer.data(kf + 73 * ncomps + c);
        const auto *kf_74 = buffer.data(kf + 74 * ncomps + c);
        const auto *kf_75 = buffer.data(kf + 75 * ncomps + c);
        const auto *kf_76 = buffer.data(kf + 76 * ncomps + c);
        const auto *kf_77 = buffer.data(kf + 77 * ncomps + c);
        const auto *kf_78 = buffer.data(kf + 78 * ncomps + c);
        const auto *kf_79 = buffer.data(kf + 79 * ncomps + c);
        const auto *kf_80 = buffer.data(kf + 80 * ncomps + c);
        const auto *kf_81 = buffer.data(kf + 81 * ncomps + c);
        const auto *kf_82 = buffer.data(kf + 82 * ncomps + c);
        const auto *kf_83 = buffer.data(kf + 83 * ncomps + c);
        const auto *kf_84 = buffer.data(kf + 84 * ncomps + c);
        const auto *kf_85 = buffer.data(kf + 85 * ncomps + c);
        const auto *kf_86 = buffer.data(kf + 86 * ncomps + c);
        const auto *kf_87 = buffer.data(kf + 87 * ncomps + c);
        const auto *kf_88 = buffer.data(kf + 88 * ncomps + c);
        const auto *kf_89 = buffer.data(kf + 89 * ncomps + c);
        const auto *kf_90 = buffer.data(kf + 90 * ncomps + c);
        const auto *kf_91 = buffer.data(kf + 91 * ncomps + c);
        const auto *kf_92 = buffer.data(kf + 92 * ncomps + c);
        const auto *kf_93 = buffer.data(kf + 93 * ncomps + c);
        const auto *kf_94 = buffer.data(kf + 94 * ncomps + c);
        const auto *kf_95 = buffer.data(kf + 95 * ncomps + c);
        const auto *kf_96 = buffer.data(kf + 96 * ncomps + c);
        const auto *kf_97 = buffer.data(kf + 97 * ncomps + c);
        const auto *kf_98 = buffer.data(kf + 98 * ncomps + c);
        const auto *kf_99 = buffer.data(kf + 99 * ncomps + c);

        const auto *lf_0 = buffer.data(lf + 0 * ncomps + c);
        const auto *lf_1 = buffer.data(lf + 1 * ncomps + c);
        const auto *lf_2 = buffer.data(lf + 2 * ncomps + c);
        const auto *lf_3 = buffer.data(lf + 3 * ncomps + c);
        const auto *lf_4 = buffer.data(lf + 4 * ncomps + c);
        const auto *lf_5 = buffer.data(lf + 5 * ncomps + c);
        const auto *lf_6 = buffer.data(lf + 6 * ncomps + c);
        const auto *lf_7 = buffer.data(lf + 7 * ncomps + c);
        const auto *lf_8 = buffer.data(lf + 8 * ncomps + c);
        const auto *lf_9 = buffer.data(lf + 9 * ncomps + c);
        const auto *lf_10 = buffer.data(lf + 10 * ncomps + c);
        const auto *lf_11 = buffer.data(lf + 11 * ncomps + c);
        const auto *lf_12 = buffer.data(lf + 12 * ncomps + c);
        const auto *lf_13 = buffer.data(lf + 13 * ncomps + c);
        const auto *lf_14 = buffer.data(lf + 14 * ncomps + c);
        const auto *lf_15 = buffer.data(lf + 15 * ncomps + c);
        const auto *lf_16 = buffer.data(lf + 16 * ncomps + c);
        const auto *lf_17 = buffer.data(lf + 17 * ncomps + c);
        const auto *lf_18 = buffer.data(lf + 18 * ncomps + c);
        const auto *lf_19 = buffer.data(lf + 19 * ncomps + c);
        const auto *lf_20 = buffer.data(lf + 20 * ncomps + c);
        const auto *lf_21 = buffer.data(lf + 21 * ncomps + c);
        const auto *lf_22 = buffer.data(lf + 22 * ncomps + c);
        const auto *lf_23 = buffer.data(lf + 23 * ncomps + c);
        const auto *lf_24 = buffer.data(lf + 24 * ncomps + c);
        const auto *lf_25 = buffer.data(lf + 25 * ncomps + c);
        const auto *lf_26 = buffer.data(lf + 26 * ncomps + c);
        const auto *lf_27 = buffer.data(lf + 27 * ncomps + c);
        const auto *lf_28 = buffer.data(lf + 28 * ncomps + c);
        const auto *lf_29 = buffer.data(lf + 29 * ncomps + c);
        const auto *lf_30 = buffer.data(lf + 30 * ncomps + c);
        const auto *lf_31 = buffer.data(lf + 31 * ncomps + c);
        const auto *lf_32 = buffer.data(lf + 32 * ncomps + c);
        const auto *lf_33 = buffer.data(lf + 33 * ncomps + c);
        const auto *lf_34 = buffer.data(lf + 34 * ncomps + c);
        const auto *lf_35 = buffer.data(lf + 35 * ncomps + c);
        const auto *lf_36 = buffer.data(lf + 36 * ncomps + c);
        const auto *lf_37 = buffer.data(lf + 37 * ncomps + c);
        const auto *lf_38 = buffer.data(lf + 38 * ncomps + c);
        const auto *lf_39 = buffer.data(lf + 39 * ncomps + c);
        const auto *lf_40 = buffer.data(lf + 40 * ncomps + c);
        const auto *lf_41 = buffer.data(lf + 41 * ncomps + c);
        const auto *lf_42 = buffer.data(lf + 42 * ncomps + c);
        const auto *lf_43 = buffer.data(lf + 43 * ncomps + c);
        const auto *lf_44 = buffer.data(lf + 44 * ncomps + c);
        const auto *lf_45 = buffer.data(lf + 45 * ncomps + c);
        const auto *lf_46 = buffer.data(lf + 46 * ncomps + c);
        const auto *lf_47 = buffer.data(lf + 47 * ncomps + c);
        const auto *lf_48 = buffer.data(lf + 48 * ncomps + c);
        const auto *lf_49 = buffer.data(lf + 49 * ncomps + c);
        const auto *lf_50 = buffer.data(lf + 50 * ncomps + c);
        const auto *lf_51 = buffer.data(lf + 51 * ncomps + c);
        const auto *lf_52 = buffer.data(lf + 52 * ncomps + c);
        const auto *lf_53 = buffer.data(lf + 53 * ncomps + c);
        const auto *lf_54 = buffer.data(lf + 54 * ncomps + c);
        const auto *lf_55 = buffer.data(lf + 55 * ncomps + c);
        const auto *lf_56 = buffer.data(lf + 56 * ncomps + c);
        const auto *lf_57 = buffer.data(lf + 57 * ncomps + c);
        const auto *lf_58 = buffer.data(lf + 58 * ncomps + c);
        const auto *lf_59 = buffer.data(lf + 59 * ncomps + c);
        const auto *lf_60 = buffer.data(lf + 60 * ncomps + c);
        const auto *lf_61 = buffer.data(lf + 61 * ncomps + c);
        const auto *lf_62 = buffer.data(lf + 62 * ncomps + c);
        const auto *lf_63 = buffer.data(lf + 63 * ncomps + c);
        const auto *lf_64 = buffer.data(lf + 64 * ncomps + c);
        const auto *lf_65 = buffer.data(lf + 65 * ncomps + c);
        const auto *lf_66 = buffer.data(lf + 66 * ncomps + c);
        const auto *lf_67 = buffer.data(lf + 67 * ncomps + c);
        const auto *lf_68 = buffer.data(lf + 68 * ncomps + c);
        const auto *lf_69 = buffer.data(lf + 69 * ncomps + c);
        const auto *lf_70 = buffer.data(lf + 70 * ncomps + c);
        const auto *lf_71 = buffer.data(lf + 71 * ncomps + c);
        const auto *lf_72 = buffer.data(lf + 72 * ncomps + c);
        const auto *lf_73 = buffer.data(lf + 73 * ncomps + c);
        const auto *lf_74 = buffer.data(lf + 74 * ncomps + c);
        const auto *lf_75 = buffer.data(lf + 75 * ncomps + c);
        const auto *lf_76 = buffer.data(lf + 76 * ncomps + c);
        const auto *lf_77 = buffer.data(lf + 77 * ncomps + c);
        const auto *lf_78 = buffer.data(lf + 78 * ncomps + c);
        const auto *lf_79 = buffer.data(lf + 79 * ncomps + c);
        const auto *lf_80 = buffer.data(lf + 80 * ncomps + c);
        const auto *lf_81 = buffer.data(lf + 81 * ncomps + c);
        const auto *lf_82 = buffer.data(lf + 82 * ncomps + c);
        const auto *lf_83 = buffer.data(lf + 83 * ncomps + c);
        const auto *lf_84 = buffer.data(lf + 84 * ncomps + c);
        const auto *lf_85 = buffer.data(lf + 85 * ncomps + c);
        const auto *lf_86 = buffer.data(lf + 86 * ncomps + c);
        const auto *lf_87 = buffer.data(lf + 87 * ncomps + c);
        const auto *lf_88 = buffer.data(lf + 88 * ncomps + c);
        const auto *lf_89 = buffer.data(lf + 89 * ncomps + c);
        const auto *lf_90 = buffer.data(lf + 90 * ncomps + c);
        const auto *lf_91 = buffer.data(lf + 91 * ncomps + c);
        const auto *lf_92 = buffer.data(lf + 92 * ncomps + c);
        const auto *lf_93 = buffer.data(lf + 93 * ncomps + c);
        const auto *lf_94 = buffer.data(lf + 94 * ncomps + c);
        const auto *lf_95 = buffer.data(lf + 95 * ncomps + c);
        const auto *lf_96 = buffer.data(lf + 96 * ncomps + c);
        const auto *lf_97 = buffer.data(lf + 97 * ncomps + c);
        const auto *lf_98 = buffer.data(lf + 98 * ncomps + c);
        const auto *lf_99 = buffer.data(lf + 99 * ncomps + c);
        const auto *lf_106 = buffer.data(lf + 106 * ncomps + c);
        const auto *lf_107 = buffer.data(lf + 107 * ncomps + c);
        const auto *lf_108 = buffer.data(lf + 108 * ncomps + c);
        const auto *lf_109 = buffer.data(lf + 109 * ncomps + c);
        const auto *lf_116 = buffer.data(lf + 116 * ncomps + c);
        const auto *lf_117 = buffer.data(lf + 117 * ncomps + c);
        const auto *lf_118 = buffer.data(lf + 118 * ncomps + c);
        const auto *lf_119 = buffer.data(lf + 119 * ncomps + c);
        const auto *lf_126 = buffer.data(lf + 126 * ncomps + c);
        const auto *lf_127 = buffer.data(lf + 127 * ncomps + c);
        const auto *lf_128 = buffer.data(lf + 128 * ncomps + c);
        const auto *lf_129 = buffer.data(lf + 129 * ncomps + c);
        const auto *lf_139 = buffer.data(lf + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, kf_0, kf_1, kf_2, kf_3, kf_4, lf_0, \
                         lf_1, lf_2, lf_3, lf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * kf_0[k]
                     + lf_0[k];

            t_1[k] = ab_x[k] * kf_1[k]
                     + lf_1[k];

            t_2[k] = ab_x[k] * kf_2[k]
                     + lf_2[k];

            t_3[k] = ab_x[k] * kf_3[k]
                     + lf_3[k];

            t_4[k] = ab_x[k] * kf_4[k]
                     + lf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, kf_5, kf_6, kf_7, kf_8, kf_9, lf_5, \
                         lf_6, lf_7, lf_8, lf_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * kf_5[k]
                     + lf_5[k];

            t_6[k] = ab_x[k] * kf_6[k]
                     + lf_6[k];

            t_7[k] = ab_x[k] * kf_7[k]
                     + lf_7[k];

            t_8[k] = ab_x[k] * kf_8[k]
                     + lf_8[k];

            t_9[k] = ab_x[k] * kf_9[k]
                     + lf_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, kf_6, kf_7, kf_8, kf_9, \
                         lf_16, lf_17, lf_18, lf_19, lf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * kf_6[k]
                      + lf_16[k];

            t_11[k] = ab_y[k] * kf_7[k]
                      + lf_17[k];

            t_12[k] = ab_y[k] * kf_8[k]
                      + lf_18[k];

            t_13[k] = ab_y[k] * kf_9[k]
                      + lf_19[k];

            t_14[k] = ab_z[k] * kf_9[k]
                      + lf_29[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, kf_10, kf_11, kf_12, kf_13, \
                         kf_14, lf_10, lf_11, lf_12, lf_13, lf_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * kf_10[k]
                      + lf_10[k];

            t_16[k] = ab_x[k] * kf_11[k]
                      + lf_11[k];

            t_17[k] = ab_x[k] * kf_12[k]
                      + lf_12[k];

            t_18[k] = ab_x[k] * kf_13[k]
                      + lf_13[k];

            t_19[k] = ab_x[k] * kf_14[k]
                      + lf_14[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, kf_15, kf_16, kf_17, kf_18, \
                         kf_19, lf_15, lf_16, lf_17, lf_18, lf_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * kf_15[k]
                      + lf_15[k];

            t_21[k] = ab_x[k] * kf_16[k]
                      + lf_16[k];

            t_22[k] = ab_x[k] * kf_17[k]
                      + lf_17[k];

            t_23[k] = ab_x[k] * kf_18[k]
                      + lf_18[k];

            t_24[k] = ab_x[k] * kf_19[k]
                      + lf_19[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, kf_16, kf_17, kf_18, kf_19, \
                         lf_36, lf_37, lf_38, lf_39, lf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_y[k] * kf_16[k]
                      + lf_36[k];

            t_26[k] = ab_y[k] * kf_17[k]
                      + lf_37[k];

            t_27[k] = ab_y[k] * kf_18[k]
                      + lf_38[k];

            t_28[k] = ab_y[k] * kf_19[k]
                      + lf_39[k];

            t_29[k] = ab_z[k] * kf_19[k]
                      + lf_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, kf_20, kf_21, kf_22, kf_23, \
                         kf_24, lf_20, lf_21, lf_22, lf_23, lf_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * kf_20[k]
                      + lf_20[k];

            t_31[k] = ab_x[k] * kf_21[k]
                      + lf_21[k];

            t_32[k] = ab_x[k] * kf_22[k]
                      + lf_22[k];

            t_33[k] = ab_x[k] * kf_23[k]
                      + lf_23[k];

            t_34[k] = ab_x[k] * kf_24[k]
                      + lf_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, kf_25, kf_26, kf_27, kf_28, \
                         kf_29, lf_25, lf_26, lf_27, lf_28, lf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * kf_25[k]
                      + lf_25[k];

            t_36[k] = ab_x[k] * kf_26[k]
                      + lf_26[k];

            t_37[k] = ab_x[k] * kf_27[k]
                      + lf_27[k];

            t_38[k] = ab_x[k] * kf_28[k]
                      + lf_28[k];

            t_39[k] = ab_x[k] * kf_29[k]
                      + lf_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, kf_26, kf_27, kf_28, kf_29, \
                         lf_46, lf_47, lf_48, lf_49, lf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * kf_26[k]
                      + lf_46[k];

            t_41[k] = ab_y[k] * kf_27[k]
                      + lf_47[k];

            t_42[k] = ab_y[k] * kf_28[k]
                      + lf_48[k];

            t_43[k] = ab_y[k] * kf_29[k]
                      + lf_49[k];

            t_44[k] = ab_z[k] * kf_29[k]
                      + lf_59[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, kf_30, kf_31, kf_32, kf_33, \
                         kf_34, lf_30, lf_31, lf_32, lf_33, lf_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * kf_30[k]
                      + lf_30[k];

            t_46[k] = ab_x[k] * kf_31[k]
                      + lf_31[k];

            t_47[k] = ab_x[k] * kf_32[k]
                      + lf_32[k];

            t_48[k] = ab_x[k] * kf_33[k]
                      + lf_33[k];

            t_49[k] = ab_x[k] * kf_34[k]
                      + lf_34[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, kf_35, kf_36, kf_37, kf_38, \
                         kf_39, lf_35, lf_36, lf_37, lf_38, lf_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * kf_35[k]
                      + lf_35[k];

            t_51[k] = ab_x[k] * kf_36[k]
                      + lf_36[k];

            t_52[k] = ab_x[k] * kf_37[k]
                      + lf_37[k];

            t_53[k] = ab_x[k] * kf_38[k]
                      + lf_38[k];

            t_54[k] = ab_x[k] * kf_39[k]
                      + lf_39[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, kf_36, kf_37, kf_38, kf_39, \
                         lf_66, lf_67, lf_68, lf_69, lf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_y[k] * kf_36[k]
                      + lf_66[k];

            t_56[k] = ab_y[k] * kf_37[k]
                      + lf_67[k];

            t_57[k] = ab_y[k] * kf_38[k]
                      + lf_68[k];

            t_58[k] = ab_y[k] * kf_39[k]
                      + lf_69[k];

            t_59[k] = ab_z[k] * kf_39[k]
                      + lf_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, kf_40, kf_41, kf_42, kf_43, \
                         kf_44, lf_40, lf_41, lf_42, lf_43, lf_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * kf_40[k]
                      + lf_40[k];

            t_61[k] = ab_x[k] * kf_41[k]
                      + lf_41[k];

            t_62[k] = ab_x[k] * kf_42[k]
                      + lf_42[k];

            t_63[k] = ab_x[k] * kf_43[k]
                      + lf_43[k];

            t_64[k] = ab_x[k] * kf_44[k]
                      + lf_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, kf_45, kf_46, kf_47, kf_48, \
                         kf_49, lf_45, lf_46, lf_47, lf_48, lf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * kf_45[k]
                      + lf_45[k];

            t_66[k] = ab_x[k] * kf_46[k]
                      + lf_46[k];

            t_67[k] = ab_x[k] * kf_47[k]
                      + lf_47[k];

            t_68[k] = ab_x[k] * kf_48[k]
                      + lf_48[k];

            t_69[k] = ab_x[k] * kf_49[k]
                      + lf_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, kf_46, kf_47, kf_48, kf_49, \
                         lf_76, lf_77, lf_78, lf_79, lf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * kf_46[k]
                      + lf_76[k];

            t_71[k] = ab_y[k] * kf_47[k]
                      + lf_77[k];

            t_72[k] = ab_y[k] * kf_48[k]
                      + lf_78[k];

            t_73[k] = ab_y[k] * kf_49[k]
                      + lf_79[k];

            t_74[k] = ab_z[k] * kf_49[k]
                      + lf_89[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, kf_50, kf_51, kf_52, kf_53, \
                         kf_54, lf_50, lf_51, lf_52, lf_53, lf_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * kf_50[k]
                      + lf_50[k];

            t_76[k] = ab_x[k] * kf_51[k]
                      + lf_51[k];

            t_77[k] = ab_x[k] * kf_52[k]
                      + lf_52[k];

            t_78[k] = ab_x[k] * kf_53[k]
                      + lf_53[k];

            t_79[k] = ab_x[k] * kf_54[k]
                      + lf_54[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, kf_55, kf_56, kf_57, kf_58, \
                         kf_59, lf_55, lf_56, lf_57, lf_58, lf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * kf_55[k]
                      + lf_55[k];

            t_81[k] = ab_x[k] * kf_56[k]
                      + lf_56[k];

            t_82[k] = ab_x[k] * kf_57[k]
                      + lf_57[k];

            t_83[k] = ab_x[k] * kf_58[k]
                      + lf_58[k];

            t_84[k] = ab_x[k] * kf_59[k]
                      + lf_59[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, kf_56, kf_57, kf_58, kf_59, \
                         lf_86, lf_87, lf_88, lf_89, lf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_y[k] * kf_56[k]
                      + lf_86[k];

            t_86[k] = ab_y[k] * kf_57[k]
                      + lf_87[k];

            t_87[k] = ab_y[k] * kf_58[k]
                      + lf_88[k];

            t_88[k] = ab_y[k] * kf_59[k]
                      + lf_89[k];

            t_89[k] = ab_z[k] * kf_59[k]
                      + lf_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, kf_60, kf_61, kf_62, kf_63, \
                         kf_64, lf_60, lf_61, lf_62, lf_63, lf_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * kf_60[k]
                      + lf_60[k];

            t_91[k] = ab_x[k] * kf_61[k]
                      + lf_61[k];

            t_92[k] = ab_x[k] * kf_62[k]
                      + lf_62[k];

            t_93[k] = ab_x[k] * kf_63[k]
                      + lf_63[k];

            t_94[k] = ab_x[k] * kf_64[k]
                      + lf_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, kf_65, kf_66, kf_67, kf_68, \
                         kf_69, lf_65, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * kf_65[k]
                      + lf_65[k];

            t_96[k] = ab_x[k] * kf_66[k]
                      + lf_66[k];

            t_97[k] = ab_x[k] * kf_67[k]
                      + lf_67[k];

            t_98[k] = ab_x[k] * kf_68[k]
                      + lf_68[k];

            t_99[k] = ab_x[k] * kf_69[k]
                      + lf_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, kf_66, kf_67, kf_68, \
                         kf_69, lf_106, lf_107, lf_108, lf_109, \
                         lf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * kf_66[k]
                       + lf_106[k];

            t_101[k] = ab_y[k] * kf_67[k]
                       + lf_107[k];

            t_102[k] = ab_y[k] * kf_68[k]
                       + lf_108[k];

            t_103[k] = ab_y[k] * kf_69[k]
                       + lf_109[k];

            t_104[k] = ab_z[k] * kf_69[k]
                       + lf_119[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, kf_70, kf_71, kf_72, kf_73, \
                         kf_74, lf_70, lf_71, lf_72, lf_73, lf_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * kf_70[k]
                       + lf_70[k];

            t_106[k] = ab_x[k] * kf_71[k]
                       + lf_71[k];

            t_107[k] = ab_x[k] * kf_72[k]
                       + lf_72[k];

            t_108[k] = ab_x[k] * kf_73[k]
                       + lf_73[k];

            t_109[k] = ab_x[k] * kf_74[k]
                       + lf_74[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, kf_75, kf_76, kf_77, kf_78, \
                         kf_79, lf_75, lf_76, lf_77, lf_78, lf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * kf_75[k]
                       + lf_75[k];

            t_111[k] = ab_x[k] * kf_76[k]
                       + lf_76[k];

            t_112[k] = ab_x[k] * kf_77[k]
                       + lf_77[k];

            t_113[k] = ab_x[k] * kf_78[k]
                       + lf_78[k];

            t_114[k] = ab_x[k] * kf_79[k]
                       + lf_79[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, kf_76, kf_77, kf_78, \
                         kf_79, lf_116, lf_117, lf_118, lf_119, \
                         lf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_y[k] * kf_76[k]
                       + lf_116[k];

            t_116[k] = ab_y[k] * kf_77[k]
                       + lf_117[k];

            t_117[k] = ab_y[k] * kf_78[k]
                       + lf_118[k];

            t_118[k] = ab_y[k] * kf_79[k]
                       + lf_119[k];

            t_119[k] = ab_z[k] * kf_79[k]
                       + lf_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, kf_80, kf_81, kf_82, kf_83, \
                         kf_84, lf_80, lf_81, lf_82, lf_83, lf_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * kf_80[k]
                       + lf_80[k];

            t_121[k] = ab_x[k] * kf_81[k]
                       + lf_81[k];

            t_122[k] = ab_x[k] * kf_82[k]
                       + lf_82[k];

            t_123[k] = ab_x[k] * kf_83[k]
                       + lf_83[k];

            t_124[k] = ab_x[k] * kf_84[k]
                       + lf_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, kf_85, kf_86, kf_87, kf_88, \
                         kf_89, lf_85, lf_86, lf_87, lf_88, lf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * kf_85[k]
                       + lf_85[k];

            t_126[k] = ab_x[k] * kf_86[k]
                       + lf_86[k];

            t_127[k] = ab_x[k] * kf_87[k]
                       + lf_87[k];

            t_128[k] = ab_x[k] * kf_88[k]
                       + lf_88[k];

            t_129[k] = ab_x[k] * kf_89[k]
                       + lf_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, kf_86, kf_87, kf_88, \
                         kf_89, lf_126, lf_127, lf_128, lf_129, \
                         lf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * kf_86[k]
                       + lf_126[k];

            t_131[k] = ab_y[k] * kf_87[k]
                       + lf_127[k];

            t_132[k] = ab_y[k] * kf_88[k]
                       + lf_128[k];

            t_133[k] = ab_y[k] * kf_89[k]
                       + lf_129[k];

            t_134[k] = ab_z[k] * kf_89[k]
                       + lf_139[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, kf_90, kf_91, kf_92, kf_93, \
                         kf_94, lf_90, lf_91, lf_92, lf_93, lf_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * kf_90[k]
                       + lf_90[k];

            t_136[k] = ab_x[k] * kf_91[k]
                       + lf_91[k];

            t_137[k] = ab_x[k] * kf_92[k]
                       + lf_92[k];

            t_138[k] = ab_x[k] * kf_93[k]
                       + lf_93[k];

            t_139[k] = ab_x[k] * kf_94[k]
                       + lf_94[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, kf_95, kf_96, kf_97, kf_98, \
                         kf_99, lf_95, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * kf_95[k]
                       + lf_95[k];

            t_141[k] = ab_x[k] * kf_96[k]
                       + lf_96[k];

            t_142[k] = ab_x[k] * kf_97[k]
                       + lf_97[k];

            t_143[k] = ab_x[k] * kf_98[k]
                       + lf_98[k];

            t_144[k] = ab_x[k] * kf_99[k]
                       + lf_99[k];
        }
    }
}

static auto
compute_hrr_kg_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t kf, const size_t lf, const size_t ncomps,
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

        const auto *kf_96 = buffer.data(kf + 96 * ncomps + c);
        const auto *kf_97 = buffer.data(kf + 97 * ncomps + c);
        const auto *kf_98 = buffer.data(kf + 98 * ncomps + c);
        const auto *kf_99 = buffer.data(kf + 99 * ncomps + c);
        const auto *kf_100 = buffer.data(kf + 100 * ncomps + c);
        const auto *kf_101 = buffer.data(kf + 101 * ncomps + c);
        const auto *kf_102 = buffer.data(kf + 102 * ncomps + c);
        const auto *kf_103 = buffer.data(kf + 103 * ncomps + c);
        const auto *kf_104 = buffer.data(kf + 104 * ncomps + c);
        const auto *kf_105 = buffer.data(kf + 105 * ncomps + c);
        const auto *kf_106 = buffer.data(kf + 106 * ncomps + c);
        const auto *kf_107 = buffer.data(kf + 107 * ncomps + c);
        const auto *kf_108 = buffer.data(kf + 108 * ncomps + c);
        const auto *kf_109 = buffer.data(kf + 109 * ncomps + c);
        const auto *kf_110 = buffer.data(kf + 110 * ncomps + c);
        const auto *kf_111 = buffer.data(kf + 111 * ncomps + c);
        const auto *kf_112 = buffer.data(kf + 112 * ncomps + c);
        const auto *kf_113 = buffer.data(kf + 113 * ncomps + c);
        const auto *kf_114 = buffer.data(kf + 114 * ncomps + c);
        const auto *kf_115 = buffer.data(kf + 115 * ncomps + c);
        const auto *kf_116 = buffer.data(kf + 116 * ncomps + c);
        const auto *kf_117 = buffer.data(kf + 117 * ncomps + c);
        const auto *kf_118 = buffer.data(kf + 118 * ncomps + c);
        const auto *kf_119 = buffer.data(kf + 119 * ncomps + c);
        const auto *kf_120 = buffer.data(kf + 120 * ncomps + c);
        const auto *kf_121 = buffer.data(kf + 121 * ncomps + c);
        const auto *kf_122 = buffer.data(kf + 122 * ncomps + c);
        const auto *kf_123 = buffer.data(kf + 123 * ncomps + c);
        const auto *kf_124 = buffer.data(kf + 124 * ncomps + c);
        const auto *kf_125 = buffer.data(kf + 125 * ncomps + c);
        const auto *kf_126 = buffer.data(kf + 126 * ncomps + c);
        const auto *kf_127 = buffer.data(kf + 127 * ncomps + c);
        const auto *kf_128 = buffer.data(kf + 128 * ncomps + c);
        const auto *kf_129 = buffer.data(kf + 129 * ncomps + c);
        const auto *kf_130 = buffer.data(kf + 130 * ncomps + c);
        const auto *kf_131 = buffer.data(kf + 131 * ncomps + c);
        const auto *kf_132 = buffer.data(kf + 132 * ncomps + c);
        const auto *kf_133 = buffer.data(kf + 133 * ncomps + c);
        const auto *kf_134 = buffer.data(kf + 134 * ncomps + c);
        const auto *kf_135 = buffer.data(kf + 135 * ncomps + c);
        const auto *kf_136 = buffer.data(kf + 136 * ncomps + c);
        const auto *kf_137 = buffer.data(kf + 137 * ncomps + c);
        const auto *kf_138 = buffer.data(kf + 138 * ncomps + c);
        const auto *kf_139 = buffer.data(kf + 139 * ncomps + c);
        const auto *kf_140 = buffer.data(kf + 140 * ncomps + c);
        const auto *kf_141 = buffer.data(kf + 141 * ncomps + c);
        const auto *kf_142 = buffer.data(kf + 142 * ncomps + c);
        const auto *kf_143 = buffer.data(kf + 143 * ncomps + c);
        const auto *kf_144 = buffer.data(kf + 144 * ncomps + c);
        const auto *kf_145 = buffer.data(kf + 145 * ncomps + c);
        const auto *kf_146 = buffer.data(kf + 146 * ncomps + c);
        const auto *kf_147 = buffer.data(kf + 147 * ncomps + c);
        const auto *kf_148 = buffer.data(kf + 148 * ncomps + c);
        const auto *kf_149 = buffer.data(kf + 149 * ncomps + c);
        const auto *kf_150 = buffer.data(kf + 150 * ncomps + c);
        const auto *kf_151 = buffer.data(kf + 151 * ncomps + c);
        const auto *kf_152 = buffer.data(kf + 152 * ncomps + c);
        const auto *kf_153 = buffer.data(kf + 153 * ncomps + c);
        const auto *kf_154 = buffer.data(kf + 154 * ncomps + c);
        const auto *kf_155 = buffer.data(kf + 155 * ncomps + c);
        const auto *kf_156 = buffer.data(kf + 156 * ncomps + c);
        const auto *kf_157 = buffer.data(kf + 157 * ncomps + c);
        const auto *kf_158 = buffer.data(kf + 158 * ncomps + c);
        const auto *kf_159 = buffer.data(kf + 159 * ncomps + c);
        const auto *kf_160 = buffer.data(kf + 160 * ncomps + c);
        const auto *kf_161 = buffer.data(kf + 161 * ncomps + c);
        const auto *kf_162 = buffer.data(kf + 162 * ncomps + c);
        const auto *kf_163 = buffer.data(kf + 163 * ncomps + c);
        const auto *kf_164 = buffer.data(kf + 164 * ncomps + c);
        const auto *kf_165 = buffer.data(kf + 165 * ncomps + c);
        const auto *kf_166 = buffer.data(kf + 166 * ncomps + c);
        const auto *kf_167 = buffer.data(kf + 167 * ncomps + c);
        const auto *kf_168 = buffer.data(kf + 168 * ncomps + c);
        const auto *kf_169 = buffer.data(kf + 169 * ncomps + c);
        const auto *kf_170 = buffer.data(kf + 170 * ncomps + c);
        const auto *kf_171 = buffer.data(kf + 171 * ncomps + c);
        const auto *kf_172 = buffer.data(kf + 172 * ncomps + c);
        const auto *kf_173 = buffer.data(kf + 173 * ncomps + c);
        const auto *kf_174 = buffer.data(kf + 174 * ncomps + c);
        const auto *kf_175 = buffer.data(kf + 175 * ncomps + c);
        const auto *kf_176 = buffer.data(kf + 176 * ncomps + c);
        const auto *kf_177 = buffer.data(kf + 177 * ncomps + c);
        const auto *kf_178 = buffer.data(kf + 178 * ncomps + c);
        const auto *kf_179 = buffer.data(kf + 179 * ncomps + c);
        const auto *kf_180 = buffer.data(kf + 180 * ncomps + c);
        const auto *kf_181 = buffer.data(kf + 181 * ncomps + c);
        const auto *kf_182 = buffer.data(kf + 182 * ncomps + c);
        const auto *kf_183 = buffer.data(kf + 183 * ncomps + c);
        const auto *kf_184 = buffer.data(kf + 184 * ncomps + c);
        const auto *kf_185 = buffer.data(kf + 185 * ncomps + c);
        const auto *kf_186 = buffer.data(kf + 186 * ncomps + c);
        const auto *kf_187 = buffer.data(kf + 187 * ncomps + c);
        const auto *kf_188 = buffer.data(kf + 188 * ncomps + c);
        const auto *kf_189 = buffer.data(kf + 189 * ncomps + c);
        const auto *kf_190 = buffer.data(kf + 190 * ncomps + c);
        const auto *kf_191 = buffer.data(kf + 191 * ncomps + c);
        const auto *kf_192 = buffer.data(kf + 192 * ncomps + c);
        const auto *kf_193 = buffer.data(kf + 193 * ncomps + c);
        const auto *kf_194 = buffer.data(kf + 194 * ncomps + c);

        const auto *lf_100 = buffer.data(lf + 100 * ncomps + c);
        const auto *lf_101 = buffer.data(lf + 101 * ncomps + c);
        const auto *lf_102 = buffer.data(lf + 102 * ncomps + c);
        const auto *lf_103 = buffer.data(lf + 103 * ncomps + c);
        const auto *lf_104 = buffer.data(lf + 104 * ncomps + c);
        const auto *lf_105 = buffer.data(lf + 105 * ncomps + c);
        const auto *lf_106 = buffer.data(lf + 106 * ncomps + c);
        const auto *lf_107 = buffer.data(lf + 107 * ncomps + c);
        const auto *lf_108 = buffer.data(lf + 108 * ncomps + c);
        const auto *lf_109 = buffer.data(lf + 109 * ncomps + c);
        const auto *lf_110 = buffer.data(lf + 110 * ncomps + c);
        const auto *lf_111 = buffer.data(lf + 111 * ncomps + c);
        const auto *lf_112 = buffer.data(lf + 112 * ncomps + c);
        const auto *lf_113 = buffer.data(lf + 113 * ncomps + c);
        const auto *lf_114 = buffer.data(lf + 114 * ncomps + c);
        const auto *lf_115 = buffer.data(lf + 115 * ncomps + c);
        const auto *lf_116 = buffer.data(lf + 116 * ncomps + c);
        const auto *lf_117 = buffer.data(lf + 117 * ncomps + c);
        const auto *lf_118 = buffer.data(lf + 118 * ncomps + c);
        const auto *lf_119 = buffer.data(lf + 119 * ncomps + c);
        const auto *lf_120 = buffer.data(lf + 120 * ncomps + c);
        const auto *lf_121 = buffer.data(lf + 121 * ncomps + c);
        const auto *lf_122 = buffer.data(lf + 122 * ncomps + c);
        const auto *lf_123 = buffer.data(lf + 123 * ncomps + c);
        const auto *lf_124 = buffer.data(lf + 124 * ncomps + c);
        const auto *lf_125 = buffer.data(lf + 125 * ncomps + c);
        const auto *lf_126 = buffer.data(lf + 126 * ncomps + c);
        const auto *lf_127 = buffer.data(lf + 127 * ncomps + c);
        const auto *lf_128 = buffer.data(lf + 128 * ncomps + c);
        const auto *lf_129 = buffer.data(lf + 129 * ncomps + c);
        const auto *lf_130 = buffer.data(lf + 130 * ncomps + c);
        const auto *lf_131 = buffer.data(lf + 131 * ncomps + c);
        const auto *lf_132 = buffer.data(lf + 132 * ncomps + c);
        const auto *lf_133 = buffer.data(lf + 133 * ncomps + c);
        const auto *lf_134 = buffer.data(lf + 134 * ncomps + c);
        const auto *lf_135 = buffer.data(lf + 135 * ncomps + c);
        const auto *lf_136 = buffer.data(lf + 136 * ncomps + c);
        const auto *lf_137 = buffer.data(lf + 137 * ncomps + c);
        const auto *lf_138 = buffer.data(lf + 138 * ncomps + c);
        const auto *lf_139 = buffer.data(lf + 139 * ncomps + c);
        const auto *lf_140 = buffer.data(lf + 140 * ncomps + c);
        const auto *lf_141 = buffer.data(lf + 141 * ncomps + c);
        const auto *lf_142 = buffer.data(lf + 142 * ncomps + c);
        const auto *lf_143 = buffer.data(lf + 143 * ncomps + c);
        const auto *lf_144 = buffer.data(lf + 144 * ncomps + c);
        const auto *lf_145 = buffer.data(lf + 145 * ncomps + c);
        const auto *lf_146 = buffer.data(lf + 146 * ncomps + c);
        const auto *lf_147 = buffer.data(lf + 147 * ncomps + c);
        const auto *lf_148 = buffer.data(lf + 148 * ncomps + c);
        const auto *lf_149 = buffer.data(lf + 149 * ncomps + c);
        const auto *lf_150 = buffer.data(lf + 150 * ncomps + c);
        const auto *lf_151 = buffer.data(lf + 151 * ncomps + c);
        const auto *lf_152 = buffer.data(lf + 152 * ncomps + c);
        const auto *lf_153 = buffer.data(lf + 153 * ncomps + c);
        const auto *lf_154 = buffer.data(lf + 154 * ncomps + c);
        const auto *lf_155 = buffer.data(lf + 155 * ncomps + c);
        const auto *lf_156 = buffer.data(lf + 156 * ncomps + c);
        const auto *lf_157 = buffer.data(lf + 157 * ncomps + c);
        const auto *lf_158 = buffer.data(lf + 158 * ncomps + c);
        const auto *lf_159 = buffer.data(lf + 159 * ncomps + c);
        const auto *lf_160 = buffer.data(lf + 160 * ncomps + c);
        const auto *lf_161 = buffer.data(lf + 161 * ncomps + c);
        const auto *lf_162 = buffer.data(lf + 162 * ncomps + c);
        const auto *lf_163 = buffer.data(lf + 163 * ncomps + c);
        const auto *lf_164 = buffer.data(lf + 164 * ncomps + c);
        const auto *lf_165 = buffer.data(lf + 165 * ncomps + c);
        const auto *lf_166 = buffer.data(lf + 166 * ncomps + c);
        const auto *lf_167 = buffer.data(lf + 167 * ncomps + c);
        const auto *lf_168 = buffer.data(lf + 168 * ncomps + c);
        const auto *lf_169 = buffer.data(lf + 169 * ncomps + c);
        const auto *lf_170 = buffer.data(lf + 170 * ncomps + c);
        const auto *lf_171 = buffer.data(lf + 171 * ncomps + c);
        const auto *lf_172 = buffer.data(lf + 172 * ncomps + c);
        const auto *lf_173 = buffer.data(lf + 173 * ncomps + c);
        const auto *lf_174 = buffer.data(lf + 174 * ncomps + c);
        const auto *lf_175 = buffer.data(lf + 175 * ncomps + c);
        const auto *lf_176 = buffer.data(lf + 176 * ncomps + c);
        const auto *lf_177 = buffer.data(lf + 177 * ncomps + c);
        const auto *lf_178 = buffer.data(lf + 178 * ncomps + c);
        const auto *lf_179 = buffer.data(lf + 179 * ncomps + c);
        const auto *lf_180 = buffer.data(lf + 180 * ncomps + c);
        const auto *lf_181 = buffer.data(lf + 181 * ncomps + c);
        const auto *lf_182 = buffer.data(lf + 182 * ncomps + c);
        const auto *lf_183 = buffer.data(lf + 183 * ncomps + c);
        const auto *lf_184 = buffer.data(lf + 184 * ncomps + c);
        const auto *lf_185 = buffer.data(lf + 185 * ncomps + c);
        const auto *lf_186 = buffer.data(lf + 186 * ncomps + c);
        const auto *lf_187 = buffer.data(lf + 187 * ncomps + c);
        const auto *lf_188 = buffer.data(lf + 188 * ncomps + c);
        const auto *lf_189 = buffer.data(lf + 189 * ncomps + c);
        const auto *lf_190 = buffer.data(lf + 190 * ncomps + c);
        const auto *lf_191 = buffer.data(lf + 191 * ncomps + c);
        const auto *lf_192 = buffer.data(lf + 192 * ncomps + c);
        const auto *lf_193 = buffer.data(lf + 193 * ncomps + c);
        const auto *lf_194 = buffer.data(lf + 194 * ncomps + c);
        const auto *lf_196 = buffer.data(lf + 196 * ncomps + c);
        const auto *lf_197 = buffer.data(lf + 197 * ncomps + c);
        const auto *lf_198 = buffer.data(lf + 198 * ncomps + c);
        const auto *lf_199 = buffer.data(lf + 199 * ncomps + c);
        const auto *lf_209 = buffer.data(lf + 209 * ncomps + c);
        const auto *lf_216 = buffer.data(lf + 216 * ncomps + c);
        const auto *lf_217 = buffer.data(lf + 217 * ncomps + c);
        const auto *lf_218 = buffer.data(lf + 218 * ncomps + c);
        const auto *lf_219 = buffer.data(lf + 219 * ncomps + c);
        const auto *lf_226 = buffer.data(lf + 226 * ncomps + c);
        const auto *lf_227 = buffer.data(lf + 227 * ncomps + c);
        const auto *lf_228 = buffer.data(lf + 228 * ncomps + c);
        const auto *lf_229 = buffer.data(lf + 229 * ncomps + c);
        const auto *lf_236 = buffer.data(lf + 236 * ncomps + c);
        const auto *lf_237 = buffer.data(lf + 237 * ncomps + c);
        const auto *lf_238 = buffer.data(lf + 238 * ncomps + c);
        const auto *lf_239 = buffer.data(lf + 239 * ncomps + c);
        const auto *lf_246 = buffer.data(lf + 246 * ncomps + c);
        const auto *lf_247 = buffer.data(lf + 247 * ncomps + c);
        const auto *lf_248 = buffer.data(lf + 248 * ncomps + c);
        const auto *lf_249 = buffer.data(lf + 249 * ncomps + c);
        const auto *lf_259 = buffer.data(lf + 259 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, kf_96, kf_97, kf_98, \
                         kf_99, lf_136, lf_137, lf_138, lf_139, \
                         lf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_y[k] * kf_96[k]
                       + lf_136[k];

            t_146[k] = ab_y[k] * kf_97[k]
                       + lf_137[k];

            t_147[k] = ab_y[k] * kf_98[k]
                       + lf_138[k];

            t_148[k] = ab_y[k] * kf_99[k]
                       + lf_139[k];

            t_149[k] = ab_z[k] * kf_99[k]
                       + lf_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, kf_100, kf_101, kf_102, \
                         kf_103, kf_104, lf_100, lf_101, lf_102, lf_103, \
                         lf_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * kf_100[k]
                       + lf_100[k];

            t_151[k] = ab_x[k] * kf_101[k]
                       + lf_101[k];

            t_152[k] = ab_x[k] * kf_102[k]
                       + lf_102[k];

            t_153[k] = ab_x[k] * kf_103[k]
                       + lf_103[k];

            t_154[k] = ab_x[k] * kf_104[k]
                       + lf_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, kf_105, kf_106, kf_107, \
                         kf_108, kf_109, lf_105, lf_106, lf_107, lf_108, \
                         lf_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * kf_105[k]
                       + lf_105[k];

            t_156[k] = ab_x[k] * kf_106[k]
                       + lf_106[k];

            t_157[k] = ab_x[k] * kf_107[k]
                       + lf_107[k];

            t_158[k] = ab_x[k] * kf_108[k]
                       + lf_108[k];

            t_159[k] = ab_x[k] * kf_109[k]
                       + lf_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, kf_106, kf_107, \
                         kf_108, kf_109, lf_156, lf_157, lf_158, lf_159, \
                         lf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * kf_106[k]
                       + lf_156[k];

            t_161[k] = ab_y[k] * kf_107[k]
                       + lf_157[k];

            t_162[k] = ab_y[k] * kf_108[k]
                       + lf_158[k];

            t_163[k] = ab_y[k] * kf_109[k]
                       + lf_159[k];

            t_164[k] = ab_z[k] * kf_109[k]
                       + lf_169[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, kf_110, kf_111, kf_112, \
                         kf_113, kf_114, lf_110, lf_111, lf_112, lf_113, \
                         lf_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * kf_110[k]
                       + lf_110[k];

            t_166[k] = ab_x[k] * kf_111[k]
                       + lf_111[k];

            t_167[k] = ab_x[k] * kf_112[k]
                       + lf_112[k];

            t_168[k] = ab_x[k] * kf_113[k]
                       + lf_113[k];

            t_169[k] = ab_x[k] * kf_114[k]
                       + lf_114[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, kf_115, kf_116, kf_117, \
                         kf_118, kf_119, lf_115, lf_116, lf_117, lf_118, \
                         lf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * kf_115[k]
                       + lf_115[k];

            t_171[k] = ab_x[k] * kf_116[k]
                       + lf_116[k];

            t_172[k] = ab_x[k] * kf_117[k]
                       + lf_117[k];

            t_173[k] = ab_x[k] * kf_118[k]
                       + lf_118[k];

            t_174[k] = ab_x[k] * kf_119[k]
                       + lf_119[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, kf_116, kf_117, \
                         kf_118, kf_119, lf_166, lf_167, lf_168, lf_169, \
                         lf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_y[k] * kf_116[k]
                       + lf_166[k];

            t_176[k] = ab_y[k] * kf_117[k]
                       + lf_167[k];

            t_177[k] = ab_y[k] * kf_118[k]
                       + lf_168[k];

            t_178[k] = ab_y[k] * kf_119[k]
                       + lf_169[k];

            t_179[k] = ab_z[k] * kf_119[k]
                       + lf_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, kf_120, kf_121, kf_122, \
                         kf_123, kf_124, lf_120, lf_121, lf_122, lf_123, \
                         lf_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * kf_120[k]
                       + lf_120[k];

            t_181[k] = ab_x[k] * kf_121[k]
                       + lf_121[k];

            t_182[k] = ab_x[k] * kf_122[k]
                       + lf_122[k];

            t_183[k] = ab_x[k] * kf_123[k]
                       + lf_123[k];

            t_184[k] = ab_x[k] * kf_124[k]
                       + lf_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, kf_125, kf_126, kf_127, \
                         kf_128, kf_129, lf_125, lf_126, lf_127, lf_128, \
                         lf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * kf_125[k]
                       + lf_125[k];

            t_186[k] = ab_x[k] * kf_126[k]
                       + lf_126[k];

            t_187[k] = ab_x[k] * kf_127[k]
                       + lf_127[k];

            t_188[k] = ab_x[k] * kf_128[k]
                       + lf_128[k];

            t_189[k] = ab_x[k] * kf_129[k]
                       + lf_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, kf_126, kf_127, \
                         kf_128, kf_129, lf_176, lf_177, lf_178, lf_179, \
                         lf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * kf_126[k]
                       + lf_176[k];

            t_191[k] = ab_y[k] * kf_127[k]
                       + lf_177[k];

            t_192[k] = ab_y[k] * kf_128[k]
                       + lf_178[k];

            t_193[k] = ab_y[k] * kf_129[k]
                       + lf_179[k];

            t_194[k] = ab_z[k] * kf_129[k]
                       + lf_189[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, kf_130, kf_131, kf_132, \
                         kf_133, kf_134, lf_130, lf_131, lf_132, lf_133, \
                         lf_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * kf_130[k]
                       + lf_130[k];

            t_196[k] = ab_x[k] * kf_131[k]
                       + lf_131[k];

            t_197[k] = ab_x[k] * kf_132[k]
                       + lf_132[k];

            t_198[k] = ab_x[k] * kf_133[k]
                       + lf_133[k];

            t_199[k] = ab_x[k] * kf_134[k]
                       + lf_134[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, kf_135, kf_136, kf_137, \
                         kf_138, kf_139, lf_135, lf_136, lf_137, lf_138, \
                         lf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * kf_135[k]
                       + lf_135[k];

            t_201[k] = ab_x[k] * kf_136[k]
                       + lf_136[k];

            t_202[k] = ab_x[k] * kf_137[k]
                       + lf_137[k];

            t_203[k] = ab_x[k] * kf_138[k]
                       + lf_138[k];

            t_204[k] = ab_x[k] * kf_139[k]
                       + lf_139[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, kf_136, kf_137, \
                         kf_138, kf_139, lf_186, lf_187, lf_188, lf_189, \
                         lf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_y[k] * kf_136[k]
                       + lf_186[k];

            t_206[k] = ab_y[k] * kf_137[k]
                       + lf_187[k];

            t_207[k] = ab_y[k] * kf_138[k]
                       + lf_188[k];

            t_208[k] = ab_y[k] * kf_139[k]
                       + lf_189[k];

            t_209[k] = ab_z[k] * kf_139[k]
                       + lf_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, kf_140, kf_141, kf_142, \
                         kf_143, kf_144, lf_140, lf_141, lf_142, lf_143, \
                         lf_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * kf_140[k]
                       + lf_140[k];

            t_211[k] = ab_x[k] * kf_141[k]
                       + lf_141[k];

            t_212[k] = ab_x[k] * kf_142[k]
                       + lf_142[k];

            t_213[k] = ab_x[k] * kf_143[k]
                       + lf_143[k];

            t_214[k] = ab_x[k] * kf_144[k]
                       + lf_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, kf_145, kf_146, kf_147, \
                         kf_148, kf_149, lf_145, lf_146, lf_147, lf_148, \
                         lf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * kf_145[k]
                       + lf_145[k];

            t_216[k] = ab_x[k] * kf_146[k]
                       + lf_146[k];

            t_217[k] = ab_x[k] * kf_147[k]
                       + lf_147[k];

            t_218[k] = ab_x[k] * kf_148[k]
                       + lf_148[k];

            t_219[k] = ab_x[k] * kf_149[k]
                       + lf_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, kf_146, kf_147, \
                         kf_148, kf_149, lf_196, lf_197, lf_198, lf_199, \
                         lf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * kf_146[k]
                       + lf_196[k];

            t_221[k] = ab_y[k] * kf_147[k]
                       + lf_197[k];

            t_222[k] = ab_y[k] * kf_148[k]
                       + lf_198[k];

            t_223[k] = ab_y[k] * kf_149[k]
                       + lf_199[k];

            t_224[k] = ab_z[k] * kf_149[k]
                       + lf_209[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, kf_150, kf_151, kf_152, \
                         kf_153, kf_154, lf_150, lf_151, lf_152, lf_153, \
                         lf_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * kf_150[k]
                       + lf_150[k];

            t_226[k] = ab_x[k] * kf_151[k]
                       + lf_151[k];

            t_227[k] = ab_x[k] * kf_152[k]
                       + lf_152[k];

            t_228[k] = ab_x[k] * kf_153[k]
                       + lf_153[k];

            t_229[k] = ab_x[k] * kf_154[k]
                       + lf_154[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, kf_155, kf_156, kf_157, \
                         kf_158, kf_159, lf_155, lf_156, lf_157, lf_158, \
                         lf_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * kf_155[k]
                       + lf_155[k];

            t_231[k] = ab_x[k] * kf_156[k]
                       + lf_156[k];

            t_232[k] = ab_x[k] * kf_157[k]
                       + lf_157[k];

            t_233[k] = ab_x[k] * kf_158[k]
                       + lf_158[k];

            t_234[k] = ab_x[k] * kf_159[k]
                       + lf_159[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, ab_z, kf_156, kf_157, \
                         kf_158, kf_159, lf_216, lf_217, lf_218, lf_219, \
                         lf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_y[k] * kf_156[k]
                       + lf_216[k];

            t_236[k] = ab_y[k] * kf_157[k]
                       + lf_217[k];

            t_237[k] = ab_y[k] * kf_158[k]
                       + lf_218[k];

            t_238[k] = ab_y[k] * kf_159[k]
                       + lf_219[k];

            t_239[k] = ab_z[k] * kf_159[k]
                       + lf_229[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, kf_160, kf_161, kf_162, \
                         kf_163, kf_164, lf_160, lf_161, lf_162, lf_163, \
                         lf_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * kf_160[k]
                       + lf_160[k];

            t_241[k] = ab_x[k] * kf_161[k]
                       + lf_161[k];

            t_242[k] = ab_x[k] * kf_162[k]
                       + lf_162[k];

            t_243[k] = ab_x[k] * kf_163[k]
                       + lf_163[k];

            t_244[k] = ab_x[k] * kf_164[k]
                       + lf_164[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, kf_165, kf_166, kf_167, \
                         kf_168, kf_169, lf_165, lf_166, lf_167, lf_168, \
                         lf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * kf_165[k]
                       + lf_165[k];

            t_246[k] = ab_x[k] * kf_166[k]
                       + lf_166[k];

            t_247[k] = ab_x[k] * kf_167[k]
                       + lf_167[k];

            t_248[k] = ab_x[k] * kf_168[k]
                       + lf_168[k];

            t_249[k] = ab_x[k] * kf_169[k]
                       + lf_169[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, ab_z, kf_166, kf_167, \
                         kf_168, kf_169, lf_226, lf_227, lf_228, lf_229, \
                         lf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_y[k] * kf_166[k]
                       + lf_226[k];

            t_251[k] = ab_y[k] * kf_167[k]
                       + lf_227[k];

            t_252[k] = ab_y[k] * kf_168[k]
                       + lf_228[k];

            t_253[k] = ab_y[k] * kf_169[k]
                       + lf_229[k];

            t_254[k] = ab_z[k] * kf_169[k]
                       + lf_239[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, kf_170, kf_171, kf_172, \
                         kf_173, kf_174, lf_170, lf_171, lf_172, lf_173, \
                         lf_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * kf_170[k]
                       + lf_170[k];

            t_256[k] = ab_x[k] * kf_171[k]
                       + lf_171[k];

            t_257[k] = ab_x[k] * kf_172[k]
                       + lf_172[k];

            t_258[k] = ab_x[k] * kf_173[k]
                       + lf_173[k];

            t_259[k] = ab_x[k] * kf_174[k]
                       + lf_174[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, kf_175, kf_176, kf_177, \
                         kf_178, kf_179, lf_175, lf_176, lf_177, lf_178, \
                         lf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * kf_175[k]
                       + lf_175[k];

            t_261[k] = ab_x[k] * kf_176[k]
                       + lf_176[k];

            t_262[k] = ab_x[k] * kf_177[k]
                       + lf_177[k];

            t_263[k] = ab_x[k] * kf_178[k]
                       + lf_178[k];

            t_264[k] = ab_x[k] * kf_179[k]
                       + lf_179[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, ab_z, kf_176, kf_177, \
                         kf_178, kf_179, lf_236, lf_237, lf_238, lf_239, \
                         lf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_y[k] * kf_176[k]
                       + lf_236[k];

            t_266[k] = ab_y[k] * kf_177[k]
                       + lf_237[k];

            t_267[k] = ab_y[k] * kf_178[k]
                       + lf_238[k];

            t_268[k] = ab_y[k] * kf_179[k]
                       + lf_239[k];

            t_269[k] = ab_z[k] * kf_179[k]
                       + lf_249[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, kf_180, kf_181, kf_182, \
                         kf_183, kf_184, lf_180, lf_181, lf_182, lf_183, \
                         lf_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * kf_180[k]
                       + lf_180[k];

            t_271[k] = ab_x[k] * kf_181[k]
                       + lf_181[k];

            t_272[k] = ab_x[k] * kf_182[k]
                       + lf_182[k];

            t_273[k] = ab_x[k] * kf_183[k]
                       + lf_183[k];

            t_274[k] = ab_x[k] * kf_184[k]
                       + lf_184[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, kf_185, kf_186, kf_187, \
                         kf_188, kf_189, lf_185, lf_186, lf_187, lf_188, \
                         lf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * kf_185[k]
                       + lf_185[k];

            t_276[k] = ab_x[k] * kf_186[k]
                       + lf_186[k];

            t_277[k] = ab_x[k] * kf_187[k]
                       + lf_187[k];

            t_278[k] = ab_x[k] * kf_188[k]
                       + lf_188[k];

            t_279[k] = ab_x[k] * kf_189[k]
                       + lf_189[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, ab_z, kf_186, kf_187, \
                         kf_188, kf_189, lf_246, lf_247, lf_248, lf_249, \
                         lf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_y[k] * kf_186[k]
                       + lf_246[k];

            t_281[k] = ab_y[k] * kf_187[k]
                       + lf_247[k];

            t_282[k] = ab_y[k] * kf_188[k]
                       + lf_248[k];

            t_283[k] = ab_y[k] * kf_189[k]
                       + lf_249[k];

            t_284[k] = ab_z[k] * kf_189[k]
                       + lf_259[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, kf_190, kf_191, kf_192, \
                         kf_193, kf_194, lf_190, lf_191, lf_192, lf_193, \
                         lf_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * kf_190[k]
                       + lf_190[k];

            t_286[k] = ab_x[k] * kf_191[k]
                       + lf_191[k];

            t_287[k] = ab_x[k] * kf_192[k]
                       + lf_192[k];

            t_288[k] = ab_x[k] * kf_193[k]
                       + lf_193[k];

            t_289[k] = ab_x[k] * kf_194[k]
                       + lf_194[k];
        }
    }
}

static auto
compute_hrr_kg_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t kf, const size_t lf, const size_t ncomps,
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
        const auto *ab_z = coordinates.data(8);

        const auto *kf_195 = buffer.data(kf + 195 * ncomps + c);
        const auto *kf_196 = buffer.data(kf + 196 * ncomps + c);
        const auto *kf_197 = buffer.data(kf + 197 * ncomps + c);
        const auto *kf_198 = buffer.data(kf + 198 * ncomps + c);
        const auto *kf_199 = buffer.data(kf + 199 * ncomps + c);
        const auto *kf_200 = buffer.data(kf + 200 * ncomps + c);
        const auto *kf_201 = buffer.data(kf + 201 * ncomps + c);
        const auto *kf_202 = buffer.data(kf + 202 * ncomps + c);
        const auto *kf_203 = buffer.data(kf + 203 * ncomps + c);
        const auto *kf_204 = buffer.data(kf + 204 * ncomps + c);
        const auto *kf_205 = buffer.data(kf + 205 * ncomps + c);
        const auto *kf_206 = buffer.data(kf + 206 * ncomps + c);
        const auto *kf_207 = buffer.data(kf + 207 * ncomps + c);
        const auto *kf_208 = buffer.data(kf + 208 * ncomps + c);
        const auto *kf_209 = buffer.data(kf + 209 * ncomps + c);
        const auto *kf_210 = buffer.data(kf + 210 * ncomps + c);
        const auto *kf_211 = buffer.data(kf + 211 * ncomps + c);
        const auto *kf_212 = buffer.data(kf + 212 * ncomps + c);
        const auto *kf_213 = buffer.data(kf + 213 * ncomps + c);
        const auto *kf_214 = buffer.data(kf + 214 * ncomps + c);
        const auto *kf_215 = buffer.data(kf + 215 * ncomps + c);
        const auto *kf_216 = buffer.data(kf + 216 * ncomps + c);
        const auto *kf_217 = buffer.data(kf + 217 * ncomps + c);
        const auto *kf_218 = buffer.data(kf + 218 * ncomps + c);
        const auto *kf_219 = buffer.data(kf + 219 * ncomps + c);
        const auto *kf_220 = buffer.data(kf + 220 * ncomps + c);
        const auto *kf_221 = buffer.data(kf + 221 * ncomps + c);
        const auto *kf_222 = buffer.data(kf + 222 * ncomps + c);
        const auto *kf_223 = buffer.data(kf + 223 * ncomps + c);
        const auto *kf_224 = buffer.data(kf + 224 * ncomps + c);
        const auto *kf_225 = buffer.data(kf + 225 * ncomps + c);
        const auto *kf_226 = buffer.data(kf + 226 * ncomps + c);
        const auto *kf_227 = buffer.data(kf + 227 * ncomps + c);
        const auto *kf_228 = buffer.data(kf + 228 * ncomps + c);
        const auto *kf_229 = buffer.data(kf + 229 * ncomps + c);
        const auto *kf_230 = buffer.data(kf + 230 * ncomps + c);
        const auto *kf_231 = buffer.data(kf + 231 * ncomps + c);
        const auto *kf_232 = buffer.data(kf + 232 * ncomps + c);
        const auto *kf_233 = buffer.data(kf + 233 * ncomps + c);
        const auto *kf_234 = buffer.data(kf + 234 * ncomps + c);
        const auto *kf_235 = buffer.data(kf + 235 * ncomps + c);
        const auto *kf_236 = buffer.data(kf + 236 * ncomps + c);
        const auto *kf_237 = buffer.data(kf + 237 * ncomps + c);
        const auto *kf_238 = buffer.data(kf + 238 * ncomps + c);
        const auto *kf_239 = buffer.data(kf + 239 * ncomps + c);
        const auto *kf_240 = buffer.data(kf + 240 * ncomps + c);
        const auto *kf_241 = buffer.data(kf + 241 * ncomps + c);
        const auto *kf_242 = buffer.data(kf + 242 * ncomps + c);
        const auto *kf_243 = buffer.data(kf + 243 * ncomps + c);
        const auto *kf_244 = buffer.data(kf + 244 * ncomps + c);
        const auto *kf_245 = buffer.data(kf + 245 * ncomps + c);
        const auto *kf_246 = buffer.data(kf + 246 * ncomps + c);
        const auto *kf_247 = buffer.data(kf + 247 * ncomps + c);
        const auto *kf_248 = buffer.data(kf + 248 * ncomps + c);
        const auto *kf_249 = buffer.data(kf + 249 * ncomps + c);
        const auto *kf_250 = buffer.data(kf + 250 * ncomps + c);
        const auto *kf_251 = buffer.data(kf + 251 * ncomps + c);
        const auto *kf_252 = buffer.data(kf + 252 * ncomps + c);
        const auto *kf_253 = buffer.data(kf + 253 * ncomps + c);
        const auto *kf_254 = buffer.data(kf + 254 * ncomps + c);
        const auto *kf_255 = buffer.data(kf + 255 * ncomps + c);
        const auto *kf_256 = buffer.data(kf + 256 * ncomps + c);
        const auto *kf_257 = buffer.data(kf + 257 * ncomps + c);
        const auto *kf_258 = buffer.data(kf + 258 * ncomps + c);
        const auto *kf_259 = buffer.data(kf + 259 * ncomps + c);
        const auto *kf_260 = buffer.data(kf + 260 * ncomps + c);
        const auto *kf_261 = buffer.data(kf + 261 * ncomps + c);
        const auto *kf_262 = buffer.data(kf + 262 * ncomps + c);
        const auto *kf_263 = buffer.data(kf + 263 * ncomps + c);
        const auto *kf_264 = buffer.data(kf + 264 * ncomps + c);
        const auto *kf_265 = buffer.data(kf + 265 * ncomps + c);
        const auto *kf_266 = buffer.data(kf + 266 * ncomps + c);
        const auto *kf_267 = buffer.data(kf + 267 * ncomps + c);
        const auto *kf_268 = buffer.data(kf + 268 * ncomps + c);
        const auto *kf_269 = buffer.data(kf + 269 * ncomps + c);
        const auto *kf_270 = buffer.data(kf + 270 * ncomps + c);
        const auto *kf_271 = buffer.data(kf + 271 * ncomps + c);
        const auto *kf_272 = buffer.data(kf + 272 * ncomps + c);
        const auto *kf_273 = buffer.data(kf + 273 * ncomps + c);
        const auto *kf_274 = buffer.data(kf + 274 * ncomps + c);
        const auto *kf_275 = buffer.data(kf + 275 * ncomps + c);
        const auto *kf_276 = buffer.data(kf + 276 * ncomps + c);
        const auto *kf_277 = buffer.data(kf + 277 * ncomps + c);
        const auto *kf_278 = buffer.data(kf + 278 * ncomps + c);
        const auto *kf_279 = buffer.data(kf + 279 * ncomps + c);
        const auto *kf_280 = buffer.data(kf + 280 * ncomps + c);
        const auto *kf_281 = buffer.data(kf + 281 * ncomps + c);
        const auto *kf_282 = buffer.data(kf + 282 * ncomps + c);
        const auto *kf_283 = buffer.data(kf + 283 * ncomps + c);
        const auto *kf_284 = buffer.data(kf + 284 * ncomps + c);
        const auto *kf_285 = buffer.data(kf + 285 * ncomps + c);
        const auto *kf_286 = buffer.data(kf + 286 * ncomps + c);
        const auto *kf_287 = buffer.data(kf + 287 * ncomps + c);
        const auto *kf_288 = buffer.data(kf + 288 * ncomps + c);
        const auto *kf_289 = buffer.data(kf + 289 * ncomps + c);

        const auto *lf_195 = buffer.data(lf + 195 * ncomps + c);
        const auto *lf_196 = buffer.data(lf + 196 * ncomps + c);
        const auto *lf_197 = buffer.data(lf + 197 * ncomps + c);
        const auto *lf_198 = buffer.data(lf + 198 * ncomps + c);
        const auto *lf_199 = buffer.data(lf + 199 * ncomps + c);
        const auto *lf_200 = buffer.data(lf + 200 * ncomps + c);
        const auto *lf_201 = buffer.data(lf + 201 * ncomps + c);
        const auto *lf_202 = buffer.data(lf + 202 * ncomps + c);
        const auto *lf_203 = buffer.data(lf + 203 * ncomps + c);
        const auto *lf_204 = buffer.data(lf + 204 * ncomps + c);
        const auto *lf_205 = buffer.data(lf + 205 * ncomps + c);
        const auto *lf_206 = buffer.data(lf + 206 * ncomps + c);
        const auto *lf_207 = buffer.data(lf + 207 * ncomps + c);
        const auto *lf_208 = buffer.data(lf + 208 * ncomps + c);
        const auto *lf_209 = buffer.data(lf + 209 * ncomps + c);
        const auto *lf_210 = buffer.data(lf + 210 * ncomps + c);
        const auto *lf_211 = buffer.data(lf + 211 * ncomps + c);
        const auto *lf_212 = buffer.data(lf + 212 * ncomps + c);
        const auto *lf_213 = buffer.data(lf + 213 * ncomps + c);
        const auto *lf_214 = buffer.data(lf + 214 * ncomps + c);
        const auto *lf_215 = buffer.data(lf + 215 * ncomps + c);
        const auto *lf_216 = buffer.data(lf + 216 * ncomps + c);
        const auto *lf_217 = buffer.data(lf + 217 * ncomps + c);
        const auto *lf_218 = buffer.data(lf + 218 * ncomps + c);
        const auto *lf_219 = buffer.data(lf + 219 * ncomps + c);
        const auto *lf_220 = buffer.data(lf + 220 * ncomps + c);
        const auto *lf_221 = buffer.data(lf + 221 * ncomps + c);
        const auto *lf_222 = buffer.data(lf + 222 * ncomps + c);
        const auto *lf_223 = buffer.data(lf + 223 * ncomps + c);
        const auto *lf_224 = buffer.data(lf + 224 * ncomps + c);
        const auto *lf_225 = buffer.data(lf + 225 * ncomps + c);
        const auto *lf_226 = buffer.data(lf + 226 * ncomps + c);
        const auto *lf_227 = buffer.data(lf + 227 * ncomps + c);
        const auto *lf_228 = buffer.data(lf + 228 * ncomps + c);
        const auto *lf_229 = buffer.data(lf + 229 * ncomps + c);
        const auto *lf_230 = buffer.data(lf + 230 * ncomps + c);
        const auto *lf_231 = buffer.data(lf + 231 * ncomps + c);
        const auto *lf_232 = buffer.data(lf + 232 * ncomps + c);
        const auto *lf_233 = buffer.data(lf + 233 * ncomps + c);
        const auto *lf_234 = buffer.data(lf + 234 * ncomps + c);
        const auto *lf_235 = buffer.data(lf + 235 * ncomps + c);
        const auto *lf_236 = buffer.data(lf + 236 * ncomps + c);
        const auto *lf_237 = buffer.data(lf + 237 * ncomps + c);
        const auto *lf_238 = buffer.data(lf + 238 * ncomps + c);
        const auto *lf_239 = buffer.data(lf + 239 * ncomps + c);
        const auto *lf_240 = buffer.data(lf + 240 * ncomps + c);
        const auto *lf_241 = buffer.data(lf + 241 * ncomps + c);
        const auto *lf_242 = buffer.data(lf + 242 * ncomps + c);
        const auto *lf_243 = buffer.data(lf + 243 * ncomps + c);
        const auto *lf_244 = buffer.data(lf + 244 * ncomps + c);
        const auto *lf_245 = buffer.data(lf + 245 * ncomps + c);
        const auto *lf_246 = buffer.data(lf + 246 * ncomps + c);
        const auto *lf_247 = buffer.data(lf + 247 * ncomps + c);
        const auto *lf_248 = buffer.data(lf + 248 * ncomps + c);
        const auto *lf_249 = buffer.data(lf + 249 * ncomps + c);
        const auto *lf_250 = buffer.data(lf + 250 * ncomps + c);
        const auto *lf_251 = buffer.data(lf + 251 * ncomps + c);
        const auto *lf_252 = buffer.data(lf + 252 * ncomps + c);
        const auto *lf_253 = buffer.data(lf + 253 * ncomps + c);
        const auto *lf_254 = buffer.data(lf + 254 * ncomps + c);
        const auto *lf_255 = buffer.data(lf + 255 * ncomps + c);
        const auto *lf_256 = buffer.data(lf + 256 * ncomps + c);
        const auto *lf_257 = buffer.data(lf + 257 * ncomps + c);
        const auto *lf_258 = buffer.data(lf + 258 * ncomps + c);
        const auto *lf_259 = buffer.data(lf + 259 * ncomps + c);
        const auto *lf_260 = buffer.data(lf + 260 * ncomps + c);
        const auto *lf_261 = buffer.data(lf + 261 * ncomps + c);
        const auto *lf_262 = buffer.data(lf + 262 * ncomps + c);
        const auto *lf_263 = buffer.data(lf + 263 * ncomps + c);
        const auto *lf_264 = buffer.data(lf + 264 * ncomps + c);
        const auto *lf_265 = buffer.data(lf + 265 * ncomps + c);
        const auto *lf_266 = buffer.data(lf + 266 * ncomps + c);
        const auto *lf_267 = buffer.data(lf + 267 * ncomps + c);
        const auto *lf_268 = buffer.data(lf + 268 * ncomps + c);
        const auto *lf_269 = buffer.data(lf + 269 * ncomps + c);
        const auto *lf_270 = buffer.data(lf + 270 * ncomps + c);
        const auto *lf_271 = buffer.data(lf + 271 * ncomps + c);
        const auto *lf_272 = buffer.data(lf + 272 * ncomps + c);
        const auto *lf_273 = buffer.data(lf + 273 * ncomps + c);
        const auto *lf_274 = buffer.data(lf + 274 * ncomps + c);
        const auto *lf_275 = buffer.data(lf + 275 * ncomps + c);
        const auto *lf_276 = buffer.data(lf + 276 * ncomps + c);
        const auto *lf_277 = buffer.data(lf + 277 * ncomps + c);
        const auto *lf_278 = buffer.data(lf + 278 * ncomps + c);
        const auto *lf_279 = buffer.data(lf + 279 * ncomps + c);
        const auto *lf_280 = buffer.data(lf + 280 * ncomps + c);
        const auto *lf_281 = buffer.data(lf + 281 * ncomps + c);
        const auto *lf_282 = buffer.data(lf + 282 * ncomps + c);
        const auto *lf_283 = buffer.data(lf + 283 * ncomps + c);
        const auto *lf_284 = buffer.data(lf + 284 * ncomps + c);
        const auto *lf_285 = buffer.data(lf + 285 * ncomps + c);
        const auto *lf_286 = buffer.data(lf + 286 * ncomps + c);
        const auto *lf_287 = buffer.data(lf + 287 * ncomps + c);
        const auto *lf_288 = buffer.data(lf + 288 * ncomps + c);
        const auto *lf_289 = buffer.data(lf + 289 * ncomps + c);
        const auto *lf_296 = buffer.data(lf + 296 * ncomps + c);
        const auto *lf_297 = buffer.data(lf + 297 * ncomps + c);
        const auto *lf_298 = buffer.data(lf + 298 * ncomps + c);
        const auto *lf_299 = buffer.data(lf + 299 * ncomps + c);
        const auto *lf_306 = buffer.data(lf + 306 * ncomps + c);
        const auto *lf_307 = buffer.data(lf + 307 * ncomps + c);
        const auto *lf_308 = buffer.data(lf + 308 * ncomps + c);
        const auto *lf_309 = buffer.data(lf + 309 * ncomps + c);
        const auto *lf_316 = buffer.data(lf + 316 * ncomps + c);
        const auto *lf_317 = buffer.data(lf + 317 * ncomps + c);
        const auto *lf_318 = buffer.data(lf + 318 * ncomps + c);
        const auto *lf_319 = buffer.data(lf + 319 * ncomps + c);
        const auto *lf_326 = buffer.data(lf + 326 * ncomps + c);
        const auto *lf_327 = buffer.data(lf + 327 * ncomps + c);
        const auto *lf_328 = buffer.data(lf + 328 * ncomps + c);
        const auto *lf_329 = buffer.data(lf + 329 * ncomps + c);
        const auto *lf_336 = buffer.data(lf + 336 * ncomps + c);
        const auto *lf_337 = buffer.data(lf + 337 * ncomps + c);
        const auto *lf_338 = buffer.data(lf + 338 * ncomps + c);
        const auto *lf_339 = buffer.data(lf + 339 * ncomps + c);
        const auto *lf_346 = buffer.data(lf + 346 * ncomps + c);
        const auto *lf_347 = buffer.data(lf + 347 * ncomps + c);
        const auto *lf_348 = buffer.data(lf + 348 * ncomps + c);
        const auto *lf_349 = buffer.data(lf + 349 * ncomps + c);
        const auto *lf_359 = buffer.data(lf + 359 * ncomps + c);
        const auto *lf_366 = buffer.data(lf + 366 * ncomps + c);
        const auto *lf_367 = buffer.data(lf + 367 * ncomps + c);
        const auto *lf_368 = buffer.data(lf + 368 * ncomps + c);
        const auto *lf_369 = buffer.data(lf + 369 * ncomps + c);
        const auto *lf_379 = buffer.data(lf + 379 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, kf_195, kf_196, kf_197, \
                         kf_198, kf_199, lf_195, lf_196, lf_197, lf_198, \
                         lf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * kf_195[k]
                       + lf_195[k];

            t_291[k] = ab_x[k] * kf_196[k]
                       + lf_196[k];

            t_292[k] = ab_x[k] * kf_197[k]
                       + lf_197[k];

            t_293[k] = ab_x[k] * kf_198[k]
                       + lf_198[k];

            t_294[k] = ab_x[k] * kf_199[k]
                       + lf_199[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, ab_z, kf_196, kf_197, \
                         kf_198, kf_199, lf_256, lf_257, lf_258, lf_259, \
                         lf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_y[k] * kf_196[k]
                       + lf_256[k];

            t_296[k] = ab_y[k] * kf_197[k]
                       + lf_257[k];

            t_297[k] = ab_y[k] * kf_198[k]
                       + lf_258[k];

            t_298[k] = ab_y[k] * kf_199[k]
                       + lf_259[k];

            t_299[k] = ab_z[k] * kf_199[k]
                       + lf_269[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, kf_200, kf_201, kf_202, \
                         kf_203, kf_204, lf_200, lf_201, lf_202, lf_203, \
                         lf_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * kf_200[k]
                       + lf_200[k];

            t_301[k] = ab_x[k] * kf_201[k]
                       + lf_201[k];

            t_302[k] = ab_x[k] * kf_202[k]
                       + lf_202[k];

            t_303[k] = ab_x[k] * kf_203[k]
                       + lf_203[k];

            t_304[k] = ab_x[k] * kf_204[k]
                       + lf_204[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, kf_205, kf_206, kf_207, \
                         kf_208, kf_209, lf_205, lf_206, lf_207, lf_208, \
                         lf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * kf_205[k]
                       + lf_205[k];

            t_306[k] = ab_x[k] * kf_206[k]
                       + lf_206[k];

            t_307[k] = ab_x[k] * kf_207[k]
                       + lf_207[k];

            t_308[k] = ab_x[k] * kf_208[k]
                       + lf_208[k];

            t_309[k] = ab_x[k] * kf_209[k]
                       + lf_209[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, ab_z, kf_206, kf_207, \
                         kf_208, kf_209, lf_266, lf_267, lf_268, lf_269, \
                         lf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_y[k] * kf_206[k]
                       + lf_266[k];

            t_311[k] = ab_y[k] * kf_207[k]
                       + lf_267[k];

            t_312[k] = ab_y[k] * kf_208[k]
                       + lf_268[k];

            t_313[k] = ab_y[k] * kf_209[k]
                       + lf_269[k];

            t_314[k] = ab_z[k] * kf_209[k]
                       + lf_279[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, kf_210, kf_211, kf_212, \
                         kf_213, kf_214, lf_210, lf_211, lf_212, lf_213, \
                         lf_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * kf_210[k]
                       + lf_210[k];

            t_316[k] = ab_x[k] * kf_211[k]
                       + lf_211[k];

            t_317[k] = ab_x[k] * kf_212[k]
                       + lf_212[k];

            t_318[k] = ab_x[k] * kf_213[k]
                       + lf_213[k];

            t_319[k] = ab_x[k] * kf_214[k]
                       + lf_214[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, kf_215, kf_216, kf_217, \
                         kf_218, kf_219, lf_215, lf_216, lf_217, lf_218, \
                         lf_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * kf_215[k]
                       + lf_215[k];

            t_321[k] = ab_x[k] * kf_216[k]
                       + lf_216[k];

            t_322[k] = ab_x[k] * kf_217[k]
                       + lf_217[k];

            t_323[k] = ab_x[k] * kf_218[k]
                       + lf_218[k];

            t_324[k] = ab_x[k] * kf_219[k]
                       + lf_219[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, ab_z, kf_216, kf_217, \
                         kf_218, kf_219, lf_286, lf_287, lf_288, lf_289, \
                         lf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_y[k] * kf_216[k]
                       + lf_286[k];

            t_326[k] = ab_y[k] * kf_217[k]
                       + lf_287[k];

            t_327[k] = ab_y[k] * kf_218[k]
                       + lf_288[k];

            t_328[k] = ab_y[k] * kf_219[k]
                       + lf_289[k];

            t_329[k] = ab_z[k] * kf_219[k]
                       + lf_299[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, kf_220, kf_221, kf_222, \
                         kf_223, kf_224, lf_220, lf_221, lf_222, lf_223, \
                         lf_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * kf_220[k]
                       + lf_220[k];

            t_331[k] = ab_x[k] * kf_221[k]
                       + lf_221[k];

            t_332[k] = ab_x[k] * kf_222[k]
                       + lf_222[k];

            t_333[k] = ab_x[k] * kf_223[k]
                       + lf_223[k];

            t_334[k] = ab_x[k] * kf_224[k]
                       + lf_224[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, kf_225, kf_226, kf_227, \
                         kf_228, kf_229, lf_225, lf_226, lf_227, lf_228, \
                         lf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * kf_225[k]
                       + lf_225[k];

            t_336[k] = ab_x[k] * kf_226[k]
                       + lf_226[k];

            t_337[k] = ab_x[k] * kf_227[k]
                       + lf_227[k];

            t_338[k] = ab_x[k] * kf_228[k]
                       + lf_228[k];

            t_339[k] = ab_x[k] * kf_229[k]
                       + lf_229[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, ab_z, kf_226, kf_227, \
                         kf_228, kf_229, lf_296, lf_297, lf_298, lf_299, \
                         lf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_y[k] * kf_226[k]
                       + lf_296[k];

            t_341[k] = ab_y[k] * kf_227[k]
                       + lf_297[k];

            t_342[k] = ab_y[k] * kf_228[k]
                       + lf_298[k];

            t_343[k] = ab_y[k] * kf_229[k]
                       + lf_299[k];

            t_344[k] = ab_z[k] * kf_229[k]
                       + lf_309[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, kf_230, kf_231, kf_232, \
                         kf_233, kf_234, lf_230, lf_231, lf_232, lf_233, \
                         lf_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * kf_230[k]
                       + lf_230[k];

            t_346[k] = ab_x[k] * kf_231[k]
                       + lf_231[k];

            t_347[k] = ab_x[k] * kf_232[k]
                       + lf_232[k];

            t_348[k] = ab_x[k] * kf_233[k]
                       + lf_233[k];

            t_349[k] = ab_x[k] * kf_234[k]
                       + lf_234[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, kf_235, kf_236, kf_237, \
                         kf_238, kf_239, lf_235, lf_236, lf_237, lf_238, \
                         lf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * kf_235[k]
                       + lf_235[k];

            t_351[k] = ab_x[k] * kf_236[k]
                       + lf_236[k];

            t_352[k] = ab_x[k] * kf_237[k]
                       + lf_237[k];

            t_353[k] = ab_x[k] * kf_238[k]
                       + lf_238[k];

            t_354[k] = ab_x[k] * kf_239[k]
                       + lf_239[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, ab_z, kf_236, kf_237, \
                         kf_238, kf_239, lf_306, lf_307, lf_308, lf_309, \
                         lf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_y[k] * kf_236[k]
                       + lf_306[k];

            t_356[k] = ab_y[k] * kf_237[k]
                       + lf_307[k];

            t_357[k] = ab_y[k] * kf_238[k]
                       + lf_308[k];

            t_358[k] = ab_y[k] * kf_239[k]
                       + lf_309[k];

            t_359[k] = ab_z[k] * kf_239[k]
                       + lf_319[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, kf_240, kf_241, kf_242, \
                         kf_243, kf_244, lf_240, lf_241, lf_242, lf_243, \
                         lf_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * kf_240[k]
                       + lf_240[k];

            t_361[k] = ab_x[k] * kf_241[k]
                       + lf_241[k];

            t_362[k] = ab_x[k] * kf_242[k]
                       + lf_242[k];

            t_363[k] = ab_x[k] * kf_243[k]
                       + lf_243[k];

            t_364[k] = ab_x[k] * kf_244[k]
                       + lf_244[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, kf_245, kf_246, kf_247, \
                         kf_248, kf_249, lf_245, lf_246, lf_247, lf_248, \
                         lf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * kf_245[k]
                       + lf_245[k];

            t_366[k] = ab_x[k] * kf_246[k]
                       + lf_246[k];

            t_367[k] = ab_x[k] * kf_247[k]
                       + lf_247[k];

            t_368[k] = ab_x[k] * kf_248[k]
                       + lf_248[k];

            t_369[k] = ab_x[k] * kf_249[k]
                       + lf_249[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, ab_z, kf_246, kf_247, \
                         kf_248, kf_249, lf_316, lf_317, lf_318, lf_319, \
                         lf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_y[k] * kf_246[k]
                       + lf_316[k];

            t_371[k] = ab_y[k] * kf_247[k]
                       + lf_317[k];

            t_372[k] = ab_y[k] * kf_248[k]
                       + lf_318[k];

            t_373[k] = ab_y[k] * kf_249[k]
                       + lf_319[k];

            t_374[k] = ab_z[k] * kf_249[k]
                       + lf_329[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, kf_250, kf_251, kf_252, \
                         kf_253, kf_254, lf_250, lf_251, lf_252, lf_253, \
                         lf_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * kf_250[k]
                       + lf_250[k];

            t_376[k] = ab_x[k] * kf_251[k]
                       + lf_251[k];

            t_377[k] = ab_x[k] * kf_252[k]
                       + lf_252[k];

            t_378[k] = ab_x[k] * kf_253[k]
                       + lf_253[k];

            t_379[k] = ab_x[k] * kf_254[k]
                       + lf_254[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, kf_255, kf_256, kf_257, \
                         kf_258, kf_259, lf_255, lf_256, lf_257, lf_258, \
                         lf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * kf_255[k]
                       + lf_255[k];

            t_381[k] = ab_x[k] * kf_256[k]
                       + lf_256[k];

            t_382[k] = ab_x[k] * kf_257[k]
                       + lf_257[k];

            t_383[k] = ab_x[k] * kf_258[k]
                       + lf_258[k];

            t_384[k] = ab_x[k] * kf_259[k]
                       + lf_259[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, ab_z, kf_256, kf_257, \
                         kf_258, kf_259, lf_326, lf_327, lf_328, lf_329, \
                         lf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_y[k] * kf_256[k]
                       + lf_326[k];

            t_386[k] = ab_y[k] * kf_257[k]
                       + lf_327[k];

            t_387[k] = ab_y[k] * kf_258[k]
                       + lf_328[k];

            t_388[k] = ab_y[k] * kf_259[k]
                       + lf_329[k];

            t_389[k] = ab_z[k] * kf_259[k]
                       + lf_339[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, kf_260, kf_261, kf_262, \
                         kf_263, kf_264, lf_260, lf_261, lf_262, lf_263, \
                         lf_264 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * kf_260[k]
                       + lf_260[k];

            t_391[k] = ab_x[k] * kf_261[k]
                       + lf_261[k];

            t_392[k] = ab_x[k] * kf_262[k]
                       + lf_262[k];

            t_393[k] = ab_x[k] * kf_263[k]
                       + lf_263[k];

            t_394[k] = ab_x[k] * kf_264[k]
                       + lf_264[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, kf_265, kf_266, kf_267, \
                         kf_268, kf_269, lf_265, lf_266, lf_267, lf_268, \
                         lf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * kf_265[k]
                       + lf_265[k];

            t_396[k] = ab_x[k] * kf_266[k]
                       + lf_266[k];

            t_397[k] = ab_x[k] * kf_267[k]
                       + lf_267[k];

            t_398[k] = ab_x[k] * kf_268[k]
                       + lf_268[k];

            t_399[k] = ab_x[k] * kf_269[k]
                       + lf_269[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, ab_z, kf_266, kf_267, \
                         kf_268, kf_269, lf_336, lf_337, lf_338, lf_339, \
                         lf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_y[k] * kf_266[k]
                       + lf_336[k];

            t_401[k] = ab_y[k] * kf_267[k]
                       + lf_337[k];

            t_402[k] = ab_y[k] * kf_268[k]
                       + lf_338[k];

            t_403[k] = ab_y[k] * kf_269[k]
                       + lf_339[k];

            t_404[k] = ab_z[k] * kf_269[k]
                       + lf_349[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, kf_270, kf_271, kf_272, \
                         kf_273, kf_274, lf_270, lf_271, lf_272, lf_273, \
                         lf_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * kf_270[k]
                       + lf_270[k];

            t_406[k] = ab_x[k] * kf_271[k]
                       + lf_271[k];

            t_407[k] = ab_x[k] * kf_272[k]
                       + lf_272[k];

            t_408[k] = ab_x[k] * kf_273[k]
                       + lf_273[k];

            t_409[k] = ab_x[k] * kf_274[k]
                       + lf_274[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, kf_275, kf_276, kf_277, \
                         kf_278, kf_279, lf_275, lf_276, lf_277, lf_278, \
                         lf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * kf_275[k]
                       + lf_275[k];

            t_411[k] = ab_x[k] * kf_276[k]
                       + lf_276[k];

            t_412[k] = ab_x[k] * kf_277[k]
                       + lf_277[k];

            t_413[k] = ab_x[k] * kf_278[k]
                       + lf_278[k];

            t_414[k] = ab_x[k] * kf_279[k]
                       + lf_279[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, ab_z, kf_276, kf_277, \
                         kf_278, kf_279, lf_346, lf_347, lf_348, lf_349, \
                         lf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_y[k] * kf_276[k]
                       + lf_346[k];

            t_416[k] = ab_y[k] * kf_277[k]
                       + lf_347[k];

            t_417[k] = ab_y[k] * kf_278[k]
                       + lf_348[k];

            t_418[k] = ab_y[k] * kf_279[k]
                       + lf_349[k];

            t_419[k] = ab_z[k] * kf_279[k]
                       + lf_359[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, kf_280, kf_281, kf_282, \
                         kf_283, kf_284, lf_280, lf_281, lf_282, lf_283, \
                         lf_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * kf_280[k]
                       + lf_280[k];

            t_421[k] = ab_x[k] * kf_281[k]
                       + lf_281[k];

            t_422[k] = ab_x[k] * kf_282[k]
                       + lf_282[k];

            t_423[k] = ab_x[k] * kf_283[k]
                       + lf_283[k];

            t_424[k] = ab_x[k] * kf_284[k]
                       + lf_284[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, kf_285, kf_286, kf_287, \
                         kf_288, kf_289, lf_285, lf_286, lf_287, lf_288, \
                         lf_289 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * kf_285[k]
                       + lf_285[k];

            t_426[k] = ab_x[k] * kf_286[k]
                       + lf_286[k];

            t_427[k] = ab_x[k] * kf_287[k]
                       + lf_287[k];

            t_428[k] = ab_x[k] * kf_288[k]
                       + lf_288[k];

            t_429[k] = ab_x[k] * kf_289[k]
                       + lf_289[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_y, ab_z, kf_286, kf_287, \
                         kf_288, kf_289, lf_366, lf_367, lf_368, lf_369, \
                         lf_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_y[k] * kf_286[k]
                       + lf_366[k];

            t_431[k] = ab_y[k] * kf_287[k]
                       + lf_367[k];

            t_432[k] = ab_y[k] * kf_288[k]
                       + lf_368[k];

            t_433[k] = ab_y[k] * kf_289[k]
                       + lf_369[k];

            t_434[k] = ab_z[k] * kf_289[k]
                       + lf_379[k];
        }
    }
}

static auto
compute_hrr_kg_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t kf, const size_t lf, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kf_290 = buffer.data(kf + 290 * ncomps + c);
        const auto *kf_291 = buffer.data(kf + 291 * ncomps + c);
        const auto *kf_292 = buffer.data(kf + 292 * ncomps + c);
        const auto *kf_293 = buffer.data(kf + 293 * ncomps + c);
        const auto *kf_294 = buffer.data(kf + 294 * ncomps + c);
        const auto *kf_295 = buffer.data(kf + 295 * ncomps + c);
        const auto *kf_296 = buffer.data(kf + 296 * ncomps + c);
        const auto *kf_297 = buffer.data(kf + 297 * ncomps + c);
        const auto *kf_298 = buffer.data(kf + 298 * ncomps + c);
        const auto *kf_299 = buffer.data(kf + 299 * ncomps + c);
        const auto *kf_300 = buffer.data(kf + 300 * ncomps + c);
        const auto *kf_301 = buffer.data(kf + 301 * ncomps + c);
        const auto *kf_302 = buffer.data(kf + 302 * ncomps + c);
        const auto *kf_303 = buffer.data(kf + 303 * ncomps + c);
        const auto *kf_304 = buffer.data(kf + 304 * ncomps + c);
        const auto *kf_305 = buffer.data(kf + 305 * ncomps + c);
        const auto *kf_306 = buffer.data(kf + 306 * ncomps + c);
        const auto *kf_307 = buffer.data(kf + 307 * ncomps + c);
        const auto *kf_308 = buffer.data(kf + 308 * ncomps + c);
        const auto *kf_309 = buffer.data(kf + 309 * ncomps + c);
        const auto *kf_310 = buffer.data(kf + 310 * ncomps + c);
        const auto *kf_311 = buffer.data(kf + 311 * ncomps + c);
        const auto *kf_312 = buffer.data(kf + 312 * ncomps + c);
        const auto *kf_313 = buffer.data(kf + 313 * ncomps + c);
        const auto *kf_314 = buffer.data(kf + 314 * ncomps + c);
        const auto *kf_315 = buffer.data(kf + 315 * ncomps + c);
        const auto *kf_316 = buffer.data(kf + 316 * ncomps + c);
        const auto *kf_317 = buffer.data(kf + 317 * ncomps + c);
        const auto *kf_318 = buffer.data(kf + 318 * ncomps + c);
        const auto *kf_319 = buffer.data(kf + 319 * ncomps + c);
        const auto *kf_320 = buffer.data(kf + 320 * ncomps + c);
        const auto *kf_321 = buffer.data(kf + 321 * ncomps + c);
        const auto *kf_322 = buffer.data(kf + 322 * ncomps + c);
        const auto *kf_323 = buffer.data(kf + 323 * ncomps + c);
        const auto *kf_324 = buffer.data(kf + 324 * ncomps + c);
        const auto *kf_325 = buffer.data(kf + 325 * ncomps + c);
        const auto *kf_326 = buffer.data(kf + 326 * ncomps + c);
        const auto *kf_327 = buffer.data(kf + 327 * ncomps + c);
        const auto *kf_328 = buffer.data(kf + 328 * ncomps + c);
        const auto *kf_329 = buffer.data(kf + 329 * ncomps + c);
        const auto *kf_330 = buffer.data(kf + 330 * ncomps + c);
        const auto *kf_331 = buffer.data(kf + 331 * ncomps + c);
        const auto *kf_332 = buffer.data(kf + 332 * ncomps + c);
        const auto *kf_333 = buffer.data(kf + 333 * ncomps + c);
        const auto *kf_334 = buffer.data(kf + 334 * ncomps + c);
        const auto *kf_335 = buffer.data(kf + 335 * ncomps + c);
        const auto *kf_336 = buffer.data(kf + 336 * ncomps + c);
        const auto *kf_337 = buffer.data(kf + 337 * ncomps + c);
        const auto *kf_338 = buffer.data(kf + 338 * ncomps + c);
        const auto *kf_339 = buffer.data(kf + 339 * ncomps + c);
        const auto *kf_340 = buffer.data(kf + 340 * ncomps + c);
        const auto *kf_341 = buffer.data(kf + 341 * ncomps + c);
        const auto *kf_342 = buffer.data(kf + 342 * ncomps + c);
        const auto *kf_343 = buffer.data(kf + 343 * ncomps + c);
        const auto *kf_344 = buffer.data(kf + 344 * ncomps + c);
        const auto *kf_345 = buffer.data(kf + 345 * ncomps + c);
        const auto *kf_346 = buffer.data(kf + 346 * ncomps + c);
        const auto *kf_347 = buffer.data(kf + 347 * ncomps + c);
        const auto *kf_348 = buffer.data(kf + 348 * ncomps + c);
        const auto *kf_349 = buffer.data(kf + 349 * ncomps + c);
        const auto *kf_350 = buffer.data(kf + 350 * ncomps + c);
        const auto *kf_351 = buffer.data(kf + 351 * ncomps + c);
        const auto *kf_352 = buffer.data(kf + 352 * ncomps + c);
        const auto *kf_353 = buffer.data(kf + 353 * ncomps + c);
        const auto *kf_354 = buffer.data(kf + 354 * ncomps + c);
        const auto *kf_355 = buffer.data(kf + 355 * ncomps + c);
        const auto *kf_356 = buffer.data(kf + 356 * ncomps + c);
        const auto *kf_357 = buffer.data(kf + 357 * ncomps + c);
        const auto *kf_358 = buffer.data(kf + 358 * ncomps + c);
        const auto *kf_359 = buffer.data(kf + 359 * ncomps + c);

        const auto *lf_290 = buffer.data(lf + 290 * ncomps + c);
        const auto *lf_291 = buffer.data(lf + 291 * ncomps + c);
        const auto *lf_292 = buffer.data(lf + 292 * ncomps + c);
        const auto *lf_293 = buffer.data(lf + 293 * ncomps + c);
        const auto *lf_294 = buffer.data(lf + 294 * ncomps + c);
        const auto *lf_295 = buffer.data(lf + 295 * ncomps + c);
        const auto *lf_296 = buffer.data(lf + 296 * ncomps + c);
        const auto *lf_297 = buffer.data(lf + 297 * ncomps + c);
        const auto *lf_298 = buffer.data(lf + 298 * ncomps + c);
        const auto *lf_299 = buffer.data(lf + 299 * ncomps + c);
        const auto *lf_300 = buffer.data(lf + 300 * ncomps + c);
        const auto *lf_301 = buffer.data(lf + 301 * ncomps + c);
        const auto *lf_302 = buffer.data(lf + 302 * ncomps + c);
        const auto *lf_303 = buffer.data(lf + 303 * ncomps + c);
        const auto *lf_304 = buffer.data(lf + 304 * ncomps + c);
        const auto *lf_305 = buffer.data(lf + 305 * ncomps + c);
        const auto *lf_306 = buffer.data(lf + 306 * ncomps + c);
        const auto *lf_307 = buffer.data(lf + 307 * ncomps + c);
        const auto *lf_308 = buffer.data(lf + 308 * ncomps + c);
        const auto *lf_309 = buffer.data(lf + 309 * ncomps + c);
        const auto *lf_310 = buffer.data(lf + 310 * ncomps + c);
        const auto *lf_311 = buffer.data(lf + 311 * ncomps + c);
        const auto *lf_312 = buffer.data(lf + 312 * ncomps + c);
        const auto *lf_313 = buffer.data(lf + 313 * ncomps + c);
        const auto *lf_314 = buffer.data(lf + 314 * ncomps + c);
        const auto *lf_315 = buffer.data(lf + 315 * ncomps + c);
        const auto *lf_316 = buffer.data(lf + 316 * ncomps + c);
        const auto *lf_317 = buffer.data(lf + 317 * ncomps + c);
        const auto *lf_318 = buffer.data(lf + 318 * ncomps + c);
        const auto *lf_319 = buffer.data(lf + 319 * ncomps + c);
        const auto *lf_320 = buffer.data(lf + 320 * ncomps + c);
        const auto *lf_321 = buffer.data(lf + 321 * ncomps + c);
        const auto *lf_322 = buffer.data(lf + 322 * ncomps + c);
        const auto *lf_323 = buffer.data(lf + 323 * ncomps + c);
        const auto *lf_324 = buffer.data(lf + 324 * ncomps + c);
        const auto *lf_325 = buffer.data(lf + 325 * ncomps + c);
        const auto *lf_326 = buffer.data(lf + 326 * ncomps + c);
        const auto *lf_327 = buffer.data(lf + 327 * ncomps + c);
        const auto *lf_328 = buffer.data(lf + 328 * ncomps + c);
        const auto *lf_329 = buffer.data(lf + 329 * ncomps + c);
        const auto *lf_330 = buffer.data(lf + 330 * ncomps + c);
        const auto *lf_331 = buffer.data(lf + 331 * ncomps + c);
        const auto *lf_332 = buffer.data(lf + 332 * ncomps + c);
        const auto *lf_333 = buffer.data(lf + 333 * ncomps + c);
        const auto *lf_334 = buffer.data(lf + 334 * ncomps + c);
        const auto *lf_335 = buffer.data(lf + 335 * ncomps + c);
        const auto *lf_336 = buffer.data(lf + 336 * ncomps + c);
        const auto *lf_337 = buffer.data(lf + 337 * ncomps + c);
        const auto *lf_338 = buffer.data(lf + 338 * ncomps + c);
        const auto *lf_339 = buffer.data(lf + 339 * ncomps + c);
        const auto *lf_340 = buffer.data(lf + 340 * ncomps + c);
        const auto *lf_341 = buffer.data(lf + 341 * ncomps + c);
        const auto *lf_342 = buffer.data(lf + 342 * ncomps + c);
        const auto *lf_343 = buffer.data(lf + 343 * ncomps + c);
        const auto *lf_344 = buffer.data(lf + 344 * ncomps + c);
        const auto *lf_345 = buffer.data(lf + 345 * ncomps + c);
        const auto *lf_346 = buffer.data(lf + 346 * ncomps + c);
        const auto *lf_347 = buffer.data(lf + 347 * ncomps + c);
        const auto *lf_348 = buffer.data(lf + 348 * ncomps + c);
        const auto *lf_349 = buffer.data(lf + 349 * ncomps + c);
        const auto *lf_350 = buffer.data(lf + 350 * ncomps + c);
        const auto *lf_351 = buffer.data(lf + 351 * ncomps + c);
        const auto *lf_352 = buffer.data(lf + 352 * ncomps + c);
        const auto *lf_353 = buffer.data(lf + 353 * ncomps + c);
        const auto *lf_354 = buffer.data(lf + 354 * ncomps + c);
        const auto *lf_355 = buffer.data(lf + 355 * ncomps + c);
        const auto *lf_356 = buffer.data(lf + 356 * ncomps + c);
        const auto *lf_357 = buffer.data(lf + 357 * ncomps + c);
        const auto *lf_358 = buffer.data(lf + 358 * ncomps + c);
        const auto *lf_359 = buffer.data(lf + 359 * ncomps + c);
        const auto *lf_376 = buffer.data(lf + 376 * ncomps + c);
        const auto *lf_377 = buffer.data(lf + 377 * ncomps + c);
        const auto *lf_378 = buffer.data(lf + 378 * ncomps + c);
        const auto *lf_379 = buffer.data(lf + 379 * ncomps + c);
        const auto *lf_386 = buffer.data(lf + 386 * ncomps + c);
        const auto *lf_387 = buffer.data(lf + 387 * ncomps + c);
        const auto *lf_388 = buffer.data(lf + 388 * ncomps + c);
        const auto *lf_389 = buffer.data(lf + 389 * ncomps + c);
        const auto *lf_396 = buffer.data(lf + 396 * ncomps + c);
        const auto *lf_397 = buffer.data(lf + 397 * ncomps + c);
        const auto *lf_398 = buffer.data(lf + 398 * ncomps + c);
        const auto *lf_399 = buffer.data(lf + 399 * ncomps + c);
        const auto *lf_406 = buffer.data(lf + 406 * ncomps + c);
        const auto *lf_407 = buffer.data(lf + 407 * ncomps + c);
        const auto *lf_408 = buffer.data(lf + 408 * ncomps + c);
        const auto *lf_409 = buffer.data(lf + 409 * ncomps + c);
        const auto *lf_416 = buffer.data(lf + 416 * ncomps + c);
        const auto *lf_417 = buffer.data(lf + 417 * ncomps + c);
        const auto *lf_418 = buffer.data(lf + 418 * ncomps + c);
        const auto *lf_419 = buffer.data(lf + 419 * ncomps + c);
        const auto *lf_426 = buffer.data(lf + 426 * ncomps + c);
        const auto *lf_427 = buffer.data(lf + 427 * ncomps + c);
        const auto *lf_428 = buffer.data(lf + 428 * ncomps + c);
        const auto *lf_429 = buffer.data(lf + 429 * ncomps + c);
        const auto *lf_436 = buffer.data(lf + 436 * ncomps + c);
        const auto *lf_437 = buffer.data(lf + 437 * ncomps + c);
        const auto *lf_438 = buffer.data(lf + 438 * ncomps + c);
        const auto *lf_439 = buffer.data(lf + 439 * ncomps + c);
        const auto *lf_449 = buffer.data(lf + 449 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, kf_290, kf_291, kf_292, \
                         kf_293, kf_294, lf_290, lf_291, lf_292, lf_293, \
                         lf_294 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * kf_290[k]
                       + lf_290[k];

            t_436[k] = ab_x[k] * kf_291[k]
                       + lf_291[k];

            t_437[k] = ab_x[k] * kf_292[k]
                       + lf_292[k];

            t_438[k] = ab_x[k] * kf_293[k]
                       + lf_293[k];

            t_439[k] = ab_x[k] * kf_294[k]
                       + lf_294[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, kf_295, kf_296, kf_297, \
                         kf_298, kf_299, lf_295, lf_296, lf_297, lf_298, \
                         lf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * kf_295[k]
                       + lf_295[k];

            t_441[k] = ab_x[k] * kf_296[k]
                       + lf_296[k];

            t_442[k] = ab_x[k] * kf_297[k]
                       + lf_297[k];

            t_443[k] = ab_x[k] * kf_298[k]
                       + lf_298[k];

            t_444[k] = ab_x[k] * kf_299[k]
                       + lf_299[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_y, ab_z, kf_296, kf_297, \
                         kf_298, kf_299, lf_376, lf_377, lf_378, lf_379, \
                         lf_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = ab_y[k] * kf_296[k]
                       + lf_376[k];

            t_446[k] = ab_y[k] * kf_297[k]
                       + lf_377[k];

            t_447[k] = ab_y[k] * kf_298[k]
                       + lf_378[k];

            t_448[k] = ab_y[k] * kf_299[k]
                       + lf_379[k];

            t_449[k] = ab_z[k] * kf_299[k]
                       + lf_389[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_x, kf_300, kf_301, kf_302, \
                         kf_303, kf_304, lf_300, lf_301, lf_302, lf_303, \
                         lf_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = ab_x[k] * kf_300[k]
                       + lf_300[k];

            t_451[k] = ab_x[k] * kf_301[k]
                       + lf_301[k];

            t_452[k] = ab_x[k] * kf_302[k]
                       + lf_302[k];

            t_453[k] = ab_x[k] * kf_303[k]
                       + lf_303[k];

            t_454[k] = ab_x[k] * kf_304[k]
                       + lf_304[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_x, kf_305, kf_306, kf_307, \
                         kf_308, kf_309, lf_305, lf_306, lf_307, lf_308, \
                         lf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = ab_x[k] * kf_305[k]
                       + lf_305[k];

            t_456[k] = ab_x[k] * kf_306[k]
                       + lf_306[k];

            t_457[k] = ab_x[k] * kf_307[k]
                       + lf_307[k];

            t_458[k] = ab_x[k] * kf_308[k]
                       + lf_308[k];

            t_459[k] = ab_x[k] * kf_309[k]
                       + lf_309[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, ab_z, kf_306, kf_307, \
                         kf_308, kf_309, lf_386, lf_387, lf_388, lf_389, \
                         lf_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = ab_y[k] * kf_306[k]
                       + lf_386[k];

            t_461[k] = ab_y[k] * kf_307[k]
                       + lf_387[k];

            t_462[k] = ab_y[k] * kf_308[k]
                       + lf_388[k];

            t_463[k] = ab_y[k] * kf_309[k]
                       + lf_389[k];

            t_464[k] = ab_z[k] * kf_309[k]
                       + lf_399[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_x, kf_310, kf_311, kf_312, \
                         kf_313, kf_314, lf_310, lf_311, lf_312, lf_313, \
                         lf_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = ab_x[k] * kf_310[k]
                       + lf_310[k];

            t_466[k] = ab_x[k] * kf_311[k]
                       + lf_311[k];

            t_467[k] = ab_x[k] * kf_312[k]
                       + lf_312[k];

            t_468[k] = ab_x[k] * kf_313[k]
                       + lf_313[k];

            t_469[k] = ab_x[k] * kf_314[k]
                       + lf_314[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_x, kf_315, kf_316, kf_317, \
                         kf_318, kf_319, lf_315, lf_316, lf_317, lf_318, \
                         lf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = ab_x[k] * kf_315[k]
                       + lf_315[k];

            t_471[k] = ab_x[k] * kf_316[k]
                       + lf_316[k];

            t_472[k] = ab_x[k] * kf_317[k]
                       + lf_317[k];

            t_473[k] = ab_x[k] * kf_318[k]
                       + lf_318[k];

            t_474[k] = ab_x[k] * kf_319[k]
                       + lf_319[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, ab_z, kf_316, kf_317, \
                         kf_318, kf_319, lf_396, lf_397, lf_398, lf_399, \
                         lf_409 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = ab_y[k] * kf_316[k]
                       + lf_396[k];

            t_476[k] = ab_y[k] * kf_317[k]
                       + lf_397[k];

            t_477[k] = ab_y[k] * kf_318[k]
                       + lf_398[k];

            t_478[k] = ab_y[k] * kf_319[k]
                       + lf_399[k];

            t_479[k] = ab_z[k] * kf_319[k]
                       + lf_409[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_x, kf_320, kf_321, kf_322, \
                         kf_323, kf_324, lf_320, lf_321, lf_322, lf_323, \
                         lf_324 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = ab_x[k] * kf_320[k]
                       + lf_320[k];

            t_481[k] = ab_x[k] * kf_321[k]
                       + lf_321[k];

            t_482[k] = ab_x[k] * kf_322[k]
                       + lf_322[k];

            t_483[k] = ab_x[k] * kf_323[k]
                       + lf_323[k];

            t_484[k] = ab_x[k] * kf_324[k]
                       + lf_324[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_x, kf_325, kf_326, kf_327, \
                         kf_328, kf_329, lf_325, lf_326, lf_327, lf_328, \
                         lf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = ab_x[k] * kf_325[k]
                       + lf_325[k];

            t_486[k] = ab_x[k] * kf_326[k]
                       + lf_326[k];

            t_487[k] = ab_x[k] * kf_327[k]
                       + lf_327[k];

            t_488[k] = ab_x[k] * kf_328[k]
                       + lf_328[k];

            t_489[k] = ab_x[k] * kf_329[k]
                       + lf_329[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, ab_z, kf_326, kf_327, \
                         kf_328, kf_329, lf_406, lf_407, lf_408, lf_409, \
                         lf_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = ab_y[k] * kf_326[k]
                       + lf_406[k];

            t_491[k] = ab_y[k] * kf_327[k]
                       + lf_407[k];

            t_492[k] = ab_y[k] * kf_328[k]
                       + lf_408[k];

            t_493[k] = ab_y[k] * kf_329[k]
                       + lf_409[k];

            t_494[k] = ab_z[k] * kf_329[k]
                       + lf_419[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_x, kf_330, kf_331, kf_332, \
                         kf_333, kf_334, lf_330, lf_331, lf_332, lf_333, \
                         lf_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = ab_x[k] * kf_330[k]
                       + lf_330[k];

            t_496[k] = ab_x[k] * kf_331[k]
                       + lf_331[k];

            t_497[k] = ab_x[k] * kf_332[k]
                       + lf_332[k];

            t_498[k] = ab_x[k] * kf_333[k]
                       + lf_333[k];

            t_499[k] = ab_x[k] * kf_334[k]
                       + lf_334[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_x, kf_335, kf_336, kf_337, \
                         kf_338, kf_339, lf_335, lf_336, lf_337, lf_338, \
                         lf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_x[k] * kf_335[k]
                       + lf_335[k];

            t_501[k] = ab_x[k] * kf_336[k]
                       + lf_336[k];

            t_502[k] = ab_x[k] * kf_337[k]
                       + lf_337[k];

            t_503[k] = ab_x[k] * kf_338[k]
                       + lf_338[k];

            t_504[k] = ab_x[k] * kf_339[k]
                       + lf_339[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_y, ab_z, kf_336, kf_337, \
                         kf_338, kf_339, lf_416, lf_417, lf_418, lf_419, \
                         lf_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = ab_y[k] * kf_336[k]
                       + lf_416[k];

            t_506[k] = ab_y[k] * kf_337[k]
                       + lf_417[k];

            t_507[k] = ab_y[k] * kf_338[k]
                       + lf_418[k];

            t_508[k] = ab_y[k] * kf_339[k]
                       + lf_419[k];

            t_509[k] = ab_z[k] * kf_339[k]
                       + lf_429[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_x, kf_340, kf_341, kf_342, \
                         kf_343, kf_344, lf_340, lf_341, lf_342, lf_343, \
                         lf_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = ab_x[k] * kf_340[k]
                       + lf_340[k];

            t_511[k] = ab_x[k] * kf_341[k]
                       + lf_341[k];

            t_512[k] = ab_x[k] * kf_342[k]
                       + lf_342[k];

            t_513[k] = ab_x[k] * kf_343[k]
                       + lf_343[k];

            t_514[k] = ab_x[k] * kf_344[k]
                       + lf_344[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_x, kf_345, kf_346, kf_347, \
                         kf_348, kf_349, lf_345, lf_346, lf_347, lf_348, \
                         lf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = ab_x[k] * kf_345[k]
                       + lf_345[k];

            t_516[k] = ab_x[k] * kf_346[k]
                       + lf_346[k];

            t_517[k] = ab_x[k] * kf_347[k]
                       + lf_347[k];

            t_518[k] = ab_x[k] * kf_348[k]
                       + lf_348[k];

            t_519[k] = ab_x[k] * kf_349[k]
                       + lf_349[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_y, ab_z, kf_346, kf_347, \
                         kf_348, kf_349, lf_426, lf_427, lf_428, lf_429, \
                         lf_439 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = ab_y[k] * kf_346[k]
                       + lf_426[k];

            t_521[k] = ab_y[k] * kf_347[k]
                       + lf_427[k];

            t_522[k] = ab_y[k] * kf_348[k]
                       + lf_428[k];

            t_523[k] = ab_y[k] * kf_349[k]
                       + lf_429[k];

            t_524[k] = ab_z[k] * kf_349[k]
                       + lf_439[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_x, kf_350, kf_351, kf_352, \
                         kf_353, kf_354, lf_350, lf_351, lf_352, lf_353, \
                         lf_354 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = ab_x[k] * kf_350[k]
                       + lf_350[k];

            t_526[k] = ab_x[k] * kf_351[k]
                       + lf_351[k];

            t_527[k] = ab_x[k] * kf_352[k]
                       + lf_352[k];

            t_528[k] = ab_x[k] * kf_353[k]
                       + lf_353[k];

            t_529[k] = ab_x[k] * kf_354[k]
                       + lf_354[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_x, kf_355, kf_356, kf_357, \
                         kf_358, kf_359, lf_355, lf_356, lf_357, lf_358, \
                         lf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = ab_x[k] * kf_355[k]
                       + lf_355[k];

            t_531[k] = ab_x[k] * kf_356[k]
                       + lf_356[k];

            t_532[k] = ab_x[k] * kf_357[k]
                       + lf_357[k];

            t_533[k] = ab_x[k] * kf_358[k]
                       + lf_358[k];

            t_534[k] = ab_x[k] * kf_359[k]
                       + lf_359[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_y, ab_z, kf_356, kf_357, \
                         kf_358, kf_359, lf_436, lf_437, lf_438, lf_439, \
                         lf_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = ab_y[k] * kf_356[k]
                       + lf_436[k];

            t_536[k] = ab_y[k] * kf_357[k]
                       + lf_437[k];

            t_537[k] = ab_y[k] * kf_358[k]
                       + lf_438[k];

            t_538[k] = ab_y[k] * kf_359[k]
                       + lf_439[k];

            t_539[k] = ab_z[k] * kf_359[k]
                       + lf_449[k];
        }
    }
}

auto
compute_hrr_kg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t kf, const size_t lf, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_kg_piece0(buffer, coordinates, target, kf, lf, ncomps, nmax);

    compute_hrr_kg_piece1(buffer, coordinates, target, kf, lf, ncomps, nmax);

    compute_hrr_kg_piece2(buffer, coordinates, target, kf, lf, ncomps, nmax);

    compute_hrr_kg_piece3(buffer, coordinates, target, kf, lf, ncomps, nmax);
}

}  // namespace simdtrf
