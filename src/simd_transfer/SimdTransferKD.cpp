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


#include "SimdTransferKD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_kd_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kp, const size_t lp,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kp_0 = buffer.data(kp + 0 * ncomps + c);
        const auto *kp_1 = buffer.data(kp + 1 * ncomps + c);
        const auto *kp_2 = buffer.data(kp + 2 * ncomps + c);
        const auto *kp_3 = buffer.data(kp + 3 * ncomps + c);
        const auto *kp_4 = buffer.data(kp + 4 * ncomps + c);
        const auto *kp_5 = buffer.data(kp + 5 * ncomps + c);
        const auto *kp_6 = buffer.data(kp + 6 * ncomps + c);
        const auto *kp_7 = buffer.data(kp + 7 * ncomps + c);
        const auto *kp_8 = buffer.data(kp + 8 * ncomps + c);
        const auto *kp_9 = buffer.data(kp + 9 * ncomps + c);
        const auto *kp_10 = buffer.data(kp + 10 * ncomps + c);
        const auto *kp_11 = buffer.data(kp + 11 * ncomps + c);
        const auto *kp_12 = buffer.data(kp + 12 * ncomps + c);
        const auto *kp_13 = buffer.data(kp + 13 * ncomps + c);
        const auto *kp_14 = buffer.data(kp + 14 * ncomps + c);
        const auto *kp_15 = buffer.data(kp + 15 * ncomps + c);
        const auto *kp_16 = buffer.data(kp + 16 * ncomps + c);
        const auto *kp_17 = buffer.data(kp + 17 * ncomps + c);
        const auto *kp_18 = buffer.data(kp + 18 * ncomps + c);
        const auto *kp_19 = buffer.data(kp + 19 * ncomps + c);
        const auto *kp_20 = buffer.data(kp + 20 * ncomps + c);
        const auto *kp_21 = buffer.data(kp + 21 * ncomps + c);
        const auto *kp_22 = buffer.data(kp + 22 * ncomps + c);
        const auto *kp_23 = buffer.data(kp + 23 * ncomps + c);
        const auto *kp_24 = buffer.data(kp + 24 * ncomps + c);
        const auto *kp_25 = buffer.data(kp + 25 * ncomps + c);
        const auto *kp_26 = buffer.data(kp + 26 * ncomps + c);
        const auto *kp_27 = buffer.data(kp + 27 * ncomps + c);
        const auto *kp_28 = buffer.data(kp + 28 * ncomps + c);
        const auto *kp_29 = buffer.data(kp + 29 * ncomps + c);
        const auto *kp_30 = buffer.data(kp + 30 * ncomps + c);
        const auto *kp_31 = buffer.data(kp + 31 * ncomps + c);
        const auto *kp_32 = buffer.data(kp + 32 * ncomps + c);
        const auto *kp_33 = buffer.data(kp + 33 * ncomps + c);
        const auto *kp_34 = buffer.data(kp + 34 * ncomps + c);
        const auto *kp_35 = buffer.data(kp + 35 * ncomps + c);
        const auto *kp_36 = buffer.data(kp + 36 * ncomps + c);
        const auto *kp_37 = buffer.data(kp + 37 * ncomps + c);
        const auto *kp_38 = buffer.data(kp + 38 * ncomps + c);
        const auto *kp_39 = buffer.data(kp + 39 * ncomps + c);
        const auto *kp_40 = buffer.data(kp + 40 * ncomps + c);
        const auto *kp_41 = buffer.data(kp + 41 * ncomps + c);
        const auto *kp_42 = buffer.data(kp + 42 * ncomps + c);
        const auto *kp_43 = buffer.data(kp + 43 * ncomps + c);
        const auto *kp_44 = buffer.data(kp + 44 * ncomps + c);
        const auto *kp_45 = buffer.data(kp + 45 * ncomps + c);
        const auto *kp_46 = buffer.data(kp + 46 * ncomps + c);
        const auto *kp_47 = buffer.data(kp + 47 * ncomps + c);
        const auto *kp_48 = buffer.data(kp + 48 * ncomps + c);
        const auto *kp_49 = buffer.data(kp + 49 * ncomps + c);
        const auto *kp_50 = buffer.data(kp + 50 * ncomps + c);
        const auto *kp_51 = buffer.data(kp + 51 * ncomps + c);
        const auto *kp_52 = buffer.data(kp + 52 * ncomps + c);
        const auto *kp_53 = buffer.data(kp + 53 * ncomps + c);
        const auto *kp_54 = buffer.data(kp + 54 * ncomps + c);
        const auto *kp_55 = buffer.data(kp + 55 * ncomps + c);
        const auto *kp_56 = buffer.data(kp + 56 * ncomps + c);
        const auto *kp_57 = buffer.data(kp + 57 * ncomps + c);
        const auto *kp_58 = buffer.data(kp + 58 * ncomps + c);
        const auto *kp_59 = buffer.data(kp + 59 * ncomps + c);
        const auto *kp_60 = buffer.data(kp + 60 * ncomps + c);
        const auto *kp_61 = buffer.data(kp + 61 * ncomps + c);
        const auto *kp_62 = buffer.data(kp + 62 * ncomps + c);
        const auto *kp_63 = buffer.data(kp + 63 * ncomps + c);
        const auto *kp_64 = buffer.data(kp + 64 * ncomps + c);
        const auto *kp_65 = buffer.data(kp + 65 * ncomps + c);
        const auto *kp_66 = buffer.data(kp + 66 * ncomps + c);
        const auto *kp_67 = buffer.data(kp + 67 * ncomps + c);
        const auto *kp_68 = buffer.data(kp + 68 * ncomps + c);
        const auto *kp_69 = buffer.data(kp + 69 * ncomps + c);
        const auto *kp_70 = buffer.data(kp + 70 * ncomps + c);
        const auto *kp_71 = buffer.data(kp + 71 * ncomps + c);

        const auto *lp_0 = buffer.data(lp + 0 * ncomps + c);
        const auto *lp_1 = buffer.data(lp + 1 * ncomps + c);
        const auto *lp_2 = buffer.data(lp + 2 * ncomps + c);
        const auto *lp_3 = buffer.data(lp + 3 * ncomps + c);
        const auto *lp_4 = buffer.data(lp + 4 * ncomps + c);
        const auto *lp_5 = buffer.data(lp + 5 * ncomps + c);
        const auto *lp_6 = buffer.data(lp + 6 * ncomps + c);
        const auto *lp_7 = buffer.data(lp + 7 * ncomps + c);
        const auto *lp_8 = buffer.data(lp + 8 * ncomps + c);
        const auto *lp_9 = buffer.data(lp + 9 * ncomps + c);
        const auto *lp_10 = buffer.data(lp + 10 * ncomps + c);
        const auto *lp_11 = buffer.data(lp + 11 * ncomps + c);
        const auto *lp_12 = buffer.data(lp + 12 * ncomps + c);
        const auto *lp_13 = buffer.data(lp + 13 * ncomps + c);
        const auto *lp_14 = buffer.data(lp + 14 * ncomps + c);
        const auto *lp_15 = buffer.data(lp + 15 * ncomps + c);
        const auto *lp_16 = buffer.data(lp + 16 * ncomps + c);
        const auto *lp_17 = buffer.data(lp + 17 * ncomps + c);
        const auto *lp_18 = buffer.data(lp + 18 * ncomps + c);
        const auto *lp_19 = buffer.data(lp + 19 * ncomps + c);
        const auto *lp_20 = buffer.data(lp + 20 * ncomps + c);
        const auto *lp_21 = buffer.data(lp + 21 * ncomps + c);
        const auto *lp_22 = buffer.data(lp + 22 * ncomps + c);
        const auto *lp_23 = buffer.data(lp + 23 * ncomps + c);
        const auto *lp_24 = buffer.data(lp + 24 * ncomps + c);
        const auto *lp_25 = buffer.data(lp + 25 * ncomps + c);
        const auto *lp_26 = buffer.data(lp + 26 * ncomps + c);
        const auto *lp_27 = buffer.data(lp + 27 * ncomps + c);
        const auto *lp_28 = buffer.data(lp + 28 * ncomps + c);
        const auto *lp_29 = buffer.data(lp + 29 * ncomps + c);
        const auto *lp_30 = buffer.data(lp + 30 * ncomps + c);
        const auto *lp_31 = buffer.data(lp + 31 * ncomps + c);
        const auto *lp_32 = buffer.data(lp + 32 * ncomps + c);
        const auto *lp_33 = buffer.data(lp + 33 * ncomps + c);
        const auto *lp_34 = buffer.data(lp + 34 * ncomps + c);
        const auto *lp_35 = buffer.data(lp + 35 * ncomps + c);
        const auto *lp_36 = buffer.data(lp + 36 * ncomps + c);
        const auto *lp_37 = buffer.data(lp + 37 * ncomps + c);
        const auto *lp_38 = buffer.data(lp + 38 * ncomps + c);
        const auto *lp_39 = buffer.data(lp + 39 * ncomps + c);
        const auto *lp_40 = buffer.data(lp + 40 * ncomps + c);
        const auto *lp_41 = buffer.data(lp + 41 * ncomps + c);
        const auto *lp_42 = buffer.data(lp + 42 * ncomps + c);
        const auto *lp_43 = buffer.data(lp + 43 * ncomps + c);
        const auto *lp_44 = buffer.data(lp + 44 * ncomps + c);
        const auto *lp_45 = buffer.data(lp + 45 * ncomps + c);
        const auto *lp_46 = buffer.data(lp + 46 * ncomps + c);
        const auto *lp_47 = buffer.data(lp + 47 * ncomps + c);
        const auto *lp_48 = buffer.data(lp + 48 * ncomps + c);
        const auto *lp_49 = buffer.data(lp + 49 * ncomps + c);
        const auto *lp_50 = buffer.data(lp + 50 * ncomps + c);
        const auto *lp_51 = buffer.data(lp + 51 * ncomps + c);
        const auto *lp_52 = buffer.data(lp + 52 * ncomps + c);
        const auto *lp_53 = buffer.data(lp + 53 * ncomps + c);
        const auto *lp_54 = buffer.data(lp + 54 * ncomps + c);
        const auto *lp_55 = buffer.data(lp + 55 * ncomps + c);
        const auto *lp_56 = buffer.data(lp + 56 * ncomps + c);
        const auto *lp_57 = buffer.data(lp + 57 * ncomps + c);
        const auto *lp_58 = buffer.data(lp + 58 * ncomps + c);
        const auto *lp_59 = buffer.data(lp + 59 * ncomps + c);
        const auto *lp_60 = buffer.data(lp + 60 * ncomps + c);
        const auto *lp_61 = buffer.data(lp + 61 * ncomps + c);
        const auto *lp_62 = buffer.data(lp + 62 * ncomps + c);
        const auto *lp_63 = buffer.data(lp + 63 * ncomps + c);
        const auto *lp_64 = buffer.data(lp + 64 * ncomps + c);
        const auto *lp_65 = buffer.data(lp + 65 * ncomps + c);
        const auto *lp_66 = buffer.data(lp + 66 * ncomps + c);
        const auto *lp_67 = buffer.data(lp + 67 * ncomps + c);
        const auto *lp_68 = buffer.data(lp + 68 * ncomps + c);
        const auto *lp_69 = buffer.data(lp + 69 * ncomps + c);
        const auto *lp_70 = buffer.data(lp + 70 * ncomps + c);
        const auto *lp_71 = buffer.data(lp + 71 * ncomps + c);
        const auto *lp_73 = buffer.data(lp + 73 * ncomps + c);
        const auto *lp_74 = buffer.data(lp + 74 * ncomps + c);
        const auto *lp_76 = buffer.data(lp + 76 * ncomps + c);
        const auto *lp_77 = buffer.data(lp + 77 * ncomps + c);
        const auto *lp_79 = buffer.data(lp + 79 * ncomps + c);
        const auto *lp_80 = buffer.data(lp + 80 * ncomps + c);
        const auto *lp_83 = buffer.data(lp + 83 * ncomps + c);
        const auto *lp_85 = buffer.data(lp + 85 * ncomps + c);
        const auto *lp_86 = buffer.data(lp + 86 * ncomps + c);
        const auto *lp_88 = buffer.data(lp + 88 * ncomps + c);
        const auto *lp_89 = buffer.data(lp + 89 * ncomps + c);
        const auto *lp_91 = buffer.data(lp + 91 * ncomps + c);
        const auto *lp_92 = buffer.data(lp + 92 * ncomps + c);
        const auto *lp_95 = buffer.data(lp + 95 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, kp_0, kp_1, kp_2, lp_0, lp_1, \
                         lp_2, lp_4, lp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * kp_0[k]
                     + lp_0[k];

            t_1[k] = ab_x[k] * kp_1[k]
                     + lp_1[k];

            t_2[k] = ab_x[k] * kp_2[k]
                     + lp_2[k];

            t_3[k] = ab_y[k] * kp_1[k]
                     + lp_4[k];

            t_4[k] = ab_y[k] * kp_2[k]
                     + lp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, kp_2, kp_3, kp_4, kp_5, lp_3, lp_4, \
                         lp_5, lp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * kp_2[k]
                     + lp_8[k];

            t_6[k] = ab_x[k] * kp_3[k]
                     + lp_3[k];

            t_7[k] = ab_x[k] * kp_4[k]
                     + lp_4[k];

            t_8[k] = ab_x[k] * kp_5[k]
                     + lp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, kp_4, kp_5, kp_6, lp_6, \
                         lp_10, lp_11, lp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * kp_4[k]
                     + lp_10[k];

            t_10[k] = ab_y[k] * kp_5[k]
                      + lp_11[k];

            t_11[k] = ab_z[k] * kp_5[k]
                      + lp_14[k];

            t_12[k] = ab_x[k] * kp_6[k]
                      + lp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, kp_7, kp_8, lp_7, \
                         lp_8, lp_13, lp_14, lp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * kp_7[k]
                      + lp_7[k];

            t_14[k] = ab_x[k] * kp_8[k]
                      + lp_8[k];

            t_15[k] = ab_y[k] * kp_7[k]
                      + lp_13[k];

            t_16[k] = ab_y[k] * kp_8[k]
                      + lp_14[k];

            t_17[k] = ab_z[k] * kp_8[k]
                      + lp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, kp_9, kp_10, kp_11, lp_9, \
                         lp_10, lp_11, lp_19, lp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * kp_9[k]
                      + lp_9[k];

            t_19[k] = ab_x[k] * kp_10[k]
                      + lp_10[k];

            t_20[k] = ab_x[k] * kp_11[k]
                      + lp_11[k];

            t_21[k] = ab_y[k] * kp_10[k]
                      + lp_19[k];

            t_22[k] = ab_y[k] * kp_11[k]
                      + lp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, kp_11, kp_12, kp_13, kp_14, \
                         lp_12, lp_13, lp_14, lp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * kp_11[k]
                      + lp_23[k];

            t_24[k] = ab_x[k] * kp_12[k]
                      + lp_12[k];

            t_25[k] = ab_x[k] * kp_13[k]
                      + lp_13[k];

            t_26[k] = ab_x[k] * kp_14[k]
                      + lp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, kp_13, kp_14, kp_15, lp_15, \
                         lp_22, lp_23, lp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * kp_13[k]
                      + lp_22[k];

            t_28[k] = ab_y[k] * kp_14[k]
                      + lp_23[k];

            t_29[k] = ab_z[k] * kp_14[k]
                      + lp_26[k];

            t_30[k] = ab_x[k] * kp_15[k]
                      + lp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, kp_16, kp_17, lp_16, \
                         lp_17, lp_25, lp_26, lp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * kp_16[k]
                      + lp_16[k];

            t_32[k] = ab_x[k] * kp_17[k]
                      + lp_17[k];

            t_33[k] = ab_y[k] * kp_16[k]
                      + lp_25[k];

            t_34[k] = ab_y[k] * kp_17[k]
                      + lp_26[k];

            t_35[k] = ab_z[k] * kp_17[k]
                      + lp_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, kp_18, kp_19, kp_20, lp_18, \
                         lp_19, lp_20, lp_31, lp_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * kp_18[k]
                      + lp_18[k];

            t_37[k] = ab_x[k] * kp_19[k]
                      + lp_19[k];

            t_38[k] = ab_x[k] * kp_20[k]
                      + lp_20[k];

            t_39[k] = ab_y[k] * kp_19[k]
                      + lp_31[k];

            t_40[k] = ab_y[k] * kp_20[k]
                      + lp_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, kp_20, kp_21, kp_22, kp_23, \
                         lp_21, lp_22, lp_23, lp_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * kp_20[k]
                      + lp_35[k];

            t_42[k] = ab_x[k] * kp_21[k]
                      + lp_21[k];

            t_43[k] = ab_x[k] * kp_22[k]
                      + lp_22[k];

            t_44[k] = ab_x[k] * kp_23[k]
                      + lp_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, kp_22, kp_23, kp_24, lp_24, \
                         lp_34, lp_35, lp_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * kp_22[k]
                      + lp_34[k];

            t_46[k] = ab_y[k] * kp_23[k]
                      + lp_35[k];

            t_47[k] = ab_z[k] * kp_23[k]
                      + lp_38[k];

            t_48[k] = ab_x[k] * kp_24[k]
                      + lp_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, kp_25, kp_26, lp_25, \
                         lp_26, lp_37, lp_38, lp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * kp_25[k]
                      + lp_25[k];

            t_50[k] = ab_x[k] * kp_26[k]
                      + lp_26[k];

            t_51[k] = ab_y[k] * kp_25[k]
                      + lp_37[k];

            t_52[k] = ab_y[k] * kp_26[k]
                      + lp_38[k];

            t_53[k] = ab_z[k] * kp_26[k]
                      + lp_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, kp_27, kp_28, kp_29, lp_27, \
                         lp_28, lp_29, lp_40, lp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * kp_27[k]
                      + lp_27[k];

            t_55[k] = ab_x[k] * kp_28[k]
                      + lp_28[k];

            t_56[k] = ab_x[k] * kp_29[k]
                      + lp_29[k];

            t_57[k] = ab_y[k] * kp_28[k]
                      + lp_40[k];

            t_58[k] = ab_y[k] * kp_29[k]
                      + lp_41[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, kp_29, kp_30, kp_31, kp_32, \
                         lp_30, lp_31, lp_32, lp_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * kp_29[k]
                      + lp_44[k];

            t_60[k] = ab_x[k] * kp_30[k]
                      + lp_30[k];

            t_61[k] = ab_x[k] * kp_31[k]
                      + lp_31[k];

            t_62[k] = ab_x[k] * kp_32[k]
                      + lp_32[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, kp_31, kp_32, kp_33, lp_33, \
                         lp_46, lp_47, lp_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = ab_y[k] * kp_31[k]
                      + lp_46[k];

            t_64[k] = ab_y[k] * kp_32[k]
                      + lp_47[k];

            t_65[k] = ab_z[k] * kp_32[k]
                      + lp_50[k];

            t_66[k] = ab_x[k] * kp_33[k]
                      + lp_33[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, kp_34, kp_35, lp_34, \
                         lp_35, lp_49, lp_50, lp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = ab_x[k] * kp_34[k]
                      + lp_34[k];

            t_68[k] = ab_x[k] * kp_35[k]
                      + lp_35[k];

            t_69[k] = ab_y[k] * kp_34[k]
                      + lp_49[k];

            t_70[k] = ab_y[k] * kp_35[k]
                      + lp_50[k];

            t_71[k] = ab_z[k] * kp_35[k]
                      + lp_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, kp_36, kp_37, kp_38, lp_36, \
                         lp_37, lp_38, lp_52, lp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * kp_36[k]
                      + lp_36[k];

            t_73[k] = ab_x[k] * kp_37[k]
                      + lp_37[k];

            t_74[k] = ab_x[k] * kp_38[k]
                      + lp_38[k];

            t_75[k] = ab_y[k] * kp_37[k]
                      + lp_52[k];

            t_76[k] = ab_y[k] * kp_38[k]
                      + lp_53[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, kp_38, kp_39, kp_40, kp_41, \
                         lp_39, lp_40, lp_41, lp_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * kp_38[k]
                      + lp_56[k];

            t_78[k] = ab_x[k] * kp_39[k]
                      + lp_39[k];

            t_79[k] = ab_x[k] * kp_40[k]
                      + lp_40[k];

            t_80[k] = ab_x[k] * kp_41[k]
                      + lp_41[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, kp_40, kp_41, kp_42, lp_42, \
                         lp_55, lp_56, lp_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = ab_y[k] * kp_40[k]
                      + lp_55[k];

            t_82[k] = ab_y[k] * kp_41[k]
                      + lp_56[k];

            t_83[k] = ab_z[k] * kp_41[k]
                      + lp_59[k];

            t_84[k] = ab_x[k] * kp_42[k]
                      + lp_42[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, kp_43, kp_44, lp_43, \
                         lp_44, lp_58, lp_59, lp_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * kp_43[k]
                      + lp_43[k];

            t_86[k] = ab_x[k] * kp_44[k]
                      + lp_44[k];

            t_87[k] = ab_y[k] * kp_43[k]
                      + lp_58[k];

            t_88[k] = ab_y[k] * kp_44[k]
                      + lp_59[k];

            t_89[k] = ab_z[k] * kp_44[k]
                      + lp_62[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, kp_45, kp_46, kp_47, lp_45, \
                         lp_46, lp_47, lp_64, lp_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * kp_45[k]
                      + lp_45[k];

            t_91[k] = ab_x[k] * kp_46[k]
                      + lp_46[k];

            t_92[k] = ab_x[k] * kp_47[k]
                      + lp_47[k];

            t_93[k] = ab_y[k] * kp_46[k]
                      + lp_64[k];

            t_94[k] = ab_y[k] * kp_47[k]
                      + lp_65[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, kp_47, kp_48, kp_49, kp_50, \
                         lp_48, lp_49, lp_50, lp_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_z[k] * kp_47[k]
                      + lp_68[k];

            t_96[k] = ab_x[k] * kp_48[k]
                      + lp_48[k];

            t_97[k] = ab_x[k] * kp_49[k]
                      + lp_49[k];

            t_98[k] = ab_x[k] * kp_50[k]
                      + lp_50[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, kp_49, kp_50, kp_51, \
                         lp_51, lp_67, lp_68, lp_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * kp_49[k]
                      + lp_67[k];

            t_100[k] = ab_y[k] * kp_50[k]
                       + lp_68[k];

            t_101[k] = ab_z[k] * kp_50[k]
                       + lp_71[k];

            t_102[k] = ab_x[k] * kp_51[k]
                       + lp_51[k];
        }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, kp_52, kp_53, \
                         lp_52, lp_53, lp_70, lp_71, lp_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_103[k] = ab_x[k] * kp_52[k]
                       + lp_52[k];

            t_104[k] = ab_x[k] * kp_53[k]
                       + lp_53[k];

            t_105[k] = ab_y[k] * kp_52[k]
                       + lp_70[k];

            t_106[k] = ab_y[k] * kp_53[k]
                       + lp_71[k];

            t_107[k] = ab_z[k] * kp_53[k]
                       + lp_74[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, kp_54, kp_55, kp_56, \
                         lp_54, lp_55, lp_56, lp_73, lp_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * kp_54[k]
                       + lp_54[k];

            t_109[k] = ab_x[k] * kp_55[k]
                       + lp_55[k];

            t_110[k] = ab_x[k] * kp_56[k]
                       + lp_56[k];

            t_111[k] = ab_y[k] * kp_55[k]
                       + lp_73[k];

            t_112[k] = ab_y[k] * kp_56[k]
                       + lp_74[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, kp_56, kp_57, kp_58, kp_59, \
                         lp_57, lp_58, lp_59, lp_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_z[k] * kp_56[k]
                       + lp_77[k];

            t_114[k] = ab_x[k] * kp_57[k]
                       + lp_57[k];

            t_115[k] = ab_x[k] * kp_58[k]
                       + lp_58[k];

            t_116[k] = ab_x[k] * kp_59[k]
                       + lp_59[k];
        }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, kp_58, kp_59, kp_60, \
                         lp_60, lp_76, lp_77, lp_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_117[k] = ab_y[k] * kp_58[k]
                       + lp_76[k];

            t_118[k] = ab_y[k] * kp_59[k]
                       + lp_77[k];

            t_119[k] = ab_z[k] * kp_59[k]
                       + lp_80[k];

            t_120[k] = ab_x[k] * kp_60[k]
                       + lp_60[k];
        }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, kp_61, kp_62, \
                         lp_61, lp_62, lp_79, lp_80, lp_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_121[k] = ab_x[k] * kp_61[k]
                       + lp_61[k];

            t_122[k] = ab_x[k] * kp_62[k]
                       + lp_62[k];

            t_123[k] = ab_y[k] * kp_61[k]
                       + lp_79[k];

            t_124[k] = ab_y[k] * kp_62[k]
                       + lp_80[k];

            t_125[k] = ab_z[k] * kp_62[k]
                       + lp_83[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ab_y, kp_63, kp_64, kp_65, \
                         lp_63, lp_64, lp_65, lp_85, lp_86 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = ab_x[k] * kp_63[k]
                       + lp_63[k];

            t_127[k] = ab_x[k] * kp_64[k]
                       + lp_64[k];

            t_128[k] = ab_x[k] * kp_65[k]
                       + lp_65[k];

            t_129[k] = ab_y[k] * kp_64[k]
                       + lp_85[k];

            t_130[k] = ab_y[k] * kp_65[k]
                       + lp_86[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, ab_x, ab_z, kp_65, kp_66, kp_67, kp_68, \
                         lp_66, lp_67, lp_68, lp_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = ab_z[k] * kp_65[k]
                       + lp_89[k];

            t_132[k] = ab_x[k] * kp_66[k]
                       + lp_66[k];

            t_133[k] = ab_x[k] * kp_67[k]
                       + lp_67[k];

            t_134[k] = ab_x[k] * kp_68[k]
                       + lp_68[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, kp_67, kp_68, kp_69, \
                         lp_69, lp_88, lp_89, lp_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_y[k] * kp_67[k]
                       + lp_88[k];

            t_136[k] = ab_y[k] * kp_68[k]
                       + lp_89[k];

            t_137[k] = ab_z[k] * kp_68[k]
                       + lp_92[k];

            t_138[k] = ab_x[k] * kp_69[k]
                       + lp_69[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, kp_70, kp_71, \
                         lp_70, lp_71, lp_91, lp_92, lp_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_x[k] * kp_70[k]
                       + lp_70[k];

            t_140[k] = ab_x[k] * kp_71[k]
                       + lp_71[k];

            t_141[k] = ab_y[k] * kp_70[k]
                       + lp_91[k];

            t_142[k] = ab_y[k] * kp_71[k]
                       + lp_92[k];

            t_143[k] = ab_z[k] * kp_71[k]
                       + lp_95[k];
        }
    }
}

static auto
compute_hrr_kd_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kp, const size_t lp,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_144 = buffer.data(target + 144 * ncomps + c);
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kp_72 = buffer.data(kp + 72 * ncomps + c);
        const auto *kp_73 = buffer.data(kp + 73 * ncomps + c);
        const auto *kp_74 = buffer.data(kp + 74 * ncomps + c);
        const auto *kp_75 = buffer.data(kp + 75 * ncomps + c);
        const auto *kp_76 = buffer.data(kp + 76 * ncomps + c);
        const auto *kp_77 = buffer.data(kp + 77 * ncomps + c);
        const auto *kp_78 = buffer.data(kp + 78 * ncomps + c);
        const auto *kp_79 = buffer.data(kp + 79 * ncomps + c);
        const auto *kp_80 = buffer.data(kp + 80 * ncomps + c);
        const auto *kp_81 = buffer.data(kp + 81 * ncomps + c);
        const auto *kp_82 = buffer.data(kp + 82 * ncomps + c);
        const auto *kp_83 = buffer.data(kp + 83 * ncomps + c);
        const auto *kp_84 = buffer.data(kp + 84 * ncomps + c);
        const auto *kp_85 = buffer.data(kp + 85 * ncomps + c);
        const auto *kp_86 = buffer.data(kp + 86 * ncomps + c);
        const auto *kp_87 = buffer.data(kp + 87 * ncomps + c);
        const auto *kp_88 = buffer.data(kp + 88 * ncomps + c);
        const auto *kp_89 = buffer.data(kp + 89 * ncomps + c);
        const auto *kp_90 = buffer.data(kp + 90 * ncomps + c);
        const auto *kp_91 = buffer.data(kp + 91 * ncomps + c);
        const auto *kp_92 = buffer.data(kp + 92 * ncomps + c);
        const auto *kp_93 = buffer.data(kp + 93 * ncomps + c);
        const auto *kp_94 = buffer.data(kp + 94 * ncomps + c);
        const auto *kp_95 = buffer.data(kp + 95 * ncomps + c);
        const auto *kp_96 = buffer.data(kp + 96 * ncomps + c);
        const auto *kp_97 = buffer.data(kp + 97 * ncomps + c);
        const auto *kp_98 = buffer.data(kp + 98 * ncomps + c);
        const auto *kp_99 = buffer.data(kp + 99 * ncomps + c);
        const auto *kp_100 = buffer.data(kp + 100 * ncomps + c);
        const auto *kp_101 = buffer.data(kp + 101 * ncomps + c);
        const auto *kp_102 = buffer.data(kp + 102 * ncomps + c);
        const auto *kp_103 = buffer.data(kp + 103 * ncomps + c);
        const auto *kp_104 = buffer.data(kp + 104 * ncomps + c);
        const auto *kp_105 = buffer.data(kp + 105 * ncomps + c);
        const auto *kp_106 = buffer.data(kp + 106 * ncomps + c);
        const auto *kp_107 = buffer.data(kp + 107 * ncomps + c);

        const auto *lp_72 = buffer.data(lp + 72 * ncomps + c);
        const auto *lp_73 = buffer.data(lp + 73 * ncomps + c);
        const auto *lp_74 = buffer.data(lp + 74 * ncomps + c);
        const auto *lp_75 = buffer.data(lp + 75 * ncomps + c);
        const auto *lp_76 = buffer.data(lp + 76 * ncomps + c);
        const auto *lp_77 = buffer.data(lp + 77 * ncomps + c);
        const auto *lp_78 = buffer.data(lp + 78 * ncomps + c);
        const auto *lp_79 = buffer.data(lp + 79 * ncomps + c);
        const auto *lp_80 = buffer.data(lp + 80 * ncomps + c);
        const auto *lp_81 = buffer.data(lp + 81 * ncomps + c);
        const auto *lp_82 = buffer.data(lp + 82 * ncomps + c);
        const auto *lp_83 = buffer.data(lp + 83 * ncomps + c);
        const auto *lp_84 = buffer.data(lp + 84 * ncomps + c);
        const auto *lp_85 = buffer.data(lp + 85 * ncomps + c);
        const auto *lp_86 = buffer.data(lp + 86 * ncomps + c);
        const auto *lp_87 = buffer.data(lp + 87 * ncomps + c);
        const auto *lp_88 = buffer.data(lp + 88 * ncomps + c);
        const auto *lp_89 = buffer.data(lp + 89 * ncomps + c);
        const auto *lp_90 = buffer.data(lp + 90 * ncomps + c);
        const auto *lp_91 = buffer.data(lp + 91 * ncomps + c);
        const auto *lp_92 = buffer.data(lp + 92 * ncomps + c);
        const auto *lp_93 = buffer.data(lp + 93 * ncomps + c);
        const auto *lp_94 = buffer.data(lp + 94 * ncomps + c);
        const auto *lp_95 = buffer.data(lp + 95 * ncomps + c);
        const auto *lp_96 = buffer.data(lp + 96 * ncomps + c);
        const auto *lp_97 = buffer.data(lp + 97 * ncomps + c);
        const auto *lp_98 = buffer.data(lp + 98 * ncomps + c);
        const auto *lp_99 = buffer.data(lp + 99 * ncomps + c);
        const auto *lp_100 = buffer.data(lp + 100 * ncomps + c);
        const auto *lp_101 = buffer.data(lp + 101 * ncomps + c);
        const auto *lp_102 = buffer.data(lp + 102 * ncomps + c);
        const auto *lp_103 = buffer.data(lp + 103 * ncomps + c);
        const auto *lp_104 = buffer.data(lp + 104 * ncomps + c);
        const auto *lp_105 = buffer.data(lp + 105 * ncomps + c);
        const auto *lp_106 = buffer.data(lp + 106 * ncomps + c);
        const auto *lp_107 = buffer.data(lp + 107 * ncomps + c);
        const auto *lp_109 = buffer.data(lp + 109 * ncomps + c);
        const auto *lp_110 = buffer.data(lp + 110 * ncomps + c);
        const auto *lp_112 = buffer.data(lp + 112 * ncomps + c);
        const auto *lp_113 = buffer.data(lp + 113 * ncomps + c);
        const auto *lp_115 = buffer.data(lp + 115 * ncomps + c);
        const auto *lp_116 = buffer.data(lp + 116 * ncomps + c);
        const auto *lp_118 = buffer.data(lp + 118 * ncomps + c);
        const auto *lp_119 = buffer.data(lp + 119 * ncomps + c);
        const auto *lp_121 = buffer.data(lp + 121 * ncomps + c);
        const auto *lp_122 = buffer.data(lp + 122 * ncomps + c);
        const auto *lp_124 = buffer.data(lp + 124 * ncomps + c);
        const auto *lp_125 = buffer.data(lp + 125 * ncomps + c);
        const auto *lp_127 = buffer.data(lp + 127 * ncomps + c);
        const auto *lp_128 = buffer.data(lp + 128 * ncomps + c);
        const auto *lp_130 = buffer.data(lp + 130 * ncomps + c);
        const auto *lp_131 = buffer.data(lp + 131 * ncomps + c);
        const auto *lp_134 = buffer.data(lp + 134 * ncomps + c);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_x, ab_y, kp_72, kp_73, kp_74, \
                         lp_72, lp_73, lp_74, lp_94, lp_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * kp_72[k]
                       + lp_72[k];

            t_145[k] = ab_x[k] * kp_73[k]
                       + lp_73[k];

            t_146[k] = ab_x[k] * kp_74[k]
                       + lp_74[k];

            t_147[k] = ab_y[k] * kp_73[k]
                       + lp_94[k];

            t_148[k] = ab_y[k] * kp_74[k]
                       + lp_95[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, ab_x, ab_z, kp_74, kp_75, kp_76, kp_77, \
                         lp_75, lp_76, lp_77, lp_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = ab_z[k] * kp_74[k]
                       + lp_98[k];

            t_150[k] = ab_x[k] * kp_75[k]
                       + lp_75[k];

            t_151[k] = ab_x[k] * kp_76[k]
                       + lp_76[k];

            t_152[k] = ab_x[k] * kp_77[k]
                       + lp_77[k];
        }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, ab_x, ab_y, ab_z, kp_76, kp_77, kp_78, \
                         lp_78, lp_97, lp_98, lp_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_153[k] = ab_y[k] * kp_76[k]
                       + lp_97[k];

            t_154[k] = ab_y[k] * kp_77[k]
                       + lp_98[k];

            t_155[k] = ab_z[k] * kp_77[k]
                       + lp_101[k];

            t_156[k] = ab_x[k] * kp_78[k]
                       + lp_78[k];
        }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, kp_79, kp_80, \
                         lp_79, lp_80, lp_100, lp_101, lp_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_157[k] = ab_x[k] * kp_79[k]
                       + lp_79[k];

            t_158[k] = ab_x[k] * kp_80[k]
                       + lp_80[k];

            t_159[k] = ab_y[k] * kp_79[k]
                       + lp_100[k];

            t_160[k] = ab_y[k] * kp_80[k]
                       + lp_101[k];

            t_161[k] = ab_z[k] * kp_80[k]
                       + lp_104[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, kp_81, kp_82, kp_83, \
                         lp_81, lp_82, lp_83, lp_103, lp_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * kp_81[k]
                       + lp_81[k];

            t_163[k] = ab_x[k] * kp_82[k]
                       + lp_82[k];

            t_164[k] = ab_x[k] * kp_83[k]
                       + lp_83[k];

            t_165[k] = ab_y[k] * kp_82[k]
                       + lp_103[k];

            t_166[k] = ab_y[k] * kp_83[k]
                       + lp_104[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, ab_x, ab_z, kp_83, kp_84, kp_85, kp_86, \
                         lp_84, lp_85, lp_86, lp_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = ab_z[k] * kp_83[k]
                       + lp_107[k];

            t_168[k] = ab_x[k] * kp_84[k]
                       + lp_84[k];

            t_169[k] = ab_x[k] * kp_85[k]
                       + lp_85[k];

            t_170[k] = ab_x[k] * kp_86[k]
                       + lp_86[k];
        }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, ab_x, ab_y, ab_z, kp_85, kp_86, kp_87, \
                         lp_87, lp_109, lp_110, lp_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_171[k] = ab_y[k] * kp_85[k]
                       + lp_109[k];

            t_172[k] = ab_y[k] * kp_86[k]
                       + lp_110[k];

            t_173[k] = ab_z[k] * kp_86[k]
                       + lp_113[k];

            t_174[k] = ab_x[k] * kp_87[k]
                       + lp_87[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, kp_88, kp_89, \
                         lp_88, lp_89, lp_112, lp_113, lp_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * kp_88[k]
                       + lp_88[k];

            t_176[k] = ab_x[k] * kp_89[k]
                       + lp_89[k];

            t_177[k] = ab_y[k] * kp_88[k]
                       + lp_112[k];

            t_178[k] = ab_y[k] * kp_89[k]
                       + lp_113[k];

            t_179[k] = ab_z[k] * kp_89[k]
                       + lp_116[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, kp_90, kp_91, kp_92, \
                         lp_90, lp_91, lp_92, lp_115, lp_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * kp_90[k]
                       + lp_90[k];

            t_181[k] = ab_x[k] * kp_91[k]
                       + lp_91[k];

            t_182[k] = ab_x[k] * kp_92[k]
                       + lp_92[k];

            t_183[k] = ab_y[k] * kp_91[k]
                       + lp_115[k];

            t_184[k] = ab_y[k] * kp_92[k]
                       + lp_116[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, ab_x, ab_z, kp_92, kp_93, kp_94, kp_95, \
                         lp_93, lp_94, lp_95, lp_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_z[k] * kp_92[k]
                       + lp_119[k];

            t_186[k] = ab_x[k] * kp_93[k]
                       + lp_93[k];

            t_187[k] = ab_x[k] * kp_94[k]
                       + lp_94[k];

            t_188[k] = ab_x[k] * kp_95[k]
                       + lp_95[k];
        }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, ab_x, ab_y, ab_z, kp_94, kp_95, kp_96, \
                         lp_96, lp_118, lp_119, lp_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_189[k] = ab_y[k] * kp_94[k]
                       + lp_118[k];

            t_190[k] = ab_y[k] * kp_95[k]
                       + lp_119[k];

            t_191[k] = ab_z[k] * kp_95[k]
                       + lp_122[k];

            t_192[k] = ab_x[k] * kp_96[k]
                       + lp_96[k];
        }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_x, ab_y, ab_z, kp_97, kp_98, \
                         lp_97, lp_98, lp_121, lp_122, lp_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_193[k] = ab_x[k] * kp_97[k]
                       + lp_97[k];

            t_194[k] = ab_x[k] * kp_98[k]
                       + lp_98[k];

            t_195[k] = ab_y[k] * kp_97[k]
                       + lp_121[k];

            t_196[k] = ab_y[k] * kp_98[k]
                       + lp_122[k];

            t_197[k] = ab_z[k] * kp_98[k]
                       + lp_125[k];
        }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, ab_x, ab_y, kp_99, kp_100, kp_101, \
                         lp_99, lp_100, lp_101, lp_124, lp_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_198[k] = ab_x[k] * kp_99[k]
                       + lp_99[k];

            t_199[k] = ab_x[k] * kp_100[k]
                       + lp_100[k];

            t_200[k] = ab_x[k] * kp_101[k]
                       + lp_101[k];

            t_201[k] = ab_y[k] * kp_100[k]
                       + lp_124[k];

            t_202[k] = ab_y[k] * kp_101[k]
                       + lp_125[k];
        }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, ab_x, ab_z, kp_101, kp_102, kp_103, \
                         kp_104, lp_102, lp_103, lp_104, lp_128 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_203[k] = ab_z[k] * kp_101[k]
                       + lp_128[k];

            t_204[k] = ab_x[k] * kp_102[k]
                       + lp_102[k];

            t_205[k] = ab_x[k] * kp_103[k]
                       + lp_103[k];

            t_206[k] = ab_x[k] * kp_104[k]
                       + lp_104[k];
        }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, ab_x, ab_y, ab_z, kp_103, kp_104, kp_105, \
                         lp_105, lp_127, lp_128, lp_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_207[k] = ab_y[k] * kp_103[k]
                       + lp_127[k];

            t_208[k] = ab_y[k] * kp_104[k]
                       + lp_128[k];

            t_209[k] = ab_z[k] * kp_104[k]
                       + lp_131[k];

            t_210[k] = ab_x[k] * kp_105[k]
                       + lp_105[k];
        }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ab_x, ab_y, ab_z, kp_106, kp_107, \
                         lp_106, lp_107, lp_130, lp_131, lp_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_211[k] = ab_x[k] * kp_106[k]
                       + lp_106[k];

            t_212[k] = ab_x[k] * kp_107[k]
                       + lp_107[k];

            t_213[k] = ab_y[k] * kp_106[k]
                       + lp_130[k];

            t_214[k] = ab_y[k] * kp_107[k]
                       + lp_131[k];

            t_215[k] = ab_z[k] * kp_107[k]
                       + lp_134[k];
        }
    }
}

auto
compute_hrr_kd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t kp, const size_t lp,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_kd_out_of_first_piece0(buffer, coordinates, target, kp, lp, ncomps, nmax);

    compute_hrr_kd_out_of_first_piece1(buffer, coordinates, target, kp, lp, ncomps, nmax);
}

static auto
compute_hrr_kd_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t kp, const size_t lp, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kp_0 = buffer.data(kp + 0 * ncomps + c);
        const auto *kp_1 = buffer.data(kp + 1 * ncomps + c);
        const auto *kp_2 = buffer.data(kp + 2 * ncomps + c);
        const auto *kp_3 = buffer.data(kp + 3 * ncomps + c);
        const auto *kp_4 = buffer.data(kp + 4 * ncomps + c);
        const auto *kp_5 = buffer.data(kp + 5 * ncomps + c);
        const auto *kp_6 = buffer.data(kp + 6 * ncomps + c);
        const auto *kp_7 = buffer.data(kp + 7 * ncomps + c);
        const auto *kp_8 = buffer.data(kp + 8 * ncomps + c);
        const auto *kp_9 = buffer.data(kp + 9 * ncomps + c);
        const auto *kp_10 = buffer.data(kp + 10 * ncomps + c);
        const auto *kp_11 = buffer.data(kp + 11 * ncomps + c);
        const auto *kp_12 = buffer.data(kp + 12 * ncomps + c);
        const auto *kp_13 = buffer.data(kp + 13 * ncomps + c);
        const auto *kp_14 = buffer.data(kp + 14 * ncomps + c);
        const auto *kp_15 = buffer.data(kp + 15 * ncomps + c);
        const auto *kp_16 = buffer.data(kp + 16 * ncomps + c);
        const auto *kp_17 = buffer.data(kp + 17 * ncomps + c);
        const auto *kp_18 = buffer.data(kp + 18 * ncomps + c);
        const auto *kp_19 = buffer.data(kp + 19 * ncomps + c);
        const auto *kp_20 = buffer.data(kp + 20 * ncomps + c);
        const auto *kp_21 = buffer.data(kp + 21 * ncomps + c);
        const auto *kp_22 = buffer.data(kp + 22 * ncomps + c);
        const auto *kp_23 = buffer.data(kp + 23 * ncomps + c);
        const auto *kp_24 = buffer.data(kp + 24 * ncomps + c);
        const auto *kp_25 = buffer.data(kp + 25 * ncomps + c);
        const auto *kp_26 = buffer.data(kp + 26 * ncomps + c);
        const auto *kp_27 = buffer.data(kp + 27 * ncomps + c);
        const auto *kp_28 = buffer.data(kp + 28 * ncomps + c);
        const auto *kp_29 = buffer.data(kp + 29 * ncomps + c);
        const auto *kp_30 = buffer.data(kp + 30 * ncomps + c);
        const auto *kp_31 = buffer.data(kp + 31 * ncomps + c);
        const auto *kp_32 = buffer.data(kp + 32 * ncomps + c);
        const auto *kp_33 = buffer.data(kp + 33 * ncomps + c);
        const auto *kp_34 = buffer.data(kp + 34 * ncomps + c);
        const auto *kp_35 = buffer.data(kp + 35 * ncomps + c);
        const auto *kp_36 = buffer.data(kp + 36 * ncomps + c);
        const auto *kp_37 = buffer.data(kp + 37 * ncomps + c);
        const auto *kp_38 = buffer.data(kp + 38 * ncomps + c);
        const auto *kp_39 = buffer.data(kp + 39 * ncomps + c);
        const auto *kp_40 = buffer.data(kp + 40 * ncomps + c);
        const auto *kp_41 = buffer.data(kp + 41 * ncomps + c);
        const auto *kp_42 = buffer.data(kp + 42 * ncomps + c);
        const auto *kp_43 = buffer.data(kp + 43 * ncomps + c);
        const auto *kp_44 = buffer.data(kp + 44 * ncomps + c);
        const auto *kp_45 = buffer.data(kp + 45 * ncomps + c);
        const auto *kp_46 = buffer.data(kp + 46 * ncomps + c);
        const auto *kp_47 = buffer.data(kp + 47 * ncomps + c);
        const auto *kp_48 = buffer.data(kp + 48 * ncomps + c);
        const auto *kp_49 = buffer.data(kp + 49 * ncomps + c);
        const auto *kp_50 = buffer.data(kp + 50 * ncomps + c);
        const auto *kp_51 = buffer.data(kp + 51 * ncomps + c);
        const auto *kp_52 = buffer.data(kp + 52 * ncomps + c);
        const auto *kp_53 = buffer.data(kp + 53 * ncomps + c);
        const auto *kp_54 = buffer.data(kp + 54 * ncomps + c);
        const auto *kp_55 = buffer.data(kp + 55 * ncomps + c);
        const auto *kp_56 = buffer.data(kp + 56 * ncomps + c);
        const auto *kp_57 = buffer.data(kp + 57 * ncomps + c);
        const auto *kp_58 = buffer.data(kp + 58 * ncomps + c);
        const auto *kp_59 = buffer.data(kp + 59 * ncomps + c);
        const auto *kp_60 = buffer.data(kp + 60 * ncomps + c);
        const auto *kp_61 = buffer.data(kp + 61 * ncomps + c);
        const auto *kp_62 = buffer.data(kp + 62 * ncomps + c);
        const auto *kp_63 = buffer.data(kp + 63 * ncomps + c);
        const auto *kp_64 = buffer.data(kp + 64 * ncomps + c);
        const auto *kp_65 = buffer.data(kp + 65 * ncomps + c);
        const auto *kp_66 = buffer.data(kp + 66 * ncomps + c);
        const auto *kp_67 = buffer.data(kp + 67 * ncomps + c);
        const auto *kp_68 = buffer.data(kp + 68 * ncomps + c);
        const auto *kp_69 = buffer.data(kp + 69 * ncomps + c);
        const auto *kp_70 = buffer.data(kp + 70 * ncomps + c);
        const auto *kp_71 = buffer.data(kp + 71 * ncomps + c);

        const auto *lp_0 = buffer.data(lp + 0 * ncomps + c);
        const auto *lp_1 = buffer.data(lp + 1 * ncomps + c);
        const auto *lp_2 = buffer.data(lp + 2 * ncomps + c);
        const auto *lp_3 = buffer.data(lp + 3 * ncomps + c);
        const auto *lp_4 = buffer.data(lp + 4 * ncomps + c);
        const auto *lp_5 = buffer.data(lp + 5 * ncomps + c);
        const auto *lp_6 = buffer.data(lp + 6 * ncomps + c);
        const auto *lp_7 = buffer.data(lp + 7 * ncomps + c);
        const auto *lp_8 = buffer.data(lp + 8 * ncomps + c);
        const auto *lp_9 = buffer.data(lp + 9 * ncomps + c);
        const auto *lp_10 = buffer.data(lp + 10 * ncomps + c);
        const auto *lp_11 = buffer.data(lp + 11 * ncomps + c);
        const auto *lp_12 = buffer.data(lp + 12 * ncomps + c);
        const auto *lp_13 = buffer.data(lp + 13 * ncomps + c);
        const auto *lp_14 = buffer.data(lp + 14 * ncomps + c);
        const auto *lp_15 = buffer.data(lp + 15 * ncomps + c);
        const auto *lp_16 = buffer.data(lp + 16 * ncomps + c);
        const auto *lp_17 = buffer.data(lp + 17 * ncomps + c);
        const auto *lp_18 = buffer.data(lp + 18 * ncomps + c);
        const auto *lp_19 = buffer.data(lp + 19 * ncomps + c);
        const auto *lp_20 = buffer.data(lp + 20 * ncomps + c);
        const auto *lp_21 = buffer.data(lp + 21 * ncomps + c);
        const auto *lp_22 = buffer.data(lp + 22 * ncomps + c);
        const auto *lp_23 = buffer.data(lp + 23 * ncomps + c);
        const auto *lp_24 = buffer.data(lp + 24 * ncomps + c);
        const auto *lp_25 = buffer.data(lp + 25 * ncomps + c);
        const auto *lp_26 = buffer.data(lp + 26 * ncomps + c);
        const auto *lp_27 = buffer.data(lp + 27 * ncomps + c);
        const auto *lp_28 = buffer.data(lp + 28 * ncomps + c);
        const auto *lp_29 = buffer.data(lp + 29 * ncomps + c);
        const auto *lp_30 = buffer.data(lp + 30 * ncomps + c);
        const auto *lp_31 = buffer.data(lp + 31 * ncomps + c);
        const auto *lp_32 = buffer.data(lp + 32 * ncomps + c);
        const auto *lp_33 = buffer.data(lp + 33 * ncomps + c);
        const auto *lp_34 = buffer.data(lp + 34 * ncomps + c);
        const auto *lp_35 = buffer.data(lp + 35 * ncomps + c);
        const auto *lp_36 = buffer.data(lp + 36 * ncomps + c);
        const auto *lp_37 = buffer.data(lp + 37 * ncomps + c);
        const auto *lp_38 = buffer.data(lp + 38 * ncomps + c);
        const auto *lp_39 = buffer.data(lp + 39 * ncomps + c);
        const auto *lp_40 = buffer.data(lp + 40 * ncomps + c);
        const auto *lp_41 = buffer.data(lp + 41 * ncomps + c);
        const auto *lp_42 = buffer.data(lp + 42 * ncomps + c);
        const auto *lp_43 = buffer.data(lp + 43 * ncomps + c);
        const auto *lp_44 = buffer.data(lp + 44 * ncomps + c);
        const auto *lp_45 = buffer.data(lp + 45 * ncomps + c);
        const auto *lp_46 = buffer.data(lp + 46 * ncomps + c);
        const auto *lp_47 = buffer.data(lp + 47 * ncomps + c);
        const auto *lp_48 = buffer.data(lp + 48 * ncomps + c);
        const auto *lp_49 = buffer.data(lp + 49 * ncomps + c);
        const auto *lp_50 = buffer.data(lp + 50 * ncomps + c);
        const auto *lp_51 = buffer.data(lp + 51 * ncomps + c);
        const auto *lp_52 = buffer.data(lp + 52 * ncomps + c);
        const auto *lp_53 = buffer.data(lp + 53 * ncomps + c);
        const auto *lp_54 = buffer.data(lp + 54 * ncomps + c);
        const auto *lp_55 = buffer.data(lp + 55 * ncomps + c);
        const auto *lp_56 = buffer.data(lp + 56 * ncomps + c);
        const auto *lp_57 = buffer.data(lp + 57 * ncomps + c);
        const auto *lp_58 = buffer.data(lp + 58 * ncomps + c);
        const auto *lp_59 = buffer.data(lp + 59 * ncomps + c);
        const auto *lp_60 = buffer.data(lp + 60 * ncomps + c);
        const auto *lp_61 = buffer.data(lp + 61 * ncomps + c);
        const auto *lp_62 = buffer.data(lp + 62 * ncomps + c);
        const auto *lp_63 = buffer.data(lp + 63 * ncomps + c);
        const auto *lp_64 = buffer.data(lp + 64 * ncomps + c);
        const auto *lp_65 = buffer.data(lp + 65 * ncomps + c);
        const auto *lp_66 = buffer.data(lp + 66 * ncomps + c);
        const auto *lp_67 = buffer.data(lp + 67 * ncomps + c);
        const auto *lp_68 = buffer.data(lp + 68 * ncomps + c);
        const auto *lp_69 = buffer.data(lp + 69 * ncomps + c);
        const auto *lp_70 = buffer.data(lp + 70 * ncomps + c);
        const auto *lp_71 = buffer.data(lp + 71 * ncomps + c);
        const auto *lp_73 = buffer.data(lp + 73 * ncomps + c);
        const auto *lp_74 = buffer.data(lp + 74 * ncomps + c);
        const auto *lp_76 = buffer.data(lp + 76 * ncomps + c);
        const auto *lp_77 = buffer.data(lp + 77 * ncomps + c);
        const auto *lp_79 = buffer.data(lp + 79 * ncomps + c);
        const auto *lp_80 = buffer.data(lp + 80 * ncomps + c);
        const auto *lp_83 = buffer.data(lp + 83 * ncomps + c);
        const auto *lp_85 = buffer.data(lp + 85 * ncomps + c);
        const auto *lp_86 = buffer.data(lp + 86 * ncomps + c);
        const auto *lp_88 = buffer.data(lp + 88 * ncomps + c);
        const auto *lp_89 = buffer.data(lp + 89 * ncomps + c);
        const auto *lp_91 = buffer.data(lp + 91 * ncomps + c);
        const auto *lp_92 = buffer.data(lp + 92 * ncomps + c);
        const auto *lp_95 = buffer.data(lp + 95 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, kp_0, kp_1, kp_2, lp_0, lp_1, \
                         lp_2, lp_4, lp_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * kp_0[k]
                     + lp_0[k];

            t_1[k] = ab_x[k] * kp_1[k]
                     + lp_1[k];

            t_2[k] = ab_x[k] * kp_2[k]
                     + lp_2[k];

            t_3[k] = ab_y[k] * kp_1[k]
                     + lp_4[k];

            t_4[k] = ab_y[k] * kp_2[k]
                     + lp_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, kp_2, kp_3, kp_4, kp_5, lp_3, lp_4, \
                         lp_5, lp_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * kp_2[k]
                     + lp_8[k];

            t_6[k] = ab_x[k] * kp_3[k]
                     + lp_3[k];

            t_7[k] = ab_x[k] * kp_4[k]
                     + lp_4[k];

            t_8[k] = ab_x[k] * kp_5[k]
                     + lp_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, kp_4, kp_5, kp_6, lp_6, \
                         lp_10, lp_11, lp_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * kp_4[k]
                     + lp_10[k];

            t_10[k] = ab_y[k] * kp_5[k]
                      + lp_11[k];

            t_11[k] = ab_z[k] * kp_5[k]
                      + lp_14[k];

            t_12[k] = ab_x[k] * kp_6[k]
                      + lp_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, kp_7, kp_8, lp_7, \
                         lp_8, lp_13, lp_14, lp_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * kp_7[k]
                      + lp_7[k];

            t_14[k] = ab_x[k] * kp_8[k]
                      + lp_8[k];

            t_15[k] = ab_y[k] * kp_7[k]
                      + lp_13[k];

            t_16[k] = ab_y[k] * kp_8[k]
                      + lp_14[k];

            t_17[k] = ab_z[k] * kp_8[k]
                      + lp_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, kp_9, kp_10, kp_11, lp_9, \
                         lp_10, lp_11, lp_19, lp_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * kp_9[k]
                      + lp_9[k];

            t_19[k] = ab_x[k] * kp_10[k]
                      + lp_10[k];

            t_20[k] = ab_x[k] * kp_11[k]
                      + lp_11[k];

            t_21[k] = ab_y[k] * kp_10[k]
                      + lp_19[k];

            t_22[k] = ab_y[k] * kp_11[k]
                      + lp_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, kp_11, kp_12, kp_13, kp_14, \
                         lp_12, lp_13, lp_14, lp_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * kp_11[k]
                      + lp_23[k];

            t_24[k] = ab_x[k] * kp_12[k]
                      + lp_12[k];

            t_25[k] = ab_x[k] * kp_13[k]
                      + lp_13[k];

            t_26[k] = ab_x[k] * kp_14[k]
                      + lp_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, kp_13, kp_14, kp_15, lp_15, \
                         lp_22, lp_23, lp_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * kp_13[k]
                      + lp_22[k];

            t_28[k] = ab_y[k] * kp_14[k]
                      + lp_23[k];

            t_29[k] = ab_z[k] * kp_14[k]
                      + lp_26[k];

            t_30[k] = ab_x[k] * kp_15[k]
                      + lp_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, kp_16, kp_17, lp_16, \
                         lp_17, lp_25, lp_26, lp_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * kp_16[k]
                      + lp_16[k];

            t_32[k] = ab_x[k] * kp_17[k]
                      + lp_17[k];

            t_33[k] = ab_y[k] * kp_16[k]
                      + lp_25[k];

            t_34[k] = ab_y[k] * kp_17[k]
                      + lp_26[k];

            t_35[k] = ab_z[k] * kp_17[k]
                      + lp_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, kp_18, kp_19, kp_20, lp_18, \
                         lp_19, lp_20, lp_31, lp_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * kp_18[k]
                      + lp_18[k];

            t_37[k] = ab_x[k] * kp_19[k]
                      + lp_19[k];

            t_38[k] = ab_x[k] * kp_20[k]
                      + lp_20[k];

            t_39[k] = ab_y[k] * kp_19[k]
                      + lp_31[k];

            t_40[k] = ab_y[k] * kp_20[k]
                      + lp_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, kp_20, kp_21, kp_22, kp_23, \
                         lp_21, lp_22, lp_23, lp_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * kp_20[k]
                      + lp_35[k];

            t_42[k] = ab_x[k] * kp_21[k]
                      + lp_21[k];

            t_43[k] = ab_x[k] * kp_22[k]
                      + lp_22[k];

            t_44[k] = ab_x[k] * kp_23[k]
                      + lp_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, kp_22, kp_23, kp_24, lp_24, \
                         lp_34, lp_35, lp_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * kp_22[k]
                      + lp_34[k];

            t_46[k] = ab_y[k] * kp_23[k]
                      + lp_35[k];

            t_47[k] = ab_z[k] * kp_23[k]
                      + lp_38[k];

            t_48[k] = ab_x[k] * kp_24[k]
                      + lp_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, kp_25, kp_26, lp_25, \
                         lp_26, lp_37, lp_38, lp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * kp_25[k]
                      + lp_25[k];

            t_50[k] = ab_x[k] * kp_26[k]
                      + lp_26[k];

            t_51[k] = ab_y[k] * kp_25[k]
                      + lp_37[k];

            t_52[k] = ab_y[k] * kp_26[k]
                      + lp_38[k];

            t_53[k] = ab_z[k] * kp_26[k]
                      + lp_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, kp_27, kp_28, kp_29, lp_27, \
                         lp_28, lp_29, lp_40, lp_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * kp_27[k]
                      + lp_27[k];

            t_55[k] = ab_x[k] * kp_28[k]
                      + lp_28[k];

            t_56[k] = ab_x[k] * kp_29[k]
                      + lp_29[k];

            t_57[k] = ab_y[k] * kp_28[k]
                      + lp_40[k];

            t_58[k] = ab_y[k] * kp_29[k]
                      + lp_41[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, kp_29, kp_30, kp_31, kp_32, \
                         lp_30, lp_31, lp_32, lp_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * kp_29[k]
                      + lp_44[k];

            t_60[k] = ab_x[k] * kp_30[k]
                      + lp_30[k];

            t_61[k] = ab_x[k] * kp_31[k]
                      + lp_31[k];

            t_62[k] = ab_x[k] * kp_32[k]
                      + lp_32[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, kp_31, kp_32, kp_33, lp_33, \
                         lp_46, lp_47, lp_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = ab_y[k] * kp_31[k]
                      + lp_46[k];

            t_64[k] = ab_y[k] * kp_32[k]
                      + lp_47[k];

            t_65[k] = ab_z[k] * kp_32[k]
                      + lp_50[k];

            t_66[k] = ab_x[k] * kp_33[k]
                      + lp_33[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, kp_34, kp_35, lp_34, \
                         lp_35, lp_49, lp_50, lp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = ab_x[k] * kp_34[k]
                      + lp_34[k];

            t_68[k] = ab_x[k] * kp_35[k]
                      + lp_35[k];

            t_69[k] = ab_y[k] * kp_34[k]
                      + lp_49[k];

            t_70[k] = ab_y[k] * kp_35[k]
                      + lp_50[k];

            t_71[k] = ab_z[k] * kp_35[k]
                      + lp_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, kp_36, kp_37, kp_38, lp_36, \
                         lp_37, lp_38, lp_52, lp_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * kp_36[k]
                      + lp_36[k];

            t_73[k] = ab_x[k] * kp_37[k]
                      + lp_37[k];

            t_74[k] = ab_x[k] * kp_38[k]
                      + lp_38[k];

            t_75[k] = ab_y[k] * kp_37[k]
                      + lp_52[k];

            t_76[k] = ab_y[k] * kp_38[k]
                      + lp_53[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, kp_38, kp_39, kp_40, kp_41, \
                         lp_39, lp_40, lp_41, lp_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * kp_38[k]
                      + lp_56[k];

            t_78[k] = ab_x[k] * kp_39[k]
                      + lp_39[k];

            t_79[k] = ab_x[k] * kp_40[k]
                      + lp_40[k];

            t_80[k] = ab_x[k] * kp_41[k]
                      + lp_41[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, kp_40, kp_41, kp_42, lp_42, \
                         lp_55, lp_56, lp_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = ab_y[k] * kp_40[k]
                      + lp_55[k];

            t_82[k] = ab_y[k] * kp_41[k]
                      + lp_56[k];

            t_83[k] = ab_z[k] * kp_41[k]
                      + lp_59[k];

            t_84[k] = ab_x[k] * kp_42[k]
                      + lp_42[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, kp_43, kp_44, lp_43, \
                         lp_44, lp_58, lp_59, lp_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * kp_43[k]
                      + lp_43[k];

            t_86[k] = ab_x[k] * kp_44[k]
                      + lp_44[k];

            t_87[k] = ab_y[k] * kp_43[k]
                      + lp_58[k];

            t_88[k] = ab_y[k] * kp_44[k]
                      + lp_59[k];

            t_89[k] = ab_z[k] * kp_44[k]
                      + lp_62[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, kp_45, kp_46, kp_47, lp_45, \
                         lp_46, lp_47, lp_64, lp_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * kp_45[k]
                      + lp_45[k];

            t_91[k] = ab_x[k] * kp_46[k]
                      + lp_46[k];

            t_92[k] = ab_x[k] * kp_47[k]
                      + lp_47[k];

            t_93[k] = ab_y[k] * kp_46[k]
                      + lp_64[k];

            t_94[k] = ab_y[k] * kp_47[k]
                      + lp_65[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, kp_47, kp_48, kp_49, kp_50, \
                         lp_48, lp_49, lp_50, lp_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_z[k] * kp_47[k]
                      + lp_68[k];

            t_96[k] = ab_x[k] * kp_48[k]
                      + lp_48[k];

            t_97[k] = ab_x[k] * kp_49[k]
                      + lp_49[k];

            t_98[k] = ab_x[k] * kp_50[k]
                      + lp_50[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, kp_49, kp_50, kp_51, \
                         lp_51, lp_67, lp_68, lp_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * kp_49[k]
                      + lp_67[k];

            t_100[k] = ab_y[k] * kp_50[k]
                       + lp_68[k];

            t_101[k] = ab_z[k] * kp_50[k]
                       + lp_71[k];

            t_102[k] = ab_x[k] * kp_51[k]
                       + lp_51[k];
        }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, kp_52, kp_53, \
                         lp_52, lp_53, lp_70, lp_71, lp_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_103[k] = ab_x[k] * kp_52[k]
                       + lp_52[k];

            t_104[k] = ab_x[k] * kp_53[k]
                       + lp_53[k];

            t_105[k] = ab_y[k] * kp_52[k]
                       + lp_70[k];

            t_106[k] = ab_y[k] * kp_53[k]
                       + lp_71[k];

            t_107[k] = ab_z[k] * kp_53[k]
                       + lp_74[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, kp_54, kp_55, kp_56, \
                         lp_54, lp_55, lp_56, lp_73, lp_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * kp_54[k]
                       + lp_54[k];

            t_109[k] = ab_x[k] * kp_55[k]
                       + lp_55[k];

            t_110[k] = ab_x[k] * kp_56[k]
                       + lp_56[k];

            t_111[k] = ab_y[k] * kp_55[k]
                       + lp_73[k];

            t_112[k] = ab_y[k] * kp_56[k]
                       + lp_74[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, kp_56, kp_57, kp_58, kp_59, \
                         lp_57, lp_58, lp_59, lp_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_z[k] * kp_56[k]
                       + lp_77[k];

            t_114[k] = ab_x[k] * kp_57[k]
                       + lp_57[k];

            t_115[k] = ab_x[k] * kp_58[k]
                       + lp_58[k];

            t_116[k] = ab_x[k] * kp_59[k]
                       + lp_59[k];
        }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, kp_58, kp_59, kp_60, \
                         lp_60, lp_76, lp_77, lp_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_117[k] = ab_y[k] * kp_58[k]
                       + lp_76[k];

            t_118[k] = ab_y[k] * kp_59[k]
                       + lp_77[k];

            t_119[k] = ab_z[k] * kp_59[k]
                       + lp_80[k];

            t_120[k] = ab_x[k] * kp_60[k]
                       + lp_60[k];
        }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, kp_61, kp_62, \
                         lp_61, lp_62, lp_79, lp_80, lp_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_121[k] = ab_x[k] * kp_61[k]
                       + lp_61[k];

            t_122[k] = ab_x[k] * kp_62[k]
                       + lp_62[k];

            t_123[k] = ab_y[k] * kp_61[k]
                       + lp_79[k];

            t_124[k] = ab_y[k] * kp_62[k]
                       + lp_80[k];

            t_125[k] = ab_z[k] * kp_62[k]
                       + lp_83[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ab_y, kp_63, kp_64, kp_65, \
                         lp_63, lp_64, lp_65, lp_85, lp_86 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = ab_x[k] * kp_63[k]
                       + lp_63[k];

            t_127[k] = ab_x[k] * kp_64[k]
                       + lp_64[k];

            t_128[k] = ab_x[k] * kp_65[k]
                       + lp_65[k];

            t_129[k] = ab_y[k] * kp_64[k]
                       + lp_85[k];

            t_130[k] = ab_y[k] * kp_65[k]
                       + lp_86[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, ab_x, ab_z, kp_65, kp_66, kp_67, kp_68, \
                         lp_66, lp_67, lp_68, lp_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = ab_z[k] * kp_65[k]
                       + lp_89[k];

            t_132[k] = ab_x[k] * kp_66[k]
                       + lp_66[k];

            t_133[k] = ab_x[k] * kp_67[k]
                       + lp_67[k];

            t_134[k] = ab_x[k] * kp_68[k]
                       + lp_68[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, kp_67, kp_68, kp_69, \
                         lp_69, lp_88, lp_89, lp_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_y[k] * kp_67[k]
                       + lp_88[k];

            t_136[k] = ab_y[k] * kp_68[k]
                       + lp_89[k];

            t_137[k] = ab_z[k] * kp_68[k]
                       + lp_92[k];

            t_138[k] = ab_x[k] * kp_69[k]
                       + lp_69[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, kp_70, kp_71, \
                         lp_70, lp_71, lp_91, lp_92, lp_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_x[k] * kp_70[k]
                       + lp_70[k];

            t_140[k] = ab_x[k] * kp_71[k]
                       + lp_71[k];

            t_141[k] = ab_y[k] * kp_70[k]
                       + lp_91[k];

            t_142[k] = ab_y[k] * kp_71[k]
                       + lp_92[k];

            t_143[k] = ab_z[k] * kp_71[k]
                       + lp_95[k];
        }
    }
}

static auto
compute_hrr_kd_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t kp, const size_t lp, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_144 = buffer.data(target + 144 * ncomps + c);
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kp_72 = buffer.data(kp + 72 * ncomps + c);
        const auto *kp_73 = buffer.data(kp + 73 * ncomps + c);
        const auto *kp_74 = buffer.data(kp + 74 * ncomps + c);
        const auto *kp_75 = buffer.data(kp + 75 * ncomps + c);
        const auto *kp_76 = buffer.data(kp + 76 * ncomps + c);
        const auto *kp_77 = buffer.data(kp + 77 * ncomps + c);
        const auto *kp_78 = buffer.data(kp + 78 * ncomps + c);
        const auto *kp_79 = buffer.data(kp + 79 * ncomps + c);
        const auto *kp_80 = buffer.data(kp + 80 * ncomps + c);
        const auto *kp_81 = buffer.data(kp + 81 * ncomps + c);
        const auto *kp_82 = buffer.data(kp + 82 * ncomps + c);
        const auto *kp_83 = buffer.data(kp + 83 * ncomps + c);
        const auto *kp_84 = buffer.data(kp + 84 * ncomps + c);
        const auto *kp_85 = buffer.data(kp + 85 * ncomps + c);
        const auto *kp_86 = buffer.data(kp + 86 * ncomps + c);
        const auto *kp_87 = buffer.data(kp + 87 * ncomps + c);
        const auto *kp_88 = buffer.data(kp + 88 * ncomps + c);
        const auto *kp_89 = buffer.data(kp + 89 * ncomps + c);
        const auto *kp_90 = buffer.data(kp + 90 * ncomps + c);
        const auto *kp_91 = buffer.data(kp + 91 * ncomps + c);
        const auto *kp_92 = buffer.data(kp + 92 * ncomps + c);
        const auto *kp_93 = buffer.data(kp + 93 * ncomps + c);
        const auto *kp_94 = buffer.data(kp + 94 * ncomps + c);
        const auto *kp_95 = buffer.data(kp + 95 * ncomps + c);
        const auto *kp_96 = buffer.data(kp + 96 * ncomps + c);
        const auto *kp_97 = buffer.data(kp + 97 * ncomps + c);
        const auto *kp_98 = buffer.data(kp + 98 * ncomps + c);
        const auto *kp_99 = buffer.data(kp + 99 * ncomps + c);
        const auto *kp_100 = buffer.data(kp + 100 * ncomps + c);
        const auto *kp_101 = buffer.data(kp + 101 * ncomps + c);
        const auto *kp_102 = buffer.data(kp + 102 * ncomps + c);
        const auto *kp_103 = buffer.data(kp + 103 * ncomps + c);
        const auto *kp_104 = buffer.data(kp + 104 * ncomps + c);
        const auto *kp_105 = buffer.data(kp + 105 * ncomps + c);
        const auto *kp_106 = buffer.data(kp + 106 * ncomps + c);
        const auto *kp_107 = buffer.data(kp + 107 * ncomps + c);

        const auto *lp_72 = buffer.data(lp + 72 * ncomps + c);
        const auto *lp_73 = buffer.data(lp + 73 * ncomps + c);
        const auto *lp_74 = buffer.data(lp + 74 * ncomps + c);
        const auto *lp_75 = buffer.data(lp + 75 * ncomps + c);
        const auto *lp_76 = buffer.data(lp + 76 * ncomps + c);
        const auto *lp_77 = buffer.data(lp + 77 * ncomps + c);
        const auto *lp_78 = buffer.data(lp + 78 * ncomps + c);
        const auto *lp_79 = buffer.data(lp + 79 * ncomps + c);
        const auto *lp_80 = buffer.data(lp + 80 * ncomps + c);
        const auto *lp_81 = buffer.data(lp + 81 * ncomps + c);
        const auto *lp_82 = buffer.data(lp + 82 * ncomps + c);
        const auto *lp_83 = buffer.data(lp + 83 * ncomps + c);
        const auto *lp_84 = buffer.data(lp + 84 * ncomps + c);
        const auto *lp_85 = buffer.data(lp + 85 * ncomps + c);
        const auto *lp_86 = buffer.data(lp + 86 * ncomps + c);
        const auto *lp_87 = buffer.data(lp + 87 * ncomps + c);
        const auto *lp_88 = buffer.data(lp + 88 * ncomps + c);
        const auto *lp_89 = buffer.data(lp + 89 * ncomps + c);
        const auto *lp_90 = buffer.data(lp + 90 * ncomps + c);
        const auto *lp_91 = buffer.data(lp + 91 * ncomps + c);
        const auto *lp_92 = buffer.data(lp + 92 * ncomps + c);
        const auto *lp_93 = buffer.data(lp + 93 * ncomps + c);
        const auto *lp_94 = buffer.data(lp + 94 * ncomps + c);
        const auto *lp_95 = buffer.data(lp + 95 * ncomps + c);
        const auto *lp_96 = buffer.data(lp + 96 * ncomps + c);
        const auto *lp_97 = buffer.data(lp + 97 * ncomps + c);
        const auto *lp_98 = buffer.data(lp + 98 * ncomps + c);
        const auto *lp_99 = buffer.data(lp + 99 * ncomps + c);
        const auto *lp_100 = buffer.data(lp + 100 * ncomps + c);
        const auto *lp_101 = buffer.data(lp + 101 * ncomps + c);
        const auto *lp_102 = buffer.data(lp + 102 * ncomps + c);
        const auto *lp_103 = buffer.data(lp + 103 * ncomps + c);
        const auto *lp_104 = buffer.data(lp + 104 * ncomps + c);
        const auto *lp_105 = buffer.data(lp + 105 * ncomps + c);
        const auto *lp_106 = buffer.data(lp + 106 * ncomps + c);
        const auto *lp_107 = buffer.data(lp + 107 * ncomps + c);
        const auto *lp_109 = buffer.data(lp + 109 * ncomps + c);
        const auto *lp_110 = buffer.data(lp + 110 * ncomps + c);
        const auto *lp_112 = buffer.data(lp + 112 * ncomps + c);
        const auto *lp_113 = buffer.data(lp + 113 * ncomps + c);
        const auto *lp_115 = buffer.data(lp + 115 * ncomps + c);
        const auto *lp_116 = buffer.data(lp + 116 * ncomps + c);
        const auto *lp_118 = buffer.data(lp + 118 * ncomps + c);
        const auto *lp_119 = buffer.data(lp + 119 * ncomps + c);
        const auto *lp_121 = buffer.data(lp + 121 * ncomps + c);
        const auto *lp_122 = buffer.data(lp + 122 * ncomps + c);
        const auto *lp_124 = buffer.data(lp + 124 * ncomps + c);
        const auto *lp_125 = buffer.data(lp + 125 * ncomps + c);
        const auto *lp_127 = buffer.data(lp + 127 * ncomps + c);
        const auto *lp_128 = buffer.data(lp + 128 * ncomps + c);
        const auto *lp_130 = buffer.data(lp + 130 * ncomps + c);
        const auto *lp_131 = buffer.data(lp + 131 * ncomps + c);
        const auto *lp_134 = buffer.data(lp + 134 * ncomps + c);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_x, ab_y, kp_72, kp_73, kp_74, \
                         lp_72, lp_73, lp_74, lp_94, lp_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * kp_72[k]
                       + lp_72[k];

            t_145[k] = ab_x[k] * kp_73[k]
                       + lp_73[k];

            t_146[k] = ab_x[k] * kp_74[k]
                       + lp_74[k];

            t_147[k] = ab_y[k] * kp_73[k]
                       + lp_94[k];

            t_148[k] = ab_y[k] * kp_74[k]
                       + lp_95[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, ab_x, ab_z, kp_74, kp_75, kp_76, kp_77, \
                         lp_75, lp_76, lp_77, lp_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = ab_z[k] * kp_74[k]
                       + lp_98[k];

            t_150[k] = ab_x[k] * kp_75[k]
                       + lp_75[k];

            t_151[k] = ab_x[k] * kp_76[k]
                       + lp_76[k];

            t_152[k] = ab_x[k] * kp_77[k]
                       + lp_77[k];
        }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, ab_x, ab_y, ab_z, kp_76, kp_77, kp_78, \
                         lp_78, lp_97, lp_98, lp_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_153[k] = ab_y[k] * kp_76[k]
                       + lp_97[k];

            t_154[k] = ab_y[k] * kp_77[k]
                       + lp_98[k];

            t_155[k] = ab_z[k] * kp_77[k]
                       + lp_101[k];

            t_156[k] = ab_x[k] * kp_78[k]
                       + lp_78[k];
        }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, kp_79, kp_80, \
                         lp_79, lp_80, lp_100, lp_101, lp_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_157[k] = ab_x[k] * kp_79[k]
                       + lp_79[k];

            t_158[k] = ab_x[k] * kp_80[k]
                       + lp_80[k];

            t_159[k] = ab_y[k] * kp_79[k]
                       + lp_100[k];

            t_160[k] = ab_y[k] * kp_80[k]
                       + lp_101[k];

            t_161[k] = ab_z[k] * kp_80[k]
                       + lp_104[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, kp_81, kp_82, kp_83, \
                         lp_81, lp_82, lp_83, lp_103, lp_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * kp_81[k]
                       + lp_81[k];

            t_163[k] = ab_x[k] * kp_82[k]
                       + lp_82[k];

            t_164[k] = ab_x[k] * kp_83[k]
                       + lp_83[k];

            t_165[k] = ab_y[k] * kp_82[k]
                       + lp_103[k];

            t_166[k] = ab_y[k] * kp_83[k]
                       + lp_104[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, ab_x, ab_z, kp_83, kp_84, kp_85, kp_86, \
                         lp_84, lp_85, lp_86, lp_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = ab_z[k] * kp_83[k]
                       + lp_107[k];

            t_168[k] = ab_x[k] * kp_84[k]
                       + lp_84[k];

            t_169[k] = ab_x[k] * kp_85[k]
                       + lp_85[k];

            t_170[k] = ab_x[k] * kp_86[k]
                       + lp_86[k];
        }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, ab_x, ab_y, ab_z, kp_85, kp_86, kp_87, \
                         lp_87, lp_109, lp_110, lp_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_171[k] = ab_y[k] * kp_85[k]
                       + lp_109[k];

            t_172[k] = ab_y[k] * kp_86[k]
                       + lp_110[k];

            t_173[k] = ab_z[k] * kp_86[k]
                       + lp_113[k];

            t_174[k] = ab_x[k] * kp_87[k]
                       + lp_87[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, kp_88, kp_89, \
                         lp_88, lp_89, lp_112, lp_113, lp_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * kp_88[k]
                       + lp_88[k];

            t_176[k] = ab_x[k] * kp_89[k]
                       + lp_89[k];

            t_177[k] = ab_y[k] * kp_88[k]
                       + lp_112[k];

            t_178[k] = ab_y[k] * kp_89[k]
                       + lp_113[k];

            t_179[k] = ab_z[k] * kp_89[k]
                       + lp_116[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, kp_90, kp_91, kp_92, \
                         lp_90, lp_91, lp_92, lp_115, lp_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * kp_90[k]
                       + lp_90[k];

            t_181[k] = ab_x[k] * kp_91[k]
                       + lp_91[k];

            t_182[k] = ab_x[k] * kp_92[k]
                       + lp_92[k];

            t_183[k] = ab_y[k] * kp_91[k]
                       + lp_115[k];

            t_184[k] = ab_y[k] * kp_92[k]
                       + lp_116[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, ab_x, ab_z, kp_92, kp_93, kp_94, kp_95, \
                         lp_93, lp_94, lp_95, lp_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_z[k] * kp_92[k]
                       + lp_119[k];

            t_186[k] = ab_x[k] * kp_93[k]
                       + lp_93[k];

            t_187[k] = ab_x[k] * kp_94[k]
                       + lp_94[k];

            t_188[k] = ab_x[k] * kp_95[k]
                       + lp_95[k];
        }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, ab_x, ab_y, ab_z, kp_94, kp_95, kp_96, \
                         lp_96, lp_118, lp_119, lp_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_189[k] = ab_y[k] * kp_94[k]
                       + lp_118[k];

            t_190[k] = ab_y[k] * kp_95[k]
                       + lp_119[k];

            t_191[k] = ab_z[k] * kp_95[k]
                       + lp_122[k];

            t_192[k] = ab_x[k] * kp_96[k]
                       + lp_96[k];
        }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_x, ab_y, ab_z, kp_97, kp_98, \
                         lp_97, lp_98, lp_121, lp_122, lp_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_193[k] = ab_x[k] * kp_97[k]
                       + lp_97[k];

            t_194[k] = ab_x[k] * kp_98[k]
                       + lp_98[k];

            t_195[k] = ab_y[k] * kp_97[k]
                       + lp_121[k];

            t_196[k] = ab_y[k] * kp_98[k]
                       + lp_122[k];

            t_197[k] = ab_z[k] * kp_98[k]
                       + lp_125[k];
        }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, ab_x, ab_y, kp_99, kp_100, kp_101, \
                         lp_99, lp_100, lp_101, lp_124, lp_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_198[k] = ab_x[k] * kp_99[k]
                       + lp_99[k];

            t_199[k] = ab_x[k] * kp_100[k]
                       + lp_100[k];

            t_200[k] = ab_x[k] * kp_101[k]
                       + lp_101[k];

            t_201[k] = ab_y[k] * kp_100[k]
                       + lp_124[k];

            t_202[k] = ab_y[k] * kp_101[k]
                       + lp_125[k];
        }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, ab_x, ab_z, kp_101, kp_102, kp_103, \
                         kp_104, lp_102, lp_103, lp_104, lp_128 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_203[k] = ab_z[k] * kp_101[k]
                       + lp_128[k];

            t_204[k] = ab_x[k] * kp_102[k]
                       + lp_102[k];

            t_205[k] = ab_x[k] * kp_103[k]
                       + lp_103[k];

            t_206[k] = ab_x[k] * kp_104[k]
                       + lp_104[k];
        }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, ab_x, ab_y, ab_z, kp_103, kp_104, kp_105, \
                         lp_105, lp_127, lp_128, lp_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_207[k] = ab_y[k] * kp_103[k]
                       + lp_127[k];

            t_208[k] = ab_y[k] * kp_104[k]
                       + lp_128[k];

            t_209[k] = ab_z[k] * kp_104[k]
                       + lp_131[k];

            t_210[k] = ab_x[k] * kp_105[k]
                       + lp_105[k];
        }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ab_x, ab_y, ab_z, kp_106, kp_107, \
                         lp_106, lp_107, lp_130, lp_131, lp_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_211[k] = ab_x[k] * kp_106[k]
                       + lp_106[k];

            t_212[k] = ab_x[k] * kp_107[k]
                       + lp_107[k];

            t_213[k] = ab_y[k] * kp_106[k]
                       + lp_130[k];

            t_214[k] = ab_y[k] * kp_107[k]
                       + lp_131[k];

            t_215[k] = ab_z[k] * kp_107[k]
                       + lp_134[k];
        }
    }
}

auto
compute_hrr_kd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t kp, const size_t lp, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_kd_piece0(buffer, coordinates, target, kp, lp, ncomps, nmax);

    compute_hrr_kd_piece1(buffer, coordinates, target, kp, lp, ncomps, nmax);
}

}  // namespace simdtrf
