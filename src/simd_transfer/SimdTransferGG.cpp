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


#include "SimdTransferGG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_gg_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t gf, const size_t hf,
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

        const auto *gf_0 = buffer.data(gf + 0 * ncomps + c);
        const auto *gf_1 = buffer.data(gf + 1 * ncomps + c);
        const auto *gf_2 = buffer.data(gf + 2 * ncomps + c);
        const auto *gf_3 = buffer.data(gf + 3 * ncomps + c);
        const auto *gf_4 = buffer.data(gf + 4 * ncomps + c);
        const auto *gf_5 = buffer.data(gf + 5 * ncomps + c);
        const auto *gf_6 = buffer.data(gf + 6 * ncomps + c);
        const auto *gf_7 = buffer.data(gf + 7 * ncomps + c);
        const auto *gf_8 = buffer.data(gf + 8 * ncomps + c);
        const auto *gf_9 = buffer.data(gf + 9 * ncomps + c);
        const auto *gf_10 = buffer.data(gf + 10 * ncomps + c);
        const auto *gf_11 = buffer.data(gf + 11 * ncomps + c);
        const auto *gf_12 = buffer.data(gf + 12 * ncomps + c);
        const auto *gf_13 = buffer.data(gf + 13 * ncomps + c);
        const auto *gf_14 = buffer.data(gf + 14 * ncomps + c);
        const auto *gf_15 = buffer.data(gf + 15 * ncomps + c);
        const auto *gf_16 = buffer.data(gf + 16 * ncomps + c);
        const auto *gf_17 = buffer.data(gf + 17 * ncomps + c);
        const auto *gf_18 = buffer.data(gf + 18 * ncomps + c);
        const auto *gf_19 = buffer.data(gf + 19 * ncomps + c);
        const auto *gf_20 = buffer.data(gf + 20 * ncomps + c);
        const auto *gf_21 = buffer.data(gf + 21 * ncomps + c);
        const auto *gf_22 = buffer.data(gf + 22 * ncomps + c);
        const auto *gf_23 = buffer.data(gf + 23 * ncomps + c);
        const auto *gf_24 = buffer.data(gf + 24 * ncomps + c);
        const auto *gf_25 = buffer.data(gf + 25 * ncomps + c);
        const auto *gf_26 = buffer.data(gf + 26 * ncomps + c);
        const auto *gf_27 = buffer.data(gf + 27 * ncomps + c);
        const auto *gf_28 = buffer.data(gf + 28 * ncomps + c);
        const auto *gf_29 = buffer.data(gf + 29 * ncomps + c);
        const auto *gf_30 = buffer.data(gf + 30 * ncomps + c);
        const auto *gf_31 = buffer.data(gf + 31 * ncomps + c);
        const auto *gf_32 = buffer.data(gf + 32 * ncomps + c);
        const auto *gf_33 = buffer.data(gf + 33 * ncomps + c);
        const auto *gf_34 = buffer.data(gf + 34 * ncomps + c);
        const auto *gf_35 = buffer.data(gf + 35 * ncomps + c);
        const auto *gf_36 = buffer.data(gf + 36 * ncomps + c);
        const auto *gf_37 = buffer.data(gf + 37 * ncomps + c);
        const auto *gf_38 = buffer.data(gf + 38 * ncomps + c);
        const auto *gf_39 = buffer.data(gf + 39 * ncomps + c);
        const auto *gf_40 = buffer.data(gf + 40 * ncomps + c);
        const auto *gf_41 = buffer.data(gf + 41 * ncomps + c);
        const auto *gf_42 = buffer.data(gf + 42 * ncomps + c);
        const auto *gf_43 = buffer.data(gf + 43 * ncomps + c);
        const auto *gf_44 = buffer.data(gf + 44 * ncomps + c);
        const auto *gf_45 = buffer.data(gf + 45 * ncomps + c);
        const auto *gf_46 = buffer.data(gf + 46 * ncomps + c);
        const auto *gf_47 = buffer.data(gf + 47 * ncomps + c);
        const auto *gf_48 = buffer.data(gf + 48 * ncomps + c);
        const auto *gf_49 = buffer.data(gf + 49 * ncomps + c);
        const auto *gf_50 = buffer.data(gf + 50 * ncomps + c);
        const auto *gf_51 = buffer.data(gf + 51 * ncomps + c);
        const auto *gf_52 = buffer.data(gf + 52 * ncomps + c);
        const auto *gf_53 = buffer.data(gf + 53 * ncomps + c);
        const auto *gf_54 = buffer.data(gf + 54 * ncomps + c);
        const auto *gf_55 = buffer.data(gf + 55 * ncomps + c);
        const auto *gf_56 = buffer.data(gf + 56 * ncomps + c);
        const auto *gf_57 = buffer.data(gf + 57 * ncomps + c);
        const auto *gf_58 = buffer.data(gf + 58 * ncomps + c);
        const auto *gf_59 = buffer.data(gf + 59 * ncomps + c);
        const auto *gf_60 = buffer.data(gf + 60 * ncomps + c);
        const auto *gf_61 = buffer.data(gf + 61 * ncomps + c);
        const auto *gf_62 = buffer.data(gf + 62 * ncomps + c);
        const auto *gf_63 = buffer.data(gf + 63 * ncomps + c);
        const auto *gf_64 = buffer.data(gf + 64 * ncomps + c);
        const auto *gf_65 = buffer.data(gf + 65 * ncomps + c);
        const auto *gf_66 = buffer.data(gf + 66 * ncomps + c);
        const auto *gf_67 = buffer.data(gf + 67 * ncomps + c);
        const auto *gf_68 = buffer.data(gf + 68 * ncomps + c);
        const auto *gf_69 = buffer.data(gf + 69 * ncomps + c);
        const auto *gf_70 = buffer.data(gf + 70 * ncomps + c);
        const auto *gf_71 = buffer.data(gf + 71 * ncomps + c);
        const auto *gf_72 = buffer.data(gf + 72 * ncomps + c);
        const auto *gf_73 = buffer.data(gf + 73 * ncomps + c);
        const auto *gf_74 = buffer.data(gf + 74 * ncomps + c);
        const auto *gf_75 = buffer.data(gf + 75 * ncomps + c);
        const auto *gf_76 = buffer.data(gf + 76 * ncomps + c);
        const auto *gf_77 = buffer.data(gf + 77 * ncomps + c);
        const auto *gf_78 = buffer.data(gf + 78 * ncomps + c);
        const auto *gf_79 = buffer.data(gf + 79 * ncomps + c);
        const auto *gf_80 = buffer.data(gf + 80 * ncomps + c);
        const auto *gf_81 = buffer.data(gf + 81 * ncomps + c);
        const auto *gf_82 = buffer.data(gf + 82 * ncomps + c);
        const auto *gf_83 = buffer.data(gf + 83 * ncomps + c);
        const auto *gf_84 = buffer.data(gf + 84 * ncomps + c);
        const auto *gf_85 = buffer.data(gf + 85 * ncomps + c);
        const auto *gf_86 = buffer.data(gf + 86 * ncomps + c);
        const auto *gf_87 = buffer.data(gf + 87 * ncomps + c);
        const auto *gf_88 = buffer.data(gf + 88 * ncomps + c);
        const auto *gf_89 = buffer.data(gf + 89 * ncomps + c);
        const auto *gf_90 = buffer.data(gf + 90 * ncomps + c);
        const auto *gf_91 = buffer.data(gf + 91 * ncomps + c);
        const auto *gf_92 = buffer.data(gf + 92 * ncomps + c);
        const auto *gf_93 = buffer.data(gf + 93 * ncomps + c);
        const auto *gf_94 = buffer.data(gf + 94 * ncomps + c);
        const auto *gf_95 = buffer.data(gf + 95 * ncomps + c);
        const auto *gf_96 = buffer.data(gf + 96 * ncomps + c);
        const auto *gf_97 = buffer.data(gf + 97 * ncomps + c);
        const auto *gf_98 = buffer.data(gf + 98 * ncomps + c);
        const auto *gf_99 = buffer.data(gf + 99 * ncomps + c);

        const auto *hf_0 = buffer.data(hf + 0 * ncomps + c);
        const auto *hf_1 = buffer.data(hf + 1 * ncomps + c);
        const auto *hf_2 = buffer.data(hf + 2 * ncomps + c);
        const auto *hf_3 = buffer.data(hf + 3 * ncomps + c);
        const auto *hf_4 = buffer.data(hf + 4 * ncomps + c);
        const auto *hf_5 = buffer.data(hf + 5 * ncomps + c);
        const auto *hf_6 = buffer.data(hf + 6 * ncomps + c);
        const auto *hf_7 = buffer.data(hf + 7 * ncomps + c);
        const auto *hf_8 = buffer.data(hf + 8 * ncomps + c);
        const auto *hf_9 = buffer.data(hf + 9 * ncomps + c);
        const auto *hf_10 = buffer.data(hf + 10 * ncomps + c);
        const auto *hf_11 = buffer.data(hf + 11 * ncomps + c);
        const auto *hf_12 = buffer.data(hf + 12 * ncomps + c);
        const auto *hf_13 = buffer.data(hf + 13 * ncomps + c);
        const auto *hf_14 = buffer.data(hf + 14 * ncomps + c);
        const auto *hf_15 = buffer.data(hf + 15 * ncomps + c);
        const auto *hf_16 = buffer.data(hf + 16 * ncomps + c);
        const auto *hf_17 = buffer.data(hf + 17 * ncomps + c);
        const auto *hf_18 = buffer.data(hf + 18 * ncomps + c);
        const auto *hf_19 = buffer.data(hf + 19 * ncomps + c);
        const auto *hf_20 = buffer.data(hf + 20 * ncomps + c);
        const auto *hf_21 = buffer.data(hf + 21 * ncomps + c);
        const auto *hf_22 = buffer.data(hf + 22 * ncomps + c);
        const auto *hf_23 = buffer.data(hf + 23 * ncomps + c);
        const auto *hf_24 = buffer.data(hf + 24 * ncomps + c);
        const auto *hf_25 = buffer.data(hf + 25 * ncomps + c);
        const auto *hf_26 = buffer.data(hf + 26 * ncomps + c);
        const auto *hf_27 = buffer.data(hf + 27 * ncomps + c);
        const auto *hf_28 = buffer.data(hf + 28 * ncomps + c);
        const auto *hf_29 = buffer.data(hf + 29 * ncomps + c);
        const auto *hf_30 = buffer.data(hf + 30 * ncomps + c);
        const auto *hf_31 = buffer.data(hf + 31 * ncomps + c);
        const auto *hf_32 = buffer.data(hf + 32 * ncomps + c);
        const auto *hf_33 = buffer.data(hf + 33 * ncomps + c);
        const auto *hf_34 = buffer.data(hf + 34 * ncomps + c);
        const auto *hf_35 = buffer.data(hf + 35 * ncomps + c);
        const auto *hf_36 = buffer.data(hf + 36 * ncomps + c);
        const auto *hf_37 = buffer.data(hf + 37 * ncomps + c);
        const auto *hf_38 = buffer.data(hf + 38 * ncomps + c);
        const auto *hf_39 = buffer.data(hf + 39 * ncomps + c);
        const auto *hf_40 = buffer.data(hf + 40 * ncomps + c);
        const auto *hf_41 = buffer.data(hf + 41 * ncomps + c);
        const auto *hf_42 = buffer.data(hf + 42 * ncomps + c);
        const auto *hf_43 = buffer.data(hf + 43 * ncomps + c);
        const auto *hf_44 = buffer.data(hf + 44 * ncomps + c);
        const auto *hf_45 = buffer.data(hf + 45 * ncomps + c);
        const auto *hf_46 = buffer.data(hf + 46 * ncomps + c);
        const auto *hf_47 = buffer.data(hf + 47 * ncomps + c);
        const auto *hf_48 = buffer.data(hf + 48 * ncomps + c);
        const auto *hf_49 = buffer.data(hf + 49 * ncomps + c);
        const auto *hf_50 = buffer.data(hf + 50 * ncomps + c);
        const auto *hf_51 = buffer.data(hf + 51 * ncomps + c);
        const auto *hf_52 = buffer.data(hf + 52 * ncomps + c);
        const auto *hf_53 = buffer.data(hf + 53 * ncomps + c);
        const auto *hf_54 = buffer.data(hf + 54 * ncomps + c);
        const auto *hf_55 = buffer.data(hf + 55 * ncomps + c);
        const auto *hf_56 = buffer.data(hf + 56 * ncomps + c);
        const auto *hf_57 = buffer.data(hf + 57 * ncomps + c);
        const auto *hf_58 = buffer.data(hf + 58 * ncomps + c);
        const auto *hf_59 = buffer.data(hf + 59 * ncomps + c);
        const auto *hf_60 = buffer.data(hf + 60 * ncomps + c);
        const auto *hf_61 = buffer.data(hf + 61 * ncomps + c);
        const auto *hf_62 = buffer.data(hf + 62 * ncomps + c);
        const auto *hf_63 = buffer.data(hf + 63 * ncomps + c);
        const auto *hf_64 = buffer.data(hf + 64 * ncomps + c);
        const auto *hf_65 = buffer.data(hf + 65 * ncomps + c);
        const auto *hf_66 = buffer.data(hf + 66 * ncomps + c);
        const auto *hf_67 = buffer.data(hf + 67 * ncomps + c);
        const auto *hf_68 = buffer.data(hf + 68 * ncomps + c);
        const auto *hf_69 = buffer.data(hf + 69 * ncomps + c);
        const auto *hf_70 = buffer.data(hf + 70 * ncomps + c);
        const auto *hf_71 = buffer.data(hf + 71 * ncomps + c);
        const auto *hf_72 = buffer.data(hf + 72 * ncomps + c);
        const auto *hf_73 = buffer.data(hf + 73 * ncomps + c);
        const auto *hf_74 = buffer.data(hf + 74 * ncomps + c);
        const auto *hf_75 = buffer.data(hf + 75 * ncomps + c);
        const auto *hf_76 = buffer.data(hf + 76 * ncomps + c);
        const auto *hf_77 = buffer.data(hf + 77 * ncomps + c);
        const auto *hf_78 = buffer.data(hf + 78 * ncomps + c);
        const auto *hf_79 = buffer.data(hf + 79 * ncomps + c);
        const auto *hf_80 = buffer.data(hf + 80 * ncomps + c);
        const auto *hf_81 = buffer.data(hf + 81 * ncomps + c);
        const auto *hf_82 = buffer.data(hf + 82 * ncomps + c);
        const auto *hf_83 = buffer.data(hf + 83 * ncomps + c);
        const auto *hf_84 = buffer.data(hf + 84 * ncomps + c);
        const auto *hf_85 = buffer.data(hf + 85 * ncomps + c);
        const auto *hf_86 = buffer.data(hf + 86 * ncomps + c);
        const auto *hf_87 = buffer.data(hf + 87 * ncomps + c);
        const auto *hf_88 = buffer.data(hf + 88 * ncomps + c);
        const auto *hf_89 = buffer.data(hf + 89 * ncomps + c);
        const auto *hf_90 = buffer.data(hf + 90 * ncomps + c);
        const auto *hf_91 = buffer.data(hf + 91 * ncomps + c);
        const auto *hf_92 = buffer.data(hf + 92 * ncomps + c);
        const auto *hf_93 = buffer.data(hf + 93 * ncomps + c);
        const auto *hf_94 = buffer.data(hf + 94 * ncomps + c);
        const auto *hf_95 = buffer.data(hf + 95 * ncomps + c);
        const auto *hf_96 = buffer.data(hf + 96 * ncomps + c);
        const auto *hf_97 = buffer.data(hf + 97 * ncomps + c);
        const auto *hf_98 = buffer.data(hf + 98 * ncomps + c);
        const auto *hf_99 = buffer.data(hf + 99 * ncomps + c);
        const auto *hf_106 = buffer.data(hf + 106 * ncomps + c);
        const auto *hf_107 = buffer.data(hf + 107 * ncomps + c);
        const auto *hf_108 = buffer.data(hf + 108 * ncomps + c);
        const auto *hf_109 = buffer.data(hf + 109 * ncomps + c);
        const auto *hf_116 = buffer.data(hf + 116 * ncomps + c);
        const auto *hf_117 = buffer.data(hf + 117 * ncomps + c);
        const auto *hf_118 = buffer.data(hf + 118 * ncomps + c);
        const auto *hf_119 = buffer.data(hf + 119 * ncomps + c);
        const auto *hf_126 = buffer.data(hf + 126 * ncomps + c);
        const auto *hf_127 = buffer.data(hf + 127 * ncomps + c);
        const auto *hf_128 = buffer.data(hf + 128 * ncomps + c);
        const auto *hf_129 = buffer.data(hf + 129 * ncomps + c);
        const auto *hf_139 = buffer.data(hf + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gf_0, gf_1, gf_2, gf_3, gf_4, hf_0, \
                         hf_1, hf_2, hf_3, hf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gf_0[k]
                     + hf_0[k];

            t_1[k] = ab_x[k] * gf_1[k]
                     + hf_1[k];

            t_2[k] = ab_x[k] * gf_2[k]
                     + hf_2[k];

            t_3[k] = ab_x[k] * gf_3[k]
                     + hf_3[k];

            t_4[k] = ab_x[k] * gf_4[k]
                     + hf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gf_5, gf_6, gf_7, gf_8, gf_9, hf_5, \
                         hf_6, hf_7, hf_8, hf_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * gf_5[k]
                     + hf_5[k];

            t_6[k] = ab_x[k] * gf_6[k]
                     + hf_6[k];

            t_7[k] = ab_x[k] * gf_7[k]
                     + hf_7[k];

            t_8[k] = ab_x[k] * gf_8[k]
                     + hf_8[k];

            t_9[k] = ab_x[k] * gf_9[k]
                     + hf_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, gf_6, gf_7, gf_8, gf_9, \
                         hf_16, hf_17, hf_18, hf_19, hf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * gf_6[k]
                      + hf_16[k];

            t_11[k] = ab_y[k] * gf_7[k]
                      + hf_17[k];

            t_12[k] = ab_y[k] * gf_8[k]
                      + hf_18[k];

            t_13[k] = ab_y[k] * gf_9[k]
                      + hf_19[k];

            t_14[k] = ab_z[k] * gf_9[k]
                      + hf_29[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, gf_10, gf_11, gf_12, gf_13, \
                         gf_14, hf_10, hf_11, hf_12, hf_13, hf_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * gf_10[k]
                      + hf_10[k];

            t_16[k] = ab_x[k] * gf_11[k]
                      + hf_11[k];

            t_17[k] = ab_x[k] * gf_12[k]
                      + hf_12[k];

            t_18[k] = ab_x[k] * gf_13[k]
                      + hf_13[k];

            t_19[k] = ab_x[k] * gf_14[k]
                      + hf_14[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gf_15, gf_16, gf_17, gf_18, \
                         gf_19, hf_15, hf_16, hf_17, hf_18, hf_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * gf_15[k]
                      + hf_15[k];

            t_21[k] = ab_x[k] * gf_16[k]
                      + hf_16[k];

            t_22[k] = ab_x[k] * gf_17[k]
                      + hf_17[k];

            t_23[k] = ab_x[k] * gf_18[k]
                      + hf_18[k];

            t_24[k] = ab_x[k] * gf_19[k]
                      + hf_19[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, gf_16, gf_17, gf_18, gf_19, \
                         hf_36, hf_37, hf_38, hf_39, hf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_y[k] * gf_16[k]
                      + hf_36[k];

            t_26[k] = ab_y[k] * gf_17[k]
                      + hf_37[k];

            t_27[k] = ab_y[k] * gf_18[k]
                      + hf_38[k];

            t_28[k] = ab_y[k] * gf_19[k]
                      + hf_39[k];

            t_29[k] = ab_z[k] * gf_19[k]
                      + hf_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gf_20, gf_21, gf_22, gf_23, \
                         gf_24, hf_20, hf_21, hf_22, hf_23, hf_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * gf_20[k]
                      + hf_20[k];

            t_31[k] = ab_x[k] * gf_21[k]
                      + hf_21[k];

            t_32[k] = ab_x[k] * gf_22[k]
                      + hf_22[k];

            t_33[k] = ab_x[k] * gf_23[k]
                      + hf_23[k];

            t_34[k] = ab_x[k] * gf_24[k]
                      + hf_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gf_25, gf_26, gf_27, gf_28, \
                         gf_29, hf_25, hf_26, hf_27, hf_28, hf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * gf_25[k]
                      + hf_25[k];

            t_36[k] = ab_x[k] * gf_26[k]
                      + hf_26[k];

            t_37[k] = ab_x[k] * gf_27[k]
                      + hf_27[k];

            t_38[k] = ab_x[k] * gf_28[k]
                      + hf_28[k];

            t_39[k] = ab_x[k] * gf_29[k]
                      + hf_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, gf_26, gf_27, gf_28, gf_29, \
                         hf_46, hf_47, hf_48, hf_49, hf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * gf_26[k]
                      + hf_46[k];

            t_41[k] = ab_y[k] * gf_27[k]
                      + hf_47[k];

            t_42[k] = ab_y[k] * gf_28[k]
                      + hf_48[k];

            t_43[k] = ab_y[k] * gf_29[k]
                      + hf_49[k];

            t_44[k] = ab_z[k] * gf_29[k]
                      + hf_59[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, gf_30, gf_31, gf_32, gf_33, \
                         gf_34, hf_30, hf_31, hf_32, hf_33, hf_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * gf_30[k]
                      + hf_30[k];

            t_46[k] = ab_x[k] * gf_31[k]
                      + hf_31[k];

            t_47[k] = ab_x[k] * gf_32[k]
                      + hf_32[k];

            t_48[k] = ab_x[k] * gf_33[k]
                      + hf_33[k];

            t_49[k] = ab_x[k] * gf_34[k]
                      + hf_34[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gf_35, gf_36, gf_37, gf_38, \
                         gf_39, hf_35, hf_36, hf_37, hf_38, hf_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * gf_35[k]
                      + hf_35[k];

            t_51[k] = ab_x[k] * gf_36[k]
                      + hf_36[k];

            t_52[k] = ab_x[k] * gf_37[k]
                      + hf_37[k];

            t_53[k] = ab_x[k] * gf_38[k]
                      + hf_38[k];

            t_54[k] = ab_x[k] * gf_39[k]
                      + hf_39[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, gf_36, gf_37, gf_38, gf_39, \
                         hf_66, hf_67, hf_68, hf_69, hf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_y[k] * gf_36[k]
                      + hf_66[k];

            t_56[k] = ab_y[k] * gf_37[k]
                      + hf_67[k];

            t_57[k] = ab_y[k] * gf_38[k]
                      + hf_68[k];

            t_58[k] = ab_y[k] * gf_39[k]
                      + hf_69[k];

            t_59[k] = ab_z[k] * gf_39[k]
                      + hf_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gf_40, gf_41, gf_42, gf_43, \
                         gf_44, hf_40, hf_41, hf_42, hf_43, hf_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * gf_40[k]
                      + hf_40[k];

            t_61[k] = ab_x[k] * gf_41[k]
                      + hf_41[k];

            t_62[k] = ab_x[k] * gf_42[k]
                      + hf_42[k];

            t_63[k] = ab_x[k] * gf_43[k]
                      + hf_43[k];

            t_64[k] = ab_x[k] * gf_44[k]
                      + hf_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gf_45, gf_46, gf_47, gf_48, \
                         gf_49, hf_45, hf_46, hf_47, hf_48, hf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * gf_45[k]
                      + hf_45[k];

            t_66[k] = ab_x[k] * gf_46[k]
                      + hf_46[k];

            t_67[k] = ab_x[k] * gf_47[k]
                      + hf_47[k];

            t_68[k] = ab_x[k] * gf_48[k]
                      + hf_48[k];

            t_69[k] = ab_x[k] * gf_49[k]
                      + hf_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, gf_46, gf_47, gf_48, gf_49, \
                         hf_76, hf_77, hf_78, hf_79, hf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * gf_46[k]
                      + hf_76[k];

            t_71[k] = ab_y[k] * gf_47[k]
                      + hf_77[k];

            t_72[k] = ab_y[k] * gf_48[k]
                      + hf_78[k];

            t_73[k] = ab_y[k] * gf_49[k]
                      + hf_79[k];

            t_74[k] = ab_z[k] * gf_49[k]
                      + hf_89[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, gf_50, gf_51, gf_52, gf_53, \
                         gf_54, hf_50, hf_51, hf_52, hf_53, hf_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * gf_50[k]
                      + hf_50[k];

            t_76[k] = ab_x[k] * gf_51[k]
                      + hf_51[k];

            t_77[k] = ab_x[k] * gf_52[k]
                      + hf_52[k];

            t_78[k] = ab_x[k] * gf_53[k]
                      + hf_53[k];

            t_79[k] = ab_x[k] * gf_54[k]
                      + hf_54[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gf_55, gf_56, gf_57, gf_58, \
                         gf_59, hf_55, hf_56, hf_57, hf_58, hf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * gf_55[k]
                      + hf_55[k];

            t_81[k] = ab_x[k] * gf_56[k]
                      + hf_56[k];

            t_82[k] = ab_x[k] * gf_57[k]
                      + hf_57[k];

            t_83[k] = ab_x[k] * gf_58[k]
                      + hf_58[k];

            t_84[k] = ab_x[k] * gf_59[k]
                      + hf_59[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, gf_56, gf_57, gf_58, gf_59, \
                         hf_86, hf_87, hf_88, hf_89, hf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_y[k] * gf_56[k]
                      + hf_86[k];

            t_86[k] = ab_y[k] * gf_57[k]
                      + hf_87[k];

            t_87[k] = ab_y[k] * gf_58[k]
                      + hf_88[k];

            t_88[k] = ab_y[k] * gf_59[k]
                      + hf_89[k];

            t_89[k] = ab_z[k] * gf_59[k]
                      + hf_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gf_60, gf_61, gf_62, gf_63, \
                         gf_64, hf_60, hf_61, hf_62, hf_63, hf_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * gf_60[k]
                      + hf_60[k];

            t_91[k] = ab_x[k] * gf_61[k]
                      + hf_61[k];

            t_92[k] = ab_x[k] * gf_62[k]
                      + hf_62[k];

            t_93[k] = ab_x[k] * gf_63[k]
                      + hf_63[k];

            t_94[k] = ab_x[k] * gf_64[k]
                      + hf_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gf_65, gf_66, gf_67, gf_68, \
                         gf_69, hf_65, hf_66, hf_67, hf_68, hf_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * gf_65[k]
                      + hf_65[k];

            t_96[k] = ab_x[k] * gf_66[k]
                      + hf_66[k];

            t_97[k] = ab_x[k] * gf_67[k]
                      + hf_67[k];

            t_98[k] = ab_x[k] * gf_68[k]
                      + hf_68[k];

            t_99[k] = ab_x[k] * gf_69[k]
                      + hf_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, gf_66, gf_67, gf_68, \
                         gf_69, hf_106, hf_107, hf_108, hf_109, \
                         hf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * gf_66[k]
                       + hf_106[k];

            t_101[k] = ab_y[k] * gf_67[k]
                       + hf_107[k];

            t_102[k] = ab_y[k] * gf_68[k]
                       + hf_108[k];

            t_103[k] = ab_y[k] * gf_69[k]
                       + hf_109[k];

            t_104[k] = ab_z[k] * gf_69[k]
                       + hf_119[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, gf_70, gf_71, gf_72, gf_73, \
                         gf_74, hf_70, hf_71, hf_72, hf_73, hf_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * gf_70[k]
                       + hf_70[k];

            t_106[k] = ab_x[k] * gf_71[k]
                       + hf_71[k];

            t_107[k] = ab_x[k] * gf_72[k]
                       + hf_72[k];

            t_108[k] = ab_x[k] * gf_73[k]
                       + hf_73[k];

            t_109[k] = ab_x[k] * gf_74[k]
                       + hf_74[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gf_75, gf_76, gf_77, gf_78, \
                         gf_79, hf_75, hf_76, hf_77, hf_78, hf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * gf_75[k]
                       + hf_75[k];

            t_111[k] = ab_x[k] * gf_76[k]
                       + hf_76[k];

            t_112[k] = ab_x[k] * gf_77[k]
                       + hf_77[k];

            t_113[k] = ab_x[k] * gf_78[k]
                       + hf_78[k];

            t_114[k] = ab_x[k] * gf_79[k]
                       + hf_79[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, gf_76, gf_77, gf_78, \
                         gf_79, hf_116, hf_117, hf_118, hf_119, \
                         hf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_y[k] * gf_76[k]
                       + hf_116[k];

            t_116[k] = ab_y[k] * gf_77[k]
                       + hf_117[k];

            t_117[k] = ab_y[k] * gf_78[k]
                       + hf_118[k];

            t_118[k] = ab_y[k] * gf_79[k]
                       + hf_119[k];

            t_119[k] = ab_z[k] * gf_79[k]
                       + hf_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gf_80, gf_81, gf_82, gf_83, \
                         gf_84, hf_80, hf_81, hf_82, hf_83, hf_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * gf_80[k]
                       + hf_80[k];

            t_121[k] = ab_x[k] * gf_81[k]
                       + hf_81[k];

            t_122[k] = ab_x[k] * gf_82[k]
                       + hf_82[k];

            t_123[k] = ab_x[k] * gf_83[k]
                       + hf_83[k];

            t_124[k] = ab_x[k] * gf_84[k]
                       + hf_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gf_85, gf_86, gf_87, gf_88, \
                         gf_89, hf_85, hf_86, hf_87, hf_88, hf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * gf_85[k]
                       + hf_85[k];

            t_126[k] = ab_x[k] * gf_86[k]
                       + hf_86[k];

            t_127[k] = ab_x[k] * gf_87[k]
                       + hf_87[k];

            t_128[k] = ab_x[k] * gf_88[k]
                       + hf_88[k];

            t_129[k] = ab_x[k] * gf_89[k]
                       + hf_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, gf_86, gf_87, gf_88, \
                         gf_89, hf_126, hf_127, hf_128, hf_129, \
                         hf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * gf_86[k]
                       + hf_126[k];

            t_131[k] = ab_y[k] * gf_87[k]
                       + hf_127[k];

            t_132[k] = ab_y[k] * gf_88[k]
                       + hf_128[k];

            t_133[k] = ab_y[k] * gf_89[k]
                       + hf_129[k];

            t_134[k] = ab_z[k] * gf_89[k]
                       + hf_139[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, gf_90, gf_91, gf_92, gf_93, \
                         gf_94, hf_90, hf_91, hf_92, hf_93, hf_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * gf_90[k]
                       + hf_90[k];

            t_136[k] = ab_x[k] * gf_91[k]
                       + hf_91[k];

            t_137[k] = ab_x[k] * gf_92[k]
                       + hf_92[k];

            t_138[k] = ab_x[k] * gf_93[k]
                       + hf_93[k];

            t_139[k] = ab_x[k] * gf_94[k]
                       + hf_94[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gf_95, gf_96, gf_97, gf_98, \
                         gf_99, hf_95, hf_96, hf_97, hf_98, hf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * gf_95[k]
                       + hf_95[k];

            t_141[k] = ab_x[k] * gf_96[k]
                       + hf_96[k];

            t_142[k] = ab_x[k] * gf_97[k]
                       + hf_97[k];

            t_143[k] = ab_x[k] * gf_98[k]
                       + hf_98[k];

            t_144[k] = ab_x[k] * gf_99[k]
                       + hf_99[k];
        }
    }
}

static auto
compute_hrr_gg_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t gf, const size_t hf,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gf_96 = buffer.data(gf + 96 * ncomps + c);
        const auto *gf_97 = buffer.data(gf + 97 * ncomps + c);
        const auto *gf_98 = buffer.data(gf + 98 * ncomps + c);
        const auto *gf_99 = buffer.data(gf + 99 * ncomps + c);
        const auto *gf_100 = buffer.data(gf + 100 * ncomps + c);
        const auto *gf_101 = buffer.data(gf + 101 * ncomps + c);
        const auto *gf_102 = buffer.data(gf + 102 * ncomps + c);
        const auto *gf_103 = buffer.data(gf + 103 * ncomps + c);
        const auto *gf_104 = buffer.data(gf + 104 * ncomps + c);
        const auto *gf_105 = buffer.data(gf + 105 * ncomps + c);
        const auto *gf_106 = buffer.data(gf + 106 * ncomps + c);
        const auto *gf_107 = buffer.data(gf + 107 * ncomps + c);
        const auto *gf_108 = buffer.data(gf + 108 * ncomps + c);
        const auto *gf_109 = buffer.data(gf + 109 * ncomps + c);
        const auto *gf_110 = buffer.data(gf + 110 * ncomps + c);
        const auto *gf_111 = buffer.data(gf + 111 * ncomps + c);
        const auto *gf_112 = buffer.data(gf + 112 * ncomps + c);
        const auto *gf_113 = buffer.data(gf + 113 * ncomps + c);
        const auto *gf_114 = buffer.data(gf + 114 * ncomps + c);
        const auto *gf_115 = buffer.data(gf + 115 * ncomps + c);
        const auto *gf_116 = buffer.data(gf + 116 * ncomps + c);
        const auto *gf_117 = buffer.data(gf + 117 * ncomps + c);
        const auto *gf_118 = buffer.data(gf + 118 * ncomps + c);
        const auto *gf_119 = buffer.data(gf + 119 * ncomps + c);
        const auto *gf_120 = buffer.data(gf + 120 * ncomps + c);
        const auto *gf_121 = buffer.data(gf + 121 * ncomps + c);
        const auto *gf_122 = buffer.data(gf + 122 * ncomps + c);
        const auto *gf_123 = buffer.data(gf + 123 * ncomps + c);
        const auto *gf_124 = buffer.data(gf + 124 * ncomps + c);
        const auto *gf_125 = buffer.data(gf + 125 * ncomps + c);
        const auto *gf_126 = buffer.data(gf + 126 * ncomps + c);
        const auto *gf_127 = buffer.data(gf + 127 * ncomps + c);
        const auto *gf_128 = buffer.data(gf + 128 * ncomps + c);
        const auto *gf_129 = buffer.data(gf + 129 * ncomps + c);
        const auto *gf_130 = buffer.data(gf + 130 * ncomps + c);
        const auto *gf_131 = buffer.data(gf + 131 * ncomps + c);
        const auto *gf_132 = buffer.data(gf + 132 * ncomps + c);
        const auto *gf_133 = buffer.data(gf + 133 * ncomps + c);
        const auto *gf_134 = buffer.data(gf + 134 * ncomps + c);
        const auto *gf_135 = buffer.data(gf + 135 * ncomps + c);
        const auto *gf_136 = buffer.data(gf + 136 * ncomps + c);
        const auto *gf_137 = buffer.data(gf + 137 * ncomps + c);
        const auto *gf_138 = buffer.data(gf + 138 * ncomps + c);
        const auto *gf_139 = buffer.data(gf + 139 * ncomps + c);
        const auto *gf_140 = buffer.data(gf + 140 * ncomps + c);
        const auto *gf_141 = buffer.data(gf + 141 * ncomps + c);
        const auto *gf_142 = buffer.data(gf + 142 * ncomps + c);
        const auto *gf_143 = buffer.data(gf + 143 * ncomps + c);
        const auto *gf_144 = buffer.data(gf + 144 * ncomps + c);
        const auto *gf_145 = buffer.data(gf + 145 * ncomps + c);
        const auto *gf_146 = buffer.data(gf + 146 * ncomps + c);
        const auto *gf_147 = buffer.data(gf + 147 * ncomps + c);
        const auto *gf_148 = buffer.data(gf + 148 * ncomps + c);
        const auto *gf_149 = buffer.data(gf + 149 * ncomps + c);

        const auto *hf_100 = buffer.data(hf + 100 * ncomps + c);
        const auto *hf_101 = buffer.data(hf + 101 * ncomps + c);
        const auto *hf_102 = buffer.data(hf + 102 * ncomps + c);
        const auto *hf_103 = buffer.data(hf + 103 * ncomps + c);
        const auto *hf_104 = buffer.data(hf + 104 * ncomps + c);
        const auto *hf_105 = buffer.data(hf + 105 * ncomps + c);
        const auto *hf_106 = buffer.data(hf + 106 * ncomps + c);
        const auto *hf_107 = buffer.data(hf + 107 * ncomps + c);
        const auto *hf_108 = buffer.data(hf + 108 * ncomps + c);
        const auto *hf_109 = buffer.data(hf + 109 * ncomps + c);
        const auto *hf_110 = buffer.data(hf + 110 * ncomps + c);
        const auto *hf_111 = buffer.data(hf + 111 * ncomps + c);
        const auto *hf_112 = buffer.data(hf + 112 * ncomps + c);
        const auto *hf_113 = buffer.data(hf + 113 * ncomps + c);
        const auto *hf_114 = buffer.data(hf + 114 * ncomps + c);
        const auto *hf_115 = buffer.data(hf + 115 * ncomps + c);
        const auto *hf_116 = buffer.data(hf + 116 * ncomps + c);
        const auto *hf_117 = buffer.data(hf + 117 * ncomps + c);
        const auto *hf_118 = buffer.data(hf + 118 * ncomps + c);
        const auto *hf_119 = buffer.data(hf + 119 * ncomps + c);
        const auto *hf_120 = buffer.data(hf + 120 * ncomps + c);
        const auto *hf_121 = buffer.data(hf + 121 * ncomps + c);
        const auto *hf_122 = buffer.data(hf + 122 * ncomps + c);
        const auto *hf_123 = buffer.data(hf + 123 * ncomps + c);
        const auto *hf_124 = buffer.data(hf + 124 * ncomps + c);
        const auto *hf_125 = buffer.data(hf + 125 * ncomps + c);
        const auto *hf_126 = buffer.data(hf + 126 * ncomps + c);
        const auto *hf_127 = buffer.data(hf + 127 * ncomps + c);
        const auto *hf_128 = buffer.data(hf + 128 * ncomps + c);
        const auto *hf_129 = buffer.data(hf + 129 * ncomps + c);
        const auto *hf_130 = buffer.data(hf + 130 * ncomps + c);
        const auto *hf_131 = buffer.data(hf + 131 * ncomps + c);
        const auto *hf_132 = buffer.data(hf + 132 * ncomps + c);
        const auto *hf_133 = buffer.data(hf + 133 * ncomps + c);
        const auto *hf_134 = buffer.data(hf + 134 * ncomps + c);
        const auto *hf_135 = buffer.data(hf + 135 * ncomps + c);
        const auto *hf_136 = buffer.data(hf + 136 * ncomps + c);
        const auto *hf_137 = buffer.data(hf + 137 * ncomps + c);
        const auto *hf_138 = buffer.data(hf + 138 * ncomps + c);
        const auto *hf_139 = buffer.data(hf + 139 * ncomps + c);
        const auto *hf_140 = buffer.data(hf + 140 * ncomps + c);
        const auto *hf_141 = buffer.data(hf + 141 * ncomps + c);
        const auto *hf_142 = buffer.data(hf + 142 * ncomps + c);
        const auto *hf_143 = buffer.data(hf + 143 * ncomps + c);
        const auto *hf_144 = buffer.data(hf + 144 * ncomps + c);
        const auto *hf_145 = buffer.data(hf + 145 * ncomps + c);
        const auto *hf_146 = buffer.data(hf + 146 * ncomps + c);
        const auto *hf_147 = buffer.data(hf + 147 * ncomps + c);
        const auto *hf_148 = buffer.data(hf + 148 * ncomps + c);
        const auto *hf_149 = buffer.data(hf + 149 * ncomps + c);
        const auto *hf_156 = buffer.data(hf + 156 * ncomps + c);
        const auto *hf_157 = buffer.data(hf + 157 * ncomps + c);
        const auto *hf_158 = buffer.data(hf + 158 * ncomps + c);
        const auto *hf_159 = buffer.data(hf + 159 * ncomps + c);
        const auto *hf_166 = buffer.data(hf + 166 * ncomps + c);
        const auto *hf_167 = buffer.data(hf + 167 * ncomps + c);
        const auto *hf_168 = buffer.data(hf + 168 * ncomps + c);
        const auto *hf_169 = buffer.data(hf + 169 * ncomps + c);
        const auto *hf_176 = buffer.data(hf + 176 * ncomps + c);
        const auto *hf_177 = buffer.data(hf + 177 * ncomps + c);
        const auto *hf_178 = buffer.data(hf + 178 * ncomps + c);
        const auto *hf_179 = buffer.data(hf + 179 * ncomps + c);
        const auto *hf_186 = buffer.data(hf + 186 * ncomps + c);
        const auto *hf_187 = buffer.data(hf + 187 * ncomps + c);
        const auto *hf_188 = buffer.data(hf + 188 * ncomps + c);
        const auto *hf_189 = buffer.data(hf + 189 * ncomps + c);
        const auto *hf_196 = buffer.data(hf + 196 * ncomps + c);
        const auto *hf_197 = buffer.data(hf + 197 * ncomps + c);
        const auto *hf_198 = buffer.data(hf + 198 * ncomps + c);
        const auto *hf_199 = buffer.data(hf + 199 * ncomps + c);
        const auto *hf_209 = buffer.data(hf + 209 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, gf_96, gf_97, gf_98, \
                         gf_99, hf_136, hf_137, hf_138, hf_139, \
                         hf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_y[k] * gf_96[k]
                       + hf_136[k];

            t_146[k] = ab_y[k] * gf_97[k]
                       + hf_137[k];

            t_147[k] = ab_y[k] * gf_98[k]
                       + hf_138[k];

            t_148[k] = ab_y[k] * gf_99[k]
                       + hf_139[k];

            t_149[k] = ab_z[k] * gf_99[k]
                       + hf_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gf_100, gf_101, gf_102, \
                         gf_103, gf_104, hf_100, hf_101, hf_102, hf_103, \
                         hf_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * gf_100[k]
                       + hf_100[k];

            t_151[k] = ab_x[k] * gf_101[k]
                       + hf_101[k];

            t_152[k] = ab_x[k] * gf_102[k]
                       + hf_102[k];

            t_153[k] = ab_x[k] * gf_103[k]
                       + hf_103[k];

            t_154[k] = ab_x[k] * gf_104[k]
                       + hf_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gf_105, gf_106, gf_107, \
                         gf_108, gf_109, hf_105, hf_106, hf_107, hf_108, \
                         hf_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * gf_105[k]
                       + hf_105[k];

            t_156[k] = ab_x[k] * gf_106[k]
                       + hf_106[k];

            t_157[k] = ab_x[k] * gf_107[k]
                       + hf_107[k];

            t_158[k] = ab_x[k] * gf_108[k]
                       + hf_108[k];

            t_159[k] = ab_x[k] * gf_109[k]
                       + hf_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, gf_106, gf_107, \
                         gf_108, gf_109, hf_156, hf_157, hf_158, hf_159, \
                         hf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * gf_106[k]
                       + hf_156[k];

            t_161[k] = ab_y[k] * gf_107[k]
                       + hf_157[k];

            t_162[k] = ab_y[k] * gf_108[k]
                       + hf_158[k];

            t_163[k] = ab_y[k] * gf_109[k]
                       + hf_159[k];

            t_164[k] = ab_z[k] * gf_109[k]
                       + hf_169[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, gf_110, gf_111, gf_112, \
                         gf_113, gf_114, hf_110, hf_111, hf_112, hf_113, \
                         hf_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * gf_110[k]
                       + hf_110[k];

            t_166[k] = ab_x[k] * gf_111[k]
                       + hf_111[k];

            t_167[k] = ab_x[k] * gf_112[k]
                       + hf_112[k];

            t_168[k] = ab_x[k] * gf_113[k]
                       + hf_113[k];

            t_169[k] = ab_x[k] * gf_114[k]
                       + hf_114[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, gf_115, gf_116, gf_117, \
                         gf_118, gf_119, hf_115, hf_116, hf_117, hf_118, \
                         hf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * gf_115[k]
                       + hf_115[k];

            t_171[k] = ab_x[k] * gf_116[k]
                       + hf_116[k];

            t_172[k] = ab_x[k] * gf_117[k]
                       + hf_117[k];

            t_173[k] = ab_x[k] * gf_118[k]
                       + hf_118[k];

            t_174[k] = ab_x[k] * gf_119[k]
                       + hf_119[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, gf_116, gf_117, \
                         gf_118, gf_119, hf_166, hf_167, hf_168, hf_169, \
                         hf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_y[k] * gf_116[k]
                       + hf_166[k];

            t_176[k] = ab_y[k] * gf_117[k]
                       + hf_167[k];

            t_177[k] = ab_y[k] * gf_118[k]
                       + hf_168[k];

            t_178[k] = ab_y[k] * gf_119[k]
                       + hf_169[k];

            t_179[k] = ab_z[k] * gf_119[k]
                       + hf_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gf_120, gf_121, gf_122, \
                         gf_123, gf_124, hf_120, hf_121, hf_122, hf_123, \
                         hf_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * gf_120[k]
                       + hf_120[k];

            t_181[k] = ab_x[k] * gf_121[k]
                       + hf_121[k];

            t_182[k] = ab_x[k] * gf_122[k]
                       + hf_122[k];

            t_183[k] = ab_x[k] * gf_123[k]
                       + hf_123[k];

            t_184[k] = ab_x[k] * gf_124[k]
                       + hf_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gf_125, gf_126, gf_127, \
                         gf_128, gf_129, hf_125, hf_126, hf_127, hf_128, \
                         hf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * gf_125[k]
                       + hf_125[k];

            t_186[k] = ab_x[k] * gf_126[k]
                       + hf_126[k];

            t_187[k] = ab_x[k] * gf_127[k]
                       + hf_127[k];

            t_188[k] = ab_x[k] * gf_128[k]
                       + hf_128[k];

            t_189[k] = ab_x[k] * gf_129[k]
                       + hf_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, gf_126, gf_127, \
                         gf_128, gf_129, hf_176, hf_177, hf_178, hf_179, \
                         hf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * gf_126[k]
                       + hf_176[k];

            t_191[k] = ab_y[k] * gf_127[k]
                       + hf_177[k];

            t_192[k] = ab_y[k] * gf_128[k]
                       + hf_178[k];

            t_193[k] = ab_y[k] * gf_129[k]
                       + hf_179[k];

            t_194[k] = ab_z[k] * gf_129[k]
                       + hf_189[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, gf_130, gf_131, gf_132, \
                         gf_133, gf_134, hf_130, hf_131, hf_132, hf_133, \
                         hf_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * gf_130[k]
                       + hf_130[k];

            t_196[k] = ab_x[k] * gf_131[k]
                       + hf_131[k];

            t_197[k] = ab_x[k] * gf_132[k]
                       + hf_132[k];

            t_198[k] = ab_x[k] * gf_133[k]
                       + hf_133[k];

            t_199[k] = ab_x[k] * gf_134[k]
                       + hf_134[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, gf_135, gf_136, gf_137, \
                         gf_138, gf_139, hf_135, hf_136, hf_137, hf_138, \
                         hf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * gf_135[k]
                       + hf_135[k];

            t_201[k] = ab_x[k] * gf_136[k]
                       + hf_136[k];

            t_202[k] = ab_x[k] * gf_137[k]
                       + hf_137[k];

            t_203[k] = ab_x[k] * gf_138[k]
                       + hf_138[k];

            t_204[k] = ab_x[k] * gf_139[k]
                       + hf_139[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, gf_136, gf_137, \
                         gf_138, gf_139, hf_186, hf_187, hf_188, hf_189, \
                         hf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_y[k] * gf_136[k]
                       + hf_186[k];

            t_206[k] = ab_y[k] * gf_137[k]
                       + hf_187[k];

            t_207[k] = ab_y[k] * gf_138[k]
                       + hf_188[k];

            t_208[k] = ab_y[k] * gf_139[k]
                       + hf_189[k];

            t_209[k] = ab_z[k] * gf_139[k]
                       + hf_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gf_140, gf_141, gf_142, \
                         gf_143, gf_144, hf_140, hf_141, hf_142, hf_143, \
                         hf_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * gf_140[k]
                       + hf_140[k];

            t_211[k] = ab_x[k] * gf_141[k]
                       + hf_141[k];

            t_212[k] = ab_x[k] * gf_142[k]
                       + hf_142[k];

            t_213[k] = ab_x[k] * gf_143[k]
                       + hf_143[k];

            t_214[k] = ab_x[k] * gf_144[k]
                       + hf_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gf_145, gf_146, gf_147, \
                         gf_148, gf_149, hf_145, hf_146, hf_147, hf_148, \
                         hf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * gf_145[k]
                       + hf_145[k];

            t_216[k] = ab_x[k] * gf_146[k]
                       + hf_146[k];

            t_217[k] = ab_x[k] * gf_147[k]
                       + hf_147[k];

            t_218[k] = ab_x[k] * gf_148[k]
                       + hf_148[k];

            t_219[k] = ab_x[k] * gf_149[k]
                       + hf_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, gf_146, gf_147, \
                         gf_148, gf_149, hf_196, hf_197, hf_198, hf_199, \
                         hf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * gf_146[k]
                       + hf_196[k];

            t_221[k] = ab_y[k] * gf_147[k]
                       + hf_197[k];

            t_222[k] = ab_y[k] * gf_148[k]
                       + hf_198[k];

            t_223[k] = ab_y[k] * gf_149[k]
                       + hf_199[k];

            t_224[k] = ab_z[k] * gf_149[k]
                       + hf_209[k];
        }
    }
}

auto
compute_hrr_gg_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t gf, const size_t hf,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_gg_out_of_first_piece0(buffer, coordinates, target, gf, hf, ncomps, nmax);

    compute_hrr_gg_out_of_first_piece1(buffer, coordinates, target, gf, hf, ncomps, nmax);
}

static auto
compute_hrr_gg_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fg, const size_t fh, const size_t ncomps,
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

        const auto *fg_0 = buffer.data(fg + 0 * ncomps + c);
        const auto *fg_1 = buffer.data(fg + 1 * ncomps + c);
        const auto *fg_2 = buffer.data(fg + 2 * ncomps + c);
        const auto *fg_3 = buffer.data(fg + 3 * ncomps + c);
        const auto *fg_4 = buffer.data(fg + 4 * ncomps + c);
        const auto *fg_5 = buffer.data(fg + 5 * ncomps + c);
        const auto *fg_6 = buffer.data(fg + 6 * ncomps + c);
        const auto *fg_7 = buffer.data(fg + 7 * ncomps + c);
        const auto *fg_8 = buffer.data(fg + 8 * ncomps + c);
        const auto *fg_9 = buffer.data(fg + 9 * ncomps + c);
        const auto *fg_10 = buffer.data(fg + 10 * ncomps + c);
        const auto *fg_11 = buffer.data(fg + 11 * ncomps + c);
        const auto *fg_12 = buffer.data(fg + 12 * ncomps + c);
        const auto *fg_13 = buffer.data(fg + 13 * ncomps + c);
        const auto *fg_14 = buffer.data(fg + 14 * ncomps + c);
        const auto *fg_15 = buffer.data(fg + 15 * ncomps + c);
        const auto *fg_16 = buffer.data(fg + 16 * ncomps + c);
        const auto *fg_17 = buffer.data(fg + 17 * ncomps + c);
        const auto *fg_18 = buffer.data(fg + 18 * ncomps + c);
        const auto *fg_19 = buffer.data(fg + 19 * ncomps + c);
        const auto *fg_20 = buffer.data(fg + 20 * ncomps + c);
        const auto *fg_21 = buffer.data(fg + 21 * ncomps + c);
        const auto *fg_22 = buffer.data(fg + 22 * ncomps + c);
        const auto *fg_23 = buffer.data(fg + 23 * ncomps + c);
        const auto *fg_24 = buffer.data(fg + 24 * ncomps + c);
        const auto *fg_25 = buffer.data(fg + 25 * ncomps + c);
        const auto *fg_26 = buffer.data(fg + 26 * ncomps + c);
        const auto *fg_27 = buffer.data(fg + 27 * ncomps + c);
        const auto *fg_28 = buffer.data(fg + 28 * ncomps + c);
        const auto *fg_29 = buffer.data(fg + 29 * ncomps + c);
        const auto *fg_30 = buffer.data(fg + 30 * ncomps + c);
        const auto *fg_31 = buffer.data(fg + 31 * ncomps + c);
        const auto *fg_32 = buffer.data(fg + 32 * ncomps + c);
        const auto *fg_33 = buffer.data(fg + 33 * ncomps + c);
        const auto *fg_34 = buffer.data(fg + 34 * ncomps + c);
        const auto *fg_35 = buffer.data(fg + 35 * ncomps + c);
        const auto *fg_36 = buffer.data(fg + 36 * ncomps + c);
        const auto *fg_37 = buffer.data(fg + 37 * ncomps + c);
        const auto *fg_38 = buffer.data(fg + 38 * ncomps + c);
        const auto *fg_39 = buffer.data(fg + 39 * ncomps + c);
        const auto *fg_40 = buffer.data(fg + 40 * ncomps + c);
        const auto *fg_41 = buffer.data(fg + 41 * ncomps + c);
        const auto *fg_42 = buffer.data(fg + 42 * ncomps + c);
        const auto *fg_43 = buffer.data(fg + 43 * ncomps + c);
        const auto *fg_44 = buffer.data(fg + 44 * ncomps + c);
        const auto *fg_45 = buffer.data(fg + 45 * ncomps + c);
        const auto *fg_46 = buffer.data(fg + 46 * ncomps + c);
        const auto *fg_47 = buffer.data(fg + 47 * ncomps + c);
        const auto *fg_48 = buffer.data(fg + 48 * ncomps + c);
        const auto *fg_49 = buffer.data(fg + 49 * ncomps + c);
        const auto *fg_50 = buffer.data(fg + 50 * ncomps + c);
        const auto *fg_51 = buffer.data(fg + 51 * ncomps + c);
        const auto *fg_52 = buffer.data(fg + 52 * ncomps + c);
        const auto *fg_53 = buffer.data(fg + 53 * ncomps + c);
        const auto *fg_54 = buffer.data(fg + 54 * ncomps + c);
        const auto *fg_55 = buffer.data(fg + 55 * ncomps + c);
        const auto *fg_56 = buffer.data(fg + 56 * ncomps + c);
        const auto *fg_57 = buffer.data(fg + 57 * ncomps + c);
        const auto *fg_58 = buffer.data(fg + 58 * ncomps + c);
        const auto *fg_59 = buffer.data(fg + 59 * ncomps + c);
        const auto *fg_60 = buffer.data(fg + 60 * ncomps + c);
        const auto *fg_61 = buffer.data(fg + 61 * ncomps + c);
        const auto *fg_62 = buffer.data(fg + 62 * ncomps + c);
        const auto *fg_63 = buffer.data(fg + 63 * ncomps + c);
        const auto *fg_64 = buffer.data(fg + 64 * ncomps + c);
        const auto *fg_65 = buffer.data(fg + 65 * ncomps + c);
        const auto *fg_66 = buffer.data(fg + 66 * ncomps + c);
        const auto *fg_67 = buffer.data(fg + 67 * ncomps + c);
        const auto *fg_68 = buffer.data(fg + 68 * ncomps + c);
        const auto *fg_69 = buffer.data(fg + 69 * ncomps + c);
        const auto *fg_70 = buffer.data(fg + 70 * ncomps + c);
        const auto *fg_71 = buffer.data(fg + 71 * ncomps + c);
        const auto *fg_72 = buffer.data(fg + 72 * ncomps + c);
        const auto *fg_73 = buffer.data(fg + 73 * ncomps + c);
        const auto *fg_74 = buffer.data(fg + 74 * ncomps + c);
        const auto *fg_75 = buffer.data(fg + 75 * ncomps + c);
        const auto *fg_76 = buffer.data(fg + 76 * ncomps + c);
        const auto *fg_77 = buffer.data(fg + 77 * ncomps + c);
        const auto *fg_78 = buffer.data(fg + 78 * ncomps + c);
        const auto *fg_79 = buffer.data(fg + 79 * ncomps + c);
        const auto *fg_80 = buffer.data(fg + 80 * ncomps + c);
        const auto *fg_81 = buffer.data(fg + 81 * ncomps + c);
        const auto *fg_82 = buffer.data(fg + 82 * ncomps + c);
        const auto *fg_83 = buffer.data(fg + 83 * ncomps + c);
        const auto *fg_84 = buffer.data(fg + 84 * ncomps + c);
        const auto *fg_85 = buffer.data(fg + 85 * ncomps + c);
        const auto *fg_86 = buffer.data(fg + 86 * ncomps + c);
        const auto *fg_87 = buffer.data(fg + 87 * ncomps + c);
        const auto *fg_88 = buffer.data(fg + 88 * ncomps + c);
        const auto *fg_89 = buffer.data(fg + 89 * ncomps + c);
        const auto *fg_90 = buffer.data(fg + 90 * ncomps + c);
        const auto *fg_91 = buffer.data(fg + 91 * ncomps + c);
        const auto *fg_92 = buffer.data(fg + 92 * ncomps + c);
        const auto *fg_93 = buffer.data(fg + 93 * ncomps + c);
        const auto *fg_94 = buffer.data(fg + 94 * ncomps + c);
        const auto *fg_95 = buffer.data(fg + 95 * ncomps + c);
        const auto *fg_96 = buffer.data(fg + 96 * ncomps + c);
        const auto *fg_97 = buffer.data(fg + 97 * ncomps + c);
        const auto *fg_98 = buffer.data(fg + 98 * ncomps + c);
        const auto *fg_99 = buffer.data(fg + 99 * ncomps + c);
        const auto *fg_100 = buffer.data(fg + 100 * ncomps + c);
        const auto *fg_101 = buffer.data(fg + 101 * ncomps + c);
        const auto *fg_102 = buffer.data(fg + 102 * ncomps + c);
        const auto *fg_103 = buffer.data(fg + 103 * ncomps + c);
        const auto *fg_104 = buffer.data(fg + 104 * ncomps + c);
        const auto *fg_105 = buffer.data(fg + 105 * ncomps + c);
        const auto *fg_106 = buffer.data(fg + 106 * ncomps + c);
        const auto *fg_107 = buffer.data(fg + 107 * ncomps + c);
        const auto *fg_108 = buffer.data(fg + 108 * ncomps + c);
        const auto *fg_109 = buffer.data(fg + 109 * ncomps + c);
        const auto *fg_110 = buffer.data(fg + 110 * ncomps + c);
        const auto *fg_111 = buffer.data(fg + 111 * ncomps + c);
        const auto *fg_112 = buffer.data(fg + 112 * ncomps + c);
        const auto *fg_113 = buffer.data(fg + 113 * ncomps + c);
        const auto *fg_114 = buffer.data(fg + 114 * ncomps + c);
        const auto *fg_115 = buffer.data(fg + 115 * ncomps + c);
        const auto *fg_116 = buffer.data(fg + 116 * ncomps + c);
        const auto *fg_117 = buffer.data(fg + 117 * ncomps + c);
        const auto *fg_118 = buffer.data(fg + 118 * ncomps + c);
        const auto *fg_119 = buffer.data(fg + 119 * ncomps + c);
        const auto *fg_120 = buffer.data(fg + 120 * ncomps + c);
        const auto *fg_121 = buffer.data(fg + 121 * ncomps + c);
        const auto *fg_122 = buffer.data(fg + 122 * ncomps + c);
        const auto *fg_123 = buffer.data(fg + 123 * ncomps + c);
        const auto *fg_124 = buffer.data(fg + 124 * ncomps + c);
        const auto *fg_125 = buffer.data(fg + 125 * ncomps + c);
        const auto *fg_126 = buffer.data(fg + 126 * ncomps + c);
        const auto *fg_127 = buffer.data(fg + 127 * ncomps + c);
        const auto *fg_128 = buffer.data(fg + 128 * ncomps + c);
        const auto *fg_129 = buffer.data(fg + 129 * ncomps + c);
        const auto *fg_130 = buffer.data(fg + 130 * ncomps + c);
        const auto *fg_131 = buffer.data(fg + 131 * ncomps + c);
        const auto *fg_132 = buffer.data(fg + 132 * ncomps + c);
        const auto *fg_133 = buffer.data(fg + 133 * ncomps + c);
        const auto *fg_134 = buffer.data(fg + 134 * ncomps + c);
        const auto *fg_135 = buffer.data(fg + 135 * ncomps + c);
        const auto *fg_136 = buffer.data(fg + 136 * ncomps + c);
        const auto *fg_137 = buffer.data(fg + 137 * ncomps + c);
        const auto *fg_138 = buffer.data(fg + 138 * ncomps + c);
        const auto *fg_139 = buffer.data(fg + 139 * ncomps + c);
        const auto *fg_140 = buffer.data(fg + 140 * ncomps + c);
        const auto *fg_141 = buffer.data(fg + 141 * ncomps + c);
        const auto *fg_142 = buffer.data(fg + 142 * ncomps + c);
        const auto *fg_143 = buffer.data(fg + 143 * ncomps + c);
        const auto *fg_144 = buffer.data(fg + 144 * ncomps + c);

        const auto *fh_0 = buffer.data(fh + 0 * ncomps + c);
        const auto *fh_1 = buffer.data(fh + 1 * ncomps + c);
        const auto *fh_2 = buffer.data(fh + 2 * ncomps + c);
        const auto *fh_3 = buffer.data(fh + 3 * ncomps + c);
        const auto *fh_4 = buffer.data(fh + 4 * ncomps + c);
        const auto *fh_5 = buffer.data(fh + 5 * ncomps + c);
        const auto *fh_6 = buffer.data(fh + 6 * ncomps + c);
        const auto *fh_7 = buffer.data(fh + 7 * ncomps + c);
        const auto *fh_8 = buffer.data(fh + 8 * ncomps + c);
        const auto *fh_9 = buffer.data(fh + 9 * ncomps + c);
        const auto *fh_10 = buffer.data(fh + 10 * ncomps + c);
        const auto *fh_11 = buffer.data(fh + 11 * ncomps + c);
        const auto *fh_12 = buffer.data(fh + 12 * ncomps + c);
        const auto *fh_13 = buffer.data(fh + 13 * ncomps + c);
        const auto *fh_14 = buffer.data(fh + 14 * ncomps + c);
        const auto *fh_21 = buffer.data(fh + 21 * ncomps + c);
        const auto *fh_22 = buffer.data(fh + 22 * ncomps + c);
        const auto *fh_23 = buffer.data(fh + 23 * ncomps + c);
        const auto *fh_24 = buffer.data(fh + 24 * ncomps + c);
        const auto *fh_25 = buffer.data(fh + 25 * ncomps + c);
        const auto *fh_26 = buffer.data(fh + 26 * ncomps + c);
        const auto *fh_27 = buffer.data(fh + 27 * ncomps + c);
        const auto *fh_28 = buffer.data(fh + 28 * ncomps + c);
        const auto *fh_29 = buffer.data(fh + 29 * ncomps + c);
        const auto *fh_30 = buffer.data(fh + 30 * ncomps + c);
        const auto *fh_31 = buffer.data(fh + 31 * ncomps + c);
        const auto *fh_32 = buffer.data(fh + 32 * ncomps + c);
        const auto *fh_33 = buffer.data(fh + 33 * ncomps + c);
        const auto *fh_34 = buffer.data(fh + 34 * ncomps + c);
        const auto *fh_35 = buffer.data(fh + 35 * ncomps + c);
        const auto *fh_42 = buffer.data(fh + 42 * ncomps + c);
        const auto *fh_43 = buffer.data(fh + 43 * ncomps + c);
        const auto *fh_44 = buffer.data(fh + 44 * ncomps + c);
        const auto *fh_45 = buffer.data(fh + 45 * ncomps + c);
        const auto *fh_46 = buffer.data(fh + 46 * ncomps + c);
        const auto *fh_47 = buffer.data(fh + 47 * ncomps + c);
        const auto *fh_48 = buffer.data(fh + 48 * ncomps + c);
        const auto *fh_49 = buffer.data(fh + 49 * ncomps + c);
        const auto *fh_50 = buffer.data(fh + 50 * ncomps + c);
        const auto *fh_51 = buffer.data(fh + 51 * ncomps + c);
        const auto *fh_52 = buffer.data(fh + 52 * ncomps + c);
        const auto *fh_53 = buffer.data(fh + 53 * ncomps + c);
        const auto *fh_54 = buffer.data(fh + 54 * ncomps + c);
        const auto *fh_55 = buffer.data(fh + 55 * ncomps + c);
        const auto *fh_56 = buffer.data(fh + 56 * ncomps + c);
        const auto *fh_63 = buffer.data(fh + 63 * ncomps + c);
        const auto *fh_64 = buffer.data(fh + 64 * ncomps + c);
        const auto *fh_65 = buffer.data(fh + 65 * ncomps + c);
        const auto *fh_66 = buffer.data(fh + 66 * ncomps + c);
        const auto *fh_67 = buffer.data(fh + 67 * ncomps + c);
        const auto *fh_68 = buffer.data(fh + 68 * ncomps + c);
        const auto *fh_69 = buffer.data(fh + 69 * ncomps + c);
        const auto *fh_70 = buffer.data(fh + 70 * ncomps + c);
        const auto *fh_71 = buffer.data(fh + 71 * ncomps + c);
        const auto *fh_72 = buffer.data(fh + 72 * ncomps + c);
        const auto *fh_73 = buffer.data(fh + 73 * ncomps + c);
        const auto *fh_74 = buffer.data(fh + 74 * ncomps + c);
        const auto *fh_75 = buffer.data(fh + 75 * ncomps + c);
        const auto *fh_76 = buffer.data(fh + 76 * ncomps + c);
        const auto *fh_77 = buffer.data(fh + 77 * ncomps + c);
        const auto *fh_84 = buffer.data(fh + 84 * ncomps + c);
        const auto *fh_85 = buffer.data(fh + 85 * ncomps + c);
        const auto *fh_86 = buffer.data(fh + 86 * ncomps + c);
        const auto *fh_87 = buffer.data(fh + 87 * ncomps + c);
        const auto *fh_88 = buffer.data(fh + 88 * ncomps + c);
        const auto *fh_89 = buffer.data(fh + 89 * ncomps + c);
        const auto *fh_90 = buffer.data(fh + 90 * ncomps + c);
        const auto *fh_91 = buffer.data(fh + 91 * ncomps + c);
        const auto *fh_92 = buffer.data(fh + 92 * ncomps + c);
        const auto *fh_93 = buffer.data(fh + 93 * ncomps + c);
        const auto *fh_94 = buffer.data(fh + 94 * ncomps + c);
        const auto *fh_95 = buffer.data(fh + 95 * ncomps + c);
        const auto *fh_96 = buffer.data(fh + 96 * ncomps + c);
        const auto *fh_97 = buffer.data(fh + 97 * ncomps + c);
        const auto *fh_98 = buffer.data(fh + 98 * ncomps + c);
        const auto *fh_105 = buffer.data(fh + 105 * ncomps + c);
        const auto *fh_106 = buffer.data(fh + 106 * ncomps + c);
        const auto *fh_107 = buffer.data(fh + 107 * ncomps + c);
        const auto *fh_108 = buffer.data(fh + 108 * ncomps + c);
        const auto *fh_109 = buffer.data(fh + 109 * ncomps + c);
        const auto *fh_110 = buffer.data(fh + 110 * ncomps + c);
        const auto *fh_111 = buffer.data(fh + 111 * ncomps + c);
        const auto *fh_112 = buffer.data(fh + 112 * ncomps + c);
        const auto *fh_113 = buffer.data(fh + 113 * ncomps + c);
        const auto *fh_114 = buffer.data(fh + 114 * ncomps + c);
        const auto *fh_115 = buffer.data(fh + 115 * ncomps + c);
        const auto *fh_116 = buffer.data(fh + 116 * ncomps + c);
        const auto *fh_117 = buffer.data(fh + 117 * ncomps + c);
        const auto *fh_118 = buffer.data(fh + 118 * ncomps + c);
        const auto *fh_119 = buffer.data(fh + 119 * ncomps + c);
        const auto *fh_126 = buffer.data(fh + 126 * ncomps + c);
        const auto *fh_127 = buffer.data(fh + 127 * ncomps + c);
        const auto *fh_128 = buffer.data(fh + 128 * ncomps + c);
        const auto *fh_129 = buffer.data(fh + 129 * ncomps + c);
        const auto *fh_130 = buffer.data(fh + 130 * ncomps + c);
        const auto *fh_131 = buffer.data(fh + 131 * ncomps + c);
        const auto *fh_132 = buffer.data(fh + 132 * ncomps + c);
        const auto *fh_133 = buffer.data(fh + 133 * ncomps + c);
        const auto *fh_134 = buffer.data(fh + 134 * ncomps + c);
        const auto *fh_135 = buffer.data(fh + 135 * ncomps + c);
        const auto *fh_136 = buffer.data(fh + 136 * ncomps + c);
        const auto *fh_137 = buffer.data(fh + 137 * ncomps + c);
        const auto *fh_138 = buffer.data(fh + 138 * ncomps + c);
        const auto *fh_139 = buffer.data(fh + 139 * ncomps + c);
        const auto *fh_140 = buffer.data(fh + 140 * ncomps + c);
        const auto *fh_147 = buffer.data(fh + 147 * ncomps + c);
        const auto *fh_148 = buffer.data(fh + 148 * ncomps + c);
        const auto *fh_149 = buffer.data(fh + 149 * ncomps + c);
        const auto *fh_150 = buffer.data(fh + 150 * ncomps + c);
        const auto *fh_151 = buffer.data(fh + 151 * ncomps + c);
        const auto *fh_152 = buffer.data(fh + 152 * ncomps + c);
        const auto *fh_153 = buffer.data(fh + 153 * ncomps + c);
        const auto *fh_154 = buffer.data(fh + 154 * ncomps + c);
        const auto *fh_155 = buffer.data(fh + 155 * ncomps + c);
        const auto *fh_156 = buffer.data(fh + 156 * ncomps + c);
        const auto *fh_157 = buffer.data(fh + 157 * ncomps + c);
        const auto *fh_158 = buffer.data(fh + 158 * ncomps + c);
        const auto *fh_159 = buffer.data(fh + 159 * ncomps + c);
        const auto *fh_160 = buffer.data(fh + 160 * ncomps + c);
        const auto *fh_161 = buffer.data(fh + 161 * ncomps + c);
        const auto *fh_168 = buffer.data(fh + 168 * ncomps + c);
        const auto *fh_169 = buffer.data(fh + 169 * ncomps + c);
        const auto *fh_170 = buffer.data(fh + 170 * ncomps + c);
        const auto *fh_171 = buffer.data(fh + 171 * ncomps + c);
        const auto *fh_172 = buffer.data(fh + 172 * ncomps + c);
        const auto *fh_173 = buffer.data(fh + 173 * ncomps + c);
        const auto *fh_174 = buffer.data(fh + 174 * ncomps + c);
        const auto *fh_175 = buffer.data(fh + 175 * ncomps + c);
        const auto *fh_176 = buffer.data(fh + 176 * ncomps + c);
        const auto *fh_177 = buffer.data(fh + 177 * ncomps + c);
        const auto *fh_178 = buffer.data(fh + 178 * ncomps + c);
        const auto *fh_179 = buffer.data(fh + 179 * ncomps + c);
        const auto *fh_180 = buffer.data(fh + 180 * ncomps + c);
        const auto *fh_181 = buffer.data(fh + 181 * ncomps + c);
        const auto *fh_182 = buffer.data(fh + 182 * ncomps + c);
        const auto *fh_189 = buffer.data(fh + 189 * ncomps + c);
        const auto *fh_190 = buffer.data(fh + 190 * ncomps + c);
        const auto *fh_191 = buffer.data(fh + 191 * ncomps + c);
        const auto *fh_192 = buffer.data(fh + 192 * ncomps + c);
        const auto *fh_193 = buffer.data(fh + 193 * ncomps + c);
        const auto *fh_194 = buffer.data(fh + 194 * ncomps + c);
        const auto *fh_195 = buffer.data(fh + 195 * ncomps + c);
        const auto *fh_196 = buffer.data(fh + 196 * ncomps + c);
        const auto *fh_197 = buffer.data(fh + 197 * ncomps + c);
        const auto *fh_198 = buffer.data(fh + 198 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fg_0, fg_1, fg_2, fg_3, fg_4, fh_0, \
                         fh_1, fh_2, fh_3, fh_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * fg_0[k]
                     + fh_0[k];

            t_1[k] = -ab_x[k] * fg_1[k]
                     + fh_1[k];

            t_2[k] = -ab_x[k] * fg_2[k]
                     + fh_2[k];

            t_3[k] = -ab_x[k] * fg_3[k]
                     + fh_3[k];

            t_4[k] = -ab_x[k] * fg_4[k]
                     + fh_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fg_5, fg_6, fg_7, fg_8, fg_9, fh_5, \
                         fh_6, fh_7, fh_8, fh_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * fg_5[k]
                     + fh_5[k];

            t_6[k] = -ab_x[k] * fg_6[k]
                     + fh_6[k];

            t_7[k] = -ab_x[k] * fg_7[k]
                     + fh_7[k];

            t_8[k] = -ab_x[k] * fg_8[k]
                     + fh_8[k];

            t_9[k] = -ab_x[k] * fg_9[k]
                     + fh_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fg_10, fg_11, fg_12, fg_13, \
                         fg_14, fh_10, fh_11, fh_12, fh_13, fh_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * fg_10[k]
                      + fh_10[k];

            t_11[k] = -ab_x[k] * fg_11[k]
                      + fh_11[k];

            t_12[k] = -ab_x[k] * fg_12[k]
                      + fh_12[k];

            t_13[k] = -ab_x[k] * fg_13[k]
                      + fh_13[k];

            t_14[k] = -ab_x[k] * fg_14[k]
                      + fh_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fg_15, fg_16, fg_17, fg_18, \
                         fg_19, fh_21, fh_22, fh_23, fh_24, fh_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * fg_15[k]
                      + fh_21[k];

            t_16[k] = -ab_x[k] * fg_16[k]
                      + fh_22[k];

            t_17[k] = -ab_x[k] * fg_17[k]
                      + fh_23[k];

            t_18[k] = -ab_x[k] * fg_18[k]
                      + fh_24[k];

            t_19[k] = -ab_x[k] * fg_19[k]
                      + fh_25[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fg_20, fg_21, fg_22, fg_23, \
                         fg_24, fh_26, fh_27, fh_28, fh_29, fh_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * fg_20[k]
                      + fh_26[k];

            t_21[k] = -ab_x[k] * fg_21[k]
                      + fh_27[k];

            t_22[k] = -ab_x[k] * fg_22[k]
                      + fh_28[k];

            t_23[k] = -ab_x[k] * fg_23[k]
                      + fh_29[k];

            t_24[k] = -ab_x[k] * fg_24[k]
                      + fh_30[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fg_25, fg_26, fg_27, fg_28, \
                         fg_29, fh_31, fh_32, fh_33, fh_34, fh_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * fg_25[k]
                      + fh_31[k];

            t_26[k] = -ab_x[k] * fg_26[k]
                      + fh_32[k];

            t_27[k] = -ab_x[k] * fg_27[k]
                      + fh_33[k];

            t_28[k] = -ab_x[k] * fg_28[k]
                      + fh_34[k];

            t_29[k] = -ab_x[k] * fg_29[k]
                      + fh_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fg_30, fg_31, fg_32, fg_33, \
                         fg_34, fh_42, fh_43, fh_44, fh_45, fh_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * fg_30[k]
                      + fh_42[k];

            t_31[k] = -ab_x[k] * fg_31[k]
                      + fh_43[k];

            t_32[k] = -ab_x[k] * fg_32[k]
                      + fh_44[k];

            t_33[k] = -ab_x[k] * fg_33[k]
                      + fh_45[k];

            t_34[k] = -ab_x[k] * fg_34[k]
                      + fh_46[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, fg_35, fg_36, fg_37, fg_38, \
                         fg_39, fh_47, fh_48, fh_49, fh_50, fh_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * fg_35[k]
                      + fh_47[k];

            t_36[k] = -ab_x[k] * fg_36[k]
                      + fh_48[k];

            t_37[k] = -ab_x[k] * fg_37[k]
                      + fh_49[k];

            t_38[k] = -ab_x[k] * fg_38[k]
                      + fh_50[k];

            t_39[k] = -ab_x[k] * fg_39[k]
                      + fh_51[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fg_40, fg_41, fg_42, fg_43, \
                         fg_44, fh_52, fh_53, fh_54, fh_55, fh_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * fg_40[k]
                      + fh_52[k];

            t_41[k] = -ab_x[k] * fg_41[k]
                      + fh_53[k];

            t_42[k] = -ab_x[k] * fg_42[k]
                      + fh_54[k];

            t_43[k] = -ab_x[k] * fg_43[k]
                      + fh_55[k];

            t_44[k] = -ab_x[k] * fg_44[k]
                      + fh_56[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, fg_45, fg_46, fg_47, fg_48, \
                         fg_49, fh_63, fh_64, fh_65, fh_66, fh_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * fg_45[k]
                      + fh_63[k];

            t_46[k] = -ab_x[k] * fg_46[k]
                      + fh_64[k];

            t_47[k] = -ab_x[k] * fg_47[k]
                      + fh_65[k];

            t_48[k] = -ab_x[k] * fg_48[k]
                      + fh_66[k];

            t_49[k] = -ab_x[k] * fg_49[k]
                      + fh_67[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fg_50, fg_51, fg_52, fg_53, \
                         fg_54, fh_68, fh_69, fh_70, fh_71, fh_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * fg_50[k]
                      + fh_68[k];

            t_51[k] = -ab_x[k] * fg_51[k]
                      + fh_69[k];

            t_52[k] = -ab_x[k] * fg_52[k]
                      + fh_70[k];

            t_53[k] = -ab_x[k] * fg_53[k]
                      + fh_71[k];

            t_54[k] = -ab_x[k] * fg_54[k]
                      + fh_72[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, fg_55, fg_56, fg_57, fg_58, \
                         fg_59, fh_73, fh_74, fh_75, fh_76, fh_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * fg_55[k]
                      + fh_73[k];

            t_56[k] = -ab_x[k] * fg_56[k]
                      + fh_74[k];

            t_57[k] = -ab_x[k] * fg_57[k]
                      + fh_75[k];

            t_58[k] = -ab_x[k] * fg_58[k]
                      + fh_76[k];

            t_59[k] = -ab_x[k] * fg_59[k]
                      + fh_77[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fg_60, fg_61, fg_62, fg_63, \
                         fg_64, fh_84, fh_85, fh_86, fh_87, fh_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * fg_60[k]
                      + fh_84[k];

            t_61[k] = -ab_x[k] * fg_61[k]
                      + fh_85[k];

            t_62[k] = -ab_x[k] * fg_62[k]
                      + fh_86[k];

            t_63[k] = -ab_x[k] * fg_63[k]
                      + fh_87[k];

            t_64[k] = -ab_x[k] * fg_64[k]
                      + fh_88[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, fg_65, fg_66, fg_67, fg_68, \
                         fg_69, fh_89, fh_90, fh_91, fh_92, fh_93 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * fg_65[k]
                      + fh_89[k];

            t_66[k] = -ab_x[k] * fg_66[k]
                      + fh_90[k];

            t_67[k] = -ab_x[k] * fg_67[k]
                      + fh_91[k];

            t_68[k] = -ab_x[k] * fg_68[k]
                      + fh_92[k];

            t_69[k] = -ab_x[k] * fg_69[k]
                      + fh_93[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fg_70, fg_71, fg_72, fg_73, \
                         fg_74, fh_94, fh_95, fh_96, fh_97, fh_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * fg_70[k]
                      + fh_94[k];

            t_71[k] = -ab_x[k] * fg_71[k]
                      + fh_95[k];

            t_72[k] = -ab_x[k] * fg_72[k]
                      + fh_96[k];

            t_73[k] = -ab_x[k] * fg_73[k]
                      + fh_97[k];

            t_74[k] = -ab_x[k] * fg_74[k]
                      + fh_98[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, fg_75, fg_76, fg_77, fg_78, \
                         fg_79, fh_105, fh_106, fh_107, fh_108, \
                         fh_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * fg_75[k]
                      + fh_105[k];

            t_76[k] = -ab_x[k] * fg_76[k]
                      + fh_106[k];

            t_77[k] = -ab_x[k] * fg_77[k]
                      + fh_107[k];

            t_78[k] = -ab_x[k] * fg_78[k]
                      + fh_108[k];

            t_79[k] = -ab_x[k] * fg_79[k]
                      + fh_109[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fg_80, fg_81, fg_82, fg_83, \
                         fg_84, fh_110, fh_111, fh_112, fh_113, \
                         fh_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * fg_80[k]
                      + fh_110[k];

            t_81[k] = -ab_x[k] * fg_81[k]
                      + fh_111[k];

            t_82[k] = -ab_x[k] * fg_82[k]
                      + fh_112[k];

            t_83[k] = -ab_x[k] * fg_83[k]
                      + fh_113[k];

            t_84[k] = -ab_x[k] * fg_84[k]
                      + fh_114[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, fg_85, fg_86, fg_87, fg_88, \
                         fg_89, fh_115, fh_116, fh_117, fh_118, \
                         fh_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * fg_85[k]
                      + fh_115[k];

            t_86[k] = -ab_x[k] * fg_86[k]
                      + fh_116[k];

            t_87[k] = -ab_x[k] * fg_87[k]
                      + fh_117[k];

            t_88[k] = -ab_x[k] * fg_88[k]
                      + fh_118[k];

            t_89[k] = -ab_x[k] * fg_89[k]
                      + fh_119[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fg_90, fg_91, fg_92, fg_93, \
                         fg_94, fh_126, fh_127, fh_128, fh_129, \
                         fh_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * fg_90[k]
                      + fh_126[k];

            t_91[k] = -ab_x[k] * fg_91[k]
                      + fh_127[k];

            t_92[k] = -ab_x[k] * fg_92[k]
                      + fh_128[k];

            t_93[k] = -ab_x[k] * fg_93[k]
                      + fh_129[k];

            t_94[k] = -ab_x[k] * fg_94[k]
                      + fh_130[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, fg_95, fg_96, fg_97, fg_98, \
                         fg_99, fh_131, fh_132, fh_133, fh_134, \
                         fh_135 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * fg_95[k]
                      + fh_131[k];

            t_96[k] = -ab_x[k] * fg_96[k]
                      + fh_132[k];

            t_97[k] = -ab_x[k] * fg_97[k]
                      + fh_133[k];

            t_98[k] = -ab_x[k] * fg_98[k]
                      + fh_134[k];

            t_99[k] = -ab_x[k] * fg_99[k]
                      + fh_135[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, fg_100, fg_101, fg_102, \
                         fg_103, fg_104, fh_136, fh_137, fh_138, fh_139, \
                         fh_140 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * fg_100[k]
                       + fh_136[k];

            t_101[k] = -ab_x[k] * fg_101[k]
                       + fh_137[k];

            t_102[k] = -ab_x[k] * fg_102[k]
                       + fh_138[k];

            t_103[k] = -ab_x[k] * fg_103[k]
                       + fh_139[k];

            t_104[k] = -ab_x[k] * fg_104[k]
                       + fh_140[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, fg_105, fg_106, fg_107, \
                         fg_108, fg_109, fh_147, fh_148, fh_149, fh_150, \
                         fh_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * fg_105[k]
                       + fh_147[k];

            t_106[k] = -ab_x[k] * fg_106[k]
                       + fh_148[k];

            t_107[k] = -ab_x[k] * fg_107[k]
                       + fh_149[k];

            t_108[k] = -ab_x[k] * fg_108[k]
                       + fh_150[k];

            t_109[k] = -ab_x[k] * fg_109[k]
                       + fh_151[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, fg_110, fg_111, fg_112, \
                         fg_113, fg_114, fh_152, fh_153, fh_154, fh_155, \
                         fh_156 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * fg_110[k]
                       + fh_152[k];

            t_111[k] = -ab_x[k] * fg_111[k]
                       + fh_153[k];

            t_112[k] = -ab_x[k] * fg_112[k]
                       + fh_154[k];

            t_113[k] = -ab_x[k] * fg_113[k]
                       + fh_155[k];

            t_114[k] = -ab_x[k] * fg_114[k]
                       + fh_156[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, fg_115, fg_116, fg_117, \
                         fg_118, fg_119, fh_157, fh_158, fh_159, fh_160, \
                         fh_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * fg_115[k]
                       + fh_157[k];

            t_116[k] = -ab_x[k] * fg_116[k]
                       + fh_158[k];

            t_117[k] = -ab_x[k] * fg_117[k]
                       + fh_159[k];

            t_118[k] = -ab_x[k] * fg_118[k]
                       + fh_160[k];

            t_119[k] = -ab_x[k] * fg_119[k]
                       + fh_161[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, fg_120, fg_121, fg_122, \
                         fg_123, fg_124, fh_168, fh_169, fh_170, fh_171, \
                         fh_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * fg_120[k]
                       + fh_168[k];

            t_121[k] = -ab_x[k] * fg_121[k]
                       + fh_169[k];

            t_122[k] = -ab_x[k] * fg_122[k]
                       + fh_170[k];

            t_123[k] = -ab_x[k] * fg_123[k]
                       + fh_171[k];

            t_124[k] = -ab_x[k] * fg_124[k]
                       + fh_172[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, fg_125, fg_126, fg_127, \
                         fg_128, fg_129, fh_173, fh_174, fh_175, fh_176, \
                         fh_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * fg_125[k]
                       + fh_173[k];

            t_126[k] = -ab_x[k] * fg_126[k]
                       + fh_174[k];

            t_127[k] = -ab_x[k] * fg_127[k]
                       + fh_175[k];

            t_128[k] = -ab_x[k] * fg_128[k]
                       + fh_176[k];

            t_129[k] = -ab_x[k] * fg_129[k]
                       + fh_177[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, fg_130, fg_131, fg_132, \
                         fg_133, fg_134, fh_178, fh_179, fh_180, fh_181, \
                         fh_182 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * fg_130[k]
                       + fh_178[k];

            t_131[k] = -ab_x[k] * fg_131[k]
                       + fh_179[k];

            t_132[k] = -ab_x[k] * fg_132[k]
                       + fh_180[k];

            t_133[k] = -ab_x[k] * fg_133[k]
                       + fh_181[k];

            t_134[k] = -ab_x[k] * fg_134[k]
                       + fh_182[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, fg_135, fg_136, fg_137, \
                         fg_138, fg_139, fh_189, fh_190, fh_191, fh_192, \
                         fh_193 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * fg_135[k]
                       + fh_189[k];

            t_136[k] = -ab_x[k] * fg_136[k]
                       + fh_190[k];

            t_137[k] = -ab_x[k] * fg_137[k]
                       + fh_191[k];

            t_138[k] = -ab_x[k] * fg_138[k]
                       + fh_192[k];

            t_139[k] = -ab_x[k] * fg_139[k]
                       + fh_193[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, fg_140, fg_141, fg_142, \
                         fg_143, fg_144, fh_194, fh_195, fh_196, fh_197, \
                         fh_198 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * fg_140[k]
                       + fh_194[k];

            t_141[k] = -ab_x[k] * fg_141[k]
                       + fh_195[k];

            t_142[k] = -ab_x[k] * fg_142[k]
                       + fh_196[k];

            t_143[k] = -ab_x[k] * fg_143[k]
                       + fh_197[k];

            t_144[k] = -ab_x[k] * fg_144[k]
                       + fh_198[k];
        }
    }
}

static auto
compute_hrr_gg_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t fg, const size_t fh, const size_t ncomps,
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

        const auto *fg_90 = buffer.data(fg + 90 * ncomps + c);
        const auto *fg_91 = buffer.data(fg + 91 * ncomps + c);
        const auto *fg_92 = buffer.data(fg + 92 * ncomps + c);
        const auto *fg_93 = buffer.data(fg + 93 * ncomps + c);
        const auto *fg_94 = buffer.data(fg + 94 * ncomps + c);
        const auto *fg_95 = buffer.data(fg + 95 * ncomps + c);
        const auto *fg_96 = buffer.data(fg + 96 * ncomps + c);
        const auto *fg_97 = buffer.data(fg + 97 * ncomps + c);
        const auto *fg_98 = buffer.data(fg + 98 * ncomps + c);
        const auto *fg_99 = buffer.data(fg + 99 * ncomps + c);
        const auto *fg_100 = buffer.data(fg + 100 * ncomps + c);
        const auto *fg_101 = buffer.data(fg + 101 * ncomps + c);
        const auto *fg_102 = buffer.data(fg + 102 * ncomps + c);
        const auto *fg_103 = buffer.data(fg + 103 * ncomps + c);
        const auto *fg_104 = buffer.data(fg + 104 * ncomps + c);
        const auto *fg_105 = buffer.data(fg + 105 * ncomps + c);
        const auto *fg_106 = buffer.data(fg + 106 * ncomps + c);
        const auto *fg_107 = buffer.data(fg + 107 * ncomps + c);
        const auto *fg_108 = buffer.data(fg + 108 * ncomps + c);
        const auto *fg_109 = buffer.data(fg + 109 * ncomps + c);
        const auto *fg_110 = buffer.data(fg + 110 * ncomps + c);
        const auto *fg_111 = buffer.data(fg + 111 * ncomps + c);
        const auto *fg_112 = buffer.data(fg + 112 * ncomps + c);
        const auto *fg_113 = buffer.data(fg + 113 * ncomps + c);
        const auto *fg_114 = buffer.data(fg + 114 * ncomps + c);
        const auto *fg_115 = buffer.data(fg + 115 * ncomps + c);
        const auto *fg_116 = buffer.data(fg + 116 * ncomps + c);
        const auto *fg_117 = buffer.data(fg + 117 * ncomps + c);
        const auto *fg_118 = buffer.data(fg + 118 * ncomps + c);
        const auto *fg_119 = buffer.data(fg + 119 * ncomps + c);
        const auto *fg_120 = buffer.data(fg + 120 * ncomps + c);
        const auto *fg_121 = buffer.data(fg + 121 * ncomps + c);
        const auto *fg_122 = buffer.data(fg + 122 * ncomps + c);
        const auto *fg_123 = buffer.data(fg + 123 * ncomps + c);
        const auto *fg_124 = buffer.data(fg + 124 * ncomps + c);
        const auto *fg_125 = buffer.data(fg + 125 * ncomps + c);
        const auto *fg_126 = buffer.data(fg + 126 * ncomps + c);
        const auto *fg_127 = buffer.data(fg + 127 * ncomps + c);
        const auto *fg_128 = buffer.data(fg + 128 * ncomps + c);
        const auto *fg_129 = buffer.data(fg + 129 * ncomps + c);
        const auto *fg_130 = buffer.data(fg + 130 * ncomps + c);
        const auto *fg_131 = buffer.data(fg + 131 * ncomps + c);
        const auto *fg_132 = buffer.data(fg + 132 * ncomps + c);
        const auto *fg_133 = buffer.data(fg + 133 * ncomps + c);
        const auto *fg_134 = buffer.data(fg + 134 * ncomps + c);
        const auto *fg_135 = buffer.data(fg + 135 * ncomps + c);
        const auto *fg_136 = buffer.data(fg + 136 * ncomps + c);
        const auto *fg_137 = buffer.data(fg + 137 * ncomps + c);
        const auto *fg_138 = buffer.data(fg + 138 * ncomps + c);
        const auto *fg_139 = buffer.data(fg + 139 * ncomps + c);
        const auto *fg_140 = buffer.data(fg + 140 * ncomps + c);
        const auto *fg_141 = buffer.data(fg + 141 * ncomps + c);
        const auto *fg_142 = buffer.data(fg + 142 * ncomps + c);
        const auto *fg_143 = buffer.data(fg + 143 * ncomps + c);
        const auto *fg_144 = buffer.data(fg + 144 * ncomps + c);
        const auto *fg_145 = buffer.data(fg + 145 * ncomps + c);
        const auto *fg_146 = buffer.data(fg + 146 * ncomps + c);
        const auto *fg_147 = buffer.data(fg + 147 * ncomps + c);
        const auto *fg_148 = buffer.data(fg + 148 * ncomps + c);
        const auto *fg_149 = buffer.data(fg + 149 * ncomps + c);

        const auto *fh_127 = buffer.data(fh + 127 * ncomps + c);
        const auto *fh_129 = buffer.data(fh + 129 * ncomps + c);
        const auto *fh_130 = buffer.data(fh + 130 * ncomps + c);
        const auto *fh_132 = buffer.data(fh + 132 * ncomps + c);
        const auto *fh_133 = buffer.data(fh + 133 * ncomps + c);
        const auto *fh_134 = buffer.data(fh + 134 * ncomps + c);
        const auto *fh_136 = buffer.data(fh + 136 * ncomps + c);
        const auto *fh_137 = buffer.data(fh + 137 * ncomps + c);
        const auto *fh_138 = buffer.data(fh + 138 * ncomps + c);
        const auto *fh_139 = buffer.data(fh + 139 * ncomps + c);
        const auto *fh_141 = buffer.data(fh + 141 * ncomps + c);
        const auto *fh_142 = buffer.data(fh + 142 * ncomps + c);
        const auto *fh_143 = buffer.data(fh + 143 * ncomps + c);
        const auto *fh_144 = buffer.data(fh + 144 * ncomps + c);
        const auto *fh_145 = buffer.data(fh + 145 * ncomps + c);
        const auto *fh_148 = buffer.data(fh + 148 * ncomps + c);
        const auto *fh_150 = buffer.data(fh + 150 * ncomps + c);
        const auto *fh_151 = buffer.data(fh + 151 * ncomps + c);
        const auto *fh_153 = buffer.data(fh + 153 * ncomps + c);
        const auto *fh_154 = buffer.data(fh + 154 * ncomps + c);
        const auto *fh_155 = buffer.data(fh + 155 * ncomps + c);
        const auto *fh_157 = buffer.data(fh + 157 * ncomps + c);
        const auto *fh_158 = buffer.data(fh + 158 * ncomps + c);
        const auto *fh_159 = buffer.data(fh + 159 * ncomps + c);
        const auto *fh_160 = buffer.data(fh + 160 * ncomps + c);
        const auto *fh_162 = buffer.data(fh + 162 * ncomps + c);
        const auto *fh_163 = buffer.data(fh + 163 * ncomps + c);
        const auto *fh_164 = buffer.data(fh + 164 * ncomps + c);
        const auto *fh_165 = buffer.data(fh + 165 * ncomps + c);
        const auto *fh_166 = buffer.data(fh + 166 * ncomps + c);
        const auto *fh_169 = buffer.data(fh + 169 * ncomps + c);
        const auto *fh_171 = buffer.data(fh + 171 * ncomps + c);
        const auto *fh_172 = buffer.data(fh + 172 * ncomps + c);
        const auto *fh_174 = buffer.data(fh + 174 * ncomps + c);
        const auto *fh_175 = buffer.data(fh + 175 * ncomps + c);
        const auto *fh_176 = buffer.data(fh + 176 * ncomps + c);
        const auto *fh_178 = buffer.data(fh + 178 * ncomps + c);
        const auto *fh_179 = buffer.data(fh + 179 * ncomps + c);
        const auto *fh_180 = buffer.data(fh + 180 * ncomps + c);
        const auto *fh_181 = buffer.data(fh + 181 * ncomps + c);
        const auto *fh_183 = buffer.data(fh + 183 * ncomps + c);
        const auto *fh_184 = buffer.data(fh + 184 * ncomps + c);
        const auto *fh_185 = buffer.data(fh + 185 * ncomps + c);
        const auto *fh_186 = buffer.data(fh + 186 * ncomps + c);
        const auto *fh_187 = buffer.data(fh + 187 * ncomps + c);
        const auto *fh_190 = buffer.data(fh + 190 * ncomps + c);
        const auto *fh_191 = buffer.data(fh + 191 * ncomps + c);
        const auto *fh_192 = buffer.data(fh + 192 * ncomps + c);
        const auto *fh_193 = buffer.data(fh + 193 * ncomps + c);
        const auto *fh_194 = buffer.data(fh + 194 * ncomps + c);
        const auto *fh_195 = buffer.data(fh + 195 * ncomps + c);
        const auto *fh_196 = buffer.data(fh + 196 * ncomps + c);
        const auto *fh_197 = buffer.data(fh + 197 * ncomps + c);
        const auto *fh_198 = buffer.data(fh + 198 * ncomps + c);
        const auto *fh_199 = buffer.data(fh + 199 * ncomps + c);
        const auto *fh_200 = buffer.data(fh + 200 * ncomps + c);
        const auto *fh_201 = buffer.data(fh + 201 * ncomps + c);
        const auto *fh_202 = buffer.data(fh + 202 * ncomps + c);
        const auto *fh_203 = buffer.data(fh + 203 * ncomps + c);
        const auto *fh_204 = buffer.data(fh + 204 * ncomps + c);
        const auto *fh_205 = buffer.data(fh + 205 * ncomps + c);
        const auto *fh_206 = buffer.data(fh + 206 * ncomps + c);
        const auto *fh_207 = buffer.data(fh + 207 * ncomps + c);
        const auto *fh_208 = buffer.data(fh + 208 * ncomps + c);
        const auto *fh_209 = buffer.data(fh + 209 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, fg_145, fg_146, fg_147, \
                         fg_148, fg_149, fh_199, fh_200, fh_201, fh_202, \
                         fh_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * fg_145[k]
                       + fh_199[k];

            t_146[k] = -ab_x[k] * fg_146[k]
                       + fh_200[k];

            t_147[k] = -ab_x[k] * fg_147[k]
                       + fh_201[k];

            t_148[k] = -ab_x[k] * fg_148[k]
                       + fh_202[k];

            t_149[k] = -ab_x[k] * fg_149[k]
                       + fh_203[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_y, fg_90, fg_91, fg_92, fg_93, \
                         fg_94, fh_127, fh_129, fh_130, fh_132, \
                         fh_133 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_y[k] * fg_90[k]
                       + fh_127[k];

            t_151[k] = -ab_y[k] * fg_91[k]
                       + fh_129[k];

            t_152[k] = -ab_y[k] * fg_92[k]
                       + fh_130[k];

            t_153[k] = -ab_y[k] * fg_93[k]
                       + fh_132[k];

            t_154[k] = -ab_y[k] * fg_94[k]
                       + fh_133[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_y, fg_95, fg_96, fg_97, fg_98, \
                         fg_99, fh_134, fh_136, fh_137, fh_138, \
                         fh_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_y[k] * fg_95[k]
                       + fh_134[k];

            t_156[k] = -ab_y[k] * fg_96[k]
                       + fh_136[k];

            t_157[k] = -ab_y[k] * fg_97[k]
                       + fh_137[k];

            t_158[k] = -ab_y[k] * fg_98[k]
                       + fh_138[k];

            t_159[k] = -ab_y[k] * fg_99[k]
                       + fh_139[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, fg_100, fg_101, fg_102, \
                         fg_103, fg_104, fh_141, fh_142, fh_143, fh_144, \
                         fh_145 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_y[k] * fg_100[k]
                       + fh_141[k];

            t_161[k] = -ab_y[k] * fg_101[k]
                       + fh_142[k];

            t_162[k] = -ab_y[k] * fg_102[k]
                       + fh_143[k];

            t_163[k] = -ab_y[k] * fg_103[k]
                       + fh_144[k];

            t_164[k] = -ab_y[k] * fg_104[k]
                       + fh_145[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_y, fg_105, fg_106, fg_107, \
                         fg_108, fg_109, fh_148, fh_150, fh_151, fh_153, \
                         fh_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_y[k] * fg_105[k]
                       + fh_148[k];

            t_166[k] = -ab_y[k] * fg_106[k]
                       + fh_150[k];

            t_167[k] = -ab_y[k] * fg_107[k]
                       + fh_151[k];

            t_168[k] = -ab_y[k] * fg_108[k]
                       + fh_153[k];

            t_169[k] = -ab_y[k] * fg_109[k]
                       + fh_154[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_y, fg_110, fg_111, fg_112, \
                         fg_113, fg_114, fh_155, fh_157, fh_158, fh_159, \
                         fh_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_y[k] * fg_110[k]
                       + fh_155[k];

            t_171[k] = -ab_y[k] * fg_111[k]
                       + fh_157[k];

            t_172[k] = -ab_y[k] * fg_112[k]
                       + fh_158[k];

            t_173[k] = -ab_y[k] * fg_113[k]
                       + fh_159[k];

            t_174[k] = -ab_y[k] * fg_114[k]
                       + fh_160[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, fg_115, fg_116, fg_117, \
                         fg_118, fg_119, fh_162, fh_163, fh_164, fh_165, \
                         fh_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_y[k] * fg_115[k]
                       + fh_162[k];

            t_176[k] = -ab_y[k] * fg_116[k]
                       + fh_163[k];

            t_177[k] = -ab_y[k] * fg_117[k]
                       + fh_164[k];

            t_178[k] = -ab_y[k] * fg_118[k]
                       + fh_165[k];

            t_179[k] = -ab_y[k] * fg_119[k]
                       + fh_166[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_y, fg_120, fg_121, fg_122, \
                         fg_123, fg_124, fh_169, fh_171, fh_172, fh_174, \
                         fh_175 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_y[k] * fg_120[k]
                       + fh_169[k];

            t_181[k] = -ab_y[k] * fg_121[k]
                       + fh_171[k];

            t_182[k] = -ab_y[k] * fg_122[k]
                       + fh_172[k];

            t_183[k] = -ab_y[k] * fg_123[k]
                       + fh_174[k];

            t_184[k] = -ab_y[k] * fg_124[k]
                       + fh_175[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_y, fg_125, fg_126, fg_127, \
                         fg_128, fg_129, fh_176, fh_178, fh_179, fh_180, \
                         fh_181 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_y[k] * fg_125[k]
                       + fh_176[k];

            t_186[k] = -ab_y[k] * fg_126[k]
                       + fh_178[k];

            t_187[k] = -ab_y[k] * fg_127[k]
                       + fh_179[k];

            t_188[k] = -ab_y[k] * fg_128[k]
                       + fh_180[k];

            t_189[k] = -ab_y[k] * fg_129[k]
                       + fh_181[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, fg_130, fg_131, fg_132, \
                         fg_133, fg_134, fh_183, fh_184, fh_185, fh_186, \
                         fh_187 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_y[k] * fg_130[k]
                       + fh_183[k];

            t_191[k] = -ab_y[k] * fg_131[k]
                       + fh_184[k];

            t_192[k] = -ab_y[k] * fg_132[k]
                       + fh_185[k];

            t_193[k] = -ab_y[k] * fg_133[k]
                       + fh_186[k];

            t_194[k] = -ab_y[k] * fg_134[k]
                       + fh_187[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_y, fg_135, fg_136, fg_137, \
                         fg_138, fg_139, fh_190, fh_192, fh_193, fh_195, \
                         fh_196 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_y[k] * fg_135[k]
                       + fh_190[k];

            t_196[k] = -ab_y[k] * fg_136[k]
                       + fh_192[k];

            t_197[k] = -ab_y[k] * fg_137[k]
                       + fh_193[k];

            t_198[k] = -ab_y[k] * fg_138[k]
                       + fh_195[k];

            t_199[k] = -ab_y[k] * fg_139[k]
                       + fh_196[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_y, fg_140, fg_141, fg_142, \
                         fg_143, fg_144, fh_197, fh_199, fh_200, fh_201, \
                         fh_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_y[k] * fg_140[k]
                       + fh_197[k];

            t_201[k] = -ab_y[k] * fg_141[k]
                       + fh_199[k];

            t_202[k] = -ab_y[k] * fg_142[k]
                       + fh_200[k];

            t_203[k] = -ab_y[k] * fg_143[k]
                       + fh_201[k];

            t_204[k] = -ab_y[k] * fg_144[k]
                       + fh_202[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, fg_145, fg_146, fg_147, \
                         fg_148, fg_149, fh_204, fh_205, fh_206, fh_207, \
                         fh_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_y[k] * fg_145[k]
                       + fh_204[k];

            t_206[k] = -ab_y[k] * fg_146[k]
                       + fh_205[k];

            t_207[k] = -ab_y[k] * fg_147[k]
                       + fh_206[k];

            t_208[k] = -ab_y[k] * fg_148[k]
                       + fh_207[k];

            t_209[k] = -ab_y[k] * fg_149[k]
                       + fh_208[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_z, fg_135, fg_136, fg_137, \
                         fg_138, fg_139, fh_191, fh_193, fh_194, fh_196, \
                         fh_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_z[k] * fg_135[k]
                       + fh_191[k];

            t_211[k] = -ab_z[k] * fg_136[k]
                       + fh_193[k];

            t_212[k] = -ab_z[k] * fg_137[k]
                       + fh_194[k];

            t_213[k] = -ab_z[k] * fg_138[k]
                       + fh_196[k];

            t_214[k] = -ab_z[k] * fg_139[k]
                       + fh_197[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_z, fg_140, fg_141, fg_142, \
                         fg_143, fg_144, fh_198, fh_200, fh_201, fh_202, \
                         fh_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_z[k] * fg_140[k]
                       + fh_198[k];

            t_216[k] = -ab_z[k] * fg_141[k]
                       + fh_200[k];

            t_217[k] = -ab_z[k] * fg_142[k]
                       + fh_201[k];

            t_218[k] = -ab_z[k] * fg_143[k]
                       + fh_202[k];

            t_219[k] = -ab_z[k] * fg_144[k]
                       + fh_203[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_z, fg_145, fg_146, fg_147, \
                         fg_148, fg_149, fh_205, fh_206, fh_207, fh_208, \
                         fh_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_z[k] * fg_145[k]
                       + fh_205[k];

            t_221[k] = -ab_z[k] * fg_146[k]
                       + fh_206[k];

            t_222[k] = -ab_z[k] * fg_147[k]
                       + fh_207[k];

            t_223[k] = -ab_z[k] * fg_148[k]
                       + fh_208[k];

            t_224[k] = -ab_z[k] * fg_149[k]
                       + fh_209[k];
        }
    }
}

auto
compute_hrr_gg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fg, const size_t fh, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_gg_piece0(buffer, coordinates, target, fg, fh, ncomps, nmax);

    compute_hrr_gg_piece1(buffer, coordinates, target, fg, fh, ncomps, nmax);
}

}  // namespace simdtrf
