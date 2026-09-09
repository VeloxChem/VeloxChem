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


#include "SimdTransferLF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_lf_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ld, const size_t md,
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

        const auto *ld_0 = buffer.data(ld + 0 * ncomps + c);
        const auto *ld_1 = buffer.data(ld + 1 * ncomps + c);
        const auto *ld_2 = buffer.data(ld + 2 * ncomps + c);
        const auto *ld_3 = buffer.data(ld + 3 * ncomps + c);
        const auto *ld_4 = buffer.data(ld + 4 * ncomps + c);
        const auto *ld_5 = buffer.data(ld + 5 * ncomps + c);
        const auto *ld_6 = buffer.data(ld + 6 * ncomps + c);
        const auto *ld_7 = buffer.data(ld + 7 * ncomps + c);
        const auto *ld_8 = buffer.data(ld + 8 * ncomps + c);
        const auto *ld_9 = buffer.data(ld + 9 * ncomps + c);
        const auto *ld_10 = buffer.data(ld + 10 * ncomps + c);
        const auto *ld_11 = buffer.data(ld + 11 * ncomps + c);
        const auto *ld_12 = buffer.data(ld + 12 * ncomps + c);
        const auto *ld_13 = buffer.data(ld + 13 * ncomps + c);
        const auto *ld_14 = buffer.data(ld + 14 * ncomps + c);
        const auto *ld_15 = buffer.data(ld + 15 * ncomps + c);
        const auto *ld_16 = buffer.data(ld + 16 * ncomps + c);
        const auto *ld_17 = buffer.data(ld + 17 * ncomps + c);
        const auto *ld_18 = buffer.data(ld + 18 * ncomps + c);
        const auto *ld_19 = buffer.data(ld + 19 * ncomps + c);
        const auto *ld_20 = buffer.data(ld + 20 * ncomps + c);
        const auto *ld_21 = buffer.data(ld + 21 * ncomps + c);
        const auto *ld_22 = buffer.data(ld + 22 * ncomps + c);
        const auto *ld_23 = buffer.data(ld + 23 * ncomps + c);
        const auto *ld_24 = buffer.data(ld + 24 * ncomps + c);
        const auto *ld_25 = buffer.data(ld + 25 * ncomps + c);
        const auto *ld_26 = buffer.data(ld + 26 * ncomps + c);
        const auto *ld_27 = buffer.data(ld + 27 * ncomps + c);
        const auto *ld_28 = buffer.data(ld + 28 * ncomps + c);
        const auto *ld_29 = buffer.data(ld + 29 * ncomps + c);
        const auto *ld_30 = buffer.data(ld + 30 * ncomps + c);
        const auto *ld_31 = buffer.data(ld + 31 * ncomps + c);
        const auto *ld_32 = buffer.data(ld + 32 * ncomps + c);
        const auto *ld_33 = buffer.data(ld + 33 * ncomps + c);
        const auto *ld_34 = buffer.data(ld + 34 * ncomps + c);
        const auto *ld_35 = buffer.data(ld + 35 * ncomps + c);
        const auto *ld_36 = buffer.data(ld + 36 * ncomps + c);
        const auto *ld_37 = buffer.data(ld + 37 * ncomps + c);
        const auto *ld_38 = buffer.data(ld + 38 * ncomps + c);
        const auto *ld_39 = buffer.data(ld + 39 * ncomps + c);
        const auto *ld_40 = buffer.data(ld + 40 * ncomps + c);
        const auto *ld_41 = buffer.data(ld + 41 * ncomps + c);
        const auto *ld_42 = buffer.data(ld + 42 * ncomps + c);
        const auto *ld_43 = buffer.data(ld + 43 * ncomps + c);
        const auto *ld_44 = buffer.data(ld + 44 * ncomps + c);
        const auto *ld_45 = buffer.data(ld + 45 * ncomps + c);
        const auto *ld_46 = buffer.data(ld + 46 * ncomps + c);
        const auto *ld_47 = buffer.data(ld + 47 * ncomps + c);
        const auto *ld_48 = buffer.data(ld + 48 * ncomps + c);
        const auto *ld_49 = buffer.data(ld + 49 * ncomps + c);
        const auto *ld_50 = buffer.data(ld + 50 * ncomps + c);
        const auto *ld_51 = buffer.data(ld + 51 * ncomps + c);
        const auto *ld_52 = buffer.data(ld + 52 * ncomps + c);
        const auto *ld_53 = buffer.data(ld + 53 * ncomps + c);
        const auto *ld_54 = buffer.data(ld + 54 * ncomps + c);
        const auto *ld_55 = buffer.data(ld + 55 * ncomps + c);
        const auto *ld_56 = buffer.data(ld + 56 * ncomps + c);
        const auto *ld_57 = buffer.data(ld + 57 * ncomps + c);
        const auto *ld_58 = buffer.data(ld + 58 * ncomps + c);
        const auto *ld_59 = buffer.data(ld + 59 * ncomps + c);
        const auto *ld_60 = buffer.data(ld + 60 * ncomps + c);
        const auto *ld_61 = buffer.data(ld + 61 * ncomps + c);
        const auto *ld_62 = buffer.data(ld + 62 * ncomps + c);
        const auto *ld_63 = buffer.data(ld + 63 * ncomps + c);
        const auto *ld_64 = buffer.data(ld + 64 * ncomps + c);
        const auto *ld_65 = buffer.data(ld + 65 * ncomps + c);
        const auto *ld_66 = buffer.data(ld + 66 * ncomps + c);
        const auto *ld_67 = buffer.data(ld + 67 * ncomps + c);
        const auto *ld_68 = buffer.data(ld + 68 * ncomps + c);
        const auto *ld_69 = buffer.data(ld + 69 * ncomps + c);
        const auto *ld_70 = buffer.data(ld + 70 * ncomps + c);
        const auto *ld_71 = buffer.data(ld + 71 * ncomps + c);
        const auto *ld_72 = buffer.data(ld + 72 * ncomps + c);
        const auto *ld_73 = buffer.data(ld + 73 * ncomps + c);
        const auto *ld_74 = buffer.data(ld + 74 * ncomps + c);
        const auto *ld_75 = buffer.data(ld + 75 * ncomps + c);
        const auto *ld_76 = buffer.data(ld + 76 * ncomps + c);
        const auto *ld_77 = buffer.data(ld + 77 * ncomps + c);
        const auto *ld_78 = buffer.data(ld + 78 * ncomps + c);
        const auto *ld_79 = buffer.data(ld + 79 * ncomps + c);
        const auto *ld_80 = buffer.data(ld + 80 * ncomps + c);
        const auto *ld_81 = buffer.data(ld + 81 * ncomps + c);
        const auto *ld_82 = buffer.data(ld + 82 * ncomps + c);
        const auto *ld_83 = buffer.data(ld + 83 * ncomps + c);
        const auto *ld_84 = buffer.data(ld + 84 * ncomps + c);
        const auto *ld_85 = buffer.data(ld + 85 * ncomps + c);
        const auto *ld_86 = buffer.data(ld + 86 * ncomps + c);
        const auto *ld_87 = buffer.data(ld + 87 * ncomps + c);
        const auto *ld_88 = buffer.data(ld + 88 * ncomps + c);

        const auto *md_0 = buffer.data(md + 0 * ncomps + c);
        const auto *md_1 = buffer.data(md + 1 * ncomps + c);
        const auto *md_2 = buffer.data(md + 2 * ncomps + c);
        const auto *md_3 = buffer.data(md + 3 * ncomps + c);
        const auto *md_4 = buffer.data(md + 4 * ncomps + c);
        const auto *md_5 = buffer.data(md + 5 * ncomps + c);
        const auto *md_6 = buffer.data(md + 6 * ncomps + c);
        const auto *md_7 = buffer.data(md + 7 * ncomps + c);
        const auto *md_8 = buffer.data(md + 8 * ncomps + c);
        const auto *md_9 = buffer.data(md + 9 * ncomps + c);
        const auto *md_10 = buffer.data(md + 10 * ncomps + c);
        const auto *md_11 = buffer.data(md + 11 * ncomps + c);
        const auto *md_12 = buffer.data(md + 12 * ncomps + c);
        const auto *md_13 = buffer.data(md + 13 * ncomps + c);
        const auto *md_14 = buffer.data(md + 14 * ncomps + c);
        const auto *md_15 = buffer.data(md + 15 * ncomps + c);
        const auto *md_16 = buffer.data(md + 16 * ncomps + c);
        const auto *md_17 = buffer.data(md + 17 * ncomps + c);
        const auto *md_18 = buffer.data(md + 18 * ncomps + c);
        const auto *md_19 = buffer.data(md + 19 * ncomps + c);
        const auto *md_20 = buffer.data(md + 20 * ncomps + c);
        const auto *md_21 = buffer.data(md + 21 * ncomps + c);
        const auto *md_22 = buffer.data(md + 22 * ncomps + c);
        const auto *md_23 = buffer.data(md + 23 * ncomps + c);
        const auto *md_24 = buffer.data(md + 24 * ncomps + c);
        const auto *md_25 = buffer.data(md + 25 * ncomps + c);
        const auto *md_26 = buffer.data(md + 26 * ncomps + c);
        const auto *md_27 = buffer.data(md + 27 * ncomps + c);
        const auto *md_28 = buffer.data(md + 28 * ncomps + c);
        const auto *md_29 = buffer.data(md + 29 * ncomps + c);
        const auto *md_30 = buffer.data(md + 30 * ncomps + c);
        const auto *md_31 = buffer.data(md + 31 * ncomps + c);
        const auto *md_32 = buffer.data(md + 32 * ncomps + c);
        const auto *md_33 = buffer.data(md + 33 * ncomps + c);
        const auto *md_34 = buffer.data(md + 34 * ncomps + c);
        const auto *md_35 = buffer.data(md + 35 * ncomps + c);
        const auto *md_36 = buffer.data(md + 36 * ncomps + c);
        const auto *md_37 = buffer.data(md + 37 * ncomps + c);
        const auto *md_38 = buffer.data(md + 38 * ncomps + c);
        const auto *md_39 = buffer.data(md + 39 * ncomps + c);
        const auto *md_40 = buffer.data(md + 40 * ncomps + c);
        const auto *md_41 = buffer.data(md + 41 * ncomps + c);
        const auto *md_42 = buffer.data(md + 42 * ncomps + c);
        const auto *md_43 = buffer.data(md + 43 * ncomps + c);
        const auto *md_44 = buffer.data(md + 44 * ncomps + c);
        const auto *md_45 = buffer.data(md + 45 * ncomps + c);
        const auto *md_46 = buffer.data(md + 46 * ncomps + c);
        const auto *md_47 = buffer.data(md + 47 * ncomps + c);
        const auto *md_48 = buffer.data(md + 48 * ncomps + c);
        const auto *md_49 = buffer.data(md + 49 * ncomps + c);
        const auto *md_50 = buffer.data(md + 50 * ncomps + c);
        const auto *md_51 = buffer.data(md + 51 * ncomps + c);
        const auto *md_52 = buffer.data(md + 52 * ncomps + c);
        const auto *md_53 = buffer.data(md + 53 * ncomps + c);
        const auto *md_54 = buffer.data(md + 54 * ncomps + c);
        const auto *md_55 = buffer.data(md + 55 * ncomps + c);
        const auto *md_56 = buffer.data(md + 56 * ncomps + c);
        const auto *md_57 = buffer.data(md + 57 * ncomps + c);
        const auto *md_58 = buffer.data(md + 58 * ncomps + c);
        const auto *md_59 = buffer.data(md + 59 * ncomps + c);
        const auto *md_60 = buffer.data(md + 60 * ncomps + c);
        const auto *md_61 = buffer.data(md + 61 * ncomps + c);
        const auto *md_62 = buffer.data(md + 62 * ncomps + c);
        const auto *md_63 = buffer.data(md + 63 * ncomps + c);
        const auto *md_64 = buffer.data(md + 64 * ncomps + c);
        const auto *md_65 = buffer.data(md + 65 * ncomps + c);
        const auto *md_66 = buffer.data(md + 66 * ncomps + c);
        const auto *md_67 = buffer.data(md + 67 * ncomps + c);
        const auto *md_68 = buffer.data(md + 68 * ncomps + c);
        const auto *md_69 = buffer.data(md + 69 * ncomps + c);
        const auto *md_70 = buffer.data(md + 70 * ncomps + c);
        const auto *md_71 = buffer.data(md + 71 * ncomps + c);
        const auto *md_72 = buffer.data(md + 72 * ncomps + c);
        const auto *md_73 = buffer.data(md + 73 * ncomps + c);
        const auto *md_74 = buffer.data(md + 74 * ncomps + c);
        const auto *md_75 = buffer.data(md + 75 * ncomps + c);
        const auto *md_76 = buffer.data(md + 76 * ncomps + c);
        const auto *md_77 = buffer.data(md + 77 * ncomps + c);
        const auto *md_78 = buffer.data(md + 78 * ncomps + c);
        const auto *md_79 = buffer.data(md + 79 * ncomps + c);
        const auto *md_80 = buffer.data(md + 80 * ncomps + c);
        const auto *md_81 = buffer.data(md + 81 * ncomps + c);
        const auto *md_82 = buffer.data(md + 82 * ncomps + c);
        const auto *md_83 = buffer.data(md + 83 * ncomps + c);
        const auto *md_84 = buffer.data(md + 84 * ncomps + c);
        const auto *md_85 = buffer.data(md + 85 * ncomps + c);
        const auto *md_86 = buffer.data(md + 86 * ncomps + c);
        const auto *md_87 = buffer.data(md + 87 * ncomps + c);
        const auto *md_88 = buffer.data(md + 88 * ncomps + c);
        const auto *md_89 = buffer.data(md + 89 * ncomps + c);
        const auto *md_93 = buffer.data(md + 93 * ncomps + c);
        const auto *md_94 = buffer.data(md + 94 * ncomps + c);
        const auto *md_95 = buffer.data(md + 95 * ncomps + c);
        const auto *md_99 = buffer.data(md + 99 * ncomps + c);
        const auto *md_100 = buffer.data(md + 100 * ncomps + c);
        const auto *md_101 = buffer.data(md + 101 * ncomps + c);
        const auto *md_105 = buffer.data(md + 105 * ncomps + c);
        const auto *md_106 = buffer.data(md + 106 * ncomps + c);
        const auto *md_107 = buffer.data(md + 107 * ncomps + c);
        const auto *md_111 = buffer.data(md + 111 * ncomps + c);
        const auto *md_112 = buffer.data(md + 112 * ncomps + c);
        const auto *md_113 = buffer.data(md + 113 * ncomps + c);
        const auto *md_119 = buffer.data(md + 119 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ld_0, ld_1, ld_2, ld_3, ld_4, md_0, \
                         md_1, md_2, md_3, md_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ld_0[k]
                     + md_0[k];

            t_1[k] = ab_x[k] * ld_1[k]
                     + md_1[k];

            t_2[k] = ab_x[k] * ld_2[k]
                     + md_2[k];

            t_3[k] = ab_x[k] * ld_3[k]
                     + md_3[k];

            t_4[k] = ab_x[k] * ld_4[k]
                     + md_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, ld_3, ld_4, ld_5, md_5, \
                         md_9, md_10, md_11, md_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * ld_5[k]
                     + md_5[k];

            t_6[k] = ab_y[k] * ld_3[k]
                     + md_9[k];

            t_7[k] = ab_y[k] * ld_4[k]
                     + md_10[k];

            t_8[k] = ab_y[k] * ld_5[k]
                     + md_11[k];

            t_9[k] = ab_z[k] * ld_5[k]
                     + md_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ld_6, ld_7, ld_8, ld_9, ld_10, \
                         md_6, md_7, md_8, md_9, md_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * ld_6[k]
                      + md_6[k];

            t_11[k] = ab_x[k] * ld_7[k]
                      + md_7[k];

            t_12[k] = ab_x[k] * ld_8[k]
                      + md_8[k];

            t_13[k] = ab_x[k] * ld_9[k]
                      + md_9[k];

            t_14[k] = ab_x[k] * ld_10[k]
                      + md_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, ld_9, ld_10, ld_11, \
                         md_11, md_21, md_22, md_23, md_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * ld_11[k]
                      + md_11[k];

            t_16[k] = ab_y[k] * ld_9[k]
                      + md_21[k];

            t_17[k] = ab_y[k] * ld_10[k]
                      + md_22[k];

            t_18[k] = ab_y[k] * ld_11[k]
                      + md_23[k];

            t_19[k] = ab_z[k] * ld_11[k]
                      + md_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, ld_12, ld_13, ld_14, ld_15, \
                         ld_16, md_12, md_13, md_14, md_15, md_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * ld_12[k]
                      + md_12[k];

            t_21[k] = ab_x[k] * ld_13[k]
                      + md_13[k];

            t_22[k] = ab_x[k] * ld_14[k]
                      + md_14[k];

            t_23[k] = ab_x[k] * ld_15[k]
                      + md_15[k];

            t_24[k] = ab_x[k] * ld_16[k]
                      + md_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, ld_15, ld_16, ld_17, \
                         md_17, md_27, md_28, md_29, md_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * ld_17[k]
                      + md_17[k];

            t_26[k] = ab_y[k] * ld_15[k]
                      + md_27[k];

            t_27[k] = ab_y[k] * ld_16[k]
                      + md_28[k];

            t_28[k] = ab_y[k] * ld_17[k]
                      + md_29[k];

            t_29[k] = ab_z[k] * ld_17[k]
                      + md_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, ld_18, ld_19, ld_20, ld_21, \
                         ld_22, md_18, md_19, md_20, md_21, md_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * ld_18[k]
                      + md_18[k];

            t_31[k] = ab_x[k] * ld_19[k]
                      + md_19[k];

            t_32[k] = ab_x[k] * ld_20[k]
                      + md_20[k];

            t_33[k] = ab_x[k] * ld_21[k]
                      + md_21[k];

            t_34[k] = ab_x[k] * ld_22[k]
                      + md_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, ld_21, ld_22, ld_23, \
                         md_23, md_39, md_40, md_41, md_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * ld_23[k]
                      + md_23[k];

            t_36[k] = ab_y[k] * ld_21[k]
                      + md_39[k];

            t_37[k] = ab_y[k] * ld_22[k]
                      + md_40[k];

            t_38[k] = ab_y[k] * ld_23[k]
                      + md_41[k];

            t_39[k] = ab_z[k] * ld_23[k]
                      + md_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, ld_24, ld_25, ld_26, ld_27, \
                         ld_28, md_24, md_25, md_26, md_27, md_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * ld_24[k]
                      + md_24[k];

            t_41[k] = ab_x[k] * ld_25[k]
                      + md_25[k];

            t_42[k] = ab_x[k] * ld_26[k]
                      + md_26[k];

            t_43[k] = ab_x[k] * ld_27[k]
                      + md_27[k];

            t_44[k] = ab_x[k] * ld_28[k]
                      + md_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, ld_27, ld_28, ld_29, \
                         md_29, md_45, md_46, md_47, md_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * ld_29[k]
                      + md_29[k];

            t_46[k] = ab_y[k] * ld_27[k]
                      + md_45[k];

            t_47[k] = ab_y[k] * ld_28[k]
                      + md_46[k];

            t_48[k] = ab_y[k] * ld_29[k]
                      + md_47[k];

            t_49[k] = ab_z[k] * ld_29[k]
                      + md_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, ld_30, ld_31, ld_32, ld_33, \
                         ld_34, md_30, md_31, md_32, md_33, md_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * ld_30[k]
                      + md_30[k];

            t_51[k] = ab_x[k] * ld_31[k]
                      + md_31[k];

            t_52[k] = ab_x[k] * ld_32[k]
                      + md_32[k];

            t_53[k] = ab_x[k] * ld_33[k]
                      + md_33[k];

            t_54[k] = ab_x[k] * ld_34[k]
                      + md_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ld_33, ld_34, ld_35, \
                         md_35, md_51, md_52, md_53, md_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * ld_35[k]
                      + md_35[k];

            t_56[k] = ab_y[k] * ld_33[k]
                      + md_51[k];

            t_57[k] = ab_y[k] * ld_34[k]
                      + md_52[k];

            t_58[k] = ab_y[k] * ld_35[k]
                      + md_53[k];

            t_59[k] = ab_z[k] * ld_35[k]
                      + md_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ld_36, ld_37, ld_38, ld_39, \
                         ld_40, md_36, md_37, md_38, md_39, md_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ld_36[k]
                      + md_36[k];

            t_61[k] = ab_x[k] * ld_37[k]
                      + md_37[k];

            t_62[k] = ab_x[k] * ld_38[k]
                      + md_38[k];

            t_63[k] = ab_x[k] * ld_39[k]
                      + md_39[k];

            t_64[k] = ab_x[k] * ld_40[k]
                      + md_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, ld_39, ld_40, ld_41, \
                         md_41, md_63, md_64, md_65, md_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * ld_41[k]
                      + md_41[k];

            t_66[k] = ab_y[k] * ld_39[k]
                      + md_63[k];

            t_67[k] = ab_y[k] * ld_40[k]
                      + md_64[k];

            t_68[k] = ab_y[k] * ld_41[k]
                      + md_65[k];

            t_69[k] = ab_z[k] * ld_41[k]
                      + md_71[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, ld_42, ld_43, ld_44, ld_45, \
                         ld_46, md_42, md_43, md_44, md_45, md_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_x[k] * ld_42[k]
                      + md_42[k];

            t_71[k] = ab_x[k] * ld_43[k]
                      + md_43[k];

            t_72[k] = ab_x[k] * ld_44[k]
                      + md_44[k];

            t_73[k] = ab_x[k] * ld_45[k]
                      + md_45[k];

            t_74[k] = ab_x[k] * ld_46[k]
                      + md_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, ld_45, ld_46, ld_47, \
                         md_47, md_69, md_70, md_71, md_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * ld_47[k]
                      + md_47[k];

            t_76[k] = ab_y[k] * ld_45[k]
                      + md_69[k];

            t_77[k] = ab_y[k] * ld_46[k]
                      + md_70[k];

            t_78[k] = ab_y[k] * ld_47[k]
                      + md_71[k];

            t_79[k] = ab_z[k] * ld_47[k]
                      + md_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, ld_48, ld_49, ld_50, ld_51, \
                         ld_52, md_48, md_49, md_50, md_51, md_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * ld_48[k]
                      + md_48[k];

            t_81[k] = ab_x[k] * ld_49[k]
                      + md_49[k];

            t_82[k] = ab_x[k] * ld_50[k]
                      + md_50[k];

            t_83[k] = ab_x[k] * ld_51[k]
                      + md_51[k];

            t_84[k] = ab_x[k] * ld_52[k]
                      + md_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, ld_51, ld_52, ld_53, \
                         md_53, md_75, md_76, md_77, md_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * ld_53[k]
                      + md_53[k];

            t_86[k] = ab_y[k] * ld_51[k]
                      + md_75[k];

            t_87[k] = ab_y[k] * ld_52[k]
                      + md_76[k];

            t_88[k] = ab_y[k] * ld_53[k]
                      + md_77[k];

            t_89[k] = ab_z[k] * ld_53[k]
                      + md_83[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ld_54, ld_55, ld_56, ld_57, \
                         ld_58, md_54, md_55, md_56, md_57, md_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * ld_54[k]
                      + md_54[k];

            t_91[k] = ab_x[k] * ld_55[k]
                      + md_55[k];

            t_92[k] = ab_x[k] * ld_56[k]
                      + md_56[k];

            t_93[k] = ab_x[k] * ld_57[k]
                      + md_57[k];

            t_94[k] = ab_x[k] * ld_58[k]
                      + md_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, ld_57, ld_58, ld_59, \
                         md_59, md_81, md_82, md_83, md_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * ld_59[k]
                      + md_59[k];

            t_96[k] = ab_y[k] * ld_57[k]
                      + md_81[k];

            t_97[k] = ab_y[k] * ld_58[k]
                      + md_82[k];

            t_98[k] = ab_y[k] * ld_59[k]
                      + md_83[k];

            t_99[k] = ab_z[k] * ld_59[k]
                      + md_89[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, ld_60, ld_61, ld_62, ld_63, \
                         ld_64, md_60, md_61, md_62, md_63, md_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_x[k] * ld_60[k]
                       + md_60[k];

            t_101[k] = ab_x[k] * ld_61[k]
                       + md_61[k];

            t_102[k] = ab_x[k] * ld_62[k]
                       + md_62[k];

            t_103[k] = ab_x[k] * ld_63[k]
                       + md_63[k];

            t_104[k] = ab_x[k] * ld_64[k]
                       + md_64[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, ld_63, ld_64, \
                         ld_65, md_65, md_93, md_94, md_95, md_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * ld_65[k]
                       + md_65[k];

            t_106[k] = ab_y[k] * ld_63[k]
                       + md_93[k];

            t_107[k] = ab_y[k] * ld_64[k]
                       + md_94[k];

            t_108[k] = ab_y[k] * ld_65[k]
                       + md_95[k];

            t_109[k] = ab_z[k] * ld_65[k]
                       + md_101[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, ld_66, ld_67, ld_68, ld_69, \
                         ld_70, md_66, md_67, md_68, md_69, md_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * ld_66[k]
                       + md_66[k];

            t_111[k] = ab_x[k] * ld_67[k]
                       + md_67[k];

            t_112[k] = ab_x[k] * ld_68[k]
                       + md_68[k];

            t_113[k] = ab_x[k] * ld_69[k]
                       + md_69[k];

            t_114[k] = ab_x[k] * ld_70[k]
                       + md_70[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, ld_69, ld_70, \
                         ld_71, md_71, md_99, md_100, md_101, md_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_x[k] * ld_71[k]
                       + md_71[k];

            t_116[k] = ab_y[k] * ld_69[k]
                       + md_99[k];

            t_117[k] = ab_y[k] * ld_70[k]
                       + md_100[k];

            t_118[k] = ab_y[k] * ld_71[k]
                       + md_101[k];

            t_119[k] = ab_z[k] * ld_71[k]
                       + md_107[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, ld_72, ld_73, ld_74, ld_75, \
                         ld_76, md_72, md_73, md_74, md_75, md_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * ld_72[k]
                       + md_72[k];

            t_121[k] = ab_x[k] * ld_73[k]
                       + md_73[k];

            t_122[k] = ab_x[k] * ld_74[k]
                       + md_74[k];

            t_123[k] = ab_x[k] * ld_75[k]
                       + md_75[k];

            t_124[k] = ab_x[k] * ld_76[k]
                       + md_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, ld_75, ld_76, \
                         ld_77, md_77, md_105, md_106, md_107, md_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * ld_77[k]
                       + md_77[k];

            t_126[k] = ab_y[k] * ld_75[k]
                       + md_105[k];

            t_127[k] = ab_y[k] * ld_76[k]
                       + md_106[k];

            t_128[k] = ab_y[k] * ld_77[k]
                       + md_107[k];

            t_129[k] = ab_z[k] * ld_77[k]
                       + md_113[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, ld_78, ld_79, ld_80, ld_81, \
                         ld_82, md_78, md_79, md_80, md_81, md_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_x[k] * ld_78[k]
                       + md_78[k];

            t_131[k] = ab_x[k] * ld_79[k]
                       + md_79[k];

            t_132[k] = ab_x[k] * ld_80[k]
                       + md_80[k];

            t_133[k] = ab_x[k] * ld_81[k]
                       + md_81[k];

            t_134[k] = ab_x[k] * ld_82[k]
                       + md_82[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, ld_81, ld_82, \
                         ld_83, md_83, md_111, md_112, md_113, md_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * ld_83[k]
                       + md_83[k];

            t_136[k] = ab_y[k] * ld_81[k]
                       + md_111[k];

            t_137[k] = ab_y[k] * ld_82[k]
                       + md_112[k];

            t_138[k] = ab_y[k] * ld_83[k]
                       + md_113[k];

            t_139[k] = ab_z[k] * ld_83[k]
                       + md_119[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, ld_84, ld_85, ld_86, ld_87, \
                         ld_88, md_84, md_85, md_86, md_87, md_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * ld_84[k]
                       + md_84[k];

            t_141[k] = ab_x[k] * ld_85[k]
                       + md_85[k];

            t_142[k] = ab_x[k] * ld_86[k]
                       + md_86[k];

            t_143[k] = ab_x[k] * ld_87[k]
                       + md_87[k];

            t_144[k] = ab_x[k] * ld_88[k]
                       + md_88[k];
        }
    }
}

static auto
compute_hrr_lf_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ld, const size_t md,
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

        const auto *ld_87 = buffer.data(ld + 87 * ncomps + c);
        const auto *ld_88 = buffer.data(ld + 88 * ncomps + c);
        const auto *ld_89 = buffer.data(ld + 89 * ncomps + c);
        const auto *ld_90 = buffer.data(ld + 90 * ncomps + c);
        const auto *ld_91 = buffer.data(ld + 91 * ncomps + c);
        const auto *ld_92 = buffer.data(ld + 92 * ncomps + c);
        const auto *ld_93 = buffer.data(ld + 93 * ncomps + c);
        const auto *ld_94 = buffer.data(ld + 94 * ncomps + c);
        const auto *ld_95 = buffer.data(ld + 95 * ncomps + c);
        const auto *ld_96 = buffer.data(ld + 96 * ncomps + c);
        const auto *ld_97 = buffer.data(ld + 97 * ncomps + c);
        const auto *ld_98 = buffer.data(ld + 98 * ncomps + c);
        const auto *ld_99 = buffer.data(ld + 99 * ncomps + c);
        const auto *ld_100 = buffer.data(ld + 100 * ncomps + c);
        const auto *ld_101 = buffer.data(ld + 101 * ncomps + c);
        const auto *ld_102 = buffer.data(ld + 102 * ncomps + c);
        const auto *ld_103 = buffer.data(ld + 103 * ncomps + c);
        const auto *ld_104 = buffer.data(ld + 104 * ncomps + c);
        const auto *ld_105 = buffer.data(ld + 105 * ncomps + c);
        const auto *ld_106 = buffer.data(ld + 106 * ncomps + c);
        const auto *ld_107 = buffer.data(ld + 107 * ncomps + c);
        const auto *ld_108 = buffer.data(ld + 108 * ncomps + c);
        const auto *ld_109 = buffer.data(ld + 109 * ncomps + c);
        const auto *ld_110 = buffer.data(ld + 110 * ncomps + c);
        const auto *ld_111 = buffer.data(ld + 111 * ncomps + c);
        const auto *ld_112 = buffer.data(ld + 112 * ncomps + c);
        const auto *ld_113 = buffer.data(ld + 113 * ncomps + c);
        const auto *ld_114 = buffer.data(ld + 114 * ncomps + c);
        const auto *ld_115 = buffer.data(ld + 115 * ncomps + c);
        const auto *ld_116 = buffer.data(ld + 116 * ncomps + c);
        const auto *ld_117 = buffer.data(ld + 117 * ncomps + c);
        const auto *ld_118 = buffer.data(ld + 118 * ncomps + c);
        const auto *ld_119 = buffer.data(ld + 119 * ncomps + c);
        const auto *ld_120 = buffer.data(ld + 120 * ncomps + c);
        const auto *ld_121 = buffer.data(ld + 121 * ncomps + c);
        const auto *ld_122 = buffer.data(ld + 122 * ncomps + c);
        const auto *ld_123 = buffer.data(ld + 123 * ncomps + c);
        const auto *ld_124 = buffer.data(ld + 124 * ncomps + c);
        const auto *ld_125 = buffer.data(ld + 125 * ncomps + c);
        const auto *ld_126 = buffer.data(ld + 126 * ncomps + c);
        const auto *ld_127 = buffer.data(ld + 127 * ncomps + c);
        const auto *ld_128 = buffer.data(ld + 128 * ncomps + c);
        const auto *ld_129 = buffer.data(ld + 129 * ncomps + c);
        const auto *ld_130 = buffer.data(ld + 130 * ncomps + c);
        const auto *ld_131 = buffer.data(ld + 131 * ncomps + c);
        const auto *ld_132 = buffer.data(ld + 132 * ncomps + c);
        const auto *ld_133 = buffer.data(ld + 133 * ncomps + c);
        const auto *ld_134 = buffer.data(ld + 134 * ncomps + c);
        const auto *ld_135 = buffer.data(ld + 135 * ncomps + c);
        const auto *ld_136 = buffer.data(ld + 136 * ncomps + c);
        const auto *ld_137 = buffer.data(ld + 137 * ncomps + c);
        const auto *ld_138 = buffer.data(ld + 138 * ncomps + c);
        const auto *ld_139 = buffer.data(ld + 139 * ncomps + c);
        const auto *ld_140 = buffer.data(ld + 140 * ncomps + c);
        const auto *ld_141 = buffer.data(ld + 141 * ncomps + c);
        const auto *ld_142 = buffer.data(ld + 142 * ncomps + c);
        const auto *ld_143 = buffer.data(ld + 143 * ncomps + c);
        const auto *ld_144 = buffer.data(ld + 144 * ncomps + c);
        const auto *ld_145 = buffer.data(ld + 145 * ncomps + c);
        const auto *ld_146 = buffer.data(ld + 146 * ncomps + c);
        const auto *ld_147 = buffer.data(ld + 147 * ncomps + c);
        const auto *ld_148 = buffer.data(ld + 148 * ncomps + c);
        const auto *ld_149 = buffer.data(ld + 149 * ncomps + c);
        const auto *ld_150 = buffer.data(ld + 150 * ncomps + c);
        const auto *ld_151 = buffer.data(ld + 151 * ncomps + c);
        const auto *ld_152 = buffer.data(ld + 152 * ncomps + c);
        const auto *ld_153 = buffer.data(ld + 153 * ncomps + c);
        const auto *ld_154 = buffer.data(ld + 154 * ncomps + c);
        const auto *ld_155 = buffer.data(ld + 155 * ncomps + c);
        const auto *ld_156 = buffer.data(ld + 156 * ncomps + c);
        const auto *ld_157 = buffer.data(ld + 157 * ncomps + c);
        const auto *ld_158 = buffer.data(ld + 158 * ncomps + c);
        const auto *ld_159 = buffer.data(ld + 159 * ncomps + c);
        const auto *ld_160 = buffer.data(ld + 160 * ncomps + c);
        const auto *ld_161 = buffer.data(ld + 161 * ncomps + c);
        const auto *ld_162 = buffer.data(ld + 162 * ncomps + c);
        const auto *ld_163 = buffer.data(ld + 163 * ncomps + c);
        const auto *ld_164 = buffer.data(ld + 164 * ncomps + c);
        const auto *ld_165 = buffer.data(ld + 165 * ncomps + c);
        const auto *ld_166 = buffer.data(ld + 166 * ncomps + c);
        const auto *ld_167 = buffer.data(ld + 167 * ncomps + c);
        const auto *ld_168 = buffer.data(ld + 168 * ncomps + c);
        const auto *ld_169 = buffer.data(ld + 169 * ncomps + c);
        const auto *ld_170 = buffer.data(ld + 170 * ncomps + c);
        const auto *ld_171 = buffer.data(ld + 171 * ncomps + c);
        const auto *ld_172 = buffer.data(ld + 172 * ncomps + c);
        const auto *ld_173 = buffer.data(ld + 173 * ncomps + c);

        const auto *md_89 = buffer.data(md + 89 * ncomps + c);
        const auto *md_90 = buffer.data(md + 90 * ncomps + c);
        const auto *md_91 = buffer.data(md + 91 * ncomps + c);
        const auto *md_92 = buffer.data(md + 92 * ncomps + c);
        const auto *md_93 = buffer.data(md + 93 * ncomps + c);
        const auto *md_94 = buffer.data(md + 94 * ncomps + c);
        const auto *md_95 = buffer.data(md + 95 * ncomps + c);
        const auto *md_96 = buffer.data(md + 96 * ncomps + c);
        const auto *md_97 = buffer.data(md + 97 * ncomps + c);
        const auto *md_98 = buffer.data(md + 98 * ncomps + c);
        const auto *md_99 = buffer.data(md + 99 * ncomps + c);
        const auto *md_100 = buffer.data(md + 100 * ncomps + c);
        const auto *md_101 = buffer.data(md + 101 * ncomps + c);
        const auto *md_102 = buffer.data(md + 102 * ncomps + c);
        const auto *md_103 = buffer.data(md + 103 * ncomps + c);
        const auto *md_104 = buffer.data(md + 104 * ncomps + c);
        const auto *md_105 = buffer.data(md + 105 * ncomps + c);
        const auto *md_106 = buffer.data(md + 106 * ncomps + c);
        const auto *md_107 = buffer.data(md + 107 * ncomps + c);
        const auto *md_108 = buffer.data(md + 108 * ncomps + c);
        const auto *md_109 = buffer.data(md + 109 * ncomps + c);
        const auto *md_110 = buffer.data(md + 110 * ncomps + c);
        const auto *md_111 = buffer.data(md + 111 * ncomps + c);
        const auto *md_112 = buffer.data(md + 112 * ncomps + c);
        const auto *md_113 = buffer.data(md + 113 * ncomps + c);
        const auto *md_114 = buffer.data(md + 114 * ncomps + c);
        const auto *md_115 = buffer.data(md + 115 * ncomps + c);
        const auto *md_116 = buffer.data(md + 116 * ncomps + c);
        const auto *md_117 = buffer.data(md + 117 * ncomps + c);
        const auto *md_118 = buffer.data(md + 118 * ncomps + c);
        const auto *md_119 = buffer.data(md + 119 * ncomps + c);
        const auto *md_120 = buffer.data(md + 120 * ncomps + c);
        const auto *md_121 = buffer.data(md + 121 * ncomps + c);
        const auto *md_122 = buffer.data(md + 122 * ncomps + c);
        const auto *md_123 = buffer.data(md + 123 * ncomps + c);
        const auto *md_124 = buffer.data(md + 124 * ncomps + c);
        const auto *md_125 = buffer.data(md + 125 * ncomps + c);
        const auto *md_126 = buffer.data(md + 126 * ncomps + c);
        const auto *md_127 = buffer.data(md + 127 * ncomps + c);
        const auto *md_128 = buffer.data(md + 128 * ncomps + c);
        const auto *md_129 = buffer.data(md + 129 * ncomps + c);
        const auto *md_130 = buffer.data(md + 130 * ncomps + c);
        const auto *md_131 = buffer.data(md + 131 * ncomps + c);
        const auto *md_132 = buffer.data(md + 132 * ncomps + c);
        const auto *md_133 = buffer.data(md + 133 * ncomps + c);
        const auto *md_134 = buffer.data(md + 134 * ncomps + c);
        const auto *md_135 = buffer.data(md + 135 * ncomps + c);
        const auto *md_136 = buffer.data(md + 136 * ncomps + c);
        const auto *md_137 = buffer.data(md + 137 * ncomps + c);
        const auto *md_138 = buffer.data(md + 138 * ncomps + c);
        const auto *md_139 = buffer.data(md + 139 * ncomps + c);
        const auto *md_140 = buffer.data(md + 140 * ncomps + c);
        const auto *md_141 = buffer.data(md + 141 * ncomps + c);
        const auto *md_142 = buffer.data(md + 142 * ncomps + c);
        const auto *md_143 = buffer.data(md + 143 * ncomps + c);
        const auto *md_144 = buffer.data(md + 144 * ncomps + c);
        const auto *md_145 = buffer.data(md + 145 * ncomps + c);
        const auto *md_146 = buffer.data(md + 146 * ncomps + c);
        const auto *md_147 = buffer.data(md + 147 * ncomps + c);
        const auto *md_148 = buffer.data(md + 148 * ncomps + c);
        const auto *md_149 = buffer.data(md + 149 * ncomps + c);
        const auto *md_150 = buffer.data(md + 150 * ncomps + c);
        const auto *md_151 = buffer.data(md + 151 * ncomps + c);
        const auto *md_152 = buffer.data(md + 152 * ncomps + c);
        const auto *md_153 = buffer.data(md + 153 * ncomps + c);
        const auto *md_154 = buffer.data(md + 154 * ncomps + c);
        const auto *md_155 = buffer.data(md + 155 * ncomps + c);
        const auto *md_156 = buffer.data(md + 156 * ncomps + c);
        const auto *md_157 = buffer.data(md + 157 * ncomps + c);
        const auto *md_158 = buffer.data(md + 158 * ncomps + c);
        const auto *md_159 = buffer.data(md + 159 * ncomps + c);
        const auto *md_160 = buffer.data(md + 160 * ncomps + c);
        const auto *md_161 = buffer.data(md + 161 * ncomps + c);
        const auto *md_162 = buffer.data(md + 162 * ncomps + c);
        const auto *md_163 = buffer.data(md + 163 * ncomps + c);
        const auto *md_164 = buffer.data(md + 164 * ncomps + c);
        const auto *md_165 = buffer.data(md + 165 * ncomps + c);
        const auto *md_166 = buffer.data(md + 166 * ncomps + c);
        const auto *md_167 = buffer.data(md + 167 * ncomps + c);
        const auto *md_168 = buffer.data(md + 168 * ncomps + c);
        const auto *md_169 = buffer.data(md + 169 * ncomps + c);
        const auto *md_170 = buffer.data(md + 170 * ncomps + c);
        const auto *md_171 = buffer.data(md + 171 * ncomps + c);
        const auto *md_172 = buffer.data(md + 172 * ncomps + c);
        const auto *md_173 = buffer.data(md + 173 * ncomps + c);
        const auto *md_177 = buffer.data(md + 177 * ncomps + c);
        const auto *md_178 = buffer.data(md + 178 * ncomps + c);
        const auto *md_179 = buffer.data(md + 179 * ncomps + c);
        const auto *md_183 = buffer.data(md + 183 * ncomps + c);
        const auto *md_184 = buffer.data(md + 184 * ncomps + c);
        const auto *md_185 = buffer.data(md + 185 * ncomps + c);
        const auto *md_189 = buffer.data(md + 189 * ncomps + c);
        const auto *md_190 = buffer.data(md + 190 * ncomps + c);
        const auto *md_191 = buffer.data(md + 191 * ncomps + c);
        const auto *md_195 = buffer.data(md + 195 * ncomps + c);
        const auto *md_196 = buffer.data(md + 196 * ncomps + c);
        const auto *md_197 = buffer.data(md + 197 * ncomps + c);
        const auto *md_201 = buffer.data(md + 201 * ncomps + c);
        const auto *md_202 = buffer.data(md + 202 * ncomps + c);
        const auto *md_203 = buffer.data(md + 203 * ncomps + c);
        const auto *md_207 = buffer.data(md + 207 * ncomps + c);
        const auto *md_208 = buffer.data(md + 208 * ncomps + c);
        const auto *md_209 = buffer.data(md + 209 * ncomps + c);
        const auto *md_215 = buffer.data(md + 215 * ncomps + c);
        const auto *md_219 = buffer.data(md + 219 * ncomps + c);
        const auto *md_220 = buffer.data(md + 220 * ncomps + c);
        const auto *md_221 = buffer.data(md + 221 * ncomps + c);
        const auto *md_227 = buffer.data(md + 227 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, ld_87, ld_88, \
                         ld_89, md_89, md_117, md_118, md_119, md_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * ld_89[k]
                       + md_89[k];

            t_146[k] = ab_y[k] * ld_87[k]
                       + md_117[k];

            t_147[k] = ab_y[k] * ld_88[k]
                       + md_118[k];

            t_148[k] = ab_y[k] * ld_89[k]
                       + md_119[k];

            t_149[k] = ab_z[k] * ld_89[k]
                       + md_125[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, ld_90, ld_91, ld_92, ld_93, \
                         ld_94, md_90, md_91, md_92, md_93, md_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ld_90[k]
                       + md_90[k];

            t_151[k] = ab_x[k] * ld_91[k]
                       + md_91[k];

            t_152[k] = ab_x[k] * ld_92[k]
                       + md_92[k];

            t_153[k] = ab_x[k] * ld_93[k]
                       + md_93[k];

            t_154[k] = ab_x[k] * ld_94[k]
                       + md_94[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, ld_93, ld_94, \
                         ld_95, md_95, md_129, md_130, md_131, md_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * ld_95[k]
                       + md_95[k];

            t_156[k] = ab_y[k] * ld_93[k]
                       + md_129[k];

            t_157[k] = ab_y[k] * ld_94[k]
                       + md_130[k];

            t_158[k] = ab_y[k] * ld_95[k]
                       + md_131[k];

            t_159[k] = ab_z[k] * ld_95[k]
                       + md_137[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, ld_96, ld_97, ld_98, ld_99, \
                         ld_100, md_96, md_97, md_98, md_99, md_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * ld_96[k]
                       + md_96[k];

            t_161[k] = ab_x[k] * ld_97[k]
                       + md_97[k];

            t_162[k] = ab_x[k] * ld_98[k]
                       + md_98[k];

            t_163[k] = ab_x[k] * ld_99[k]
                       + md_99[k];

            t_164[k] = ab_x[k] * ld_100[k]
                       + md_100[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, ld_99, ld_100, \
                         ld_101, md_101, md_135, md_136, md_137, \
                         md_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * ld_101[k]
                       + md_101[k];

            t_166[k] = ab_y[k] * ld_99[k]
                       + md_135[k];

            t_167[k] = ab_y[k] * ld_100[k]
                       + md_136[k];

            t_168[k] = ab_y[k] * ld_101[k]
                       + md_137[k];

            t_169[k] = ab_z[k] * ld_101[k]
                       + md_143[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, ld_102, ld_103, ld_104, \
                         ld_105, ld_106, md_102, md_103, md_104, md_105, \
                         md_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * ld_102[k]
                       + md_102[k];

            t_171[k] = ab_x[k] * ld_103[k]
                       + md_103[k];

            t_172[k] = ab_x[k] * ld_104[k]
                       + md_104[k];

            t_173[k] = ab_x[k] * ld_105[k]
                       + md_105[k];

            t_174[k] = ab_x[k] * ld_106[k]
                       + md_106[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, ld_105, ld_106, \
                         ld_107, md_107, md_141, md_142, md_143, \
                         md_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * ld_107[k]
                       + md_107[k];

            t_176[k] = ab_y[k] * ld_105[k]
                       + md_141[k];

            t_177[k] = ab_y[k] * ld_106[k]
                       + md_142[k];

            t_178[k] = ab_y[k] * ld_107[k]
                       + md_143[k];

            t_179[k] = ab_z[k] * ld_107[k]
                       + md_149[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, ld_108, ld_109, ld_110, \
                         ld_111, ld_112, md_108, md_109, md_110, md_111, \
                         md_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * ld_108[k]
                       + md_108[k];

            t_181[k] = ab_x[k] * ld_109[k]
                       + md_109[k];

            t_182[k] = ab_x[k] * ld_110[k]
                       + md_110[k];

            t_183[k] = ab_x[k] * ld_111[k]
                       + md_111[k];

            t_184[k] = ab_x[k] * ld_112[k]
                       + md_112[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, ld_111, ld_112, \
                         ld_113, md_113, md_147, md_148, md_149, \
                         md_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * ld_113[k]
                       + md_113[k];

            t_186[k] = ab_y[k] * ld_111[k]
                       + md_147[k];

            t_187[k] = ab_y[k] * ld_112[k]
                       + md_148[k];

            t_188[k] = ab_y[k] * ld_113[k]
                       + md_149[k];

            t_189[k] = ab_z[k] * ld_113[k]
                       + md_155[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, ld_114, ld_115, ld_116, \
                         ld_117, ld_118, md_114, md_115, md_116, md_117, \
                         md_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_x[k] * ld_114[k]
                       + md_114[k];

            t_191[k] = ab_x[k] * ld_115[k]
                       + md_115[k];

            t_192[k] = ab_x[k] * ld_116[k]
                       + md_116[k];

            t_193[k] = ab_x[k] * ld_117[k]
                       + md_117[k];

            t_194[k] = ab_x[k] * ld_118[k]
                       + md_118[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, ld_117, ld_118, \
                         ld_119, md_119, md_153, md_154, md_155, \
                         md_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * ld_119[k]
                       + md_119[k];

            t_196[k] = ab_y[k] * ld_117[k]
                       + md_153[k];

            t_197[k] = ab_y[k] * ld_118[k]
                       + md_154[k];

            t_198[k] = ab_y[k] * ld_119[k]
                       + md_155[k];

            t_199[k] = ab_z[k] * ld_119[k]
                       + md_161[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, ld_120, ld_121, ld_122, \
                         ld_123, ld_124, md_120, md_121, md_122, md_123, \
                         md_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * ld_120[k]
                       + md_120[k];

            t_201[k] = ab_x[k] * ld_121[k]
                       + md_121[k];

            t_202[k] = ab_x[k] * ld_122[k]
                       + md_122[k];

            t_203[k] = ab_x[k] * ld_123[k]
                       + md_123[k];

            t_204[k] = ab_x[k] * ld_124[k]
                       + md_124[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, ld_123, ld_124, \
                         ld_125, md_125, md_159, md_160, md_161, \
                         md_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_x[k] * ld_125[k]
                       + md_125[k];

            t_206[k] = ab_y[k] * ld_123[k]
                       + md_159[k];

            t_207[k] = ab_y[k] * ld_124[k]
                       + md_160[k];

            t_208[k] = ab_y[k] * ld_125[k]
                       + md_161[k];

            t_209[k] = ab_z[k] * ld_125[k]
                       + md_167[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, ld_126, ld_127, ld_128, \
                         ld_129, ld_130, md_126, md_127, md_128, md_129, \
                         md_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * ld_126[k]
                       + md_126[k];

            t_211[k] = ab_x[k] * ld_127[k]
                       + md_127[k];

            t_212[k] = ab_x[k] * ld_128[k]
                       + md_128[k];

            t_213[k] = ab_x[k] * ld_129[k]
                       + md_129[k];

            t_214[k] = ab_x[k] * ld_130[k]
                       + md_130[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ab_y, ab_z, ld_129, ld_130, \
                         ld_131, md_131, md_171, md_172, md_173, \
                         md_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * ld_131[k]
                       + md_131[k];

            t_216[k] = ab_y[k] * ld_129[k]
                       + md_171[k];

            t_217[k] = ab_y[k] * ld_130[k]
                       + md_172[k];

            t_218[k] = ab_y[k] * ld_131[k]
                       + md_173[k];

            t_219[k] = ab_z[k] * ld_131[k]
                       + md_179[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, ld_132, ld_133, ld_134, \
                         ld_135, ld_136, md_132, md_133, md_134, md_135, \
                         md_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * ld_132[k]
                       + md_132[k];

            t_221[k] = ab_x[k] * ld_133[k]
                       + md_133[k];

            t_222[k] = ab_x[k] * ld_134[k]
                       + md_134[k];

            t_223[k] = ab_x[k] * ld_135[k]
                       + md_135[k];

            t_224[k] = ab_x[k] * ld_136[k]
                       + md_136[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, ab_y, ab_z, ld_135, ld_136, \
                         ld_137, md_137, md_177, md_178, md_179, \
                         md_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * ld_137[k]
                       + md_137[k];

            t_226[k] = ab_y[k] * ld_135[k]
                       + md_177[k];

            t_227[k] = ab_y[k] * ld_136[k]
                       + md_178[k];

            t_228[k] = ab_y[k] * ld_137[k]
                       + md_179[k];

            t_229[k] = ab_z[k] * ld_137[k]
                       + md_185[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, ld_138, ld_139, ld_140, \
                         ld_141, ld_142, md_138, md_139, md_140, md_141, \
                         md_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * ld_138[k]
                       + md_138[k];

            t_231[k] = ab_x[k] * ld_139[k]
                       + md_139[k];

            t_232[k] = ab_x[k] * ld_140[k]
                       + md_140[k];

            t_233[k] = ab_x[k] * ld_141[k]
                       + md_141[k];

            t_234[k] = ab_x[k] * ld_142[k]
                       + md_142[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, ab_y, ab_z, ld_141, ld_142, \
                         ld_143, md_143, md_183, md_184, md_185, \
                         md_191 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_x[k] * ld_143[k]
                       + md_143[k];

            t_236[k] = ab_y[k] * ld_141[k]
                       + md_183[k];

            t_237[k] = ab_y[k] * ld_142[k]
                       + md_184[k];

            t_238[k] = ab_y[k] * ld_143[k]
                       + md_185[k];

            t_239[k] = ab_z[k] * ld_143[k]
                       + md_191[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, ld_144, ld_145, ld_146, \
                         ld_147, ld_148, md_144, md_145, md_146, md_147, \
                         md_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * ld_144[k]
                       + md_144[k];

            t_241[k] = ab_x[k] * ld_145[k]
                       + md_145[k];

            t_242[k] = ab_x[k] * ld_146[k]
                       + md_146[k];

            t_243[k] = ab_x[k] * ld_147[k]
                       + md_147[k];

            t_244[k] = ab_x[k] * ld_148[k]
                       + md_148[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, ab_y, ab_z, ld_147, ld_148, \
                         ld_149, md_149, md_189, md_190, md_191, \
                         md_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * ld_149[k]
                       + md_149[k];

            t_246[k] = ab_y[k] * ld_147[k]
                       + md_189[k];

            t_247[k] = ab_y[k] * ld_148[k]
                       + md_190[k];

            t_248[k] = ab_y[k] * ld_149[k]
                       + md_191[k];

            t_249[k] = ab_z[k] * ld_149[k]
                       + md_197[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, ld_150, ld_151, ld_152, \
                         ld_153, ld_154, md_150, md_151, md_152, md_153, \
                         md_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_x[k] * ld_150[k]
                       + md_150[k];

            t_251[k] = ab_x[k] * ld_151[k]
                       + md_151[k];

            t_252[k] = ab_x[k] * ld_152[k]
                       + md_152[k];

            t_253[k] = ab_x[k] * ld_153[k]
                       + md_153[k];

            t_254[k] = ab_x[k] * ld_154[k]
                       + md_154[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, ab_y, ab_z, ld_153, ld_154, \
                         ld_155, md_155, md_195, md_196, md_197, \
                         md_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * ld_155[k]
                       + md_155[k];

            t_256[k] = ab_y[k] * ld_153[k]
                       + md_195[k];

            t_257[k] = ab_y[k] * ld_154[k]
                       + md_196[k];

            t_258[k] = ab_y[k] * ld_155[k]
                       + md_197[k];

            t_259[k] = ab_z[k] * ld_155[k]
                       + md_203[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, ld_156, ld_157, ld_158, \
                         ld_159, ld_160, md_156, md_157, md_158, md_159, \
                         md_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * ld_156[k]
                       + md_156[k];

            t_261[k] = ab_x[k] * ld_157[k]
                       + md_157[k];

            t_262[k] = ab_x[k] * ld_158[k]
                       + md_158[k];

            t_263[k] = ab_x[k] * ld_159[k]
                       + md_159[k];

            t_264[k] = ab_x[k] * ld_160[k]
                       + md_160[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, ld_159, ld_160, \
                         ld_161, md_161, md_201, md_202, md_203, \
                         md_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_x[k] * ld_161[k]
                       + md_161[k];

            t_266[k] = ab_y[k] * ld_159[k]
                       + md_201[k];

            t_267[k] = ab_y[k] * ld_160[k]
                       + md_202[k];

            t_268[k] = ab_y[k] * ld_161[k]
                       + md_203[k];

            t_269[k] = ab_z[k] * ld_161[k]
                       + md_209[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, ld_162, ld_163, ld_164, \
                         ld_165, ld_166, md_162, md_163, md_164, md_165, \
                         md_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * ld_162[k]
                       + md_162[k];

            t_271[k] = ab_x[k] * ld_163[k]
                       + md_163[k];

            t_272[k] = ab_x[k] * ld_164[k]
                       + md_164[k];

            t_273[k] = ab_x[k] * ld_165[k]
                       + md_165[k];

            t_274[k] = ab_x[k] * ld_166[k]
                       + md_166[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, ab_y, ab_z, ld_165, ld_166, \
                         ld_167, md_167, md_207, md_208, md_209, \
                         md_215 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * ld_167[k]
                       + md_167[k];

            t_276[k] = ab_y[k] * ld_165[k]
                       + md_207[k];

            t_277[k] = ab_y[k] * ld_166[k]
                       + md_208[k];

            t_278[k] = ab_y[k] * ld_167[k]
                       + md_209[k];

            t_279[k] = ab_z[k] * ld_167[k]
                       + md_215[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, ld_168, ld_169, ld_170, \
                         ld_171, ld_172, md_168, md_169, md_170, md_171, \
                         md_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_x[k] * ld_168[k]
                       + md_168[k];

            t_281[k] = ab_x[k] * ld_169[k]
                       + md_169[k];

            t_282[k] = ab_x[k] * ld_170[k]
                       + md_170[k];

            t_283[k] = ab_x[k] * ld_171[k]
                       + md_171[k];

            t_284[k] = ab_x[k] * ld_172[k]
                       + md_172[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, ab_y, ab_z, ld_171, ld_172, \
                         ld_173, md_173, md_219, md_220, md_221, \
                         md_227 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * ld_173[k]
                       + md_173[k];

            t_286[k] = ab_y[k] * ld_171[k]
                       + md_219[k];

            t_287[k] = ab_y[k] * ld_172[k]
                       + md_220[k];

            t_288[k] = ab_y[k] * ld_173[k]
                       + md_221[k];

            t_289[k] = ab_z[k] * ld_173[k]
                       + md_227[k];
        }
    }
}

static auto
compute_hrr_lf_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ld, const size_t md,
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

        const auto *ld_174 = buffer.data(ld + 174 * ncomps + c);
        const auto *ld_175 = buffer.data(ld + 175 * ncomps + c);
        const auto *ld_176 = buffer.data(ld + 176 * ncomps + c);
        const auto *ld_177 = buffer.data(ld + 177 * ncomps + c);
        const auto *ld_178 = buffer.data(ld + 178 * ncomps + c);
        const auto *ld_179 = buffer.data(ld + 179 * ncomps + c);
        const auto *ld_180 = buffer.data(ld + 180 * ncomps + c);
        const auto *ld_181 = buffer.data(ld + 181 * ncomps + c);
        const auto *ld_182 = buffer.data(ld + 182 * ncomps + c);
        const auto *ld_183 = buffer.data(ld + 183 * ncomps + c);
        const auto *ld_184 = buffer.data(ld + 184 * ncomps + c);
        const auto *ld_185 = buffer.data(ld + 185 * ncomps + c);
        const auto *ld_186 = buffer.data(ld + 186 * ncomps + c);
        const auto *ld_187 = buffer.data(ld + 187 * ncomps + c);
        const auto *ld_188 = buffer.data(ld + 188 * ncomps + c);
        const auto *ld_189 = buffer.data(ld + 189 * ncomps + c);
        const auto *ld_190 = buffer.data(ld + 190 * ncomps + c);
        const auto *ld_191 = buffer.data(ld + 191 * ncomps + c);
        const auto *ld_192 = buffer.data(ld + 192 * ncomps + c);
        const auto *ld_193 = buffer.data(ld + 193 * ncomps + c);
        const auto *ld_194 = buffer.data(ld + 194 * ncomps + c);
        const auto *ld_195 = buffer.data(ld + 195 * ncomps + c);
        const auto *ld_196 = buffer.data(ld + 196 * ncomps + c);
        const auto *ld_197 = buffer.data(ld + 197 * ncomps + c);
        const auto *ld_198 = buffer.data(ld + 198 * ncomps + c);
        const auto *ld_199 = buffer.data(ld + 199 * ncomps + c);
        const auto *ld_200 = buffer.data(ld + 200 * ncomps + c);
        const auto *ld_201 = buffer.data(ld + 201 * ncomps + c);
        const auto *ld_202 = buffer.data(ld + 202 * ncomps + c);
        const auto *ld_203 = buffer.data(ld + 203 * ncomps + c);
        const auto *ld_204 = buffer.data(ld + 204 * ncomps + c);
        const auto *ld_205 = buffer.data(ld + 205 * ncomps + c);
        const auto *ld_206 = buffer.data(ld + 206 * ncomps + c);
        const auto *ld_207 = buffer.data(ld + 207 * ncomps + c);
        const auto *ld_208 = buffer.data(ld + 208 * ncomps + c);
        const auto *ld_209 = buffer.data(ld + 209 * ncomps + c);
        const auto *ld_210 = buffer.data(ld + 210 * ncomps + c);
        const auto *ld_211 = buffer.data(ld + 211 * ncomps + c);
        const auto *ld_212 = buffer.data(ld + 212 * ncomps + c);
        const auto *ld_213 = buffer.data(ld + 213 * ncomps + c);
        const auto *ld_214 = buffer.data(ld + 214 * ncomps + c);
        const auto *ld_215 = buffer.data(ld + 215 * ncomps + c);
        const auto *ld_216 = buffer.data(ld + 216 * ncomps + c);
        const auto *ld_217 = buffer.data(ld + 217 * ncomps + c);
        const auto *ld_218 = buffer.data(ld + 218 * ncomps + c);
        const auto *ld_219 = buffer.data(ld + 219 * ncomps + c);
        const auto *ld_220 = buffer.data(ld + 220 * ncomps + c);
        const auto *ld_221 = buffer.data(ld + 221 * ncomps + c);
        const auto *ld_222 = buffer.data(ld + 222 * ncomps + c);
        const auto *ld_223 = buffer.data(ld + 223 * ncomps + c);
        const auto *ld_224 = buffer.data(ld + 224 * ncomps + c);
        const auto *ld_225 = buffer.data(ld + 225 * ncomps + c);
        const auto *ld_226 = buffer.data(ld + 226 * ncomps + c);
        const auto *ld_227 = buffer.data(ld + 227 * ncomps + c);
        const auto *ld_228 = buffer.data(ld + 228 * ncomps + c);
        const auto *ld_229 = buffer.data(ld + 229 * ncomps + c);
        const auto *ld_230 = buffer.data(ld + 230 * ncomps + c);
        const auto *ld_231 = buffer.data(ld + 231 * ncomps + c);
        const auto *ld_232 = buffer.data(ld + 232 * ncomps + c);
        const auto *ld_233 = buffer.data(ld + 233 * ncomps + c);
        const auto *ld_234 = buffer.data(ld + 234 * ncomps + c);
        const auto *ld_235 = buffer.data(ld + 235 * ncomps + c);
        const auto *ld_236 = buffer.data(ld + 236 * ncomps + c);
        const auto *ld_237 = buffer.data(ld + 237 * ncomps + c);
        const auto *ld_238 = buffer.data(ld + 238 * ncomps + c);
        const auto *ld_239 = buffer.data(ld + 239 * ncomps + c);
        const auto *ld_240 = buffer.data(ld + 240 * ncomps + c);
        const auto *ld_241 = buffer.data(ld + 241 * ncomps + c);
        const auto *ld_242 = buffer.data(ld + 242 * ncomps + c);
        const auto *ld_243 = buffer.data(ld + 243 * ncomps + c);
        const auto *ld_244 = buffer.data(ld + 244 * ncomps + c);
        const auto *ld_245 = buffer.data(ld + 245 * ncomps + c);
        const auto *ld_246 = buffer.data(ld + 246 * ncomps + c);
        const auto *ld_247 = buffer.data(ld + 247 * ncomps + c);
        const auto *ld_248 = buffer.data(ld + 248 * ncomps + c);
        const auto *ld_249 = buffer.data(ld + 249 * ncomps + c);
        const auto *ld_250 = buffer.data(ld + 250 * ncomps + c);
        const auto *ld_251 = buffer.data(ld + 251 * ncomps + c);
        const auto *ld_252 = buffer.data(ld + 252 * ncomps + c);
        const auto *ld_253 = buffer.data(ld + 253 * ncomps + c);
        const auto *ld_254 = buffer.data(ld + 254 * ncomps + c);
        const auto *ld_255 = buffer.data(ld + 255 * ncomps + c);
        const auto *ld_256 = buffer.data(ld + 256 * ncomps + c);
        const auto *ld_257 = buffer.data(ld + 257 * ncomps + c);
        const auto *ld_258 = buffer.data(ld + 258 * ncomps + c);
        const auto *ld_259 = buffer.data(ld + 259 * ncomps + c);
        const auto *ld_260 = buffer.data(ld + 260 * ncomps + c);
        const auto *ld_261 = buffer.data(ld + 261 * ncomps + c);
        const auto *ld_262 = buffer.data(ld + 262 * ncomps + c);

        const auto *md_174 = buffer.data(md + 174 * ncomps + c);
        const auto *md_175 = buffer.data(md + 175 * ncomps + c);
        const auto *md_176 = buffer.data(md + 176 * ncomps + c);
        const auto *md_177 = buffer.data(md + 177 * ncomps + c);
        const auto *md_178 = buffer.data(md + 178 * ncomps + c);
        const auto *md_179 = buffer.data(md + 179 * ncomps + c);
        const auto *md_180 = buffer.data(md + 180 * ncomps + c);
        const auto *md_181 = buffer.data(md + 181 * ncomps + c);
        const auto *md_182 = buffer.data(md + 182 * ncomps + c);
        const auto *md_183 = buffer.data(md + 183 * ncomps + c);
        const auto *md_184 = buffer.data(md + 184 * ncomps + c);
        const auto *md_185 = buffer.data(md + 185 * ncomps + c);
        const auto *md_186 = buffer.data(md + 186 * ncomps + c);
        const auto *md_187 = buffer.data(md + 187 * ncomps + c);
        const auto *md_188 = buffer.data(md + 188 * ncomps + c);
        const auto *md_189 = buffer.data(md + 189 * ncomps + c);
        const auto *md_190 = buffer.data(md + 190 * ncomps + c);
        const auto *md_191 = buffer.data(md + 191 * ncomps + c);
        const auto *md_192 = buffer.data(md + 192 * ncomps + c);
        const auto *md_193 = buffer.data(md + 193 * ncomps + c);
        const auto *md_194 = buffer.data(md + 194 * ncomps + c);
        const auto *md_195 = buffer.data(md + 195 * ncomps + c);
        const auto *md_196 = buffer.data(md + 196 * ncomps + c);
        const auto *md_197 = buffer.data(md + 197 * ncomps + c);
        const auto *md_198 = buffer.data(md + 198 * ncomps + c);
        const auto *md_199 = buffer.data(md + 199 * ncomps + c);
        const auto *md_200 = buffer.data(md + 200 * ncomps + c);
        const auto *md_201 = buffer.data(md + 201 * ncomps + c);
        const auto *md_202 = buffer.data(md + 202 * ncomps + c);
        const auto *md_203 = buffer.data(md + 203 * ncomps + c);
        const auto *md_204 = buffer.data(md + 204 * ncomps + c);
        const auto *md_205 = buffer.data(md + 205 * ncomps + c);
        const auto *md_206 = buffer.data(md + 206 * ncomps + c);
        const auto *md_207 = buffer.data(md + 207 * ncomps + c);
        const auto *md_208 = buffer.data(md + 208 * ncomps + c);
        const auto *md_209 = buffer.data(md + 209 * ncomps + c);
        const auto *md_210 = buffer.data(md + 210 * ncomps + c);
        const auto *md_211 = buffer.data(md + 211 * ncomps + c);
        const auto *md_212 = buffer.data(md + 212 * ncomps + c);
        const auto *md_213 = buffer.data(md + 213 * ncomps + c);
        const auto *md_214 = buffer.data(md + 214 * ncomps + c);
        const auto *md_215 = buffer.data(md + 215 * ncomps + c);
        const auto *md_216 = buffer.data(md + 216 * ncomps + c);
        const auto *md_217 = buffer.data(md + 217 * ncomps + c);
        const auto *md_218 = buffer.data(md + 218 * ncomps + c);
        const auto *md_219 = buffer.data(md + 219 * ncomps + c);
        const auto *md_220 = buffer.data(md + 220 * ncomps + c);
        const auto *md_221 = buffer.data(md + 221 * ncomps + c);
        const auto *md_222 = buffer.data(md + 222 * ncomps + c);
        const auto *md_223 = buffer.data(md + 223 * ncomps + c);
        const auto *md_224 = buffer.data(md + 224 * ncomps + c);
        const auto *md_225 = buffer.data(md + 225 * ncomps + c);
        const auto *md_226 = buffer.data(md + 226 * ncomps + c);
        const auto *md_227 = buffer.data(md + 227 * ncomps + c);
        const auto *md_228 = buffer.data(md + 228 * ncomps + c);
        const auto *md_229 = buffer.data(md + 229 * ncomps + c);
        const auto *md_230 = buffer.data(md + 230 * ncomps + c);
        const auto *md_231 = buffer.data(md + 231 * ncomps + c);
        const auto *md_232 = buffer.data(md + 232 * ncomps + c);
        const auto *md_233 = buffer.data(md + 233 * ncomps + c);
        const auto *md_234 = buffer.data(md + 234 * ncomps + c);
        const auto *md_235 = buffer.data(md + 235 * ncomps + c);
        const auto *md_236 = buffer.data(md + 236 * ncomps + c);
        const auto *md_237 = buffer.data(md + 237 * ncomps + c);
        const auto *md_238 = buffer.data(md + 238 * ncomps + c);
        const auto *md_239 = buffer.data(md + 239 * ncomps + c);
        const auto *md_240 = buffer.data(md + 240 * ncomps + c);
        const auto *md_241 = buffer.data(md + 241 * ncomps + c);
        const auto *md_242 = buffer.data(md + 242 * ncomps + c);
        const auto *md_243 = buffer.data(md + 243 * ncomps + c);
        const auto *md_244 = buffer.data(md + 244 * ncomps + c);
        const auto *md_245 = buffer.data(md + 245 * ncomps + c);
        const auto *md_246 = buffer.data(md + 246 * ncomps + c);
        const auto *md_247 = buffer.data(md + 247 * ncomps + c);
        const auto *md_248 = buffer.data(md + 248 * ncomps + c);
        const auto *md_249 = buffer.data(md + 249 * ncomps + c);
        const auto *md_250 = buffer.data(md + 250 * ncomps + c);
        const auto *md_251 = buffer.data(md + 251 * ncomps + c);
        const auto *md_252 = buffer.data(md + 252 * ncomps + c);
        const auto *md_253 = buffer.data(md + 253 * ncomps + c);
        const auto *md_254 = buffer.data(md + 254 * ncomps + c);
        const auto *md_255 = buffer.data(md + 255 * ncomps + c);
        const auto *md_256 = buffer.data(md + 256 * ncomps + c);
        const auto *md_257 = buffer.data(md + 257 * ncomps + c);
        const auto *md_258 = buffer.data(md + 258 * ncomps + c);
        const auto *md_259 = buffer.data(md + 259 * ncomps + c);
        const auto *md_260 = buffer.data(md + 260 * ncomps + c);
        const auto *md_261 = buffer.data(md + 261 * ncomps + c);
        const auto *md_262 = buffer.data(md + 262 * ncomps + c);
        const auto *md_263 = buffer.data(md + 263 * ncomps + c);
        const auto *md_269 = buffer.data(md + 269 * ncomps + c);
        const auto *md_273 = buffer.data(md + 273 * ncomps + c);
        const auto *md_274 = buffer.data(md + 274 * ncomps + c);
        const auto *md_275 = buffer.data(md + 275 * ncomps + c);
        const auto *md_279 = buffer.data(md + 279 * ncomps + c);
        const auto *md_280 = buffer.data(md + 280 * ncomps + c);
        const auto *md_281 = buffer.data(md + 281 * ncomps + c);
        const auto *md_285 = buffer.data(md + 285 * ncomps + c);
        const auto *md_286 = buffer.data(md + 286 * ncomps + c);
        const auto *md_287 = buffer.data(md + 287 * ncomps + c);
        const auto *md_291 = buffer.data(md + 291 * ncomps + c);
        const auto *md_292 = buffer.data(md + 292 * ncomps + c);
        const auto *md_293 = buffer.data(md + 293 * ncomps + c);
        const auto *md_297 = buffer.data(md + 297 * ncomps + c);
        const auto *md_298 = buffer.data(md + 298 * ncomps + c);
        const auto *md_299 = buffer.data(md + 299 * ncomps + c);
        const auto *md_303 = buffer.data(md + 303 * ncomps + c);
        const auto *md_304 = buffer.data(md + 304 * ncomps + c);
        const auto *md_305 = buffer.data(md + 305 * ncomps + c);
        const auto *md_309 = buffer.data(md + 309 * ncomps + c);
        const auto *md_310 = buffer.data(md + 310 * ncomps + c);
        const auto *md_311 = buffer.data(md + 311 * ncomps + c);
        const auto *md_317 = buffer.data(md + 317 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, ld_174, ld_175, ld_176, \
                         ld_177, ld_178, md_174, md_175, md_176, md_177, \
                         md_178 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * ld_174[k]
                       + md_174[k];

            t_291[k] = ab_x[k] * ld_175[k]
                       + md_175[k];

            t_292[k] = ab_x[k] * ld_176[k]
                       + md_176[k];

            t_293[k] = ab_x[k] * ld_177[k]
                       + md_177[k];

            t_294[k] = ab_x[k] * ld_178[k]
                       + md_178[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, ab_y, ab_z, ld_177, ld_178, \
                         ld_179, md_179, md_225, md_226, md_227, \
                         md_233 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_x[k] * ld_179[k]
                       + md_179[k];

            t_296[k] = ab_y[k] * ld_177[k]
                       + md_225[k];

            t_297[k] = ab_y[k] * ld_178[k]
                       + md_226[k];

            t_298[k] = ab_y[k] * ld_179[k]
                       + md_227[k];

            t_299[k] = ab_z[k] * ld_179[k]
                       + md_233[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, ld_180, ld_181, ld_182, \
                         ld_183, ld_184, md_180, md_181, md_182, md_183, \
                         md_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * ld_180[k]
                       + md_180[k];

            t_301[k] = ab_x[k] * ld_181[k]
                       + md_181[k];

            t_302[k] = ab_x[k] * ld_182[k]
                       + md_182[k];

            t_303[k] = ab_x[k] * ld_183[k]
                       + md_183[k];

            t_304[k] = ab_x[k] * ld_184[k]
                       + md_184[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, ab_y, ab_z, ld_183, ld_184, \
                         ld_185, md_185, md_231, md_232, md_233, \
                         md_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * ld_185[k]
                       + md_185[k];

            t_306[k] = ab_y[k] * ld_183[k]
                       + md_231[k];

            t_307[k] = ab_y[k] * ld_184[k]
                       + md_232[k];

            t_308[k] = ab_y[k] * ld_185[k]
                       + md_233[k];

            t_309[k] = ab_z[k] * ld_185[k]
                       + md_239[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, ld_186, ld_187, ld_188, \
                         ld_189, ld_190, md_186, md_187, md_188, md_189, \
                         md_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_x[k] * ld_186[k]
                       + md_186[k];

            t_311[k] = ab_x[k] * ld_187[k]
                       + md_187[k];

            t_312[k] = ab_x[k] * ld_188[k]
                       + md_188[k];

            t_313[k] = ab_x[k] * ld_189[k]
                       + md_189[k];

            t_314[k] = ab_x[k] * ld_190[k]
                       + md_190[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, ab_y, ab_z, ld_189, ld_190, \
                         ld_191, md_191, md_237, md_238, md_239, \
                         md_245 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * ld_191[k]
                       + md_191[k];

            t_316[k] = ab_y[k] * ld_189[k]
                       + md_237[k];

            t_317[k] = ab_y[k] * ld_190[k]
                       + md_238[k];

            t_318[k] = ab_y[k] * ld_191[k]
                       + md_239[k];

            t_319[k] = ab_z[k] * ld_191[k]
                       + md_245[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, ld_192, ld_193, ld_194, \
                         ld_195, ld_196, md_192, md_193, md_194, md_195, \
                         md_196 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * ld_192[k]
                       + md_192[k];

            t_321[k] = ab_x[k] * ld_193[k]
                       + md_193[k];

            t_322[k] = ab_x[k] * ld_194[k]
                       + md_194[k];

            t_323[k] = ab_x[k] * ld_195[k]
                       + md_195[k];

            t_324[k] = ab_x[k] * ld_196[k]
                       + md_196[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, ab_y, ab_z, ld_195, ld_196, \
                         ld_197, md_197, md_243, md_244, md_245, \
                         md_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_x[k] * ld_197[k]
                       + md_197[k];

            t_326[k] = ab_y[k] * ld_195[k]
                       + md_243[k];

            t_327[k] = ab_y[k] * ld_196[k]
                       + md_244[k];

            t_328[k] = ab_y[k] * ld_197[k]
                       + md_245[k];

            t_329[k] = ab_z[k] * ld_197[k]
                       + md_251[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, ld_198, ld_199, ld_200, \
                         ld_201, ld_202, md_198, md_199, md_200, md_201, \
                         md_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * ld_198[k]
                       + md_198[k];

            t_331[k] = ab_x[k] * ld_199[k]
                       + md_199[k];

            t_332[k] = ab_x[k] * ld_200[k]
                       + md_200[k];

            t_333[k] = ab_x[k] * ld_201[k]
                       + md_201[k];

            t_334[k] = ab_x[k] * ld_202[k]
                       + md_202[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, ab_y, ab_z, ld_201, ld_202, \
                         ld_203, md_203, md_249, md_250, md_251, \
                         md_257 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * ld_203[k]
                       + md_203[k];

            t_336[k] = ab_y[k] * ld_201[k]
                       + md_249[k];

            t_337[k] = ab_y[k] * ld_202[k]
                       + md_250[k];

            t_338[k] = ab_y[k] * ld_203[k]
                       + md_251[k];

            t_339[k] = ab_z[k] * ld_203[k]
                       + md_257[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, ld_204, ld_205, ld_206, \
                         ld_207, ld_208, md_204, md_205, md_206, md_207, \
                         md_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_x[k] * ld_204[k]
                       + md_204[k];

            t_341[k] = ab_x[k] * ld_205[k]
                       + md_205[k];

            t_342[k] = ab_x[k] * ld_206[k]
                       + md_206[k];

            t_343[k] = ab_x[k] * ld_207[k]
                       + md_207[k];

            t_344[k] = ab_x[k] * ld_208[k]
                       + md_208[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, ab_y, ab_z, ld_207, ld_208, \
                         ld_209, md_209, md_255, md_256, md_257, \
                         md_263 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * ld_209[k]
                       + md_209[k];

            t_346[k] = ab_y[k] * ld_207[k]
                       + md_255[k];

            t_347[k] = ab_y[k] * ld_208[k]
                       + md_256[k];

            t_348[k] = ab_y[k] * ld_209[k]
                       + md_257[k];

            t_349[k] = ab_z[k] * ld_209[k]
                       + md_263[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, ld_210, ld_211, ld_212, \
                         ld_213, ld_214, md_210, md_211, md_212, md_213, \
                         md_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * ld_210[k]
                       + md_210[k];

            t_351[k] = ab_x[k] * ld_211[k]
                       + md_211[k];

            t_352[k] = ab_x[k] * ld_212[k]
                       + md_212[k];

            t_353[k] = ab_x[k] * ld_213[k]
                       + md_213[k];

            t_354[k] = ab_x[k] * ld_214[k]
                       + md_214[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, ab_y, ab_z, ld_213, ld_214, \
                         ld_215, md_215, md_261, md_262, md_263, \
                         md_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_x[k] * ld_215[k]
                       + md_215[k];

            t_356[k] = ab_y[k] * ld_213[k]
                       + md_261[k];

            t_357[k] = ab_y[k] * ld_214[k]
                       + md_262[k];

            t_358[k] = ab_y[k] * ld_215[k]
                       + md_263[k];

            t_359[k] = ab_z[k] * ld_215[k]
                       + md_269[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, ld_216, ld_217, ld_218, \
                         ld_219, ld_220, md_216, md_217, md_218, md_219, \
                         md_220 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * ld_216[k]
                       + md_216[k];

            t_361[k] = ab_x[k] * ld_217[k]
                       + md_217[k];

            t_362[k] = ab_x[k] * ld_218[k]
                       + md_218[k];

            t_363[k] = ab_x[k] * ld_219[k]
                       + md_219[k];

            t_364[k] = ab_x[k] * ld_220[k]
                       + md_220[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, ab_y, ab_z, ld_219, ld_220, \
                         ld_221, md_221, md_273, md_274, md_275, \
                         md_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * ld_221[k]
                       + md_221[k];

            t_366[k] = ab_y[k] * ld_219[k]
                       + md_273[k];

            t_367[k] = ab_y[k] * ld_220[k]
                       + md_274[k];

            t_368[k] = ab_y[k] * ld_221[k]
                       + md_275[k];

            t_369[k] = ab_z[k] * ld_221[k]
                       + md_281[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, ld_222, ld_223, ld_224, \
                         ld_225, ld_226, md_222, md_223, md_224, md_225, \
                         md_226 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * ld_222[k]
                       + md_222[k];

            t_371[k] = ab_x[k] * ld_223[k]
                       + md_223[k];

            t_372[k] = ab_x[k] * ld_224[k]
                       + md_224[k];

            t_373[k] = ab_x[k] * ld_225[k]
                       + md_225[k];

            t_374[k] = ab_x[k] * ld_226[k]
                       + md_226[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, ab_y, ab_z, ld_225, ld_226, \
                         ld_227, md_227, md_279, md_280, md_281, \
                         md_287 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * ld_227[k]
                       + md_227[k];

            t_376[k] = ab_y[k] * ld_225[k]
                       + md_279[k];

            t_377[k] = ab_y[k] * ld_226[k]
                       + md_280[k];

            t_378[k] = ab_y[k] * ld_227[k]
                       + md_281[k];

            t_379[k] = ab_z[k] * ld_227[k]
                       + md_287[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, ld_228, ld_229, ld_230, \
                         ld_231, ld_232, md_228, md_229, md_230, md_231, \
                         md_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * ld_228[k]
                       + md_228[k];

            t_381[k] = ab_x[k] * ld_229[k]
                       + md_229[k];

            t_382[k] = ab_x[k] * ld_230[k]
                       + md_230[k];

            t_383[k] = ab_x[k] * ld_231[k]
                       + md_231[k];

            t_384[k] = ab_x[k] * ld_232[k]
                       + md_232[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, ab_y, ab_z, ld_231, ld_232, \
                         ld_233, md_233, md_285, md_286, md_287, \
                         md_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_x[k] * ld_233[k]
                       + md_233[k];

            t_386[k] = ab_y[k] * ld_231[k]
                       + md_285[k];

            t_387[k] = ab_y[k] * ld_232[k]
                       + md_286[k];

            t_388[k] = ab_y[k] * ld_233[k]
                       + md_287[k];

            t_389[k] = ab_z[k] * ld_233[k]
                       + md_293[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, ld_234, ld_235, ld_236, \
                         ld_237, ld_238, md_234, md_235, md_236, md_237, \
                         md_238 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * ld_234[k]
                       + md_234[k];

            t_391[k] = ab_x[k] * ld_235[k]
                       + md_235[k];

            t_392[k] = ab_x[k] * ld_236[k]
                       + md_236[k];

            t_393[k] = ab_x[k] * ld_237[k]
                       + md_237[k];

            t_394[k] = ab_x[k] * ld_238[k]
                       + md_238[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, ab_y, ab_z, ld_237, ld_238, \
                         ld_239, md_239, md_291, md_292, md_293, \
                         md_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * ld_239[k]
                       + md_239[k];

            t_396[k] = ab_y[k] * ld_237[k]
                       + md_291[k];

            t_397[k] = ab_y[k] * ld_238[k]
                       + md_292[k];

            t_398[k] = ab_y[k] * ld_239[k]
                       + md_293[k];

            t_399[k] = ab_z[k] * ld_239[k]
                       + md_299[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, ld_240, ld_241, ld_242, \
                         ld_243, ld_244, md_240, md_241, md_242, md_243, \
                         md_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_x[k] * ld_240[k]
                       + md_240[k];

            t_401[k] = ab_x[k] * ld_241[k]
                       + md_241[k];

            t_402[k] = ab_x[k] * ld_242[k]
                       + md_242[k];

            t_403[k] = ab_x[k] * ld_243[k]
                       + md_243[k];

            t_404[k] = ab_x[k] * ld_244[k]
                       + md_244[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, ab_y, ab_z, ld_243, ld_244, \
                         ld_245, md_245, md_297, md_298, md_299, \
                         md_305 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * ld_245[k]
                       + md_245[k];

            t_406[k] = ab_y[k] * ld_243[k]
                       + md_297[k];

            t_407[k] = ab_y[k] * ld_244[k]
                       + md_298[k];

            t_408[k] = ab_y[k] * ld_245[k]
                       + md_299[k];

            t_409[k] = ab_z[k] * ld_245[k]
                       + md_305[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, ld_246, ld_247, ld_248, \
                         ld_249, ld_250, md_246, md_247, md_248, md_249, \
                         md_250 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * ld_246[k]
                       + md_246[k];

            t_411[k] = ab_x[k] * ld_247[k]
                       + md_247[k];

            t_412[k] = ab_x[k] * ld_248[k]
                       + md_248[k];

            t_413[k] = ab_x[k] * ld_249[k]
                       + md_249[k];

            t_414[k] = ab_x[k] * ld_250[k]
                       + md_250[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, ab_y, ab_z, ld_249, ld_250, \
                         ld_251, md_251, md_303, md_304, md_305, \
                         md_311 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_x[k] * ld_251[k]
                       + md_251[k];

            t_416[k] = ab_y[k] * ld_249[k]
                       + md_303[k];

            t_417[k] = ab_y[k] * ld_250[k]
                       + md_304[k];

            t_418[k] = ab_y[k] * ld_251[k]
                       + md_305[k];

            t_419[k] = ab_z[k] * ld_251[k]
                       + md_311[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, ld_252, ld_253, ld_254, \
                         ld_255, ld_256, md_252, md_253, md_254, md_255, \
                         md_256 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * ld_252[k]
                       + md_252[k];

            t_421[k] = ab_x[k] * ld_253[k]
                       + md_253[k];

            t_422[k] = ab_x[k] * ld_254[k]
                       + md_254[k];

            t_423[k] = ab_x[k] * ld_255[k]
                       + md_255[k];

            t_424[k] = ab_x[k] * ld_256[k]
                       + md_256[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, ab_y, ab_z, ld_255, ld_256, \
                         ld_257, md_257, md_309, md_310, md_311, \
                         md_317 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * ld_257[k]
                       + md_257[k];

            t_426[k] = ab_y[k] * ld_255[k]
                       + md_309[k];

            t_427[k] = ab_y[k] * ld_256[k]
                       + md_310[k];

            t_428[k] = ab_y[k] * ld_257[k]
                       + md_311[k];

            t_429[k] = ab_z[k] * ld_257[k]
                       + md_317[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, ld_258, ld_259, ld_260, \
                         ld_261, ld_262, md_258, md_259, md_260, md_261, \
                         md_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * ld_258[k]
                       + md_258[k];

            t_431[k] = ab_x[k] * ld_259[k]
                       + md_259[k];

            t_432[k] = ab_x[k] * ld_260[k]
                       + md_260[k];

            t_433[k] = ab_x[k] * ld_261[k]
                       + md_261[k];

            t_434[k] = ab_x[k] * ld_262[k]
                       + md_262[k];
        }
    }
}

static auto
compute_hrr_lf_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ld, const size_t md,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ld_261 = buffer.data(ld + 261 * ncomps + c);
        const auto *ld_262 = buffer.data(ld + 262 * ncomps + c);
        const auto *ld_263 = buffer.data(ld + 263 * ncomps + c);
        const auto *ld_264 = buffer.data(ld + 264 * ncomps + c);
        const auto *ld_265 = buffer.data(ld + 265 * ncomps + c);
        const auto *ld_266 = buffer.data(ld + 266 * ncomps + c);
        const auto *ld_267 = buffer.data(ld + 267 * ncomps + c);
        const auto *ld_268 = buffer.data(ld + 268 * ncomps + c);
        const auto *ld_269 = buffer.data(ld + 269 * ncomps + c);

        const auto *md_263 = buffer.data(md + 263 * ncomps + c);
        const auto *md_264 = buffer.data(md + 264 * ncomps + c);
        const auto *md_265 = buffer.data(md + 265 * ncomps + c);
        const auto *md_266 = buffer.data(md + 266 * ncomps + c);
        const auto *md_267 = buffer.data(md + 267 * ncomps + c);
        const auto *md_268 = buffer.data(md + 268 * ncomps + c);
        const auto *md_269 = buffer.data(md + 269 * ncomps + c);
        const auto *md_315 = buffer.data(md + 315 * ncomps + c);
        const auto *md_316 = buffer.data(md + 316 * ncomps + c);
        const auto *md_317 = buffer.data(md + 317 * ncomps + c);
        const auto *md_321 = buffer.data(md + 321 * ncomps + c);
        const auto *md_322 = buffer.data(md + 322 * ncomps + c);
        const auto *md_323 = buffer.data(md + 323 * ncomps + c);
        const auto *md_329 = buffer.data(md + 329 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, ab_y, ab_z, ld_261, ld_262, \
                         ld_263, md_263, md_315, md_316, md_317, \
                         md_323 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * ld_263[k]
                       + md_263[k];

            t_436[k] = ab_y[k] * ld_261[k]
                       + md_315[k];

            t_437[k] = ab_y[k] * ld_262[k]
                       + md_316[k];

            t_438[k] = ab_y[k] * ld_263[k]
                       + md_317[k];

            t_439[k] = ab_z[k] * ld_263[k]
                       + md_323[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, ld_264, ld_265, ld_266, \
                         ld_267, ld_268, md_264, md_265, md_266, md_267, \
                         md_268 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * ld_264[k]
                       + md_264[k];

            t_441[k] = ab_x[k] * ld_265[k]
                       + md_265[k];

            t_442[k] = ab_x[k] * ld_266[k]
                       + md_266[k];

            t_443[k] = ab_x[k] * ld_267[k]
                       + md_267[k];

            t_444[k] = ab_x[k] * ld_268[k]
                       + md_268[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_x, ab_y, ab_z, ld_267, ld_268, \
                         ld_269, md_269, md_321, md_322, md_323, \
                         md_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = ab_x[k] * ld_269[k]
                       + md_269[k];

            t_446[k] = ab_y[k] * ld_267[k]
                       + md_321[k];

            t_447[k] = ab_y[k] * ld_268[k]
                       + md_322[k];

            t_448[k] = ab_y[k] * ld_269[k]
                       + md_323[k];

            t_449[k] = ab_z[k] * ld_269[k]
                       + md_329[k];
        }
    }
}

auto
compute_hrr_lf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t ld, const size_t md,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_lf_out_of_first_piece0(buffer, coordinates, target, ld, md, ncomps, nmax);

    compute_hrr_lf_out_of_first_piece1(buffer, coordinates, target, ld, md, ncomps, nmax);

    compute_hrr_lf_out_of_first_piece2(buffer, coordinates, target, ld, md, ncomps, nmax);

    compute_hrr_lf_out_of_first_piece3(buffer, coordinates, target, ld, md, ncomps, nmax);
}

static auto
compute_hrr_lf_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ld, const size_t md, const size_t ncomps,
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

        const auto *ld_0 = buffer.data(ld + 0 * ncomps + c);
        const auto *ld_1 = buffer.data(ld + 1 * ncomps + c);
        const auto *ld_2 = buffer.data(ld + 2 * ncomps + c);
        const auto *ld_3 = buffer.data(ld + 3 * ncomps + c);
        const auto *ld_4 = buffer.data(ld + 4 * ncomps + c);
        const auto *ld_5 = buffer.data(ld + 5 * ncomps + c);
        const auto *ld_6 = buffer.data(ld + 6 * ncomps + c);
        const auto *ld_7 = buffer.data(ld + 7 * ncomps + c);
        const auto *ld_8 = buffer.data(ld + 8 * ncomps + c);
        const auto *ld_9 = buffer.data(ld + 9 * ncomps + c);
        const auto *ld_10 = buffer.data(ld + 10 * ncomps + c);
        const auto *ld_11 = buffer.data(ld + 11 * ncomps + c);
        const auto *ld_12 = buffer.data(ld + 12 * ncomps + c);
        const auto *ld_13 = buffer.data(ld + 13 * ncomps + c);
        const auto *ld_14 = buffer.data(ld + 14 * ncomps + c);
        const auto *ld_15 = buffer.data(ld + 15 * ncomps + c);
        const auto *ld_16 = buffer.data(ld + 16 * ncomps + c);
        const auto *ld_17 = buffer.data(ld + 17 * ncomps + c);
        const auto *ld_18 = buffer.data(ld + 18 * ncomps + c);
        const auto *ld_19 = buffer.data(ld + 19 * ncomps + c);
        const auto *ld_20 = buffer.data(ld + 20 * ncomps + c);
        const auto *ld_21 = buffer.data(ld + 21 * ncomps + c);
        const auto *ld_22 = buffer.data(ld + 22 * ncomps + c);
        const auto *ld_23 = buffer.data(ld + 23 * ncomps + c);
        const auto *ld_24 = buffer.data(ld + 24 * ncomps + c);
        const auto *ld_25 = buffer.data(ld + 25 * ncomps + c);
        const auto *ld_26 = buffer.data(ld + 26 * ncomps + c);
        const auto *ld_27 = buffer.data(ld + 27 * ncomps + c);
        const auto *ld_28 = buffer.data(ld + 28 * ncomps + c);
        const auto *ld_29 = buffer.data(ld + 29 * ncomps + c);
        const auto *ld_30 = buffer.data(ld + 30 * ncomps + c);
        const auto *ld_31 = buffer.data(ld + 31 * ncomps + c);
        const auto *ld_32 = buffer.data(ld + 32 * ncomps + c);
        const auto *ld_33 = buffer.data(ld + 33 * ncomps + c);
        const auto *ld_34 = buffer.data(ld + 34 * ncomps + c);
        const auto *ld_35 = buffer.data(ld + 35 * ncomps + c);
        const auto *ld_36 = buffer.data(ld + 36 * ncomps + c);
        const auto *ld_37 = buffer.data(ld + 37 * ncomps + c);
        const auto *ld_38 = buffer.data(ld + 38 * ncomps + c);
        const auto *ld_39 = buffer.data(ld + 39 * ncomps + c);
        const auto *ld_40 = buffer.data(ld + 40 * ncomps + c);
        const auto *ld_41 = buffer.data(ld + 41 * ncomps + c);
        const auto *ld_42 = buffer.data(ld + 42 * ncomps + c);
        const auto *ld_43 = buffer.data(ld + 43 * ncomps + c);
        const auto *ld_44 = buffer.data(ld + 44 * ncomps + c);
        const auto *ld_45 = buffer.data(ld + 45 * ncomps + c);
        const auto *ld_46 = buffer.data(ld + 46 * ncomps + c);
        const auto *ld_47 = buffer.data(ld + 47 * ncomps + c);
        const auto *ld_48 = buffer.data(ld + 48 * ncomps + c);
        const auto *ld_49 = buffer.data(ld + 49 * ncomps + c);
        const auto *ld_50 = buffer.data(ld + 50 * ncomps + c);
        const auto *ld_51 = buffer.data(ld + 51 * ncomps + c);
        const auto *ld_52 = buffer.data(ld + 52 * ncomps + c);
        const auto *ld_53 = buffer.data(ld + 53 * ncomps + c);
        const auto *ld_54 = buffer.data(ld + 54 * ncomps + c);
        const auto *ld_55 = buffer.data(ld + 55 * ncomps + c);
        const auto *ld_56 = buffer.data(ld + 56 * ncomps + c);
        const auto *ld_57 = buffer.data(ld + 57 * ncomps + c);
        const auto *ld_58 = buffer.data(ld + 58 * ncomps + c);
        const auto *ld_59 = buffer.data(ld + 59 * ncomps + c);
        const auto *ld_60 = buffer.data(ld + 60 * ncomps + c);
        const auto *ld_61 = buffer.data(ld + 61 * ncomps + c);
        const auto *ld_62 = buffer.data(ld + 62 * ncomps + c);
        const auto *ld_63 = buffer.data(ld + 63 * ncomps + c);
        const auto *ld_64 = buffer.data(ld + 64 * ncomps + c);
        const auto *ld_65 = buffer.data(ld + 65 * ncomps + c);
        const auto *ld_66 = buffer.data(ld + 66 * ncomps + c);
        const auto *ld_67 = buffer.data(ld + 67 * ncomps + c);
        const auto *ld_68 = buffer.data(ld + 68 * ncomps + c);
        const auto *ld_69 = buffer.data(ld + 69 * ncomps + c);
        const auto *ld_70 = buffer.data(ld + 70 * ncomps + c);
        const auto *ld_71 = buffer.data(ld + 71 * ncomps + c);
        const auto *ld_72 = buffer.data(ld + 72 * ncomps + c);
        const auto *ld_73 = buffer.data(ld + 73 * ncomps + c);
        const auto *ld_74 = buffer.data(ld + 74 * ncomps + c);
        const auto *ld_75 = buffer.data(ld + 75 * ncomps + c);
        const auto *ld_76 = buffer.data(ld + 76 * ncomps + c);
        const auto *ld_77 = buffer.data(ld + 77 * ncomps + c);
        const auto *ld_78 = buffer.data(ld + 78 * ncomps + c);
        const auto *ld_79 = buffer.data(ld + 79 * ncomps + c);
        const auto *ld_80 = buffer.data(ld + 80 * ncomps + c);
        const auto *ld_81 = buffer.data(ld + 81 * ncomps + c);
        const auto *ld_82 = buffer.data(ld + 82 * ncomps + c);
        const auto *ld_83 = buffer.data(ld + 83 * ncomps + c);
        const auto *ld_84 = buffer.data(ld + 84 * ncomps + c);
        const auto *ld_85 = buffer.data(ld + 85 * ncomps + c);
        const auto *ld_86 = buffer.data(ld + 86 * ncomps + c);
        const auto *ld_87 = buffer.data(ld + 87 * ncomps + c);
        const auto *ld_88 = buffer.data(ld + 88 * ncomps + c);

        const auto *md_0 = buffer.data(md + 0 * ncomps + c);
        const auto *md_1 = buffer.data(md + 1 * ncomps + c);
        const auto *md_2 = buffer.data(md + 2 * ncomps + c);
        const auto *md_3 = buffer.data(md + 3 * ncomps + c);
        const auto *md_4 = buffer.data(md + 4 * ncomps + c);
        const auto *md_5 = buffer.data(md + 5 * ncomps + c);
        const auto *md_6 = buffer.data(md + 6 * ncomps + c);
        const auto *md_7 = buffer.data(md + 7 * ncomps + c);
        const auto *md_8 = buffer.data(md + 8 * ncomps + c);
        const auto *md_9 = buffer.data(md + 9 * ncomps + c);
        const auto *md_10 = buffer.data(md + 10 * ncomps + c);
        const auto *md_11 = buffer.data(md + 11 * ncomps + c);
        const auto *md_12 = buffer.data(md + 12 * ncomps + c);
        const auto *md_13 = buffer.data(md + 13 * ncomps + c);
        const auto *md_14 = buffer.data(md + 14 * ncomps + c);
        const auto *md_15 = buffer.data(md + 15 * ncomps + c);
        const auto *md_16 = buffer.data(md + 16 * ncomps + c);
        const auto *md_17 = buffer.data(md + 17 * ncomps + c);
        const auto *md_18 = buffer.data(md + 18 * ncomps + c);
        const auto *md_19 = buffer.data(md + 19 * ncomps + c);
        const auto *md_20 = buffer.data(md + 20 * ncomps + c);
        const auto *md_21 = buffer.data(md + 21 * ncomps + c);
        const auto *md_22 = buffer.data(md + 22 * ncomps + c);
        const auto *md_23 = buffer.data(md + 23 * ncomps + c);
        const auto *md_24 = buffer.data(md + 24 * ncomps + c);
        const auto *md_25 = buffer.data(md + 25 * ncomps + c);
        const auto *md_26 = buffer.data(md + 26 * ncomps + c);
        const auto *md_27 = buffer.data(md + 27 * ncomps + c);
        const auto *md_28 = buffer.data(md + 28 * ncomps + c);
        const auto *md_29 = buffer.data(md + 29 * ncomps + c);
        const auto *md_30 = buffer.data(md + 30 * ncomps + c);
        const auto *md_31 = buffer.data(md + 31 * ncomps + c);
        const auto *md_32 = buffer.data(md + 32 * ncomps + c);
        const auto *md_33 = buffer.data(md + 33 * ncomps + c);
        const auto *md_34 = buffer.data(md + 34 * ncomps + c);
        const auto *md_35 = buffer.data(md + 35 * ncomps + c);
        const auto *md_36 = buffer.data(md + 36 * ncomps + c);
        const auto *md_37 = buffer.data(md + 37 * ncomps + c);
        const auto *md_38 = buffer.data(md + 38 * ncomps + c);
        const auto *md_39 = buffer.data(md + 39 * ncomps + c);
        const auto *md_40 = buffer.data(md + 40 * ncomps + c);
        const auto *md_41 = buffer.data(md + 41 * ncomps + c);
        const auto *md_42 = buffer.data(md + 42 * ncomps + c);
        const auto *md_43 = buffer.data(md + 43 * ncomps + c);
        const auto *md_44 = buffer.data(md + 44 * ncomps + c);
        const auto *md_45 = buffer.data(md + 45 * ncomps + c);
        const auto *md_46 = buffer.data(md + 46 * ncomps + c);
        const auto *md_47 = buffer.data(md + 47 * ncomps + c);
        const auto *md_48 = buffer.data(md + 48 * ncomps + c);
        const auto *md_49 = buffer.data(md + 49 * ncomps + c);
        const auto *md_50 = buffer.data(md + 50 * ncomps + c);
        const auto *md_51 = buffer.data(md + 51 * ncomps + c);
        const auto *md_52 = buffer.data(md + 52 * ncomps + c);
        const auto *md_53 = buffer.data(md + 53 * ncomps + c);
        const auto *md_54 = buffer.data(md + 54 * ncomps + c);
        const auto *md_55 = buffer.data(md + 55 * ncomps + c);
        const auto *md_56 = buffer.data(md + 56 * ncomps + c);
        const auto *md_57 = buffer.data(md + 57 * ncomps + c);
        const auto *md_58 = buffer.data(md + 58 * ncomps + c);
        const auto *md_59 = buffer.data(md + 59 * ncomps + c);
        const auto *md_60 = buffer.data(md + 60 * ncomps + c);
        const auto *md_61 = buffer.data(md + 61 * ncomps + c);
        const auto *md_62 = buffer.data(md + 62 * ncomps + c);
        const auto *md_63 = buffer.data(md + 63 * ncomps + c);
        const auto *md_64 = buffer.data(md + 64 * ncomps + c);
        const auto *md_65 = buffer.data(md + 65 * ncomps + c);
        const auto *md_66 = buffer.data(md + 66 * ncomps + c);
        const auto *md_67 = buffer.data(md + 67 * ncomps + c);
        const auto *md_68 = buffer.data(md + 68 * ncomps + c);
        const auto *md_69 = buffer.data(md + 69 * ncomps + c);
        const auto *md_70 = buffer.data(md + 70 * ncomps + c);
        const auto *md_71 = buffer.data(md + 71 * ncomps + c);
        const auto *md_72 = buffer.data(md + 72 * ncomps + c);
        const auto *md_73 = buffer.data(md + 73 * ncomps + c);
        const auto *md_74 = buffer.data(md + 74 * ncomps + c);
        const auto *md_75 = buffer.data(md + 75 * ncomps + c);
        const auto *md_76 = buffer.data(md + 76 * ncomps + c);
        const auto *md_77 = buffer.data(md + 77 * ncomps + c);
        const auto *md_78 = buffer.data(md + 78 * ncomps + c);
        const auto *md_79 = buffer.data(md + 79 * ncomps + c);
        const auto *md_80 = buffer.data(md + 80 * ncomps + c);
        const auto *md_81 = buffer.data(md + 81 * ncomps + c);
        const auto *md_82 = buffer.data(md + 82 * ncomps + c);
        const auto *md_83 = buffer.data(md + 83 * ncomps + c);
        const auto *md_84 = buffer.data(md + 84 * ncomps + c);
        const auto *md_85 = buffer.data(md + 85 * ncomps + c);
        const auto *md_86 = buffer.data(md + 86 * ncomps + c);
        const auto *md_87 = buffer.data(md + 87 * ncomps + c);
        const auto *md_88 = buffer.data(md + 88 * ncomps + c);
        const auto *md_89 = buffer.data(md + 89 * ncomps + c);
        const auto *md_93 = buffer.data(md + 93 * ncomps + c);
        const auto *md_94 = buffer.data(md + 94 * ncomps + c);
        const auto *md_95 = buffer.data(md + 95 * ncomps + c);
        const auto *md_99 = buffer.data(md + 99 * ncomps + c);
        const auto *md_100 = buffer.data(md + 100 * ncomps + c);
        const auto *md_101 = buffer.data(md + 101 * ncomps + c);
        const auto *md_105 = buffer.data(md + 105 * ncomps + c);
        const auto *md_106 = buffer.data(md + 106 * ncomps + c);
        const auto *md_107 = buffer.data(md + 107 * ncomps + c);
        const auto *md_111 = buffer.data(md + 111 * ncomps + c);
        const auto *md_112 = buffer.data(md + 112 * ncomps + c);
        const auto *md_113 = buffer.data(md + 113 * ncomps + c);
        const auto *md_119 = buffer.data(md + 119 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ld_0, ld_1, ld_2, ld_3, ld_4, md_0, \
                         md_1, md_2, md_3, md_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ld_0[k]
                     + md_0[k];

            t_1[k] = ab_x[k] * ld_1[k]
                     + md_1[k];

            t_2[k] = ab_x[k] * ld_2[k]
                     + md_2[k];

            t_3[k] = ab_x[k] * ld_3[k]
                     + md_3[k];

            t_4[k] = ab_x[k] * ld_4[k]
                     + md_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, ld_3, ld_4, ld_5, md_5, \
                         md_9, md_10, md_11, md_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * ld_5[k]
                     + md_5[k];

            t_6[k] = ab_y[k] * ld_3[k]
                     + md_9[k];

            t_7[k] = ab_y[k] * ld_4[k]
                     + md_10[k];

            t_8[k] = ab_y[k] * ld_5[k]
                     + md_11[k];

            t_9[k] = ab_z[k] * ld_5[k]
                     + md_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ld_6, ld_7, ld_8, ld_9, ld_10, \
                         md_6, md_7, md_8, md_9, md_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * ld_6[k]
                      + md_6[k];

            t_11[k] = ab_x[k] * ld_7[k]
                      + md_7[k];

            t_12[k] = ab_x[k] * ld_8[k]
                      + md_8[k];

            t_13[k] = ab_x[k] * ld_9[k]
                      + md_9[k];

            t_14[k] = ab_x[k] * ld_10[k]
                      + md_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, ld_9, ld_10, ld_11, \
                         md_11, md_21, md_22, md_23, md_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * ld_11[k]
                      + md_11[k];

            t_16[k] = ab_y[k] * ld_9[k]
                      + md_21[k];

            t_17[k] = ab_y[k] * ld_10[k]
                      + md_22[k];

            t_18[k] = ab_y[k] * ld_11[k]
                      + md_23[k];

            t_19[k] = ab_z[k] * ld_11[k]
                      + md_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, ld_12, ld_13, ld_14, ld_15, \
                         ld_16, md_12, md_13, md_14, md_15, md_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * ld_12[k]
                      + md_12[k];

            t_21[k] = ab_x[k] * ld_13[k]
                      + md_13[k];

            t_22[k] = ab_x[k] * ld_14[k]
                      + md_14[k];

            t_23[k] = ab_x[k] * ld_15[k]
                      + md_15[k];

            t_24[k] = ab_x[k] * ld_16[k]
                      + md_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, ld_15, ld_16, ld_17, \
                         md_17, md_27, md_28, md_29, md_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * ld_17[k]
                      + md_17[k];

            t_26[k] = ab_y[k] * ld_15[k]
                      + md_27[k];

            t_27[k] = ab_y[k] * ld_16[k]
                      + md_28[k];

            t_28[k] = ab_y[k] * ld_17[k]
                      + md_29[k];

            t_29[k] = ab_z[k] * ld_17[k]
                      + md_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, ld_18, ld_19, ld_20, ld_21, \
                         ld_22, md_18, md_19, md_20, md_21, md_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * ld_18[k]
                      + md_18[k];

            t_31[k] = ab_x[k] * ld_19[k]
                      + md_19[k];

            t_32[k] = ab_x[k] * ld_20[k]
                      + md_20[k];

            t_33[k] = ab_x[k] * ld_21[k]
                      + md_21[k];

            t_34[k] = ab_x[k] * ld_22[k]
                      + md_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, ld_21, ld_22, ld_23, \
                         md_23, md_39, md_40, md_41, md_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * ld_23[k]
                      + md_23[k];

            t_36[k] = ab_y[k] * ld_21[k]
                      + md_39[k];

            t_37[k] = ab_y[k] * ld_22[k]
                      + md_40[k];

            t_38[k] = ab_y[k] * ld_23[k]
                      + md_41[k];

            t_39[k] = ab_z[k] * ld_23[k]
                      + md_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, ld_24, ld_25, ld_26, ld_27, \
                         ld_28, md_24, md_25, md_26, md_27, md_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * ld_24[k]
                      + md_24[k];

            t_41[k] = ab_x[k] * ld_25[k]
                      + md_25[k];

            t_42[k] = ab_x[k] * ld_26[k]
                      + md_26[k];

            t_43[k] = ab_x[k] * ld_27[k]
                      + md_27[k];

            t_44[k] = ab_x[k] * ld_28[k]
                      + md_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, ld_27, ld_28, ld_29, \
                         md_29, md_45, md_46, md_47, md_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * ld_29[k]
                      + md_29[k];

            t_46[k] = ab_y[k] * ld_27[k]
                      + md_45[k];

            t_47[k] = ab_y[k] * ld_28[k]
                      + md_46[k];

            t_48[k] = ab_y[k] * ld_29[k]
                      + md_47[k];

            t_49[k] = ab_z[k] * ld_29[k]
                      + md_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, ld_30, ld_31, ld_32, ld_33, \
                         ld_34, md_30, md_31, md_32, md_33, md_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * ld_30[k]
                      + md_30[k];

            t_51[k] = ab_x[k] * ld_31[k]
                      + md_31[k];

            t_52[k] = ab_x[k] * ld_32[k]
                      + md_32[k];

            t_53[k] = ab_x[k] * ld_33[k]
                      + md_33[k];

            t_54[k] = ab_x[k] * ld_34[k]
                      + md_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ld_33, ld_34, ld_35, \
                         md_35, md_51, md_52, md_53, md_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * ld_35[k]
                      + md_35[k];

            t_56[k] = ab_y[k] * ld_33[k]
                      + md_51[k];

            t_57[k] = ab_y[k] * ld_34[k]
                      + md_52[k];

            t_58[k] = ab_y[k] * ld_35[k]
                      + md_53[k];

            t_59[k] = ab_z[k] * ld_35[k]
                      + md_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ld_36, ld_37, ld_38, ld_39, \
                         ld_40, md_36, md_37, md_38, md_39, md_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ld_36[k]
                      + md_36[k];

            t_61[k] = ab_x[k] * ld_37[k]
                      + md_37[k];

            t_62[k] = ab_x[k] * ld_38[k]
                      + md_38[k];

            t_63[k] = ab_x[k] * ld_39[k]
                      + md_39[k];

            t_64[k] = ab_x[k] * ld_40[k]
                      + md_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, ld_39, ld_40, ld_41, \
                         md_41, md_63, md_64, md_65, md_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * ld_41[k]
                      + md_41[k];

            t_66[k] = ab_y[k] * ld_39[k]
                      + md_63[k];

            t_67[k] = ab_y[k] * ld_40[k]
                      + md_64[k];

            t_68[k] = ab_y[k] * ld_41[k]
                      + md_65[k];

            t_69[k] = ab_z[k] * ld_41[k]
                      + md_71[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, ld_42, ld_43, ld_44, ld_45, \
                         ld_46, md_42, md_43, md_44, md_45, md_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_x[k] * ld_42[k]
                      + md_42[k];

            t_71[k] = ab_x[k] * ld_43[k]
                      + md_43[k];

            t_72[k] = ab_x[k] * ld_44[k]
                      + md_44[k];

            t_73[k] = ab_x[k] * ld_45[k]
                      + md_45[k];

            t_74[k] = ab_x[k] * ld_46[k]
                      + md_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, ld_45, ld_46, ld_47, \
                         md_47, md_69, md_70, md_71, md_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * ld_47[k]
                      + md_47[k];

            t_76[k] = ab_y[k] * ld_45[k]
                      + md_69[k];

            t_77[k] = ab_y[k] * ld_46[k]
                      + md_70[k];

            t_78[k] = ab_y[k] * ld_47[k]
                      + md_71[k];

            t_79[k] = ab_z[k] * ld_47[k]
                      + md_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, ld_48, ld_49, ld_50, ld_51, \
                         ld_52, md_48, md_49, md_50, md_51, md_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * ld_48[k]
                      + md_48[k];

            t_81[k] = ab_x[k] * ld_49[k]
                      + md_49[k];

            t_82[k] = ab_x[k] * ld_50[k]
                      + md_50[k];

            t_83[k] = ab_x[k] * ld_51[k]
                      + md_51[k];

            t_84[k] = ab_x[k] * ld_52[k]
                      + md_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, ld_51, ld_52, ld_53, \
                         md_53, md_75, md_76, md_77, md_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * ld_53[k]
                      + md_53[k];

            t_86[k] = ab_y[k] * ld_51[k]
                      + md_75[k];

            t_87[k] = ab_y[k] * ld_52[k]
                      + md_76[k];

            t_88[k] = ab_y[k] * ld_53[k]
                      + md_77[k];

            t_89[k] = ab_z[k] * ld_53[k]
                      + md_83[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ld_54, ld_55, ld_56, ld_57, \
                         ld_58, md_54, md_55, md_56, md_57, md_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * ld_54[k]
                      + md_54[k];

            t_91[k] = ab_x[k] * ld_55[k]
                      + md_55[k];

            t_92[k] = ab_x[k] * ld_56[k]
                      + md_56[k];

            t_93[k] = ab_x[k] * ld_57[k]
                      + md_57[k];

            t_94[k] = ab_x[k] * ld_58[k]
                      + md_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, ld_57, ld_58, ld_59, \
                         md_59, md_81, md_82, md_83, md_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * ld_59[k]
                      + md_59[k];

            t_96[k] = ab_y[k] * ld_57[k]
                      + md_81[k];

            t_97[k] = ab_y[k] * ld_58[k]
                      + md_82[k];

            t_98[k] = ab_y[k] * ld_59[k]
                      + md_83[k];

            t_99[k] = ab_z[k] * ld_59[k]
                      + md_89[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, ld_60, ld_61, ld_62, ld_63, \
                         ld_64, md_60, md_61, md_62, md_63, md_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_x[k] * ld_60[k]
                       + md_60[k];

            t_101[k] = ab_x[k] * ld_61[k]
                       + md_61[k];

            t_102[k] = ab_x[k] * ld_62[k]
                       + md_62[k];

            t_103[k] = ab_x[k] * ld_63[k]
                       + md_63[k];

            t_104[k] = ab_x[k] * ld_64[k]
                       + md_64[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, ld_63, ld_64, \
                         ld_65, md_65, md_93, md_94, md_95, md_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * ld_65[k]
                       + md_65[k];

            t_106[k] = ab_y[k] * ld_63[k]
                       + md_93[k];

            t_107[k] = ab_y[k] * ld_64[k]
                       + md_94[k];

            t_108[k] = ab_y[k] * ld_65[k]
                       + md_95[k];

            t_109[k] = ab_z[k] * ld_65[k]
                       + md_101[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, ld_66, ld_67, ld_68, ld_69, \
                         ld_70, md_66, md_67, md_68, md_69, md_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * ld_66[k]
                       + md_66[k];

            t_111[k] = ab_x[k] * ld_67[k]
                       + md_67[k];

            t_112[k] = ab_x[k] * ld_68[k]
                       + md_68[k];

            t_113[k] = ab_x[k] * ld_69[k]
                       + md_69[k];

            t_114[k] = ab_x[k] * ld_70[k]
                       + md_70[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, ld_69, ld_70, \
                         ld_71, md_71, md_99, md_100, md_101, md_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_x[k] * ld_71[k]
                       + md_71[k];

            t_116[k] = ab_y[k] * ld_69[k]
                       + md_99[k];

            t_117[k] = ab_y[k] * ld_70[k]
                       + md_100[k];

            t_118[k] = ab_y[k] * ld_71[k]
                       + md_101[k];

            t_119[k] = ab_z[k] * ld_71[k]
                       + md_107[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, ld_72, ld_73, ld_74, ld_75, \
                         ld_76, md_72, md_73, md_74, md_75, md_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * ld_72[k]
                       + md_72[k];

            t_121[k] = ab_x[k] * ld_73[k]
                       + md_73[k];

            t_122[k] = ab_x[k] * ld_74[k]
                       + md_74[k];

            t_123[k] = ab_x[k] * ld_75[k]
                       + md_75[k];

            t_124[k] = ab_x[k] * ld_76[k]
                       + md_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, ld_75, ld_76, \
                         ld_77, md_77, md_105, md_106, md_107, md_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * ld_77[k]
                       + md_77[k];

            t_126[k] = ab_y[k] * ld_75[k]
                       + md_105[k];

            t_127[k] = ab_y[k] * ld_76[k]
                       + md_106[k];

            t_128[k] = ab_y[k] * ld_77[k]
                       + md_107[k];

            t_129[k] = ab_z[k] * ld_77[k]
                       + md_113[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, ld_78, ld_79, ld_80, ld_81, \
                         ld_82, md_78, md_79, md_80, md_81, md_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_x[k] * ld_78[k]
                       + md_78[k];

            t_131[k] = ab_x[k] * ld_79[k]
                       + md_79[k];

            t_132[k] = ab_x[k] * ld_80[k]
                       + md_80[k];

            t_133[k] = ab_x[k] * ld_81[k]
                       + md_81[k];

            t_134[k] = ab_x[k] * ld_82[k]
                       + md_82[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, ld_81, ld_82, \
                         ld_83, md_83, md_111, md_112, md_113, md_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * ld_83[k]
                       + md_83[k];

            t_136[k] = ab_y[k] * ld_81[k]
                       + md_111[k];

            t_137[k] = ab_y[k] * ld_82[k]
                       + md_112[k];

            t_138[k] = ab_y[k] * ld_83[k]
                       + md_113[k];

            t_139[k] = ab_z[k] * ld_83[k]
                       + md_119[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, ld_84, ld_85, ld_86, ld_87, \
                         ld_88, md_84, md_85, md_86, md_87, md_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * ld_84[k]
                       + md_84[k];

            t_141[k] = ab_x[k] * ld_85[k]
                       + md_85[k];

            t_142[k] = ab_x[k] * ld_86[k]
                       + md_86[k];

            t_143[k] = ab_x[k] * ld_87[k]
                       + md_87[k];

            t_144[k] = ab_x[k] * ld_88[k]
                       + md_88[k];
        }
    }
}

static auto
compute_hrr_lf_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ld, const size_t md, const size_t ncomps,
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

        const auto *ld_87 = buffer.data(ld + 87 * ncomps + c);
        const auto *ld_88 = buffer.data(ld + 88 * ncomps + c);
        const auto *ld_89 = buffer.data(ld + 89 * ncomps + c);
        const auto *ld_90 = buffer.data(ld + 90 * ncomps + c);
        const auto *ld_91 = buffer.data(ld + 91 * ncomps + c);
        const auto *ld_92 = buffer.data(ld + 92 * ncomps + c);
        const auto *ld_93 = buffer.data(ld + 93 * ncomps + c);
        const auto *ld_94 = buffer.data(ld + 94 * ncomps + c);
        const auto *ld_95 = buffer.data(ld + 95 * ncomps + c);
        const auto *ld_96 = buffer.data(ld + 96 * ncomps + c);
        const auto *ld_97 = buffer.data(ld + 97 * ncomps + c);
        const auto *ld_98 = buffer.data(ld + 98 * ncomps + c);
        const auto *ld_99 = buffer.data(ld + 99 * ncomps + c);
        const auto *ld_100 = buffer.data(ld + 100 * ncomps + c);
        const auto *ld_101 = buffer.data(ld + 101 * ncomps + c);
        const auto *ld_102 = buffer.data(ld + 102 * ncomps + c);
        const auto *ld_103 = buffer.data(ld + 103 * ncomps + c);
        const auto *ld_104 = buffer.data(ld + 104 * ncomps + c);
        const auto *ld_105 = buffer.data(ld + 105 * ncomps + c);
        const auto *ld_106 = buffer.data(ld + 106 * ncomps + c);
        const auto *ld_107 = buffer.data(ld + 107 * ncomps + c);
        const auto *ld_108 = buffer.data(ld + 108 * ncomps + c);
        const auto *ld_109 = buffer.data(ld + 109 * ncomps + c);
        const auto *ld_110 = buffer.data(ld + 110 * ncomps + c);
        const auto *ld_111 = buffer.data(ld + 111 * ncomps + c);
        const auto *ld_112 = buffer.data(ld + 112 * ncomps + c);
        const auto *ld_113 = buffer.data(ld + 113 * ncomps + c);
        const auto *ld_114 = buffer.data(ld + 114 * ncomps + c);
        const auto *ld_115 = buffer.data(ld + 115 * ncomps + c);
        const auto *ld_116 = buffer.data(ld + 116 * ncomps + c);
        const auto *ld_117 = buffer.data(ld + 117 * ncomps + c);
        const auto *ld_118 = buffer.data(ld + 118 * ncomps + c);
        const auto *ld_119 = buffer.data(ld + 119 * ncomps + c);
        const auto *ld_120 = buffer.data(ld + 120 * ncomps + c);
        const auto *ld_121 = buffer.data(ld + 121 * ncomps + c);
        const auto *ld_122 = buffer.data(ld + 122 * ncomps + c);
        const auto *ld_123 = buffer.data(ld + 123 * ncomps + c);
        const auto *ld_124 = buffer.data(ld + 124 * ncomps + c);
        const auto *ld_125 = buffer.data(ld + 125 * ncomps + c);
        const auto *ld_126 = buffer.data(ld + 126 * ncomps + c);
        const auto *ld_127 = buffer.data(ld + 127 * ncomps + c);
        const auto *ld_128 = buffer.data(ld + 128 * ncomps + c);
        const auto *ld_129 = buffer.data(ld + 129 * ncomps + c);
        const auto *ld_130 = buffer.data(ld + 130 * ncomps + c);
        const auto *ld_131 = buffer.data(ld + 131 * ncomps + c);
        const auto *ld_132 = buffer.data(ld + 132 * ncomps + c);
        const auto *ld_133 = buffer.data(ld + 133 * ncomps + c);
        const auto *ld_134 = buffer.data(ld + 134 * ncomps + c);
        const auto *ld_135 = buffer.data(ld + 135 * ncomps + c);
        const auto *ld_136 = buffer.data(ld + 136 * ncomps + c);
        const auto *ld_137 = buffer.data(ld + 137 * ncomps + c);
        const auto *ld_138 = buffer.data(ld + 138 * ncomps + c);
        const auto *ld_139 = buffer.data(ld + 139 * ncomps + c);
        const auto *ld_140 = buffer.data(ld + 140 * ncomps + c);
        const auto *ld_141 = buffer.data(ld + 141 * ncomps + c);
        const auto *ld_142 = buffer.data(ld + 142 * ncomps + c);
        const auto *ld_143 = buffer.data(ld + 143 * ncomps + c);
        const auto *ld_144 = buffer.data(ld + 144 * ncomps + c);
        const auto *ld_145 = buffer.data(ld + 145 * ncomps + c);
        const auto *ld_146 = buffer.data(ld + 146 * ncomps + c);
        const auto *ld_147 = buffer.data(ld + 147 * ncomps + c);
        const auto *ld_148 = buffer.data(ld + 148 * ncomps + c);
        const auto *ld_149 = buffer.data(ld + 149 * ncomps + c);
        const auto *ld_150 = buffer.data(ld + 150 * ncomps + c);
        const auto *ld_151 = buffer.data(ld + 151 * ncomps + c);
        const auto *ld_152 = buffer.data(ld + 152 * ncomps + c);
        const auto *ld_153 = buffer.data(ld + 153 * ncomps + c);
        const auto *ld_154 = buffer.data(ld + 154 * ncomps + c);
        const auto *ld_155 = buffer.data(ld + 155 * ncomps + c);
        const auto *ld_156 = buffer.data(ld + 156 * ncomps + c);
        const auto *ld_157 = buffer.data(ld + 157 * ncomps + c);
        const auto *ld_158 = buffer.data(ld + 158 * ncomps + c);
        const auto *ld_159 = buffer.data(ld + 159 * ncomps + c);
        const auto *ld_160 = buffer.data(ld + 160 * ncomps + c);
        const auto *ld_161 = buffer.data(ld + 161 * ncomps + c);
        const auto *ld_162 = buffer.data(ld + 162 * ncomps + c);
        const auto *ld_163 = buffer.data(ld + 163 * ncomps + c);
        const auto *ld_164 = buffer.data(ld + 164 * ncomps + c);
        const auto *ld_165 = buffer.data(ld + 165 * ncomps + c);
        const auto *ld_166 = buffer.data(ld + 166 * ncomps + c);
        const auto *ld_167 = buffer.data(ld + 167 * ncomps + c);
        const auto *ld_168 = buffer.data(ld + 168 * ncomps + c);
        const auto *ld_169 = buffer.data(ld + 169 * ncomps + c);
        const auto *ld_170 = buffer.data(ld + 170 * ncomps + c);
        const auto *ld_171 = buffer.data(ld + 171 * ncomps + c);
        const auto *ld_172 = buffer.data(ld + 172 * ncomps + c);
        const auto *ld_173 = buffer.data(ld + 173 * ncomps + c);

        const auto *md_89 = buffer.data(md + 89 * ncomps + c);
        const auto *md_90 = buffer.data(md + 90 * ncomps + c);
        const auto *md_91 = buffer.data(md + 91 * ncomps + c);
        const auto *md_92 = buffer.data(md + 92 * ncomps + c);
        const auto *md_93 = buffer.data(md + 93 * ncomps + c);
        const auto *md_94 = buffer.data(md + 94 * ncomps + c);
        const auto *md_95 = buffer.data(md + 95 * ncomps + c);
        const auto *md_96 = buffer.data(md + 96 * ncomps + c);
        const auto *md_97 = buffer.data(md + 97 * ncomps + c);
        const auto *md_98 = buffer.data(md + 98 * ncomps + c);
        const auto *md_99 = buffer.data(md + 99 * ncomps + c);
        const auto *md_100 = buffer.data(md + 100 * ncomps + c);
        const auto *md_101 = buffer.data(md + 101 * ncomps + c);
        const auto *md_102 = buffer.data(md + 102 * ncomps + c);
        const auto *md_103 = buffer.data(md + 103 * ncomps + c);
        const auto *md_104 = buffer.data(md + 104 * ncomps + c);
        const auto *md_105 = buffer.data(md + 105 * ncomps + c);
        const auto *md_106 = buffer.data(md + 106 * ncomps + c);
        const auto *md_107 = buffer.data(md + 107 * ncomps + c);
        const auto *md_108 = buffer.data(md + 108 * ncomps + c);
        const auto *md_109 = buffer.data(md + 109 * ncomps + c);
        const auto *md_110 = buffer.data(md + 110 * ncomps + c);
        const auto *md_111 = buffer.data(md + 111 * ncomps + c);
        const auto *md_112 = buffer.data(md + 112 * ncomps + c);
        const auto *md_113 = buffer.data(md + 113 * ncomps + c);
        const auto *md_114 = buffer.data(md + 114 * ncomps + c);
        const auto *md_115 = buffer.data(md + 115 * ncomps + c);
        const auto *md_116 = buffer.data(md + 116 * ncomps + c);
        const auto *md_117 = buffer.data(md + 117 * ncomps + c);
        const auto *md_118 = buffer.data(md + 118 * ncomps + c);
        const auto *md_119 = buffer.data(md + 119 * ncomps + c);
        const auto *md_120 = buffer.data(md + 120 * ncomps + c);
        const auto *md_121 = buffer.data(md + 121 * ncomps + c);
        const auto *md_122 = buffer.data(md + 122 * ncomps + c);
        const auto *md_123 = buffer.data(md + 123 * ncomps + c);
        const auto *md_124 = buffer.data(md + 124 * ncomps + c);
        const auto *md_125 = buffer.data(md + 125 * ncomps + c);
        const auto *md_126 = buffer.data(md + 126 * ncomps + c);
        const auto *md_127 = buffer.data(md + 127 * ncomps + c);
        const auto *md_128 = buffer.data(md + 128 * ncomps + c);
        const auto *md_129 = buffer.data(md + 129 * ncomps + c);
        const auto *md_130 = buffer.data(md + 130 * ncomps + c);
        const auto *md_131 = buffer.data(md + 131 * ncomps + c);
        const auto *md_132 = buffer.data(md + 132 * ncomps + c);
        const auto *md_133 = buffer.data(md + 133 * ncomps + c);
        const auto *md_134 = buffer.data(md + 134 * ncomps + c);
        const auto *md_135 = buffer.data(md + 135 * ncomps + c);
        const auto *md_136 = buffer.data(md + 136 * ncomps + c);
        const auto *md_137 = buffer.data(md + 137 * ncomps + c);
        const auto *md_138 = buffer.data(md + 138 * ncomps + c);
        const auto *md_139 = buffer.data(md + 139 * ncomps + c);
        const auto *md_140 = buffer.data(md + 140 * ncomps + c);
        const auto *md_141 = buffer.data(md + 141 * ncomps + c);
        const auto *md_142 = buffer.data(md + 142 * ncomps + c);
        const auto *md_143 = buffer.data(md + 143 * ncomps + c);
        const auto *md_144 = buffer.data(md + 144 * ncomps + c);
        const auto *md_145 = buffer.data(md + 145 * ncomps + c);
        const auto *md_146 = buffer.data(md + 146 * ncomps + c);
        const auto *md_147 = buffer.data(md + 147 * ncomps + c);
        const auto *md_148 = buffer.data(md + 148 * ncomps + c);
        const auto *md_149 = buffer.data(md + 149 * ncomps + c);
        const auto *md_150 = buffer.data(md + 150 * ncomps + c);
        const auto *md_151 = buffer.data(md + 151 * ncomps + c);
        const auto *md_152 = buffer.data(md + 152 * ncomps + c);
        const auto *md_153 = buffer.data(md + 153 * ncomps + c);
        const auto *md_154 = buffer.data(md + 154 * ncomps + c);
        const auto *md_155 = buffer.data(md + 155 * ncomps + c);
        const auto *md_156 = buffer.data(md + 156 * ncomps + c);
        const auto *md_157 = buffer.data(md + 157 * ncomps + c);
        const auto *md_158 = buffer.data(md + 158 * ncomps + c);
        const auto *md_159 = buffer.data(md + 159 * ncomps + c);
        const auto *md_160 = buffer.data(md + 160 * ncomps + c);
        const auto *md_161 = buffer.data(md + 161 * ncomps + c);
        const auto *md_162 = buffer.data(md + 162 * ncomps + c);
        const auto *md_163 = buffer.data(md + 163 * ncomps + c);
        const auto *md_164 = buffer.data(md + 164 * ncomps + c);
        const auto *md_165 = buffer.data(md + 165 * ncomps + c);
        const auto *md_166 = buffer.data(md + 166 * ncomps + c);
        const auto *md_167 = buffer.data(md + 167 * ncomps + c);
        const auto *md_168 = buffer.data(md + 168 * ncomps + c);
        const auto *md_169 = buffer.data(md + 169 * ncomps + c);
        const auto *md_170 = buffer.data(md + 170 * ncomps + c);
        const auto *md_171 = buffer.data(md + 171 * ncomps + c);
        const auto *md_172 = buffer.data(md + 172 * ncomps + c);
        const auto *md_173 = buffer.data(md + 173 * ncomps + c);
        const auto *md_177 = buffer.data(md + 177 * ncomps + c);
        const auto *md_178 = buffer.data(md + 178 * ncomps + c);
        const auto *md_179 = buffer.data(md + 179 * ncomps + c);
        const auto *md_183 = buffer.data(md + 183 * ncomps + c);
        const auto *md_184 = buffer.data(md + 184 * ncomps + c);
        const auto *md_185 = buffer.data(md + 185 * ncomps + c);
        const auto *md_189 = buffer.data(md + 189 * ncomps + c);
        const auto *md_190 = buffer.data(md + 190 * ncomps + c);
        const auto *md_191 = buffer.data(md + 191 * ncomps + c);
        const auto *md_195 = buffer.data(md + 195 * ncomps + c);
        const auto *md_196 = buffer.data(md + 196 * ncomps + c);
        const auto *md_197 = buffer.data(md + 197 * ncomps + c);
        const auto *md_201 = buffer.data(md + 201 * ncomps + c);
        const auto *md_202 = buffer.data(md + 202 * ncomps + c);
        const auto *md_203 = buffer.data(md + 203 * ncomps + c);
        const auto *md_207 = buffer.data(md + 207 * ncomps + c);
        const auto *md_208 = buffer.data(md + 208 * ncomps + c);
        const auto *md_209 = buffer.data(md + 209 * ncomps + c);
        const auto *md_215 = buffer.data(md + 215 * ncomps + c);
        const auto *md_219 = buffer.data(md + 219 * ncomps + c);
        const auto *md_220 = buffer.data(md + 220 * ncomps + c);
        const auto *md_221 = buffer.data(md + 221 * ncomps + c);
        const auto *md_227 = buffer.data(md + 227 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, ld_87, ld_88, \
                         ld_89, md_89, md_117, md_118, md_119, md_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * ld_89[k]
                       + md_89[k];

            t_146[k] = ab_y[k] * ld_87[k]
                       + md_117[k];

            t_147[k] = ab_y[k] * ld_88[k]
                       + md_118[k];

            t_148[k] = ab_y[k] * ld_89[k]
                       + md_119[k];

            t_149[k] = ab_z[k] * ld_89[k]
                       + md_125[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, ld_90, ld_91, ld_92, ld_93, \
                         ld_94, md_90, md_91, md_92, md_93, md_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ld_90[k]
                       + md_90[k];

            t_151[k] = ab_x[k] * ld_91[k]
                       + md_91[k];

            t_152[k] = ab_x[k] * ld_92[k]
                       + md_92[k];

            t_153[k] = ab_x[k] * ld_93[k]
                       + md_93[k];

            t_154[k] = ab_x[k] * ld_94[k]
                       + md_94[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, ld_93, ld_94, \
                         ld_95, md_95, md_129, md_130, md_131, md_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * ld_95[k]
                       + md_95[k];

            t_156[k] = ab_y[k] * ld_93[k]
                       + md_129[k];

            t_157[k] = ab_y[k] * ld_94[k]
                       + md_130[k];

            t_158[k] = ab_y[k] * ld_95[k]
                       + md_131[k];

            t_159[k] = ab_z[k] * ld_95[k]
                       + md_137[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, ld_96, ld_97, ld_98, ld_99, \
                         ld_100, md_96, md_97, md_98, md_99, md_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * ld_96[k]
                       + md_96[k];

            t_161[k] = ab_x[k] * ld_97[k]
                       + md_97[k];

            t_162[k] = ab_x[k] * ld_98[k]
                       + md_98[k];

            t_163[k] = ab_x[k] * ld_99[k]
                       + md_99[k];

            t_164[k] = ab_x[k] * ld_100[k]
                       + md_100[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, ld_99, ld_100, \
                         ld_101, md_101, md_135, md_136, md_137, \
                         md_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * ld_101[k]
                       + md_101[k];

            t_166[k] = ab_y[k] * ld_99[k]
                       + md_135[k];

            t_167[k] = ab_y[k] * ld_100[k]
                       + md_136[k];

            t_168[k] = ab_y[k] * ld_101[k]
                       + md_137[k];

            t_169[k] = ab_z[k] * ld_101[k]
                       + md_143[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, ld_102, ld_103, ld_104, \
                         ld_105, ld_106, md_102, md_103, md_104, md_105, \
                         md_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * ld_102[k]
                       + md_102[k];

            t_171[k] = ab_x[k] * ld_103[k]
                       + md_103[k];

            t_172[k] = ab_x[k] * ld_104[k]
                       + md_104[k];

            t_173[k] = ab_x[k] * ld_105[k]
                       + md_105[k];

            t_174[k] = ab_x[k] * ld_106[k]
                       + md_106[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, ld_105, ld_106, \
                         ld_107, md_107, md_141, md_142, md_143, \
                         md_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * ld_107[k]
                       + md_107[k];

            t_176[k] = ab_y[k] * ld_105[k]
                       + md_141[k];

            t_177[k] = ab_y[k] * ld_106[k]
                       + md_142[k];

            t_178[k] = ab_y[k] * ld_107[k]
                       + md_143[k];

            t_179[k] = ab_z[k] * ld_107[k]
                       + md_149[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, ld_108, ld_109, ld_110, \
                         ld_111, ld_112, md_108, md_109, md_110, md_111, \
                         md_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * ld_108[k]
                       + md_108[k];

            t_181[k] = ab_x[k] * ld_109[k]
                       + md_109[k];

            t_182[k] = ab_x[k] * ld_110[k]
                       + md_110[k];

            t_183[k] = ab_x[k] * ld_111[k]
                       + md_111[k];

            t_184[k] = ab_x[k] * ld_112[k]
                       + md_112[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, ld_111, ld_112, \
                         ld_113, md_113, md_147, md_148, md_149, \
                         md_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * ld_113[k]
                       + md_113[k];

            t_186[k] = ab_y[k] * ld_111[k]
                       + md_147[k];

            t_187[k] = ab_y[k] * ld_112[k]
                       + md_148[k];

            t_188[k] = ab_y[k] * ld_113[k]
                       + md_149[k];

            t_189[k] = ab_z[k] * ld_113[k]
                       + md_155[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, ld_114, ld_115, ld_116, \
                         ld_117, ld_118, md_114, md_115, md_116, md_117, \
                         md_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_x[k] * ld_114[k]
                       + md_114[k];

            t_191[k] = ab_x[k] * ld_115[k]
                       + md_115[k];

            t_192[k] = ab_x[k] * ld_116[k]
                       + md_116[k];

            t_193[k] = ab_x[k] * ld_117[k]
                       + md_117[k];

            t_194[k] = ab_x[k] * ld_118[k]
                       + md_118[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, ld_117, ld_118, \
                         ld_119, md_119, md_153, md_154, md_155, \
                         md_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * ld_119[k]
                       + md_119[k];

            t_196[k] = ab_y[k] * ld_117[k]
                       + md_153[k];

            t_197[k] = ab_y[k] * ld_118[k]
                       + md_154[k];

            t_198[k] = ab_y[k] * ld_119[k]
                       + md_155[k];

            t_199[k] = ab_z[k] * ld_119[k]
                       + md_161[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, ld_120, ld_121, ld_122, \
                         ld_123, ld_124, md_120, md_121, md_122, md_123, \
                         md_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * ld_120[k]
                       + md_120[k];

            t_201[k] = ab_x[k] * ld_121[k]
                       + md_121[k];

            t_202[k] = ab_x[k] * ld_122[k]
                       + md_122[k];

            t_203[k] = ab_x[k] * ld_123[k]
                       + md_123[k];

            t_204[k] = ab_x[k] * ld_124[k]
                       + md_124[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, ld_123, ld_124, \
                         ld_125, md_125, md_159, md_160, md_161, \
                         md_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_x[k] * ld_125[k]
                       + md_125[k];

            t_206[k] = ab_y[k] * ld_123[k]
                       + md_159[k];

            t_207[k] = ab_y[k] * ld_124[k]
                       + md_160[k];

            t_208[k] = ab_y[k] * ld_125[k]
                       + md_161[k];

            t_209[k] = ab_z[k] * ld_125[k]
                       + md_167[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, ld_126, ld_127, ld_128, \
                         ld_129, ld_130, md_126, md_127, md_128, md_129, \
                         md_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * ld_126[k]
                       + md_126[k];

            t_211[k] = ab_x[k] * ld_127[k]
                       + md_127[k];

            t_212[k] = ab_x[k] * ld_128[k]
                       + md_128[k];

            t_213[k] = ab_x[k] * ld_129[k]
                       + md_129[k];

            t_214[k] = ab_x[k] * ld_130[k]
                       + md_130[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ab_y, ab_z, ld_129, ld_130, \
                         ld_131, md_131, md_171, md_172, md_173, \
                         md_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * ld_131[k]
                       + md_131[k];

            t_216[k] = ab_y[k] * ld_129[k]
                       + md_171[k];

            t_217[k] = ab_y[k] * ld_130[k]
                       + md_172[k];

            t_218[k] = ab_y[k] * ld_131[k]
                       + md_173[k];

            t_219[k] = ab_z[k] * ld_131[k]
                       + md_179[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, ld_132, ld_133, ld_134, \
                         ld_135, ld_136, md_132, md_133, md_134, md_135, \
                         md_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * ld_132[k]
                       + md_132[k];

            t_221[k] = ab_x[k] * ld_133[k]
                       + md_133[k];

            t_222[k] = ab_x[k] * ld_134[k]
                       + md_134[k];

            t_223[k] = ab_x[k] * ld_135[k]
                       + md_135[k];

            t_224[k] = ab_x[k] * ld_136[k]
                       + md_136[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, ab_y, ab_z, ld_135, ld_136, \
                         ld_137, md_137, md_177, md_178, md_179, \
                         md_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * ld_137[k]
                       + md_137[k];

            t_226[k] = ab_y[k] * ld_135[k]
                       + md_177[k];

            t_227[k] = ab_y[k] * ld_136[k]
                       + md_178[k];

            t_228[k] = ab_y[k] * ld_137[k]
                       + md_179[k];

            t_229[k] = ab_z[k] * ld_137[k]
                       + md_185[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, ld_138, ld_139, ld_140, \
                         ld_141, ld_142, md_138, md_139, md_140, md_141, \
                         md_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * ld_138[k]
                       + md_138[k];

            t_231[k] = ab_x[k] * ld_139[k]
                       + md_139[k];

            t_232[k] = ab_x[k] * ld_140[k]
                       + md_140[k];

            t_233[k] = ab_x[k] * ld_141[k]
                       + md_141[k];

            t_234[k] = ab_x[k] * ld_142[k]
                       + md_142[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, ab_y, ab_z, ld_141, ld_142, \
                         ld_143, md_143, md_183, md_184, md_185, \
                         md_191 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_x[k] * ld_143[k]
                       + md_143[k];

            t_236[k] = ab_y[k] * ld_141[k]
                       + md_183[k];

            t_237[k] = ab_y[k] * ld_142[k]
                       + md_184[k];

            t_238[k] = ab_y[k] * ld_143[k]
                       + md_185[k];

            t_239[k] = ab_z[k] * ld_143[k]
                       + md_191[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, ld_144, ld_145, ld_146, \
                         ld_147, ld_148, md_144, md_145, md_146, md_147, \
                         md_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * ld_144[k]
                       + md_144[k];

            t_241[k] = ab_x[k] * ld_145[k]
                       + md_145[k];

            t_242[k] = ab_x[k] * ld_146[k]
                       + md_146[k];

            t_243[k] = ab_x[k] * ld_147[k]
                       + md_147[k];

            t_244[k] = ab_x[k] * ld_148[k]
                       + md_148[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, ab_y, ab_z, ld_147, ld_148, \
                         ld_149, md_149, md_189, md_190, md_191, \
                         md_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * ld_149[k]
                       + md_149[k];

            t_246[k] = ab_y[k] * ld_147[k]
                       + md_189[k];

            t_247[k] = ab_y[k] * ld_148[k]
                       + md_190[k];

            t_248[k] = ab_y[k] * ld_149[k]
                       + md_191[k];

            t_249[k] = ab_z[k] * ld_149[k]
                       + md_197[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, ld_150, ld_151, ld_152, \
                         ld_153, ld_154, md_150, md_151, md_152, md_153, \
                         md_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_x[k] * ld_150[k]
                       + md_150[k];

            t_251[k] = ab_x[k] * ld_151[k]
                       + md_151[k];

            t_252[k] = ab_x[k] * ld_152[k]
                       + md_152[k];

            t_253[k] = ab_x[k] * ld_153[k]
                       + md_153[k];

            t_254[k] = ab_x[k] * ld_154[k]
                       + md_154[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, ab_y, ab_z, ld_153, ld_154, \
                         ld_155, md_155, md_195, md_196, md_197, \
                         md_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * ld_155[k]
                       + md_155[k];

            t_256[k] = ab_y[k] * ld_153[k]
                       + md_195[k];

            t_257[k] = ab_y[k] * ld_154[k]
                       + md_196[k];

            t_258[k] = ab_y[k] * ld_155[k]
                       + md_197[k];

            t_259[k] = ab_z[k] * ld_155[k]
                       + md_203[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, ld_156, ld_157, ld_158, \
                         ld_159, ld_160, md_156, md_157, md_158, md_159, \
                         md_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * ld_156[k]
                       + md_156[k];

            t_261[k] = ab_x[k] * ld_157[k]
                       + md_157[k];

            t_262[k] = ab_x[k] * ld_158[k]
                       + md_158[k];

            t_263[k] = ab_x[k] * ld_159[k]
                       + md_159[k];

            t_264[k] = ab_x[k] * ld_160[k]
                       + md_160[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, ld_159, ld_160, \
                         ld_161, md_161, md_201, md_202, md_203, \
                         md_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_x[k] * ld_161[k]
                       + md_161[k];

            t_266[k] = ab_y[k] * ld_159[k]
                       + md_201[k];

            t_267[k] = ab_y[k] * ld_160[k]
                       + md_202[k];

            t_268[k] = ab_y[k] * ld_161[k]
                       + md_203[k];

            t_269[k] = ab_z[k] * ld_161[k]
                       + md_209[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, ld_162, ld_163, ld_164, \
                         ld_165, ld_166, md_162, md_163, md_164, md_165, \
                         md_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * ld_162[k]
                       + md_162[k];

            t_271[k] = ab_x[k] * ld_163[k]
                       + md_163[k];

            t_272[k] = ab_x[k] * ld_164[k]
                       + md_164[k];

            t_273[k] = ab_x[k] * ld_165[k]
                       + md_165[k];

            t_274[k] = ab_x[k] * ld_166[k]
                       + md_166[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, ab_y, ab_z, ld_165, ld_166, \
                         ld_167, md_167, md_207, md_208, md_209, \
                         md_215 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * ld_167[k]
                       + md_167[k];

            t_276[k] = ab_y[k] * ld_165[k]
                       + md_207[k];

            t_277[k] = ab_y[k] * ld_166[k]
                       + md_208[k];

            t_278[k] = ab_y[k] * ld_167[k]
                       + md_209[k];

            t_279[k] = ab_z[k] * ld_167[k]
                       + md_215[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, ld_168, ld_169, ld_170, \
                         ld_171, ld_172, md_168, md_169, md_170, md_171, \
                         md_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_x[k] * ld_168[k]
                       + md_168[k];

            t_281[k] = ab_x[k] * ld_169[k]
                       + md_169[k];

            t_282[k] = ab_x[k] * ld_170[k]
                       + md_170[k];

            t_283[k] = ab_x[k] * ld_171[k]
                       + md_171[k];

            t_284[k] = ab_x[k] * ld_172[k]
                       + md_172[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, ab_y, ab_z, ld_171, ld_172, \
                         ld_173, md_173, md_219, md_220, md_221, \
                         md_227 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * ld_173[k]
                       + md_173[k];

            t_286[k] = ab_y[k] * ld_171[k]
                       + md_219[k];

            t_287[k] = ab_y[k] * ld_172[k]
                       + md_220[k];

            t_288[k] = ab_y[k] * ld_173[k]
                       + md_221[k];

            t_289[k] = ab_z[k] * ld_173[k]
                       + md_227[k];
        }
    }
}

static auto
compute_hrr_lf_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ld, const size_t md, const size_t ncomps,
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

        const auto *ld_174 = buffer.data(ld + 174 * ncomps + c);
        const auto *ld_175 = buffer.data(ld + 175 * ncomps + c);
        const auto *ld_176 = buffer.data(ld + 176 * ncomps + c);
        const auto *ld_177 = buffer.data(ld + 177 * ncomps + c);
        const auto *ld_178 = buffer.data(ld + 178 * ncomps + c);
        const auto *ld_179 = buffer.data(ld + 179 * ncomps + c);
        const auto *ld_180 = buffer.data(ld + 180 * ncomps + c);
        const auto *ld_181 = buffer.data(ld + 181 * ncomps + c);
        const auto *ld_182 = buffer.data(ld + 182 * ncomps + c);
        const auto *ld_183 = buffer.data(ld + 183 * ncomps + c);
        const auto *ld_184 = buffer.data(ld + 184 * ncomps + c);
        const auto *ld_185 = buffer.data(ld + 185 * ncomps + c);
        const auto *ld_186 = buffer.data(ld + 186 * ncomps + c);
        const auto *ld_187 = buffer.data(ld + 187 * ncomps + c);
        const auto *ld_188 = buffer.data(ld + 188 * ncomps + c);
        const auto *ld_189 = buffer.data(ld + 189 * ncomps + c);
        const auto *ld_190 = buffer.data(ld + 190 * ncomps + c);
        const auto *ld_191 = buffer.data(ld + 191 * ncomps + c);
        const auto *ld_192 = buffer.data(ld + 192 * ncomps + c);
        const auto *ld_193 = buffer.data(ld + 193 * ncomps + c);
        const auto *ld_194 = buffer.data(ld + 194 * ncomps + c);
        const auto *ld_195 = buffer.data(ld + 195 * ncomps + c);
        const auto *ld_196 = buffer.data(ld + 196 * ncomps + c);
        const auto *ld_197 = buffer.data(ld + 197 * ncomps + c);
        const auto *ld_198 = buffer.data(ld + 198 * ncomps + c);
        const auto *ld_199 = buffer.data(ld + 199 * ncomps + c);
        const auto *ld_200 = buffer.data(ld + 200 * ncomps + c);
        const auto *ld_201 = buffer.data(ld + 201 * ncomps + c);
        const auto *ld_202 = buffer.data(ld + 202 * ncomps + c);
        const auto *ld_203 = buffer.data(ld + 203 * ncomps + c);
        const auto *ld_204 = buffer.data(ld + 204 * ncomps + c);
        const auto *ld_205 = buffer.data(ld + 205 * ncomps + c);
        const auto *ld_206 = buffer.data(ld + 206 * ncomps + c);
        const auto *ld_207 = buffer.data(ld + 207 * ncomps + c);
        const auto *ld_208 = buffer.data(ld + 208 * ncomps + c);
        const auto *ld_209 = buffer.data(ld + 209 * ncomps + c);
        const auto *ld_210 = buffer.data(ld + 210 * ncomps + c);
        const auto *ld_211 = buffer.data(ld + 211 * ncomps + c);
        const auto *ld_212 = buffer.data(ld + 212 * ncomps + c);
        const auto *ld_213 = buffer.data(ld + 213 * ncomps + c);
        const auto *ld_214 = buffer.data(ld + 214 * ncomps + c);
        const auto *ld_215 = buffer.data(ld + 215 * ncomps + c);
        const auto *ld_216 = buffer.data(ld + 216 * ncomps + c);
        const auto *ld_217 = buffer.data(ld + 217 * ncomps + c);
        const auto *ld_218 = buffer.data(ld + 218 * ncomps + c);
        const auto *ld_219 = buffer.data(ld + 219 * ncomps + c);
        const auto *ld_220 = buffer.data(ld + 220 * ncomps + c);
        const auto *ld_221 = buffer.data(ld + 221 * ncomps + c);
        const auto *ld_222 = buffer.data(ld + 222 * ncomps + c);
        const auto *ld_223 = buffer.data(ld + 223 * ncomps + c);
        const auto *ld_224 = buffer.data(ld + 224 * ncomps + c);
        const auto *ld_225 = buffer.data(ld + 225 * ncomps + c);
        const auto *ld_226 = buffer.data(ld + 226 * ncomps + c);
        const auto *ld_227 = buffer.data(ld + 227 * ncomps + c);
        const auto *ld_228 = buffer.data(ld + 228 * ncomps + c);
        const auto *ld_229 = buffer.data(ld + 229 * ncomps + c);
        const auto *ld_230 = buffer.data(ld + 230 * ncomps + c);
        const auto *ld_231 = buffer.data(ld + 231 * ncomps + c);
        const auto *ld_232 = buffer.data(ld + 232 * ncomps + c);
        const auto *ld_233 = buffer.data(ld + 233 * ncomps + c);
        const auto *ld_234 = buffer.data(ld + 234 * ncomps + c);
        const auto *ld_235 = buffer.data(ld + 235 * ncomps + c);
        const auto *ld_236 = buffer.data(ld + 236 * ncomps + c);
        const auto *ld_237 = buffer.data(ld + 237 * ncomps + c);
        const auto *ld_238 = buffer.data(ld + 238 * ncomps + c);
        const auto *ld_239 = buffer.data(ld + 239 * ncomps + c);
        const auto *ld_240 = buffer.data(ld + 240 * ncomps + c);
        const auto *ld_241 = buffer.data(ld + 241 * ncomps + c);
        const auto *ld_242 = buffer.data(ld + 242 * ncomps + c);
        const auto *ld_243 = buffer.data(ld + 243 * ncomps + c);
        const auto *ld_244 = buffer.data(ld + 244 * ncomps + c);
        const auto *ld_245 = buffer.data(ld + 245 * ncomps + c);
        const auto *ld_246 = buffer.data(ld + 246 * ncomps + c);
        const auto *ld_247 = buffer.data(ld + 247 * ncomps + c);
        const auto *ld_248 = buffer.data(ld + 248 * ncomps + c);
        const auto *ld_249 = buffer.data(ld + 249 * ncomps + c);
        const auto *ld_250 = buffer.data(ld + 250 * ncomps + c);
        const auto *ld_251 = buffer.data(ld + 251 * ncomps + c);
        const auto *ld_252 = buffer.data(ld + 252 * ncomps + c);
        const auto *ld_253 = buffer.data(ld + 253 * ncomps + c);
        const auto *ld_254 = buffer.data(ld + 254 * ncomps + c);
        const auto *ld_255 = buffer.data(ld + 255 * ncomps + c);
        const auto *ld_256 = buffer.data(ld + 256 * ncomps + c);
        const auto *ld_257 = buffer.data(ld + 257 * ncomps + c);
        const auto *ld_258 = buffer.data(ld + 258 * ncomps + c);
        const auto *ld_259 = buffer.data(ld + 259 * ncomps + c);
        const auto *ld_260 = buffer.data(ld + 260 * ncomps + c);
        const auto *ld_261 = buffer.data(ld + 261 * ncomps + c);
        const auto *ld_262 = buffer.data(ld + 262 * ncomps + c);

        const auto *md_174 = buffer.data(md + 174 * ncomps + c);
        const auto *md_175 = buffer.data(md + 175 * ncomps + c);
        const auto *md_176 = buffer.data(md + 176 * ncomps + c);
        const auto *md_177 = buffer.data(md + 177 * ncomps + c);
        const auto *md_178 = buffer.data(md + 178 * ncomps + c);
        const auto *md_179 = buffer.data(md + 179 * ncomps + c);
        const auto *md_180 = buffer.data(md + 180 * ncomps + c);
        const auto *md_181 = buffer.data(md + 181 * ncomps + c);
        const auto *md_182 = buffer.data(md + 182 * ncomps + c);
        const auto *md_183 = buffer.data(md + 183 * ncomps + c);
        const auto *md_184 = buffer.data(md + 184 * ncomps + c);
        const auto *md_185 = buffer.data(md + 185 * ncomps + c);
        const auto *md_186 = buffer.data(md + 186 * ncomps + c);
        const auto *md_187 = buffer.data(md + 187 * ncomps + c);
        const auto *md_188 = buffer.data(md + 188 * ncomps + c);
        const auto *md_189 = buffer.data(md + 189 * ncomps + c);
        const auto *md_190 = buffer.data(md + 190 * ncomps + c);
        const auto *md_191 = buffer.data(md + 191 * ncomps + c);
        const auto *md_192 = buffer.data(md + 192 * ncomps + c);
        const auto *md_193 = buffer.data(md + 193 * ncomps + c);
        const auto *md_194 = buffer.data(md + 194 * ncomps + c);
        const auto *md_195 = buffer.data(md + 195 * ncomps + c);
        const auto *md_196 = buffer.data(md + 196 * ncomps + c);
        const auto *md_197 = buffer.data(md + 197 * ncomps + c);
        const auto *md_198 = buffer.data(md + 198 * ncomps + c);
        const auto *md_199 = buffer.data(md + 199 * ncomps + c);
        const auto *md_200 = buffer.data(md + 200 * ncomps + c);
        const auto *md_201 = buffer.data(md + 201 * ncomps + c);
        const auto *md_202 = buffer.data(md + 202 * ncomps + c);
        const auto *md_203 = buffer.data(md + 203 * ncomps + c);
        const auto *md_204 = buffer.data(md + 204 * ncomps + c);
        const auto *md_205 = buffer.data(md + 205 * ncomps + c);
        const auto *md_206 = buffer.data(md + 206 * ncomps + c);
        const auto *md_207 = buffer.data(md + 207 * ncomps + c);
        const auto *md_208 = buffer.data(md + 208 * ncomps + c);
        const auto *md_209 = buffer.data(md + 209 * ncomps + c);
        const auto *md_210 = buffer.data(md + 210 * ncomps + c);
        const auto *md_211 = buffer.data(md + 211 * ncomps + c);
        const auto *md_212 = buffer.data(md + 212 * ncomps + c);
        const auto *md_213 = buffer.data(md + 213 * ncomps + c);
        const auto *md_214 = buffer.data(md + 214 * ncomps + c);
        const auto *md_215 = buffer.data(md + 215 * ncomps + c);
        const auto *md_216 = buffer.data(md + 216 * ncomps + c);
        const auto *md_217 = buffer.data(md + 217 * ncomps + c);
        const auto *md_218 = buffer.data(md + 218 * ncomps + c);
        const auto *md_219 = buffer.data(md + 219 * ncomps + c);
        const auto *md_220 = buffer.data(md + 220 * ncomps + c);
        const auto *md_221 = buffer.data(md + 221 * ncomps + c);
        const auto *md_222 = buffer.data(md + 222 * ncomps + c);
        const auto *md_223 = buffer.data(md + 223 * ncomps + c);
        const auto *md_224 = buffer.data(md + 224 * ncomps + c);
        const auto *md_225 = buffer.data(md + 225 * ncomps + c);
        const auto *md_226 = buffer.data(md + 226 * ncomps + c);
        const auto *md_227 = buffer.data(md + 227 * ncomps + c);
        const auto *md_228 = buffer.data(md + 228 * ncomps + c);
        const auto *md_229 = buffer.data(md + 229 * ncomps + c);
        const auto *md_230 = buffer.data(md + 230 * ncomps + c);
        const auto *md_231 = buffer.data(md + 231 * ncomps + c);
        const auto *md_232 = buffer.data(md + 232 * ncomps + c);
        const auto *md_233 = buffer.data(md + 233 * ncomps + c);
        const auto *md_234 = buffer.data(md + 234 * ncomps + c);
        const auto *md_235 = buffer.data(md + 235 * ncomps + c);
        const auto *md_236 = buffer.data(md + 236 * ncomps + c);
        const auto *md_237 = buffer.data(md + 237 * ncomps + c);
        const auto *md_238 = buffer.data(md + 238 * ncomps + c);
        const auto *md_239 = buffer.data(md + 239 * ncomps + c);
        const auto *md_240 = buffer.data(md + 240 * ncomps + c);
        const auto *md_241 = buffer.data(md + 241 * ncomps + c);
        const auto *md_242 = buffer.data(md + 242 * ncomps + c);
        const auto *md_243 = buffer.data(md + 243 * ncomps + c);
        const auto *md_244 = buffer.data(md + 244 * ncomps + c);
        const auto *md_245 = buffer.data(md + 245 * ncomps + c);
        const auto *md_246 = buffer.data(md + 246 * ncomps + c);
        const auto *md_247 = buffer.data(md + 247 * ncomps + c);
        const auto *md_248 = buffer.data(md + 248 * ncomps + c);
        const auto *md_249 = buffer.data(md + 249 * ncomps + c);
        const auto *md_250 = buffer.data(md + 250 * ncomps + c);
        const auto *md_251 = buffer.data(md + 251 * ncomps + c);
        const auto *md_252 = buffer.data(md + 252 * ncomps + c);
        const auto *md_253 = buffer.data(md + 253 * ncomps + c);
        const auto *md_254 = buffer.data(md + 254 * ncomps + c);
        const auto *md_255 = buffer.data(md + 255 * ncomps + c);
        const auto *md_256 = buffer.data(md + 256 * ncomps + c);
        const auto *md_257 = buffer.data(md + 257 * ncomps + c);
        const auto *md_258 = buffer.data(md + 258 * ncomps + c);
        const auto *md_259 = buffer.data(md + 259 * ncomps + c);
        const auto *md_260 = buffer.data(md + 260 * ncomps + c);
        const auto *md_261 = buffer.data(md + 261 * ncomps + c);
        const auto *md_262 = buffer.data(md + 262 * ncomps + c);
        const auto *md_263 = buffer.data(md + 263 * ncomps + c);
        const auto *md_269 = buffer.data(md + 269 * ncomps + c);
        const auto *md_273 = buffer.data(md + 273 * ncomps + c);
        const auto *md_274 = buffer.data(md + 274 * ncomps + c);
        const auto *md_275 = buffer.data(md + 275 * ncomps + c);
        const auto *md_279 = buffer.data(md + 279 * ncomps + c);
        const auto *md_280 = buffer.data(md + 280 * ncomps + c);
        const auto *md_281 = buffer.data(md + 281 * ncomps + c);
        const auto *md_285 = buffer.data(md + 285 * ncomps + c);
        const auto *md_286 = buffer.data(md + 286 * ncomps + c);
        const auto *md_287 = buffer.data(md + 287 * ncomps + c);
        const auto *md_291 = buffer.data(md + 291 * ncomps + c);
        const auto *md_292 = buffer.data(md + 292 * ncomps + c);
        const auto *md_293 = buffer.data(md + 293 * ncomps + c);
        const auto *md_297 = buffer.data(md + 297 * ncomps + c);
        const auto *md_298 = buffer.data(md + 298 * ncomps + c);
        const auto *md_299 = buffer.data(md + 299 * ncomps + c);
        const auto *md_303 = buffer.data(md + 303 * ncomps + c);
        const auto *md_304 = buffer.data(md + 304 * ncomps + c);
        const auto *md_305 = buffer.data(md + 305 * ncomps + c);
        const auto *md_309 = buffer.data(md + 309 * ncomps + c);
        const auto *md_310 = buffer.data(md + 310 * ncomps + c);
        const auto *md_311 = buffer.data(md + 311 * ncomps + c);
        const auto *md_317 = buffer.data(md + 317 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, ld_174, ld_175, ld_176, \
                         ld_177, ld_178, md_174, md_175, md_176, md_177, \
                         md_178 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * ld_174[k]
                       + md_174[k];

            t_291[k] = ab_x[k] * ld_175[k]
                       + md_175[k];

            t_292[k] = ab_x[k] * ld_176[k]
                       + md_176[k];

            t_293[k] = ab_x[k] * ld_177[k]
                       + md_177[k];

            t_294[k] = ab_x[k] * ld_178[k]
                       + md_178[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, ab_y, ab_z, ld_177, ld_178, \
                         ld_179, md_179, md_225, md_226, md_227, \
                         md_233 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_x[k] * ld_179[k]
                       + md_179[k];

            t_296[k] = ab_y[k] * ld_177[k]
                       + md_225[k];

            t_297[k] = ab_y[k] * ld_178[k]
                       + md_226[k];

            t_298[k] = ab_y[k] * ld_179[k]
                       + md_227[k];

            t_299[k] = ab_z[k] * ld_179[k]
                       + md_233[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, ld_180, ld_181, ld_182, \
                         ld_183, ld_184, md_180, md_181, md_182, md_183, \
                         md_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * ld_180[k]
                       + md_180[k];

            t_301[k] = ab_x[k] * ld_181[k]
                       + md_181[k];

            t_302[k] = ab_x[k] * ld_182[k]
                       + md_182[k];

            t_303[k] = ab_x[k] * ld_183[k]
                       + md_183[k];

            t_304[k] = ab_x[k] * ld_184[k]
                       + md_184[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, ab_y, ab_z, ld_183, ld_184, \
                         ld_185, md_185, md_231, md_232, md_233, \
                         md_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * ld_185[k]
                       + md_185[k];

            t_306[k] = ab_y[k] * ld_183[k]
                       + md_231[k];

            t_307[k] = ab_y[k] * ld_184[k]
                       + md_232[k];

            t_308[k] = ab_y[k] * ld_185[k]
                       + md_233[k];

            t_309[k] = ab_z[k] * ld_185[k]
                       + md_239[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, ld_186, ld_187, ld_188, \
                         ld_189, ld_190, md_186, md_187, md_188, md_189, \
                         md_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_x[k] * ld_186[k]
                       + md_186[k];

            t_311[k] = ab_x[k] * ld_187[k]
                       + md_187[k];

            t_312[k] = ab_x[k] * ld_188[k]
                       + md_188[k];

            t_313[k] = ab_x[k] * ld_189[k]
                       + md_189[k];

            t_314[k] = ab_x[k] * ld_190[k]
                       + md_190[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, ab_y, ab_z, ld_189, ld_190, \
                         ld_191, md_191, md_237, md_238, md_239, \
                         md_245 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * ld_191[k]
                       + md_191[k];

            t_316[k] = ab_y[k] * ld_189[k]
                       + md_237[k];

            t_317[k] = ab_y[k] * ld_190[k]
                       + md_238[k];

            t_318[k] = ab_y[k] * ld_191[k]
                       + md_239[k];

            t_319[k] = ab_z[k] * ld_191[k]
                       + md_245[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, ld_192, ld_193, ld_194, \
                         ld_195, ld_196, md_192, md_193, md_194, md_195, \
                         md_196 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * ld_192[k]
                       + md_192[k];

            t_321[k] = ab_x[k] * ld_193[k]
                       + md_193[k];

            t_322[k] = ab_x[k] * ld_194[k]
                       + md_194[k];

            t_323[k] = ab_x[k] * ld_195[k]
                       + md_195[k];

            t_324[k] = ab_x[k] * ld_196[k]
                       + md_196[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, ab_y, ab_z, ld_195, ld_196, \
                         ld_197, md_197, md_243, md_244, md_245, \
                         md_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_x[k] * ld_197[k]
                       + md_197[k];

            t_326[k] = ab_y[k] * ld_195[k]
                       + md_243[k];

            t_327[k] = ab_y[k] * ld_196[k]
                       + md_244[k];

            t_328[k] = ab_y[k] * ld_197[k]
                       + md_245[k];

            t_329[k] = ab_z[k] * ld_197[k]
                       + md_251[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, ld_198, ld_199, ld_200, \
                         ld_201, ld_202, md_198, md_199, md_200, md_201, \
                         md_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * ld_198[k]
                       + md_198[k];

            t_331[k] = ab_x[k] * ld_199[k]
                       + md_199[k];

            t_332[k] = ab_x[k] * ld_200[k]
                       + md_200[k];

            t_333[k] = ab_x[k] * ld_201[k]
                       + md_201[k];

            t_334[k] = ab_x[k] * ld_202[k]
                       + md_202[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, ab_y, ab_z, ld_201, ld_202, \
                         ld_203, md_203, md_249, md_250, md_251, \
                         md_257 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * ld_203[k]
                       + md_203[k];

            t_336[k] = ab_y[k] * ld_201[k]
                       + md_249[k];

            t_337[k] = ab_y[k] * ld_202[k]
                       + md_250[k];

            t_338[k] = ab_y[k] * ld_203[k]
                       + md_251[k];

            t_339[k] = ab_z[k] * ld_203[k]
                       + md_257[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, ld_204, ld_205, ld_206, \
                         ld_207, ld_208, md_204, md_205, md_206, md_207, \
                         md_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_x[k] * ld_204[k]
                       + md_204[k];

            t_341[k] = ab_x[k] * ld_205[k]
                       + md_205[k];

            t_342[k] = ab_x[k] * ld_206[k]
                       + md_206[k];

            t_343[k] = ab_x[k] * ld_207[k]
                       + md_207[k];

            t_344[k] = ab_x[k] * ld_208[k]
                       + md_208[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, ab_y, ab_z, ld_207, ld_208, \
                         ld_209, md_209, md_255, md_256, md_257, \
                         md_263 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * ld_209[k]
                       + md_209[k];

            t_346[k] = ab_y[k] * ld_207[k]
                       + md_255[k];

            t_347[k] = ab_y[k] * ld_208[k]
                       + md_256[k];

            t_348[k] = ab_y[k] * ld_209[k]
                       + md_257[k];

            t_349[k] = ab_z[k] * ld_209[k]
                       + md_263[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, ld_210, ld_211, ld_212, \
                         ld_213, ld_214, md_210, md_211, md_212, md_213, \
                         md_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * ld_210[k]
                       + md_210[k];

            t_351[k] = ab_x[k] * ld_211[k]
                       + md_211[k];

            t_352[k] = ab_x[k] * ld_212[k]
                       + md_212[k];

            t_353[k] = ab_x[k] * ld_213[k]
                       + md_213[k];

            t_354[k] = ab_x[k] * ld_214[k]
                       + md_214[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, ab_y, ab_z, ld_213, ld_214, \
                         ld_215, md_215, md_261, md_262, md_263, \
                         md_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_x[k] * ld_215[k]
                       + md_215[k];

            t_356[k] = ab_y[k] * ld_213[k]
                       + md_261[k];

            t_357[k] = ab_y[k] * ld_214[k]
                       + md_262[k];

            t_358[k] = ab_y[k] * ld_215[k]
                       + md_263[k];

            t_359[k] = ab_z[k] * ld_215[k]
                       + md_269[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, ld_216, ld_217, ld_218, \
                         ld_219, ld_220, md_216, md_217, md_218, md_219, \
                         md_220 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * ld_216[k]
                       + md_216[k];

            t_361[k] = ab_x[k] * ld_217[k]
                       + md_217[k];

            t_362[k] = ab_x[k] * ld_218[k]
                       + md_218[k];

            t_363[k] = ab_x[k] * ld_219[k]
                       + md_219[k];

            t_364[k] = ab_x[k] * ld_220[k]
                       + md_220[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, ab_y, ab_z, ld_219, ld_220, \
                         ld_221, md_221, md_273, md_274, md_275, \
                         md_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * ld_221[k]
                       + md_221[k];

            t_366[k] = ab_y[k] * ld_219[k]
                       + md_273[k];

            t_367[k] = ab_y[k] * ld_220[k]
                       + md_274[k];

            t_368[k] = ab_y[k] * ld_221[k]
                       + md_275[k];

            t_369[k] = ab_z[k] * ld_221[k]
                       + md_281[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, ld_222, ld_223, ld_224, \
                         ld_225, ld_226, md_222, md_223, md_224, md_225, \
                         md_226 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * ld_222[k]
                       + md_222[k];

            t_371[k] = ab_x[k] * ld_223[k]
                       + md_223[k];

            t_372[k] = ab_x[k] * ld_224[k]
                       + md_224[k];

            t_373[k] = ab_x[k] * ld_225[k]
                       + md_225[k];

            t_374[k] = ab_x[k] * ld_226[k]
                       + md_226[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, ab_y, ab_z, ld_225, ld_226, \
                         ld_227, md_227, md_279, md_280, md_281, \
                         md_287 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * ld_227[k]
                       + md_227[k];

            t_376[k] = ab_y[k] * ld_225[k]
                       + md_279[k];

            t_377[k] = ab_y[k] * ld_226[k]
                       + md_280[k];

            t_378[k] = ab_y[k] * ld_227[k]
                       + md_281[k];

            t_379[k] = ab_z[k] * ld_227[k]
                       + md_287[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, ld_228, ld_229, ld_230, \
                         ld_231, ld_232, md_228, md_229, md_230, md_231, \
                         md_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * ld_228[k]
                       + md_228[k];

            t_381[k] = ab_x[k] * ld_229[k]
                       + md_229[k];

            t_382[k] = ab_x[k] * ld_230[k]
                       + md_230[k];

            t_383[k] = ab_x[k] * ld_231[k]
                       + md_231[k];

            t_384[k] = ab_x[k] * ld_232[k]
                       + md_232[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, ab_y, ab_z, ld_231, ld_232, \
                         ld_233, md_233, md_285, md_286, md_287, \
                         md_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_x[k] * ld_233[k]
                       + md_233[k];

            t_386[k] = ab_y[k] * ld_231[k]
                       + md_285[k];

            t_387[k] = ab_y[k] * ld_232[k]
                       + md_286[k];

            t_388[k] = ab_y[k] * ld_233[k]
                       + md_287[k];

            t_389[k] = ab_z[k] * ld_233[k]
                       + md_293[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, ld_234, ld_235, ld_236, \
                         ld_237, ld_238, md_234, md_235, md_236, md_237, \
                         md_238 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * ld_234[k]
                       + md_234[k];

            t_391[k] = ab_x[k] * ld_235[k]
                       + md_235[k];

            t_392[k] = ab_x[k] * ld_236[k]
                       + md_236[k];

            t_393[k] = ab_x[k] * ld_237[k]
                       + md_237[k];

            t_394[k] = ab_x[k] * ld_238[k]
                       + md_238[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, ab_y, ab_z, ld_237, ld_238, \
                         ld_239, md_239, md_291, md_292, md_293, \
                         md_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * ld_239[k]
                       + md_239[k];

            t_396[k] = ab_y[k] * ld_237[k]
                       + md_291[k];

            t_397[k] = ab_y[k] * ld_238[k]
                       + md_292[k];

            t_398[k] = ab_y[k] * ld_239[k]
                       + md_293[k];

            t_399[k] = ab_z[k] * ld_239[k]
                       + md_299[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, ld_240, ld_241, ld_242, \
                         ld_243, ld_244, md_240, md_241, md_242, md_243, \
                         md_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_x[k] * ld_240[k]
                       + md_240[k];

            t_401[k] = ab_x[k] * ld_241[k]
                       + md_241[k];

            t_402[k] = ab_x[k] * ld_242[k]
                       + md_242[k];

            t_403[k] = ab_x[k] * ld_243[k]
                       + md_243[k];

            t_404[k] = ab_x[k] * ld_244[k]
                       + md_244[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, ab_y, ab_z, ld_243, ld_244, \
                         ld_245, md_245, md_297, md_298, md_299, \
                         md_305 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * ld_245[k]
                       + md_245[k];

            t_406[k] = ab_y[k] * ld_243[k]
                       + md_297[k];

            t_407[k] = ab_y[k] * ld_244[k]
                       + md_298[k];

            t_408[k] = ab_y[k] * ld_245[k]
                       + md_299[k];

            t_409[k] = ab_z[k] * ld_245[k]
                       + md_305[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, ld_246, ld_247, ld_248, \
                         ld_249, ld_250, md_246, md_247, md_248, md_249, \
                         md_250 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * ld_246[k]
                       + md_246[k];

            t_411[k] = ab_x[k] * ld_247[k]
                       + md_247[k];

            t_412[k] = ab_x[k] * ld_248[k]
                       + md_248[k];

            t_413[k] = ab_x[k] * ld_249[k]
                       + md_249[k];

            t_414[k] = ab_x[k] * ld_250[k]
                       + md_250[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, ab_y, ab_z, ld_249, ld_250, \
                         ld_251, md_251, md_303, md_304, md_305, \
                         md_311 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_x[k] * ld_251[k]
                       + md_251[k];

            t_416[k] = ab_y[k] * ld_249[k]
                       + md_303[k];

            t_417[k] = ab_y[k] * ld_250[k]
                       + md_304[k];

            t_418[k] = ab_y[k] * ld_251[k]
                       + md_305[k];

            t_419[k] = ab_z[k] * ld_251[k]
                       + md_311[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, ld_252, ld_253, ld_254, \
                         ld_255, ld_256, md_252, md_253, md_254, md_255, \
                         md_256 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * ld_252[k]
                       + md_252[k];

            t_421[k] = ab_x[k] * ld_253[k]
                       + md_253[k];

            t_422[k] = ab_x[k] * ld_254[k]
                       + md_254[k];

            t_423[k] = ab_x[k] * ld_255[k]
                       + md_255[k];

            t_424[k] = ab_x[k] * ld_256[k]
                       + md_256[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, ab_y, ab_z, ld_255, ld_256, \
                         ld_257, md_257, md_309, md_310, md_311, \
                         md_317 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * ld_257[k]
                       + md_257[k];

            t_426[k] = ab_y[k] * ld_255[k]
                       + md_309[k];

            t_427[k] = ab_y[k] * ld_256[k]
                       + md_310[k];

            t_428[k] = ab_y[k] * ld_257[k]
                       + md_311[k];

            t_429[k] = ab_z[k] * ld_257[k]
                       + md_317[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, ld_258, ld_259, ld_260, \
                         ld_261, ld_262, md_258, md_259, md_260, md_261, \
                         md_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * ld_258[k]
                       + md_258[k];

            t_431[k] = ab_x[k] * ld_259[k]
                       + md_259[k];

            t_432[k] = ab_x[k] * ld_260[k]
                       + md_260[k];

            t_433[k] = ab_x[k] * ld_261[k]
                       + md_261[k];

            t_434[k] = ab_x[k] * ld_262[k]
                       + md_262[k];
        }
    }
}

static auto
compute_hrr_lf_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ld, const size_t md, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ld_261 = buffer.data(ld + 261 * ncomps + c);
        const auto *ld_262 = buffer.data(ld + 262 * ncomps + c);
        const auto *ld_263 = buffer.data(ld + 263 * ncomps + c);
        const auto *ld_264 = buffer.data(ld + 264 * ncomps + c);
        const auto *ld_265 = buffer.data(ld + 265 * ncomps + c);
        const auto *ld_266 = buffer.data(ld + 266 * ncomps + c);
        const auto *ld_267 = buffer.data(ld + 267 * ncomps + c);
        const auto *ld_268 = buffer.data(ld + 268 * ncomps + c);
        const auto *ld_269 = buffer.data(ld + 269 * ncomps + c);

        const auto *md_263 = buffer.data(md + 263 * ncomps + c);
        const auto *md_264 = buffer.data(md + 264 * ncomps + c);
        const auto *md_265 = buffer.data(md + 265 * ncomps + c);
        const auto *md_266 = buffer.data(md + 266 * ncomps + c);
        const auto *md_267 = buffer.data(md + 267 * ncomps + c);
        const auto *md_268 = buffer.data(md + 268 * ncomps + c);
        const auto *md_269 = buffer.data(md + 269 * ncomps + c);
        const auto *md_315 = buffer.data(md + 315 * ncomps + c);
        const auto *md_316 = buffer.data(md + 316 * ncomps + c);
        const auto *md_317 = buffer.data(md + 317 * ncomps + c);
        const auto *md_321 = buffer.data(md + 321 * ncomps + c);
        const auto *md_322 = buffer.data(md + 322 * ncomps + c);
        const auto *md_323 = buffer.data(md + 323 * ncomps + c);
        const auto *md_329 = buffer.data(md + 329 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, ab_y, ab_z, ld_261, ld_262, \
                         ld_263, md_263, md_315, md_316, md_317, \
                         md_323 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * ld_263[k]
                       + md_263[k];

            t_436[k] = ab_y[k] * ld_261[k]
                       + md_315[k];

            t_437[k] = ab_y[k] * ld_262[k]
                       + md_316[k];

            t_438[k] = ab_y[k] * ld_263[k]
                       + md_317[k];

            t_439[k] = ab_z[k] * ld_263[k]
                       + md_323[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, ld_264, ld_265, ld_266, \
                         ld_267, ld_268, md_264, md_265, md_266, md_267, \
                         md_268 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * ld_264[k]
                       + md_264[k];

            t_441[k] = ab_x[k] * ld_265[k]
                       + md_265[k];

            t_442[k] = ab_x[k] * ld_266[k]
                       + md_266[k];

            t_443[k] = ab_x[k] * ld_267[k]
                       + md_267[k];

            t_444[k] = ab_x[k] * ld_268[k]
                       + md_268[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_x, ab_y, ab_z, ld_267, ld_268, \
                         ld_269, md_269, md_321, md_322, md_323, \
                         md_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = ab_x[k] * ld_269[k]
                       + md_269[k];

            t_446[k] = ab_y[k] * ld_267[k]
                       + md_321[k];

            t_447[k] = ab_y[k] * ld_268[k]
                       + md_322[k];

            t_448[k] = ab_y[k] * ld_269[k]
                       + md_323[k];

            t_449[k] = ab_z[k] * ld_269[k]
                       + md_329[k];
        }
    }
}

auto
compute_hrr_lf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ld, const size_t md, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_lf_piece0(buffer, coordinates, target, ld, md, ncomps, nmax);

    compute_hrr_lf_piece1(buffer, coordinates, target, ld, md, ncomps, nmax);

    compute_hrr_lf_piece2(buffer, coordinates, target, ld, md, ncomps, nmax);

    compute_hrr_lf_piece3(buffer, coordinates, target, ld, md, ncomps, nmax);
}

}  // namespace simdtrf
