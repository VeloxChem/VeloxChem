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


#include "SimdTransferGeom100YGG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_100y_gg_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t gf_1,
                                             const size_t gf_0, const size_t hf_1,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gf_1_0 = buffer.data(gf_1 + 0 * ncomps + c);
        const auto *gf_1_1 = buffer.data(gf_1 + 1 * ncomps + c);
        const auto *gf_1_2 = buffer.data(gf_1 + 2 * ncomps + c);
        const auto *gf_1_3 = buffer.data(gf_1 + 3 * ncomps + c);
        const auto *gf_1_4 = buffer.data(gf_1 + 4 * ncomps + c);
        const auto *gf_1_5 = buffer.data(gf_1 + 5 * ncomps + c);
        const auto *gf_1_6 = buffer.data(gf_1 + 6 * ncomps + c);
        const auto *gf_1_7 = buffer.data(gf_1 + 7 * ncomps + c);
        const auto *gf_1_8 = buffer.data(gf_1 + 8 * ncomps + c);
        const auto *gf_1_9 = buffer.data(gf_1 + 9 * ncomps + c);
        const auto *gf_1_10 = buffer.data(gf_1 + 10 * ncomps + c);
        const auto *gf_1_11 = buffer.data(gf_1 + 11 * ncomps + c);
        const auto *gf_1_12 = buffer.data(gf_1 + 12 * ncomps + c);
        const auto *gf_1_13 = buffer.data(gf_1 + 13 * ncomps + c);
        const auto *gf_1_14 = buffer.data(gf_1 + 14 * ncomps + c);
        const auto *gf_1_15 = buffer.data(gf_1 + 15 * ncomps + c);
        const auto *gf_1_16 = buffer.data(gf_1 + 16 * ncomps + c);
        const auto *gf_1_17 = buffer.data(gf_1 + 17 * ncomps + c);
        const auto *gf_1_18 = buffer.data(gf_1 + 18 * ncomps + c);
        const auto *gf_1_19 = buffer.data(gf_1 + 19 * ncomps + c);
        const auto *gf_1_20 = buffer.data(gf_1 + 20 * ncomps + c);
        const auto *gf_1_21 = buffer.data(gf_1 + 21 * ncomps + c);
        const auto *gf_1_22 = buffer.data(gf_1 + 22 * ncomps + c);
        const auto *gf_1_23 = buffer.data(gf_1 + 23 * ncomps + c);
        const auto *gf_1_24 = buffer.data(gf_1 + 24 * ncomps + c);
        const auto *gf_1_25 = buffer.data(gf_1 + 25 * ncomps + c);
        const auto *gf_1_26 = buffer.data(gf_1 + 26 * ncomps + c);
        const auto *gf_1_27 = buffer.data(gf_1 + 27 * ncomps + c);
        const auto *gf_1_28 = buffer.data(gf_1 + 28 * ncomps + c);
        const auto *gf_1_29 = buffer.data(gf_1 + 29 * ncomps + c);
        const auto *gf_1_30 = buffer.data(gf_1 + 30 * ncomps + c);
        const auto *gf_1_31 = buffer.data(gf_1 + 31 * ncomps + c);
        const auto *gf_1_32 = buffer.data(gf_1 + 32 * ncomps + c);
        const auto *gf_1_33 = buffer.data(gf_1 + 33 * ncomps + c);
        const auto *gf_1_34 = buffer.data(gf_1 + 34 * ncomps + c);
        const auto *gf_1_35 = buffer.data(gf_1 + 35 * ncomps + c);
        const auto *gf_1_36 = buffer.data(gf_1 + 36 * ncomps + c);
        const auto *gf_1_37 = buffer.data(gf_1 + 37 * ncomps + c);
        const auto *gf_1_38 = buffer.data(gf_1 + 38 * ncomps + c);
        const auto *gf_1_39 = buffer.data(gf_1 + 39 * ncomps + c);
        const auto *gf_1_40 = buffer.data(gf_1 + 40 * ncomps + c);
        const auto *gf_1_41 = buffer.data(gf_1 + 41 * ncomps + c);
        const auto *gf_1_42 = buffer.data(gf_1 + 42 * ncomps + c);
        const auto *gf_1_43 = buffer.data(gf_1 + 43 * ncomps + c);
        const auto *gf_1_44 = buffer.data(gf_1 + 44 * ncomps + c);
        const auto *gf_1_45 = buffer.data(gf_1 + 45 * ncomps + c);
        const auto *gf_1_46 = buffer.data(gf_1 + 46 * ncomps + c);
        const auto *gf_1_47 = buffer.data(gf_1 + 47 * ncomps + c);
        const auto *gf_1_48 = buffer.data(gf_1 + 48 * ncomps + c);
        const auto *gf_1_49 = buffer.data(gf_1 + 49 * ncomps + c);
        const auto *gf_1_50 = buffer.data(gf_1 + 50 * ncomps + c);
        const auto *gf_1_51 = buffer.data(gf_1 + 51 * ncomps + c);
        const auto *gf_1_52 = buffer.data(gf_1 + 52 * ncomps + c);
        const auto *gf_1_53 = buffer.data(gf_1 + 53 * ncomps + c);
        const auto *gf_1_54 = buffer.data(gf_1 + 54 * ncomps + c);
        const auto *gf_1_55 = buffer.data(gf_1 + 55 * ncomps + c);
        const auto *gf_1_56 = buffer.data(gf_1 + 56 * ncomps + c);
        const auto *gf_1_57 = buffer.data(gf_1 + 57 * ncomps + c);
        const auto *gf_1_58 = buffer.data(gf_1 + 58 * ncomps + c);
        const auto *gf_1_59 = buffer.data(gf_1 + 59 * ncomps + c);
        const auto *gf_1_60 = buffer.data(gf_1 + 60 * ncomps + c);
        const auto *gf_1_61 = buffer.data(gf_1 + 61 * ncomps + c);
        const auto *gf_1_62 = buffer.data(gf_1 + 62 * ncomps + c);
        const auto *gf_1_63 = buffer.data(gf_1 + 63 * ncomps + c);
        const auto *gf_1_64 = buffer.data(gf_1 + 64 * ncomps + c);
        const auto *gf_1_65 = buffer.data(gf_1 + 65 * ncomps + c);
        const auto *gf_1_66 = buffer.data(gf_1 + 66 * ncomps + c);
        const auto *gf_1_67 = buffer.data(gf_1 + 67 * ncomps + c);
        const auto *gf_1_68 = buffer.data(gf_1 + 68 * ncomps + c);
        const auto *gf_1_69 = buffer.data(gf_1 + 69 * ncomps + c);
        const auto *gf_1_70 = buffer.data(gf_1 + 70 * ncomps + c);
        const auto *gf_1_71 = buffer.data(gf_1 + 71 * ncomps + c);
        const auto *gf_1_72 = buffer.data(gf_1 + 72 * ncomps + c);
        const auto *gf_1_73 = buffer.data(gf_1 + 73 * ncomps + c);
        const auto *gf_1_74 = buffer.data(gf_1 + 74 * ncomps + c);
        const auto *gf_1_75 = buffer.data(gf_1 + 75 * ncomps + c);
        const auto *gf_1_76 = buffer.data(gf_1 + 76 * ncomps + c);
        const auto *gf_1_77 = buffer.data(gf_1 + 77 * ncomps + c);
        const auto *gf_1_78 = buffer.data(gf_1 + 78 * ncomps + c);
        const auto *gf_1_79 = buffer.data(gf_1 + 79 * ncomps + c);
        const auto *gf_1_80 = buffer.data(gf_1 + 80 * ncomps + c);
        const auto *gf_1_81 = buffer.data(gf_1 + 81 * ncomps + c);
        const auto *gf_1_82 = buffer.data(gf_1 + 82 * ncomps + c);
        const auto *gf_1_83 = buffer.data(gf_1 + 83 * ncomps + c);
        const auto *gf_1_84 = buffer.data(gf_1 + 84 * ncomps + c);
        const auto *gf_1_85 = buffer.data(gf_1 + 85 * ncomps + c);
        const auto *gf_1_86 = buffer.data(gf_1 + 86 * ncomps + c);
        const auto *gf_1_87 = buffer.data(gf_1 + 87 * ncomps + c);
        const auto *gf_1_88 = buffer.data(gf_1 + 88 * ncomps + c);
        const auto *gf_1_89 = buffer.data(gf_1 + 89 * ncomps + c);

        const auto *gf_0_6 = buffer.data(gf_0 + 6 * ncomps + c);
        const auto *gf_0_7 = buffer.data(gf_0 + 7 * ncomps + c);
        const auto *gf_0_8 = buffer.data(gf_0 + 8 * ncomps + c);
        const auto *gf_0_9 = buffer.data(gf_0 + 9 * ncomps + c);
        const auto *gf_0_16 = buffer.data(gf_0 + 16 * ncomps + c);
        const auto *gf_0_17 = buffer.data(gf_0 + 17 * ncomps + c);
        const auto *gf_0_18 = buffer.data(gf_0 + 18 * ncomps + c);
        const auto *gf_0_19 = buffer.data(gf_0 + 19 * ncomps + c);
        const auto *gf_0_26 = buffer.data(gf_0 + 26 * ncomps + c);
        const auto *gf_0_27 = buffer.data(gf_0 + 27 * ncomps + c);
        const auto *gf_0_28 = buffer.data(gf_0 + 28 * ncomps + c);
        const auto *gf_0_29 = buffer.data(gf_0 + 29 * ncomps + c);
        const auto *gf_0_36 = buffer.data(gf_0 + 36 * ncomps + c);
        const auto *gf_0_37 = buffer.data(gf_0 + 37 * ncomps + c);
        const auto *gf_0_38 = buffer.data(gf_0 + 38 * ncomps + c);
        const auto *gf_0_39 = buffer.data(gf_0 + 39 * ncomps + c);
        const auto *gf_0_46 = buffer.data(gf_0 + 46 * ncomps + c);
        const auto *gf_0_47 = buffer.data(gf_0 + 47 * ncomps + c);
        const auto *gf_0_48 = buffer.data(gf_0 + 48 * ncomps + c);
        const auto *gf_0_49 = buffer.data(gf_0 + 49 * ncomps + c);
        const auto *gf_0_56 = buffer.data(gf_0 + 56 * ncomps + c);
        const auto *gf_0_57 = buffer.data(gf_0 + 57 * ncomps + c);
        const auto *gf_0_58 = buffer.data(gf_0 + 58 * ncomps + c);
        const auto *gf_0_59 = buffer.data(gf_0 + 59 * ncomps + c);
        const auto *gf_0_66 = buffer.data(gf_0 + 66 * ncomps + c);
        const auto *gf_0_67 = buffer.data(gf_0 + 67 * ncomps + c);
        const auto *gf_0_68 = buffer.data(gf_0 + 68 * ncomps + c);
        const auto *gf_0_69 = buffer.data(gf_0 + 69 * ncomps + c);
        const auto *gf_0_76 = buffer.data(gf_0 + 76 * ncomps + c);
        const auto *gf_0_77 = buffer.data(gf_0 + 77 * ncomps + c);
        const auto *gf_0_78 = buffer.data(gf_0 + 78 * ncomps + c);
        const auto *gf_0_79 = buffer.data(gf_0 + 79 * ncomps + c);
        const auto *gf_0_86 = buffer.data(gf_0 + 86 * ncomps + c);
        const auto *gf_0_87 = buffer.data(gf_0 + 87 * ncomps + c);
        const auto *gf_0_88 = buffer.data(gf_0 + 88 * ncomps + c);

        const auto *hf_1_0 = buffer.data(hf_1 + 0 * ncomps + c);
        const auto *hf_1_1 = buffer.data(hf_1 + 1 * ncomps + c);
        const auto *hf_1_2 = buffer.data(hf_1 + 2 * ncomps + c);
        const auto *hf_1_3 = buffer.data(hf_1 + 3 * ncomps + c);
        const auto *hf_1_4 = buffer.data(hf_1 + 4 * ncomps + c);
        const auto *hf_1_5 = buffer.data(hf_1 + 5 * ncomps + c);
        const auto *hf_1_6 = buffer.data(hf_1 + 6 * ncomps + c);
        const auto *hf_1_7 = buffer.data(hf_1 + 7 * ncomps + c);
        const auto *hf_1_8 = buffer.data(hf_1 + 8 * ncomps + c);
        const auto *hf_1_9 = buffer.data(hf_1 + 9 * ncomps + c);
        const auto *hf_1_10 = buffer.data(hf_1 + 10 * ncomps + c);
        const auto *hf_1_11 = buffer.data(hf_1 + 11 * ncomps + c);
        const auto *hf_1_12 = buffer.data(hf_1 + 12 * ncomps + c);
        const auto *hf_1_13 = buffer.data(hf_1 + 13 * ncomps + c);
        const auto *hf_1_14 = buffer.data(hf_1 + 14 * ncomps + c);
        const auto *hf_1_15 = buffer.data(hf_1 + 15 * ncomps + c);
        const auto *hf_1_16 = buffer.data(hf_1 + 16 * ncomps + c);
        const auto *hf_1_17 = buffer.data(hf_1 + 17 * ncomps + c);
        const auto *hf_1_18 = buffer.data(hf_1 + 18 * ncomps + c);
        const auto *hf_1_19 = buffer.data(hf_1 + 19 * ncomps + c);
        const auto *hf_1_20 = buffer.data(hf_1 + 20 * ncomps + c);
        const auto *hf_1_21 = buffer.data(hf_1 + 21 * ncomps + c);
        const auto *hf_1_22 = buffer.data(hf_1 + 22 * ncomps + c);
        const auto *hf_1_23 = buffer.data(hf_1 + 23 * ncomps + c);
        const auto *hf_1_24 = buffer.data(hf_1 + 24 * ncomps + c);
        const auto *hf_1_25 = buffer.data(hf_1 + 25 * ncomps + c);
        const auto *hf_1_26 = buffer.data(hf_1 + 26 * ncomps + c);
        const auto *hf_1_27 = buffer.data(hf_1 + 27 * ncomps + c);
        const auto *hf_1_28 = buffer.data(hf_1 + 28 * ncomps + c);
        const auto *hf_1_29 = buffer.data(hf_1 + 29 * ncomps + c);
        const auto *hf_1_30 = buffer.data(hf_1 + 30 * ncomps + c);
        const auto *hf_1_31 = buffer.data(hf_1 + 31 * ncomps + c);
        const auto *hf_1_32 = buffer.data(hf_1 + 32 * ncomps + c);
        const auto *hf_1_33 = buffer.data(hf_1 + 33 * ncomps + c);
        const auto *hf_1_34 = buffer.data(hf_1 + 34 * ncomps + c);
        const auto *hf_1_35 = buffer.data(hf_1 + 35 * ncomps + c);
        const auto *hf_1_36 = buffer.data(hf_1 + 36 * ncomps + c);
        const auto *hf_1_37 = buffer.data(hf_1 + 37 * ncomps + c);
        const auto *hf_1_38 = buffer.data(hf_1 + 38 * ncomps + c);
        const auto *hf_1_39 = buffer.data(hf_1 + 39 * ncomps + c);
        const auto *hf_1_40 = buffer.data(hf_1 + 40 * ncomps + c);
        const auto *hf_1_41 = buffer.data(hf_1 + 41 * ncomps + c);
        const auto *hf_1_42 = buffer.data(hf_1 + 42 * ncomps + c);
        const auto *hf_1_43 = buffer.data(hf_1 + 43 * ncomps + c);
        const auto *hf_1_44 = buffer.data(hf_1 + 44 * ncomps + c);
        const auto *hf_1_45 = buffer.data(hf_1 + 45 * ncomps + c);
        const auto *hf_1_46 = buffer.data(hf_1 + 46 * ncomps + c);
        const auto *hf_1_47 = buffer.data(hf_1 + 47 * ncomps + c);
        const auto *hf_1_48 = buffer.data(hf_1 + 48 * ncomps + c);
        const auto *hf_1_49 = buffer.data(hf_1 + 49 * ncomps + c);
        const auto *hf_1_50 = buffer.data(hf_1 + 50 * ncomps + c);
        const auto *hf_1_51 = buffer.data(hf_1 + 51 * ncomps + c);
        const auto *hf_1_52 = buffer.data(hf_1 + 52 * ncomps + c);
        const auto *hf_1_53 = buffer.data(hf_1 + 53 * ncomps + c);
        const auto *hf_1_54 = buffer.data(hf_1 + 54 * ncomps + c);
        const auto *hf_1_55 = buffer.data(hf_1 + 55 * ncomps + c);
        const auto *hf_1_56 = buffer.data(hf_1 + 56 * ncomps + c);
        const auto *hf_1_57 = buffer.data(hf_1 + 57 * ncomps + c);
        const auto *hf_1_58 = buffer.data(hf_1 + 58 * ncomps + c);
        const auto *hf_1_59 = buffer.data(hf_1 + 59 * ncomps + c);
        const auto *hf_1_60 = buffer.data(hf_1 + 60 * ncomps + c);
        const auto *hf_1_61 = buffer.data(hf_1 + 61 * ncomps + c);
        const auto *hf_1_62 = buffer.data(hf_1 + 62 * ncomps + c);
        const auto *hf_1_63 = buffer.data(hf_1 + 63 * ncomps + c);
        const auto *hf_1_64 = buffer.data(hf_1 + 64 * ncomps + c);
        const auto *hf_1_65 = buffer.data(hf_1 + 65 * ncomps + c);
        const auto *hf_1_66 = buffer.data(hf_1 + 66 * ncomps + c);
        const auto *hf_1_67 = buffer.data(hf_1 + 67 * ncomps + c);
        const auto *hf_1_68 = buffer.data(hf_1 + 68 * ncomps + c);
        const auto *hf_1_69 = buffer.data(hf_1 + 69 * ncomps + c);
        const auto *hf_1_70 = buffer.data(hf_1 + 70 * ncomps + c);
        const auto *hf_1_71 = buffer.data(hf_1 + 71 * ncomps + c);
        const auto *hf_1_72 = buffer.data(hf_1 + 72 * ncomps + c);
        const auto *hf_1_73 = buffer.data(hf_1 + 73 * ncomps + c);
        const auto *hf_1_74 = buffer.data(hf_1 + 74 * ncomps + c);
        const auto *hf_1_75 = buffer.data(hf_1 + 75 * ncomps + c);
        const auto *hf_1_76 = buffer.data(hf_1 + 76 * ncomps + c);
        const auto *hf_1_77 = buffer.data(hf_1 + 77 * ncomps + c);
        const auto *hf_1_78 = buffer.data(hf_1 + 78 * ncomps + c);
        const auto *hf_1_79 = buffer.data(hf_1 + 79 * ncomps + c);
        const auto *hf_1_80 = buffer.data(hf_1 + 80 * ncomps + c);
        const auto *hf_1_81 = buffer.data(hf_1 + 81 * ncomps + c);
        const auto *hf_1_82 = buffer.data(hf_1 + 82 * ncomps + c);
        const auto *hf_1_83 = buffer.data(hf_1 + 83 * ncomps + c);
        const auto *hf_1_84 = buffer.data(hf_1 + 84 * ncomps + c);
        const auto *hf_1_85 = buffer.data(hf_1 + 85 * ncomps + c);
        const auto *hf_1_86 = buffer.data(hf_1 + 86 * ncomps + c);
        const auto *hf_1_87 = buffer.data(hf_1 + 87 * ncomps + c);
        const auto *hf_1_88 = buffer.data(hf_1 + 88 * ncomps + c);
        const auto *hf_1_89 = buffer.data(hf_1 + 89 * ncomps + c);
        const auto *hf_1_99 = buffer.data(hf_1 + 99 * ncomps + c);
        const auto *hf_1_106 = buffer.data(hf_1 + 106 * ncomps + c);
        const auto *hf_1_107 = buffer.data(hf_1 + 107 * ncomps + c);
        const auto *hf_1_108 = buffer.data(hf_1 + 108 * ncomps + c);
        const auto *hf_1_109 = buffer.data(hf_1 + 109 * ncomps + c);
        const auto *hf_1_116 = buffer.data(hf_1 + 116 * ncomps + c);
        const auto *hf_1_117 = buffer.data(hf_1 + 117 * ncomps + c);
        const auto *hf_1_118 = buffer.data(hf_1 + 118 * ncomps + c);
        const auto *hf_1_119 = buffer.data(hf_1 + 119 * ncomps + c);
        const auto *hf_1_126 = buffer.data(hf_1 + 126 * ncomps + c);
        const auto *hf_1_127 = buffer.data(hf_1 + 127 * ncomps + c);
        const auto *hf_1_128 = buffer.data(hf_1 + 128 * ncomps + c);
        const auto *hf_1_129 = buffer.data(hf_1 + 129 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gf_1_0, gf_1_1, gf_1_2, gf_1_3, \
                         gf_1_4, hf_1_0, hf_1_1, hf_1_2, hf_1_3, \
                         hf_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gf_1_0[k]
                     + hf_1_0[k];

            t_1[k] = ab_x[k] * gf_1_1[k]
                     + hf_1_1[k];

            t_2[k] = ab_x[k] * gf_1_2[k]
                     + hf_1_2[k];

            t_3[k] = ab_x[k] * gf_1_3[k]
                     + hf_1_3[k];

            t_4[k] = ab_x[k] * gf_1_4[k]
                     + hf_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gf_1_5, gf_1_6, gf_1_7, gf_1_8, \
                         gf_1_9, hf_1_5, hf_1_6, hf_1_7, hf_1_8, \
                         hf_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * gf_1_5[k]
                     + hf_1_5[k];

            t_6[k] = ab_x[k] * gf_1_6[k]
                     + hf_1_6[k];

            t_7[k] = ab_x[k] * gf_1_7[k]
                     + hf_1_7[k];

            t_8[k] = ab_x[k] * gf_1_8[k]
                     + hf_1_8[k];

            t_9[k] = ab_x[k] * gf_1_9[k]
                     + hf_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, ab_y, gf_1_6, gf_1_7, gf_1_8, gf_0_6, gf_0_7, \
                         gf_0_8, hf_1_16, hf_1_17, hf_1_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * gf_1_6[k]
                      + gf_0_6[k]
                      + hf_1_16[k];

            t_11[k] = ab_y[k] * gf_1_7[k]
                      + gf_0_7[k]
                      + hf_1_17[k];

            t_12[k] = ab_y[k] * gf_1_8[k]
                      + gf_0_8[k]
                      + hf_1_18[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, gf_1_9, gf_1_10, gf_1_11, \
                         gf_0_9, hf_1_10, hf_1_11, hf_1_19, hf_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_y[k] * gf_1_9[k]
                      + gf_0_9[k]
                      + hf_1_19[k];

            t_14[k] = ab_z[k] * gf_1_9[k]
                      + hf_1_29[k];

            t_15[k] = ab_x[k] * gf_1_10[k]
                      + hf_1_10[k];

            t_16[k] = ab_x[k] * gf_1_11[k]
                      + hf_1_11[k];
        }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, gf_1_12, gf_1_13, gf_1_14, \
                         gf_1_15, gf_1_16, hf_1_12, hf_1_13, hf_1_14, hf_1_15, \
                         hf_1_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_17[k] = ab_x[k] * gf_1_12[k]
                      + hf_1_12[k];

            t_18[k] = ab_x[k] * gf_1_13[k]
                      + hf_1_13[k];

            t_19[k] = ab_x[k] * gf_1_14[k]
                      + hf_1_14[k];

            t_20[k] = ab_x[k] * gf_1_15[k]
                      + hf_1_15[k];

            t_21[k] = ab_x[k] * gf_1_16[k]
                      + hf_1_16[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, ab_x, ab_y, gf_1_16, gf_1_17, gf_1_18, \
                         gf_1_19, gf_0_16, hf_1_17, hf_1_18, hf_1_19, \
                         hf_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_x[k] * gf_1_17[k]
                      + hf_1_17[k];

            t_23[k] = ab_x[k] * gf_1_18[k]
                      + hf_1_18[k];

            t_24[k] = ab_x[k] * gf_1_19[k]
                      + hf_1_19[k];

            t_25[k] = ab_y[k] * gf_1_16[k]
                      + gf_0_16[k]
                      + hf_1_36[k];
        }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, ab_y, ab_z, gf_1_17, gf_1_18, gf_1_19, \
                         gf_0_17, gf_0_18, gf_0_19, hf_1_37, hf_1_38, hf_1_39, \
                         hf_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_26[k] = ab_y[k] * gf_1_17[k]
                      + gf_0_17[k]
                      + hf_1_37[k];

            t_27[k] = ab_y[k] * gf_1_18[k]
                      + gf_0_18[k]
                      + hf_1_38[k];

            t_28[k] = ab_y[k] * gf_1_19[k]
                      + gf_0_19[k]
                      + hf_1_39[k];

            t_29[k] = ab_z[k] * gf_1_19[k]
                      + hf_1_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gf_1_20, gf_1_21, gf_1_22, \
                         gf_1_23, gf_1_24, hf_1_20, hf_1_21, hf_1_22, hf_1_23, \
                         hf_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * gf_1_20[k]
                      + hf_1_20[k];

            t_31[k] = ab_x[k] * gf_1_21[k]
                      + hf_1_21[k];

            t_32[k] = ab_x[k] * gf_1_22[k]
                      + hf_1_22[k];

            t_33[k] = ab_x[k] * gf_1_23[k]
                      + hf_1_23[k];

            t_34[k] = ab_x[k] * gf_1_24[k]
                      + hf_1_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gf_1_25, gf_1_26, gf_1_27, \
                         gf_1_28, gf_1_29, hf_1_25, hf_1_26, hf_1_27, hf_1_28, \
                         hf_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * gf_1_25[k]
                      + hf_1_25[k];

            t_36[k] = ab_x[k] * gf_1_26[k]
                      + hf_1_26[k];

            t_37[k] = ab_x[k] * gf_1_27[k]
                      + hf_1_27[k];

            t_38[k] = ab_x[k] * gf_1_28[k]
                      + hf_1_28[k];

            t_39[k] = ab_x[k] * gf_1_29[k]
                      + hf_1_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, ab_y, gf_1_26, gf_1_27, gf_1_28, gf_0_26, gf_0_27, \
                         gf_0_28, hf_1_46, hf_1_47, hf_1_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * gf_1_26[k]
                      + gf_0_26[k]
                      + hf_1_46[k];

            t_41[k] = ab_y[k] * gf_1_27[k]
                      + gf_0_27[k]
                      + hf_1_47[k];

            t_42[k] = ab_y[k] * gf_1_28[k]
                      + gf_0_28[k]
                      + hf_1_48[k];
        }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, ab_x, ab_y, ab_z, gf_1_29, gf_1_30, gf_1_31, \
                         gf_0_29, hf_1_30, hf_1_31, hf_1_49, hf_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_43[k] = ab_y[k] * gf_1_29[k]
                      + gf_0_29[k]
                      + hf_1_49[k];

            t_44[k] = ab_z[k] * gf_1_29[k]
                      + hf_1_59[k];

            t_45[k] = ab_x[k] * gf_1_30[k]
                      + hf_1_30[k];

            t_46[k] = ab_x[k] * gf_1_31[k]
                      + hf_1_31[k];
        }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, ab_x, gf_1_32, gf_1_33, gf_1_34, \
                         gf_1_35, gf_1_36, hf_1_32, hf_1_33, hf_1_34, hf_1_35, \
                         hf_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_47[k] = ab_x[k] * gf_1_32[k]
                      + hf_1_32[k];

            t_48[k] = ab_x[k] * gf_1_33[k]
                      + hf_1_33[k];

            t_49[k] = ab_x[k] * gf_1_34[k]
                      + hf_1_34[k];

            t_50[k] = ab_x[k] * gf_1_35[k]
                      + hf_1_35[k];

            t_51[k] = ab_x[k] * gf_1_36[k]
                      + hf_1_36[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, ab_x, ab_y, gf_1_36, gf_1_37, gf_1_38, \
                         gf_1_39, gf_0_36, hf_1_37, hf_1_38, hf_1_39, \
                         hf_1_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * gf_1_37[k]
                      + hf_1_37[k];

            t_53[k] = ab_x[k] * gf_1_38[k]
                      + hf_1_38[k];

            t_54[k] = ab_x[k] * gf_1_39[k]
                      + hf_1_39[k];

            t_55[k] = ab_y[k] * gf_1_36[k]
                      + gf_0_36[k]
                      + hf_1_66[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, ab_y, ab_z, gf_1_37, gf_1_38, gf_1_39, \
                         gf_0_37, gf_0_38, gf_0_39, hf_1_67, hf_1_68, hf_1_69, \
                         hf_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = ab_y[k] * gf_1_37[k]
                      + gf_0_37[k]
                      + hf_1_67[k];

            t_57[k] = ab_y[k] * gf_1_38[k]
                      + gf_0_38[k]
                      + hf_1_68[k];

            t_58[k] = ab_y[k] * gf_1_39[k]
                      + gf_0_39[k]
                      + hf_1_69[k];

            t_59[k] = ab_z[k] * gf_1_39[k]
                      + hf_1_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gf_1_40, gf_1_41, gf_1_42, \
                         gf_1_43, gf_1_44, hf_1_40, hf_1_41, hf_1_42, hf_1_43, \
                         hf_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * gf_1_40[k]
                      + hf_1_40[k];

            t_61[k] = ab_x[k] * gf_1_41[k]
                      + hf_1_41[k];

            t_62[k] = ab_x[k] * gf_1_42[k]
                      + hf_1_42[k];

            t_63[k] = ab_x[k] * gf_1_43[k]
                      + hf_1_43[k];

            t_64[k] = ab_x[k] * gf_1_44[k]
                      + hf_1_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gf_1_45, gf_1_46, gf_1_47, \
                         gf_1_48, gf_1_49, hf_1_45, hf_1_46, hf_1_47, hf_1_48, \
                         hf_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * gf_1_45[k]
                      + hf_1_45[k];

            t_66[k] = ab_x[k] * gf_1_46[k]
                      + hf_1_46[k];

            t_67[k] = ab_x[k] * gf_1_47[k]
                      + hf_1_47[k];

            t_68[k] = ab_x[k] * gf_1_48[k]
                      + hf_1_48[k];

            t_69[k] = ab_x[k] * gf_1_49[k]
                      + hf_1_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, ab_y, gf_1_46, gf_1_47, gf_1_48, gf_0_46, gf_0_47, \
                         gf_0_48, hf_1_76, hf_1_77, hf_1_78 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * gf_1_46[k]
                      + gf_0_46[k]
                      + hf_1_76[k];

            t_71[k] = ab_y[k] * gf_1_47[k]
                      + gf_0_47[k]
                      + hf_1_77[k];

            t_72[k] = ab_y[k] * gf_1_48[k]
                      + gf_0_48[k]
                      + hf_1_78[k];
        }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, gf_1_49, gf_1_50, gf_1_51, \
                         gf_0_49, hf_1_50, hf_1_51, hf_1_79, hf_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_73[k] = ab_y[k] * gf_1_49[k]
                      + gf_0_49[k]
                      + hf_1_79[k];

            t_74[k] = ab_z[k] * gf_1_49[k]
                      + hf_1_89[k];

            t_75[k] = ab_x[k] * gf_1_50[k]
                      + hf_1_50[k];

            t_76[k] = ab_x[k] * gf_1_51[k]
                      + hf_1_51[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, ab_x, gf_1_52, gf_1_53, gf_1_54, \
                         gf_1_55, gf_1_56, hf_1_52, hf_1_53, hf_1_54, hf_1_55, \
                         hf_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_x[k] * gf_1_52[k]
                      + hf_1_52[k];

            t_78[k] = ab_x[k] * gf_1_53[k]
                      + hf_1_53[k];

            t_79[k] = ab_x[k] * gf_1_54[k]
                      + hf_1_54[k];

            t_80[k] = ab_x[k] * gf_1_55[k]
                      + hf_1_55[k];

            t_81[k] = ab_x[k] * gf_1_56[k]
                      + hf_1_56[k];
        }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, ab_x, ab_y, gf_1_56, gf_1_57, gf_1_58, \
                         gf_1_59, gf_0_56, hf_1_57, hf_1_58, hf_1_59, \
                         hf_1_86 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_82[k] = ab_x[k] * gf_1_57[k]
                      + hf_1_57[k];

            t_83[k] = ab_x[k] * gf_1_58[k]
                      + hf_1_58[k];

            t_84[k] = ab_x[k] * gf_1_59[k]
                      + hf_1_59[k];

            t_85[k] = ab_y[k] * gf_1_56[k]
                      + gf_0_56[k]
                      + hf_1_86[k];
        }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, ab_y, ab_z, gf_1_57, gf_1_58, gf_1_59, \
                         gf_0_57, gf_0_58, gf_0_59, hf_1_87, hf_1_88, hf_1_89, \
                         hf_1_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_86[k] = ab_y[k] * gf_1_57[k]
                      + gf_0_57[k]
                      + hf_1_87[k];

            t_87[k] = ab_y[k] * gf_1_58[k]
                      + gf_0_58[k]
                      + hf_1_88[k];

            t_88[k] = ab_y[k] * gf_1_59[k]
                      + gf_0_59[k]
                      + hf_1_89[k];

            t_89[k] = ab_z[k] * gf_1_59[k]
                      + hf_1_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gf_1_60, gf_1_61, gf_1_62, \
                         gf_1_63, gf_1_64, hf_1_60, hf_1_61, hf_1_62, hf_1_63, \
                         hf_1_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * gf_1_60[k]
                      + hf_1_60[k];

            t_91[k] = ab_x[k] * gf_1_61[k]
                      + hf_1_61[k];

            t_92[k] = ab_x[k] * gf_1_62[k]
                      + hf_1_62[k];

            t_93[k] = ab_x[k] * gf_1_63[k]
                      + hf_1_63[k];

            t_94[k] = ab_x[k] * gf_1_64[k]
                      + hf_1_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gf_1_65, gf_1_66, gf_1_67, \
                         gf_1_68, gf_1_69, hf_1_65, hf_1_66, hf_1_67, hf_1_68, \
                         hf_1_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * gf_1_65[k]
                      + hf_1_65[k];

            t_96[k] = ab_x[k] * gf_1_66[k]
                      + hf_1_66[k];

            t_97[k] = ab_x[k] * gf_1_67[k]
                      + hf_1_67[k];

            t_98[k] = ab_x[k] * gf_1_68[k]
                      + hf_1_68[k];

            t_99[k] = ab_x[k] * gf_1_69[k]
                      + hf_1_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, ab_y, gf_1_66, gf_1_67, gf_1_68, gf_0_66, \
                         gf_0_67, gf_0_68, hf_1_106, hf_1_107, \
                         hf_1_108 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * gf_1_66[k]
                       + gf_0_66[k]
                       + hf_1_106[k];

            t_101[k] = ab_y[k] * gf_1_67[k]
                       + gf_0_67[k]
                       + hf_1_107[k];

            t_102[k] = ab_y[k] * gf_1_68[k]
                       + gf_0_68[k]
                       + hf_1_108[k];
        }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, ab_x, ab_y, ab_z, gf_1_69, gf_1_70, \
                         gf_1_71, gf_0_69, hf_1_70, hf_1_71, hf_1_109, \
                         hf_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_103[k] = ab_y[k] * gf_1_69[k]
                       + gf_0_69[k]
                       + hf_1_109[k];

            t_104[k] = ab_z[k] * gf_1_69[k]
                       + hf_1_119[k];

            t_105[k] = ab_x[k] * gf_1_70[k]
                       + hf_1_70[k];

            t_106[k] = ab_x[k] * gf_1_71[k]
                       + hf_1_71[k];
        }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, ab_x, gf_1_72, gf_1_73, gf_1_74, \
                         gf_1_75, gf_1_76, hf_1_72, hf_1_73, hf_1_74, hf_1_75, \
                         hf_1_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_107[k] = ab_x[k] * gf_1_72[k]
                       + hf_1_72[k];

            t_108[k] = ab_x[k] * gf_1_73[k]
                       + hf_1_73[k];

            t_109[k] = ab_x[k] * gf_1_74[k]
                       + hf_1_74[k];

            t_110[k] = ab_x[k] * gf_1_75[k]
                       + hf_1_75[k];

            t_111[k] = ab_x[k] * gf_1_76[k]
                       + hf_1_76[k];
        }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, ab_x, ab_y, gf_1_76, gf_1_77, gf_1_78, \
                         gf_1_79, gf_0_76, hf_1_77, hf_1_78, hf_1_79, \
                         hf_1_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_112[k] = ab_x[k] * gf_1_77[k]
                       + hf_1_77[k];

            t_113[k] = ab_x[k] * gf_1_78[k]
                       + hf_1_78[k];

            t_114[k] = ab_x[k] * gf_1_79[k]
                       + hf_1_79[k];

            t_115[k] = ab_y[k] * gf_1_76[k]
                       + gf_0_76[k]
                       + hf_1_116[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, ab_y, ab_z, gf_1_77, gf_1_78, gf_1_79, \
                         gf_0_77, gf_0_78, gf_0_79, hf_1_117, hf_1_118, hf_1_119, \
                         hf_1_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_y[k] * gf_1_77[k]
                       + gf_0_77[k]
                       + hf_1_117[k];

            t_117[k] = ab_y[k] * gf_1_78[k]
                       + gf_0_78[k]
                       + hf_1_118[k];

            t_118[k] = ab_y[k] * gf_1_79[k]
                       + gf_0_79[k]
                       + hf_1_119[k];

            t_119[k] = ab_z[k] * gf_1_79[k]
                       + hf_1_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gf_1_80, gf_1_81, gf_1_82, \
                         gf_1_83, gf_1_84, hf_1_80, hf_1_81, hf_1_82, hf_1_83, \
                         hf_1_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * gf_1_80[k]
                       + hf_1_80[k];

            t_121[k] = ab_x[k] * gf_1_81[k]
                       + hf_1_81[k];

            t_122[k] = ab_x[k] * gf_1_82[k]
                       + hf_1_82[k];

            t_123[k] = ab_x[k] * gf_1_83[k]
                       + hf_1_83[k];

            t_124[k] = ab_x[k] * gf_1_84[k]
                       + hf_1_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gf_1_85, gf_1_86, gf_1_87, \
                         gf_1_88, gf_1_89, hf_1_85, hf_1_86, hf_1_87, hf_1_88, \
                         hf_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * gf_1_85[k]
                       + hf_1_85[k];

            t_126[k] = ab_x[k] * gf_1_86[k]
                       + hf_1_86[k];

            t_127[k] = ab_x[k] * gf_1_87[k]
                       + hf_1_87[k];

            t_128[k] = ab_x[k] * gf_1_88[k]
                       + hf_1_88[k];

            t_129[k] = ab_x[k] * gf_1_89[k]
                       + hf_1_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, ab_y, gf_1_86, gf_1_87, gf_1_88, gf_0_86, \
                         gf_0_87, gf_0_88, hf_1_126, hf_1_127, \
                         hf_1_128 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * gf_1_86[k]
                       + gf_0_86[k]
                       + hf_1_126[k];

            t_131[k] = ab_y[k] * gf_1_87[k]
                       + gf_0_87[k]
                       + hf_1_127[k];

            t_132[k] = ab_y[k] * gf_1_88[k]
                       + gf_0_88[k]
                       + hf_1_128[k];
        }
    }
}

static auto
compute_hrr_geom_100y_gg_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t gf_1,
                                             const size_t gf_0, const size_t hf_1,
                                             const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *gf_1_89 = buffer.data(gf_1 + 89 * ncomps + c);
        const auto *gf_1_90 = buffer.data(gf_1 + 90 * ncomps + c);
        const auto *gf_1_91 = buffer.data(gf_1 + 91 * ncomps + c);
        const auto *gf_1_92 = buffer.data(gf_1 + 92 * ncomps + c);
        const auto *gf_1_93 = buffer.data(gf_1 + 93 * ncomps + c);
        const auto *gf_1_94 = buffer.data(gf_1 + 94 * ncomps + c);
        const auto *gf_1_95 = buffer.data(gf_1 + 95 * ncomps + c);
        const auto *gf_1_96 = buffer.data(gf_1 + 96 * ncomps + c);
        const auto *gf_1_97 = buffer.data(gf_1 + 97 * ncomps + c);
        const auto *gf_1_98 = buffer.data(gf_1 + 98 * ncomps + c);
        const auto *gf_1_99 = buffer.data(gf_1 + 99 * ncomps + c);
        const auto *gf_1_100 = buffer.data(gf_1 + 100 * ncomps + c);
        const auto *gf_1_101 = buffer.data(gf_1 + 101 * ncomps + c);
        const auto *gf_1_102 = buffer.data(gf_1 + 102 * ncomps + c);
        const auto *gf_1_103 = buffer.data(gf_1 + 103 * ncomps + c);
        const auto *gf_1_104 = buffer.data(gf_1 + 104 * ncomps + c);
        const auto *gf_1_105 = buffer.data(gf_1 + 105 * ncomps + c);
        const auto *gf_1_106 = buffer.data(gf_1 + 106 * ncomps + c);
        const auto *gf_1_107 = buffer.data(gf_1 + 107 * ncomps + c);
        const auto *gf_1_108 = buffer.data(gf_1 + 108 * ncomps + c);
        const auto *gf_1_109 = buffer.data(gf_1 + 109 * ncomps + c);
        const auto *gf_1_110 = buffer.data(gf_1 + 110 * ncomps + c);
        const auto *gf_1_111 = buffer.data(gf_1 + 111 * ncomps + c);
        const auto *gf_1_112 = buffer.data(gf_1 + 112 * ncomps + c);
        const auto *gf_1_113 = buffer.data(gf_1 + 113 * ncomps + c);
        const auto *gf_1_114 = buffer.data(gf_1 + 114 * ncomps + c);
        const auto *gf_1_115 = buffer.data(gf_1 + 115 * ncomps + c);
        const auto *gf_1_116 = buffer.data(gf_1 + 116 * ncomps + c);
        const auto *gf_1_117 = buffer.data(gf_1 + 117 * ncomps + c);
        const auto *gf_1_118 = buffer.data(gf_1 + 118 * ncomps + c);
        const auto *gf_1_119 = buffer.data(gf_1 + 119 * ncomps + c);
        const auto *gf_1_120 = buffer.data(gf_1 + 120 * ncomps + c);
        const auto *gf_1_121 = buffer.data(gf_1 + 121 * ncomps + c);
        const auto *gf_1_122 = buffer.data(gf_1 + 122 * ncomps + c);
        const auto *gf_1_123 = buffer.data(gf_1 + 123 * ncomps + c);
        const auto *gf_1_124 = buffer.data(gf_1 + 124 * ncomps + c);
        const auto *gf_1_125 = buffer.data(gf_1 + 125 * ncomps + c);
        const auto *gf_1_126 = buffer.data(gf_1 + 126 * ncomps + c);
        const auto *gf_1_127 = buffer.data(gf_1 + 127 * ncomps + c);
        const auto *gf_1_128 = buffer.data(gf_1 + 128 * ncomps + c);
        const auto *gf_1_129 = buffer.data(gf_1 + 129 * ncomps + c);
        const auto *gf_1_130 = buffer.data(gf_1 + 130 * ncomps + c);
        const auto *gf_1_131 = buffer.data(gf_1 + 131 * ncomps + c);
        const auto *gf_1_132 = buffer.data(gf_1 + 132 * ncomps + c);
        const auto *gf_1_133 = buffer.data(gf_1 + 133 * ncomps + c);
        const auto *gf_1_134 = buffer.data(gf_1 + 134 * ncomps + c);
        const auto *gf_1_135 = buffer.data(gf_1 + 135 * ncomps + c);
        const auto *gf_1_136 = buffer.data(gf_1 + 136 * ncomps + c);
        const auto *gf_1_137 = buffer.data(gf_1 + 137 * ncomps + c);
        const auto *gf_1_138 = buffer.data(gf_1 + 138 * ncomps + c);
        const auto *gf_1_139 = buffer.data(gf_1 + 139 * ncomps + c);
        const auto *gf_1_140 = buffer.data(gf_1 + 140 * ncomps + c);
        const auto *gf_1_141 = buffer.data(gf_1 + 141 * ncomps + c);
        const auto *gf_1_142 = buffer.data(gf_1 + 142 * ncomps + c);
        const auto *gf_1_143 = buffer.data(gf_1 + 143 * ncomps + c);
        const auto *gf_1_144 = buffer.data(gf_1 + 144 * ncomps + c);
        const auto *gf_1_145 = buffer.data(gf_1 + 145 * ncomps + c);
        const auto *gf_1_146 = buffer.data(gf_1 + 146 * ncomps + c);
        const auto *gf_1_147 = buffer.data(gf_1 + 147 * ncomps + c);
        const auto *gf_1_148 = buffer.data(gf_1 + 148 * ncomps + c);
        const auto *gf_1_149 = buffer.data(gf_1 + 149 * ncomps + c);

        const auto *gf_0_89 = buffer.data(gf_0 + 89 * ncomps + c);
        const auto *gf_0_96 = buffer.data(gf_0 + 96 * ncomps + c);
        const auto *gf_0_97 = buffer.data(gf_0 + 97 * ncomps + c);
        const auto *gf_0_98 = buffer.data(gf_0 + 98 * ncomps + c);
        const auto *gf_0_99 = buffer.data(gf_0 + 99 * ncomps + c);
        const auto *gf_0_106 = buffer.data(gf_0 + 106 * ncomps + c);
        const auto *gf_0_107 = buffer.data(gf_0 + 107 * ncomps + c);
        const auto *gf_0_108 = buffer.data(gf_0 + 108 * ncomps + c);
        const auto *gf_0_109 = buffer.data(gf_0 + 109 * ncomps + c);
        const auto *gf_0_116 = buffer.data(gf_0 + 116 * ncomps + c);
        const auto *gf_0_117 = buffer.data(gf_0 + 117 * ncomps + c);
        const auto *gf_0_118 = buffer.data(gf_0 + 118 * ncomps + c);
        const auto *gf_0_119 = buffer.data(gf_0 + 119 * ncomps + c);
        const auto *gf_0_126 = buffer.data(gf_0 + 126 * ncomps + c);
        const auto *gf_0_127 = buffer.data(gf_0 + 127 * ncomps + c);
        const auto *gf_0_128 = buffer.data(gf_0 + 128 * ncomps + c);
        const auto *gf_0_129 = buffer.data(gf_0 + 129 * ncomps + c);
        const auto *gf_0_136 = buffer.data(gf_0 + 136 * ncomps + c);
        const auto *gf_0_137 = buffer.data(gf_0 + 137 * ncomps + c);
        const auto *gf_0_138 = buffer.data(gf_0 + 138 * ncomps + c);
        const auto *gf_0_139 = buffer.data(gf_0 + 139 * ncomps + c);
        const auto *gf_0_146 = buffer.data(gf_0 + 146 * ncomps + c);
        const auto *gf_0_147 = buffer.data(gf_0 + 147 * ncomps + c);
        const auto *gf_0_148 = buffer.data(gf_0 + 148 * ncomps + c);
        const auto *gf_0_149 = buffer.data(gf_0 + 149 * ncomps + c);

        const auto *hf_1_90 = buffer.data(hf_1 + 90 * ncomps + c);
        const auto *hf_1_91 = buffer.data(hf_1 + 91 * ncomps + c);
        const auto *hf_1_92 = buffer.data(hf_1 + 92 * ncomps + c);
        const auto *hf_1_93 = buffer.data(hf_1 + 93 * ncomps + c);
        const auto *hf_1_94 = buffer.data(hf_1 + 94 * ncomps + c);
        const auto *hf_1_95 = buffer.data(hf_1 + 95 * ncomps + c);
        const auto *hf_1_96 = buffer.data(hf_1 + 96 * ncomps + c);
        const auto *hf_1_97 = buffer.data(hf_1 + 97 * ncomps + c);
        const auto *hf_1_98 = buffer.data(hf_1 + 98 * ncomps + c);
        const auto *hf_1_99 = buffer.data(hf_1 + 99 * ncomps + c);
        const auto *hf_1_100 = buffer.data(hf_1 + 100 * ncomps + c);
        const auto *hf_1_101 = buffer.data(hf_1 + 101 * ncomps + c);
        const auto *hf_1_102 = buffer.data(hf_1 + 102 * ncomps + c);
        const auto *hf_1_103 = buffer.data(hf_1 + 103 * ncomps + c);
        const auto *hf_1_104 = buffer.data(hf_1 + 104 * ncomps + c);
        const auto *hf_1_105 = buffer.data(hf_1 + 105 * ncomps + c);
        const auto *hf_1_106 = buffer.data(hf_1 + 106 * ncomps + c);
        const auto *hf_1_107 = buffer.data(hf_1 + 107 * ncomps + c);
        const auto *hf_1_108 = buffer.data(hf_1 + 108 * ncomps + c);
        const auto *hf_1_109 = buffer.data(hf_1 + 109 * ncomps + c);
        const auto *hf_1_110 = buffer.data(hf_1 + 110 * ncomps + c);
        const auto *hf_1_111 = buffer.data(hf_1 + 111 * ncomps + c);
        const auto *hf_1_112 = buffer.data(hf_1 + 112 * ncomps + c);
        const auto *hf_1_113 = buffer.data(hf_1 + 113 * ncomps + c);
        const auto *hf_1_114 = buffer.data(hf_1 + 114 * ncomps + c);
        const auto *hf_1_115 = buffer.data(hf_1 + 115 * ncomps + c);
        const auto *hf_1_116 = buffer.data(hf_1 + 116 * ncomps + c);
        const auto *hf_1_117 = buffer.data(hf_1 + 117 * ncomps + c);
        const auto *hf_1_118 = buffer.data(hf_1 + 118 * ncomps + c);
        const auto *hf_1_119 = buffer.data(hf_1 + 119 * ncomps + c);
        const auto *hf_1_120 = buffer.data(hf_1 + 120 * ncomps + c);
        const auto *hf_1_121 = buffer.data(hf_1 + 121 * ncomps + c);
        const auto *hf_1_122 = buffer.data(hf_1 + 122 * ncomps + c);
        const auto *hf_1_123 = buffer.data(hf_1 + 123 * ncomps + c);
        const auto *hf_1_124 = buffer.data(hf_1 + 124 * ncomps + c);
        const auto *hf_1_125 = buffer.data(hf_1 + 125 * ncomps + c);
        const auto *hf_1_126 = buffer.data(hf_1 + 126 * ncomps + c);
        const auto *hf_1_127 = buffer.data(hf_1 + 127 * ncomps + c);
        const auto *hf_1_128 = buffer.data(hf_1 + 128 * ncomps + c);
        const auto *hf_1_129 = buffer.data(hf_1 + 129 * ncomps + c);
        const auto *hf_1_130 = buffer.data(hf_1 + 130 * ncomps + c);
        const auto *hf_1_131 = buffer.data(hf_1 + 131 * ncomps + c);
        const auto *hf_1_132 = buffer.data(hf_1 + 132 * ncomps + c);
        const auto *hf_1_133 = buffer.data(hf_1 + 133 * ncomps + c);
        const auto *hf_1_134 = buffer.data(hf_1 + 134 * ncomps + c);
        const auto *hf_1_135 = buffer.data(hf_1 + 135 * ncomps + c);
        const auto *hf_1_136 = buffer.data(hf_1 + 136 * ncomps + c);
        const auto *hf_1_137 = buffer.data(hf_1 + 137 * ncomps + c);
        const auto *hf_1_138 = buffer.data(hf_1 + 138 * ncomps + c);
        const auto *hf_1_139 = buffer.data(hf_1 + 139 * ncomps + c);
        const auto *hf_1_140 = buffer.data(hf_1 + 140 * ncomps + c);
        const auto *hf_1_141 = buffer.data(hf_1 + 141 * ncomps + c);
        const auto *hf_1_142 = buffer.data(hf_1 + 142 * ncomps + c);
        const auto *hf_1_143 = buffer.data(hf_1 + 143 * ncomps + c);
        const auto *hf_1_144 = buffer.data(hf_1 + 144 * ncomps + c);
        const auto *hf_1_145 = buffer.data(hf_1 + 145 * ncomps + c);
        const auto *hf_1_146 = buffer.data(hf_1 + 146 * ncomps + c);
        const auto *hf_1_147 = buffer.data(hf_1 + 147 * ncomps + c);
        const auto *hf_1_148 = buffer.data(hf_1 + 148 * ncomps + c);
        const auto *hf_1_149 = buffer.data(hf_1 + 149 * ncomps + c);
        const auto *hf_1_156 = buffer.data(hf_1 + 156 * ncomps + c);
        const auto *hf_1_157 = buffer.data(hf_1 + 157 * ncomps + c);
        const auto *hf_1_158 = buffer.data(hf_1 + 158 * ncomps + c);
        const auto *hf_1_159 = buffer.data(hf_1 + 159 * ncomps + c);
        const auto *hf_1_166 = buffer.data(hf_1 + 166 * ncomps + c);
        const auto *hf_1_167 = buffer.data(hf_1 + 167 * ncomps + c);
        const auto *hf_1_168 = buffer.data(hf_1 + 168 * ncomps + c);
        const auto *hf_1_169 = buffer.data(hf_1 + 169 * ncomps + c);
        const auto *hf_1_176 = buffer.data(hf_1 + 176 * ncomps + c);
        const auto *hf_1_177 = buffer.data(hf_1 + 177 * ncomps + c);
        const auto *hf_1_178 = buffer.data(hf_1 + 178 * ncomps + c);
        const auto *hf_1_179 = buffer.data(hf_1 + 179 * ncomps + c);
        const auto *hf_1_186 = buffer.data(hf_1 + 186 * ncomps + c);
        const auto *hf_1_187 = buffer.data(hf_1 + 187 * ncomps + c);
        const auto *hf_1_188 = buffer.data(hf_1 + 188 * ncomps + c);
        const auto *hf_1_189 = buffer.data(hf_1 + 189 * ncomps + c);
        const auto *hf_1_196 = buffer.data(hf_1 + 196 * ncomps + c);
        const auto *hf_1_197 = buffer.data(hf_1 + 197 * ncomps + c);
        const auto *hf_1_198 = buffer.data(hf_1 + 198 * ncomps + c);
        const auto *hf_1_199 = buffer.data(hf_1 + 199 * ncomps + c);
        const auto *hf_1_209 = buffer.data(hf_1 + 209 * ncomps + c);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, ab_x, ab_y, ab_z, gf_1_89, gf_1_90, \
                         gf_1_91, gf_0_89, hf_1_90, hf_1_91, hf_1_129, \
                         hf_1_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_133[k] = ab_y[k] * gf_1_89[k]
                       + gf_0_89[k]
                       + hf_1_129[k];

            t_134[k] = ab_z[k] * gf_1_89[k]
                       + hf_1_139[k];

            t_135[k] = ab_x[k] * gf_1_90[k]
                       + hf_1_90[k];

            t_136[k] = ab_x[k] * gf_1_91[k]
                       + hf_1_91[k];
        }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, ab_x, gf_1_92, gf_1_93, gf_1_94, \
                         gf_1_95, gf_1_96, hf_1_92, hf_1_93, hf_1_94, hf_1_95, \
                         hf_1_96 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_137[k] = ab_x[k] * gf_1_92[k]
                       + hf_1_92[k];

            t_138[k] = ab_x[k] * gf_1_93[k]
                       + hf_1_93[k];

            t_139[k] = ab_x[k] * gf_1_94[k]
                       + hf_1_94[k];

            t_140[k] = ab_x[k] * gf_1_95[k]
                       + hf_1_95[k];

            t_141[k] = ab_x[k] * gf_1_96[k]
                       + hf_1_96[k];
        }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, ab_x, ab_y, gf_1_96, gf_1_97, gf_1_98, \
                         gf_1_99, gf_0_96, hf_1_97, hf_1_98, hf_1_99, \
                         hf_1_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_142[k] = ab_x[k] * gf_1_97[k]
                       + hf_1_97[k];

            t_143[k] = ab_x[k] * gf_1_98[k]
                       + hf_1_98[k];

            t_144[k] = ab_x[k] * gf_1_99[k]
                       + hf_1_99[k];

            t_145[k] = ab_y[k] * gf_1_96[k]
                       + gf_0_96[k]
                       + hf_1_136[k];
        }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, ab_y, ab_z, gf_1_97, gf_1_98, gf_1_99, \
                         gf_0_97, gf_0_98, gf_0_99, hf_1_137, hf_1_138, hf_1_139, \
                         hf_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_146[k] = ab_y[k] * gf_1_97[k]
                       + gf_0_97[k]
                       + hf_1_137[k];

            t_147[k] = ab_y[k] * gf_1_98[k]
                       + gf_0_98[k]
                       + hf_1_138[k];

            t_148[k] = ab_y[k] * gf_1_99[k]
                       + gf_0_99[k]
                       + hf_1_139[k];

            t_149[k] = ab_z[k] * gf_1_99[k]
                       + hf_1_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gf_1_100, gf_1_101, \
                         gf_1_102, gf_1_103, gf_1_104, hf_1_100, hf_1_101, hf_1_102, hf_1_103, \
                         hf_1_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * gf_1_100[k]
                       + hf_1_100[k];

            t_151[k] = ab_x[k] * gf_1_101[k]
                       + hf_1_101[k];

            t_152[k] = ab_x[k] * gf_1_102[k]
                       + hf_1_102[k];

            t_153[k] = ab_x[k] * gf_1_103[k]
                       + hf_1_103[k];

            t_154[k] = ab_x[k] * gf_1_104[k]
                       + hf_1_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gf_1_105, gf_1_106, \
                         gf_1_107, gf_1_108, gf_1_109, hf_1_105, hf_1_106, hf_1_107, hf_1_108, \
                         hf_1_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * gf_1_105[k]
                       + hf_1_105[k];

            t_156[k] = ab_x[k] * gf_1_106[k]
                       + hf_1_106[k];

            t_157[k] = ab_x[k] * gf_1_107[k]
                       + hf_1_107[k];

            t_158[k] = ab_x[k] * gf_1_108[k]
                       + hf_1_108[k];

            t_159[k] = ab_x[k] * gf_1_109[k]
                       + hf_1_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, ab_y, gf_1_106, gf_1_107, gf_1_108, gf_0_106, \
                         gf_0_107, gf_0_108, hf_1_156, hf_1_157, \
                         hf_1_158 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * gf_1_106[k]
                       + gf_0_106[k]
                       + hf_1_156[k];

            t_161[k] = ab_y[k] * gf_1_107[k]
                       + gf_0_107[k]
                       + hf_1_157[k];

            t_162[k] = ab_y[k] * gf_1_108[k]
                       + gf_0_108[k]
                       + hf_1_158[k];
        }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, ab_x, ab_y, ab_z, gf_1_109, gf_1_110, \
                         gf_1_111, gf_0_109, hf_1_110, hf_1_111, hf_1_159, \
                         hf_1_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_163[k] = ab_y[k] * gf_1_109[k]
                       + gf_0_109[k]
                       + hf_1_159[k];

            t_164[k] = ab_z[k] * gf_1_109[k]
                       + hf_1_169[k];

            t_165[k] = ab_x[k] * gf_1_110[k]
                       + hf_1_110[k];

            t_166[k] = ab_x[k] * gf_1_111[k]
                       + hf_1_111[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, ab_x, gf_1_112, gf_1_113, \
                         gf_1_114, gf_1_115, gf_1_116, hf_1_112, hf_1_113, hf_1_114, hf_1_115, \
                         hf_1_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = ab_x[k] * gf_1_112[k]
                       + hf_1_112[k];

            t_168[k] = ab_x[k] * gf_1_113[k]
                       + hf_1_113[k];

            t_169[k] = ab_x[k] * gf_1_114[k]
                       + hf_1_114[k];

            t_170[k] = ab_x[k] * gf_1_115[k]
                       + hf_1_115[k];

            t_171[k] = ab_x[k] * gf_1_116[k]
                       + hf_1_116[k];
        }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, ab_x, ab_y, gf_1_116, gf_1_117, gf_1_118, \
                         gf_1_119, gf_0_116, hf_1_117, hf_1_118, hf_1_119, \
                         hf_1_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_172[k] = ab_x[k] * gf_1_117[k]
                       + hf_1_117[k];

            t_173[k] = ab_x[k] * gf_1_118[k]
                       + hf_1_118[k];

            t_174[k] = ab_x[k] * gf_1_119[k]
                       + hf_1_119[k];

            t_175[k] = ab_y[k] * gf_1_116[k]
                       + gf_0_116[k]
                       + hf_1_166[k];
        }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, ab_y, ab_z, gf_1_117, gf_1_118, gf_1_119, \
                         gf_0_117, gf_0_118, gf_0_119, hf_1_167, hf_1_168, hf_1_169, \
                         hf_1_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_176[k] = ab_y[k] * gf_1_117[k]
                       + gf_0_117[k]
                       + hf_1_167[k];

            t_177[k] = ab_y[k] * gf_1_118[k]
                       + gf_0_118[k]
                       + hf_1_168[k];

            t_178[k] = ab_y[k] * gf_1_119[k]
                       + gf_0_119[k]
                       + hf_1_169[k];

            t_179[k] = ab_z[k] * gf_1_119[k]
                       + hf_1_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gf_1_120, gf_1_121, \
                         gf_1_122, gf_1_123, gf_1_124, hf_1_120, hf_1_121, hf_1_122, hf_1_123, \
                         hf_1_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * gf_1_120[k]
                       + hf_1_120[k];

            t_181[k] = ab_x[k] * gf_1_121[k]
                       + hf_1_121[k];

            t_182[k] = ab_x[k] * gf_1_122[k]
                       + hf_1_122[k];

            t_183[k] = ab_x[k] * gf_1_123[k]
                       + hf_1_123[k];

            t_184[k] = ab_x[k] * gf_1_124[k]
                       + hf_1_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gf_1_125, gf_1_126, \
                         gf_1_127, gf_1_128, gf_1_129, hf_1_125, hf_1_126, hf_1_127, hf_1_128, \
                         hf_1_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * gf_1_125[k]
                       + hf_1_125[k];

            t_186[k] = ab_x[k] * gf_1_126[k]
                       + hf_1_126[k];

            t_187[k] = ab_x[k] * gf_1_127[k]
                       + hf_1_127[k];

            t_188[k] = ab_x[k] * gf_1_128[k]
                       + hf_1_128[k];

            t_189[k] = ab_x[k] * gf_1_129[k]
                       + hf_1_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, ab_y, gf_1_126, gf_1_127, gf_1_128, gf_0_126, \
                         gf_0_127, gf_0_128, hf_1_176, hf_1_177, \
                         hf_1_178 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * gf_1_126[k]
                       + gf_0_126[k]
                       + hf_1_176[k];

            t_191[k] = ab_y[k] * gf_1_127[k]
                       + gf_0_127[k]
                       + hf_1_177[k];

            t_192[k] = ab_y[k] * gf_1_128[k]
                       + gf_0_128[k]
                       + hf_1_178[k];
        }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, ab_x, ab_y, ab_z, gf_1_129, gf_1_130, \
                         gf_1_131, gf_0_129, hf_1_130, hf_1_131, hf_1_179, \
                         hf_1_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_193[k] = ab_y[k] * gf_1_129[k]
                       + gf_0_129[k]
                       + hf_1_179[k];

            t_194[k] = ab_z[k] * gf_1_129[k]
                       + hf_1_189[k];

            t_195[k] = ab_x[k] * gf_1_130[k]
                       + hf_1_130[k];

            t_196[k] = ab_x[k] * gf_1_131[k]
                       + hf_1_131[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_x, gf_1_132, gf_1_133, \
                         gf_1_134, gf_1_135, gf_1_136, hf_1_132, hf_1_133, hf_1_134, hf_1_135, \
                         hf_1_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_x[k] * gf_1_132[k]
                       + hf_1_132[k];

            t_198[k] = ab_x[k] * gf_1_133[k]
                       + hf_1_133[k];

            t_199[k] = ab_x[k] * gf_1_134[k]
                       + hf_1_134[k];

            t_200[k] = ab_x[k] * gf_1_135[k]
                       + hf_1_135[k];

            t_201[k] = ab_x[k] * gf_1_136[k]
                       + hf_1_136[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, ab_x, ab_y, gf_1_136, gf_1_137, gf_1_138, \
                         gf_1_139, gf_0_136, hf_1_137, hf_1_138, hf_1_139, \
                         hf_1_186 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_x[k] * gf_1_137[k]
                       + hf_1_137[k];

            t_203[k] = ab_x[k] * gf_1_138[k]
                       + hf_1_138[k];

            t_204[k] = ab_x[k] * gf_1_139[k]
                       + hf_1_139[k];

            t_205[k] = ab_y[k] * gf_1_136[k]
                       + gf_0_136[k]
                       + hf_1_186[k];
        }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, ab_y, ab_z, gf_1_137, gf_1_138, gf_1_139, \
                         gf_0_137, gf_0_138, gf_0_139, hf_1_187, hf_1_188, hf_1_189, \
                         hf_1_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_206[k] = ab_y[k] * gf_1_137[k]
                       + gf_0_137[k]
                       + hf_1_187[k];

            t_207[k] = ab_y[k] * gf_1_138[k]
                       + gf_0_138[k]
                       + hf_1_188[k];

            t_208[k] = ab_y[k] * gf_1_139[k]
                       + gf_0_139[k]
                       + hf_1_189[k];

            t_209[k] = ab_z[k] * gf_1_139[k]
                       + hf_1_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gf_1_140, gf_1_141, \
                         gf_1_142, gf_1_143, gf_1_144, hf_1_140, hf_1_141, hf_1_142, hf_1_143, \
                         hf_1_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * gf_1_140[k]
                       + hf_1_140[k];

            t_211[k] = ab_x[k] * gf_1_141[k]
                       + hf_1_141[k];

            t_212[k] = ab_x[k] * gf_1_142[k]
                       + hf_1_142[k];

            t_213[k] = ab_x[k] * gf_1_143[k]
                       + hf_1_143[k];

            t_214[k] = ab_x[k] * gf_1_144[k]
                       + hf_1_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gf_1_145, gf_1_146, \
                         gf_1_147, gf_1_148, gf_1_149, hf_1_145, hf_1_146, hf_1_147, hf_1_148, \
                         hf_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * gf_1_145[k]
                       + hf_1_145[k];

            t_216[k] = ab_x[k] * gf_1_146[k]
                       + hf_1_146[k];

            t_217[k] = ab_x[k] * gf_1_147[k]
                       + hf_1_147[k];

            t_218[k] = ab_x[k] * gf_1_148[k]
                       + hf_1_148[k];

            t_219[k] = ab_x[k] * gf_1_149[k]
                       + hf_1_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, ab_y, gf_1_146, gf_1_147, gf_1_148, gf_0_146, \
                         gf_0_147, gf_0_148, hf_1_196, hf_1_197, \
                         hf_1_198 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * gf_1_146[k]
                       + gf_0_146[k]
                       + hf_1_196[k];

            t_221[k] = ab_y[k] * gf_1_147[k]
                       + gf_0_147[k]
                       + hf_1_197[k];

            t_222[k] = ab_y[k] * gf_1_148[k]
                       + gf_0_148[k]
                       + hf_1_198[k];
        }

#pragma omp simd aligned(t_223, t_224, ab_y, ab_z, gf_1_149, gf_0_149, hf_1_199, \
                         hf_1_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_223[k] = ab_y[k] * gf_1_149[k]
                       + gf_0_149[k]
                       + hf_1_199[k];

            t_224[k] = ab_z[k] * gf_1_149[k]
                       + hf_1_209[k];
        }
    }
}

auto
compute_hrr_geom_100y_gg_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t gf_1, const size_t gf_0,
                                      const size_t hf_1, const size_t ncomps,
                                      const size_t nmax) -> void
{
    compute_hrr_geom_100y_gg_out_of_first_piece0(buffer, coordinates, target, gf_1, gf_0, hf_1,
                                                 ncomps, nmax);

    compute_hrr_geom_100y_gg_out_of_first_piece1(buffer, coordinates, target, gf_1, gf_0, hf_1,
                                                 ncomps, nmax);
}

}  // namespace simdtrf
