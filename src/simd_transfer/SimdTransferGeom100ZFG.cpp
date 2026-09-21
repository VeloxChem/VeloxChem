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


#include "SimdTransferGeom100ZFG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_100z_fg_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t ff_1,
                                             const size_t ff_0, const size_t gf_1,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ff_1_0 = buffer.data(ff_1 + 0 * ncomps + c);
        const auto *ff_1_1 = buffer.data(ff_1 + 1 * ncomps + c);
        const auto *ff_1_2 = buffer.data(ff_1 + 2 * ncomps + c);
        const auto *ff_1_3 = buffer.data(ff_1 + 3 * ncomps + c);
        const auto *ff_1_4 = buffer.data(ff_1 + 4 * ncomps + c);
        const auto *ff_1_5 = buffer.data(ff_1 + 5 * ncomps + c);
        const auto *ff_1_6 = buffer.data(ff_1 + 6 * ncomps + c);
        const auto *ff_1_7 = buffer.data(ff_1 + 7 * ncomps + c);
        const auto *ff_1_8 = buffer.data(ff_1 + 8 * ncomps + c);
        const auto *ff_1_9 = buffer.data(ff_1 + 9 * ncomps + c);
        const auto *ff_1_10 = buffer.data(ff_1 + 10 * ncomps + c);
        const auto *ff_1_11 = buffer.data(ff_1 + 11 * ncomps + c);
        const auto *ff_1_12 = buffer.data(ff_1 + 12 * ncomps + c);
        const auto *ff_1_13 = buffer.data(ff_1 + 13 * ncomps + c);
        const auto *ff_1_14 = buffer.data(ff_1 + 14 * ncomps + c);
        const auto *ff_1_15 = buffer.data(ff_1 + 15 * ncomps + c);
        const auto *ff_1_16 = buffer.data(ff_1 + 16 * ncomps + c);
        const auto *ff_1_17 = buffer.data(ff_1 + 17 * ncomps + c);
        const auto *ff_1_18 = buffer.data(ff_1 + 18 * ncomps + c);
        const auto *ff_1_19 = buffer.data(ff_1 + 19 * ncomps + c);
        const auto *ff_1_20 = buffer.data(ff_1 + 20 * ncomps + c);
        const auto *ff_1_21 = buffer.data(ff_1 + 21 * ncomps + c);
        const auto *ff_1_22 = buffer.data(ff_1 + 22 * ncomps + c);
        const auto *ff_1_23 = buffer.data(ff_1 + 23 * ncomps + c);
        const auto *ff_1_24 = buffer.data(ff_1 + 24 * ncomps + c);
        const auto *ff_1_25 = buffer.data(ff_1 + 25 * ncomps + c);
        const auto *ff_1_26 = buffer.data(ff_1 + 26 * ncomps + c);
        const auto *ff_1_27 = buffer.data(ff_1 + 27 * ncomps + c);
        const auto *ff_1_28 = buffer.data(ff_1 + 28 * ncomps + c);
        const auto *ff_1_29 = buffer.data(ff_1 + 29 * ncomps + c);
        const auto *ff_1_30 = buffer.data(ff_1 + 30 * ncomps + c);
        const auto *ff_1_31 = buffer.data(ff_1 + 31 * ncomps + c);
        const auto *ff_1_32 = buffer.data(ff_1 + 32 * ncomps + c);
        const auto *ff_1_33 = buffer.data(ff_1 + 33 * ncomps + c);
        const auto *ff_1_34 = buffer.data(ff_1 + 34 * ncomps + c);
        const auto *ff_1_35 = buffer.data(ff_1 + 35 * ncomps + c);
        const auto *ff_1_36 = buffer.data(ff_1 + 36 * ncomps + c);
        const auto *ff_1_37 = buffer.data(ff_1 + 37 * ncomps + c);
        const auto *ff_1_38 = buffer.data(ff_1 + 38 * ncomps + c);
        const auto *ff_1_39 = buffer.data(ff_1 + 39 * ncomps + c);
        const auto *ff_1_40 = buffer.data(ff_1 + 40 * ncomps + c);
        const auto *ff_1_41 = buffer.data(ff_1 + 41 * ncomps + c);
        const auto *ff_1_42 = buffer.data(ff_1 + 42 * ncomps + c);
        const auto *ff_1_43 = buffer.data(ff_1 + 43 * ncomps + c);
        const auto *ff_1_44 = buffer.data(ff_1 + 44 * ncomps + c);
        const auto *ff_1_45 = buffer.data(ff_1 + 45 * ncomps + c);
        const auto *ff_1_46 = buffer.data(ff_1 + 46 * ncomps + c);
        const auto *ff_1_47 = buffer.data(ff_1 + 47 * ncomps + c);
        const auto *ff_1_48 = buffer.data(ff_1 + 48 * ncomps + c);
        const auto *ff_1_49 = buffer.data(ff_1 + 49 * ncomps + c);
        const auto *ff_1_50 = buffer.data(ff_1 + 50 * ncomps + c);
        const auto *ff_1_51 = buffer.data(ff_1 + 51 * ncomps + c);
        const auto *ff_1_52 = buffer.data(ff_1 + 52 * ncomps + c);
        const auto *ff_1_53 = buffer.data(ff_1 + 53 * ncomps + c);
        const auto *ff_1_54 = buffer.data(ff_1 + 54 * ncomps + c);
        const auto *ff_1_55 = buffer.data(ff_1 + 55 * ncomps + c);
        const auto *ff_1_56 = buffer.data(ff_1 + 56 * ncomps + c);
        const auto *ff_1_57 = buffer.data(ff_1 + 57 * ncomps + c);
        const auto *ff_1_58 = buffer.data(ff_1 + 58 * ncomps + c);
        const auto *ff_1_59 = buffer.data(ff_1 + 59 * ncomps + c);
        const auto *ff_1_60 = buffer.data(ff_1 + 60 * ncomps + c);
        const auto *ff_1_61 = buffer.data(ff_1 + 61 * ncomps + c);
        const auto *ff_1_62 = buffer.data(ff_1 + 62 * ncomps + c);
        const auto *ff_1_63 = buffer.data(ff_1 + 63 * ncomps + c);
        const auto *ff_1_64 = buffer.data(ff_1 + 64 * ncomps + c);
        const auto *ff_1_65 = buffer.data(ff_1 + 65 * ncomps + c);
        const auto *ff_1_66 = buffer.data(ff_1 + 66 * ncomps + c);
        const auto *ff_1_67 = buffer.data(ff_1 + 67 * ncomps + c);
        const auto *ff_1_68 = buffer.data(ff_1 + 68 * ncomps + c);
        const auto *ff_1_69 = buffer.data(ff_1 + 69 * ncomps + c);
        const auto *ff_1_70 = buffer.data(ff_1 + 70 * ncomps + c);
        const auto *ff_1_71 = buffer.data(ff_1 + 71 * ncomps + c);
        const auto *ff_1_72 = buffer.data(ff_1 + 72 * ncomps + c);
        const auto *ff_1_73 = buffer.data(ff_1 + 73 * ncomps + c);
        const auto *ff_1_74 = buffer.data(ff_1 + 74 * ncomps + c);
        const auto *ff_1_75 = buffer.data(ff_1 + 75 * ncomps + c);
        const auto *ff_1_76 = buffer.data(ff_1 + 76 * ncomps + c);
        const auto *ff_1_77 = buffer.data(ff_1 + 77 * ncomps + c);
        const auto *ff_1_78 = buffer.data(ff_1 + 78 * ncomps + c);
        const auto *ff_1_79 = buffer.data(ff_1 + 79 * ncomps + c);
        const auto *ff_1_80 = buffer.data(ff_1 + 80 * ncomps + c);
        const auto *ff_1_81 = buffer.data(ff_1 + 81 * ncomps + c);
        const auto *ff_1_82 = buffer.data(ff_1 + 82 * ncomps + c);
        const auto *ff_1_83 = buffer.data(ff_1 + 83 * ncomps + c);
        const auto *ff_1_84 = buffer.data(ff_1 + 84 * ncomps + c);
        const auto *ff_1_85 = buffer.data(ff_1 + 85 * ncomps + c);
        const auto *ff_1_86 = buffer.data(ff_1 + 86 * ncomps + c);
        const auto *ff_1_87 = buffer.data(ff_1 + 87 * ncomps + c);
        const auto *ff_1_88 = buffer.data(ff_1 + 88 * ncomps + c);
        const auto *ff_1_89 = buffer.data(ff_1 + 89 * ncomps + c);
        const auto *ff_1_90 = buffer.data(ff_1 + 90 * ncomps + c);
        const auto *ff_1_91 = buffer.data(ff_1 + 91 * ncomps + c);
        const auto *ff_1_92 = buffer.data(ff_1 + 92 * ncomps + c);

        const auto *ff_0_9 = buffer.data(ff_0 + 9 * ncomps + c);
        const auto *ff_0_19 = buffer.data(ff_0 + 19 * ncomps + c);
        const auto *ff_0_29 = buffer.data(ff_0 + 29 * ncomps + c);
        const auto *ff_0_39 = buffer.data(ff_0 + 39 * ncomps + c);
        const auto *ff_0_49 = buffer.data(ff_0 + 49 * ncomps + c);
        const auto *ff_0_59 = buffer.data(ff_0 + 59 * ncomps + c);
        const auto *ff_0_69 = buffer.data(ff_0 + 69 * ncomps + c);
        const auto *ff_0_79 = buffer.data(ff_0 + 79 * ncomps + c);
        const auto *ff_0_89 = buffer.data(ff_0 + 89 * ncomps + c);

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
        const auto *gf_1_90 = buffer.data(gf_1 + 90 * ncomps + c);
        const auto *gf_1_91 = buffer.data(gf_1 + 91 * ncomps + c);
        const auto *gf_1_92 = buffer.data(gf_1 + 92 * ncomps + c);
        const auto *gf_1_99 = buffer.data(gf_1 + 99 * ncomps + c);
        const auto *gf_1_106 = buffer.data(gf_1 + 106 * ncomps + c);
        const auto *gf_1_107 = buffer.data(gf_1 + 107 * ncomps + c);
        const auto *gf_1_108 = buffer.data(gf_1 + 108 * ncomps + c);
        const auto *gf_1_109 = buffer.data(gf_1 + 109 * ncomps + c);
        const auto *gf_1_116 = buffer.data(gf_1 + 116 * ncomps + c);
        const auto *gf_1_117 = buffer.data(gf_1 + 117 * ncomps + c);
        const auto *gf_1_118 = buffer.data(gf_1 + 118 * ncomps + c);
        const auto *gf_1_119 = buffer.data(gf_1 + 119 * ncomps + c);
        const auto *gf_1_126 = buffer.data(gf_1 + 126 * ncomps + c);
        const auto *gf_1_127 = buffer.data(gf_1 + 127 * ncomps + c);
        const auto *gf_1_128 = buffer.data(gf_1 + 128 * ncomps + c);
        const auto *gf_1_129 = buffer.data(gf_1 + 129 * ncomps + c);
        const auto *gf_1_139 = buffer.data(gf_1 + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ff_1_0, ff_1_1, ff_1_2, ff_1_3, \
                         ff_1_4, gf_1_0, gf_1_1, gf_1_2, gf_1_3, \
                         gf_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ff_1_0[k]
                     + gf_1_0[k];

            t_1[k] = ab_x[k] * ff_1_1[k]
                     + gf_1_1[k];

            t_2[k] = ab_x[k] * ff_1_2[k]
                     + gf_1_2[k];

            t_3[k] = ab_x[k] * ff_1_3[k]
                     + gf_1_3[k];

            t_4[k] = ab_x[k] * ff_1_4[k]
                     + gf_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ff_1_5, ff_1_6, ff_1_7, ff_1_8, \
                         ff_1_9, gf_1_5, gf_1_6, gf_1_7, gf_1_8, \
                         gf_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * ff_1_5[k]
                     + gf_1_5[k];

            t_6[k] = ab_x[k] * ff_1_6[k]
                     + gf_1_6[k];

            t_7[k] = ab_x[k] * ff_1_7[k]
                     + gf_1_7[k];

            t_8[k] = ab_x[k] * ff_1_8[k]
                     + gf_1_8[k];

            t_9[k] = ab_x[k] * ff_1_9[k]
                     + gf_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, ab_y, ff_1_6, ff_1_7, ff_1_8, ff_1_9, \
                         gf_1_16, gf_1_17, gf_1_18, gf_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * ff_1_6[k]
                      + gf_1_16[k];

            t_11[k] = ab_y[k] * ff_1_7[k]
                      + gf_1_17[k];

            t_12[k] = ab_y[k] * ff_1_8[k]
                      + gf_1_18[k];

            t_13[k] = ab_y[k] * ff_1_9[k]
                      + gf_1_19[k];
        }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, ab_x, ab_z, ff_1_9, ff_1_10, ff_1_11, \
                         ff_1_12, ff_0_9, gf_1_10, gf_1_11, gf_1_12, \
                         gf_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_14[k] = ab_z[k] * ff_1_9[k]
                      + ff_0_9[k]
                      + gf_1_29[k];

            t_15[k] = ab_x[k] * ff_1_10[k]
                      + gf_1_10[k];

            t_16[k] = ab_x[k] * ff_1_11[k]
                      + gf_1_11[k];

            t_17[k] = ab_x[k] * ff_1_12[k]
                      + gf_1_12[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ff_1_13, ff_1_14, ff_1_15, \
                         ff_1_16, ff_1_17, gf_1_13, gf_1_14, gf_1_15, gf_1_16, \
                         gf_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * ff_1_13[k]
                      + gf_1_13[k];

            t_19[k] = ab_x[k] * ff_1_14[k]
                      + gf_1_14[k];

            t_20[k] = ab_x[k] * ff_1_15[k]
                      + gf_1_15[k];

            t_21[k] = ab_x[k] * ff_1_16[k]
                      + gf_1_16[k];

            t_22[k] = ab_x[k] * ff_1_17[k]
                      + gf_1_17[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, ab_x, ab_y, ff_1_16, ff_1_17, ff_1_18, \
                         ff_1_19, gf_1_18, gf_1_19, gf_1_36, gf_1_37, \
                         gf_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_x[k] * ff_1_18[k]
                      + gf_1_18[k];

            t_24[k] = ab_x[k] * ff_1_19[k]
                      + gf_1_19[k];

            t_25[k] = ab_y[k] * ff_1_16[k]
                      + gf_1_36[k];

            t_26[k] = ab_y[k] * ff_1_17[k]
                      + gf_1_37[k];

            t_27[k] = ab_y[k] * ff_1_18[k]
                      + gf_1_38[k];
        }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ff_1_19, ff_1_20, ff_1_21, \
                         ff_0_19, gf_1_20, gf_1_21, gf_1_39, gf_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_28[k] = ab_y[k] * ff_1_19[k]
                      + gf_1_39[k];

            t_29[k] = ab_z[k] * ff_1_19[k]
                      + ff_0_19[k]
                      + gf_1_49[k];

            t_30[k] = ab_x[k] * ff_1_20[k]
                      + gf_1_20[k];

            t_31[k] = ab_x[k] * ff_1_21[k]
                      + gf_1_21[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, ab_x, ff_1_22, ff_1_23, ff_1_24, \
                         ff_1_25, ff_1_26, gf_1_22, gf_1_23, gf_1_24, gf_1_25, \
                         gf_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_x[k] * ff_1_22[k]
                      + gf_1_22[k];

            t_33[k] = ab_x[k] * ff_1_23[k]
                      + gf_1_23[k];

            t_34[k] = ab_x[k] * ff_1_24[k]
                      + gf_1_24[k];

            t_35[k] = ab_x[k] * ff_1_25[k]
                      + gf_1_25[k];

            t_36[k] = ab_x[k] * ff_1_26[k]
                      + gf_1_26[k];
        }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, ab_x, ab_y, ff_1_26, ff_1_27, ff_1_28, \
                         ff_1_29, gf_1_27, gf_1_28, gf_1_29, gf_1_46, \
                         gf_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_37[k] = ab_x[k] * ff_1_27[k]
                      + gf_1_27[k];

            t_38[k] = ab_x[k] * ff_1_28[k]
                      + gf_1_28[k];

            t_39[k] = ab_x[k] * ff_1_29[k]
                      + gf_1_29[k];

            t_40[k] = ab_y[k] * ff_1_26[k]
                      + gf_1_46[k];

            t_41[k] = ab_y[k] * ff_1_27[k]
                      + gf_1_47[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, ab_x, ab_y, ab_z, ff_1_28, ff_1_29, ff_1_30, \
                         ff_0_29, gf_1_30, gf_1_48, gf_1_49, gf_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_y[k] * ff_1_28[k]
                      + gf_1_48[k];

            t_43[k] = ab_y[k] * ff_1_29[k]
                      + gf_1_49[k];

            t_44[k] = ab_z[k] * ff_1_29[k]
                      + ff_0_29[k]
                      + gf_1_59[k];

            t_45[k] = ab_x[k] * ff_1_30[k]
                      + gf_1_30[k];
        }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, ab_x, ff_1_31, ff_1_32, ff_1_33, \
                         ff_1_34, ff_1_35, gf_1_31, gf_1_32, gf_1_33, gf_1_34, \
                         gf_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_46[k] = ab_x[k] * ff_1_31[k]
                      + gf_1_31[k];

            t_47[k] = ab_x[k] * ff_1_32[k]
                      + gf_1_32[k];

            t_48[k] = ab_x[k] * ff_1_33[k]
                      + gf_1_33[k];

            t_49[k] = ab_x[k] * ff_1_34[k]
                      + gf_1_34[k];

            t_50[k] = ab_x[k] * ff_1_35[k]
                      + gf_1_35[k];
        }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, ab_x, ab_y, ff_1_36, ff_1_37, ff_1_38, \
                         ff_1_39, gf_1_36, gf_1_37, gf_1_38, gf_1_39, \
                         gf_1_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_51[k] = ab_x[k] * ff_1_36[k]
                      + gf_1_36[k];

            t_52[k] = ab_x[k] * ff_1_37[k]
                      + gf_1_37[k];

            t_53[k] = ab_x[k] * ff_1_38[k]
                      + gf_1_38[k];

            t_54[k] = ab_x[k] * ff_1_39[k]
                      + gf_1_39[k];

            t_55[k] = ab_y[k] * ff_1_36[k]
                      + gf_1_66[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, ab_y, ab_z, ff_1_37, ff_1_38, ff_1_39, \
                         ff_0_39, gf_1_67, gf_1_68, gf_1_69, gf_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = ab_y[k] * ff_1_37[k]
                      + gf_1_67[k];

            t_57[k] = ab_y[k] * ff_1_38[k]
                      + gf_1_68[k];

            t_58[k] = ab_y[k] * ff_1_39[k]
                      + gf_1_69[k];

            t_59[k] = ab_z[k] * ff_1_39[k]
                      + ff_0_39[k]
                      + gf_1_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ff_1_40, ff_1_41, ff_1_42, \
                         ff_1_43, ff_1_44, gf_1_40, gf_1_41, gf_1_42, gf_1_43, \
                         gf_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ff_1_40[k]
                      + gf_1_40[k];

            t_61[k] = ab_x[k] * ff_1_41[k]
                      + gf_1_41[k];

            t_62[k] = ab_x[k] * ff_1_42[k]
                      + gf_1_42[k];

            t_63[k] = ab_x[k] * ff_1_43[k]
                      + gf_1_43[k];

            t_64[k] = ab_x[k] * ff_1_44[k]
                      + gf_1_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ff_1_45, ff_1_46, ff_1_47, \
                         ff_1_48, ff_1_49, gf_1_45, gf_1_46, gf_1_47, gf_1_48, \
                         gf_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * ff_1_45[k]
                      + gf_1_45[k];

            t_66[k] = ab_x[k] * ff_1_46[k]
                      + gf_1_46[k];

            t_67[k] = ab_x[k] * ff_1_47[k]
                      + gf_1_47[k];

            t_68[k] = ab_x[k] * ff_1_48[k]
                      + gf_1_48[k];

            t_69[k] = ab_x[k] * ff_1_49[k]
                      + gf_1_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, ab_y, ff_1_46, ff_1_47, ff_1_48, ff_1_49, \
                         gf_1_76, gf_1_77, gf_1_78, gf_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * ff_1_46[k]
                      + gf_1_76[k];

            t_71[k] = ab_y[k] * ff_1_47[k]
                      + gf_1_77[k];

            t_72[k] = ab_y[k] * ff_1_48[k]
                      + gf_1_78[k];

            t_73[k] = ab_y[k] * ff_1_49[k]
                      + gf_1_79[k];
        }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, ab_x, ab_z, ff_1_49, ff_1_50, ff_1_51, \
                         ff_1_52, ff_0_49, gf_1_50, gf_1_51, gf_1_52, \
                         gf_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_74[k] = ab_z[k] * ff_1_49[k]
                      + ff_0_49[k]
                      + gf_1_89[k];

            t_75[k] = ab_x[k] * ff_1_50[k]
                      + gf_1_50[k];

            t_76[k] = ab_x[k] * ff_1_51[k]
                      + gf_1_51[k];

            t_77[k] = ab_x[k] * ff_1_52[k]
                      + gf_1_52[k];
        }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, ab_x, ff_1_53, ff_1_54, ff_1_55, \
                         ff_1_56, ff_1_57, gf_1_53, gf_1_54, gf_1_55, gf_1_56, \
                         gf_1_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_78[k] = ab_x[k] * ff_1_53[k]
                      + gf_1_53[k];

            t_79[k] = ab_x[k] * ff_1_54[k]
                      + gf_1_54[k];

            t_80[k] = ab_x[k] * ff_1_55[k]
                      + gf_1_55[k];

            t_81[k] = ab_x[k] * ff_1_56[k]
                      + gf_1_56[k];

            t_82[k] = ab_x[k] * ff_1_57[k]
                      + gf_1_57[k];
        }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ff_1_56, ff_1_57, ff_1_58, \
                         ff_1_59, gf_1_58, gf_1_59, gf_1_86, gf_1_87, \
                         gf_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_83[k] = ab_x[k] * ff_1_58[k]
                      + gf_1_58[k];

            t_84[k] = ab_x[k] * ff_1_59[k]
                      + gf_1_59[k];

            t_85[k] = ab_y[k] * ff_1_56[k]
                      + gf_1_86[k];

            t_86[k] = ab_y[k] * ff_1_57[k]
                      + gf_1_87[k];

            t_87[k] = ab_y[k] * ff_1_58[k]
                      + gf_1_88[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, ab_x, ab_y, ab_z, ff_1_59, ff_1_60, ff_1_61, \
                         ff_0_59, gf_1_60, gf_1_61, gf_1_89, gf_1_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * ff_1_59[k]
                      + gf_1_89[k];

            t_89[k] = ab_z[k] * ff_1_59[k]
                      + ff_0_59[k]
                      + gf_1_99[k];

            t_90[k] = ab_x[k] * ff_1_60[k]
                      + gf_1_60[k];

            t_91[k] = ab_x[k] * ff_1_61[k]
                      + gf_1_61[k];
        }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, ab_x, ff_1_62, ff_1_63, ff_1_64, \
                         ff_1_65, ff_1_66, gf_1_62, gf_1_63, gf_1_64, gf_1_65, \
                         gf_1_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_92[k] = ab_x[k] * ff_1_62[k]
                      + gf_1_62[k];

            t_93[k] = ab_x[k] * ff_1_63[k]
                      + gf_1_63[k];

            t_94[k] = ab_x[k] * ff_1_64[k]
                      + gf_1_64[k];

            t_95[k] = ab_x[k] * ff_1_65[k]
                      + gf_1_65[k];

            t_96[k] = ab_x[k] * ff_1_66[k]
                      + gf_1_66[k];
        }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, ab_x, ab_y, ff_1_66, ff_1_67, \
                         ff_1_68, ff_1_69, gf_1_67, gf_1_68, gf_1_69, gf_1_106, \
                         gf_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_97[k] = ab_x[k] * ff_1_67[k]
                      + gf_1_67[k];

            t_98[k] = ab_x[k] * ff_1_68[k]
                      + gf_1_68[k];

            t_99[k] = ab_x[k] * ff_1_69[k]
                      + gf_1_69[k];

            t_100[k] = ab_y[k] * ff_1_66[k]
                       + gf_1_106[k];

            t_101[k] = ab_y[k] * ff_1_67[k]
                       + gf_1_107[k];
        }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, ab_x, ab_y, ab_z, ff_1_68, ff_1_69, \
                         ff_1_70, ff_0_69, gf_1_70, gf_1_108, gf_1_109, \
                         gf_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_102[k] = ab_y[k] * ff_1_68[k]
                       + gf_1_108[k];

            t_103[k] = ab_y[k] * ff_1_69[k]
                       + gf_1_109[k];

            t_104[k] = ab_z[k] * ff_1_69[k]
                       + ff_0_69[k]
                       + gf_1_119[k];

            t_105[k] = ab_x[k] * ff_1_70[k]
                       + gf_1_70[k];
        }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, ab_x, ff_1_71, ff_1_72, ff_1_73, \
                         ff_1_74, ff_1_75, gf_1_71, gf_1_72, gf_1_73, gf_1_74, \
                         gf_1_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_106[k] = ab_x[k] * ff_1_71[k]
                       + gf_1_71[k];

            t_107[k] = ab_x[k] * ff_1_72[k]
                       + gf_1_72[k];

            t_108[k] = ab_x[k] * ff_1_73[k]
                       + gf_1_73[k];

            t_109[k] = ab_x[k] * ff_1_74[k]
                       + gf_1_74[k];

            t_110[k] = ab_x[k] * ff_1_75[k]
                       + gf_1_75[k];
        }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, ab_x, ab_y, ff_1_76, ff_1_77, \
                         ff_1_78, ff_1_79, gf_1_76, gf_1_77, gf_1_78, gf_1_79, \
                         gf_1_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_111[k] = ab_x[k] * ff_1_76[k]
                       + gf_1_76[k];

            t_112[k] = ab_x[k] * ff_1_77[k]
                       + gf_1_77[k];

            t_113[k] = ab_x[k] * ff_1_78[k]
                       + gf_1_78[k];

            t_114[k] = ab_x[k] * ff_1_79[k]
                       + gf_1_79[k];

            t_115[k] = ab_y[k] * ff_1_76[k]
                       + gf_1_116[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, ab_y, ab_z, ff_1_77, ff_1_78, ff_1_79, \
                         ff_0_79, gf_1_117, gf_1_118, gf_1_119, \
                         gf_1_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_y[k] * ff_1_77[k]
                       + gf_1_117[k];

            t_117[k] = ab_y[k] * ff_1_78[k]
                       + gf_1_118[k];

            t_118[k] = ab_y[k] * ff_1_79[k]
                       + gf_1_119[k];

            t_119[k] = ab_z[k] * ff_1_79[k]
                       + ff_0_79[k]
                       + gf_1_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, ff_1_80, ff_1_81, ff_1_82, \
                         ff_1_83, ff_1_84, gf_1_80, gf_1_81, gf_1_82, gf_1_83, \
                         gf_1_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * ff_1_80[k]
                       + gf_1_80[k];

            t_121[k] = ab_x[k] * ff_1_81[k]
                       + gf_1_81[k];

            t_122[k] = ab_x[k] * ff_1_82[k]
                       + gf_1_82[k];

            t_123[k] = ab_x[k] * ff_1_83[k]
                       + gf_1_83[k];

            t_124[k] = ab_x[k] * ff_1_84[k]
                       + gf_1_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ff_1_85, ff_1_86, ff_1_87, \
                         ff_1_88, ff_1_89, gf_1_85, gf_1_86, gf_1_87, gf_1_88, \
                         gf_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * ff_1_85[k]
                       + gf_1_85[k];

            t_126[k] = ab_x[k] * ff_1_86[k]
                       + gf_1_86[k];

            t_127[k] = ab_x[k] * ff_1_87[k]
                       + gf_1_87[k];

            t_128[k] = ab_x[k] * ff_1_88[k]
                       + gf_1_88[k];

            t_129[k] = ab_x[k] * ff_1_89[k]
                       + gf_1_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, ab_y, ff_1_86, ff_1_87, ff_1_88, ff_1_89, \
                         gf_1_126, gf_1_127, gf_1_128, gf_1_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * ff_1_86[k]
                       + gf_1_126[k];

            t_131[k] = ab_y[k] * ff_1_87[k]
                       + gf_1_127[k];

            t_132[k] = ab_y[k] * ff_1_88[k]
                       + gf_1_128[k];

            t_133[k] = ab_y[k] * ff_1_89[k]
                       + gf_1_129[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, ab_x, ab_z, ff_1_89, ff_1_90, ff_1_91, \
                         ff_1_92, ff_0_89, gf_1_90, gf_1_91, gf_1_92, \
                         gf_1_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = ab_z[k] * ff_1_89[k]
                       + ff_0_89[k]
                       + gf_1_139[k];

            t_135[k] = ab_x[k] * ff_1_90[k]
                       + gf_1_90[k];

            t_136[k] = ab_x[k] * ff_1_91[k]
                       + gf_1_91[k];

            t_137[k] = ab_x[k] * ff_1_92[k]
                       + gf_1_92[k];
        }
    }
}

static auto
compute_hrr_geom_100z_fg_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t ff_1,
                                             const size_t ff_0, const size_t gf_1,
                                             const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ff_1_93 = buffer.data(ff_1 + 93 * ncomps + c);
        const auto *ff_1_94 = buffer.data(ff_1 + 94 * ncomps + c);
        const auto *ff_1_95 = buffer.data(ff_1 + 95 * ncomps + c);
        const auto *ff_1_96 = buffer.data(ff_1 + 96 * ncomps + c);
        const auto *ff_1_97 = buffer.data(ff_1 + 97 * ncomps + c);
        const auto *ff_1_98 = buffer.data(ff_1 + 98 * ncomps + c);
        const auto *ff_1_99 = buffer.data(ff_1 + 99 * ncomps + c);

        const auto *ff_0_99 = buffer.data(ff_0 + 99 * ncomps + c);

        const auto *gf_1_93 = buffer.data(gf_1 + 93 * ncomps + c);
        const auto *gf_1_94 = buffer.data(gf_1 + 94 * ncomps + c);
        const auto *gf_1_95 = buffer.data(gf_1 + 95 * ncomps + c);
        const auto *gf_1_96 = buffer.data(gf_1 + 96 * ncomps + c);
        const auto *gf_1_97 = buffer.data(gf_1 + 97 * ncomps + c);
        const auto *gf_1_98 = buffer.data(gf_1 + 98 * ncomps + c);
        const auto *gf_1_99 = buffer.data(gf_1 + 99 * ncomps + c);
        const auto *gf_1_136 = buffer.data(gf_1 + 136 * ncomps + c);
        const auto *gf_1_137 = buffer.data(gf_1 + 137 * ncomps + c);
        const auto *gf_1_138 = buffer.data(gf_1 + 138 * ncomps + c);
        const auto *gf_1_139 = buffer.data(gf_1 + 139 * ncomps + c);
        const auto *gf_1_149 = buffer.data(gf_1 + 149 * ncomps + c);

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, ab_x, ff_1_93, ff_1_94, ff_1_95, \
                         ff_1_96, ff_1_97, gf_1_93, gf_1_94, gf_1_95, gf_1_96, \
                         gf_1_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_138[k] = ab_x[k] * ff_1_93[k]
                       + gf_1_93[k];

            t_139[k] = ab_x[k] * ff_1_94[k]
                       + gf_1_94[k];

            t_140[k] = ab_x[k] * ff_1_95[k]
                       + gf_1_95[k];

            t_141[k] = ab_x[k] * ff_1_96[k]
                       + gf_1_96[k];

            t_142[k] = ab_x[k] * ff_1_97[k]
                       + gf_1_97[k];
        }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, ab_x, ab_y, ff_1_96, ff_1_97, \
                         ff_1_98, ff_1_99, gf_1_98, gf_1_99, gf_1_136, gf_1_137, \
                         gf_1_138 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_143[k] = ab_x[k] * ff_1_98[k]
                       + gf_1_98[k];

            t_144[k] = ab_x[k] * ff_1_99[k]
                       + gf_1_99[k];

            t_145[k] = ab_y[k] * ff_1_96[k]
                       + gf_1_136[k];

            t_146[k] = ab_y[k] * ff_1_97[k]
                       + gf_1_137[k];

            t_147[k] = ab_y[k] * ff_1_98[k]
                       + gf_1_138[k];
        }

#pragma omp simd aligned(t_148, t_149, ab_y, ab_z, ff_1_99, ff_0_99, gf_1_139, \
                         gf_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_148[k] = ab_y[k] * ff_1_99[k]
                       + gf_1_139[k];

            t_149[k] = ab_z[k] * ff_1_99[k]
                       + ff_0_99[k]
                       + gf_1_149[k];
        }
    }
}

auto
compute_hrr_geom_100z_fg_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t ff_1, const size_t ff_0,
                                      const size_t gf_1, const size_t ncomps,
                                      const size_t nmax) -> void
{
    compute_hrr_geom_100z_fg_out_of_first_piece0(buffer, coordinates, target, ff_1, ff_0, gf_1,
                                                 ncomps, nmax);

    compute_hrr_geom_100z_fg_out_of_first_piece1(buffer, coordinates, target, ff_1, ff_0, gf_1,
                                                 ncomps, nmax);
}

}  // namespace simdtrf
