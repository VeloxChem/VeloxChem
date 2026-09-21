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


#include "SimdTransferGeom100ZGF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_100z_gf_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t gd_1,
                                             const size_t gd_0, const size_t hd_1,
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

        const auto *gd_1_0 = buffer.data(gd_1 + 0 * ncomps + c);
        const auto *gd_1_1 = buffer.data(gd_1 + 1 * ncomps + c);
        const auto *gd_1_2 = buffer.data(gd_1 + 2 * ncomps + c);
        const auto *gd_1_3 = buffer.data(gd_1 + 3 * ncomps + c);
        const auto *gd_1_4 = buffer.data(gd_1 + 4 * ncomps + c);
        const auto *gd_1_5 = buffer.data(gd_1 + 5 * ncomps + c);
        const auto *gd_1_6 = buffer.data(gd_1 + 6 * ncomps + c);
        const auto *gd_1_7 = buffer.data(gd_1 + 7 * ncomps + c);
        const auto *gd_1_8 = buffer.data(gd_1 + 8 * ncomps + c);
        const auto *gd_1_9 = buffer.data(gd_1 + 9 * ncomps + c);
        const auto *gd_1_10 = buffer.data(gd_1 + 10 * ncomps + c);
        const auto *gd_1_11 = buffer.data(gd_1 + 11 * ncomps + c);
        const auto *gd_1_12 = buffer.data(gd_1 + 12 * ncomps + c);
        const auto *gd_1_13 = buffer.data(gd_1 + 13 * ncomps + c);
        const auto *gd_1_14 = buffer.data(gd_1 + 14 * ncomps + c);
        const auto *gd_1_15 = buffer.data(gd_1 + 15 * ncomps + c);
        const auto *gd_1_16 = buffer.data(gd_1 + 16 * ncomps + c);
        const auto *gd_1_17 = buffer.data(gd_1 + 17 * ncomps + c);
        const auto *gd_1_18 = buffer.data(gd_1 + 18 * ncomps + c);
        const auto *gd_1_19 = buffer.data(gd_1 + 19 * ncomps + c);
        const auto *gd_1_20 = buffer.data(gd_1 + 20 * ncomps + c);
        const auto *gd_1_21 = buffer.data(gd_1 + 21 * ncomps + c);
        const auto *gd_1_22 = buffer.data(gd_1 + 22 * ncomps + c);
        const auto *gd_1_23 = buffer.data(gd_1 + 23 * ncomps + c);
        const auto *gd_1_24 = buffer.data(gd_1 + 24 * ncomps + c);
        const auto *gd_1_25 = buffer.data(gd_1 + 25 * ncomps + c);
        const auto *gd_1_26 = buffer.data(gd_1 + 26 * ncomps + c);
        const auto *gd_1_27 = buffer.data(gd_1 + 27 * ncomps + c);
        const auto *gd_1_28 = buffer.data(gd_1 + 28 * ncomps + c);
        const auto *gd_1_29 = buffer.data(gd_1 + 29 * ncomps + c);
        const auto *gd_1_30 = buffer.data(gd_1 + 30 * ncomps + c);
        const auto *gd_1_31 = buffer.data(gd_1 + 31 * ncomps + c);
        const auto *gd_1_32 = buffer.data(gd_1 + 32 * ncomps + c);
        const auto *gd_1_33 = buffer.data(gd_1 + 33 * ncomps + c);
        const auto *gd_1_34 = buffer.data(gd_1 + 34 * ncomps + c);
        const auto *gd_1_35 = buffer.data(gd_1 + 35 * ncomps + c);
        const auto *gd_1_36 = buffer.data(gd_1 + 36 * ncomps + c);
        const auto *gd_1_37 = buffer.data(gd_1 + 37 * ncomps + c);
        const auto *gd_1_38 = buffer.data(gd_1 + 38 * ncomps + c);
        const auto *gd_1_39 = buffer.data(gd_1 + 39 * ncomps + c);
        const auto *gd_1_40 = buffer.data(gd_1 + 40 * ncomps + c);
        const auto *gd_1_41 = buffer.data(gd_1 + 41 * ncomps + c);
        const auto *gd_1_42 = buffer.data(gd_1 + 42 * ncomps + c);
        const auto *gd_1_43 = buffer.data(gd_1 + 43 * ncomps + c);
        const auto *gd_1_44 = buffer.data(gd_1 + 44 * ncomps + c);
        const auto *gd_1_45 = buffer.data(gd_1 + 45 * ncomps + c);
        const auto *gd_1_46 = buffer.data(gd_1 + 46 * ncomps + c);
        const auto *gd_1_47 = buffer.data(gd_1 + 47 * ncomps + c);
        const auto *gd_1_48 = buffer.data(gd_1 + 48 * ncomps + c);
        const auto *gd_1_49 = buffer.data(gd_1 + 49 * ncomps + c);
        const auto *gd_1_50 = buffer.data(gd_1 + 50 * ncomps + c);
        const auto *gd_1_51 = buffer.data(gd_1 + 51 * ncomps + c);
        const auto *gd_1_52 = buffer.data(gd_1 + 52 * ncomps + c);
        const auto *gd_1_53 = buffer.data(gd_1 + 53 * ncomps + c);
        const auto *gd_1_54 = buffer.data(gd_1 + 54 * ncomps + c);
        const auto *gd_1_55 = buffer.data(gd_1 + 55 * ncomps + c);
        const auto *gd_1_56 = buffer.data(gd_1 + 56 * ncomps + c);
        const auto *gd_1_57 = buffer.data(gd_1 + 57 * ncomps + c);
        const auto *gd_1_58 = buffer.data(gd_1 + 58 * ncomps + c);
        const auto *gd_1_59 = buffer.data(gd_1 + 59 * ncomps + c);
        const auto *gd_1_60 = buffer.data(gd_1 + 60 * ncomps + c);
        const auto *gd_1_61 = buffer.data(gd_1 + 61 * ncomps + c);
        const auto *gd_1_62 = buffer.data(gd_1 + 62 * ncomps + c);
        const auto *gd_1_63 = buffer.data(gd_1 + 63 * ncomps + c);
        const auto *gd_1_64 = buffer.data(gd_1 + 64 * ncomps + c);
        const auto *gd_1_65 = buffer.data(gd_1 + 65 * ncomps + c);
        const auto *gd_1_66 = buffer.data(gd_1 + 66 * ncomps + c);
        const auto *gd_1_67 = buffer.data(gd_1 + 67 * ncomps + c);
        const auto *gd_1_68 = buffer.data(gd_1 + 68 * ncomps + c);
        const auto *gd_1_69 = buffer.data(gd_1 + 69 * ncomps + c);
        const auto *gd_1_70 = buffer.data(gd_1 + 70 * ncomps + c);
        const auto *gd_1_71 = buffer.data(gd_1 + 71 * ncomps + c);
        const auto *gd_1_72 = buffer.data(gd_1 + 72 * ncomps + c);
        const auto *gd_1_73 = buffer.data(gd_1 + 73 * ncomps + c);
        const auto *gd_1_74 = buffer.data(gd_1 + 74 * ncomps + c);
        const auto *gd_1_75 = buffer.data(gd_1 + 75 * ncomps + c);
        const auto *gd_1_76 = buffer.data(gd_1 + 76 * ncomps + c);
        const auto *gd_1_77 = buffer.data(gd_1 + 77 * ncomps + c);
        const auto *gd_1_78 = buffer.data(gd_1 + 78 * ncomps + c);
        const auto *gd_1_79 = buffer.data(gd_1 + 79 * ncomps + c);
        const auto *gd_1_80 = buffer.data(gd_1 + 80 * ncomps + c);
        const auto *gd_1_81 = buffer.data(gd_1 + 81 * ncomps + c);
        const auto *gd_1_82 = buffer.data(gd_1 + 82 * ncomps + c);
        const auto *gd_1_83 = buffer.data(gd_1 + 83 * ncomps + c);

        const auto *gd_0_5 = buffer.data(gd_0 + 5 * ncomps + c);
        const auto *gd_0_11 = buffer.data(gd_0 + 11 * ncomps + c);
        const auto *gd_0_17 = buffer.data(gd_0 + 17 * ncomps + c);
        const auto *gd_0_23 = buffer.data(gd_0 + 23 * ncomps + c);
        const auto *gd_0_29 = buffer.data(gd_0 + 29 * ncomps + c);
        const auto *gd_0_35 = buffer.data(gd_0 + 35 * ncomps + c);
        const auto *gd_0_41 = buffer.data(gd_0 + 41 * ncomps + c);
        const auto *gd_0_47 = buffer.data(gd_0 + 47 * ncomps + c);
        const auto *gd_0_53 = buffer.data(gd_0 + 53 * ncomps + c);
        const auto *gd_0_59 = buffer.data(gd_0 + 59 * ncomps + c);
        const auto *gd_0_65 = buffer.data(gd_0 + 65 * ncomps + c);
        const auto *gd_0_71 = buffer.data(gd_0 + 71 * ncomps + c);
        const auto *gd_0_77 = buffer.data(gd_0 + 77 * ncomps + c);

        const auto *hd_1_0 = buffer.data(hd_1 + 0 * ncomps + c);
        const auto *hd_1_1 = buffer.data(hd_1 + 1 * ncomps + c);
        const auto *hd_1_2 = buffer.data(hd_1 + 2 * ncomps + c);
        const auto *hd_1_3 = buffer.data(hd_1 + 3 * ncomps + c);
        const auto *hd_1_4 = buffer.data(hd_1 + 4 * ncomps + c);
        const auto *hd_1_5 = buffer.data(hd_1 + 5 * ncomps + c);
        const auto *hd_1_6 = buffer.data(hd_1 + 6 * ncomps + c);
        const auto *hd_1_7 = buffer.data(hd_1 + 7 * ncomps + c);
        const auto *hd_1_8 = buffer.data(hd_1 + 8 * ncomps + c);
        const auto *hd_1_9 = buffer.data(hd_1 + 9 * ncomps + c);
        const auto *hd_1_10 = buffer.data(hd_1 + 10 * ncomps + c);
        const auto *hd_1_11 = buffer.data(hd_1 + 11 * ncomps + c);
        const auto *hd_1_12 = buffer.data(hd_1 + 12 * ncomps + c);
        const auto *hd_1_13 = buffer.data(hd_1 + 13 * ncomps + c);
        const auto *hd_1_14 = buffer.data(hd_1 + 14 * ncomps + c);
        const auto *hd_1_15 = buffer.data(hd_1 + 15 * ncomps + c);
        const auto *hd_1_16 = buffer.data(hd_1 + 16 * ncomps + c);
        const auto *hd_1_17 = buffer.data(hd_1 + 17 * ncomps + c);
        const auto *hd_1_18 = buffer.data(hd_1 + 18 * ncomps + c);
        const auto *hd_1_19 = buffer.data(hd_1 + 19 * ncomps + c);
        const auto *hd_1_20 = buffer.data(hd_1 + 20 * ncomps + c);
        const auto *hd_1_21 = buffer.data(hd_1 + 21 * ncomps + c);
        const auto *hd_1_22 = buffer.data(hd_1 + 22 * ncomps + c);
        const auto *hd_1_23 = buffer.data(hd_1 + 23 * ncomps + c);
        const auto *hd_1_24 = buffer.data(hd_1 + 24 * ncomps + c);
        const auto *hd_1_25 = buffer.data(hd_1 + 25 * ncomps + c);
        const auto *hd_1_26 = buffer.data(hd_1 + 26 * ncomps + c);
        const auto *hd_1_27 = buffer.data(hd_1 + 27 * ncomps + c);
        const auto *hd_1_28 = buffer.data(hd_1 + 28 * ncomps + c);
        const auto *hd_1_29 = buffer.data(hd_1 + 29 * ncomps + c);
        const auto *hd_1_30 = buffer.data(hd_1 + 30 * ncomps + c);
        const auto *hd_1_31 = buffer.data(hd_1 + 31 * ncomps + c);
        const auto *hd_1_32 = buffer.data(hd_1 + 32 * ncomps + c);
        const auto *hd_1_33 = buffer.data(hd_1 + 33 * ncomps + c);
        const auto *hd_1_34 = buffer.data(hd_1 + 34 * ncomps + c);
        const auto *hd_1_35 = buffer.data(hd_1 + 35 * ncomps + c);
        const auto *hd_1_36 = buffer.data(hd_1 + 36 * ncomps + c);
        const auto *hd_1_37 = buffer.data(hd_1 + 37 * ncomps + c);
        const auto *hd_1_38 = buffer.data(hd_1 + 38 * ncomps + c);
        const auto *hd_1_39 = buffer.data(hd_1 + 39 * ncomps + c);
        const auto *hd_1_40 = buffer.data(hd_1 + 40 * ncomps + c);
        const auto *hd_1_41 = buffer.data(hd_1 + 41 * ncomps + c);
        const auto *hd_1_42 = buffer.data(hd_1 + 42 * ncomps + c);
        const auto *hd_1_43 = buffer.data(hd_1 + 43 * ncomps + c);
        const auto *hd_1_44 = buffer.data(hd_1 + 44 * ncomps + c);
        const auto *hd_1_45 = buffer.data(hd_1 + 45 * ncomps + c);
        const auto *hd_1_46 = buffer.data(hd_1 + 46 * ncomps + c);
        const auto *hd_1_47 = buffer.data(hd_1 + 47 * ncomps + c);
        const auto *hd_1_48 = buffer.data(hd_1 + 48 * ncomps + c);
        const auto *hd_1_49 = buffer.data(hd_1 + 49 * ncomps + c);
        const auto *hd_1_50 = buffer.data(hd_1 + 50 * ncomps + c);
        const auto *hd_1_51 = buffer.data(hd_1 + 51 * ncomps + c);
        const auto *hd_1_52 = buffer.data(hd_1 + 52 * ncomps + c);
        const auto *hd_1_53 = buffer.data(hd_1 + 53 * ncomps + c);
        const auto *hd_1_54 = buffer.data(hd_1 + 54 * ncomps + c);
        const auto *hd_1_55 = buffer.data(hd_1 + 55 * ncomps + c);
        const auto *hd_1_56 = buffer.data(hd_1 + 56 * ncomps + c);
        const auto *hd_1_57 = buffer.data(hd_1 + 57 * ncomps + c);
        const auto *hd_1_58 = buffer.data(hd_1 + 58 * ncomps + c);
        const auto *hd_1_59 = buffer.data(hd_1 + 59 * ncomps + c);
        const auto *hd_1_60 = buffer.data(hd_1 + 60 * ncomps + c);
        const auto *hd_1_61 = buffer.data(hd_1 + 61 * ncomps + c);
        const auto *hd_1_62 = buffer.data(hd_1 + 62 * ncomps + c);
        const auto *hd_1_63 = buffer.data(hd_1 + 63 * ncomps + c);
        const auto *hd_1_64 = buffer.data(hd_1 + 64 * ncomps + c);
        const auto *hd_1_65 = buffer.data(hd_1 + 65 * ncomps + c);
        const auto *hd_1_66 = buffer.data(hd_1 + 66 * ncomps + c);
        const auto *hd_1_67 = buffer.data(hd_1 + 67 * ncomps + c);
        const auto *hd_1_68 = buffer.data(hd_1 + 68 * ncomps + c);
        const auto *hd_1_69 = buffer.data(hd_1 + 69 * ncomps + c);
        const auto *hd_1_70 = buffer.data(hd_1 + 70 * ncomps + c);
        const auto *hd_1_71 = buffer.data(hd_1 + 71 * ncomps + c);
        const auto *hd_1_72 = buffer.data(hd_1 + 72 * ncomps + c);
        const auto *hd_1_73 = buffer.data(hd_1 + 73 * ncomps + c);
        const auto *hd_1_74 = buffer.data(hd_1 + 74 * ncomps + c);
        const auto *hd_1_75 = buffer.data(hd_1 + 75 * ncomps + c);
        const auto *hd_1_76 = buffer.data(hd_1 + 76 * ncomps + c);
        const auto *hd_1_77 = buffer.data(hd_1 + 77 * ncomps + c);
        const auto *hd_1_78 = buffer.data(hd_1 + 78 * ncomps + c);
        const auto *hd_1_79 = buffer.data(hd_1 + 79 * ncomps + c);
        const auto *hd_1_80 = buffer.data(hd_1 + 80 * ncomps + c);
        const auto *hd_1_81 = buffer.data(hd_1 + 81 * ncomps + c);
        const auto *hd_1_82 = buffer.data(hd_1 + 82 * ncomps + c);
        const auto *hd_1_83 = buffer.data(hd_1 + 83 * ncomps + c);
        const auto *hd_1_89 = buffer.data(hd_1 + 89 * ncomps + c);
        const auto *hd_1_93 = buffer.data(hd_1 + 93 * ncomps + c);
        const auto *hd_1_94 = buffer.data(hd_1 + 94 * ncomps + c);
        const auto *hd_1_95 = buffer.data(hd_1 + 95 * ncomps + c);
        const auto *hd_1_99 = buffer.data(hd_1 + 99 * ncomps + c);
        const auto *hd_1_100 = buffer.data(hd_1 + 100 * ncomps + c);
        const auto *hd_1_101 = buffer.data(hd_1 + 101 * ncomps + c);
        const auto *hd_1_105 = buffer.data(hd_1 + 105 * ncomps + c);
        const auto *hd_1_106 = buffer.data(hd_1 + 106 * ncomps + c);
        const auto *hd_1_107 = buffer.data(hd_1 + 107 * ncomps + c);
        const auto *hd_1_111 = buffer.data(hd_1 + 111 * ncomps + c);
        const auto *hd_1_112 = buffer.data(hd_1 + 112 * ncomps + c);
        const auto *hd_1_113 = buffer.data(hd_1 + 113 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gd_1_0, gd_1_1, gd_1_2, gd_1_3, \
                         gd_1_4, hd_1_0, hd_1_1, hd_1_2, hd_1_3, \
                         hd_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gd_1_0[k]
                     + hd_1_0[k];

            t_1[k] = ab_x[k] * gd_1_1[k]
                     + hd_1_1[k];

            t_2[k] = ab_x[k] * gd_1_2[k]
                     + hd_1_2[k];

            t_3[k] = ab_x[k] * gd_1_3[k]
                     + hd_1_3[k];

            t_4[k] = ab_x[k] * gd_1_4[k]
                     + hd_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_y, gd_1_3, gd_1_4, gd_1_5, hd_1_5, \
                         hd_1_9, hd_1_10, hd_1_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * gd_1_5[k]
                     + hd_1_5[k];

            t_6[k] = ab_y[k] * gd_1_3[k]
                     + hd_1_9[k];

            t_7[k] = ab_y[k] * gd_1_4[k]
                     + hd_1_10[k];

            t_8[k] = ab_y[k] * gd_1_5[k]
                     + hd_1_11[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_z, gd_1_5, gd_1_6, gd_1_7, gd_1_8, \
                         gd_0_5, hd_1_6, hd_1_7, hd_1_8, hd_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_z[k] * gd_1_5[k]
                     + gd_0_5[k]
                     + hd_1_17[k];

            t_10[k] = ab_x[k] * gd_1_6[k]
                      + hd_1_6[k];

            t_11[k] = ab_x[k] * gd_1_7[k]
                      + hd_1_7[k];

            t_12[k] = ab_x[k] * gd_1_8[k]
                      + hd_1_8[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, gd_1_9, gd_1_10, gd_1_11, \
                         hd_1_9, hd_1_10, hd_1_11, hd_1_21, hd_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * gd_1_9[k]
                      + hd_1_9[k];

            t_14[k] = ab_x[k] * gd_1_10[k]
                      + hd_1_10[k];

            t_15[k] = ab_x[k] * gd_1_11[k]
                      + hd_1_11[k];

            t_16[k] = ab_y[k] * gd_1_9[k]
                      + hd_1_21[k];

            t_17[k] = ab_y[k] * gd_1_10[k]
                      + hd_1_22[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, gd_1_11, gd_1_12, gd_1_13, \
                         gd_0_11, hd_1_12, hd_1_13, hd_1_23, hd_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_y[k] * gd_1_11[k]
                      + hd_1_23[k];

            t_19[k] = ab_z[k] * gd_1_11[k]
                      + gd_0_11[k]
                      + hd_1_29[k];

            t_20[k] = ab_x[k] * gd_1_12[k]
                      + hd_1_12[k];

            t_21[k] = ab_x[k] * gd_1_13[k]
                      + hd_1_13[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, gd_1_14, gd_1_15, gd_1_16, \
                         gd_1_17, hd_1_14, hd_1_15, hd_1_16, hd_1_17, \
                         hd_1_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_x[k] * gd_1_14[k]
                      + hd_1_14[k];

            t_23[k] = ab_x[k] * gd_1_15[k]
                      + hd_1_15[k];

            t_24[k] = ab_x[k] * gd_1_16[k]
                      + hd_1_16[k];

            t_25[k] = ab_x[k] * gd_1_17[k]
                      + hd_1_17[k];

            t_26[k] = ab_y[k] * gd_1_15[k]
                      + hd_1_27[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, gd_1_16, gd_1_17, gd_1_18, \
                         gd_0_17, hd_1_18, hd_1_28, hd_1_29, hd_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * gd_1_16[k]
                      + hd_1_28[k];

            t_28[k] = ab_y[k] * gd_1_17[k]
                      + hd_1_29[k];

            t_29[k] = ab_z[k] * gd_1_17[k]
                      + gd_0_17[k]
                      + hd_1_35[k];

            t_30[k] = ab_x[k] * gd_1_18[k]
                      + hd_1_18[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, gd_1_19, gd_1_20, gd_1_21, \
                         gd_1_22, gd_1_23, hd_1_19, hd_1_20, hd_1_21, hd_1_22, \
                         hd_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * gd_1_19[k]
                      + hd_1_19[k];

            t_32[k] = ab_x[k] * gd_1_20[k]
                      + hd_1_20[k];

            t_33[k] = ab_x[k] * gd_1_21[k]
                      + hd_1_21[k];

            t_34[k] = ab_x[k] * gd_1_22[k]
                      + hd_1_22[k];

            t_35[k] = ab_x[k] * gd_1_23[k]
                      + hd_1_23[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, ab_y, ab_z, gd_1_21, gd_1_22, gd_1_23, \
                         gd_0_23, hd_1_39, hd_1_40, hd_1_41, hd_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_y[k] * gd_1_21[k]
                      + hd_1_39[k];

            t_37[k] = ab_y[k] * gd_1_22[k]
                      + hd_1_40[k];

            t_38[k] = ab_y[k] * gd_1_23[k]
                      + hd_1_41[k];

            t_39[k] = ab_z[k] * gd_1_23[k]
                      + gd_0_23[k]
                      + hd_1_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gd_1_24, gd_1_25, gd_1_26, \
                         gd_1_27, gd_1_28, hd_1_24, hd_1_25, hd_1_26, hd_1_27, \
                         hd_1_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * gd_1_24[k]
                      + hd_1_24[k];

            t_41[k] = ab_x[k] * gd_1_25[k]
                      + hd_1_25[k];

            t_42[k] = ab_x[k] * gd_1_26[k]
                      + hd_1_26[k];

            t_43[k] = ab_x[k] * gd_1_27[k]
                      + hd_1_27[k];

            t_44[k] = ab_x[k] * gd_1_28[k]
                      + hd_1_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, gd_1_27, gd_1_28, gd_1_29, \
                         hd_1_29, hd_1_45, hd_1_46, hd_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * gd_1_29[k]
                      + hd_1_29[k];

            t_46[k] = ab_y[k] * gd_1_27[k]
                      + hd_1_45[k];

            t_47[k] = ab_y[k] * gd_1_28[k]
                      + hd_1_46[k];

            t_48[k] = ab_y[k] * gd_1_29[k]
                      + hd_1_47[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, ab_x, ab_z, gd_1_29, gd_1_30, gd_1_31, \
                         gd_1_32, gd_0_29, hd_1_30, hd_1_31, hd_1_32, \
                         hd_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_z[k] * gd_1_29[k]
                      + gd_0_29[k]
                      + hd_1_53[k];

            t_50[k] = ab_x[k] * gd_1_30[k]
                      + hd_1_30[k];

            t_51[k] = ab_x[k] * gd_1_31[k]
                      + hd_1_31[k];

            t_52[k] = ab_x[k] * gd_1_32[k]
                      + hd_1_32[k];
        }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, ab_x, ab_y, gd_1_33, gd_1_34, gd_1_35, \
                         hd_1_33, hd_1_34, hd_1_35, hd_1_51, hd_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_53[k] = ab_x[k] * gd_1_33[k]
                      + hd_1_33[k];

            t_54[k] = ab_x[k] * gd_1_34[k]
                      + hd_1_34[k];

            t_55[k] = ab_x[k] * gd_1_35[k]
                      + hd_1_35[k];

            t_56[k] = ab_y[k] * gd_1_33[k]
                      + hd_1_51[k];

            t_57[k] = ab_y[k] * gd_1_34[k]
                      + hd_1_52[k];
        }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, ab_x, ab_y, ab_z, gd_1_35, gd_1_36, gd_1_37, \
                         gd_0_35, hd_1_36, hd_1_37, hd_1_53, hd_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_58[k] = ab_y[k] * gd_1_35[k]
                      + hd_1_53[k];

            t_59[k] = ab_z[k] * gd_1_35[k]
                      + gd_0_35[k]
                      + hd_1_59[k];

            t_60[k] = ab_x[k] * gd_1_36[k]
                      + hd_1_36[k];

            t_61[k] = ab_x[k] * gd_1_37[k]
                      + hd_1_37[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, ab_x, ab_y, gd_1_38, gd_1_39, gd_1_40, \
                         gd_1_41, hd_1_38, hd_1_39, hd_1_40, hd_1_41, \
                         hd_1_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_x[k] * gd_1_38[k]
                      + hd_1_38[k];

            t_63[k] = ab_x[k] * gd_1_39[k]
                      + hd_1_39[k];

            t_64[k] = ab_x[k] * gd_1_40[k]
                      + hd_1_40[k];

            t_65[k] = ab_x[k] * gd_1_41[k]
                      + hd_1_41[k];

            t_66[k] = ab_y[k] * gd_1_39[k]
                      + hd_1_63[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, gd_1_40, gd_1_41, gd_1_42, \
                         gd_0_41, hd_1_42, hd_1_64, hd_1_65, hd_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = ab_y[k] * gd_1_40[k]
                      + hd_1_64[k];

            t_68[k] = ab_y[k] * gd_1_41[k]
                      + hd_1_65[k];

            t_69[k] = ab_z[k] * gd_1_41[k]
                      + gd_0_41[k]
                      + hd_1_71[k];

            t_70[k] = ab_x[k] * gd_1_42[k]
                      + hd_1_42[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ab_x, gd_1_43, gd_1_44, gd_1_45, \
                         gd_1_46, gd_1_47, hd_1_43, hd_1_44, hd_1_45, hd_1_46, \
                         hd_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_x[k] * gd_1_43[k]
                      + hd_1_43[k];

            t_72[k] = ab_x[k] * gd_1_44[k]
                      + hd_1_44[k];

            t_73[k] = ab_x[k] * gd_1_45[k]
                      + hd_1_45[k];

            t_74[k] = ab_x[k] * gd_1_46[k]
                      + hd_1_46[k];

            t_75[k] = ab_x[k] * gd_1_47[k]
                      + hd_1_47[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_y, ab_z, gd_1_45, gd_1_46, gd_1_47, \
                         gd_0_47, hd_1_69, hd_1_70, hd_1_71, hd_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_y[k] * gd_1_45[k]
                      + hd_1_69[k];

            t_77[k] = ab_y[k] * gd_1_46[k]
                      + hd_1_70[k];

            t_78[k] = ab_y[k] * gd_1_47[k]
                      + hd_1_71[k];

            t_79[k] = ab_z[k] * gd_1_47[k]
                      + gd_0_47[k]
                      + hd_1_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gd_1_48, gd_1_49, gd_1_50, \
                         gd_1_51, gd_1_52, hd_1_48, hd_1_49, hd_1_50, hd_1_51, \
                         hd_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * gd_1_48[k]
                      + hd_1_48[k];

            t_81[k] = ab_x[k] * gd_1_49[k]
                      + hd_1_49[k];

            t_82[k] = ab_x[k] * gd_1_50[k]
                      + hd_1_50[k];

            t_83[k] = ab_x[k] * gd_1_51[k]
                      + hd_1_51[k];

            t_84[k] = ab_x[k] * gd_1_52[k]
                      + hd_1_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, ab_x, ab_y, gd_1_51, gd_1_52, gd_1_53, \
                         hd_1_53, hd_1_75, hd_1_76, hd_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * gd_1_53[k]
                      + hd_1_53[k];

            t_86[k] = ab_y[k] * gd_1_51[k]
                      + hd_1_75[k];

            t_87[k] = ab_y[k] * gd_1_52[k]
                      + hd_1_76[k];

            t_88[k] = ab_y[k] * gd_1_53[k]
                      + hd_1_77[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, ab_x, ab_z, gd_1_53, gd_1_54, gd_1_55, \
                         gd_1_56, gd_0_53, hd_1_54, hd_1_55, hd_1_56, \
                         hd_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_z[k] * gd_1_53[k]
                      + gd_0_53[k]
                      + hd_1_83[k];

            t_90[k] = ab_x[k] * gd_1_54[k]
                      + hd_1_54[k];

            t_91[k] = ab_x[k] * gd_1_55[k]
                      + hd_1_55[k];

            t_92[k] = ab_x[k] * gd_1_56[k]
                      + hd_1_56[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, ab_x, ab_y, gd_1_57, gd_1_58, gd_1_59, \
                         hd_1_57, hd_1_58, hd_1_59, hd_1_81, hd_1_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = ab_x[k] * gd_1_57[k]
                      + hd_1_57[k];

            t_94[k] = ab_x[k] * gd_1_58[k]
                      + hd_1_58[k];

            t_95[k] = ab_x[k] * gd_1_59[k]
                      + hd_1_59[k];

            t_96[k] = ab_y[k] * gd_1_57[k]
                      + hd_1_81[k];

            t_97[k] = ab_y[k] * gd_1_58[k]
                      + hd_1_82[k];
        }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, ab_x, ab_y, ab_z, gd_1_59, gd_1_60, \
                         gd_1_61, gd_0_59, hd_1_60, hd_1_61, hd_1_83, \
                         hd_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_98[k] = ab_y[k] * gd_1_59[k]
                      + hd_1_83[k];

            t_99[k] = ab_z[k] * gd_1_59[k]
                      + gd_0_59[k]
                      + hd_1_89[k];

            t_100[k] = ab_x[k] * gd_1_60[k]
                       + hd_1_60[k];

            t_101[k] = ab_x[k] * gd_1_61[k]
                       + hd_1_61[k];
        }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, ab_x, ab_y, gd_1_62, gd_1_63, \
                         gd_1_64, gd_1_65, hd_1_62, hd_1_63, hd_1_64, hd_1_65, \
                         hd_1_93 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_102[k] = ab_x[k] * gd_1_62[k]
                       + hd_1_62[k];

            t_103[k] = ab_x[k] * gd_1_63[k]
                       + hd_1_63[k];

            t_104[k] = ab_x[k] * gd_1_64[k]
                       + hd_1_64[k];

            t_105[k] = ab_x[k] * gd_1_65[k]
                       + hd_1_65[k];

            t_106[k] = ab_y[k] * gd_1_63[k]
                       + hd_1_93[k];
        }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, ab_x, ab_y, ab_z, gd_1_64, gd_1_65, \
                         gd_1_66, gd_0_65, hd_1_66, hd_1_94, hd_1_95, \
                         hd_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_107[k] = ab_y[k] * gd_1_64[k]
                       + hd_1_94[k];

            t_108[k] = ab_y[k] * gd_1_65[k]
                       + hd_1_95[k];

            t_109[k] = ab_z[k] * gd_1_65[k]
                       + gd_0_65[k]
                       + hd_1_101[k];

            t_110[k] = ab_x[k] * gd_1_66[k]
                       + hd_1_66[k];
        }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, ab_x, gd_1_67, gd_1_68, gd_1_69, \
                         gd_1_70, gd_1_71, hd_1_67, hd_1_68, hd_1_69, hd_1_70, \
                         hd_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_111[k] = ab_x[k] * gd_1_67[k]
                       + hd_1_67[k];

            t_112[k] = ab_x[k] * gd_1_68[k]
                       + hd_1_68[k];

            t_113[k] = ab_x[k] * gd_1_69[k]
                       + hd_1_69[k];

            t_114[k] = ab_x[k] * gd_1_70[k]
                       + hd_1_70[k];

            t_115[k] = ab_x[k] * gd_1_71[k]
                       + hd_1_71[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, ab_y, ab_z, gd_1_69, gd_1_70, gd_1_71, \
                         gd_0_71, hd_1_99, hd_1_100, hd_1_101, \
                         hd_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_y[k] * gd_1_69[k]
                       + hd_1_99[k];

            t_117[k] = ab_y[k] * gd_1_70[k]
                       + hd_1_100[k];

            t_118[k] = ab_y[k] * gd_1_71[k]
                       + hd_1_101[k];

            t_119[k] = ab_z[k] * gd_1_71[k]
                       + gd_0_71[k]
                       + hd_1_107[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gd_1_72, gd_1_73, gd_1_74, \
                         gd_1_75, gd_1_76, hd_1_72, hd_1_73, hd_1_74, hd_1_75, \
                         hd_1_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * gd_1_72[k]
                       + hd_1_72[k];

            t_121[k] = ab_x[k] * gd_1_73[k]
                       + hd_1_73[k];

            t_122[k] = ab_x[k] * gd_1_74[k]
                       + hd_1_74[k];

            t_123[k] = ab_x[k] * gd_1_75[k]
                       + hd_1_75[k];

            t_124[k] = ab_x[k] * gd_1_76[k]
                       + hd_1_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, ab_x, ab_y, gd_1_75, gd_1_76, gd_1_77, \
                         hd_1_77, hd_1_105, hd_1_106, hd_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * gd_1_77[k]
                       + hd_1_77[k];

            t_126[k] = ab_y[k] * gd_1_75[k]
                       + hd_1_105[k];

            t_127[k] = ab_y[k] * gd_1_76[k]
                       + hd_1_106[k];

            t_128[k] = ab_y[k] * gd_1_77[k]
                       + hd_1_107[k];
        }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, ab_x, ab_z, gd_1_77, gd_1_78, gd_1_79, \
                         gd_1_80, gd_0_77, hd_1_78, hd_1_79, hd_1_80, \
                         hd_1_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_129[k] = ab_z[k] * gd_1_77[k]
                       + gd_0_77[k]
                       + hd_1_113[k];

            t_130[k] = ab_x[k] * gd_1_78[k]
                       + hd_1_78[k];

            t_131[k] = ab_x[k] * gd_1_79[k]
                       + hd_1_79[k];

            t_132[k] = ab_x[k] * gd_1_80[k]
                       + hd_1_80[k];
        }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, ab_x, ab_y, gd_1_81, gd_1_82, \
                         gd_1_83, hd_1_81, hd_1_82, hd_1_83, hd_1_111, \
                         hd_1_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_133[k] = ab_x[k] * gd_1_81[k]
                       + hd_1_81[k];

            t_134[k] = ab_x[k] * gd_1_82[k]
                       + hd_1_82[k];

            t_135[k] = ab_x[k] * gd_1_83[k]
                       + hd_1_83[k];

            t_136[k] = ab_y[k] * gd_1_81[k]
                       + hd_1_111[k];

            t_137[k] = ab_y[k] * gd_1_82[k]
                       + hd_1_112[k];
        }
    }
}

static auto
compute_hrr_geom_100z_gf_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t gd_1,
                                             const size_t gd_0, const size_t hd_1,
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

        const auto *gd_1_83 = buffer.data(gd_1 + 83 * ncomps + c);
        const auto *gd_1_84 = buffer.data(gd_1 + 84 * ncomps + c);
        const auto *gd_1_85 = buffer.data(gd_1 + 85 * ncomps + c);
        const auto *gd_1_86 = buffer.data(gd_1 + 86 * ncomps + c);
        const auto *gd_1_87 = buffer.data(gd_1 + 87 * ncomps + c);
        const auto *gd_1_88 = buffer.data(gd_1 + 88 * ncomps + c);
        const auto *gd_1_89 = buffer.data(gd_1 + 89 * ncomps + c);

        const auto *gd_0_83 = buffer.data(gd_0 + 83 * ncomps + c);
        const auto *gd_0_89 = buffer.data(gd_0 + 89 * ncomps + c);

        const auto *hd_1_84 = buffer.data(hd_1 + 84 * ncomps + c);
        const auto *hd_1_85 = buffer.data(hd_1 + 85 * ncomps + c);
        const auto *hd_1_86 = buffer.data(hd_1 + 86 * ncomps + c);
        const auto *hd_1_87 = buffer.data(hd_1 + 87 * ncomps + c);
        const auto *hd_1_88 = buffer.data(hd_1 + 88 * ncomps + c);
        const auto *hd_1_89 = buffer.data(hd_1 + 89 * ncomps + c);
        const auto *hd_1_113 = buffer.data(hd_1 + 113 * ncomps + c);
        const auto *hd_1_117 = buffer.data(hd_1 + 117 * ncomps + c);
        const auto *hd_1_118 = buffer.data(hd_1 + 118 * ncomps + c);
        const auto *hd_1_119 = buffer.data(hd_1 + 119 * ncomps + c);
        const auto *hd_1_125 = buffer.data(hd_1 + 125 * ncomps + c);

#pragma omp simd aligned(t_138, t_139, t_140, t_141, ab_x, ab_y, ab_z, gd_1_83, gd_1_84, \
                         gd_1_85, gd_0_83, hd_1_84, hd_1_85, hd_1_113, \
                         hd_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_138[k] = ab_y[k] * gd_1_83[k]
                       + hd_1_113[k];

            t_139[k] = ab_z[k] * gd_1_83[k]
                       + gd_0_83[k]
                       + hd_1_119[k];

            t_140[k] = ab_x[k] * gd_1_84[k]
                       + hd_1_84[k];

            t_141[k] = ab_x[k] * gd_1_85[k]
                       + hd_1_85[k];
        }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, ab_x, ab_y, gd_1_86, gd_1_87, \
                         gd_1_88, gd_1_89, hd_1_86, hd_1_87, hd_1_88, hd_1_89, \
                         hd_1_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_142[k] = ab_x[k] * gd_1_86[k]
                       + hd_1_86[k];

            t_143[k] = ab_x[k] * gd_1_87[k]
                       + hd_1_87[k];

            t_144[k] = ab_x[k] * gd_1_88[k]
                       + hd_1_88[k];

            t_145[k] = ab_x[k] * gd_1_89[k]
                       + hd_1_89[k];

            t_146[k] = ab_y[k] * gd_1_87[k]
                       + hd_1_117[k];
        }

#pragma omp simd aligned(t_147, t_148, t_149, ab_y, ab_z, gd_1_88, gd_1_89, gd_0_89, hd_1_118, \
                         hd_1_119, hd_1_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_147[k] = ab_y[k] * gd_1_88[k]
                       + hd_1_118[k];

            t_148[k] = ab_y[k] * gd_1_89[k]
                       + hd_1_119[k];

            t_149[k] = ab_z[k] * gd_1_89[k]
                       + gd_0_89[k]
                       + hd_1_125[k];
        }
    }
}

auto
compute_hrr_geom_100z_gf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t gd_1, const size_t gd_0,
                                      const size_t hd_1, const size_t ncomps,
                                      const size_t nmax) -> void
{
    compute_hrr_geom_100z_gf_out_of_first_piece0(buffer, coordinates, target, gd_1, gd_0, hd_1,
                                                 ncomps, nmax);

    compute_hrr_geom_100z_gf_out_of_first_piece1(buffer, coordinates, target, gd_1, gd_0, hd_1,
                                                 ncomps, nmax);
}

}  // namespace simdtrf
