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


#include "SimdTransferGeom010YFH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_010y_fh_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t dh_1, const size_t dh_0,
                                const size_t di_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);

        const auto *dh_1_0 = buffer.data(dh_1 + 0 * ncomps + c);
        const auto *dh_1_1 = buffer.data(dh_1 + 1 * ncomps + c);
        const auto *dh_1_2 = buffer.data(dh_1 + 2 * ncomps + c);
        const auto *dh_1_3 = buffer.data(dh_1 + 3 * ncomps + c);
        const auto *dh_1_4 = buffer.data(dh_1 + 4 * ncomps + c);
        const auto *dh_1_5 = buffer.data(dh_1 + 5 * ncomps + c);
        const auto *dh_1_6 = buffer.data(dh_1 + 6 * ncomps + c);
        const auto *dh_1_7 = buffer.data(dh_1 + 7 * ncomps + c);
        const auto *dh_1_8 = buffer.data(dh_1 + 8 * ncomps + c);
        const auto *dh_1_9 = buffer.data(dh_1 + 9 * ncomps + c);
        const auto *dh_1_10 = buffer.data(dh_1 + 10 * ncomps + c);
        const auto *dh_1_11 = buffer.data(dh_1 + 11 * ncomps + c);
        const auto *dh_1_12 = buffer.data(dh_1 + 12 * ncomps + c);
        const auto *dh_1_13 = buffer.data(dh_1 + 13 * ncomps + c);
        const auto *dh_1_14 = buffer.data(dh_1 + 14 * ncomps + c);
        const auto *dh_1_15 = buffer.data(dh_1 + 15 * ncomps + c);
        const auto *dh_1_16 = buffer.data(dh_1 + 16 * ncomps + c);
        const auto *dh_1_17 = buffer.data(dh_1 + 17 * ncomps + c);
        const auto *dh_1_18 = buffer.data(dh_1 + 18 * ncomps + c);
        const auto *dh_1_19 = buffer.data(dh_1 + 19 * ncomps + c);
        const auto *dh_1_20 = buffer.data(dh_1 + 20 * ncomps + c);
        const auto *dh_1_21 = buffer.data(dh_1 + 21 * ncomps + c);
        const auto *dh_1_22 = buffer.data(dh_1 + 22 * ncomps + c);
        const auto *dh_1_23 = buffer.data(dh_1 + 23 * ncomps + c);
        const auto *dh_1_24 = buffer.data(dh_1 + 24 * ncomps + c);
        const auto *dh_1_25 = buffer.data(dh_1 + 25 * ncomps + c);
        const auto *dh_1_26 = buffer.data(dh_1 + 26 * ncomps + c);
        const auto *dh_1_27 = buffer.data(dh_1 + 27 * ncomps + c);
        const auto *dh_1_28 = buffer.data(dh_1 + 28 * ncomps + c);
        const auto *dh_1_29 = buffer.data(dh_1 + 29 * ncomps + c);
        const auto *dh_1_30 = buffer.data(dh_1 + 30 * ncomps + c);
        const auto *dh_1_31 = buffer.data(dh_1 + 31 * ncomps + c);
        const auto *dh_1_32 = buffer.data(dh_1 + 32 * ncomps + c);
        const auto *dh_1_33 = buffer.data(dh_1 + 33 * ncomps + c);
        const auto *dh_1_34 = buffer.data(dh_1 + 34 * ncomps + c);
        const auto *dh_1_35 = buffer.data(dh_1 + 35 * ncomps + c);
        const auto *dh_1_36 = buffer.data(dh_1 + 36 * ncomps + c);
        const auto *dh_1_37 = buffer.data(dh_1 + 37 * ncomps + c);
        const auto *dh_1_38 = buffer.data(dh_1 + 38 * ncomps + c);
        const auto *dh_1_39 = buffer.data(dh_1 + 39 * ncomps + c);
        const auto *dh_1_40 = buffer.data(dh_1 + 40 * ncomps + c);
        const auto *dh_1_41 = buffer.data(dh_1 + 41 * ncomps + c);
        const auto *dh_1_42 = buffer.data(dh_1 + 42 * ncomps + c);
        const auto *dh_1_43 = buffer.data(dh_1 + 43 * ncomps + c);
        const auto *dh_1_44 = buffer.data(dh_1 + 44 * ncomps + c);
        const auto *dh_1_45 = buffer.data(dh_1 + 45 * ncomps + c);
        const auto *dh_1_46 = buffer.data(dh_1 + 46 * ncomps + c);
        const auto *dh_1_47 = buffer.data(dh_1 + 47 * ncomps + c);
        const auto *dh_1_48 = buffer.data(dh_1 + 48 * ncomps + c);
        const auto *dh_1_49 = buffer.data(dh_1 + 49 * ncomps + c);
        const auto *dh_1_50 = buffer.data(dh_1 + 50 * ncomps + c);
        const auto *dh_1_51 = buffer.data(dh_1 + 51 * ncomps + c);
        const auto *dh_1_52 = buffer.data(dh_1 + 52 * ncomps + c);
        const auto *dh_1_53 = buffer.data(dh_1 + 53 * ncomps + c);
        const auto *dh_1_54 = buffer.data(dh_1 + 54 * ncomps + c);
        const auto *dh_1_55 = buffer.data(dh_1 + 55 * ncomps + c);
        const auto *dh_1_56 = buffer.data(dh_1 + 56 * ncomps + c);
        const auto *dh_1_57 = buffer.data(dh_1 + 57 * ncomps + c);
        const auto *dh_1_58 = buffer.data(dh_1 + 58 * ncomps + c);
        const auto *dh_1_59 = buffer.data(dh_1 + 59 * ncomps + c);
        const auto *dh_1_60 = buffer.data(dh_1 + 60 * ncomps + c);
        const auto *dh_1_61 = buffer.data(dh_1 + 61 * ncomps + c);
        const auto *dh_1_62 = buffer.data(dh_1 + 62 * ncomps + c);
        const auto *dh_1_63 = buffer.data(dh_1 + 63 * ncomps + c);
        const auto *dh_1_64 = buffer.data(dh_1 + 64 * ncomps + c);
        const auto *dh_1_65 = buffer.data(dh_1 + 65 * ncomps + c);
        const auto *dh_1_66 = buffer.data(dh_1 + 66 * ncomps + c);
        const auto *dh_1_67 = buffer.data(dh_1 + 67 * ncomps + c);
        const auto *dh_1_68 = buffer.data(dh_1 + 68 * ncomps + c);
        const auto *dh_1_69 = buffer.data(dh_1 + 69 * ncomps + c);
        const auto *dh_1_70 = buffer.data(dh_1 + 70 * ncomps + c);
        const auto *dh_1_71 = buffer.data(dh_1 + 71 * ncomps + c);
        const auto *dh_1_72 = buffer.data(dh_1 + 72 * ncomps + c);
        const auto *dh_1_73 = buffer.data(dh_1 + 73 * ncomps + c);
        const auto *dh_1_74 = buffer.data(dh_1 + 74 * ncomps + c);
        const auto *dh_1_75 = buffer.data(dh_1 + 75 * ncomps + c);
        const auto *dh_1_76 = buffer.data(dh_1 + 76 * ncomps + c);
        const auto *dh_1_77 = buffer.data(dh_1 + 77 * ncomps + c);
        const auto *dh_1_78 = buffer.data(dh_1 + 78 * ncomps + c);
        const auto *dh_1_79 = buffer.data(dh_1 + 79 * ncomps + c);
        const auto *dh_1_80 = buffer.data(dh_1 + 80 * ncomps + c);
        const auto *dh_1_81 = buffer.data(dh_1 + 81 * ncomps + c);
        const auto *dh_1_82 = buffer.data(dh_1 + 82 * ncomps + c);
        const auto *dh_1_83 = buffer.data(dh_1 + 83 * ncomps + c);
        const auto *dh_1_84 = buffer.data(dh_1 + 84 * ncomps + c);
        const auto *dh_1_85 = buffer.data(dh_1 + 85 * ncomps + c);
        const auto *dh_1_86 = buffer.data(dh_1 + 86 * ncomps + c);
        const auto *dh_1_87 = buffer.data(dh_1 + 87 * ncomps + c);
        const auto *dh_1_88 = buffer.data(dh_1 + 88 * ncomps + c);
        const auto *dh_1_89 = buffer.data(dh_1 + 89 * ncomps + c);
        const auto *dh_1_90 = buffer.data(dh_1 + 90 * ncomps + c);
        const auto *dh_1_91 = buffer.data(dh_1 + 91 * ncomps + c);
        const auto *dh_1_92 = buffer.data(dh_1 + 92 * ncomps + c);
        const auto *dh_1_93 = buffer.data(dh_1 + 93 * ncomps + c);
        const auto *dh_1_94 = buffer.data(dh_1 + 94 * ncomps + c);
        const auto *dh_1_95 = buffer.data(dh_1 + 95 * ncomps + c);
        const auto *dh_1_96 = buffer.data(dh_1 + 96 * ncomps + c);
        const auto *dh_1_97 = buffer.data(dh_1 + 97 * ncomps + c);
        const auto *dh_1_98 = buffer.data(dh_1 + 98 * ncomps + c);
        const auto *dh_1_99 = buffer.data(dh_1 + 99 * ncomps + c);
        const auto *dh_1_100 = buffer.data(dh_1 + 100 * ncomps + c);
        const auto *dh_1_101 = buffer.data(dh_1 + 101 * ncomps + c);
        const auto *dh_1_102 = buffer.data(dh_1 + 102 * ncomps + c);
        const auto *dh_1_103 = buffer.data(dh_1 + 103 * ncomps + c);
        const auto *dh_1_104 = buffer.data(dh_1 + 104 * ncomps + c);
        const auto *dh_1_105 = buffer.data(dh_1 + 105 * ncomps + c);
        const auto *dh_1_106 = buffer.data(dh_1 + 106 * ncomps + c);
        const auto *dh_1_107 = buffer.data(dh_1 + 107 * ncomps + c);
        const auto *dh_1_108 = buffer.data(dh_1 + 108 * ncomps + c);
        const auto *dh_1_109 = buffer.data(dh_1 + 109 * ncomps + c);
        const auto *dh_1_110 = buffer.data(dh_1 + 110 * ncomps + c);
        const auto *dh_1_111 = buffer.data(dh_1 + 111 * ncomps + c);
        const auto *dh_1_112 = buffer.data(dh_1 + 112 * ncomps + c);
        const auto *dh_1_113 = buffer.data(dh_1 + 113 * ncomps + c);
        const auto *dh_1_114 = buffer.data(dh_1 + 114 * ncomps + c);
        const auto *dh_1_115 = buffer.data(dh_1 + 115 * ncomps + c);
        const auto *dh_1_116 = buffer.data(dh_1 + 116 * ncomps + c);
        const auto *dh_1_117 = buffer.data(dh_1 + 117 * ncomps + c);
        const auto *dh_1_118 = buffer.data(dh_1 + 118 * ncomps + c);
        const auto *dh_1_119 = buffer.data(dh_1 + 119 * ncomps + c);
        const auto *dh_1_120 = buffer.data(dh_1 + 120 * ncomps + c);
        const auto *dh_1_121 = buffer.data(dh_1 + 121 * ncomps + c);
        const auto *dh_1_122 = buffer.data(dh_1 + 122 * ncomps + c);
        const auto *dh_1_123 = buffer.data(dh_1 + 123 * ncomps + c);
        const auto *dh_1_124 = buffer.data(dh_1 + 124 * ncomps + c);
        const auto *dh_1_125 = buffer.data(dh_1 + 125 * ncomps + c);

        const auto *dh_0_63 = buffer.data(dh_0 + 63 * ncomps + c);
        const auto *dh_0_64 = buffer.data(dh_0 + 64 * ncomps + c);
        const auto *dh_0_65 = buffer.data(dh_0 + 65 * ncomps + c);
        const auto *dh_0_66 = buffer.data(dh_0 + 66 * ncomps + c);
        const auto *dh_0_67 = buffer.data(dh_0 + 67 * ncomps + c);
        const auto *dh_0_68 = buffer.data(dh_0 + 68 * ncomps + c);
        const auto *dh_0_69 = buffer.data(dh_0 + 69 * ncomps + c);
        const auto *dh_0_70 = buffer.data(dh_0 + 70 * ncomps + c);
        const auto *dh_0_71 = buffer.data(dh_0 + 71 * ncomps + c);
        const auto *dh_0_72 = buffer.data(dh_0 + 72 * ncomps + c);
        const auto *dh_0_73 = buffer.data(dh_0 + 73 * ncomps + c);
        const auto *dh_0_74 = buffer.data(dh_0 + 74 * ncomps + c);
        const auto *dh_0_75 = buffer.data(dh_0 + 75 * ncomps + c);
        const auto *dh_0_76 = buffer.data(dh_0 + 76 * ncomps + c);

        const auto *di_1_0 = buffer.data(di_1 + 0 * ncomps + c);
        const auto *di_1_1 = buffer.data(di_1 + 1 * ncomps + c);
        const auto *di_1_2 = buffer.data(di_1 + 2 * ncomps + c);
        const auto *di_1_3 = buffer.data(di_1 + 3 * ncomps + c);
        const auto *di_1_4 = buffer.data(di_1 + 4 * ncomps + c);
        const auto *di_1_5 = buffer.data(di_1 + 5 * ncomps + c);
        const auto *di_1_6 = buffer.data(di_1 + 6 * ncomps + c);
        const auto *di_1_7 = buffer.data(di_1 + 7 * ncomps + c);
        const auto *di_1_8 = buffer.data(di_1 + 8 * ncomps + c);
        const auto *di_1_9 = buffer.data(di_1 + 9 * ncomps + c);
        const auto *di_1_10 = buffer.data(di_1 + 10 * ncomps + c);
        const auto *di_1_11 = buffer.data(di_1 + 11 * ncomps + c);
        const auto *di_1_12 = buffer.data(di_1 + 12 * ncomps + c);
        const auto *di_1_13 = buffer.data(di_1 + 13 * ncomps + c);
        const auto *di_1_14 = buffer.data(di_1 + 14 * ncomps + c);
        const auto *di_1_15 = buffer.data(di_1 + 15 * ncomps + c);
        const auto *di_1_16 = buffer.data(di_1 + 16 * ncomps + c);
        const auto *di_1_17 = buffer.data(di_1 + 17 * ncomps + c);
        const auto *di_1_18 = buffer.data(di_1 + 18 * ncomps + c);
        const auto *di_1_19 = buffer.data(di_1 + 19 * ncomps + c);
        const auto *di_1_20 = buffer.data(di_1 + 20 * ncomps + c);
        const auto *di_1_28 = buffer.data(di_1 + 28 * ncomps + c);
        const auto *di_1_29 = buffer.data(di_1 + 29 * ncomps + c);
        const auto *di_1_30 = buffer.data(di_1 + 30 * ncomps + c);
        const auto *di_1_31 = buffer.data(di_1 + 31 * ncomps + c);
        const auto *di_1_32 = buffer.data(di_1 + 32 * ncomps + c);
        const auto *di_1_33 = buffer.data(di_1 + 33 * ncomps + c);
        const auto *di_1_34 = buffer.data(di_1 + 34 * ncomps + c);
        const auto *di_1_35 = buffer.data(di_1 + 35 * ncomps + c);
        const auto *di_1_36 = buffer.data(di_1 + 36 * ncomps + c);
        const auto *di_1_37 = buffer.data(di_1 + 37 * ncomps + c);
        const auto *di_1_38 = buffer.data(di_1 + 38 * ncomps + c);
        const auto *di_1_39 = buffer.data(di_1 + 39 * ncomps + c);
        const auto *di_1_40 = buffer.data(di_1 + 40 * ncomps + c);
        const auto *di_1_41 = buffer.data(di_1 + 41 * ncomps + c);
        const auto *di_1_42 = buffer.data(di_1 + 42 * ncomps + c);
        const auto *di_1_43 = buffer.data(di_1 + 43 * ncomps + c);
        const auto *di_1_44 = buffer.data(di_1 + 44 * ncomps + c);
        const auto *di_1_45 = buffer.data(di_1 + 45 * ncomps + c);
        const auto *di_1_46 = buffer.data(di_1 + 46 * ncomps + c);
        const auto *di_1_47 = buffer.data(di_1 + 47 * ncomps + c);
        const auto *di_1_48 = buffer.data(di_1 + 48 * ncomps + c);
        const auto *di_1_56 = buffer.data(di_1 + 56 * ncomps + c);
        const auto *di_1_57 = buffer.data(di_1 + 57 * ncomps + c);
        const auto *di_1_58 = buffer.data(di_1 + 58 * ncomps + c);
        const auto *di_1_59 = buffer.data(di_1 + 59 * ncomps + c);
        const auto *di_1_60 = buffer.data(di_1 + 60 * ncomps + c);
        const auto *di_1_61 = buffer.data(di_1 + 61 * ncomps + c);
        const auto *di_1_62 = buffer.data(di_1 + 62 * ncomps + c);
        const auto *di_1_63 = buffer.data(di_1 + 63 * ncomps + c);
        const auto *di_1_64 = buffer.data(di_1 + 64 * ncomps + c);
        const auto *di_1_65 = buffer.data(di_1 + 65 * ncomps + c);
        const auto *di_1_66 = buffer.data(di_1 + 66 * ncomps + c);
        const auto *di_1_67 = buffer.data(di_1 + 67 * ncomps + c);
        const auto *di_1_68 = buffer.data(di_1 + 68 * ncomps + c);
        const auto *di_1_69 = buffer.data(di_1 + 69 * ncomps + c);
        const auto *di_1_70 = buffer.data(di_1 + 70 * ncomps + c);
        const auto *di_1_71 = buffer.data(di_1 + 71 * ncomps + c);
        const auto *di_1_72 = buffer.data(di_1 + 72 * ncomps + c);
        const auto *di_1_73 = buffer.data(di_1 + 73 * ncomps + c);
        const auto *di_1_74 = buffer.data(di_1 + 74 * ncomps + c);
        const auto *di_1_75 = buffer.data(di_1 + 75 * ncomps + c);
        const auto *di_1_76 = buffer.data(di_1 + 76 * ncomps + c);
        const auto *di_1_84 = buffer.data(di_1 + 84 * ncomps + c);
        const auto *di_1_85 = buffer.data(di_1 + 85 * ncomps + c);
        const auto *di_1_86 = buffer.data(di_1 + 86 * ncomps + c);
        const auto *di_1_87 = buffer.data(di_1 + 87 * ncomps + c);
        const auto *di_1_88 = buffer.data(di_1 + 88 * ncomps + c);
        const auto *di_1_89 = buffer.data(di_1 + 89 * ncomps + c);
        const auto *di_1_90 = buffer.data(di_1 + 90 * ncomps + c);
        const auto *di_1_91 = buffer.data(di_1 + 91 * ncomps + c);
        const auto *di_1_92 = buffer.data(di_1 + 92 * ncomps + c);
        const auto *di_1_93 = buffer.data(di_1 + 93 * ncomps + c);
        const auto *di_1_94 = buffer.data(di_1 + 94 * ncomps + c);
        const auto *di_1_95 = buffer.data(di_1 + 95 * ncomps + c);
        const auto *di_1_96 = buffer.data(di_1 + 96 * ncomps + c);
        const auto *di_1_97 = buffer.data(di_1 + 97 * ncomps + c);
        const auto *di_1_98 = buffer.data(di_1 + 98 * ncomps + c);
        const auto *di_1_99 = buffer.data(di_1 + 99 * ncomps + c);
        const auto *di_1_100 = buffer.data(di_1 + 100 * ncomps + c);
        const auto *di_1_101 = buffer.data(di_1 + 101 * ncomps + c);
        const auto *di_1_102 = buffer.data(di_1 + 102 * ncomps + c);
        const auto *di_1_103 = buffer.data(di_1 + 103 * ncomps + c);
        const auto *di_1_104 = buffer.data(di_1 + 104 * ncomps + c);
        const auto *di_1_112 = buffer.data(di_1 + 112 * ncomps + c);
        const auto *di_1_113 = buffer.data(di_1 + 113 * ncomps + c);
        const auto *di_1_114 = buffer.data(di_1 + 114 * ncomps + c);
        const auto *di_1_115 = buffer.data(di_1 + 115 * ncomps + c);
        const auto *di_1_116 = buffer.data(di_1 + 116 * ncomps + c);
        const auto *di_1_117 = buffer.data(di_1 + 117 * ncomps + c);
        const auto *di_1_118 = buffer.data(di_1 + 118 * ncomps + c);
        const auto *di_1_119 = buffer.data(di_1 + 119 * ncomps + c);
        const auto *di_1_120 = buffer.data(di_1 + 120 * ncomps + c);
        const auto *di_1_121 = buffer.data(di_1 + 121 * ncomps + c);
        const auto *di_1_122 = buffer.data(di_1 + 122 * ncomps + c);
        const auto *di_1_123 = buffer.data(di_1 + 123 * ncomps + c);
        const auto *di_1_124 = buffer.data(di_1 + 124 * ncomps + c);
        const auto *di_1_125 = buffer.data(di_1 + 125 * ncomps + c);
        const auto *di_1_126 = buffer.data(di_1 + 126 * ncomps + c);
        const auto *di_1_127 = buffer.data(di_1 + 127 * ncomps + c);
        const auto *di_1_128 = buffer.data(di_1 + 128 * ncomps + c);
        const auto *di_1_129 = buffer.data(di_1 + 129 * ncomps + c);
        const auto *di_1_130 = buffer.data(di_1 + 130 * ncomps + c);
        const auto *di_1_131 = buffer.data(di_1 + 131 * ncomps + c);
        const auto *di_1_132 = buffer.data(di_1 + 132 * ncomps + c);
        const auto *di_1_140 = buffer.data(di_1 + 140 * ncomps + c);
        const auto *di_1_141 = buffer.data(di_1 + 141 * ncomps + c);
        const auto *di_1_142 = buffer.data(di_1 + 142 * ncomps + c);
        const auto *di_1_143 = buffer.data(di_1 + 143 * ncomps + c);
        const auto *di_1_144 = buffer.data(di_1 + 144 * ncomps + c);
        const auto *di_1_145 = buffer.data(di_1 + 145 * ncomps + c);
        const auto *di_1_146 = buffer.data(di_1 + 146 * ncomps + c);
        const auto *di_1_147 = buffer.data(di_1 + 147 * ncomps + c);
        const auto *di_1_148 = buffer.data(di_1 + 148 * ncomps + c);
        const auto *di_1_149 = buffer.data(di_1 + 149 * ncomps + c);
        const auto *di_1_150 = buffer.data(di_1 + 150 * ncomps + c);
        const auto *di_1_151 = buffer.data(di_1 + 151 * ncomps + c);
        const auto *di_1_152 = buffer.data(di_1 + 152 * ncomps + c);
        const auto *di_1_153 = buffer.data(di_1 + 153 * ncomps + c);
        const auto *di_1_154 = buffer.data(di_1 + 154 * ncomps + c);
        const auto *di_1_155 = buffer.data(di_1 + 155 * ncomps + c);
        const auto *di_1_156 = buffer.data(di_1 + 156 * ncomps + c);
        const auto *di_1_157 = buffer.data(di_1 + 157 * ncomps + c);
        const auto *di_1_158 = buffer.data(di_1 + 158 * ncomps + c);
        const auto *di_1_159 = buffer.data(di_1 + 159 * ncomps + c);
        const auto *di_1_160 = buffer.data(di_1 + 160 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dh_1_0, dh_1_1, dh_1_2, dh_1_3, \
                         dh_1_4, di_1_0, di_1_1, di_1_2, di_1_3, \
                         di_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dh_1_0[k]
                     + di_1_0[k];

            t_1[k] = -ab_x[k] * dh_1_1[k]
                     + di_1_1[k];

            t_2[k] = -ab_x[k] * dh_1_2[k]
                     + di_1_2[k];

            t_3[k] = -ab_x[k] * dh_1_3[k]
                     + di_1_3[k];

            t_4[k] = -ab_x[k] * dh_1_4[k]
                     + di_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dh_1_5, dh_1_6, dh_1_7, dh_1_8, \
                         dh_1_9, di_1_5, di_1_6, di_1_7, di_1_8, \
                         di_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dh_1_5[k]
                     + di_1_5[k];

            t_6[k] = -ab_x[k] * dh_1_6[k]
                     + di_1_6[k];

            t_7[k] = -ab_x[k] * dh_1_7[k]
                     + di_1_7[k];

            t_8[k] = -ab_x[k] * dh_1_8[k]
                     + di_1_8[k];

            t_9[k] = -ab_x[k] * dh_1_9[k]
                     + di_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dh_1_10, dh_1_11, dh_1_12, \
                         dh_1_13, dh_1_14, di_1_10, di_1_11, di_1_12, di_1_13, \
                         di_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dh_1_10[k]
                      + di_1_10[k];

            t_11[k] = -ab_x[k] * dh_1_11[k]
                      + di_1_11[k];

            t_12[k] = -ab_x[k] * dh_1_12[k]
                      + di_1_12[k];

            t_13[k] = -ab_x[k] * dh_1_13[k]
                      + di_1_13[k];

            t_14[k] = -ab_x[k] * dh_1_14[k]
                      + di_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dh_1_15, dh_1_16, dh_1_17, \
                         dh_1_18, dh_1_19, di_1_15, di_1_16, di_1_17, di_1_18, \
                         di_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dh_1_15[k]
                      + di_1_15[k];

            t_16[k] = -ab_x[k] * dh_1_16[k]
                      + di_1_16[k];

            t_17[k] = -ab_x[k] * dh_1_17[k]
                      + di_1_17[k];

            t_18[k] = -ab_x[k] * dh_1_18[k]
                      + di_1_18[k];

            t_19[k] = -ab_x[k] * dh_1_19[k]
                      + di_1_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dh_1_20, dh_1_21, dh_1_22, \
                         dh_1_23, dh_1_24, di_1_20, di_1_28, di_1_29, di_1_30, \
                         di_1_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * dh_1_20[k]
                      + di_1_20[k];

            t_21[k] = -ab_x[k] * dh_1_21[k]
                      + di_1_28[k];

            t_22[k] = -ab_x[k] * dh_1_22[k]
                      + di_1_29[k];

            t_23[k] = -ab_x[k] * dh_1_23[k]
                      + di_1_30[k];

            t_24[k] = -ab_x[k] * dh_1_24[k]
                      + di_1_31[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dh_1_25, dh_1_26, dh_1_27, \
                         dh_1_28, dh_1_29, di_1_32, di_1_33, di_1_34, di_1_35, \
                         di_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * dh_1_25[k]
                      + di_1_32[k];

            t_26[k] = -ab_x[k] * dh_1_26[k]
                      + di_1_33[k];

            t_27[k] = -ab_x[k] * dh_1_27[k]
                      + di_1_34[k];

            t_28[k] = -ab_x[k] * dh_1_28[k]
                      + di_1_35[k];

            t_29[k] = -ab_x[k] * dh_1_29[k]
                      + di_1_36[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dh_1_30, dh_1_31, dh_1_32, \
                         dh_1_33, dh_1_34, di_1_37, di_1_38, di_1_39, di_1_40, \
                         di_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dh_1_30[k]
                      + di_1_37[k];

            t_31[k] = -ab_x[k] * dh_1_31[k]
                      + di_1_38[k];

            t_32[k] = -ab_x[k] * dh_1_32[k]
                      + di_1_39[k];

            t_33[k] = -ab_x[k] * dh_1_33[k]
                      + di_1_40[k];

            t_34[k] = -ab_x[k] * dh_1_34[k]
                      + di_1_41[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dh_1_35, dh_1_36, dh_1_37, \
                         dh_1_38, dh_1_39, di_1_42, di_1_43, di_1_44, di_1_45, \
                         di_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * dh_1_35[k]
                      + di_1_42[k];

            t_36[k] = -ab_x[k] * dh_1_36[k]
                      + di_1_43[k];

            t_37[k] = -ab_x[k] * dh_1_37[k]
                      + di_1_44[k];

            t_38[k] = -ab_x[k] * dh_1_38[k]
                      + di_1_45[k];

            t_39[k] = -ab_x[k] * dh_1_39[k]
                      + di_1_46[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dh_1_40, dh_1_41, dh_1_42, \
                         dh_1_43, dh_1_44, di_1_47, di_1_48, di_1_56, di_1_57, \
                         di_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * dh_1_40[k]
                      + di_1_47[k];

            t_41[k] = -ab_x[k] * dh_1_41[k]
                      + di_1_48[k];

            t_42[k] = -ab_x[k] * dh_1_42[k]
                      + di_1_56[k];

            t_43[k] = -ab_x[k] * dh_1_43[k]
                      + di_1_57[k];

            t_44[k] = -ab_x[k] * dh_1_44[k]
                      + di_1_58[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dh_1_45, dh_1_46, dh_1_47, \
                         dh_1_48, dh_1_49, di_1_59, di_1_60, di_1_61, di_1_62, \
                         di_1_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * dh_1_45[k]
                      + di_1_59[k];

            t_46[k] = -ab_x[k] * dh_1_46[k]
                      + di_1_60[k];

            t_47[k] = -ab_x[k] * dh_1_47[k]
                      + di_1_61[k];

            t_48[k] = -ab_x[k] * dh_1_48[k]
                      + di_1_62[k];

            t_49[k] = -ab_x[k] * dh_1_49[k]
                      + di_1_63[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dh_1_50, dh_1_51, dh_1_52, \
                         dh_1_53, dh_1_54, di_1_64, di_1_65, di_1_66, di_1_67, \
                         di_1_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * dh_1_50[k]
                      + di_1_64[k];

            t_51[k] = -ab_x[k] * dh_1_51[k]
                      + di_1_65[k];

            t_52[k] = -ab_x[k] * dh_1_52[k]
                      + di_1_66[k];

            t_53[k] = -ab_x[k] * dh_1_53[k]
                      + di_1_67[k];

            t_54[k] = -ab_x[k] * dh_1_54[k]
                      + di_1_68[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dh_1_55, dh_1_56, dh_1_57, \
                         dh_1_58, dh_1_59, di_1_69, di_1_70, di_1_71, di_1_72, \
                         di_1_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * dh_1_55[k]
                      + di_1_69[k];

            t_56[k] = -ab_x[k] * dh_1_56[k]
                      + di_1_70[k];

            t_57[k] = -ab_x[k] * dh_1_57[k]
                      + di_1_71[k];

            t_58[k] = -ab_x[k] * dh_1_58[k]
                      + di_1_72[k];

            t_59[k] = -ab_x[k] * dh_1_59[k]
                      + di_1_73[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dh_1_60, dh_1_61, dh_1_62, \
                         dh_1_63, dh_1_64, di_1_74, di_1_75, di_1_76, di_1_84, \
                         di_1_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * dh_1_60[k]
                      + di_1_74[k];

            t_61[k] = -ab_x[k] * dh_1_61[k]
                      + di_1_75[k];

            t_62[k] = -ab_x[k] * dh_1_62[k]
                      + di_1_76[k];

            t_63[k] = -ab_x[k] * dh_1_63[k]
                      + di_1_84[k];

            t_64[k] = -ab_x[k] * dh_1_64[k]
                      + di_1_85[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dh_1_65, dh_1_66, dh_1_67, \
                         dh_1_68, dh_1_69, di_1_86, di_1_87, di_1_88, di_1_89, \
                         di_1_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * dh_1_65[k]
                      + di_1_86[k];

            t_66[k] = -ab_x[k] * dh_1_66[k]
                      + di_1_87[k];

            t_67[k] = -ab_x[k] * dh_1_67[k]
                      + di_1_88[k];

            t_68[k] = -ab_x[k] * dh_1_68[k]
                      + di_1_89[k];

            t_69[k] = -ab_x[k] * dh_1_69[k]
                      + di_1_90[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dh_1_70, dh_1_71, dh_1_72, \
                         dh_1_73, dh_1_74, di_1_91, di_1_92, di_1_93, di_1_94, \
                         di_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * dh_1_70[k]
                      + di_1_91[k];

            t_71[k] = -ab_x[k] * dh_1_71[k]
                      + di_1_92[k];

            t_72[k] = -ab_x[k] * dh_1_72[k]
                      + di_1_93[k];

            t_73[k] = -ab_x[k] * dh_1_73[k]
                      + di_1_94[k];

            t_74[k] = -ab_x[k] * dh_1_74[k]
                      + di_1_95[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dh_1_75, dh_1_76, dh_1_77, \
                         dh_1_78, dh_1_79, di_1_96, di_1_97, di_1_98, di_1_99, \
                         di_1_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * dh_1_75[k]
                      + di_1_96[k];

            t_76[k] = -ab_x[k] * dh_1_76[k]
                      + di_1_97[k];

            t_77[k] = -ab_x[k] * dh_1_77[k]
                      + di_1_98[k];

            t_78[k] = -ab_x[k] * dh_1_78[k]
                      + di_1_99[k];

            t_79[k] = -ab_x[k] * dh_1_79[k]
                      + di_1_100[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dh_1_80, dh_1_81, dh_1_82, \
                         dh_1_83, dh_1_84, di_1_101, di_1_102, di_1_103, di_1_104, \
                         di_1_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * dh_1_80[k]
                      + di_1_101[k];

            t_81[k] = -ab_x[k] * dh_1_81[k]
                      + di_1_102[k];

            t_82[k] = -ab_x[k] * dh_1_82[k]
                      + di_1_103[k];

            t_83[k] = -ab_x[k] * dh_1_83[k]
                      + di_1_104[k];

            t_84[k] = -ab_x[k] * dh_1_84[k]
                      + di_1_112[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dh_1_85, dh_1_86, dh_1_87, \
                         dh_1_88, dh_1_89, di_1_113, di_1_114, di_1_115, di_1_116, \
                         di_1_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * dh_1_85[k]
                      + di_1_113[k];

            t_86[k] = -ab_x[k] * dh_1_86[k]
                      + di_1_114[k];

            t_87[k] = -ab_x[k] * dh_1_87[k]
                      + di_1_115[k];

            t_88[k] = -ab_x[k] * dh_1_88[k]
                      + di_1_116[k];

            t_89[k] = -ab_x[k] * dh_1_89[k]
                      + di_1_117[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, dh_1_90, dh_1_91, dh_1_92, \
                         dh_1_93, dh_1_94, di_1_118, di_1_119, di_1_120, di_1_121, \
                         di_1_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * dh_1_90[k]
                      + di_1_118[k];

            t_91[k] = -ab_x[k] * dh_1_91[k]
                      + di_1_119[k];

            t_92[k] = -ab_x[k] * dh_1_92[k]
                      + di_1_120[k];

            t_93[k] = -ab_x[k] * dh_1_93[k]
                      + di_1_121[k];

            t_94[k] = -ab_x[k] * dh_1_94[k]
                      + di_1_122[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, dh_1_95, dh_1_96, dh_1_97, \
                         dh_1_98, dh_1_99, di_1_123, di_1_124, di_1_125, di_1_126, \
                         di_1_127 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * dh_1_95[k]
                      + di_1_123[k];

            t_96[k] = -ab_x[k] * dh_1_96[k]
                      + di_1_124[k];

            t_97[k] = -ab_x[k] * dh_1_97[k]
                      + di_1_125[k];

            t_98[k] = -ab_x[k] * dh_1_98[k]
                      + di_1_126[k];

            t_99[k] = -ab_x[k] * dh_1_99[k]
                      + di_1_127[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, dh_1_100, dh_1_101, \
                         dh_1_102, dh_1_103, dh_1_104, di_1_128, di_1_129, di_1_130, di_1_131, \
                         di_1_132 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * dh_1_100[k]
                       + di_1_128[k];

            t_101[k] = -ab_x[k] * dh_1_101[k]
                       + di_1_129[k];

            t_102[k] = -ab_x[k] * dh_1_102[k]
                       + di_1_130[k];

            t_103[k] = -ab_x[k] * dh_1_103[k]
                       + di_1_131[k];

            t_104[k] = -ab_x[k] * dh_1_104[k]
                       + di_1_132[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, dh_1_105, dh_1_106, \
                         dh_1_107, dh_1_108, dh_1_109, di_1_140, di_1_141, di_1_142, di_1_143, \
                         di_1_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * dh_1_105[k]
                       + di_1_140[k];

            t_106[k] = -ab_x[k] * dh_1_106[k]
                       + di_1_141[k];

            t_107[k] = -ab_x[k] * dh_1_107[k]
                       + di_1_142[k];

            t_108[k] = -ab_x[k] * dh_1_108[k]
                       + di_1_143[k];

            t_109[k] = -ab_x[k] * dh_1_109[k]
                       + di_1_144[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, dh_1_110, dh_1_111, \
                         dh_1_112, dh_1_113, dh_1_114, di_1_145, di_1_146, di_1_147, di_1_148, \
                         di_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * dh_1_110[k]
                       + di_1_145[k];

            t_111[k] = -ab_x[k] * dh_1_111[k]
                       + di_1_146[k];

            t_112[k] = -ab_x[k] * dh_1_112[k]
                       + di_1_147[k];

            t_113[k] = -ab_x[k] * dh_1_113[k]
                       + di_1_148[k];

            t_114[k] = -ab_x[k] * dh_1_114[k]
                       + di_1_149[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, dh_1_115, dh_1_116, \
                         dh_1_117, dh_1_118, dh_1_119, di_1_150, di_1_151, di_1_152, di_1_153, \
                         di_1_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * dh_1_115[k]
                       + di_1_150[k];

            t_116[k] = -ab_x[k] * dh_1_116[k]
                       + di_1_151[k];

            t_117[k] = -ab_x[k] * dh_1_117[k]
                       + di_1_152[k];

            t_118[k] = -ab_x[k] * dh_1_118[k]
                       + di_1_153[k];

            t_119[k] = -ab_x[k] * dh_1_119[k]
                       + di_1_154[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, dh_1_120, dh_1_121, \
                         dh_1_122, dh_1_123, dh_1_124, di_1_155, di_1_156, di_1_157, di_1_158, \
                         di_1_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * dh_1_120[k]
                       + di_1_155[k];

            t_121[k] = -ab_x[k] * dh_1_121[k]
                       + di_1_156[k];

            t_122[k] = -ab_x[k] * dh_1_122[k]
                       + di_1_157[k];

            t_123[k] = -ab_x[k] * dh_1_123[k]
                       + di_1_158[k];

            t_124[k] = -ab_x[k] * dh_1_124[k]
                       + di_1_159[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, ab_x, ab_y, dh_1_63, dh_1_64, dh_1_125, dh_0_63, \
                         dh_0_64, di_1_85, di_1_87, di_1_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * dh_1_125[k]
                       + di_1_160[k];

            t_126[k] = -ab_y[k] * dh_1_63[k]
                       + dh_0_63[k]
                       + di_1_85[k];

            t_127[k] = -ab_y[k] * dh_1_64[k]
                       + dh_0_64[k]
                       + di_1_87[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, ab_y, dh_1_65, dh_1_66, dh_1_67, dh_0_65, \
                         dh_0_66, dh_0_67, di_1_88, di_1_90, di_1_91 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = -ab_y[k] * dh_1_65[k]
                       + dh_0_65[k]
                       + di_1_88[k];

            t_129[k] = -ab_y[k] * dh_1_66[k]
                       + dh_0_66[k]
                       + di_1_90[k];

            t_130[k] = -ab_y[k] * dh_1_67[k]
                       + dh_0_67[k]
                       + di_1_91[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, ab_y, dh_1_68, dh_1_69, dh_1_70, dh_0_68, \
                         dh_0_69, dh_0_70, di_1_92, di_1_94, di_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = -ab_y[k] * dh_1_68[k]
                       + dh_0_68[k]
                       + di_1_92[k];

            t_132[k] = -ab_y[k] * dh_1_69[k]
                       + dh_0_69[k]
                       + di_1_94[k];

            t_133[k] = -ab_y[k] * dh_1_70[k]
                       + dh_0_70[k]
                       + di_1_95[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, ab_y, dh_1_71, dh_1_72, dh_1_73, dh_0_71, \
                         dh_0_72, dh_0_73, di_1_96, di_1_97, di_1_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = -ab_y[k] * dh_1_71[k]
                       + dh_0_71[k]
                       + di_1_96[k];

            t_135[k] = -ab_y[k] * dh_1_72[k]
                       + dh_0_72[k]
                       + di_1_97[k];

            t_136[k] = -ab_y[k] * dh_1_73[k]
                       + dh_0_73[k]
                       + di_1_99[k];
        }

#pragma omp simd aligned(t_137, t_138, t_139, ab_y, dh_1_74, dh_1_75, dh_1_76, dh_0_74, \
                         dh_0_75, dh_0_76, di_1_100, di_1_101, \
                         di_1_102 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_137[k] = -ab_y[k] * dh_1_74[k]
                       + dh_0_74[k]
                       + di_1_100[k];

            t_138[k] = -ab_y[k] * dh_1_75[k]
                       + dh_0_75[k]
                       + di_1_101[k];

            t_139[k] = -ab_y[k] * dh_1_76[k]
                       + dh_0_76[k]
                       + di_1_102[k];
        }
    }
}

static auto
compute_hrr_geom_010y_fh_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t dh_1, const size_t dh_0,
                                const size_t di_1, const size_t ncomps,
                                const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dh_1_77 = buffer.data(dh_1 + 77 * ncomps + c);
        const auto *dh_1_78 = buffer.data(dh_1 + 78 * ncomps + c);
        const auto *dh_1_79 = buffer.data(dh_1 + 79 * ncomps + c);
        const auto *dh_1_80 = buffer.data(dh_1 + 80 * ncomps + c);
        const auto *dh_1_81 = buffer.data(dh_1 + 81 * ncomps + c);
        const auto *dh_1_82 = buffer.data(dh_1 + 82 * ncomps + c);
        const auto *dh_1_83 = buffer.data(dh_1 + 83 * ncomps + c);
        const auto *dh_1_84 = buffer.data(dh_1 + 84 * ncomps + c);
        const auto *dh_1_85 = buffer.data(dh_1 + 85 * ncomps + c);
        const auto *dh_1_86 = buffer.data(dh_1 + 86 * ncomps + c);
        const auto *dh_1_87 = buffer.data(dh_1 + 87 * ncomps + c);
        const auto *dh_1_88 = buffer.data(dh_1 + 88 * ncomps + c);
        const auto *dh_1_89 = buffer.data(dh_1 + 89 * ncomps + c);
        const auto *dh_1_90 = buffer.data(dh_1 + 90 * ncomps + c);
        const auto *dh_1_91 = buffer.data(dh_1 + 91 * ncomps + c);
        const auto *dh_1_92 = buffer.data(dh_1 + 92 * ncomps + c);
        const auto *dh_1_93 = buffer.data(dh_1 + 93 * ncomps + c);
        const auto *dh_1_94 = buffer.data(dh_1 + 94 * ncomps + c);
        const auto *dh_1_95 = buffer.data(dh_1 + 95 * ncomps + c);
        const auto *dh_1_96 = buffer.data(dh_1 + 96 * ncomps + c);
        const auto *dh_1_97 = buffer.data(dh_1 + 97 * ncomps + c);
        const auto *dh_1_98 = buffer.data(dh_1 + 98 * ncomps + c);
        const auto *dh_1_99 = buffer.data(dh_1 + 99 * ncomps + c);
        const auto *dh_1_100 = buffer.data(dh_1 + 100 * ncomps + c);
        const auto *dh_1_101 = buffer.data(dh_1 + 101 * ncomps + c);
        const auto *dh_1_102 = buffer.data(dh_1 + 102 * ncomps + c);
        const auto *dh_1_103 = buffer.data(dh_1 + 103 * ncomps + c);
        const auto *dh_1_104 = buffer.data(dh_1 + 104 * ncomps + c);
        const auto *dh_1_105 = buffer.data(dh_1 + 105 * ncomps + c);
        const auto *dh_1_106 = buffer.data(dh_1 + 106 * ncomps + c);
        const auto *dh_1_107 = buffer.data(dh_1 + 107 * ncomps + c);
        const auto *dh_1_108 = buffer.data(dh_1 + 108 * ncomps + c);
        const auto *dh_1_109 = buffer.data(dh_1 + 109 * ncomps + c);
        const auto *dh_1_110 = buffer.data(dh_1 + 110 * ncomps + c);
        const auto *dh_1_111 = buffer.data(dh_1 + 111 * ncomps + c);
        const auto *dh_1_112 = buffer.data(dh_1 + 112 * ncomps + c);
        const auto *dh_1_113 = buffer.data(dh_1 + 113 * ncomps + c);
        const auto *dh_1_114 = buffer.data(dh_1 + 114 * ncomps + c);
        const auto *dh_1_115 = buffer.data(dh_1 + 115 * ncomps + c);
        const auto *dh_1_116 = buffer.data(dh_1 + 116 * ncomps + c);
        const auto *dh_1_117 = buffer.data(dh_1 + 117 * ncomps + c);
        const auto *dh_1_118 = buffer.data(dh_1 + 118 * ncomps + c);
        const auto *dh_1_119 = buffer.data(dh_1 + 119 * ncomps + c);
        const auto *dh_1_120 = buffer.data(dh_1 + 120 * ncomps + c);
        const auto *dh_1_121 = buffer.data(dh_1 + 121 * ncomps + c);
        const auto *dh_1_122 = buffer.data(dh_1 + 122 * ncomps + c);
        const auto *dh_1_123 = buffer.data(dh_1 + 123 * ncomps + c);
        const auto *dh_1_124 = buffer.data(dh_1 + 124 * ncomps + c);
        const auto *dh_1_125 = buffer.data(dh_1 + 125 * ncomps + c);

        const auto *dh_0_77 = buffer.data(dh_0 + 77 * ncomps + c);
        const auto *dh_0_78 = buffer.data(dh_0 + 78 * ncomps + c);
        const auto *dh_0_79 = buffer.data(dh_0 + 79 * ncomps + c);
        const auto *dh_0_80 = buffer.data(dh_0 + 80 * ncomps + c);
        const auto *dh_0_81 = buffer.data(dh_0 + 81 * ncomps + c);
        const auto *dh_0_82 = buffer.data(dh_0 + 82 * ncomps + c);
        const auto *dh_0_83 = buffer.data(dh_0 + 83 * ncomps + c);
        const auto *dh_0_84 = buffer.data(dh_0 + 84 * ncomps + c);
        const auto *dh_0_85 = buffer.data(dh_0 + 85 * ncomps + c);
        const auto *dh_0_86 = buffer.data(dh_0 + 86 * ncomps + c);
        const auto *dh_0_87 = buffer.data(dh_0 + 87 * ncomps + c);
        const auto *dh_0_88 = buffer.data(dh_0 + 88 * ncomps + c);
        const auto *dh_0_89 = buffer.data(dh_0 + 89 * ncomps + c);
        const auto *dh_0_90 = buffer.data(dh_0 + 90 * ncomps + c);
        const auto *dh_0_91 = buffer.data(dh_0 + 91 * ncomps + c);
        const auto *dh_0_92 = buffer.data(dh_0 + 92 * ncomps + c);
        const auto *dh_0_93 = buffer.data(dh_0 + 93 * ncomps + c);
        const auto *dh_0_94 = buffer.data(dh_0 + 94 * ncomps + c);
        const auto *dh_0_95 = buffer.data(dh_0 + 95 * ncomps + c);
        const auto *dh_0_96 = buffer.data(dh_0 + 96 * ncomps + c);
        const auto *dh_0_97 = buffer.data(dh_0 + 97 * ncomps + c);
        const auto *dh_0_98 = buffer.data(dh_0 + 98 * ncomps + c);
        const auto *dh_0_99 = buffer.data(dh_0 + 99 * ncomps + c);
        const auto *dh_0_100 = buffer.data(dh_0 + 100 * ncomps + c);
        const auto *dh_0_101 = buffer.data(dh_0 + 101 * ncomps + c);
        const auto *dh_0_102 = buffer.data(dh_0 + 102 * ncomps + c);
        const auto *dh_0_103 = buffer.data(dh_0 + 103 * ncomps + c);
        const auto *dh_0_104 = buffer.data(dh_0 + 104 * ncomps + c);
        const auto *dh_0_105 = buffer.data(dh_0 + 105 * ncomps + c);
        const auto *dh_0_106 = buffer.data(dh_0 + 106 * ncomps + c);
        const auto *dh_0_107 = buffer.data(dh_0 + 107 * ncomps + c);
        const auto *dh_0_108 = buffer.data(dh_0 + 108 * ncomps + c);
        const auto *dh_0_109 = buffer.data(dh_0 + 109 * ncomps + c);
        const auto *dh_0_110 = buffer.data(dh_0 + 110 * ncomps + c);
        const auto *dh_0_111 = buffer.data(dh_0 + 111 * ncomps + c);
        const auto *dh_0_112 = buffer.data(dh_0 + 112 * ncomps + c);
        const auto *dh_0_113 = buffer.data(dh_0 + 113 * ncomps + c);
        const auto *dh_0_114 = buffer.data(dh_0 + 114 * ncomps + c);
        const auto *dh_0_115 = buffer.data(dh_0 + 115 * ncomps + c);
        const auto *dh_0_116 = buffer.data(dh_0 + 116 * ncomps + c);
        const auto *dh_0_117 = buffer.data(dh_0 + 117 * ncomps + c);
        const auto *dh_0_118 = buffer.data(dh_0 + 118 * ncomps + c);
        const auto *dh_0_119 = buffer.data(dh_0 + 119 * ncomps + c);
        const auto *dh_0_120 = buffer.data(dh_0 + 120 * ncomps + c);
        const auto *dh_0_121 = buffer.data(dh_0 + 121 * ncomps + c);
        const auto *dh_0_122 = buffer.data(dh_0 + 122 * ncomps + c);
        const auto *dh_0_123 = buffer.data(dh_0 + 123 * ncomps + c);
        const auto *dh_0_124 = buffer.data(dh_0 + 124 * ncomps + c);
        const auto *dh_0_125 = buffer.data(dh_0 + 125 * ncomps + c);

        const auto *di_1_103 = buffer.data(di_1 + 103 * ncomps + c);
        const auto *di_1_105 = buffer.data(di_1 + 105 * ncomps + c);
        const auto *di_1_106 = buffer.data(di_1 + 106 * ncomps + c);
        const auto *di_1_107 = buffer.data(di_1 + 107 * ncomps + c);
        const auto *di_1_108 = buffer.data(di_1 + 108 * ncomps + c);
        const auto *di_1_109 = buffer.data(di_1 + 109 * ncomps + c);
        const auto *di_1_110 = buffer.data(di_1 + 110 * ncomps + c);
        const auto *di_1_113 = buffer.data(di_1 + 113 * ncomps + c);
        const auto *di_1_115 = buffer.data(di_1 + 115 * ncomps + c);
        const auto *di_1_116 = buffer.data(di_1 + 116 * ncomps + c);
        const auto *di_1_118 = buffer.data(di_1 + 118 * ncomps + c);
        const auto *di_1_119 = buffer.data(di_1 + 119 * ncomps + c);
        const auto *di_1_120 = buffer.data(di_1 + 120 * ncomps + c);
        const auto *di_1_122 = buffer.data(di_1 + 122 * ncomps + c);
        const auto *di_1_123 = buffer.data(di_1 + 123 * ncomps + c);
        const auto *di_1_124 = buffer.data(di_1 + 124 * ncomps + c);
        const auto *di_1_125 = buffer.data(di_1 + 125 * ncomps + c);
        const auto *di_1_127 = buffer.data(di_1 + 127 * ncomps + c);
        const auto *di_1_128 = buffer.data(di_1 + 128 * ncomps + c);
        const auto *di_1_129 = buffer.data(di_1 + 129 * ncomps + c);
        const auto *di_1_130 = buffer.data(di_1 + 130 * ncomps + c);
        const auto *di_1_131 = buffer.data(di_1 + 131 * ncomps + c);
        const auto *di_1_133 = buffer.data(di_1 + 133 * ncomps + c);
        const auto *di_1_134 = buffer.data(di_1 + 134 * ncomps + c);
        const auto *di_1_135 = buffer.data(di_1 + 135 * ncomps + c);
        const auto *di_1_136 = buffer.data(di_1 + 136 * ncomps + c);
        const auto *di_1_137 = buffer.data(di_1 + 137 * ncomps + c);
        const auto *di_1_138 = buffer.data(di_1 + 138 * ncomps + c);
        const auto *di_1_141 = buffer.data(di_1 + 141 * ncomps + c);
        const auto *di_1_142 = buffer.data(di_1 + 142 * ncomps + c);
        const auto *di_1_143 = buffer.data(di_1 + 143 * ncomps + c);
        const auto *di_1_144 = buffer.data(di_1 + 144 * ncomps + c);
        const auto *di_1_145 = buffer.data(di_1 + 145 * ncomps + c);
        const auto *di_1_146 = buffer.data(di_1 + 146 * ncomps + c);
        const auto *di_1_147 = buffer.data(di_1 + 147 * ncomps + c);
        const auto *di_1_148 = buffer.data(di_1 + 148 * ncomps + c);
        const auto *di_1_149 = buffer.data(di_1 + 149 * ncomps + c);
        const auto *di_1_150 = buffer.data(di_1 + 150 * ncomps + c);
        const auto *di_1_151 = buffer.data(di_1 + 151 * ncomps + c);
        const auto *di_1_152 = buffer.data(di_1 + 152 * ncomps + c);
        const auto *di_1_153 = buffer.data(di_1 + 153 * ncomps + c);
        const auto *di_1_154 = buffer.data(di_1 + 154 * ncomps + c);
        const auto *di_1_155 = buffer.data(di_1 + 155 * ncomps + c);
        const auto *di_1_156 = buffer.data(di_1 + 156 * ncomps + c);
        const auto *di_1_157 = buffer.data(di_1 + 157 * ncomps + c);
        const auto *di_1_158 = buffer.data(di_1 + 158 * ncomps + c);
        const auto *di_1_159 = buffer.data(di_1 + 159 * ncomps + c);
        const auto *di_1_160 = buffer.data(di_1 + 160 * ncomps + c);
        const auto *di_1_161 = buffer.data(di_1 + 161 * ncomps + c);
        const auto *di_1_162 = buffer.data(di_1 + 162 * ncomps + c);
        const auto *di_1_163 = buffer.data(di_1 + 163 * ncomps + c);
        const auto *di_1_164 = buffer.data(di_1 + 164 * ncomps + c);
        const auto *di_1_165 = buffer.data(di_1 + 165 * ncomps + c);
        const auto *di_1_166 = buffer.data(di_1 + 166 * ncomps + c);
        const auto *di_1_167 = buffer.data(di_1 + 167 * ncomps + c);

#pragma omp simd aligned(t_140, t_141, t_142, ab_y, dh_1_77, dh_1_78, dh_1_79, dh_0_77, \
                         dh_0_78, dh_0_79, di_1_103, di_1_105, \
                         di_1_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_y[k] * dh_1_77[k]
                       + dh_0_77[k]
                       + di_1_103[k];

            t_141[k] = -ab_y[k] * dh_1_78[k]
                       + dh_0_78[k]
                       + di_1_105[k];

            t_142[k] = -ab_y[k] * dh_1_79[k]
                       + dh_0_79[k]
                       + di_1_106[k];
        }

#pragma omp simd aligned(t_143, t_144, t_145, ab_y, dh_1_80, dh_1_81, dh_1_82, dh_0_80, \
                         dh_0_81, dh_0_82, di_1_107, di_1_108, \
                         di_1_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_143[k] = -ab_y[k] * dh_1_80[k]
                       + dh_0_80[k]
                       + di_1_107[k];

            t_144[k] = -ab_y[k] * dh_1_81[k]
                       + dh_0_81[k]
                       + di_1_108[k];

            t_145[k] = -ab_y[k] * dh_1_82[k]
                       + dh_0_82[k]
                       + di_1_109[k];
        }

#pragma omp simd aligned(t_146, t_147, t_148, ab_y, dh_1_83, dh_1_84, dh_1_85, dh_0_83, \
                         dh_0_84, dh_0_85, di_1_110, di_1_113, \
                         di_1_115 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_146[k] = -ab_y[k] * dh_1_83[k]
                       + dh_0_83[k]
                       + di_1_110[k];

            t_147[k] = -ab_y[k] * dh_1_84[k]
                       + dh_0_84[k]
                       + di_1_113[k];

            t_148[k] = -ab_y[k] * dh_1_85[k]
                       + dh_0_85[k]
                       + di_1_115[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, ab_y, dh_1_86, dh_1_87, dh_1_88, dh_0_86, \
                         dh_0_87, dh_0_88, di_1_116, di_1_118, \
                         di_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = -ab_y[k] * dh_1_86[k]
                       + dh_0_86[k]
                       + di_1_116[k];

            t_150[k] = -ab_y[k] * dh_1_87[k]
                       + dh_0_87[k]
                       + di_1_118[k];

            t_151[k] = -ab_y[k] * dh_1_88[k]
                       + dh_0_88[k]
                       + di_1_119[k];
        }

#pragma omp simd aligned(t_152, t_153, t_154, ab_y, dh_1_89, dh_1_90, dh_1_91, dh_0_89, \
                         dh_0_90, dh_0_91, di_1_120, di_1_122, \
                         di_1_123 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_152[k] = -ab_y[k] * dh_1_89[k]
                       + dh_0_89[k]
                       + di_1_120[k];

            t_153[k] = -ab_y[k] * dh_1_90[k]
                       + dh_0_90[k]
                       + di_1_122[k];

            t_154[k] = -ab_y[k] * dh_1_91[k]
                       + dh_0_91[k]
                       + di_1_123[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, ab_y, dh_1_92, dh_1_93, dh_1_94, dh_0_92, \
                         dh_0_93, dh_0_94, di_1_124, di_1_125, \
                         di_1_127 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_y[k] * dh_1_92[k]
                       + dh_0_92[k]
                       + di_1_124[k];

            t_156[k] = -ab_y[k] * dh_1_93[k]
                       + dh_0_93[k]
                       + di_1_125[k];

            t_157[k] = -ab_y[k] * dh_1_94[k]
                       + dh_0_94[k]
                       + di_1_127[k];
        }

#pragma omp simd aligned(t_158, t_159, t_160, ab_y, dh_1_95, dh_1_96, dh_1_97, dh_0_95, \
                         dh_0_96, dh_0_97, di_1_128, di_1_129, \
                         di_1_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_158[k] = -ab_y[k] * dh_1_95[k]
                       + dh_0_95[k]
                       + di_1_128[k];

            t_159[k] = -ab_y[k] * dh_1_96[k]
                       + dh_0_96[k]
                       + di_1_129[k];

            t_160[k] = -ab_y[k] * dh_1_97[k]
                       + dh_0_97[k]
                       + di_1_130[k];
        }

#pragma omp simd aligned(t_161, t_162, t_163, ab_y, dh_1_98, dh_1_99, dh_1_100, dh_0_98, \
                         dh_0_99, dh_0_100, di_1_131, di_1_133, \
                         di_1_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_161[k] = -ab_y[k] * dh_1_98[k]
                       + dh_0_98[k]
                       + di_1_131[k];

            t_162[k] = -ab_y[k] * dh_1_99[k]
                       + dh_0_99[k]
                       + di_1_133[k];

            t_163[k] = -ab_y[k] * dh_1_100[k]
                       + dh_0_100[k]
                       + di_1_134[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, ab_y, dh_1_101, dh_1_102, dh_1_103, dh_0_101, \
                         dh_0_102, dh_0_103, di_1_135, di_1_136, \
                         di_1_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = -ab_y[k] * dh_1_101[k]
                       + dh_0_101[k]
                       + di_1_135[k];

            t_165[k] = -ab_y[k] * dh_1_102[k]
                       + dh_0_102[k]
                       + di_1_136[k];

            t_166[k] = -ab_y[k] * dh_1_103[k]
                       + dh_0_103[k]
                       + di_1_137[k];
        }

#pragma omp simd aligned(t_167, t_168, t_169, ab_y, dh_1_104, dh_1_105, dh_1_106, dh_0_104, \
                         dh_0_105, dh_0_106, di_1_138, di_1_141, \
                         di_1_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = -ab_y[k] * dh_1_104[k]
                       + dh_0_104[k]
                       + di_1_138[k];

            t_168[k] = -ab_y[k] * dh_1_105[k]
                       + dh_0_105[k]
                       + di_1_141[k];

            t_169[k] = -ab_y[k] * dh_1_106[k]
                       + dh_0_106[k]
                       + di_1_143[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, ab_y, dh_1_107, dh_1_108, dh_1_109, dh_0_107, \
                         dh_0_108, dh_0_109, di_1_144, di_1_146, \
                         di_1_147 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_y[k] * dh_1_107[k]
                       + dh_0_107[k]
                       + di_1_144[k];

            t_171[k] = -ab_y[k] * dh_1_108[k]
                       + dh_0_108[k]
                       + di_1_146[k];

            t_172[k] = -ab_y[k] * dh_1_109[k]
                       + dh_0_109[k]
                       + di_1_147[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, ab_y, dh_1_110, dh_1_111, dh_1_112, dh_0_110, \
                         dh_0_111, dh_0_112, di_1_148, di_1_150, \
                         di_1_151 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = -ab_y[k] * dh_1_110[k]
                       + dh_0_110[k]
                       + di_1_148[k];

            t_174[k] = -ab_y[k] * dh_1_111[k]
                       + dh_0_111[k]
                       + di_1_150[k];

            t_175[k] = -ab_y[k] * dh_1_112[k]
                       + dh_0_112[k]
                       + di_1_151[k];
        }

#pragma omp simd aligned(t_176, t_177, t_178, ab_y, dh_1_113, dh_1_114, dh_1_115, dh_0_113, \
                         dh_0_114, dh_0_115, di_1_152, di_1_153, \
                         di_1_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_176[k] = -ab_y[k] * dh_1_113[k]
                       + dh_0_113[k]
                       + di_1_152[k];

            t_177[k] = -ab_y[k] * dh_1_114[k]
                       + dh_0_114[k]
                       + di_1_153[k];

            t_178[k] = -ab_y[k] * dh_1_115[k]
                       + dh_0_115[k]
                       + di_1_155[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, ab_y, dh_1_116, dh_1_117, dh_1_118, dh_0_116, \
                         dh_0_117, dh_0_118, di_1_156, di_1_157, \
                         di_1_158 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = -ab_y[k] * dh_1_116[k]
                       + dh_0_116[k]
                       + di_1_156[k];

            t_180[k] = -ab_y[k] * dh_1_117[k]
                       + dh_0_117[k]
                       + di_1_157[k];

            t_181[k] = -ab_y[k] * dh_1_118[k]
                       + dh_0_118[k]
                       + di_1_158[k];
        }

#pragma omp simd aligned(t_182, t_183, t_184, ab_y, dh_1_119, dh_1_120, dh_1_121, dh_0_119, \
                         dh_0_120, dh_0_121, di_1_159, di_1_161, \
                         di_1_162 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_182[k] = -ab_y[k] * dh_1_119[k]
                       + dh_0_119[k]
                       + di_1_159[k];

            t_183[k] = -ab_y[k] * dh_1_120[k]
                       + dh_0_120[k]
                       + di_1_161[k];

            t_184[k] = -ab_y[k] * dh_1_121[k]
                       + dh_0_121[k]
                       + di_1_162[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, ab_y, dh_1_122, dh_1_123, dh_1_124, dh_0_122, \
                         dh_0_123, dh_0_124, di_1_163, di_1_164, \
                         di_1_165 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_y[k] * dh_1_122[k]
                       + dh_0_122[k]
                       + di_1_163[k];

            t_186[k] = -ab_y[k] * dh_1_123[k]
                       + dh_0_123[k]
                       + di_1_164[k];

            t_187[k] = -ab_y[k] * dh_1_124[k]
                       + dh_0_124[k]
                       + di_1_165[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, ab_y, ab_z, dh_1_105, dh_1_106, dh_1_107, \
                         dh_1_125, dh_0_125, di_1_142, di_1_144, di_1_145, \
                         di_1_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = -ab_y[k] * dh_1_125[k]
                       + dh_0_125[k]
                       + di_1_166[k];

            t_189[k] = -ab_z[k] * dh_1_105[k]
                       + di_1_142[k];

            t_190[k] = -ab_z[k] * dh_1_106[k]
                       + di_1_144[k];

            t_191[k] = -ab_z[k] * dh_1_107[k]
                       + di_1_145[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ab_z, dh_1_108, dh_1_109, \
                         dh_1_110, dh_1_111, dh_1_112, di_1_147, di_1_148, di_1_149, di_1_151, \
                         di_1_152 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = -ab_z[k] * dh_1_108[k]
                       + di_1_147[k];

            t_193[k] = -ab_z[k] * dh_1_109[k]
                       + di_1_148[k];

            t_194[k] = -ab_z[k] * dh_1_110[k]
                       + di_1_149[k];

            t_195[k] = -ab_z[k] * dh_1_111[k]
                       + di_1_151[k];

            t_196[k] = -ab_z[k] * dh_1_112[k]
                       + di_1_152[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_z, dh_1_113, dh_1_114, \
                         dh_1_115, dh_1_116, dh_1_117, di_1_153, di_1_154, di_1_156, di_1_157, \
                         di_1_158 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = -ab_z[k] * dh_1_113[k]
                       + di_1_153[k];

            t_198[k] = -ab_z[k] * dh_1_114[k]
                       + di_1_154[k];

            t_199[k] = -ab_z[k] * dh_1_115[k]
                       + di_1_156[k];

            t_200[k] = -ab_z[k] * dh_1_116[k]
                       + di_1_157[k];

            t_201[k] = -ab_z[k] * dh_1_117[k]
                       + di_1_158[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, ab_z, dh_1_118, dh_1_119, \
                         dh_1_120, dh_1_121, dh_1_122, di_1_159, di_1_160, di_1_162, di_1_163, \
                         di_1_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = -ab_z[k] * dh_1_118[k]
                       + di_1_159[k];

            t_203[k] = -ab_z[k] * dh_1_119[k]
                       + di_1_160[k];

            t_204[k] = -ab_z[k] * dh_1_120[k]
                       + di_1_162[k];

            t_205[k] = -ab_z[k] * dh_1_121[k]
                       + di_1_163[k];

            t_206[k] = -ab_z[k] * dh_1_122[k]
                       + di_1_164[k];
        }

#pragma omp simd aligned(t_207, t_208, t_209, ab_z, dh_1_123, dh_1_124, dh_1_125, di_1_165, \
                         di_1_166, di_1_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_207[k] = -ab_z[k] * dh_1_123[k]
                       + di_1_165[k];

            t_208[k] = -ab_z[k] * dh_1_124[k]
                       + di_1_166[k];

            t_209[k] = -ab_z[k] * dh_1_125[k]
                       + di_1_167[k];
        }
    }
}

auto
compute_hrr_geom_010y_fh(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t dh_1, const size_t dh_0,
                         const size_t di_1, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_geom_010y_fh_piece0(buffer, coordinates, target, dh_1, dh_0, di_1, ncomps,
                                    nmax);

    compute_hrr_geom_010y_fh_piece1(buffer, coordinates, target, dh_1, dh_0, di_1, ncomps,
                                    nmax);
}

}  // namespace simdtrf
